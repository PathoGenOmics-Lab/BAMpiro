#!/usr/bin/env python3
import sys
import json
import argparse
import datetime
import gzip
import os

# List of keys expected by the legacy system
# Added LIN_COUNTS at the end for the specific counts
ORDERED_KEYS = [
    'DAT_OUT', 'GID_IN', 'FQ_IN', 'PHE_OUT', 'FV_IN', 'FV_OUT',
    'FQC_IN', 'TRD_OUT', 'RDL_OUT', 'LNG_OUT', 'TRM_IN', 'TRM_OUT',
    'ASV_OUT', 'FSV_OUT', 'RSV_OUT', 'ADP_OUT', 'SQP_IN', 'SIN_OUT',
    'PM_OUT', 'PA_OUT', 'PD_OUT', 'MG_OUT', 'BWA_IN', 'SM1_IN',
    'MD_IN', 'SM2_IN', 'IDX1_IN', 'RTC_IN', 'IDR_IN', 'BAMF_IN',
    'IDX2_IN', 'IDX3_IN', 'SFS_IN', 'RDC_OUT', 'DUP_OUT', 'PP_OUT',
    'MAP_OUT', 'MAP_PCT', 'BQC_IN', 'COV_OUT', 'COV_BREADTH',
    'COV_MEDIAN', 'COV_BREADTH5', 'COV_BREADTH10', 'COV_EVENNESS', 'GENOME_LEN',
    'UNMAP_PCT', 'SINGLE_PCT', 'MQ0_PCT', 'ERR_RATE', 'ISIZE_MEAN', 'ISIZE_SD',
    'Q20_PCT', 'Q30_PCT', 'GC_PCT', 'STD_OUT',
    'NC_OUT', 'N7_OUT',
    'VC_IN', 'BPF_OUT', 'VFL_OUT', 'VAR_OUT', 'SNP_OUT', 'HOM_OUT',
    'HET_OUT', 'HOMO_IND_OUT', 'TITV_OUT', 'SNP_PROFILE', 'HET_PROFILE', 'INDEL_PROFILE',
    'ANN_IN',
    'ANN_PRESENT', 'ANN_ANNOTATED_FRAC', 'ANN_TOTAL', 'ANN_ANNOTATED',
    'ANN_IMPACT_HIGH', 'ANN_IMPACT_MODERATE', 'ANN_IMPACT_LOW', 'ANN_IMPACT_MODIFIER',
    'ANN_MISSENSE', 'ANN_SYNON', 'ANN_STOP_GAINED', 'ANN_STOP_LOST', 'ANN_START_LOST',
    'ANN_FRAMESHIFT', 'ANN_INFRAME_INDEL', 'ANN_SPLICE', 'ANN_INTERGENIC', 'ANN_REGULATORY',
    'ANN_LOF', 'ANN_CODING', 'ANN_MISSENSE_SILENT',
    'ANN_WARN', 'ANN_DB_ERROR', 'ANN_GENE_BURDEN',
    'LIN_OUT', 'LIN_COUNTS',
    'DR_OUT', 'DR_COUNTS', 'ANC_OUT'
]

# Transition pairs (unordered). Everything else between distinct bases is a transversion.
_TRANSITIONS = {frozenset(('A', 'G')), frozenset(('C', 'T'))}

# --- snpEff ANN parsing (functional annotation features): fixed subfield layout, 0-indexed ---
_ANN_EFFECT = 1     # Annotation (SO effect terms, '&'-joined)
_ANN_IMPACT = 2     # Annotation_Impact: HIGH|MODERATE|LOW|MODIFIER
_ANN_GENE   = 3     # Gene_Name
_ANN_ERRORS = 15    # ERRORS/WARNINGS/INFO subfield (may be empty)
_IMPACTS = ('HIGH', 'MODERATE', 'LOW', 'MODIFIER')
# Whole-VCF DB failures snpEff stamps onto the ERRORS field of every record when the VCF contig name
# does not match the built database (wrong reference / bad GFF->database build).
_DB_ERROR_TOKENS = ('ERROR_CHROMOSOME_NOT_FOUND', 'ERROR_OUT_OF_CHROMOSOME_RANGE')
_GENE_BURDEN_CAP = 40   # bound the per-sample gene list so a hypermutated sample never blows up the log line
_EFF_CLASSES = ('missense', 'synonymous', 'stop_gained', 'stop_lost', 'start_lost',
                'frameshift', 'inframe_indel', 'splice', 'intergenic', 'regulatory', 'other')


def _effect_class(e):
    """Map a lowercased, first-'&'-token snpEff effect string to a compact CLASS.
    Order matters: most specific first. Organism-agnostic (pure SO-term string rules)."""
    if 'frameshift' in e:                                       return 'frameshift'
    if 'stop_gained' in e:                                      return 'stop_gained'
    if 'stop_lost' in e:                                        return 'stop_lost'
    if 'start_lost' in e:                                       return 'start_lost'
    if 'missense' in e:                                         return 'missense'
    if 'stop_retained' in e or 'start_retained' in e:           return 'synonymous'
    if 'synonymous' in e:                                       return 'synonymous'
    if ('inframe_insertion' in e or 'inframe_deletion' in e
            or 'disruptive_inframe' in e):                      return 'inframe_indel'
    if 'splice' in e:                                           return 'splice'
    if 'intergenic' in e or 'intragenic' in e:                  return 'intergenic'
    if ('upstream' in e or 'downstream' in e or 'utr' in e
            or 'regulatory' in e or 'tf_binding' in e):         return 'regulatory'
    return 'other'


def _parse_ann_field(info_str):
    """From one record's INFO string return (class, impact, gene, has_error, db_error) for the FIRST
    (most-severe) ANN entry, or None if the record has no usable ANN= (un-annotated / malformed)."""
    ann_val = None
    for field in info_str.split(';'):
        if field.startswith('ANN='):
            ann_val = field[4:]
            break
    if not ann_val:
        return None
    first = ann_val.split(',', 1)[0]        # snpEff pre-sorts by severity; first entry is primary
    sub = first.split('|')
    if len(sub) <= _ANN_IMPACT:
        return None                          # malformed -> treat as unannotated
    primary = sub[_ANN_EFFECT].split('&', 1)[0].strip().lower()
    cls = _effect_class(primary)
    impact = sub[_ANN_IMPACT].strip().upper()
    if impact not in _IMPACTS:
        impact = 'MODIFIER'
    gene = sub[_ANN_GENE] if len(sub) > _ANN_GENE else ''
    err = sub[_ANN_ERRORS] if len(sub) > _ANN_ERRORS else ''
    has_error = err.strip() != ''
    db_error = any(tok in err for tok in _DB_ERROR_TOKENS)
    return (cls, impact, gene, has_error, db_error)


def _format_gene_burden(ann_gene):
    """ann_gene: {gene: [high, mod, {cls: count}]} -> 'gene:HIGH:MOD:topEff;...' (top-40 by HIGH then MOD,
    HIGH/MOD-carrying genes only, delimiter-sanitised). Read only by the Phase-B cohort aggregation."""
    if not ann_gene:
        return "NA"
    items = []
    for gene, (high, mod, effs) in ann_gene.items():
        if high == 0 and mod == 0:
            continue
        top_eff = max(effs.items(), key=lambda kv: kv[1])[0] if effs else 'other'
        g = gene.replace(';', '_').replace(':', '_').replace('\t', '_')
        items.append((high, mod, g, top_eff))
    if not items:
        return "NA"
    items.sort(key=lambda t: (-t[0], -t[1], t[2]))
    items = items[:_GENE_BURDEN_CAP]
    return ';'.join(f"{g}:{high}:{mod}:{eff}" for (high, mod, g, eff) in items)

def parse_args():
    parser = argparse.ArgumentParser(description="Generate Legacy TB Pipeline Log")
    parser.add_argument('--sample', required=True, help="Sample ID")
    parser.add_argument('--fastp-json', required=True, help="Path to fastp.json")
    parser.add_argument('--bam-stats', required=True, help="Path to samtools stats output")
    parser.add_argument('--vcf', required=True, help="Path to final VCF")
    parser.add_argument('--ref-fai', required=True, help="Path to reference .fai index")
    parser.add_argument('--pathotypr-report', required=False, help="Path to Pathotypr lineage summary TSV output")
    parser.add_argument('--pathotypr-dr-report', required=False, help="Path to Pathotypr DR summary TSV output")
    parser.add_argument('--ancient', default='no',
                        help="'yes' if this is an ancient (aDNA) sample; carried into QC as the ancient flag.")
    return parser.parse_args()

def get_genome_size(fai_file):
    total = 0
    try:
        with open(fai_file, 'r') as f:
            for line in f:
                parts = line.split('\t')
                if len(parts) >= 2:
                    total += int(parts[1])
    except Exception as e:
        sys.stderr.write(f"[stats_to_legacy] WARN get_genome_size({fai_file}): {e}\n")
    return total   # 0 signals unknown; every consumer guards g_size > 0 and emits NA (no organism-specific fallback)

def get_bam_stats(stats_file):
    """Parse samtools stats output for mapping metrics."""
    result = {
        'total': 0,
        'mapped': 0,
        'mapped_paired': 0,
        'properly_paired': 0,
        'duplicates': 0,
        'avg_quality': 0.0,
        'insert_size_avg': 0.0,
        'insert_size_sd': 0.0,
        'bases_mapped': 0,
        'reads_unmapped': 0,
        'reads_mq0': 0,
        'error_rate': 0.0,
    }
    try:
        with open(stats_file, 'r') as f:
            for line in f:
                if not line.startswith('SN\t'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 3:
                    continue
                key = parts[1].rstrip(':')
                val = parts[2].strip()

                if key == 'raw total sequences':
                    result['total'] = int(val)
                elif key == 'reads mapped':
                    result['mapped'] = int(val)
                elif key == 'reads mapped and paired':
                    result['mapped_paired'] = int(val)
                elif key == 'reads properly paired':
                    result['properly_paired'] = int(val)
                elif key == 'reads duplicated':
                    result['duplicates'] = int(val)
                elif key == 'reads unmapped':
                    result['reads_unmapped'] = int(val)
                elif key == 'reads MQ0':
                    result['reads_mq0'] = int(val)
                elif key == 'average quality':
                    result['avg_quality'] = float(val)
                elif key == 'insert size average':
                    result['insert_size_avg'] = float(val)
                elif key == 'insert size standard deviation':
                    result['insert_size_sd'] = float(val)
                elif key == 'error rate':
                    result['error_rate'] = float(val)
                elif key == 'bases mapped (cigar)':
                    result['bases_mapped'] = int(val)
    except Exception as e:
        sys.stderr.write(f"[stats_to_legacy] WARN get_bam_stats({stats_file}): {e}\n")
    return result

def calc_breadth_of_coverage(stats_file, genome_size, min_dp=1):
    """
    Calculate breadth of coverage from samtools stats coverage histogram (COV lines).
    Returns fraction of genome covered at >= min_dp depth.
    """
    if genome_size <= 0:
        return 0.0

    # COV lines: COV [range] bases_at_this_depth fraction
    # We need to sum bases with depth >= min_dp
    bases_covered = 0
    try:
        with open(stats_file, 'r') as f:
            for line in f:
                if not line.startswith('COV\t'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 4:
                    continue
                # samtools stats COV format: COV  [depth_range]  depth  count
                # e.g. COV	[1-1]	1	1234   (parts[2]=depth, parts[3]=#bases). The overflow bin is
                # [1000<]; strip the '<'/'>'/'=' so those deep-covered bases are not dropped as ValueError.
                depth_range = parts[1].strip('[]').strip('<>=')
                try:
                    if '-' in depth_range:
                        low = int(depth_range.split('-')[0])
                    else:
                        low = int(depth_range)
                except ValueError:
                    continue

                if low >= min_dp:
                    bases_covered += int(parts[3])
    except Exception as e:
        sys.stderr.write(f"[stats_to_legacy] WARN calc_breadth_of_coverage({stats_file}): {e}\n")

    return bases_covered / genome_size if genome_size > 0 else 0.0

def cov_histogram_stats(stats_file, genome_size):
    """From the samtools-stats COV histogram derive: median depth, breadth at >=5x / >=10x,
    and coverage evenness = coefficient of variation of depth (sd/mean over covered+uncovered bases).
    The COV block bins are [low-high] count; we take the bin midpoint as the depth for that many bases.
    genome_size anchors the zero-depth bases the histogram omits (breadth/CoV need the full denominator)."""
    out = {'median_depth': 0.0, 'breadth5': 0.0, 'breadth10': 0.0, 'evenness': 0.0}
    if genome_size <= 0:
        return out
    bins = []  # (depth_midpoint, n_bases)
    try:
        with open(stats_file, 'r') as f:
            for line in f:
                if not line.startswith('COV\t'):
                    continue
                p = line.rstrip('\n').split('\t')
                if len(p) < 4:
                    continue
                rng = p[1].strip('[]').strip('<>=')   # '[1000<]' overflow bin -> '1000' (not a ValueError)
                try:
                    if '-' in rng:
                        lo, hi = rng.split('-', 1)
                        d = (int(lo) + int(hi)) / 2.0
                    else:
                        d = float(rng)
                    n = int(p[3])   # samtools COV: [range] depth count -> count is col 4
                except ValueError:
                    continue
                bins.append((d, n))
    except Exception as e:
        sys.stderr.write(f"[stats_to_legacy] WARN cov_histogram_stats({stats_file}): {e}\n")
        return out
    covered = sum(n for _, n in bins)
    zero = max(0, genome_size - covered)         # bases at depth 0 the COV block omits
    total = covered + zero
    if total <= 0:
        return out
    out['breadth5'] = 100.0 * sum(n for d, n in bins if d >= 5) / total
    out['breadth10'] = 100.0 * sum(n for d, n in bins if d >= 10) / total
    # median depth over ALL bases (including zero-depth) via a running cumulative count
    half = total / 2.0
    run = zero
    med = 0.0
    if run < half:
        for d, n in sorted(bins):
            run += n
            if run >= half:
                med = d
                break
    out['median_depth'] = med
    # coverage evenness = CoV = sd/mean over the full genome (zero bases included)
    ssum = sum(d * n for d, n in bins)
    mean = ssum / total
    if mean > 0:
        var = (sum((d * d) * n for d, n in bins) + 0.0) / total - mean * mean
        var = max(0.0, var)
        out['evenness'] = (var ** 0.5) / mean
    return out

BIN_N = 200  # bins along the reference for the per-position variant-density profiles (must match qc_report.NBINS)

def analyze_vcf(vcf_file, g_size=0, nbins=BIN_N):
    stats = {'snps': 0, 'indels': 0, 'homo_total': 0, 'het_total': 0, 'homo_indels': 0,
             'ti': 0, 'tv': 0,
             'snp_prof': [0] * nbins, 'het_prof': [0] * nbins, 'indel_prof': [0] * nbins,
             'ann_present': False, 'ann_records': 0, 'ann_total': 0,
             'ann_impact': {k: 0 for k in _IMPACTS},
             'ann_class': {c: 0 for c in _EFF_CLASSES},
             'ann_warn': 0, 'ann_db_error': False, 'ann_gene': {}}
    opener = gzip.open if vcf_file.endswith('.gz') else open
    try:
        with opener(vcf_file, 'rt') as f:
            for line in f:
                if line.startswith('#'): continue
                parts = line.strip().split('\t')
                if len(parts) < 10: continue
                ref = parts[3]
                alt = parts[4]
                if alt in ['.', '*']: continue
                fmt = parts[8].split(':')
                if 'GT' not in fmt: continue
                gt_idx = fmt.index('GT')
                sample = parts[9].split(':')
                if len(sample) <= gt_idx: continue
                gt = sample[gt_idx]
                if gt in ['./.', '0/0', '0|0']: continue

                # snpEff functional annotation of this real variant (no-op when the VCF has no ANN=)
                info_str = parts[7] if len(parts) > 7 else ''
                stats['ann_total'] += 1
                ann = _parse_ann_field(info_str)
                if ann is not None:
                    cls, impact, gene, has_error, db_error = ann
                    stats['ann_present'] = True
                    stats['ann_records'] += 1
                    stats['ann_impact'][impact] += 1
                    stats['ann_class'][cls] += 1
                    if has_error:
                        stats['ann_warn'] += 1
                    if db_error:
                        stats['ann_db_error'] = True
                    if gene:
                        g = stats['ann_gene'].get(gene)
                        if g is None:
                            g = [0, 0, {}]
                            stats['ann_gene'][gene] = g
                        if impact == 'HIGH':
                            g[0] += 1
                        elif impact == 'MODERATE':
                            g[1] += 1
                        g[2][cls] = g[2].get(cls, 0) + 1

                is_snp = (len(ref) == 1 and len(alt) == 1)
                is_homo = False
                alleles = gt.replace('|', '/').split('/')
                if len(alleles) == 2 and alleles[0] == alleles[1]:
                    is_homo = True

                # positional bin (POS over the reference) for the genome-landscape variant tracks
                b = None
                if g_size > 0:
                    try:
                        b = min(nbins - 1, (int(parts[1]) * nbins) // g_size)
                    except ValueError:
                        b = None

                if is_snp:
                    stats['snps'] += 1
                    r = ref.upper(); a = alt.upper()
                    if r in 'ACGT' and a in 'ACGT' and r != a:
                        if frozenset((r, a)) in _TRANSITIONS:
                            stats['ti'] += 1
                        else:
                            stats['tv'] += 1
                    if b is not None and is_homo:
                        stats['snp_prof'][b] += 1
                else:
                    stats['indels'] += 1
                    if is_homo: stats['homo_indels'] += 1
                    if b is not None:
                        stats['indel_prof'][b] += 1
                if is_homo:
                    stats['homo_total'] += 1
                else:
                    stats['het_total'] += 1
                    if b is not None:
                        stats['het_prof'][b] += 1
    except Exception as e:
        sys.stderr.write(f"[stats_to_legacy] WARN analyze_vcf({vcf_file}): {e} (Ti/Tv & SNP counts may be partial)\n")
    return stats

def get_pathotypr_data(report_file):
    """
    Parses the Pathotypr summary TSV.
    Header usually: genome <tab> lineage:count <tab> major_lineage
    Returns a dict with 'major' and 'counts'.
    """
    results = {'major': 'NA', 'counts': 'NA'}

    if not report_file or not os.path.exists(report_file) or "NO_LINEAGE_FILE" in report_file:
        return results

    try:
        with open(report_file, 'r') as f:
            # Read header
            header_line = f.readline().strip()
            if not header_line: return results
            headers = header_line.split('\t')

            # Identify columns
            try:
                idx_major = headers.index('major_lineage')
            except ValueError:
                idx_major = -1

            try:
                idx_count = headers.index('lineage:count')
            except ValueError:
                idx_count = -1

            # Read data line
            data_line = f.readline().strip()
            if data_line:
                cols = data_line.split('\t')

                # Extract Major Lineage
                if idx_major != -1 and len(cols) > idx_major:
                    results['major'] = cols[idx_major]

                # Extract Lineage Counts
                if idx_count != -1 and len(cols) > idx_count:
                    results['counts'] = cols[idx_count]

    except Exception as e:
        sys.stderr.write(f"Warning: Error parsing Pathotypr file: {e}\n")

    return results

def main():
    args = parse_args()
    data = {k: "NA" for k in ORDERED_KEYS}

    # 1. Metadata & Fastp
    data['DAT_OUT'] = str(datetime.date.today())
    data['GID_IN'] = args.sample
    data['FQ_IN'] = f"{args.sample}_R1/R2"
    data['ANC_OUT'] = 'yes' if str(args.ancient).strip().lower() in ('yes', 'y', 'true', '1', 'ancient') else 'no'

    g_size = get_genome_size(args.ref_fai)
    if g_size > 0:
        data['GENOME_LEN'] = str(g_size)   # reference length -> SNP density works even when make_consensus is off
    r_len = 0.0

    try:
        with open(args.fastp_json, 'r') as f:
            j = json.load(f)
            before = j.get('summary', {}).get('before_filtering', {})
            after = j.get('summary', {}).get('after_filtering', {})
            dup = j.get('duplication', {})

            data['TRD_OUT'] = str(before.get('total_reads', 'NA'))
            data['TRM_OUT'] = str(after.get('total_reads', 'NA'))
            r_len = float(after.get('read1_mean_length', 0))
            data['RDL_OUT'] = f"{r_len:.1f}"
            data['DUP_OUT'] = f"{float(dup.get('rate', 0))*100:.2f}"
            if 'q20_rate' in after: data['Q20_PCT'] = f"{float(after['q20_rate'])*100:.2f}"
            if 'q30_rate' in after: data['Q30_PCT'] = f"{float(after['q30_rate'])*100:.2f}"
            if 'gc_content' in after: data['GC_PCT'] = f"{float(after['gc_content'])*100:.2f}"
            data['PHE_OUT'] = '33'
            data['FV_OUT'] = 'FASTQ_SUCCESS'

            if g_size > 0:
                data['LNG_OUT'] = f"{float(after.get('total_bases', 0))/g_size:.2f}"
    except Exception as e:
        sys.stderr.write(f"[stats_to_legacy] WARN fastp parse ({args.fastp_json}): {e}\n")

    # 2. Mapping
    bam = get_bam_stats(args.bam_stats)
    data['MAP_OUT'] = str(bam['mapped'])

    # Mapping percentage
    if bam['total'] > 0:
        data['MAP_PCT'] = f"{(bam['mapped'] / bam['total']) * 100:.2f}"

    # Properly paired percentage
    if bam['total'] > 0:
        data['PP_OUT'] = f"{(bam['properly_paired'] / bam['total']) * 100:.2f}"

    # Average base quality
    if bam['avg_quality'] > 0:
        data['STD_OUT'] = f"{bam['avg_quality']:.1f}"

    # Coverage: depth from CIGAR-based mapped bases (more accurate than reads * read_length)
    if g_size > 0:
        if bam['bases_mapped'] > 0:
            data['COV_OUT'] = f"{bam['bases_mapped'] / g_size:.2f}"
        elif r_len > 0 and bam['mapped'] > 0:
            data['COV_OUT'] = f"{(bam['mapped'] * r_len) / g_size:.2f}"

    # Breadth of coverage (% genome covered at >= 1x)
    breadth = calc_breadth_of_coverage(args.bam_stats, g_size, min_dp=1)
    if breadth > 0:
        data['COV_BREADTH'] = f"{breadth * 100:.2f}"

    # Extended coverage metrics from the same COV histogram (no extra tool)
    ch = cov_histogram_stats(args.bam_stats, g_size)
    if g_size > 0:
        if ch['median_depth'] > 0:
            data['COV_MEDIAN'] = f"{ch['median_depth']:.2f}"
        data['COV_BREADTH5'] = f"{ch['breadth5']:.2f}"
        data['COV_BREADTH10'] = f"{ch['breadth10']:.2f}"
        if ch['evenness'] > 0:
            data['COV_EVENNESS'] = f"{ch['evenness']:.3f}"

    # Extended mapping metrics from the samtools stats SN block
    if bam['total'] > 0:
        data['UNMAP_PCT'] = f"{(bam['reads_unmapped'] / bam['total']) * 100:.2f}"
        data['MQ0_PCT'] = f"{(bam['reads_mq0'] / bam['total']) * 100:.2f}"
        singletons = max(0, bam['mapped'] - bam['mapped_paired'])   # approx (no flagstat passed)
        data['SINGLE_PCT'] = f"{(singletons / bam['total']) * 100:.2f}"
    if bam['insert_size_avg'] > 0:
        data['ISIZE_MEAN'] = f"{bam['insert_size_avg']:.1f}"
    if bam['insert_size_sd'] > 0:
        data['ISIZE_SD'] = f"{bam['insert_size_sd']:.1f}"
    if bam['error_rate'] > 0:
        data['ERR_RATE'] = f"{bam['error_rate']*100:.3f}"   # emit as a percent

    # 3. Variants
    v = analyze_vcf(args.vcf, g_size, BIN_N)
    if g_size > 0:
        data['SNP_PROFILE'] = ','.join(map(str, v['snp_prof']))       # homozygous SNP counts per reference bin
        data['HET_PROFILE'] = ','.join(map(str, v['het_prof']))       # heterozygous variant counts per bin
        data['INDEL_PROFILE'] = ','.join(map(str, v['indel_prof']))   # indel counts per bin
    data['SNP_OUT'] = str(v['snps'])
    data['HOM_OUT'] = str(v['homo_total'])
    data['HET_OUT'] = str(v['het_total'])
    data['VAR_OUT'] = str(v['snps'] + v['indels'])
    data['HOMO_IND_OUT'] = str(v['homo_indels'])
    # Ti/Tv ratio: an aDNA damage / mapping-noise sentinel. Deamination pushes it up (transition-heavy);
    # random mapping/base-call noise pushes it down (transversion-heavy, as on the excluded long-branch mummy).
    data['TITV_OUT'] = f"{(v['ti'] / v['tv']):.3f}" if v['tv'] > 0 else "NA"

    # 3b. Functional annotation (snpEff ANN); an un-annotated VCF leaves every ANN_* at the default "NA"
    if v['ann_present']:
        imp = v['ann_impact']; cl = v['ann_class']
        denom = v['ann_total']
        data['ANN_PRESENT'] = 'yes'
        data['ANN_TOTAL'] = str(denom)
        data['ANN_ANNOTATED'] = str(v['ann_records'])
        data['ANN_ANNOTATED_FRAC'] = f"{v['ann_records']/denom:.4f}" if denom > 0 else "NA"
        data['ANN_IMPACT_HIGH']     = str(imp['HIGH'])
        data['ANN_IMPACT_MODERATE'] = str(imp['MODERATE'])
        data['ANN_IMPACT_LOW']      = str(imp['LOW'])
        data['ANN_IMPACT_MODIFIER'] = str(imp['MODIFIER'])
        data['ANN_MISSENSE']      = str(cl['missense'])
        data['ANN_SYNON']         = str(cl['synonymous'])
        data['ANN_STOP_GAINED']   = str(cl['stop_gained'])
        data['ANN_STOP_LOST']     = str(cl['stop_lost'])
        data['ANN_START_LOST']    = str(cl['start_lost'])
        data['ANN_FRAMESHIFT']    = str(cl['frameshift'])
        data['ANN_INFRAME_INDEL'] = str(cl['inframe_indel'])
        data['ANN_SPLICE']        = str(cl['splice'])
        data['ANN_INTERGENIC']    = str(cl['intergenic'])
        data['ANN_REGULATORY']    = str(cl['regulatory'])
        # loss-of-function = stop_gained + frameshift + start_lost + stop_lost + splice (single authoritative def)
        data['ANN_LOF'] = str(cl['stop_gained'] + cl['frameshift'] + cl['start_lost'] + cl['stop_lost'] + cl['splice'])
        # CODING is the superset of LoF (splice included) so that lof_pct = LOF/CODING stays in [0, 100]
        data['ANN_CODING'] = str(cl['missense'] + cl['synonymous'] + cl['stop_gained'] + cl['stop_lost'] +
                                 cl['start_lost'] + cl['frameshift'] + cl['inframe_indel'] + cl['splice'])
        syn = cl['synonymous']
        # missense:silent VARIANT-count ratio = a per-sample pN/pS PROXY (not site-corrected, not dN/dS)
        data['ANN_MISSENSE_SILENT'] = f"{cl['missense']/syn:.3f}" if syn > 0 else "NA"
        data['ANN_WARN'] = str(v['ann_warn'])
        data['ANN_DB_ERROR'] = 'yes' if v['ann_db_error'] else 'no'
        data['ANN_GENE_BURDEN'] = _format_gene_burden(v['ann_gene'])

    # 4. Pathotypr Lineage Data
    if args.pathotypr_report:
        p_data = get_pathotypr_data(args.pathotypr_report)
        data['LIN_OUT'] = p_data['major']
        data['LIN_COUNTS'] = p_data['counts']

    # 5. Pathotypr DR Data
    if args.pathotypr_dr_report:
        dr_data = get_pathotypr_data(args.pathotypr_dr_report)
        data['DR_OUT'] = dr_data['major']
        data['DR_COUNTS'] = dr_data['counts']

    # Write legacy log
    with open(f"{args.sample}.log", 'w') as out:
        for key in ORDERED_KEYS:
            out.write(f"{key}\t{data[key]}\n")

if __name__ == "__main__":
    main()
