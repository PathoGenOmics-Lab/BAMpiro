#!/usr/bin/env python3
"""
Consolidates individual sample .log files into a single summary TSV.
Selects the most useful QC columns for a quick overview table.
"""
import sys
import os
import argparse

# Columns to include in the summary (human-readable names)
SUMMARY_COLUMNS = {
    'GID_IN':       'sample_id',
    'TRD_OUT':      'raw_reads',
    'TRM_OUT':      'trimmed_reads',
    'RDL_OUT':      'read_length',
    'DUP_OUT':      'duplication_pct',
    'MAP_OUT':      'mapped_reads',
    'MAP_PCT':      'mapped_pct',
    'PP_OUT':       'properly_paired_pct',
    'COV_OUT':      'mean_depth',
    'COV_MEDIAN':   'median_depth',
    'COV_BREADTH':  'breadth_pct',
    'COV_BREADTH5': 'breadth5x_pct',
    'COV_BREADTH10':'breadth10x_pct',
    'COV_EVENNESS': 'coverage_cv',
    'GENOME_LEN':   'genome_len',
    'UNMAP_PCT':    'unmapped_pct',
    'SINGLE_PCT':   'singleton_pct',
    'MQ0_PCT':      'mapq0_pct',
    'ERR_RATE':     'error_rate_pct',
    'ISIZE_MEAN':   'insert_size',
    'ISIZE_SD':     'insert_size_sd',
    'Q20_PCT':      'q20_pct',
    'Q30_PCT':      'q30_pct',
    'GC_PCT':       'gc_pct',
    'STD_OUT':      'avg_base_quality',
    'LNG_OUT':      'est_genome_cov',
    'VAR_OUT':      'total_variants',
    'SNP_OUT':      'snps',
    'HOM_OUT':      'homo_variants',
    'HET_OUT':      'het_variants',
    'HOMO_IND_OUT': 'homo_indels',
    'TITV_OUT':     'ti_tv',
    'SNP_PROFILE':  'snp_profile',
    'HET_PROFILE':  'het_profile',
    'INDEL_PROFILE':'indel_profile',
    'ANN_PRESENT':          'ann_present',
    'ANN_ANNOTATED_FRAC':   'ann_annotated_frac',
    'ANN_TOTAL':            'ann_total',
    'ANN_ANNOTATED':        'ann_annotated',
    'ANN_IMPACT_HIGH':      'ann_high',
    'ANN_IMPACT_MODERATE':  'ann_moderate',
    'ANN_IMPACT_LOW':       'ann_low',
    'ANN_IMPACT_MODIFIER':  'ann_modifier',
    'ANN_MISSENSE':         'eff_missense',
    'ANN_SYNON':            'eff_synonymous',
    'ANN_STOP_GAINED':      'eff_stop_gained',
    'ANN_STOP_LOST':        'eff_stop_lost',
    'ANN_START_LOST':       'eff_start_lost',
    'ANN_FRAMESHIFT':       'eff_frameshift',
    'ANN_INFRAME_INDEL':    'eff_inframe_indel',
    'ANN_SPLICE':           'eff_splice',
    'ANN_INTERGENIC':       'eff_intergenic',
    'ANN_REGULATORY':       'eff_regulatory',
    'ANN_LOF':              'lof_count',
    'ANN_CODING':           'coding_count',
    'ANN_MISSENSE_SILENT':  'missense_silent',
    'ANN_WARN':             'snpeff_warn',
    'ANN_DB_ERROR':         'ann_db_error',
    'LIN_OUT':      'lineage',
    'LIN_COUNTS':   'lineage_counts',
    'DR_OUT':       'drug_resistance',
    'DR_COUNTS':    'dr_counts',
    'DAT_OUT':      'date',
    'ANC_OUT':      'ancient',
}

# Order for the output table
COLUMN_ORDER = [
    'GID_IN', 'TRD_OUT', 'TRM_OUT', 'RDL_OUT', 'Q20_PCT', 'Q30_PCT', 'GC_PCT', 'DUP_OUT',
    'MAP_OUT', 'MAP_PCT', 'PP_OUT', 'UNMAP_PCT', 'SINGLE_PCT', 'MQ0_PCT', 'ERR_RATE',
    'ISIZE_MEAN', 'ISIZE_SD',
    'COV_OUT', 'COV_MEDIAN', 'COV_BREADTH', 'COV_BREADTH5', 'COV_BREADTH10', 'COV_EVENNESS', 'GENOME_LEN', 'STD_OUT',
    'LNG_OUT', 'VAR_OUT', 'SNP_OUT', 'HOM_OUT', 'HET_OUT', 'HOMO_IND_OUT', 'TITV_OUT',
    'ANN_PRESENT', 'ANN_ANNOTATED_FRAC', 'ANN_TOTAL', 'ANN_ANNOTATED',
    'ANN_IMPACT_HIGH', 'ANN_IMPACT_MODERATE', 'ANN_IMPACT_LOW', 'ANN_IMPACT_MODIFIER',
    'ANN_MISSENSE', 'ANN_SYNON', 'ANN_STOP_GAINED', 'ANN_STOP_LOST', 'ANN_START_LOST',
    'ANN_FRAMESHIFT', 'ANN_INFRAME_INDEL', 'ANN_SPLICE', 'ANN_INTERGENIC', 'ANN_REGULATORY',
    'ANN_LOF', 'ANN_CODING', 'ANN_MISSENSE_SILENT', 'ANN_WARN', 'ANN_DB_ERROR',
    'SNP_PROFILE', 'HET_PROFILE', 'INDEL_PROFILE',
    'LIN_OUT', 'LIN_COUNTS', 'DR_OUT', 'DR_COUNTS', 'DAT_OUT', 'ANC_OUT',
]

# Percentage columns the report reads directly; DERIVED here (once, in Python) so qc_flags.tsv and the HTML agree.
DERIVED_COLUMNS = ['annotated_pct', 'coding_pct', 'lof_pct']

def _to_int(v):
    try:
        return int(v)
    except (TypeError, ValueError):
        return None

def _derive_ann(data):
    """(annotated_pct, coding_pct, lof_pct) as strings, or ('NA','NA','NA') for an un-annotated sample."""
    if data.get('ANN_PRESENT') != 'yes':
        return ('NA', 'NA', 'NA')
    tot  = _to_int(data.get('ANN_TOTAL'))
    ann  = _to_int(data.get('ANN_ANNOTATED'))
    lofc = _to_int(data.get('ANN_LOF'))
    codc = _to_int(data.get('ANN_CODING'))
    hi   = _to_int(data.get('ANN_IMPACT_HIGH'))
    mod  = _to_int(data.get('ANN_IMPACT_MODERATE'))
    low  = _to_int(data.get('ANN_IMPACT_LOW'))
    modf = _to_int(data.get('ANN_IMPACT_MODIFIER'))
    ann_pct = f"{100.0*ann/tot:.2f}" if (ann is not None and tot) else 'NA'
    # coding % = variants in a transcript (HIGH+MODERATE+LOW) / all impact-annotated (MODIFIER = intergenic/UTR bucket)
    denom_imp = sum(x for x in (hi, mod, low, modf) if x is not None)
    coding_pct = f"{100.0*((hi or 0)+(mod or 0)+(low or 0))/denom_imp:.2f}" if denom_imp else 'NA'
    lof_pct = f"{100.0*lofc/codc:.2f}" if (lofc is not None and codc) else 'NA'
    return (ann_pct, coding_pct, lof_pct)

def _accumulate_gene_burden(agg, gene_burden_str):
    """agg: {gene: {'high':int, 'mod':int, 'n_samples':int, 'eff':{cls:count}}}. Reads the compact
    ANN_GENE_BURDEN string ('gene:HIGH:MOD:topEff;...') from one sample. O(<=40) per sample."""
    if not gene_burden_str or gene_burden_str == 'NA':
        return
    for tok in gene_burden_str.split(';'):
        parts = tok.split(':')
        if len(parts) != 4:
            continue
        gene, high_s, mod_s, top_eff = parts
        try:
            high = int(high_s); mod = int(mod_s)
        except ValueError:
            continue
        g = agg.get(gene)
        if g is None:
            g = {'high': 0, 'mod': 0, 'n_samples': 0, 'eff': {}}
            agg[gene] = g
        g['high'] += high
        g['mod'] += mod
        g['n_samples'] += 1
        g['eff'][top_eff] = g['eff'].get(top_eff, 0) + 1


def _write_gene_burden(agg, path):
    """Cohort gene-burden aux TSV, most-burdened first. low/modifier/start/end are emitted empty
    (the compact per-sample key carries only HIGH/MOD counts, no coordinates)."""
    header = ['gene', 'high', 'moderate', 'low', 'modifier',
              'dominant_effect', 'n_samples', 'total_impactful', 'start', 'end']
    rows = []
    for gene, g in agg.items():
        dom = max(g['eff'].items(), key=lambda kv: kv[1])[0] if g['eff'] else 'NA'
        rows.append((gene, g['high'], g['mod'], '', '', dom, g['n_samples'], g['high'] + g['mod'], '', ''))
    rows.sort(key=lambda r: (-r[1], -r[2], r[0]))
    with open(path, 'w') as out:
        out.write('\t'.join(header) + '\n')
        for r in rows:
            out.write('\t'.join(str(x) for x in r) + '\n')


def parse_log(log_file):
    """Parse a single sample .log file into a dict."""
    data = {}
    try:
        with open(log_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or '\t' not in line:
                    continue
                key, val = line.split('\t', 1)
                data[key] = val
    except Exception as e:
        sys.stderr.write(f"Warning: Could not parse {log_file}: {e}\n")
    return data

def main():
    parser = argparse.ArgumentParser(description="Collect sample logs into summary TSV")
    parser.add_argument('logs', nargs='+', help="Sample .log files")
    parser.add_argument('-o', '--output', default='pipeline_summary.tsv', help="Output TSV file")
    parser.add_argument('--gene-burden-out', default=None,
                        help="Optional: also write a cohort per-gene functional-burden TSV here")
    args = parser.parse_args()

    # Write header
    header = [SUMMARY_COLUMNS[k] for k in COLUMN_ORDER] + DERIVED_COLUMNS

    rows = []
    gene_agg = {}
    for log_file in sorted(args.logs):
        if not os.path.exists(log_file):
            continue
        data = parse_log(log_file)
        if not data:
            continue
        row = [data.get(k, 'NA') for k in COLUMN_ORDER]
        row.extend(_derive_ann(data))          # annotated_pct, coding_pct, lof_pct
        rows.append(row)
        _accumulate_gene_burden(gene_agg, data.get('ANN_GENE_BURDEN', 'NA'))

    with open(args.output, 'w') as out:
        out.write('\t'.join(header) + '\n')
        for row in rows:
            out.write('\t'.join(row) + '\n')

    if args.gene_burden_out and gene_agg:
        _write_gene_burden(gene_agg, args.gene_burden_out)
        sys.stderr.write(f"Gene burden: {len(gene_agg)} genes written to {args.gene_burden_out}\n")

    sys.stderr.write(f"Summary: {len(rows)} samples written to {args.output}\n")

if __name__ == "__main__":
    main()
