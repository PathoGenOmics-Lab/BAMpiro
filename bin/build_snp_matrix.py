#!/usr/bin/env python3
"""Master SNP matrix.

Aggregate the per-sample annotated VCFs of a run into a single wide TSV:
  rows    = every SNP site seen in any sample (sorted by contig, position)
  columns = the reference / annotation, then TWO columns per sample: allele
            frequency (|AF) and depth (|DP)

Allele frequency comes from AD (or AO+RO) when present, else the GT dosage
(hom-alt 1.0 / het 0.5). Depth comes from the FORMAT DP field, else AO+RO / sum(AD).

A sample with no call at a site is not the same as a sample nobody read there, and a blank
cell said both. With the per-sample all-positions VCFs (--depth-vcfs), which hold a record for
every position of the sample's reference, the cell says which it was:

  AF 0 and the site's depth   the sample was read there and carries no alternate allele
  AF blank and DP 0           no read covers the site
  both blank                  the file has no record there: a site of another reference, in a
                              run mapped against several, or no all-positions VCF for the sample

That depth is the pileup the consensus is built from, which does not apply the caller's
mapping-quality floor, so in a repeat it can read higher than the depth of a called cell.

The same files also carry the calls the variant VCFs spell as part of a longer allele (an MNP
or a complex record, which this matrix skips): atomised, they are single-base SNPs at the site,
and the cell takes their allele fraction rather than a 0 that would deny them.

Only SNPs (single-base REF and ALT) are kept. The '--reference' value labels the
reference the samples were mapped against.
"""
from __future__ import annotations
import argparse
import gzip
import os
import sys
import zlib
from array import array
from concurrent.futures import ProcessPoolExecutor


def _open(p):
    return gzip.open(p, 'rt', encoding='utf-8', errors='replace') if str(p).endswith('.gz') \
        else open(p, encoding='utf-8', errors='replace')


def _af_dp(fmt, val):
    """(allele_fraction_of_first_alt, total_depth) from one FORMAT/sample pair."""
    d = dict(zip(fmt.split(':'), val.split(':')))
    dp = None
    if d.get('DP', '.') not in ('.', ''):
        try:
            dp = int(d['DP'])
        except ValueError:
            dp = None
    ad = None
    if 'AD' in d:
        try:
            ad = [int(x) for x in d['AD'].split(',') if x not in ('.', '')]
        except ValueError:
            ad = None
    ao = ro = None
    if 'AO' in d and 'RO' in d:
        try:
            ao = sum(int(x) for x in d['AO'].split(',') if x not in ('.', ''))
            ro = int(d['RO']) if d['RO'] not in ('.', '') else 0
        except ValueError:
            ao = ro = None
    af = None
    if ad and len(ad) >= 2 and sum(ad) > 0:
        af = sum(ad[1:]) / sum(ad)
        if dp is None:
            dp = sum(ad)
    elif ao is not None and ro is not None and (ao + ro) > 0:
        af = ao / (ao + ro)
        if dp is None:
            dp = ao + ro
    else:
        alleles = [a for a in d.get('GT', './.').replace('|', '/').split('/') if a not in ('.', '')]
        if alleles:
            af = sum(1 for a in alleles if a != '0') / len(alleles)
    # A single-element AD cannot give an allele fraction, but it still states a depth: without this
    # the matrix emits an empty DP cell for a site whose depth is right there in the record.
    if dp is None and ad:
        dp = sum(ad)
    return af, dp


# Read once per worker process from _init_worker, so the site index crosses to each worker once
# rather than once per file.
_INDEX = {}
_NSITES = 0


def _init_worker(index, n_sites):
    global _INDEX, _NSITES
    _INDEX, _NSITES = index, n_sites


def _site_index(keys):
    """{contig: {position as written in a VCF: row}} for the matrix's sorted sites.

    Keyed on the position's text so that matching a four-million-line all-positions VCF costs a
    split and two dictionary lookups per line, with no number parsed for the lines it skips.
    """
    index = {}
    for i, (contig, pos) in enumerate(keys):
        index.setdefault(contig, {})[str(pos)] = i
    return index


def _read_depths(path, index=None, n_sites=None):
    """(sample, depths, calls) from one all-positions VCF, or None when it cannot be read.

    `depths` holds one entry per matrix row: the sample's depth at the site, or -1 where the
    file has no record for it. `calls` maps a row to (allele, AF, DP) where the record there is
    a variant rather than the reference backbone. A file that fails part-way yields nothing:
    half a sample's depths would read as the other half being uncovered.
    """
    index = _INDEX if index is None else index
    n_sites = _NSITES if n_sites is None else n_sites
    depths = array('i', [-1]) * n_sites
    calls = {}
    sample = None
    try:
        with _open(path) as fh:
            for line in fh:
                if line[0] == '#':
                    if line.startswith('#CHROM'):
                        cols = line.rstrip('\n').split('\t')
                        sample = cols[9] if len(cols) > 9 else None
                    continue
                parts = line.split('\t', 2)
                if len(parts) < 3:
                    continue
                rows = index.get(parts[0])
                i = rows.get(parts[1]) if rows is not None else None
                if i is None:
                    continue
                c = line.rstrip('\n').split('\t')
                if len(c) < 10:
                    continue
                af, dp = _af_dp(c[8], c[9])
                if dp is None:
                    continue
                depths[i] = dp
                alt = c[4].split(',')[0]
                if alt not in ('.', '') and len(alt) == 1 and len(c[3]) == 1 and af is not None:
                    calls[i] = (alt, af, dp)
    except (OSError, EOFError, zlib.error) as e:
        sys.stderr.write(f"[snp_matrix] WARN could not read depths from {path}: {e}\n")
        return None
    if sample is None:
        sys.stderr.write(f"[snp_matrix] WARN {path} names no sample; its depths are not used\n")
        return None
    return sample, depths, calls


def _ann(info):
    """(gene, effect, HGVS.p) from the first snpEff ANN entry."""
    for f in info.split(';'):
        if f.startswith('ANN='):
            a = f[4:].split(',')[0].split('|')
            return (a[3] if len(a) > 3 else ''), (a[1] if len(a) > 1 else ''), (a[10] if len(a) > 10 else '')
    return '', '', ''


def main():
    ap = argparse.ArgumentParser(description="Build a master SNP matrix from per-sample annotated VCFs.")
    ap.add_argument('--vcfs', nargs='+', required=True, help="per-sample (annotated) VCFs")
    ap.add_argument('--reference', default='', help="reference id the samples were mapped against")
    ap.add_argument('--depth-vcfs', nargs='*', default=[],
                    help="per-sample all-positions VCFs (one record per reference position). A sample "
                         "with no call at a site then reads AF 0 with the site's depth instead of a "
                         "blank, and a blank is left only where nothing was read")
    ap.add_argument('--threads', type=int, default=1,
                    help="all-positions VCFs read in parallel")
    ap.add_argument('-o', '--output', required=True)
    args = ap.parse_args()

    samples = []           # preserve first-seen order
    sites = {}             # (contig, pos) -> {ref, alt:set, gene, eff, aa, cells:{sample:(af,dp)}}
    for p in args.vcfs:
        if not p or not os.path.exists(p):
            continue
        try:
            with _open(p) as fh:
                sample = None
                for line in fh:
                    if line.startswith('##'):
                        continue
                    if line.startswith('#CHROM'):
                        cols = line.rstrip('\n').split('\t')
                        sample = cols[9] if len(cols) > 9 else os.path.basename(p).split('.')[0]
                        if sample not in samples:
                            samples.append(sample)
                        continue
                    if sample is None:
                        continue
                    c = line.rstrip('\n').split('\t')
                    if len(c) < 8:
                        continue
                    contig, pos, ref, alt = c[0], c[1], c[3], c[4]
                    if alt in ('.', '') or len(ref) != 1:
                        continue
                    alt1 = alt.split(',')[0]
                    if len(alt1) != 1:
                        continue  # SNPs only
                    af, dp = _af_dp(c[8], c[9]) if len(c) >= 10 else (1.0, None)
                    if af is None:
                        continue
                    gene, eff, aa = _ann(c[7])
                    try:
                        key = (contig, int(pos))
                    except ValueError:
                        continue
                    st = sites.setdefault(key, {'ref': ref, 'alt': set(), 'gene': '', 'eff': '', 'aa': '', 'cells': {}})
                    st['alt'].add(alt1)
                    if gene and not st['gene']:
                        st['gene'] = gene
                    if eff and not st['eff']:
                        st['eff'] = eff
                    if aa and not st['aa']:
                        st['aa'] = aa
                    st['cells'][sample] = (af, dp)
        except OSError as e:
            sys.stderr.write(f"[snp_matrix] WARN could not read {p}: {e}\n")

    keys = sorted(sites)
    depths, n_recovered = {}, 0
    depth_vcfs = [p for p in args.depth_vcfs if p and os.path.exists(p)]
    if depth_vcfs and keys:
        index = _site_index(keys)
        if args.threads > 1 and len(depth_vcfs) > 1:
            with ProcessPoolExecutor(max_workers=args.threads, initializer=_init_worker,
                                     initargs=(index, len(keys))) as pool:
                results = list(pool.map(_read_depths, depth_vcfs))
        else:
            results = [_read_depths(p, index, len(keys)) for p in depth_vcfs]
        for res in results:
            if res is None:
                continue
            sample, d, calls = res
            if sample not in samples:
                continue   # the columns are the samples whose variants were given
            have = depths.get(sample)
            if have is None:
                depths[sample] = have = d
            else:          # one sample mapped against two references: two files, disjoint contigs
                for i, v in enumerate(d):
                    if v >= 0 and have[i] < 0:
                        have[i] = v
            for i, (alt, af, dp) in calls.items():
                st = sites[keys[i]]
                if sample not in st['cells']:
                    st['alt'].add(alt)
                    st['cells'][sample] = (af, dp)
                    n_recovered += 1

    header = ['reference', 'contig', 'pos', 'ref_allele', 'alt_allele', 'gene', 'effect', 'aa_change']
    for s in samples:
        header += [f'{s}|AF', f'{s}|DP']
    n_absent = n_unread = 0
    with open(args.output, 'w') as out:
        out.write('\t'.join(header) + '\n')
        for i, (contig, pos) in enumerate(keys):
            st = sites[(contig, pos)]
            row = [args.reference or contig, contig, str(pos), st['ref'],
                   ','.join(sorted(st['alt'])), st['gene'], st['eff'], st['aa']]
            for s in samples:
                if s in st['cells']:
                    af, dp = st['cells'][s]
                    row += [f'{af:.4f}', ('' if dp is None else str(dp))]
                    continue
                d = depths.get(s)
                dp = d[i] if d is not None else -1
                if dp > 0:
                    row += ['0', str(dp)]    # read there, and no alternate allele was called
                    n_absent += 1
                elif dp == 0:
                    row += ['', '0']         # nobody read it: absence says nothing
                    n_unread += 1
                else:
                    row += ['', '']          # no record: another reference, or no depths given
            out.write('\t'.join(row) + '\n')
    sys.stderr.write(f"[snp_matrix] {len(sites)} SNP site(s) x {len(samples)} sample(s) -> {args.output}\n")
    if args.depth_vcfs and keys:
        missing = [s for s in samples if s not in depths]
        more = f" and {len(missing) - 5} more" if len(missing) > 5 else ""
        sys.stderr.write(
            f"[snp_matrix] depths from {len(depths)} of {len(samples)} sample(s): {n_absent} cell(s) "
            f"read without the allele (AF 0), {n_unread} with no read, {n_recovered} call(s) "
            f"recovered from an MNP or complex record"
            + (f"; no depths for {', '.join(missing[:5])}{more}" if missing else "") + "\n")


if __name__ == '__main__':
    main()
