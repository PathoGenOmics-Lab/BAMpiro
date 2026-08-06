#!/usr/bin/env python3
"""Master SNP matrix.

Aggregate the per-sample annotated VCFs of a run into a single wide TSV:
  rows    = every SNP site seen in any sample (sorted by contig, position)
  columns = the reference / annotation, then TWO columns per sample: allele
            frequency (|AF) and depth (|DP)

Allele frequency comes from AD (or AO+RO) when present, else the GT dosage
(hom-alt 1.0 / het 0.5). Depth comes from the FORMAT DP field, else AO+RO / sum(AD).
A sample with no call at a site is left blank (not called in that sample).

Only SNPs (single-base REF and ALT) are kept. The '--reference' value labels the
reference the samples were mapped against.
"""
from __future__ import annotations
import argparse
import gzip
import os
import sys


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

    header = ['reference', 'contig', 'pos', 'ref_allele', 'alt_allele', 'gene', 'effect', 'aa_change']
    for s in samples:
        header += [f'{s}|AF', f'{s}|DP']
    with open(args.output, 'w') as out:
        out.write('\t'.join(header) + '\n')
        for contig, pos in sorted(sites):
            st = sites[(contig, pos)]
            row = [args.reference or contig, contig, str(pos), st['ref'],
                   ','.join(sorted(st['alt'])), st['gene'], st['eff'], st['aa']]
            for s in samples:
                if s in st['cells']:
                    af, dp = st['cells'][s]
                    row += [f'{af:.4f}', ('' if dp is None else str(dp))]
                else:
                    row += ['', '']   # not called in this sample
            out.write('\t'.join(row) + '\n')
    sys.stderr.write(f"[snp_matrix] {len(sites)} SNP site(s) x {len(samples)} sample(s) -> {args.output}\n")


if __name__ == '__main__':
    main()
