#!/usr/bin/env python3
"""Master indel matrix.

The indels of a run in one wide TSV, as the SNP matrix holds its SNPs:
  rows    = every indel (contig, position, REF, ALT) that at least one sample calls with a PASS:
            the hom or het rule of bin/vcf_filter_rules.py, the depth, strand support and allele
            fraction a SNP needs
  columns = the reference and the annotation (the indel's length, gene, effect, HGVS.c, HGVS.p),
            how many samples call it with a PASS, then THREE columns per sample: allele fraction
            (|AF), depth (|DP) and whether that sample's own call passes (|FT: PASS or LowSupport)

It reads the per-sample indel VCFs of CALL_FREEBAYES (every normalised indel call, soft-filtered),
annotated by snpEff when the run annotates. A sample that calls a row's indel below the rule keeps
its fraction, marked LowSupport: a minority indel is where resistance through a frameshift often
starts, and a blank would say the sample does not have it.

A sample with no call at all reads as in the SNP matrix, from its all-positions VCF (--depth-vcfs):
AF 0 with the depth at the indel's anchor base where it was read, AF blank and DP 0 where no read
covers it, both blank where the file has no record there (a site of another reference). That depth
counts the reads that show a base at the anchor, so for a long indel it is a floor: a read that
ends inside the event is not counted.
"""
from __future__ import annotations

import argparse
import os
import sys
from concurrent.futures import ProcessPoolExecutor

from build_snp_matrix import _af_dp, _init_worker, _open, _read_depths, _site_index

PASSING = ("PASS", ".", "")        # an unfiltered VCF (no soft filter written) reads as passing


def _ann(info):
    """(gene, effect, HGVS.c, HGVS.p) from the first snpEff ANN entry."""
    for f in info.split(';'):
        if f.startswith('ANN='):
            a = f[4:].split(',')[0].split('|')
            get = lambda k: a[k] if len(a) > k else ''
            return get(3), get(1), get(9), get(10)
    return '', '', '', ''


def read_indels(paths):
    """(samples in first-seen order, {(contig, pos, ref, alt): site}) from per-sample VCFs."""
    samples, sites = [], {}
    for p in paths:
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
                    contig, ref, alt = c[0], c[3], c[4].split(',')[0]
                    if alt in ('.', '', '*') or alt.startswith('<') or len(ref) == len(alt):
                        continue        # SNPs and MNPs are the SNP matrix's
                    try:
                        pos = int(c[1])
                    except ValueError:
                        continue
                    af, dp = _af_dp(c[8], c[9]) if len(c) >= 10 else (1.0, None)
                    if af is None:
                        continue
                    st = sites.setdefault((contig, pos, ref, alt), {'gene': '', 'eff': '', 'hgvs_c': '',
                                                                    'hgvs_p': '', 'cells': {}})
                    gene, eff, hc, hp = _ann(c[7])
                    for k, v in (('gene', gene), ('eff', eff), ('hgvs_c', hc), ('hgvs_p', hp)):
                        if v and not st[k]:
                            st[k] = v
                    if sample not in st['cells']:
                        st['cells'][sample] = (af, dp, c[6] if len(c) > 6 else '')
        except OSError as e:
            sys.stderr.write(f"[indel_matrix] WARN could not read {p}: {e}\n")
    return samples, sites


def main(argv=None):
    ap = argparse.ArgumentParser(description="Build a master indel matrix from per-sample indel VCFs.")
    ap.add_argument('--vcfs', nargs='+', required=True, help="per-sample indel VCFs (annotated or not)")
    ap.add_argument('--reference', default='', help="reference id the samples were mapped against")
    ap.add_argument('--depth-vcfs', nargs='*', default=[],
                    help="per-sample all-positions VCFs: the depth at the anchor base of an indel a "
                         "sample does not call")
    ap.add_argument('--threads', type=int, default=1, help="all-positions VCFs read in parallel")
    ap.add_argument('-o', '--output', required=True)
    args = ap.parse_args(argv)

    samples, sites = read_indels(args.vcfs)
    keys = sorted(k for k, st in sites.items()
                  if any(ft in PASSING for _, _, ft in st['cells'].values()))

    # the depth at each row's anchor base, for the samples that do not call it
    anchors = sorted({(contig, pos) for contig, pos, _, _ in keys})
    at = {k: i for i, k in enumerate(anchors)}
    depths = {}
    depth_vcfs = [p for p in args.depth_vcfs if p and os.path.exists(p)]
    if depth_vcfs and anchors:
        index = _site_index(anchors)
        if args.threads > 1 and len(depth_vcfs) > 1:
            with ProcessPoolExecutor(max_workers=args.threads, initializer=_init_worker,
                                     initargs=(index, len(anchors))) as pool:
                results = list(pool.map(_read_depths, depth_vcfs))
        else:
            results = [_read_depths(p, index, len(anchors)) for p in depth_vcfs]
        for res in results:
            if res is None:
                continue
            sample, d, _, _ = res
            if sample not in samples:
                continue
            have = depths.get(sample)
            if have is None:
                depths[sample] = d
            else:          # one sample mapped against two references: disjoint contigs
                for i, v in enumerate(d):
                    if v >= 0 and have[i] < 0:
                        have[i] = v

    header = ['reference', 'contig', 'pos', 'ref_allele', 'alt_allele', 'length', 'gene', 'effect',
              'hgvs_c', 'hgvs_p', 'n_pass']
    for s in samples:
        header += [f'{s}|AF', f'{s}|DP', f'{s}|FT']
    with open(args.output, 'w') as out:
        out.write('\t'.join(header) + '\n')
        for contig, pos, ref, alt in keys:
            st = sites[(contig, pos, ref, alt)]
            n_pass = sum(1 for _, _, ft in st['cells'].values() if ft in PASSING)
            row = [args.reference or contig, contig, str(pos), ref, alt, str(len(alt) - len(ref)),
                   st['gene'], st['eff'], st['hgvs_c'], st['hgvs_p'], str(n_pass)]
            i = at[(contig, pos)]
            for s in samples:
                if s in st['cells']:
                    af, dp, ft = st['cells'][s]
                    row += [f'{af:.4f}', '' if dp is None else str(dp), 'PASS' if ft in PASSING else ft]
                    continue
                d = depths.get(s)
                dp = d[i] if d is not None else -1
                if dp > 0:
                    row += ['0', str(dp), '']      # read at the anchor, and no such indel called
                elif dp == 0:
                    row += ['', '0', '']           # nobody read it: absence says nothing
                else:
                    row += ['', '', '']            # no record: another reference, or no depths given
            out.write('\t'.join(row) + '\n')
    n_low = sum(1 for k in keys for _, _, ft in sites[k]['cells'].values() if ft not in PASSING)
    sys.stderr.write(f"[indel_matrix] {len(keys)} indel(s) called with a PASS in at least one of "
                     f"{len(samples)} sample(s) ({len(sites) - len(keys)} more only below the rule, left "
                     f"out; {n_low} cell(s) kept below it) -> {args.output}\n")
    return 0


if __name__ == '__main__':
    sys.exit(main())
