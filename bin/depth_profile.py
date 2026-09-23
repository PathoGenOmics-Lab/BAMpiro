#!/usr/bin/env python3
"""What one sample's reads cover, from its all-positions VCF.

The all-positions VCF has a record for every reference position, so it holds the depth of the
whole genome. This reduces it to three small tables for the cohort report:

  <prefix>.depth_windows.tsv  one row per window: mean depth, the share of positions no read
                              covers, and the share deep enough for the consensus to call a
                              base (depth of --callable-dp or more)
  <prefix>.zero_depth.tsv     every stretch of at least --min-run positions no read covers
  <prefix>.gene_depth.tsv     per gene of the GFF: the share of its positions read at all, the
                              share the consensus can call, its mean depth, and that depth
                              against the sample's genome-wide median (with --gff)

A stretch no read covers is a deletion only if the rest of the cohort reads it; otherwise it is
a part of the reference none of these genomes has, or one no read can be placed on. That
comparison needs the cohort, so it is the report's to make, not this script's.

Each table opens with '#' lines naming the sample, the reference and the sample's genome-wide
depth, so a table can be read without its file name.
"""
from __future__ import annotations

import argparse
import gzip
import os
import sys

import numpy as np

from qcreport.parsers import _gff_attr


def _open(p):
    return gzip.open(p, 'rt', encoding='utf-8', errors='replace') if str(p).endswith('.gz') \
        else open(p, encoding='utf-8', errors='replace')


def read_depths(path):
    """(sample, {contig: depth array}, {contig: read-at-all mask}) from an all-positions VCF.

    Arrays are sized from the ##contig lengths and grown when a record lies beyond them. A
    position without a record is not a position of depth 0, so the mask keeps the two apart.
    """
    lengths, order, depth, known = {}, [], {}, {}
    sample = None
    dp_index = {}                  # FORMAT string -> index of DP in it (-1 when it has none)

    def arrays(contig, need):
        if contig not in depth:
            n = max(lengths.get(contig, 0), need)
            depth[contig] = np.zeros(n, dtype=np.int32)
            known[contig] = np.zeros(n, dtype=bool)
            order.append(contig)
        elif need > len(depth[contig]):
            grow = max(need, 2 * len(depth[contig]))
            depth[contig] = np.concatenate([depth[contig], np.zeros(grow - len(depth[contig]), np.int32)])
            known[contig] = np.concatenate([known[contig], np.zeros(grow - len(known[contig]), bool)])
        return depth[contig], known[contig]

    last = (None, None, None)      # (contig, depth array, mask): most lines stay on one contig
    with _open(path) as fh:
        for line in fh:
            if line[0] == '#':
                if line.startswith('##contig='):
                    body = line.strip()[len('##contig=<'):-1]
                    fields = dict(f.split('=', 1) for f in body.split(',') if '=' in f)
                    try:
                        lengths[fields['ID']] = int(fields.get('length', 0))
                    except (KeyError, ValueError):
                        pass
                elif line.startswith('#CHROM'):
                    cols = line.rstrip('\n').split('\t')
                    sample = cols[9] if len(cols) > 9 else None
                continue
            c = line.rstrip('\n').split('\t')
            if len(c) < 10:
                continue
            k = dp_index.get(c[8])
            if k is None:
                keys = c[8].split(':')
                k = dp_index[c[8]] = keys.index('DP') if 'DP' in keys else -1
            try:
                pos = int(c[1])
            except ValueError:
                continue
            vals = c[9].split(':')
            raw = vals[k] if 0 <= k < len(vals) else '.'
            try:
                dp = int(raw)
            except ValueError:
                continue                   # a depth the record does not state is not a depth of 0
            contig = c[0]
            if last[0] != contig or pos > len(last[1]):
                d, m = arrays(contig, pos)
                last = (contig, d, m)
            last[1][pos - 1] = dp
            last[2][pos - 1] = True
    # Trim what growing over-allocated past the last record of a contig with no stated length.
    for contig in order:
        if not lengths.get(contig):
            idx = np.flatnonzero(known[contig])
            n = int(idx[-1]) + 1 if idx.size else 0
            depth[contig], known[contig] = depth[contig][:n], known[contig][:n]
    return sample, {c: depth[c] for c in order}, {c: known[c] for c in order}


def windows(depth, known, size, callable_dp):
    """[(contig, start, end, mean depth, zero fraction, callable fraction)], 1-based inclusive.

    Fractions are of the positions with a record; a window with none has no values.
    """
    out = []
    for contig, d in depth.items():
        m = known[contig]
        for s in range(0, len(d), size):
            e = min(len(d), s + size)
            dm, km = d[s:e], m[s:e]
            n = int(km.sum())
            if n == 0:
                out.append((contig, s + 1, e, None, None, None))
                continue
            dk = dm[km]
            out.append((contig, s + 1, e, float(dk.mean()), float((dk == 0).mean()),
                        float((dk >= callable_dp).mean())))
    return out


def zero_runs(depth, known, min_run):
    """[(contig, start, end)] of the stretches of at least `min_run` read positions at depth 0."""
    out = []
    for contig, d in depth.items():
        z = ((d == 0) & known[contig]).astype(np.int8)
        edges = np.diff(np.concatenate([[0], z, [0]]))
        starts, ends = np.flatnonzero(edges == 1), np.flatnonzero(edges == -1)
        for s, e in zip(starts, ends):
            if e - s >= min_run:
                out.append((contig, int(s) + 1, int(e)))
    return out


def read_genes(path):
    """[(contig, start, end, name, locus_tag)] of the gene features of a GFF3, CDS where there
    are no genes (a GFF from an annotator that writes only CDS). Empty without a readable file."""
    if not path or not os.path.exists(path) or os.path.basename(path).startswith('NO_FILE'):
        return []
    genes, cds = [], []
    try:
        with _open(path) as fh:
            for line in fh:
                if line.startswith('#') or '\t' not in line:
                    continue
                c = line.rstrip('\n').split('\t')
                if len(c) < 9 or c[2] not in ('gene', 'CDS'):
                    continue
                try:
                    start, end = int(c[3]), int(c[4])
                except ValueError:
                    continue
                name = None
                for key in ('Name', 'gene', 'locus_tag', 'ID'):
                    name = _gff_attr(c[8], key)
                    if name:
                        break
                row = (c[0], start, end, name or f'{start}-{end}', _gff_attr(c[8], 'locus_tag') or '')
                (genes if c[2] == 'gene' else cds).append(row)
    except OSError:
        return []
    return sorted(genes or cds, key=lambda g: (g[0], g[1]))


def gene_depths(depth, known, genes, callable_dp, median):
    """[(contig, start, end, name, locus_tag, length, breadth, callable, mean depth, relative)].

    breadth is the share of the gene's read positions at depth 1 or more, callable the share at
    `callable_dp` or more, relative the mean depth against the genome-wide median (None at
    median 0).
    """
    out = []
    for contig, start, end, name, locus in genes:
        d = depth.get(contig)
        if d is None:
            continue
        s, e = max(0, start - 1), min(len(d), end)
        if e <= s:
            continue
        km = known[contig][s:e]
        dk = d[s:e][km]
        if dk.size == 0:
            out.append((contig, start, end, name, locus, end - start + 1, None, None, None, None))
            continue
        mean = float(dk.mean())
        out.append((contig, start, end, name, locus, end - start + 1, float((dk > 0).mean()),
                    float((dk >= callable_dp).mean()), mean, (mean / median) if median else None))
    return out


def _fmt(v, nd=4):
    return '' if v is None else (f'{v:.{nd}f}' if isinstance(v, float) else str(v))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--vcf', required=True, help="the sample's all-positions VCF")
    ap.add_argument('--gff', default=None, help="the reference's GFF3, for the per-gene table")
    ap.add_argument('--reference', default='', help="the reference id, written into each table")
    ap.add_argument('--window', type=int, default=1000, help="window size in bp")
    # The consensus writes a reference base only where the backbone calls a confident reference,
    # which it does from --allpos_min_cov reads (30); between --consensus_min_dp and that it
    # writes N. A position read at 20x is not callable, and counting it so would make every
    # SNP-per-callable-kb density read low.
    ap.add_argument('--callable-dp', type=int, default=30,
                    help="depth from which the consensus can call a base (--allpos_min_cov)")
    ap.add_argument('--min-run', type=int, default=50,
                    help="shortest stretch without reads that is written out")
    ap.add_argument('--out-prefix', required=True)
    a = ap.parse_args(argv)

    sample, depth, known = read_depths(a.vcf)
    if sample is None:
        sys.stderr.write(f"[depth_profile] {a.vcf} names no sample\n")
        return 1
    read = np.concatenate([depth[c][known[c]] for c in depth]) if depth else np.zeros(0, np.int32)
    median = float(np.median(read)) if read.size else 0.0
    mean = float(read.mean()) if read.size else 0.0
    head = [f"# sample={sample}", f"# reference={a.reference}", f"# genome_len={int(read.size)}",
            f"# median_dp={median:g}", f"# mean_dp={mean:.2f}",
            f"# zero_frac={float((read == 0).mean()) if read.size else 0:.4f}",
            f"# callable_frac={float((read >= a.callable_dp).mean()) if read.size else 0:.4f}",
            f"# callable_dp={a.callable_dp}", f"# window={a.window}", f"# min_run={a.min_run}",
            "# contigs=" + ",".join(f"{c}:{len(depth[c])}" for c in depth)]

    with open(f"{a.out_prefix}.depth_windows.tsv", 'w') as fh:
        fh.write("\n".join(head) + "\ncontig\tstart\tend\tmean_dp\tzero_frac\tcallable_frac\n")
        for row in windows(depth, known, a.window, a.callable_dp):
            fh.write("\t".join(_fmt(v, 2 if i == 3 else 4) for i, v in enumerate(row)) + "\n")
    runs = zero_runs(depth, known, a.min_run)
    with open(f"{a.out_prefix}.zero_depth.tsv", 'w') as fh:
        fh.write("\n".join(head) + "\ncontig\tstart\tend\tlength\n")
        for contig, s, e in runs:
            fh.write(f"{contig}\t{s}\t{e}\t{e - s + 1}\n")
    genes = read_genes(a.gff)
    if genes:
        with open(f"{a.out_prefix}.gene_depth.tsv", 'w') as fh:
            fh.write("\n".join(head) + "\ncontig\tstart\tend\tgene\tlocus_tag\tlength\tbreadth\t"
                                       "callable\tmean_dp\trel_depth\n")
            for row in gene_depths(depth, known, genes, a.callable_dp, median):
                fh.write("\t".join(_fmt(v, 2 if i == 8 else 4) for i, v in enumerate(row)) + "\n")
    sys.stderr.write(f"[depth_profile] {sample}: {int(read.size)} positions, median depth {median:g}, "
                     f"{len(runs)} stretch(es) of {a.min_run}+ bp without reads"
                     + (f", {len(genes)} gene(s)" if genes else "") + "\n")
    return 0


if __name__ == '__main__':
    sys.exit(main())
