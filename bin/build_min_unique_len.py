#!/usr/bin/env python3
"""Collapse a genmap (k,E)-mappability k-sweep into per-contig min_unique_len tracks.

min_unique_len[pos] = smallest probed read length K at which the K-mer STARTING at pos
is genome-unique (genmap mappability == 1; both strands; <= E mismatches). A position
never unique across the sweep -> 65535 (infinity / always-repetitive). Contig-end 'tail'
positions genmap could not score (pos > L-K for every K) are distinguished from
always-repetitive and handled by --tail-policy.

Output: <prefix>.min_unique_len.npz -- numpy .npz, one uint16 array per contig, key =
contig name exactly as in the FASTA/BAM header; value 65535 = infinity.
Coordinates: 0-based half-open, identical to genmap bedgraph 'start' and to BAM
reference_start, so the read filter looks up min_unique_len[read.reference_start] directly.
"""
import argparse
import numpy as np


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fai", required=True, help="reference.fa.fai")
    ap.add_argument("--bedgraphs", nargs="+", required=True, help="'K:path.bedgraph' tokens")
    ap.add_argument("--sentinel", type=int, default=65535, help="uint16 value = never unique (infinity)")
    ap.add_argument("--tail-policy", choices=["keep", "drop"], default="keep",
                    help="contig-end positions genmap never scored: keep (min_unique_len=0) or drop (=infinity)")
    ap.add_argument("--unique-eps", type=float, default=0.999999, help="mappability >= this counts as unique (==1)")
    ap.add_argument("--out-prefix", required=True)
    a = ap.parse_args()

    contigs = [(l.split("\t")[0], int(l.split("\t")[1])) for l in open(a.fai) if l.strip()]
    SENT = np.uint16(a.sentinel)
    tracks = {n: np.full(L, SENT, np.uint16) for n, L in contigs}
    covered = {n: np.zeros(L, dtype=bool) for n, L in contigs}

    items = sorted((int(t.split(":", 1)[0]), t.split(":", 1)[1]) for t in a.bedgraphs)
    kmin = items[0][0]
    for K, path in items:                                  # ascending K -> first (smallest) unique wins
        kv = np.uint16(min(K, a.sentinel - 1))             # never store the sentinel itself as a real length
        for line in open(path):
            if not line or line[0] == "#":
                continue
            c, s, e, v = line.rstrip("\n").split("\t")
            arr = tracks.get(c)
            if arr is None:
                continue
            s, e = int(s), int(e)
            if K == kmin:
                covered[c][s:e] = True                     # smallest K = widest genmap coverage
            if float(v) >= a.unique_eps:
                seg = arr[s:e]
                seg[seg == SENT] = kv                       # only fill positions not yet unique at a smaller K
                arr[s:e] = seg

    if a.tail_policy == "keep":                            # tail = never scored (contig end), NOT repetitive
        for n, L in contigs:
            tail = (tracks[n] == SENT) & (~covered[n])
            tracks[n][tail] = 0                            # span >= 0 always -> such reads are kept

    np.savez(a.out_prefix + ".min_unique_len.npz", **tracks)


if __name__ == "__main__":
    main()
