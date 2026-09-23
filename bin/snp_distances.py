#!/usr/bin/env python3
"""Pairwise SNP distances between consensus sequences, as a phylogeny built on them would see.

Two samples differ at a position when both called a base there (A, C, G or T) and the bases
differ. A position either one did not call, masked, or wrote as an ambiguity code (a mixed site)
says nothing about the pair and is left out of that comparison, so the distance never counts a
gap as a difference. The consensus already carries every masking decision of the pipeline, which
is why the distances are taken from it rather than from the VCFs.

Only samples mapped against the same reference are compared: their consensus sequences share
the reference's coordinates, and sequences of different references do not.

Two outputs:
  -o       square matrix, samples x samples, NA between samples of different references
  --pairs  one row per pair on a reference: the SNPs between them and how many of that reference's
           variable positions both called, which says how much the distance rests on

  snp_distances.py --manifest consensus.tsv --dir cons -o snp_distances.tsv --pairs pairs.tsv

The manifest is 'sample<TAB>reference<TAB>fasta' (file name, looked up under --dir). Without one,
--consensus takes the files directly, names each sample after its first record (or the file name
up to the first dot) and groups them by sequence length.
"""
from __future__ import annotations

import argparse
import gzip
import os
import sys

import numpy as np

ACGT = np.frombuffer(b"ACGT", dtype=np.uint8)


def read_fasta(path):
    """(record names, the records' sequences concatenated as upper-case uint8)."""
    op = gzip.open if str(path).endswith(".gz") else open
    names, chunks = [], []
    with op(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                names.append(line[1:].split()[0] if line[1:].split() else "")
                continue
            chunks.append(line.strip())
    seq = "".join(chunks).upper().encode("ascii", "replace")
    return names, np.frombuffer(seq, dtype=np.uint8)


def entries(args):
    """[(sample, reference, path)] from --manifest or --consensus."""
    out = []
    if args.manifest:
        with open(args.manifest, encoding="utf-8") as fh:
            for line in fh:
                c = line.rstrip("\n").split("\t")
                if len(c) >= 3 and c[0] and not line.startswith("#"):
                    out.append((c[0], c[1], os.path.join(args.dir, c[2]) if args.dir else c[2]))
        return out
    for path in args.consensus or []:
        names, seq = read_fasta(path)
        # The pipeline names a one-record consensus after the sample and prefixes each record of
        # a multi-contig one with it, so only a single record's name is the sample's.
        sample = names[0] if len(names) == 1 and names[0] else os.path.basename(path).split(".")[0]
        out.append((sample, f"L{len(seq)}", path))
    return out


def variable_positions(paths):
    """(indices of the positions where two samples called different bases, length) over `paths`.

    A running 'first base called here' is kept instead of the reference, which this script is
    never given: a later sample calling another base at a position marks it variable.
    """
    first, variable, length = None, None, None
    for path in paths:
        _, seq = read_fasta(path)
        if first is None:
            length = len(seq)
            first = np.zeros(length, dtype=np.uint8)
            variable = np.zeros(length, dtype=bool)
        if len(seq) != length:
            raise ValueError(f"{path}: {len(seq)} bp where the other consensus sequences of its "
                             f"reference have {length}")
        called = np.isin(seq, ACGT)
        variable |= called & (first != 0) & (seq != first)
        new = called & (first == 0)
        first[new] = seq[new]
    return (np.flatnonzero(variable) if variable is not None else np.zeros(0, dtype=np.int64)), length


def codes_at(paths, positions):
    """samples x positions uint8: the base called (A, C, G or T as ASCII) or 0 where none was."""
    out = np.zeros((len(paths), len(positions)), dtype=np.uint8)
    for i, path in enumerate(paths):
        _, seq = read_fasta(path)
        c = seq[positions]
        out[i] = np.where(np.isin(c, ACGT), c, 0)
    return out


def distances(codes):
    """(snps, compared): samples x samples matrices. `compared` counts the positions both called."""
    n = codes.shape[0]
    snps = np.zeros((n, n), dtype=np.int64)
    compared = np.zeros((n, n), dtype=np.int64)
    called = codes != 0
    for i in range(n):
        both = called[i] & called[i:]
        diff = both & (codes[i:] != codes[i])
        snps[i, i:] = snps[i:, i] = diff.sum(axis=1)
        compared[i, i:] = compared[i:, i] = both.sum(axis=1)
    return snps, compared


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", help="sample<TAB>reference<TAB>fasta")
    ap.add_argument("--dir", default="", help="folder the manifest's fasta names are under")
    ap.add_argument("--consensus", nargs="*", default=[], help="consensus FASTAs, without a manifest")
    ap.add_argument("-o", "--output", required=True, help="square matrix TSV")
    ap.add_argument("--pairs", help="one row per pair TSV")
    a = ap.parse_args(argv)

    items = entries(a)
    groups = {}
    for sample, ref, path in items:
        if not os.path.exists(path):
            sys.stderr.write(f"[snp_distances] WARN {sample}: {path} not found\n")
            continue
        groups.setdefault(ref, []).append((sample, path))
    samples = [s for ref in sorted(groups) for s, _ in groups[ref]]
    index = {s: i for i, s in enumerate(samples)}
    full = [["NA"] * len(samples) for _ in samples]
    rows = []
    for ref in sorted(groups):
        members = groups[ref]
        paths = [p for _, p in members]
        pos, length = variable_positions(paths)
        codes = codes_at(paths, pos)
        snps, compared = distances(codes)
        for i, (si, _) in enumerate(members):
            for j, (sj, _) in enumerate(members):
                full[index[si]][index[sj]] = str(int(snps[i, j]))
                if j > i:
                    rows.append((si, sj, ref, int(snps[i, j]), int(compared[i, j]), len(pos)))
        sys.stderr.write(f"[snp_distances] {ref}: {len(members)} sample(s), {length} bp, "
                         f"{len(pos)} variable position(s)\n")

    with open(a.output, "w", encoding="utf-8") as fh:
        fh.write("sample\t" + "\t".join(samples) + "\n")
        for s in samples:
            fh.write(s + "\t" + "\t".join(full[index[s]]) + "\n")
    if a.pairs:
        with open(a.pairs, "w", encoding="utf-8") as fh:
            fh.write("sample_a\tsample_b\treference\tsnps\tcompared\tvariable_sites\n")
            for r in rows:
                fh.write("\t".join(map(str, r)) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
