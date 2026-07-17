#!/usr/bin/env python3
"""K-mer coordinate liftover between two MTBC references, driven by `pathotypr classify`.

pathotypr is alignment-free: for a marker at position P on a --ref_fasta it builds the k-mer of
flanking context around P and locates that k-mer in a target genome, reporting the target position
(`snp_position`) next to the reference one (`ref_position`). We exploit that as a per-position liftover
that needs NO whole-genome alignment and makes NO global-synteny assumption:

    positions on ref A  --(markers TSV)-->  pathotypr classify --ref_fasta A --fasta_genomes B
                        --(classify out)-->  position map / lifted BED in ref B coordinates

Two subcommands wrap the two halves (pathotypr runs in between, in the container):

  markers  <positions> <ref.fasta>          -> the classify --tsv_pos file (context k-mers use the REF base,
                                               so a position lifts iff its k-mer is conserved+unique in B)
  apply    <classify_out>                    -> src->tgt position map (+ optional lifted BED), with lift rate

Limitations (report them, do not hide them): a position does NOT lift if its flanking k-mer is not unique
in B (repeats -> exactly the blind-spot regions) or is broken by a nearby A/B difference. Calibrate the
coordinate convention once with an A->A (identity) run: snp_position should equal the input position;
pass any constant delta as --offset.
"""
from __future__ import annotations
import argparse
import sys


def _read_first_contig(fasta):
    """Uppercase sequence of the first FASTA record (MTBC is single-contig)."""
    seq = []
    started = False
    with open(fasta, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if started:
                    break
                started = True
                continue
            if started:
                seq.append(line.strip())
    return "".join(seq).upper()


def _positions(path):
    """Yield 1-based positions from either a BED (contig start end ...; 0-based half-open, expanded) or a
    plain one-position-per-line list. Deduplicated, sorted."""
    pos = set()
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#") or line.lower().startswith(("chrom", "track")):
                continue
            c = line.replace(",", "\t").split("\t")
            try:
                if len(c) >= 3 and c[1].lstrip("-").isdigit() and c[2].lstrip("-").isdigit():
                    a, b = int(c[1]), int(c[2])          # BED 0-based half-open -> 1-based positions a+1..b
                    for p in range(a + 1, b + 1):
                        pos.add(p)
                else:
                    pos.add(int(c[0]))                    # plain 1-based position
            except ValueError:
                continue
    return sorted(pos)


def cmd_markers(args):
    seq = _read_first_contig(args.fasta)
    n = len(seq)
    written = skipped = 0
    with open(args.out, "w", encoding="utf-8") as w:
        w.write("pos\talt\tlevel\n")                      # header: pathotypr's csv reader skips row 1
        for p in _positions(args.positions):
            if p < 1 or p > n:
                skipped += 1
                continue
            base = seq[p - 1]
            if base not in "ACGT":                        # N / IUPAC -> no reliable context k-mer
                skipped += 1
                continue
            # alt == ref base -> classify's context k-mer is the pure reference window (position liftover,
            # not variant detection); the position lifts where that exact context recurs uniquely in B.
            w.write("%d\t%s\t%s\n" % (p, base, args.level))
            written += 1
    sys.stderr.write("[liftover] markers: %d written, %d skipped (out-of-range / non-ACGT) -> %s\n"
                     % (written, skipped, args.out))


def cmd_apply(args):
    # A k-mer that recurs in the target genome makes the lift AMBIGUOUS: pathotypr's index keeps the last
    # occurrence, so it would report a WRONG coordinate (seen on real MTBC: a duplicated 51-mer at
    # H37Rv 2,300,000 also sits at 2,306,069). Guard against it: if --target-fasta is given, drop any marker
    # whose k-mer is not unique in the target, so we only ever emit a correct coordinate (or none).
    def _kmer_counts(fasta, k):
        seq = _read_first_contig(fasta)
        d = {}
        for i in range(len(seq) - k + 1):
            km = seq[i:i + k]
            d[km] = d.get(km, 0) + 1
        return d
    # A k-mer must be UNIQUE in both genomes to give a correct coordinate: non-unique in the target means
    # pathotypr's index kept the wrong (last) occurrence; non-unique in the source means the marker itself
    # collapsed onto the wrong source position. Both happen on real MTBC (duplicated segments). Drop either.
    tgt_counts = _kmer_counts(args.target_fasta, args.kmer_size) if args.target_fasta else None
    src_counts = _kmer_counts(args.source_fasta, args.kmer_size) if args.source_fasta else None
    # classify output columns: genome  k-mer  k-merPOS  SNPgenome  SNPreference  lineage
    mapping = {}   # src (SNPreference) -> tgt (SNPgenome + offset)
    nonuniq = 0
    with open(args.classify, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) < 5:
                continue
            try:
                tgt = int(c[3]) + args.offset              # SNPgenome: position in the target genome
                src = int(c[4])                            # SNPreference: the marker's source coordinate
            except ValueError:
                continue
            if (tgt_counts is not None and tgt_counts.get(c[1], 0) != 1) or \
               (src_counts is not None and src_counts.get(c[1], 0) != 1):
                nonuniq += 1                               # k-mer recurs in source or target -> ambiguous, drop
                continue
            # keep the first / drop ambiguous multi-hits (a source that maps to >1 target)
            if src in mapping and mapping[src] != tgt:
                mapping[src] = None
            elif src not in mapping:
                mapping[src] = tgt
    good = {s: t for s, t in mapping.items() if t is not None}
    if args.out_map:
        with open(args.out_map, "w", encoding="utf-8") as w:
            w.write("src_pos\ttgt_pos\n")
            for s in sorted(good):
                w.write("%d\t%d\n" % (s, good[s]))
    if args.out_bed:
        tgts = sorted(set(good.values()))
        with open(args.out_bed, "w", encoding="utf-8") as w:
            if tgts:
                start = prev = tgts[0]
                for q in tgts[1:]:
                    if q == prev + 1:
                        prev = q
                    else:
                        w.write("%s\t%d\t%d\t%s\n" % (args.contig, start - 1, prev, args.name))
                        start = prev = q
                w.write("%s\t%d\t%d\t%s\n" % (args.contig, start - 1, prev, args.name))
    amb = sum(1 for t in mapping.values() if t is None)
    sys.stderr.write("[liftover] apply: %d lifted, %d dropped (k-mer recurs in target), %d ambiguous%s\n"
                     % (len(good), nonuniq, amb, (" -> " + args.out_bed) if args.out_bed else ""))


def main():
    ap = argparse.ArgumentParser(description="k-mer coordinate liftover via pathotypr classify")
    sub = ap.add_subparsers(dest="cmd", required=True)

    m = sub.add_parser("markers", help="positions + ref FASTA -> a pathotypr classify --tsv_pos file")
    m.add_argument("positions", help="BED (0-based) or one-1-based-position-per-line list")
    m.add_argument("fasta", help="the reference the positions are defined on (classify --ref_fasta)")
    m.add_argument("-o", "--out", required=True)
    m.add_argument("--level", default="lift")
    m.set_defaults(func=cmd_markers)

    a = sub.add_parser("apply", help="pathotypr classify output -> src->tgt map (+ lifted BED)")
    a.add_argument("classify", help="the classify per-marker output TSV")
    a.add_argument("--out-map", default=None, help="write a src_pos<TAB>tgt_pos table")
    a.add_argument("--out-bed", default=None, help="collapse lifted positions into a BED in target coords")
    a.add_argument("--contig", default="target", help="contig name for --out-bed")
    a.add_argument("--name", default="blindspot", help="BED feature name")
    a.add_argument("--offset", type=int, default=0, help="constant added to snp_position (calibrate with an identity run)")
    a.add_argument("--target-fasta", default=None,
                   help="the genome that was searched (classify --tsv_genomes); drop positions whose k-mer is "
                        "not unique in it, so a lifted coordinate is never wrong (only correct or absent).")
    a.add_argument("--source-fasta", default=None,
                   help="the reference the markers were built on (classify --ref_fasta); drop positions whose "
                        "k-mer is not unique in it (the marker collapsed onto the wrong source position).")
    a.add_argument("--kmer-size", type=int, default=21, help="k-mer size used in classify (for the uniqueness check)")
    a.set_defaults(func=cmd_apply)

    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
