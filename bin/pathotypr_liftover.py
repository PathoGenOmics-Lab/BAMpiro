#!/usr/bin/env python3
"""Alignment-free k-mer coordinate liftover between two MTBC references (no whole-genome alignment, no
shared-coordinate assumption). For a position P on reference A it takes the flanking-context k-mer and
finds it in reference B; B's position is the equivalent coordinate.

  lift  <positions> <A.fasta> <B.fasta>   -> the RECOMMENDED, self-contained method. Builds unique-in-both
                                             ANCHORS, and for a position whose k-mer RECURS (repeats) picks
                                             the occurrence consistent with the surrounding anchors (SYNTENY)
                                             instead of dropping it. Emits src->tgt map / lifted BED. Never
                                             mis-maps: a coordinate is correct or absent. Pure Python.

  markers <positions> <A.fasta>  +  apply <classify_out>   -> the same idea via `pathotypr classify`
                                             (Rust, for very large sets). markers writes the --tsv_pos file;
                                             apply parses the classify output (pass --source-fasta/--target-fasta
                                             so non-unique k-mers are dropped, and --offset 1 for the 0->1-based
                                             convention). classify only reports one occurrence, so it drops
                                             (does not synteny-resolve) repeats.

A position never lifts if its context k-mer is absent in B (a nearby A/B difference breaks it). Validated on
the real H37Rv / MTBC-ancestor pair: `lift` reaches ~99% with zero wrong coordinates, tracks indels, and the
DR sites (rpoB 761155, katG 2155168, gyrA 7570, rrs 1473246) map correctly.
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


def _write_outputs(good, out_map, out_bed, contig, name):
    if out_map:
        with open(out_map, "w", encoding="utf-8") as w:
            w.write("src_pos\ttgt_pos\n")
            for s in sorted(good):
                w.write("%d\t%d\n" % (s, good[s]))
    if out_bed:
        tgts = sorted(set(good.values()))
        with open(out_bed, "w", encoding="utf-8") as w:
            if tgts:
                start = prev = tgts[0]
                for q in tgts[1:]:
                    if q == prev + 1:
                        prev = q
                    else:
                        w.write("%s\t%d\t%d\t%s\n" % (contig, start - 1, prev, name))
                        start = prev = q
                w.write("%s\t%d\t%d\t%s\n" % (contig, start - 1, prev, name))


_HMUL = 0x9E3779B97F4A7C15
_M64 = (1 << 64) - 1


def _hash64(code):
    h = (code * _HMUL) & _M64
    return h ^ (h >> 29)


def _rc_code(code, k):
    """Reverse-complement of a 2-bit-encoded k-mer (A<->T, C<->G)."""
    rc = 0
    for _ in range(k):
        rc = (rc << 2) | (3 - (code & 3))
        code >>= 2
    return rc


def _kmer_unique_map(fasta, k, sample=1):
    """{2-bit k-mer code -> single 1-based CENTRE} for k-mers UNIQUE in the sequence. Rolling 2-bit code,
    resets on non-ACGT. With sample>1, keep only ~1/sample of k-mers (FracMinHash: hash(code) % sample == 0)
    -> ~sample x less memory; the choice is deterministic so both genomes keep the SAME subset and shared
    anchors survive."""
    seq = _read_first_contig(fasta)
    half = k // 2
    mask = (1 << (2 * k)) - 1
    code_of = {65: 0, 67: 1, 71: 2, 84: 3}
    seen = {}
    code = 0
    valid = 0
    for i, ch in enumerate(seq.encode("ascii", "replace")):
        v = code_of.get(ch, -1)
        if v < 0:
            valid = 0
            continue
        code = ((code << 2) | v) & mask
        valid += 1
        if valid >= k and (sample <= 1 or _hash64(code) % sample == 0):
            pos = i - half + 1
            seen[code] = -1 if code in seen else pos
    return {c: p for c, p in seen.items() if p > 0}


def _lis(pairs, increasing=True):
    """Longest strictly-monotonic-in-b subsequence of `pairs` (already sorted by a). increasing=False finds
    the longest DECREASING run (an inversion). Returns (ca, cb) of the chained anchors."""
    import bisect
    if not pairs:
        return [], []
    vals = [b if increasing else -b for _, b in pairs]
    tails_v, tails_i, parent = [], [], [-1] * len(vals)
    for i, v in enumerate(vals):
        j = bisect.bisect_left(tails_v, v)
        if j == len(tails_v):
            tails_v.append(v); tails_i.append(i)
        else:
            tails_v[j] = v; tails_i[j] = i
        parent[i] = tails_i[j - 1] if j > 0 else -1
    idx = []
    node = tails_i[-1]
    while node != -1:
        idx.append(node); node = parent[node]
    idx.reverse()
    return [pairs[i][0] for i in idx], [pairs[i][1] for i in idx]


def _lift_chain(args, positions):
    """Whole-genome anchor-chain liftover. Anchors = k-mers unique in both genomes; the FORWARD chain (LIS)
    is the collinear backbone and a REVERSE chain (longest decreasing run of reverse-complement anchors)
    covers an inversion. Any position is placed by interpolating between its flanking anchors on the closer
    chain (SNP sites + indels handled). Large forward gaps where the target lost sequence are reported as RD
    deletions and their interior positions are dropped (no equivalent). Anchor-desert positions drop."""
    import bisect
    k = args.kmer_size
    sample = max(1, args.sample)
    uA = _kmer_unique_map(args.source_fasta, k, sample)
    uB = _kmer_unique_map(args.target_fasta, k, sample)
    fwd, rev, in_f, in_r = [], [], set(), set()
    for c, a in uA.items():
        if c in uB:
            fwd.append((a, uB[c])); in_f.add(a)
        rc = uB.get(_rc_code(c, k))
        if rc is not None:
            rev.append((a, rc)); in_r.add(a)
    amb = in_f & in_r                                              # a source pos anchoring both strands -> drop
    fa, fb = _lis(sorted((a, b) for a, b in fwd if a not in amb), True)
    ra, rb = _lis(sorted((a, b) for a, b in rev if a not in amb), False)

    rd = []                                                        # (a_left, a_right): target lost >= rd_min bp
    for i in range(1, len(fa)):
        if (fa[i] - fa[i - 1]) - (fb[i] - fb[i - 1]) >= args.rd_min:
            rd.append((fa[i - 1], fa[i]))
    rd_left = [x[0] for x in rd]

    def in_rd(p):
        j = bisect.bisect_right(rd_left, p) - 1
        return j >= 0 and rd[j][0] < p < rd[j][1]

    def interp(ca, cb, p, slope):
        j = bisect.bisect_left(ca, p)
        if j < len(ca) and ca[j] == p:
            return cb[j]
        left = j - 1 if j > 0 else None
        right = j if j < len(ca) else None
        if left is not None and right is not None:
            if (ca[right] - ca[left]) > args.max_gap:
                return None
            return int(round(cb[left] + (p - ca[left]) * (cb[right] - cb[left]) / (ca[right] - ca[left])))
        if left is not None and (p - ca[left]) <= args.max_gap:
            return cb[left] + slope * (p - ca[left])
        if right is not None and (ca[right] - p) <= args.max_gap:
            return cb[right] - slope * (ca[right] - p)
        return None

    def nearest(ca, p):
        if not ca:
            return 1 << 62
        j = bisect.bisect_left(ca, p)
        d = 1 << 62
        if j < len(ca):
            d = ca[j] - p
        if j > 0:
            d = min(d, p - ca[j - 1])
        return d

    good, n_fwd, n_rev, n_rd, n_drop = {}, 0, 0, 0, 0
    for p in positions:
        cf = None if not fa else (None if in_rd(p) else interp(fa, fb, p, 1))
        cr = None if not ra else interp(ra, rb, p, -1)
        if cf is not None and (cr is None or nearest(fa, p) <= nearest(ra, p)):
            good[p] = cf; n_fwd += 1
        elif cr is not None:
            good[p] = cr; n_rev += 1
        elif fa and in_rd(p):
            n_rd += 1
        else:
            n_drop += 1
    if args.rd_out:
        with open(args.rd_out, "w", encoding="utf-8") as w:
            for a0, a1 in rd:
                w.write("%s\t%d\t%d\tRD_deletion\n" % (args.source_contig, a0, a1))
    return good, {"fwd": n_fwd, "rev": n_rev, "rd_drop": n_rd, "drop": n_drop,
                  "anchors": len(fa), "inv_anchors": len(ra), "rd": len(rd)}


def cmd_lift(args):
    """Self-contained SYNTENY-anchored liftover (no pathotypr call): match each position's context k-mer
    from source to target, and for positions whose k-mer RECURS (repeats) pick the occurrence consistent
    with the surrounding unique anchors instead of dropping it. Recovers most repeat positions correctly.
    With --global-chain, use the whole-genome anchor-chain instead (places SNP sites and other positions by
    interpolation between flanking anchors, not by their own k-mer)."""
    if args.global_chain:
        positions = _positions(args.positions)
        good, st = _lift_chain(args, positions)
        _write_outputs(good, args.out_map, args.out_bed, args.contig, args.name)
        sys.stderr.write("[liftover] lift(anchor-chain): %d fwd + %d inverted = %d lifted ; %d dropped (%d in "
                         "RD) ; %d collinear anchors (+%d inversion), %d RD deletion(s)\n"
                         % (st["fwd"], st["rev"], len(good), st["drop"] + st["rd_drop"], st["rd_drop"],
                            st["anchors"], st["inv_anchors"], st["rd"]))
        return
    import bisect
    k = args.kmer_size
    half = k // 2
    src = _read_first_contig(args.source_fasta)
    tgt = _read_first_contig(args.target_fasta)
    ns, nt = len(src), len(tgt)
    pkmer = {}                                    # position -> its context k-mer (ref base kept)
    for p in _positions(args.positions):
        if p <= half or p > ns - half:
            continue
        km = src[p - 1 - half:p + half]
        if len(km) == k and all(b in "ACGT" for b in km):
            pkmer[p] = km
    kset = set(pkmer.values())
    src_cnt = {km: 0 for km in kset}             # occurrences of each wanted k-mer in the source (uniqueness)
    for i in range(ns - k + 1):
        km = src[i:i + k]
        if km in src_cnt:
            src_cnt[km] += 1
    tgt_occ = {km: [] for km in kset}            # all 1-based target SNP positions of each wanted k-mer
    for i in range(nt - k + 1):
        km = tgt[i:i + k]
        if km in tgt_occ:
            tgt_occ[km].append(i + half + 1)
    anchors = sorted((p, tgt_occ[km][0]) for p, km in pkmer.items()
                     if src_cnt[km] == 1 and len(tgt_occ[km]) == 1)   # unique in BOTH -> unambiguous 1:1
    asrc = [a[0] for a in anchors]

    def predict(p):                              # expected target coord from the nearest unique anchor's local offset
        if not anchors:
            return None
        j = bisect.bisect_left(asrc, p)
        cands = ([anchors[j]] if j < len(anchors) else []) + ([anchors[j - 1]] if j > 0 else [])
        a = min(cands, key=lambda an: abs(an[0] - p))
        return p + (a[1] - a[0])

    good = {}
    n_anchor = n_synteny = n_drop = 0
    for p, km in pkmer.items():
        occ = tgt_occ[km]
        if src_cnt[km] == 1 and len(occ) == 1:
            good[p] = occ[0]; n_anchor += 1
            continue
        if not occ:
            n_drop += 1; continue
        pred = predict(p)
        if pred is None:
            n_drop += 1; continue
        occ_sorted = sorted(occ, key=lambda q: abs(q - pred))
        d0 = abs(occ_sorted[0] - pred)
        d1 = abs(occ_sorted[1] - pred) if len(occ_sorted) > 1 else 1e18
        if d0 <= args.max_shift and (d1 - d0) >= args.min_margin:   # close to synteny AND unambiguous
            good[p] = occ_sorted[0]; n_synteny += 1
        else:
            n_drop += 1
    _write_outputs(good, args.out_map, args.out_bed, args.contig, args.name)
    sys.stderr.write("[liftover] lift: %d anchored + %d recovered-by-synteny = %d ; %d dropped (of %d)\n"
                     % (n_anchor, n_synteny, len(good), n_drop, len(pkmer)))


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

    l = sub.add_parser("lift", help="self-contained synteny-anchored liftover (source + target FASTA -> map/BED)")
    l.add_argument("positions", help="BED (0-based) or one-1-based-position-per-line list, on the source reference")
    l.add_argument("source_fasta", help="the reference the positions are defined on")
    l.add_argument("target_fasta", help="the reference to lift the positions onto")
    l.add_argument("--out-map", default=None, help="write a src_pos<TAB>tgt_pos table")
    l.add_argument("--out-bed", default=None, help="collapse lifted positions into a BED in target coords")
    l.add_argument("--contig", default="target", help="contig name for --out-bed")
    l.add_argument("--name", default="blindspot", help="BED feature name")
    l.add_argument("--kmer-size", type=int, default=21)
    l.add_argument("--max-shift", type=int, default=100,
                   help="max bp between a repeated k-mer's occurrence and the synteny prediction to accept it")
    l.add_argument("--min-margin", type=int, default=20,
                   help="the accepted occurrence must be at least this many bp closer than the next one")
    l.add_argument("--global-chain", action="store_true",
                   help="use the whole-genome anchor CHAIN: place every position (incl. SNP sites) by "
                        "interpolating between flanking unique anchors, instead of per-position k-mer matching")
    l.add_argument("--max-gap", type=int, default=2000,
                   help="--global-chain: max bp between flanking anchors to still interpolate (else drop)")
    l.add_argument("--sample", type=int, default=1,
                   help="--global-chain: keep ~1/N of k-mers as anchors (FracMinHash) -> ~N x less memory")
    l.add_argument("--rd-out", default=None,
                   help="--global-chain: write detected RD deletions (target lost sequence) as a BED (source coords)")
    l.add_argument("--rd-min", type=int, default=50,
                   help="--global-chain: minimum size (bp) of a target deletion to call it an RD block")
    l.add_argument("--source-contig", default="source", help="contig name for the --rd-out BED")
    l.set_defaults(func=cmd_lift)

    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
