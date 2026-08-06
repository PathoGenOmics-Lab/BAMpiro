#!/usr/bin/env python3
"""Length-aware per-read-PAIR mappability filter (stdlib+numpy; no pysam).

Keep a fragment iff AT LEAST ONE mate's reference span reaches the length at which
its start position becomes genome-unique (ref_span >= min_unique_len[start]).
Mate rescue ONLY for concordant (proper) pairs. Drops-only -> coordinate order
preserved -> caller just re-indexes.

Fixes from adversarial review:
 - FIX1 (soft-clip below k-floor): a mapped read whose ref span < genmap_min_k at a
   locus that IS resolvable (mul != inf) is kept (short span is clipping, not ambiguity).
 - FIX3 (concordance): rescue requires the proper-pair flag 0x2 (else a discordant read
   sitting in the wrong repeat copy gets re-admitted -> the exact false positive).
 - Finding2 (contig mismatch): self on an unmodelled contig fails OPEN (keep) + a loud
   startup check that every BAM @SQ contig is in the track.
Usage: samtools view -h in.bam | filter_..._sam.py --mul MUL.npz --kmin 35 | samtools view -b -o out.bam -
"""
import sys, argparse
import numpy as np

def span(cig):
    if cig == "*":
        return 0
    s = n = 0
    for ch in cig:
        if ch.isdigit():
            n = n * 10 + (ord(ch) - 48)
        else:
            if ch in "MDN=X":
                s += n
            n = 0
    return s

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mul", required=True)
    ap.add_argument("--kmin", type=int, default=35, help="smallest read length probed by genmap")
    ap.add_argument("--sentinel", type=int, default=65535)
    ap.add_argument("--keep-unmapped", action="store_true")
    ap.add_argument("--strict-contigs", action="store_true",
                    help="fail loudly if any BAM @SQ contig is absent from the track (else fail-open with a warning)")
    a = ap.parse_args()
    npz = np.load(a.mul)
    MUL = {k: npz[k] for k in npz.files}
    SENT = a.sentinel
    KMIN = a.kmin
    seen_sq = set()

    def anchored(A, idx, sp):
        if not (0 <= idx < len(A)):
            return False
        mul = int(A[idx])
        if sp >= mul:
            return True
        if sp < KMIN and mul != SENT:     # FIX1: short (clipped) span at a resolvable locus
            return True
        return False

    w = sys.stdout.write
    for line in sys.stdin:
        if not line.strip():                               # a blank line has no fields to index
            continue
        if line[0] == "@":
            if line.startswith("@SQ"):
                for tok in line.split("\t"):
                    if tok.startswith("SN:"):
                        sn = tok[3:].strip()
                        if sn not in MUL:
                            if a.strict_contigs:
                                sys.exit("[filter_reads] ERROR: BAM contig %r absent from the mappability track "
                                         "(wrong reference/track?). Aborting (--strict-contigs)." % sn)
                            sys.stderr.write("[filter_reads] WARNING: contig %r not in track -> fail-open (masking disabled there)\n" % sn)
                        seen_sq.add(sn)
            w(line)
            continue
        f = line.split("\t", 11)
        flag = int(f[1])
        if flag & 0xF04:                                   # unmapped|secondary|qcfail|dup|supplementary
            if a.keep_unmapped and (flag & 0x4) and not (flag & 0x900):
                w(line)
            continue
        A = MUL.get(f[2])
        if A is None:                                      # self on unmodelled contig -> fail-open KEEP
            w(line); continue
        idx = int(f[3]) - 1
        self_anchor = anchored(A, idx, span(f[5]))
        if not (flag & 0x1):                               # SE -> judged alone
            if self_anchor:
                w(line)
            continue
        if self_anchor:                                    # short-circuit before rescue
            w(line); continue
        if not (flag & 0x2) or (flag & 0x8):               # FIX3: rescue only concordant, mate-mapped pairs
            continue
        mrn = f[2] if f[6] == "=" else f[6]
        mA = MUL.get(mrn)
        if mA is None:                                     # mate on unmodelled contig -> fail-open KEEP
            w(line); continue
        midx = int(f[7]) - 1
        if not (0 <= midx < len(mA)):
            continue
        mc = ""
        for tok in (f[11] if len(f) > 11 else "").rstrip("\n").split("\t"):
            if tok.startswith("MC:Z:"):
                mc = tok[5:]; break
        if mc and anchored(mA, midx, span(mc)):
            w(line)

if __name__ == "__main__":
    main()
