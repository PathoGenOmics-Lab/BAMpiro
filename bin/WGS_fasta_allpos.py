#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Consensus FASTA from BAMpiro all.pos VCF (multiallelic-aware) with "original-like" X rules.

- Keep multiallelic support:
  * 3 distinct A/C/G/T bases at a site -> IUPAC (B/D/H/V)
  * 4 distinct bases -> N
- X (mask_char) in as many "original-like" cases as possible:
  1) any position within --exclude intervals
  2) any position listed in --mask-sites (e.g., "str10-equivalent" sites from pipeline)
  3) any VCF record with FILTER containing "str10" (VarScan compatibility)
- '-' (nocall_char) ONLY if DP <= --min-dp (default 0)
  (DP read from FORMAT/DP, else INFO/DP, else INFO/ADP, else 0)

Notes:
- all.pos can contain multiple records at the same CHROM/POS (bcftools norm -m -).
  We group records per position and compute the allele set using GT per record.
- We do NOT treat INFO/NC as gap. (NC in BAMpiro backbone can represent DP<MINCOV;
  gaps only when DP = 0.)
"""

from __future__ import annotations

import argparse
import gzip
import sys
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Set, Tuple

DNA_BASES: Set[str] = {"A", "C", "G", "T"}

_IUPAC: Dict[frozenset, str] = {
    # 2-base
    frozenset({"A", "G"}): "R",
    frozenset({"C", "T"}): "Y",
    frozenset({"G", "C"}): "S",
    frozenset({"A", "T"}): "W",
    frozenset({"G", "T"}): "K",
    frozenset({"A", "C"}): "M",
    # 3-base
    frozenset({"C", "G", "T"}): "B",
    frozenset({"A", "G", "T"}): "D",
    frozenset({"A", "C", "T"}): "H",
    frozenset({"A", "C", "G"}): "V",
    # 4-base
    frozenset({"A", "C", "G", "T"}): "N",
}


def iupac_for_bases(bases: Iterable[str]) -> str:
    norm: Set[str] = set()
    for b in bases:
        if not b:
            continue
        b = b.upper()
        if b in DNA_BASES:
            norm.add(b)
    if not norm:
        return "N"
    if len(norm) == 1:
        return next(iter(norm))
    return _IUPAC.get(frozenset(norm), "N")


def parse_info_field(info_field: str) -> Dict[str, str]:
    info: Dict[str, str] = {}
    if not info_field or info_field == ".":
        return info
    for item in info_field.split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" in item:
            k, v = item.split("=", 1)
            info[k] = v
        else:
            info[item] = "true"
    return info


def parse_format_sample(fmt: str, sample: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    if not fmt or fmt == ".":
        return out
    keys = fmt.split(":")
    vals = sample.split(":") if sample else []
    for i, k in enumerate(keys):
        out[k] = vals[i] if i < len(vals) else "."
    return out


def get_dp(fmt: str, sample: str, info: str) -> int:
    fs = parse_format_sample(fmt, sample)
    dp_s = fs.get("DP")
    if dp_s and dp_s not in {".", ""}:
        try:
            return int(float(dp_s))
        except ValueError:
            pass

    d = parse_info_field(info)
    for key in ("DP", "ADP"):
        v = d.get(key)
        if v and v not in {".", ""}:
            try:
                return int(float(v))
            except ValueError:
                continue

    return 0


def get_gt(fmt: str, sample: str) -> Optional[str]:
    fs = parse_format_sample(fmt, sample)
    gt = fs.get("GT")
    if not gt or gt in {".", "./.", ".|.", "./", ".|"}:
        return None
    return gt


def parse_gt_alleles(gt: str) -> Optional[List[int]]:
    if not gt or gt in {".", "./.", ".|."}:
        return None
    sep = "/" if "/" in gt else ("|" if "|" in gt else None)
    parts = gt.split(sep) if sep else [gt]
    alleles: List[int] = []
    for p in parts:
        p = p.strip()
        if not p or p == ".":
            return None
        try:
            alleles.append(int(p))
        except ValueError:
            return None
    return alleles


@dataclass(frozen=True)
class VcfRecord:
    chrom: str
    pos: int
    ref: str
    alt: str
    qual: str
    flt: str
    info: str
    fmt: str
    sample: str


class IntervalMasker:
    def __init__(self, exclude_file: Optional[str]):
        self.intervals: Dict[str, List[Tuple[int, int]]] = {}
        self._ptr: Dict[str, int] = {}
        if exclude_file:
            self._load(exclude_file)

    def _load(self, path: str) -> None:
        try:
            with open(path, "r", encoding="utf-8") as fh:
                for line in fh:
                    line = line.strip()
                    if not line or line.startswith("#"):
                        continue
                    if line.lower().startswith("chrom"):
                        continue
                    parts = line.split("\t")
                    if len(parts) < 3:
                        continue
                    chrom = parts[0].strip()
                    try:
                        start = int(parts[1]); end = int(parts[2])
                    except ValueError:
                        continue
                    if start > end:
                        start, end = end, start
                    self.intervals.setdefault(chrom, []).append((start, end))
        except FileNotFoundError:
            sys.stderr.write(f"[WARN] Exclude file not found: {path}\n")
            return

        for chrom, ivs in self.intervals.items():
            ivs.sort()
            merged: List[Tuple[int, int]] = []
            for s, e in ivs:
                if not merged or s > merged[-1][1] + 1:
                    merged.append((s, e))
                else:
                    merged[-1] = (merged[-1][0], max(merged[-1][1], e))
            self.intervals[chrom] = merged
            self._ptr[chrom] = 0

    def masked(self, chrom: str, pos: int) -> bool:
        ivs = self.intervals.get(chrom)
        if not ivs:
            return False
        i = self._ptr.get(chrom, 0)
        while i < len(ivs) and ivs[i][1] < pos:
            i += 1
        self._ptr[chrom] = i
        if i >= len(ivs):
            return False
        s, e = ivs[i]
        return s <= pos <= e


class PointMasker:
    def __init__(self, mask_file: Optional[str]):
        self.points: Dict[str, Set[int]] = {}
        if mask_file:
            self._load(mask_file)

    def _load(self, path: str) -> None:
        try:
            with open(path, "r", encoding="utf-8") as fh:
                for line in fh:
                    line = line.strip()
                    if not line or line.startswith("#"):
                        continue
                    parts = line.split("\t")
                    if len(parts) < 2:
                        continue
                    chrom = parts[0].strip()
                    try:
                        pos = int(parts[1])
                    except ValueError:
                        continue
                    self.points.setdefault(chrom, set()).add(pos)
        except FileNotFoundError:
            sys.stderr.write(f"[WARN] Mask-sites file not found: {path}\n")

    def masked(self, chrom: str, pos: int) -> bool:
        s = self.points.get(chrom)
        return (pos in s) if s else False


class FastaReader:
    def __init__(self, fasta_path: str):
        self.seqs: Dict[str, str] = {}
        self._load(fasta_path)

    def _load(self, fasta_path: str) -> None:
        op = gzip.open if fasta_path.endswith(".gz") else open
        name: Optional[str] = None
        parts: List[str] = []
        with op(fasta_path, "rt", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                line = line.rstrip("\n")
                if line.startswith(">"):
                    if name is not None:
                        self.seqs[name] = "".join(parts).upper()
                    name = line[1:].split()[0]
                    parts = []
                else:
                    parts.append(line.strip())
            if name is not None:
                self.seqs[name] = "".join(parts).upper()


class FastaWriter:
    def __init__(self, out_path: str, wrap: int):
        self.wrap = max(0, int(wrap))
        self._fh = open(out_path, "w", encoding="utf-8")
        self._line_len = 0
        self._started = False

    def start(self, header: str) -> None:
        if self._started and self._line_len != 0:
            self._fh.write("\n")
        self._fh.write(f">{header}\n")
        self._started = True
        self._line_len = 0

    def write_base(self, base: str) -> None:
        b = (base or "N")[0].upper()
        self._fh.write(b)
        self._line_len += 1
        if self.wrap and self._line_len >= self.wrap:
            self._fh.write("\n")
            self._line_len = 0

    def close(self) -> None:
        if self._started and self._line_len != 0:
            self._fh.write("\n")
        self._fh.close()


def parse_vcf_record(line: str) -> Optional[VcfRecord]:
    fields = line.rstrip("\n").split("\t")
    if len(fields) < 10:
        return None
    try:
        pos = int(fields[1])
    except ValueError:
        return None
    return VcfRecord(
        chrom=fields[0],
        pos=pos,
        ref=fields[3],
        alt=fields[4],
        qual=fields[5],
        flt=fields[6],
        info=fields[7],
        fmt=fields[8],
        sample=fields[9],
    )


def allele_bases_from_record(rec: VcfRecord, ref_base: str) -> Optional[Set[str]]:
    gt = get_gt(rec.fmt, rec.sample)
    if gt is None:
        return None
    alleles = parse_gt_alleles(gt)
    if alleles is None:
        return None

    # If ALT is '.', treat as reference
    if rec.alt == "." or rec.alt == "":
        out = set()
        for aidx in alleles:
            if aidx == 0 and ref_base in DNA_BASES:
                out.add(ref_base)
        return out if out else None

    alt_alleles = [a.strip().upper() for a in rec.alt.split(",") if a.strip()]
    out: Set[str] = set()
    for aidx in alleles:
        if aidx == 0:
            if ref_base in DNA_BASES:
                out.add(ref_base)
        elif 1 <= aidx <= len(alt_alleles):
            a = alt_alleles[aidx - 1]
            if len(a) != 1 or a not in DNA_BASES:
                return {"N"}
            out.add(a)
        else:
            return {"N"}
    return out


def consensus_for_position(
    recs: List[VcfRecord],
    ref_seqs: Dict[str, str],
    interval_mask: IntervalMasker,
    point_mask: PointMasker,
    mask_char: str,
    nocall_char: str,
    min_dp: int,
) -> str:
    chrom = recs[0].chrom
    pos = recs[0].pos

    # X priority
    if interval_mask.masked(chrom, pos) or point_mask.masked(chrom, pos):
        return mask_char

    for r in recs:
        if r.flt and "str10" in r.flt:
            return mask_char

    # DP threshold for gap: use max DP at site
    dp_pos = max(get_dp(r.fmt, r.sample, r.info) for r in recs) if recs else 0
    if dp_pos <= min_dp:
        return nocall_char

    # Reference base
    ref_seq = ref_seqs.get(chrom)
    ref_base = "N"
    if ref_seq and 1 <= pos <= len(ref_seq):
        ref_base = ref_seq[pos - 1].upper()
    elif recs[0].ref:
        ref_base = recs[0].ref[0].upper()

    # Pure backbone record
    if len(recs) == 1 and (recs[0].alt == "." or recs[0].alt == ""):
        return ref_base if ref_base in DNA_BASES else "N"

    bases_union: Set[str] = set()
    for r in recs:
        bset = allele_bases_from_record(r, ref_base)
        if bset is None:
            return mask_char
        if "N" in bset:
            return "N"
        bases_union |= bset

    return iupac_for_bases(bases_union)


def build_consensus(
    vcf_path: str,
    reference_path: str,
    exclude_path: Optional[str],
    mask_sites_path: Optional[str],
    out_path: str,
    wrap: int,
    mask_char: str,
    nocall_char: str,
    min_dp: int,
    strict: bool,
) -> None:
    ref_reader = FastaReader(reference_path)
    ref_seqs = ref_reader.seqs

    interval_mask = IntervalMasker(exclude_path)
    point_mask = PointMasker(mask_sites_path)

    op = gzip.open if vcf_path.endswith(".gz") else open
    writer = FastaWriter(out_path, wrap=wrap)

    current_chrom: Optional[str] = None
    expected_pos = 1
    pending: List[VcfRecord] = []

    def start_contig(chrom: str) -> None:
        nonlocal current_chrom, expected_pos
        current_chrom = chrom
        expected_pos = 1
        writer.start(chrom)

    def fill_gap(chrom: str, start: int, end_excl: int) -> None:
        nonlocal expected_pos
        for p in range(start, end_excl):
            if interval_mask.masked(chrom, p) or point_mask.masked(chrom, p):
                writer.write_base(mask_char)
            else:
                writer.write_base(nocall_char)
        expected_pos = end_excl

    def flush_pending() -> None:
        nonlocal pending, current_chrom, expected_pos
        if not pending:
            return
        chrom = pending[0].chrom
        pos = pending[0].pos

        if current_chrom != chrom:
            start_contig(chrom)

        if chrom not in ref_seqs:
            if strict:
                raise RuntimeError(f"Contig '{chrom}' not found in reference")
            pending = []
            return

        if pos > expected_pos:
            fill_gap(chrom, expected_pos, pos)
        elif pos < expected_pos:
            if strict:
                raise RuntimeError(f"Unexpected duplicate/out-of-order {chrom}:{pos}")
            pending = []
            return

        base = consensus_for_position(pending, ref_seqs, interval_mask, point_mask, mask_char, nocall_char, min_dp)
        writer.write_base(base)
        expected_pos = pos + 1
        pending = []

    with op(vcf_path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            rec = parse_vcf_record(line)
            if rec is None:
                continue

            if not pending:
                pending = [rec]
                continue

            if rec.chrom == pending[0].chrom and rec.pos == pending[0].pos:
                pending.append(rec)
            else:
                flush_pending()
                pending = [rec]

    flush_pending()
    writer.close()


def main() -> None:
    ap = argparse.ArgumentParser(description="Consensus FASTA from BAMpiro all.pos (multiallelic + original-like X, DP<=min_dp -> '-').")
    ap.add_argument("--vcf", required=True, help="Input all.pos.vcf(.gz)")
    ap.add_argument("--reference", required=True, help="Reference FASTA (.fa or .fa.gz)")
    ap.add_argument("--exclude", default=None, help="Exclude intervals file (Chrom\\tStart\\tEnd...)")
    ap.add_argument("--mask-sites", default=None, help="Mask sites file (CHROM\\tPOS per line) -> X (str10-equivalent)")
    ap.add_argument("--output", required=True, help="Output consensus FASTA")
    ap.add_argument("--wrap", type=int, default=80, help="FASTA wrap length (0 disables wrapping)")
    ap.add_argument("--mask-char", default="X", help="Character for masked sites")
    ap.add_argument("--nocall-char", default="-", help="Character for no-call sites")
    ap.add_argument("--min-dp", type=int, default=0, help="Write nocall_char when DP <= min-dp (default 0)")
    ap.add_argument("--strict", action="store_true", help="Fail on unsorted VCF / unexpected duplicates")
    args = ap.parse_args()

    build_consensus(
        vcf_path=args.vcf,
        reference_path=args.reference,
        exclude_path=args.exclude,
        mask_sites_path=args.mask_sites,
        out_path=args.output,
        wrap=args.wrap,
        mask_char=(args.mask_char or "X")[0],
        nocall_char=(args.nocall_char or "-")[0],
        min_dp=int(args.min_dp),
        strict=args.strict,
    )


if __name__ == "__main__":
    main()

