#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Consensus FASTA from BAMpiro all.pos VCF (multiallelic-aware) with "original-like" X rules.
NOW INCLUDES: Special CSV Masking support.

- Keep multiallelic support:
  * 3 distinct A/C/G/T bases at a site -> IUPAC (B/D/H/V)
  * 4 distinct bases -> N
- X (mask_char) priority logic:
  1) Position within --exclude intervals.
  2) Position listed in --mask-sites.
  3) Position marked as 'repetitive' (1) or 'blindspot' (1) in --special-mask-csv.
  4) VCF record with any non-PASS FILTER value (e.g., str10, baq_dropout, caller soft filters).
- '-' (nocall_char) ONLY if DP <= --min-dp (default 0).
"""

from __future__ import annotations

import argparse
import gzip
import sys
import csv
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


def get_ro_ao(fmt: str, sample: str) -> Tuple[Optional[int], Optional[int]]:
    """Reference- and alternate-supporting read counts (freebayes RO / AO, or backbone AD). None if unavailable."""
    fs = parse_format_sample(fmt, sample)

    def _int(key: str) -> Optional[int]:
        v = fs.get(key)
        if v and v not in {".", ""}:
            try:
                return int(float(v.split(",")[0]))
            except ValueError:
                return None
        return None

    ro = _int("RO")
    ao = _int("AO")
    if ro is None or ao is None:                       # fall back to AD = ref,alt
        ad = fs.get("AD")
        if ad and ad not in {".", ""}:
            parts = ad.split(",")
            try:
                ro = int(parts[0])
                ao = sum(int(x) for x in parts[1:]) if len(parts) > 1 else 0
            except ValueError:
                pass
    return ro, ao


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


class CsvMasker:
    """
    Parses a CSV with columns: rv_position,pal,homopolymer,GCrich,repetitive,blindspot
    Masks the position if 'repetitive' or 'blindspot' == 1.
    """
    def __init__(self, csv_path: Optional[str]):
        self.points: Set[int] = set()
        if csv_path:
            self._load(csv_path)

    def _load(self, path: str) -> None:
        try:
            with open(path, "r", encoding="utf-8") as fh:
                reader = csv.DictReader(fh)
                # Check for required 'rv_position' header
                if not reader.fieldnames or 'rv_position' not in reader.fieldnames:
                     sys.stderr.write(f"[WARN] CSV missing 'rv_position' header: {path}\n")
                     return
                
                for row in reader:
                    try:
                        pos = int(row['rv_position'])
                        # Parse flags (handle strings '1' or integers)
                        rep = int(row.get('repetitive', 0))
                        blind = int(row.get('blindspot', 0))
                        
                        if rep == 1 or blind == 1:
                            self.points.add(pos)
                            
                    except (ValueError, KeyError):
                        continue
        except FileNotFoundError:
            sys.stderr.write(f"[WARN] Special Mask CSV file not found: {path}\n")

    def masked(self, pos: int) -> bool:
        return pos in self.points


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
        self.counts: Dict[str, int] = {}   # per-contig bases written (truncation guard)
        self._cur: Optional[str] = None

    def start(self, header: str) -> None:
        if self._started and self._line_len != 0:
            self._fh.write("\n")
        self._fh.write(f">{header}\n")
        self._started = True
        self._line_len = 0
        self._cur = header
        self.counts[header] = 0

    def write_base(self, base: str) -> None:
        b = (base or "N")[0].upper()
        self._fh.write(b)
        self._line_len += 1
        if self._cur is not None:
            self.counts[self._cur] += 1
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
            if a == "*":
                # Spanning-deletion allele: skip it and use the remaining allele(s)
                continue
            if len(a) != 1 or a not in DNA_BASES:
                # Non-SNP allele (MNP / indel / complex): signal for masking
                return None
            out.add(a)
        else:
            # Allele index out of range: signal for masking
            return None
    return out if out else None


def consensus_for_position(
    recs: List[VcfRecord],
    ref_seqs: Dict[str, str],
    interval_mask: IntervalMasker,
    point_mask: PointMasker,
    csv_mask: CsvMasker,
    mask_char: str,
    nocall_char: str,
    min_dp: int,
    ref_min_dp: int,
    max_ref_altfrac: float,
    max_ref_min_alt: int,
) -> str:
    chrom = recs[0].chrom
    pos = recs[0].pos

    # CHECK MASKING PRIORITY:
    # 1. Interval Mask (Exclusions)
    # 2. Point Mask (Mask Sites)
    # 3. CSV Mask (Repetitive/Blindspot)
    if interval_mask.masked(chrom, pos) or point_mask.masked(chrom, pos) or csv_mask.masked(pos):
        return mask_char

    for r in recs:
        if r.flt and r.flt not in {".", "PASS", ""}:
            return mask_char

    dp_pos = max(get_dp(r.fmt, r.sample, r.info) for r in recs) if recs else 0
    if dp_pos <= min_dp:
        return nocall_char

    ref_seq = ref_seqs.get(chrom)
    ref_base = "N"
    if ref_seq and 1 <= pos <= len(ref_seq):
        ref_base = ref_seq[pos - 1].upper()
    elif recs[0].ref:
        ref_base = recs[0].ref[0].upper()

    if len(recs) == 1 and (recs[0].alt == "." or recs[0].alt == ""):
        # Monomorphic-reference site: emit the reference base ONLY if the reference call is
        # CONFIDENT, else N. Without this gate, uncertain / low-coverage / deletion-spanning
        # positions default to the reference base (false-ancestral) and bias a phylogeny.
        if ref_base not in DNA_BASES:
            return "N"
        r0 = recs[0]
        gt_alleles = parse_gt_alleles(get_gt(r0.fmt, r0.sample) or "")
        if gt_alleles is None or any(a != 0 for a in gt_alleles):   # GT missing or not homozygous-ref
            return "N"
        ro, ao = get_ro_ao(r0.fmt, r0.sample)
        # Gate on RO+AO (real base-calling depth), NOT raw DP: samtools counts deletion-spanning
        # '*' and ref-skip reads into DP but they carry no base evidence, so a deletion site would
        # otherwise satisfy ref_min_dp on DP alone and emit a false-ancestral reference call.
        base_depth = (ro + ao) if (ro is not None and ao is not None) else dp_pos
        if base_depth < ref_min_dp:
            return "N"
        # Enough alt reads to doubt the ref call: require BOTH a high alt fraction AND >= max_ref_min_alt alt reads.
        if (ro is not None and ao is not None and base_depth > 0
                and ao >= max_ref_min_alt and ao / base_depth >= max_ref_altfrac):
            return "N"
        return ref_base

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
    special_csv_path: Optional[str],
    out_path: str,
    wrap: int,
    mask_char: str,
    nocall_char: str,
    min_dp: int,
    ref_min_dp: int,
    max_ref_altfrac: float,
    max_ref_min_alt: int,
    strict: bool,
) -> None:
    ref_reader = FastaReader(reference_path)
    ref_seqs = ref_reader.seqs

    interval_mask = IntervalMasker(exclude_path)
    point_mask = PointMasker(mask_sites_path)
    csv_mask = CsvMasker(special_csv_path)

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
            # Also check CSV mask for gaps
            if interval_mask.masked(chrom, p) or point_mask.masked(chrom, p) or csv_mask.masked(p):
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

        base = consensus_for_position(pending, ref_seqs, interval_mask, point_mask, csv_mask, mask_char, nocall_char, min_dp, ref_min_dp, max_ref_altfrac, max_ref_min_alt)
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

    # Fail loud on a truncated / interrupted stream: every reference contig must be emitted at
    # full length. A mid-stream failure can still leave a valid-gzip but SHORT FASTA that is
    # plausible-but-wrong for a phylogeny; assert here rather than let it flow downstream.
    missing = [c for c in ref_seqs if c not in writer.counts]
    if missing:
        raise RuntimeError(f"consensus is missing contig(s) entirely: {missing} (empty / truncated VCF?)")
    for chrom, seq in ref_seqs.items():
        got = writer.counts.get(chrom, 0)
        if got != len(seq):
            raise RuntimeError(
                f"consensus length mismatch for {chrom}: wrote {got} != reference {len(seq)} "
                f"(truncated VCF / interrupted stream?)")


def main() -> None:
    ap = argparse.ArgumentParser(description="Consensus FASTA from BAMpiro all.pos.")
    ap.add_argument("--vcf", required=True, help="Input all.pos.vcf(.gz)")
    ap.add_argument("--reference", required=True, help="Reference FASTA (.fa or .fa.gz)")
    ap.add_argument("--exclude", default=None, help="Exclude intervals file")
    ap.add_argument("--mask-sites", default=None, help="Mask sites file")
    
    # New argument for the special CSV mask
    ap.add_argument("--special-mask-csv", default=None, help="CSV with rv_position,repetitive,blindspot to mask")
    
    ap.add_argument("--output", required=True, help="Output consensus FASTA")
    ap.add_argument("--wrap", type=int, default=80, help="FASTA wrap length")
    ap.add_argument("--mask-char", default="X", help="Character for masked sites")
    ap.add_argument("--nocall-char", default="-", help="Character for no-call sites")
    ap.add_argument("--min-dp", type=int, default=0, help="Write nocall_char when DP <= min-dp")
    ap.add_argument("--ref-min-dp", type=int, default=10,
                    help="Min base-supporting depth (RO+AO) to emit the REFERENCE base at a monomorphic site, else N. Prevents reference bias. 0 = disable.")
    ap.add_argument("--max-ref-altfrac", type=float, default=0.10,
                    help="If a monomorphic-reference call has alt-read fraction >= this (and >= --max-ref-min-alt alt reads), write N. 1.0 = disable.")
    ap.add_argument("--max-ref-min-alt", type=int, default=2,
                    help="Only apply --max-ref-altfrac if there are >= this many alt reads (stops a single low-cov read N-ing a ref call).")
    ap.add_argument("--strict", action="store_true", help="Fail on unsorted VCF")
    args = ap.parse_args()

    build_consensus(
        vcf_path=args.vcf,
        reference_path=args.reference,
        exclude_path=args.exclude,
        mask_sites_path=args.mask_sites,
        special_csv_path=args.special_mask_csv,
        out_path=args.output,
        wrap=args.wrap,
        mask_char=(args.mask_char or "X")[0],
        nocall_char=(args.nocall_char or "-")[0],
        min_dp=int(args.min_dp),
        ref_min_dp=int(args.ref_min_dp),
        max_ref_altfrac=float(args.max_ref_altfrac),
        max_ref_min_alt=int(args.max_ref_min_alt),
        strict=args.strict,
    )

if __name__ == "__main__":
    main()
