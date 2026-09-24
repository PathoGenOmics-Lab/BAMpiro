#!/usr/bin/env python3
"""The codon-level changes of a run: one row per sample and codon where reads carry an MNV.

get_MNV groups a sample's SNVs that fall in one codon and reads the combined codon off the reads
that span it (the BAM). Two SNVs in one codon change one amino acid, and annotating them one base
at a time names two changes that are neither: CAC>GAG is His>Glu, where each base alone says
His>Asp and His>Gln, and GAG>TTG is Glu>Leu where the first base alone says a stop.

Its TSV writes the combined amino acid for every codon that holds more than one SNV, also when no
read carries them together: two SNVs on different molecules of a mixed population (type SNP/MNV,
0 MNV reads). A codon change is only real when molecules carry it, so a row is kept here when at
least --min-reads reads carry every change of the codon at once.

    collect_mnv.py --manifest manifest.tsv --min-reads 2 -o <samplesheet>_mnv.tsv

The manifest has one line per sample: sample, reference, get_MNV TSV (the file names stay free).
Columns: sample, reference, contig, gene, positions, ref_bases, alt_bases, ref_codon, mnv_codon,
aa_change (the codon read whole), snp_aa_changes (each SNV alone: what a per-base annotation says),
change_type, consequence_shift (MNV-masked: the SNVs alone name another consequence, such as a stop;
MNV-gained; Concordant), mnv_reads, total_reads, mnv_frequency, phasing_support (the share of the
reads spanning the codon with its leading change that carry the whole MNV).
"""
from __future__ import annotations

import argparse
import csv
import os
import sys

COLUMNS = ["sample", "reference", "contig", "gene", "positions", "ref_bases", "alt_bases", "ref_codon",
           "mnv_codon", "aa_change", "snp_aa_changes", "change_type", "consequence_shift", "mnv_reads",
           "total_reads", "mnv_frequency", "phasing_support"]


def _int(v):
    try:
        return int(float(v))
    except (TypeError, ValueError):
        return None


def _clean(v):
    v = (v or "").strip()
    return "" if v in ("-", ".") else v


def _list(v):
    """'6139, 6141' -> '6139,6141'."""
    return ",".join(x.strip() for x in _clean(v).split(",") if x.strip())


def read_mnv(path, sample, reference, min_reads):
    """The rows of one get_MNV TSV whose MNV at least min_reads reads carry."""
    out = []
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if "MNV" not in (r.get("Variant Type") or ""):
                continue
            reads = _int(r.get("MNV Reads"))
            if reads is None or reads < min_reads:
                continue
            out.append({
                "sample": sample, "reference": reference, "contig": r.get("Chromosome", ""),
                "gene": _clean(r.get("Gene")), "positions": _list(r.get("Positions")),
                "ref_bases": _list(r.get("Reference Bases")), "alt_bases": _list(r.get("Base Changes")),
                "ref_codon": _clean(r.get("Reference Codon")), "mnv_codon": _clean(r.get("MNV Codon")),
                "aa_change": _clean(r.get("AA Changes")), "snp_aa_changes": _clean(r.get("SNP AA Changes")),
                "change_type": _clean(r.get("Change Type")),
                "consequence_shift": _clean(r.get("MNV Consequence Shift")),
                "mnv_reads": str(reads), "total_reads": _clean(r.get("Total Reads")),
                "mnv_frequency": _clean(r.get("MNV Frequencies")),
                "phasing_support": _clean(r.get("MNV Phasing Support")),
            })
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", required=True, help="TSV: sample, reference, get_MNV TSV (one line per sample)")
    ap.add_argument("--min-reads", type=int, default=2, help="reads that must carry every change of the codon")
    ap.add_argument("-o", "--output", required=True)
    args = ap.parse_args(argv)

    rows, n_files, n_skipped = [], 0, 0
    with open(args.manifest, encoding="utf-8") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3 or not parts[2]:
                continue
            sample, reference, path = parts[0], parts[1], parts[2]
            if not os.path.exists(path) or os.path.getsize(path) == 0:
                n_skipped += 1
                continue
            n_files += 1
            rows += read_mnv(path, sample, reference, args.min_reads)
    rows.sort(key=lambda r: (r["contig"], _int(r["positions"].split(",")[0]) or 0, r["sample"]))
    with open(args.output, "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    masked = sum(1 for r in rows if r["consequence_shift"] == "MNV-masked")
    sys.stderr.write(f"[mnv] {len(rows)} phased codon change(s) with >= {args.min_reads} reads in {n_files} "
                     f"sample(s)" + (f", {masked} of them masked by a per-base annotation" if masked else "")
                     + (f"; {n_skipped} empty or missing file(s)" if n_skipped else "") + f" -> {args.output}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
