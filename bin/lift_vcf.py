#!/usr/bin/env python3
"""A sample's SNPs moved to the canonical reference's coordinates, to be annotated against its database.

A variant VCF is in the coordinates of the reference the sample was mapped against. Annotated as it stands
against the canonical (H37Rv) snpEff database, each position is read as an H37Rv one: on any other reference
that names another gene and codon, and on every reference it names nothing at all, since the database calls
its chromosome 'Chromosome' and a mapping reference does not. This writes each SNP where the liftover map of
LIFT_VARIANTS puts it:

  - CHROM is the canonical contig, or --chrom when the canonical genome has one contig (the name the
    database gives it), and POS the lifted position;
  - REF is the canonical base there, and ALT the sample's allele, complemented where the mapping
    reference lies reversed against the canonical genome;
  - INFO/OPOS keeps the mapping coordinate, which is how the report pairs the record with its variant;
  - a record whose allele IS the canonical base is left out: in canonical numbering nothing changed.

Records the lift dropped, and anything but a SNP (the report reads SNPs only), are left out too. Every ANN,
which described the mapping reference, is stripped, and the records are sorted by their new position, as
tabix needs.
"""
from __future__ import annotations

import argparse
import gzip
import sys

COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}
_OLD_ANNOTATION = ("ANN=", "LOF=", "NMD=", "EFF=")


def read_map(path):
    """{(contig, position): (canonical contig, position, strand)} from a liftover map
    (src_contig src_pos tgt_contig tgt_pos strand)."""
    out = {}
    try:
        with open(path, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                c = line.rstrip("\n").split("\t")
                if len(c) >= 5 and c[1].isdigit() and c[3].isdigit():
                    out[(c[0], int(c[1]))] = (c[2], int(c[3]), c[4])
    except OSError:
        pass
    return out


def read_fasta(path):
    """[(contig, uppercase sequence)] in file order."""
    out, name, parts = [], None, []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    out.append((name, "".join(parts).upper()))
                head = line[1:].strip()
                name, parts = (head.split()[0] if head else ""), []
            elif name is not None:
                parts.append(line.strip())
    if name is not None:
        out.append((name, "".join(parts).upper()))
    return out


def _open(path):
    return gzip.open(path, "rt", encoding="utf-8", errors="replace") if str(path).endswith(".gz") \
        else open(path, encoding="utf-8", errors="replace")


def lift(lines, lmap, contigs, chrom=None):
    """(header lines, record lines sorted by canonical position, counts) for the lines of a VCF."""
    seq = dict(contigs)
    order = {n: i for i, (n, _) in enumerate(contigs)}
    rename = chrom if (chrom and len(contigs) == 1) else None
    head, recs, used, chrom_line = [], [], set(), None
    n = {"snps": 0, "lifted": 0, "unlifted": 0, "canonical_base": 0, "not_snp": 0}
    for line in lines:
        if line.startswith("##"):
            if not line.startswith(("##contig=", "##SnpEff", "##INFO=<ID=ANN,", "##INFO=<ID=LOF,",
                                    "##INFO=<ID=NMD,", "##INFO=<ID=EFF,", "##INFO=<ID=OPOS,")):
                head.append(line.rstrip("\n"))
            continue
        if line.startswith("#"):
            chrom_line = line.rstrip("\n")
            continue
        c = line.rstrip("\n").split("\t")
        if len(c) < 8 or not c[1].isdigit():
            continue
        ref, alts = c[3].upper(), c[4].upper().split(",")
        if ref not in COMP or any(a not in COMP for a in alts):
            n["not_snp"] += 1
            continue
        n["snps"] += 1
        hit = lmap.get((c[0], int(c[1])))
        if hit is None or hit[0] not in seq or not 1 <= hit[1] <= len(seq[hit[0]]):
            n["unlifted"] += 1
            continue
        tc, tp, strand = hit
        cref = seq[tc][tp - 1]
        calts = alts if strand == "+" else [COMP[a] for a in alts]
        if cref not in COMP or calts[0] == cref:
            n["canonical_base"] += 1
            continue
        info = [f for f in c[7].split(";") if f and f != "." and not f.startswith(_OLD_ANNOTATION)
                and not f.startswith("OPOS=")]
        c[7] = ";".join(["OPOS=%s:%s" % (c[0], c[1])] + info)
        c[0], c[1], c[3], c[4] = rename or tc, str(tp), cref, ",".join(calts)
        used.add(tc)
        recs.append(((order[tc], tp), "\t".join(c)))
        n["lifted"] += 1
    for tc in sorted(used, key=order.get):
        head.append("##contig=<ID=%s,length=%d>" % (rename or tc, len(seq[tc])))
    head.append('##INFO=<ID=OPOS,Number=1,Type=String,Description="The contig:position of the reference the '
                'sample was mapped against, which this record was lifted from">')
    head.append(chrom_line or "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
    return head, [r for _, r in sorted(recs, key=lambda x: x[0])], n


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vcf", required=True, help="the sample's VCF, in mapping-reference coordinates (may be .gz)")
    ap.add_argument("--map", required=True, help="the liftover map (src_contig src_pos tgt_contig tgt_pos strand)")
    ap.add_argument("--fasta", required=True, help="the canonical genome the map lifts onto (for its bases)")
    ap.add_argument("--chrom", default=None,
                    help="the chromosome name of the canonical annotation database, written in place of the "
                         "canonical FASTA's contig name when that genome has one contig")
    ap.add_argument("-o", "--output", required=True)
    a = ap.parse_args(argv)
    with _open(a.vcf) as fh:
        head, recs, n = lift(fh, read_map(a.map), read_fasta(a.fasta), a.chrom)
    with open(a.output, "w", encoding="utf-8") as out:
        out.write("\n".join(head + recs) + "\n")
    sys.stderr.write("[lift_vcf] %s: %d SNP(s), %d lifted, %d where the lift has no coordinate, %d carrying the "
                     "canonical base, %d other record(s) left out\n"
                     % (a.vcf, n["snps"], n["lifted"], n["unlifted"], n["canonical_base"], n["not_snp"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
