"""Unit tests for bin/lift_vcf.py: a sample's SNPs moved to the canonical coordinates before annotation."""

from __future__ import annotations

import argparse
import gzip

from conftest import load_script

lv = load_script("lift_vcf")
qc = load_script("qc_report")

HEAD = ("##fileformat=VCFv4.2\n##contig=<ID=E1,length=1000>\n"
        "##INFO=<ID=ANN,Number=.,Type=String,Description=\"old\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n")


def _rec(chrom, pos, ref, alt, info="DP=40;ANN=T|missense_variant|MODERATE|E1_00001|x"):
    return f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t50\tPASS\t{info}\tGT:AD:DP\t1:2,38:40\n"


CANON = [("MTB_anc", "ACGTACGTAC" * 10)]          # canonical base at p is "ACGTACGTAC"[(p - 1) % 10]


def _lift(records, lmap, contigs=CANON, chrom="Chromosome"):
    return lv.lift((HEAD + "".join(records)).splitlines(True), lmap, contigs, chrom)


def test_a_snp_moves_to_its_canonical_position_with_the_canonical_base():
    head, recs, n = _lift([_rec("E1", 500, "C", "T")], {("E1", 500): ("MTB_anc", 42, "+")})

    c = recs[0].split("\t")
    assert c[:5] == ["Chromosome", "42", ".", "C", "T"], "CHROM renamed for the database, REF the canonical base"
    assert c[7] == "OPOS=E1:500;DP=40", "the mapping coordinate kept, the old annotation stripped"
    assert "##contig=<ID=Chromosome,length=100>" in head
    assert not any(h.startswith("##INFO=<ID=ANN,") for h in head)
    assert n["lifted"] == 1


def test_a_reference_lying_reversed_complements_the_allele():
    _, recs, _ = _lift([_rec("E1", 500, "G", "A,C")], {("E1", 500): ("MTB_anc", 42, "-")})

    assert recs[0].split("\t")[3:5] == ["C", "T,G"]


def test_a_sample_carrying_the_canonical_base_has_nothing_to_annotate():
    """A lineage SNP of the mapping reference that the sample reverts: in canonical numbering, no change."""
    _, recs, n = _lift([_rec("E1", 500, "T", "C")], {("E1", 500): ("MTB_anc", 42, "+")})

    assert recs == [] and n["canonical_base"] == 1


def test_records_the_lift_dropped_and_non_snps_are_left_out():
    _, recs, n = _lift([_rec("E1", 7, "C", "T"), _rec("E1", 8, "CA", "C"), _rec("E1", 9, "C", "*")], {})

    assert recs == []
    assert (n["unlifted"], n["not_snp"]) == (1, 2)


def test_records_are_sorted_by_their_canonical_position():
    lmap = {("E1", 100): ("MTB_anc", 90, "+"), ("E1", 200): ("MTB_anc", 10, "+")}
    _, recs, _ = _lift([_rec("E1", 100, "C", "T"), _rec("E1", 200, "C", "A")], lmap)

    assert [r.split("\t")[1] for r in recs] == ["10", "90"]


def test_a_multi_contig_canonical_genome_keeps_its_own_names():
    contigs = [("chr1", "ACGT" * 25), ("chr2", "TTTT" * 25)]
    head, recs, _ = _lift([_rec("E1", 5, "C", "A")], {("E1", 5): ("chr2", 3, "+")}, contigs)

    assert recs[0].split("\t")[:2] == ["chr2", "3"]
    assert "##contig=<ID=chr2,length=100>" in head


def test_the_lifted_record_is_paired_with_its_variant_by_the_report(tmp_path):
    """End to end: lift_vcf writes OPOS, the report pairs by it and takes the canonical coordinate."""
    vcf = tmp_path / "s.vcf.gz"
    with gzip.open(vcf, "wt") as fh:
        fh.write(HEAD + _rec("E1", 500, "C", "T"))
    lmap = tmp_path / "lift.tsv"
    lmap.write_text("src_contig\tsrc_pos\ttgt_contig\ttgt_pos\tstrand\nE1\t500\tMTB_anc\t42\t+\n")
    fa = tmp_path / "canon.fa"
    fa.write_text(">MTB_anc\n" + CANON[0][1] + "\n")
    lifted = tmp_path / "lifted.vcf"

    assert lv.main(["--vcf", str(vcf), "--map", str(lmap), "--fasta", str(fa), "--chrom", "Chromosome",
                    "-o", str(lifted)]) == 0
    # what snpEff adds against the canonical database
    annotated = tmp_path / "canonical.vcf"
    annotated.write_text(lifted.read_text().replace(
        "OPOS=E1:500;DP=40", "OPOS=E1:500;DP=40;ANN=T|missense_variant|MODERATE|rpoB|g|t|t|pc|1/1|c.1A>T|p.Ser450Leu"))
    original = tmp_path / "s.vcf"
    original.write_text(HEAD + _rec("E1", 500, "C", "T"))

    v = qc.load_variants(argparse.Namespace(vcfs=[str(original)], vcfs_h37rv=[str(annotated)],
                                            pos_liftover=None, aa2_label="H37Rv"))["S1"]["E1:500"]

    assert (v["gene_h37rv"], v["aa_h37rv"], v["pos_h37rv"]) == ("rpoB", "p.Ser450Leu", "H37Rv:42")
