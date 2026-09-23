"""The canonical (H37Rv) coordinate, amino acid and gene the report attaches to each variant (bin/qc_report.py).

A cohort can be mapped against several references, none of them H37Rv. Each variant's canonical coordinate
comes from the liftover map of ITS reference, looked up by contig and position: the same position on another
contig or reference is another site. A canonical record lends a coordinate only when it was lifted to the
canonical reference, which it says by carrying the mapping coordinate in OPOS.
"""

from __future__ import annotations

import argparse

from conftest import load_script

qc = load_script("qc_report")


def _vcf(path, sample, records):
    """records: (chrom, pos, ref, alt, info)."""
    head = ("##fileformat=VCFv4.2\n"
            f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n")
    body = "".join(f"{c}\t{p}\t.\t{r}\t{a}\t50\tPASS\t{i}\tGT:AD:DP\t1:2,38:40\n" for c, p, r, a, i in records)
    path.write_text(head + body)
    return str(path)


def _ann(gene, aa):
    return f"ANN=T|missense_variant|MODERATE|{gene}|g|transcript|t|protein_coding|1/1|c.1A>T|{aa}"


def _args(vcfs, vcfs_h37rv=(), pos_liftover=None):
    return argparse.Namespace(vcfs=list(vcfs), vcfs_h37rv=list(vcfs_h37rv), pos_liftover=pos_liftover,
                              aa2_label="H37Rv")


def _map(path, rows):
    path.write_text("src_contig\tsrc_pos\ttgt_contig\ttgt_pos\tstrand\n"
                    + "".join(f"{c}\t{p}\tMTB_anc\t{t}\t{s}\n" for c, p, t, s in rows))
    return str(path)


def test_read_lift_map_reads_the_contig_aware_map_and_the_older_one(tmp_path):
    new = _map(tmp_path / "new.tsv", [("E1", 100, 761155, "+")])
    old = tmp_path / "old.tsv"
    old.write_text("src_pos\ttgt_pos\n100\t761155\n")

    assert qc.read_lift_map(new) == {("E1", 100): 761155}
    assert qc.read_lift_map(str(old)) == {(None, 100): 761155}
    assert qc.read_lift_map(str(tmp_path / "NO_FILE_LIFTOVER")) == {}


def test_each_reference_gets_its_own_canonical_coordinate(tmp_path):
    """Two references, the same position: each variant takes the coordinate its own reference lifts to.
    One map for the first reference, applied by position alone, gave both the first one's."""
    a = _vcf(tmp_path / "a.vcf", "A", [("E1", 100, "C", "T", ".")])
    b = _vcf(tmp_path / "b.vcf", "B", [("E2", 100, "C", "T", "."), ("E2", 200, "G", "A", ".")])
    lift = _map(tmp_path / "lift.tsv", [("E1", 100, 761155, "+"), ("E2", 100, 900001, "-")])

    v = qc.load_variants(_args([a, b], pos_liftover=lift))

    assert v["A"]["E1:100"]["pos_h37rv"] == "H37Rv:761155"
    assert v["B"]["E2:100"]["pos_h37rv"] == "H37Rv:900001"
    assert not v["B"]["E2:200"].get("pos_h37rv"), "a position the lift dropped has no canonical coordinate"


def test_a_map_without_contigs_is_applied_only_to_a_single_contig_cohort(tmp_path):
    old = tmp_path / "old.tsv"
    old.write_text("src_pos\ttgt_pos\n100\t761155\n")
    a = _vcf(tmp_path / "a.vcf", "A", [("E1", 100, "C", "T", ".")])
    b = _vcf(tmp_path / "b.vcf", "B", [("E2", 100, "C", "T", ".")])

    one = qc.load_variants(_args([a], pos_liftover=str(old)))
    two = qc.load_variants(_args([a, b], pos_liftover=str(old)))

    assert one["A"]["E1:100"]["pos_h37rv"] == "H37Rv:761155"
    assert not two["A"]["E1:100"].get("pos_h37rv") and not two["B"]["E2:100"].get("pos_h37rv")


def test_a_lifted_canonical_record_lends_its_coordinate_amino_acid_and_gene(tmp_path):
    a = _vcf(tmp_path / "a.vcf", "A", [("E1", 100, "C", "T", _ann("E1_00667", "p.Ser451Leu"))])
    h = _vcf(tmp_path / "h.vcf", "A", [("Chromosome", 761155, "C", "T",
                                        "OPOS=E1:100;" + _ann("rpoB", "p.Ser450Leu"))])

    v = qc.load_variants(_args([a], vcfs_h37rv=[h]))["A"]["E1:100"]

    assert (v["aa_h37rv"], v["gene_h37rv"], v["pos_h37rv"]) == ("p.Ser450Leu", "rpoB", "H37Rv:761155")


def test_a_canonical_record_annotated_in_place_lends_no_coordinate(tmp_path):
    """A VCF annotated against the canonical database without being lifted sits at the mapping coordinate.
    Its position is not a canonical one, so it must not stand in for the lift, as it used to."""
    a = _vcf(tmp_path / "a.vcf", "A", [("E1", 100, "C", "T", ".")])
    h = _vcf(tmp_path / "h.vcf", "A", [("E1", 100, "C", "T", _ann("rpoB", "p.Ser450Leu"))])
    lift = _map(tmp_path / "lift.tsv", [("E1", 100, 761155, "+")])

    v = qc.load_variants(_args([a], vcfs_h37rv=[h], pos_liftover=lift))["A"]["E1:100"]

    assert v["pos_h37rv"] == "H37Rv:761155", "the lift, not the in-place record's own position"
    assert v["aa_h37rv"] == "p.Ser450Leu"


def test_an_in_place_record_pairs_by_contig_as_well_as_position(tmp_path):
    a = _vcf(tmp_path / "a.vcf", "A", [("E1", 100, "C", "T", ".")])
    h = _vcf(tmp_path / "h.vcf", "A", [("E2", 100, "C", "T", _ann("rpoB", "p.Ser450Leu"))])

    v = qc.load_variants(_args([a], vcfs_h37rv=[h]))["A"]["E1:100"]

    assert not v.get("aa_h37rv")
