"""The report's indel view and its codon-level marks.

parse_indel_matrix turns the run's indel matrix into the matrix panel's second view, in the SNP view's
column order; parse_mnv_table and build_snp_matrix mark the SNPs that are one amino-acid change with
another one of their codon on the same reads.
"""

from __future__ import annotations

from qcreport import panels as qc_panels
from qcreport import parsers as qc_parsers

INDEL_HEADER = ["reference", "contig", "pos", "ref_allele", "alt_allele", "length", "gene", "effect", "hgvs_c",
                "hgvs_p", "n_pass", "A|AF", "A|DP", "A|FT", "B|AF", "B|DP", "B|FT"]


def write(path, lines):
    path.write_text("".join("\t".join(map(str, ln)) + "\n" for ln in lines), encoding="utf-8")
    return str(path)


def test_the_indel_view_keeps_only_called_cells_and_their_filter(tmp_path):
    p = write(tmp_path / "ind.tsv", [INDEL_HEADER,
              ["r", "c1", 100, "T", "TC", 1, "Rv0678", "frameshift_variant", "c.144dupC", "p.Glu49fs", 1,
               "1.0000", 40, "PASS", "0.1000", 40, "LowSupport"],
              ["r", "c1", 200, "CGA", "C", -2, "g", "frameshift_variant", "", "", 1, "0.9000", 30, "PASS", "0", 50, ""]])
    m = qc_parsers.parse_indel_matrix(p)
    assert m["samples"] == ["A", "B"]
    first, second = m["rows"]
    assert (first["gene"], first["len"], first["aa"], first["hgvs_c"]) == ("Rv0678", 1, "p.Glu49fs", "c.144dupC")
    assert first["cells"] == {0: [1.0, 40, "PASS"], 1: [0.1, 40, "LowSupport"]}
    assert second["cells"] == {0: [0.9, 30, "PASS"]}, "AF 0 is read without the indel, not a call"
    assert first["n"] == 2 and second["n"] == 1


def test_the_indel_view_follows_the_snp_view_column_order(tmp_path):
    p = write(tmp_path / "ind.tsv", [INDEL_HEADER,
              ["r", "c1", 100, "T", "TC", 1, "", "", "", "", 1, "1.0000", 40, "PASS", "1.0000", 40, "PASS"]])
    m = qc_parsers.parse_indel_matrix(p, order=["B", "Z", "A"])
    assert m["samples"] == ["B", "A"]
    assert m["rows"][0]["cells"] == {0: [1.0, 40, "PASS"], 1: [1.0, 40, "PASS"]}


def test_an_indel_matrix_without_calls_or_a_placeholder_hides_the_view(tmp_path):
    empty = write(tmp_path / "ind.tsv", [INDEL_HEADER])
    assert qc_parsers.parse_indel_matrix(empty) is None
    ph = tmp_path / "NO_FILE_INDELS"
    ph.write_text("")
    assert qc_parsers.parse_indel_matrix(str(ph)) is None
    assert qc_parsers.parse_indel_matrix(None) is None


def test_the_codon_table_maps_every_position_of_a_codon(tmp_path):
    p = write(tmp_path / "mnv.tsv", [["sample", "reference", "contig", "gene", "positions", "aa_change",
                                      "snp_aa_changes", "consequence_shift", "mnv_frequency"],
                                     ["S1", "r", "c1", "rpoB", "139,141", "His445Glu", "His445Asp, His445Gln",
                                      "Concordant", "0.9904"]])
    m = qc_parsers.parse_mnv_table(p)
    entry = ["His445Glu", 0.9904, "His445Asp, His445Gln", "Concordant"]
    assert m == {"S1": {"c1:139": entry, "c1:141": entry}}
    assert qc_parsers.parse_mnv_table(None) == {}


def test_the_snp_view_marks_only_the_sample_whose_reads_carry_the_codon_change():
    variants = {"S1": {"c1:139": {"ref": "C", "alt": "G", "af": 1.0, "dp": 40, "aa": "p.His445Asp"}},
                "S2": {"c1:139": {"ref": "C", "alt": "G", "af": 1.0, "dp": 38, "aa": "p.His445Asp"}}}
    mnv = {"S1": {"c1:139": ["His445Glu", 0.99, "His445Asp, His445Gln", "Concordant"]}}
    row = qc_panels.build_snp_matrix(variants, mnv=mnv)["rows"][0]
    assert row["mnv"] == {0: ["His445Glu", 0.99, "His445Asp, His445Gln", "Concordant"]}
    assert row["mnv_aa"] == ["His445Glu"]


def test_a_snp_view_without_codon_changes_carries_no_marks():
    variants = {"S1": {"c1:139": {"ref": "C", "alt": "G", "af": 1.0, "dp": 40}}}
    row = qc_panels.build_snp_matrix(variants)["rows"][0]
    assert "mnv" not in row and "mnv_aa" not in row


def test_an_mnp_record_is_every_snp_it_changes(tmp_path):
    """FreeBayes writes two changes of one codon on the same reads as one record (CAC>GAG): each base
    it changes is a SNP at its own position, with the record's fraction and its codon-level annotation."""
    vcf = tmp_path / "s.vcf"
    vcf.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
                   "c1\t139\t.\tCAC\tGAG\t50\t.\tANN=GAG|missense_variant|MODERATE|rpoB|g|transcript|t|protein_coding|1/1|c.1333_1335delCACinsGAG|"
                   "p.His445Glu|||||\tGT:AD:DP\t1/1:1,103:104\n"
                   "c1\t200\t.\tAT\tA\t50\t.\t.\tGT:AD:DP\t1/1:0,40:40\n")
    v = qc_parsers.parse_vcfs([str(vcf)])["S1"]
    assert sorted(v) == ["c1:139", "c1:141"], "the unchanged middle base is no SNP, and the indel is not one"
    assert (v["c1:139"]["ref"], v["c1:139"]["alt"], v["c1:141"]["ref"], v["c1:141"]["alt"]) == ("C", "G", "C", "G")
    assert v["c1:139"]["aa"] == v["c1:141"]["aa"] == "p.His445Glu"
    assert v["c1:141"]["af"] == round(103 / 104, 4)
