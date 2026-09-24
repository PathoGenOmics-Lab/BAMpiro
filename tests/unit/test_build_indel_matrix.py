"""Unit tests for bin/build_indel_matrix.py.

The run's indels in one wide matrix: one row per indel that at least one sample calls with a PASS,
three columns per sample (AF, DP, FT). A call below the rule keeps its fraction, marked LowSupport;
a sample without a call reads AF 0 at the depth of the indel's anchor base where it was read.
"""

from __future__ import annotations

import subprocess
import sys

from conftest import BIN, load_script

load_script("build_snp_matrix")          # build_indel_matrix reads VCFs with its helpers
bim = load_script("build_indel_matrix")

HEADER = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=chr1,length=1000>\n"
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n"
    "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"d\">\n"
    "##FORMAT=<ID=RO,Number=1,Type=Integer,Description=\"r\">\n"
    "##FORMAT=<ID=AO,Number=A,Type=Integer,Description=\"a\">\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n")
ANN = ("ANN=TC|frameshift_variant|HIGH|Rv0678|gene-Rv0678|transcript|t|protein_coding|1/1|"
       "c.144dupC|p.Glu49fs|145/498|145/498|49/165||")


def _vcf(path, sample, rows):
    """rows: (pos, ref, alt, filter, ro, ao, info)."""
    body = "".join(f"chr1\t{pos}\t.\t{ref}\t{alt}\t50\t{flt}\t{info or '.'}\tGT:DP:RO:AO\t0/1:{ro + ao}:{ro}:{ao}\n"
                   for pos, ref, alt, flt, ro, ao, info in rows)
    path.write_text(HEADER.format(sample=sample) + body)
    return path


def _allpos(path, sample, depth):
    """An all-positions VCF with the given depth at every position 1..300."""
    body = "".join(f"chr1\t{p}\t.\tA\t.\t.\t.\t.\tGT:DP:AD\t0:{depth}:{depth},0\n" for p in range(1, 301))
    path.write_text(HEADER.format(sample=sample).replace(
        "##FORMAT=<ID=RO", "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"a\">\n##FORMAT=<ID=RO") + body)
    return path


def _run(tmp_path, vcfs, depth_vcfs=()):
    out = tmp_path / "indel_matrix.tsv"
    args = [sys.executable, str(BIN / "build_indel_matrix.py"), "--vcfs", *map(str, vcfs), "--reference", "ref",
            "-o", str(out)]
    if depth_vcfs:
        args += ["--depth-vcfs", *map(str, depth_vcfs)]
    proc = subprocess.run(args, capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    lines = out.read_text().splitlines()
    header = lines[0].split("\t")
    return header, [dict(zip(header, ln.split("\t"))) for ln in lines[1:]], proc.stderr


def test_a_row_needs_a_pass_in_at_least_one_sample(tmp_path):
    a = _vcf(tmp_path / "a.vcf", "A", [(100, "T", "TC", "PASS", 0, 40, ANN),
                                       (200, "CG", "C", "LowSupport", 38, 2, None)])
    header, rows, err = _run(tmp_path, [a])
    assert [r["pos"] for r in rows] == ["100"], "an indel no sample calls with a PASS is left out"
    assert "1 more only below the rule" in err


def test_a_call_below_the_rule_keeps_its_fraction_marked_low_support(tmp_path):
    """A minority frameshift starts below the rule: blanking it would say the sample does not have it."""
    a = _vcf(tmp_path / "a.vcf", "A", [(100, "T", "TC", "PASS", 0, 40, ANN)])
    b = _vcf(tmp_path / "b.vcf", "B", [(100, "T", "TC", "LowSupport", 36, 4, ANN)])
    _, rows, _ = _run(tmp_path, [a, b])
    assert (rows[0]["B|AF"], rows[0]["B|DP"], rows[0]["B|FT"]) == ("0.1000", "40", "LowSupport")
    assert (rows[0]["A|AF"], rows[0]["A|FT"]) == ("1.0000", "PASS")
    assert rows[0]["n_pass"] == "1"


def test_a_sample_without_the_indel_reads_zero_at_the_anchor_depth(tmp_path):
    a = _vcf(tmp_path / "a.vcf", "A", [(100, "T", "TC", "PASS", 0, 40, ANN)])
    b = _vcf(tmp_path / "b.vcf", "B", [])
    _, rows, _ = _run(tmp_path, [a, b], [_allpos(tmp_path / "a.all.vcf", "A", 40), _allpos(tmp_path / "b.all.vcf", "B", 55)])
    assert (rows[0]["B|AF"], rows[0]["B|DP"], rows[0]["B|FT"]) == ("0", "55", "")


def test_a_sample_nobody_read_there_is_not_zero(tmp_path):
    a = _vcf(tmp_path / "a.vcf", "A", [(100, "T", "TC", "PASS", 0, 40, ANN)])
    b = _vcf(tmp_path / "b.vcf", "B", [])
    _, rows, _ = _run(tmp_path, [a, b], [_allpos(tmp_path / "b.all.vcf", "B", 0)])
    assert (rows[0]["B|AF"], rows[0]["B|DP"]) == ("", "0"), "absence where nothing was read says nothing"


def test_without_depths_an_absent_cell_stays_blank(tmp_path):
    a = _vcf(tmp_path / "a.vcf", "A", [(100, "T", "TC", "PASS", 0, 40, ANN)])
    b = _vcf(tmp_path / "b.vcf", "B", [])
    _, rows, _ = _run(tmp_path, [a, b])
    assert (rows[0]["B|AF"], rows[0]["B|DP"], rows[0]["B|FT"]) == ("", "", "")


def test_snps_are_left_to_the_snp_matrix(tmp_path):
    a = _vcf(tmp_path / "a.vcf", "A", [(50, "A", "G", "PASS", 0, 40, None), (60, "AC", "GT", "PASS", 0, 40, None),
                                       (100, "T", "TC", "PASS", 0, 40, None)])
    _, rows, _ = _run(tmp_path, [a])
    assert [r["pos"] for r in rows] == ["100"]


def test_the_annotation_and_the_length_are_carried(tmp_path):
    a = _vcf(tmp_path / "a.vcf", "A", [(100, "T", "TC", "PASS", 0, 40, ANN), (200, "CGAT", "C", "PASS", 0, 40, None)])
    header, rows, _ = _run(tmp_path, [a])
    assert header[:11] == ["reference", "contig", "pos", "ref_allele", "alt_allele", "length", "gene", "effect",
                           "hgvs_c", "hgvs_p", "n_pass"]
    ins, dele = rows
    assert (ins["length"], ins["gene"], ins["effect"], ins["hgvs_c"], ins["hgvs_p"]) == \
           ("1", "Rv0678", "frameshift_variant", "c.144dupC", "p.Glu49fs")
    assert dele["length"] == "-3"


def test_an_unfiltered_vcf_reads_as_passing(tmp_path):
    """An indel VCF written without the soft filter (FILTER '.') is every call the caller kept."""
    a = _vcf(tmp_path / "a.vcf", "A", [(100, "T", "TC", ".", 0, 40, None)])
    _, rows, _ = _run(tmp_path, [a])
    assert rows[0]["A|FT"] == "PASS" and rows[0]["n_pass"] == "1"


def test_two_indels_at_one_position_are_two_rows(tmp_path):
    a = _vcf(tmp_path / "a.vcf", "A", [(100, "T", "TC", "PASS", 0, 40, None)])
    b = _vcf(tmp_path / "b.vcf", "B", [(100, "T", "TCC", "PASS", 0, 40, None)])
    _, rows, _ = _run(tmp_path, [a, b])
    assert [(r["alt_allele"], r["A|AF"], r["B|AF"]) for r in rows] == [("TC", "1.0000", ""), ("TCC", "", "1.0000")]


def test_ann_reads_gene_effect_and_both_hgvs():
    assert bim._ann(ANN) == ("Rv0678", "frameshift_variant", "c.144dupC", "p.Glu49fs")
    assert bim._ann("DP=10") == ("", "", "", "")
