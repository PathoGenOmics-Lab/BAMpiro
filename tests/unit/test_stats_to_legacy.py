"""Unit tests for bin/stats_to_legacy.py.

The script turns fastp JSON, samtools stats, a VCF and an optional pathotypr
report into the flat legacy `<sample>.log`. Almost every parser here is written
to swallow exceptions and hand back a zeroed default, so a large part of these
tests exists to PIN that "no input file" never propagates an exception.
"""

from __future__ import annotations

import gzip
import subprocess
import sys

import pytest

from conftest import load_script

stl = load_script("stats_to_legacy")


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------

# snpEff ANN has a fixed 16-subfield layout; index 15 is ERRORS/WARNINGS/INFO.
_ANN_TEMPLATE = [
    "G", "missense_variant", "MODERATE", "geneA", "GENE_ID", "transcript",
    "TX_ID", "protein_coding", "1/1", "c.100A>G", "p.Lys34Arg", "100/1000",
    "100/900", "34/300", "", "",
]


def _ann_entry(allele="G", effect="missense_variant", impact="MODERATE", gene="geneA", errors=""):
    sub = list(_ANN_TEMPLATE)
    sub[0], sub[1], sub[2], sub[3], sub[15] = allele, effect, impact, gene, errors
    return "|".join(sub)


def _write_fai(tmp_path, rows, name="ref.fa.fai"):
    p = tmp_path / name
    p.write_text("".join(f"{contig}\t{length}\t0\t60\t61\n" for contig, length in rows))
    return p


_VCF_HEADER = (
    "##fileformat=VCFv4.2\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
)


def _rec(pos, ref, alt, gt, info="DP=10", contig="chr1", fmt="GT:DP", sample_extra="10"):
    sample = gt if fmt == "GT" else f"{gt}:{sample_extra}"
    return "\t".join([contig, str(pos), ".", ref, alt, "50", "PASS", info, fmt, sample])


def _write_vcf(tmp_path, records, name="calls.vcf", header=_VCF_HEADER):
    p = tmp_path / name
    p.write_text(header + "".join(r + "\n" for r in records))
    return p


# --------------------------------------------------------------------------
# _effect_class
# --------------------------------------------------------------------------

_EFFECT_CASES = [
    ("frameshift_variant", "frameshift"),
    # ORDER MATTERS: 'frameshift' is the first test in the cascade, so a term
    # that also mentions stop_gained still resolves to the frameshift class.
    ("frameshift_variant&stop_gained", "frameshift"),
    ("stop_gained", "stop_gained"),
    # stop_gained is checked before splice, so the compound term is a stop_gained.
    ("stop_gained&splice_region_variant", "stop_gained"),
    ("stop_lost", "stop_lost"),
    ("stop_lost&splice_region_variant", "stop_lost"),
    ("start_lost", "start_lost"),
    ("missense_variant", "missense"),
    # missense is checked before splice for the same reason.
    ("missense_variant&splice_region_variant", "missense"),
    ("stop_retained_variant", "synonymous"),
    ("start_retained_variant", "synonymous"),
    ("synonymous_variant", "synonymous"),
    ("inframe_insertion", "inframe_indel"),
    ("inframe_deletion", "inframe_indel"),
    ("disruptive_inframe_deletion", "inframe_indel"),
    ("conservative_inframe_insertion", "inframe_indel"),
    ("splice_acceptor_variant", "splice"),
    ("splice_donor_variant", "splice"),
    ("splice_region_variant", "splice"),
    ("intergenic_region", "intergenic"),
    ("intragenic_variant", "intergenic"),
    ("upstream_gene_variant", "regulatory"),
    ("downstream_gene_variant", "regulatory"),
    ("5_prime_utr_variant", "regulatory"),
    ("3_prime_utr_premature_start_codon_gain_variant", "regulatory"),
    ("regulatory_region_variant", "regulatory"),
    ("tf_binding_site_variant", "regulatory"),
    ("non_coding_transcript_exon_variant", "other"),
    ("sequence_feature", "other"),
    ("", "other"),
]


@pytest.mark.parametrize("effect,expected", _EFFECT_CASES)
def test_effect_class(effect, expected):
    assert stl._effect_class(effect) == expected


def test_effect_class_cases_cover_every_declared_class():
    """The parametrised table above must exercise all of _EFF_CLASSES."""
    seen = {expected for _, expected in _EFFECT_CASES}
    assert seen == set(stl._EFF_CLASSES)


def test_effect_class_is_substring_based_not_exact():
    """The cascade uses `in`, so any term containing the token is classified."""
    assert stl._effect_class("some_prefix_missense_thing") == "missense"


# --------------------------------------------------------------------------
# _parse_ann_field
# --------------------------------------------------------------------------

def test_parse_ann_field_no_ann_returns_none():
    assert stl._parse_ann_field("DP=10;AF=0.5;MQ=60") is None


def test_parse_ann_field_empty_info_returns_none():
    assert stl._parse_ann_field("") is None


def test_parse_ann_field_empty_ann_value_returns_none():
    """`ANN=` with nothing after it is falsy and treated as un-annotated."""
    assert stl._parse_ann_field("DP=10;ANN=") is None


def test_parse_ann_field_too_few_subfields_returns_none():
    """Fewer subfields than the impact index (2) means malformed -> unannotated."""
    assert stl._parse_ann_field("ANN=G|missense_variant") is None
    assert stl._parse_ann_field("ANN=G") is None


def test_parse_ann_field_exactly_three_subfields_is_usable():
    """Three subfields is the minimum; gene and errors default to empty."""
    cls, impact, gene, has_error, db_error = stl._parse_ann_field("ANN=G|missense_variant|HIGH")
    assert (cls, impact, gene, has_error, db_error) == ("missense", "HIGH", "", False, False)


def test_parse_ann_field_full_entry():
    info = f"DP=10;ANN={_ann_entry(effect='stop_gained', impact='HIGH', gene='rpoB')};MQ=60"
    assert stl._parse_ann_field(info) == ("stop_gained", "HIGH", "rpoB", False, False)


def test_parse_ann_field_uses_first_ann_entry_only():
    """snpEff pre-sorts entries by severity, so only the first one is read."""
    info = "ANN=" + ",".join([
        _ann_entry(effect="stop_gained", impact="HIGH", gene="geneA"),
        _ann_entry(effect="synonymous_variant", impact="LOW", gene="geneB"),
    ])
    assert stl._parse_ann_field(info) == ("stop_gained", "HIGH", "geneA", False, False)


def test_parse_ann_field_splits_effect_on_ampersand():
    info = f"ANN={_ann_entry(effect='frameshift_variant&stop_gained', impact='HIGH', gene='katG')}"
    assert stl._parse_ann_field(info)[0] == "frameshift"


@pytest.mark.parametrize("raw", ["high", " HIGH ", "HiGh"])
def test_parse_ann_field_normalises_impact_case_and_space(raw):
    assert stl._parse_ann_field(f"ANN={_ann_entry(impact=raw)}")[1] == "HIGH"


@pytest.mark.parametrize("raw", ["BOGUS", "", "moderate_ish"])
def test_parse_ann_field_unknown_impact_is_coerced_to_modifier(raw):
    """Anything outside HIGH/MODERATE/LOW/MODIFIER falls back to MODIFIER."""
    assert stl._parse_ann_field(f"ANN={_ann_entry(impact=raw)}")[1] == "MODIFIER"


def test_parse_ann_field_flags_warnings():
    info = f"ANN={_ann_entry(errors='WARNING_TRANSCRIPT_NO_START_CODON')}"
    _, _, _, has_error, db_error = stl._parse_ann_field(info)
    assert has_error is True
    assert db_error is False


@pytest.mark.parametrize("token", ["ERROR_CHROMOSOME_NOT_FOUND", "ERROR_OUT_OF_CHROMOSOME_RANGE"])
def test_parse_ann_field_flags_db_errors(token):
    _, _, _, has_error, db_error = stl._parse_ann_field(f"ANN={_ann_entry(errors=token)}")
    assert has_error is True
    assert db_error is True


def test_parse_ann_field_whitespace_only_errors_is_not_an_error():
    _, _, _, has_error, _ = stl._parse_ann_field(f"ANN={_ann_entry(errors='   ')}")
    assert has_error is False


# --------------------------------------------------------------------------
# _format_gene_burden
# --------------------------------------------------------------------------

def test_format_gene_burden_empty_is_na():
    assert stl._format_gene_burden({}) == "NA"


def test_format_gene_burden_no_high_or_moderate_is_na():
    """Genes carrying only LOW/MODIFIER variants are dropped, leaving nothing."""
    agg = {"geneA": [0, 0, {"synonymous": 5}], "geneB": [0, 0, {"regulatory": 3}]}
    assert stl._format_gene_burden(agg) == "NA"


def test_format_gene_burden_format_and_sort_order():
    agg = {
        "geneLow": [0, 1, {"missense": 1}],
        "geneTop": [3, 0, {"stop_gained": 2, "frameshift": 1}],
        "geneMid": [1, 9, {"missense": 4}],
        "geneSkip": [0, 0, {"synonymous": 7}],
    }
    out = stl._format_gene_burden(agg)
    assert out == "geneTop:3:0:stop_gained;geneMid:1:9:missense;geneLow:0:1:missense"


def test_format_gene_burden_ties_break_on_gene_name():
    agg = {"zeta": [2, 1, {"missense": 1}], "alpha": [2, 1, {"missense": 1}]}
    assert stl._format_gene_burden(agg) == "alpha:2:1:missense;zeta:2:1:missense"


def test_format_gene_burden_empty_effect_map_falls_back_to_other():
    assert stl._format_gene_burden({"geneA": [1, 0, {}]}) == "geneA:1:0:other"


def test_format_gene_burden_sanitises_delimiters_in_gene_names():
    agg = {"we;ird:gene\tname": [1, 0, {"missense": 1}]}
    assert stl._format_gene_burden(agg) == "we_ird_gene_name:1:0:missense"


def test_format_gene_burden_is_capped_at_40_genes():
    # Distinct HIGH counts make the ordering fully deterministic.
    agg = {f"gene{i:03d}": [100 - i, 0, {"missense": 1}] for i in range(45)}
    out = stl._format_gene_burden(agg)
    items = out.split(";")
    assert stl._GENE_BURDEN_CAP == 40
    assert len(items) == 40
    assert items[0] == "gene000:100:0:missense"
    assert items[-1] == "gene039:61:0:missense"


# --------------------------------------------------------------------------
# get_genome_size / get_contig_offsets
# --------------------------------------------------------------------------

def test_get_genome_size_sums_second_column(tmp_path):
    fai = _write_fai(tmp_path, [("chr1", 1000), ("chr2", 500)])
    assert stl.get_genome_size(str(fai)) == 1500


def test_get_genome_size_missing_file_returns_zero(tmp_path):
    """Never raises: a missing .fai degrades to an unknown (0) genome size."""
    assert stl.get_genome_size(str(tmp_path / "absent.fai")) == 0


def test_get_genome_size_skips_short_lines(tmp_path):
    fai = tmp_path / "ref.fai"
    fai.write_text("chr1\t1000\t0\t60\t61\nnotabhere\n\nchr2\t7\t0\t60\t61\n")
    assert stl.get_genome_size(str(fai)) == 1007


def test_get_genome_size_non_numeric_length_aborts_and_returns_partial(tmp_path):
    """The try/except wraps the whole loop, so a bad row truncates the total
    instead of skipping just that row. Pinned as observed behaviour."""
    fai = tmp_path / "ref.fai"
    fai.write_text("chr1\t1000\t0\t60\t61\nchr2\tNOT_A_NUMBER\t0\t60\t61\nchr3\t99\t0\t60\t61\n")
    assert stl.get_genome_size(str(fai)) == 1000


def test_get_contig_offsets_is_cumulative(tmp_path):
    fai = _write_fai(tmp_path, [("chr1", 1000), ("chr2", 500), ("chr3", 250)])
    assert stl.get_contig_offsets(str(fai)) == {"chr1": 0, "chr2": 1000, "chr3": 1500}


def test_get_contig_offsets_missing_file_returns_empty(tmp_path):
    assert stl.get_contig_offsets(str(tmp_path / "absent.fai")) == {}


# --------------------------------------------------------------------------
# get_bam_stats
# --------------------------------------------------------------------------

_SN_BLOCK = """# This file was produced by samtools stats
CHK\t0000\t0000\t0000
SN\traw total sequences:\t1000\t# excluding supplementary
SN\treads mapped:\t900
SN\treads mapped and paired:\t880
SN\treads properly paired:\t850
SN\treads duplicated:\t50
SN\treads unmapped:\t100
SN\treads MQ0:\t10
SN\taverage quality:\t35.5
SN\tinsert size average:\t300.0
SN\tinsert size standard deviation:\t50.0
SN\terror rate:\t0.001\t# mismatches / bases mapped
SN\tbases mapped (cigar):\t90000
SN\tsome unknown key:\t7
SN\ttruncated line
"""


def test_get_bam_stats_parses_sn_block(tmp_path):
    p = tmp_path / "stats.txt"
    p.write_text(_SN_BLOCK)
    got = stl.get_bam_stats(str(p))
    assert got == {
        "total": 1000,
        "mapped": 900,
        "mapped_paired": 880,
        "properly_paired": 850,
        "duplicates": 50,
        "avg_quality": 35.5,
        "insert_size_avg": 300.0,
        "insert_size_sd": 50.0,
        "bases_mapped": 90000,
        "reads_unmapped": 100,
        "reads_mq0": 10,
        "error_rate": 0.001,
    }


def test_get_bam_stats_missing_file_returns_zeroed_defaults(tmp_path):
    got = stl.get_bam_stats(str(tmp_path / "absent.txt"))
    assert got["total"] == 0
    assert got["mapped"] == 0
    assert got["error_rate"] == 0.0
    assert set(got) == {
        "total", "mapped", "mapped_paired", "properly_paired", "duplicates",
        "avg_quality", "insert_size_avg", "insert_size_sd", "bases_mapped",
        "reads_unmapped", "reads_mq0", "error_rate",
    }


def test_get_bam_stats_ignores_non_sn_lines(tmp_path):
    p = tmp_path / "stats.txt"
    p.write_text("FFQ\t1\t2\t3\nCOV\t[1-1]\t1\t500\nSN\treads mapped:\t7\n")
    assert stl.get_bam_stats(str(p))["mapped"] == 7


# --------------------------------------------------------------------------
# calc_breadth_of_coverage
# --------------------------------------------------------------------------

_COV_BLOCK = """SN\traw total sequences:\t1000
COV\t[1-1]\t1\t200
COV\t[2-2]\t2\t200
COV\t[6-6]\t6\t100
COV\t[11-11]\t11\t100
"""


def test_calc_breadth_of_coverage_fraction(tmp_path):
    p = tmp_path / "stats.txt"
    p.write_text(_COV_BLOCK)
    assert stl.calc_breadth_of_coverage(str(p), 1000, min_dp=1) == pytest.approx(0.6)


def test_calc_breadth_of_coverage_respects_min_dp(tmp_path):
    p = tmp_path / "stats.txt"
    p.write_text(_COV_BLOCK)
    # Only the >=5x bins (100 + 100 bases) qualify.
    assert stl.calc_breadth_of_coverage(str(p), 1000, min_dp=5) == pytest.approx(0.2)


def test_calc_breadth_of_coverage_counts_the_overflow_bin(tmp_path):
    """samtools emits a final '[1000<]' bin; the '<' is stripped so those deep
    bases are counted rather than dropped as an unparsable range."""
    p = tmp_path / "stats.txt"
    p.write_text("COV\t[1-1]\t1\t100\nCOV\t[1000<]\t1000\t400\n")
    assert stl.calc_breadth_of_coverage(str(p), 1000, min_dp=1) == pytest.approx(0.5)
    assert stl.calc_breadth_of_coverage(str(p), 1000, min_dp=500) == pytest.approx(0.4)


def test_calc_breadth_of_coverage_zero_genome_short_circuits(tmp_path, capsys):
    """genome_size <= 0 returns before the file is even opened, so a missing
    file produces no warning and no division."""
    assert stl.calc_breadth_of_coverage(str(tmp_path / "absent.txt"), 0) == 0.0
    assert stl.calc_breadth_of_coverage(str(tmp_path / "absent.txt"), -5) == 0.0
    assert capsys.readouterr().err == ""


def test_calc_breadth_of_coverage_missing_file_returns_zero(tmp_path, capsys):
    assert stl.calc_breadth_of_coverage(str(tmp_path / "absent.txt"), 1000) == 0.0
    assert "WARN calc_breadth_of_coverage" in capsys.readouterr().err


def test_calc_breadth_of_coverage_skips_unparsable_ranges(tmp_path):
    p = tmp_path / "stats.txt"
    p.write_text("COV\t[junk]\tx\t900\nCOV\t[1-1]\t1\t100\nCOV\t[2-2]\n")
    assert stl.calc_breadth_of_coverage(str(p), 1000, min_dp=1) == pytest.approx(0.1)


# --------------------------------------------------------------------------
# cov_histogram_stats
# --------------------------------------------------------------------------

def test_cov_histogram_stats_zero_genome_short_circuits(tmp_path, capsys):
    """A genome size of 0 returns the zeroed default before any file read or
    division by the (zero) denominator."""
    out = stl.cov_histogram_stats(str(tmp_path / "absent.txt"), 0)
    assert out == {"median_depth": 0.0, "breadth5": 0.0, "breadth10": 0.0, "evenness": 0.0}
    assert capsys.readouterr().err == ""


def test_cov_histogram_stats_missing_file_returns_zeroed_defaults(tmp_path, capsys):
    out = stl.cov_histogram_stats(str(tmp_path / "absent.txt"), 1000)
    assert out == {"median_depth": 0.0, "breadth5": 0.0, "breadth10": 0.0, "evenness": 0.0}
    assert "WARN cov_histogram_stats" in capsys.readouterr().err


def test_cov_histogram_stats_full_calculation(tmp_path):
    p = tmp_path / "stats.txt"
    p.write_text(_COV_BLOCK)
    out = stl.cov_histogram_stats(str(p), 1000)
    # covered = 600 bases, zero-depth = 400, total = 1000
    assert out["breadth5"] == pytest.approx(20.0)    # 200 bases at depth >= 5
    assert out["breadth10"] == pytest.approx(10.0)   # 100 bases at depth >= 10
    # half = 500; the running count starts at 400 (zero bin) and crosses at [1-1]
    assert out["median_depth"] == pytest.approx(1.0)
    # mean = 2300/1000 = 2.3; var = 16700/1000 - 2.3^2 = 11.41
    assert out["evenness"] == pytest.approx((11.41 ** 0.5) / 2.3)


def test_cov_histogram_stats_uses_bin_midpoint(tmp_path):
    p = tmp_path / "stats.txt"
    p.write_text("COV\t[10-20]\t10\t1000\n")
    out = stl.cov_histogram_stats(str(p), 1000)
    assert out["median_depth"] == pytest.approx(15.0)
    assert out["breadth10"] == pytest.approx(100.0)
    assert out["evenness"] == pytest.approx(0.0)


def test_cov_histogram_stats_counts_the_overflow_bin(tmp_path):
    p = tmp_path / "stats.txt"
    p.write_text("COV\t[1000<]\t1000\t1000\n")
    out = stl.cov_histogram_stats(str(p), 1000)
    assert out["median_depth"] == pytest.approx(1000.0)
    assert out["breadth10"] == pytest.approx(100.0)


def test_cov_histogram_stats_median_stays_zero_when_half_the_genome_is_uncovered(tmp_path):
    """The zero-depth bin is entered into the running count first. When it
    already reaches half the genome the median loop never runs and the median
    is reported as 0.0 (the true median depth IS 0 here)."""
    p = tmp_path / "stats.txt"
    p.write_text("COV\t[5-5]\t5\t400\n")
    out = stl.cov_histogram_stats(str(p), 1000)   # 600 of 1000 bases at depth 0
    assert out["median_depth"] == 0.0
    assert out["breadth5"] == pytest.approx(40.0)


def test_cov_histogram_stats_exactly_half_uncovered_also_stays_zero(tmp_path):
    """Boundary: run == half is not < half, so the loop is skipped."""
    p = tmp_path / "stats.txt"
    p.write_text("COV\t[5-5]\t5\t500\n")
    out = stl.cov_histogram_stats(str(p), 1000)
    assert out["median_depth"] == 0.0


def test_cov_histogram_stats_skips_unparsable_and_truncated_bins(tmp_path):
    p = tmp_path / "stats.txt"
    p.write_text("COV\t[junk-bin]\tx\t900\nCOV\t[9-9]\n\nCOV\t[4-4]\t4\t1000\n")
    out = stl.cov_histogram_stats(str(p), 1000)
    assert out["median_depth"] == pytest.approx(4.0)


# --------------------------------------------------------------------------
# analyze_vcf
# --------------------------------------------------------------------------

def test_analyze_vcf_missing_file_returns_zeroed_defaults(tmp_path, capsys):
    v = stl.analyze_vcf(str(tmp_path / "absent.vcf"), 1000, 10, {})
    assert v["snps"] == 0
    assert v["indels"] == 0
    assert v["ann_present"] is False
    assert v["snp_prof"] == [0] * 10
    assert "WARN analyze_vcf" in capsys.readouterr().err


def test_analyze_vcf_haploid_one_counts_as_homozygous(tmp_path):
    """Deliberate fix: homozygosity is `len(set(alleles)) == 1` over non-'.'
    alleles, so a haploid '1' (freebayes ploidy=1) is homozygous, not het.
    The old 'len == 2 and equal' test miscounted every haploid call."""
    vcf = _write_vcf(tmp_path, [_rec(100, "A", "G", "1")])
    v = stl.analyze_vcf(str(vcf))
    assert v["homo_total"] == 1
    assert v["het_total"] == 0
    assert v["snps"] == 1


@pytest.mark.parametrize("gt,homo", [
    ("1", True),
    ("1/1", True),
    ("1|1", True),
    ("2/2", True),
    ("1/2", False),
    ("0/1", False),
    ("0|1", False),
    ("./1", True),      # the missing allele is dropped, leaving a single '1'
    ("1/.", True),
])
def test_analyze_vcf_zygosity(tmp_path, gt, homo):
    vcf = _write_vcf(tmp_path, [_rec(100, "A", "G", gt)], name=f"z{gt.replace('/', '_').replace('|', 'p')}.vcf")
    v = stl.analyze_vcf(str(vcf))
    assert v["homo_total"] == (1 if homo else 0)
    assert v["het_total"] == (0 if homo else 1)


def test_analyze_vcf_haploid_ref_call_is_counted_as_a_homozygous_variant(tmp_path):
    """BUG PIN (reported, not fixed): the skip list only covers './.', '0/0'
    and '0|0'. A HAPLOID reference call ('0') survives, and because
    len(set(['0'])) == 1 it is counted as a homozygous SNP."""
    vcf = _write_vcf(tmp_path, [_rec(100, "A", "G", "0")])
    v = stl.analyze_vcf(str(vcf))
    assert v["snps"] == 1
    assert v["homo_total"] == 1


def test_analyze_vcf_haploid_missing_call_is_counted_as_heterozygous(tmp_path):
    """BUG PIN (reported, not fixed): a haploid no-call '.' is not in the skip
    list; it yields an empty allele list, so is_homo is False and the record
    lands in the heterozygous bucket."""
    vcf = _write_vcf(tmp_path, [_rec(100, "A", "G", ".")])
    v = stl.analyze_vcf(str(vcf))
    assert v["snps"] == 1
    assert v["het_total"] == 1


@pytest.mark.parametrize("gt", ["./.", "0/0", "0|0"])
def test_analyze_vcf_skips_non_calls_and_ref_calls(tmp_path, gt):
    vcf = _write_vcf(tmp_path, [_rec(100, "A", "G", gt)], name=f"s{gt.replace('/', '_').replace('|', 'p')}.vcf")
    v = stl.analyze_vcf(str(vcf))
    assert v["snps"] == 0
    assert v["ann_total"] == 0


@pytest.mark.parametrize("alt", [".", "*"])
def test_analyze_vcf_skips_placeholder_alt(tmp_path, alt):
    vcf = _write_vcf(tmp_path, [_rec(100, "A", alt, "1/1")], name=f"alt{'dot' if alt == '.' else 'star'}.vcf")
    assert stl.analyze_vcf(str(vcf))["snps"] == 0


def test_analyze_vcf_skips_records_with_fewer_than_ten_columns(tmp_path):
    p = tmp_path / "short.vcf"
    p.write_text(_VCF_HEADER + "chr1\t100\t.\tA\tG\t50\tPASS\tDP=10\tGT:DP\n")
    assert stl.analyze_vcf(str(p))["snps"] == 0


def test_analyze_vcf_skips_records_without_gt_in_format(tmp_path):
    p = tmp_path / "nogt.vcf"
    p.write_text(_VCF_HEADER + "chr1\t100\t.\tA\tG\t50\tPASS\tDP=10\tDP:AD\t10:5,5\n")
    assert stl.analyze_vcf(str(p))["snps"] == 0


def test_analyze_vcf_skips_when_sample_column_is_shorter_than_the_gt_index(tmp_path):
    p = tmp_path / "shortsample.vcf"
    p.write_text(_VCF_HEADER + "chr1\t100\t.\tA\tG\t50\tPASS\tDP=10\tDP:GT\t10\n")
    assert stl.analyze_vcf(str(p))["snps"] == 0


def test_analyze_vcf_snp_indel_split_and_titv(tmp_path):
    vcf = _write_vcf(tmp_path, [
        _rec(100, "A", "G", "1/1"),     # transition
        _rec(200, "C", "T", "1/1"),     # transition
        _rec(300, "G", "T", "0/1"),     # transversion, heterozygous
        _rec(400, "AT", "A", "1/1"),    # homozygous deletion
        _rec(500, "A", "ATT", "0/1"),   # heterozygous insertion
    ])
    v = stl.analyze_vcf(str(vcf))
    assert v["snps"] == 3
    assert v["indels"] == 2
    assert v["homo_total"] == 3
    assert v["het_total"] == 2
    assert v["homo_indels"] == 1
    assert v["ti"] == 2
    assert v["tv"] == 1


def test_analyze_vcf_titv_ignores_non_acgt_and_identical_alleles(tmp_path):
    vcf = _write_vcf(tmp_path, [
        _rec(100, "N", "G", "1/1"),
        _rec(200, "A", "N", "1/1"),
        _rec(300, "A", "A", "1/1"),
    ])
    v = stl.analyze_vcf(str(vcf))
    assert v["snps"] == 3
    assert (v["ti"], v["tv"]) == (0, 0)


def test_analyze_vcf_binning_formula(tmp_path):
    """b = min(nbins - 1, (contig_offset + POS) * nbins // genome_size)."""
    vcf = _write_vcf(tmp_path, [
        _rec(1, "A", "G", "1/1"),      # 1*10//1000   -> bin 0
        _rec(500, "A", "G", "1/1"),    # 500*10//1000 -> bin 5
        _rec(1000, "A", "G", "1/1"),   # 1000*10//1000 = 10 -> clamped to 9
    ])
    v = stl.analyze_vcf(str(vcf), g_size=1000, nbins=10, contig_offsets={"chr1": 0})
    assert v["snp_prof"] == [1, 0, 0, 0, 0, 1, 0, 0, 0, 1]


def test_analyze_vcf_binning_applies_contig_offsets(tmp_path):
    """VCF POS restarts at 1 per contig, so the cumulative offset places the
    variant on a single whole-genome axis."""
    vcf = _write_vcf(tmp_path, [
        _rec(500, "A", "G", "1/1", contig="chr1"),   # gpos 500  -> 500*10//2000  = 2
        _rec(500, "A", "G", "1/1", contig="chr2"),   # gpos 1500 -> 1500*10//2000 = 7
    ])
    v = stl.analyze_vcf(str(vcf), g_size=2000, nbins=10, contig_offsets={"chr1": 0, "chr2": 1000})
    assert v["snp_prof"] == [0, 0, 1, 0, 0, 0, 0, 1, 0, 0]


def test_analyze_vcf_unknown_contig_offsets_to_zero(tmp_path):
    vcf = _write_vcf(tmp_path, [_rec(500, "A", "G", "1/1", contig="chrUnknown")])
    v = stl.analyze_vcf(str(vcf), g_size=1000, nbins=10, contig_offsets={"chr1": 0})
    assert v["snp_prof"][5] == 1


def test_analyze_vcf_no_binning_when_genome_size_is_unknown(tmp_path):
    vcf = _write_vcf(tmp_path, [_rec(500, "A", "G", "1/1"), _rec(600, "AT", "A", "0/1")])
    v = stl.analyze_vcf(str(vcf), g_size=0, nbins=10)
    assert v["snps"] == 1
    assert v["snp_prof"] == [0] * 10
    assert v["het_prof"] == [0] * 10
    assert v["indel_prof"] == [0] * 10


def test_analyze_vcf_non_numeric_position_is_not_binned(tmp_path):
    p = tmp_path / "badpos.vcf"
    p.write_text(_VCF_HEADER + "chr1\tNOT_A_POS\t.\tA\tG\t50\tPASS\tDP=10\tGT:DP\t1/1:10\n")
    v = stl.analyze_vcf(str(p), g_size=1000, nbins=10)
    assert v["snps"] == 1
    assert v["snp_prof"] == [0] * 10


def test_analyze_vcf_only_homozygous_snps_reach_the_snp_profile(tmp_path):
    """Heterozygous SNPs go to het_prof only; indels go to indel_prof
    regardless of zygosity, and a het indel is additionally in het_prof."""
    vcf = _write_vcf(tmp_path, [
        _rec(100, "A", "G", "1/1"),     # homozygous SNP  -> snp_prof[1]
        _rec(300, "A", "G", "0/1"),     # het SNP         -> het_prof[3] only
        _rec(500, "AT", "A", "1/1"),    # homozygous indel-> indel_prof[5] only
        _rec(700, "AT", "A", "0/1"),    # het indel       -> indel_prof[7] + het_prof[7]
    ])
    v = stl.analyze_vcf(str(vcf), g_size=1000, nbins=10, contig_offsets={"chr1": 0})
    assert v["snp_prof"] == [0, 1, 0, 0, 0, 0, 0, 0, 0, 0]
    assert v["het_prof"] == [0, 0, 0, 1, 0, 0, 0, 1, 0, 0]
    assert v["indel_prof"] == [0, 0, 0, 0, 0, 1, 0, 1, 0, 0]


def test_analyze_vcf_reads_gzipped_input(tmp_path):
    p = tmp_path / "calls.vcf.gz"
    with gzip.open(p, "wt") as fh:
        fh.write(_VCF_HEADER + _rec(100, "A", "G", "1/1") + "\n")
    assert stl.analyze_vcf(str(p))["snps"] == 1


def test_analyze_vcf_without_ann_leaves_annotation_absent(tmp_path):
    vcf = _write_vcf(tmp_path, [_rec(100, "A", "G", "1/1", info="DP=10")])
    v = stl.analyze_vcf(str(vcf))
    assert v["ann_present"] is False
    assert v["ann_total"] == 1        # counted as a candidate record...
    assert v["ann_records"] == 0      # ...but not as an annotated one


def test_analyze_vcf_accumulates_annotation_counts(tmp_path):
    vcf = _write_vcf(tmp_path, [
        _rec(100, "A", "G", "1/1", info="ANN=" + _ann_entry(effect="missense_variant", impact="MODERATE", gene="geneA")),
        _rec(200, "C", "T", "1/1", info="ANN=" + _ann_entry(effect="synonymous_variant", impact="LOW", gene="geneA")),
        _rec(300, "G", "T", "0/1", info="ANN=" + _ann_entry(effect="stop_gained", impact="HIGH", gene="geneB",
                                                            errors="WARNING_TRANSCRIPT_INCOMPLETE")),
        _rec(400, "AT", "A", "1/1", info="DP=10"),
    ])
    v = stl.analyze_vcf(str(vcf))
    assert v["ann_present"] is True
    assert v["ann_total"] == 4
    assert v["ann_records"] == 3
    assert v["ann_impact"] == {"HIGH": 1, "MODERATE": 1, "LOW": 1, "MODIFIER": 0}
    assert v["ann_class"]["missense"] == 1
    assert v["ann_class"]["synonymous"] == 1
    assert v["ann_class"]["stop_gained"] == 1
    assert v["ann_warn"] == 1
    assert v["ann_db_error"] is False
    assert v["ann_gene"] == {
        "geneA": [0, 1, {"missense": 1, "synonymous": 1}],
        "geneB": [1, 0, {"stop_gained": 1}],
    }


def test_analyze_vcf_flags_db_error(tmp_path):
    vcf = _write_vcf(tmp_path, [
        _rec(100, "A", "G", "1/1", info="ANN=" + _ann_entry(errors="ERROR_CHROMOSOME_NOT_FOUND")),
    ])
    assert stl.analyze_vcf(str(vcf))["ann_db_error"] is True


def test_analyze_vcf_ignores_genes_with_an_empty_name(tmp_path):
    vcf = _write_vcf(tmp_path, [_rec(100, "A", "G", "1/1", info="ANN=G|intergenic_region|MODIFIER|")])
    v = stl.analyze_vcf(str(vcf))
    assert v["ann_records"] == 1
    assert v["ann_gene"] == {}


# --------------------------------------------------------------------------
# get_pathotypr_data
# --------------------------------------------------------------------------

_PATHO_TSV = "genome\tlineage:count\tmajor_lineage\nS1\tL4:10;L4.3:5\tL4.3\n"


def test_get_pathotypr_data_parses_report(tmp_path):
    p = tmp_path / "lineage.tsv"
    p.write_text(_PATHO_TSV)
    assert stl.get_pathotypr_data(str(p)) == {"major": "L4.3", "counts": "L4:10;L4.3:5"}


@pytest.mark.parametrize("path", [None, "", "some/NO_LINEAGE_FILE/report.tsv"])
def test_get_pathotypr_data_sentinel_paths_return_na(path):
    assert stl.get_pathotypr_data(path) == {"major": "NA", "counts": "NA"}


def test_get_pathotypr_data_missing_file_returns_na(tmp_path):
    assert stl.get_pathotypr_data(str(tmp_path / "absent.tsv")) == {"major": "NA", "counts": "NA"}


def test_get_pathotypr_data_empty_file_returns_na(tmp_path):
    p = tmp_path / "empty.tsv"
    p.write_text("")
    assert stl.get_pathotypr_data(str(p)) == {"major": "NA", "counts": "NA"}


def test_get_pathotypr_data_header_only_returns_na(tmp_path):
    p = tmp_path / "hdr.tsv"
    p.write_text("genome\tlineage:count\tmajor_lineage\n")
    assert stl.get_pathotypr_data(str(p)) == {"major": "NA", "counts": "NA"}


def test_get_pathotypr_data_unknown_columns_return_na(tmp_path):
    p = tmp_path / "other.tsv"
    p.write_text("genome\tsomething\tsomething_else\nS1\ta\tb\n")
    assert stl.get_pathotypr_data(str(p)) == {"major": "NA", "counts": "NA"}


def test_get_pathotypr_data_partial_columns(tmp_path):
    """Only the columns actually present in the header are filled in."""
    p = tmp_path / "partial.tsv"
    p.write_text("genome\tmajor_lineage\nS1\tL2.2.1\n")
    assert stl.get_pathotypr_data(str(p)) == {"major": "L2.2.1", "counts": "NA"}


def test_get_pathotypr_data_short_data_row_returns_na(tmp_path):
    p = tmp_path / "short.tsv"
    p.write_text("genome\tlineage:count\tmajor_lineage\nS1\n")
    assert stl.get_pathotypr_data(str(p)) == {"major": "NA", "counts": "NA"}


def test_get_pathotypr_data_unreadable_path_returns_na(tmp_path, capsys):
    """A directory passes the os.path.exists() check and then fails to open;
    the exception is swallowed and the defaults are returned."""
    d = tmp_path / "a_directory"
    d.mkdir()
    assert stl.get_pathotypr_data(str(d)) == {"major": "NA", "counts": "NA"}
    assert "Error parsing Pathotypr file" in capsys.readouterr().err


# --------------------------------------------------------------------------
# end to end
# --------------------------------------------------------------------------

def _fastp_json():
    return (
        '{"summary": {"before_filtering": {"total_reads": 1000},'
        ' "after_filtering": {"total_reads": 900, "total_bases": 90000,'
        ' "read1_mean_length": 100, "q20_rate": 0.98, "q30_rate": 0.95,'
        ' "gc_content": 0.65}}, "duplication": {"rate": 0.1}}'
    )


def test_end_to_end_writes_the_legacy_log(tmp_path, repo_root):
    """main() has NO output flag: it writes '<sample>.log' into the CURRENT
    WORKING DIRECTORY, hence cwd=tmp_path."""
    (tmp_path / "fastp.json").write_text(_fastp_json())
    (tmp_path / "stats.txt").write_text(_SN_BLOCK + _COV_BLOCK.split("\n", 1)[1])
    _write_fai(tmp_path, [("chr1", 1000)])
    _write_vcf(tmp_path, [
        _rec(100, "A", "G", "1", info="ANN=" + _ann_entry(effect="missense_variant", impact="MODERATE", gene="geneA")),
        _rec(200, "C", "T", "1/1", info="ANN=" + _ann_entry(effect="synonymous_variant", impact="LOW", gene="geneA")),
        _rec(300, "G", "T", "0/1", info="ANN=" + _ann_entry(effect="stop_gained", impact="HIGH", gene="geneB")),
        _rec(400, "AT", "A", "1/1", info="ANN=" + _ann_entry(effect="frameshift_variant", impact="HIGH", gene="geneB")),
    ])
    (tmp_path / "lineage.tsv").write_text(_PATHO_TSV)
    (tmp_path / "dr.tsv").write_text("genome\tlineage:count\tmajor_lineage\nS1\trpoB:3\tRIF\n")

    proc = subprocess.run(
        [
            sys.executable, str(repo_root / "bin" / "stats_to_legacy.py"),
            "--sample", "S1",
            "--fastp-json", str(tmp_path / "fastp.json"),
            "--bam-stats", str(tmp_path / "stats.txt"),
            "--vcf", str(tmp_path / "calls.vcf"),
            "--ref-fai", str(tmp_path / "ref.fa.fai"),
            "--pathotypr-report", str(tmp_path / "lineage.tsv"),
            "--pathotypr-dr-report", str(tmp_path / "dr.tsv"),
        ],
        cwd=tmp_path, capture_output=True, text=True, timeout=120,
    )
    assert proc.returncode == 0, proc.stderr

    lines = (tmp_path / "S1.log").read_text().splitlines()
    keys = [ln.split("\t", 1)[0] for ln in lines]
    assert keys == stl.ORDERED_KEYS
    data = dict(ln.split("\t", 1) for ln in lines)

    # DAT_OUT is today's date, so it is excluded from the golden comparison.
    assert data.pop("DAT_OUT")
    profiles = {k: data.pop(k) for k in ("SNP_PROFILE", "HET_PROFILE", "INDEL_PROFILE")}

    golden = {
        "GID_IN": "S1",
        "FQ_IN": "S1_R1/R2",
        "ANC_OUT": "no",
        "GENOME_LEN": "1000",
        # fastp
        "PHE_OUT": "33",
        "FV_OUT": "FASTQ_SUCCESS",
        "TRD_OUT": "1000",
        "TRM_OUT": "900",
        "RDL_OUT": "100.0",
        "DUP_OUT": "10.00",
        "Q20_PCT": "98.00",
        "Q30_PCT": "95.00",
        "GC_PCT": "65.00",
        "LNG_OUT": "90.00",       # 90000 after-filter bases / 1000 bp genome
        # samtools stats
        "MAP_OUT": "900",
        "MAP_PCT": "90.00",
        "PP_OUT": "85.00",
        "UNMAP_PCT": "10.00",
        "SINGLE_PCT": "2.00",     # (900 mapped - 880 mapped-and-paired) / 1000
        "MQ0_PCT": "1.00",
        "ERR_RATE": "0.100",      # 0.001 emitted as a percent
        "ISIZE_MEAN": "300.0",
        "ISIZE_SD": "50.0",
        "STD_OUT": "35.5",
        "COV_OUT": "90.00",       # 90000 cigar bases / 1000 bp
        "COV_BREADTH": "60.00",
        "COV_MEDIAN": "1.00",
        "COV_BREADTH5": "20.00",
        "COV_BREADTH10": "10.00",
        "COV_EVENNESS": "1.469",
        # VCF
        "VAR_OUT": "4",
        "SNP_OUT": "3",
        "HOM_OUT": "3",
        "HET_OUT": "1",
        "HOMO_IND_OUT": "1",
        "TITV_OUT": "2.000",
        # snpEff ANN
        "ANN_PRESENT": "yes",
        "ANN_TOTAL": "4",
        "ANN_ANNOTATED": "4",
        "ANN_ANNOTATED_FRAC": "1.0000",
        "ANN_IMPACT_HIGH": "2",
        "ANN_IMPACT_MODERATE": "1",
        "ANN_IMPACT_LOW": "1",
        "ANN_IMPACT_MODIFIER": "0",
        "ANN_MISSENSE": "1",
        "ANN_SYNON": "1",
        "ANN_STOP_GAINED": "1",
        "ANN_STOP_LOST": "0",
        "ANN_START_LOST": "0",
        "ANN_FRAMESHIFT": "1",
        "ANN_INFRAME_INDEL": "0",
        "ANN_SPLICE": "0",
        "ANN_INTERGENIC": "0",
        "ANN_REGULATORY": "0",
        "ANN_LOF": "2",
        "ANN_CODING": "4",
        "ANN_MISSENSE_SILENT": "1.000",
        "ANN_WARN": "0",
        "ANN_DB_ERROR": "no",
        "ANN_GENE_BURDEN": "geneB:2:0:stop_gained;geneA:0:1:missense",
        # pathotypr
        "LIN_OUT": "L4.3",
        "LIN_COUNTS": "L4:10;L4.3:5",
        "DR_OUT": "RIF",
        "DR_COUNTS": "rpoB:3",
    }
    for key, expected in golden.items():
        assert data.pop(key) == expected, key
    # Every remaining legacy field is untouched.
    assert set(data.values()) == {"NA"}

    for name, profile in profiles.items():
        assert len(profile.split(",")) == stl.BIN_N, name
    # 200 bins over 1000 bp -> bin = POS // 5
    assert profiles["SNP_PROFILE"].split(",")[20] == "1"    # homozygous SNP at 100
    assert profiles["SNP_PROFILE"].split(",")[40] == "1"    # homozygous SNP at 200
    assert profiles["HET_PROFILE"].split(",")[60] == "1"    # het SNP at 300
    assert profiles["INDEL_PROFILE"].split(",")[80] == "1"  # homozygous indel at 400
