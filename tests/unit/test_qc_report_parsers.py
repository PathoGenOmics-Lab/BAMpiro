"""Unit tests for the file parsers of bin/qcreport/parsers.py.

Each parser gets a happy path built from a small inline fixture in tmp_path, plus
the empty / missing-file behaviour, which differs per function: most return [], {}
or None, but parse_summary and consensus_stats deliberately raise (the pipeline
always hands them a file it just produced). The GFF and samplesheet parsers are
additionally exercised against the committed fixture cohort in tests/data/, so a
change to those files breaks a test instead of silently changing a report.
"""

from __future__ import annotations

import gzip

import pytest
from qcreport import metrics as qc_metrics       # build_gene_map: a metric helper over a parsed GFF
from qcreport import parsers as qc_parsers


def write(path, text):
    """Write a fixture file and return its path as a string (the parsers take str paths)."""
    path.write_text(text, encoding="utf-8")
    return str(path)


def write_gz(path, text):
    with gzip.open(path, "wt", encoding="utf-8") as fh:
        fh.write(text)
    return str(path)


MISSING = "/nonexistent/directory/definitely-not-here.tsv"


# --------------------------------------------------------------------------- parse_summary


SUMMARY = (
    "sample_id\tmean_depth\tbreadth_pct\tsnps\tlineage\n"
    "S1\t42.5\t99.1\t1200\tL4\n"
    "S2\t8.0\t70.0\t900\tL2\n"
)


def test_parse_summary_reads_one_dict_per_sample(tmp_path):
    rows = qc_parsers.parse_summary(write(tmp_path / "summary.tsv", SUMMARY))
    assert sorted(rows) == ["S1", "S2"]
    assert rows["S1"] == {"sample_id": "S1", "mean_depth": "42.5", "breadth_pct": "99.1",
                          "snps": "1200", "lineage": "L4"}


def test_parse_summary_pads_a_short_row_with_empty_strings(tmp_path):
    path = write(tmp_path / "summary.tsv", SUMMARY + "S3\t10.0\n")
    rows = qc_parsers.parse_summary(path)
    assert rows["S3"] == {"sample_id": "S3", "mean_depth": "10.0", "breadth_pct": "",
                          "snps": "", "lineage": ""}


def test_parse_summary_skips_blank_lines_and_rows_without_a_sample_id(tmp_path):
    path = write(tmp_path / "summary.tsv", SUMMARY + "\n   \n\t5.0\t5.0\t5\tL1\n")
    assert sorted(qc_parsers.parse_summary(path)) == ["S1", "S2"]


def test_parse_summary_falls_back_to_the_first_column_without_a_sample_id_header(tmp_path):
    path = write(tmp_path / "summary.tsv", "name\tmean_depth\nS9\t30\n")
    rows = qc_parsers.parse_summary(path)
    assert list(rows) == ["S9"]
    assert rows["S9"]["mean_depth"] == "30"


def test_parse_summary_keeps_the_last_row_of_a_duplicated_sample(tmp_path):
    path = write(tmp_path / "summary.tsv", SUMMARY + "S1\t99.9\t50.0\t1\tL9\n")
    assert qc_parsers.parse_summary(path)["S1"]["mean_depth"] == "99.9"


def test_parse_summary_of_an_empty_file_is_an_empty_cohort(tmp_path):
    assert qc_parsers.parse_summary(write(tmp_path / "summary.tsv", "")) == {}


def test_parse_summary_of_a_missing_file_raises(tmp_path):
    # the summary is a required input the pipeline just produced: absence is a bug, not a state
    with pytest.raises(OSError):
        qc_parsers.parse_summary(str(tmp_path / "nope.tsv"))


# --------------------------------------------------------------------------- consensus_stats


# 80 callable bases, 10 N, 5 gap, 5 IUPAC -> callable 80%, missing 15%, IUPAC 5%
CONSENSUS_SEQ = "ACGT" * 20 + "N" * 10 + "-" * 5 + "R" * 5
CONSENSUS_FASTA = ">S1 some description\n" + "\n".join(
    CONSENSUS_SEQ[i:i + 60] for i in range(0, len(CONSENSUS_SEQ), 60)) + "\n"


def assert_consensus_fixture(stats):
    assert stats["consensus_id"] == "S1"
    assert stats["length"] == 100
    assert stats["callable_pct"] == pytest.approx(80.0)
    assert stats["missing_pct"] == pytest.approx(15.0)
    assert stats["iupac_pct"] == pytest.approx(5.0)
    assert stats["longest_n_run"] == 15   # the N run and the gap run are contiguous
    assert stats["n_gaps"] == 1
    assert len(stats["miss_profile"]) == qc_parsers.NBINS


def test_consensus_stats_measures_a_plain_fasta(tmp_path):
    assert_consensus_fixture(qc_parsers.consensus_stats(write(tmp_path / "S1.fa", CONSENSUS_FASTA)))


def test_consensus_stats_reads_a_gzipped_fasta_identically(tmp_path):
    assert_consensus_fixture(qc_parsers.consensus_stats(write_gz(tmp_path / "S1.fa.gz", CONSENSUS_FASTA)))


def test_consensus_stats_is_case_insensitive(tmp_path):
    lower = qc_parsers.consensus_stats(write(tmp_path / "S1.fa", CONSENSUS_FASTA.lower()))
    assert lower["callable_pct"] == pytest.approx(80.0)
    assert lower["missing_pct"] == pytest.approx(15.0)


def test_consensus_stats_measures_every_record_but_keeps_the_first_id(tmp_path):
    path = write(tmp_path / "multi.fa", ">contig1\nACGT\n>contig2\nNNNN\n")
    stats = qc_parsers.consensus_stats(path)
    assert stats["consensus_id"] == "contig1"
    assert stats["length"] == 8
    assert stats["callable_pct"] == pytest.approx(50.0)
    assert stats["missing_pct"] == pytest.approx(50.0)


def test_consensus_stats_counts_separate_gap_runs(tmp_path):
    path = write(tmp_path / "S1.fa", ">S1\nACNNGT--AC\n")
    stats = qc_parsers.consensus_stats(path)
    assert stats["n_gaps"] == 2
    assert stats["longest_n_run"] == 2


def test_consensus_stats_of_a_zero_length_sequence_does_not_divide_by_zero(tmp_path):
    # a header with no bases: the `len(seq) or 1` guard keeps every percentage at 0
    stats = qc_parsers.consensus_stats(write(tmp_path / "empty.fa", ">S1 empty\n"))
    assert stats["length"] == 0
    assert stats["missing_pct"] == 0.0
    assert stats["iupac_pct"] == 0.0
    assert stats["callable_pct"] == 0.0
    assert stats["longest_n_run"] == 0
    assert stats["n_gaps"] == 0
    assert stats["miss_profile"] == [0.0] * qc_parsers.NBINS


def test_consensus_stats_of_a_completely_empty_file_has_no_id(tmp_path):
    stats = qc_parsers.consensus_stats(write(tmp_path / "empty.fa", ""))
    assert stats["consensus_id"] is None
    assert stats["length"] == 0


def test_consensus_stats_of_a_missing_file_raises(tmp_path):
    with pytest.raises(OSError):
        qc_parsers.consensus_stats(str(tmp_path / "nope.fa"))


# --------------------------------------------------------------------------- parse_lineage_colors


def test_parse_lineage_colors_reads_a_palette_tsv(tmp_path):
    path = write(tmp_path / "colors.tsv",
                 "# lineage palette\n"
                 "L1\t#ff0000\tIndo-Oceanic\n"
                 "L2\t#0f0\tEast-Asian\n"
                 "\n"
                 "L3\t#00ff00\n")
    assert qc_parsers.parse_lineage_colors(path) == {"L1": "#ff0000", "L2": "#0f0", "L3": "#00ff00"}


def test_parse_lineage_colors_rejects_anything_that_is_not_a_hex_colour(tmp_path):
    path = write(tmp_path / "colors.tsv",
                 "L1\tred\n"            # not a hex colour
                 "L2\t#ff00\n"          # wrong length
                 "L3\t#abcdef\n"
                 "L4\n")                # no colour column
    assert qc_parsers.parse_lineage_colors(path) == {"L3": "#abcdef"}


@pytest.mark.parametrize("path", [None, "", MISSING])
def test_parse_lineage_colors_of_an_absent_file_is_empty(path):
    assert qc_parsers.parse_lineage_colors(path) == {}


# --------------------------------------------------------------------------- parse_bed


BED = ("# a masking BED\n"
       "track name=mask\n"
       "browser hide all\n"
       "chr\t0\t100\tPE_PGRS\n"
       "chr\t500\t900\n"
       "\n"
       "chr\tbad\t900\n"
       "chr\t950\n")


def test_parse_bed_reads_the_half_open_intervals(tmp_path):
    assert qc_parsers.parse_bed(write(tmp_path / "mask.bed", BED)) == [(0, 100), (500, 900)]


def test_parse_bed_reads_a_gzipped_bed(tmp_path):
    assert qc_parsers.parse_bed(write_gz(tmp_path / "mask.bed.gz", BED)) == [(0, 100), (500, 900)]


@pytest.mark.parametrize("path", [None, "", MISSING])
def test_parse_bed_of_an_absent_file_is_empty(path):
    assert qc_parsers.parse_bed(path) == []


def test_parse_bed_of_an_empty_file_is_empty(tmp_path):
    assert qc_parsers.parse_bed(write(tmp_path / "mask.bed", "")) == []


# --------------------------------------------------------------------------- parse_gene_burden


GENE_BURDEN = (
    "gene\thigh\tmoderate\tlow\tmodifier\tdominant_effect\tn_samples\ttotal_impactful\tstart\tend\n"
    "katG\t3\t1\t0\t5\tstop_gained\t4\t4\t3301\t4700\n"
    "\n"
    "rpoB\t1\t5\t2\t1\tmissense_variant\t3\t6\t1601\t3100\n"
)


def test_parse_gene_burden_sorts_by_descending_impactful_count(tmp_path):
    rows = qc_parsers.parse_gene_burden(write(tmp_path / "burden.tsv", GENE_BURDEN))
    assert [r["gene"] for r in rows] == ["rpoB", "katG"]
    assert rows[0] == {"gene": "rpoB", "high": 1, "moderate": 5, "low": 2, "modifier": 1,
                       "dominant_effect": "missense_variant", "n_samples": 3,
                       "total_impactful": 6, "start": 1601, "end": 3100}


def test_parse_gene_burden_derives_a_missing_total_from_high_plus_moderate(tmp_path):
    path = write(tmp_path / "burden.tsv",
                 "gene\thigh\tmoderate\n"
                 "katG\t3\t4\n"
                 "rpoB\t\t\n")
    rows = qc_parsers.parse_gene_burden(path)
    assert rows[0]["total_impactful"] == 7
    assert rows[1] == {"gene": "rpoB", "high": 0, "moderate": 0, "low": None, "modifier": None,
                       "dominant_effect": "NA", "n_samples": None, "total_impactful": 0,
                       "start": None, "end": None}


def test_parse_gene_burden_skips_rows_without_a_gene(tmp_path):
    path = write(tmp_path / "burden.tsv", GENE_BURDEN + "\t9\t9\t0\t0\tx\t1\t18\t1\t2\n")
    assert [r["gene"] for r in qc_parsers.parse_gene_burden(path)] == ["rpoB", "katG"]


def test_parse_gene_burden_is_capped_at_two_hundred_genes(tmp_path):
    body = "".join(f"gene{i:04d}\t1\t1\t0\t0\teff\t1\t{i}\t1\t2\n" for i in range(250))
    path = write(tmp_path / "burden.tsv", GENE_BURDEN.splitlines()[0] + "\n" + body)
    rows = qc_parsers.parse_gene_burden(path)
    assert len(rows) == 200
    assert rows[0]["gene"] == "gene0249"     # the highest total_impactful survives the cap


@pytest.mark.parametrize("path", [None, "", MISSING])
def test_parse_gene_burden_of_an_absent_file_is_empty(path):
    assert qc_parsers.parse_gene_burden(path) == []


def test_parse_gene_burden_of_an_empty_file_is_empty(tmp_path):
    assert qc_parsers.parse_gene_burden(write(tmp_path / "burden.tsv", "")) == []


# --------------------------------------------------------------------------- parse_pnps


PNPS = (
    "gene\tgene_name\tn_samples\tn_pairs\tmean_dN\tmean_dS\tdNdS\teffect_dominant\n"
    "TEST_0001\tdnaA\t4\t6\t0.0012\t0.0009\t1.33\tmissense\n"
    "\n"
    "TEST_0002\trpoB\t4\t6\t0.0000\t0.0000\tNA\tsynonymous\n"
    "TEST_0003\tkatG\t4\t6\t0.0020\t0.0000\tinf\tmissense\n"
)


def test_parse_pnps_reads_the_eskaks_aggregate_columns(tmp_path):
    rows = qc_parsers.parse_pnps(write(tmp_path / "pnps.tsv", PNPS))
    assert [r["gene"] for r in rows] == ["dnaA", "rpoB", "katG"]
    assert rows[0] == {"gene": "dnaA", "pn": pytest.approx(0.0012), "ps": pytest.approx(0.0009),
                       "pnps": pytest.approx(1.33), "n": 6.0, "effect": "missense"}


@pytest.mark.parametrize("row_index", [1, 2])
def test_parse_pnps_reads_the_non_finite_sentinels_as_none(tmp_path, row_index):
    # NA and inf must never be blindly float()-ed: an inf dN/dS would poison every downstream axis
    rows = qc_parsers.parse_pnps(write(tmp_path / "pnps.tsv", PNPS))
    assert rows[row_index]["pnps"] is None


def test_parse_pnps_falls_back_to_the_gene_column(tmp_path):
    path = write(tmp_path / "pnps.tsv", "gene\tdNdS\nTEST_0001\t0.5\n")
    assert qc_parsers.parse_pnps(path) == [{"gene": "TEST_0001", "pn": None, "ps": None,
                                    "pnps": 0.5, "n": None, "effect": None}]


def test_parse_pnps_skips_rows_without_a_gene(tmp_path):
    path = write(tmp_path / "pnps.tsv", "gene_name\tdNdS\n\t0.5\ndnaA\t0.7\n")
    assert [r["gene"] for r in qc_parsers.parse_pnps(path)] == ["dnaA"]


def test_parse_pnps_is_capped_at_three_hundred_genes(tmp_path):
    body = "".join(f"g{i:04d}\t0.5\n" for i in range(320))
    assert len(qc_parsers.parse_pnps(write(tmp_path / "pnps.tsv", "gene\tdNdS\n" + body))) == 300


@pytest.mark.parametrize("path", [None, "", MISSING])
def test_parse_pnps_of_an_absent_file_is_empty(path):
    assert qc_parsers.parse_pnps(path) == []


def test_parse_pnps_of_an_empty_file_is_empty(tmp_path):
    assert qc_parsers.parse_pnps(write(tmp_path / "pnps.tsv", "")) == []


# --------------------------------------------------------------------------- parse_dr


DR = (
    "# pathotypr drug-resistance calls\n"
    "sample\tdrug\tgene\tmutation\tgrade\tmarker_name\taf\tdp\n"
    "S1\tINH\tkatG\tS315T\t1) Assoc w R\tkatG_S315T\t0.98\t40\n"
    "S1\tRIF\trpoB\tS450L\t2) Assoc w R - Interim\trpoB_S450L\t0.55\t22\n"
    "S2\tEMB\tembB\tM306V\tUncertain\t\tNA\t.\n"
)


def test_parse_dr_reads_the_calls_samples_and_drug_order(tmp_path):
    dr = qc_parsers.parse_dr(write(tmp_path / "dr.tsv", DR))
    assert dr["samples"] == ["S1", "S2"]
    # drugs follow the curated first/second-line order, not alphabetical order
    assert dr["drugs"] == ["RIF", "INH", "EMB"]
    assert len(dr["calls"]) == 3
    assert dr["calls"][0] == {"s": "S1", "drug": "INH", "gene": "katG", "mutation": "S315T",
                              "grade": "1) Assoc w R", "gn": 1, "marker": "katG_S315T",
                              "af": pytest.approx(0.98), "dp": 40}


def test_parse_dr_extracts_the_leading_who_grade_number(tmp_path):
    calls = qc_parsers.parse_dr(write(tmp_path / "dr.tsv", DR))["calls"]
    assert [c["gn"] for c in calls] == [1, 2, None]


def test_parse_dr_reads_the_missing_sentinels_as_none(tmp_path):
    last = qc_parsers.parse_dr(write(tmp_path / "dr.tsv", DR))["calls"][-1]
    assert last["af"] is None
    assert last["dp"] is None
    assert last["marker"] == ""


def test_parse_dr_falls_back_to_the_marker_column(tmp_path):
    path = write(tmp_path / "dr.tsv", "sample\tdrug\tmarker\nS1\tRIF\tm1\n")
    assert qc_parsers.parse_dr(path)["calls"][0]["marker"] == "m1"


def test_parse_dr_puts_an_unknown_drug_last(tmp_path):
    path = write(tmp_path / "dr.tsv", "sample\tdrug\nS1\tZZZ\nS1\tRIF\n")
    assert qc_parsers.parse_dr(path)["drugs"] == ["RIF", "ZZZ"]


def test_parse_dr_skips_blank_lines_and_rows_without_a_sample(tmp_path):
    path = write(tmp_path / "dr.tsv", "sample\tdrug\n\nS1\tRIF\n\t INH\n")
    dr = qc_parsers.parse_dr(path)
    assert dr["samples"] == ["S1"]
    assert len(dr["calls"]) == 1


def test_parse_dr_reads_a_float_formatted_depth(tmp_path):
    """A depth written as '40.0' is still a depth of 40, not a parse failure."""
    path = write(tmp_path / "dr.tsv", "sample\tdrug\tdp\nS1\tRIF\t40.0\n")
    assert qc_parsers.parse_dr(path)["calls"][0]["dp"] == 40


def test_parse_dr_of_a_non_numeric_depth_is_none(tmp_path):
    path = write(tmp_path / "dr.tsv", "sample\tdrug\tdp\nS1\tRIF\tdeep\n")
    assert qc_parsers.parse_dr(path)["calls"][0]["dp"] is None


@pytest.mark.parametrize("path", [None, "", MISSING])
def test_parse_dr_of_an_absent_file_is_none(path):
    assert qc_parsers.parse_dr(path) is None


def test_parse_dr_without_any_call_is_none(tmp_path):
    assert qc_parsers.parse_dr(write(tmp_path / "dr.tsv", "")) is None
    assert qc_parsers.parse_dr(write(tmp_path / "dr2.tsv", "sample\tdrug\n")) is None


# --------------------------------------------------------------------------- parse_gene_conversion


GCONV_HEADER = (
    "sample\tpair_id\tcontig\tdonor\tverdict\treason\tstart\tend\tspan_bp\t"
    "post_conv\tlog10_bf\tlog10_bf_vs_null\tmismap_frac\tstart_ci\tend_ci\t"
    "n_sites\tn_sites_outside\tn_undetermined\tdonor_af_in\tdonor_af_outside\t"
    "min_depth\tcis_reads\tbreakpoint_reads\tdonor_only_reads\n"
)
GCONV = GCONV_HEADER + (
    "S1\t7\tPPE34\tPPE12\tgene_conversion\tlog10 Bayes factor 18.4 over the best alternative"
    "\t1000\t1400\t401\t1.0\t18.4\t612.5\t0.0031\t1000-1000\t1400-1400\t6\t9\t0\t0.97\t0.01\t18\t11\t3\t0\n"
    "S1\t9\tPE_PGRS4\tPE_PGRS5\tmismapping\tthe locus is explained by reads from the donor at 0.48"
    "\t2000\t2300\t301\t0.02\t0.4\t2.1\t0.4812\t2000-2140\t2210-2300\t4\t12\t1\t0.55\t0.48\t9\t5\t0\t7\n"
    # the whole-locus case: no site outside the tract, so donor_af_outside cannot be measured at all
    "S2\t7\tPPE34\tPPE12\tambiguous\tthe tract covers every diagnostic site\t1000\t1600\t601"
    "\t0.31\t1.1\t1.1\t0.9701\t1000-1000\t1600-1600\t15\t0\t0\t0.91\t\t12\t8\t0\t8\n"
)


def test_parse_gene_conversion_reads_the_tracts_samples_and_verdict_counts(tmp_path):
    g = qc_parsers.parse_gene_conversion(write(tmp_path / "gconv.tsv", GCONV))
    assert g["samples"] == ["S1", "S2"]
    assert g["counts"] == {"gene_conversion": 1, "mismapping": 1, "ambiguous": 1,
                           "coverage_shift": 0}
    assert g["tracts"][0] == {
        "s": "S1", "pair": "7", "contig": "PPE34", "donor": "PPE12", "verdict": "gene_conversion",
        "reason": "log10 Bayes factor 18.4 over the best alternative",
        "start": 1000, "end": 1400, "span": 401,
        "post": pytest.approx(1.0), "bf": pytest.approx(18.4), "bf_null": pytest.approx(612.5),
        "mismap": pytest.approx(0.0031), "start_ci": "1000-1000", "end_ci": "1400-1400",
        "n_sites": 6, "n_out": 9, "n_undet": 0,
        "af_in": pytest.approx(0.97), "af_out": pytest.approx(0.01),
        "depth": 18, "cis_reads": 11, "bp_reads": 3, "donor_only": 0}


def test_parse_gene_conversion_keeps_the_model_numbers_and_the_checkable_ones(tmp_path):
    """The Bayes factor is the conclusion and the allele fractions are what it can be checked
    against. Both travel to the panel: a verdict nobody can audit is not much use."""
    tracts = qc_parsers.parse_gene_conversion(write(tmp_path / "gconv.tsv", GCONV))["tracts"]

    assert [t["bf"] for t in tracts] == [pytest.approx(18.4), pytest.approx(0.4),
                                         pytest.approx(1.1)]
    assert [t["mismap"] for t in tracts] == [pytest.approx(0.0031), pytest.approx(0.4812),
                                             pytest.approx(0.9701)]
    assert [t["af_in"] for t in tracts] == [pytest.approx(0.97), pytest.approx(0.55),
                                            pytest.approx(0.91)]


def test_parse_gene_conversion_keeps_the_three_pieces_of_deciding_evidence(tmp_path):
    """donor_af_in, donor_af_outside and breakpoint_reads are what separate a conversion from a
    mismapping, so the panel reads them per row; none of them may be collapsed or dropped."""
    tracts = qc_parsers.parse_gene_conversion(write(tmp_path / "gconv.tsv", GCONV))["tracts"]
    assert [t["af_in"] for t in tracts] == [pytest.approx(0.97), pytest.approx(0.55), pytest.approx(0.91)]
    assert [t["af_out"] for t in tracts] == [pytest.approx(0.01), pytest.approx(0.48), None]
    assert [t["bp_reads"] for t in tracts] == [3, 0, 0]


def test_parse_gene_conversion_reads_an_unmeasurable_value_as_none_not_zero(tmp_path):
    """A tract covering every diagnostic site has no donor_af_outside. Reading the empty field as
    0.0 would make the least-bounded case look like the best-bounded one."""
    last = qc_parsers.parse_gene_conversion(write(tmp_path / "gconv.tsv", GCONV))["tracts"][-1]
    assert last["af_out"] is None
    assert last["n_out"] == 0


def _gconv_row(**cells):
    """One row laid out against GCONV_HEADER, so a column added upstream cannot silently shift
    the fields of these fixtures into each other."""
    columns = GCONV_HEADER.rstrip("\n").split("\t")
    defaults = dict(sample="S1", pair_id="1", contig="c", donor="d", verdict="ambiguous",
                    reason="r", start="1", end="2", span_bp="2", post_conv="0.5", log10_bf="1.0",
                    log10_bf_vs_null="2.0", mismap_frac="0.01", start_ci="1-1", end_ci="2-2",
                    n_sites="3", n_sites_outside="3", n_undetermined="0", donor_af_in="1.0",
                    donor_af_outside="0.0", min_depth="9", cis_reads="3", breakpoint_reads="1",
                    donor_only_reads="0")
    defaults.update(cells)
    return "\t".join(defaults[c] for c in columns) + "\n"


def test_parse_gene_conversion_normalises_the_verdict_case(tmp_path):
    row = _gconv_row(verdict="Gene_Conversion")
    g = qc_parsers.parse_gene_conversion(write(tmp_path / "gconv.tsv", GCONV_HEADER + row))
    assert g["tracts"][0]["verdict"] == "gene_conversion"
    assert g["counts"]["gene_conversion"] == 1


def test_parse_gene_conversion_counts_an_unrecognised_verdict_without_losing_the_known_ones(tmp_path):
    row = _gconv_row(verdict="something_else")
    counts = qc_parsers.parse_gene_conversion(write(tmp_path / "gconv.tsv", GCONV_HEADER + row))["counts"]
    assert counts == {"gene_conversion": 0, "mismapping": 0, "ambiguous": 0,
                      "coverage_shift": 0, "something_else": 1}


def test_parse_gene_conversion_reads_a_float_formatted_count(tmp_path):
    row = _gconv_row(start="1.0", end="2.0", span_bp="2.0", n_sites="3.0", n_sites_outside="3.0",
                     min_depth="9.0", cis_reads="3.0", breakpoint_reads="1.0",
                     donor_only_reads="0.0")
    t = qc_parsers.parse_gene_conversion(write(tmp_path / "gconv.tsv", GCONV_HEADER + row))["tracts"][0]
    assert (t["start"], t["n_sites"], t["bp_reads"], t["depth"]) == (1, 3, 1, 9)


def test_parse_gene_conversion_reads_the_missing_sentinels_as_none(tmp_path):
    row = _gconv_row(start="NA", end=".", span_bp="", log10_bf="NA", mismap_frac=".",
                     donor_af_in="NA", donor_af_outside=".", min_depth="", cis_reads="NA",
                     breakpoint_reads=".", donor_only_reads="")
    t = qc_parsers.parse_gene_conversion(write(tmp_path / "gconv.tsv", GCONV_HEADER + row))["tracts"][0]
    assert [t[k] for k in ("start", "end", "span", "bf", "mismap", "af_in", "af_out", "depth",
                           "cis_reads", "bp_reads", "donor_only")] == [None] * 11


def test_parse_gene_conversion_skips_blank_lines_and_rows_without_a_sample(tmp_path):
    body = "\n" + _gconv_row(sample="")
    g = qc_parsers.parse_gene_conversion(write(tmp_path / "gconv.tsv", GCONV + body))
    assert g["samples"] == ["S1", "S2"]
    assert len(g["tracts"]) == 3


@pytest.mark.parametrize("path", [None, "", MISSING])
def test_parse_gene_conversion_of_an_absent_file_is_none(path):
    assert qc_parsers.parse_gene_conversion(path) is None


def test_parse_gene_conversion_without_any_tract_is_none(tmp_path):
    """The pipeline always writes the file, header-only when the stage found nothing, so an empty
    result must hide the panel rather than render an empty one."""
    assert qc_parsers.parse_gene_conversion(write(tmp_path / "a.tsv", "")) is None
    assert qc_parsers.parse_gene_conversion(write(tmp_path / "b.tsv", GCONV_HEADER)) is None


# --------------------------------------------------------------------------- parse_kraken


KRAKEN = (
    " 12.00\t1200\t1200\tU\t0\tunclassified\n"
    " 88.00\t8800\t100\tR\t1\troot\n"
    " 70.00\t7000\t7000\tS\t1773\t    Mycobacterium tuberculosis\n"
    " 10.00\t1000\t1000\tS\t1280\t    Staphylococcus aureus\n"
    "short line\n"
    "not-a-number\t1\t1\tS\t9\t  Bogus\n"
)


def test_parse_kraken_summarises_the_species_composition(tmp_path):
    path = write(tmp_path / "S1__RUN1.kraken.report", KRAKEN)
    sample = qc_parsers.parse_kraken([path])["samples"][0]
    assert sample["s"] == "S1"
    assert sample["unclassified"] == 12.0
    assert sample["classified"] == 88.0
    assert sample["primary"] == {"name": "Mycobacterium tuberculosis", "pct": 70.0}
    assert sample["secondary"] == {"name": "Staphylococcus aureus", "pct": 10.0}
    assert len(sample["top"]) == 2


def test_parse_kraken_derives_the_sample_id_from_the_file_name(tmp_path):
    a = write(tmp_path / "S1__RUN1.kraken.report", KRAKEN)
    b = write(tmp_path / "S2.kraken.report", KRAKEN)
    assert [s["s"] for s in qc_parsers.parse_kraken([a, b])["samples"]] == ["S1", "S2"]


def test_parse_kraken_keeps_only_the_top_six_species(tmp_path):
    body = "".join(f"{i}.00\t10\t10\tS\t{i}\ttaxon{i}\n" for i in range(1, 12))
    path = write(tmp_path / "S1.kraken.report", body)
    sample = qc_parsers.parse_kraken([path])["samples"][0]
    assert len(sample["top"]) == 6
    assert sample["primary"]["name"] == "taxon11"      # sorted by descending percent


def test_parse_kraken_keeps_a_fully_unclassified_sample(tmp_path):
    path = write(tmp_path / "S1.kraken.report", "100.00\t10\t10\tU\t0\tunclassified\n")
    sample = qc_parsers.parse_kraken([path])["samples"][0]
    assert sample["unclassified"] == 100.0
    assert sample["primary"] is None
    assert sample["secondary"] is None


def test_parse_kraken_skips_an_uninformative_report(tmp_path):
    path = write(tmp_path / "S1.kraken.report", "88.00\t8800\t100\tR\t1\troot\n")
    assert qc_parsers.parse_kraken([path]) is None


@pytest.mark.parametrize("paths", [None, [], [MISSING], ["", None]])
def test_parse_kraken_without_any_readable_report_is_none(paths):
    assert qc_parsers.parse_kraken(paths) is None


# --------------------------------------------------------------------------- parse_gff / parse_gene_locus


def test_parse_gff_reads_the_committed_reference_annotation(data_dir):
    genes = qc_parsers.parse_gff(str(data_dir / "reference.gff3"))
    assert genes == [{"name": "dnaA", "start": 201, "end": 1400},
                     {"name": "rpoB", "start": 1601, "end": 3100},
                     {"name": "katG", "start": 3301, "end": 4700}]


def test_parse_gff_reads_a_gzipped_gff(tmp_path, data_dir):
    text = (data_dir / "reference.gff3").read_text(encoding="utf-8")
    assert qc_parsers.parse_gff(write_gz(tmp_path / "ref.gff3.gz", text)) == \
        qc_parsers.parse_gff(str(data_dir / "reference.gff3"))


def test_parse_gff_prefers_a_gene_feature_over_a_cds_of_the_same_name(tmp_path):
    path = write(tmp_path / "ref.gff3",
                 "chr\tt\tCDS\t50\t60\t.\t+\t0\tName=abc\n"
                 "chr\tt\tgene\t10\t100\t.\t+\t.\tName=abc\n")
    assert qc_parsers.parse_gff(path) == [{"name": "abc", "start": 10, "end": 100}]


def test_parse_gff_falls_back_through_the_attribute_keys_then_to_the_coordinates(tmp_path):
    path = write(tmp_path / "ref.gff3",
                 "chr\tt\tgene\t10\t20\t.\t+\t.\tlocus_tag=LT1\n"
                 "chr\tt\tgene\t30\t40\t.\t+\t.\tID=only-an-id\n"
                 "chr\tt\tgene\t50\t60\t.\t+\t.\t.\n")
    assert [g["name"] for g in qc_parsers.parse_gff(path)] == ["LT1", "only-an-id", "50-60"]


def test_parse_gff_skips_non_gene_features_headers_and_bad_coordinates(tmp_path):
    path = write(tmp_path / "ref.gff3",
                 "##gff-version 3\n"
                 "chr\tt\tregion\t1\t5000\t.\t+\t.\tName=chr\n"
                 "chr\tt\tgene\tx\t20\t.\t+\t.\tName=bad\n"
                 "a line with no tabs\n"
                 "chr\tt\tgene\t10\t20\t.\t+\t.\tName=ok\n")
    assert [g["name"] for g in qc_parsers.parse_gff(path)] == ["ok"]


@pytest.mark.parametrize("path", [None, "", MISSING])
def test_parse_gff_of_an_absent_file_is_empty(path):
    assert qc_parsers.parse_gff(path) == []


def test_parse_gff_attribute_lookup_matches_the_whole_key_not_a_suffix(tmp_path):
    """The attribute key is matched whole, so 'locus_tag' is not read out of 'old_locus_tag'
    and 'gene' is not read out of 'pseudogene'. Both attributes are legal GFF3 and routine in
    RefSeq / Prokka output, so a substring search named features after an obsolete tag or after
    a pseudogene TYPE ('unitary')."""
    path = write(tmp_path / "ref.gff3",
                 "chr\tt\tgene\t1\t100\t.\t+\t.\tID=g1;old_locus_tag=OLDTAG;locus_tag=NEWTAG\n"
                 "chr\tt\tgene\t200\t300\t.\t+\t.\tID=g2;pseudogene=unitary;locus_tag=LT2\n")
    assert [g["name"] for g in qc_parsers.parse_gff(path)] == ["NEWTAG", "LT2"]


def test_parse_gff_tolerates_spaces_around_the_attribute_separator(tmp_path):
    path = write(tmp_path / "ref.gff3", "chr\tt\tgene\t1\t100\t.\t+\t.\tID=g1; locus_tag=LT1 ; note=x\n")
    assert [g["name"] for g in qc_parsers.parse_gff(path)] == ["LT1"]


def test_parse_gff_keeps_an_equals_sign_inside_an_attribute_value(tmp_path):
    """Only the FIRST '=' splits key from value, so a value containing '=' survives intact."""
    path = write(tmp_path / "ref.gff3", "chr\tt\tgene\t1\t100\t.\t+\t.\tName=a=b;locus_tag=LT1\n")
    assert [g["name"] for g in qc_parsers.parse_gff(path)] == ["a=b"]


def test_parse_gene_locus_maps_the_committed_gene_names_to_locus_tags(data_dir):
    assert qc_parsers.parse_gene_locus(str(data_dir / "reference.gff3")) == {
        "dnaA": "TEST_0001", "rpoB": "TEST_0002", "katG": "TEST_0003"}


def test_parse_gene_locus_keeps_the_first_locus_tag_seen_for_a_gene(tmp_path):
    path = write(tmp_path / "ref.gff3",
                 "chr\tt\tgene\t10\t20\t.\t+\t.\tgene=rpoB;locus_tag=FIRST\n"
                 "chr\tt\tCDS\t10\t20\t.\t+\t0\tgene=rpoB;locus_tag=SECOND\n")
    assert qc_parsers.parse_gene_locus(path) == {"rpoB": "FIRST"}


def test_parse_gene_locus_needs_both_a_name_and_a_locus_tag(tmp_path):
    path = write(tmp_path / "ref.gff3",
                 "chr\tt\tgene\t10\t20\t.\t+\t.\tlocus_tag=NO_NAME\n"
                 "chr\tt\tgene\t30\t40\t.\t+\t.\tName=no_locus\n"
                 "chr\tt\tgene\t50\t60\t.\t+\t.\tName=both;locus_tag=LT9\n")
    assert qc_parsers.parse_gene_locus(path) == {"both": "LT9"}


@pytest.mark.parametrize("path", [None, "", MISSING])
def test_parse_gene_locus_of_an_absent_file_is_empty(path):
    assert qc_parsers.parse_gene_locus(path) == {}


def test_parse_gene_locus_attribute_lookup_matches_the_whole_key_not_a_suffix(tmp_path):
    """Same whole-key match as parse_gff: 'old_locus_tag' must not answer for 'locus_tag'
    (that mapped the gene to its OBSOLETE tag) and 'pseudogene' must not answer for 'gene'
    (that turned the pseudogene TYPE into a gene name)."""
    path = write(tmp_path / "ref.gff3",
                 "chr\tt\tgene\t1\t100\t.\t+\t.\tID=g1;old_locus_tag=MTV_OLD;locus_tag=Rv0667;gene=rpoB\n"
                 "chr\tt\tgene\t200\t300\t.\t+\t.\tID=g2;pseudogene=unitary;locus_tag=LT2\n")
    assert qc_parsers.parse_gene_locus(path) == {"rpoB": "Rv0667"}


def test_build_gene_map_takes_the_real_locus_tag_not_the_obsolete_one(tmp_path):
    """build_gene_map overlays H37Rv-style tags from the GFF onto the curated map. An Rv-style
    'old_locus_tag' read in place of 'locus_tag' passed the _RV_LOCUS_RE filter and silently
    overrode the curated locus with a retired one."""
    path = write(tmp_path / "ref.gff3",
                 "chr\tt\tgene\t1\t100\t.\t+\t.\tgene=rpoB;old_locus_tag=Rv9999;locus_tag=Rv0667\n")
    assert qc_metrics.build_gene_map(path)["rpoB"] == "Rv0667"


# --------------------------------------------------------------------------- samplesheet parsers


ANNOTATED_SHEET = (
    "sampleId\trunId\tr1\tr2\tpassage\tpatient\ttreatment\tdose\n"
    "S1\tRUN1\ta_R1.fq.gz\ta_R2.fq.gz\tP1\tPAT1\tdrugA\t10\n"
    "S2\tRUN1\tb_R1.fq.gz\tb_R2.fq.gz\tP2\tPAT1\tdrugA\t20\n"
    "S3\tRUN1\tc_R1.fq.gz\tc_R2.fq.gz\tP1\tPAT2\tcontrol\t\n"
)


def test_parse_sample_meta_keeps_only_the_user_annotation_columns(tmp_path):
    meta = qc_parsers.parse_sample_meta(write(tmp_path / "sheet.tsv", ANNOTATED_SHEET))
    # the sample id, the pipeline's file/reference columns and the numeric dose are all excluded
    assert meta["fields"] == ["passage", "patient", "treatment"]
    assert meta["rows"]["S1"] == {"passage": "P1", "patient": "PAT1", "treatment": "drugA"}
    assert sorted(meta["rows"]) == ["S1", "S2", "S3"]


def test_parse_sample_meta_tags_the_time_group_and_treatment_columns(tmp_path):
    meta = qc_parsers.parse_sample_meta(write(tmp_path / "sheet.tsv", ANNOTATED_SHEET))
    assert meta["time_field"] == "passage"
    assert meta["group_field"] == "patient"
    assert meta["tx_field"] == "treatment"


def test_parse_sample_meta_keeps_the_first_row_of_a_repeated_sample(tmp_path):
    path = write(tmp_path / "sheet.tsv", ANNOTATED_SHEET +
                 "S1\tRUN2\td_R1.fq.gz\td_R2.fq.gz\tP9\tPAT9\tother\t99\n")
    assert qc_parsers.parse_sample_meta(path)["rows"]["S1"]["patient"] == "PAT1"


def test_parse_sample_meta_of_the_committed_samplesheet_is_none(data_dir):
    # the fixture samplesheet carries only pipeline columns (runId / reads / reference / taxId),
    # so there is no user annotation to show and the panel self-hides
    assert qc_parsers.parse_sample_meta(str(data_dir / "samplesheet.tsv")) is None


@pytest.mark.parametrize("path", [None, "", MISSING])
def test_parse_sample_meta_of_an_absent_file_is_none(path):
    assert qc_parsers.parse_sample_meta(path) is None


def test_parse_sample_meta_without_rows_or_header_is_none(tmp_path):
    assert qc_parsers.parse_sample_meta(write(tmp_path / "sheet.tsv", "")) is None
    assert qc_parsers.parse_sample_meta(write(tmp_path / "h.tsv", "sampleId\tpassage\n")) is None


def test_parse_dose_reads_the_numeric_dose_column(tmp_path):
    # S3 has an empty dose and is simply absent from the map
    assert qc_parsers.parse_dose(write(tmp_path / "sheet.tsv", ANNOTATED_SHEET)) == {"S1": 10.0, "S2": 20.0}


def test_parse_dose_accepts_the_spanish_column_name(tmp_path):
    path = write(tmp_path / "sheet.tsv", "sampleId\tdosis\nS1\t2.5\n")
    assert qc_parsers.parse_dose(path) == {"S1": 2.5}


def test_parse_dose_without_a_dose_column_is_empty(tmp_path, data_dir):
    assert qc_parsers.parse_dose(write(tmp_path / "sheet.tsv", "sampleId\tpassage\nS1\tP1\n")) == {}
    assert qc_parsers.parse_dose(str(data_dir / "samplesheet.tsv")) == {}


@pytest.mark.parametrize("path", [None, "", MISSING])
def test_parse_dose_of_an_absent_file_is_empty(path):
    assert qc_parsers.parse_dose(path) == {}


def test_parse_metadata_auto_detects_the_time_and_group_columns(tmp_path):
    meta = qc_parsers.parse_metadata(write(tmp_path / "sheet.tsv", ANNOTATED_SHEET))
    assert meta == {"S1": {"time": "P1", "tnum": 1.0, "group": "PAT1"},
                    "S2": {"time": "P2", "tnum": 2.0, "group": "PAT1"},
                    "S3": {"time": "P1", "tnum": 1.0, "group": "PAT2"}}


def test_parse_metadata_of_the_committed_samplesheet_has_no_time_or_group(data_dir):
    # one entry per sample (the two SAMPLE_A runs collapse to the first), but nothing to plot
    meta = qc_parsers.parse_metadata(str(data_dir / "samplesheet.tsv"))
    assert sorted(meta) == ["SAMPLE_A", "SAMPLE_B", "SAMPLE_C"]
    assert all(v == {"time": None, "tnum": None, "group": None} for v in meta.values())


def test_parse_metadata_skips_comments_and_blank_lines(tmp_path):
    path = write(tmp_path / "sheet.tsv",
                 "# a comment\n\nsample\ttimepoint\tpatient\nS1\tday 3\tPAT1\n")
    assert qc_parsers.parse_metadata(path) == {"S1": {"time": "day 3", "tnum": 3.0, "group": "PAT1"}}


@pytest.mark.parametrize("path", [None, "", MISSING])
def test_parse_metadata_of_an_absent_file_is_empty(path):
    assert qc_parsers.parse_metadata(path) == {}


def test_parse_metadata_of_an_empty_file_is_empty(tmp_path):
    assert qc_parsers.parse_metadata(write(tmp_path / "sheet.tsv", "")) == {}


@pytest.mark.parametrize(("parser", "empty"), [
    (qc_parsers.parse_sample_meta, None), (qc_parsers.parse_dose, {}), (qc_parsers.parse_metadata, {})])
def test_samplesheet_parsers_skip_a_row_shorter_than_the_sample_column(tmp_path, parser, empty):
    path = write(tmp_path / "sheet.tsv", "# header comment\nrunId\tsampleId\tdose\n\nRUN1\n")
    assert parser(path) == empty


def test_parse_sample_meta_skips_a_row_shorter_than_the_sample_column(tmp_path):
    path = write(tmp_path / "sheet.tsv", "runId\tsampleId\tpassage\nRUN1\nRUN2\tS2\tP2\n")
    assert list(qc_parsers.parse_sample_meta(path)["rows"]) == ["S2"]


def test_parse_dose_of_an_empty_file_is_empty(tmp_path):
    assert qc_parsers.parse_dose(write(tmp_path / "sheet.tsv", "")) == {}


def test_parse_dose_keeps_the_first_row_of_a_repeated_sample(tmp_path):
    # the samplesheet has one row per run, so a multi-run sample appears more than once
    path = write(tmp_path / "sheet.tsv",
                 "sampleId\trunId\tdose\nS1\tRUN1\t10\nS1\tRUN2\t99\n\tRUN3\t5\n")
    assert qc_parsers.parse_dose(path) == {"S1": 10.0}


def test_parse_dose_skips_a_row_without_a_dose_cell(tmp_path):
    path = write(tmp_path / "sheet.tsv", "sampleId\tpassage\tdose\nS1\tP1\nS2\tP2\t5\n")
    assert qc_parsers.parse_dose(path) == {"S2": 5.0}


# --------------------------------------------------------------------------- parse_vcfs


ANN = ("ANN=T|missense_variant|MODERATE|rpoB|GENE_rpoB|transcript|TX1|protein_coding|"
       "1/1|c.1349C>T|p.Ser450Leu|450/3519|")
VCF = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=chr,length=5000>\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
    f"chr\t100\t.\tC\tT\t60\tPASS\t{ANN}\tGT:AD\t0/1:2,8\n"
    "chr\t200\t.\tA\tAT\t60\tPASS\t.\tGT:AD\t1/1:0,10\n"      # an insertion: not a SNP
    "chr\t300\t.\tAC\tA\t60\tPASS\t.\tGT:AD\t1/1:0,10\n"      # a deletion: not a SNP
    "chr\t400\t.\tA\t.\t60\tPASS\t.\tGT:AD\t0/0:10,0\n"       # no alt allele
    "chr\t500\t.\tG\tA\t60\tPASS\t.\tGT:DP\t./.:30\n"         # no genotype -> no frequency
    "chr\t600\t.\tG\tA\t60\tPASS\t.\tGT:DP:AD\t1/1:33:1,32\n"
)


def test_parse_vcfs_keeps_only_snps_with_a_usable_frequency(tmp_path):
    variants = qc_parsers.parse_vcfs([write(tmp_path / "S1.vcf", VCF)])
    assert list(variants) == ["S1"]
    assert sorted(variants["S1"]) == ["chr:100", "chr:600"]


def test_parse_vcfs_decodes_the_frequency_depth_and_annotation(tmp_path):
    variants = qc_parsers.parse_vcfs([write(tmp_path / "S1.vcf", VCF)])
    assert variants["S1"]["chr:100"] == {
        "ref": "C", "alt": "T", "gene": "rpoB", "eff": "missense_variant", "imp": "MODERATE",
        "aa": "p.Ser450Leu", "aa_h37rv": "", "af": 0.8, "dp": 10, "opos": None}
    # AD wins over the 1/1 genotype dosage, and the explicit DP wins over sum(AD)
    assert variants["S1"]["chr:600"]["af"] == 0.9697
    assert variants["S1"]["chr:600"]["dp"] == 33
    assert variants["S1"]["chr:600"]["gene"] == ""


def test_parse_vcfs_takes_the_sample_from_the_chrom_header_not_the_file_name(tmp_path):
    path = write(tmp_path / "unrelated_name.vcf", VCF.replace("\tS1\n", "\tREAL_ID\n"))
    assert list(qc_parsers.parse_vcfs([path])) == ["REAL_ID"]


def test_parse_vcfs_falls_back_to_the_file_stem_for_a_sampleless_vcf(tmp_path):
    path = write(tmp_path / "SAMPLE_Z.vcf",
                 "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                 "chr\t100\t.\tC\tT\t60\tPASS\t.\n")
    variants = qc_parsers.parse_vcfs([path])
    assert list(variants) == ["SAMPLE_Z"]
    # with no FORMAT/sample column the record is taken at face value: fixed, unknown depth
    assert variants["SAMPLE_Z"]["chr:100"]["af"] == 1.0
    assert variants["SAMPLE_Z"]["chr:100"]["dp"] is None


def test_parse_vcfs_reads_a_gzipped_vcf(tmp_path):
    assert qc_parsers.parse_vcfs([write_gz(tmp_path / "S1.vcf.gz", VCF)]) == \
        qc_parsers.parse_vcfs([write(tmp_path / "S1.vcf", VCF)])


def test_parse_vcfs_records_the_liftover_original_position(tmp_path):
    path = write(tmp_path / "S1.vcf",
                 "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
                 "chr\t100\t.\tC\tT\t60\tPASS\tOriginalContig=used;OriginalStart=77\tGT\t1/1\n")
    assert qc_parsers.parse_vcfs([path])["S1"]["chr:100"]["opos"] == "used:77"


def test_parse_vcfs_keeps_the_first_alt_of_a_multiallelic_snp(tmp_path):
    path = write(tmp_path / "S1.vcf",
                 "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
                 "chr\t100\t.\tC\tT,G\t60\tPASS\t.\tGT:AD\t1/2:1,5,4\n")
    record = qc_parsers.parse_vcfs([path])["S1"]["chr:100"]
    assert record["alt"] == "T"
    assert record["af"] == 0.9


def test_parse_vcfs_merges_several_files(tmp_path):
    a = write(tmp_path / "a.vcf", VCF)
    b = write(tmp_path / "b.vcf", VCF.replace("\tS1\n", "\tS2\n"))
    assert sorted(qc_parsers.parse_vcfs([a, b])) == ["S1", "S2"]


def test_parse_vcfs_without_a_chrom_header_yields_nothing(tmp_path):
    path = write(tmp_path / "S1.vcf", "##fileformat=VCFv4.2\nchr\t100\t.\tC\tT\t60\tPASS\t.\n")
    assert qc_parsers.parse_vcfs([path]) == {}


def test_parse_vcfs_skips_a_truncated_data_line(tmp_path):
    path = write(tmp_path / "S1.vcf",
                 "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
                 "chr\t100\t.\tC\n"
                 "chr\t200\t.\tC\tT\t60\tPASS\t.\tGT\t1/1\n")
    assert list(qc_parsers.parse_vcfs([path])["S1"]) == ["chr:200"]


@pytest.mark.parametrize("paths", [None, [], [MISSING], ["", None]])
def test_parse_vcfs_without_any_readable_file_is_empty(paths):
    assert qc_parsers.parse_vcfs(paths) == {}


# --------------------------------------------------------------------------- mapDamage readers


FREQ = "pos\t5pC>T\n1\t0.25\n2\t0.12\n3\t0.05\n"
LGDIST = "# table produced by mapDamage\nStd\tLength\tOccurences\n-\t50\t10\n-\t100\t10\n"


def test_read_freq_file_reads_the_terminal_damage_profile(tmp_path):
    freq = qc_parsers._read_freq_file(write(tmp_path / "5pCtoT_freq.txt", FREQ))
    assert freq["pos1"] == pytest.approx(0.25)
    assert freq["profile"] == [[1, 0.25], [2, 0.12], [3, 0.05]]


def test_read_freq_file_accepts_comma_and_space_separated_rows(tmp_path):
    comma = qc_parsers._read_freq_file(write(tmp_path / "a.txt", "1,0.25\n2,0.12\n"))
    spaced = qc_parsers._read_freq_file(write(tmp_path / "b.txt", "1 0.25\n2 0.12\n"))
    assert comma == spaced
    assert comma["pos1"] == pytest.approx(0.25)


def test_read_freq_file_sorts_by_position_and_caps_the_profile(tmp_path):
    body = "".join(f"{p}\t0.1\n" for p in range(30, 0, -1))
    freq = qc_parsers._read_freq_file(write(tmp_path / "freq.txt", body))
    assert len(freq["profile"]) == 25
    assert freq["profile"][0][0] == 1


def test_read_freq_file_honours_max_pos(tmp_path):
    freq = qc_parsers._read_freq_file(write(tmp_path / "freq.txt", FREQ), max_pos=2)
    assert len(freq["profile"]) == 2


def test_read_freq_file_falls_back_to_the_first_row_without_a_position_one(tmp_path):
    freq = qc_parsers._read_freq_file(write(tmp_path / "freq.txt", "3\t0.05\n2\t0.12\n"))
    assert freq["pos1"] == pytest.approx(0.12)


def test_read_freq_file_skips_malformed_rows(tmp_path):
    freq = qc_parsers._read_freq_file(write(tmp_path / "freq.txt",
                                    "1\t0.25\nonlyonefield\nx\ty\n2\tNA\n2\t0.12\n"))
    assert freq["profile"] == [[1, 0.25], [2, 0.12]]


def test_read_freq_file_of_an_empty_or_missing_file_is_none(tmp_path):
    assert qc_parsers._read_freq_file(write(tmp_path / "freq.txt", "")) is None
    assert qc_parsers._read_freq_file(write(tmp_path / "h.txt", "pos\tfreq\n")) is None
    assert qc_parsers._read_freq_file(str(tmp_path / "nope.txt")) is None


def test_mean_fraglen_is_the_occurrence_weighted_mean_length(tmp_path):
    assert qc_parsers._mean_fraglen(write(tmp_path / "lgdistribution.txt", LGDIST)) == pytest.approx(75.0)


def test_mean_fraglen_skips_comments_headers_and_short_rows(tmp_path):
    path = write(tmp_path / "lg.txt", LGDIST + "-\t999\n" + "# trailing comment\n")
    assert qc_parsers._mean_fraglen(path) == pytest.approx(75.0)


def test_mean_fraglen_skips_a_row_with_a_non_numeric_length_or_count(tmp_path):
    path = write(tmp_path / "lg.txt", LGDIST + "-\tlong\tmany\n-\t50\tNA\n")
    assert qc_parsers._mean_fraglen(path) == pytest.approx(75.0)


def test_mean_fraglen_without_any_occurrence_is_none(tmp_path):
    assert qc_parsers._mean_fraglen(write(tmp_path / "lg.txt", "-\t50\t0\n")) is None
    assert qc_parsers._mean_fraglen(write(tmp_path / "empty.txt", "")) is None
    assert qc_parsers._mean_fraglen(str(tmp_path / "nope.txt")) is None


def test_mapdamage_stats_reads_a_complete_sample_directory(tmp_path):
    d = tmp_path / "S1"
    d.mkdir()
    write(d / "5pCtoT_freq.txt", FREQ)
    write(d / "3pGtoA_freq.txt", "pos\t3pG>A\n1\t0.18\n2\t0.09\n")
    write(d / "lgdistribution.txt", LGDIST)
    stats = qc_parsers.mapdamage_stats(str(d))
    assert stats["ct1"] == pytest.approx(0.25)
    assert stats["ga1"] == pytest.approx(0.18)
    assert stats["fraglen"] == pytest.approx(75.0)
    assert stats["ct_profile"] == [[1, 0.25], [2, 0.12], [3, 0.05]]


def test_mapdamage_stats_tolerates_a_partial_directory(tmp_path):
    d = tmp_path / "S1"
    d.mkdir()
    write(d / "lgdistribution.txt", LGDIST)
    assert qc_parsers.mapdamage_stats(str(d)) == {"ct1": None, "ga1": None,
                                          "fraglen": pytest.approx(75.0), "ct_profile": None}


def test_mapdamage_stats_of_an_empty_or_absent_directory_is_none(tmp_path):
    empty = tmp_path / "empty"
    empty.mkdir()
    assert qc_parsers.mapdamage_stats(str(empty)) is None
    assert qc_parsers.mapdamage_stats(str(tmp_path / "nope")) is None
    assert qc_parsers.mapdamage_stats(None) is None
    # a file where a directory is expected is not a mapDamage output
    assert qc_parsers.mapdamage_stats(write(tmp_path / "a-file.txt", FREQ)) is None


# --------------------------------------------------------------------------- unreadable optional inputs


@pytest.mark.parametrize(("parser", "empty"), [
    (qc_parsers.parse_lineage_colors, {}),
    (qc_parsers.parse_bed, []),
    (qc_parsers.parse_gene_burden, []),
    (qc_parsers.parse_pnps, []),
    (qc_parsers.parse_dr, None),
    (qc_parsers.parse_gene_conversion, None),
    (qc_parsers.parse_gff, []),
    (qc_parsers.parse_gene_locus, {}),
    (qc_parsers.parse_sample_meta, None),
    (qc_parsers.parse_dose, {}),
    (qc_parsers.parse_metadata, {}),
])
def test_optional_parsers_degrade_quietly_on_an_unreadable_path(tmp_path, parser, empty):
    """A path that exists but cannot be opened (here a directory) hides its panel rather
    than crashing the whole report. Only the two required inputs, parse_summary and
    consensus_stats, are allowed to raise."""
    assert parser(str(tmp_path)) == empty


@pytest.mark.parametrize(("parser", "empty"), [(qc_parsers.parse_kraken, None), (qc_parsers.parse_vcfs, {})])
def test_multi_file_parsers_skip_an_unreadable_path(tmp_path, parser, empty):
    assert parser([str(tmp_path)]) == empty
