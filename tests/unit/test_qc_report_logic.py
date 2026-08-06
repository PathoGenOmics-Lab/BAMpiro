"""Unit tests for the QC verdict engine and the derived statistics of bin/qc_report.py.

The file parsers live in test_qc_report_parsers.py; this module covers the pure
in-memory logic: flag_sample, the robust statistics, the lineage/het derivations,
the small VCF field decoders and the three payload builders.

Every expectation here is pinned against the CURRENT behaviour of the script, and
the thresholds are read from qc_report.DEF / ANC_DEF rather than hardcoded, so a
deliberate change to a default shows up as a behaviour change and not as a
spurious failure.
"""

from __future__ import annotations

import pytest
from conftest import load_script

# qc_report reads bin/report_assets/** at import time, so it must be imported in place.
qc = load_script("qc_report")

THR = dict(qc.DEF)
ANC = dict(qc.ANC_DEF)


def metrics(**over):
    """A summary row that PASSes every default check, overridable per test."""
    m = {
        "mean_depth": "50",
        "breadth_pct": "99",
        "missing_pct": "1",
        "mapped_pct": "95",
        "duplication_pct": "5",
        "iupac_pct": "0.5",
        "ti_tv": "1.8",
        "snps": "1000",
    }
    m.update(over)
    return m


# --------------------------------------------------------------------------- flag_sample: verdicts


def test_flag_sample_clean_row_passes_with_no_flags():
    assert qc.flag_sample(metrics(), THR, None, None) == ("PASS", [])


def test_flag_sample_all_missing_metrics_is_a_hard_no_data_fail():
    assert qc.flag_sample({}, THR, None, None) == ("FAIL", ["NO_DATA"])


def test_flag_sample_na_sentinels_count_as_missing_for_the_no_data_guard():
    row = {"mean_depth": "NA", "breadth_pct": "", "missing_pct": ".", "mapped_pct": None}
    assert qc.flag_sample(row, THR, None, None) == ("FAIL", ["NO_DATA"])


def test_flag_sample_no_data_guard_only_looks_at_its_four_keys():
    # duplication_pct alone is not evidence that the sample produced anything, so the
    # guard still fires and the HIGH_DUP check is never reached.
    assert qc.flag_sample({"duplication_pct": "99"}, THR, None, None) == ("FAIL", ["NO_DATA"])


@pytest.mark.parametrize(
    ("key", "value", "flag"),
    [
        ("mean_depth", qc.DEF["depth_min"] - 1, "LOW_DEPTH"),
        ("breadth_pct", qc.DEF["breadth_min"] - 1, "LOW_BREADTH"),
        ("missing_pct", qc.DEF["missing_max"] + 1, "HIGH_MISSING"),
    ],
)
def test_flag_sample_fail_flags_drive_a_fail_verdict(key, value, flag):
    verdict, flags = qc.flag_sample(metrics(**{key: value}), THR, None, None)
    assert flags == [flag]
    assert flag in qc.FAIL_FLAGS
    assert verdict == "FAIL"


@pytest.mark.parametrize(
    ("key", "value", "flag"),
    [
        ("mapped_pct", qc.DEF["mapping_min"] - 1, "MAPPING_LOW"),
        ("duplication_pct", qc.DEF["dup_max"] + 1, "HIGH_DUP"),
        ("iupac_pct", qc.DEF["iupac_max"] + 1, "HIGH_IUPAC"),
    ],
)
def test_flag_sample_soft_flags_only_warn(key, value, flag):
    verdict, flags = qc.flag_sample(metrics(**{key: value}), THR, None, None)
    assert flags == [flag]
    assert flag not in qc.FAIL_FLAGS
    assert verdict == "WARN"


@pytest.mark.parametrize(
    ("key", "value"),
    [
        ("mean_depth", qc.DEF["depth_min"]),      # strict <, so exactly at the floor passes
        ("breadth_pct", qc.DEF["breadth_min"]),   # strict <
        ("missing_pct", qc.DEF["missing_max"]),   # strict >, so exactly at the ceiling passes
        ("mapped_pct", qc.DEF["mapping_min"]),    # strict <
        ("duplication_pct", qc.DEF["dup_max"]),   # strict >
        ("iupac_pct", qc.DEF["iupac_max"]),       # strict >
    ],
)
def test_flag_sample_threshold_boundaries_are_inclusive_of_the_threshold(key, value):
    assert qc.flag_sample(metrics(**{key: value}), THR, None, None) == ("PASS", [])


def test_flag_sample_multiple_problems_accumulate_and_a_single_fail_flag_wins():
    verdict, flags = qc.flag_sample(
        metrics(mean_depth="1", duplication_pct="99", iupac_pct="50"), THR, None, None)
    assert flags == ["LOW_DEPTH", "HIGH_DUP", "HIGH_IUPAC"]
    assert verdict == "FAIL"


def test_flag_sample_missing_single_metric_is_simply_not_checked():
    row = metrics()
    del row["duplication_pct"]
    assert qc.flag_sample(row, THR, None, None) == ("PASS", [])


# --------------------------------------------------------------------------- flag_sample: Ti/Tv


def test_flag_sample_titv_check_is_disabled_at_the_default_threshold_of_zero():
    assert qc.DEF["titv_min"] == 0.0
    assert qc.flag_sample(metrics(ti_tv="0.3"), THR, None, None) == ("PASS", [])


def test_flag_sample_titv_low_warns_once_the_threshold_is_raised():
    thr = dict(THR, titv_min=1.5)
    assert qc.flag_sample(metrics(ti_tv="1.2"), thr, None, None) == ("WARN", ["TITV_LOW"])
    assert qc.flag_sample(metrics(ti_tv="1.5"), thr, None, None) == ("PASS", [])


# --------------------------------------------------------------------------- flag_sample: SNP outliers


def test_flag_sample_snp_low_and_high_are_cohort_z_outliers():
    med, sig = 1000.0, 10.0
    span = THR["snp_z"] * sig
    assert qc.flag_sample(metrics(snps=med - span - 1), THR, med, sig) == ("WARN", ["SNP_LOW"])
    assert qc.flag_sample(metrics(snps=med + span + 1), THR, med, sig) == ("WARN", ["SNP_HIGH"])
    # both comparisons are strict, so a sample exactly z sigma out is not an outlier
    assert qc.flag_sample(metrics(snps=med - span), THR, med, sig) == ("PASS", [])
    assert qc.flag_sample(metrics(snps=med + span), THR, med, sig) == ("PASS", [])


def test_flag_sample_snp_check_is_skipped_without_a_cohort_scale():
    assert qc.flag_sample(metrics(snps="999999"), THR, 1000.0, None) == ("PASS", [])
    assert qc.flag_sample(metrics(snps="999999"), THR, 1000.0, 0.0) == ("PASS", [])
    assert qc.flag_sample(metrics(snps="999999"), THR, None, 10.0) == ("PASS", [])


# --------------------------------------------------------------------------- flag_sample: het fraction


def test_flag_sample_het_high_uses_the_het_fraction_as_a_percent():
    # 5 het / 100 variant genotypes = 5%, above the 3% default
    row = metrics(het_variants="5", homo_variants="95")
    assert qc.flag_sample(row, THR, None, None) == ("WARN", ["HET_HIGH"])
    assert qc.flag_sample(metrics(het_variants="2", homo_variants="98"), THR, None, None) == ("PASS", [])


def test_flag_sample_het_check_is_skipped_below_the_minimum_variant_count():
    # 50% het, but only 10 variant genotypes -> het_frac() returns None, so no flag
    row = metrics(het_variants="5", homo_variants="5")
    assert qc.flag_sample(row, THR, None, None) == ("PASS", [])


# --------------------------------------------------------------------------- flag_sample: MIXED


def test_flag_sample_mixed_needs_two_lineages_over_both_the_fraction_and_the_raw_count():
    row = metrics(lineage_counts="L4:100;L2:10")
    assert qc.flag_sample(row, THR, None, None) == ("WARN", ["MIXED"])


def test_flag_sample_mixed_ignores_a_minor_lineage_below_three_markers():
    # 2 markers is 1.96% of 102 and also below the >= 3 raw-count guard
    assert qc.flag_sample(metrics(lineage_counts="L4:100;L2:2"), THR, None, None) == ("PASS", [])


def test_flag_sample_mixed_ignores_a_minor_lineage_below_the_fraction_cutoff():
    # 3 markers clears the raw count but is 0.03% of the total, below mixed_min_frac
    assert qc.flag_sample(metrics(lineage_counts="L4:10000;L2:3"), THR, None, None) == ("PASS", [])


def test_flag_sample_mixed_is_not_raised_by_two_sublineages_of_the_same_lineage():
    assert qc.flag_sample(metrics(lineage_counts="L4.1:100;L4.3:100"), THR, None, None) == ("PASS", [])


def test_flag_sample_mixed_is_skipped_when_there_are_no_lineage_counts():
    assert qc.flag_sample(metrics(lineage_counts="NA"), THR, None, None) == ("PASS", [])


# --------------------------------------------------------------------------- flag_sample: ancient view


def test_flag_sample_ancient_row_is_judged_against_the_adna_thresholds():
    # 5x depth / 40% breadth / 50% missing FAILs the modern gate but is fine for aDNA
    row = metrics(mean_depth="5", breadth_pct="40", missing_pct="50", mapped_pct="25",
                  duplication_pct="55", iupac_pct="4")
    assert qc.flag_sample(row, THR, None, None)[0] == "FAIL"
    assert qc.flag_sample(row, THR, None, None, ancient=True, anc_thr=ANC) == ("PASS", [])


def test_flag_sample_ancient_falls_back_to_the_modern_thresholds_without_anc_thr():
    row = metrics(mean_depth="5")
    assert qc.flag_sample(row, THR, None, None, ancient=True, anc_thr=None) == ("FAIL", ["LOW_DEPTH"])


def test_flag_sample_ancient_skips_the_snp_outlier_check():
    row = metrics(snps="1")
    assert qc.flag_sample(row, THR, 1000.0, 10.0) == ("WARN", ["SNP_LOW"])
    assert qc.flag_sample(row, THR, 1000.0, 10.0, ancient=True, anc_thr=ANC) == ("PASS", [])


def test_flag_sample_ancient_skips_the_het_check_because_deamination_inflates_it():
    row = metrics(het_variants="40", homo_variants="60")
    assert qc.flag_sample(row, THR, None, None) == ("WARN", ["HET_HIGH"])
    assert qc.flag_sample(row, THR, None, None, ancient=True, anc_thr=ANC) == ("PASS", [])


def test_flag_sample_damage_low_flags_an_ancient_sample_under_the_deamination_floor():
    dmg = {"ct1": ANC["damage_min_ct"] / 2.0}
    verdict, flags = qc.flag_sample(metrics(), THR, None, None, ancient=True, anc_thr=ANC, dmg=dmg)
    assert (verdict, flags) == ("WARN", ["DAMAGE_LOW"])


def test_flag_sample_damage_low_is_not_raised_at_or_above_the_floor():
    for ct1 in (ANC["damage_min_ct"], 0.40):
        assert qc.flag_sample(metrics(), THR, None, None,
                              ancient=True, anc_thr=ANC, dmg={"ct1": ct1}) == ("PASS", [])


@pytest.mark.parametrize("dmg", [None, {}, {"ct1": None}])
def test_flag_sample_damage_check_needs_a_measured_ct1(dmg):
    assert qc.flag_sample(metrics(), THR, None, None, ancient=True, anc_thr=ANC, dmg=dmg) == ("PASS", [])


def test_flag_sample_damage_check_is_not_applied_to_a_modern_sample():
    dmg = {"ct1": 0.0}
    assert qc.flag_sample(metrics(), THR, None, None, ancient=False, anc_thr=ANC, dmg=dmg) == ("PASS", [])


# --------------------------------------------------------------------------- flag_sample: threshold dict contract


@pytest.mark.parametrize("missing_key", sorted(qc.DEF))
def test_flag_sample_raises_key_error_when_the_threshold_dict_is_incomplete(missing_key):
    """flag_sample does no .get() on thr: a caller-built partial threshold dict is a hard error."""
    thr = {k: v for k, v in THR.items() if k != missing_key}
    row = metrics(het_variants="5", homo_variants="95", lineage_counts="L4:100;L2:10")
    with pytest.raises(KeyError) as exc:
        qc.flag_sample(row, thr, 1000.0, 10.0)
    assert missing_key in str(exc.value)


# --------------------------------------------------------------------------- robust


def test_robust_of_an_empty_or_all_none_sequence_is_a_pair_of_nones():
    assert qc.robust([]) == (None, None)
    assert qc.robust([None, None]) == (None, None)


def test_robust_of_a_single_value_falls_all_the_way_back_to_sigma_one():
    assert qc.robust([7.0]) == (7.0, 1.0)


def test_robust_sigma_is_never_zero_for_an_all_ties_cohort():
    # a zero sigma would make every downstream z-score divide by zero
    med, sigma = qc.robust([5.0] * 4)
    assert med == 5.0
    assert sigma == 1.0


def test_robust_uses_the_mad_scaled_by_the_normal_consistency_constant():
    med, sigma = qc.robust([1.0, 2.0, 3.0])
    assert med == 2.0
    assert sigma == pytest.approx(1.4826)


def test_robust_falls_back_to_the_population_sd_when_the_mad_is_zero():
    # [1, 1, 1, 5]: median 1, MAD 0 -> population sd of the raw values
    med, sigma = qc.robust([1.0, 1.0, 1.0, 5.0])
    assert med == 1.0
    assert sigma == pytest.approx(3.0 ** 0.5)


def test_robust_uses_the_even_length_midpoint_and_ignores_nones():
    med, sigma = qc.robust([3.0, None, 1.0, 2.0, None, 4.0])
    assert med == 2.5
    assert sigma == pytest.approx(1.4826)


def test_robust_median_is_resistant_to_a_large_outlier():
    med, _ = qc.robust([1.0, 2.0, 3.0, 1000.0])
    assert med == 2.5


# --------------------------------------------------------------------------- het_frac


def test_het_frac_returns_the_het_share_of_variant_genotypes():
    assert qc.het_frac({"het_variants": "10", "homo_variants": "90"}) == pytest.approx(0.1)


def test_het_frac_floor_is_inclusive_of_the_minimum_variant_count():
    n = qc.MIN_HET_VARIANTS
    assert qc.het_frac({"het_variants": "1", "homo_variants": n - 1}) is not None
    assert qc.het_frac({"het_variants": "1", "homo_variants": n - 2}) is None


@pytest.mark.parametrize("row", [
    {},
    {"het_variants": "10"},
    {"homo_variants": "90"},
    {"het_variants": "NA", "homo_variants": "90"},
])
def test_het_frac_needs_both_counts(row):
    assert qc.het_frac(row) is None


# --------------------------------------------------------------------------- lineage counts / fractions


def test_lineage_counts_parsed_collapses_sublineages_on_the_first_dot():
    assert qc.lineage_counts_parsed("L4.1:3200;L4.3:120;L2:40") == {"L4": 3320.0, "L2": 40.0}


def test_lineage_counts_parsed_treats_commas_as_semicolons():
    assert qc.lineage_counts_parsed("L4:10,L2:20") == {"L4": 10.0, "L2": 20.0}


@pytest.mark.parametrize("value", [None, "", "NA", ".", 0])
def test_lineage_counts_parsed_of_an_absent_value_is_empty(value):
    assert qc.lineage_counts_parsed(value) == {}


def test_lineage_counts_parsed_drops_malformed_and_non_numeric_entries():
    assert qc.lineage_counts_parsed("L4;L2:20;L3:abc") == {"L2": 20.0}


def test_lineage_fracs_normalises_the_collapsed_counts():
    assert qc.lineage_fracs("L4.1:30;L4.2:10;L2:10") == {"L4": 0.8, "L2": 0.2}


def test_lineage_fracs_of_a_zero_total_is_empty_rather_than_a_division_by_zero():
    assert qc.lineage_fracs("L4:0") == {}
    assert qc.lineage_fracs("NA") == {}


# --------------------------------------------------------------------------- is_ancient


@pytest.mark.parametrize("value", ["yes", "Yes", " YES ", "y", "true", "True", "1", "ancient"])
def test_is_ancient_accepts_the_documented_truthy_spellings(value):
    assert qc.is_ancient({"ancient": value}) is True


@pytest.mark.parametrize("value", ["no", "", "NA", "0", "false", "modern", None])
def test_is_ancient_is_false_for_anything_else(value):
    assert qc.is_ancient({"ancient": value}) is False


def test_is_ancient_of_a_row_without_the_column_is_modern():
    assert qc.is_ancient({}) is False


# --------------------------------------------------------------------------- infer_metric / discover_extra_metrics


@pytest.mark.parametrize("key", ["foo_pct", "error_rate", "some_percent", "alt_fraction", "af_frac"])
def test_infer_metric_recognises_a_percentage_column_by_name(key):
    _, kind, direction = qc.infer_metric(key, [1.5])
    assert kind == "pct"
    assert direction == "neu"


def test_infer_metric_recognises_an_integer_column_by_name_and_by_value():
    assert qc.infer_metric("custom_reads", [10.0, 20.0]) == ("Custom reads", "int", "neu")
    # the same name with a non-integral value is not an int column
    assert qc.infer_metric("custom_reads", [10.5]) == ("Custom reads", "float", "neu")


def test_infer_metric_defaults_to_float_and_prettifies_the_label():
    assert qc.infer_metric("weird_thing", [1.5]) == ("Weird thing", "float", "neu")
    assert qc.infer_metric("foo_pct", [1.5])[0] == "Foo %"


def test_infer_metric_never_claims_a_direction():
    # an auto-detected column must not silently start flagging samples
    for key in ("foo_pct", "custom_reads", "weird_thing"):
        assert qc.infer_metric(key, [1.0])[2] == "neu"


def test_discover_extra_metrics_of_an_empty_cohort_is_empty():
    assert qc.discover_extra_metrics({}) == []
    assert qc.discover_extra_metrics(None) == []


def test_discover_extra_metrics_picks_up_only_unknown_numeric_columns():
    summ = {
        "S1": {"sample_id": "S1", "mean_depth": "30", "lineage": "L4",
               "custom_reads": "100", "note": "hello"},
        "S2": {"sample_id": "S2", "mean_depth": "40", "lineage": "L2",
               "custom_reads": "200", "note": "world"},
    }
    assert qc.discover_extra_metrics(summ) == [
        {"key": "custom_reads", "label": "Custom reads", "kind": "int", "dir": "neu"}]


def test_discover_extra_metrics_needs_only_one_numeric_sample():
    summ = {"S1": {"sample_id": "S1", "odd_metric": "NA"},
            "S2": {"sample_id": "S2", "odd_metric": "3.5"}}
    assert [e["key"] for e in qc.discover_extra_metrics(summ)] == ["odd_metric"]


def test_discover_extra_metrics_is_sorted_by_key():
    summ = {"S1": {"sample_id": "S1", "zeta_metric": "1", "alpha_metric": "2"}}
    assert [e["key"] for e in qc.discover_extra_metrics(summ)] == ["alpha_metric", "zeta_metric"]


# --------------------------------------------------------------------------- to_float / clean_str


@pytest.mark.parametrize("value", [None, "", "NA", "."])
def test_to_float_maps_the_missing_sentinels_to_none(value):
    assert qc.to_float(value) is None


@pytest.mark.parametrize(("value", "expected"), [
    ("1.5", 1.5), ("0", 0.0), (7, 7.0), ("-3", -3.0), ("1e3", 1000.0), (" 2.5 ", 2.5),
])
def test_to_float_parses_a_numeric_value(value, expected):
    assert qc.to_float(value) == pytest.approx(expected)


@pytest.mark.parametrize("value", ["abc", "1,5", [1], {}])
def test_to_float_of_an_unparseable_value_is_none_rather_than_a_raise(value):
    assert qc.to_float(value) is None


@pytest.mark.parametrize("value", [None, "", "   ", "NA", ".", "-"])
def test_clean_str_maps_the_metadata_sentinels_to_none(value):
    assert qc.clean_str(value) is None


def test_clean_str_strips_a_real_value():
    assert qc.clean_str("  L4.3.4  ") == "L4.3.4"
    assert qc.clean_str(4) == "4"


# --------------------------------------------------------------------------- parse_profile


@pytest.mark.parametrize("value", [None, "", "NA", ".", 0])
def test_parse_profile_of_an_absent_value_is_none(value):
    assert qc.parse_profile(value) is None


@pytest.mark.parametrize("value", ["a,b,c", "1;2;3"])
def test_parse_profile_of_an_unparseable_value_is_none(value):
    assert qc.parse_profile(value) is None


def test_parse_profile_keeps_a_correctly_binned_profile_verbatim():
    prof = qc.parse_profile(",".join(str(i) for i in range(qc.NBINS)))
    assert prof == list(range(qc.NBINS))


def test_parse_profile_truncates_floats_to_ints():
    qc._PROFILE_WARNED.clear()
    assert qc.parse_profile("1.7,2.9,0.4")[:3] == [1, 2, 0]


def test_parse_profile_pads_a_short_profile_to_nbins():
    qc._PROFILE_WARNED.clear()
    prof = qc.parse_profile("5,6,7")
    assert len(prof) == qc.NBINS
    assert prof[:3] == [5, 6, 7]
    assert set(prof[3:]) == {0}


def test_parse_profile_truncates_a_long_profile_to_nbins():
    qc._PROFILE_WARNED.clear()
    prof = qc.parse_profile(",".join("1" for _ in range(qc.NBINS + 50)))
    assert prof == [1] * qc.NBINS


def test_parse_profile_warns_once_per_observed_bin_count(capsys):
    qc._PROFILE_WARNED.clear()
    qc.parse_profile("1,2,3")
    qc.parse_profile("4,5,6")          # same length -> already warned, stays quiet
    qc.parse_profile("1,2,3,4")        # a new length warns again
    err = capsys.readouterr().err
    assert err.count("variant profile has 3 bins") == 1
    assert err.count("variant profile has 4 bins") == 1
    qc._PROFILE_WARNED.clear()


# --------------------------------------------------------------------------- mask_profile


def test_mask_profile_without_intervals_or_genome_length_is_none():
    assert qc.mask_profile([], 100) == (None, 0.0)
    assert qc.mask_profile(None, 100) == (None, 0.0)
    assert qc.mask_profile([(0, 10)], 0) == (None, 0.0)


def test_mask_profile_marks_the_covered_bins_and_the_total_percent():
    prof, pct = qc.mask_profile([(0, 50)], 100, nbins=10)
    assert prof == [1.0] * 5 + [0.0] * 5
    assert pct == 50.0


def test_mask_profile_clamps_intervals_to_the_genome():
    prof, pct = qc.mask_profile([(-10, 200)], 100, nbins=10)
    assert prof == [1.0] * 10
    assert pct == 100.0


def test_mask_profile_skips_an_empty_or_inverted_interval():
    prof, pct = qc.mask_profile([(50, 50), (80, 20)], 100, nbins=10)
    assert prof == [0.0] * 10
    assert pct == 0.0


def test_mask_profile_reports_a_partially_covered_bin_as_a_fraction():
    prof, pct = qc.mask_profile([(0, 5)], 100, nbins=10)
    assert prof[0] == 0.5
    assert prof[1] == 0.0
    assert pct == 5.0


def test_mask_profile_defaults_to_the_module_bin_count():
    prof, _ = qc.mask_profile([(0, 1000)], 1000)
    assert len(prof) == qc.NBINS


def test_mask_profile_double_counts_overlapping_intervals_in_the_total():
    # KNOWN BEHAVIOUR (pinned, not endorsed): the per-bin profile is clamped to 1.0 but the
    # masked TOTAL sums interval lengths without merging, so an un-merged BED with overlapping
    # records reports a masked % above the true union (here 100% of the genome twice = 200%).
    prof, pct = qc.mask_profile([(0, 100), (0, 100)], 100, nbins=10)
    assert prof == [1.0] * 10
    assert pct == 200.0


# --------------------------------------------------------------------------- _pearson


def test_pearson_of_a_perfectly_correlated_pair_is_one():
    assert qc._pearson([1.0, 2.0, 3.0], [2.0, 4.0, 6.0]) == pytest.approx(1.0)


def test_pearson_of_a_perfectly_anticorrelated_pair_is_minus_one():
    assert qc._pearson([1.0, 2.0, 3.0], [6.0, 4.0, 2.0]) == pytest.approx(-1.0)


def test_pearson_needs_at_least_three_points():
    assert qc._pearson([1.0, 2.0], [1.0, 2.0]) is None
    assert qc._pearson([], []) is None


def test_pearson_of_unequal_lengths_is_none():
    assert qc._pearson([1.0, 2.0, 3.0], [1.0, 2.0, 3.0, 4.0]) is None


def test_pearson_of_a_flat_vector_is_none():
    assert qc._pearson([1.0, 1.0, 1.0], [1.0, 2.0, 3.0]) is None
    assert qc._pearson([1.0, 2.0, 3.0], [5.0, 5.0, 5.0]) is None


# --------------------------------------------------------------------------- _perm_p


PAIR = [([0.0, 0.5, 1.0, 0.5], [0.0, 0.5, 1.0, 0.5])]


def test_perm_p_is_reproducible_for_a_fixed_seed():
    assert qc._perm_p(PAIR, 1.0, B=100, seed=7) == qc._perm_p(PAIR, 1.0, B=100, seed=7)


def test_perm_p_is_the_plus_one_smoothed_fraction():
    p = qc._perm_p(PAIR, 1.0, B=100, seed=7)
    assert p == pytest.approx((round(p * 101) - 1 + 1) / 101)
    assert 1 / 101 <= p <= 1.0


def test_perm_p_bottoms_out_at_one_over_b_plus_one_when_no_null_reaches_the_observed():
    # |r| can never exceed 1, so an observed mean of 2 is unreachable
    assert qc._perm_p(PAIR, 2.0, B=100, seed=1) == pytest.approx(1 / 101)


def test_perm_p_is_one_when_every_null_reaches_the_observed():
    p = qc._perm_p([([0.0, 0.5, 1.0, 0.5], [0.1, 0.6, 0.2, 0.9])], 0.0, B=100, seed=1)
    assert p == pytest.approx(1.0)


def test_perm_p_with_undefined_correlations_reports_the_floor():
    # 2-point trajectories make every _pearson None -> no null ever counts
    assert qc._perm_p([([0.0, 1.0], [0.0, 1.0])], 0.0, B=100, seed=1) == pytest.approx(1 / 101)


# --------------------------------------------------------------------------- _dyn_af / _dyn_dp


def test_dyn_af_prefers_ad_over_ao_ro_and_over_the_gt_dosage():
    assert qc._dyn_af("GT:AD:AO:RO", "0/1:2,8:99:1") == pytest.approx(0.8)


def test_dyn_af_falls_back_to_ao_over_ao_plus_ro():
    assert qc._dyn_af("GT:AO:RO", "0/1:3:1") == pytest.approx(0.75)


def test_dyn_af_sums_every_alt_allele_of_a_multiallelic_ad():
    assert qc._dyn_af("GT:AD", "1/2:2,4,4") == pytest.approx(0.8)


@pytest.mark.parametrize(("gt", "expected"), [("1/1", 1.0), ("0/1", 0.5), ("0/0", 0.0), ("1|1", 1.0), ("1", 1.0)])
def test_dyn_af_falls_back_to_the_gt_dosage(gt, expected):
    assert qc._dyn_af("GT", gt) == pytest.approx(expected)


def test_dyn_af_of_a_no_call_genotype_is_none():
    assert qc._dyn_af("GT", "./.") is None
    assert qc._dyn_af("GT", ".|.") is None


def test_dyn_af_skips_an_unusable_ad_and_keeps_going():
    # a zero-depth AD, a single-entry AD and a non-numeric AD all fall through to the GT dosage
    assert qc._dyn_af("GT:AD", "1/1:0,0") == pytest.approx(1.0)
    assert qc._dyn_af("GT:AD", "0/1:10") == pytest.approx(0.5)
    assert qc._dyn_af("GT:AD", "0/1:x,y") == pytest.approx(0.5)


def test_dyn_af_of_a_zero_depth_ao_ro_falls_through_to_the_dosage():
    assert qc._dyn_af("GT:AO:RO", "0/1:0:0") == pytest.approx(0.5)


def test_dyn_af_of_a_non_numeric_ao_ro_falls_through_to_the_dosage():
    assert qc._dyn_af("GT:AO:RO", "1/1:x:1") == pytest.approx(1.0)


def test_dyn_dp_prefers_the_format_dp_field():
    assert qc._dyn_dp("GT:DP:AD", "0/1:37:2,8") == 37


def test_dyn_dp_falls_back_to_the_ad_sum_then_to_ao_plus_ro():
    assert qc._dyn_dp("GT:AD", "0/1:2,8") == 10
    assert qc._dyn_dp("GT:DP:AD", "0/1:.:2,8") == 10
    assert qc._dyn_dp("GT:AO:RO", "0/1:3:1") == 4


def test_dyn_dp_without_any_depth_field_is_none():
    assert qc._dyn_dp("GT", "0/1") is None
    assert qc._dyn_dp("GT:DP", "0/1:.") is None


def test_dyn_dp_of_a_non_numeric_dp_falls_through():
    assert qc._dyn_dp("GT:DP:AD", "0/1:deep:2,8") == 10
    assert qc._dyn_dp("GT:DP", "0/1:deep") is None


def test_dyn_dp_of_non_numeric_ad_or_ao_ro_is_none():
    assert qc._dyn_dp("GT:AD", "0/1:x,y") is None
    assert qc._dyn_dp("GT:AO:RO", "0/1:x:1") is None


# --------------------------------------------------------------------------- _dyn_ann / _dyn_opos / _dyn_num


ANN_INFO = ("DP=40;ANN=T|missense_variant|MODERATE|rpoB|GENE_rpoB|transcript|TX1|protein_coding|"
            "1/1|c.1349C>T|p.Ser450Leu|450/3519|;AF=0.9")


def test_dyn_ann_pulls_the_gene_effect_impact_and_protein_change():
    assert qc._dyn_ann(ANN_INFO) == ("rpoB", "missense_variant", "MODERATE", "p.Ser450Leu")


def test_dyn_ann_uses_only_the_first_of_several_annotations():
    info = "ANN=T|synonymous_variant|LOW|geneA|x|x|x|x|x|x|p.Aaa1Aaa|,T|missense_variant|MODERATE|geneB|"
    assert qc._dyn_ann(info)[0] == "geneA"


def test_dyn_ann_of_an_unannotated_or_truncated_record_is_all_empty_strings():
    assert qc._dyn_ann("DP=40;AF=0.9") == ("", "", "", "")
    assert qc._dyn_ann(".") == ("", "", "", "")
    assert qc._dyn_ann("ANN=T|missense_variant") == ("", "missense_variant", "", "")


def test_dyn_opos_reads_the_picard_liftover_original_coordinate():
    assert qc._dyn_opos("DP=40;OriginalContig=chrX;OriginalStart=1234") == "chrX:1234"


def test_dyn_opos_reads_a_plain_opos_field():
    assert qc._dyn_opos("DP=40;OPOS=chrX:1234") == "chrX:1234"


@pytest.mark.parametrize("info", [None, "", "DP=40", "OriginalContig=chrX", "OriginalStart=99", "OPOS="])
def test_dyn_opos_without_a_complete_original_coordinate_is_none(info):
    assert qc._dyn_opos(info) is None


@pytest.mark.parametrize(("value", "expected"), [
    ("P3", 3.0), ("passage 12", 12.0), ("-2.5", -2.5), (7, 7.0), ("7.", 7.0), ("2025-01-15", 2025.0),
])
def test_dyn_num_extracts_the_first_number_it_finds(value, expected):
    assert qc._dyn_num(value) == pytest.approx(expected)


@pytest.mark.parametrize("value", [None, "", "abc", "NA"])
def test_dyn_num_without_a_number_is_none(value):
    assert qc._dyn_num(value) is None


# --------------------------------------------------------------------------- build_dynamics


def variant(af, gene="rpoB", eff="missense_variant", imp="MODERATE", dp=40):
    return {"ref": "A", "alt": "T", "gene": gene, "eff": eff, "imp": imp,
            "aa": "p.Ser450Leu", "aa_h37rv": "", "af": af, "dp": dp}


def timeseries(afs_by_pos, group="G"):
    """(metadata, variants) for one connected group with one sample per timepoint."""
    n = len(next(iter(afs_by_pos.values())))
    metadata = {f"S{i + 1}": {"time": f"T{i + 1}", "tnum": float(i + 1), "group": group} for i in range(n)}
    variants = {f"S{i + 1}": {pos: variant(afs[i]) for pos, afs in afs_by_pos.items()} for i in range(n)}
    return metadata, variants


def test_build_dynamics_without_metadata_or_variants_is_none():
    metadata, variants = timeseries({"c:10": [0.05, 0.95]})
    assert qc.build_dynamics({}, variants) is None
    assert qc.build_dynamics(None, variants) is None
    assert qc.build_dynamics(metadata, {}) is None


def test_build_dynamics_builds_a_trajectory_with_its_event_flags():
    metadata, variants = timeseries({"c:10": [0.05, 0.95]})
    dyn = qc.build_dynamics(metadata, variants)
    assert dyn["thresholds"] == {"emerge": 0.25, "fix": 0.90, "loss": 0.10}
    assert len(dyn["groups"]) == 1
    group = dyn["groups"][0]
    assert group["group"] == "G"
    assert group["samples"] == ["S1", "S2"]
    assert group["times"] == ["T1", "T2"]
    assert group["tnums"] == [1.0, 2.0]
    assert group["n_flagged"] == 1
    series = group["series"][0]
    assert series["pos"] == "c:10"
    assert series["traj"] == [0.05, 0.95]
    assert series["dp"] == [40, 40]
    assert series["flags"] == ["emergence", "fixation", "nonsyn"]


def test_build_dynamics_flags_a_lost_variant():
    metadata, variants = timeseries({"c:10": [0.95, 0.05]})
    series = qc.build_dynamics(metadata, variants)["groups"][0]["series"][0]
    assert "loss" in series["flags"]
    assert "emergence" not in series["flags"]


def test_build_dynamics_flags_a_high_impact_variant():
    metadata, variants = timeseries({"c:10": [0.05, 0.95]})
    for sample in variants:
        variants[sample]["c:10"]["imp"] = "HIGH"
        variants[sample]["c:10"]["eff"] = "stop_gained"
    series = qc.build_dynamics(metadata, variants)["groups"][0]["series"][0]
    assert "high_impact" in series["flags"]


def test_build_dynamics_drops_a_variant_that_does_not_move():
    metadata, variants = timeseries({"c:10": [0.50, 0.55]})
    assert qc.build_dynamics(metadata, variants) is None


def test_build_dynamics_treats_an_uncalled_site_as_frequency_zero():
    metadata, variants = timeseries({"c:10": [0.05, 0.95]})
    del variants["S1"]["c:10"]
    series = qc.build_dynamics(metadata, variants)["groups"][0]["series"][0]
    assert series["traj"] == [0.0, 0.95]
    assert series["dp"] == [None, 40]


def test_build_dynamics_needs_more_than_one_distinct_timepoint():
    metadata, variants = timeseries({"c:10": [0.05, 0.95]})
    metadata["S2"]["time"] = "T1"
    assert qc.build_dynamics(metadata, variants) is None


def test_build_dynamics_ignores_a_sample_without_a_group():
    metadata, variants = timeseries({"c:10": [0.05, 0.95]})
    metadata["S1"]["group"] = None
    assert qc.build_dynamics(metadata, variants) is None


def test_build_dynamics_orders_samples_by_the_numeric_time():
    metadata, variants = timeseries({"c:10": [0.95, 0.05]})
    metadata["S1"]["tnum"], metadata["S2"]["tnum"] = 2.0, 1.0
    group = qc.build_dynamics(metadata, variants)["groups"][0]
    assert group["samples"] == ["S2", "S1"]
    assert group["series"][0]["traj"] == [0.05, 0.95]


def test_build_dynamics_attaches_the_group_invariant_sample_metadata():
    metadata, variants = timeseries({"c:10": [0.05, 0.95]})
    sample_meta = {"fields": ["patient", "arm"],
                   "rows": {"S1": {"patient": "P1", "arm": "A"}, "S2": {"patient": "P1", "arm": "B"}}}
    group = qc.build_dynamics(metadata, variants, sample_meta)["groups"][0]
    assert group["meta"] == {"patient": ["P1"], "arm": ["A", "B"]}


def test_build_dynamics_honours_custom_event_thresholds():
    metadata, variants = timeseries({"c:10": [0.30, 0.60]})
    # a 0.30 move is under a 0.5 min_move, so the variant is not reported at all
    assert qc.build_dynamics(metadata, variants, min_move=0.5) is None
    # at the default fixation cut-off 0.60 is not fixed; lowering it makes the same trajectory fix
    assert qc.build_dynamics(metadata, variants)["groups"][0]["series"][0]["flags"] == ["nonsyn"]
    series = qc.build_dynamics(metadata, variants, fix=0.55)["groups"][0]["series"][0]
    assert series["flags"] == ["fixation", "nonsyn"]


# --------------------------------------------------------------------------- build_epistasis


@pytest.fixture()
def four_point_dynamics():
    """One 4-timepoint series with two co-rising variants and one falling one."""
    rising = [0.05, 0.35, 0.65, 0.95]
    falling = list(reversed(rising))
    metadata, variants = timeseries({"c:10": rising, "c:20": rising, "c:30": falling})
    for sample in variants:
        variants[sample]["c:30"]["gene"] = "katG"
    return qc.build_dynamics(metadata, variants)


def test_build_epistasis_without_dynamics_is_none():
    assert qc.build_epistasis(None) is None
    assert qc.build_epistasis({}) is None
    assert qc.build_epistasis({"groups": []}) is None


def test_build_epistasis_needs_min_points_distinct_timepoints():
    metadata, variants = timeseries({"c:10": [0.05, 0.95], "c:20": [0.05, 0.95]})
    two_point = qc.build_dynamics(metadata, variants)
    assert two_point is not None
    assert qc.build_epistasis(two_point, perm=10) is None


def test_build_epistasis_reports_every_pair_over_the_correlation_threshold(four_point_dynamics):
    ep = qc.build_epistasis(four_point_dynamics, perm=100)
    assert ep["n_series"] == 1
    assert ep["n_candidates"] == 3
    assert len(ep["pairs"]) == 3
    assert {(p["posA"], p["posB"]) for p in ep["pairs"]} == {
        ("c:10", "c:20"), ("c:10", "c:30"), ("c:20", "c:30")}
    concordant = next(p for p in ep["pairs"] if p["posA"] == "c:10" and p["posB"] == "c:20")
    assert concordant["r"] == pytest.approx(1.0)
    assert concordant["direction"] == "concordant"
    assert concordant["consistent"] is True
    assert concordant["n"] == 1
    assert ep["n_concordant"] == 1
    assert ep["n_discordant"] == 2
    assert ep["pairs"][0]["geneA"] == "rpoB"


def test_build_epistasis_drops_pairs_below_min_r(four_point_dynamics):
    assert qc.build_epistasis(four_point_dynamics, min_r=1.5, perm=10) is None


def test_build_epistasis_q_values_are_bounded_and_monotone_in_p(four_point_dynamics):
    # the Benjamini-Hochberg loop is inline in build_epistasis, so it is only reachable this way
    ep = qc.build_epistasis(four_point_dynamics, perm=100)
    qs = [p["q"] for p in ep["pairs"]]
    assert all(0.0 < q <= 1.0 for q in qs)
    by_p = sorted(ep["pairs"], key=lambda p: p["p"])
    assert [p["q"] for p in by_p] == sorted(p["q"] for p in by_p)
    # the largest p keeps its own value (rank m: q = p * m / m)
    assert by_p[-1]["q"] == pytest.approx(by_p[-1]["p"], abs=1e-4)
    assert all(p["tier"] in ("strong", "moderate", "weak") for p in ep["pairs"])
    # the payload is sorted by q, then by descending |r|
    assert qs == sorted(qs)


def test_build_epistasis_is_reproducible_across_runs(four_point_dynamics):
    first = qc.build_epistasis(four_point_dynamics, perm=100)
    second = qc.build_epistasis(four_point_dynamics, perm=100)
    assert [(p["posA"], p["posB"], p["p"], p["q"]) for p in first["pairs"]] == \
           [(p["posA"], p["posB"], p["p"], p["q"]) for p in second["pairs"]]


def test_build_epistasis_ignores_a_series_whose_correlation_is_undefined():
    # a hand-made payload with a flat trajectory: _pearson returns None for every pair it
    # takes part in, so it contributes no pair and never reaches the correlation matrix
    dynamics = {"groups": [{"group": "G", "times": ["T1", "T2", "T3"], "series": [
        {"pos": "c:10", "gene": "a", "aa": "", "eff": "", "traj": [0.1, 0.5, 0.9]},
        {"pos": "c:20", "gene": "b", "aa": "", "eff": "", "traj": [0.5, 0.5, 0.5]},
        {"pos": "c:30", "gene": "c", "aa": "", "eff": "", "traj": [0.2, 0.6, 1.0]},
    ]}]}
    ep = qc.build_epistasis(dynamics, perm=20)
    assert [(p["posA"], p["posB"]) for p in ep["pairs"]] == [("c:10", "c:30")]
    assert ep["matrix"]["total_nodes"] == 2


def test_build_epistasis_matrix_holds_every_pairwise_mean_r(four_point_dynamics):
    matrix = qc.build_epistasis(four_point_dynamics, perm=10)["matrix"]
    assert matrix["total_nodes"] == 3
    assert matrix["truncated"] is False
    assert len(matrix["nodes"]) == 3
    assert sorted(matrix["cells"].values()) == [-1.0, -1.0, 1.0]


def test_build_epistasis_reports_its_scaling_caps(four_point_dynamics):
    ep = qc.build_epistasis(four_point_dynamics, perm=10)
    assert ep["vars_capped"] is False
    assert ep["pairs_capped"] is False
    assert ep["n_vars_max"] == 3
    assert ep["max_vars"] == 250
    assert ep["max_pairs_perm"] == 800


def test_build_epistasis_caps_the_variants_tested_per_series(four_point_dynamics):
    capped = qc.build_epistasis(four_point_dynamics, perm=10, max_vars=2)
    assert capped["vars_capped"] is True
    assert capped["n_vars_max"] == 3
    assert len(capped["pairs"]) == 1   # only one pair survives out of the two kept variants


def test_build_epistasis_caps_the_pairs_that_get_a_permutation_test(four_point_dynamics):
    capped = qc.build_epistasis(four_point_dynamics, perm=10, max_pairs_perm=1)
    assert capped["pairs_capped"] is True
    assert capped["n_candidates"] == 3   # all three were candidates, only the strongest was tested
    assert len(capped["pairs"]) == 1


# --------------------------------------------------------------------------- build_snp_matrix


def test_build_snp_matrix_without_variants_is_none():
    assert qc.build_snp_matrix({}) is None
    assert qc.build_snp_matrix(None) is None


def test_build_snp_matrix_of_a_sample_without_calls_is_an_empty_matrix():
    assert qc.build_snp_matrix({"S1": {}}) == {
        "samples": ["S1"], "reference": "", "rows": [], "total_sites": 0, "truncated": False}


def test_build_snp_matrix_places_one_row_per_site_with_the_called_cells_only():
    variants = {"S1": {"c:10": variant(0.9), "c:20": variant(0.5, gene="", eff="")},
                "S2": {"c:10": variant(0.4, dp=12)}}
    matrix = qc.build_snp_matrix(variants, reference="TESTREF")
    assert matrix["samples"] == ["S1", "S2"]
    assert matrix["reference"] == "TESTREF"
    assert matrix["total_sites"] == 2
    assert matrix["truncated"] is False
    shared, private = matrix["rows"]
    assert (shared["contig"], shared["pos"], shared["ref"], shared["alt"]) == ("c", 10, "A", "T")
    assert shared["gene"] == "rpoB"
    assert shared["n"] == 2
    assert shared["cells"] == {0: [0.9, 40], 1: [0.4, 12]}
    assert private["pos"] == 20
    assert private["cells"] == {0: [0.5, 40]}


def test_build_snp_matrix_skips_a_site_with_a_non_numeric_position():
    matrix = qc.build_snp_matrix({"S1": {"c:10": variant(0.9), "c:bad": variant(0.9)}})
    assert [r["pos"] for r in matrix["rows"]] == [10]


def test_build_snp_matrix_merges_alternate_alleles_of_the_same_site():
    a, b = variant(0.9), variant(0.9)
    b["alt"] = "G"
    matrix = qc.build_snp_matrix({"S1": {"c:10": a}, "S2": {"c:10": b}})
    assert matrix["rows"][0]["alt"] == "G,T"


def test_build_snp_matrix_truncation_keeps_the_most_shared_sites_in_coordinate_order():
    variants = {"S1": {"c:10": variant(0.9), "c:20": variant(0.9), "c:30": variant(0.9)},
                "S2": {"c:10": variant(0.9), "c:20": variant(0.9)}}
    matrix = qc.build_snp_matrix(variants, max_sites=2)
    assert matrix["truncated"] is True
    assert matrix["total_sites"] == 3
    assert [r["pos"] for r in matrix["rows"]] == [10, 20]


# --------------------------------------------------------------------------- build_gene_map


def test_build_gene_map_without_a_gff_is_the_curated_h37rv_map():
    assert qc.build_gene_map(None) == qc.GENE_RV
    assert qc.build_gene_map("/nonexistent/reference.gff3") == qc.GENE_RV


def test_build_gene_map_ignores_a_non_h37rv_reference(data_dir):
    # the committed fixture GFF uses TEST_000n locus tags, which are not H37Rv-style
    assert qc.build_gene_map(str(data_dir / "reference.gff3")) == qc.GENE_RV


def test_build_gene_map_overlays_h37rv_style_locus_tags_from_the_run_gff(tmp_path):
    gff = tmp_path / "ref.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "chr\ttest\tgene\t1\t100\t.\t+\t.\tID=g1;gene=rpoB;locus_tag=Rv9999\n"
        "chr\ttest\tgene\t200\t300\t.\t+\t.\tID=g2;gene=newGene;locus_tag=Rv1234\n"
        "chr\ttest\tgene\t400\t500\t.\t+\t.\tID=g3;gene=strainGene;locus_tag=ABC_0001\n",
        encoding="utf-8")
    gene_map = qc.build_gene_map(str(gff))
    assert gene_map["rpoB"] == "Rv9999"          # overrides the curated Rv0667
    assert gene_map["newGene"] == "Rv1234"       # a gene the curated map does not carry
    assert "strainGene" not in gene_map          # a strain-specific tag is ignored
    assert gene_map["katG"] == qc.GENE_RV["katG"]
