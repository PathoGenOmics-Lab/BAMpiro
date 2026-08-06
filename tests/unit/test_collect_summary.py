"""Unit tests for bin/collect_summary.py.

Phase-B aggregation: turns the per-sample legacy `.log` files written by
stats_to_legacy.py into one wide summary TSV, plus an optional cohort
gene-burden table.
"""

from __future__ import annotations

import subprocess
import sys

import pytest

from conftest import load_script

cs = load_script("collect_summary")


# --------------------------------------------------------------------------
# _to_int
# --------------------------------------------------------------------------

@pytest.mark.parametrize("raw,expected", [
    ("5", 5),
    ("0", 0),
    ("-3", -3),
    ("  7  ", 7),      # int() tolerates surrounding whitespace
    (12, 12),
    (5.9, 5),          # truncates towards zero
])
def test_to_int_accepts(raw, expected):
    assert cs._to_int(raw) == expected


@pytest.mark.parametrize("raw", [None, "", "NA", "abc", "5.5", "1e3", []])
def test_to_int_rejects(raw):
    """Anything unparsable becomes None (never raises, never 0)."""
    assert cs._to_int(raw) is None


# --------------------------------------------------------------------------
# parse_log
# --------------------------------------------------------------------------

def test_parse_log_reads_tab_separated_pairs(tmp_path):
    p = tmp_path / "S1.log"
    p.write_text("GID_IN\tS1\nSNP_OUT\t42\nTITV_OUT\tNA\n")
    assert cs.parse_log(str(p)) == {"GID_IN": "S1", "SNP_OUT": "42", "TITV_OUT": "NA"}


def test_parse_log_splits_only_on_the_first_tab(tmp_path):
    """split('\\t', 1) means a value is allowed to contain tabs of its own."""
    p = tmp_path / "S1.log"
    p.write_text("KEY\tvalue\twith\ttabs\n")
    assert cs.parse_log(str(p)) == {"KEY": "value\twith\ttabs"}


def test_parse_log_skips_blank_and_tabless_lines(tmp_path):
    p = tmp_path / "S1.log"
    p.write_text("GID_IN\tS1\n\n   \nJUST_A_KEY\nSNP_OUT\t7\n")
    assert cs.parse_log(str(p)) == {"GID_IN": "S1", "SNP_OUT": "7"}


def test_parse_log_last_duplicate_key_wins(tmp_path):
    p = tmp_path / "S1.log"
    p.write_text("GID_IN\tfirst\nGID_IN\tsecond\n")
    assert cs.parse_log(str(p)) == {"GID_IN": "second"}


def test_parse_log_missing_file_returns_empty_dict(tmp_path, capsys):
    """Never raises: an unreadable log degrades to {} and a warning."""
    assert cs.parse_log(str(tmp_path / "absent.log")) == {}
    assert "Could not parse" in capsys.readouterr().err


def test_parse_log_empty_file_returns_empty_dict(tmp_path):
    p = tmp_path / "empty.log"
    p.write_text("")
    assert cs.parse_log(str(p)) == {}


# --------------------------------------------------------------------------
# _derive_ann
# --------------------------------------------------------------------------

def _ann_data(**overrides):
    data = {
        "ANN_PRESENT": "yes",
        "ANN_TOTAL": "100",
        "ANN_ANNOTATED": "80",
        "ANN_LOF": "5",
        "ANN_CODING": "50",
        "ANN_IMPACT_HIGH": "5",
        "ANN_IMPACT_MODERATE": "20",
        "ANN_IMPACT_LOW": "25",
        "ANN_IMPACT_MODIFIER": "30",
    }
    data.update(overrides)
    return data


@pytest.mark.parametrize("present", ["no", "NA", "", "YES", None])
def test_derive_ann_requires_ann_present_yes(present):
    data = _ann_data()
    if present is None:
        data.pop("ANN_PRESENT")
    else:
        data["ANN_PRESENT"] = present
    assert cs._derive_ann(data) == ("NA", "NA", "NA")


def test_derive_ann_returns_formatted_strings():
    ann_pct, coding_pct, lof_pct = cs._derive_ann(_ann_data())
    assert (ann_pct, coding_pct, lof_pct) == ("80.00", "62.50", "10.00")
    assert all(isinstance(x, str) for x in (ann_pct, coding_pct, lof_pct))


@pytest.mark.parametrize("total", ["0", "NA", "", None])
def test_derive_ann_guard_one_zero_or_missing_total(total):
    """Guard 1: annotated_pct needs a non-zero ANN_TOTAL denominator."""
    data = _ann_data()
    if total is None:
        data.pop("ANN_TOTAL")
    else:
        data["ANN_TOTAL"] = total
    ann_pct, coding_pct, lof_pct = cs._derive_ann(data)
    assert ann_pct == "NA"
    # the other two are unaffected by this guard
    assert (coding_pct, lof_pct) == ("62.50", "10.00")


def test_derive_ann_guard_one_also_trips_on_an_unparsable_numerator():
    ann_pct, _, _ = cs._derive_ann(_ann_data(ANN_ANNOTATED="NA"))
    assert ann_pct == "NA"


def test_derive_ann_guard_two_zero_impact_denominator():
    """Guard 2: coding_pct needs the four impact buckets to sum above zero."""
    data = _ann_data(ANN_IMPACT_HIGH="0", ANN_IMPACT_MODERATE="0",
                     ANN_IMPACT_LOW="0", ANN_IMPACT_MODIFIER="0")
    ann_pct, coding_pct, lof_pct = cs._derive_ann(data)
    assert coding_pct == "NA"
    assert (ann_pct, lof_pct) == ("80.00", "10.00")


def test_derive_ann_guard_two_all_impacts_unparsable():
    """sum() over an all-None generator is 0, so the guard still catches it."""
    data = _ann_data(ANN_IMPACT_HIGH="NA", ANN_IMPACT_MODERATE="NA",
                     ANN_IMPACT_LOW="NA", ANN_IMPACT_MODIFIER="NA")
    assert cs._derive_ann(data)[1] == "NA"


def test_derive_ann_guard_two_treats_unparsable_impacts_as_zero_in_the_numerator():
    data = _ann_data(ANN_IMPACT_HIGH="NA", ANN_IMPACT_MODERATE="10",
                     ANN_IMPACT_LOW="10", ANN_IMPACT_MODIFIER="80")
    # numerator = 0 + 10 + 10, denominator = 100
    assert cs._derive_ann(data)[1] == "20.00"


@pytest.mark.parametrize("coding", ["0", "NA", "", None])
def test_derive_ann_guard_three_zero_or_missing_coding(coding):
    """Guard 3: lof_pct needs a non-zero ANN_CODING denominator."""
    data = _ann_data()
    if coding is None:
        data.pop("ANN_CODING")
    else:
        data["ANN_CODING"] = coding
    ann_pct, coding_pct, lof_pct = cs._derive_ann(data)
    assert lof_pct == "NA"
    assert (ann_pct, coding_pct) == ("80.00", "62.50")


def test_derive_ann_guard_three_also_trips_on_an_unparsable_numerator():
    assert cs._derive_ann(_ann_data(ANN_LOF="NA"))[2] == "NA"


def test_derive_ann_lof_pct_stays_within_range():
    """CODING is the superset of LOF by construction, so the ratio is bounded."""
    assert cs._derive_ann(_ann_data(ANN_LOF="50", ANN_CODING="50"))[2] == "100.00"


# --------------------------------------------------------------------------
# _accumulate_gene_burden
# --------------------------------------------------------------------------

@pytest.mark.parametrize("burden", ["NA", "", None])
def test_accumulate_gene_burden_is_a_noop_for_empty_input(burden):
    agg = {}
    cs._accumulate_gene_burden(agg, burden)
    assert agg == {}


def test_accumulate_gene_burden_mutates_the_aggregate_in_place():
    agg = {}
    assert cs._accumulate_gene_burden(agg, "rpoB:2:1:missense") is None
    assert agg == {"rpoB": {"high": 2, "mod": 1, "n_samples": 1, "eff": {"missense": 1}}}


def test_accumulate_gene_burden_sums_across_samples():
    agg = {}
    cs._accumulate_gene_burden(agg, "rpoB:2:1:missense;katG:1:0:frameshift")
    cs._accumulate_gene_burden(agg, "rpoB:3:4:stop_gained")
    assert agg["rpoB"] == {"high": 5, "mod": 5, "n_samples": 2,
                           "eff": {"missense": 1, "stop_gained": 1}}
    assert agg["katG"] == {"high": 1, "mod": 0, "n_samples": 1, "eff": {"frameshift": 1}}


@pytest.mark.parametrize("token", [
    "rpoB:2:1",             # too few fields
    "rpoB:2:1:missense:x",  # too many fields
    "rpoB",
    "",
    "rpoB:x:1:missense",    # non-integer HIGH
    "rpoB:2:y:missense",    # non-integer MOD
])
def test_accumulate_gene_burden_skips_malformed_tokens(token):
    agg = {}
    cs._accumulate_gene_burden(agg, token)
    assert agg == {}


def test_accumulate_gene_burden_skips_only_the_bad_token():
    agg = {}
    cs._accumulate_gene_burden(agg, "bad:token;katG:1:0:frameshift;also:bad")
    assert list(agg) == ["katG"]


def test_accumulate_gene_burden_counts_n_samples_once_per_token():
    """PIN: n_samples is incremented per TOKEN, not per distinct sample. A
    per-sample burden string never repeats a gene (stats_to_legacy builds it
    from a dict), but if it did, that one sample would count twice."""
    agg = {}
    cs._accumulate_gene_burden(agg, "rpoB:1:0:missense;rpoB:2:0:missense")
    assert agg["rpoB"]["n_samples"] == 2
    assert agg["rpoB"]["high"] == 3


def test_accumulate_gene_burden_preserves_existing_entries():
    agg = {"katG": {"high": 9, "mod": 9, "n_samples": 4, "eff": {"missense": 4}}}
    cs._accumulate_gene_burden(agg, "katG:1:1:frameshift")
    assert agg["katG"] == {"high": 10, "mod": 10, "n_samples": 5,
                           "eff": {"missense": 4, "frameshift": 1}}


# --------------------------------------------------------------------------
# _write_gene_burden
# --------------------------------------------------------------------------

def test_write_gene_burden_header_and_sort_order(tmp_path):
    agg = {
        "geneLow": {"high": 0, "mod": 1, "n_samples": 1, "eff": {"missense": 1}},
        "geneTop": {"high": 5, "mod": 0, "n_samples": 3, "eff": {"stop_gained": 2, "frameshift": 1}},
        "zeta": {"high": 2, "mod": 2, "n_samples": 1, "eff": {"missense": 1}},
        "alpha": {"high": 2, "mod": 2, "n_samples": 2, "eff": {"missense": 2}},
        "geneMid": {"high": 2, "mod": 7, "n_samples": 1, "eff": {"splice": 1}},
    }
    out = tmp_path / "gene_burden.tsv"
    cs._write_gene_burden(agg, str(out))
    rows = [line.split("\t") for line in out.read_text().splitlines()]

    assert rows[0] == ["gene", "high", "moderate", "low", "modifier",
                       "dominant_effect", "n_samples", "total_impactful", "start", "end"]
    # sorted by -high, then -moderate, then gene name
    assert [r[0] for r in rows[1:]] == ["geneTop", "geneMid", "alpha", "zeta", "geneLow"]
    assert rows[1] == ["geneTop", "5", "0", "", "", "stop_gained", "3", "5", "", ""]
    assert rows[2] == ["geneMid", "2", "7", "", "", "splice", "1", "9", "", ""]


def test_write_gene_burden_empty_effect_map_gives_na(tmp_path):
    out = tmp_path / "gene_burden.tsv"
    cs._write_gene_burden({"geneX": {"high": 1, "mod": 0, "n_samples": 1, "eff": {}}}, str(out))
    assert out.read_text().splitlines()[1].split("\t")[5] == "NA"


def test_write_gene_burden_empty_aggregate_writes_header_only(tmp_path):
    out = tmp_path / "gene_burden.tsv"
    cs._write_gene_burden({}, str(out))
    assert len(out.read_text().splitlines()) == 1


# --------------------------------------------------------------------------
# end to end
# --------------------------------------------------------------------------

def _log_text(sample, **overrides):
    data = {
        "GID_IN": sample,
        "TRD_OUT": "1000",
        "TRM_OUT": "900",
        "SNP_OUT": "42",
        "ANN_PRESENT": "yes",
        "ANN_TOTAL": "100",
        "ANN_ANNOTATED": "80",
        "ANN_LOF": "5",
        "ANN_CODING": "50",
        "ANN_IMPACT_HIGH": "5",
        "ANN_IMPACT_MODERATE": "20",
        "ANN_IMPACT_LOW": "25",
        "ANN_IMPACT_MODIFIER": "30",
        "ANN_GENE_BURDEN": "rpoB:2:1:missense",
    }
    data.update(overrides)
    return "".join(f"{k}\t{v}\n" for k, v in data.items())


def test_end_to_end_summary(tmp_path, repo_root):
    """Row order follows the SORTED FILE NAMES, not the sample ids inside."""
    (tmp_path / "a.log").write_text(_log_text("zeta"))
    (tmp_path / "b.log").write_text(_log_text("alpha", ANN_PRESENT="no",
                                              ANN_GENE_BURDEN="katG:3:0:frameshift"))
    out = tmp_path / "summary.tsv"
    burden = tmp_path / "gene_burden.tsv"

    proc = subprocess.run(
        [
            sys.executable, str(repo_root / "bin" / "collect_summary.py"),
            # deliberately passed out of order, and with one path that does not exist
            str(tmp_path / "b.log"), str(tmp_path / "a.log"), str(tmp_path / "absent.log"),
            "-o", str(out),
            "--gene-burden-out", str(burden),
        ],
        cwd=tmp_path, capture_output=True, text=True, timeout=120,
    )
    assert proc.returncode == 0, proc.stderr
    assert "2 samples" in proc.stderr

    rows = [line.split("\t") for line in out.read_text().splitlines()]
    header = rows[0]
    assert header == [cs.SUMMARY_COLUMNS[k] for k in cs.COLUMN_ORDER] + cs.DERIVED_COLUMNS
    assert header[0] == "sample_id"
    assert header[-3:] == ["annotated_pct", "coding_pct", "lof_pct"]

    assert len(rows) == 3
    assert [r[0] for r in rows[1:]] == ["zeta", "alpha"]   # a.log before b.log
    assert all(len(r) == len(header) for r in rows)

    first = dict(zip(header, rows[1]))
    assert first["raw_reads"] == "1000"
    assert first["snps"] == "42"
    assert first["annotated_pct"] == "80.00"
    assert first["coding_pct"] == "62.50"
    assert first["lof_pct"] == "10.00"
    # A field absent from the .log becomes 'NA' in the summary.
    assert first["median_depth"] == "NA"

    second = dict(zip(header, rows[2]))
    assert second["ann_present"] == "no"
    assert (second["annotated_pct"], second["coding_pct"], second["lof_pct"]) == ("NA", "NA", "NA")

    # The gene burden is accumulated even for the sample whose ANN_PRESENT is 'no'.
    burden_rows = [line.split("\t") for line in burden.read_text().splitlines()]
    assert [r[0] for r in burden_rows[1:]] == ["katG", "rpoB"]


def test_end_to_end_with_no_usable_logs_writes_header_only(tmp_path, repo_root):
    (tmp_path / "empty.log").write_text("")
    out = tmp_path / "summary.tsv"
    burden = tmp_path / "gene_burden.tsv"
    proc = subprocess.run(
        [
            sys.executable, str(repo_root / "bin" / "collect_summary.py"),
            str(tmp_path / "empty.log"), str(tmp_path / "absent.log"),
            "-o", str(out), "--gene-burden-out", str(burden),
        ],
        cwd=tmp_path, capture_output=True, text=True, timeout=120,
    )
    assert proc.returncode == 0, proc.stderr
    assert out.read_text().splitlines() == ["\t".join(
        [cs.SUMMARY_COLUMNS[k] for k in cs.COLUMN_ORDER] + cs.DERIVED_COLUMNS)]
    # The aux table is only written when at least one gene was aggregated.
    assert not burden.exists()
