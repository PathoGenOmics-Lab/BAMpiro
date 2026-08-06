"""Unit tests for bin/collect_dr.py.

Aggregates the per-sample pathotypr `*_mutations.tsv` drug-resistance detail
files into one run TSV. Only `sample_from_name` is importable; the row parsing
lives inline in `main()`, so it is exercised through the CLI.
"""

from __future__ import annotations

import subprocess
import sys

import pytest

from conftest import load_script

cdr = load_script("collect_dr")

_OUT_HEADER = ["sample", "drug", "gene", "mutation", "grade", "marker_name", "af", "dp"]

# pathotypr split-fastq column layout for the WHO DR marker panel
_IN_HEADER = "pos\tref_allele\talt_allele\tref_count\talt_count\talt_fraction\tlineage_path\n"


# --------------------------------------------------------------------------
# sample_from_name
# --------------------------------------------------------------------------

@pytest.mark.parametrize("path,expected", [
    # the three suffix forms the regex strips
    ("S1.dr_mutations.tsv", "S1"),
    ("S1_mutations.tsv", "S1"),
    ("S1.tsv", "S1"),
    # basename is taken first
    ("/some/deep/path/S1.dr_mutations.tsv", "S1"),
    ("./rel/S1_mutations.tsv", "S1"),
    # '__' is the pipeline's sample/run separator: everything after it is dropped
    ("S1__run2.dr_mutations.tsv", "S1"),
    ("S1__run2_mutations.tsv", "S1"),
    ("S1__run2.tsv", "S1"),
    ("S1__run2__ref3.dr_mutations.tsv", "S1"),
    ("S1__.tsv", "S1"),
    # no recognised suffix: the name is returned as-is (minus any '__' tail)
    ("S1", "S1"),
    ("S1__run2", "S1"),
    # only the trailing occurrence is stripped, and only one of the three
    ("sample.with.dots.tsv", "sample.with.dots"),
    ("S1.dr_mutations.tsv.gz", "S1.dr_mutations.tsv.gz"),
    ("S1_mutations.tsv.bak", "S1_mutations.tsv.bak"),
    ("", ""),
])
def test_sample_from_name(path, expected):
    assert cdr.sample_from_name(path) == expected


def test_sample_from_name_leading_double_underscore_yields_empty_sample():
    """PIN: a name that starts with '__' splits into an empty sample id."""
    assert cdr.sample_from_name("__run2.dr_mutations.tsv") == ""


# --------------------------------------------------------------------------
# end to end (the row parsing lives inside main())
# --------------------------------------------------------------------------

def _run(tmp_path, repo_root, files, output="dr.tsv", sample=None):
    cmd = [sys.executable, str(repo_root / "bin" / "collect_dr.py"), "-o", str(tmp_path / output)]
    if sample is not None:
        cmd += ["--sample", sample]
    cmd += [str(f) for f in files]
    proc = subprocess.run(cmd, cwd=tmp_path, capture_output=True, text=True, timeout=120)
    assert proc.returncode == 0, proc.stderr
    rows = [line.split("\t") for line in (tmp_path / output).read_text().splitlines()]
    return proc, rows


def test_end_to_end_parses_a_mutations_file(tmp_path, repo_root):
    f = tmp_path / "S1.dr_mutations.tsv"
    f.write_text(
        _IN_HEADER
        + "761155\tC\tT\t2\t38\t95.0\tRIF;R;marker1;A;rpoB;S450L\n"
        + "2155168\tG\tC\t5\t45\t90.0\tINH;R;marker2;B;katG;S315T\n"
    )
    proc, rows = _run(tmp_path, repo_root, [f])

    assert rows[0] == _OUT_HEADER
    # sample, drug, gene, mutation, grade, marker_name, af, dp
    assert rows[1] == ["S1", "RIF", "rpoB", "S450L", "A", "marker1", "0.9500", "40"]
    assert rows[2] == ["S1", "INH", "katG", "S315T", "B", "marker2", "0.9000", "50"]
    assert "2 resistance call(s) from 1 file(s)" in proc.stderr


def test_end_to_end_alt_fraction_is_a_percent(tmp_path, repo_root):
    """alt_fraction comes out of pathotypr as a PERCENT and is divided by 100,
    then formatted to four decimals."""
    f = tmp_path / "S1.dr_mutations.tsv"
    f.write_text(
        _IN_HEADER
        + "1\tA\tG\t0\t100\t100.0\tRIF;R;m;A;rpoB;X\n"
        + "2\tA\tG\t95\t5\t5.0\tRIF;R;m;A;rpoB;Y\n"
        + "3\tA\tG\t99\t1\t1.234\tRIF;R;m;A;rpoB;Z\n"
    )
    _, rows = _run(tmp_path, repo_root, [f])
    assert [r[6] for r in rows[1:]] == ["1.0000", "0.0500", "0.0123"]


def test_end_to_end_dp_is_ref_plus_alt(tmp_path, repo_root):
    f = tmp_path / "S1.dr_mutations.tsv"
    f.write_text(
        _IN_HEADER
        + "1\tA\tG\t7\t13\t65.0\tRIF;R;m;A;rpoB;X\n"
        + "2\tA\tG\t0\t0\t0.0\tRIF;R;m;A;rpoB;Y\n"
        + "3\tA\tG\t\t9\t90.0\tRIF;R;m;A;rpoB;Z\n"    # empty ref_count -> 0
    )
    _, rows = _run(tmp_path, repo_root, [f])
    assert [r[7] for r in rows[1:]] == ["20", "0", "9"]


def test_end_to_end_lineage_path_maxsplit_keeps_semicolons_in_the_mutation(tmp_path, repo_root):
    """split(';', 5) means only the first five separators are consumed, so a
    ';' inside the trailing mutation descriptor survives intact."""
    f = tmp_path / "S1.dr_mutations.tsv"
    f.write_text(
        _IN_HEADER
        + "1\tA\tG\t1\t9\t90.0\tRIF;R;marker1;A;rpoB;S450L;c.1349C>T\n"
    )
    _, rows = _run(tmp_path, repo_root, [f])
    assert rows[1][3] == "S450L;c.1349C>T"
    assert rows[1][:3] == ["S1", "RIF", "rpoB"]
    assert rows[1][4:6] == ["A", "marker1"]


def test_end_to_end_short_lineage_path_leaves_trailing_fields_empty(tmp_path, repo_root):
    f = tmp_path / "S1.dr_mutations.tsv"
    f.write_text(
        _IN_HEADER
        + "1\tA\tG\t1\t9\t90.0\tRIF;R;marker1\n"
        + "2\tA\tG\t1\t9\t90.0\t\n"
    )
    _, rows = _run(tmp_path, repo_root, [f])
    assert rows[1] == ["S1", "RIF", "", "", "", "marker1", "0.9000", "10"]
    assert rows[2] == ["S1", "", "", "", "", "", "0.9000", "10"]


@pytest.mark.parametrize("alt_fraction", ["NA", "", "not_a_number", "95%"])
def test_end_to_end_non_numeric_alt_fraction_degrades_to_empty(tmp_path, repo_root, alt_fraction):
    """A bad value produces an empty cell rather than raising."""
    f = tmp_path / "S1.dr_mutations.tsv"
    f.write_text(_IN_HEADER + f"1\tA\tG\t1\t9\t{alt_fraction}\tRIF;R;m;A;rpoB;X\n")
    _, rows = _run(tmp_path, repo_root, [f])
    assert rows[1][6] == ""
    assert rows[1][7] == "10"     # dp is unaffected


@pytest.mark.parametrize("counts", [("x", "9"), ("1", "y"), ("1.5", "9")])
def test_end_to_end_non_numeric_counts_degrade_to_empty(tmp_path, repo_root, counts):
    ref_count, alt_count = counts
    f = tmp_path / "S1.dr_mutations.tsv"
    f.write_text(_IN_HEADER + f"1\tA\tG\t{ref_count}\t{alt_count}\t90.0\tRIF;R;m;A;rpoB;X\n")
    _, rows = _run(tmp_path, repo_root, [f])
    assert rows[1][7] == ""
    assert rows[1][6] == "0.9000"  # af is unaffected


def test_end_to_end_missing_columns_degrade_to_empty(tmp_path, repo_root):
    """A header without alt_fraction / counts still yields a row, just blank."""
    f = tmp_path / "S1.dr_mutations.tsv"
    f.write_text("pos\tlineage_path\n1\tRIF;R;marker1;A;rpoB;S450L\n")
    _, rows = _run(tmp_path, repo_root, [f])
    assert rows[1] == ["S1", "RIF", "rpoB", "S450L", "A", "marker1", "", "0"]


def test_end_to_end_header_is_case_insensitive_and_blank_lines_are_skipped(tmp_path, repo_root):
    f = tmp_path / "S1.dr_mutations.tsv"
    f.write_text(
        "\nPOS\tREF_ALLELE\tALT_ALLELE\tREF_COUNT\tALT_COUNT\t ALT_FRACTION \tLineage_Path\n"
        "\n1\tA\tG\t1\t9\t90.0\tRIF;R;m;A;rpoB;X\n\n"
    )
    _, rows = _run(tmp_path, repo_root, [f])
    assert len(rows) == 2
    assert rows[1] == ["S1", "RIF", "rpoB", "X", "A", "m", "0.9000", "10"]


def test_end_to_end_infers_the_sample_per_file(tmp_path, repo_root):
    a = tmp_path / "sampleA.dr_mutations.tsv"
    b = tmp_path / "sampleB__run7.dr_mutations.tsv"
    a.write_text(_IN_HEADER + "1\tA\tG\t1\t9\t90.0\tRIF;R;m;A;rpoB;X\n")
    b.write_text(_IN_HEADER + "2\tA\tG\t1\t9\t90.0\tINH;R;m;A;katG;Y\n")
    _, rows = _run(tmp_path, repo_root, [a, b])
    assert [r[0] for r in rows[1:]] == ["sampleA", "sampleB"]


def test_end_to_end_sample_flag_overrides_every_file(tmp_path, repo_root):
    a = tmp_path / "sampleA.dr_mutations.tsv"
    b = tmp_path / "sampleB.dr_mutations.tsv"
    a.write_text(_IN_HEADER + "1\tA\tG\t1\t9\t90.0\tRIF;R;m;A;rpoB;X\n")
    b.write_text(_IN_HEADER + "2\tA\tG\t1\t9\t90.0\tINH;R;m;A;katG;Y\n")
    _, rows = _run(tmp_path, repo_root, [a, b], sample="OVERRIDE")
    assert [r[0] for r in rows[1:]] == ["OVERRIDE", "OVERRIDE"]


def test_end_to_end_skips_non_existent_inputs_silently(tmp_path, repo_root):
    good = tmp_path / "S1.dr_mutations.tsv"
    good.write_text(_IN_HEADER + "1\tA\tG\t1\t9\t90.0\tRIF;R;m;A;rpoB;X\n")
    proc, rows = _run(tmp_path, repo_root, [tmp_path / "absent.dr_mutations.tsv", good])
    assert len(rows) == 2
    assert "absent" not in proc.stderr
    # the summary line counts the files ASKED for, not the ones actually read
    assert "1 resistance call(s) from 2 file(s)" in proc.stderr


def test_end_to_end_writes_the_header_with_zero_data_rows(tmp_path, repo_root):
    """The output file is always created, so the downstream report never has a
    missing input, even when nothing resistant was detected."""
    empty = tmp_path / "S1.dr_mutations.tsv"
    empty.write_text(_IN_HEADER)               # header only
    _, rows = _run(tmp_path, repo_root, [empty])
    assert rows == [_OUT_HEADER]


def test_end_to_end_with_no_input_files_at_all(tmp_path, repo_root):
    proc, rows = _run(tmp_path, repo_root, [])
    assert rows == [_OUT_HEADER]
    assert "0 resistance call(s) from 0 file(s)" in proc.stderr


def test_end_to_end_completely_empty_file(tmp_path, repo_root):
    f = tmp_path / "S1.dr_mutations.tsv"
    f.write_text("")
    _, rows = _run(tmp_path, repo_root, [f])
    assert rows == [_OUT_HEADER]


def test_end_to_end_unreadable_input_warns_but_does_not_fail(tmp_path, repo_root):
    """A directory passes the exists() check and then fails to open; the OSError
    is reported on stderr and the remaining inputs are still processed."""
    d = tmp_path / "S9.dr_mutations.tsv"
    d.mkdir()
    good = tmp_path / "S1.dr_mutations.tsv"
    good.write_text(_IN_HEADER + "1\tA\tG\t1\t9\t90.0\tRIF;R;m;A;rpoB;X\n")
    proc, rows = _run(tmp_path, repo_root, [d, good])
    assert "[collect_dr] WARN could not read" in proc.stderr
    assert len(rows) == 2
    assert rows[1][0] == "S1"
