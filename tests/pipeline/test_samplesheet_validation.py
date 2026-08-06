"""The samplesheet validator in main.nf.

Validation runs at script level, before any process is scheduled, so these tests are fast:
Nextflow fails during script evaluation and never touches an executor. Its defining feature is
that it collects EVERY problem and reports them together, so the tests assert on the full list
rather than on the first failure.
"""

from __future__ import annotations

import pytest

from conftest import DATA, run_pipeline

pytestmark = pytest.mark.nextflow

INVALID = DATA / "invalid"


def validate(nextflow, nextflow_env, tmp_path, sheet):
    """Run the pipeline far enough to trigger samplesheet validation, and return the result."""
    return run_pipeline(
        nextflow,
        nextflow_env,
        tmp_path,
        "-profile", "test",
        "-stub-run",
        "--tsv", str(sheet),
        "--outdir", str(tmp_path / "out"),
        timeout=300,
    )


def test_the_fixture_samplesheet_is_accepted(nextflow, nextflow_env, tmp_path):
    result = validate(nextflow, nextflow_env, tmp_path, DATA / "samplesheet.tsv")
    assert result.returncode == 0, result.stdout + result.stderr


@pytest.mark.parametrize(
    ("sheet", "expected"),
    [
        ("missing_column.tsv", "missing the required column 'refGff'"),
        ("no_rows.tsv", "has a header but no data rows"),
        ("no_usable_rows.tsv", "produced no usable samples"),
        ("empty_field.tsv", "empty required field(s): sampleId"),
        ("missing_file.tsv", "file(s) not found"),
        ("ref_conflict.tsv", "points to a different refFasta"),
    ],
)
def test_each_kind_of_broken_samplesheet_is_rejected(
    nextflow, nextflow_env, tmp_path, sheet, expected
):
    result = validate(nextflow, nextflow_env, tmp_path, INVALID / sheet)
    output = result.stdout + result.stderr
    assert result.returncode != 0, f"{sheet} was accepted:\n{output}"
    assert expected in output, f"{sheet} did not explain the problem:\n{output}"


def test_every_problem_is_reported_in_one_pass(nextflow, nextflow_env, tmp_path):
    """The whole point of the validator: fix everything in one edit, not one error per rerun."""
    result = validate(nextflow, nextflow_env, tmp_path, INVALID / "multi_error.tsv")
    output = result.stdout + result.stderr

    assert result.returncode != 0
    assert "has 3 problem(s)" in output, f"expected all three problems at once:\n{output}"
    for expected in ("empty required field(s)", "file(s) not found", "points to a different refFasta"):
        assert expected in output, f"missing '{expected}' from the report:\n{output}"


def test_a_missing_samplesheet_names_the_path(nextflow, nextflow_env, tmp_path):
    result = validate(nextflow, nextflow_env, tmp_path, tmp_path / "definitely_absent.tsv")
    output = result.stdout + result.stderr
    assert result.returncode != 0
    assert "Samplesheet not found" in output, output


def test_row_numbers_refer_to_lines_in_the_file(nextflow, nextflow_env, tmp_path):
    """A row number that does not match the user's editor makes the message useless."""
    sheet = tmp_path / "second_row_bad.tsv"
    good = (DATA / "samplesheet.tsv").read_text().splitlines()
    sheet.write_text("\n".join([good[0], good[1], "\tRUN9\t\t\t\t\t\t"]) + "\n")

    result = validate(nextflow, nextflow_env, tmp_path, sheet)
    output = result.stdout + result.stderr
    assert result.returncode != 0
    # Header is line 1 and the first data row line 2, so the broken third line is row 3.
    assert "row 3:" in output, output
