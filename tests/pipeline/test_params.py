"""Parameter declaration, validation and --help.

Nextflow accepts any `--flag` into `params` whether or not the pipeline declares it, so without
validation a typo runs to completion with the default silently in place. main.nf reads the declared
set straight out of nextflow.config, which is also what these tests check against.
"""

from __future__ import annotations

import re

import pytest

from conftest import DATA, REPO_ROOT, run_pipeline

pytestmark = pytest.mark.nextflow

CONFIG = REPO_ROOT / "nextflow.config"


def declared_params() -> set[str]:
    """The names inside the top-level `params { ... }` block of nextflow.config."""
    text = CONFIG.read_text()
    start = text.index("params {")
    depth, end = 0, None
    for i in range(text.index("{", start), len(text)):
        if text[i] == "{":
            depth += 1
        elif text[i] == "}":
            depth -= 1
            if depth == 0:
                end = i
                break
    assert end is not None, "unterminated params block in nextflow.config"
    body = text[start:end]
    return set(re.findall(r"^\s*([A-Za-z_][A-Za-z0-9_]*)\s*=", body, re.MULTILINE))


def test_every_param_referenced_in_the_pipeline_is_declared():
    """`params.foo` with no declaration is null at runtime, which usually reads as "feature off"."""
    declared = declared_params()
    # Set by Nextflow itself or by a profile rather than by the params block.
    builtin = {"outdir"}

    used = {}
    for source in [REPO_ROOT / "main.nf", *sorted((REPO_ROOT / "modules").glob("*.nf"))]:
        text = source.read_text()
        for match in re.finditer(r"params\.([A-Za-z_][A-Za-z0-9_]*)", text):
            # Skip method calls on the params map itself, e.g. params.keySet(). Checked after the
            # match rather than with a lookahead, which would just backtrack to a shorter name.
            if text[match.end():].lstrip().startswith("("):
                continue
            used.setdefault(match.group(1), source.name)

    undeclared = {name: where for name, where in used.items()
                  if name not in declared and name not in builtin}
    assert undeclared == {}, f"referenced but never declared in nextflow.config: {undeclared}"


def test_help_lists_the_parameters_and_exits_cleanly(nextflow, nextflow_env, tmp_path):
    result = run_pipeline(nextflow, nextflow_env, tmp_path, "--help", timeout=300)
    output = result.stdout + result.stderr

    assert result.returncode == 0, output
    # A representative parameter from three different sections.
    for expected in ("--tsv", "--kraken2_db", "--freebayes_ploidy"):
        assert expected in output, f"--help does not mention {expected}:\n{output}"
    assert "Profiles:" in output


def test_help_stays_in_step_with_the_config(nextflow, nextflow_env, tmp_path):
    """--help is generated from nextflow.config, so a new parameter appears without being listed twice."""
    result = run_pipeline(nextflow, nextflow_env, tmp_path, "--help", timeout=300)
    listed = set(re.findall(r"^\s*--(\w+)", result.stdout + result.stderr, re.MULTILINE))
    missing = declared_params() - listed
    assert missing == set(), f"declared but absent from --help: {sorted(missing)}"


@pytest.mark.parametrize(
    ("typo", "suggestion"),
    [
        ("treads", "--threads"),
        ("outdirr", "--outdir"),
        ("kraken_db", "--kraken2_db"),
    ],
)
def test_a_misspelled_parameter_is_rejected_with_a_suggestion(
    nextflow, nextflow_env, tmp_path, typo, suggestion
):
    result = run_pipeline(
        nextflow, nextflow_env, tmp_path,
        "-profile", "test", "-stub-run",
        "--outdir", str(tmp_path / "out"),
        f"--{typo}", "1",
        timeout=300,
    )
    output = result.stdout + result.stderr
    assert result.returncode != 0, f"--{typo} was accepted:\n{output}"
    assert f"--{typo}" in output
    assert suggestion in output, f"no suggestion for --{typo}:\n{output}"


def test_a_declared_parameter_is_accepted(nextflow, nextflow_env, tmp_path):
    """The validator must not be so eager that it rejects a real flag."""
    result = run_pipeline(
        nextflow, nextflow_env, tmp_path,
        "-profile", "test", "-stub-run",
        "--outdir", str(tmp_path / "out"),
        "--mappability_dir", str(tmp_path / "mappability"),
        "--threads", "2",
        "--freebayes_ploidy", "1",
        timeout=900,
    )
    assert result.returncode == 0, result.stdout + result.stderr


def test_a_missing_samplesheet_parameter_says_what_to_pass(nextflow, nextflow_env, tmp_path):
    result = run_pipeline(nextflow, nextflow_env, tmp_path, "-stub-run", timeout=300)
    output = result.stdout + result.stderr
    assert result.returncode != 0
    assert "--tsv is required" in output, output
    assert "sampleId" in output, "the error should name the samplesheet columns"


def test_the_startup_banner_reports_where_the_run_will_execute(nextflow, nextflow_env, tmp_path):
    """Running on a cluster login node by accident is the easiest mistake to make."""
    result = run_pipeline(
        nextflow, nextflow_env, tmp_path,
        "-profile", "test", "-stub-run",
        "--outdir", str(tmp_path / "out"),
        "--mappability_dir", str(tmp_path / "mappability"),
        timeout=900,
    )
    output = result.stdout + result.stderr
    assert result.returncode == 0, output
    assert "Executor         : local" in output, output
    assert "Profile(s)       : test" in output, output
    assert str(DATA / "samplesheet.tsv") in output


def test_resource_caps_let_an_oversized_request_run(nextflow, nextflow_env, tmp_path):
    """Kraken2 asks for 80 GB, which no laptop has and the local executor will never schedule.

    Only the success direction is asserted: whether an UNCAPPED run fails depends on how much
    memory the host happens to have, so testing that would pass on a laptop and fail on a big
    workstation. `-profile test_full` enables Kraken, and 1 GB is available anywhere.
    """
    result = run_pipeline(
        nextflow, nextflow_env, tmp_path,
        "-profile", "test_full", "-stub-run",
        "--outdir", str(tmp_path / "out"),
        "--mappability_dir", str(tmp_path / "mappability"),
        "--max_cpus", "1",
        "--max_memory", "1.GB",
        timeout=900,
    )
    output = result.stdout + result.stderr
    assert result.returncode == 0, output
    assert "exceeds available" not in output, output
