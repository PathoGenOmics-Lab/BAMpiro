"""End-to-end DAG tests via `-stub-run`.

A stub run executes each process's `stub:` block instead of the real tool, so the entire workflow
is exercised in seconds with no container and no reference data. What it proves is the wiring:
that every channel joins, every declared output is produced and collected, and no process is left
starved or duplicated. It proves nothing about the science, which is what the unit tests are for.

`-profile test` covers the DEFAULT feature set. `-profile test_full` turns on every optional
branch so a single run instantiates every process in the pipeline.
"""

from __future__ import annotations

import csv
import re

import pytest

from conftest import REPO_ROOT, run_pipeline

pytestmark = pytest.mark.nextflow

# Declared in modules/*.nf. CALL_BACKBONE, MERGE_VCFS and CONSENSUS_FASTA are also instantiated
# a second time under a _RAW alias for the virgin-consensus path.
ALIASES = {"CALL_BACKBONE_RAW", "MERGE_VCFS_RAW", "CONSENSUS_FASTA_RAW"}


def declared_processes():
    names = set()
    for module in sorted((REPO_ROOT / "modules").glob("*.nf")):
        names.update(re.findall(r"^process\s+(\w+)", module.read_text(), re.MULTILINE))
    return names


def stub_run(nextflow, nextflow_env, tmp_path, profile):
    trace = tmp_path / "trace.txt"
    result = run_pipeline(
        nextflow,
        nextflow_env,
        tmp_path,
        "-profile", profile,
        "-stub-run",
        "--outdir", str(tmp_path / "out"),
        # BUILD_MAPPABILITY uses storeDir, so a cache shared between tests makes every run after
        # the first skip it. Give each test its own so the coverage assertions stay honest.
        "--mappability_dir", str(tmp_path / "mappability"),
        "-with-trace", str(trace),
        timeout=900,
    )
    return result, trace


def executed_processes(trace):
    with trace.open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    # trace "name" is like "MAPPING_PE (SAMPLE_A__RUN1)"; "process" is the fully qualified name.
    return {row["process"].split(":")[-1].strip() for row in rows}


def test_there_is_a_stub_block_for_every_process():
    """A process without a stub aborts the whole stub run, so this is checked statically too."""
    missing = []
    for module in sorted((REPO_ROOT / "modules").glob("*.nf")):
        text = module.read_text()
        for block in re.split(r"^process\s+", text, flags=re.MULTILINE)[1:]:
            name = block.split(None, 1)[0]
            if "stub:" not in block:
                missing.append(f"{module.name}:{name}")
    assert missing == [], f"processes with no stub: block: {missing}"


@pytest.mark.parametrize("profile", ["test", "test_full"])
def test_the_dag_runs_end_to_end(nextflow, nextflow_env, tmp_path, profile):
    result, trace = stub_run(nextflow, nextflow_env, tmp_path, profile)
    output = result.stdout + result.stderr
    assert result.returncode == 0, f"-profile {profile} failed:\n{output}"
    assert trace.exists(), "no trace file was written"

    with trace.open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert rows, "the run scheduled no tasks at all"
    failed = [r["name"] for r in rows if r["status"] != "COMPLETED"]
    assert failed == [], f"tasks did not complete: {failed}"


def test_the_full_profile_reaches_every_process(nextflow, nextflow_env, tmp_path):
    """Guards against an optional branch silently dropping out of the coverage of the stub run."""
    result, trace = stub_run(nextflow, nextflow_env, tmp_path, "test_full")
    assert result.returncode == 0, result.stdout + result.stderr

    ran = executed_processes(trace)
    never_ran = declared_processes() - ran
    assert never_ran == set(), f"-profile test_full never instantiates: {sorted(never_ran)}"
    assert ALIASES <= ran, f"the virgin-consensus aliases did not run: {sorted(ALIASES - ran)}"


def test_the_default_profile_produces_the_report(nextflow, nextflow_env, tmp_path):
    """The QC report is the pipeline's headline output, and it is on by default.

    This is a regression test: with pathotypr, canonical annotation and the liftover all off,
    which is the default, QC_REPORT used to receive three inputs sharing the placeholder name
    NO_FILE and Nextflow refused to stage them (`input file name collision`).
    """
    result, trace = stub_run(nextflow, nextflow_env, tmp_path, "test")
    assert result.returncode == 0, result.stdout + result.stderr
    assert "QC_REPORT" in executed_processes(trace), "QC_REPORT never ran on the default profile"

    outdir = tmp_path / "out"
    reports = list(outdir.glob("*_qc_report.html"))
    flags = list(outdir.glob("*_qc_flags.tsv"))
    assert reports, f"no QC report was published into {outdir}"
    assert flags, f"no qc_flags.tsv was published into {outdir}"


def test_multi_run_samples_are_merged_once(nextflow, nextflow_env, tmp_path):
    """SAMPLE_A has two runs, so it must map twice but merge and call variants only once."""
    result, trace = stub_run(nextflow, nextflow_env, tmp_path, "test")
    assert result.returncode == 0, result.stdout + result.stderr

    with trace.open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    def tasks_for(process):
        return [r["name"] for r in rows if r["process"].split(":")[-1].strip() == process]

    mapping = [t for t in tasks_for("MAPPING_PE") if "SAMPLE_A" in t]
    merges = [t for t in tasks_for("MERGE_AND_MARKDUP") if "SAMPLE_A" in t]
    assert len(mapping) == 2, f"expected SAMPLE_A to map once per run, got {mapping}"
    assert len(merges) == 1, f"expected SAMPLE_A to merge exactly once, got {merges}"


def test_single_end_and_paired_end_samples_both_flow_through(nextflow, nextflow_env, tmp_path):
    """SAMPLE_C is single-end; the SE branch is easy to break without noticing."""
    result, trace = stub_run(nextflow, nextflow_env, tmp_path, "test")
    assert result.returncode == 0, result.stdout + result.stderr

    ran = executed_processes(trace)
    assert "MAPPING_SE" in ran, "the single-end mapping branch never ran"
    assert "MAPPING_PE" in ran, "the paired-end mapping branch never ran"
    assert "FASTP_SE" in ran and "FASTP_PE" in ran
