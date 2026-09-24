"""A real run of the pipeline on a cohort whose truth is known.

Every other leg of the suite stops short of the science: `-stub-run` replaces each script block,
and the unit tests exercise `bin/` one script at a time. This one runs BAMpiro for real, in its
containers, on the cohort tests/e2e/cohort.py simulates, and holds the deliverables against what
the cohort planted: the SNPs and their fractions, the protein change each makes, the codon changes
the reads carry and the one they do not, the indels with their filters, the consensus over a
repeat and over a lost stretch, the distances between the samples, the deletion, and the report's
own copy of the matrices.

It needs nextflow, Java 17+ and Docker, pulls the pinned images and takes minutes, so it runs
only when asked:

    BAMPIRO_E2E=1 python -m pytest tests/e2e -q

BAMPIRO_E2E_DIR keeps the cohort, the results, the trace and the Nextflow log in a directory of
your choice (CI uploads it when the job fails); by default they go to pytest's temporary directory.
"""

from __future__ import annotations

import base64
import csv
import gzip
import json
import os
import re
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

import pytest

import cohort
from conftest import REPO_ROOT, run_pipeline

pytestmark = [pytest.mark.e2e, pytest.mark.nextflow]

PROFILE = "docker,test_e2e"
SAMPLES = ("S1", "S2", "S3")


@dataclass
class Run:
    base: Path
    out: Path
    truth: dict
    trace: Path


def docker_available():
    if not shutil.which("docker"):
        return False
    return subprocess.run(["docker", "info"], capture_output=True).returncode == 0


def pull_images(nextflow, env):
    """Pull the pinned images first, so a slow registry is a slow step and not a failed task."""
    out = subprocess.run([nextflow, "config", "-profile", PROFILE, str(REPO_ROOT)], capture_output=True,
                         text=True, cwd=REPO_ROOT, env=env, timeout=300)
    assert out.returncode == 0, out.stderr[-2000:]
    images = {m for m in re.findall(r"^\s*container\s*=\s*'([^']+)'", out.stdout, re.M)
              if not m.startswith("docker://")}   # the params block keeps Singularity's scheme
    assert images, "no container resolved under the docker profile"
    for image in sorted(images):
        subprocess.run(["docker", "pull", "-q", image], check=True, timeout=1800)


@pytest.fixture(scope="module")
def e2e(nextflow, nextflow_env, tmp_path_factory):
    if os.environ.get("BAMPIRO_E2E") != "1":
        pytest.skip("a real pipeline run: set BAMPIRO_E2E=1 (needs Docker and takes minutes)")
    if not docker_available():
        pytest.skip("docker is not available")
    base = Path(os.environ.get("BAMPIRO_E2E_DIR") or tmp_path_factory.mktemp("e2e")).resolve()
    base.mkdir(parents=True, exist_ok=True)
    truth = cohort.build(base / "cohort")
    pull_images(nextflow, nextflow_env)

    out, trace = base / "results", base / "trace.txt"
    result = run_pipeline(
        nextflow, nextflow_env, base,
        "-profile", PROFILE,
        "--tsv", str(base / "cohort" / "e2e.tsv"),
        "--outdir", str(out),
        # BUILD_MAPPABILITY keeps a storeDir cache: inside the run's directory, so it is built here.
        "--mappability_dir", str(base / "mappability"),
        "-with-trace", str(trace),
        timeout=3600,
    )
    output = result.stdout + result.stderr
    (base / "nextflow.stdout").write_text(output)
    if result.returncode != 0:
        # From Nextflow's own account of the failure: a script that does not compile ("Error <file>:<line>")
        # or a process that failed ("ERROR ~", then its cause and command), whichever comes first.
        m = re.search(r"^(ERROR ~|Error )", output, re.M)
        pytest.fail(f"the pipeline failed (exit {result.returncode}):\n"
                    + (output[m.start():m.start() + 6000] if m else output[-6000:]))
    return Run(base, out, truth, trace)


# ---------------------------------------------------------------------------------------------
# Readers
# ---------------------------------------------------------------------------------------------

def one(root: Path, name: str) -> Path:
    found = sorted(root.rglob(name))
    assert len(found) == 1, f"expected one {name} under {root}, found {[str(p) for p in found]}"
    return found[0]


def tsv(path: Path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def consensus(run: Run, sample: str) -> str:
    text = one(run.out, f"{sample}.{cohort.REF_ID}.consensus.fasta").read_text()
    return "".join(ln.strip() for ln in text.splitlines() if not ln.startswith(">"))


def vcf_records(path: Path):
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if not line.startswith("#"):
                yield line.rstrip("\n").split("\t")


def report_payload(run: Run) -> dict:
    html = one(run.out, "e2e_qc_report.html").read_text()
    m = re.search(r'var REPORT_GZ="([^"]+)"', html)
    assert m, "the report carries no payload"
    return json.loads(gzip.decompress(base64.b64decode(m.group(1))))


def close(value, expected, tol=0.15):
    """A fixed call reads at least 0.9; a minority within tol of the fraction planted."""
    return value >= 0.9 if expected == 1.0 else abs(value - expected) <= tol


# ---------------------------------------------------------------------------------------------
# The run
# ---------------------------------------------------------------------------------------------

def test_every_task_completes(e2e):
    with e2e.trace.open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    failed = [r["name"] for r in rows if r["status"] not in ("COMPLETED", "CACHED")]
    assert failed == [], f"tasks that did not complete: {failed}"
    ran = {r["process"].split(":")[-1] for r in rows}
    needed = {"FASTP_PE", "FASTP_SE", "MAPPING_PE", "MAPPING_SE", "MERGE_AND_MARKDUP", "CALL_FREEBAYES",
              "GET_MNV", "MNV_TABLE", "SNP_MATRIX", "INDEL_MATRIX", "CONSENSUS_FASTA", "SNP_DISTANCES",
              "DEPTH_PROFILE", "QC_REPORT", "DUMP_VERSIONS"}
    assert needed <= ran, f"never ran: {sorted(needed - ran)}"


def test_the_run_records_its_software(e2e):
    versions = one(e2e.out, "software_versions.txt").read_text()
    for tool in ("samtools", "bcftools", "freebayes", "bwa-mem2", "fastp", "genmap", "snpEff"):
        line = next((ln for ln in versions.splitlines() if ln.startswith(f"{tool}: ")), "")
        assert re.search(r"\d+\.\d+", line), f"no version recorded for {tool}: {line!r}"
    for s in SAMPLES:
        manifest = json.loads(one(e2e.out, f"{s}.{cohort.REF_ID}.mnv.manifest.json").read_text())
        assert manifest["tool_version"] == "1.1.5", f"{s}: get_MNV {manifest['tool_version']} ran"


def test_each_sample_gets_a_verdict(e2e):
    flags = {r["sample"]: r["verdict"] for r in tsv(one(e2e.out, "e2e_qc_flags.tsv"))}
    assert set(flags) == set(SAMPLES), flags
    assert set(flags.values()) <= {"PASS", "WARN", "FAIL"}, flags


# ---------------------------------------------------------------------------------------------
# SNPs
# ---------------------------------------------------------------------------------------------

def snp_matrix(run: Run):
    return {int(r["pos"]): r for r in tsv(one(run.out, "e2e_snp_matrix.tsv"))}


def test_the_snp_matrix_holds_the_planted_snps_and_nothing_else(e2e):
    rows = snp_matrix(e2e)
    want = {e["pos"] for e in e2e.truth["snps"]}
    assert set(rows) - want == set(), f"SNPs nobody carries: {sorted(set(rows) - want)}"
    assert want - set(rows) == set(), f"planted SNPs not called: {sorted(want - set(rows))}"
    repeat_snp = e2e.truth["uncalled"][0]["pos"]
    assert repeat_snp not in rows, "a SNP no read can be placed at was called"


def test_every_cell_says_what_the_sample_carries(e2e):
    rows = snp_matrix(e2e)
    wrong = []
    for e in e2e.truth["snps"]:
        r = rows[e["pos"]]
        if (r["ref_allele"], r["alt_allele"]) != (e["ref"], e["alt"]):
            wrong.append((e["pos"], "alleles", r["ref_allele"], r["alt_allele"]))
        for s in SAMPLES:
            af, dp = r[f"{s}|AF"], r[f"{s}|DP"]
            if s in e["af"]:
                if not af or not close(float(af), e["af"][s]):
                    wrong.append((e["pos"], s, f"AF {af}, planted {e['af'][s]}"))
            elif af != "0" or int(dp or 0) < 30:
                wrong.append((e["pos"], s, f"AF {af!r} DP {dp!r}: read without the SNP should be 0 at depth"))
    assert wrong == []


def test_each_snp_names_the_protein_change_it_makes(e2e):
    rows = snp_matrix(e2e)
    wrong = [(e["pos"], rows[e["pos"]]["gene"], rows[e["pos"]]["aa_change"], e["gene"], e["hgvs_p"])
             for e in e2e.truth["snps"] if e["gene"]
             and (rows[e["pos"]]["gene"], rows[e["pos"]]["aa_change"]) != (e["gene"], e["hgvs_p"])]
    assert wrong == [], "pos, gene and change written, then planted"


# ---------------------------------------------------------------------------------------------
# Codon-level changes
# ---------------------------------------------------------------------------------------------

def test_a_codon_change_is_reported_where_reads_carry_it_whole(e2e):
    rows = tsv(one(e2e.out, "e2e_mnv.tsv"))
    got = {(r["sample"], r["positions"], r["aa_change"]) for r in rows}
    want = {(m["sample"], ",".join(map(str, m["positions"])), m["aa_change"]) for m in e2e.truth["mnvs"]}
    assert got == want
    unphased = ",".join(map(str, e2e.truth["unphased_codon"]["positions"]))
    assert all(r["positions"] != unphased for r in rows), "two SNPs no read carries together were merged"


def test_a_codon_change_carries_its_fraction_and_what_the_bases_alone_would_say(e2e):
    by_sample = {r["sample"]: r for r in tsv(one(e2e.out, "e2e_mnv.tsv"))}
    for m in e2e.truth["mnvs"]:
        r = by_sample[m["sample"]]
        assert close(float(r["mnv_frequency"]), m["af"]), (m["sample"], r["mnv_frequency"])
        assert (r["consequence_shift"] == "MNV-masked") == m["masked"], (m["sample"], r["consequence_shift"])


def test_the_snp_matrix_names_the_codon_change_of_its_snps(e2e):
    rows = snp_matrix(e2e)
    for m in e2e.truth["mnvs"]:
        for p in m["positions"]:
            assert rows[p]["codon_change"] == m["aa_change"], (p, rows[p]["codon_change"])
            assert rows[p]["codon_change_samples"].startswith(f"{m['aa_change']}={m['sample']}"), \
                (p, rows[p]["codon_change_samples"])
    for p in e2e.truth["unphased_codon"]["positions"]:
        assert rows[p]["codon_change"] == "", (p, rows[p]["codon_change"])


# ---------------------------------------------------------------------------------------------
# Indels
# ---------------------------------------------------------------------------------------------

def test_the_indel_matrix_holds_each_indel_with_its_fraction_and_filter(e2e):
    rows = {(int(r["pos"]), r["ref_allele"], r["alt_allele"]): r for r in tsv(one(e2e.out, "e2e_indel_matrix.tsv"))}
    want = {(i["pos"], i["ref"], i["alt"]): i for i in e2e.truth["indels"]}
    assert set(rows) == set(want)
    wrong = []
    for key, i in want.items():
        r = rows[key]
        if (r["gene"], r["effect"], r["hgvs_p"]) != (i["gene"], i["effect"], i["hgvs_p"]):
            wrong.append((i["name"], r["gene"], r["effect"], r["hgvs_p"]))
        for s in SAMPLES:
            af, dp, ft = r[f"{s}|AF"], r[f"{s}|DP"], r[f"{s}|FT"]
            if s in i["af"]:
                if not af or not close(float(af), i["af"][s]) or ft != i["filter"][s]:
                    wrong.append((i["name"], s, af, ft))
            elif af != "0" or int(dp or 0) < 30 or ft:
                wrong.append((i["name"], s, af, dp, ft))
        if int(r["n_pass"]) != sum(v == "PASS" for v in i["filter"].values()):
            wrong.append((i["name"], "n_pass", r["n_pass"]))
    assert wrong == []


def test_a_call_below_the_rule_stays_in_the_samples_vcf_marked(e2e):
    """S1 carries S2's frameshift on forward reads only: kept, with its fraction, as LowSupport."""
    ins = e2e.truth["indels"][0]
    recs = [f for f in vcf_records(one(e2e.out, f"S1.{cohort.REF_ID}.indels.ann.vcf.gz"))
            if (int(f[1]), f[3], f[4]) == (ins["pos"], ins["ref"], ins["alt"])]
    assert [f[6] for f in recs] == ["LowSupport"]


# ---------------------------------------------------------------------------------------------
# Consensus, distances and the lost stretch
# ---------------------------------------------------------------------------------------------

def test_the_consensus_writes_what_each_sample_carries(e2e):
    t = e2e.truth
    for s in SAMPLES:
        seq = consensus(e2e, s)
        assert len(seq) == t["genome_length"], (s, len(seq))
        for e in t["snps"]:
            base = seq[e["pos"] - 1]
            af = e["af"].get(s)
            if af == 1.0:
                assert base == e["alt"], (s, e["pos"], base, "a fixed SNP")
            elif af:
                assert base not in "ACGT", (s, e["pos"], base, "a minority written as a base")
            else:
                assert base == e["ref"], (s, e["pos"], base, "read without the SNP")


def test_the_consensus_masks_the_repeat_and_leaves_the_lost_stretch_uncalled(e2e):
    t = e2e.truth
    for s in SAMPLES:
        seq = consensus(e2e, s)
        for a, b in t["repeat"]:
            assert set(seq[a - 1:b]) == {"X"}, (s, a, b, sorted(set(seq[a - 1:b])))
    # The ends of a lost stretch are only as sharp as the sequence either side of the break allows:
    # bases that repeat across the junction align on either side of it (in this cohort, 3 at one
    # end and 1 at the other), so the interior is what must be uncalled.
    lost = t["lost"]
    interior = consensus(e2e, "S3")[lost["start"] - 1 + 10:lost["end"] - 10]
    assert set(interior) == {"-"}, sorted(set(interior))
    # The codon S1 has lost is read by deletion-carrying reads only: never the reference base.
    in_frame = next(i for i in t["indels"] if i["name"] == "in-frame deletion")
    gone = consensus(e2e, "S1")[in_frame["pos"]:in_frame["pos"] + len(in_frame["ref"]) - 1]
    assert not set(gone) & set("ACGT"), gone


def test_distances_count_the_fixed_differences_only(e2e):
    matrix = {r["sample"]: r for r in tsv(one(e2e.out, "e2e_snp_distances.tsv"))}
    got = {f"{a}|{b}": int(matrix[a][b]) for a, b in (k.split("|") for k in e2e.truth["distances"])}
    assert got == e2e.truth["distances"]


def test_the_lost_stretch_is_one_samples_deletion(e2e):
    rows = tsv(one(e2e.out, "e2e_deletions.tsv"))
    private = [r for r in rows if r["class"] == "private"]
    assert len(private) == 1, private
    assert not [r for r in rows if r["class"] == "shared"], rows
    r, lost = private[0], e2e.truth["lost"]
    who, _, span = r["samples"].partition(":")
    start, end = map(int, span.split("-"))
    overlap = min(end, lost["end"]) - max(start, lost["start"]) + 1
    assert who == lost["sample"] and overlap >= 0.9 * (lost["end"] - lost["start"] + 1), r
    assert set(lost["genes"]) <= set(r["genes"].split(",")), r["genes"]


# ---------------------------------------------------------------------------------------------
# The report
# ---------------------------------------------------------------------------------------------

def test_the_report_carries_the_same_matrices(e2e):
    payload = report_payload(e2e)
    snps = payload["snp_matrix"]
    assert set(snps["samples"]) == set(SAMPLES), snps["samples"]
    assert {r["pos"] for r in snps["rows"]} == {e["pos"] for e in e2e.truth["snps"]}
    by_pos = {r["pos"]: r for r in snps["rows"]}
    for m in e2e.truth["mnvs"]:
        for p in m["positions"]:
            assert by_pos[p].get("mnv_aa") == [m["aa_change"]], (p, by_pos[p].get("mnv_aa"))
    indels = payload["indel_matrix"]
    assert {(r["pos"], r["ref"], r["alt"]) for r in indels["rows"]} == \
           {(i["pos"], i["ref"], i["alt"]) for i in e2e.truth["indels"]}
