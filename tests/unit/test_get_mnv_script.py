"""The GET_MNV script hands get_MNV a GFF it can read.

get_MNV reads the GFF as plain text and fails on a gzipped one ("stream did not contain valid
UTF-8"), while the samplesheet takes either: the SnpEff build and every other GFF reader of the
pipeline open both. A stub run skips the script block, so the script is rendered here and run in
bash, with a get_mnv on PATH that records what it was given.
"""

from __future__ import annotations

import gzip
import os
import re
import stat
import subprocess

from conftest import REPO_ROOT

MODULE = REPO_ROOT / "modules" / "mnv.nf"
GFF = "##gff-version 3\nc1\tRefSeq\tCDS\t1\t300\t.\t+\t0\tID=cds-1;gene=rpoB\n"


def get_mnv_script():
    text = MODULE.read_text()
    body = text[text.index("process GET_MNV"):]
    m = re.search(r'script:\s*"""(.*?)"""', body, re.S)
    assert m, "GET_MNV has no script block"
    return m.group(1)


def render(block, values):
    """What Nextflow runs: ${...} interpolated, \\$ and \\\\ unescaped."""
    def sub(m):
        t = m.group(0)
        if t == "\\\\":
            return "\\"
        if t == "\\$":
            return "$"
        return str(values[t[2:-1]])
    return re.sub(r"\\\\|\\\$|\$\{[^}]+\}", sub, block)


def run_with(tmp_path, gff_name):
    for f in ("codon_calls.vcf.gz", "codon_calls.vcf.gz.tbi", "S1.bam", "S1.bam.bai", "reference.fa"):
        (tmp_path / f).write_text("")
    bindir = tmp_path / "bin"
    bindir.mkdir()
    fake = bindir / "get_mnv"
    fake.write_text('#!/usr/bin/env bash\nprintf \'%s\\n\' "$@" > get_mnv.args\n')
    fake.chmod(fake.stat().st_mode | stat.S_IEXEC)
    values = {"calls": "codon_calls.vcf.gz", "calls_tbi": "codon_calls.vcf.gz.tbi", "sampleId": "S1", "refId": "r",
              "bam": "S1.bam", "ref_fa": "reference.fa", "gff": gff_name, "task.cpus": 2,
              "params.mnv_gff_features": "CDS", "params.mnv_translation_table": 11,
              "params.freebayes_min_base_qual": 20, "params.freebayes_min_map_qual": 30}
    env = dict(os.environ, PATH=f"{bindir}{os.pathsep}{os.environ['PATH']}")
    proc = subprocess.run(["bash", "-c", render(get_mnv_script(), values)], cwd=tmp_path, env=env,
                          capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    args = (tmp_path / "get_mnv.args").read_text().splitlines()
    return args, tmp_path / args[args.index("--gff") + 1]


def test_a_gzipped_gff_reaches_get_mnv_as_plain_text(tmp_path):
    with gzip.open(tmp_path / "genes.gff.gz", "wt") as fh:
        fh.write(GFF)
    args, gff = run_with(tmp_path, "genes.gff.gz")
    assert gff.name != "genes.gff.gz", "get_MNV was handed the gzipped GFF, which it cannot read"
    assert gff.read_text() == GFF


def test_a_plain_gff_is_handed_over_as_given(tmp_path):
    (tmp_path / "genes.gff").write_text(GFF)
    args, gff = run_with(tmp_path, "genes.gff")
    assert gff.name == "genes.gff"
    assert args[args.index("--vcf") + 1] == "S1.r.vcf.gz", "get_MNV names its outputs after the calls it reads"
