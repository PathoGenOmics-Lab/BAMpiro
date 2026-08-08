"""The gene-conversion stage end to end: paralog map, tract caller, cohort pass.

Every other test file takes one of the three scripts and hands it inputs made for it. This one
runs the chain the pipeline runs, so each script reads what the previous one actually wrote:
MUMmer output goes into `paralog_map.py`, its sites file into `gene_conversion.py`, and both of
its outputs into `gconv_cohort.py`. That is where the failures no unit test can see live, because
a column renamed at one end and read by name at the other passes both files' own tests.

The alignments are built rather than aligned, so this needs only samtools. The MUMmer tools are
stood in for with committed text, which is what `--show-coords` and `--show-snps` exist for.

The genome behind it: a 3000 bp contig with a 600 bp paralog pair at 401-1000 and 1601-2200,
differing at four positions and with one base the donor does not have. A conversion is implanted
over the middle of the acceptor copy.
"""

from __future__ import annotations

import shutil
import subprocess

import pytest

from conftest import load_script

pm = load_script("paralog_map")
gc = load_script("gene_conversion")
gcc = load_script("gconv_cohort")

CONTIG = "chr"
# Acceptor position -> (its own base, the donor's). The middle three are the tract.
DIFFS = [(500, "A", "G"), (600, "T", "C"), (700, "G", "A"), (800, "C", "T")]
DELETION = (750, "A")          # a base at 750 the donor does not have
TRACT = (600, 750, 700)        # what the conversion copies: two substitutions and the deletion


@pytest.fixture
def samtools():
    exe = shutil.which("samtools")
    if not exe:
        pytest.skip("samtools is not on PATH")
    return exe


def _tool(tmp_path, name, text):
    """A stand-in for a MUMmer binary that prints committed output."""
    path = tmp_path / name
    path.write_text("#!/bin/sh\ncat <<'OUT'\n" + text + "OUT\n")
    path.chmod(0o755)
    return str(path)


def _coords():
    row = "\t".join(["401", "1000", "1601", "2200", "600", "600", "97.00",
                     "3000", "3000", "20.00", "20.00", CONTIG, CONTIG])
    return row + "\n"


def _snps():
    out = []
    for pos, acc, don in DIFFS:
        out.append("\t".join([str(pos), acc, don, str(pos + 1200), "50", "50", "1", "1",
                              "3000", "3000", "1", "1", CONTIG, CONTIG]))
    out.append("\t".join([str(DELETION[0]), DELETION[1], ".", str(DELETION[0] + 1200),
                          "50", "50", "1", "1", "3000", "3000", "1", "1", CONTIG, CONTIG]))
    return "\n".join(out) + "\n"


def _reads(converted, n=25, length=400, start=450):
    """Molecules over the whole locus, carrying the converted alleles when `converted`."""
    records = []
    for i in range(n):
        seq = []
        # Walk the acceptor: a base per position, and at the deletion site either the base or,
        # once converted, nothing at all, which the CIGAR has to say with a D.
        ops = []
        for p in range(start, start + length):
            base = "N"
            for pos, acc, don in DIFFS:
                if p == pos:
                    base = don if (converted and TRACT[0] <= p <= TRACT[1]) else acc
            if p == DELETION[0]:
                if converted:
                    ops.append("D")
                    continue
                base = DELETION[1]
            ops.append("M")
            seq.append(base)
        # collapse the operation run-lengths into a CIGAR
        cig, run, prev = "", 0, ops[0]
        for op in ops + [None]:
            if op == prev:
                run += 1
            else:
                cig += f"{run}{prev}"
                prev, run = op, 1
        records.append("\t".join([f"r{i}", "0", CONTIG, str(start), "0", cig, "*", "0", "0",
                                  "".join(seq), "I" * len(seq)]))
    return records


def _bam(tmp_path, records, name):
    sam = tmp_path / f"{name}.sam"
    sam.write_text(f"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:{CONTIG}\tLN:3000\n"
                   + "\n".join(records) + "\n")
    bam = tmp_path / f"{name}.bam"
    subprocess.run(["samtools", "view", "-b", "-o", str(bam), str(sam)], check=True,
                   capture_output=True)
    subprocess.run(["samtools", "index", str(bam)], check=True, capture_output=True)
    return str(bam)


def _run_chain(tmp_path, samtools, samples):
    """paralog_map -> gene_conversion per sample -> gconv_cohort, as the pipeline runs them."""
    pairs, sites = tmp_path / "pairs.tsv", tmp_path / "sites.tsv"
    assert pm.main(["--delta", str(tmp_path / "x.delta"),
                    "--out-pairs", str(pairs), "--out-sites", str(sites),
                    "--show-coords", _tool(tmp_path, "coords.sh", _coords()),
                    "--show-snps", _tool(tmp_path, "snps.sh", _snps())]) == 0

    tracts, loci = [], []
    for name, converted in samples:
        bam = _bam(tmp_path, _reads(converted), name)
        t, lo = tmp_path / f"{name}.tracts.tsv", tmp_path / f"{name}.loci.tsv"
        assert gc.main(["--sites", str(sites), "--bam", bam, "--sample", name,
                        "-o", str(t), "--output-loci", str(lo),
                        "--samtools", samtools, "--min-sites", "2"]) == 0
        tracts.append(str(t))
        loci.append(str(lo))

    cohort = tmp_path / "cohort.tsv"
    assert gcc.main(["--tracts", *tracts, "--loci", *loci, "-o", str(cohort)]) == 0
    return sites, tracts, cohort


def _rows(path):
    text = open(path).read().splitlines()
    body = [ln for ln in text if not ln.startswith("#")]
    header = body[0].split("\t")
    return [dict(zip(header, ln.split("\t"))) for ln in body[1:]]


def test_the_chain_finds_the_implanted_tract(tmp_path, samtools):
    """One converted sample among four, run through all three scripts in order."""
    samples = [("s1", True)] + [(f"s{i}", False) for i in range(2, 6)]
    sites, tracts, cohort = _run_chain(tmp_path, samtools, samples)

    called = [r for r in _rows(cohort)
              if r["sample"] == "s1" and r["cohort_verdict"] == "gene_conversion"]
    assert called, "the converted sample produced no call"
    assert (int(called[0]["start"]), int(called[0]["end"])) == (TRACT[0], TRACT[1])

    others = [r for r in _rows(cohort) if r["sample"] != "s1"
              and r["cohort_verdict"] == "gene_conversion"]
    assert not others, "the unconverted samples must stay quiet"


def test_the_deletion_between_the_copies_reaches_the_caller_as_a_marker(tmp_path, samtools):
    """The whole point of the map emitting deletions: the site has to survive being written to a
    file, read back, matched against a CIGAR `D`, and priced by the model."""
    sites, tracts, cohort = _run_chain(tmp_path, samtools, [("s1", True), ("s2", False)])

    assert any(r["kind"] == "del" for r in _rows(sites)), "the map did not emit the deletion"
    call = next(r for r in _rows(cohort) if r["sample"] == "s1"
                and r["cohort_verdict"] == "gene_conversion")
    # the tract spans the deletion, so the marker is one of the ones supporting it
    assert int(call["start"]) <= DELETION[0] <= int(call["end"])
    assert int(call["n_sites"]) >= 3


def test_an_event_the_whole_cohort_shows_is_demoted_by_the_chain(tmp_path, samtools):
    """The cohort pass reading real per-sample output, not a fixture written for it."""
    samples = [(f"s{i}", True) for i in range(1, 7)]
    _, _, cohort = _run_chain(tmp_path, samtools, samples)

    verdicts = {r["cohort_verdict"] for r in _rows(cohort)}
    assert verdicts == {"reference_artifact"}


def test_every_column_the_cohort_pass_adds_survives_the_round_trip(tmp_path, samtools):
    """A column renamed at one end and read by name at the other passes both files' own tests."""
    _, _, cohort = _run_chain(tmp_path, samtools, [("s1", True), ("s2", False)])
    rows = _rows(cohort)

    assert rows
    for column in gc.COLUMNS + gcc.COHORT_COLUMNS:
        assert column in rows[0], f"{column} did not survive the chain"


def test_the_donor_coordinates_reach_the_cohort_pass_and_are_used(tmp_path, samtools):
    """Asserting on the BEHAVIOUR, because asserting on the column names cannot catch this.

    `gconv_cohort` reads the donor span by name and falls back to "no span" when it is absent,
    which is right for a file written before the column existed and wrong as a way to notice the
    column was renamed. Renaming it at both ends passed the round-trip check above, and the
    donor resolution silently put every relationship in its own candidate. So this asks the
    question the other test cannot: did the coordinates actually arrive and do their job?
    """
    _, tracts, cohort = _run_chain(tmp_path, samtools, [("s1", True), ("s2", False)])

    called = next(r for r in _rows(cohort) if r["sample"] == "s1"
                  and r["cohort_verdict"] == "gene_conversion")
    assert called["don_start"] and called["don_end"]
    # the donor stretch sits 1200 bases along, which is the offset the fixture alignment has
    assert int(called["don_start"]) - int(called["start"]) == 1200
    assert int(called["don_end"]) - int(called["end"]) == 1200
    # and it is a STRETCH: collapsing it to a point is what the cohort pass then merges donors on
    assert int(called["don_end"]) > int(called["don_start"])
    # and having arrived, it was used: one relationship, so one candidate and nothing to choose
    assert called["n_donors"] == "1"
    assert called["donor_call"] == "only candidate"
    assert called["is_representative"] == "1"


def test_the_genes_a_tract_lands_on_reach_the_cohort_file(tmp_path, samtools):
    """The annotation has to survive the chain, not merely work when called directly.

    A tract reported as coordinates is a finding nobody can act on, and the columns that make it
    one are written per sample and read back twice: once by the cohort pass and once by the
    report. A GFF read into a shape the row builder does not use fails silently, with every
    annotation column empty and no error anywhere.
    """
    gff = tmp_path / "ref.gff"
    gff.write_text("##gff-version 3\n"
                   f"{CONTIG}\t.\tCDS\t401\t1000\t.\t+\t0\tID=c1;gene=ppeA;locus_tag=Rv0001\n")
    fasta = tmp_path / "ref.fa"
    # The acceptor copy runs 401-1000, so the CDS covers it and every diagnostic site is in frame.
    fasta.write_text(f">{CONTIG}\n" + "ACG" * 1000 + "\n")

    pairs, sites = tmp_path / "pairs.tsv", tmp_path / "sites.tsv"
    assert pm.main(["--delta", str(tmp_path / "x.delta"),
                    "--out-pairs", str(pairs), "--out-sites", str(sites),
                    "--show-coords", _tool(tmp_path, "coords.sh", _coords()),
                    "--show-snps", _tool(tmp_path, "snps.sh", _snps())]) == 0

    bam = _bam(tmp_path, _reads(True), "s1")
    out = tmp_path / "s1.tracts.tsv"
    assert gc.main(["--sites", str(sites), "--bam", bam, "--sample", "s1", "-o", str(out),
                    "--samtools", samtools, "--min-sites", "2",
                    "--gff", str(gff), "--reference", str(fasta)]) == 0

    called = next(r for r in _rows(out) if r["verdict"] == "gene_conversion")
    assert called["genes"] == "Rv0001", "the tract did not pick up the gene it sits in"
    assert int(called["n_syn"]) + int(called["n_nonsyn"]) > 0, "no codon was scored"


def test_the_settings_of_every_step_reach_the_cohort_file(tmp_path, samtools):
    _, _, cohort = _run_chain(tmp_path, samtools, [("s1", True), ("s2", False)])
    notes = [ln for ln in open(cohort).read().splitlines() if ln.startswith("#")]

    assert any(n.startswith("# gene_conversion.py") for n in notes)
    assert any(n.startswith("# gconv_cohort.py") for n in notes)
