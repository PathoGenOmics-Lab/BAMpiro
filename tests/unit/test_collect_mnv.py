"""Unit tests for bin/collect_mnv.py: the run's codon-level changes from get_MNV.

get_MNV writes a combined amino acid for every codon holding more than one SNV, also when no read
carries them together (two SNVs on different molecules). A codon change is kept only where reads
carry every change of the codon at once.
"""

from __future__ import annotations

import csv
import subprocess
import sys

from conftest import BIN, load_script

cm = load_script("collect_mnv")

COLS = ["Chromosome", "Gene", "Positions", "Reference Bases", "Base Changes", "AA Changes", "SNP AA Changes",
        "Variant Type", "Change Type", "Reference Codon", "SNP Codon", "MNV Codon", "SNP Reads", "MNV Reads",
        "Total Reads", "SNP Frequencies", "MNV Frequencies", "MNV Consequence Shift", "MNV Phasing Support"]
# Rows as get_MNV 1.1.5 writes them for a simulated sample (tests/unit data from a real run of the tool).
ROWS = [
    ["win", "rpoB_Rv0667", "6155", "C", "T", "Ser450Leu", "Ser450Leu", "SNP", "Non-synonymous", "TCG", "TTG", "TTG",
     "102", "0", "103", "0.9903", "0.0000", "-", "-"],
    # on the same reads: His445Glu, where each base alone says His445Asp and His445Gln
    ["win", "rpoB_Rv0667", "6139, 6141", "C, C", "G, G", "His445Glu", "His445Asp, His445Gln", "SNP/MNV",
     "Non-synonymous", "CAC", "GAC, CAG", "GAG", "2, 1", "103", "104", "0.0189, 0.0096", "0.9904", "Concordant",
     "1.0000"],
    # two SNVs of one codon on different molecules: get_MNV still writes a combined change, carried by no read
    ["win", "rpoB_Rv0667", "6094, 6096", "C, G", "A, C", "Leu430Ile", "Leu430Met, Leu430Leu", "SNP/MNV",
     "Non-synonymous", "CTG", "ATG, CTC", "ATC", "89, 26", "0", "115", "0.7739, 0.2261", "0.0000", "Concordant",
     "0.0000"],
    # a minority MNV whose bases alone would name a stop
    ["win", "rpoB_Rv0667", "6184, 6185", "G, A", "T, T", "Glu460Leu", "Glu460*, Glu460Val", "SNP/MNV",
     "Non-synonymous", "GAG", "TAG, GTG", "TTG", "0, 0", "23", "99", "0.0000, 0.0000", "0.2323", "MNV-masked",
     "1.0000"],
    ["win", "NP_215192.1_Rv0678", "24130", "T", "TC", "Glu49Argfs", "-", "INDEL", "Frameshift Indel", "GAT", "-",
     "GAT", "-", "-", "-", "-", "-", "-", "-"],
]


def _tsv(path, rows=ROWS):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(COLS)
        w.writerows(rows)
    return path


def _run(tmp_path, manifest_lines, min_reads=2):
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text("".join("\t".join(x) + "\n" for x in manifest_lines))
    out = tmp_path / "mnv.tsv"
    proc = subprocess.run([sys.executable, str(BIN / "collect_mnv.py"), "--manifest", str(manifest),
                           "--min-reads", str(min_reads), "-o", str(out)], capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    with open(out, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t")), proc.stderr


def test_only_codons_that_reads_carry_whole_are_kept(tmp_path):
    rows, err = _run(tmp_path, [("S1", "win", str(_tsv(tmp_path / "S1.win.MNV.tsv")))])
    assert [r["aa_change"] for r in rows] == ["His445Glu", "Glu460Leu"]
    assert "Leu430Ile" not in {r["aa_change"] for r in rows}, "no read carries it: a change nobody has"
    assert "1 of them masked" in err


def test_a_row_says_what_the_bases_say_alone_and_how_many_reads_carry_it(tmp_path):
    rows, _ = _run(tmp_path, [("S1", "win", str(_tsv(tmp_path / "S1.win.MNV.tsv")))])
    his = rows[0]
    assert (his["sample"], his["reference"], his["gene"], his["positions"]) == ("S1", "win", "rpoB_Rv0667", "6139,6141")
    assert (his["ref_codon"], his["mnv_codon"], his["snp_aa_changes"]) == ("CAC", "GAG", "His445Asp, His445Gln")
    assert (his["mnv_reads"], his["total_reads"], his["mnv_frequency"], his["phasing_support"]) == \
           ("103", "104", "0.9904", "1.0000")
    assert rows[1]["consequence_shift"] == "MNV-masked"


def test_the_read_floor_is_the_one_given(tmp_path):
    rows, _ = _run(tmp_path, [("S1", "win", str(_tsv(tmp_path / "S1.win.MNV.tsv")))], min_reads=50)
    assert [r["aa_change"] for r in rows] == ["His445Glu"]


def test_samples_and_references_come_from_the_manifest_not_the_file_name(tmp_path):
    a, b = _tsv(tmp_path / "x.tsv"), _tsv(tmp_path / "y.tsv", ROWS[1:2])
    rows, _ = _run(tmp_path, [("A.1", "refA", str(a)), ("B", "refB", str(b))])
    assert sorted({(r["sample"], r["reference"]) for r in rows}) == [("A.1", "refA"), ("B", "refB")]


def test_an_empty_or_missing_file_is_skipped_and_counted(tmp_path):
    empty = tmp_path / "empty.tsv"
    empty.write_text("")
    rows, err = _run(tmp_path, [("S1", "win", str(empty)), ("S2", "win", str(tmp_path / "nope.tsv"))])
    assert rows == []
    assert "2 empty or missing file(s)" in err


def test_rows_are_ordered_by_position_then_sample(tmp_path):
    rows, _ = _run(tmp_path, [("S2", "win", str(_tsv(tmp_path / "a.tsv"))), ("S1", "win", str(_tsv(tmp_path / "b.tsv")))])
    assert [(r["positions"], r["sample"]) for r in rows] == [("6139,6141", "S1"), ("6139,6141", "S2"),
                                                             ("6184,6185", "S1"), ("6184,6185", "S2")]


def test_helpers():
    assert cm._list("6139, 6141") == "6139,6141"
    assert cm._clean("-") == "" and cm._clean(".") == "" and cm._clean(" x ") == "x"
    assert cm._int("12") == 12 and cm._int("-") is None
