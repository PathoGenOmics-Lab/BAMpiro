"""Pairwise SNP distances between consensus sequences (bin/snp_distances.py).

Two samples differ at a position only when both called a base there and the bases differ: a gap,
a masked position or an ambiguity code is not a difference, whichever sample it is in.
"""

from __future__ import annotations

from conftest import load_script

sd = load_script("snp_distances")


def fasta(path, seq, name=None, wrap=5):
    records = seq if isinstance(seq, list) else [(name or path.stem.split(".")[0], seq)]
    path.write_text("".join(f">{n}\n" + "\n".join(s[i:i + wrap] for i in range(0, len(s), wrap)) + "\n"
                            for n, s in records))
    return str(path)


def run(tmp_path, seqs, refs=None):
    """{(a, b): (snps, compared)} and the square matrix rows, from {sample: sequence}."""
    lines = []
    for s, seq in seqs.items():
        fasta(tmp_path / f"{s}.fa", seq, s)
        lines.append(f"{s}\t{(refs or {}).get(s, 'REF')}\t{s}.fa")
    (tmp_path / "manifest.tsv").write_text("\n".join(lines) + "\n")
    assert sd.main(["--manifest", str(tmp_path / "manifest.tsv"), "--dir", str(tmp_path),
                    "-o", str(tmp_path / "square.tsv"), "--pairs", str(tmp_path / "pairs.tsv")]) == 0
    pairs = {}
    for line in (tmp_path / "pairs.tsv").read_text().splitlines()[1:]:
        a, b, _, snps, comp, _ = line.split("\t")
        pairs[(a, b)] = (int(snps), int(comp))
    square = [r.split("\t") for r in (tmp_path / "square.tsv").read_text().splitlines()]
    return pairs, square


def test_bases_that_differ_where_both_called_are_snps(tmp_path):
    pairs, _ = run(tmp_path, {"A": "ACGTACGTAC", "B": "ACGTTCGTAA", "C": "ACGTACGTAC"})

    assert pairs[("A", "B")][0] == 2
    assert pairs[("A", "C")][0] == 0


def test_a_gap_mask_or_ambiguity_code_is_not_a_difference(tmp_path):
    """The consensus writes '-' where depth was too low, 'X' where the pipeline masked, and an
    IUPAC code at a mixed site. None of them says the samples differ."""
    pairs, _ = run(tmp_path, {"A": "ACGTACGTAC", "B": "-CXTRCGTAC", "C": "TCGTACGTAC"})

    # Position 1 is the only variable one, and B did not call it: B differs from nobody there.
    assert pairs[("A", "B")] == (0, 0)
    assert pairs[("B", "C")] == (0, 0)
    assert pairs[("A", "C")] == (1, 1)


def test_compared_counts_the_variable_positions_both_samples_called(tmp_path):
    pairs, _ = run(tmp_path, {"A": "AAAAAAAAAA", "B": "CAAAAAAAAC", "C": "A-AAAAAAA-"})

    assert pairs[("A", "B")] == (2, 2)
    assert pairs[("A", "C")] == (0, 1), "C did not call position 10, so only position 1 compares"


def test_samples_of_different_references_are_not_compared(tmp_path):
    pairs, square = run(tmp_path, {"A": "ACGTACGTAC", "B": "ACGTACGTAA", "X": "TTTTTTTTTTTT"},
                        refs={"X": "OTHER"})

    assert ("A", "X") not in pairs and ("B", "X") not in pairs
    head = square[0]
    row_a = dict(zip(head[1:], next(r for r in square[1:] if r[0] == "A")[1:]))
    assert row_a == {"A": "0", "B": "1", "X": "NA"}


def test_a_multi_contig_consensus_is_compared_whole(tmp_path):
    """A reference with a plasmid gives a consensus of two records, compared end to end."""
    fasta(tmp_path / "A.fa", [("A_chr", "ACGTA"), ("A_plasmid", "CCCC")])
    fasta(tmp_path / "B.fa", [("B_chr", "ACGTA"), ("B_plasmid", "CCCA")])

    assert sd.main(["--consensus", str(tmp_path / "A.fa"), str(tmp_path / "B.fa"),
                    "-o", str(tmp_path / "sq.tsv"), "--pairs", str(tmp_path / "p.tsv")]) == 0

    row = (tmp_path / "p.tsv").read_text().splitlines()[1].split("\t")
    assert row[:2] == ["A", "B"], "a multi-record consensus is named after its file"
    assert row[3:5] == ["1", "1"]


def test_a_single_record_names_its_sample(tmp_path):
    fasta(tmp_path / "x.consensus.fasta", "ACGT", name="S1")
    fasta(tmp_path / "y.consensus.fasta", "ACGA", name="S2")

    sd.main(["--consensus", str(tmp_path / "x.consensus.fasta"), str(tmp_path / "y.consensus.fasta"),
             "-o", str(tmp_path / "sq.tsv")])

    assert (tmp_path / "sq.tsv").read_text().splitlines()[0] == "sample\tS1\tS2"


def test_the_matrix_is_symmetric_with_a_zero_diagonal(tmp_path):
    _, square = run(tmp_path, {"A": "ACGTACGTAC", "B": "TCGTACGTAA", "C": "ACCTACGTAC"})

    head = square[0][1:]
    m = {r[0]: dict(zip(head, r[1:])) for r in square[1:]}
    assert all(m[s][s] == "0" for s in head)
    assert all(m[a][b] == m[b][a] for a in head for b in head)


def test_a_sample_on_two_references_is_one_row_per_reference(tmp_path):
    """Named once per reference, with each reference's distances in its own row and column."""
    fasta(tmp_path / "A1.fa", "ACGTACGTAC", "A")
    fasta(tmp_path / "B.fa", "ACGTACGTAA", "B")
    fasta(tmp_path / "A2.fa", "TTTTTTTT", "A")
    fasta(tmp_path / "C.fa", "TTTTTTTA", "C")
    (tmp_path / "m.tsv").write_text("A\tR1\tA1.fa\nB\tR1\tB.fa\nA\tR2\tA2.fa\nC\tR2\tC.fa\n")

    assert sd.main(["--manifest", str(tmp_path / "m.tsv"), "--dir", str(tmp_path),
                    "-o", str(tmp_path / "sq.tsv")]) == 0

    rows = [r.split("\t") for r in (tmp_path / "sq.tsv").read_text().splitlines()]
    assert rows[0] == ["sample", "A@R1", "B", "A@R2", "C"]
    m = {r[0]: dict(zip(rows[0][1:], r[1:])) for r in rows[1:]}
    assert m["A@R1"]["B"] == "1" and m["A@R2"]["C"] == "1" and m["A@R1"]["C"] == "NA"


def test_a_truncated_consensus_is_left_out_by_name_and_the_rest_compared(tmp_path, capsys):
    fasta(tmp_path / "A.fa", "ACGTAC", "A")
    fasta(tmp_path / "B.fa", "ACGTAA", "B")
    fasta(tmp_path / "C.fa", "ACG", "C")
    (tmp_path / "m.tsv").write_text("C\tR\tC.fa\nA\tR\tA.fa\nB\tR\tB.fa\n")

    assert sd.main(["--manifest", str(tmp_path / "m.tsv"), "--dir", str(tmp_path),
                    "-o", str(tmp_path / "sq.tsv"), "--pairs", str(tmp_path / "p.tsv")]) == 0

    err = capsys.readouterr().err
    assert "C: " in err and "left out" in err, "the short file, not the next one, is named"
    assert (tmp_path / "p.tsv").read_text().splitlines()[1].split("\t")[:4] == ["A", "B", "R", "1"]
