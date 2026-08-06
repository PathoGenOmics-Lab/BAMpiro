"""Unit tests for bin/pathotypr_liftover.py (alignment-free k-mer coordinate liftover).

The script is a CLI with three subcommands (markers / apply / lift) built on a set of
pure helpers. The helpers are tested directly; `_lift_chain` is the integration target
because it takes a plain namespace and RETURNS its result instead of printing it.

Sequences are generated from a seeded `random.Random`, so every test is deterministic
while still being long and varied enough that a 21-mer is unique by construction.
"""

from __future__ import annotations

import argparse
import random

import pytest

from conftest import load_script

lift = load_script("pathotypr_liftover")

COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}
CODE = {"A": 0, "C": 1, "G": 2, "T": 3}


def _rc(seq: str) -> str:
    """Reverse complement of an ACGT string."""
    return "".join(COMP[b] for b in reversed(seq))


def _random_seq(n: int, seed: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def _encode(kmer: str) -> int:
    """The same 2-bit packing `_kmer_unique_map` uses (first base in the high bits)."""
    code = 0
    for base in kmer:
        code = (code << 2) | CODE[base]
    return code


def _write_fasta(path, name: str, seq: str, width: int = 60) -> str:
    lines = [">%s" % name]
    lines += [seq[i:i + width] for i in range(0, len(seq), width)]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return str(path)


def _chain_args(source_fasta: str, target_fasta: str, **overrides) -> argparse.Namespace:
    """Every attribute `_lift_chain` reads, with the CLI defaults from `main()`.

    `max_gap` is deliberately None here: the code has an adaptive default
    (2*k + 2*sample + 2) that only kicks in when the attribute is None.
    """
    kwargs = dict(
        source_fasta=source_fasta,
        target_fasta=target_fasta,
        kmer_size=21,
        sample=1,
        indel_tol=0,
        max_gap=None,
        min_chain=10,
        min_density=0.2,
        min_identity=0.9,
        rd_min=50,
        rd_out=None,
        source_contig="source",
    )
    kwargs.update(overrides)
    return argparse.Namespace(**kwargs)


# --------------------------------------------------------------------------- _hash64


def test_hash64_is_deterministic_and_64_bit():
    assert lift._hash64(12345) == lift._hash64(12345)
    assert 0 <= lift._hash64(12345) < (1 << 64)
    assert lift._hash64(0) != lift._hash64(1)


def test_hash64_spreads_small_consecutive_codes():
    """FracMinHash sampling relies on the hash mixing, not on the raw code order."""
    values = [lift._hash64(c) % 7 for c in range(200)]
    assert len(set(values)) == 7


# --------------------------------------------------------------------------- _rc_code


@pytest.mark.parametrize("k", [1, 2, 3, 4, 8, 21])
def test_rc_code_is_an_involution(k):
    """Reverse-complementing twice must return the original 2-bit code."""
    rng = random.Random(k)
    codes = [rng.randrange(1 << (2 * k)) for _ in range(50)]
    for code in codes:
        assert lift._rc_code(lift._rc_code(code, k), k) == code


def test_rc_code_matches_string_reverse_complement():
    rng = random.Random(3)
    for _ in range(50):
        kmer = "".join(rng.choice("ACGT") for _ in range(11))
        assert lift._rc_code(_encode(kmer), 11) == _encode(_rc(kmer))


def test_rc_code_known_values():
    assert lift._rc_code(_encode("AAAA"), 4) == _encode("TTTT")
    assert lift._rc_code(_encode("ACGT"), 4) == _encode("ACGT")   # self-reverse-complement
    assert lift._rc_code(_encode("A"), 1) == _encode("T")
    assert lift._rc_code(_encode("C"), 1) == _encode("G")


# -------------------------------------------------------------------- _kmer_unique_map


def test_kmer_unique_map_k_longer_than_sequence_is_empty():
    assert lift._kmer_unique_map("ACGT", 5) == {}


def test_kmer_unique_map_all_n_is_empty():
    assert lift._kmer_unique_map("NNNNNNNNNN", 3) == {}


def test_kmer_unique_map_centres_are_one_based():
    """k=3 over ACGTA: the three k-mers ACG/CGT/GTA have centres 2, 3 and 4."""
    got = lift._kmer_unique_map("ACGTA", 3)
    assert got == {_encode("ACG"): 2, _encode("CGT"): 3, _encode("GTA"): 4}


def test_kmer_unique_map_drops_a_repeated_kmer():
    """ACGT occurs twice in ACGTTACGT, so it must not be an anchor."""
    got = lift._kmer_unique_map("ACGTTACGT", 4)
    assert _encode("ACGT") not in got
    assert _encode("CGTT") in got


def test_kmer_unique_map_encoding_resets_on_non_acgt():
    """The rolling code must not bridge an N: ACG and TTT survive, no chimera in between."""
    got = lift._kmer_unique_map("ACGNTTT", 3)
    assert got == {_encode("ACG"): 2, _encode("TTT"): 6}


def test_kmer_unique_map_lowercase_is_not_recognised():
    """Callers pass uppercase (via `_read_first_contig`); lowercase resets the code."""
    assert lift._kmer_unique_map("acgta", 3) == {}


def test_kmer_unique_map_sampling_keeps_a_subset():
    seq = _random_seq(4000, seed=42)
    full = lift._kmer_unique_map(seq, 21, sample=1)
    sampled = lift._kmer_unique_map(seq, 21, sample=5)
    assert 0 < len(sampled) < len(full)
    assert set(sampled) <= set(full)
    assert all(full[code] == centre for code, centre in sampled.items())


# ------------------------------------------------------------------------------- _lis


def test_lis_empty_input():
    assert lift._lis([]) == ([], [])


def test_lis_strictly_decreasing_input_gives_a_length_one_chain():
    ca, cb = lift._lis([(1, 5), (2, 4), (3, 3)], increasing=True)
    assert len(ca) == 1 and len(cb) == 1


def test_lis_decreasing_mode_takes_the_whole_run():
    assert lift._lis([(1, 5), (2, 4), (3, 3)], increasing=False) == ([1, 2, 3], [5, 4, 3])


def test_lis_worked_example():
    pairs = [(1, 1), (2, 5), (3, 2), (4, 3), (5, 9)]
    assert lift._lis(pairs, increasing=True) == ([1, 3, 4, 5], [1, 2, 3, 9])


def test_lis_is_strict_so_ties_do_not_extend_a_chain():
    ca, _ = lift._lis([(1, 5), (2, 5), (3, 5)], increasing=True)
    assert len(ca) == 1


# ------------------------------------------------------------------------ _homologous


def test_homologous_identity_context_passes():
    seq = _random_seq(400, seed=1)
    assert lift._homologous(seq, seq, 200, 200, 1, 10, 0.9) is True


def test_homologous_rejects_a_shifted_coordinate():
    seq = _random_seq(400, seed=1)
    assert lift._homologous(seq, seq, 200, 240, 1, 10, 0.9) is False


def test_homologous_checks_each_side_independently():
    """Right half at 0.4 identity, left half at 1.0: the aggregate (0.7) clears 0.6 but
    the per-side rule must still reject it."""
    seq = _random_seq(400, seed=2)
    mutated = list(seq)
    for i in range(200, 206):                       # centre p=200 -> right side is 0-based 200..209
        mutated[i] = "A" if seq[i] != "A" else "C"
    mutated = "".join(mutated)
    assert lift._homologous(seq, mutated, 200, 200, 1, 10, 0.6) is False
    assert lift._homologous(seq, mutated, 200, 200, 1, 10, 0.3) is True


def test_homologous_reverse_orientation_uses_the_reverse_complement():
    seq = _random_seq(400, seed=3)
    rc_seq = _rc(seq)
    p = 200
    t = len(seq) - p + 1                            # reflected 1-based coordinate
    assert lift._homologous(seq, rc_seq, p, t, -1, 10, 0.9) is True
    assert lift._homologous(seq, rc_seq, p, t, 1, 10, 0.9) is False


def test_homologous_refuses_a_trivial_pass_at_a_sequence_edge():
    """Only 4 in-range comparisons against a half-window of 10 -> too little context."""
    short = "ACGTA"
    assert lift._homologous(short, short, 3, 3, 1, 10, 0.9) is False


def test_homologous_accepts_a_genome_edge_when_one_full_side_is_in_range():
    """The width rule is `left + right < w`, so a centre at position 1 of a long sequence
    still passes on its right side alone."""
    seq = _random_seq(400, seed=4)
    assert lift._homologous(seq, seq, 1, 1, 1, 10, 0.9) is True


# ---------------------------------------------------------------------------- _chains


def test_chains_extracts_successive_blocks():
    pairs = sorted([(1, 10), (2, 9), (3, 8), (4, 7), (10, 100), (11, 99), (12, 98)])
    assert lift._chains(pairs, increasing=False, min_anchors=3) == [
        ([1, 2, 3, 4], [10, 9, 8, 7]),
        ([10, 11, 12], [100, 99, 98]),
    ]


def test_chains_stops_below_min_anchors():
    pairs = sorted([(1, 10), (2, 9), (3, 8), (4, 7), (10, 100), (11, 99), (12, 98)])
    assert lift._chains(pairs, increasing=False, min_anchors=5) == []


def test_chains_respects_max_chains():
    pairs = sorted([(1, 10), (2, 9), (3, 8), (4, 7), (10, 100), (11, 99), (12, 98)])
    assert lift._chains(pairs, increasing=False, min_anchors=3, max_chains=1) == [
        ([1, 2, 3, 4], [10, 9, 8, 7])
    ]


def test_chains_on_empty_input():
    assert lift._chains([], increasing=True, min_anchors=2) == []


# --------------------------------------------------------------------- _write_outputs


def test_write_outputs_writes_a_sorted_map(tmp_path):
    good = {20: 500, 5: 100, 6: 101}
    out_map = tmp_path / "lift.map"
    lift._write_outputs(good, str(out_map), None, "ctg", "feat")
    assert out_map.read_text().splitlines() == [
        "src_pos\ttgt_pos",
        "5\t100",
        "6\t101",
        "20\t500",
    ]


def test_write_outputs_bed_collapses_runs_and_is_zero_based_start(tmp_path):
    """Targets 100-102 and 500 collapse into two rows. The BED start is `start - 1`
    (0-based) while the end is the last 1-based position, i.e. an inclusive end that
    doubles as the 0-based half-open end."""
    good = {5: 100, 6: 101, 7: 102, 20: 500}
    out_bed = tmp_path / "lift.bed"
    lift._write_outputs(good, None, str(out_bed), "ctg", "feat")
    assert out_bed.read_text() == "ctg\t99\t102\tfeat\nctg\t499\t500\tfeat\n"


def test_write_outputs_bed_single_position_run(tmp_path):
    out_bed = tmp_path / "one.bed"
    lift._write_outputs({7: 42}, None, str(out_bed), "chrom", "name")
    assert out_bed.read_text() == "chrom\t41\t42\tname\n"


def test_write_outputs_deduplicates_targets(tmp_path):
    """Two sources landing on the same target give one BED position, not two."""
    out_bed = tmp_path / "dup.bed"
    lift._write_outputs({1: 10, 2: 10, 3: 11}, None, str(out_bed), "ctg", "feat")
    assert out_bed.read_text() == "ctg\t9\t11\tfeat\n"


def test_write_outputs_empty_mapping(tmp_path):
    out_map = tmp_path / "empty.map"
    out_bed = tmp_path / "empty.bed"
    lift._write_outputs({}, str(out_map), str(out_bed), "ctg", "feat")
    assert out_map.read_text() == "src_pos\ttgt_pos\n"     # header only
    assert out_bed.read_text() == ""                       # no rows at all


def test_write_outputs_skips_files_that_were_not_requested(tmp_path):
    lift._write_outputs({1: 2}, None, None, "ctg", "feat")
    assert list(tmp_path.iterdir()) == []


# ----------------------------------------------------------------- _read_first_contig


def test_read_first_contig_reads_only_the_first_record(tmp_path):
    fasta = tmp_path / "two.fa"
    fasta.write_text(">one\nacgt\nAACC\n>two\nTTTTTTTT\n")
    assert lift._read_first_contig(fasta) == "ACGTAACC"


def test_read_first_contig_ignores_leading_content_without_a_header(tmp_path):
    fasta = tmp_path / "noheader.fa"
    fasta.write_text("ACGT\nACGT\n")
    assert lift._read_first_contig(fasta) == ""


def test_read_first_contig_strips_whitespace(tmp_path):
    fasta = tmp_path / "ws.fa"
    fasta.write_text(">c\n  ACGT \n\nAC GT\n")
    # Only leading/trailing whitespace is stripped; an interior space is kept verbatim.
    assert lift._read_first_contig(fasta) == "ACGTAC GT"


# ------------------------------------------------------------------------- _positions


def test_positions_plain_list_is_sorted_and_deduplicated(tmp_path):
    path = tmp_path / "pos.txt"
    path.write_text("7\n5\n7\n\n1\n")
    assert lift._positions(path) == [1, 5, 7]


def test_positions_expands_a_bed_interval_to_one_based(tmp_path):
    """BED 0-based half-open [10, 13) becomes 1-based 11, 12, 13."""
    path = tmp_path / "regions.bed"
    path.write_text("chr1\t10\t13\n")
    assert lift._positions(path) == [11, 12, 13]


def test_positions_accepts_extra_bed_columns(tmp_path):
    path = tmp_path / "named.bed"
    path.write_text("chr1\t10\t13\tname\t0\t+\n")
    assert lift._positions(path) == [11, 12, 13]


def test_positions_treats_commas_as_tabs(tmp_path):
    path = tmp_path / "csv.bed"
    path.write_text("chr1,20,22\n")
    assert lift._positions(path) == [21, 22]


def test_positions_skips_header_and_comment_lines(tmp_path):
    path = tmp_path / "hdr.bed"
    path.write_text("# a comment\ntrack name=x\nChrom\tStart\tEnd\nchr1\t0\t2\n")
    assert lift._positions(path) == [1, 2]


def test_positions_skips_unparseable_rows(tmp_path):
    path = tmp_path / "junk.txt"
    path.write_text("chr1\tfoo\tbar\nnot-a-number\n9\n")
    assert lift._positions(path) == [9]


def test_positions_empty_file(tmp_path):
    path = tmp_path / "none.txt"
    path.write_text("")
    assert lift._positions(path) == []


def test_positions_negative_bed_start_is_expanded_verbatim(tmp_path):
    """Pinning current behaviour: `-2` passes the `lstrip("-").isdigit()` test, so a
    malformed BED with a negative start yields non-positive positions instead of being
    rejected. Downstream (`cmd_markers`, `_lift_chain`) drops them anyway."""
    path = tmp_path / "neg.bed"
    path.write_text("chr1\t-2\t1\n")
    assert lift._positions(path) == [-1, 0, 1]


# ------------------------------------------------------- _lift_chain (integration)


def test_lift_chain_identity_maps_positions_onto_themselves(tmp_path):
    seq = _random_seq(600, seed=11)
    src = _write_fasta(tmp_path / "src.fa", "A", seq)
    tgt = _write_fasta(tmp_path / "tgt.fa", "B", seq)
    args = _chain_args(src, tgt)

    positions = [50, 100, 200, 300, 500]
    good, stats = lift._lift_chain(args, positions)

    assert good == {p: p for p in positions}
    assert stats["fwd"] == len(positions)
    assert stats["rev"] == 0
    assert stats["drop"] == 0
    assert stats["rd"] == 0
    assert stats["anchors"] > 0


def test_lift_chain_drops_positions_outside_the_anchored_span(tmp_path):
    """No extrapolation: the last k/2 bases have no k-mer centre, so they are past the
    chain's last anchor and drop instead of receiving a smeared coordinate."""
    seq = _random_seq(600, seed=11)
    src = _write_fasta(tmp_path / "src.fa", "A", seq)
    tgt = _write_fasta(tmp_path / "tgt.fa", "B", seq)
    good, stats = lift._lift_chain(_chain_args(src, tgt), [1, 600])
    assert good == {}
    assert stats["drop"] == 2


def test_lift_chain_tracks_an_insertion_in_the_target(tmp_path):
    """Positions before the insertion keep their coordinate; positions after it shift by
    exactly the inserted length."""
    seq = _random_seq(600, seed=11)
    insert = "GGATCCTTAAGGCCTTAACCGGTTAACCGGTTAACC"
    src = _write_fasta(tmp_path / "src.fa", "A", seq)
    tgt = _write_fasta(tmp_path / "tgt.fa", "B", seq[:300] + insert + seq[300:])
    args = _chain_args(src, tgt)

    good, stats = lift._lift_chain(args, [50, 100, 250, 350, 400, 500])

    assert good[50] == 50
    assert good[100] == 100
    assert good[250] == 250
    assert good[350] == 350 + len(insert)
    assert good[400] == 400 + len(insert)
    assert good[500] == 500 + len(insert)
    assert stats["fwd"] == 6
    assert stats["drop"] == 0


def test_lift_chain_maps_an_inversion_to_the_reflected_coordinate(tmp_path):
    """The target carries seq[200:400] reverse-complemented, so a source position p in
    201..400 lifts to 601 - p via a reverse chain."""
    seq = _random_seq(600, seed=11)
    inverted = seq[:200] + _rc(seq[200:400]) + seq[400:]
    src = _write_fasta(tmp_path / "src.fa", "A", seq)
    tgt = _write_fasta(tmp_path / "tgt.fa", "B", inverted)

    good, stats = lift._lift_chain(_chain_args(src, tgt), [50, 250, 300, 350, 500])

    assert good[50] == 50 and good[500] == 500            # collinear flanks
    assert good[250] == 601 - 250
    assert good[300] == 601 - 300
    assert good[350] == 601 - 350
    assert stats["inv_chains"] == 1
    assert stats["rev"] == 3
    assert stats["fwd"] == 2


def test_lift_chain_reports_an_rd_deletion(tmp_path):
    """100 bp missing from the target is >= rd_min, so it is written to --rd-out in
    SOURCE coordinates as the two flanking anchor centres."""
    seq = _random_seq(600, seed=11)
    src = _write_fasta(tmp_path / "src.fa", "A", seq)
    tgt = _write_fasta(tmp_path / "tgt.fa", "B", seq[:200] + seq[300:])
    rd_out = tmp_path / "rd.bed"
    args = _chain_args(src, tgt, rd_out=str(rd_out), source_contig="H37Rv")

    good, stats = lift._lift_chain(args, [50, 100, 250])

    assert stats["rd"] == 1
    rows = [line.split("\t") for line in rd_out.read_text().splitlines()]
    assert len(rows) == 1
    contig, left, right, tag = rows[0]
    assert contig == "H37Rv" and tag == "RD_deletion"
    assert int(left) < 200 < 300 < int(right)             # the anchors bracket the lost block
    assert good[50] == 50 and good[100] == 100            # flank unaffected
    assert 250 not in good                                # inside the deleted block -> dropped


def test_lift_chain_does_not_write_rd_out_when_unset(tmp_path):
    seq = _random_seq(400, seed=5)
    src = _write_fasta(tmp_path / "src.fa", "A", seq)
    tgt = _write_fasta(tmp_path / "tgt.fa", "B", seq)
    lift._lift_chain(_chain_args(src, tgt), [100])
    assert sorted(p.name for p in tmp_path.iterdir()) == ["src.fa", "tgt.fa"]


def test_lift_chain_sampling_keeps_the_same_coordinates(tmp_path):
    """FracMinHash sampling thins the anchors but must not move a lifted coordinate."""
    seq = _random_seq(600, seed=11)
    insert = "GGATCCTTAAGGCCTTAACCGGTTAACCGGTTAACC"
    src = _write_fasta(tmp_path / "src.fa", "A", seq)
    tgt = _write_fasta(tmp_path / "tgt.fa", "B", seq[:300] + insert + seq[300:])
    positions = [50, 100, 250, 400, 500]

    dense, dense_stats = lift._lift_chain(_chain_args(src, tgt), positions)
    sparse, sparse_stats = lift._lift_chain(_chain_args(src, tgt, sample=4), positions)

    assert sparse == dense
    assert sparse_stats["anchors"] < dense_stats["anchors"]


def test_lift_chain_on_unrelated_sequences_lifts_nothing(tmp_path):
    src = _write_fasta(tmp_path / "src.fa", "A", _random_seq(600, seed=11))
    tgt = _write_fasta(tmp_path / "tgt.fa", "B", _random_seq(600, seed=99))
    good, stats = lift._lift_chain(_chain_args(src, tgt), [100, 200, 300])
    assert good == {}
    assert stats["drop"] == 3


# ------------------------------------------------------------------------- cmd_lift


def _lift_args(tmp_path, positions_file, source_fasta, target_fasta, **overrides):
    """Every attribute `cmd_lift` reads, with the CLI defaults from `main()`."""
    kwargs = dict(
        positions=positions_file,
        source_fasta=source_fasta,
        target_fasta=target_fasta,
        out_map=str(tmp_path / "lift.map"),
        out_bed=str(tmp_path / "lift.bed"),
        contig="target",
        name="blindspot",
        kmer_size=21,
        max_shift=100,
        min_margin=20,
        global_chain=False,
    )
    kwargs.update(overrides)
    return argparse.Namespace(**kwargs)


def test_cmd_lift_anchors_unique_context_kmers(tmp_path, capsys):
    seq = _random_seq(600, seed=11)
    src = _write_fasta(tmp_path / "src.fa", "A", seq)
    tgt = _write_fasta(tmp_path / "tgt.fa", "B", seq)
    positions = tmp_path / "pos.txt"
    positions.write_text("100\n101\n300\n")

    args = _lift_args(tmp_path, str(positions), src, tgt)
    lift.cmd_lift(args)

    assert (tmp_path / "lift.map").read_text() == "src_pos\ttgt_pos\n100\t100\n101\t101\n300\t300\n"
    assert (tmp_path / "lift.bed").read_text() == "target\t99\t101\tblindspot\ntarget\t299\t300\tblindspot\n"
    assert "3 anchored + 0 recovered-by-synteny" in capsys.readouterr().err


def test_cmd_lift_recovers_a_repeated_kmer_by_synteny(tmp_path, capsys):
    """The target carries a second copy of seq[100:200], so the context k-mer of 150 is
    no longer unique there. The nearest anchor's offset predicts the syntenic occurrence
    and the position is recovered instead of dropped.

    Note that 90 and 250 are in the request only to serve as anchors: the anchor set is
    built from the REQUESTED positions, not from the whole genome (see the next test)."""
    seq = _random_seq(600, seed=11)
    src = _write_fasta(tmp_path / "src.fa", "A", seq)
    tgt = _write_fasta(tmp_path / "tgt.fa", "B", seq + seq[100:200])
    positions = tmp_path / "pos.txt"
    positions.write_text("90\n150\n250\n")

    lift.cmd_lift(_lift_args(tmp_path, str(positions), src, tgt))

    assert (tmp_path / "lift.map").read_text() == "src_pos\ttgt_pos\n90\t90\n150\t150\n250\t250\n"
    assert "2 anchored + 1 recovered-by-synteny" in capsys.readouterr().err


def test_cmd_lift_synteny_needs_other_requested_positions_as_anchors(tmp_path, capsys):
    """Pinning current behaviour: `predict()` only sees anchors derived from the
    requested position list, so asking for a single repeated position on its own drops it
    even though the genome is full of usable anchors. Lifting a whole BED is unaffected;
    lifting one coordinate at a time is not equivalent to lifting them together."""
    seq = _random_seq(600, seed=11)
    src = _write_fasta(tmp_path / "src.fa", "A", seq)
    tgt = _write_fasta(tmp_path / "tgt.fa", "B", seq + seq[100:200])
    positions = tmp_path / "pos.txt"
    positions.write_text("150\n")

    lift.cmd_lift(_lift_args(tmp_path, str(positions), src, tgt))

    assert (tmp_path / "lift.map").read_text() == "src_pos\ttgt_pos\n"
    assert "0 anchored + 0 recovered-by-synteny" in capsys.readouterr().err


def test_cmd_lift_drops_a_position_whose_kmer_is_absent_from_the_target(tmp_path, capsys):
    src = _write_fasta(tmp_path / "src.fa", "A", _random_seq(600, seed=11))
    tgt = _write_fasta(tmp_path / "tgt.fa", "B", _random_seq(600, seed=99))
    positions = tmp_path / "pos.txt"
    positions.write_text("300\n")

    lift.cmd_lift(_lift_args(tmp_path, str(positions), src, tgt))

    assert (tmp_path / "lift.map").read_text() == "src_pos\ttgt_pos\n"
    assert (tmp_path / "lift.bed").read_text() == ""
    assert "1 dropped (of 1)" in capsys.readouterr().err


def test_cmd_lift_skips_positions_without_a_full_context_window(tmp_path, capsys):
    """A position closer than k//2 to either end has no centred k-mer, so it never even
    enters the candidate set (it is not counted as a drop either)."""
    seq = _random_seq(600, seed=11)
    src = _write_fasta(tmp_path / "src.fa", "A", seq)
    tgt = _write_fasta(tmp_path / "tgt.fa", "B", seq)
    positions = tmp_path / "pos.txt"
    positions.write_text("1\n5\n600\n300\n")

    lift.cmd_lift(_lift_args(tmp_path, str(positions), src, tgt))

    assert (tmp_path / "lift.map").read_text() == "src_pos\ttgt_pos\n300\t300\n"
    assert "1 dropped (of 1)" not in capsys.readouterr().err


def test_cmd_lift_global_chain_dispatches_to_lift_chain(tmp_path, capsys):
    """--global-chain places positions by anchor interpolation, so it also lifts a
    position whose own context k-mer was destroyed by a nearby insertion."""
    seq = _random_seq(600, seed=11)
    insert = "GGATCCTTAAGGCCTTAACCGGTTAACCGGTTAACC"
    src = _write_fasta(tmp_path / "src.fa", "A", seq)
    tgt = _write_fasta(tmp_path / "tgt.fa", "B", seq[:300] + insert + seq[300:])
    positions = tmp_path / "pos.txt"
    positions.write_text("100\n400\n")

    args = _lift_args(tmp_path, str(positions), src, tgt, global_chain=True,
                      sample=1, indel_tol=0, max_gap=None, min_chain=10, min_density=0.2,
                      min_identity=0.9, rd_min=50, rd_out=None, source_contig="A")
    lift.cmd_lift(args)

    assert (tmp_path / "lift.map").read_text() == "src_pos\ttgt_pos\n100\t100\n400\t436\n"
    assert "lift(anchor-chain): 2 fwd + 0 inverted = 2 lifted" in capsys.readouterr().err


# --------------------------------------------------------------------- cmd_markers


def test_cmd_markers_writes_header_and_reference_bases(tmp_path, capsys):
    seq = _random_seq(100, seed=7)
    fasta = _write_fasta(tmp_path / "ref.fa", "REF", seq)
    positions = tmp_path / "pos.txt"
    positions.write_text("5\n50\n")
    out = tmp_path / "markers.tsv"

    lift.cmd_markers(argparse.Namespace(fasta=fasta, positions=str(positions),
                                        out=str(out), level="lift"))

    assert out.read_text() == "pos\talt\tlevel\n5\t%s\tlift\n50\t%s\tlift\n" % (seq[4], seq[49])
    assert "2 written, 0 skipped" in capsys.readouterr().err


def test_cmd_markers_skips_out_of_range_and_non_acgt(tmp_path, capsys):
    seq = "N" + _random_seq(99, seed=7)[1:]
    fasta = _write_fasta(tmp_path / "ref.fa", "REF", seq)
    positions = tmp_path / "pos.txt"
    positions.write_text("0\n1\n5\n1000\n")
    out = tmp_path / "markers.tsv"

    lift.cmd_markers(argparse.Namespace(fasta=fasta, positions=str(positions),
                                        out=str(out), level="lineage4"))

    rows = out.read_text().splitlines()
    assert rows[0] == "pos\talt\tlevel"
    assert rows[1:] == ["5\t%s\tlineage4" % seq[4]]        # 0 and 1000 out of range, 1 is an N
    assert "1 written, 3 skipped" in capsys.readouterr().err


# ------------------------------------------------------------------------ cmd_apply


def _classify_rows(rows) -> str:
    return "".join("\t".join(r) + "\n" for r in rows)


def test_cmd_apply_writes_map_and_bed_with_offset(tmp_path, capsys):
    classify = tmp_path / "classify.tsv"
    classify.write_text(_classify_rows([
        ("genome", "ACGT", "0", "1000", "111", "lin"),
        ("genome", "ACGT", "0", "1001", "112", "lin"),
        ("genome", "TTTT", "0", "2000", "500", "lin"),
        ("truncated", "row"),                              # < 5 columns -> ignored
        ("genome", "ACGT", "0", "notanint", "113", "lin"),  # unparseable -> ignored
    ]))
    out_map = tmp_path / "out.map"
    out_bed = tmp_path / "out.bed"

    lift.cmd_apply(argparse.Namespace(
        classify=str(classify), target_fasta=None, source_fasta=None, kmer_size=4,
        offset=1, out_map=str(out_map), out_bed=str(out_bed), contig="ANCESTOR",
        name="blindspot"))

    assert out_map.read_text() == "src_pos\ttgt_pos\n111\t1001\n112\t1002\n500\t2001\n"
    # 1001 and 1002 are consecutive -> one run; 2001 is its own run.
    assert out_bed.read_text() == "ANCESTOR\t1000\t1002\tblindspot\nANCESTOR\t2000\t2001\tblindspot\n"
    assert "3 lifted" in capsys.readouterr().err


def test_cmd_apply_drops_a_source_mapping_to_two_targets(tmp_path, capsys):
    classify = tmp_path / "classify.tsv"
    classify.write_text(_classify_rows([
        ("genome", "ACGT", "0", "1000", "111", "lin"),
        ("genome", "ACGT", "0", "2000", "111", "lin"),     # same src, different tgt -> ambiguous
        ("genome", "TTTT", "0", "3000", "222", "lin"),
        ("genome", "TTTT", "0", "3000", "222", "lin"),     # repeated but consistent -> kept
    ]))
    out_map = tmp_path / "out.map"

    lift.cmd_apply(argparse.Namespace(
        classify=str(classify), target_fasta=None, source_fasta=None, kmer_size=4,
        offset=0, out_map=str(out_map), out_bed=None, contig="ctg", name="feat"))

    assert out_map.read_text() == "src_pos\ttgt_pos\n222\t3000\n"
    assert "1 ambiguous" in capsys.readouterr().err


def test_cmd_apply_drops_kmers_that_recur_in_the_target(tmp_path, capsys):
    seq = _random_seq(600, seed=13)
    unique_kmer = seq[100:121]
    duplicated_kmer = seq[200:221]
    source = _write_fasta(tmp_path / "src.fa", "A", seq)
    # The target carries a second copy of `duplicated_kmer`, so its coordinate is ambiguous.
    target = _write_fasta(tmp_path / "tgt.fa", "B", seq + duplicated_kmer)
    classify = tmp_path / "classify.tsv"
    classify.write_text(_classify_rows([
        ("genome", unique_kmer, "0", "111", "111", "lin"),
        ("genome", duplicated_kmer, "0", "211", "211", "lin"),
    ]))
    out_map = tmp_path / "out.map"

    lift.cmd_apply(argparse.Namespace(
        classify=str(classify), target_fasta=target, source_fasta=source, kmer_size=21,
        offset=0, out_map=str(out_map), out_bed=None, contig="ctg", name="feat"))

    assert out_map.read_text() == "src_pos\ttgt_pos\n111\t111\n"
    assert "1 lifted, 1 dropped" in capsys.readouterr().err


def test_cmd_apply_drops_kmers_that_recur_in_the_source(tmp_path, capsys):
    seq = _random_seq(600, seed=13)
    duplicated_kmer = seq[300:321]
    source = _write_fasta(tmp_path / "src.fa", "A", seq + duplicated_kmer)
    target = _write_fasta(tmp_path / "tgt.fa", "B", seq)
    classify = tmp_path / "classify.tsv"
    classify.write_text(_classify_rows([
        ("genome", duplicated_kmer, "0", "311", "311", "lin"),
    ]))
    out_map = tmp_path / "out.map"

    lift.cmd_apply(argparse.Namespace(
        classify=str(classify), target_fasta=target, source_fasta=source, kmer_size=21,
        offset=0, out_map=str(out_map), out_bed=None, contig="ctg", name="feat"))

    assert out_map.read_text() == "src_pos\ttgt_pos\n"
    assert "0 lifted, 1 dropped" in capsys.readouterr().err


def test_cmd_apply_with_no_outputs_requested_writes_nothing(tmp_path, capsys):
    classify = tmp_path / "classify.tsv"
    classify.write_text(_classify_rows([("genome", "ACGT", "0", "10", "20", "lin")]))

    lift.cmd_apply(argparse.Namespace(
        classify=str(classify), target_fasta=None, source_fasta=None, kmer_size=4,
        offset=0, out_map=None, out_bed=None, contig="ctg", name="feat"))

    assert sorted(p.name for p in tmp_path.iterdir()) == ["classify.tsv"]
    assert "1 lifted" in capsys.readouterr().err
