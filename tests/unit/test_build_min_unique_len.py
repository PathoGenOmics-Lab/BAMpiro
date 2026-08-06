"""Unit tests for bin/build_min_unique_len.py (genmap k-sweep -> min_unique_len tracks).

The whole script is one `main()` with no importable helpers, so every test drives it as
a subprocess: a small `.fa.fai`, one bedgraph per probed read length K, and assertions
on the resulting `.npz` (and on the optional repeat BED).

Conventions under test:
  * bedgraph coordinates are 0-based half-open and carry through to the arrays;
  * the array value is the SMALLEST probed K at which the k-mer starting there is
    genome-unique, with the sentinel (default 65535) meaning "never unique";
  * the repeat BED is 1-based inclusive, matching the nucmer / bcftools -T contract.
"""

from __future__ import annotations

import subprocess
import sys

import numpy as np
import pytest

SENTINEL = 65535


# --------------------------------------------------------------------------- helpers


def _fai(tmp_path, contigs) -> str:
    """A minimal samtools .fa.fai; only columns 1 (name) and 2 (length) are read."""
    path = tmp_path / "ref.fa.fai"
    path.write_text("".join("%s\t%d\t6\t60\t61\n" % (name, length) for name, length in contigs))
    return str(path)


def _bedgraph(tmp_path, name, rows) -> str:
    path = tmp_path / name
    path.write_text("".join("%s\t%d\t%d\t%s\n" % row for row in rows))
    return str(path)


def _build(repo_root, tmp_path, fai, bedgraphs, *extra):
    """Run the script; return (CompletedProcess, {contig: array})."""
    prefix = str(tmp_path / "out")
    cmd = [sys.executable, str(repo_root / "bin" / "build_min_unique_len.py"),
           "--fai", fai, "--bedgraphs", *bedgraphs, "--out-prefix", prefix, *extra]
    res = subprocess.run(cmd, capture_output=True, text=True)
    tracks = {}
    if res.returncode == 0:
        with np.load(prefix + ".min_unique_len.npz") as npz:
            tracks = {key: npz[key] for key in npz.files}
    return res, tracks


# ---------------------------------------------------------------- the k-sweep collapse


def test_smallest_unique_k_wins(repo_root, tmp_path):
    """Bedgraphs are processed in ascending K and only positions still at the sentinel
    are filled, so a position unique at both 20 and 50 keeps 20."""
    fai = _fai(tmp_path, [("chr1", 10)])
    b20 = _bedgraph(tmp_path, "k20.bedgraph", [("chr1", 0, 5, "1"), ("chr1", 5, 10, "0")])
    b50 = _bedgraph(tmp_path, "k50.bedgraph", [("chr1", 0, 10, "1")])

    res, tracks = _build(repo_root, tmp_path, fai, ["20:" + b20, "50:" + b50])

    assert res.returncode == 0, res.stderr
    assert list(tracks["chr1"]) == [20] * 5 + [50] * 5


def test_command_line_order_of_the_bedgraphs_does_not_matter(repo_root, tmp_path):
    """The tokens are sorted by K before use, so a caller may list them in any order."""
    fai = _fai(tmp_path, [("chr1", 10)])
    b20 = _bedgraph(tmp_path, "k20.bedgraph", [("chr1", 0, 5, "1"), ("chr1", 5, 10, "0")])
    b50 = _bedgraph(tmp_path, "k50.bedgraph", [("chr1", 0, 10, "1")])

    _, ascending = _build(repo_root, tmp_path, fai, ["20:" + b20, "50:" + b50])
    _, descending = _build(repo_root, tmp_path, fai, ["50:" + b50, "20:" + b20])

    assert list(descending["chr1"]) == list(ascending["chr1"])


def test_never_unique_positions_get_the_sentinel(repo_root, tmp_path):
    fai = _fai(tmp_path, [("chr1", 6)])
    b20 = _bedgraph(tmp_path, "k20.bedgraph", [("chr1", 0, 3, "1"), ("chr1", 3, 6, "0")])
    b50 = _bedgraph(tmp_path, "k50.bedgraph", [("chr1", 0, 6, "0")])

    _, tracks = _build(repo_root, tmp_path, fai, ["20:" + b20, "50:" + b50])

    assert list(tracks["chr1"]) == [20, 20, 20, SENTINEL, SENTINEL, SENTINEL]


def test_k_is_never_stored_as_the_sentinel(repo_root, tmp_path):
    """`kv = min(K, sentinel - 1)`: a probed length at or above the sentinel is clamped,
    so a real (unique) position can never be mistaken for "never unique"."""
    fai = _fai(tmp_path, [("chr1", 4)])
    b30 = _bedgraph(tmp_path, "k30.bedgraph", [("chr1", 0, 4, "1")])

    _, tracks = _build(repo_root, tmp_path, fai, ["30:" + b30], "--sentinel", "30")

    assert list(tracks["chr1"]) == [29, 29, 29, 29]


def test_clamping_still_lets_a_smaller_k_win(repo_root, tmp_path):
    fai = _fai(tmp_path, [("chr1", 6)])
    b20 = _bedgraph(tmp_path, "k20.bedgraph", [("chr1", 0, 2, "1"), ("chr1", 2, 6, "0")])
    b50 = _bedgraph(tmp_path, "k50.bedgraph", [("chr1", 0, 6, "1")])

    _, tracks = _build(repo_root, tmp_path, fai, ["20:" + b20, "50:" + b50], "--sentinel", "30")

    assert list(tracks["chr1"]) == [20, 20, 29, 29, 29, 29]


def test_arrays_are_uint16_and_sized_from_the_fai(repo_root, tmp_path):
    fai = _fai(tmp_path, [("chr1", 10), ("chr2", 4)])
    b20 = _bedgraph(tmp_path, "k20.bedgraph", [("chr1", 0, 10, "1"), ("chr2", 0, 4, "1")])

    _, tracks = _build(repo_root, tmp_path, fai, ["20:" + b20])

    assert sorted(tracks) == ["chr1", "chr2"]
    assert tracks["chr1"].dtype == np.uint16 and tracks["chr1"].shape == (10,)
    assert tracks["chr2"].dtype == np.uint16 and tracks["chr2"].shape == (4,)


# ------------------------------------------------------------------ uniqueness cutoff


def test_uniqueness_uses_a_greater_or_equal_test(repo_root, tmp_path):
    """`value >= --unique-eps` counts as unique, so the threshold itself passes."""
    fai = _fai(tmp_path, [("chr1", 4)])
    rows = [("chr1", 0, 1, "1"), ("chr1", 1, 2, "0.5"),
            ("chr1", 2, 3, "0.4999"), ("chr1", 3, 4, "0")]
    b20 = _bedgraph(tmp_path, "k20.bedgraph", rows)

    _, tracks = _build(repo_root, tmp_path, fai, ["20:" + b20], "--unique-eps", "0.5")

    assert list(tracks["chr1"]) == [20, 20, SENTINEL, SENTINEL]


def test_default_epsilon_rejects_almost_unique_positions(repo_root, tmp_path):
    """The default 0.999999 exists so that genmap's float 1.0 passes while a genuinely
    multi-mapping 0.99 does not."""
    fai = _fai(tmp_path, [("chr1", 4)])
    rows = [("chr1", 0, 1, "1"), ("chr1", 1, 2, "0.999999"),
            ("chr1", 2, 3, "0.9999989"), ("chr1", 3, 4, "0.99")]
    b20 = _bedgraph(tmp_path, "k20.bedgraph", rows)

    _, tracks = _build(repo_root, tmp_path, fai, ["20:" + b20])

    assert list(tracks["chr1"]) == [20, 20, SENTINEL, SENTINEL]


def test_comment_lines_in_a_bedgraph_are_skipped(repo_root, tmp_path):
    fai = _fai(tmp_path, [("chr1", 4)])
    path = tmp_path / "k20.bedgraph"
    path.write_text("#track type=bedGraph\nchr1\t0\t4\t1\n")

    _, tracks = _build(repo_root, tmp_path, fai, ["20:" + str(path)])

    assert list(tracks["chr1"]) == [20] * 4


# ----------------------------------------------------------------------- tail policy


def test_tail_policy_keep_zeroes_never_scored_contig_ends(repo_root, tmp_path):
    """genmap cannot score the last K-1 starts of a contig. Those positions are not
    repetitive, so the default policy gives them 0 (any span >= 0 keeps the read)."""
    fai = _fai(tmp_path, [("chr1", 20)])
    b10 = _bedgraph(tmp_path, "k10.bedgraph", [("chr1", 0, 11, "1")])

    _, tracks = _build(repo_root, tmp_path, fai, ["10:" + b10])

    assert list(tracks["chr1"]) == [10] * 11 + [0] * 9


def test_tail_policy_drop_leaves_the_ends_at_the_sentinel(repo_root, tmp_path):
    fai = _fai(tmp_path, [("chr1", 20)])
    b10 = _bedgraph(tmp_path, "k10.bedgraph", [("chr1", 0, 11, "1")])

    _, tracks = _build(repo_root, tmp_path, fai, ["10:" + b10], "--tail-policy", "drop")

    assert list(tracks["chr1"]) == [10] * 11 + [SENTINEL] * 9


def test_a_repetitive_tail_position_is_not_zeroed(repo_root, tmp_path):
    """Only positions the smallest K never covered are treated as tail; a scored but
    non-unique position keeps the sentinel under either policy."""
    fai = _fai(tmp_path, [("chr1", 6)])
    b10 = _bedgraph(tmp_path, "k10.bedgraph", [("chr1", 0, 4, "0")])

    _, tracks = _build(repo_root, tmp_path, fai, ["10:" + b10])

    assert list(tracks["chr1"]) == [SENTINEL] * 4 + [0, 0]


def test_coverage_comes_only_from_the_smallest_k(repo_root, tmp_path):
    """Pinning current behaviour: `covered` is filled only while processing the smallest
    K, so a region scored at a larger K but not at the smallest one still counts as an
    unscored tail. Real genmap output is contiguous per K, so the pipeline never hits it,
    but a mid-contig gap in the smallest-K bedgraph becomes min_unique_len 0 (fully
    permissive) rather than the sentinel."""
    fai = _fai(tmp_path, [("chr1", 10)])
    b10 = _bedgraph(tmp_path, "k10.bedgraph", [("chr1", 0, 3, "1"), ("chr1", 7, 10, "1")])
    b20 = _bedgraph(tmp_path, "k20.bedgraph", [("chr1", 3, 7, "0")])

    _, tracks = _build(repo_root, tmp_path, fai, ["10:" + b10, "20:" + b20])

    assert list(tracks["chr1"]) == [10, 10, 10, 0, 0, 0, 0, 10, 10, 10]


# --------------------------------------------------------------- fai / bedgraph mismatch


def test_a_contig_missing_from_the_fai_is_skipped_silently(repo_root, tmp_path):
    fai = _fai(tmp_path, [("chr1", 4)])
    b20 = _bedgraph(tmp_path, "k20.bedgraph", [
        ("chr1", 0, 4, "1"),
        ("plasmid", 0, 100, "1"),      # absent from the .fai -> ignored, no warning
    ])

    res, tracks = _build(repo_root, tmp_path, fai, ["20:" + b20])

    assert res.returncode == 0
    assert res.stderr == ""
    assert sorted(tracks) == ["chr1"]
    assert list(tracks["chr1"]) == [20] * 4


def test_a_contig_with_no_bedgraph_rows_is_all_tail(repo_root, tmp_path):
    fai = _fai(tmp_path, [("chr1", 4), ("chr2", 3)])
    b20 = _bedgraph(tmp_path, "k20.bedgraph", [("chr1", 0, 4, "1")])

    _, tracks = _build(repo_root, tmp_path, fai, ["20:" + b20])

    assert list(tracks["chr2"]) == [0, 0, 0]                       # keep policy
    _, dropped = _build(repo_root, tmp_path, fai, ["20:" + b20], "--tail-policy", "drop")
    assert list(dropped["chr2"]) == [SENTINEL] * 3


def test_blank_lines_in_the_fai_are_ignored(repo_root, tmp_path):
    path = tmp_path / "ref.fa.fai"
    path.write_text("chr1\t4\t6\t60\t61\n\n")
    b20 = _bedgraph(tmp_path, "k20.bedgraph", [("chr1", 0, 4, "1")])

    res, tracks = _build(repo_root, tmp_path, str(path), ["20:" + b20])

    assert res.returncode == 0, res.stderr
    assert sorted(tracks) == ["chr1"]


@pytest.mark.parametrize("body", [
    "chr1\t0\t4\t1\n\n",                                           # trailing blank line
    "\nchr1\t0\t4\t1\n",                                           # leading blank line
    "chr1\t0\t2\t1\n   \nchr1\t2\t4\t1\n",                         # whitespace-only line in the middle
])
def test_blank_lines_in_a_bedgraph_are_skipped(repo_root, tmp_path, body):
    """A blank line is "\\n" (never "") when iterating a file, so it has to be skipped on
    its content, not on emptiness, or the four-column unpack raises ValueError."""
    fai = _fai(tmp_path, [("chr1", 4)])
    path = tmp_path / "k20.bedgraph"
    path.write_text(body)

    res, tracks = _build(repo_root, tmp_path, fai, ["20:" + str(path)])

    assert res.returncode == 0, res.stderr
    assert list(tracks["chr1"]) == [20] * 4


# -------------------------------------------------------------------------- --mask-bed


def _mask_rows(path):
    lines = path.read_text().splitlines()
    assert lines[0] == "Chrom\tStart\tEnd\tTag\tComment"
    return lines[1:]


def test_mask_bed_is_one_based_inclusive(repo_root, tmp_path):
    """chr1 is unique at 3 bp over 0-based 0..2 and 15..19 and repetitive in between.
    With a 5 bp window a base is masked when no start in the preceding 5 bp can anchor,
    which is 0-based 7..14, written as 1-based inclusive 8..15."""
    fai = _fai(tmp_path, [("chr1", 20)])
    b3 = _bedgraph(tmp_path, "k3.bedgraph", [
        ("chr1", 0, 3, "1"), ("chr1", 3, 15, "0"), ("chr1", 15, 20, "1")])
    mask = tmp_path / "repeats.bed"

    res, tracks = _build(repo_root, tmp_path, fai, ["3:" + b3],
                         "--mask-bed", str(mask), "--mask-window", "5")

    assert res.returncode == 0, res.stderr
    assert list(tracks["chr1"]) == [3, 3, 3] + [SENTINEL] * 12 + [3] * 5
    assert _mask_rows(mask) == ["chr1\t8\t15\t\t"]


def test_mask_bed_is_empty_when_nothing_is_repetitive(repo_root, tmp_path):
    fai = _fai(tmp_path, [("chr1", 10)])
    b3 = _bedgraph(tmp_path, "k3.bedgraph", [("chr1", 0, 10, "1")])
    mask = tmp_path / "repeats.bed"

    _build(repo_root, tmp_path, fai, ["3:" + b3], "--mask-bed", str(mask), "--mask-window", "5")

    assert _mask_rows(mask) == []                      # header written, no intervals


def test_mask_bed_covers_a_fully_repetitive_contig(repo_root, tmp_path):
    """The `-W` initial value (not -1) is what makes a repeat that starts at the contig
    start mask from base 1."""
    fai = _fai(tmp_path, [("chr1", 10)])
    b3 = _bedgraph(tmp_path, "k3.bedgraph", [("chr1", 0, 10, "0")])
    mask = tmp_path / "repeats.bed"

    _build(repo_root, tmp_path, fai, ["3:" + b3], "--mask-bed", str(mask), "--mask-window", "5")

    assert _mask_rows(mask) == ["chr1\t1\t10\t\t"]


def test_mask_bed_lists_every_contig(repo_root, tmp_path):
    fai = _fai(tmp_path, [("chr1", 10), ("chr2", 10)])
    b3 = _bedgraph(tmp_path, "k3.bedgraph", [("chr1", 0, 10, "1"), ("chr2", 0, 10, "0")])
    mask = tmp_path / "repeats.bed"

    _build(repo_root, tmp_path, fai, ["3:" + b3], "--mask-bed", str(mask), "--mask-window", "5")

    assert _mask_rows(mask) == ["chr2\t1\t10\t\t"]


def test_mask_window_widens_the_masked_block(repo_root, tmp_path):
    """A longer read length considered means more bases count as unmappable."""
    fai = _fai(tmp_path, [("chr1", 20)])
    b3 = _bedgraph(tmp_path, "k3.bedgraph", [
        ("chr1", 0, 3, "1"), ("chr1", 3, 15, "0"), ("chr1", 15, 20, "1")])
    narrow = tmp_path / "narrow.bed"
    wide = tmp_path / "wide.bed"

    _build(repo_root, tmp_path, fai, ["3:" + b3], "--mask-bed", str(narrow), "--mask-window", "5")
    _build(repo_root, tmp_path, fai, ["3:" + b3], "--mask-bed", str(wide), "--mask-window", "9")

    assert _mask_rows(narrow) == ["chr1\t8\t15\t\t"]
    assert _mask_rows(wide) == ["chr1\t12\t15\t\t"]


def test_no_mask_bed_by_default(repo_root, tmp_path):
    fai = _fai(tmp_path, [("chr1", 4)])
    b20 = _bedgraph(tmp_path, "k20.bedgraph", [("chr1", 0, 4, "1")])

    _build(repo_root, tmp_path, fai, ["20:" + b20])

    assert not any(p.suffix == ".bed" for p in tmp_path.iterdir())


# -------------------------------------------------------------------------- CLI knobs


@pytest.mark.parametrize("missing", ["--fai", "--bedgraphs", "--out-prefix"])
def test_required_arguments_are_enforced(repo_root, tmp_path, missing):
    fai = _fai(tmp_path, [("chr1", 4)])
    b20 = _bedgraph(tmp_path, "k20.bedgraph", [("chr1", 0, 4, "1")])
    full = {"--fai": [fai], "--bedgraphs": ["20:" + b20],
            "--out-prefix": [str(tmp_path / "out")]}
    argv = []
    for flag, value in full.items():
        if flag != missing:
            argv += [flag] + value

    res = subprocess.run(
        [sys.executable, str(repo_root / "bin" / "build_min_unique_len.py")] + argv,
        capture_output=True, text=True)

    assert res.returncode != 0
    assert missing in res.stderr


def test_output_file_name_is_the_prefix_plus_a_fixed_suffix(repo_root, tmp_path):
    fai = _fai(tmp_path, [("chr1", 4)])
    b20 = _bedgraph(tmp_path, "k20.bedgraph", [("chr1", 0, 4, "1")])
    prefix = tmp_path / "H37Rv"

    res = subprocess.run(
        [sys.executable, str(repo_root / "bin" / "build_min_unique_len.py"),
         "--fai", fai, "--bedgraphs", "20:" + b20, "--out-prefix", str(prefix)],
        capture_output=True, text=True)

    assert res.returncode == 0, res.stderr
    assert (tmp_path / "H37Rv.min_unique_len.npz").exists()


def test_an_unscored_hole_away_from_the_contig_ends_is_reported(repo_root, tmp_path):
    """A never-scored position becomes 0, i.e. "placeable at any read length", so the read filter
    keeps everything there. At a contig end that is deliberate; in the middle it silently disables
    masking over that stretch, so it has to be visible."""
    fai = _fai(tmp_path, [("chr1", 20)])
    # Scored 0-5 and 15-20, nothing in between.
    bg = _bedgraph(tmp_path, "k35", [("chr1", 0, 5, 1.0), ("chr1", 15, 20, 1.0)])
    res, tracks = _build(repo_root, tmp_path, fai, [f"35:{bg}"])

    assert res.returncode == 0, res.stderr
    assert "unscored position(s) away from the contig ends" in res.stderr
    assert "chr1" in res.stderr


def test_unscored_contig_ends_alone_are_not_reported(repo_root, tmp_path):
    """genmap cannot score a k-mer that runs off the sequence, so a tail is expected and silent."""
    fai = _fai(tmp_path, [("chr1", 20)])
    bg = _bedgraph(tmp_path, "k35", [("chr1", 0, 15, 1.0)])
    res, tracks = _build(repo_root, tmp_path, fai, [f"35:{bg}"])

    assert res.returncode == 0, res.stderr
    assert "unscored position" not in res.stderr
