"""Unit tests for bin/filter_reads_mappability.py (length-aware mappability read filter).

`span()` is importable and tested directly. Everything else lives in a closure inside
`main()`, so the filter is exercised end to end as a subprocess: a hand-written SAM on
stdin, the filtered SAM on stdout, and a numpy `.npz` mappability track on disk.

Track convention (see bin/build_min_unique_len.py): one uint16 array per contig,
0-based, value = the shortest read length that is genome-unique at that start, with
65535 meaning "never unique" (always repetitive).
"""

from __future__ import annotations

import subprocess
import sys

import numpy as np
import pytest

from conftest import load_script

filt = load_script("filter_reads_mappability")

SENTINEL = 65535
KMIN = 35                      # the script default; also what these tests pass explicitly


# --------------------------------------------------------------------------- helpers


def _track(tmp_path, arrays, name="mul.npz") -> str:
    """Write a mappability .npz, one uint16 array per contig."""
    path = tmp_path / name
    np.savez(path, **{k: np.asarray(v, dtype=np.uint16) for k, v in arrays.items()})
    return str(path)


def _default_track(tmp_path) -> str:
    """chr1: unique at 30 bp everywhere except 1-based 51..100, which is a repeat."""
    arr = np.full(200, 30, dtype=np.uint16)
    arr[50:100] = SENTINEL
    return _track(tmp_path, {"chr1": arr})


def _sam(name, flag, rname="chr1", pos=1, cigar="50M", rnext="*", pnext=0, tags=()):
    fields = [name, str(flag), rname, str(pos), "60", cigar, rnext, str(pnext), "0",
              "A" * 10, "I" * 10]
    return "\t".join(fields + list(tags)) + "\n"


HEADER = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:200\n@PG\tID:bwa\tPN:bwa\n"


def _run(repo_root, track, sam_text, *extra):
    cmd = [sys.executable, str(repo_root / "bin" / "filter_reads_mappability.py"),
           "--mul", track, *extra]
    return subprocess.run(cmd, input=sam_text, capture_output=True, text=True)


def _kept(res):
    """Read names of the alignment records that survived, in output order.

    The exit status is asserted here so that an unnoticed crash can never masquerade as
    "every read was filtered out".
    """
    assert res.returncode == 0, res.stderr
    return [line.split("\t")[0] for line in res.stdout.splitlines() if not line.startswith("@")]


# ------------------------------------------------------------------------------ span


@pytest.mark.parametrize("cigar,expected", [
    ("10M2I5M", 15),          # I does not consume the reference
    ("5S10M5S", 10),          # neither does S
    ("10M5N10M", 25),         # N does
    ("10M2D5M", 17),          # D does
    ("5=3X2M", 10),           # = and X do
    ("5H10M", 10),            # hard clip does not
    ("10P5M", 5),             # padding does not
    ("76M", 76),
])
def test_span_sums_only_reference_consuming_ops(cigar, expected):
    assert filt.span(cigar) == expected


def test_span_of_a_missing_cigar_is_zero():
    assert filt.span("*") == 0


@pytest.mark.parametrize("cigar", ["", "10", "not-a-cigar", "MMM", "10m"])
def test_span_of_a_malformed_cigar_is_zero_without_raising(cigar):
    # No validation happens: unknown letters are simply not counted, and a run of
    # digits with no operator contributes nothing. Lowercase 'm' is NOT an operator.
    assert filt.span(cigar) == 0


def test_span_ignores_unknown_operators_but_keeps_the_valid_ones():
    assert filt.span("5Z10M") == 10


# --------------------------------------------------------------------------- headers


def test_header_lines_pass_through_unchanged(repo_root, tmp_path):
    track = _default_track(tmp_path)
    res = _run(repo_root, track, HEADER, "--kmin", str(KMIN))
    assert res.returncode == 0
    assert res.stdout == HEADER


def test_unmodelled_contig_in_the_header_only_warns(repo_root, tmp_path):
    track = _default_track(tmp_path)
    header = HEADER + "@SQ\tSN:plasmid\tLN:5000\n"
    res = _run(repo_root, track, header, "--kmin", str(KMIN))
    assert res.returncode == 0
    assert res.stdout == header
    assert "plasmid" in res.stderr
    assert "fail-open" in res.stderr


def test_strict_contigs_aborts_on_an_unmodelled_contig(repo_root, tmp_path):
    track = _default_track(tmp_path)
    header = HEADER + "@SQ\tSN:plasmid\tLN:5000\n"
    res = _run(repo_root, track, header, "--kmin", str(KMIN), "--strict-contigs")
    assert res.returncode != 0
    assert "plasmid" in res.stderr
    assert "--strict-contigs" in res.stderr


def test_strict_contigs_is_silent_when_every_contig_is_modelled(repo_root, tmp_path):
    track = _default_track(tmp_path)
    res = _run(repo_root, track, HEADER, "--strict-contigs")
    assert res.returncode == 0
    assert res.stderr == ""


# ----------------------------------------------------------------- single-end reads


def test_single_end_read_at_a_unique_locus_is_kept(repo_root, tmp_path):
    track = _default_track(tmp_path)
    sam = HEADER + _sam("unique", 0, pos=11, cigar="50M")     # min_unique_len = 30
    res = _run(repo_root, track, sam, "--kmin", str(KMIN))
    assert _kept(res) == ["unique"]


def test_single_end_read_inside_a_repeat_is_dropped(repo_root, tmp_path):
    track = _default_track(tmp_path)
    sam = HEADER + _sam("repeat", 0, pos=61, cigar="50M")     # min_unique_len = 65535
    res = _run(repo_root, track, sam, "--kmin", str(KMIN))
    assert _kept(res) == []


def test_single_end_read_shorter_than_its_min_unique_len_is_dropped(repo_root, tmp_path):
    """40 bp span at a locus that needs 50 bp: ambiguous, and 40 >= kmin so the
    soft-clip escape hatch does not apply."""
    track = _track(tmp_path, {"chr1": np.full(200, 50, dtype=np.uint16)})
    sam = HEADER + _sam("tooshort", 0, pos=11, cigar="40M")
    res = _run(repo_root, track, sam, "--kmin", str(KMIN))
    assert _kept(res) == []


def test_clipped_read_below_kmin_at_a_resolvable_locus_is_kept(repo_root, tmp_path):
    """FIX1: a 20 bp reference span below the smallest probed k is soft clipping, not
    ambiguity, so a resolvable (non-sentinel) locus keeps the read."""
    track = _default_track(tmp_path)
    sam = HEADER + _sam("clipped", 0, pos=11, cigar="20M80S")
    res = _run(repo_root, track, sam, "--kmin", str(KMIN))
    assert _kept(res) == ["clipped"]


def test_clipped_read_below_kmin_inside_a_repeat_is_still_dropped(repo_root, tmp_path):
    track = _default_track(tmp_path)
    sam = HEADER + _sam("clipped_repeat", 0, pos=61, cigar="20M80S")
    res = _run(repo_root, track, sam, "--kmin", str(KMIN))
    assert _kept(res) == []


def test_read_starting_past_the_end_of_the_track_is_dropped(repo_root, tmp_path):
    track = _default_track(tmp_path)
    sam = HEADER + _sam("offtrack", 0, pos=5000, cigar="50M")
    res = _run(repo_root, track, sam, "--kmin", str(KMIN))
    assert _kept(res) == []


# ------------------------------------------------------------------------ flag drops


@pytest.mark.parametrize("flag,label", [
    (4, "unmapped"),
    (256, "secondary"),
    (512, "qcfail"),
    (1024, "duplicate"),
    (2048, "supplementary"),
])
def test_flag_masked_records_are_dropped(repo_root, tmp_path, flag, label):
    track = _default_track(tmp_path)
    sam = HEADER + _sam(label, flag, pos=11, cigar="50M")
    res = _run(repo_root, track, sam, "--kmin", str(KMIN))
    assert _kept(res) == []


def test_keep_unmapped_readmits_only_plain_unmapped_records(repo_root, tmp_path):
    track = _default_track(tmp_path)
    sam = (HEADER
           + _sam("unmapped", 4, rname="*", pos=0, cigar="*")
           + _sam("unmapped_secondary", 4 | 256, pos=11)
           + _sam("unmapped_supplementary", 4 | 2048, pos=11)
           + _sam("duplicate", 1024, pos=11))
    res = _run(repo_root, track, sam, "--kmin", str(KMIN), "--keep-unmapped")
    assert _kept(res) == ["unmapped"]


def test_without_keep_unmapped_the_unmapped_record_is_gone(repo_root, tmp_path):
    track = _default_track(tmp_path)
    sam = HEADER + _sam("unmapped", 4, rname="*", pos=0, cigar="*")
    res = _run(repo_root, track, sam, "--kmin", str(KMIN))
    assert _kept(res) == []


# ------------------------------------------------------------------ contig fail-open


def test_read_on_a_contig_absent_from_the_track_is_kept(repo_root, tmp_path):
    """Finding2: masking is simply disabled on an unmodelled contig rather than
    silently deleting every read there."""
    track = _default_track(tmp_path)
    sam = HEADER + "@SQ\tSN:plasmid\tLN:5000\n" + _sam("onplasmid", 0, rname="plasmid", pos=1)
    res = _run(repo_root, track, sam, "--kmin", str(KMIN))
    assert res.returncode == 0
    assert _kept(res) == ["onplasmid"]


def test_read_whose_mate_is_on_an_unmodelled_contig_is_kept(repo_root, tmp_path):
    """The self anchor fails (repeat), the mate contig is unknown -> fail open."""
    track = _default_track(tmp_path)
    sam = HEADER + _sam("matefaropen", 1 | 2 | 64, pos=61, cigar="50M",
                        rnext="plasmid", pnext=100, tags=("MC:Z:50M",))
    res = _run(repo_root, track, sam, "--kmin", str(KMIN))
    assert _kept(res) == ["matefaropen"]


# --------------------------------------------------------------------- mate rescue


def test_repeat_read_is_rescued_by_an_anchored_mate(repo_root, tmp_path):
    """Proper pair (0x2), mate mapped, MC:Z gives the mate a 50 bp span at a locus that
    needs 30 -> the whole fragment is kept."""
    track = _default_track(tmp_path)
    sam = HEADER + _sam("rescued", 1 | 2 | 64, pos=61, cigar="50M",
                        rnext="=", pnext=11, tags=("MC:Z:50M",))
    res = _run(repo_root, track, sam, "--kmin", str(KMIN))
    assert _kept(res) == ["rescued"]


def test_rescue_requires_the_proper_pair_flag(repo_root, tmp_path):
    """FIX3: a discordant read sitting in the wrong repeat copy must not be re-admitted
    by its mate, so 0x2 is mandatory."""
    track = _default_track(tmp_path)
    sam = HEADER + _sam("discordant", 1 | 64, pos=61, cigar="50M",
                        rnext="=", pnext=11, tags=("MC:Z:50M",))
    res = _run(repo_root, track, sam, "--kmin", str(KMIN))
    assert _kept(res) == []


def test_rescue_requires_a_mapped_mate(repo_root, tmp_path):
    track = _default_track(tmp_path)
    sam = HEADER + _sam("mateunmapped", 1 | 2 | 8 | 64, pos=61, cigar="50M",
                        rnext="=", pnext=11, tags=("MC:Z:50M",))
    res = _run(repo_root, track, sam, "--kmin", str(KMIN))
    assert _kept(res) == []


def test_rescue_requires_the_mc_tag(repo_root, tmp_path):
    """Without MC:Z the mate's reference span is unknown, so the read cannot be rescued."""
    track = _default_track(tmp_path)
    sam = HEADER + _sam("nomc", 1 | 2 | 64, pos=61, cigar="50M", rnext="=", pnext=11)
    res = _run(repo_root, track, sam, "--kmin", str(KMIN))
    assert _kept(res) == []


def test_rescue_finds_the_mc_tag_among_other_tags(repo_root, tmp_path):
    track = _default_track(tmp_path)
    sam = HEADER + _sam("taggy", 1 | 2 | 64, pos=61, cigar="50M", rnext="=", pnext=11,
                        tags=("NM:i:0", "MD:Z:50", "MC:Z:50M", "RG:Z:s1"))
    res = _run(repo_root, track, sam, "--kmin", str(KMIN))
    assert _kept(res) == ["taggy"]


def test_rescue_fails_when_the_mate_is_also_in_the_repeat(repo_root, tmp_path):
    track = _default_track(tmp_path)
    sam = HEADER + _sam("bothrepeat", 1 | 2 | 64, pos=61, cigar="50M",
                        rnext="=", pnext=71, tags=("MC:Z:50M",))
    res = _run(repo_root, track, sam, "--kmin", str(KMIN))
    assert _kept(res) == []


def test_rescue_fails_when_the_mate_position_is_off_the_track(repo_root, tmp_path):
    track = _default_track(tmp_path)
    sam = HEADER + _sam("mateoff", 1 | 2 | 64, pos=61, cigar="50M",
                        rnext="=", pnext=9000, tags=("MC:Z:50M",))
    res = _run(repo_root, track, sam, "--kmin", str(KMIN))
    assert _kept(res) == []


def test_self_anchored_pair_is_kept_without_looking_at_the_mate(repo_root, tmp_path):
    """The self anchor short-circuits, so neither 0x2 nor MC:Z is needed."""
    track = _default_track(tmp_path)
    sam = HEADER + _sam("selfok", 1 | 64, pos=11, cigar="50M", rnext="=", pnext=61)
    res = _run(repo_root, track, sam, "--kmin", str(KMIN))
    assert _kept(res) == ["selfok"]


def test_mate_rnext_equals_sign_resolves_to_the_read_contig(repo_root, tmp_path):
    """RNEXT '=' and the explicit contig name must behave identically."""
    track = _default_track(tmp_path)
    explicit = HEADER + _sam("explicit", 1 | 2 | 64, pos=61, cigar="50M",
                             rnext="chr1", pnext=11, tags=("MC:Z:50M",))
    res = _run(repo_root, track, explicit, "--kmin", str(KMIN))
    assert _kept(res) == ["explicit"]


# -------------------------------------------------------------------------- ordering


def test_output_preserves_input_order_and_only_drops(repo_root, tmp_path):
    track = _default_track(tmp_path)
    records = [
        _sam("a_unique", 0, pos=11),
        _sam("b_repeat", 0, pos=61),
        _sam("c_unique", 0, pos=111),
        _sam("d_secondary", 256, pos=11),
        _sam("e_unique", 0, pos=151),
        _sam("f_repeat", 0, pos=99),
    ]
    res = _run(repo_root, track, HEADER + "".join(records), "--kmin", str(KMIN))
    assert _kept(res) == ["a_unique", "c_unique", "e_unique"]
    # Drops only: every surviving line is byte-identical to its input line.
    kept_lines = [ln + "\n" for ln in res.stdout.splitlines() if not ln.startswith("@")]
    assert kept_lines == [records[0], records[2], records[4]]


def test_headers_stay_ahead_of_the_records(repo_root, tmp_path):
    track = _default_track(tmp_path)
    res = _run(repo_root, track, HEADER + _sam("r", 0, pos=11), "--kmin", str(KMIN))
    lines = res.stdout.splitlines()
    assert lines[:3] == HEADER.splitlines()
    assert lines[3].startswith("r\t")


# ------------------------------------------------------------------------- CLI knobs


def test_custom_sentinel_is_honoured(repo_root, tmp_path):
    """With --sentinel 999 the value 999 means "never unique", so the soft-clip escape
    hatch must not fire there."""
    arr = np.full(200, 999, dtype=np.uint16)
    track = _track(tmp_path, {"chr1": arr})
    sam = HEADER + _sam("clipped", 0, pos=11, cigar="20M80S")
    res = _run(repo_root, track, sam, "--kmin", str(KMIN), "--sentinel", "999")
    assert _kept(res) == []
    # Without the override, 999 is an ordinary length and the clipped read is rescued.
    res2 = _run(repo_root, track, sam, "--kmin", str(KMIN))
    assert _kept(res2) == ["clipped"]


def test_kmin_controls_the_soft_clip_escape_hatch(repo_root, tmp_path):
    track = _default_track(tmp_path)
    sam = HEADER + _sam("clipped", 0, pos=11, cigar="20M80S")
    assert _kept(_run(repo_root, track, sam, "--kmin", "35")) == ["clipped"]
    assert _kept(_run(repo_root, track, sam, "--kmin", "15")) == []


def test_missing_mul_argument_is_an_error(repo_root, tmp_path):
    cmd = [sys.executable, str(repo_root / "bin" / "filter_reads_mappability.py")]
    res = subprocess.run(cmd, input="", capture_output=True, text=True)
    assert res.returncode != 0
    assert "--mul" in res.stderr


def test_a_blank_line_in_the_stream_crashes(repo_root, tmp_path):
    """Pinning current behaviour: the record parser indexes fields unconditionally, so a
    stray blank line raises IndexError instead of being skipped. samtools never emits
    one, but the filter is not defensive about it."""
    track = _default_track(tmp_path)
    res = _run(repo_root, track, HEADER + "\n" + _sam("r", 0, pos=11), "--kmin", str(KMIN))
    assert res.returncode != 0
    assert "IndexError" in res.stderr
