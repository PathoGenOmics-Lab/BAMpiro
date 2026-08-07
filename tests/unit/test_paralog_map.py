"""Unit tests for bin/paralog_map.py (paralog pairs and paralog-diagnostic sites).

The script aligns nothing. It re-reads the nucmer self-alignment PREPARE_REFERENCE already
produced and recovers what the repeat-masking collapse throws away: which copy aligns to
which, in which orientation, and at which positions the two copies actually differ. So the
whole surface under test is parsing of MUMmer text plus two small filters, and every fixture
here is literal `show-coords -r -c -l -T` / `show-snps -l -r -T -H` output.

The genome behind the fixtures: one 2300 bp contig `chr` with a 600 bp paralog pair at
401-1000 and 1301-1900, the two copies differing at 471/1371 and at 551/1451.
"""

from __future__ import annotations

import os
import subprocess

import pytest

from conftest import load_script

pm = load_script("paralog_map")


# --------------------------------------------------------------------------- fixtures

# show-coords -r -c -l -T writes two banner lines, a blank line, then the column header.
# The header line is padded to 13 fields on purpose: the length guard alone would drop it,
# so this makes the isdigit() guard the thing that has to catch it.
COORDS_BANNER = (
    "/ref/genome.fasta /ref/genome.fasta\n"
    "NUCMER\n"
    "\n"
    "[S1]\t[E1]\t[S2]\t[E2]\t[LEN 1]\t[LEN 2]\t[% IDY]\t[LEN R]\t[LEN Q]\t[COV R]\t[COV Q]\t[TAGS]\t\n"
)

# The self-diagonal first, then the paralog pair reported once per direction.
COORDS_ROWS = (
    "1\t2300\t1\t2300\t2300\t2300\t100.00\t2300\t2300\t100.00\t100.00\tchr\tchr\n"
    "401\t1000\t1301\t1900\t600\t600\t98.00\t2300\t2300\t26.09\t26.09\tchr\tchr\n"
    "1301\t1900\t401\t1000\t600\t600\t98.00\t2300\t2300\t26.09\t26.09\tchr\tchr\n"
)

COORDS = COORDS_BANNER + COORDS_ROWS

# show-snps -H emits no banner. Both directions of the self-alignment are reported, as in a
# real self-alignment, so each copy appears once as acceptor and once as donor.
SNPS = (
    "471\tA\tG\t1371\t71\t471\t1\t1\t2300\t2300\t1\t1\tchr\tchr\n"
    "551\tT\tC\t1451\t11\t551\t1\t1\t2300\t2300\t1\t1\tchr\tchr\n"
    "1371\tG\tA\t471\t71\t471\t1\t1\t2300\t2300\t1\t1\tchr\tchr\n"
    "1451\tC\tT\t551\t11\t551\t1\t1\t2300\t2300\t1\t1\tchr\tchr\n"
)


def _coords_row(s1, e1, s2, e2, idy=98.00, tag1="chr", tag2="chr", lenr=2300, lenq=2300):
    """One show-coords -T data row (13 tab-separated fields)."""
    return "\t".join([str(s1), str(e1), str(s2), str(e2),
                      str(abs(e1 - s1) + 1), str(abs(e2 - s2) + 1), f"{idy:.2f}",
                      str(lenr), str(lenq), "26.09", "26.09", tag1, tag2]) + "\n"


def _snps_row(p1, sub1, sub2, p2, frm2="1", tag1="chr", tag2="chr"):
    """One show-snps -T -H row (14 tab-separated fields)."""
    return "\t".join([str(p1), sub1, sub2, str(p2), "10", "10", "1", "1", "2300", "2300",
                      "1", frm2, tag1, tag2]) + "\n"


# ---------------------------------------------------------------------- parse_coords


def test_self_diagonal_is_dropped():
    """The contig aligned to itself over its whole length is an artefact, not a paralog.

    This is the reason parse_coords exists at all: without the guard, every base of the
    genome is reported as its own paralog and the diagnostic-site set is meaningless.
    """
    pairs = pm.parse_coords(COORDS)

    assert len(pairs) == 2
    assert all(p["length"] == 600 for p in pairs)
    assert not any(p["acc_start"] == p["don_start"] and p["acc_end"] == p["don_end"] for p in pairs)


def test_self_diagonal_guard_needs_same_contig_and_same_interval():
    """The guard is deliberately narrow: identical coordinates on DIFFERENT contigs are a
    genuine paralog pair, and the same contig aligned to itself at other coordinates is the
    signal we are after."""
    text = _coords_row(1, 500, 1, 500, tag1="chrA", tag2="chrB")
    assert len(pm.parse_coords(text)) == 1

    text = _coords_row(1, 500, 1, 500, tag1="chr", tag2="chr")
    assert pm.parse_coords(text) == []


def test_both_directions_of_a_pair_are_reported():
    """Acceptor and donor are not interchangeable downstream, so each direction is a row."""
    pairs = pm.parse_coords(COORDS)

    assert [(p["acc_start"], p["acc_end"], p["don_start"], p["don_end"]) for p in pairs] == [
        (401, 1000, 1301, 1900),
        (1301, 1900, 401, 1000),
    ]
    assert all(p["acceptor"] == "chr" and p["donor"] == "chr" for p in pairs)
    assert all(p["identity"] == 98.0 and p["strand"] == "+" for p in pairs)


def test_banner_and_header_lines_are_skipped():
    """Only the four data-free lines differ between these two inputs."""
    assert pm.parse_coords(COORDS) == pm.parse_coords(COORDS_ROWS)


@pytest.mark.parametrize("junk", [
    "",                                     # blank line
    "\n",
    "not a coords line at all\n",
    "1\t2\t3\n",                            # too few fields
    "1\t2300\t1\t2300\n",                   # numeric but truncated
    "-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\n",   # right width, nothing numeric
])
def test_short_or_garbage_lines_do_not_crash(junk):
    pairs = pm.parse_coords(junk + COORDS_ROWS)
    assert len(pairs) == 2


def test_min_identity_filters_divergent_pairs():
    text = COORDS_ROWS + _coords_row(2000, 2200, 2100, 2300, idy=88.00)

    assert len(pm.parse_coords(text, min_identity=90.0)) == 2
    assert len(pm.parse_coords(text, min_identity=85.0)) == 3
    # The boundary is inclusive: a pair exactly at the threshold is kept.
    assert len(pm.parse_coords(text, min_identity=88.0)) == 3


def test_min_length_filters_short_pairs():
    text = COORDS_ROWS + _coords_row(2000, 2100, 2150, 2250)   # 101 bp

    assert len(pm.parse_coords(text, min_length=200)) == 2
    assert len(pm.parse_coords(text, min_length=100)) == 3
    # Inclusive as well: the acceptor length is e1 - s1 + 1 = 101.
    assert len(pm.parse_coords(text, min_length=101)) == 3


def test_reverse_strand_pair_is_normalised_and_flagged():
    """nucmer reports the QUERY interval with S2 > E2 when the copies are inverted.

    Both halves matter: the interval has to come out ascending so downstream range checks
    work, and the orientation has to survive as `strand`, because on an inverted copy the
    donor position moves the other way along the alignment.
    """
    pairs = pm.parse_coords(_coords_row(401, 1000, 1900, 1301))

    assert len(pairs) == 1
    assert pairs[0]["strand"] == "-"
    assert (pairs[0]["don_start"], pairs[0]["don_end"]) == (1301, 1900)
    assert (pairs[0]["acc_start"], pairs[0]["acc_end"]) == (401, 1000)
    assert pairs[0]["length"] == 600


def test_reverse_strand_reference_interval_is_also_normalised():
    """-r normally sorts the reference interval ascending, but the swap is unconditional."""
    pairs = pm.parse_coords(_coords_row(1000, 401, 1301, 1900))

    assert (pairs[0]["acc_start"], pairs[0]["acc_end"]) == (401, 1000)
    assert pairs[0]["strand"] == "+"          # only S2 > E2 means inverted


# ------------------------------------------------------------------------ parse_snps


def test_parse_snps_reads_diagnostic_positions():
    sites, indels = pm.parse_snps(SNPS)

    assert indels == 0
    assert len(sites) == 4
    assert sites[0] == {"acceptor": "chr", "acc_pos": 471, "acc_base": "A",
                        "donor": "chr", "don_pos": 1371, "don_base": "G", "strand": "+"}
    assert (sites[2]["acc_pos"], sites[2]["acc_base"], sites[2]["don_base"]) == (1371, "G", "A")


def test_parse_snps_skips_indels_and_counts_them():
    """A '.' on either side is an indel: the position stops mapping one-to-one between the
    copies, which is exactly what the tract caller assumes. Kept out of the set, counted so
    the caller can report how many it skipped."""
    text = SNPS + _snps_row(600, ".", "T", 1500) + _snps_row(700, "A", ".", 1600)

    sites, indels = pm.parse_snps(text)

    assert len(sites) == 4
    assert indels == 2
    assert all(s["acc_base"] != "." and s["don_base"] != "." for s in sites)


def test_parse_snps_skips_rows_where_both_bases_are_equal():
    """Not diagnostic: if both copies carry the same base, no read can testify either way.

    Observed behaviour, worth pinning: such a row is added to the SAME counter as indels,
    so the "indel/ambiguous position(s) skipped" message on stderr covers both.
    """
    sites, indels = pm.parse_snps(_snps_row(800, "A", "A", 1700))

    assert sites == []
    assert indels == 1


def test_parse_snps_is_case_insensitive():
    sites, _ = pm.parse_snps(_snps_row(471, "a", "g", 1371))

    assert (sites[0]["acc_base"], sites[0]["don_base"]) == ("A", "G")


def test_parse_snps_reads_strand_from_the_query_frame():
    plus, _ = pm.parse_snps(_snps_row(471, "A", "G", 1371, frm2="1"))
    minus, _ = pm.parse_snps(_snps_row(471, "A", "G", 1371, frm2="-1"))

    assert plus[0]["strand"] == "+"
    assert minus[0]["strand"] == "-"


@pytest.mark.parametrize("junk", [
    "",
    "\n",
    "/ref/genome.fasta /ref/genome.fasta\nNUCMER\n\n",       # banner, if -H were dropped
    "[P1]\t[SUB]\t[SUB]\t[P2]\t[BUFF]\t[DIST]\t[R]\t[Q]\t[LEN R]\t[LEN Q]\t[FRM]\t[FRM]\t[TAGS]\t\n",
    "471\tA\tG\n",                                           # too few fields
])
def test_parse_snps_ignores_banner_header_and_short_lines(junk):
    sites, indels = pm.parse_snps(junk + SNPS)

    assert len(sites) == 4
    assert indels == 0


# ----------------------------------------------------------------- sites_within_pairs


def test_sites_within_pairs_keeps_only_sites_inside_a_pair_and_labels_them():
    pairs = pm.parse_coords(COORDS)
    sites, _ = pm.parse_snps(SNPS)

    kept = pm.sites_within_pairs(sites, pairs)

    assert [s["pair_id"] for s in kept] == [0, 0, 1, 1]
    assert [s["acc_pos"] for s in kept] == [471, 551, 1371, 1451]
    assert all(s["identity"] == 98.0 for s in kept)


def test_sites_within_pairs_drops_sites_outside_every_pair():
    """show-snps reports every difference in the delta, including ones from alignments too
    short or too divergent to have been reported as a pair."""
    pairs = pm.parse_coords(COORDS)
    sites, _ = pm.parse_snps(_snps_row(50, "A", "G", 2250))

    assert pm.sites_within_pairs(sites, pairs) == []


def test_sites_within_pairs_requires_both_ends_inside_the_pair():
    """The acceptor position landing inside a pair is not enough: the donor position has to
    land inside the matching interval of the same pair."""
    pairs = pm.parse_coords(COORDS)
    sites, _ = pm.parse_snps(_snps_row(471, "A", "G", 2250))    # acceptor in pair 0, donor not

    assert pm.sites_within_pairs(sites, pairs) == []


def test_sites_within_pairs_requires_matching_contigs():
    pairs = pm.parse_coords(COORDS)
    sites, _ = pm.parse_snps(_snps_row(471, "A", "G", 1371, tag1="other", tag2="other"))

    assert pm.sites_within_pairs(sites, pairs) == []


def test_a_site_belongs_to_every_pair_that_spans_it():
    """`--maxmatch --nosimplify` emits overlapping and nested alignments on purpose, so a tandem
    repeat produces several pairs over the same bases and one difference is diagnostic for all of
    them.

    This used to stop at the first pair found. On H37Rv that cost 140 sites and left 12 pairs with
    NO diagnostic site at all, among them a 5.8 kb paralog at 98% identity, which the tract caller
    then skipped for having nothing to work with. A locus the tool cannot see is worse than one it
    gets wrong.
    """
    pairs = pm.parse_coords(COORDS_ROWS + _coords_row(401, 1000, 1301, 1900, idy=97.00))
    sites, _ = pm.parse_snps(_snps_row(471, "A", "G", 1371))

    kept = pm.sites_within_pairs(sites, pairs)

    assert sorted(k["pair_id"] for k in kept) == [0, 2]


def test_the_alignment_offset_picks_which_nested_alignment_a_site_belongs_to():
    """Taking every span brings back an ambiguity that stopping early used to hide.

    Two nested alignments can put two different donor positions against the SAME acceptor
    position inside one pair's coordinate box, and only one of them is that pair's own alignment.
    The offset between the copies is constant along an alignment and drifting only with indels,
    so it is what tells them apart. On H37Rv this is 196 positions, all of them resolved.
    """
    pairs = pm.parse_coords(_coords_row(401, 1000, 1301, 1900))     # offset +900
    sites, _ = pm.parse_snps(_snps_row(471, "A", "G", 1371)         # offset +900: this alignment
                             + _snps_row(471, "A", "T", 1500))      # offset +1029: a nested one

    kept = pm.sites_within_pairs(sites, pairs)

    assert len(kept) == 1
    assert (kept[0]["don_pos"], kept[0]["don_base"]) == (1371, "G")


def test_a_site_on_the_other_strand_does_not_belong_to_the_pair():
    """A pair is one alignment in one orientation. A difference reported on the other strand
    came from a different one, however well its coordinates happen to fit."""
    pairs = pm.parse_coords(_coords_row(401, 1000, 1301, 1900))     # a plus-strand pair
    sites, _ = pm.parse_snps(_snps_row(471, "A", "G", 1371, frm2="-1"))

    assert pm.sites_within_pairs(sites, pairs) == []


# ------------------------------------------------------------------------------- _run


def _fake_tool(tmp_path, name, body):
    """A tiny executable stand-in for a MUMmer binary."""
    path = tmp_path / name
    path.write_text("#!/bin/sh\n" + body)
    os.chmod(path, 0o755)
    return str(path)


def test_run_returns_stdout(tmp_path):
    tool = _fake_tool(tmp_path, "show-coords", "printf 'hello\\n'\n")

    assert pm._run([tool]) == "hello\n"


def test_run_raises_with_the_tool_name_status_and_stderr(tmp_path):
    """A MUMmer failure must not reach the pipeline as an empty, plausible-looking TSV."""
    tool = _fake_tool(tmp_path, "show-coords", "echo 'ERROR: could not open delta' >&2\nexit 3\n")

    with pytest.raises(RuntimeError) as excinfo:
        pm._run([tool, "-r", "missing.delta"])

    message = str(excinfo.value)
    assert "show-coords" in message
    assert "(3)" in message
    assert "could not open delta" in message


def test_run_truncates_a_huge_stderr(tmp_path):
    tool = _fake_tool(tmp_path, "show-snps", "python3 -c \"print('x' * 5000)\" >&2\nexit 1\n")

    with pytest.raises(RuntimeError) as excinfo:
        pm._run([tool])

    assert len(str(excinfo.value)) < 600


# ------------------------------------------------------------------------------- main


def _cli_workspace(tmp_path, coords_text=COORDS, snps_text=SNPS):
    """A delta plus fake show-coords/show-snps that replay canned MUMmer output."""
    (tmp_path / "coords.txt").write_text(coords_text)
    (tmp_path / "snps.txt").write_text(snps_text)
    (tmp_path / "self_aln.delta").write_text("/ref/genome.fasta /ref/genome.fasta\nNUCMER\n")

    return {
        "delta": str(tmp_path / "self_aln.delta"),
        "coords": _fake_tool(tmp_path, "show-coords", f'cat "{tmp_path / "coords.txt"}"\n'),
        "snps": _fake_tool(tmp_path, "show-snps", f'cat "{tmp_path / "snps.txt"}"\n'),
        "pairs_out": str(tmp_path / "pairs.tsv"),
        "sites_out": str(tmp_path / "sites.tsv"),
    }


def _cli_argv(ws, *extra):
    return ["--delta", ws["delta"], "--out-pairs", ws["pairs_out"], "--out-sites", ws["sites_out"],
            "--show-coords", ws["coords"], "--show-snps", ws["snps"], *extra]


def _read_tsv(path):
    rows = path.read_text().rstrip("\n").split("\n")
    header = rows[0].split("\t")
    return header, [dict(zip(header, r.split("\t"))) for r in rows[1:]]


def test_main_writes_both_tables(tmp_path, capsys):
    ws = _cli_workspace(tmp_path)

    assert pm.main(_cli_argv(ws)) == 0

    pair_header, pair_rows = _read_tsv(tmp_path / "pairs.tsv")
    site_header, site_rows = _read_tsv(tmp_path / "sites.tsv")

    assert pair_header == pm.PAIR_COLS
    assert site_header == pm.SITE_COLS
    # The self-diagonal is gone; the pair survives once per direction.
    assert len(pair_rows) == 2
    assert len(site_rows) == 4

    assert [r["pair_id"] for r in pair_rows] == ["0", "1"]
    assert [r["n_diagnostic"] for r in pair_rows] == ["2", "2"]
    assert [(r["acc_start"], r["acc_end"], r["don_start"], r["don_end"]) for r in pair_rows] == [
        ("401", "1000", "1301", "1900"),
        ("1301", "1900", "401", "1000"),
    ]

    assert [r["pair_id"] for r in site_rows] == ["0", "0", "1", "1"]
    assert [r["acc_pos"] for r in site_rows] == ["471", "551", "1371", "1451"]
    assert [(r["acc_base"], r["don_base"]) for r in site_rows] == [
        ("A", "G"), ("T", "C"), ("G", "A"), ("C", "T"),
    ]

    assert "2 paralogous pair(s), 4 diagnostic site(s)" in capsys.readouterr().err


def test_main_reports_n_diagnostic_per_pair(tmp_path):
    """n_diagnostic is per pair, not a global count: a pair with no difference between the
    copies gets a zero, which is what tells the caller the locus carries no evidence."""
    coords = COORDS_ROWS + _coords_row(2000, 2299, 1, 300, idy=99.00)
    ws = _cli_workspace(tmp_path, coords_text=coords)

    assert pm.main(_cli_argv(ws)) == 0

    _, pair_rows = _read_tsv(tmp_path / "pairs.tsv")
    assert [r["n_diagnostic"] for r in pair_rows] == ["2", "2", "0"]


def test_main_counts_skipped_indels_on_stderr(tmp_path, capsys):
    snps = SNPS + _snps_row(600, ".", "T", 1500)
    ws = _cli_workspace(tmp_path, snps_text=snps)

    assert pm.main(_cli_argv(ws)) == 0

    _, site_rows = _read_tsv(tmp_path / "sites.tsv")
    assert len(site_rows) == 4
    assert "1 indel/ambiguous position(s) skipped" in capsys.readouterr().err


def test_main_passes_the_filters_through(tmp_path):
    coords = COORDS_ROWS + _coords_row(2000, 2100, 2150, 2250, idy=95.00)
    ws = _cli_workspace(tmp_path, coords_text=coords)

    assert pm.main(_cli_argv(ws, "--min-length", "100", "--min-identity", "96")) == 0

    _, pair_rows = _read_tsv(tmp_path / "pairs.tsv")
    assert len(pair_rows) == 2          # short pair now long enough, but too divergent


def test_main_calls_the_tools_with_the_expected_flags(tmp_path):
    """-C is deliberately absent from show-snps: in a self-alignment every position has an
    ambiguous mapping, so -C would return an empty file."""
    log = tmp_path / "argv.log"
    ws = _cli_workspace(tmp_path)
    for key, name in (("coords", "show-coords"), ("snps", "show-snps")):
        ws[key] = _fake_tool(tmp_path, name + ".logging",
                             f'echo "{name} $@" >> "{log}"\ncat "{tmp_path / (key + ".txt")}"\n')

    assert pm.main(_cli_argv(ws)) == 0

    lines = log.read_text().splitlines()
    assert lines[0].startswith("show-coords -r -c -l -T ")
    assert lines[1].startswith("show-snps -l -r -T -H ")
    assert "-C" not in lines[1].split()
    assert all(line.endswith(ws["delta"]) for line in lines)


def test_main_fails_loudly_when_a_mummer_tool_fails(tmp_path):
    ws = _cli_workspace(tmp_path)
    ws["coords"] = _fake_tool(tmp_path, "show-coords.broken", "echo boom >&2\nexit 2\n")

    with pytest.raises(RuntimeError, match="failed"):
        pm.main(_cli_argv(ws))


def test_script_runs_as_a_subprocess(tmp_path, repo_root):
    """The pipeline invokes it as an executable, not as an import."""
    ws = _cli_workspace(tmp_path)
    res = subprocess.run(["python3", str(repo_root / "bin" / "paralog_map.py"), *_cli_argv(ws)],
                         capture_output=True, text=True)

    assert res.returncode == 0, res.stderr
    assert (tmp_path / "pairs.tsv").exists()
    assert "2 paralogous pair(s)" in res.stderr
