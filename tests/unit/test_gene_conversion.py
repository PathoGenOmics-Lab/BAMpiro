"""Unit tests for bin/gene_conversion.py (gene conversion tract detection).

Two things are being tested, and they fail in different ways.

`read_bases_at` is arithmetic: it walks a CIGAR to line reference positions up with read
offsets. When it is wrong it is silently wrong, reading the right position off the wrong
base, so it is tested operation by operation against reads whose bases are laid out to make
an off-by-anything visible.

The detector is judgement. Reads that mismap from the donor produce the same per-site
picture as a real conversion: donor bases where the reference expects acceptor bases. The
tests for `classify` are the specification of how the two are told apart, and each of the
four verdict routes gets its own test.
"""

from __future__ import annotations

import shutil
import subprocess

import pytest

from conftest import load_script

gc = load_script("gene_conversion")

CONTIG = "chr"
HIGH_BQ = "I"                       # Phred 40
LOW_BQ = "&"                        # Phred 5, under the default floor of 13


# --------------------------------------------------------------------------- helpers


def _sam_line(name, pos, cigar, seq, qual="*", flag=0, rname=CONTIG, mapq=0):
    """One SAM alignment record. MAPQ defaults to 0 because paralogous regions are exactly
    where an aligner assigns MAPQ 0, and the tool deliberately does not filter on it."""
    return "\t".join([name, str(flag), rname, str(pos), str(mapq), cigar,
                      "*", "0", "0", seq, qual])


def _per_read(spec):
    """{read_name: {pos: 'acceptor'|'donor'|'other'}} written out literally."""
    return {name: dict(calls) for name, calls in spec.items()}


def _ev(**overrides):
    """Evidence for a textbook conversion, so each classify test changes one thing."""
    ev = {"n_sites": 5, "start": 110, "end": 150, "span_bp": 41, "n_sites_outside": 6,
          "donor_af_in": 0.98, "donor_af_outside": 0.0, "min_depth": 20,
          "cis_reads": 8, "breakpoint_reads": 3, "donor_only_reads": 0}
    ev.update(overrides)
    return ev


# -------------------------------------------------------------------------- parse_cigar


def test_parse_cigar_splits_into_length_operation_pairs():
    assert gc.parse_cigar("10S90M") == [(10, "S"), (90, "M")]
    assert gc.parse_cigar("50M5I2D3=1X4H") == [(50, "M"), (5, "I"), (2, "D"), (3, "="),
                                               (1, "X"), (4, "H")]
    assert gc.parse_cigar("*") == []
    assert gc.parse_cigar("") == []


# ------------------------------------------------------------------------ read_bases_at
#
# Every read below starts at reference position 1000 and marks the bases the test cares
# about with a distinctive letter, so a shifted read offset returns a visibly wrong base
# rather than a plausible one.


def test_read_bases_at_plain_match():
    seq = "A" * 49 + "C" + "A" * 49 + "T"          # offsets 49 and 99 are marked
    got = gc.read_bases_at(1000, "100M", seq, "I" * 100, {1000, 1049, 1099})

    assert got == {1000: "A", 1049: "C", 1099: "T"}


def test_read_bases_at_leading_soft_clip_consumes_read_not_reference():
    """POS points at the first ALIGNED base, so the clipped bases sit before it in SEQ.

    Treating S as consuming the reference would return the clipped G here.
    """
    seq = "G" * 10 + "A" * 89 + "T"
    got = gc.read_bases_at(1000, "10S90M", seq, "I" * 100, {1000, 1089})

    assert got == {1000: "A", 1089: "T"}


def test_read_bases_at_trailing_soft_clip_is_not_aligned():
    seq = "A" * 90 + "G" * 10
    got = gc.read_bases_at(1000, "90M10S", seq, "I" * 100, {1000, 1089, 1090, 1095})

    assert got == {1000: "A", 1089: "A"}          # 1090+ is past the alignment


def test_read_bases_at_insertion_shifts_the_read_but_not_the_reference():
    """Every base after an insertion moves in SEQ while the reference stands still.

    Reference 1050 is the first C. Ignoring the 5I would read a G, the inserted base.
    """
    seq = "A" * 50 + "G" * 5 + "C" * 50
    got = gc.read_bases_at(1000, "50M5I50M", seq, "I" * 105, {1049, 1050, 1099})

    assert got == {1049: "A", 1050: "C", 1099: "C"}
    assert "G" not in got.values()


def test_read_bases_at_deletion_positions_are_absent_not_guessed():
    """The read has no base inside a deletion, so those positions are omitted. Returning
    the neighbouring base instead would invent an allele call at a diagnostic site."""
    seq = "A" * 50 + "C" * 50
    got = gc.read_bases_at(1000, "50M5D50M", seq, "I" * 100,
                           {1049, 1050, 1052, 1054, 1055, 1104})

    assert got == {1049: "A", 1055: "C", 1104: "C"}
    assert 1050 not in got and 1052 not in got and 1054 not in got


def test_read_bases_at_skipped_reference_behaves_like_a_deletion():
    seq = "A" * 50 + "C" * 50
    deletion = gc.read_bases_at(1000, "50M5D50M", seq, "I" * 100, set(range(1000, 1105)))
    skipped = gc.read_bases_at(1000, "50M5N50M", seq, "I" * 100, set(range(1000, 1105)))

    assert skipped == deletion
    assert set(range(1050, 1055)).isdisjoint(skipped)


def test_read_bases_at_equal_and_mismatch_operations_behave_like_match():
    seq = "A" * 50 + "C" * 50
    explicit = gc.read_bases_at(1000, "50=50X", seq, "I" * 100, {1000, 1049, 1050, 1099})
    plain = gc.read_bases_at(1000, "100M", seq, "I" * 100, {1000, 1049, 1050, 1099})

    assert explicit == plain == {1000: "A", 1049: "A", 1050: "C", 1099: "C"}


def test_read_bases_at_hard_clip_consumes_neither():
    """Hard-clipped bases are absent from SEQ, so the first aligned base is offset 0.

    Treating H like S would skip 10 real bases and return an A here.
    """
    seq = "T" + "A" * 98 + "C"
    got = gc.read_bases_at(1000, "10H100M", seq, "I" * 100, {1000, 1099})

    assert got == {1000: "T", 1099: "C"}


def test_read_bases_at_combines_clip_insertion_and_deletion():
    """One read with everything on it, as a real spanning read tends to look."""
    seq = "G" * 5 + "A" * 20 + "T" * 3 + "C" * 20      # 5S 20M 3I 20M, with a 4D in between
    got = gc.read_bases_at(1000, "5S20M3I10M4D10M", seq, "I" * len(seq),
                           {1000, 1019, 1020, 1029, 1030, 1033, 1034, 1043})

    assert got == {1000: "A", 1019: "A", 1020: "C", 1029: "C", 1034: "C", 1043: "C"}
    assert 1030 not in got and 1033 not in got       # the 4D window
    assert "T" not in got.values()                    # the insertion is not on the reference


def test_read_bases_at_drops_bases_below_the_quality_floor():
    got = gc.read_bases_at(1000, "3M", "ACG", HIGH_BQ + LOW_BQ + HIGH_BQ, {1000, 1001, 1002})

    assert got == {1000: "A", 1002: "G"}


def test_read_bases_at_quality_floor_is_inclusive():
    at_floor = chr(33 + 13)
    below = chr(33 + 12)

    assert gc.read_bases_at(1000, "1M", "A", at_floor, {1000}) == {1000: "A"}
    assert gc.read_bases_at(1000, "1M", "A", below, {1000}) == {}
    assert gc.read_bases_at(1000, "1M", "A", below, {1000}, min_bq=12) == {1000: "A"}
    assert gc.read_bases_at(1000, "1M", "A", below, {1000}, min_bq=0) == {1000: "A"}


def test_read_bases_at_accepts_every_base_when_qualities_are_absent():
    """QUAL '*' means the record carries no qualities, which must not silence the read."""
    assert gc.read_bases_at(1000, "3M", "ACG", "*", {1000, 1001, 1002}) == {1000: "A", 1001: "C",
                                                                           1002: "G"}
    assert gc.read_bases_at(1000, "3M", "ACG", "", {1000, 1001, 1002}) == {1000: "A", 1001: "C",
                                                                          1002: "G"}


def test_read_bases_at_ignores_positions_the_read_does_not_reach():
    got = gc.read_bases_at(1000, "100M", "A" * 100, "I" * 100, {999, 1050, 1100, 5000})

    assert got == {1050: "A"}


def test_read_bases_at_stops_at_the_end_of_seq():
    """A CIGAR claiming more aligned bases than SEQ holds must not run off the sequence."""
    got = gc.read_bases_at(1000, "100M", "ACGTA", "*", {1000, 1004, 1005, 1099})

    assert got == {1000: "A", 1004: "A"}


def test_read_bases_at_uppercases_the_base():
    assert gc.read_bases_at(1000, "2M", "ac", "*", {1000, 1001}) == {1000: "A", 1001: "C"}


# ------------------------------------------------------------------- collect_read_alleles

SITES = {100: ("A", "G"), 110: ("A", "G"), 120: ("T", "C")}


def test_collect_read_alleles_labels_each_base_by_the_copy_it_belongs_to():
    lines = [
        _sam_line("acc", 100, "21M", "A" + "N" * 9 + "A" + "N" * 9 + "T"),
        _sam_line("don", 100, "21M", "G" + "N" * 9 + "G" + "N" * 9 + "C"),
        _sam_line("neither", 100, "21M", "C" + "N" * 9 + "C" + "N" * 9 + "A"),
    ]

    per_read = gc.collect_read_alleles(lines, SITES)

    assert per_read["acc"] == {100: "acceptor", 110: "acceptor", 120: "acceptor"}
    assert per_read["don"] == {100: "donor", 110: "donor", 120: "donor"}
    # A third base is neither copy's allele: sequencing error or a real third state, but not
    # evidence about conversion either way, so it is kept out of the allele fraction.
    assert per_read["neither"] == {100: "other", 110: "other", 120: "other"}


def test_collect_read_alleles_treats_two_mates_as_one_molecule():
    """Keyed by read NAME: a pair is one physical molecule, and counting it twice would
    double the apparent support for whatever it carries."""
    lines = [
        _sam_line("frag", 100, "1M", "G", flag=99),
        _sam_line("frag", 120, "1M", "T", flag=147),
    ]

    per_read = gc.collect_read_alleles(lines, SITES)

    assert list(per_read) == ["frag"]
    assert per_read["frag"] == {100: "donor", 120: "acceptor"}


@pytest.mark.parametrize("flag,what", [
    (0x4, "unmapped"),
    (0x100, "secondary"),
    (0x800, "supplementary"),
    (0x904, "all three at once"),
])
def test_collect_read_alleles_skips_unmapped_secondary_and_supplementary(flag, what):
    """A secondary or supplementary record is another copy of a molecule already counted."""
    lines = [_sam_line("r", 100, "1M", "G", flag=flag)]

    assert gc.collect_read_alleles(lines, SITES) == {}


def test_collect_read_alleles_keeps_reverse_strand_and_duplicate_records():
    """0x10 and 0x400 are not in the mask, so those reads still count."""
    lines = [_sam_line("rev", 100, "1M", "G", flag=16),
             _sam_line("dup", 100, "1M", "G", flag=1024)]

    assert set(gc.collect_read_alleles(lines, SITES)) == {"rev", "dup"}


def test_collect_read_alleles_skips_headers_blank_and_malformed_lines():
    lines = ["@HD\tVN:1.6", "@SQ\tSN:chr\tLN:2300", "", "truncated\t0\tchr",
             _sam_line("r", 100, "1M", "G")]

    assert gc.collect_read_alleles(lines, SITES) == {"r": {100: "donor"}}


def test_collect_read_alleles_skips_records_without_an_alignment():
    lines = [_sam_line("nocigar", 100, "*", "G"), _sam_line("nopos", 0, "1M", "G")]

    assert gc.collect_read_alleles(lines, SITES) == {}


def test_collect_read_alleles_applies_the_default_base_quality_floor():
    """Observed behaviour worth pinning: the floor is `read_bases_at`'s default of 13.

    `collect_read_alleles` takes no min_bq, and main() never forwards its own --min-bq, so
    that option currently has no effect. See the report accompanying these tests.
    """
    seq = "G" + "N" * 9 + "G" + "N" * 9 + "C"           # donor allele at all three sites
    qual = HIGH_BQ * 10 + LOW_BQ + HIGH_BQ * 10         # but site 110 is a bad base call
    lines = [_sam_line("r", 100, "21M", seq, qual=qual)]

    per_read = gc.collect_read_alleles(lines, SITES)

    assert per_read["r"] == {100: "donor", 120: "donor"}
    assert 110 not in per_read["r"]


def test_collect_read_alleles_reads_only_the_diagnostic_positions():
    lines = [_sam_line("r", 90, "40M", "G" * 40)]

    assert set(gc.collect_read_alleles(lines, SITES)["r"]) == set(SITES)


# ----------------------------------------------------------------------------- pileup


def test_pileup_computes_the_donor_fraction_over_informative_reads_only():
    """`other` calls count towards depth but not towards the fraction: they say nothing
    about which copy the molecule came from."""
    per_read = _per_read({
        "d1": {100: "donor"}, "d2": {100: "donor"}, "d3": {100: "donor"},
        "a1": {100: "acceptor"},
        "o1": {100: "other"}, "o2": {100: "other"},
    })

    counts = gc.pileup(per_read, [100])

    assert counts[100]["donor"] == 3
    assert counts[100]["acceptor"] == 1
    assert counts[100]["other"] == 2
    assert counts[100]["depth"] == 6
    assert counts[100]["donor_af"] == pytest.approx(0.75)


def test_pileup_leaves_the_fraction_undefined_without_informative_reads():
    """None, not 0.0: no informative read is not the same as no donor allele, and the
    tract caller must not read the second out of the first."""
    per_read = _per_read({"o1": {100: "other"}, "o2": {100: "other"}})

    counts = gc.pileup(per_read, [100, 200])

    assert counts[100]["donor_af"] is None and counts[100]["depth"] == 2
    assert counts[200]["donor_af"] is None and counts[200]["depth"] == 0


def test_pileup_ignores_calls_outside_the_requested_positions():
    per_read = _per_read({"r": {100: "donor", 999: "donor"}})

    counts = gc.pileup(per_read, [100])

    assert set(counts) == {100}
    assert counts[100]["donor"] == 1


# ------------------------------------------------------------------------- call_tracts


def _counts(spec, depth=10):
    """Build a pileup-shaped table straight from {position: donor_af}."""
    out = {}
    for pos, af in spec.items():
        d = depth[pos] if isinstance(depth, dict) else depth
        donor = 0 if af is None else round(af * d)
        out[pos] = {"acceptor": d - donor, "donor": donor, "other": 0,
                    "depth": d, "donor_af": af}
    return out


def test_call_tracts_finds_a_maximal_run_of_donor_sites():
    positions = [10, 20, 30, 40, 50, 60]
    counts = _counts({10: 0.0, 20: 1.0, 30: 0.95, 40: 1.0, 50: 0.0, 60: 0.0})

    assert gc.call_tracts(positions, counts) == [[20, 30, 40]]


def test_call_tracts_means_consecutive_in_the_site_list_not_in_coordinates():
    """The bases between two diagnostic sites are identical in both copies, so they cannot
    testify either way. A 4 kb gap between adjacent diagnostic sites does not break a run."""
    positions = [100, 200, 5000, 5001]
    counts = _counts(dict.fromkeys(positions, 1.0))

    assert gc.call_tracts(positions, counts) == [[100, 200, 5000, 5001]]


def test_call_tracts_finds_several_runs_and_flushes_the_last_one():
    positions = [10, 20, 30, 40, 50, 60, 70]
    counts = _counts({10: 1.0, 20: 1.0, 30: 1.0, 40: 0.0, 50: 1.0, 60: 1.0, 70: 1.0})

    assert gc.call_tracts(positions, counts) == [[10, 20, 30], [50, 60, 70]]


def test_call_tracts_requires_min_sites():
    positions = [10, 20, 30]
    counts = _counts({10: 1.0, 20: 1.0, 30: 0.0})

    assert gc.call_tracts(positions, counts, min_sites=3) == []
    assert gc.call_tracts(positions, counts, min_sites=2) == [[10, 20]]


def test_call_tracts_requires_min_af():
    positions = [10, 20, 30, 40]
    counts = _counts({10: 1.0, 20: 0.6, 30: 1.0, 40: 1.0})

    assert gc.call_tracts(positions, counts, min_af=0.7, min_sites=2) == [[30, 40]]
    assert gc.call_tracts(positions, counts, min_af=0.6, min_sites=2) == [[10, 20, 30, 40]]
    assert gc.call_tracts(positions, counts, min_af=0.61, min_sites=2) == [[30, 40]]


def test_a_shallow_site_bridges_a_tract_instead_of_breaking_it():
    """A donor fraction computed from two reads is not evidence of anything, in EITHER direction.

    It used to break the run, which meant one ordinary coverage dip inside a clean tract split it
    into fragments that each fell below min_sites, and the whole tract vanished. Found on a
    simulated cohort: a real four-site conversion went undetected because one site sat at depth 2.
    The shallow site is skipped, so it neither joins the tract nor ends it.
    """
    positions = [10, 20, 30, 40]
    counts = _counts(dict.fromkeys(positions, 1.0), depth={10: 10, 20: 10, 30: 2, 40: 10})

    # The shallow site is not a member, but the sites either side still form one tract.
    assert gc.call_tracts(positions, counts, min_sites=2, min_depth=5) == [[10, 20, 40]]
    # With the floor low enough to genotype it, it joins like any other site.
    assert gc.call_tracts(positions, counts, min_sites=2, min_depth=2) == [[10, 20, 30, 40]]


def test_a_site_carrying_the_acceptor_allele_does_end_a_tract():
    """The counterpart: a MEASURED acceptor site is evidence, and it bounds the tract."""
    positions = [10, 20, 30, 40]
    counts = _counts({10: 1.0, 20: 1.0, 30: 0.0, 40: 1.0}, depth=20)

    assert gc.call_tracts(positions, counts, min_sites=2, min_depth=5) == [[10, 20]]


def test_a_site_with_no_informative_reads_also_bridges():
    """No informative read at all is the same situation as too few: undetermined, not acceptor."""
    positions = [10, 20, 30, 40, 50]
    counts = _counts({10: 1.0, 20: 1.0, 30: None, 40: 1.0, 50: 1.0})

    assert gc.call_tracts(positions, counts, min_sites=2, min_depth=1) == [[10, 20, 40, 50]]


# ---------------------------------------------------------------------- tract_evidence


def test_tract_evidence_summarises_a_bounded_tract():
    positions = [10, 20, 30, 40, 50]
    counts = _counts({10: 0.0, 20: 1.0, 30: 0.9, 40: 1.0, 50: 0.0}, depth=8)
    per_read = _per_read({"r": {20: "donor", 30: "donor"}})

    ev = gc.tract_evidence([20, 30, 40], positions, counts, per_read)

    assert ev["n_sites"] == 3
    assert (ev["start"], ev["end"], ev["span_bp"]) == (20, 40, 21)
    assert ev["n_sites_outside"] == 2
    assert ev["donor_af_in"] == pytest.approx(0.9667, abs=1e-4)
    assert ev["donor_af_outside"] == 0.0
    assert ev["min_depth"] == 8


def test_tract_evidence_counts_cis_breakpoint_and_donor_only_reads():
    """The three read-level counters, one molecule of each kind.

    Only `breakpoint` is unfakeable by a mismapping: donor alleles inside the tract and
    acceptor alleles outside it, on the same physical fragment.
    """
    positions = [10, 20, 30, 40, 50]
    counts = _counts(dict.fromkeys(positions, 0.5), depth=6)
    per_read = _per_read({
        "breakpoint": {20: "donor", 30: "donor", 50: "acceptor"},
        "cis_only": {20: "donor", 40: "donor"},
        "donor_only": {20: "donor", 30: "donor", 50: "donor"},
        "one_site": {20: "donor"},
        "acceptor_read": {20: "acceptor", 50: "acceptor"},
    })

    ev = gc.tract_evidence([20, 30, 40], positions, counts, per_read)

    assert ev["cis_reads"] == 3                 # breakpoint, cis_only, donor_only
    assert ev["breakpoint_reads"] == 1
    assert ev["donor_only_reads"] == 1


def test_tract_evidence_does_not_count_a_donor_only_read_that_also_reads_back():
    """A single acceptor call outside the tract disqualifies a read from being donor-only:
    a molecule that came wholesale from the donor never crosses back."""
    positions = [10, 20, 30]
    counts = _counts(dict.fromkeys(positions, 0.5), depth=6)
    per_read = _per_read({"r": {20: "donor", 30: "donor", 10: "acceptor"}})

    ev = gc.tract_evidence([20], positions, counts, per_read)

    assert ev["donor_only_reads"] == 0
    assert ev["breakpoint_reads"] == 1


def test_tract_evidence_averages_the_donor_fraction_outside_the_tract():
    """donor_af_outside is the mean over sites outside, and it is the direct test for
    mismapping: a conversion is bounded, so it leaves those sites alone."""
    positions = [10, 20, 30, 40]
    counts = _counts({10: 0.4, 20: 1.0, 30: 1.0, 40: 0.2})

    ev = gc.tract_evidence([20, 30], positions, counts, {})

    assert ev["donor_af_outside"] == pytest.approx(0.3)
    assert ev["n_sites_outside"] == 2


def test_an_undetermined_site_counts_as_outside_in_neither_direction():
    """A site too shallow to genotype must not land in the "outside" average either.

    Excluding it from the tract but still averaging its fraction into donor_af_outside lets a
    number measured off a couple of reads push the mean past the mismapping threshold and flip a
    real conversion. Found on a simulated cohort, where exactly that turned a true positive into
    a mismapping call.
    """
    positions = [10, 20, 30, 40]
    counts = _counts({10: 0.4, 20: 1.0, 30: 1.0, 40: None})

    ev = gc.tract_evidence([20, 30], positions, counts, {}, min_depth=5)

    assert ev["donor_af_outside"] == pytest.approx(0.4)
    assert ev["n_sites_outside"] == 1, "the undetermined site must not be counted as outside"


def test_a_shallow_site_does_not_drag_the_outside_average():
    """The regression in full: one shallow donor-looking site outside a clean bounded tract."""
    positions = [10, 20, 30, 40, 50]
    counts = _counts({10: 0.0, 20: 1.0, 30: 1.0, 40: 1.0, 50: 0.0},
                     depth={10: 30, 20: 30, 30: 30, 40: 2, 50: 30})

    ev = gc.tract_evidence([20, 30], positions, counts, {}, min_depth=5)

    assert ev["donor_af_outside"] == pytest.approx(0.0)
    assert gc.classify(ev)[0] != "mismapping"


def test_two_tracts_in_one_locus_shadow_each_other():
    """Known limitation, pinned rather than asserted as desirable.

    `outside` is every diagnostic site of the locus that is not in THIS tract, including
    the sites of a second tract. Two genuine conversion tracts in the same paralog pair
    therefore inflate each other's `donor_af_outside` and both come out as mismapping,
    even with breakpoint reads, because the outside test runs first. Reported alongside
    these tests; the conservative direction, but a false negative.
    """
    positions = [10, 20, 30, 40, 50, 60, 70, 80]
    counts = _counts({10: 0.0, 20: 1.0, 30: 1.0, 40: 1.0,
                      50: 0.0, 60: 1.0, 70: 1.0, 80: 1.0})
    per_read = _per_read({"spans": {10: "acceptor", 20: "donor", 30: "donor"}})

    first = gc.tract_evidence([20, 30, 40], positions, counts, per_read)

    assert first["breakpoint_reads"] == 1
    assert first["donor_af_outside"] == pytest.approx(0.6)     # the second tract, mostly
    assert gc.classify(first)[0] == "mismapping"


def test_tract_evidence_reports_no_outside_fraction_for_a_whole_locus_tract():
    positions = [10, 20, 30]
    counts = _counts(dict.fromkeys(positions, 1.0))

    ev = gc.tract_evidence(positions, positions, counts, {})

    assert ev["n_sites_outside"] == 0
    assert ev["donor_af_outside"] is None


def test_tract_evidence_min_depth_is_the_worst_site_inside_the_tract():
    positions = [10, 20, 30]
    counts = _counts(dict.fromkeys(positions, 1.0), depth={10: 40, 20: 7, 30: 30})

    assert gc.tract_evidence(positions, positions, counts, {})["min_depth"] == 7


# --------------------------------------------------------------------------- classify
#
# The four routes through classify are the specification of the tool. One test each.


def test_classify_calls_a_bounded_near_fixed_tract_with_a_spanning_read_a_conversion():
    """Case 1, a real conversion: donor alleles near-fixed inside a bounded tract, acceptor
    alleles outside, and at least one molecule carrying donor INSIDE and acceptor OUTSIDE.

    That last one is the piece a mismapping cannot fake, so it decides the verdict and the
    reason says so.
    """
    verdict, reason = gc.classify(_ev(donor_af_outside=0.02, donor_af_in=0.97,
                                      n_sites_outside=6, cis_reads=9, breakpoint_reads=4))

    assert verdict == "gene_conversion"
    assert reason == "4 read(s) cross a breakpoint in cis"


def test_classify_calls_donor_alleles_outside_the_tract_mismapping():
    """Case 2: mismapped donor reads carry donor alleles at EVERY diagnostic site of the
    locus, including outside the candidate tract. A conversion is bounded, so it does not."""
    verdict, reason = gc.classify(_ev(donor_af_outside=0.45, donor_only_reads=30,
                                      breakpoint_reads=0, cis_reads=30))

    assert verdict == "mismapping"
    assert "outside the tract" in reason


def test_classify_puts_the_outside_test_before_the_breakpoint_test():
    """Observed precedence: donor alleles all over the locus outrank breakpoint reads, so a
    handful of chimeric or misassembled reads cannot rescue a wholesale mismapping."""
    assert gc.classify(_ev(donor_af_outside=0.45, breakpoint_reads=5))[0] == "mismapping"


def test_classify_mismapping_threshold_is_strict():
    assert gc.classify(_ev(donor_af_outside=0.25, breakpoint_reads=0, cis_reads=5))[0] != "mismapping"
    assert gc.classify(_ev(donor_af_outside=0.2501, breakpoint_reads=0))[0] == "mismapping"
    assert gc.classify(_ev(donor_af_outside=0.6, breakpoint_reads=0),
                       max_outside_af=0.7)[0] != "mismapping"


def test_classify_calls_an_unbounded_whole_locus_tract_ambiguous():
    """Case 3, and a regression test: the tract covers EVERY diagnostic site of the locus,
    so `n_sites_outside == 0`, and no read spans a breakpoint.

    This used to be called `gene_conversion`, which is a false positive. With no site
    outside the tract there is nothing the donor alleles are bounded BY: a locus wholly
    replaced by its paralog and a locus whose reads all arrived from its paralog are
    identical from site data alone. `donor_af_outside` is None here, not low, and absence
    of contrary evidence is not evidence.
    """
    ev = _ev(n_sites_outside=0, donor_af_outside=None, donor_af_in=1.0,
             cis_reads=25, breakpoint_reads=0, donor_only_reads=25)

    verdict, reason = gc.classify(ev)

    assert verdict == "ambiguous"
    assert "not bounded" in reason
    assert "mismapping" in reason


def test_classify_still_believes_a_whole_locus_tract_with_a_breakpoint_read():
    """The one thing that can rescue case 3: a molecule that reads back to acceptor alleles
    at a site the tract does not cover. It cannot come from a donor read."""
    ev = _ev(n_sites_outside=0, donor_af_outside=None, breakpoint_reads=1)

    assert gc.classify(ev)[0] == "gene_conversion"


def test_classify_accepts_a_bounded_tract_without_a_breakpoint_read_only_in_cis():
    """Case 4: bounded and near-fixed, but no molecule spans a breakpoint.

    Then the fallback is read-level support in cis. With none, the tract is still reported,
    never hidden, but it is called ambiguous and the reason names the test it failed.
    """
    bounded = dict(breakpoint_reads=0, n_sites_outside=4, donor_af_in=0.95, donor_af_outside=0.01)

    assert gc.classify(_ev(**bounded, cis_reads=3)) == (
        "gene_conversion", "donor allele near-fixed within a bounded tract")

    verdict, reason = gc.classify(_ev(**bounded, cis_reads=0))
    assert verdict == "ambiguous"
    assert "no read spans a breakpoint" in reason


def test_classify_requires_the_donor_allele_to_be_near_fixed_without_a_breakpoint_read():
    """A fraction sitting in between is what an aligner splitting reads between two copies
    produces, not what a clonal sample's converted locus looks like."""
    middling = _ev(breakpoint_reads=0, n_sites_outside=4, donor_af_in=0.6,
                   donor_af_outside=0.1, cis_reads=12)

    assert gc.classify(middling)[0] == "ambiguous"
    assert gc.classify(middling, min_af_in=0.5)[0] == "gene_conversion"


def test_classify_handles_missing_fractions():
    """Neither None can crash the classifier or be silently read as a number."""
    assert gc.classify(_ev(donor_af_outside=None, breakpoint_reads=2))[0] == "gene_conversion"
    assert gc.classify(_ev(donor_af_in=None, breakpoint_reads=0, cis_reads=5))[0] == "ambiguous"


# ------------------------------------------------------------------------- load_sites


def _sites_tsv(path, loci, acc="A", don="G"):
    """A sites.tsv in paralog_map.py's format: {pair_id: [acceptor positions]}."""
    lines = ["\t".join(["pair_id", "acceptor", "acc_pos", "acc_base",
                        "donor", "don_pos", "don_base", "strand"])]
    for pair_id, positions in loci.items():
        for p in positions:
            lines.append("\t".join([str(pair_id), CONTIG, str(p), acc,
                                    CONTIG, str(p + 1200), don, "+"]))
    path.write_text("\n".join(lines) + "\n")
    return str(path)


def test_load_sites_groups_by_locus(tmp_path):
    path = _sites_tsv(tmp_path / "sites.tsv", {0: [100, 110], 1: [1500]})

    loci = gc.load_sites(path)

    assert set(loci) == {"0", "1"}
    assert loci["0"]["sites"] == {100: ("A", "G"), 110: ("A", "G")}
    assert loci["0"]["contig"] == CONTIG and loci["0"]["donor"] == CONTIG
    assert loci["1"]["sites"] == {1500: ("A", "G")}


def test_load_sites_uppercases_bases_and_reads_columns_by_name(tmp_path):
    """Column order is read from the header, not assumed."""
    path = tmp_path / "sites.tsv"
    path.write_text("acc_pos\tacc_base\tdon_base\tdonor\tacceptor\tpair_id\n"
                    "100\ta\tg\tchrB\tchrA\t7\n")

    loci = gc.load_sites(str(path))

    assert loci["7"]["sites"] == {100: ("A", "G")}
    assert loci["7"]["contig"] == "chrA" and loci["7"]["donor"] == "chrB"


def test_load_sites_skips_truncated_rows(tmp_path):
    path = tmp_path / "sites.tsv"
    _sites_tsv(path, {0: [100, 110]})
    path.write_text(path.read_text() + "0\tchr\t120\n")

    assert set(gc.load_sites(str(path))["0"]["sites"]) == {100, 110}


# ------------------------------------------------------------------------------- main


def _bam(tmp_path, records, contig_length=2300):
    """Sort a hand-written SAM into an indexed BAM with samtools."""
    sam = tmp_path / "aln.sam"
    sam.write_text("@HD\tVN:1.6\tSO:unsorted\n"
                   f"@SQ\tSN:{CONTIG}\tLN:{contig_length}\n"
                   + "\n".join(records) + "\n")
    bam = tmp_path / "aln.bam"
    subprocess.run(["samtools", "sort", "-o", str(bam), str(sam)], check=True,
                   capture_output=True)
    subprocess.run(["samtools", "index", str(bam)], check=True, capture_output=True)
    return str(bam)


def _read_at(name, pos, alleles, length=60):
    """A read whose bases at the given reference positions are set, everything else filler."""
    seq = list("N" * length)
    for p, base in alleles.items():
        seq[p - pos] = base
    return _sam_line(name, pos, f"{length}M", "".join(seq), qual=HIGH_BQ * length)


@pytest.fixture
def samtools():
    exe = shutil.which("samtools")
    if not exe:
        pytest.skip("samtools is not on PATH")
    return exe


def test_main_end_to_end(tmp_path, samtools, capsys):
    """Three loci through the real CLI: one converted, one clean, one unbounded.

    Every read is MAPQ 0, which is what an aligner gives a paralogous region. If the tool
    ever grows a MAPQ filter, this test goes to zero rows.
    """
    converted = [100, 110, 120, 130, 140, 150]          # tract over the middle four
    clean = [1500, 1510, 1520, 1530]
    whole_locus = [2000, 2010, 2020, 2030]

    records = []
    for i in range(6):
        # Donor bases at 110-140, the locus's own alleles at 100 and 150: each of these
        # molecules crosses both breakpoints.
        alleles = {p: ("G" if 110 <= p <= 140 else "A") for p in converted}
        records.append(_read_at(f"conv{i}", 100, alleles))
    for i in range(6):
        records.append(_read_at(f"clean{i}", 1500, dict.fromkeys(clean, "A")))
    for i in range(8):
        records.append(_read_at(f"whole{i}", 2000, dict.fromkeys(whole_locus, "G")))

    bam = _bam(tmp_path, records)
    sites = _sites_tsv(tmp_path / "sites.tsv",
                       {0: converted, 1: clean, 2: whole_locus})
    out = tmp_path / "tracts.tsv"

    assert gc.main(["--sites", sites, "--bam", bam, "--sample", "S1", "-o", str(out),
                    "--samtools", samtools]) == 0

    lines = out.read_text().rstrip("\n").split("\n")
    header = lines[0].split("\t")
    rows = [dict(zip(header, line.split("\t"))) for line in lines[1:]]

    assert header == gc.COLUMNS
    # The clean locus produces no candidate tract at all, so it produces no row.
    assert [r["pair_id"] for r in rows] == ["0", "2"]

    converted_row = rows[0]
    assert converted_row["sample"] == "S1"
    assert converted_row["contig"] == CONTIG and converted_row["donor"] == CONTIG
    assert converted_row["verdict"] == "gene_conversion"
    assert "cross a breakpoint in cis" in converted_row["reason"]
    assert (converted_row["start"], converted_row["end"], converted_row["span_bp"]) == \
        ("110", "140", "31")
    assert (converted_row["n_sites"], converted_row["n_sites_outside"]) == ("4", "2")
    assert float(converted_row["donor_af_in"]) == 1.0
    assert float(converted_row["donor_af_outside"]) == 0.0
    assert converted_row["min_depth"] == "6"
    assert (converted_row["cis_reads"], converted_row["breakpoint_reads"]) == ("6", "6")
    assert converted_row["donor_only_reads"] == "0"

    unbounded_row = rows[1]
    assert unbounded_row["verdict"] == "ambiguous"
    assert "not bounded" in unbounded_row["reason"]
    assert unbounded_row["n_sites_outside"] == "0"
    assert unbounded_row["donor_af_outside"] == ""      # None is written as an empty cell
    # With no site outside the tract, neither counter that discriminates can fire: no read
    # can cross a breakpoint, and no read can look donor-only either. That symmetry is the
    # whole reason the verdict has to be ambiguous.
    assert unbounded_row["cis_reads"] == "8"
    assert (unbounded_row["breakpoint_reads"], unbounded_row["donor_only_reads"]) == ("0", "0")

    err = capsys.readouterr().err
    assert "S1: 2 candidate tract(s), 1 called as conversion" in err


def test_main_skips_a_locus_with_too_few_diagnostic_sites(tmp_path, samtools):
    """Fewer diagnostic sites than min_sites means no tract is reachable, so the locus is
    dropped before samtools is invoked at all."""
    records = [_read_at(f"r{i}", 100, {100: "G", 110: "G"}, length=30) for i in range(6)]
    bam = _bam(tmp_path, records)
    sites = _sites_tsv(tmp_path / "sites.tsv", {0: [100, 110]})
    out = tmp_path / "tracts.tsv"

    assert gc.main(["--sites", sites, "--bam", bam, "--sample", "S1", "-o", str(out),
                    "--samtools", samtools]) == 0

    assert out.read_text() == "\t".join(gc.COLUMNS) + "\n"


def test_main_thresholds_reach_the_tract_caller(tmp_path, samtools):
    """Four reads is under the default informative depth of 5, so nothing is called until
    --min-depth is lowered."""
    positions = [100, 110, 120]
    records = [_read_at(f"r{i}", 100, dict.fromkeys(positions, "G"), length=30) for i in range(4)]
    bam = _bam(tmp_path, records)
    sites = _sites_tsv(tmp_path / "sites.tsv", {0: positions})
    out = tmp_path / "tracts.tsv"
    argv = ["--sites", sites, "--bam", bam, "--sample", "S1", "-o", str(out),
            "--samtools", samtools]

    assert gc.main(argv) == 0
    assert len(out.read_text().splitlines()) == 1

    assert gc.main(argv + ["--min-depth", "4"]) == 0
    assert len(out.read_text().splitlines()) == 2


def test_main_fails_loudly_when_samtools_fails(tmp_path, samtools):
    sites = _sites_tsv(tmp_path / "sites.tsv", {0: [100, 110, 120]})

    with pytest.raises(RuntimeError, match="samtools view failed"):
        gc.main(["--sites", sites, "--bam", str(tmp_path / "missing.bam"), "--sample", "S1",
                 "-o", str(tmp_path / "out.tsv"), "--samtools", samtools])


def test_script_runs_as_a_subprocess(tmp_path, samtools, repo_root):
    """The pipeline invokes it as an executable, not as an import."""
    records = [_read_at(f"r{i}", 100, {100: "G", 110: "G", 120: "G"}, length=30) for i in range(6)]
    bam = _bam(tmp_path, records)
    sites = _sites_tsv(tmp_path / "sites.tsv", {0: [100, 110, 120]})
    out = tmp_path / "tracts.tsv"

    res = subprocess.run(["python3", str(repo_root / "bin" / "gene_conversion.py"),
                          "--sites", sites, "--bam", bam, "--sample", "S1", "-o", str(out)],
                         capture_output=True, text=True)

    assert res.returncode == 0, res.stderr
    assert "1 candidate tract(s)" in res.stderr
    assert out.read_text().splitlines()[1].startswith("S1\t0\t")


# --------------------------------------------------------------------------- #
# Regressions: two defects the first version of these tools shipped with
# --------------------------------------------------------------------------- #

def _flat_counts(positions, donor_positions, depth=30):
    """Per-site counts where the given positions are pure donor and the rest pure acceptor."""
    counts = {}
    for p in positions:
        donor = p in donor_positions
        counts[p] = {"acceptor": 0 if donor else depth, "donor": depth if donor else 0,
                     "other": 0, "depth": depth, "donor_af": 1.0 if donor else 0.0}
    return counts


def test_two_tracts_in_one_locus_do_not_shadow_each_other():
    """A second genuine tract is not "outside" the first one.

    Regression: `outside` used to be every site not in THIS tract, so a locus with two real
    conversions had each one counting the other's donor alleles as evidence against itself.
    Because the outside test runs before the breakpoint test, both came out as mismapping. Two
    conversions in one locus is an ordinary outcome, not a pathological input.
    """
    positions = list(range(1, 13))
    counts = _flat_counts(positions, donor_positions=set(range(1, 5)) | set(range(9, 13)))
    per_read = {"r1": {1: "donor", 2: "donor", 5: "acceptor"}}

    tracts = gc.call_tracts(positions, counts, 0.7, 3, 5)
    assert [(min(t), max(t)) for t in tracts] == [(1, 4), (9, 12)]

    in_any = {p for t in tracts for p in t}
    for tract in tracts:
        ev = gc.tract_evidence(tract, positions, counts, per_read, in_any)
        assert ev["donor_af_outside"] == 0.0, "the other tract was counted as outside"
        assert gc.classify(ev)[0] != "mismapping"


def test_the_base_quality_floor_reaches_the_base_reader():
    """Regression: --min-bq was parsed and then never forwarded, so it did nothing at all."""
    sites = {10: ("A", "G")}
    # One read whose base at position 10 is Phred 30 ('?').
    sam = "\t".join(["r1", "0", "chr", "10", "0", "1M", "*", "0", "0", "G", "?"])

    permissive = gc.collect_read_alleles([sam], sites, min_bq=13)
    strict = gc.collect_read_alleles([sam], sites, min_bq=40)

    assert permissive["r1"] == {10: "donor"}
    assert strict == {}, "a base below the floor must not be read"


def test_a_quality_string_shorter_than_the_sequence_does_not_crash():
    """Malformed SAM should skip the base, not raise IndexError out of a CIGAR walk."""
    got = gc.read_bases_at(10, "5M", "ACGTA", "II", {10, 11, 12, 13, 14}, min_bq=13)
    assert got == {10: "A", 11: "C"}, "positions past the end of QUAL must be dropped, not crash"
