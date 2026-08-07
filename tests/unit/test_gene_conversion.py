"""Unit tests for bin/gene_conversion.py (gene conversion tract detection).

Two things are being tested, and they fail in different ways.

`read_bases_at` is arithmetic: it walks a CIGAR to line reference positions up with read
offsets. When it is wrong it is silently wrong, reading the right position off the wrong
base, so it is tested operation by operation against reads whose bases are laid out to make
an off-by-anything visible.

The rest is plumbing around a model that lives in `gconv_model` and is specified in
tests/unit/test_gconv_model.py. What is tested here is that the reads reach it intact and that
what it decides reaches the output file intact: the CIGAR arithmetic above it, the descriptive
statistics reported beside its verdict, and the coverage check that catches the tracts whose
reads have moved to the donor, which is the one thing no allele model can see.
"""

from __future__ import annotations

import math
import re
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


def test_a_deletion_site_is_read_as_presence_or_absence():
    """Where the donor has no base, the question is whether the read has one, not which.

    Inverting that reads every molecule backwards at exactly the sites the deletion markers were
    added for, and it does so without changing anything else: the columns still fill in, the
    verdict still comes out, and the tract is supported by the reads that contradict it.
    """
    sites = {100: ("A", gc.GAP), 110: ("A", "G")}
    lines = [
        # 13 aligned bases from 98, so 100 carries the base the donor does not have
        _sam_line("kept", 98, "13M", "NN" + "A" + "N" * 9 + "G"),
        # the same stretch with 100 deleted: no base there at all, which is the donor's state
        _sam_line("lost", 98, "2M1D10M", "NN" + "N" * 9 + "G"),
    ]

    per_read = gc.collect_read_alleles(lines, sites)

    assert per_read["kept"][100] == "acceptor"
    assert per_read["lost"][100] == "donor"
    assert per_read["lost"][110] == "donor", "the rest of the read still reads normally"


@pytest.mark.parametrize("phred,kept", [(13, True), (12, False), (40, True)])
def test_a_base_at_exactly_the_quality_floor_is_kept(phred, kept):
    """`--min-bq` is the quality a base must reach, so reaching it is enough. The same floor is
    applied in two places and only one of them was pinned."""
    lines = [_sam_line("r", 100, "1M", "G", qual=chr(33 + phred))]

    per_read = gc.collect_read_alleles(lines, SITES)

    assert (per_read.get("r", {}).get(100) == "donor") is kept


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
             # A header line is recognised by its @, not by being short. Read groups carry a tag
             # per sequencing run and go past eleven fields easily, and one that reaches the
             # record parser meets int() on a tag name rather than on a FLAG.
             "@RG\tID:1\tSM:S1\tLB:l\tPL:ILLUMINA\tPU:u\tCN:c\tDS:d\tDT:2026-01-01\tPI:300\tPM:m\tPG:p",
             _sam_line("r", 100, "1M", "G")]

    assert gc.collect_read_alleles(lines, SITES) == {"r": {100: "donor"}}


def test_collect_read_alleles_skips_records_without_an_alignment():
    """Both records reach the diagnostic sites if they are let through, which is the point.

    A record with POS 0 is one the aligner declined to place, and reading it anyway lays the
    read down from the start of the contig: the calls it produces are at real diagnostic sites
    and look like any other, which is the failure that reports a tract nobody sequenced.
    """
    lines = [_sam_line("nocigar", 100, "*", "G" * 200),
             _sam_line("nopos", 0, "200M", "G" * 200, qual="I" * 200)]

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


# ---------------------------------------------------------------------- tract_evidence


def _counts(spec, depth=10):
    """Build a pileup-shaped table straight from {position: donor_af}."""
    out = {}
    for pos, af in spec.items():
        d = depth[pos] if isinstance(depth, dict) else depth
        donor = 0 if af is None else round(af * d)
        out[pos] = {"acceptor": d - donor, "donor": donor, "other": 0,
                    "depth": d, "donor_af": af}
    return out


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


def test_a_site_at_exactly_the_depth_floor_is_determined():
    """`--min-depth` is the depth a site must reach to be genotyped, and the same floor decides
    the two sides of one split: what counts as outside, and what counts as unmeasured inside.

    Moving either comparison alone tears the split apart, and a site then counts on both sides or
    on neither, which no total in the row would reveal.
    """
    positions = [10, 20, 30, 40]
    # 10 sits at the floor and is outside; 20 sits at the floor and is inside; 30 is below it.
    counts = _counts({10: 0.0, 20: 1.0, 30: 1.0, 40: 0.0},
                     depth={10: 5, 20: 5, 30: 4, 40: 30})

    ev = gc.tract_evidence([20, 30], positions, counts, {}, min_depth=5)

    assert ev["n_sites_outside"] == 2, "a site at exactly the floor is measured"
    assert ev["n_undetermined"] == 1, "and inside the tract, only the one below it is not"


def test_n_undetermined_counts_the_unmeasured_sites_inside_the_tract():
    """A tract resting on measured sites and empty ones looks solid in the coordinates.

    The model weights each site by the reads actually there, so an empty site simply adds
    nothing. The count is reported next to the tract so a reader can see how much of its span
    was never measured at all.
    """
    positions = [10, 20, 30, 40]
    counts = _counts({10: 1.0, 20: None, 30: 1.0, 40: 0.0},
                     depth={10: 30, 20: 0, 30: 30, 40: 30})

    ev = gc.tract_evidence([10, 20, 30], positions, counts, {}, min_depth=5)

    assert ev["n_sites"] == 3
    assert ev["n_undetermined"] == 1
    assert ev["min_depth"] == 0


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


# ------------------------------------------------------------------------- load_sites


def _sites_tsv(path, loci, acc="A", don="G", kinds=None):
    """A sites.tsv in paralog_map.py's format: {pair_id: [acceptor positions]}.

    `kinds` maps an acceptor position to "del", which is how the map records a base one copy has
    and the other does not.
    """
    kinds = kinds or {}
    lines = ["\t".join(["pair_id", "acceptor", "acc_pos", "acc_base",
                        "donor", "don_pos", "don_base", "strand", "kind"])]
    for pair_id, positions in loci.items():
        for p in positions:
            kind = kinds.get(p, "snp")
            lines.append("\t".join([str(pair_id), CONTIG, str(p), acc, CONTIG, str(p + 1200),
                                    gc.GAP if kind == "del" else don, "+", kind]))
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

    lines = [ln for ln in out.read_text().rstrip("\n").split("\n") if not ln.startswith("#")]
    header = lines[0].split("\t")
    rows = [dict(zip(header, line.split("\t"))) for line in lines[1:]]

    assert header == gc.COLUMNS
    # The clean locus produces no candidate tract at all, so it produces no row.
    assert [r["pair_id"] for r in rows] == ["0", "2"]

    converted_row = rows[0]
    assert converted_row["sample"] == "S1"
    assert converted_row["contig"] == CONTIG and converted_row["donor"] == CONTIG
    assert converted_row["verdict"] == "gene_conversion"
    assert (converted_row["start"], converted_row["end"], converted_row["span_bp"]) == \
        ("110", "140", "31")
    # Both breakpoints are pinned to a single site, because every molecule crosses them.
    assert (converted_row["start_ci"], converted_row["end_ci"]) == ("110-110", "140-140")
    assert float(converted_row["log10_bf"]) > 3.0
    assert float(converted_row["post_conv"]) > 0.99
    assert (converted_row["n_sites"], converted_row["n_sites_outside"]) == ("4", "2")
    assert float(converted_row["donor_af_in"]) == 1.0
    assert float(converted_row["donor_af_outside"]) == 0.0
    assert converted_row["min_depth"] == "6"
    assert (converted_row["cis_reads"], converted_row["breakpoint_reads"]) == ("6", "6")
    assert converted_row["donor_only_reads"] == "0"

    # The third locus is entirely donor bases, with no site outside the tract to bound them.
    # A wholly converted locus and a locus whose reads all arrived from the donor have the same
    # likelihood, so the Bayes factor collapses to the ratio of their priors on its own.
    unbounded_row = rows[1]
    assert unbounded_row["verdict"] == "ambiguous"
    assert "every diagnostic site" in unbounded_row["reason"]
    assert float(unbounded_row["log10_bf"]) < 3.0
    assert float(unbounded_row["mismap_frac"]) > 0.5, "the model should reach for mismapping here"
    assert unbounded_row["n_sites_outside"] == "0"
    assert unbounded_row["donor_af_outside"] == ""      # None is written as an empty cell
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

    body = [ln for ln in out.read_text().splitlines() if not ln.startswith("#")]
    assert body == ["\t".join(gc.COLUMNS)]


def test_main_model_thresholds_reach_the_model(tmp_path, samtools):
    """The Bayes factor thresholds are what decide, and both of them are reachable.

    Four sites of six carry the donor base on every read. That is far past the default
    --min-bf, and raising the bar past the evidence turns the same tract into `ambiguous`
    rather than making it disappear: a tract that fails a test is reported with the test it
    failed, never hidden.
    """
    positions = [100, 110, 120, 130, 140, 150]
    records = [_read_at(f"r{i}", 100, {p: ("G" if 110 <= p <= 140 else "A") for p in positions})
               for i in range(6)]
    bam = _bam(tmp_path, records)
    sites = _sites_tsv(tmp_path / "sites.tsv", {0: positions})
    out = tmp_path / "tracts.tsv"
    argv = ["--sites", sites, "--bam", bam, "--sample", "S1", "-o", str(out),
            "--samtools", samtools]

    def body():
        return [ln for ln in out.read_text().splitlines() if not ln.startswith("#")]

    assert gc.main(argv) == 0
    assert len(body()) == 2 and "\tgene_conversion\t" in body()[1]

    assert gc.main(argv + ["--min-bf", "99"]) == 0
    assert len(body()) == 2 and "\tambiguous\t" in body()[1]

    # Past the reporting threshold the row goes away entirely.
    assert gc.main(argv + ["--report-bf", "99"]) == 0
    assert len(body()) == 1


def test_the_nextflow_module_writes_the_same_header_the_tool_does(repo_root):
    """A stub run is how the DAG is tested, and an empty cohort is written by the module rather
    than by the tool. A header that has drifted from this file's is worse than none: the run
    still produces a file, the report still parses it, and every new column reads as missing.
    """
    module = (repo_root / "modules" / "gene_conversion.nf").read_text()

    body = module.split("def gconvHeader()", 1)[1].split("}", 1)[0]
    assert re.findall(r"'([a-z0-9_]+)'", body) == gc.COLUMNS
    assert "printf 'sample" not in module, "the header belongs in gconvHeader(), not restated"


def test_the_report_download_writes_the_same_header_the_tool_does(repo_root):
    """The panel's "download tracts (TSV)" button rebuilds the file from the parsed payload.

    Its column list is hand written, so a column added here and not there silently mislabels
    every field after the gap in a file someone will open in a spreadsheet and believe.
    """
    js = (repo_root / "bin" / "report_assets" / "js" / "13_drug_kraken.js").read_text()

    # Anchored on the file the button writes, because the drug-resistance panel in the same
    # module builds its own download the same way and matching the first one found would test
    # that instead and pass for the wrong reason.
    before = js.split("'gene_conversion.tsv'", 1)[0]
    block = before.rsplit("var hdr=[", 1)[1].split("]", 1)[0]
    cohort = load_script("gconv_cohort")
    assert re.findall(r"'([a-z0-9_]+)'", block) == gc.COLUMNS + cohort.COHORT_COLUMNS


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
    body = [ln for ln in out.read_text().splitlines() if not ln.startswith("#")]
    assert body[1].startswith("S1\t0\t")


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
    conversions had each one counting the other's donor alleles as evidence against itself, and
    both came out as mismapping. Two conversions in one locus is an ordinary outcome, not a
    pathological input.
    """
    positions = list(range(1, 13))
    counts = _flat_counts(positions, donor_positions=set(range(1, 5)) | set(range(9, 13)))
    per_read = {"r1": {1: "donor", 2: "donor", 5: "acceptor"}}
    tracts = [[1, 2, 3, 4], [9, 10, 11, 12]]

    in_any = {p for t in tracts for p in t}
    for tract in tracts:
        ev = gc.tract_evidence(tract, positions, counts, per_read, in_any)
        assert ev["donor_af_outside"] == 0.0, "the other tract was counted as outside"


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


# --------------------------------------------------------------------------- #
# The high-identity blind spot
# --------------------------------------------------------------------------- #

def test_a_depleted_run_is_reported_when_the_donor_holds_the_reads():
    """A conversion longer than the library insert loses its reads to the donor.

    Once the tract makes the acceptor identical to the donor, a read pair falling entirely inside
    it has no unique anchor and the aligner assigns it arbitrarily. Measured on a simulated 99.7%
    paralog pair: the acceptor's tract went from 84x to 37x while the donor's matching region went
    from 83x to 127x, and every allele-based signal vanished. Reporting the depletion is the only
    thing that keeps that from being a silent zero.
    """
    positions = [10, 20, 30, 40, 50]
    counts = _counts(dict.fromkeys(positions, None),
                     depth={10: 40, 20: 1, 30: 0, 40: 2, 50: 40})
    donor = {10: 40, 20: 80, 30: 90, 40: 85, 50: 40}

    runs = gc.depleted_runs(positions, counts, donor, min_depth=5, min_sites=3)
    assert runs == [[20, 30, 40]]


def test_a_depletion_running_to_the_last_site_is_still_reported():
    """A run is closed when a healthy site ends it, and the last site of the locus ends nothing.

    A depletion reaching the end of the diagnostic sites has no site after it to trigger the
    close, so it is only reported if the loop flushes what it is still holding when it stops.
    Miss that and the run is dropped in silence: not a wrong verdict, no verdict, and the case is
    not exotic. A conversion running off the end of the aligned stretch is exactly where the
    reads are hardest to place.
    """
    positions = [10, 20, 30, 40, 50]
    counts = _counts(dict.fromkeys(positions, None),
                     depth={10: 40, 20: 40, 30: 1, 40: 0, 50: 2})
    donor = {10: 40, 20: 40, 30: 90, 40: 95, 50: 85}

    assert gc.depleted_runs(positions, counts, donor, min_depth=5, min_sites=3) == [[30, 40, 50]]


def test_a_run_at_the_end_still_has_to_be_long_enough():
    """The flush at the end is not a way around `min_sites`."""
    positions = [10, 20, 30, 40]
    counts = _counts(dict.fromkeys(positions, None), depth={10: 40, 20: 40, 30: 1, 40: 0})
    donor = {10: 40, 20: 40, 30: 90, 40: 95}

    assert gc.depleted_runs(positions, counts, donor, min_depth=5, min_sites=3) == []


# Every constant in `depleted_runs` is a threshold, and each one reads perfectly well when it is
# a step out. The cases below sit exactly ON each of them, in the direction that has to pass.
#
# They all share a shape: healthy flanks either side of a depleted middle. The flanks are not
# decoration. The donor has to have GAINED over its own level ELSEWHERE, so a locus whose every
# site is in the run has no elsewhere to be compared with and is never reported at all, which is
# how the first version of these tests managed to assert nothing four times over.

FLANKS, MIDDLE = [10, 20, 70, 80], [30, 40, 50]


def _shift(acc, don, **kw):
    """A locus described by its two depth profiles: {position: acceptor}, {position: donor}."""
    positions = sorted(acc)
    counts = {p: {"acceptor": 0, "donor": 0, "other": 0, "depth": acc[p], "donor_af": None}
              for p in positions}
    return gc.depleted_runs(positions, counts, don, **kw)


def _locus(acc_mid, don_mid, acc_flank=40, don_flank=40):
    acc = {p: acc_flank for p in FLANKS} | {p: acc_mid for p in MIDDLE}
    don = {p: don_flank for p in FLANKS} | {p: don_mid for p in MIDDLE}
    return acc, don


def test_a_depletion_over_the_minimum_number_of_sites_is_reported():
    """The base case the rest are nudged away from, and it carries two of the checks on its own.

    Both marks of a shift are required, not either: loosen the acceptor test and the healthy
    flanks become part of the run, which leaves no baseline and reports nothing at all.
    """
    assert _shift(*_locus(0, 100), min_depth=5, min_sites=3) == [MIDDLE]


@pytest.mark.parametrize("acc_mid,depleted", [(4, True), (5, False)])
def test_a_site_at_exactly_the_depth_floor_is_not_depleted(acc_mid, depleted):
    """`min_depth` is the depth a site needs to count as covered, so at it, it is covered."""
    runs = _shift(*_locus(acc_mid, 100), min_depth=5, min_sites=3)

    assert (runs == [MIDDLE]) is depleted


@pytest.mark.parametrize("don_mid,enough", [(10, True), (9, False)])
def test_the_two_copies_need_twice_the_floor_between_them_before_anything_is_said(don_mid, enough):
    """Under `2 * min_depth` across both copies there is not enough sequence to tell reads that
    moved from reads that never arrived, and exactly that much is enough.

    The donor's level in the flanks is kept low here so that the enrichment check clears in both
    arms and only the one condition under test decides the answer.
    """
    runs = _shift(*_locus(0, don_mid, don_flank=4), min_depth=5, min_sites=3)

    assert (runs == [MIDDLE]) is enough


@pytest.mark.parametrize("max_ratio,lopsided", [(0.4, True), (0.39, False)])
def test_a_site_holding_exactly_the_permitted_share_still_counts_as_lopsided(max_ratio, lopsided):
    """`max_ratio` is the largest share of the reads the acceptor may still hold: 4 of 10 here."""
    runs = _shift(*_locus(4, 6, don_flank=4), min_depth=5, max_ratio=max_ratio, min_sites=3)

    assert (runs == [MIDDLE]) is lopsided


def test_a_well_covered_acceptor_is_not_depleted_however_deep_the_donor_is():
    """The share of the reads is not on its own a statement that any were lost.

    A tandem repeat can pull an order of magnitude more reads than its paralog while the paralog
    itself is perfectly covered. Nothing moved: there is simply more of the other copy.
    """
    assert _shift(*_locus(40, 400), min_depth=5, min_sites=3) == []


@pytest.mark.parametrize("max_local,inside", [(0.6, True), (0.59, False)])
def test_a_site_at_exactly_the_local_ceiling_is_part_of_the_depletion(max_local, inside):
    """The run is cut where the acceptor is back up to `max_local` of its own level elsewhere,
    and a site sitting exactly at that fraction is still down. 24 is 0.6 of 40."""
    runs = _shift(*_locus(24, 400), min_depth=25, max_local=max_local, min_sites=3)

    assert (runs == [MIDDLE]) is inside


def test_two_healthy_sites_are_enough_to_know_what_the_locus_runs_at():
    """With two the baseline is a measurement and the local ceiling applies; with fewer it is not
    and the ceiling is skipped, which lets through a site that is down but not down far.

    Here the first site of the would-be run sits at 30 against a baseline of 40, well above the
    ceiling, so it breaks the run and what is left is too short to report. Refusing to measure a
    baseline from two sites reports it instead.
    """
    acc = {10: 40, 30: 30, 40: 0, 50: 0, 80: 40}
    don = {10: 40, 30: 400, 40: 400, 50: 400, 80: 40}

    assert _shift(acc, don, min_depth=35, min_sites=3) == []


@pytest.mark.parametrize("don_mid,gained", [(50, True), (49, False)])
def test_a_donor_gaining_exactly_the_required_amount_has_gained(don_mid, gained):
    """`min_enrichment` is the rise the donor must show over its own level elsewhere: 50 over 40
    is exactly 1.25."""
    runs = _shift(*_locus(0, don_mid), min_depth=5, min_enrichment=1.25, min_sites=3)

    assert (runs == [MIDDLE]) is gained


def test_a_donor_with_no_reads_outside_the_run_has_not_gained_anything():
    """There is no level to have risen above, and zero clears any ratio you care to name.

    This is the shape a paralog on a contig edge has: the donor copy is only covered where the
    reads piled up, and nowhere else. Comparing against nothing and calling the result a rise
    turns every one of those into a reported shift.
    """
    acc = {p: 40 for p in FLANKS} | {p: 0 for p in MIDDLE}
    don = {p: 0 for p in FLANKS} | {p: 100 for p in MIDDLE}

    assert _shift(acc, don, min_depth=5, min_sites=3) == []


def test_a_locus_covered_at_neither_copy_says_nothing():
    positions = [10, 20, 30]
    counts = {p: {"acceptor": 0, "donor": 0, "other": 0, "depth": 0, "donor_af": None}
              for p in positions}

    assert gc.depleted_runs(positions, counts, {}, min_depth=5, min_sites=3) == []


def test_a_locus_covered_on_both_sides_reports_no_depletion():
    """Ordinary coverage on both copies is not a shift, however low it is overall."""
    positions = [10, 20, 30, 40]
    counts = _counts(dict.fromkeys(positions, 0.0), depth=40)
    donor = dict.fromkeys(positions, 40)

    assert gc.depleted_runs(positions, counts, donor, min_depth=5, min_sites=3) == []


def _depth_case(acc, don):
    positions = sorted(acc)
    counts = {p: {"acceptor": 0, "donor": 0, "other": 0, "depth": acc[p], "donor_af": None}
              for p in positions}
    return positions, counts, don


def test_a_low_coverage_locus_is_not_a_shift_of_reads():
    """Regression, and the reason the test above is not enough.

    The check used to be an absolute one: the acceptor under `min_depth` while the donor held
    most of the reads. That is true of every low-coverage paralog whether or not anything moved,
    and on a clean negative control over 411 real pairs it reported two shifts on a locus running
    at 4x throughout, where the acceptor was barely dipping (0.75 of its own level) and the donor
    barely gaining (1.12 of its own). Both numbers now have to say the reads MOVED.
    """
    positions, counts, donor = _depth_case(
        {10: 5, 20: 5, 30: 4, 40: 4, 50: 4, 60: 4, 70: 5, 80: 5},
        {10: 14, 20: 14, 30: 16, 40: 16, 50: 16, 60: 16, 70: 14, 80: 14})

    assert gc.depleted_runs(positions, counts, donor, min_depth=5, min_sites=3) == []


def test_an_acceptor_dropout_the_donor_did_not_absorb_is_not_a_shift():
    """A hole in the coverage is a hole. Reads that moved arrive somewhere."""
    positions, counts, donor = _depth_case(
        {10: 80, 20: 80, 30: 1, 40: 0, 50: 1, 60: 2, 70: 80, 80: 80},
        dict.fromkeys([10, 20, 30, 40, 50, 60, 70, 80], 80))

    assert gc.depleted_runs(positions, counts, donor, min_depth=5, min_sites=3) == []


def test_the_measured_migration_is_still_reported():
    """The case the check exists for, at the depths measured on the 99.7% paralog pair: the
    acceptor's tract fell from 84x to 37x while the donor's matching region rose from 83x to
    127x, and the reads were conserved between them."""
    positions, counts, donor = _depth_case(
        {10: 84, 20: 84, 30: 37, 40: 37, 50: 37, 60: 37, 70: 84, 80: 84},
        {10: 83, 20: 83, 30: 127, 40: 127, 50: 127, 60: 127, 70: 83, 80: 83})

    assert gc.depleted_runs(positions, counts, donor, min_depth=50, min_sites=3) == [[30, 40, 50, 60]]


def test_the_baseline_comes_from_the_sites_that_are_not_suspect():
    """A depletion covering half the locus must not drag down the level it is compared against.

    Taken over every site, the median is inside the depletion itself, and the run then breaks up
    on its own deepest site.
    """
    positions, counts, donor = _depth_case(
        {10: 40, 20: 1, 30: 0, 40: 2, 50: 40},
        {10: 40, 20: 80, 30: 90, 40: 85, 50: 40})

    assert gc.depleted_runs(positions, counts, donor, min_depth=5, min_sites=3) == [[20, 30, 40]]


def test_depletion_needs_the_reads_to_be_somewhere():
    """A locus with no coverage on either copy is missing data, not a shift of reads."""
    positions = [10, 20, 30, 40]
    counts = _counts(dict.fromkeys(positions, None), depth=0)
    donor = dict.fromkeys(positions, 0)

    assert gc.depleted_runs(positions, counts, donor, min_depth=5, min_sites=3) == []


def test_the_output_records_the_settings_that_produced_it(tmp_path, samtools):
    """A results table whose verdicts depend on seventeen settings, and which does not say what
    they were, cannot be checked against another run or reproduced a year later.

    Written as a `#` header so the file stays a plain TSV: every reader here skips those, and so
    does the report's parser.
    """
    records = [_read_at(f"r{i}", 100, {100: "G", 110: "G", 120: "G"}, length=30) for i in range(6)]
    bam = _bam(tmp_path, records)
    sites = _sites_tsv(tmp_path / "sites.tsv", {0: [100, 110, 120]})
    out = tmp_path / "tracts.tsv"

    assert gc.main(["--sites", sites, "--bam", bam, "--sample", "S1", "-o", str(out),
                    "--samtools", samtools, "--min-bf", "4.5", "--mut-rate", "0.002"]) == 0

    head = out.read_text().splitlines()[0]
    assert head.startswith("# gene_conversion.py")
    assert "min_bf=4.5" in head and "mut_rate=0.002" in head
    for name in gc.PROVENANCE_ARGS:
        assert f"{name}=" in head, f"{name} moves the numbers but is not recorded"
    # and the table underneath is still exactly a TSV
    assert out.read_text().splitlines()[1] == "\t".join(gc.COLUMNS)


# --------------------------------------------------------------------------- #
# The settings are checked before the work starts                              #
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("flag,value,says", [
    ("--mut-rate", "0", "--mut-rate"),
    ("--mut-rate", "1", "--mut-rate"),
    ("--mut-rate", "-0.1", "--mut-rate"),
    ("--prior", "1.5", "--prior"),
    ("--prior", "-0.1", "--prior"),
    ("--mean-tract-bp", "0", "--mean-tract-bp"),
    ("--mean-tract-bp", "-500", "--mean-tract-bp"),
    ("--max-tracts", "0", "--max-tracts"),
])
def test_a_setting_the_model_cannot_use_is_refused_before_any_reading_is_done(flag, value, says):
    """These guards exist to fail while the message can still name the flag.

    The model raises on all of them, but it raises after the BAM has been read, from inside an
    expression, naming a variable nobody typed. Every one of these was uncovered: the guards were
    written and nothing checked that they fire, so a wrong comparison in any of them would have
    turned a clean refusal back into a crash halfway through a sample.
    """
    with pytest.raises(SystemExit) as e:
        gc.parse_args(["--sites", "s.tsv", "--bam", "b.bam", "--sample", "S1", "-o", "o.tsv",
                       flag, value])
    assert e.value.code != 0


@pytest.mark.parametrize("flag,value", [("--prior", "0"), ("--prior", "1"), ("--max-tracts", "1")])
def test_a_setting_at_the_edge_of_what_is_allowed_is_allowed(flag, value):
    """The other half of the guards above, and the half that stays broken quietly.

    A guard one step too strict refuses a setting the tool documents, and the user reads a
    refusal for a value the help told them to use. A prior of 0 and a prior of 1 are both
    meaningful instructions, and one tract is the smallest number of tracts to look for.
    """
    a = gc.parse_args(["--sites", "s.tsv", "--bam", "b.bam", "--sample", "S1", "-o", "o.tsv",
                       flag, value])

    assert getattr(a, flag.lstrip("-").replace("-", "_")) == float(value)


def test_a_deletion_marker_is_priced_as_an_indel_not_as_a_substitution(tmp_path, samtools):
    """The KIND of each site has to reach the model, not just the position.

    A base one copy has and the other does not is far rarer than a substitution, so a tract
    carrying one has that much more going for it. Reading the column and then ignoring it, or
    using it the wrong way round, changes no column in the output and no verdict here: only the
    evidence moves, and nothing was asking about the evidence.

    The tract has to leave sites outside it, because on a locus a tract covers entirely the
    binding alternative is mismapping and the pricing of substitutions never comes into it. The
    first version of this test made exactly that mistake and reported the same number twice.
    """
    pos = [100, 110, 120, 130, 140, 150, 160, 170]
    # 25M covers 95-119, the 1D is 120, and 55M covers 121-175
    read_at = lambda ref: ref - 95 if ref < 120 else 25 + (ref - 121)
    seq = ["A"] * 80
    for p in (110, 130):
        seq[read_at(p)] = "G"                       # the donor base, either side of the deletion
    records = [_sam_line(f"r{i}", 95, "25M1D55M", "".join(seq), qual="I" * 80) for i in range(14)]
    bam = _bam(tmp_path, records)

    def evidence(kinds, name):
        out = tmp_path / f"{name}.tsv"
        sites = _sites_tsv(tmp_path / f"{name}.sites", {0: pos}, kinds=kinds)
        assert gc.main(["--sites", sites, "--bam", bam, "--sample", "S1", "-o", str(out),
                        "--samtools", samtools, "--min-sites", "2", "--report-bf", "-99"]) == 0
        body = [ln for ln in out.read_text().splitlines() if not ln.startswith("#")]
        row = dict(zip(body[0].split("\t"), body[1].split("\t")))
        assert (row["start"], row["end"]) == ("110", "130")
        return float(row["log10_bf"])

    as_indel = evidence({120: "del"}, "indel")
    as_snp = evidence({}, "snp")

    # An indel is priced at INDEL_MUT_FACTOR of a substitution, so the tract that carries one is
    # ahead by exactly that ratio in log10. Asserting the amount rather than the direction is
    # what makes this a test of the pricing instead of a test that the two runs differ.
    assert as_indel - as_snp == pytest.approx(-math.log10(gc.gm.INDEL_MUT_FACTOR), abs=0.02)


def test_a_coverage_shift_row_reports_the_stretch_it_actually_covers(tmp_path, samtools):
    """The one verdict built by hand rather than by the model, so its columns are its own.

    Nothing else checks them: the run's own coordinates, the span they imply and the shallowest
    depth in it are all assembled here, and each is a min or a max that reads perfectly well when
    it is the wrong one. A row saying the acceptor is depleted over a single base, or over a
    stretch running backwards, is what a reader would be handed.
    """
    positions = [1000, 1010, 1020, 1030]
    sites = _sites_tsv(tmp_path / "sites.tsv", {0: positions})
    # The acceptor runs deep at 1000 and is empty over the other three. The donor, at 1200 along,
    # has to be ENRICHED over those three against its own baseline, not merely present: an
    # acceptor that is simply shallow everywhere is not a shift of reads, and the check says so.
    records = [_sam_line(f"acc{i}", 995, "10M", "A" * 10, qual="I" * 10) for i in range(30)]
    # Two stragglers at 1010, so the run's depths are not all the same and the shallowest site
    # in it is a different number from the deepest.
    records += [_sam_line(f"thin{i}", 1010, "1M", "A", qual="I") for i in range(2)]
    records += [_sam_line(f"base{i}", 2195, "50M", "A" * 50, qual="I" * 50) for i in range(5)]
    records += [_sam_line(f"don{i}", 2205, "40M", "A" * 40, qual="I" * 40) for i in range(40)]
    bam = _bam(tmp_path, records)
    out = tmp_path / "tracts.tsv"

    assert gc.main(["--sites", sites, "--bam", bam, "--sample", "S1", "-o", str(out),
                    "--samtools", samtools, "--min-sites", "3", "--min-depth", "5"]) == 0

    body = [ln for ln in out.read_text().splitlines() if not ln.startswith("#")]
    header = body[0].split("\t")
    shifts = [dict(zip(header, ln.split("\t"))) for ln in body[1:]]
    shifts = [r for r in shifts if r["verdict"] == "coverage_shift"]

    assert shifts, "the depletion was not reported at all"
    row = shifts[0]
    assert (int(row["start"]), int(row["end"])) == (1010, 1030)
    assert int(row["span_bp"]) == 1030 - 1010 + 1
    assert int(row["n_sites"]) == 3 and int(row["n_sites_outside"]) == 1
    assert int(row["min_depth"]) == 0, "the shallowest site in the run, not the deepest"


def _one_locus(tmp_path, tract, name, positions=(100, 110, 120, 130, 140, 150, 160, 170)):
    """A BAM and a sites file for one locus whose molecules carry the donor base over `tract`."""
    seq = ["A"] * 90
    for p in tract:
        seq[p - 95] = "G"
    records = [_sam_line(f"r{i}", 95, "90M", "".join(seq), qual="I" * 90) for i in range(14)]
    return (_bam(tmp_path, records),
            _sites_tsv(tmp_path / f"{name}.sites", {0: list(positions)}))


def _only_row(path):
    body = [ln for ln in path.read_text().splitlines() if not ln.startswith("#")]
    return dict(zip(body[0].split("\t"), body[1].split("\t")))


def test_a_tract_touching_one_end_of_a_locus_does_not_cover_it(tmp_path, samtools):
    """Covering the locus means covering ALL of it, and it is not a detail of wording.

    A tract over every diagnostic site has exactly the likelihood of every read having come from
    the donor, so the model cannot separate the two and says so. A tract that merely STARTS at
    the first site is nothing of the kind, and treating it as such takes a clean call and returns
    ambiguous, which is the answer that looks like caution rather than like a defect.
    """
    bam, sites = _one_locus(tmp_path, {100, 110, 120}, "edge")
    out = tmp_path / "edge.tsv"

    assert gc.main(["--sites", sites, "--bam", bam, "--sample", "S1", "-o", str(out),
                    "--samtools", samtools, "--min-sites", "2", "--report-bf", "-99"]) == 0

    row = _only_row(out)
    assert (row["start"], row["end"]) == ("100", "120")
    assert row["verdict"] == "gene_conversion"

    # Covering the locus only decides anything once the evidence is short of a call, so the
    # threshold is put out of reach and the reason is read instead of the verdict.
    weak = tmp_path / "edge_weak.tsv"
    assert gc.main(["--sites", sites, "--bam", bam, "--sample", "S1", "-o", str(weak),
                    "--samtools", samtools, "--min-sites", "2", "--report-bf", "-99",
                    "--min-bf", "99"]) == 0

    assert "covers every diagnostic site" not in _only_row(weak)["reason"]


def test_a_locus_searched_exhaustively_is_not_told_it_was_not(tmp_path, samtools):
    """The note about coarse breakpoints belongs on the loci that got them.

    At stride 1 every interval was evaluated, and appending the note anyway tells a reader their
    breakpoints were placed every 1 diagnostic sites, which is both untrue and unreadable.
    """
    bam, sites = _one_locus(tmp_path, {110, 120, 130}, "fine")
    out = tmp_path / "fine.tsv"

    assert gc.main(["--sites", sites, "--bam", bam, "--sample", "S1", "-o", str(out),
                    "--samtools", samtools, "--min-sites", "2", "--report-bf", "-99"]) == 0

    assert "resolved to every" not in _only_row(out)["reason"]


@pytest.mark.parametrize("bf,mismap,kept", [
    (3.0, 0.0, True),        # exactly the reporting threshold
    (2.9999, 0.0, False),
    (0.0, 0.2, True),        # exactly the mismapping rate
    (0.0, 0.1999, False),
    (2.9999, 0.1999, False), # neither, which is the only way to be dropped
])
def test_a_fit_is_written_out_at_the_thresholds_not_past_them(bf, mismap, kept):
    """`--report-bf` and `--min-mismap` decide whether a row exists at all.

    Both were unreachable from outside: the only view of either number is a column rounded to two
    decimals, where a fit sitting on the threshold and a fit just past it read the same, so the
    first attempt at this test asserted nothing twice over.
    """
    fit = {"log10_bf": bf, "mismap_frac": mismap}

    assert gc.worth_reporting(fit, report_bf=3.0, min_mismap=0.2) is kept
