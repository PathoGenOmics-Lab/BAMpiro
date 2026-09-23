"""Unit tests for bin/gconv_cohort.py (reading a cohort's conversion calls together).

Two revisions are made here and nowhere else, in opposite directions, and each is only possible
because the other samples are present:

* an event in nearly every sample is demoted to `reference_artifact`. The same conversion did not
  happen independently in fifty isolates; the reference is wrong there, or the aligner does the
  same thing to everybody. A single sample cannot tell those apart from a real tract, because
  mismapped reads look identical either way.
* a tract too weak for one sample to commit to is promoted when the SAME tract, at the SAME
  coordinates, is called outright in another. Contamination and index hopping do not reproduce a
  specific tract across independent libraries.

The second is the dangerous one, so most of what follows is about the ways it must NOT fire.
"""

from __future__ import annotations

import pytest

from conftest import load_script

gcc = load_script("gconv_cohort")


def tract(sample, pair="1", contig="chr", start=1000, end=1200, verdict="gene_conversion",
          bf="12.0", don_start=5000, don_end=5200, bf_null="100.0", n_sites="5", **extra):
    row = {"sample": sample, "pair_id": pair, "contig": contig, "donor": "chr",
           "verdict": verdict, "reason": "because", "start": str(start), "end": str(end),
           "don_start": str(don_start), "don_end": str(don_end),
           "log10_bf": bf, "log10_bf_vs_null": bf_null, "n_sites": n_sites,
           "tract_af": "0.9", "mismap_frac": "0.01"}
    row.update(extra)
    return row


def locus(sample, pair="1", mismap="0.01", bf="1.0", contig="chr"):
    return {"sample": sample, "pair_id": pair, "contig": contig, "donor": contig,
            "log10_bf": bf, "mismap_frac": mismap, "mut_rate": "0.0003"}


def cohort(n, **kw):
    """Per-locus rows for `n` samples, which is how the cohort size is known."""
    return [locus(f"S{i}", **kw) for i in range(n)]


# ------------------------------------------------------------------ grouping


def test_the_same_stretch_in_different_samples_is_one_event():
    rows = [tract("S1", start=1000, end=1200), tract("S2", start=1000, end=1200)]

    out, _ = gcc.annotate(rows, cohort(10))

    assert out[0]["event_id"] == out[1]["event_id"]
    assert out[0]["event_samples"] == 2


def test_the_same_stretch_reported_through_two_paralog_pairs_is_still_one_event():
    """A gene family reports the same converted stretch once per relationship, measured at 1.2
    rows per event on the real genome and up to 3. Grouping on the acceptor's coordinates rather
    than on the pair collapses them, so counting samples counts samples and not relationships."""
    rows = [tract("S1", pair="1", start=1000, end=1200),
            tract("S1", pair="7", start=1000, end=1200)]

    out, _ = gcc.annotate(rows, cohort(10))

    assert out[0]["event_id"] == out[1]["event_id"]
    assert out[0]["event_samples"] == 1, "one sample, two relationships, one sample"


def test_tracts_that_do_not_overlap_are_different_events():
    rows = [tract("S1", start=1000, end=1200), tract("S1", start=5000, end=5200)]

    out, _ = gcc.annotate(rows, cohort(10))

    assert out[0]["event_id"] != out[1]["event_id"]


def test_the_same_coordinates_on_another_contig_are_a_different_event():
    rows = [tract("S1", contig="ctg1"), tract("S2", contig="ctg2")]

    out, _ = gcc.annotate(rows, cohort(10))

    assert out[0]["event_id"] != out[1]["event_id"]


# ------------------------------------------------- demoting what everyone has


def test_an_event_in_the_whole_cohort_is_the_reference_not_the_isolates():
    rows = [tract(f"S{i}") for i in range(10)]

    out, _ = gcc.annotate(rows, cohort(10))

    assert all(r["cohort_verdict"] == "reference_artifact" for r in out)
    assert "10 of 10 samples" in out[0]["reason"]


def test_an_event_in_a_few_samples_is_left_alone():
    rows = [tract("S0"), tract("S1")]

    out, _ = gcc.annotate(rows, cohort(10))

    assert all(r["cohort_verdict"] == "gene_conversion" for r in out)


def test_a_cohort_too_small_to_argue_from_recurrence_is_not_argued_from():
    """With three samples, "in all of them" is three. That is not evidence about the reference,
    and demoting a real finding because the cohort is small would be worse than saying nothing."""
    rows = [tract(f"S{i}") for i in range(3)]

    out, _ = gcc.annotate(rows, cohort(3), min_samples=5)

    assert all(r["cohort_verdict"] == "gene_conversion" for r in out)


def test_a_tract_inside_another_does_not_shorten_the_event():
    """The running end of an event is the furthest end seen, not the latest one.

    Tracts are walked in order of where they start, so a short one nested inside a long one comes
    second and would pull the event's end back to its own. The next tract then falls outside an
    event it belongs to and starts a new one, which splits a single recurrent event into two
    that each look rare.
    """
    rows = [tract("S0", start=1000, end=2000),
            tract("S1", start=1100, end=1200),        # nested inside the first
            tract("S2", start=1900, end=2100)]        # still inside the first, not the second

    out, _ = gcc.annotate(rows, cohort(10))

    assert len({r["event_id"] for r in out}) == 1, "one stretch, one event"


def test_merging_two_views_of_a_donor_widens_the_candidate_to_cover_both():
    """The merged span is the union of the two, because it is what later rows are tested against.

    Anything narrower and a third relationship pointing at the same stretch falls outside it and
    is counted as a separate source, which turns one donor into two and makes a resolved call
    ambiguous.
    """
    rows = [tract("S0", pair="1", don_start=5000, don_end=5100),
            tract("S0", pair="2", don_start=5090, don_end=5300),
            tract("S0", pair="3", don_start=5200, don_end=5250)]

    out, _ = gcc.annotate(rows, cohort(10))

    assert all(r["n_donors"] == 1 for r in out), "one stretch of donor, one candidate"


def test_a_row_with_no_markers_does_not_erase_the_evidence_of_the_one_it_joins():
    """A tract row can carry no usable count, and then it has no evidence per marker to offer.

    Joining a candidate has to leave that candidate's own evidence alone. Letting the empty row
    through instead replaces a measured number with nothing, and the source that was resolved
    stops being resolvable.
    """
    rows = [tract("S0", pair="1", don_start=5000, don_end=5200, bf_null="200.0", n_sites="20"),
            tract("S0", pair="2", don_start=5100, don_end=5300, bf_null="", n_sites="0"),
            tract("S0", pair="3", don_start=9000, don_end=9200, bf_null="12.0", n_sites="9")]

    out, _ = gcc.annotate(rows, cohort(10))

    assert out[0]["donor_rank"] == 1, "the well-supported stretch is still the best candidate"
    assert out[0]["donor_call"] == "resolved"


@pytest.mark.parametrize("gap,together", [(5, True), (6, False), (0, True)])
def test_slack_is_the_tolerance_it_says_it_is(gap, together):
    """A gap of exactly `--slack` bases still makes two tracts the same event.

    Both directions matter. Too tight and one event in two samples becomes two events, so it can
    never look recurrent and never be demoted; too loose and neighbouring events merge and a real
    tract inherits another one's verdict.
    """
    rows = [tract("S0", start=1000, end=1200),
            tract("S1", start=1200 + gap, end=1400 + gap)]

    out, _ = gcc.annotate(rows, cohort(10), slack=5)

    assert (out[0]["event_id"] == out[1]["event_id"]) is together


def test_the_samples_with_nothing_to_report_still_count_towards_the_cohort():
    """The denominator of the ubiquity rule is how many samples were RUN, not how many had
    something to say. A clean sample writes no tract row, so counting the samples named in the
    tract files alone turns an event in 5 of 8 into an event in 5 of 5 and demotes a real
    conversion to a reference artifact, which is this pass doing the opposite of its job."""
    rows = [tract(f"S{i}") for i in range(5)]

    out, n = gcc.annotate(rows, cohort(8))

    assert n == 8
    assert all(r["cohort_verdict"] == "gene_conversion" for r in out)
    assert out[0]["event_frac"] == 0.625


def test_a_cohort_of_unknown_size_is_not_a_cohort_of_the_size_that_spoke():
    """Without the per-locus rows there is no census, and the quiet samples cannot be counted at
    all. An unknown denominator is not a small one, so recurrence is left unapplied rather than
    computed against the only number to hand."""
    rows = [tract(f"S{i}") for i in range(5)]

    out, n = gcc.annotate(rows, [])

    assert n is None
    assert all(r["cohort_verdict"] == "gene_conversion" for r in out)
    assert all(r["event_frac"] == "" for r in out)

    stated, n_stated = gcc.annotate(rows, [], cohort_size=5)
    assert n_stated == 5
    assert all(r["cohort_verdict"] == "reference_artifact" for r in stated)


def test_a_stated_cohort_size_cannot_be_smaller_than_the_samples_seen():
    rows = [tract(f"S{i}") for i in range(6)]

    out, n = gcc.annotate(rows, [], cohort_size=2)

    assert n == 6, "a size that contradicts the rows read would make the fraction exceed 1"
    assert out[0]["event_frac"] == 1.0


@pytest.mark.parametrize("n,demoted", [(5, True), (4, False), (6, True)])
def test_a_cohort_of_exactly_the_minimum_is_big_enough_to_argue_from(n, demoted):
    """`--gconv_cohort_min_samples` is the size at which recurrence starts being evidence."""
    rows = [tract(f"S{i}") for i in range(n)]

    out, _ = gcc.annotate(rows, cohort(n), min_samples=5)

    assert all(r["cohort_verdict"] == "reference_artifact" for r in out) is demoted


@pytest.mark.parametrize("carrying,demoted", [(9, True), (8, False), (10, True)])
def test_an_event_at_exactly_the_ubiquity_fraction_is_demoted(carrying, demoted):
    """`--gconv_ubiquitous` is the fraction AT WHICH an event becomes a reference artifact."""
    rows = [tract(f"S{i}") for i in range(carrying)]

    out, _ = gcc.annotate(rows, cohort(10), ubiquitous=0.9)

    assert all(r["cohort_verdict"] == "reference_artifact" for r in out) is demoted


@pytest.mark.parametrize("bf,rescued", [("2.0", True), ("1.99", False), ("2.5", True)])
def test_a_tract_at_exactly_the_corroboration_threshold_is_corroborated(bf, rescued):
    """`--corroborated-bf` is the evidence a sub-threshold tract needs before another sample's
    outright call is allowed to speak for it."""
    rows = [tract("S0", verdict="gene_conversion"),
            tract("S1", verdict="ambiguous", bf=bf)]

    out, _ = gcc.annotate(rows, cohort(10), corroborated_bf=2.0)

    assert (out[1]["cohort_verdict"] == "gene_conversion") is rescued


@pytest.mark.parametrize("bf_null,resolved", [("30.0", True), ("29.0", False), ("40.0", True)])
def test_a_donor_ahead_by_exactly_the_margin_is_resolved(bf_null, resolved):
    """`--donor-margin` is the lead a source needs, in evidence per marker, to be named.

    The runner-up sits at 2.0 per marker over 10 markers, so the leader clears the margin of 1.0
    exactly at 30.0 over 10. Below it the two are compatible with the reads and the call is
    ambiguous, which is the honest answer short reads allow.
    """
    rows = [tract("S0", pair="1", don_start=5000, don_end=5200, bf_null=bf_null, n_sites="10"),
            tract("S0", pair="2", don_start=9000, don_end=9200, bf_null="20.0", n_sites="10")]

    out, _ = gcc.annotate(rows, cohort(10), donor_margin=1.0)

    assert (out[0]["donor_call"] == "resolved") is resolved


@pytest.mark.parametrize("first,second", [((5000, 5200), (5200, 5400)),
                                          ((5200, 5400), (5000, 5200))])
def test_donor_spans_that_touch_are_the_same_stretch(first, second):
    """Adjacent self-alignments meeting end to start describe one stretch of donor, not two.

    Treating them as two makes the source ambiguous where it is not, and the tract then reports
    a donor it could have named. Both orders are asked because the comparison is not symmetric in
    the code: rows arrive in whatever order the pair ids fall, and only the second order here
    reaches the test against the candidate's start.
    """
    rows = [tract("S0", pair="1", don_start=first[0], don_end=first[1]),
            tract("S0", pair="2", don_start=second[0], don_end=second[1])]

    out, _ = gcc.annotate(rows, cohort(10))

    assert all(r["n_donors"] == 1 for r in out)


def test_a_merged_candidate_reaches_back_to_the_earlier_of_the_two():
    """The union runs from the earliest start, so a third row to the left of it still belongs.

    Keeping the later start instead leaves the candidate covering only part of what it merged,
    and the next relationship pointing into the part left out is counted as another source.
    """
    rows = [tract("S0", pair="1", don_start=5200, don_end=5400),
            tract("S0", pair="2", don_start=5000, don_end=5250),
            tract("S0", pair="3", don_start=5050, don_end=5100)]

    out, _ = gcc.annotate(rows, cohort(10))

    assert all(r["n_donors"] == 1 for r in out), "one stretch, however it was assembled"


def test_a_row_without_a_verdict_column_is_read_as_having_none():
    """The same case as the missing reason, on the column every later decision is compared with."""
    row = tract("S0")
    del row["verdict"]

    out, _ = gcc.annotate([row], cohort(10))

    assert out[0]["cohort_verdict"] == ""


def test_a_row_without_a_reason_column_still_gets_a_verdict():
    """Files written before a column existed are the normal case for a tool used across runs, and
    a missing reason is a blank one, not something to fall over on the way to the verdict.

    The row's own verdict stands here, which is the path that reads the reason back, so it is
    also the one that meets the absent column.
    """
    row = tract("S0")
    del row["reason"]

    out, _ = gcc.annotate([row], cohort(10))

    assert out[0]["cohort_verdict"] == "gene_conversion"
    assert out[0].get("reason", "") == "", "a reason was invented for a row that had none"


def test_the_recurrence_at_which_an_event_is_demoted_is_a_parameter():
    rows = [tract(f"S{i}") for i in range(7)]

    lenient = gcc.annotate(rows, cohort(10), ubiquitous=0.9)[0]
    strict = gcc.annotate(rows, cohort(10), ubiquitous=0.6)[0]

    assert all(r["cohort_verdict"] == "gene_conversion" for r in lenient)
    assert all(r["cohort_verdict"] == "reference_artifact" for r in strict)


# ---------------------------------------------- promoting what one sample saw


def test_a_weak_tract_is_corroborated_by_a_sample_that_called_it_outright():
    """The floor under `tract_af` cost two of twenty-three benchmark tracts, at 25% and 40%
    frequency, because a minority signal in ONE sample cannot be told from contamination. The
    same tract at the same coordinates in another sample is a different proposition."""
    rows = [tract("S0", verdict="gene_conversion", bf="12.0"),
            tract("S1", verdict="ambiguous", bf="2.4")]

    out, _ = gcc.annotate(rows, cohort(10), min_bf=3.0, corroborated_bf=2.0)

    assert out[1]["cohort_verdict"] == "gene_conversion"
    assert "called outright in 1 other sample" in out[1]["reason"]


def test_a_weak_tract_nobody_else_called_stays_weak():
    rows = [tract("S0", verdict="ambiguous", bf="2.4"),
            tract("S1", verdict="ambiguous", bf="2.4")]

    out, _ = gcc.annotate(rows, cohort(10))

    assert all(r["cohort_verdict"] == "ambiguous" for r in out)


def test_a_sample_does_not_corroborate_itself():
    """Two relationships in one sample report the same stretch twice. Letting one vouch for the
    other would turn a single sample's family echo into cohort evidence, which is circular."""
    rows = [tract("S0", pair="1", verdict="gene_conversion", bf="12.0"),
            tract("S0", pair="7", verdict="ambiguous", bf="2.4")]

    out, _ = gcc.annotate(rows, cohort(10))

    assert out[1]["cohort_verdict"] == "ambiguous"


def test_evidence_far_too_thin_is_not_rescued_by_anyone():
    rows = [tract("S0", verdict="gene_conversion", bf="12.0"),
            tract("S1", verdict="ambiguous", bf="0.4")]

    out, _ = gcc.annotate(rows, cohort(10), corroborated_bf=2.0)

    assert out[1]["cohort_verdict"] == "ambiguous"


def test_corroboration_cannot_resurrect_a_reference_artifact():
    """The two revisions pull in opposite directions and the demotion has to win: an event
    everybody has would otherwise corroborate itself into a cohort-wide finding."""
    rows = [tract(f"S{i}", verdict="gene_conversion") for i in range(9)]
    rows.append(tract("S9", verdict="ambiguous", bf="2.5"))

    out, _ = gcc.annotate(rows, cohort(10))

    assert all(r["cohort_verdict"] == "reference_artifact" for r in out)


def test_mismapping_is_not_promoted_by_corroboration():
    """Only `ambiguous` is eligible. A locus the model explained as donor reads was not a weak
    conversion call, and another sample's tract says nothing about it."""
    rows = [tract("S0", verdict="gene_conversion", bf="12.0"),
            tract("S1", verdict="mismapping", bf="2.5")]

    out, _ = gcc.annotate(rows, cohort(10))

    assert out[1]["cohort_verdict"] == "mismapping"


# ------------------------------------------- what a sample's own reads can carry
#
# On a 185-sample cohort 1,520 of the 1,579 calls were promotions. Each rested on 2 to 5% of the
# reads carrying the donor's bases, fitted at the model's 20% floor, and each was vouched for by a
# single "outright" call from a sample at 0.2x to 2x depth whose tract was read by one to four
# molecules. Reads from a third copy of the family reproduce the same stretch in every library
# mapped to the same reference, so recurrence at that level is the artefact, not the confirmation.


def test_a_call_resting_on_a_handful_of_reads_corroborates_nobody():
    rows = [tract("S0", verdict="gene_conversion", n_sites="5", n_undetermined="4"),
            tract("S1", verdict="ambiguous", bf="2.4")]

    out, _ = gcc.annotate(rows, cohort(10))

    assert out[1]["cohort_verdict"] == "ambiguous"


def test_a_call_resting_on_a_handful_of_reads_is_not_a_call_either():
    """Per-sample files written before the check still hold these as calls; the cohort pass
    applies it too, so re-running only this step over them is enough to correct a run."""
    rows = [tract("S0", verdict="gene_conversion", n_sites="5", n_undetermined="4")]

    out, _ = gcc.annotate(rows, cohort(10))

    assert out[0]["cohort_verdict"] == "ambiguous"
    assert "only 1 of the tract's 5 sites are read by enough molecules" in out[0]["reason"]


def test_a_trickle_is_not_promoted_however_many_samples_carry_it():
    rows = [tract("S0", verdict="gene_conversion", donor_af_in="0.9")]
    rows += [tract(f"S{i}", verdict="ambiguous", bf="6.0", tract_af="0.2", donor_af_in="0.03",
                   reason="only 20% of the reads here carry the tract") for i in range(1, 8)]

    out, _ = gcc.annotate(rows, cohort(40))

    assert all(r["cohort_verdict"] == "ambiguous" for r in out[1:])
    assert "only 3% of the reads over the tract" in out[1]["reason"], \
        "the reason quotes the reads, not the floor the fit was parked on"


def test_corroboration_says_what_the_sample_fell_short_of():
    """The old reason quoted the Bayes factor as short of the threshold on rows at log10 BF 5.6
    against a threshold of 3: what they lacked was the read fraction."""
    rows = [tract("S0", verdict="gene_conversion", bf="12.0"),
            tract("S1", verdict="ambiguous", bf="6.0", tract_af="0.22", donor_af_in="0.25")]

    out, _ = gcc.annotate(rows, cohort(10), min_tract_af=0.25)

    assert out[1]["cohort_verdict"] == "gene_conversion"
    assert out[1]["reason"].startswith("carried by 22% of the reads, under the 25% one sample "
                                       "needs on its own, but the same tract is called outright")
    assert "Bayes factor" not in out[1]["reason"]


# ------------------------------------------------ the samples of one reference


def test_an_event_in_every_sample_of_its_reference_is_the_reference():
    """Measured against the whole cohort, an artefact of a reference carrying 113 of 185 samples
    could never pass 61%, and the ubiquity rule never fired on a two-reference cohort."""
    rows = [tract(f"A{i}", contig="refA") for i in range(10)]
    loci = ([locus(f"A{i}", contig="refA") for i in range(10)]
            + [locus(f"B{i}", contig="refB") for i in range(10)])

    out, _ = gcc.annotate(rows, loci)

    assert all(r["cohort_verdict"] == "reference_artifact" for r in out)
    assert "10 of 10 samples mapped to this reference" in out[0]["reason"]
    assert out[0]["event_frac"] == 1.0


def test_the_same_pair_number_on_two_references_is_two_loci():
    """Pair ids number one reference's paralog map, so pair 7 of one genome is not pair 7 of the
    other, and pooling their mismapping rates described neither."""
    loci = ([locus(f"A{i}", pair="7", mismap="0.30", contig="refA") for i in range(5)]
            + [locus(f"B{i}", pair="7", mismap="0.01", contig="refB") for i in range(5)])
    rows = [tract("A0", pair="7", contig="refA"), tract("B0", pair="7", contig="refB")]

    out, _ = gcc.annotate(rows, loci)

    assert out[0]["cohort_mismap"] == pytest.approx(0.30)
    assert out[1]["cohort_mismap"] == pytest.approx(0.01)


def _diverged_cohort():
    """30 stretches with something reported, a sample calling tracts in ten of them (diverged
    from its reference), and a real event at 500-600 in two clean samples and the diverged one."""
    rows = [tract("F", pair=str(p), start=10_000 + 1_000 * p, end=10_000 + 1_000 * p + 50,
                  verdict="ambiguous") for p in range(30)]
    rows += [tract("D", pair=str(p), start=10_000 + 1_000 * p, end=10_000 + 1_000 * p + 50)
             for p in range(10)]
    rows += [tract(s, pair="99", start=500, end=600) for s in ("S1", "S2", "D")]
    return rows


def test_a_sample_diverged_from_its_reference_does_not_count_towards_an_event():
    """It carries the donor's base at every paralogous locus by inheritance, so it would join
    every event on its reference and say nothing about how often the reference misleads."""
    out, _ = gcc.annotate(_diverged_cohort(), cohort(20))

    event = [r for r in out if r["pair_id"] == "99"]
    assert {r["cohort_verdict"] for r in event if r["sample"] == "D"} == {"divergent_sample"}
    assert all(r["event_samples"] == 2 for r in event)


def test_one_stretch_seen_through_many_relatives_is_one_locus():
    """A gene family reports one converted stretch through a pair per relative. Counted in pairs,
    one conversion seen through eleven relatives put a clean sample on the divergence threshold."""
    rows = [tract("F", pair=str(p), start=10_000 + 1_000 * p, end=10_000 + 1_000 * p + 50,
                  verdict="ambiguous") for p in range(100)]
    rows += [tract("S", pair=str(p), start=500, end=600) for p in range(100, 112)]

    assert gcc.divergent_samples(rows, 0.1, event_of=gcc.group_events(rows)) == {}


# ------------------------------------------------------------- the background


def test_the_locus_statistics_are_pooled_across_the_cohort():
    """The mismapping rate is a property of the reference and the aligner, so the same locus
    troubles everybody the same way. One sample measures it noisily; the cohort median does not."""
    rows = [tract("S0")]
    loci = [locus(f"S{i}", mismap=m) for i, m in enumerate(["0.30", "0.32", "0.31", "0.05",
                                                            "0.33", "0.31", "0.30", "0.32"])]

    out, _ = gcc.annotate(rows, loci)

    assert out[0]["cohort_mismap"] == pytest.approx(0.31, abs=0.01)


def test_samples_that_reported_nothing_still_count_towards_the_cohort():
    """A sample with no tract row is a sample that looked and found nothing, and it has to be in
    the denominator. Counting only the samples that reported would make every event look
    universal, and demote every real finding to a reference artifact."""
    rows = [tract("S0"), tract("S1")]

    out, _ = gcc.annotate(rows, cohort(20))

    assert out[0]["event_frac"] == pytest.approx(0.1)
    assert all(r["cohort_verdict"] == "gene_conversion" for r in out)


def test_a_cohort_with_no_tract_at_all_is_not_an_error():
    out, n = gcc.annotate([], cohort(10))

    assert out == []
    assert n == 10


# ----------------------------------------------------------- which relative it came from


def test_relationships_naming_the_same_stretch_of_donor_are_one_candidate():
    """Overlapping self-alignments of a tandem repeat report one event several times and all
    point at the same donor sequence. There is no source to choose between; the extra rows are
    the family saying the same thing, and on the real genome their evidence ties to two decimals.
    """
    rows = [tract("S0", pair="5", don_start=5000, don_end=5200),
            tract("S0", pair="6", don_start=5010, don_end=5210)]

    out, _ = gcc.annotate(rows, cohort(10))

    assert all(r["n_donors"] == 1 for r in out)
    assert all(r["donor_call"] == "only candidate" for r in out)
    assert sum(r["is_representative"] for r in out) == 1, "one row stands for the event"


def test_a_merged_candidate_is_judged_on_its_best_view_not_its_first():
    """Merging two relationships that name the same donor stretch has to keep the better of them.

    The rows arrive in whatever order the pair ids happen to fall, so a candidate whose first row
    is a poor view of the source would carry that weak number into the ranking against a genuinely
    different candidate and lose to it. Here the same stretch is seen twice, badly then well, and
    a third place in the genome sits between the two readings.
    """
    rows = [tract("S0", pair="1", don_start=5000, don_end=5200, bf_null="10.0", n_sites="10"),
            tract("S0", pair="2", don_start=5100, don_end=5300, bf_null="100.0", n_sites="10"),
            tract("S0", pair="3", don_start=9000, don_end=9200, bf_null="50.0", n_sites="10")]

    out, _ = gcc.annotate(rows, cohort(10))

    assert out[0]["n_donors"] == 2, "the two overlapping views are one candidate"
    winner = next(r for r in out if r["is_representative"])
    assert winner["pair_id"] == "2", "the better view of the merged pair has to win the ranking"
    assert out[2]["donor_rank"] == 2, "the unrelated place comes second on its own evidence"


def test_the_relative_whose_sequence_fits_best_is_the_source():
    """Where the candidates are genuinely different places in the genome, the evidence per marker
    picks between them: how well that donor's sequence accounts for the reads, divided by how
    many markers it had to work with. On the real genome the true source beat a bystander by
    19.1 against 1.4 per site."""
    rows = [tract("S0", pair="333", don_start=9000, don_end=9200, bf_null="12.0", n_sites="9"),
            tract("S0", pair="334", don_start=5000, don_end=5200, bf_null="382.0", n_sites="20")]

    out, _ = gcc.annotate(rows, cohort(10))

    winner = next(r for r in out if r["is_representative"])
    assert winner["pair_id"] == "334"
    assert winner["n_donors"] == 2
    assert winner["donor_call"] == "resolved"
    assert winner["donor_margin"] > 1.0


def test_two_relatives_that_fit_equally_well_are_not_chosen_between():
    """The honest failure. In a family, several relatives can explain a tract equally well and
    short reads do not carry what would settle it. Measured on the real genome: four candidates
    separated by 0.17 per marker, where guessing would be wrong half the time."""
    rows = [tract("S0", pair="124", don_start=1160539, don_end=1161165, bf_null="695.8", n_sites="20"),
            tract("S0", pair="127", don_start=4059984, don_end=4060591, bf_null="692.3", n_sites="20")]

    out, _ = gcc.annotate(rows, cohort(10))

    assert all(r["donor_call"] == "ambiguous" for r in out)
    assert all(r["n_donors"] == 2 for r in out)
    assert sum(r["is_representative"] for r in out) == 1, "still one row per event to read"


def test_the_margin_the_source_has_to_win_by_is_a_parameter():
    rows = [tract("S0", pair="1", don_start=9000, don_end=9200, bf_null="90.0", n_sites="10"),
            tract("S0", pair="2", don_start=5000, don_end=5200, bf_null="100.0", n_sites="10")]

    lenient = gcc.annotate(rows, cohort(10), donor_margin=0.5)[0]
    strict = gcc.annotate(rows, cohort(10), donor_margin=5.0)[0]

    assert next(r for r in lenient if r["is_representative"])["donor_call"] == "resolved"
    assert all(r["donor_call"] == "ambiguous" for r in strict)


def test_each_sample_reads_the_donor_for_itself():
    """Two samples carrying one event each get their own reading of where it came from, because
    the evidence is theirs. Pooling them would let a well-covered sample speak for a thin one."""
    rows = [tract("S0", pair="1", don_start=5000, don_end=5200),
            tract("S1", pair="1", don_start=5000, don_end=5200)]

    out, _ = gcc.annotate(rows, cohort(10))

    assert sum(r["is_representative"] for r in out) == 2
    assert {r["sample"] for r in out if r["is_representative"]} == {"S0", "S1"}


def test_a_row_with_no_donor_coordinates_still_gets_a_verdict():
    """A file written before the donor span was recorded must not lose its rows."""
    rows = [tract("S0"), tract("S1")]
    for r in rows:
        del r["don_start"], r["don_end"]

    out, _ = gcc.annotate(rows, cohort(10))

    assert len(out) == 2
    assert all(r["cohort_verdict"] == "gene_conversion" for r in out)


# ------------------------------------------------------------------ provenance


def test_the_cohort_file_carries_every_sample_settings_and_its_own(tmp_path):
    """A cohort file combines numbers produced by per-sample runs, so it has to say under which
    settings, and a cohort assembled from samples run DIFFERENTLY is a thing a reader must see."""
    a = tmp_path / "a.tsv"
    b = tmp_path / "b.tsv"
    header = "\t".join(["sample", "pair_id", "contig", "donor", "verdict", "reason",
                        "start", "end", "log10_bf"])
    a.write_text("# gene_conversion.py min_bf=3.0\n" + header +
                 "\nS0\t1\tchr\tchr\tgene_conversion\tx\t100\t200\t9.0\n")
    b.write_text("# gene_conversion.py min_bf=9.9\n" + header +
                 "\nS1\t1\tchr\tchr\tgene_conversion\tx\t100\t200\t9.0\n")
    out = tmp_path / "cohort.tsv"

    assert gcc.main(["--tracts", str(a), str(b), "-o", str(out)]) == 0

    notes = [ln for ln in out.read_text().splitlines() if ln.startswith("#")]
    assert "# gene_conversion.py min_bf=3.0" in notes
    assert "# gene_conversion.py min_bf=9.9" in notes, "a mixed cohort must show it was mixed"
    assert any(n.startswith("# gconv_cohort.py") for n in notes)
    assert any("ubiquitous=" in n for n in notes)


def test_a_repeated_settings_line_is_recorded_once(tmp_path):
    """Fifty samples run the same way should not put fifty identical lines in the header."""
    files = []
    header = "\t".join(["sample", "pair_id", "contig", "donor", "verdict", "reason",
                        "start", "end", "log10_bf"])
    for i in range(5):
        f = tmp_path / f"s{i}.tsv"
        f.write_text("# gene_conversion.py min_bf=3.0\n" + header +
                     f"\nS{i}\t1\tchr\tchr\tgene_conversion\tx\t100\t200\t9.0\n")
        files.append(str(f))
    out = tmp_path / "cohort.tsv"

    assert gcc.main(["--tracts", *files, "-o", str(out)]) == 0

    notes = [ln for ln in out.read_text().splitlines() if ln.startswith("#")]
    assert sum(1 for n in notes if n.startswith("# gene_conversion.py")) == 1


def _cohort_files(tmp_path, rows, loci_samples):
    """One tracts file per sample plus the per-locus census, as the pipeline writes them."""
    header = "\t".join(["sample", "pair_id", "contig", "donor", "verdict", "reason",
                        "start", "end", "don_start", "don_end", "log10_bf",
                        "log10_bf_vs_null", "n_sites", "tract_af", "mismap_frac"])
    tracts = []
    by_sample = {}
    for r in rows:
        by_sample.setdefault(r["sample"], []).append(r)
    for sample, rs in by_sample.items():
        f = tmp_path / f"{sample}.tracts.tsv"
        f.write_text("# gene_conversion.py min_bf=3.0\n" + header + "\n"
                     + "".join("\t".join(r[c] for c in header.split("\t")) + "\n" for r in rs))
        tracts.append(str(f))
    lo = tmp_path / "loci.tsv"
    lo.write_text("sample\tpair_id\tlog10_bf\tmismap_frac\tmut_rate\n"
                  + "".join(f"S{i}\t1\t0.1\t0.02\t0.0003\n" for i in range(loci_samples)))
    return tracts, [str(lo)]


def test_the_summary_line_counts_what_the_run_actually_did(tmp_path, capsys):
    """The line on stderr is what a reader sees without opening the file, and every count in it
    is a filter that reads the same when it selects the opposite set.

    The cohort here is deliberately mixed: one event demoted to a reference artifact, one tract
    the cohort talks up from ambiguous, and one left exactly as its sample called it.
    """
    rows = [tract(f"S{i}", start=1000, end=1200) for i in range(9)]     # 9 of 10: demoted
    # S0 reports the second event through two relationships onto the same stretch of donor, so
    # one of those two rows stands for it and the other does not. Without that every row is its
    # own representative and the count says nothing.
    rows += [tract("S0", pair="2", start=5000, end=5200),
             tract("S0", pair="3", start=5000, end=5200, don_start=5100, don_end=5300),
             tract("S1", pair="2", start=5000, end=5200, verdict="ambiguous", bf="2.5")]
    tracts, loci = _cohort_files(tmp_path, rows, loci_samples=10)
    out = tmp_path / "cohort.tsv"

    assert gcc.main(["--tracts", *tracts, "--loci", *loci, "-o", str(out)]) == 0

    line = capsys.readouterr().err
    assert "12 tract row(s)" in line
    assert "over 10 sample(s)" in line
    assert "2 distinct event(s)" in line
    assert "10 verdict(s) revised" in line, "nine demotions and one corroboration"
    assert "(9 as reference artifacts, 1 corroborated across samples)" in line
    assert "collapsing to 11 event/donor call(s)" in line
    # and the warning about an unknown cohort has no business on a run that had the census
    assert "no --loci" not in line


def test_the_summary_says_when_the_cohort_size_was_never_known(tmp_path, capsys):
    rows = [tract("S0"), tract("S1")]
    tracts, _ = _cohort_files(tmp_path, rows, loci_samples=0)
    out = tmp_path / "cohort.tsv"

    assert gcc.main(["--tracts", *tracts, "-o", str(out)]) == 0

    err = capsys.readouterr().err
    assert "over an unknown number of sample(s)" in err
    assert "no --loci and no --cohort-size" in err
    notes = [ln for ln in out.read_text().splitlines() if ln.startswith("#")]
    assert any("cohort_size=unknown" in n for n in notes)


def test_the_settings_line_records_the_cohort_size_it_used(tmp_path):
    rows = [tract("S0")]
    tracts, loci = _cohort_files(tmp_path, rows, loci_samples=7)
    out = tmp_path / "cohort.tsv"

    assert gcc.main(["--tracts", *tracts, "--loci", *loci, "-o", str(out)]) == 0

    notes = [ln for ln in out.read_text().splitlines() if ln.startswith("#")]
    assert any("cohort_size=7" in n for n in notes)
