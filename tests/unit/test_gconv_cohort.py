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
          bf="12.0", **extra):
    row = {"sample": sample, "pair_id": pair, "contig": contig, "donor": "chr",
           "verdict": verdict, "reason": "because", "start": str(start), "end": str(end),
           "log10_bf": bf, "tract_af": "0.9", "mismap_frac": "0.01"}
    row.update(extra)
    return row


def locus(sample, pair="1", mismap="0.01", bf="1.0"):
    return {"sample": sample, "pair_id": pair, "contig": "chr", "donor": "chr",
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
