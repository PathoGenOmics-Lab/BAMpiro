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
