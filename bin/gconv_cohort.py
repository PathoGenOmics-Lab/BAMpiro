#!/usr/bin/env python3
"""Read a cohort's gene-conversion calls together instead of one sample at a time.

Every sample so far has been judged on its own, and that leaves the two hardest questions
unanswerable. Both are answerable the moment the other samples are in the room.

**Is this the isolate, or is it the reference?**

A tract reported in one sample of fifty is a finding. The same tract at the same coordinates in
all fifty is not a conversion that happened fifty times: it is the reference being wrong there,
or the aligner doing the same thing to everybody. Nothing in a single sample distinguishes those,
and the per-sample caller cannot help: mismapped reads look exactly like a conversion whether the
cause is this sample or the reference. Recurrence is the discriminator, and it needs the cohort.

**Is a weak signal a minority conversion or contamination?**

The per-sample caller refuses to speak below `--gconv_min_tract_af`, because under about a fifth
of the reads a minority conversion, index hopping and a contaminating sample are the same picture.
That floor cost two of the twenty-three tracts in the benchmark, at 25% and 40% frequency. A weak
signal is a different proposition when the SAME tract, at the SAME coordinates, is unambiguous in
another sample: contamination and index hopping do not reproduce a specific tract across
independent libraries. So a sub-threshold call can be corroborated, and the recurrence check above
is what stops that from rescuing an artefact instead. Reads from a third copy of a gene family DO
reproduce a stretch in every library mapped to the same reference, at a few percent of the reads,
so neither a trickle like that nor a call read by a handful of molecules takes part, on either
side (`gene_conversion.unreadable`).

    gconv_cohort.py --tracts s1.tsv s2.tsv ... --loci s1.loci.tsv ... -o cohort.tsv

The per-locus files are what make the quiet samples visible. A sample that reported nothing at a
locus could have found nothing there or merely have fallen under the reporting threshold, and
those mean opposite things when counting how many samples carry an event.
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
from collections import defaultdict

from gene_conversion import COLUMNS as TRACT_COLUMNS
from gene_conversion import unreadable

# Added to every row of the cohort file. The per-sample columns pass through untouched.
COHORT_COLUMNS = ["event_id", "event_samples", "event_frac", "cohort_verdict",
                  "cohort_mismap", "cohort_bf_median",
                  "donor_rank", "n_donors", "donor_margin", "donor_call", "is_representative"]


def _file_unit(path):
    """The sample and reference a per-sample file belongs to, from its name
    (`<sample>.<ref>.gene_conversion.tsv` and `..._loci.tsv`), or None for another name."""
    base = os.path.basename(str(path))
    for suffix in (".gene_conversion_loci.tsv", ".gene_conversion.tsv"):
        if base.endswith(suffix):
            return base[:-len(suffix)]
    return None


def unit_of(row):
    """What a divergence is judged on: the sample on one reference. A sample mapped against two
    references is two genomes compared with two references, and one can be diverged from its
    reference while the other is not."""
    return row.get("_unit") or row.get("sample")


def read_tsv(paths):
    """Rows, and the `#` provenance lines the per-sample files carry.

    Those lines say which settings produced the numbers being combined here. Dropping them would
    leave a cohort file that cannot be reproduced or compared with another run, and a cohort
    assembled from samples run with DIFFERENT settings is a thing a reader has to be able to see.
    Each row remembers the sample and reference its file was written for (`_unit`), which the
    output does not carry.
    """
    rows, notes = [], []
    for path in paths or []:
        unit = _file_unit(path)
        try:
            with open(path) as fh:
                for line in fh:
                    if line.startswith("#"):
                        note = line.rstrip("\n")
                        if note not in notes:
                            notes.append(note)
                    else:
                        break
                fh.seek(0)
                for r in csv.DictReader((ln for ln in fh if not ln.startswith("#")), delimiter="\t"):
                    if unit:
                        r["_unit"] = unit
                    rows.append(r)
        except OSError:
            continue
    return rows, notes


def num(row, key):
    v = (row.get(key) or "").strip()
    try:
        return float(v)
    except ValueError:
        return None


def median(values):
    v = sorted(x for x in values if x is not None)
    if not v:
        return None
    mid = len(v) // 2
    return v[mid] if len(v) % 2 else (v[mid - 1] + v[mid]) / 2.0


def group_events(rows, slack=0):
    """Group tract rows into events: the same stretch of the same acceptor, across samples.

    Grouped on the acceptor's COORDINATES rather than on the paralog pair, so the several
    relationships a gene family reports the same converted stretch through collapse into one
    event. That was measured at 1.2 rows per event on the real genome, up to 3.
    """
    by_contig = defaultdict(list)
    for i, r in enumerate(rows):
        try:
            start, end = int(r["start"]), int(r["end"])
        except (KeyError, ValueError):
            continue
        by_contig[r.get("contig", "")].append((start, end, i))

    event_of = {}
    next_id = 0
    for contig in sorted(by_contig):
        spans = sorted(by_contig[contig])
        cur_end, cur_id = None, None
        for start, end, i in spans:
            if cur_end is None or start > cur_end + slack:
                cur_id, next_id = next_id, next_id + 1
                cur_end = end
            else:
                cur_end = max(cur_end, end)
            event_of[i] = f"{contig}:{cur_id}" if contig else str(cur_id)
    return event_of


def rank_donors(rows, indices, margin=1.0):
    """Which relative the converted stretch actually came from, where the reads can say.

    A gene family reports one event once per relationship, and those rows are two different
    things wearing the same shape. Sometimes they point at the SAME stretch of donor, in which
    case there is no choice to make and the extra rows are redundancy: overlapping self-alignments
    of a tandem repeat do this, and their evidence comes out tied to two decimal places. Sometimes
    they point at genuinely different places in the genome, and then one of them is the source.

    The evidence is `log10_bf_vs_null` per marker: how well that relationship's donor sequence
    accounts for what the reads carry, divided by how many markers it had to work with. Measured on
    the real genome it separates a true source from a bystander by more than an order of magnitude
    per site (19.1 against 1.4). It also honestly fails to separate two relatives that fit equally
    well, which is why `donor_margin` is reported rather than a bare winner: below `margin` the
    tract is compatible with either, and short reads do not carry what would settle it.
    """
    # Candidates are grouped by where in the DONOR the copied stretch sits, so the several
    # relationships that name the same stretch count once.
    cands = []
    for i in indices:
        r = rows[i]
        try:
            span = (int(r["don_start"]), int(r["don_end"]))
        except (KeyError, TypeError, ValueError):
            span = None
        per = None
        n = num(r, "n_sites")
        bf = num(r, "log10_bf_vs_null")
        if bf is not None and n:
            per = bf / n
        placed = False
        for c in cands:
            if span and c["span"] and not (span[1] < c["span"][0] or span[0] > c["span"][1]):
                c["rows"].append(i)
                c["span"] = (min(span[0], c["span"][0]), max(span[1], c["span"][1]))
                if per is not None and (c["per"] is None or per > c["per"]):
                    c["per"] = per
                placed = True
                break
        if not placed:
            cands.append({"span": span, "rows": [i], "per": per})

    cands.sort(key=lambda c: (c["per"] is None, -(c["per"] or 0.0)))
    gap = None
    if len(cands) > 1 and cands[0]["per"] is not None and cands[1]["per"] is not None:
        gap = cands[0]["per"] - cands[1]["per"]

    out = {}
    for rank, c in enumerate(cands, start=1):
        # Inside one candidate the rows describe the same source, so the best-supported of them
        # represents it and the rest are the family saying it again.
        best = max(c["rows"], key=lambda i: (num(rows[i], "log10_bf_vs_null") or -1e9))
        for i in c["rows"]:
            out[i] = {"donor_rank": rank, "n_donors": len(cands),
                      "donor_margin": "" if gap is None else round(gap, 2),
                      "is_representative": int(rank == 1 and i == best),
                      "donor_call": ("resolved" if rank == 1 and gap is not None and gap >= margin
                                     else "only candidate" if len(cands) == 1
                                     else "ambiguous")}
    return out


def locus_key(row):
    """A paralog pair is named by its contig as well as its number.

    `pair_id` numbers the pairs of ONE reference's paralog map. A cohort mapped against two
    references has a pair 61 in each, and they are different loci: keyed on the number alone, the
    medians below mixed a locus of one genome with an unrelated locus of the other.
    """
    return (row.get("contig") or "", row.get("pair_id") or "")


def readable_call(row):
    """A per-sample conversion call whose reads can carry it (see gene_conversion.unreadable).

    Checked here as well as where the call is made, so a cohort pass over per-sample files written
    before that check still refuses to count a call resting on a handful of reads.
    """
    return ((row.get("verdict") or "").strip() == "gene_conversion"
            and not unreadable(row.get("n_sites"), row.get("n_undetermined"), row.get("donor_af_in")))


def locus_background(locus_rows):
    """Per locus (contig and pair) cohort statistics from the per-locus files.

    The mismapping rate and the substitution rate are properties of the reference and the
    aligner, not of the sample: the same locus produces the same trouble in everybody. Measured
    once per sample they are noisy; the cohort median is the same quantity with the noise
    averaged out, and it is what says whether one sample's number is unusual.
    """
    by_pair = defaultdict(list)
    for r in locus_rows:
        by_pair[locus_key(r)].append(r)
    out = {}
    for pair, rs in by_pair.items():
        out[pair] = {
            "n_samples": len({r.get("sample") for r in rs}),
            "mismap": median([num(r, "mismap_frac") for r in rs]),
            "bf": median([num(r, "log10_bf") for r in rs]),
            "mut_rate": median([num(r, "mut_rate") for r in rs]),
        }
    return out


def divergent_samples(tract_rows, max_locus_frac, min_loci=20, min_sample_loci=5, event_of=None,
                      locus_rows=()):
    """Samples calling tracts in so many loci at once that the reference is the problem.

    Gene conversion is a local event. A sample converting a quarter of every paralogous locus it
    has is not converting: it is a different genotype from the one the reads were mapped to, and
    at every paralogous locus it carries the donor's base by inheritance, which is exactly the
    pattern the model scores as a tract.

    This is not hypothetical. In a 185-sample cohort, 22 samples labelled as one lineage were
    really another, were therefore mapped against a reference 2,000 SNPs away, and produced 507 of
    the 522 tracts. Tracts per sample against genome-wide SNP count correlated at r = 0.98. The
    separation is not subtle: those samples called tracts in 16 to 29 percent of the stretches
    their reference reports anything in, and every other sample in at most 1 percent.

    Counted in distinct STRETCHES of the acceptor (the events of `event_of`), not in rows and not
    in pairs. One stretch scored against several donors yields a row per donor, and a gene family
    reports one stretch through a pair per relative: counted in pairs, a single conversion seen
    through eleven relatives put a clean sample exactly on the threshold. Rows without
    coordinates fall back to their pair. The fraction is taken over the stretches of every contig
    of the sample's own reference, which its per-locus rows name: counted over only the contigs
    it has tracts on, five tracts on a plasmid made a sample with a clean chromosome divergent.

    Keyed on the sample and its reference (see unit_of), not the sample alone.
    """
    event_of = event_of or {}
    units, contigs_of, by_sample = defaultdict(set), defaultdict(set), {}
    for r in locus_rows or ():
        contigs_of[unit_of(r)].add(r.get("contig") or "")
    for i, r in enumerate(tract_rows):
        pair = r.get("pair_id")
        if pair in (None, ""):
            continue
        contig = r.get("contig") or ""
        unit = event_of.get(i) or ("pair", contig, pair)
        units[contig].add(unit)
        contigs_of[unit_of(r)].add(contig)
        if readable_call(r):
            by_sample.setdefault(unit_of(r), set()).add(unit)
    # A fraction of a handful of loci says nothing: one tract out of three is 33% and means only
    # that the reference has three paralogous loci. Both floors are absolute for that reason, and
    # a cohort below them is left alone rather than judged on a ratio it cannot support.
    out = {}
    for s, called in by_sample.items():
        n = sum(len(units[c]) for c in contigs_of[s])
        if n >= min_loci and len(called) >= min_sample_loci and len(called) / n > max_locus_frac:
            out[s] = len(called) / n
    return out


def cohort_verdict(row, event, background, n_samples, min_bf, corroborated_bf,
                   ubiquitous, min_samples, divergent=None, min_tract_af=0.25, scoped=False):
    """Revise a per-sample verdict in the light of the rest of the cohort.

    Only two revisions are made, in opposite directions, and both need the cohort to be possible.
    Everything else is left exactly as the per-sample caller decided it, except a call its own
    reads cannot carry, which per-sample files written before that check still hold.

    `n_samples` is the cohort the event could have appeared in: the samples mapped to the same
    reference (`scoped`), or the whole cohort when the per-locus files cannot say which those are.
    """
    verdict = (row.get("verdict") or "").strip()
    bf = num(row, "log10_bf")
    why_not = unreadable(row.get("n_sites"), row.get("n_undetermined"), row.get("donor_af_in"))

    # Before anything else: if the SAMPLE is diverged from the reference, none of its tracts mean
    # what they say, however good each one looks on its own. This one comes first because the
    # corroboration rule below would otherwise let these samples vouch for each other: they share
    # a genotype, so they share their artefacts at identical coordinates, which reads as exactly
    # the independent recurrence that rule is looking for.
    frac_div = (divergent or {}).get(unit_of(row))
    if frac_div is not None:
        return "divergent_sample", (
            f"this sample calls tracts in {frac_div:.0%} of its paralogous loci. Conversion is "
            "local; a sample converting that many loci at once is a different genotype from the "
            "reference it was mapped to, and carries the donor base at each of them by "
            "inheritance. Check the sample's lineage against the reference before reading this")
    # None means the cohort size is unknown, which is not the same as small. A fraction cannot be
    # formed at all, so the recurrence rule below is skipped rather than fed a stand-in.
    frac = event["n_samples"] / n_samples if n_samples else 0.0

    # An event in nearly every sample of the cohort is usually not an event: the reference
    # carries the wrong base there, or the aligner puts the same reads in the same wrong place for
    # everybody. This is the finding a single sample can never make.
    #
    # It is an inference from recurrence, not a proof, and there is a third explanation the data
    # cannot rule out: a clonal cohort really does share an ancestral conversion. The reason says
    # the count so a reader who knows their isolates are related can read it that way.
    if n_samples is not None and n_samples >= min_samples and frac >= ubiquitous:
        return "reference_artifact", (
            f"present in {event['n_samples']} of {n_samples} samples"
            f"{' mapped to this reference' if scoped else ''} ({frac:.0%}): at that "
            "recurrence the reference or the aligner explains it more simply than the same "
            "conversion arising in every isolate. In a clonal cohort it may instead be shared "
            "ancestry, which recurrence alone cannot distinguish")

    # Both verdicts: an ambiguous row resting on a trickle still carries the reason the model
    # gave, "only 20% of the reads here carry the tract", which quotes the floor the fit was
    # parked on rather than the 3% the reads actually hold.
    # The same verdicts the per-sample caller demotes, reciprocal_exchange included, so re-running
    # only this step corrects files written before the check. A reason that already starts with
    # the check is such a file's own, with its copy-number and breakpoint notes after it: kept.
    # Compared up to "read by", since this pass does not know --gconv_min_depth and writes
    # "enough molecules" where the caller wrote "5 or more molecules".
    if why_not and verdict in ("gene_conversion", "reciprocal_exchange", "ambiguous"):
        reason = (row.get("reason") or "").strip()
        return "ambiguous", reason if reason.startswith(why_not.split(" read by ")[0]) else why_not

    # A tract the sample alone could not commit to, corroborated by the same tract at the same
    # coordinates in a sample that could. Contamination and index hopping do not reproduce a
    # specific tract across independent libraries; a real minority conversion does.
    #
    # Reads from a third copy of the family DO reproduce one, in every library mapped to the same
    # reference, because they are a property of the reference and the aligner. That is why a row
    # its own reads cannot carry is never promoted: on a 185-sample cohort, 1,520 of 1,579 calls
    # were rows of 2 to 5% donor reads that a single thinly read sample had "called outright".
    if verdict == "ambiguous" and not why_not and bf is not None and bf >= corroborated_bf \
            and event["n_called"] >= 1:
        taf = num(row, "tract_af")
        short = []
        if bf < min_bf:
            short.append(f"log10 Bayes factor {bf:.1f} on its own, below the {min_bf:.1f} needed")
        if taf is not None and taf < min_tract_af:
            short.append(f"carried by {taf:.0%} of the reads, under the {min_tract_af:.0%} one "
                         "sample needs on its own")
        return "gene_conversion", (
            ("; ".join(short) or "not called on its own") +
            f", but the same tract is called outright in {event['n_called']} other sample(s) of "
            "the cohort")

    return verdict, (row.get("reason") or "").strip()


def annotate(tract_rows, locus_rows, min_bf=3.0, corroborated_bf=2.0,
             ubiquitous=0.9, min_samples=5, slack=0, donor_margin=1.0, cohort_size=None,
             max_locus_frac=0.1, min_tract_af=0.25):
    """Add the cohort columns to every tract row, revising the verdict where the cohort speaks.

    `cohort_size` is how many samples were RUN. It matters because the ubiquity rule is a
    fraction, and the samples a tract file can name are only the ones that had something to
    report: a sample whose genome is clean writes no tract row and is invisible to the numerator
    and the denominator alike. Counting only those leaves an event in 5 of 8 samples looking like
    an event in 5 of 5, which demotes a real conversion to a reference artifact, the exact
    inversion of what this pass exists to do.

    The per-locus rows are the census that fixes it, since every sample writes one per pair it
    evaluated whether or not anything came of it. Without them, and without an explicit size, the
    cohort is genuinely unknown and the rule is left unapplied rather than applied to a number
    that only looks like an answer.

    The census is kept per contig. A cohort mapped against two references splits into the samples
    of each, and an event on one reference can only ever appear in that reference's samples:
    measured against the whole cohort, an artefact of a reference carrying 113 of 185 samples
    could never reach more than 61% and was never demoted.

    Samples diverged from their reference count towards neither side of the fraction. They carry
    the donor's base at every paralogous locus by inheritance, so they would add themselves to
    every event on the reference and nothing about how often the reference misleads.
    """
    event_of = group_events(tract_rows, slack)
    background = locus_background(locus_rows)
    divergent = divergent_samples(tract_rows, max_locus_frac, event_of=event_of, locus_rows=locus_rows)
    # The samples a divergence was found in, and on which contigs, since a sample on two
    # references counts as diverged only on the one it is diverged from.
    div_on = {(r.get("sample"), r.get("contig") or "") for r in list(tract_rows) + list(locus_rows)
              if unit_of(r) in divergent}

    samples = {r.get("sample") for r in tract_rows} | {r.get("sample") for r in locus_rows}
    samples.discard(None)
    if cohort_size is not None:
        n_samples = max(int(cohort_size), len(samples))
    elif locus_rows:
        n_samples = len(samples)
    else:
        n_samples = None

    census = defaultdict(set)
    for r in locus_rows:
        census[r.get("contig") or ""].add(r.get("sample"))
    for r in tract_rows:                         # a sample with a tract on a contig was run on it
        contig = r.get("contig") or ""
        if contig in census:
            census[contig].add(r.get("sample"))
    unscoped = None if n_samples is None else n_samples - len({s for s, _ in div_on} & samples)

    def cohort_of(contig):
        """(number of samples the event could have appeared in, whether that is per reference)"""
        if contig in census:
            return len({s for s in census[contig] if (s, contig) not in div_on}), True
        return unscoped, False

    events = defaultdict(lambda: {"samples": set(), "called": set(), "bf": []})
    for i, r in enumerate(tract_rows):
        eid = event_of.get(i)
        if eid is None or unit_of(r) in divergent:
            continue
        e = events[eid]
        e["samples"].add(r.get("sample"))
        e["bf"].append(num(r, "log10_bf"))
        if readable_call(r):
            e["called"].add(r.get("sample"))

    # Donor resolution is per event PER SAMPLE: the several relationships are one sample's
    # several views of one stretch, and two samples carrying the same event each get their own
    # reading of where it came from.
    by_event_sample = defaultdict(list)
    for i, r in enumerate(tract_rows):
        if event_of.get(i) is not None:
            by_event_sample[(event_of[i], r.get("sample"))].append(i)
    donors = {}
    for idxs in by_event_sample.values():
        donors.update(rank_donors(tract_rows, idxs, donor_margin))

    out = []
    for i, r in enumerate(tract_rows):
        eid = event_of.get(i)
        e = events.get(eid, {"samples": set(), "called": set(), "bf": []})
        # A sample does not corroborate itself.
        info = {"n_samples": len(e["samples"]),
                "n_called": len(e["called"] - {r.get("sample")})}
        n_here, scoped = cohort_of(r.get("contig") or "")
        verdict, reason = cohort_verdict(r, info, background, n_here, min_bf,
                                         corroborated_bf, ubiquitous, min_samples, divergent,
                                         min_tract_af=min_tract_af, scoped=scoped)
        bg = background.get(locus_key(r), {})
        row = dict(r)
        row.update({
            "event_id": eid or "",
            "event_samples": info["n_samples"],
            "event_frac": round(info["n_samples"] / n_here, 3) if n_here else "",
            "cohort_verdict": verdict,
            "cohort_mismap": "" if bg.get("mismap") is None else round(bg["mismap"], 4),
            "cohort_bf_median": "" if bg.get("bf") is None else round(bg["bf"], 2),
        })
        row.update(donors.get(i, {"donor_rank": "", "n_donors": "", "donor_margin": "",
                                  "donor_call": "", "is_representative": ""}))
        if verdict != (r.get("verdict") or "").strip() or reason != (r.get("reason") or "").strip():
            row["reason"] = reason
        out.append(row)
    return out, n_samples


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--tracts", nargs="*", default=[], help="per-sample gene_conversion TSVs")
    p.add_argument("--loci", nargs="*", default=[],
                   help="per-sample per-locus TSVs, so the quiet samples are visible too")
    p.add_argument("-o", "--output", required=True)
    p.add_argument("--min-bf", type=float, default=3.0,
                   help="the per-sample calling threshold, quoted back in the reasons")
    p.add_argument("--min-tract-af", type=float, default=0.25,
                   help="the per-sample read-fraction floor, quoted back in the reasons")
    p.add_argument("--corroborated-bf", type=float, default=2.0,
                   help="Bayes factor a sub-threshold tract needs before another sample's "
                        "outright call is allowed to corroborate it")
    p.add_argument("--max-locus-frac", type=float, default=0.1,
                   help="A sample calling tracts in more than this FRACTION of its paralogous "
                        "loci is diverged from the reference rather than converting, and all of "
                        "its tracts are marked divergent_sample. Counted in distinct stretches of "
                        "its reference. Measured on a real cohort the two groups do not overlap: "
                        "mislabelled samples called tracts in 16-29%% of the stretches and every "
                        "correctly mapped sample in at most 1%%.")
    p.add_argument("--ubiquitous", type=float, default=0.9,
                   help="fraction of the samples mapped to the same reference at which an event "
                        "is called a reference artifact")
    p.add_argument("--min-samples", type=int, default=5,
                   help="cohort size below which recurrence says too little to act on")
    p.add_argument("--donor-margin", type=float, default=1.0,
                   help="evidence per marker, in log10, by which the best donor must beat the "
                        "next distinct one before the source is called resolved")
    p.add_argument("--slack", type=int, default=0,
                   help="bases of tolerance when deciding two tracts are the same event")
    p.add_argument("--cohort-size", type=int, default=None,
                   help="how many samples were run, when --loci is not available to say so. "
                        "It is counted against every reference alike, which is right only when "
                        "the run has one; without either, recurrence has no denominator and is "
                        "left alone")
    return p.parse_args(argv)


def main(argv=None) -> int:
    a = parse_args(argv)
    tracts, notes = read_tsv(a.tracts)
    loci, _ = read_tsv(a.loci)

    # Taken from the tool that wrote them when there are rows, and from its declared column list
    # when there are none, so an empty cohort still gets a file the report can parse rather than
    # a header made only of the columns this script adds.
    columns = (list(tracts[0].keys()) if tracts else list(TRACT_COLUMNS)) + COHORT_COLUMNS
    rows, n_samples = annotate(tracts, loci, a.min_bf, a.corroborated_bf,
                               a.ubiquitous, a.min_samples, a.slack, a.donor_margin,
                               a.cohort_size, max_locus_frac=a.max_locus_frac,
                               min_tract_af=a.min_tract_af)

    if n_samples is None and tracts:
        sys.stderr.write(
            "[gconv_cohort] no --loci and no --cohort-size, so the number of samples run is not "
            "known: a sample with nothing to report writes no tract row and cannot be counted. "
            "Recurrence is left unapplied, and no event will be called a reference artifact\n")

    with open(a.output, "w") as fh:
        for note in notes:
            fh.write(note + "\n")
        fh.write(f"# gconv_cohort.py min_bf={a.min_bf} min_tract_af={a.min_tract_af} "
                 f"corroborated_bf={a.corroborated_bf} "
                 f"ubiquitous={a.ubiquitous} min_samples={a.min_samples} "
                 f"donor_margin={a.donor_margin} slack={a.slack} "
                 f"cohort_size={n_samples if n_samples is not None else 'unknown'}\n")
        fh.write("\t".join(columns) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")) for c in columns) + "\n")

    changed = sum(1 for r in rows if r["cohort_verdict"] != (r.get("verdict") or "").strip())
    artifacts = sum(1 for r in rows if r["cohort_verdict"] == "reference_artifact")
    rescued = sum(1 for r in rows if r["cohort_verdict"] == "gene_conversion"
                  and (r.get("verdict") or "").strip() != "gene_conversion")
    reps = sum(1 for r in rows if r.get("is_representative") == 1)
    sys.stderr.write(
        f"[gconv_cohort] {len(rows)} tract row(s) collapsing to {reps} event/donor call(s) "
        f"over {n_samples if n_samples is not None else 'an unknown number of'} sample(s), "
        f"{len({r['event_id'] for r in rows})} distinct event(s); {changed} verdict(s) revised "
        f"({artifacts} as reference artifacts, {rescued} corroborated across samples)\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
