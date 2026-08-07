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
another sample: the alternatives do not reproduce a specific tract across independent libraries.
So a sub-threshold call can be corroborated, and the recurrence check above is what stops that
from rescuing an artefact instead.

    gconv_cohort.py --tracts s1.tsv s2.tsv ... --loci s1.loci.tsv ... -o cohort.tsv

The per-locus files are what make the quiet samples visible. A sample that reported nothing at a
locus could have found nothing there or merely have fallen under the reporting threshold, and
those mean opposite things when counting how many samples carry an event.
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict

from gene_conversion import COLUMNS as TRACT_COLUMNS

# Added to every row of the cohort file. The per-sample columns pass through untouched.
COHORT_COLUMNS = ["event_id", "event_samples", "event_frac", "cohort_verdict",
                  "cohort_mismap", "cohort_bf_median"]


def read_tsv(paths):
    rows = []
    for path in paths or []:
        try:
            with open(path) as fh:
                rows.extend(csv.DictReader(fh, delimiter="\t"))
        except OSError:
            continue
    return rows


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


def locus_background(locus_rows):
    """Per (pair_id) cohort statistics from the per-locus files.

    The mismapping rate and the substitution rate are properties of the reference and the
    aligner, not of the sample: the same locus produces the same trouble in everybody. Measured
    once per sample they are noisy; the cohort median is the same quantity with the noise
    averaged out, and it is what says whether one sample's number is unusual.
    """
    by_pair = defaultdict(list)
    for r in locus_rows:
        by_pair[r.get("pair_id", "")].append(r)
    out = {}
    for pair, rs in by_pair.items():
        out[pair] = {
            "n_samples": len({r.get("sample") for r in rs}),
            "mismap": median([num(r, "mismap_frac") for r in rs]),
            "bf": median([num(r, "log10_bf") for r in rs]),
            "mut_rate": median([num(r, "mut_rate") for r in rs]),
        }
    return out


def cohort_verdict(row, event, background, n_samples, min_bf, corroborated_bf,
                   ubiquitous, min_samples):
    """Revise a per-sample verdict in the light of the rest of the cohort.

    Only two revisions are made, in opposite directions, and both need the cohort to be possible.
    Everything else is left exactly as the per-sample caller decided it.
    """
    verdict = (row.get("verdict") or "").strip()
    bf = num(row, "log10_bf")
    frac = event["n_samples"] / n_samples if n_samples else 0.0

    # An event in nearly every sample of the cohort is usually not an event: the reference
    # carries the wrong base there, or the aligner puts the same reads in the same wrong place for
    # everybody. This is the finding a single sample can never make.
    #
    # It is an inference from recurrence, not a proof, and there is a third explanation the data
    # cannot rule out: a clonal cohort really does share an ancestral conversion. The reason says
    # the count so a reader who knows their isolates are related can read it that way.
    if n_samples >= min_samples and frac >= ubiquitous:
        return "reference_artifact", (
            f"present in {event['n_samples']} of {n_samples} samples ({frac:.0%}): at that "
            "recurrence the reference or the aligner explains it more simply than the same "
            "conversion arising in every isolate. In a clonal cohort it may instead be shared "
            "ancestry, which recurrence alone cannot distinguish")

    # A tract the sample alone could not commit to, corroborated by the same tract at the same
    # coordinates in a sample that could. Contamination and index hopping do not reproduce a
    # specific tract across independent libraries; a real minority conversion does.
    if verdict == "ambiguous" and bf is not None and bf >= corroborated_bf \
            and event["n_called"] >= 1:
        return "gene_conversion", (
            f"log10 Bayes factor {bf:.1f} on its own, below the {min_bf:.1f} needed, but the same "
            f"tract is called outright in {event['n_called']} other sample(s) of the cohort")

    return verdict, (row.get("reason") or "").strip()


def annotate(tract_rows, locus_rows, min_bf=3.0, corroborated_bf=2.0,
             ubiquitous=0.9, min_samples=5, slack=0):
    """Add the cohort columns to every tract row, revising the verdict where the cohort speaks."""
    event_of = group_events(tract_rows, slack)
    background = locus_background(locus_rows)

    samples = {r.get("sample") for r in tract_rows} | {r.get("sample") for r in locus_rows}
    samples.discard(None)
    n_samples = len(samples)

    events = defaultdict(lambda: {"samples": set(), "called": set(), "bf": []})
    for i, r in enumerate(tract_rows):
        eid = event_of.get(i)
        if eid is None:
            continue
        e = events[eid]
        e["samples"].add(r.get("sample"))
        e["bf"].append(num(r, "log10_bf"))
        if (r.get("verdict") or "").strip() == "gene_conversion":
            e["called"].add(r.get("sample"))

    out = []
    for i, r in enumerate(tract_rows):
        eid = event_of.get(i)
        e = events.get(eid, {"samples": set(), "called": set(), "bf": []})
        # A sample does not corroborate itself.
        info = {"n_samples": len(e["samples"]),
                "n_called": len(e["called"] - {r.get("sample")})}
        verdict, reason = cohort_verdict(r, info, background, n_samples, min_bf,
                                         corroborated_bf, ubiquitous, min_samples)
        bg = background.get(r.get("pair_id"), {})
        row = dict(r)
        row.update({
            "event_id": eid or "",
            "event_samples": info["n_samples"],
            "event_frac": round(info["n_samples"] / n_samples, 3) if n_samples else "",
            "cohort_verdict": verdict,
            "cohort_mismap": "" if bg.get("mismap") is None else round(bg["mismap"], 4),
            "cohort_bf_median": "" if bg.get("bf") is None else round(bg["bf"], 2),
        })
        if verdict != (r.get("verdict") or "").strip():
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
    p.add_argument("--corroborated-bf", type=float, default=2.0,
                   help="Bayes factor a sub-threshold tract needs before another sample's "
                        "outright call is allowed to corroborate it")
    p.add_argument("--ubiquitous", type=float, default=0.9,
                   help="fraction of the cohort at which an event is called a reference artifact")
    p.add_argument("--min-samples", type=int, default=5,
                   help="cohort size below which recurrence says too little to act on")
    p.add_argument("--slack", type=int, default=0,
                   help="bases of tolerance when deciding two tracts are the same event")
    return p.parse_args(argv)


def main(argv=None) -> int:
    a = parse_args(argv)
    tracts = read_tsv(a.tracts)
    loci = read_tsv(a.loci)

    # Taken from the tool that wrote them when there are rows, and from its declared column list
    # when there are none, so an empty cohort still gets a file the report can parse rather than
    # a header made only of the columns this script adds.
    columns = (list(tracts[0].keys()) if tracts else list(TRACT_COLUMNS)) + COHORT_COLUMNS
    rows, n_samples = annotate(tracts, loci, a.min_bf, a.corroborated_bf,
                               a.ubiquitous, a.min_samples, a.slack)

    with open(a.output, "w") as fh:
        fh.write("\t".join(columns) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")) for c in columns) + "\n")

    changed = sum(1 for r in rows if r["cohort_verdict"] != (r.get("verdict") or "").strip())
    artifacts = sum(1 for r in rows if r["cohort_verdict"] == "reference_artifact")
    rescued = sum(1 for r in rows if r["cohort_verdict"] == "gene_conversion"
                  and (r.get("verdict") or "").strip() != "gene_conversion")
    sys.stderr.write(
        f"[gconv_cohort] {len(rows)} tract row(s) over {n_samples} sample(s), "
        f"{len({r['event_id'] for r in rows})} distinct event(s); {changed} verdict(s) revised "
        f"({artifacts} as reference artifacts, {rescued} corroborated across samples)\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
