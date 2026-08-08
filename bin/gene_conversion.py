#!/usr/bin/env python3
"""Find gene conversion tracts from short reads.

Gene conversion copies a stretch of one paralog onto another without reciprocal exchange, so the
acceptor locus stops carrying its own alleles and starts carrying the donor's over a bounded
tract. On a reference that has never seen the event, that reads as a run of "variants" which are
not random: they are exactly the donor's bases, in order, between two sharp breakpoints.

    gene_conversion.py --sites sites.tsv --bam sample.bam --sample S1 -o tracts.tsv

`sites.tsv` comes from paralog_map.py: the positions where the two copies differ. Everywhere else
the copies are identical and a read is uninformative by construction, so only those are examined.

THE HARD PART IS NOT FINDING TRACTS, IT IS NOT BELIEVING THEM.

There are three ways an acceptor site can show the donor's base, and only one of them is a
conversion:

* it was converted, along with its neighbours, as one tract;
* it mutated to that base on its own, which for a single site is entirely ordinary;
* the read carrying it came from the donor and mismapped, which in a paralogous region is what
  aligners do.

Deciding between them is `gconv_model`, which evaluates every possible tract against the other
two explanations, weighting each base by its quality and treating each MOLECULE as the unit of
evidence rather than each site. What comes out is a Bayes factor and a posterior over the
breakpoints, not a threshold crossing.

What this file adds around it is the descriptive statistics a person needs to check the call
against the BAM: how fixed the donor allele is inside the tract (`donor_af_in`), whether donor
alleles also turn up outside it (`donor_af_outside`), and how many single molecules carry donor
alleles on one side of a breakpoint and acceptor alleles on the other, in cis (`breakpoint_reads`,
the one thing a mismapped read cannot fake).

MAPQ is deliberately NOT filtered on. Paralogous regions are exactly where an aligner assigns
MAPQ 0, so the usual quality gate would discard the entire signal this tool exists to find.
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from collections import defaultdict

import gconv_annotate as ga
import gconv_model as gm

# CIGAR operations that consume the reference, the read, or both.
CONSUMES_REF = set("MDN=X")
CONSUMES_READ = set("MIS=X")

# SAM lets the quality string be absent ("*"). There is then no evidence about how reliable a
# base is, so it passes any quality floor, and the model caps it at its own ceiling rather than
# taking silence for certainty.
NO_QUAL_PHRED = 60

# A diagnostic site where the DONOR is missing bases the acceptor has. The observation is not
# which base a read carries but whether it carries one at all: a read from the unconverted
# acceptor has a base there, a read from the converted acceptor or from the donor spans it with a
# deletion. Written into the site table as the donor's "base" so one code path serves both kinds.
GAP = "-"

# A deletion has no per-base quality to weigh it by, so one is assigned. Deliberately below what a
# good base is worth: short-read aligners place indels less reliably than substitutions, and this
# is an alignment-level call rather than a base call.
INDEL_PHRED = 20


def parse_cigar(cigar):
    """CIGAR string to [(length, op), ...]."""
    ops, num = [], ""
    for ch in cigar:
        if ch.isdigit():
            num += ch
        elif num:
            ops.append((int(num), ch))
            num = ""
    return ops


def read_calls_at(pos, cigar, seq, qual, wanted):
    """Base and Phred score this read carries at the reference positions in `wanted`.

    Walks the CIGAR so an insertion or a soft clip does not shift every base after it, which is
    the classic way to read the right position off the wrong read.

    Returns {ref_pos (1-based): (base, phred)}. A position inside a deletion, or past the end of
    a truncated quality string, is omitted rather than guessed at. Nothing is filtered here: the
    model weights each base by its quality instead of discarding it, so the score travels with
    the base.
    """
    out = {}
    ref, qry = pos, 0
    for length, op in parse_cigar(cigar):
        if op in CONSUMES_REF and op in CONSUMES_READ:          # M, =, X
            for i in range(length):
                p = ref + i
                if p in wanted:
                    q = qry + i
                    if q >= len(seq):
                        continue
                    if not qual or qual == "*":
                        phred = NO_QUAL_PHRED
                    elif q < len(qual):
                        phred = ord(qual[q]) - 33
                    else:
                        continue                                 # quality string ran out
                    out[p] = (seq[q].upper(), phred)
            ref += length
            qry += length
        elif op == "D":
            # No base to read, and that IS the reading: at a site where the donor is missing
            # bases, a deletion here is the donor's allele. N is left out on purpose, being a
            # skipped region rather than a claim that the bases are absent.
            for i in range(length):
                if ref + i in wanted:
                    out[ref + i] = (GAP, INDEL_PHRED)
            ref += length
        elif op in CONSUMES_REF:                                 # N: no base, and no claim
            ref += length
        elif op in CONSUMES_READ:                                # I, S: consumes no reference
            qry += length
    return out


def read_bases_at(pos, cigar, seq, qual, wanted, min_bq=13):
    """Bases this read carries at `wanted`, with everything below the quality floor dropped.

    Deletions are not bases and are left out, so this keeps the meaning it always had.
    """
    return {p: base for p, (base, phred) in read_calls_at(pos, cigar, seq, qual, wanted).items()
            if phred >= min_bq and base != GAP}


def collect_read_observations(sam_lines, sites):
    """Per-read base and Phred score at the diagnostic sites.

    Keyed by read NAME so the two mates of a pair contribute to one observation of one molecule,
    which is what the model treats as the unit of evidence.
    """
    wanted = set(sites)
    obs = defaultdict(dict)
    for line in sam_lines:
        if not line or line.startswith("@"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 11:
            continue
        flag = int(f[1])
        if flag & 0x904:                                        # unmapped, secondary, supplementary
            continue
        name, pos, cigar, seq, qual = f[0], int(f[3]), f[5], f[9], f[10]
        if cigar == "*" or pos <= 0:
            continue
        obs[name].update(read_calls_at(pos, cigar, seq, qual, wanted))
    return obs


def allele_view(observations, sites, min_bq=13):
    """Which copy each observed base belongs to: {read: {pos: 'acceptor'|'donor'|'other'}}.

    The categorical view is what the descriptive columns are built from. The model does not use
    it: it works off the log-likelihood ratios, where a third base contributes nothing instead of
    being counted as depth and dropped from the numerator.
    """
    per_read = defaultdict(dict)
    for name, calls in observations.items():
        for p, (base, phred) in calls.items():
            if phred < min_bq:
                continue
            acc, don = sites[p]
            if don == GAP:                    # presence or absence, not which base
                per_read[name][p] = "donor" if base == GAP else "acceptor"
            else:
                per_read[name][p] = ("acceptor" if base == acc else
                                     "donor" if base == don else "other")
    return per_read


def collect_read_alleles(sam_lines, sites, min_bq=13):
    """Per-read allele calls at the diagnostic sites, straight from SAM records."""
    return allele_view(collect_read_observations(sam_lines, sites), sites, min_bq)


def pileup(per_read, positions):
    """Per-site counts across reads: how many molecules carry each copy's allele."""
    counts = {p: {"acceptor": 0, "donor": 0, "other": 0} for p in positions}
    for calls in per_read.values():
        for p, which in calls.items():
            if p in counts:
                counts[p][which] += 1
    for p, c in counts.items():
        informative = c["acceptor"] + c["donor"]
        c["depth"] = informative + c["other"]
        c["donor_af"] = (c["donor"] / informative) if informative else None
    return counts


def locus_copies(counts, positions, genome_depth, min_depth=5):
    """How many copies of this locus the reads look like: (ratio, copies, expected tract AF).

    A paralog family with more members than the reference pair puts the extra copies' reads here
    too, and they carry the acceptor's alleles because they were never converted. A clonal
    conversion of ONE copy then shows up in one copy's worth of the reads, not in all of them,
    and the read fraction that would otherwise read as a minority event is exactly what a clonal
    one looks like at that copy number.

    Measured against the sample's own genome-wide depth rather than an absolute, because half the
    reads of a shallow run and half the reads of a deep one are the same fraction and not the
    same thing.

    (None, None, None) without a genome depth to compare against. The number is not guessed at:
    an assumed copy number of 1 is exactly the assumption that makes a clonal conversion in a
    three-copy family look like contamination.
    """
    if not genome_depth or genome_depth <= 0:
        return None, None, None
    seen = [counts[p]["depth"] for p in positions if counts[p]["depth"] >= min_depth]
    if not seen:
        return None, None, None
    ratio = _median(seen) / genome_depth
    copies = max(1, int(round(ratio)))
    return round(ratio, 2), copies, round(1.0 / copies, 3)


def dilution_notes(tract_af, min_tract_af, cn_ratio, copies, expected_af, het_fraction):
    """What else fits a tract carried by less than `min_tract_af` of the reads.

    Copy number EXPLAINS a diluted fraction, it does not discriminate one: at a locus with five
    copies, a clonal conversion of one and a contamination at 20% are the same fraction and no
    amount of depth separates them. A mixed sample is weaker still, saying only that another
    strain is present somewhere, not that it is present here. So neither moves the verdict, and
    both are said, which is the difference between a bucket and a reading.

    A function rather than the condition inline, because the only view of `tract_af` from outside
    is a column rounded to three decimals, where a value on the floor and one just under it read
    the same. That is the second threshold in this file to have hidden behind its own rounding.
    """
    if tract_af is None or tract_af >= min_tract_af:
        return ""
    notes = []
    if copies and copies > 1:
        notes.append(
            f"the locus runs at {cn_ratio} times the genome's depth, so it looks like {copies} "
            f"copies and a clonal conversion of one would reach about {expected_af:.0%} of the "
            "reads")
    if het_fraction is not None:
        notes.append(
            f"{het_fraction:.0%} of this sample's variant sites are heterozygous genome-wide, so "
            "a diluted tract here may belong to another strain. That is context, not a reading "
            "of this locus")
    return "".join("; " + n for n in notes)


def donor_side(args, loc, positions):
    """What the DONOR copy carries: (depth per site, fraction of reads carrying the ACCEPTOR's).

    Both keyed by the ACCEPTOR position, so a caller can line them up with everything else.

    The depth is what catches a tract whose reads moved next door. The alleles are what say
    whether this was a conversion at all. Gene conversion is NON-reciprocal by definition: the
    donor hands over a copy of its sequence and keeps its own. If the donor has also taken on the
    acceptor's bases over the same stretch, the two copies swapped, and a swap is an unequal
    crossover, not a conversion.

    The site pair is given to the reader the other way round on purpose. At the donor's own
    position, the donor's base is the resident one and the acceptor's is the foreign one, so
    `donor_af` there reads as "how much of the donor now looks like the acceptor".
    """
    don_pos = loc.get("donor_pos", {})
    if not don_pos:
        return {}, {}
    dvals = sorted(don_pos.values())
    proc = subprocess.run(
        [args.samtools, "view", args.bam, f"{loc['donor']}:{dvals[0]}-{dvals[-1]}"],
        capture_output=True, text=True)
    if proc.returncode != 0:
        return {}, {}

    # A deletion marker cannot be read from the donor's side: the donor is the copy that HAS no
    # base there, so there is nothing at that position for a read to carry either way.
    sites = {}
    for acc_pos, dp in don_pos.items():
        acc_base, don_base = loc["sites"][acc_pos]
        if don_base != GAP:
            sites[dp] = (don_base, acc_base)
    if not sites:
        return {}, {}

    reads = collect_read_alleles(proc.stdout.splitlines(), sites, args.min_bq)
    counts = pileup(reads, sorted(sites))
    depth = {ap: counts[dp]["depth"] for ap, dp in don_pos.items() if dp in counts}
    swapped = {ap: counts[dp]["donor_af"] for ap, dp in don_pos.items()
               if dp in counts and counts[dp]["donor_af"] is not None}
    return depth, swapped


def reciprocity(tract, positions, donor_swapped, min_sites=2):
    """How much of the donor carries the acceptor's bases INSIDE the tract and outside it.

    Both, because the fraction inside on its own says nothing. Reads from the unconverted
    acceptor that the aligner placed at the donor carry the acceptor's bases at every site it
    reaches, which is the same picture at a single site as a donor that genuinely swapped. What
    tells them apart is the one thing this whole tool leans on: a real event is BOUNDED. An
    exchange shows the acceptor's bases over the tract and the donor's own outside it; reads that
    arrived from next door show the acceptor's everywhere.

    (None, None) when too few sites on either side could be read. Silence is not evidence that
    the donor stayed put, and it is not evidence that it moved.
    """
    inside = [donor_swapped[p] for p in tract if p in donor_swapped]
    tract_set = set(tract)
    outside = [donor_swapped[p] for p in positions
               if p not in tract_set and p in donor_swapped]
    if len(inside) < min_sites or len(outside) < min_sites:
        return None, None
    return sum(inside) / len(inside), sum(outside) / len(outside)


def is_exchange(inside, outside, reciprocal_af):
    """Whether the donor swapped with the acceptor over this stretch and only over it.

    The second half is what keeps wholesale mismapping out. A donor locus whose reads came from
    the acceptor carries the acceptor's bases outside the tract too, and calling that an exchange
    both invents a finding and buries the conversion that is really there.
    """
    if inside is None or outside is None:
        return False
    return inside >= reciprocal_af and outside < reciprocal_af


def worth_reporting(fit, report_bf, min_mismap):
    """Whether a fit is written out at all.

    Either it clears the reporting threshold, or the locus is explained by reads arriving from
    the paralog, which is a result too and not a silence: a reader looking for a conversion there
    needs to be told the reads plainly came from next door.

    A function rather than the condition inline, because both thresholds are settings a user
    picks and the only way their edges were visible from outside was through a column rounded to
    two decimals, where a value on the threshold and a value just past it read the same.
    """
    return fit["log10_bf"] >= report_bf or fit["mismap_frac"] >= min_mismap


def tract_evidence(tract, positions, counts, per_read, in_any_tract=None, min_depth=5):
    """Descriptive statistics for a tract the model has already located.

    None of this decides anything any more: the verdict comes from `gconv_model`, which weighs
    the same reads by base quality and marginalises the mismapping rate out instead of testing
    an average against a threshold. What survives here is what a person reads to see WHY, in
    terms they can check against the BAM: how fixed the donor allele is inside the tract, whether
    donor alleles also turn up outside it, and how many single molecules carry both.

    `in_any_tract` is every site belonging to ANY tract found in this locus. Sites of a second,
    genuine tract are not "outside" this one, and counting them there would make each of two real
    conversions look like the other's contamination.
    """
    inside = set(tract)
    excluded = set(in_any_tract) if in_any_tract else inside
    # A site too shallow to genotype is undetermined, and undetermined testifies in NEITHER
    # direction. Excluding it from the tract but still counting it as "outside" lets an unreliable
    # fraction, measured off a handful of reads, drag the outside average over the threshold and
    # flip a real conversion to mismapping.
    determined = {p for p in positions
                  if counts[p]["depth"] >= min_depth and counts[p]["donor_af"] is not None}
    outside = [p for p in positions if p not in excluded and p in determined]

    outside_af = [counts[p]["donor_af"] for p in outside if counts[p]["donor_af"] is not None]
    donor_af_outside = sum(outside_af) / len(outside_af) if outside_af else None

    cis_reads = 0          # one molecule carrying the donor allele at 2+ sites inside the tract
    breakpoint_reads = 0   # one molecule carrying donor INSIDE and acceptor OUTSIDE: the signature
    donor_only_reads = 0   # carries donor everywhere it reaches: what a mismapped read looks like
    for calls in per_read.values():
        in_donor = [p for p, w in calls.items() if p in inside and w == "donor"]
        out_acc = [p for p, w in calls.items()
                   if p not in excluded and p in determined and w == "acceptor"]
        out_donor = [p for p, w in calls.items()
                     if p not in excluded and p in determined and w == "donor"]
        if len(in_donor) >= 2:
            cis_reads += 1
        if in_donor and out_acc:
            breakpoint_reads += 1
        if in_donor and out_donor and not out_acc:
            donor_only_reads += 1

    # Sites inside the tract that were too shallow to genotype. The model does not need them to
    # be: a site with two reads contributes what two reads are worth and no more. But a tract
    # resting on a handful of measured sites and a stretch of empty ones looks more solid in the
    # coordinates than it is, so the count is reported next to it.
    n_undetermined = sum(1 for p in tract
                         if counts[p]["depth"] < min_depth or counts[p]["donor_af"] is None)

    afs = [counts[p]["donor_af"] for p in tract if counts[p]["donor_af"] is not None]
    return {
        "n_sites": len(tract),
        "n_undetermined": n_undetermined,
        "start": min(tract),
        "end": max(tract),
        "span_bp": max(tract) - min(tract) + 1,
        # How many diagnostic sites of this locus sit OUTSIDE the tract. Zero means the tract
        # covers the whole locus, and boundedness - the core claim - cannot be tested at all.
        "n_sites_outside": len(outside),
        "donor_af_in": round(sum(afs) / len(afs), 4) if afs else None,
        "donor_af_outside": round(donor_af_outside, 4) if donor_af_outside is not None else None,
        "min_depth": min(counts[p]["depth"] for p in tract),
        "cis_reads": cis_reads,
        "breakpoint_reads": breakpoint_reads,
        "donor_only_reads": donor_only_reads,
    }


def depleted_runs(positions, counts, donor_counts, min_depth=5, max_ratio=0.4, min_sites=3,
                  exclude=None, max_local=0.6, min_enrichment=1.25):
    """Runs of diagnostic sites where the acceptor has lost its reads to the donor.

    A conversion makes the acceptor identical to the donor over the tract. Once the tract is
    LONGER than the library insert, a read pair falling inside it has no unique anchor left, and
    the aligner assigns it to one copy arbitrarily. In practice the acceptor's tract is depleted
    and the donor's matching region is enriched by the same reads, so the allele evidence this
    tool is built on evaporates exactly when the conversion is most complete.

    That failure is silent: the sites simply go undetermined and no tract is called. This reports
    the depletion instead, as a candidate needing a different kind of evidence (longer reads, or
    read-depth analysis over the pair). It is NOT a conversion call: a deletion of the acceptor
    produces the same picture.

    `exclude` holds the sites of tracts already called. A long conversion often keeps enough
    reads at one end to be called there and loses them everywhere else, so the two findings sit
    side by side in the same locus. An excluded site ends the run it is in rather than
    disqualifying it, which would throw away the depleted part of exactly the tract that was
    hardest to see.

    Everything here is measured against the LOCUS's own coverage, not against an absolute depth.
    It used to fire whenever the acceptor sat under `min_depth` while the donor held most of the
    reads, which is true of every low-coverage paralog whether or not anything moved: on a clean
    negative control it reported a shift over a locus running at 4x throughout, where the acceptor
    was barely dipping and the donor barely gaining. Reads that MOVED leave two marks, and both
    are now required: the acceptor falls well below its own level elsewhere, and the donor rises
    above its own. On the measured case those were 0.46 and 1.53; on the false one, 0.75 and 1.12.
    """
    skip = set(exclude) if exclude else set()
    usable = [p for p in positions if p not in skip]
    if len(usable) < min_sites:
        return []

    def lopsided(p):
        """The donor holds the reads at this site."""
        acc_dp, don_dp = counts[p]["depth"], donor_counts.get(p, 0)
        total = acc_dp + don_dp
        return (acc_dp < min_depth and total >= 2 * min_depth
                and (acc_dp / total if total else 1.0) <= max_ratio)

    # What the acceptor runs at where it is NOT suspected of having lost anything. Taking the
    # median over the whole locus instead would be circular: a depletion covering half the sites
    # drags the very number it is being compared against down with it, and the run then breaks up
    # on its own deepest site.
    healthy = [counts[p]["depth"] for p in usable if not lopsided(p)]
    acc_base = _median(healthy) if len(healthy) >= 2 else None

    runs, run = [], []
    for p in positions:
        if p not in skip and lopsided(p) and (acc_base is None
                                              or counts[p]["depth"] <= max_local * acc_base):
            run.append(p)
        else:
            if len(run) >= min_sites:
                runs.append(run)
            run = []
    if len(run) >= min_sites:
        runs.append(run)

    # And the donor has to have GAINED. Without this, an ordinary coverage hole reads exactly
    # like a tract whose reads walked next door, and a locus simply running at 4x throughout
    # reported a shift on a clean negative control.
    kept = []
    for r in runs:
        inside = set(r)
        elsewhere = [donor_counts.get(p, 0) for p in usable if p not in inside]
        don_base = _median(elsewhere) if elsewhere else 0.0
        don_here = _median([donor_counts.get(p, 0) for p in r])
        if don_base > 0 and don_here >= min_enrichment * don_base:
            kept.append(r)
    return kept


def _median(values):
    v = sorted(values)
    if not v:
        return 0.0
    mid = len(v) // 2
    return float(v[mid]) if len(v) % 2 else (v[mid - 1] + v[mid]) / 2.0


def load_sites(path):
    """Diagnostic sites grouped by acceptor locus, from paralog_map.py."""
    loci = defaultdict(lambda: {"sites": {}, "donor": None, "contig": None})
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < len(header):
                continue
            key = f[idx["pair_id"]]
            loc = loci[key]
            loc["contig"] = f[idx["acceptor"]]
            loc["donor"] = f[idx["donor"]]
            loc["sites"][int(f[idx["acc_pos"]])] = (f[idx["acc_base"]].upper(),
                                                    f[idx["don_base"]].upper())
            # Which markers are deletions, so the model can price them differently: two copies
            # sharing an indel is a far longer coincidence than two copies sharing a base.
            if "kind" in idx:
                loc.setdefault("kind", {})[int(f[idx["acc_pos"]])] = f[idx["kind"]].strip()
            # Optional: only the depletion check needs it, and a hand-made sites file may omit it.
            if "don_pos" in idx and f[idx["don_pos"]].strip().isdigit():
                loc.setdefault("donor_pos", {})[int(f[idx["acc_pos"]])] = int(f[idx["don_pos"]])
            # Which copy the difference is on, when an outgroup was given to say so.
            if "polarity" in idx:
                loc.setdefault("polarity", {})[int(f[idx["acc_pos"]])] = \
                    f[idx["polarity"]].strip()
    return loci


def polarity_counts(tract, polarity):
    """How the sites under a tract are polarised: (derived, ancestral, unreadable).

    `derived` sites are the ones that can support a conversion: the reference's acceptor carries
    the ancestral base there, so a sample carrying the donor's base has changed. `ancestral` sites
    say the opposite, that the REFERENCE carries the derived base and a sample carrying the
    donor's has changed nothing.

    A tract made mostly of the second kind is a real observation about the reference and not
    about the sample, and without an outgroup the two are the same picture.
    """
    if not polarity:
        return 0, 0, 0
    kinds = [polarity.get(p, "") for p in tract]
    derived = sum(1 for k in kinds if k == "derived")
    ancestral = sum(1 for k in kinds if k == "ancestral")
    return derived, ancestral, len(kinds) - derived - ancestral


# One row per paralog pair analysed, whatever came of it. Feeds the cohort pass, which needs to
# see the quiet loci as well as the loud ones.
LOCUS_COLUMNS = ["sample", "pair_id", "contig", "donor", "n_sites", "n_reads", "start", "end",
                 "n_tract_sites", "log10_bf", "log10_bf_vs_null", "tract_af", "mismap_frac",
                 "mut_rate"]

COLUMNS = ["sample", "pair_id", "contig", "donor", "verdict", "reason", "start", "end", "span_bp",
           "don_start", "don_end",
           "post_conv", "log10_bf", "log10_bf_vs_null", "tract_af", "mismap_frac", "mut_rate",
           "start_ci", "end_ci",
           "n_sites", "n_sites_outside", "n_undetermined", "donor_af_in", "donor_af_outside",
           "min_depth", "cis_reads", "breakpoint_reads", "donor_only_reads",
           "n_derived", "n_ancestral", "n_unpolarised",
           "donor_swap_af", "donor_swap_af_outside",
           "locus_cn", "expected_af"] + ga.ANNOTATION_COLUMNS


# Parameters that change the numbers in the output. Recorded in the file itself, because a
# results table whose verdicts depend on seventeen settings and does not say what they were
# cannot be checked, compared against another run, or reproduced a year later.
PROVENANCE_ARGS = ["min_bf", "report_bf", "prior", "mean_tract_bp", "max_tract_bp", "max_tracts",
                   "mut_rate", "indel_factor", "min_tract_af", "min_mismap", "min_sites",
                   "min_depth", "min_bq"]


def provenance_line(args):
    """A `#` header naming the tool and every setting that moved a number in this file."""
    kv = " ".join(f"{k}={getattr(args, k)}" for k in PROVENANCE_ARGS if hasattr(args, k))
    return f"# gene_conversion.py sites={os.path.basename(args.sites)} {kv}"


def write_tsv(path, columns, rows, header=None):
    with open(path, "w") as fh:
        if header:
            fh.write(header + "\n")
        fh.write("\t".join(columns) + "\n")
        for r in rows:
            fh.write("\t".join("" if r.get(c) is None else str(r.get(c, "")) for c in columns) + "\n")


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--sites", required=True, help="diagnostic sites TSV from paralog_map.py")
    p.add_argument("--bam", required=True,
                   help="alignment BAM. Use the PRE-FILTER one: the mappability filter removes "
                        "exactly the reads this analysis depends on")
    p.add_argument("--sample", required=True)
    p.add_argument("-o", "--output", required=True)
    p.add_argument("--output-loci",
                   help="one row per paralog pair analysed, reported or not. This is what the "
                        "cohort pass reads: a locus that says nothing in one sample only means "
                        "something next to the same locus in the others")
    p.add_argument("--min-bf", type=float, default=3.0,
                   help="log10 Bayes factor over the best alternative before a tract is called "
                        "a conversion")
    p.add_argument("--report-bf", type=float, default=1.0,
                   help="log10 Bayes factor below which a tract is not written out at all")
    p.add_argument("--prior", type=float, default=0.01,
                   help="prior probability that a given paralog pair carries a conversion tract")
    p.add_argument("--mean-tract-bp", type=float, default=1000.0,
                   help="mean of the exponential prior on tract length")
    p.add_argument("--max-tract-bp", type=int, default=10000,
                   help="longest tract considered; beyond this the reads have moved to the donor "
                        "and the coverage_shift check is what applies")
    p.add_argument("--max-tracts", type=int, default=2,
                   help="conversion tracts looked for per paralog pair")
    p.add_argument("--mut-rate", type=float, default=gm.MUT_RATE,
                   help="probability that a diagnostic site carries the donor base by independent "
                        "substitution. This is what sets how many sites a tract needs")
    p.add_argument("--indel-factor", type=float, default=gm.INDEL_MUT_FACTOR,
                   help="how much less likely a shared DELETION is to have arisen twice by "
                        "chance than a shared substitution. This is the whole value of an indel "
                        "marker: it is worth more, not merely one more")
    p.add_argument("--min-tract-af", type=float, default=0.25,
                   help="fraction of the reads that must carry a tract before it is called. "
                        "Under this a minority conversion cannot be told from contamination")
    p.add_argument("--min-mismap", type=float, default=0.2,
                   help="fitted donor-read fraction at which a locus is reported as mismapping")
    p.add_argument("--min-sites", type=int, default=3,
                   help="diagnostic sites a paralog pair must have before it is analysed, and the "
                        "length a depleted run must reach")
    p.add_argument("--min-depth", type=int, default=5,
                   help="informative depth below which a site counts as undetermined in the "
                        "descriptive columns and in the depletion check")
    p.add_argument("--min-bq", type=int, default=13, help="base-quality floor")
    p.add_argument("--genome-depth", type=float, default=None,
                   help="the sample's genome-wide mean depth. With it, each locus reports how "
                        "many copies its reads look like, and the read fraction a clonal "
                        "conversion is expected to reach at that copy number")
    p.add_argument("--het-fraction", type=float, default=None,
                   help="fraction of the sample's variant sites that are heterozygous, as a "
                        "genome-wide signal that the sample is mixed. Reported beside a diluted "
                        "tract, never used to decide one")
    p.add_argument("--reciprocal-af", type=float, default=0.5,
                   help="fraction of the DONOR's reads carrying the acceptor's bases at which "
                        "the event is an exchange between the copies rather than a conversion "
                        "of one by the other (default: %(default)s)")
    p.add_argument("--gff", help="GFF3 for the reference. With it, each tract reports the genes "
                                 "it covers and what its copied bases do to their proteins")
    p.add_argument("--reference", help="reference FASTA, needed with --gff to read the codons")
    p.add_argument("--samtools", default="samtools")
    a = p.parse_args(argv)
    # Check the ranges here rather than letting the model take log(0) on the first locus: by then
    # the run has already spent a samtools call on it, and the error names a variable rather than
    # the flag that set it.
    if not 0.0 < a.mut_rate < 1.0:
        p.error(f"--mut-rate must be between 0 and 1, exclusive (got {a.mut_rate})")
    if not 0.0 <= a.prior <= 1.0:
        p.error(f"--prior must be between 0 and 1 (got {a.prior})")
    if a.mean_tract_bp <= 0:
        p.error(f"--mean-tract-bp must be positive (got {a.mean_tract_bp})")
    if a.max_tracts < 1:
        p.error(f"--max-tracts must be at least 1 (got {a.max_tracts})")
    return a


def main(argv=None) -> int:
    a = parse_args(argv)
    loci = load_sites(a.sites)
    # Read once for the whole sample. Absent, every annotation column comes out empty rather than
    # guessed at, which is the default: a reference has a GFF or it does not.
    features = ga.parse_cds(a.gff) if a.gff else []
    seqs = ga.read_fasta(a.reference) if a.gff and a.reference else {}

    rows, locus_rows = [], []
    for pair_id, loc in sorted(loci.items(), key=lambda kv: int(kv[0])):
        positions = sorted(loc["sites"])
        if len(positions) < a.min_sites:
            continue
        region = f"{loc['contig']}:{positions[0]}-{positions[-1]}"
        proc = subprocess.run([a.samtools, "view", a.bam, region],
                              capture_output=True, text=True)
        if proc.returncode != 0:
            raise RuntimeError(f"samtools view failed on {region}: {proc.stderr.strip()[:300]}")

        don_pos_of = loc.get("donor_pos", {})
        obs = collect_read_observations(proc.stdout.splitlines(), loc["sites"])
        per_read = allele_view(obs, loc["sites"], a.min_bq)
        counts = pileup(per_read, positions)

        # The DONOR side of the same pair. Read before the tracts are judged, because what the
        # donor carries decides what kind of event this is: a conversion leaves it alone, and a
        # reciprocal exchange does not.
        donor_counts, donor_swapped = donor_side(a, loc, positions)
        cn_ratio, copies, expected_af = locus_copies(counts, positions, a.genome_depth,
                                                     a.min_depth)

        # The model decides where the tracts are and whether they mean anything. It sees the
        # bases and their qualities, not the pileup summary, and it compares a tract against the
        # two other ways donor bases turn up here: independent substitution, and reads that came
        # from the donor in the first place.
        _, delta = gm.delta_matrix(obs, positions, loc["sites"], a.min_bq)
        kinds = loc.get("kind", {})
        is_del = [kinds.get(p) == "del" for p in positions]
        fits = gm.segment(delta, positions, max_tracts=a.max_tracts, min_report_bf=a.report_bf,
                          prior=a.prior, mean_span_bp=a.mean_tract_bp, mut_rate=a.mut_rate,
                          max_span_bp=a.max_tract_bp, is_del=is_del,
                          indel_factor=a.indel_factor)
        # What the best tract at this locus looked like, whatever it came to. The tracts file
        # holds findings and is what a person reads; this holds the evidence at EVERY locus,
        # including the loci where there was none, and it is what the cohort pass needs. A locus
        # silent in this sample is only informative next to the same locus in the others, and
        # "no row" cannot tell "nothing there" from "not written out".
        if fits:
            f0 = fits[0]
            locus_rows.append({
                "sample": a.sample, "pair_id": pair_id, "contig": loc["contig"],
                "donor": loc["donor"], "n_sites": len(positions), "n_reads": f0["n_reads"],
                "start": f0["start"], "end": f0["end"], "n_tract_sites": f0["map_j"] - f0["map_i"] + 1,
                "log10_bf": round(f0["log10_bf"], 2),
                "log10_bf_vs_null": round(f0["log10_bf_null"], 2),
                "tract_af": round(f0["tract_af"], 3),
                "mismap_frac": round(f0["mismap_frac"], 4),
                "mut_rate": round(f0["mut_rate"], 5)})

        keep = [f for f in fits if worth_reporting(f, a.report_bf, a.min_mismap)]
        tracts = [positions[f["map_i"]:f["map_j"] + 1] for f in keep]
        in_any = {p for t in tracts for p in t}
        for fit, tract in zip(keep, tracts):
            ev = tract_evidence(tract, positions, counts, per_read, in_any, a.min_depth)
            # Where in the DONOR the copied stretch came from. A gene family reports one event
            # through several relationships, and this is what says whether they are pointing at
            # the same piece of donor or at different ones, which is the difference between
            # redundancy and a genuine choice of source.
            dpos = [don_pos_of[p] for p in tract if p in don_pos_of]
            covers_locus = fit["map_i"] == 0 and fit["map_j"] == len(positions) - 1
            verdict, reason = gm.verdict(fit, covers_locus, a.min_bf, a.min_mismap,
                                         a.min_tract_af)
            reason += dilution_notes(fit.get("tract_af"), a.min_tract_af, cn_ratio, copies,
                                     expected_af, a.het_fraction)
            # Which copy the change is on, where an outgroup was given to say so. The model has
            # no view on this: it sees the acceptor carrying the donor's base and that is the
            # same picture whether the sample changed or the reference did.
            derived, ancestral, unread = polarity_counts(tract, loc.get("polarity"))
            # Non-reciprocal is not a detail of the definition, it is the definition. If the donor
            # has taken on the acceptor's bases over the same stretch, the two copies swapped, and
            # a swap is an unequal crossover: one event, both copies changed, and the consequences
            # for the family are not the consequences of a conversion.
            swapped, swapped_out = reciprocity(tract, positions, donor_swapped)
            if is_exchange(swapped, swapped_out, a.reciprocal_af):
                verdict = "reciprocal_exchange"
                reason = (
                    f"the donor carries the acceptor's bases over the same stretch "
                    f"({swapped:.0%} of its reads there against {swapped_out:.0%} outside it), so "
                    "both copies changed and only over this stretch. That is an exchange between "
                    "them rather than one copy being overwritten, and gene conversion is "
                    "non-reciprocal by definition")
            elif ancestral > derived:
                verdict = "reference_derived"
                reason = (
                    f"{ancestral} of the {derived + ancestral} polarised site(s) here carry the "
                    "ancestral base in the reads and a derived one in the reference, so it is the "
                    "REFERENCE's copy that was converted or mutated and this sample retains what "
                    "the outgroup has. Not a conversion in this sample")
            if fit["stride"] > 1:
                reason += (f"; breakpoints resolved to every {fit['stride']} diagnostic sites "
                           "because the locus has too many to search exhaustively")
            rows.append({"sample": a.sample, "pair_id": pair_id, "contig": loc["contig"],
                         "donor": loc["donor"], "verdict": verdict, "reason": reason,
                         "post_conv": round(fit["post_conv"], 4),
                         "log10_bf": round(fit["log10_bf"], 2),
                         "log10_bf_vs_null": round(fit["log10_bf_null"], 2),
                         "tract_af": round(fit["tract_af"], 3),
                         "mismap_frac": round(fit["mismap_frac"], 4),
                         "mut_rate": round(fit["mut_rate"], 5),
                         "start_ci": "{}-{}".format(*fit["start_ci"]),
                         "end_ci": "{}-{}".format(*fit["end_ci"]),
                         "don_start": min(dpos) if dpos else None,
                         "don_end": max(dpos) if dpos else None,
                         "n_derived": derived, "n_ancestral": ancestral,
                         "n_unpolarised": unread,
                         "donor_swap_af": None if swapped is None else round(swapped, 3),
                         "donor_swap_af_outside": (None if swapped_out is None
                                                   else round(swapped_out, 3)),
                         "locus_cn": cn_ratio, "expected_af": expected_af,
                         # Only the diagnostic sites change: everywhere else the two copies are
                         # identical, so a conversion there is invisible and inconsequential.
                         # Deletion markers are left out, being a frameshift question rather than
                         # a codon one, and half an answer there is worse than none.
                         **ga.annotate_tract(
                             features, seqs, loc["contig"], tract[0], tract[-1],
                             {p: loc["sites"][p][1] for p in tract
                              if loc["sites"][p][1] != GAP}),
                         **ev})

        for run in depleted_runs(positions, counts, donor_counts, a.min_depth,
                                 min_sites=a.min_sites, exclude=in_any):
            acc = sum(counts[p]["depth"] for p in run)
            don = sum(donor_counts.get(p, 0) for p in run)
            rows.append({
                "sample": a.sample, "pair_id": pair_id, "contig": loc["contig"],
                "donor": loc["donor"], "verdict": "coverage_shift",
                "reason": (f"acceptor depleted ({acc} reads) while the donor carries {don} over the "
                           "same sites; a tract longer than the insert loses its reads to the donor, "
                           "but so does a deletion"),
                "n_sites": len(run), "start": min(run), "end": max(run),
                "span_bp": max(run) - min(run) + 1, "n_sites_outside": len(positions) - len(run),
                "n_undetermined": 0, "donor_af_in": None, "donor_af_outside": None,
                "n_derived": 0, "n_ancestral": 0, "n_unpolarised": len(run),
                "donor_swap_af": None, "donor_swap_af_outside": None, "locus_cn": cn_ratio, "expected_af": expected_af,
                **ga.annotate_tract(features, seqs, loc["contig"], run[0], run[-1], {}),
                "min_depth": min(counts[p]["depth"] for p in run),
                "cis_reads": 0, "breakpoint_reads": 0, "donor_only_reads": 0})

    prov = provenance_line(a)
    write_tsv(a.output, COLUMNS, rows, prov)
    if a.output_loci:
        write_tsv(a.output_loci, LOCUS_COLUMNS, locus_rows, prov)

    called = sum(1 for r in rows if r["verdict"] == "gene_conversion")
    sys.stderr.write(f"[gene_conversion] {a.sample}: {len(rows)} candidate tract(s), "
                     f"{called} called as conversion, {len(locus_rows)} locus/loci measured\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
