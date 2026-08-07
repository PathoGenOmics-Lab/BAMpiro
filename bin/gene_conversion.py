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

Reads from the donor that mismap onto the acceptor produce the same per-site picture: donor bases
where the reference expects acceptor bases. Three things separate the two, and this tool reports
all three rather than collapsing them into one number:

* **Where the donor alleles are.** A conversion is bounded. Mismapping is not: mismapped reads
  carry donor alleles at every diagnostic site of the locus, including outside the candidate
  tract. `donor_af_outside` is the direct test.
* **How fixed they are.** In a clonal sample a converted locus is homozygous for the donor allele
  (allele fraction near 1). Mismapping mixes donor and acceptor reads at whatever ratio the
  aligner produced, so the fraction sits in between and, tellingly, is much the same at every site.
* **What single reads carry.** This is the one a mismapping cannot fake. A read that spans a
  breakpoint carries donor alleles on one side and acceptor alleles on the other, in cis, on the
  same physical molecule. A mismapped read is a donor read: it carries donor alleles at every site
  it covers, and never crosses back. `breakpoint_reads` counts the former.

MAPQ is deliberately NOT filtered on. Paralogous regions are exactly where an aligner assigns
MAPQ 0, so the usual quality gate would discard the entire signal this tool exists to find.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from collections import defaultdict

# CIGAR operations that consume the reference, the read, or both.
CONSUMES_REF = set("MDN=X")
CONSUMES_READ = set("MIS=X")


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


def read_bases_at(pos, cigar, seq, qual, wanted, min_bq=13):
    """Bases this read carries at the reference positions in `wanted`.

    Walks the CIGAR so an insertion or a soft clip does not shift every base after it, which is
    the classic way to read the right position off the wrong read.

    Returns {ref_pos (1-based): base}. A position inside a deletion or below the base-quality
    floor is omitted rather than guessed at.
    """
    out = {}
    ref, qry = pos, 0
    for length, op in parse_cigar(cigar):
        if op in CONSUMES_REF and op in CONSUMES_READ:          # M, =, X
            for i in range(length):
                p = ref + i
                if p in wanted:
                    q = qry + i
                    if q < len(seq) and (not qual or qual == "*"
                                         or (q < len(qual) and (ord(qual[q]) - 33) >= min_bq)):
                        out[p] = seq[q].upper()
            ref += length
            qry += length
        elif op in CONSUMES_REF:                                 # D, N: no base to read
            ref += length
        elif op in CONSUMES_READ:                                # I, S: consumes no reference
            qry += length
    return out


def collect_read_alleles(sam_lines, sites, min_bq=13):
    """Per-read allele calls at the diagnostic sites.

    `sites` maps a reference position to (acceptor_base, donor_base). Returns
    {read_name: {pos: 'acceptor' | 'donor' | 'other'}}, keyed by read NAME so the two mates of a
    pair contribute to one observation of one molecule.
    """
    wanted = set(sites)
    per_read = defaultdict(dict)
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
        for p, base in read_bases_at(pos, cigar, seq, qual, wanted, min_bq).items():
            acc, don = sites[p]
            per_read[name][p] = "acceptor" if base == acc else "donor" if base == don else "other"
    return per_read


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


def call_tracts(positions, counts, min_af=0.7, min_sites=3, min_depth=5):
    """Maximal runs of consecutive diagnostic sites that carry the donor allele.

    Consecutive means adjacent in the diagnostic-site list, not in base coordinates: the sites
    between them are identical in both copies and cannot testify either way.

    A site with too little depth to genotype is UNDETERMINED, and undetermined is not the same as
    "the acceptor allele is here". It bridges rather than breaks: a single coverage dip inside an
    otherwise clean tract would otherwise shatter it into fragments that each fall below
    min_sites, and the whole tract disappears. Coverage dips are ordinary, so that failure mode is
    ordinary too. The bridged sites do not count towards min_sites and do not become tract
    members; tract_evidence reports how many were skipped.
    """
    tracts, run = [], []
    for p in positions:
        c = counts[p]
        if c["depth"] < min_depth or c["donor_af"] is None:
            continue                                    # undetermined: no evidence either way
        if c["donor_af"] >= min_af:
            run.append(p)
        else:
            if len(run) >= min_sites:
                tracts.append(run)
            run = []
    if len(run) >= min_sites:
        tracts.append(run)
    return tracts


def tract_evidence(tract, positions, counts, per_read, in_any_tract=None, min_depth=5):
    """The evidence that separates a real conversion from reads arriving from the donor.

    `in_any_tract` is every site belonging to ANY tract called in this locus. Sites of a second,
    genuine tract are not "outside" this one: counting them there inflates donor_af_outside, and
    since the outside test runs first, two real tracts in one locus would each report the other as
    mismapping. Two conversions in one locus is an ordinary outcome, not a pathological input.
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

    # Diagnostic sites inside the tract's span that were too shallow to genotype. They were
    # bridged rather than allowed to break the tract, so say how many, or the tract looks more
    # solid than the data supports.
    n_undetermined = sum(1 for p in positions
                         if min(tract) < p < max(tract) and p not in inside)

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


def classify(ev, max_outside_af=0.25, min_af_in=0.85):
    """Turn the evidence into a verdict, and say which piece of it decided.

    Deliberately conservative: a tract with no read-level support is reported, not hidden, but it
    is not called a conversion, and the reason says which test it failed.
    """
    if ev["donor_af_outside"] is not None and ev["donor_af_outside"] > max_outside_af:
        return "mismapping", "donor alleles are present outside the tract as well"

    # A read crossing a breakpoint carries donor and acceptor alleles on one molecule. Reads
    # arriving from the donor cannot produce that, so this outranks everything else.
    if ev["breakpoint_reads"] > 0:
        return "gene_conversion", f"{ev['breakpoint_reads']} read(s) cross a breakpoint in cis"

    # No site outside the tract means there is nothing the donor alleles are bounded BY. A whole
    # locus replaced by its paralog and a locus whose reads all arrived from its paralog look
    # identical from site data alone, so absence of contrary evidence is not evidence.
    if ev["n_sites_outside"] == 0:
        return "ambiguous", ("tract covers every diagnostic site, so it is not bounded: "
                             "indistinguishable from wholesale mismapping without a breakpoint read")

    if ev["donor_af_in"] is not None and ev["donor_af_in"] >= min_af_in and ev["cis_reads"] > 0:
        return "gene_conversion", "donor allele near-fixed within a bounded tract"
    return "ambiguous", "no read spans a breakpoint; tract is bounded but unconfirmed"


def depleted_runs(positions, counts, donor_counts, min_depth=5, max_ratio=0.4, min_sites=3):
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
    """
    runs, run = [], []
    for p in positions:
        acc_dp = counts[p]["depth"]
        don_dp = donor_counts.get(p, 0)
        total = acc_dp + don_dp
        depleted = (acc_dp < min_depth and total >= 2 * min_depth
                    and (acc_dp / total if total else 1.0) <= max_ratio)
        if depleted:
            run.append(p)
        else:
            if len(run) >= min_sites:
                runs.append(run)
            run = []
    if len(run) >= min_sites:
        runs.append(run)
    return runs


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
            # Optional: only the depletion check needs it, and a hand-made sites file may omit it.
            if "don_pos" in idx and f[idx["don_pos"]].strip().isdigit():
                loc.setdefault("donor_pos", {})[int(f[idx["acc_pos"]])] = int(f[idx["don_pos"]])
    return loci


COLUMNS = ["sample", "pair_id", "contig", "donor", "verdict", "reason", "start", "end", "span_bp",
           "n_sites", "n_sites_outside", "n_undetermined", "donor_af_in", "donor_af_outside", "min_depth", "cis_reads",
           "breakpoint_reads", "donor_only_reads"]


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--sites", required=True, help="diagnostic sites TSV from paralog_map.py")
    p.add_argument("--bam", required=True,
                   help="alignment BAM. Use the PRE-FILTER one: the mappability filter removes "
                        "exactly the reads this analysis depends on")
    p.add_argument("--sample", required=True)
    p.add_argument("-o", "--output", required=True)
    p.add_argument("--min-af", type=float, default=0.7, help="donor fraction for a site to join a tract")
    p.add_argument("--min-sites", type=int, default=3, help="diagnostic sites needed to call a tract")
    p.add_argument("--min-depth", type=int, default=5, help="informative depth needed at a site")
    p.add_argument("--min-bq", type=int, default=13, help="base-quality floor")
    p.add_argument("--samtools", default="samtools")
    return p.parse_args(argv)


def main(argv=None) -> int:
    a = parse_args(argv)
    loci = load_sites(a.sites)

    rows = []
    for pair_id, loc in sorted(loci.items(), key=lambda kv: int(kv[0])):
        positions = sorted(loc["sites"])
        if len(positions) < a.min_sites:
            continue
        region = f"{loc['contig']}:{positions[0]}-{positions[-1]}"
        proc = subprocess.run([a.samtools, "view", a.bam, region],
                              capture_output=True, text=True)
        if proc.returncode != 0:
            raise RuntimeError(f"samtools view failed on {region}: {proc.stderr.strip()[:300]}")

        per_read = collect_read_alleles(proc.stdout.splitlines(), loc["sites"], a.min_bq)
        counts = pileup(per_read, positions)
        tracts = call_tracts(positions, counts, a.min_af, a.min_sites, a.min_depth)
        in_any = {p for t in tracts for p in t}
        for tract in tracts:
            ev = tract_evidence(tract, positions, counts, per_read, in_any, a.min_depth)
            verdict, reason = classify(ev)
            rows.append({"sample": a.sample, "pair_id": pair_id, "contig": loc["contig"],
                         "donor": loc["donor"], "verdict": verdict, "reason": reason, **ev})

        # Depth on the DONOR side of the same pair, to catch the tracts whose reads moved there.
        don_pos = loc.get("donor_pos", {})
        donor_counts = {}
        if don_pos:
            dvals = sorted(don_pos.values())
            dproc = subprocess.run(
                [a.samtools, "view", a.bam, f"{loc['donor']}:{dvals[0]}-{dvals[-1]}"],
                capture_output=True, text=True)
            if dproc.returncode == 0:
                dsites = {dp: ("N", "N") for dp in dvals}     # only depth matters here
                dreads = collect_read_alleles(dproc.stdout.splitlines(), dsites, a.min_bq)
                dcounts = pileup(dreads, dvals)
                donor_counts = {ap: dcounts[dp]["depth"] for ap, dp in don_pos.items()}

        for run in depleted_runs(positions, counts, donor_counts, a.min_depth,
                                 min_sites=a.min_sites):
            if set(run) & in_any:
                continue                                     # already reported as a tract
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
                "min_depth": min(counts[p]["depth"] for p in run),
                "cis_reads": 0, "breakpoint_reads": 0, "donor_only_reads": 0})

    with open(a.output, "w") as fh:
        fh.write("\t".join(COLUMNS) + "\n")
        for r in rows:
            fh.write("\t".join("" if r.get(c) is None else str(r.get(c, "")) for c in COLUMNS) + "\n")

    called = sum(1 for r in rows if r["verdict"] == "gene_conversion")
    sys.stderr.write(f"[gene_conversion] {a.sample}: {len(rows)} candidate tract(s), "
                     f"{called} called as conversion\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
