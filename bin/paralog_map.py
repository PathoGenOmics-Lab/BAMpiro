#!/usr/bin/env python3
"""Recover the paralog map the repeat-masking step already computes and throws away.

PREPARE_REFERENCE aligns the reference against itself with `nucmer --maxmatch --nosimplify`
and then collapses the result to a flat list of repetitive intervals. That collapse discards
the one thing gene-conversion detection needs: WHICH copy aligns to WHICH, at what offset and
in which orientation, and at which positions the two copies actually differ.

This reads the same `.delta` and emits both halves of that map:

    paralog_map.py --delta self_aln.delta --out-pairs pairs.tsv --out-sites sites.tsv

* pairs.tsv - one row per paralogous interval pair (the donor/acceptor graph of the reference)
* sites.tsv - one row per PARALOG-DIAGNOSTIC position: a base where the two copies differ.

Only the diagnostic sites can carry evidence of conversion. Everywhere else the two copies are
identical, so a read is uninformative there by construction, and the tract caller works over
this list rather than over the whole locus.

Both files come from MUMmer's own `show-coords` and `show-snps`, so no alignment is redone.
"""

from __future__ import annotations

import argparse
import bisect
import os
import pathlib
import subprocess
import sys

# show-coords -r -c -l -T: S1 E1 S2 E2 LEN1 LEN2 %IDY LENR LENQ COVR COVQ TAG1 TAG2
COORDS_FIELDS = 13
# show-snps -l -r -T -H: P1 SUB1 SUB2 P2 BUFF DIST R Q LENR LENQ FRM1 FRM2 TAG1 TAG2
SNPS_FIELDS = 14

# Coordinate bucket for the site-to-pair lookup. Paralogous alignments run to a few kb, so at
# this size a pair lands in one or two buckets and a site only ever looks at its neighbours.
BIN_SIZE = 100_000


def rehomed_delta(delta, ref_fasta, query_fasta=None, out=None):
    """A copy of `delta` whose header names files that exist here.

    nucmer writes the ABSOLUTE paths of its two inputs into the first line of the delta, and
    `show-snps` opens them to report the bases around each difference. Those paths are the ones
    the alignment ran under. Under a scheduler that gives each task its own scratch directory,
    that path belonged to a different task on possibly a different node and is gone by the time
    this runs, so show-snps fails with "Could not open file" on a delta that is otherwise intact.

    Rewriting one line is the whole fix. Everything after it is offsets, which do not care where
    the sequence is read from.
    """
    lines = pathlib.Path(delta).read_text().splitlines()
    if not lines:
        return delta
    q = query_fasta or ref_fasta
    lines[0] = f"{os.path.abspath(ref_fasta)} {os.path.abspath(q)}"
    out = out or (delta + ".rehomed")
    pathlib.Path(out).write_text("\n".join(lines) + "\n")
    return out


def _run(cmd):
    """Run a MUMmer command and return its stdout, failing loudly rather than silently empty."""
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"{cmd[0]} failed ({proc.returncode}): {proc.stderr.strip()[:400]}")
    return proc.stdout


def parse_coords(text, min_identity=90.0, min_length=200):
    """Paralogous interval pairs from `show-coords -r -c -l -T` output.

    The self-diagonal (every contig aligned to itself over its whole length) is dropped: it is an
    artefact of aligning a sequence against itself, not a paralog. A pair is reported once per
    direction, because acceptor and donor are not interchangeable downstream.
    """
    pairs = []
    for line in text.splitlines():
        f = line.rstrip("\n").split("\t")
        if len(f) < COORDS_FIELDS:
            continue
        if not (f[0].isdigit() and f[1].isdigit() and f[2].isdigit() and f[3].isdigit()):
            continue                                    # header and banner lines
        s1, e1, s2, e2 = (int(f[0]), int(f[1]), int(f[2]), int(f[3]))
        idy = float(f[6])
        tag1, tag2 = f[11], f[12]

        # nucmer reports the query interval reversed when the match is on the other strand.
        strand = "-" if s2 > e2 else "+"
        if s1 > e1:
            s1, e1 = e1, s1
        if s2 > e2:
            s2, e2 = e2, s2

        if tag1 == tag2 and s1 == s2 and e1 == e2:
            continue                                    # self-diagonal
        if idy < min_identity or (e1 - s1 + 1) < min_length:
            continue
        pairs.append({
            "acceptor": tag1, "acc_start": s1, "acc_end": e1,
            "donor": tag2, "don_start": s2, "don_end": e2,
            "strand": strand, "identity": idy, "length": e1 - s1 + 1,
        })
    return pairs


def parse_snps(text):
    """Paralog-diagnostic sites from `show-snps -l -r -T -H` output.

    Two kinds of difference come out of this, and both are diagnostic.

    A **substitution** is one position where the copies carry different bases, and it is what the
    tract caller was built on.

    An **indel** is reported with '.' on one side, one row per base. Where the DONOR is the one
    missing bases, the acceptor's own copy still has a coordinate, so a read either has a base
    there or spans it with a deletion, and that is as readable as any substitution. Consecutive
    rows are collapsed into ONE site: a five-base deletion is one event that happened once, not
    five independent observations, and counting it five times is the same overcounting the caller
    avoids by treating molecules rather than sites as its unit.

    Rows where the ACCEPTOR is the one missing bases are skipped. The donor's extra sequence has
    no acceptor coordinate to hang on, and every pair is emitted in both directions, so the same
    difference is available as a donor gap from the other side. They are counted, not used.
    """
    sites, skipped = [], 0
    gaps = []
    for line in text.splitlines():
        f = line.rstrip("\n").split("\t")
        if len(f) < SNPS_FIELDS or not f[0].isdigit():
            continue
        acc_base, don_base = f[1].upper(), f[2].upper()
        strand = "-" if f[11] == "-1" else "+"
        if acc_base == "." or acc_base == don_base:
            skipped += 1
            continue
        if don_base == ".":
            gaps.append((f[12], f[13], strand, int(f[0]), int(f[3]), acc_base))
            continue
        sites.append({
            "acceptor": f[12], "acc_pos": int(f[0]), "acc_base": acc_base,
            "donor": f[13], "don_pos": int(f[3]), "don_base": don_base,
            "strand": strand, "kind": "snp", "length": 1,
        })

    # One site per run of consecutive acceptor positions the donor lacks. show-snps holds the
    # donor position still across a run, which is what makes the run recognisable.
    gaps.sort()
    i = 0
    while i < len(gaps):
        tag1, tag2, strand, acc_pos, don_pos, base = gaps[i]
        j, bases = i, [base]
        while (j + 1 < len(gaps) and gaps[j + 1][:3] == (tag1, tag2, strand)
               and gaps[j + 1][3] == gaps[j][3] + 1 and gaps[j + 1][4] == don_pos):
            j += 1
            bases.append(gaps[j][5])
        sites.append({
            "acceptor": tag1, "acc_pos": acc_pos, "acc_base": "".join(bases),
            "donor": tag2, "don_pos": don_pos, "don_base": "-",
            "strand": strand, "kind": "del", "length": len(bases),
        })
        i = j + 1
    sites.sort(key=lambda s: (s["acceptor"], s["acc_pos"]))
    return sites, skipped


def aligned_intervals(text):
    """Reference intervals an outgroup alignment covers, as {contig: [(start, end)]}.

    Needed because "no difference was reported here" and "nothing was aligned here" are the same
    silence in a SNP file and mean opposite things. A position the outgroup does not reach is not
    a position where the outgroup agrees.
    """
    out = {}
    for line in text.splitlines():
        f = line.rstrip("\n").split("\t")
        if len(f) < COORDS_FIELDS or not (f[0].isdigit() and f[1].isdigit()):
            continue
        s, e = int(f[0]), int(f[1])
        out.setdefault(f[11], []).append((min(s, e), max(s, e)))

    # Merged, and that is not tidying. An outgroup is aligned with the same `--maxmatch` that
    # produces the paralog map, so it emits nested and overlapping alignments on purpose. The
    # lookup below finds the last interval starting at or before a position and asks only that
    # one, so a position inside a long alignment that also contains a short nested one was
    # answered by the short one and came back uncovered. Whole stretches of the outgroup then
    # went unpolarised, and an unpolarised site is indistinguishable from one the outgroup
    # genuinely does not reach.
    for contig, spans in out.items():
        spans.sort()
        merged = []
        for start, end in spans:
            if merged and start <= merged[-1][1] + 1:
                merged[-1] = (merged[-1][0], max(merged[-1][1], end))
            else:
                merged.append((start, end))
        out[contig] = merged
    return out


def ancestral_bases(snps_text, coords_text):
    """What an outgroup carries at each reference position: (differences, covered intervals).

    `differences` holds only the positions where the outgroup differs from the reference, which is
    what `show-snps` reports. Everywhere inside a covered interval the two agree, so the reference
    base IS the ancestral base there and no entry is needed.

    Outgroup indels are recorded as None rather than skipped. A position the outgroup has deleted
    has no ancestral base to speak of, and treating it as agreement would hand the reference's own
    allele the authority of the ancestor.
    """
    diffs = {}
    for line in snps_text.splitlines():
        f = line.rstrip("\n").split("\t")
        if len(f) < SNPS_FIELDS or not f[0].isdigit():
            continue
        ref_base, anc_base = f[1].upper(), f[2].upper()
        if ref_base == ".":
            continue                      # the outgroup has extra sequence; no reference position
        diffs[(f[12], int(f[0]))] = None if anc_base == "." else anc_base
    return diffs, aligned_intervals(coords_text)


def _covered(intervals, contig, pos):
    spans = intervals.get(contig)
    if not spans:
        return False
    i = bisect.bisect_right(spans, (pos, float("inf")))
    return i > 0 and spans[i - 1][1] >= pos


def ancestral_at(diffs, intervals, contig, pos, ref_base):
    """The outgroup's base at one reference position, or None when it cannot be read there."""
    if (contig, pos) in diffs:
        return diffs[(contig, pos)]
    return ref_base if _covered(intervals, contig, pos) else None


def polarise(sites, diffs, intervals):
    """Say, per diagnostic site, which copy the change is on.

    A sample whose acceptor carries the DONOR's base at a diagnostic site has either changed or
    not, and the site alone cannot tell which. The reference is one genome among many: where its
    acceptor copy carries a derived allele, a sample carrying the donor's base is holding the
    ANCESTRAL state and has changed nothing at all. Read as a conversion, that is the reference's
    history reported as the sample's.

    An outgroup settles it, per site:

    * `derived` - the reference's acceptor is ancestral, so a sample carrying the donor's base
      there has changed. This is the site that can support a conversion.
    * `ancestral` - the reference's acceptor is derived and the ancestral allele is the donor's.
      A sample carrying the donor's base has retained it, and the finding belongs to the
      reference rather than to the sample.
    * `third` - the outgroup carries neither base. Something happened here, but not this.
    * `` - the outgroup does not reach the position, or has deleted it.

    `don_derived` is separate and stronger where it holds: when the DONOR's allele is itself an
    innovation, a sample cannot be carrying it by retention. It had to be copied.
    """
    for s in sites:
        # A deletion site is not polarised, so it is not looked up either. Its `acc_base` is the
        # whole run of bases the donor lacks, and handing that to a per-position lookup gets it
        # back unchanged: the column would then claim the outgroup carries all of them, on the
        # strength of having checked one position.
        if s["kind"] != "snp":
            s["anc_acc"] = s["anc_don"] = ""
            s["polarity"] = ""
            s["don_derived"] = 0
            continue
        anc_acc = ancestral_at(diffs, intervals, s["acceptor"], s["acc_pos"], s["acc_base"])
        anc_don = ancestral_at(diffs, intervals, s["donor"], s["don_pos"], s["don_base"])
        s["anc_acc"] = anc_acc or ""
        s["anc_don"] = anc_don or ""
        if anc_acc is None:
            s["polarity"] = ""
        elif anc_acc == s["acc_base"]:
            s["polarity"] = "derived"
        elif anc_acc == s["don_base"]:
            s["polarity"] = "ancestral"
        else:
            s["polarity"] = "third"
        s["don_derived"] = int(anc_don is not None and anc_don != s["don_base"])
    return sites


def alignment_offset(start, end, acc_start, strand):
    """Where the donor sits relative to the acceptor along an alignment.

    Constant along a gap-free alignment and drifting only with indels, which makes it the thing
    that tells two nested alignments between the same two intervals apart when nothing else can.
    """
    return start - acc_start if strand == "+" else end + acc_start


def sites_within_pairs(sites, pairs):
    """Attach each diagnostic site to every paralog pair it belongs to.

    show-snps reports every difference in the delta, including ones from alignments too short or
    too divergent to be a credible paralog. Tying each site to a pair is what lets the caller
    group sites by locus, and it is less simple than it looks in a SELF-alignment.

    A site belongs to EVERY pair that spans it, not just the first one found. `--maxmatch
    --nosimplify` emits overlapping and nested alignments on purpose, so a tandem repeat produces
    several pairs covering the same bases and one difference is diagnostic for all of them.
    Stopping at the first match cost 140 sites on H37Rv and, worse, left 12 pairs with NO sites at
    all, among them a 5.8 kb paralog at 98% identity that the tract caller then skipped entirely
    for having nothing to work with.

    Taking every span brings back an ambiguity, though: two nested alignments can put two
    different donor positions against the same acceptor position inside one pair's coordinate box,
    and only one of them is that pair's own alignment. The offset picks it out. On H37Rv that is
    196 positions, every one of them resolved.
    """
    # Bucketed by acceptor contig and coordinate, because without the early exit every site
    # otherwise scans every pair. On a finished genome that is 411 pairs and does not matter; on
    # a draft assembly, which is where the repeats are, 3000 pairs against 300k show-snps rows
    # took a minute of quadratic scanning for an answer that touches a handful of them.
    index = {}
    for i, p in enumerate(pairs):
        for b in range(p["acc_start"] // BIN_SIZE, p["acc_end"] // BIN_SIZE + 1):
            index.setdefault((p["acceptor"], b), []).append(i)

    chosen = {}
    for site in sites:
        for i in index.get((site["acceptor"], site["acc_pos"] // BIN_SIZE), ()):
            p = pairs[i]
            if site["donor"] != p["donor"]:
                continue
            if not p["acc_start"] <= site["acc_pos"] <= p["acc_end"]:
                continue
            if not p["don_start"] <= site["don_pos"] <= p["don_end"]:
                continue
            if site["strand"] != p["strand"]:
                continue
            drift = abs(alignment_offset(site["don_pos"], site["don_pos"], site["acc_pos"],
                                         site["strand"])
                        - alignment_offset(p["don_start"], p["don_end"], p["acc_start"],
                                           p["strand"]))
            key = (i, site["acc_pos"])
            if key not in chosen or drift < chosen[key][1]:
                chosen[key] = (site, drift)
    return [{**site, "pair_id": i, "identity": pairs[i]["identity"]}
            for (i, _), (site, _) in sorted(chosen.items())]


PAIR_COLS = ["pair_id", "acceptor", "acc_start", "acc_end", "donor", "don_start", "don_end",
             "strand", "identity", "length", "n_diagnostic"]
SITE_COLS = ["pair_id", "acceptor", "acc_pos", "acc_base", "donor", "don_pos", "don_base",
             "strand", "kind", "length", "anc_acc", "anc_don", "polarity", "don_derived"]


def write_tsv(path, columns, rows):
    with open(path, "w") as fh:
        fh.write("\t".join(columns) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")) for c in columns) + "\n")


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--delta", required=True, help="nucmer .delta from the reference self-alignment")
    p.add_argument("--out-pairs", required=True, help="paralogous interval pairs TSV")
    p.add_argument("--out-sites", required=True, help="paralog-diagnostic sites TSV")
    p.add_argument("--min-identity", type=float, default=90.0,
                   help="ignore pairs below this %% identity (default: %(default)s)")
    p.add_argument("--min-length", type=int, default=200,
                   help="ignore pairs shorter than this many bp (default: %(default)s)")
    p.add_argument("--ancestor-delta",
                   help="nucmer .delta of an OUTGROUP aligned against this reference, as "
                        "`nucmer reference.fasta outgroup.fasta`. With it, each diagnostic site "
                        "says which copy the difference is on, so a sample carrying the donor's "
                        "base can be told from a reference that carries a derived one")
    p.add_argument("--reference",
                   help="the FASTA the delta was built from. nucmer records the absolute path of "
                        "its inputs in the delta header and show-snps opens them, so under a "
                        "scheduler that scratches each task that path is gone and the delta is "
                        "unreadable. With this, the header is rewritten to a file that is here")
    p.add_argument("--outgroup",
                   help="the outgroup FASTA, for the same reason, when --ancestor-delta is given")
    p.add_argument("--show-coords", default="show-coords")
    p.add_argument("--show-snps", default="show-snps")
    return p.parse_args(argv)


def main(argv=None) -> int:
    a = parse_args(argv)

    delta = rehomed_delta(a.delta, a.reference) if a.reference else a.delta
    pairs = parse_coords(_run([a.show_coords, "-r", "-c", "-l", "-T", delta]),
                         min_identity=a.min_identity, min_length=a.min_length)
    # -C is deliberately absent: it drops SNPs from alignments with an ambiguous mapping, and in a
    # SELF-alignment every position is ambiguous by definition, so -C returns an empty file.
    all_sites, indels = parse_snps(_run([a.show_snps, "-l", "-r", "-T", "-H", delta]))
    sites = sites_within_pairs(all_sites, pairs)

    if a.ancestor_delta:
        anc = (rehomed_delta(a.ancestor_delta, a.reference, a.outgroup,
                             out=a.ancestor_delta + ".rehomed")
               if a.reference else a.ancestor_delta)
        diffs, intervals = ancestral_bases(
            _run([a.show_snps, "-l", "-r", "-T", "-H", anc]),
            _run([a.show_coords, "-r", "-c", "-l", "-T", anc]))
        polarise(sites, diffs, intervals)

    counts = {}
    for s in sites:
        counts[s["pair_id"]] = counts.get(s["pair_id"], 0) + 1
    for i, p in enumerate(pairs):
        p["pair_id"] = i
        p["n_diagnostic"] = counts.get(i, 0)

    write_tsv(a.out_pairs, PAIR_COLS, pairs)
    write_tsv(a.out_sites, SITE_COLS, sites)

    sys.stderr.write(
        f"[paralog_map] {len(pairs)} paralogous pair(s), {len(sites)} diagnostic site(s) "
        f"({sum(1 for s in sites if s['kind'] == 'del')} of them deletions)"
        f"{f', {indels} position(s) skipped' if indels else ''}\n")
    if a.ancestor_delta:
        by = {}
        for s in sites:
            by[s["polarity"]] = by.get(s["polarity"], 0) + 1
        sys.stderr.write(
            f"[paralog_map] polarised against the outgroup: {by.get('derived', 0)} site(s) where "
            f"the reference is ancestral, {by.get('ancestral', 0)} where it is not and the donor's "
            f"base is the ancestral one, {by.get('third', 0)} where the outgroup carries neither, "
            f"{by.get('', 0)} it cannot reach\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
