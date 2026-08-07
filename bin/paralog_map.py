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
import subprocess
import sys

# show-coords -r -c -l -T: S1 E1 S2 E2 LEN1 LEN2 %IDY LENR LENQ COVR COVQ TAG1 TAG2
COORDS_FIELDS = 13
# show-snps -l -r -T -H: P1 SUB1 SUB2 P2 BUFF DIST R Q LENR LENQ FRM1 FRM2 TAG1 TAG2
SNPS_FIELDS = 14


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

    An indel is reported with '.' on one side. Those are kept out of the diagnostic set: a read
    carrying one is still evidence, but the position no longer maps one-to-one between the copies,
    which the tract caller assumes. They are counted so the caller can say how many it skipped.
    """
    sites, indels = [], 0
    for line in text.splitlines():
        f = line.rstrip("\n").split("\t")
        if len(f) < SNPS_FIELDS or not f[0].isdigit():
            continue
        acc_base, don_base = f[1].upper(), f[2].upper()
        if acc_base == "." or don_base == "." or acc_base == don_base:
            indels += 1
            continue
        sites.append({
            "acceptor": f[12], "acc_pos": int(f[0]), "acc_base": acc_base,
            "donor": f[13], "don_pos": int(f[3]), "don_base": don_base,
            "strand": "-" if f[11] == "-1" else "+",
        })
    return sites, indels


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
    chosen = {}
    for site in sites:
        for i, p in enumerate(pairs):
            if site["acceptor"] != p["acceptor"] or site["donor"] != p["donor"]:
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
SITE_COLS = ["pair_id", "acceptor", "acc_pos", "acc_base", "donor", "don_pos", "don_base", "strand"]


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
    p.add_argument("--show-coords", default="show-coords")
    p.add_argument("--show-snps", default="show-snps")
    return p.parse_args(argv)


def main(argv=None) -> int:
    a = parse_args(argv)

    pairs = parse_coords(_run([a.show_coords, "-r", "-c", "-l", "-T", a.delta]),
                         min_identity=a.min_identity, min_length=a.min_length)
    # -C is deliberately absent: it drops SNPs from alignments with an ambiguous mapping, and in a
    # SELF-alignment every position is ambiguous by definition, so -C returns an empty file.
    all_sites, indels = parse_snps(_run([a.show_snps, "-l", "-r", "-T", "-H", a.delta]))
    sites = sites_within_pairs(all_sites, pairs)

    counts = {}
    for s in sites:
        counts[s["pair_id"]] = counts.get(s["pair_id"], 0) + 1
    for i, p in enumerate(pairs):
        p["pair_id"] = i
        p["n_diagnostic"] = counts.get(i, 0)

    write_tsv(a.out_pairs, PAIR_COLS, pairs)
    write_tsv(a.out_sites, SITE_COLS, sites)

    sys.stderr.write(
        f"[paralog_map] {len(pairs)} paralogous pair(s), {len(sites)} diagnostic site(s)"
        f"{f', {indels} indel/ambiguous position(s) skipped' if indels else ''}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
