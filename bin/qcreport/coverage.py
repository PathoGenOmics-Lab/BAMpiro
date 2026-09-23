"""Deletions and callable depth along the genome, from the per-sample depth profiles.

bin/depth_profile.py reduces each sample's all-positions VCF to windows, the stretches no read
covers, and a per-gene table. What those stretches mean needs the cohort:

- a stretch almost no sample reads is a part of the reference none of these genomes has, or one
  no read can be placed on (a repeat). It is not a deletion of any sample;
- a stretch some samples read and others do not is a deletion in the others, shared by a
  lineage when several carry it and private to one sample when only it does;
- in a sample read too thinly, stretches without reads turn up by chance, so its stretches are
  not assessed at all rather than reported as deletions.

Stretches of different samples are the same deletion when each covers at least half of the
other (reciprocal overlap), measured against the stretch that opened the region: overlap alone
would chain a run of neighbouring small deletions into one, and let one long private deletion
swallow a short gap every sample shares.

The windows also give, per bin of the genome landscape, how much of each sample could be called,
which is what turns a count of SNPs per bin into a density: a bin half of which was not read
holds half the SNPs whatever the genome looks like there.

A sample mapped against two references has a profile for each, so profiles are keyed on the
sample and the reference together.
"""
from __future__ import annotations

import math
import os
from array import array

from .parsers import NBINS, to_float


def _kind(header):
    """Which of depth_profile.py's tables a column header belongs to. The gene table also has a
    mean_dp column, so it is recognised first."""
    if "gene" in header:
        return "genes"
    if "zero_frac" in header:
        return "windows"
    if "length" in header:
        return "zero"
    return None


def _read_table(path):
    """(meta, header, rows as lists) of one depth_profile.py table, or None."""
    meta, header, rows = {}, None, []
    try:
        with open(path, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                line = line.rstrip("\n")
                if line.startswith("#"):
                    k, sep, v = line[1:].strip().partition("=")
                    if sep:
                        meta[k.strip()] = v.strip()
                    continue
                if header is None:
                    header = line.split("\t")
                    continue
                if line:
                    rows.append(line.split("\t"))
    except OSError:
        return None
    return meta, header or [], rows


def parse_depth_profiles(paths, min_len=0, nbins=NBINS):
    """(profiles, gene lists) from depth_profile.py's tables, any of the three per sample.

    profiles maps (sample, reference) to what the report needs and no more, so a cohort of a
    thousand samples fits: the header values, callable positions per landscape bin (the windows
    are not kept), the stretches without reads of at least `min_len` bp, and the breadth of each
    gene as a compact array aligned with its reference's gene list (gene lists maps the reference
    to [(contig, start, end, gene, locus_tag)], the same for every sample on it).
    """
    profiles, gene_lists = {}, {}
    for path in paths or []:
        if not path or not os.path.exists(path) or os.path.basename(path).startswith("NO_FILE"):
            continue
        table = _read_table(path)
        if table is None:
            continue
        meta, header, rows = table
        sample, kind = meta.get("sample"), _kind(header)
        if not sample or not kind:
            continue
        ref = meta.get("reference") or ""
        p = profiles.setdefault((sample, ref), {"sample": sample, "reference": ref, "median_dp": 0.0,
                                                "contigs": [], "callable": None, "zero": [], "genes": None})
        if not p["contigs"]:
            p["contigs"] = _contigs(meta.get("contigs"))
            p["median_dp"] = to_float(meta.get("median_dp")) or 0.0
        col = {h: i for i, h in enumerate(header)}

        def cell(r, name):
            i = col.get(name)
            return r[i] if i is not None and i < len(r) else ""

        rows = [r for r in rows if cell(r, "start").isdigit() and cell(r, "end").isdigit()]
        if kind == "windows":
            spans = [(cell(r, "contig"), int(cell(r, "start")), int(cell(r, "end")),
                      to_float(cell(r, "callable_frac"))) for r in rows]
            spans = [s for s in spans if s[3] is not None]
            p["callable"] = _spread(spans, p["contigs"], nbins) if spans and p["contigs"] else None
        elif kind == "zero":
            p["zero"] = [(cell(r, "contig"), int(cell(r, "start")), int(cell(r, "end"))) for r in rows
                         if int(cell(r, "end")) - int(cell(r, "start")) + 1 >= min_len]
        else:
            keys = [(cell(r, "contig"), int(cell(r, "start")), int(cell(r, "end")), cell(r, "gene"),
                     cell(r, "locus_tag")) for r in rows]
            breadth = [to_float(cell(r, "breadth")) for r in rows]
            known = gene_lists.setdefault(ref, keys)
            if keys != known:              # another GFF version for the same reference: align by key
                at = {k: i for i, k in enumerate(keys)}
                breadth = [breadth[at[k]] if k in at else None for k in known]
            p["genes"] = array("f", [-1.0 if b is None else b for b in breadth])
    return profiles, gene_lists


def _contigs(text):
    """[(name, length)] from '# contigs=a:100,b:20', in the order of the reference."""
    out = []
    for item in (text or "").split(","):
        name, sep, n = item.rpartition(":")
        if sep and n.isdigit():
            out.append((name, int(n)))
    return out


def _offsets(contigs):
    """({contig: offset}, genome length): the contigs laid end to end, as the landscape does."""
    offs, total = {}, 0
    for name, n in contigs:
        offs[name] = total
        total += n
    return offs, total


def _spread(intervals, contigs, nbins=NBINS):
    """Per-bin amounts from [(contig, start, end, amount per position)], 1-based inclusive.

    Uses the landscape's own binning, bin = global position * nbins // genome length (see
    stats_to_legacy.analyze_vcf), so a value lands in the bin the SNPs of those positions do.
    """
    offs, total = _offsets(contigs)
    out = [0.0] * nbins
    if total <= 0:
        return out
    for contig, s, e, per in intervals:
        if contig not in offs or not per:
            continue
        g0, g1 = offs[contig] + s, offs[contig] + e
        b = min(nbins - 1, g0 * nbins // total)
        while g0 <= g1:
            nxt = min(g1, ((b + 1) * total - 1) // nbins)   # last position of bin b
            if nxt < g0:
                nxt = g0
            out[b] += per * (nxt - g0 + 1)
            g0, b = nxt + 1, min(nbins - 1, b + 1)
    return out


def _bin_len(contigs, nbins=NBINS):
    return _spread([(c, 1, n, 1.0) for c, n in contigs], contigs, nbins)


def snp_per_kb(snp_profile, profile, nbins=NBINS):
    """SNPs per callable kb in each bin, or None without both profiles.

    A bin with fewer than a tenth of its positions callable has no value: a few called bases
    make any count look dense or empty.
    """
    cb = (profile or {}).get("callable")
    if not snp_profile or not cb or not profile.get("contigs"):
        return None
    full = _bin_len(profile["contigs"], nbins)
    return [round(1000.0 * (snp_profile[b] or 0) / cb[b], 3)
            if cb[b] > 0 and full[b] and cb[b] >= 0.1 * full[b] else None
            for b in range(min(nbins, len(snp_profile)))]


def _overlap(a0, a1, b0, b1):
    return max(0, min(a1, b1) - max(a0, b0) + 1)


def _cluster(stretches, frac=0.5):
    """Regions of [(start, end, carrier)] sorted by start: a stretch joins the first open region
    whose founding stretch and itself each overlap by at least `frac` of their length."""
    regions, open_ = [], []
    for s, e, who in stretches:
        open_ = [r for r in open_ if r["e0"] >= s]
        home = None
        for r in open_:
            ov = _overlap(r["s0"], r["e0"], s, e)
            if ov >= frac * (r["e0"] - r["s0"] + 1) and ov >= frac * (e - s + 1):
                home = r
                break
        if home is None:
            home = {"s0": s, "e0": e, "members": {}}
            regions.append(home)
            open_.append(home)
        home["members"].setdefault(who, []).append((s, e))
    return regions


def _median(v):
    v = sorted(v)
    return v[len(v) // 2]


def build_coverage(profiles, gene_lists=None, min_len=200, min_depth=10.0, cohort_frac=0.9,
                   gene_breadth=0.5, nbins=NBINS):
    """(report section, {(sample, reference): per-bin deleted fraction}, rows for the TSV).

    Stretches without reads of at least `min_len` bp, in samples whose genome-wide median depth
    is at least `min_depth`, are grouped across samples into regions (see _cluster) and classed
    by how many of the assessed samples on the same reference lack them:

      cohort   at least `cohort_frac` of them, and at least two: not a deletion of anyone
      shared   two or more, fewer than that: a deletion several samples carry
      private  one sample, with at least one other assessed sample reading the region
      alone    the only assessed sample on its reference, so nothing to compare it with

    Genes are named for a region when a carrier reads less than `gene_breadth` of them.
    """
    gene_lists = gene_lists or {}
    by_ref = {}
    for key, p in profiles.items():
        by_ref.setdefault(p.get("reference") or ",".join(c for c, _ in p.get("contigs", [])), []).append(key)

    # The genes each profile barely reads, found once: naming a region's genes then looks at a
    # handful of genes per carrier instead of every gene of the reference.
    thin_genes = {}
    for key, p in profiles.items():
        genes, b = gene_lists.get(p.get("reference") or "", []), p.get("genes")
        thin_genes[key] = ([genes[i] for i in range(min(len(genes), len(b))) if 0 <= b[i] < gene_breadth]
                           if b is not None else [])

    regions, not_assessed = [], []
    for ref, keys in sorted(by_ref.items()):
        assessed = [k for k in keys if profiles[k]["median_dp"] >= min_depth]
        not_assessed += [k[0] for k in keys if k not in assessed]
        n = len(assessed)
        by_contig = {}
        for k in assessed:
            for c, s, e in profiles[k]["zero"]:
                if e - s + 1 >= min_len:
                    by_contig.setdefault(c, []).append((s, e, k))
        for contig, stretches in sorted(by_contig.items()):
            for r in _cluster(sorted(stretches)):
                k = len(r["members"])
                if n < 2:
                    cls = "alone"
                elif k >= max(2, math.ceil(cohort_frac * n)):
                    cls = "cohort"
                elif k >= 2:
                    cls = "shared"
                else:
                    cls = "private"
                own = {who: (min(a for a, _ in iv), max(b for _, b in iv)) for who, iv in r["members"].items()}
                start, end = _median([a for a, _ in own.values()]), _median([b for _, b in own.values()])
                names = {}
                for who in r["members"]:
                    for gc, gs, ge, gene, locus in thin_genes.get(who, []):
                        if gc == contig and _overlap(gs, ge, start, end):
                            names.setdefault(gene, locus)
                regions.append({
                    "ref": ref, "contig": contig, "start": start, "end": end, "len": end - start + 1,
                    "cls": cls, "n": k, "of": n, "genes": sorted(names),
                    "samples": sorted([who[0], a, b] for who, (a, b) in own.items()),
                    "_members": r["members"],
                })

    # Per-bin share of each profile's positions inside its own stretches of the regions it
    # carries. A stretch nearly nobody reads is nobody's deletion, and a sample alone on its
    # reference has nothing to be compared with, so neither lights up the track.
    per = {}
    for reg in regions:
        if reg["cls"] in ("cohort", "alone"):
            continue
        for who, iv in reg["_members"].items():
            per.setdefault(who, []).extend((reg["contig"], s, e, 1.0) for s, e in iv)
    tracks = {}
    for key, p in profiles.items():
        if p["median_dp"] < min_depth or not p.get("contigs"):
            continue
        full = _bin_len(p["contigs"], nbins)
        lost = _spread(per.get(key, []), p["contigs"], nbins)
        tracks[key] = [round(min(1.0, lost[b] / full[b]), 3) if full[b] else 0.0 for b in range(nbins)]

    for r in regions:
        del r["_members"]
    lost_genes = _lost_genes(profiles, gene_lists, by_ref, min_depth, gene_breadth)
    regions.sort(key=lambda r: (r["cls"] == "cohort", -r["n"] if r["cls"] == "shared" else 0,
                                r["ref"], r["contig"], r["start"]))
    for i, r in enumerate(regions, 1):
        r["id"] = f"D{i}"
    # Each reference's contigs in order, so the report can place a region on the landscape, which
    # lays the contigs of a sample's genome end to end.
    refs = {ref: [list(c) for c in profiles[keys[0]].get("contigs", [])] for ref, keys in by_ref.items()}
    section = {"min_len": min_len, "min_depth": min_depth, "cohort_frac": cohort_frac, "refs": refs,
               "regions": regions, "genes_lost": lost_genes, "not_assessed": sorted(set(not_assessed)),
               "n_assessed": sum(1 for p in profiles.values() if p["median_dp"] >= min_depth)}
    return section, tracks, regions


def _lost_genes(profiles, gene_lists, by_ref, min_depth, gene_breadth, cohort_breadth=0.9):
    """Genes a sample reads less than `gene_breadth` of while at least half of the assessed
    samples on its reference read `cohort_breadth` of them: gene loss, rather than a gene
    nobody here carries."""
    out = []
    for ref, keys in sorted(by_ref.items()):
        genes = gene_lists.get(ref) or []
        assessed = [k for k in keys if profiles[k]["median_dp"] >= min_depth and profiles[k].get("genes")]
        if len(assessed) < 2 or not genes:
            continue
        for i, (contig, s, e, gene, locus) in enumerate(genes):
            vals = [(k[0], profiles[k]["genes"][i]) for k in assessed
                    if i < len(profiles[k]["genes"]) and profiles[k]["genes"][i] >= 0]
            lost = sorted([sid, round(b, 3)] for sid, b in vals if b < gene_breadth)
            read = sum(1 for _, b in vals if b >= cohort_breadth)
            if lost and read >= 0.5 * len(vals):
                out.append({"ref": ref, "contig": contig, "start": s, "end": e, "gene": gene,
                            "locus": locus, "samples": lost, "of": len(vals)})
    return out


def write_deletions(path, regions):
    """The regions as a TSV: one row per region, each carrier with its own stretch."""
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("region\treference\tcontig\tstart\tend\tlength\tclass\tn_samples\tn_assessed\t"
                 "samples\tgenes\n")
        for r in regions:
            fh.write("\t".join(str(v) for v in (
                r["id"], r["ref"], r["contig"], r["start"], r["end"], r["len"], r["cls"], r["n"], r["of"],
                ",".join(f"{sid}:{s}-{e}" for sid, s, e in r["samples"]), ",".join(r["genes"]))) + "\n")
