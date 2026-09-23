"""Deletions and callable depth along the genome, from the per-sample depth profiles.

bin/depth_profile.py reduces each sample's all-positions VCF to windows, the stretches no read
covers, and a per-gene table. What those stretches mean needs the cohort:

- a stretch no sample reads is a part of the reference none of these genomes has, or one no
  read can be placed on (a repeat). It is not a deletion of any sample;
- a stretch some samples read and others do not is a deletion in the others, shared by a
  lineage when several carry it and private to one sample when only it does;
- in a sample read too thinly, stretches without reads turn up by chance, so its stretches are
  not assessed at all rather than reported as deletions.

The windows also give, per bin of the genome landscape, how much of each sample could be called,
which is what turns a count of SNPs per bin into a density: a bin half of which was not read
holds half the SNPs whatever the genome looks like there.
"""
from __future__ import annotations

import math
import os

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


def parse_depth_profiles(paths):
    """{sample: profile} from depth_profile.py's tables, any of the three per sample.

    The table is recognised by its columns and the sample by its '# sample=' line, so neither
    depends on a file name. A profile holds the header values (reference, median_dp, contigs...)
    and the rows of whichever tables were given.
    """
    out = {}
    for path in paths or []:
        if not path or not os.path.exists(path) or os.path.basename(path).startswith("NO_FILE"):
            continue
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
                        rows.append(dict(zip(header, line.split("\t"))))
        except OSError:
            continue
        sample = meta.get("sample")
        kind = _kind(header or [])
        if not sample or not kind:
            continue
        p = out.setdefault(sample, {"windows": [], "zero": [], "genes": []})
        for k in ("reference", "median_dp", "mean_dp", "genome_len", "zero_frac", "callable_frac",
                  "min_dp", "window", "min_run", "contigs"):
            if k in meta and k not in p:
                p[k] = meta[k]
        if kind == "windows":
            p["windows"] = [(r["contig"], int(r["start"]), int(r["end"]), to_float(r.get("mean_dp")),
                             to_float(r.get("zero_frac")), to_float(r.get("callable_frac")))
                            for r in rows if r.get("start", "").isdigit() and r.get("end", "").isdigit()]
        elif kind == "zero":
            p["zero"] = [(r["contig"], int(r["start"]), int(r["end"]))
                         for r in rows if r.get("start", "").isdigit() and r.get("end", "").isdigit()]
        else:
            p["genes"] = [{"contig": r["contig"], "start": int(r["start"]), "end": int(r["end"]),
                           "gene": r.get("gene", ""), "locus": r.get("locus_tag", ""),
                           "breadth": to_float(r.get("breadth")), "rel": to_float(r.get("rel_depth"))}
                          for r in rows if r.get("start", "").isdigit() and r.get("end", "").isdigit()]
    for p in out.values():
        p["median_dp"] = to_float(p.get("median_dp")) or 0.0
        p["contigs"] = _contigs(p.get("contigs"))
    return out


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


def callable_bp(profile, nbins=NBINS):
    """Per-bin number of positions deep enough for the consensus to call, or None."""
    rows = [(c, s, e, f) for c, s, e, _, _, f in profile.get("windows", []) if f is not None]
    if not rows or not profile.get("contigs"):
        return None
    return _spread(rows, profile["contigs"], nbins)


def _bin_len(contigs, nbins=NBINS):
    return _spread([(c, 1, n, 1.0) for c, n in contigs], contigs, nbins)


def snp_per_kb(snp_profile, profile, nbins=NBINS):
    """SNPs per callable kb in each bin, or None without both profiles.

    A bin with fewer than a tenth of its positions callable has no value: a few called bases
    make any count look dense or empty.
    """
    cb = callable_bp(profile, nbins)
    if not snp_profile or cb is None:
        return None
    full = _bin_len(profile["contigs"], nbins)
    return [round(1000.0 * (snp_profile[b] or 0) / cb[b], 3)
            if cb[b] > 0 and full[b] and cb[b] >= 0.1 * full[b] else None for b in range(nbins)]


def _overlap(a0, a1, b0, b1):
    return max(0, min(a1, b1) - max(a0, b0) + 1)


def build_coverage(profiles, min_len=200, min_depth=10.0, cohort_frac=0.9, gene_breadth=0.5,
                   nbins=NBINS):
    """(report section, {sample: per-bin deleted fraction}, rows for the deletions TSV).

    Stretches without reads of at least `min_len` bp, in samples whose genome-wide median depth
    is at least `min_depth`, are merged across samples into regions and classed by how many of
    the assessed samples on the same reference lack them:

      cohort   at least `cohort_frac` of them, and at least two: not a deletion of anyone
      shared   two or more, fewer than that: a deletion several samples carry
      private  one sample, with at least one other assessed sample reading the region
      alone    the only assessed sample on its reference, so nothing to compare it with

    Genes are named for a region when a carrier reads less than `gene_breadth` of them.
    """
    by_ref = {}
    for sid, p in profiles.items():
        by_ref.setdefault(p.get("reference") or ",".join(c for c, _ in p.get("contigs", [])), []).append(sid)

    regions, tracks, not_assessed = [], {}, []
    for ref, samples in sorted(by_ref.items()):
        assessed = [s for s in samples if profiles[s]["median_dp"] >= min_depth]
        not_assessed += [s for s in samples if s not in assessed]
        runs = sorted((c, s, e, sid) for sid in assessed for c, s, e in profiles[sid]["zero"]
                      if e - s + 1 >= min_len)
        merged = []
        for c, s, e, sid in runs:
            if merged and merged[-1]["contig"] == c and s <= merged[-1]["end"]:
                r = merged[-1]
                r["end"] = max(r["end"], e)
                r["carriers"].setdefault(sid, []).append((s, e))
            else:
                merged.append({"contig": c, "start": s, "end": e, "carriers": {sid: [(s, e)]}})
        n = len(assessed)
        for r in merged:
            k = len(r["carriers"])
            if n < 2:
                cls = "alone"
            elif k >= max(2, math.ceil(cohort_frac * n)):
                cls = "cohort"
            elif k >= 2:
                cls = "shared"
            else:
                cls = "private"
            genes = {}
            for sid in r["carriers"]:
                for g in profiles[sid]["genes"]:
                    if g["contig"] == r["contig"] and _overlap(g["start"], g["end"], r["start"], r["end"]) \
                            and g["breadth"] is not None and g["breadth"] < gene_breadth:
                        genes.setdefault(g["gene"], g["locus"])
            regions.append({
                "ref": ref, "contig": r["contig"], "start": r["start"], "end": r["end"],
                "len": r["end"] - r["start"] + 1, "cls": cls, "n": k, "of": n,
                "genes": sorted(genes),
                "samples": sorted([sid, min(a for a, _ in iv), max(b for _, b in iv)]
                                  for sid, iv in r["carriers"].items()),
            })

    # Per-bin share of each sample's positions inside a deletion it carries (cohort-wide gaps
    # excluded, since they are nobody's deletion).
    per_sample = {}
    for reg in regions:
        if reg["cls"] == "cohort":
            continue
        for sid, s, e in reg["samples"]:
            per_sample.setdefault(sid, []).append((reg["contig"], s, e, 1.0))
    for sid, ivs in per_sample.items():
        contigs = profiles[sid].get("contigs") or []
        full = _bin_len(contigs, nbins)
        lost = _spread(ivs, contigs, nbins)
        tracks[sid] = [round(min(1.0, lost[b] / full[b]), 3) if full[b] else 0.0 for b in range(nbins)]
    for sid in profiles:
        if sid not in tracks and profiles[sid]["median_dp"] >= min_depth and profiles[sid].get("contigs"):
            tracks[sid] = [0.0] * nbins

    lost_genes = _lost_genes(profiles, by_ref, min_depth, gene_breadth)
    regions.sort(key=lambda r: (r["cls"] == "cohort", -r["n"] if r["cls"] == "shared" else 0,
                                r["ref"], r["contig"], r["start"]))
    for i, r in enumerate(regions, 1):
        r["id"] = f"D{i}"
    # Each reference's contigs in order, so the report can place a region on the landscape, which
    # lays the contigs of a sample's genome end to end.
    refs = {ref: [list(c) for c in profiles[samples[0]].get("contigs", [])] for ref, samples in by_ref.items()}
    section = {"min_len": min_len, "min_depth": min_depth, "cohort_frac": cohort_frac, "refs": refs,
               "regions": regions, "genes_lost": lost_genes, "not_assessed": sorted(not_assessed),
               "n_assessed": sum(1 for p in profiles.values() if p["median_dp"] >= min_depth)}
    return section, tracks, regions


def _lost_genes(profiles, by_ref, min_depth, gene_breadth, cohort_breadth=0.9):
    """Genes a sample reads less than `gene_breadth` of while at least half of the assessed
    samples on its reference read `cohort_breadth` of them: gene loss, rather than a gene
    nobody here carries."""
    out = []
    for ref, samples in sorted(by_ref.items()):
        assessed = [s for s in samples if profiles[s]["median_dp"] >= min_depth and profiles[s]["genes"]]
        if len(assessed) < 2:
            continue
        table = {}
        for sid in assessed:
            for g in profiles[sid]["genes"]:
                if g["breadth"] is not None:
                    table.setdefault((g["contig"], g["start"], g["end"], g["gene"], g["locus"]), {})[sid] = g["breadth"]
        for (contig, s, e, gene, locus), vals in sorted(table.items()):
            lost = sorted([sid, round(b, 3)] for sid, b in vals.items() if b < gene_breadth)
            read = sum(1 for b in vals.values() if b >= cohort_breadth)
            if lost and read >= 0.5 * len(vals):
                out.append({"ref": ref, "contig": contig, "start": s, "end": e, "gene": gene,
                            "locus": locus, "samples": lost, "of": len(vals)})
    return out


def write_deletions(path, regions):
    """The regions as a TSV: one row per region, its carriers with their own coordinates."""
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("region\treference\tcontig\tstart\tend\tlength\tclass\tn_samples\tn_assessed\t"
                 "samples\tgenes\n")
        for r in regions:
            fh.write("\t".join(str(v) for v in (
                r["id"], r["ref"], r["contig"], r["start"], r["end"], r["len"], r["cls"], r["n"], r["of"],
                ",".join(f"{sid}:{s}-{e}" for sid, s, e in r["samples"]), ",".join(r["genes"]))) + "\n")
