"""How low an allele frequency can be trusted: the variant calls below fixation.

Every sample: how many calls sit below fixation (allele fraction under 0.9), their spread of
fractions, and how many of them rest on three alternate reads or fewer, where a sequencing error
and a real minority look the same.

When the samplesheet says which samples are libraries of the same DNA (a column such as dna_id,
extract or biosample), the question has an empirical answer: a real minority variant is in the
DNA, so another library of it calls it too, while an error is not reproduced. The share of calls
one library makes that another library of the same DNA also makes, in each band of allele
fraction, is where the noise ends. Two samples that share a FASTQ file share reads, so their
agreement proves nothing and they are not compared; a merged sample and the runs it was merged
from are the usual case. A call is only counted as missing from the other library where that
library was read at the site (the SNP matrix's depth); a site it did not read says nothing.
"""
from __future__ import annotations

import csv
import os
import re

FIXED = 0.9
EDGES = [0.0, 0.1, 0.2, 0.3, 0.5, 0.7, FIXED]
_REP_RE = re.compile(r"^(dna[_ ]?id|dna|extract(ion)?([_ ]?id)?|biosample([_ ]?id)?|specimen([_ ]?id)?|"
                     r"isolate[_ ]?id|library[_ ]?of|replicate[_ ]?of|tech(nical)?[_ ]?rep(licate)?([_ ]?of)?)$", re.I)
_SAMPLE_RE = re.compile(r"^(sample[_ ]?id|sample|id|name)$", re.I)


def _bin(af):
    """The band of EDGES an allele fraction falls in, (lo, hi]; None at or above fixation."""
    for i in range(len(EDGES) - 1):
        if EDGES[i] < af <= EDGES[i + 1]:
            return i
    return None


def replicate_sets(path):
    """(column, {sample: value}, {sample: {fastq paths}}) from the samplesheet, or (None, {}, {}).

    The column is the first one named like a DNA or extract identifier. FASTQ paths come from
    r1/r2 and are kept per sample, since a merged sample lists every run it was merged from.
    """
    if not path or not os.path.exists(path):
        return None, {}, {}
    try:
        with open(path, encoding="utf-8", errors="replace") as fh:
            rows = list(csv.DictReader((line for line in fh if not line.startswith("#")), delimiter="\t"))
    except OSError:
        return None, {}, {}
    if not rows:
        return None, {}, {}
    header = list(rows[0].keys())
    scol = next((h for h in header if h and _SAMPLE_RE.match(h.strip())), header[0])
    rcol = next((h for h in header if h and _REP_RE.match(h.strip())), None)
    value, reads = {}, {}
    for r in rows:
        s = (r.get(scol) or "").strip()
        if not s:
            continue
        for k in ("r1", "r2", "R1", "R2", "fastq_1", "fastq_2"):
            v = (r.get(k) or "").strip()
            if v:
                reads.setdefault(s, set()).add(v)
        if rcol and (r.get(rcol) or "").strip():
            value.setdefault(s, r[rcol].strip())
    return rcol, value, reads


def replicate_pairs(value, reads, samples, ref=None):
    """([(a, b)], n left out for sharing reads): independent libraries of the same DNA.

    With `ref` ({sample: reference}), libraries of one DNA mapped against different references
    are not paired either: their coordinates do not correspond, and the pairing is itself a
    sign one of them went to the wrong genome.
    """
    by = {}
    for s, v in value.items():
        if s in samples:
            by.setdefault(v, []).append(s)
    pairs, shared = [], 0
    for members in by.values():
        members.sort()
        for i, a in enumerate(members):
            for b in members[i + 1:]:
                if ref and ref.get(a) != ref.get(b):
                    continue
                if reads.get(a, set()) & reads.get(b, set()):
                    shared += 1
                else:
                    pairs.append((a, b))
    return pairs, shared


def needed_cells(pairs, variants, fixed=FIXED):
    """{site: {samples}}: each library at the sites the other of its pair called."""
    need = {}
    for a, b in pairs:
        for x, y in ((a, b), (b, a)):
            for site, v in variants.get(x, {}).items():
                if v.get("af") and site not in variants.get(y, {}):
                    need.setdefault(site, set()).add(y)
    return need


def build_minority(variants, pairs=(), shared=0, cells=None, column=None, min_dp=7, fixed=FIXED,
                   checked=True):
    """The payload section, or None without variant calls. Without the SNP matrix (`checked`
    False) nothing tells a call the other library missed from a site it did not read, so the
    replicate comparison is left out rather than biased."""
    if not variants:
        return None
    samples, hist_all = {}, [0] * (len(EDGES) - 1)
    low_all = total_all = 0
    for s, calls in variants.items():
        hist, alts, low = [0] * (len(EDGES) - 1), [], 0
        for v in calls.values():
            af = v.get("af")
            if af is None or af <= 0 or af >= fixed:
                continue
            b = _bin(af)
            if b is None:
                continue
            hist[b] += 1
            if v.get("dp"):
                alt = round(af * v["dp"])
                alts.append(alt)
                low += alt <= 3
        n = sum(hist)
        for i, c in enumerate(hist):
            hist_all[i] += c
        low_all += low
        total_all += len(alts)
        alts.sort()
        samples[s] = {"n": n, "hist": hist, "le3": low, "with_dp": len(alts),
                      "median_alt": alts[len(alts) // 2] if alts else None}
    section = {"fixed": fixed, "edges": EDGES, "hist": hist_all, "samples": samples,
               "le3": low_all, "with_dp": total_all, "replicates": None}
    if column:
        section["replicates"] = (_reproducibility(variants, pairs, shared, {} if cells is None else cells,
                                                  column, min_dp, fixed)
                                 if checked else {"column": column, "pairs": len(pairs), "shared": shared,
                                                  "tested": None})
    return section


def _reproducibility(variants, pairs, shared, cells, column, min_dp, fixed):
    """Per band of allele fraction, the calls of one library the other library of its pair
    also makes, among those the other library was read at; and the same for fixed calls."""
    tested, repro = [0] * (len(EDGES) - 1), [0] * (len(EDGES) - 1)
    fx_t = fx_r = 0
    for a, b in pairs:
        for x, y in ((a, b), (b, a)):
            vy = variants.get(y, {})
            for site, v in variants.get(x, {}).items():
                af = v.get("af")
                if not af:
                    continue
                if site in vy and (vy[site].get("af") or 0) > 0:
                    hit = True
                else:
                    cell = cells.get((site, y))
                    if cell is None or cell[1] is None or cell[1] < min_dp:
                        continue          # the other library was not read there: nothing to learn
                    hit = bool(cell[0])
                if af >= fixed:
                    fx_t += 1
                    fx_r += hit
                    continue
                i = _bin(af)
                if i is not None:
                    tested[i] += 1
                    repro[i] += hit
    return {"column": column, "pairs": len(pairs), "shared": shared, "tested": tested,
            "reproduced": repro, "fixed_tested": fx_t, "fixed_reproduced": fx_r}
