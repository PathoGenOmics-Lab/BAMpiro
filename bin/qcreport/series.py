"""What each series gained since its first time point.

A series is the samples the samplesheet puts in one group (a patient, a passage line) at several
time points. For every later sample, the SNPs it carries fixed (allele fraction 0.9 or more) are
read against the series' first time point:

  new      absent there, at a site read deeply enough to say so
  risen    present there already, below fixation (a minority that took over)
  unknown  the first time point was not read there, so nothing says whether it was absent
  lost     fixed at the first time point, read and absent in the later sample: a reversion, or a
           sample that does not belong to the series

Telling 'absent' from 'not read' needs the depth of the sites a sample has no call at, which is
what the SNP matrix holds. Without it every site not called at the start counts as new, and the
section says so rather than passing that off as checked. A site counts as read above
--consensus_min_dp reads, where the consensus stops leaving it uncalled.

The first time point is the earliest one with a sample the QC keeps and places in its series:
a failed or swapped first sample would make every later sample look as if it had gained
everything that sample lacks. A sample whose time cannot be read as a number has no place in
the series and is left out of it.
"""
from __future__ import annotations

import os

FIXED = 0.9


def _order(meta, samples):
    return sorted(samples, key=lambda s: (meta[s].get("tnum") is None,
                                          meta[s].get("tnum") if meta[s].get("tnum") is not None else 0,
                                          str(meta[s].get("time") or "")))


def series_groups(meta, variants, excluded=()):
    """{group: (samples from the first time point on, in time order; first-time-point samples)}
    for groups with two or more times among the samples not `excluded`."""
    groups = {}
    for s, md in (meta or {}).items():
        if md.get("group") and s in variants:
            groups.setdefault(md["group"], []).append(s)
    out = {}
    for g, samples in groups.items():
        timed = [s for s in samples if meta[s].get("tnum") is not None]
        samples = _order(meta, timed if timed else samples)
        kept = [s for s in samples if s not in excluded]
        times = [meta[s].get("time") for s in kept]
        if len(set(times)) < 2:
            continue
        base = [s for s in kept if meta[s].get("time") == times[0]]
        start = samples.index(base[0])
        out[g] = ([s for s in samples[start:] if s in base or meta[s].get("time") != times[0]], base)
    return out


def needed_cells(meta, variants, fixed=FIXED, excluded=()):
    """{site: {samples}} whose matrix cells the comparison reads: the first time point at every
    site a later sample carries fixed, and every later sample at the sites fixed at the start."""
    need = {}
    for samples, base in series_groups(meta, variants, excluded).values():
        later = [s for s in samples if s not in base]
        for s in later:
            for site, v in variants[s].items():
                if v.get("af") is not None and v["af"] >= fixed:
                    need.setdefault(site, set()).update(base)
        for b in base:
            for site, v in variants[b].items():
                if v.get("af") is not None and v["af"] >= fixed:
                    need.setdefault(site, set()).update(later)
    return need


def read_matrix_cells(path, need):
    """{(site, sample): (af, dp)} of the SNP matrix TSV (build_snp_matrix.py) for the cells in
    `need`. af and dp are None where the cell is empty. Reads only the rows it needs."""
    out = {}
    if not need or not path or not os.path.exists(path) or os.path.basename(path).startswith("NO_FILE"):
        return out
    try:
        with open(path, encoding="utf-8", errors="replace") as fh:
            header = fh.readline().rstrip("\n").split("\t")
            col = {h[:-3]: i for i, h in enumerate(header) if h.endswith("|AF")}
            for line in fh:
                parts = line.split("\t", 3)
                if len(parts) < 4:
                    continue
                site = f"{parts[1]}:{parts[2]}"
                wanted = need.get(site)
                if not wanted:
                    continue
                c = line.rstrip("\n").split("\t")
                for s in wanted:
                    i = col.get(s)
                    if i is None or i + 1 >= len(c):
                        continue
                    af = float(c[i]) if c[i] not in ("", "NA") else None
                    dp = int(float(c[i + 1])) if c[i + 1] not in ("", "NA") else None
                    out[(site, s)] = (af, dp)
    except (OSError, ValueError):
        return {}
    return out


def _state(site, sample, variants, cells, min_dp, fixed):
    """'fixed', 'minority', 'absent', 'unread' or 'unchecked' for one sample at one site."""
    v = variants.get(sample, {}).get(site)
    if v is not None and v.get("af") is not None:
        return "fixed" if v["af"] >= fixed else ("minority" if v["af"] > 0 else "absent")
    cell = cells.get((site, sample))
    if cell is None:
        return "unchecked"
    af, dp = cell
    if af is not None and af > 0:
        return "fixed" if af >= fixed else "minority"
    if dp is None:
        return "unchecked"
    return "absent" if dp > min_dp else "unread"


_RANK = ("fixed", "minority", "absent", "unread", "unchecked")


def build_series(meta, variants, cells, sample_meta=None, min_dp=7, fixed=FIXED, checked=True,
                 excluded=()):
    """The payload section, or None without a series of two or more time points. `excluded`
    are the samples that cannot be a series' first time point (failed, or placed outside it)."""
    groups = series_groups(meta, variants, excluded)
    if not groups:
        return None
    tx = (sample_meta or {}).get("tx_field")
    rows_meta = (sample_meta or {}).get("rows") or {}
    out, sites = [], {}
    for g, (samples, base) in sorted(groups.items()):
        memo = {}

        def start(site):
            if site not in memo:
                memo[site] = min((_state(site, b, variants, cells, min_dp, fixed) for b in base),
                                 key=_RANK.index)
            return memo[site]
        base_fixed = {site for b in base for site, v in variants[b].items()
                      if v.get("af") is not None and v["af"] >= fixed}
        rows = []
        for s in samples:
            if s in base:
                continue
            counts = {"new": [], "risen": [], "unknown": []}
            for site, v in variants[s].items():
                if v.get("af") is None or v["af"] < fixed:
                    continue
                st = start(site)
                if st == "fixed":
                    continue
                # Without the matrix, nothing tells a site absent at the start from one not read
                # there: it counts as new, and the section carries checked=False to say so.
                new = st == "absent" or (st == "unchecked" and not checked)
                key = "new" if new else ("risen" if st == "minority" else "unknown")
                counts[key].append(site)
            lost = [site for site in base_fixed
                    if _state(site, s, variants, cells, min_dp, fixed) == "absent"]
            for site in counts["new"] + counts["risen"] + lost:
                v = next((variants[x][site] for x in samples + base if site in variants.get(x, {})), {})
                sites.setdefault(site, {"gene": v.get("gene", ""), "aa": v.get("aa", ""),
                                        "eff": v.get("eff", ""), "alt": v.get("alt", "")})
            rows.append({"s": s, "time": meta[s].get("time"), "tnum": meta[s].get("tnum"),
                         "new": sorted(counts["new"]), "risen": sorted(counts["risen"]),
                         "unknown": len(counts["unknown"]), "lost": sorted(lost)})
        treat = [(rows_meta.get(s) or {}).get(tx) for s in samples] if tx else []
        treat = [t for t in treat if t]
        out.append({"group": g, "tx": max(set(treat), key=treat.count) if treat else None,
                    "base": base, "t0": meta[base[0]].get("time"), "tnum0": meta[base[0]].get("tnum"),
                    "rows": rows})
    return {"groups": out, "sites": sites, "fixed": fixed, "min_dp": min_dp, "checked": checked,
            "tx_field": tx}
