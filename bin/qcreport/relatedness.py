"""How close the samples are to each other, from the pairwise SNP distances (bin/snp_distances.py).

The distances are between consensus sequences, over the positions both samples called, so a gap
never reads as a difference. Two uses here:

- the relatedness page: every pair of samples on a reference, which the report clusters at a
  threshold the reader can move;
- GROUP_MISMATCH: a sample the samplesheet puts in a group (a patient, a line, an outbreak) but
  that sits further than the threshold from every other member of it, in a group whose members
  are normally within it of each other. That is a swapped or mislabelled sample, a contaminated
  or mixed one, or a reinfection, and each of them changes how the sample should be read. The
  sample nearest to it, and its group, go with the flag: a swap usually sits next to a sample of
  another group, a contaminant next to nobody.

The check only runs in groups that are genetically tight to begin with: a column like 'hospital'
also matches the samplesheet's group names, and its members need not be related at all. It needs
three members, since with two there is no telling which one is out of place.
"""
from __future__ import annotations

import os


def parse_pairs(path):
    """{"ref": {sample: reference}, "variable": {reference: n}, "pairs": {(a, b): (snps, compared)}}."""
    out = {"ref": {}, "variable": {}, "pairs": {}}
    if not path or not os.path.exists(path) or os.path.basename(path).startswith("NO_FILE"):
        return out
    try:
        with open(path, encoding="utf-8", errors="replace") as fh:
            header = fh.readline().rstrip("\n").split("\t")
            col = {c: i for i, c in enumerate(header)}
            need = ("sample_a", "sample_b", "reference", "snps", "compared", "variable_sites")
            if not all(c in col for c in need):
                return out
            for line in fh:
                c = line.rstrip("\n").split("\t")
                if len(c) < len(header):
                    continue
                try:
                    snps, comp, var = int(c[col["snps"]]), int(c[col["compared"]]), int(c[col["variable_sites"]])
                except ValueError:
                    continue
                a, b, ref = c[col["sample_a"]], c[col["sample_b"]], c[col["reference"]]
                out["ref"][a] = out["ref"][b] = ref
                out["variable"][ref] = var
                out["pairs"][(a, b)] = (snps, comp)
    except OSError:
        return out
    return out


def _dist(pairs, a, b):
    return pairs.get((a, b)) or pairs.get((b, a))


def build_relatedness(dist, groups=None, threshold=12, min_compared=0.5):
    """The payload section for the relatedness page, or None without distances.

    Per reference: its samples, and the upper triangle of the SNP matrix and of the share of the
    reference's variable positions both samples called (percent), row by row.
    """
    if not dist["pairs"]:
        return None
    by_ref = {}
    for s, ref in dist["ref"].items():
        by_ref.setdefault(ref, []).append(s)
    refs = {}
    for ref, samples in sorted(by_ref.items()):
        samples.sort()
        var = dist["variable"].get(ref) or 0
        snps, cmp = [], []
        for i, a in enumerate(samples):
            for b in samples[i + 1:]:
                d = _dist(dist["pairs"], a, b)
                snps.append(d[0] if d else None)
                cmp.append(round(100.0 * d[1] / var) if d and var else None)
        refs[ref] = {"samples": samples, "variable": var, "snps": snps, "cmp": cmp}
    return {"threshold": threshold, "min_compared": min_compared, "refs": refs,
            "group": {s: g for s, g in (groups or {}).items() if g and s in dist["ref"]}}


def group_mismatches(dist, groups, threshold=12, min_compared=0.5):
    """{sample: detail} of the samples far from the rest of their own group.

    A pair counts only when both samples called at least `min_compared` of their reference's
    variable positions: a thinly called sample is close to everyone. A group is judged only when
    it has three or more members with a distance and its members' median distance to their
    nearest groupmate is within `threshold`. A member of it is then flagged when its nearest
    groupmate is further than `threshold`. The detail names that groupmate, and the nearest
    sample of any other group with its group (None when there is none on the reference).
    """
    pairs, var = dist["pairs"], dist["variable"]
    groups = {s: g for s, g in (groups or {}).items() if g and s in dist["ref"]}
    if not pairs or not groups:
        return {}

    def ok(a, b):
        d = _dist(pairs, a, b)
        v = var.get(dist["ref"].get(a)) or 0
        return d if d and v and d[1] >= min_compared * v else None

    near_in, near_out = {}, {}
    samples = sorted(dist["ref"])
    for a in samples:
        if a not in groups:
            continue
        for b in samples:
            if b == a or dist["ref"].get(b) != dist["ref"].get(a):
                continue
            d = ok(a, b)
            if not d:
                continue
            same = groups.get(b) == groups[a]
            best = (near_in if same else near_out).get(a)
            if best is None or d[0] < best[0]:
                (near_in if same else near_out)[a] = (d[0], b)

    by_group = {}
    for s, (d, _) in near_in.items():
        by_group.setdefault(groups[s], []).append(d)
    tight = {g for g, ds in by_group.items() if len(ds) >= 3 and sorted(ds)[len(ds) // 2] <= threshold}

    out = {}
    for s, (d_in, t_in) in near_in.items():
        if groups[s] not in tight or d_in <= threshold:
            continue
        d_out, t_out = near_out.get(s, (None, None))
        out[s] = {"g": groups[s], "din": d_in, "nin": t_in, "dout": d_out, "nout": t_out,
                  "gout": groups.get(t_out, "") if t_out else ""}
    return out
