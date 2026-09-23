"""File parsers for the QC report: pipeline artefacts in, plain Python values out.

Every parser is defensive by design. The optional inputs (BED, GFF, gene burden,
pN/pS, drug resistance, gene conversion, Kraken, mapDamage, samplesheet, VCF) return an empty
value instead of raising, so a missing or malformed optional file hides its panel
rather than failing the run. The two REQUIRED inputs, parse_summary and
consensus_stats, deliberately raise: the pipeline always hands them a file it has
just produced, so a problem there is a bug and must be loud.

The low-level `_dyn_*` helpers live here rather than in panels.py because they are
parsers too (of a VCF FORMAT/INFO field), and the `parse_*` functions above them
are their only callers.
"""
from __future__ import annotations

import gzip
import os
import re
import sys

NBINS = 200  # bins along the reference for the genome callability-landscape heatmap


def parse_summary(path):
    rows = {}
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < len(header):
                parts += [""] * (len(header) - len(parts))
            d = dict(zip(header, parts))
            sid = (d.get("sample_id") or parts[0]).strip()
            if not sid:
                continue
            rows[sid] = d
    return rows


def consensus_stats(path):
    op = gzip.open if str(path).endswith(".gz") else open
    name, chunks = None, []
    with op(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is None:            # keep the first record's id, but keep reading ALL records so a
                    name = line[1:].split()[0]  # multi-contig consensus is measured whole-genome, not contig 1 only
            else:
                chunks.append(line.strip())
    seq = "".join(chunks).upper()
    L = len(seq) or 1
    gap, n = seq.count("-"), seq.count("N")
    iup = sum(seq.count(c) for c in "RYSWKMBDHV")
    acgt = sum(seq.count(c) for c in "ACGT")
    # gap-length view: one big deletion vs many scattered holes are very different for a phylogeny
    runs = re.findall(r"[-N]+", seq)
    longest_run = max((len(r) for r in runs), default=0)
    # binned missingness profile ALONG the reference (consensus is reference-aligned) -> the genome landscape heatmap.
    # A callability proxy for coverage: per bin, the fraction of positions that are N / gap (no confident base).
    prof = []
    for b in range(NBINS):
        s0, e0 = (b * L) // NBINS, ((b + 1) * L) // NBINS
        sub = seq[s0:e0]
        ln = len(sub) or 1
        prof.append(round((sub.count("-") + sub.count("N")) / ln, 3))
    return {"consensus_id": name, "length": len(seq), "missing_pct": 100.0 * (gap + n) / L,
            "iupac_pct": 100.0 * iup / L, "callable_pct": 100.0 * acgt / L,
            "longest_n_run": longest_run, "n_gaps": len(runs), "miss_profile": prof}


def to_float(x):
    try:
        if x in (None, "", "NA", "."):
            return None
        return float(x)
    except (ValueError, TypeError):
        return None


def clean_str(x):
    """A metadata string, or None for the NA/empty sentinels (so 'NA' lineage does not read as a real lineage)."""
    if x is None:
        return None
    s = str(x).strip()
    return s if s and s not in ("NA", ".", "-") else None


_PROFILE_WARNED = set()


def parse_profile(s):
    """A comma-joined per-bin variant-count profile string -> [int]*NBINS, or None if absent/unparseable.
    Coerced (pad with 0 / truncate) to NBINS: the genome-landscape JS assumes every profile has exactly
    NBINS bins so that a brushed span maps 1:1 onto the hotspot table. Upstream (stats_to_legacy.BIN_N)
    already emits NBINS, so a mismatch means a foreign / re-binned summary; warn once per observed length."""
    if not s or str(s) in ("NA", ".", ""):
        return None
    try:
        v = [int(float(x)) for x in str(s).split(",")]
    except (ValueError, TypeError):
        return None
    if not v:
        return None
    if len(v) != NBINS:
        if len(v) not in _PROFILE_WARNED:
            sys.stderr.write(
                "[qc_report] warning: variant profile has %d bins, expected %d; "
                "padding/truncating so the genome landscape stays aligned.\n" % (len(v), NBINS))
            _PROFILE_WARNED.add(len(v))
        v = (v + [0] * NBINS)[:NBINS]
    return v


def parse_lineage_colors(path):
    """Read a 'lineage<TAB>#hex<TAB>desc' palette (e.g. mycolorsTB) -> {lineage: '#hex'}. Empty if absent."""
    if not path or not os.path.exists(path):
        return {}
    out = {}
    try:
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                c = line.rstrip("\n").split("\t")
                if len(c) >= 2 and c[1].strip().startswith("#") and len(c[1].strip()) in (4, 7):
                    out[c[0].strip()] = c[1].strip()
    except OSError:
        return {}
    return out


def parse_bed(path):
    """Parse a BED (0-based half-open) of regions to mask -> [(start, end)]. Empty if absent."""
    if not path or not os.path.exists(path):
        return []
    iv = []
    try:
        op = gzip.open if str(path).endswith(".gz") else open
        with op(path, "rt", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if not line.strip() or line.startswith(("#", "track", "browser")):
                    continue
                c = line.split("\t")
                if len(c) >= 3:
                    try:
                        iv.append((int(c[1]), int(c[2])))
                    except ValueError:
                        pass
    except OSError:
        return []
    return iv


def mask_profile(intervals, genome_len, nbins=NBINS):
    """Per-bin masked fraction (0-1) along the reference + the total masked %."""
    if not intervals or genome_len <= 0:
        return None, 0.0
    binbp = genome_len / nbins
    cov = [0.0] * nbins
    kept = []
    for s, e in intervals:
        s = max(0, s); e = min(genome_len, e)
        if e <= s:
            continue
        kept.append((s, e))
        b0 = int(s // binbp); b1 = min(nbins - 1, int((e - 1) // binbp))
        for b in range(b0, b1 + 1):
            bs = b * binbp; be = (b + 1) * binbp
            cov[b] += max(0.0, min(e, be) - max(s, bs))
    # the total is the UNION of the intervals: an un-merged BED with overlapping records must not
    # count the same base twice (that reported >100% masked), so sweep them in coordinate order
    total = 0
    reach = 0   # rightmost end already counted
    for s, e in sorted(kept):
        if e > reach:
            total += e - max(s, reach)
            reach = e
    prof = [round(min(1.0, cov[b] / binbp), 3) if binbp > 0 else 0.0 for b in range(nbins)]
    return prof, round(100.0 * total / genome_len, 2)


def parse_gene_burden(path):
    """Cohort gene-burden aux TSV (from collect_summary --gene-burden-out) -> list of gene dicts sorted by
    total_impactful desc. Optional file; empty on any problem. Bounded to the top 200 genes."""
    if not path or not os.path.exists(path):
        return []
    out = []
    try:
        with open(path, encoding="utf-8", errors="replace") as fh:
            header = fh.readline().rstrip("\n").split("\t")
            idx = {h: i for i, h in enumerate(header)}

            def sget(c, k):
                i = idx.get(k)
                return c[i] if (i is not None and i < len(c) and c[i] != "") else None

            def nget(c, k):
                v = sget(c, k)
                try:
                    return int(float(v))
                except (TypeError, ValueError):
                    return None
            for line in fh:
                if not line.strip():
                    continue
                c = line.rstrip("\n").split("\t")
                gene = sget(c, "gene")
                if not gene:
                    continue
                hi, mod = nget(c, "high") or 0, nget(c, "moderate") or 0
                tot = nget(c, "total_impactful")
                out.append({"gene": gene, "high": hi, "moderate": mod,
                            "low": nget(c, "low"), "modifier": nget(c, "modifier"),
                            "dominant_effect": sget(c, "dominant_effect") or sget(c, "dominant") or "NA",
                            "n_samples": nget(c, "n_samples"),
                            "total_impactful": tot if tot is not None else hi + mod,
                            "start": nget(c, "start"), "end": nget(c, "end")})
    except Exception:
        return []
    out.sort(key=lambda r: (-(r.get("total_impactful") or 0), r["gene"]))
    return out[:200]


def parse_pnps(path):
    """Optional cohort dN/dS (eskaks) per-gene TSV -> [{gene, pn, ps, pnps, n, effect}]. Empty on any problem
    (panel hidden). Real eskaks-aggregate columns: gene, gene_name, n_samples, n_pairs, mean_dN, mean_dS, dNdS,
    effect_dominant. inf/nan/empty are read as None, never blindly float()-ed."""
    if not path or not os.path.exists(path):
        return []

    def num(v):
        if v is None:
            return None
        v = str(v).strip().lower()
        if v in ("", "na", "nan", "inf", "-inf", "+inf"):
            return None
        return to_float(v)
    rows = []
    try:
        with open(path, encoding="utf-8", errors="replace") as fh:
            header = fh.readline().rstrip("\n").split("\t")
            idx = {h.strip().lower(): i for i, h in enumerate(header)}

            def col(c, p):
                i = idx.get(c)
                return p[i] if (i is not None and i < len(p)) else None
            for line in fh:
                if not line.strip():
                    continue
                p = line.rstrip("\n").split("\t")
                g = col("gene_name", p) or col("gene", p)
                if not g:
                    continue
                rows.append({"gene": g, "pn": num(col("mean_dn", p)), "ps": num(col("mean_ds", p)),
                             "pnps": num(col("dnds", p)), "n": num(col("n_pairs", p)),
                             "effect": col("effect_dominant", p)})
    except OSError:
        return []
    return rows[:300]


# Roughly clinical order: first line, then the injectables and the rest. AMI is the catalogue's
# spelling of amikacin; AMK is kept beside it because other panels use that one. No composite label
# appears here, because `_dr_components` splits every composite into the drugs it names before the
# matrix is built.
_DR_DRUG_ORDER = ["RIF", "INH", "EMB", "PZA", "STR", "STM", "FQ", "LFX", "MFX", "OFX", "AMI", "AMK",
                  "KAN", "CAP", "ETH", "PTO", "LZD", "BDQ", "CFZ", "DLM", "PAS", "CS", "OTHER"]

# The single-drug codes the WHO catalogue uses, which is what makes a composite recognisable.
_DR_KNOWN = {"RIF", "INH", "EMB", "PZA", "STR", "STM", "FQ", "LFX", "MFX", "OFX", "AMI", "AMK",
             "KAN", "CAP", "ETH", "PTO", "LZD", "BDQ", "CFZ", "DLM", "PAS", "CS"}


def _dr_components(label):
    """The drugs a label names. `AMI_KAN_CAP` is three drugs, `RIF` is one.

    The WHO catalogue grades every variant-drug pair separately, so one variant can be an
    established marker for several drugs at once. Catalogue v1.0.2 writes that as an
    underscore-joined label, and a variant graded 1 for amikacin and kanamycin arrives as
    `AMI_KAN`. Left whole it would get a column of its own and contribute nothing to either drug,
    so a sample whose only isoniazid evidence is an inhA promoter variant, now labelled `INH_ETH`,
    would read as clean in the INH column. The calls table and the TSV keep the label as written;
    only the matrix expands it.

    Split only when EVERY part is a drug code we recognise. A label from somebody else's panel
    that happens to contain an underscore is left alone rather than chopped into nonsense.
    """
    if not label:
        return []
    if label in _DR_KNOWN or "_" not in label:
        return [label]
    parts = label.split("_")
    return parts if all(p in _DR_KNOWN for p in parts) else [label]


def parse_dr(path):
    """Per-sample drug-resistance calls -> {'samples','drugs','calls':[{s,drug,gene,mutation,grade,gn,marker,
    af,dp}]} or None. Tab-separated with a header; recognised columns (case-insensitive): sample, drug, gene,
    mutation, grade, marker_name/marker, af, dp. `grade` may be a WHO string like '1) Assoc w R' -> gn=1."""
    if not path or not os.path.exists(path):
        return None
    calls, samples = [], []
    try:
        with open(path, encoding="utf-8", errors="replace") as fh:
            header = None
            for line in fh:
                if not line.strip() or line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if header is None:
                    header = [h.strip().lower() for h in parts]
                    continue
                d = dict(zip(header, parts))
                s = (d.get("sample") or (parts[0] if parts else "")).strip()
                if not s:
                    continue
                grade = (d.get("grade") or "").strip()
                mg = re.match(r"\s*(\d)", grade)
                dp = d.get("dp")
                try:   # via float(): a depth written as '40.0' is still an integer depth, not a parse failure
                    dp = int(float(dp)) if dp not in (None, "", "NA", ".") else None
                except ValueError:
                    dp = None
                drug = (d.get("drug") or "").strip()
                calls.append({"s": s, "drug": drug, "dr": _dr_components(drug),
                              "gene": (d.get("gene") or "").strip(),
                              "mutation": (d.get("mutation") or "").strip(), "grade": grade,
                              "gn": (int(mg.group(1)) if mg else None),
                              "marker": (d.get("marker_name") or d.get("marker") or "").strip(),
                              "af": to_float(d.get("af")), "dp": dp})
                if s not in samples:
                    samples.append(s)
    except OSError:
        return None
    if not calls:
        return None
    drugs = sorted({d for c in calls for d in c["dr"] if d},
                   key=lambda x: (_DR_DRUG_ORDER.index(x) if x in _DR_DRUG_ORDER else 99, x))
    return {"samples": samples, "drugs": drugs, "calls": calls}


_GCONV_VERDICTS = ("gene_conversion", "mismapping", "ambiguous", "coverage_shift",
                   "reference_artifact", "reference_derived", "reciprocal_exchange")


def parse_gene_conversion(path):
    """Cohort gene-conversion tracts (COLLECT_GENE_CONVERSION) -> {'samples','counts','tracts':[...]} or None.

    One row per candidate tract. Two kinds of number are carried through, and they answer
    different questions.

    `log10_bf` and `post_conv` are what the model concluded: a conversion weighed against an
    independent substitution and against reads that arrived from the donor. `mismap_frac` is the
    donor-read rate it had to assume to say so.

    `donor_af_in`, `donor_af_outside` and `breakpoint_reads` are what a person checks that
    against in the BAM: how fixed the donor allele is inside the tract, whether it also turns up
    outside it, and how many single molecules carry both in cis. The panel shows the judgement
    and the evidence side by side, because a conclusion nobody can check is not much use.

    Optional input: missing, empty or header-only returns None and the panel self-hides.
    """
    if not path or not os.path.exists(path):
        return None

    def ival(v):
        try:   # via float(): a count written as '3.0' is still 3, not a parse failure
            return int(float(v)) if v not in (None, "", "NA", ".") else None
        except (TypeError, ValueError):
            return None
    tracts, samples = [], []
    try:
        with open(path, encoding="utf-8", errors="replace") as fh:
            header = None
            for line in fh:
                if not line.strip() or line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if header is None:
                    header = [h.strip().lower() for h in parts]
                    continue
                d = dict(zip(header, parts))
                s = (d.get("sample") or "").strip()
                if not s:
                    continue
                tracts.append({"s": s, "pair": (d.get("pair_id") or "").strip(),
                               "contig": (d.get("contig") or "").strip(),
                               "donor": (d.get("donor") or "").strip(),
                               "verdict": ((d.get("cohort_verdict") or d.get("verdict") or "")
                                           .strip().lower()),
                               "sample_verdict": (d.get("verdict") or "").strip().lower(),
                               "reason": (d.get("reason") or "").strip(),
                               "start": ival(d.get("start")), "end": ival(d.get("end")),
                               "span": ival(d.get("span_bp")), "n_sites": ival(d.get("n_sites")),
                               "n_out": ival(d.get("n_sites_outside")),
                               "n_undet": ival(d.get("n_undetermined")),
                               "bf": to_float(d.get("log10_bf")),
                               "bf_null": to_float(d.get("log10_bf_vs_null")),
                               "post": to_float(d.get("post_conv")),
                               "mismap": to_float(d.get("mismap_frac")),
                               "tract_af": to_float(d.get("tract_af")),
                               "mut_rate": to_float(d.get("mut_rate")),
                               "event": (d.get("event_id") or "").strip(),
                               "n_ev": ival(d.get("event_samples")),
                               "ev_frac": to_float(d.get("event_frac")),
                               "co_mismap": to_float(d.get("cohort_mismap")),
                               "co_bf": to_float(d.get("cohort_bf_median")),
                               "don_start": ival(d.get("don_start")),
                               "don_end": ival(d.get("don_end")),
                               "n_don": ival(d.get("n_donors")),
                               "don_call": (d.get("donor_call") or "").strip(),
                               "rep": ival(d.get("is_representative")),
                               "don_rank": ival(d.get("donor_rank")),
                               "don_margin": to_float(d.get("donor_margin")),
                               "start_ci": (d.get("start_ci") or "").strip(),
                               "end_ci": (d.get("end_ci") or "").strip(),
                               "af_in": to_float(d.get("donor_af_in")),
                               "af_out": to_float(d.get("donor_af_outside")),
                               "depth": ival(d.get("min_depth")),
                               "cis_reads": ival(d.get("cis_reads")),
                               "bp_reads": ival(d.get("breakpoint_reads")),
                               "donor_only": ival(d.get("donor_only_reads")),
                               "n_derived": ival(d.get("n_derived")),
                               "n_ancestral": ival(d.get("n_ancestral")),
                               "n_unpol": ival(d.get("n_unpolarised")),
                               "donor_swap": to_float(d.get("donor_swap_af")),
                               "donor_swap_out": to_float(d.get("donor_swap_af_outside")),
                               "locus_cn": to_float(d.get("locus_cn")),
                               "expected_af": to_float(d.get("expected_af")),
                               "genes": (d.get("genes") or "").strip(),
                               "n_syn": ival(d.get("n_syn")),
                               "n_nonsyn": ival(d.get("n_nonsyn")),
                               "aa_changes": (d.get("aa_changes") or "").strip()})
                if s not in samples:
                    samples.append(s)
    except OSError:
        return None
    if not tracts:
        return None
    counts = dict.fromkeys(_GCONV_VERDICTS, 0)
    for t in tracts:
        counts[t["verdict"]] = counts.get(t["verdict"], 0) + 1
    return {"samples": samples, "counts": counts, "tracts": tracts}


# A child clade holding at least this share of its parent's reads is where the reads really sit, so
# the dominant-clade walk steps into it; below it the reads are spread over several children (or
# assigned to the parent itself) and the walk stops.
KRAKEN_DOMINANT = 0.5
KRAKEN_TOP = 6


def _kraken_read(path, nodes, children):
    """Add one Kraken2 report to a sample's tree -> the unclassified reads it carried, or None when
    the file cannot be read. `nodes` is {taxid: {name, rank, reads, pct}} and `children` is
    {parent taxid (None for the top level): [taxid]}; both are shared by the reports of one sample,
    so the reads of several runs add up on the same taxa.

    The report records the tree only through the indentation of its name column (two spaces per
    level), so a node's parent is the nearest line above it that sits one level shallower."""
    unclass, stack = 0, []
    try:
        with open(path, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                c = line.rstrip("\n").split("\t")
                if len(c) < 6:
                    continue
                try:
                    pct, reads = float(c[0]), int(float(c[1]))
                except ValueError:
                    continue
                rank, raw = c[3].strip(), c[5]
                if rank == "U":
                    unclass += reads
                    continue
                name = raw.strip()
                depth = (len(raw) - len(raw.lstrip(" "))) // 2
                key = c[4].strip() or name
                while stack and stack[-1][0] >= depth:
                    stack.pop()
                parent = stack[-1][1] if stack else None
                stack.append((depth, key))
                if key not in nodes:
                    nodes[key] = {"name": name, "rank": rank, "reads": 0, "pct": 0.0}
                    children.setdefault(parent, []).append(key)
                nodes[key]["reads"] += reads
                nodes[key]["pct"] += pct
    except OSError:
        return None
    return unclass


def _kraken_walk(start, nodes, children, parent_reads):
    """Follow the most abundant child down from `start` while it holds at least KRAKEN_DOMINANT of its
    parent's reads -> the taxon the reads actually sit in.

    For an M. tuberculosis culture that is the complex, not the species: Kraken2 leaves most reads on
    the complex node (the members are too alike for k-mers to split), so the species clade holds a
    few percent and a species-level reading would call every clean sample contaminated. The walk
    always steps out of the top level and through the root ranks ('root', 'cellular organisms'),
    which say nothing about what is in the tube, and never goes below a species: Kraken2's strain
    assignments are too unreliable to name what a sample is."""
    key = start
    while children.get(key):
        if key is not None and nodes[key]["rank"].startswith("S"):
            break
        best = max(children[key], key=lambda k: (nodes[k]["reads"], nodes[k]["pct"]))
        here = parent_reads if key is None else nodes[key]["reads"]
        forced = key is None or nodes[key]["rank"].startswith("R")
        if not forced and nodes[best]["reads"] < KRAKEN_DOMINANT * here:
            break
        key = best
    return key


def parse_kraken(paths):
    """Per-sample Kraken2 `.report` files (6 columns: clade%, clade reads, taxon reads, rank, taxid,
    name) -> {'samples': [...], 'target': {...}} or None.

    The sample id is the filename up to the first '__' (the pipeline writes one
    '<sampleId>__<runId>.kraken.report' per run), and every report of a sample is merged by summing
    its reads, so a sample sequenced over five runs is one row weighted by what each run produced,
    not five rows. Per sample:
      primary     the dominant clade (see _kraken_walk) and its share of all reads
      secondary   the most abundant taxon outside the primary's lineage, e.g. a contaminant or host
      top         the primary and the next taxa outside its lineage, for the composition bar
      target_pct  share of the CLASSIFIED reads that fall in the cohort's target clade
    The target is the clade most samples have as their primary, so a culture of the wrong organism
    reads as ~0% target rather than as a clean sample of something else."""
    trees, order = {}, []
    for p in paths or []:
        if not p or not os.path.exists(p):
            continue
        sid = os.path.basename(p).split("__")[0].split(".")[0]
        tree = trees.get(sid) or {"nodes": {}, "children": {}, "unclass": 0, "runs": 0}
        unclass = _kraken_read(p, tree["nodes"], tree["children"])
        if unclass is None:
            continue
        tree["unclass"] += unclass
        tree["runs"] += 1
        if sid not in trees:
            trees[sid] = tree
            order.append(sid)

    out = []
    for sid in order:
        t = trees[sid]
        nodes, children = t["nodes"], t["children"]
        classified = sum(nodes[k]["reads"] for k in children.get(None, []))
        total = classified + t["unclass"]
        below_root = any(children.get(k) for k in children.get(None, [])) or \
            any(not nodes[k]["rank"].startswith("R") for k in children.get(None, []))
        if not total or (not below_root and not t["unclass"]):
            continue

        def share(k, _total=total):
            return round(100.0 * nodes[k]["reads"] / _total, 2) if _total else 0.0

        def taxon(k, _share=share):
            return {"name": nodes[k]["name"], "rank": nodes[k]["rank"], "pct": _share(k)}

        prim = _kraken_walk(None, nodes, children, classified) if below_root else None
        # Everything off the primary's lineage: at each ancestor, the sibling clades not on the path.
        path, k = [], prim
        parent_of = {c: par for par, cs in children.items() for c in cs}
        while k is not None:
            path.append(k)
            k = parent_of.get(k)
        on_path = set(path)
        others = [c for anc in [None] + path[1:] for c in children.get(anc, []) if c not in on_path]
        others.sort(key=lambda c: (nodes[c]["reads"], nodes[c]["pct"]), reverse=True)
        described = []
        for c in others:
            if not nodes[c]["reads"] and not nodes[c]["pct"]:
                continue
            d = _kraken_walk(c, nodes, children, nodes[c]["reads"])
            described.append(taxon(d))
            if len(described) >= KRAKEN_TOP - 1:
                break
        top = ([taxon(prim)] if prim is not None else []) + described
        out.append({"s": sid, "runs": t["runs"], "reads": total, "taxid": prim,
                    "unclassified": round(100.0 * t["unclass"] / total, 2),
                    "classified": round(100.0 * classified / total, 2),
                    "primary": taxon(prim) if prim is not None else None,
                    "secondary": described[0] if described else None,
                    "top": top, "_tree": t, "_classified": classified})
    if not out:
        return None

    # The target: the clade most samples are dominated by (ties go to the most reads).
    votes = {}
    for s in out:
        if s["taxid"] is not None:
            v = votes.setdefault(s["taxid"], [0, 0, s["primary"]])
            v[0] += 1
            v[1] += s["_tree"]["nodes"][s["taxid"]]["reads"]
    target = max(votes.items(), key=lambda kv: (kv[1][0], kv[1][1]))[0] if votes else None
    for s in out:
        nodes, classified = s["_tree"]["nodes"], s.pop("_classified")
        s.pop("_tree")
        s.pop("taxid")
        hit = nodes.get(target)
        s["target_pct"] = (round(100.0 * hit["reads"] / classified, 2) if hit and classified else
                           (0.0 if classified else None))
    tgt = None
    if target is not None:
        tgt = {"name": votes[target][2]["name"], "rank": votes[target][2]["rank"], "n": votes[target][0]}
    return {"samples": out, "target": tgt}


def _gff_attr(attrs, key):
    """Value of GFF3 attribute `key` from the column-9 string, or None when absent. The key is matched
    WHOLE against each ';'-separated field (a substring search would read 'locus_tag' out of
    'old_locus_tag' and 'gene' out of 'pseudogene', both legal GFF3 and routine in RefSeq/Prokka).
    Only the first '=' splits, so a value that itself contains '=' survives."""
    for field in attrs.split(";"):
        k, sep, v = field.partition("=")
        if sep and k.strip() == key:
            return v.strip()
    return None


def parse_gff(path):
    """Best-effort GFF3 parse of gene/CDS features -> [{name, start, end}] (1-based). Empty on any problem."""
    if not path or not os.path.exists(path):
        return []
    genes = {}
    try:
        op = gzip.open if str(path).endswith(".gz") else open
        with op(path, "rt", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith("#") or "\t" not in line:
                    continue
                c = line.rstrip("\n").split("\t")
                if len(c) < 9 or c[2] not in ("gene", "CDS"):
                    continue
                try:
                    start, end = int(c[3]), int(c[4])
                except ValueError:
                    continue
                name = None
                for key in ("Name", "gene", "locus_tag", "ID"):
                    name = _gff_attr(c[8], key)
                    if name is not None:
                        break
                name = name or f"{start}-{end}"
                # prefer a 'gene' feature over a 'CDS' of the same name
                if name not in genes or (c[2] == "gene" and genes[name].get("type") != "gene"):
                    genes[name] = {"name": name, "start": start, "end": end, "type": c[2]}
    except (OSError, ValueError):
        return []
    return sorted(({"name": g["name"], "start": g["start"], "end": g["end"]} for g in genes.values()),
                  key=lambda g: g["start"])


def parse_gene_locus(path):
    """{gene_display_name: locus_tag} from a GFF3 (e.g. rpoB -> Rv0667). Empty on any problem."""
    out = {}
    if not path or not os.path.exists(path):
        return out
    try:
        op = gzip.open if str(path).endswith(".gz") else open
        with op(path, "rt", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith("#") or "\t" not in line:
                    continue
                c = line.rstrip("\n").split("\t")
                if len(c) < 9 or c[2] not in ("gene", "CDS"):
                    continue
                attrs = c[8]
                locus = _gff_attr(attrs, "locus_tag")
                name = _gff_attr(attrs, "gene") or _gff_attr(attrs, "Name")
                if name and locus and name not in out:
                    out[name] = locus
    except (OSError, ValueError):
        return out
    return out


def _read_freq_file(path, max_pos=25):
    """Parse a mapDamage 5pCtoT_freq.txt / 3pGtoA_freq.txt ('pos<TAB>freq'). Returns
    {'pos1': float|None, 'profile': [[pos,freq],...] up to max_pos} or None. Robust to header + blanks."""
    prof = []
    try:
        with open(path) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.lower().startswith("pos"):
                    continue
                parts = line.replace(",", "\t").split("\t")
                if len(parts) < 2:
                    parts = line.split()
                if len(parts) < 2:
                    continue
                p = to_float(parts[0]); v = to_float(parts[1])
                if p is None or v is None:
                    continue
                prof.append((int(p), v))
    except (OSError, ValueError):
        return None
    if not prof:
        return None
    prof.sort(key=lambda x: x[0])
    pos1 = next((v for p, v in prof if p == 1), prof[0][1])
    return {"pos1": pos1, "profile": [[p, round(v, 5)] for p, v in prof[:max_pos]]}


def _mean_fraglen(path):
    """Mean fragment length from lgdistribution.txt (Std<TAB>Length<TAB>Occurences, '#'-comment lines)."""
    num = den = 0.0
    try:
        with open(path) as fh:
            for line in fh:
                if not line.strip() or line.startswith("#"):
                    continue
                parts = line.split()
                if len(parts) < 3 or parts[1].lower() == "length":
                    continue
                L = to_float(parts[1]); occ = to_float(parts[2])
                if L is None or occ is None:
                    continue
                num += L * occ; den += occ
    except (OSError, ValueError):
        return None
    return (num / den) if den > 0 else None


def mapdamage_stats(sample_dir):
    """Read a mapDamage2 output directory for one sample. Returns a dict or None if nothing parseable.
    keys: ct1 (5' C>T pos1), ga1 (3' G>A pos1), fraglen (mean), ct_profile [[pos,freq],...]."""
    if not sample_dir or not os.path.isdir(sample_dir):
        return None
    five = _read_freq_file(os.path.join(sample_dir, "5pCtoT_freq.txt"))
    three = _read_freq_file(os.path.join(sample_dir, "3pGtoA_freq.txt"))
    fl = _mean_fraglen(os.path.join(sample_dir, "lgdistribution.txt"))
    ct1 = five["pos1"] if five else None
    if five is None and three is None and fl is None and ct1 is None:
        return None
    return {"ct1": ct1, "ga1": (three["pos1"] if three else None), "fraglen": fl,
            "ct_profile": five["profile"] if five else None}


# ---------------------------------------------------------------- samplesheet + VCF (dynamics inputs)
_DYN_SAMPLE_RE = re.compile(r'^(sample_?id|sampleid|sample|name|gid|strain|isolate)$', re.I)
_DYN_TIME_RE   = re.compile(r'(passage|pase|timepoint|time_?point|^time$|^day$|date|week|month|hour|generation|^tp$|visit|^t\d*$)', re.I)
_DYN_GROUP_RE  = re.compile(r'(group|series|patient|host|subject|cluster|experiment|^line$|replicate|chain|pair|lineage_?id|donor|case|animal|^samples$)', re.I)
_DOSE_RE       = re.compile(r'^(dose|dosis)$', re.I)                           # a numeric annotation -> a metric, not a categorical level
_TX_RE         = re.compile(r'(treatment|tratamiento|regimen|therapy|^arm$)', re.I)   # the categorical column the dose x treatment test groups by


def _dyn_open(path):
    return gzip.open(path, 'rt', encoding='utf-8', errors='replace') if str(path).endswith('.gz') \
        else open(path, 'r', encoding='utf-8', errors='replace')


def _dyn_num(x):
    if x is None:
        return None
    m = re.search(r'-?\d+\.?\d*', str(x))
    return float(m.group()) if m else None


_META_SKIP = {'r1', 'r2', 'reffasta', 'refgff', 'refid', 'taxid', 'runid'}


def parse_sample_meta(path):
    """{'fields':[col...], 'rows':{sample:{col:value}}} for the user's annotation columns of the
    samplesheet (everything except the sample id and the pipeline's file/reference columns), in
    samplesheet order -> one column-header level per field. None if nothing usable."""
    if not path or not os.path.exists(path):
        return None
    header, rows = None, []
    try:
        with _dyn_open(path) as fh:
            for line in fh:
                if not line.strip() or line.startswith('#'):
                    continue
                cells = line.rstrip('\n').split('\t')
                if header is None:
                    header = [c.strip() for c in cells]
                    continue
                rows.append(cells)
    except OSError:
        return None
    if not header:
        return None
    si = next((i for i, h in enumerate(header) if _DYN_SAMPLE_RE.match(h.strip())), 0)
    # the dose column (if any) is a numeric metric, not a categorical annotation -> keep it out of fields
    fields = [(i, h) for i, h in enumerate(header)
              if i != si and h.strip() and h.strip().lower() not in _META_SKIP
              and not _DOSE_RE.match(h.strip())]
    if not fields:
        return None
    out = {}
    for r in rows:
        if si >= len(r):
            continue
        s = r[si].strip()
        if not s or s in out:
            continue
        out[s] = {h: (r[i].strip() if i < len(r) else '') for i, h in fields}
    if not out:
        return None
    # tag the time / group columns (dynamics axes) and the treatment column (dose x treatment test)
    # so the report can route each correctly; they all stay as SNP-matrix header levels either way.
    time_field = next((h for _, h in fields if _DYN_TIME_RE.search(h.replace(' ', '_'))), None)
    group_field = next((h for _, h in fields if _DYN_GROUP_RE.search(h.replace(' ', '_'))), None)
    tx_field = next((h for _, h in fields if _TX_RE.search(h.replace(' ', '_'))), None)
    return {'fields': [h for _, h in fields], 'rows': out,
            'time_field': time_field, 'group_field': group_field, 'tx_field': tx_field}


_DATE_RE = re.compile(r'^(collection|sampling|isolation|sample)?[_ ]?(date|year)$|^fecha$', re.I)


def parse_collection_dates(path):
    """{sample: date string} from the samplesheet column that says when a sample was COLLECTED
    ('collection_date', 'sampling_date', 'date', 'year'...). {} when there is no such column.

    This is the only source of a sampling date. The summary TSV carries a 'date' too, but it is
    DAT_OUT, the day the pipeline processed the sample: read as a collection date it put a whole
    cohort in one year and the temporal panel reported a sampling span of zero years."""
    if not path or not os.path.exists(path):
        return {}
    header, out = None, {}
    try:
        with _dyn_open(path) as fh:
            for line in fh:
                if not line.strip() or line.startswith('#'):
                    continue
                cells = line.rstrip('\n').split('\t')
                if header is None:
                    header = [c.strip() for c in cells]
                    si = next((i for i, h in enumerate(header) if _DYN_SAMPLE_RE.match(h)), 0)
                    di = next((i for i, h in enumerate(header) if _DATE_RE.match(h)), None)
                    if di is None:
                        return {}
                    continue
                if max(si, di) >= len(cells):
                    continue
                s, d = cells[si].strip(), clean_str(cells[di])
                if s and d and s not in out:
                    out[s] = d
    except OSError:
        return {}
    return out


def parse_dose(path):
    """{sample: float} for a numeric samplesheet column named 'dose'/'dosis' (a quantitative
    annotation exposed as a report metric). {} if the file or column is absent / non-numeric."""
    if not path or not os.path.exists(path):
        return {}
    header, rows = None, []
    try:
        with _dyn_open(path) as fh:
            for line in fh:
                if not line.strip() or line.startswith('#'):
                    continue
                cells = line.rstrip('\n').split('\t')
                if header is None:
                    header = [c.strip() for c in cells]
                    continue
                rows.append(cells)
    except OSError:
        return {}
    if not header:
        return {}
    si = next((i for i, h in enumerate(header) if _DYN_SAMPLE_RE.match(h.strip())), 0)
    di = next((i for i, h in enumerate(header) if _DOSE_RE.match(h.strip())), None)
    if di is None:
        return {}
    out = {}
    for r in rows:
        if si >= len(r) or di >= len(r):
            continue
        s = r[si].strip()
        if not s or s in out:
            continue
        v = to_float(r[di].strip())
        if v is not None:
            out[s] = v
    return out


def parse_metadata(path):
    """{sample: {'time','tnum','group'}} auto-detecting sample/time/group columns; {} if unusable."""
    if not path or not os.path.exists(path):
        return {}
    header, rows = None, []
    try:
        with _dyn_open(path) as fh:
            for line in fh:
                if not line.strip() or line.startswith('#'):
                    continue
                cells = line.rstrip('\n').split('\t')
                if header is None:
                    header = [c.strip() for c in cells]
                    continue
                rows.append(cells)
    except OSError:
        return {}
    if not header:
        return {}
    def find(rx):
        for i, h in enumerate(header):
            if rx.search(h.replace(' ', '_')):
                return i
        return None
    si = next((i for i, h in enumerate(header) if _DYN_SAMPLE_RE.match(h.strip())), 0)
    ti, gi = find(_DYN_TIME_RE), find(_DYN_GROUP_RE)
    out = {}
    for r in rows:
        if si >= len(r):
            continue
        s = r[si].strip()
        if not s or s in out:
            continue
        out[s] = {'time':  r[ti].strip() if (ti is not None and ti < len(r)) else None,
                  'tnum':  _dyn_num(r[ti]) if (ti is not None and ti < len(r)) else None,
                  'group': r[gi].strip() if (gi is not None and gi < len(r)) else None}
    return out


def _dyn_af(fmt, val):
    """First-alt allele fraction: AD/AO+RO (true AF) else GT dosage (hom 1.0 / het 0.5)."""
    d = dict(zip(fmt.split(':'), val.split(':')))
    if 'AD' in d:
        try:
            ad = [int(x) for x in d['AD'].split(',') if x not in ('.', '')]
            if len(ad) >= 2 and sum(ad) > 0:
                return sum(ad[1:]) / sum(ad)
        except ValueError:
            pass
    if 'AO' in d and 'RO' in d:
        try:
            ao = sum(int(x) for x in d['AO'].split(',') if x not in ('.', ''))
            ro = int(d['RO']) if d['RO'] not in ('.', '') else 0
            if ao + ro > 0:
                return ao / (ao + ro)
        except ValueError:
            pass
    alleles = [a for a in d.get('GT', './.').replace('|', '/').split('/') if a not in ('.', '')]
    if not alleles:
        return None
    return sum(1 for a in alleles if a != '0') / len(alleles)


def _dyn_dp(fmt, val):
    """Total depth at the site: FORMAT DP, else sum(AD), else AO+RO."""
    d = dict(zip(fmt.split(':'), val.split(':')))
    if d.get('DP', '.') not in ('.', ''):
        try:
            return int(d['DP'])
        except ValueError:
            pass
    if 'AD' in d:
        try:
            ad = [int(x) for x in d['AD'].split(',') if x not in ('.', '')]
            if ad:
                return sum(ad)
        except ValueError:
            pass
    if 'AO' in d and 'RO' in d:
        try:
            return sum(int(x) for x in d['AO'].split(',') if x not in ('.', '')) + \
                (int(d['RO']) if d['RO'] not in ('.', '') else 0)
        except ValueError:
            pass
    return None


def _dyn_ann(info):
    # snpEff ANN fields: 3=gene, 1=effect, 2=impact, 10=HGVS.p (protein change, e.g. p.Ser315Thr)
    for field in info.split(';'):
        if field.startswith('ANN='):
            a = field[4:].split(',')[0].split('|')
            return ((a[3] if len(a) > 3 else ''), (a[1] if len(a) > 1 else ''),
                    (a[2] if len(a) > 2 else ''), (a[10] if len(a) > 10 else ''))
    return '', '', '', ''


def _dyn_opos(info):
    """The ORIGINAL (used-reference) coordinate a liftover stamped onto a canonical-reference VCF record,
    so a used-ref variant can be paired with its reference-of-interest position WITHOUT assuming the two
    references share coordinates. Reads Picard LiftoverVcf's OriginalContig/OriginalStart, else a plain
    OPOS=<contig:pos>. Returns 'contig:pos' (used-ref coordinates) or None."""
    if not info:
        return None
    oc = ost = None
    for field in info.split(';'):
        if field.startswith('OriginalContig='):
            oc = field[len('OriginalContig='):]
        elif field.startswith('OriginalStart='):
            ost = field[len('OriginalStart='):]
        elif field.startswith('OPOS='):
            return field[len('OPOS='):].strip() or None
    return (oc + ':' + ost) if (oc and ost) else None


def parse_vcfs(paths):
    """{sample: {'chrom:pos': {ref,alt,gene,eff,imp,af}}} from annotated per-sample VCFs (SNPs only).
    Sample = the VCF #CHROM last column, so filenames are irrelevant."""
    out = {}
    for p in paths or []:
        if not p or not os.path.exists(p):
            continue
        try:
            with _dyn_open(p) as fh:
                sample = None
                for line in fh:
                    if line.startswith('##'):
                        continue
                    if line.startswith('#CHROM'):
                        cols = line.rstrip('\n').split('\t')
                        sample = cols[9] if len(cols) > 9 else os.path.basename(p).split('.')[0]
                        out.setdefault(sample, {})
                        continue
                    if sample is None:
                        continue
                    c = line.rstrip('\n').split('\t')
                    if len(c) < 8:
                        continue
                    chrom, pos, ref, alt = c[0], c[1], c[3], c[4]
                    if alt in ('.', '') or len(ref) != 1 or any(len(a) != 1 for a in alt.split(',')):
                        continue
                    af = _dyn_af(c[8], c[9]) if len(c) >= 10 else 1.0
                    if af is None:
                        continue
                    dp = _dyn_dp(c[8], c[9]) if len(c) >= 10 else None
                    gene, eff, imp, aa = _dyn_ann(c[7])
                    out[sample][f'{chrom}:{pos}'] = {'ref': ref, 'alt': alt.split(',')[0], 'gene': gene,
                                                     'eff': eff, 'imp': imp, 'aa': aa, 'aa_h37rv': '',
                                                     'af': round(af, 4), 'dp': dp, 'opos': _dyn_opos(c[7])}
        except OSError:
            continue
    return out
