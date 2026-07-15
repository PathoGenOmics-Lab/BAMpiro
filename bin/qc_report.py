#!/usr/bin/env python3
"""Consolidated interactive QC report for sBAMpiro (organism-agnostic).

Reads the per-sample summary TSV (collect_summary output) + the consensus FASTAs, computes consensus completeness
(missing / N / IUPAC / callable), embeds everything as JSON, and writes a SELF-CONTAINED, data-driven qc_report.html
(one inline vanilla-JS module, no CDN) plus a qc_flags.tsv (PASS/WARN/FAIL). Features:
  - General Statistics table: value-coloured cells, sort, sample filter, column show/hide, CSV export, sticky header
    + frozen sample column.
  - LIVE thresholds: adjust the flag cut-offs in the report and everything re-flags instantly.
  - Distributions: beeswarm / bar / histogram (switchable) per metric, with hover tooltips.
  - Correlations: a scatter with X/Y metric pickers, coloured by verdict.
  - Consensus completeness: stacked callable / IUPAC / missing bar per sample.
  - Click any row, dot or point to highlight a sample everywhere; click a flag chip to filter.
With --gate it exits non-zero if any sample FAILs (at the default thresholds), so a downstream target can block.
Pure standard library (no numpy / matplotlib / CDN). Study-specific signatures live in a separate downstream gate.
"""
from __future__ import annotations

import argparse
import gzip
import html
import json
import os
import re
import sys
from datetime import datetime, timezone

DEF = dict(depth_min=10.0, breadth_min=90.0, missing_max=10.0, dup_max=40.0,
           iupac_max=2.0, mapping_min=80.0, titv_min=0.0, snp_z=3.0,
           het_max_frac=3.0, mixed_min_frac=2.0)  # het_max_frac + mixed_min_frac are PERCENTS

# Ancient (aDNA) sample threshold view: looser cov/breadth/missing/mapping/dup/iupac (a 5x mummy must not FAIL the
# modern gate), plus damage_min_ct = the terminal 5' C>T deamination floor for the DAMAGE_LOW authentication screen.
ANC_DEF = dict(depth_min=3.0, breadth_min=30.0, missing_max=70.0, dup_max=60.0,
               iupac_max=5.0, mapping_min=20.0, damage_min_ct=0.05)
FAIL_FLAGS = {"LOW_DEPTH", "LOW_BREADTH", "HIGH_MISSING", "NO_DATA"}
NBINS = 200  # bins along the reference for the genome callability-landscape heatmap

# metric definitions + typical acceptable ranges (surfaced as in-report help)
DEFS = {
    "raw_reads": ["Read pairs/reads before trimming.", "context-dependent"],
    "trimmed_reads": ["Reads after fastp quality/adapter trimming.", ""],
    "read_length": ["Mean read length after trimming.", ""],
    "duplication_pct": ["% reads marked as duplicates (PCR/optical).", "typically < 30-40%"],
    "mapped_pct": ["% reads mapping to the reference.", "> 90% for a clean library"],
    "properly_paired_pct": ["% read pairs mapping in the expected orientation/insert.", "> 90%"],
    "mean_depth": ["Mean sequencing depth over the genome (sensitive to a few high-copy repeats; see Median depth).", ">= 10-20x for confident calls"],
    "median_depth": ["Median sequencing depth over the genome (robust to high-copy repeats that inflate the mean).", ">= 10-20x"],
    "breadth_pct": ["% of the genome covered at >= 1x. Note: 1x is not enough to CALL a base; see Callable % for usable genome.", ">= 90% at 1x to enter a phylogeny"],
    "breadth5x_pct": ["% of the genome covered at >= 5x.", "high for a confident consensus"],
    "breadth10x_pct": ["% of the genome covered at >= 10x.", ">= 90% for confident SNP calls"],
    "coverage_cv": ["Coverage evenness = coefficient of variation of depth (sd / mean). Lower is more even; spikes flag repeats/contamination.", "< ~0.5-1 for an even library"],
    "avg_base_quality": ["Mean per-base quality of mapped reads.", ">= 30"],
    "unmapped_pct": ["% of reads that did not map to the reference (an off-target / contamination proxy).", "low for a clean, on-target library"],
    "singleton_pct": ["% of reads mapped without their mate (approx: reads mapped minus reads mapped-and-paired).", "low"],
    "mapq0_pct": ["% of reads with MAPQ 0 (ambiguous multi-mapping, e.g. repeats / IS elements).", "low"],
    "error_rate_pct": ["samtools mismatch (error) rate, as a percent; elevated with damage or contamination.", "< ~1%"],
    "insert_size": ["Mean fragment/insert size (bp) of properly paired reads. For aDNA a short distribution is expected.", "library-dependent"],
    "insert_size_sd": ["Standard deviation of insert size; a very wide spread signals chimeras/contamination.", ""],
    "q20_pct": ["% of trimmed bases with quality >= Q20 (fastp).", "> 90%"],
    "q30_pct": ["% of trimmed bases with quality >= Q30 (fastp).", "> 80-90%"],
    "gc_pct": ["GC content of the trimmed reads (fastp); a large deviation from the reference GC flags contamination.", "match the reference GC"],
    "est_genome_cov": ["Estimated fold coverage = trimmed bases / genome length (fastp-based, pre-mapping).", "tracks mean_depth"],
    "total_variants": ["Total variants (SNPs + indels) called vs the reference.", "context-dependent"],
    "snps": ["SNPs called vs the reference. Dominated by true phylogenetic distance, NOT quality; the SNP outlier flag is a cohort distributional signal, not a defect.", "expected to be higher for divergent lineages"],
    "snp_density": ["SNP calls per callable megabase = snps / (callable_pct/100 * genome_len/1e6). Normalises SNP burden for how much of the genome was callable; compare within a lineage. Suppressed for a toy genome (< 100 kb).", "compare within a lineage"],
    "het_variants": ["Heterozygous variant calls. In a clonal/haploid organism these are expected ~0; a nonzero rate is consistent with mixture, contamination, or (for aDNA) post-mortem deamination.", "~0"],
    "homo_indels": ["Homozygous indels called.", ""],
    "ti_tv": ["Transition/transversion ratio (a spectrum sanity check, not a rate). Typical genome-wide MTBC ~1.7-2.0; M. bovis (A4) runs higher (~2.6) from its transition-biased spectrum, and aDNA deamination inflates it. LOW Ti/Tv (< ~1.5) is the QC concern (mapping/base-call noise is transversion-heavy); high is biology or damage, not a defect.", "~1.5-2.1 for human lineages"],
    "missing_pct": ["% of the consensus that is '-' or N (no confident base).", "< 10%"],
    "iupac_pct": ["% of the consensus that is an IUPAC ambiguity (mixed call).", "< 2%"],
    "callable_pct": ["% of the consensus that is a confident A/C/G/T (the honest, depth-aware usable-genome measure).", "> 90%"],
    "longest_n_run": ["Longest single run of consecutive N/gap in the consensus (one big deletion vs many scattered holes).", "context-dependent"],
    "n_gaps": ["Number of separate N/gap runs in the consensus.", "context-dependent"],
    "het_frac": ["Fraction of VARIANT genotypes that are heterozygous (het / (het+hom)). An illustrative within-sample mixture signal, not the textbook heterozygosity rate; the denominator is variant-count, so read it with care. Suppressed below a minimum variant count, and NOT applied to aDNA (deamination inflates it).", "~0 for a clonal isolate"],
    "n_lineages": ["Number of lineages (collapsed to main lineage) above the mixture cut-off AND above a minimum marker count; > 1 is consistent with a mixed infection or cross-contamination.", "1"],
    "damage_ct": ["Terminal 5' C>T deamination frequency (mapDamage2). Elevated in authentic aDNA; consistent with post-mortem damage but not proof of authenticity or endogeneity on its own.", "aDNA typically >= 0.05-0.40"],
    "damage_ga": ["Terminal 3' G>A deamination frequency (mapDamage2).", "aDNA typically >= 0.05-0.40"],
    "fraglen": ["Mean fragment (read) length from mapDamage2. Short (~30-70 bp) is expected for aDNA.", "aDNA ~30-70 bp"],
    "ann_high": ["snpEff HIGH-impact variants (stop_gained, frameshift, start/stop_lost, splice). A cohort-relative excess is a contamination / misassembly / wrong-reference-annotation signal, not a per-sample defect on its own.", "cohort-relative robust-z; high flags a look-here sample"],
    "missense_silent": ["Missense / synonymous variant-count ratio: a within-sample pN/pS PROXY (a spectrum sanity check, NOT dN/dS and NOT a selection test; no site-degeneracy correction, entangled with divergence and with how snpEff bins effects). Read within a lineage. For a real pN/pS use the eskaks panel.", "compare within a lineage; robust-z high tail only"],
    "lof_pct": ["Fraction of protein-coding-effect variants that are loss-of-function (stop_gained + frameshift + start/stop_lost + splice). Elevated LoF can be real pseudogenisation OR a frameshift storm from an indel-calling / assembly / annotation artifact.", "cohort-relative robust-z"],
    "coding_pct": ["Percent of annotated variants falling in a transcript (HIGH+MODERATE+LOW; MODIFIER is the intergenic / intron / UTR bucket). Tracks how gene-dense the reference annotation is; not a quality gate.", "reference-dependent, descriptive"],
    "annotated_pct": ["Percent of this sample's real variants that received an ANN annotation from snpEff. A low value flags a GFF-to-database mismatch: the annotation may be silently garbage.", "cohort-relative low + absolute floor"],
    "snpeff_warn": ["Variants carrying a snpEff ERROR/WARNING/INFO subfield (e.g. WARNING_TRANSCRIPT_INCOMPLETE, ERROR_CHROMOSOME_NOT_FOUND). A value near the total variant count means the reference database was built wrong for this contig; treat every annotation with suspicion.", "low; near-total = broken DB"],
}
# Contamination/mixture cut-offs (het_max_frac, mixed_min_frac) and the aDNA damage floor are ILLUSTRATIVE sensitivity
# defaults, not validated gates; they are WARN-only and live-adjustable. Tune to your depth and organism.

# (key, label, kind, direction). kind: int|float|pct ; direction: hi_good|hi_bad|neu
METRICS = [
    ("raw_reads", "Raw reads", "int", "neu"),
    ("trimmed_reads", "Trimmed", "int", "neu"),
    ("read_length", "Read len", "float", "neu"),
    ("q20_pct", "Q20 %", "pct", "hi_good"),
    ("q30_pct", "Q30 %", "pct", "hi_good"),
    ("gc_pct", "GC %", "pct", "neu"),
    ("duplication_pct", "Dup %", "pct", "hi_bad"),
    ("mapped_pct", "Mapped %", "pct", "hi_good"),
    ("properly_paired_pct", "Proper %", "pct", "hi_good"),
    ("unmapped_pct", "Unmapped %", "pct", "hi_bad"),
    ("singleton_pct", "Singleton %", "pct", "hi_bad"),
    ("mapq0_pct", "MAPQ0 %", "pct", "hi_bad"),
    ("error_rate_pct", "Error %", "pct", "hi_bad"),
    ("insert_size", "Insert bp", "float", "neu"),
    ("insert_size_sd", "Insert sd", "float", "neu"),
    ("mean_depth", "Depth", "float", "hi_good"),
    ("median_depth", "Median depth", "float", "hi_good"),
    ("breadth_pct", "Breadth %", "pct", "hi_good"),
    ("breadth5x_pct", "Breadth 5x %", "pct", "hi_good"),
    ("breadth10x_pct", "Breadth 10x %", "pct", "hi_good"),
    ("coverage_cv", "Cov CV", "float", "hi_bad"),
    ("avg_base_quality", "Base Q", "float", "hi_good"),
    ("est_genome_cov", "Est cov", "float", "hi_good"),
    ("total_variants", "Variants", "int", "neu"),
    ("snps", "SNPs", "int", "neu"),
    ("het_variants", "Het", "int", "hi_bad"),
    ("homo_indels", "Indels", "int", "neu"),
    ("ti_tv", "Ti/Tv", "float", "neu"),
    ("missing_pct", "Missing %", "pct", "hi_bad"),
    ("iupac_pct", "IUPAC %", "pct", "hi_bad"),
    ("callable_pct", "Callable %", "pct", "hi_good"),
    ("longest_n_run", "Longest N-run", "int", "hi_bad"),
    ("n_gaps", "N-gaps", "int", "neu"),
    ("ann_high", "HIGH impact", "int", "hi_bad"),
    ("ann_moderate", "MODERATE", "int", "neu"),
    ("ann_low", "LOW impact", "int", "neu"),
    ("ann_modifier", "MODIFIER", "int", "neu"),
    ("missense_silent", "Missense/silent", "float", "neu"),
    ("lof_pct", "LoF %", "pct", "hi_bad"),
    ("coding_pct", "Coding %", "pct", "neu"),
    ("annotated_pct", "Annotated %", "pct", "hi_good"),
    ("snpeff_warn", "snpEff warn", "int", "hi_bad"),
]
DIST = ["mean_depth", "median_depth", "breadth_pct", "breadth10x_pct", "coverage_cv", "mapped_pct",
        "unmapped_pct", "properly_paired_pct", "duplication_pct", "mapq0_pct", "error_rate_pct",
        "missing_pct", "callable_pct", "iupac_pct", "snps", "ti_tv",
        "ann_high", "missense_silent", "annotated_pct"]

# columns that are identifiers / metadata / already-surfaced -> NOT auto-detected as "extra" metrics
KNOWN_KEYS = {k for k, *_ in METRICS} | {
    "sample_id", "lineage", "lineage_counts", "drug_resistance", "dr_counts", "date", "ancient",
    "consensus_id", "length", "genome_len", "mapped_reads", "homo_variants",
    "snp_profile", "het_profile", "indel_profile",   # per-position variant tracks (huge strings, not metrics)
    "eff_missense", "eff_synonymous", "eff_stop_gained", "eff_stop_lost", "eff_start_lost",
    "eff_frameshift", "eff_inframe_indel", "eff_splice", "eff_intergenic", "eff_regulatory",
    "ann_present", "ann_annotated_frac", "ann_total", "ann_annotated", "coding_count", "lof_count", "ann_db_error",
}
_PCT_HINTS = ("_pct", "_rate", "percent", "fraction", "frac")
_INT_HINTS = ("reads", "count", "_n", "num_", "variants", "snps", "indels", "bases")


def parse_summary(path):
    rows = {}
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < len(header):
                parts += [""] * (len(header) - len(parts))
            d = dict(zip(header, parts))
            rows[d.get("sample_id") or parts[0]] = d
    return rows


def consensus_stats(path):
    op = gzip.open if str(path).endswith(".gz") else open
    name, chunks = None, []
    with op(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    break
                name = line[1:].split()[0]
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
    total = 0
    for s, e in intervals:
        s = max(0, s); e = min(genome_len, e)
        if e <= s:
            continue
        total += e - s
        b0 = int(s // binbp); b1 = min(nbins - 1, int((e - 1) // binbp))
        for b in range(b0, b1 + 1):
            bs = b * binbp; be = (b + 1) * binbp
            cov[b] += max(0.0, min(e, be) - max(s, bs))
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
                for key in ("Name=", "gene=", "locus_tag=", "ID="):
                    i = c[8].find(key)
                    if i >= 0:
                        name = c[8][i + len(key):].split(";")[0].strip()
                        break
                name = name or f"{start}-{end}"
                # prefer a 'gene' feature over a 'CDS' of the same name
                if name not in genes or (c[2] == "gene" and genes[name].get("type") != "gene"):
                    genes[name] = {"name": name, "start": start, "end": end, "type": c[2]}
    except (OSError, ValueError):
        return []
    return sorted(({"name": g["name"], "start": g["start"], "end": g["end"]} for g in genes.values()),
                  key=lambda g: g["start"])


def robust(vals):
    v = sorted(x for x in vals if x is not None)
    if not v:
        return None, None
    n = len(v)
    med = v[n // 2] if n % 2 else (v[n // 2 - 1] + v[n // 2]) / 2
    devs = sorted(abs(x - med) for x in v)
    mad = devs[n // 2] if n % 2 else (devs[n // 2 - 1] + devs[n // 2]) / 2
    sigma = 1.4826 * mad
    if sigma <= 0:
        mean = sum(v) / n
        sigma = (sum((x - mean) ** 2 for x in v) / n) ** 0.5 or 1.0
    return med, sigma


MIN_HET_VARIANTS = 50  # below this variant count the het fraction is too noisy to interpret -> None


def het_frac(m):
    """Fraction of VARIANT genotypes that are heterozygous (het/(het+hom)). Illustrative mixture signal; the
    denominator is variant-count (not callable sites), so it is entangled with divergence -> read with care.
    Returns None below a minimum variant count (too noisy to interpret)."""
    het = to_float(m.get("het_variants"))
    hom = to_float(m.get("homo_variants"))
    if het is None or hom is None or (het + hom) < MIN_HET_VARIANTS:
        return None
    return het / (het + hom)


def lineage_counts_parsed(lc):
    """Parse 'L4.1:3200;L4.3:120;L2:40' into {main_lineage: total_count}, COLLAPSING sub-lineages (L4.1 -> L4)
    so two sub-lineages of the same main lineage do not spuriously read as a mixed infection."""
    if not lc or str(lc) in ("NA", ".", ""):
        return {}
    out = {}
    for part in str(lc).replace(",", ";").split(";"):
        if ":" in part:
            k, v = part.rsplit(":", 1)
            n = to_float(v)
            if n is not None:
                main = k.strip().split(".")[0]   # collapse sub-lineage to main lineage
                out[main] = out.get(main, 0.0) + n
    return out


def lineage_fracs(lc):
    """{main_lineage: fraction} from the collapsed counts."""
    counts = lineage_counts_parsed(lc)
    tot = sum(counts.values())
    return {k: v / tot for k, v in counts.items()} if tot > 0 else {}


def is_ancient(m):
    """True if the summary row marks this sample ancient (aDNA). Missing column -> False (modern)."""
    return str(m.get("ancient", "")).strip().lower() in ("yes", "y", "true", "1", "ancient")


def infer_metric(key, values):
    """Infer (label, kind, dir) for an un-registered numeric summary column. Neutral dir = no auto-flag / no claim."""
    kl = key.lower()
    if any(h in kl for h in _PCT_HINTS):
        kind = "pct"
    elif all(v is None or float(v).is_integer() for v in values) and any(h in kl for h in _INT_HINTS):
        kind = "int"
    else:
        kind = "float"
    label = key.replace("_pct", " %").replace("_", " ").strip()
    label = label[:1].upper() + label[1:]
    return label, kind, "neu"


def discover_extra_metrics(summ):
    """Any summary column not a KNOWN_KEY and numeric for >=1 sample -> an auto 'extra' metric (future-proofing)."""
    if not summ:
        return []
    seen = set()
    for mm in summ.values():
        for k in mm.keys():
            if k not in KNOWN_KEYS and k != "sample_id":
                seen.add(k)
    extras = []
    for k in sorted(seen):
        vals = [to_float(mm.get(k)) for mm in summ.values()]
        if any(v is not None for v in vals):
            label, kind, direction = infer_metric(k, vals)
            extras.append({"key": k, "label": label, "kind": kind, "dir": direction})
    return extras


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


def flag_sample(m, thr, snp_med, snp_sig, ancient=False, anc_thr=None, dmg=None):
    flags = []
    at = anc_thr if (ancient and anc_thr) else thr   # ancient samples judged against the aDNA threshold view
    g = lambda k: to_float(m.get(k))
    # produced-nothing guard: a sample with no consensus AND no mapping metrics is a FAIL, not a quiet PASS-by-absence
    if all(g(k) is None for k in ("mean_depth", "breadth_pct", "missing_pct", "mapped_pct")):
        return "FAIL", ["NO_DATA"]
    if (v := g("mean_depth")) is not None and v < at["depth_min"]:
        flags.append("LOW_DEPTH")
    if (v := g("breadth_pct")) is not None and v < at["breadth_min"]:
        flags.append("LOW_BREADTH")
    if (v := g("missing_pct")) is not None and v > at["missing_max"]:
        flags.append("HIGH_MISSING")
    if (v := g("mapped_pct")) is not None and v < at["mapping_min"]:
        flags.append("MAPPING_LOW")
    if (v := g("duplication_pct")) is not None and v > at["dup_max"]:
        flags.append("HIGH_DUP")
    if (v := g("iupac_pct")) is not None and v > at["iupac_max"]:
        flags.append("HIGH_IUPAC")
    if thr["titv_min"] > 0 and (v := g("ti_tv")) is not None and v < thr["titv_min"]:
        flags.append("TITV_LOW")
    # SNP count is dominated by phylogenetic divergence, not quality -> a cohort distributional outlier (WARN only).
    # snp_med/sig are computed over MODERN samples; ancient SNP counts are low by design, so skip the z-flag for aDNA.
    if not ancient:
        snp = g("snps")
        if snp is not None and snp_med is not None and snp_sig:
            if snp < snp_med - thr["snp_z"] * snp_sig:
                flags.append("SNP_LOW")
            elif snp > snp_med + thr["snp_z"] * snp_sig:
                flags.append("SNP_HIGH")
    # het fraction: NOT applied to aDNA (post-mortem deamination inflates apparent het -> would false-flag mummies)
    if not ancient:
        hf = het_frac(m)
        if hf is not None and hf * 100.0 > thr["het_max_frac"]:
            flags.append("HET_HIGH")
    # MIXED: >1 MAIN lineage above the fraction cut-off AND above a minimum marker count (guards single-read noise)
    counts = lineage_counts_parsed(m.get("lineage_counts"))
    tot = sum(counts.values())
    if tot > 0 and sum(1 for c in counts.values()
                       if (c / tot) * 100.0 >= thr["mixed_min_frac"] and c >= 3) > 1:
        flags.append("MIXED")
    # aDNA authentication: an ANCIENT sample whose terminal 5' C>T is below the floor -> possible modern contaminant
    if ancient and anc_thr is not None and dmg and dmg.get("ct1") is not None \
            and dmg["ct1"] < anc_thr["damage_min_ct"]:
        flags.append("DAMAGE_LOW")
    verdict = "FAIL" if any(f in FAIL_FLAGS for f in flags) else ("WARN" if flags else "PASS")
    return verdict, flags


# ============================================================================= CSS / JS / SHELL
CSS = r"""
:root{--bg:#eef2f7;--panel:#fff;--line:#e6ebf2;--ink:#15202e;--mut:#6a7889;--soft:#f6f8fb;
 --accent:#0e8ba8;--accent-soft:#e3f3f7;
 --pass:#0f9d6b;--warn:#dd8a1a;--fail:#e23a4a;--good:#16a37a;--bad:#e5615c;--neu:#4a90b8;
 --bandfill:rgba(15,157,107,.09);--bandedge:rgba(15,157,107,.42);
 --sh:0 1px 2px rgba(16,24,40,.04),0 3px 8px rgba(16,24,40,.05);--r:15px}
*{box-sizing:border-box} html{scroll-behavior:smooth}
body{margin:0;font-family:ui-sans-serif,system-ui,-apple-system,"Segoe UI",Roboto,Helvetica,sans-serif;color:var(--ink);font-size:14px;line-height:1.5;-webkit-font-smoothing:antialiased;
 background-color:#eef1f6;
 background-image:radial-gradient(1200px 520px at 50% -260px,#ffffff,transparent 70%),linear-gradient(180deg,#f4f6fa 0%,#e9edf3 100%);
 background-attachment:fixed,fixed;background-repeat:no-repeat,no-repeat}
.num,td{font-variant-numeric:tabular-nums} :focus-visible{outline:2px solid var(--accent);outline-offset:1px}
header{position:sticky;top:0;z-index:40;background:rgba(255,255,255,.82);-webkit-backdrop-filter:saturate(1.15) blur(10px);backdrop-filter:saturate(1.15) blur(10px);color:var(--ink);padding:13px 22px;display:flex;align-items:center;gap:14px;border-bottom:1px solid var(--line);box-shadow:0 1px 0 rgba(16,24,40,.02),0 8px 24px rgba(16,24,40,.035)}
header .logo{font-weight:700;letter-spacing:-.2px;font-size:16px;display:flex;align-items:center;gap:9px;color:var(--ink)}
header .logo b{color:var(--accent);font-weight:700}
header .logo .dot{width:9px;height:9px;border-radius:50%;background:var(--accent);box-shadow:0 0 0 3px rgba(14,139,168,.16)}
header .meta{color:var(--mut);font-size:12px}
nav{margin-left:auto;display:flex;gap:2px;overflow:auto} nav::-webkit-scrollbar{display:none}
nav a{color:#5b6b7e;text-decoration:none;font-size:12.5px;padding:5px 10px;border-radius:8px;white-space:nowrap;font-weight:500} nav a:hover{color:var(--accent);background:var(--accent-soft)}
.wrap{max-width:1180px;margin:0 auto;padding:24px 22px 90px}
section{margin-top:34px} section:first-of-type{margin-top:24px}
h2{font-size:11.5px;text-transform:uppercase;letter-spacing:.09em;color:var(--mut);font-weight:600;margin:0 0 14px;display:flex;align-items:center;gap:12px;flex-wrap:wrap}
h2 .c{text-transform:none;letter-spacing:0;font-weight:400;font-size:12.5px;color:#95a3b4}
.panel{background:var(--panel);border:1px solid var(--line);border-radius:var(--r);box-shadow:var(--sh);overflow:hidden}
.pad{padding:16px 18px}
/* hero */
.hero{display:flex;gap:16px;flex-wrap:wrap;align-items:stretch}
.summary{display:flex;align-items:center;gap:22px;background:var(--panel);border:1px solid var(--line);border-radius:var(--r);box-shadow:var(--sh);padding:16px 24px 16px 20px}
.summary svg{flex:none} .counts{display:flex;gap:24px}
.counts .c .n{font-size:27px;font-weight:700;line-height:1;letter-spacing:-.5px} .counts .c .l{font-size:10px;text-transform:uppercase;letter-spacing:.08em;color:var(--mut);margin-top:4px}
.counts .all .n{color:var(--ink)} .counts .pass .n{color:var(--pass)} .counts .warn .n{color:var(--warn)} .counts .fail .n{color:var(--fail)}
.chips{display:flex;gap:8px;flex-wrap:wrap;align-content:center;background:var(--panel);border:1px solid var(--line);border-radius:var(--r);box-shadow:var(--sh);padding:14px 16px;flex:1;min-width:220px}
.chips .t{font-size:10px;color:var(--mut);text-transform:uppercase;letter-spacing:.08em;align-self:center;margin-right:2px}
.chip{font-size:12px;padding:3px 12px;border-radius:20px;border:1px solid var(--line);background:var(--soft);color:#516074;cursor:pointer;user-select:none;transition:.12s;font-weight:500}
.chip:hover{border-color:var(--accent);color:var(--accent);background:var(--accent-soft)} .chip.on{background:var(--accent);color:#fff;border-color:var(--accent)}
.chip .k{opacity:.65;margin-left:5px;font-variant-numeric:tabular-nums}
/* controls */
.controls{display:flex;gap:9px;align-items:center;flex-wrap:wrap;margin-bottom:12px}
.controls input[type=search]{padding:7px 12px;border:1px solid var(--line);border-radius:9px;font-size:13px;min-width:200px;background:#fff;box-shadow:var(--sh)}
.controls input[type=search]:focus{border-color:var(--accent)}
.gsearch{padding:4px 10px;border:1px solid var(--line);border-radius:8px;font-size:12px;min-width:110px;max-width:170px;background:#fff;box-shadow:var(--sh)}
.gsearch:focus{border-color:var(--accent);outline:none}
.panel.expanded .gtable{max-height:82vh}
.controls label{font-size:12.5px;color:var(--mut);display:flex;align-items:center;gap:6px;cursor:pointer}
.btn{font-size:12.5px;color:#33465c;border:1px solid var(--line);border-radius:9px;padding:7px 13px;background:#fff;cursor:pointer;box-shadow:var(--sh)}
.btn:hover{border-color:var(--accent);color:var(--accent)}
details.dd{position:relative} details.dd>summary{cursor:pointer;font-size:12.5px;color:#33465c;list-style:none;border:1px solid var(--line);border-radius:9px;padding:7px 13px;background:#fff;box-shadow:var(--sh)}
details.dd>summary::-webkit-details-marker{display:none} details.dd>summary::marker{content:""}
details.dd .menu{position:absolute;z-index:50;background:#fff;border:1px solid var(--line);border-radius:12px;padding:12px 16px;box-shadow:0 14px 40px rgba(16,24,40,.16);max-height:340px;overflow:auto;margin-top:6px;min-width:170px}
details.dd .menu label{display:block;padding:3px 0;font-size:12.5px;color:#33465c}
.seg{display:inline-flex;border:1px solid var(--line);border-radius:9px;overflow:hidden;box-shadow:var(--sh)}
.seg button{border:0;background:#fff;color:#516074;font-size:12.5px;padding:6px 13px;cursor:pointer;border-right:1px solid var(--line)}
.seg button:last-child{border-right:0} .seg button.on{background:var(--accent);color:#fff}
select.msel{font-size:12.5px;border:1px solid var(--line);border-radius:9px;padding:6px 10px;background:#fff;color:#33465c;box-shadow:var(--sh)}
.hint{margin-left:auto;color:var(--mut);font-size:12px}
/* table */
.gtable{overflow:auto;max-height:76vh;border-radius:var(--r)}
.gtable::-webkit-scrollbar{height:10px;width:10px} .gtable::-webkit-scrollbar-thumb{background:#cfd8e3;border-radius:6px;border:2px solid #fff}
table{border-collapse:separate;border-spacing:0;width:100%;font-size:12.5px}
th,td{padding:8px 12px;white-space:nowrap;border-bottom:1px solid #eef2f6;text-align:right}
th{background:#f7f9fc;color:#556579;cursor:pointer;user-select:none;position:sticky;top:0;font-weight:600;font-size:11px;text-transform:uppercase;letter-spacing:.03em;z-index:4;box-shadow:0 1px 0 var(--line)}
th:hover{background:#eef3f8} td.na{color:#c4ccd7}
tbody tr{cursor:pointer;transition:background .1s} tbody tr:hover td{background:#fafcfe}
tr.hl td{background:#fff8e1!important}
th.s,td.s{text-align:left;position:sticky;left:0;background:var(--panel);border-right:1px solid var(--line);z-index:3;font-weight:600}
th.s{z-index:6;background:#f7f9fc} tr.hl td.s{background:#fff3ce!important}
.v{font-weight:600;padding:2px 10px;border-radius:20px;font-size:11px;display:inline-block;letter-spacing:.02em}
.v.PASS{color:#0b7350;background:#dff5ec} .v.WARN{color:#95560d;background:#fdefd6} .v.FAIL{color:#a01f2d;background:#fde3e6}
.flags{color:var(--mut);font-size:11.5px;text-align:left;white-space:normal}
/* plots */
.bee{display:flex;align-items:center;border-bottom:1px solid #f2f5f9;height:40px} .bee:last-child{border:0}
.bl{flex:0 0 150px;padding:0 14px;font-size:12px;color:#516074;text-align:right;font-weight:500} .nd{color:#c4ccd7;font-size:12px;padding-left:14px}
.plotwrap{flex:1} .plotwrap svg{display:block} .plotwrap circle,.plotwrap rect.hit{cursor:pointer}
.stack{display:flex;align-items:center;border-bottom:1px solid #f4f7fa;height:26px} .stack:last-child{border:0}
.stack .sl{flex:0 0 150px;padding:0 14px;font-size:11.5px;color:#516074;text-align:right;overflow:hidden;text-overflow:ellipsis;white-space:nowrap}
.stack .sb{flex:1;height:13px;display:flex;border-radius:4px;overflow:hidden;background:#eef2f7;box-shadow:inset 0 0 0 1px rgba(16,24,40,.03)}
.stack .sv{flex:0 0 52px;text-align:right;padding-right:14px;font-size:11px;color:#9aa7b6;font-variant-numeric:tabular-nums}
.legend{display:flex;gap:18px;font-size:11.5px;color:var(--mut);padding:10px 16px;border-top:1px solid var(--line);flex-wrap:wrap}
.legend i{display:inline-block;width:10px;height:10px;border-radius:50%;margin-right:6px;vertical-align:middle}
.footer{color:var(--mut);font-size:11.5px;margin-top:30px;border-top:1px solid var(--line);padding-top:14px;line-height:1.8}
#tt{position:fixed;pointer-events:none;background:#0f2431;color:#fff;font-size:12px;line-height:1.5;padding:7px 11px;border-radius:9px;opacity:0;transition:opacity .08s;z-index:80;white-space:nowrap;box-shadow:0 8px 24px rgba(15,36,49,.32)}
/* curation basket */
.cbx{vertical-align:middle;margin-right:7px;accent-color:var(--accent);cursor:pointer;width:14px;height:14px}
.sname{cursor:pointer;border-bottom:1px dashed transparent;transition:.12s} .sname:hover{color:var(--accent);border-bottom-color:var(--accent)}
th.s .hlab{font-weight:600}
.curation{display:flex;align-items:center;gap:14px;flex-wrap:wrap;margin-top:12px;padding:12px 16px;background:linear-gradient(180deg,#fff,var(--soft));border:1px solid var(--line);border-radius:12px;box-shadow:var(--sh)}
.cur-intro{flex-basis:100%;font-size:11.5px;color:var(--mut);line-height:1.55;margin-bottom:2px} .cur-intro b{color:#33465c}
.cur-read{font-size:13px;color:#33465c} .cur-read b{color:var(--ink);font-variant-numeric:tabular-nums;font-size:15px} .cur-read .arw{color:var(--accent);margin:0 3px;font-weight:700}
/* expand-to-fill (fullscreen-within-window) */
.exp-h{display:inline-flex;align-items:center;gap:5px;font-size:11.5px;color:#33465c;border:1px solid var(--line);border-radius:8px;padding:4px 9px;background:#fff;cursor:pointer;box-shadow:var(--sh)}
.exp-h:hover{border-color:var(--accent);color:var(--accent)}
.panel.expanded{position:fixed;inset:14px;z-index:96;overflow:auto;box-shadow:0 24px 80px rgba(16,24,40,.34);max-height:none}
.panel.expanded .gtable{max-height:calc(100vh - 60px)}
body.has-expanded{overflow:hidden} body.has-expanded::after{content:"";position:fixed;inset:0;background:rgba(15,36,49,.42);backdrop-filter:blur(1px);z-index:95}
.exp-close{display:none;position:fixed;top:20px;right:24px;z-index:98;align-items:center;gap:6px;font-size:13px;font-weight:600;color:#fff;background:#0f2431;border:0;border-radius:10px;padding:9px 15px;cursor:pointer;box-shadow:0 8px 24px rgba(15,36,49,.4)}
.exp-close:hover{background:var(--accent)} body.has-expanded .exp-close{display:inline-flex}
.cur-btns{display:flex;gap:7px;flex-wrap:wrap;margin-left:auto}
.btn.prim{background:var(--accent);color:#fff;border-color:var(--accent)} .btn.prim:hover{background:#0b7387;color:#fff}
.btn.on{background:var(--accent);color:#fff;border-color:var(--accent)}
/* per-sample detail modal */
.modal{position:fixed;inset:0;background:rgba(15,36,49,.42);backdrop-filter:blur(2px);display:none;align-items:flex-start;justify-content:center;z-index:90;padding:5vh 16px;overflow:auto}
.modal.open{display:flex}
.modalcard{background:var(--panel);border-radius:16px;box-shadow:0 24px 70px rgba(15,36,49,.4);max-width:560px;width:100%;padding:22px 24px;position:relative;animation:pop .16s ease-out}
@keyframes pop{from{transform:translateY(8px);opacity:.4}to{transform:none;opacity:1}}
.modalx{position:absolute;top:12px;right:14px;border:0;background:#eef2f7;color:#516074;width:28px;height:28px;border-radius:50%;font-size:18px;line-height:1;cursor:pointer} .modalx:hover{background:#e0e6ee}
.dhead{display:flex;align-items:flex-start;justify-content:space-between;gap:12px;margin-bottom:12px} .dtitle{font-size:18px;font-weight:700;letter-spacing:-.3px} .dmeta{font-size:12px;color:var(--mut);margin-top:2px}
.dflags{display:flex;gap:6px;flex-wrap:wrap;margin-bottom:6px} .dflags .chip{cursor:default} .dflags .chip.failc{background:#fde3e6;color:#a01f2d;border-color:#f5c2c8}
.dsub{font-size:10.5px;text-transform:uppercase;letter-spacing:.08em;color:var(--mut);font-weight:600;margin:16px 0 8px}
.lcomp{display:flex;flex-direction:column;gap:5px;margin-bottom:2px}
.lrow{display:flex;align-items:center;gap:9px;font-size:12px} .lk{flex:0 0 58px;color:#33465c;font-weight:600} .lbarw{flex:1;height:9px;background:#eef2f7;border-radius:5px;overflow:hidden} .lbar{height:100%;background:var(--accent)} .lv{flex:0 0 46px;text-align:right;color:#8895a6;font-variant-numeric:tabular-nums}
.drow{display:flex;align-items:center;gap:10px;padding:3px 0;font-size:12px;border-bottom:1px solid #f4f7fa} .drow:last-of-type{border:0}
.dk{flex:0 0 120px;color:#516074} .dbarwrap{flex:1;height:7px;background:#eef2f7;border-radius:4px;overflow:hidden} .dbar{height:100%;border-radius:4px} .dv{flex:0 0 74px;text-align:right;font-weight:600;font-variant-numeric:tabular-nums} .dp{flex:0 0 34px;text-align:right;color:#9aa7b6;font-size:11px}
.dbtns{margin-top:16px;display:flex;justify-content:flex-end}
/* metric help */
#helpmenu{max-width:340px;white-space:normal}
.hitem{padding:6px 0;border-bottom:1px solid #f0f3f7} .hitem:last-child{border:0} .hitem b{font-size:12px} .hk{color:#aab6c4;font-size:10.5px} .hd{font-size:11.5px;color:#516074;margin-top:2px} .hr{color:var(--accent)}
/* provenance bar */
.provbar{display:flex;align-items:center;gap:12px;margin-top:12px}
.prov{display:flex;gap:8px;flex-wrap:wrap;flex:1;font-size:11px;color:var(--mut)}
.prov span{background:var(--panel);border:1px solid var(--line);border-radius:7px;padding:3px 9px;font-variant-numeric:tabular-nums}
/* beeswarm cohort sub-label + bands */
.bl .stat{display:block;font-weight:400;font-size:9.5px;color:#9aa7b6;font-variant-numeric:tabular-nums;margin-top:1px;line-height:1.1;overflow:hidden;text-overflow:ellipsis}
/* lineage */
.ldot{display:inline-block;width:9px;height:9px;border-radius:50%;margin-right:8px;vertical-align:middle}
tr.lingrp td{background:#f0f5f9;color:#33465c;font-weight:600;font-size:11px;letter-spacing:.03em;position:sticky;left:0} tr.lingrp{cursor:default}
.abadge{margin-left:7px;font-size:9px;font-weight:700;letter-spacing:.04em;color:#8a5a12;background:#fdefd6;border:1px solid #f0d9a8;border-radius:6px;padding:1px 5px;vertical-align:middle}
.lincomp-bar{display:flex;height:20px;border-radius:6px;overflow:hidden;margin:14px 16px 6px;box-shadow:inset 0 0 0 1px rgba(16,24,40,.05)}
.lincomp-bar .lseg{height:100%;min-width:2px;transition:filter .1s} .lincomp-bar .lseg:hover{filter:brightness(1.08)}
.lincomp-lab{display:flex;gap:12px;flex-wrap:wrap;padding:0 16px 12px;font-size:11.5px;color:var(--mut)}
.lincomp-lab .lchip i{display:inline-block;width:9px;height:9px;border-radius:2px;margin-right:5px;vertical-align:middle}
.lincomp-lab .lchip b{color:var(--ink);font-variant-numeric:tabular-nums;margin-left:2px}
#linsumtable tbody tr{cursor:pointer} #linsumtable tbody tr:hover td{background:#fafcfe}
/* genome landscape */
.genome-scroll{overflow-x:auto;overflow-y:hidden} .genome-scroll svg{display:block;cursor:crosshair} .genome-scroll rect{shape-rendering:crispEdges}
#scatter svg{cursor:crosshair} #corr_body{overflow-x:auto} #corrsvg text{pointer-events:none}
.hot-note{font-size:11.5px;color:var(--mut);line-height:1.5;padding:12px 16px 6px}
#hottable tbody tr{cursor:pointer} #hottable tbody tr:hover td{background:#fafcfe}
/* aDNA panel */
.anote{font-size:12px;color:#516074;line-height:1.6;background:var(--soft);border:1px solid var(--line);border-radius:10px;padding:10px 13px;margin-bottom:12px} .anote b{color:var(--ink)}
.acards{display:grid;grid-template-columns:repeat(auto-fill,minmax(230px,1fr));gap:11px}
.acard{border:1px solid var(--line);border-radius:12px;padding:11px 13px;background:linear-gradient(180deg,#fff,var(--soft))}
.acard.low{border-color:#f0c2be;background:linear-gradient(180deg,#fff,#fdf1f0)}
.acard-h{display:flex;align-items:center;justify-content:space-between;gap:8px;margin-bottom:7px} .acard-h .sname{font-weight:600;font-size:12.5px}
.acard-row{display:flex;align-items:center;gap:12px} .aspark{flex:0 0 auto} .astats{display:flex;flex-direction:column;gap:3px;flex:1}
.astat{display:flex;justify-content:space-between;font-size:11.5px} .ak{color:#8895a6} .av{font-weight:600;font-variant-numeric:tabular-nums}
.av.good{color:#0b7350} .av.bad{color:#a01f2d} .av.na{color:#c4ccd7}
.alow{margin-top:8px;font-size:11px;color:#a01f2d;line-height:1.45}
/* modal 'why' block */
.dwhy{font-size:11.5px;color:#516074;background:var(--soft);border:1px solid var(--line);border-radius:9px;padding:8px 11px;margin-bottom:6px;line-height:1.55} .dwhy b{color:#a01f2d}
@media print{
  header,nav,.controls,.provbar #printBtn,#thbox,#athbox,.dd,.chips{display:none!important}
  body{background:#fff} .wrap{max-width:none;padding:0} .gtable{max-height:none!important;overflow:visible!important}
  section{break-inside:avoid} .panel{box-shadow:none}
}
@media (max-width:760px){
  header{padding:11px 16px} .wrap{padding:16px 14px 70px}
  .hero{flex-direction:column} .summary{justify-content:space-between;padding:14px 18px} .counts{gap:16px}
  .counts .c .n{font-size:22px} .chips{min-width:0}
  section{margin-top:26px}
  th,td{padding:7px 9px} .bl,.stack .sl{flex-basis:96px;font-size:10.5px;padding:0 8px} .stack .sv{flex-basis:42px;padding-right:8px}
  .hint{margin-left:0;flex-basis:100%}
  .cur-btns{margin-left:0} .curation{gap:10px} .modalcard{padding:18px 15px} .dk{flex-basis:84px} .dv{flex-basis:58px} .dp{flex-basis:28px}
  details.dd .menu{min-width:0;max-width:calc(100vw - 28px);box-sizing:border-box} #thbox{min-width:0!important}
}
@media (prefers-reduced-motion:reduce){*{transition:none!important;scroll-behavior:auto!important}}
/* SNP dynamics panel */
.dyngrp{margin:12px 0;padding:10px 12px;border:1px solid var(--line);border-radius:10px;background:var(--soft)}
.dynhd{margin-bottom:6px}
.dynleg{display:flex;gap:14px;flex-wrap:wrap;align-items:center;font-size:11px;color:var(--mut);margin-bottom:8px}
.dynleg i{display:inline-block;width:11px;height:11px;border-radius:3px;margin-right:5px;vertical-align:-1px}
.dyntbl{width:100%;border-collapse:collapse;font-size:11px;margin-top:8px}
.dyntbl th,.dyntbl td{text-align:left;padding:3px 8px;border-bottom:1px solid var(--line)}
.dyntbl th{color:var(--mut);font-weight:600}
.dtraj{font-variant-numeric:tabular-nums;color:var(--ink)}
.dchip{display:inline-block;color:#fff;border-radius:5px;padding:1px 6px;font-size:10px;margin-right:3px}
"""

JS = r"""
(function(){
var R=REPORT;
var DIST=R.dist;
var VCOL={PASS:'#94a3b8',WARN:'#d97706',FAIL:'#dc2626'}, VFILL={PASS:'#16a34a',WARN:'#d97706',FAIL:'#dc2626'};
var BAR={hi_good:'#22a06b',hi_bad:'#e0544f',neu:'#4f83c2'};
// ---- lineage palette: deterministic, colour-blind-safe, self-contained (Okabe-Ito + Tol-muted; golden-angle overflow)
var LINPAL_BASE=['#4477aa','#ee6677','#228833','#ccbb44','#66ccee','#aa3377','#e69f00','#0072b2','#d55e00','#009e73','#cc79a7','#882255'];
R.lin_colors=R.lin_colors||{};   // canonical palette (e.g. mycolorsTB) keyed by lineage label
var LINCOL={};(function(){(R.lineages||[]).forEach(function(lab,i){
  if(R.lin_colors[lab]){LINCOL[lab]=R.lin_colors[lab];}                       // canonical colour if provided
  else if(i<LINPAL_BASE.length){LINCOL[lab]=LINPAL_BASE[i];}                  // else colour-blind-safe fallback
  else{var h=(i*137.508)%360;LINCOL[lab]='hsl('+h.toFixed(0)+',58%,52%)';}});})();
function linColor(lab){return (lab!=null&&LINCOL[lab])?LINCOL[lab]:'#b8c2cf';}
var thr=Object.assign({},R.thresholds);
var athr=Object.assign({},R.anc_thresholds||{});      // ancient (aDNA) threshold view
function actv(s){return (s.anc&&R.n_ancient)?athr:thr;}   // active threshold set for a sample
R.defs=R.defs||{}; R.extra=R.extra||[]; R.provenance=R.provenance||{};
// ---- derived metrics registered as first-class metrics ----
function nlin(s){if(!s.linc)return null;var n=0,tot=0,k;for(k in s.linc)tot+=s.linc[k];if(tot<=0)return null;for(k in s.linc){if((s.linc[k]/tot)*100>=thr.mixed_min_frac&&s.linc[k]>=3)n++;}return n;}
R.samples.forEach(function(s){s.m.het_frac=(s.hetf!=null)?s.hetf*100:null;s.m.n_lineages=nlin(s);});
R.metrics.push({key:'het_frac',label:'Het %',kind:'pct',dir:'hi_bad'});
R.metrics.push({key:'n_lineages',label:'# lin',kind:'int',dir:'hi_bad'});
if(DIST.indexOf('het_frac')<0)DIST.push('het_frac');
// ---- SNP density per callable Mb (guarded: only for a plausible genome length) ----
if(R.snp_density_ok){var GMB=R.genome_len/1e6;R.samples.forEach(function(s){var snp=s.m.snps,cal=s.m.callable_pct;s.m.snp_density=(snp!=null&&cal!=null&&cal>0)?(snp/(cal/100*GMB)):null;});R.metrics.push({key:'snp_density',label:'SNP/Mb',kind:'float',dir:'neu'});if(DIST.indexOf('snp_density')<0)DIST.push('snp_density');}
// ===================== SIGNATURE ANALYSIS LAYER (QC-space PCA + robust Mahalanobis + divergence/completeness + temporal) =====================
R.defs=R.defs||{};
var PCA_CAPTION='QC-metric space. Samples ordinated by their quality-metric profiles (each metric robustly z-scored across the samples in view, then projected onto the top two principal components of the metric covariance). This is a descriptive map of how QC profiles co-vary, not a map of biological relatedness and not a test. PC1 usually captures a coverage-and-completeness composite because depth, breadth and callability are strongly correlated; see the loadings for what drives each axis. A sample sitting apart can reflect genuine sequence divergence, a library or coverage artifact, or contamination, and the ordination alone cannot separate these, so cross-read the divergence-vs-completeness panel below. Axes are cohort-relative and recompute when you filter. The shaded ellipse marks the data spread (2 SD along each PC), not a confidence or significance region. Rings mark robust-Mahalanobis outliers.';
var REFBIAS_CAPTION='Reference divergence vs completeness. Each sample is plotted by consensus incompleteness (x) against a proxy for how far it sits from the reference (y = SNPs, per callable Mb when available to avoid coupling the two axes). This is a screen that separates two ways a sample can look reference-like. Bottom-left (few SNPs, little missing) is the reference-bias suspect zone: low apparent divergence that is not explained by missing data, consistent with reference-guided calling collapsing a divergent sample onto the reference, but a genuinely reference-close or smaller-genome sample looks identical here, so this is a list to check, not a verdict. Bottom-right (few SNPs, lots missing) is honest low coverage: the low count is just missing data, not bias. High divergence is not automatically real either, since contamination or a wrong reference also inflates SNP counts, so cross-read the mapping and error-rate metrics. Guide lines are the live missing gate (x) and a robust low-divergence cut (y = cohort median minus 2.5 robust SD); only samples below that cut are called out, so a normally-divergent sample is never flagged here. Confirm any suspected reference bias downstream; this panel only points.';
var TEMPORAL_CAPTION='Temporal sampling overview. Per group, the span and count of samples carrying a parseable collection year, alongside a rough sequence-information proxy (median SNPs vs reference). {DATED} of {TOT} samples carry a parseable year; overall span {SPAN}. This describes whether the cohort has the dated spread a downstream molecular-dating analysis would need, and it does not test for clock signal. Wide sampling in time is necessary but not sufficient: a cohort can span decades and still show no root-to-tip temporal signal, which only a proper date-randomization or regression test downstream can establish. The SNP median is an information proxy, not a count of clock-informative sites.';
// core quality/completeness metrics for the ordination + Mahalanobis. Divergence (snps/snp_density) is DELIBERATELY
// excluded here (a divergent lineage is not a QC outlier); it is the y-axis of the divergence-vs-completeness panel.
var PCA_CORE=['mean_depth','breadth_pct','callable_pct','missing_pct','iupac_pct','mapped_pct','duplication_pct','coverage_cv','mapq0_pct','error_rate_pct','q30_pct','het_frac','ti_tv','insert_size'];
function pcaColUsable(pool,k){var v=[],i;for(i=0;i<pool.length;i++){var x=pool[i].m[k];if(x!=null&&isFinite(x))v.push(x);}
  if(v.length<Math.max(4,Math.ceil(0.5*pool.length)))return null;   // require >50% present AND >=4 observations
  v.sort(function(a,b){return a-b;});var md=med2(v);
  var d=v.map(function(x){return Math.abs(x-md);}).sort(function(a,b){return a-b;});var sig=1.4826*med2(d);
  if(!(sig>0)){var mn=v.reduce(function(a,b){return a+b;},0)/v.length;sig=Math.sqrt(v.reduce(function(a,b){return a+(b-mn)*(b-mn);},0)/v.length);}
  if(!(sig>0))return null;return {med:md,sig:sig};}                 // truly constant -> drop the column
function buildZ(pool){var keys=[],meds=[],sigs=[],c;
  for(c=0;c<PCA_CORE.length;c++){var k=PCA_CORE[c],u=pcaColUsable(pool,k);if(u){keys.push(k);meds.push(u.med);sigs.push(u.sig);}}
  var K=keys.length,Z=[],samples=[],i,j;
  for(i=0;i<pool.length;i++){var s=pool[i],row=new Array(K);
    for(j=0;j<K;j++){var x=s.m[keys[j]];if(x==null||!isFinite(x))x=meds[j];var z=(x-meds[j])/sigs[j];if(z>8)z=8;else if(z<-8)z=-8;row[j]=z;}
    Z.push(row);samples.push(s);}
  return {Z:Z,keys:keys,samples:samples,k:K,n:Z.length};}
function covMatrix(Z,k){var n=Z.length,C=[],i,j,r;
  for(i=0;i<k;i++){C[i]=new Array(k);for(j=0;j<k;j++)C[i][j]=0;}
  for(r=0;r<n;r++){var row=Z[r];for(i=0;i<k;i++){var zi=row[i];for(j=i;j<k;j++)C[i][j]+=zi*row[j];}}
  var d=(n>1)?(n-1):1;for(i=0;i<k;i++)for(j=i;j<k;j++){C[i][j]/=d;C[j][i]=C[i][j];}return C;}
function topEigs(C,k,howMany){
  function matVec(M,x){var y=new Array(k),i,j;for(i=0;i<k;i++){var s=0;for(j=0;j<k;j++)s+=M[i][j]*x[j];y[i]=s;}return y;}
  function nrm(x){var s=0,i;for(i=0;i<k;i++)s+=x[i]*x[i];return Math.sqrt(s);}
  var A=[],a;for(a=0;a<k;a++)A[a]=C[a].slice();var vecs=[],vals=[],e,i,j;
  for(e=0;e<howMany;e++){var x=new Array(k);for(i=0;i<k;i++)x[i]=(((i*2654435761)>>>0)%1000)/1000-0.5+1e-6*(e+1);
    var nx=nrm(x)||1;for(i=0;i<k;i++)x[i]/=nx;var lam=0,it;
    for(it=0;it<500;it++){var y=matVec(A,x),ny=nrm(y);if(ny<1e-14)break;for(i=0;i<k;i++)y[i]/=ny;
      var Ay=matVec(A,y),nl=0;for(i=0;i<k;i++)nl+=y[i]*Ay[i];var dot=0;for(i=0;i<k;i++)dot+=y[i]*x[i];if(dot<0)for(i=0;i<k;i++)y[i]=-y[i];
      x=y;if(Math.abs(nl-lam)<1e-9*(Math.abs(nl)+1e-12)){lam=nl;break;}lam=nl;}
    if(lam<0)lam=0;vecs.push(x);vals.push(lam);
    for(i=0;i<k;i++)for(j=0;j<k;j++)A[i][j]-=lam*x[i]*x[j];}     // Hotelling deflation
  return {vecs:vecs,vals:vals};}
function pca2(pool){var B=buildZ(pool);
  if(B.n<4)return {ok:false,reason:'need at least 4 samples with a shared core metric in view',B:B};
  if(B.k<2)return {ok:false,reason:'need at least 2 core QC metrics in view',B:B};
  var C=covMatrix(B.Z,B.k),E=topEigs(C,B.k,2),tot=0,i;for(i=0;i<B.k;i++)tot+=C[i][i];
  if(!(tot>0))return {ok:false,reason:'no QC variance (samples identical on the core metrics)',B:B};
  var pev=[100*E.vals[0]/tot,100*E.vals[1]/tot];
  var scores=B.Z.map(function(row){var p1=0,p2=0,j;for(j=0;j<B.k;j++){p1+=row[j]*E.vecs[0][j];p2+=row[j]*E.vecs[1][j];}return [p1,p2];});
  function loadTop(vec){return B.keys.map(function(k,j){return {k:k,label:(MET[k]||{}).label||k,w:vec[j]};}).sort(function(a,b){return Math.abs(b.w)-Math.abs(a.w);}).slice(0,5);}
  return {ok:true,B:B,scores:scores,pev:pev,vecs:E.vecs,load:[loadTop(E.vecs[0]),loadTop(E.vecs[1])]};}
function invRidge(C,k,lambda){var A=[],I=[],i,j,r;
  for(i=0;i<k;i++){A[i]=C[i].slice();A[i][i]+=lambda;I[i]=new Array(k);for(j=0;j<k;j++)I[i][j]=(i==j)?1:0;}
  for(i=0;i<k;i++){var p=i,mx=Math.abs(A[i][i]);for(r=i+1;r<k;r++){var av=Math.abs(A[r][i]);if(av>mx){mx=av;p=r;}}
    if(mx<1e-12)return null;if(p!=i){var t=A[p];A[p]=A[i];A[i]=t;t=I[p];I[p]=I[i];I[i]=t;}
    var piv=A[i][i];for(j=0;j<k;j++){A[i][j]/=piv;I[i][j]/=piv;}
    for(r=0;r<k;r++){if(r==i)continue;var f=A[r][i];if(f===0)continue;for(j=0;j<k;j++){A[r][j]-=f*A[i][j];I[r][j]-=f*I[i][j];}}}
  return I;}
function robustCovInv(C,k){var tr=0,i;for(i=0;i<k;i++)tr+=C[i][i];var base=(tr/k)||1;var fac=[1e-3,1e-2,1e-1,1,10],t;
  for(t=0;t<fac.length;t++){var Iv=invRidge(C,k,fac[t]*base);if(Iv)return {inv:Iv,lambda:fac[t]*base};}return null;}
function chi2q975(k){var zp=1.959964,a=2/(9*k),tt=1-a+zp*Math.sqrt(a);return k*tt*tt*tt;}   // Wilson-Hilferty; nominal reference only
function mahalanobis(pool){var B=buildZ(pool);if(B.k<2||B.n<3)return {ok:false,B:B};
  var C=covMatrix(B.Z,B.k),R2=robustCovInv(C,B.k);if(!R2)return {ok:false,B:B};var Iv=R2.inv,rows=[],r;
  for(r=0;r<B.n;r++){var z=B.Z[r],acc=0,i,j;for(i=0;i<B.k;i++){var Iz=0;for(j=0;j<B.k;j++)Iz+=Iv[i][j]*z[j];acc+=z[i]*Iz;}if(acc<0)acc=0;
    var drv=B.keys.map(function(kk,ki){return {k:kk,label:(MET[kk]||{}).label||kk,z:z[ki]};}).sort(function(a,b){return Math.abs(b.z)-Math.abs(a.z);}).slice(0,3);
    rows.push({s:B.samples[r],d2:acc,drivers:drv});}
  return {ok:true,B:B,rows:rows};}
function spreadEllipse(pts){if(pts.length<4)return null;var n=pts.length,mx=0,my=0,i;   // 2-SD spread ellipse in pixel space
  for(i=0;i<n;i++){mx+=pts[i][0];my+=pts[i][1];}mx/=n;my/=n;
  var sxx=0,syy=0,sxy=0;for(i=0;i<n;i++){var dx=pts[i][0]-mx,dy=pts[i][1]-my;sxx+=dx*dx;syy+=dy*dy;sxy+=dx*dy;}
  var d=(n-1)||1;sxx/=d;syy/=d;sxy/=d;var tr=sxx+syy,det=sxx*syy-sxy*sxy,disc=Math.sqrt(Math.max(0,tr*tr/4-det));
  var l1=tr/2+disc,l2=tr/2-disc;if(l1<0)l1=0;if(l2<0)l2=0;
  var ang=(Math.abs(sxy)<1e-12&&sxx>=syy)?0:Math.atan2(l1-sxx,sxy);
  return {cx:mx,cy:my,rx:2*Math.sqrt(l1),ry:2*Math.sqrt(l2),angle:ang};}
function yearOf(dstr){if(dstr==null)return null;var m=String(dstr).match(/(1[5-9]\d\d|20\d\d|2100)/);if(!m)return null;var y=+m[1];return (y>=1500&&y<=2100)?y:null;}
// register d^2 as a first-class metric + a WARN-only flag (never in FAILF -> never a FAIL, never auto-excluded)
R.metrics.push({key:'qc_mahal',label:'QC outlier d²',kind:'float',dir:'hi_bad'});
if(DIST.indexOf('qc_mahal')<0)DIST.push('qc_mahal');
R.defs.qc_mahal=['Robust Mahalanobis distance squared in standardized QC-metric space (ridge-regularized robust covariance, median center). A multivariate rarity rank: high = an unusual COMBINATION of QC metrics even when each metric is individually in range. A heuristic composite, not a hypothesis test; missing metrics are imputed to the cohort median so a sparse sample reads more central than it may be.','robust flag QC_OUTLIER when d2 > median + 3*MAD of the cohort d2 distribution'];
R.mahal_cut=null;R.mahal_k=null;R.mahal_ref=null;R.mahal_rows=null;
// compute d^2 over the WHOLE cohort (stable flag under filtering); called inside recompute() so it stays live.
function applyMahal(){
  var pool=R.samples.filter(function(s){return !s.anc&&(!s.f||s.f.indexOf('NO_DATA')<0);});   // modern, produced-something (aDNA is legitimately extreme; scope like the SNP-z)
  R.samples.forEach(function(s){s.m.qc_mahal=null;s._mdrv=null;});                             // default: not scored (aDNA / NO_DATA -> NA, no ring, no flag)
  var M=mahalanobis(pool);
  if(!M.ok){R.mahal_cut=null;R.mahal_k=null;R.mahal_ref=null;R.mahal_rows=null;return;}
  var d2s=M.rows.map(function(o){return o.d2;}).sort(function(a,b){return a-b;});
  var md=med2(d2s),mad=med2(d2s.map(function(v){return Math.abs(v-md);}).sort(function(a,b){return a-b;}));
  var scale=1.4826*mad;   // MAD=0 (ties from median-imputation) -> fall back to SD; genuinely zero spread -> Infinity (flag nothing), like the SNP-z guard
  if(!(scale>0)){var mn=d2s.reduce(function(a,b){return a+b;},0)/d2s.length;scale=Math.sqrt(d2s.reduce(function(a,b){return a+(b-mn)*(b-mn);},0)/d2s.length);}
  R.mahal_cut=(scale>0)?md+3*scale:Infinity;R.mahal_k=M.B.k;R.mahal_ref=chi2q975(M.B.k);R.mahal_rows=M.rows;
  var by={};M.rows.forEach(function(o){by[o.s.s]=o;});
  R.samples.forEach(function(s){var o=by[s.s];if(o){s.m.qc_mahal=o.d2;s._mdrv=o.drivers;}});
  // refresh this live-derived metric's display scaling (RANGES/MED/IQR are built once at load, before d2 exists)
  var vs=R.samples.map(function(s){return s.m.qc_mahal;}).filter(function(v){return v!=null;});
  if(vs.length){RANGES.qc_mahal=[Math.min.apply(null,vs),Math.max.apply(null,vs)];var so=vs.slice().sort(function(a,b){return a-b;});MED.qc_mahal=so[Math.floor(so.length/2)];IQR.qc_mahal=[quant(so,0.25),quant(so,0.5),quant(so,0.75)];}}
// ---- Panel 1: QC-metric-space PCA biplot (over the live cohort) ----
function renderQCspace(){var host=el('qcpca_body'),cap=el('qcpca_caption'),ot=el('pca_outtable');if(!host)return;
  var P=pca2(visible());
  if(!P.ok){host.innerHTML='<div class="nd">'+esc(P.reason)+'</div>';if(cap)cap.innerHTML='';if(ot)ot.innerHTML='';return;}
  var pad=44,W=Math.min(host.clientWidth||620,760),H=Math.min(W,420);
  var xs=P.scores.map(function(p){return p[0];}),ys=P.scores.map(function(p){return p[1];});
  var xr=[Math.min.apply(null,xs),Math.max.apply(null,xs)],yr=[Math.min.apply(null,ys),Math.max.apply(null,ys)];
  var mgx=(xr[1]-xr[0])*0.08||1,mgy=(yr[1]-yr[0])*0.08||1;xr=[xr[0]-mgx,xr[1]+mgx];yr=[yr[0]-mgy,yr[1]+mgy];
  function sx(v){return pad+(v-xr[0])/(xr[1]-xr[0]||1)*(W-pad-14);}
  function sy(v){return H-pad-(v-yr[0])/(yr[1]-yr[0]||1)*(H-pad-16);}
  var groups={};P.B.samples.forEach(function(s,i){var g=(st.colorBy=='lineage')?(s.lineage||'NA'):s.v;(groups[g]=groups[g]||[]).push([sx(P.scores[i][0]),sy(P.scores[i][1])]);});
  var ellSVG='',gn;for(gn in groups){var E=spreadEllipse(groups[gn]);if(!E)continue;var col=(st.colorBy=='lineage')?linColor(gn):(VCOL[gn]||'#94a3b8');
    ellSVG+='<ellipse cx="'+E.cx.toFixed(1)+'" cy="'+E.cy.toFixed(1)+'" rx="'+E.rx.toFixed(1)+'" ry="'+E.ry.toFixed(1)+'" transform="rotate('+(E.angle*180/Math.PI).toFixed(1)+' '+E.cx.toFixed(1)+' '+E.cy.toFixed(1)+')" fill="'+col+'" fill-opacity="0.06" stroke="'+col+'" stroke-opacity="0.4" stroke-width="1"/>';}
  var dots=P.B.samples.map(function(s,i){var big=(st.hi==s.s),cx=sx(P.scores[i][0]).toFixed(1),cy=sy(P.scores[i][1]).toFixed(1);
    var ring=(R.mahal_cut!=null&&s.m.qc_mahal!=null&&s.m.qc_mahal>R.mahal_cut)?'<circle cx="'+cx+'" cy="'+cy+'" r="7" fill="none" stroke="'+VCOL.WARN+'" stroke-width="1.3"/>':'';
    return ring+'<circle cx="'+cx+'" cy="'+cy+'" r="'+(big?5.4:3.4)+'" fill="'+dotColor(s)+'" opacity="0.85"'+((s.v=='FAIL'||big)?' stroke="#0f1c29" stroke-width="'+(big?1.4:0.6)+'"':'')+' data-s="'+esc(s.s)+'" data-lin="'+esc(s.lineage||'')+'" data-x="'+P.scores[i][0]+'" data-y="'+P.scores[i][1]+'" data-xl="PC1" data-yl="PC2" data-xk="float" data-yk="float"/>';}).join('');
  var frame='<line x1="'+pad+'" y1="'+(H-pad)+'" x2="'+(W-14)+'" y2="'+(H-pad)+'" stroke="#cbd5e1"/><line x1="'+pad+'" y1="14" x2="'+pad+'" y2="'+(H-pad)+'" stroke="#cbd5e1"/>';
  var zero='';if(xr[0]<0&&xr[1]>0)zero+='<line x1="'+sx(0).toFixed(1)+'" y1="14" x2="'+sx(0).toFixed(1)+'" y2="'+(H-pad)+'" stroke="#eef2f6"/>';if(yr[0]<0&&yr[1]>0)zero+='<line x1="'+pad+'" y1="'+sy(0).toFixed(1)+'" x2="'+(W-14)+'" y2="'+sy(0).toFixed(1)+'" stroke="#eef2f6"/>';
  var xt='<text x="'+((pad+W-14)/2)+'" y="'+(H-6)+'" text-anchor="middle" font-size="11" fill="#475569">PC1 ('+P.pev[0].toFixed(1)+'%)</text>';
  var yt='<text transform="rotate(-90 13 '+((14+H-pad)/2)+')" x="13" y="'+((14+H-pad)/2)+'" text-anchor="middle" font-size="11" fill="#475569">PC2 ('+(P.pev[1]<1e-3?'~0%, rank-deficient':P.pev[1].toFixed(1)+'%')+')</text>';
  var load='<div style="display:flex;gap:24px;flex-wrap:wrap;margin-top:8px">'+[0,1].map(function(pc){return '<div style="flex:1;min-width:170px"><div class="dsub" style="margin:2px 0 6px">PC'+(pc+1)+' loadings</div>'+P.load[pc].map(function(l){var w=Math.abs(l.w),col=l.w>=0?'#22a06b':'#e0544f';return '<div class="drow"><span class="dk">'+esc(l.label)+'</span><div class="dbarwrap"><div class="dbar" style="width:'+Math.round(w*100)+'%;background:'+col+'"></div></div><span class="dv">'+(l.w>=0?'+':'')+l.w.toFixed(2)+'</span></div>';}).join('')+'</div>';}).join('')+'</div>';
  host.innerHTML='<svg viewBox="0 0 '+W+' '+H+'" style="width:100%;height:auto;display:block;cursor:crosshair">'+ellSVG+frame+zero+xt+yt+dots+'</svg>'+load;
  if(cap)cap.innerHTML=PCA_CAPTION;
  if(ot){var uv=visible().filter(function(s){return s.m.qc_mahal!=null;}).sort(function(a,b){return b.m.qc_mahal-a.m.qc_mahal;}).slice(0,8);
    ot.innerHTML='<thead><tr><th class="s" style="text-align:left">Sample</th><th>d²</th><th>QC</th><th style="text-align:left">Top deviating metrics</th></tr></thead><tbody>'+
     (uv.length?uv.map(function(s){var drv=(s._mdrv||[]).map(function(dd){var col=dd.z>=0?'#a01f2d':'#0b7350';return '<span style="color:'+col+'">'+esc(dd.label)+' '+(dd.z>=0?'+':'')+dd.z.toFixed(1)+'σ</span>';}).join(', ');
       return '<tr data-s="'+esc(s.s)+'"'+(st.hi==s.s?' class="hl"':'')+'><td class="s"><span class="sname" data-s="'+esc(s.s)+'">'+esc(s.s)+'</span></td><td>'+shortv(s.m.qc_mahal,'float')+'</td><td><span class="v '+s.v+'">'+s.v+'</span></td><td style="text-align:left;white-space:normal">'+drv+'</td></tr>';}).join('')
      :'<tr><td colspan="4" style="text-align:left;color:#94a3b8;padding:8px">no samples in view</td></tr>')+'</tbody>';
    Array.prototype.forEach.call(ot.querySelectorAll('tbody tr[data-s]'),function(tr){tr.onclick=function(){setHi(tr.getAttribute('data-s'));};});}}
// ---- Panel 2: divergence vs completeness (reference-bias screen) ----
function renderRefBias(){var host=el('divcomp_body'),cap=el('divcomp_caption'),qn=el('divcomp_quadn');if(!host)return;
  st.divx=st.divx||'missing_pct';
  var yk=(R.snp_density_ok?'snp_density':'snps'),ylabel=(R.snp_density_ok?'SNPs per callable Mb':'SNPs vs reference');
  function xof(s){return st.divx=='callable_inv'?(s.m.callable_pct!=null?100-s.m.callable_pct:null):(s.m.missing_pct!=null?s.m.missing_pct:(s.m.callable_pct!=null?100-s.m.callable_pct:null));}
  var pool=visible().filter(function(s){return xof(s)!=null&&s.m[yk]!=null;});
  if(pool.length<3){host.innerHTML='<div class="nd">needs missing/callable and SNP data for at least 3 samples in view</div>';if(cap)cap.innerHTML='';if(qn)qn.innerHTML='';return;}
  var pad=44,W=Math.min(host.clientWidth||620,760),H=Math.min(W*0.72,360);
  var xs=pool.map(xof),ys=pool.map(function(s){return s.m[yk];});
  var xhi=Math.max.apply(null,xs)*1.06||1,yhi=Math.max.apply(null,ys)*1.06||1;
  function sx(v){return pad+v/(xhi||1)*(W-pad-14);}function sy(v){return H-pad-v/(yhi||1)*(H-pad-16);}
  var gx=thr.missing_max,divs=ys.slice().sort(function(a,b){return a-b;});
  var dmed=med2(divs),dmad=med2(divs.map(function(v){return Math.abs(v-dmed);}).sort(function(a,b){return a-b;}));
  var ylo=dmed-2.5*1.4826*dmad;if(!(dmad>0)||!(ylo>0))ylo=0.5*dmed;   // robust low cut (guard MAD=0 ties + wide spread): only ANOMALOUSLY low divergence is called out (not just below median)
  var GX=sx(Math.min(gx,xhi)),GY=sy(Math.max(0,ylo));
  var rects='<rect x="'+pad+'" y="14" width="'+Math.max(0,W-14-pad).toFixed(1)+'" height="'+Math.max(0,GY-14).toFixed(1)+'" fill="var(--pass)" opacity="0.04"/>'
    +'<rect x="'+pad+'" y="'+GY.toFixed(1)+'" width="'+Math.max(0,GX-pad).toFixed(1)+'" height="'+Math.max(0,H-pad-GY).toFixed(1)+'" fill="var(--fail)" opacity="0.07"/>'
    +'<rect x="'+GX.toFixed(1)+'" y="'+GY.toFixed(1)+'" width="'+Math.max(0,W-14-GX).toFixed(1)+'" height="'+Math.max(0,H-pad-GY).toFixed(1)+'" fill="var(--warn)" opacity="0.06"/>';
  var guides='<line x1="'+GX.toFixed(1)+'" y1="14" x2="'+GX.toFixed(1)+'" y2="'+(H-pad)+'" stroke="#94a3b8" stroke-dasharray="4 3"/><line x1="'+pad+'" y1="'+GY.toFixed(1)+'" x2="'+(W-14)+'" y2="'+GY.toFixed(1)+'" stroke="#94a3b8" stroke-dasharray="4 3"/>';
  var quad={refbias:0,lowcov:0,typical:0};
  var dots=pool.map(function(s){var xv=xof(s),yv=s.m[yk],lowMiss=xv<=gx,lowDiv=yv<ylo;
    var cls=lowDiv?(lowMiss?'refbias':'lowcov'):'typical';quad[cls]++;
    var big=(st.hi==s.s),cx=sx(xv).toFixed(1),cy=sy(yv).toFixed(1);
    var ring=(cls=='refbias')?'<circle cx="'+cx+'" cy="'+cy+'" r="6.5" fill="none" stroke="#a01f2d" stroke-width="1.2"/>':'';
    return ring+'<circle cx="'+cx+'" cy="'+cy+'" r="'+(big?5.4:3.4)+'" fill="'+dotColor(s)+'" opacity="0.82"'+((s.v=='FAIL'||big)?' stroke="#0f1c29" stroke-width="'+(big?1.4:0.6)+'"':'')+' data-s="'+esc(s.s)+'" data-lin="'+esc(s.lineage||'')+'" data-x="'+xv+'" data-y="'+yv+'" data-xl="'+(st.divx=='callable_inv'?'100 - callable %':'Missing %')+'" data-yl="'+esc(ylabel)+'" data-xk="pct" data-yk="'+(R.snp_density_ok?'float':'int')+'"/>';}).join('');
  var frame='<line x1="'+pad+'" y1="'+(H-pad)+'" x2="'+(W-14)+'" y2="'+(H-pad)+'" stroke="#cbd5e1"/><line x1="'+pad+'" y1="14" x2="'+pad+'" y2="'+(H-pad)+'" stroke="#cbd5e1"/>';
  var titles='<text x="'+((pad+W-14)/2)+'" y="'+(H-6)+'" text-anchor="middle" font-size="11" fill="#475569">'+(st.divx=='callable_inv'?'100 - callable % (incompleteness)':'Missing % (incompleteness)')+'</text><text transform="rotate(-90 13 '+((14+H-pad)/2)+')" x="13" y="'+((14+H-pad)/2)+'" text-anchor="middle" font-size="11" fill="#475569">'+esc(ylabel)+'</text>';
  var labs='<text x="'+(pad+6)+'" y="'+(H-pad-6)+'" font-size="8.5" font-weight="600" fill="var(--fail)">reference-bias suspect</text>'
    +'<text x="'+(pad+6)+'" y="24" font-size="8.5" font-weight="600" fill="#3f7d55">typical divergence</text>'
    +'<text x="'+(W-16)+'" y="'+(H-pad-6)+'" text-anchor="end" font-size="8.5" font-weight="600" fill="var(--warn)">low coverage</text>';
  host.innerHTML='<svg viewBox="0 0 '+W+' '+H+'" style="width:100%;height:auto;display:block;cursor:crosshair">'+rects+guides+frame+titles+labs+dots+'</svg>';
  if(qn)qn.innerHTML=[['refbias','reference-bias suspect','var(--fail)'],['lowcov','honest low-coverage','var(--warn)'],['typical','typical divergence','var(--pass)']].map(function(q){return '<span><i style="background:'+q[2]+'"></i>'+q[1]+' <b>'+quad[q[0]]+'</b></span>';}).join('');
  if(cap)cap.innerHTML=REFBIAS_CAPTION;
  var dxb=el('divx');if(dxb)Array.prototype.forEach.call(dxb.querySelectorAll('button'),function(b){b.onclick=function(){st.divx=b.getAttribute('data-x');Array.prototype.forEach.call(dxb.querySelectorAll('button'),function(x){x.classList.toggle('on',x==b);});renderRefBias();};});}
// ---- Panel 3: temporal sampling overview (per lineage) ----
function renderTemporal(){var host=el('temporal_body'),cap=el('temporal_caption');if(!host)return;
  var pool=visible();
  var dated=pool.map(function(s){return {s:s,y:yearOf(s.date)};}).filter(function(o){return o.y!=null;});
  if(!dated.length){host.innerHTML='<div class="anote">Temporal overview needs sample collection years (a <b>date</b> column); none parseable in this run. A time-resolved / clock analysis is a downstream step.</div>';if(cap)cap.innerHTML='';return;}
  var ymin=Math.min.apply(null,dated.map(function(o){return o.y;})),ymax=Math.max.apply(null,dated.map(function(o){return o.y;})),span=Math.max(1,ymax-ymin);
  var g={};pool.forEach(function(s){var L=R.lin_present?(s.lineage||'NA'):'all samples';(g[L]=g[L]||[]).push(s);});
  var order=R.lin_present?(R.lineages||[]).concat(['NA']):['all samples'],seen={},rows=[];
  order.forEach(function(L){if(g[L]&&!seen[L]){seen[L]=1;rows.push(L);}});Object.keys(g).forEach(function(L){if(!seen[L]){seen[L]=1;rows.push(L);}});
  var gut=140,W=Math.min(host.clientWidth||620,900),plotW=Math.max(80,W-gut-70);
  function ax(y){return gut+(y-ymin)/span*plotW;}
  var ticks='';for(var t=0;t<=4;t++){var yr=Math.round(ymin+span*t/4),X=ax(yr);ticks+='<line x1="'+X.toFixed(1)+'" y1="14" x2="'+X.toFixed(1)+'" y2="18" stroke="#cbd5e1"/><text x="'+X.toFixed(1)+'" y="11" text-anchor="middle" font-size="8.5" fill="#94a3b8">'+yr+'</text>';}
  var axisSVG='<svg viewBox="0 0 '+W+' 22" style="width:100%;height:auto;display:block"><line x1="'+gut+'" y1="18" x2="'+(gut+plotW)+'" y2="18" stroke="#e6ebf1"/>'+ticks+'</svg>';
  var body=rows.map(function(L){var grp=g[L],ys=grp.map(function(s){return yearOf(s.date);}).filter(function(y){return y!=null;});
    var snps=grp.map(function(s){return s.m.snps;}).filter(function(v){return v!=null;}).sort(function(a,b){return a-b;});
    var titv=grp.map(function(s){return s.m.ti_tv;}).filter(function(v){return v!=null;}).sort(function(a,b){return a-b;});
    var mn=ys.length?Math.min.apply(null,ys):null,mxx=ys.length?Math.max.apply(null,ys):null;
    var col=R.lin_present?linColor(L):'var(--accent)';
    var rug=grp.map(function(s,i){var y=yearOf(s.date);if(y==null)return '';var jx=(((i*2654435761)>>>0)%997)/997-0.5,big=(st.hi==s.s);
      return '<circle cx="'+ax(y).toFixed(1)+'" cy="'+(13+jx*8).toFixed(1)+'" r="'+(big?4.5:2.4)+'" fill="'+col+'" opacity="0.72"'+(big?' stroke="#0f1c29" stroke-width="1"':'')+' data-s="'+esc(s.s)+'" data-lin="'+esc(s.lineage||'')+'" data-x="'+y+'" data-y="'+(s.m.snps!=null?s.m.snps:0)+'" data-xl="year" data-yl="SNPs" data-xk="int" data-yk="int"/>';}).join('');
    var bar=(mn!=null)?'<rect x="'+ax(mn).toFixed(1)+'" y="10" width="'+Math.max(2,ax(mxx)-ax(mn)).toFixed(1)+'" height="6" rx="3" fill="'+col+'" opacity="0.18"/>':'';
    var lineSVG='<svg viewBox="0 0 '+W+' 26" style="width:100%;height:auto;display:block">'+bar+rug+'</svg>';
    var proxy=snps.length?med2(snps):null,mtitv=titv.length?med2(titv):null;
    return '<div class="bee" style="height:auto;min-height:34px"><div class="bl"><span class="ldot" style="background:'+col+'"></span>'+esc(L)+'<span class="stat">n='+grp.length+(mn!=null?' · '+mn+'-'+mxx:' · no dates')+(mtitv!=null?' · Ti/Tv '+mtitv.toFixed(2):'')+'</span></div><div class="plotwrap">'+lineSVG+'</div><div class="sv" style="flex:0 0 62px" title="informative-site proxy (median SNPs)">'+(proxy!=null?Math.round(proxy).toLocaleString('en-US'):'-')+'</div></div>';}).join('');
  host.innerHTML=axisSVG+body;
  var cohTitv=(function(){var tv=pool.map(function(s){return s.m.ti_tv;}).filter(function(v){return v!=null;}).sort(function(a,b){return a-b;});return tv.length?med2(tv):null;})();
  if(cap)cap.innerHTML=TEMPORAL_CAPTION.replace('{DATED}',dated.length).replace('{TOT}',pool.length).replace('{SPAN}',ymin+' - '+ymax+' ('+span+' yr)')+(cohTitv!=null?' Cohort median Ti/Tv = '+cohTitv.toFixed(2)+' (a ratio); a ratio drifting toward parity indicates mutational saturation where the SNP proxy over-counts, and a high ratio is spectrum or damage, not saturation.':'');}
// ===============================================================================
// ---- auto-discovered extra metrics: merge into the registry, hidden by default ----
var extraSet={}; R.extra.forEach(function(e){extraSet[e.key]=1;R.metrics.push(e);R.defs[e.key]=R.defs[e.key]||['Auto-detected metric from the summary TSV (not a named QC metric).',''];});
var st={sortKey:'s',asc:true,q:'',onlyFlagged:false,hidden:{},hi:null,flagFilter:null,ptype:'beeswarm',
        sx:'mean_depth',sy:'breadth_pct',excl:{},detail:null,colorBy:'qc',groupLin:false,ancOnly:null,linFilter:null,gtrack:'missing',maskOn:false,gsel:null,gzoom:null,gbq:'',hotq:'',pnpsq:''};
var SGEO=null, GGEO=null;   // scatter + genome brush geometry caches (for inverse-mapping the rubber-band)
R.extra.forEach(function(e){st.hidden[e.key]=1;});
// narrow (embedded panel / mobile): show only the essentials by default; the long tail is one click away in "columns"
if(window.innerWidth<860){['raw_reads','trimmed_reads','read_length','q20_pct','q30_pct','gc_pct','properly_paired_pct','avg_base_quality','het_variants','homo_indels','iupac_pct','het_frac','n_lineages','est_genome_cov','total_variants','median_depth','breadth5x_pct','coverage_cv','unmapped_pct','singleton_pct','mapq0_pct','error_rate_pct','insert_size','insert_size_sd','longest_n_run','n_gaps','ann_moderate','ann_low','ann_modifier','coding_pct','lof_pct','missense_silent','snpeff_warn'].forEach(function(k){st.hidden[k]=1;});}
var MET={}; R.metrics.forEach(function(m){MET[m.key]=m;});
function quant(so,p){if(!so.length)return null;var i=(so.length-1)*p,lo=Math.floor(i),hi=Math.ceil(i);return lo==hi?so[lo]:so[lo]+(so[hi]-so[lo])*(i-lo);}
var RANGES={},MED={},IQR={}; R.metrics.forEach(function(m){
  var vs=R.samples.map(function(s){return s.m[m.key];}).filter(function(v){return v!=null;});
  RANGES[m.key]=vs.length?[Math.min.apply(null,vs),Math.max.apply(null,vs)]:[0,1];
  var so=vs.slice().sort(function(a,b){return a-b;}); MED[m.key]=so.length?so[Math.floor(so.length/2)]:null;
  IQR[m.key]=so.length?[quant(so,0.25),quant(so,0.5),quant(so,0.75)]:null;
});
// robust snp median/mad over MODERN samples (matches the Python cohort scoping)
// even-n midpoint average to match the Python robust() convention (so qc_flags.tsv and the HTML agree)
function med2(a){var n=a.length;return n%2?a[(n-1)/2]:(a[n/2-1]+a[n/2])/2;}
var SNPMED=null,SNPSIG=null;(function(){var vs=R.samples.filter(function(s){return !s.anc;}).map(function(s){return s.m.snps;}).filter(function(v){return v!=null;}).sort(function(a,b){return a-b;});
  if(!vs.length)return; var n=vs.length,med=med2(vs);var d=vs.map(function(v){return Math.abs(v-med);}).sort(function(a,b){return a-b;});
  var mad=med2(d),sig=1.4826*mad; if(sig<=0){var mean=vs.reduce(function(a,b){return a+b;},0)/n; sig=Math.sqrt(vs.reduce(function(a,b){return a+(b-mean)*(b-mean);},0)/n)||1;} SNPMED=med;SNPSIG=sig;})();
// ---- snpEff functional-annotation cohort fits (over NON-ANCIENT samples; robust median + MAD like the SNP-z) ----
function robustMedSig(vals){var v=vals.filter(function(x){return x!=null;}).sort(function(a,b){return a-b;});
  if(v.length<4)return null;var md=med2(v);
  var d=v.map(function(x){return Math.abs(x-md);}).sort(function(a,b){return a-b;});var sig=1.4826*med2(d);
  if(sig<=0){var mn=v.reduce(function(a,b){return a+b;},0)/v.length;sig=Math.sqrt(v.reduce(function(a,b){return a+(b-mn)*(b-mn);},0)/v.length)||0;}
  return sig>0?{med:md,sig:sig}:null;}
var ANNHI=null,LOFHI=null,ANNPCT=null,MSIL=null;
(function(){var mod=R.samples.filter(function(s){return !s.anc;});
  ANNHI =robustMedSig(mod.map(function(s){return s.m.ann_high;}));
  LOFHI =robustMedSig(mod.map(function(s){return s.m.lof_pct;}));
  ANNPCT=robustMedSig(mod.map(function(s){return s.m.annotated_pct;}));
  MSIL  =robustMedSig(mod.map(function(s){return s.m.missense_silent;}));})();
var ANNOT_FLOOR=(typeof thr.annot_floor=='number')?thr.annot_floor:50;   // absolute annotated-% sanity floor
var FAILF={LOW_DEPTH:1,LOW_BREADTH:1,HIGH_MISSING:1,NO_DATA:1};
// ---- expected-range bands: acceptable [lo,hi] per metric, derived LIVE from the ACTIVE threshold set ----
var TITV_WINDOW=[1.5,2.1];  // Ti/Tv sanity window (a ratio); A4/M. bovis legitimately runs higher (~2.6)
function bandFor(pk,T){switch(pk){
  case 'mean_depth':return [T.depth_min,Infinity];
  case 'breadth_pct':return [T.breadth_min,Infinity];
  case 'mapped_pct':return [T.mapping_min,Infinity];
  case 'missing_pct':return [-Infinity,T.missing_max];
  case 'duplication_pct':return [-Infinity,T.dup_max];
  case 'iupac_pct':return [-Infinity,T.iupac_max];
  case 'het_frac':return [-Infinity,thr.het_max_frac];
  case 'ti_tv':return [(thr.titv_min>0?Math.max(thr.titv_min,TITV_WINDOW[0]):TITV_WINDOW[0]),TITV_WINDOW[1]];
  case 'snps':return (SNPMED==null||!SNPSIG)?null:[SNPMED-thr.snp_z*SNPSIG,SNPMED+thr.snp_z*SNPSIG];
  case 'ann_high':return (ANNHI==null)?null:[-Infinity,ANNHI.med+thr.snp_z*ANNHI.sig];
  case 'annotated_pct':return (ANNPCT==null)?null:[ANNPCT.med-thr.snp_z*ANNPCT.sig,Infinity];
  case 'missense_silent':return null;   // deliberately no band: it is a proxy, high-tail only
  default:return null;}}
var RULES=[['LOW_DEPTH','mean_depth','<','depth_min','FAIL'],['LOW_BREADTH','breadth_pct','<','breadth_min','FAIL'],
 ['HIGH_MISSING','missing_pct','>','missing_max','FAIL'],['MAPPING_LOW','mapped_pct','<','mapping_min','WARN'],
 ['HIGH_DUP','duplication_pct','>','dup_max','WARN'],['HIGH_IUPAC','iupac_pct','>','iupac_max','WARN'],
 ['TITV_LOW','ti_tv','<','titv_min','WARN']];
function recompute(){
  R.samples.forEach(function(s){
    if(['mean_depth','breadth_pct','missing_pct','mapped_pct'].every(function(k){return s.m[k]==null;})){s.f=['NO_DATA'];s.v='FAIL';return;}
    var f=[],T=actv(s);
    RULES.forEach(function(r){var v=s.m[r[1]];if(v==null)return;
      var t=(r[0]=='TITV_LOW')?thr[r[3]]:(T[r[3]]!=null?T[r[3]]:thr[r[3]]);
      if(r[0]=='TITV_LOW'&&!(t>0))return; if(r[2]=='<'?v<t:v>t)f.push(r[0]);});
    if(!s.anc){var snp=s.m.snps; if(snp!=null&&SNPMED!=null){if(snp<SNPMED-thr.snp_z*SNPSIG)f.push('SNP_LOW');else if(snp>SNPMED+thr.snp_z*SNPSIG)f.push('SNP_HIGH');}
      if(s.m.het_frac!=null&&s.m.het_frac>thr.het_max_frac)f.push('HET_HIGH');
      // snpEff functional-annotation flags (WARN-only, cohort-relative robust-z; aDNA-guarded: deamination inflates missense/stop)
      if(ANNHI&&s.m.ann_high!=null&&s.m.ann_high>ANNHI.med+thr.snp_z*ANNHI.sig)f.push('HIGH_IMPACT_EXCESS');
      if(LOFHI&&s.m.lof_pct!=null&&s.m.lof_pct>LOFHI.med+thr.snp_z*LOFHI.sig)f.push('LOF_EXCESS');
      if(MSIL&&s.m.missense_silent!=null&&s.m.missense_silent>MSIL.med+thr.snp_z*MSIL.sig)f.push('PNPS_PROXY_HIGH');
      if(s.m.annotated_pct!=null&&(s.m.annotated_pct<ANNOT_FLOOR||(ANNPCT&&s.m.annotated_pct<90&&s.m.annotated_pct<ANNPCT.med-thr.snp_z*ANNPCT.sig)))f.push('ANNOTATION_POOR');}
    s.m.n_lineages=nlin(s); if(s.m.n_lineages!=null&&s.m.n_lineages>1)f.push('MIXED');
    if(s.anc&&s.dmg&&s.dmg.ct1!=null&&athr.damage_min_ct!=null&&s.dmg.ct1<athr.damage_min_ct)f.push('DAMAGE_LOW');
    s.f=f; s.v=f.some(function(x){return FAILF[x];})?'FAIL':(f.length?'WARN':'PASS');});
  // multivariate QC-outlier flag (cohort-wide fit; WARN-only, never FAIL, never auto-excluded)
  applyMahal();
  if(R.mahal_cut!=null)R.samples.forEach(function(s){
    if(s.f.indexOf('NO_DATA')>=0)return;                                  // a produced-nothing FAIL is not a "combination outlier"
    if(s.m.qc_mahal!=null&&s.m.qc_mahal>R.mahal_cut){if(s.f.indexOf('QC_OUTLIER')<0)s.f.push('QC_OUTLIER');if(s.v=='PASS')s.v='WARN';}});
  R.counts={PASS:0,WARN:0,FAIL:0}; R.samples.forEach(function(s){R.counts[s.v]++;});
}
// flag margin ("why"): human-readable reason for one flag on a sample (for the modal + exclusion reason)
function flagWhy(s,fl){var T=actv(s),M=MET;function u(k){return (M[k]&&M[k].kind=='pct')?'%':'';}
  switch(fl){
    case 'LOW_DEPTH':return 'depth '+shortv(s.m.mean_depth,'float')+' < '+T.depth_min;
    case 'LOW_BREADTH':return 'breadth '+shortv(s.m.breadth_pct,'pct')+'% < '+T.breadth_min+'%';
    case 'HIGH_MISSING':return 'missing '+shortv(s.m.missing_pct,'pct')+'% > '+T.missing_max+'%';
    case 'MAPPING_LOW':return 'mapped '+shortv(s.m.mapped_pct,'pct')+'% < '+T.mapping_min+'%';
    case 'HIGH_DUP':return 'dup '+shortv(s.m.duplication_pct,'pct')+'% > '+T.dup_max+'%';
    case 'HIGH_IUPAC':return 'IUPAC '+shortv(s.m.iupac_pct,'pct')+'% > '+T.iupac_max+'%';
    case 'TITV_LOW':return 'Ti/Tv '+shortv(s.m.ti_tv,'float')+' < '+thr.titv_min;
    case 'SNP_LOW':return 'SNPs '+shortv(s.m.snps,'int')+' below cohort (z<-'+thr.snp_z+')';
    case 'SNP_HIGH':return 'SNPs '+shortv(s.m.snps,'int')+' above cohort (z>'+thr.snp_z+')';
    case 'HET_HIGH':return 'het fraction '+shortv(s.m.het_frac,'pct')+'% > '+thr.het_max_frac+'%';
    case 'MIXED':return s.m.n_lineages+' lineages above '+thr.mixed_min_frac+'%';
    case 'DAMAGE_LOW':return "5' C>T "+(s.dmg&&s.dmg.ct1!=null?(s.dmg.ct1*100).toFixed(1)+'%':'NA')+' < '+(athr.damage_min_ct*100).toFixed(0)+'% (aDNA auth)';
    case 'NO_DATA':return 'no consensus / mapping output produced';
    case 'QC_OUTLIER':return 'unusual combination of QC metrics - Mahalanobis d² '+shortv(s.m.qc_mahal,'float')+' > cohort robust cut '+shortv(R.mahal_cut,'float')+(s._mdrv?' (driven by '+s._mdrv.map(function(d){return d.label+' '+(d.z>=0?'+':'')+d.z.toFixed(1)+'σ';}).join(', ')+')':'')+'; WARN only, no single metric fails';
    case 'HIGH_IMPACT_EXCESS':return 'HIGH-impact variants '+shortv(s.m.ann_high,'int')+' above cohort (robust z > '+thr.snp_z+'); possible contamination, misassembly, or a wrong reference annotation, or a genuinely reduced/divergent genome, not confirmed';
    case 'LOF_EXCESS':return 'loss-of-function fraction '+shortv(s.m.lof_pct,'pct')+'% above cohort (robust z > '+thr.snp_z+'); possible frameshift storm from indel-calling / assembly / annotation, or real pseudogenisation, not confirmed';
    case 'PNPS_PROXY_HIGH':return 'missense/silent ratio '+shortv(s.m.missense_silent,'float')+' above cohort (robust z > '+thr.snp_z+'); a spectrum proxy only, possible base-call error or contamination; not dN/dS, see the pN/pS panel';
    case 'ANNOTATION_POOR':return 'annotated '+shortv(s.m.annotated_pct,'pct')+'% below floor '+ANNOT_FLOOR+' or below cohort; the GFF-to-snpEff database may be mismatched for this reference';
    default:return fl;}}
function esc(s){return String(s).replace(/[&<>"]/g,function(c){return{'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;'}[c];});}
function fmt(v,k){if(v==null)return'NA';if(k=='int')return Math.round(v).toLocaleString('en-US');if(k=='pct')return v.toFixed(1);return v.toFixed(2);}
function shortv(v,k){if(v==null)return'';if(k=='int')return Math.round(v).toLocaleString('en-US');return (+v).toPrecision(3);}
function el(id){return document.getElementById(id);}
function visible(){return R.samples.filter(function(s){
  if(st.onlyFlagged&&s.v=='PASS')return false;
  if(st.flagFilter&&s.f.indexOf(st.flagFilter)<0)return false;
  if(st.ancOnly=='anc'&&!s.anc)return false;
  if(st.ancOnly=='mod'&&s.anc)return false;
  if(st.linFilter&&s.lineage!=st.linFilter)return false;
  if(st.q&&s.s.toLowerCase().indexOf(st.q)<0)return false; return true;});}
function dotColor(s){return st.colorBy=='lineage'?linColor(s.lineage):VCOL[s.v];}
// genome landscape: missing-fraction (0 callable -> 1 missing) mapped to a pale->red heat colour
function heatCol(mv){if(mv==null)return '#e9edf2';var a=[238,244,240],b=[214,64,58];
  return 'rgb('+Math.round(a[0]+(b[0]-a[0])*mv)+','+Math.round(a[1]+(b[1]-a[1])*mv)+','+Math.round(a[2]+(b[2]-a[2])*mv)+')';}
function fmtpos(p){return p>=1e6?(p/1e6).toFixed(2)+' Mb':p>=1e3?Math.round(p/1e3)+' kb':(''+p)+' bp';}
function tip(h,x,y){var t=el('tt'); if(!h){t.style.opacity=0;return;} t.innerHTML=h;
  var w=window.innerWidth; t.style.left=Math.min(x+13,w-t.offsetWidth-10)+'px'; t.style.top=(y+13)+'px'; t.style.opacity=1;}
function setHi(s){st.hi=(st.hi==s?null:s); renderAll();}

function donut(c){var t=(c.PASS+c.WARN+c.FAIL)||1,R0=38,C=2*Math.PI*R0,off=0,segs='';
  segs+='<circle cx="46" cy="46" r="'+R0+'" fill="none" stroke="#edf1f6" stroke-width="13"/>';
  [['PASS',VFILL.PASS],['WARN',VFILL.WARN],['FAIL',VFILL.FAIL]].forEach(function(p){var frac=c[p[0]]/t,len=frac*C;
    if(len<=0)return;
    segs+='<circle cx="46" cy="46" r="'+R0+'" fill="none" stroke="'+p[1]+'" stroke-width="13" stroke-dasharray="'+len.toFixed(2)+' '+(C-len).toFixed(2)+'" stroke-dashoffset="'+(-off).toFixed(2)+'" transform="rotate(-90 46 46)"/>'; off+=len;});
  var pct=Math.round(100*c.PASS/t);
  return '<svg width="92" height="92" viewBox="0 0 92 92">'+segs+'<text x="46" y="42" text-anchor="middle" font-size="22" font-weight="700" fill="#15202e" letter-spacing="-.5">'+pct+'%</text><text x="46" y="58" text-anchor="middle" font-size="9.5" fill="#8895a6" letter-spacing=".1em">PASS</text></svg>';}

function renderOverview(){
  var c=R.counts;
  el('summary').innerHTML=donut(c)+'<div class="counts">'+
    '<div class="c all"><div class="n">'+R.samples.length+'</div><div class="l">samples</div></div>'+
    ['PASS','WARN','FAIL'].map(function(v){return '<div class="c '+v.toLowerCase()+'"><div class="n">'+c[v]+'</div><div class="l">'+v.toLowerCase()+'</div></div>';}).join('')+'</div>';
  var freq={}; R.samples.forEach(function(s){s.f.forEach(function(f){freq[f]=(freq[f]||0)+1;});});
  var keys=Object.keys(freq).sort(function(a,b){return freq[b]-freq[a];});
  el('chips').innerHTML='<span class="t">flags</span>'+(keys.length?keys.map(function(f){return '<span class="chip'+(st.flagFilter==f?' on':'')+'" data-f="'+f+'">'+f+'<span class="k">'+freq[f]+'</span></span>';}).join(''):'<span style="color:#94a3b8;font-size:12px">none - every sample clear ✓</span>');
  Array.prototype.forEach.call(document.querySelectorAll('#chips .chip'),function(ch){ch.onclick=function(){var f=ch.getAttribute('data-f');st.flagFilter=(st.flagFilter==f?null:f);st.onlyFlagged=false;renderAll();};});
}

function renderTable(){
  var mets=R.metrics.filter(function(m){return !st.hidden[m.key];});
  var head='<tr><th class="s" data-k="s"><input type="checkbox" id="cbAll" title="exclude all shown samples"><span class="hlab"> Sample</span></th><th data-k="v">QC</th>'+
    mets.map(function(m){var d=(R.defs[m.key]||[''])[0];return '<th data-k="'+m.key+'" title="'+esc(d)+'">'+esc(m.label)+(st.sortKey==m.key?(st.asc?' ▲':' ▼'):'')+'</th>';}).join('')+
    '<th data-k="lineage" style="text-align:left">Lineage</th></tr>';
  var linRank={}; (R.lineages||[]).forEach(function(l,i){linRank[l]=i;});
  function lr(s){return (s.lineage&&linRank[s.lineage]!=null)?linRank[s.lineage]:9999;}
  var rows=visible().slice().sort(function(a,b){
    if(st.groupLin){var la=lr(a),lb=lr(b);if(la!=lb)return la-lb;}
    var k=st.sortKey,x=(k=='s')?a.s:(k=='v'?a.v:a.m[k]),y=(k=='s')?b.s:(k=='v'?b.v:b.m[k]),c;
    if(typeof x=='number'&&typeof y=='number')c=x-y;else c=String(x==null?'':x).localeCompare(String(y==null?'':y));return st.asc?c:-c;});
  var lastLin=null, ncol=mets.length+3;
  var body=rows.map(function(s){
    var pre='';
    if(st.groupLin){var lk=s.lineage||'NA'; if(lk!==lastLin){lastLin=lk;
      pre='<tr class="lingrp"><td class="s" colspan="'+ncol+'" style="text-align:left"><span class="ldot" style="background:'+linColor(s.lineage)+'"></span>'+esc(lk)+'</td></tr>';}}
    var badge=s.anc?'<span class="abadge" title="ancient (aDNA) sample">aDNA</span>':'';
    var tds='<td class="s" data-s="'+esc(s.s)+'"><input type="checkbox" class="cbx" data-s="'+esc(s.s)+'"'+(st.excl[s.s]?' checked':'')+'><span class="sname" data-s="'+esc(s.s)+'">'+esc(s.s)+'</span>'+badge+'</td><td data-v="'+s.v+'"><span class="v '+s.v+'">'+s.v+'</span></td>';
    mets.forEach(function(m){var v=s.m[m.key];
      if(v==null){tds+='<td class="na" data-v="">NA</td>';return;}
      var r=RANGES[m.key],nn=r[1]>r[0]?(v-r[0])/(r[1]-r[0]):0;nn=Math.max(0,Math.min(1,nn));var p=(nn*100).toFixed(1);
      tds+='<td data-v="'+v+'" style="background:linear-gradient(90deg,'+BAR[m.dir]+'2b 0 '+p+'%,#0000 '+p+'%)">'+fmt(v,m.kind)+'</td>';});
    tds+='<td data-v="'+esc(s.lineage||'')+'" style="text-align:left">'+esc(s.lineage||'NA')+'</td>';
    return pre+'<tr class="'+(st.hi==s.s?'hl':'')+'" data-s="'+esc(s.s)+'">'+tds+'</tr>';}).join('');
  var t=el('gstable'); t.innerHTML='<thead>'+head+'</thead><tbody>'+body+'</tbody>';
  el('nshown').textContent=rows.length+' / '+R.samples.length+' shown';
  Array.prototype.forEach.call(t.querySelectorAll('th'),function(th){th.onclick=function(){var k=th.getAttribute('data-k');if(st.sortKey==k)st.asc=!st.asc;else{st.sortKey=k;st.asc=(k=='s');}renderTable();};});
  Array.prototype.forEach.call(t.querySelectorAll('tbody tr[data-s]'),function(tr){tr.onclick=function(){setHi(tr.getAttribute('data-s'));};});
  // curation basket: exclusion checkboxes (stopPropagation so they don't sort/highlight)
  var cbAll=el('cbAll');
  if(cbAll){var vis=visible();cbAll.checked=vis.length>0&&vis.every(function(s){return st.excl[s.s];});
    cbAll.onclick=function(e){e.stopPropagation();var on=cbAll.checked;visible().forEach(function(s){if(on)st.excl[s.s]=1;else delete st.excl[s.s];});renderTable();renderCuration();};}
  Array.prototype.forEach.call(t.querySelectorAll('.cbx'),function(cb){
    cb.onclick=function(e){e.stopPropagation();};
    cb.onchange=function(){var s=cb.getAttribute('data-s');if(cb.checked)st.excl[s]=1;else delete st.excl[s];renderCuration();
      var cba=el('cbAll');if(cba){var vv=visible();cba.checked=vv.length>0&&vv.every(function(x){return st.excl[x.s];});}};});
  Array.prototype.forEach.call(t.querySelectorAll('.sname'),function(sp){sp.onclick=function(e){e.stopPropagation();openDetail(sp.getAttribute('data-s'));};});
}

function renderPlots(){
  var host=el('plots'); host.innerHTML='';
  var W=host.clientWidth||900, labelW=156, svgW=Math.max(240,W-labelW-4), padL=8, rightPad=64, pw=svgW-padL-rightPad, H=30, cy=H/2;
  var vis={}; visible().forEach(function(s){vis[s.s]=1;});
  DIST.forEach(function(pk){
    var mt=MET[pk]||{label:pk,kind:'float'};
    var rows=R.samples.filter(function(s){return s.m[pk]!=null;});
    var row=document.createElement('div'); row.className='bee';
    if(!rows.length){row.innerHTML='<span class="bl">'+esc(mt.label)+'</span><span class="nd">no data</span>';host.appendChild(row);return;}
    var lo=Math.min.apply(null,rows.map(function(s){return s.m[pk];})),hi=Math.max.apply(null,rows.map(function(s){return s.m[pk];}));
    if(hi<=lo)hi=lo+(Math.abs(lo)||1)*1e-3+1e-9;
    // acceptable-range band (moves with live thresholds); drawn behind marks. bar mode = ranked X, no band.
    var bandSVG='';(function(){var b=bandFor(pk,thr);if(!b||st.ptype=='bar')return;var blo=b[0],bhi=b[1];
      var x0=padL+(Math.max(blo,lo)-lo)/(hi-lo)*pw, x1=padL+(Math.min(bhi,hi)-lo)/(hi-lo)*pw, xw=x1-x0;
      if(xw>0.5){bandSVG='<rect x="'+x0.toFixed(1)+'" y="2" width="'+xw.toFixed(1)+'" height="'+(H-4)+'" fill="var(--bandfill)"/>';
        if(blo>lo)bandSVG+='<line x1="'+x0.toFixed(1)+'" y1="1" x2="'+x0.toFixed(1)+'" y2="'+(H-1)+'" stroke="var(--bandedge)" stroke-width="1"/>';
        if(bhi<hi)bandSVG+='<line x1="'+x1.toFixed(1)+'" y1="1" x2="'+x1.toFixed(1)+'" y2="'+(H-1)+'" stroke="var(--bandedge)" stroke-width="1"/>';}})();
    var inner='';
    if(st.ptype=='beeswarm'){
      inner=rows.map(function(s,i){var x=padL+(s.m[pk]-lo)/(hi-lo)*pw,j=((i*2654435761)%997)/997-0.5,y=cy+j*(H-9),
        big=(st.hi==s.s),dim=(st.q||st.onlyFlagged||st.flagFilter||st.ancOnly||st.linFilter)&&!vis[s.s],rr=big?4.7:(s.v!='PASS'?3.1:2.3),op=dim?0.1:(s.v!='PASS'?0.95:0.5),
        stk=(s.v=='FAIL'||big)?' stroke="#0f1c29" stroke-width="'+(big?1.3:0.6)+'"':'';
        return '<circle cx="'+x.toFixed(1)+'" cy="'+y.toFixed(1)+'" r="'+rr+'" fill="'+dotColor(s)+'" opacity="'+op+'"'+stk+' data-s="'+esc(s.s)+'" data-lin="'+esc(s.lineage||'')+'" data-pk="'+esc(pk)+'" data-val="'+s.m[pk]+'" data-lab="'+esc(mt.label)+'" data-kind="'+mt.kind+'"/>';}).join('');
      var mx=padL+(MED[pk]-lo)/(hi-lo)*pw;
      inner=bandSVG+'<line x1="'+padL+'" y1="'+cy+'" x2="'+(padL+pw)+'" y2="'+cy+'" stroke="#e6ebf1"/><line x1="'+mx.toFixed(1)+'" y1="4" x2="'+mx.toFixed(1)+'" y2="'+(H-4)+'" stroke="#64748b" stroke-dasharray="2 2"/>'+inner;
    }else if(st.ptype=='bar'){
      var sr=rows.slice().sort(function(a,b){return b.m[pk]-a.m[pk];}); var bw=pw/sr.length;
      inner=sr.map(function(s,i){var h=(s.m[pk]-Math.min(lo,0))/(hi-Math.min(lo,0))*(H-6),x=padL+i*bw,dim=(st.q||st.onlyFlagged||st.flagFilter||st.ancOnly||st.linFilter)&&!vis[s.s],big=(st.hi==s.s);
        return '<rect class="hit" x="'+x.toFixed(1)+'" y="'+(H-3-h).toFixed(1)+'" width="'+Math.max(bw-0.5,0.6).toFixed(1)+'" height="'+Math.max(h,0.5).toFixed(1)+'" fill="'+(big?'#0f1c29':dotColor(s))+'" opacity="'+(dim?0.12:(s.v!='PASS'?0.95:0.62))+'" data-s="'+esc(s.s)+'" data-lin="'+esc(s.lineage||'')+'" data-pk="'+esc(pk)+'" data-val="'+s.m[pk]+'" data-lab="'+esc(mt.label)+'" data-kind="'+mt.kind+'"/>';}).join('');
    }else{ // histogram
      var nb=Math.min(30,Math.max(8,Math.round(Math.sqrt(rows.length)))),cnt=new Array(nb).fill(0);
      rows.forEach(function(s){var b=Math.floor((s.m[pk]-lo)/(hi-lo)*nb);if(b>=nb)b=nb-1;if(b<0)b=0;cnt[b]++;});
      var cm=Math.max.apply(null,cnt)||1,bw=pw/nb;
      inner=bandSVG+cnt.map(function(c,i){var h=c/cm*(H-6),x=padL+i*bw;return '<rect x="'+x.toFixed(1)+'" y="'+(H-3-h).toFixed(1)+'" width="'+Math.max(bw-1,0.6).toFixed(1)+'" height="'+Math.max(h,0.4).toFixed(1)+'" fill="'+BAR[mt.dir]+'" opacity="0.8"/>';}).join('');
    }
    var iq=IQR[pk],statTxt=iq?('med '+shortv(iq[1],mt.kind)+' · IQR '+shortv(iq[0],mt.kind)+' to '+shortv(iq[2],mt.kind)):'';
    row.innerHTML='<span class="bl">'+esc(mt.label)+(statTxt?'<span class="stat">'+esc(statTxt)+'</span>':'')+'</span><div class="plotwrap"><svg width="'+svgW+'" height="'+H+'">'+inner+
      '<text x="'+(padL+pw+6)+'" y="'+(cy+3)+'" font-size="9" fill="#94a3b8">'+shortv(hi,mt.kind)+'</text></svg></div>';
    host.appendChild(row);
  });
}

function renderScatter(){
  var host=el('scatter'); var W=Math.min(host.classList.contains('expanded')?880:560,(host.clientWidth||520)); var S=Math.max(300,W); var pad=42, plot=S-pad-14, H=(host.classList.contains('expanded')?Math.min(660,S):340), ph=H-pad-14;
  var xk=st.sx,yk=st.sy,xm=MET[xk],ym=MET[yk];
  var rows=R.samples.filter(function(s){return s.m[xk]!=null&&s.m[yk]!=null;});
  if(!rows.length){host.innerHTML='<div class="pad nd">no data for these axes</div>';return;}
  var xr=RANGES[xk],yr=RANGES[yk];
  function sx(v){return pad+(xr[1]>xr[0]?(v-xr[0])/(xr[1]-xr[0]):0.5)*plot;}
  function sy(v){return H-pad-(yr[1]>yr[0]?(v-yr[0])/(yr[1]-yr[0]):0.5)*ph;}
  var vis={}; visible().forEach(function(s){vis[s.s]=1;});
  var dots=rows.map(function(s){var big=(st.hi==s.s),dim=(st.q||st.onlyFlagged||st.flagFilter||st.ancOnly||st.linFilter)&&!vis[s.s];
    return '<circle cx="'+sx(s.m[xk]).toFixed(1)+'" cy="'+sy(s.m[yk]).toFixed(1)+'" r="'+(big?5.4:3.4)+'" fill="'+dotColor(s)+'" opacity="'+(dim?0.12:0.82)+'"'+((s.v=='FAIL'||big)?' stroke="#0f1c29" stroke-width="'+(big?1.4:0.6)+'"':'')+' data-s="'+esc(s.s)+'" data-lin="'+esc(s.lineage||'')+'" data-x="'+s.m[xk]+'" data-y="'+s.m[yk]+'" data-xl="'+esc(xm.label)+'" data-yl="'+esc(ym.label)+'" data-xk="'+xm.kind+'" data-yk="'+ym.kind+'"/>';}).join('');
  var ticks='';[0,0.5,1].forEach(function(t){var gx=pad+t*plot,gy=H-pad-t*ph;
    ticks+='<line x1="'+gx+'" y1="'+pad+'" x2="'+gx+'" y2="'+(H-pad)+'" stroke="#f0f3f7"/><line x1="'+pad+'" y1="'+gy+'" x2="'+(pad+plot)+'" y2="'+gy+'" stroke="#f0f3f7"/>'+
    '<text x="'+gx+'" y="'+(H-pad+13)+'" font-size="9" fill="#94a3b8" text-anchor="middle">'+shortv(xr[0]+t*(xr[1]-xr[0]),xm.kind)+'</text>'+
    '<text x="'+(pad-6)+'" y="'+(gy+3)+'" font-size="9" fill="#94a3b8" text-anchor="end">'+shortv(yr[0]+t*(yr[1]-yr[0]),ym.kind)+'</text>';});
  host.innerHTML='<svg width="'+S+'" height="'+H+'" id="scsvg">'+
    '<line x1="'+pad+'" y1="'+(H-pad)+'" x2="'+(pad+plot)+'" y2="'+(H-pad)+'" stroke="#cbd5e1"/><line x1="'+pad+'" y1="'+pad+'" x2="'+pad+'" y2="'+(H-pad)+'" stroke="#cbd5e1"/>'+
    ticks+dots+
    '<text x="'+(pad+plot/2)+'" y="'+(H-6)+'" font-size="11" fill="#475569" text-anchor="middle">'+esc(xm.label)+'</text>'+
    '<text x="12" y="'+(pad+ph/2)+'" font-size="11" fill="#475569" text-anchor="middle" transform="rotate(-90 12 '+(pad+ph/2)+')">'+esc(ym.label)+'</text></svg>';
  SGEO={pad:pad,plot:plot,ph:ph,H:H,xr:xr,yr:yr,xk:xk,yk:yk};   // for the rubber-band select inverse-mapping
  var scsvg=el('scsvg'); if(scsvg){var ov=document.createElementNS('http://www.w3.org/2000/svg','rect');
    ov.setAttribute('id','scbrush');ov.setAttribute('fill','rgba(14,139,168,.12)');ov.setAttribute('stroke','#0e8ba8');
    ov.setAttribute('stroke-dasharray','3 2');ov.setAttribute('pointer-events','none');ov.style.display='none';scsvg.appendChild(ov);}
}

// ---- Spearman correlation matrix across the core metrics (click a cell -> loads that pair into the scatter) ----
function rankvec(a){var idx=a.map(function(v,i){return [v,i];}).sort(function(x,y){return x[0]-y[0];});
  var r=new Array(a.length),k=0;
  while(k<idx.length){var j=k;while(j+1<idx.length&&idx[j+1][0]===idx[k][0])j++;
    var avg=(k+j)/2+1;for(var t=k;t<=j;t++)r[idx[t][1]]=avg;k=j+1;}
  return r;}
function spearman(x,y){var n=x.length; if(n<4)return null;
  var rx=rankvec(x),ry=rankvec(y),mx=0,my=0,i;
  for(i=0;i<n;i++){mx+=rx[i];my+=ry[i];} mx/=n;my/=n;
  var sxy=0,sxx=0,syy=0;
  for(i=0;i<n;i++){var dx=rx[i]-mx,dy=ry[i]-my;sxy+=dx*dy;sxx+=dx*dx;syy+=dy*dy;}
  return (sxx>0&&syy>0)?sxy/Math.sqrt(sxx*syy):null;}
function corrCol(r){if(r==null)return '#f0f3f7';var a=Math.abs(r),base=r>=0?[224,84,79]:[79,131,194],w=[247,249,252];
  return 'rgb('+w.map(function(c,i){return Math.round(c+(base[i]-c)*a);}).join(',')+')';}
function renderCorr(){
  var host=el('corr_body'); if(!host)return;
  var keys=DIST.filter(function(k){return MET[k];}).slice(0,16);
  var pool=visible();
  if(pool.length<4||keys.length<2){host.innerHTML='<div class="pad nd">not enough data for a correlation matrix (need >= 4 samples in view).</div>';return;}
  var vals={}; keys.forEach(function(k){vals[k]=pool.map(function(s){return s.m[k];});});
  var n=keys.length, cell=Math.max(15,Math.min(30,Math.floor((Math.min(host.clientWidth||520,host.classList.contains('expanded')?900:560)-100)/n)));
  var padL=94,padT=8, W=padL+n*cell+8, H=padT+n*cell+96;
  var svg='<svg width="'+W+'" height="'+H+'" id="corrsvg" style="max-width:100%;display:block">';
  keys.forEach(function(k,j){var cx=padL+j*cell+cell/2;
    svg+='<text x="'+cx+'" y="'+(padT+n*cell+12)+'" font-size="8.5" fill="#67788b" text-anchor="end" transform="rotate(-55 '+cx+' '+(padT+n*cell+12)+')">'+esc(MET[k].label)+'</text>';});
  keys.forEach(function(k,i){var cy=padT+i*cell+cell/2;
    svg+='<text x="'+(padL-6)+'" y="'+(cy+3)+'" font-size="8.5" fill="#67788b" text-anchor="end">'+esc(MET[k].label)+'</text>';});
  for(var i=0;i<n;i++)for(var j=0;j<n;j++){
    var cx=padL+j*cell,cy=padT+i*cell,r;
    if(i===j){r=vals[keys[i]].some(function(v){return v!=null;})?1:null;}   // grey out an all-NA metric
    else{var xs=[],ys=[]; for(var t=0;t<pool.length;t++){var xv=vals[keys[j]][t],yv=vals[keys[i]][t];
        if(xv!=null&&yv!=null){xs.push(xv);ys.push(yv);}} r=spearman(xs,ys);}
    svg+='<rect x="'+cx+'" y="'+cy+'" width="'+(cell-1)+'" height="'+(cell-1)+'" rx="2" fill="'+corrCol(r)+'"'+
      ' data-xk="'+keys[j]+'" data-yk="'+keys[i]+'" data-r="'+(r==null?'':r.toFixed(2))+'"'+(i!==j?' style="cursor:pointer"':'')+'/>';
    if(cell>=22&&r!=null)svg+='<text x="'+(cx+(cell-1)/2)+'" y="'+(cy+(cell-1)/2+3)+'" font-size="7.5" fill="'+(Math.abs(r)>0.55?'#fff':'#5b6b7e')+'" text-anchor="middle" pointer-events="none">'+(r>0?'':'-')+Math.abs(r).toFixed(1).replace('0.','.')+'</text>';
  }
  var ly=padT+n*cell+64;
  svg+='<text x="'+padL+'" y="'+(ly-4)+'" font-size="8.5" fill="#94a3b8">Spearman rho</text>';
  for(var g=0;g<=20;g++){var rr=-1+g/10; svg+='<rect x="'+(padL+g*7)+'" y="'+ly+'" width="7" height="9" fill="'+corrCol(rr)+'"/>';}
  svg+='<text x="'+padL+'" y="'+(ly+20)+'" font-size="8" fill="#94a3b8">-1</text>'+
       '<text x="'+(padL+70)+'" y="'+(ly+20)+'" font-size="8" fill="#94a3b8" text-anchor="middle">0</text>'+
       '<text x="'+(padL+140)+'" y="'+(ly+20)+'" font-size="8" fill="#94a3b8" text-anchor="end">+1</text>';
  host.innerHTML=svg+'</svg>';
}

function renderStacks(){
  var host=el('stacks');
  var rows=R.samples.filter(function(s){return s.m.callable_pct!=null||s.m.missing_pct!=null;})
    .slice().sort(function(a,b){return (b.m.missing_pct||0)-(a.m.missing_pct||0);});
  if(!rows.length){host.innerHTML='<div class="pad nd">no consensus data</div>';return;}
  var vis={}; visible().forEach(function(s){vis[s.s]=1;});
  host.innerHTML=rows.map(function(s){var cal=s.m.callable_pct||0,iup=s.m.iupac_pct||0,mis=Math.max(0,100-cal-iup);
    var dim=(st.q||st.onlyFlagged||st.flagFilter||st.ancOnly||st.linFilter)&&!vis[s.s];
    return '<div class="stack" data-s="'+esc(s.s)+'" style="opacity:'+(dim?0.25:1)+(st.hi==s.s?';background:#fff6d6':'')+'">'+
     '<span class="sl">'+esc(s.s)+'</span><div class="sb">'+
     '<div style="width:'+cal.toFixed(2)+'%;background:#22a06b" title="callable"></div>'+
     '<div style="width:'+iup.toFixed(2)+'%;background:#e6b25a" title="IUPAC"></div>'+
     '<div style="width:'+mis.toFixed(2)+'%;background:#cbd5e1" title="missing"></div></div>'+
     '<span class="sv">'+mis.toFixed(1)+'%</span></div>';}).join('');
  Array.prototype.forEach.call(host.querySelectorAll('.stack'),function(d){d.onclick=function(){setHi(d.getAttribute('data-s'));};});
}

function qcScore(s){var fails=s.f.filter(function(f){return FAILF[f];}).length;return fails*100+(s.f.length-fails)*10;}
function renderFlags(){
  var fl=R.samples.filter(function(s){return s.v!='PASS';}).sort(function(a,b){return qcScore(b)-qcScore(a)||a.s.localeCompare(b.s);});
  var body=fl.length?fl.map(function(s){return '<tr data-s="'+esc(s.s)+'">'+
    '<td class="s">'+esc(s.s)+(s.anc?'<span class="abadge">aDNA</span>':'')+'</td><td><span class="v '+s.v+'">'+s.v+'</span></td>'+
    '<td class="flags">'+s.f.map(function(f){return '<span class="chip'+(FAILF[f]?' failc':'')+'" data-f="'+f+'" title="'+esc(flagWhy(s,f))+'" style="cursor:pointer">'+f+'</span>';}).join(' ')+'</td></tr>';}).join('')
    :'<tr><td colspan="3" style="text-align:left;color:#16a34a;padding:10px">All samples pass at the current thresholds.</td></tr>';
  var t=el('flagtable');
  t.innerHTML='<thead><tr><th class="s">Sample (worst first)</th><th>QC</th><th style="text-align:left">Flags - hover for the margin</th></tr></thead><tbody>'+body+'</tbody>';
  Array.prototype.forEach.call(t.querySelectorAll('tbody tr[data-s]'),function(tr){tr.onclick=function(e){
    if(e.target.classList.contains('chip')){var f=e.target.getAttribute('data-f');st.flagFilter=(st.flagFilter==f?null:f);renderAll();el('gstats').scrollIntoView();return;}
    setHi(tr.getAttribute('data-s'));el('gstats').scrollIntoView();};});
  el('nflag').textContent=fl.length;
  var bw=el('basketFlagged'); if(bw){bw.onclick=function(){fl.forEach(function(s){st.excl[s.s]=1;});renderTable();renderCuration();};}
}

// ---- curation basket + exclusion exports ----
function nExcl(){return R.samples.filter(function(s){return st.excl[s.s];}).length;}
function renderCuration(){var ex=nExcl(),keep=R.samples.length-ex;
  el('curation').innerHTML='<div class="cur-intro"><b>Exclusion basket</b> - the set of samples you are dropping from the downstream analysis (phylogeny / clock). Tick a sample&#39;s box in the table (or use the buttons), then export the drop list (<b>exclusion.tsv</b>) or the survivors (<b>keep_list.txt</b>).</div>'+
    '<div class="cur-read"><b>'+ex+'</b> to exclude <span class="arw">→</span> <b>'+keep+'</b> kept for downstream</div>'+
   '<div class="cur-btns"><button class="btn" data-cur="fail">exclude FAILs</button><button class="btn" data-cur="flagged">exclude all flagged</button>'+
   '<button class="btn" data-cur="clear">clear</button><button class="btn" data-cur="invert">invert (shown)</button>'+
   '<button class="btn prim" data-cur="excl">↓ exclusion.tsv</button><button class="btn prim" data-cur="keep">↓ keep_list.txt</button></div>';
  Array.prototype.forEach.call(el('curation').querySelectorAll('[data-cur]'),function(b){b.onclick=function(){curAction(b.getAttribute('data-cur'));};});
  if(typeof saveState=='function')saveState();}
function curAction(a){
  if(a=='fail')R.samples.forEach(function(s){if(s.v=='FAIL')st.excl[s.s]=1;});
  else if(a=='flagged')R.samples.forEach(function(s){if(s.v!='PASS')st.excl[s.s]=1;});
  else if(a=='clear')st.excl={};
  else if(a=='invert')visible().forEach(function(s){if(st.excl[s.s])delete st.excl[s.s];else st.excl[s.s]=1;});
  else if(a=='excl'){exportExcl();return;} else if(a=='keep'){exportKeep();return;}
  renderTable();renderCuration();}
function exportExcl(){var lines=['sample\tverdict\tancient\tflags\treason'];
  R.samples.filter(function(s){return st.excl[s.s];}).sort(function(a,b){return a.s.localeCompare(b.s);}).forEach(function(s){
    var reason=s.f.length?s.f.map(function(f){return flagWhy(s,f);}).join('; '):'manual_qc_exclusion';
    lines.push([s.s,s.v,(s.anc?'yes':'no'),(s.f.join(';')||'.'),reason].join('\t'));});
  dl(lines.join('\n')+'\n','exclusion.tsv','text/tab-separated-values');}
function exportKeep(){var lines=R.samples.filter(function(s){return !st.excl[s.s];}).map(function(s){return s.s;}).sort();
  dl(lines.join('\n')+'\n','keep_list.txt','text/plain');}
function dl(txt,name,type){var blob=new Blob([txt],{type:type}),a=document.createElement('a');a.href=URL.createObjectURL(blob);a.download=name;a.click();URL.revokeObjectURL(a.href);}

// ---- per-sample detail modal ----
function pctRank(key,val){if(val==null)return null;var vs=R.samples.map(function(s){return s.m[key];}).filter(function(v){return v!=null;});if(!vs.length)return null;var b=0;vs.forEach(function(v){if(v<val)b++;});return Math.round(100*b/vs.length);}
function openDetail(sid){var s=null;R.samples.forEach(function(x){if(x.s==sid)s=x;});if(!s)return;st.detail=sid;
  var rows=R.metrics.map(function(m){var v=s.m[m.key],r=RANGES[m.key],nn=(v==null||r[1]<=r[0])?0:Math.max(0,Math.min(1,(v-r[0])/(r[1]-r[0]))),pr=pctRank(m.key,v);
    return '<div class="drow"><div class="dk" title="'+esc((R.defs[m.key]||[''])[0])+'">'+esc(m.label)+'</div>'+
      '<div class="dbarwrap"><div class="dbar" style="width:'+(nn*100).toFixed(1)+'%;background:'+BAR[m.dir]+'"></div></div>'+
      '<div class="dv">'+(v==null?'<span class="na">NA</span>':fmt(v,m.kind))+'</div><div class="dp">'+(pr==null?'':('p'+pr))+'</div></div>';}).join('');
  var lin='';if(s.linf){var ks=Object.keys(s.linf).sort(function(a,b){return s.linf[b]-s.linf[a];});
    lin='<div class="dsub">Lineage composition</div><div class="lcomp">'+ks.map(function(k){return '<div class="lrow"><span class="lk">'+esc(k)+'</span><div class="lbarw"><div class="lbar" style="width:'+(s.linf[k]*100).toFixed(1)+'%"></div></div><span class="lv">'+(s.linf[k]*100).toFixed(1)+'%</span></div>';}).join('')+'</div>';}
  var dmg='';if(s.anc){var d=s.dmg;
    if(d){dmg='<div class="dsub">aDNA damage (mapDamage2) - consistent with ancient DNA, not proof of authenticity</div>'+
      '<div class="drow"><div class="dk">5&#39; C&gt;T (pos 1)</div><div class="dbarwrap"><div class="dbar" style="width:'+Math.min(100,(d.ct1||0)*100/0.3).toFixed(1)+'%;background:'+(s.f.indexOf('DAMAGE_LOW')>=0?'#e0544f':'#22a06b')+'"></div></div><div class="dv">'+(d.ct1==null?'NA':(d.ct1*100).toFixed(1)+'%')+'</div><div class="dp"></div></div>'+
      '<div class="drow"><div class="dk">3&#39; G&gt;A (pos 1)</div><div class="dbarwrap"><div class="dbar" style="width:'+Math.min(100,(d.ga1||0)*100/0.3).toFixed(1)+'%;background:#4f83c2"></div></div><div class="dv">'+(d.ga1==null?'NA':(d.ga1*100).toFixed(1)+'%')+'</div><div class="dp"></div></div>'+
      '<div class="drow"><div class="dk">mean frag len</div><div class="dv" style="flex:1;text-align:left;color:#516074">'+(d.fraglen==null?'NA':d.fraglen.toFixed(0)+' bp')+'</div></div>';
    }else{dmg='<div class="dsub">aDNA damage</div><div class="nd" style="padding:4px 0">no mapDamage2 output found.</div>';}}
  var fl=s.f.length?s.f.map(function(f){return '<span class="chip'+(FAILF[f]?' failc':'')+'" title="'+esc(flagWhy(s,f))+'">'+f+'</span>';}).join(' '):'<span style="color:#16a34a">no flags ✓</span>';
  var why=s.f.length?'<div class="dwhy">'+s.f.map(function(f){return '<div><b>'+esc(f)+'</b> &middot; '+esc(flagWhy(s,f))+'</div>';}).join('')+'</div>':'';
  var annb='';
  (function(){if(s.m.ann_high==null&&s.m.annotated_pct==null)return;
    var warn=s.m.snpeff_warn, tot=s.m.total_variants, badFrac=(warn!=null&&tot)?warn/tot:null;
    var annBad=(s.m.annotated_pct!=null&&s.m.annotated_pct<50)||s.ann_db_error;
    var cls=(annBad||(badFrac!=null&&badFrac>0.5))?' low':'';
    annb='<div class="dsub">Functional annotation (snpEff) <span style="font-weight:400;color:#94a3b8">- impact + annotation QC</span></div>'+
      '<div class="acards">'+
      '<div class="acard'+cls+'"><div class="acard-h"><span class="sname">Impact</span></div>'+
        '<div class="astats">'+
        '<div class="astat"><span class="ak">HIGH</span><span class="av'+(s.m.ann_high>0?' bad':'')+'">'+(s.m.ann_high==null?'NA':Math.round(s.m.ann_high).toLocaleString("en-US"))+'</span></div>'+
        '<div class="astat"><span class="ak">MODERATE</span><span class="av">'+(s.m.ann_moderate==null?"NA":Math.round(s.m.ann_moderate).toLocaleString("en-US"))+'</span></div>'+
        '<div class="astat"><span class="ak">LoF %</span><span class="av">'+(s.m.lof_pct==null?"NA":s.m.lof_pct.toFixed(1)+"%")+'</span></div>'+
        '<div class="astat"><span class="ak">missense/silent</span><span class="av">'+(s.m.missense_silent==null?"NA":s.m.missense_silent.toFixed(2))+'</span></div>'+
        '</div></div>'+
      '<div class="acard'+cls+'"><div class="acard-h"><span class="sname">Annotation QC</span></div>'+
        '<div class="astats">'+
        '<div class="astat"><span class="ak">annotated</span><span class="av'+(annBad?' bad':' good')+'">'+(s.m.annotated_pct==null?"NA":s.m.annotated_pct.toFixed(1)+"%")+'</span></div>'+
        '<div class="astat"><span class="ak">coding</span><span class="av">'+(s.m.coding_pct==null?"NA":s.m.coding_pct.toFixed(1)+"%")+'</span></div>'+
        '<div class="astat"><span class="ak">snpEff warnings</span><span class="av'+(badFrac!=null&&badFrac>0.5?' bad':'')+'">'+(s.m.snpeff_warn==null?"NA":Math.round(s.m.snpeff_warn).toLocaleString("en-US"))+'</span></div>'+
        '</div>'+((annBad||(badFrac!=null&&badFrac>0.5))?'<div class="alow">Low annotation coverage, a snpEff database error, or many snpEff warnings: the reference GFF-to-database build may be wrong for this contig; treat the impact counts and the missense/silent proxy for this sample with suspicion.</div>':'')+
      '</div></div>';})();
  el('modalbody').innerHTML='<div class="dhead"><div><div class="dtitle">'+esc(s.s)+'</div><div class="dmeta">'+esc(s.lineage||'lineage NA')+(s.anc?' &middot; <b style="color:#8a5a12">aDNA</b>':'')+(s.dr?' &middot; DR: '+esc(s.dr):'')+(s.date?' &middot; '+esc(s.date):'')+'</div></div>'+
    '<span class="v '+s.v+'" style="font-size:12px">'+s.v+'</span></div><div class="dflags">'+fl+'</div>'+why+genomeSpark(s)+lin+dmg+annb+
    '<div class="dsub">All metrics <span style="font-weight:400;color:#94a3b8">(bar = position in cohort range · p = percentile)</span></div>'+rows+
    '<div class="dbtns"><button class="btn" id="dexcl"></button></div>';
  var dx=el('dexcl');function setlbl(){dx.textContent=st.excl[s.s]?'✓ excluded - click to keep':'exclude this sample';dx.classList.toggle('prim',!st.excl[s.s]);}
  setlbl();dx.onclick=function(){if(st.excl[s.s])delete st.excl[s.s];else st.excl[s.s]=1;setlbl();renderTable();renderCuration();};
  el('modal').classList.add('open');}
function closeDetail(){st.detail=null;el('modal').classList.remove('open');}

// ---- genome landscape: samples x reference-position heatmap, multi-track (missing / SNPs / het / indels) ----
var GTRACKS=[{k:'missing',lab:'Missing',base:[214,64,58]},{k:'snp',lab:'SNPs',base:[14,139,168]},
             {k:'het',lab:'Het',base:[221,138,26]},{k:'indel',lab:'Indels',base:[124,92,191]}];
function gtrackHas(k){return R.samples.some(function(s){return k=='missing'?!!s.miss:(s.trk&&s.trk[k]);});}
function gtrackMax(k){var mx=0;R.samples.forEach(function(s){var p=k=='missing'?s.miss:(s.trk&&s.trk[k]);if(p)p.forEach(function(v){if(v!=null&&v>mx)mx=v;});});return mx;}
function gcol(base,mv){if(mv==null)return '#e9edf2';if(mv>1)mv=1;var a=[238,244,240];
  return 'rgb('+Math.round(a[0]+(base[0]-a[0])*mv)+','+Math.round(a[1]+(base[1]-a[1])*mv)+','+Math.round(a[2]+(base[2]-a[2])*mv)+')';}
// ---- Functional annotation (snpEff impact + effect classes per sample) ----
var IMP=[['ann_high','HIGH','#e0544f'],['ann_moderate','MODERATE','#e6b25a'],
         ['ann_low','LOW','#4f83c2'],['ann_modifier','MODIFIER','#cbd5e1']];
var EFFCLS=[['eff_missense','missense','#4f83c2'],['eff_synonymous','synonymous','#22a06b'],
            ['eff_stop_gained','stop gained','#e0544f'],['eff_frameshift','frameshift','#c0392b'],
            ['eff_start_lost','start lost','#b5651d'],['eff_stop_lost','stop lost','#d98c2b'],
            ['eff_inframe_indel','inframe indel','#7c5cbf'],['eff_splice','splice','#a0508a'],
            ['eff_intergenic','intergenic','#9aa7b6'],['eff_regulatory','regulatory','#c4ccd7']];
var FUNCTION_CAPTION='snpEff functional annotation per sample. The bar splits each sample&#39;s annotated variants into snpEff impact classes (HIGH, MODERATE, LOW, MODIFIER; most severe first). Impact is snpEff&#39;s own severity call from the reference gene model, so these counts inherit any error in the reference annotation, and a HIGH-impact excess can be genuine loss-of-function, a wrong reference database, contamination, or an indel-calling artifact rather than biology. The missense/silent ratio shown is a pN/pS PROXY (a spectrum sanity check within a lineage), NOT dN/dS and NOT a test of selection; it is entangled with divergence and with how snpEff bins effects, so read it alongside the Ti/Tv ratio and use the eskaks pN/pS panel for a real estimate. Counts are absolute (a divergent lineage has more of every class); the excess-HIGH and excess-LoF flags are cohort-relative robust-z, not absolute cut-offs.';
function renderFunction(){
  var host=el('fn_stacks'),cap=el('fn_caption'),leg=el('fn_legend'),coh=el('fn_cohort');
  if(!host)return;
  var rows=R.samples.filter(function(s){return IMP.some(function(p){return s.m[p[0]]!=null;});})
    .slice().sort(function(a,b){return (b.m.ann_high||0)-(a.m.ann_high||0);});
  if(!rows.length){host.innerHTML='<div class="pad nd">no snpEff annotation in the summary (pass annotated VCFs upstream).</div>';if(cap)cap.innerHTML='';if(leg)leg.innerHTML='';if(coh)coh.innerHTML='';return;}
  var vis={}; visible().forEach(function(s){vis[s.s]=1;});
  host.innerHTML=rows.map(function(s){
    var tot=IMP.reduce(function(a,p){return a+(s.m[p[0]]||0);},0)||1;
    var dim=(st.q||st.onlyFlagged||st.flagFilter||st.ancOnly||st.linFilter)&&!vis[s.s];
    var segs=IMP.map(function(p){var v=s.m[p[0]]||0,w=100*v/tot;if(w<=0)return '';
      return '<div style="width:'+w.toFixed(2)+'%;background:'+p[2]+'" title="'+p[1]+': '+Math.round(v).toLocaleString('en-US')+'"></div>';}).join('');
    var hi=s.m.ann_high!=null?Math.round(s.m.ann_high).toLocaleString('en-US'):'-';
    return '<div class="stack" data-s="'+esc(s.s)+'" style="opacity:'+(dim?0.25:1)+(st.hi==s.s?';background:#fff6d6':'')+'">'+
      '<span class="sl">'+esc(s.s)+'</span><div class="sb">'+segs+'</div>'+
      '<span class="sv" title="HIGH-impact count">'+hi+'</span></div>';}).join('');
  Array.prototype.forEach.call(host.querySelectorAll('.stack'),function(d){d.onclick=function(){setHi(d.getAttribute('data-s'));};});
  if(leg)leg.innerHTML=IMP.map(function(p){return '<span><i style="background:'+p[2]+'"></i>'+p[1]+'</span>';}).join('')+'<span style="margin-left:auto">right number = HIGH-impact count; sorted worst first</span>';
  var pool=visible().filter(function(s){return EFFCLS.some(function(p){return s.m[p[0]]!=null;});});
  var sums={},grand=0;EFFCLS.forEach(function(p){var t=0;pool.forEach(function(s){if(s.m[p[0]]!=null)t+=s.m[p[0]];});sums[p[0]]=t;grand+=t;});
  var mx=1;EFFCLS.forEach(function(p){if(sums[p[0]]>mx)mx=sums[p[0]];});
  var effRows=EFFCLS.filter(function(p){return sums[p[0]]>0;}).map(function(p){var v=sums[p[0]];
    return '<div class="drow"><div class="dk">'+esc(p[1])+'</div><div class="dbarwrap"><div class="dbar" style="width:'+(100*v/mx).toFixed(1)+'%;background:'+p[2]+'"></div></div>'+
      '<div class="dv">'+Math.round(v).toLocaleString('en-US')+'</div><div class="dp">'+(grand?Math.round(100*v/grand):0)+'%</div></div>';}).join('');
  var mis=sums.eff_missense||0,syn=sums.eff_synonymous||0,ratio=(syn>0)?(mis/syn):null;
  var strip='<div class="cur-read" style="margin:6px 0 10px"><b>'+(ratio==null?'NA':ratio.toFixed(2))+'</b> cohort missense/silent ratio '+
    '<span style="color:#94a3b8;font-weight:400">(pN/pS proxy; '+Math.round(mis).toLocaleString('en-US')+' missense / '+Math.round(syn).toLocaleString('en-US')+' synonymous; not a selection test)</span></div>';
  if(coh)coh.innerHTML='<div style="flex-basis:100%"><div class="dsub" style="margin:2px 0 4px">Cohort effect classes <span style="font-weight:400;color:#94a3b8">(summed over samples in view)</span></div>'+strip+effRows+'</div>';
  if(cap)cap.innerHTML=FUNCTION_CAPTION;
}
function renderGenome(){
  var host=el('genome_body'); if(!host)return;
  var samp=R.samples.filter(function(s){return s.miss||s.trk;});
  if(!samp.length){host.innerHTML='<span class="nd" style="padding:0">no consensus/variant data for a genome landscape.</span>';return;}
  var tk=st.gtrack||'missing'; if(!gtrackHas(tk)){tk='missing';st.gtrack='missing';}
  var meta=GTRACKS[0]; GTRACKS.forEach(function(g){if(g.k==tk)meta=g;}); var base=meta.base;
  var rows0=samp.filter(function(s){return (tk=='missing')?s.miss:(s.trk&&s.trk[tk]);});
  if(!rows0.length){host.innerHTML='<span class="nd" style="padding:0">no data for this track.</span>';return;}
  var nb=(tk=='missing')?rows0[0].miss.length:rows0[0].trk[tk].length, gl=R.genome_len||nb;
  var mx=(tk=='missing')?1:(gtrackMax(tk)||1);
  var vis={}; visible().forEach(function(s){vis[s.s]=1;});
  var rows=rows0.slice().sort(function(a,b){var la=(R.lineages||[]).indexOf(a.lineage),lb=(R.lineages||[]).indexOf(b.lineage);
    if(la<0)la=999; if(lb<0)lb=999; if(la!=lb)return la-lb; return a.s.localeCompare(b.s);});
  var W=Math.max(280,host.clientWidth||900), gut=54, profH=54, rowH=Math.max(3,Math.min(9,Math.floor(300/rows.length))), plotW=W-gut-8;
  var z0=0,z1=nb-1; if(st.gzoom){z0=Math.max(0,Math.min(nb-1,st.gzoom.b0|0));z1=Math.max(z0,Math.min(nb-1,st.gzoom.b1|0));}
  var winN=z1-z0+1;                                            // zoom window (bins); the whole reference when st.gzoom is null
  function x(i){return gut+(i-z0)/winN*plotW;} var bw=plotW/winN+0.6;
  function raw(s,i){if(tk=='missing')return s.miss?s.miss[i]:null;var p=s.trk&&s.trk[tk];return p?p[i]:null;}
  var visRows=rows.filter(function(s){return vis[s.s];}); if(!visRows.length)visRows=rows;
  var agg=[];for(var i=0;i<nb;i++){var sum=0,c=0;visRows.forEach(function(s){var vv=raw(s,i);if(vv!=null){sum+=vv;c++;}});agg.push(c?sum/c:null);}
  var mb=(st.maskOn&&R.mask_bins)?R.mask_bins:null;
  if(mb)for(var mi=0;mi<nb;mi++)if(mb[mi]>=0.5)agg[mi]=null;   // masked bins drop out of the density profile
  var y0=profH+14, totH=y0+rows.length*rowH+18;
  // cohort DENSITY PROFILE (area + line) of the selected track along the reference
  var pmax=Math.max.apply(null,agg.slice(z0,z1+1).map(function(v){return v||0;}).concat([1e-9]));
  var pts=[];for(var pi=z0;pi<=z1;pi++){var pv=agg[pi]||0;pts.push(x(pi).toFixed(1)+','+(profH-(pv/pmax)*(profH-10)).toFixed(1));}
  var bs='rgb('+base[0]+','+base[1]+','+base[2]+')',bf='rgba('+base[0]+','+base[1]+','+base[2]+',.15)';
  var TLAB={missing:'missing %',snp:'SNP density',het:'het density',indel:'indel density'};
  var svg='<svg width="'+W+'" height="'+totH+'">';
  svg+='<line x1="'+gut+'" y1="'+profH+'" x2="'+(gut+plotW).toFixed(1)+'" y2="'+profH+'" stroke="#e6ebf1"/>';
  svg+='<path d="M '+gut.toFixed(1)+' '+profH+' L '+pts.join(' L ')+' L '+(gut+plotW).toFixed(1)+' '+profH+' Z" fill="'+bf+'"/>';
  svg+='<polyline points="'+pts.join(' ')+'" fill="none" stroke="'+bs+'" stroke-width="1.3" stroke-linejoin="round"/>';
  svg+='<text x="0" y="11" font-size="9" font-weight="600" fill="#67788b">'+TLAB[tk]+'</text>';
  svg+='<text x="'+(gut+plotW).toFixed(1)+'" y="11" font-size="8.5" fill="#94a3b8" text-anchor="end">peak '+(tk=='missing'?(pmax*100).toFixed(0)+'%':(Math.round(pmax*10)/10)+'/bin')+'</text>';
  for(var pj=z0;pj<=z1;pj++)svg+='<rect x="'+x(pj).toFixed(1)+'" y="0" width="'+bw.toFixed(1)+'" height="'+profH+'" fill="transparent" data-bin="'+pj+'" data-v="'+(agg[pj]==null?'':(tk=='missing'?(agg[pj]*100).toFixed(0):(Math.round(agg[pj]*10)/10)))+'"/>';
  rows.forEach(function(s,r){var dim=(st.q||st.onlyFlagged||st.flagFilter||st.ancOnly||st.linFilter)&&!vis[s.s],yy=y0+r*rowH;
    svg+='<rect x="'+(gut-9)+'" y="'+yy+'" width="5" height="'+(rowH-0.4).toFixed(1)+'" fill="'+linColor(s.lineage)+'"'+(dim?' opacity="0.25"':'')+'/>';
    for(var i=z0;i<=z1;i++){var vv=raw(s,i),mv=vv==null?null:(tk=='missing'?vv:vv/mx);
      svg+='<rect x="'+x(i).toFixed(1)+'" y="'+yy+'" width="'+bw.toFixed(1)+'" height="'+(rowH-0.4).toFixed(1)+'" fill="'+gcol(base,mv)+'"'+(dim?' opacity="0.3"':'')+' data-s="'+esc(s.s)+'" data-bin="'+i+'" data-v="'+(vv==null?'':(tk=='missing'?(vv*100).toFixed(0):vv))+'"/>';}
  });
  var yb=y0+rows.length*rowH;
  if(mb)for(var mk=z0;mk<=z1;mk++){var mf=mb[mk]||0;if(mf<=0)continue;   // grey out the masked zones
    svg+='<rect x="'+x(mk).toFixed(1)+'" y="0" width="'+bw.toFixed(1)+'" height="'+yb+'" fill="#5b6b7e" opacity="'+(0.12+0.42*mf).toFixed(2)+'"/>';}
  [0,0.25,0.5,0.75,1].forEach(function(t){var px=gut+t*plotW,pos=Math.round((z0+t*winN)/nb*gl);
    svg+='<line x1="'+px.toFixed(1)+'" y1="'+yb+'" x2="'+px.toFixed(1)+'" y2="'+(yb+4)+'" stroke="#cbd5e1"/>'+
      '<text x="'+px.toFixed(1)+'" y="'+(yb+14)+'" font-size="8.5" fill="#94a3b8" text-anchor="'+(t==0?'start':t==1?'end':'middle')+'">'+fmtpos(pos)+'</text>';});
  if(tk=='snp'&&st.geneMark&&st.geneMark.b1>=z0&&st.geneMark.b0<=z1){var gm=st.geneMark,mx0=Math.max(gut,x(gm.b0)),mx1=Math.min(gut+plotW,x(gm.b1+1));
    svg+='<rect x="'+mx0.toFixed(1)+'" y="0" width="'+Math.max(2,mx1-mx0).toFixed(1)+'" height="'+yb+'" fill="none" stroke="#0f1c29" stroke-width="1.2" stroke-dasharray="3 2"/>'+
      '<text x="'+Math.min(W-2,(mx0+mx1)/2).toFixed(1)+'" y="'+(profH+11)+'" font-size="9" font-weight="600" fill="#0f1c29" text-anchor="middle">'+esc(gm.name||'')+'</text>';}
  svg+='<rect id="gbrush" x="0" y="0" width="0" height="'+yb+'" fill="rgba(14,139,168,.12)" stroke="#0e8ba8" stroke-width="1" stroke-dasharray="3 2" pointer-events="none" style="display:none"/>';
  GGEO={gut:gut,plotW:plotW,nb:nb,gl:gl,yb:yb,z0:z0,winN:winN};
  host.innerHTML=svg+'</svg>';
  // selection overlay only when it is a SUB-region of the current view (when zoomed to exactly the selection, the whole plot is it)
  if(st.gsel&&!(st.gsel.b0<=z0&&st.gsel.b1>=z1)){var gbx=el('gbrush');if(gbx){var bx0=Math.max(gut,x(st.gsel.b0)),bx1=Math.min(gut+plotW,x(st.gsel.b1+1));
    if(bx1>bx0){gbx.setAttribute('x',bx0.toFixed(1));gbx.setAttribute('width',(bx1-bx0).toFixed(1));gbx.style.display='';}else gbx.style.display='none';}}
}

// ---- SNP-dense gene / region detection (cohort SNP density along the reference) ----
function cohortSnp(){var nb=R.nbins||200,agg=new Array(nb).fill(0),has=false,mb=(st.maskOn&&R.mask_bins)?R.mask_bins:null;
  R.samples.forEach(function(s){var p=s.trk&&s.trk.snp;if(p){has=true;for(var i=0;i<nb&&i<p.length;i++){if(mb&&mb[i]>=0.5)continue;agg[i]+=p[i]||0;}}});
  return has?agg:null;}
function robustMS(a){var v=a.slice().sort(function(x,y){return x-y;}),n=v.length;if(!n)return[0,1];
  var med=n%2?v[(n-1)/2]:(v[n/2-1]+v[n/2])/2,d=a.map(function(x){return Math.abs(x-med);}).sort(function(x,y){return x-y;});
  var mad=n%2?d[(n-1)/2]:(d[n/2-1]+d[n/2])/2,sig=1.4826*mad;
  if(sig<=0){var m=a.reduce(function(p,q){return p+q;},0)/n;sig=Math.sqrt(a.reduce(function(p,q){return p+(q-m)*(q-m);},0)/n)||1;}
  return[med,sig];}
function renderHotspots(){
  var host=el('hot_body'),tbl=el('hottable'),note=el('hot_note'); if(!host||!tbl)return;
  var agg=cohortSnp();
  if(!agg){host.style.display='none';return;} host.style.display='';
  var nb=R.nbins||agg.length, gl=R.genome_len||nb, binbp=gl/nb, genes=R.genes||[], head, body;
  function b0of(p){return Math.max(0,Math.min(nb-1,Math.floor(p/gl*nb)));}
  function geneMaskedFrac(s,e){var iv=R.mask_iv;if(!iv)return 0;var len=(e-s+1)||1,ov=0;
    for(var i=0;i<iv.length;i++){var a=iv[i][0],b=iv[i][1];if(b<=s||a>=e)continue;ov+=Math.max(0,Math.min(e,b)-Math.max(s,a));}
    return ov/len;}
  if(genes.length){
    var items=genes.map(function(g){var s=Math.max(1,g.start),e=Math.max(s,g.end),len=(e-s+1)||1,b0=b0of(s),b1=b0of(e),tot=0;
      for(var b=b0;b<=b1;b++){var binS=b*binbp,binE=(b+1)*binbp,ov=Math.max(0,Math.min(e,binE)-Math.max(s,binS));tot+=agg[b]*(binbp>0?ov/binbp:1);}
      return {name:g.name,start:s,end:e,len:len,snp:tot,dens:tot/(len/1000),b0:b0,b1:b1};});
    if(st.maskOn&&R.mask_iv)items=items.filter(function(x){return geneMaskedFrac(x.start,x.end)<0.5;});   // drop masked genes
    if(st.gsel)items=items.filter(function(x){return x.b1>=st.gsel.b0&&x.b0<=st.gsel.b1;});   // brushed span only
    if(st.hotq){var _hq=st.hotq.toLowerCase();items=items.filter(function(x){return (x.name||'').toLowerCase().indexOf(_hq)>=0;});}   // gene search
    var ms=robustMS(items.map(function(x){return x.dens;}));
    items.forEach(function(x){x.z=ms[1]?(x.dens-ms[0])/ms[1]:0;});
    items.sort(function(a,b){return b.dens-a.dens;});
    var _hcap=st.hotq?300:((el('hotspotsPanel')&&el('hotspotsPanel').classList.contains('expanded'))?150:15);
    head='<tr><th class="s" style="text-align:left">Gene</th><th style="text-align:left">Position</th><th>Length</th><th>Cohort SNPs</th><th>SNPs/kb</th><th>robust z</th></tr>';
    body=items.slice(0,_hcap).map(function(x){return '<tr data-b0="'+x.b0+'" data-b1="'+x.b1+'" data-name="'+esc(x.name)+'">'+
      '<td class="s" style="text-align:left">'+esc(x.name)+'</td><td style="text-align:left">'+fmtpos(x.start)+' - '+fmtpos(x.end)+'</td>'+
      '<td>'+Math.round(x.len).toLocaleString('en-US')+'</td><td>'+Math.round(x.snp).toLocaleString('en-US')+'</td><td>'+x.dens.toFixed(1)+'</td>'+
      '<td'+(x.z>=3?' style="color:var(--fail);font-weight:600"':'')+'>'+(x.z>0?'+':'')+x.z.toFixed(1)+'</td></tr>';}).join('');
    note.textContent=(st.hotq?('Showing '+Math.min(items.length,_hcap)+' of '+items.length+' matching genes. '):'')+'Genes ranked by cohort SNP density (SNPs per kb, summed over all samples); z = robust outlier vs the gene distribution. Resolution is limited to the '+(binbp/1000).toFixed(0)+' kb profile bins.';
    if(!items.length)body='<tr><td class="na" colspan="6" style="text-align:left;padding:8px">no gene matches "'+esc(st.hotq||'')+'".</td></tr>';
  }else{
    var ms2=robustMS(agg),thr=ms2[0]+3*ms2[1],hot=[];
    for(var i=0;i<nb;i++)if(agg[i]>thr){if(hot.length&&hot[hot.length-1].b1==i-1)hot[hot.length-1].b1=i;else hot.push({b0:i,b1:i});}
    hot.forEach(function(h){h.snp=0;for(var b=h.b0;b<=h.b1;b++)h.snp+=agg[b];});
    if(st.gsel)hot=hot.filter(function(h){return h.b1>=st.gsel.b0&&h.b0<=st.gsel.b1;});
    hot.sort(function(a,b){return b.snp-a.snp;});
    head='<tr><th class="s" style="text-align:left">Region</th><th>Bins</th><th>Cohort SNPs</th></tr>';
    body=hot.slice(0,15).map(function(h){return '<tr data-b0="'+h.b0+'" data-b1="'+h.b1+'"><td class="s" style="text-align:left">'+fmtpos(Math.round(h.b0*binbp))+' - '+fmtpos(Math.round((h.b1+1)*binbp))+'</td><td>'+(h.b1-h.b0+1)+'</td><td>'+Math.round(h.snp).toLocaleString('en-US')+'</td></tr>';}).join('');
    note.textContent='SNP-dense regions (robust-z > 3 over the '+(binbp/1000).toFixed(0)+' kb bins). Pass a gene GFF (--gff) to name the genes.';
  }
  if(st.maskOn&&R.mask_bins)note.textContent+=' Masked regions excluded.';
  if(st.gsel)note.textContent+=' Restricted to the brushed span.';
  tbl.innerHTML='<thead>'+head+'</thead><tbody>'+(body||'<tr><td colspan="'+(genes.length?6:3)+'" style="text-align:left;color:#16a34a;padding:8px">no SNP-dense outliers.</td></tr>')+'</tbody>';
  Array.prototype.forEach.call(tbl.querySelectorAll('tbody tr[data-b0]'),function(tr){tr.onclick=function(){
    var _gb0=+tr.getAttribute('data-b0'),_gb1=+tr.getAttribute('data-b1'),_nb=R.nbins||200,_gp=Math.max(3,Math.round((_gb1-_gb0)*0.6)+2);
    st.gtrack='snp'; st.geneMark={b0:_gb0,b1:_gb1,name:tr.getAttribute('data-name')||''}; st.gzoom={b0:Math.max(0,_gb0-_gp),b1:Math.min(_nb-1,_gb1+_gp)};   // zoom to the clicked gene
    Array.prototype.forEach.call(document.querySelectorAll('#gtrack button'),function(x){x.classList.toggle('on',x.getAttribute('data-gt')=='snp');});
    renderGenome(); el('genome').scrollIntoView();};});
}

// ---- per-lineage summary + cohort composition bar ----
function med(a){if(!a.length)return null;var s=a.slice().sort(function(x,y){return x-y;});return s[Math.floor(s.length/2)];}
function renderLineages(){
  if(!R.lin_present)return;
  var groups={}; R.samples.forEach(function(s){var k=s.lineage||'NA';(groups[k]=groups[k]||[]).push(s);});
  var order=R.lineages.slice(); if(groups['NA'])order.push('NA');
  var tot=R.samples.length||1;
  var comp='<div class="lincomp-bar">'+order.map(function(k){var n=(groups[k]||[]).length;if(!n)return '';var w=100*n/tot;
    return '<div class="lseg" style="width:'+w.toFixed(3)+'%;background:'+linColor(k=='NA'?null:k)+'" title="'+esc(k)+': '+n+' ('+w.toFixed(1)+'%)"></div>';}).join('')+'</div>'+
    '<div class="lincomp-lab">'+order.map(function(k){var n=(groups[k]||[]).length;if(!n)return '';
      return '<span class="lchip"><i style="background:'+linColor(k=='NA'?null:k)+'"></i>'+esc(k)+' <b>'+n+'</b></span>';}).join('')+'</div>';
  el('lincomp').innerHTML=comp;
  var head='<tr><th class="s" style="text-align:left">Lineage</th><th>n</th><th>%PASS</th><th>med depth</th><th>med breadth</th><th>med SNPs</th><th># MIXED</th></tr>';
  var body=order.map(function(k){var g=groups[k]||[];if(!g.length)return '';
    var np=g.filter(function(s){return s.v=='PASS';}).length,pct=g.length?100*np/g.length:0;
    var md=med(g.map(function(s){return s.m.mean_depth;}).filter(function(v){return v!=null;}));
    var mb=med(g.map(function(s){return s.m.breadth_pct;}).filter(function(v){return v!=null;}));
    var ms=med(g.map(function(s){return s.m.snps;}).filter(function(v){return v!=null;}));
    var nmix=g.filter(function(s){return s.f.indexOf('MIXED')>=0;}).length;
    return '<tr data-lin="'+esc(k)+'"><td class="s" style="text-align:left"><span class="ldot" style="background:'+linColor(k=='NA'?null:k)+'"></span>'+esc(k)+'</td>'+
      '<td>'+g.length+'</td><td>'+pct.toFixed(0)+'</td>'+
      '<td>'+(md==null?'<span class="na">NA</span>':md.toFixed(1))+'</td><td>'+(mb==null?'<span class="na">NA</span>':mb.toFixed(1))+'</td>'+
      '<td>'+(ms==null?'<span class="na">NA</span>':Math.round(ms).toLocaleString('en-US'))+'</td>'+
      '<td'+(nmix>0?' style="color:var(--fail);font-weight:600"':'')+'>'+nmix+'</td></tr>';}).join('');
  el('linsumtable').innerHTML='<thead>'+head+'</thead><tbody>'+body+'</tbody>';
  Array.prototype.forEach.call(el('linsumtable').querySelectorAll('tbody tr'),function(tr){tr.onclick=function(){
    var k=tr.getAttribute('data-lin'); st.linFilter=(st.linFilter==(k=='NA'?null:k))?null:(k=='NA'?null:k);
    renderAll();el('gstats').scrollIntoView();};});
}

// per-sample genome sparkline for the detail modal: missing fraction (area) + SNP density (line) along the reference
function genomeSpark(s){
  var miss=s.miss, snp=s.trk&&s.trk.snp;
  if(!miss&&!snp)return '';
  var nb=(miss?miss.length:snp.length), gl=R.genome_len||nb;
  var W=512,H=46,pl=4,pr=4,pt=6,pb=12,iw=W-pl-pr,ih=H-pt-pb;
  function x(i){return pl+i/nb*iw;}
  var svg='<svg width="'+W+'" height="'+H+'" viewBox="0 0 '+W+' '+H+'" preserveAspectRatio="none" style="width:100%;height:'+H+'px">';
  svg+='<line x1="'+pl+'" y1="'+(pt+ih)+'" x2="'+(W-pr)+'" y2="'+(pt+ih)+'" stroke="#e6ebf1"/>';
  if(R.mask_bins)for(var mk=0;mk<nb&&mk<R.mask_bins.length;mk++){var mf=R.mask_bins[mk]||0;if(mf<0.5)continue;
    svg+='<rect x="'+x(mk).toFixed(1)+'" y="'+pt+'" width="'+(iw/nb+0.6).toFixed(1)+'" height="'+ih+'" fill="#5b6b7e" opacity="0.10"/>';}
  if(miss){var pts=[];for(var i=0;i<nb;i++){var mv=miss[i]==null?0:miss[i];pts.push(x(i).toFixed(1)+','+(pt+ih-mv*ih).toFixed(1));}
    svg+='<path d="M '+pl+' '+(pt+ih)+' L '+pts.join(' L ')+' L '+(W-pr)+' '+(pt+ih)+' Z" fill="rgba(214,64,58,.16)"/>'+
         '<polyline points="'+pts.join(' ')+'" fill="none" stroke="#d6403a" stroke-width="1"/>';}
  if(snp){var smax=Math.max.apply(null,snp.concat([1]));var sp=[];for(var j=0;j<nb;j++){sp.push(x(j).toFixed(1)+','+(pt+ih-(snp[j]||0)/smax*ih).toFixed(1));}
    svg+='<polyline points="'+sp.join(' ')+'" fill="none" stroke="#0e8ba8" stroke-width="1" opacity="0.75"/>';}
  [0,0.5,1].forEach(function(t){svg+='<text x="'+x(nb*t).toFixed(1)+'" y="'+(H-2)+'" font-size="7.5" fill="#94a3b8" text-anchor="'+(t==0?'start':t==1?'end':'middle')+'">'+fmtpos(Math.round(t*gl))+'</text>';});
  svg+='<text x="'+(W-pr)+'" y="'+(pt+6)+'" font-size="7.5" fill="#94a3b8" text-anchor="end">'+(snp?'missing (red) / SNP density (blue)':'missing along reference')+'</text>';
  return '<div class="dsub">Genome profile <span style="font-weight:400;color:#94a3b8">(is the loss one block or scattered?)</span></div><div style="margin-bottom:6px">'+svg+'</svg></div>';
}
// ---- aDNA damage authentication panel ----
function ctSpark(prof,lim){
  var W=132,H=34,pl=3,pr=3,pt=4,pb=8,iw=W-pl-pr,ih=H-pt-pb;
  if(!prof||!prof.length)return '<span class="nd" style="padding:0">no C&gt;T profile</span>';
  var n=Math.min(prof.length,25),xs=function(i){return pl+(n<=1?0:i/(n-1)*iw);};
  var ymax=Math.max(0.30,Math.max.apply(null,prof.slice(0,n).map(function(p){return p[1];})));
  var ys=function(v){return pt+ih-(v/ymax)*ih;};
  var pts=prof.slice(0,n).map(function(p,i){return xs(i).toFixed(1)+','+ys(p[1]).toFixed(1);}).join(' ');
  var tl=(lim!=null)?'<line x1="'+pl+'" y1="'+ys(lim).toFixed(1)+'" x2="'+(W-pr)+'" y2="'+ys(lim).toFixed(1)+'" stroke="#e0544f" stroke-width="1" stroke-dasharray="2 2"/>':'';
  return '<svg width="'+W+'" height="'+H+'" style="overflow:visible">'+
    '<line x1="'+pl+'" y1="'+(pt+ih)+'" x2="'+(W-pr)+'" y2="'+(pt+ih)+'" stroke="#e6ebf1"/>'+tl+
    '<polyline points="'+pts+'" fill="none" stroke="var(--accent)" stroke-width="1.5" stroke-linejoin="round"/>'+
    '<circle cx="'+xs(0).toFixed(1)+'" cy="'+ys(prof[0][1]).toFixed(1)+'" r="2.2" fill="var(--accent)"/>'+
    '<text x="'+pl+'" y="'+H+'" font-size="7.5" fill="#94a3b8">5&#39; pos</text>'+
    '<text x="'+(W-pr)+'" y="'+(pt-1)+'" font-size="7.5" fill="#94a3b8" text-anchor="end">'+(ymax*100).toFixed(0)+'%</text></svg>';
}
// ---- Functional gene burden (cohort HIGH+MODERATE snpEff burden per gene; optional --gene-burden aux TSV) ----
var GENEBURDEN_CAPTION='Genes ranked by cohort HIGH + MODERATE snpEff burden (summed across all samples); dominant effect = the most frequent effect class in that gene. HIGH + MODERATE counts are absolute, so a longer or more divergent gene ranks higher; this is a functional companion to the positional Variable genes panel below, not a selection test.';
function renderGeneBurden(){
  var host=el('gb_body'),tbl=el('gbtable'),note=el('gb_note'); if(!host||!tbl)return;
  var gb=R.gene_burden; if(!(gb&&gb.length)){host.style.display='none';return;} host.style.display='';
  var q=(st.gbq||'').toLowerCase(),exp=(el('geneburdenPanel')&&el('geneburdenPanel').classList.contains('expanded'));
  var flt=q?gb.filter(function(g){return (g.gene||'').toLowerCase().indexOf(q)>=0;}):gb;
  var rows=flt.slice(0,q?300:(exp?150:20));
  tbl.innerHTML='<thead><tr><th class="s" style="text-align:left">Gene</th><th>HIGH</th><th>MODERATE</th><th style="text-align:left">Dominant effect</th><th>Samples</th><th>Burden</th></tr></thead><tbody>'+
    rows.map(function(g){var mark=(g.start!=null&&g.end!=null);var burd=(g.total_impactful!=null?g.total_impactful:(g.high||0)+(g.moderate||0));
      return '<tr'+(mark?' data-start="'+g.start+'" data-end="'+g.end+'" data-name="'+esc(g.gene)+'" style="cursor:pointer"':'')+'>'+
        '<td class="s" style="text-align:left">'+esc(g.gene)+'</td>'+
        '<td'+((g.high||0)>0?' style="color:var(--fail);font-weight:600"':'')+'>'+(g.high||0)+'</td>'+
        '<td>'+(g.moderate||0)+'</td>'+
        '<td style="text-align:left">'+esc(g.dominant_effect||'NA')+'</td>'+
        '<td>'+(g.n_samples==null?'-':g.n_samples)+'</td>'+
        '<td>'+Math.round(burd).toLocaleString('en-US')+'</td></tr>';}).join('')+'</tbody>';
  var gbcount=q?('Showing '+rows.length+' of '+flt.length+' matching genes. '):(flt.length>rows.length?('Top '+rows.length+' of '+flt.length+' genes (search or expand for more). '):'');
  note.textContent=gbcount+GENEBURDEN_CAPTION+(rows.some(function(g){return g.start!=null;})?'':' Click-to-mark is available when gene coordinates are joined from a GFF upstream.');
  if(!rows.length)tbl.innerHTML='<tbody><tr><td class="na" colspan="6" style="text-align:left;padding:8px">no gene matches "'+esc(st.gbq||'')+'".</td></tr></tbody>';
  Array.prototype.forEach.call(tbl.querySelectorAll('tbody tr[data-start]'),function(tr){tr.onclick=function(){
    var gl=R.genome_len||0,nb=R.nbins||200; if(!gl)return;
    var b0=Math.max(0,Math.min(nb-1,Math.floor(+tr.getAttribute('data-start')/gl*nb)));
    var b1=Math.max(0,Math.min(nb-1,Math.floor(+tr.getAttribute('data-end')/gl*nb)));
    st.gtrack='snp'; st.geneMark={b0:b0,b1:b1,name:tr.getAttribute('data-name')||''};
    renderGenome(); var gs=el('genome');if(gs)gs.scrollIntoView();};});
}
// ---- Selection: cohort per-gene dN/dS from eskaks (optional --pnps; population genetics, NOT per-sample QC) ----
var PNPS_CAPTION='Cohort selection analysis (pairwise dN/dS across samples per gene, eskaks Nei/Li). This is population genetics, NOT per-sample QC: a gene here is not a flag on any one sample. dN/dS > 1 is a signal of positive / diversifying selection at cohort scale but is noisy per gene, sensitive to alignment and to the codon model, and saturates at low divergence. dN, dS and the ratio can be undefined at low divergence (shown as NA when dS is near zero). This panel appears only when a cohort pN/pS table is provided upstream; it is not derivable from one annotated VCF.';
function renderPnps(){
  var host=el('pnps_body'),tbl=el('pnpstable'),note=el('pnps_note'); if(!host||!tbl)return;
  var pp=R.pnps; if(!(pp&&pp.length)){host.style.display='none';return;} host.style.display='';
  var q=(st.pnpsq||'').toLowerCase(),exp=(el('pnpsPanel')&&el('pnpsPanel').classList.contains('expanded'));
  var flt=q?pp.filter(function(g){return (g.gene||'').toLowerCase().indexOf(q)>=0;}):pp;
  var rows=flt.slice(0,q?300:(exp?150:40));
  tbl.innerHTML='<thead><tr><th class="s" style="text-align:left">Gene</th><th>dN</th><th>dS</th><th>dN/dS</th><th>pairs</th><th style="text-align:left">Effect</th></tr></thead><tbody>'+
    rows.map(function(g){var pos=(g.pnps!=null&&g.pnps>1);
      return '<tr><td class="s" style="text-align:left">'+esc(g.gene)+'</td>'+
        '<td>'+(g.pn==null?'NA':g.pn.toFixed(3))+'</td>'+
        '<td>'+(g.ps==null?'NA':g.ps.toFixed(3))+'</td>'+
        '<td'+(pos?' style="color:var(--warn);font-weight:600"':'')+'>'+(g.pnps==null?'NA':g.pnps.toFixed(2))+'</td>'+
        '<td>'+(g.n==null?'-':Math.round(g.n))+'</td>'+
        '<td style="text-align:left">'+esc(g.effect||'NA')+'</td></tr>';}).join('')+'</tbody>';
  if(note)note.textContent=(q?('Showing '+rows.length+' of '+flt.length+' matching genes. '):'')+PNPS_CAPTION;
  if(!rows.length&&tbl)tbl.innerHTML='<tbody><tr><td class="na" colspan="6" style="text-align:left;padding:8px">no gene matches "'+esc(st.pnpsq||'')+'".</td></tr></tbody>';
}
function renderADNA(){
  var host=el('adna_body'); if(!host)return;
  var anc=R.samples.filter(function(s){return s.anc;});
  if(!R.n_ancient||!anc.length){host.innerHTML='<span class="nd" style="padding:0">no ancient (aDNA) samples in this run.</span>';return;}
  var lim=(athr&&athr.damage_min_ct!=null)?athr.damage_min_ct:null;
  var pctv=function(v){return v==null?'NA':(v*100).toFixed(1)+'%';};
  var cards=anc.map(function(s){var d=s.dmg,low=s.f.indexOf('DAMAGE_LOW')>=0;
    var ctcls=(d&&d.ct1!=null)?(low?'bad':'good'):'na';
    var head='<div class="acard-h"><span class="sname" data-s="'+esc(s.s)+'">'+esc(s.s)+'</span><span class="v '+s.v+'" style="font-size:10px">'+s.v+'</span></div>';
    if(!d)return '<div class="acard'+(low?' low':'')+'">'+head+'<div class="nd" style="padding:6px 0">damage data not found (run mapDamage2 / pass --mapdamage-dir)</div></div>';
    var body='<div class="acard-row"><div class="aspark">'+ctSpark(d.ct_profile,lim)+'</div><div class="astats">'+
      '<div class="astat"><span class="ak">5&#39; C&gt;T</span><span class="av '+ctcls+'">'+pctv(d.ct1)+'</span></div>'+
      '<div class="astat"><span class="ak">3&#39; G&gt;A</span><span class="av">'+pctv(d.ga1)+'</span></div>'+
      '<div class="astat"><span class="ak">frag len</span><span class="av">'+(d.fraglen==null?'NA':d.fraglen.toFixed(0)+' bp')+'</span></div></div></div>'+
      (low?'<div class="alow">⚠ terminal C&gt;T below the '+(lim*100).toFixed(0)+'% authentication floor - possible modern contamination; verify before using as a calibration tip.</div>':'');
    return '<div class="acard'+(low?' low':'')+'">'+head+body+'</div>';}).join('');
  host.innerHTML='<div class="anote">Elevated terminal C&gt;T (5&#39;) / G&gt;A (3&#39;) deamination is <b>consistent with</b> post-mortem damage and is a necessary authentication signal. It does <b>not</b> by itself prove the DNA is ancient (deaminated modern contaminant DNA can mimic it) or endogenous, and its absence (below the floor) is a red flag for a modern sample mislabelled ancient. Treat this as a screen, not a proof.</div>'+
    '<div class="acards">'+cards+'</div>';
  Array.prototype.forEach.call(host.querySelectorAll('.sname'),function(sp){sp.onclick=function(){openDetail(sp.getAttribute('data-s'));};});
}

function renderDynamics(){
  var host=el('dyn_body'), sec=el('dynamics'); if(!host)return;
  var D=R.dynamics;
  if(!(D&&D.groups&&D.groups.length)){ if(sec)sec.style.display='none'; var nv=el('nav-dyn'); if(nv)nv.style.display='none'; return; }
  if(sec)sec.style.display='';
  var th=D.thresholds||{emerge:0.25,fix:0.9,loss:0.1};
  var FCOL={fixation:'#2f6fed',emergence:'#22a06b',loss:'#e6893a',nonsyn:'#d1495b',high_impact:'#7c3aed'};
  function serColor(f){ if(f.indexOf('fixation')>=0)return FCOL.fixation; if(f.indexOf('emergence')>=0)return FCOL.emergence; if(f.indexOf('loss')>=0)return FCOL.loss; return '#c3ccda'; }
  var MAXSER=60;
  var body=D.groups.map(function(g){
    var n=g.times.length, W=680,H=250, ml=42,mr=14,mt=12,mb=34, pw=W-ml-mr, ph=H-mt-mb;
    function X(i){ return ml + (n<=1? pw/2 : (i/(n-1))*pw); }
    function Y(a){ return mt + (1-a)*ph; }
    var svg='<svg viewBox="0 0 '+W+' '+H+'" width="100%" style="max-width:'+W+'px;display:block;font:11px system-ui">';
    [0,0.25,0.5,0.75,1].forEach(function(a){ svg+='<line x1="'+ml+'" y1="'+Y(a).toFixed(1)+'" x2="'+(W-mr)+'" y2="'+Y(a).toFixed(1)+'" stroke="#eef2f7"/><text x="'+(ml-6)+'" y="'+(Y(a)+3).toFixed(1)+'" text-anchor="end" fill="#8a97a8">'+a.toFixed(2)+'</text>'; });
    [[th.emerge,FCOL.emergence],[th.fix,FCOL.fixation],[th.loss,FCOL.loss]].forEach(function(t){ svg+='<line x1="'+ml+'" y1="'+Y(t[0]).toFixed(1)+'" x2="'+(W-mr)+'" y2="'+Y(t[0]).toFixed(1)+'" stroke="'+t[1]+'" stroke-dasharray="3 3" opacity="0.35"/>'; });
    g.times.forEach(function(t,i){ svg+='<text x="'+X(i).toFixed(1)+'" y="'+(H-14)+'" text-anchor="middle" fill="#4a5768">'+esc(t==null?i:t)+'</text>'; });
    svg+='<text x="'+ml+'" y="'+(H-1)+'" fill="#8a97a8" font-size="10">time / passage &#8594;</text>';
    g.series.slice(0,MAXSER).forEach(function(s){
      var col=serColor(s.flags), flagged=s.flags.length>0, nonsyn=s.flags.indexOf('nonsyn')>=0;
      var pts=s.traj.map(function(a,i){return X(i).toFixed(1)+','+Y(a).toFixed(1);}).join(' ');
      var tip=esc(s.pos+(s.gene?(' '+s.gene):'')+(s.eff?(' '+s.eff):'')+(s.flags.length?(' ['+s.flags.join(',')+']'):''));
      svg+='<polyline points="'+pts+'" fill="none" stroke="'+col+'" stroke-width="'+(flagged?2:1)+'" opacity="'+(flagged?0.95:0.4)+'"><title>'+tip+'</title></polyline>';
      s.traj.forEach(function(a,i){ svg+='<circle cx="'+X(i).toFixed(1)+'" cy="'+Y(a).toFixed(1)+'" r="'+(flagged?3:2)+'" fill="'+col+'" stroke="'+(nonsyn?FCOL.nonsyn:'#fff')+'" stroke-width="'+(nonsyn?1.6:0.6)+'"><title>'+tip+' | AF='+a.toFixed(2)+'</title></circle>'; });
    });
    svg+='</svg>';
    var flagged=g.series.filter(function(s){return s.flags.length;});
    var tbl='';
    if(flagged.length){
      tbl='<table class="dyntbl"><thead><tr><th>position</th><th>gene</th><th>effect</th><th>events</th><th>trajectory</th></tr></thead><tbody>'+
        flagged.slice(0,40).map(function(s){
          var chips=s.flags.map(function(f){return '<span class="dchip" style="background:'+(FCOL[f]||'#8895a6')+'">'+f+'</span>';}).join(' ');
          return '<tr><td>'+esc(s.pos)+'</td><td>'+esc(s.gene||'-')+'</td><td>'+esc(s.eff||'-')+'</td><td>'+chips+'</td><td class="dtraj">'+s.traj.map(function(a){return a.toFixed(2);}).join(' &#8594; ')+'</td></tr>';
        }).join('')+'</tbody></table>';
    }
    var extra=g.series.length>MAXSER?('<span class="c"> (showing '+MAXSER+' most dynamic of '+g.series.length+')</span>'):'';
    return '<div class="dyngrp"><div class="dynhd"><b>'+esc(g.group)+'</b> <span class="c">- '+g.samples.length+' samples, '+g.n_flagged+' flagged SNP(s)'+extra+'</span></div>'+svg+tbl+'</div>';
  }).join('');
  var legend='<div class="dynleg"><span><i style="background:'+FCOL.emergence+'"></i>emergence</span><span><i style="background:'+FCOL.fixation+'"></i>fixation</span><span><i style="background:'+FCOL.loss+'"></i>loss</span><span><i style="border:2px solid '+FCOL.nonsyn+';background:#fff"></i>non-synonymous</span><span class="c" style="margin-left:auto">dashed guides = thresholds; each line = one SNP</span></div>';
  host.innerHTML=legend+body;
}

function renderAll(){renderOverview();renderTable();renderLineages();renderPlots();renderScatter();renderCorr();renderQCspace();renderRefBias();renderStacks();renderGenome();renderFunction();renderGeneBurden();renderHotspots();renderTemporal();renderPnps();renderADNA();renderDynamics();renderFlags();renderCuration();}

// ---- static wiring ----
el('meta').textContent=R.samples.length+' samples · '+R.generated;
el('foot').innerHTML='Generated '+R.generated+' · thresholds are adjustable live above; the pipeline gate uses the defaults ('+
  Object.keys(R.thresholds).map(function(k){return k+'='+R.thresholds[k];}).join(', ')+'). Values scale within each column; NA = not reported.';
el('colmenu').innerHTML=R.metrics.map(function(m){return '<label><input type="checkbox" data-k="'+m.key+'"'+(st.hidden[m.key]?'':' checked')+'> '+esc(m.label)+'</label>';}).join('');
Array.prototype.forEach.call(document.querySelectorAll('#colmenu input'),function(cb){cb.onchange=function(){if(cb.checked)delete st.hidden[cb.getAttribute('data-k')];else st.hidden[cb.getAttribute('data-k')]=1;renderTable();};});
el('q').oninput=function(e){st.q=e.target.value.toLowerCase().trim();renderTable();renderPlots();renderScatter();renderCorr();renderQCspace();renderRefBias();renderStacks();renderGenome();renderFunction();renderTemporal();};
// per-panel gene search (Functional gene burden / Variable genes / pN-pS): filter each gene table by gene name
[['gbq','gbq',renderGeneBurden],['hotq','hotq',renderHotspots],['pnpsq','pnpsq',renderPnps]].forEach(function(w){var inp=el(w[0]);if(inp)inp.oninput=function(e){st[w[1]]=e.target.value.trim();w[2]();};});
el('of').onchange=function(e){st.onlyFlagged=e.target.checked;renderAll();};
Array.prototype.forEach.call(document.querySelectorAll('#ptype button'),function(b){b.onclick=function(){st.ptype=b.getAttribute('data-t');
  Array.prototype.forEach.call(document.querySelectorAll('#ptype button'),function(x){x.classList.toggle('on',x==b);});renderPlots();};});
// scatter axis pickers
['sx','sy'].forEach(function(ax){var sel=el(ax);sel.innerHTML=R.metrics.map(function(m){return '<option value="'+m.key+'"'+(st[ax]==m.key?' selected':'')+'>'+esc(m.label)+'</option>';}).join('');
  sel.onchange=function(){st[ax]=sel.value;renderScatter();};});
// live thresholds + presets
var TH=[['depth_min','Depth min'],['breadth_min','Breadth min %'],['missing_max','Missing max %'],['mapping_min','Mapped min %'],['dup_max','Dup max %'],['iupac_max','IUPAC max %'],['titv_min','Ti/Tv min'],['snp_z','SNP z'],['het_max_frac','Het % max'],['mixed_min_frac','Mixed lin % min']];
var PRESETS=[
  {id:'gate',label:'gate defaults',th:Object.assign({},R.thresholds),note:'The cut-offs the Snakemake qc_gate uses (config report_* keys).'},
  {id:'strict',label:'strict (modern WGS)',th:{depth_min:20,breadth_min:95,missing_max:5,mapping_min:90,dup_max:30,iupac_max:2,titv_min:1.5,snp_z:3,het_max_frac:1.5,mixed_min_frac:2},note:'Confident modern Illumina isolate: 20x, 95% breadth, <5% missing, Ti/Tv >=1.5.'},
  {id:'lenient',label:'lenient (aDNA / low-cov)',th:{depth_min:3,breadth_min:60,missing_max:40,mapping_min:50,dup_max:80,iupac_max:5,titv_min:0,snp_z:4,het_max_frac:8,mixed_min_frac:5},note:'Degraded / low-coverage library at the 3x calling floor: rescues calibration tips.'}
];
function applyThr(next){Object.keys(next).forEach(function(k){if(k in thr)thr[k]=next[k];});
  Array.prototype.forEach.call(document.querySelectorAll('#thbox input'),function(inp){var k=inp.getAttribute('data-t');if(k in thr)inp.value=thr[k];});
  recompute();renderAll();saveState();}
el('thbox').innerHTML='<label style="display:inline-flex;flex-direction:column;font-size:10px;color:#64748b;gap:2px">preset<select id="thpreset" class="msel" style="padding:4px 6px"><option value="">custom…</option>'+
  PRESETS.map(function(p){return '<option value="'+p.id+'">'+esc(p.label)+'</option>';}).join('')+'</select></label>'+
  '<span id="thnote" style="flex-basis:100%;font-size:10.5px;color:#8895a6;margin-top:-2px"></span>'+
  TH.map(function(t){return '<label style="display:inline-flex;flex-direction:column;font-size:10px;color:#64748b;gap:2px">'+t[1]+
  '<input type="number" step="any" data-t="'+t[0]+'" value="'+thr[t[0]]+'" style="width:78px;padding:4px 6px;border:1px solid var(--line);border-radius:6px;font-size:12px"></label>';}).join('')+
  '<button class="btn" id="threset" style="align-self:flex-end">reset</button>';
var thpre=el('thpreset'),thnote=el('thnote');
if(thpre)thpre.onchange=function(){var p=null;PRESETS.forEach(function(x){if(x.id==thpre.value)p=x;});if(!p){thnote.textContent='';return;}thnote.textContent=p.note;applyThr(p.th);};
Array.prototype.forEach.call(document.querySelectorAll('#thbox input'),function(inp){inp.oninput=function(){var v=parseFloat(inp.value);if(!isNaN(v)){thr[inp.getAttribute('data-t')]=v;if(thpre)thpre.value='';if(thnote)thnote.textContent='';recompute();renderAll();saveState();}};});
el('threset').onclick=function(){if(thpre)thpre.value='';if(thnote)thnote.textContent='';applyThr(Object.assign({},R.thresholds));};
// ancient (aDNA) live thresholds
if(R.n_ancient){var ATH=[['depth_min','aDNA depth min'],['breadth_min','aDNA breadth %'],['missing_max','aDNA missing %'],['mapping_min','aDNA mapped %'],['dup_max','aDNA dup %'],['iupac_max','aDNA IUPAC %'],['damage_min_ct',"5′ C>T min (0-1)"]];
  el('athbox').innerHTML='<div style="flex-basis:100%;font-size:10px;color:#8a5a12;font-weight:600;text-transform:uppercase;letter-spacing:.06em">Ancient (aDNA) thresholds &middot; a 5x mummy is judged here, not against the modern gate</div>'+
    ATH.map(function(t){return '<label style="display:inline-flex;flex-direction:column;font-size:10px;color:#64748b;gap:2px">'+t[1]+
    '<input type="number" step="any" data-t="'+t[0]+'" value="'+(athr[t[0]]!=null?athr[t[0]]:'')+'" style="width:78px;padding:4px 6px;border:1px solid var(--line);border-radius:6px;font-size:12px"></label>';}).join('');
  Array.prototype.forEach.call(document.querySelectorAll('#athbox input'),function(inp){inp.oninput=function(){var v=parseFloat(inp.value);if(!isNaN(v)){athr[inp.getAttribute('data-t')]=v;recompute();renderAll();saveState();}};});}
// CSV export of the current (filtered rows, visible columns) table
el('csv').onclick=function(){var mets=R.metrics.filter(function(m){return !st.hidden[m.key];});
  var head=['sample','verdict'].concat(mets.map(function(m){return m.key;})).concat(['lineage','flags']);
  var lines=[head.join('\t')]; visible().forEach(function(s){lines.push([s.s,s.v].concat(mets.map(function(m){return s.m[m.key]==null?'':s.m[m.key];})).concat([s.lineage||'',s.f.join(';')]).join('\t'));});
  var blob=new Blob([lines.join('\n')],{type:'text/tab-separated-values'}),a=document.createElement('a');
  a.href=URL.createObjectURL(blob);a.download='qc_table.tsv';a.click();URL.revokeObjectURL(a.href);};
// tooltips + click on plots/scatter
function bandTip(pk,val){var b=bandFor(pk,thr);if(!b)return'';var mk=(MET[pk]||{}).kind,inb=val>=b[0]&&val<=b[1];
  var lab=(b[0]==-Infinity)?('≤ '+shortv(b[1],mk)):(b[1]==Infinity)?('≥ '+shortv(b[0],mk)):(shortv(b[0],mk)+' to '+shortv(b[1],mk));
  return '<br><span style="color:'+(inb?'#8fe3c0':'#f2b8bc')+'">band '+lab+(inb?' ✓ in':' ✗ out')+'</span>';}
function wireHover(host){host.addEventListener('mousemove',function(e){var t=e.target,lin;
  if(t.tagName=='circle'&&t.hasAttribute('data-val')){lin=t.getAttribute('data-lin');var pk=t.getAttribute('data-pk'),vv=+t.getAttribute('data-val');
    tip('<b>'+esc(t.getAttribute('data-s'))+'</b>'+(lin?' · '+esc(lin):'')+'<br>'+esc(t.getAttribute('data-lab'))+': '+shortv(vv,t.getAttribute('data-kind'))+(pk?bandTip(pk,vv):''),e.clientX,e.clientY);}
  else if(t.tagName=='rect'&&t.hasAttribute('data-val')){lin=t.getAttribute('data-lin');var pk2=t.getAttribute('data-pk'),vv2=+t.getAttribute('data-val');
    tip('<b>'+esc(t.getAttribute('data-s'))+'</b>'+(lin?' · '+esc(lin):'')+'<br>'+esc(t.getAttribute('data-lab'))+': '+shortv(vv2,t.getAttribute('data-kind'))+(pk2?bandTip(pk2,vv2):''),e.clientX,e.clientY);}
  else if(t.tagName=='circle'&&t.hasAttribute('data-x'))tip('<b>'+esc(t.getAttribute('data-s'))+'</b>'+(t.getAttribute('data-lin')?' · '+esc(t.getAttribute('data-lin')):'')+'<br>'+esc(t.getAttribute('data-xl'))+': '+shortv(+t.getAttribute('data-x'),t.getAttribute('data-xk'))+'<br>'+esc(t.getAttribute('data-yl'))+': '+shortv(+t.getAttribute('data-y'),t.getAttribute('data-yk')),e.clientX,e.clientY);
  else tip('');});
  host.addEventListener('mouseleave',function(){tip('');});
  host.addEventListener('click',function(e){if((e.target.tagName=='circle'||e.target.tagName=='rect')&&e.target.getAttribute('data-s'))setHi(e.target.getAttribute('data-s'));});}
wireHover(el('plots')); wireHover(el('scatter')); wireHover(el('qcpca_body')); wireHover(el('divcomp_body')); wireHover(el('temporal_body'));
// correlation matrix: hover tip + click a cell -> load that metric pair into the scatter
(function(){var c=el('corr_body'); if(!c)return;
  c.addEventListener('mousemove',function(e){var t=e.target;
    if(t.tagName=='rect'&&t.hasAttribute('data-r')&&t.getAttribute('data-xk')){var r=t.getAttribute('data-r');
      tip('<b>'+esc((MET[t.getAttribute('data-yk')]||{}).label||'')+'</b> vs <b>'+esc((MET[t.getAttribute('data-xk')]||{}).label||'')+'</b><br>Spearman rho '+(r===''?'NA':r),e.clientX,e.clientY);}
    else tip('');});
  c.addEventListener('mouseleave',function(){tip('');});
  c.addEventListener('click',function(e){var t=e.target;
    if(t.tagName=='rect'&&t.getAttribute('data-xk')&&t.getAttribute('data-xk')!=t.getAttribute('data-yk')){
      st.sx=t.getAttribute('data-xk');st.sy=t.getAttribute('data-yk');
      var sxs=el('sx'),sys=el('sy'); if(sxs)sxs.value=st.sx; if(sys)sys.value=st.sy;
      renderScatter();el('corr').scrollIntoView({behavior:'smooth'});}});})();
// scatter rubber-band select -> add the enclosed samples to the exclusion basket (rAF-gated overlay, hit-test on release)
(function(){var host=el('scatter'); if(!host)return;
  var dragging=false,rectL=0,rectT=0,x0=0,y0=0,x1=0,y1=0,raf=0;
  function paint(){raf=0;var ov=el('scbrush');if(!ov)return;
    ov.setAttribute('x',Math.min(x0,x1));ov.setAttribute('y',Math.min(y0,y1));ov.setAttribute('width',Math.abs(x1-x0));ov.setAttribute('height',Math.abs(y1-y0));ov.style.display='';}
  function move(e){x1=e.clientX-rectL;y1=e.clientY-rectT;if(!raf)raf=requestAnimationFrame(paint);}
  function up(){if(!dragging)return;dragging=false;
    document.removeEventListener('mousemove',move);document.removeEventListener('mouseup',up);
    var ov=el('scbrush');if(ov)ov.style.display='none';
    if(!SGEO||(Math.abs(x1-x0)<4&&Math.abs(y1-y0)<4))return;   // a click, not a drag
    var G=SGEO,l=Math.min(x0,x1),rr=Math.max(x0,x1),tp=Math.min(y0,y1),bt=Math.max(y0,y1);
    function sx(v){return G.pad+(G.xr[1]>G.xr[0]?(v-G.xr[0])/(G.xr[1]-G.xr[0]):0.5)*G.plot;}
    function sy(v){return G.H-G.pad-(G.yr[1]>G.yr[0]?(v-G.yr[0])/(G.yr[1]-G.yr[0]):0.5)*G.ph;}
    function inbox(s){var xv=s.m[G.xk],yv=s.m[G.yk];if(xv==null||yv==null)return false;
      var px=sx(xv),py=sy(yv);return px>=l&&px<=rr&&py>=tp&&py<=bt;}
    var vis=visible(),visSet={}; vis.forEach(function(s){visSet[s.s]=1;});
    var picked=vis.filter(inbox);
    var skipped=R.samples.filter(function(s){return !visSet[s.s]&&inbox(s);}).length;   // dimmed, out-of-filter dots inside the box
    var dim=skipped?' <span style="color:#94a3b8">('+skipped+' filtered-out dot'+(skipped>1?'s':'')+' ignored)</span>':'';
    var rd=el('scbrushinfo'); if(!picked.length){if(rd)rd.innerHTML=skipped?'no in-view samples in the box'+dim:'';return;}
    var added=picked.filter(function(s){return !st.excl[s.s];});   // undo only reverts what this brush added
    added.forEach(function(s){st.excl[s.s]=1;}); renderTable();renderCuration();
    if(rd){rd.innerHTML='basketed <b>'+added.length+'</b> of '+picked.length+' &#8594; exclusion'+dim+(added.length?' <button class="btn" id="scbrushundo" style="padding:2px 8px;font-size:11px">undo</button>':'');
      var ub=el('scbrushundo'); if(ub)ub.onclick=function(){added.forEach(function(s){delete st.excl[s.s];});renderTable();renderCuration();rd.innerHTML='';};}}
  host.addEventListener('mousedown',function(e){if(e.button!==0)return;var svg=el('scsvg');if(!svg)return;
    var rc=svg.getBoundingClientRect();rectL=rc.left;rectT=rc.top;
    dragging=true;x0=x1=e.clientX-rectL;y0=y1=e.clientY-rectT;
    document.addEventListener('mousemove',move);document.addEventListener('mouseup',up);e.preventDefault();});})();
// genome landscape hover + click
(function(){var g=el('genome_body'); if(!g)return;
  var LAB={missing:'missing',snp:'homozygous SNPs',het:'het variants',indel:'indels'};
  g.addEventListener('mousemove',function(e){var t=e.target;
    if(t.tagName=='rect'&&t.hasAttribute('data-v')){var bin=+t.getAttribute('data-bin'),nb=R.nbins||200,gl=R.genome_len||nb,tk=st.gtrack||'missing';
      var p0=Math.round(bin/nb*gl),p1=Math.round((bin+1)/nb*gl),v=t.getAttribute('data-v');
      var val=(v==='')?'NA':(tk=='missing'?v+'% missing':v+' '+LAB[tk]);
      tip('<b>'+esc(t.getAttribute('data-s')||'cohort (mean)')+'</b><br>'+fmtpos(p0)+' - '+fmtpos(p1)+'<br>'+val,e.clientX,e.clientY);}
    else tip('');});
  g.addEventListener('mouseleave',function(){tip('');});
  g.addEventListener('click',function(e){if(e.target.getAttribute('data-s'))setHi(e.target.getAttribute('data-s'));});})();
// genome brush: drag to ZOOM the plot into a reference span (and filter the Variable-genes table); nested drags zoom further
function genomeReadout(b0,b1,zoomed){var nb=R.nbins||200,gl=R.genome_len||nb;
  var p0=Math.round(b0/nb*gl),p1=Math.round((b1+1)/nb*gl);
  var ng=(R.genes||[]).filter(function(ge){return ge.end>=p0&&ge.start<=p1;}).length;
  var rd=el('gselreadout');if(rd)rd.innerHTML=(zoomed?'zoomed ':'')+'<b>'+fmtpos(p0)+' - '+fmtpos(p1)+'</b> &middot; '+(p1-p0).toLocaleString('en-US')+' bp'+(R.genes&&R.genes.length?' &middot; '+ng+' gene'+(ng==1?'':'s'):'')+' <button class="btn" id="gselclear" style="padding:2px 8px;font-size:11px">'+(zoomed?'reset zoom':'release to zoom')+'</button>';}
function genomeResetZoom(){st.gsel=null;st.gzoom=null;st.geneMark=null;var gg=el('genegoto');if(gg)gg.value='';renderGenome();renderHotspots();var rd=el('gselreadout');if(rd)rd.innerHTML='';}
(function(){var g=el('genome_body'); if(!g)return;
  var dragging=false,rectL=0,startBin=0,curBin=0,raf=0;
  function binAt(cx){if(!GGEO)return 0;var b=GGEO.z0+Math.floor((cx-rectL-GGEO.gut)/GGEO.plotW*GGEO.winN);return Math.max(GGEO.z0,Math.min(GGEO.z0+GGEO.winN-1,b));}
  function paint(){raf=0;var gb=el('gbrush');if(!gb||!GGEO)return;
    var b0=Math.min(startBin,curBin),b1=Math.max(startBin,curBin);
    var x0=GGEO.gut+(b0-GGEO.z0)/GGEO.winN*GGEO.plotW,x1=GGEO.gut+(b1+1-GGEO.z0)/GGEO.winN*GGEO.plotW;
    gb.setAttribute('x',x0.toFixed(1));gb.setAttribute('width',(x1-x0).toFixed(1));gb.style.display='';
    genomeReadout(b0,b1,false);}
  function move(e){curBin=binAt(e.clientX);if(!raf)raf=requestAnimationFrame(paint);}
  function up(){if(!dragging)return;dragging=false;
    document.removeEventListener('mousemove',move);document.removeEventListener('mouseup',up);
    var b0=Math.min(startBin,curBin),b1=Math.max(startBin,curBin);
    if(b1>b0){st.gsel={b0:b0,b1:b1};st.gzoom={b0:b0,b1:b1};renderGenome();renderHotspots();genomeReadout(b0,b1,true);
      var cb=el('gselclear');if(cb)cb.onclick=function(ev){ev.stopPropagation();genomeResetZoom();};}}
  g.addEventListener('mousedown',function(e){if(e.button!==0)return;var svg=g.querySelector('svg');if(!svg)return;
    rectL=svg.getBoundingClientRect().left;dragging=true;startBin=curBin=binAt(e.clientX);
    document.addEventListener('mousemove',move);document.addEventListener('mouseup',up);e.preventDefault();});})();
// go-to-gene: type a gene name -> zoom the plot to that gene and mark it (needs a GFF, i.e. R.genes)
(function(){var gg=el('genegoto');if(!gg)return;
  if(!(R.genes&&R.genes.length)){gg.style.display='none';return;}
  gg.oninput=function(e){var q=e.target.value.trim().toLowerCase();
    if(!q){genomeResetZoom();return;}
    var g=null,i;for(i=0;i<R.genes.length;i++){if((R.genes[i].name||'').toLowerCase().indexOf(q)>=0){g=R.genes[i];break;}}
    if(!g)return;
    var nb=R.nbins||200,gl=R.genome_len||nb;
    var b0=Math.max(0,Math.min(nb-1,Math.floor(g.start/gl*nb))),b1=Math.max(0,Math.min(nb-1,Math.floor(g.end/gl*nb)));
    var pad=Math.max(3,Math.round((b1-b0)*0.6)+2);
    st.gzoom={b0:Math.max(0,b0-pad),b1:Math.min(nb-1,b1+pad)}; st.geneMark={b0:b0,b1:b1,name:g.name};
    if(gtrackHas('snp')){st.gtrack='snp';Array.prototype.forEach.call(document.querySelectorAll('#gtrack button'),function(b){b.classList.toggle('on',b.getAttribute('data-gt')=='snp');});}
    renderGenome();renderHotspots();
    var rd=el('gselreadout');if(rd)rd.innerHTML='gene <b>'+esc(g.name)+'</b> &middot; '+fmtpos(g.start)+' - '+fmtpos(g.end)+' <button class="btn" id="gselclear" style="padding:2px 8px;font-size:11px">reset zoom</button>';
    var cb=el('gselclear');if(cb)cb.onclick=function(ev){ev.stopPropagation();genomeResetZoom();};};})();
if(!R.samples.some(function(s){return s.miss||s.trk;})){var gs=el('genome');if(gs)gs.style.display='none';var ng=el('nav-genome');if(ng)ng.style.display='none';}
if(!R.samples.some(function(s){return s.trk&&s.trk.snp;})){var hsx=el('hotspots');if(hsx)hsx.style.display='none';var nhx=el('nav-hot');if(nhx)nhx.style.display='none';}
// signature layer gates: hide a panel when its data type is absent cohort-wide (organism-agnostic; render fns also degrade per-view)
if(!pca2(R.samples).ok){var qp=el('qcpca');if(qp)qp.style.display='none';var nqp=el('nav-pca');if(nqp)nqp.style.display='none';}
if(!R.samples.some(function(s){return s.m.snps!=null||s.m.snp_density!=null;})){var dc=el('divcomp');if(dc)dc.style.display='none';var ndc=el('nav-divcomp');if(ndc)ndc.style.display='none';}
if(!R.samples.some(function(s){return yearOf(s.date)!=null;})){var tpx=el('temporal');if(tpx)tpx.style.display='none';var ntpx=el('nav-temporal');if(ntpx)ntpx.style.display='none';}
if(!R.samples.some(function(s){return s.m.ann_high!=null||s.m.ann_moderate!=null||s.m.ann_modifier!=null;})){var fnx=el('function');if(fnx)fnx.style.display='none';var nfx=el('nav-function');if(nfx)nfx.style.display='none';}
if(!(R.gene_burden&&R.gene_burden.length)){var gbx=el('geneburden');if(gbx)gbx.style.display='none';var ngb=el('nav-geneburden');if(ngb)ngb.style.display='none';}
if(!(R.pnps&&R.pnps.length)){var ppx=el('pnps');if(ppx)ppx.style.display='none';var npp=el('nav-pnps');if(npp)npp.style.display='none';}
// genome track selector (Missing / SNPs / Het / Indels) - only offer tracks that have data
(function(){var host=el('gtrack'); if(!host)return;var avail=GTRACKS.filter(function(g){return gtrackHas(g.k);});
  if(avail.length<=1){host.style.display='none';return;}
  host.innerHTML=avail.map(function(g){return '<button'+(g.k==st.gtrack?' class="on"':'')+' data-gt="'+g.k+'">'+g.lab+'</button>';}).join('');
  Array.prototype.forEach.call(host.querySelectorAll('button'),function(b){b.onclick=function(){st.gtrack=b.getAttribute('data-gt');
    Array.prototype.forEach.call(host.querySelectorAll('button'),function(x){x.classList.toggle('on',x==b);});renderGenome();};});})();
// mask-regions toggle (mtbc_mask etc.): grey the masked zones + exclude them from the SNP density + variable-gene ranking
(function(){var mb=el('maskbtn'); if(!mb)return;
  if(!R.mask_bins){mb.style.display='none';return;}
  mb.title=(R.mask_pct||0)+'% of the reference masked (PE/PPE, IS, DR, repeats); toggle to exclude these zones';
  mb.onclick=function(){st.maskOn=!st.maskOn;mb.classList.toggle('on',st.maskOn);renderGenome();renderHotspots();};})();
var rz;window.addEventListener('resize',function(){clearTimeout(rz);rz=setTimeout(function(){renderPlots();renderScatter();renderCorr();renderQCspace();renderRefBias();renderGenome();renderTemporal();},120);});
// metric help panel
el('helpmenu').innerHTML=R.metrics.map(function(m){var d=R.defs[m.key]||['',''];return '<div class="hitem"><b>'+esc(m.label)+'</b> <span class="hk">'+esc(m.key)+'</span><div class="hd">'+esc(d[0]||'')+(d[1]?' <span class="hr">('+esc(d[1])+')</span>':'')+'</div></div>';}).join('');
// ---- colour-by toggle (beeswarm + scatter); hides the whole lineage UI when there is no lineage data ----
(function(){var cb=el('colorby');
  if(!R.lin_present){if(cb)cb.style.display='none';['leg-lin','linsum','nav-lin','groupui'].forEach(function(id){var e=el(id);if(e)e.style.display='none';});return;}
  if(cb){Array.prototype.forEach.call(cb.querySelectorAll('button'),function(b){b.onclick=function(){st.colorBy=b.getAttribute('data-cb');
    Array.prototype.forEach.call(cb.querySelectorAll('button'),function(x){x.classList.toggle('on',x==b);});
    el('leg-qc').style.display=(st.colorBy=='lineage')?'none':'flex';el('leg-lin').style.display=(st.colorBy=='lineage')?'flex':'none';
    var pc=el('pcacb');if(pc)Array.prototype.forEach.call(pc.querySelectorAll('button'),function(x){x.classList.toggle('on',x.getAttribute('data-cb')==st.colorBy);});
    renderPlots();renderScatter();renderQCspace();};});}
  el('leg-lin').innerHTML=R.lineages.map(function(l){return '<span><i style="background:'+linColor(l)+'"></i>'+esc(l)+'</span>';}).join('')+'<span style="margin-left:auto">beeswarm dashed line = median</span>';
})();
// PCA colour-by (mirrors #colorby into st.colorBy; works even with no lineage data by dropping only the lineage button)
(function(){var cb=el('pcacb');if(!cb)return;
  if(!R.lin_present){var lb=cb.querySelector('[data-cb="lineage"]');if(lb)lb.style.display='none';}
  Array.prototype.forEach.call(cb.querySelectorAll('button'),function(b){b.onclick=function(){st.colorBy=b.getAttribute('data-cb');
    Array.prototype.forEach.call(cb.querySelectorAll('button'),function(x){x.classList.toggle('on',x==b);});
    var main=el('colorby');if(main)Array.prototype.forEach.call(main.querySelectorAll('button'),function(x){x.classList.toggle('on',x.getAttribute('data-cb')==st.colorBy);});
    var lq=el('leg-qc'),ll=el('leg-lin');if(lq)lq.style.display=(st.colorBy=='lineage')?'none':'flex';if(ll)ll.style.display=(st.colorBy=='lineage')?'flex':'none';
    renderPlots();renderScatter();renderQCspace();};});})();
// group-by-lineage checkbox
var glin=el('glin'); if(glin)glin.onchange=function(){st.groupLin=glin.checked;renderTable();};
// ancient / modern filter
if(R.n_ancient){el('ancfilter').innerHTML='<span class="seg" id="ancseg"><button class="on" data-a="">all</button><button data-a="mod">modern</button><button data-a="anc">aDNA <span class="k">'+R.n_ancient+'</span></button></span>';
  Array.prototype.forEach.call(document.querySelectorAll('#ancseg button'),function(b){b.onclick=function(){st.ancOnly=b.getAttribute('data-a')||null;
    Array.prototype.forEach.call(document.querySelectorAll('#ancseg button'),function(x){x.classList.toggle('on',x==b);});renderAll();};});}
// provenance / run-manifest header (self-documenting for a citable exclusion set)
(function(){var p=R.provenance||{},items=[];
  items.push('reference '+(p.reference||'NA')+(R.genome_len?' ('+R.genome_len.toLocaleString('en-US')+' bp)':''));
  if(p.container&&p.container!='none')items.push('container '+p.container);
  if(p.commit)items.push('commit '+p.commit);
  items.push(R.samples.length+' samples'+(R.n_ancient?' · '+R.n_ancient+' aDNA':''));
  var pe=el('prov'); if(pe)pe.innerHTML=items.map(function(t){return '<span>'+esc(t)+'</span>';}).join('');})();
// print
var pbtn=el('printBtn'); if(pbtn)pbtn.onclick=function(){window.print();};
// ---- persistence of the curated view (basket + thresholds), namespaced per sample-set so two reports don't bleed ----
var SKEY='sbampiro_qc_v2:'+R.samples.length+':'+(R.samples[0]?R.samples[0].s:'')+':'+(R.samples.length?R.samples[R.samples.length-1].s:'');
function saveState(){try{localStorage.setItem(SKEY,JSON.stringify({excl:st.excl,thr:thr,athr:athr}));}catch(e){}}
function loadState(){try{var s=JSON.parse(localStorage.getItem(SKEY)||'null');if(!s)return false;
  if(s.thr)Object.keys(s.thr).forEach(function(k){if(k in thr)thr[k]=s.thr[k];});
  if(s.athr)Object.keys(s.athr).forEach(function(k){if(k in athr)athr[k]=s.athr[k];});
  if(s.excl&&typeof s.excl=='object'){var have={};R.samples.forEach(function(x){have[x.s]=1;});
    st.excl={};Object.keys(s.excl).forEach(function(k){if(have[k])st.excl[k]=1;});}  // drop unknown sample ids
  return true;}catch(e){return false;}}
// ---- expand-to-fill (fullscreen within the window) for the big panels ----
function collapseExpanded(){var ex=document.querySelector('.panel.expanded');if(!ex)return;ex.classList.remove('expanded');document.body.classList.remove('has-expanded');
  Array.prototype.forEach.call(document.querySelectorAll('.exp-h'),function(b){b.textContent='⤢ full';});
  renderGenome();renderPlots();renderScatter();renderTable();renderQCspace();renderRefBias();renderFunction();renderGeneBurden();renderHotspots();renderPnps();}
Array.prototype.forEach.call(document.querySelectorAll('.exp-h'),function(b){b.onclick=function(){
  var panel=el(b.getAttribute('data-panel')); if(!panel)return;
  var willExpand=!panel.classList.contains('expanded'); collapseExpanded();
  if(willExpand){panel.classList.add('expanded');document.body.classList.add('has-expanded');b.textContent='⤡ close';}
  var rn=b.getAttribute('data-render');
  setTimeout(function(){if(rn=='genome')renderGenome();else if(rn=='plots')renderPlots();else if(rn=='scatter')renderScatter();else if(rn=='corr')renderCorr();else if(rn=='pca')renderQCspace();else if(rn=='divcomp')renderRefBias();else if(rn=='function')renderFunction();else if(rn=='geneburden')renderGeneBurden();else if(rn=='hotspots')renderHotspots();else if(rn=='pnps')renderPnps();else if(rn=='table')renderTable();},20);};});
if(el('expClose'))el('expClose').onclick=collapseExpanded;
document.addEventListener('keydown',function(e){if(e.key=='Escape')collapseExpanded();});
// per-sample detail modal close (button, backdrop, Esc)
el('modalx').onclick=closeDetail; el('modal').addEventListener('click',function(e){if(e.target==el('modal'))closeDetail();});
document.addEventListener('keydown',function(e){if(e.key=='Escape'&&st.detail)closeDetail();});

var hadSaved=loadState();
Array.prototype.forEach.call(document.querySelectorAll('#thbox input'),function(inp){var k=inp.getAttribute('data-t');if(k in thr)inp.value=thr[k];});
if(R.n_ancient)Array.prototype.forEach.call(document.querySelectorAll('#athbox input'),function(inp){var k=inp.getAttribute('data-t');if(k in athr)inp.value=athr[k];});
recompute();
if(!(hadSaved&&Object.keys(st.excl).length))R.samples.forEach(function(s){if(s.v=='FAIL')st.excl[s.s]=1;});  // preselect FAILs unless a saved basket exists
renderAll();
})();
"""

SHELL = """<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>__TITLE__</title>
<style>__CSS__</style></head><body>
<header><span class="logo">sBAM<b>piro</b> QC</span><span class="meta" id="meta"></span>
<nav><a href="#gstats">Stats</a><a href="#linsum" id="nav-lin">Lineages</a><a href="#dist">Distributions</a><a href="#corr">Correlations</a><a href="#corrmatrix">Corr matrix</a><a href="#qcpca" id="nav-pca">QC space</a><a href="#divcomp" id="nav-divcomp">Divergence</a><a href="#cons">Consensus</a><a href="#genome" id="nav-genome">Genome</a><a href="#function" id="nav-function">Function</a><a href="#geneburden" id="nav-geneburden">Gene burden</a><a href="#hotspots" id="nav-hot">Variable genes</a><a href="#temporal" id="nav-temporal">Temporal</a><a href="#pnps" id="nav-pnps">pN/pS</a><a href="#adna" id="nav-adna">aDNA</a><a href="#flagged">Flagged</a></nav></header>
<div class="wrap">
<section class="hero"><div class="summary" id="summary"></div><div class="chips" id="chips"></div></section>
<div class="provbar"><div class="prov" id="prov"></div><button class="btn" id="printBtn" title="expand + print / save as PDF">⎙ print</button></div>
<section><details class="dd" style="display:inline-block"><summary>⚙ Live thresholds &amp; presets - adjust and everything re-flags</summary>
<div class="menu" style="display:flex;gap:12px;flex-wrap:wrap;align-items:flex-end;min-width:520px">
  <div style="display:flex;gap:12px;flex-wrap:wrap;align-items:flex-end;flex-basis:100%" id="thbox"></div>
  <div style="display:flex;gap:12px;flex-wrap:wrap;align-items:flex-end;flex-basis:100%" id="athbox"></div>
</div></details></section>
<section id="linsum"><h2>Per-lineage summary <span class="c">- medians per lineage; # MIXED = samples with &gt;1 lineage above the mixture cut-off; click a row to filter</span></h2>
<div class="panel"><div id="lincomp"></div><div class="gtable"><table id="linsumtable"></table></div></div></section>
<section id="gstats"><h2>General statistics <span class="c">- tick a box to basket a sample for exclusion; click a sample name for its full profile; a header to sort</span><button class="exp-h" data-panel="gstatsPanel" data-render="table" style="margin-left:auto">⤢ full</button></h2>
<div class="controls">
  <input id="q" type="search" placeholder="filter samples…">
  <label><input id="of" type="checkbox"> only flagged</label>
  <label id="groupui"><input id="glin" type="checkbox"> group by lineage</label>
  <span id="ancfilter"></span>
  <details class="dd"><summary>columns ▾</summary><div class="menu" id="colmenu"></div></details>
  <details class="dd"><summary>? metric help</summary><div class="menu" id="helpmenu"></div></details>
  <button class="btn" id="csv">↓ export TSV</button>
  <span class="hint"><span id="nshown"></span></span>
</div>
<div class="panel gtable" id="gstatsPanel"><table id="gstable"></table></div>
<div class="curation" id="curation"></div></section>
<section id="dist"><h2>Distributions <span class="c">- one mark per sample; shaded band = acceptable range (modern gate); hover for detail</span>
<span class="seg" id="colorby" style="margin-left:auto"><button class="on" data-cb="qc">colour: QC</button><button data-cb="lineage">lineage</button></span>
<span class="seg" id="ptype"><button class="on" data-t="beeswarm">beeswarm</button><button data-t="bar">bar</button><button data-t="histogram">histogram</button></span>
<button class="exp-h" data-panel="plotsPanel" data-render="plots">⤢ full</button></h2>
<div class="panel" id="plotsPanel"><div id="plots" class="pad" style="padding-top:6px;padding-bottom:6px"></div>
<div class="legend" id="leg-qc"><span><i style="background:#94a3b8"></i>PASS</span><span><i style="background:#d97706"></i>WARN</span><span><i style="background:#dc2626"></i>FAIL</span><span style="margin-left:auto">beeswarm dashed line = median</span></div>
<div class="legend" id="leg-lin" style="display:none"></div></div></section>
<section id="corr"><h2>Correlations <span class="c">- pick two metrics, coloured by QC; drag a box to basket the enclosed samples</span>
<span id="scbrushinfo" style="margin-left:auto;font-size:11px;color:var(--accent)"></span>
<select class="msel" id="sx"></select><span style="color:#94a3b8">vs</span><select class="msel" id="sy"></select>
<button class="exp-h" data-panel="scatter" data-render="scatter">⤢ full</button></h2>
<div class="panel pad" id="scatter"></div></section>
<section id="corrmatrix"><h2>Metric correlation <span class="c">- Spearman rank correlation across the core metrics, over the samples in view; click a cell to load that pair into the scatter above</span>
<button class="exp-h" data-panel="corrmatrix" data-render="corr" style="margin-left:auto">⤢ full</button></h2>
<div class="panel pad" id="corr_body"></div></section>
<section id="qcpca"><h2>QC-metric space <span class="c">- samples ordinated in standardized QC-metric space (PCA on robust z-scores); a descriptive map of how quality profiles co-vary, not phylogeny and not a test</span>
<span class="seg" id="pcacb" style="margin-left:auto"><button class="on" data-cb="qc">colour: QC</button><button data-cb="lineage">lineage</button></span>
<button class="exp-h" data-panel="qcpcaPanel" data-render="pca">⤢ full</button></h2>
<div class="panel" id="qcpcaPanel"><div id="qcpca_body" class="pad"></div>
<div class="hot-note" id="qcpca_caption"></div>
<div class="dsub" style="margin:14px 16px 4px">Most unusual samples <span style="font-weight:400;color:#94a3b8">(robust Mahalanobis rank + top deviating metrics)</span></div>
<div class="gtable" style="max-height:34vh"><table id="pca_outtable"></table></div></div></section>
<section id="divcomp"><h2>Divergence vs completeness <span class="c">- separates reference-bias suspects (low divergence at low missingness) from honest low-coverage (low divergence explained by high missingness); a screen to inspect, not an automatic flag</span>
<span class="seg" id="divx" style="margin-left:auto"><button class="on" data-x="missing_pct">x: missing %</button><button data-x="callable_inv">100 - callable %</button></span>
<button class="exp-h" data-panel="divcompPanel" data-render="divcomp">⤢ full</button></h2>
<div class="panel" id="divcompPanel"><div id="divcomp_body" class="pad"></div>
<div class="legend" id="divcomp_quadn"></div>
<div class="hot-note" id="divcomp_caption"></div></div></section>
<section id="cons"><h2>Consensus completeness <span class="c">- callable / IUPAC / missing per sample, worst first</span></h2>
<div class="panel"><div id="stacks" style="max-height:60vh;overflow:auto"></div>
<div class="legend"><span><i style="background:#22a06b"></i>callable</span><span><i style="background:#e6b25a"></i>IUPAC</span><span><i style="background:#cbd5e1"></i>missing (- / N)</span></div></div></section>
<section id="genome"><h2>Genome landscape <span class="c">- a signal along the reference; rows = samples (grouped by lineage). A contiguous block = a localised feature (e.g. an RD deletion, a variant cluster). Missing = consensus callability (coverage proxy); SNPs/Het/Indels = VCF variant density</span>
<input class="gsearch" id="genegoto" type="search" placeholder="go to gene" style="margin-left:auto"><span id="gselreadout" style="font-size:11px;color:var(--accent)"></span><button class="btn" id="maskbtn">⛆ mask regions</button><span class="seg" id="gtrack"></span>
<button class="exp-h" data-panel="genomePanel" data-render="genome">⤢ full</button></h2>
<div class="panel" id="genomePanel"><div class="genome-scroll"><div id="genome_body" class="pad"></div></div>
<div class="legend"><span><i style="background:rgb(238,244,240)"></i>none / callable</span><span><i style="background:rgb(160,120,150)"></i>&#8594;</span><span><i style="background:rgb(120,70,90)"></i>high</span><span style="margin-left:auto">top strip = cohort mean; hover a cell for the position + value</span></div></div></section>
<section id="function"><h2>Functional annotation <span class="c">- snpEff impact + effect classes per sample; the missense/silent ratio is a pN/pS PROXY, not a selection test</span>
<button class="exp-h" data-panel="functionPanel" data-render="function" style="margin-left:auto">⤢ full</button></h2>
<div class="panel" id="functionPanel">
  <div class="provbar" style="margin:12px 16px 4px"><div class="prov" id="fn_cohort"></div></div>
  <div id="fn_stacks" style="max-height:56vh;overflow:auto"></div>
  <div class="legend" id="fn_legend"></div>
  <div class="hot-note" id="fn_caption"></div>
</div></section>
<section id="geneburden"><h2>Functional gene burden <span class="c">- genes carrying the most HIGH+MODERATE variants across the cohort, with the dominant effect class; a functional companion to Variable genes below, not a selection test</span>
<input class="gsearch" id="gbq" type="search" placeholder="search gene" style="margin-left:auto"><button class="exp-h" data-panel="geneburdenPanel" data-render="geneburden">⤢ full</button></h2>
<div class="panel" id="geneburdenPanel"><div id="gb_body"><div class="hot-note" id="gb_note"></div><div class="gtable" style="max-height:44vh"><table id="gbtable"></table></div></div></div></section>
<section id="hotspots"><h2>Variable genes <span class="c">- genes/regions with the most SNPs across the cohort; click a row to mark it on the SNP track</span>
<input class="gsearch" id="hotq" type="search" placeholder="search gene" style="margin-left:auto"><button class="exp-h" data-panel="hotspotsPanel" data-render="hotspots">⤢ full</button></h2>
<div class="panel" id="hotspotsPanel"><div id="hot_body"><div class="hot-note" id="hot_note"></div><div class="gtable" style="max-height:44vh"><table id="hottable"></table></div></div></div></section>
<section id="dynamics"><h2>SNP dynamics <span class="c">- allele-frequency trajectories over time per connected series; points flag emergence / fixation / loss / non-synonymous. Appears only when the metadata carries time + group columns.</span></h2>
<div class="panel pad" id="dyn_body"></div></section>
<section id="adna"><h2>aDNA damage authentication <span class="c">- terminal deamination per ancient sample; a screen, not a proof of authenticity</span></h2>
<div class="panel"><div id="adna_body" class="pad"></div></div></section>
<section id="temporal"><h2>Temporal sampling overview <span class="c">- per-lineage year span parsed from dates + an informative-site proxy; a readiness check for downstream time-resolved analysis, NOT a clock estimate</span></h2>
<div class="panel"><div id="temporal_body" class="pad"></div>
<div class="legend"><span><i style="background:var(--accent)"></i>sample year</span><span>bar = sampling span</span><span style="margin-left:auto">informative-site proxy = median SNPs</span></div>
<div class="hot-note" id="temporal_caption"></div></div></section>
<section id="pnps"><h2>Selection: pN/pS and dN/dS <span class="c">- alignment-based per-gene dN/dS from eskaks; a cohort selection screen, not a per-sample QC metric</span>
<input class="gsearch" id="pnpsq" type="search" placeholder="search gene" style="margin-left:auto"><button class="exp-h" data-panel="pnpsPanel" data-render="pnps">⤢ full</button></h2>
<div class="panel" id="pnpsPanel"><div id="pnps_body"><div class="hot-note" id="pnps_note"></div><div class="gtable" style="max-height:48vh"><table id="pnpstable"></table></div></div></div></section>
<section id="flagged"><h2>Flagged samples <span class="c">- <span id="nflag"></span> to review; click a flag to filter, hover for the margin</span>
<button class="btn" id="basketFlagged" style="margin-left:auto">↧ basket all flagged</button></h2>
<div class="panel gtable" style="max-height:50vh"><table id="flagtable"></table></div></section>
<div class="footer" id="foot"></div>
</div>
<button id="expClose" class="exp-close" aria-label="exit fullscreen">✕ close (Esc)</button>
<div id="modal" class="modal"><div class="modalcard"><button class="modalx" id="modalx" aria-label="close">×</button><div id="modalbody"></div></div></div>
<div id="tt"></div>
<script>const REPORT=__JSON__;</script>
<script>__JS__</script>
</body></html>"""


def build_html(title, payload):
    data = json.dumps(payload, separators=(",", ":")).replace("</", "<\\/")
    return (SHELL.replace("__TITLE__", html.escape(title)).replace("__CSS__", CSS)
                 .replace("__JS__", JS).replace("__JSON__", data))


# ============================================================================
#  SNP DYNAMICS (optional): flexible metadata + per-sample VCFs -> per-connected-group
#  allele-frequency trajectories over time with flagged events. Self-hides if the metadata
#  has no time+group columns or no VCFs are given.
# ============================================================================
_DYN_SAMPLE_RE = re.compile(r'^(sample_?id|sampleid|sample|name|gid|strain|isolate)$', re.I)
_DYN_TIME_RE   = re.compile(r'(passage|pase|timepoint|time_?point|^time$|^day$|date|week|month|hour|generation|^tp$|visit|^t\d*$)', re.I)
_DYN_GROUP_RE  = re.compile(r'(group|series|patient|host|subject|cluster|experiment|^line$|replicate|chain|pair|lineage_?id|donor|case|animal)', re.I)
_DYN_NONSYN    = re.compile(r'missense|stop_gained|stop_lost|start_lost|frameshift|inframe|splice|initiator', re.I)


def _dyn_open(path):
    return gzip.open(path, 'rt', encoding='utf-8', errors='replace') if str(path).endswith('.gz') \
        else open(path, 'r', encoding='utf-8', errors='replace')


def _dyn_num(x):
    if x is None:
        return None
    m = re.search(r'-?\d+\.?\d*', str(x))
    return float(m.group()) if m else None


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


def _dyn_ann(info):
    for field in info.split(';'):
        if field.startswith('ANN='):
            a = field[4:].split(',')[0].split('|')
            return (a[3] if len(a) > 3 else ''), (a[1] if len(a) > 1 else ''), (a[2] if len(a) > 2 else '')
    return '', '', ''


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
                    gene, eff, imp = _dyn_ann(c[7])
                    out[sample][f'{chrom}:{pos}'] = {'ref': ref, 'alt': alt.split(',')[0], 'gene': gene,
                                                     'eff': eff, 'imp': imp, 'af': round(af, 4)}
        except OSError:
            continue
    return out


def build_dynamics(metadata, variants, emerge=0.25, fix=0.90, loss=0.10, min_points=2, min_move=0.15):
    """Per-connected-group AF trajectories over time, only for SNPs that move. None if nothing to show."""
    if not metadata:
        return None
    groups = {}
    for s, md in metadata.items():
        g = md.get('group')
        if g and s in variants:
            groups.setdefault(g, []).append(s)
    out_groups = []
    for g, samples in sorted(groups.items()):
        samples = sorted(samples, key=lambda s: (metadata[s]['tnum'] is None,
                                                 metadata[s]['tnum'] if metadata[s]['tnum'] is not None else 0,
                                                 str(metadata[s]['time'])))
        times = [metadata[s]['time'] for s in samples]
        if len(set(times)) < min_points:
            continue
        allpos = set()
        for s in samples:
            allpos.update(variants[s].keys())
        series = []
        for pos in allpos:
            traj, meta = [], None
            for s in samples:
                v = variants[s].get(pos)
                traj.append(v['af'] if v else 0.0)
                if v and meta is None:
                    meta = v
            if max(traj) - min(traj) < min_move:
                continue
            flags = []
            if traj[0] <= loss and max(traj) >= emerge:
                flags.append('emergence')
            if traj[-1] >= fix and traj[0] < fix:
                flags.append('fixation')
            if traj[0] >= emerge and traj[-1] <= loss:
                flags.append('loss')
            if meta and _DYN_NONSYN.search(meta.get('eff', '')):
                flags.append('nonsyn')
            if meta and meta.get('imp', '') == 'HIGH':
                flags.append('high_impact')
            series.append({'pos': pos, 'gene': (meta or {}).get('gene', ''), 'eff': (meta or {}).get('eff', ''),
                           'imp': (meta or {}).get('imp', ''), 'alt': (meta or {}).get('alt', ''),
                           'traj': [round(x, 4) for x in traj], 'flags': flags})
        if not series:
            continue
        series.sort(key=lambda x: (-(max(x['traj']) - min(x['traj'])), -len(x['flags'])))
        out_groups.append({'group': g, 'samples': samples, 'times': times,
                           'tnums': [metadata[s]['tnum'] for s in samples],
                           'series': series, 'n_flagged': sum(1 for x in series if x['flags'])})
    if not out_groups:
        return None
    return {'groups': out_groups, 'thresholds': {'emerge': emerge, 'fix': fix, 'loss': loss}}


def main():
    ap = argparse.ArgumentParser(description="sBAMpiro interactive QC report + flags.")
    ap.add_argument("--summary", required=True)
    ap.add_argument("--consensus", nargs="*", default=[])
    ap.add_argument("--out-html", required=True)
    ap.add_argument("--out-flags", required=True)
    ap.add_argument("--title", default="sBAMpiro QC report")
    for k, v in DEF.items():
        ap.add_argument("--" + k.replace("_", "-"), type=float, default=v)
    for k, v in ANC_DEF.items():
        ap.add_argument("--anc-" + k.replace("_", "-"), type=float, default=v)
    ap.add_argument("--mapdamage-dir", nargs="*", default=[],
                    help="mapDamage2 output dir(s) for ancient samples (per-sample subdir or a parent to walk).")
    ap.add_argument("--gff", default=None,
                    help="GFF3 of gene coordinates -> per-gene SNP-density hotspot detection (optional).")
    ap.add_argument("--lineage-colors", default=None,
                    help="TSV 'lineage<TAB>#hex' (e.g. mycolorsTB) to colour lineages with a canonical palette (optional).")
    ap.add_argument("--mask-bed", default=None,
                    help="BED of regions to mask (e.g. mtbc_mask.bed: PE/PPE, IS, DR, repeats) -> a 'mask regions' toggle.")
    ap.add_argument("--gene-burden", default=None,
                    help="Cohort gene-burden TSV (collect_summary --gene-burden-out) -> the Functional gene burden panel (optional).")
    ap.add_argument("--pnps", default=None,
                    help="Cohort per-gene dN/dS TSV (eskaks) -> the Selection pN/pS panel (optional; a cohort analysis, not per-sample QC).")
    ap.add_argument("--provenance", nargs="*", default=[],
                    help="key=value provenance pairs surfaced in the report header (reference, container, commit...).")
    ap.add_argument("--metadata", default=None,
                    help="Optional TSV (e.g. the samplesheet) with a sample column + time (passage/timepoint) "
                         "and group (patient/series/cluster) columns -> the SNP dynamics panel. Auto-detected.")
    ap.add_argument("--vcfs", nargs="*", default=[],
                    help="Optional per-sample annotated VCFs -> per-SNP allele frequencies for the dynamics panel.")
    ap.add_argument("--gate", action="store_true")
    args = ap.parse_args()
    thr = {k: getattr(args, k) for k in DEF}
    anc_thr = {k: getattr(args, "anc_" + k) for k in ANC_DEF}

    summ = parse_summary(args.summary)
    miss_by_sample = {}    # per-sample binned missingness profile along the reference (genome landscape)
    for path in args.consensus:
        try:
            cs = consensus_stats(path)
        except Exception as e:
            sys.stderr.write(f"[qc_report] WARN consensus {path}: {e}\n")
            continue
        sid = cs.get("consensus_id")
        if sid not in summ:
            stem = os.path.basename(path).split(".")[0]
            sid = stem if stem in summ else sid
        keys = ("length", "missing_pct", "iupac_pct", "callable_pct", "longest_n_run", "n_gaps")
        if sid in summ:
            summ[sid].update({k: cs[k] for k in keys})
        elif sid:
            summ[sid] = {"sample_id": sid, **{k: cs[k] for k in keys}}
        if sid:
            miss_by_sample[sid] = cs.get("miss_profile")

    # SNP-count robust-z is computed over MODERN samples only (ancient counts are low by design and would skew it)
    snp_med, snp_sig = robust([to_float(m.get("snps")) for m in summ.values() if not is_ancient(m)])
    # representative genome length: consensus length, else the reference length column (works when make_consensus is off)
    genome_len = max([int(v) for v in (to_float(m.get("length")) or to_float(m.get("genome_len"))
                                       for m in summ.values()) if v], default=0)
    snp_density_ok = genome_len > 100000

    # ---- aDNA damage authentication: map each mapDamage dir to a sample id ----
    ancient_ids = {sid for sid, m in summ.items() if is_ancient(m)}
    dmg_by_sample = {}

    def _index_dir(d):
        if not os.path.isdir(d):
            return
        base = os.path.basename(os.path.normpath(d))
        if base in ancient_ids and base not in dmg_by_sample:
            st = mapdamage_stats(d)
            if st:
                dmg_by_sample[base] = st
        for aid in ancient_ids:
            sub = os.path.join(d, aid)
            if aid not in dmg_by_sample and os.path.isdir(sub):
                st = mapdamage_stats(sub)
                if st:
                    dmg_by_sample[aid] = st
        if any(os.path.exists(os.path.join(d, f)) for f in ("5pCtoT_freq.txt", "3pGtoA_freq.txt")):
            comps = os.path.normpath(d).split(os.sep)   # match a whole path COMPONENT, not a loose substring
            for aid in ancient_ids:
                if aid not in dmg_by_sample and aid in comps:
                    st = mapdamage_stats(d)
                    if st:
                        dmg_by_sample[aid] = st

    for d in (args.mapdamage_dir or []):
        _index_dir(d)
        if os.path.isdir(d):
            for name in sorted(os.listdir(d)):
                p = os.path.join(d, name)
                _index_dir(p)
                sub = os.path.join(p, "mapDamage")   # published layout: {OUTDIR}/{sample_dir}/mapDamage
                if os.path.isdir(sub):
                    _index_dir(sub)
    if ancient_ids and not dmg_by_sample:
        sys.stderr.write(f"[qc_report] {len(ancient_ids)} ancient sample(s) but no mapDamage output matched; "
                         f"DAMAGE authentication skipped.\n")

    metric_keys = [k for k, _, _, _ in METRICS]
    # fine-grained snpEff effect-class columns: carried into s.m for the Functional-annotation panel, not registered metrics
    EFF_KEYS = ["eff_missense", "eff_synonymous", "eff_stop_gained", "eff_stop_lost", "eff_start_lost",
                "eff_frameshift", "eff_inframe_indel", "eff_splice", "eff_intergenic", "eff_regulatory"]
    extra_metrics = discover_extra_metrics(summ)
    extra_keys = [e["key"] for e in extra_metrics]
    jsamples, counts = [], {"PASS": 0, "WARN": 0, "FAIL": 0}
    for sid, m in summ.items():
        anc = is_ancient(m)
        dmg = dmg_by_sample.get(sid)
        verdict, flags = flag_sample(m, thr, snp_med, snp_sig, ancient=anc, anc_thr=anc_thr, dmg=dmg)
        counts[verdict] += 1
        jsamples.append({"s": sid, "v": verdict, "f": flags,
                         "lineage": clean_str(m.get("lineage")),
                         "hetf": het_frac(m),
                         "linf": lineage_fracs(m.get("lineage_counts")) or None,
                         "linc": lineage_counts_parsed(m.get("lineage_counts")) or None,
                         "dr": clean_str(m.get("drug_resistance")),
                         "anc": anc,
                         "dmg": dmg,
                         "date": clean_str(m.get("date")),
                         "miss": miss_by_sample.get(sid),
                         "trk": ({k: v for k, v in (("snp", parse_profile(m.get("snp_profile"))),
                                                    ("het", parse_profile(m.get("het_profile"))),
                                                    ("indel", parse_profile(m.get("indel_profile")))) if v is not None}
                                 or None),
                         "ann_db_error": (clean_str(m.get("ann_db_error")) == "yes"),
                         "m": {k: to_float(m.get(k)) for k in metric_keys + extra_keys + EFF_KEYS}})

    lineages = sorted({s["lineage"] for s in jsamples if s["lineage"]})
    n_ancient = sum(1 for s in jsamples if s["anc"])
    mask_iv = parse_bed(args.mask_bed)
    mask_bins, mask_pct = mask_profile(mask_iv, genome_len)
    provenance = {}
    for kv in (args.provenance or []):
        if "=" in kv:
            k, v = kv.split("=", 1)
            provenance[k.strip()] = v.strip()

    now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    payload = {"generated": now, "counts": counts, "thresholds": thr, "dist": DIST,
               "genome_len": genome_len, "snp_density_ok": snp_density_ok, "defs": DEFS, "nbins": NBINS,
               "genes": parse_gff(args.gff), "lin_colors": parse_lineage_colors(args.lineage_colors),
               "gene_burden": parse_gene_burden(args.gene_burden) or None,
               "pnps": parse_pnps(args.pnps) or None,
               "mask_bins": mask_bins, "mask_pct": mask_pct,
               "mask_iv": (mask_iv[:5000] if mask_iv else None),
               "lineages": lineages, "lin_present": len(lineages) > 0,
               "n_ancient": n_ancient, "anc_thresholds": anc_thr, "provenance": provenance,
               "dynamics": build_dynamics(parse_metadata(args.metadata), parse_vcfs(args.vcfs)),
               "metrics": [{"key": k, "label": l, "kind": kind, "dir": d} for k, l, kind, d in METRICS],
               "extra": extra_metrics,
               "samples": jsamples}

    os.makedirs(os.path.dirname(os.path.abspath(args.out_html)) or ".", exist_ok=True)
    with open(args.out_html, "w", encoding="utf-8") as fh:
        fh.write(build_html(args.title, payload))
    with open(args.out_flags, "w", encoding="utf-8") as fh:
        fh.write("sample\tverdict\tmean_depth\tbreadth_pct\tmapped_pct\tmissing_pct\tiupac_pct\tsnps\tti_tv\tflags\n")
        for s in sorted(jsamples, key=lambda x: (x["v"] != "FAIL", x["v"] != "WARN", x["s"])):
            m = s["m"]
            def g(k):
                v = m.get(k)
                return "" if v is None else f"{v:.4g}"
            fh.write("\t".join([s["s"], s["v"], g("mean_depth"), g("breadth_pct"), g("mapped_pct"),
                                g("missing_pct"), g("iupac_pct"), g("snps"), g("ti_tv"),
                                ",".join(s["f"]) or "."]) + "\n")

    sys.stderr.write(f"[qc_report] {len(jsamples)} samples -> {counts['PASS']} PASS, {counts['WARN']} WARN, "
                     f"{counts['FAIL']} FAIL. Wrote {args.out_html} + {args.out_flags}\n")
    for s in sorted([x for x in jsamples if x["v"] == "FAIL"], key=lambda x: x["s"]):
        sys.stderr.write(f"    FAIL  {s['s']:<40s} [{','.join(s['f'])}]\n")
    if args.gate and counts["FAIL"]:
        sys.stderr.write(f"[qc_report] GATE: {counts['FAIL']} sample(s) FAIL -> blocking.\n")
        sys.exit(1)


if __name__ == "__main__":
    main()
