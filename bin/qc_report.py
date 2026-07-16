#!/usr/bin/env python3
"""Consolidated interactive QC report for BAMpiro (organism-agnostic).

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
import random
import re
import sys
import zlib
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
REPO_URL = "https://github.com/PathoGenOmics-Lab/BAMpiro"   # surfaced in the report header + footer

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


_DR_DRUG_ORDER = ["RIF", "INH", "EMB", "PZA", "STR", "STM", "FQ", "LFX", "MFX", "OFX", "KAN", "AMK",
                  "CAP", "ETH", "PTO", "LZD", "BDQ", "CFZ", "BDQ_CFZ", "DLM", "PAS", "CS"]


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
                try:
                    dp = int(dp) if dp not in (None, "", "NA", ".") else None
                except ValueError:
                    dp = None
                calls.append({"s": s, "drug": (d.get("drug") or "").strip(), "gene": (d.get("gene") or "").strip(),
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
    drugs = sorted({c["drug"] for c in calls if c["drug"]},
                   key=lambda x: (_DR_DRUG_ORDER.index(x) if x in _DR_DRUG_ORDER else 99, x))
    return {"samples": samples, "drugs": drugs, "calls": calls}


def parse_kraken(paths):
    """Per-sample Kraken2 `.report` files (6 columns: clade%, clade reads, taxon reads, rank, taxid,
    name) -> {'samples':[{s, unclassified, classified, primary, secondary, top}]} or None. The sample id
    is the filename up to the first '__' (the pipeline writes '<sampleId>__<runId>.kraken.report'). The
    'primary' taxon is the top species-level clade; contamination shows as a low primary% / large
    secondary taxon / high unclassified%."""
    out = []
    for p in paths or []:
        if not p or not os.path.exists(p):
            continue
        sid = os.path.basename(p).split("__")[0].split(".")[0]
        unclass, species = 0.0, []
        try:
            with open(p, encoding="utf-8", errors="replace") as fh:
                for line in fh:
                    c = line.rstrip("\n").split("\t")
                    if len(c) < 6:
                        continue
                    try:
                        pct = float(c[0])
                    except ValueError:
                        continue
                    rank, name = c[3].strip(), c[5].strip()
                    if rank == "U":
                        unclass = pct
                    elif rank == "S":            # species-level clades = the interpretable composition
                        species.append((pct, name))
        except OSError:
            continue
        if not species and unclass == 0.0:
            continue
        species.sort(key=lambda x: -x[0])
        top = [{"name": n, "pct": round(p, 2)} for p, n in species[:6]]
        out.append({"s": sid, "unclassified": round(unclass, 2),
                    "classified": round(100.0 - unclass, 2),
                    "primary": (top[0] if top else None),
                    "secondary": (top[1] if len(top) > 1 else None), "top": top})
    if not out:
        return None
    return {"samples": out}


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


# Curated H37Rv locus tags (Mycobrowser) for the common MTB resistance genes, used as a fallback when
# the run's reference GFF carries no H37Rv-style locus_tag for them (e.g. a non-H37Rv reference). These
# are the WHO-catalogue first/second-line resistance genes.
GENE_RV = {
    "rpoB": "Rv0667", "rpoC": "Rv0668", "katG": "Rv1908c", "inhA": "Rv1484", "fabG1": "Rv1483",
    "ahpC": "Rv2428", "gyrA": "Rv0006", "gyrB": "Rv0005", "embA": "Rv3794", "embB": "Rv3795",
    "embC": "Rv3793", "pncA": "Rv2043c", "rpsL": "Rv0682", "eis": "Rv2416c", "ethA": "Rv3854c",
    "ethR": "Rv3855", "gid": "Rv3919c", "gidB": "Rv3919c", "tlyA": "Rv1694", "Rv0678": "Rv0678",
    "atpE": "Rv1305", "pepQ": "Rv2535c", "folC": "Rv2447c", "thyA": "Rv2764c", "alr": "Rv3423c",
    "ddn": "Rv3547", "fgd1": "Rv0407", "fbiA": "Rv3261", "fbiB": "Rv3262", "fbiC": "Rv1173",
    "rrs": "MTB000019", "rrl": "MTB000020", "mmpL5": "Rv0676c", "mmpS5": "Rv0677c", "whiB7": "Rv3197A",
    "clpC1": "Rv3596c", "panD": "Rv3601c", "ddlA": "Rv2981c",
}
_RV_LOCUS_RE = re.compile(r'^(Rv\d|MTB\d)', re.I)   # H37Rv (Mycobrowser) locus-tag schemes


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

                def attr(key):
                    i = attrs.find(key)
                    return attrs[i + len(key):].split(";")[0].strip() if i >= 0 else ""
                locus = attr("locus_tag=")
                name = attr("gene=") or attr("Name=")
                if name and locus and name not in out:
                    out[name] = locus
    except (OSError, ValueError):
        return out
    return out


def build_gene_map(gff_path):
    """{gene: Mycobrowser (H37Rv) locus tag}. Curated resistance-gene map, overlaid with any H37Rv-style
    locus tags from the run's GFF (so an H37Rv reference contributes every gene; a non-H37Rv reference's
    strain-specific tags are ignored in favour of the curated H37Rv equivalents)."""
    m = dict(GENE_RV)
    for gene, locus in parse_gene_locus(gff_path).items():
        if _RV_LOCUS_RE.match(locus):
            m[gene] = locus
    return m


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
 --accent:#0b7f97;--accent-soft:#e3f3f7;   /* 4.67:1 on white -> WCAG AA for the action buttons / active chips (was #0e8ba8 at 3.98:1) */
 --pass:#0f9d6b;--warn:#dd8a1a;--fail:#e23a4a;--good:#16a37a;--bad:#e5615c;--neu:#4a90b8;
 --bandfill:rgba(15,157,107,.09);--bandedge:rgba(15,157,107,.42);
 --sh:0 1px 2px rgba(16,24,40,.04),0 3px 8px rgba(16,24,40,.05);--r:15px;
 --label:#33465c;--txt2:#516074;--track:#eef2f7}
/* ---- Dark theme (html.dark, set by the early head script from the saved pref or the OS) ---- */
html.dark{color-scheme:dark;
 --bg:#0f1720;--panel:#18232f;--line:#2a3543;--ink:#e7edf4;--mut:#94a4b6;--soft:#1e2a38;
 --label:#c1cede;--txt2:#a3b2c4;--track:#243141;
 --accent:#3bb8d8;--accent-soft:#123640;
 --pass:#2fbb8b;--warn:#eaa63c;--fail:#f26274;--good:#34c294;--bad:#f07b87;--neu:#5ba8d9;
 --bandfill:rgba(47,187,139,.13);--bandedge:rgba(47,187,139,.46);
 --sh:0 1px 2px rgba(0,0,0,.3),0 4px 14px rgba(0,0,0,.36)}
html.dark body{background-color:#0f1720;background-image:radial-gradient(1200px 520px at 50% -260px,#18242f,transparent 70%),linear-gradient(180deg,#131d27 0%,#0d151e 100%)}
html.dark header{background:rgba(20,30,40,.82);box-shadow:0 1px 0 rgba(0,0,0,.25),0 8px 24px rgba(0,0,0,.32)}
html.dark #toc{background:rgba(19,28,38,.97)}
html.dark #toc-toggle{background:var(--panel);color:var(--mut)}
html.dark .toc-link{color:#9fb0c2}
html.dark h2 .c{color:#8798ab}
html.dark .controls input[type=search],html.dark .menu input,html.dark #thbox input,html.dark #athbox input,html.dark tr.colfilt .cfx{background:var(--soft);color:var(--ink);border-color:var(--line)}
html.dark th{background:#1b2836;color:#9fb0c2;box-shadow:0 1px 0 var(--line)}
html.dark th:hover{background:#213142} html.dark th.s{background:#1b2836}
html.dark tr.colfilt th{background:#172431}
html.dark td.na{color:#5a6b7d}
html.dark tbody tr:hover td{background:#1e2b3a}
html.dark tr.hl td{background:#33371c!important} html.dark tr.hl td.s{background:#33371c!important}
html.dark .v.PASS{color:#7fe4be;background:#123727} html.dark .v.WARN{color:#f1c179;background:#3a2c10} html.dark .v.FAIL{color:#f5a1ac;background:#3b1a22}
html.dark .chip{color:#aebccb} html.dark .chip.on{color:#0f1720}
html.dark .modal{background:rgba(4,8,12,.62)}
html.dark #toc-toggle:hover{color:var(--accent);border-color:var(--accent)}
html.dark #themeToggle,html.dark #ghlink{color:var(--mut)}
html.dark .btn,html.dark details.dd>summary,html.dark .seg button,html.dark .exp-h{background:var(--panel);color:#b7c4d3}
html.dark .btn:hover,html.dark .exp-h:hover,html.dark details.dd>summary:hover{border-color:var(--accent);color:var(--accent)}
html.dark details.dd .menu{background:#1b2836;box-shadow:0 16px 44px rgba(0,0,0,.55)}
html.dark .seg button.on{background:var(--accent);color:#0f1720}
html.dark .btn.prim{background:var(--accent);color:#0f1720;border-color:var(--accent)}
html.dark .acard{background:linear-gradient(180deg,var(--soft),var(--panel))}
html.dark .curation{background:linear-gradient(180deg,var(--soft),var(--panel))}
html.dark #nbasket{color:var(--accent)}
html.dark code,html.dark .kbd{background:#243141;color:#cdd8e4}
html.dark .abadge{color:#e8c37a;background:#3a2c10;border-color:#5a4a24}
html.dark .chip.failc{background:#3b1a22;color:#f5a1ac;border-color:#5a2b33}
html.dark tr.lingrp td{background:#1b2836}
html.dark .acard.low{border-color:#6a3540;background:linear-gradient(180deg,var(--soft),#2c1a20)}
html.dark .epi-moderate{background:#1c3350;color:#8fb3e0} html.dark .epi-weak{background:#232f3d;color:#93a0b0}
html.dark .modalx:hover{background:#2a3a4c}
html.dark .hitem{border-color:var(--line)}
html.dark .infoi{background:#2a3a4c;color:#9fb0c2} html.dark th .infoi,html.dark .dyn-legend .infoi{background:#2a3a4c}
html.dark .av.good{color:#4fd0a0} html.dark .av.bad{color:#f0808c} html.dark .av.na{color:#5a6b7d}
html.dark .alow,html.dark .dwhy b{color:#f0808c}
html.dark .epimx-cell,html.dark .drmx-cell{border-color:#0f1720}
html.dark table.snpmx th,html.dark table.snpmx td{border-bottom-color:#0f1720}
html.dark table.snpmx .snpmx-metacell{border-bottom-color:#0f1720}
html.dark .epimx-diag{background:repeating-linear-gradient(45deg,#2a3543,#2a3543 3px,#222e3c 3px,#222e3c 6px)}
html.dark .epimx-grad{background:linear-gradient(90deg,#a24a8f,#2a3543,#2f8f5b)}
html.dark table.snpmx td.snpmx-empty{background:repeating-linear-gradient(45deg,#1e2a38,#1e2a38 3px,#243141 3px,#243141 6px)}
html.dark .gtable::-webkit-scrollbar-thumb{background:#3a485a;border-color:#18232f}
#themeToggle,#ghlink{background:none;border:1px solid var(--line);border-radius:9px;width:32px;height:32px;cursor:pointer;color:var(--mut);font-size:15px;display:flex;align-items:center;justify-content:center;line-height:1;box-shadow:var(--sh);text-decoration:none;flex:none}
#themeToggle:hover,#ghlink:hover{color:var(--accent);border-color:var(--accent)}
.ver{font-size:11px;font-weight:700;color:var(--accent);background:var(--accent-soft);border-radius:20px;padding:2px 9px;letter-spacing:.2px}
html.dark .ver{background:var(--accent-soft);color:var(--accent)}
*{box-sizing:border-box} html{scroll-behavior:smooth}
body{margin:0;font-family:ui-sans-serif,system-ui,-apple-system,"Segoe UI",Roboto,Helvetica,sans-serif;color:var(--ink);font-size:15px;line-height:1.55;-webkit-font-smoothing:antialiased;
 --tocw:232px;padding-left:var(--tocw);transition:padding-left .2s ease;
 background-color:#eef1f6;
 background-image:radial-gradient(1200px 520px at 50% -260px,#ffffff,transparent 70%),linear-gradient(180deg,#f4f6fa 0%,#e9edf3 100%);
 background-attachment:fixed,fixed;background-repeat:no-repeat,no-repeat}
body.toc-collapsed{padding-left:0}
.num,td{font-variant-numeric:tabular-nums} :focus-visible{outline:2px solid var(--accent);outline-offset:1px}
header{position:sticky;top:0;z-index:40;background:rgba(255,255,255,.82);-webkit-backdrop-filter:saturate(1.15) blur(10px);backdrop-filter:saturate(1.15) blur(10px);color:var(--ink);padding:13px 22px;display:flex;align-items:center;gap:14px;border-bottom:1px solid var(--line);box-shadow:0 1px 0 rgba(16,24,40,.02),0 8px 24px rgba(16,24,40,.035)}
header .logo{font-weight:700;letter-spacing:-.2px;font-size:16px;display:flex;align-items:center;gap:9px;color:var(--ink)}
header .logo b{color:var(--accent);font-weight:700}
header .logo .dot{width:9px;height:9px;border-radius:50%;background:var(--accent);box-shadow:0 0 0 3px rgba(14,139,168,.16)}
header .logo .brandlogo{height:30px;width:auto;flex:0 0 auto}
.toc-brand .brandlogo{height:24px;width:auto;vertical-align:middle;margin-right:7px}
header .meta{color:var(--mut);font-size:12px}
/* Left table-of-contents sidebar (collapsible, grouped, scroll-spy) */
#toc-toggle{position:fixed;top:11px;left:11px;z-index:70;width:34px;height:34px;border:1px solid var(--line);border-radius:9px;background:var(--panel);color:var(--label);font-size:16px;cursor:pointer;box-shadow:var(--sh);display:flex;align-items:center;justify-content:center;line-height:1}
#toc-toggle:hover{color:var(--accent);border-color:var(--accent)}
#tocscrim{display:none;position:fixed;inset:0;z-index:55;background:rgba(8,12,18,.44);-webkit-backdrop-filter:blur(1px);backdrop-filter:blur(1px)}
#toc{position:fixed;left:0;top:0;bottom:0;width:var(--tocw);overflow-y:auto;background:rgba(255,255,255,.97);-webkit-backdrop-filter:blur(8px);backdrop-filter:blur(8px);border-right:1px solid var(--line);z-index:60;padding:54px 12px 26px;transition:transform .2s ease}
body.toc-collapsed #toc{transform:translateX(-100%)}
body.toc-collapsed header{padding-left:56px}
.toc-brand{font-weight:700;font-size:14px;color:var(--ink);padding:0 10px 12px;letter-spacing:-.2px}
.toc-brand b{color:var(--accent)}
.toc-group{margin-bottom:4px}
.toc-gh{display:flex;align-items:center;justify-content:space-between;font-size:11px;font-weight:700;letter-spacing:.5px;text-transform:uppercase;color:var(--mut);padding:9px 10px 4px;cursor:pointer;user-select:none;border-radius:7px}
.toc-gh:hover{color:var(--ink)}
.toc-chev{transition:transform .15s;opacity:.7;display:inline-flex}
.toc-group.closed .toc-chev{transform:rotate(-90deg)}
/* homogeneous icon family (inline Lucide SVG, inherits text colour) */
.ic{width:1.05em;height:1.05em;flex:none;vertical-align:-.16em;stroke-width:2;margin-right:.36em}
.ic.sort{width:.85em;height:.85em;margin:0 0 0 .25em;opacity:.85;vertical-align:-.05em}
.toc-chev .ic{width:13px;height:13px;margin:0}
#toc-toggle .ic,#themeToggle .ic,.modalx .ic,#ghlink .ic{margin:0}
/* first-run orientation lede + expandable-dropdown caret */
.lede{color:var(--mut);font-size:13px;line-height:1.55;margin:0 0 16px;max-width:74ch}
.lede a{color:var(--accent);text-decoration:none} .lede a:hover{text-decoration:underline}
.ddcaret{opacity:.55;transition:transform .15s;margin-left:.3em} details[open]>summary .ddcaret{transform:rotate(180deg)}
.toc-group.closed .toc-items{display:none}
.toc-items{display:flex;flex-direction:column;gap:1px}
.toc-link{color:#5b6b7e;text-decoration:none;font-size:13.5px;padding:6px 12px;border-radius:8px;font-weight:500;border-left:2px solid transparent}
.toc-link:hover{color:var(--accent);background:var(--accent-soft)}
.toc-link.active{color:var(--accent);background:var(--accent-soft);border-left-color:var(--accent);font-weight:600}
section[id]{scroll-margin-top:64px}
@media (max-width:860px){ body{padding-left:0} #toc{box-shadow:0 10px 40px rgba(16,24,40,.18)} }
.wrap{max-width:1180px;margin:0 auto;padding:24px 22px 90px}
section{margin-top:34px} section:first-of-type{margin-top:24px}
h2{font-size:13px;text-transform:uppercase;letter-spacing:.07em;color:var(--mut);font-weight:600;margin:0 0 14px;display:flex;align-items:center;gap:12px;flex-wrap:wrap}
h2 .c{text-transform:none;letter-spacing:0;font-weight:400;font-size:13.5px;color:#5f6f81}
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
.chip{font-size:12.5px;padding:4px 13px;border-radius:20px;border:1px solid var(--line);background:var(--soft);color:var(--txt2);cursor:pointer;user-select:none;transition:.12s;font-weight:500}
.chip:hover{border-color:var(--accent);color:var(--accent);background:var(--accent-soft)} .chip.on{background:var(--accent);color:#fff;border-color:var(--accent)}
.chip .k{opacity:.65;margin-left:5px;font-variant-numeric:tabular-nums}
/* controls */
.controls{display:flex;gap:9px;align-items:center;flex-wrap:wrap;margin-bottom:12px}
.controls input[type=search]{padding:7px 12px;border:1px solid var(--line);border-radius:9px;font-size:13px;min-width:200px;background:var(--panel);box-shadow:var(--sh)}
.controls input[type=search]:focus{border-color:var(--accent)}
.gsearch{padding:4px 10px;border:1px solid var(--line);border-radius:8px;font-size:12px;min-width:110px;max-width:170px;background:var(--panel);box-shadow:var(--sh)}
.gsearch:focus{border-color:var(--accent)} .gsearch:focus:not(:focus-visible){outline:none}
.panel.expanded .gtable{max-height:82vh}
.controls label{font-size:12.5px;color:var(--mut);display:flex;align-items:center;gap:6px;cursor:pointer}
.btn{font-size:12.5px;color:var(--label);border:1px solid var(--line);border-radius:9px;padding:7px 13px;background:var(--panel);cursor:pointer;box-shadow:var(--sh)}
.btn:hover{border-color:var(--accent);color:var(--accent)}
details.dd{position:relative} details.dd>summary{cursor:pointer;font-size:12.5px;color:var(--label);list-style:none;border:1px solid var(--line);border-radius:9px;padding:7px 13px;background:var(--panel);box-shadow:var(--sh)}
details.dd>summary::-webkit-details-marker{display:none} details.dd>summary::marker{content:""}
details.dd .menu{position:absolute;z-index:50;background:var(--panel);border:1px solid var(--line);border-radius:12px;padding:12px 16px;box-shadow:0 14px 40px rgba(16,24,40,.16);max-height:340px;overflow:auto;margin-top:6px;min-width:170px}
details.dd .menu label{display:block;padding:3px 0;font-size:12.5px;color:var(--label)}
.seg{display:inline-flex;border:1px solid var(--line);border-radius:9px;overflow:hidden;box-shadow:var(--sh)}
.seg button{border:0;background:var(--panel);color:var(--txt2);font-size:12.5px;padding:6px 13px;cursor:pointer;border-right:1px solid var(--line)}
.seg button:last-child{border-right:0} .seg button.on{background:var(--accent);color:#fff}
select.msel{font-size:12.5px;border:1px solid var(--line);border-radius:9px;padding:6px 10px;background:var(--panel);color:var(--label);box-shadow:var(--sh)}
.hint{margin-left:auto;color:var(--mut);font-size:12px}
/* table */
.gtable{overflow:auto;max-height:76vh;border-radius:var(--r)}
.gtable::-webkit-scrollbar{height:10px;width:10px} .gtable::-webkit-scrollbar-thumb{background:#cfd8e3;border-radius:6px;border:2px solid #fff}
table{border-collapse:separate;border-spacing:0;width:100%;font-size:13.5px}
th,td{padding:8px 12px;white-space:nowrap;border-bottom:1px solid #eef2f6;text-align:right}
th{background:var(--soft);color:var(--txt2);cursor:pointer;user-select:none;position:sticky;top:0;font-weight:600;font-size:12px;text-transform:uppercase;letter-spacing:.03em;z-index:4;box-shadow:0 1px 0 var(--line)}
th:hover{background:#eef3f8} td.na{color:#c4ccd7}
tbody tr{cursor:pointer;transition:background .1s} tbody tr:hover td{background:var(--soft)}
tr.hl td{background:#fff8e1!important}
th.s,td.s{text-align:left;position:sticky;left:0;background:var(--panel);border-right:1px solid var(--line);z-index:3;font-weight:600}
th.s{z-index:6;background:var(--soft)} tr.hl td.s{background:#fff3ce!important}
.v{font-weight:600;padding:3px 11px;border-radius:20px;font-size:12px;display:inline-block;letter-spacing:.02em}
.v.PASS{color:#0b7350;background:#dff5ec} .v.WARN{color:#95560d;background:#fdefd6} .v.FAIL{color:#a01f2d;background:#fde3e6}
.flags{color:var(--mut);font-size:12.5px;text-align:left;white-space:normal}
/* per-column filter row under the header */
tr.colfilt th{position:sticky;top:33px;background:var(--soft);cursor:auto;text-transform:none;letter-spacing:0;padding:4px 7px;z-index:4;box-shadow:0 1px 0 var(--line)}
tr.colfilt th.s{left:0;z-index:6;background:var(--soft)}
tr.colfilt .cfx{width:100%;min-width:56px;box-sizing:border-box;padding:3px 6px;border:1px solid var(--line);border-radius:6px;font-size:11.5px;font-weight:400;text-transform:none;background:var(--panel);color:var(--ink)}
tr.colfilt .cfx:focus{border-color:var(--accent)} tr.colfilt .cfx:focus:not(:focus-visible){outline:none}
tr.colfilt .cfx::placeholder{color:#aeb8c6}
/* plots */
.bee{display:flex;align-items:center;border-bottom:1px solid #f2f5f9;height:40px} .bee:last-child{border:0}
.bl{flex:0 0 150px;padding:0 14px;font-size:12px;color:var(--txt2);text-align:right;font-weight:500} .nd{color:#c4ccd7;font-size:12px;padding-left:14px}
.plotwrap{flex:1} .plotwrap svg{display:block} .plotwrap circle,.plotwrap rect.hit{cursor:pointer}
.stack{display:flex;align-items:center;border-bottom:1px solid #f4f7fa;height:26px} .stack:last-child{border:0}
.stack .sl{flex:0 0 150px;padding:0 14px;font-size:11.5px;color:var(--txt2);text-align:right;overflow:hidden;text-overflow:ellipsis;white-space:nowrap}
.stack .sb{flex:1;height:13px;display:flex;border-radius:4px;overflow:hidden;background:var(--track);box-shadow:inset 0 0 0 1px rgba(16,24,40,.03)}
.stack .sv{flex:0 0 52px;text-align:right;padding-right:14px;font-size:11px;color:#9aa7b6;font-variant-numeric:tabular-nums}
.legend{display:flex;gap:18px;font-size:11.5px;color:var(--mut);padding:10px 16px;border-top:1px solid var(--line);flex-wrap:wrap}
.legend i{display:inline-block;width:10px;height:10px;border-radius:50%;margin-right:6px;vertical-align:middle}
.footer{color:var(--mut);font-size:11.5px;margin-top:30px;border-top:1px solid var(--line);padding-top:14px;line-height:1.8}
.footer a{color:var(--accent);text-decoration:none} .footer a:hover{text-decoration:underline} .foot-brand{white-space:nowrap}
#tt{position:fixed;pointer-events:none;background:#0f2431;color:#fff;font-size:12px;line-height:1.5;padding:7px 11px;border-radius:9px;opacity:0;transition:opacity .08s;z-index:80;white-space:nowrap;box-shadow:0 8px 24px rgba(15,36,49,.32)}
#toast{position:fixed;left:50%;bottom:26px;transform:translateX(-50%) translateY(8px);background:#0f2431;color:#fff;font-size:12.5px;padding:10px 16px;border-radius:10px;box-shadow:0 12px 34px rgba(15,36,49,.4);opacity:0;pointer-events:none;transition:opacity .18s,transform .18s;z-index:120}
#toast.show{opacity:1;transform:translateX(-50%) translateY(0)}
#infopop{position:fixed;display:none;max-width:320px;background:#0f2431;color:#fff;font-size:13px;line-height:1.55;padding:11px 14px;border-radius:11px;z-index:200;box-shadow:0 12px 34px rgba(15,36,49,.42)}
/* curation basket */
.cbx{vertical-align:middle;margin-right:7px;accent-color:var(--accent);cursor:pointer;width:14px;height:14px}
.sname{cursor:pointer;border-bottom:1px dashed transparent;transition:.12s} .sname:hover{color:var(--accent);border-bottom-color:var(--accent)}
th.s .hlab{font-weight:600}
.curation{display:flex;align-items:center;gap:14px;flex-wrap:wrap;margin-top:12px;padding:12px 16px;background:linear-gradient(180deg,#fff,var(--soft));border:1px solid var(--line);border-radius:12px;box-shadow:var(--sh);scroll-margin-top:64px}
#nbasket:hover{text-decoration:underline!important}
.cur-intro{flex-basis:100%;font-size:11.5px;color:var(--mut);line-height:1.55;margin-bottom:2px} .cur-intro b{color:var(--label)}
.cur-read{font-size:13px;color:var(--label)} .cur-read b{color:var(--ink);font-variant-numeric:tabular-nums;font-size:15px} .cur-read .arw{color:var(--accent);margin:0 3px;font-weight:700}
/* expand-to-fill (fullscreen-within-window) */
.exp-h{display:inline-flex;align-items:center;gap:5px;font-size:11.5px;color:var(--label);border:1px solid var(--line);border-radius:8px;padding:4px 9px;background:var(--panel);cursor:pointer;box-shadow:var(--sh)}
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
.modalx{position:absolute;top:12px;right:14px;border:0;background:var(--track);color:var(--txt2);width:28px;height:28px;border-radius:50%;font-size:18px;line-height:1;cursor:pointer} .modalx:hover{background:#e0e6ee}
.dhead{display:flex;align-items:flex-start;justify-content:space-between;gap:12px;margin-bottom:12px} .dtitle{font-size:18px;font-weight:700;letter-spacing:-.3px} .dmeta{font-size:12px;color:var(--mut);margin-top:2px}
.dflags{display:flex;gap:6px;flex-wrap:wrap;margin-bottom:6px} .dflags .chip{cursor:default}
.chip.failc{background:#fde3e6;color:#a01f2d;border-color:#f5c2c8}   /* FAIL-tier flags read red everywhere (flagged panel + table + modal), not just the modal */
.fchips{display:flex;gap:6px;flex-wrap:wrap}
.chip-hint{font-size:10.5px;color:var(--mut);align-self:center;margin-left:2px;white-space:nowrap}
/* Kraken taxonomic-composition panel */
.krktable td{vertical-align:middle} .krktable .krk-barcell{width:34%;min-width:150px}
.krk-bar{display:flex;height:13px;border-radius:4px;overflow:hidden;box-shadow:inset 0 0 0 1px rgba(16,24,40,.06)}
.krk-seg{height:100%;transition:filter .1s} .krk-seg:hover{filter:brightness(1.1)}
.krk-legend{display:flex;gap:16px;flex-wrap:wrap;align-items:center;font-size:11.5px;color:var(--mut);margin-bottom:12px}
.krk-legend i{display:inline-block;width:10px;height:10px;border-radius:2px;margin-right:5px;vertical-align:middle}
.krk-mut{color:var(--mut);font-size:11px}
.krk-chart{margin-bottom:4px} .krk-plotscroll{overflow-x:auto;overflow-y:hidden}
.krk-plotscroll svg rect{transition:opacity .1s} .krk-plotscroll svg rect:hover{opacity:.82}
/* Executive summary */
.exec-narr{font-size:14.5px;line-height:1.65;color:var(--ink);margin:0 0 18px;max-width:82ch}
.exec-narr b{font-weight:700} .tone-bad{color:var(--fail)} .tone-warn{color:var(--warn)}
.exec-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(155px,1fr));gap:12px;margin-bottom:20px}
.kpi{background:var(--panel);border:1px solid var(--line);border-radius:13px;padding:13px 15px 12px;box-shadow:var(--sh);position:relative;overflow:hidden}
.kpi::before{content:"";position:absolute;left:0;top:0;bottom:0;width:3px;background:var(--line)}
.kpi.good::before{background:var(--pass)} .kpi.warn::before{background:var(--warn)} .kpi.bad::before{background:var(--fail)}
.kpi-l{font-size:10px;font-weight:700;letter-spacing:.07em;text-transform:uppercase;color:var(--mut)}
.kpi-n{font-size:27px;font-weight:800;letter-spacing:-.6px;color:var(--ink);margin:2px 0 3px;font-variant-numeric:tabular-nums;line-height:1.05}
.kpi.good .kpi-n{color:var(--pass)} .kpi.warn .kpi-n{color:var(--warn)} .kpi.bad .kpi-n{color:var(--fail)}
.kpi-s{font-size:11.5px;color:var(--txt2);line-height:1.45} .v.xs{font-size:9.5px;padding:1px 6px;vertical-align:middle}
.exec-cols{display:grid;grid-template-columns:1.25fr 1fr;gap:26px}
.exec-block .exec-h{font-size:11.5px;font-weight:700;letter-spacing:.05em;text-transform:uppercase;color:var(--mut);margin-bottom:11px}
.qcp-row{display:flex;align-items:center;gap:11px;margin-bottom:7px}
.qcp-lab{flex:0 0 128px;font-size:11.5px;color:var(--txt2);font-weight:600;font-variant-numeric:tabular-nums}
.qcp-bar{flex:1;height:8px;background:var(--track);border-radius:5px;overflow:hidden}
.qcp-bar span{display:block;height:100%;border-radius:5px}
.qcp-n{flex:0 0 22px;text-align:right;font-size:12px;font-weight:700;font-variant-numeric:tabular-nums;color:var(--ink)}
@media(max-width:760px){.exec-cols{grid-template-columns:1fr;gap:18px}}
.flagrsn{display:flex;flex-direction:column;gap:1px;margin-top:5px;font-size:11px;color:var(--txt2);font-variant-numeric:tabular-nums}   /* always-visible flag margins (touch-safe, not hover-only) */
.dsub{font-size:10.5px;text-transform:uppercase;letter-spacing:.08em;color:var(--mut);font-weight:600;margin:16px 0 8px}
.lcomp{display:flex;flex-direction:column;gap:5px;margin-bottom:2px}
.lrow{display:flex;align-items:center;gap:9px;font-size:12px} .lk{flex:0 0 58px;color:var(--label);font-weight:600} .lbarw{flex:1;height:9px;background:var(--track);border-radius:5px;overflow:hidden} .lbar{height:100%;background:var(--accent)} .lv{flex:0 0 46px;text-align:right;color:#8895a6;font-variant-numeric:tabular-nums}
.drow{display:flex;align-items:center;gap:10px;padding:3px 0;font-size:12px;border-bottom:1px solid #f4f7fa} .drow:last-of-type{border:0}
.dk{flex:0 0 120px;color:var(--txt2)} .dbarwrap{flex:1;height:7px;background:var(--track);border-radius:4px;overflow:hidden} .dbar{height:100%;border-radius:4px} .dv{flex:0 0 74px;text-align:right;font-weight:600;font-variant-numeric:tabular-nums} .dp{flex:0 0 34px;text-align:right;color:#9aa7b6;font-size:11px}
.dbtns{margin-top:16px;display:flex;justify-content:flex-end}
/* metric help */
#helpmenu{max-width:340px;white-space:normal}
.hitem{padding:6px 0;border-bottom:1px solid #f0f3f7} .hitem:last-child{border:0} .hitem b{font-size:12px} .hk{color:#aab6c4;font-size:10.5px} .hd{font-size:11.5px;color:var(--txt2);margin-top:2px} .hr{color:var(--accent)}
/* provenance bar */
.provbar{display:flex;align-items:center;gap:12px;margin-top:12px}
.prov{display:flex;gap:8px;flex-wrap:wrap;flex:1;font-size:11px;color:var(--mut)}
.prov span{background:var(--panel);border:1px solid var(--line);border-radius:7px;padding:3px 9px;font-variant-numeric:tabular-nums}
/* beeswarm cohort sub-label + bands */
.bl .stat{display:block;font-weight:400;font-size:9.5px;color:#9aa7b6;font-variant-numeric:tabular-nums;margin-top:1px;line-height:1.1;overflow:hidden;text-overflow:ellipsis}
/* lineage */
.ldot{display:inline-block;width:9px;height:9px;border-radius:50%;margin-right:8px;vertical-align:middle}
tr.lingrp td{background:#f0f5f9;color:var(--label);font-weight:600;font-size:11px;letter-spacing:.03em;position:sticky;left:0} tr.lingrp{cursor:default}
.abadge{margin-left:7px;font-size:9px;font-weight:700;letter-spacing:.04em;color:#8a5a12;background:#fdefd6;border:1px solid #f0d9a8;border-radius:6px;padding:1px 5px;vertical-align:middle}
.lincomp-bar{display:flex;height:20px;border-radius:6px;overflow:hidden;margin:14px 16px 6px;box-shadow:inset 0 0 0 1px rgba(16,24,40,.05)}
.lincomp-bar .lseg{height:100%;min-width:2px;transition:filter .1s} .lincomp-bar .lseg:hover{filter:brightness(1.08)}
.lincomp-lab{display:flex;gap:12px;flex-wrap:wrap;padding:0 16px 12px;font-size:11.5px;color:var(--mut)}
.lincomp-lab .lchip i{display:inline-block;width:9px;height:9px;border-radius:2px;margin-right:5px;vertical-align:middle}
.lincomp-lab .lchip b{color:var(--ink);font-variant-numeric:tabular-nums;margin-left:2px}
#linsumtable tbody tr{cursor:pointer} #linsumtable tbody tr:hover td{background:var(--soft)}
/* genome landscape */
.genome-scroll{overflow-x:auto;overflow-y:hidden} .genome-scroll svg{display:block;cursor:crosshair} .genome-scroll rect{shape-rendering:crispEdges}
#scatter svg{cursor:crosshair} #corr_body{overflow-x:auto} #corrsvg text{pointer-events:none}
.hot-note{font-size:11.5px;color:var(--mut);line-height:1.5;padding:12px 16px 6px}
#hottable tbody tr{cursor:pointer} #hottable tbody tr:hover td{background:var(--soft)}
/* aDNA panel */
.anote{font-size:12px;color:var(--txt2);line-height:1.6;background:var(--soft);border:1px solid var(--line);border-radius:10px;padding:10px 13px;margin-bottom:12px} .anote b{color:var(--ink)}
.acards{display:grid;grid-template-columns:repeat(auto-fill,minmax(230px,1fr));gap:11px}
.acard{border:1px solid var(--line);border-radius:12px;padding:11px 13px;background:linear-gradient(180deg,#fff,var(--soft))}
.acard.low{border-color:#f0c2be;background:linear-gradient(180deg,#fff,#fdf1f0)}
.acard-h{display:flex;align-items:center;justify-content:space-between;gap:8px;margin-bottom:7px} .acard-h .sname{font-weight:600;font-size:12.5px}
.acard-row{display:flex;align-items:center;gap:12px} .aspark{flex:0 0 auto} .astats{display:flex;flex-direction:column;gap:3px;flex:1}
.astat{display:flex;justify-content:space-between;font-size:11.5px} .ak{color:#8895a6} .av{font-weight:600;font-variant-numeric:tabular-nums}
.av.good{color:#0b7350} .av.bad{color:#a01f2d} .av.na{color:#c4ccd7}
.alow{margin-top:8px;font-size:11px;color:#a01f2d;line-height:1.45}
/* modal 'why' block */
.dwhy{font-size:11.5px;color:var(--txt2);background:var(--soft);border:1px solid var(--line);border-radius:9px;padding:8px 11px;margin-bottom:6px;line-height:1.55} .dwhy b{color:#a01f2d}
@media print{
  *,*::before,*::after{-webkit-print-color-adjust:exact!important;print-color-adjust:exact!important;color-adjust:exact!important}
  header,nav,#toc,#toc-toggle,.controls,.provbar #printBtn,#thbox,#athbox,.dd,.chips{display:none!important}
  #modal{display:none!important}   /* never capture an open detail dialog over the printed page */
  body{background:#ffffff;padding-left:0} .wrap{max-width:none;padding:0}
  .gtable,.snpmx-wrap,.epimx-wrap,.dr-mxwrap,.epitbl-wrap,#stacks,#fn_stacks{max-height:none!important;overflow:visible!important}
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
  header{flex-wrap:wrap;row-gap:6px}
  /* larger touch targets on phones */
  .controls label{padding:4px 2px} .controls input[type=checkbox],tr.colfilt input[type=checkbox]{width:17px;height:17px}
  .chip{padding:8px 13px} .seg button{padding:8px 13px} .exp-h{padding:7px 11px}
  .modalx{width:38px;height:38px;font-size:22px} .btn{padding:8px 13px}
  #themeToggle,#ghlink{width:44px;height:44px} #toc-toggle{width:44px;height:44px}
  /* the frozen (sticky left:0) sample column is opaque and covers the viewport; long sample
     names would otherwise push the metric columns off-screen and out of reach, so cap it */
  td.s .sname,th.s .hlab{display:inline-block;max-width:40vw;overflow:hidden;text-overflow:ellipsis;white-space:nowrap;vertical-align:middle}
  #flagtable td.s{white-space:normal;overflow-wrap:anywhere;word-break:break-word;max-width:48vw}
  .dr-row{max-width:150px;overflow:hidden;text-overflow:ellipsis}
}
@media (max-width:860px){ body:not(.toc-collapsed) #tocscrim{display:block} }
@media (prefers-reduced-motion:reduce){*{transition:none!important;scroll-behavior:auto!important;animation:none!important}}
/* visible info icon signalling a hover tooltip */
.infoi{display:inline-flex;align-items:center;justify-content:center;width:15px;height:15px;border-radius:50%;background:#d5deea;color:var(--txt2);font-size:10px;font-weight:700;font-style:italic;font-family:Georgia,'Times New Roman',serif;text-transform:none;margin-left:5px;cursor:help;vertical-align:middle;line-height:1;transition:.12s;user-select:none}
.infoi:hover{background:var(--accent);color:#fff}
th .infoi,.dyn-legend .infoi{background:#dde5f0}
/* SNP dynamics panel (gene-centric, searchable) */
.dyn-controls{display:flex;gap:10px;align-items:center;flex-wrap:wrap;margin-bottom:12px}
.dyn-search{padding:9px 14px;border:1px solid var(--line);border-radius:10px;font-size:14px;min-width:240px;background:var(--panel);box-shadow:var(--sh)}
.dyn-btn{font-size:13px;color:var(--label);border:1px solid var(--line);border-radius:9px;padding:8px 14px;background:var(--panel);cursor:pointer;font-weight:500}
.dyn-btn:hover{background:var(--accent-soft);color:var(--accent);border-color:var(--accent)}
.dyn-count{font-size:13px;color:var(--mut);margin-left:auto}
.dyn-legend{display:flex;gap:18px;flex-wrap:wrap;align-items:center;font-size:13px;color:var(--mut);margin-bottom:14px}
.dyn-legend i{display:inline-block;width:12px;height:12px;border-radius:3px;margin-right:6px;vertical-align:-1px}
.dyn-genechips{display:flex;gap:9px;flex-wrap:wrap;margin-bottom:20px;max-height:140px;overflow:auto;padding:2px}
.dyn-chip{font-size:13.5px;padding:6px 13px;border-radius:22px;border:1px solid var(--line);background:var(--soft);color:var(--txt2);cursor:pointer;display:inline-flex;align-items:center;gap:6px;transition:.12s}
.dyn-chip:hover{border-color:var(--accent);color:var(--accent)}
.dyn-chip b{font-weight:700;color:#98a6b8}
.dyn-chip.sel{background:var(--accent);border-color:var(--accent);color:#fff} .dyn-chip.sel:hover{color:#fff} .dyn-chip.sel b{color:#d7e6ff}
.dyn-dot{width:8px;height:8px;border-radius:50%;background:#d1495b;display:inline-block}
.dyn-grid{display:grid;grid-template-columns:repeat(auto-fill,minmax(var(--dyncw,250px),1fr));gap:14px;align-items:start}
.dyn-card{border:1px solid var(--line);border-radius:13px;padding:11px 13px;background:var(--soft);display:flex;flex-direction:column}
.dyn-card.flagged{border-color:#cdd9ea;box-shadow:0 1px 0 rgba(31,120,180,.04)}
.dyn-cardgene-dot{width:7px;height:7px;border-radius:50%;background:#d1495b;display:inline-block;margin-left:5px;vertical-align:1px}
.dyn-card-h{display:flex;justify-content:space-between;align-items:center;gap:8px}
.dyn-cardgene{font-weight:800;font-size:13.5px;color:var(--ink)}
.dyn-pos{font-weight:600;font-size:11.5px;color:var(--mut);font-variant-numeric:tabular-nums}
.dyn-grp{font-size:11.5px;color:var(--txt2);background:var(--panel);border:1px solid var(--line);border-radius:7px;padding:1px 8px}
.dyn-eff{font-size:12px;color:var(--txt2);margin:3px 0 7px;line-height:1.35}
.dyn-zoom{display:flex;align-items:center;gap:7px;font-size:12.5px;color:var(--mut)}
.dyn-zoom input[type=range]{cursor:pointer;accent-color:var(--accent);width:118px}
.dyn-toggle{display:flex;align-items:center;gap:5px;font-size:12.5px;color:var(--mut);cursor:pointer}
.dyn-toggle input{accent-color:#5b8fc9;cursor:pointer}
.dyn-aa{color:var(--accent);font-weight:700;font-variant-numeric:tabular-nums}
.rvtag{font-size:10.5px;color:#8a97a8;text-decoration:none;font-weight:600;white-space:nowrap}
.rvtag:hover{color:var(--accent);text-decoration:underline}
.aah37{color:#8a97a8;font-weight:600;font-size:.9em}
.dyn-card-f{display:flex;flex-wrap:wrap;gap:6px;align-items:center;margin-top:8px}
.dyn-fchip{color:#fff;border-radius:7px;padding:2px 9px;font-size:12px;font-weight:500}
.dyn-traj{font-size:12.5px;color:var(--txt2);font-variant-numeric:tabular-nums;margin-left:auto}
.dyn-empty{padding:36px;text-align:center;color:var(--mut);font-size:14.5px;border:1px dashed var(--line);border-radius:14px}
/* Epistasis (co-dynamics pairs) */
.epi-controls{display:flex;gap:10px;align-items:center;flex-wrap:wrap;margin-bottom:12px}
.dyn-btn.on{background:var(--accent-soft);color:var(--accent);border-color:var(--accent)}
.epi-legend{display:flex;gap:20px;flex-wrap:wrap;align-items:center;font-size:13px;color:var(--mut);margin-bottom:14px}
.epi-legend i{display:inline-block;width:12px;height:12px;border-radius:3px;margin-right:6px;vertical-align:-1px}
.epi-card{border:1px solid var(--line);border-radius:13px;padding:11px 13px;background:var(--soft);display:flex;flex-direction:column}
.epi-card.discordant{border-color:#e7d2e6}
.epi-card-h{display:flex;align-items:center;gap:8px;margin-bottom:5px}
.epi-badge{font-size:12px;font-weight:700;border-radius:7px;padding:2px 9px;color:#fff;font-variant-numeric:tabular-nums}
.epi-badge.concordant{background:#2f8f5b}
.epi-badge.discordant{background:#a24a8f}
.epi-dir{font-size:11.5px;color:var(--mut)}
.epi-pair{font-size:12.5px;color:var(--txt2);margin-bottom:3px;line-height:1.45}
.epi-vs{color:var(--mut);font-weight:700;margin:0 4px}
.epi-card-f{display:flex;justify-content:space-between;align-items:center;gap:8px;font-size:11.5px;color:var(--mut);margin-top:6px}
.epi-flabel{font-size:12px;color:var(--mut);font-weight:600;margin-left:6px}
.epi-tier{font-size:10.5px;font-weight:700;border-radius:6px;padding:2px 8px;text-transform:uppercase;letter-spacing:.3px;margin-left:auto}
.epi-strong{background:#1f7a4d;color:#fff}
.epi-moderate{background:#e7eef7;color:#3f6fa8}
.epi-weak{background:#eef1f5;color:#93a0b0}
.epi-recur{color:#3f6fa8;font-weight:600}
.epi-views{display:inline-flex;border:1px solid var(--line);border-radius:10px;overflow:hidden;margin-bottom:12px}
.epi-viewbtn{font-size:13px;padding:7px 16px;background:var(--panel);border:0;border-right:1px solid var(--line);color:var(--txt2);cursor:pointer;font-weight:500}
.epi-viewbtn:last-child{border-right:0}
.epi-viewbtn.on{background:var(--accent);color:#fff}
.epimx-note{font-size:12px;color:var(--mut);margin-bottom:8px}
.epimx-wrap{overflow:auto;max-height:72vh;border:1px solid var(--line);border-radius:12px}
table.epimx{border-collapse:separate;border-spacing:0;font-size:11px}
table.epimx th{position:sticky;background:var(--soft);z-index:2}
.epimx-corner{left:0;top:0;z-index:4}
.epimx-hcell{top:0;height:78px;vertical-align:bottom;padding:3px 0;z-index:3}
.epimx-h{writing-mode:vertical-rl;transform:rotate(180deg);white-space:nowrap;font-size:10.5px;font-weight:600;color:var(--label)}
.epimx-row{left:0;text-align:right;padding:2px 9px;font-size:10.5px;font-weight:600;color:var(--label);white-space:nowrap}
.epimx-cell{min-width:26px;height:22px;text-align:center;border-bottom:1px solid #fff;border-right:1px solid #fff;font-size:9px;color:var(--label);font-variant-numeric:tabular-nums}
.epimx-diag{background:repeating-linear-gradient(45deg,#e6ebf2,#e6ebf2 3px,#eef2f7 3px,#eef2f7 6px)}
.epimx-scale{display:flex;align-items:center;gap:8px;font-size:12px;color:var(--mut);margin-top:12px}
.epimx-grad{width:170px;height:12px;border-radius:3px;background:linear-gradient(90deg,#a24a8f,#f3f5f8,#2f8f5b);display:inline-block}
.epitbl-top{display:flex;align-items:center;gap:12px;margin-bottom:10px}
.epitbl-wrap{overflow-x:auto;border:1px solid var(--line);border-radius:12px}
table.epitbl{border-collapse:collapse;width:100%;font-size:12.5px}
table.epitbl th{text-align:left;padding:8px 11px;background:var(--soft);border-bottom:2px solid var(--line);cursor:pointer;color:var(--label);white-space:nowrap;position:sticky;top:0}
table.epitbl td{padding:6px 11px;border-bottom:1px solid #eef2f6;white-space:nowrap}
table.epitbl tbody tr:hover td{background:var(--soft)}
.epitbl-r{font-weight:700;font-variant-numeric:tabular-nums}
/* Drug resistance panel */
.dr-legend{display:flex;gap:18px;flex-wrap:wrap;align-items:center;font-size:13px;color:var(--mut);margin-bottom:12px}
.dr-legend i{display:inline-block;width:12px;height:12px;border-radius:3px;margin-right:6px;vertical-align:-1px}
.dr-mxwrap{overflow:auto;max-height:70vh;border:1px solid var(--line);border-radius:12px;margin-bottom:14px}
table.drmx{border-collapse:separate;border-spacing:0;font-size:12px}
table.drmx th{position:sticky;background:var(--soft);z-index:2}
.dr-corner{left:0;top:0;z-index:4;text-align:right;padding:2px 9px;color:var(--mut)}
.dr-hcell{top:0;height:74px;vertical-align:bottom;padding:3px 0;z-index:3}
.dr-h{writing-mode:vertical-rl;transform:rotate(180deg);white-space:nowrap;font-size:11px;font-weight:700;color:var(--label)}
.dr-row{left:0;text-align:right;padding:3px 10px;font-size:11.5px;font-weight:600;color:var(--label);white-space:nowrap}
.drmx-cell{min-width:26px;height:24px;text-align:center;border-bottom:1px solid #fff;border-right:1px solid #fff}
.drmx-cell b{color:#fff;font-size:11px}
.dr-controls{display:flex;gap:10px;align-items:center;flex-wrap:wrap;margin-bottom:10px}
.dr-badge{color:#fff;border-radius:6px;padding:2px 8px;font-size:11px;font-weight:600;white-space:nowrap}
/* SNP matrix (explorable heatmap) */
.snpmx-controls{display:flex;gap:12px;align-items:center;flex-wrap:wrap;margin-bottom:10px}
.snpmx-toggle{font-size:12.5px;color:var(--mut);display:flex;align-items:center;gap:5px;cursor:pointer}
.snpmx-metanote{font-size:12px;color:var(--mut);margin-bottom:10px}
.snpmx-filters,.dyn-filters{display:flex;gap:12px;align-items:center;flex-wrap:wrap;margin-bottom:10px;padding:9px 13px;background:var(--soft);border:1px solid var(--line);border-radius:10px}
.dyn-filters{margin-bottom:14px}
.snpmx-flabel{font-size:12.5px;color:var(--mut);font-weight:600}
.snpmx-fsel{font-size:12.5px;color:var(--label);display:flex;align-items:center;gap:5px}
.snpmx-fsel select{font-size:12.5px;border:1px solid var(--line);border-radius:8px;padding:4px 9px;background:var(--panel);color:var(--label);cursor:pointer}
.snpmx-wrap{overflow:auto;max-height:74vh;border:1px solid var(--line);border-radius:12px}
table.snpmx{border-collapse:separate;border-spacing:0;font-size:12px;width:auto;margin:0 auto}
table.snpmx th,table.snpmx td{border-bottom:1px solid #eef2f6}
table.snpmx thead th{position:sticky;background:var(--soft);z-index:5}   /* top offset set inline per header row */
table.snpmx .snpmx-metacell{font-size:10px;text-align:center;color:var(--label);padding:2px 4px;height:22px;white-space:nowrap;overflow:hidden;text-overflow:ellipsis;max-width:70px;border-bottom:1px solid #fff}
table.snpmx .snpmx-metalabel{font-size:10.5px;font-weight:600;color:var(--txt2);text-align:right;padding:2px 12px;height:22px}
table.snpmx .snpmx-hcell{padding:4px 1px;vertical-align:bottom;height:104px}
table.snpmx .snpmx-h{writing-mode:vertical-rl;transform:rotate(180deg);font-size:11px;color:var(--txt2);font-weight:600;white-space:nowrap;display:inline-block;max-height:96px;overflow:hidden;text-overflow:ellipsis}
table.snpmx .snpmx-info{position:sticky;left:0;background:var(--panel);z-index:4;text-align:left;padding:5px 13px;border-right:1px solid var(--line);white-space:nowrap;font-size:12.5px}
table.snpmx thead .snpmx-info{z-index:6;background:var(--soft)}
.snpmx-aa{color:var(--accent);font-weight:700}
table.snpmx td.snpmx-cell{min-width:32px;text-align:center;color:#0f2431;font-variant-numeric:tabular-nums;padding:3px 2px}
table.snpmx td.snpmx-cell.wdp{min-width:42px}
.snpmx-af{display:block;font-weight:600;font-size:10px}
.snpmx-dp{display:block;font-size:8.5px;color:var(--txt2);line-height:1.15}
table.snpmx td.snpmx-empty{background:repeating-linear-gradient(45deg,#f6f8fb,#f6f8fb 3px,#eef2f7 3px,#eef2f7 6px)}
table.snpmx tbody tr:hover td.snpmx-info{background:var(--soft)}
/* leading search glyph on search fields (background SVG). Placed last so it wins over the
   `background:var(--panel)` shorthands above, which would otherwise reset background-image. */
.dyn-search,.gsearch,.controls input[type=search]{
  background-image:url("data:image/svg+xml,%3Csvg xmlns='http://www.w3.org/2000/svg' width='16' height='16' viewBox='0 0 24 24' fill='none' stroke='%238895a6' stroke-width='2' stroke-linecap='round' stroke-linejoin='round'%3E%3Ccircle cx='11' cy='11' r='8'/%3E%3Cpath d='m21 21-4.3-4.3'/%3E%3C/svg%3E");
  background-repeat:no-repeat;background-position:11px 50%;background-size:15px;padding-left:33px}
html.dark .dyn-search,html.dark .gsearch,html.dark .controls input[type=search]{
  background-image:url("data:image/svg+xml,%3Csvg xmlns='http://www.w3.org/2000/svg' width='16' height='16' viewBox='0 0 24 24' fill='none' stroke='%236f7f92' stroke-width='2' stroke-linecap='round' stroke-linejoin='round'%3E%3Ccircle cx='11' cy='11' r='8'/%3E%3Cpath d='m21 21-4.3-4.3'/%3E%3C/svg%3E");
  background-repeat:no-repeat;background-position:11px 50%;background-size:15px}   /* re-assert: the html.dark `background:var(--soft)` shorthand above resets repeat, tiling the glyph */
"""

JS = r"""
(function(){
var R=REPORT;
function assign(t,s){for(var _k in s){if(Object.prototype.hasOwnProperty.call(s,_k))t[_k]=s[_k];}return t;}  // ES5 shallow copy
function zeros(n){var _a=new Array(n);for(var _i=0;_i<n;_i++)_a[_i]=0;return _a;}                              // ES5 zero-filled array
var DIST=R.dist;
var VCOL={PASS:'#94a3b8',WARN:'#d97706',FAIL:'#dc2626'}, VFILL={PASS:'#16a34a',WARN:'#d97706',FAIL:'#dc2626'};
// Theme colours for the JS-drawn SVG panels (plots, heatmaps): swapped when the dark theme toggles.
function isDark(){return document.documentElement.classList.contains('dark');}
var TH_LIGHT={ink:'#15202e',mut:'#6a7889',panel:'#ffffff',soft:'#f6f8fb',line:'#e6ebf2',grid:'#e6ebf1',axis:'#8895a6',cellnull:'#e9edf2',track:'#edf1f6',hl:'#fff8e1',faint:'#5b6b7e'};
var TH_DARK ={ink:'#e7edf4',mut:'#94a4b6',panel:'#18232f',soft:'#1e2a38',line:'#2a3543',grid:'#2a3543',axis:'#6f7f92',cellnull:'#222e3c',track:'#232f3d',hl:'#33371c',faint:'#9fb0c2'};
var TH=isDark()?TH_DARK:TH_LIGHT;
var BAR={hi_good:'#22a06b',hi_bad:'#e0544f',neu:'#4f83c2'};
// ---- lineage palette: deterministic, colour-blind-safe, self-contained (Okabe-Ito + Tol-muted; golden-angle overflow)
var LINPAL_BASE=['#4477aa','#ee6677','#228833','#ccbb44','#66ccee','#aa3377','#e69f00','#0072b2','#d55e00','#009e73','#cc79a7','#882255'];
R.lin_colors=R.lin_colors||{};   // canonical palette (e.g. mycolorsTB) keyed by lineage label
var LINCOL={};(function(){(R.lineages||[]).forEach(function(lab,i){
  if(R.lin_colors[lab]){LINCOL[lab]=R.lin_colors[lab];}                       // canonical colour if provided
  else if(i<LINPAL_BASE.length){LINCOL[lab]=LINPAL_BASE[i];}                  // else colour-blind-safe fallback
  else{var h=(i*137.508)%360;LINCOL[lab]='hsl('+h.toFixed(0)+',58%,52%)';}});})();
function linColor(lab){return (lab!=null&&LINCOL[lab])?LINCOL[lab]:'#b8c2cf';}
R.gene_map=R.gene_map||{};
function geneRv(g){ return (g&&R.gene_map[g])||''; }   // Mycobrowser (H37Rv) locus tag for a gene, or ''
function geneRvTag(g){ var rv=geneRv(g); return rv?(' <a class="rvtag" href="https://mycobrowser.epfl.ch/genes/'+esc(rv)+'" target="_blank" rel="noopener" title="Mycobrowser locus tag of '+esc(g)+' (opens mycobrowser.epfl.ch)">'+esc(rv)+'</a>'):''; }
// amino-acid change in the used-reference numbering, plus the H37Rv/Mycobrowser one in brackets when it differs
var AA2LBL=R.aa2_label||'H37Rv';   // label for the canonical-reference amino-acid numbering
function aaDual(aa,aaH){ if(!aa) return ''; return esc(aa)+((aaH&&aaH!==aa)?(' <span class="aah37" title="same variant in the '+esc(AA2LBL)+' reference numbering">['+esc(AA2LBL)+' '+esc(aaH)+']</span>'):''); }
var thr=assign({},R.thresholds);
var athr=assign({},R.anc_thresholds||{});      // ancient (aDNA) threshold view
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
    return ring+'<circle cx="'+cx+'" cy="'+cy+'" r="'+(big?5.4:3.4)+'" fill="'+dotColor(s)+'" opacity="0.85"'+((s.v=='FAIL'||big)?' stroke="'+TH.ink+'" stroke-width="'+(big?1.4:0.6)+'"':'')+' data-s="'+esc(s.s)+'" data-lin="'+esc(s.lineage||'')+'" data-x="'+P.scores[i][0]+'" data-y="'+P.scores[i][1]+'" data-xl="PC1" data-yl="PC2" data-xk="float" data-yk="float"/>';}).join('');
  var frame='<line x1="'+pad+'" y1="'+(H-pad)+'" x2="'+(W-14)+'" y2="'+(H-pad)+'" stroke="'+TH.axis+'"/><line x1="'+pad+'" y1="14" x2="'+pad+'" y2="'+(H-pad)+'" stroke="'+TH.axis+'"/>';
  var zero='';if(xr[0]<0&&xr[1]>0)zero+='<line x1="'+sx(0).toFixed(1)+'" y1="14" x2="'+sx(0).toFixed(1)+'" y2="'+(H-pad)+'" stroke="#eef2f6"/>';if(yr[0]<0&&yr[1]>0)zero+='<line x1="'+pad+'" y1="'+sy(0).toFixed(1)+'" x2="'+(W-14)+'" y2="'+sy(0).toFixed(1)+'" stroke="#eef2f6"/>';
  var xt='<text x="'+((pad+W-14)/2)+'" y="'+(H-6)+'" text-anchor="middle" font-size="11" fill="'+TH.mut+'">PC1 ('+P.pev[0].toFixed(1)+'%)</text>';
  var yt='<text transform="rotate(-90 13 '+((14+H-pad)/2)+')" x="13" y="'+((14+H-pad)/2)+'" text-anchor="middle" font-size="11" fill="'+TH.mut+'">PC2 ('+(P.pev[1]<1e-3?'~0%, rank-deficient':P.pev[1].toFixed(1)+'%')+')</text>';
  var load='<div style="display:flex;gap:24px;flex-wrap:wrap;margin-top:8px">'+[0,1].map(function(pc){return '<div style="flex:1;min-width:170px"><div class="dsub" style="margin:2px 0 6px">PC'+(pc+1)+' loadings</div>'+P.load[pc].map(function(l){var w=Math.abs(l.w),col=l.w>=0?'#22a06b':'#e0544f';return '<div class="drow"><span class="dk">'+esc(l.label)+'</span><div class="dbarwrap"><div class="dbar" style="width:'+Math.round(w*100)+'%;background:'+col+'"></div></div><span class="dv">'+(l.w>=0?'+':'')+l.w.toFixed(2)+'</span></div>';}).join('')+'</div>';}).join('')+'</div>';
  host.innerHTML='<svg viewBox="0 0 '+W+' '+H+'" style="width:100%;height:auto;display:block;cursor:crosshair">'+ellSVG+frame+zero+xt+yt+dots+'</svg>'+colorLegend()+load;
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
    return ring+'<circle cx="'+cx+'" cy="'+cy+'" r="'+(big?5.4:3.4)+'" fill="'+dotColor(s)+'" opacity="0.82"'+((s.v=='FAIL'||big)?' stroke="'+TH.ink+'" stroke-width="'+(big?1.4:0.6)+'"':'')+' data-s="'+esc(s.s)+'" data-lin="'+esc(s.lineage||'')+'" data-x="'+xv+'" data-y="'+yv+'" data-xl="'+(st.divx=='callable_inv'?'100 - callable %':'Missing %')+'" data-yl="'+esc(ylabel)+'" data-xk="pct" data-yk="'+(R.snp_density_ok?'float':'int')+'"/>';}).join('');
  var frame='<line x1="'+pad+'" y1="'+(H-pad)+'" x2="'+(W-14)+'" y2="'+(H-pad)+'" stroke="'+TH.axis+'"/><line x1="'+pad+'" y1="14" x2="'+pad+'" y2="'+(H-pad)+'" stroke="'+TH.axis+'"/>';
  var titles='<text x="'+((pad+W-14)/2)+'" y="'+(H-6)+'" text-anchor="middle" font-size="11" fill="'+TH.mut+'">'+(st.divx=='callable_inv'?'100 - callable % (incompleteness)':'Missing % (incompleteness)')+'</text><text transform="rotate(-90 13 '+((14+H-pad)/2)+')" x="13" y="'+((14+H-pad)/2)+'" text-anchor="middle" font-size="11" fill="'+TH.mut+'">'+esc(ylabel)+'</text>';
  var labs='<text x="'+(pad+6)+'" y="'+(H-pad-6)+'" font-size="8.5" font-weight="600" fill="var(--fail)">reference-bias suspect</text>'
    +'<text x="'+(pad+6)+'" y="24" font-size="8.5" font-weight="600" fill="#3f7d55">typical divergence</text>'
    +'<text x="'+(W-16)+'" y="'+(H-pad-6)+'" text-anchor="end" font-size="8.5" font-weight="600" fill="var(--warn)">low coverage</text>';
  host.innerHTML='<svg viewBox="0 0 '+W+' '+H+'" style="width:100%;height:auto;display:block;cursor:crosshair">'+rects+guides+frame+titles+labs+dots+'</svg>'+colorLegend();
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
  var ticks='';for(var t=0;t<=4;t++){var yr=Math.round(ymin+span*t/4),X=ax(yr);ticks+='<line x1="'+X.toFixed(1)+'" y1="14" x2="'+X.toFixed(1)+'" y2="18" stroke="'+TH.axis+'"/><text x="'+X.toFixed(1)+'" y="11" text-anchor="middle" font-size="8.5" fill="#94a3b8">'+yr+'</text>';}
  var axisSVG='<svg viewBox="0 0 '+W+' 22" style="width:100%;height:auto;display:block"><line x1="'+gut+'" y1="18" x2="'+(gut+plotW)+'" y2="18" stroke="'+TH.grid+'"/>'+ticks+'</svg>';
  var body=rows.map(function(L){var grp=g[L],ys=grp.map(function(s){return yearOf(s.date);}).filter(function(y){return y!=null;});
    var snps=grp.map(function(s){return s.m.snps;}).filter(function(v){return v!=null;}).sort(function(a,b){return a-b;});
    var titv=grp.map(function(s){return s.m.ti_tv;}).filter(function(v){return v!=null;}).sort(function(a,b){return a-b;});
    var mn=ys.length?Math.min.apply(null,ys):null,mxx=ys.length?Math.max.apply(null,ys):null;
    var col=R.lin_present?linColor(L):'var(--accent)';
    var rug=grp.map(function(s,i){var y=yearOf(s.date);if(y==null)return '';var jx=(((i*2654435761)>>>0)%997)/997-0.5,big=(st.hi==s.s);
      return '<circle cx="'+ax(y).toFixed(1)+'" cy="'+(13+jx*8).toFixed(1)+'" r="'+(big?4.5:2.4)+'" fill="'+col+'" opacity="0.72"'+(big?' stroke="'+TH.ink+'" stroke-width="1"':'')+' data-s="'+esc(s.s)+'" data-lin="'+esc(s.lineage||'')+'" data-x="'+y+'" data-y="'+(s.m.snps!=null?s.m.snps:0)+'" data-xl="year" data-yl="SNPs" data-xk="int" data-yk="int"/>';}).join('');
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
var st={sortKey:'v',asc:false,q:'',onlyFlagged:false,hidden:{},hi:null,flagFilter:null,ptype:'beeswarm',   // default view: worst QC first (a returning user's saved sort overrides this)
        sx:'mean_depth',sy:'breadth_pct',excl:{},detail:null,colorBy:'qc',groupLin:false,ancOnly:null,linFilter:null,gtrack:'missing',maskOn:false,gsel:null,gzoom:null,gbq:'',hotq:'',pnpsq:'',colf:{},showColF:false};
var SGEO=null, GGEO=null;   // scatter + genome brush geometry caches (for inverse-mapping the rubber-band)
var _thdb,_qdb,_mxdb,_cfdb;  // debounce timers: keep live inputs snappy at cohort scale (defer heavy re-renders)
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
// ---- Icon set (Lucide, MIT): one homogeneous stroke family; currentColor -> auto light/dark ----
var IC={
  printer:'<path d="M6 18H4a2 2 0 0 1-2-2v-5a2 2 0 0 1 2-2h16a2 2 0 0 1 2 2v5a2 2 0 0 1-2 2h-2"/><path d="M6 9V3a1 1 0 0 1 1-1h10a1 1 0 0 1 1 1v6"/><rect x="6" y="14" width="12" height="8" rx="1"/>',
  sun:'<circle cx="12" cy="12" r="4"/><path d="M12 2v2"/><path d="M12 20v2"/><path d="m4.9 4.9 1.4 1.4"/><path d="m17.7 17.7 1.4 1.4"/><path d="M2 12h2"/><path d="M20 12h2"/><path d="m6.3 17.7-1.4 1.4"/><path d="m19.1 4.9-1.4 1.4"/>',
  moon:'<path d="M12 3a6 6 0 0 0 9 9 9 9 0 1 1-9-9Z"/>',
  panel:'<rect width="18" height="18" x="3" y="3" rx="2"/><path d="M9 3v18"/>',
  sliders:'<line x1="21" x2="14" y1="4" y2="4"/><line x1="10" x2="3" y1="4" y2="4"/><line x1="21" x2="12" y1="12" y2="12"/><line x1="8" x2="3" y1="12" y2="12"/><line x1="21" x2="16" y1="20" y2="20"/><line x1="12" x2="3" y1="20" y2="20"/><line x1="14" x2="14" y1="2" y2="6"/><line x1="8" x2="8" y1="10" y2="14"/><line x1="16" x2="16" y1="18" y2="22"/>',
  maximize:'<polyline points="15 3 21 3 21 9"/><polyline points="9 21 3 21 3 15"/><line x1="21" x2="14" y1="3" y2="10"/><line x1="3" x2="10" y1="21" y2="14"/>',
  minimize:'<polyline points="4 14 10 14 10 20"/><polyline points="20 10 14 10 14 4"/><line x1="14" x2="21" y1="10" y2="3"/><line x1="3" x2="10" y1="21" y2="14"/>',
  basket:'<path d="m15 11-1 9"/><path d="m19 11-4-7"/><path d="M2 11h20"/><path d="m3.5 11 1.6 7.4a2 2 0 0 0 2 1.6h9.8a2 2 0 0 0 2-1.6l1.6-7.4"/><path d="m5 11 4-7"/><path d="m9 11 1 9"/>',
  download:'<path d="M21 15v4a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2v-4"/><polyline points="7 10 12 15 17 10"/><line x1="12" x2="12" y1="15" y2="3"/>',
  search:'<circle cx="11" cy="11" r="8"/><path d="m21 21-4.3-4.3"/>',
  x:'<path d="M18 6 6 18"/><path d="m6 6 12 12"/>',
  chevronDown:'<path d="m6 9 6 6 6-6"/>',
  chevronUp:'<path d="m18 15-6-6-6 6"/>',
  help:'<circle cx="12" cy="12" r="10"/><path d="M9.09 9a3 3 0 0 1 5.83 1c0 2-3 3-3 3"/><path d="M12 17h.01"/>',
  alert:'<path d="m21.73 18-8-14a2 2 0 0 0-3.48 0l-8 14A2 2 0 0 0 4 21h16a2 2 0 0 0 1.73-3Z"/><path d="M12 9v4"/><path d="M12 17h.01"/>',
  ban:'<circle cx="12" cy="12" r="10"/><path d="m4.9 4.9 14.2 14.2"/>',
  github:'<path d="M15 22v-4a4.8 4.8 0 0 0-1-3.5c3 0 6-2 6-5.5.08-1.25-.27-2.48-1-3.5.28-1.15.28-2.35 0-3.5 0 0-1 0-3 1.5-2.64-.5-5.36-.5-8 0C6 2 5 2 5 2c-.3 1.15-.3 2.35 0 3.5A5.4 5.4 0 0 0 4 9c0 3.5 3 5.5 6 5.5-.39.49-.68 1.05-.85 1.65-.17.6-.22 1.23-.15 1.85v4"/><path d="M9 18c-4.51 2-5-2-7-2"/>',
  ext:'<path d="M15 3h6v6"/><path d="M10 14 21 3"/><path d="M18 13v6a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2V8a2 2 0 0 1 2-2h6"/>'
};
function icon(n,cls){return '<svg class="ic'+(cls?' '+cls:'')+'" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true" focusable="false">'+(IC[n]||'')+'</svg>';}
function fillIcons(root){Array.prototype.forEach.call((root||document).querySelectorAll('[data-ic]'),function(e){e.innerHTML=icon(e.getAttribute('data-ic'),e.getAttribute('data-ic-cls')||'');e.removeAttribute('data-ic');});}
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
// Per-column filters for the General Statistics table (table-scoped; do NOT touch the global visible()
// so the plots stay driven by the global filters). A filter string is a numeric operator/range on numeric
// columns (>50, >=50, <10, 5-9, 5..9), otherwise a case-insensitive substring on the displayed cell text.
function colMatchOne(raw,val,txt){
  var q=(raw||'').trim(); if(!q) return true;
  if(typeof val=='number'&&!isNaN(val)){
    var m=q.match(/^(>=|<=|>|<|=)\s*(-?\d+(?:\.\d+)?)$/);
    if(m){var n=parseFloat(m[2]),o=m[1];
      if(o=='>')return val>n; if(o=='>=')return val>=n; if(o=='<')return val<n; if(o=='<=')return val<=n; return val==n;}
    m=q.match(/^(-?\d+(?:\.\d+)?)\s*(?:\.\.|-|to)\s*(-?\d+(?:\.\d+)?)$/);
    if(m){var a=parseFloat(m[1]),b=parseFloat(m[2]); if(a>b){var t=a;a=b;b=t;} return val>=a&&val<=b;}
  }
  return String(txt==null?'':txt).toLowerCase().indexOf(q.toLowerCase())>=0;
}
function colFilterVal(s,k){   // -> [numericValueOrNull, displayText] for column key k
  if(k=='s')return [null,s.s];
  if(k=='v')return [null,s.v];
  if(k=='lineage')return [null,s.lineage||'NA'];
  var v=s.m[k]; return [v,(v==null?'NA':fmt(v,(MET[k]||{}).kind))];
}
function colMatch(s){if(!st.showColF)return true;   // filters apply only while the filter row is shown
  for(var k in st.colf){var raw=st.colf[k]; if(!raw||!raw.trim())continue;
  var pv=colFilterVal(s,k); if(!colMatchOne(raw,pv[0],pv[1]))return false;} return true;}
function colAnyActive(){if(!st.showColF)return false; for(var k in st.colf){if(st.colf[k]&&st.colf[k].trim())return true;} return false;}
// any row filter active (used to surface a "clear filters" affordance so users never lose track of why rows vanished)
function anyFilterActive(){return !!(st.q||st.onlyFlagged||st.flagFilter||st.linFilter||st.ancOnly||colAnyActive());}
function clearAllFilters(){
  st.q=''; st.onlyFlagged=false; st.flagFilter=null; st.linFilter=null; st.ancOnly=null; st.colf={};
  var q=el('q'); if(q)q.value=''; var of=el('of'); if(of)of.checked=false;
  Array.prototype.forEach.call(document.querySelectorAll('#ancseg button'),function(x){x.classList.toggle('on',(x.getAttribute('data-a')||'')=='');});
  renderAll();
}
function dotColor(s){return st.colorBy=='lineage'?linColor(s.lineage):VCOL[s.v];}
// shared colour key for every dot plot (adapts to the QC/lineage colour toggle)
function colorLegend(){
  var sw='display:inline-block;width:11px;height:11px;border-radius:3px;margin-right:5px;vertical-align:-1px;background:';
  var items=(st.colorBy=='lineage'
    ? (R.lineages||[]).map(function(l){return '<span><i style="'+sw+linColor(l)+'"></i>'+esc(l)+'</span>';}).join('')||'<span class="c">no lineage assigned</span>'
    : '<span><i style="'+sw+VCOL.PASS+'"></i>PASS</span><span><i style="'+sw+VCOL.WARN+'"></i>WARN</span><span><i style="'+sw+VCOL.FAIL+'"></i>FAIL</span>');
  return '<div class="legend" style="justify-content:center">'+items+'</div>';
}
// genome landscape: missing-fraction (0 callable -> 1 missing) mapped to a pale->red heat colour
function heatCol(mv){if(mv==null)return TH.cellnull;var a=isDark()?[34,46,60]:[238,244,240],b=[214,64,58];
  return 'rgb('+Math.round(a[0]+(b[0]-a[0])*mv)+','+Math.round(a[1]+(b[1]-a[1])*mv)+','+Math.round(a[2]+(b[2]-a[2])*mv)+')';}
function fmtpos(p){return p>=1e6?(p/1e6).toFixed(2)+' Mb':p>=1e3?Math.round(p/1e3)+' kb':(''+p)+' bp';}
function tip(h,x,y){var t=el('tt'); if(!h){t.style.opacity=0;return;} t.innerHTML=h;
  var w=window.innerWidth; t.style.left=Math.min(x+13,w-t.offsetWidth-10)+'px'; t.style.top=(y+13)+'px'; t.style.opacity=1;}
function setHi(s){st.hi=(st.hi==s?null:s); renderAll();}

function donut(c){var t=(c.PASS+c.WARN+c.FAIL)||1,R0=38,C=2*Math.PI*R0,off=0,segs='';
  segs+='<circle cx="46" cy="46" r="'+R0+'" fill="none" stroke="'+TH.track+'" stroke-width="13"/>';
  [['PASS',VFILL.PASS],['WARN',VFILL.WARN],['FAIL',VFILL.FAIL]].forEach(function(p){var frac=c[p[0]]/t,len=frac*C;
    if(len<=0)return;
    segs+='<circle cx="46" cy="46" r="'+R0+'" fill="none" stroke="'+p[1]+'" stroke-width="13" stroke-dasharray="'+len.toFixed(2)+' '+(C-len).toFixed(2)+'" stroke-dashoffset="'+(-off).toFixed(2)+'" transform="rotate(-90 46 46)"/>'; off+=len;});
  var pct=Math.round(100*c.PASS/t);
  return '<svg width="92" height="92" viewBox="0 0 92 92">'+segs+'<text x="46" y="42" text-anchor="middle" font-size="22" font-weight="700" fill="'+TH.ink+'" letter-spacing="-.5">'+pct+'%</text><text x="46" y="58" text-anchor="middle" font-size="9.5" fill="'+TH.axis+'" letter-spacing=".1em">PASS</text></svg>';}

// ---- Executive summary: cohort KPIs, quality profile and headline findings (for the PI receiving the file)
function _median(vals){var a=vals.filter(function(v){return v!=null;}).sort(function(x,y){return x-y;}); if(!a.length)return null; var m=Math.floor(a.length/2); return a.length%2?a[m]:(a[m-1]+a[m])/2;}
function _range(vals){var a=vals.filter(function(v){return v!=null;}); return a.length?[Math.min.apply(null,a),Math.max.apply(null,a)]:null;}
function renderExec(){
  var host=el('exec_body'); if(!host)return;
  var S=R.samples, N=S.length, c=R.counts, toEx=c.FAIL;
  function med(k){return _median(S.map(function(s){return s.m[k];}));}
  function rg(k){return _range(S.map(function(s){return s.m[k];}));}
  var mDepth=med('mean_depth'), mBreadth=med('breadth_pct'), rDepth=rg('mean_depth'), rBreadth=rg('breadth_pct');
  var contam=(R.kraken&&R.kraken.samples)?R.kraken.samples.filter(function(k){return k.primary&&k.primary.pct<90;}).length:null;
  var resSamp=null, resDrugs=[];
  if(R.dr&&R.dr.calls){var rs={},dd={}; R.dr.calls.forEach(function(cl){if(cl.gn===1||cl.gn===2){rs[cl.s]=1; if(cl.drug)dd[cl.drug]=(dd[cl.drug]||0)+1;}});
    resSamp=Object.keys(rs).length; resDrugs=Object.keys(dd).sort(function(a,b){return dd[b]-dd[a];}).slice(0,4);}
  function card(l,n,s,tone){return '<div class="kpi'+(tone?' '+tone:'')+'"><div class="kpi-l">'+l+'</div><div class="kpi-n">'+n+'</div><div class="kpi-s">'+(s||'')+'</div></div>';}
  var passRate=N?Math.round(c.PASS/N*100):0;
  var cards=[
    card('Samples', N, '<span class="v PASS xs">'+c.PASS+' pass</span> <span class="v WARN xs">'+c.WARN+' warn</span> <span class="v FAIL xs">'+c.FAIL+' fail</span>'),
    card('Pass rate', passRate+'%', toEx?('<b>'+toEx+'</b> to exclude'):'all usable', passRate>=80?'good':(passRate>=50?'warn':'bad')),
    card('Median depth', (mDepth!=null?fmt(mDepth,'float')+'&#215;':'NA'), rDepth?(fmt(rDepth[0],'float')+'-'+fmt(rDepth[1],'float')+'&#215; range'):''),
    card('Median breadth', (mBreadth!=null?mBreadth.toFixed(1)+'%':'NA'), rBreadth?(rBreadth[0].toFixed(0)+'-'+rBreadth[1].toFixed(0)+'% range'):'')];
  if(contam!=null) cards.push(card('Contamination', contam, contam?'sample(s) &lt; 90% primary':'none flagged', contam?'warn':'good'));
  if(resSamp!=null) cards.push(card('Resistance', resSamp, resDrugs.length?esc(resDrugs.join(' · ')):'no R calls', resSamp?'warn':''));
  var flagc={}; S.forEach(function(s){(s.f||[]).forEach(function(f){flagc[f]=(flagc[f]||0)+1;});});
  var flags=Object.keys(flagc).sort(function(a,b){return flagc[b]-flagc[a];});
  var prof=flags.length?flags.map(function(f){var n=flagc[f],w=Math.round(n/N*100),fatal=FAILF[f];
    return '<div class="qcp-row"><span class="qcp-lab">'+f+'</span><span class="qcp-bar"><span style="width:'+Math.max(4,w)+'%;background:'+(fatal?'var(--fail)':'var(--warn)')+'"></span></span><span class="qcp-n">'+n+'</span></div>';}).join(''):'<div class="krk-mut">No sample trips any check at the current thresholds.</div>';
  var linc={}; S.forEach(function(s){if(s.lineage)linc[s.lineage]=(linc[s.lineage]||0)+1;});
  var lins=Object.keys(linc).sort();
  var strip=lins.length?'<div class="lincomp-bar" style="margin:0 0 8px">'+lins.map(function(l){return '<div class="lseg" style="width:'+(linc[l]/N*100).toFixed(2)+'%;background:'+linColor(l)+'" title="'+esc(l)+' n='+linc[l]+'"></div>';}).join('')+'</div><div class="lincomp-lab" style="padding:0">'+lins.map(function(l){return '<span class="lchip"><i style="background:'+linColor(l)+'"></i>'+esc(l)+' <b>'+linc[l]+'</b></span>';}).join('')+'</div>':'<div class="krk-mut">no lineage calls</div>';
  var narr='<b>'+N+'</b> samples analysed against '+(R.provenance&&R.provenance.reference?esc(R.provenance.reference):'the reference')+' &#183; <b>'+c.PASS+'</b> pass, '+(toEx?'<b class="tone-bad">'+toEx+'</b> recommended for exclusion':'0 to exclude')+' &#183; median depth <b>'+(mDepth!=null?fmt(mDepth,'float')+'&#215;':'NA')+'</b>, breadth <b>'+(mBreadth!=null?mBreadth.toFixed(1)+'%':'NA')+'</b>'+(contam?' &#183; <b class="tone-warn">'+contam+'</b> possibly contaminated':'')+(resSamp?' &#183; drug resistance in <b class="tone-warn">'+resSamp+'</b> sample(s)':'')+'.';
  host.innerHTML='<div class="exec-narr">'+narr+'</div><div class="exec-grid">'+cards.join('')+'</div>'+
    '<div class="exec-cols"><div class="exec-block"><div class="exec-h">QC quality profile <span class="krk-mut">- samples tripping each check (red = gate-failing)</span></div>'+prof+'</div>'+
    '<div class="exec-block"><div class="exec-h">Cohort lineages</div>'+strip+'</div></div>';
}
function renderOverview(){
  var c=R.counts;
  el('summary').innerHTML=donut(c)+'<div class="counts">'+
    '<div class="c all"><div class="n">'+R.samples.length+'</div><div class="l">samples</div></div>'+
    ['PASS','WARN','FAIL'].map(function(v){return '<div class="c '+v.toLowerCase()+'"><div class="n">'+c[v]+'</div><div class="l">'+v.toLowerCase()+'</div></div>';}).join('')+'</div>';
  var freq={}; R.samples.forEach(function(s){s.f.forEach(function(f){freq[f]=(freq[f]||0)+1;});});
  var keys=Object.keys(freq).sort(function(a,b){return freq[b]-freq[a];});
  el('chips').innerHTML='<span class="t">flags</span>'+(keys.length?keys.map(function(f){return '<span class="chip'+(st.flagFilter==f?' on':'')+'" data-f="'+f+'" role="button" tabindex="0" aria-pressed="'+(st.flagFilter==f?'true':'false')+'" aria-label="filter by '+f+'">'+f+'<span class="k">'+freq[f]+'</span></span>';}).join('')+'<span class="chip-hint">click to filter</span>':'<span style="color:#94a3b8;font-size:12px">none - every sample clear ✓</span>');
  Array.prototype.forEach.call(document.querySelectorAll('#chips .chip'),function(ch){function tog(){var f=ch.getAttribute('data-f');st.flagFilter=(st.flagFilter==f?null:f);st.onlyFlagged=false;renderAll();}
    ch.onclick=tog; ch.onkeydown=function(e){if(e.key=='Enter'||e.key==' '||e.key=='Spacebar'){e.preventDefault();tog();}};});
}

function renderTable(){
  var mets=R.metrics.filter(function(m){return !st.hidden[m.key];});
  function hsa(k){return ' tabindex="0" aria-sort="'+(st.sortKey==k?(st.asc?'ascending':'descending'):'none')+'"';}  // sortable-header a11y
  function sarr(k){return st.sortKey==k?(st.asc?icon('chevronUp','sort'):icon('chevronDown','sort')):'';}  // active-sort direction caret
  var head='<tr><th class="s" data-k="s"'+hsa('s')+'><input type="checkbox" id="cbAll" title="exclude all shown samples"><span class="hlab"> Sample</span>'+sarr('s')+'</th><th data-k="v"'+hsa('v')+'>QC'+sarr('v')+'</th>'+
    mets.map(function(m){var d=(R.defs[m.key]||[''])[0];return '<th data-k="'+m.key+'"'+hsa(m.key)+' title="'+esc(d)+'">'+esc(m.label)+(d?'<span class="infoi" title="'+esc(d)+'">i</span>':'')+sarr(m.key)+'</th>';}).join('')+
    '<th data-k="lineage"'+hsa('lineage')+' style="text-align:left">Lineage</th></tr>';
  // optional per-column filter row: one input per column (numeric ops on metric columns, substring otherwise)
  function cfIn(k,ph,lab,num){return '<input class="cfx" type="search" data-fk="'+k+'" value="'+esc(st.colf[k]||'')+'" placeholder="'+esc(ph)+'" aria-label="Filter '+esc(lab||k)+'"'+(num?' title="operators: &gt; &gt;= &lt; &lt;= = , a range 5-9 or 5..9; otherwise matches the text"':'')+'>';}
  var filtRow = st.showColF ? ('<tr class="colfilt"><th class="s">'+cfIn('s','name…','Sample')+'</th><th>'+cfIn('v','PASS/WARN…','QC status')+'</th>'+
    mets.map(function(m){return '<th>'+cfIn(m.key,'>50  5-9…',m.label,1)+'</th>';}).join('')+'<th>'+cfIn('lineage','L4…','Lineage')+'</th></tr>') : '';
  var linRank={}; (R.lineages||[]).forEach(function(l,i){linRank[l]=i;});
  function lr(s){return (s.lineage&&linRank[s.lineage]!=null)?linRank[s.lineage]:9999;}
  var shown=visible().filter(colMatch);
  var rows=shown.slice().sort(function(a,b){
    if(st.groupLin){var la=lr(a),lb=lr(b);if(la!=lb)return la-lb;}
    var k=st.sortKey,x=(k=='s')?a.s:(k=='v'?qcScore(a):a.m[k]),y=(k=='s')?b.s:(k=='v'?qcScore(b):b.m[k]),c;   // QC column sorts by severity, not the verdict string
    if(typeof x=='number'&&typeof y=='number')c=x-y;else c=String(x==null?'':x).localeCompare(String(y==null?'':y));return st.asc?c:-c;});
  var lastLin=null, ncol=mets.length+3;
  var TBL_CAP=400, total=rows.length, capped=total>TBL_CAP, draw=capped?rows.slice(0,TBL_CAP):rows;
  var body=draw.map(function(s){
    var pre='';
    if(st.groupLin){var lk=s.lineage||'NA'; if(lk!==lastLin){lastLin=lk;
      pre='<tr class="lingrp"><td class="s" colspan="'+ncol+'" style="text-align:left"><span class="ldot" style="background:'+linColor(s.lineage)+'"></span>'+esc(lk)+'</td></tr>';}}
    var badge=s.anc?'<span class="abadge" title="ancient (aDNA) sample">aDNA</span>':'';
    var tds='<td class="s" data-s="'+esc(s.s)+'"><input type="checkbox" class="cbx" data-s="'+esc(s.s)+'"'+(st.excl[s.s]?' checked':'')+' aria-label="basket '+esc(s.s)+'"><span class="sname" data-s="'+esc(s.s)+'" role="button" tabindex="0" aria-label="Open profile for '+esc(s.s)+'">'+esc(s.s)+'</span>'+badge+'</td><td data-v="'+s.v+'"><span class="v '+s.v+'">'+s.v+'</span></td>';
    mets.forEach(function(m){var v=s.m[m.key];
      if(v==null){tds+='<td class="na" data-v="">NA</td>';return;}
      var r=RANGES[m.key],nn=r[1]>r[0]?(v-r[0])/(r[1]-r[0]):0;nn=Math.max(0,Math.min(1,nn));var p=(nn*100).toFixed(1);
      tds+='<td data-v="'+v+'" style="background:linear-gradient(90deg,'+BAR[m.dir]+'2b 0 '+p+'%,#0000 '+p+'%)">'+fmt(v,m.kind)+'</td>';});
    tds+='<td data-v="'+esc(s.lineage||'')+'" style="text-align:left">'+(s.lineage?'<span class="ldot" style="background:'+linColor(s.lineage)+'"></span>':'')+esc(s.lineage||'NA')+'</td>';
    return pre+'<tr class="'+(st.hi==s.s?'hl':'')+'" data-s="'+esc(s.s)+'">'+tds+'</tr>';}).join('');
  var bodyOut=draw.length?body:'<tr><td colspan="'+ncol+'" style="text-align:left;color:#5f6f81;padding:14px 12px">No samples match the current filters.</td></tr>';
  var t=el('gstable'); t.innerHTML='<thead>'+head+filtRow+'</thead><tbody>'+bodyOut+'</tbody>';
  var af=anyFilterActive(), cntTxt=capped?('first '+TBL_CAP+' of '+total):(total+' / '+R.samples.length);
  el('nshown').innerHTML=cntTxt+' shown'+(af?' <a href="#" id="clrfilt" style="color:var(--accent);cursor:pointer;margin-left:7px;text-decoration:none">clear filters '+icon('x','sort')+'</a>':'');
  var cf=el('clrfilt'); if(cf)cf.onclick=function(e){e.preventDefault();clearAllFilters();};
  Array.prototype.forEach.call(t.querySelectorAll('th[data-k]'),function(th){function srt(){var k=th.getAttribute('data-k');if(st.sortKey==k)st.asc=!st.asc;else{st.sortKey=k;st.asc=(k=='s');}renderTable();saveState();}
    th.onclick=srt; th.onkeydown=function(e){if(e.key=='Enter'||e.key==' '||e.key=='Spacebar'){e.preventDefault();srt();}};});
  Array.prototype.forEach.call(t.querySelectorAll('.cfx'),function(inp){
    inp.onclick=function(e){e.stopPropagation();};
    inp.oninput=function(){var k=inp.getAttribute('data-fk'),pos=inp.selectionStart;st.colf[k]=inp.value;clearTimeout(_cfdb);_cfdb=setTimeout(function(){renderTable();
      var again=el('gstable').querySelector('.cfx[data-fk="'+k+'"]');if(again){again.focus();try{again.setSelectionRange(pos,pos);}catch(e){}}},140);};});
  Array.prototype.forEach.call(t.querySelectorAll('tbody tr[data-s]'),function(tr){tr.onclick=function(){setHi(tr.getAttribute('data-s'));};});
  // curation basket: exclusion checkboxes (stopPropagation so they don't sort/highlight)
  var cbAll=el('cbAll');
  if(cbAll){cbAll.checked=shown.length>0&&shown.every(function(s){return st.excl[s.s];});
    cbAll.onclick=function(e){e.stopPropagation();var on=cbAll.checked;shown.forEach(function(s){if(on)st.excl[s.s]=1;else delete st.excl[s.s];});renderTable();renderCuration();};}
  Array.prototype.forEach.call(t.querySelectorAll('.cbx'),function(cb){
    cb.onclick=function(e){e.stopPropagation();};
    cb.onchange=function(){var s=cb.getAttribute('data-s');if(cb.checked)st.excl[s]=1;else delete st.excl[s];renderCuration();
      var cba=el('cbAll');if(cba){cba.checked=shown.length>0&&shown.every(function(x){return st.excl[x.s];});}};});
  Array.prototype.forEach.call(t.querySelectorAll('.sname'),function(sp){sp.onclick=function(e){e.stopPropagation();openDetail(sp.getAttribute('data-s'));};
    sp.onkeydown=function(e){if(e.key=='Enter'||e.key==' '||e.key=='Spacebar'){e.preventDefault();e.stopPropagation();openDetail(sp.getAttribute('data-s'));}};});
}

function renderPlots(){
  var host=el('plots'); host.innerHTML='';
  var W=host.clientWidth||900, narrow=W<520, labelW=narrow?96:156, svgW=Math.max(narrow?190:240,W-labelW-4), padL=8, rightPad=64, pw=svgW-padL-rightPad, H=30, cy=H/2;
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
        stk=(s.v=='FAIL'||big)?' stroke="'+TH.ink+'" stroke-width="'+(big?1.3:0.6)+'"':'';
        return '<circle cx="'+x.toFixed(1)+'" cy="'+y.toFixed(1)+'" r="'+rr+'" fill="'+dotColor(s)+'" opacity="'+op+'"'+stk+' data-s="'+esc(s.s)+'" data-lin="'+esc(s.lineage||'')+'" data-pk="'+esc(pk)+'" data-val="'+s.m[pk]+'" data-lab="'+esc(mt.label)+'" data-kind="'+mt.kind+'"/>';}).join('');
      var mx=padL+(MED[pk]-lo)/(hi-lo)*pw;
      inner=bandSVG+'<line x1="'+padL+'" y1="'+cy+'" x2="'+(padL+pw)+'" y2="'+cy+'" stroke="'+TH.grid+'"/><line x1="'+mx.toFixed(1)+'" y1="4" x2="'+mx.toFixed(1)+'" y2="'+(H-4)+'" stroke="#64748b" stroke-dasharray="2 2"/>'+inner;
    }else if(st.ptype=='bar'){
      var sr=rows.slice().sort(function(a,b){return b.m[pk]-a.m[pk];}); var bw=pw/sr.length;
      inner=sr.map(function(s,i){var h=(s.m[pk]-Math.min(lo,0))/(hi-Math.min(lo,0))*(H-6),x=padL+i*bw,dim=(st.q||st.onlyFlagged||st.flagFilter||st.ancOnly||st.linFilter)&&!vis[s.s],big=(st.hi==s.s);
        return '<rect class="hit" x="'+x.toFixed(1)+'" y="'+(H-3-h).toFixed(1)+'" width="'+Math.max(bw-0.5,0.6).toFixed(1)+'" height="'+Math.max(h,0.5).toFixed(1)+'" fill="'+(big?''+TH.ink+'':dotColor(s))+'" opacity="'+(dim?0.12:(s.v!='PASS'?0.95:0.62))+'" data-s="'+esc(s.s)+'" data-lin="'+esc(s.lineage||'')+'" data-pk="'+esc(pk)+'" data-val="'+s.m[pk]+'" data-lab="'+esc(mt.label)+'" data-kind="'+mt.kind+'"/>';}).join('');
    }else{ // histogram
      var nb=Math.min(30,Math.max(8,Math.round(Math.sqrt(rows.length)))),cnt=zeros(nb);
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
  var host=el('scatter'); var cs=getComputedStyle(host); var avail=(host.clientWidth||560)-(parseFloat(cs.paddingLeft)||0)-(parseFloat(cs.paddingRight)||0); if(!(avail>0))avail=520; var cap=host.classList.contains('expanded')?880:700; var S=Math.max(240,Math.min(cap,avail)); var pad=42, plot=S-pad-14, H=(host.classList.contains('expanded')?Math.min(660,S):340), ph=H-pad-14;
  var xk=st.sx,yk=st.sy,xm=MET[xk],ym=MET[yk];
  var rows=R.samples.filter(function(s){return s.m[xk]!=null&&s.m[yk]!=null;});
  if(!rows.length){host.innerHTML='<div class="pad nd">no data for these axes</div>';return;}
  var xr=RANGES[xk],yr=RANGES[yk];
  function sx(v){return pad+(xr[1]>xr[0]?(v-xr[0])/(xr[1]-xr[0]):0.5)*plot;}
  function sy(v){return H-pad-(yr[1]>yr[0]?(v-yr[0])/(yr[1]-yr[0]):0.5)*ph;}
  var vis={}; visible().forEach(function(s){vis[s.s]=1;});
  var dots=rows.map(function(s){var big=(st.hi==s.s),dim=(st.q||st.onlyFlagged||st.flagFilter||st.ancOnly||st.linFilter)&&!vis[s.s];
    return '<circle cx="'+sx(s.m[xk]).toFixed(1)+'" cy="'+sy(s.m[yk]).toFixed(1)+'" r="'+(big?5.4:3.4)+'" fill="'+dotColor(s)+'" opacity="'+(dim?0.12:0.82)+'"'+((s.v=='FAIL'||big)?' stroke="'+TH.ink+'" stroke-width="'+(big?1.4:0.6)+'"':'')+' data-s="'+esc(s.s)+'" data-lin="'+esc(s.lineage||'')+'" data-x="'+s.m[xk]+'" data-y="'+s.m[yk]+'" data-xl="'+esc(xm.label)+'" data-yl="'+esc(ym.label)+'" data-xk="'+xm.kind+'" data-yk="'+ym.kind+'"/>';}).join('');
  var ticks='';[0,0.5,1].forEach(function(t){var gx=pad+t*plot,gy=H-pad-t*ph;
    ticks+='<line x1="'+gx+'" y1="'+pad+'" x2="'+gx+'" y2="'+(H-pad)+'" stroke="'+TH.grid+'"/><line x1="'+pad+'" y1="'+gy+'" x2="'+(pad+plot)+'" y2="'+gy+'" stroke="'+TH.grid+'"/>'+
    '<text x="'+gx+'" y="'+(H-pad+13)+'" font-size="9" fill="#94a3b8" text-anchor="middle">'+shortv(xr[0]+t*(xr[1]-xr[0]),xm.kind)+'</text>'+
    '<text x="'+(pad-6)+'" y="'+(gy+3)+'" font-size="9" fill="#94a3b8" text-anchor="end">'+shortv(yr[0]+t*(yr[1]-yr[0]),ym.kind)+'</text>';});
  host.innerHTML='<svg width="'+S+'" height="'+H+'" id="scsvg" style="display:block;max-width:100%;margin:0 auto">'+
    '<line x1="'+pad+'" y1="'+(H-pad)+'" x2="'+(pad+plot)+'" y2="'+(H-pad)+'" stroke="'+TH.axis+'"/><line x1="'+pad+'" y1="'+pad+'" x2="'+pad+'" y2="'+(H-pad)+'" stroke="'+TH.axis+'"/>'+
    ticks+dots+
    '<text x="'+(pad+plot/2)+'" y="'+(H-6)+'" font-size="11" fill="'+TH.mut+'" text-anchor="middle">'+esc(xm.label)+'</text>'+
    '<text x="12" y="'+(pad+ph/2)+'" font-size="11" fill="'+TH.mut+'" text-anchor="middle" transform="rotate(-90 12 '+(pad+ph/2)+')">'+esc(ym.label)+'</text></svg>'+
    colorLegend();
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
function corrCol(r){if(r==null)return TH.cellnull;var a=Math.abs(r),base=r>=0?[224,84,79]:[79,131,194],w=isDark()?[30,42,56]:[247,249,252];
  return 'rgb('+w.map(function(c,i){return Math.round(c+(base[i]-c)*a);}).join(',')+')';}
function renderCorr(){
  var host=el('corr_body'); if(!host)return;
  var keys=DIST.filter(function(k){return MET[k];}).slice(0,16);
  var pool=visible();
  if(pool.length<4||keys.length<2){host.innerHTML='<div class="pad nd">not enough data for a correlation matrix (need >= 4 samples in view).</div>';return;}
  var vals={}; keys.forEach(function(k){vals[k]=pool.map(function(s){return s.m[k];});});
  var n=keys.length, cell=Math.max(16,Math.min(34,Math.floor((Math.min(host.clientWidth||560,host.classList.contains('expanded')?1100:940)-110)/n)));
  var padL=96,padT=8, W=padL+n*cell+8, H=padT+n*cell+128;
  var svg='<svg width="'+W+'" height="'+H+'" id="corrsvg" style="display:block">';   // no max-width: let #corr_body{overflow-x:auto} scroll instead of clipping the right columns
  keys.forEach(function(k,j){var cx=padL+j*cell+cell/2;
    svg+='<text x="'+cx+'" y="'+(padT+n*cell+12)+'" font-size="8.5" fill="'+TH.mut+'" text-anchor="end" transform="rotate(-55 '+cx+' '+(padT+n*cell+12)+')">'+esc(MET[k].label)+'</text>';});
  keys.forEach(function(k,i){var cy=padT+i*cell+cell/2;
    svg+='<text x="'+(padL-6)+'" y="'+(cy+3)+'" font-size="8.5" fill="'+TH.mut+'" text-anchor="end">'+esc(MET[k].label)+'</text>';});
  for(var i=0;i<n;i++)for(var j=0;j<n;j++){
    var cx=padL+j*cell,cy=padT+i*cell,r;
    if(i===j){r=vals[keys[i]].some(function(v){return v!=null;})?1:null;}   // grey out an all-NA metric
    else{var xs=[],ys=[]; for(var t=0;t<pool.length;t++){var xv=vals[keys[j]][t],yv=vals[keys[i]][t];
        if(xv!=null&&yv!=null){xs.push(xv);ys.push(yv);}} r=spearman(xs,ys);}
    svg+='<rect x="'+cx+'" y="'+cy+'" width="'+(cell-1)+'" height="'+(cell-1)+'" rx="2" fill="'+corrCol(r)+'"'+
      ' data-xk="'+keys[j]+'" data-yk="'+keys[i]+'" data-r="'+(r==null?'':r.toFixed(2))+'"'+(i!==j?' style="cursor:pointer"':'')+'/>';
    if(cell>=22&&r!=null)svg+='<text x="'+(cx+(cell-1)/2)+'" y="'+(cy+(cell-1)/2+3)+'" font-size="7.5" fill="'+(Math.abs(r)>0.55?'#fff':'#5b6b7e')+'" text-anchor="middle" pointer-events="none">'+(r>0?'':'-')+Math.abs(r).toFixed(1).replace('0.','.')+'</text>';
  }
  var ly=padT+n*cell+94;
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
    return '<div class="stack" data-s="'+esc(s.s)+'" style="opacity:'+(dim?0.25:1)+(st.hi==s.s?';background:'+TH.hl:'')+'">'+
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
    '<td class="flags"><div class="fchips">'+s.f.map(function(f){return '<span class="chip'+(FAILF[f]?' failc':'')+'" data-f="'+f+'" title="'+esc(flagWhy(s,f))+'" style="cursor:pointer">'+f+'</span>';}).join(' ')+'</div>'+
      '<div class="flagrsn">'+s.f.map(function(f){return '<span>'+esc(flagWhy(s,f))+'</span>';}).join('')+'</div></td></tr>';}).join('')
    :'<tr><td colspan="3" style="text-align:left;color:#16a34a;padding:10px">All samples pass at the current thresholds.</td></tr>';
  var t=el('flagtable');
  t.innerHTML='<thead><tr><th class="s">Sample (worst first)</th><th>QC</th><th style="text-align:left">Flags &amp; reason</th></tr></thead><tbody>'+body+'</tbody>';
  Array.prototype.forEach.call(t.querySelectorAll('tbody tr[data-s]'),function(tr){tr.onclick=function(e){
    if(e.target.classList.contains('chip')){var f=e.target.getAttribute('data-f');st.flagFilter=(st.flagFilter==f?null:f);renderAll();el('gstats').scrollIntoView();return;}
    setHi(tr.getAttribute('data-s'));el('gstats').scrollIntoView();};});
  el('nflag').textContent=fl.length;
  var bw=el('basketFlagged'); if(bw){bw.onclick=function(){fl.forEach(function(s){st.excl[s.s]=1;});renderTable();renderCuration();};}
}

// ---- curation basket + exclusion exports ----
function nExcl(){return R.samples.filter(function(s){return st.excl[s.s];}).length;}
function renderCuration(){var ex=nExcl(),keep=R.samples.length-ex;
  var nb=el('nbasket'); if(nb)nb.innerHTML=icon('basket')+ex+' basketed';   // always-visible toolbar mirror of the basket
  el('curation').innerHTML='<div class="cur-intro"><b>Exclusion basket</b> - the set of samples you are dropping from the downstream analysis (phylogeny / clock). <b>FAIL samples start pre-selected</b> (their boxes are ticked); tick or untick any box in the table (or use the buttons), then export the drop list (<b>exclusion.tsv</b>) or the survivors (<b>keep_list.txt</b>).</div>'+
    '<div class="cur-read"><b>'+ex+'</b> to exclude <span class="arw">→</span> <b>'+keep+'</b> kept for downstream</div>'+
   '<div class="cur-btns"><button class="btn" data-cur="fail">exclude FAILs</button><button class="btn" data-cur="flagged">exclude all flagged</button>'+
   '<button class="btn" data-cur="clear">clear</button><button class="btn" data-cur="invert">invert (shown)</button>'+
   '<button class="btn prim" data-cur="excl">'+icon('download')+'exclusion.tsv</button><button class="btn prim" data-cur="keep">'+icon('download')+'keep_list.txt</button></div>';
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
function dl(txt,name,type){var blob=new Blob([txt],{type:type}),a=document.createElement('a');a.href=URL.createObjectURL(blob);a.download=name;a.click();URL.revokeObjectURL(a.href);toast('Saved '+name);}
var _toastT;
function toast(msg){var e0=el('toast'); if(!e0)return; e0.textContent=msg; e0.className='show'; clearTimeout(_toastT); _toastT=setTimeout(function(){e0.className='';},2200);}
// keep a screen reader / keyboard user inside the open dialog: the page behind is made inert
function setBgInert(on){['toc','toc-toggle'].forEach(function(id){var e=el(id); if(e){if(on){e.setAttribute('inert','');e.setAttribute('aria-hidden','true');}else{e.removeAttribute('inert');e.removeAttribute('aria-hidden');}}});
  Array.prototype.forEach.call(document.querySelectorAll('header,.wrap'),function(e){if(on){e.setAttribute('inert','');e.setAttribute('aria-hidden','true');}else{e.removeAttribute('inert');e.removeAttribute('aria-hidden');}});}

// ---- per-sample detail modal ----
function pctRank(key,val){if(val==null)return null;var vs=R.samples.map(function(s){return s.m[key];}).filter(function(v){return v!=null;});if(!vs.length)return null;var b=0;vs.forEach(function(v){if(v<val)b++;});return Math.round(100*b/vs.length);}
function openDetail(sid){var s=null;R.samples.forEach(function(x){if(x.s==sid)s=x;});if(!s)return;st.detail=sid;
  var rows=R.metrics.map(function(m){var v=s.m[m.key],r=RANGES[m.key],nn=(v==null||r[1]<=r[0])?0:Math.max(0,Math.min(1,(v-r[0])/(r[1]-r[0]))),pr=pctRank(m.key,v);
    var _dd=(R.defs[m.key]||[''])[0];return '<div class="drow"><div class="dk" title="'+esc(_dd)+'">'+esc(m.label)+(_dd?'<span class="infoi" title="'+esc(_dd)+'">i</span>':'')+'</div>'+
      '<div class="dbarwrap"><div class="dbar" style="width:'+(nn*100).toFixed(1)+'%;background:'+BAR[m.dir]+'"></div></div>'+
      '<div class="dv">'+(v==null?'<span class="na">NA</span>':fmt(v,m.kind))+'</div><div class="dp">'+(pr==null?'':('p'+pr))+'</div></div>';}).join('');
  var lin='';if(s.linf){var ks=Object.keys(s.linf).sort(function(a,b){return s.linf[b]-s.linf[a];});
    lin='<div class="dsub">Lineage composition</div><div class="lcomp">'+ks.map(function(k){return '<div class="lrow"><span class="lk">'+esc(k)+'</span><div class="lbarw"><div class="lbar" style="width:'+(s.linf[k]*100).toFixed(1)+'%"></div></div><span class="lv">'+(s.linf[k]*100).toFixed(1)+'%</span></div>';}).join('')+'</div>';}
  var dmg='';if(s.anc){var d=s.dmg;
    if(d){dmg='<div class="dsub">aDNA damage (mapDamage2) - consistent with ancient DNA, not proof of authenticity</div>'+
      '<div class="drow"><div class="dk">5&#39; C&gt;T (pos 1)</div><div class="dbarwrap"><div class="dbar" style="width:'+Math.min(100,(d.ct1||0)*100/0.3).toFixed(1)+'%;background:'+(s.f.indexOf('DAMAGE_LOW')>=0?'#e0544f':'#22a06b')+'"></div></div><div class="dv">'+(d.ct1==null?'NA':(d.ct1*100).toFixed(1)+'%')+'</div><div class="dp"></div></div>'+
      '<div class="drow"><div class="dk">3&#39; G&gt;A (pos 1)</div><div class="dbarwrap"><div class="dbar" style="width:'+Math.min(100,(d.ga1||0)*100/0.3).toFixed(1)+'%;background:#4f83c2"></div></div><div class="dv">'+(d.ga1==null?'NA':(d.ga1*100).toFixed(1)+'%')+'</div><div class="dp"></div></div>'+
      '<div class="drow"><div class="dk">mean frag len</div><div class="dv" style="flex:1;text-align:left;color:var(--txt2)">'+(d.fraglen==null?'NA':d.fraglen.toFixed(0)+' bp')+'</div></div>';
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
  st._opener=document.activeElement; el('modal').classList.add('open'); setBgInert(true); var mx=el('modalx'); if(mx)mx.focus();}
function closeDetail(){st.detail=null;el('modal').classList.remove('open');setBgInert(false);
  if(st._opener&&st._opener.focus){try{st._opener.focus();}catch(e){}} st._opener=null;}

// ---- genome landscape: samples x reference-position heatmap, multi-track (missing / SNPs / het / indels) ----
var GTRACKS=[{k:'missing',lab:'Missing',base:[214,64,58]},{k:'snp',lab:'SNPs',base:[14,139,168]},
             {k:'het',lab:'Het',base:[221,138,26]},{k:'indel',lab:'Indels',base:[124,92,191]}];
function gtrackHas(k){return R.samples.some(function(s){return k=='missing'?!!s.miss:(s.trk&&s.trk[k]);});}
function gtrackMax(k){var mx=0;R.samples.forEach(function(s){var p=k=='missing'?s.miss:(s.trk&&s.trk[k]);if(p)p.forEach(function(v){if(v!=null&&v>mx)mx=v;});});return mx;}
function gcol(base,mv){if(mv==null)return TH.cellnull;if(mv>1)mv=1;var a=isDark()?[30,42,56]:[238,244,240];
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
    return '<div class="stack" data-s="'+esc(s.s)+'" style="opacity:'+(dim?0.25:1)+(st.hi==s.s?';background:'+TH.hl:'')+'">'+
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
  svg+='<line x1="'+gut+'" y1="'+profH+'" x2="'+(gut+plotW).toFixed(1)+'" y2="'+profH+'" stroke="'+TH.grid+'"/>';
  svg+='<path d="M '+gut.toFixed(1)+' '+profH+' L '+pts.join(' L ')+' L '+(gut+plotW).toFixed(1)+' '+profH+' Z" fill="'+bf+'"/>';
  svg+='<polyline points="'+pts.join(' ')+'" fill="none" stroke="'+bs+'" stroke-width="1.3" stroke-linejoin="round"/>';
  svg+='<text x="0" y="11" font-size="9" font-weight="600" fill="'+TH.mut+'">'+TLAB[tk]+'</text>';
  svg+='<text x="'+(gut+plotW).toFixed(1)+'" y="11" font-size="8.5" fill="#94a3b8" text-anchor="end">peak '+(tk=='missing'?(pmax*100).toFixed(0)+'%':(Math.round(pmax*10)/10)+'/bin')+'</text>';
  for(var pj=z0;pj<=z1;pj++)svg+='<rect x="'+x(pj).toFixed(1)+'" y="0" width="'+bw.toFixed(1)+'" height="'+profH+'" fill="transparent" data-bin="'+pj+'" data-v="'+(agg[pj]==null?'':(tk=='missing'?(agg[pj]*100).toFixed(0):(Math.round(agg[pj]*10)/10)))+'"/>';
  rows.forEach(function(s,r){var dim=(st.q||st.onlyFlagged||st.flagFilter||st.ancOnly||st.linFilter)&&!vis[s.s],yy=y0+r*rowH;
    svg+='<rect x="'+(gut-9)+'" y="'+yy+'" width="5" height="'+(rowH-0.4).toFixed(1)+'" fill="'+linColor(s.lineage)+'"'+(dim?' opacity="0.25"':'')+'/>';
    for(var i=z0;i<=z1;i++){var vv=raw(s,i),mv=vv==null?null:(tk=='missing'?vv:vv/mx);
      svg+='<rect x="'+x(i).toFixed(1)+'" y="'+yy+'" width="'+bw.toFixed(1)+'" height="'+(rowH-0.4).toFixed(1)+'" fill="'+gcol(base,mv)+'"'+(dim?' opacity="0.3"':'')+' data-s="'+esc(s.s)+'" data-bin="'+i+'" data-v="'+(vv==null?'':(tk=='missing'?(vv*100).toFixed(0):vv))+'"/>';}
  });
  var yb=y0+rows.length*rowH;
  if(mb)for(var mk=z0;mk<=z1;mk++){var mf=mb[mk]||0;if(mf<=0)continue;   // grey out the masked zones
    svg+='<rect x="'+x(mk).toFixed(1)+'" y="0" width="'+bw.toFixed(1)+'" height="'+yb+'" fill="'+TH.faint+'" opacity="'+(0.12+0.42*mf).toFixed(2)+'"/>';}
  [0,0.25,0.5,0.75,1].forEach(function(t){var px=gut+t*plotW,pos=Math.round((z0+t*winN)/nb*gl);
    svg+='<line x1="'+px.toFixed(1)+'" y1="'+yb+'" x2="'+px.toFixed(1)+'" y2="'+(yb+4)+'" stroke="'+TH.axis+'"/>'+
      '<text x="'+px.toFixed(1)+'" y="'+(yb+14)+'" font-size="8.5" fill="#94a3b8" text-anchor="'+(t==0?'start':t==1?'end':'middle')+'">'+fmtpos(pos)+'</text>';});
  if(tk=='snp'&&st.geneMark&&st.geneMark.b1>=z0&&st.geneMark.b0<=z1){var gm=st.geneMark,mx0=Math.max(gut,x(gm.b0)),mx1=Math.min(gut+plotW,x(gm.b1+1));
    svg+='<rect x="'+mx0.toFixed(1)+'" y="0" width="'+Math.max(2,mx1-mx0).toFixed(1)+'" height="'+yb+'" fill="none" stroke="'+TH.ink+'" stroke-width="1.2" stroke-dasharray="3 2"/>'+
      '<text x="'+Math.min(W-2,(mx0+mx1)/2).toFixed(1)+'" y="'+(profH+11)+'" font-size="9" font-weight="600" fill="'+TH.ink+'" text-anchor="middle">'+esc(gm.name||'')+'</text>';}
  svg+='<rect id="gbrush" x="0" y="0" width="0" height="'+yb+'" fill="rgba(14,139,168,.12)" stroke="#0e8ba8" stroke-width="1" stroke-dasharray="3 2" pointer-events="none" style="display:none"/>';
  GGEO={gut:gut,plotW:plotW,nb:nb,gl:gl,yb:yb,z0:z0,winN:winN};
  host.innerHTML=svg+'</svg>';
  // selection overlay only when it is a SUB-region of the current view (when zoomed to exactly the selection, the whole plot is it)
  if(st.gsel&&!(st.gsel.b0<=z0&&st.gsel.b1>=z1)){var gbx=el('gbrush');if(gbx){var bx0=Math.max(gut,x(st.gsel.b0)),bx1=Math.min(gut+plotW,x(st.gsel.b1+1));
    if(bx1>bx0){gbx.setAttribute('x',bx0.toFixed(1));gbx.setAttribute('width',(bx1-bx0).toFixed(1));gbx.style.display='';}else gbx.style.display='none';}}
}

// ---- SNP-dense gene / region detection (cohort SNP density along the reference) ----
function cohortSnp(){var nb=R.nbins||200,agg=zeros(nb),has=false,mb=(st.maskOn&&R.mask_bins)?R.mask_bins:null;
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
      '<td class="s" style="text-align:left">'+esc(x.name)+geneRvTag(x.name)+'</td><td style="text-align:left">'+fmtpos(x.start)+' - '+fmtpos(x.end)+'</td>'+
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
  svg+='<line x1="'+pl+'" y1="'+(pt+ih)+'" x2="'+(W-pr)+'" y2="'+(pt+ih)+'" stroke="'+TH.grid+'"/>';
  if(R.mask_bins)for(var mk=0;mk<nb&&mk<R.mask_bins.length;mk++){var mf=R.mask_bins[mk]||0;if(mf<0.5)continue;
    svg+='<rect x="'+x(mk).toFixed(1)+'" y="'+pt+'" width="'+(iw/nb+0.6).toFixed(1)+'" height="'+ih+'" fill="'+TH.faint+'" opacity="0.10"/>';}
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
    '<line x1="'+pl+'" y1="'+(pt+ih)+'" x2="'+(W-pr)+'" y2="'+(pt+ih)+'" stroke="'+TH.grid+'"/>'+tl+
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
        '<td class="s" style="text-align:left">'+esc(g.gene)+geneRvTag(g.gene)+'</td>'+
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
      return '<tr><td class="s" style="text-align:left">'+esc(g.gene)+geneRvTag(g.gene)+'</td>'+
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
      (low?'<div class="alow">'+icon('alert','sort')+'terminal C&gt;T below the '+(lim*100).toFixed(0)+'% authentication floor - possible modern contamination; verify before using as a calibration tip.</div>':'');
    return '<div class="acard'+(low?' low':'')+'">'+head+body+'</div>';}).join('');
  host.innerHTML='<div class="anote">Elevated terminal C&gt;T (5&#39;) / G&gt;A (3&#39;) deamination is <b>consistent with</b> post-mortem damage and is a necessary authentication signal. It does <b>not</b> by itself prove the DNA is ancient (deaminated modern contaminant DNA can mimic it) or endogenous, and its absence (below the floor) is a red flag for a modern sample mislabelled ancient. Treat this as a screen, not a proof.</div>'+
    '<div class="acards">'+cards+'</div>';
  Array.prototype.forEach.call(host.querySelectorAll('.sname'),function(sp){sp.onclick=function(){openDetail(sp.getAttribute('data-s'));};});
}

var DYNCOL={fixation:'#2f6fed',emergence:'#1f9d6b',loss:'#e6893a',nonsyn:'#d1495b',high_impact:'#7c3aed'};
var DYNHELP={emergence:'Emergence: the variant is (near-)absent at the first timepoint, then rises above the emergence threshold - a new allele appearing in this series.',fixation:'Fixation: the allele frequency reaches near 1.0 by the last timepoint - the variant has (almost) taken over.',loss:'Loss: the variant is present early then falls back toward 0 - an allele being lost from the series.',nonsyn:'Non-synonymous: the variant changes the protein (missense / stop / frameshift / splice / inframe indel), per snpEff - potentially functional.',high_impact:'High impact: snpEff predicts a HIGH-impact effect (frameshift, stop gained/lost...) - likely to disrupt the gene.'};
var dynState={sel:null,q:''};
var dynFilter={};
var dynZoom=250;   // trajectory-card width in px (zoom slider); smaller -> more charts per row
var dynShowDP=true;   // draw the per-timepoint read-depth (DP) bars behind each trajectory
function dynHasFlag(v){return v.flags&&v.flags.length;}
function dynColor(flags){ if(!flags)return '#9fb0c3'; if(flags.indexOf('fixation')>=0)return DYNCOL.fixation; if(flags.indexOf('emergence')>=0)return DYNCOL.emergence; if(flags.indexOf('loss')>=0)return DYNCOL.loss; if(flags.indexOf('high_impact')>=0)return DYNCOL.high_impact; return '#5b6b7e'; }
function dynMiniChart(v,th,showDP){
  var n=v.times.length,W=250,H=150,ml=30,mt=10,mb=26;
  var dps=v.dp||[], hasDP=showDP&&dps.some(function(d){return d!=null;});
  var maxDP=1; if(hasDP){ dps.forEach(function(d){ if(d!=null&&d>maxDP)maxDP=d; }); }
  var mr=hasDP?16:10, pw=W-ml-mr, ph=H-mt-mb;
  function X(i){ return ml+(n<=1?pw/2:(i/(n-1))*pw); }
  function Y(a){ return mt+(1-a)*ph; }                  // allele frequency (left axis)
  function YD(d){ return mt+ph-(d/maxDP)*ph; }          // read depth (right axis)
  var col=dynColor(v.flags), nonsyn=v.flags&&v.flags.indexOf('nonsyn')>=0;
  var svg='<svg viewBox="0 0 '+W+' '+H+'" width="100%" style="display:block"><title>Allele frequency (0-1, left axis, line) across timepoints'+(hasDP?'; read depth DP as bars with the value on top':'')+'. Hover for exact values.</title>';
  [0,0.5,1].forEach(function(a){ svg+='<line x1="'+ml+'" y1="'+Y(a).toFixed(1)+'" x2="'+(W-mr)+'" y2="'+Y(a).toFixed(1)+'" stroke="'+TH.grid+'"/><text x="'+(ml-5)+'" y="'+(Y(a)+3.5).toFixed(1)+'" text-anchor="end" font-size="10.5" fill="'+TH.mut+'">'+a.toFixed(1)+'</text>'; });
  if(hasDP){   // depth bars behind the AF line; each bar carries its DP value on top (see the pass after the line)
    var bw=Math.min(n<=1?18:(pw/n)*0.5, 16);
    dps.forEach(function(d,i){ if(d==null)return; var x=X(i), y=YD(d), h=(mt+ph)-y; svg+='<rect x="'+(x-bw/2).toFixed(1)+'" y="'+y.toFixed(1)+'" width="'+bw.toFixed(1)+'" height="'+Math.max(0,h).toFixed(1)+'" fill="#7ea8d6" opacity="0.45" rx="1.5"><title>t='+esc(v.times[i]==null?i:v.times[i])+'  DP='+d+'</title></rect>'; });
  }
  if(th){ [[th.emerge,DYNCOL.emergence],[th.fix,DYNCOL.fixation]].forEach(function(t){ svg+='<line x1="'+ml+'" y1="'+Y(t[0]).toFixed(1)+'" x2="'+(W-mr)+'" y2="'+Y(t[0]).toFixed(1)+'" stroke="'+t[1]+'" stroke-dasharray="3 3" opacity="0.3"/>'; }); }
  var ptsA=v.traj.map(function(a,i){return X(i).toFixed(1)+','+Y(a).toFixed(1);}), pts=ptsA.join(' '), y0=Y(0).toFixed(1), lastI=v.traj.length-1;
  svg+='<path d="M'+X(0).toFixed(1)+','+y0+' L'+ptsA.join(' L')+' L'+X(lastI).toFixed(1)+','+y0+' Z" fill="'+col+'" opacity="0.13"/>';   // soft area fill under the trajectory
  svg+='<polyline points="'+pts+'" fill="none" stroke="'+col+'" stroke-width="2.4" stroke-linecap="round" stroke-linejoin="round"/>';
  v.traj.forEach(function(a,i){ var last=(i===lastI); svg+='<circle cx="'+X(i).toFixed(1)+'" cy="'+Y(a).toFixed(1)+'" r="'+(last?4.4:3.4)+'" fill="'+col+'" stroke="'+(nonsyn?DYNCOL.nonsyn:TH.panel)+'" stroke-width="'+(nonsyn?1.8:(last?1.7:1))+'"><title>t='+esc(v.times[i]==null?i:v.times[i])+'  AF='+a.toFixed(3)+(hasDP&&dps[i]!=null?('  DP='+dps[i]):'')+'</title></circle>'; });
  if(hasDP){ dps.forEach(function(d,i){ if(d==null)return; svg+='<text x="'+X(i).toFixed(1)+'" y="'+(YD(d)-3).toFixed(1)+'" text-anchor="middle" font-size="8.5" font-weight="600" fill="'+TH.ink+'" stroke="'+TH.panel+'" stroke-width="2.6" paint-order="stroke" style="paint-order:stroke">'+d+'</text>'; }); }   // DP value on top of each bar
  v.times.forEach(function(t,i){ svg+='<text x="'+X(i).toFixed(1)+'" y="'+(H-8)+'" text-anchor="middle" font-size="10.5" fill="'+TH.mut+'">'+esc(t==null?i:t)+'</text>'; });
  svg+='</svg>';
  return svg;
}
function renderDynamics(){
  var host=el('dyn_body'), sec=el('dynamics'); if(!host)return;
  var D=R.dynamics;
  if(!(D&&D.groups&&D.groups.length)){ if(sec)sec.style.display='none'; return; }
  if(sec)sec.style.display='';
  var th=D.thresholds||{emerge:0.25,fix:0.9,loss:0.1};
  var vars=[];
  D.groups.forEach(function(g){ g.series.forEach(function(s){ vars.push({gene:s.gene||'(intergenic)',pos:s.pos,group:g.group,times:g.times,traj:s.traj,dp:s.dp,flags:s.flags,eff:s.eff,imp:s.imp,alt:s.alt,aa:s.aa,aa_h37rv:s.aa_h37rv}); }); });
  // per-series (group) metadata + the fields usable as a series filter: group-invariant (one value per
  // series, so timepoint/date -> the trajectory axis -> excluded) and with >1 value across series.
  var groupMeta={}; D.groups.forEach(function(g){ groupMeta[g.group]=g.meta||{}; });
  var mfields=(R.sample_meta&&R.sample_meta.fields)||[];
  var dynFields=[], dynFieldVals={};
  mfields.forEach(function(f){
    var invariant=true, all={};
    D.groups.forEach(function(g){ var vals=(g.meta&&g.meta[f])||[]; if(vals.length>1) invariant=false; vals.forEach(function(v){all[v]=1;}); });
    var distinct=Object.keys(all);
    if(invariant && distinct.length>1){ dynFields.push(f); dynFieldVals[f]=distinct.sort(); }
  });
  function passFilter(grp){ var m=groupMeta[grp]||{}; for(var f in dynFilter){ if(dynFilter[f]){ var vals=m[f]||[]; if(vals.indexOf(dynFilter[f])<0) return false; } } return true; }
  function filterActive(){ for(var f in dynFilter){ if(dynFilter[f]) return true; } return false; }
  var genes, geneList, singleGroup=false, visGroupName='';
  function recompute(){
    genes={}; var gset={};
    vars.forEach(function(v){ if(passFilter(v.group)){ (genes[v.gene]=genes[v.gene]||[]).push(v); gset[v.group]=1; } });
    var gkeys=Object.keys(gset);
    singleGroup=gkeys.length<=1;   // one visible series -> the per-card group tag is redundant
    visGroupName=gkeys.length===1?gkeys[0]:'';
    geneList=Object.keys(genes).sort(function(a,b){
      var fa=genes[a].filter(dynHasFlag).length, fb=genes[b].filter(dynHasFlag).length;
      return fb-fa || genes[b].length-genes[a].length || a.localeCompare(b);
    });
  }
  recompute();
  if(dynState.sel===null){
    dynState.sel={};
    var fg=geneList.filter(function(g){return genes[g].some(dynHasFlag);});
    (fg.length?fg:geneList.slice(0,6)).forEach(function(g){dynState.sel[g]=1;});
  }
  var filterUI=dynFields.length?('<div class="dyn-filters"><span class="snpmx-flabel" title="Show only the connected series (patient / passage line...) matching these metadata values. The time axis of each trajectory is unchanged.">filter series:</span>'+
    dynFields.map(function(f){ return '<label class="snpmx-fsel">'+esc(f)+' <select data-df="'+esc(f)+'"><option value="">all</option>'+
      dynFieldVals[f].map(function(v){return '<option value="'+esc(v)+'"'+(dynFilter[f]===v?' selected':'')+'>'+esc(v)+'</option>';}).join('')+'</select></label>'; }).join('')+
    '<button class="dyn-btn" id="dynfclear">clear</button></div>'):'';
  host.innerHTML=
    '<div class="dyn-controls">'+
      '<input id="dynsearch" class="dyn-search" type="search" title="Type a gene name to filter the gene chips and the grid below" placeholder="search gene..." value="'+esc(dynState.q)+'">'+
      '<button class="dyn-btn" id="dynFlag" title="Show only genes that have at least one flagged variant">flagged genes</button>'+
      '<button class="dyn-btn" id="dynAll" title="Select every gene that has a moving variant">all</button>'+
      '<button class="dyn-btn" id="dynNone" title="Deselect all genes">clear</button>'+
      '<label class="dyn-toggle" title="Show a per-timepoint read-depth (DP) bar behind each trajectory"><input type="checkbox" id="dynDP"'+(dynShowDP?' checked':'')+'> depth bars</label>'+
      '<label class="dyn-zoom" title="Resize the trajectory cards - drag left to fit more charts per row"><span>'+icon('search','sort')+'&#8211;/+</span><input type="range" id="dynzoom" min="165" max="360" step="5" value="'+dynZoom+'"></label>'+
      '<span class="dyn-count" id="dynCount"></span></div>'+
    filterUI+
    '<div class="dyn-legend"><span title="'+DYNHELP.emergence+'"><i style="background:'+DYNCOL.emergence+'"></i>emergence <span class="infoi">i</span></span><span title="'+DYNHELP.fixation+'"><i style="background:'+DYNCOL.fixation+'"></i>fixation <span class="infoi">i</span></span><span title="'+DYNHELP.loss+'"><i style="background:'+DYNCOL.loss+'"></i>loss <span class="infoi">i</span></span><span title="'+DYNHELP.nonsyn+'"><i style="border:2px solid '+DYNCOL.nonsyn+';background:var(--panel)"></i>non-synonymous <span class="infoi">i</span></span></div>'+
    '<div class="dyn-genechips" id="dynchips"></div>'+
    '<div class="dyn-grid" id="dyngrid"></div>';
  function paintChips(){
    var q=dynState.q.toLowerCase();
    var list=geneList.filter(function(g){return !q||g.toLowerCase().indexOf(q)>=0;});
    el('dynchips').innerHTML=list.length?list.map(function(g){
      var nf=genes[g].filter(dynHasFlag).length;
      return '<button class="dyn-chip'+(dynState.sel[g]?' sel':'')+'" data-g="'+esc(g)+'" title="Click to toggle. '+genes[g].length+' variant trajectory(ies)'+(nf?(', '+nf+' flagged'):'')+'">'+esc(g)+' <b>'+genes[g].length+'</b>'+(nf?'<i class="dyn-dot" title="'+nf+' flagged variant(s) in this gene"></i>':'')+'</button>';
    }).join(''):'<span class="c">'+(filterActive()?'no gene matches the current series filter':'no gene matches "'+esc(dynState.q)+'"')+'</span>';
    Array.prototype.forEach.call(el('dynchips').querySelectorAll('.dyn-chip'),function(b){ b.onclick=function(){ var g=b.getAttribute('data-g'); if(dynState.sel[g])delete dynState.sel[g]; else dynState.sel[g]=1; paintChips(); paintGrid(); }; });
  }
  function paintGrid(){
    var q=dynState.q.toLowerCase();
    var sel=geneList.filter(function(g){return dynState.sel[g] && (!q||g.toLowerCase().indexOf(q)>=0);});
    var nvar=0; sel.forEach(function(g){nvar+=genes[g].length;});
    el('dynCount').innerHTML=sel.length+' of '+geneList.length+' genes'+(sel.length?(' &#183; '+nvar+' trajectories'):'')+(singleGroup&&visGroupName?(' &#183; series <b>'+esc(visGroupName)+'</b>'):'')+(filterActive()?' (filtered)':'');
    var grid=el('dyngrid');
    grid.style.setProperty('--dyncw', dynZoom+'px');
    if(!sel.length){ grid.innerHTML='<div class="dyn-empty" style="grid-column:1/-1">&#128204; '+(geneList.length?(dynState.q?('no selected gene matches &quot;'+esc(dynState.q)+'&quot;'):'Search and select one or more genes above to see the allele-frequency trajectories of their variants.'):'no variant trajectory matches the current series filter.')+'</div>'; return; }
    var cards=[];
    sel.forEach(function(g){
      genes[g].slice().sort(function(a,b){ return ((dynHasFlag(b)?1:0)-(dynHasFlag(a)?1:0)) || (String(a.pos)>String(b.pos)?1:-1); }).forEach(function(v){
        var flagged=dynHasFlag(v);
        var chips=(v.flags||[]).map(function(f){return '<span class="dyn-fchip" style="background:'+(DYNCOL[f]||'#8895a6')+'" title="'+(DYNHELP[f]||f)+'">'+f+'</span>';}).join('');
        var posNum=String(v.pos).split(':').pop();
        cards.push('<div class="dyn-card'+(flagged?' flagged':'')+'">'+
          '<div class="dyn-card-h"><span class="dyn-cardgene" title="Gene (click its chip above to toggle)">'+esc(v.gene||'(intergenic)')+'</span>'+geneRvTag(v.gene)+(singleGroup?'':'<span class="dyn-grp" title="Connected series this variant belongs to (the metadata group column, e.g. patient / passage line)">'+esc(v.group)+'</span>')+'<span class="dyn-pos" title="Genomic position (contig:position) of this SNP: '+esc(v.pos)+'">'+esc(posNum)+'</span></div>'+
          '<div class="dyn-eff" title="Predicted effect (snpEff) and protein change HGVS.p: ref amino acid, codon position, alt amino acid">'+esc(v.eff||'variant')+(v.aa?(' &#183; <b class="dyn-aa">'+aaDual(v.aa,v.aa_h37rv)+'</b>'):(v.alt?(' &#183; &#8594;'+esc(v.alt)):''))+'</div>'+
          dynMiniChart(v,th,dynShowDP)+
          '<div class="dyn-card-f">'+(chips||'<span class="c" title="no emergence / fixation / loss / non-synonymous event for this variant">no event</span>')+'<span class="dyn-traj" title="Allele frequency at each timepoint, in chronological order">'+v.traj.map(function(a){return a.toFixed(2);}).join(' &#8594; ')+'</span></div>'+
        '</div>');
      });
    });
    grid.innerHTML=cards.join('');
  }
  el('dynsearch').oninput=function(){ dynState.q=this.value; paintChips(); paintGrid(); };
  el('dynFlag').onclick=function(){ dynState.sel={}; geneList.filter(function(g){return genes[g].some(dynHasFlag);}).forEach(function(g){dynState.sel[g]=1;}); paintChips(); paintGrid(); };
  el('dynAll').onclick=function(){ dynState.sel={}; geneList.forEach(function(g){dynState.sel[g]=1;}); paintChips(); paintGrid(); };
  el('dynNone').onclick=function(){ dynState.sel={}; paintChips(); paintGrid(); };
  el('dynzoom').oninput=function(){ dynZoom=+this.value; el('dyngrid').style.setProperty('--dyncw', dynZoom+'px'); };
  el('dynDP').onchange=function(){ dynShowDP=this.checked; paintGrid(); };
  if(dynFields.length){
    Array.prototype.forEach.call(host.querySelectorAll('.dyn-filters select'),function(sel){ sel.onchange=function(){ var f=sel.getAttribute('data-df'); if(sel.value)dynFilter[f]=sel.value; else delete dynFilter[f]; recompute(); paintChips(); paintGrid(); }; });
    el('dynfclear').onclick=function(){ dynFilter={}; Array.prototype.forEach.call(host.querySelectorAll('.dyn-filters select'),function(s){s.value='';}); recompute(); paintChips(); paintGrid(); };
  }
  paintChips(); paintGrid();
}


var EPICOL={A:'#2f6fed',B:'#e6893a'};
var epiState={dir:'all',q:'',minr:null,conf:'all',view:'cards',tsort:{k:'q',asc:true}};
function epiMiniChart(p){
  var times=p.times||[], A=p.trajA||[], B=p.trajB||[], n=times.length;
  var W=250,H=150,ml=30,mr=12,mt=10,mb=26,pw=W-ml-mr,ph=H-mt-mb;
  function X(i){ return ml+(n<=1?pw/2:(i/(n-1))*pw); }
  function Y(a){ return mt+(1-a)*ph; }
  var svg='<svg viewBox="0 0 '+W+' '+H+'" width="100%" style="display:block"><title>Two allele-frequency trajectories over time; parallel lines = concordant, mirrored = discordant. Hover a point for its value.</title>';
  [0,0.5,1].forEach(function(a){ svg+='<line x1="'+ml+'" y1="'+Y(a).toFixed(1)+'" x2="'+(W-mr)+'" y2="'+Y(a).toFixed(1)+'" stroke="'+TH.grid+'"/><text x="'+(ml-5)+'" y="'+(Y(a)+3.5).toFixed(1)+'" text-anchor="end" font-size="10.5" fill="'+TH.mut+'">'+a.toFixed(1)+'</text>'; });
  [[A,EPICOL.A],[B,EPICOL.B]].forEach(function(pr){ var t=pr[0],c=pr[1],lastI=t.length-1; var pts=t.map(function(a,i){return X(i).toFixed(1)+','+Y(a).toFixed(1);}).join(' '); svg+='<polyline points="'+pts+'" fill="none" stroke="'+c+'" stroke-width="2.4" stroke-linecap="round" stroke-linejoin="round"/>'; t.forEach(function(a,i){ var last=(i===lastI); svg+='<circle cx="'+X(i).toFixed(1)+'" cy="'+Y(a).toFixed(1)+'" r="'+(last?4.2:3)+'" fill="'+c+'" stroke="'+TH.panel+'" stroke-width="'+(last?1.4:0)+'"><title>t='+esc(times[i]==null?i:times[i])+'  AF='+a.toFixed(3)+'</title></circle>'; }); });
  times.forEach(function(t,i){ svg+='<text x="'+X(i).toFixed(1)+'" y="'+(H-8)+'" text-anchor="middle" font-size="10.5" fill="'+TH.mut+'">'+esc(t==null?i:t)+'</text>'; });
  svg+='</svg>';
  return svg;
}
function epiColor(r){
  if(r==null) return '#f3f5f8';
  var a=Math.min(1,Math.abs(r));
  if(r>=0) return 'rgba(47,143,91,'+(0.10+a*0.82).toFixed(2)+')';   // concordant (green)
  return 'rgba(162,74,143,'+(0.10+a*0.82).toFixed(2)+')';           // discordant (purple)
}
function epiFmtR(v){ return (v<0?'-':(v>0?'+':''))+Math.abs(v).toFixed(2).replace(/^0/,''); }
function epiNodeLbl(nd){ return (nd.gene||'(intergenic)')+' '+String(nd.pos).split(':').pop(); }
function renderEpistasis(){
  var host=el('epi_body'), sec=el('epistasis'); if(!host)return;
  var E=R.epistasis;
  if(!(E&&E.pairs&&E.pairs.length)){ if(sec)sec.style.display='none'; var nv=el('nav-epi'); if(nv)nv.style.display='none'; return; }
  if(sec)sec.style.display='';
  if(epiState.minr==null) epiState.minr=E.min_r||0.8;
  var DIRS=[['all','all pairs'],['concordant','same dynamics'],['discordant','opposite dynamics']];
  var CONF=[['all','all pairs'],['mod','moderate or better (permutation p &#8804; 0.05)'],['strong','strong only (FDR q &#8804; 0.05)']];
  var VIEWS=[['cards','cards'],['matrix','matrix'],['table','table']];
  host.innerHTML=
    '<div class="epi-views">'+VIEWS.map(function(v){return '<button class="epi-viewbtn'+(epiState.view===v[0]?' on':'')+'" data-v="'+v[0]+'">'+v[1]+'</button>';}).join('')+'</div>'+
    '<div class="epi-controls">'+
      '<input id="epiq" class="dyn-search" type="search" title="Filter the pairs by gene name or position" placeholder="filter by gene / position..." value="'+esc(epiState.q)+'">'+
      '<span class="epi-flabel">dynamics</span>'+
      DIRS.map(function(d){return '<button class="dyn-btn epi-dirbtn'+(epiState.dir===d[0]?' on':'')+'" data-d="'+d[0]+'" title="Show '+d[1]+'">'+d[0]+'</button>';}).join('')+
      '<span class="epi-flabel">confidence</span>'+
      CONF.map(function(c){return '<button class="dyn-btn epi-confbtn'+(epiState.conf===c[0]?' on':'')+'" data-c="'+c[0]+'" title="'+c[1]+'">'+c[0]+'</button>';}).join('')+
      '<label class="dyn-zoom" title="Minimum |Pearson r| for a pair to be shown (cards / table)"><span>|r| &#8805;</span><input type="range" id="epir" min="'+(E.min_r||0.8)+'" max="0.99" step="0.01" value="'+epiState.minr+'"><b id="epirv">'+epiState.minr.toFixed(2)+'</b></label>'+
      '<span class="dyn-count" id="epicount"></span></div>'+
    '<div class="epi-legend">'+
      '<span title="The two variants rise and fall together across the series - candidate linkage or co-selection."><i style="background:#2f8f5b"></i>concordant / same dynamics <span class="infoi">i</span></span>'+
      '<span title="One variant rises as the other falls across the series - competing lineages / clonal interference."><i style="background:#a24a8f"></i>discordant / opposite dynamics <span class="infoi">i</span></span>'+
      '<span title="Permutation p-value: how often shuffling the timepoints of one trajectory reaches this |r| by chance. q = Benjamini-Hochberg FDR across every reported pair. strong = q &#8804; 0.05, moderate = p &#8804; 0.05, weak otherwise. With few timepoints a single series cannot beat ~1/n! by chance, so a pattern RECURRING across independent series is what drives a pair to strong."><b>strong</b> &#183; moderate &#183; weak = permutation p &amp; FDR q <span class="infoi">i</span></span>'+
      '<span class="c">mean Pearson r within a series (&#8805; '+E.min_points+' timepoints), across '+E.n_series+' series; '+E.perm+' permutations.</span>'+
    '</div>'+
    '<div class="dyn-grid" id="epi-cards"></div>'+
    '<div id="epi-matrix"></div>'+
    '<div id="epi-table"></div>';
  function posn(x){ return String(x).split(':').pop(); }
  function confok(p){ if(epiState.conf==='strong')return p.tier==='strong'; if(epiState.conf==='mod')return p.tier==='strong'||p.tier==='moderate'; return true; }
  function epiFilter(){
    var q=epiState.q.toLowerCase();
    return E.pairs.filter(function(p){
      if(Math.abs(p.r)<epiState.minr) return false;
      if(epiState.dir!=='all'&&p.direction!==epiState.dir) return false;
      if(!confok(p)) return false;
      if(q && !((p.geneA&&p.geneA.toLowerCase().indexOf(q)>=0)||(p.geneB&&p.geneB.toLowerCase().indexOf(q)>=0)||String(p.posA).indexOf(q)>=0||String(p.posB).indexOf(q)>=0)) return false;
      return true;
    });
  }
  function epiCards(list){
    var grid=el('epi-cards');
    if(!list.length){ grid.innerHTML='<div class="dyn-empty" style="grid-column:1/-1">&#128204; no variant pair matches the current filter.</div>'; return; }
    grid.innerHTML=list.map(function(p){
      var arrow=p.direction==='concordant'?'&#8596;':'&#8646;';
      var rec=(p.n>1)?('<span class="epi-recur" title="seen in '+p.n+' independent series'+(p.consistent?' with the same sign - recurrent':'')+'">&#8635; '+p.n+' series</span>'):('<span title="from a single series">series '+esc(p.group)+'</span>');
      return '<div class="epi-card '+p.direction+'">'+
        '<div class="epi-card-h"><span class="epi-badge '+p.direction+'">r = '+(p.r>0?'+':'')+p.r.toFixed(2)+'</span><span class="epi-tier epi-'+p.tier+'" title="permutation p = '+p.p+', FDR q = '+p.q+'">'+p.tier+'</span></div>'+
        '<div class="epi-pair"><span style="color:'+EPICOL.A+'"><b>'+esc(p.geneA||'(intergenic)')+'</b>'+geneRvTag(p.geneA)+' '+posn(p.posA)+(p.aaA?(' '+aaDual(p.aaA,p.aaA_h37rv)):'')+'</span><span class="epi-vs">'+arrow+'</span><span style="color:'+EPICOL.B+'"><b>'+esc(p.geneB||'(intergenic)')+'</b>'+geneRvTag(p.geneB)+' '+posn(p.posB)+(p.aaB?(' '+aaDual(p.aaB,p.aaB_h37rv)):'')+'</span></div>'+
        epiMiniChart(p)+
        '<div class="epi-card-f"><span title="permutation p-value / Benjamini-Hochberg FDR q-value">p '+p.p.toFixed(3)+' &#183; q '+p.q.toFixed(3)+'</span>'+rec+'</div>'+
      '</div>';
    }).join('');
  }
  function epiMatrix(){
    var box=el('epi-matrix'), M=E.matrix;
    if(!(M&&M.nodes&&M.nodes.length)){ box.innerHTML='<div class="dyn-empty">not enough correlated variants for a matrix.</div>'; return; }
    var nodes=M.nodes, N=nodes.length, showVal=N<=18;
    function cell(i,j){ if(i===j)return 1; var lo=Math.min(i,j),hi=Math.max(i,j); var v=M.cells[lo+','+hi]; return (v==null)?null:v; }
    var h='<div class="epimx-note">'+N+' variant(s)'+(M.truncated?(' (top '+N+' of '+M.total_nodes+' by connectivity)'):'')+' &#183; every pairwise mean r, including sub-threshold cells &#183; hover a cell for the pair.</div>';
    h+='<div class="epimx-wrap"><table class="epimx"><thead><tr><th class="epimx-corner"></th>';
    nodes.forEach(function(nd){ h+='<th class="epimx-hcell" title="'+esc(epiNodeLbl(nd))+(nd.aa?(' '+esc(nd.aa)):'')+'"><span class="epimx-h">'+esc(epiNodeLbl(nd))+'</span></th>'; });
    h+='</tr></thead><tbody>';
    nodes.forEach(function(nd,i){
      h+='<tr><th class="epimx-row" title="'+esc(epiNodeLbl(nd))+(nd.aa?(' '+esc(nd.aa)):'')+'">'+esc(epiNodeLbl(nd))+'</th>';
      nodes.forEach(function(nd2,j){
        if(i===j){ h+='<td class="epimx-cell epimx-diag" title="'+esc(epiNodeLbl(nd))+' (self)"></td>'; return; }
        var v=cell(i,j);
        var t=esc(epiNodeLbl(nd))+' '+(v!=null&&v<0?'&#8646;':'&#8596;')+' '+esc(epiNodeLbl(nd2))+(v==null?': no shared series':': r = '+epiFmtR(v));
        h+='<td class="epimx-cell" style="background:'+epiColor(v)+'" title="'+t+'">'+((showVal&&v!=null)?('<span>'+epiFmtR(v)+'</span>'):'')+'</td>';
      });
      h+='</tr>';
    });
    h+='</tbody></table></div>';
    h+='<div class="epimx-scale"><span>discordant &#8722;1</span><i class="epimx-grad"></i><span>+1 concordant</span></div>';
    box.innerHTML=h;
  }
  function tsortVal(p,k){ if(k==='pair')return (p.geneA||'')+p.posA; if(k==='r')return p.r; if(k==='ar')return Math.abs(p.r); if(k==='n')return p.n; if(k==='p')return p.p; if(k==='q')return p.q; if(k==='tier')return ({strong:0,moderate:1,weak:2})[p.tier]; return p.direction; }
  function epiTable(list){
    var box=el('epi-table');
    var COLS=[['pair','Variant A &#8596; Variant B'],['direction','dynamics'],['ar','|r|'],['r','r'],['n','series'],['p','p'],['q','q (FDR)'],['tier','confidence']];
    var k=epiState.tsort.k, asc=epiState.tsort.asc;
    var rows=list.slice().sort(function(a,b){ var x=tsortVal(a,k),y=tsortVal(b,k),c; if(typeof x==='number'&&typeof y==='number')c=x-y; else c=String(x).localeCompare(String(y)); return asc?c:-c; });
    var h='<div class="epitbl-top"><button class="dyn-btn" id="epidl" title="Download every reported pair as a TSV">'+icon('download')+'download pairs (TSV)</button><span class="dyn-count">'+rows.length+' pair(s)</span></div>';
    h+='<div class="epitbl-wrap"><table class="epitbl"><thead><tr>'+COLS.map(function(c){return '<th data-k="'+c[0]+'">'+c[1]+(k===c[0]?(asc?icon('chevronUp','sort'):icon('chevronDown','sort')):'')+'</th>';}).join('')+'</tr></thead><tbody>';
    if(!rows.length){ h+='<tr><td colspan="'+COLS.length+'" class="c" style="padding:20px;text-align:center">no variant pair matches the current filter.</td></tr>'; }
    rows.forEach(function(p){
      var a='<b style="color:'+EPICOL.A+'">'+esc(p.geneA||'(intergenic)')+'</b>'+geneRvTag(p.geneA)+' '+String(p.posA).split(':').pop()+(p.aaA?(' '+aaDual(p.aaA,p.aaA_h37rv)):'');
      var b='<b style="color:'+EPICOL.B+'">'+esc(p.geneB||'(intergenic)')+'</b>'+geneRvTag(p.geneB)+' '+String(p.posB).split(':').pop()+(p.aaB?(' '+aaDual(p.aaB,p.aaB_h37rv)):'');
      h+='<tr><td>'+a+' <span class="epi-vs">'+(p.direction==='concordant'?'&#8596;':'&#8646;')+'</span> '+b+'</td>'+
        '<td>'+p.direction+'</td>'+
        '<td class="epitbl-r" style="color:'+(p.r>=0?'#2f8f5b':'#a24a8f')+'">'+(p.r>0?'+':'')+p.r.toFixed(2)+'</td>'+
        '<td>'+(p.r>0?'+':'')+p.r.toFixed(2)+'</td>'+
        '<td>'+p.n+(p.n>1?'&#8635;':'')+'</td>'+
        '<td>'+p.p.toFixed(3)+'</td>'+
        '<td>'+p.q.toFixed(3)+'</td>'+
        '<td><span class="epi-tier epi-'+p.tier+'">'+p.tier+'</span></td></tr>';
    });
    h+='</tbody></table></div>';
    box.innerHTML=h;
    Array.prototype.forEach.call(box.querySelectorAll('th[data-k]'),function(th){ th.onclick=function(){ var kk=th.getAttribute('data-k'); if(epiState.tsort.k===kk)epiState.tsort.asc=!epiState.tsort.asc; else{epiState.tsort.k=kk;epiState.tsort.asc=(kk==='pair'||kk==='direction'||kk==='tier');} epiTable(epiFilter()); }; });
    el('epidl').onclick=function(){
      var hdr=['geneA','posA','aa_A','geneB','posB','aa_B','dynamics','mean_r','n_series','r_min','r_max','perm_p','fdr_q','confidence'];
      var lines=[hdr.join('\t')];
      E.pairs.forEach(function(p){ lines.push([p.geneA,p.posA,p.aaA,p.geneB,p.posB,p.aaB,p.direction,p.r,p.n,p.rmin,p.rmax,p.p,p.q,p.tier].join('\t')); });
      dl(lines.join('\n')+'\n','epistasis_pairs.tsv','text/tab-separated-values');
    };
  }
  function draw(){
    var list=epiFilter();
    el('epicount').innerHTML=(epiState.view==='matrix')?((E.matrix?E.matrix.nodes.length:0)+' variant(s) in the matrix'):(list.length+' pair(s)'+((epiState.dir==='all'&&epiState.conf==='all')?(' &#183; '+E.n_concordant+' concordant / '+E.n_discordant+' discordant &#183; '+E.n_strong+' strong'):''));
    el('epi-cards').style.display=epiState.view==='cards'?'':'none';
    el('epi-matrix').style.display=epiState.view==='matrix'?'':'none';
    el('epi-table').style.display=epiState.view==='table'?'':'none';
    if(epiState.view==='cards') epiCards(list);
    else if(epiState.view==='matrix') epiMatrix();
    else epiTable(list);
  }
  el('epiq').oninput=function(){ epiState.q=this.value; draw(); };
  Array.prototype.forEach.call(host.querySelectorAll('.epi-viewbtn'),function(b){ b.onclick=function(){ epiState.view=b.getAttribute('data-v'); Array.prototype.forEach.call(host.querySelectorAll('.epi-viewbtn'),function(x){x.className='epi-viewbtn'+(x.getAttribute('data-v')===epiState.view?' on':'');}); draw(); }; });
  Array.prototype.forEach.call(host.querySelectorAll('.epi-dirbtn'),function(b){ b.onclick=function(){ epiState.dir=b.getAttribute('data-d'); Array.prototype.forEach.call(host.querySelectorAll('.epi-dirbtn'),function(x){x.className='dyn-btn epi-dirbtn'+(x.getAttribute('data-d')===epiState.dir?' on':'');}); draw(); }; });
  Array.prototype.forEach.call(host.querySelectorAll('.epi-confbtn'),function(b){ b.onclick=function(){ epiState.conf=b.getAttribute('data-c'); Array.prototype.forEach.call(host.querySelectorAll('.epi-confbtn'),function(x){x.className='dyn-btn epi-confbtn'+(x.getAttribute('data-c')===epiState.conf?' on':'');}); draw(); }; });
  el('epir').oninput=function(){ epiState.minr=+this.value; el('epirv').textContent=epiState.minr.toFixed(2); draw(); };
  draw();
}

function snpAfColor(af){return 'rgba(31,120,180,'+(0.16+af*0.8).toFixed(2)+')';}
var SNPMX_PAL=['#bcd0ea','#f3d1b0','#c3e0c9','#f0c4cf','#d6c9ec','#b8e0dd','#eadfb0','#dfe4ea','#f2c4c4','#cdd1a8','#e6c3e0','#b9d6ee'];
// dark-theme categorical palette for the SNP-matrix metadata rows (same hues, muted/dark so they don't glare as light bands)
var SNPMX_PAL_DARK=['#2f4a6b','#6b4c2f','#2f5a40','#6b3a4a','#4a3a6b','#2f5a55','#5c4a2f','#3a4450','#6b3838','#4c4c2f','#5a2f5a','#35506b'];
var snpmxFilter={};
function renderSnpMatrix(){
  var host=el('snpmx_body'), sec=el('snpmatrix'); if(!host)return;
  var M=R.snp_matrix;
  if(!(M&&M.rows&&M.rows.length)){ if(sec)sec.style.display='none'; var nv=el('nav-snpmx'); if(nv)nv.style.display='none'; return; }
  if(sec)sec.style.display='';
  var samples=M.samples, MAXR=400;
  var meta=(R.sample_meta&&R.sample_meta.fields&&R.sample_meta.fields.length)?R.sample_meta:null;
  var metaMaps={}, metaVals={};
  var PAL=isDark()?SNPMX_PAL_DARK:SNPMX_PAL;
  if(meta){ meta.fields.forEach(function(f){ var m={},vals=[],k=0; samples.forEach(function(s){var v=(meta.rows[s]||{})[f]; if(v&&!(v in m)){m[v]=PAL[k%PAL.length];k++;vals.push(v);}}); metaMaps[f]=m; metaVals[f]=vals.sort(); }); }
  // a field is "the lineage field" only if every (non-NA) value is a known lineage -> reuse the canonical
  // report palette (mycolorsTB / linColor) for it; all other fields get the neutral pastels. Keying on the
  // field (not just the value) stops an unrelated column whose value happens to equal a lineage label from
  // being mis-coloured.
  var linField={};
  if(meta){ meta.fields.forEach(function(f){ var vs=(metaVals[f]||[]).filter(function(v){return v&&v!=='NA'&&v!=='.'&&v!=='-';}); if(vs.length&&vs.every(function(v){return LINCOL[v];})) linField[f]=1; }); }
  function metaColor(f,v){ if(v&&linField[f]&&LINCOL[v]) return LINCOL[v]; return (v&&metaMaps[f]&&metaMaps[f][v])?metaMaps[f][v]:TH.cellnull; }
  function metaText(bg){
    if(bg&&bg.charAt(0)==='#'){ var h=bg.length===4?('#'+bg.charAt(1)+bg.charAt(1)+bg.charAt(2)+bg.charAt(2)+bg.charAt(3)+bg.charAt(3)):bg;
      var L=(0.299*parseInt(h.substr(1,2),16)+0.587*parseInt(h.substr(3,2),16)+0.114*parseInt(h.substr(5,2),16))/255; return L<0.62?'#fff':'#1c2b3a'; }
    var m=/hsl\(\s*[\d.]+\s*,\s*[\d.]+%\s*,\s*([\d.]+)%/i.exec(bg||'');   // algorithmic lineage fallback is hsl(h,58%,52%)
    if(m) return (+m[1])<62?'#fff':'#1c2b3a';
    return '#1c2b3a';
  }
  function visIdx(){ var idx=[]; samples.forEach(function(s,i){ var ok=true; if(meta){ for(var f in snpmxFilter){ if(snpmxFilter[f] && (meta.rows[s]||{})[f]!==snpmxFilter[f]){ ok=false; break; } } } if(ok) idx.push(i); }); return idx; }
  var filterUI=meta?('<div class="snpmx-filters"><span class="snpmx-flabel">filter columns:</span>'+
    meta.fields.map(function(f){ return '<label class="snpmx-fsel">'+esc(f)+' <select data-f="'+esc(f)+'"><option value="">all</option>'+
      metaVals[f].map(function(v){return '<option value="'+esc(v)+'"'+(snpmxFilter[f]===v?' selected':'')+'>'+esc(v)+'</option>';}).join('')+'</select></label>'; }).join('')+
    '<button class="dyn-btn" id="snpmxfclear">clear</button></div>'):'';
  host.innerHTML=
    '<div class="snpmx-controls">'+
      '<input id="snpmxq" class="dyn-search" type="search" placeholder="filter by gene / position / amino acid...">'+
      '<label class="snpmx-toggle"><input type="checkbox" id="snpmxdp" checked> show depth</label>'+
      '<button class="dyn-btn" id="snpmxdl" title="Download the full matrix (all samples) as a wide TSV">'+icon('download')+'download matrix (TSV)</button>'+
      '<span class="dyn-count" id="snpmxcount"></span></div>'+
    filterUI+
    (meta?('<div class="snpmx-metanote">column levels from the samplesheet: '+meta.fields.map(function(f){return '<b>'+esc(f)+'</b>';}).join(' &#183; ')+' &#183; hover a header cell for its value</div>'):'')+
    '<div class="snpmx-wrap"><table class="snpmx" id="snpmxtable"></table></div>';
  function draw(){
    var q=(el('snpmxq').value||'').toLowerCase(), showDP=el('snpmxdp').checked;
    var vi=visIdx();
    var rows=M.rows.filter(function(r){
      if(q && !((r.gene&&r.gene.toLowerCase().indexOf(q)>=0)||String(r.pos).indexOf(q)>=0||(r.aa&&r.aa.toLowerCase().indexOf(q)>=0))) return false;
      return vi.some(function(i){return r.cells[i];});   // only SNPs seen in the visible columns
    });
    var nfilt=0; for(var kf in snpmxFilter){ if(snpmxFilter[kf]) nfilt++; }
    el('snpmxcount').innerHTML=rows.length+' SNP site(s) &#215; '+vi.length+' sample(s)'+(nfilt?' (filtered)':'')+(rows.length>MAXR?(' - showing '+MAXR):'')+(M.truncated?' &#183; full matrix in the TSV':'');
    var shown=rows.slice(0,MAXR), mh=22, nf=meta?meta.fields.length:0;
    var metaRows=meta?meta.fields.map(function(f,k){
      return '<tr>'+'<th class="snpmx-info snpmx-metalabel" style="top:'+(k*mh)+'px">'+esc(f)+'</th>'+
        vi.map(function(i){var s=samples[i],v=(meta.rows[s]||{})[f]||'',bg=metaColor(f,v); return '<th class="snpmx-metacell" style="top:'+(k*mh)+'px;background:'+bg+';color:'+metaText(bg)+'" title="'+esc(f)+': '+esc(v||'-')+'">'+esc(v)+'</th>';}).join('')+'</tr>';
    }).join(''):'';
    var stop=nf*mh;
    var nameRow='<tr><th class="snpmx-info snpmx-corner" style="top:'+stop+'px">SNP '+esc(M.reference?('('+M.reference+')'):'')+'</th>'+
      vi.map(function(i){var s=samples[i];return '<th class="snpmx-hcell" style="top:'+stop+'px" title="'+esc(s)+'"><span class="snpmx-h">'+esc(s)+'</span></th>';}).join('')+'</tr>';
    var body='<tbody>'+shown.map(function(r){
      var lbl='<b>'+esc(r.gene||r.contig)+'</b>'+geneRvTag(r.gene)+' '+r.pos+' '+esc(r.ref)+'&#8594;'+esc(r.alt)+(r.aa?(' <span class="snpmx-aa">'+aaDual(r.aa,r.aa_h37rv)+'</span>'):'');
      var cells=vi.map(function(i){var c=r.cells[i], s=samples[i];
        if(!c) return '<td class="snpmx-cell snpmx-empty" title="'+esc(s)+' - not called"></td>';
        var afTxt=c[0].toFixed(2).replace(/^0/,'').replace(/^1\.00$/,'1');
        var dpTxt=(showDP&&c[1]!=null)?('<span class="snpmx-dp">'+c[1]+'</span>'):'';
        return '<td class="snpmx-cell'+(showDP?' wdp':'')+'" style="background:'+snpAfColor(c[0])+'" title="'+esc(s)+'  AF='+c[0].toFixed(3)+(c[1]!=null?('  DP='+c[1]):'')+'"><span class="snpmx-af">'+afTxt+'</span>'+dpTxt+'</td>';
      }).join('');
      return '<tr><td class="snpmx-info">'+lbl+'</td>'+cells+'</tr>';
    }).join('')+'</tbody>';
    el('snpmxtable').innerHTML='<thead>'+metaRows+nameRow+'</thead>'+body;
  }
  el('snpmxq').oninput=function(){clearTimeout(_mxdb);_mxdb=setTimeout(draw,160);};
  el('snpmxdp').onchange=draw;
  if(meta){
    Array.prototype.forEach.call(host.querySelectorAll('.snpmx-filters select'),function(sel){ sel.onchange=function(){ var f=sel.getAttribute('data-f'); if(sel.value) snpmxFilter[f]=sel.value; else delete snpmxFilter[f]; draw(); }; });
    el('snpmxfclear').onclick=function(){ snpmxFilter={}; Array.prototype.forEach.call(host.querySelectorAll('.snpmx-filters select'),function(s){s.value='';}); draw(); };
  }
  el('snpmxdl').onclick=function(){
    var hdr=['reference','contig','pos','ref_allele','alt_allele','gene','effect','aa_change'];
    samples.forEach(function(s){hdr.push(s+'|AF');hdr.push(s+'|DP');});
    var lines=[hdr.join('\t')];
    if(meta){ meta.fields.forEach(function(f){ var row=['# '+f,'','','','','','','']; samples.forEach(function(s){row.push((meta.rows[s]||{})[f]||'');row.push('');}); lines.push(row.join('\t')); }); }
    M.rows.forEach(function(r){var row=[M.reference||r.contig,r.contig,r.pos,r.ref,r.alt,r.gene,r.eff,r.aa];
      samples.forEach(function(s,i){var c=r.cells[i]; if(c){row.push(c[0].toFixed(4));row.push(c[1]==null?'':c[1]);}else{row.push('');row.push('');}});
      lines.push(row.join('\t'));});
    dl(lines.join('\n')+'\n','snp_matrix.tsv','text/tab-separated-values');
  };
  draw();
}

function drGColor(gn){ return (gn===1||gn===2)?'#dc2626':(gn===3?'#d97706':((gn===4||gn===5)?'#94a3b8':'#b8c2cf')); }
function drStatus(gns){ var r=false,u=false,n=false; for(var i=0;i<gns.length;i++){var g=gns[i]; if(g===1||g===2)r=true; else if(g===3)u=true; else if(g===4||g===5)n=true;} return r?{t:'R',c:'#dc2626'}:(u?{t:'?',c:'#d97706'}:(n?{t:'&#183;',c:'#94a3b8'}:{t:'',c:'#eef2f7'})); }
// ---- Kraken2 taxonomic composition (contamination / host check on the reads before mapping) ----
function renderKraken(){
  var host=el('kraken_body'), sec=el('kraken'); if(!host)return;
  var K=R.kraken;
  if(!(K&&K.samples&&K.samples.length)){ if(sec)sec.style.display='none'; var nv=el('nav-kraken'); if(nv)nv.style.display='none'; return; }
  var rows=K.samples.slice().sort(function(a,b){ return (a.primary?a.primary.pct:0)-(b.primary?b.primary.pct:0); });  // most-contaminated (lowest primary %) first
  // ---- cohort stacked-composition plot: one bar per sample, full taxa breakdown ----
  var agg={}; rows.forEach(function(s){ (s.top||[]).forEach(function(t){ agg[t.name]=(agg[t.name]||0)+t.pct; }); });
  var taxa=Object.keys(agg).sort(function(a,b){ return agg[b]-agg[a]; }).slice(0,8);
  var TAXPAL=['#3b7dd8','#e0544f','#2ea36b','#e0a11f','#8a63c9','#26a0a0','#d06fae','#c98a3b'], OTHERC='#9aa7b6';
  var tcol={}; taxa.forEach(function(t,i){ tcol[t]=TAXPAL[i%TAXPAL.length]; });
  var kW=Math.max(320,(host.clientWidth||760)), kLab=Math.min(150,Math.round(kW*0.28)), kRP=10, kBx=kLab, kBw=kW-kLab-kRP,
      kRowH=Math.max(15,Math.min(22,Math.floor(340/rows.length))), kTop=6, kH=kTop+rows.length*kRowH+20;
  var ksvg='<svg width="'+kW+'" height="'+kH+'" style="display:block;max-width:100%">';
  [0,0.5,1].forEach(function(f){ var x=kBx+f*kBw; ksvg+='<line x1="'+x.toFixed(1)+'" y1="'+kTop+'" x2="'+x.toFixed(1)+'" y2="'+(kTop+rows.length*kRowH).toFixed(1)+'" stroke="'+TH.grid+'"/><text x="'+x.toFixed(1)+'" y="'+(kH-6)+'" text-anchor="'+(f===0?'start':f===1?'end':'middle')+'" font-size="9.5" fill="'+TH.mut+'">'+(f*100)+'%</text>'; });
  rows.forEach(function(s,i){
    var y=kTop+i*kRowH, cx=kBx, used=0, arr=s.top||[];
    ksvg+='<g data-s="'+esc(s.s)+'" style="cursor:pointer"><text x="'+(kLab-6)+'" y="'+(y+kRowH/2+3).toFixed(1)+'" text-anchor="end" font-size="10" fill="'+TH.ink+'">'+esc(s.s.length>20?s.s.slice(0,19)+'…':s.s)+'</text>';
    taxa.forEach(function(t){ var m=null; for(var j=0;j<arr.length;j++){ if(arr[j].name===t){ m=arr[j]; break; } }
      if(m&&m.pct>0){ var w=m.pct/100*kBw; ksvg+='<rect x="'+cx.toFixed(1)+'" y="'+(y+2).toFixed(1)+'" width="'+Math.max(0.4,w).toFixed(1)+'" height="'+(kRowH-4)+'" fill="'+tcol[t]+'"><title>'+esc(s.s)+' · '+esc(t)+' '+m.pct.toFixed(1)+'%</title></rect>'; cx+=w; used+=m.pct; } });
    var other=Math.max(0,(s.classified||0)-used); if(other>0.05){ var wo=other/100*kBw; ksvg+='<rect x="'+cx.toFixed(1)+'" y="'+(y+2).toFixed(1)+'" width="'+wo.toFixed(1)+'" height="'+(kRowH-4)+'" fill="'+OTHERC+'"><title>'+esc(s.s)+' · other classified '+other.toFixed(1)+'%</title></rect>'; cx+=wo; }
    var unc=s.unclassified||0; if(unc>0.05){ var wu=unc/100*kBw; ksvg+='<rect x="'+cx.toFixed(1)+'" y="'+(y+2).toFixed(1)+'" width="'+wu.toFixed(1)+'" height="'+(kRowH-4)+'" fill="'+TH.track+'"><title>'+esc(s.s)+' · unclassified '+unc.toFixed(1)+'%</title></rect>'; }
    ksvg+='</g>';
  });
  ksvg+='</svg>';
  var kleg='<div class="krk-legend">'+taxa.map(function(t){ return '<span><i style="background:'+tcol[t]+'"></i>'+esc(t)+'</span>'; }).join('')+
    '<span><i style="background:'+OTHERC+'"></i>other classified</span><span><i style="background:'+TH.track+'"></i>unclassified</span> <span class="krk-mut">- one bar per sample, worst first; hover a segment for the %. Click a bar/row to highlight the sample everywhere.</span></div>';
  var body=rows.map(function(s){
    var pri=s.primary?s.primary.pct:0, unc=s.unclassified||0, other=Math.max(0,100-pri-unc), lowPri=pri<90, hiUnc=unc>15;
    var bar='<div class="krk-bar">'+
      '<div class="krk-seg" style="width:'+pri.toFixed(1)+'%;background:'+(lowPri?'#e0a11f':'#2ea36b')+'" title="primary: '+esc(s.primary?s.primary.name:'-')+' '+pri.toFixed(1)+'%"></div>'+
      '<div class="krk-seg" style="width:'+other.toFixed(1)+'%;background:#c0704f" title="other classified '+other.toFixed(1)+'%"></div>'+
      '<div class="krk-seg" style="width:'+unc.toFixed(1)+'%;background:var(--track)" title="unclassified '+unc.toFixed(1)+'%"></div></div>';
    return '<tr class="hit" data-s="'+esc(s.s)+'"'+(st.hi==s.s?' style="background:'+TH.hl+'"':'')+'><td class="s">'+esc(s.s)+'</td>'+
      '<td style="text-align:left"><b>'+esc(s.primary?s.primary.name:'-')+'</b></td>'+
      '<td'+(lowPri?' style="color:var(--fail);font-weight:600"':'')+'>'+pri.toFixed(1)+'</td>'+
      '<td style="text-align:left">'+((s.secondary&&s.secondary.pct>=1)?esc(s.secondary.name)+' <span class="krk-mut">'+s.secondary.pct.toFixed(1)+'%</span>':'<span class="krk-mut">-</span>')+'</td>'+
      '<td'+(hiUnc?' style="color:var(--warn)"':'')+'>'+unc.toFixed(1)+'</td>'+
      '<td class="krk-barcell">'+bar+'</td></tr>';
  }).join('');
  host.innerHTML='<div class="krk-chart">'+kleg+'<div class="krk-plotscroll">'+ksvg+'</div></div>'+
    '<div class="gtable" style="margin-top:14px"><table class="krktable"><thead><tr><th class="s">Sample</th><th style="text-align:left">Primary taxon</th><th>Primary %</th><th style="text-align:left">Top other</th><th>Unclass. %</th><th style="text-align:left">Composition</th></tr></thead><tbody>'+body+'</tbody></table></div>';
  Array.prototype.forEach.call(host.querySelectorAll('[data-s]'),function(e){e.onclick=function(){setHi(e.getAttribute('data-s'));};});
}
var drState={q:''};
function renderDrug(){
  var host=el('drug_body'), sec=el('drug'); if(!host)return;
  var D=R.dr;
  if(!(D&&D.calls&&D.calls.length)){ if(sec)sec.style.display='none'; var nv=el('nav-drug'); if(nv)nv.style.display='none'; return; }
  if(sec)sec.style.display='';
  var samples=D.samples, drugs=D.drugs, calls=D.calls;
  var cell={}, cmut={};
  calls.forEach(function(c){ (cell[c.s]=cell[c.s]||{}); (cell[c.s][c.drug]=cell[c.s][c.drug]||[]).push(c.gn); (cmut[c.s]=cmut[c.s]||{}); (cmut[c.s][c.drug]=cmut[c.s][c.drug]||[]).push(c); });
  var mx='<div class="dr-mxwrap"><table class="drmx"><thead><tr><th class="dr-corner">sample \\ drug</th>'+
    drugs.map(function(dr){return '<th class="dr-hcell" title="'+esc(dr)+'"><span class="dr-h">'+esc(dr)+'</span></th>';}).join('')+'</tr></thead><tbody>'+
    samples.map(function(s){ return '<tr><th class="dr-row" title="'+esc(s)+'">'+esc(s)+'</th>'+drugs.map(function(dr){
      var gns=(cell[s]||{})[dr];
      if(!gns) return '<td class="drmx-cell" title="'+esc(s)+' &#183; '+esc(dr)+': no mutation detected"></td>';
      var st=drStatus(gns);
      var muts=((cmut[s]||{})[dr]||[]).map(function(c){return c.gene+' '+c.mutation+(c.gn?(' (WHO '+c.gn+')'):'');}).join('; ');
      return '<td class="drmx-cell" style="background:'+st.c+'" title="'+esc(s)+' &#183; '+esc(dr)+' &#8212; '+esc(muts)+'"><b>'+st.t+'</b></td>';
    }).join('')+'</tr>'; }).join('')+'</tbody></table></div>';
  host.innerHTML=
    '<div class="dr-legend">'+
      '<span title="WHO groups 1-2: associated with resistance"><i style="background:#dc2626"></i>1&#8211;2 associated with R</span>'+
      '<span title="WHO group 3: uncertain significance"><i style="background:#d97706"></i>3 uncertain</span>'+
      '<span title="WHO groups 4-5: not associated with resistance"><i style="background:#94a3b8"></i>4&#8211;5 not associated</span>'+
      '<span class="c">Cell = worst grade per drug (R / ? / &#183;). A genomic screen (pathotypr, WHO catalogue, H37Rv numbering), not a clinical DST result.</span></div>'+
    mx+
    '<div class="dr-controls"><input id="drq" class="dyn-search" type="search" placeholder="filter by sample / drug / gene / mutation..." value="'+esc(drState.q)+'"><button class="dyn-btn" id="drdl" title="Download every resistance call as a TSV">'+icon('download')+'download calls (TSV)</button><span class="dyn-count" id="drcount"></span></div>'+
    '<div class="epitbl-wrap"><table class="epitbl" id="drtable"></table></div>';
  function draw(){
    var q=drState.q.toLowerCase();
    var rows=calls.filter(function(c){ return !q||(c.s.toLowerCase().indexOf(q)>=0)||(c.drug.toLowerCase().indexOf(q)>=0)||(c.gene.toLowerCase().indexOf(q)>=0)||(c.mutation.toLowerCase().indexOf(q)>=0); });
    rows=rows.slice().sort(function(a,b){ if(a.s!==b.s) return a.s<b.s?-1:1; return (a.gn||9)-(b.gn||9); });
    el('drcount').textContent=rows.length+' call(s)';
    var h='<thead><tr><th>Sample</th><th>Drug</th><th>Gene</th><th>Mutation (H37Rv)</th><th>WHO grade</th><th>AF</th><th>DP</th></tr></thead><tbody>';
    if(!rows.length) h+='<tr><td colspan="7" class="c" style="padding:18px;text-align:center">no call matches the filter.</td></tr>';
    h+=rows.slice(0,600).map(function(c){
      return '<tr><td>'+esc(c.s)+'</td><td><b>'+esc(c.drug)+'</b></td><td>'+esc(c.gene)+geneRvTag(c.gene)+'</td><td class="epitbl-r">'+esc(c.mutation)+'</td>'+
        '<td><span class="dr-badge" style="background:'+drGColor(c.gn)+'" title="'+esc(c.marker||'')+'">'+esc(c.grade||'?')+'</span></td>'+
        '<td>'+(c.af==null?'':c.af.toFixed(2))+'</td><td>'+(c.dp==null?'':c.dp)+'</td></tr>';
    }).join('')+'</tbody>';
    el('drtable').innerHTML=h;
  }
  el('drq').oninput=function(){ drState.q=this.value; draw(); };
  el('drdl').onclick=function(){
    var hdr=['sample','drug','gene','mutation_h37rv','who_grade','marker','af','dp'];
    var lines=[hdr.join('\t')];
    calls.forEach(function(c){ lines.push([c.s,c.drug,c.gene,c.mutation,c.grade,c.marker,(c.af==null?'':c.af),(c.dp==null?'':c.dp)].join('\t')); });
    dl(lines.join('\n')+'\n','drug_resistance.tsv','text/tab-separated-values');
  };
  draw();
}

function renderAll(){renderExec();renderOverview();renderTable();renderLineages();renderPlots();renderScatter();renderCorr();renderQCspace();renderRefBias();renderStacks();renderGenome();renderFunction();renderGeneBurden();renderHotspots();renderTemporal();renderPnps();renderADNA();renderDynamics();renderEpistasis();renderSnpMatrix();renderDrug();renderKraken();renderFlags();renderCuration();}

// ---- static wiring ----
el('meta').textContent=R.samples.length+' samples · '+R.generated;
(function(){var hv=el('hver'); if(hv){ if(R.version){hv.textContent='v'+R.version; hv.title='BAMpiro version '+R.version;} else hv.style.display='none'; }})();
fillIcons();   // swap every static data-ic placeholder (header, TOC chevrons, buttons, modal) for its inline SVG
// click an (i) info icon -> show a persistent popover with its definition (capture phase so it
// beats the column-sort handler); click anywhere / Esc / scroll to dismiss.
(function(){
  var pop=el('infopop');
  document.addEventListener('click',function(e){
    var ic=(e.target&&e.target.closest)?e.target.closest('.infoi'):null;
    if(ic){
      e.stopPropagation(); e.preventDefault();
      if(pop._for===ic&&pop.style.display==='block'){ pop.style.display='none'; pop._for=null; return; }
      var _h=ic.closest('[title]'); pop.textContent=ic.getAttribute('data-info')||(_h?_h.getAttribute('title'):'')||'';
      var pw=Math.min(320,window.innerWidth-24); pop.style.maxWidth=pw+'px'; pop.style.display='block';
      var r=ic.getBoundingClientRect();
      var left=Math.min(Math.max(8,r.left-4),window.innerWidth-pw-8), top=r.bottom+8;
      if(top+pop.offsetHeight+8>window.innerHeight){ var up=r.top-8-pop.offsetHeight; if(up>4)top=up; }
      pop.style.left=left+'px'; pop.style.top=top+'px'; pop._for=ic;
    } else if(pop.style.display==='block'){ pop.style.display='none'; pop._for=null; }
  },true);
  document.addEventListener('keydown',function(e){ if(e.key==='Escape'&&pop.style.display==='block'){pop.style.display='none';pop._for=null;} });
  window.addEventListener('scroll',function(){ if(pop.style.display==='block'){pop.style.display='none';pop._for=null;} },true);
})();
el('foot').innerHTML='<span class="foot-brand">BAMpiro'+(R.version?' <b>v'+esc(R.version)+'</b>':'')+' · <a href="'+esc(R.repo_url)+'" target="_blank" rel="noopener noreferrer">'+icon('github','sort')+'source on GitHub'+icon('ext','sort')+'</a></span> · Generated '+R.generated+' · thresholds are adjustable live above; the pipeline gate uses the defaults ('+
  Object.keys(R.thresholds).map(function(k){return k+'='+R.thresholds[k];}).join(', ')+'). Values scale within each column; NA = not reported.';
el('colmenu').innerHTML='<div style="display:flex;gap:12px;margin-bottom:6px;padding-bottom:6px;border-bottom:1px solid var(--line);font-size:12px"><a href="#" id="colall" style="color:var(--accent)">show all</a><a href="#" id="colnone" style="color:var(--accent)">hide all</a></div>'+R.metrics.map(function(m){return '<label><input type="checkbox" data-k="'+m.key+'"'+(st.hidden[m.key]?'':' checked')+'> '+esc(m.label)+'</label>';}).join('');
Array.prototype.forEach.call(document.querySelectorAll('#colmenu input'),function(cb){cb.onchange=function(){if(cb.checked)delete st.hidden[cb.getAttribute('data-k')];else st.hidden[cb.getAttribute('data-k')]=1;renderTable();saveState();};});
(function(){var ca=el('colall'),cn=el('colnone');
  if(ca)ca.onclick=function(e){e.preventDefault();st.hidden={};Array.prototype.forEach.call(document.querySelectorAll('#colmenu input'),function(cb){cb.checked=true;});renderTable();saveState();};
  if(cn)cn.onclick=function(e){e.preventDefault();Array.prototype.forEach.call(document.querySelectorAll('#colmenu input'),function(cb){cb.checked=false;st.hidden[cb.getAttribute('data-k')]=1;});renderTable();saveState();};})();
el('q').oninput=function(e){st.q=e.target.value.toLowerCase().trim();renderTable();clearTimeout(_qdb);_qdb=setTimeout(function(){renderPlots();renderScatter();renderCorr();renderQCspace();renderRefBias();renderStacks();renderGenome();renderFunction();renderTemporal();},160);};
// per-panel gene search (Functional gene burden / Variable genes / pN-pS): filter each gene table by gene name
[['gbq','gbq',renderGeneBurden],['hotq','hotq',renderHotspots],['pnpsq','pnpsq',renderPnps]].forEach(function(w){var inp=el(w[0]);if(inp)inp.oninput=function(e){st[w[1]]=e.target.value.trim();w[2]();};});
el('of').onchange=function(e){st.onlyFlagged=e.target.checked;renderAll();};
Array.prototype.forEach.call(document.querySelectorAll('#ptype button'),function(b){b.onclick=function(){st.ptype=b.getAttribute('data-t');
  Array.prototype.forEach.call(document.querySelectorAll('#ptype button'),function(x){x.classList.toggle('on',x==b);});renderPlots();};});
// scatter axis pickers
['sx','sy'].forEach(function(ax){var sel=el(ax);sel.innerHTML=R.metrics.map(function(m){return '<option value="'+m.key+'"'+(st[ax]==m.key?' selected':'')+'>'+esc(m.label)+'</option>';}).join('');
  sel.onchange=function(){st[ax]=sel.value;renderScatter();};});
// live thresholds + presets
var THL=[['depth_min','Depth min'],['breadth_min','Breadth min %'],['missing_max','Missing max %'],['mapping_min','Mapped min %'],['dup_max','Dup max %'],['iupac_max','IUPAC max %'],['titv_min','Ti/Tv min'],['snp_z','SNP z'],['het_max_frac','Het % max'],['mixed_min_frac','Mixed lin % min']];
var PRESETS=[
  {id:'gate',label:'gate defaults',th:assign({},R.thresholds),note:'The cut-offs the Snakemake qc_gate uses (config report_* keys).'},
  {id:'strict',label:'strict (modern WGS)',th:{depth_min:20,breadth_min:95,missing_max:5,mapping_min:90,dup_max:30,iupac_max:2,titv_min:1.5,snp_z:3,het_max_frac:1.5,mixed_min_frac:2},note:'Confident modern Illumina isolate: 20x, 95% breadth, <5% missing, Ti/Tv >=1.5.'},
  {id:'lenient',label:'lenient (aDNA / low-cov)',th:{depth_min:3,breadth_min:60,missing_max:40,mapping_min:50,dup_max:80,iupac_max:5,titv_min:0,snp_z:4,het_max_frac:8,mixed_min_frac:5},note:'Degraded / low-coverage library at the 3x calling floor: rescues calibration tips.'}
];
function applyThr(next){Object.keys(next).forEach(function(k){if(k in thr)thr[k]=next[k];});
  Array.prototype.forEach.call(document.querySelectorAll('#thbox input'),function(inp){var k=inp.getAttribute('data-t');if(k in thr)inp.value=thr[k];});
  recompute();renderAll();saveState();}
el('thbox').innerHTML='<label style="display:inline-flex;flex-direction:column;font-size:10px;color:#64748b;gap:2px">preset<select id="thpreset" class="msel" style="padding:4px 6px"><option value="">custom…</option>'+
  PRESETS.map(function(p){return '<option value="'+p.id+'">'+esc(p.label)+'</option>';}).join('')+'</select></label>'+
  '<span id="thnote" style="flex-basis:100%;font-size:10.5px;color:#8895a6;margin-top:-2px"></span>'+
  THL.map(function(t){return '<label style="display:inline-flex;flex-direction:column;font-size:10px;color:#64748b;gap:2px">'+t[1]+
  '<input type="number" step="any" data-t="'+t[0]+'" value="'+thr[t[0]]+'" style="width:78px;padding:4px 6px;border:1px solid var(--line);border-radius:6px;font-size:12px"></label>';}).join('')+
  '<button class="btn" id="threset" style="align-self:flex-end">reset</button>';
var thpre=el('thpreset'),thnote=el('thnote');
if(thpre)thpre.onchange=function(){var p=null;PRESETS.forEach(function(x){if(x.id==thpre.value)p=x;});if(!p){thnote.textContent='';return;}thnote.textContent=p.note;applyThr(p.th);};
Array.prototype.forEach.call(document.querySelectorAll('#thbox input'),function(inp){inp.oninput=function(){var v=parseFloat(inp.value);if(!isNaN(v)){thr[inp.getAttribute('data-t')]=v;if(thpre)thpre.value='';if(thnote)thnote.textContent='';clearTimeout(_thdb);_thdb=setTimeout(function(){recompute();renderAll();saveState();},180);}};});
el('threset').onclick=function(){if(thpre)thpre.value='';if(thnote)thnote.textContent='';applyThr(assign({},R.thresholds));};
// ancient (aDNA) live thresholds
if(R.n_ancient){var ATH=[['depth_min','aDNA depth min'],['breadth_min','aDNA breadth %'],['missing_max','aDNA missing %'],['mapping_min','aDNA mapped %'],['dup_max','aDNA dup %'],['iupac_max','aDNA IUPAC %'],['damage_min_ct',"5′ C>T min (0-1)"]];
  el('athbox').innerHTML='<div style="flex-basis:100%;font-size:10px;color:#8a5a12;font-weight:600;text-transform:uppercase;letter-spacing:.06em">Ancient (aDNA) thresholds &middot; a 5x mummy is judged here, not against the modern gate</div>'+
    ATH.map(function(t){return '<label style="display:inline-flex;flex-direction:column;font-size:10px;color:#64748b;gap:2px">'+t[1]+
    '<input type="number" step="any" data-t="'+t[0]+'" value="'+(athr[t[0]]!=null?athr[t[0]]:'')+'" style="width:78px;padding:4px 6px;border:1px solid var(--line);border-radius:6px;font-size:12px"></label>';}).join('');
  Array.prototype.forEach.call(document.querySelectorAll('#athbox input'),function(inp){inp.oninput=function(){var v=parseFloat(inp.value);if(!isNaN(v)){athr[inp.getAttribute('data-t')]=v;clearTimeout(_thdb);_thdb=setTimeout(function(){recompute();renderAll();saveState();},180);}};});}
// CSV export of the current (filtered rows, visible columns) table
el('csv').onclick=function(){var mets=R.metrics.filter(function(m){return !st.hidden[m.key];});
  var head=['sample','verdict'].concat(mets.map(function(m){return m.key;})).concat(['lineage','flags']);
  var lines=[head.join('\t')]; visible().filter(colMatch).forEach(function(s){lines.push([s.s,s.v].concat(mets.map(function(m){return s.m[m.key]==null?'':s.m[m.key];})).concat([s.lineage||'',s.f.join(';')]).join('\t'));});
  var blob=new Blob([lines.join('\n')],{type:'text/tab-separated-values'}),a=document.createElement('a');
  a.href=URL.createObjectURL(blob);a.download='qc_table.tsv';a.click();URL.revokeObjectURL(a.href);toast('Saved qc_table.tsv ('+visible().filter(colMatch).length+' rows)');};
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
// aDNA panel: only when there are ancient samples (no placeholder otherwise)
if(!R.n_ancient){var adx=el('adna');if(adx)adx.style.display='none';var nadx=el('nav-adna');if(nadx)nadx.style.display='none';}
// SNP dynamics: only when the metadata gave connected time-series (else no section at all)
if(!(R.dynamics&&R.dynamics.groups&&R.dynamics.groups.length)){var dyx=el('dynamics');if(dyx)dyx.style.display='none';var ndyx=el('nav-dyn');if(ndyx)ndyx.style.display='none';}
if(!(R.epistasis&&R.epistasis.pairs&&R.epistasis.pairs.length)){var epx=el('epistasis');if(epx)epx.style.display='none';var nepx=el('nav-epi');if(nepx)nepx.style.display='none';}
if(!(R.dr&&R.dr.calls&&R.dr.calls.length)){var drx=el('drug');if(drx)drx.style.display='none';var ndrx=el('nav-drug');if(ndrx)ndrx.style.display='none';}
if(!(R.kraken&&R.kraken.samples&&R.kraken.samples.length)){var kkx=el('kraken');if(kkx)kkx.style.display='none';var nkkx=el('nav-kraken');if(nkkx)nkkx.style.display='none';}
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
var colfCb=el('colf'); if(colfCb)colfCb.onchange=function(){st.showColF=colfCb.checked;renderTable();};
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
var SKEY='bampiro_qc_v2:'+R.samples.length+':'+(R.samples[0]?R.samples[0].s:'')+':'+(R.samples.length?R.samples[R.samples.length-1].s:'');
function saveState(){try{localStorage.setItem(SKEY,JSON.stringify({excl:st.excl,thr:thr,athr:athr,hidden:st.hidden,sortKey:st.sortKey,asc:st.asc}));}catch(e){}}
function loadState(){try{var s=JSON.parse(localStorage.getItem(SKEY)||'null');if(!s)return false;
  if(s.thr)Object.keys(s.thr).forEach(function(k){if(k in thr)thr[k]=s.thr[k];});
  if(s.athr)Object.keys(s.athr).forEach(function(k){if(k in athr)athr[k]=s.athr[k];});
  if(s.excl&&typeof s.excl=='object'){var have={};R.samples.forEach(function(x){have[x.s]=1;});
    st.excl={};Object.keys(s.excl).forEach(function(k){if(have[k])st.excl[k]=1;});}  // drop unknown sample ids
  if(s.hidden&&typeof s.hidden=='object')st.hidden=s.hidden;        // restore the chosen visible columns
  if(s.sortKey){st.sortKey=s.sortKey;st.asc=!!s.asc;}              // and the sort order (SKEY is per-cohort)
  return true;}catch(e){return false;}}
// ---- expand-to-fill (fullscreen within the window) for the big panels ----
function collapseExpanded(){var ex=document.querySelector('.panel.expanded');if(!ex)return;ex.classList.remove('expanded');document.body.classList.remove('has-expanded');
  Array.prototype.forEach.call(document.querySelectorAll('.exp-h'),function(b){b.innerHTML=icon('maximize')+'full';});
  renderGenome();renderPlots();renderScatter();renderTable();renderQCspace();renderRefBias();renderFunction();renderGeneBurden();renderHotspots();renderPnps();}
Array.prototype.forEach.call(document.querySelectorAll('.exp-h'),function(b){b.onclick=function(){
  var panel=el(b.getAttribute('data-panel')); if(!panel)return;
  var willExpand=!panel.classList.contains('expanded'); collapseExpanded();
  if(willExpand){panel.classList.add('expanded');document.body.classList.add('has-expanded');b.innerHTML=icon('minimize')+'close';}
  var rn=b.getAttribute('data-render');
  setTimeout(function(){if(rn=='genome')renderGenome();else if(rn=='plots')renderPlots();else if(rn=='scatter')renderScatter();else if(rn=='corr')renderCorr();else if(rn=='pca')renderQCspace();else if(rn=='divcomp')renderRefBias();else if(rn=='function')renderFunction();else if(rn=='geneburden')renderGeneBurden();else if(rn=='hotspots')renderHotspots();else if(rn=='pnps')renderPnps();else if(rn=='table')renderTable();},20);};});
if(el('expClose'))el('expClose').onclick=collapseExpanded;
document.addEventListener('keydown',function(e){if(e.key=='Escape')collapseExpanded();});
// per-sample detail modal close (button, backdrop, Esc)
el('modalx').onclick=closeDetail; el('modal').addEventListener('click',function(e){if(e.target==el('modal'))closeDetail();});
el('modal').addEventListener('keydown',function(e){   // trap Tab focus inside the open dialog
  if(e.key!=='Tab')return;
  var f=Array.prototype.filter.call(el('modal').querySelectorAll('a[href],button:not([disabled]),input:not([disabled]),select,textarea,[tabindex]:not([tabindex="-1"])'),function(x){return x.offsetParent!==null;});
  if(!f.length)return; var first=f[0],last=f[f.length-1];
  if(e.shiftKey&&document.activeElement===first){e.preventDefault();last.focus();}
  else if(!e.shiftKey&&document.activeElement===last){e.preventDefault();first.focus();}
});
document.addEventListener('keydown',function(e){if(e.key=='Escape'&&st.detail)closeDetail();});

var hadSaved=loadState();
Array.prototype.forEach.call(document.querySelectorAll('#thbox input'),function(inp){var k=inp.getAttribute('data-t');if(k in thr)inp.value=thr[k];});
if(R.n_ancient)Array.prototype.forEach.call(document.querySelectorAll('#athbox input'),function(inp){var k=inp.getAttribute('data-t');if(k in athr)inp.value=athr[k];});
Array.prototype.forEach.call(document.querySelectorAll('#colmenu input'),function(cb){cb.checked=!st.hidden[cb.getAttribute('data-k')];});  // sync the column checkboxes to any restored/hidden set
recompute();
if(!(hadSaved&&Object.keys(st.excl).length))R.samples.forEach(function(s){if(s.v=='FAIL')st.excl[s.s]=1;});  // preselect FAILs unless a saved basket exists
renderAll();
(function(){  // dark / light theme toggle (the early head script set the initial class from the saved pref or OS)
  var tb=el('themeToggle'); if(!tb)return;
  function setIcon(){tb.innerHTML=icon(isDark()?'sun':'moon');}
  setIcon();
  tb.onclick=function(){var d=!isDark();document.documentElement.classList.toggle('dark',d);
    try{localStorage.setItem('bampiro_theme',d?'dark':'light');}catch(e){}
    TH=d?TH_DARK:TH_LIGHT; setIcon(); renderAll();};
})();
(function(){  // left contents sidebar: collapse toggle, collapsible groups, scroll-spy highlight
  var toc=el('toc'), tg=el('toc-toggle'); if(!toc||!tg)return;
  tg.onclick=function(){ document.body.classList.toggle('toc-collapsed'); };
  var scrim=el('tocscrim'); if(scrim)scrim.onclick=function(){ document.body.classList.add('toc-collapsed'); };  // tap outside the drawer to dismiss (mobile)
  if(window.innerWidth&&window.innerWidth<860) document.body.classList.add('toc-collapsed');   // start collapsed on small screens
  Array.prototype.forEach.call(toc.querySelectorAll('.toc-gh'),function(gh){ gh.onclick=function(){ gh.parentNode.classList.toggle('closed'); }; });
  Array.prototype.forEach.call(toc.querySelectorAll('.toc-group'),function(g){   // hide a whole group if every section in it was self-hidden
    var any=false; Array.prototype.forEach.call(g.querySelectorAll('.toc-link'),function(a){ if(a.style.display!=='none') any=true; });
    if(!any) g.style.display='none';
  });
  var links=Array.prototype.slice.call(toc.querySelectorAll('.toc-link'));
  var items=links.map(function(a){ return {a:a, el:el(a.getAttribute('href').slice(1))}; }).filter(function(x){return x.el;});
  links.forEach(function(a){ a.onclick=function(){ if(window.innerWidth&&window.innerWidth<860) document.body.classList.add('toc-collapsed'); }; });
  function spy(){
    var se=document.scrollingElement||document.documentElement, y=se.scrollTop+92, cur=null;
    var best=-1;   // pick the section with the GREATEST offsetTop<=y (TOC order need not match DOM order)
    items.forEach(function(s){ if(s.el.style.display!=='none' && s.el.offsetTop<=y && s.el.offsetTop>=best){ best=s.el.offsetTop; cur=s; } });
    if(!cur){ for(var i=0;i<items.length;i++){ if(items[i].el.style.display!=='none'){ cur=items[i]; break; } } }  // above the 1st section -> highlight it
    links.forEach(function(a){ a.className='toc-link'; });
    if(cur) cur.a.className='toc-link active';
  }
  window.addEventListener('scroll',spy);
  window.addEventListener('resize',spy);
  spy();
})();
(function(){  // drop a clickable (i) into every section heading; the explanation text comes from R.section_info
  var info=R.section_info||{};
  Array.prototype.forEach.call(document.querySelectorAll('section[id]'),function(sec){
    var txt=info[sec.id]; if(!txt) return;
    var h2=sec.querySelector('h2'); if(!h2||h2.querySelector('.sec-infoi')) return;
    var ic=document.createElement('span');
    ic.className='infoi sec-infoi'; ic.textContent='i';
    ic.setAttribute('data-info',txt);
    ic.setAttribute('role','button'); ic.setAttribute('aria-label','About this analysis');
    var cap=h2.querySelector('.c');
    if(cap) h2.insertBefore(ic,cap); else h2.appendChild(ic);
  });
})();
})();
"""

# BAMpiro header logo (bampiro2.png resized to ~90x100 and embedded so the report stays
# self-contained - no external image request). Regenerate from .github/bampiro2.png if the logo changes.
LOGO_DATA_URI = "data:image/png;base64,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"

SHELL = """<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>__TITLE__</title>
<script>(function(){try{var t=localStorage.getItem('bampiro_theme');if(t=='dark'||(!t&&window.matchMedia&&matchMedia('(prefers-color-scheme:dark)').matches))document.documentElement.className+=' dark';}catch(e){}})();</script>
<style>__CSS__</style></head><body>
<button id="toc-toggle" title="Show / hide the contents sidebar" aria-label="Toggle contents"><span data-ic="panel"></span></button>
<div id="tocscrim" aria-hidden="true"></div>
<nav id="toc" aria-label="Contents">
<div class="toc-brand"><img class="brandlogo" src="__LOGO__" alt="BAMpiro logo"><b>BAMpiro</b> QC</div>
<div class="toc-group"><div class="toc-gh">Overview<span class="toc-chev" data-ic="chevronDown"></span></div><div class="toc-items"><a class="toc-link" href="#exec">Summary</a><a class="toc-link" href="#gstats">Stats</a><a class="toc-link" href="#flagged">Flagged</a><a class="toc-link" href="#linsum" id="nav-lin">Lineages</a><a class="toc-link" href="#dist">Distributions</a><a class="toc-link" href="#kraken" id="nav-kraken">Kraken</a></div></div>
<div class="toc-group"><div class="toc-gh">Correlation &amp; structure<span class="toc-chev" data-ic="chevronDown"></span></div><div class="toc-items"><a class="toc-link" href="#corr">Correlations</a><a class="toc-link" href="#corrmatrix">Corr matrix</a><a class="toc-link" href="#qcpca" id="nav-pca">QC space</a><a class="toc-link" href="#divcomp" id="nav-divcomp">Divergence</a></div></div>
<div class="toc-group"><div class="toc-gh">Genome &amp; genes<span class="toc-chev" data-ic="chevronDown"></span></div><div class="toc-items"><a class="toc-link" href="#cons">Consensus</a><a class="toc-link" href="#genome" id="nav-genome">Genome</a><a class="toc-link" href="#function" id="nav-function">Function</a><a class="toc-link" href="#geneburden" id="nav-geneburden">Gene burden</a><a class="toc-link" href="#hotspots" id="nav-hot">Variable genes</a></div></div>
<div class="toc-group"><div class="toc-gh">Evolution<span class="toc-chev" data-ic="chevronDown"></span></div><div class="toc-items"><a class="toc-link" href="#temporal" id="nav-temporal">Temporal</a><a class="toc-link" href="#pnps" id="nav-pnps">pN/pS</a><a class="toc-link" href="#adna" id="nav-adna">aDNA</a></div></div>
<div class="toc-group"><div class="toc-gh">Variants over time<span class="toc-chev" data-ic="chevronDown"></span></div><div class="toc-items"><a class="toc-link" href="#dynamics" id="nav-dyn">SNP dynamics</a><a class="toc-link" href="#epistasis" id="nav-epi">Epistasis</a><a class="toc-link" href="#snpmatrix" id="nav-snpmx">SNP matrix</a><a class="toc-link" href="#drug" id="nav-drug">Drug resistance</a></div></div>
</nav>
<header><span class="logo"><img class="brandlogo" src="__LOGO__" alt="BAMpiro logo"><b>BAMpiro</b> QC</span><span class="ver" id="hver"></span><span class="meta" id="meta"></span><a id="ghlink" class="hbtn" href="__REPO__" target="_blank" rel="noopener noreferrer" title="BAMpiro source on GitHub" aria-label="BAMpiro source on GitHub" style="margin-left:auto"><span data-ic="github"></span></a><button id="themeToggle" class="hbtn" title="Toggle dark / light theme" aria-label="Toggle dark / light theme"></button></header>
<div class="wrap">
<p class="lede">Short-read bacterial / MTBC cohort QC. Review <a href="#flagged">flagged samples</a>, tick any to drop, then export <b>keep_list.txt</b> / <b>exclusion.tsv</b>. Thresholds below are live; the pipeline gate itself is unchanged.</p>
<section class="hero"><div class="summary" id="summary"></div><div class="chips" id="chips"></div></section>
<section id="exec"><h2>Executive summary <span class="c">- cohort health and headline findings at a glance</span></h2>
<div class="panel pad" id="exec_body"></div></section>
<div class="provbar"><div class="prov" id="prov"></div><button class="btn" id="printBtn" title="expand + print / save as PDF"><span data-ic="printer"></span> print</button></div>
<section><details class="dd" style="display:inline-block"><summary><span data-ic="sliders"></span> Live thresholds &amp; presets - adjust and everything re-flags<span class="ddcaret" data-ic="chevronDown" data-ic-cls="sort"></span></summary>
<div class="menu" style="display:flex;gap:12px;flex-wrap:wrap;align-items:flex-end;min-width:min(520px,calc(100vw - 28px))">
  <div style="display:flex;gap:12px;flex-wrap:wrap;align-items:flex-end;flex-basis:100%" id="thbox"></div>
  <div style="display:flex;gap:12px;flex-wrap:wrap;align-items:flex-end;flex-basis:100%" id="athbox"></div>
</div></details></section>
<section id="gstats"><h2>General statistics <span class="c">- tick a box to basket a sample for exclusion; click a sample name for its full profile; a header to sort</span><button class="exp-h" data-panel="gstatsPanel" data-render="table" style="margin-left:auto"><span data-ic="maximize"></span>full</button></h2>
<div class="controls">
  <input id="q" type="search" placeholder="filter samples…">
  <label><input id="of" type="checkbox"> only flagged</label>
  <label id="groupui"><input id="glin" type="checkbox"> group by lineage</label>
  <label title="show a search / filter box under every column header"><input id="colf" type="checkbox"> column filters</label>
  <span id="ancfilter"></span>
  <details class="dd"><summary>columns <span data-ic="chevronDown" data-ic-cls="sort"></span></summary><div class="menu" id="colmenu"></div></details>
  <details class="dd"><summary><span data-ic="help"></span> metric help</summary><div class="menu" id="helpmenu"></div></details>
  <button class="btn" id="csv"><span data-ic="download"></span> export TSV</button>
  <a class="hint" id="nbasket" href="#curation" title="samples in the exclusion basket - click to jump to it" style="text-decoration:none;color:var(--accent);cursor:pointer"></a>
  <span class="hint"><span id="nshown"></span></span>
</div>
<div class="panel gtable" id="gstatsPanel"><table id="gstable"></table></div></section>
<section id="flagged"><h2>Flagged samples <span class="c">- <span id="nflag"></span> to review; the margin is shown under each flag; click a flag to filter the cohort</span>
<button class="btn" id="basketFlagged" style="margin-left:auto"><span data-ic="basket"></span> basket all flagged</button></h2>
<div class="panel gtable" style="max-height:50vh"><table id="flagtable"></table></div></section>
<section id="curation-sec"><div class="curation" id="curation"></div></section>
<section id="linsum"><h2>Per-lineage summary <span class="c">- medians per lineage; # MIXED = samples with &gt;1 lineage above the mixture cut-off; click a row to filter</span></h2>
<div class="panel"><div id="lincomp"></div><div class="gtable"><table id="linsumtable"></table></div></div></section>
<section id="dist"><h2>Distributions <span class="c">- one mark per sample; shaded band = acceptable range (modern gate); hover for detail</span>
<span class="seg" id="colorby" style="margin-left:auto"><button class="on" data-cb="qc">colour: QC</button><button data-cb="lineage">lineage</button></span>
<span class="seg" id="ptype"><button class="on" data-t="beeswarm">beeswarm</button><button data-t="bar">bar</button><button data-t="histogram">histogram</button></span>
<button class="exp-h" data-panel="plotsPanel" data-render="plots"><span data-ic="maximize"></span>full</button></h2>
<div class="panel" id="plotsPanel"><div id="plots" class="pad" style="padding-top:6px;padding-bottom:6px"></div>
<div class="legend" id="leg-qc"><span><i style="background:#94a3b8"></i>PASS</span><span><i style="background:#d97706"></i>WARN</span><span><i style="background:#dc2626"></i>FAIL</span><span style="margin-left:auto">beeswarm dashed line = median</span></div>
<div class="legend" id="leg-lin" style="display:none"></div></div></section>
<section id="kraken"><h2>Taxonomic composition <span class="c">- Kraken2 classification of the reads before mapping: the dominant (primary) taxon, the top other taxon, and the unclassified fraction. A low primary % or a large secondary taxon flags possible contamination / a mixed sample.</span></h2>
<div class="panel pad" id="kraken_body"></div></section>
<section id="corr"><h2>Correlations <span class="c">- pick two metrics, coloured by QC; drag a box to basket the enclosed samples</span>
<span id="scbrushinfo" style="margin-left:auto;font-size:11px;color:var(--accent)"></span>
<select class="msel" id="sx"></select><span style="color:#94a3b8">vs</span><select class="msel" id="sy"></select>
<button class="exp-h" data-panel="scatter" data-render="scatter"><span data-ic="maximize"></span>full</button></h2>
<div class="panel pad" id="scatter"></div></section>
<section id="corrmatrix"><h2>Metric correlation <span class="c">- Spearman rank correlation across the core metrics, over the samples in view; click a cell to load that pair into the scatter above</span>
<button class="exp-h" data-panel="corrmatrix" data-render="corr" style="margin-left:auto"><span data-ic="maximize"></span>full</button></h2>
<div class="panel pad" id="corr_body"></div></section>
<section id="qcpca"><h2>QC-metric space <span class="c">- samples ordinated in standardized QC-metric space (PCA on robust z-scores); a descriptive map of how quality profiles co-vary, not phylogeny and not a test</span>
<span class="seg" id="pcacb" style="margin-left:auto"><button class="on" data-cb="qc">colour: QC</button><button data-cb="lineage">lineage</button></span>
<button class="exp-h" data-panel="qcpcaPanel" data-render="pca"><span data-ic="maximize"></span>full</button></h2>
<div class="panel" id="qcpcaPanel"><div id="qcpca_body" class="pad"></div>
<div class="hot-note" id="qcpca_caption"></div>
<div class="dsub" style="margin:14px 16px 4px">Most unusual samples <span style="font-weight:400;color:#94a3b8">(robust Mahalanobis rank + top deviating metrics)</span></div>
<div class="gtable" style="max-height:34vh"><table id="pca_outtable"></table></div></div></section>
<section id="divcomp"><h2>Divergence vs completeness <span class="c">- separates reference-bias suspects (low divergence at low missingness) from honest low-coverage (low divergence explained by high missingness); a screen to inspect, not an automatic flag</span>
<span class="seg" id="divx" style="margin-left:auto"><button class="on" data-x="missing_pct">x: missing %</button><button data-x="callable_inv">100 - callable %</button></span>
<button class="exp-h" data-panel="divcompPanel" data-render="divcomp"><span data-ic="maximize"></span>full</button></h2>
<div class="panel" id="divcompPanel"><div id="divcomp_body" class="pad"></div>
<div class="legend" id="divcomp_quadn"></div>
<div class="hot-note" id="divcomp_caption"></div></div></section>
<section id="cons"><h2>Consensus completeness <span class="c">- callable / IUPAC / missing per sample, worst first</span></h2>
<div class="panel"><div id="stacks" style="max-height:60vh;overflow:auto"></div>
<div class="legend"><span><i style="background:#22a06b"></i>callable</span><span><i style="background:#e6b25a"></i>IUPAC</span><span><i style="background:#cbd5e1"></i>missing (- / N)</span></div></div></section>
<section id="genome"><h2>Genome landscape <span class="c">- a signal along the reference; rows = samples (grouped by lineage). A contiguous block = a localised feature (e.g. an RD deletion, a variant cluster). Missing = consensus callability (coverage proxy); SNPs/Het/Indels = VCF variant density. Drag across the track to zoom a region</span>
<input class="gsearch" id="genegoto" type="search" placeholder="go to gene" style="margin-left:auto"><span id="gselreadout" style="font-size:11px;color:var(--accent)"></span><button class="btn" id="maskbtn"><span data-ic="ban"></span>mask regions</button><span class="seg" id="gtrack"></span>
<button class="exp-h" data-panel="genomePanel" data-render="genome"><span data-ic="maximize"></span>full</button></h2>
<div class="panel" id="genomePanel"><div class="genome-scroll"><div id="genome_body" class="pad"></div></div>
<div class="legend"><span><i style="background:rgb(238,244,240)"></i>none / callable</span><span><i style="background:rgb(160,120,150)"></i>&#8594;</span><span><i style="background:rgb(120,70,90)"></i>high</span><span style="margin-left:auto">top strip = cohort mean; hover a cell for the position + value</span></div></div></section>
<section id="function"><h2>Functional annotation <span class="c">- snpEff impact + effect classes per sample; the missense/silent ratio is a pN/pS PROXY, not a selection test</span>
<button class="exp-h" data-panel="functionPanel" data-render="function" style="margin-left:auto"><span data-ic="maximize"></span>full</button></h2>
<div class="panel" id="functionPanel">
  <div class="provbar" style="margin:12px 16px 4px"><div class="prov" id="fn_cohort"></div></div>
  <div id="fn_stacks" style="max-height:56vh;overflow:auto"></div>
  <div class="legend" id="fn_legend"></div>
  <div class="hot-note" id="fn_caption"></div>
</div></section>
<section id="geneburden"><h2>Functional gene burden <span class="c">- genes carrying the most HIGH+MODERATE variants across the cohort, with the dominant effect class; a functional companion to Variable genes below, not a selection test</span>
<input class="gsearch" id="gbq" type="search" placeholder="search gene" style="margin-left:auto"><button class="exp-h" data-panel="geneburdenPanel" data-render="geneburden"><span data-ic="maximize"></span>full</button></h2>
<div class="panel" id="geneburdenPanel"><div id="gb_body"><div class="hot-note" id="gb_note"></div><div class="gtable" style="max-height:44vh"><table id="gbtable"></table></div></div></div></section>
<section id="hotspots"><h2>Variable genes <span class="c">- genes/regions with the most SNPs across the cohort; click a row to mark it on the SNP track</span>
<input class="gsearch" id="hotq" type="search" placeholder="search gene" style="margin-left:auto"><button class="exp-h" data-panel="hotspotsPanel" data-render="hotspots"><span data-ic="maximize"></span>full</button></h2>
<div class="panel" id="hotspotsPanel"><div id="hot_body"><div class="hot-note" id="hot_note"></div><div class="gtable" style="max-height:44vh"><table id="hottable"></table></div></div></div></section>
<section id="dynamics"><h2>SNP dynamics <span class="c">- search &amp; select genes to see the allele-frequency trajectories of their variants over time; events flag emergence / fixation / loss / non-synonymous</span></h2>
<div class="panel pad" id="dyn_body"></div></section>
<section id="epistasis"><h2>Epistasis <span class="c">- pairs of variants whose allele-frequency trajectories move together (concordant) or in opposition (discordant) within a connected series; candidate linked / co-selected / competing SNPs</span></h2>
<div class="panel pad" id="epi_body"></div></section>
<section id="snpmatrix"><h2>SNP matrix <span class="c">- every SNP site (rows) &#215; sample (columns); each cell is the allele frequency (hover for AF &amp; depth), a striped cell = not called. Filter by gene / position, then download the full matrix as a TSV.</span></h2>
<div class="panel pad" id="snpmx_body"></div></section>
<section id="drug"><h2>Drug resistance <span class="c">- resistance-associated mutations detected by pathotypr (WHO catalogue, H37Rv numbering). A genomic screen, NOT a clinical result; grade 1&#8211;2 = associated with resistance, 3 = uncertain, 4&#8211;5 = not associated.</span></h2>
<div class="panel pad" id="drug_body"></div></section>
<section id="adna"><h2>aDNA damage authentication <span class="c">- terminal deamination per ancient sample; a screen, not a proof of authenticity</span></h2>
<div class="panel"><div id="adna_body" class="pad"></div></div></section>
<section id="temporal"><h2>Temporal sampling overview <span class="c">- per-lineage year span parsed from dates + an informative-site proxy; a readiness check for downstream time-resolved analysis, NOT a clock estimate</span></h2>
<div class="panel"><div id="temporal_body" class="pad"></div>
<div class="legend"><span><i style="background:var(--accent)"></i>sample year</span><span>bar = sampling span</span><span style="margin-left:auto">informative-site proxy = median SNPs</span></div>
<div class="hot-note" id="temporal_caption"></div></div></section>
<section id="pnps"><h2>Selection: pN/pS and dN/dS <span class="c">- alignment-based per-gene dN/dS from eskaks; a cohort selection screen, not a per-sample QC metric</span>
<input class="gsearch" id="pnpsq" type="search" placeholder="search gene" style="margin-left:auto"><button class="exp-h" data-panel="pnpsPanel" data-render="pnps"><span data-ic="maximize"></span>full</button></h2>
<div class="panel" id="pnpsPanel"><div id="pnps_body"><div class="hot-note" id="pnps_note"></div><div class="gtable" style="max-height:48vh"><table id="pnpstable"></table></div></div></div></section>
<div class="footer" id="foot"></div>
</div>
<button id="expClose" class="exp-close" aria-label="exit fullscreen"><span data-ic="x"></span>close (Esc)</button>
<div id="modal" class="modal"><div class="modalcard" role="dialog" aria-modal="true" aria-label="Sample detail"><button class="modalx" id="modalx" aria-label="close"><span data-ic="x"></span></button><div id="modalbody"></div></div></div>
<div id="tt"></div><div id="toast" role="status" aria-live="polite"></div>
<div id="infopop"></div>
<script>var REPORT=__JSON__;</script>
<script>__JS__</script>
</body></html>"""


def build_html(title, payload):
    data = json.dumps(payload, separators=(",", ":")).replace("</", "<\\/")
    # Inject content first, then substitute __TITLE__ LAST so a user --title that happens to
    # contain a placeholder token (e.g. "__JS__") can never pull in the CSS/JS/logo/JSON blob.
    return (SHELL.replace("__CSS__", CSS).replace("__LOGO__", LOGO_DATA_URI)
                 .replace("__JS__", JS).replace("__JSON__", data)
                 .replace("__REPO__", html.escape(REPO_URL))
                 .replace("__TITLE__", html.escape(title)))


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
    fields = [(i, h) for i, h in enumerate(header)
              if i != si and h.strip() and h.strip().lower() not in _META_SKIP]
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
    return {'fields': [h for _, h in fields], 'rows': out}


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
                                                     'af': round(af, 4), 'dp': dp}
        except OSError:
            continue
    return out


def build_dynamics(metadata, variants, sample_meta=None, emerge=0.25, fix=0.90, loss=0.10, min_points=2, min_move=0.15):
    """Per-connected-group AF trajectories over time, only for SNPs that move. None if nothing to show.
    sample_meta (parse_sample_meta output) attaches each group's metadata signature so the report can
    filter the series by any group-invariant field (e.g. patient / lineage)."""
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
        for pos in sorted(allpos):   # deterministic order (a set iterates in a PYTHONHASHSEED-dependent order)
            traj, dps, meta = [], [], None
            for s in samples:
                v = variants[s].get(pos)
                traj.append(v['af'] if v else 0.0)
                dps.append(v.get('dp') if v else None)   # per-timepoint depth (None where not called)
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
                           'aa': (meta or {}).get('aa', ''), 'aa_h37rv': (meta or {}).get('aa_h37rv', ''),
                           'traj': [round(x, 4) for x in traj], 'dp': dps, 'flags': flags})
        if not series:
            continue
        series.sort(key=lambda x: (-(max(x['traj']) - min(x['traj'])), -len(x['flags'])))
        gmeta = {}
        if sample_meta and sample_meta.get('fields'):
            rows = sample_meta.get('rows', {})
            for f in sample_meta['fields']:
                vals = sorted({v for v in ((rows.get(s) or {}).get(f) for s in samples) if v})
                if vals:
                    gmeta[f] = vals
        out_groups.append({'group': g, 'samples': samples, 'times': times,
                           'tnums': [metadata[s]['tnum'] for s in samples], 'meta': gmeta,
                           'series': series, 'n_flagged': sum(1 for x in series if x['flags'])})
    if not out_groups:
        return None
    return {'groups': out_groups, 'thresholds': {'emerge': emerge, 'fix': fix, 'loss': loss}}


def _pearson(x, y):
    """Pearson correlation of two equal-length trajectories, or None if undefined (n<3 or a flat one)."""
    n = len(x)
    if n < 3 or len(y) != n:
        return None
    mx, my = sum(x) / n, sum(y) / n
    sxx = sum((a - mx) ** 2 for a in x)
    syy = sum((b - my) ** 2 for b in y)
    if sxx <= 1e-9 or syy <= 1e-9:
        return None
    sxy = sum((a - mx) * (b - my) for a, b in zip(x, y))
    return sxy / ((sxx * syy) ** 0.5)


def _perm_p(series_traj, obs_mean, B=2000, seed=0):
    """Permutation p-value for a pair's mean correlation: shuffle one trajectory's values within EACH
    series independently (breaking the time linkage) and recompute the mean |r|; p = fraction of the B
    nulls reaching the observed |r| (+1 smoothing). Deterministic (fixed per-pair seed) -> reproducible.
    With few timepoints the resolution is inherently coarse (a 4-point series has only 4!=24 orders), so
    a lone short series lands near p~0.04 while a pattern recurring across series drives p far lower."""
    obs = abs(obs_mean)
    rng = random.Random(seed)
    ge = 0
    for _ in range(B):
        tot, cnt = 0.0, 0
        for ta, tb in series_traj:
            perm = list(ta)
            rng.shuffle(perm)
            r = _pearson(perm, tb)
            if r is not None:
                tot += r
                cnt += 1
        if cnt and abs(tot / cnt) >= obs:
            ge += 1
    return (ge + 1) / (B + 1)


def build_epistasis(dynamics, min_r=0.8, min_points=3, top=300, perm=2000):
    """Candidate epistatic / linked variant pairs: two SNPs whose allele-frequency trajectories co-vary
    WITHIN a connected series - concordant (rise/fall together) or discordant (one rises as the other
    falls). Score = mean Pearson r of the two trajectories across every series where both move and the
    series has >= min_points timepoints. Significance = a permutation p-value (per pair), then a
    Benjamini-Hochberg FDR q-value across all reported pairs. Recurrence across independent series is
    tracked so a pattern seen in several series (low q) is separable from a lucky single short series.
    None if nothing qualifies. Built from the dynamics payload (moving variants per series)."""
    if not dynamics or not dynamics.get('groups'):
        return None
    pairs, n_series = {}, 0
    for g in dynamics['groups']:
        times = g.get('times') or []
        if len(set(times)) < min_points:
            continue
        n_series += 1
        series = g.get('series') or []
        for i in range(len(series)):
            for j in range(i + 1, len(series)):
                a, b = series[i], series[j]
                r = _pearson(a['traj'], b['traj'])
                if r is None:
                    continue
                lo, hi = (a, b) if str(a['pos']) <= str(b['pos']) else (b, a)  # order-independent key
                key = (lo['pos'], hi['pos'])
                rec = pairs.get(key)
                if rec is None:
                    rec = pairs[key] = {'A': lo, 'B': hi, 'rs': [], 'straj': [], 'bestn': -1,
                                        'times': times, 'trajA': lo['traj'], 'trajB': hi['traj'], 'group': g['group']}
                rec['rs'].append(r)
                rec['straj'].append((lo['traj'], hi['traj']))   # every series' pair, for the permutation test
                if len(set(times)) > rec['bestn']:   # keep the richest series for the mini-chart
                    rec['bestn'] = len(set(times))
                    rec['times'], rec['trajA'], rec['trajB'], rec['group'] = times, lo['traj'], hi['traj'], g['group']
    out = []
    for key, rec in pairs.items():
        rs = rec['rs']
        mean_r = sum(rs) / len(rs)
        if abs(mean_r) < min_r:
            continue
        A, B = rec['A'], rec['B']
        n_pos = sum(1 for r in rs if r > 0)
        # seed from the STABLE pair key (positions), not the iteration index, so the permutation p-value is
        # reproducible regardless of dict/set iteration order. crc32 is deterministic (unlike hash()).
        seed = zlib.crc32(('%s|%s' % (key[0], key[1])).encode('utf-8'))
        p = _perm_p(rec['straj'], mean_r, B=perm, seed=seed)
        out.append({'geneA': A.get('gene', ''), 'posA': A['pos'], 'aaA': A.get('aa', ''), 'aaA_h37rv': A.get('aa_h37rv', ''), 'effA': A.get('eff', ''),
                    'geneB': B.get('gene', ''), 'posB': B['pos'], 'aaB': B.get('aa', ''), 'aaB_h37rv': B.get('aa_h37rv', ''), 'effB': B.get('eff', ''),
                    'r': round(mean_r, 3), 'n': len(rs), 'rmin': round(min(rs), 3), 'rmax': round(max(rs), 3),
                    'direction': 'concordant' if mean_r > 0 else 'discordant',
                    'p': p, 'consistent': (n_pos == len(rs) or n_pos == 0),
                    'times': rec['times'], 'trajA': rec['trajA'], 'trajB': rec['trajB'], 'group': rec['group']})
    if not out:
        return None
    # Benjamini-Hochberg FDR across every reported pair (multiple-testing correction)
    m = len(out)
    order = sorted(range(m), key=lambda i: out[i]['p'])
    qmin = 1.0
    for rank in range(m, 0, -1):
        i = order[rank - 1]
        qmin = min(qmin, out[i]['p'] * m / rank)
        out[i]['q'] = round(min(qmin, 1.0), 4)
    for pr in out:
        pr['tier'] = 'strong' if pr['q'] <= 0.05 else ('moderate' if pr['p'] <= 0.05 else 'weak')
        pr['p'] = round(pr['p'], 4)
    out.sort(key=lambda p: (p['q'], -abs(p['r'])))
    out = out[:top]
    # all-vs-all correlation matrix over the most-connected variants (INCLUDES sub-threshold cells, so it
    # is the full co-dynamics landscape, not only the reported pairs). Capped to bound the embedded size.
    allmeta, deg, meanr = {}, {}, {}
    for (a, b), rec in pairs.items():
        meanr[(a, b)] = sum(rec['rs']) / len(rec['rs'])
        deg[a] = deg.get(a, 0) + 1
        deg[b] = deg.get(b, 0) + 1
        for v in (rec['A'], rec['B']):
            allmeta.setdefault(v['pos'], {'pos': v['pos'], 'gene': v.get('gene', ''), 'aa': v.get('aa', '')})
    MAXNODES = 50
    node_pos = sorted(deg, key=lambda x: (-deg[x], x))[:MAXNODES]          # keep the best-connected variants
    node_pos.sort(key=lambda x: (allmeta[x]['gene'] or 'zzz', x))         # display order: gene, then position
    nidx = {p: i for i, p in enumerate(node_pos)}
    cells = {}
    for (a, b), mr in meanr.items():
        if a in nidx and b in nidx:
            i, j = nidx[a], nidx[b]
            lo, hi = (i, j) if i <= j else (j, i)
            cells['%d,%d' % (lo, hi)] = round(mr, 3)
    matrix = {'nodes': [allmeta[p] for p in node_pos], 'cells': cells,
              'truncated': len(deg) > MAXNODES, 'total_nodes': len(deg)}
    return {'pairs': out, 'min_r': min_r, 'min_points': min_points, 'n_series': n_series, 'perm': perm,
            'matrix': matrix,
            'n_concordant': sum(1 for p in out if p['direction'] == 'concordant'),
            'n_discordant': sum(1 for p in out if p['direction'] == 'discordant'),
            'n_strong': sum(1 for p in out if p['tier'] == 'strong')}


def build_snp_matrix(variants, reference='', max_sites=8000):
    """Sparse SNP matrix for the report: samples + one row per SNP site with its ref/alt/annotation
    and only the cells (sample index -> [af, dp]) actually called. Capped to the most-shared sites to
    bound the embedded HTML size; the pipeline writes the complete matrix TSV separately."""
    if not variants:
        return None
    samples = list(variants.keys())
    sidx = {s: i for i, s in enumerate(samples)}
    sites = {}
    for s, posmap in variants.items():
        for key, v in posmap.items():
            contig, _, pos = key.rpartition(':')
            try:
                ipos = int(pos)
            except ValueError:
                continue
            st = sites.setdefault(key, {'contig': contig, 'pos': ipos, 'ref': v.get('ref', ''),
                                        'alt': set(), 'gene': '', 'eff': '', 'aa': '', 'aa_h37rv': '', 'cells': {}})
            if v.get('alt'):
                st['alt'].add(v['alt'])
            for k in ('gene', 'eff', 'aa', 'aa_h37rv'):
                if v.get(k) and not st[k]:
                    st[k] = v[k]
            st['cells'][sidx[s]] = [v['af'], v.get('dp')]
    ordered = sorted(sites.values(), key=lambda x: (x['contig'], x['pos']))
    truncated = len(ordered) > max_sites
    if truncated:
        ordered = sorted(sites.values(), key=lambda x: -len(x['cells']))[:max_sites]
        ordered.sort(key=lambda x: (x['contig'], x['pos']))
    rows = [{'contig': x['contig'], 'pos': x['pos'], 'ref': x['ref'], 'alt': ','.join(sorted(x['alt'])),
             'gene': x['gene'], 'eff': x['eff'], 'aa': x['aa'], 'aa_h37rv': x['aa_h37rv'], 'n': len(x['cells']), 'cells': x['cells']}
            for x in ordered]
    return {'samples': samples, 'reference': reference, 'rows': rows,
            'total_sites': len(sites), 'truncated': truncated}


# Per-panel explanation shown by the (i) icon in each section heading: what the analysis is, how to read
# it, and its caveats. Keyed by the section id in the HTML shell.
SECTION_INFO = {
    "gstats": "Per-sample sequencing and mapping QC: reads, duplication, mapped %, mean/median depth, "
              "breadth, variant counts and Ti/Tv. Each numeric cell is shaded by where the sample sits in "
              "that column's range; the QC column flags PASS/WARN/FAIL against the live thresholds. Tick a "
              "box to basket a sample for exclusion; click a name for its full profile.",
    "linsum": "Median QC metrics grouped by the assigned MTBC lineage, with the number of samples, %PASS and "
              "how many are 'mixed' (more than one lineage above the mixture cut-off). Lineage comes from "
              "pathotypr; a mixed call can indicate co-infection or contamination.",
    "dist": "Per-metric value distributions across the cohort (a strip/bee plot, or a ranked bar). The shaded "
            "band is the acceptable range for the current thresholds; points outside it are outliers. Useful "
            "to see cohort spread and choose sensible cut-offs.",
    "corr": "A scatter of any two QC metrics, one point per sample, coloured by QC status or lineage. Drag a "
            "box to basket the enclosed samples. Reveals metric relationships (e.g. depth vs breadth) and "
            "samples that are joint outliers.",
    "corrmatrix": "Pairwise Pearson correlation between the QC metrics across the cohort, as a heatmap. Shows "
                  "which metrics are redundant (strongly correlated) and which capture independent axes of quality.",
    "qcpca": "PCA of the standardised QC metrics: each point is a sample placed by its overall quality "
             "profile. Samples that cluster share a profile; a lone point is a multivariate outlier that no "
             "single metric flags.",
    "divcomp": "Each sample's genomic divergence (SNPs vs the reference) against its consensus completeness. "
               "Separates genuinely divergent samples from those that only look divergent because they are "
               "incomplete or contaminated.",
    "cons": "How much of the reference each sample's consensus recovers: callable %, missing % (N/gap), IUPAC "
            "(ambiguous) % and the longest gap. A phylogeny-oriented view of completeness, not just average depth.",
    "genome": "Per-position callability / variant density along the reference, binned into a heatmap per "
              "sample. Brush a region to list the genes under it. Reveals systematically low-callability "
              "regions (repeats, deletions) shared across samples.",
    "function": "The snpEff functional class of each sample's variants (HIGH/MODERATE/LOW/MODIFIER impact and "
                "effect types such as missense / synonymous). A per-sample mutational-impact profile.",
    "geneburden": "Genes carrying the most impactful (HIGH/MODERATE) variants across the cohort, with the "
                  "dominant effect and how many samples are hit. A cohort-level view of which genes accumulate "
                  "functional change.",
    "hotspots": "Genes with the highest SNP density across the cohort (variants per gene length). Flags the "
                "most variable loci - often repeats, antigens, or genes under diversifying selection.",
    "temporal": "The distribution of sampling dates, when the samplesheet provides them. Context for temporal "
                "analyses (e.g. SNP dynamics); not a QC metric itself.",
    "pnps": "Cohort-level pairwise dN/dS per gene (eskaks, Nei/Li). dN/dS > 1 suggests positive / diversifying "
            "selection at cohort scale. This is population genetics, NOT per-sample QC: noisy per gene, "
            "sensitive to alignment, and saturating at low divergence (NA when dS is near zero).",
    "adna": "For samples marked ancient: the characteristic post-mortem damage signal (elevated terminal 5' "
            "C to T). A sample below the damage floor may be a modern contaminant rather than authentic "
            "ancient DNA.",
    "dynamics": "For connected time-series (patient / passage line), each variant's allele-frequency "
                "trajectory over time, with per-timepoint read depth. Events flag emergence, fixation, loss "
                "and non-synonymous changes. Needs a time column and a group column in the samplesheet.",
    "epistasis": "Pairs of variants whose allele-frequency trajectories co-vary within a series: concordant "
                 "(rise/fall together) or discordant (one rises as the other falls). Scored by the Pearson "
                 "correlation, a permutation p-value and a Benjamini-Hochberg FDR q; a pattern recurring across "
                 "independent series is what makes a pair 'strong'. Candidate linked / co-selected / competing "
                 "SNPs, NOT proof of a functional interaction.",
    "snpmatrix": "Every SNP site (rows) by sample (columns); each cell is the allele frequency with its depth, "
                 "plus a reference column and samplesheet metadata as column-header levels. Filter by gene / "
                 "position or by metadata, and download the full matrix as a TSV.",
    "drug": "Resistance-associated mutations detected by pathotypr against the WHO catalogue (H37Rv numbering), "
            "as a sample x drug matrix (worst grade per drug) and a per-mutation table with the WHO confidence "
            "grade. Alignment-free (k-mer), so it works regardless of the mapping reference. A genomic screen, "
            "NOT a clinical drug-susceptibility result.",
    "flagged": "The samples the current thresholds flag as WARN / FAIL and the specific reason for each. This "
               "is the actionable QC summary; adjust the thresholds above to re-flag the whole report.",
}


def main():
    ap = argparse.ArgumentParser(description="BAMpiro interactive QC report + flags.")
    ap.add_argument("--summary", required=True)
    ap.add_argument("--consensus", nargs="*", default=[])
    ap.add_argument("--out-html", required=True)
    ap.add_argument("--out-flags", required=True)
    ap.add_argument("--title", default="BAMpiro QC report")
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
    ap.add_argument("--dr-report", default=None,
                    help="Per-sample drug-resistance calls TSV (from pathotypr DR markers) -> the Drug resistance panel (optional).")
    ap.add_argument("--kraken", nargs="*", default=[],
                    help="Per-sample Kraken2 .report files -> the Taxonomic composition / contamination panel (optional).")
    ap.add_argument("--pnps", default=None,
                    help="Cohort per-gene dN/dS TSV (eskaks) -> the Selection pN/pS panel (optional; a cohort analysis, not per-sample QC).")
    ap.add_argument("--provenance", nargs="*", default=[],
                    help="key=value provenance pairs surfaced in the report header (reference, container, commit...).")
    ap.add_argument("--version", default="",
                    help="BAMpiro version string (workflow.manifest.version) shown in the header + footer.")
    ap.add_argument("--metadata", default=None,
                    help="Optional TSV (e.g. the samplesheet) with a sample column + time (passage/timepoint) "
                         "and group (patient/series/cluster) columns -> the SNP dynamics panel. Auto-detected.")
    ap.add_argument("--vcfs", nargs="*", default=[],
                    help="Optional per-sample annotated VCFs -> per-SNP allele frequencies for the dynamics panel.")
    ap.add_argument("--vcfs-h37rv", nargs="*", default=[],
                    help="Optional per-sample VCFs annotated against a canonical reference -> the canonical "
                         "amino-acid position shown alongside the used-reference one (matched by sample + contig:pos).")
    ap.add_argument("--aa2-label", default="H37Rv",
                    help="Label for the canonical-reference amino-acid numbering shown by --vcfs-h37rv (default H37Rv).")
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

    _variants = parse_vcfs(args.vcfs)   # parsed once, feeds both the dynamics panel and the SNP matrix
    if args.vcfs_h37rv:   # attach the H37Rv/Mycobrowser amino-acid change per variant
        _h37 = parse_vcfs(args.vcfs_h37rv)
        for _s, _pm in _variants.items():
            _hs = _h37.get(_s, {})
            # also index by bare position: MTB is single-contig, so a differing contig NAME between the
            # used-reference and H37Rv VCFs should still match (only the coordinate needs to line up).
            _hpos = {_k.rpartition(':')[2]: _hv for _k, _hv in _hs.items()}
            for _key, _v in _pm.items():
                _hv = _hs.get(_key) or _hpos.get(_key.rpartition(':')[2])
                if _hv and _hv.get('aa'):
                    _v['aa_h37rv'] = _hv['aa']
    _sample_meta = parse_sample_meta(args.metadata)   # shared by the dynamics filter and the SNP matrix header
    _dynamics = build_dynamics(parse_metadata(args.metadata), _variants, _sample_meta)   # feeds dynamics + epistasis
    now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    payload = {"generated": now, "version": clean_str(args.version) or "", "repo_url": REPO_URL,
               "counts": counts, "thresholds": thr, "dist": DIST,
               "genome_len": genome_len, "snp_density_ok": snp_density_ok, "defs": DEFS, "nbins": NBINS,
               "genes": parse_gff(args.gff), "lin_colors": parse_lineage_colors(args.lineage_colors),
               "gene_map": build_gene_map(args.gff), "aa2_label": args.aa2_label,
               "section_info": SECTION_INFO,
               "gene_burden": parse_gene_burden(args.gene_burden) or None,
               "dr": parse_dr(args.dr_report),
               "kraken": parse_kraken(args.kraken),
               "pnps": parse_pnps(args.pnps) or None,
               "mask_bins": mask_bins, "mask_pct": mask_pct,
               "mask_iv": (mask_iv[:5000] if mask_iv else None),
               "lineages": lineages, "lin_present": len(lineages) > 0,
               "n_ancient": n_ancient, "anc_thresholds": anc_thr, "provenance": provenance,
               "dynamics": _dynamics,
               "epistasis": build_epistasis(_dynamics),
               "snp_matrix": build_snp_matrix(_variants, provenance.get('reference', '')),
               "sample_meta": _sample_meta,
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
