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
import base64
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
    "dose": ["Per-sample dose from the samplesheet (a numeric annotation, not a QC metric). Feeds the scatter/correlation axes and the dose x treatment test.", ""],
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
# ---- report front-end assets (CSS / JS / HTML shell / logo) live next to this script in
# report_assets/, one file per language so each stays editable with its own tooling; build_html()
# splices them into the single self-contained output HTML. ----
_ASSET_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "report_assets")


def _asset(name):
    with open(os.path.join(_ASSET_DIR, name), encoding="utf-8") as fh:
        return fh.read()


# CSS split into ordered modules under report_assets/css/ (base tokens -> layout -> components ->
# panels); concatenated in cascade order, so the result is identical to a single stylesheet.
_CSS_MODULES = ["css/01_base.css", "css/02_layout.css", "css/03_components.css", "css/04_panels.css"]
CSS = "".join(_asset(m) for m in _CSS_MODULES)

# The report front-end is one ES5 IIFE split into ordered modules under report_assets/js/ (one per
# area) purely for editability; they are concatenated verbatim here, so the running code is identical
# to a single file. Order matters: 01 opens the IIFE and sets up shared state; 14 wires events + closes it.
_JS_MODULES = [
    "js/01_prelude.js", "js/02_qcspace.js", "js/03_state.js", "js/04_helpers.js",
    "js/05_insights_qc.js", "js/06_insights_genome.js", "js/07_render_core.js", "js/08_curation.js",
    "js/09_genome_genes.js", "js/10_dynamics.js", "js/11_epistasis.js", "js/12_snpmatrix.js",
    "js/13_drug_kraken.js", "js/14_boot.js",
]
JS = "".join(_asset(m) for m in _JS_MODULES)

# BAMpiro header logo (bampiro2.png resized to ~90x100 and embedded so the report stays
# self-contained - no external image request). Regenerate from .github/bampiro2.png if the logo changes.
LOGO_DATA_URI = _asset("logo.b64")

SHELL = _asset("shell.html")   # named shell.html (not report_*.html) so .gitignore's report*.html rule can't swallow it


def build_html(title, payload):
    # Embed the payload gzip-compressed + base64 so the whole cohort travels in a light, self-contained
    # HTML (JSON deflates ~5-10x); the browser inflates it natively at load. base64 has no "</" so it
    # cannot break out of the <script>. Nothing is dropped for size - compression is what lets us keep it all.
    raw = json.dumps(payload, separators=(",", ":")).encode("utf-8")
    data_gz = base64.b64encode(gzip.compress(raw, 9)).decode("ascii")
    sys.stderr.write("[qc_report] payload %.2f MB JSON -> %.2f MB embedded (gzip+base64, %.0f%% smaller)\n"
                     % (len(raw) / 1048576.0, len(data_gz) / 1048576.0,
                        100.0 * (1.0 - len(data_gz) / max(1, len(raw)))))
    # Inject content first, then substitute __TITLE__ LAST so a user --title that happens to
    # contain a placeholder token (e.g. "__JS__") can never pull in the CSS/JS/logo/JSON blob.
    return (SHELL.replace("__CSS__", CSS).replace("__LOGO__", LOGO_DATA_URI)
                 .replace("__JS__", JS).replace("__JSON_GZ__", data_gz)
                 .replace("__REPO__", html.escape(REPO_URL))
                 .replace("__TITLE__", html.escape(title)))


# ============================================================================
#  SNP DYNAMICS (optional): flexible metadata + per-sample VCFs -> per-connected-group
#  allele-frequency trajectories over time with flagged events. Self-hides if the metadata
#  has no time+group columns or no VCFs are given.
# ============================================================================
_DYN_SAMPLE_RE = re.compile(r'^(sample_?id|sampleid|sample|name|gid|strain|isolate)$', re.I)
_DYN_TIME_RE   = re.compile(r'(passage|pase|timepoint|time_?point|^time$|^day$|date|week|month|hour|generation|^tp$|visit|^t\d*$)', re.I)
_DYN_GROUP_RE  = re.compile(r'(group|series|patient|host|subject|cluster|experiment|^line$|replicate|chain|pair|lineage_?id|donor|case|animal|^samples$)', re.I)
_DYN_NONSYN    = re.compile(r'missense|stop_gained|stop_lost|start_lost|frameshift|inframe|splice|initiator', re.I)
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
                           'pos_h37rv': (meta or {}).get('pos_h37rv', ''),
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


def _amp(traj):
    """Amplitude (max-min) of a trajectory's non-null allele frequencies; 0 if flat/empty. Used to pick the
    most-variable variants when a series has too many to test all pairs of."""
    vals = [v for v in (traj or []) if v is not None]
    return (max(vals) - min(vals)) if len(vals) >= 2 else 0.0


def build_epistasis(dynamics, min_r=0.8, min_points=3, top=300, perm=2000, max_vars=250, max_pairs_perm=800):
    """Candidate epistatic / linked variant pairs: two SNPs whose allele-frequency trajectories co-vary
    WITHIN a connected series - concordant (rise/fall together) or discordant (one rises as the other
    falls). Score = mean Pearson r of the two trajectories across every series where both move and the
    series has >= min_points timepoints. Significance = a permutation p-value (per pair), then a
    Benjamini-Hochberg FDR q-value across all reported pairs. Recurrence across independent series is
    tracked so a pattern seen in several series (low q) is separable from a lucky single short series.
    None if nothing qualifies. Built from the dynamics payload (moving variants per series).

    Scaling: the pairwise step is O(variants^2) per series and the permutation test is the dominant cost,
    so two caps keep it bounded for large cohorts (many variants and/or many samples -> long trajectories):
    each series contributes at most max_vars variants (the highest-amplitude / most-variable ones), and only
    the max_pairs_perm strongest-|r| candidate pairs get the (expensive) permutation p-value. Both caps are
    reported in the payload so the UI can say what was covered; the full trajectories stay in the dynamics/TSV."""
    if not dynamics or not dynamics.get('groups'):
        return None
    pairs, n_series = {}, 0
    vars_capped, n_vars_max = False, 0
    for g in dynamics['groups']:
        times = g.get('times') or []
        if len(set(times)) < min_points:
            continue
        n_series += 1
        series = g.get('series') or []
        n_vars_max = max(n_vars_max, len(series))
        if len(series) > max_vars:   # too many variants to test every pair: keep the most-variable ones
            vars_capped = True
            series = sorted(series, key=lambda s: -_amp(s.get('traj')))[:max_vars]
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
    # Collect the pairs that clear the |r| threshold, then permutation-test only the strongest ones: the
    # permutation p-value is by far the dominant cost, so with many candidate pairs (common when trajectories
    # are short and spurious high correlations abound) we cap it to the top max_pairs_perm by |mean r|.
    cands = []
    for key, rec in pairs.items():
        rs = rec['rs']
        mean_r = sum(rs) / len(rs)
        if abs(mean_r) < min_r:
            continue
        cands.append((abs(mean_r), key, rec, mean_r))
    n_candidates = len(cands)
    pairs_capped = n_candidates > max_pairs_perm
    if pairs_capped:
        cands.sort(key=lambda c: -c[0])
        cands = cands[:max_pairs_perm]
    out = []
    for _absr, key, rec, mean_r in cands:
        rs = rec['rs']
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
            'vars_capped': vars_capped, 'max_vars': max_vars, 'n_vars_max': n_vars_max,
            'pairs_capped': pairs_capped, 'n_candidates': n_candidates, 'max_pairs_perm': max_pairs_perm,
            'n_concordant': sum(1 for p in out if p['direction'] == 'concordant'),
            'n_discordant': sum(1 for p in out if p['direction'] == 'discordant'),
            'n_strong': sum(1 for p in out if p['tier'] == 'strong')}


def build_snp_matrix(variants, reference='', max_sites=50000):
    """Sparse SNP matrix for the report: samples + one row per SNP site with its ref/alt/annotation
    and only the cells (sample index -> [af, dp]) actually called. The payload is gzip-compressed in the
    HTML, so this cap is a high safety backstop (most-shared sites kept) rather than a size limit; the
    pipeline also writes the complete matrix TSV separately."""
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
                                        'alt': set(), 'gene': '', 'eff': '', 'aa': '', 'aa_h37rv': '', 'pos_h37rv': '', 'cells': {}})
            if v.get('alt'):
                st['alt'].add(v['alt'])
            for k in ('gene', 'eff', 'aa', 'aa_h37rv', 'pos_h37rv'):
                if v.get(k) and not st[k]:
                    st[k] = v[k]
            st['cells'][sidx[s]] = [v['af'], v.get('dp')]
    ordered = sorted(sites.values(), key=lambda x: (x['contig'], x['pos']))
    truncated = len(ordered) > max_sites
    if truncated:
        ordered = sorted(sites.values(), key=lambda x: -len(x['cells']))[:max_sites]
        ordered.sort(key=lambda x: (x['contig'], x['pos']))
    rows = [{'contig': x['contig'], 'pos': x['pos'], 'ref': x['ref'], 'alt': ','.join(sorted(x['alt'])),
             'gene': x['gene'], 'eff': x['eff'], 'aa': x['aa'], 'aa_h37rv': x['aa_h37rv'],
             'pos_h37rv': x['pos_h37rv'], 'n': len(x['cells']), 'cells': x['cells']}
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
    "dosetx": "The per-sample dose (a numeric samplesheet column) split by treatment group, drawn as a box (IQR "
              "+ median) with the individual samples, plus a Kruskal-Wallis rank test of whether dose differs "
              "across the treatment groups (a Mann-Whitney-equivalent when there are two groups). The test runs "
              "over the whole cohort, independent of the live filters; groups with fewer than two dosed samples "
              "are drawn but not tested. Click a point to highlight that sample across the report. Descriptive, "
              "not a claim about efficacy.",
    "vardose": "An association scan: for every variant site, the per-sample allele frequency (0 where the site is "
               "reference) is rank-correlated with dose across the dosed samples (Spearman rho + two-sided p), and "
               "a Benjamini-Hochberg FDR q is computed across all tested variants so the multiple testing is "
               "controlled. The table ranks the variants; click a row to plot its allele-frequency-vs-dose scatter, "
               "click a column header to re-sort, and click a point to highlight that sample. Sites with fewer than "
               "three carriers or no allele-frequency variation are skipped. A screen for dose-associated variants, "
               "NOT proof of causation.",
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
    ap.add_argument("--pos-liftover", default=None,
                    help="TSV 'mapping_pos<TAB>canonical_pos' (e.g. pathotypr_liftover.py apply --out-map) -> the "
                         "reference-of-interest COORDINATE per variant, alignment-free. Fills pos_h37rv for the SNP "
                         "tables; an alternative to --vcfs-h37rv when the references do not share coordinates.")
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
    # a numeric 'dose' samplesheet column becomes a first-class (hidden-by-default) metric so it
    # can go on the scatter axes + correlation matrix and drive the dose x treatment test.
    dose_map = parse_dose(args.metadata)
    if dose_map:
        extra_metrics = extra_metrics + [{"key": "dose", "label": "Dose", "kind": "float", "dir": "neu"}]
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
                         "m": dict({k: to_float(m.get(k)) for k in metric_keys + extra_keys + EFF_KEYS},
                                   **({"dose": dose_map.get(sid)} if dose_map else {}))})

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
    if args.vcfs_h37rv:   # attach the reference-of-interest (H37Rv) amino-acid change AND coordinate per variant
        _h37 = parse_vcfs(args.vcfs_h37rv)
        for _s, _pm in _variants.items():
            _hs = _h37.get(_s, {})
            # Pair each used-ref variant with its canonical-reference record. Prefer the ORIGINAL used-ref
            # coordinate a liftover stamped (OriginalContig/Start or OPOS) so the two references need NOT
            # share coordinates; else fall back to the bare position (references that DO share H37Rv coords).
            _by_opos, _by_pos = {}, {}
            for _ck, _cv in _hs.items():
                if _cv.get('opos'):
                    _by_opos[_cv['opos']] = (_ck, _cv)
                _by_pos.setdefault(_ck.rpartition(':')[2], (_ck, _cv))
            for _key, _v in _pm.items():
                _hit = _by_opos.get(_key) or _by_pos.get(_key.rpartition(':')[2])
                if _hit:
                    _ck, _cv = _hit
                    if _cv.get('aa'):
                        _v['aa_h37rv'] = _cv['aa']
                    _v['pos_h37rv'] = _ck   # reference-of-interest coordinate (contig:pos), possibly != the mapping one
    if args.pos_liftover:   # alignment-free canonical coordinate per variant (pathotypr k-mer liftover map)
        _lift = {}
        try:
            with open(args.pos_liftover, encoding="utf-8", errors="replace") as _fh:
                for _line in _fh:
                    _c = _line.rstrip("\n").split("\t")
                    if len(_c) >= 2 and _c[0].strip().isdigit() and _c[1].strip().isdigit():
                        _lift[_c[0].strip()] = _c[1].strip()
        except OSError as _e:
            sys.stderr.write("[qc_report] WARN pos-liftover (%s): %s\n" % (args.pos_liftover, _e))
        for _s, _pm in _variants.items():
            for _key, _v in _pm.items():
                _p = _key.rpartition(":")[2]
                if not _v.get("pos_h37rv") and _p in _lift:
                    _v["pos_h37rv"] = "%s:%s" % (args.aa2_label, _lift[_p])
    _sample_meta = parse_sample_meta(args.metadata)   # shared by the dynamics filter and the SNP matrix header
    _dynamics = build_dynamics(parse_metadata(args.metadata), _variants, _sample_meta)   # feeds dynamics + epistasis
    # front-load 'dose' so it survives the correlation matrix's top-N view
    dist_keys = (["dose"] + DIST) if dose_map else DIST
    now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    payload = {"generated": now, "version": clean_str(args.version) or "", "repo_url": REPO_URL,
               "counts": counts, "thresholds": thr, "dist": dist_keys,
               "genome_len": genome_len, "snp_density_ok": snp_density_ok, "defs": DEFS, "nbins": NBINS,
               "genes": parse_gff(args.gff), "lin_colors": parse_lineage_colors(args.lineage_colors),
               "gene_map": build_gene_map(args.gff), "aa2_label": args.aa2_label,
               "ref_name": provenance.get('reference', ''),   # mapping reference name for the SNP tables' headers
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
