"""The metric registry, the flag thresholds and the QC verdict engine.

DEF / ANC_DEF are the default cut-offs (the modern view and the looser ancient-DNA
view); METRICS is the ordered table of every metric the report knows about, DEFS
its in-report documentation, and DIST the subset drawn as distributions. Around
them sit the derived statistics (robust, het_frac, the lineage helpers) and
flag_sample, which turns one summary row into a PASS / WARN / FAIL verdict plus
the list of reasons.

Everything here is pure in-memory logic over an already-parsed summary row.
"""
from __future__ import annotations

import re

from .parsers import parse_gene_locus, to_float

DEF = dict(depth_min=10.0, breadth_min=90.0, missing_max=10.0, dup_max=40.0,
           iupac_max=2.0, mapping_min=80.0, titv_min=0.0, snp_z=3.0,
           het_max_frac=3.0, mixed_min_frac=2.0)  # het_max_frac + mixed_min_frac are PERCENTS

# Ancient (aDNA) sample threshold view: looser cov/breadth/missing/mapping/dup/iupac (a 5x mummy must not FAIL the
# modern gate), plus damage_min_ct = the terminal 5' C>T deamination floor for the DAMAGE_LOW authentication screen.
ANC_DEF = dict(depth_min=3.0, breadth_min=30.0, missing_max=70.0, dup_max=60.0,
               iupac_max=5.0, mapping_min=20.0, damage_min_ct=0.05)
FAIL_FLAGS = {"LOW_DEPTH", "LOW_BREADTH", "HIGH_MISSING", "NO_DATA"}

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
