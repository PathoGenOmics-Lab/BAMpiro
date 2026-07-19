# Reading the QC report

Every run drops a single self-contained `<samplesheet>_qc_report.html` next to your outputs — the whole cohort folded into one interactive dashboard. This tutorial is a guided tour: by the end you can open it, triage a cohort worst-first, tune thresholds on the fly, and build an exclusion list for downstream work. Allow about **15 minutes**, ideally with the [interactive demo report](../examples/qc_report_demo.html){ target="_blank" } (a full 17-sample MTBC cohort) open in another tab so you can click along.

!!! note "Where it comes from"

    The report is produced by default (`--make_qc_report true`) alongside a machine-readable
    `<samplesheet>_qc_flags.tsv` of per-sample **PASS / WARN / FAIL** verdicts. See
    [Interactive QC Report](../qc-report.md) for the exhaustive panel-by-panel reference —
    this page just teaches the workflow.

## 1. Start at the top: the executive summary

Open the HTML file in any browser (no internet needed). The **Executive summary** is the first panel: cohort health and the headline findings at a glance. Read it before anything else — it tells you whether you have a clean cohort or a triage job on your hands.

Right below the general **Stats**, the sidebar puts the **Flagged samples** panel near the very top on purpose. The layout is triage-first: the sample table opens sorted **worst-QC first**, and the failing samples are the first thing you reach.

## 2. Triage the flagged samples

Work the **Flagged** panel top-down — it is ordered worst-first, and the **margin** behind every flag is shown inline, so you can see *how far* a sample missed each cut-off without hovering.

To decide what to drop, use the **exclusion basket**:

1. **Tick** the samples you want to exclude — from the table, by drawing a **box-select** in any scatter plot, or with **basket all flagged** in one click. FAIL samples start **pre-selected**.
2. When the basket reflects your call, **export** it. You get two files:
    - `exclusion.tsv` — the samples you are dropping, with reasons.
    - `keep_list.txt` — the complement, ready to feed the phylogeny step.

!!! tip "The basket is yours, not the pipeline's"

    Basketing and exporting are a **human decision layer**. They never alter the run —
    the machine-readable `_qc_flags.tsv` gate verdicts stay exactly as computed.

## 3. Re-flag live with thresholds & presets

The report re-computes flags **in your browser**. Open the **thresholds** controls and nudge any cut-off — depth, breadth, missing, duplication, mapping, IUPAC, Ti/Tv, SNP-z, heteroplasmy, mixed-lineage — and the **whole report re-flags instantly**. Or apply a **preset**:

=== "gate defaults"

    The pipeline's own cut-offs — what produced `_qc_flags.tsv`. Your starting point.

=== "strict (modern WGS)"

    Tighter thresholds for high-coverage modern sequencing.

=== "lenient (aDNA / low-cov)"

    Relaxed cut-offs for ancient-DNA or low-coverage cohorts.

!!! warning "Live threshold ≠ pipeline gate"

    Editing thresholds re-flags the **view** so you can explore "what if". It does **not**
    change the run or the exported `_qc_flags.tsv`. To change the actual gate, re-run with
    the relevant parameters (see [Configuring a run](configuring-a-run.md)).

## 4. Narrow the view with the toolbar

Three toolbar controls restrict what you are looking at:

- **Search** — jump to a sample by name.
- **Only flagged** — hide the passing samples and focus on problems.
- **Metadata cohort filter** — when your samplesheet carries categorical columns (e.g. `treatment`, `site`, `ward`), each gets a **dropdown** that restricts the *entire* view — sample table, distribution and QC-space plots, and the genome & gene panels — to one value. (Time, group and lineage columns are handled separately; lineage has its own filter.)

## 5. A quick pass over the panel groups

The sidebar groups every panel; panels with no data hide themselves. In brief:

- **Overview** — executive summary, general statistics, flagged samples, per-lineage summary, distributions, taxonomic composition.
- **Correlation & structure** — metric-pair scatter, correlation heatmap, QC-space PCA, divergence vs completeness, **Dose × treatment**.
- **Genome & genes** — consensus completeness, genome landscape, functional annotation & gene burden, variable genes.
- **Evolution** — temporal sampling, selection pN/pS, aDNA damage authentication.
- **Variants over time** — SNP dynamics, epistasis, SNP matrix, **Variant × dose**, drug resistance.

Every panel carries an **(i) info popover** explaining its method and caveats — click it whenever you are unsure what a plot means. The full catalogue lives in the [Interactive QC Report](../qc-report.md) reference.

!!! tip "Longitudinal cohorts"

    If your samples span timepoints, **SNP dynamics** and **Epistasis** turn into a within-host
    evolution story — see [Longitudinal & within-host analysis](longitudinal-analysis.md).

## 6. Spotting contamination (Kraken2)

Before mapping, BAMpiro screens every sample's reads with **Kraken2**; the **Taxonomic composition** panel (under *Overview*) shows the breakdown per sample, sorted **worst-first**: the **primary** taxon %, the largest **contaminants**, and the **unclassified** fraction.

A clean sample is almost entirely one taxon. Warning signs:

- **Low primary %** — the executive summary flags samples below **90 % primary** as *possibly contaminated*.
- **A large secondary taxon** — co-infection, cross-contamination, or the wrong reference for these reads.
- **High unclassified %** — material absent from the Kraken2 database (a divergent or novel organism), or a mismatched reference.

Basket obviously contaminated samples for exclusion alongside the QC failures. If a whole batch shows the *same* contaminant, suspect a shared lab / reagent source rather than the biology.

## 7. The dose analyses

If your metadata has a numeric **`dose`** column, two extra panels light up.

**Dose × treatment** (under *Correlation & structure*) draws a per-treatment dose distribution (box + points) and runs a **Kruskal–Wallis** rank test of whether dose differs across the treatment groups (a Mann–Whitney-equivalent for two groups). Groups with fewer than two dosed samples are drawn but not tested. **Click a group name** to filter the whole report to that treatment — the same cohort filter as the toolbar dropdown — or click a point to highlight that sample everywhere.

**Variant × dose** (under *Variants over time*) needs per-sample VCFs as well. For every variant site it correlates the per-sample **allele frequency** with dose (**Spearman** ρ + two-sided *p*), applies a **Benjamini–Hochberg FDR** across all tested variants, and ranks them — so real dose-tracking mutations surface while multiple-testing keeps incidental hits in check. Sites carried by fewer than three samples (or with no AF variation) are skipped. The table is **sortable** (click a header to sort by ρ, *p*, *q* or carrier count) and **searchable** (gene / position / amino acid); **click a row** to plot that variant's AF-vs-dose scatter, and **click a point** to highlight that sample across the report.

## What's next

You can now triage a cohort and hand a clean `keep_list.txt` downstream. Next, dig into what the samples *are*:

[:octicons-arrow-right-24: Lineage & drug-resistance typing](lineage-and-dr-typing.md)
