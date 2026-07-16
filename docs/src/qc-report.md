# Interactive QC Report

Every run produces a **single self-contained HTML file**
(`<samplesheet>_qc_report.html` - no internet or CDN needed) that folds the whole
cohort into one interactive dashboard, plus a machine-readable
`<samplesheet>_qc_flags.tsv` of per-sample **PASS/WARN/FAIL** verdicts. It is built
by `bin/qc_report.py` and controlled by `--make_qc_report` (default `true`).

## Live & interactive

- **Live thresholds & presets** - edit any QC cut-off (depth, breadth, missing,
  duplication, mapping, IUPAC, Ti/Tv, SNP-z, heteroplasmy, mixed-lineage) and the
  whole report re-flags instantly. Presets: *gate defaults*, *strict (modern WGS)*,
  *lenient (aDNA / low-cov)*.
- **Collapsible table-of-contents sidebar** with scroll-spy; panels with no data
  hide themselves (and their nav link).
- **Per-section (i) info popovers** explaining each analysis and its caveats.
- **Exclusion basket** - tick samples (via the table, a drag-box in the scatter, or
  *basket all flagged*); FAIL samples start pre-selected, and you export
  `exclusion.tsv` / `keep_list.txt` for downstream phylogeny. The margin behind every
  flag is shown inline in the Flagged panel (no hover needed).
- **Triage-first layout** - the sample table opens sorted worst-QC first, and the
  **Flagged** panel sits directly under the general statistics (and at the top of
  the sidebar) so the failing samples are the first thing you reach.
- **Dark / light theme** - a header toggle that follows your OS preference and is
  remembered per viewer.
- **Responsive** - reflows to a phone: the frozen sample column is capped, wide
  tables and matrices scroll horizontally, and the contents sidebar becomes a
  dismissable drawer.
- The header shows the **pipeline version** and links to the **source on GitHub**;
  most panels have a fullscreen (expand) view; the whole report prints / saves to PDF.

## Panels

Grouped as in the sidebar:

| Group | Panels |
| :--- | :--- |
| **Overview** | General statistics (value-coloured, sortable, filterable, TSV export) · **Flagged samples** (worst-first, each failing margin shown inline) · Per-lineage summary (canonical *mycolorsTB* palette) · Distributions (beeswarm / bar / histogram) · **Taxonomic composition** (Kraken2: primary taxon / contaminants / unclassified, worst-first) |
| **Correlation & structure** | Metric-pair scatter (box-select to basket) · Metric correlation heatmap · QC-space PCA (+ most-unusual-samples table) · Divergence vs completeness |
| **Genome & genes** | Consensus completeness · Genome landscape (per-position callability / variant heatmap, gene search, mask-region toggle) · Functional annotation (snpEff classes) · Functional gene burden · Variable genes (SNP-density hotspots) |
| **Evolution** | Temporal sampling overview · Selection pN/pS (dN/dS, eskaks) · aDNA damage authentication (mapDamage) |
| **Variants over time** | **SNP dynamics** (allele-frequency trajectories over time, per-timepoint DP bars, zoom, series filter) · **Epistasis** (co-varying variant pairs, permutation *p* + BH-FDR *q*, cards / matrix / table views) · **SNP matrix** (site × sample AF matrix, metadata column filter, TSV export) · **Drug resistance** (sample × drug WHO-grade matrix) |

Optional panels appear only when their input is present: SNP dynamics / epistasis /
SNP matrix need `--metadata` + the per-sample VCFs; gene burden needs the cohort
burden TSV; pN/pS needs an [eskaks](https://github.com/PathoGenOmics-Lab/eskaks)
TSV; the [drug resistance](pathotypr.md) panel needs the pathotypr DR calls; aDNA
needs mapDamage output; the taxonomic-composition panel needs the Kraken2 `.report`
files (i.e. Kraken decontamination ran).
