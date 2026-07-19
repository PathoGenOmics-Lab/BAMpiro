# Interactive QC Report

Every run produces a **single self-contained HTML file**
(`<samplesheet>_qc_report.html` - no internet or CDN needed) that folds the whole
cohort into one interactive dashboard, plus a machine-readable
`<samplesheet>_qc_flags.tsv` of per-sample **PASS/WARN/FAIL** verdicts. It is built
by `bin/qc_report.py` and controlled by `--make_qc_report` (default `true`).

!!! tip "See it live"

    Explore a full example report built from a 17-sample demo cohort — every panel
    populated, fully interactive (live thresholds, the exclusion basket, dark mode):

    [:octicons-play-16: Open the interactive demo report](examples/qc_report_demo.html){ .md-button .md-button--primary target="_blank" rel="noopener" }

    *Synthetic demo data. It opens as a standalone page — the very same self-contained
    HTML file each real run produces.*

## Live & interactive

- **Live thresholds & presets** - edit any QC cut-off (depth, breadth, missing,
  duplication, mapping, IUPAC, Ti/Tv, SNP-z, heteroplasmy, mixed-lineage) and the
  whole report re-flags instantly. Presets: *gate defaults*, *strict (modern WGS)*,
  *lenient (aDNA / low-cov)*.
- **Filter by metadata** - when the samplesheet carries categorical annotation columns
  (`treatment`, `site`, `ward`…), the toolbar shows a dropdown per column that restricts
  the whole QC view - the sample table, the distribution / QC-space plots, and the genome
  & gene panels - to one value, exactly like the lineage and flag filters. The time and
  group (dynamics) columns and the lineage column are left out (lineage has its own filter).
- **Collapsible table-of-contents sidebar** with scroll-spy; panels with no data
  hide themselves (and their nav link).
- **Per-section (i) info popovers** explaining each analysis and its caveats.
- **Exclusion basket** - tick samples (via the table, a drag-box in the scatter, or
  *basket all flagged*); FAIL samples start pre-selected, and you export
  `exclusion.tsv` / `keep_list.txt` for downstream phylogeny. The margin behind every
  flag is shown inline in the Flagged panel (no hover needed).
- **Triage-first layout** - the sample table opens sorted worst-QC first, and the
  **Flagged** panel sits directly under the general statistics (near the top of the
  sidebar, right under Summary and Stats) so the failing samples are the first thing
  you reach.
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
| **Overview** | **Executive summary** (cohort health & headline findings at a glance — the first panel) · General statistics (value-coloured, sortable, filterable, TSV export) · **Flagged samples** (worst-first, each failing margin shown inline) · Per-lineage summary (canonical *mycolorsTB* palette) · Distributions (beeswarm / bar / histogram) · **Taxonomic composition** (Kraken2: primary taxon / contaminants / unclassified, worst-first) |
| **Correlation & structure** | Metric-pair scatter (box-select to basket, with a Spearman *r* + *p* read-out) · Metric correlation heatmap · QC-space PCA (+ most-unusual-samples table) · Divergence vs completeness · **Dose × treatment** (per-group dose distribution + Kruskal–Wallis test) |
| **Genome & genes** | Consensus completeness · Genome landscape (per-position callability / variant heatmap, gene search, mask-region toggle) · Functional annotation (snpEff classes) · Functional gene burden · Variable genes (SNP-density hotspots) |
| **Evolution** | Temporal sampling overview · Selection pN/pS (dN/dS, eskaks) · aDNA damage authentication (mapDamage) |
| **Variants over time** | **SNP dynamics** (allele-frequency trajectories over time, per-timepoint DP bars, zoom, series filter) · **Epistasis** (co-varying variant pairs, permutation *p* + BH-FDR *q*, cards / matrix / table views) · **SNP matrix** (site × sample AF matrix, metadata column filter, TSV export) · **Variant × dose** (per-variant AF ~ dose Spearman scan with BH-FDR) · **Drug resistance** (sample × drug WHO-grade matrix) |

## Optional inputs

Nothing below is required — each optional input just lights up the matching panel or
feature, and any panel with no data hides itself and its nav link. The pipeline wires
these from the run's own outputs; you only get what you have.

| Input | Adds |
| :--- | :--- |
| `--metadata <samplesheet.tsv>` | SNP-dynamics grouping + the SNP-matrix column-header levels (see **Optional metadata** below) |
| `--vcfs <sample.vcf …>` | Per-SNP allele frequencies → SNP dynamics, epistasis, SNP matrix |
| `--gff <genes.gff3>` | Variable-genes hotspots, genome-landscape gene search, and the **gene → Mycobrowser (H37Rv) locus-tag** links |
| `--dr-report <dr.tsv>` | **Drug resistance** panel (WHO-grade sample × drug matrix) — the [pathotypr](pathotypr.md) DR calls |
| `--kraken <sample.report …>` | Taxonomic composition / contamination panel (Kraken2) |
| `--gene-burden <burden.tsv>` | Functional gene-burden panel |
| `--pnps <dnds.tsv>` | Selection pN/pS panel ([eskaks](https://github.com/PathoGenOmics-Lab/eskaks)) |
| `--mapdamage-dir <dir>` | aDNA damage-authentication panel (mapDamage) |
| `--mask-bed <mask.bed>` | Mask-region toggle in the genome-landscape panel |
| `--consensus <sample.fa …>` | Consensus-completeness panel |
| `--lineage-colors <palette.tsv>` | Canonical lineage palette (`lineage<TAB>#hex`, e.g. *mycolorsTB*) |
| `--vcfs-h37rv <sample.vcf …>` | Canonical (H37Rv) **amino-acid** numbering shown beside the used-reference one |
| `--pos-liftover <map.tsv>` | Canonical (H37Rv) **coordinate** per variant, alignment-free (`pathotypr_liftover.py`) |
| `--aa2-label <name>` | Label for the canonical numbering (default `H37Rv`) |

## Optional metadata

The report reads richer context from the run **samplesheet** (or any TSV passed as
`--metadata`). Every column is optional and matched **by name** (case-insensitive), so
there is nothing to configure — add a column and the matching panel reacts.

| Column (matched by name) | Examples | What it drives |
| :--- | :--- | :--- |
| **sample id** | `sample`, `sample_id`, `name`, `gid`, `strain`, `isolate`… | Keys the metadata to each sample's VCF / stats (falls back to the first column). |
| **time** | `timepoint`, `day`, `date`, `week`, `month`, `hour`, `passage`, `generation`, `visit`, `tp`, `t0`… | The x-axis of the **SNP dynamics** trajectories. |
| **group / series** | `group`, `patient`, `series`, `host`, `subject`, `cluster`, `experiment`, `donor`, `case`, `replicate`, `chain`, `samples`… | Connects samples into one longitudinal series (a trajectory set per group) for **SNP dynamics** and **epistasis**. |
| **any other column** | `treatment`, `site`, `region`, `ward`, `batch`… | A categorical annotation - becomes a **cohort filter** dropdown in the toolbar (restrict the whole QC view to one value) and a SNP-matrix header level. |
| **dose** | `dose`, `dosis` | A **numeric** column - becomes a first-class metric (selectable on the scatter axes + the correlation matrix, with a Spearman *r* + *p* read-out) and drives the **Dose × treatment** test. |

**Every** annotation column — including the time and group ones — also becomes a
**column-header level** in the **SNP matrix** (filter the matrix by it, hover a header
for its value); the time and group columns *additionally* drive the dynamics, and each
plain categorical column *additionally* becomes a **cohort filter** in the toolbar. The
pipeline's own file/reference columns (`r1`, `r2`, `reffasta`, `refgff`, `refid`,
`taxid`, `runid`) are ignored, so a normal BAMpiro samplesheet works unchanged. The
SNP-dynamics / epistasis panels appear only when the metadata has **both** a time and
a group column *and* per-sample VCFs are present; the SNP matrix needs the VCFs alone.

Example `--metadata` TSV (tab-separated; the demo cohort):

```tsv
sample      samples   timepoint   site      treatment     dose   lineage
TB-P1-d0    TB-P1     0           Madrid    HRZE          360    L2
TB-P1-d60   TB-P1     60          Madrid    HRZE          360    L2
TB-P1-d180  TB-P1     180         Madrid    HRZE          450    L2
TB-2020-C   .         .           Sevilla   MDR regimen   740    L4
```

Here `samples` is the group (the longitudinal patient series *TB-P1*), `timepoint` is
the time axis, `treatment` and `site` become cohort filters (and SNP-matrix header
levels), `dose` becomes a numeric metric, and `lineage` gets its canonical palette; the
singleton `TB-2020-C` (no series) simply has no dynamics trajectory.

### Dose × treatment and quantitative metadata

A numeric **`dose`** column is treated as a quantitative variable rather than a category.
It becomes a first-class **metric** — pick it on either **Correlations** scatter axis (or
in the **Metric correlation** matrix) to get a Spearman *r* with a two-sided *p*-value
against any QC or genomic metric — and it powers a dedicated **Dose × treatment** panel:
a per-treatment dose distribution (box + points) with a **Kruskal–Wallis** rank test of
whether dose differs across the treatment groups (a Mann–Whitney-equivalent when there
are two groups). The test runs over the whole cohort; groups with fewer than two dosed
samples are drawn but not tested, and the panel hides itself when there is no `dose`
column or no treatment column.

When per-sample **VCFs** are also present, a **Variant × dose** panel (under *Variants
over time*) runs an association scan: for every variant site it correlates the per-sample
**allele frequency** (0 where the site is reference) with dose across the dosed samples
(**Spearman** ρ + two-sided *p*), applies a **Benjamini–Hochberg FDR** across all tested
variants, and ranks them — so you can see which mutations track the dose while the
multiple-testing correction keeps incidental hits in check. Click any row to plot that
variant's allele-frequency-vs-dose scatter, click a column header to re-sort (by ρ, *p*,
*q* or carrier count), and click a point to highlight that sample across the whole report.
Sites carried by fewer than three samples (or with no allele-frequency variation) are
skipped. Both this panel and **Dose × treatment** carry an info **(i)** with the full
method note (as every analysis panel does).
