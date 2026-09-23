# Interactive QC Report

Every run produces a **single self-contained HTML file**
(`<samplesheet>_qc_report.html` - no internet or CDN needed) that reads the whole cohort
back as a short report, plus a machine-readable `<samplesheet>_qc_flags.tsv` of per-sample
**PASS/WARN/FAIL** verdicts. It is built by `bin/qc_report.py` and controlled by
`--make_qc_report` (default `true`).

!!! tip "See it live"

    Explore a full example report built from a 17-sample demo cohort - every page
    populated, fully interactive (live thresholds, the exclusion list, dark mode):

    [:octicons-play-16: Open the interactive demo report](examples/qc_report_demo.html){ .md-button .md-button--primary target="_blank" rel="noopener" }

    *Synthetic demo data. It opens as a standalone page - the very same self-contained
    HTML file each real run produces.*

## How it reads

The report is eight pages, listed in the sidebar in the order a reader needs them. A page
whose data the run did not produce disappears with its entry, and the rest are renumbered.
The sidebar badge next to a page says whether it needs attention: the samples to exclude
on *Sample QC*, the samples with resistance mutations their lineage does not share on
*Resistance*.

| Page | What it answers |
| :--- | :--- |
| **1 · Summary** | What the run found, as sentences with their numbers: an *In short* paragraph, then one card per finding (sample QC, identity, resistance, variants over time, lineages, gene conversion, coverage). Every card links to the page that holds its evidence and says what it does not prove. |
| **2 · Sample QC** | Which samples can be trusted: the verdicts, the live thresholds, the **flagged samples** with the value behind every flag, the **exclusion list** and its exports, the table of all samples, lineages, **contamination** (Kraken2), distributions and aDNA damage. |
| **3 · Genome & genes** | Consensus completeness, the genome landscape (missing calls, **deletions**, SNPs per bin and **per callable kb**, het, indels), the **deletions** found by comparing each sample's stretches without reads with the rest of the cohort, functional impact (snpEff), gene burden, variable genes and dN/dS. |
| **4 · Relatedness** | SNP distances between the consensus sequences, over the positions both samples called: a heatmap in the order that keeps each cluster together, the **clusters** at a threshold you can move, and the samples far from the rest of their own group (patient, line, series). |
| **5 · Variants over time** | SNP dynamics, co-varying pairs (epistasis), the SNP matrix and variant &#215; dose. |
| **6 · Resistance** | WHO-catalogue mutations grouped by mutation, the sample &#215; drug matrix and every call. |
| **7 · Gene conversion** | Candidate tracts and the evidence behind each verdict. |
| **8 · Diagnostics** | Metric pairs, the correlation matrix, QC space, divergence vs completeness, dose &#215; treatment and sampling dates: views for digging into a problem, none of which flags a sample. |

Printing (or saving as PDF) lays out every page one after another.

## Reading the results

- **Flags in words.** A flag is shown as what went wrong and the value against the rule it
  broke: *Low depth 2.2&#215;, needs 10&#215;*; *Incomplete consensus 93.7% of the consensus
  missing, at most 10%*. The codes (`LOW_DEPTH`, `HIGH_MISSING`...) stay in the tooltips and in
  every export, so `exclusion.tsv` and `qc_flags.tsv` do not change.
- **Flags the pipeline decided stay.** Changing a threshold re-flags the page, but a flag that
  needs evidence the page does not hold is kept rather than recomputed away:
  `LINEAGE_MISMATCH` compares a sample's lineage with the other samples on its reference, and
  the flagged list says both (*types as A4; the rest of E1ASM0035 types as L7*).
- **The table of all samples** opens on the fourteen metrics a verdict is made of, in reading
  order; the other metrics are under *columns*. A cell is coloured only when its value tripped
  a flag for that sample, red when the flag fails the sample and amber when it asks for a look.
- **Contamination is measured against what the cohort is.** Every Kraken2 report of a sample is
  merged (summing reads), so a sample sequenced over five runs is one row. Each sample is
  placed in the clade its reads actually sit in rather than at species level: most
  *M. tuberculosis* reads stop at the complex, so a species-level reading put every clean
  culture at about 5% and called it contaminated. The *target* is the clade most samples are
  dominated by; below 90% of the classified reads in it a sample is **mixed**, below 50% it is
  **another organism**.
- **Resistance a lineage shares is set apart.** A mutation carried by at least 90% of the samples
  of a lineage is that lineage's own (pncA H57D in every *M. bovis* is why *M. bovis* resists
  PZA), so it is listed as a lineage marker and not counted among the resistance the cohort
  acquired.
- **Big panels open small.** The SNP dynamics draw 24 trajectories at a time, list the first 40
  genes and the 12 genes rising in most series; the resistance panel opens on grade 1&#8211;2 calls,
  grouped by mutation; the gene-conversion table lists the called events first. Everything else
  is one click away.
- **Written read-outs** at the top of each panel state its result in a sentence; the header
  *read-outs* button hides them all. They list at most eight samples or genes and say how many
  more there are.

## Live & interactive

- **Live thresholds & presets** - edit any QC cut-off (depth, breadth, missing,
  duplication, mapping, IUPAC, Ti/Tv, SNP-z, heteroplasmy, mixed-lineage) on the *Sample QC*
  page and the whole report re-flags instantly, the summary included. Presets: *gate
  defaults*, *strict (modern WGS)*, *lenient (aDNA / low-cov)*.
- **Filter by metadata** - when the samplesheet carries categorical annotation columns
  (`treatment`, `site`, `ward`…), the toolbar shows a dropdown per column that restricts
  the whole QC view - the sample table, the distribution / QC-space plots, and the genome
  & gene panels - to one value, exactly like the lineage and flag filters. The time and
  group (dynamics) columns and the lineage column are left out (lineage has its own filter).
- **Per-section (i) info popovers** explaining each analysis and its caveats.
- **Exclusion list** - tick samples (in the flagged list, the table, or with a drag-box in the
  metric-pairs plot); FAIL samples start in it, and you export `exclusion.tsv` /
  `keep_list.txt` for the downstream analysis.
- **Sample profile** - click a sample name for every metric, its flags in words and its genome
  profile.
- **Dark / light theme** - a header toggle that follows your OS preference and is
  remembered per viewer.
- **Responsive** - reflows to a phone: the sidebar becomes a drawer, wide tables and matrices
  scroll horizontally.
- The header shows the **pipeline version** and links to the **source on GitHub**;
  most panels have a fullscreen (expand) view.

## Optional inputs

Nothing below is required - each optional input just lights up the matching panel or
feature, and any panel with no data hides itself and its nav link. The pipeline wires
these from the run's own outputs; you only get what you have.

| Input | Adds |
| :--- | :--- |
| `--metadata <samplesheet.tsv>` | SNP-dynamics grouping + the SNP-matrix column-header levels (see **Optional metadata** below) |
| `--vcfs <sample.vcf …>` | Per-SNP allele frequencies → SNP dynamics, epistasis, SNP matrix |
| `--gff <genes.gff3>` | Variable-genes hotspots, genome-landscape gene search, and the **gene → Mycobrowser (H37Rv) locus-tag** links |
| `--dr-report <dr.tsv>` | **Resistance** page (mutations grouped across samples, lineage markers set apart, the sample × drug matrix) - the [pathotypr](pathotypr.md) DR calls |
| `--kraken <sample.report …>` | **Contamination** panel (Kraken2): one row per sample, its runs merged, measured against the cohort's target clade. The pipeline only has these reports to pass when you set `--kraken2_db`, which has no default, so the panel hides itself otherwise |
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
there is nothing to configure - add a column and the matching panel reacts.

| Column (matched by name) | Examples | What it drives |
| :--- | :--- | :--- |
| **sample id** | `sample`, `sample_id`, `name`, `gid`, `strain`, `isolate`… | Keys the metadata to each sample's VCF / stats (falls back to the first column). |
| **time** | `timepoint`, `day`, `date`, `week`, `month`, `hour`, `passage`, `generation`, `visit`, `tp`, `t0`… | The x-axis of the **SNP dynamics** trajectories. |
| **group / series** | `group`, `patient`, `series`, `host`, `subject`, `cluster`, `experiment`, `donor`, `case`, `replicate`, `chain`, `samples`… | Connects samples into one longitudinal series (a trajectory set per group) for **SNP dynamics** and **epistasis**. |
| **any other column** | `treatment`, `site`, `region`, `ward`, `batch`… | A categorical annotation - becomes a **cohort filter** dropdown in the toolbar (restrict the whole QC view to one value) and a SNP-matrix header level. |
| **collection date** | `collection_date`, `sampling_date`, `isolation_date`, `date`, `year`, `fecha` | The *Sampling dates* panel (Diagnostics) and the date in each sample's profile. The day the pipeline processed a sample is never read as one. |
| **dose** | `dose`, `dosis` | A **numeric** column - becomes a first-class metric (selectable on the scatter axes + the correlation matrix, with a Spearman *r* + *p* read-out) and drives the **Dose × treatment** test. |

**Every** annotation column - including the time and group ones - also becomes a
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
It becomes a first-class **metric** - pick it on either axis of **Metric pairs** (or
in the **Correlation matrix**, both on the *Diagnostics* page) to get a Spearman *r* with a two-sided *p*-value
against any QC or genomic metric - and it powers a dedicated **Dose × treatment** panel:
a per-treatment dose distribution (box + points) with a **Kruskal–Wallis** rank test of
whether dose differs across the treatment groups (a Mann–Whitney-equivalent when there
are two groups). The test runs over the whole cohort; groups with fewer than two dosed
samples are drawn but not tested, and the panel hides itself when there is no `dose`
column or no treatment column. Click a group name to filter the whole report to that
treatment (it drives the same cohort filter as the toolbar dropdown), or click a point to
highlight that sample everywhere.

When per-sample **VCFs** are also present, a **Variant × dose** panel (under *Variants
over time*) runs an association scan: for every variant site it correlates the per-sample
**allele frequency** (0 where the site is reference) with dose across the dosed samples
(**Spearman** ρ + two-sided *p*), applies a **Benjamini–Hochberg FDR** across all tested
variants, and ranks them - so you can see which mutations track the dose while the
multiple-testing correction keeps incidental hits in check. Click any row to plot that
variant's allele-frequency-vs-dose scatter, click a column header to re-sort (by ρ, *p*,
*q* or carrier count), search by gene / position / amino acid, and click a point to
highlight that sample across the whole report. Sites carried by fewer than three samples
(or with no allele-frequency variation) are skipped. Both this panel and **Dose × treatment** carry an info **(i)** with the full
method note (as every analysis panel does).
