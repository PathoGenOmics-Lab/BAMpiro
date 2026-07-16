

# BAMpiro 🧛‍♂️🧬
### *General Bacterial Short Read Mapping, Variant Calling & Lineage/DR Typing Pipeline*
__Paula Ruiz-Rodriguez<sup>1</sup>__ 
__and Mireia Coscolla<sup>1</sup>__
<br>
<sub> 1. I<sup>2</sup>SysBio, University of Valencia-CSIC, FISABIO Joint Research Unit Infection and Public Health, Valencia, Spain </sub>  

<table>
  <tr>
    <td width="300">
      <img src=".github/bampiro2.png" title="BAMpiro logo" style="width:300px; height: auto;">
    </td>
    <td style="padding-left: 20px;">
      <p>
        <strong>BAMpiro</strong> is a modular, containerized bioinformatics pipeline built with <strong>Nextflow (DSL2)</strong>. While optimized by default for <em>Mycobacterium tuberculosis</em> (TB), its architecture is <strong>agnostic</strong> and can be used to analyze <strong>any bacterial genome</strong> (e.g., <em>E. coli</em>, <em>Salmonella</em>, <em>Staphylococcus</em>) by adjusting a few parameters.
      </p>
      <p>The pipeline automates the workflow from raw reads to annotated variants, consensus sequences, and comprehensive quality control reports.</p>
    </td>
  </tr>
</table>

## Key Features

* **Universal Bacterial Support:** Works with any reference genome and GFF annotation.
* **Automated Reference Prep:** Automatically indexes genomes and builds SnpEff databases on the fly.
* **Repeat Masking:** Uses `nucmer` to auto-detect and exclude repetitive regions from variant calling (crucial for accurate bacterial genomics).
* **Robust QC:** `FastP` for cleaning and `Kraken2` for taxonomic contamination checks.
* **Variant Calling:** `FreeBayes` with customizable ploidy (1 or 2) and strict filtering.
* **Backbone Generation:** Creates "All-sites" VCFs (WT + Variants) suitable for phylogenetic tree construction.
* **Lineage & Drug-Resistance Typing:** Alignment-free (k-mer) MTBC lineage + WHO drug-resistance genotyping with **Pathotypr**, bundled in the container and **reference-agnostic** (it types straight from the reads, so calls are correct even when samples are mapped to a non-H37Rv reference).
* **Interactive QC Report:** A single self-contained HTML dashboard (20 panels: QC stats, lineages, genome landscape, SNP dynamics, epistasis, SNP matrix, drug resistance, …) plus a machine-readable per-sample `qc_flags.tsv`, alongside the classic **MultiQC** report.
* **Dual Amino-Acid Numbering:** Optionally re-annotates variants against a canonical reference (H37Rv by default) so the report shows each protein change in both numberings, with links to the Mycobrowser locus.

---

## Workflow Summary

1.  **Reference:** Indexing + Repeat masking + SnpEff DB building.
2.  **QC:** Read validation -> Kraken2 (Taxonomy) -> FastP (Trimming).
3.  **Pathotypr:** (Optional) Alignment-free MTBC lineage (nested sub-lineage) **and** WHO drug-resistance typing from reads.
4.  **Mapping:** `bwa-mem2` alignment -> Merge runs -> Mark Duplicates (`samtools`) -> length-aware read masking (`genmap`).
5.  **Variants:** `FreeBayes` calling + `mpileup` for backbone generation.
6.  **Consensus:** Fasta generation masking low-coverage/low-quality sites.
7.  **Annotation:** `SnpEff` annotation of main and legacy VCFs (+ optional canonical/H37Rv pass for dual amino-acid numbering).
8.  **Report:** Aggregation into the classic **MultiQC** report **and** the interactive **HTML QC report** (`qc_report.py`) with per-sample PASS/WARN/FAIL flags.

---

## 🛠 Prerequisites

* **Nextflow** (`>=24.04.2`)
* **Singularity** or **Docker**
* **Java** (version 11 or later)

The pipeline automatically pulls the container `docker://paururo/bambard:latest`, which contains all necessary tools (BWA, Samtools, FreeBayes, Python, etc.).

---

## Input Format

Create a Tab-Separated Value (TSV) file (e.g., `samples.tsv`) with the following columns.

| Column | Description |
| :--- | :--- |
| `sampleId` | Unique identifier for the sample (e.g., `Sample_A`). |
| `runId` | (Optional) Sequencing run ID. |
| `r1` | Path to Read 1 (FastQ). |
| `r2` | Path to Read 2 (FastQ). Leave empty for Single-End. |
| `refId` | Identifier for the reference genome (e.g., `H37Rv`, `Ecoli_K12`). |
| `refFasta` | Path to the reference FASTA file. |
| `refGff` | Path to the reference GFF file (for SnpEff annotation). |
| `taxId` | (Optional) NCBI TaxID for Kraken2 filtering (e.g., `1773` for TB, `562` for E. coli). |

**Example `samples.tsv`:**
```tsv
sampleId	runId	r1	r2	refId	refFasta	refGff	taxId
T00001	RUN1	/data/T1_R1.fq.gz	/data/T1_R2.fq.gz	H37Rv	/refs/tb.fa	/refs/tb.gff	1773
ECOLI_1	RUN2	/data/EC_R1.fq.gz	/data/EC_R2.fq.gz	K12	/refs/ecoli.fa	/refs/ecoli.gff	562
```
## Usage
By default, the pipeline assumes M. tuberculosis settings (Ploidy=2 to detect mixed infections). Lineage/DR typing is **off** by default — enable it with `--run_pathotypr true`.
```
nextflow run main.nf \
    --tsv samples.tsv \
    --outdir results_bampiro \
    -profile slurm
```
### Handling Multiple Runs per Sample (Merging)

BAMpiro automatically handles multiple sequencing runs (e.g., different lanes or re-sequencing) for the same biological sample. 

* **How to trigger merging:** Simply assign the **same `sampleId`** to multiple rows in your TSV file.
* **The Logic:** The pipeline will process QC and Mapping for each run independently (in parallel) and then **merge** all BAM files associated with that `sampleId` before the Deduplication and Variant Calling steps.

**Example of merging 2 runs into 1 sample:**
```tsv
sampleId    runId   r1                 r2                 refId ...
Sample_A    Run_L1  /data/A_L1_R1.fq   /data/A_L1_R2.fq   H37Rv ...
Sample_A    Run_L2  /data/A_L2_R1.fq   /data/A_L2_R2.fq   H37Rv ...
```
## Configuration Parameters

You can customize the pipeline execution by providing parameters via the command line (e.g., `--threads 16`) or by modifying a config file.

| Category | Parameter | Default | Description |
| :--- | :--- | :--- | :--- |
| **Input/Output** | `--tsv` | `samples_legio.tsv` | Path to the input sample sheet (TSV). |
| | `--outdir` | `results_bampiro` | Directory where results will be saved. |
| | `--threads` | `8` | Max CPUs per process (where applicable). |
| | `--container` | *(pinned digest)* | Container image. Defaults to a pinned `paururo/bambard` digest for reproducibility. |
| | `--nested_output` | `true` | Nest per-sample folders (e.g. `MP001` → `MP/00/1`). |
| | `--publish_mode` | `copy` | `copy` duplicates outputs into `outdir`; `link` hardlinks them to the work dir. |
| | `--output_cram` | `false` | Publish the alignment as CRAM (~40–50% smaller) instead of BAM. |
| **Lineage & DR (Pathotypr)** | `--run_pathotypr` | `false` | Enable alignment-free MTBC lineage + WHO drug-resistance typing. |
| | `--pathotypr_bin` | `pathotypr` | Executable name (on `PATH` inside the container). |
| | `--pathotypr_ref` | `/opt/pathotypr/reference.fasta` | MTBC-ancestor FASTA the markers are defined on (bundled). |
| | `--pathotypr_markers` | `/opt/pathotypr/lineage_markers.tsv` | Zenodo lineage markers (bundled). |
| | `--pathotypr_dr_markers` | `/opt/pathotypr/dr_markers.tsv` | Zenodo WHO drug-resistance markers (bundled). |
| | `--pathotypr_rf_model` | `/opt/pathotypr/rf_model.pathotypr` | Zenodo pre-trained RF lineage model (bundled). |
| | `--pathotypr_min_alt` | `95` | Min alt-allele % for a DR call. Lower it (10–25) to catch heteroresistant / minority alleles. |
| **Dual AA numbering** | `--annotate_canonical` | `false` | Re-annotate variants against a canonical snpEff DB for dual (H37Rv) numbering. |
| | `--canonical_snpeff_db` | `Mycobacterium_tuberculosis_h37rv` | Canonical snpEff genome for the second annotation. |
| | `--canonical_label` | `H37Rv` | How that numbering is labelled in the report. |
| **QC Report** | `--make_qc_report` | `true` | Build the interactive HTML QC report + `qc_flags.tsv`. |
| | `--make_snp_matrix` | `true` | Also emit the master SNP-matrix TSV. |
| | `--report_gate` | `false` | Fail the run if any sample is flagged **FAIL**. |
| | `--report_depth_min` | `10` | Gate: min mean depth (all `report_*` cut-offs are editable live in the report). |
| | `--report_breadth_min` | `90` | Gate: min breadth %. |
| | `--report_missing_max` | `10` | Gate: max missing %. |
| | `--report_dup_max` | `40` | Gate: max duplication %. |
| | `--report_iupac_max` | `2` | Gate: max IUPAC (ambiguous) %. |
| | `--report_mapping_min` | `80` | Gate: min mapped %. |
| **QC & Filter** | `--kraken2_db` | *(path)* | Path to the Kraken2 database directory. |
| | `--fastp_min_length` | `35` | Discard reads shorter than this length. |
| **Mapping/Backbone**| `--allpos_min_cov` | `30` | Minimum coverage to call a site "WT" (otherwise "NC"). |
| | `--allpos_max_depth` | `10000` | Max depth for mpileup to avoid memory issues. |
| | `--allpos_min_bq` | `20` | Minimum base quality for backbone calling. |
| **Variant Calling** | `--freebayes_ploidy` | `2` | Ploidy (1 for haploid, 2 for mixed/diploid). |
| | `--freebayes_min_map_qual`| `30` | Min mapping quality to use a read. |
| | `--freebayes_min_base_qual`| `20` | Min base quality to use a base. |
| | `--freebayes_min_alt_fraction`| `0.05` | Min fraction of alt reads to propose a variant. |
| | `--freebayes_min_alt_count`| `2` | Min alt-supporting reads to propose a variant. |
| **VAF & Filters** | `--hom_threshold` | `0.90` | Frequency ≥ 0.90 is called **Homozygous**. |
| | `--het_min_frac` | `0.10` | Frequency between 0.10 and 0.90 is **Heterozygous**. |
| | `--filter_min_dp` | `30` | Minimum depth required to call a variant. |
| | `--min_alt_fwd/rev` | `2` | Min variant supporting reads in FWD and REV strands. |
| **Consensus** | `--make_consensus` | `true` | Generate a consensus FASTA for each sample. |
| | `--consensus_min_dp` | `7` | Depth threshold below which a base becomes "No Call". |
| | `--consensus_mask_char` | `X` | Character for masked/low-quality sites. |
| | `--consensus_nocall_char`| `-` | Character for no-coverage sites (gaps). |
| **Flags** | `--exclude_repeats` | `true` | Mask self-aligned repetitive regions from reference. |
| | `--annotate_main_vcf` | `true` | Run SnpEff on the main VCF. |
| | `--annotate_legacy_vcfs`| `true` | Run SnpEff on split VCFs (homo/het/indel). |

> **More knobs.** `nextflow.config` also exposes the length-aware read filter (`genmap_*`, `dynamic_read_filter`, `mappability_*`), the reference-bias consensus gate (`consensus_ref_*`), the virgin/unmasked consensus (`keep_virgin_consensus`, `virgin_full_freebayes`), the region-parallel FreeBayes (`freebayes_parallel*`), and the backbone/masking options (`allpos_*`, `mask_baq_dropouts`). See the commented `params { … }` block for the full, annotated list.

## 📊 Interactive QC Report

Every run produces a **single self-contained HTML file** (`<samplesheet>_qc_report.html` — no internet or CDN needed) that folds the whole cohort into one interactive dashboard, plus a machine-readable `<samplesheet>_qc_flags.tsv` of per-sample **PASS/WARN/FAIL** verdicts. It is built by `bin/qc_report.py` and controlled by `--make_qc_report` (default `true`).

**Live & interactive**

* **Live thresholds & presets** — edit any QC cut-off (depth, breadth, missing, duplication, mapping, IUPAC, Ti/Tv, SNP-z, heteroplasmy, mixed-lineage) and the whole report re-flags instantly. Presets: *gate defaults*, *strict (modern WGS)*, *lenient (aDNA / low-cov)*.
* **Collapsible table-of-contents sidebar** with scroll-spy; panels with no data hide themselves (and their nav link).
* **Per-section (i) info popovers** explaining each analysis and its caveats.
* **Exclusion basket** — tick samples (via the table, a drag-box in the scatter, or *basket all flagged*) and export `exclusion.tsv` / `keep_list.txt` for downstream phylogeny.
* Most panels have a fullscreen (⤢) view; the whole report prints / saves to PDF.

**Panels** (grouped as in the sidebar):

| Group | Panels |
| :--- | :--- |
| **Overview** | General statistics (value-coloured, sortable, filterable, TSV export) · Per-lineage summary · Distributions (beeswarm / bar / histogram) |
| **Correlation & structure** | Metric-pair scatter (box-select to basket) · Metric correlation heatmap · QC-space PCA (+ most-unusual-samples table) · Divergence vs completeness |
| **Genome & genes** | Consensus completeness · Genome landscape (per-position callability / variant heatmap, gene search, mask-region toggle) · Functional annotation (snpEff classes) · Functional gene burden · Variable genes (SNP-density hotspots) |
| **Evolution** | Temporal sampling overview · Selection pN/pS (dN/dS, eskaks) · aDNA damage authentication (mapDamage) |
| **Variants over time** | **SNP dynamics** (allele-frequency trajectories over time, per-timepoint DP bars, zoom, series filter) · **Epistasis** (co-varying variant pairs, permutation *p* + BH-FDR *q*, cards / matrix / table views) · **SNP matrix** (site × sample AF matrix, metadata column filter, TSV export) · **Drug resistance** |
| **Quality** | Flagged samples (with the exact failing margin) |

Optional panels appear only when their input is present: SNP dynamics / epistasis / SNP matrix need `--metadata` + the per-sample VCFs; gene burden needs the cohort burden TSV; pN/pS needs an eskaks TSV; drug resistance needs the pathotypr DR calls; aDNA needs mapDamage output.

## 🧬 Lineage & Drug-Resistance Typing (Pathotypr)

`--run_pathotypr true` runs [**Pathotypr**](https://github.com/PathoGenOmics-Lab/pathotypr) — alignment-free (k-mer) MTBC lineage + WHO drug-resistance genotyping — **directly from the container** (a bioconda build; the old hard-coded cluster binary is gone). Because it types straight from the reads using diagnostic k-mers, it is **reference-agnostic**: lineage and DR calls are correct even when your samples were mapped to a non-H37Rv reference.

Each sample gets two `split-fastq` passes:

1. **Lineage** — nested sub-lineage classification against the bundled lineage markers → `<sample>.pathotypr.lineage_summary.tsv` (feeds the per-lineage panel and the lineage colours).
2. **Drug resistance** — WHO-catalogue markers, gated by `--pathotypr_min_alt` (default `95` = near-fixed only; lower to ~10–25 to also catch heteroresistant / minority alleles) → per-sample DR mutations, aggregated into `<samplesheet>_dr.tsv` and shown in the **Drug resistance** panel (a sample × drug matrix with the worst WHO grade per drug, plus a per-mutation table).

The marker panels and pre-trained model (Zenodo v1.0.0, DOI [10.5281/zenodo.19210044](https://doi.org/10.5281/zenodo.19210044)) and the MTBC-ancestor reference are **bundled in the image** under `/opt/pathotypr/`; override any of them with `--pathotypr_ref` / `--pathotypr_markers` / `--pathotypr_dr_markers` / `--pathotypr_rf_model` if you supply your own.

### Dual amino-acid numbering (H37Rv / Mycobrowser)

With `--annotate_canonical true`, each sample's variants are re-annotated with a canonical snpEff database (`--canonical_snpeff_db`, default H37Rv) *in addition to* your mapping reference. The report then shows every amino-acid change in **both** numberings (the used reference + `--canonical_label`, default `H37Rv`) and links each gene to its Mycobrowser locus (`Rv…`). Turn it off for non-MTB organisms, or point `--canonical_snpeff_db` at another snpEff genome. This is only meaningful when the mapping reference shares the canonical reference's coordinates.

## Output Structure
The pipeline organizes results by `sampleId`. Below is a detailed breakdown of the output files using a sample named `MP00091` mapped against reference `LENS`.
Run/cohort-level deliverables land at the top of `outdir` and are prefixed with the sample-sheet name (here `samples_legio`); per-sample files live under their `sampleId` folder.
```text
results_bampiro/
├── samples_legio_qc_report.html          # 🌟 INTERACTIVE consolidated QC report (self-contained HTML)
├── samples_legio_qc_flags.tsv            # Per-sample PASS/WARN/FAIL verdicts
├── samples_legio_summary.tsv             # Cohort metrics table (feeds the report)
├── samples_legio_snp_matrix.tsv          # Master SNP matrix (site × sample AF & depth)
├── samples_legio_gene_burden.tsv         # Per-gene functional burden (cohort)
├── samples_legio_dr.tsv                  # 💊 Drug-resistance calls (pathotypr; only if --run_pathotypr)
│
├── multiqc/
│   └── samples_legio_multiqc_report.html   # 📊 Aggregate Report (QC, Mapping, Variants summary)
│
├── pipeline_info/                        # Nextflow execution timeline / report / trace + software versions
│
├── references/
│   └── LENS/                               # Processed Reference indices & SnpEff DB
│
└── MP00091/                                # 📁 Per-Sample Results Directory
    │
    ├── MP00091.LENS.final.bam              # 🧬 Merged, Coordinate-sorted, Deduplicated BAM
    ├── MP00091.LENS.final.bam.bai          # BAM Index  (or .cram with --output_cram)
    │
    ├── MP00091.LENS.ann.vcf.gz             # 🎯 MAIN OUTPUT: Annotated Variants (SNPs/Indels)
    ├── MP00091.LENS.ann.vcf.gz.tbi         # Index for the main VCF
    ├── MP00091.LENS.canonical.ann.vcf.gz   # 🧬 Variants re-annotated vs the canonical ref (only if --annotate_canonical)
    │
    ├── MP00091.LENS.all.pos.vcf.gz         # 🦴 BACKBONE: VCF containing ALL positions (WT + Variants)
    │                                       # (Ideal for phylogenetic supermatrices)
    │
    ├── MP00091.LENS.consensus.fasta        # 📝 Consensus Sequence (Fasta generated from VCF)
    │
    ├── MP00091.LENS.freebayes.raw...       # 🧪 RAW VCF: Unfiltered calls (debug/comparison)
    ├── MP00091.LENS.var.homo.SNPs.ann...   # 📂 Split VCFs: Subset of Homozygous SNPs (Annotated)
    ├── MP00091.LENS.var.het.SNPs.ann...    # 📂 Split VCFs: Subset of Heterozygous SNPs (Annotated)
    ├── MP00091.LENS.var.homo.indel...      # 📂 SPLIT VCF: Homozygous Indels only
    │
    └── stats/                              # 📉 Statistics & Logs Folder
        ├── MP00091.log                     # -> LEGACY summary log (Tab-separated metrics)
        ├── MP00091.LENS.dedup.stats        # -> Samtools stats (reads mapped, coverage, etc.)
        ├── MP00091.LENS.mask_sites.tsv     # -> Specific positions masked due to low confidence
        ├── MP00091...fastp.html/.json      # -> Trimming quality reports
        ├── MP00091...kraken.report         # -> Taxonomic classification report
        ├── MP00091.LENS.snpeff.csv         # -> Variant effect statistics
        ├── MP00091.pathotypr.lineage_summary.tsv  # -> Pathotypr lineage/sub-lineage call (only if --run_pathotypr)
        ├── MP00091.dr_mutations.tsv        # -> Pathotypr per-sample DR mutations (only if --run_pathotypr)
        └── Locus_to_exclude_LENS.txt       # -> List of repetitive regions excluded from calling
```
## Directory Layout
```
BAMpiro/
├── bin/                     # Python helpers
│   ├── stats_to_legacy.py         # Per-sample metrics -> legacy log (parses the pathotypr lineage call)
│   ├── collect_summary.py         # Aggregate per-sample logs -> cohort summary + gene-burden TSVs
│   ├── collect_dr.py              # Aggregate pathotypr DR calls -> run drug-resistance TSV
│   ├── build_snp_matrix.py        # Build the master SNP matrix (site × sample)
│   ├── qc_report.py               # Build the interactive self-contained HTML QC report
│   ├── WGS_fasta_allpos.py        # Consensus FASTA from the all-positions VCF
│   ├── build_min_unique_len.py    # Per-reference mappability track (genmap)
│   └── filter_reads_mappability.py  # Length-aware read filter
├── assets/
│   └── mycolorsTB_nature.tsv     # Canonical MTBC lineage colour palette (report)
├── modules/                 # Nextflow DSL2 Modules
│   ├── qc.nf                # FastP, Kraken, MultiQC, software versions
│   ├── mapping.nf           # BWA-MEM2, MarkDup, length-aware read filter
│   ├── variants.nf          # FreeBayes, Backbone, Merge
│   ├── annotation.nf        # SnpEff (main / legacy / canonical), legacy stats
│   ├── consensus.nf         # Consensus Fasta
│   ├── report.nf            # Cohort summary, SNP matrix, DR collection, QC report
│   ├── pathotypr.nf         # Lineage + drug-resistance typing
│   ├── reference.nf         # Reference Prep
│   └── utils.nf             # Publish-path routing / clean publish dir
├── .github/dockerfile/      # Container recipe (bundles pathotypr + Zenodo panels + H37Rv snpEff DB)
├── nextflow.config          # Global configuration & params
├── main.nf                  # Main workflow entry point
└── README.md                # This file
```

## 🧛‍♂️ Why "BAMpiro"?
The name is a play on words (a pun) combining bioinformatics and folklore:

- BAM: Stands for Binary Alignment Map. It is the standard file format for storing sequence data aligned to a reference genome. It is the "heart" of this pipeline (mapping -> variant calling).
- Piro: Combined with "BAM", it sounds like "Vampiro" (the Spanish/Portuguese word for Vampire).

The Metaphor: Just as a vampire seeks blood to sustain itself, BAMpiro seeks BAM files (and FASTQ data) to extract vital information (variants, lineages, and stats). It is a "creature" that lives in your cluster and processes bacterial genomes.

---
<h2 id="contributors" align="center">

✨ <a href="https://github.com/PathoGenOmics-Lab/BAMpiro/graphs/contributors">Contributors</a>
</h2>

<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->
<div align="center">
BAMpiro is developed with ❤️ by:
<table>
  <tr>
    <td align="center">
      <a href="https://github.com/paururo">
        <img src="https://avatars.githubusercontent.com/u/50167687?v=4&s=100" width="100px;" alt=""/>
        <br />
        <sub><b>Paula Ruiz-Rodriguez</b></sub>
      </a>
      <br />
      <a href="" title="Code">💻</a>
      <a href="" title="Research">🔬</a>
      <a href="" title="Ideas">🤔</a>
      <a href="" title="Data">🔣</a>
      <a href="" title="Desing">🎨</a>
      <a href="" title="Tool">🔧</a>
    </td> 
    <td align="center">
      <a href="https://github.com/mireiacoscolla">
        <img src="https://avatars.githubusercontent.com/u/29301737?v=4&s=100" width="100px;" alt=""/>
        <br />
        <sub><b>Mireia Coscolla</b></sub>
      </a>
      <br />
      <a href="https://www.uv.es/instituto-biologia-integrativa-sistemas-i2sysbio/es/investigacion/proyectos/proyectos-actuales/mol-tb-host-1286169137294/ProjecteInves.html?id=1286289780236" title="Funding/Grant Finders">🔍</a>
      <a href="" title="Ideas">🤔</a>
      <a href="" title="Mentoring">🧑‍🏫</a>
      <a href="" title="Research">🔬</a>
      <a href="" title="User Testing">📓</a>
    </td> 
  </tr>
</table>

This project follows the [all-contributors](https://github.com/all-contributors/all-contributors) specification ([emoji key](https://allcontributors.org/docs/en/emoji-key)).

<!-- markdownlint-restore -->
<!-- prettier-ignore-end -->

<!-- ALL-CONTRIBUTORS-LIST:END -->
---  
<h2 id="contributors" align="center">
</h2>

<div align="justify">

___
## 📦 Container Specifications (Software Versions)

The Docker container (`paururo/bambard:latest`) includes the following tools:

| Tool | Version | Purpose |
| :--- | :--- | :--- |
| **Nextflow** | `25.10.2` | Workflow management engine |
| **Python** | `3.14.2` | Scripting and orchestration |
| **Java (OpenJDK)** | `23.0.2` | Runtime for Nextflow, SnpEff & FastQC |
| **BWA-MEM2** | `2.3` | High-performance read alignment |
| **Samtools** | `1.23` | BAM/SAM processing and stats |
| **BCFtools** | `1.23` | Variant manipulation and filtering |
| **HTSlib** | `1.23` | C library for high-throughput sequencing data |
| **FreeBayes** | `1.3.10` | Haplotype-based variant caller |
| **SnpEff** | `5.4.0a` | Variant annotation and effect prediction |
| **FastP** | `1.0.1` | Fast all-in-one read pre-processing |
| **Kraken2** | `2.17.1` | Taxonomic classification |
| **Genmap** | *(bioconda)* | Genome mappability (length-aware read filter) |
| **Pathotypr** | *(bioconda)* | Alignment-free MTBC lineage + WHO drug-resistance typing |
| **MultiQC** | `1.33` | Aggregate results reporting |
| **Bedtools** | `2.31.1` | Genome arithmetic |
| **BLAST** | `2.17.0` | Sequence alignment search |
| **MUMmer4** | `4.0.1` | Efficient sequence alignment (used for repeat masking) |
| **Biopython** | `1.86` | Biological computation library |
| **Pandas** | `2.3.3` | Data analysis library |

The image also bundles, under `/opt/pathotypr/`, Pathotypr's marker panels + pre-trained RF model (Zenodo v1.0.0, DOI [10.5281/zenodo.19210044](https://doi.org/10.5281/zenodo.19210044)) and the MTBC-ancestor reference, plus the pre-downloaded **H37Rv snpEff database** (`Mycobacterium_tuberculosis_h37rv`) — so lineage/DR typing and dual amino-acid annotation run offline and reproducibly. `nextflow.config` pins the image to a specific digest; override it with `--container`.
