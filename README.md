

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
* **Lineage/DR SNPs Classification:** Optional integration with **Pathotypr** (Needed custom input files).
* **Reporting:** Generates statistic logs and a dynamic **MultiQC** report.

---

## Workflow Summary

1.  **Reference:** Indexing + Repeat masking + SnpEff DB building.
2.  **QC:** Read validation -> Kraken2 (Taxonomy) -> FastP (Trimming).
3.  **Pathotypr:** (Optional) Rapid lineage classification from reads.
4.  **Mapping:** `bwa-mem2` alignment -> Merge runs -> Mark Duplicates (`samtools`).
5.  **Variants:** `FreeBayes` calling + `mpileup` for backbone generation.
6.  **Consensus:** Fasta generation masking low-coverage/low-quality sites.
7.  **Annotation:** `SnpEff` annotation of main and legacy VCFs.
8.  **Stats:** Aggregation of all metrics into a single MultiQC report.

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
By default, the pipeline assumes M. tuberculosis settings (Ploidy=2 to detect mixed infections, Pathotypr enabled).
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
| **Pathotypr** | `--run_pathotypr` | `false` | Set to `true` to enable lineage classification. |
| | `--pathotypr_bin` | *(path)* | Path to the Pathotypr binary executable. |
| | `--pathotypr_markers` | *(path)* | Path to the lineage markers file. |
| | `--pathotypr_ref` | *(path)* | Reference fasta used for Pathotypr. |
| **QC & Filter** | `--kraken2_db` | *(path)* | Path to the Kraken2 database directory. |
| | `--fastp_min_length` | `35` | Discard reads shorter than this length. |
| **Mapping/Backbone**| `--allpos_min_cov` | `30` | Minimum coverage to call a site "WT" (otherwise "NC"). |
| | `--allpos_max_depth` | `10000` | Max depth for mpileup to avoid memory issues. |
| | `--allpos_min_bq` | `20` | Minimum base quality for backbone calling. |
| **Variant Calling** | `--freebayes_ploidy` | `2` | Ploidy (1 for haploid, 2 for mixed/diploid). |
| | `--freebayes_min_map_qual`| `30` | Min mapping quality to use a read. |
| | `--freebayes_min_base_qual`| `20` | Min base quality to use a base. |
| | `--freebayes_min_alt_frac`| `0.05` | Min fraction of alt reads to propose a variant. |
| **VAF & Filters** | `--hom_threshold` | `0.90` | Frequency ≥ 0.90 is called **Homozygous**. |
| | `--het_min_frac` | `0.10` | Frequency between 0.10 and 0.90 is **Heterozygous**. |
| | `--filter_min_dp` | `30` | Minimum depth required to call a variant. |
| | `--min_alt_fwd/rev` | `2` | Min variant supporting reads in FWD and REV strands. |
| **Consensus** | `--make_consensus` | `true` | Generate a consensus FASTA for each sample. |
| | `--consensus_min_dp` | `7` | Depth threshold below which a base becomes "No Call". |
| | `--consensus_mask_char` | `X` | Character for masked/low-quality sites. |
| | `--consensus_nocall_char`| `-` | Character for no-coverage sites (gaps). |
| **Flags** | `--exclude_repeats` | `true` | Mask self-aligned repetitive regions from reference. |
| | `--annotate_legacy_vcfs`| `true` | Run SnpEff on split VCFs (homo/het/indel). |

## Output Structure
The pipeline organizes results by `sampleId`. Below is a detailed breakdown of the output files using a sample named `MP00091` mapped against reference `LENS`.
```text
results_bampiro/
├── multiqc/
│   └── samples_legio_multiqc_report.html   # 📊 Aggregate Report (QC, Mapping, Variants summary)
│
├── references/
│   └── LENS/                               # Processed Reference indices & SnpEff DB
│
└── MP00091/                                # 📁 Per-Sample Results Directory
    │
    ├── MP00091.LENS.final.bam              # 🧬 Merged, Coordinate-sorted, Deduplicated BAM
    ├── MP00091.LENS.final.bam.bai          # BAM Index
    │
    ├── MP00091.LENS.ann.vcf.gz             # 🎯 MAIN OUTPUT: Annotated Variants (SNPs/Indels)
    ├── MP00091.LENS.ann.vcf.gz.tbi         # Index for the main VCF
    │
    ├── MP00091.LENS.all.pos.vcf.gz         # 🦴 BACKBONE: VCF containing ALL positions (WT + Variants)
    │                                       # (Ideal for phylogenetic supermatrices)
    │
    ├── MP00091.LENS.consensus.fasta        # 📝 Consensus Sequence (Fasta generated from VCF)
    │
    ├── MP00091.LENS.var.homo.SNPs.ann...   # 📂 Split VCFs: Subset of Homozygous SNPs (Annotated)
    ├── MP00091.LENS.var.het.SNPs.ann...    # 📂 Split VCFs: Subset of Heterozygous SNPs (Annotated)
    │
    └── stats/                              # 📉 Statistics & Logs Folder
        ├── MP00091.log                     # -> LEGACY summary log (Tab-separated metrics)
        ├── MP00091.LENS.dedup.stats        # -> Samtools stats (reads mapped, coverage, etc.)
        ├── MP00091...fastp.html/.json      # -> Trimming quality reports
        ├── MP00091...kraken.report         # -> Taxonomic classification report
        ├── MP00091.LENS.snpeff.csv         # -> Variant effect statistics
        └── Locus_to_exclude_LENS.txt       # -> List of repetitive regions excluded from calling
```
## Directory Layout
```
BAMpiro/
├── bin/                     # Python scripts (stats_to_legacy.py, WGS_fasta_allpos.py)
├── modules/                 # Nextflow DSL2 Modules
│   ├── qc.nf                # FastP, Kraken, MultiQC
│   ├── mapping.nf           # BWA, MarkDup
│   ├── variants.nf          # FreeBayes, Backbone, Merge
│   ├── annotation.nf        # SnpEff, Stats Legacy
│   ├── consensus.nf         # Consensus Fasta
│   ├── pathotypr.nf         # Pathotypr Logic
│   └── reference.nf         # Reference Prep
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

✨ [Contributors]((https://github.com/PathoGenOmics-Lab/AMAP/graphs/contributors))
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
