

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

## Output Structure
After a successful run, the results_bampiro/ directory will look like this:
```
results_bampiro/
├── multiqc/
│   ├── samples_multiqc_report.html   # Final interactive Summary (FastP, Kraken, Alignment, Variants)
│   └── samples_multiqc_report_data/  # Raw data for the report
├── references/
│   └── H37Rv/                        # Processed reference indices & SnpEff DB
├── T00001/                           # Per-sample results
│   ├── stats/
│   │   └── T00001.log                # Legacy statistics log (tab-separated)
│   ├── lineage/
│   │   └── T00001.pathotypr...       # Lineage classification results (if enabled)
│   ├── T00001.H37Rv.final.bam        # Final Deduped BAM
│   ├── T00001.H37Rv.vcf.gz           # Main Variant VCF (Annotated)
│   ├── T00001.H37Rv.all.pos.vcf.gz   # All-sites VCF (Backbone)
│   └── T00001.H37Rv.consensus.fasta  # Consensus Sequence
└── ...
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
