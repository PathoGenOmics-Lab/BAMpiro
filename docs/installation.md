# Installation

## Prerequisites

- **Nextflow** (`>=24.04.2`)
- **Singularity** or **Docker**
- **Java** (version 17 or later, as required by Nextflow >= 24.04.2)

The pipeline automatically pulls the container `docker://paururo/bampiro:1.0.1`,
which contains every tool it needs (BWA-MEM2, Samtools, FreeBayes, SnpEff,
Pathotypr, Python, …). `nextflow.config` pins the image to a specific digest for
reproducibility; override it with `--container`.

```bash
nextflow run main.nf --tsv samples.tsv --outdir results_bampiro -profile slurm
```

See [Quick Start](quickstart.md) for a full example and the
[samplesheet format](quickstart.md#samplesheet).

## Container contents (software versions)

The Docker container (`paururo/bampiro:1.0.1`) bundles the following tools:

| Tool | Version | Purpose |
| :--- | :--- | :--- |
| **Nextflow** | `25.10.2` | Workflow management engine |
| **Python** | `3.14.2` | Scripting and orchestration |
| **Java (OpenJDK)** | `23.0.2` | Runtime for Nextflow & SnpEff |
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

The image also bundles, under `/opt/pathotypr/`, Pathotypr's marker panels +
pre-trained RF model (Zenodo v1.0.0, DOI
[10.5281/zenodo.19210044](https://doi.org/10.5281/zenodo.19210044)) and the
MTBC-ancestor reference, plus the pre-downloaded **H37Rv snpEff database**
(`Mycobacterium_tuberculosis_h37rv`) - so [lineage/DR typing](pathotypr.md) and
[dual amino-acid annotation](pathotypr.md#dual-amino-acid-numbering-h37rv--mycobrowser)
run offline and reproducibly.
