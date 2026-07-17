# Installation

## Prerequisites

- **Nextflow** (`>=24.04.2`)
- **Singularity** or **Docker**
- **Java** (version 17 or later, as required by Nextflow >= 24.04.2)

The pipeline automatically pulls the container `docker://paururo/bampiro:1.0.1`,
which contains every tool it needs (BWA-MEM2, Samtools, FreeBayes, SnpEff,
Pathotypr, Python, …). `nextflow.config` references the image by its mutable `1.0.1`
tag; for byte-for-byte reproducibility re-pin it to the digest the image was built
from (`docker://paururo/bampiro@sha256:<digest>`). Override either with `--container`.

```bash
nextflow run main.nf --tsv samples.tsv --outdir results_bampiro -profile standard
```

See [Quick Start](quickstart.md) for a full example and the
[samplesheet format](quickstart.md#samplesheet).

## Running without SLURM

The **default** (`-profile standard`) targets the authors' SLURM cluster — it uses the `slurm` executor,
`module load singularity`, and site-specific Singularity `--bind` paths. On any other machine those fail
before the pipeline runs, so use the `local` (and, for Docker, `docker`) profiles, which drop all of that:

```bash
# Singularity/Apptainer on a single machine (no SLURM, no environment-modules):
nextflow run main.nf --tsv samples.tsv --outdir results -profile local

# Docker on a laptop:
nextflow run main.nf --tsv samples.tsv --outdir results -profile local,docker
```

`local` runs every task on the current host with the executor set to `local` and the cluster-only bind
paths cleared. Point `--kraken2_db` at your own Kraken2 database — if it lives outside the launch directory,
Singularity may need it bound explicitly, e.g. `-profile local` plus
`--kraken2_db /data/kraken2` and `export NXF_SINGULARITY_RUN_OPTIONS="--bind /data/kraken2"`.

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
| **MUMmer4** | `4.0.1` | `nucmer` / `show-coords` self-alignment for repeat masking |
| **Biopython** | `1.86` | Biological computation library |
| **Pandas** | `2.3.3` | Data analysis library |

The image also bundles, under `/opt/pathotypr/`, Pathotypr's marker panels +
pre-trained RF model (Zenodo v1.0.0, DOI
[10.5281/zenodo.19210044](https://doi.org/10.5281/zenodo.19210044)) and the
MTBC-ancestor reference, plus the pre-downloaded **H37Rv snpEff database**
(`Mycobacterium_tuberculosis_h37rv`) - so [lineage/DR typing](pathotypr.md) and
[dual amino-acid annotation](pathotypr.md#dual-amino-acid-numbering-h37rv--mycobrowser)
run offline and reproducibly.
