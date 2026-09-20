# Installation

## Prerequisites

- **Nextflow** (`>=24.04.2`)
- **Singularity** or **Docker**
- **Java** (version 17 or later, as required by Nextflow >= 24.04.2)

The pipeline is parsed by both the classic and the strict Nextflow parser, so it runs
unchanged on 24.04.2 and on 25.10 or newer without setting `NXF_SYNTAX_PARSER=v1`.

A run needs two things: a samplesheet (`--tsv`, required) and a profile that says where to
execute (`-profile`, defaulting to the current host).

```bash
nextflow run main.nf --tsv samples.tsv --outdir results_bampiro -profile local,docker
```

See [Quick Start](quickstart.md) for a full example and the
[samplesheet format](quickstart.md#samplesheet). `nextflow run main.nf --help` lists every
parameter with its default and exits.

## The container

The pipeline pulls `paururo/bampiro` automatically; it contains every tool it needs
(BWA-MEM2, Samtools, FreeBayes, SnpEff, Pathotypr, Python, …). `nextflow.config` pins the
image **by digest**, not by a tag:

```groovy
container = "docker://paururo/bampiro@sha256:78355aa4909dad87c44e2157ccb8032c86d24da3a8e10f37b5075a2b9ba57626"
```

A tag can be repointed at a rebuilt image, so two runs a month apart could silently use
different tool versions; a digest cannot. That digest is the 1.1.0-rc1 image. Override it with
`--container`; either way the value actually used is recorded in the report's provenance
footer.

The same image is published to Docker Hub (`paururo/bampiro`, public) and to GHCR
(`ghcr.io/pathogenomics-lab/bampiro`, currently private) with the same digest.

## Choosing where it runs

`-profile` decides **where** each task is executed. Nothing site-specific is baked into
the defaults: with no `-profile` Nextflow applies `standard`, which runs every task on the
current host, so a fresh clone works anywhere without being edited.

| Profile | Executor | Use it for |
| :--- | :--- | :--- |
| `standard` *(applied when you pass no `-profile`)* | `local` | Whatever machine you launched from |
| `local` | `local` | The same thing, spelled out in the command |
| `slurm` | `slurm` | Any SLURM cluster. No queue, account or module system assumed |
| `garnatxa` | `slurm` | The I2SysBio / University of Valencia cluster, configured end to end |
| `docker` | *(leaves it alone)* | Docker instead of Singularity (combine it: `local,docker`) |
| `generic` | *(leaves it alone)* | A non-tuberculosis organism (ploidy 1, the MTBC-only features off) |
| `test` · `test_full` | `local` | The bundled fixture cohort, for `-stub-run` |

Pick one executor profile and add as many option profiles as you like, comma-separated
and with no spaces: `-profile local,docker`, `-profile slurm,generic`.

!!! warning "On a cluster, choose the executor profile deliberately"

    The default runs on the machine you launched from. On a cluster login node that means
    the entire pipeline runs **on the login node**: nothing is ever submitted to the
    scheduler, and your admins will notice. Pass `-profile slurm` (or `-profile garnatxa`)
    for real work. The startup banner prints the profile, executor, container and Kraken
    state it resolved, so read the `Executor` line before walking away.

=== "Local (Singularity/Apptainer)"

    On a single machine, no scheduler and no environment modules:

    ```bash
    nextflow run main.nf --tsv samples.tsv --outdir results -profile local
    ```

=== "Docker (laptop)"

    ```bash
    nextflow run main.nf --tsv samples.tsv --outdir results -profile local,docker
    ```

=== "SLURM cluster"

    ```bash
    nextflow run main.nf --tsv samples.tsv --outdir results -profile slurm
    ```

    `slurm` submits tasks with `sbatch` and assumes nothing else: no queue, no account, no
    environment modules. Add your site's QoS / partition and your Singularity bind paths in
    a small `-c` config, as shown in
    [Configuring a run](tutorials/configuring-a-run.md#on-an-hpc-cluster-make-your-files-visible-to-the-container).

=== "Garnatxa (I2SysBio)"

    ```bash
    nextflow run main.nf --tsv samples.tsv --outdir results -profile garnatxa
    ```

    Everything that cluster needs is already in `conf/garnatxa.config`: the SLURM QoS per
    step, `module load singularity`, the Singularity bind paths, the shared Kraken2
    database, a reusable mappability cache and a work directory on scratch. It also works
    when the pipeline is pulled straight from GitHub
    (`nextflow run PathoGenOmics-Lab/BAMpiro -profile garnatxa`), so there is nothing to
    write or copy first.

    Submit the driver rather than running it on the login node, which is what that cluster's
    own documentation asks for:

    ```bash
    sbatch conf/garnatxa.sbatch samples.tsv results
    ```

    `conf/garnatxa.sbatch` asks for one core and 4 GB for Nextflow itself, which submits every
    task as its own job and waits. Its `--time` is set explicitly and deliberately: **every QoS
    on Garnatxa defaults to six hours**, and a driver killed at its limit orphans whatever it
    was waiting on.

    !!! tip "The Kraken2 database is already set there"

        `-profile garnatxa` points `--kraken2_db` at the shared copy under
        `/storage/shared_datasets`, so it needs no flag. What it does need is a `taxId` column
        in the samplesheet, or there is nothing to filter against. Kraken2 asks for about 80 GB,
        which is the memory ceiling of the whole run.

## Kraken2 is opt-in

`--kraken2_db` has no default: without one the contamination screen is skipped and the
report's taxonomy panel hides itself. To run it, point `--kraken2_db` at your own database
(and give the samplesheet a `taxId` column). If the database lives outside the launch
directory, Singularity needs it bound explicitly, e.g. `--kraken2_db /data/kraken2` plus
`export NXF_SINGULARITY_RUN_OPTIONS="--bind /data/kraken2"`. `-profile garnatxa` already
sets both the database and its bind.

## Container contents (software versions)

The container (`paururo/bampiro`, the 1.1.0-rc1 image pinned above) bundles the following tools:

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

## Resource requirements

Nextflow schedules each process with its own CPU / RAM request (retrying with more RAM
on an out-of-memory kill). The peak driver is **Kraken2**, so size your machine / queue
for it whenever you set `--kraken2_db`:

| Step | CPUs | Memory | Notes |
| :--- | :--- | :--- | :--- |
| **Kraken2** (contamination) | - | **~80 GB** | loads the whole DB; the run's memory ceiling. `--kraken_memory_mapping` (default on) mmaps it so parallel tasks share RAM |
| Mappability track (`genmap`) | 8 | 16 GB | once per reference (cached across runs) |
| Reference prep / SnpEff DB build | 4 | 8 GB | once per reference |
| Mapping (`bwa-mem2`) | `--threads` (8) | scales with genome | the only step `--threads` controls |
| Pathotypr typing | 4 | 8 GB | only with `--run_pathotypr` |
| Variant calling / consensus / report | 1–4 | 2–4 GB | - |

!!! warning "Kraken2 sets the memory ceiling (~80 GB)"

    Kraken2 loads the entire database into RAM, so it is the single biggest memory
    request in the run - an under-sized host kills it with an **OOM (exit 137)**. On a
    laptop, use a smaller Kraken2 DB, or leave `--kraken2_db` unset and the whole screen
    is skipped. Disk is roughly a few GB per sample (BAM/CRAM + VCFs + consensus).

    To make the whole pipeline fit a small machine, cap every request at once with
    `--max_cpus` / `--max_memory` / `--max_time` (e.g. `--max_cpus 4 --max_memory 8.GB`).
    Uncapped is the right choice on a cluster, where the per-process numbers are the point.
