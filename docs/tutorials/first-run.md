# Your first run

This tutorial takes you from an empty directory to an open, interactive QC report in a single sitting. Budget roughly **15–20 minutes** of hands-on time, plus a one-off container pull the first time you run.

By the end you'll have launched BAMpiro on a tiny samplesheet and opened the self-contained HTML report it produces.

## Before you start

You only need two things on your machine - the container carries every bioinformatics tool (BWA-MEM2, Samtools, FreeBayes, SnpEff, Pathotypr, …), so there is nothing else to install.

- **Nextflow** `>= 24.04.2` (which needs **Java 17+**)
- **Docker** *or* **Singularity / Apptainer**

!!! tip "Installing Nextflow in a conda environment"

    The cleanest way to get Nextflow - **and the Java it needs** - without touching your system is a
    dedicated **conda / mamba** environment:

    ```bash
    conda create -n bampiro -c conda-forge -c bioconda nextflow
    conda activate bampiro
    nextflow -version   # confirm it is >= 24.04.2
    ```

    Bioconda pulls in a compatible Java for you. Activate this environment whenever you run BAMpiro.
    (The official installer `curl -s https://get.nextflow.io | bash` also works, but you must already
    have Java 17+ on your `PATH`.)

!!! note "The container does the heavy lifting"

    On the first run, BAMpiro automatically pulls its image with every tool built in.
    The image is pinned **by digest** (`docker://paururo/bampiro@sha256:bef4375…`), not by
    a tag, so the run you do today and the one you repeat next year use the same binaries.
    That download is a few GB and happens once - later runs reuse the cached image. Full
    requirements and the bundled software versions are in
    [Installation](../installation.md#the-container).

## 1. Get BAMpiro

Clone the repository so you have `main.nf` and the bundled configs to hand:

```bash
git clone https://github.com/PathoGenOmics-Lab/BAMpiro
cd BAMpiro
```

## 2. Prepare a tiny samplesheet

The only required input is a Tab-Separated samplesheet passed with `--tsv` - there is no default. For your first run, keep it to a single sample. Each row points at your reads and a reference genome + GFF:

```tsv
sampleId	runId	r1	r2	refId	refFasta	refGff	taxId
T00001	RUN1	/data/T1_R1.fq.gz	/data/T1_R2.fq.gz	H37Rv	/refs/tb.fa	/refs/tb.gff	1773
```

Reads, `refFasta` and `refGff` may be plain or gzipped. Don't worry about the finer points yet - the full column reference and multi-run merging are covered in [The samplesheet](samplesheet.md).

!!! tip "Just want to look around first?"

    A bundled **17-sample MTBC demo cohort** and a runnable Jupyter notebook let you
    explore a finished run's outputs with `pandas` / `matplotlib` - no pipeline run
    needed. See the [Quick Start](../quickstart.md) for that shortcut.

## 3. Launch the run

Run the pipeline, pointing `--tsv` at your samplesheet and `--outdir` at where results should land. Pick the tab for your container engine - both use `-profile local`, which runs everything on the current host:

=== "Docker"

    ```bash
    nextflow run main.nf \
        --tsv samples.tsv \
        --outdir results_bampiro \
        -profile local,docker
    ```

=== "Singularity / Apptainer"

    ```bash
    nextflow run main.nf \
        --tsv samples.tsv \
        --outdir results_bampiro \
        -profile local
    ```

!!! warning "On a cluster, say where you want the work to run"

    `local` is also what you get with no `-profile` at all, because the default profile
    (`standard`) runs on the machine you launched from. That is what makes a fresh clone
    work anywhere, but it cuts both ways: launch it from a cluster **login node** and the
    whole pipeline runs on the login node instead of being submitted to the scheduler.

    On a cluster, choose the executor deliberately: `-profile slurm` for a generic SLURM
    site, or `-profile garnatxa` on the I2SysBio cluster. The startup banner prints the
    `Executor` it resolved, so you can confirm it in the first seconds of a run. See
    [Choosing where it runs](../installation.md#choosing-where-it-runs).

    On **your own cluster** you will also need to bind your files into the Singularity
    container so tasks can read them - see
    [On an HPC cluster](configuring-a-run.md#on-an-hpc-cluster-make-your-files-visible-to-the-container).

Nextflow prints a live table of processes as they run. Behind that table it QCs and trims your reads, screens for contamination, aligns to the reference, masks repeats and low-mappability regions, calls and annotates variants, builds a consensus, and folds everything into one report.

If a step goes red, the [Troubleshooting](../troubleshooting.md) page covers the common first-run errors (a frequent one on small machines is Kraken2 running out of memory).

## 4. Open the QC report

When the run finishes, everything is under your `--outdir`. The one file to open is the consolidated, self-contained report:

```
results_bampiro/<samplesheet>_qc_report.html
```

Open it in any browser - double-click it, or:

```bash
open results_bampiro/*_qc_report.html   # macOS  (Linux: xdg-open)
```

It needs no internet and no server: every panel, threshold and theme toggle is baked into that single HTML file, so it works even on an HPC login node. The full result tree is described in [Outputs](../outputs.md).

## What's next

You've completed an end-to-end run. Next, learn how to describe **any** cohort - multiple samples, multiple runs, single-end reads and several references - in [**The samplesheet**](samplesheet.md).
