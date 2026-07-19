# Your first run

This tutorial takes you from an empty directory to an open, interactive QC report in a single sitting. Budget roughly **15–20 minutes** of hands-on time, plus a one-off container pull the first time you run.

By the end you'll have launched BAMpiro on a tiny samplesheet and opened the self-contained HTML report it produces.

## Before you start

You only need two things on your machine — the container carries every bioinformatics tool (BWA-MEM2, Samtools, FreeBayes, SnpEff, Pathotypr, …), so there is nothing else to install.

- **Nextflow** `>= 24.04.2` (which needs **Java 17+**)
- **Docker** *or* **Singularity / Apptainer**

!!! note "The container does the heavy lifting"

    On the first run, BAMpiro automatically pulls its pinned image
    (`docker://paururo/bampiro:1.0.1`) with every tool built in. That download is a
    few GB and happens once — later runs reuse the cached image. Full requirements
    and the bundled software versions are in [Installation](../installation.md).

## 1. Get BAMpiro

Clone the repository so you have `main.nf` and the bundled configs to hand:

```bash
git clone https://github.com/PathoGenOmics-Lab/BAMpiro
cd BAMpiro
```

## 2. Prepare a tiny samplesheet

The only required input is a Tab-Separated samplesheet passed with `--tsv` — there is no default. For your first run, keep it to a single sample. Each row points at your reads and a reference genome + GFF:

```tsv
sampleId	runId	r1	r2	refId	refFasta	refGff	taxId
T00001	RUN1	/data/T1_R1.fq.gz	/data/T1_R2.fq.gz	H37Rv	/refs/tb.fa	/refs/tb.gff	1773
```

Reads, `refFasta` and `refGff` may be plain or gzipped. Don't worry about the finer points yet — the full column reference and multi-run merging are covered in [The samplesheet](samplesheet.md).

!!! tip "Just want to look around first?"

    A bundled **17-sample MTBC demo cohort** and a runnable Jupyter notebook let you
    explore a finished run's outputs with `pandas` / `matplotlib` — no pipeline run
    needed. See the [Quick Start](../quickstart.md) for that shortcut.

## 3. Launch the run

Run the pipeline, pointing `--tsv` at your samplesheet and `--outdir` at where results should land. Pick the tab for your container engine — both use `-profile local`, which runs everything on the current host:

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

!!! warning "`-profile standard` targets a SLURM cluster"

    The default `standard` profile (as shown in the [Quick Start](../quickstart.md))
    uses the `slurm` executor and site-specific paths. On a laptop or any non-SLURM
    host, use `local` (or `local,docker`) as above, or the run fails before it
    starts. See [Running without SLURM](../installation.md#running-without-slurm).

    Running on **your own cluster**? You'll also need to bind your files into the
    Singularity container so tasks can read them — see
    [On an HPC cluster](configuring-a-run.md#on-an-hpc-cluster-make-your-files-visible-to-the-container).

Nextflow prints a live table of processes as they run. Behind that table it QCs and trims your reads, screens for contamination, aligns to the reference, masks repeats and low-mappability regions, calls and annotates variants, builds a consensus, and folds everything into one report.

If a step goes red, the [Troubleshooting](../troubleshooting.md) page covers the common first-run errors (a frequent one on small machines is Kraken2 running out of memory).

## 4. Open the QC report

When the run finishes, everything is under your `--outdir`. The one file to open is the consolidated, self-contained report:

```
results_bampiro/<samplesheet>_qc_report.html
```

Open it in any browser — double-click it, or:

```bash
open results_bampiro/*_qc_report.html   # macOS  (Linux: xdg-open)
```

It needs no internet and no server: every panel, threshold and theme toggle is baked into that single HTML file, so it works even on an HPC login node. The full result tree is described in [Outputs](../outputs.md).

## What's next

You've completed an end-to-end run. Next, learn how to describe **any** cohort — multiple samples, multiple runs, single-end reads and several references — in [**The samplesheet**](samplesheet.md).
