# Configuring a run

By now you've completed a run and built a samplesheet. This page shows you how to point BAMpiro at *your* machine and *your* organism: picking an execution profile, setting parameters, and switching on the optional MTBC features. Budget about 15 minutes to read and adapt to your setup.

For the exhaustive parameter list (every flag, its default and description), keep [Configuration](../configuration.md) open alongside this page — here we cover only the choices that matter most for a first tuning.

## 1. Choose an execution profile

A `-profile` tells Nextflow *where* and *how* to run each task. The default (`standard`) targets the authors' SLURM cluster, so on any other machine you'll want a different one.

=== "Laptop / workstation (Docker)"

    ```bash
    nextflow run main.nf --tsv samples.tsv --outdir results -profile local,docker
    ```

=== "Single machine (Singularity/Apptainer)"

    ```bash
    nextflow run main.nf --tsv samples.tsv --outdir results -profile local
    ```

=== "SLURM cluster (default)"

    ```bash
    nextflow run main.nf --tsv samples.tsv --outdir results -profile standard
    ```

`local` runs every task on the current host and clears the cluster-only `module load` and Singularity bind paths; add `docker` if you'd rather use Docker than Singularity. Profiles compose with a comma (e.g. `local,docker`).

!!! warning "Kraken2 sets the memory ceiling (~80 GB)"

    Kraken2 loads its whole database into RAM and is the single biggest request in the run — an under-sized host kills it with an OOM (exit 137). On a laptop, point `--kraken2_db` at a smaller database. See [Resource requirements](../installation.md#resource-requirements).

## 2. Set parameters

Any parameter can be set on the command line with a `--` prefix, for example:

```bash
nextflow run main.nf --tsv samples.tsv --outdir results \
  -profile local,docker --threads 16 --freebayes_ploidy 1
```

For anything you want to reuse, collect your overrides in a config file and pass it with `-c` — command-line `--flags` still win over the file:

```groovy title="my_run.config"
params {
    threads         = 16
    run_pathotypr   = true
    report_gate     = true
}
```

```bash
nextflow run main.nf --tsv samples.tsv --outdir results -profile local,docker -c my_run.config
```

!!! tip "Single dash vs double dash"

    `-profile` and `-c` are **Nextflow** options (one dash). `--tsv`, `--threads`, `--run_pathotypr` and friends are **pipeline** parameters (two dashes).

## 3. Run against another organism

BAMpiro is organism-agnostic; only its *defaults* are TB-tuned (diploid calling for mixed infections, plus the MTBC-only features). To type a different species, point the samplesheet's `refFasta` / `refGff` columns at your genome, set `--kraken2_db` and the `taxId` column, and add the built-in **`generic`** profile:

```bash
nextflow run main.nf --tsv samples.tsv --outdir results -profile local,docker,generic
```

The `generic` profile sets `freebayes_ploidy 1` (a clonal bacterium is haploid) and forces the *M. tuberculosis*-specific features off. The report's TB-only panels (lineage, drug resistance) simply self-hide when there's no such data. The commented template [`conf/organism.config`](https://github.com/PathoGenOmics-Lab/BAMpiro/blob/indel-mask/conf/organism.config) covers the rest.

## 4. Turn on the MTBC features

For a TB cohort like the bundled 17-sample demo, the built-in defaults already give you a valid run. To add the full recommended MTBC feature set, switch on lineage/DR typing and dual amino-acid numbering:

```bash
nextflow run main.nf --tsv samples.tsv --outdir results -profile local,docker \
  --run_pathotypr true --annotate_canonical true
```

- `--run_pathotypr true` enables alignment-free MTBC lineage + WHO drug-resistance typing ([Lineage & Drug-Resistance Typing](../pathotypr.md)).
- `--annotate_canonical true` adds dual **amino-acid** numbering — re-annotating variants against H37Rv, with Mycobrowser gene links.

!!! tip "Or use the ready-made TB config"

    [`conf/tuberculosis.config`](https://github.com/PathoGenOmics-Lab/BAMpiro/blob/indel-mask/conf/tuberculosis.config) bundles those two plus the H37Rv coordinate liftover (`--variant_liftover`) and blind-spot masking (`--mask_blindspots`). Apply it with `-c conf/tuberculosis.config`.

## 5. A word on the QC gate

Two things are easy to confuse. The interactive QC report lets you drag the cut-offs and re-flags samples **live** in your browser — that never changes any file. The pipeline's own pass/fail **gate** is separate: it uses the `--report_*` thresholds from your config, and only stops the run when you opt in with `--report_gate true`.

```bash
# Fail the whole run if any sample is flagged FAIL, using a stricter depth cut-off
--report_gate true --report_depth_min 20
```

The gate thresholds (`--report_depth_min`, `--report_breadth_min`, `--report_missing_max`, and so on) are listed in full under [Configuration](../configuration.md). You'll see exactly how each one flags a sample in the next tutorial.

## What's next

You've tuned a run for your environment and organism. Next, learn to read what it produced → [Reading the QC report](reading-the-qc-report.md).
