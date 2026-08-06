# Configuring a run

By now you've completed a run and built a samplesheet. This page shows you how to point BAMpiro at *your* machine and *your* organism: picking an execution profile, setting parameters, and switching on the optional MTBC features. Budget about 15 minutes to read and adapt to your setup.

For the exhaustive parameter list (every flag, its default and description), keep [Configuration](../configuration.md) open alongside this page - here we cover only the choices that matter most for a first tuning.

## 1. Choose an execution profile

A `-profile` tells Nextflow *where* and *how* to run each task. Nothing site-specific is baked into the defaults, so the choice is entirely yours: with no `-profile`, Nextflow applies `standard`, which runs everything on the machine you launched from.

=== "Laptop / workstation (Docker)"

    ```bash
    nextflow run main.nf --tsv samples.tsv --outdir results -profile local,docker
    ```

=== "Single machine (Singularity/Apptainer)"

    ```bash
    nextflow run main.nf --tsv samples.tsv --outdir results -profile local
    ```

=== "SLURM cluster"

    ```bash
    nextflow run main.nf --tsv samples.tsv --outdir results -profile slurm
    ```

    Submits every task with `sbatch`. It assumes nothing about your site: no queue, no account, no environment modules. Add those in a small `-c` config, [as shown below](#on-an-hpc-cluster-make-your-files-visible-to-the-container).

=== "Garnatxa (I2SysBio)"

    ```bash
    nextflow run main.nf --tsv samples.tsv --outdir results -profile garnatxa
    ```

    That cluster's QoS, module load, bind paths, Kraken2 database and scratch work directory are already set. Nothing else to write.

`local` runs every task on the current host; add `docker` if you'd rather use Docker than Singularity. Profiles compose with a comma and no spaces (e.g. `local,docker`, `slurm,generic`).

!!! warning "The default runs where you launched it"

    Because `standard` means "this host", launching from a cluster **login node** with no `-profile` runs the entire pipeline on the login node: nothing is ever submitted to the scheduler, and it will be noticed. Pick `slurm` or `garnatxa` deliberately for real work.

    You do not have to guess which one you got. The startup banner prints the resolved profile, executor, container and Kraken state:

    ```text
    Profile(s)       : slurm
    Executor         : slurm
    Container        : docker://paururo/bampiro@sha256:c3bc851…
    Kraken2          : disabled
    ```

    `Kraken2` reads `disabled` until you point `--kraken2_db` at a database, which has no default. While it is disabled the contamination screen is skipped and the report's taxonomy panel hides itself.

!!! warning "Kraken2 sets the memory ceiling (~80 GB)"

    Kraken2 loads its whole database into RAM and is the single biggest request in the run - an under-sized host kills it with an OOM (exit 137). On a laptop, point `--kraken2_db` at a smaller database, or leave it unset and the screen is skipped entirely. To make every step fit a small machine in one go, cap the requests with `--max_cpus 4 --max_memory 8.GB`. See [Resource requirements](../installation.md#resource-requirements).

### On an HPC cluster: make your files visible to the container

On a cluster BAMpiro runs each task through **Singularity / Apptainer**, and a container only sees the host paths that are **bound** into it. Nextflow auto-mounts the work and launch directories, but your **reads, references and Kraken2 database usually live elsewhere on a shared filesystem** - so the container can't read them and tasks fail with *"No such file or directory"* even though the paths are correct on the login node.

Binds are therefore per-site, like the queue names and the way Singularity is provided. BAMpiro keeps all of that in **site profiles** instead of in the defaults, so there are two cases.

=== "Garnatxa (I2SysBio)"

    Nothing to write. `-profile garnatxa` applies [`conf/garnatxa.config`](https://github.com/PathoGenOmics-Lab/BAMpiro/blob/indel-mask/conf/garnatxa.config), which already carries that cluster's bind paths, its `--qos` per step, `module load singularity`, the shared Kraken2 database, a reusable mappability cache and a work directory on scratch:

    ```bash
    nextflow run main.nf --tsv samples.tsv --outdir results -profile garnatxa
    ```

    It also works when the pipeline is pulled straight from GitHub, with no clone to configure:

    ```bash
    nextflow run PathoGenOmics-Lab/BAMpiro --tsv samples.tsv --outdir results -profile garnatxa
    ```

=== "Any other cluster"

    Start from `-profile slurm`, which submits with `sbatch` and assumes nothing else, and add your site's specifics in a small config passed with `-c`. Copy `conf/garnatxa.config` as the template and replace the paths and queue names with yours:

    ```groovy title="my_site.config"
    process {
        clusterOptions = '--qos=normal'            // your QoS / partition / account
        beforeScript   = 'module load singularity' // however your site provides Singularity

        // the steps that outgrow a short queue
        withName: 'MERGE_AND_MARKDUP|CALL_FREEBAYES|KRAKEN_FILTER_.*' {
            clusterOptions = '--qos=long'
        }
    }

    singularity {
        // every host directory your samplesheet paths + --kraken2_db sit under,
        // comma-separated, no spaces (parent dirs are fine: /data covers /data/*).
        // Add $HOME too if you let Nextflow pull the pipeline into ~/.nextflow/assets.
        runOptions = "--bind /scratch/me/reads,/data/refs,/shared/kraken2,${System.getenv('HOME') ?: '/home'}"
    }

    // keep work/ off a quota-limited home; -resume needs it to survive
    workDir = '/scratch/me/bampiro/work'
    ```

    ```bash
    nextflow run main.nf --tsv samples.tsv --outdir results -profile slurm -c my_site.config
    ```

    If your site is one other people use, that file is worth contributing back as a profile of its own.

!!! tip "Binding without a config file"

    For a one-off you can set the binds straight from the environment instead:

    ```bash
    export NXF_SINGULARITY_RUN_OPTIONS="--bind /scratch/me/reads,/data/refs,/shared/kraken2"
    ```

    Bind **parent** directories so every file underneath is visible.

## 2. Set parameters

Any parameter can be set on the command line with a `--` prefix, for example:

```bash
nextflow run main.nf --tsv samples.tsv --outdir results \
  -profile local,docker --threads 16 --freebayes_ploidy 1
```

To see what is available, ask the pipeline rather than the docs:

```bash
nextflow run main.nf --help    # every parameter, with its default, then exits
```

Typos are caught before any work starts. A `--flag` BAMpiro does not know stops the run with the nearest match instead of being silently ignored:

```console
$ nextflow run main.nf --tsv samples.tsv --freebayes_ploydi 1
Unrecognised parameter(s): --freebayes_ploydi
  --freebayes_ploydi: did you mean --freebayes_ploidy?
Run with --help to list every parameter.
```

For anything you want to reuse, collect your overrides in a config file and pass it with `-c` - command-line `--flags` still win over the file:

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

!!! tip "Not working with TB?"

    The dedicated tutorial [**Running a non-TB organism**](other-organisms.md) walks through
    the whole thing end to end - your reference, the `generic` profile, keeping the typing off,
    tuning thresholds, and a worked *E. coli* example.

## 4. Turn on the MTBC features

For a TB cohort like the bundled 17-sample demo, the built-in defaults already give you a valid run. To add the full recommended MTBC feature set, switch on lineage/DR typing and dual amino-acid numbering:

```bash
nextflow run main.nf --tsv samples.tsv --outdir results -profile local,docker \
  --run_pathotypr true --annotate_canonical true
```

- `--run_pathotypr true` enables alignment-free MTBC lineage + WHO drug-resistance typing ([Lineage & Drug-Resistance Typing](../pathotypr.md)).
- `--annotate_canonical true` adds dual **amino-acid** numbering - re-annotating variants against H37Rv, with Mycobrowser gene links.

!!! tip "Or use the ready-made TB config"

    [`conf/tuberculosis.config`](https://github.com/PathoGenOmics-Lab/BAMpiro/blob/indel-mask/conf/tuberculosis.config) bundles those two plus the H37Rv coordinate liftover (`--variant_liftover`) and blind-spot masking (`--mask_blindspots`). Apply it with `-c conf/tuberculosis.config`.

## 5. A word on the QC gate

Two things are easy to confuse. The interactive QC report lets you drag the cut-offs and re-flags samples **live** in your browser - that never changes any file. The pipeline's own pass/fail **gate** is separate: it uses the `--report_*` thresholds from your config, and only stops the run when you opt in with `--report_gate true`.

```bash
# Fail the whole run if any sample is flagged FAIL, using a stricter depth cut-off
--report_gate true --report_depth_min 20
```

The gate thresholds (`--report_depth_min`, `--report_breadth_min`, `--report_missing_max`, and so on) are listed in full under [Configuration](../configuration.md). You'll see exactly how each one flags a sample in the next tutorial.

## 6. Resuming & monitoring a run

Runs are long, and clusters kill jobs. Nextflow caches every finished task in the `work/` directory, so you rarely have to start over.

- **Resume** - if a run fails or you cancel it, re-run the **exact same command** with `-resume` added. Nextflow reuses the completed tasks and only recomputes what's actually needed:

    ```bash
    nextflow run main.nf --tsv samples.tsv --outdir results -profile local,docker -resume
    ```

- **Execution reports** - every run already writes a resource-usage report, a per-task trace and a timeline under `<outdir>/pipeline_info/`; they are what you need to right-size CPU / memory / time requests on a cluster. To send them somewhere else, name them yourself:

    ```bash
    nextflow run main.nf --tsv samples.tsv --outdir results -profile slurm \
      -with-report report.html -with-trace trace.txt -with-timeline timeline.html
    ```

    `-with-report` is an HTML resource-usage summary, `-with-trace` a per-task TSV, and `-with-timeline` a Gantt-style timeline.

!!! note "Keep `work/`… until you're done"

    `-resume` only works while the `work/` directory survives - clean it (`nextflow clean` or `rm -rf work`) **after** a run finishes and you've collected your `--outdir` results, not before.

## What's next

You've tuned a run for your environment and organism. Next, learn to read what it produced → [Reading the QC report](reading-the-qc-report.md).
