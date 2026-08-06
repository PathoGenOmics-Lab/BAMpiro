# `conf/` - ready-made configs

## Site configs (applied by a profile)

- **[`garnatxa.config`](garnatxa.config)** - the I2SysBio / University of Valencia cluster, applied by
  `-profile garnatxa`: SLURM QoS per step, `module load singularity`, the Singularity bind paths, the
  shared Kraken2 database, a reusable mappability cache and a work directory on scratch. Nothing to pass
  with `-c`, and it works when the pipeline is pulled straight from GitHub.

  ```bash
  nextflow run PathoGenOmics-Lab/BAMpiro --tsv samples.tsv --outdir results -profile garnatxa
  ```

  On any other cluster, start from the generic `-profile slurm` and copy this file as the template for
  your own site config. `test.config` / `test_full.config` back `-profile test` / `-profile test_full`.

## Organism configs (applied with `-c`)

Layer either on top of BAMpiro's defaults with `-c`:

- **[`tuberculosis.config`](tuberculosis.config)** - the recommended *M. tuberculosis* run (the **default**
  organism): turns on the full MTBC feature set, that is lineage + drug-resistance typing, dual H37Rv
  amino-acid numbering, the H37Rv coordinate liftover, and blind-spot masking.

  ```bash
  nextflow run main.nf --tsv samples.tsv --outdir results -profile local,docker -c conf/tuberculosis.config
  ```

- **[`organism.config`](organism.config)** - a commented, copy-and-edit template for **any other organism**:
  ploidy 1, the MTBC-only features off, and placeholders for a Kraken DB / canonical numbering / gate
  thresholds. Pair it with the built-in `generic` profile.

  ```bash
  nextflow run main.nf --tsv samples.tsv --outdir results -profile local,docker,generic -c my_organism.config
  ```

Swap `local,docker` for the executor profile you actually want (`slurm`, `garnatxa`, …); profiles compose
with a comma.

See [docs/configuration.md → Example configs](../docs/configuration.md#example-configs) for the details.
