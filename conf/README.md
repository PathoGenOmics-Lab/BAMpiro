# `conf/` — ready-made example configs

Layer either on top of BAMpiro's defaults with `-c`:

- **[`tuberculosis.config`](tuberculosis.config)** — the recommended *M. tuberculosis* run (the **default**
  organism): turns on the full MTBC feature set — lineage + drug-resistance typing, dual H37Rv amino-acid
  numbering, the H37Rv coordinate liftover, and blind-spot masking.

  ```bash
  nextflow run main.nf --tsv samples.tsv --outdir results -profile standard -c conf/tuberculosis.config
  ```

- **[`organism.config`](organism.config)** — a commented, copy-and-edit template for **any other organism**:
  ploidy 1, the MTBC-only features off, and placeholders for a Kraken DB / canonical numbering / gate
  thresholds. Pair it with the built-in `generic` profile.

  ```bash
  nextflow run main.nf --tsv samples.tsv --outdir results -profile standard,generic -c my_organism.config
  ```

See [docs/configuration.md → Example configs](../docs/configuration.md#example-configs) for the details.
