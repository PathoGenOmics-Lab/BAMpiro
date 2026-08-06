# `NO_FILE*` placeholders

Nextflow has no notion of an optional process input: a `path` input must always
receive something. When one of BAMpiro's optional features is off, `main.nf`
stages one of these empty files instead, and the process recognises it by the
`NO_FILE*` prefix and skips the corresponding argument.

| File | Stands in for | Off when |
| :--- | :--- | :--- |
| `NO_FILE` | pathotypr lineage summary (`GENERATE_LEGACY_STATS`) | pathotypr did not type that sample |
| `NO_FILE_DR` | drug-resistance calls (`QC_REPORT`) | `--run_pathotypr false` |
| `NO_FILE_H37RV` | canonically annotated VCFs (`QC_REPORT`) | `--annotate_canonical false` |
| `NO_FILE_LIFTOVER` | variant coordinate map (`QC_REPORT`) | `--variant_liftover false` |

Two properties matter, and both are load-bearing:

- **The names are distinct.** `QC_REPORT` takes three of these at once and all
  three features are off by default. Nextflow refuses to stage two inputs of the
  same task under the same filename, so sharing one `NO_FILE` fails the run with
  `input file name collision`.
- **The files exist and are empty.** A path that does not exist stages as a
  dangling symlink, which happens to work with the default staging mode but fails
  outright under `stageInMode 'copy'`. Empty keeps the `[ -s ... ]` guards false.

They are referenced as `${projectDir}/assets/...`, so a file called `NO_FILE`
sitting in someone's launch directory is never picked up by mistake.
