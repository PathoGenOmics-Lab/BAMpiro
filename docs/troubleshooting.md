# Troubleshooting

When a run fails, Nextflow names the failed process, its exit code, and the task's work
directory - `cd` there and read `.command.err` / `.command.log` / `.command.sh`. The
most common first-run problems:

| Symptom (error) | Cause | Fix |
| :--- | :--- | :--- |
| `Unknown configuration profile: '<x>'` | that profile isn't defined | the profiles are `standard` (the default, this host), `local`, `slurm`, `garnatxa`, `docker`, `generic`, `test` - combine as needed (`-profile local,docker`) |
| The run executes on the login node, nothing reaches the queue | no executor profile was chosen; the default runs on the current host | add `-profile slurm` (or `-profile garnatxa`) - see [Choosing where it runs](installation.md#choosing-where-it-runs). The startup banner's `Executor` line tells you which one you got |
| `Unrecognised parameter(s): --…` | a misspelled or invented `--flag` | fix the spelling (the error suggests the nearest match); `nextflow run main.nf --help` lists every parameter |
| `module: command not found` in `.command.run` | a site profile ran `module load singularity`, but this machine has no module system | that comes from `-profile garnatxa`; on another cluster use `-profile slurm` plus your own `-c` config |
| `FATAL … mount source /scr/… no such file or directory` | `-profile garnatxa` binds paths that exist only on that cluster | use `-profile slurm` and bind your own paths (`-c` config, or `NXF_SINGULARITY_RUN_OPTIONS="--bind /your/path"`) |
| `manifest for paururo/bampiro@sha256:… not found` / auth denied | the pinned image isn't pullable on your setup | pull it from the registry the maintainers published (Docker Hub or GHCR), or point `--container` at your own copy |
| Kraken2 task killed (OOM, exit 137) | Kraken2 loads the whole DB - it asks for **~80 GB** | run on a host / partition with enough RAM, use a smaller DB, or drop `--kraken2_db` and skip the screen; `--kraken_memory_mapping true` (default) mmaps the DB so parallel tasks share it |
| `--tsv is required: a Tab-Separated samplesheet …` | `--tsv` was omitted | `--tsv` has no default - pass your samplesheet (the message lists the columns) |
| `Missing required column(s): refGff …` / SnpEff DB build fails | `refGff` is a **required** samplesheet column, and its seqid must match the FASTA | give every sample a GFF3 whose first column equals the FASTA header id - see [reference requirements](quickstart.md#reference-requirements) |
| Report has no **lineage** / **drug-resistance** panel | those need `--run_pathotypr true` | expected without pathotypr; panels self-hide when there is no lineage/DR data |
| Report has no **taxonomic composition** panel | `--kraken2_db` has no default, so the screen was skipped | point `--kraken2_db` at a database (and give the samplesheet a `taxId` column) |
| Report has no **SNP dynamics** / **epistasis** panel | need per-sample VCFs + `--metadata` with a time+group column | add the metadata columns (see [Optional metadata](qc-report.md#optional-metadata)) |
| A sample failed with `SLURM … TIMEOUT` | the task exceeded the queue time limit | give that step a longer QoS / partition through `process.clusterOptions` in your own `-c` config; `conf/garnatxa.config` shows the pattern (a `withName:` block for the slow steps) |

!!! question "Still stuck?"

    Re-run with `-resume` (skips finished work) and open an
    [issue](https://github.com/PathoGenOmics-Lab/BAMpiro/issues) with the failed process
    name and the contents of its `.command.err`.
