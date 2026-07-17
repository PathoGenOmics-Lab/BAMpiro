# Troubleshooting

When a run fails, Nextflow names the failed process, its exit code, and the task's work
directory — `cd` there and read `.command.err` / `.command.log` / `.command.sh`. The
most common first-run problems:

| Symptom (error) | Cause | Fix |
| :--- | :--- | :--- |
| `Unknown configuration profile: '<x>'` | that profile isn't defined | the profiles are `standard` (SLURM), `local`, `docker`, `generic` — combine as needed (`-profile local,docker`) |
| `Cannot run program "sbatch"` / jobs never start | default targets SLURM | you're not on a SLURM cluster — add `-profile local` (see [Running without SLURM](installation.md#running-without-slurm)) |
| `module: command not found` in `.command.run` | base config runs `module load singularity` | use `-profile local` (or `docker`), which drops the module-load |
| `FATAL … mount source /scr/… no such file or directory` | the default binds site-specific Singularity paths | use `-profile local` (clears those binds); bind your own paths with `NXF_SINGULARITY_RUN_OPTIONS="--bind /your/path"` |
| `manifest for paururo/bampiro:1.0.1 not found` / auth denied | the image isn't pullable on your setup | pull it from the registry the maintainers published (Docker Hub or GHCR), or point `--container` at your own copy |
| Kraken2 task killed (OOM, exit 137) | Kraken2 loads the whole DB — it asks for **~80 GB** | run on a host / partition with enough RAM, or use a smaller DB; `--kraken_memory_mapping true` (default) mmaps the DB so parallel tasks share it |
| `Samplesheet not found: samples_legio.tsv` | `--tsv` was omitted | `--tsv` is required (no usable default) — pass your samplesheet |
| `Missing required column(s): refGff …` / SnpEff DB build fails | `refGff` is a **required** samplesheet column, and its seqid must match the FASTA | give every sample a GFF3 whose first column equals the FASTA header id — see [reference requirements](quickstart.md#reference-requirements) |
| Report has no **lineage** / **drug-resistance** panel | those need `--run_pathotypr true` | expected without pathotypr; panels self-hide when there is no lineage/DR data |
| Report has no **SNP dynamics** / **epistasis** panel | need per-sample VCFs + `--metadata` with a time+group column | add the metadata columns (see [Optional metadata](qc-report.md#optional-metadata)) |
| A sample failed with `SLURM … TIMEOUT` | the task exceeded the queue time limit | raise the time / use a longer `--qos` / partition for that step |

Still stuck? Re-run with `-resume` (skips finished work) and open an
[issue](https://github.com/PathoGenOmics-Lab/BAMpiro/issues) with the failed process
name and the contents of its `.command.err`.
