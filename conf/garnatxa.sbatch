#!/bin/bash
# Launch BAMpiro on Garnatxa (I2SysBio / Universitat de Valencia).
#
#   sbatch conf/garnatxa.sbatch samples.tsv results
#
# What this submits is the DRIVER, not the work. Nextflow itself uses one core and a couple of
# gigabytes: it submits every task as its own job and waits. Garnatxa's own documentation does it
# this way rather than leaving the driver on the login node, where a long run is both against the
# usage rules and liable to be killed under you.
#
# The per-task resources are NOT here. They come from `-profile garnatxa`, which sets the QoS per
# step, the Singularity bind paths, the shared Kraken2 database, the mappability cache and a work
# directory on scratch. Change those in conf/garnatxa.config, not in this file.

#SBATCH --job-name=bampiro
#SBATCH --qos=medium
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
# Explicit, and not optional: every QoS on Garnatxa defaults to SIX HOURS. The driver has to
# outlive the whole cohort, and when it is killed every task it was waiting on is orphaned.
# medium allows up to 7 days; raise this if the cohort is large, and use --qos=long past 7 days.
#SBATCH --time=3-00:00:00
#SBATCH --output=bampiro_driver_%j.log
#SBATCH --error=bampiro_driver_%j.log

set -euo pipefail

# ---------------------------------------------------------------- what to run

SAMPLESHEET="${1:-samples.tsv}"
OUTDIR="${2:-results_bampiro}"

# The pipeline revision. A branch, a tag or a commit: whatever you pass is what runs, and it is
# recorded in the report's provenance footer.
REVISION="feat/gene-conversion"

# Off by default. Needs --exclude_repeats true, which is also the default.
FIND_GENE_CONVERSION="false"

# The Kraken2 database is ALREADY SET by -profile garnatxa, to the shared copy at
# /storage/shared_datasets/kraken2. Leave this empty to use it. Point it somewhere else only if
# you want a different database; an empty string here means "whatever the profile says".
#
# Two things to know either way:
#   - the samplesheet needs a `taxId` column, or there is nothing to filter against
#   - Kraken2 loads the whole database and asks for ~80 GB, which is the memory ceiling of the run
KRAKEN2_DB=""

# ---------------------------------------------------------------- preflight

if [ ! -s "${SAMPLESHEET}" ]; then
    echo "samplesheet not found or empty: ${SAMPLESHEET}" >&2
    echo "usage: sbatch $0 <samplesheet.tsv> [outdir]" >&2
    exit 2
fi

module load nextflow
if ! command -v nextflow >/dev/null; then
    echo "nextflow did not load. Check \`module avail nextflow\` on this cluster" >&2
    exit 2
fi

echo "BAMpiro on $(hostname), SLURM job ${SLURM_JOB_ID:-none}"
echo "  samplesheet : ${SAMPLESHEET}"
echo "  outdir      : ${OUTDIR}"
echo "  revision    : ${REVISION}"
echo "  nextflow    : $(nextflow -version 2>&1 | awk '/version/{print $3; exit}')"

# ---------------------------------------------------------------- run

# Built as a string rather than an array: expanding an EMPTY array under `set -u` is an
# unbound-variable error on bash 3.x, and the usual workaround reads worse than this does.
EXTRA=""
if [ -n "${KRAKEN2_DB}" ]; then
    EXTRA="${EXTRA} --kraken2_db ${KRAKEN2_DB}"
fi
if [ "${FIND_GENE_CONVERSION}" = "true" ]; then
    EXTRA="${EXTRA} --find_gene_conversion true"
fi

# -resume reuses finished work, so a driver that hit its time limit costs only the tasks that
# were still running. -ansi-log false keeps the log readable in a file rather than redrawing it.
nextflow run PathoGenOmics-Lab/BAMpiro \
    -r "${REVISION}" \
    -profile garnatxa \
    -resume \
    -ansi-log false \
    --tsv "${SAMPLESHEET}" \
    --outdir "${OUTDIR}" \
    ${EXTRA}

# The pipeline writes its own execution report, timeline and trace under
# <outdir>/pipeline_info/, so there is nothing to add with -with-report here.
echo "done. QC report: ${OUTDIR}/$(basename "${SAMPLESHEET}" .tsv)_qc_report.html"
