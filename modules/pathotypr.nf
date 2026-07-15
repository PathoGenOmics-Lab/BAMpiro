nextflow.enable.dsl=2

// Import centralized functions for path generation and file classification
include { getSavePath; getSampleDir } from './utils'

/* ====================================================================
    PATHOTYPR MODULES
    Alignment-free (k-mer) MTBC lineage + WHO drug-resistance genotyping.
    Reference-agnostic: the diagnostic k-mers are built from the MTBC-ancestor
    reference the markers are defined on, so calls do NOT depend on the reference
    the sample was mapped against. Two split-fastq passes per sample:
      - lineage markers (nested sub-lineage)  -> emit: summary
      - WHO drug-resistance markers           -> emit: dr_summary

    NOTE: confirm the split-fastq flags and output filenames of bioconda pathotypr
    v1.0.0 with `pathotypr split-fastq --help`; the cp normalisation below adapts
    whatever it writes to the fixed names the stats/report steps read.
==================================================================== */

process RUN_PATHOTYPR_PE {
    tag "PathotyprPE: ${sampleId}"
    publishDir "${params.outdir}/${getSampleDir(sampleId, params)}", mode: params.publish_mode, saveAs: { filename -> getSavePath(filename, params) }
    cpus 4
    memory '8 GB'

    input:
    tuple val(refId), val(sampleId), val(runId), path(r1), path(r2)
    path ref_fasta_pathotypr
    path lineage_markers
    path dr_markers
    val pathotypr_bin

    output:
    tuple val(sampleId), path("${sampleId}.pathotypr.lineage.summary.tsv"), emit: summary
    tuple val(sampleId), path("${sampleId}.pathotypr.dr.summary.tsv"),      emit: dr_summary
    path("${sampleId}.pathotypr.*"), emit: results

    script:
    """
    set -euo pipefail
    # Lineage (nested sub-lineage classification)
    ${pathotypr_bin} split-fastq -i ${r1} -i ${r2} --paired \\
        --reference ${ref_fasta_pathotypr} --markers ${lineage_markers} \\
        --nested-classification --output-prefix ${sampleId}.pathotypr.lineage --threads ${task.cpus}
    # Drug resistance (WHO catalogue markers)
    ${pathotypr_bin} split-fastq -i ${r1} -i ${r2} --paired \\
        --reference ${ref_fasta_pathotypr} --markers ${dr_markers} \\
        --output-prefix ${sampleId}.pathotypr.dr --threads ${task.cpus}
    # Normalise to the fixed names the stats/report steps read (adapt to pathotypr's actual outputs)
    [ -f ${sampleId}.pathotypr.lineage.summary.tsv ] || cp "\$(ls ${sampleId}.pathotypr.lineage*.tsv | head -1)" ${sampleId}.pathotypr.lineage.summary.tsv
    [ -f ${sampleId}.pathotypr.dr.summary.tsv ]      || cp "\$(ls ${sampleId}.pathotypr.dr*.tsv | head -1)"      ${sampleId}.pathotypr.dr.summary.tsv
    """
}

process RUN_PATHOTYPR_SE {
    tag "PathotyprSE: ${sampleId}"
    publishDir "${params.outdir}/${getSampleDir(sampleId, params)}", mode: params.publish_mode, saveAs: { filename -> getSavePath(filename, params) }
    cpus 4
    memory '8 GB'

    input:
    tuple val(refId), val(sampleId), val(runId), path(r1)
    path ref_fasta_pathotypr
    path lineage_markers
    path dr_markers
    val pathotypr_bin

    output:
    tuple val(sampleId), path("${sampleId}.pathotypr.lineage.summary.tsv"), emit: summary
    tuple val(sampleId), path("${sampleId}.pathotypr.dr.summary.tsv"),      emit: dr_summary
    path("${sampleId}.pathotypr.*"), emit: results

    script:
    """
    set -euo pipefail
    ${pathotypr_bin} split-fastq -i ${r1} \\
        --reference ${ref_fasta_pathotypr} --markers ${lineage_markers} \\
        --nested-classification --output-prefix ${sampleId}.pathotypr.lineage --threads ${task.cpus}
    ${pathotypr_bin} split-fastq -i ${r1} \\
        --reference ${ref_fasta_pathotypr} --markers ${dr_markers} \\
        --output-prefix ${sampleId}.pathotypr.dr --threads ${task.cpus}
    [ -f ${sampleId}.pathotypr.lineage.summary.tsv ] || cp "\$(ls ${sampleId}.pathotypr.lineage*.tsv | head -1)" ${sampleId}.pathotypr.lineage.summary.tsv
    [ -f ${sampleId}.pathotypr.dr.summary.tsv ]      || cp "\$(ls ${sampleId}.pathotypr.dr*.tsv | head -1)"      ${sampleId}.pathotypr.dr.summary.tsv
    """
}
