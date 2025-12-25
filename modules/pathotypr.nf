nextflow.enable.dsl=2

/* ====================================================================
    PATHOTYPR MODULES
    Contains: Lineage classification processes (PE and SE)
==================================================================== */

process RUN_PATHOTYPR_PE {
    tag "PathotyprPE: ${sampleId}"
    // Publish results to the 'lineage' subdirectory
    publishDir "${params.outdir}/${sampleId}/lineage", mode: 'copy'
    cpus 4
    memory '8 GB'

    input:
    tuple val(refId), val(sampleId), val(runId), path(r1), path(r2) 
    path ref_fasta_pathotypr 
    path markers_file
    val pathotypr_bin

    output:
    // Specific capture of the summary file for downstream statistics
    tuple val(sampleId), path("*_summary.tsv"), emit: summary
    // Capture all other output files (lineage report, kmer counts, etc.)
    path("${sampleId}*"), emit: results

    script:
    """
    set -euo pipefail
    
    # Run Pathotypr in Paired-End mode
    ${pathotypr_bin} split-fastq \\
        -i ${r1} -i ${r2} --paired \\
        --reference ${ref_fasta_pathotypr} \\
        --markers ${markers_file} \\
        --nested-classification \\
        --output-prefix ${sampleId}.pathotypr.pe \\
        --threads ${task.cpus}
    """
}

process RUN_PATHOTYPR_SE {
    tag "PathotyprSE: ${sampleId}"
    // Publish results to the 'lineage' subdirectory
    publishDir "${params.outdir}/${sampleId}/lineage", mode: 'copy'
    cpus 4
    memory '8 GB'

    input:
    tuple val(refId), val(sampleId), val(runId), path(r1)
    path ref_fasta_pathotypr 
    path markers_file
    val pathotypr_bin

    output:
    // Specific capture of the summary file for downstream statistics
    tuple val(sampleId), path("*_summary.tsv"), emit: summary
    // Capture all other output files
    path("${sampleId}*"), emit: results

    script:
    """
    set -euo pipefail
    
    # Run Pathotypr in Single-End mode
    ${pathotypr_bin} split-fastq \\
        -i ${r1} \\
        --reference ${ref_fasta_pathotypr} \\
        --markers ${markers_file} \\
        --nested-classification \\
        --output-prefix ${sampleId}.pathotypr.se \\
        --threads ${task.cpus}
    """
}
