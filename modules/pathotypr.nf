nextflow.enable.dsl=2

// Import centralized functions for path generation and file classification
include { getSavePath; getSampleDir } from './utils'

/* ====================================================================
    PATHOTYPR MODULES
    Contains: Lineage classification processes (PE and SE)
==================================================================== */

process RUN_PATHOTYPR_PE {
    tag "PathotyprPE: ${sampleId}"
    
    // Use getSampleDir for nested output support.
    // getSavePath will automatically place these files into the 'lineage/' subdirectory.
    publishDir "${params.outdir}/${getSampleDir(sampleId, params)}", mode: 'copy', saveAs: { filename -> getSavePath(filename, params) }
    
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
    
    // Use getSampleDir for nested output support
    publishDir "${params.outdir}/${getSampleDir(sampleId, params)}", mode: 'copy', saveAs: { filename -> getSavePath(filename, params) }
    
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
