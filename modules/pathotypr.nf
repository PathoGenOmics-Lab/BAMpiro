nextflow.enable.dsl=2

// Import centralized functions for path generation and file classification
include { getSavePath; getSampleDir } from './utils'

/* ====================================================================
    PATHOTYPR MODULES
    Alignment-free (k-mer) MTBC lineage + WHO drug-resistance genotyping.
    Reference-agnostic: the diagnostic k-mers are built from the MTBC-ancestor
    reference the markers are defined on, so calls do NOT depend on the reference
    the sample was mapped against. Two split-fastq passes per sample:
      - lineage markers (nested sub-lineage) -> emit: summary   (${prefix}_summary.tsv:
        genome / lineage:count / major_lineage; parsed by stats_to_legacy.py)
      - WHO drug-resistance markers          -> emit: dr_mutations (${prefix}_<s>_mutations.tsv:
        pos/ref/alt/ref_count/alt_count/alt_fraction/lineage_path[=drug;resistance;marker;grade;gene;mutation];
        aggregated by collect_dr.py for the report's Drug-resistance panel)
    pathotypr v1.0.0 CLI verified: split-fastq -i -r -m -o --paired --nested-classification --threads.
==================================================================== */

process RUN_PATHOTYPR_PE {
    tag "PathotyprPE: ${sampleId}"
    publishDir path: { "${params.outdir}/${getSampleDir(sampleId, params)}" }, mode: params.publish_mode, saveAs: { filename -> getSavePath(filename, params) }
    cpus 4
    memory '8 GB'

    input:
    tuple val(refId), val(sampleId), val(runId), path(r1), path(r2)
    path ref_fasta_pathotypr
    path lineage_markers
    path dr_markers
    val pathotypr_bin

    output:
    // runId in the file names so a multi-lane (or multi-reference) sample does not emit identically-named
    // files that collide when COLLECT_DR .collect()s them; the tuple key stays sampleId for the joins.
    tuple val(sampleId), path("${sampleId}__${runId}.pathotypr.lineage_summary.tsv"), emit: summary
    tuple val(sampleId), path("${sampleId}__${runId}.dr_mutations.tsv"),              emit: dr_mutations
    path("${sampleId}__${runId}.pathotypr.*"), emit: results

    script:
    """
    set -euo pipefail
    # Lineage (nested sub-lineage classification) -> \${prefix}_summary.tsv
    ${pathotypr_bin} split-fastq -i ${r1} -i ${r2} --paired \\
        --reference ${ref_fasta_pathotypr} --markers ${lineage_markers} \\
        --nested-classification --output-prefix ${sampleId}__${runId}.pathotypr.lineage --threads ${task.cpus}
    # Drug resistance (WHO catalogue markers) -> \${prefix}_<sample>_mutations.tsv
    ${pathotypr_bin} split-fastq -i ${r1} -i ${r2} --paired \\
        --reference ${ref_fasta_pathotypr} --markers ${dr_markers} \\
        --output-prefix ${sampleId}__${runId}.pathotypr.dr --min-alt-percent ${params.pathotypr_min_alt} --threads ${task.cpus}
    # Fix the DR detail file name to carry our sampleId+runId (pathotypr names the sample after the FASTQ)
    cp "\$(ls ${sampleId}__${runId}.pathotypr.dr_*_mutations.tsv | head -1)" ${sampleId}__${runId}.dr_mutations.tsv
    """

    // The lineage summary also satisfies the pathotypr.* glob (emit: results), so two files cover all three outputs.
    stub:
    """
    touch ${sampleId}__${runId}.pathotypr.lineage_summary.tsv
    touch ${sampleId}__${runId}.dr_mutations.tsv
    """
}

process RUN_PATHOTYPR_SE {
    tag "PathotyprSE: ${sampleId}"
    publishDir path: { "${params.outdir}/${getSampleDir(sampleId, params)}" }, mode: params.publish_mode, saveAs: { filename -> getSavePath(filename, params) }
    cpus 4
    memory '8 GB'

    input:
    tuple val(refId), val(sampleId), val(runId), path(r1)
    path ref_fasta_pathotypr
    path lineage_markers
    path dr_markers
    val pathotypr_bin

    output:
    // runId in the file names (see RUN_PATHOTYPR_PE) so multi-lane / multi-ref samples don't collide in COLLECT_DR
    tuple val(sampleId), path("${sampleId}__${runId}.pathotypr.lineage_summary.tsv"), emit: summary
    tuple val(sampleId), path("${sampleId}__${runId}.dr_mutations.tsv"),              emit: dr_mutations
    path("${sampleId}__${runId}.pathotypr.*"), emit: results

    script:
    """
    set -euo pipefail
    ${pathotypr_bin} split-fastq -i ${r1} \\
        --reference ${ref_fasta_pathotypr} --markers ${lineage_markers} \\
        --nested-classification --output-prefix ${sampleId}__${runId}.pathotypr.lineage --threads ${task.cpus}
    ${pathotypr_bin} split-fastq -i ${r1} \\
        --reference ${ref_fasta_pathotypr} --markers ${dr_markers} \\
        --output-prefix ${sampleId}__${runId}.pathotypr.dr --min-alt-percent ${params.pathotypr_min_alt} --threads ${task.cpus}
    cp "\$(ls ${sampleId}__${runId}.pathotypr.dr_*_mutations.tsv | head -1)" ${sampleId}__${runId}.dr_mutations.tsv
    """

    stub:
    """
    touch ${sampleId}__${runId}.pathotypr.lineage_summary.tsv
    touch ${sampleId}__${runId}.dr_mutations.tsv
    """
}
