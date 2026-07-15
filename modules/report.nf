nextflow.enable.dsl=2

/* ====================================================================
    REPORT MODULES
    Cohort-level consolidated QC report (ported from sBAMpiro):
      COLLECT_SUMMARY -> one summary TSV from all per-sample legacy logs
      QC_REPORT       -> self-contained interactive HTML + qc_flags.tsv
==================================================================== */

process COLLECT_SUMMARY {
    tag "Summary"
    publishDir "${params.outdir}", mode: params.publish_mode
    cpus 1
    memory '2 GB'

    input:
    path(logs)              // every ${sampleId}.log (collected)
    val(basename)

    output:
    path("${basename}_summary.tsv"),     emit: summary
    path("${basename}_gene_burden.tsv"), emit: gene_burden

    script:
    """
    set -euo pipefail
    python3 ${projectDir}/bin/collect_summary.py \\
        -o ${basename}_summary.tsv \\
        --gene-burden-out ${basename}_gene_burden.tsv \\
        ${logs}
    # collect_summary only writes the gene-burden file when snpEff ANN is present in the logs;
    # guarantee it exists (header only) so the channel/QC_REPORT always have a file.
    [ -f ${basename}_gene_burden.tsv ] || printf 'gene\\thigh\\tmoderate\\tlow\\tmodifier\\tdominant_effect\\tn_samples\\ttotal_impactful\\tstart\\tend\\n' > ${basename}_gene_burden.tsv
    """
}

process QC_REPORT {
    tag "QC report"
    publishDir "${params.outdir}", mode: params.publish_mode
    cpus 1
    memory '4 GB'

    input:
    path(summary)
    path(gene_burden)
    path(consensus)         // all masked consensus FASTAs (may be empty)
    val(basename)

    output:
    path("${basename}_qc_report.html"), emit: html
    path("${basename}_qc_flags.tsv"),   emit: flags

    script:
    def gate_arg = params.report_gate ? "--gate" : ""
    def cons_arg = consensus ? "--consensus ${consensus}" : ""
    """
    set -euo pipefail
    python3 ${projectDir}/bin/qc_report.py \\
        --summary ${summary} \\
        ${cons_arg} \\
        --gene-burden ${gene_burden} \\
        --out-html ${basename}_qc_report.html \\
        --out-flags ${basename}_qc_flags.tsv \\
        --title "BAMpiro QC report" \\
        --depth-min ${params.report_depth_min} \\
        --breadth-min ${params.report_breadth_min} \\
        --missing-max ${params.report_missing_max} \\
        --dup-max ${params.report_dup_max} \\
        --iupac-max ${params.report_iupac_max} \\
        --mapping-min ${params.report_mapping_min} \\
        ${gate_arg}
    """
}
