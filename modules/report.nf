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

process COLLECT_DR {
    tag "DR calls"
    publishDir "${params.outdir}", mode: params.publish_mode
    cpus 1
    memory '2 GB'

    input:
    path(dr_mutations)      // every sample's ${sampleId}.dr_mutations.tsv (pathotypr DR run)
    val(basename)

    output:
    path("${basename}_dr.tsv"), emit: dr

    script:
    """
    set -euo pipefail
    python3 ${projectDir}/bin/collect_dr.py -o ${basename}_dr.tsv ${dr_mutations}
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
    path(gff)               // reference GFF3 -> per-gene SNP-density hotspots panel
    path(mask_bed)          // reference repeat/exclude BED -> masked-regions / callability panel
    path(metadata)          // samplesheet/metadata TSV -> SNP dynamics panel (auto-detects time+group)
    path(vcfs)              // per-sample annotated VCFs -> per-SNP allele frequencies for dynamics
    path(vcfs_h37rv)        // per-sample VCFs annotated vs H37Rv -> dual amino-acid numbering (may be NO_FILE)
    path(dr_report)         // run drug-resistance calls TSV (collect_dr) -> Drug resistance panel (may be NO_FILE)
    path(kraken_reports)    // per-sample Kraken2 .report files -> Taxonomic composition panel (may be empty)
    val(provenance)         // pre-quoted provenance tokens (container=..., reference=...)
    val(basename)

    output:
    path("${basename}_qc_report.html"), emit: html
    path("${basename}_qc_flags.tsv"),   emit: flags

    script:
    def gate_arg = params.report_gate ? "--gate" : ""
    def cons_arg = consensus ? "--consensus ${consensus}" : ""
    def palette  = "${projectDir}/assets/mycolorsTB_nature.tsv"
    """
    set -euo pipefail
    # Optional inputs self-hide their panel when absent/empty (parse_gff & parse_bed are tolerant).
    LC_ARG=""; [ -f "${palette}" ] && LC_ARG="--lineage-colors ${palette}"
    MASK_ARG=""; [ -s "${mask_bed}" ] && MASK_ARG="--mask-bed ${mask_bed}"
    MD_ARG=""; [ -s "${metadata}" ] && MD_ARG="--metadata ${metadata}"
    VCF_ARG=""; [ -n "${vcfs}" ] && VCF_ARG="--vcfs ${vcfs}"
    VH_ARG="";  case "${vcfs_h37rv}" in ""|NO_FILE) ;; *) VH_ARG="--vcfs-h37rv ${vcfs_h37rv}";; esac
    DR_ARG="";  [ -s "${dr_report}" ] && [ "${dr_report}" != "NO_FILE" ] && DR_ARG="--dr-report ${dr_report}"
    KRK_ARG=""; [ -n "${kraken_reports}" ] && KRK_ARG="--kraken ${kraken_reports}"
    python3 ${projectDir}/bin/qc_report.py \\
        --summary ${summary} \\
        ${cons_arg} \\
        --gene-burden ${gene_burden} \\
        --gff ${gff} \\
        \$MASK_ARG \$LC_ARG \$MD_ARG \$VCF_ARG \$VH_ARG \$DR_ARG \$KRK_ARG \\
        --aa2-label "${params.canonical_label}" \\
        --provenance ${provenance} \\
        --version "${workflow.manifest.version}" \\
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

process SNP_MATRIX {
    tag "SNP matrix"
    publishDir "${params.outdir}", mode: params.publish_mode
    cpus 1
    memory { 4.GB * task.attempt }

    input:
    path(vcfs)              // all per-sample annotated VCFs
    val(reference)          // reference id the samples were mapped against
    val(basename)

    output:
    path("${basename}_snp_matrix.tsv"), emit: matrix

    script:
    """
    set -euo pipefail
    python3 ${projectDir}/bin/build_snp_matrix.py \\
        --vcfs ${vcfs} \\
        --reference "${reference}" \\
        -o ${basename}_snp_matrix.tsv
    """
}
