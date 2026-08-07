// Gene conversion detection.
//
// Both processes deliberately work OUTSIDE the pipeline's masking. Repeat exclusion and the
// length-aware read filter remove paralogous regions and multi-mapping reads, which is correct for
// variant calling and is exactly the signal a conversion lives in. See subworkflows/gene_conversion.nf.

include { getSampleDir; getSavePath } from './utils'

// The output columns of bin/gene_conversion.py, in one place. Three shell blocks need to write
// this header when there is nothing to write, and a copy that has drifted from the tool's own is
// worse than no header at all: a stub run then tests the wrong shape and says it passed.
// tests/unit/test_gene_conversion.py checks this list against the tool.
def gconvHeader() {
    return ['sample', 'pair_id', 'contig', 'donor', 'verdict', 'reason', 'start', 'end', 'span_bp',
            'post_conv', 'log10_bf', 'log10_bf_vs_null', 'tract_af', 'mismap_frac', 'mut_rate',
            'start_ci', 'end_ci',
            'n_sites', 'n_sites_outside', 'n_undetermined', 'donor_af_in', 'donor_af_outside',
            'min_depth', 'cis_reads', 'breakpoint_reads', 'donor_only_reads'].join('\\t')
}

def gconvLocusHeader() {
    return ['sample', 'pair_id', 'contig', 'donor', 'n_sites', 'n_reads', 'start', 'end',
            'n_tract_sites', 'log10_bf', 'log10_bf_vs_null', 'tract_af', 'mismap_frac',
            'mut_rate'].join('\\t')
}

process PARALOG_MAP {
    tag "Paralogs: ${refId}"
    publishDir path: { "${params.outdir}/references/${refId}" }, mode: params.publish_mode
    cpus 1
    memory { 4.GB * task.attempt }

    input:
    tuple val(refId), path(delta)

    output:
    tuple val(refId), path("${refId}.paralog_pairs.tsv"), emit: pairs
    tuple val(refId), path("${refId}.paralog_sites.tsv"), emit: sites

    script:
    """
    set -euo pipefail
    python3 ${projectDir}/bin/paralog_map.py \\
        --delta ${delta} \\
        --out-pairs ${refId}.paralog_pairs.tsv \\
        --out-sites ${refId}.paralog_sites.tsv \\
        --min-identity ${params.gconv_min_identity} \\
        --min-length ${params.gconv_min_paralog_length}
    """

    stub:
    """
    printf 'pair_id\\tacceptor\\tacc_start\\tacc_end\\tdonor\\tdon_start\\tdon_end\\tstrand\\tidentity\\tlength\\tn_diagnostic\\n' > ${refId}.paralog_pairs.tsv
    printf 'pair_id\\tacceptor\\tacc_pos\\tacc_base\\tdonor\\tdon_pos\\tdon_base\\tstrand\\tkind\\tlength\\n' > ${refId}.paralog_sites.tsv
    """
}

process FIND_GENE_CONVERSION {
    tag "Gene conversion: ${sampleId}"
    publishDir path: { "${params.outdir}/${getSampleDir(sampleId, params)}" }, mode: params.publish_mode, saveAs: { filename -> getSavePath(filename, params) }
    cpus 1
    memory { 4.GB * task.attempt }

    input:
    // The BAM is the deduplicated, PRE-FILTER one on purpose: the mappability filter drops the
    // multi-mapping reads that carry the evidence.
    tuple val(sampleId), val(refId), path(bam), path(bai), path(sites)

    output:
    tuple val(sampleId), val(refId), path("${sampleId}.${refId}.gene_conversion.tsv"), emit: tracts
    // One row per paralog pair analysed, reported or not. Only the cohort pass reads it: a locus
    // that says nothing in one sample means something only next to the same locus in the others.
    tuple val(sampleId), val(refId), path("${sampleId}.${refId}.gene_conversion_loci.tsv"), emit: loci

    script:
    """
    set -euo pipefail
    python3 ${projectDir}/bin/gene_conversion.py \\
        --sites ${sites} \\
        --bam ${bam} \\
        --sample ${sampleId} \\
        --output ${sampleId}.${refId}.gene_conversion.tsv \\
        --output-loci ${sampleId}.${refId}.gene_conversion_loci.tsv \\
        --min-bf ${params.gconv_min_bf} \\
        --report-bf ${params.gconv_report_bf} \\
        --prior ${params.gconv_prior} \\
        --mean-tract-bp ${params.gconv_mean_tract_bp} \\
        --max-tract-bp ${params.gconv_max_tract_bp} \\
        --max-tracts ${params.gconv_max_tracts} \\
        --mut-rate ${params.gconv_mut_rate} \\
        --indel-factor ${params.gconv_indel_factor} \\
        --min-tract-af ${params.gconv_min_tract_af} \\
        --min-mismap ${params.gconv_min_mismap} \\
        --min-sites ${params.gconv_min_sites} \\
        --min-depth ${params.gconv_min_depth} \\
        --min-bq ${params.gconv_min_bq}
    """

    stub:
    """
    printf '${gconvHeader()}\\n' > ${sampleId}.${refId}.gene_conversion.tsv
    printf '${gconvLocusHeader()}\\n' > ${sampleId}.${refId}.gene_conversion_loci.tsv
    """
}

process COLLECT_GENE_CONVERSION {
    // Not a concatenation. Two questions a sample cannot answer about itself are answerable here:
    // whether a tract is this isolate or the reference (recurrence across the cohort), and whether
    // a signal too weak to call alone is the same tract another sample calls outright.
    tag "Gene conversion cohort"
    publishDir path: { "${params.outdir}" }, mode: params.publish_mode
    cpus 1
    memory '2 GB'

    input:
    path(tracts)
    path(loci)
    val(basename)

    output:
    path("${basename}_gene_conversion.tsv"), emit: cohort

    script:
    """
    set -euo pipefail
    python3 ${projectDir}/bin/gconv_cohort.py \\
        --tracts ${tracts} \\
        --loci ${loci} \\
        --output ${basename}_gene_conversion.tsv \\
        --min-bf ${params.gconv_min_bf} \\
        --corroborated-bf ${params.gconv_corroborated_bf} \\
        --ubiquitous ${params.gconv_ubiquitous} \\
        --min-samples ${params.gconv_cohort_min_samples}
    """

    stub:
    """
    printf '${gconvHeader()}\\tevent_id\\tevent_samples\\tevent_frac\\tcohort_verdict\\tcohort_mismap\\tcohort_bf_median\\n' > ${basename}_gene_conversion.tsv
    """
}
