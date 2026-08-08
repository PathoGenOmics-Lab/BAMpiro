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
            'don_start', 'don_end',
            'post_conv', 'log10_bf', 'log10_bf_vs_null', 'tract_af', 'mismap_frac', 'mut_rate',
            'start_ci', 'end_ci',
            'n_sites', 'n_sites_outside', 'n_undetermined', 'donor_af_in', 'donor_af_outside',
            'min_depth', 'cis_reads', 'breakpoint_reads', 'donor_only_reads',
            'n_derived', 'n_ancestral', 'n_unpolarised', 'donor_swap_af',
            'locus_cn', 'expected_af',
            'genes', 'n_syn', 'n_nonsyn', 'aa_changes'].join('\\t')
}

def gconvLocusHeader() {
    return ['sample', 'pair_id', 'contig', 'donor', 'n_sites', 'n_reads', 'start', 'end',
            'n_tract_sites', 'log10_bf', 'log10_bf_vs_null', 'tract_af', 'mismap_frac',
            'mut_rate'].join('\\t')
}

process ALIGN_OUTGROUP {
    // An outgroup aligned against the reference, which is what says WHICH copy a difference is
    // on. Without it a sample carrying the donor's base and a reference carrying a derived one
    // are the same picture, and the reference's history gets reported as the sample's.
    tag "Outgroup: ${refId}"
    cpus 2
    memory { 8.GB * task.attempt }

    input:
    tuple val(refId), path(ref_fa), path(outgroup)

    output:
    tuple val(refId), path("${refId}.outgroup.delta"), emit: delta

    script:
    """
    set -euo pipefail
    nucmer --maxmatch --prefix ${refId}.outgroup ${ref_fa} ${outgroup} --threads ${task.cpus}
    """

    stub:
    """
    touch ${refId}.outgroup.delta
    """
}

process PARALOG_MAP {
    tag "Paralogs: ${refId}"
    publishDir path: { "${params.outdir}/references/${refId}" }, mode: params.publish_mode
    cpus 1
    memory { 4.GB * task.attempt }

    input:
    // The outgroup delta is optional and arrives as assets/NO_FILE_OUTGROUP when absent. With it
    // each diagnostic site says which copy the difference is on, so a sample carrying the donor's
    // base can be told from a reference that carries a derived one.
    tuple val(refId), path(delta), path(outgroup_delta)

    output:
    tuple val(refId), path("${refId}.paralog_pairs.tsv"), emit: pairs
    tuple val(refId), path("${refId}.paralog_sites.tsv"), emit: sites

    script:
    """
    set -euo pipefail
    ANC_ARG=""; case "${outgroup_delta}" in ""|NO_FILE*) ;; *) ANC_ARG="--ancestor-delta ${outgroup_delta}";; esac
    python3 ${projectDir}/bin/paralog_map.py \\
        --delta ${delta} \\
        --out-pairs ${refId}.paralog_pairs.tsv \\
        --out-sites ${refId}.paralog_sites.tsv \\
        --min-identity ${params.gconv_min_identity} \\
        --min-length ${params.gconv_min_paralog_length} \\
        \$ANC_ARG
    """

    stub:
    """
    printf 'pair_id\\tacceptor\\tacc_start\\tacc_end\\tdonor\\tdon_start\\tdon_end\\tstrand\\tidentity\\tlength\\tn_diagnostic\\n' > ${refId}.paralog_pairs.tsv
    printf 'pair_id\\tacceptor\\tacc_pos\\tacc_base\\tdonor\\tdon_pos\\tdon_base\\tstrand\\tkind\\tlength\\tanc_acc\\tanc_don\\tpolarity\\tdon_derived\\n' > ${refId}.paralog_sites.tsv
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
    tuple val(sampleId), val(refId), path(bam), path(bai), path(sites), path(ref_fa), path(gff)

    output:
    tuple val(sampleId), val(refId), path("${sampleId}.${refId}.gene_conversion.tsv"), emit: tracts
    // One row per paralog pair analysed, reported or not. Only the cohort pass reads it: a locus
    // that says nothing in one sample means something only next to the same locus in the others.
    tuple val(sampleId), val(refId), path("${sampleId}.${refId}.gene_conversion_loci.tsv"), emit: loci

    script:
    """
    set -euo pipefail
    # Without a GFF the tract is coordinates and nothing else, which is a finding nobody can act
    # on. With one it names the genes it landed on and what its copied bases do to their proteins.
    GFF_ARG=""; case "${gff}" in ""|NO_FILE*) ;; *) GFF_ARG="--gff ${gff} --reference ${ref_fa}";; esac
    # The sample's own genome-wide depth, so a locus can be compared with it. One extra pass
    # over an already indexed BAM, and it is what turns "a diluted tract" into "a clonal
    # conversion of one copy out of three".
    GDEPTH=\$(samtools coverage ${bam} \\
        | awk 'NR>1 && \$3>0 {len+=\$3; sum+=\$3*\$7} END {print (len>0) ? sum/len : 0}')
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
        --min-bq ${params.gconv_min_bq} \\
        --reciprocal-af ${params.gconv_reciprocal_af} \\
        --genome-depth "\$GDEPTH" \\
        \$GFF_ARG
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
        --min-samples ${params.gconv_cohort_min_samples} \\
        --donor-margin ${params.gconv_donor_margin}
    """

    stub:
    """
    printf '${gconvHeader()}\\tevent_id\\tevent_samples\\tevent_frac\\tcohort_verdict\\tcohort_mismap\\tcohort_bf_median\\tdonor_rank\\tn_donors\\tdonor_margin\\tdonor_call\\tis_representative\\n' > ${basename}_gene_conversion.tsv
    """
}
