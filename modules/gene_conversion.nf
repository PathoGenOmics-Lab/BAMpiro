// Gene conversion detection.
//
// Both processes deliberately work OUTSIDE the pipeline's masking. Repeat exclusion and the
// length-aware read filter remove paralogous regions and multi-mapping reads, which is correct for
// variant calling and is exactly the signal a conversion lives in. See subworkflows/gene_conversion.nf.

include { getSampleDir; getSavePath } from './utils'

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
    printf 'pair_id\\tacceptor\\tacc_pos\\tacc_base\\tdonor\\tdon_pos\\tdon_base\\tstrand\\n' > ${refId}.paralog_sites.tsv
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

    script:
    """
    set -euo pipefail
    python3 ${projectDir}/bin/gene_conversion.py \\
        --sites ${sites} \\
        --bam ${bam} \\
        --sample ${sampleId} \\
        --output ${sampleId}.${refId}.gene_conversion.tsv \\
        --min-af ${params.gconv_min_af} \\
        --min-sites ${params.gconv_min_sites} \\
        --min-depth ${params.gconv_min_depth} \\
        --min-bq ${params.gconv_min_bq}
    """

    stub:
    """
    printf 'sample\\tpair_id\\tcontig\\tdonor\\tverdict\\treason\\tstart\\tend\\tspan_bp\\tn_sites\\tn_sites_outside\\tdonor_af_in\\tdonor_af_outside\\tmin_depth\\tcis_reads\\tbreakpoint_reads\\tdonor_only_reads\\n' > ${sampleId}.${refId}.gene_conversion.tsv
    """
}

process COLLECT_GENE_CONVERSION {
    tag "Gene conversion cohort"
    publishDir path: { "${params.outdir}" }, mode: params.publish_mode
    cpus 1
    memory '2 GB'

    input:
    path(tracts)
    val(basename)

    output:
    path("${basename}_gene_conversion.tsv"), emit: cohort

    script:
    """
    set -euo pipefail
    out=${basename}_gene_conversion.tsv
    first=1
    for f in ${tracts}; do
        if [ "\$first" = "1" ]; then head -1 "\$f" > "\$out"; first=0; fi
        tail -n +2 "\$f" >> "\$out"
    done
    [ -s "\$out" ] || printf 'sample\\tpair_id\\tcontig\\tdonor\\tverdict\\treason\\tstart\\tend\\tspan_bp\\tn_sites\\tn_sites_outside\\tdonor_af_in\\tdonor_af_outside\\tmin_depth\\tcis_reads\\tbreakpoint_reads\\tdonor_only_reads\\n' > "\$out"
    """

    stub:
    """
    printf 'sample\\tpair_id\\tcontig\\tdonor\\tverdict\\treason\\tstart\\tend\\tspan_bp\\tn_sites\\tn_sites_outside\\tdonor_af_in\\tdonor_af_outside\\tmin_depth\\tcis_reads\\tbreakpoint_reads\\tdonor_only_reads\\n' > ${basename}_gene_conversion.tsv
    """
}
