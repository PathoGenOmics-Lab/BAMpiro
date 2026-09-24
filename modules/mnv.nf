nextflow.enable.dsl=2

include { getSavePath; getSampleDir } from './utils'

/* ====================================================================
    CODON-LEVEL CHANGES (get_MNV)
    Two SNPs in one codon change one amino acid. Annotated a base at a time they name two changes that
    are not there: CAC>GAG is His>Glu, where each base alone says His>Asp and His>Gln, and GAG>TTG is
    Glu>Leu where the first base alone says a stop. get_MNV groups a sample's calls by codon and reads
    the codon whole on the reads that span it, so two SNPs on different molecules of a mixed population
    are not merged into a change no molecule carries.
==================================================================== */

process GET_MNV {
    tag "MNV: ${sampleId}"
    publishDir path: { "${params.outdir}/${getSampleDir(sampleId, params)}" }, mode: params.publish_mode, saveAs: { filename -> getSavePath(filename, params) }
    // Reads the BAM only around the called SNPs and indels: light.
    cpus 2
    memory { 2.GB * task.attempt }

    input:
    // The SNPs and indels that pass (CALL_FREEBAYES codon_calls), the alignment they were called on,
    // and the reference with its GFF
    tuple val(sampleId), val(refId), path(calls), path(calls_tbi), path(bam), path(bai), path(ref_fa), path(gff)

    output:
    tuple val(sampleId), val(refId), path("${sampleId}.${refId}.MNV.tsv"), emit: tsv
    // BGZF-compressed but not indexed: get_MNV indexes with tabix, which its image does not carry
    // (`tabix -p vcf` indexes the file as it is).
    tuple val(sampleId), val(refId), path("${sampleId}.${refId}.MNV.vcf.gz"), emit: vcf
    path("${sampleId}.${refId}.mnv.summary.json"),  emit: summary
    path("${sampleId}.${refId}.mnv.manifest.json"), emit: manifest

    script:
    """
    set -euo pipefail
    # get_MNV names its outputs after its input: give the calls the sample's name.
    ln -s ${calls} ${sampleId}.${refId}.vcf.gz
    ln -s ${calls_tbi} ${sampleId}.${refId}.vcf.gz.tbi
    # get_MNV reads the GFF as plain text, and the samplesheet may give it gzipped.
    GFF=${gff}
    if [[ "${gff}" == *.gz ]]; then
        gzip -cd "${gff}" > get_mnv.gff
        GFF=get_mnv.gff
    fi
    # The same base and mapping quality floors FreeBayes called with, so a codon's reads are the ones
    # its calls rest on. No read or frequency floor here: the calls already passed theirs, and the
    # run's table (MNV_TABLE) keeps a codon change only where reads carry it.
    get_mnv --vcf ${sampleId}.${refId}.vcf.gz --bam ${bam} --fasta ${ref_fa} --gff \$GFF \\
        --gff-features ${params.mnv_gff_features} --translation-table ${params.mnv_translation_table} \\
        --quality ${params.freebayes_min_base_qual} --min-mapq ${params.freebayes_min_map_qual} \\
        --threads ${task.cpus} --both --vcf-gz \\
        --summary-json ${sampleId}.${refId}.mnv.summary.json \\
        --run-manifest ${sampleId}.${refId}.mnv.manifest.json
    """

    stub:
    """
    touch ${sampleId}.${refId}.MNV.tsv
    touch ${sampleId}.${refId}.MNV.vcf.gz
    touch ${sampleId}.${refId}.mnv.summary.json
    touch ${sampleId}.${refId}.mnv.manifest.json
    """
}

process MNV_TABLE {
    tag "MNV table"
    publishDir "${params.outdir}", mode: params.publish_mode
    cpus 1
    memory { 2.GB * task.attempt }

    input:
    path(manifest)                          // sample<TAB>reference<TAB>get_MNV TSV, one line per sample
    path(tsvs, stageAs: 'mnv/*')            // every sample's get_MNV TSV
    val(basename)

    output:
    path("${basename}_mnv.tsv"), emit: table

    script:
    """
    set -euo pipefail
    awk -F'\\t' 'BEGIN { OFS = "\\t" } NF >= 3 { print \$1, \$2, "mnv/" \$3 }' ${manifest} > manifest.staged.tsv
    python3 ${projectDir}/bin/collect_mnv.py --manifest manifest.staged.tsv \\
        --min-reads ${params.mnv_min_reads} -o ${basename}_mnv.tsv
    """

    stub:
    """
    touch ${basename}_mnv.tsv
    """
}
