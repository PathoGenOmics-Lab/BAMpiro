/*
===============================================================================
 Stage 1: Reference Preparation
===============================================================================
 Index every reference in the samplesheet and build its snpEff database.
*/

include { PREPARE_REFERENCE; SNPEFF_BUILD_DB } from '../modules/reference'

workflow PREPARE_REFERENCES {

    take:
    ref_in_ch       // (refId, refFasta)
    snpeff_db_in    // (refId, refFasta, refGff)

    main:

    // 1. Reference Preparation
    def ref_out    = PREPARE_REFERENCE(ref_in_ch)
    def snpeff_out = SNPEFF_BUILD_DB(snpeff_db_in)

    emit:
    bundle = ref_out.bundle       // (refId, ref_fa, indices, exclude_txt)
    snpeff = snpeff_out.db        // (refId, snpEff.config, data)
    // The nucmer self-alignment, already computed for repeat masking. Only the gene-conversion
    // stage reads it, and it is absent when exclude_repeats is off.
    delta  = ref_out.delta        // (refId, self_aln.delta)
}
