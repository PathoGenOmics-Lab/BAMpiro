/*
===============================================================================
 Stage 6b: Codon-level changes
===============================================================================
 get_MNV reads the SNPs that share a codon whole, on the reads that carry them:
 the amino-acid change a per-base annotation misses.
*/

include { GET_MNV } from '../modules/mnv'
include { asBool } from '../modules/utils'

workflow CODON_CHANGES {

    take:
    codon_calls   // (sId, rId, codon_calls.vcf.gz, tbi): the passing SNPs + indels
    vbase         // (sId, rId, bam, bai, ref_fa, exclude_txt): the alignment they were called on
    refGffMap     // refId -> refGff

    main:
    def mnv_tsv = Channel.empty()
    if (asBool(params.run_mnv)) {
        def mnv_in = codon_calls.join(vbase, by: [0, 1])
            .map { sId, rId, vcf, tbi, bam, bai, ref_fa, excl -> tuple(sId, rId, vcf, tbi, bam, bai, ref_fa, file(refGffMap[rId])) }
        mnv_tsv = GET_MNV(mnv_in).tsv
    }

    emit:
    tsv = mnv_tsv // (sId, rId, <sample>.<ref>.MNV.tsv); empty when run_mnv is off
}
