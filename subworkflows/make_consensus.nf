/*
===============================================================================
 Stage 7: Consensus Generation (Optional)
===============================================================================
 The masked consensus, plus the parallel virgin (unmasked) path built from the
 original dedup BAM.
*/

include { CONSENSUS_FASTA } from '../modules/consensus'
include { CALL_FREEBAYES_RAW } from '../modules/variants'
// Aliases for the parallel virgin (unmasked) consensus path (a DSL2 process runs once per name)
include { CALL_BACKBONE as CALL_BACKBONE_RAW; MERGE_VCFS as MERGE_VCFS_RAW } from '../modules/variants'
include { CONSENSUS_FASTA as CONSENSUS_FASTA_RAW } from '../modules/consensus'
include { asBool } from '../modules/utils'

workflow MAKE_CONSENSUS {

    take:
    vbase          // (sId, rId, bam, bai, ref_fa, exclude_txt): the BAM variants were called on
    dedup_bam      // the ORIGINAL dedup BAM, before any read filtering
    allpos_vcf     // (sId, rId, all.pos.vcf.gz, tbi)
    fb_snps        // (sId, rId, valid_snps.vcf.gz, tbi)
    fb_mask_sites  // (sId, rId, mask_sites.tsv)

    main:

    // 7. Consensus Generation (Optional)
    def masked_consensus = Channel.empty()   // captured for the cohort QC report (section 11)
    if (asBool(params.make_consensus)) {
        // Prepare inputs: VCF + Reference + Mask sites
        // Note: Script is called from bin/ directly in the module
        // refmeta from vbase so the masked consensus gets the combined (nucmer + genmap) exclude
        def refmeta = vbase.map { sId, rId, bam, bai, ref_fa, exclude_txt -> tuple(sId, rId, ref_fa, exclude_txt) }
        def allpos_mask = allpos_vcf.join(fb_mask_sites, by: [0,1])
        def cons_in = allpos_mask.join(refmeta, by: [0,1]).map { sId, rId, vcf_gz, tbi, mask, ref_fa, exclude_txt -> tuple(sId, rId, vcf_gz, tbi, mask, ref_fa, exclude_txt, "") }

        masked_consensus = CONSENSUS_FASTA(cons_in).fasta

        // 7b. Virgin (unmasked) consensus from the ORIGINAL dedup BAM, in parallel with the masked one.
        if (asBool(params.keep_virgin_consensus) && asBool(params.dynamic_read_filter)) {
            def raw_bb_in = dedup_bam.map { sId, rId, bam, bai, ref_fa, excl -> tuple(sId, rId, bam, bai, ref_fa) }
            def bb_raw = CALL_BACKBONE_RAW(raw_bb_in)
            // Virgin SNPs: re-call FreeBayes on the raw bam (shows pre-filter variants) or reuse masked SNPs.
            def raw_snps = asBool(params.virgin_full_freebayes) ? CALL_FREEBAYES_RAW(dedup_bam).snps : fb_snps
            def merge_raw = raw_snps.join(bb_raw.backbone, by: [0,1]).join(bb_raw.header, by: [0,1])
                .map { sId, rId, sv, st, bv, bt, hdr -> tuple(sId, rId, sv, st, bv, bt, hdr, ".raw") }
            def vcf_raw = MERGE_VCFS_RAW(merge_raw)
            def refmeta_raw = dedup_bam.map { sId, rId, bam, bai, ref_fa, excl -> tuple(sId, rId, ref_fa, excl) }
            def cons_raw = vcf_raw.allpos.join(fb_mask_sites, by: [0,1]).join(refmeta_raw, by: [0,1])
                .map { sId, rId, vcf_gz, tbi, mask, ref_fa, excl -> tuple(sId, rId, vcf_gz, tbi, mask, ref_fa, excl, ".raw") }
            CONSENSUS_FASTA_RAW(cons_raw)
        }
    }

    emit:
    fasta = masked_consensus   // (sId, rId, consensus fasta); empty when make_consensus is off
}
