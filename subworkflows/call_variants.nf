/*
===============================================================================
 Stage 6: Variant Calling
===============================================================================
 FreeBayes for the variant sites, the backbone for every other position; the
 two are merged into the all-positions VCF.
*/

include { CALL_FREEBAYES; CALL_BACKBONE; MERGE_VCFS } from '../modules/variants'

workflow CALL_VARIANTS {

    take:
    vbase   // (sId, rId, bam, bai, ref_fa, exclude_txt)

    main:

    // 6. Variant Calling
    // Run FreeBayes and Backbone parallel
    def fb_out = CALL_FREEBAYES(vbase)
    def bb_in  = vbase.map { sId, rId, bam, bai, ref_fa, exclude_txt -> tuple(sId, rId, bam, bai, ref_fa) }
    def bb_out = CALL_BACKBONE(bb_in)

    // Merge specific variants with backbone for All-Positions VCF (outLabel "" = masked path)
    def merge_in = fb_out.snps.join(bb_out.backbone, by: [0,1]).join(bb_out.header, by: [0,1])
        .map { sId, rId, sv, st, bv, bt, hdr -> tuple(sId, rId, sv, st, bv, bt, hdr, "") }
    def vcf_ch = MERGE_VCFS(merge_in)

    emit:
    snps       = fb_out.snps         // (sId, rId, valid_snps.vcf.gz, tbi)
    mask_sites = fb_out.mask_sites   // (sId, rId, mask_sites.tsv)
    raw_fb     = fb_out.raw_fb       // (sId, rId, freebayes.raw.vcf.gz, tbi)
    homo_snp   = fb_out.homo_snp     // legacy split VCFs
    het_snp    = fb_out.het_snp
    homo_indel = fb_out.homo_indel
    allpos     = vcf_ch.allpos       // (sId, rId, all.pos.vcf.gz, tbi)
    main_vcf   = vcf_ch.main_vcf     // (sId, rId, vcf.gz, tbi)
}
