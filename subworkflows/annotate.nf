/*
===============================================================================
 Stage 8: Annotation
===============================================================================
 snpEff over the legacy split VCFs and over the main per-sample VCF.
*/

include { ANNOTATE_LEGACY_VCF; ANNOTATE_MAIN_VCF } from '../modules/annotation'

workflow ANNOTATE {

    take:
    homo_snp     // (sId, rId, var.homo.SNPs.vcf)
    het_snp      // (sId, rId, var.het.SNPs.vcf)
    homo_indel   // (sId, rId, var.homo.indel.vcf)
    raw_fb       // (sId, rId, freebayes.raw.vcf.gz, tbi)
    main_vcf     // (sId, rId, vcf.gz, tbi)
    allpos_vcf   // (sId, rId, all.pos.vcf.gz, tbi)
    snpeff_db    // (refId, snpEff.config, data)

    main:

    // 8. Annotation
    // Annotate Legacy split VCFs (Homo/Het/Indel/Raw)
    def freebayes_ann = Channel.empty()   // annotated freebayes.raw (AD + snpEff ANN) -> SNP dynamics
    if (params.annotate_legacy_vcfs) {
        def leg_inputs = Channel.empty()
            .mix(homo_snp.map   { sId, rId, vcf -> tuple(rId, sId, "var.homo.SNPs", vcf) })
            .mix(het_snp.map    { sId, rId, vcf -> tuple(rId, sId, "var.het.SNPs", vcf) })
            .mix(homo_indel.map { sId, rId, vcf -> tuple(rId, sId, "var.homo.indel", vcf) })
            .mix(raw_fb.map     { sId, rId, vcf, tbi -> tuple(rId, sId, "freebayes.raw", vcf) })

        def ann_leg_in = leg_inputs.combine(snpeff_db, by: 0)
            .map { rId, sId, lbl, vcf, cfg, dat -> tuple(sId, rId, lbl, vcf, cfg, dat) }

        def ann_leg = ANNOTATE_LEGACY_VCF(ann_leg_in)
        // Keep the annotated freebayes VCF: it carries per-alt AD (true AF) AND gene/effect (ANN),
        // so the dynamics panel gets continuous allele frequencies without losing gene annotation.
        freebayes_ann = ann_leg.out
            .filter { sId, rId, vcf, tbi -> vcf.name.contains('freebayes.raw') }
            .map { sId, rId, vcf, tbi -> vcf }
    }

    // Annotate Main VCF
    def vcf_for_stats = Channel.empty()
    def ch_snpeff_stats = Channel.empty()

    if (params.annotate_main_vcf) {
        def main_pre = main_vcf.map { sId, rId, vcf, tbi -> tuple(rId, sId, vcf, tbi) }
        def ann_main_in = main_pre.combine(snpeff_db, by: 0)
            .map { rId, sId, vcf, tbi, cfg, dat -> tuple(sId, rId, vcf, tbi, cfg, dat) }

        def ann_out = ANNOTATE_MAIN_VCF(ann_main_in)

        vcf_for_stats = ann_out.vcf_ann.map { sId, rId, vcf, tbi -> tuple(sId, rId, vcf) }
        ch_snpeff_stats = ann_out.stats // Collect SnpEff stats CSV for MultiQC
    } else {
        vcf_for_stats = allpos_vcf.map { sId, rId, vcf, tbi -> tuple(sId, rId, vcf) }
    }

    emit:
    freebayes_vcf = freebayes_ann     // annotated freebayes.raw; empty when annotate_legacy_vcfs is off
    stats_vcf     = vcf_for_stats     // (sId, rId, vcf) feeding the legacy stats and the cohort report
    snpeff_stats  = ch_snpeff_stats   // SnpEff CSV for MultiQC; empty when annotate_main_vcf is off
}
