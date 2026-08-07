/*
===============================================================================
 Stage 11: Cohort-level outputs
===============================================================================
 The consolidated interactive QC report and the master SNP matrix.
*/

include { ANNOTATE_CANONICAL } from '../modules/annotation'
include { COLLECT_SUMMARY; COLLECT_DR; LIFT_VARIANTS; QC_REPORT; SNP_MATRIX } from '../modules/report'

workflow COHORT_REPORT {

    take:
    legacy_log         // per-sample GENERATE_LEGACY_STATS logs
    masked_consensus   // (sId, rId, consensus fasta); empty when make_consensus is off
    ref_bundle         // (refId, ref_fa, indices, exclude_txt)
    kraken_reports     // per-sample Kraken2 reports; empty when Kraken is off
    patho_dr_results   // (sId, dr mutations); empty when pathotypr is off
    vcf_for_stats      // (sId, rId, vcf)
    freebayes_ann      // annotated freebayes.raw VCFs; empty when annotate_legacy_vcfs is off
    gconv_cohort       // cohort gene-conversion TSV; empty when --find_gene_conversion is off
    refMap             // refId -> refFasta
    refGffMap          // refId -> refGff
    tsv_name           // samplesheet basename, used to name the cohort-level outputs

    main:

    // 11. Cohort-level outputs (QC report + master SNP matrix).
    // Per-sample annotated VCFs, shared by both: prefer the annotated freebayes VCF (per-alt AD ->
    // continuous AF + gene/effect); fall back to the annotated main VCF (GT-based AF) if legacy
    // annotation is off. reference name(s) the samples were mapped against.
    def report_vcfs = (params.annotate_legacy_vcfs ? freebayes_ann
                                                   : vcf_for_stats.map { sId, rId, vcf -> vcf })
                      .collect().ifEmpty([])
    def ref_name = refMap.keySet().join(',')

    // 11a. Consolidated QC report: aggregate every per-sample legacy log into one summary TSV, then
    // render the self-contained interactive HTML + PASS/WARN/FAIL flags.
    if (params.make_qc_report) {
        def all_logs = legacy_log.collect()
        def summ = COLLECT_SUMMARY(all_logs, tsv_name)
        def cons_files = masked_consensus.map { sId, rId, fa -> fa }.collect().ifEmpty([])
        // Reference-level extras (cohort report -> take the reference bundle; single-ref is the norm):
        // GFF enables the per-gene SNP-density panel, the nucmer/repeat BED the masked-regions panel.
        def report_ref  = refGffMap.keySet().toList().first()   // deterministic: the first reference
        def report_gff  = file(refGffMap[report_ref])
        def report_mask = ref_bundle.filter { rId, fa, idx, excl -> rId == report_ref }
                                    .map { rId, fa, idx, excl -> excl }.first()
        // Provenance footer: pinned container digest + reference(s).
        def provenance  = (["container=${params.container}"] + refMap.keySet().collect { "reference=${it}" })
                          .collect { "\"${it}\"" }.join(' ')
        // SNP dynamics: the samplesheet is the metadata source (auto-detects time/group columns; the
        // panel self-hides if absent).
        def report_meta = file(params.tsv)
        // Optional QC_REPORT inputs use a per-input placeholder from assets/. They must have DISTINCT
        // names: Nextflow refuses to stage two inputs of the same task under one filename, and all three
        // of these are off by default, so a single shared "NO_FILE" fails the run with an input file name
        // collision. They are real empty files (a dangling symlink breaks under stageInMode 'copy') and
        // are referenced through projectDir so a stray NO_FILE in the launch directory is never picked up.
        // QC_REPORT recognises them by the NO_FILE* prefix.
        // Drug-resistance calls (pathotypr DR run -> one run TSV); placeholder when pathotypr is off.
        def dr_report = params.run_pathotypr
            ? COLLECT_DR(patho_dr_results.map { sId, f -> f }.collect().ifEmpty([]), tsv_name).dr
            : file("${projectDir}/assets/NO_FILE_DR")
        // Canonical-reference-annotated VCFs for the dual amino-acid numbering (off unless annotate_canonical).
        // vcf_for_stats is the per-sample (sId,rId,vcf) channel; its positions match report_vcfs, so the
        // report's sample+position merge finds each variant's canonical amino-acid change.
        def report_vcfs_h37rv = params.annotate_canonical
            ? ANNOTATE_CANONICAL(vcf_for_stats, params.canonical_snpeff_db).out
                             .map { sId, vcf -> vcf }.collect().ifEmpty([])
            : file("${projectDir}/assets/NO_FILE_H37RV")
        // Kraken2 per-sample reports (deduped) -> Taxonomic composition panel; empty when Kraken is off.
        def report_kraken = kraken_reports.unique { it.name }.collect().ifEmpty([])
        // Gene conversion tracts (stage 12, opt-in --find_gene_conversion) -> the Gene conversion panel.
        // Empty when the stage is off, and QC_REPORT only passes --gene-conversion for a non-empty
        // file, so the panel hides itself rather than rendering an empty table.
        def report_gconv = gconv_cohort.ifEmpty([])
        // Alignment-free canonical COORDINATE per variant via pathotypr (alternative to the --vcfs-h37rv path):
        // lift the run's variant positions (mapping-reference coords) onto H37Rv and hand the map to the report.
        def report_ref_fa = ref_bundle.filter { rId, fa, idx, excl -> rId == report_ref }
                                      .map { rId, fa, idx, excl -> fa }.first()
        def pos_liftover  = params.variant_liftover
            ? LIFT_VARIANTS(report_vcfs, report_ref_fa, report_ref, tsv_name).map
            : file("${projectDir}/assets/NO_FILE_LIFTOVER")
        QC_REPORT(summ.summary, summ.gene_burden, cons_files, report_gff, report_mask,
                  report_meta, report_vcfs, report_vcfs_h37rv, pos_liftover, dr_report, report_gconv,
                  report_kraken, provenance, tsv_name)
    }

    // 11b. Master SNP matrix: rows = SNP sites, columns = reference/annotation + per-sample AF & depth.
    if (params.make_snp_matrix) {
        SNP_MATRIX(report_vcfs, ref_name, tsv_name)
    }
}
