/*
===============================================================================
 Stage 11: Cohort-level outputs
===============================================================================
 The consolidated interactive QC report and the master SNP matrix.
*/

include { ANNOTATE_CANONICAL } from '../modules/annotation'
include { COLLECT_SUMMARY; COLLECT_DR; DEPTH_PROFILE; INDEL_MATRIX; LIFT_VARIANTS; QC_REPORT; SNP_DISTANCES; SNP_MATRIX } from '../modules/report'
include { MNV_TABLE } from '../modules/mnv'
include { asBool } from '../modules/utils'

workflow COHORT_REPORT {

    take:
    legacy_log         // per-sample GENERATE_LEGACY_STATS logs
    masked_consensus   // (sId, rId, consensus fasta); empty when make_consensus is off
    ref_bundle         // (refId, ref_fa, indices, exclude_txt)
    kraken_reports     // per-sample Kraken2 reports; empty when Kraken is off
    patho_dr_results   // (sId, dr mutations); empty when pathotypr is off
    vcf_for_stats      // (sId, rId, vcf)
    freebayes_ann      // annotated freebayes.raw VCFs; empty when annotate_legacy_vcfs is off
    allpos_vcf         // (sId, rId, all.pos.vcf.gz, tbi): the depth of every site, for the SNP matrix
    indels_vcf         // (sId, rId, indels(.ann).vcf.gz): every indel call, PASS or LowSupport
    mnv_tsv            // (sId, rId, <sample>.<ref>.MNV.tsv): get_MNV's codons; empty when run_mnv is off
    gconv_cohort       // cohort gene-conversion TSV; empty when --find_gene_conversion is off
    refMap             // refId -> refFasta
    refGffMap          // refId -> refGff
    tsv_name           // samplesheet basename, used to name the cohort-level outputs

    main:

    // 11. Cohort-level outputs (QC report + master SNP matrix).
    // Per-sample annotated VCFs, shared by both: prefer the annotated freebayes VCF (per-alt AD ->
    // continuous AF + gene/effect); fall back to the annotated main VCF (GT-based AF) if legacy
    // annotation is off. reference name(s) the samples were mapped against.
    def report_vcfs = (asBool(params.annotate_legacy_vcfs) ? freebayes_ann
                                                   : vcf_for_stats.map { sId, rId, vcf -> vcf })
                      .collect().ifEmpty([])
    def ref_name = refMap.keySet().join(',')

    // 11a0. The run's codon-level changes: every sample's get_MNV codons that reads carry whole. The
    // manifest names each file's sample and reference, so neither is read back out of a file name.
    def mnv_table = asBool(params.run_mnv)
        ? MNV_TABLE(mnv_tsv.map { sId, rId, f -> "${sId}\t${rId}\t${f.name}" }
                           .collectFile(name: 'mnv_manifest.tsv', newLine: true),
                    mnv_tsv.map { sId, rId, f -> f }.collect(), tsv_name).table
              .ifEmpty(file("${projectDir}/assets/NO_FILE_MNV"))
        : file("${projectDir}/assets/NO_FILE_MNV")

    // 11a. Master SNP matrix: rows = SNP sites, columns = reference/annotation + per-sample AF & depth.
    // The all-positions VCFs say what a sample's reads show at a site it has no call at. Built before
    // the report, which reads it to tell a site absent at a series' first time point from one not read.
    // The codon-level table says which SNPs are one amino-acid change together (codon_change columns).
    def snp_matrix = asBool(params.make_snp_matrix)
        ? SNP_MATRIX(report_vcfs, allpos_vcf.map { sId, rId, vcf, tbi -> vcf }.collect().ifEmpty([]),
                     mnv_table, ref_name, tsv_name).matrix
        : file("${projectDir}/assets/NO_FILE_MATRIX")

    // 11a'. Master indel matrix: every indel a sample calls with a PASS, each sample's AF, depth and filter.
    def indel_matrix = asBool(params.make_indel_matrix)
        ? INDEL_MATRIX(indels_vcf.map { sId, rId, vcf -> vcf }.collect().ifEmpty([]),
                       allpos_vcf.map { sId, rId, vcf, tbi -> vcf }.collect().ifEmpty([]),
                       ref_name, tsv_name).matrix
        : file("${projectDir}/assets/NO_FILE_INDELS")

    // 11b. What each sample's reads cover, from its all-positions VCF and its reference's GFF: the
    // report reads them against each other for deletions and SNPs per callable kb, and the tables
    // are published per sample whether or not a report is made.
    def depth = DEPTH_PROFILE(allpos_vcf.map { sId, rId, vcf, tbi -> tuple(sId, rId, vcf, file(refGffMap[rId])) })
    def depth_profiles = depth.profile.map { sId, rId, win, zero -> [win, zero] }
                              .mix(depth.genes.map { sId, rId, genes -> [genes] })
                              .flatten().collect().ifEmpty([])

    // 11c. Pairwise SNP distances between the consensus sequences, a deliverable of their own. The
    // manifest names each file's sample and reference, so neither has to be read back out of a file
    // name. Without consensus sequences SNP_DISTANCES never runs, and ifEmpty hands the report the
    // placeholder instead of leaving it waiting for an input that never comes.
    def cons_files = masked_consensus.map { sId, rId, fa -> fa }.collect().ifEmpty([])
    def snp_dist = asBool(params.make_snp_distances)
        ? SNP_DISTANCES(masked_consensus.map { sId, rId, fa -> "${sId}\t${rId}\t${fa.name}" }
                                        .collectFile(name: 'consensus_manifest.tsv', newLine: true),
                        cons_files, tsv_name).pairs
              .ifEmpty(file("${projectDir}/assets/NO_FILE_DISTANCES"))
        : file("${projectDir}/assets/NO_FILE_DISTANCES")

    // 11d. Consolidated QC report: aggregate every per-sample legacy log into one summary TSV, then
    // render the self-contained interactive HTML + PASS/WARN/FAIL flags.
    if (asBool(params.make_qc_report)) {
        def all_logs = legacy_log.collect()
        def summ = COLLECT_SUMMARY(all_logs, tsv_name)
        // Every reference's GFF, repeat/exclude list and FASTA index, in reference order: the report draws
        // each reference's own genes and masked regions under its own samples, placed on the axis their
        // profiles are binned on (the index's contig order). Taking the first reference's for all drew one
        // reference's genes under the other's samples.
        def ref_files = ref_bundle
            .map { rId, fa, idx, excl -> tuple(rId, file(refGffMap[rId]), excl, idx.find { it.name.endsWith('.fai') }) }
            .toSortedList { a, b -> a[0] <=> b[0] }
        def report_ref_ids = ref_files.map { l -> l.collect { it[0] } }
        def report_gffs    = ref_files.map { l -> l.collect { it[1] } }
        def report_masks   = ref_files.map { l -> l.collect { it[2] } }
        def report_fais    = ref_files.map { l -> l.collect { it[3] } }
        // Provenance footer: pinned container digest + reference(s).
        // The footer records what produced the numbers. When the gene-conversion panel is
        // showing, its settings belong there too: its verdicts depend on them, and a report
        // that displays a verdict without saying under which rules cannot be checked.
        def prov_items = ["container=${params.container}"] + refMap.keySet().collect { "reference=${it}" }
        if (asBool(params.find_gene_conversion)) {
            prov_items += ["gconv_min_bf=${params.gconv_min_bf}",
                           "gconv_prior=${params.gconv_prior}",
                           "gconv_mut_rate=${params.gconv_mut_rate}",
                           "gconv_indel_factor=${params.gconv_indel_factor}",
                           "gconv_min_tract_af=${params.gconv_min_tract_af}",
                           "gconv_ubiquitous=${params.gconv_ubiquitous}",
                           "gconv_donor_margin=${params.gconv_donor_margin}"]
        }
        def provenance  = prov_items.collect { "\"${it}\"" }.join(' ')
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
        def dr_report = asBool(params.run_pathotypr)
            ? COLLECT_DR(patho_dr_results.map { sId, f -> f }.collect().ifEmpty([]), tsv_name).dr
            : file("${projectDir}/assets/NO_FILE_DR")
        // Canonical (H37Rv) coordinates of the variant positions: one liftover per reference, each from its
        // own sequence and for the positions on its own contigs, joined into one map keyed by contig. The
        // report's SNP tables read it, and so does the canonical annotation below, which needs each record's
        // canonical position. Lifting every reference with the first one's sequence put 96% of another
        // reference's positions at a wrong coordinate.
        def lift_needed  = asBool(params.variant_liftover) || asBool(params.annotate_canonical)
        def pos_liftover = lift_needed
            ? LIFT_VARIANTS(ref_bundle.map { rId, fa, idx, excl -> tuple(rId, fa) }, report_vcfs).map
                  .collectFile(name: "${tsv_name}_pos_liftover.tsv", keepHeader: true)
                  .first()
            : file("${projectDir}/assets/NO_FILE_LIFTOVER")
        // Canonical-reference-annotated VCFs for the dual amino-acid numbering (off unless annotate_canonical).
        // Each record is moved to its canonical position before it is annotated and keeps its mapping
        // coordinate in OPOS, which is how the report pairs it with its variant.
        def report_vcfs_h37rv = asBool(params.annotate_canonical)
            ? ANNOTATE_CANONICAL(vcf_for_stats, pos_liftover, params.canonical_snpeff_db).out
                             .map { sId, vcf -> vcf }.collect().ifEmpty([])
            : file("${projectDir}/assets/NO_FILE_H37RV")
        // Kraken2 per-sample reports (deduped) -> Taxonomic composition panel; empty when Kraken is off.
        def report_kraken = kraken_reports.unique { it.name }.collect().ifEmpty([])
        // Gene conversion tracts (stage 12, opt-in --find_gene_conversion) -> the Gene conversion panel.
        // Empty when the stage is off, and QC_REPORT only passes --gene-conversion for a non-empty
        // file, so the panel hides itself rather than rendering an empty table.
        def report_gconv = gconv_cohort.ifEmpty([])
        QC_REPORT(summ.summary, summ.gene_burden, cons_files, report_gffs, report_masks, report_fais,
                  report_ref_ids, report_meta, report_vcfs, report_vcfs_h37rv, pos_liftover, dr_report,
                  report_gconv, report_kraken, depth_profiles, snp_dist, snp_matrix, indel_matrix, mnv_table,
                  provenance, tsv_name)
    }

}
