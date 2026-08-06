/*
===============================================================================
 Stage 10: MultiQC Report
===============================================================================
*/

include { MULTIQC; DUMP_VERSIONS } from '../modules/qc'

workflow MULTIQC_REPORT {

    take:
    pe_json          // (sId, fastp json) from FASTP_PE
    se_json          // (sId, fastp json) from FASTP_SE
    bam_stats        // (sId, rId, dedup.stats)
    kraken_reports   // per-sample Kraken2 reports; empty when Kraken is off
    snpeff_stats     // SnpEff CSV; empty when annotate_main_vcf is off
    filter_stats     // read-filter drop rates; empty when dynamic_read_filter is off
    report_filename  // "<samplesheet>_multiqc_report"

    main:

    // 10. MultiQC Report
    // Collect all relevant metrics from previous processes.
    // FastP JSONs and Kraken reports are reference-independent (QC of the raw reads), so a sample
    // mapped to >1 reference emits identically-named copies from parallel tasks. Deduplicate by
    // filename before collecting, otherwise MultiQC hits a fatal input-name collision.
    def fastp_json_mqc = pe_json.mix(se_json).map { sId, json -> json }.unique { it.name }
    def kraken_mqc     = kraken_reports.unique { it.name }
    // Provenance: dump the container's tool versions once and surface them in the report.
    def versions_mqc   = DUMP_VERSIONS().mqc
    def qc_collection = Channel.empty()
        .mix(fastp_json_mqc)                     // FastP
        .mix(bam_stats.map { s, r, st -> st })   // Samtools (drop the (sId,rId) key -> bare path)
        .mix(kraken_mqc)                         // Kraken
        .mix(snpeff_stats)                       // SnpEff
        .mix(filter_stats)                       // Length-aware read filter (drop rate)
        .mix(versions_mqc)                       // Software versions
        .collect()

    MULTIQC(qc_collection, report_filename)
}
