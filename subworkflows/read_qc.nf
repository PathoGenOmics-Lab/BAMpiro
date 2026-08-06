/*
===============================================================================
 Stages 2-3: Read Validation, optional Kraken filtering and FastP QC
===============================================================================
*/

include { VALIDATE_RAW_READS_PE; VALIDATE_RAW_READS_SE; KRAKEN_FILTER_PE; KRAKEN_FILTER_SE; FASTP_PE; FASTP_SE } from '../modules/qc'

workflow READ_QC {

    take:
    pe_reads_ch      // (sId, runId, r1, r2, refId, taxId)
    se_reads_ch      // (sId, runId, r1, refId, taxId)
    KRAKEN_ENABLED   // boolean: a usable Kraken DB AND at least one taxId in the samplesheet

    main:

    // 2. Read Validation
    def checked_pe = VALIDATE_RAW_READS_PE(pe_reads_ch)
    def checked_se = VALIDATE_RAW_READS_SE(se_reads_ch)
    def pe_final = checked_pe.reads
    def se_final = checked_se.reads

    // Channel to collect Kraken reports if enabled
    def ch_kraken_reports = Channel.empty()

    if (KRAKEN_ENABLED) {
        def pe_with_tax = checked_pe.reads.filter { it[5] != null }
        def pe_no_tax   = checked_pe.reads.filter { it[5] == null }
        def se_with_tax = checked_se.reads.filter { it[4] != null }
        def se_no_tax   = checked_se.reads.filter { it[4] == null }

        def krk_pe = KRAKEN_FILTER_PE(pe_with_tax, params.kraken2_db)
        def krk_se = KRAKEN_FILTER_SE(se_with_tax, params.kraken2_db)

        pe_final = krk_pe.reads.mix(pe_no_tax)
        se_final = krk_se.reads.mix(se_no_tax)

        // Accumulate reports for MultiQC
        ch_kraken_reports = krk_pe.report.mix(krk_se.report)
    }

    // 3. FastP QC
    def fastp_pe = FASTP_PE(pe_final)
    def fastp_se = FASTP_SE(se_final)

    emit:
    pe_reads       = fastp_pe.pe_reads   // (sId, runId, r1, r2, refId, taxId)
    pe_orphans     = fastp_pe.se_reads   // (sId, runId, rSE, refId, taxId): PE mates orphaned by FastP
    se_reads       = fastp_se.se_reads   // (sId, runId, rSE, refId, taxId)
    pe_json        = fastp_pe.json       // (sId, fastp json)
    se_json        = fastp_se.json       // (sId, fastp json)
    kraken_reports = ch_kraken_reports   // empty when Kraken is disabled
}
