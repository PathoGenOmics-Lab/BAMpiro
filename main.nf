#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
===============================================================================
 BAMpiro Pipeline - Short Read Mapping, Variant Calling, Lineage & Reporting
===============================================================================
*/

// --- MODULE IMPORTS ---
include { PREPARE_REFERENCE; SNPEFF_BUILD_DB; BUILD_MAPPABILITY } from './modules/reference'
include { VALIDATE_RAW_READS_PE; VALIDATE_RAW_READS_SE; KRAKEN_FILTER_PE; KRAKEN_FILTER_SE; FASTP_PE; FASTP_SE; MULTIQC } from './modules/qc'
include { RUN_PATHOTYPR_PE; RUN_PATHOTYPR_SE } from './modules/pathotypr'
include { MAPPING_PE; MAPPING_SE; MERGE_AND_MARKDUP; FILTER_READS } from './modules/mapping'
include { CALL_FREEBAYES; CALL_BACKBONE; MERGE_VCFS; CALL_FREEBAYES_RAW } from './modules/variants'
include { CONSENSUS_FASTA } from './modules/consensus'
// Aliases for the parallel virgin (unmasked) consensus path (a DSL2 process runs once per name)
include { CALL_BACKBONE as CALL_BACKBONE_RAW; MERGE_VCFS as MERGE_VCFS_RAW } from './modules/variants'
include { CONSENSUS_FASTA as CONSENSUS_FASTA_RAW } from './modules/consensus'
include { ANNOTATE_LEGACY_VCF; ANNOTATE_MAIN_VCF; GENERATE_LEGACY_STATS } from './modules/annotation'

/* ----------------------------- Configuration Logic ----------------------------- */

def final_outdir = file(params.outdir).toAbsolutePath().toString()

// Logic to extract the base name of the TSV (e.g., "samples.tsv" -> "samples")
// This is used to name the final MultiQC report dynamically
def tsv_name = file(params.tsv).simpleName
def multiqc_report_filename = "${tsv_name}_multiqc_report"

/* ----------------------------- Helpers ----------------------------- */

def cleanStr = { v ->
    if (v == null) return null
    v.toString().replace('\r','').trim()
}

def nullish = { v ->
    if (v == null) return true
    def s = cleanStr(v)
    if (!s) return true
    def sl = s.toLowerCase()
    return (s == "." || sl == "na" || sl == "n/a" || sl == "null")
}

def sanitizeId = { v ->
    def s = cleanStr(v)
    if (!s) return s
    s.replaceAll(/[^A-Za-z0-9_.-]+/, "_")
}

def inferRunId = { r1 ->
    def name = new File(r1).getName()
    name = name.replaceAll(/(\.fastq|\.fq)(\.gz|\.bz2)?$/,'')
    name = name.replaceAll(/(\.gz|\.bz2)$/,'')
    sanitizeId(name)
}

// 1. Output Directory Strategy
def SAMPLE_PUBLISH_DIR = { sid -> "${final_outdir}/${sid}" }

/* ----------------------------- TSV Parsing ----------------------------- */

def tsvFile = new File(params.tsv as String)
if (!tsvFile.exists()) throw new RuntimeException("TSV not found: ${params.tsv}")

def lines = tsvFile.readLines()
if (lines.size() < 2) throw new RuntimeException("TSV has no data rows: ${params.tsv}")

def header = lines[0].replace('\r','').split('\t')*.trim()
def col = [:]
header.eachWithIndex { h,i -> col[h]=i }

['sampleId','r1','refId','refFasta','refGff'].each { req ->
    if (!col.containsKey(req)) throw new RuntimeException("TSV missing required column: ${req}")
}

def refMap    = [:]
def refGffMap = [:] 
def expectedMap = [:].withDefault{0}
def seenRuns = [:].withDefault{0}
def peList = []
def seList = []
def hasTax = false

lines.drop(1).each { raw ->
    def line = raw.replace('\r','')
    if (!line.trim()) return
    if (line.startsWith('#')) return

    def p = line.split('\t', -1)

    def sampleId = sanitizeId(p[col.sampleId])
    def r1Str    = cleanStr(p[col.r1])
    def r2Str    = col.containsKey('r2') ? cleanStr(p[col.r2]) : null
    def refId    = sanitizeId(p[col.refId])
    def refFasta = cleanStr(p[col.refFasta])
    def refGff   = cleanStr(p[col.refGff])
    def taxId    = col.containsKey('taxId') ? cleanStr(p[col.taxId]) : null
    def runId    = col.containsKey('runId') ? cleanStr(p[col.runId]) : null

    if (!sampleId) throw new RuntimeException("Empty sampleId in TSV line: ${raw}")
    if (!r1Str)    throw new RuntimeException("Empty r1 for sampleId=${sampleId}")
    if (!refId)    throw new RuntimeException("Empty refId for sampleId=${sampleId}")
    if (!refFasta) throw new RuntimeException("Empty refFasta for sampleId=${sampleId}")
    if (!refGff)   throw new RuntimeException("Empty refGff for sampleId=${sampleId}")

    if (nullish(r2Str)) r2Str = null
    def mode = (r2Str ? "PE" : "SE")

    if (taxId != null && (taxId == "" || taxId == ".")) taxId = null
    if (taxId != null) hasTax = true

    if (!runId) runId = inferRunId(r1Str)
    runId = sanitizeId(runId)
    if (!runId) runId = "run"

    // Guarantee a unique run token per (sampleId, refId): lane-split inputs with no runId column
    // infer the same basename for every lane -> identical BAM names collide in the merge group.
    def runKey = "${sampleId}||${refId}||${runId}"
    def runN = (seenRuns[runKey] = seenRuns[runKey] + 1)
    if (runN > 1) runId = "${runId}_${runN}"

    if (!new File(r1Str).exists()) throw new RuntimeException("R1 not found: ${r1Str}")
    if (mode == "PE" && !new File(r2Str).exists()) throw new RuntimeException("R2 not found: ${r2Str}")
    if (!new File(refFasta).exists()) throw new RuntimeException("refFasta not found: ${refFasta}")
    if (!new File(refGff).exists())   throw new RuntimeException("refGff not found: ${refGff}")

    if (!refMap.containsKey(refId)) {
        refMap[refId] = refFasta
        refGffMap[refId] = refGff
    } else {
        if (refMap[refId] != refFasta) throw new RuntimeException("refId ${refId} has multiple refFasta paths")
    }

    def key = "${sampleId}||${refId}"
    expectedMap[key] = expectedMap[key] + (mode == "PE" ? 2 : 1)

    if (mode == "PE") {
        peList << [sampleId, runId, file(r1Str), file(r2Str), refId, taxId]
    } else {
        seList << [sampleId, runId, file(r1Str), refId, taxId]
    }
}

// Check if Kraken DB exists
def KRAKEN_ENABLED = false
if (params.kraken2_db && hasTax) {
    def db = new File(params.kraken2_db as String)
    if (db.exists()) KRAKEN_ENABLED = true
    else log.warn "Kraken DB not found at: ${params.kraken2_db} -> Kraken disabled"
}

// Check consensus script
if (params.make_consensus) {
    // Note: Script is now looked for in 'bin/', handled by Nextflow automatically
}

log.info """
================================================================
 BAMpiro Pipeline ðŸ§›â€â™‚ï¸
================================================================
TSV              : ${params.tsv}
Output Absolute  : ${final_outdir}
Report Name      : ${multiqc_report_filename}.html
Annotate Legacy  : ${params.annotate_legacy_vcfs}
Run Pathotypr    : ${params.run_pathotypr}
================================================================
"""

/* ----------------------------- Channels ----------------------------- */

def uniqueRefsList = refMap.collect { k,v -> [k, file(v)] }
def ref_in_ch   = Channel.fromList(uniqueRefsList).map { refId, fasta -> tuple(refId, fasta) }

def snpeff_db_in = Channel.fromList(refMap.keySet().collect { rid -> 
    [rid, file(refMap[rid]), file(refGffMap[rid])] 
}).map { rid, fa, gff -> tuple(rid, fa, gff) }

def pe_reads_ch = Channel.fromList(peList).map { sId, runId, r1, r2, refId, taxId -> tuple(sId, runId, r1, r2, refId, taxId) }
def se_reads_ch = Channel.fromList(seList).map { sId, runId, r1, refId, taxId -> tuple(sId, runId, r1, refId, taxId) }

/* ----------------------------- Main Workflow ----------------------------- */

workflow {

    // 1. Reference Preparation
    def ref_bundle = PREPARE_REFERENCE(ref_in_ch)
    def snpeff_db  = SNPEFF_BUILD_DB(snpeff_db_in)

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
    
    // 4. Pathotypr (Optional Lineage Classification)
    def patho_results = Channel.empty()

    if (params.run_pathotypr) {
        // Reshape channels for Pathotypr
        def patho_pe_ch = fastp_pe.pe_reads.map { sId, runId, r1, r2, refId, taxId -> tuple(refId, sId, runId, r1, r2) }
        def patho_se_ch = fastp_se.se_reads.map { sId, runId, rSE, refId, taxId -> tuple(refId, sId, runId, rSE) }

        def run_pe = RUN_PATHOTYPR_PE(
            patho_pe_ch,
            file(params.pathotypr_ref),
            file(params.pathotypr_markers),
            params.pathotypr_bin
        )

        def run_se = RUN_PATHOTYPR_SE(
            patho_se_ch,
            file(params.pathotypr_ref),
            file(params.pathotypr_markers),
            params.pathotypr_bin
        )
        
        // Combine summaries for later statistics
        patho_results = run_pe.summary.mix(run_se.summary)
    }

    // 5. Mapping (BWA)
    def map_pe  = fastp_pe.pe_reads.map { sId, runId, r1, r2, refId, taxId -> tuple(refId, "PE", sId, runId, r1, r2, taxId) }
    def map_se1 = fastp_pe.se_reads.map { sId, runId, rSE, refId, taxId -> tuple(refId, "SE", sId, runId, rSE, null, taxId) }
    def map_se2 = fastp_se.se_reads.map { sId, runId, rSE, refId, taxId -> tuple(refId, "SE", sId, runId, rSE, null, taxId) }
    def map_all = map_pe.mix(map_se1).mix(map_se2)

    // Join reads with references (matching RefID)
    def joined = map_all.combine(ref_bundle.bundle).filter { it[0] == it[7] }.map { refId, mode, sId, runId, r1, r2, taxId, refId2, ref_fa, indices, exclude_txt -> tuple(refId, mode, sId, runId, r1, r2, taxId, ref_fa, indices, exclude_txt) }
    def branched = joined.branch { pe: it[1] == "PE"; se: it[1] == "SE" }
    
    def map_pe_in = branched.pe.map { refId, mode, sId, runId, r1, r2, taxId, ref_fa, indices, exclude_txt -> tuple(refId, sId, runId, r1, r2, taxId, ref_fa, indices, exclude_txt) }
    def map_se_in = branched.se.map { refId, mode, sId, runId, r1, r2, taxId, ref_fa, indices, exclude_txt -> tuple(refId, sId, runId, r1, taxId, ref_fa, indices, exclude_txt) }

    def mapped_pe = MAPPING_PE(map_pe_in)
    def mapped_se = MAPPING_SE(map_se_in)
    def all_bams = mapped_pe.bam.mix(mapped_se.bam)

    // Group BAMs by SampleID + RefID for merging
    def bams_grouped = all_bams.map { sId, rId, bam, ref_fa, exclude_txt ->
            def key = "${sId}||${rId}"
            def nExp = expectedMap[key]
            tuple(groupKey(key, nExp as int), sId, rId, bam, ref_fa, exclude_txt)
        }.groupTuple(by: 0).map { gk, sIds, rIds, bams, ref_fas, excls -> tuple(sIds[0], rIds[0], bams, ref_fas[0], excls[0]) }

    def final_bams = MERGE_AND_MARKDUP(bams_grouped)

    // 6. Variant Calling
    // Optional length-aware read filter: drop reads too short to map uniquely at their locus
    // (short-read false positives in near-repeats), then call/consensus on the filtered BAM.
    // Feature off -> call on the dedup BAM unchanged.
    def vbase
    def filter_stats = Channel.empty()
    if (params.dynamic_read_filter) {
        def mapp = BUILD_MAPPABILITY( ref_bundle.bundle.map { rId, fa, idx, excl -> tuple(rId, fa) } )
        def track_bed = mapp.track.join(mapp.repeat_bed, by: 0)                    // (rId, npz, repeat_bed)
        def filt_in = final_bams.final_bam
            .map { sId, rId, bam, bai, fa, excl -> tuple(rId, sId, bam, bai, fa, excl) }
            .combine(track_bed, by: 0)                                             // fan the per-reference track+bed to each sample
            .map { rId, sId, bam, bai, fa, excl, npz, bed -> tuple(sId, rId, bam, bai, fa, excl, npz, bed) }
        def filt = FILTER_READS(filt_in)
        vbase = filt.filtered_bam                                                  // exclude_txt now = nucmer + genmap repeats
        filter_stats = filt.stats
    } else {
        vbase = final_bams.final_bam
    }

    // Run FreeBayes and Backbone parallel
    def fb_out = CALL_FREEBAYES(vbase)
    def bb_in  = vbase.map { sId, rId, bam, bai, ref_fa, exclude_txt -> tuple(sId, rId, bam, bai, ref_fa) }
    def bb_out = CALL_BACKBONE(bb_in)
    
    // Merge specific variants with backbone for All-Positions VCF (outLabel "" = masked path)
    def merge_in = fb_out.snps.join(bb_out.backbone, by: [0,1]).join(bb_out.header, by: [0,1])
        .map { sId, rId, sv, st, bv, bt, hdr -> tuple(sId, rId, sv, st, bv, bt, hdr, "") }
    def vcf_ch = MERGE_VCFS(merge_in)

    // 7. Consensus Generation (Optional)
    if (params.make_consensus) {
        // Prepare inputs: VCF + Reference + Mask sites
        // Note: Script is called from bin/ directly in the module
        // refmeta from vbase so the masked consensus gets the combined (nucmer + genmap) exclude
        def refmeta = vbase.map { sId, rId, bam, bai, ref_fa, exclude_txt -> tuple(sId, rId, ref_fa, exclude_txt) }
        def allpos_mask = vcf_ch.allpos.join(fb_out.mask_sites, by: [0,1])
        def cons_in = allpos_mask.join(refmeta, by: [0,1]).map { sId, rId, vcf_gz, tbi, mask, ref_fa, exclude_txt -> tuple(sId, rId, vcf_gz, tbi, mask, ref_fa, exclude_txt, "") }

        CONSENSUS_FASTA(cons_in)

        // 7b. Virgin (unmasked) consensus from the ORIGINAL dedup BAM, in parallel with the masked one.
        if (params.keep_virgin_consensus && params.dynamic_read_filter) {
            def raw_bb_in = final_bams.final_bam.map { sId, rId, bam, bai, ref_fa, excl -> tuple(sId, rId, bam, bai, ref_fa) }
            def bb_raw = CALL_BACKBONE_RAW(raw_bb_in)
            // Virgin SNPs: re-call FreeBayes on the raw bam (shows pre-filter variants) or reuse masked SNPs.
            def raw_snps = params.virgin_full_freebayes ? CALL_FREEBAYES_RAW(final_bams.final_bam).snps : fb_out.snps
            def merge_raw = raw_snps.join(bb_raw.backbone, by: [0,1]).join(bb_raw.header, by: [0,1])
                .map { sId, rId, sv, st, bv, bt, hdr -> tuple(sId, rId, sv, st, bv, bt, hdr, ".raw") }
            def vcf_raw = MERGE_VCFS_RAW(merge_raw)
            def refmeta_raw = final_bams.final_bam.map { sId, rId, bam, bai, ref_fa, excl -> tuple(sId, rId, ref_fa, excl) }
            def cons_raw = vcf_raw.allpos.join(fb_out.mask_sites, by: [0,1]).join(refmeta_raw, by: [0,1])
                .map { sId, rId, vcf_gz, tbi, mask, ref_fa, excl -> tuple(sId, rId, vcf_gz, tbi, mask, ref_fa, excl, ".raw") }
            CONSENSUS_FASTA_RAW(cons_raw)
        }
    }

    // 8. Annotation
    // Annotate Legacy split VCFs (Homo/Het/Indel/Raw)
    if (params.annotate_legacy_vcfs) {
        def leg_inputs = Channel.empty()
            .mix(fb_out.homo_snp.map   { sId, rId, vcf -> tuple(rId, sId, "var.homo.SNPs", vcf) })
            .mix(fb_out.het_snp.map    { sId, rId, vcf -> tuple(rId, sId, "var.het.SNPs", vcf) })
            .mix(fb_out.homo_indel.map { sId, rId, vcf -> tuple(rId, sId, "var.homo.indel", vcf) })
            .mix(fb_out.raw_fb.map     { sId, rId, vcf, tbi -> tuple(rId, sId, "freebayes.raw", vcf) })
        
        def ann_leg_in = leg_inputs.combine(snpeff_db.db, by: 0)
            .map { rId, sId, lbl, vcf, cfg, dat -> tuple(sId, rId, lbl, vcf, cfg, dat) }
        
        ANNOTATE_LEGACY_VCF(ann_leg_in)
    }

    // Annotate Main VCF
    def vcf_for_stats = Channel.empty()
    def ch_snpeff_stats = Channel.empty()
    
    if (params.annotate_main_vcf) {
        def main_pre = vcf_ch.main_vcf.map { sId, rId, vcf, tbi -> tuple(rId, sId, vcf, tbi) }
        def ann_main_in = main_pre.combine(snpeff_db.db, by: 0)
            .map { rId, sId, vcf, tbi, cfg, dat -> tuple(sId, rId, vcf, tbi, cfg, dat) }
        
        def ann_out = ANNOTATE_MAIN_VCF(ann_main_in)
        
        vcf_for_stats = ann_out.vcf_ann.map { sId, rId, vcf, tbi -> tuple(sId, rId, vcf) }
        ch_snpeff_stats = ann_out.stats // Collect SnpEff stats CSV for MultiQC
    } else {
        vcf_for_stats = vcf_ch.allpos.map { sId, rId, vcf, tbi -> tuple(sId, rId, vcf) }
    }

    // 9. Legacy Stats Generation
    // Prepare Reference Index channel
    def ref_fai_ch = ref_bundle.bundle.map { rId, fa, indices, excl -> 
        def fai_file = indices.find { it.name.endsWith('.fai') }
        tuple(rId, fai_file)
    }

    // Prepare FastP JSONs (grouped by sample)
    def json_ch = fastp_pe.json.mix(fastp_se.json)
        .map { json -> 
            def meta = json.name.split('__') 
            tuple(meta[0], json) // [sId, json]
        }
        .groupTuple()
        .map { sId, jsons -> tuple(sId, jsons[0]) } // Take first JSON if multiple

    // BAM stats already carry (sampleId, refId) as vals -- do NOT re-parse the filename with
    // tokenize('.'), which truncates any dotted refId (e.g. NC_000962.3) and breaks the join.
    def bam_stats_ch = final_bams.stats // [sId, rId, stats]

    // Join logic for Stats
    def vcf_bam_joined = vcf_for_stats.join(bam_stats_ch, by: [0,1]) // [sId, rId, vcf, stats]
    // combine (not join) on sampleId: one json per sample must fan out to EVERY (sample,ref)
    // row, else a sample mapped to >1 reference silently loses all but one reference's stats.
    def vcf_bam_fastp = vcf_bam_joined.combine(json_ch, by: 0) // [sId, rId, vcf, stats, json]

    // Join with Reference Index (Using COMBINE to reuse reference)
    def ready_no_patho = vcf_bam_fastp
        .map { sId, rId, vcf, stats, json -> tuple(rId, sId, vcf, stats, json) }
        .combine(ref_fai_ch)
        .filter { it[0] == it[5] } // Filter where sample RefID == FAI RefID
        .map { rId, sId, vcf, stats, json, rId2, fai -> tuple(sId, rId, json, stats, vcf, fai) }

    // Join with Pathotypr (Using remainder to handle missing files safely)
    def final_stats_input = ready_no_patho
        .join(patho_results, by: 0, remainder: true) 
        .map { list ->
            // list size < 7 means missing left-side data (e.g. mapping failed)
            if (list.size() < 7) return null
            
            def sId   = list[0]
            def rId   = list[1]
            def json  = list[2]
            def stats = list[3]
            def vcf   = list[4]
            def fai   = list[5]
            def patho = list[6] // Can be null

            // Fallback placeholder file
            def real_patho = patho ? patho : file("NO_FILE")
            tuple(sId, rId, json, stats, vcf, fai, real_patho)
        }
        .filter { it != null }

    
    GENERATE_LEGACY_STATS(final_stats_input)

    // 10. MultiQC Report
    // Collect all relevant metrics from previous processes.
    // FastP JSONs and Kraken reports are reference-independent (QC of the raw reads), so a sample
    // mapped to >1 reference emits identically-named copies from parallel tasks. Deduplicate by
    // filename before collecting, otherwise MultiQC hits a fatal input-name collision.
    def fastp_json_mqc = fastp_pe.json.mix(fastp_se.json).unique { it.name }
    def kraken_mqc     = ch_kraken_reports.unique { it.name }
    def qc_collection = Channel.empty()
        .mix(fastp_json_mqc)                          // FastP
        .mix(final_bams.stats.map { s, r, st -> st }) // Samtools (drop the (sId,rId) key -> bare path)
        .mix(kraken_mqc)                              // Kraken
        .mix(ch_snpeff_stats)                         // SnpEff
        .mix(filter_stats)                            // Length-aware read filter (drop rate)
        .collect()

    MULTIQC(qc_collection, multiqc_report_filename)
}
