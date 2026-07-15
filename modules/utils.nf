// modules/utils.nf

/**
 * Helper function to clean strings (trim and remove carriage returns)
 */
def cleanStr(v) {
    if (v == null) return null
    v.toString().replace('\r','').trim()
}

/**
 * Helper function to check for "null-like" values (NA, ., null)
 */
def nullish(v) {
    if (v == null) return true
    def s = cleanStr(v)
    if (!s) return true
    def sl = s.toLowerCase()
    return (s == "." || sl == "na" || sl == "n/a" || sl == "null")
}

/**
 * Sanitizes IDs by replacing non-alphanumeric characters with underscores
 */
def sanitizeId(v) {
    def s = cleanStr(v)
    if (!s) return s
    s.replaceAll(/[^A-Za-z0-9_.-]+/, "_")
}

/**
 * Logic to calculate the Sample Directory path.
 * If params.nested_output is true, it splits the ID into a nested structure.
 * * Logic:
 * 1. Separates Letter Prefix from Numeric Suffix.
 * 2. Isolates the LAST digit of the number into its own folder.
 * 3. Splits the remaining numbers into chunks of 2.
 * * Example: MP00001  -> MP/00/00/1
 * Example: MIP00123 -> MIP/00/12/3
 */
def getSampleDir(sampleId, params) {
    if (params.nested_output) {
        // Regex to separate alphabetical prefix from numeric suffix
        def match = sampleId =~ /^([A-Za-z]+)(\d+)$/
        
        if (match.matches()) {
            def prefix = match[0][1] // e.g., "MIP"
            def number = match[0][2] // e.g., "00123"
            
            if (number.length() > 0) {
                // 1. Isolate the last digit
                def lastDigit = number[-1]
                
                // 2. Take the rest (everything except the last digit)
                def rest = (number.length() > 1) ? number[0..-2] : ""
                
                // 3. If there is a rest, split it into chunks of 2
                def nestedRest = ""
                if (rest) {
                    nestedRest = rest.replaceAll("(.{2})", "\$1/")
                    // Ensure it ends with / to join correctly
                    if (!nestedRest.endsWith("/")) nestedRest += "/"
                }
                
                // Result: Prefix / Nested_Rest / Last_Digit
                return "${prefix}/${nestedRest}${lastDigit}"
            }
        }
        
        // Fallback: If ID does not follow Letters+Numbers pattern (e.g., "Sample_A")
        // Just isolate the last character as a folder to maintain depth consistency
        if (sampleId.length() > 1) {
            def last = sampleId[-1]
            def rest = sampleId[0..-2]
            def nestedRest = rest.replaceAll("(.{2})", "\$1/")
            if (!nestedRest.endsWith("/")) nestedRest += "/"
            return "${nestedRest}${last}"
        }
    }
    // Default behavior: Return original ID
    return sampleId
}

/**
 * Centralized logic for 'saveAs' in publishDir.
 * Filters intermediate files and organizes output into stats/ and lineage/ folders.
 */
def getSavePath(filename, params) {
    def name  = filename.tokenize('/').last()
    def lower = name.toLowerCase()

    // A. INTERMEDIATE FILES TO ALWAYS BLOCK
    if (lower.endsWith('.fastq')    || lower.endsWith('.fastq.gz') || 
        lower.endsWith('.fq')       || lower.endsWith('.fq.gz') ||
        lower.contains('backbone')  ||
        lower.contains('header_template') ||
        lower.contains('valid_snps_formatted') ||
        lower == 'reference.fa' ||          
        lower.endsWith('.fai')  ||          
        lower.endsWith('.amb')  ||          
        lower.endsWith('.ann')  ||
        lower.endsWith('.bwt')  ||
        lower.endsWith('.pac')  ||
        lower.endsWith('.sa')) {
        return null
    }

    // A1b. With output_cram, CRAM is the published alignment form -> suppress the BAM/BAI
    //      (analysis still runs on the BAM in the work dir; only publishing changes).
    if (params.output_cram && (lower.endsWith('.bam') || lower.endsWith('.bam.bai'))) {
        return null
    }

    // A2. Pre-filter dedup alignment: once the length-aware filter runs, filtered.* is the analysis
    //     alignment and final.* is a ~redundant second full copy. Drop it from the outdir by
    //     default (still kept in the work dir, so the virgin/raw branch is unaffected). When the
    //     filter is off, final.* is the ONLY alignment and this gate does not trigger. Matches both
    //     .final.bam(.bai) and .final.cram(.crai).
    if ((lower.contains('.final.bam') || lower.contains('.final.cram'))
        && params.dynamic_read_filter && !params.publish_prefilter_bam) {
        return null
    }

    // A3. All-positions (per-position) VCFs are the consensus substrate. Gate publishing here
    //     (hoisted OUT of the annotate block so the flags work regardless of annotate_main_vcf).
    if (lower.contains('all.pos') && (lower.endsWith('.vcf.gz') || lower.endsWith('.vcf.gz.tbi'))) {
        if (!params.publish_allpos_vcf) return null
        if (lower.contains('.raw.') && !params.publish_virgin_allpos_vcf) return null
        return name
    }

    // B. ANNOTATION LOGIC (Filter Raw VCFs if annotation is enabled)
    if (params.annotate_legacy_vcfs) {
        if ((lower.contains('.var.') || lower.contains('freebayes.raw')) && !lower.contains('.ann.')) {
            return null
        }
    }

    if (params.annotate_main_vcf) {
        // all.pos is already handled above; here we only drop the un-annotated main/virgin VCFs
        // (their .ann. versions are the deliverables).
        if ((lower.endsWith('.vcf.gz') || lower.endsWith('.vcf.gz.tbi')) &&
            !lower.contains('.ann.') &&
            !lower.contains('.var.') &&
            !lower.contains('freebayes.raw')) {
            return null
        }
    }

    // C. STATS FOLDER
    if (lower.endsWith('.log')    || lower.endsWith('.stats') || 
        lower.endsWith('.json')   || lower.endsWith('.html')  || 
        lower.endsWith('.report') || lower.endsWith('.csv')   || 
        lower.endsWith('.tsv')    || lower.endsWith('.txt')) {
        return "stats/${name}"
    }
    
    // D. LINEAGE FOLDER (PATHOTYPR)
    if (lower.contains('pathotypr') || lower.endsWith('.lineage.tsv') || lower.endsWith('.split_kmer.tsv')) {
        return "lineage/${name}"
    }

    // E. DEFAULT PUBLISH
    return name 
}
