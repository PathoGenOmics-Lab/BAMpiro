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

    // B. ANNOTATION LOGIC (Filter Raw VCFs if annotation is enabled)
    if (params.annotate_legacy_vcfs) {
        if ((lower.contains('.var.') || lower.contains('freebayes.raw')) && !lower.contains('.ann.')) {
            return null
        }
    }

    if (params.annotate_main_vcf) {
        if ((lower.endsWith('.vcf.gz') || lower.endsWith('.vcf.gz.tbi')) && 
            !lower.contains('.ann.') && 
            !lower.contains('.var.') && 
            !lower.contains('freebayes.raw')) { 
            
            // Keep the All Positions VCF (Backbone)
            if (lower.contains('all.pos')) {
                return name
            }
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
