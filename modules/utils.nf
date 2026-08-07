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
 * A parameter read as a real boolean, whatever type it arrived as.
 *
 * `--make_qc_report false` on the command line does not reach the workflow as `false`. Up to
 * Nextflow 24 it was coerced to the type of the config default; from 26 it stays the STRING
 * "false", and a non-empty String is truthy in Groovy. So `if (params.make_qc_report)` was true
 * and the report was produced by a run that had explicitly asked for no report. Every boolean
 * flag in this pipeline had the same hole, and it opened silently on a Nextflow upgrade rather
 * than on any change here: the run does the opposite of what was asked and says nothing.
 *
 * Anything that is not recognisably true is false, so a typo turns a feature off rather than
 * quietly on.
 */
def asBool(v) {
    if (v instanceof Boolean) return v
    if (v == null) return false
    return cleanStr(v).toLowerCase() in ['true', 'yes', 'on', '1']
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
    if (asBool(params.nested_output)) {
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
    if (asBool(params.output_cram) && (lower.endsWith('.bam') || lower.endsWith('.bam.bai'))) {
        return null
    }

    // A2. Pre-filter dedup alignment: once the length-aware filter runs, filtered.* is the analysis
    //     alignment and final.* is a ~redundant second full copy. Drop it from the outdir by
    //     default (still kept in the work dir, so the virgin/raw branch is unaffected). When the
    //     filter is off, final.* is the ONLY alignment and this gate does not trigger. Matches both
    //     .final.bam(.bai) and .final.cram(.crai).
    if ((lower.contains('.final.bam') || lower.contains('.final.cram'))
        && asBool(params.dynamic_read_filter) && !asBool(params.publish_prefilter_bam)) {
        return null
    }

    // A3. All-positions (per-position) VCFs are the consensus substrate. Gate publishing here
    //     (hoisted OUT of the annotate block so the flags work regardless of annotate_main_vcf).
    if (lower.contains('all.pos') && (lower.endsWith('.vcf.gz') || lower.endsWith('.vcf.gz.tbi'))) {
        if (!asBool(params.publish_allpos_vcf)) return null
        if (lower.contains('.raw.') && !asBool(params.publish_virgin_allpos_vcf)) return null
        return name
    }

    // B. ANNOTATION LOGIC (Filter Raw VCFs if annotation is enabled)
    if (asBool(params.annotate_legacy_vcfs)) {
        if ((lower.contains('.var.') || lower.contains('freebayes.raw')) && !lower.contains('.ann.')) {
            return null
        }
    }

    if (asBool(params.annotate_main_vcf)) {
        // all.pos is already handled above; here we only drop the un-annotated main/virgin VCFs
        // (their .ann. versions are the deliverables).
        if ((lower.endsWith('.vcf.gz') || lower.endsWith('.vcf.gz.tbi')) &&
            !lower.contains('.ann.') &&
            !lower.contains('.var.') &&
            !lower.contains('freebayes.raw')) {
            return null
        }
    }

    // C. LINEAGE FOLDER (PATHOTYPR) - must precede the generic stats/.tsv rule below, otherwise the
    //    pathotypr *.tsv deliverables end in .tsv and get routed to stats/ (this block never fires).
    if (lower.contains('pathotypr') || lower.endsWith('.lineage.tsv') || lower.endsWith('.split_kmer.tsv')) {
        return "lineage/${name}"
    }

    // D. STATS FOLDER
    if (lower.endsWith('.log')    || lower.endsWith('.stats') ||
        lower.endsWith('.json')   || lower.endsWith('.html')  ||
        lower.endsWith('.report') || lower.endsWith('.csv')   ||
        lower.endsWith('.tsv')    || lower.endsWith('.txt')) {
        return "stats/${name}"
    }

    // E. DEFAULT PUBLISH
    return name 
}

/**
 * Parse the declared parameters straight out of nextflow.config.
 *
 * Reading the config rather than keeping a second list of names means the two can never drift.
 * Returns a list of [name, default, description, section] in declaration order.
 */
def declaredParams(String configPath) {
    def f = new File(configPath)
    if (!f.exists()) return []

    def out = []
    def section = 'General'
    def depth = 0
    def inside = false

    f.readLines().each { raw ->
        def line = raw.trim()
        if (!inside) {
            if (line ==~ /^params\s*\{.*/) { inside = true; depth = 1 }
            return
        }
        // Track nesting by line so the matching close ends the block and a later one is ignored.
        // A balanced "${...}" inside a value nets to zero, which is what we want.
        def net = depth + line.count('{') - line.count('}')
        if (net <= 0) { inside = false; return }
        depth = net

        def head = (line =~ /^\/\/\s*---\s*(.+?)\s*-*\s*$/)
        if (head) { section = head[0][1]; return }
        if (line.startsWith('//') || !line) return
        // A quoted value is matched whole first, so a "//" inside it (a URL, a docker:// ref) is
        // not mistaken for the start of a trailing comment.
        def m = (line =~ /^([A-Za-z_][A-Za-z0-9_]*)\s*=\s*("[^"]*"|'[^']*'|.+?)\s*(?:\/\/\s*(.*))?$/)
        if (m) out << [name: m[0][1], value: m[0][2], description: (m[0][3] ?: '').trim(), section: section]
    }
    return out
}

/**
 * Reject an unknown --parameter instead of accepting it and doing nothing.
 *
 * Nextflow puts every --flag into params whether or not the pipeline declares it, so a typo like
 * `--treads 16` is silently ignored and the run completes with the default. Anything Nextflow
 * itself injects is allowed through.
 */
def validateParams(Map params, String configPath) {
    def declared = declaredParams(configPath)*.name as Set
    if (!declared) return   // config unreadable: do not block the run over a validation helper

    // Nextflow-injected and convention keys that are never declared in the params block.
    def allowed = ['help', 'validate_params', 'monochrome_logs', 'igenomes_base'] as Set

    def unknown = params.keySet().findAll { !declared.contains(it) && !allowed.contains(it) }.sort()
    if (!unknown) return

    def msg = new StringBuilder("Unrecognised parameter(s): ${unknown.collect{ "--${it}" }.join(', ')}\n")
    unknown.each { u ->
        // Cheap edit-distance suggestion: most real cases are one transposed or dropped letter.
        def near = declared.findAll { d ->
            d.size() > 3 && (d.toLowerCase().contains(u.toLowerCase()) || u.toLowerCase().contains(d.toLowerCase())
                             || (d.size() - u.size()).abs() <= 2 && (d.toSet().intersect(u.toSet()).size() >= u.size() - 1))
        }.sort().take(3)
        if (near) msg << "  --${u}: did you mean ${near.collect{ "--${it}" }.join(' or ')}?\n"
    }
    msg << "Run with --help to list every parameter."
    throw new RuntimeException(msg.toString())
}

/**
 * Print every declared parameter with its default, grouped by the section headings in the config.
 */
def paramsHelp(String configPath, String version) {
    def declared = declaredParams(configPath)
    def out = new StringBuilder("\nBAMpiro ${version}\n\n")
    out << "  nextflow run main.nf --tsv samples.tsv --outdir results -profile local,docker\n\n"
    out << "Profiles: standard (default, this host), local, slurm, garnatxa, docker, generic, test\n"
    declared.groupBy{ it.section }.each { section, entries ->
        out << "\n${section}\n"
        entries.each { p ->
            def val = p.value.replaceAll(/\s+/, ' ')
            if (val.size() > 46) val = val.take(43) + '...'
            out << String.format("  --%-26s %s\n", p.name, val)
        }
    }
    return out.toString()
}
