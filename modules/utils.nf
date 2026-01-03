// modules/utils.nf

def cleanStr(v) {
    if (v == null) return null
    v.toString().replace('\r','').trim()
}

def nullish(v) {
    if (v == null) return true
    def s = cleanStr(v)
    if (!s) return true
    def sl = s.toLowerCase()
    return (s == "." || sl == "na" || sl == "n/a" || sl == "null")
}

def sanitizeId(v) {
    def s = cleanStr(v)
    if (!s) return s
    s.replaceAll(/[^A-Za-z0-9_.-]+/, "_")
}

def inferRunId(r1) {
    def name = new File(r1).getName()
    name = name.replaceAll(/(\.fastq|\.fq)(\.gz|\.bz2)?$/,'')
    name = name.replaceAll(/(\.gz|\.bz2)$/,'')
    sanitizeId(name)
}


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

    // B. ANNOTATION LOGIC
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
