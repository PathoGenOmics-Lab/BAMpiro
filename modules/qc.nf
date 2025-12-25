nextflow.enable.dsl=2

process VALIDATE_RAW_READS_PE {
    tag "Validate PE: ${sampleId}"
    cpus 1
    
    input:
    tuple val(sampleId), val(runId), path(r1), path(r2), val(refId), val(taxId)
    
    output:
    tuple val(sampleId), val(runId), path(r1), path(r2), val(refId), val(taxId), emit: reads
    
    shell:
    '''
    set -euo pipefail
    # Check if files are valid gzip or non-empty files
    if [[ "!{r1}" == *.gz ]]; then gzip -t "!{r1}" >/dev/null; else [[ -s "!{r1}" ]]; fi
    if [[ "!{r2}" == *.gz ]]; then gzip -t "!{r2}" >/dev/null; else [[ -s "!{r2}" ]]; fi
    '''
}

process VALIDATE_RAW_READS_SE {
    tag "Validate SE: ${sampleId}"
    cpus 1
    
    input:
    tuple val(sampleId), val(runId), path(r1), val(refId), val(taxId)
    
    output:
    tuple val(sampleId), val(runId), path(r1), val(refId), val(taxId), emit: reads
    
    shell:
    '''
    set -euo pipefail
    if [[ "!{r1}" == *.gz ]]; then gzip -t "!{r1}" >/dev/null; else [[ -s "!{r1}" ]]; fi
    '''
}
