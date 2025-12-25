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

process KRAKEN_FILTER_PE {
    tag "Kraken PE: ${sampleId}"
    cpus 12
    memory '80 GB'
    publishDir "${params.outdir}/${sampleId}", mode: 'copy'
    
    input:
    tuple val(sampleId), val(runId), path(r1), path(r2), val(refId), val(taxId)
    val kraken_db
    
    output:
    tuple val(sampleId), val(runId), path("${sampleId}__${runId}_R1.kraken.fq.gz"), path("${sampleId}__${runId}_R2.kraken.fq.gz"), val(refId), val(taxId), emit: reads
    path("${sampleId}__${runId}.kraken.report"), emit: report
    
    shell:
    '''
    set -euo pipefail
    prefix="!{sampleId}__!{runId}"
    
    # Download the KrakenTools script dynamically
    wget -qO extract_kraken_reads.py "!{params.krakentools_url}"
    chmod +x extract_kraken_reads.py

    # Run Kraken2
    kraken2 --db !{kraken_db} --threads !{task.cpus} --paired !{r1} !{r2} \
      --output kraken.output --report ${prefix}.kraken.report

    # Filter reads based on TaxID
    python3 extract_kraken_reads.py -k kraken.output -r ${prefix}.kraken.report \
      -s !{r1} -s2 !{r2} -t !{taxId} --include-children --fastq-output \
      -o ${prefix}_R1.kraken.fq -o2 ${prefix}_R2.kraken.fq

    # Compress outputs
    gzip -f ${prefix}_R1.kraken.fq
    gzip -f ${prefix}_R2.kraken.fq
    '''
}

process KRAKEN_FILTER_SE {
    tag "Kraken SE: ${sampleId}"
    cpus 12
    memory '80 GB'
    publishDir "${params.outdir}/${sampleId}", mode: 'copy'
    
    input:
    tuple val(sampleId), val(runId), path(r1), val(refId), val(taxId)
    val kraken_db
    
    output:
    tuple val(sampleId), val(runId), path("${sampleId}__${runId}.kraken.fq.gz"), val(refId), val(taxId), emit: reads
    path("${sampleId}__${runId}.kraken.report"), emit: report
    
    shell:
    '''
    set -euo pipefail
    prefix="!{sampleId}__!{runId}"
    
    wget -qO extract_kraken_reads.py "!{params.krakentools_url}"
    chmod +x extract_kraken_reads.py

    kraken2 --db !{kraken_db} --threads !{task.cpus} !{r1} \
      --output kraken.output --report ${prefix}.kraken.report

    python3 extract_kraken_reads.py -k kraken.output -r ${prefix}.kraken.report \
      -s !{r1} -t !{taxId} --include-children --fastq-output \
      -o ${prefix}.kraken.fq

    gzip -f ${prefix}.kraken.fq
    '''
}
