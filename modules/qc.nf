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

process FASTP_PE {
    tag "fastp PE: ${sampleId}"
    cpus 4
    memory '8 GB'
    publishDir "${params.outdir}/${sampleId}", mode: 'copy', saveAs: { filename ->
        if (filename.endsWith('.json') || filename.endsWith('.html')) return "stats/${filename}"
        return filename
    }
    
    input:
    tuple val(sampleId), val(runId), path(r1), path(r2), val(refId), val(taxId)
    
    output:
    tuple val(sampleId), val(runId), path("${sampleId}__${runId}_R1.clean.fq.gz"), path("${sampleId}__${runId}_R2.clean.fq.gz"), val(refId), val(taxId), emit: pe_reads
    tuple val(sampleId), val(runId), path("${sampleId}__${runId}_se_combined.fq.gz"), val(refId), val(taxId), emit: se_reads
    path("${sampleId}__${runId}_fastp.json"), emit: json
    path("${sampleId}__${runId}_fastp.html"), emit: html
    
    shell:
    '''
    set -euo pipefail
    prefix="!{sampleId}__!{runId}"
    
    # Run FastP
    fastp -i !{r1} -I !{r2} \
      --out1 ${prefix}_R1.clean.fq.gz --out2 ${prefix}_R2.clean.fq.gz \
      --merge --merged_out ${prefix}_merged.fq.gz \
      --unpaired1 ${prefix}_u1.fq.gz --unpaired2 ${prefix}_u2.fq.gz \
      --detect_adapter_for_pe \
      --thread !{task.cpus} \
      --length_required !{params.fastp_min_length} \
      --json ${prefix}_fastp.json --html ${prefix}_fastp.html

    # Handle optional merged/unpaired outputs
    for f in ${prefix}_merged.fq.gz ${prefix}_u1.fq.gz ${prefix}_u2.fq.gz; do
      if [[ ! -f "$f" ]]; then gzip -c /dev/null > "$f"; fi
    done
    
    # Concatenate orphans and merged reads into a single SE file
    cat ${prefix}_merged.fq.gz ${prefix}_u1.fq.gz ${prefix}_u2.fq.gz > ${prefix}_se_combined.fq.gz
    rm -f ${prefix}_merged.fq.gz ${prefix}_u1.fq.gz ${prefix}_u2.fq.gz
    '''
}

process FASTP_SE {
    tag "fastp SE: ${sampleId}"
    cpus 4
    memory '8 GB'
    publishDir "${params.outdir}/${sampleId}", mode: 'copy', saveAs: { filename ->
        if (filename.endsWith('.json') || filename.endsWith('.html')) return "stats/${filename}"
        return filename
    }
    
    input:
    tuple val(sampleId), val(runId), path(r1), val(refId), val(taxId)
    
    output:
    tuple val(sampleId), val(runId), path("${sampleId}__${runId}_SE.clean.fq.gz"), val(refId), val(taxId), emit: se_reads
    path("${sampleId}__${runId}_fastp.json"), emit: json
    path("${sampleId}__${runId}_fastp.html"), emit: html
    
    shell:
    '''
    set -euo pipefail
    prefix="!{sampleId}__!{runId}"
    
    fastp -i !{r1} \
      --out1 ${prefix}_SE.clean.fq.gz \
      --thread !{task.cpus} \
      --length_required !{params.fastp_min_length} \
      --json ${prefix}_fastp.json --html ${prefix}_fastp.html
    '''
}

process MULTIQC {
    tag "MultiQC"
    publishDir "${params.outdir}/multiqc", mode: 'copy'
    cpus 2
    memory '4 GB'

    input:
    path multiqc_files
    val report_name

    output:
    path "${report_name}.html", emit: report
    path "${report_name}_data", emit: data

    script:
    """
    # Create configuration file to ORDER the report logically.
    # We do NOT hide 'generalstats' to keep the summary table.
    cat <<EOF > multiqc_config.yaml
    title: "BAMpiro Report 🧛‍♂️"
    module_order:
        - fastp
        - samtools
        - kraken
        - snpeff
    EOF

    # Run MultiQC
    # -c : Use custom config
    # -n : Set dynamic output name
    multiqc . -c multiqc_config.yaml -n ${report_name}
    """
}
