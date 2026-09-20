nextflow.enable.dsl=2

// Import centralized functions for path generation and file classification
include { getSavePath; getSampleDir } from './utils'

/* ====================================================================
    QC MODULES
    Contains: Validation, Kraken2, FastP, and MultiQC processes
====================================================================
*/

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

    stub:
    """
    # the outputs re-emit the staged input reads, so there is nothing to create
    true
    """
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

    stub:
    """
    true
    """
}

process KRAKEN_FILTER_PE {
    tag "Kraken PE: ${sampleId}"
    cpus 12
    memory '80 GB'
    
    // Use getSampleDir for nested output support
    publishDir path: { "${params.outdir}/${getSampleDir(sampleId, params)}" }, mode: params.publish_mode, saveAs: { filename -> getSavePath(filename, params) }
    
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

    # KrakenTools helper: reuse the vendored copy in bin/ if present, else fetch it (no repeat
    # download per task once vendored, and no network dependency).
    if [ -f "!{projectDir}/bin/extract_kraken_reads.py" ]; then
      cp "!{projectDir}/bin/extract_kraken_reads.py" extract_kraken_reads.py
    else
      wget -qO extract_kraken_reads.py "!{params.krakentools_url}"
    fi
    chmod +x extract_kraken_reads.py

    # --memory-mapping (optional): share the DB in RAM across parallel tasks instead of each
    # loading the whole DB -> much lower peak memory.
    MM=""
    if [[ "!{params.kraken_memory_mapping}" == "true" ]]; then MM="--memory-mapping"; fi

    # Run Kraken2
    kraken2 --db !{kraken_db} --threads !{task.cpus} $MM --paired !{r1} !{r2} \
      --output kraken.output --report ${prefix}.kraken.report

    # Filter reads based on TaxID
    python3 extract_kraken_reads.py -k kraken.output -r ${prefix}.kraken.report \
      -s !{r1} -s2 !{r2} -t !{taxId} --include-children --fastq-output \
      -o ${prefix}_R1.kraken.fq -o2 ${prefix}_R2.kraken.fq

    # Compress outputs
    gzip -f ${prefix}_R1.kraken.fq
    gzip -f ${prefix}_R2.kraken.fq
    '''

    stub:
    """
    touch ${sampleId}__${runId}_R1.kraken.fq.gz
    touch ${sampleId}__${runId}_R2.kraken.fq.gz
    touch ${sampleId}__${runId}.kraken.report
    """
}

process KRAKEN_FILTER_SE {
    tag "Kraken SE: ${sampleId}"
    cpus 12
    memory '80 GB'
    
    // Use getSampleDir for nested output support
    publishDir path: { "${params.outdir}/${getSampleDir(sampleId, params)}" }, mode: params.publish_mode, saveAs: { filename -> getSavePath(filename, params) }
    
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

    # KrakenTools helper: reuse the vendored copy in bin/ if present, else fetch it.
    if [ -f "!{projectDir}/bin/extract_kraken_reads.py" ]; then
      cp "!{projectDir}/bin/extract_kraken_reads.py" extract_kraken_reads.py
    else
      wget -qO extract_kraken_reads.py "!{params.krakentools_url}"
    fi
    chmod +x extract_kraken_reads.py

    MM=""
    if [[ "!{params.kraken_memory_mapping}" == "true" ]]; then MM="--memory-mapping"; fi

    kraken2 --db !{kraken_db} --threads !{task.cpus} $MM !{r1} \
      --output kraken.output --report ${prefix}.kraken.report

    python3 extract_kraken_reads.py -k kraken.output -r ${prefix}.kraken.report \
      -s !{r1} -t !{taxId} --include-children --fastq-output \
      -o ${prefix}.kraken.fq

    gzip -f ${prefix}.kraken.fq
    '''

    stub:
    """
    touch ${sampleId}__${runId}.kraken.fq.gz
    touch ${sampleId}__${runId}.kraken.report
    """
}

process FASTP_PE {
    tag "fastp PE: ${sampleId}"
    cpus 4
    memory { 4.GB * task.attempt }
    
    // Use getSampleDir for nested output support
    publishDir path: { "${params.outdir}/${getSampleDir(sampleId, params)}" }, mode: params.publish_mode, saveAs: { filename -> getSavePath(filename, params) }
    
    input:
    tuple val(sampleId), val(runId), path(r1), path(r2), val(refId), val(taxId)
    
    output:
    tuple val(sampleId), val(runId), path("${sampleId}__${runId}_R1.clean.fq.gz"), path("${sampleId}__${runId}_R2.clean.fq.gz"), val(refId), val(taxId), emit: pe_reads
    tuple val(sampleId), val(runId), path("${sampleId}__${runId}_se_combined.fq.gz"), val(refId), val(taxId), emit: se_reads
    tuple val(sampleId), path("${sampleId}__${runId}_fastp.json"), emit: json
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

    stub:
    """
    touch ${sampleId}__${runId}_R1.clean.fq.gz
    touch ${sampleId}__${runId}_R2.clean.fq.gz
    touch ${sampleId}__${runId}_se_combined.fq.gz
    touch ${sampleId}__${runId}_fastp.json
    touch ${sampleId}__${runId}_fastp.html
    """
}

process FASTP_SE {
    tag "fastp SE: ${sampleId}"
    cpus 4
    memory { 4.GB * task.attempt }
    
    // Use getSampleDir for nested output support
    publishDir path: { "${params.outdir}/${getSampleDir(sampleId, params)}" }, mode: params.publish_mode, saveAs: { filename -> getSavePath(filename, params) }
    
    input:
    tuple val(sampleId), val(runId), path(r1), val(refId), val(taxId)
    
    output:
    tuple val(sampleId), val(runId), path("${sampleId}__${runId}_SE.clean.fq.gz"), val(refId), val(taxId), emit: se_reads
    tuple val(sampleId), path("${sampleId}__${runId}_fastp.json"), emit: json
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

    stub:
    """
    touch ${sampleId}__${runId}_SE.clean.fq.gz
    touch ${sampleId}__${runId}_fastp.json
    touch ${sampleId}__${runId}_fastp.html
    """
}

process MULTIQC {
    tag "MultiQC"
    // MultiQC has its own dedicated path and does not use per-sample logic
    publishDir "${params.outdir}/multiqc", mode: params.publish_mode
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
    title: "BAMpiro Report"
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

    stub:
    """
    touch ${report_name}.html
    mkdir -p ${report_name}_data
    """
}

process DUMP_VERSIONS {
    tag "versions"
    cpus 1
    memory '1 GB'
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_mode

    output:
    path "software_versions_mqc.yml", emit: mqc
    path "software_versions.txt",     emit: txt

    shell:
    '''
    # Provenance: record the exact tool versions from the (pinned) container plus the
    # pipeline/Nextflow versions, so every run documents its own software environment.
    ver() {
        # $1 = label, $2 = prefix to strip, rest = version command
        local raw line
        raw=$("${@:3}" 2>&1) || raw=""
        line=$(head -n1 <<< "$raw")
        if [ -n "$2" ]; then line=${line#"$2"}; fi
        line=${line#"${line%%[![:space:]]*}"}   # trim leading whitespace
        echo "$1: ${line:-NA}"
    }

    {
      echo "BAMpiro: !{workflow.manifest.version}"
      echo "Nextflow: !{workflow.nextflow.version}"
      echo "container: !{params.container}"
      ver samtools  'samtools '          samtools  --version
      ver bcftools  'bcftools '          bcftools  --version
      ver freebayes 'version: '          freebayes --version
      ver bwa-mem2  ''                   bwa-mem2  version
      ver fastp     'fastp '             fastp     --version
      ver kraken2   'Kraken version '    kraken2   --version
      ver genmap    ''                   genmap    --version
      ver snpEff    ''                   snpEff    -version
      ver python    'Python '            python3   --version
      ver multiqc   'multiqc, version '  multiqc   --version
      ver pathotypr 'pathotypr '         pathotypr --version

      # The marker catalogue is provenance too, and it is the half that cannot be recovered from
      # the results. Catalogue v1.0.0 assigned each variant a single drug inherited from its gene;
      # v1.0.2 grades per variant-drug pair and disagrees with it on 15,969 rows, so two runs with
      # identical calls can mean different things. The checksum covers the case of a catalogue
      # supplied with --pathotypr_dr_markers, which carries no version string at all.
      if [ -s "!{params.pathotypr_dr_markers}" ]; then
          echo "pathotypr markers: $(cat "$(dirname "!{params.pathotypr_dr_markers}")/VERSION" 2>/dev/null || echo 'not recorded')"
          echo "pathotypr dr_markers sha256: $(sha256sum "!{params.pathotypr_dr_markers}" | cut -d' ' -f1)"
      else
          echo "pathotypr markers: NA"
      fi
    } > software_versions.txt

    # MultiQC custom-content section (files ending in _mqc.yml are auto-detected)
    {
      echo 'id: "software_versions"'
      echo 'section_name: "Software Versions"'
      echo 'section_href: "https://github.com/PathoGenOmics-Lab/BAMpiro"'
      echo 'plot_type: "html"'
      echo 'description: "Captured at runtime from the pipeline container."'
      echo 'data: |'
      echo '    <dl class="dl-horizontal">'
      while IFS= read -r kv; do
          k=${kv%%:*}
          v=${kv#*: }
          echo "        <dt>${k}</dt><dd><samp>${v}</samp></dd>"
      done < software_versions.txt
      echo '    </dl>'
    } > software_versions_mqc.yml
    '''

    stub:
    """
    touch software_versions.txt
    touch software_versions_mqc.yml
    """
}
