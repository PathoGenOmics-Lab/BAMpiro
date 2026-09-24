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
    // Sized from a measured run rather than from caution. Three of these tasks averaged 102% CPU
    // over ten minutes on a 12-CPU reservation: kraken2 threads well but spends that time faulting
    // the database in, and the two steps after it were single-threaded until bgzip. 8 is what bgzip
    // actually scales to here, and a smaller reservation is scheduled sooner on a shared queue,
    // which is most of what a run waits for.
    cpus 8
    // Sized from params.kraken_memory, because the right number is a property of the DATABASE and
    // not of this pipeline: with --memory-mapping the resident set is the part of the index
    // actually probed, measured at 34.7-39.7 GB against a 133 GB standard database and well under
    // 20 GB against the capped 16 GB one. A fixed 56 GB here is what put 61 of these behind
    // QOSMaxMemoryPerUser on a real cohort. Scales with the attempt, like every other memory
    // directive in this pipeline, so an underestimate costs a retry and not the run.
    memory { (params.kraken_memory as nextflow.util.MemoryUnit) * task.attempt }
    
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

    # Compress outputs. bgzip, not gzip, and with the task's own threads: gzip is single-threaded,
    # so on a 1.5M-pair sample it spent 72 s compressing 480 MB per mate while the other eleven
    # reserved cores sat idle. bgzip -@ 8 does the same work in 2.8 s. BGZF is valid gzip, every
    # downstream reader is unaffected, and the output is byte-for-byte identical after decompression
    # and slightly smaller. Measured against this container, not assumed.
    bgzip -@ !{task.cpus} -f ${prefix}_R1.kraken.fq
    bgzip -@ !{task.cpus} -f ${prefix}_R2.kraken.fq
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
    cpus 8              // as KRAKEN_FILTER_PE above
    memory { (params.kraken_memory as nextflow.util.MemoryUnit) * task.attempt }
    
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

    # bgzip for the reason given in KRAKEN_FILTER_PE above.
    bgzip -@ !{task.cpus} -f ${prefix}.kraken.fq
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
    # A container run as the invoking user (the docker profile) has HOME=/, which it cannot write:
    # genmap then prints a mkdir error where its version goes.
    [ -w "${HOME:-/}" ] || export HOME="$PWD"

    ver() {
        # $1 = label, $2 = prefix to strip, rest = version command. The first line that carries a
        # version number, not the first line: bwa-mem2 opens with the SIMD build it launches.
        local raw line
        raw=$("${@:3}" 2>&1) || raw=""
        line=$(grep -m1 -E '[0-9]+\.[0-9]+' <<< "$raw" || true)
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
      # get_MNV runs in its own image until the pipeline image carries it: NA here then, and each
      # sample's mnv/<sample>.<ref>.mnv.manifest.json records the version that ran.
      ver get_mnv   'get_mnv '           get_mnv   --version
      echo "get_mnv container: !{params.mnv_container ?: params.container}"

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
