nextflow.enable.dsl=2

// Import centralized functions for path generation and file classification
include { getSavePath; getSampleDir } from './utils'

/* ====================================================================
    VARIANTS MODULES
    Contains: FreeBayes Calling, Backbone Generation, and VCF Merging
==================================================================== */

process CALL_FREEBAYES {
    tag "FreeBayes: ${sampleId}"
    
    // Use getSampleDir for nested output support.
    // getSavePath handles the filtering of intermediate files if needed.
    publishDir "${params.outdir}/${getSampleDir(sampleId, params)}", mode: params.publish_mode, saveAs: { filename -> getSavePath(filename, params) }
    
    // FreeBayes is single-threaded unless region-parallel is enabled; only bgzip uses extra cores.
    cpus { params.freebayes_parallel ? (params.freebayes_parallel_jobs as int) : 2 }
    memory { 8.GB * task.attempt }

    input:
    tuple val(sampleId), val(refId), path(bam), path(bai), path(ref_fa), path(exclude_txt)

    output:
    tuple val(sampleId), val(refId), path("valid_snps_formatted.vcf.gz"), path("valid_snps_formatted.vcf.gz.tbi"), emit: snps
    tuple val(sampleId), val(refId), path("${sampleId}.${refId}.mask_sites.tsv"), emit: mask_sites
    tuple val(sampleId), val(refId), path("${sampleId}.${refId}.freebayes.raw.vcf.gz"), path("${sampleId}.${refId}.freebayes.raw.vcf.gz.tbi"), emit: raw_fb
    // Splits for specific annotation
    tuple val(sampleId), val(refId), path("${sampleId}.${refId}.var.homo.SNPs.vcf"),  emit: homo_snp
    tuple val(sampleId), val(refId), path("${sampleId}.${refId}.var.het.SNPs.vcf"),   emit: het_snp
    tuple val(sampleId), val(refId), path("${sampleId}.${refId}.var.homo.indel.vcf"), emit: homo_indel

    shell:
    """
    set -euo pipefail
    
    # Set up a temporary directory
    if [ -n "\${SLURM_TMPDIR:-}" ]; then export TMPDIR="\$SLURM_TMPDIR"; else export TMPDIR="./tmp_fb"; fi
    mkdir -p "\$TMPDIR"
    
    # Helper function for safe tabix indexing
    safe_tabix () {
      local gz="\$1"; local idx="\${gz}.tbi"
      set +e; tabix -f -p vcf "\$gz"; st=\$?; set -e
      if [ \$st -ne 0 ]; then
        # Only a genuinely header-only VCF gets an empty index. bcftools reads the whole file,
        # so a truncated/malformed VCF errors under set -e (fails loud) instead of a fake index.
        local nrec; nrec=\$(bcftools view -H "\$gz" | wc -l)
        if [ "\$nrec" -ne 0 ]; then
          echo "ERROR: tabix failed on \$gz with \$nrec records (malformed/truncated) -> refusing a fake index" >&2
          exit \$st
        fi
        : > "\$idx"
      fi
    }

    # 1. Handle Excluded Regions (Repeats)
    # Keep only numeric-coordinate rows -> strips EVERY header (the combined nucmer+genmap exclude has two)
    awk 'tolower(\$1)!="chrom" && \$2 ~ /^[0-9]+\$/ && \$3 ~ /^[0-9]+\$/ {print \$1"\\t"\$2"\\t"\$3}' !{exclude_txt} > exclude.regions.tsv || true
    EXCL_ARG=""
    if [[ "!{params.exclude_repeats}" == "true" && -s exclude.regions.tsv ]]; then 
        EXCL_ARG="-T ^exclude.regions.tsv"
    fi

    # 2. Run FreeBayes (optionally region-parallel; params.freebayes_parallel)
    if [[ "!{params.freebayes_parallel}" == "true" ]]; then
      samtools faidx !{ref_fa}
      fasta_generate_regions.py !{ref_fa}.fai !{params.freebayes_parallel_chunk} > fb_regions.txt
      freebayes-parallel fb_regions.txt !{task.cpus} \\
        -f !{ref_fa} -p !{params.freebayes_ploidy} \\
        --min-alternate-count !{params.freebayes_min_alt_count} \\
        --min-alternate-fraction !{params.freebayes_min_alt_fraction} \\
        --min-mapping-quality !{params.freebayes_min_map_qual} \\
        --min-base-quality !{params.freebayes_min_base_qual} \\
        !{bam} > raw_freebayes.vcf
    else
      freebayes -f !{ref_fa} -p !{params.freebayes_ploidy} \\
        --min-alternate-count !{params.freebayes_min_alt_count} \\
        --min-alternate-fraction !{params.freebayes_min_alt_fraction} \\
        --min-mapping-quality !{params.freebayes_min_map_qual} \\
        --min-base-quality !{params.freebayes_min_base_qual} \\
        !{bam} > raw_freebayes.vcf
    fi

    # Compress Raw Output
    bgzip -@ !{task.cpus} -c raw_freebayes.vcf > !{sampleId}.!{refId}.freebayes.raw.vcf.gz
    safe_tabix !{sampleId}.!{refId}.freebayes.raw.vcf.gz

    # 3. Normalize Variants
    bcftools norm -m - -a -f !{ref_fa} raw_freebayes.vcf | \
    bcftools view -e 'GT="0/0" || GT="0"' -Ov -o normalized.vcf
    
    # 4. Generate Mask Sites (Positions with low allelic balance support)
    ( bcftools view \$EXCL_ARG --types snps,mnps \\
        -i "(INFO/SAF[0] < !{params.min_alt_fwd} || INFO/SAR[0] < !{params.min_alt_rev})" \\
        normalized.vcf \\
      | bcftools query -f '%CHROM\\t%POS\\n' \\
      | tr -d '\\r' | awk 'NF' | sort -u > !{sampleId}.!{refId}.mask_sites.tsv ) || true
    
    # Ensure file exists even if empty
    [ -s !{sampleId}.!{refId}.mask_sites.tsv ] || : > !{sampleId}.!{refId}.mask_sites.tsv

    # 5. Define Filter Rules (Bash variables for bcftools expressions)
    base_filt="(FMT/DP >= !{params.filter_min_dp} || INFO/DP >= !{params.filter_min_dp}) && INFO/SAF[0] >= !{params.min_alt_fwd} && INFO/SAR[0] >= !{params.min_alt_rev}"
    af_guard="(INFO/AO[0] + INFO/RO) > 0"
    af_calc="((INFO/AO[0] * 1.0) / (INFO/AO[0] + INFO/RO))"
    
    # Homozygous Rule: High Allele Fraction
    hom_rule="\$base_filt && \$af_guard && \$af_calc >= !{params.hom_threshold}"
    # Heterozygous Rule: Mid Allele Fraction
    het_rule="\$base_filt && \$af_guard && \$af_calc >= !{params.het_min_frac} && \$af_calc < !{params.hom_threshold}"

    # 6. Apply Filters
    bcftools view \$EXCL_ARG --types snps -i "\$hom_rule || \$het_rule" normalized.vcf > valid_snps.vcf || true
    [ -s valid_snps.vcf ] || bcftools view -h normalized.vcf > valid_snps.vcf

    # 7. Output Split Files (Legacy format required)
    bcftools view \$EXCL_ARG --types snps,mnps -i "\$hom_rule" -Ov -o !{sampleId}.!{refId}.var.homo.SNPs.vcf  normalized.vcf
    bcftools view \$EXCL_ARG --types snps,mnps -i "\$het_rule" -Ov -o !{sampleId}.!{refId}.var.het.SNPs.vcf   normalized.vcf
    bcftools view \$EXCL_ARG --types indels    -i "\$hom_rule" -Ov -o !{sampleId}.!{refId}.var.homo.indel.vcf normalized.vcf

    # 8. Format for Backbone Integration
    # Adds ADP, WT, HET, HOM, NC tags to INFO/FORMAT for the consensus step. INFO is tokenised on ';'
    # (anchored ^DP=/^RO=/^AO=) so DP is never grabbed from a look-alike field, and the genotype is
    # re-validated with RO/AO.
    awk -v HETMIN=!{params.het_min_frac} 'BEGIN{OFS="\\t"}
      /^##/ { print; next }
      /^#CHROM/ { print; next }
      {
        dp=0; ro=0; ao=0
        ni=split(\$8,I,";")
        for(k=1;k<=ni;k++){
          if(I[k] ~ /^DP=/){ v=I[k]; sub(/^DP=/,"",v); dp=v+0 }
          else if(I[k] ~ /^RO=/){ v=I[k]; sub(/^RO=/,"",v); ro=v+0 }
          else if(I[k] ~ /^AO=/){ v=I[k]; sub(/^AO=/,"",v); split(v,ax,","); ao=ax[1]+0 }
        }
        if (dp==0) { n=split(\$9,fmt,":"); m=split(\$10,dat,":"); for(i=1;i<=n;i++) if(fmt[i]=="DP") dp=dat[i] }
        if (dp==0) dp=ro+ao
        if (dp==0) dp=1
        split(\$10,b,":"); gtype=b[1]
        # Re-validate het: 'bcftools norm -m -' splits a multiallelic 1/2 into biallelic records whose
        # reference index 0 is an ARTIFACT (RO~0). A het whose reference allele is essentially
        # unsupported is not a real het -> homozygous-alt (matches the af-based hom/het rule and stops
        # the reference base leaking into the consensus IUPAC).
        if ((gtype=="0/1" || gtype=="1/0") && (ro+ao)>0 && ro/(ro+ao) < HETMIN) gtype="1/1"
        het=0; hom=0
        if (gtype=="0/1" || gtype=="1/0") het=1
        else if (gtype=="1/1") hom=1
        \$8="ADP="dp";WT=0;HET="het";HOM="hom";NC=0"
        \$9="GT:DP"
        \$10=gtype":"dp
        print
      }' valid_snps.vcf | bgzip -@ !{task.cpus} -c > valid_snps_formatted.vcf.gz
    
    safe_tabix valid_snps_formatted.vcf.gz
    """

    // Stubs create only the NEW outputs; note valid_snps_formatted.* is a fixed, non-interpolated name.
    stub:
    """
    touch valid_snps_formatted.vcf.gz
    touch valid_snps_formatted.vcf.gz.tbi
    touch ${sampleId}.${refId}.mask_sites.tsv
    touch ${sampleId}.${refId}.freebayes.raw.vcf.gz
    touch ${sampleId}.${refId}.freebayes.raw.vcf.gz.tbi
    touch ${sampleId}.${refId}.var.homo.SNPs.vcf
    touch ${sampleId}.${refId}.var.het.SNPs.vcf
    touch ${sampleId}.${refId}.var.homo.indel.vcf
    """
}

process CALL_BACKBONE {
    tag "Backbone: ${sampleId}"
    // samtools mpileup is single-threaded; only bgzip uses extra cores.
    cpus 2
    memory { 6.GB * task.attempt }
    
    input:
    tuple val(sampleId), val(refId), path(bam), path(bai), path(ref_fa)
    
    output:
    tuple val(sampleId), val(refId), path("backbone.vcf.gz"), path("backbone.vcf.gz.tbi"), emit: backbone
    tuple val(sampleId), val(refId), path("header_template.txt"), emit: header
    
    shell:
    """
    set -euo pipefail
    
    if [ -n "\${SLURM_TMPDIR:-}" ]; then export TMPDIR="\$SLURM_TMPDIR"; else export TMPDIR="./tmp_bb"; fi
    mkdir -p "\$TMPDIR"
    
    safe_tabix () {
      local gz="\$1"; local idx="\${gz}.tbi"
      set +e; tabix -f -p vcf "\$gz"; st=\$?; set -e
      if [ \$st -ne 0 ]; then
        local nrec; nrec=\$(bcftools view -H "\$gz" | wc -l)
        if [ "\$nrec" -ne 0 ]; then
          echo "ERROR: tabix failed on \$gz with \$nrec records (malformed/truncated) -> refusing a fake index" >&2
          exit \$st
        fi
        : > "\$idx"
      fi
    }
    
    # Ensure index exists
    samtools faidx !{ref_fa} >/dev/null 2>&1 || true

    # 1. Create VCF Header Template
    echo "##fileformat=VCFv4.2" > header_template.txt
    awk 'BEGIN{OFS=""} {print "##contig=<ID=" \$1 ",length=" \$2 ">"}' !{ref_fa}.fai >> header_template.txt
    echo '##INFO=<ID=ADP,Number=1,Type=Integer,Description="AverageDepth">' >> header_template.txt
    echo '##INFO=<ID=WT,Number=1,Type=Integer,Description="Legacy_WT">' >> header_template.txt
    echo '##INFO=<ID=HET,Number=1,Type=Integer,Description="Legacy_HET">' >> header_template.txt
    echo '##INFO=<ID=HOM,Number=1,Type=Integer,Description="Legacy_HOM">' >> header_template.txt
    echo '##INFO=<ID=NC,Number=1,Type=Integer,Description="Legacy_NC">' >> header_template.txt
    echo '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' >> header_template.txt
    echo '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">' >> header_template.txt
    echo '##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths ref,alt (backbone: pileup-derived; variants: freebayes RO,AO)">' >> header_template.txt
    echo '##FILTER=<ID=baq_dropout,Description="Well-covered site whose depth collapses below consensus_min_dp only under BAQ (indel-adjacent homopolymer); masked in consensus">' >> header_template.txt
    echo -e "#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\tFORMAT\\t!{sampleId}" >> header_template.txt

    # 2. Run Mpileup for ALL positions (-aa)
    # ADP/DP is the real (no-BAQ) depth. A site is flagged 'baq_dropout' when it is
    # well covered but BAQ collapses it to <= consensus_min_dp (indel-adjacent
    # homopolymer artifact), so the consensus masks it (X) instead of leaving a gap.
    if [[ "!{params.mask_baq_dropouts}" == "true" ]]; then
        # Materialise both pileups to files (NOT process substitution): a mpileup that dies mid-stream
        # inside <(...) is invisible to 'set -euo pipefail', and paste would then pad the truncated
        # BAQ column with empties -> dp_baq=0 -> every tail position falsely flagged 'baq_dropout'.
        # Writing to files makes the failure fatal, and we assert equal line counts before pasting.
        samtools mpileup -aa -B -f !{ref_fa} -Q !{params.allpos_min_bq} -d !{params.allpos_max_depth} !{bam} > mpileup_nobaq.txt
        samtools mpileup -aa    -f !{ref_fa} -Q !{params.allpos_min_bq} -d !{params.allpos_max_depth} !{bam} > mpileup_baq.txt
        n1=\$(wc -l < mpileup_nobaq.txt); n2=\$(wc -l < mpileup_baq.txt)
        if [ "\$n1" -ne "\$n2" ]; then echo "ERROR: mpileup line counts differ (\$n1 vs \$n2) -- truncated stream" >&2; exit 1; fi
        paste mpileup_nobaq.txt mpileup_baq.txt \\
        | awk -v MINCOV=!{params.allpos_min_cov} -v GAPDP=!{params.consensus_min_dp} -v OFS="\\t" '
            function count_bases(s,   i,n,c,num){
              RO_G=0; AO_G=0; i=1; n=length(s)
              while(i<=n){
                c=substr(s,i,1)
                if(c=="^"){ i+=2 }                                    # ^ + mapping-quality char (read start)
                else if(c=="\$"){ i++ }                               # read end
                else if(c=="+"||c=="-"){                              # indel: skip sign + count + that many bases
                  i++; num=""
                  while(i<=n && substr(s,i,1) ~ /[0-9]/){ num=num substr(s,i,1); i++ }
                  i+=(num+0)
                }
                else if(c=="."||c==","){ RO_G++; i++ }                # match to reference
                else if(c ~ /[ACGTacgt]/){ AO_G++; i++ }             # mismatch = alt observation
                else{ i++ }                                          # * (del), <>, N, etc. -> ignore
              }
            }
            {
              chrom=\$1; pos=\$2; ref=\$3; dp=(\$4+0); dp_baq=(\$10+0);
              count_bases(\$5)                                        # RO/AO from the no-BAQ pileup (real base evidence)
              wt=(dp>=MINCOV?1:0);
              nc=(dp>=MINCOV?0:1);
              gt=(dp>=MINCOV?"0":"./.");
              flt=(dp>GAPDP && dp_baq<=GAPDP)?"baq_dropout":".";
              info="ADP="dp";WT="wt";HET=0;HOM=0;NC="nc;
              print chrom, pos, ".", ref, ".", ".", flt, info, "GT:DP:AD", gt":"dp":"RO_G","AO_G
            }' | cat header_template.txt - | bgzip -@ !{task.cpus} -c > backbone.vcf.gz
    else
        samtools mpileup -aa -f !{ref_fa} -Q !{params.allpos_min_bq} -d !{params.allpos_max_depth} !{bam} \\
            | awk -v MINCOV=!{params.allpos_min_cov} -v OFS="\\t" '
                function count_bases(s,   i,n,c,num){
                  RO_G=0; AO_G=0; i=1; n=length(s)
                  while(i<=n){
                    c=substr(s,i,1)
                    if(c=="^"){ i+=2 }
                    else if(c=="\$"){ i++ }
                    else if(c=="+"||c=="-"){
                      i++; num=""
                      while(i<=n && substr(s,i,1) ~ /[0-9]/){ num=num substr(s,i,1); i++ }
                      i+=(num+0)
                    }
                    else if(c=="."||c==","){ RO_G++; i++ }
                    else if(c ~ /[ACGTacgt]/){ AO_G++; i++ }
                    else{ i++ }
                  }
                }
                {
                  chrom=\$1; pos=\$2; ref=\$3; dp=(\$4+0);
                  count_bases(\$5)
                  wt=(dp>=MINCOV?1:0);
                  nc=(dp>=MINCOV?0:1);
                  gt=(dp>=MINCOV?"0":"./.");
                  info="ADP="dp";WT="wt";HET=0;HOM=0;NC="nc;
                  print chrom, pos, ".", ref, ".", ".", ".", info, "GT:DP:AD", gt":"dp":"RO_G","AO_G
                }' | cat header_template.txt - | bgzip -@ !{task.cpus} -c > backbone.vcf.gz
    fi

    safe_tabix backbone.vcf.gz
    rm -f mpileup_nobaq.txt mpileup_baq.txt
    """

    stub:
    """
    touch backbone.vcf.gz
    touch backbone.vcf.gz.tbi
    touch header_template.txt
    """
}

process MERGE_VCFS {
    tag "Merge: ${sampleId}"
    
    // Use getSampleDir for nested output support
    publishDir "${params.outdir}/${getSampleDir(sampleId, params)}", mode: params.publish_mode, saveAs: { filename -> getSavePath(filename, params) }
    
    // bcftools (concat/sort/view) on a bacterial VCF is light and mostly single-threaded.
    cpus 2
    memory { 4.GB * task.attempt }

    input:
    // outLabel is "" for the normal (masked) all.pos and ".raw" for the virgin one.
    tuple val(sampleId), val(refId), path(snps_vcf), path(snps_tbi), path(back_vcf), path(back_tbi), path(hdr_template), val(outLabel)

    output:
    tuple val(sampleId), val(refId), path("${sampleId}.${refId}${outLabel}.all.pos.vcf.gz"), path("${sampleId}.${refId}${outLabel}.all.pos.vcf.gz.tbi"), emit: allpos
    tuple val(sampleId), val(refId), path("${sampleId}.${refId}${outLabel}.vcf.gz"), path("${sampleId}.${refId}${outLabel}.vcf.gz.tbi"), emit: main_vcf
    
    shell:
    """
    set -euo pipefail
    if [ -n "\${SLURM_TMPDIR:-}" ]; then export TMPDIR="\$SLURM_TMPDIR"; else export TMPDIR="./tmp_merge"; fi
    mkdir -p "\$TMPDIR"
    
    safe_tabix () {
      local gz="\$1"; local idx="\${gz}.tbi"
      set +e; tabix -f -p vcf "\$gz"; st=\$?; set -e
      if [ \$st -ne 0 ]; then
        # Only a genuinely header-only VCF gets an empty index. bcftools reads the whole
        # file, so a truncated/malformed VCF errors under set -e (fails loud) instead of
        # getting a fake index that would let a short VCF flow into the consensus.
        local nrec; nrec=\$(bcftools view -H "\$gz" | wc -l)
        if [ "\$nrec" -ne 0 ]; then
          echo "ERROR: tabix failed on \$gz with \$nrec records (malformed/truncated) -> refusing a fake index" >&2
          exit \$st
        fi
        : > "\$idx"
      fi
    }

    # 1. Apply proper header to SNP VCF
    bcftools reheader -h !{hdr_template} !{snps_vcf} > clean_snps.vcf.gz
    safe_tabix clean_snps.vcf.gz
    
    # Save a copy as the main variant VCF
    cp clean_snps.vcf.gz !{sampleId}.!{refId}!{outLabel}.vcf.gz
    cp clean_snps.vcf.gz.tbi !{sampleId}.!{refId}!{outLabel}.vcf.gz.tbi

    # 2. Merge Backbone + SNPs
    # Count with bcftools (reads the whole file): a truncated clean_snps errors under set -e
    # instead of the old 'zgrep | head | wc || true', where SIGPIPE+pipefail could swallow a
    # truncation and yield n_vars=0 -> backbone-only branch -> every SNP silently dropped and
    # the consensus collapses to the reference.
    n_vars=\$(bcftools view -H clean_snps.vcf.gz | wc -l)
    echo "DEBUG: Number of variants found: \$n_vars" >&2

    if [[ "\$n_vars" -gt 0 ]]; then

      # A. Identify Variant Positions
      bcftools query -f '%CHROM\\t%POS\\n' clean_snps.vcf.gz | sort -u > variant_positions.txt
      
      # B. Filter those positions OUT of the backbone (to avoid duplication)
      bcftools view -T ^variant_positions.txt !{back_vcf} -Oz -o backbone_clean.vcf.gz
      safe_tabix backbone_clean.vcf.gz
      
      # C. Concatenate and Sort
      bcftools concat -a backbone_clean.vcf.gz clean_snps.vcf.gz \\
        | bcftools sort -T "\$TMPDIR" -m !{params.bcftools_sort_mem} -Oz -o !{sampleId}.!{refId}!{outLabel}.all.pos.vcf.gz
    else
      # If no variants, the all.pos VCF is identical to the backbone
      cp !{back_vcf} !{sampleId}.!{refId}!{outLabel}.all.pos.vcf.gz
    fi

    # The all.pos VCF is the direct consensus substrate and always has records; verify BGZF
    # integrity and index with plain tabix so any truncation/corruption fails the task loudly
    # (never a fake empty index here).
    bgzip -t !{sampleId}.!{refId}!{outLabel}.all.pos.vcf.gz
    tabix -f -p vcf !{sampleId}.!{refId}!{outLabel}.all.pos.vcf.gz
    """

    stub:
    """
    touch ${sampleId}.${refId}${outLabel}.all.pos.vcf.gz
    touch ${sampleId}.${refId}${outLabel}.all.pos.vcf.gz.tbi
    touch ${sampleId}.${refId}${outLabel}.vcf.gz
    touch ${sampleId}.${refId}${outLabel}.vcf.gz.tbi
    """
}

process CALL_FREEBAYES_RAW {
    // Lean FreeBayes for the VIRGIN path (virgin_full_freebayes=true): the SAME SNP calls but on
    // the UNfiltered bam, with no legacy outputs and no publishing -> feeds only the virgin
    // consensus so it shows the pre-filter variant set. Mirrors CALL_FREEBAYES steps 2/3/5/6/8.
    tag "FreeBayes(raw): ${sampleId}"
    // FreeBayes is single-threaded unless region-parallel is enabled.
    cpus { params.freebayes_parallel ? (params.freebayes_parallel_jobs as int) : 2 }
    memory { 8.GB * task.attempt }

    input:
    tuple val(sampleId), val(refId), path(bam), path(bai), path(ref_fa), path(exclude_txt)

    output:
    tuple val(sampleId), val(refId), path("valid_snps_formatted.vcf.gz"), path("valid_snps_formatted.vcf.gz.tbi"), emit: snps

    shell:
    """
    set -euo pipefail
    if [ -n "\${SLURM_TMPDIR:-}" ]; then export TMPDIR="\$SLURM_TMPDIR"; else export TMPDIR="./tmp_fbraw"; fi
    mkdir -p "\$TMPDIR"

    # Keep only numeric-coordinate rows -> strips EVERY header (the combined nucmer+genmap exclude has two)
    awk 'tolower(\$1)!="chrom" && \$2 ~ /^[0-9]+\$/ && \$3 ~ /^[0-9]+\$/ {print \$1"\\t"\$2"\\t"\$3}' !{exclude_txt} > exclude.regions.tsv || true
    EXCL_ARG=""
    if [[ "!{params.exclude_repeats}" == "true" && -s exclude.regions.tsv ]]; then EXCL_ARG="-T ^exclude.regions.tsv"; fi

    if [[ "!{params.freebayes_parallel}" == "true" ]]; then
      samtools faidx !{ref_fa}
      fasta_generate_regions.py !{ref_fa}.fai !{params.freebayes_parallel_chunk} > fb_regions.txt
      freebayes-parallel fb_regions.txt !{task.cpus} \\
        -f !{ref_fa} -p !{params.freebayes_ploidy} \\
        --min-alternate-count !{params.freebayes_min_alt_count} \\
        --min-alternate-fraction !{params.freebayes_min_alt_fraction} \\
        --min-mapping-quality !{params.freebayes_min_map_qual} \\
        --min-base-quality !{params.freebayes_min_base_qual} \\
        !{bam} > raw_freebayes.vcf
    else
      freebayes -f !{ref_fa} -p !{params.freebayes_ploidy} \\
        --min-alternate-count !{params.freebayes_min_alt_count} \\
        --min-alternate-fraction !{params.freebayes_min_alt_fraction} \\
        --min-mapping-quality !{params.freebayes_min_map_qual} \\
        --min-base-quality !{params.freebayes_min_base_qual} \\
        !{bam} > raw_freebayes.vcf
    fi

    bcftools norm -m - -a -f !{ref_fa} raw_freebayes.vcf | bcftools view -e 'GT="0/0" || GT="0"' -Ov -o normalized.vcf

    base_filt="(FMT/DP >= !{params.filter_min_dp} || INFO/DP >= !{params.filter_min_dp}) && INFO/SAF[0] >= !{params.min_alt_fwd} && INFO/SAR[0] >= !{params.min_alt_rev}"
    af_guard="(INFO/AO[0] + INFO/RO) > 0"
    af_calc="((INFO/AO[0] * 1.0) / (INFO/AO[0] + INFO/RO))"
    hom_rule="\$base_filt && \$af_guard && \$af_calc >= !{params.hom_threshold}"
    het_rule="\$base_filt && \$af_guard && \$af_calc >= !{params.het_min_frac} && \$af_calc < !{params.hom_threshold}"

    bcftools view \$EXCL_ARG --types snps -i "\$hom_rule || \$het_rule" normalized.vcf > valid_snps.vcf || true
    [ -s valid_snps.vcf ] || bcftools view -h normalized.vcf > valid_snps.vcf

    awk -v HETMIN=!{params.het_min_frac} 'BEGIN{OFS="\\t"}
      /^##/ { print; next }
      /^#CHROM/ { print; next }
      {
        dp=0; ro=0; ao=0
        ni=split(\$8,I,";")
        for(k=1;k<=ni;k++){
          if(I[k] ~ /^DP=/){ v=I[k]; sub(/^DP=/,"",v); dp=v+0 }
          else if(I[k] ~ /^RO=/){ v=I[k]; sub(/^RO=/,"",v); ro=v+0 }
          else if(I[k] ~ /^AO=/){ v=I[k]; sub(/^AO=/,"",v); split(v,ax,","); ao=ax[1]+0 }
        }
        if (dp==0) { n=split(\$9,fmt,":"); m=split(\$10,dat,":"); for(i=1;i<=n;i++) if(fmt[i]=="DP") dp=dat[i] }
        if (dp==0) dp=ro+ao
        if (dp==0) dp=1
        split(\$10,b,":"); gtype=b[1]
        if ((gtype=="0/1" || gtype=="1/0") && (ro+ao)>0 && ro/(ro+ao) < HETMIN) gtype="1/1"
        het=0; hom=0
        if (gtype=="0/1" || gtype=="1/0") het=1
        else if (gtype=="1/1") hom=1
        \$8="ADP="dp";WT=0;HET="het";HOM="hom";NC=0"
        \$9="GT:DP"
        \$10=gtype":"dp
        print
      }' valid_snps.vcf | bgzip -@ !{task.cpus} -c > valid_snps_formatted.vcf.gz

    set +e; tabix -f -p vcf valid_snps_formatted.vcf.gz; st=\$?; set -e
    if [ \$st -ne 0 ]; then : > valid_snps_formatted.vcf.gz.tbi; fi
    """

    stub:
    """
    touch valid_snps_formatted.vcf.gz
    touch valid_snps_formatted.vcf.gz.tbi
    """
}
