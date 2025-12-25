nextflow.enable.dsl=2

/* ====================================================================
    VARIANTS MODULES
    Contains: FreeBayes Calling, Backbone Generation, and VCF Merging
==================================================================== */

process CALL_FREEBAYES {
    tag "FreeBayes: ${sampleId}"
    publishDir "${params.outdir}/${sampleId}", mode: 'copy', saveAs: { filename ->
        // We save specific VCFs but might discard intermediate temps if needed
        return filename
    }
    
    cpus 4
    memory '16 GB'

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
        # If indexing failed, check if file is truly empty or just has header
        if zgrep -v '^#' "\$gz" | grep -q "."; then exit \$st; else : > "\$idx"; fi
      fi
    }

    # 1. Handle Excluded Regions (Repeats)
    awk 'NR>1 && \$1!="" && \$2!="" && \$3!="" {print \$1"\\t"\$2"\\t"\$3}' !{exclude_txt} > exclude.regions.tsv || true
    EXCL_ARG=""
    if [[ "!{params.exclude_repeats}" == "true" && -s exclude.regions.tsv ]]; then 
        EXCL_ARG="-T ^exclude.regions.tsv"
    fi

    # 2. Run FreeBayes
    freebayes -f !{ref_fa} -p !{params.freebayes_ploidy} \\
      --min-alternate-count !{params.freebayes_min_alt_count} \\
      --min-alternate-fraction !{params.freebayes_min_alt_fraction} \\
      --min-mapping-quality !{params.freebayes_min_map_qual} \\
      --min-base-quality !{params.freebayes_min_base_qual} \\
      !{bam} > raw_freebayes.vcf

    # Compress Raw Output
    bgzip -@ !{task.cpus} -c raw_freebayes.vcf > !{sampleId}.!{refId}.freebayes.raw.vcf.gz
    safe_tabix !{sampleId}.!{refId}.freebayes.raw.vcf.gz

    # 3. Normalize Variants
    bcftools norm -m - -a -f !{ref_fa} -Ov -o normalized.vcf raw_freebayes.vcf
    
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
    # Adds ADP, WT, HET, HOM, NC tags to INFO/FORMAT for the consensus step
    awk 'BEGIN{OFS="\\t"}
      /^##/ { print; next }
      /^#CHROM/ { print; next }
      {
        dp=0
        if (match(\$8, /DP=[0-9]+/)) { s=substr(\$8,RSTART,RLENGTH); split(s,a,"="); dp=a[2] }
        if (dp==0) { n=split(\$9,fmt,":"); m=split(\$10,dat,":"); for(i=1;i<=n;i++) if(fmt[i]=="DP") dp=dat[i] }
        if (dp==0) dp=1
        split(\$10,b,":"); gtype=b[1]
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
}

process CALL_BACKBONE {
    tag "Backbone: ${sampleId}"
    cpus 4
    memory '16 GB'
    
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
      if [ \$st -ne 0 ]; then if zgrep -v '^#' "\$gz" | grep -q "."; then exit \$st; else : > "\$idx"; fi; fi
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
    echo -e "#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\tFORMAT\\t!{sampleId}" >> header_template.txt

    # 2. Run Mpileup for ALL positions (-aa)
    # Convert mpileup output to VCF lines using awk
    samtools mpileup -aa -f !{ref_fa} -Q !{params.allpos_min_bq} -d !{params.allpos_max_depth} !{bam} \\
        | awk -v MINCOV=!{params.allpos_min_cov} -v OFS="\\t" '{
            chrom=\$1; pos=\$2; ref=\$3; dp=(\$4+0);
            
            # Determine status based on coverage
            wt=(dp>=MINCOV?1:0); 
            nc=(dp>=MINCOV?0:1); 
            gt=(dp>=MINCOV?"0":"./.");
            
            info="ADP="dp";WT="wt";HET=0;HOM=0;NC="nc;
            
            # Print VCF Record (ALT is dot, QUAL is dot)
            print chrom, pos, ".", ref, ".", ".", ".", info, "GT:DP", gt":"dp
          }' | cat header_template.txt - | bgzip -@ !{task.cpus} -c > backbone.vcf.gz
    
    safe_tabix backbone.vcf.gz
    """
}

process MERGE_VCFS {
    tag "Merge: ${sampleId}"
    publishDir "${params.outdir}/${sampleId}", mode: 'copy', saveAs: { filename ->
        // Save the main VCF and the All-Positions VCF
        return filename
    }
    cpus 4
    memory '8 GB'
    
    input:
    tuple val(sampleId), val(refId), path(snps_vcf), path(snps_tbi), path(back_vcf), path(back_tbi), path(hdr_template)
    
    output:
    tuple val(sampleId), val(refId), path("${sampleId}.${refId}.all.pos.vcf.gz"), path("${sampleId}.${refId}.all.pos.vcf.gz.tbi"), emit: allpos
    tuple val(sampleId), val(refId), path("${sampleId}.${refId}.vcf.gz"), path("${sampleId}.${refId}.vcf.gz.tbi"), emit: main_vcf
    
    shell:
    """
    set -euo pipefail
    if [ -n "\${SLURM_TMPDIR:-}" ]; then export TMPDIR="\$SLURM_TMPDIR"; else export TMPDIR="./tmp_merge"; fi
    mkdir -p "\$TMPDIR"
    
    safe_tabix () {
      local gz="\$1"; local idx="\${gz}.tbi"
      set +e; tabix -f -p vcf "\$gz"; st=\$?; set -e
      if [ \$st -ne 0 ]; then : > "\$idx"; fi
    }

    # 1. Apply proper header to SNP VCF
    bcftools reheader -h !{hdr_template} !{snps_vcf} > clean_snps.vcf.gz
    safe_tabix clean_snps.vcf.gz
    
    # Save a copy as the main variant VCF
    cp clean_snps.vcf.gz !{sampleId}.!{refId}.vcf.gz
    cp clean_snps.vcf.gz.tbi !{sampleId}.!{refId}.vcf.gz.tbi

    # 2. Merge Backbone + SNPs
    has_snps=0
    if zgrep -v '^#' clean_snps.vcf.gz | grep -q "."; then has_snps=1; fi
    
    if [[ \$has_snps -eq 1 ]]; then
      # A. Identify Variant Positions
      bcftools query -f '%CHROM\\t%POS\\n' clean_snps.vcf.gz | sort -u > variant_positions.txt
      
      # B. Filter those positions OUT of the backbone (to avoid duplication)
      bcftools view -T ^variant_positions.txt !{back_vcf} -Oz -o backbone_clean.vcf.gz
      safe_tabix backbone_clean.vcf.gz
      
      # C. Concatenate and Sort
      bcftools concat -a backbone_clean.vcf.gz clean_snps.vcf.gz \\
        | bcftools sort -T "\$TMPDIR" -m !{params.bcftools_sort_mem} -Oz -o !{sampleId}.!{refId}.all.pos.vcf.gz
    else
      # If no variants, the all.pos VCF is identical to the backbone
      cp !{back_vcf} !{sampleId}.!{refId}.all.pos.vcf.gz
    fi
    
    safe_tabix !{sampleId}.!{refId}.all.pos.vcf.gz
    """
}
