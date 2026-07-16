# Use minimal Debian image (amd64 architecture enforced)
FROM --platform=linux/amd64 debian:bullseye-slim

LABEL maintainer="paula.ruiz.rodriguez@csic.es"
LABEL description="Self-contained Docker image for the BAMpiro Nextflow pipeline (Short-Reads)"

WORKDIR /app

# 1. Install system dependencies
# procps is needed for Nextflow monitoring
# git/wget/tar/bzip2 are standard utilities
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    wget \
    curl \
    bzip2 \
    ca-certificates \
    procps \
    git \
    g++ \
    make \
    zlib1g-dev \
    && apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# 2. Install Micromamba (Faster alternative to Conda)
RUN wget -qO- https://micromamba.snakepit.net/api/micromamba/linux-64/latest \
    | tar -xvj -C /usr/local/bin --strip-components=1 bin/micromamba

ENV MAMBA_ROOT_PREFIX=/opt/conda

# 3. Create Conda environment with BAMpiro tools
# Channels: conda-forge (general), bioconda (bioinformatics)
# Tools:
# - nextflow: Workflow manager
# - bwa: Mapper
# - fastp: QC and trimming
# - samtools: BAM manipulation
# - freebayes: Variant calling (bacteria/ploidy 1)
# - snpeff: Annotation
# - kraken2: Decontamination
# - multiqc: Reporting
# - pathotypr: alignment-free MTBC lineage + WHO drug-resistance genotyping (k-mer, reference-agnostic)
# - python=3.12: Latest stable python requested
RUN micromamba create -n bampiro-env -c conda-forge -c bioconda -y \
    'python' \
    'openjdk' \
    'nextflow' \
    'fastp' \
    'bwa-mem2' \
    'samtools' \
    'freebayes' \
    'bcftools' \
    'genmap' \
    'snpeff' \
    'kraken2' \
    'multiqc' \
    'pathotypr' \
    'pandas' \
    'biopython' \
    && micromamba clean --all --yes

# Put the env binaries on PATH up front so the RUN steps below (snpEff, pathotypr) can find them.
ENV PATH=/opt/conda/envs/bampiro-env/bin:/opt/conda/bin:$PATH

# 3b. Bundle the pathotypr marker panels + pre-trained model from Zenodo (v1.0.0, pinned) so runs are
# reproducible and offline. Positions are in MTBC-ancestor coordinates; pathotypr is k-mer based so it is
# reference-agnostic. Override the paths with --pathotypr_* if you supply your own.
ENV PATHOTYPR_DATA=/opt/pathotypr
RUN mkdir -p ${PATHOTYPR_DATA} && \
    wget -q -O ${PATHOTYPR_DATA}/lineage_markers.tsv "https://zenodo.org/records/19210044/files/pathotypr_lineage_markers_v1.0.0.tsv?download=1" && \
    wget -q -O ${PATHOTYPR_DATA}/dr_markers.tsv       "https://zenodo.org/records/19210044/files/pathotypr_dr_markers_v1.0.0.tsv?download=1" && \
    wget -q -O ${PATHOTYPR_DATA}/rf_model.pathotypr   "https://zenodo.org/records/19210044/files/pathotypr_rf_model_v1.0.0.pathotypr?download=1" && \
    wget -q -O ${PATHOTYPR_DATA}/reference.fasta      "https://zenodo.org/records/3497110/files/MTB_ancestor_reference.fasta?download=1"

# 3c. Pre-download the H37Rv snpEff database so variants can also be annotated in H37Rv / Mycobrowser
# numbering for the report's dual amino-acid display (see qc_report.py --vcfs-h37rv). Adjust the genome
# name if snpEff renames it. This needs network at build time only.
RUN snpEff download -v Mycobacterium_tuberculosis_h37rv || \
    echo "WARN: could not pre-fetch the H37Rv snpEff DB; build it manually or set the correct genome name"

# 4. Pipeline Setup
# Create directory structure
WORKDIR /pipeline

# Copy pipeline scripts into the container
# (Assumes you build this dockerfile from the root of your project)
COPY . .

# 6. Runtime Configuration
# Set working directory for data mounting (where user binds their data)
WORKDIR /data

# Default command
CMD ["/bin/bash"]
