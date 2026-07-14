# Use minimal Debian image (amd64 architecture enforced)
FROM --platform=linux/amd64 debian:bullseye-slim

LABEL maintainer="paula.ruiz.rodriguez@csic.es"
LABEL description="Self-contained Docker image for the BAMbard Nextflow pipeline (Short-Reads)"

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

# 3. Create Conda environment with BAMbard tools
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
# - python=3.12: Latest stable python requested
RUN micromamba create -n bambard-env -c conda-forge -c bioconda -y \
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
    'pandas' \
    'biopython' \
    && micromamba clean --all --yes

# 4. Set PATH to include the new environment binaries globally
# This ensures that calling 'bwa' or 'nextflow' works directly
ENV PATH=/opt/conda/envs/bambard-env/bin:/opt/conda/bin:$PATH

# 5. Pipeline Setup
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
