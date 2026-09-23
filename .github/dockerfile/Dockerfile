# Use minimal Debian image (amd64 architecture enforced).
# bookworm, not bullseye: bullseye left LTS on 2026-08-31 and its security pool is being pruned
# while the package index still advertises what was in it, so apt resolves versions whose .deb
# files 404. The build failed on that, not on anything in this file.
FROM --platform=linux/amd64 debian:bookworm-slim

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
# From the project's own release channel rather than micromamba.snakepit.net, whose certificate
# expired on 2026-09-15: that host is an alias not covered by the *.mamba.pm wildcard it is served,
# so the download fails TLS verification. This removes a third-party CDN from the build entirely.
RUN wget -qO /usr/local/bin/micromamba \
      https://github.com/mamba-org/micromamba-releases/releases/latest/download/micromamba-linux-64 \
    && chmod +x /usr/local/bin/micromamba

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
    # freebayes-parallel ships inside the freebayes package but its dependencies do not:
    # it pipes `parallel` into vcffirstheader, vcfstreamsort and vcfuniq, and none of the
    # four were in the image. --freebayes_parallel is documented and off by default, so
    # turning it on failed inside the container with a command not found.
    'parallel' \
    'vcflib' \
    'bcftools' \
    'genmap' \
    'mummer4' \
    'snpeff' \
    'kraken2' \
    'multiqc' \
    # Pinned, like its marker files below. The catalogue and the binary carry the same
    # version on purpose: v1.0.2 renumbered the data files so the two match again.
    'pathotypr=1.0.2' \
    'pandas' \
    # pandas already pulls numpy in, but bin/gconv_model.py and bin/build_min_unique_len.py import
    # it directly, and a direct import has no business resting on somebody else's dependency.
    'numpy' \
    'biopython' \
    && micromamba clean --all --yes

# Put the env binaries on PATH up front so the RUN steps below (snpEff, pathotypr) can find them.
ENV PATH=/opt/conda/envs/bampiro-env/bin:/opt/conda/bin:$PATH

# 3b. Bundle the pathotypr marker panels + pre-trained model from Zenodo (v1.0.2, pinned) so runs are
# reproducible and offline. Positions are in MTBC-ancestor coordinates; pathotypr is k-mer based so it is
# reference-agnostic. Override the paths with --pathotypr_* if you supply your own.
#
# v1.0.2 is a CORRECTIVE release of the drug-resistance catalogue and replaces v1.0.0, which assigned
# each variant a single drug inherited from its gene and took the confidence grade from whichever
# catalogue row was found first. The WHO catalogue grades every variant-drug pair separately, so that
# was wrong in both directions: amikacin did not occur anywhere in v1.0.0, and 15,969 rows carried a
# grade belonging to another drug. Upgrading moves 9,439 rows into grade 1-2 and none out of it.
#
# v1.0.2 ships two coordinate frames. Take the ANCESTOR one: the reference bundled below is the MTBC
# ancestor, all 102,216 REF alleles of the ancestor file match it, and the H37Rv file mismatches 609
# of them at identical coordinates, so the wrong choice would drop those markers with no error at all.
# The lineage markers and the RF model are byte-identical to v1.0.0 (sha256 verified); only the
# resistance catalogue changes.
ENV PATHOTYPR_DATA=/opt/pathotypr
RUN mkdir -p ${PATHOTYPR_DATA} && cd ${PATHOTYPR_DATA} && \
    wget -nv --tries=3 --timeout=60 -O lineage_markers.tsv "https://zenodo.org/records/21915539/files/pathotypr_lineage_markers_v1.0.2.tsv?download=1" && \
    wget -nv --tries=3 --timeout=60 -O dr_markers.tsv       "https://zenodo.org/records/21915539/files/pathotypr_dr_markers_ancestor_v1.0.2.tsv?download=1" && \
    wget -nv --tries=3 --timeout=60 -O rf_model.pathotypr   "https://zenodo.org/records/21915539/files/pathotypr_rf_model_v1.0.2.pathotypr?download=1" && \
    wget -nv --tries=3 --timeout=60 -O reference.fasta      "https://zenodo.org/records/3497110/files/MTB_ancestor_reference.fasta?download=1" && \
    # Fail the build LOUDLY if an asset came down truncated or was replaced upstream. `test -s` used to
    # stand here and only caught an empty file, which is the one failure wget already reports itself.
    printf '%s\n' \
      '923214cc73bc10d9b4bb772a81e355f89a60921fe46efd0925b67a4966bc1a6c  lineage_markers.tsv' \
      '2de9cba89d23b6a8be1bf88aa3f3a3d52d32ecf9f81af647e0c08628e2cea773  dr_markers.tsv' \
      'f2eea62ae60164f796bf6d754cf498a79af9fadd23e7e4d94f716efa821fce65  rf_model.pathotypr' \
      'b7943e21f84b2129737c1badd6514b6a63d09523db722bb08df4d40e938bf833  reference.fasta' \
      > SHA256SUMS && \
    # --strict, because without it sha256sum treats a malformed digest line as a warning, skips
    # it and still exits 0 as long as the well-formed lines pass. A typo in one of the constants
    # above would then silently stop checking that file. Verified by building with a 65-character
    # digest: it passes without --strict and fails with it.
    sha256sum --strict -c SHA256SUMS && \
    { head -c1 reference.fasta | grep -q '>' || { echo "ERROR: reference.fasta is not a FASTA"; exit 1; }; } && \
    # Which catalogue produced a set of resistance calls is not recoverable from the calls
    # themselves, and v1.0.0 and v1.0.2 disagree on 15,969 rows. DUMP_VERSIONS reads this file so
    # every run records it, which is what the question "is this the current list?" needed and did
    # not have.
    printf '%s\n' 'v1.0.2 (Zenodo 10.5281/zenodo.21915539, MTBC-ancestor frame)' > VERSION

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
