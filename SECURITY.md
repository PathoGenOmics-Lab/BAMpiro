# Security Policy

BAMpiro is a research pipeline, not a network service. It runs on a laptop or an HPC cluster, under
your own account, on data you point it at. It listens on no port, has no users, no authentication and
no persistent state beyond the files it writes to `--outdir` and the Nextflow work directory. That
rules out most of what "vulnerability" usually means, so this policy is short and specific to what is
actually plausible here.

## Supported versions

Only the latest release is supported. Fixes go into the next release; there are no backports to
earlier tags.

| Version | Supported |
| :--- | :--- |
| Latest release (see [CHANGELOG.md](CHANGELOG.md) and the [releases page](https://github.com/PathoGenOmics-Lab/BAMpiro/releases)) | Yes |
| Anything older | No |

## Reporting a vulnerability

Please report privately rather than opening a public issue:

- **GitHub private vulnerability reporting**: the *Security* tab of this repository, then
  *Report a vulnerability*. This keeps the report visible only to the maintainers.
- **Email**: paula.ruiz.rodriguez@csic.es

Include what you found, how to reproduce it, and what an attacker would gain. If it involves an input
file, a minimal example helps more than a real dataset (and please do not send us patient data).

This is a small academic team with no on-call rotation, so we cannot promise a response time. We will
acknowledge your report as soon as we see it and tell you what we intend to do about it. If you would
like credit in the changelog entry for the fix, say so.

## What is realistically in scope

**The container's pinned dependencies.** Runs execute inside `paururo/bampiro`, which bundles
BWA-MEM2, Samtools, BCFtools, FreeBayes, SnpEff, FastP, Kraken2, MUMmer, Python and a Java runtime,
all pinned (versions are listed in [docs/installation.md](docs/installation.md)). A vulnerability in
one of those reaches BAMpiro users through the image. Tell us if a bundled tool needs an urgent
version bump.

Two related notes. `nextflow.config` references the image by its mutable version **tag**, so the same
tag can be rebuilt; re-pin to `docker://paururo/bampiro@sha256:<digest>` if you need byte-for-byte
reproducibility, as documented in the installation page. And `modules/qc.nf` falls back to
downloading `extract_kraken_reads.py` from `params.krakentools_url` over the network when the
vendored copy under `bin/` is missing, which is a supply-chain surface worth knowing about; the
vendored copy is there so the normal path never touches the network.

**Untrusted input files.** The pipeline parses whatever you hand it: the TSV samplesheet, FASTQ
reads, a reference FASTA, a GFF3, and optionally a Kraken2 database. Paths from the samplesheet reach
shell commands inside process scripts. Running BAMpiro on files from a source you do not trust is
therefore closer to running an untrusted script than to opening a data file. Treat a crafted
samplesheet or reference that causes something other than a clean parse error as a bug worth
reporting.

**The QC report is your data.** Each run writes a single self-contained
`<samplesheet>_qc_report.html` that embeds sample identifiers, per-sample metrics and any metadata
columns you supplied. It contains no remote references and phones nothing home, but it is not
anonymised. Publishing or sharing one shares whatever was in your samplesheet.

## What is out of scope

- Denial of service through resource exhaustion. Kraken2 alone asks for around 80 GB of RAM by
  design; an under-sized host being killed by the OOM reaper is a sizing problem, covered in
  [Troubleshooting](docs/troubleshooting.md).
- Anything requiring an attacker who already has your shell account, your cluster credentials or
  write access to the work directory.
- The security of Nextflow, Docker, Singularity or your scheduler. Report those upstream.

Ordinary bugs, wrong results and crashes are not security issues. Open a
[normal issue](https://github.com/PathoGenOmics-Lab/BAMpiro/issues/new/choose) for those.
