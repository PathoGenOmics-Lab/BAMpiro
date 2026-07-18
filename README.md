<p align="center">
  <img src=".github/bampiro2.png" height="200" alt="BAMpiro logo" />
</p>

<div align="center">

[![License: GPL v3](https://img.shields.io/badge/license-GPL%20v3-%23af64d1?style=flat-square)](LICENSE)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A524.04.2-%2323aa62?style=flat-square)](https://www.nextflow.io/)
[![Version](https://img.shields.io/badge/version-1.0.1-%23149389?style=flat-square)](CHANGELOG.md)
[![Container](https://img.shields.io/badge/container-paururo%2Fbampiro-%232496ed?style=flat-square)](https://hub.docker.com/r/paururo/bampiro)
[![PGO](https://img.shields.io/badge/PathoGenOmics-lab-%23E52421?style=flat-square)](https://github.com/PathoGenOmics-Lab)
[![Docs](https://img.shields.io/badge/docs-online-%23149389?style=flat-square)](https://pathogenomics-lab.github.io/BAMpiro/)
[![Live QC report](https://img.shields.io/badge/QC%20report-live%20demo-%23af64d1?style=flat-square)](https://pathogenomics-lab.github.io/BAMpiro/examples/qc_report_demo.html)

**General bacterial short-read mapping, variant calling & lineage/DR typing.**
**Nextflow (DSL2) · containerized · self-contained interactive QC reports.**

[Docs site](https://pathogenomics-lab.github.io/BAMpiro/) · [Live QC report](https://pathogenomics-lab.github.io/BAMpiro/examples/qc_report_demo.html) · [Quick Start](#quick-start) · [Configuration](docs/configuration.md) · [Citation](#citation)

</div>

__Paula Ruiz-Rodriguez<sup>1</sup>__
__and Mireia Coscolla<sup>1</sup>__
<br>
<sub> 1. I<sup>2</sup>SysBio, University of Valencia-CSIC, FISABIO Joint Research Unit Infection and Public Health, Valencia, Spain </sub>

---

## What is BAMpiro?

> [!TIP]
> **New here?** Start with the [Quick Start](docs/quickstart.md) - a run command, the
> samplesheet format, and where to find the report. Every parameter is listed in the
> [configuration reference](docs/configuration.md).

**BAMpiro** is a modular, containerized **Nextflow (DSL2)** pipeline that takes raw
bacterial short reads all the way to annotated variants, consensus sequences, and a
single **interactive QC report**. It is tuned by default for *Mycobacterium
tuberculosis* but is **organism-agnostic** - point it at any reference genome + GFF.

**Main features:**

- Reference-agnostic mapping, variant calling (`FreeBayes`), and consensus
- Alignment-free MTBC lineage + WHO drug-resistance typing ([Pathotypr](docs/pathotypr.md))
- A self-contained, interactive [HTML QC report](docs/qc-report.md) with 21 panels
- Dual amino-acid numbering (used reference + H37Rv / Mycobrowser)
- One pinned container with every tool and marker panel built in

## Features

| Feature | Description |
|---|---|
| 🧬 Any bacterial genome | Reference-agnostic mapping + variant calling; TB-tuned defaults, works on any species |
| 🧹 Repeat & mappability masking | `nucmer` repeat exclusion plus a length-aware `genmap` read filter |
| 🧪 Variants & backbone | `FreeBayes` (ploidy 1/2) + "all-sites" VCFs for phylogenetic supermatrices |
| 🩺 Lineage & drug resistance | Alignment-free MTBC lineage + WHO DR typing ([Pathotypr](docs/pathotypr.md)), reference-agnostic |
| 📊 Interactive QC report | Self-contained HTML dashboard, [21 panels](docs/qc-report.md) + per-sample `qc_flags.tsv` |
| 🔤 Dual amino-acid numbering | Protein changes in both the used reference and H37Rv/Mycobrowser numbering |
| 📦 Containerized & reproducible | A single pinned image with every tool + bundled marker panels ([details](docs/installation.md)) |
| ⚙️ Fully configurable | Every step exposed as a Nextflow parameter ([reference](docs/configuration.md)) |

## Installation

Requires **Nextflow ≥ 24.04.2** and **Docker** or **Singularity**. The pipeline
pulls a pinned `paururo/bampiro` image with every tool built in - nothing else to
install.

```bash
nextflow run main.nf --tsv samples.tsv --outdir results_bampiro -profile standard
```

Full requirements and the bundled software versions: [Installation](docs/installation.md).

## Quick Start

BAMpiro assumes *M. tuberculosis* settings by default (ploidy = 2 for mixed
infections). Lineage/DR typing is **off** by default - enable it (and dual
amino-acid numbering) with:

```bash
nextflow run main.nf \
    --tsv samples.tsv --outdir results_bampiro -profile standard \
    --run_pathotypr true --annotate_canonical true
```

The samplesheet is a TSV (`sampleId`, `r1`, `r2`, `refId`, `refFasta`, `refGff`, …).
See the [Quick Start guide](docs/quickstart.md) for the full format and multi-run
merging, and [Outputs](docs/outputs.md) for the result layout.

**New to BAMpiro?** Follow the hands-on
[**Jupyter tutorial**](docs/tutorial/bampiro_tutorial.ipynb) — install → samplesheet →
run → explore the outputs with `pandas` / `matplotlib` on a bundled 17-sample example
cohort (runnable without running the pipeline first).

## Interactive QC Report

Every run writes a single self-contained `<samplesheet>_qc_report.html` (no internet,
no CDN) that folds the whole cohort into one dashboard: **live-adjustable QC
thresholds**, a **dark / light theme**, a collapsible sidebar, and **21 linked
panels** - general statistics, flagged samples, per-lineage summary (canonical
*mycolorsTB* palette), QC-space PCA, genome landscape, SNP dynamics, epistasis, a
full SNP matrix, and drug resistance - plus a machine-readable per-sample
`qc_flags.tsv`. It is organism-agnostic, mobile-responsive, and works offline on an
HPC login node.

> [!TIP]
> **Try it live:** open the
> [**interactive example report**](https://pathogenomics-lab.github.io/BAMpiro/examples/qc_report_demo.html)
> — a full 17-sample demo cohort, every panel populated, right in your browser (the
> very same self-contained HTML each run produces).

Full panel list and interactive features: [Interactive QC Report](docs/qc-report.md).

## Documentation

| Document | Description |
| :--- | :--- |
| [Introduction](docs/introduction.md) | What BAMpiro is, key features, and the workflow at a glance |
| [Installation](docs/installation.md) | Requirements, the container, and bundled software versions |
| [Quick Start](docs/quickstart.md) | Run commands, the samplesheet format, and multi-run merging |
| [Configuration](docs/configuration.md) | The full parameter reference |
| [Interactive QC Report](docs/qc-report.md) | The self-contained HTML dashboard and its 21 panels |
| [Lineage & Drug-Resistance Typing](docs/pathotypr.md) | Pathotypr typing and dual amino-acid numbering |
| [Outputs](docs/outputs.md) | The result file tree and the repository layout |
| [Changelog](CHANGELOG.md) | Version history |

Browse the whole set under [`docs/`](docs/), read it online at
[pathogenomics-lab.github.io/BAMpiro](https://pathogenomics-lab.github.io/BAMpiro/), or
build the searchable site locally with `make docs-serve`
(needs [MkDocs Material](https://squidfunk.github.io/mkdocs-material/):
`pip install -r docs/requirements.txt`).

## Why "BAMpiro"?

The name is a play on words combining bioinformatics and folklore:

- **BAM** - Binary Alignment Map, the standard format for reads aligned to a
  reference genome; the "heart" of this pipeline (mapping → variant calling).
- **Piro** - combined with "BAM" it sounds like *Vampiro* (Spanish/Portuguese for
  vampire).

Just as a vampire seeks blood, BAMpiro seeks BAM files (and FASTQ data) to extract
vital information - variants, lineages, and stats. A creature that lives in your
cluster and processes bacterial genomes.

## Citation

If you use BAMpiro in your research, please cite:

> Ruiz-Rodriguez P, Coscollá M. **BAMpiro: a Nextflow pipeline for bacterial
> short-read mapping, variant calling and lineage/drug-resistance typing.**
> https://github.com/PathoGenOmics-Lab/BAMpiro

```bibtex
@software{ruiz-rodriguez_bampiro,
  title   = {BAMpiro: bacterial short-read mapping, variant calling and lineage/drug-resistance typing},
  author  = {Ruiz-Rodriguez, Paula and Coscoll{\'a}, Mireia},
  url      = {https://github.com/PathoGenOmics-Lab/BAMpiro},
  version = {1.0.1},
  license = {GPL-3.0}
}
```

## License

[GNU General Public License v3.0](LICENSE)

---

<h2 id="contributors" align="center">✨ <a href="https://github.com/PathoGenOmics-Lab/BAMpiro/graphs/contributors">Contributors</a></h2>

<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->
<div align="center">
BAMpiro is developed with ❤️ by:
<table>
  <tr>
    <td align="center">
      <a href="https://github.com/paururo">
        <img src="https://avatars.githubusercontent.com/u/50167687?v=4&s=100" width="100px;" alt=""/>
        <br />
        <sub><b>Paula Ruiz-Rodriguez</b></sub>
      </a>
      <br />
      <a href="" title="Code">💻</a>
      <a href="" title="Research">🔬</a>
      <a href="" title="Ideas">🤔</a>
      <a href="" title="Data">🔣</a>
      <a href="" title="Desing">🎨</a>
      <a href="" title="Tool">🔧</a>
    </td> 
    <td align="center">
      <a href="https://github.com/mireiacoscolla">
        <img src="https://avatars.githubusercontent.com/u/29301737?v=4&s=100" width="100px;" alt=""/>
        <br />
        <sub><b>Mireia Coscolla</b></sub>
      </a>
      <br />
      <a href="https://www.uv.es/instituto-biologia-integrativa-sistemas-i2sysbio/es/investigacion/proyectos/proyectos-actuales/mol-tb-host-1286169137294/ProjecteInves.html?id=1286289780236" title="Funding/Grant Finders">🔍</a>
      <a href="" title="Ideas">🤔</a>
      <a href="" title="Mentoring">🧑‍🏫</a>
      <a href="" title="Research">🔬</a>
      <a href="" title="User Testing">📓</a>
    </td> 
  </tr>
</table>

This project follows the [all-contributors](https://github.com/all-contributors/all-contributors) specification ([emoji key](https://allcontributors.org/docs/en/emoji-key)).

<!-- markdownlint-restore -->
<!-- prettier-ignore-end -->

<!-- ALL-CONTRIBUTORS-LIST:END -->
