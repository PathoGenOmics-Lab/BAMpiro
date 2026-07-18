# BAMpiro Documentation

Detailed documentation for **BAMpiro** - a general bacterial short-read mapping,
variant-calling and lineage / drug-resistance typing pipeline. Read it online at
**[pathogenomics-lab.github.io/BAMpiro](https://pathogenomics-lab.github.io/BAMpiro/)**,
browse the pages here on GitHub, or build the site locally (see below).

## Contents

| Document | Description |
| :--- | :--- |
| [Introduction](introduction.md) | What BAMpiro is, key features, and the workflow at a glance |
| [Installation](installation.md) | Requirements, the container, and bundled software versions |
| [Quick Start](quickstart.md) | Run commands, the samplesheet format, and multi-run merging |
| [Configuration](configuration.md) | The full parameter reference and feature toggles |
| [Troubleshooting](troubleshooting.md) | Common first-run errors and how to fix them |
| [Interactive QC Report](qc-report.md) | The self-contained HTML dashboard and its panels |
| [Lineage & Drug-Resistance Typing](pathotypr.md) | Pathotypr typing and dual amino-acid numbering |
| [Tutorial](tutorial/bampiro_tutorial.ipynb) | An end-to-end walkthrough (Jupyter notebook) |
| [Outputs](outputs.md) | The result file tree and the repository layout |
| [From outputs to a phylogeny](downstream.md) | Building a tree from the pipeline outputs |
| [Changelog](https://github.com/PathoGenOmics-Lab/BAMpiro/blob/indel-mask/CHANGELOG.md) | Version history |

## Building the docs locally

The site is built with [MkDocs Material](https://squidfunk.github.io/mkdocs-material/).
The Markdown lives in `docs/*.md` (readable straight from GitHub) and the site
configuration is [`mkdocs.yml`](https://github.com/PathoGenOmics-Lab/BAMpiro/blob/indel-mask/mkdocs.yml)
at the repository root.

```bash
pip install -r docs/requirements.txt
make docs-serve    # live preview at http://127.0.0.1:8000
make docs          # build the static site into ./site
```

Every push to `main` / `indel-mask` that touches the docs rebuilds and validates the
site via
[`.github/workflows/docs.yml`](https://github.com/PathoGenOmics-Lab/BAMpiro/blob/indel-mask/.github/workflows/docs.yml)
(and uploads it as a downloadable artifact). Publishing to GitHub Pages activates
automatically once the repo is public and the `ENABLE_PAGES` variable is set - see the
comments at the top of that workflow.
