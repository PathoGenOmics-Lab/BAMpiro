# BAMpiro Documentation

Detailed documentation for **BAMpiro** - a general bacterial short-read mapping,
variant-calling and lineage / drug-resistance typing pipeline. Read it online at
**[pathogenomics-lab.github.io/BAMpiro](https://pathogenomics-lab.github.io/BAMpiro/)**,
browse the pages here on GitHub, or build the site locally (see below).

## Contents

<div class="grid cards" markdown>

-   :material-book-open-variant:{ .lg .middle } &nbsp; **[Introduction](introduction.md)**

    ---

    What BAMpiro is, key features, and the workflow at a glance.

-   :material-download:{ .lg .middle } &nbsp; **[Installation](installation.md)**

    ---

    Requirements, the container, and bundled software versions.

-   :material-rocket-launch:{ .lg .middle } &nbsp; **[Quick Start](quickstart.md)**

    ---

    Run commands, the samplesheet format, and multi-run merging.

-   :material-tune:{ .lg .middle } &nbsp; **[Configuration](configuration.md)**

    ---

    The full parameter reference and feature toggles.

-   :material-lifebuoy:{ .lg .middle } &nbsp; **[Troubleshooting](troubleshooting.md)**

    ---

    Common first-run errors and how to fix them.

-   :material-chart-box:{ .lg .middle } &nbsp; **[Interactive QC Report](qc-report.md)**

    ---

    The self-contained HTML dashboard and its panels.

-   :material-dna:{ .lg .middle } &nbsp; **[Lineage & Drug-Resistance Typing](pathotypr.md)**

    ---

    Pathotypr typing and dual amino-acid numbering.

-   :material-notebook:{ .lg .middle } &nbsp; **[Tutorial](tutorial/bampiro_tutorial.ipynb)**

    ---

    An end-to-end walkthrough (Jupyter notebook).

-   :material-file-tree:{ .lg .middle } &nbsp; **[Outputs](outputs.md)**

    ---

    The result file tree and the repository layout.

-   :material-family-tree:{ .lg .middle } &nbsp; **[From outputs to a phylogeny](downstream.md)**

    ---

    Building a tree from the pipeline outputs.

</div>

See the [Changelog](https://github.com/PathoGenOmics-Lab/BAMpiro/blob/indel-mask/CHANGELOG.md)
for the version history.

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
