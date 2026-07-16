# Documentation

Detailed documentation for BAMpiro. Browse it here on GitHub, or build the
[mdbook](https://rust-lang.github.io/mdBook/) locally for a searchable docs site.

## Contents

| Document | Description |
| :--- | :--- |
| [Introduction](introduction.md) | What BAMpiro is, key features, and the workflow at a glance |
| [Installation](installation.md) | Requirements, the container, and bundled software versions |
| [Quick Start](quickstart.md) | Run commands, the samplesheet format, and multi-run merging |
| [Configuration](configuration.md) | The full parameter reference |
| [Interactive QC Report](qc-report.md) | The self-contained HTML dashboard and its 20 panels |
| [Lineage & Drug-Resistance Typing](pathotypr.md) | Pathotypr typing and dual amino-acid numbering |
| [Outputs](outputs.md) | The result file tree and the repository layout |

## Building locally

The mdbook sources live in [`src/`](src/). With
[mdbook](https://rust-lang.github.io/mdBook/) installed:

```bash
make docs          # build HTML → docs/book/
make docs-serve    # build + serve with live reload
```
