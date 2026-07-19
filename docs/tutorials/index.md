# Tutorials

A guided, hands-on path through **BAMpiro** — from your first run to a finished
phylogeny. Each tutorial is short and task-oriented; follow them in order, or jump to the
one you need. For exhaustive detail, every tutorial links back to the reference pages.

!!! tip "New to BAMpiro?"
    Start with **[Your first run](first-run.md)** — you'll go from an empty folder to an
    open interactive QC report in one sitting.

<div class="grid cards" markdown>

-   :material-rocket-launch:{ .lg .middle } &nbsp; **[1 · Your first run](first-run.md)**

    ---

    Install the prerequisites, run the pipeline on a small cohort, and open the QC report.

-   :material-table:{ .lg .middle } &nbsp; **[2 · The samplesheet](samplesheet.md)**

    ---

    The required columns, multi-run merging, and the metadata columns that unlock report
    features (`treatment`, `dose`, `site`, timepoints…).

-   :material-tune:{ .lg .middle } &nbsp; **[3 · Configuring a run](configuring-a-run.md)**

    ---

    Profiles, key parameters, running on any organism, and turning on lineage / DR typing.

-   :material-chart-box:{ .lg .middle } &nbsp; **[4 · Reading the QC report](reading-the-qc-report.md)**

    ---

    A guided tour: triage, live thresholds, the metadata cohort filter, and the
    dose analyses (Dose × treatment, Variant × dose).

-   :material-dna:{ .lg .middle } &nbsp; **[5 · Lineage & drug-resistance typing](lineage-and-dr-typing.md)**

    ---

    Enable pathotypr, read the WHO-graded drug-resistance calls, and use dual amino-acid
    numbering.

-   :material-family-tree:{ .lg .middle } &nbsp; **[6 · Building a phylogeny](building-a-phylogeny.md)**

    ---

    Take the pipeline outputs and the QC exclusion set through to a tree.

-   :material-notebook:{ .lg .middle } &nbsp; **[Analysing outputs (Jupyter)](../tutorial/bampiro_tutorial.ipynb)**

    ---

    A fully runnable notebook: explore a real cohort's outputs with `pandas` /
    `matplotlib` — no pipeline run required.

</div>

!!! note "Following along without running the pipeline"
    The tutorials use a bundled **17-sample MTBC demo cohort**. You can read the whole path
    without a cluster; when you're ready, point the same commands at your own data. There's
    also a [live interactive demo report](../examples/qc_report_demo.html) to click through.
