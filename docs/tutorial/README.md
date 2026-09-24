# BAMpiro tutorial (Jupyter notebook)

A hands-on, end-to-end walkthrough: install → samplesheet → run → outputs → **exploring
the results with `pandas` / `matplotlib`** → the interactive QC report.

- **[`bampiro_tutorial.ipynb`](bampiro_tutorial.ipynb)** — open on GitHub to read it, or run it locally:

  ```bash
  pip install jupyterlab pandas matplotlib
  jupyter lab bampiro_tutorial.ipynb
  ```

The analysis cells are **fully runnable** against the small example cohort in
[`example_outputs/`](https://github.com/PathoGenOmics-Lab/BAMpiro/tree/main/docs/tutorial/example_outputs)
(17 samples, real pipeline-format TSVs — cohort summary, QC flags, drug-resistance
calls, SNP matrix, samplesheet metadata), so you can
follow the whole analysis without running the pipeline first. Point the same code at
your own run's `results_bampiro/<samplesheet>_*.tsv` to analyse a real cohort.

The example data is synthetic demo output; see the [main docs](../introduction.md) for the
real pipeline.
