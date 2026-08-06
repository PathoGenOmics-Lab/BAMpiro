"""Internals of the BAMpiro QC report.

`bin/qc_report.py` is the CLI; everything it does lives here:

  parsers.py  read the pipeline's files (summary TSV, consensus FASTA, GFF, BED,
              VCF, mapDamage, Kraken, samplesheet) into plain Python values
  metrics.py  the metric registry and the QC verdict engine (flag_sample)
  panels.py   the optional analysis panels built from parsed VCFs
  render.py   the front-end assets and the self-contained HTML assembly

Nothing is re-exported here on purpose: each name is imported from the module that
owns it, so the ownership stays visible at every call site.
"""
