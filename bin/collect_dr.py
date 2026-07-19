#!/usr/bin/env python3
"""Collect pathotypr split-fastq drug-resistance detail files into one run TSV for the QC report.

Each per-sample pathotypr `*_mutations.tsv` (from `pathotypr split-fastq -m <WHO DR markers>`) has columns:
    pos  ref_allele  alt_allele  ref_count  alt_count  alt_fraction  lineage_path  [extra_annotations]
The WHO drug-resistance marker panel has no empty column, so pathotypr reads its annotation columns as the
marker "levels" and joins them into `lineage_path` (';'-separated) as:
    drug ; resistance ; marker_name ; grade ; gene ; mutation

Output (consumed by qc_report.py --dr-report):
    sample  drug  gene  mutation  grade  marker_name  af  dp
Only markers whose alt_fraction passed pathotypr's --min-alt-percent are present (i.e. detected mutations).
The sample name is inferred from the file name (<sample>.dr_mutations.tsv) unless --sample is given.
"""
from __future__ import annotations
import argparse
import os
import re
import sys

_LP = ["drug", "resistance", "marker_name", "grade", "gene", "mutation"]


def sample_from_name(path):
    b = os.path.basename(path)
    b = re.sub(r'\.dr_mutations\.tsv$|_mutations\.tsv$|\.tsv$', '', b)
    # files are named <sample>.dr_mutations.tsv or <sample>__<runId>.dr_mutations.tsv (multi-lane / multi-ref);
    # '__' is the pipeline's sample/run separator and never appears inside a sample id, so take the part before it.
    return b.split('__')[0]


def main():
    ap = argparse.ArgumentParser(description="Aggregate pathotypr DR mutation files into a run DR TSV.")
    ap.add_argument('-o', '--output', required=True)
    ap.add_argument('--sample', default=None, help="sample name (else inferred from each file name)")
    ap.add_argument('files', nargs='*', help="per-sample pathotypr *_mutations.tsv from the DR marker run")
    args = ap.parse_args()

    rows = []
    for p in args.files:
        if not p or not os.path.exists(p):
            continue
        s = args.sample or sample_from_name(p)
        try:
            with open(p, encoding='utf-8', errors='replace') as fh:
                header = None
                for line in fh:
                    if not line.strip():
                        continue
                    c = line.rstrip('\n').split('\t')
                    if header is None:
                        header = [h.strip().lower() for h in c]
                        continue
                    d = dict(zip(header, c))
                    # maxsplit=5 keeps a ';' inside the last field (mutation, e.g. a compound descriptor)
                    parts = (d.get('lineage_path', '') or '').split(';', 5)
                    lp = dict(zip(_LP, parts))
                    try:
                        af = '%.4f' % (float(d.get('alt_fraction', '')) / 100.0)
                    except ValueError:
                        af = ''
                    try:
                        dp = str(int(d.get('ref_count', '0') or 0) + int(d.get('alt_count', '0') or 0))
                    except ValueError:
                        dp = ''
                    rows.append([s, lp.get('drug', ''), lp.get('gene', ''), lp.get('mutation', ''),
                                 lp.get('grade', ''), lp.get('marker_name', ''), af, dp])
        except OSError as e:
            sys.stderr.write("[collect_dr] WARN could not read %s: %s\n" % (p, e))

    with open(args.output, 'w') as out:
        out.write('\t'.join(['sample', 'drug', 'gene', 'mutation', 'grade', 'marker_name', 'af', 'dp']) + '\n')
        for r in rows:
            out.write('\t'.join(r) + '\n')
    sys.stderr.write("[collect_dr] %d resistance call(s) from %d file(s) -> %s\n" % (len(rows), len(args.files), args.output))


if __name__ == '__main__':
    main()
