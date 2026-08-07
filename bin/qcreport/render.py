"""The report front-end: the CSS/JS/HTML assets and the self-contained HTML assembly.

The assets are read from bin/report_assets/ at IMPORT time, so importing this
module touches the disk. build_html() splices them into shell.html together with
the gzip+base64 payload, producing the single file the pipeline publishes.

SECTION_INFO lives here rather than with the metrics because it is report copy:
the text behind the (i) icon of each section of the HTML shell.
"""
from __future__ import annotations

import base64
import gzip
import html
import json
import os
import sys

REPO_URL = "https://github.com/PathoGenOmics-Lab/BAMpiro"   # surfaced in the report header + footer


# ============================================================================= CSS / JS / SHELL
# ---- report front-end assets (CSS / JS / HTML shell / logo) live next to bin/qc_report.py in
# report_assets/, one file per language so each stays editable with its own tooling; build_html()
# splices them into the single self-contained output HTML. ----
# NOTE: report_assets/ sits beside bin/qc_report.py, i.e. one level ABOVE this package, so the
# path below walks bin/qcreport/render.py -> bin/qcreport -> bin -> bin/report_assets. The extra
# dirname() is what keeps it resolving after the move out of bin/qc_report.py; drop it and every
# asset read fails at import time.
_ASSET_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "report_assets")


def _asset(name):
    with open(os.path.join(_ASSET_DIR, name), encoding="utf-8") as fh:
        return fh.read()


# CSS split into ordered modules under report_assets/css/ (base tokens -> layout -> components ->
# panels); concatenated in cascade order, so the result is identical to a single stylesheet.
_CSS_MODULES = ["css/01_base.css", "css/02_layout.css", "css/03_components.css", "css/04_panels.css"]
CSS = "".join(_asset(m) for m in _CSS_MODULES)

# The report front-end is one ES5 IIFE split into ordered modules under report_assets/js/ (one per
# area) purely for editability; they are concatenated verbatim here, so the running code is identical
# to a single file. Order matters: 01 opens the IIFE and sets up shared state; 14 wires events + closes it.
_JS_MODULES = [
    "js/01_prelude.js", "js/02_qcspace.js", "js/03_state.js", "js/04_helpers.js",
    "js/05_insights_qc.js", "js/06_insights_genome.js", "js/07_render_core.js", "js/08_curation.js",
    "js/09_genome_genes.js", "js/10_dynamics.js", "js/11_epistasis.js", "js/12_snpmatrix.js",
    "js/13_drug_kraken.js", "js/14_boot.js",
]
JS = "".join(_asset(m) for m in _JS_MODULES)

# BAMpiro header logo (bampiro2.png resized to ~90x100 and embedded so the report stays
# self-contained - no external image request). Regenerate from .github/bampiro2.png if the logo changes.
LOGO_DATA_URI = _asset("logo.b64")

SHELL = _asset("shell.html")   # named shell.html (not report_*.html) so .gitignore's report*.html rule can't swallow it


def build_html(title, payload):
    # Embed the payload gzip-compressed + base64 so the whole cohort travels in a light, self-contained
    # HTML (JSON deflates ~5-10x); the browser inflates it natively at load. base64 has no "</" so it
    # cannot break out of the <script>. Nothing is dropped for size - compression is what lets us keep it all.
    raw = json.dumps(payload, separators=(",", ":")).encode("utf-8")
    data_gz = base64.b64encode(gzip.compress(raw, 9)).decode("ascii")
    sys.stderr.write("[qc_report] payload %.2f MB JSON -> %.2f MB embedded (gzip+base64, %.0f%% smaller)\n"
                     % (len(raw) / 1048576.0, len(data_gz) / 1048576.0,
                        100.0 * (1.0 - len(data_gz) / max(1, len(raw)))))
    # Inject content first, then substitute __TITLE__ LAST so a user --title that happens to
    # contain a placeholder token (e.g. "__JS__") can never pull in the CSS/JS/logo/JSON blob.
    return (SHELL.replace("__CSS__", CSS).replace("__LOGO__", LOGO_DATA_URI)
                 .replace("__JS__", JS).replace("__JSON_GZ__", data_gz)
                 .replace("__REPO__", html.escape(REPO_URL))
                 .replace("__TITLE__", html.escape(title)))


# Per-panel explanation shown by the (i) icon in each section heading: what the analysis is, how to read
# it, and its caveats. Keyed by the section id in the HTML shell.
SECTION_INFO = {
    "gstats": "Per-sample sequencing and mapping QC: reads, duplication, mapped %, mean/median depth, "
              "breadth, variant counts and Ti/Tv. Each numeric cell is shaded by where the sample sits in "
              "that column's range; the QC column flags PASS/WARN/FAIL against the live thresholds. Tick a "
              "box to basket a sample for exclusion; click a name for its full profile.",
    "linsum": "Median QC metrics grouped by the assigned MTBC lineage, with the number of samples, %PASS and "
              "how many are 'mixed' (more than one lineage above the mixture cut-off). Lineage comes from "
              "pathotypr; a mixed call can indicate co-infection or contamination.",
    "dist": "Per-metric value distributions across the cohort (a strip/bee plot, or a ranked bar). The shaded "
            "band is the acceptable range for the current thresholds; points outside it are outliers. Useful "
            "to see cohort spread and choose sensible cut-offs.",
    "corr": "A scatter of any two QC metrics, one point per sample, coloured by QC status or lineage. Drag a "
            "box to basket the enclosed samples. Reveals metric relationships (e.g. depth vs breadth) and "
            "samples that are joint outliers.",
    "corrmatrix": "Pairwise Pearson correlation between the QC metrics across the cohort, as a heatmap. Shows "
                  "which metrics are redundant (strongly correlated) and which capture independent axes of quality.",
    "qcpca": "PCA of the standardised QC metrics: each point is a sample placed by its overall quality "
             "profile. Samples that cluster share a profile; a lone point is a multivariate outlier that no "
             "single metric flags.",
    "divcomp": "Each sample's genomic divergence (SNPs vs the reference) against its consensus completeness. "
               "Separates genuinely divergent samples from those that only look divergent because they are "
               "incomplete or contaminated.",
    "cons": "How much of the reference each sample's consensus recovers: callable %, missing % (N/gap), IUPAC "
            "(ambiguous) % and the longest gap. A phylogeny-oriented view of completeness, not just average depth.",
    "genome": "Per-position callability / variant density along the reference, binned into a heatmap per "
              "sample. Brush a region to list the genes under it. Reveals systematically low-callability "
              "regions (repeats, deletions) shared across samples.",
    "function": "The snpEff functional class of each sample's variants (HIGH/MODERATE/LOW/MODIFIER impact and "
                "effect types such as missense / synonymous). A per-sample mutational-impact profile.",
    "geneburden": "Genes carrying the most impactful (HIGH/MODERATE) variants across the cohort, with the "
                  "dominant effect and how many samples are hit. A cohort-level view of which genes accumulate "
                  "functional change.",
    "hotspots": "Genes with the highest SNP density across the cohort (variants per gene length). Flags the "
                "most variable loci - often repeats, antigens, or genes under diversifying selection.",
    "temporal": "The distribution of sampling dates, when the samplesheet provides them. Context for temporal "
                "analyses (e.g. SNP dynamics); not a QC metric itself.",
    "pnps": "Cohort-level pairwise dN/dS per gene (eskaks, Nei/Li). dN/dS > 1 suggests positive / diversifying "
            "selection at cohort scale. This is population genetics, NOT per-sample QC: noisy per gene, "
            "sensitive to alignment, and saturating at low divergence (NA when dS is near zero).",
    "adna": "For samples marked ancient: the characteristic post-mortem damage signal (elevated terminal 5' "
            "C to T). A sample below the damage floor may be a modern contaminant rather than authentic "
            "ancient DNA.",
    "dynamics": "For connected time-series (patient / passage line), each variant's allele-frequency "
                "trajectory over time, with per-timepoint read depth. Events flag emergence, fixation, loss "
                "and non-synonymous changes. Needs a time column and a group column in the samplesheet.",
    "epistasis": "Pairs of variants whose allele-frequency trajectories co-vary within a series: concordant "
                 "(rise/fall together) or discordant (one rises as the other falls). Scored by the Pearson "
                 "correlation, a permutation p-value and a Benjamini-Hochberg FDR q; a pattern recurring across "
                 "independent series is what makes a pair 'strong'. Candidate linked / co-selected / competing "
                 "SNPs, NOT proof of a functional interaction.",
    "snpmatrix": "Every SNP site (rows) by sample (columns); each cell is the allele frequency with its depth, "
                 "plus a reference column and samplesheet metadata as column-header levels. Filter by gene / "
                 "position or by metadata, and download the full matrix as a TSV.",
    "drug": "Resistance-associated mutations detected by pathotypr against the WHO catalogue (H37Rv numbering), "
            "as a sample x drug matrix (worst grade per drug) and a per-mutation table with the WHO confidence "
            "grade. Alignment-free (k-mer), so it works regardless of the mapping reference. A genomic screen, "
            "NOT a clinical drug-susceptibility result.",
    "gconv": "Stretches where one paralog appears to have been copied onto another: the acceptor locus stops "
             "carrying its own alleles and carries the donor's between two breakpoints. Finding a candidate is "
             "easy; deciding whether to believe it is the whole problem. There are three ways an acceptor site "
             "can show the donor's base and only one of them is a conversion: it was converted along with its "
             "neighbours, it mutated to that base on its own, or the read carrying it came from the donor. All "
             "three are weighed against each other read by read, with the base qualities and the fraction of "
             "reads that arrived from the donor fitted rather than assumed. BF is the log10 Bayes factor for a "
             "conversion over the best of the other two, and above 3 is decisive; the fraction next to it is how "
             "much of the locus the model had to write off as donor reads. A tract covering every diagnostic "
             "site of its locus stays 'ambiguous' by arithmetic rather than by rule: it predicts exactly the "
             "same bases as every read having come from the donor. The allele fractions in the plot and the "
             "table are descriptive, and they are what you can go and check in the BAM. These regions are "
             "excluded from variant calling and the consensus by design, so a tract will not appear in the SNP "
             "matrix. Candidates to inspect, NOT confirmed recombination events, and breakpoints are located "
             "only to diagnostic-site resolution.",
    "flagged": "The samples the current thresholds flag as WARN / FAIL and the specific reason for each. This "
               "is the actionable QC summary; adjust the thresholds above to re-flag the whole report.",
    "dosetx": "The per-sample dose (a numeric samplesheet column) split by treatment group, drawn as a box (IQR "
              "+ median) with the individual samples, plus a Kruskal-Wallis rank test of whether dose differs "
              "across the treatment groups (a Mann-Whitney-equivalent when there are two groups). The test runs "
              "over the whole cohort, independent of the live filters; groups with fewer than two dosed samples "
              "are drawn but not tested. Click a point to highlight that sample across the report. Descriptive, "
              "not a claim about efficacy.",
    "vardose": "An association scan: for every variant site, the per-sample allele frequency (0 where the site is "
               "reference) is rank-correlated with dose across the dosed samples (Spearman rho + two-sided p), and "
               "a Benjamini-Hochberg FDR q is computed across all tested variants so the multiple testing is "
               "controlled. The table ranks the variants; click a row to plot its allele-frequency-vs-dose scatter, "
               "click a column header to re-sort, and click a point to highlight that sample. Sites with fewer than "
               "three carriers or no allele-frequency variation are skipped. A screen for dose-associated variants, "
               "NOT proof of causation.",
}
