"""Two guards against the same failure: a sample mapped to a genome it does not belong to.

In a 185-sample cohort, 22 samples annotated as one lineage were really another. They were routed
to that lineage's reference, carried about 2,000 SNPs where their line mates carried 5, and
produced 507 of the run's 522 gene-conversion tracts: at every paralogous locus a diverged
genotype carries the donor's base by inheritance, which is exactly what the model scores as a
converted tract. Tracts per sample against genome-wide SNP count correlated at r = 0.98.

Nothing in the run failed. The reference was a valid genome, the reads mapped to it, the metrics
were computed correctly. They were simply answering a question nobody asked.

The first guard is at the report: k-mer lineage typing does not use the reference, so it is
independent evidence, and a sample disagreeing with the rest of its reference's samples is
flagged. The second is in the cohort step, which refuses to read a sample's tracts as conversion
when there are too many of them to be local events.
"""

from __future__ import annotations

import sys

import pytest

from conftest import REPO_ROOT

sys.path.insert(0, str(REPO_ROOT / "bin"))
from qcreport.metrics import flag_sample                 # noqa: E402
import gconv_cohort as gc                                # noqa: E402

THR = {"depth_min": 10, "breadth_min": 90, "missing_max": 10, "mapping_min": 80, "dup_max": 50,
       "iupac_max": 1, "titv_min": 0, "snp_z": 3, "het_max_frac": 1, "mixed_min_frac": 10}
CLEAN = {"mean_depth": 50, "breadth_pct": 99.9, "missing_pct": 1, "mapped_pct": 99}


def flags(lineage, ref_lineage):
    return flag_sample({**CLEAN, "lineage": lineage}, THR, None, None, ref_lineage=ref_lineage)[1]


def test_a_sample_typed_as_another_lineage_is_flagged():
    assert "LINEAGE_MISMATCH" in flags("A4;M_bovis;A4.7", "L7")


def test_agreement_is_not_flagged():
    assert "LINEAGE_MISMATCH" not in flags("L7", "L7")


def test_sublineages_of_the_same_lineage_agree():
    """The comparison is on the top level. A4.7.1 is still A4 and must not be flagged."""
    assert "LINEAGE_MISMATCH" not in flags("A4;M_bovis;A4.7;A4.7.1", "A4")


def test_nothing_is_flagged_without_a_reference_majority():
    """A cohort too small to have a majority must not invent one."""
    assert "LINEAGE_MISMATCH" not in flags("A4", None)


def test_an_untyped_sample_is_not_accused():
    """No lineage call is missing evidence, not conflicting evidence."""
    assert "LINEAGE_MISMATCH" not in flags("", "L7")


# --------------------------------------------------------------------------- the cohort guard

def rows(spec, n_loci=100):
    """spec is {sample: number of loci it calls a tract in}, over a cohort of n_loci loci."""
    out = []
    for pair in range(n_loci):
        out.append({"sample": "filler", "pair_id": str(pair), "verdict": "ambiguous"})
    for sample, k in spec.items():
        for pair in range(k):
            out.append({"sample": sample, "pair_id": str(pair), "verdict": "gene_conversion"})
    return out


def test_a_sample_converting_a_third_of_its_loci_is_divergent():
    div = gc.divergent_samples(rows({"bad": 34, "good": 3}), 0.1)

    assert set(div) == {"bad"}
    assert div["bad"] == pytest.approx(0.34)


def test_the_measured_separation_is_respected():
    """The real cohort put mislabelled samples at 0.20 to 0.39 and every correct one at most at
    0.082. The default threshold has to sit in that gap, and it has to sit there with room."""
    div = gc.divergent_samples(rows({"mislabelled": 20, "correct": 8}), 0.1)

    assert set(div) == {"mislabelled"}


def test_one_locus_scored_against_several_donors_counts_once():
    """Counting rows rather than loci would push a clean sample over the line on duplication
    alone: a single tract is emitted once per candidate donor."""
    dup = [{"sample": "s", "pair_id": "7", "verdict": "gene_conversion"} for _ in range(40)]
    dup += [{"sample": "filler", "pair_id": str(i), "verdict": "ambiguous"} for i in range(100)]

    assert gc.divergent_samples(dup, 0.1) == {}


def test_the_divergent_verdict_beats_corroboration():
    """The rule has to come first. Diverged samples share a genotype, so they share their
    artefacts at identical coordinates, which is exactly what the corroboration rule rewards."""
    row = {"sample": "bad", "verdict": "ambiguous", "log10_bf": "5.0", "reason": ""}
    event = {"n_samples": 20, "n_called": 19}

    verdict, reason = gc.cohort_verdict(row, event, {}, 185, 3.0, 1.0, 0.9, 5,
                                        divergent={"bad": 0.33})

    assert verdict == "divergent_sample"
    assert "33%" in reason


def test_a_clean_sample_is_still_corroborated():
    row = {"sample": "ok", "verdict": "ambiguous", "log10_bf": "5.0", "reason": ""}
    event = {"n_samples": 20, "n_called": 19}

    verdict, _ = gc.cohort_verdict(row, event, {}, 185, 3.0, 1.0, 0.9, 5, divergent={"bad": 0.33})

    assert verdict == "gene_conversion"


def test_a_handful_of_loci_cannot_support_a_fraction():
    """One tract out of three loci is 33 percent and means only that the reference has three
    paralogous loci. Firing there broke 61 existing tests before the floors went in."""
    assert gc.divergent_samples(rows({"s": 1}, n_loci=3), 0.1) == {}


def test_a_few_tracts_in_a_big_locus_set_are_not_divergence():
    """25 loci, 3 tracts: 12 percent clears the ratio but three local events are just three
    local events."""
    assert gc.divergent_samples(rows({"s": 3}, n_loci=25), 0.1) == {}


def test_the_floors_do_not_rescue_a_really_divergent_sample():
    assert set(gc.divergent_samples(rows({"s": 30}, n_loci=100), 0.1)) == {"s"}


# --------------------------------------------------------------------------- the wiring
#
# The flag above shipped working and never fired. flag_sample was right; qc_report fed it a
# majority computed from summ[sid]["reference"], a key the summary has never had, so the majority
# was always empty. On the cohort it was written for it flagged nothing while 18 samples sat on
# the wrong reference. Every test above passes ref_lineage in by hand, which is exactly why none
# of them saw it. These read the reference the way qc_report does.

import qc_report as qr                                   # noqa: E402


def test_the_reference_comes_from_the_samplesheet(tmp_path):
    ss = tmp_path / "ss.tsv"
    ss.write_text("sampleId\trunId\tr1\tr2\trefId\n"
                  "S1\tR1\ta\tb\tREF_A\n"
                  "S1\tR2\tc\td\tREF_A\n"
                  "S2\tR1\te\tf\tREF_B\n")

    assert qr.sample_references(str(ss)) == {"S1": "REF_A", "S2": "REF_B"}


def test_no_samplesheet_is_no_reference_rather_than_an_error():
    assert qr.sample_references(None) == {}


def test_the_summary_has_no_reference_column_to_rely_on():
    """The regression itself: if a future summary gains a 'reference' column this is harmless, but
    the flag must not DEPEND on it, because today it does not exist."""
    src = (REPO_ROOT / "bin" / "qc_report.py").read_text()
    assert 'm.get("reference")' not in src, "the majority is being read from the summary again"
    assert "sample_references(" in src
