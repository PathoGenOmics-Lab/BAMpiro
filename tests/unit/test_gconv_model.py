"""Unit tests for bin/gconv_model.py (the gene conversion changepoint model).

These are the specification of what replaced the threshold cascade, and most of them are
statements the old version could not have made at all. Each one is written as a comparison
between two datasets rather than as a magic number, because what matters is not that a Bayes
factor is 8.4, it is that the evidence moves in the right direction when the data change:

* more depth at the same allele fraction is more evidence;
* a high-quality base is worth more than a low-quality one;
* one molecule carrying a donor base and an acceptor base is worth more than two molecules each
  carrying one of them, which is the whole reason the likelihood is a product over READS;
* one site carrying the donor base is better explained by an ordinary substitution, and three
  sites are not, which is where the old hand-picked `min_sites = 3` comes from.

The arithmetic that is checked against exact numbers is only the part that has a closed form:
the per-observation log-likelihood ratio, and the candidate enumeration.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from conftest import load_script

gm = load_script("gconv_model")


# --------------------------------------------------------------------------- helpers


def _positions(n, step=40, first=1000):
    return [first + step * i for i in range(n)]


def _sites(positions, acc="A", don="C"):
    return {p: (acc, don) for p in positions}


def _obs(reads, positions, phred=30):
    """{name: {pos: (base, phred)}} from {name: {site index: base}}."""
    return {name: {positions[i]: (base, phred) for i, base in calls.items()}
            for name, calls in reads.items()}


def _tract_reads(n, span, n_sites, tract, prefix="r", phred=30):
    """`n` molecules per starting offset, each covering `span` consecutive sites."""
    reads = {}
    for start in range(n_sites - span + 1):
        for k in range(n):
            calls = {s: ("C" if s in tract else "A") for s in range(start, start + span)}
            reads[f"{prefix}{start}_{k}"] = calls
    return reads


def _fit(reads, positions, phred=30, **kw):
    sites = _sites(positions)
    _, delta = gm.delta_matrix(_obs(reads, positions, phred), positions, sites)
    return gm.fit_locus(delta, positions, **kw)


# ----------------------------------------------------------------- delta_matrix


def test_delta_is_the_log_likelihood_ratio_of_donor_over_acceptor():
    """The one piece with a closed form, checked against it.

    A base equal to the donor's is (1-e) likely under conversion and e/3 under the acceptor, so
    the ratio is log(3(1-e)/e). A base equal to the acceptor's is the same size with the sign
    flipped.
    """
    positions = _positions(2)
    e = 10.0 ** -3.0                                    # Phred 30
    expected = math.log(1 - e) - math.log(e / 3)

    names, delta = gm.delta_matrix(_obs({"donor": {0: "C"}, "acceptor": {1: "A"}}, positions),
                                   positions, _sites(positions))

    assert names == ["donor", "acceptor"]
    assert delta[0, 0] == pytest.approx(expected)
    assert delta[1, 1] == pytest.approx(-expected)


def test_a_third_base_contributes_nothing_in_either_direction():
    """Neither copy expects a G here, so it is exactly as unlikely under both hypotheses.

    The allele-fraction version could not say this: a third base counted towards depth and
    vanished from the numerator, which silently moved the fraction it was excluded from.
    """
    positions = _positions(1)

    _, delta = gm.delta_matrix(_obs({"r": {0: "G"}}, positions), positions, _sites(positions))

    assert delta.shape == (0, 1), "a read with nothing informative is not an observation at all"


def test_a_better_base_quality_asserts_more():
    positions = _positions(1)
    sites = _sites(positions)

    _, low = gm.delta_matrix(_obs({"r": {0: "C"}}, positions, phred=15), positions, sites)
    _, high = gm.delta_matrix(_obs({"r": {0: "C"}}, positions, phred=30), positions, sites)

    assert high[0, 0] > low[0, 0] > 0


def test_the_quality_cap_bounds_what_one_base_can_claim():
    """A Q40 base claims one error in ten thousand. In a paralogous region an apparent error is
    often a real base from the other copy, so the score is capped before it is believed."""
    positions = _positions(1)
    sites = _sites(positions)

    _, capped = gm.delta_matrix(_obs({"r": {0: "C"}}, positions, phred=40), positions, sites)
    _, at_cap = gm.delta_matrix(_obs({"r": {0: "C"}}, positions, phred=gm.MAX_PHRED),
                                positions, sites)

    assert capped[0, 0] == pytest.approx(at_cap[0, 0])


def test_bases_below_the_quality_floor_are_dropped():
    positions = _positions(2)
    obs = {"r": {positions[0]: ("C", 5), positions[1]: ("C", 30)}}

    _, delta = gm.delta_matrix(obs, positions, _sites(positions), min_bq=13)

    assert delta[0, 0] == 0.0 and delta[0, 1] > 0


# ------------------------------------------------------------------ candidates


def test_the_search_is_exhaustive_at_stride_one():
    positions = _positions(5)
    starts, ends = gm._candidates(positions, np.ones(5, bool))

    assert len(starts) == 5 * 6 // 2
    assert set(zip(starts.tolist(), ends.tolist())) == {(i, j) for i in range(5)
                                                        for j in range(i, 5)}


def test_sites_already_claimed_by_a_tract_are_not_searched_through():
    """A second tract is looked for around the first, not through it. Over an already converted
    site the two hypotheses agree, so an interval crossing one would be tied with the interval
    that stops short of it and the model would have no reason to prefer either."""
    free = np.array([True, True, False, True, True])
    starts, ends = gm._candidates(_positions(5), free)

    pairs = set(zip(starts.tolist(), ends.tolist()))
    assert pairs == {(0, 0), (0, 1), (1, 1), (3, 3), (3, 4), (4, 4)}


def test_a_coarse_grid_still_keeps_every_short_interval():
    """Above the size limit the endpoints are spread out, which costs breakpoint resolution.

    Short tracts are exempt: a two-site tract is exactly what a coarse grid would mangle, by
    offering only intervals several sites longer than the tract really is.
    """
    positions = _positions(40)
    starts, ends = gm._candidates(positions, np.ones(40, bool), stride=5, fine_len=3)

    pairs = set(zip(starts.tolist(), ends.tolist()))
    assert len(pairs) < 40 * 41 // 2
    for i in range(40):
        assert (i, i) in pairs
        if i + 2 < 40:
            assert (i, i + 2) in pairs, "short intervals survive the stride"
    assert (0, 35) in pairs and (0, 34) not in pairs, "long intervals land on the grid"


def test_the_last_site_of_a_locus_can_still_end_a_tract_on_a_coarse_grid():
    """The grid steps from the first site and stops before the last unless the last is added.

    A run whose length is not a multiple of the stride leaves its final site off the grid, and
    then the only intervals ending there are the short ones the fine window contributes. A long
    tract running to the end of the aligned stretch would have nothing to be found by, which is
    both the hardest place to place reads and the easiest failure to mistake for a quiet locus.
    """
    n = 302                                     # not a multiple of the stride below
    positions = _positions(n)
    starts, ends = gm._candidates(positions, np.ones(n, bool), stride=3, fine_len=8)

    pairs = set(zip(starts.tolist(), ends.tolist()))
    assert any(j == n - 1 for _, j in pairs), "nothing can end at the last site"
    assert any(j == n - 1 and i < n - 9 for i, j in pairs), \
        "the last site is only reachable from inside the fine window"


def test_the_span_cap_removes_the_tracts_that_are_out_of_scope():
    positions = _positions(6, step=1000)
    starts, ends = gm._candidates(positions, np.ones(6, bool), max_span_bp=2500)

    spans = [positions[j] - positions[i] + 1 for i, j in zip(starts, ends)]
    assert spans and max(spans) <= 2500


def test_the_span_cap_is_a_limit_not_a_boundary_to_stay_under():
    """A tract of exactly the cap is in scope. Off by one here silently narrows the search by a
    site at every locus, and the tracts it drops are the longest ones considered."""
    positions = [0, 500, 999]                   # first to last is exactly 1000 bases
    free = np.ones(3, bool)

    at_cap = set(zip(*(x.tolist() for x in gm._candidates(positions, free, max_span_bp=1000))))
    under = set(zip(*(x.tolist() for x in gm._candidates(positions, free, max_span_bp=999))))

    assert (0, 2) in at_cap, "a tract of exactly the cap is within it"
    assert (0, 2) not in under, "and one base over is not"


# -------------------------------------------------------------------- fit_locus


def test_a_clean_tract_is_located_exactly_and_called():
    positions = _positions(10)
    res = _fit(_tract_reads(6, 3, 10, tract={4, 5, 6}), positions)

    assert (res["start"], res["end"]) == (positions[4], positions[6])
    assert res["log10_bf"] > 3.0
    assert res["post_conv"] > 0.99
    assert res["start_ci"] == (positions[4], positions[4])
    assert res["end_ci"] == (positions[6], positions[6])
    assert gm.verdict(res, covers_locus=False)[0] == "gene_conversion"


def _res(**kw):
    """A model result, for asking `verdict` about one number at a time."""
    out = {"log10_bf": 9.0, "log10_bf_null": 9.0, "log10_bf_mut": 9.0, "tract_af": 1.0,
           "mismap_frac": 0.01, "map_i": 0, "map_j": 3}
    out.update(kw)
    return out


@pytest.mark.parametrize("bf,called", [(3.0, True), (2.9999, False), (3.1, True)])
def test_a_tract_at_exactly_the_calling_threshold_is_called(bf, called):
    """`--min-bf` is the evidence a tract needs, and needing it means reaching it is enough.

    This is the headline setting of the whole stage, and the comparison could be moved by one
    without a single test failing.
    """
    got = gm.verdict(_res(log10_bf=bf, log10_bf_mut=0.0, log10_bf_null=0.0),
                     covers_locus=False, min_bf=3.0)[0]

    assert (got == "gene_conversion") is called


@pytest.mark.parametrize("frac,blamed", [(0.2, True), (0.1999, False), (0.3, True)])
def test_a_locus_at_exactly_the_mismapping_rate_is_reported_as_mismapping(frac, blamed):
    """`--min-mismap` is the rate AT WHICH the locus is blamed on mismapped reads, not past it."""
    got = gm.verdict(_res(log10_bf=0.5, log10_bf_mut=9.0, log10_bf_null=0.0, mismap_frac=frac),
                     covers_locus=False, min_bf=3.0, min_mismap=0.2)[0]

    assert (got == "mismapping") is blamed


@pytest.mark.parametrize("af,called", [(0.25, True), (0.2499, False), (0.26, True)])
def test_a_tract_carried_by_exactly_the_required_fraction_is_called(af, called):
    """`--min-tract-af` is the fraction a tract must reach, so reaching it is enough.

    An off-by-one at a user-facing threshold is invisible: the setting still appears to work,
    every tract at exactly the documented value quietly becomes `ambiguous`, and no failure ever
    points at the comparison. It is worth pinning for that reason rather than for the odds of a
    fraction landing on 0.25 exactly.
    """
    got = gm.verdict(_res(tract_af=af), covers_locus=False, min_tract_af=0.25)[0]

    assert (got == "gene_conversion") is called


def test_a_locus_with_no_conversion_is_not_called():
    positions = _positions(10)
    res = _fit(_tract_reads(6, 3, 10, tract=set()), positions)

    assert res["log10_bf"] < 0
    assert res["post_conv"] < 0.01
    assert gm.verdict(res, covers_locus=False)[0] != "gene_conversion"


def test_the_reported_bayes_factor_is_the_weakest_link():
    """`log10_bf` is a conversion against the best of the other two explanations, so it is
    limited by whichever question is still open.

    Thin data leave the sites themselves in doubt and the comparison against no conversion
    binds. Deep data settle the sites, and what remains is whether a run of donor bases is a
    conversion or a coincidence of substitutions, which no amount of extra depth can answer.
    """
    positions = _positions(10)

    reads = _tract_reads(1, 3, 10, tract={4, 5, 6})
    thin = _fit(reads, positions, phred=13)             # one molecule per offset, at the floor
    deep = _fit(_tract_reads(30, 3, 10, tract={4, 5, 6}), positions)

    assert thin["log10_bf"] == pytest.approx(thin["log10_bf_null"], abs=0.5)
    assert thin["log10_bf"] < thin["log10_bf_mut"]
    assert deep["log10_bf"] == pytest.approx(deep["log10_bf_mut"], abs=1e-6)
    assert deep["log10_bf_null"] > thin["log10_bf_null"] + 100


def _mixed(scale, frac, positions):
    """`frac` molecules in ten carry the donor allele at sites 4 to 6, at any depth."""
    reads = {}
    for start in range(6):
        for k in range(10 * scale):
            donor = k % 10 < frac
            reads[f"r{start}_{k}"] = {s: ("C" if (4 <= s <= 6 and donor) else "A")
                                      for s in range(start, start + 3)}
    return reads


def test_depth_is_evidence_and_an_allele_fraction_is_not():
    """The headline defect of the threshold version: 7 donor reads out of 10 and 700 out of
    1000 both give a fraction of 0.7 and were treated identically.

    It shows up in the comparison against no conversion. The headline number does not move,
    because once the sites are settled the open question is whether a run of donor bases is a
    conversion or a coincidence, and depth has nothing to say about that.
    """
    positions = _positions(8)

    shallow = _fit(_mixed(1, 7, positions), positions)
    deep = _fit(_mixed(10, 7, positions), positions)

    assert deep["log10_bf_null"] > shallow["log10_bf_null"] + 100
    assert deep["log10_bf"] == pytest.approx(shallow["log10_bf"], abs=1e-6)


def test_the_tract_fraction_is_recovered_rather_than_held_against_the_tract():
    """How much of the read pool carries the tract is a quantity, not a disqualification.

    Two ordinary things pull it below 1 and from allele data they are the same thing: a mixed
    infection, and a gene family with more than two members whose unconverted copies contribute
    reads to the acceptor. In H37Rv the second is not rare, so a model that scored those reads
    as evidence against the tract would throw away perfectly clonal conversions.
    """
    positions = _positions(8)

    for frac in (10, 7, 4):
        res = _fit(_mixed(5, frac, positions), positions)
        assert res["tract_af"] == pytest.approx(frac / 10, abs=0.05)
        assert (res["start"], res["end"]) == (positions[4], positions[6])


def test_a_diluted_tract_costs_evidence_but_is_still_reachable():
    """The fraction is not free: it has a prior, and dilution costs against no conversion.

    It costs nothing against independent substitution, and that is correct rather than an
    oversight. A substitution can sit at 30% of the population just as a conversion can, so the
    dilution says nothing about which of the two produced these bases and cancels out of that
    comparison. What it does say is how sure we are the bases are there at all.
    """
    positions = _positions(8)

    clonal = _fit(_mixed(5, 10, positions), positions)
    diluted = _fit(_mixed(5, 3, positions), positions)
    noise = _fit(_mixed(5, 0, positions), positions)

    assert clonal["log10_bf_null"] > diluted["log10_bf_null"] + 50
    assert diluted["log10_bf"] > 3.0
    assert noise["log10_bf"] < 0


def test_a_signal_too_thin_to_be_a_frequency_is_not_given_one():
    """The floor under the tract fraction, which a genome-wide run put there.

    Before it, a handful of donor bases at 4% of the reads on a real paralog pair came out as a
    called conversion at the grid's minimum frequency. Below roughly a fifth of the reads there
    is nothing to separate a minority conversion from index hopping or a contaminating sample,
    so the model is not allowed to reach down there and pick one.
    """
    positions = _positions(8)
    reads = {}
    for start in range(6):
        for k in range(250):
            donor = k % 25 == 0                              # 4% of the molecules
            reads[f"r{start}_{k}"] = {s: ("C" if (4 <= s <= 6 and donor) else "A")
                                      for s in range(start, start + 3)}

    res = _fit(reads, positions)

    assert res["tract_af"] >= min(gm.TRACT_AF_GRID)
    assert gm.verdict(res, covers_locus=False)[0] != "gene_conversion"


def test_one_molecule_spanning_a_breakpoint_beats_two_molecules_that_do_not():
    """Why the likelihood is a product over READS rather than over sites.

    Both datasets give exactly the same pileup at every site. In the first, each molecule carries
    an acceptor base outside the tract AND a donor base inside it, which is something a read that
    arrived from the donor cannot do. In the second the same observations are spread over twice
    as many molecules, and each one on its own is consistent with a donor read. A model working
    off per-site counts cannot tell them apart at all.

    The difference lands in `log10_bf_null`, the comparison against no conversion, where reads
    arriving from the donor are the competing explanation. The headline `log10_bf` cannot show
    it: for a site set the model is already certain about, that number is bounded by how
    implausible independent substitution is, and both datasets name the same two sites. Which is
    itself the right behaviour, and the reason both are reported.
    """
    positions = _positions(6)
    tract = {2, 3}

    linked, split = {}, {}
    for k in range(12):
        linked[f"L{k}"] = {1: "A", 2: "C", 3: "C", 4: "A"}
        split[f"Sout{k}"] = {1: "A", 4: "A"}
        split[f"Sin{k}"] = {2: "C", 3: "C"}
    for k in range(12):                                  # identical flanks in both datasets
        linked[f"F{k}"] = split[f"F{k}"] = {0: "A", 5: "A"}

    a = _fit(linked, positions)
    b = _fit(split, positions)

    assert (a["start"], a["end"]) == (positions[min(tract)], positions[max(tract)])
    assert (b["start"], b["end"]) == (a["start"], a["end"])
    assert a["log10_bf_null"] > b["log10_bf_null"] + 50
    assert a["log10_bf"] == pytest.approx(b["log10_bf"], abs=1e-3)


def test_wholesale_mismapping_is_a_fitted_rate_and_not_a_tract():
    """Reads arriving from the donor carry its bases at EVERY site, including outside any
    candidate tract. That raises the no-conversion likelihood everywhere at once, so there is
    nothing left for a tract to explain, and the rate itself is what gets reported."""
    positions = _positions(10)
    reads = {}
    for start in range(8):
        for k in range(10):
            donor = k < 8                                # 80% of molecules, everywhere
            reads[f"r{start}_{k}"] = {s: ("C" if donor else "A") for s in range(start, start + 3)}

    res = _fit(reads, positions)

    assert res["mismap_frac"] > 0.5
    assert res["log10_bf"] < 3.0
    assert gm.verdict(res, covers_locus=False)[0] == "mismapping"


def test_a_tract_covering_the_whole_locus_is_degenerate_with_a_mismapping_rate_of_one():
    """Every read carrying the donor base everywhere is exactly what both hypotheses predict.

    The threshold version needed a special case for this. Here the two likelihoods are equal by
    construction, so the Bayes factor collapses to the ratio of their priors on its own.
    """
    positions = _positions(6)
    reads = {f"r{start}_{k}": {s: "C" for s in range(start, start + 3)}
             for start in range(4) for k in range(10)}

    res = _fit(reads, positions)

    assert (res["map_i"], res["map_j"]) == (0, len(positions) - 1)
    assert res["mismap_frac"] > 0.9
    assert res["log10_bf"] < 3.0
    verdict, reason = gm.verdict(res, covers_locus=True)
    assert verdict == "ambiguous" and "every diagnostic site" in reason


def test_one_site_is_a_substitution_and_three_sites_are_a_conversion():
    """Where `min_sites` went.

    For the same set of sites, "converted" and "mutated independently" have identical
    likelihoods, so what separates them is only their priors: one conversion event against k
    independent substitutions. At one site the substitution wins and at three it does not, and
    nothing in the code says the number three.
    """
    positions = _positions(10)

    one = _fit(_tract_reads(8, 3, 10, tract={5}), positions)
    three = _fit(_tract_reads(8, 3, 10, tract={4, 5, 6}), positions)

    # Both are certain about which sites carry the donor base: that is not what is in question.
    assert one["log10_bf_null"] > 50 and three["log10_bf_null"] > 50

    assert one["log10_bf"] < 3.0
    assert three["log10_bf"] > 3.0
    verdict, reason = gm.verdict(one, covers_locus=False)
    assert verdict == "ambiguous" and "independent substitution" in reason


def test_the_mutation_rate_moves_the_number_of_sites_a_tract_needs():
    """The bar is derived from a rate, so it can be argued with instead of just overruled."""
    positions = _positions(10)
    reads = _tract_reads(8, 3, 10, tract={5})

    plausible = _fit(reads, positions, mut_rate=3e-4)
    implausible = _fit(reads, positions, mut_rate=1e-9)

    assert plausible["log10_bf"] < 3.0 < implausible["log10_bf"]


def test_an_uncertain_boundary_widens_the_credible_interval():
    """A breakpoint is only as sharp as the molecules that cross it.

    With no read reaching from outside the tract to inside it, the model cannot tell where
    between two diagnostic sites the boundary is, and says so instead of putting it wherever a
    threshold was crossed.
    """
    positions = _positions(8)
    tract = {3, 4}

    crossed = _fit(_tract_reads(8, 4, 8, tract=tract), positions)
    # Only single-site molecules: nothing links a site inside the tract to a site outside it,
    # and the sites either side of the tract are never measured at all.
    isolated = {}
    for s in [3, 4]:
        for k in range(8):
            isolated[f"in{s}_{k}"] = {s: "C"}
    for s in [0, 7]:
        for k in range(8):
            isolated[f"out{s}_{k}"] = {s: "A"}
    loose = _fit(isolated, positions)

    def width(res, key):
        lo, hi = res[key]
        return hi - lo

    assert width(crossed, "start_ci") == 0 and width(crossed, "end_ci") == 0
    assert width(loose, "start_ci") > 0 or width(loose, "end_ci") > 0


def test_a_locus_with_no_informative_read_returns_nothing_rather_than_a_guess():
    positions = _positions(5)

    assert gm.fit_locus(np.zeros((0, 5)), positions) is None


def test_the_grid_is_coarsened_only_when_the_locus_is_too_big_to_search_exactly():
    positions = _positions(40)
    reads = _tract_reads(2, 4, 40, tract={10, 11, 12})

    exact = _fit(reads, positions, max_grid=150)
    coarse = _fit(reads, positions, max_grid=10)

    assert exact["stride"] == 1
    assert coarse["stride"] == 4
    assert coarse["n_candidates"] < exact["n_candidates"]
    # The tract is short, so it survives the stride intact.
    assert (coarse["start"], coarse["end"]) == (exact["start"], exact["end"])


# ------------------------------------------------- the locus substitution rate


def _scatter(positions, tract, lone_sites, n=8, span=3):
    """Reads over a locus with a tract plus isolated substituted sites elsewhere."""
    reads = {}
    carry = set(tract) | set(lone_sites)
    for start in range(len(positions) - span + 1):
        for k in range(n):
            reads[f"r{start}_{k}"] = {s: ("C" if s in carry else "A")
                                      for s in range(start, start + span)}
    return reads


def test_a_locus_that_substitutes_freely_makes_a_short_run_ordinary():
    """The substitution rate is measured at the locus, not taken from the genome average.

    A gene that carries a scattering of donor-matching bases on its own is a gene where a run of
    two or three of them is unremarkable. Assuming the genome average there says such a run is a
    one-in-ten-million coincidence, and a run over the real H37Rv measured that assumption
    calling 39% of hypervariable loci a conversion, against 8% once the locus is asked.
    """
    positions = _positions(24)
    tract = {10, 11}
    # The scattered sites are kept clear of the tract. Put one right beside it and the two
    # intervals tie in likelihood exactly, because the substituted site and the unsubstituted one
    # between them cancel, and the winner is then decided by a prior term worth half a nat. Real
    # data does not tie like that, and a test should not pin an arbitrary tie-break.
    clean = _fit(_scatter(positions, tract, []), positions)
    peppered = _fit(_scatter(positions, tract, [1, 4, 7, 16, 19, 22]), positions)

    assert clean["mut_rate"] == pytest.approx(gm.MUT_RATE)
    assert peppered["mut_rate"] > 0.2
    assert clean["log10_bf"] > 3.0 > peppered["log10_bf"]
    assert (peppered["start"], peppered["end"]) == (positions[10], positions[11])


def test_a_second_tract_is_not_mistaken_for_a_high_substitution_rate():
    """Only ISOLATED substituted sites count towards the rate.

    A run of them is what a conversion looks like, so counting runs would let a second genuine
    tract in the same locus talk the first one down, and two conversions in one paralog pair is
    an ordinary outcome rather than a pathological input.
    """
    positions = _positions(24)
    res = _fit(_scatter(positions, {4, 5, 6}, [15, 16, 17, 18]), positions)

    assert res["mut_rate"] == pytest.approx(gm.MUT_RATE)
    assert res["n_substituted_outside"] == 0


def test_the_locus_rate_stays_between_the_floor_asked_for_and_the_ceiling():
    positions = _positions(24)
    reads = _scatter(positions, {10, 11}, [2, 5, 8, 15, 18, 21])

    assert _fit(reads, positions, mut_rate=0.4)["mut_rate"] == pytest.approx(0.4)
    assert _fit(reads, positions)["mut_rate"] <= gm.MAX_MUT_RATE


# ---------------------------------------------------------------------- segment


def test_two_tracts_in_one_locus_are_both_found():
    """Two conversions in one paralog pair is an ordinary outcome. The threshold version had
    each one counting the other's donor bases as evidence against itself."""
    positions = _positions(14)
    fits = gm.segment(*_delta(_tract_reads(6, 3, 14, tract={2, 3, 4} | {9, 10, 11}), positions),
                      max_tracts=2)

    found = sorted((f["start"], f["end"]) for f in fits if f["log10_bf"] > 3.0)
    assert found == [(positions[2], positions[4]), (positions[9], positions[11])]


def test_the_second_search_is_a_test_of_another_tract_on_top_of_the_first():
    """Once a tract is accepted its sites become the background, so the second fit is not free
    to simply restate the first one in a wider interval."""
    positions = _positions(12)
    fits = gm.segment(*_delta(_tract_reads(6, 3, 12, tract={4, 5, 6}), positions), max_tracts=2)

    assert len(fits) == 2
    assert (fits[0]["start"], fits[0]["end"]) == (positions[4], positions[6])
    assert fits[1]["log10_bf"] < fits[0]["log10_bf"]
    first = set(range(fits[0]["map_i"], fits[0]["map_j"] + 1))
    second = set(range(fits[1]["map_i"], fits[1]["map_j"] + 1))
    assert not (first & second)


def test_segment_stops_once_a_fit_is_not_worth_reporting():
    positions = _positions(10)
    fits = gm.segment(*_delta(_tract_reads(6, 3, 10, tract=set()), positions), max_tracts=4)

    assert len(fits) == 1, "a locus with nothing in it is searched once, not four times"


def test_a_fit_exactly_at_the_reporting_threshold_is_worth_reporting():
    """`--report-bf` is the factor below which a fit stops the search, so landing on it does not.

    The threshold is asked instead of assumed: the tract is fitted once to learn what it scores,
    and the cap is then set to exactly that. Choosing a number and hoping the fit lands on it
    would test the arithmetic of the fixture rather than the comparison.
    """
    positions = _positions(16)
    delta, pos = _delta(_tract_reads(6, 3, 16, tract={4, 5, 6}), positions)
    at = gm.segment(delta, pos, max_tracts=1)[0]["log10_bf"]

    assert len(gm.segment(delta, pos, max_tracts=2, min_report_bf=at)) == 2
    assert len(gm.segment(delta, pos, max_tracts=2, min_report_bf=at + 1e-9)) == 1


def _delta(reads, positions):
    """(delta, positions) ready for segment()."""
    _, delta = gm.delta_matrix(_obs(reads, positions), positions, _sites(positions))
    return delta, positions


# ------------------------------------------------------------- invariants
#
# Things that must hold for any input at all. A numerical model can be wrong in ways no single
# example reveals, and these are the statements that would catch it.


def test_the_order_the_reads_arrive_in_does_not_change_the_answer():
    """The likelihood is a product over molecules, so it is symmetric in them. If this ever
    fails, something is accumulating in a way that depends on iteration order."""
    positions = _positions(12)
    reads = _tract_reads(6, 3, 12, tract={4, 5, 6})
    shuffled = dict(sorted(reads.items(), key=lambda kv: kv[0][::-1]))

    a, b = _fit(reads, positions), _fit(shuffled, positions)

    assert a["log10_bf"] == pytest.approx(b["log10_bf"])
    assert (a["map_i"], a["map_j"]) == (b["map_i"], b["map_j"])


@pytest.mark.parametrize("depths", [(1, 2, 4, 8, 16)])
def test_more_agreeing_molecules_never_weaken_the_case_against_no_conversion(depths):
    positions = _positions(12)
    seen = [_fit(_tract_reads(n, 3, 12, tract={4, 5, 6}), positions)["log10_bf_null"]
            for n in depths]

    assert seen == sorted(seen), f"evidence fell as depth rose: {seen}"


@pytest.mark.parametrize("quals", [(13, 20, 25, 30)])
def test_better_base_quality_never_weakens_it_either(quals):
    positions = _positions(12)
    reads = _tract_reads(6, 3, 12, tract={4, 5, 6})
    seen = [_fit(reads, positions, phred=q)["log10_bf_null"] for q in quals]

    assert seen == sorted(seen), f"evidence fell as quality rose: {seen}"


def test_molecules_that_contradict_the_tract_never_strengthen_it():
    positions = _positions(12)
    reads = _tract_reads(6, 3, 12, tract={4, 5, 6})
    contradicted = dict(reads)
    for k in range(40):
        contradicted[f"against{k}"] = {4: "A", 5: "A", 6: "A"}

    assert _fit(contradicted, positions)["log10_bf_null"] < _fit(reads, positions)["log10_bf_null"]


def test_a_locus_read_end_to_end_gives_the_mirrored_tract():
    """No hidden preference for one side of the locus."""
    positions = _positions(12)
    reads = _tract_reads(6, 3, 12, tract={4, 5, 6})
    mirrored = {n: {11 - i: b for i, b in c.items()} for n, c in reads.items()}

    res = _fit(mirrored, positions)

    assert (res["map_i"], res["map_j"]) == (11 - 6, 11 - 4)


def test_moving_the_whole_locus_along_the_genome_changes_nothing_but_coordinates():
    positions = _positions(12)
    reads = _tract_reads(6, 3, 12, tract={4, 5, 6})
    sites = _sites(positions)
    _, delta = gm.delta_matrix(_obs(reads, positions), positions, sites)

    here = gm.fit_locus(delta, positions)
    far = gm.fit_locus(delta, [p + 1_000_000 for p in positions])

    assert here["log10_bf"] == pytest.approx(far["log10_bf"])
    assert (here["map_i"], here["map_j"]) == (far["map_i"], far["map_j"])


def test_the_prior_does_not_move_the_bayes_factor():
    """The Bayes factor is a statement about the data. Only the posterior may follow the prior,
    and a run over the real genome checks this on every reported row."""
    positions = _positions(12)
    reads = _tract_reads(6, 3, 12, tract={4, 5, 6})

    seen = [_fit(reads, positions, prior=p)["log10_bf"] for p in (0.001, 0.01, 0.5, 0.99)]
    posts = [_fit(reads, positions, prior=p)["post_conv"] for p in (0.001, 0.01, 0.5, 0.99)]

    assert all(v == pytest.approx(seen[0]) for v in seen)
    assert posts == sorted(posts)


# ------------------------------------------------------- degenerate inputs
#
# Each of these was a real failure found by running the model over random inputs. They are
# grouped because they share a shape: a value that is legal to pass, produces no error where it
# is passed, and goes wrong somewhere it cannot be traced back from.


@pytest.mark.parametrize("prior", [0.0, 1.0])
def test_a_prior_of_zero_or_one_is_answered_rather_than_crashed(prior):
    """log(0) is a domain error, and "I do not believe this happens" is a fair thing to say."""
    positions = _positions(8)
    res = _fit(_tract_reads(6, 3, 8, tract={3, 4}), positions, prior=prior)

    assert res is not None
    assert math.isfinite(res["log10_bf"])
    assert res["post_conv"] == (0.0 if prior == 0.0 else 1.0)


@pytest.mark.parametrize("kw", [{"mut_rate": 0.0}, {"mut_rate": 1.0}, {"mut_rate": -0.1},
                                {"prior": -0.1}, {"prior": 1.5}, {"mean_span_bp": 0}])
def test_a_degenerate_parameter_is_refused_where_it_was_passed(kw):
    """Not deep inside the arithmetic, where the message names a local variable instead."""
    positions = _positions(6)
    with pytest.raises(ValueError):
        _fit(_tract_reads(4, 3, 6, tract={2, 3}), positions, **kw)


@pytest.mark.parametrize("positions", [[40, 30, 20, 10], [7, 7, 7, 7], [1, 5, 5, 9]])
def test_positions_out_of_order_are_refused_rather_than_answered(positions):
    """Descending positions give negative spans, which turns the length prior upside down and
    makes the model REWARD long tracts, and two sites at one coordinate is not a thing. Both
    produced a confident, plausible-looking, wrong answer."""
    with pytest.raises(ValueError, match="strictly increasing"):
        gm.fit_locus(np.zeros((4, len(positions))) + 8.0, positions)


@pytest.mark.parametrize("delta,positions", [
    (np.zeros((20, 10)), _positions(10)),                    # no information anywhere
    (np.full((10, 6), -1e3), _positions(6)),                 # every site emphatically acceptor
    (np.full((10, 6), 1e3), _positions(6)),                  # every site emphatically donor
])
def test_the_tract_is_always_inside_its_own_credible_interval(delta, positions):
    """The boundary comes from the joint posterior over intervals and the interval from a
    marginal, so on a flat or bimodal posterior the two disagreed and the tract was reported
    outside the range it was supposed to lie in. Nobody can read that."""
    res = gm.fit_locus(delta, positions)

    assert res["start_ci"][0] <= res["start"] <= res["start_ci"][1]
    assert res["end_ci"][0] <= res["end"] <= res["end_ci"][1]


# ---------------------------------------------------------------------- verdict


def test_every_verdict_says_which_alternative_came_closest():
    """"Not called" is only useful when it names what the data look like instead."""
    base = {"log10_bf": 0.0, "log10_bf_null": 0.0, "log10_bf_mut": 0.0,
            "mismap_frac": 0.0, "map_i": 2, "map_j": 4}

    called = gm.verdict({**base, "log10_bf": 8.0}, covers_locus=False)
    assert called[0] == "gene_conversion" and "Bayes factor 8.0" in called[1]

    mismapped = gm.verdict({**base, "log10_bf": 1.0, "log10_bf_mut": 5.0, "mismap_frac": 0.4},
                           covers_locus=False)
    assert mismapped[0] == "mismapping" and "0.40" in mismapped[1]

    mutated = gm.verdict({**base, "log10_bf": 1.0, "log10_bf_null": 40.0, "log10_bf_mut": 1.0},
                         covers_locus=False)
    assert mutated[0] == "ambiguous" and "independent substitution" in mutated[1]

    unbounded = gm.verdict({**base, "log10_bf": 1.0}, covers_locus=True)
    assert unbounded[0] == "ambiguous" and "every diagnostic site" in unbounded[1]


def test_the_calling_threshold_is_a_parameter_and_not_a_constant():
    res = {"log10_bf": 4.0, "log10_bf_null": 40.0, "log10_bf_mut": 4.0,
           "mismap_frac": 0.0, "map_i": 0, "map_j": 3}

    assert gm.verdict(res, covers_locus=False, min_bf=3.0)[0] == "gene_conversion"
    assert gm.verdict(res, covers_locus=False, min_bf=5.0)[0] == "ambiguous"


# ------------------------------------------------------------ deletion markers


def _del_sites(positions, del_at):
    """Site table where the positions in `del_at` are bases the donor does not have."""
    return {p: (("T", gm.GAP) if i in del_at else ("A", "C"))
            for i, p in enumerate(positions)}


def _del_obs(reads, positions, phred=30):
    return {n: {positions[i]: (b, phred) for i, b in c.items()} for n, c in reads.items()}


def test_a_deletion_site_reads_presence_or_absence_not_which_base():
    """There is no donor base to compare against, so the question changes: a read either carries
    a base here, which is the acceptor's, or spans it with a deletion, which is the donor's."""
    positions = _positions(3)
    sites = _del_sites(positions, {1})
    obs = _del_obs({"present": {1: "T"}, "absent": {1: gm.GAP}, "other": {1: "G"}}, positions)

    names, delta = gm.delta_matrix(obs, positions, sites)

    assert delta[names.index("present"), 1] == pytest.approx(-gm.INDEL_LR)
    assert delta[names.index("absent"), 1] == pytest.approx(gm.INDEL_LR)
    # Any base at all is the acceptor's: it is presence that is being read, not identity.
    assert delta[names.index("other"), 1] == pytest.approx(-gm.INDEL_LR)


def test_a_deletion_carries_a_fixed_weight_rather_than_a_base_quality():
    """A base that is absent has no Phred score to be weighed by, and an aligner places indels
    less reliably than substitutions, so the weight is fixed and deliberately modest."""
    positions = _positions(3)
    sites = _del_sites(positions, {1})

    low = gm.delta_matrix(_del_obs({"r": {1: gm.GAP}}, positions, phred=13), positions, sites)[1]
    high = gm.delta_matrix(_del_obs({"r": {1: gm.GAP}}, positions, phred=41), positions, sites)[1]

    assert low[0, 1] == pytest.approx(high[0, 1]) == pytest.approx(gm.INDEL_LR)


def test_a_tract_carrying_a_deletion_beats_coincidence_by_more():
    """Where the whole value of an indel marker is.

    The headline Bayes factor is capped by how implausible it is that the same markers arose
    independently, so once a tract's sites are settled, more depth cannot move it. Two copies
    losing the same bases at the same place is far longer odds than two copies mutating to the
    same base, so a tract that includes a deletion has more to say against that coincidence.
    """
    positions = _positions(10)
    reads = _tract_reads(8, 3, 10, tract={4, 5, 6})

    all_snp = _fit(reads, positions)
    with_del = gm.fit_locus(
        gm.delta_matrix(_obs(reads, positions), positions, _sites(positions))[1],
        positions, is_del=[i == 5 for i in range(10)])

    assert with_del["log10_bf"] > all_snp["log10_bf"] + 0.5
    assert with_del["log10_bf_mut"] > all_snp["log10_bf_mut"]


def test_the_weight_given_to_a_deletion_marker_is_a_parameter():
    positions = _positions(10)
    reads = _tract_reads(8, 3, 10, tract={4, 5, 6})
    delta = gm.delta_matrix(_obs(reads, positions), positions, _sites(positions))[1]
    is_del = [i == 5 for i in range(10)]

    modest = gm.fit_locus(delta, positions, is_del=is_del, indel_factor=0.5)
    strong = gm.fit_locus(delta, positions, is_del=is_del, indel_factor=0.01)

    assert strong["log10_bf"] > modest["log10_bf"]


def test_a_locus_with_no_deletion_marker_is_unaffected():
    positions = _positions(10)
    reads = _tract_reads(8, 3, 10, tract={4, 5, 6})
    delta = gm.delta_matrix(_obs(reads, positions), positions, _sites(positions))[1]

    without = gm.fit_locus(delta, positions)
    with_flags = gm.fit_locus(delta, positions, is_del=[False] * 10)

    assert without["log10_bf"] == pytest.approx(with_flags["log10_bf"])
