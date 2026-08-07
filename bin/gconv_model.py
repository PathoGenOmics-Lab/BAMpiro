#!/usr/bin/env python3
"""A changepoint model for gene conversion, evaluated over reads instead of over sites.

The first version of this detector was a cascade of thresholds: call a site converted if the
donor allele fraction clears 0.7, join adjacent calls into runs, keep runs of three or more,
then reject the run if the average donor fraction outside it clears 0.25. Every one of those
numbers is a decision the data were never asked about, and the cascade throws away three things
that matter:

* **Depth.** 7 donor reads out of 10 and 700 out of 1000 both give a fraction of 0.7. They are
  not remotely the same evidence.
* **Base quality.** A Q40 base and a Q14 base are worth the same to a hard floor, and a base
  under the floor is worth nothing at all rather than a little.
* **Read linkage.** A pileup treats each site as an independent observation, but at 100 bp reads
  and diagnostic sites 30 bp apart it is the SAME molecules voting at both. The effective sample
  size is the number of molecules, not the sum of the depths.

This module replaces the cascade with an explicit model of how the data came to look that way.

THE MODEL
---------

A locus has S diagnostic sites. Each read is either **native** to the acceptor, in which case it
carries whatever the acceptor genome carries, or it **mismapped** from the donor, in which case
it carries the donor base everywhere it reaches. The fraction that mismapped, `m`, is a property
of the locus and of the aligner, not of the hypothesis being tested, so it is a nuisance
parameter and it is marginalised over rather than thresholded on.

The acceptor genome is described by a single interval [i, j] of diagnostic sites, the conversion
tract, which carries donor bases; everything outside carries acceptor bases. The null is the
empty interval.

For one read r and one site s, with observed base b and error probability e from the Phred score,
what separates the two hypotheses is a log-likelihood ratio:

    delta[r, s] = log P(b | donor base) - log P(b | acceptor base)
                = +log(3(1-e)/e)   if b is the donor base
                  -log(3(1-e)/e)   if b is the acceptor base
                   0               otherwise

A third base is exactly as unlikely under both, so it contributes nothing. That is the correct
statement, and it is one the allele-fraction version could not make: there, a third base counted
towards depth and vanished from the numerator, which quietly moved the fraction.

The likelihood of the locus is then a product over READS of a two-component mixture:

    L = prod_r [ m * L_donor(r) + (1-m) * L_native(r | i, j) ]

and this is the point of the whole exercise. A read that spans a breakpoint carries donor bases
on one side and acceptor bases on the other, so its native likelihood is large under the right
[i, j] and tiny under any other, while its donor-only likelihood is tiny everywhere. Read-level
evidence is not a separate counter bolted onto the side of the score: it IS the score. And
because the product is over reads rather than over sites, a molecule covering four sites is one
observation with four correlated parts, not four independent votes.

WHY THE SEARCH IS EXACT AND NOT A HEURISTIC
-------------------------------------------

Writing the native log-likelihood with the acceptor as the baseline,

    log L_native(r | i, j) = base_r + sum_{s in [i,j]} delta[r, s]

the baseline term is shared by every hypothesis and cancels, and the sum over the interval is a
difference of two prefix sums. So the likelihood of ANY interval costs one subtraction per read.
Every one of the S(S+1)/2 intervals is evaluated exactly, on every grid point of `m`. There is no
optimiser, no local maximum to get stuck in, and no segmentation heuristic to defend.

WHAT COMES OUT
--------------

Marginalising over the intervals with a prior on tract length, and over `m`, gives a Bayes factor
for "a tract exists" against "no tract", which is the number to read. Normalising over intervals
gives a posterior over breakpoints, so the boundaries come with a credible interval instead of
being wherever a threshold happened to be crossed.

Two things fall out of the model that previously needed special cases:

* **Wholesale mismapping** is no longer a competing verdict. It is `m` being large, which raises
  the null's likelihood everywhere at once and leaves the tract nothing to explain.
* **A tract covering every site** is genuinely indistinguishable from `m = 1`: both say every read
  carries donor bases. The two hypotheses have the same likelihood, so the Bayes factor collapses
  to the prior ratio on its own. The model reports the ambiguity instead of having to be told
  about it.

WHAT A RUN OVER A REAL GENOME ADDED
-----------------------------------

Two components exist because a genome-wide run on H37Rv, with its 411 real paralog pairs and an
isolate simulated with real error profiles, indels, duplicates and contamination, showed the model
without them getting the wrong answer for a reason no toy dataset would have shown.

**The tract fraction.** Not every read that should carry a tract does. A mixed infection is the
obvious reason and the boring one; the interesting one is that a gene family usually has more than
two members. 21% of H37Rv's paralog pairs have a relative CLOSER than their own donor, and that
relative's reads land on the acceptor carrying unconverted bases. Scoring them as evidence against
the tract is how a perfectly clonal conversion comes out at an allele fraction of 0.5 and gets
thrown away, which is exactly what happened. The fraction is now fitted and reported, with a prior
that keeps a stray blip from becoming "a tract at some convenient minority frequency", and a floor
under it because below about a fifth of the reads a minority conversion and a contaminating sample
are the same picture.

**The substitution rate, measured at the locus.** The rate a run of donor-matching sites has to
beat is not a property of the genome, it is a property of the gene. PE/PPE genes are hypervariable
and they are hypervariable in exactly the positions where the copies already differ, so a locus
that carries a scattering of donor-matching bases on its own is a locus where a run of two or three
is unremarkable. Every diagnostic site outside the candidate tract is a direct observation of that
rate, so it is estimated rather than assumed. Measured on a simulated hypervariable isolate: 39% of
those loci were called conversions with the genome-average rate, 8% with the locus's own, and the
same 10 of 10 real tracts were found either way.
"""

from __future__ import annotations

import math

import numpy as np

LN10 = math.log(10.0)

# A discrete prior over the fraction of reads at a locus that arrived from its paralog. Log
# spaced because the interesting range spans orders of magnitude: a clean locus sits near 0.001
# and a badly collapsed repeat near 0.5. Uniform weight over the grid points, so no value of the
# mismapping rate is assumed.
MISMAP_GRID = (0.0, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05,
               0.1, 0.15, 0.2, 0.3, 0.5, 0.7, 0.9, 0.95, 0.98, 0.995)

# Fraction of the reads that are native to the acceptor AND carry the tract. A clonal conversion
# in a two-copy family sits at 1. Two ordinary situations pull it down, and from allele data they
# are the same situation:
#
#   * a mixed infection, where only part of the population carries the conversion;
#   * a gene family with more than two members, where the copies that were not converted are
#     similar enough that their reads land on the acceptor as well. In H37Rv this is not an edge
#     case: 21% of paralog pairs have a relative CLOSER than their own donor.
#
# Without this the model scores every one of those reads as evidence against the tract, which is
# how a perfectly clonal conversion in a three-copy family comes out at an allele fraction of 0.5
# and gets thrown away.
# The grid stops at 0.2 on purpose. Below about a fifth of the reads there is nothing left to
# tell a subclonal conversion from index hopping, cross-sample contamination or an aligner having
# a bad day, and a floor of 0.1 was measured to turn a 4% blip into a called conversion on the
# real genome. What sits under the floor belongs in donor_af_in for a human to look at, not in a
# verdict.
TRACT_AF_GRID = (1.0, 0.95, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2)

# How much a diluted tract has to pay for it. An exponential prior on the shortfall from 1, so a
# clonal tract is free, a half-strength one costs about an order of magnitude, and one at the
# floor costs nearly two.
TRACT_AF_SCALE = 0.18

# Probability that a given diagnostic site independently carries the donor base in this sample
# because it mutated to it, rather than because it was converted. Roughly the per-site SNP rate
# of a clonal genome against its reference (about 1e-3 for a tuberculosis isolate) divided by the
# three bases it could have mutated to. This is what stops the model calling every SNP that
# happens to match the paralog a one-site conversion tract.
MUT_RATE = 3e-4

# Ceiling on the rate estimated from a locus, and the log-likelihood ratio at which a site counts
# as carrying the donor base at all. Without the ceiling a locus dense in substitutions could
# explain away a long, clean tract as a run of coincidences.
MAX_MUT_RATE = 0.25
SITE_CALL_LR = 10.0

# Phred scores from a real sequencer are optimistic in exactly the situation this tool works in:
# a paralogous region, where a "sequencing error" is often a real base from the other copy. Left
# uncapped, a Q40 base claims one error in ten thousand and a handful of them can outvote the
# rest of the pileup. Capping the score bounds how much any single observation can assert.
MAX_PHRED = 30

# Cap on the (candidates x reads) working array, in elements. The exact search is quadratic in
# the number of diagnostic sites, and a 5 kb repeat at 95% identity has a few hundred of them, so
# the array is chunked rather than allocated whole.
CHUNK_ELEMENTS = 4_000_000

# Diagnostic sites above which breakpoints are placed on a coarser grid. The exhaustive search
# costs O(sites^2 x reads), so a 400-site locus takes a minute and a 1000-site one is out of the
# question. Above this, endpoints are spaced out and breakpoint resolution drops; nothing else
# about the model changes. Set high enough that ordinary paralog pairs are searched exactly.
MAX_GRID = 150


def _logsumexp(x):
    """log(sum(exp(x))) without overflowing. Local so this module needs no scipy."""
    x = np.asarray(x, dtype=float)
    if x.size == 0:
        return -np.inf
    top = np.max(x)
    if not np.isfinite(top):
        return float(top)
    return float(top + np.log(np.sum(np.exp(x - top))))


def _logsumexp_axis(x):
    """_logsumexp down the rows of a 2-D array, one result per column."""
    top = np.max(x, axis=0)
    out = np.array(top, dtype=float)
    finite = np.isfinite(top)
    if finite.any():
        out[finite] = top[finite] + np.log(np.sum(np.exp(x[:, finite] - top[finite]), axis=0))
    return out


def delta_matrix(observations, positions, sites, min_bq=13, max_phred=MAX_PHRED):
    """Per-read, per-site log-likelihood ratio of the donor allele over the acceptor allele.

    `observations` is {read_name: {ref_pos: (base, phred)}}, keyed by read NAME so the two mates
    of a pair are one molecule. `sites` maps a position to (acceptor_base, donor_base).

    Returns (read_names, delta) with delta of shape (n_reads, n_sites). Reads carrying no
    informative base are dropped: they contribute an identical constant to every hypothesis.
    """
    index = {p: s for s, p in enumerate(positions)}
    names, rows = [], []
    for name, calls in observations.items():
        row = np.zeros(len(positions))
        informative = False
        for p, (base, phred) in calls.items():
            s = index.get(p)
            if s is None or phred < min_bq:
                continue
            acc, don = sites[p]
            if base == acc:
                sign = -1.0
            elif base == don:
                sign = 1.0
            else:
                continue                    # equally unlikely under both: it separates nothing
            e = 10.0 ** (-min(phred, max_phred) / 10.0)
            row[s] = sign * (math.log(1.0 - e) - math.log(e / 3.0))
            informative = True
        if informative:
            names.append(name)
            rows.append(row)
    delta = np.vstack(rows) if rows else np.zeros((0, len(positions)))
    return names, delta


def _free_runs(free):
    """Maximal runs of sites not already claimed by a tract, as (first, last) index pairs.

    A second tract is searched for AROUND the first rather than through it. Over an
    already-converted site the two hypotheses agree, so an interval crossing one would be tied
    with the interval that stops short of it, and the model would have no reason to prefer either.
    """
    runs = []
    n = len(free)
    s = 0
    while s < n:
        if not free[s]:
            s += 1
            continue
        e = s
        while e + 1 < n and free[e + 1]:
            e += 1
        runs.append((s, e))
        s = e + 1
    return runs


def _candidates(positions, free, max_span_bp=None, stride=1, fine_len=8):
    """Candidate tracts, as parallel (start index, end index) arrays.

    At stride 1 this is every interval of consecutive free sites, which is the exact search. A
    locus with several hundred diagnostic sites has tens of thousands of intervals and the cost is
    quadratic, so above a size limit the endpoints are placed every `stride` sites instead. That
    costs breakpoint resolution, and only that: the likelihood of each interval evaluated is still
    exact.

    Short intervals are always kept at full resolution regardless of stride. A two-site tract is
    exactly the case a coarse grid would mangle, by offering only intervals several sites longer
    than the tract really is, and short tracts are the ones whose evidence is thinnest to begin
    with.
    """
    pairs = set()
    for a, b in _free_runs(free):
        for i in range(a, b + 1):
            for j in range(i, min(b, i + fine_len - 1) + 1):
                pairs.add((i, j))
        if stride > 1:
            grid = list(range(a, b + 1, stride))
            if grid[-1] != b:
                grid.append(b)
            for gi in grid:
                for gj in grid:
                    if gj >= gi:
                        pairs.add((gi, gj))
        else:
            for i in range(a, b + 1):
                for j in range(i, b + 1):
                    pairs.add((i, j))
    if max_span_bp:
        pairs = {(i, j) for i, j in pairs
                 if positions[j] - positions[i] + 1 <= max_span_bp}
    ordered = sorted(pairs)
    return (np.array([i for i, _ in ordered], dtype=int),
            np.array([j for _, j in ordered], dtype=int))


def _credible_interval(weights, mass=0.95):
    """Smallest index range holding `mass` of the posterior, as (lo, hi) indices.

    Highest-density first: take candidates in order of decreasing weight until enough mass is
    covered, then report the span of what was taken. A sharp posterior gives a single index and a
    diffuse one gives a wide range, which is the honest statement about where a breakpoint is.
    """
    order = np.argsort(weights)[::-1]
    taken, acc = [], 0.0
    for k in order:
        taken.append(int(k))
        acc += float(weights[k])
        if acc >= mass:
            break
    return min(taken), max(taken)


def fit_locus(delta, positions, free=None, prior=0.01, mean_span_bp=1000.0,
              m_grid=MISMAP_GRID, max_span_bp=None, mut_rate=MUT_RATE, max_grid=MAX_GRID,
              af_grid=TRACT_AF_GRID):
    """Evaluate every conversion tract against everything else that could produce donor bases.

    There are three ways an acceptor site can show the donor's base, and all three are in the
    comparison:

    * it was converted, together with its neighbours, as one tract;
    * it mutated to that base on its own, which for a single site is a perfectly ordinary thing
      to happen and needs no conversion at all;
    * the read carrying it came from the donor.

    The second is what sets the size a tract has to reach before it means anything. For the SAME
    set of sites the conversion and mutation hypotheses have identical likelihoods, so the
    comparison between them is a pure ratio of priors: one conversion event against k independent
    substitutions each of probability `mut_rate`. At one site the mutation wins, at two it is
    close, and by three the conversion wins by orders of magnitude. That is where the old
    `min_sites = 3` came from, except now it is derived rather than chosen, and it moves on its
    own with the density of diagnostic sites.

    The mutation hypothesis is enumerated over the same intervals as the conversion hypothesis
    rather than over arbitrary subsets. That is exact where it matters, at the one and two site
    sets that a conversion has to beat, and understates it only for long runs, where the prior
    ratio is astronomical in the conversion's favour anyway.

    Returns None when there is nothing to evaluate (no informative read, or no free site). The
    caller decides what to report; this function never applies a threshold.
    """
    n_sites = len(positions)
    if free is None:
        free = np.ones(n_sites, dtype=bool)
    n_free = int(free.sum())
    stride = max(1, -(-n_free // max_grid)) if max_grid else 1
    starts, ends = _candidates(positions, free, max_span_bp, stride)
    n_reads = delta.shape[0]
    if n_reads == 0 or starts.size == 0:
        return None

    pos = np.asarray(positions, dtype=float)
    # prefix[:, k] is the sum of delta over sites 0..k-1, so any interval is one subtraction.
    prefix = np.concatenate([np.zeros((n_reads, 1)), np.cumsum(delta, axis=1)], axis=1)
    # A mismapped read carries the donor base at EVERY site, which is the sum of the whole row.
    total = prefix[:, -1]

    # Prior over how much of the acceptor's own read pool carries the tract. Weighted, unlike the
    # mismapping grid: a clonal tract is the expected thing and a diluted one has to earn it.
    af = np.asarray(af_grid, dtype=float)
    log_w_af = -(1.0 - af) / TRACT_AF_SCALE
    log_w_af -= _logsumexp(log_w_af)
    log_w_m = -math.log(len(m_grid))
    chunk = max(1, CHUNK_ELEMENTS // max(1, n_reads))
    null_by_m = np.empty(len(m_grid))
    cand_by_m = np.empty((len(m_grid), starts.size))
    # Posterior weight of each tract fraction, accumulated across candidates, so the reported
    # fraction is the one the winning tract needed rather than a global average.
    cand_by_af = np.full((len(af_grid), starts.size), -np.inf)
    for k, m in enumerate(m_grid):
        log_m = math.log(m) if m > 0 else -np.inf
        log_1m = math.log1p(-m)
        # The null: no tract, so a native read follows the all-acceptor background and its
        # native term is the baseline, which is zero on this scale. The tract fraction does not
        # enter: with no tract to carry, carrying it at any frequency is the same thing.
        null_by_m[k] = float(np.logaddexp(log_m + total, log_1m).sum())
        for a in range(0, starts.size, chunk):
            b = min(a + chunk, starts.size)
            gain = prefix[:, ends[a:b] + 1] - prefix[:, starts[a:b]]
            acc = np.full(b - a, -np.inf)
            for q, p in enumerate(af):
                # A native read either carries the tract or does not; a read that does not looks
                # exactly like the unconverted acceptor, which is the baseline.
                inner = (np.logaddexp(math.log(p) + gain, math.log1p(-p)) if p < 1.0 else gain)
                ll = np.logaddexp(log_m + total[:, None], log_1m + inner).sum(axis=0)
                acc = np.logaddexp(acc, ll + log_w_af[q])
                cand_by_af[q, a:b] = np.logaddexp(cand_by_af[q, a:b], ll + log_w_m)
            cand_by_m[k, a:b] = acc

    # Marginal over the mismapping rate.
    null_ll = _logsumexp(null_by_m + log_w_m)
    cand_ll = _logsumexp_axis(cand_by_m + log_w_m)

    # Prior over tracts: exponential in span, so the model prefers a short tract to a long one at
    # equal likelihood. Bacterial conversion tracts run from tens of bases to a few kb.
    span = pos[ends] - pos[starts] + 1.0
    log_prior_span = -span / float(mean_span_bp)
    log_prior_span -= _logsumexp(log_prior_span)

    joint = cand_ll + log_prior_span
    log_ev_conv = _logsumexp(joint)
    # Posterior over tracts, given that one exists. Only the span prior enters, so this is
    # settled before the substitution rate is, and the tract it picks is what that rate is then
    # estimated around.
    weights = np.exp(joint - log_ev_conv)
    best = int(np.argmax(weights))

    # How often a site at THIS locus carries the donor base on its own, measured here rather than
    # assumed. A locus is not a random stretch of genome: PE/PPE genes are hypervariable, and
    # they are hypervariable in exactly the positions where the copies already differ. A rate
    # taken from the genome average says two adjacent donor-matching sites are a one-in-ten-
    # million coincidence, and a genome-wide run showed that calling 22 of 60 isolated
    # substitutions a conversion. Every diagnostic site outside the candidate tract is a direct
    # observation of that rate, so the locus is asked instead of told.
    #
    # Only ISOLATED substituted sites count. A run of them is what a conversion looks like, so
    # counting runs would let a second genuine tract in the same locus inflate the rate and talk
    # the first one down, and two conversions in one paralog pair is an ordinary outcome.
    site_lr = delta.sum(axis=0)
    outside = np.ones(n_sites, dtype=bool)
    outside[starts[best]:ends[best] + 1] = False
    measured = outside & (np.abs(site_lr) > SITE_CALL_LR)
    carries = measured & (site_lr > 0)
    lone = carries.copy()
    lone[:-1] &= ~carries[1:]
    lone[1:] &= ~carries[:-1]
    substituted = int(lone.sum())
    local_rate = substituted / int(measured.sum()) if measured.any() else 0.0
    # The ceiling bounds the ESTIMATE, not the caller: a locus that happens to be noisy should
    # not be able to explain away a long clean tract, but an explicit rate is an instruction.
    mut_rate = max(mut_rate, min(MAX_MUT_RATE, local_rate))

    # Prior over the same site sets arising as independent substitutions: one factor of mut_rate
    # per site that carries the donor base, and one of (1 - mut_rate) per site that does not.
    n_in = (ends - starts + 1).astype(float)
    log_prior_mut = n_in * math.log(mut_rate) + (n_sites - n_in) * math.log1p(-mut_rate)

    log_ev_mut = _logsumexp(cand_ll + log_prior_mut)
    log_ev_null = null_ll + n_sites * math.log1p(-mut_rate)

    # Priors between the three families. A conversion is the rare event; not-a-conversion carries
    # the rest, and inside it the mutation and no-mutation cases are already weighted by mut_rate.
    log_pi, log_1pi = math.log(prior), math.log1p(-prior)
    log_ev_alt = _logsumexp([log_1pi + log_ev_null, log_1pi + log_ev_mut])
    log_bf = log_ev_conv + log_pi - log_ev_alt
    # The Bayes factor itself, prior on the conversion divided back out, so it stays a statement
    # about the data: how much better a tract explains them than anything else could.
    log10_bf = (log_ev_conv - _logsumexp([log_ev_null, log_ev_mut])) / LN10
    log10_bf_null = (log_ev_conv - log_ev_null) / LN10
    log10_bf_mut = (log_ev_conv - log_ev_mut) / LN10

    post_conv = 1.0 / (1.0 + math.exp(-log_bf)) if log_bf > -700 else 0.0

    start_w = np.zeros(n_sites)
    np.add.at(start_w, starts, weights)
    end_w = np.zeros(n_sites)
    np.add.at(end_w, ends, weights)
    # Posterior that each site is inside the tract: add at the start, subtract past the end, and
    # a running sum turns the intervals into per-site mass in one pass.
    edges = np.zeros(n_sites + 1)
    np.add.at(edges, starts, weights)
    np.add.at(edges, ends + 1, -weights)
    site_post = np.cumsum(edges)[:n_sites]

    s_lo, s_hi = _credible_interval(start_w)
    e_lo, e_hi = _credible_interval(end_w)

    # Posterior over the mismapping rate under the full model, both hypotheses weighted by their
    # priors, so the number reported is what the locus as a whole implies rather than what the
    # winning hypothesis alone would like it to be.
    conv_by_m = np.array([_logsumexp(cand_by_m[k] + log_prior_span) for k in range(len(m_grid))])
    m_ev = np.logaddexp(math.log1p(-prior) + null_by_m,
                        (math.log(prior) if prior > 0 else -np.inf) + conv_by_m)
    m_post = np.exp(m_ev - _logsumexp(m_ev))
    mismap = float(np.dot(m_post, np.array(m_grid, dtype=float)))

    # How much of the acceptor's read pool the MAP tract needed to carry it. Read it as an allele
    # fraction: 1 is a clonal conversion in a two-copy family, and well under 1 means either a
    # mixed infection or a third copy of the family contributing unconverted reads. Nothing in
    # the reads tells those two apart.
    af_ev = cand_by_af[:, best] + log_w_af
    af_post = np.exp(af_ev - _logsumexp(af_ev))
    tract_af = float(np.dot(af_post, af))

    return {
        "map_i": int(starts[best]),
        "map_j": int(ends[best]),
        "start": int(pos[starts[best]]),
        "end": int(pos[ends[best]]),
        "post_tract": float(weights[best]),
        "post_conv": float(post_conv),
        "log10_bf": float(log10_bf),
        "log10_bf_null": float(log10_bf_null),
        "log10_bf_mut": float(log10_bf_mut),
        "mismap_frac": mismap,
        "tract_af": tract_af,
        "start_ci": (int(pos[s_lo]), int(pos[s_hi])),
        "end_ci": (int(pos[e_lo]), int(pos[e_hi])),
        "site_post": site_post,
        "n_free": n_free,
        "n_sites_total": n_sites,
        "n_reads": n_reads,
        "n_candidates": int(starts.size),
        # The substitution rate actually used, after the locus was allowed to raise it.
        "mut_rate": float(mut_rate),
        "n_substituted_outside": substituted,
        # 1 when every interval was evaluated. Above that, breakpoints were placed every `stride`
        # diagnostic sites and the reported boundaries are that much coarser.
        "stride": stride,
    }


def segment(delta, positions, max_tracts=2, min_report_bf=1.0, **kw):
    """Find tracts one at a time, each judged against everything already found.

    The first fit is returned whatever its Bayes factor: a locus whose data are explained by
    mismapping is a result worth reporting, not a silence. Further fits only run while the
    previous one was worth reporting, and each is a test of "another tract ON TOP of the ones
    already accepted", because the accepted sites are folded into the background first.
    """
    free = np.ones(len(positions), dtype=bool)
    work = delta
    out = []
    for _ in range(max(1, max_tracts)):
        res = fit_locus(work, positions, free=free, **kw)
        if res is None:
            break
        out.append(res)
        if res["log10_bf"] < min_report_bf:
            break
        i, j = res["map_i"], res["map_j"]
        # Fold the accepted tract into the background: over its sites the donor base is now what
        # the acceptor genome carries, so a further tract there would change nothing.
        work = work.copy()
        work[:, i:j + 1] = 0.0
        free = free.copy()
        free[i:j + 1] = False
        if not free.any():
            break
    return out


def verdict(res, covers_locus, min_bf=3.0, min_mismap=0.15, min_tract_af=0.25):
    """Turn the model's output into the same vocabulary the threshold version used.

    The reason always names the alternative that came closest, because "not called" is only
    useful when it says what else the data look like.
    """
    bf, bf_null, bf_mut = res["log10_bf"], res["log10_bf_null"], res["log10_bf_mut"]
    # A tract the model can only fit at the bottom of its frequency range is one it cannot speak
    # about. Below roughly a fifth of the reads, a minority conversion, index hopping and a
    # contaminating sample are the same picture, and the fitted fraction is pinned at the floor
    # rather than measured. Reported with its numbers, not called.
    if res.get("tract_af") is not None and res["tract_af"] < min_tract_af:
        return "ambiguous", (
            f"only {res['tract_af']:.0%} of the reads here carry the tract, which is under what "
            "separates a minority conversion from contamination or index hopping")
    if bf >= min_bf:
        return "gene_conversion", (
            f"log10 Bayes factor {bf:.1f} over the best alternative, with the mismapping rate "
            f"marginalised out (posterior mean {res['mismap_frac']:.3f})"
            + (f"; carried by {res['tract_af']:.0%} of the reads"
               if res.get("tract_af", 1.0) < 0.9 else ""))

    if covers_locus:
        return "ambiguous", (
            "the tract covers every diagnostic site, which has the same likelihood as every read "
            f"having arrived from the donor: log10 Bayes factor only {bf:.1f}")
    if bf_mut <= bf_null:
        n = res["map_j"] - res["map_i"] + 1
        return "ambiguous", (
            f"{n} site(s) carrying the donor base is better explained by independent "
            f"substitution than by a conversion: log10 Bayes factor {bf_mut:.1f}")
    if res["mismap_frac"] >= min_mismap:
        return "mismapping", (
            f"the locus is explained by reads arriving from the donor at rate "
            f"{res['mismap_frac']:.2f}; log10 Bayes factor for a tract is only {bf:.1f}")
    return "ambiguous", (
        f"log10 Bayes factor {bf:.1f} is below the {min_bf:.1f} needed to call a conversion")
