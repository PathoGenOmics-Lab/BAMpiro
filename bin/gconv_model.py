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

# Probability that a given diagnostic site independently carries the donor base in this sample
# because it mutated to it, rather than because it was converted. Roughly the per-site SNP rate
# of a clonal genome against its reference (about 1e-3 for a tuberculosis isolate) divided by the
# three bases it could have mutated to. This is what stops the model calling every SNP that
# happens to match the paralog a one-site conversion tract.
MUT_RATE = 3e-4

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
              m_grid=MISMAP_GRID, max_span_bp=None, mut_rate=MUT_RATE, max_grid=MAX_GRID):
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

    log_w = -math.log(len(m_grid))
    chunk = max(1, CHUNK_ELEMENTS // max(1, n_reads))
    null_by_m = np.empty(len(m_grid))
    cand_by_m = np.empty((len(m_grid), starts.size))
    for k, m in enumerate(m_grid):
        log_m = math.log(m) if m > 0 else -np.inf
        log_1m = math.log1p(-m)
        # The null: no tract, so a native read follows the all-acceptor background and its
        # native term is the baseline, which is zero on this scale.
        null_by_m[k] = float(np.logaddexp(log_m + total, log_1m).sum())
        for a in range(0, starts.size, chunk):
            b = min(a + chunk, starts.size)
            gain = prefix[:, ends[a:b] + 1] - prefix[:, starts[a:b]]
            cand_by_m[k, a:b] = np.logaddexp(log_m + total[:, None], log_1m + gain).sum(axis=0)

    # Marginal over the mismapping rate.
    null_ll = _logsumexp(null_by_m + log_w)
    cand_ll = _logsumexp_axis(cand_by_m + log_w)

    # Prior over tracts: exponential in span, so the model prefers a short tract to a long one at
    # equal likelihood. Bacterial conversion tracts run from tens of bases to a few kb.
    span = pos[ends] - pos[starts] + 1.0
    log_prior_span = -span / float(mean_span_bp)
    log_prior_span -= _logsumexp(log_prior_span)

    # Prior over the same site sets arising as independent substitutions: one factor of mut_rate
    # per site that carries the donor base, and one of (1 - mut_rate) per site that does not.
    n_in = (ends - starts + 1).astype(float)
    log_prior_mut = n_in * math.log(mut_rate) + (n_sites - n_in) * math.log1p(-mut_rate)

    joint = cand_ll + log_prior_span
    log_ev_conv = _logsumexp(joint)
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

    # Posterior over tracts, given that one exists.
    weights = np.exp(joint - log_ev_conv)
    best = int(np.argmax(weights))

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
        "start_ci": (int(pos[s_lo]), int(pos[s_hi])),
        "end_ci": (int(pos[e_lo]), int(pos[e_hi])),
        "site_post": site_post,
        "n_free": n_free,
        "n_sites_total": n_sites,
        "n_reads": n_reads,
        "n_candidates": int(starts.size),
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


def verdict(res, covers_locus, min_bf=3.0, min_mismap=0.15):
    """Turn the model's output into the same vocabulary the threshold version used.

    The reason always names the alternative that came closest, because "not called" is only
    useful when it says what else the data look like.
    """
    bf, bf_null, bf_mut = res["log10_bf"], res["log10_bf_null"], res["log10_bf_mut"]
    if bf >= min_bf:
        return "gene_conversion", (
            f"log10 Bayes factor {bf:.1f} over the best alternative, with the mismapping rate "
            f"marginalised out (posterior mean {res['mismap_frac']:.3f})")

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
