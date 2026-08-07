"""Randomised checks over the gene-conversion code, kept rather than run once and thrown away.

Every defect found in this stage so far came from running something like what is in this file:
a prior of 0 taking log(0), positions out of order producing a confident wrong answer, a rate of
0 or 1 dying deep in the arithmetic, a tract reported outside its own credible interval. Each was
turned into a fixed-case regression test, and the GENERATOR that found it was thrown away, which
means the next defect of the same shape waits for someone to go looking again.

So the generators live here. They are seeded, so a failure names the seed and is reproducible,
and they explore, so they are not limited to the cases already known to be broken. What they
assert is not an expected output but the properties any correct answer has to have, which is the
only thing that can be checked without knowing the answer.

None of this needs a real sample. Real data is what would test whether the MODEL describes
biology; these test whether the CODE does what the model says.
"""

from __future__ import annotations

import math
import random
import shutil
import subprocess
import warnings

import numpy as np
import pytest

from conftest import load_script

gm = load_script("gconv_model")
gc = load_script("gene_conversion")
gcc = load_script("gconv_cohort")
pm = load_script("paralog_map")

CASES = 60          # per property; the whole file stays under a few seconds


# --------------------------------------------------------------------------- helpers


def _random_locus(rng, max_sites=18, max_reads=40):
    """A locus of random size, depth, signal strength and sparsity.

    Half of them carry a tract, because a property only checked on noise is only checked on the
    easy half of the input space, and the interesting failures are where there is something to
    get wrong.
    """
    n = int(rng.integers(1, max_sites))
    r = int(rng.integers(1, max_reads))
    scale = float(rng.choice([0.5, 8.0, 200.0]))

    # A fifth of the loci are degenerate on purpose. Random noise almost never produces a FLAT
    # or perfectly tied posterior, and that is exactly where the boundary a joint maximum picks
    # and the interval a marginal covers can disagree. Leaving these out let the tract be
    # reported outside its own credible interval and the suite still pass.
    shape = rng.random()
    degenerate = shape < 0.25
    if shape < 0.09:
        delta = np.zeros((r, n))                       # no information anywhere
        n = max(n, 8)                                  # a flat posterior needs room to be flat
        delta = np.zeros((r, n))
    elif shape < 0.17:
        delta = np.full((r, n), scale)                 # every site emphatically the donor's
    elif shape < 0.25:
        delta = np.full((r, n), -scale)                # every site emphatically the acceptor's
    else:
        delta = rng.normal(0, scale, (r, n)) * (rng.random((r, n)) > rng.uniform(0.2, 0.9))
    if n >= 4 and r >= 4 and rng.random() < 0.5:
        i = int(rng.integers(0, n - 2))
        j = int(rng.integers(i + 1, n))
        covered = rng.random(r) < 0.8
        delta[np.ix_(covered, np.arange(i, j + 1))] = scale
    # Evenly spaced positions half the time. Random spacing gives every interval a different
    # span, and the length prior then breaks every tie on its own: 0 of 200 degenerate loci with
    # random spacing reach the tied posterior where the boundary a joint maximum picks and the
    # interval a marginal covers disagree, while evenly spaced ones do it at once. The values
    # were being randomised and the STRUCTURE was not.
    # A degenerate delta ALWAYS gets even spacing, because the two only produce a tie together:
    # equal likelihoods need equal spans before the length prior stops deciding for them.
    if degenerate or rng.random() < 0.5:
        first = int(rng.integers(1, 100000))
        step = int(rng.choice([1, 40, 500]))
        positions = [first + step * i for i in range(n)]
    else:
        positions = sorted(rng.choice(200000, n, replace=False).tolist())
    return delta, positions


def _assert_sane(res, positions, label):
    """What any answer must satisfy, whatever the input was."""
    if res is None:
        return
    for k, v in res.items():
        if isinstance(v, float):
            assert math.isfinite(v), f"{label}: {k} is {v}"
    assert 0.0 <= res["post_conv"] <= 1.0, f"{label}: post_conv {res['post_conv']}"
    assert 0.0 <= res["post_tract"] <= 1.0, f"{label}: post_tract {res['post_tract']}"
    assert 0 <= res["map_i"] <= res["map_j"] < len(positions), f"{label}: MAP out of range"
    assert res["mut_rate"] > 0, f"{label}: mut_rate {res['mut_rate']}"
    assert min(gm.TRACT_AF_GRID) - 1e-9 <= res["tract_af"] <= 1 + 1e-9, f"{label}: tract_af"
    lo, hi = res["start_ci"]
    assert lo <= res["start"] <= hi, f"{label}: start {res['start']} outside {res['start_ci']}"
    lo, hi = res["end_ci"]
    assert lo <= res["end"] <= hi, f"{label}: end {res['end']} outside {res['end_ci']}"
    assert np.all(np.isfinite(res["site_post"])), f"{label}: site_post has non-finite entries"


# ------------------------------------------------------------------- the model


@pytest.mark.parametrize("seed", range(CASES))
def test_the_model_answers_sanely_whatever_it_is_given(seed):
    """The check that found the crash on a prior of 0, the silent wrong answer on unsorted
    positions, and the tract reported outside its own credible interval."""
    rng = np.random.default_rng(seed)
    delta, positions = _random_locus(rng)
    kw = {}
    if rng.random() < 0.3:
        kw["prior"] = float(rng.choice([0.0, 1e-6, 0.01, 0.5, 1.0]))
    if rng.random() < 0.3:
        kw["mut_rate"] = float(rng.choice([1e-9, 3e-4, 0.05, 0.24]))
    if rng.random() < 0.3:
        kw["mean_span_bp"] = float(rng.choice([1.0, 100.0, 1e6]))
    if rng.random() < 0.3:
        kw["max_grid"] = int(rng.choice([1, 5, 150]))
    if rng.random() < 0.3:
        kw["is_del"] = (rng.random(len(positions)) < 0.3).tolist()

    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)   # an overflow here is a defect, not noise
        res = gm.fit_locus(delta, positions, **kw)
    _assert_sane(res, positions, f"seed {seed} kw={kw}")


@pytest.mark.parametrize("seed", range(CASES))
def test_a_degenerate_setting_is_refused_and_never_answered_wrongly(seed):
    """The other half: inputs that have no right answer must raise where they were passed,
    not produce a plausible number somewhere downstream."""
    rng = np.random.default_rng(1000 + seed)
    delta, positions = _random_locus(rng)
    bad = rng.choice(["mut_rate_0", "mut_rate_1", "prior_neg", "span_0", "unsorted", "duplicate"])

    kw, pos = {}, list(positions)
    if bad == "mut_rate_0":
        kw["mut_rate"] = 0.0
    elif bad == "mut_rate_1":
        kw["mut_rate"] = 1.0
    elif bad == "prior_neg":
        kw["prior"] = -0.5
    elif bad == "span_0":
        kw["mean_span_bp"] = 0.0
    elif bad == "unsorted" and len(pos) > 1:
        pos[0], pos[-1] = pos[-1], pos[0]
    elif bad == "duplicate" and len(pos) > 1:
        pos[1] = pos[0]
    else:
        pytest.skip("locus too small for this mutation")

    with pytest.raises(ValueError):
        gm.fit_locus(delta, pos, **kw)


@pytest.mark.parametrize("seed", range(CASES))
def test_doubling_the_molecules_argues_harder_for_whatever_they_already_said(seed):
    """More of the same data makes the same case more strongly, in whichever direction.

    Not "more evidence for a tract": on random input the evidence points AGAINST one, and
    doubling it correctly makes the case against stronger. Writing the property the naive way
    failed on 11 of 60 random loci and the model was right every time, which is the sort of
    wrong assumption a generator catches and a hand-written example does not.

    Only loci that already say something clearly are checked. Near zero the prior terms and a
    switch of the winning tract dominate, and there is nothing to be monotone about.
    """
    rng = np.random.default_rng(2000 + seed)
    delta, positions = _random_locus(rng, max_sites=10, max_reads=20)
    if delta.shape[0] == 0:
        pytest.skip("no reads")

    one = gm.fit_locus(delta, positions)
    two = gm.fit_locus(np.vstack([delta, delta]), positions)
    if one is None or two is None or abs(one["log10_bf_null"]) < 1.0:
        pytest.skip("the locus says nothing clearly enough to be monotone about")

    assert math.copysign(1, two["log10_bf_null"]) == math.copysign(1, one["log10_bf_null"]), \
        f"seed {seed}: doubling the data changed which way the evidence points"
    assert abs(two["log10_bf_null"]) >= abs(one["log10_bf_null"]) - 1e-6, \
        f"seed {seed}: doubling the data weakened the case it already made"


@pytest.mark.parametrize("seed", range(CASES))
def test_the_prior_never_moves_the_bayes_factor(seed):
    """It is a statement about the data. Only the posterior may follow the prior."""
    rng = np.random.default_rng(3000 + seed)
    delta, positions = _random_locus(rng)

    a = gm.fit_locus(delta, positions, prior=0.001)
    b = gm.fit_locus(delta, positions, prior=0.9)
    if a is None or b is None:
        pytest.skip("nothing to compare")

    assert a["log10_bf"] == pytest.approx(b["log10_bf"], abs=1e-9), f"seed {seed}"
    assert b["post_conv"] >= a["post_conv"] - 1e-12


@pytest.mark.parametrize("seed", range(CASES))
def test_the_order_of_the_molecules_is_irrelevant(seed):
    """The likelihood is a product over reads, so it is symmetric in them."""
    rng = np.random.default_rng(4000 + seed)
    delta, positions = _random_locus(rng)
    if delta.shape[0] < 2:
        pytest.skip("needs two reads")

    order = rng.permutation(delta.shape[0])
    a = gm.fit_locus(delta, positions)
    b = gm.fit_locus(delta[order], positions)
    if a is None or b is None:
        pytest.skip("nothing to compare")

    assert a["log10_bf"] == pytest.approx(b["log10_bf"], abs=1e-9), f"seed {seed}"
    assert (a["map_i"], a["map_j"]) == (b["map_i"], b["map_j"]), f"seed {seed}"


# --------------------------------------------------------- the CIGAR arithmetic


def _oracle(pos, cigar, seq, qual, wanted):
    """An independent, deliberately slow reimplementation: build the whole map, then look up.

    Written literally rather than as a walk, so it cannot share an arithmetic bug with the code
    under test. The op semantics are shared knowledge; the offset bookkeeping is not, and that is
    where a walker goes wrong.
    """
    ops, num = [], ""
    for ch in cigar:
        if ch.isdigit():
            num += ch
        elif num:
            ops.append((int(num), ch))
            num = ""
    ref_to_read, deleted = {}, set()
    ref, qry = pos, 0
    for length, op in ops:
        for _ in range(length):
            if op in "M=X":
                ref_to_read[ref] = qry
                ref += 1
                qry += 1
            elif op == "D":
                deleted.add(ref)
                ref += 1
            elif op == "N":
                ref += 1
            elif op in "IS":
                qry += 1
    out = {}
    for p in wanted:
        if p in deleted:
            out[p] = (gc.GAP, gc.INDEL_PHRED)
            continue
        q = ref_to_read.get(p)
        if q is None or q >= len(seq):
            continue
        if not qual or qual == "*":
            out[p] = (seq[q].upper(), gc.NO_QUAL_PHRED)
        elif q < len(qual):
            out[p] = (seq[q].upper(), ord(qual[q]) - 33)
    return out


def test_the_cigar_walker_agrees_with_an_independent_reimplementation():
    """Five thousand random alignments, including the operations that shift one coordinate and
    not the other. When this is wrong it is silently wrong: the right position off the wrong
    base, which no output ever looks strange enough to catch."""
    rng = random.Random(9)
    for trial in range(5000):
        cig, seqlen = "", 0
        for _ in range(rng.randint(1, 7)):
            op = rng.choice("MIDNSHP=X")
            ln = rng.randint(1, 12)
            cig += f"{ln}{op}"
            if op in "MIS=X":
                seqlen += ln
        seq = "".join(rng.choice("ACGTNacgt") for _ in range(max(0, seqlen + rng.randint(-3, 1))))
        mode = rng.random()
        qual = ("*" if mode < 0.15 else "" if mode < 0.25 else
                "".join(chr(33 + rng.randint(0, 41)) for _ in range(len(seq) + rng.randint(-4, 1))))
        pos = rng.randint(1, 50)
        wanted = set(rng.sample(range(1, 90), rng.randint(1, 15)))

        got = gc.read_calls_at(pos, cig, seq, qual, wanted)
        assert got == _oracle(pos, cig, seq, qual, wanted), \
            f"trial {trial}: pos={pos} cigar={cig} seq={seq!r} qual={qual!r}"


# ----------------------------------------------------------------- end to end


@pytest.fixture
def samtools():
    exe = shutil.which("samtools")
    if not exe:
        pytest.skip("samtools is not on PATH")
    return exe


def _readlen(cigar):
    import re
    return sum(int(n) for n, o in re.findall(r"(\d+)([A-Z])", cigar) if o in "MIS=X")


@pytest.mark.parametrize("seed", range(12))
def test_the_tool_produces_a_well_formed_file_from_any_alignment(seed, samtools, tmp_path):
    """Random sites and random reads through the real CLI. Nothing here knows the right answer;
    it checks that whatever comes out is a file someone can read: the declared header, verdicts
    from the declared set, coordinates the right way round, intervals containing their point."""
    rng = random.Random(seed)
    lines = ["\t".join(["pair_id", "acceptor", "acc_pos", "acc_base", "donor", "don_pos",
                        "don_base", "strand", "kind", "length"])]
    allpos = []
    for pid in range(rng.randint(1, 3)):
        base = 1000 + pid * 5000
        pos = sorted(rng.sample(range(base, base + 900), rng.randint(1, 12)))
        allpos += pos
        for p in pos:
            a, d = rng.sample("ACGT", 2)
            if rng.random() < 0.2:
                lines.append("\t".join([str(pid), "chr", str(p), a, "chr", str(p + 20000),
                                        "-", rng.choice("+-"), "del", "1"]))
            else:
                lines.append("\t".join([str(pid), "chr", str(p), a, "chr", str(p + 20000),
                                        d, rng.choice("+-"), "snp", "1"]))
    recs = []
    for i in range(rng.randint(0, 40)):
        start = max(1, rng.choice(allpos) - rng.randint(0, 120))
        ln = rng.randint(20, 150)
        cig = rng.choice([f"{ln}M", f"5S{ln}M3D7M", f"{ln}M4I9M", f"3H{ln}M", f"{ln}M2N5M"])
        L = _readlen(cig)
        seq = "".join(rng.choice("ACGTN") for _ in range(L))
        qual = "*" if rng.random() < 0.2 else "".join(chr(33 + rng.randint(0, 41)) for _ in range(L))
        recs.append("\t".join([f"r{i}", str(rng.choice([0, 16, 99, 147, 0x100, 0x400])), "chr",
                               str(start), "0", cig, "*", "0", "0", seq, qual]))

    sam = tmp_path / "x.sam"
    sam.write_text("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr\tLN:100000\n"
                   + "\n".join(sorted(recs, key=lambda x: int(x.split("\t")[3]))) + "\n")
    bam = tmp_path / "t.bam"
    subprocess.run([samtools, "view", "-b", "-o", str(bam), str(sam)], check=True,
                   capture_output=True)
    subprocess.run([samtools, "index", str(bam)], check=True, capture_output=True)
    sites = tmp_path / "sites.tsv"
    sites.write_text("\n".join(lines) + "\n")
    out = tmp_path / "out.tsv"

    assert gc.main(["--sites", str(sites), "--bam", str(bam), "--sample", "S", "-o", str(out),
                    "--samtools", samtools, "--min-sites", "1"]) == 0

    body = [ln for ln in out.read_text().splitlines() if not ln.startswith("#")]
    assert body[0].split("\t") == gc.COLUMNS
    header = body[0].split("\t")
    for line in body[1:]:
        r = dict(zip(header, line.split("\t")))
        assert r["verdict"] in {"gene_conversion", "ambiguous", "mismapping", "coverage_shift"}
        assert int(r["start"]) <= int(r["end"])
        for key, point in (("start_ci", "start"), ("end_ci", "end")):
            if r[key]:
                lo, hi = (int(x) for x in r[key].split("-"))
                assert lo <= int(r[point]) <= hi, f"seed {seed}: {key} excludes its estimate"
        if r["post_conv"]:
            assert 0.0 <= float(r["post_conv"]) <= 1.0
        assert "\t" not in r["reason"] and "\n" not in r["reason"]


# --------------------------------------------------------------------------------------------
# The map parsers and the cohort pass, which had no randomised coverage at all.

VERDICTS = ["gene_conversion", "ambiguous", "mismapping", "coverage_shift"]


def _snp_text(rng):
    """`show-snps` output, including deliberate multi-base deletion runs.

    Drawing every position independently, which is the obvious way to write this, means a run of
    consecutive acceptor positions sharing one donor position essentially never comes up, so the
    branch that collapses a run into one site is never reached however the rows are checked
    afterwards. The runs are therefore built on purpose.
    """
    def row(p1, b1, b2, p2, strand):
        return "\t".join([str(p1), b1, b2, str(p2), "10", "10", "1", "1", "3000", "3000", "1",
                          strand, "chr", "chr"])

    out = []
    for _ in range(rng.randint(0, 10)):
        out.append(row(rng.randint(1, 5000), rng.choice(["A", "C", ".", "N"]),
                       rng.choice(["A", "G", ".", "N"]), rng.randint(1, 5000),
                       rng.choice(["1", "-1"])))
    for _ in range(rng.randint(0, 3)):
        p1, p2 = rng.randint(1, 4900), rng.randint(1, 5000)
        strand = rng.choice(["1", "-1"])
        for k in range(rng.randint(1, 6)):
            out.append(row(p1 + k, rng.choice("ACGT"), ".", p2, strand))
        if rng.random() < 0.5:          # a run broken by a gap in the acceptor positions
            out.append(row(p1 + 20, rng.choice("ACGT"), ".", p2, strand))
    return out


def _expected_deletions(rows):
    """What the deletion runs in `rows` have to come back as, worked out independently.

    Deliberately not an invariant over the output. A mistake in collapsing a run corrupts the
    length and the bases it was built from together, so the two stay consistent with each other
    and every check that compares one against the other passes. Only an expectation built from
    the INPUT can see it.
    """
    gaps = sorted((r[12], r[13], "-" if r[11] == "-1" else "+", int(r[0]), int(r[3]), r[1].upper())
                  for r in rows if r[2].upper() == "." and r[1].upper() != ".")
    want, i = {}, 0
    while i < len(gaps):
        t1, t2, strand, acc, don, base = gaps[i]
        j, bases = i, [base]
        while (j + 1 < len(gaps) and gaps[j + 1][:3] == (t1, t2, strand)
               and gaps[j + 1][3] == gaps[j][3] + 1 and gaps[j + 1][4] == don):
            j += 1
            bases.append(gaps[j][5])
        want[(t1, t2, strand, acc)] = "".join(bases)
        i = j + 1
    return want


@pytest.mark.parametrize("seed", range(CASES))
def test_the_map_parsers_survive_whatever_mummer_prints(seed):
    rng = random.Random(9000 + seed)
    lines = []
    for _ in range(rng.randint(0, 8)):
        if rng.random() < 0.15:
            lines.append(rng.choice(["", "\t", "NUCMER", "junk", "1\t2\t3"]))
        else:
            lines.append("\t".join(
                [str(rng.randint(1, 5000)) for _ in range(4)]
                + ["600", "600", f"{rng.uniform(80, 100):.2f}", "3000", "3000", "20", "20",
                   rng.choice(["chr", "ctg2"]), rng.choice(["chr", "ctg2"])]))
    pairs = pm.parse_coords("\n".join(lines) + "\n", min_identity=rng.choice([0.0, 95.0]),
                            min_length=rng.choice([0, 300]))
    for pair in pairs:
        assert pair["acc_start"] <= pair["acc_end"] and pair["don_start"] <= pair["don_end"]
        assert pair["length"] == pair["acc_end"] - pair["acc_start"] + 1
        assert pair["strand"] in "+-"

    snp = _snp_text(rng)
    sites, _skipped = pm.parse_snps("\n".join(snp) + "\n")
    pm.sites_within_pairs(sites, pairs)

    found = {(s["acceptor"], s["donor"], s["strand"], s["acc_pos"]): s
             for s in sites if s["kind"] == "del"}
    want = _expected_deletions([ln.split("\t") for ln in snp])
    assert set(found) == set(want), f"seed {seed}: the deletion sites are not the ones in the text"
    for key, bases in want.items():
        site = found[key]
        assert (site["acc_base"], site["length"]) == (bases, len(bases)), f"seed {seed}: {key}"

    for site in sites:
        assert site["kind"] in ("snp", "del")
        assert site["acc_pos"] >= 1 and site["don_pos"] >= 1 and site["strand"] in "+-"
        if site["kind"] == "del":
            assert site["don_base"] == gm.GAP
        else:
            assert site["length"] == 1 and len(site["acc_base"]) == 1
            assert site["acc_base"] != site["don_base"]


@pytest.mark.parametrize("seed", range(CASES))
def test_the_cohort_pass_annotates_every_row_it_is_given(seed):
    rng = random.Random(4000 + seed)
    rows = [{"sample": rng.choice(["A", "B", "C", ""]), "pair_id": str(rng.randint(0, 3)),
             "contig": rng.choice(["chr", "ctg2", ""]), "donor": "chr",
             "verdict": rng.choice(VERDICTS), "reason": "r",
             "start": str(rng.randint(1, 500)), "end": str(rng.randint(1, 500)),
             "don_start": rng.choice([str(rng.randint(1, 500)), "", "NA"]),
             "don_end": rng.choice([str(rng.randint(1, 500)), "", "NA"]),
             "log10_bf": rng.choice([str(round(rng.uniform(-9, 40), 2)), "", "NA"]),
             "log10_bf_vs_null": rng.choice([str(round(rng.uniform(-9, 900), 2)), ""]),
             "n_sites": rng.choice([str(rng.randint(0, 30)), "", "0"]),
             "tract_af": "0.9", "mismap_frac": rng.choice(["0.01", "", "NA"])}
            for _ in range(rng.randint(0, 15))]
    loci = [{"sample": rng.choice(["A", "B", "C"]), "pair_id": str(rng.randint(0, 3)),
             "log10_bf": rng.choice(["1.0", ""]), "mismap_frac": rng.choice(["0.02", ""]),
             "mut_rate": "0.0003"} for _ in range(rng.randint(0, 10))]

    out, n_samples = gcc.annotate(rows, loci, ubiquitous=rng.choice([0.1, 0.9]),
                                  min_samples=rng.choice([0, 5]),
                                  donor_margin=rng.choice([0.0, 1.0]))

    assert len(out) == len(rows), f"seed {seed}: the cohort pass must not drop or invent rows"
    for r in out:
        assert r["cohort_verdict"] in VERDICTS + ["reference_artifact"]
        for column in gcc.COHORT_COLUMNS:
            assert column in r, f"seed {seed}: {column} missing"
        assert r["is_representative"] in ("", 0, 1)
        assert r["donor_call"] in ("", "resolved", "ambiguous", "only candidate")
        if r["donor_rank"] != "":
            assert 1 <= r["donor_rank"] <= r["n_donors"], f"seed {seed}: rank outside the field"
        # None is the cohort size being unknown, which happens when no per-locus rows were given:
        # the samples with nothing to report are then uncountable, so there is no denominator.
        if n_samples is None:
            assert r["event_frac"] == "", f"seed {seed}: a fraction of an unknown cohort"
            assert r["cohort_verdict"] != "reference_artifact", \
                f"seed {seed}: recurrence applied without a cohort to be recurrent in"
        elif r["event_samples"]:
            assert r["event_samples"] <= n_samples, f"seed {seed}: more samples than the cohort has"
