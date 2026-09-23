"""Relatedness and GROUP_MISMATCH (bin/qcreport/relatedness.py).

A sample the samplesheet puts in a group but that sits far from all of it, in a group whose
members are normally close, is swapped, mislabelled, contaminated or reinfected. These tests pin
when that is said, and when it is not: groups that are not genetic units, groups of two, and
pairs that rest on too few positions.
"""

from __future__ import annotations

from qcreport import relatedness as rel


def dist(pairs, variable=100, ref="R"):
    """A parse_pairs() result from {(a, b): snps}, every pair comparing all variable positions."""
    out = {"ref": {}, "variable": {ref: variable}, "pairs": {}}
    for (a, b), v in pairs.items():
        snps, comp = v if isinstance(v, tuple) else (v, variable)
        out["pairs"][(a, b)] = (snps, comp)
        out["ref"][a] = out["ref"][b] = ref
    return out


def full(groups_of, within=2, between=2000, override=None):
    """Every pair of the samples in `groups_of`: `within` SNPs inside a group, `between` across."""
    names = sorted(groups_of)
    pairs = {}
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            pairs[(a, b)] = within if groups_of[a] == groups_of[b] else between
    pairs.update(override or {})
    return dist(pairs)


LINES = {"a1": "A", "a2": "A", "a3": "A", "a4": "A", "b1": "B", "b2": "B", "b3": "B"}


def test_a_sample_far_from_its_group_and_close_to_another_is_flagged():
    """a4 is labelled line A but carries line B's genome: 2,000 SNPs from its line, 3 from b1."""
    over = {("a1", "a4"): 2000, ("a2", "a4"): 2000, ("a3", "a4"): 2000, ("a4", "b1"): 3,
            ("a4", "b2"): 5, ("a4", "b3"): 4}

    out = rel.group_mismatches(full(LINES, override=over), LINES, threshold=12)

    assert list(out) == ["a4"]
    assert out["a4"] == {"g": "A", "din": 2000, "nin": "a1", "dout": 3, "nout": "b1", "gout": "B"}


def test_a_sample_far_from_its_group_and_from_everyone_else_is_flagged_too():
    """A contaminant, or a culture of another strain: close to no sample at all. On a real
    cohort the samples typed as another lineage sat 47 to 72 SNPs from each other and 2,000 from
    their line mates, so requiring a close sample of another group missed every one of them."""
    over = {("a1", "a4"): 2000, ("a2", "a4"): 2000, ("a3", "a4"): 2000,
            ("a4", "b1"): 60, ("a4", "b2"): 64, ("a4", "b3"): 70}

    out = rel.group_mismatches(full(LINES, override=over), LINES, threshold=12)

    assert out["a4"]["din"] == 2000 and (out["a4"]["dout"], out["a4"]["nout"]) == (60, "b1")


def test_a_consistent_cohort_flags_nobody():
    assert rel.group_mismatches(full(LINES), LINES) == {}


def test_a_group_that_is_not_a_genetic_unit_is_not_judged():
    """A column like 'hospital' also reads as a group; its members need not be related."""
    hospitals = {"s1": "H1", "s2": "H1", "s3": "H1", "s4": "H2", "s5": "H2", "s6": "H2"}
    pairs = {(a, b): 500 for a in hospitals for b in hospitals if a < b}
    pairs[("s1", "s4")] = 2                     # a transmission across hospitals

    assert rel.group_mismatches(dist(pairs), hospitals) == {}


def test_a_group_of_two_cannot_say_which_member_is_out_of_place():
    groups = {"a1": "A", "a2": "A", "b1": "B", "b2": "B", "b3": "B"}
    over = {("a1", "a2"): 2000, ("a2", "b1"): 3}

    assert rel.group_mismatches(full(groups, override=over), groups) == {}


def test_a_pair_resting_on_few_positions_does_not_count():
    """A distance over a handful of positions says nothing, near or far."""
    over = {("a1", "a4"): (2000, 10), ("a2", "a4"): (2000, 10), ("a3", "a4"): (2000, 10)}

    out = rel.group_mismatches(full(LINES, override=over), LINES, min_compared=0.5)

    assert out == {}, "a4's distances to its group compared 10 of 100 variable positions"


def test_samples_of_another_reference_are_never_compared():
    d = full(LINES, override={("a1", "a4"): 2000, ("a2", "a4"): 2000, ("a3", "a4"): 2000,
                              ("a4", "b1"): 900, ("a4", "b2"): 900, ("a4", "b3"): 900})
    d["ref"]["x"] = "OTHER"
    d["pairs"][("a4", "x")] = (1, 100)
    groups = dict(LINES, x="B")

    assert rel.group_mismatches(d, groups)["a4"]["nout"] == "b1", "x is on another reference"


def test_parse_pairs_reads_the_distances_tsv(tmp_path):
    p = tmp_path / "pairs.tsv"
    p.write_text("sample_a\tsample_b\treference\tsnps\tcompared\tvariable_sites\n"
                 "S1\tS2\tR\t4\t90\t100\nS1\tS3\tR\t7\t88\t100\n")

    d = rel.parse_pairs(str(p))

    assert d["pairs"][("S1", "S2")] == (4, 90)
    assert d["ref"] == {"S1": "R", "S2": "R", "S3": "R"} and d["variable"] == {"R": 100}


def test_parse_pairs_of_a_placeholder_is_empty(tmp_path):
    p = tmp_path / "NO_FILE_DISTANCES"
    p.write_text("")

    assert rel.parse_pairs(str(p)) == {"ref": {}, "variable": {}, "pairs": {}}


def test_build_relatedness_lays_out_each_references_upper_triangle():
    d = dist({("S1", "S2"): (4, 90), ("S1", "S3"): (7, 50), ("S2", "S3"): (9, 100)})

    sec = rel.build_relatedness(d, {"S1": "P1", "S2": "P1"}, threshold=5)

    r = sec["refs"]["R"]
    assert r["samples"] == ["S1", "S2", "S3"]
    assert r["snps"] == [4, 7, 9] and r["cmp"] == [90, 50, 100]
    assert sec["threshold"] == 5 and sec["group"] == {"S1": "P1", "S2": "P1"}


def test_build_relatedness_without_distances_is_none():
    assert rel.build_relatedness({"ref": {}, "variable": {}, "pairs": {}}) is None
