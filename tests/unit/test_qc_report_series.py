"""What each series gained since its first time point (bin/qcreport/series.py).

A fixed SNP in a later sample is new only if the first time point was read at that site and did
not carry it. These tests pin that reading, the losses, and what happens without the depth the
SNP matrix provides.
"""

from __future__ import annotations

from qcreport import series as ser


def meta(**samples):
    """{sample: (group, time)} -> the parse_metadata() shape."""
    return {s: {"group": g, "time": str(t), "tnum": float(t)} for s, (g, t) in samples.items()}


def var(af, gene="g", aa=""):
    return {"af": af, "dp": 50, "gene": gene, "aa": aa, "eff": "missense_variant", "alt": "T"}


M = meta(p0=("line1", 0), p3=("line1", 3), p6=("line1", 6))


def test_a_site_read_at_the_start_without_the_allele_is_new():
    variants = {"p0": {}, "p3": {"c:10": var(1.0, "atpE", "p.Ile66Met")}, "p6": {}}
    cells = {("c:10", "p0"): (0.0, 40)}

    sec = ser.build_series(M, variants, cells)

    row = sec["groups"][0]["rows"][0]
    assert (row["s"], row["new"], row["risen"], row["unknown"]) == ("p3", ["c:10"], [], 0)
    assert sec["sites"]["c:10"]["gene"] == "atpE" and sec["checked"] is True


def test_a_site_the_start_did_not_read_is_unknown_not_new():
    variants = {"p0": {}, "p3": {"c:10": var(1.0)}, "p6": {}}
    cells = {("c:10", "p0"): (None, 2)}

    row = ser.build_series(M, variants, cells, min_dp=7)["groups"][0]["rows"][0]

    assert (row["new"], row["unknown"]) == ([], 1)


def test_a_minority_at_the_start_that_became_fixed_has_risen():
    variants = {"p0": {"c:10": var(0.2)}, "p3": {"c:10": var(0.95)}, "p6": {}}

    row = ser.build_series(M, variants, {})["groups"][0]["rows"][0]

    assert (row["new"], row["risen"]) == ([], ["c:10"])


def test_a_site_fixed_from_the_start_is_inherited_not_gained():
    variants = {"p0": {"c:10": var(1.0)}, "p3": {"c:10": var(1.0)}, "p6": {"c:10": var(1.0)}}

    rows = ser.build_series(M, variants, {})["groups"][0]["rows"]

    assert all(r["new"] == [] and r["risen"] == [] for r in rows)


def test_a_site_fixed_at_the_start_and_read_without_it_later_is_lost():
    variants = {"p0": {"c:10": var(1.0)}, "p3": {"c:10": var(1.0)}, "p6": {}}
    cells = {("c:10", "p6"): (0.0, 55)}

    rows = {r["s"]: r for r in ser.build_series(M, variants, cells)["groups"][0]["rows"]}

    assert rows["p6"]["lost"] == ["c:10"] and rows["p3"]["lost"] == []


def test_without_the_matrix_a_site_not_called_at_the_start_counts_as_new_and_says_so():
    variants = {"p0": {"c:10": var(1.0)}, "p3": {"c:20": var(1.0)}, "p6": {}}

    sec = ser.build_series(M, variants, {}, checked=False)

    rows = {r["s"]: r for r in sec["groups"][0]["rows"]}
    assert rows["p3"]["new"] == ["c:20"] and sec["checked"] is False
    assert rows["p6"]["lost"] == [], "a loss needs the depth to be told from a site not read"


def test_every_sample_of_the_first_time_point_is_the_start():
    m = meta(a=("L", 0), b=("L", 0), c=("L", 5))
    variants = {"a": {}, "b": {"c:10": var(1.0)}, "c": {"c:10": var(1.0)}}

    g = ser.build_series(m, variants, {})["groups"][0]

    assert g["base"] == ["a", "b"] and g["rows"][0]["new"] == []


def test_a_group_with_a_single_time_point_is_not_a_series():
    m = meta(a=("L", 0), b=("L", 0))

    assert ser.build_series(m, {"a": {}, "b": {}}, {}) is None


def test_the_series_carries_its_treatment():
    m = meta(a=("L", 0), b=("L", 3), c=("L", 6))
    sm = {"tx_field": "treatment", "rows": {"a": {"treatment": "CTRL"}, "b": {"treatment": "BDQ"},
                                            "c": {"treatment": "BDQ"}}}

    g = ser.build_series(m, {"a": {}, "b": {}, "c": {}}, {}, sample_meta=sm)["groups"][0]

    assert g["tx"] == "BDQ"


def test_needed_cells_asks_for_the_start_at_later_fixed_sites_and_back():
    variants = {"p0": {"c:5": var(1.0)}, "p3": {"c:10": var(1.0), "c:11": var(0.3)}, "p6": {}}

    need = ser.needed_cells(M, variants)

    assert need == {"c:10": {"p0"}, "c:5": {"p3", "p6"}}


def test_read_matrix_cells_reads_only_what_is_needed(tmp_path):
    p = tmp_path / "m.tsv"
    p.write_text("reference\tcontig\tpos\tref_allele\talt_allele\tgene\teffect\taa_change\t"
                 "p0|AF\tp0|DP\tp3|AF\tp3|DP\n"
                 "R\tc\t10\tA\tT\t\t\t\t0\t40\t1.0000\t38\n"
                 "R\tc\t20\tA\tT\t\t\t\t\t0\t1.0000\t30\n")

    cells = ser.read_matrix_cells(str(p), {"c:10": {"p0"}, "c:20": {"p0"}})

    assert cells == {("c:10", "p0"): (0.0, 40), ("c:20", "p0"): (None, 0)}
