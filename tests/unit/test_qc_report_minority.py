"""How low an allele frequency can be trusted (bin/qcreport/minority.py).

Calls below fixation are counted per sample with the reads they rest on, and, where the
samplesheet names libraries of the same DNA, measured against each other: a real minority is in
the DNA and another library calls it too; an error is not reproduced.
"""

from __future__ import annotations

from qcreport import minority as mn


def call(af, dp=100):
    return {"af": af, "dp": dp}


def test_calls_below_fixation_are_counted_with_the_reads_they_rest_on():
    variants = {"S": {"c:1": call(0.05, 60), "c:2": call(0.15, 100), "c:3": call(0.5, 40),
                      "c:4": call(1.0, 50)}}

    sec = mn.build_minority(variants)

    s = sec["samples"]["S"]
    assert s["n"] == 3, "the fixed call is not a minority"
    assert s["hist"] == [1, 1, 0, 1, 0, 0]
    assert s["le3"] == 1, "0.05 x 60 is 3 alternate reads"
    assert s["median_alt"] == 15


def _sheet(tmp_path, rows):
    p = tmp_path / "samples.tsv"
    p.write_text("sampleId\trunId\tr1\tr2\tdna_id\n" + "".join("\t".join(r) + "\n" for r in rows))
    return str(p)


def test_libraries_of_one_dna_are_found_and_those_sharing_reads_left_out(tmp_path):
    """A merged sample lists the runs it was merged from: it shares their reads, so its
    agreement with them proves nothing."""
    sheet = _sheet(tmp_path, [
        ["A", "r1", "a_1.fq", "a_2.fq", "DE1"],
        ["B", "r1", "b_1.fq", "b_2.fq", "DE1"],
        ["M", "r1", "a_1.fq", "a_2.fq", "DE1"],
        ["M", "r2", "b_1.fq", "b_2.fq", "DE1"],
        ["C", "r1", "c_1.fq", "c_2.fq", "DE2"],
    ])

    col, value, reads = mn.replicate_sets(sheet)
    pairs, shared = mn.replicate_pairs(value, reads, {"A", "B", "M", "C"})

    assert col == "dna_id"
    assert reads["M"] == {"a_1.fq", "a_2.fq", "b_1.fq", "b_2.fq"}
    assert pairs == [("A", "B")] and shared == 2


def test_libraries_of_one_dna_on_different_references_are_not_paired():
    """Their coordinates do not correspond; one of them went to the wrong genome."""
    value = {"A": "DE1", "B": "DE1", "C": "DE1"}
    reads = {"A": {"a.fq"}, "B": {"b.fq"}, "C": {"c.fq"}}

    pairs, _ = mn.replicate_pairs(value, reads, {"A", "B", "C"}, ref={"A": "R1", "B": "R1", "C": "R2"})

    assert pairs == [("A", "B")]


def test_a_samplesheet_without_a_dna_column_has_no_replicates(tmp_path):
    p = tmp_path / "s.tsv"
    p.write_text("sampleId\tr1\nA\ta.fq\n")

    assert mn.replicate_sets(str(p))[0] is None


def test_reproducibility_counts_only_where_the_other_library_could_have_called_it():
    variants = {"A": {"c:1": call(0.05), "c:2": call(0.06), "c:3": call(0.4), "c:4": call(0.05),
                      "c:5": call(0.3), "c:9": call(1.0)},
                "B": {"c:3": call(0.35), "c:9": call(1.0)}}
    cells = {("c:1", "B"): (0.0, 200),     # read deeply, and not there: an error
             ("c:2", "B"): (None, 2),      # not read: says nothing
             ("c:4", "B"): (0.0, 60),      # read, but 5% of 60 is 3 reads: too thin to tell
             ("c:5", "B"): (0.0, 7)}       # at the consensus minimum: not read

    rep = mn.build_minority(variants, [("A", "B")], 0, cells, "dna_id", min_dp=7)["replicates"]

    assert rep["tested"][0] == 1 and rep["reproduced"][0] == 0
    assert rep["tested"][3] == 2 and rep["reproduced"][3] == 2, "c:3 in both directions"
    assert (rep["fixed_tested"], rep["fixed_reproduced"]) == (2, 2)


def test_an_na_cell_of_the_other_library_is_a_call_it_made():
    """NA in the matrix is a call whose fraction the files do not keep, not a site read without it."""
    from qcreport.series import CALLED
    variants = {"A": {"c:1": call(0.4)}, "B": {}}

    rep = mn.build_minority(variants, [("A", "B")], 0, {("c:1", "B"): (CALLED, 80)}, "dna_id",
                            min_dp=7)["replicates"]

    assert rep["tested"][3] == 1 and rep["reproduced"][3] == 1


def test_without_the_matrix_the_replicates_are_not_measured():
    variants = {"A": {"c:1": call(0.05)}, "B": {}}

    rep = mn.build_minority(variants, [("A", "B")], 0, {}, "dna_id", checked=False)["replicates"]

    assert rep["tested"] is None and rep["pairs"] == 1


def test_needed_cells_asks_each_library_about_the_other_ones_calls():
    variants = {"A": {"c:1": call(0.05), "c:3": call(0.4)}, "B": {"c:3": call(0.4), "c:5": call(0.1)}}

    assert mn.needed_cells([("A", "B")], variants) == {"c:1": {"B"}, "c:5": {"A"}}
