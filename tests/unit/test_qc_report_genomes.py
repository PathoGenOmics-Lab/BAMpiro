"""Each reference's own genome for the report (bin/qcreport/genomes.py).

A sample's landscape profile is binned over its own reference, contigs end to end in the order of the
FASTA index, so the genes and masked regions drawn under it must be its reference's and placed the same way.
"""

from __future__ import annotations

from qcreport import genomes as gm


def _fai(path, contigs):
    path.write_text("".join(f"{c}\t{n}\t0\t60\t61\n" for c, n in contigs))
    return str(path)


def _gff(path, rows):
    path.write_text("##gff-version 3\n" + "".join(
        f"{c}\t.\t{t}\t{s}\t{e}\t.\t+\t.\t{a}\n" for c, t, s, e, a in rows))
    return str(path)


def test_a_gene_of_the_second_contig_sits_after_the_whole_first_one(tmp_path):
    fai = _fai(tmp_path / "r.fa.fai", [("chr", 1000), ("plasmid", 200)])
    gff = _gff(tmp_path / "r.gff", [("chr", "gene", 10, 90, "Name=dnaA"), ("plasmid", "gene", 5, 50, "Name=repA")])

    g = gm.build_genome(gff, None, fai)

    assert g["len"] == 1200
    assert g["genes"] == [{"name": "dnaA", "start": 10, "end": 90}, {"name": "repA", "start": 1005, "end": 1050}]


def test_a_single_contig_reference_takes_its_genes_whatever_the_gff_calls_the_chromosome(tmp_path):
    fai = _fai(tmp_path / "r.fa.fai", [("NC_000962.3", 4411532)])
    gff = _gff(tmp_path / "r.gff", [("Chromosome", "gene", 759807, 763325, "Name=rpoB")])

    assert gm.build_genome(gff, None, fai)["genes"] == [{"name": "rpoB", "start": 759807, "end": 763325}]


def test_a_gene_on_a_contig_the_reference_lacks_is_left_out(tmp_path):
    fai = _fai(tmp_path / "r.fa.fai", [("chr", 1000), ("plasmid", 200)])
    gff = _gff(tmp_path / "r.gff", [("other", "gene", 10, 90, "Name=x")])

    assert gm.build_genome(gff, None, fai)["genes"] == []


def test_a_gene_is_preferred_to_its_cds(tmp_path):
    gff = _gff(tmp_path / "r.gff", [("chr", "CDS", 12, 88, "Name=dnaA"), ("chr", "gene", 10, 90, "Name=dnaA")])

    assert gm.build_genome(gff, None, None, length=1000)["genes"] == [{"name": "dnaA", "start": 10, "end": 90}]


def test_masked_regions_are_placed_on_the_same_axis_and_binned_over_it(tmp_path):
    fai = _fai(tmp_path / "r.fa.fai", [("chr", 1000), ("plasmid", 1000)])
    mask = tmp_path / "Locus_to_exclude_R.txt"
    mask.write_text("plasmid\t1\t500\trepeat\t\n")

    g = gm.build_genome(None, str(mask), fai, nbins=4)

    assert g["mask_iv"] == [(1001, 1500)]
    assert g["mask_bins"][:2] == [0.0, 0.0] and g["mask_bins"][2] > 0.99
    assert round(g["mask_pct"]) == 25


def test_each_reference_gets_its_own_genes_and_length(tmp_path):
    fa = _fai(tmp_path / "a.fa.fai", [("A", 1000)])
    ga = _gff(tmp_path / "a.gff", [("A", "gene", 10, 90, "Name=atpE")])
    gb = _gff(tmp_path / "b.gff", [("B", "gene", 20, 80, "Name=E1_01312")])

    out = gm.build_genomes({"RA": ga, "RB": gb}, {}, {"RA": fa}, lengths={"RB": 900})

    assert set(out) == {"RA", "RB"}
    assert out["RA"]["len"] == 1000 and out["RA"]["genes"][0]["name"] == "atpE"
    assert out["RB"]["len"] == 900, "without an index, the length the samples report"
    assert out["RB"]["genes"][0]["name"] == "E1_01312"


def test_placeholders_and_missing_files_give_an_empty_genome(tmp_path):
    (tmp_path / "NO_FILE_GFF").write_text("")

    g = gm.build_genome(str(tmp_path / "NO_FILE_GFF"), str(tmp_path / "absent.bed"), None, length=500)

    assert (g["len"], g["genes"], g["mask_bins"], g["mask_iv"]) == (500, [], None, None)


def test_pairs_reads_ref_equals_path_and_skips_the_rest():
    assert gm.pairs(["E1=a.gff", "E2=dir/b=c.gff", "bad", "=x", "E3="]) == {"E1": "a.gff", "E2": "dir/b=c.gff"}
