"""Deletions and callable depth along the genome (bin/qcreport/coverage.py).

A stretch no read covers means different things depending on the rest of the cohort: nobody reads
it (a part of the reference these genomes lack, or a repeat), several samples lack it (a shared
deletion), or one sample does (a private one). These tests pin that reading.
"""

from __future__ import annotations

import pytest

from conftest import load_script
from qcreport import coverage as cov

dp = load_script("depth_profile")


def profile(zero=(), median=40.0, ref="REF", length=10_000, genes=(), windows=None):
    return {"reference": ref, "median_dp": median, "contigs": [("chr", length)],
            "zero": [("chr", s, e) for s, e in zero],
            "genes": [{"contig": "chr", "start": s, "end": e, "gene": g, "locus": "", "breadth": b, "rel": None}
                      for g, s, e, b in genes],
            "windows": windows or []}


def cohort(n, **kw):
    return {f"S{i}": profile(**kw) for i in range(n)}


# ------------------------------------------------------------------------------------ classes


def test_a_stretch_one_sample_lacks_is_its_private_deletion():
    profiles = cohort(5)
    profiles["S0"] = profile(zero=[(1001, 3000)])

    section, _, _ = cov.build_coverage(profiles, min_len=200)

    [r] = section["regions"]
    assert (r["cls"], r["n"], r["of"], r["start"], r["end"]) == ("private", 1, 5, 1001, 3000)
    assert r["samples"] == [["S0", 1001, 3000]]


def test_a_stretch_several_samples_lack_is_a_shared_deletion():
    profiles = cohort(6)
    for s in ("S0", "S1"):
        profiles[s] = profile(zero=[(1001, 3000)])

    [r] = cov.build_coverage(profiles)[0]["regions"]

    assert (r["cls"], r["n"]) == ("shared", 2)


def test_a_stretch_nobody_reads_is_not_anybodys_deletion():
    """A repeat, or a part of the reference none of these genomes has."""
    profiles = cohort(10, zero=[(5001, 5800)])

    section, tracks, _ = cov.build_coverage(profiles)

    [r] = section["regions"]
    assert (r["cls"], r["n"], r["of"]) == ("cohort", 10, 10)
    assert all(v == 0 for v in tracks["S0"]), "nobody's deletion lights up nobody's track"


def test_the_only_sample_on_a_reference_has_nothing_to_compare_with():
    profiles = cohort(4)
    profiles["X"] = profile(zero=[(1001, 3000)], ref="OTHER")

    regions = cov.build_coverage(profiles)[0]["regions"]

    assert [(r["ref"], r["cls"]) for r in regions] == [("OTHER", "alone")]


def test_a_thinly_read_sample_is_not_assessed():
    """At 3x, stretches without reads turn up by chance: they would all read as deletions."""
    profiles = cohort(4)
    profiles["thin"] = profile(zero=[(1001, 3000), (6001, 7000)], median=3.0)

    section, tracks, _ = cov.build_coverage(profiles, min_depth=10)

    assert section["regions"] == []
    assert section["not_assessed"] == ["thin"]
    assert "thin" not in tracks


def test_a_stretch_shorter_than_the_minimum_is_not_a_deletion():
    profiles = cohort(4)
    profiles["S0"] = profile(zero=[(1001, 1150)])

    assert cov.build_coverage(profiles, min_len=200)[0]["regions"] == []


def test_overlapping_stretches_of_different_samples_are_one_region():
    profiles = cohort(6)
    profiles["S0"] = profile(zero=[(1001, 3000)])
    profiles["S1"] = profile(zero=[(1200, 3400)])

    [r] = cov.build_coverage(profiles)[0]["regions"]

    assert (r["start"], r["end"], r["n"]) == (1001, 3400, 2)
    assert r["samples"] == [["S0", 1001, 3000], ["S1", 1200, 3400]]


def test_a_region_names_the_genes_its_carriers_do_not_read():
    profiles = cohort(4, genes=[("katG", 1101, 2500, 1.0), ("rpoB", 6001, 7000, 1.0)])
    profiles["S0"] = profile(zero=[(1001, 3000)], genes=[("katG", 1101, 2500, 0.0), ("rpoB", 6001, 7000, 1.0)])

    [r] = cov.build_coverage(profiles)[0]["regions"]

    assert r["genes"] == ["katG"]


# ------------------------------------------------------------------------------------ genes


def test_a_gene_most_samples_read_and_one_does_not_is_lost_in_that_one():
    profiles = cohort(4, genes=[("katG", 1101, 2500, 0.99)])
    profiles["S0"] = profile(genes=[("katG", 1101, 2500, 0.2)])

    [g] = cov.build_coverage(profiles)[0]["genes_lost"]

    assert (g["gene"], g["samples"], g["of"]) == ("katG", [["S0", 0.2]], 4)


def test_a_gene_nobody_reads_is_not_lost_by_anyone():
    profiles = cohort(4, genes=[("absent", 1101, 2500, 0.0)])

    assert cov.build_coverage(profiles)[0]["genes_lost"] == []


# ------------------------------------------------------------------------------------ binning


def test_spread_uses_the_landscapes_own_binning():
    """stats_to_legacy puts a variant at global position g in bin g * nbins // genome length."""
    contigs = [("a", 60), ("b", 40)]
    for contig, pos in (("a", 1), ("a", 9), ("a", 10), ("a", 60), ("b", 1), ("b", 40)):
        g = dict(a=0, b=60)[contig] + pos
        bins = cov._spread([(contig, pos, pos, 1.0)], contigs, nbins=10)
        assert bins.index(1.0) == min(9, g * 10 // 100), (contig, pos)


def test_spread_splits_an_interval_across_bins():
    bins = cov._spread([("a", 5, 25, 1.0)], [("a", 100)], nbins=10)

    assert sum(bins) == 21
    assert bins[:3] == [5.0, 10.0, 6.0]


def test_snp_per_kb_divides_by_the_callable_positions_of_each_bin():
    windows = [("chr", 1, 500, 40.0, 0.0, 1.0), ("chr", 501, 1000, 40.0, 0.5, 0.5)]
    p = profile(length=1000, windows=windows)

    dens = cov.snp_per_kb([2, 1], p, nbins=2)

    # The landscape's bins are 1-based (499 and 501 positions here), hence the approximation.
    assert dens == [pytest.approx(4.0, rel=0.01), pytest.approx(4.0, rel=0.01)], \
        "1 SNP in 250 callable bp is 4 per kb"


def test_snp_per_kb_is_empty_where_almost_nothing_could_be_called():
    windows = [("chr", 1, 1000, 1.0, 0.95, 0.05)]

    assert cov.snp_per_kb([3], profile(length=1000, windows=windows), nbins=1) == [None]


def test_the_deletion_track_is_the_share_of_each_bin_a_sample_lacks():
    profiles = cohort(4, length=1000)
    profiles["S0"] = profile(zero=[(1, 250)], length=1000)

    _, tracks, _ = cov.build_coverage(profiles, nbins=2)

    assert tracks["S0"] == [pytest.approx(0.5, abs=0.01), 0.0]


# ------------------------------------------------------------------------------------ files


def _write_allpos(path, sample, depths):
    head = ("##fileformat=VCFv4.2\n##contig=<ID=chr,length=%d>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\t"
            "INFO\tFORMAT\t%s\n" % (len(depths), sample))
    body = "".join(f"chr\t{i + 1}\t.\tA\t.\t.\t.\t.\tGT:DP:AD\t0:{d}:{d},0\n" for i, d in enumerate(depths))
    path.write_text(head + body)


def test_the_tables_of_depth_profile_are_read_back_by_the_report(tmp_path):
    gff = tmp_path / "ref.gff"
    gff.write_text("chr\t.\tgene\t400\t600\t.\t+\t.\tName=lost\n")
    paths = []
    for i in range(4):
        depths = [30] * 1000 if i else [30] * 300 + [0] * 400 + [30] * 300
        _write_allpos(tmp_path / f"S{i}.vcf", f"S{i}", depths)
        assert dp.main(["--vcf", str(tmp_path / f"S{i}.vcf"), "--gff", str(gff), "--reference", "R",
                        "--window", "100", "--out-prefix", str(tmp_path / f"S{i}")]) == 0
        paths += [str(tmp_path / f"S{i}.{k}.tsv") for k in ("depth_windows", "zero_depth", "gene_depth")]

    profiles = cov.parse_depth_profiles(paths)

    assert sorted(profiles) == ["S0", "S1", "S2", "S3"]
    assert profiles["S0"]["median_dp"] == 30 and profiles["S0"]["contigs"] == [("chr", 1000)]
    assert len(profiles["S0"]["windows"]) == 10 and profiles["S0"]["zero"] == [("chr", 301, 700)]
    section, _, regions = cov.build_coverage(profiles)
    assert [(r["cls"], r["genes"], r["samples"]) for r in regions] == [("private", ["lost"], [["S0", 301, 700]])]
    assert section["genes_lost"][0]["gene"] == "lost"

    out = tmp_path / "deletions.tsv"
    cov.write_deletions(out, regions)
    assert out.read_text().splitlines()[1].split("\t") == [
        "D1", "R", "chr", "301", "700", "400", "private", "1", "4", "S0:301-700", "lost"]
