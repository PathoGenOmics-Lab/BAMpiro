"""Deletions and callable depth along the genome (bin/qcreport/coverage.py).

A stretch no read covers means different things depending on the rest of the cohort: nobody reads
it (a part of the reference these genomes lack, or a repeat), several samples lack it (a shared
deletion), or one sample does (a private one). These tests pin that reading.
"""

from __future__ import annotations

from array import array

import pytest

from conftest import load_script
from qcreport import coverage as cov

dp = load_script("depth_profile")

GENES = [("chr", 1101, 2500, "katG", "T1"), ("chr", 6001, 7000, "rpoB", "T2")]


def profile(zero=(), median=40.0, ref="REF", length=10_000, breadth=None, callable_bins=None, sample=None):
    return {"sample": sample, "reference": ref, "median_dp": median, "contigs": [("chr", length)],
            "zero": [("chr", s, e) for s, e in zero], "callable": callable_bins,
            "genes": array("f", breadth) if breadth is not None else None}


def cohort(n, ref="REF", **kw):
    return {(f"S{i}", ref): profile(sample=f"S{i}", ref=ref, **kw) for i in range(n)}


def run(profiles, **kw):
    return cov.build_coverage(profiles, {"REF": GENES}, **kw)


# ------------------------------------------------------------------------------------ classes


def test_a_stretch_one_sample_lacks_is_its_private_deletion():
    profiles = cohort(5)
    profiles[("S0", "REF")] = profile(zero=[(1001, 3000)])

    section, _, _ = run(profiles, min_len=200)

    [r] = section["regions"]
    assert (r["cls"], r["n"], r["of"], r["start"], r["end"]) == ("private", 1, 5, 1001, 3000)
    assert r["samples"] == [["S0", 1001, 3000]]


def test_a_stretch_several_samples_lack_is_a_shared_deletion():
    profiles = cohort(6)
    for s in ("S0", "S1"):
        profiles[(s, "REF")] = profile(zero=[(1001, 3000)])

    [r] = run(profiles)[0]["regions"]

    assert (r["cls"], r["n"]) == ("shared", 2)


def test_a_stretch_nobody_reads_is_not_anybodys_deletion():
    """A repeat, or a part of the reference none of these genomes has."""
    profiles = cohort(10, zero=[(5001, 5800)])

    section, tracks, _ = run(profiles)

    [r] = section["regions"]
    assert (r["cls"], r["n"], r["of"]) == ("cohort", 10, 10)
    assert all(v == 0 for v in tracks[("S0", "REF")]), "nobody's deletion lights up nobody's track"


def test_the_only_sample_on_a_reference_has_nothing_to_compare_with():
    profiles = cohort(4)
    profiles[("X", "OTHER")] = profile(zero=[(1001, 3000)], ref="OTHER")

    section, tracks, _ = run(profiles)

    assert [(r["ref"], r["cls"]) for r in section["regions"]] == [("OTHER", "alone")]
    assert not any(tracks[("X", "OTHER")]), "a stretch nothing compares with is nobody's deletion either"


def test_a_thinly_read_sample_is_not_assessed():
    """At 3x, stretches without reads turn up by chance: they would all read as deletions."""
    profiles = cohort(4)
    profiles[("thin", "REF")] = profile(zero=[(1001, 3000), (6001, 7000)], median=3.0)

    section, tracks, _ = run(profiles, min_depth=10)

    assert section["regions"] == []
    assert section["not_assessed"] == ["thin"]
    assert ("thin", "REF") not in tracks


def test_a_stretch_shorter_than_the_minimum_is_not_a_deletion():
    profiles = cohort(4)
    profiles[("S0", "REF")] = profile(zero=[(1001, 1150)])

    assert run(profiles, min_len=200)[0]["regions"] == []


def test_stretches_that_overlap_each_other_by_half_are_one_deletion():
    profiles = cohort(6)
    profiles[("S0", "REF")] = profile(zero=[(1001, 3000)])
    profiles[("S1", "REF")] = profile(zero=[(1200, 3400)])

    [r] = run(profiles)[0]["regions"]

    assert r["n"] == 2 and r["cls"] == "shared"
    assert r["samples"] == [["S0", 1001, 3000], ["S1", 1200, 3400]], "each carrier keeps its own stretch"


def test_a_private_deletion_is_not_swallowed_by_a_gap_every_sample_shares():
    """Nine of ten samples lack 4000-4500 (a repeat); one sample alone lacks 2000-8000."""
    profiles = cohort(10, zero=[(4000, 4500)])
    profiles[("S9", "REF")] = profile(zero=[(2000, 8000)])

    regions = run(profiles)[0]["regions"]

    by = {(r["start"], r["end"]): r for r in regions}
    assert by[(2000, 8000)]["cls"] == "private" and by[(2000, 8000)]["samples"] == [["S9", 2000, 8000]]
    assert by[(4000, 4500)]["cls"] == "cohort" and by[(4000, 4500)]["n"] == 9


def test_a_chain_of_small_neighbouring_deletions_is_not_one_cohort_gap():
    """Ten samples each lack 700 bp, each overlapping the next by 250: no position is lacked by
    more than two of them, so none of it is a stretch nobody reads."""
    profiles = {(f"S{i}", "REF"): profile(zero=[(1000 + 450 * i, 1699 + 450 * i)]) for i in range(10)}

    regions = run(profiles)[0]["regions"]

    assert all(r["cls"] in ("private", "shared") for r in regions)
    assert max(r["n"] for r in regions) <= 2


def test_a_carriers_separate_stretches_are_not_stretched_across_a_region():
    """A lacks 1000-1300 and 5000-5300, B lacks 1200-5100: A lacks 602 bp, not 4,301."""
    profiles = cohort(6)
    profiles[("A", "REF")] = profile(zero=[(1000, 1300), (5000, 5300)])
    profiles[("B", "REF")] = profile(zero=[(1200, 5100)])

    section, tracks, _ = run(profiles, min_len=200, nbins=10)

    for r in section["regions"]:
        for sid, s, e in r["samples"]:
            if sid == "A":
                assert e - s + 1 == 301
    assert sum(tracks[("A", "REF")]) * 1000 == pytest.approx(602, abs=2)


def test_a_region_names_the_genes_its_carriers_do_not_read():
    profiles = cohort(4, breadth=[1.0, 1.0])
    profiles[("S0", "REF")] = profile(zero=[(1001, 3000)], breadth=[0.0, 1.0])

    [r] = run(profiles)[0]["regions"]

    assert r["genes"] == ["katG"]


def test_a_sample_on_two_references_is_two_profiles():
    """Its deletion on one reference is compared with that reference's samples, not merged
    with its tables of the other."""
    profiles = cohort(3, ref="REF")
    profiles.update(cohort(3, ref="OTHER"))
    profiles[("S0", "OTHER")] = profile(zero=[(5001, 7000)], ref="OTHER")

    section, tracks, _ = run(profiles)

    [r] = section["regions"]
    assert (r["ref"], r["cls"], r["of"]) == ("OTHER", "private", 3)
    assert not any(tracks[("S0", "REF")]) and any(tracks[("S0", "OTHER")])


# ------------------------------------------------------------------------------------ genes


def test_a_gene_most_samples_read_and_one_does_not_is_lost_in_that_one():
    profiles = cohort(4, breadth=[0.99, 1.0])
    profiles[("S0", "REF")] = profile(breadth=[0.2, 1.0])

    [g] = run(profiles)[0]["genes_lost"]

    assert (g["gene"], g["samples"], g["of"]) == ("katG", [["S0", 0.2]], 4)


def test_a_gene_nobody_reads_is_not_lost_by_anyone():
    profiles = cohort(4, breadth=[0.0, 1.0])

    assert run(profiles)[0]["genes_lost"] == []


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
    p = profile(length=1000, callable_bins=cov._spread(
        [("chr", 1, 500, 1.0), ("chr", 501, 1000, 0.5)], [("chr", 1000)], nbins=2))

    dens = cov.snp_per_kb([2, 1], p, nbins=2)

    # The landscape's bins are 1-based (499 and 501 positions here), hence the approximation.
    assert dens == [pytest.approx(4.0, rel=0.01), pytest.approx(4.0, rel=0.01)], \
        "1 SNP in 250 callable bp is 4 per kb"


def test_snp_per_kb_is_empty_where_almost_nothing_could_be_called():
    p = profile(length=1000, callable_bins=[50.0])

    assert cov.snp_per_kb([3], p, nbins=1) == [None]


def test_the_deletion_track_is_the_share_of_each_bin_a_sample_lacks():
    profiles = cohort(4, length=1000)
    profiles[("S0", "REF")] = profile(zero=[(1, 250)], length=1000)

    _, tracks, _ = run(profiles, nbins=2)

    assert tracks[("S0", "REF")] == [pytest.approx(0.5, abs=0.01), 0.0]


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

    profiles, genes = cov.parse_depth_profiles(paths, min_len=200)

    assert sorted(profiles) == [("S0", "R"), ("S1", "R"), ("S2", "R"), ("S3", "R")]
    p = profiles[("S0", "R")]
    assert p["median_dp"] == 30 and p["contigs"] == [("chr", 1000)]
    assert p["zero"] == [("chr", 301, 700)] and len(p["callable"]) == 200
    assert genes["R"] == [("chr", 400, 600, "lost", "")] and list(p["genes"]) == [0.0]
    section, _, regions = cov.build_coverage(profiles, genes)
    assert [(r["cls"], r["genes"], r["samples"]) for r in regions] == [("private", ["lost"], [["S0", 301, 700]])]
    assert section["genes_lost"][0]["gene"] == "lost"

    out = tmp_path / "deletions.tsv"
    cov.write_deletions(out, regions)
    assert out.read_text().splitlines()[1].split("\t") == [
        "D1", "R", "chr", "301", "700", "400", "private", "1", "4", "S0:301-700", "lost"]


def test_a_thousand_samples_are_parsed_into_little_memory(tmp_path):
    """The windows are reduced to the landscape's bins as they are read, and gene breadths kept
    as one small array per sample next to one gene list per reference: a thousand samples used
    to take about 3 GB. What each further sample costs is what has to stay small."""
    import tracemalloc

    def tables(sample):
        win = tmp_path / f"{sample}.depth_windows.tsv"
        win.write_text(f"# sample={sample}\n# reference=R\n# median_dp=40\n# contigs=chr:4400000\n"
                       "contig\tstart\tend\tmean_dp\tzero_frac\tcallable_frac\n"
                       + "".join(f"chr\t{s + 1}\t{s + 1000}\t40.00\t0.0000\t1.0000\n" for s in range(0, 4_400_000, 1000)))
        genes = tmp_path / f"{sample}.gene_depth.tsv"
        genes.write_text(f"# sample={sample}\n# reference=R\n# median_dp=40\n# contigs=chr:4400000\n"
                         "contig\tstart\tend\tgene\tlocus_tag\tlength\tbreadth\tcallable\tmean_dp\trel_depth\n"
                         + "".join(f"chr\t{s}\t{s + 900}\tg{s}\t\t901\t1.0000\t1.0000\t40.00\t1.0000\n"
                                   for s in range(1, 4_000_000, 1000)))
        return [str(win), str(genes)]

    paths = [p for i in range(20) for p in tables(f"S{i}")]
    tracemalloc.start()
    profiles, _ = cov.parse_depth_profiles(paths)
    kept, _ = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    # Twenty samples of a 4.4 Mb genome with 4,000 genes kept about 72 MB before; the reference's
    # gene list, shared by all of them, is about 1 MB of what is left.
    assert kept < 3_000_000, f"twenty samples kept {kept} bytes"
    assert len(profiles[("S19", "R")]["genes"]) == 4000
