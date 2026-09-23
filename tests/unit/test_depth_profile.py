"""Unit tests for bin/depth_profile.py: what one sample's reads cover, from its all-positions VCF."""

from __future__ import annotations

import gzip

import pytest

from conftest import load_script

dp = load_script("depth_profile")


def _backbone(contig, pos, depth):
    gt = "0" if depth >= 30 else "./."
    return "\t".join([contig, str(pos), ".", "A", ".", ".", ".", f"ADP={depth};WT=1;HET=0;HOM=0;NC=0",
                      "GT:DP:AD", f"{gt}:{depth}:{depth},0"])


def _allpos(path, depths, sample="S1", contigs=None, extra=()):
    """An all-positions VCF: `depths` maps contig -> list of depths from position 1 (None skips
    the record), `contigs` maps contig -> the length its ##contig line states."""
    contigs = contigs if contigs is not None else {c: len(v) for c, v in depths.items()}
    head = "##fileformat=VCFv4.2\n" + "".join(
        f"##contig=<ID={c},length={n}>\n" if n else f"##contig=<ID={c}>\n" for c, n in contigs.items())
    head += f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n"
    body = [_backbone(c, i + 1, d) for c, v in depths.items() for i, d in enumerate(v) if d is not None]
    text = head + "\n".join(list(body) + list(extra)) + "\n"
    if str(path).endswith(".gz"):
        with gzip.open(path, "wt") as fh:
            fh.write(text)
    else:
        path.write_text(text)
    return path


# ------------------------------------------------------------------------------------ reading


def test_read_depths_keeps_a_missing_record_apart_from_depth_zero(tmp_path):
    p = _allpos(tmp_path / "a.vcf", {"chr": [5, 0, None, 7]})

    sample, depth, known = dp.read_depths(str(p))

    assert sample == "S1"
    assert list(depth["chr"]) == [5, 0, 0, 7]
    assert list(known["chr"]) == [True, True, False, True], "no record is not the same as no reads"


def test_read_depths_reads_dp_wherever_the_format_puts_it(tmp_path):
    """A variant record keeps the caller's FORMAT, which need not put DP second."""
    variant = "\t".join(["chr", "2", ".", "C", "T", "50", ".", ".", "GT:AD:DP", "1:2,38:40"])
    p = _allpos(tmp_path / "a.vcf", {"chr": [5, None, 7]}, extra=[variant])

    _, depth, _ = dp.read_depths(str(p))

    assert list(depth["chr"]) == [5, 40, 7]


def test_read_depths_grows_a_contig_whose_length_is_not_stated(tmp_path):
    p = _allpos(tmp_path / "a.vcf.gz", {"chr": [3] * 1500}, contigs={"chr": 0})

    _, depth, known = dp.read_depths(str(p))

    assert len(depth["chr"]) == 1500
    assert known["chr"].all()


def test_read_depths_does_not_take_an_unstated_depth_for_zero(tmp_path):
    record = "\t".join(["chr", "2", ".", "C", ".", ".", ".", ".", "GT:DP:AD", "./.:.:0,0"])
    p = _allpos(tmp_path / "a.vcf", {"chr": [5, None, 7]}, extra=[record])

    _, _, known = dp.read_depths(str(p))

    assert not known["chr"][1]


# ------------------------------------------------------------------------------------ tables


def _arrays(values):
    import numpy as np
    d = np.array([0 if v is None else v for v in values], dtype=np.int32)
    return {"chr": d}, {"chr": np.array([v is not None for v in values])}


def test_windows_measure_only_the_positions_with_a_record():
    depth, known = _arrays([10, 0, None, 20] + [8] * 4)

    rows = dp.windows(depth, known, 4, min_dp=7)

    assert rows[0][:3] == ("chr", 1, 4)
    assert rows[0][3] == pytest.approx(10.0)            # (10 + 0 + 20) / 3
    assert rows[0][4] == pytest.approx(1 / 3)
    assert rows[0][5] == pytest.approx(2 / 3)           # 10 and 20 are above 7; 0 is not
    assert rows[1][5] == pytest.approx(1.0)


def test_a_window_without_records_has_no_values():
    depth, known = _arrays([None, None, 5, 5])

    assert dp.windows(depth, known, 2, 7)[0][3:] == (None, None, None)


def test_zero_runs_keep_the_stretches_long_enough_with_one_based_ends():
    depth, known = _arrays([5, 0, 0, 0, 5, 0, 5, 0, 0, 0])

    assert dp.zero_runs(depth, known, 3) == [("chr", 2, 4), ("chr", 8, 10)]


def test_a_missing_record_breaks_a_stretch_without_reads():
    depth, known = _arrays([0, 0, None, 0, 0])

    assert dp.zero_runs(depth, known, 3) == []


# ------------------------------------------------------------------------------------ genes


def _gff(path, rows):
    path.write_text("##gff-version 3\n" + "".join("\t".join(r) + "\n" for r in rows))
    return str(path)


def test_read_genes_prefers_gene_features_and_keeps_the_contig(tmp_path):
    p = _gff(tmp_path / "a.gff", [
        ["chr", ".", "gene", "1", "4", ".", "+", ".", "ID=g1;Name=dnaA;locus_tag=T_0001"],
        ["chr", ".", "CDS", "1", "4", ".", "+", ".", "ID=c1;Parent=g1"],
        ["plasmid", ".", "gene", "10", "20", ".", "+", ".", "ID=g2;locus_tag=T_0002"],
    ])

    assert dp.read_genes(p) == [("chr", 1, 4, "dnaA", "T_0001"), ("plasmid", 10, 20, "T_0002", "T_0002")]


def test_read_genes_falls_back_to_cds_when_there_are_no_genes(tmp_path):
    p = _gff(tmp_path / "a.gff", [["chr", ".", "CDS", "5", "9", ".", "+", ".", "ID=cds1;Name=katG"]])

    assert dp.read_genes(p) == [("chr", 5, 9, "katG", "")]


@pytest.mark.parametrize("name", ["NO_FILE", "NO_FILE_GFF", "absent.gff"])
def test_read_genes_without_a_gff_is_empty(tmp_path, name):
    p = tmp_path / name
    if name.startswith("NO_FILE"):
        p.write_text("")
    assert dp.read_genes(str(p)) == []


def test_gene_depths_measure_breadth_callable_depth_and_depth_against_the_median():
    depth, known = _arrays([0, 0, 10, 30, 20, 20])

    row = dp.gene_depths(depth, known, [("chr", 1, 4, "g", "")], min_dp=7, median=20.0)[0]

    assert row[5] == 4
    assert row[6] == pytest.approx(0.5)                 # two of four positions read at all
    assert row[7] == pytest.approx(0.5)
    assert row[8] == pytest.approx(10.0)
    assert row[9] == pytest.approx(0.5)


def test_gene_depths_of_a_gene_on_a_contig_the_sample_lacks_are_skipped():
    depth, known = _arrays([5, 5])

    assert dp.gene_depths(depth, known, [("plasmid", 1, 2, "g", "")], 7, 5.0) == []


# ------------------------------------------------------------------------------------ end to end


def test_end_to_end_writes_the_three_tables_with_their_provenance(tmp_path):
    vcf = _allpos(tmp_path / "S1.all.pos.vcf.gz", {"chr": [20] * 100 + [0] * 60 + [20] * 40})
    gff = _gff(tmp_path / "ref.gff", [["chr", ".", "gene", "101", "160", ".", "+", ".", "Name=lost"]])

    assert dp.main(["--vcf", str(vcf), "--gff", gff, "--reference", "REF", "--window", "50",
                    "--out-prefix", str(tmp_path / "S1")]) == 0

    win = (tmp_path / "S1.depth_windows.tsv").read_text().splitlines()
    assert "# sample=S1" in win and "# reference=REF" in win and "# median_dp=20" in win
    assert "# contigs=chr:200" in win
    rows = [r.split("\t") for r in win if not r.startswith("#")]
    assert rows[0] == ["contig", "start", "end", "mean_dp", "zero_frac", "callable_frac"]
    assert rows[3] == ["chr", "101", "150", "0.00", "1.0000", "0.0000"]
    zero = [r for r in (tmp_path / "S1.zero_depth.tsv").read_text().splitlines() if not r.startswith("#")]
    assert zero == ["contig\tstart\tend\tlength", "chr\t101\t160\t60"]
    genes = [r.split("\t") for r in (tmp_path / "S1.gene_depth.tsv").read_text().splitlines()
             if not r.startswith("#")]
    assert genes[1][3] == "lost" and genes[1][6] == "0.0000" and genes[1][9] == "0.0000"


def test_end_to_end_without_a_gff_writes_no_gene_table(tmp_path):
    vcf = _allpos(tmp_path / "S1.all.pos.vcf", {"chr": [20] * 10})

    assert dp.main(["--vcf", str(vcf), "--out-prefix", str(tmp_path / "S1")]) == 0

    assert not (tmp_path / "S1.gene_depth.tsv").exists()
    assert (tmp_path / "S1.zero_depth.tsv").exists()
