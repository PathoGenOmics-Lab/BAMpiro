"""Unit tests for bin/build_snp_matrix.py.

Aggregates the per-sample annotated VCFs of a run into one wide SNP matrix:
one row per site, two columns (AF and DP) per sample.
"""

from __future__ import annotations

import subprocess
import sys

import pytest

from conftest import load_script

bsm = load_script("build_snp_matrix")


# --------------------------------------------------------------------------
# _af_dp: allele frequency
# --------------------------------------------------------------------------

def test_af_dp_ad_is_preferred_over_ao_ro_and_gt():
    """AD wins even when AO/RO would give a different answer."""
    af, dp = bsm._af_dp("GT:AD:AO:RO", "0/1:10,30:1:99")
    assert af == pytest.approx(0.75)
    assert dp == 40


def test_af_dp_ad_sums_every_alternate_allele():
    af, dp = bsm._af_dp("GT:AD", "1/2:10,15,25")
    assert af == pytest.approx(0.8)      # (15 + 25) / 50
    assert dp == 50


def test_af_dp_ad_summing_to_zero_falls_through_without_dividing():
    """`sum(ad) > 0` is what keeps this out of a ZeroDivisionError. The AD branch is skipped
    and the GT dosage supplies the fraction, but AD still states a depth of zero."""
    af, dp = bsm._af_dp("GT:AD", "0/1:0,0")
    assert af == pytest.approx(0.5)      # from the GT dosage, not from AD
    assert dp == 0, "a measured depth of zero is not the same as an unknown depth"


def test_af_dp_single_element_ad_is_ignored_for_af_but_still_gives_a_depth():
    """A one-entry AD (no alternate depth) cannot give an allele fraction, so the GT dosage is
    used. The depth it does state must not be thrown away: the matrix would emit an empty DP
    cell for a site whose depth is right there in the record."""
    af, dp = bsm._af_dp("GT:AD", "1:5")
    assert af == pytest.approx(1.0)
    assert dp == 5


@pytest.mark.parametrize("ad", [".,.", ".", ""])
def test_af_dp_ad_of_only_missing_values_is_ignored(ad):
    af, dp = bsm._af_dp("GT:AD", f"1/1:{ad}")
    assert af == pytest.approx(1.0)
    assert dp is None


def test_af_dp_unparsable_ad_is_ignored():
    af, dp = bsm._af_dp("GT:AD", "0/1:a,b")
    assert af == pytest.approx(0.5)
    assert dp is None


def test_af_dp_ao_ro_when_ad_is_absent():
    af, dp = bsm._af_dp("GT:AO:RO", "0/1:8:2")
    assert af == pytest.approx(0.8)
    assert dp == 10


def test_af_dp_ao_ro_sums_multi_allelic_ao():
    af, dp = bsm._af_dp("GT:AO:RO", "1/2:5,5:10")
    assert af == pytest.approx(0.5)
    assert dp == 20


def test_af_dp_ro_missing_value_counts_as_zero():
    af, dp = bsm._af_dp("GT:AO:RO", "1/1:8:.")
    assert af == pytest.approx(1.0)
    assert dp == 8


def test_af_dp_ao_without_ro_is_ignored():
    """Both keys are required; AO alone falls through to the GT dosage."""
    af, dp = bsm._af_dp("GT:AO", "0/1:8")
    assert af == pytest.approx(0.5)
    assert dp is None


def test_af_dp_ao_plus_ro_of_zero_falls_through():
    af, dp = bsm._af_dp("GT:AO:RO", "1/1:0:0")
    assert af == pytest.approx(1.0)
    assert dp is None


def test_af_dp_unparsable_ao_is_ignored():
    af, dp = bsm._af_dp("GT:AO:RO", "0/1:x:2")
    assert af == pytest.approx(0.5)
    assert dp is None


@pytest.mark.parametrize("gt,expected", [
    ("1/1", 1.0),
    ("1|1", 1.0),
    ("0/1", 0.5),
    ("0|1", 0.5),
    ("1", 1.0),
    ("0", 0.0),
    ("1/2", 1.0),        # every non-reference allele counts
    ("0/0", 0.0),
    ("./1", 1.0),        # missing alleles are dropped before the dosage
    ("0/./1", 0.5),
])
def test_af_dp_gt_dosage(gt, expected):
    af, dp = bsm._af_dp("GT", gt)
    assert af == pytest.approx(expected)
    # PIN: the GT branch computes a fraction only - it never fills in a depth.
    assert dp is None


@pytest.mark.parametrize("gt", ["./.", ".", ".|."])
def test_af_dp_no_call_returns_none(gt):
    assert bsm._af_dp("GT", gt) == (None, None)


def test_af_dp_absent_gt_key_defaults_to_no_call():
    """No AD, no AO/RO and no GT at all: DP is still read, but AF is None and
    main() drops the record."""
    assert bsm._af_dp("DP", "30") == (None, 30)


# --------------------------------------------------------------------------
# _af_dp: depth source order
# --------------------------------------------------------------------------

def test_af_dp_format_dp_wins_over_ad_sum():
    _, dp = bsm._af_dp("GT:AD:DP", "0/1:10,30:99")
    assert dp == 99


def test_af_dp_format_dp_wins_over_ao_plus_ro():
    _, dp = bsm._af_dp("GT:AO:RO:DP", "0/1:8:2:99")
    assert dp == 99


@pytest.mark.parametrize("raw_dp", [".", "", "abc", "1.5"])
def test_af_dp_unusable_format_dp_falls_back_to_the_ad_sum(raw_dp):
    _, dp = bsm._af_dp("GT:DP:AD", f"0/1:{raw_dp}:10,30")
    assert dp == 40


def test_af_dp_ad_sum_is_used_when_dp_is_absent():
    _, dp = bsm._af_dp("GT:AD", "0/1:10,30")
    assert dp == 40


def test_af_dp_extra_format_keys_without_values_are_dropped():
    """zip() truncates to the shorter of FORMAT and the sample column."""
    af, dp = bsm._af_dp("GT:AD:DP", "0/1:10,30")
    assert af == pytest.approx(0.75)
    assert dp == 40


# --------------------------------------------------------------------------
# _ann
# --------------------------------------------------------------------------

_FULL_ANN = ("ANN=G|missense_variant|MODERATE|rpoB|GENE_ID|transcript|TX_ID|"
             "protein_coding|1/1|c.1349C>T|p.Ser450Leu|1349/3519|1349/3519|450/1172||")


def test_ann_extracts_gene_effect_and_protein_change():
    assert bsm._ann(f"DP=30;{_FULL_ANN};MQ=60") == ("rpoB", "missense_variant", "p.Ser450Leu")


@pytest.mark.parametrize("info", ["", "DP=30;MQ=60", "AF=0.5", "XANN=a|b|c|d"])
def test_ann_returns_three_empty_strings_without_ann(info):
    assert bsm._ann(info) == ("", "", "")


def test_ann_uses_the_first_entry_only():
    info = _FULL_ANN + ",T|synonymous_variant|LOW|katG|GID2|transcript|TX2|protein_coding|1/1|c.1A>G|p.Lys1Lys"
    assert bsm._ann(info) == ("rpoB", "missense_variant", "p.Ser450Leu")


@pytest.mark.parametrize("ann,expected", [
    ("ANN=G|missense_variant|MODERATE|rpoB", ("rpoB", "missense_variant", "")),
    ("ANN=G|missense_variant", ("", "missense_variant", "")),
    ("ANN=G", ("", "", "")),
    ("ANN=", ("", "", "")),
])
def test_ann_tolerates_truncated_entries(ann, expected):
    assert bsm._ann(ann) == expected


# --------------------------------------------------------------------------
# end to end
# --------------------------------------------------------------------------

def _vcf(path, sample, records):
    header = (
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1>\n"
        f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n"
    )
    path.write_text(header + "".join(r + "\n" for r in records))
    return path


def _row(contig, pos, ref, alt, info, fmt, val):
    return "\t".join([contig, str(pos), ".", ref, alt, "50", "PASS", info, fmt, val])


def test_end_to_end_matrix(tmp_path, repo_root):
    a = _vcf(tmp_path / "SA.vcf", "SA", [
        _row("chr1", 100, "A", "G", _FULL_ANN, "GT:AD:DP", "1/1:0,30:30"),
        _row("chr1", 300, "C", "T", ".", "GT:AD:DP", "0/1:10,10:20"),
        _row("chr1", 500, "AT", "A", ".", "GT:AD:DP", "1/1:0,30:30"),   # deletion, dropped
        _row("chr1", 700, "A", "GG", ".", "GT:AD:DP", "1/1:0,30:30"),   # insertion, dropped
        _row("chr2", 50, "G", "A", ".", "GT:AD:DP", "1/1:0,20:20"),
    ])
    b = _vcf(tmp_path / "SB.vcf", "SB", [
        _row("chr1", 90, "T", "C", ".", "GT:AO:RO", "0/1:6:2"),
        _row("chr1", 100, "A", "T", ".", "GT:AD:DP", "1/1:0,40:40"),
        _row("chr1", 300, "C", "T", ".", "GT", "./."),                  # no call, dropped
    ])
    out = tmp_path / "snp_matrix.tsv"

    proc = subprocess.run(
        [
            sys.executable, str(repo_root / "bin" / "build_snp_matrix.py"),
            "--vcfs", str(a), str(b), str(tmp_path / "absent.vcf"),
            "--reference", "REF1",
            "-o", str(out),
        ],
        cwd=tmp_path, capture_output=True, text=True, timeout=120,
    )
    assert proc.returncode == 0, proc.stderr
    assert "4 SNP site(s) x 2 sample(s)" in proc.stderr

    rows = [line.split("\t") for line in out.read_text().splitlines()]
    assert rows[0] == [
        "reference", "contig", "pos", "ref_allele", "alt_allele", "gene", "effect", "aa_change",
        "SA|AF", "SA|DP", "SB|AF", "SB|DP",
    ]

    # Only SNPs survive; sites sort by contig then NUMERIC position (90 < 100 < 300).
    assert [(r[1], r[2]) for r in rows[1:]] == [
        ("chr1", "90"), ("chr1", "100"), ("chr1", "300"), ("chr2", "50"),
    ]
    assert all(len(r) == len(rows[0]) for r in rows)

    # chr1:90 - called only in SB, so SA gets two EMPTY cells (not zeros).
    assert rows[1] == ["REF1", "chr1", "90", "T", "C", "", "", "", "", "", "0.7500", "8"]
    # chr1:100 - both samples; the ALT set is merged and sorted, the annotation
    # is taken from the first VCF that carried one.
    assert rows[2] == ["REF1", "chr1", "100", "A", "G,T", "rpoB", "missense_variant",
                       "p.Ser450Leu", "1.0000", "30", "1.0000", "40"]
    # chr1:300 - SB's './.' record was dropped, so SB is blank here too.
    assert rows[3] == ["REF1", "chr1", "300", "C", "T", "", "", "", "0.5000", "20", "", ""]
    assert rows[4] == ["REF1", "chr2", "50", "G", "A", "", "", "", "1.0000", "20", "", ""]


def test_end_to_end_reference_defaults_to_the_contig_name(tmp_path, repo_root):
    a = _vcf(tmp_path / "SA.vcf", "SA", [_row("chr1", 100, "A", "G", ".", "GT:DP", "1/1:30")])
    out = tmp_path / "snp_matrix.tsv"
    proc = subprocess.run(
        [sys.executable, str(repo_root / "bin" / "build_snp_matrix.py"),
         "--vcfs", str(a), "-o", str(out)],
        cwd=tmp_path, capture_output=True, text=True, timeout=120,
    )
    assert proc.returncode == 0, proc.stderr
    row = out.read_text().splitlines()[1].split("\t")
    assert row[0] == "chr1"
    assert row[8:] == ["1.0000", "30"]


def test_end_to_end_multi_allelic_alt_keeps_only_the_first_allele(tmp_path, repo_root):
    a = _vcf(tmp_path / "SA.vcf", "SA", [
        _row("chr1", 100, "A", "G,T", ".", "GT:AD:DP", "1/2:0,20,10:30"),
        _row("chr1", 200, "A", "GG,T", ".", "GT:AD:DP", "1/2:0,20,10:30"),   # first ALT is not a SNP
    ])
    out = tmp_path / "snp_matrix.tsv"
    proc = subprocess.run(
        [sys.executable, str(repo_root / "bin" / "build_snp_matrix.py"),
         "--vcfs", str(a), "-o", str(out)],
        cwd=tmp_path, capture_output=True, text=True, timeout=120,
    )
    assert proc.returncode == 0, proc.stderr
    rows = [line.split("\t") for line in out.read_text().splitlines()]
    assert len(rows) == 2
    assert rows[1][2] == "100"
    assert rows[1][4] == "G"


def test_end_to_end_skips_malformed_records(tmp_path, repo_root):
    """Records before #CHROM, records with fewer than 8 columns and records
    with a non-numeric POS are dropped; an unreadable input only warns."""
    p = tmp_path / "SA.vcf"
    p.write_text(
        "##fileformat=VCFv4.2\n"
        + _row("chr1", 10, "A", "G", ".", "GT:DP", "1/1:30") + "\n"   # before #CHROM
        + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSA\n"
        + "chr1\t20\t.\tA\tG\t50\n"                                    # too few columns
        + _row("chr1", "NOT_A_POS", "A", "G", ".", "GT:DP", "1/1:30") + "\n"
        + _row("chr1", 100, "A", "G", ".", "GT:DP", "1/1:30") + "\n"
    )
    unreadable = tmp_path / "SB.vcf"
    unreadable.mkdir()
    out = tmp_path / "snp_matrix.tsv"
    proc = subprocess.run(
        [sys.executable, str(repo_root / "bin" / "build_snp_matrix.py"),
         "--vcfs", str(p), str(unreadable), "-o", str(out)],
        cwd=tmp_path, capture_output=True, text=True, timeout=120,
    )
    assert proc.returncode == 0, proc.stderr
    assert "[snp_matrix] WARN could not read" in proc.stderr
    rows = out.read_text().splitlines()
    assert len(rows) == 2
    assert rows[1].split("\t")[2] == "100"


def test_end_to_end_with_no_sites_writes_the_header_only(tmp_path, repo_root):
    a = _vcf(tmp_path / "SA.vcf", "SA", [])
    out = tmp_path / "snp_matrix.tsv"
    proc = subprocess.run(
        [sys.executable, str(repo_root / "bin" / "build_snp_matrix.py"),
         "--vcfs", str(a), str(tmp_path / "absent.vcf"), "--reference", "REF1", "-o", str(out)],
        cwd=tmp_path, capture_output=True, text=True, timeout=120,
    )
    assert proc.returncode == 0, proc.stderr
    lines = out.read_text().splitlines()
    assert lines == ["\t".join(["reference", "contig", "pos", "ref_allele", "alt_allele",
                                "gene", "effect", "aa_change", "SA|AF", "SA|DP"])]


# --------------------------------------------------------------------------
# depths from the all-positions VCFs: absent is not unknown
# --------------------------------------------------------------------------
#
# A blank cell said both "this sample was read here and carries no alternate allele" and
# "nobody read this sample here". The all-positions VCF holds a record for every position of the
# sample's reference, so it can tell them apart.

def _allpos(path, sample, records, contigs=("chr1",)):
    header = ("##fileformat=VCFv4.2\n" + "".join(f"##contig=<ID={c}>\n" for c in contigs)
              + f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n")
    text = header + "".join(r + "\n" for r in records)
    if str(path).endswith(".gz"):
        import gzip
        with gzip.open(path, "wt") as fh:
            fh.write(text)
    else:
        path.write_text(text)
    return path


def _backbone(contig, pos, dp, ao=0):
    """A backbone record as bin/backbone_allpos.awk writes it: no ALT, depth and RO,AO."""
    gt = "0" if dp >= 30 else "./."
    return "\t".join([contig, str(pos), ".", "A", ".", ".", ".",
                      f"ADP={dp};WT=1;HET=0;HOM=0;NC=0", "GT:DP:AD", f"{gt}:{dp}:{dp - ao},{ao}"])


def _index(*keys):
    return bsm._site_index(sorted(keys)), len(keys)


def test_read_depths_reads_the_matrix_sites_and_nothing_else(tmp_path):
    p = _allpos(tmp_path / "SA.all.pos.vcf", "SA", [
        _backbone("chr1", 99, 50), _backbone("chr1", 100, 45), _backbone("chr1", 101, 50),
        _backbone("chr1", 300, 0),
    ])
    index, n = _index(("chr1", 100), ("chr1", 300), ("chr2", 50))

    sample, depths, calls = bsm._read_depths(str(p), index, n)

    assert sample == "SA"
    assert list(depths) == [45, 0, -1], "a site of another reference has no record, not depth 0"
    assert calls == {}


def test_read_depths_takes_a_variant_record_as_a_call(tmp_path):
    """Atomised, an MNP of the variant VCF is a single-base SNP here, which the matrix would
    otherwise report as absent at a site where the sample carries it."""
    p = _allpos(tmp_path / "SA.all.pos.vcf", "SA", [
        "\t".join(["chr1", "200", ".", "C", "T", "50", ".", ".", "GT:DP:AD:RO:AO", "1:40:2,38:2:38"]),
    ])
    index, n = _index(("chr1", 200))

    _, depths, calls = bsm._read_depths(str(p), index, n)

    assert list(depths) == [40]
    assert calls[0][0] == "T"
    assert calls[0][1] == pytest.approx(0.95)


def test_read_depths_of_a_bgzipped_file(tmp_path):
    p = _allpos(tmp_path / "SA.all.pos.vcf.gz", "SA", [_backbone("chr1", 100, 45)])
    index, n = _index(("chr1", 100))

    assert list(bsm._read_depths(str(p), index, n)[1]) == [45]


def test_read_depths_of_a_truncated_file_is_nothing_rather_than_half(tmp_path, capsys):
    """Half a sample's depths would read as the other half never having been sequenced."""
    p = _allpos(tmp_path / "SA.all.pos.vcf.gz", "SA",
                [_backbone("chr1", pos, 45) for pos in range(1, 20_000)])
    p.write_bytes(p.read_bytes()[:-2_000])
    index, n = _index(("chr1", 10), ("chr1", 19_000))

    assert bsm._read_depths(str(p), index, n) is None
    assert "could not read depths" in capsys.readouterr().err


def test_read_depths_skips_a_blank_line(tmp_path):
    p = tmp_path / "SA.all.pos.vcf"
    _allpos(p, "SA", [_backbone("chr1", 100, 45)])
    p.write_text(p.read_text() + "\n")
    index, n = _index(("chr1", 100))

    assert list(bsm._read_depths(str(p), index, n)[1]) == [45]


def _run(tmp_path, repo_root, *args):
    out = tmp_path / "snp_matrix.tsv"
    proc = subprocess.run(
        [sys.executable, str(repo_root / "bin" / "build_snp_matrix.py"), *map(str, args), "-o", str(out)],
        cwd=tmp_path, capture_output=True, text=True, timeout=120,
    )
    assert proc.returncode == 0, proc.stderr
    return {r.split("\t")[1] + ":" + r.split("\t")[2]: r.split("\t")
            for r in out.read_text().splitlines()[1:]}, out.read_text(), proc.stderr


def _two_reference_cohort(tmp_path):
    """SA and SB mapped against chr1, SC against chr2, SD with no all-positions VCF."""
    vcfs = [
        _vcf(tmp_path / "SA.vcf", "SA", [_row("chr1", 100, "A", "G", ".", "GT:AD:DP", "1:0,30:30"),
                                         _row("chr1", 300, "A", "G", ".", "GT:AD:DP", "1:0,30:30")]),
        # SB carries a dinucleotide change the caller wrote as one MNP, which the matrix skips.
        _vcf(tmp_path / "SB.vcf", "SB", [_row("chr1", 199, "GC", "AT", ".", "GT:AD:DP", "1:2,38:40")]),
        _vcf(tmp_path / "SC.vcf", "SC", [_row("chr2", 50, "G", "A", ".", "GT:AD:DP", "1:0,20:20")]),
        _vcf(tmp_path / "SD.vcf", "SD", [_row("chr1", 200, "C", "T", ".", "GT:AD:DP", "1:0,25:25")]),
    ]
    depth = [
        _allpos(tmp_path / "SA.all.pos.vcf.gz", "SA",
                [_backbone("chr1", 200, 38), _backbone("chr1", 300, 30)]),
        _allpos(tmp_path / "SB.all.pos.vcf.gz", "SB", [
            _backbone("chr1", 100, 45),
            "\t".join(["chr1", "200", ".", "C", "T", "50", ".", ".", "GT:DP:AD", "1:40:2,38"]),
            _backbone("chr1", 300, 0),
        ]),
        _allpos(tmp_path / "SC.all.pos.vcf.gz", "SC", [_backbone("chr2", 50, 20)], contigs=("chr2",)),
    ]
    return vcfs, depth


def test_end_to_end_a_sample_read_without_the_allele_is_zero_with_the_depth(tmp_path, repo_root):
    vcfs, depth = _two_reference_cohort(tmp_path)

    rows, _, err = _run(tmp_path, repo_root, "--vcfs", *vcfs, "--depth-vcfs", *depth)

    sa, sb, sc, sd = slice(8, 10), slice(10, 12), slice(12, 14), slice(14, 16)
    assert rows["chr1:100"][sb] == ["0", "45"], "read 45 times without the allele"
    assert rows["chr1:300"][sb] == ["", "0"], "no read: absence says nothing, so no 0"
    assert rows["chr1:100"][sc] == ["", ""], "SC was mapped against the other reference"
    assert rows["chr2:50"][sa] == ["", ""]
    assert rows["chr1:100"][sd] == ["", ""], "no all-positions VCF, so still unknown"
    assert rows["chr1:100"][sa] == ["1.0000", "30"], "a called cell is left as the caller wrote it"
    assert "depths from 3 of 4 sample(s)" in err
    assert "no depths for SD" in err


def test_end_to_end_a_call_spelled_as_a_longer_allele_is_not_denied(tmp_path, repo_root):
    """SB's MNP at 199-200 is not in the matrix as written, but its all-positions VCF holds the
    atomised SNP at 200: the cell takes its fraction, not a 0 the reads contradict."""
    vcfs, depth = _two_reference_cohort(tmp_path)

    rows, _, err = _run(tmp_path, repo_root, "--vcfs", *vcfs, "--depth-vcfs", *depth)

    assert rows["chr1:200"][10:12] == ["0.9500", "40"]
    assert rows["chr1:200"][4] == "T"
    assert rows["chr1:200"][8:10] == ["0", "38"]
    assert "1 call(s) recovered from an MNP or complex record" in err


def test_end_to_end_without_depths_the_cells_stay_blank(tmp_path, repo_root):
    vcfs, _ = _two_reference_cohort(tmp_path)

    rows, _, err = _run(tmp_path, repo_root, "--vcfs", *vcfs)

    assert rows["chr1:100"][10:12] == ["", ""]
    assert "depths from" not in err


def test_end_to_end_reading_the_depths_in_parallel_changes_nothing(tmp_path, repo_root):
    vcfs, depth = _two_reference_cohort(tmp_path)

    _, one, _ = _run(tmp_path, repo_root, "--vcfs", *vcfs, "--depth-vcfs", *depth, "--threads", "1")
    _, many, _ = _run(tmp_path, repo_root, "--vcfs", *vcfs, "--depth-vcfs", *depth, "--threads", "3")

    assert one == many


def test_end_to_end_an_unreadable_depth_file_leaves_that_sample_unknown(tmp_path, repo_root):
    vcfs, depth = _two_reference_cohort(tmp_path)
    depth[1].write_bytes(b"not gzip at all")

    rows, _, err = _run(tmp_path, repo_root, "--vcfs", *vcfs, "--depth-vcfs", *depth)

    assert "could not read depths" in err
    assert rows["chr1:100"][10:12] == ["", ""]
    assert rows["chr1:200"][8:10] == ["0", "38"], "the other samples are unaffected"
