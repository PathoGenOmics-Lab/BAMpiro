"""bin/backbone_allpos.awk - one VCF record per reference position.

The backbone is what makes an all-sites VCF, and therefore a phylogenetic supermatrix, possible.
Its `count_bases` reads an mpileup base string, which is not one character per read: a read start
is '^' plus a mapping-quality character that must not be counted, and an indel is a sign, a length,
and that many bases belonging to the previous position. Miscount any of that and the allelic depth
the consensus relies on is quietly wrong.

Both branches of CALL_BACKBONE carried their own copy of this until it was pulled out, and neither
had any coverage: a stub run skips the script block.
"""

from __future__ import annotations

import subprocess

import pytest

from conftest import BIN

AWK = BIN / "backbone_allpos.awk"


def run(rows, mincov=30, gapdp=7, baq=0, tmp_path=None):
    """Feed mpileup-shaped rows through the awk program and return the VCF fields per line."""
    src = tmp_path / "pileup.txt"
    src.write_text("".join(rows))
    cmd = ["awk", "-v", f"MINCOV={mincov}", "-v", f"BAQ={baq}", "-v", "OFS=\t"]
    if baq:
        cmd += ["-v", f"GAPDP={gapdp}"]
    cmd += ["-f", str(AWK), str(src)]
    out = subprocess.run(cmd, capture_output=True, text=True, check=True)
    return [line.split("\t") for line in out.stdout.splitlines()]


def pileup(bases, depth=None, pos=100, ref="A", baq_depth=None):
    """One mpileup row. Column 10 carries the BAQ-adjusted depth when the two pileups are pasted."""
    dp = len(bases) if depth is None else depth
    row = f"c\t{pos}\t{ref}\t{dp}\t{bases}\tIIII"
    if baq_depth is not None:
        row += f"\tc\t{pos}\t{ref}\t{baq_depth}\t{bases}\tIIII"
    return row + "\n"


def ad_of(fields):
    """The ref,alt observation counts out of the GT:DP:AD sample column."""
    return fields[9].split(":")[2]


def info_of(fields):
    return dict(kv.split("=", 1) for kv in fields[7].split(";"))


class TestBaseCounting:
    def test_matches_to_the_reference_count_as_reference_observations(self, tmp_path):
        [f] = run([pileup("....,,,,")], tmp_path=tmp_path)
        assert ad_of(f) == "8,0"

    def test_mismatches_count_as_alt_observations_in_either_case(self, tmp_path):
        [f] = run([pileup("..AAtt")], tmp_path=tmp_path)
        assert ad_of(f) == "2,4"

    def test_a_read_start_consumes_its_mapping_quality_character(self, tmp_path):
        """'^I.' is ONE reference observation: the 'I' is a quality score, not a base."""
        [f] = run([pileup("^I.^I.^I.")], tmp_path=tmp_path)
        assert ad_of(f) == "3,0", "the mapping-quality character was counted as a base"

    def test_a_read_start_whose_quality_char_looks_like_a_base_is_not_counted(self, tmp_path):
        """The quality character can legitimately be 'A' or 'T'; it must still be skipped."""
        [f] = run([pileup("^A.^T.")], tmp_path=tmp_path)
        assert ad_of(f) == "2,0", "a quality character was mistaken for an alt observation"

    def test_a_read_end_marker_is_not_a_base(self, tmp_path):
        [f] = run([pileup("..$..$")], tmp_path=tmp_path)
        assert ad_of(f) == "4,0"

    def test_an_insertion_and_its_bases_are_skipped(self, tmp_path):
        """'+2AG' belongs to the previous position, so neither the length nor the bases count here."""
        [f] = run([pileup(".+2AG.")], tmp_path=tmp_path)
        assert ad_of(f) == "2,0", "inserted bases were counted at this position"

    def test_a_multi_digit_indel_length_is_parsed_whole(self, tmp_path):
        """A 12-base insertion must consume 12 bases, not 1 then treat '2' as a length."""
        [f] = run([pileup(".+12ACGTACGTACGT.")], tmp_path=tmp_path)
        assert ad_of(f) == "2,0"

    def test_a_deletion_and_its_bases_are_skipped(self, tmp_path):
        [f] = run([pileup(".-3ACG.")], tmp_path=tmp_path)
        assert ad_of(f) == "2,0"

    def test_a_spanning_deletion_placeholder_is_ignored(self, tmp_path):
        """'*' marks a position inside a deletion: neither reference nor alt evidence."""
        [f] = run([pileup("..**")], tmp_path=tmp_path)
        assert ad_of(f) == "2,0"

    def test_an_n_base_is_ignored(self, tmp_path):
        [f] = run([pileup("..NN")], tmp_path=tmp_path)
        assert ad_of(f) == "2,0"

    def test_a_position_with_no_reads_counts_nothing(self, tmp_path):
        [f] = run([pileup("*", depth=0)], tmp_path=tmp_path)
        assert ad_of(f) == "0,0"


class TestRecordShape:
    def test_the_record_is_a_non_variant_site(self, tmp_path):
        """The backbone asserts coverage, it does not call variants: ALT stays '.'."""
        [f] = run([pileup("." * 40, pos=4242, ref="G")], tmp_path=tmp_path)
        assert f[0] == "c" and f[1] == "4242" and f[3] == "G"
        assert f[4] == "." and f[5] == "."
        assert f[8] == "GT:DP:AD"

    def test_depth_comes_from_the_pileup_column_not_the_base_string(self, tmp_path):
        [f] = run([pileup("....", depth=99)], tmp_path=tmp_path)
        assert info_of(f)["ADP"] == "99"
        assert f[9].split(":")[1] == "99"


class TestCoverageCall:
    def test_a_well_covered_position_is_a_confident_reference_call(self, tmp_path):
        [f] = run([pileup("." * 40, depth=40)], mincov=30, tmp_path=tmp_path)
        assert f[9].startswith("0:")
        assert info_of(f)["WT"] == "1" and info_of(f)["NC"] == "0"

    def test_a_shallow_position_is_a_no_call(self, tmp_path):
        [f] = run([pileup("." * 5, depth=5)], mincov=30, tmp_path=tmp_path)
        assert f[9].startswith("./.:")
        assert info_of(f)["WT"] == "0" and info_of(f)["NC"] == "1"

    def test_the_threshold_is_inclusive(self, tmp_path):
        at = run([pileup(".", depth=30)], mincov=30, tmp_path=tmp_path)[0]
        below = run([pileup(".", depth=29)], mincov=30, tmp_path=tmp_path)[0]
        assert at[9].startswith("0:")
        assert below[9].startswith("./.:")

    def test_the_backbone_never_calls_a_heterozygote(self, tmp_path):
        [f] = run([pileup("..AATT", depth=40)], tmp_path=tmp_path)
        assert info_of(f)["HET"] == "0" and info_of(f)["HOM"] == "0"


class TestBaqDropout:
    """A well-covered site whose depth collapses only under BAQ is an indel-adjacent artefact.

    Flagging it lets the consensus mask it with X rather than leave a gap, which is the difference
    between "we could not read here" and "there is nothing here".
    """

    def test_a_collapsed_site_is_flagged(self, tmp_path):
        [f] = run([pileup("." * 40, depth=40, baq_depth=3)], gapdp=7, baq=1, tmp_path=tmp_path)
        assert f[6] == "baq_dropout"

    def test_a_site_that_holds_up_under_baq_is_not_flagged(self, tmp_path):
        [f] = run([pileup("." * 40, depth=40, baq_depth=38)], gapdp=7, baq=1, tmp_path=tmp_path)
        assert f[6] == "."

    def test_a_genuinely_shallow_site_is_not_flagged(self, tmp_path):
        """Depth at or below the gap threshold to begin with is a gap, not a BAQ artefact."""
        [f] = run([pileup("..", depth=2, baq_depth=1)], gapdp=7, baq=1, tmp_path=tmp_path)
        assert f[6] == "."

    def test_without_the_baq_column_nothing_is_ever_flagged(self, tmp_path):
        """The single-pileup branch has nothing to compare against, so the filter stays empty."""
        [f] = run([pileup("." * 40, depth=40)], baq=0, tmp_path=tmp_path)
        assert f[6] == "."

    def test_the_depth_reported_is_the_real_one_not_the_baq_one(self, tmp_path):
        """ADP must stay the no-BAQ depth; BAQ is only used to decide the flag."""
        [f] = run([pileup("." * 40, depth=40, baq_depth=3)], gapdp=7, baq=1, tmp_path=tmp_path)
        assert info_of(f)["ADP"] == "40"


def test_both_backbone_branches_share_this_one_program():
    """The BAQ and non-BAQ branches used to carry separate copies of count_bases."""
    variants = (BIN.parent / "modules" / "variants.nf").read_text()
    assert variants.count("backbone_allpos.awk") == 2
    assert "count_bases" not in variants, "the pileup parser is back inline in the module"


@pytest.mark.parametrize("bases,expected", [
    ("^].^].,,", "4,0"),                 # read starts with ']' as the quality char
    (".$.+1A.", "3,0"),                  # read end, then an insertion, then a match
    ("A^AA", "0,2"),                     # alt, then a read start whose quality is 'A', then alt
    ("", "0,0"),                         # empty base string
])
def test_mixed_pileup_strings(bases, expected, tmp_path):
    [f] = run([pileup(bases, depth=10)], tmp_path=tmp_path)
    assert ad_of(f) == expected
