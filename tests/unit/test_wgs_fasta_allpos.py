"""Unit tests for bin/WGS_fasta_allpos.py, the all.pos consensus builder.

The script is a standalone file under bin/, so it is loaded through the shared
`load_script` helper in tests/conftest.py rather than imported as a package.

Every assertion here pins the CURRENT, observed behaviour of the script. Where
that behaviour looks questionable the test carries a `QUIRK:` comment naming the
issue; those tests are deliberately written to pass against the code as it is.
"""

from __future__ import annotations

import dataclasses
import gzip

import pytest

from conftest import load_script

wgs = load_script("WGS_fasta_allpos")


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #

# The argparse defaults from main(), so a test only states the knob it exercises.
DEFAULT_OPTS = {
    "mask_char": "X",
    "nocall_char": "-",
    "min_dp": 0,
    "ref_min_dp": 10,
    "max_ref_altfrac": 0.10,
    "max_ref_min_alt": 2,
    "max_ref_del_frac": 0.5,
}

BUILD_DEFAULTS = {
    "exclude_path": None,
    "mask_sites_path": None,
    "special_csv_path": None,
    "wrap": 60,
    "strict": False,
    **DEFAULT_OPTS,
}

VCF_HEADER = (
    "##fileformat=VCFv4.2\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
)


def make_rec(chrom="chr1", pos=1, ref="A", alt=".", qual="60", flt="PASS",
             info=".", fmt="GT:DP:AD", sample="0/0:30:30,0"):
    """A VcfRecord that, with the default options, yields a confident reference call."""
    return wgs.VcfRecord(chrom=chrom, pos=pos, ref=ref, alt=alt, qual=qual,
                         flt=flt, info=info, fmt=fmt, sample=sample)


def call_consensus(recs, ref_seqs=None, interval_mask=None, point_mask=None,
                   csv_mask=None, **overrides):
    opts = dict(DEFAULT_OPTS)
    opts.update(overrides)
    return wgs.consensus_for_position(
        recs,
        {} if ref_seqs is None else ref_seqs,
        interval_mask if interval_mask is not None else wgs.IntervalMasker(None),
        point_mask if point_mask is not None else wgs.PointMasker(None),
        csv_mask if csv_mask is not None else wgs.CsvMasker(None),
        **opts,
    )


def vcf_line(chrom, pos, ref, alt, flt="PASS", info=".", fmt="GT:DP:AD",
             sample="0/0:30:30,0"):
    return "\t".join([chrom, str(pos), ".", ref, alt, "60", flt, info, fmt, sample])


def write_text(path, text, gz=False):
    if gz:
        with gzip.open(path, "wt", encoding="utf-8") as fh:
            fh.write(text)
    else:
        path.write_text(text, encoding="utf-8")
    return str(path)


def write_vcf(path, lines, gz=False):
    return write_text(path, VCF_HEADER + "".join(f"{ln}\n" for ln in lines), gz=gz)


def write_fasta(path, seqs, wrap=60, gz=False):
    out = []
    for name, seq in seqs.items():
        out.append(f">{name}\n")
        for i in range(0, max(len(seq), 1), wrap):
            out.append(seq[i:i + wrap] + "\n")
    return write_text(path, "".join(out), gz=gz)


def read_fasta(path):
    """Minimal FASTA parser for the produced output (independent of FastaReader)."""
    seqs = {}
    name = None
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            name = line[1:]
            seqs[name] = ""
        elif name is not None:
            seqs[name] += line
    return seqs


def run_build(tmp_path, vcf_lines, ref_seqs, **overrides):
    kwargs = dict(BUILD_DEFAULTS)
    kwargs.update(overrides)
    out = tmp_path / "consensus.fa"
    wgs.build_consensus(
        vcf_path=write_vcf(tmp_path / "in.vcf", vcf_lines),
        reference_path=write_fasta(tmp_path / "ref.fa", ref_seqs),
        out_path=str(out),
        **kwargs,
    )
    return out


# --------------------------------------------------------------------------- #
# iupac_for_bases
# --------------------------------------------------------------------------- #

class TestIupacForBases:
    @pytest.mark.parametrize("bases,expected", [
        (["A", "C", "G"], "V"),
        (["A", "C", "T"], "H"),
        (["A", "G", "T"], "D"),
        (["C", "G", "T"], "B"),
        (["A", "G"], "R"),
        (["C", "T"], "Y"),
        (["G", "C"], "S"),
        (["A", "T"], "W"),
        (["G", "T"], "K"),
        (["A", "C"], "M"),
    ])
    def test_ambiguity_codes(self, bases, expected):
        assert wgs.iupac_for_bases(bases) == expected

    def test_all_four_bases_collapse_to_n(self):
        assert wgs.iupac_for_bases(["A", "C", "G", "T"]) == "N"

    def test_empty_input_is_n(self):
        assert wgs.iupac_for_bases([]) == "N"

    def test_only_unknown_characters_is_n(self):
        assert wgs.iupac_for_bases(["X", "Z", "-", "N"]) == "N", \
            "characters outside ACGT contribute nothing, so the set is empty -> N"

    @pytest.mark.parametrize("bases,expected", [
        (["A"], "A"),
        (["t"], "T"),
        (["G", "G", "G"], "G"),
    ])
    def test_single_distinct_base_returns_itself(self, bases, expected):
        assert wgs.iupac_for_bases(bases) == expected

    def test_non_dna_members_are_dropped_before_lookup(self):
        # 'N' and '' are discarded, leaving {A, G} -> R rather than N.
        assert wgs.iupac_for_bases(["A", "", "N", "g"]) == "R"

    def test_accepts_any_iterable_of_characters(self):
        assert wgs.iupac_for_bases("ACG") == "V"


# --------------------------------------------------------------------------- #
# parse_info_field
# --------------------------------------------------------------------------- #

class TestParseInfoField:
    @pytest.mark.parametrize("info", [".", ""])
    def test_missing_info_is_empty_dict(self, info):
        assert wgs.parse_info_field(info) == {}

    def test_key_value_pairs(self):
        assert wgs.parse_info_field("DP=30;ADP=28") == {"DP": "30", "ADP": "28"}

    def test_flag_without_value_becomes_true(self):
        assert wgs.parse_info_field("INDEL;DP=30") == {"INDEL": "true", "DP": "30"}

    def test_items_are_stripped_and_blanks_skipped(self):
        assert wgs.parse_info_field(" DP=30 ;; INDEL ") == {"DP": "30", "INDEL": "true"}

    def test_only_the_first_equals_sign_splits(self):
        assert wgs.parse_info_field("ANN=a=b=c") == {"ANN": "a=b=c"}


# --------------------------------------------------------------------------- #
# parse_format_sample
# --------------------------------------------------------------------------- #

class TestParseFormatSample:
    @pytest.mark.parametrize("fmt", ["", "."])
    def test_missing_format_is_empty_dict(self, fmt):
        assert wgs.parse_format_sample(fmt, "0/1:30") == {}

    def test_zips_keys_to_values(self):
        assert wgs.parse_format_sample("GT:DP:AD", "0/1:30:15,15") == {
            "GT": "0/1", "DP": "30", "AD": "15,15",
        }

    def test_short_sample_column_pads_with_dot(self):
        assert wgs.parse_format_sample("GT:DP:AD:RO", "0/1:30") == {
            "GT": "0/1", "DP": "30", "AD": ".", "RO": ".",
        }, "keys past the end of the sample column must be filled with '.'"

    def test_empty_sample_column_pads_every_key(self):
        assert wgs.parse_format_sample("GT:DP", "") == {"GT": ".", "DP": "."}

    def test_extra_sample_values_are_ignored(self):
        assert wgs.parse_format_sample("GT", "0/1:30:99") == {"GT": "0/1"}


# --------------------------------------------------------------------------- #
# get_dp
# --------------------------------------------------------------------------- #

class TestGetDp:
    def test_format_dp_wins_over_info(self):
        assert wgs.get_dp("GT:DP", "0/0:42", "DP=7;ADP=8") == 42

    def test_falls_back_to_info_dp_when_format_dp_is_missing(self):
        assert wgs.get_dp("GT", "0/0", "DP=7") == 7

    @pytest.mark.parametrize("dp_value", [".", ""])
    def test_falls_back_to_info_dp_when_format_dp_is_blank(self, dp_value):
        assert wgs.get_dp("GT:DP", f"0/0:{dp_value}", "DP=7") == 7

    def test_falls_back_to_info_dp_when_format_dp_is_not_numeric(self):
        assert wgs.get_dp("GT:DP", "0/0:high", "DP=7") == 7

    def test_falls_back_to_info_adp_when_info_dp_is_absent(self):
        assert wgs.get_dp("GT", "0/0", "ADP=9") == 9

    def test_falls_back_to_info_adp_when_info_dp_is_not_numeric(self):
        assert wgs.get_dp("GT", "0/0", "DP=deep;ADP=9") == 9

    def test_zero_when_no_depth_anywhere(self):
        assert wgs.get_dp("GT", "0/0", ".") == 0

    def test_float_depth_is_truncated_to_int(self):
        assert wgs.get_dp("GT:DP", "0/0:12.9", ".") == 12


# --------------------------------------------------------------------------- #
# get_gt
# --------------------------------------------------------------------------- #

class TestGetGt:
    @pytest.mark.parametrize("gt", [".", "./.", ".|.", "./", ".|", ""])
    def test_missing_genotypes_are_none(self, gt):
        assert wgs.get_gt("GT:DP", f"{gt}:30") is None

    def test_absent_gt_key_is_none(self):
        assert wgs.get_gt("DP:AD", "30:30,0") is None

    @pytest.mark.parametrize("gt", ["0/0", "0/1", "1|2", "1"])
    def test_present_genotypes_are_returned_verbatim(self, gt):
        assert wgs.get_gt("GT:DP", f"{gt}:30") == gt

    def test_half_missing_genotype_is_not_filtered_here(self):
        # QUIRK: '0/.' is not in the missing-genotype set, so get_gt passes it
        # through; parse_gt_alleles is what ultimately rejects it.
        assert wgs.get_gt("GT", "0/.") == "0/."
        assert wgs.parse_gt_alleles("0/.") is None


# --------------------------------------------------------------------------- #
# get_ro_ao
# --------------------------------------------------------------------------- #

class TestGetRoAo:
    def test_reads_ro_and_ao(self):
        assert wgs.get_ro_ao("GT:RO:AO", "0/1:80:20") == (80, 20)

    def test_ao_with_multiple_alts_uses_only_the_first_value(self):
        assert wgs.get_ro_ao("GT:RO:AO", "1/2:5:7,9") == (5, 7)

    def test_falls_back_to_ad_when_ro_and_ao_are_absent(self):
        assert wgs.get_ro_ao("GT:AD", "0/1:12,8") == (12, 8)

    def test_ad_alt_counts_are_summed(self):
        assert wgs.get_ro_ao("GT:AD", "1/2:4,3,2") == (4, 5), \
            "ro is AD[0] and ao is the sum of every remaining AD field"

    def test_ad_with_a_single_field_gives_zero_alt_reads(self):
        assert wgs.get_ro_ao("GT:AD", "0/0:30") == (30, 0)

    @pytest.mark.parametrize("fmt,sample", [
        ("GT:DP", "0/0:30"),
        ("GT:RO:AO:AD", "0/0:.:.:."),
    ])
    def test_returns_none_none_when_nothing_is_available(self, fmt, sample):
        assert wgs.get_ro_ao(fmt, sample) == (None, None)

    def test_ad_fallback_overrides_a_present_ro_when_ao_is_missing(self):
        # QUIRK: the AD fallback triggers when EITHER RO or AO is missing and it
        # reassigns BOTH, so the explicit RO=99 is silently discarded.
        assert wgs.get_ro_ao("GT:RO:AD", "0/1:99:12,8") == (12, 8)

    def test_ad_fallback_rejects_float_formatted_counts(self):
        # QUIRK: RO/AO go through int(float(...)) but the AD fallback uses a bare
        # int(), so a float-formatted AD is dropped entirely.
        assert wgs.get_ro_ao("GT:AD", "0/0:30.0,0") == (None, None)

    def test_malformed_ad_can_leave_a_half_assigned_result(self):
        # BUG (pinned, not fixed): ro is assigned before the ao sum raises, and the
        # except branch only passes, so a malformed AD returns (ro, None) instead of
        # (None, None). Callers guard with `ro is not None and ao is not None`, so
        # today this only degrades to DP-based gating rather than corrupting output.
        assert wgs.get_ro_ao("GT:AD", "0/1:5,x") == (5, None)


# --------------------------------------------------------------------------- #
# parse_gt_alleles
# --------------------------------------------------------------------------- #

class TestParseGtAlleles:
    @pytest.mark.parametrize("gt,expected", [
        ("1", [1]),
        ("0", [0]),
        ("0/0", [0, 0]),
        ("0/1", [0, 1]),
        ("0|1", [0, 1]),
        ("1/2/3", [1, 2, 3]),
        (" 0 / 1 ", [0, 1]),
    ])
    def test_parses_allele_indices(self, gt, expected):
        assert wgs.parse_gt_alleles(gt) == expected

    @pytest.mark.parametrize("gt", ["", ".", "./.", ".|.", "0/.", "./1", "A/T", "0/1x"])
    def test_missing_or_non_integer_genotypes_are_none(self, gt):
        assert wgs.parse_gt_alleles(gt) is None


# --------------------------------------------------------------------------- #
# parse_vcf_record
# --------------------------------------------------------------------------- #

class TestParseVcfRecord:
    def test_parses_a_ten_column_record(self):
        rec = wgs.parse_vcf_record(
            "chr1\t42\trs1\tA\tG\t60\tPASS\tDP=30\tGT:DP\t0/1:30\n"
        )
        assert (rec.chrom, rec.pos, rec.ref, rec.alt) == ("chr1", 42, "A", "G")
        assert (rec.qual, rec.flt, rec.info) == ("60", "PASS", "DP=30")
        assert (rec.fmt, rec.sample) == ("GT:DP", "0/1:30")

    def test_fewer_than_ten_columns_is_none(self):
        assert wgs.parse_vcf_record("chr1\t42\trs1\tA\tG\t60\tPASS\tDP=30\tGT:DP") is None

    def test_non_integer_position_is_none(self):
        assert wgs.parse_vcf_record(
            "chr1\tNA\trs1\tA\tG\t60\tPASS\t.\tGT\t0/1"
        ) is None

    def test_space_separated_line_is_none(self):
        # The split is tab-only, so a space-delimited line collapses to one field.
        assert wgs.parse_vcf_record("chr1 42 rs1 A G 60 PASS . GT 0/1") is None

    def test_extra_sample_columns_are_ignored(self):
        rec = wgs.parse_vcf_record(
            "chr1\t7\t.\tA\tG\t60\tPASS\t.\tGT\t0/1\t1/1\t0/0\n"
        )
        assert rec.sample == "0/1", "only the first sample column is read"

    def test_record_is_frozen(self):
        rec = wgs.parse_vcf_record("chr1\t7\t.\tA\tG\t60\tPASS\t.\tGT\t0/1")
        with pytest.raises(dataclasses.FrozenInstanceError):
            rec.pos = 8


# --------------------------------------------------------------------------- #
# allele_bases_from_record
# --------------------------------------------------------------------------- #

class TestAlleleBasesFromRecord:
    def test_homozygous_alt_snp(self):
        rec = make_rec(ref="A", alt="T", sample="1/1:30:0,30")
        assert wgs.allele_bases_from_record(rec, "A") == {"T"}

    def test_heterozygous_snp_includes_the_reference_base(self):
        rec = make_rec(ref="A", alt="T", sample="0/1:30:15,15")
        assert wgs.allele_bases_from_record(rec, "A") == {"A", "T"}

    def test_ignore_ref_drops_the_reference_allele(self):
        rec = make_rec(ref="A", alt="T", sample="0/1:30:15,15")
        assert wgs.allele_bases_from_record(rec, "A", ignore_ref=True) == {"T"}

    def test_multiallelic_genotype_collects_both_alts(self):
        rec = make_rec(ref="A", alt="C,G", sample="1/2:30:0,15,15")
        assert wgs.allele_bases_from_record(rec, "A") == {"C", "G"}

    def test_monomorphic_record_yields_the_reference_base(self):
        rec = make_rec(ref="A", alt=".", sample="0/0:30:30,0")
        assert wgs.allele_bases_from_record(rec, "A") == {"A"}

    def test_monomorphic_record_with_ignore_ref_signals_mask(self):
        # Nothing is left to report once the reference allele is dropped.
        rec = make_rec(ref="A", alt=".", sample="0/0:30:30,0")
        assert wgs.allele_bases_from_record(rec, "A", ignore_ref=True) is None

    def test_monomorphic_record_with_non_dna_reference_signals_mask(self):
        rec = make_rec(ref="N", alt=".", sample="0/0:30:30,0")
        assert wgs.allele_bases_from_record(rec, "N") is None

    def test_reference_allele_is_skipped_when_ref_base_is_not_dna(self):
        rec = make_rec(ref="N", alt="T", sample="0/1:30:15,15")
        assert wgs.allele_bases_from_record(rec, "N") == {"T"}

    @pytest.mark.parametrize("alt,sample", [
        ("AT", "1/1:30:0,30"),      # insertion
        ("ATG", "1/1:30:0,30"),     # MNP
        ("<DEL>", "1/1:30:0,30"),   # symbolic allele
    ])
    def test_non_snp_allele_signals_mask(self, alt, sample):
        rec = make_rec(ref="A", alt=alt, sample=sample)
        assert wgs.allele_bases_from_record(rec, "A") is None, \
            "a non-SNP allele must return None, the mask signal"

    def test_allele_index_out_of_range_signals_mask(self):
        rec = make_rec(ref="A", alt="T", sample="2/2:30:0,30")
        assert wgs.allele_bases_from_record(rec, "A") is None

    def test_negative_allele_index_signals_mask(self):
        rec = make_rec(ref="A", alt="T", sample="-1/-1:30:0,30")
        assert wgs.allele_bases_from_record(rec, "A") is None

    def test_missing_genotype_signals_mask(self):
        rec = make_rec(ref="A", alt="T", sample="./.:30:15,15")
        assert wgs.allele_bases_from_record(rec, "A") is None

    def test_spanning_deletion_star_allele_is_skipped(self):
        rec = make_rec(ref="A", alt="*,T", sample="1/2:30:0,15,15")
        assert wgs.allele_bases_from_record(rec, "A") == {"T"}

    def test_only_star_allele_signals_mask(self):
        rec = make_rec(ref="A", alt="*", sample="1/1:30:0,30")
        assert wgs.allele_bases_from_record(rec, "A") is None

    def test_alt_alleles_are_uppercased(self):
        rec = make_rec(ref="A", alt="t", sample="1/1:30:0,30")
        assert wgs.allele_bases_from_record(rec, "A") == {"T"}


# --------------------------------------------------------------------------- #
# consensus_for_position: decision order and every branch
# --------------------------------------------------------------------------- #

def interval_masker(tmp_path, body):
    return wgs.IntervalMasker(write_text(tmp_path / "exclude.tsv", body))


def point_masker(tmp_path, body):
    return wgs.PointMasker(write_text(tmp_path / "mask_sites.tsv", body))


def csv_masker(tmp_path, body, contigs=None):
    return wgs.CsvMasker(write_text(tmp_path / "special.csv", body), contigs)


class TestConsensusMasking:
    def test_interval_mask_returns_mask_char(self, tmp_path):
        rec = make_rec(pos=10)
        assert call_consensus([rec], interval_mask=interval_masker(tmp_path, "chr1\t5\t15\n")) == "X"

    def test_point_mask_returns_mask_char(self, tmp_path):
        rec = make_rec(pos=10)
        assert call_consensus([rec], point_mask=point_masker(tmp_path, "chr1\t10\n")) == "X"

    def test_csv_mask_returns_mask_char(self, tmp_path):
        body = "rv_position,pal,homopolymer,GCrich,repetitive,blindspot\n10,0,0,0,1,0\n"
        rec = make_rec(pos=10)
        assert call_consensus([rec], csv_mask=csv_masker(tmp_path, body, ["chr1"])) == "X"

    def test_mask_char_is_configurable(self, tmp_path):
        rec = make_rec(pos=10)
        masker = point_masker(tmp_path, "chr1\t10\n")
        assert call_consensus([rec], point_mask=masker, mask_char="?") == "?"

    def test_masking_beats_the_depth_gate(self, tmp_path):
        # Masks are checked before the DP gate, so a zero-depth masked position is
        # X and not the no-call char. (Mask vs FILTER ordering is not observable:
        # both branches return mask_char.)
        rec = make_rec(pos=10, fmt="GT:DP", sample="0/0:0")
        masker = point_masker(tmp_path, "chr1\t10\n")
        assert call_consensus([rec], point_mask=masker, nocall_char="-") == "X"

    def test_position_outside_the_mask_is_not_masked(self, tmp_path):
        rec = make_rec(pos=20, ref="A")
        assert call_consensus([rec], interval_mask=interval_masker(tmp_path, "chr1\t5\t15\n")) == "A"


class TestConsensusFilterAndDepth:
    @pytest.mark.parametrize("flt", ["str10", "baq_dropout", "LowQual", "q10;str10"])
    def test_non_pass_filter_is_masked(self, flt):
        assert call_consensus([make_rec(flt=flt)]) == "X"

    @pytest.mark.parametrize("flt", ["PASS", ".", ""])
    def test_pass_like_filters_are_accepted(self, flt):
        assert call_consensus([make_rec(ref="A", flt=flt)]) == "A"

    def test_a_single_failing_record_masks_the_whole_position(self):
        recs = [make_rec(pos=5, alt="C", sample="1/0:30:15,15"),
                make_rec(pos=5, alt="T", flt="str10", sample="0/1:30:15,15")]
        assert call_consensus(recs) == "X"

    def test_filter_is_checked_before_depth(self):
        # DP=0 alone would give the no-call char; the failing FILTER wins.
        rec = make_rec(flt="str10", fmt="GT:DP", sample="0/0:0")
        assert call_consensus([rec]) == "X", "FILTER is evaluated before the depth gate"

    def test_depth_at_or_below_min_dp_is_a_no_call(self):
        rec = make_rec(fmt="GT:DP", sample="0/0:5")
        assert call_consensus([rec], min_dp=5) == "-", "the gate is `dp <= min_dp`, not `<`"

    def test_depth_just_above_min_dp_proceeds(self):
        rec = make_rec(ref="A", fmt="GT:DP:AD", sample="0/0:6:6,0")
        assert call_consensus([rec], min_dp=5, ref_min_dp=0) == "A"

    def test_zero_depth_is_a_no_call_with_default_min_dp(self):
        assert call_consensus([make_rec(fmt="GT:DP", sample="0/0:0")]) == "-"

    def test_nocall_char_is_configurable(self):
        rec = make_rec(fmt="GT:DP", sample="0/0:0")
        assert call_consensus([rec], nocall_char="?") == "?"

    def test_depth_is_the_maximum_across_records_at_the_position(self):
        recs = [make_rec(pos=5, alt="C", fmt="GT:DP:AD", sample="1/0:0:0,0"),
                make_rec(pos=5, alt="T", fmt="GT:DP:AD", sample="0/1:50:25,25")]
        assert call_consensus(recs, min_dp=10) == "Y", \
            "dp_pos is max() over the records, so one deep record rescues the site"


class TestConsensusReferenceBaseSource:
    def test_reference_fasta_wins_over_the_record_ref_column(self):
        rec = make_rec(pos=2, ref="A", alt=".", sample="0/0:30:30,0")
        assert call_consensus([rec], ref_seqs={"chr1": "GGGG"}) == "G"

    def test_record_ref_is_used_when_the_contig_is_absent(self):
        rec = make_rec(pos=2, ref="C", alt=".", sample="0/0:30:30,0")
        assert call_consensus([rec], ref_seqs={"other": "GGGG"}) == "C"

    def test_record_ref_is_used_when_the_position_is_past_the_contig_end(self):
        rec = make_rec(pos=99, ref="T", alt=".", sample="0/0:30:30,0")
        assert call_consensus([rec], ref_seqs={"chr1": "GG"}) == "T"

    def test_only_the_first_character_of_a_multi_base_ref_is_used(self):
        rec = make_rec(pos=1, ref="CTT", alt=".", sample="0/0:30:30,0")
        assert call_consensus([rec]) == "C"

    def test_reference_base_is_uppercased(self):
        rec = make_rec(pos=1, ref="a", alt=".", sample="0/0:30:30,0")
        assert call_consensus([rec], ref_seqs={"chr1": "acgt"}) == "A"

    def test_falls_back_to_n_when_no_reference_is_available(self):
        rec = make_rec(pos=1, ref="", alt=".", sample="0/0:100:100,0")
        assert call_consensus([rec]) == "N", \
            "with no FASTA and an empty REF the reference base defaults to N"

    def test_non_dna_reference_base_is_n_at_a_monomorphic_site(self):
        rec = make_rec(pos=1, ref="N", alt=".", sample="0/0:100:100,0")
        assert call_consensus([rec]) == "N"


class TestConsensusMonomorphicGates:
    """The confidence gates that guard a monomorphic-reference (ALT '.') call."""

    def test_confident_reference_call_emits_the_reference_base(self):
        rec = make_rec(ref="A", alt=".", sample="0/0:30:30,0")
        assert call_consensus([rec]) == "A"

    def test_empty_alt_column_is_treated_as_monomorphic(self):
        rec = make_rec(ref="A", alt="", sample="0/0:30:30,0")
        assert call_consensus([rec]) == "A"

    def test_missing_genotype_is_n(self):
        rec = make_rec(ref="A", alt=".", sample="./.:30:30,0")
        assert call_consensus([rec]) == "N"

    def test_non_homozygous_reference_genotype_is_n(self):
        rec = make_rec(ref="A", alt=".", sample="0/1:30:30,0")
        assert call_consensus([rec]) == "N", \
            "a non-ref allele index with no ALT to resolve it cannot be trusted"

    def test_haploid_reference_genotype_is_accepted(self):
        rec = make_rec(ref="A", alt=".", sample="0:30:30,0")
        assert call_consensus([rec]) == "A"

    def test_base_depth_below_ref_min_dp_is_n(self):
        rec = make_rec(ref="A", alt=".", sample="0/0:30:5,0")
        assert call_consensus([rec], max_ref_del_frac=1.0) == "N", \
            "RO+AO (5) is below ref_min_dp (10) even though DP (30) is not"

    def test_base_depth_uses_ro_plus_ao_not_dp(self):
        rec = make_rec(ref="A", alt=".", fmt="GT:DP:RO:AO", sample="0/0:100:8:1")
        assert call_consensus([rec], max_ref_del_frac=1.0) == "N", \
            "the gate must read base-calling depth (RO+AO=9), not raw DP=100"

    def test_base_depth_falls_back_to_dp_without_ro_ao_or_ad(self):
        rec = make_rec(ref="A", alt=".", fmt="GT:DP", sample="0/0:30")
        assert call_consensus([rec]) == "A", \
            "with no AD/RO/AO the raw DP stands in for base depth"

    def test_ref_min_dp_zero_disables_the_depth_gate(self):
        rec = make_rec(ref="A", alt=".", sample="0/0:3:3,0")
        assert call_consensus([rec], ref_min_dp=0) == "A"

    def test_deletion_fraction_at_or_above_the_threshold_is_n(self):
        # DP counts 100 reads but only 40 carry a base: 60% span a deletion.
        rec = make_rec(ref="A", alt=".", sample="0/0:100:40,0")
        assert call_consensus([rec]) == "N"

    def test_deletion_fraction_below_the_threshold_emits_the_reference(self):
        rec = make_rec(ref="A", alt=".", sample="0/0:100:60,0")
        assert call_consensus([rec]) == "A", "(100-60)/100 = 0.4 is below the 0.5 default"

    def test_deletion_gate_boundary_is_inclusive(self):
        rec = make_rec(ref="A", alt=".", sample="0/0:100:50,0")
        assert call_consensus([rec]) == "N", "the comparison is `>=`, so exactly 0.5 masks"

    def test_max_ref_del_frac_one_disables_the_deletion_gate(self):
        rec = make_rec(ref="A", alt=".", sample="0/0:100:10,0")
        assert call_consensus([rec], max_ref_del_frac=1.0) == "A"

    def test_alt_fraction_at_or_above_the_threshold_is_n(self):
        rec = make_rec(ref="A", alt=".", sample="0/0:100:90,10")
        assert call_consensus([rec]) == "N", "10/100 = 0.10 reaches max_ref_altfrac"

    def test_alt_fraction_below_the_threshold_emits_the_reference(self):
        rec = make_rec(ref="A", alt=".", sample="0/0:100:91,9")
        assert call_consensus([rec]) == "A"

    def test_alt_fraction_gate_needs_enough_alt_reads(self):
        # 1/10 = 0.10 clears max_ref_altfrac but only one alt read is below
        # max_ref_min_alt, so a single stray read cannot N a reference call.
        rec = make_rec(ref="A", alt=".", sample="0/0:10:9,1")
        assert call_consensus([rec]) == "A"

    def test_alt_fraction_gate_fires_at_exactly_max_ref_min_alt(self):
        rec = make_rec(ref="A", alt=".", sample="0/0:20:18,2")
        assert call_consensus([rec]) == "N"

    def test_alt_fraction_gate_works_from_ro_and_ao(self):
        rec = make_rec(ref="A", alt=".", fmt="GT:DP:RO:AO", sample="0/0:100:90:10")
        assert call_consensus([rec]) == "N"

    def test_max_ref_altfrac_one_disables_the_alt_fraction_gate(self):
        rec = make_rec(ref="A", alt=".", sample="0/0:100:50,50")
        assert call_consensus([rec], max_ref_altfrac=1.0, max_ref_del_frac=1.0) == "A"

    def test_alt_fraction_gate_is_skipped_without_ro_ao_evidence(self):
        # QUIRK: with no AD/RO/AO the alt-fraction gate cannot run at all, so a
        # DP-only record always passes it.
        rec = make_rec(ref="A", alt=".", fmt="GT:DP", sample="0/0:100")
        assert call_consensus([rec]) == "A"


class TestConsensusVariantSites:
    def test_homozygous_alt_snp(self):
        rec = make_rec(ref="A", alt="T", sample="1/1:30:0,30")
        assert call_consensus([rec]) == "T"

    def test_heterozygous_snp_becomes_an_iupac_code(self):
        rec = make_rec(ref="A", alt="G", sample="0/1:30:15,15")
        assert call_consensus([rec]) == "R"

    def test_three_way_multiallelic_in_one_record(self):
        rec = make_rec(ref="A", alt="C,G", sample="1/2:30:0,15,15")
        assert call_consensus([rec], ref_seqs={"chr1": "A"}) == "S", \
            "{C,G} without the reference allele is S"

    def test_indel_allele_is_masked(self):
        rec = make_rec(ref="A", alt="AGG", sample="1/1:30:0,30")
        assert call_consensus([rec]) == "X"

    def test_out_of_range_allele_index_is_masked(self):
        rec = make_rec(ref="A", alt="T", sample="3/3:30:0,30")
        assert call_consensus([rec]) == "X"

    def test_missing_genotype_at_a_variant_site_is_masked(self):
        rec = make_rec(ref="A", alt="T", sample="./.:30:15,15")
        assert call_consensus([rec]) == "X"

    def test_variant_sites_skip_the_monomorphic_confidence_gates(self):
        # A shallow ALT call is emitted as-is; ref_min_dp only guards ALT '.' sites.
        rec = make_rec(ref="A", alt="T", fmt="GT:DP:AD", sample="1/1:2:0,2")
        assert call_consensus([rec], ref_min_dp=100) == "T"


class TestConsensusMultiRecordUnion:
    def test_norm_split_multiallelic_unions_both_alts(self):
        recs = [make_rec(pos=5, ref="A", alt="C", sample="1/0:30:15,15"),
                make_rec(pos=5, ref="A", alt="G", sample="0/1:30:15,15")]
        assert call_consensus(recs, ref_seqs={"chr1": "AAAAA"}) == "S", \
            "the per-record 0 is a split artifact, so the reference base is dropped"

    def test_single_record_keeps_the_reference_allele(self):
        # Same genotype shape, one record: here 0 is a real reference call.
        rec = make_rec(pos=5, ref="A", alt="C", sample="0/1:30:15,15")
        assert call_consensus([rec], ref_seqs={"chr1": "AAAAA"}) == "M"

    def test_three_records_union_to_a_three_base_code(self):
        recs = [make_rec(pos=5, ref="T", alt="A", sample="1/0:30:10,10"),
                make_rec(pos=5, ref="T", alt="C", sample="0/1:30:10,10"),
                make_rec(pos=5, ref="T", alt="G", sample="0/1:30:10,10")]
        assert call_consensus(recs, ref_seqs={"chr1": "TTTTT"}) == "V"

    def test_union_of_all_four_bases_is_n(self):
        recs = [make_rec(pos=5, ref="T", alt="A", sample="1/0:30:10,10"),
                make_rec(pos=5, ref="T", alt="C", sample="0/1:30:10,10"),
                make_rec(pos=5, ref="T", alt="G,T", sample="1/2:30:5,5,5")]
        assert call_consensus(recs, ref_seqs={"chr1": "TTTTT"}) == "N"

    def test_one_non_snp_sibling_masks_the_position(self):
        recs = [make_rec(pos=5, ref="A", alt="C", sample="1/0:30:15,15"),
                make_rec(pos=5, ref="A", alt="ACC", sample="0/1:30:15,15")]
        assert call_consensus(recs, ref_seqs={"chr1": "AAAAA"}) == "X"

    def test_a_monomorphic_sibling_masks_the_position(self):
        # QUIRK: with ignore_ref the ALT '.' record contributes nothing, which
        # allele_bases_from_record reports as None, i.e. the mask signal. A VCF
        # that carries both an ALT '.' row and a variant row at one position
        # therefore masks rather than reporting the variant.
        recs = [make_rec(pos=5, ref="A", alt=".", sample="0/0:30:30,0"),
                make_rec(pos=5, ref="A", alt="C", sample="0/1:30:15,15")]
        assert call_consensus(recs, ref_seqs={"chr1": "AAAAA"}) == "X"

    def test_a_star_only_sibling_masks_the_position(self):
        # QUIRK: the '*' allele is skipped and, with ignore_ref, the sibling's 0 is
        # dropped too, so that record contributes an empty set. allele_bases_from_record
        # reports empty as None, i.e. mask, so the site is X rather than C.
        recs = [make_rec(pos=5, ref="A", alt="*", sample="1/0:30:15,15"),
                make_rec(pos=5, ref="A", alt="C", sample="0/1:30:15,15")]
        assert call_consensus(recs, ref_seqs={"chr1": "AAAAA"}) == "X"


# --------------------------------------------------------------------------- #
# IntervalMasker
# --------------------------------------------------------------------------- #

class TestIntervalMasker:
    def test_no_file_masks_nothing(self):
        assert wgs.IntervalMasker(None).masked("chr1", 10) is False

    def test_missing_file_warns_without_raising(self, tmp_path, capsys):
        masker = wgs.IntervalMasker(str(tmp_path / "nope.tsv"))
        assert masker.masked("chr1", 10) is False
        assert "Exclude file not found" in capsys.readouterr().err

    @pytest.mark.parametrize("pos,expected", [
        (4, False), (5, True), (10, True), (15, True), (16, False),
    ])
    def test_interval_bounds_are_inclusive(self, tmp_path, pos, expected):
        masker = interval_masker(tmp_path, "chr1\t5\t15\n")
        assert masker.masked("chr1", pos) is expected

    def test_overlapping_intervals_are_merged(self, tmp_path):
        masker = interval_masker(tmp_path, "chr1\t1\t10\nchr1\t5\t20\n")
        assert masker.intervals["chr1"] == [(1, 20)]

    def test_adjacent_intervals_are_merged(self, tmp_path):
        masker = interval_masker(tmp_path, "chr1\t1\t5\nchr1\t6\t10\n")
        assert masker.intervals["chr1"] == [(1, 10)], \
            "intervals touching end+1 are joined into one"

    def test_intervals_with_a_one_base_gap_stay_separate(self, tmp_path):
        masker = interval_masker(tmp_path, "chr1\t1\t5\nchr1\t7\t10\n")
        assert masker.intervals["chr1"] == [(1, 5), (7, 10)]

    def test_reversed_interval_is_swapped(self, tmp_path):
        masker = interval_masker(tmp_path, "chr1\t20\t10\n")
        assert masker.intervals["chr1"] == [(10, 20)]
        assert masker.masked("chr1", 15) is True

    def test_contained_interval_is_absorbed(self, tmp_path):
        masker = interval_masker(tmp_path, "chr1\t1\t100\nchr1\t20\t30\n")
        assert masker.intervals["chr1"] == [(1, 100)]

    def test_comment_header_and_short_lines_are_skipped(self, tmp_path):
        body = (
            "# a comment\n"
            "chrom\tstart\tend\n"
            "\n"
            "chr1\t5\n"
            "chr1\tstart\tend\n"
            "chr1\t5\t15\n"
        )
        masker = interval_masker(tmp_path, body)
        assert masker.intervals == {"chr1": [(5, 15)]}

    def test_contig_named_like_the_header_is_dropped(self, tmp_path):
        # BUG (pinned, not fixed): the header check is a case-insensitive
        # startswith("chrom"), so a real contig called e.g. "chromosome_1" is
        # silently discarded and never masked.
        masker = interval_masker(tmp_path, "chromosome_1\t5\t15\n")
        assert masker.intervals == {}
        assert masker.masked("chromosome_1", 10) is False

    def test_intervals_are_tracked_per_contig(self, tmp_path):
        masker = interval_masker(tmp_path, "chr1\t5\t15\nchr2\t100\t200\n")
        assert masker.masked("chr1", 10) is True
        assert masker.masked("chr2", 150) is True
        assert masker.masked("chr2", 10) is False

    def test_ascending_queries_across_several_intervals(self, tmp_path):
        masker = interval_masker(tmp_path, "chr1\t5\t10\nchr1\t20\t25\n")
        assert [masker.masked("chr1", p) for p in (1, 7, 12, 22, 30)] == \
            [False, True, False, True, False]

    def test_forward_only_pointer_misses_a_descending_query(self, tmp_path):
        # DOCUMENTED LIMITATION (pinned): masked() advances a per-contig pointer
        # and never rewinds, so it is only correct for ascending positions. That
        # holds for build_consensus, which walks the VCF in order. This test pins
        # the limitation rather than asserting the "right" answer.
        masker = interval_masker(tmp_path, "chr1\t10\t20\n")
        assert masker.masked("chr1", 15) is True
        assert masker.masked("chr1", 25) is False
        assert masker.masked("chr1", 15) is False, \
            "after querying past the interval the pointer cannot go back"


# --------------------------------------------------------------------------- #
# PointMasker
# --------------------------------------------------------------------------- #

class TestPointMasker:
    def test_no_file_masks_nothing(self):
        assert wgs.PointMasker(None).masked("chr1", 10) is False

    def test_missing_file_warns_without_raising(self, tmp_path, capsys):
        masker = wgs.PointMasker(str(tmp_path / "nope.tsv"))
        assert masker.masked("chr1", 10) is False
        assert "Mask-sites file not found" in capsys.readouterr().err

    def test_masks_listed_positions_only(self, tmp_path):
        masker = point_masker(tmp_path, "chr1\t10\nchr1\t12\nchr2\t10\n")
        assert masker.masked("chr1", 10) is True
        assert masker.masked("chr1", 11) is False
        assert masker.masked("chr2", 10) is True
        assert masker.masked("chr3", 10) is False

    def test_lookup_is_stateless_and_order_independent(self, tmp_path):
        masker = point_masker(tmp_path, "chr1\t10\nchr1\t20\n")
        assert masker.masked("chr1", 20) is True
        assert masker.masked("chr1", 10) is True, \
            "unlike IntervalMasker this is a plain set lookup, so order is irrelevant"

    def test_comments_blanks_short_and_non_numeric_lines_are_skipped(self, tmp_path):
        body = "# comment\n\nchr1\nchr1\tNA\nchr1\t10\textra\n"
        masker = point_masker(tmp_path, body)
        assert masker.points == {"chr1": {10}}, "extra columns past POS are ignored"


# --------------------------------------------------------------------------- #
# CsvMasker
# --------------------------------------------------------------------------- #

CSV_HEADER = "rv_position,pal,homopolymer,GCrich,repetitive,blindspot\n"


class TestCsvMasker:
    def test_no_file_masks_nothing(self):
        assert wgs.CsvMasker(None).masked("chr1", 10) is False

    def test_missing_file_warns_without_raising(self, tmp_path, capsys):
        masker = wgs.CsvMasker(str(tmp_path / "nope.csv"), ["chr1"])
        assert masker.masked("chr1", 10) is False
        assert "Special Mask CSV file not found" in capsys.readouterr().err

    def test_masks_on_repetitive_flag(self, tmp_path):
        masker = csv_masker(tmp_path, CSV_HEADER + "10,0,0,0,1,0\n", ["chr1"])
        assert masker.masked("chr1", 10) is True

    def test_masks_on_blindspot_flag(self, tmp_path):
        masker = csv_masker(tmp_path, CSV_HEADER + "10,0,0,0,0,1\n", ["chr1"])
        assert masker.masked("chr1", 10) is True

    def test_other_flags_do_not_mask(self, tmp_path):
        masker = csv_masker(tmp_path, CSV_HEADER + "10,1,1,1,0,0\n", ["chr1"])
        assert masker.masked("chr1", 10) is False, \
            "pal / homopolymer / GCrich are read but never mask"

    @pytest.mark.parametrize("flag,masked", [
        ("1", True), ("1.0", True), ("0", False), ("", False),
        ("2", False), ("yes", False),
    ])
    def test_flag_parsing(self, tmp_path, flag, masked):
        masker = csv_masker(tmp_path, CSV_HEADER + f"10,0,0,0,{flag},0\n", ["chr1"])
        assert masker.masked("chr1", 10) is masked, \
            "only a flag that parses to exactly 1 masks"

    def test_rows_with_a_bad_position_are_skipped_individually(self, tmp_path):
        body = CSV_HEADER + "NA,0,0,0,1,0\n11,0,0,0,1,0\n"
        masker = csv_masker(tmp_path, body, ["chr1"])
        assert masker.points == {11}, "a malformed row must not drop the rest of the file"

    def test_missing_rv_position_header_warns_and_masks_nothing(self, tmp_path, capsys):
        body = "position,repetitive,blindspot\n10,1,0\n"
        masker = csv_masker(tmp_path, body, ["chr1"])
        assert masker.points == set()
        assert "missing 'rv_position' header" in capsys.readouterr().err

    def test_mask_is_bound_to_the_single_reference_contig(self, tmp_path):
        masker = csv_masker(tmp_path, CSV_HEADER + "10,0,0,0,1,0\n", ["chr1"])
        assert masker.contig == "chr1"
        assert masker.masked("chr2", 10) is False, \
            "single-contig rv_position coordinates must not leak onto another contig"

    def test_without_contigs_the_mask_applies_to_any_contig(self, tmp_path):
        masker = csv_masker(tmp_path, CSV_HEADER + "10,0,0,0,1,0\n", None)
        assert masker.contig is None
        assert masker.masked("anything", 10) is True

    def test_multiple_contigs_raise(self, tmp_path):
        with pytest.raises(RuntimeError, match="single-contig"):
            csv_masker(tmp_path, CSV_HEADER + "10,0,0,0,1,0\n", ["chr1", "chr2"])

    def test_multiple_contigs_without_a_csv_do_not_raise(self):
        assert wgs.CsvMasker(None, ["chr1", "chr2"]).masked("chr1", 10) is False


# --------------------------------------------------------------------------- #
# FastaReader
# --------------------------------------------------------------------------- #

class TestFastaReader:
    def test_reads_multiple_contigs(self, tmp_path):
        path = write_fasta(tmp_path / "ref.fa", {"chr1": "ACGT", "chr2": "TTTT"}, wrap=2)
        assert wgs.FastaReader(path).seqs == {"chr1": "ACGT", "chr2": "TTTT"}

    def test_sequence_is_uppercased_and_unwrapped(self, tmp_path):
        path = write_text(tmp_path / "ref.fa", ">chr1\nac\n  gt  \nAC\n")
        assert wgs.FastaReader(path).seqs == {"chr1": "ACGTAC"}, \
            "line whitespace is stripped and the sequence is joined and uppercased"

    def test_header_is_truncated_at_the_first_whitespace(self, tmp_path):
        path = write_text(tmp_path / "ref.fa", ">chr1 Mycobacterium tuberculosis\nACGT\n")
        assert list(wgs.FastaReader(path).seqs) == ["chr1"]

    def test_reads_a_gzipped_reference(self, tmp_path):
        path = write_text(tmp_path / "ref.fa.gz", ">chr1\nACGT\n", gz=True)
        assert wgs.FastaReader(path).seqs == {"chr1": "ACGT"}

    def test_empty_file_yields_no_sequences(self, tmp_path):
        assert wgs.FastaReader(write_text(tmp_path / "ref.fa", "")).seqs == {}

    def test_bare_greater_than_header_raises_index_error(self, tmp_path):
        # BUG (pinned, not fixed): `line[1:].split()[0]` has no guard, so a header
        # line that is just ">" fails with an opaque IndexError instead of a clear
        # "malformed FASTA" message.
        path = write_text(tmp_path / "ref.fa", ">\nACGT\n")
        with pytest.raises(IndexError):
            wgs.FastaReader(path)


# --------------------------------------------------------------------------- #
# FastaWriter
# --------------------------------------------------------------------------- #

class TestFastaWriter:
    def test_writes_header_and_bases(self, tmp_path):
        out = tmp_path / "o.fa"
        w = wgs.FastaWriter(str(out), wrap=0)
        w.start("chr1")
        for b in "ACGT":
            w.write_base(b)
        w.close()
        assert out.read_text(encoding="utf-8") == ">chr1\nACGT\n"

    def test_wraps_at_the_requested_width(self, tmp_path):
        out = tmp_path / "o.fa"
        w = wgs.FastaWriter(str(out), wrap=3)
        w.start("chr1")
        for b in "ABCDEFG":
            w.write_base(b)
        w.close()
        assert out.read_text(encoding="utf-8") == ">chr1\nABC\nDEF\nG\n"

    @pytest.mark.parametrize("wrap", [0, -5])
    def test_non_positive_wrap_disables_wrapping(self, tmp_path, wrap):
        out = tmp_path / "o.fa"
        w = wgs.FastaWriter(str(out), wrap=wrap)
        w.start("chr1")
        for b in "ABCDEFG":
            w.write_base(b)
        w.close()
        assert out.read_text(encoding="utf-8") == ">chr1\nABCDEFG\n"

    def test_empty_base_is_written_as_n(self, tmp_path):
        out = tmp_path / "o.fa"
        w = wgs.FastaWriter(str(out), wrap=0)
        w.start("chr1")
        w.write_base("")
        w.close()
        assert out.read_text(encoding="utf-8") == ">chr1\nN\n"

    def test_only_the_first_character_is_written_and_uppercased(self, tmp_path):
        out = tmp_path / "o.fa"
        w = wgs.FastaWriter(str(out), wrap=0)
        w.start("chr1")
        w.write_base("acgt")
        w.close()
        assert out.read_text(encoding="utf-8") == ">chr1\nA\n"

    def test_counts_bases_per_contig(self, tmp_path):
        w = wgs.FastaWriter(str(tmp_path / "o.fa"), wrap=0)
        w.start("chr1")
        for _ in range(3):
            w.write_base("A")
        w.start("chr2")
        w.write_base("C")
        w.close()
        assert w.counts == {"chr1": 3, "chr2": 1}

    def test_a_new_contig_closes_the_previous_line(self, tmp_path):
        out = tmp_path / "o.fa"
        w = wgs.FastaWriter(str(out), wrap=0)
        w.start("chr1")
        w.write_base("A")
        w.start("chr2")
        w.write_base("C")
        w.close()
        assert out.read_text(encoding="utf-8") == ">chr1\nA\n>chr2\nC\n"

    def test_no_blank_line_when_the_wrap_already_closed_the_line(self, tmp_path):
        out = tmp_path / "o.fa"
        w = wgs.FastaWriter(str(out), wrap=2)
        w.start("chr1")
        w.write_base("A")
        w.write_base("C")
        w.start("chr2")
        w.write_base("G")
        w.close()
        assert out.read_text(encoding="utf-8") == ">chr1\nAC\n>chr2\nG\n"


# --------------------------------------------------------------------------- #
# build_consensus: end to end
# --------------------------------------------------------------------------- #

# A 13 bp reference exercised by every branch of the writer in one pass.
E2E_REF = "ACGTACGTACGTA"

E2E_VCF = [
    # pos 1: confident monomorphic reference call -> A
    vcf_line("chr1", 1, "A", ".", sample="0/0:30:30,0"),
    # pos 2: absent from the VCF and unmasked -> gap filled with the no-call char
    # pos 3: absent from the VCF but inside an exclude interval -> gap filled with X
    # pos 4: base depth 5 < ref_min_dp 10 -> N
    vcf_line("chr1", 4, "T", ".", sample="0/0:5:5,0"),
    # pos 5: non-PASS FILTER -> X
    vcf_line("chr1", 5, "A", "G", flt="str10", sample="1/1:30:0,30"),
    # pos 6: DP 0 <= min_dp 0 -> no-call char
    vcf_line("chr1", 6, "C", ".", fmt="GT:DP", sample="0/0:0"),
    # pos 7: heterozygous G/A -> R (and present in the CSV with both flags 0)
    vcf_line("chr1", 7, "G", "A", sample="0/1:30:15,15"),
    # pos 8: a confident reference call that --mask-sites overrides -> X
    vcf_line("chr1", 8, "T", ".", sample="0/0:40:40,0"),
    # pos 9: homozygous alt -> C
    vcf_line("chr1", 9, "A", "C", sample="1/1:40:0,40"),
    # pos 10: norm-split multiallelic A + T with the artifact ref allele dropped -> W
    vcf_line("chr1", 10, "C", "A", sample="1/0:40:20,20"),
    vcf_line("chr1", 10, "C", "T", sample="0/1:40:20,20"),
    # pos 11: monomorphic but 20% alt reads -> N
    vcf_line("chr1", 11, "G", ".", sample="0/0:100:80,20"),
    # pos 12: confident monomorphic reference call -> T
    vcf_line("chr1", 12, "T", ".", sample="0/0:30:30,0"),
    # pos 13: a confident reference call that the special CSV mask overrides -> X
    vcf_line("chr1", 13, "A", ".", sample="0/0:30:30,0"),
]

E2E_EXPECTED = "A-XNX-RXCWNTX"


class TestBuildConsensusEndToEnd:
    @pytest.fixture()
    def e2e_out(self, tmp_path):
        exclude = write_text(
            tmp_path / "exclude.tsv",
            "# regions to drop\nchrom\tstart\tend\nchr1\t3\t3\n",
        )
        mask_sites = write_text(tmp_path / "mask_sites.tsv", "# sites\nchr1\t8\n")
        special = write_text(
            tmp_path / "special.csv",
            CSV_HEADER + "13,0,0,0,0,1\n7,0,0,0,0,0\n",
        )
        return run_build(
            tmp_path, E2E_VCF, {"chr1": E2E_REF},
            exclude_path=exclude, mask_sites_path=mask_sites,
            special_csv_path=special, wrap=5,
        )

    def test_consensus_sequence(self, e2e_out):
        assert read_fasta(e2e_out) == {"chr1": E2E_EXPECTED}

    def test_consensus_length_matches_the_reference(self, e2e_out):
        assert len(read_fasta(e2e_out)["chr1"]) == len(E2E_REF)

    def test_output_is_wrapped(self, e2e_out):
        lines = e2e_out.read_text(encoding="utf-8").splitlines()
        assert lines[0] == ">chr1"
        assert lines[1:] == ["A-XNX", "-RXCW", "NTX"]

    def test_two_contigs_are_written_in_vcf_order(self, tmp_path):
        vcf = [vcf_line("chr1", 1, "A", "."), vcf_line("chr1", 2, "C", "."),
               vcf_line("chr2", 1, "G", "."), vcf_line("chr2", 2, "T", ".")]
        out = run_build(tmp_path, vcf, {"chr1": "AC", "chr2": "GT"})
        assert read_fasta(out) == {"chr1": "AC", "chr2": "GT"}

    def test_reads_gzipped_inputs(self, tmp_path):
        vcf = write_vcf(tmp_path / "in.vcf.gz",
                        [vcf_line("chr1", 1, "A", "."), vcf_line("chr1", 2, "C", ".")],
                        gz=True)
        ref = write_fasta(tmp_path / "ref.fa.gz", {"chr1": "AC"}, gz=True)
        out = tmp_path / "consensus.fa"
        wgs.build_consensus(vcf_path=vcf, reference_path=ref, out_path=str(out),
                            **BUILD_DEFAULTS)
        assert read_fasta(out) == {"chr1": "AC"}

    def test_leading_gap_is_filled_from_position_one(self, tmp_path):
        out = run_build(tmp_path, [vcf_line("chr1", 3, "G", ".")], {"chr1": "ACG"})
        assert read_fasta(out) == {"chr1": "--G"}

    def test_malformed_and_header_lines_are_skipped(self, tmp_path):
        vcf = ["##contig=<ID=chr1>", "not\ta\tvcf\tline",
               vcf_line("chr1", 1, "A", "."), vcf_line("chr1", 2, "C", ".")]
        out = run_build(tmp_path, vcf, {"chr1": "AC"})
        assert read_fasta(out) == {"chr1": "AC"}

    def test_custom_mask_and_nocall_characters(self, tmp_path):
        mask_sites = write_text(tmp_path / "mask_sites.tsv", "chr1\t1\n")
        vcf = [vcf_line("chr1", 1, "A", "."),
               vcf_line("chr1", 3, "G", ".", fmt="GT:DP", sample="0/0:0")]
        out = run_build(tmp_path, vcf, {"chr1": "ACG"},
                        mask_sites_path=mask_sites, mask_char="?", nocall_char="~")
        assert read_fasta(out) == {"chr1": "?~~"}, \
            "position 2 is a filled gap and position 3 is a zero-depth record"


class TestBuildConsensusErrors:
    def test_contig_missing_from_the_reference_raises_under_strict(self, tmp_path):
        vcf = [vcf_line("chr1", 1, "A", "."), vcf_line("chrZ", 1, "G", ".")]
        with pytest.raises(RuntimeError, match="Contig 'chrZ' not found in reference"):
            run_build(tmp_path, vcf, {"chr1": "A"}, strict=True)

    def test_contig_missing_from_the_reference_is_skipped_without_strict(self, tmp_path):
        vcf = [vcf_line("chr1", 1, "A", "."), vcf_line("chrZ", 1, "G", ".")]
        out = run_build(tmp_path, vcf, {"chr1": "A"}, strict=False)
        assert read_fasta(out) == {"chr1": "A"}

    def test_short_output_raises_a_length_mismatch(self, tmp_path):
        vcf = [vcf_line("chr1", 1, "A", "."), vcf_line("chr1", 2, "C", ".")]
        with pytest.raises(RuntimeError, match=r"consensus length mismatch for chr1: wrote 2 != reference 5"):
            run_build(tmp_path, vcf, {"chr1": "ACGTA"})

    def test_a_reference_contig_with_no_records_raises(self, tmp_path):
        vcf = [vcf_line("chr1", 1, "A", "."), vcf_line("chr1", 2, "C", ".")]
        with pytest.raises(RuntimeError, match=r"missing contig\(s\) entirely: \['chr2'\]"):
            run_build(tmp_path, vcf, {"chr1": "AC", "chr2": "GT"})

    def test_empty_vcf_raises(self, tmp_path):
        with pytest.raises(RuntimeError, match=r"missing contig\(s\) entirely"):
            run_build(tmp_path, [], {"chr1": "AC"})

    def test_out_of_order_position_raises_under_strict(self, tmp_path):
        vcf = [vcf_line("chr1", 1, "A", "."), vcf_line("chr1", 3, "G", "."),
               vcf_line("chr1", 2, "C", ".")]
        with pytest.raises(RuntimeError, match=r"Unexpected duplicate/out-of-order chr1:2"):
            run_build(tmp_path, vcf, {"chr1": "ACG"}, strict=True)

    def test_out_of_order_position_is_dropped_without_strict(self, tmp_path):
        # The late record is discarded and position 2 keeps the gap-filled no-call.
        vcf = [vcf_line("chr1", 1, "A", "."), vcf_line("chr1", 3, "G", "."),
               vcf_line("chr1", 2, "C", ".")]
        out = run_build(tmp_path, vcf, {"chr1": "ACG"}, strict=False)
        assert read_fasta(out) == {"chr1": "A-G"}

    @pytest.mark.parametrize("strict", [True, False])
    def test_a_contig_that_reappears_raises_regardless_of_strict(self, tmp_path, strict):
        vcf = [vcf_line("chr1", 1, "A", "."), vcf_line("chr2", 1, "G", "."),
               vcf_line("chr1", 2, "C", ".")]
        with pytest.raises(RuntimeError, match=r"Contig 'chr1' reappears out of order"):
            run_build(tmp_path, vcf, {"chr1": "AC", "chr2": "G"}, strict=strict)
