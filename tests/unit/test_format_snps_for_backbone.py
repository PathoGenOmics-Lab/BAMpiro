"""bin/format_snps_for_backbone.awk - the last rewrite a SNP gets before the consensus.

It decides the depth the consensus divides by and, more importantly, whether a call stays
heterozygous. Both FreeBayes paths pipe through it, and until it was pulled out of the process
script it had no coverage: a stub run skips the script block entirely.

The interesting rule is the het re-validation. `bcftools norm -m -` splits a multiallelic 1/2 into
biallelic records whose reference index 0 is an artefact with RO near zero. Left alone, that reads
as a heterozygote and the reference base leaks into the consensus IUPAC code.
"""

from __future__ import annotations

import subprocess

import pytest

from conftest import BIN

AWK = BIN / "format_snps_for_backbone.awk"

HEADER = (
    "##fileformat=VCFv4.2\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
)


def run(rows, hetmin="0.10", tmp_path=None):
    """Run the awk program over a small VCF and return the data rows, split into fields."""
    vcf = tmp_path / "in.vcf"
    vcf.write_text(HEADER + "".join(rows))
    out = subprocess.run(["awk", "-v", f"HETMIN={hetmin}", "-f", str(AWK), str(vcf)],
                         capture_output=True, text=True, check=True)
    return [line.split("\t") for line in out.stdout.splitlines() if not line.startswith("#")]


def row(info, gt, dp_field=None, pos=100):
    fmt_val = gt if dp_field is None else f"{gt}:{dp_field}"
    fmt = "GT" if dp_field is None else "GT:DP"
    return f"c\t{pos}\t.\tA\tT\t50\t.\t{info}\t{fmt}\t{fmt_val}\n"


def info_of(fields):
    return dict(kv.split("=", 1) for kv in fields[7].split(";"))


class TestHeaders:
    def test_the_header_passes_through_unchanged(self, tmp_path):
        vcf = tmp_path / "in.vcf"
        vcf.write_text(HEADER + row("DP=10;RO=1;AO=9", "1/1", "10"))
        out = subprocess.run(["awk", "-v", "HETMIN=0.10", "-f", str(AWK), str(vcf)],
                             capture_output=True, text=True, check=True).stdout
        assert out.startswith("##fileformat=VCFv4.2\n#CHROM\t")

    def test_a_header_only_vcf_produces_no_data_rows(self, tmp_path):
        assert run([], tmp_path=tmp_path) == []


class TestOutputShape:
    def test_info_is_replaced_with_the_backbone_tags(self, tmp_path):
        [f] = run([row("DP=30;RO=2;AO=28", "1/1", "30")], tmp_path=tmp_path)
        assert f[7] == "ADP=30;WT=0;HET=0;HOM=1;NC=0"

    def test_format_is_reduced_to_genotype_and_depth(self, tmp_path):
        [f] = run([row("DP=30;RO=2;AO=28", "1/1", "30")], tmp_path=tmp_path)
        assert f[8] == "GT:DP"
        assert f[9] == "1/1:30"

    def test_the_leading_columns_are_untouched(self, tmp_path):
        [f] = run([row("DP=30;RO=2;AO=28", "1/1", "30", pos=4411532)], tmp_path=tmp_path)
        assert f[:7] == ["c", "4411532", ".", "A", "T", "50", "."]


class TestDepth:
    def test_depth_comes_from_info_dp_first(self, tmp_path):
        [f] = run([row("DP=42;RO=2;AO=28", "1/1", "99")], tmp_path=tmp_path)
        assert info_of(f)["ADP"] == "42"

    def test_a_lookalike_info_key_is_not_mistaken_for_dp(self, tmp_path):
        """Keys are anchored, so MQMDP= and ADP= must not be read as the depth."""
        [f] = run([row("MQMDP=999;ADP=888;RO=4;AO=6", "1/1", "12")], tmp_path=tmp_path)
        assert info_of(f)["ADP"] == "12", "depth was taken from a look-alike INFO key"

    def test_depth_falls_back_to_the_format_column(self, tmp_path):
        [f] = run([row("RO=4;AO=6", "1/1", "17")], tmp_path=tmp_path)
        assert info_of(f)["ADP"] == "17"

    def test_depth_falls_back_to_ro_plus_ao(self, tmp_path):
        [f] = run([row("RO=4;AO=6", "1/1")], tmp_path=tmp_path)
        assert info_of(f)["ADP"] == "10"

    def test_depth_never_ends_up_zero(self, tmp_path):
        """The consensus divides by this, so the last fallback is 1 rather than 0."""
        [f] = run([row("RO=0;AO=0", "1/1")], tmp_path=tmp_path)
        assert info_of(f)["ADP"] == "1"

    def test_only_the_first_alt_observation_count_is_used(self, tmp_path):
        """AO is per-alt; after `bcftools norm -m -` the record is biallelic."""
        [f] = run([row("RO=4;AO=6,3", "1/1")], tmp_path=tmp_path)
        assert info_of(f)["ADP"] == "10"

    def test_a_missing_format_depth_lands_in_adp_verbatim(self, tmp_path):
        """BUG PIN: a '.' in FORMAT/DP is copied straight into ADP.

        `dat[i]` is a string, so the following `dp == 0` is a string comparison that '.' fails,
        and all three numeric fallbacks are skipped. Downstream `get_dp` reads the '.' as 0 and
        the position becomes a no-call, so a called SNP silently turns into a gap. Probably
        unreachable in production because the upstream filter requires FMT/DP >= filter_min_dp.
        """
        [f] = run([row("RO=3;AO=7", "0/1", ".")], tmp_path=tmp_path)
        assert info_of(f)["ADP"] == "."


class TestGenotype:
    def test_a_balanced_heterozygote_stays_heterozygous(self, tmp_path):
        [f] = run([row("DP=40;RO=20;AO=20", "0/1", "40")], tmp_path=tmp_path)
        assert f[9].startswith("0/1")
        assert info_of(f)["HET"] == "1" and info_of(f)["HOM"] == "0"

    def test_a_heterozygote_with_no_reference_support_becomes_homozygous(self, tmp_path):
        """The multiallelic-split artefact: reference index 0 exists but nothing supports it."""
        [f] = run([row("DP=30;RO=2;AO=28", "0/1", "30")], tmp_path=tmp_path)
        assert f[9].startswith("1/1"), "an unsupported reference allele must not read as a het"
        assert info_of(f)["HET"] == "0" and info_of(f)["HOM"] == "1"

    @pytest.mark.parametrize("gt", ["0/1", "1/0"])
    def test_both_spellings_of_a_heterozygote_are_re_validated(self, gt, tmp_path):
        [f] = run([row("DP=30;RO=1;AO=29", gt, "30")], tmp_path=tmp_path)
        assert f[9].startswith("1/1")

    def test_the_threshold_is_hetmin_and_it_is_exclusive(self, tmp_path):
        """RO fraction exactly at HETMIN stays a het; below it is rewritten."""
        at = run([row("DP=100;RO=10;AO=90", "0/1", "100")], tmp_path=tmp_path)[0]
        below = run([row("DP=100;RO=9;AO=91", "0/1", "100")], tmp_path=tmp_path)[0]
        assert at[9].startswith("0/1"), "0.10 is not below HETMIN, so it stays heterozygous"
        assert below[9].startswith("1/1")

    def test_the_threshold_follows_the_hetmin_argument(self, tmp_path):
        [f] = run([row("DP=100;RO=30;AO=70", "0/1", "100")], hetmin="0.40", tmp_path=tmp_path)
        assert f[9].startswith("1/1")

    def test_an_already_homozygous_call_is_left_alone(self, tmp_path):
        [f] = run([row("DP=30;RO=15;AO=15", "1/1", "30")], tmp_path=tmp_path)
        assert f[9].startswith("1/1")
        assert info_of(f)["HET"] == "0" and info_of(f)["HOM"] == "1"

    def test_a_genotype_that_is_neither_counts_as_neither(self, tmp_path):
        """A '.' or a 1/2 that survived is flagged as neither het nor hom rather than guessed at."""
        [f] = run([row("DP=30;RO=10;AO=20", "1/2", "30")], tmp_path=tmp_path)
        assert info_of(f)["HET"] == "0" and info_of(f)["HOM"] == "0"

    def test_a_heterozygote_with_no_observations_is_not_rewritten(self, tmp_path):
        """RO+AO of zero leaves the fraction undefined, so the call is left as it was."""
        [f] = run([row("DP=30;RO=0;AO=0", "0/1", "30")], tmp_path=tmp_path)
        assert f[9].startswith("0/1")


def test_the_two_freebayes_processes_share_this_one_program():
    """Both paths must pipe through the same file, or the het rule can drift between them."""
    variants = (BIN.parent / "modules" / "variants.nf").read_text()
    assert variants.count("format_snps_for_backbone.awk") == 2
    assert "gtype=\"1/1\"" not in variants, "the genotype logic is back inline in the module"
