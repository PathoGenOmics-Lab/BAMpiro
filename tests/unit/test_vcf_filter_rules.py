"""The bcftools expressions that define a homozygous or heterozygous SNP call.

These are the scientific definition of a call, not plumbing. They lived inline in two FreeBayes
processes, written out twice, and had no coverage at all: a stub run skips the script block and
nothing else could reach them. The expressions are also checked against a real bcftools here, so a
change that produces a syntactically valid but wrong filter still fails.
"""

from __future__ import annotations

import shutil
import subprocess
import sys
import textwrap

import pytest

from conftest import BIN, load_script

rules = load_script("vcf_filter_rules")

# Production defaults, from nextflow.config.
DEFAULTS = dict(min_dp=30, min_alt_fwd=2, min_alt_rev=2, hom_threshold=0.90, het_min_frac=0.10)


def hom(**kw):
    a = {**DEFAULTS, **kw}
    return rules.hom_rule(a["min_dp"], a["min_alt_fwd"], a["min_alt_rev"], a["hom_threshold"])


def het(**kw):
    a = {**DEFAULTS, **kw}
    return rules.het_rule(a["min_dp"], a["min_alt_fwd"], a["min_alt_rev"], a["hom_threshold"],
                          a["het_min_frac"])


class TestExpressionShape:
    def test_both_rules_require_depth_from_either_source(self):
        # A caller may record depth per-sample or per-site; either satisfies the threshold.
        for expr in (hom(), het()):
            assert "FMT/DP >= 30" in expr
            assert "INFO/DP >= 30" in expr

    def test_both_rules_require_support_on_both_strands(self):
        """A variant seen on one strand only is the classic short-read artefact."""
        for expr in (hom(), het()):
            assert "INFO/SAF[0] >= 2" in expr
            assert "INFO/SAR[0] >= 2" in expr

    def test_both_rules_guard_the_allele_fraction_division(self):
        for expr in (hom(), het()):
            assert "(INFO/AO[0] + INFO/RO) > 0" in expr

    def test_the_het_band_is_closed_below_and_open_above(self):
        """The boundary must belong to exactly one rule, or a call at 0.90 is both or neither."""
        assert ">= 0.1" in het() and "< 0.9" in het()
        assert ">= 0.9" in hom()

    def test_the_thresholds_come_from_the_arguments(self):
        assert ">= 0.75" in hom(hom_threshold=0.75)
        assert "FMT/DP >= 5" in het(min_dp=5)


@pytest.mark.skipif(not shutil.which("bcftools"), reason="bcftools is not installed")
class TestAgainstRealBcftools:
    """The expression has to be one bcftools accepts AND one that selects the right records."""

    @staticmethod
    def _vcf(tmp_path, rows):
        header = textwrap.dedent("""\
            ##fileformat=VCFv4.2
            ##contig=<ID=c,length=1000>
            ##INFO=<ID=DP,Number=1,Type=Integer,Description="d">
            ##INFO=<ID=RO,Number=1,Type=Integer,Description="r">
            ##INFO=<ID=AO,Number=A,Type=Integer,Description="a">
            ##INFO=<ID=SAF,Number=A,Type=Integer,Description="f">
            ##INFO=<ID=SAR,Number=A,Type=Integer,Description="r">
            ##FORMAT=<ID=GT,Number=1,Type=String,Description="g">
            ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="d">
            #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
            """)
        path = tmp_path / "in.vcf"
        path.write_text(header + "".join(rows))
        return path

    @staticmethod
    def _row(pos, dp, ro, ao, saf, sar, gt="0/1"):
        return (f"c\t{pos}\t.\tA\tT\t50\t.\t"
                f"DP={dp};RO={ro};AO={ao};SAF={saf};SAR={sar}\tGT:DP\t{gt}:{dp}\n")

    def _select(self, vcf, expr):
        out = subprocess.run(["bcftools", "view", "-H", "-i", expr, str(vcf)],
                             capture_output=True, text=True, check=True)
        return [line.split("\t")[1] for line in out.stdout.splitlines()]

    def test_the_rules_partition_calls_by_allele_fraction(self, tmp_path):
        vcf = self._vcf(tmp_path, [
            self._row(100, dp=50, ro=1, ao=49, saf=25, sar=24),    # AF 0.98 -> hom
            self._row(200, dp=50, ro=25, ao=25, saf=13, sar=12),   # AF 0.50 -> het
            self._row(300, dp=50, ro=48, ao=2, saf=1, sar=1),      # AF 0.04 -> neither
        ])
        assert self._select(vcf, hom()) == ["100"]
        assert self._select(vcf, het()) == ["200"]

    def test_a_shallow_site_passes_neither_rule(self, tmp_path):
        vcf = self._vcf(tmp_path, [self._row(100, dp=5, ro=0, ao=5, saf=3, sar=2)])
        assert self._select(vcf, hom()) == []
        assert self._select(vcf, het()) == []

    def test_a_single_strand_variant_passes_neither_rule(self, tmp_path):
        """Full depth and a clean allele fraction, but every alt read on the forward strand."""
        vcf = self._vcf(tmp_path, [self._row(100, dp=50, ro=1, ao=49, saf=49, sar=0)])
        assert self._select(vcf, hom()) == []
        assert self._select(vcf, het()) == []

    def test_a_site_with_no_observations_does_not_divide_by_zero(self, tmp_path):
        vcf = self._vcf(tmp_path, [self._row(100, dp=40, ro=0, ao=0, saf=0, sar=0)])
        assert self._select(vcf, hom()) == []
        assert self._select(vcf, het()) == []

    def test_the_boundary_fraction_is_homozygous_and_not_heterozygous(self, tmp_path):
        """Exactly 0.90: hom is inclusive, het is exclusive, so it must land in hom only."""
        vcf = self._vcf(tmp_path, [self._row(100, dp=50, ro=5, ao=45, saf=23, sar=22)])
        assert self._select(vcf, hom()) == ["100"]
        assert self._select(vcf, het()) == []


@pytest.mark.skipif(not shutil.which("bcftools"), reason="bcftools is not installed")
class TestIndelSoftFilter:
    """CALL_FREEBAYES keeps every indel and marks the ones short of the rule instead of dropping them:
    `bcftools view --types indels | bcftools filter -s LowSupport -i "<hom> || <het>"`, as in the process."""

    def test_every_indel_is_kept_and_only_the_short_ones_are_marked(self, tmp_path):
        header = TestAgainstRealBcftools._vcf(tmp_path, []).read_text()
        rows = [
            "c\t100\t.\tA\tAT\t50\t.\tDP=50;RO=1;AO=49;SAF=25;SAR=24\tGT:DP\t1/1:50\n",    # hom indel
            "c\t200\t.\tAT\tA\t50\t.\tDP=50;RO=40;AO=10;SAF=5;SAR=5\tGT:DP\t0/1:50\n",    # 20% indel
            "c\t300\t.\tAT\tA\t50\t.\tDP=50;RO=46;AO=4;SAF=4;SAR=0\tGT:DP\t0/1:50\n",     # one strand
            "c\t400\t.\tA\tT\t50\t.\tDP=50;RO=1;AO=49;SAF=25;SAR=24\tGT:DP\t1/1:50\n",     # a SNP
        ]
        vcf = tmp_path / "calls.vcf"
        vcf.write_text(header + "".join(rows))
        rule = f"{hom()} || {het()}"
        view = subprocess.run(["bcftools", "view", "--types", "indels", "-Ou", str(vcf)],
                              capture_output=True, check=True)
        out = subprocess.run(["bcftools", "filter", "-s", "LowSupport", "-i", rule, "-Ov", "-"],
                             input=view.stdout, capture_output=True, check=True)
        got = [(ln.split(b"\t")[1].decode(), ln.split(b"\t")[6].decode())
               for ln in out.stdout.splitlines() if not ln.startswith(b"#")]
        assert got == [("100", "PASS"), ("200", "PASS"), ("300", "LowSupport")]


class TestCli:
    """The processes call this as a command and substitute the result into bcftools."""

    def _run(self, *args):
        cmd = [sys.executable, str(BIN / "vcf_filter_rules.py"), *args,
               "--min-dp", "30", "--min-alt-fwd", "2", "--min-alt-rev", "2",
               "--hom-threshold", "0.90", "--het-min-frac", "0.10"]
        return subprocess.run(cmd, capture_output=True, text=True, check=True).stdout.strip()

    def test_each_rule_is_printed_on_one_line(self):
        for rule in ("hom", "het", "either", "base"):
            out = self._run(rule)
            assert out and "\n" not in out, f"{rule} did not print a single line"

    def test_either_is_the_disjunction_the_processes_use(self):
        assert self._run("either") == f"{self._run('hom')} || {self._run('het')}"

    def test_an_unknown_rule_is_rejected(self):
        cmd = [sys.executable, str(BIN / "vcf_filter_rules.py"), "nonsense",
               "--min-dp", "1", "--min-alt-fwd", "1", "--min-alt-rev", "1",
               "--hom-threshold", "0.9", "--het-min-frac", "0.1"]
        assert subprocess.run(cmd, capture_output=True).returncode != 0
