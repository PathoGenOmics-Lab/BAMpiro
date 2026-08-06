#!/usr/bin/env python3
"""Print the bcftools filter expression that defines a homozygous or heterozygous SNP call.

    vcf_filter_rules.py hom --min-dp 30 --min-alt-fwd 2 --min-alt-rev 2 \
                            --hom-threshold 0.90 --het-min-frac 0.10

This is the scientific definition of a call, not plumbing: which depth and strand support a site
needs, and where the allele-fraction boundary between heterozygous and homozygous sits. It lived
inline in two FreeBayes processes, written out twice, which is how the two drift apart. Keeping it
here means one definition and tests/unit/test_vcf_filter_rules.py covering it.

The expressions are consumed as `bcftools view -i "<expr>"`.
"""

from __future__ import annotations

import argparse
import sys

# Depth may be recorded per-sample or per-site depending on the caller, so either satisfies it.
# SAF/SAR are the forward and reverse alt-supporting reads: requiring both rejects a variant seen
# on one strand only, the classic short-read artefact.
BASE = ("(FMT/DP >= {min_dp} || INFO/DP >= {min_dp})"
        " && INFO/SAF[0] >= {min_alt_fwd} && INFO/SAR[0] >= {min_alt_rev}")

# Guard the division: a site with neither ref nor alt observations has no defined allele fraction.
AF_GUARD = "(INFO/AO[0] + INFO/RO) > 0"
AF_CALC = "((INFO/AO[0] * 1.0) / (INFO/AO[0] + INFO/RO))"


def base_filter(min_dp, min_alt_fwd, min_alt_rev) -> str:
    return BASE.format(min_dp=min_dp, min_alt_fwd=min_alt_fwd, min_alt_rev=min_alt_rev)


def hom_rule(min_dp, min_alt_fwd, min_alt_rev, hom_threshold) -> str:
    """Alt allele at or above the homozygous fraction."""
    return (f"{base_filter(min_dp, min_alt_fwd, min_alt_rev)} && {AF_GUARD}"
            f" && {AF_CALC} >= {hom_threshold}")


def het_rule(min_dp, min_alt_fwd, min_alt_rev, hom_threshold, het_min_frac) -> str:
    """Alt allele between the heterozygous floor and the homozygous threshold, upper bound exclusive."""
    return (f"{base_filter(min_dp, min_alt_fwd, min_alt_rev)} && {AF_GUARD}"
            f" && {AF_CALC} >= {het_min_frac} && {AF_CALC} < {hom_threshold}")


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("rule", choices=["hom", "het", "either", "base"])
    p.add_argument("--min-dp", required=True)
    p.add_argument("--min-alt-fwd", required=True)
    p.add_argument("--min-alt-rev", required=True)
    p.add_argument("--hom-threshold", required=True)
    p.add_argument("--het-min-frac", required=True)
    return p.parse_args(argv)


def main(argv=None) -> int:
    a = parse_args(argv)
    hom = hom_rule(a.min_dp, a.min_alt_fwd, a.min_alt_rev, a.hom_threshold)
    het = het_rule(a.min_dp, a.min_alt_fwd, a.min_alt_rev, a.hom_threshold, a.het_min_frac)
    out = {
        "hom": hom,
        "het": het,
        "either": f"{hom} || {het}",
        "base": base_filter(a.min_dp, a.min_alt_fwd, a.min_alt_rev),
    }[a.rule]
    sys.stdout.write(out + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
