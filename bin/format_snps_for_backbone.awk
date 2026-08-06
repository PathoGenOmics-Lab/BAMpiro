# Rewrite a FreeBayes SNP VCF into the reduced form the backbone merge and the consensus step
# expect: INFO becomes ADP/WT/HET/HOM/NC and the sample column becomes GT:DP.
#
#   awk -v HETMIN=0.10 -f format_snps_for_backbone.awk valid_snps.vcf
#
# HETMIN is params.het_min_frac. Both FreeBayes processes pipe through this, so the genotype
# decision below lives in exactly one place; tests/unit/test_format_snps_for_backbone.py covers it.

BEGIN { OFS = "\t" }

/^##/    { print; next }
/^#CHROM/ { print; next }

{
  dp = 0; ro = 0; ao = 0

  # INFO is tokenised on ';' and each key anchored (^DP=/^RO=/^AO=) so DP is never grabbed from a
  # look-alike field such as MQMDP= or ADP=. AO is per-alt; after `bcftools norm -m -` the record is
  # biallelic, so the first value is the one that matters.
  ni = split($8, I, ";")
  for (k = 1; k <= ni; k++) {
    if (I[k] ~ /^DP=/)      { v = I[k]; sub(/^DP=/, "", v); dp = v + 0 }
    else if (I[k] ~ /^RO=/) { v = I[k]; sub(/^RO=/, "", v); ro = v + 0 }
    else if (I[k] ~ /^AO=/) { v = I[k]; sub(/^AO=/, "", v); split(v, ax, ","); ao = ax[1] + 0 }
  }

  # Depth falls back through FORMAT/DP, then RO+AO, then 1 so the consensus never divides by zero.
  if (dp == 0) { n = split($9, fmt, ":"); m = split($10, dat, ":"); for (i = 1; i <= n; i++) if (fmt[i] == "DP") dp = dat[i] }
  if (dp == 0) dp = ro + ao
  if (dp == 0) dp = 1

  split($10, b, ":"); gtype = b[1]

  # Re-validate het: `bcftools norm -m -` splits a multiallelic 1/2 into biallelic records whose
  # reference index 0 is an ARTIFACT (RO~0). A het whose reference allele is essentially
  # unsupported is not a real het -> homozygous-alt (matches the af-based hom/het rule and stops
  # the reference base leaking into the consensus IUPAC).
  if ((gtype == "0/1" || gtype == "1/0") && (ro + ao) > 0 && ro / (ro + ao) < HETMIN) gtype = "1/1"

  het = 0; hom = 0
  if (gtype == "0/1" || gtype == "1/0") het = 1
  else if (gtype == "1/1") hom = 1

  $8  = "ADP=" dp ";WT=0;HET=" het ";HOM=" hom ";NC=0"
  $9  = "GT:DP"
  $10 = gtype ":" dp
  print
}
