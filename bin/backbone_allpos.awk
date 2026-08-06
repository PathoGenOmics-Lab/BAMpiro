# Turn a samtools mpileup stream into the all-positions backbone VCF body.
#
#   samtools mpileup -aa ... | awk -v MINCOV=30 -v GAPDP=7 -v BAQ=1 -v OFS='\t' \
#                                  -f backbone_allpos.awk
#
#   MINCOV  params.allpos_min_cov - depth at or above this is a confident reference call
#   GAPDP   params.consensus_min_dp - only read when BAQ=1
#   BAQ     1 when the input carries a second, BAQ-adjusted depth in column 10 (the two pileups
#           pasted together), 0 when it is a plain single mpileup
#
# The backbone supplies one record per reference position, which is what makes an all-sites VCF and
# therefore a phylogenetic supermatrix possible. Both branches of CALL_BACKBONE used to carry their
# own copy of count_bases; tests/unit/test_backbone_allpos.py covers this one.

# Count reference and alternate observations in an mpileup base string.
#
# The format is not one character per read: a read start is '^' followed by a mapping-quality
# character that must not be counted, and an indel is a sign, a length, and that many bases which
# belong to the PREVIOUS position. Getting this wrong silently inflates or deflates the depth
# evidence the consensus relies on.
function count_bases(s,   i, n, c, num) {
  RO_G = 0; AO_G = 0; i = 1; n = length(s)
  while (i <= n) {
    c = substr(s, i, 1)
    if (c == "^") { i += 2 }                                  # ^ + mapping-quality char (read start)
    else if (c == "$") { i++ }                                # read end
    else if (c == "+" || c == "-") {                          # indel: skip sign + count + that many bases
      i++; num = ""
      while (i <= n && substr(s, i, 1) ~ /[0-9]/) { num = num substr(s, i, 1); i++ }
      i += (num + 0)
    }
    else if (c == "." || c == ",") { RO_G++; i++ }            # match to reference
    else if (c ~ /[ACGTacgt]/)     { AO_G++; i++ }            # mismatch = alt observation
    else { i++ }                                              # * (del), <>, N, etc. -> ignore
  }
}

{
  chrom = $1; pos = $2; ref = $3; dp = ($4 + 0)
  count_bases($5)                                             # RO/AO from the no-BAQ pileup (real base evidence)

  wt = (dp >= MINCOV ? 1 : 0)
  nc = (dp >= MINCOV ? 0 : 1)
  gt = (dp >= MINCOV ? "0" : "./.")

  # A site is flagged baq_dropout when it is well covered but BAQ collapses it to at or below
  # consensus_min_dp (the indel-adjacent homopolymer artefact), so the consensus masks it with X
  # rather than leaving a gap. Without the BAQ column there is nothing to compare and no flag.
  flt = "."
  if (BAQ == 1) { dp_baq = ($10 + 0); if (dp > GAPDP && dp_baq <= GAPDP) flt = "baq_dropout" }

  info = "ADP=" dp ";WT=" wt ";HET=0;HOM=0;NC=" nc
  print chrom, pos, ".", ref, ".", ".", flt, info, "GT:DP:AD", gt ":" dp ":" RO_G "," AO_G
}
