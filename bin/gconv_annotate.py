#!/usr/bin/env python3
"""What a conversion tract actually does to the genes it lands on.

A tract reported as coordinates is a finding nobody can act on. The questions asked of one are
which gene it hit, whether the copied bases change the protein, and whether it moved an allele
into a gene where alleles matter. All of that is a GFF and the reference sequence away, and the
pipeline already carries both.

The tract is a run of DIAGNOSTIC sites, and only those change: everywhere else the two copies are
identical, so a conversion there is invisible and also inconsequential. So the consequence of a
tract is the consequence of substituting the donor's base at each diagnostic site under it.

Each is scored ALONE against the reference codon, not compounded with its neighbours. Two changes
in one codon are rare at the density paralogs differ, and reporting a joint effect would mean
claiming a phase for the conversion that the breakpoints do not resolve. One site, one codon, one
statement a reader can check.
"""

from __future__ import annotations

import gzip
import re

# The genetic code in its canonical TCAG ordering, which is the form every reference prints it in
# and therefore the form that can be checked by eye. Bacterial table 11 differs from this only in
# which codons may START a protein; the sense codons and the stops are identical, and a conversion
# lands in the middle of a gene.
_BASES = "TCAG"
_AAS = ("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRR"
        "VVVVAAAADDEEGGGG")
CODONS = {a + b + c: aa for (a, b, c), aa in
          zip(((x, y, z) for x in _BASES for y in _BASES for z in _BASES), _AAS)}

COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def _attr(field, key):
    m = re.search(rf"(?:^|;)\s*{re.escape(key)}=([^;]*)", field)
    return m.group(1).strip() if m else None


def parse_cds(path):
    """Coding features from a GFF3: [{contig, start, end, strand, phase, name, locus_tag}].

    Only CDS, because only a CDS has a reading frame and the question here is what happens to the
    protein. `gene` features would give a name over a wider span and no way to place a codon.

    Kept as separate records rather than merged by name: a gene split across several CDS lines is
    a gene with several segments, and a site falls in one of them.
    """
    out = []
    if not path:
        return out
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#") or "\t" not in line:
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9 or c[2] != "CDS":
                continue
            try:
                start, end = int(c[3]), int(c[4])
            except ValueError:
                continue
            phase = c[7] if c[7] in ("0", "1", "2") else "0"
            locus = _attr(c[8], "locus_tag")
            name = _attr(c[8], "gene") or _attr(c[8], "Name") or locus or _attr(c[8], "ID")
            out.append({"contig": c[0], "start": min(start, end), "end": max(start, end),
                        "strand": "-" if c[6] == "-" else "+", "phase": int(phase),
                        "name": name or f"{c[0]}:{start}", "locus_tag": locus or ""})
    out.sort(key=lambda f: (f["contig"], f["start"]))
    return out


def read_fasta(path):
    """{contig: sequence}, keyed on the first word of the header as every other tool does."""
    seqs, name, chunks = {}, None, []
    if not path:
        return seqs
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks)
                name, chunks = line[1:].split()[0] if line[1:].split() else "", []
            elif name is not None:
                chunks.append(line.strip())
    if name is not None:
        seqs[name] = "".join(chunks)
    return seqs


def features_at(features, contig, start, end):
    """Every CDS overlapping [start, end]. Linear: a tract spans a handful of genes at most."""
    return [f for f in features
            if f["contig"] == contig and f["start"] <= end and f["end"] >= start]


def codon_change(seq, feature, pos, alt):
    """What substituting `alt` at reference position `pos` does to this CDS.

    Returns (codon_number, ref_aa, alt_aa) with the codon numbered from 1 at the start of the
    feature, or None where the question does not arise: outside the feature, off the end of the
    contig, or a codon running past it.

    The phase is honoured, so a CDS whose first base is not the first base of a codon is read in
    the frame the annotation declares rather than the one its coordinates suggest.
    """
    if not (feature["start"] <= pos <= feature["end"]) or len(alt) != 1:
        return None
    if feature["strand"] == "+":
        offset = pos - (feature["start"] + feature["phase"])
    else:
        offset = (feature["end"] - feature["phase"]) - pos
    if offset < 0:
        return None                                   # inside the phase, before the first codon
    index = offset % 3                                # which base of its codon this position is

    if feature["strand"] == "+":
        first = pos - index
        codon = seq[first - 1:first + 2]
        alt_base = alt
    else:
        first = pos + index                           # the codon's first base, in genome order
        codon = seq[first - 3:first].translate(COMPLEMENT)[::-1]
        alt_base = alt.translate(COMPLEMENT)
    if len(codon) != 3:
        return None

    codon = codon.upper()
    mutated = codon[:index] + alt_base.upper() + codon[index + 1:]
    ref_aa, alt_aa = CODONS.get(codon), CODONS.get(mutated)
    if ref_aa is None or alt_aa is None:
        return None                                   # an N or another ambiguity in the codon
    return offset // 3 + 1, ref_aa, alt_aa


def annotate_tract(features, seqs, contig, start, end, changes):
    """The genes a tract covers and what its copied bases do to them.

    `changes` is {acceptor position: donor base} for the diagnostic sites under the tract, which
    is exactly what a conversion writes there.

    A site in more than one CDS is reported once per CDS, because overlapping reading frames are
    real and the same base is a different codon in each.
    """
    genes, syn, nonsyn, effects = [], 0, 0, []
    seq = seqs.get(contig, "")
    for f in features_at(features, contig, start, end):
        label = f["locus_tag"] or f["name"]
        if label not in genes:
            genes.append(label)
        if not seq:
            continue
        for pos in sorted(changes):
            got = codon_change(seq, f, pos, changes[pos])
            if got is None:
                continue
            n, ref_aa, alt_aa = got
            if ref_aa == alt_aa:
                syn += 1
            else:
                nonsyn += 1
                effects.append(f"{label}:{ref_aa}{n}{alt_aa}")
    return {"genes": ",".join(genes), "n_syn": syn, "n_nonsyn": nonsyn,
            "aa_changes": ",".join(effects)}


ANNOTATION_COLUMNS = ["genes", "n_syn", "n_nonsyn", "aa_changes"]
