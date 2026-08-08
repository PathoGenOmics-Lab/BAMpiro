"""Unit tests for bin/gconv_annotate.py (what a conversion tract does to the genes it lands on).

The arithmetic here is codon arithmetic, and it fails silently: a frame off by one, or a strand
read the wrong way, returns a perfectly plausible amino acid change at a perfectly plausible
position. So the reverse strand is tested against the same sequence read forwards, the frame is
tested against a phase the coordinates alone would get wrong, and the genetic code is checked
against a table written out independently rather than trusted because it has 64 entries.
"""

from __future__ import annotations

import pytest

from conftest import load_script

ga = load_script("gconv_annotate")

CONTIG = "chr"


@pytest.fixture
def write(tmp_path):
    """These parsers take a path, so their fixtures have to be files."""
    def _write(text, name):
        path = tmp_path / name
        path.write_text(text)
        return str(path)
    return _write


# --------------------------------------------------------------------------- the genetic code


def test_the_genetic_code_is_the_genetic_code():
    """Written out here from a different source than the one under test.

    A table built by zipping two strings is one transposition away from being wrong in a way that
    no other test would catch: every codon still maps to some amino acid, and the two or three
    that moved are the ones nobody looks up.
    """
    known = {
        "ATG": "M", "TGG": "W", "TAA": "*", "TAG": "*", "TGA": "*",
        "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTG": "L",
        "ATT": "I", "ATC": "I", "ATA": "I", "GTG": "V", "GTA": "V",
        "TCT": "S", "AGT": "S", "AGC": "S", "CCG": "P", "ACT": "T",
        "GCA": "A", "TAT": "Y", "CAT": "H", "CAA": "Q", "AAT": "N",
        "AAA": "K", "GAT": "D", "GAA": "E", "TGT": "C", "CGA": "R",
        "AGA": "R", "GGG": "G",
    }

    assert len(ga.CODONS) == 64
    for codon, aa in known.items():
        assert ga.CODONS[codon] == aa, codon


# --------------------------------------------------------------------------- GFF and FASTA


GFF = """\
##gff-version 3
chr\t.\tgene\t100\t199\t.\t+\t.\tID=gene1;Name=ppeA;locus_tag=Rv0001
chr\t.\tCDS\t100\t199\t.\t+\t0\tID=cds1;gene=ppeA;locus_tag=Rv0001
chr\t.\tCDS\t300\t399\t.\t-\t0\tID=cds2;gene=ppeB;locus_tag=Rv0002
chr\t.\ttRNA\t500\t560\t.\t+\t.\tID=t1;locus_tag=Rv0003
"""


def test_only_coding_features_are_read(write):
    """A gene line gives a name over a wider span and no reading frame to place a codon in, and a
    tRNA has no protein to change. Neither can answer the question this asks."""
    feats = ga.parse_cds(write(GFF, "x.gff"))

    assert [f["locus_tag"] for f in feats] == ["Rv0001", "Rv0002"]
    assert [f["strand"] for f in feats] == ["+", "-"]
    assert feats[0]["name"] == "ppeA"


@pytest.mark.parametrize("attrs,name", [
    ("gene=ppeA;Name=other;locus_tag=Rv1;ID=cds1", "ppeA"),
    ("Name=other;locus_tag=Rv1;ID=cds1", "other"),
    ("locus_tag=Rv1;ID=cds1", "Rv1"),
    ("ID=orphan", "orphan"),
])
def test_a_feature_is_named_by_the_best_attribute_it_carries(write, attrs, name):
    """GFFs from different sources put the readable name in different places, so the fallbacks
    are tried in order. Collapsing that order means a gene called by its ID when it had a name."""
    feats = ga.parse_cds(write(f"chr\t.\tCDS\t10\t20\t.\t+\t0\t{attrs}\n", "y.gff"))

    assert feats[0]["name"] == name


def test_a_feature_without_a_locus_tag_has_an_empty_one(write):
    feats = ga.parse_cds(write("chr\t.\tCDS\t10\t20\t.\t+\t0\tID=orphan\n", "n.gff"))

    assert feats[0]["locus_tag"] == ""


def test_a_reversed_feature_is_stored_ascending(write):
    """A range is a range whichever way the annotation wrote it, and every lookup below assumes
    start <= end."""
    feats = ga.parse_cds(write("chr\t.\tCDS\t200\t100\t.\t-\t0\tlocus_tag=Rv9\n", "z.gff"))

    assert (feats[0]["start"], feats[0]["end"]) == (100, 200)


def test_fasta_is_keyed_on_the_first_word_of_the_header(write):
    """Which is what every other tool in the pipeline does, so `chr` in a BAM finds `chr` here
    even when the FASTA header carries a description after it."""
    seqs = ga.read_fasta(write(">chr some description here\nACGT\nACGT\n>other\nTTTT\n", "r.fa"))

    assert seqs == {"chr": "ACGTACGT", "other": "TTTT"}


# --------------------------------------------------------------------------- codons


def _cds(start, end, strand="+", phase=0, tag="Rv0001"):
    return {"contig": CONTIG, "start": start, "end": end, "strand": strand,
            "phase": phase, "name": tag, "locus_tag": tag}


# ATG AAA TTT TGG TAA, starting at position 11 of the sequence below
SEQ = "NNNNNNNNNN" + "ATGAAATTTTGGTAA"


def test_a_substitution_is_placed_in_the_right_codon_and_the_right_position_in_it():
    """Every base of the second codon in turn, so a frame off by one shows up as the wrong codon
    number and a base index off by one as the wrong amino acid."""
    cds = _cds(11, 25)

    assert ga.codon_change(SEQ, cds, 14, "C") == (2, "K", "Q")     # AAA -> CAA
    assert ga.codon_change(SEQ, cds, 15, "G") == (2, "K", "R")     # AAA -> AGA
    assert ga.codon_change(SEQ, cds, 16, "G") == (2, "K", "K")     # AAA -> AAG, synonymous


def test_a_change_that_creates_a_stop_is_reported_as_one():
    assert ga.codon_change(SEQ, _cds(11, 25), 22, "A") == (4, "W", "*")   # TGG -> TAG


def test_the_reverse_strand_is_read_on_the_reverse_strand():
    """The same bases, annotated the other way round, are a different protein entirely.

    Read forwards, positions 11-25 are M K F W *. Read backwards they are the reverse complement,
    and a tool that gets this wrong still returns a plausible amino acid at a plausible codon.
    """
    forward = ga.codon_change(SEQ, _cds(11, 25, "+"), 14, "C")
    reverse = ga.codon_change(SEQ, _cds(11, 25, "-"), 14, "C")

    assert forward == (2, "K", "Q")
    # TTACCAAAATTTCAT is the reverse complement; position 14 is the 12th base of it, so the last
    # base of codon 4, and the substituted base is complemented too.
    assert reverse == (4, "F", "L")
    assert forward != reverse


@pytest.mark.parametrize("phase,expected", [(0, (2, "K", "Q")), (1, (1, "*", "C")),
                                            (2, (1, "E", "A"))])
def test_the_phase_the_annotation_declares_is_the_phase_that_is_used(phase, expected):
    """A CDS whose first base is not the first base of a codon is common in a GFF, and reading the
    frame off the coordinates instead puts every codon in the feature one position out.

    One position, three phases, three different answers. Each is a plausible amino acid change at
    a plausible codon, which is why getting this wrong is not something an eye catches downstream.
    """
    assert ga.codon_change(SEQ, _cds(11, 25, phase=phase), 14, "C") == expected


@pytest.mark.parametrize("pos,inside", [(10, False), (11, True), (25, True), (26, False)])
def test_a_feature_holds_its_own_first_and_last_base(pos, inside):
    """Inclusive at both ends. One step in either direction and the codons at the edges of every
    gene are either lost or invented, and the ends of genes are where paralogs differ most."""
    assert (ga.codon_change(SEQ, _cds(11, 25), pos, "C") is not None) is inside


def test_a_position_past_a_short_gene_is_not_read_as_a_later_codon_of_it():
    """The bases after a gene ends are somebody else's, or nobody's.

    Dropping the containment check does not fail loudly: the frame keeps counting past the stop
    codon and returns a codon number and an amino acid change that look exactly like a real one.
    Only a position with sequence AROUND it shows this, which is why the cases at the very ends
    of the contig did not.
    """
    short = _cds(11, 16)                       # ends four bases before the position asked about

    assert ga.codon_change(SEQ, short, 20, "C") is None


def test_the_first_base_of_the_first_codon_is_in_the_first_codon():
    """Offset zero is the start of the reading frame, not a position before it."""
    assert ga.codon_change(SEQ, _cds(11, 25), 11, "C") == (1, "M", "L")   # ATG -> CTG


def test_only_a_single_base_substitution_is_scored():
    """A deletion marker is not a substitution and the answer to what it does is a frameshift
    question, not a codon one. Half an answer there is worse than none."""
    assert ga.codon_change(SEQ, _cds(11, 25), 14, "AC") is None
    assert ga.codon_change(SEQ, _cds(11, 25), 14, "-") is None


def test_a_codon_running_off_the_end_of_the_contig_is_not_guessed_at():
    assert ga.codon_change("ATGAA", _cds(1, 6), 5, "C") is None


def test_an_ambiguous_base_in_the_codon_gives_no_answer():
    """Either side of the comparison. An N in the reference codon leaves nothing to change FROM,
    and one in the substituted base leaves nothing to change TO."""
    assert ga.codon_change("ATGANA", _cds(1, 6), 6, "C") is None
    assert ga.codon_change("ATGAAA", _cds(1, 6), 6, "N") is None


# --------------------------------------------------------------------------- whole tracts


def test_a_tract_reports_the_genes_it_covers_and_what_it_does_to_them():
    features = [_cds(11, 25, tag="Rv0001")]
    changes = {14: "C", 16: "G", 22: "A"}          # nonsynonymous, synonymous, a stop

    got = ga.annotate_tract(features, {CONTIG: SEQ}, CONTIG, 11, 25, changes)

    assert got["genes"] == "Rv0001"
    assert (got["n_syn"], got["n_nonsyn"]) == (1, 2)
    assert got["aa_changes"] == "Rv0001:K2Q,Rv0001:W4*"


def test_a_tract_spanning_two_genes_names_both():
    features = [_cds(11, 16, tag="Rv0001"), _cds(17, 25, tag="Rv0002")]

    got = ga.annotate_tract(features, {CONTIG: SEQ}, CONTIG, 11, 25, {14: "C"})

    assert got["genes"] == "Rv0001,Rv0002"


def test_a_tract_in_no_gene_at_all_says_so_without_failing():
    got = ga.annotate_tract([_cds(500, 600)], {CONTIG: SEQ}, CONTIG, 11, 25, {14: "C"})

    assert got == {"genes": "", "n_syn": 0, "n_nonsyn": 0, "aa_changes": ""}


def test_a_site_in_two_overlapping_reading_frames_is_counted_in_each():
    """Overlapping genes are real, and the same base is a different codon in each of them. One
    report per frame is the honest answer; picking one of the two is not."""
    features = [_cds(11, 25, tag="Rv0001"), _cds(12, 25, tag="Rv0002")]

    got = ga.annotate_tract(features, {CONTIG: SEQ}, CONTIG, 11, 25, {14: "C"})

    assert got["genes"] == "Rv0001,Rv0002"
    assert got["n_syn"] + got["n_nonsyn"] == 2


def test_without_the_reference_sequence_the_genes_are_still_named():
    """A GFF alone answers "where did it land", which is most of the question, and a codon needs
    bases nobody has here. Half an answer beats refusing to give one."""
    got = ga.annotate_tract([_cds(11, 25)], {}, CONTIG, 11, 25, {14: "C"})

    assert got["genes"] == "Rv0001"
    assert (got["n_syn"], got["n_nonsyn"], got["aa_changes"]) == (0, 0, "")


@pytest.mark.parametrize("start,end,covered", [(1, 10, False), (1, 11, True),
                                               (25, 30, True), (26, 30, False)])
def test_a_feature_is_covered_when_it_touches_the_tract_at_all(start, end, covered):
    """A tract reaching a gene's first base has reached the gene."""
    got = ga.annotate_tract([_cds(start, end)], {}, CONTIG, 11, 25, {})

    assert bool(got["genes"]) is covered


def test_a_gene_with_no_locus_tag_is_reported_by_its_name():
    """Not every annotation carries locus tags, and a blank column is not a report."""
    feature = dict(_cds(11, 25), locus_tag="", name="ppeA")

    got = ga.annotate_tract([feature], {CONTIG: SEQ}, CONTIG, 11, 25, {14: "C"})

    assert got["genes"] == "ppeA"
    assert got["aa_changes"] == "ppeA:K2Q"


def test_a_feature_on_another_contig_is_not_reported():
    features = [dict(_cds(11, 25), contig="other")]

    assert ga.annotate_tract(features, {CONTIG: SEQ}, CONTIG, 11, 25, {14: "C"})["genes"] == ""
