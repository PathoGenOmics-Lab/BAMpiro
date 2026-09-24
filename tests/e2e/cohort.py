#!/usr/bin/env python3
"""The end-to-end cohort: a genome, three samples, and everything each of them carries.

The stub runs prove the wiring and the unit tests prove each script, but neither runs a process's
script block with the real tools on reads. This cohort is simulated so that what the pipeline
must report about it is known in advance: which SNPs are fixed and which are minorities, which
two SNPs of a codon sit on the same reads and which only share the codon, which indels a sample
carries and at what fraction, a stretch one sample has lost, and a repeat no read can be placed
in. tests/e2e/test_e2e.py runs the pipeline on it and holds the outputs against that truth.

    python3 tests/e2e/cohort.py OUTDIR

writes the reference, its GFF, the reads, a samplesheet with absolute paths and truth.json.
Everything derives from one seed, so every machine gets the same files.

The samples:

    S1  single-end, two runs merged. Fixed SNPs; two SNPs of one codon of a minus-strand gene on
        the same reads (Phe>Met; each base alone says Ile or Leu); a fixed in-frame deletion of
        one codon; and S2's frameshift at 20%, read on the forward strand only, which a call must
        keep but mark LowSupport.
    S2  paired-end. Fixed SNPs; two SNPs of one codon on the same reads (His>Glu; each base alone
        says Asp or Gln); a fixed frameshift insertion; and a SNP in one copy of a repeat, which
        no read can be placed in and so must not be called.
    S3  paired-end, a mixture. Fixed SNPs and a lost stretch of 1.5 kb holding a whole gene. At
        25%, a SNP, two SNPs of one codon on the same reads (Glu>Leu, where the first base alone
        says a stop) and a frameshift deletion. At 30% and 25%, two SNPs of one codon carried by
        different molecules: no read has both, so no codon change may be reported.
"""

from __future__ import annotations

import gzip
import json
import random
import sys
from pathlib import Path

SEED = 20260924
CONTIG = "e2e_chr"
REF_ID = "E2EREF"
GENOME_LEN = 24_000
READ_LEN = 150
INSERT_MEAN, INSERT_SD = 350, 30
COVERAGE = 100
ERROR_RATE = 0.001
QUAL = "I" * READ_LEN

BASES = "ACGT"
COMPLEMENT = str.maketrans("ACGT", "TGCA")

CODON_TABLE = {
    "TTT": "Phe", "TTC": "Phe", "TTA": "Leu", "TTG": "Leu", "CTT": "Leu", "CTC": "Leu", "CTA": "Leu",
    "CTG": "Leu", "ATT": "Ile", "ATC": "Ile", "ATA": "Ile", "ATG": "Met", "GTT": "Val", "GTC": "Val",
    "GTA": "Val", "GTG": "Val", "TCT": "Ser", "TCC": "Ser", "TCA": "Ser", "TCG": "Ser", "CCT": "Pro",
    "CCC": "Pro", "CCA": "Pro", "CCG": "Pro", "ACT": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr",
    "GCT": "Ala", "GCC": "Ala", "GCA": "Ala", "GCG": "Ala", "TAT": "Tyr", "TAC": "Tyr", "TAA": "*",
    "TAG": "*", "CAT": "His", "CAC": "His", "CAA": "Gln", "CAG": "Gln", "AAT": "Asn", "AAC": "Asn",
    "AAA": "Lys", "AAG": "Lys", "GAT": "Asp", "GAC": "Asp", "GAA": "Glu", "GAG": "Glu", "TGT": "Cys",
    "TGC": "Cys", "TGA": "*", "TGG": "Trp", "CGT": "Arg", "CGC": "Arg", "CGA": "Arg", "CGG": "Arg",
    "AGT": "Ser", "AGC": "Ser", "AGA": "Arg", "AGG": "Arg", "GGT": "Gly", "GGC": "Gly", "GGA": "Gly",
    "GGG": "Gly",
}
SENSE = sorted(c for c, aa in CODON_TABLE.items() if aa != "*")

# (locus_tag, name, first base, codons including start and stop, strand)
GENES = [
    ("E2E_0001", "geneA", 1001, 800, "+"),
    ("E2E_0002", "geneB", 4001, 200, "+"),
    ("E2E_0003", "geneC", 6001, 800, "-"),
    ("E2E_0004", "geneD", 10001, 400, "+"),
    ("E2E_0005", "geneE", 20201, 300, "+"),   # inside the stretch S3 has lost
]

# Codons set by hand, so each planted change has a known meaning and an indel has only one place
# it can be written (neither end of it repeats the base beside it, so no aligner can shift it).
FIXED_CODONS = {
    "geneA": {145: "CAC", 260: "GAG", 380: "CTG"},
    "geneB": {49: "GAT", 101: "AGC", 134: "TCT", 135: "CGA", 136: "GTT"},
    "geneC": {100: "TTT", 500: "TGG"},
}

# The first copy of a repeat, copied again at REPEAT_COPY: reads from either copy map equally well
# to both, so none of them has the mapping quality a call needs.
REPEAT = (13001, 14500)
REPEAT_COPY = 17001
LOST = (20001, 21500)   # the stretch S3 does not have

# Fixed SNPs by who carries them. None falls in a repeat, in the lost stretch or near a gene edge.
FIXED_SNPS = {
    ("S1",): [3601, 10151, 15201],
    ("S2",): [5001, 10451, 12001, 19001],
    ("S3",): [8801, 22601],
    ("S1", "S2"): [9201, 10751],
    ("S1", "S2", "S3"): [11601],
}
REPEAT_SNP = REPEAT_COPY + 500   # S2's SNP in the second copy of the repeat

TRANSITION = {"A": "G", "G": "A", "C": "T", "T": "C"}


def revcomp(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


# ---------------------------------------------------------------------------------------------
# The reference
# ---------------------------------------------------------------------------------------------

def build_reference(rng: random.Random):
    """The genome and each gene's coding sequence, with the fixed codons in place."""
    genome = [rng.choice(BASES) for _ in range(GENOME_LEN)]
    cds = {}
    for locus, name, start, n_codons, strand in GENES:
        fixed = FIXED_CODONS.get(name, {})
        codons = ["ATG"] + [fixed.get(k) or rng.choice(SENSE) for k in range(2, n_codons)] + ["TAA"]
        seq = "".join(codons)
        cds[name] = seq
        genome[start - 1:start - 1 + len(seq)] = list(seq if strand == "+" else revcomp(seq))
    a, b = REPEAT
    genome[REPEAT_COPY - 1:REPEAT_COPY - 1 + (b - a + 1)] = genome[a - 1:b]
    return "".join(genome), cds


def gene(name):
    return next(g for g in GENES if g[1] == name)


def codon_positions(name, k):
    """Genomic positions (1-based) of codon k's three bases, in coding order."""
    _, _, start, n_codons, strand = gene(name)
    first = 3 * (k - 1)
    if strand == "+":
        return [start + first + j for j in range(3)]
    end = start + 3 * n_codons - 1
    return [end - first - j for j in range(3)]


def gene_at(pos):
    return next((g[1] for g in GENES if g[2] <= pos <= g[2] + 3 * g[3] - 1), None)


def translate(seq):
    return [CODON_TABLE.get(seq[i:i + 3], "?") for i in range(0, len(seq) - 2, 3)]


def protein_change(genome, pos, ref, alt):
    """(gene, HGVS.p) of one change as SnpEff writes it: the first amino acid it changes, then
    the new one, `fs` for a frameshift or `del` for a codon removed. (None, None) outside a gene."""
    name = gene_at(pos)
    if name is None:
        return None, None
    _, _, start, n_codons, strand = gene(name)

    def coding(seq):
        window = seq[start - 1:start - 1 + 3 * n_codons]
        return window if strand == "+" else revcomp(window)

    edited = genome[:pos - 1] + alt + genome[pos - 1 + len(ref):]
    before, after = translate(coding(genome)), translate(coding(edited))
    k = next((i for i, (a, b) in enumerate(zip(before, after)) if a != b), None)
    shift = len(alt) - len(ref)
    if shift == 0 and k is None:   # synonymous: the codon holding the base
        k = ((pos - start) if strand == "+" else (start + 3 * n_codons - 1 - pos)) // 3
        return name, f"p.{before[k]}{k + 1}{before[k]}"
    if shift == 0:
        return name, f"p.{before[k]}{k + 1}{after[k]}"
    return name, f"p.{before[k]}{k + 1}" + ("fs" if shift % 3 else "del")


def codon_change(genome, name, k, new_codon):
    """The SNPs turning codon k of a gene into new_codon: [(pos, ref, alt)] in genome terms."""
    strand = gene(name)[4]
    out = []
    for pos, base in zip(codon_positions(name, k), new_codon):
        alt = base if strand == "+" else base.translate(COMPLEMENT)
        if genome[pos - 1] != alt:
            out.append((pos, genome[pos - 1], alt))
    return out


# ---------------------------------------------------------------------------------------------
# What each sample carries
# ---------------------------------------------------------------------------------------------

def plan(genome, cds):
    """Every planted change, and each sample as a mixture of haplotypes that carry them."""
    snp = lambda pos: (pos, genome[pos - 1], TRANSITION[genome[pos - 1]])   # noqa: E731

    fixed = {s: [] for s in ("S1", "S2", "S3")}
    for carriers, positions in FIXED_SNPS.items():
        for pos in positions:
            for s in carriers:
                fixed[s].append(snp(pos))

    mnv = {
        "S1": ("geneC", 100, "ATG"),    # TTT>ATG Phe>Met: alone, ATT Ile and TTG Leu
        "S2": ("geneA", 145, "GAG"),    # CAC>GAG His>Glu: alone, GAC Asp and CAG Gln
        "S3": ("geneA", 260, "TTG"),    # GAG>TTG Glu>Leu: alone, TAG stop and GTG Val
    }
    mnv_snps = {s: codon_change(genome, *m) for s, m in mnv.items()}
    split_1 = codon_change(genome, "geneA", 380, "ATG")   # CTG>ATG, on some molecules
    split_2 = codon_change(genome, "geneA", 380, "CTC")   # CTG>CTC, on others
    # geneC codon 500 is TGG (Trp); TGG>CGG (Arg) is a plain missense change at one base.
    minority_snp = codon_change(genome, "geneC", 500, "CGG")

    # Indels in genome terms: (anchor position, REF, ALT) as a normalised VCF writes them.
    b49 = codon_positions("geneB", 49)       # GAT: insert C after GA -> GAC T..., a frameshift
    ins = (b49[1], genome[b49[1] - 1], genome[b49[1] - 1] + "C")
    b101 = codon_positions("geneB", 101)     # AGC: delete its G, a frameshift
    del1 = (b101[0], genome[b101[0] - 1:b101[1]], genome[b101[0] - 1])
    b135 = codon_positions("geneB", 135)     # CGA (Arg), between TCT (Ser) and GTT (Val)
    anchor = b135[0] - 1
    del3 = (anchor, genome[anchor - 1:b135[2]], genome[anchor - 1])

    s1 = fixed["S1"] + mnv_snps["S1"]
    s2 = fixed["S2"] + mnv_snps["S2"] + [snp(REPEAT_SNP)]
    s3 = fixed["S3"]
    samples = {
        # name: (layout, runs, [(weight, snps, indels, lost stretch, strands)])
        "S1": ("SE", ["RUN1", "RUN2"], [
            (0.80, s1, [del3], None, "both"),
            (0.20, s1, [del3, ins], None, "forward"),
        ]),
        "S2": ("PE", ["RUN1"], [
            (1.00, s2, [ins], None, "both"),
        ]),
        "S3": ("PE", ["RUN1"], [
            (0.20, s3, [], LOST, "both"),
            (0.25, s3 + minority_snp + mnv_snps["S3"], [del1], LOST, "both"),
            (0.30, s3 + split_1, [], LOST, "both"),
            (0.25, s3 + split_2, [], LOST, "both"),
        ]),
    }

    def aa(name, k, codon):
        return f"{CODON_TABLE[cds[name][3 * (k - 1):3 * k]]}{k}{CODON_TABLE[codon]}"

    truth = {
        "contig": CONTIG, "reference": REF_ID, "genome_length": GENOME_LEN,
        "samples": {s: {"layout": v[0], "runs": v[1]} for s, v in samples.items()},
        "genes": [{"locus": g[0], "name": g[1], "start": g[2],
                   "end": g[2] + 3 * g[3] - 1, "strand": g[4]} for g in GENES],
        "repeat": [list(REPEAT), [REPEAT_COPY, REPEAT_COPY + REPEAT[1] - REPEAT[0]]],
        "lost": {"sample": "S3", "start": LOST[0], "end": LOST[1], "genes": ["geneE"]},
        "snps": [], "mnvs": [], "indels": [], "uncalled": [],
    }

    # One entry per site, with the fraction of each sample's molecules that carry it.
    sites = {}
    for entries, s, af in [(fixed["S1"], "S1", 1.0), (fixed["S2"], "S2", 1.0), (fixed["S3"], "S3", 1.0),
                           (mnv_snps["S1"], "S1", 1.0), (mnv_snps["S2"], "S2", 1.0),
                           (mnv_snps["S3"], "S3", 0.25), (minority_snp, "S3", 0.25),
                           (split_1, "S3", 0.30), (split_2, "S3", 0.25)]:
        for site in entries:
            sites.setdefault(site, {})[s] = af
    # Two changes of one codon on the same reads are one MNP record in FreeBayes' raw VCF, annotated
    # as the codon they make; the SNP matrix and the report read that VCF and keep the codon's change
    # on each base. The main VCF splits the record first, so each base carries its own change there.
    whole = {pos: f"p.{aa(name, k, codon)}" for s, (name, k, codon) in mnv.items() for pos, _, _ in mnv_snps[s]}
    for (pos, ref, alt), af in sorted(sites.items()):
        name, hgvs = protein_change(genome, pos, ref, alt)
        truth["snps"].append({"pos": pos, "ref": ref, "alt": alt, "af": af, "gene": name, "hgvs_p": hgvs,
                              "matrix_aa": whole.get(pos, hgvs)})

    for s, (name, k, codon) in mnv.items():
        truth["mnvs"].append({
            "sample": s, "gene": name, "codon": k, "positions": sorted(p for p, _, _ in mnv_snps[s]),
            "aa_change": aa(name, k, codon), "af": 1.0 if s != "S3" else 0.25,
            "masked": s == "S3",   # its first base alone would read a stop
        })
    truth["unphased_codon"] = {"sample": "S3", "gene": "geneA", "codon": 380,
                               "positions": sorted(p for p, _, _ in split_1 + split_2)}

    for label, (pos, ref, alt), effect, af, flt in [
        ("frameshift insertion", ins, "frameshift_variant", {"S2": 1.0, "S1": 0.20},
         {"S2": "PASS", "S1": "LowSupport"}),   # S1's copy is on forward reads only
        ("minority frameshift deletion", del1, "frameshift_variant", {"S3": 0.25}, {"S3": "PASS"}),
        ("in-frame deletion", del3, "conservative_inframe_deletion", {"S1": 1.0}, {"S1": "PASS"}),
    ]:
        name, hgvs = protein_change(genome, pos, ref, alt)
        truth["indels"].append({"name": label, "pos": pos, "ref": ref, "alt": alt, "gene": name,
                                "effect": effect, "hgvs_p": hgvs, "af": af, "filter": flt})
    truth["uncalled"] = [{"pos": REPEAT_SNP, "sample": "S2", "why": "in a repeat"}]
    return samples, truth


# ---------------------------------------------------------------------------------------------
# Reads
# ---------------------------------------------------------------------------------------------

def haplotype(genome, snps, indels, lost):
    """The haplotype's sequence. Edits are applied right to left so earlier positions hold."""
    seq = list(genome)
    edits = [(pos, "snp", (ref, alt)) for pos, ref, alt in snps]
    edits += [(pos, "indel", (ref, alt)) for pos, ref, alt in indels]
    if lost:
        edits.append((lost[0], "lost", lost))
    for pos, kind, data in sorted(edits, key=lambda e: -e[0]):
        if kind == "snp":
            ref, alt = data
            assert seq[pos - 1] == ref, (pos, ref, seq[pos - 1])
            seq[pos - 1] = alt
        elif kind == "indel":
            ref, alt = data
            assert "".join(seq[pos - 1:pos - 1 + len(ref)]) == ref, (pos, ref)
            seq[pos - 1:pos - 1 + len(ref)] = list(alt)
        else:
            seq[data[0] - 1:data[1]] = []
    return "".join(seq)


def with_errors(read, rng):
    out = list(read)
    for i in range(len(out)):
        if rng.random() < ERROR_RATE:
            out[i] = rng.choice([b for b in BASES if b != out[i]])
    return "".join(out)


def simulate(seq, coverage, layout, strands, rng):
    """Reads from one haplotype at the given coverage: pairs (r1, r2) or single reads."""
    if layout == "PE":
        n = int(len(seq) * coverage / (2 * READ_LEN))
        for _ in range(n):
            size = max(2 * READ_LEN, min(len(seq), int(rng.gauss(INSERT_MEAN, INSERT_SD))))
            start = rng.randint(0, len(seq) - size)
            frag = seq[start:start + size]
            if strands == "both" and rng.random() < 0.5:
                frag = revcomp(frag)
            yield with_errors(frag[:READ_LEN], rng), with_errors(revcomp(frag[-READ_LEN:]), rng)
    else:
        n = int(len(seq) * coverage / READ_LEN)
        for _ in range(n):
            start = rng.randint(0, len(seq) - READ_LEN)
            read = seq[start:start + READ_LEN]
            if strands == "both" and rng.random() < 0.5:
                read = revcomp(read)
            yield with_errors(read, rng), None


def write_fastq(path, names, reads):
    # mtime=0 keeps the gzip header constant, so the same seed gives the same bytes.
    payload = "".join(f"@{n}\n{r}\n+\n{QUAL}\n" for n, r in zip(names, reads))
    with path.open("wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", compresslevel=6, mtime=0) as fh:
        fh.write(payload.encode())


# ---------------------------------------------------------------------------------------------
# Files
# ---------------------------------------------------------------------------------------------

def write_fasta(path, seq):
    with path.open("w") as fh:
        fh.write(f">{CONTIG}\n")
        for i in range(0, len(seq), 80):
            fh.write(seq[i:i + 80] + "\n")


def write_gff(path):
    with path.open("w") as fh:
        fh.write(f"##gff-version 3\n##sequence-region {CONTIG} 1 {GENOME_LEN}\n")
        fh.write(f"{CONTIG}\te2e\tregion\t1\t{GENOME_LEN}\t.\t+\t.\tID={CONTIG};Name={CONTIG}\n")
        for locus, name, start, n_codons, strand in GENES:
            end = start + 3 * n_codons - 1
            fh.write(f"{CONTIG}\te2e\tgene\t{start}\t{end}\t.\t{strand}\t.\t"
                     f"ID=gene-{locus};Name={name};gene={name};locus_tag={locus}\n")
            fh.write(f"{CONTIG}\te2e\tCDS\t{start}\t{end}\t.\t{strand}\t0\t"
                     f"ID=cds-{locus};Parent=gene-{locus};Name={name};gene={name};locus_tag={locus};"
                     f"product={name} e2e product\n")


def build(outdir) -> dict:
    """Write the cohort under outdir and return the truth (also written as truth.json)."""
    out = Path(outdir).resolve()
    reads_dir = out / "reads"
    reads_dir.mkdir(parents=True, exist_ok=True)
    rng = random.Random(SEED)

    genome, cds = build_reference(rng)
    samples, truth = plan(genome, cds)
    write_fasta(out / "reference.fasta", genome)
    write_gff(out / "reference.gff3")

    rows = []
    for sid, (layout, runs, haps) in samples.items():
        per_run = COVERAGE / len(runs)
        for run in runs:
            reads = []
            for weight, snps, indels, lost, strands in haps:
                seq = haplotype(genome, snps, indels, lost)
                reads += list(simulate(seq, per_run * weight, layout, strands, rng))
            rng.shuffle(reads)
            names = [f"{sid}_{run}_{i:06d}" for i in range(len(reads))]
            if layout == "PE":
                r1 = reads_dir / f"{sid}_{run}_R1.fq.gz"
                r2 = reads_dir / f"{sid}_{run}_R2.fq.gz"
                write_fastq(r1, [f"{n}/1" for n in names], [a for a, _ in reads])
                write_fastq(r2, [f"{n}/2" for n in names], [b for _, b in reads])
                rows.append((sid, run, r1, r2))
            else:
                r1 = reads_dir / f"{sid}_{run}.fq.gz"
                write_fastq(r1, names, [a for a, _ in reads])
                rows.append((sid, run, r1, ""))

    with (out / "e2e.tsv").open("w") as fh:
        fh.write("sampleId\trunId\tr1\tr2\trefId\trefFasta\trefGff\ttaxId\n")
        for sid, run, r1, r2 in rows:
            fh.write(f"{sid}\t{run}\t{r1}\t{r2}\t{REF_ID}\t{out / 'reference.fasta'}\t"
                     f"{out / 'reference.gff3'}\t\n")

    truth["distances"] = expected_distances(truth)
    (out / "truth.json").write_text(json.dumps(truth, indent=1) + "\n")
    return truth


def expected_distances(truth):
    """SNPs between two consensus sequences: the fixed SNPs one sample carries and the other not.

    A minority never counts (the consensus writes an ambiguity code or the major base), and no
    fixed SNP lies in a repeat or in the stretch S3 lost, so every fixed difference is compared.
    """
    fixed = {s: {e["pos"] for e in truth["snps"] if e["af"].get(s) == 1.0} for s in truth["samples"]}
    names = sorted(truth["samples"])
    return {f"{a}|{b}": len(fixed[a] ^ fixed[b]) for i, a in enumerate(names) for b in names[i + 1:]}


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    t = build(sys.argv[1])
    print(f"{GENOME_LEN} bp, {len(GENES)} genes, {len(t['samples'])} samples, {len(t['snps'])} SNP sites, "
          f"{len(t['mnvs'])} codon changes, {len(t['indels'])} indels -> {sys.argv[1]}")
