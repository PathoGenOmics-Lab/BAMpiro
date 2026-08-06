#!/usr/bin/env python3
"""Regenerate the tiny fixture cohort under tests/data/.

Everything is derived from a fixed seed, so re-running this script reproduces the
committed files byte for byte. The fixtures are deliberately small (a 5 kb genome,
three samples) so they can live in the repository and still exercise the pipeline's
paired-end, single-end and multi-run code paths.

    python3 tests/data/make_test_data.py

The generated cohort is what `-profile test` and the pytest suite read.
"""

from __future__ import annotations

import gzip
import random
from pathlib import Path

HERE = Path(__file__).resolve().parent

SEED = 20260806
CONTIG = "test_chr"
TAXID = 1773
GENOME_LEN = 5000
READ_LEN = 100
COVERAGE = 30
INSERT = 300
QUAL = "I" * READ_LEN

# Genes are laid out so every one of them is hit by at least one planted variant.
GENES = [
    # (locus_tag, gene_name, start, end, strand)
    ("TEST_0001", "dnaA", 201, 1400, "+"),
    ("TEST_0002", "rpoB", 1601, 3100, "+"),
    ("TEST_0003", "katG", 3301, 4700, "-"),
]

# Planted variants: (position 1-based, alt base, which samples carry it).
VARIANTS = [
    (450, "T", {"SAMPLE_A", "SAMPLE_B"}),
    (1750, "G", {"SAMPLE_B"}),
    (2100, "A", {"SAMPLE_A", "SAMPLE_C"}),
    (3500, "C", {"SAMPLE_C"}),
    (4200, "T", {"SAMPLE_A", "SAMPLE_B", "SAMPLE_C"}),
]

# (sampleId, runId, layout) - SAMPLE_A has two runs so the merge path is exercised.
SAMPLES = [
    ("SAMPLE_A", "RUN1", "PE"),
    ("SAMPLE_A", "RUN2", "PE"),
    ("SAMPLE_B", "RUN1", "PE"),
    ("SAMPLE_C", "RUN1", "SE"),
]

COMPLEMENT = str.maketrans("ACGT", "TGCA")


def revcomp(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def build_genome(rng: random.Random) -> str:
    return "".join(rng.choice("ACGT") for _ in range(GENOME_LEN))


def sample_genome(reference: str, sample_id: str) -> str:
    bases = list(reference)
    for pos, alt, carriers in VARIANTS:
        if sample_id in carriers:
            bases[pos - 1] = alt
    return "".join(bases)


def write_fasta(path: Path, header: str, seq: str, wrap: int = 60) -> None:
    with path.open("w") as fh:
        fh.write(f">{header}\n")
        for i in range(0, len(seq), wrap):
            fh.write(seq[i : i + wrap] + "\n")


def write_gff(path: Path) -> None:
    with path.open("w") as fh:
        fh.write("##gff-version 3\n")
        fh.write(f"##sequence-region {CONTIG} 1 {GENOME_LEN}\n")
        fh.write(
            f"{CONTIG}\ttest\tregion\t1\t{GENOME_LEN}\t.\t+\t.\t"
            f"ID={CONTIG};Name={CONTIG}\n"
        )
        for locus, name, start, end, strand in GENES:
            fh.write(
                f"{CONTIG}\ttest\tgene\t{start}\t{end}\t.\t{strand}\t.\t"
                f"ID=gene-{locus};Name={name};locus_tag={locus}\n"
            )
            fh.write(
                f"{CONTIG}\ttest\tCDS\t{start}\t{end}\t.\t{strand}\t0\t"
                f"ID=cds-{locus};Parent=gene-{locus};Name={name};locus_tag={locus};"
                f"product={name} test product\n"
            )


def paired_records(seq: str, sample_id: str, run_id: str, rng: random.Random):
    n_pairs = (GENOME_LEN * COVERAGE) // (READ_LEN * 2)
    for i in range(n_pairs):
        start = rng.randint(0, len(seq) - INSERT)
        fragment = seq[start : start + INSERT]
        r1 = fragment[:READ_LEN]
        r2 = revcomp(fragment[-READ_LEN:])
        name = f"{sample_id}_{run_id}_{i:05d}"
        yield name, r1, r2


def single_records(seq: str, sample_id: str, run_id: str, rng: random.Random):
    n_reads = (GENOME_LEN * COVERAGE) // READ_LEN
    for i in range(n_reads):
        start = rng.randint(0, len(seq) - READ_LEN)
        read = seq[start : start + READ_LEN]
        if rng.random() < 0.5:
            read = revcomp(read)
        yield f"{sample_id}_{run_id}_{i:05d}", read


def write_fastq(path: Path, records) -> None:
    # mtime=0 keeps the gzip header constant so regenerating gives an identical file.
    payload = "".join(f"@{name}\n{seq}\n+\n{QUAL}\n" for name, seq in records)
    with path.open("wb") as raw:
        with gzip.GzipFile(fileobj=raw, mode="wb", compresslevel=9, mtime=0) as fh:
            fh.write(payload.encode())


def main() -> None:
    rng = random.Random(SEED)
    reads_dir = HERE / "reads"
    reads_dir.mkdir(parents=True, exist_ok=True)

    reference = build_genome(rng)
    write_fasta(HERE / "reference.fasta", CONTIG, reference)
    write_gff(HERE / "reference.gff3")

    rows = []
    for sample_id, run_id, layout in SAMPLES:
        seq = sample_genome(reference, sample_id)
        if layout == "PE":
            pairs = list(paired_records(seq, sample_id, run_id, rng))
            r1 = reads_dir / f"{sample_id}_{run_id}_R1.fq.gz"
            r2 = reads_dir / f"{sample_id}_{run_id}_R2.fq.gz"
            write_fastq(r1, ((f"{n}/1", s) for n, s, _ in pairs))
            write_fastq(r2, ((f"{n}/2", s) for n, _, s in pairs))
            rows.append((sample_id, run_id, r1.name, r2.name))
        else:
            reads = list(single_records(seq, sample_id, run_id, rng))
            r1 = reads_dir / f"{sample_id}_{run_id}.fq.gz"
            write_fastq(r1, reads)
            rows.append((sample_id, run_id, r1.name, ""))

    header = "sampleId\trunId\tr1\tr2\trefId\trefFasta\trefGff\ttaxId\n"
    with (HERE / "samplesheet.tsv").open("w") as fh:
        fh.write(header)
        for sample_id, run_id, r1, r2 in rows:
            r1p = f"tests/data/reads/{r1}"
            r2p = f"tests/data/reads/{r2}" if r2 else ""
            fh.write(
                f"{sample_id}\t{run_id}\t{r1p}\t{r2p}\tTESTREF\t"
                f"tests/data/reference.fasta\ttests/data/reference.gff3\t{TAXID}\n"
            )

    write_invalid_samplesheets(rows)
    write_pathotypr_panels(reference)

    # -profile test_full points --kraken2_db here; main.nf only checks that it exists.
    kraken_dir = HERE / "kraken2_db"
    kraken_dir.mkdir(exist_ok=True)
    (kraken_dir / ".gitkeep").write_text(
        "Placeholder so git tracks this directory. -profile test_full points\n"
        "--kraken2_db here purely so the Kraken branch of the DAG is exercised;\n"
        "a stub run never opens the database.\n"
    )

    print(f"reference : {GENOME_LEN} bp, {len(GENES)} genes")
    print(f"samples   : {len({s for s, _, _ in SAMPLES})} in {len(SAMPLES)} runs")
    print(f"variants  : {len(VARIANTS)} planted")


HEADER = ["sampleId", "runId", "r1", "r2", "refId", "refFasta", "refGff", "taxId"]
REF_FA = "tests/data/reference.fasta"
REF_GFF = "tests/data/reference.gff3"


def write_invalid_samplesheets(rows) -> None:
    """Samplesheets that main.nf's validator must reject, one per failure mode.

    tests/pipeline/test_samplesheet_validation.py asserts on the exact message each of these
    produces, so the shape of every row matters. In particular the conflicting-reference case
    needs a VALID row first: the validator only notices a conflict once a refId has been bound,
    and rows that fail an earlier check return before binding anything.
    """
    out = HERE / "invalid"
    out.mkdir(exist_ok=True)
    good = rows[0]
    ok_row = [good[0], good[1], f"tests/data/reads/{good[2]}",
              f"tests/data/reads/{good[3]}", "TESTREF", REF_FA, REF_GFF, str(TAXID)]
    # An existing file that is not the reference, for the conflicting-refFasta case.
    other_file = f"tests/data/reads/{rows[-1][2]}"

    def write(name, header, data_rows):
        with (out / name).open("w") as fh:
            fh.write("\t".join(header) + "\n")
            for row in data_rows:
                fh.write("\t".join(row) + "\n")

    write("missing_column.tsv",
          [c for c in HEADER if c != "refGff"],
          [[v for i, v in enumerate(ok_row) if HEADER[i] != "refGff"]])

    write("no_rows.tsv", HEADER, [])

    with (out / "no_usable_rows.tsv").open("w") as fh:
        fh.write("\t".join(HEADER) + "\n")
        fh.write("#" + "\t".join(ok_row) + "\n")
        fh.write("\n")

    write("empty_field.tsv", HEADER, [[""] + ok_row[1:]])

    write("missing_file.tsv", HEADER,
          [[ok_row[0], ok_row[1], "tests/data/reads/DOES_NOT_EXIST_R1.fq.gz"] + ok_row[3:]])

    write("ref_conflict.tsv", HEADER,
          [ok_row, ["SAMPLE_B", "RUN1", other_file, "", "TESTREF", other_file, REF_GFF, str(TAXID)]])

    # Three DIFFERENT failures at once: the validator must report all of them in one pass.
    write("multi_error.tsv", HEADER, [
        ok_row,
        ["", "RUN1", other_file, "", "TESTREF", REF_FA, REF_GFF, str(TAXID)],
        ["SAMPLE_B", "RUN1", "tests/data/reads/NOPE_R1.fq.gz", "", "TESTREF", REF_FA, REF_GFF, str(TAXID)],
        ["SAMPLE_C", "RUN1", other_file, "", "TESTREF", other_file, REF_GFF, str(TAXID)],
    ])


def write_pathotypr_panels(reference: str) -> None:
    """Marker panels shaped like the real Zenodo v1.0.0 files.

    -profile test_full stages these so the pathotypr branch of the DAG runs; a stub
    run never parses them, but the column layout matches the real panels so they are
    also usable as fixtures for the parsing tests.
    """
    out = HERE / "pathotypr"
    out.mkdir(exist_ok=True)
    write_fasta(out / "reference.fasta", CONTIG, reference)

    with (out / "lineage_markers.tsv").open("w") as fh:
        fh.write("position\tref_base\talt_base\tlineage1\tlineage2\n")
        for pos, alt, _ in VARIANTS[:3]:
            fh.write(f"{pos}\t{reference[pos - 1]}\t{alt}\tL1\tL1.1\n")

    with (out / "dr_markers.tsv").open("w") as fh:
        fh.write("#pos\tref\talt\tdrug\tresistance\tmarker_name\tgrade\tgene\tmutation\n")
        fh.write(
            f"1750\t{reference[1749]}\tG\trifampicin\tR\trpoB_S450L\t"
            "1) Assoc w R\trpoB\tp.Ser450Leu\n"
        )
        fh.write(
            f"3500\t{reference[3499]}\tC\tisoniazid\tR\tkatG_S315T\t"
            "1) Assoc w R\tkatG\tp.Ser315Thr\n"
        )


if __name__ == "__main__":
    main()
