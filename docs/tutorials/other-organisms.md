# Running a non-TB organism

BAMpiro maps, calls variants and builds its QC report for **any** bacterium - only its *defaults* are TB-tuned. This page shows you the handful of switches to flip so a non-tuberculosis organism runs correctly, and which TB-specific features to leave off. Budget about 10 minutes.

If you have not run the pipeline before, do [Your first run](first-run.md) first - the mechanics are identical; here we only change the reference, one profile and a couple of flags.

## The key idea

BAMpiro's built-in defaults assume *M. tuberculosis*: **diploid** calling (ploidy 2, to catch mixed MTBC infections) plus the MTBC-only extras (H37Rv lineage/DR typing, H37Rv coordinate liftover, blind-spot masking, H37Rv amino-acid numbering). For another organism you point the samplesheet at your genome and add one profile that flips those TB assumptions back to sane bacterial defaults.

## 1. Point the samplesheet at your reference

Everything organism-specific lives in the samplesheet columns - no code changes. In each row set:

- **`refFasta` / `refGff`** to *your* organism's genome and annotation (GFF3),
- **`taxId`** to your species' NCBI TaxID, which drives the Kraken2 contamination screen (e.g. `562` for *E. coli*).

See [The samplesheet](samplesheet.md) for the full column reference.

!!! warning "The GFF seqid must match the FASTA header"

    `refGff`'s first column (seqid) must be identical to the `refFasta` header id (up to the first whitespace), or the per-sample SnpEff database build fails or annotates nothing. If you have no annotation, pass a minimal GFF3 stub - see [Reference requirements](../quickstart.md#reference-requirements).

## 2. Add the `generic` profile

`generic` is a built-in profile that flips the TB defaults. Combine it with your executor profile (they compose with commas):

=== "SLURM cluster"

    ```bash
    nextflow run main.nf --tsv samples.tsv --outdir results -profile slurm,generic
    ```

=== "Laptop / workstation (Docker)"

    ```bash
    nextflow run main.nf --tsv samples.tsv --outdir results -profile local,docker,generic
    ```

What `generic` changes versus the TB defaults:

| Setting | TB default | `generic` |
| :--- | :--- | :--- |
| `freebayes_ploidy` | `2` (diploid - detects mixed infections) | **`1`** (a clonal bacterium is haploid) |
| `run_pathotypr` · `variant_liftover` · `mask_blindspots` · `annotate_canonical` | optional / off | **forced off** (all MTBC-specific) |

Everything else is already generic. For the full table see [Configuration](../configuration.md).

## 3. Leave the TB typing off

The `generic` profile already handles this, but it is worth understanding:

- **Lineage / DR typing** (`--run_pathotypr`) is MTBC-specific and **off by default** - leave it off. See [Lineage & Drug-Resistance Typing](../pathotypr.md).
- **`--annotate_canonical`** stays **off**: H37Rv amino-acid numbering is meaningless for a non-MTBC organism. If you *do* want a canonical numbering for your species, set `--annotate_canonical true` with `--canonical_snpeff_db` (and `--canonical_label`) pointing at another SnpEff genome.

!!! note "Only the H37Rv SnpEff DB is bundled"

    The image ships just the `Mycobacterium_tuberculosis_h37rv` database. To use `--canonical_snpeff_db` for another organism, either add its genome to the image and rebuild, or ensure a bare database id can be downloaded by SnpEff at run time (needs internet).

## 4. Tune the QC thresholds

The pass/fail gate defaults (depth, breadth, missing, and so on) were picked for a ~4.4 Mb TB genome at typical MTBC depths. A different genome size or expected coverage may warrant different cut-offs - adjust the `--report_*` parameters. See [Configuring a run](configuring-a-run.md) for how the gate works and which knobs to turn.

## 5. The QC report is organism-agnostic

The interactive report still builds in full. The MTBC-only panels (per-lineage summary, drug resistance, dual amino-acid numbering) simply **self-hide** when there is no such data - everything else (depth/breadth, SNP matrix, genome landscape, gene burden, SNP dynamics…) behaves exactly as for TB. See the [Interactive QC Report](../qc-report.md).

## Worked example - *E. coli* K-12

A haploid *E. coli* run on a laptop with Docker, screened against TaxID `562`:

```bash
nextflow run main.nf \
    --tsv samples.tsv \
    --outdir results_ecoli \
    -profile local,docker,generic \
    --kraken2_db /data/kraken2
```

with a samplesheet row such as:

```tsv
sampleId	runId	r1	r2	refId	refFasta	refGff	taxId
ECOLI_1	RUN1	/data/EC_R1.fq.gz	/data/EC_R2.fq.gz	K12	/refs/ecoli_K12.fa	/refs/ecoli_K12.gff	562
```

That is a complete, correct non-TB run: haploid calling, the TB-only features off, and Kraken2 screening for *E. coli*.

## Where to go next

- Fine-tune profiles, parameters and the QC gate in [Configuring a run](configuring-a-run.md).
- Return to the [Tutorials overview](index.md) for the rest of the guided path.
