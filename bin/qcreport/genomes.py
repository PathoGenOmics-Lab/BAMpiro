"""Each reference's own genome for the report: its length, genes and masked regions, on the axis its samples'
profiles are binned on.

A cohort can be mapped against several references. Every sample's landscape profile is binned over its own
reference, contigs end to end in the order of the reference's FASTA index (stats_to_legacy.py), so the genes
and masked regions drawn under it have to be that reference's and placed the same way: a gene of the second
contig sits after the whole of the first one. One genome per reference, each shaped like the single one the
report has always carried (len, genes, mask_bins, mask_pct, mask_iv), so the page can switch between them.
"""
from __future__ import annotations

import gzip
import os

from .parsers import NBINS, gff_features, mask_profile

MASK_IV_CAP = 5000


def _open(path):
    return gzip.open(path, "rt", encoding="utf-8", errors="replace") if str(path).endswith(".gz") \
        else open(path, encoding="utf-8", errors="replace")


def _usable(path):
    return bool(path) and os.path.exists(path) and not os.path.basename(str(path)).startswith("NO_FILE")


def read_fai(path):
    """[(contig, length)] of a FASTA index, in its order; [] without one."""
    out = []
    if not _usable(path):
        return out
    try:
        with open(path, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                c = line.rstrip("\n").split("\t")
                if len(c) >= 2 and c[1].isdigit():
                    out.append((c[0], int(c[1])))
    except OSError:
        return []
    return out


def bed_intervals(path):
    """[(contig, start, end)] of a BED or of a reference's exclusion list; [] without one."""
    out = []
    if not _usable(path):
        return out
    try:
        with _open(path) as fh:
            for line in fh:
                if not line.strip() or line.startswith(("#", "track", "browser")):
                    continue
                c = line.rstrip("\n").split("\t")
                if len(c) >= 3 and c[1].strip().isdigit() and c[2].strip().isdigit():
                    out.append((c[0], int(c[1]), int(c[2])))
    except OSError:
        return []
    return out


def _placer(contigs):
    """A function placing (contig, position) on the genome-wide axis, or None where it cannot. With a single
    contig every name is it: a GFF may call the chromosome otherwise than the FASTA does."""
    offsets, run = {}, 0
    for name, length in contigs:
        offsets[name] = run
        run += length
    single = len(contigs) <= 1

    def place(contig, pos):
        if contig in offsets:
            return offsets[contig] + pos
        return pos if single else None
    return place


def build_genome(gff=None, mask=None, fai=None, length=0, nbins=NBINS):
    """{len, genes, mask_bins, mask_pct, mask_iv} of one reference, positions on its genome-wide axis."""
    contigs = read_fai(fai)
    place = _placer(contigs)
    genes = []
    for contig, name, start, end in gff_features(gff):
        s, e = place(contig, start), place(contig, end)
        if s is not None and e is not None:
            genes.append({"name": name, "start": s, "end": e})
    genes.sort(key=lambda g: (g["start"], g["end"]))
    glen = sum(n for _, n in contigs) or length or max((g["end"] for g in genes), default=0)
    iv = []
    for contig, start, end in bed_intervals(mask):
        s, e = place(contig, start), place(contig, end)
        if s is not None and e is not None:
            iv.append((s, e))
    bins, pct = mask_profile(iv, glen, nbins)
    return {"len": glen, "genes": genes, "mask_bins": bins, "mask_pct": pct,
            "mask_iv": (iv[:MASK_IV_CAP] if iv else None)}


def build_genomes(ref_gff, ref_mask, ref_fai, lengths=None, nbins=NBINS):
    """{reference: genome} for every reference named in any of the three {reference: path} maps. `lengths`
    ({reference: bp}, from the samples' own summary) stands in for a reference without a FASTA index."""
    lengths = lengths or {}
    refs = sorted(set(ref_gff) | set(ref_mask) | set(ref_fai))
    return {r: build_genome(ref_gff.get(r), ref_mask.get(r), ref_fai.get(r), lengths.get(r, 0), nbins)
            for r in refs}


def pairs(items):
    """{reference: path} from REF=PATH arguments."""
    out = {}
    for it in items or []:
        ref, sep, path = str(it).partition("=")
        if sep and ref and path:
            out[ref] = path
    return out
