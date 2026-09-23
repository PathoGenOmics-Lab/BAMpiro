#!/usr/bin/env python3
"""Consolidated interactive QC report for BAMpiro (organism-agnostic).

Reads the per-sample summary TSV (collect_summary output) + the consensus FASTAs, computes consensus completeness
(missing / N / IUPAC / callable), embeds everything as JSON, and writes a SELF-CONTAINED, data-driven qc_report.html
(one inline vanilla-JS module, no CDN) plus a qc_flags.tsv (PASS/WARN/FAIL). Features:
  - General Statistics table: value-coloured cells, sort, sample filter, column show/hide, CSV export, sticky header
    + frozen sample column.
  - LIVE thresholds: adjust the flag cut-offs in the report and everything re-flags instantly.
  - Distributions: beeswarm / bar / histogram (switchable) per metric, with hover tooltips.
  - Correlations: a scatter with X/Y metric pickers, coloured by verdict.
  - Consensus completeness: stacked callable / IUPAC / missing bar per sample.
  - Click any row, dot or point to highlight a sample everywhere; click a flag chip to filter.
With --gate it exits non-zero if any sample FAILs (at the default thresholds), so a downstream target can block.
Pure standard library (no numpy / matplotlib / CDN). Study-specific signatures live in a separate downstream gate.

The implementation lives in the qcreport package next to this file: qcreport.parsers
reads the pipeline's artefacts, qcreport.metrics holds the metric registry and the
verdict engine, qcreport.panels builds the optional analysis panels and qcreport.render
assembles the HTML. This file is the CLI: arguments in, payload assembled, two files out.
"""
from __future__ import annotations

import argparse
import os
import sys
from datetime import datetime, timezone

from qcreport.coverage import build_coverage, parse_depth_profiles, snp_per_kb, write_deletions
from qcreport.genomes import build_genomes, pairs
from qcreport.metrics import (ANC_DEF, DEF, DEFS, DIST, METRICS, build_gene_map, discover_extra_metrics,
                              flag_sample, het_frac, is_ancient, lineage_counts_parsed, lineage_fracs, robust)
from qcreport.panels import build_dynamics, build_epistasis, build_snp_matrix
from qcreport.parsers import (NBINS, clean_str, consensus_stats, mapdamage_stats, mask_profile, parse_bed,
                              parse_collection_dates, parse_dose, parse_dr, parse_gene_burden,
                              parse_gene_conversion, parse_gff,
                              parse_kraken, parse_lineage_colors, parse_metadata, parse_pnps, parse_profile,
                              parse_sample_meta, parse_summary, parse_vcfs, to_float)
from qcreport.relatedness import build_relatedness, group_mismatches, parse_pairs
from qcreport.minority import build_minority, replicate_pairs, replicate_sets
from qcreport.minority import needed_cells as replicate_cells
from qcreport.series import build_series, needed_cells, read_matrix_cells
from qcreport.render import REPO_URL, SECTION_INFO, build_html


def build_parser():
    """The CLI. Every threshold in DEF / ANC_DEF gets a --flag automatically, so adding a
    threshold to the metrics module is all it takes to make it settable from the command line."""
    ap = argparse.ArgumentParser(description="BAMpiro interactive QC report + flags.")
    ap.add_argument("--summary", required=True)
    ap.add_argument("--consensus", nargs="*", default=[])
    ap.add_argument("--out-html", required=True)
    ap.add_argument("--out-flags", required=True)
    ap.add_argument("--title", default="BAMpiro QC report")
    for k, v in DEF.items():
        ap.add_argument("--" + k.replace("_", "-"), type=float, default=v)
    for k, v in ANC_DEF.items():
        ap.add_argument("--anc-" + k.replace("_", "-"), type=float, default=v)
    ap.add_argument("--mapdamage-dir", nargs="*", default=[],
                    help="mapDamage2 output dir(s) for ancient samples (per-sample subdir or a parent to walk).")
    ap.add_argument("--gff", default=None,
                    help="GFF3 of gene coordinates -> per-gene SNP-density hotspot detection (optional).")
    ap.add_argument("--lineage-colors", default=None,
                    help="TSV 'lineage<TAB>#hex' (e.g. mycolorsTB) to colour lineages with a canonical palette (optional).")
    ap.add_argument("--mask-bed", default=None,
                    help="BED of regions to mask (e.g. mtbc_mask.bed: PE/PPE, IS, DR, repeats) -> a 'mask regions' toggle.")
    ap.add_argument("--ref-gff", nargs="*", default=[],
                    help="REF=GFF3 for each reference of the run -> its own genes under its own samples (genome "
                         "landscape, gene hotspots, go-to-gene). With these, --gff and --mask-bed are not used.")
    ap.add_argument("--ref-mask", nargs="*", default=[],
                    help="REF=BED for each reference: its repeat/exclude regions (the masked-regions toggle).")
    ap.add_argument("--ref-fai", nargs="*", default=[],
                    help="REF=FASTA index for each reference: its contigs' order and lengths, the axis its samples' "
                         "profiles are binned on, where genes and masked regions are placed.")
    ap.add_argument("--gene-burden", default=None,
                    help="Cohort gene-burden TSV (collect_summary --gene-burden-out) -> the Functional gene burden panel (optional).")
    ap.add_argument("--dr-report", default=None,
                    help="Per-sample drug-resistance calls TSV (from pathotypr DR markers) -> the Drug resistance panel (optional).")
    ap.add_argument("--gene-conversion", default=None,
                    help="Cohort gene-conversion tracts TSV (COLLECT_GENE_CONVERSION) -> the Gene conversion panel "
                         "(optional).")
    ap.add_argument("--kraken", nargs="*", default=[],
                    help="Per-sample Kraken2 .report files -> the Taxonomic composition / contamination panel (optional).")
    ap.add_argument("--pnps", default=None,
                    help="Cohort per-gene dN/dS TSV (eskaks) -> the Selection pN/pS panel (optional; a cohort analysis, not per-sample QC).")
    ap.add_argument("--provenance", nargs="*", default=[],
                    help="key=value provenance pairs surfaced in the report header (reference, container, commit...).")
    ap.add_argument("--version", default="",
                    help="BAMpiro version string (workflow.manifest.version) shown in the header + footer.")
    ap.add_argument("--metadata", default=None,
                    help="Optional TSV (e.g. the samplesheet) with a sample column + time (passage/timepoint) "
                         "and group (patient/series/cluster) columns -> the SNP dynamics panel. Auto-detected.")
    ap.add_argument("--vcfs", nargs="*", default=[],
                    help="Optional per-sample annotated VCFs -> per-SNP allele frequencies for the dynamics panel.")
    ap.add_argument("--vcfs-h37rv", nargs="*", default=[],
                    help="Optional per-sample VCFs annotated against a canonical reference -> the canonical "
                         "amino-acid change and gene shown alongside the used-reference ones. Records lifted to the "
                         "canonical coordinates carry the mapping coordinate in INFO/OPOS and are paired by it.")
    ap.add_argument("--aa2-label", default="H37Rv",
                    help="Label for the canonical-reference amino-acid numbering shown by --vcfs-h37rv (default H37Rv).")
    ap.add_argument("--pos-liftover", default=None,
                    help="TSV 'src_contig src_pos tgt_contig tgt_pos strand' (pathotypr_liftover.py lift "
                         "--global-chain --out-map, one per reference joined) -> the canonical COORDINATE of each "
                         "variant, looked up by contig and position. Fills pos_h37rv for the SNP tables.")
    ap.add_argument("--depth-profiles", nargs="*", default=[],
                    help="Per-sample depth tables from depth_profile.py (windows, stretches without reads, "
                         "genes) -> deletions and SNPs per callable kb along the genome (optional).")
    ap.add_argument("--deletion-min-len", type=int, default=200,
                    help="Shortest stretch without reads reported as a deletion, in bp.")
    ap.add_argument("--deletion-min-depth", type=float, default=10.0,
                    help="Median depth a sample needs before its stretches without reads are assessed: "
                         "in a thinly read sample they turn up by chance.")
    ap.add_argument("--out-deletions", default=None,
                    help="Write the deletion regions (with --depth-profiles) as a TSV.")
    ap.add_argument("--snp-distances", default=None,
                    help="Pairwise SNP distances between consensus sequences (snp_distances.py --pairs) -> "
                         "the relatedness page and the GROUP_MISMATCH flag (optional).")
    ap.add_argument("--cluster-snps", type=int, default=12,
                    help="SNPs within which two samples are drawn as one cluster, and beyond which a "
                         "sample is far from its group. The report can move it live.")
    ap.add_argument("--snp-matrix", default=None,
                    help="The master SNP matrix TSV (build_snp_matrix.py with --depth-vcfs) -> what each "
                         "series gained since its first time point, telling a site absent there from one "
                         "not read (optional).")
    ap.add_argument("--min-dp", type=int, default=7,
                    help="Depth a site needs before a sample without a call there counts as lacking the "
                         "allele (--consensus_min_dp).")
    ap.add_argument("--gate", action="store_true")
    return ap


def collect_damage(mapdamage_dirs, summ):
    """{sample: mapDamage stats} for the ancient samples, from the --mapdamage-dir path(s).
    A dir may be the sample's own mapDamage output, a parent of per-sample subdirs, or the
    published {OUTDIR}/{sample_dir}/mapDamage layout; all three are probed."""
    # ---- aDNA damage authentication: map each mapDamage dir to a sample id ----
    ancient_ids = {sid for sid, m in summ.items() if is_ancient(m)}
    dmg_by_sample = {}

    def _index_dir(d):
        if not os.path.isdir(d):
            return
        base = os.path.basename(os.path.normpath(d))
        if base in ancient_ids and base not in dmg_by_sample:
            st = mapdamage_stats(d)
            if st:
                dmg_by_sample[base] = st
        for aid in ancient_ids:
            sub = os.path.join(d, aid)
            if aid not in dmg_by_sample and os.path.isdir(sub):
                st = mapdamage_stats(sub)
                if st:
                    dmg_by_sample[aid] = st
        if any(os.path.exists(os.path.join(d, f)) for f in ("5pCtoT_freq.txt", "3pGtoA_freq.txt")):
            comps = os.path.normpath(d).split(os.sep)   # match a whole path COMPONENT, not a loose substring
            for aid in ancient_ids:
                if aid not in dmg_by_sample and aid in comps:
                    st = mapdamage_stats(d)
                    if st:
                        dmg_by_sample[aid] = st

    for d in (mapdamage_dirs or []):
        _index_dir(d)
        if os.path.isdir(d):
            for name in sorted(os.listdir(d)):
                p = os.path.join(d, name)
                _index_dir(p)
                sub = os.path.join(p, "mapDamage")   # published layout: {OUTDIR}/{sample_dir}/mapDamage
                if os.path.isdir(sub):
                    _index_dir(sub)
    if ancient_ids and not dmg_by_sample:
        sys.stderr.write(f"[qc_report] {len(ancient_ids)} ancient sample(s) but no mapDamage output matched; "
                         f"DAMAGE authentication skipped.\n")
    return dmg_by_sample


def read_lift_map(path):
    """{(contig, position): canonical position, or None where the canonical genome has no counterpart} from
    the liftover map, or {} without one.

    LIFT_VARIANTS writes src_contig src_pos tgt_contig tgt_pos strand, one map per reference joined into
    one, and '.' as the target of a position inserted relative to the canonical genome. A map in the older
    two-column form (src_pos tgt_pos) names no contig and is keyed on (None, pos): it can only be applied
    where the cohort has a single contig, since the same position of another contig, or of another reference,
    is another site."""
    out = {}
    if not path or not os.path.exists(path) or os.path.basename(path).startswith("NO_FILE"):
        return out
    try:
        with open(path, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                c = line.rstrip("\n").split("\t")
                if len(c) >= 4 and c[1].strip().isdigit() and c[3].strip().isdigit():
                    out[(c[0].strip(), int(c[1]))] = int(c[3])
                elif len(c) >= 4 and c[1].strip().isdigit() and c[3].strip() == ".":
                    out[(c[0].strip(), int(c[1]))] = None
                elif len(c) >= 2 and c[0].strip().isdigit() and c[1].strip().isdigit():
                    out[(None, int(c[0]))] = int(c[1])
    except OSError as e:
        sys.stderr.write("[qc_report] WARN pos-liftover (%s): %s\n" % (path, e))
    return out


def load_variants(args):
    """{sample: {'chrom:pos': variant}} from --vcfs, with the canonical-reference (H37Rv) amino acid, gene
    and coordinate attached where they are known.

    A canonical record is paired with a variant by the mapping coordinate it carries (OPOS, or Picard's
    OriginalContig/OriginalStart), which ANNOTATE_CANONICAL writes as it lifts each record: its position is
    then the canonical coordinate. A canonical record without one sits at the mapping coordinate itself (a
    VCF annotated in place, right only for a reference that shares the canonical coordinates): it lends its
    amino acid to the variant at the same contig and position, but not a coordinate, which it never had.
    The liftover map fills the coordinate of every variant still without one, looked up by contig."""
    _variants = parse_vcfs(args.vcfs)   # parsed once, feeds both the dynamics panel and the SNP matrix
    label = args.aa2_label
    if args.vcfs_h37rv:
        _h37 = parse_vcfs(args.vcfs_h37rv)
        for _s, _pm in _variants.items():
            _by_opos, _in_place = {}, {}
            for _ck, _cv in _h37.get(_s, {}).items():
                if _cv.get('opos'):
                    _by_opos[_cv['opos']] = (_ck, _cv)
                else:
                    _in_place[_ck] = _cv
            for _key, _v in _pm.items():
                _hit = _by_opos.get(_key)
                _cv = _hit[1] if _hit else _in_place.get(_key)
                if not _cv:
                    continue
                if _cv.get('aa'):
                    _v['aa_h37rv'] = _cv['aa']
                if _cv.get('gene'):
                    _v['gene_h37rv'] = _cv['gene']
                if _hit:
                    _v['pos_h37rv'] = "%s:%s" % (label, _hit[0].rpartition(':')[2])
    _lift = read_lift_map(args.pos_liftover)
    if _lift:
        contigs = {k.rpartition(':')[0] for _pm in _variants.values() for k in _pm}
        legacy = len(contigs) <= 1      # a map without contigs is only unambiguous on a single contig
        for _pm in _variants.values():
            for _key, _v in _pm.items():
                if _v.get("pos_h37rv"):
                    continue
                _c, _, _p = _key.rpartition(":")
                if not _p.isdigit():
                    continue
                _k = (_c, int(_p))
                if _k not in _lift and legacy:
                    _k = (None, int(_p))
                if _k in _lift:
                    # None: inserted relative to the canonical genome, which has no coordinate to give
                    _t = _lift[_k]
                    _v["pos_h37rv"] = "%s:%s" % (label, "absent" if _t is None else _t)
    return _variants


def build_payload(args, thr, anc_thr):
    """Everything the report embeds: (payload, jsamples, counts). jsamples and counts are
    returned alongside because the flags TSV and the --gate exit are built from them."""
    summ = parse_summary(args.summary)
    miss_by_sample = {}    # per-sample binned missingness profile along the reference (genome landscape)
    for path in args.consensus:
        try:
            cs = consensus_stats(path)
        except Exception as e:
            sys.stderr.write(f"[qc_report] WARN consensus {path}: {e}\n")
            continue
        sid = cs.get("consensus_id")
        if sid not in summ:
            stem = os.path.basename(path).split(".")[0]
            sid = stem if stem in summ else sid
        keys = ("length", "missing_pct", "iupac_pct", "callable_pct", "longest_n_run", "n_gaps")
        if sid in summ:
            summ[sid].update({k: cs[k] for k in keys})
        elif sid:
            summ[sid] = {"sample_id": sid, **{k: cs[k] for k in keys}}
        if sid:
            miss_by_sample[sid] = cs.get("miss_profile")

    # SNP-count robust-z is computed over MODERN samples only (ancient counts are low by design and would skew it)
    snp_med, snp_sig = robust([to_float(m.get("snps")) for m in summ.values() if not is_ancient(m)])
    # representative genome length: consensus length, else the reference length column (works when make_consensus is off)
    genome_len = max([int(v) for v in (to_float(m.get("length")) or to_float(m.get("genome_len"))
                                       for m in summ.values()) if v], default=0)
    snp_density_ok = genome_len > 100000

    dmg_by_sample = collect_damage(args.mapdamage_dir, summ)

    metric_keys = [k for k, _, _, _ in METRICS]
    # fine-grained snpEff effect-class columns: carried into s.m for the Functional-annotation panel, not registered metrics
    EFF_KEYS = ["eff_missense", "eff_synonymous", "eff_stop_gained", "eff_stop_lost", "eff_start_lost",
                "eff_frameshift", "eff_inframe_indel", "eff_splice", "eff_intergenic", "eff_regulatory"]
    extra_metrics = discover_extra_metrics(summ)
    # a numeric 'dose' samplesheet column becomes a first-class (hidden-by-default) metric so it
    # can go on the scatter axes + correlation matrix and drive the dose x treatment test.
    dose_map = parse_dose(args.metadata)
    if dose_map:
        extra_metrics = extra_metrics + [{"key": "dose", "label": "Dose", "kind": "float", "dir": "neu"}]
    extra_keys = [e["key"] for e in extra_metrics]
    # The lineage the samples sharing a reference agree on. A reference is chosen per sample, so a
    # sample whose reads type as something else was routed to the wrong one; the majority is what
    # "should have been here" without needing to know the reference's own lineage.
    # The reference comes from the SAMPLESHEET, because the summary does not carry one. It used to
    # be read from summ[sid]["reference"], a key that never exists, so the majority below was
    # always empty and LINEAGE_MISMATCH could never fire: on the cohort it was written for it
    # flagged nothing, while 18 samples sat on the wrong reference. The unit test passed the
    # majority in directly and never exercised this lookup.
    ref_of = sample_references(args.metadata)
    by_ref = {}
    for sid, m in summ.items():
        ref, lin = ref_of.get(sid), clean_str(m.get("lineage"))
        if ref and lin and lin.lower() not in ("unclassified", "nan"):
            by_ref.setdefault(ref, []).append(lin.split(";")[0])
    ref_major = {r: max(set(v), key=v.count) for r, v in by_ref.items() if len(v) >= 3}
    dates = parse_collection_dates(args.metadata)
    # What each sample's reads cover, compared across the cohort: the stretches no read covers that
    # other samples do read (deletions), and how much of each bin could be called at all.
    # How close the samples are to each other, and which of them sit far from their own group (a
    # patient, a line) while close to a sample of another one.
    series_meta = parse_metadata(args.metadata)
    groups = {s: md.get("group") for s, md in series_meta.items()}
    distances = parse_pairs(args.snp_distances)
    outside = group_mismatches(distances, groups, args.cluster_snps)
    profiles, gene_lists = parse_depth_profiles(args.depth_profiles, min_len=args.deletion_min_len)
    coverage, del_tracks = (build_coverage(profiles, gene_lists, min_len=args.deletion_min_len,
                                           min_depth=args.deletion_min_depth)[:2]
                            if profiles else (None, {}))
    # A sample mapped against two references has a profile for each; its row in the report is the
    # one of the reference the samplesheet gives it first.
    profile_keys = {}
    for key in profiles:
        profile_keys.setdefault(key[0], []).append(key)

    def profile_of(sid):
        keys = profile_keys.get(sid, [])
        return (sid, ref_of.get(sid)) if (sid, ref_of.get(sid)) in profiles else (keys[0] if len(keys) == 1 else None)

    jsamples, counts = [], {"PASS": 0, "WARN": 0, "FAIL": 0}
    for sid, m in summ.items():
        anc = is_ancient(m)
        dmg = dmg_by_sample.get(sid)
        snp_prof = parse_profile(m.get("snp_profile"))
        pkey = profile_of(sid)
        verdict, flags = flag_sample(m, thr, snp_med, snp_sig, ancient=anc, anc_thr=anc_thr, dmg=dmg,
                                     ref_lineage=ref_major.get(ref_of.get(sid)))
        if sid in outside:
            flags.append("GROUP_MISMATCH")
            verdict = "FAIL" if verdict == "FAIL" else "WARN"
        counts[verdict] += 1
        jsamples.append({"s": sid, "v": verdict, "f": flags,
                         "lineage": clean_str(m.get("lineage")),
                         "hetf": het_frac(m),
                         "linf": lineage_fracs(m.get("lineage_counts")) or None,
                         "linc": lineage_counts_parsed(m.get("lineage_counts")) or None,
                         "dr": clean_str(m.get("drug_resistance")),
                         "anc": anc,
                         "dmg": dmg,
                         # the samplesheet's collection date; the summary's own 'date' is the day
                         # the pipeline ran, which says nothing about when the sample was taken
                         "date": dates.get(sid),
                         # the reference the sample was mapped to and the lineage its other samples
                         # type as: what a LINEAGE_MISMATCH flag compares, so the report can say so
                         "ref": ref_of.get(sid),
                         "ref_lin": ref_major.get(ref_of.get(sid)),
                         "miss": miss_by_sample.get(sid),
                         # what a GROUP_MISMATCH compares: the nearest sample of its own group and the
                         # nearest of another, with their distances
                         "grpd": outside.get(sid),
                         "trk": ({k: v for k, v in (("snp", snp_prof),
                                                    ("snpkb", snp_per_kb(snp_prof, profiles[pkey])
                                                     if pkey else None),
                                                    ("het", parse_profile(m.get("het_profile"))),
                                                    ("indel", parse_profile(m.get("indel_profile"))),
                                                    ("del", del_tracks.get(pkey))) if v is not None}
                                 or None),
                         "ann_db_error": (clean_str(m.get("ann_db_error")) == "yes"),
                         "m": dict({k: to_float(m.get(k)) for k in metric_keys + extra_keys + EFF_KEYS},
                                   **({"dose": dose_map.get(sid)} if dose_map else {}))})

    lineages = sorted({s["lineage"] for s in jsamples if s["lineage"]})
    n_ancient = sum(1 for s in jsamples if s["anc"])
    # Each reference's own genome (length, genes, masked regions), placed as its samples' profiles are
    # binned. The top-level fields keep the reference most samples were mapped to, which is also what a
    # single-reference run has always had there.
    ref_gff, ref_mask, ref_fai = pairs(args.ref_gff), pairs(args.ref_mask), pairs(args.ref_fai)
    genomes = {}
    if ref_gff or ref_mask or ref_fai:
        lengths = {}
        for sid, m in summ.items():
            v = to_float(m.get("length")) or to_float(m.get("genome_len"))
            if v and ref_of.get(sid):
                lengths[ref_of[sid]] = max(lengths.get(ref_of[sid], 0), int(v))
        genomes = build_genomes(ref_gff, ref_mask, ref_fai, lengths, NBINS)
    if genomes:
        n_on = {}
        for s in jsamples:
            if s.get("ref"):
                n_on[s["ref"]] = n_on.get(s["ref"], 0) + 1
        main_ref = sorted(genomes, key=lambda r: (-n_on.get(r, 0), r))[0]
        g0 = genomes[main_ref]
        genome_len = g0["len"] or genome_len
        genes, mask_bins, mask_pct, mask_iv = g0["genes"], g0["mask_bins"], g0["mask_pct"], g0["mask_iv"]
        gene_map = build_gene_map(list(ref_gff.values()))
    else:
        genes, gene_map = parse_gff(args.gff), build_gene_map(args.gff)
        mask_iv = parse_bed(args.mask_bed)
        mask_bins, mask_pct = mask_profile(mask_iv, genome_len)
    provenance = {}
    for kv in (args.provenance or []):
        if "=" in kv:
            k, v = kv.split("=", 1)
            provenance[k.strip()] = v.strip()

    _variants = load_variants(args)
    _sample_meta = parse_sample_meta(args.metadata)   # shared by the dynamics filter and the SNP matrix header
    # What each series gained since its first time point. The matrix is read only at the cells the
    # comparison needs, and without it a site not called at the start cannot be told from one not read.
    matrix_ok = bool(args.snp_matrix and os.path.exists(args.snp_matrix)
                     and not os.path.basename(args.snp_matrix).startswith("NO_FILE"))
    # Libraries of the same DNA, for how far down a minority call reproduces. The matrix is read once
    # for the cells both comparisons need.
    rep_col, rep_val, rep_reads = replicate_sets(args.metadata)
    rep_pairs, rep_shared = (replicate_pairs(rep_val, rep_reads, set(_variants), ref_of)
                             if rep_col else ([], 0))
    # A sample the QC fails or places outside its series (another lineage, far from its group)
    # cannot be the first time point every later one is read against.
    not_in_series = {s["s"] for s in jsamples
                     if s["v"] == "FAIL" or {"GROUP_MISMATCH", "LINEAGE_MISMATCH"} & set(s["f"])}
    _need = needed_cells(series_meta, _variants, excluded=not_in_series)
    for site, who in replicate_cells(rep_pairs, _variants).items():
        _need.setdefault(site, set()).update(who)
    _cells = read_matrix_cells(args.snp_matrix, _need) if matrix_ok else {}
    _series = build_series(series_meta, _variants, _cells, _sample_meta, args.min_dp, checked=matrix_ok,
                           excluded=not_in_series)
    _minority = build_minority(_variants, rep_pairs, rep_shared, _cells, rep_col, args.min_dp, checked=matrix_ok)
    _dynamics = build_dynamics(series_meta, _variants, _sample_meta)   # feeds dynamics + epistasis
    # front-load 'dose' so it survives the correlation matrix's top-N view
    dist_keys = (["dose"] + DIST) if dose_map else DIST
    now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    payload = {"generated": now, "version": clean_str(args.version) or "", "repo_url": REPO_URL,
               "counts": counts, "thresholds": thr, "dist": dist_keys,
               "genome_len": genome_len, "snp_density_ok": snp_density_ok, "defs": DEFS, "nbins": NBINS,
               "genes": genes, "lin_colors": parse_lineage_colors(args.lineage_colors),
               "gene_map": gene_map, "aa2_label": args.aa2_label,
               "genomes": genomes or None,
               "ref_name": provenance.get('reference', ''),   # mapping reference name for the SNP tables' headers
               "section_info": SECTION_INFO,
               "gene_burden": parse_gene_burden(args.gene_burden) or None,
               "dr": parse_dr(args.dr_report),
               "gconv": parse_gene_conversion(args.gene_conversion),
               "kraken": parse_kraken(args.kraken),
               "pnps": parse_pnps(args.pnps) or None,
               "mask_bins": mask_bins, "mask_pct": mask_pct,
               "mask_iv": (mask_iv[:5000] if mask_iv else None),
               "lineages": lineages, "lin_present": len(lineages) > 0,
               "n_ancient": n_ancient, "anc_thresholds": anc_thr, "provenance": provenance,
               "dynamics": _dynamics,
               "epistasis": build_epistasis(_dynamics),
               "snp_matrix": build_snp_matrix(_variants, provenance.get('reference', '')),
               "coverage": coverage,
               "relatedness": build_relatedness(distances, groups, args.cluster_snps),
               "series": _series,
               "minority": _minority,
               "sample_meta": _sample_meta,
               "metrics": [{"key": k, "label": l, "kind": kind, "dir": d} for k, l, kind, d in METRICS],
               "extra": extra_metrics,
               "samples": jsamples}
    return payload, jsamples, counts


def write_flags(path, jsamples):
    """The machine-readable verdict table: one row per sample, FAIL first, then WARN, then PASS."""
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("sample\tverdict\tmean_depth\tbreadth_pct\tmapped_pct\tmissing_pct\tiupac_pct\tsnps\tti_tv\tflags\n")
        for s in sorted(jsamples, key=lambda x: (x["v"] != "FAIL", x["v"] != "WARN", x["s"])):
            m = s["m"]
            def g(k):
                v = m.get(k)
                return "" if v is None else f"{v:.4g}"
            fh.write("\t".join([s["s"], s["v"], g("mean_depth"), g("breadth_pct"), g("mapped_pct"),
                                g("missing_pct"), g("iupac_pct"), g("snps"), g("ti_tv"),
                                ",".join(s["f"]) or "."]) + "\n")


def sample_references(path):
    """{sampleId: refId} from the samplesheet. Empty when there is no samplesheet or no refId column.

    Read here rather than through parse_sample_meta, which drops refId on purpose because it is a
    pipeline column and not an annotation to display. A sample merged from several runs appears
    once per run with the same refId, so the first row per sample is enough.
    """
    import csv
    out = {}
    if not path:
        return out
    try:
        with open(path, encoding="utf-8", errors="replace") as fh:
            for r in csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"):
                sid = (r.get("sampleId") or r.get("sample") or "").strip()
                ref = (r.get("refId") or "").strip()
                if sid and ref:
                    out.setdefault(sid, ref)
    except OSError:
        return {}
    return out


def main():
    args = build_parser().parse_args()
    thr = {k: getattr(args, k) for k in DEF}
    anc_thr = {k: getattr(args, "anc_" + k) for k in ANC_DEF}
    payload, jsamples, counts = build_payload(args, thr, anc_thr)

    os.makedirs(os.path.dirname(os.path.abspath(args.out_html)) or ".", exist_ok=True)
    with open(args.out_html, "w", encoding="utf-8") as fh:
        fh.write(build_html(args.title, payload))
    write_flags(args.out_flags, jsamples)
    if args.out_deletions and payload.get("coverage"):
        write_deletions(args.out_deletions, payload["coverage"]["regions"])

    sys.stderr.write(f"[qc_report] {len(jsamples)} samples -> {counts['PASS']} PASS, {counts['WARN']} WARN, "
                     f"{counts['FAIL']} FAIL. Wrote {args.out_html} + {args.out_flags}\n")
    for s in sorted([x for x in jsamples if x["v"] == "FAIL"], key=lambda x: x["s"]):
        sys.stderr.write(f"    FAIL  {s['s']:<40s} [{','.join(s['f'])}]\n")
    if args.gate and counts["FAIL"]:
        sys.stderr.write(f"[qc_report] GATE: {counts['FAIL']} sample(s) FAIL -> blocking.\n")
        sys.exit(1)


if __name__ == "__main__":
    main()
