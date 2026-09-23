// The summary page and the read-outs that share its numbers. Every builder reads the payload
// alone, so each is checked here on a small cohort built for the case, the way the report sees it.

import { strict as assert } from "node:assert";
import { describe, it } from "node:test";
import { emptyReport, loadFunctions } from "./harness.mjs";

const METRICS = [
  ["mean_depth", "Depth", "float", "hi_good"], ["breadth_pct", "Breadth %", "pct", "hi_good"],
  ["missing_pct", "Missing %", "pct", "hi_bad"], ["mapped_pct", "Mapped %", "pct", "hi_good"],
  ["duplication_pct", "Dup %", "pct", "hi_bad"], ["iupac_pct", "IUPAC %", "pct", "hi_bad"],
  ["snps", "SNPs", "int", "neu"], ["ti_tv", "Ti/Tv", "float", "neu"],
].map(([key, label, kind, dir]) => ({ key, label, kind, dir }));

const CLEAN = { mean_depth: 60, breadth_pct: 99.5, missing_pct: 1, mapped_pct: 99, duplication_pct: 5, iupac_pct: 0.1, snps: 10 };

function sample(s, lineage, m = {}, extra = {}) {
  return { s, v: "PASS", f: [], lineage, hetf: null, linf: null, linc: null, anc: false, dmg: null,
           m: { ...CLEAN, ...m }, ...extra };
}

function report(samples, more = {}) {
  return {
    ...emptyReport(),
    samples,
    metrics: METRICS.map((x) => ({ ...x })),
    lineages: [...new Set(samples.map((s) => s.lineage).filter(Boolean))].sort(),
    lin_present: samples.some((s) => s.lineage),
    thresholds: { depth_min: 10, breadth_min: 90, missing_max: 10, dup_max: 40, iupac_max: 2, mapping_min: 80,
                  titv_min: 0, snp_z: 3, het_max_frac: 3, mixed_min_frac: 2 },
    counts: { PASS: samples.length, WARN: 0, FAIL: 0 },
    ...more,
  };
}

function load(names, rep = report([])) {
  return loadFunctions(names, rep);
}

// Values built inside the vm context carry its own Array.prototype; copy them out before comparing.
const host = (x) => JSON.parse(JSON.stringify(x));

describe("linLabel, a lineage call as a reader says it", () => {
  const fn = load(["linLabel", "untyped", "linMain"]);

  it("keeps the finest level and the species label", () => {
    assert.equal(fn.linLabel("A4;M_bovis;A4.7;A4.7.1"), "A4.7.1 · M. bovis");
    assert.equal(fn.linLabel("L6;L6.2;L6.2.1"), "L6.2.1");
    assert.equal(fn.linLabel("L7"), "L7");
    assert.equal(fn.linLabel(null), "NA");
  });

  it("compares on the top level, and knows an untyped call when it sees one", () => {
    assert.equal(fn.linMain("A4;M_bovis;A4.7"), "A4");
    assert.equal(fn.untyped("Unclassified"), true);
    assert.equal(fn.untyped("NA"), true);
    assert.equal(fn.untyped("L4.1"), false);
  });
});

describe("coupled, the metric pairs whose correlation is arithmetic", () => {
  const fn = load(["coupled"]);

  it("pairs a percent with its complement, and every variant count with the SNPs", () => {
    assert.equal(fn.coupled("mapped_pct", "unmapped_pct"), true);
    assert.equal(fn.coupled("missing_pct", "callable_pct"), true);
    assert.equal(fn.coupled("snps", "ann_high"), true);
  });

  it("leaves two different quantities alone", () => {
    assert.equal(fn.coupled("mean_depth", "breadth_pct"), false);
    assert.equal(fn.coupled("iupac_pct", "snps"), false);
  });
});

describe("nv, a value as a sentence quotes it", () => {
  const fn = load(["nv"]);

  it("keeps three significant digits without padding zeros", () => {
    assert.equal(fn.nv(2.2), "2.2");
    assert.equal(fn.nv(0.05), "0.05");
    assert.equal(fn.nv(93.73), "93.7");
    assert.equal(fn.nv(260.55), "261");
    assert.equal(fn.nv(2195), "2,195");
    assert.equal(fn.nv(null), "NA");
  });
});

describe("krkState, how much of a sample is the organism the cohort is about", () => {
  const fn = load(["krkState", "krkPurity"]);

  it("reads the share of the classified reads in the target clade", () => {
    assert.equal(fn.krkState({ target_pct: 1.2, unclassified: 3 }), "off");
    assert.equal(fn.krkState({ target_pct: 75, unclassified: 3 }), "mixed");
    assert.equal(fn.krkState({ target_pct: 96, unclassified: 22 }), "uncl");
    assert.equal(fn.krkState({ target_pct: 96, unclassified: 3 }), "ok");
  });
});

describe("divScreen, the reference-bias screen", () => {
  it("switches itself off when the cohort sits on its reference", () => {
    // mapped to their own ancestor, the samples carry a handful of SNPs per Mb: fewer than the rest is noise
    const fn = load(["divScreen"], report([], { snp_density_ok: true }));
    const scr = fn.divScreen([1.2, 1.5, 0.5, 2.1, 1.1, 0.9]);
    assert.equal(scr.on, false);
  });

  it("runs, with a cut below the median, when the reference is another strain", () => {
    const fn = load(["divScreen"], report([], { snp_density_ok: true }));
    const scr = fn.divScreen([210, 230, 250, 240, 220, 60]);
    assert.equal(scr.on, true);
    assert.ok(scr.lo > 0 && scr.lo < scr.med, `cut ${scr.lo} should sit below the median ${scr.med}`);
  });
});

describe("recompute, the live re-flagging", () => {
  it("keeps a flag the pipeline decided from evidence the page does not hold", () => {
    // LINEAGE_MISMATCH compares a sample with the rest of its reference. Rebuilding the flags from the
    // thresholds used to drop it, so the page hid every wrong-reference sample.
    const rep = report([
      sample("A", "L7"), sample("B", "L7"), sample("C", "L7"),
      sample("X", "A4;M_bovis", {}, { f: ["LINEAGE_MISMATCH"], v: "WARN", ref: "R7", ref_lin: "L7" }),
    ]);
    const fn = load(["recompute", "flagReason"], rep);
    fn.recompute();
    const x = rep.samples.find((s) => s.s === "X");
    assert.ok(x.f.includes("LINEAGE_MISMATCH"));
    assert.equal(x.v, "WARN");
    assert.match(fn.flagReason(x, "LINEAGE_MISMATCH")[1], /types as A4; the rest of R7 types as L7/);
  });

  it("still fails a sample on the thresholds in force", () => {
    const rep = report([sample("A", "L7", { mean_depth: 2 }), sample("B", "L7"), sample("C", "L7")]);
    const fn = load(["recompute"], rep);
    fn.recompute();
    assert.equal(rep.samples[0].v, "FAIL");
    assert.ok(rep.samples[0].f.includes("LOW_DEPTH"));
  });
});

describe("insLineages", () => {
  it("counts mixed samples from the flag, not from every typed sample's lineage breakdown", () => {
    // every typed sample carries s.linf; counting those once reported 176 mixed samples in a cohort with none
    const samples = ["A", "B", "C", "D"].map((id) => sample(id, "L7", {}, { linf: { L7: 1 }, linc: { L7: 277 } }));
    const fn = load(["insLineages"], report(samples));
    const html = fn.insLineages();
    assert.match(html, /No sample carries a second lineage/);
    assert.doesNotMatch(html, /mixed/i);
  });
});

describe("drSummary, resistance a lineage shares against what some samples carry", () => {
  function cohort() {
    const samples = [];
    for (let i = 0; i < 6; i++) samples.push(sample(`A${i}`, "A4;M_bovis;A4.7"));
    for (let i = 0; i < 6; i++) samples.push(sample(`L${i}`, "L7"));
    const call = (s, drug, gene, mutation, gn = 1) => ({ s, drug, dr: [drug], gene, mutation, gn, grade: `${gn}) Assoc w R`, af: 1, dp: 50 });
    const calls = samples.filter((s) => s.s.startsWith("A")).map((s) => call(s.s, "PZA", "pncA", "H57D"));
    calls.push(call("L0", "PZA", "pncA", "H57D"));                 // the A4 marker, called in one L7 sample
    calls.push(call("A1", "BDQ", "atpE", "I66M", 2), call("A2", "BDQ", "atpE", "I66M", 2));
    calls.push(call("L3", "RIF", "rpoB", "S450L", 3));           // grade 3: not a resistance call
    return report(samples, { dr: { samples: samples.map((s) => s.s), drugs: ["RIF", "PZA", "BDQ"], calls } });
  }

  it("sets aside a mutation every sample of a lineage carries", () => {
    const fn = load(["drSummary"], cohort());
    const d = fn.drSummary();
    assert.deepEqual(host(d.wide.map((m) => `${m.gene} ${m.mutation} in ${m.wide}`)), ["pncA H57D in A4"]);
    assert.deepEqual(host(d.acquired.map((m) => `${m.gene} ${m.mutation}`)), ["atpE I66M"]);
  });

  it("counts the samples and drugs of what the lineage does not share", () => {
    const fn = load(["drSummary", "findResistance"], cohort());
    const d = fn.drSummary();
    assert.equal(d.nCarriers, 2);
    assert.deepEqual(host(d.drugs.map((x) => [x.drug, x.n])), [["BDQ", 2]]);
    const f = fn.findResistance();
    assert.match(f.head, /2<\/b> samples carry resistance-associated mutations their lineage does not share/);
    assert.match(f.body, /pncA H57D<\/b> \(PZA\) in all 6 A4 samples, and in 1 of 6 L7/);
  });
});

describe("findQC", () => {
  it("names the rule each excluded sample broke", () => {
    const rep = report([sample("A", "L7", { mean_depth: 2 }), sample("B", "L7", { missing_pct: 40 }), sample("C", "L7")]);
    const fn = load(["recompute", "findQC"], rep);
    fn.recompute();
    const f = fn.findQC();
    assert.match(f.head, /<b>1<\/b> of 3 samples are ready to use; <b class="tone-bad">2<\/b> should be excluded/);
    assert.match(f.body, /depth below 10/);
    assert.match(f.body, /more than 10% of the consensus missing/);
  });
});

describe("findGconv", () => {
  it("counts the events of samples the QC keeps apart from those only failing samples carry", () => {
    // A mixed or contaminated culture carries the donor's bases for reasons of its own; a summary
    // counting its tracts with the rest overstates what the run found.
    const t = (s, event, extra = {}) => ({ s, event, verdict: "gene_conversion", rep: 1, contig: "chr", start: 100, bp_reads: 1, ...extra });
    const rep = report([sample("OK1", "L7"), sample("OK2", "L7"), sample("REV", "L7", {}, { v: "WARN", f: ["LOW_DEPTH"] }),
                        sample("BAD", "L7", {}, { v: "FAIL", f: ["HIGH_MISSING"] })], {
      gconv: { tracts: [t("OK1", "chr:1"), t("OK2", "chr:1"), t("REV", "chr:4", { bp_reads: 0 }), t("BAD", "chr:2"), t("BAD", "chr:3", { verdict: "ambiguous" })] },
    });
    const fn = load(["gconvSummary", "findGconv"], rep);
    const g = fn.gconvSummary();
    assert.equal(g.nPassEvents, 2, "a WARN sample is one to review, not one to exclude");
    assert.equal(g.failOnly, 1);
    const f = fn.findGconv();
    assert.match(f.head, /<b>2<\/b> gene-conversion events called in samples the QC does not fail/);
    assert.match(f.body, /1 more is only in samples the QC fails \(BAD\)/);
    assert.match(f.body, /2 of the 3 events are backed by a read crossing a breakpoint/,
      "counted in events like the headline, not in the calls of each sample");
  });

  it("does not say 'more' when no sample the QC keeps carries an event", () => {
    const rep = report([sample("OK1", "L7"), sample("BAD", "L7", {}, { v: "FAIL", f: ["HIGH_MISSING"] })], {
      gconv: { tracts: [{ s: "BAD", event: "chr:2", verdict: "gene_conversion", rep: 1, bp_reads: 1 }] },
    });
    const f = load(["gconvSummary", "findGconv"], rep).findGconv();
    assert.match(f.head, /<b>0<\/b> gene-conversion events called in samples the QC does not fail/);
    assert.match(f.body, /^1 is only in samples the QC fails \(BAD\).*That event is backed/);
  });

  it("names only the failed samples of events no kept sample shares", () => {
    const t = (s, event) => ({ s, event, verdict: "gene_conversion", rep: 1, contig: "chr", start: 100, bp_reads: 0 });
    const rep = report([sample("A", "L7"), sample("B", "L7", {}, { v: "FAIL", f: ["HIGH_MISSING"] }),
                        sample("C", "L7", {}, { v: "FAIL", f: ["HIGH_MISSING"] })], {
      gconv: { tracts: [t("A", "chr:1"), t("B", "chr:1"), t("C", "chr:2")] },
    });
    const g = load(["gconvSummary"], rep).gconvSummary();
    assert.deepEqual(host(g.failSamples), ["C"], "B shares chr:1 with a kept sample");
  });
});

describe("seriesGainSummary, what the series gained since their first time point", () => {
  it("summarises series without the samples the chart leaves out, and with a true median", () => {
    const row = (s, time, n) => ({ s, time, tnum: +time, new: Array(n).fill("c:1"), risen: [], unknown: 0, lost: [] });
    const rep = report([sample("a1", "L7"), sample("a2", "L7"), sample("bad", "L7", {}, { v: "WARN", f: ["LINEAGE_MISMATCH"] }),
                        sample("b1", "L7")], {
      series: { checked: true, groups: [
        { group: "A", rows: [row("a1", "6", 0), row("a2", "6", 10), row("bad", "9", 2462)] },
        { group: "B", rows: [row("b1", "6", 3)] },
      ] },
    });
    const g = host(load(["seriesGainSummary", "serOutside"], rep).seriesGainSummary());
    assert.equal(g.max, 5, "A's last counted time point is 6, median of 0 and 10");
    assert.equal(g.median, 4, "the median of 5 and 3");
  });
});

describe("minFloor, the allele fraction from which calls are reproduced", () => {
  it("reads the noise floor from the top band down", () => {
    // 50/85/60/70/95/97%: the first passing band is 0.1-0.2, but 0.2-0.5 fall short again.
    const tested = [100, 100, 100, 100, 100, 100];
    const rep = report([], { minority: { edges: [0, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9],
      replicates: { tested, reproduced: [50, 85, 60, 70, 95, 97] } } });
    assert.equal(load(["minFloor"], rep).minFloor(rep.minority), 4, "from 0.5 up");
  });

  it("gives no floor when the band just below fixation already fails", () => {
    const rep = report([], { minority: { edges: [0, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9],
      replicates: { tested: [100, 100, 100, 100, 100, 100], reproduced: [90, 90, 90, 90, 90, 50] } } });
    assert.equal(load(["minFloor"], rep).minFloor(rep.minority), -1);
  });

  it("gives none at all when no band was compared often enough to judge", () => {
    const M = { edges: [0, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9], hist: [4, 3, 2, 1, 1, 1], le3: 2, with_dp: 12,
                replicates: { column: "dna_id", tested: [5, 4, 3, 2, 1, 1], reproduced: [1, 1, 1, 1, 1, 1] } };
    const fn = load(["minFloor", "findMinority"], report([], { minority: M }));
    assert.equal(fn.minFloor(M), null);
    assert.match(fn.findMinority().body, /too few to say where the noise ends/);
  });
});
