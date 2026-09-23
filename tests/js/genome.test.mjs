// The genome the landscape draws when a cohort is mapped against several references (09_genome_genes.js).
//
// Each sample's profile is binned over its own reference, so the page shows one reference at a time: its
// samples, and its own length, genes and masked regions. A report without per-reference genomes keeps the
// single one it always had.

import { strict as assert } from "node:assert";
import { describe, it } from "node:test";
import { emptyReport, loadFunctions } from "./harness.mjs";

const FNS = ["genome", "genomeRef", "genomeRefs", "onGenome", "genomeSwitch", "insGenome"];

function sample(s, ref) {
  return { s, ref, v: "PASS", f: [], m: { callable_pct: 95 }, miss: [0.1, 0.2], trk: { snp: [1, 2] } };
}

function twoRefs() {
  return {
    ...emptyReport(),
    samples: [sample("a1", "L7ref"), sample("b1", "A4ref"), sample("b2", "A4ref")],
    genome_len: 4409697, genes: [{ name: "l7gene", start: 1, end: 9 }],
    genomes: {
      L7ref: { len: 4409697, genes: [{ name: "l7gene", start: 1, end: 9 }], mask_bins: null, mask_pct: 0, mask_iv: null },
      A4ref: { len: 4346880, genes: [{ name: "a4gene", start: 5, end: 50 }], mask_bins: [0, 1], mask_pct: 3, mask_iv: [[1, 10]] },
    },
  };
}

describe("the genome in view", () => {
  it("opens on the reference most samples were mapped to, with its own genes and length", () => {
    const fn = loadFunctions(FNS, twoRefs());
    assert.equal(fn.genomeRef(), "A4ref");
    assert.equal(fn.genome().len, 4346880);
    assert.equal(fn.genome().genes[0].name, "a4gene");
  });

  it("shows only the samples of the reference in view, and switches with it", () => {
    const rep = twoRefs();
    const fn = loadFunctions(FNS, rep);
    assert.deepEqual(rep.samples.filter(fn.onGenome).map((s) => s.s), ["b1", "b2"]);
    fn.genomeSwitch("L7ref");
    assert.equal(fn.genomeRef(), "L7ref");
    assert.deepEqual(rep.samples.filter(fn.onGenome).map((s) => s.s), ["a1"]);
    assert.equal(fn.genome().genes[0].name, "l7gene");
  });

  it("keeps the single genome of a report without per-reference ones", () => {
    const rep = { ...twoRefs(), genomes: null };
    const fn = loadFunctions(FNS, rep);
    assert.equal(fn.genomeRef(), null);
    assert.equal(fn.genome().genes[0].name, "l7gene");
    assert.ok(rep.samples.every(fn.onGenome));
  });

  it("describes the callability of the reference in view, by name", () => {
    const fn = loadFunctions(FNS, twoRefs());
    assert.match(fn.insGenome(), /4\.35 Mb reference A4ref/);
  });
});
