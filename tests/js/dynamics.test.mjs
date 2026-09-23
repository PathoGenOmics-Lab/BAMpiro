// Matching a trajectory's variant to the drug-resistance catalogue (10_dynamics.js).
//
// The catalogue writes a protein change as I66M (or LoF for any loss of function of the gene) and names
// genes as H37Rv does. A variant carries HGVS (p.Ile66Met), in the numbering and gene names of the
// reference it was mapped against, and in H37Rv's where the pipeline lifted it there.

import { strict as assert } from "node:assert";
import { describe, it } from "node:test";
import { emptyReport, loadFunctions } from "./harness.mjs";

const calls = [
  { gene: "atpE", mutation: "I66M", drug: "BDQ", gn: 1 },
  { gene: "katG", mutation: "c.609C>T", drug: "INH", gn: 3 },
  { gene: "pncA", mutation: "LoF", drug: "PZA", gn: 1 },
  { gene: "katG", mutation: "W191!", drug: "INH", gn: 1 },
];

function load() {
  return loadFunctions(["dynDRindex", "dynDRmatch", "aaShort", "aaLoF"], { ...emptyReport(), dr: { calls } });
}

describe("aaShort, HGVS to the catalogue's short form", () => {
  it("turns three-letter codes into one", () => {
    const fn = load();
    assert.equal(fn.aaShort("p.Ile66Met"), "I66M");
    assert.equal(fn.aaShort("p.Trp191*"), "W191*");
    assert.equal(fn.aaShort("p.Gln53Ter"), "Q53*");
    assert.equal(fn.aaShort("p.Lys12fs"), "", "a frameshift has no short form");
    assert.equal(fn.aaShort(""), "");
  });

  it("knows a loss of function", () => {
    const fn = load();
    assert.ok(fn.aaLoF("p.Lys12fs") && fn.aaLoF("p.Trp191*") && fn.aaLoF("p.Gln53Ter"));
    assert.ok(!fn.aaLoF("p.Ile66Met"));
  });
});

describe("dynDRmatch", () => {
  it("matches a missense change written in HGVS", () => {
    const fn = load();
    const hit = fn.dynDRmatch(fn.dynDRindex(), { gene: "atpE", aa: "p.Ile66Met" });
    assert.equal(hit && hit.drug, "BDQ");
  });

  it("matches in the H37Rv gene and numbering where the reference names and numbers its own way", () => {
    const fn = load();
    const v = { gene: "E1ASM0057_01312", aa: "p.Ile67Met", gene_h37rv: "atpE", aa_h37rv: "p.Ile66Met" };
    assert.equal(fn.dynDRmatch(fn.dynDRindex(), v).drug, "BDQ");
  });

  it("matches a stop written with ! and a loss of function by gene", () => {
    const fn = load();
    const idx = fn.dynDRindex();
    assert.equal(fn.dynDRmatch(idx, { gene: "katG", aa: "p.Trp191*" }).mutation, "W191!");
    assert.equal(fn.dynDRmatch(idx, { gene: "pncA", aa: "p.Gln10fs" }).drug, "PZA");
  });

  it("does not match another change of the same gene", () => {
    const fn = load();
    assert.equal(fn.dynDRmatch(fn.dynDRindex(), { gene: "atpE", aa: "p.Ala63Val" }), null);
    assert.equal(fn.dynDRmatch(fn.dynDRindex(), { gene: "atpE", aa: "" }), null);
  });
});
