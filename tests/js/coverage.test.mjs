// The deletions panel places a region on the genome landscape and searches it. Both read the
// payload alone, so they are checked on small hand-built cohorts.

import { strict as assert } from "node:assert";
import { describe, it } from "node:test";
import { emptyReport, loadFunctions } from "./harness.mjs";

function report(coverage, more = {}) {
  return { ...emptyReport(), samples: [], metrics: [], nbins: 10, genome_len: 1000, coverage, ...more };
}

describe("delBins, a region on the landscape", () => {
  it("uses the landscape's own binning over one contig", () => {
    const rep = report({ refs: { R: [["chr", 1000]] }, regions: [] });
    const fn = loadFunctions(["delBins"], rep);
    assert.deepEqual(JSON.parse(JSON.stringify(fn.delBins({ ref: "R", contig: "chr", start: 250, end: 449 }))), [2, 4]);
  });

  it("lays the contigs of a reference end to end, as the landscape does", () => {
    // 600 bp of chromosome then a 400 bp plasmid: plasmid position 1 is genome position 601.
    const rep = report({ refs: { R: [["chr", 600], ["plasmid", 400]] }, regions: [] });
    const fn = loadFunctions(["delBins"], rep);
    assert.deepEqual(JSON.parse(JSON.stringify(fn.delBins({ ref: "R", contig: "plasmid", start: 1, end: 400 }))), [6, 9]);
  });

  it("scales a reference of another length to its own genome", () => {
    const rep = report({ refs: { A: [["a", 1000]], B: [["b", 2000]] }, regions: [] });
    const fn = loadFunctions(["delBins"], rep);
    assert.deepEqual(JSON.parse(JSON.stringify(fn.delBins({ ref: "B", contig: "b", start: 1000, end: 1999 }))), [5, 9]);
  });
});

describe("delMatch, searching the deletions", () => {
  const region = { ref: "R", contig: "chr", start: 1200, end: 3400, genes: ["katG"], samples: [["S0", 1200, 3000]] };
  const fn = loadFunctions(["delMatch"], report({ refs: {}, regions: [] }));

  it("finds a region by sample, gene or coordinates", () => {
    assert.equal(fn.delMatch(region, "s0"), true);
    assert.equal(fn.delMatch(region, "katg"), true);
    assert.equal(fn.delMatch(region, "chr:1200"), true);
    assert.equal(fn.delMatch(region, "rpob"), false);
  });
});
