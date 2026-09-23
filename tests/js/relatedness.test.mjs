// The relatedness page clusters the samples and orders its heatmap from the pairwise distances
// alone, so both are checked on hand-built references.

import { strict as assert } from "node:assert";
import { describe, it } from "node:test";
import { emptyReport, loadFunctions } from "./harness.mjs";

const host = (x) => JSON.parse(JSON.stringify(x));

// Upper triangle, row by row, of a symmetric matrix given as {"a,b": snps}.
function ref(samples, d, cmp = 100) {
  const snps = [], c = [];
  for (let i = 0; i < samples.length; i++)
    for (let j = i + 1; j < samples.length; j++) {
      snps.push(d[`${samples[i]},${samples[j]}`] ?? d[`${samples[j]},${samples[i]}`] ?? 1000);
      c.push(typeof cmp === "object" ? (cmp[`${samples[i]},${samples[j]}`] ?? 100) : cmp);
    }
  return { samples, variable: 100, snps, cmp: c };
}

const fns = loadFunctions(["relIdx", "relD", "relClusters", "relOrder"], { ...emptyReport(), samples: [], metrics: [] });

describe("relClusters, single linkage at a threshold", () => {
  it("chains pairs within the threshold into one cluster", () => {
    // a-b 3 and b-c 4 join a, b and c although a-c is 7; d is far from all of them.
    const r = ref(["a", "b", "c", "d"], { "a,b": 3, "b,c": 4, "a,c": 7 });
    assert.deepEqual(host(fns.relClusters(r, 5, 50)), [[0, 1, 2], [3]]);
  });

  it("leaves out a pair that compared too few positions", () => {
    const r = ref(["a", "b"], { "a,b": 0 }, { "a,b": 10 });
    assert.deepEqual(host(fns.relClusters(r, 5, 50)), [[0], [1]], "0 SNPs over 10% of the positions says nothing");
  });
});

describe("relOrder, the heatmap's order", () => {
  it("keeps every cluster contiguous", () => {
    // Two tight pairs interleaved in the input: a,c close and b,d close.
    const r = ref(["a", "b", "c", "d"], { "a,c": 1, "b,d": 2 });
    const order = host(fns.relOrder(r, 50));
    const pos = Object.fromEntries(order.map((x, i) => [x, i]));
    assert.equal(Math.abs(pos[0] - pos[2]), 1, "a and c sit together");
    assert.equal(Math.abs(pos[1] - pos[3]), 1, "b and d sit together");
    assert.deepEqual([...order].sort(), [0, 1, 2, 3]);
  });

  it("indexes the upper triangle row by row", () => {
    const r = ref(["a", "b", "c"], { "a,b": 1, "a,c": 2, "b,c": 3 });
    assert.equal(fns.relD(r, 0, 2), 2);
    assert.equal(fns.relD(r, 2, 1), 3);
    assert.equal(fns.relD(r, 1, 1), 0);
  });
});
