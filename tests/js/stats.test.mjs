// The report computes its statistics in hand-written ES5 (no library is available inside a
// self-contained HTML file). These tests check that arithmetic against reference values from
// SciPy 1.15.2 and statsmodels, so a regression is caught against an independent implementation
// rather than against a snapshot of our own output.
//
// Reference values were produced with:
//   scipy.stats.spearmanr, scipy.stats.kruskal, scipy.stats.chi2.sf, scipy.stats.t.sf,
//   statsmodels.stats.multitest.multipletests(method="fdr_bh")

import { strict as assert } from "node:assert";
import { before, describe, it } from "node:test";
import { loadFunctions } from "./harness.mjs";

const TOL = 1e-9;

// Arrays built inside the vm context carry that context's Array.prototype, so deepStrictEqual
// rejects them as not reference-equal. Copy into a host array before structural comparison.
const arr = (x) => Array.from(x);

function closeTo(actual, expected, tol = TOL, what = "") {
  assert.equal(typeof actual, "number", `${what}: expected a number, got ${actual}`);
  const diff = Math.abs(actual - expected);
  const scale = Math.max(1, Math.abs(expected));
  assert.ok(
    diff / scale <= tol,
    `${what}: ${actual} is not within ${tol} (relative) of ${expected}`,
  );
}

let fn;
before(() => {
  fn = loadFunctions([
    "rankvec", "spearman", "spearmanP", "studentP", "betai",
    "chiSqP", "kruskalWallis", "bhFDR", "pfmt",
  ]);
});

describe("rankvec", () => {
  it("ranks from 1 and preserves input order", () => {
    assert.deepEqual(arr(fn.rankvec([3, 1, 2])), [3, 1, 2]);
    assert.deepEqual(arr(fn.rankvec([10, 20, 30])), [1, 2, 3]);
  });

  it("averages ranks across ties, like scipy's default", () => {
    assert.deepEqual(arr(fn.rankvec([5, 5, 5])), [2, 2, 2]);
    assert.deepEqual(arr(fn.rankvec([1, 2, 2, 3])), [1, 2.5, 2.5, 4]);
  });

  it("returns an empty array for empty input", () => {
    assert.deepEqual(arr(fn.rankvec([])), []);
  });
});

describe("spearman", () => {
  // rho values from scipy.stats.spearmanr
  const cases = [
    ["perfect positive", [1, 2, 3, 4], [1, 2, 3, 4], 1.0],
    ["perfect negative", [1, 2, 3, 4], [4, 3, 2, 1], -1.0],
    ["noisy", [1, 2, 3, 4, 5, 6, 7, 8], [2, 1, 4, 3, 6, 5, 8, 7], 0.9047619047619048],
    ["with ties", [1, 2, 2, 3, 4, 5], [2, 3, 3, 1, 5, 4], 0.5882352941176471],
  ];

  for (const [name, x, y, expected] of cases) {
    it(`matches scipy for ${name}`, () => {
      closeTo(fn.spearman(x, y), expected, TOL, name);
    });
  }

  it("returns null below four points", () => {
    assert.equal(fn.spearman([1, 2, 3], [1, 2, 3]), null);
    assert.equal(fn.spearman([], []), null);
  });

  it("returns null when a vector has no variance", () => {
    assert.equal(fn.spearman([1, 1, 1, 1], [1, 2, 3, 4]), null);
    assert.equal(fn.spearman([1, 2, 3, 4], [7, 7, 7, 7]), null);
  });
});

describe("spearmanP", () => {
  // scipy.stats.spearmanr p-values (scipy uses the same Student-t approximation)
  it("matches scipy for a noisy correlation", () => {
    const rho = fn.spearman([1, 2, 3, 4, 5, 6, 7, 8], [2, 1, 4, 3, 6, 5, 8, 7]);
    closeTo(fn.spearmanP(rho, 8), 0.0020082755054294677, 1e-9, "noisy p");
  });

  it("matches scipy in the presence of ties", () => {
    const rho = fn.spearman([1, 2, 2, 3, 4, 5], [2, 3, 3, 1, 5, 4]);
    closeTo(fn.spearmanP(rho, 6), 0.21941787095461004, 1e-9, "tied p");
  });

  it("returns exactly 0 for a perfect correlation, avoiding a division by zero", () => {
    assert.equal(fn.spearmanP(1, 10), 0);
    assert.equal(fn.spearmanP(-1, 10), 0);
  });

  it("returns null when rho is null or there are too few points", () => {
    assert.equal(fn.spearmanP(null, 10), null);
    assert.equal(fn.spearmanP(0.5, 2), null);
  });
});

describe("kruskalWallis", () => {
  // H and p from scipy.stats.kruskal, which applies the same tie correction
  const cases = [
    ["two groups", [[1, 2], [3, 4]], 2.3999999999999986, 0.12133525035848222],
    ["three groups", [[1, 2, 3], [4, 5, 6], [7, 8, 9]], 7.200000000000003, 0.02732372244729252],
    ["tied values", [[1, 1, 2], [2, 3, 3], [3, 4, 4]], 6.4424778761061905, 0.03990558707075906],
    [
      "uneven group sizes",
      [[2.1, 3.4, 1.9, 4.0], [5.5, 6.1, 4.9], [9.0, 8.2, 7.7, 8.8, 9.4]],
      9.692307692307693,
      0.0078585446701517,
    ],
  ];

  for (const [name, groups, h, p] of cases) {
    it(`matches scipy for ${name}`, () => {
      const got = fn.kruskalWallis(groups);
      closeTo(got.H, h, TOL, `${name} H`);
      closeTo(got.p, p, 1e-8, `${name} p`);
      assert.equal(got.k, groups.length);
      assert.equal(got.df, groups.length - 1);
    });
  }

  it("gives H = 0 and p = 1 when every value is identical", () => {
    const got = fn.kruskalWallis([[1, 1], [1, 1]]);
    assert.equal(got.H, 0);
    assert.equal(got.p, 1);
  });

  it("returns null with fewer than two groups or fewer than three values", () => {
    assert.equal(fn.kruskalWallis([[1, 2, 3]]), null);
    assert.equal(fn.kruskalWallis([[1], [2]]), null);
    assert.equal(fn.kruskalWallis([]), null);
  });
});

describe("chiSqP", () => {
  // scipy.stats.chi2.sf
  const cases = [
    [3.84, 1, 0.050043521248705085],
    [2.4, 1, 0.1213352503584821],
    [5.99, 2, 0.05003662708658629],
    [100, 2, 1.9287498479639183e-22],
    [12.5, 7, 0.0852692751582693],
  ];

  for (const [x, df, expected] of cases) {
    it(`matches scipy at x=${x}, df=${df}`, () => {
      closeTo(fn.chiSqP(x, df), expected, 1e-8, `chi2(${x},${df})`);
    });
  }

  it("returns 1 for degenerate input instead of NaN", () => {
    assert.equal(fn.chiSqP(0, 1), 1);
    assert.equal(fn.chiSqP(-1, 3), 1);
    assert.equal(fn.chiSqP(NaN, 3), 1);
    assert.equal(fn.chiSqP(5, 0), 1);
  });
});

describe("studentP", () => {
  // 2 * scipy.stats.t.sf(|t|, df)
  it("matches scipy", () => {
    closeTo(fn.studentP(2.5, 10), 0.03144684423660878, 1e-8, "t=2.5 df=10");
    closeTo(fn.studentP(0, 5), 1.0, 1e-12, "t=0");
    closeTo(fn.studentP(1.0, 1), 0.49999999999999956, 1e-8, "t=1 df=1");
  });

  it("returns 1 for a degenerate degrees of freedom", () => {
    assert.equal(fn.studentP(2.5, 0), 1);
  });
});

describe("betai", () => {
  it("is pinned at the boundaries", () => {
    assert.equal(fn.betai(1, 1, 0), 0);
    assert.equal(fn.betai(1, 1, 1), 1);
  });

  it("matches the regularized incomplete beta at an interior point", () => {
    // scipy.special.betainc(2, 3, 0.5)
    closeTo(fn.betai(2, 3, 0.5), 0.6875, 1e-12, "betai(2,3,0.5)");
  });
});

describe("bhFDR", () => {
  // statsmodels.stats.multitest.multipletests(..., method="fdr_bh")
  it("matches statsmodels", () => {
    const q = fn.bhFDR([0.01, 0.02, 0.5]);
    [0.03, 0.03, 0.5].forEach((e, i) => closeTo(q[i], e, 1e-12, `q[${i}]`));
  });

  it("returns q values in INPUT order, not sorted order", () => {
    const q = fn.bhFDR([0.5, 0.01, 0.02]);
    [0.5, 0.03, 0.03].forEach((e, i) => closeTo(q[i], e, 1e-12, `q[${i}]`));
  });

  it("matches statsmodels on a ten-hypothesis set", () => {
    const q = fn.bhFDR([0.001, 0.008, 0.039, 0.041, 0.042, 0.06, 0.074, 0.205, 0.212, 0.216]);
    const expected = [0.01, 0.04, 0.084, 0.084, 0.084, 0.1, 0.10571428571428572, 0.216, 0.216, 0.216];
    expected.forEach((e, i) => closeTo(q[i], e, 1e-12, `q[${i}]`));
  });

  it("stays monotone and never exceeds 1", () => {
    const ps = [0.9, 0.95, 0.99];
    const q = fn.bhFDR(ps);
    assert.ok(q.every((v) => v <= 1), "a q value exceeded 1");
    const sorted = ps.map((p, i) => [p, q[i]]).sort((a, b) => a[0] - b[0]).map((r) => r[1]);
    for (let i = 1; i < sorted.length; i++) {
      assert.ok(sorted[i] >= sorted[i - 1] - 1e-12, "q values are not monotone in p");
    }
  });

  it("returns an empty array for empty input", () => {
    assert.deepEqual(arr(fn.bhFDR([])), []);
  });
});

describe("pfmt", () => {
  it("formats p-values for display", () => {
    assert.equal(fn.pfmt(null), "n/a");
    assert.equal(fn.pfmt(1e-9), "< 1e-4");
    assert.equal(fn.pfmt(5e-4), "5.0e-4");
    assert.equal(fn.pfmt(0.2), "0.200");
  });
});
