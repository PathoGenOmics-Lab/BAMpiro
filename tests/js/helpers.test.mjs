// Pure helpers behind the report's interactive panels: the column filter grammar, the SNP-dynamics
// classifier, the QC-space linear algebra, and the small formatters. None of these touch the DOM.

import { strict as assert } from "node:assert";
import { before, describe, it } from "node:test";
import { loadFunctions } from "./harness.mjs";

const arr = (x) => Array.from(x);

function closeTo(actual, expected, tol, what) {
  assert.ok(
    Math.abs(actual - expected) <= tol,
    `${what}: ${actual} is not within ${tol} of ${expected}`,
  );
}

let fn;
before(() => {
  fn = loadFunctions([
    "colMatchOne", "esc", "fmt", "shortv", "fmtpos",
    "quant", "med2", "robustMedSig", "robustMS",
    "dynLogit", "dynSlope", "dynSelCls",
    "covMatrix", "topEigs", "invRidge", "chi2q975", "spreadEllipse", "yearOf",
    "epiFmtR", "epiNodeLbl", "drStatus", "drGColor", "snpAfColor", "zeros", "assign",
  ]);
});

describe("colMatchOne, the per-column filter grammar", () => {
  it("treats a blank filter as matching everything", () => {
    assert.equal(fn.colMatchOne("", 5, "anything"), true);
    assert.equal(fn.colMatchOne("   ", null, ""), true);
  });

  it("applies numeric comparisons", () => {
    assert.equal(fn.colMatchOne(">=10", 10, ""), true);
    assert.equal(fn.colMatchOne(">=10", 9.9, ""), false);
    assert.equal(fn.colMatchOne("<5", 4, ""), true);
    assert.equal(fn.colMatchOne("<5", 5, ""), false);
    assert.equal(fn.colMatchOne("=42", 42, ""), true);
    assert.equal(fn.colMatchOne(">-3", -2, ""), true);
  });

  it("applies ranges and tolerates a reversed one", () => {
    assert.equal(fn.colMatchOne("10..20", 15, ""), true);
    assert.equal(fn.colMatchOne("10..20", 21, ""), false);
    assert.equal(fn.colMatchOne("20..10", 15, ""), true, "a reversed range should be swapped");
    assert.equal(fn.colMatchOne("1 to 3", 2, ""), true);
  });

  it("falls back to a case-insensitive substring match on the rendered text", () => {
    assert.equal(fn.colMatchOne("kat", null, "katG"), true);
    assert.equal(fn.colMatchOne("KAT", null, "katG"), true);
    assert.equal(fn.colMatchOne("rpo", null, "katG"), false);
  });

  it("does not match a numeric filter against a null value", () => {
    assert.equal(fn.colMatchOne(">=10", null, "10"), false);
  });
});

describe("formatters", () => {
  it("escapes the HTML metacharacters it claims to", () => {
    assert.equal(fn.esc('<a href="x">&'), "&lt;a href=&quot;x&quot;&gt;&amp;");
  });

  it("renders NA rather than a crash for a missing value", () => {
    assert.equal(fn.fmt(null, "int"), "NA");
    assert.equal(fn.fmt(null, "pct"), "NA");
    assert.equal(fn.shortv(null, "int"), "");
  });

  it("formats by metric kind", () => {
    assert.equal(fn.fmt(1234.6, "int"), "1,235");
    assert.equal(fn.fmt(12.345, "pct"), "12.3");
    assert.equal(fn.fmt(12.345, "float"), "12.35");
  });

  it("scales genome coordinates", () => {
    assert.equal(fn.fmtpos(4411532), "4.41 Mb");
    assert.equal(fn.fmtpos(15000), "15 kb");
    assert.equal(fn.fmtpos(742), "742 bp");
  });
});

describe("quantiles and robust spread", () => {
  it("interpolates quantiles on a pre-sorted array", () => {
    const s = [1, 2, 3, 4, 5];
    assert.equal(fn.quant(s, 0), 1);
    assert.equal(fn.quant(s, 1), 5);
    assert.equal(fn.quant(s, 0.5), 3);
    closeTo(fn.quant(s, 0.25), 2, 1e-12, "q25");
    assert.equal(fn.quant([], 0.5), null);
  });

  it("averages the two middle values for an even-length median", () => {
    assert.equal(fn.med2([1, 2, 3, 4]), 2.5);
    assert.equal(fn.med2([1, 2, 3]), 2);
  });

  it("needs four non-null values before reporting a robust spread", () => {
    assert.equal(fn.robustMedSig([1, 2, 3]), null);
    assert.equal(fn.robustMedSig([1, null, 2, null]), null);
    const r = fn.robustMedSig([10, 12, 14, 16]);
    assert.equal(r.med, 13);
    assert.ok(r.sig > 0, "sigma must be positive");
  });

  it("never returns a zero sigma, which would divide by zero downstream", () => {
    // robustMS has no minimum-n and floors sigma at 1.
    const [med, sig] = arr(fn.robustMS([7, 7, 7, 7]));
    assert.equal(med, 7);
    assert.ok(sig > 0, `sigma was ${sig}`);
  });
});

describe("SNP dynamics", () => {
  it("clamps the logit so a fixed or absent allele stays finite", () => {
    assert.ok(Number.isFinite(fn.dynLogit(0)), "AF 0 produced a non-finite logit");
    assert.ok(Number.isFinite(fn.dynLogit(1)), "AF 1 produced a non-finite logit");
    assert.equal(fn.dynLogit(0.5), 0);
  });

  it("fits a positive slope to a rising trajectory", () => {
    const r = fn.dynSlope([0.05, 0.3, 0.7, 0.95], [0, 1, 2, 3]);
    assert.ok(r.s > 0, `slope was ${r.s}`);
    assert.ok(r.r2 > 0.9, `r2 was ${r.r2}`);
    assert.equal(r.n, 4);
  });

  it("returns null with fewer than two observations", () => {
    assert.equal(fn.dynSlope([null, 0.4, null], [0, 1, 2]), null);
    assert.equal(fn.dynSlope([], []), null);
  });

  it("reports a zero slope when every timepoint is the same", () => {
    const r = fn.dynSlope([0.2, 0.4, 0.6], [5, 5, 5]);
    assert.equal(r.s, 0);
  });

  it("classifies the selection patterns the panel colours by", () => {
    assert.equal(fn.dynSelCls([0.05, 0.5, 0.95]).cls, "sweep");
    assert.equal(fn.dynSelCls([0.02, 0.2, 0.4]).cls, "emerge");
    assert.equal(fn.dynSelCls([0.9, 0.5, 0.05]).cls, "lost");
    assert.equal(fn.dynSelCls([0.5, 0.52, 0.49]).cls, "stable");
    assert.equal(fn.dynSelCls([0.4]).cls, "single");
  });
});

describe("QC-space linear algebra", () => {
  it("computes a covariance matrix with the n-1 divisor", () => {
    // covMatrix does NOT re-center: buildZ has already z-scored the columns, so the input must
    // be centered. Columns here are [-1,0,1] and [-2,0,2]; var = 1 and 4, cov = 2.
    const C = fn.covMatrix([[-1, -2], [0, 0], [1, 2]], 2);
    closeTo(C[0][0], 1, 1e-12, "var x");
    closeTo(C[1][1], 4, 1e-12, "var y");
    closeTo(C[0][1], 2, 1e-12, "cov");
    closeTo(C[0][1], C[1][0], 1e-15, "symmetry");
  });

  it("recovers the eigenvalues of a diagonal matrix", () => {
    const e = fn.topEigs([[9, 0], [0, 4]], 2, 2);
    closeTo(e.vals[0], 9, 1e-6, "first eigenvalue");
    closeTo(e.vals[1], 4, 1e-6, "second eigenvalue");
  });

  it("is deterministic, since power iteration is seeded rather than random", () => {
    const C = [[5, 2], [2, 3]];
    const a = fn.topEigs(C, 2, 2);
    const b = fn.topEigs(C, 2, 2);
    assert.deepEqual(arr(a.vals), arr(b.vals));
  });

  it("inverts a ridge-regularized matrix", () => {
    const inv = fn.invRidge([[2, 0], [0, 4]], 2, 0);
    closeTo(inv[0][0], 0.5, 1e-12, "inv[0][0]");
    closeTo(inv[1][1], 0.25, 1e-12, "inv[1][1]");
  });

  it("returns null instead of dividing by a singular pivot", () => {
    assert.equal(fn.invRidge([[0, 0], [0, 0]], 2, 0), null);
  });

  it("approximates the chi-square 97.5th percentile", () => {
    // Wilson-Hilferty, so an approximation. Reference: scipy.stats.chi2.ppf(0.975, k).
    // The tolerance is the accuracy the approximation actually achieves, which improves with k.
    closeTo(fn.chi2q975(2), 7.3778, 0.05, "k=2");
    closeTo(fn.chi2q975(4), 11.1433, 0.03, "k=4");
    closeTo(fn.chi2q975(10), 20.4832, 0.02, "k=10");
  });

  it("fits a spread ellipse and needs at least four points", () => {
    assert.equal(fn.spreadEllipse([[0, 0], [1, 1], [2, 2]]), null);
    const e = fn.spreadEllipse([[0, 0], [1, 0], [0, 1], [1, 1], [2, 2]]);
    assert.ok(e.rx >= 0 && e.ry >= 0, "ellipse radii must be non-negative");
  });

  it("extracts a plausible year and rejects an implausible one", () => {
    assert.equal(fn.yearOf("collected 2018-04-02"), 2018);
    assert.equal(fn.yearOf("1998"), 1998);
    assert.equal(fn.yearOf("sample 12"), null);
    assert.equal(fn.yearOf(null), null);
  });
});

describe("panel colouring and labels", () => {
  it("marks resistance from the worst WHO grade present", () => {
    assert.equal(fn.drStatus([1, 4]).t, "R");
    assert.equal(fn.drStatus([3, 5]).t, "?");
    assert.equal(fn.drStatus([5]).t, "&#183;");
    assert.equal(fn.drStatus([]).t, "");
  });

  it("colours WHO grades 1 and 2 as resistant", () => {
    assert.equal(fn.drGColor(1), fn.drGColor(2));
    assert.notEqual(fn.drGColor(1), fn.drGColor(3));
  });

  it("formats a correlation with an explicit sign and no leading zero", () => {
    assert.equal(fn.epiFmtR(-0.85), "-.85");
    assert.equal(fn.epiFmtR(0.85), "+.85");
    assert.equal(fn.epiFmtR(0), ".00");
  });

  it("labels an intergenic epistasis node rather than leaving it blank", () => {
    assert.equal(fn.epiNodeLbl({ gene: "", pos: "chr:761155" }), "(intergenic) 761155");
    assert.equal(fn.epiNodeLbl({ gene: "rpoB", pos: 761155 }), "rpoB 761155");
  });

  it("scales the SNP-matrix cell opacity with allele frequency", () => {
    assert.notEqual(fn.snpAfColor(0), fn.snpAfColor(1));
  });
});

describe("ES5 shims", () => {
  it("zeros builds a zero-filled array", () => {
    assert.deepEqual(arr(fn.zeros(3)), [0, 0, 0]);
    assert.deepEqual(arr(fn.zeros(0)), []);
  });

  it("assign copies only own properties", () => {
    const target = fn.assign({ a: 1 }, { b: 2 });
    assert.equal(target.a, 1);
    assert.equal(target.b, 2);
    assert.equal(target.toString, Object.prototype.toString, "inherited keys must not be copied");
  });
});
