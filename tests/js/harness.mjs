// Loads the report's front-end so its pure functions can be unit-tested from Node.
//
// The report ships as ONE ES5 IIFE: bin/qcreport/render.py concatenates bin/report_assets/js/*.js
// in a fixed order and drops the result into shell.html. 01_prelude.js opens the IIFE and
// 14_boot.js closes it, so no individual fragment (and no prefix of them) is a valid
// script on its own.
//
// To reach the functions without editing the sources we concatenate everything except
// 14_boot.js - which is pure DOM wiring - and append our own epilogue that captures the
// declarations and closes the IIFE. Function declarations hoist to the top of the IIFE,
// so the epilogue sees every function regardless of which fragment declared it. Nothing
// here depends on line numbers, so the tests survive refactoring of the JS.

import { readFileSync } from "node:fs";
import { dirname, join } from "node:path";
import { fileURLToPath } from "node:url";
import vm from "node:vm";

const ROOT = join(dirname(fileURLToPath(import.meta.url)), "..", "..");
const JS_DIR = join(ROOT, "bin", "report_assets", "js");

/** The module order the report uses, read from bin/qcreport/render.py itself so the two cannot drift. */
export function moduleOrder() {
  const py = readFileSync(join(ROOT, "bin", "qcreport", "render.py"), "utf8");
  const block = py.match(/_JS_MODULES\s*=\s*\[([\s\S]*?)\]/);
  if (!block) throw new Error("could not find _JS_MODULES in bin/qcreport/render.py");
  return [...block[1].matchAll(/"js\/([^"]+)"/g)].map((m) => m[1]);
}

export function readModule(name) {
  return readFileSync(join(JS_DIR, name), "utf8");
}

/** The exact string render.py injects into shell.html (a bare concatenation, no separator). */
export function bundle() {
  return moduleOrder().map(readModule).join("");
}

// Minimal stand-ins for what the fragments touch while they are being evaluated:
// 01_prelude.js reads the REPORT global and calls isDark(); 03_state.js reads
// window.innerWidth and iterates R.samples / R.metrics.
function stubContext(report) {
  const classList = { contains: () => false, add() {}, remove() {}, toggle() {} };
  const element = {
    classList,
    style: {},
    textContent: "",
    innerHTML: "",
    appendChild() {},
    setAttribute() {},
    getAttribute: () => null,
    addEventListener() {},
    querySelectorAll: () => [],
  };
  const document = {
    documentElement: { classList },
    body: { classList },
    getElementById: () => null,
    querySelector: () => null,
    querySelectorAll: () => [],
    createElement: () => element,
    createElementNS: () => element,
    addEventListener() {},
  };
  const window = {
    innerWidth: 1280,
    addEventListener() {},
    matchMedia: () => ({ matches: false, addEventListener() {} }),
    localStorage: { getItem: () => null, setItem() {}, removeItem() {} },
  };
  const ctx = {
    REPORT: report,
    document,
    window,
    console,
    Math,
    JSON,
    setTimeout,
    clearTimeout,
    requestAnimationFrame: (fn) => fn(),
  };
  ctx.globalThis = ctx;
  return ctx;
}

/**
 * A REPORT payload just rich enough for every fragment to finish evaluating.
 *
 * The fragments run real setup code at load: 01_prelude.js derives metrics off R.samples and
 * pushes onto R.dist, 03_state.js builds the threshold and metric tables. An empty cohort keeps
 * all of that a no-op while still exercising the code path.
 */
export function emptyReport() {
  return {
    samples: [],
    metrics: [],
    dist: [],
    lineages: [],
    gene_map: {},
    defs: {},
    n_ancient: 0,
    nbins: 0,
    genome_len: 0,
    snp_density_ok: false,
    mask_bins: [],
    thresholds: {},
    version: "test",
  };
}

/**
 * Evaluate every fragment except the DOM-wiring one and return the named functions.
 *
 * @param {string[]} names function names to pull out of the IIFE scope
 * @param {object} [report] optional REPORT payload override
 */
export function loadFunctions(names, report = emptyReport()) {
  const order = moduleOrder().filter((m) => m !== "14_boot.js");
  const epilogue =
    "\n__CAPTURE__(" +
    "{" +
    names.map((n) => `${JSON.stringify(n)}: typeof ${n} === "undefined" ? undefined : ${n}`).join(",") +
    "});\n})();\n";

  const captured = {};
  const ctx = stubContext(report);
  ctx.__CAPTURE__ = (obj) => Object.assign(captured, obj);

  vm.createContext(ctx);
  new vm.Script(order.map(readModule).join("") + epilogue, { filename: "bampiro-report-bundle.js" }).runInContext(ctx);

  const missing = names.filter((n) => typeof captured[n] !== "function");
  if (missing.length) throw new Error(`not found in the report bundle: ${missing.join(", ")}`);
  return captured;
}
