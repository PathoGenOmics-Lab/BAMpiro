// Integrity of the report bundle itself.
//
// bin/qcreport/render.py builds the report by concatenating a HARD-CODED list of CSS and JS fragments
// and splicing them into shell.html. Nothing checks that list against what is actually on disk, so a
// new fragment that nobody registers ships as dead code and a syntax error in any one fragment
// only surfaces when a human opens the HTML. These tests close both gaps.

import { strict as assert } from "node:assert";
import { readdirSync, readFileSync } from "node:fs";
import { dirname, join } from "node:path";
import { fileURLToPath } from "node:url";
import { describe, it } from "node:test";
import vm from "node:vm";
import { bundle, moduleOrder, readModule } from "./harness.mjs";

const ROOT = join(dirname(fileURLToPath(import.meta.url)), "..", "..");
const ASSETS = join(ROOT, "bin", "report_assets");
const RENDER_PY = readFileSync(join(ROOT, "bin", "qcreport", "render.py"), "utf8");

function pyList(name) {
  const block = RENDER_PY.match(new RegExp(`${name}\\s*=\\s*\\[([\\s\\S]*?)\\]`));
  assert.ok(block, `could not find ${name} in bin/qcreport/render.py`);
  return [...block[1].matchAll(/"([^"]+)"/g)].map((m) => m[1]);
}

describe("asset registration", () => {
  it("ships every .js fragment on disk", () => {
    const onDisk = readdirSync(join(ASSETS, "js")).filter((f) => f.endsWith(".js")).sort();
    const registered = moduleOrder().slice().sort();
    assert.deepEqual(
      registered,
      onDisk,
      "bin/report_assets/js/ and _JS_MODULES in bin/qcreport/render.py disagree; " +
        "an unregistered fragment is never included in the report",
    );
  });

  it("ships every .css fragment on disk", () => {
    const onDisk = readdirSync(join(ASSETS, "css")).filter((f) => f.endsWith(".css")).sort();
    const registered = pyList("_CSS_MODULES").map((p) => p.replace(/^css\//, "")).sort();
    assert.deepEqual(registered, onDisk, "css/ and _CSS_MODULES disagree");
  });

  it("concatenates the JS fragments in numeric order", () => {
    const order = moduleOrder();
    const sorted = order.slice().sort();
    assert.deepEqual(order, sorted, "the fragments are numbered, so the list must stay sorted");
  });

  it("ends every fragment with a newline, since they are joined with no separator", () => {
    // render.py does "".join(...), so a fragment without a trailing newline would splice its
    // last line into the next fragment's first line.
    for (const name of moduleOrder()) {
      const text = readModule(name);
      assert.ok(text.endsWith("\n"), `${name} does not end with a newline`);
    }
  });
});

describe("bundle syntax", () => {
  it("parses as the single IIFE the report expects", () => {
    // 01_prelude.js opens the IIFE and 14_boot.js closes it, so only the whole bundle balances.
    assert.doesNotThrow(
      () => new vm.Script(bundle(), { filename: "bundle.js" }),
      "the concatenated report JS is not syntactically valid",
    );
  });

  it("parses inside the function wrapper shell.html puts it in", () => {
    const wrapped = `function __qcmain__(){\n${bundle()}\n}`;
    assert.doesNotThrow(() => new vm.Script(wrapped, { filename: "wrapped.js" }));
  });

  it("stays ES5, with no syntax an old browser would reject", () => {
    // The report is deliberately ES5: it has to open from a file:// URL on whatever browser a
    // collaborator happens to have, and past releases have shipped ES5 slips. Node parses ES2015+
    // happily, so detect the common slips by pattern. Only patterns that cannot plausibly occur in
    // prose are checked unanchored; `const`/`let` are anchored to a statement position, because
    // "constant" and words like "first-class" appear all over the comments.
    const offenders = [
      [/=>/g, "an arrow function"],
      [/`/g, "a template literal"],
      [/(?:^|[;{}()])\s*(?:const|let)\s+[A-Za-z_$]/gm, "a const/let declaration"],
    ];
    const problems = [];
    for (const name of moduleOrder()) {
      const text = readModule(name);
      for (const [re, label] of offenders) {
        for (const hit of text.matchAll(re)) {
          const line = text.slice(0, hit.index).split("\n").length;
          problems.push(`${name}:${line} uses ${label}`);
        }
      }
    }
    assert.deepEqual(problems, [], `non-ES5 syntax in the report bundle:\n  ${problems.join("\n  ")}`);
  });
});

describe("shell.html", () => {
  const shell = readFileSync(join(ASSETS, "shell.html"), "utf8");

  it("contains every placeholder render.py substitutes", () => {
    for (const token of ["__CSS__", "__JS__", "__JSON_GZ__", "__LOGO__", "__REPO__", "__TITLE__"]) {
      assert.ok(shell.includes(token), `shell.html is missing the ${token} placeholder`);
    }
  });

  it("does not itself contain a placeholder token inside the JS or CSS fragments", () => {
    // A fragment containing "__JS__" would be re-substituted and duplicate the whole bundle.
    for (const name of moduleOrder()) {
      const text = readModule(name);
      for (const token of ["__CSS__", "__JS__", "__JSON_GZ__", "__LOGO__"]) {
        assert.ok(!text.includes(token), `${name} contains the substitution token ${token}`);
      }
    }
  });
});
