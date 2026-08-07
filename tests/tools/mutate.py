"""Change the code on purpose and see whether the tests notice.

    python3 tests/tools/mutate.py [substring of a target file]

A passing suite says nothing about whether it would still pass on a broken version, which is the
only thing worth knowing about it. This flips one comparison, boundary or connective at a time and
runs the tests: a mutant that survives is a decision nothing checks.

It found, among others: a depletion running to the last site of a locus being dropped in silence,
--min-bf and --min-tract-af and --slack each movable by one, a tract nested inside a longer one
shortening the event that contains it, and a merged donor candidate forgetting the earlier of the
two starts it merged. None of it needed real data. All of it passed the suite.

Not a general mutation tester. It only touches lines that decide something, and never inside a
string, because a mutant in a message survives for reasons that say nothing and buries the ones
that do. Survivors still need reading: some are genuinely equivalent, and the ones over the
progress counters change a log line rather than a result.

Run from anywhere in the checkout. It works on a COPY of the tracked files, never on the working
tree, because a mutant sitting on disk poisons anything else that reads that file while it is
there, up to and including a backup taken during the run.
"""

import re
import shutil
import subprocess
import sys
import tempfile
import tokenize
from collections import defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
TARGETS = ["bin/gconv_model.py", "bin/gconv_cohort.py", "bin/gene_conversion.py",
           "bin/paralog_map.py"]
TESTS = ["tests/unit/test_gconv_model.py", "tests/unit/test_gconv_cohort.py",
         "tests/unit/test_gene_conversion.py", "tests/unit/test_paralog_map.py",
         "tests/unit/test_gconv_properties.py", "tests/unit/test_gene_conversion_chain.py"]

# (pattern, replacement). Applied to one occurrence on one line at a time.
RULES = [
    (r"(?<![<>=!])>=(?!=)", ">"), (r"(?<![<>=!])<=(?!=)", "<"),
    (r"(?<![<>=!])>(?![=])", ">="), (r"(?<![<>=!])<(?![=])", "<="),
    (r"(?<![<>=!])==(?!=)", "!="), (r"!=", "=="),
    (r"\band\b", "or"), (r"\bor\b", "and"),
    (r"\bnot\s+", ""),
    (r"\bmax\(", "min("), (r"\bmin\(", "max("),
]

IGNORE = {tokenize.NL, tokenize.NEWLINE, tokenize.INDENT, tokenize.DEDENT,
          tokenize.ENDMARKER, tokenize.ENCODING}
TEXT = {tokenize.STRING, getattr(tokenize, "FSTRING_MIDDLE", tokenize.STRING)}


def prose(path):
    """Line numbers carrying nothing but a string literal or a comment.

    Tokenised rather than matched, because the first attempt only recognised the line a docstring
    OPENS on and took every line of prose after it for code: this codebase explains its reasoning
    at length, so most of what came back was an `and` in a sentence.

    A line counts as prose only when EVERYTHING on it is prose. Excluding every line that merely
    contains a string was the opposite mistake, and it silently removed the comparisons that
    matter most here, since a verdict is a string and deciding about one means naming it.
    """
    lines = defaultdict(list)
    with open(path, "rb") as fh:
        for tok in tokenize.tokenize(fh.readline):
            if tok.type in IGNORE:
                continue
            for n in range(tok.start[0], tok.end[0] + 1):
                lines[n].append(tok.type)
    return {n for n, kinds in lines.items()
            if all(k in TEXT | {tokenize.COMMENT} for k in kinds)}


def literals(path):
    """(line, start, end) column spans of every string literal, to mutate around.

    A line can be code and prose at once. `raise ValueError("must be between 0 and 1")` is a
    decision, so it is kept, and then the `and` inside the message is flipped and survives for
    the obvious reason. Those survivors are noise that hides the real ones.
    """
    spans = defaultdict(list)
    with open(path, "rb") as fh:
        for tok in tokenize.tokenize(fh.readline):
            # FSTRING_MIDDLE, because an f-string is not one STRING token here: it is split
            # around its {} and the prose between them is what needs protecting.
            if (tok.type in TEXT and tok.start[0] == tok.end[0]):
                spans[tok.start[0]].append((tok.start[1], tok.end[1]))
    return spans


def decisions(path):
    """Lines that decide something, with the mutations to try on each."""
    out = []
    text, inside = prose(path), literals(path)
    for n, line in enumerate(path.read_text().splitlines(), 1):
        if n in text or not line.strip():
            continue
        code = line.split("#")[0]
        if not re.search(r"\b(if|while|elif|assert|return|and|or|max|min)\b|[<>=!]=|[<>]", code):
            continue
        for pat, rep in RULES:
            for m in re.finditer(pat, code):
                if any(a <= m.start() < b for a, b in inside.get(n, ())):
                    continue
                out.append((n, m.start(), m.end(), rep, line))
    return out


def checkout(into):
    """Every tracked file, copied. Uncommitted work included, so what runs is what you have."""
    listed = subprocess.run(["git", "-C", str(REPO), "ls-files", "-z"],
                            capture_output=True, text=True, check=True)
    for name in listed.stdout.split("\0"):
        if not name:
            continue
        src, dst = REPO / name, into / name
        if src.is_file():
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src, dst)
    return into


def run(root, only):
    survivors, killed, total = [], 0, 0
    for target in TARGETS:
        if only and only not in target:
            continue
        path = root / target
        original = path.read_text()
        backup = original
        muts = decisions(path)
        print(f"{target}: {len(muts)} mutants", flush=True)
        try:
            for n, a, b, rep, line in muts:
                lines = original.splitlines(keepends=True)
                raw = lines[n - 1]
                lines[n - 1] = raw[:a] + rep + raw[b:]
                path.write_text("".join(lines))
                try:
                    compile(path.read_text(), str(path), "exec")
                except SyntaxError:
                    continue
                total += 1
                r = subprocess.run([sys.executable, "-m", "pytest", *TESTS, "-x", "-q",
                                    "--no-header", "-p", "no:cacheprovider"],
                                   cwd=root, capture_output=True, text=True, timeout=600)
                if r.returncode == 0:
                    survivors.append((target, n, line.strip(), raw[a:b], rep))
                    print(f"  SURVIVED {target}:{n}  {raw[a:b]!r} -> {rep!r}\n"
                          f"           {line.strip()[:100]}", flush=True)
                else:
                    killed += 1
        finally:
            path.write_text(backup)

    print(f"\n{killed}/{total} killed, {len(survivors)} survived")
    return survivors


def main():
    only = sys.argv[1] if len(sys.argv) > 1 else None
    work = Path(tempfile.mkdtemp(prefix="bampiro-mutants-"))
    try:
        run(checkout(work), only)
    finally:
        shutil.rmtree(work, ignore_errors=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
