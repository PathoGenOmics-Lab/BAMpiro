"""The pinned container digest is written down in more places than one, and they must agree.

`nextflow.config` is the only copy that does anything: it is what every run pulls. The others are
documentation, and a stale one there is worse than no mention at all, because it tells a reader
they are running an image they are not.

There were six copies when this was written, and updating them by grep found the last three only
after the first three were already changed. That is the argument for a test rather than a habit.
"""

from __future__ import annotations

import re
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
DIGEST = re.compile(r"sha256:([0-9a-f]{7,64})")


def _pinned():
    """The digest nextflow.config actually pins, which is the one that counts."""
    for line in (REPO / "nextflow.config").read_text().splitlines():
        if line.strip().startswith("container") and "sha256:" in line:
            return DIGEST.search(line).group(1)
    raise AssertionError("nextflow.config pins no container digest")


def _mentions():
    """Every other file that writes a digest down, and what it says. Abbreviated forms count."""
    out = {}
    for path in sorted(REPO.glob("docs/**/*.md")) + sorted(REPO.glob("docs/**/*.ipynb")) \
            + [REPO / "README.md"]:
        if not path.is_file():
            continue
        found = {m.group(1) for m in DIGEST.finditer(path.read_text(errors="replace"))}
        if found:
            out[path.relative_to(REPO).as_posix()] = found
    return out


def test_every_written_down_digest_matches_the_one_the_pipeline_pulls():
    pinned = _pinned()
    wrong = {}
    for where, found in _mentions().items():
        for d in found:
            # the docs abbreviate, so a prefix of the real digest is the match to make
            if not (pinned.startswith(d) or d.startswith(pinned)):
                wrong.setdefault(where, set()).add(d)

    assert not wrong, (
        f"nextflow.config pins sha256:{pinned[:12]}, but these disagree: "
        + "; ".join(f"{w} says {', '.join(sorted(d)[:2])}" for w, d in sorted(wrong.items())))


def test_the_documentation_mentions_the_digest_at_all():
    """Guards the guard: a test that finds nothing to compare passes for the wrong reason."""
    assert len(_mentions()) >= 3, "the digest used to appear in six files; found almost none"
