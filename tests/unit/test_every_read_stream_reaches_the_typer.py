"""The typer must see the same reads as the mapper, and must actually read the ones it is handed.

FASTP_PE is run with `--merge`, so a pair whose mates overlap leaves as one merged fragment in the
orphan stream instead of as r1/r2. How much of a library ends up there is a property of the insert
size, not of the pipeline: in the cohort that surfaced this it was 81-92 percent, and a fully
overlapping library puts ALL of it there.

`MAP_READS` took all three streams from the start. `LINEAGE_TYPING` took two, so pathotypr typed
whatever happened to be left unmerged. Nothing failed, no count was zero in any log, and a sample
typed on nothing at all would have reported `Unclassified` exactly as a genuinely unclassifiable
sample does. Measured on a purpose-built fully overlapping library: r1 and r2 reached the typer
holding 0 reads each while the merged stream held 3,000.

Two checks, because the bug had two halves and either alone can come back. The first is that the
channels are wired: whatever read streams reach the mapper reach the typer. The second is that the
process then uses them, since an input a script never names is staged and ignored in silence.

Static on purpose. Proving it by running needs a real container and a library whose mates overlap,
which is not what the suite this lives in is for.
"""

from __future__ import annotations

import re

import pytest

from conftest import REPO_ROOT


def call_args(text, process):
    """The argument list of a call, flattened, for a call written over several lines."""
    m = re.search(rf"\b{re.escape(process)}\s*\(", text)
    if not m:
        return None
    depth, out, i = 0, [], m.end() - 1
    while i < len(text):
        if text[i] == "(":
            depth += 1
        elif text[i] == ")":
            depth -= 1
            if depth == 0:
                break
        out.append(text[i])
        i += 1
    return "".join(out[1:])


def streams_passed_to(process):
    """The `reads.<name>` channels a call in main.nf receives."""
    args = call_args((REPO_ROOT / "main.nf").read_text(), process)
    assert args is not None, f"no call to {process} in main.nf"
    return set(re.findall(r"\breads\.(\w+)", args))


def test_the_mapper_still_takes_more_than_one_read_stream():
    """The comparison below is vacuous if MAP_READS stops taking several streams, which is exactly
    the refactor that would hide a regression rather than cause one."""
    assert len(streams_passed_to("MAP_READS")) >= 2


def test_the_typer_sees_every_read_stream_the_mapper_sees():
    mapper, typer = streams_passed_to("MAP_READS"), streams_passed_to("LINEAGE_TYPING")

    missing = mapper - typer
    assert not missing, (
        "LINEAGE_TYPING is not given " + ", ".join(f"reads.{s}" for s in sorted(missing))
        + ", so pathotypr types a subset of the reads that were mapped. FastP --merge puts most "
          "of an overlapping library in reads.pe_orphans")


def test_the_pathotypr_process_reads_every_file_it_is_given():
    """Wiring a channel in is half of it. The script has to name the file too."""
    text = (REPO_ROOT / "modules" / "pathotypr.nf").read_text()
    body = text[text.index("process RUN_PATHOTYPR_PE"):text.index("process RUN_PATHOTYPR_SE")]
    inputs = re.search(r"^\s*input:\s*$(.*?)^\s*output:\s*$", body, re.M | re.S)
    assert inputs, "RUN_PATHOTYPR_PE has no input block"
    reads = re.findall(r"path\((\w+)\)", inputs.group(1))
    assert len(reads) >= 3, f"expected the pair and the orphan stream, found {reads}"

    script = body[body.index("script:"):]
    cat = [line for line in script.splitlines() if line.strip().startswith("cat ")]
    assert cat, "RUN_PATHOTYPR_PE no longer combines its inputs before typing"

    for name in reads:
        assert f"${{{name}}}" in cat[0], (
            f"input `{name}` is staged but never reaches pathotypr: {cat[0].strip()}")


@pytest.mark.parametrize("flag", ["--no-auto-paired"])
def test_the_combined_file_is_not_mistaken_for_half_of_a_pair(flag):
    """Auto-detection keys on _R1 / _1 in the file name, and the combined file is named after the
    sample. A sampleId containing that pattern would otherwise be read as one mate of a pair."""
    text = (REPO_ROOT / "modules" / "pathotypr.nf").read_text()
    body = text[text.index("process RUN_PATHOTYPR_PE"):text.index("process RUN_PATHOTYPR_SE")]
    calls = [line for line in body.splitlines() if "split-fastq" in line]

    assert calls, "no pathotypr invocation in RUN_PATHOTYPR_PE"
    for line in calls:
        assert flag in line, f"{flag} missing from: {line.strip()}"
