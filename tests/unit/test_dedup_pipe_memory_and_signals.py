"""The dedup pipe runs two sorts at once, and an OOM inside it does not look like an OOM.

`samtools sort -m` is PER THREAD, and MERGE_AND_MARKDUP has two sorts alive simultaneously:
a name sort feeding fixmate, and a coordinate sort feeding markdup. Sizing each at 80% of the
task's allocation asks for 160% of it, which at 8 GB over 4 threads was 1638M x 4 x 2 = 13.1 GB.

It hid for a long time because a single sequencing run produces a BAM too small to fill the
buffers. It surfaced on the first MERGED sample, whose four input BAMs gave a 2.8 GB stream:

    [E::bgzf_read_block] Failed to read BGZF block data at offset 2817752658
                         expected 64482 bytes; hread returned 48284
    samtools fixmate: Couldn't read from input file

The second half is worse than the first. The kernel kills one stage of a pipe, its stdout closes
mid-BGZF-block, and the NEXT stage reports a short read and exits 1. From outside, an
out-of-memory is indistinguishable from a corrupt file, so `errorStrategy` retried nothing and the
run died on something a retry with more memory would have fixed. Hence the PIPESTATUS guard: it
re-raises a signal death as 137 so the existing policy can see it.
"""

from __future__ import annotations

import re

import pytest

from conftest import REPO_ROOT

MODULE = REPO_ROOT / "modules" / "mapping.nf"


def dedup_block():
    t = MODULE.read_text()
    start = t.index("process MERGE_AND_MARKDUP")
    nxt = t.find("\nprocess ", start + 1)
    return t[start:nxt if nxt > 0 else len(t)]


def test_two_concurrent_sorts_fit_in_the_allocation():
    """The whole point: the factor has to leave room for BOTH sorts, not one."""
    block = dedup_block()
    sorts = re.findall(r"samtools sort[^|\n]*-m \$mem_per_thread", block)
    assert len(sorts) == 2, f"expected two sorts sharing the budget, found {len(sorts)}"

    m = re.search(r"mem_per_thread=.*?\* ([0-9.]+) / ", block)
    assert m, "the per-thread memory is no longer computed from a fraction of the allocation"
    factor = float(m.group(1))

    assert len(sorts) * factor <= 0.8, (
        f"{len(sorts)} sorts at {factor} of the allocation each asks for "
        f"{len(sorts) * factor:.0%} of it, and samtools -m is per thread")


@pytest.mark.parametrize("mb, cpus", [(8192, 4), (16384, 4), (4096, 2), (2048, 1)])
def test_the_arithmetic_leaves_headroom_at_every_size(mb, cpus):
    """Evaluated rather than eyeballed, because the failure is silent until a BAM is big enough."""
    block = dedup_block()
    factor = float(re.search(r"mem_per_thread=.*?\* ([0-9.]+) / ", block).group(1))
    floor = int(re.search(r"max\((\d+),", block).group(1))

    per_thread = max(floor, int(mb * factor / cpus))
    both_sorts = 2 * cpus * per_thread

    assert both_sorts < mb, f"{both_sorts} MB of sort buffers against a {mb} MB allocation"


def test_a_signal_death_inside_the_pipe_is_re_raised():
    """Without this the retry policy never fires: bash reports the short read, not the kill."""
    block = dedup_block()

    assert "PIPESTATUS" in block, "the pipe's per-stage exit codes are no longer captured"
    assert re.search(r'-ge 128', block), "a stage killed by a signal is no longer detected"
    assert re.search(r'exit "\$s"', block), "the captured status is never re-raised"


def test_the_retry_policy_still_covers_the_code_the_guard_raises():
    """The guard raises 128+signal. That is only useful if errorStrategy retries on it."""
    block = dedup_block()
    m = re.search(r"task\.exitStatus in \[([^\]]+)\]", block)

    assert m, "MERGE_AND_MARKDUP no longer declares which exit codes retry"
    codes = {int(c.strip()) for c in m.group(1).split(",")}
    assert 137 in codes, f"SIGKILL (137) is not retried; the guard would raise it for nothing: {codes}"
