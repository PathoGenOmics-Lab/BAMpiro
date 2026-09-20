"""A file that lives inside the container must never be declared as a Nextflow `path` input.

`path` means "stage this from the host". Nextflow honours it by adding a bind for the file's
parent directory, and Singularity then refuses to start a container whose bind source does not
exist on the cluster:

    FATAL: container creation failed: mount /opt/pathotypr->/opt/pathotypr error:
    while mounting /opt/pathotypr: mount source /opt/pathotypr doesn't exist

The three pathotypr assets and the canonical reference are bundled into the image at build time,
so the host has no copy of them and never will. Declared as `val` they are interpolated into the
script unchanged and resolved inside the container, where they do exist.

This cost a production run on garnatxa. Nothing caught it earlier because the stub tests replace
the script and the local Docker profile binds the work directory anyway, so the failure needs a
real Singularity task to appear. The check here is static, which is the only way it is cheap.

It derives the parameter list from `nextflow.config` rather than hard-coding it, so a new asset
bundled into the image is covered the day it is added.
"""

from __future__ import annotations

import re

import pytest

from conftest import REPO_ROOT

# `param = "/opt/..."`, which is the marker for "this path exists only inside the image".
IN_IMAGE = re.compile(r"^\s*(\w+)\s*=\s*[\"'](/opt/[^\"']+)[\"']", re.M)


def in_image_params():
    found = {}
    for cfg in [REPO_ROOT / "nextflow.config", *sorted((REPO_ROOT / "conf").glob("*.config"))]:
        for name, value in IN_IMAGE.findall(cfg.read_text()):
            found.setdefault(name, value)
    return found


def nf_sources():
    return sorted((REPO_ROOT / "modules").glob("*.nf")) + \
        sorted((REPO_ROOT / "subworkflows").glob("*.nf")) + [REPO_ROOT / "main.nf"]


def test_the_config_still_bundles_something():
    """The two tests below pass vacuously if the pattern stops matching, which it would on a
    rename of the parameters or a move of the assets. Fail loudly instead of quietly."""
    assert in_image_params(), "no in-image parameter found; has the /opt convention changed?"


@pytest.mark.parametrize("name", sorted(in_image_params()))
def test_no_in_image_param_is_wrapped_in_file(name):
    """`file(params.x)` is the other way to ask Nextflow to stage it, and it fails the same way."""
    wrapped = re.compile(rf"\bfile\s*\(\s*params\.{re.escape(name)}\b")
    for src in nf_sources():
        text = src.read_text()
        for n, line in enumerate(text.splitlines(), 1):
            assert not wrapped.search(line), \
                f"{src.relative_to(REPO_ROOT)}:{n} stages an in-image file from the host: {line.strip()}"


@pytest.mark.parametrize("name", sorted(in_image_params()))
def test_every_in_image_param_reaches_a_val_input(name):
    """Follow each call site to the process input at the same position and read its qualifier.

    Positional, because that is how Nextflow itself matches an argument to an input. A call
    spread over several lines with one argument on each is the shape used throughout; anything
    else is reported rather than guessed at, so a reformatting that defeats the parser shows up
    as a failure instead of as silence.
    """
    checked = 0
    for src in nf_sources():
        lines = src.read_text().splitlines()
        for n, line in enumerate(lines):
            if not re.search(rf"^\s*params\.{re.escape(name)}\s*,?\s*$", line):
                continue
            call = _enclosing_call(lines, n)
            assert call, f"{src.relative_to(REPO_ROOT)}:{n + 1} is not a recognisable call site"
            process, index = call
            qualifier = _input_qualifier(process, index)
            assert qualifier == "val", (
                f"{src.relative_to(REPO_ROOT)}:{n + 1} passes params.{name} as argument "
                f"{index + 1} of {process}, declared `{qualifier}`. It lives inside the image, "
                f"so it has to be `val`")
            checked += 1
    # Script interpolation (`${params.x}` inside a shell block) needs no input at all and is fine.
    assert checked >= 0


def _enclosing_call(lines, n):
    """(process name, zero-based argument index) for the call the line at `n` is an argument of."""
    for start in range(n - 1, max(-1, n - 40), -1):
        m = re.search(r"\b([A-Z][A-Z0-9_]*)\s*\(\s*$", lines[start])
        if m:
            args = [line for line in lines[start + 1:n + 1] if line.strip()]
            return m.group(1), len(args) - 1
        if ")" in lines[start]:
            return None
    return None


def _input_qualifier(process, index):
    """The qualifier of the `index`-th input of `process`, as written in modules/."""
    for src in (REPO_ROOT / "modules").glob("*.nf"):
        text = src.read_text()
        m = re.search(rf"^process\s+{re.escape(process)}\s*\{{", text, re.M)
        if not m:
            continue
        body = text[m.end():]
        block = re.search(r"^\s*input:\s*$(.*?)^\s*output:\s*$", body, re.M | re.S)
        assert block, f"{process} has no input block to read"
        decls = [line.strip() for line in block.group(1).splitlines() if line.strip()
                 and not line.strip().startswith("//")]
        assert index < len(decls), f"{process} takes {len(decls)} inputs, asked for {index + 1}"
        return decls[index].split("(")[0].split()[0]
    raise AssertionError(f"process {process} not found in modules/")
