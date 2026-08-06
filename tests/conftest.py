"""Shared fixtures for the BAMpiro test suite.

The pipeline's Python lives in `bin/` as standalone scripts, not as an installed
package, so tests load them by path. Every script guards its entry point behind
`if __name__ == "__main__"`, so importing one is side-effect free apart from
`qc_report`, which reads `bin/report_assets/` at import time.

The one exception is `bin/qcreport/`, the package behind `bin/qc_report.py`: it is
a real package and its tests import it normally (`from qcreport import parsers`),
which works because pyproject.toml puts `bin` on the pytest pythonpath. Loading
`qc_report` itself by path still works too, for the same reason.
"""

from __future__ import annotations

import importlib.util
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
BIN = REPO_ROOT / "bin"
DATA = Path(__file__).resolve().parent / "data"


def load_script(name: str):
    """Import a script from bin/ by name (without the .py extension).

    Registering the module in sys.modules before executing it is required, not
    cosmetic: WGS_fasta_allpos defines a @dataclass, and dataclasses resolve the
    module through sys.modules while the class body is being built.
    """
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, BIN / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="session")
def repo_root() -> Path:
    return REPO_ROOT


@pytest.fixture(scope="session")
def data_dir() -> Path:
    """The committed fixture cohort (see tests/data/make_test_data.py)."""
    return DATA


@pytest.fixture(scope="session")
def nextflow() -> str:
    """Path to the nextflow CLI, skipping the test when it is not installed."""
    exe = shutil.which("nextflow")
    if not exe:
        pytest.skip("nextflow is not on PATH")
    return exe


@pytest.fixture(scope="session")
def nextflow_env() -> dict:
    """Environment for invoking Nextflow from a test.

    No parser pin: main.nf is accepted by the strict (v2) parser that Nextflow
    25.10 and later use by default, as well as by the v1 parser of the minimum
    supported 24.04.2. CI runs both.
    """
    env = dict(os.environ)
    env["NXF_ANSI_LOG"] = "false"
    return env


def run_pipeline(nextflow_exe, env, workdir: Path, *args, timeout=600):
    """Run `nextflow run main.nf ...` from the repository root and capture the result."""
    cmd = [
        nextflow_exe,
        "-log", str(workdir / "nextflow.log"),
        "run", "main.nf",
        "-work-dir", str(workdir / "work"),
        *args,
    ]
    return subprocess.run(
        cmd,
        cwd=REPO_ROOT,
        env=env,
        capture_output=True,
        text=True,
        timeout=timeout,
    )
