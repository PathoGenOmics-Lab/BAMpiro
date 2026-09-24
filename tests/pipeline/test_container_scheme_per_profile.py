"""The pinned container reference carries a scheme that only one of the two runtimes accepts.

Singularity wants `docker://`. Docker refuses it: the run dies with `invalid reference format`
before it pulls anything, which is a confusing way to be told the image is fine and the prefix is
not. The default keeps the scheme, because that is the path every cluster run takes, and the
`docker` profile is the one that strips it.

Both halves are asserted here. Dropping the scheme from the default would fix a laptop and break
every cluster, which is the more expensive direction to get wrong. get_MNV runs in an image of its
own (`--mnv_container`), and the same two runtimes read that reference too.
"""

from __future__ import annotations

import re
import subprocess

import pytest

from conftest import REPO_ROOT

pytestmark = pytest.mark.nextflow

CONTAINER = re.compile(r"^\s*container\s*=\s*'([^']+)'", re.M)


def effective_container(nextflow, env, profile, image="bampiro"):
    """What a process's container resolves to under a profile, as Nextflow itself reports it."""
    out = subprocess.run(
        [nextflow, "config", "-profile", profile, str(REPO_ROOT)],
        capture_output=True, text=True, cwd=REPO_ROOT, env=env, timeout=300)
    assert out.returncode == 0, out.stderr[-2000:]
    found = [m.group(1) for m in CONTAINER.finditer(out.stdout) if image in m.group(1)]
    assert found, f"no container resolved under -profile {profile}"
    return found[-1]


@pytest.mark.parametrize("profile", ["standard", "local", "slurm", "garnatxa"])
def test_every_singularity_profile_keeps_the_scheme(nextflow, nextflow_env, profile):
    """Singularity needs it, and these are the profiles a cluster run uses."""
    got = effective_container(nextflow, nextflow_env, profile)

    assert got.startswith("docker://"), f"-profile {profile} resolved {got}"


def test_the_docker_profile_drops_the_scheme(nextflow, nextflow_env):
    """Docker rejects it outright, so the one profile that uses Docker has to remove it."""
    got = effective_container(nextflow, nextflow_env, "local,docker")

    assert not got.startswith("docker://"), got
    assert "bampiro@sha256:" in got, "the digest pin has to survive the rewrite"


def test_the_two_forms_are_the_same_image(nextflow, nextflow_env):
    """The point of stripping a scheme is that nothing else changes with it. A profile that
    quietly resolved a different tag would pass both tests above and still be wrong."""
    cluster = effective_container(nextflow, nextflow_env, "garnatxa")
    laptop = effective_container(nextflow, nextflow_env, "local,docker")

    assert cluster == "docker://" + laptop


def test_get_mnv_image_takes_the_scheme_each_runtime_wants(nextflow, nextflow_env):
    """GET_MNV's container is set on its own, so the docker profile has to strip it on its own."""
    cluster = effective_container(nextflow, nextflow_env, "garnatxa", image="get_mnv")
    laptop = effective_container(nextflow, nextflow_env, "local,docker", image="get_mnv")

    assert cluster.startswith("docker://"), cluster
    assert "get_mnv@sha256:" in laptop, "the digest pin has to survive the rewrite"
    assert cluster == "docker://" + laptop
