"""WF-13 step 1: the installed package resolves outside the repository.

Brief: docs/planning/chunks/WF-13.md, §9 item 220. ``scripts/package/wf13_setup.sh``
installs the package with ``pip install .`` into a fresh venv under /tmp/wf13;
this module, run from that venv with cwd outside /workspace and no PYTHONPATH,
asserts the package came from the venv and that the requested DolfinX build is
the one loaded (the ``sitecustomize`` trap: the wrong build would otherwise come
back silently).

``test_installed_package_is_outside_the_repo`` is skipped unless
``WF13_INSTALLED_CHECK=1`` so ordinary ``pytest tests`` runs from the repo never
go red on it. ``test_control_repo_location_is_detected`` always runs: it is the
negative control (c1) showing the checker flags a repo-resolved package.
"""

from __future__ import annotations

import os
import sys

import pytest

REPO = "/workspace"


def _under(path: str, root: str) -> bool:
    path = os.path.abspath(path)
    root = os.path.abspath(root)
    return path == root or path.startswith(root.rstrip("/") + "/")


def _location_violations(pkg_file, sys_path, environ, cwd) -> list[str]:
    """Return a named violation for each way the run could resolve to the repo."""
    out: list[str] = []
    if _under(pkg_file, REPO):
        out.append(f"package file under {REPO}: {pkg_file}")
    if environ.get("PYTHONPATH"):
        out.append(f"PYTHONPATH set: {environ['PYTHONPATH']}")
    repo_entries = [p for p in sys_path if p and _under(p, REPO)]
    if repo_entries:
        out.append(f"sys.path entries under {REPO}: {repo_entries}")
    if _under(cwd, REPO):
        out.append(f"cwd under {REPO}: {cwd}")
    return out


def test_control_repo_location_is_detected():
    """(c1) A planted repo location returns all four violations, each named."""
    v = _location_violations(
        "/workspace/src/fem_em_solver/__init__.py",
        ["/workspace/src", "/usr/lib/python3.12"],
        {"PYTHONPATH": "/workspace/src"},
        "/workspace",
    )
    assert len(v) == 4, v
    assert v[0].startswith("package file under")
    assert v[1].startswith("PYTHONPATH set")
    assert v[2].startswith("sys.path entries under")
    assert v[3].startswith("cwd under")
    # And a clean location returns none (the checker is not vacuously red).
    assert (
        _location_violations(
            "/tmp/wf13/venv/lib/python3.12/site-packages/fem_em_solver/__init__.py",
            ["/tmp/wf13/venv/lib/python3.12/site-packages", "/usr/lib/python3.12"],
            {},
            "/tmp/wf13/run",
        )
        == []
    )


@pytest.mark.skipif(
    os.environ.get("WF13_INSTALLED_CHECK") != "1",
    reason="WF-13 installed-location check runs only from the /tmp/wf13 venv",
)
def test_installed_package_is_outside_the_repo():
    import importlib.metadata

    import fem_em_solver
    import fem_em_solver.core
    import fem_em_solver.post.evaluation
    import dolfinx
    from mpi4py import MPI
    from petsc4py import PETSc

    import numpy as np

    v = _location_violations(fem_em_solver.__file__, sys.path, os.environ, os.getcwd())
    assert v == [], v

    for mod in (fem_em_solver, fem_em_solver.core, fem_em_solver.post.evaluation):
        assert _under(mod.__file__, sys.prefix), (mod.__name__, mod.__file__, sys.prefix)

    assert importlib.metadata.version("fem-em-solver") == fem_em_solver.__version__

    expect = os.environ["WF13_EXPECT_SCALAR"]
    assert np.dtype(PETSc.ScalarType) == np.dtype(expect), (PETSc.ScalarType, expect)

    if MPI.COMM_WORLD.rank == 0:
        print(
            f"\n[wf13] sys.prefix={sys.prefix} base={sys.base_prefix}"
            f"\n[wf13] fem_em_solver.__file__={fem_em_solver.__file__}"
            f"\n[wf13] fem-em-solver dist version={importlib.metadata.version('fem-em-solver')}"
            f"\n[wf13] dolfinx {dolfinx.__version__} at {dolfinx.__file__}"
            f"\n[wf13] PETSc.ScalarType={np.dtype(PETSc.ScalarType)} (expect {expect})"
            f"\n[wf13] cwd={os.getcwd()} PYTHONPATH={os.environ.get('PYTHONPATH')!r}"
        )
