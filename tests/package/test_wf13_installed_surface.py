"""`WF-13` step 5 (§9 item 432, §10 U47), venv side: the installed package's
location, version and public surface, read from the install itself.

Brief: docs/planning/chunks/WF-13.md, "§9 item 432".  The window copies this
module, ``public_surface.txt`` and ``pyproject.toml`` (and the sphere YAML) into
the rehearsal directory ``WF13_REHEARSAL_DIR`` (``/tmp/wf13/run/u47``) and runs
it there with ``/tmp/wf13/venv/bin/python`` and no ``PYTHONPATH``.  It imports
only the standard library, ``pytest`` and -- inside the node -- the installed
``fem_em_solver`` / ``petsc4py``; never a ``tests/`` name, because the
rehearsal directory is the outside user's stand-in.

(L) ``fem_em_solver.__file__`` under ``sys.prefix``, not under ``/workspace``;
    ``PYTHONPATH`` unset; cwd outside ``/workspace``; ``PETSc.ScalarType`` is
    ``complex128`` (the ``sitecustomize`` guard, item 220).
(V) ``importlib.metadata.version("fem-em-solver") == fem_em_solver.__version__
    == pyproject["project"]["version"]``.
(S) ``sorted(set(fem_em_solver.__all__))`` equals the snapshot's non-blank
    lines and every name resolves by ``getattr``.

The install node is skipped unless ``WF13_INSTALLED_CHECK=1``; on green it
writes ``out/installed.json`` for the repo-side reader
(``test_wf13_exit_rehearsal.py``).  The two control nodes (c1) planted surface
name and (c3) planted version always run.
"""

from __future__ import annotations

import hashlib
import importlib.metadata
import json
import os
import sys
import tomllib
from pathlib import Path

import pytest

REPO = "/workspace"
DIST = "fem-em-solver"
REHEARSAL_ENV = "WF13_REHEARSAL_DIR"


def _under(path: str, root: str) -> bool:
    path = os.path.abspath(path)
    root = os.path.abspath(root)
    return path == root or path.startswith(root.rstrip("/") + "/")


# Restated in stdlib from `WF-14`'s `_surface_violations`
# (tests/package/test_wf14_validation_envelope.py:86-100), surface/snapshot
# half only: the venv cannot import tests/, and the registry half is the
# repo-side envelope's job.
def _surface_violations(surface, snapshot) -> list[str]:
    surface_set, snapshot_set = set(surface), set(snapshot)
    out: list[str] = []
    for n in sorted(surface_set - snapshot_set):
        out.append(f"not in snapshot: {n}")
    for n in sorted(snapshot_set - surface_set):
        out.append(f"in snapshot but gone: {n}")
    return out


def _version_violations(meta, attr, toml) -> list[str]:
    out: list[str] = []
    if meta != toml:
        out.append(f"metadata {meta!r} != pyproject {toml!r}")
    if attr != toml:
        out.append(f"__version__ {attr!r} != pyproject {toml!r}")
    return out


def _snapshot_lines(text: str) -> list[str]:
    return [ln.strip() for ln in text.splitlines() if ln.strip()]


def test_wf13_s5_c1_planted_surface_name_is_reported():
    """(c1) A planted public name is reported exactly; a clean surface is not."""
    snapshot = ["Alpha", "Beta", "Gamma"]
    assert _surface_violations(snapshot + ["_wf13_planted"], snapshot) == [
        "not in snapshot: _wf13_planted"
    ]
    assert _surface_violations(["Alpha", "Beta"], snapshot) == ["in snapshot but gone: Gamma"]
    assert _surface_violations(list(snapshot), snapshot) == []


def test_wf13_s5_c3_planted_version_is_reported():
    """(c3) A planted metadata version is reported exactly; agreement is clean."""
    assert _version_violations("0.0.0", "0.1.0", "0.1.0") == [
        "metadata '0.0.0' != pyproject '0.1.0'"
    ]
    assert _version_violations("0.1.0", "0.1.0", "0.1.0") == []


@pytest.mark.skipif(
    os.environ.get("WF13_INSTALLED_CHECK", "") != "1", reason="WF13_INSTALLED_CHECK unset"
)
def test_wf13_s5_installed_location_version_surface():
    """(L), (V), (S) from the install; writes out/installed.json on green."""
    raw = os.environ.get(REHEARSAL_ENV)
    assert raw, f"{REHEARSAL_ENV} unset"
    rdir = Path(raw)
    surface_file = rdir / "public_surface.txt"
    toml_file = rdir / "pyproject.toml"
    assert surface_file.is_file() and toml_file.is_file(), (surface_file, toml_file)

    import fem_em_solver
    import numpy as np
    from petsc4py import PETSc

    # (L)
    pkg_file = fem_em_solver.__file__
    scalar = np.dtype(PETSc.ScalarType)
    print(
        f"\n[WF-13 s5] (L) sys.prefix {sys.prefix}; fem_em_solver.__file__ {pkg_file}\n"
        f"  cwd {os.getcwd()}; PYTHONPATH {os.environ.get('PYTHONPATH')!r}; "
        f"PETSc.ScalarType {scalar}",
        flush=True,
    )
    assert _under(pkg_file, sys.prefix), (pkg_file, sys.prefix)
    assert not _under(pkg_file, REPO), pkg_file
    assert not _under(sys.prefix, REPO), sys.prefix
    assert "PYTHONPATH" not in os.environ, os.environ.get("PYTHONPATH")
    assert not _under(os.getcwd(), REPO), os.getcwd()
    assert scalar == np.dtype(np.complex128), scalar

    # (V)
    meta = importlib.metadata.version(DIST)
    attr = fem_em_solver.__version__
    with toml_file.open("rb") as fh:
        toml = tomllib.load(fh)["project"]["version"]
    v = _version_violations(meta, attr, toml)
    print(f"  (V) metadata {meta!r}; __version__ {attr!r}; pyproject {toml!r}; violations {v}")
    assert v == [], v

    # (S)
    snap_bytes = surface_file.read_bytes()
    snapshot = _snapshot_lines(snap_bytes.decode())
    surface = sorted(set(fem_em_solver.__all__))
    for n in surface:
        getattr(fem_em_solver, n)  # a dangling __all__ entry raises here
    s = _surface_violations(surface, snapshot)
    sha = hashlib.sha256(snap_bytes).hexdigest()
    print(
        f"  (S) installed __all__ {len(surface)} names; snapshot {len(snapshot)} lines "
        f"(sha256 {sha}); violations {s}",
        flush=True,
    )
    assert s == [], s
    assert surface == sorted(snapshot)

    out = rdir / "out"
    out.mkdir(parents=True, exist_ok=True)
    record = {
        "version": meta,
        "file": pkg_file,
        "n_names": len(surface),
        "snapshot_sha256": sha,
        "sys_prefix": sys.prefix,
    }
    (out / "installed.json").write_text(json.dumps(record, indent=2) + "\n")
    print(f"  wrote {out / 'installed.json'}: {record}", flush=True)
