"""`WF-13` step 5 (§9 item 432, §10 U47): the exit-gate rehearsal, our half --
the installed `WF-1` ``femem run`` from outside the repo reproduces the `WF-15`
sphere record, with the install's version and public surface asserted.

Brief: docs/planning/chunks/WF-13.md, "§9 item 432".  The window installs the
package into the `WF-13` venv (``scripts/package/wf13_setup.sh``), copies
``sphere_mat4.yaml``, ``public_surface.txt``, ``pyproject.toml`` and
``test_wf13_installed_surface.py`` into ``/tmp/wf13/run/u47`` and, there, with
no ``PYTHONPATH``, runs ``mpiexec -n 2 /tmp/wf13/venv/bin/femem run`` and the
venv-side module.  This module (run afterwards, single process, from the
repository) only reads what those wrote, located by ``WF13_REHEARSAL_OUT``.

(V)/(S) ``installed.json``: version equals the repo's ``__version__``; the
        snapshot hash equals the committed ``public_surface.txt``'s.
(P)     summary ``source`` is the rehearsal copy, ``sha256`` the repo file's,
        ``ranks == 2``, ``degree == 1``, ``family == "sphere_in_box"``; the
        console script's shebang is the venv interpreter.
(R)     SAR within ``L0_BAND`` of ``MAT4_DEG1_RECORD``, ``.6e`` string equal to
        ``MAT4_DEG1_PRINT``, within ``MAT4_BAND`` of the closed form.
(c2)    degree confusion: ``|SAR / WF15_DEG2_RECORD - 1| > L0_BAND`` (asserted).

Every record and band is imported by name.  Missing inputs FAIL, never skip.
"""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

import pytest

import fem_em_solver

from tests.materials.test_phantom_material_model import MAT4_BAND
from tests.package.test_wf15_study_lossy_sphere import MAT4_DEG1_PRINT, MAT4_DEG1_RECORD, _closed
from tests.package.test_wf1_study_yaml import SPHERE_FILE, WF15_DEG2_RECORD
from tests.validation.test_lossy_sphere_sar import EPSILON_R_SPHERE, SIGMA_HIGH
from tests.validation.test_mat5_em_coupled_sphere import L0_BAND

HERE = Path(__file__).resolve().parent
SNAPSHOT_PATH = HERE / "public_surface.txt"
VENV = Path("/tmp/wf13/venv")
VENV_PY = VENV / "bin" / "python"
FEMEM = VENV / "bin" / "femem"
RUN_DIR = Path("/tmp/wf13/run/u47")
OUT_ENV = "WF13_REHEARSAL_OUT"
RECORD_RANKS = 2


def _root() -> Path:
    raw = os.environ.get(OUT_ENV)
    if not raw:
        pytest.fail(f"{OUT_ENV} unset: the rehearsal directory is required")
    root = Path(raw)
    assert root.is_dir(), f"{OUT_ENV}={raw} is not a directory"
    return root


def _json(path: Path) -> dict:
    assert path.is_file(), f"missing {path}: the rehearsal wrote nothing there"
    with open(path) as fh:
        return json.load(fh)


def test_wf13_s5_control_dropped_name_changes_snapshot_hash():
    """(c4, §9 item 784) The (S) hash equality can fail: the committed
    snapshot rehashed with one name dropped differs from its own hash."""
    text = SNAPSHOT_PATH.read_text()
    names = [ln for ln in text.splitlines() if ln.strip()]
    assert len(names) >= 2, names
    full = hashlib.sha256(text.encode()).hexdigest()
    dropped = "".join(f"{n}\n" for n in names[1:])
    sha_dropped = hashlib.sha256(dropped.encode()).hexdigest()
    print(
        f"\n[WF-13 s5] (c4) snapshot {len(names)} names sha256 {full}; "
        f"drop {names[0]!r} -> {sha_dropped}",
        flush=True,
    )
    assert full == hashlib.sha256(SNAPSHOT_PATH.read_bytes()).hexdigest()
    assert sha_dropped != full


def test_wf13_exit_rehearsal_reproduces_mat4():
    """(V)/(S) from installed.json, (P), (R), (c2) on the rehearsal outputs."""
    root = _root()
    inst = _json(root / "out" / "installed.json")
    summary = _json(root / "out" / "sphere" / "summary.json")

    # (V) / (S): the install checked this repo's version and committed snapshot.
    snap_sha = hashlib.sha256(SNAPSHOT_PATH.read_bytes()).hexdigest()
    n_snapshot = len([ln for ln in SNAPSHOT_PATH.read_text().splitlines() if ln.strip()])
    print(
        f"\n[WF-13 s5] installed.json {inst}\n"
        f"  repo __version__ {fem_em_solver.__version__!r}; snapshot sha256 {snap_sha} "
        f"({n_snapshot} names)",
        flush=True,
    )
    assert inst["version"] == fem_em_solver.__version__
    assert inst["snapshot_sha256"] == snap_sha
    assert inst["n_names"] == n_snapshot
    assert inst["file"].startswith(str(VENV) + "/"), inst["file"]
    assert not inst["file"].startswith("/workspace"), inst["file"]

    # (P)
    shebang = FEMEM.read_text().splitlines()[0] if FEMEM.is_file() else ""
    repo_sha = hashlib.sha256(SPHERE_FILE.read_bytes()).hexdigest()
    print(
        f"  (P) femem shebang {shebang!r}\n"
        f"  summary source {summary['source']!r}; sha256 {summary['sha256']} vs repo {repo_sha}\n"
        f"  ranks {summary['ranks']}; degree {summary['degree']}; family {summary['family']}; "
        f"cells {summary['cells']}; elapsed_s {summary['elapsed_s']:.1f}",
        flush=True,
    )
    assert shebang == f"#!{VENV_PY}", shebang
    assert summary["source"] == str(RUN_DIR / "sphere_mat4.yaml"), summary["source"]
    assert summary["sha256"] == repo_sha
    assert summary["ranks"] == RECORD_RANKS
    assert summary["degree"] == 1
    assert summary["family"] == "sphere_in_box"

    # (R) and (c2)
    sar = summary["mean_sar"]["mean_sar_w_per_kg"]
    r1 = abs(sar - MAT4_DEG1_RECORD) / MAT4_DEG1_RECORD
    r2 = f"{sar:.6e}"
    c = _closed(SIGMA_HIGH, EPSILON_R_SPHERE)
    g_err = abs(sar - c["sar_w_per_kg"]) / c["sar_w_per_kg"]
    c2 = abs(sar / WF15_DEG2_RECORD - 1.0)
    print(
        f"  (R) SAR {sar!r} rel to MAT4_DEG1_RECORD {r1:.3e} (L0_BAND {L0_BAND:.0e}); "
        f"'{r2}' vs '{MAT4_DEG1_PRINT}'\n"
        f"  closed {c['sar_w_per_kg']:.6e}: err {g_err:.3%} (MAT4_BAND {MAT4_BAND:.0%})\n"
        f"  (c2) |SAR / WF15_DEG2_RECORD - 1| = {c2:.4%} ({c2 / L0_BAND:.2e} x L0_BAND)",
        flush=True,
    )
    assert r1 <= L0_BAND
    assert r2 == MAT4_DEG1_PRINT
    assert g_err <= MAT4_BAND
    assert c2 > L0_BAND
