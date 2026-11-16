"""`WF-1` step 2 (§9 item 293, chain U18): ``femem run`` with a solve through
the *installed* console script, outside the repository.

The window (brief: docs/planning/chunks/WF-1.md, "§9 item 293") installs the
package into `WF-13`'s fresh venv (``scripts/package/wf13_setup.sh``), copies
one committed `WF-1` YAML file to ``/tmp/wf13/run/u18/`` and runs
``mpiexec -n 2 /tmp/wf13/venv/bin/femem run <file> --output-dir ...`` there,
with no ``PYTHONPATH``.  This module (run afterwards, single process, from the
repository) only reads what that run wrote, located by ``WF1_INSTALLED_OUT``.
Asserted:

(P)   provenance: the venv imports ``fem_em_solver`` from the venv, not from
      ``/workspace/src``; the console script's shebang is the venv's
      interpreter; each summary's ``source`` lies under ``/tmp/wf13/run/u18/``,
      its ``sha256`` equals the repo file's, ``ranks`` equals the record
      width, ``degree == 1``;
(R-S) sphere: SAR within ``L0_BAND`` of ``MAT4_DEG1_RECORD``, ``.6e`` string
      ``== MAT4_DEG1_PRINT``, within ``MAT4_BAND`` of the closed form;
(R-F) F-small: S within ``LEG_D0_REPRODUCTION_BAND`` of ``P1`` at 64 MHz; the
      installed-written ``.s4p`` read back by ``load_touchstone`` and
      scikit-rf equals the summary S within ``ANCHOR_TOL``; `PORT-11` gates.

Negative controls (asserted): (c2) sphere vs the degree-2 record; (c3)
F-small vs ``P2_X1``; (c4) P1<->P2 swap misses ``P1``.

Every record and band is imported by name.  Missing inputs FAIL, never skip:
a skipped node reads green in the harness footer.
"""

from __future__ import annotations

import hashlib
import json
import os
import subprocess
from pathlib import Path

import numpy as np
import pytest

from fem_em_solver.ports.touchstone import load_touchstone

from tests.io.test_touchstone_skrf_crosscheck import ANCHOR_TOL, FREQ_ATOL_HZ
from tests.materials.test_phantom_material_model import MAT4_BAND
from tests.package.test_wf15_study_lossy_sphere import MAT4_DEG1_PRINT, MAT4_DEG1_RECORD, _closed
from tests.package.test_wf16_study_four_port import _record, _rel
from tests.package.test_wf1_study_yaml import BIRDCAGE_FILE, SPHERE_FILE, WF15_DEG2_RECORD
from tests.validation.test_ans4_exact_strip_records import ANS4_EXACT_STRIP_RECORD_RANK_WIDTH
from tests.validation.test_lossy_sphere_fullwave import FREQUENCY_64_HZ
from tests.validation.test_lossy_sphere_sar import EPSILON_R_SPHERE, SIGMA_HIGH
from tests.validation.test_mat5_em_coupled_sphere import L0_BAND
from tests.validation.test_ops61_exact_strip_ports import ANS4_EXACT_STRIP_CELL_RECORD
from tests.validation.test_port_birdcage_four_port import (
    LEG_D0_REPRODUCTION_BAND,
    PASSIVITY_SIGMA_TOLERANCE,
    POOLED_SEPARATION_FLOOR,
    _circulant_classes,
    _class_spread,
)
from tests.validation.test_port_birdcage_lumped_column import ADJACENT_SPREAD_BAND
from tests.validation.test_port_lumped_sheet_sweep import RECIPROCITY_BAND

VENV = Path("/tmp/wf13/venv")
VENV_PY = VENV / "bin" / "python"
FEMEM = VENV / "bin" / "femem"
RUN_DIR = Path("/tmp/wf13/run/u18")
OUT_ENV = "WF1_INSTALLED_OUT"


def _out_root() -> Path:
    raw = os.environ.get(OUT_ENV)
    if not raw:
        pytest.fail(f"{OUT_ENV} unset: the installed `femem run` output root is required")
    root = Path(raw)
    assert root.is_dir(), f"{OUT_ENV}={raw} is not a directory"
    return root


def _summary(sub: str) -> dict:
    path = _out_root() / sub / "summary.json"
    assert path.is_file(), f"missing {path}: the installed `femem run` wrote no summary"
    with open(path) as fh:
        return json.load(fh)


def _provenance(summary: dict, repo_file: Path) -> None:
    """(P), shared by both nodes."""
    env = {k: v for k, v in os.environ.items() if k != "PYTHONPATH"}
    proc = subprocess.run(
        [str(VENV_PY), "-c", "import fem_em_solver; print(fem_em_solver.__file__)"],
        cwd="/tmp",
        env=env,
        capture_output=True,
        text=True,
        timeout=60,
    )
    where = proc.stdout.strip().splitlines()[-1] if proc.stdout.strip() else ""
    shebang = FEMEM.read_text().splitlines()[0] if FEMEM.is_file() else ""
    repo_sha = hashlib.sha256(repo_file.read_bytes()).hexdigest()
    print(
        f"\n[WF-1 s2] (P) venv fem_em_solver.__file__ = {where!r} (rc {proc.returncode})\n"
        f"  femem shebang {shebang!r}\n"
        f"  summary source {summary['source']!r}; sha256 {summary['sha256']} vs repo {repo_sha}\n"
        f"  ranks {summary['ranks']} (record width {ANS4_EXACT_STRIP_RECORD_RANK_WIDTH}); "
        f"degree {summary['degree']}; family {summary['family']}; cells {summary['cells']}; "
        f"elapsed_s {summary['elapsed_s']:.1f}",
        flush=True,
    )
    assert proc.returncode == 0, proc.stderr[-2000:]
    assert where.startswith(str(VENV) + "/"), where
    assert not where.startswith("/workspace"), where
    assert shebang == f"#!{VENV_PY}", shebang
    assert summary["femem"] == 1
    assert summary["source"].startswith(str(RUN_DIR) + "/"), summary["source"]
    assert Path(summary["source"]).name == repo_file.name
    assert summary["sha256"] == repo_sha
    assert summary["ranks"] == ANS4_EXACT_STRIP_RECORD_RANK_WIDTH
    assert summary["degree"] == 1


def test_wf1_installed_sphere_reproduces_mat4():
    """(P), (R-S), (c2) on the installed `femem run sphere_mat4.yaml` output."""
    summary = _summary("sphere")
    _provenance(summary, SPHERE_FILE)
    assert summary["family"] == "sphere_in_box"

    sar = summary["mean_sar"]["mean_sar_w_per_kg"]
    r1 = abs(sar - MAT4_DEG1_RECORD) / MAT4_DEG1_RECORD
    r2 = f"{sar:.6e}"
    c = _closed(SIGMA_HIGH, EPSILON_R_SPHERE)
    g_err = abs(sar - c["sar_w_per_kg"]) / c["sar_w_per_kg"]
    c2 = abs(sar / WF15_DEG2_RECORD - 1.0)
    print(
        f"  (R-S) SAR {sar!r} rel to MAT4_DEG1_RECORD {r1:.3e} (L0_BAND {L0_BAND:.0e}); "
        f"'{r2}' vs '{MAT4_DEG1_PRINT}'\n"
        f"  closed {c['sar_w_per_kg']:.6e}: err {g_err:.3%} (MAT4_BAND {MAT4_BAND:.0%})\n"
        f"  cells {summary['cells']} (WF-1 W1 read 74020; printed, not asserted)\n"
        f"  (c2) |SAR / WF15_DEG2_RECORD - 1| = {c2:.4%} ({c2 / L0_BAND:.2e} x L0_BAND)",
        flush=True,
    )
    assert r1 <= L0_BAND
    assert r2 == MAT4_DEG1_PRINT
    assert g_err <= MAT4_BAND
    assert c2 > L0_BAND


def test_wf1_installed_f_small_reproduces_p1():
    """(P), (R-F), (c3), (c4) on the installed `femem run f_small_four_port.yaml` output."""
    import skrf

    summary = _summary("f_small")
    _provenance(summary, BIRDCAGE_FILE)
    assert summary["family"] == "four_port_birdcage"
    assert summary["touchstone"] == "f_small_four_port.s4p"
    s = np.array([[complex(*p) for p in row] for row in summary["s_matrix"]])
    z = np.array([[complex(*p) for p in row] for row in summary["z_matrix"]])
    assert s.shape == (4, 4) and z.shape == (4, 4)

    p1 = _record(1)
    r1 = _rel(s, p1)
    c3 = _rel(s, _record(2))
    ts = _out_root() / "f_small" / summary["touchstone"]
    assert ts.is_file(), f"missing {ts}"
    freqs, mats, ports = load_touchstone(ts)
    loaded = np.asarray(mats[0], dtype=np.complex128)
    ntw = skrf.Network(str(ts))
    sk = np.asarray(ntw.s[0], dtype=np.complex128)
    d_load = float(np.max(np.abs(loaded - s)))
    d_skrf = float(np.max(np.abs(sk - s)))
    t2 = _rel(loaded, p1)
    perm = [1, 0, 2, 3]
    c4 = _rel(loaded[np.ix_(perm, perm)], p1)
    recip = float(np.linalg.norm(s - s.T) / np.linalg.norm(s))
    sigma = np.linalg.svd(s, compute_uv=False)
    power = np.sum(np.abs(s) ** 2, axis=0)
    classes = _circulant_classes(z)
    spreads = {n: _class_spread(v) for n, v in classes.items()}
    pooled = _class_spread(np.concatenate([classes["adjacent"], classes["opposite"]]))
    worst = max(spreads.values())
    separation = pooled / worst if worst > 0.0 else np.inf
    print(
        f"  cells {summary['cells']} vs ANS4_EXACT_STRIP_CELL_RECORD "
        f"{ANS4_EXACT_STRIP_CELL_RECORD} (printed, not asserted)\n"
        f"  (R-F) ||S - P1||/||P1|| {r1:.3e} (band {LEG_D0_REPRODUCTION_BAND:.0e}); "
        f"(c3) vs P2_X1 {c3:.4f}\n"
        f"  Touchstone max|dS| load_touchstone {d_load:.3e}, skrf {d_skrf:.3e} (ANCHOR_TOL "
        f"{ANCHOR_TOL:.0e}); f {float(freqs[0])!r}; ports {tuple(ports)}; nports {ntw.nports}\n"
        f"  ||S_file - P1||/||P1|| {t2:.3e}; (c4) P1<->P2 swap {c4:.4f}\n"
        f"  PORT-11: reciprocity {recip:.3e}; sigma_max {float(np.max(sigma)):.12f}; max column "
        f"power {float(np.max(power)):.12f}; spreads "
        + ", ".join(f"{k} {v * 100:.4f}%" for k, v in spreads.items())
        + f"; separation {separation:.2f}x",
        flush=True,
    )
    assert r1 <= LEG_D0_REPRODUCTION_BAND
    assert c3 > LEG_D0_REPRODUCTION_BAND
    assert d_load <= ANCHOR_TOL and d_skrf <= ANCHOR_TOL
    assert abs(float(freqs[0]) - FREQUENCY_64_HZ) <= FREQ_ATOL_HZ
    assert abs(float(ntw.f[0]) - FREQUENCY_64_HZ) <= FREQ_ATOL_HZ
    assert ntw.nports == 4 and tuple(ports) == ("P1", "P2", "P3", "P4")
    assert np.all(np.asarray(ntw.z0) == 50.0)
    assert t2 <= LEG_D0_REPRODUCTION_BAND
    assert c4 > LEG_D0_REPRODUCTION_BAND
    assert recip <= RECIPROCITY_BAND
    assert float(np.max(sigma)) <= 1.0 + PASSIVITY_SIGMA_TOLERANCE
    for v in power:
        assert v <= 1.0 + PASSIVITY_SIGMA_TOLERANCE
    for name, v in spreads.items():
        assert v <= ADJACENT_SPREAD_BAND, name
    assert separation >= POOLED_SEPARATION_FLOOR
