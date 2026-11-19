"""`WF-13` step 3 (§9 item 294, chain U19): ``femem run --ranks N``, the MPI
launch helper, run from the *installed* package.

The window (brief: docs/planning/chunks/WF-13.md, "§9 item 294") installs the
package into `WF-13`'s venv (``scripts/package/wf13_setup.sh``), copies `WF-1`'s
sphere file to ``/tmp/wf13/run/u19/`` and runs, with no ``PYTHONPATH``,
``/tmp/wf13/venv/bin/femem run --ranks N sphere_mat4.yaml`` for N = 1 and 2.
The helper re-executes ``mpiexec -n N <venv python> -m fem_em_solver run ...``
(the ``__main__`` route).  This module (single process, from the repository)
reads what those runs wrote, located by ``WF13_LAUNCH_OUT``.  Asserted:

(L) the launch happened at the asked width: ``summary["ranks"] == N``; the
    printed child command carries ``-n N`` and ``-m fem_em_solver`` and the
    launcher's exit code is passed through (stub launcher, no solve);
(P) provenance: the venv imports ``fem_em_solver`` from the venv, not
    ``/workspace/src``; each summary's ``source`` lies under
    ``/tmp/wf13/run/u19/``;
(R) at each width SAR within ``L0_BAND`` of ``MAT4_DEG1_RECORD``, ``.6e``
    string ``== MAT4_DEG1_PRINT``, within ``MAT4_BAND`` of the closed form;
    the cross-width ratio is printed, not asserted;
(X) refusals, no solve: ``--ranks 0`` returns 2; ``--ranks`` inside an
    ``mpiexec -n 2`` launch returns 2 on both ranks and names the variable;
    neither writes an output directory.

Negative controls (asserted): (c1) launcher blindness — a helper that dropped
``-n`` would read ``ranks == 1`` at N = 2 (structural: (L) at N = 2);
(c2) degree confusion — ``|SAR / WF15_DEG2_RECORD - 1| > L0_BAND`` at both
widths.

Every record and band is imported by name.  Missing inputs FAIL, never skip.
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
from pathlib import Path

import pytest

from tests.materials.test_phantom_material_model import MAT4_BAND
from tests.package.test_wf1_study_yaml import SPHERE_FILE, WF15_DEG2_RECORD
from tests.package.test_wf15_study_lossy_sphere import MAT4_DEG1_PRINT, MAT4_DEG1_RECORD, _closed
from tests.validation.test_lossy_sphere_sar import EPSILON_R_SPHERE, SIGMA_HIGH
from tests.validation.test_mat5_em_coupled_sphere import L0_BAND

VENV = Path("/tmp/wf13/venv")
VENV_PY = VENV / "bin" / "python"
FEMEM = VENV / "bin" / "femem"
RUN_DIR = Path("/tmp/wf13/run/u19")
OUT_ENV = "WF13_LAUNCH_OUT"
WIDTHS = (1, 2)
REFUSAL = "femem: error: --ranks inside an MPI launch"


def _out_root() -> Path:
    raw = os.environ.get(OUT_ENV)
    if not raw:
        pytest.fail(f"{OUT_ENV} unset: the `femem run --ranks N` output root is required")
    root = Path(raw)
    assert root.is_dir(), f"{OUT_ENV}={raw} is not a directory"
    return root


def _summary(n: int) -> dict:
    path = _out_root() / f"n{n}" / "summary.json"
    assert path.is_file(), f"missing {path}: `femem run --ranks {n}` wrote no summary"
    with open(path) as fh:
        return json.load(fh)


def _venv_env(**extra) -> dict:
    env = {k: v for k, v in os.environ.items() if k != "PYTHONPATH"}
    env.update(extra)
    return env


def _sar(summary: dict) -> float:
    return float(summary["mean_sar"]["mean_sar_w_per_kg"])


def test_wf13_launch_provenance_and_width():
    """(P), (L), (c1) on both installed launches."""
    proc = subprocess.run(
        [str(VENV_PY), "-c", "import fem_em_solver; print(fem_em_solver.__file__)"],
        cwd="/tmp",
        env=_venv_env(),
        capture_output=True,
        text=True,
        timeout=60,
    )
    where = proc.stdout.strip().splitlines()[-1] if proc.stdout.strip() else ""
    print(
        f"\n[WF-13 s3] (P) venv fem_em_solver.__file__ = {where!r} (rc {proc.returncode})",
        flush=True,
    )
    assert proc.returncode == 0, proc.stderr[-2000:]
    assert where.startswith(str(VENV) + "/"), where
    assert not where.startswith("/workspace"), where

    for n in WIDTHS:
        s = _summary(n)
        print(
            f"  N={n}: ranks {s['ranks']}; source {s['source']!r}; degree {s['degree']}; "
            f"cells {s['cells']}; elapsed_s {s['elapsed_s']:.1f}",
            flush=True,
        )
        assert s["femem"] == 1
        assert s["family"] == "sphere_in_box"
        assert s["degree"] == 1
        assert s["source"].startswith(str(RUN_DIR) + "/"), s["source"]
        assert Path(s["source"]).name == SPHERE_FILE.name
        assert s["ranks"] == n  # (L)
    # (c1) a helper that dropped `-n` launches a singleton: ranks 1 at N = 2.
    assert _summary(2)["ranks"] != 1


@pytest.mark.parametrize("n", WIDTHS)
def test_wf13_launch_reproduces_mat4_record(n):
    """(R), (c2) at width N."""
    sar = _sar(_summary(n))
    r1 = abs(sar - MAT4_DEG1_RECORD) / MAT4_DEG1_RECORD
    r2 = f"{sar:.6e}"
    c = _closed(SIGMA_HIGH, EPSILON_R_SPHERE)
    g_err = abs(sar - c["sar_w_per_kg"]) / c["sar_w_per_kg"]
    c2 = abs(sar / WF15_DEG2_RECORD - 1.0)
    print(
        f"\n[WF-13 s3] (R) N={n}: SAR {sar!r} rel to MAT4_DEG1_RECORD {r1:.3e} "
        f"(L0_BAND {L0_BAND:.0e}); '{r2}' vs '{MAT4_DEG1_PRINT}'\n"
        f"  closed {c['sar_w_per_kg']:.6e}: err {g_err:.3%} (MAT4_BAND {MAT4_BAND:.0%})\n"
        f"  (c2) |SAR / WF15_DEG2_RECORD - 1| = {c2:.4%} ({c2 / L0_BAND:.2e} x L0_BAND)",
        flush=True,
    )
    if n == WIDTHS[-1]:
        x = abs(_sar(_summary(1)) / _sar(_summary(2)) - 1.0)
        print(f"  cross-width |SAR_1 / SAR_2 - 1| = {x:.3e} (printed, not asserted)", flush=True)
    assert r1 <= L0_BAND
    assert r2 == MAT4_DEG1_PRINT
    assert g_err <= MAT4_BAND
    assert c2 > L0_BAND


def test_wf13_launch_command_line_and_exit_code(tmp_path):
    """(L) the printed child command, and exit-code pass-through, via a stub
    launcher (no solve, no MPI)."""
    true_bin, false_bin = shutil.which("true"), shutil.which("false")
    assert true_bin and false_bin
    out = tmp_path / "stub"
    for stub, want in ((true_bin, 0), (false_bin, 1)):
        proc = subprocess.run(
            [str(FEMEM), "run", "--ranks", "2", str(SPHERE_FILE), "--output-dir", str(out)],
            cwd="/tmp",
            env=_venv_env(FEMEM_MPIEXEC=stub),
            capture_output=True,
            text=True,
            timeout=60,
        )
        line = [ln for ln in proc.stdout.splitlines() if ln.startswith("femem: launching: ")]
        print(f"\n[WF-13 s3] stub {stub}: rc {proc.returncode}; {line!r}", flush=True)
        assert proc.returncode == want, proc.stderr[-2000:]
        assert len(line) == 1
        argv = line[0].split()
        assert argv[2] == stub
        assert argv[3:5] == ["-n", "2"]
        assert argv[5] == str(VENV_PY)
        assert argv[6:9] == ["-m", "fem_em_solver", "run"]
    assert not out.exists()

    proc = subprocess.run(
        [str(FEMEM), "run", "--ranks", "2", str(SPHERE_FILE)],
        cwd="/tmp",
        env=_venv_env(FEMEM_MPIEXEC="/nonexistent/mpiexec"),
        capture_output=True,
        text=True,
        timeout=60,
    )
    print(f"  absent launcher: rc {proc.returncode}; stderr {proc.stderr.strip()!r}", flush=True)
    assert proc.returncode == 2
    assert "femem: error: MPI launcher not found" in proc.stderr


def test_wf13_launch_refusals(tmp_path):
    """(X): `--ranks 0` (installed), and `--ranks` inside `mpiexec -n 2`
    (repo route); neither solves nor writes an output directory."""
    out0 = tmp_path / "zero"
    proc = subprocess.run(
        [str(FEMEM), "run", "--ranks", "0", str(SPHERE_FILE), "--output-dir", str(out0)],
        cwd="/tmp",
        env=_venv_env(),
        capture_output=True,
        text=True,
        timeout=60,
    )
    print(
        f"\n[WF-13 s3] (X) --ranks 0: rc {proc.returncode}; stderr {proc.stderr.strip()!r}",
        flush=True,
    )
    assert proc.returncode == 2
    assert "--ranks must be >= 1" in proc.stderr
    assert "femem: launching" not in proc.stdout
    assert not out0.exists()

    out2 = tmp_path / "nested"
    proc = subprocess.run(
        [
            "mpiexec",
            "-n",
            "2",
            "python3",
            "-m",
            "fem_em_solver",
            "run",
            "--ranks",
            "2",
            str(SPHERE_FILE),
            "--output-dir",
            str(out2),
        ],
        cwd="/workspace",
        capture_output=True,
        text=True,
        timeout=120,
    )
    refusals = [ln for ln in proc.stderr.splitlines() if ln.startswith(REFUSAL)]
    print(f"  nested: rc {proc.returncode}; refusals {refusals!r}", flush=True)
    assert proc.returncode == 2, proc.stderr[-2000:]
    assert len(refusals) == 2  # one per rank
    for ln in refusals:
        assert any(
            k in ln for k in ("PMI_SIZE", "PMI_RANK", "OMPI_COMM_WORLD_SIZE", "PMIX_RANK")
        ), ln
    assert "femem: launching" not in proc.stdout
    assert not out2.exists()
