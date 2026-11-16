"""`WF-1` step 1 (§9 item 256), W5: the installed ``femem`` console script.

Stdlib only — no import from ``tests/`` and none of ``fem_em_solver`` — so it
runs in `WF-13`'s fresh venv (``scripts/package/wf13_setup.sh real``) against
the wheel built from the copied ``pyproject.toml``.  Skipped unless
``WF13_INSTALLED_CHECK=1``.  (X′): the script exists, its shebang names the
venv's interpreter, ``femem check`` accepts the good files (exit 0) and
refuses the three malformed ones (exit 2) naming the field on stderr.
"""

from __future__ import annotations

import os
import subprocess
import sys
import time
from pathlib import Path

import pytest

WF1 = Path(__file__).resolve().parent / "wf1"
FEMEM = Path(sys.prefix) / "bin" / "femem"
MALFORMED = {
    "malformed_unknown_key.yaml": "geometry.sphere_raduis",
    "malformed_wrong_type.yaml": "study.frequency_hz",
    "malformed_missing_section.yaml": "excitation",
}

pytestmark = pytest.mark.skipif(
    os.environ.get("WF13_INSTALLED_CHECK", "") != "1", reason="WF13_INSTALLED_CHECK unset"
)


def _check(name):
    started = time.perf_counter()
    proc = subprocess.run(
        [str(FEMEM), "check", str(WF1 / name)], capture_output=True, text=True, timeout=60
    )
    seconds = time.perf_counter() - started
    print(f"\n[WF-1] femem check {name}: rc {proc.returncode} ({seconds:.1f} s)")
    print(f"  stdout: {proc.stdout.strip()!r}")
    print(f"  stderr: {proc.stderr.strip()[-400:]!r}")
    return proc


def test_wf1_installed_femem_refuses_the_malformed_files():
    print(f"\n[WF-1] sys.prefix {sys.prefix}; femem {FEMEM}")
    assert FEMEM.is_file(), f"{FEMEM} was not installed"
    shebang = FEMEM.read_text().splitlines()[0]
    print(f"  shebang {shebang!r}")
    assert shebang.startswith("#!")
    assert str(Path(sys.prefix) / "bin" / "python") in shebang

    good = _check("sphere_mat4.yaml")
    assert good.returncode == 0
    assert "sphere_in_box" in good.stdout and "degree: 1" in good.stdout
    assert _check("f_small_four_port.yaml").returncode == 0
    assert _check("control_sigma_ulp.yaml").returncode == 0

    for name, field in MALFORMED.items():
        proc = _check(name)
        assert proc.returncode == 2, name
        assert field in proc.stderr, (name, field)
