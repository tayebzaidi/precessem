"""(N5) controls for `tests/validation/gmsh_stop_rule.py` (§9 item 522).

Pure; run at ``-n 1`` only (W0). Each child is a bare ``python3 -c`` that
imports only the stdlib helper — never ``mpi4py`` — so it is not an MPI rank.
"""

from __future__ import annotations

import os
import subprocess
import sys
import time

import pytest

from tests.validation.gmsh_stop_rule import (
    GMSH_STOP_ENV,
    GMSH_STOP_MARKER,
    GMSH_STOP_RC,
    gmsh_stop_armed,
)

_REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def _child(code):
    env = dict(os.environ)
    env["PYTHONPATH"] = _REPO
    t0 = time.perf_counter()
    p = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True,
                       env=env, cwd=_REPO, timeout=30)
    return p, time.perf_counter() - t0


def test_n5a_fires_across_a_c_call():
    code = (
        "import ctypes\n"
        "from tests.validation.gmsh_stop_rule import _StopWatch\n"
        "w = _StopWatch(1.0, 'N5a'); w.arm()\n"
        "ctypes.CDLL(None).sleep(3)\n"
        "print('fell through')\n"
    )
    p, dt = _child(code)
    print(f"\n[N5a] rc {p.returncode}, {dt:.3f} s, stdout {p.stdout!r}", flush=True)
    assert p.returncode == GMSH_STOP_RC
    assert GMSH_STOP_MARKER in p.stdout
    assert "fell through" not in p.stdout
    assert dt < 2.9


def test_n5b_cancel_is_silent():
    code = (
        "import ctypes, time\n"
        "from tests.validation.gmsh_stop_rule import _StopWatch\n"
        "w = _StopWatch(2.0, 'N5b'); w.arm()\n"
        "ctypes.CDLL(None).usleep(100000)\n"
        "w.cancel()\n"
        "time.sleep(2.5)\n"
        "print('done')\n"
    )
    p, dt = _child(code)
    print(f"\n[N5b] rc {p.returncode}, {dt:.3f} s, stdout {p.stdout!r}", flush=True)
    assert p.returncode == 0
    assert GMSH_STOP_MARKER not in p.stdout
    assert "done" in p.stdout


def test_n5c_gmsh_lib_releases_the_gil():
    import ctypes

    import gmsh

    print(f"\n[N5c] gmsh.lib type {type(gmsh.lib).__name__}", flush=True)
    assert isinstance(gmsh.lib, ctypes.CDLL)
    assert not isinstance(gmsh.lib, ctypes.PyDLL)


def test_armed_parse(monkeypatch):
    monkeypatch.delenv(GMSH_STOP_ENV, raising=False)
    assert gmsh_stop_armed() is False
    monkeypatch.setenv(GMSH_STOP_ENV, "")
    assert gmsh_stop_armed() is False
    monkeypatch.setenv(GMSH_STOP_ENV, "1")
    assert gmsh_stop_armed() is True
    monkeypatch.setenv(GMSH_STOP_ENV, "450")
    with pytest.raises(ValueError):
        gmsh_stop_armed()
