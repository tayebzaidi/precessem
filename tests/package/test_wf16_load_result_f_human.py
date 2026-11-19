"""`WF-16` step 15 (§9 item 680, §10 U100): ``Study.load_result`` on the 32-port
copper F-human, so a saved human-scale matrix tunes with no solve.

``load_result`` accepts an ``FHumanBirdcageGeometry`` Study (port count
``2 * leg_count`` = 32) and the ring tune reads its port names from
``LoadedPortResult.port_ids``.  Fixture: step 12's copper F-human hole Study
(``_study`` of ``test_wf16_study_tune_f_human_copper.py``), degree 1, 64 MHz,
κ off, ``-n 8``.

Gates (asserted; every band imported, none set here):
  (CM)  the loaded ``C_mode`` equals ``WF12_STEP3_C_MODE_CU_F`` within
        ``REDUCTION_BAND``.
  (T=)  the loaded ``C_mode`` equals the in-memory tune's within
        ``REDUCTION_BAND`` (Touchstone text round-trips at ≈ 6e-14, so no
        ``float.hex`` equality is claimed).
  (S0)  class-level spies on ``Study.solve`` / ``Study._solve_ports`` record 0
        calls on the fresh Study.
  (Ref) ``sar`` / ``b1_plus`` / ``temperature`` on the loaded F-human raise
        ``TypeError``.
Negative control (asserted): the saved file with ports 1 and 2 swapped (U94's
``_swap12`` generalised to 32 ports: data permuted, names kept) misses the
record by more than ``LEG_D0_REPRODUCTION_BAND``; the miss is printed.

Opt-in: ``FEM_EM_WF16_S15=1`` (unset => the solve node skips; any other value
raises); width 8 asserted when set.  The no-mesh node (W0) always runs.
Scope: a saved 32-port matrix tunes with no solve.  Not claimed: a 64-port
(F-clinical) load, a saved-result field output, Route A.
"""

from __future__ import annotations

import os
import time
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from mpi4py import MPI

from fem_em_solver import study as study_mod
from fem_em_solver.ports.definitions import PortDefinition
from fem_em_solver.ports.touchstone import export_touchstone
from fem_em_solver.study import (
    FHumanBirdcageGeometry,
    FourPortBirdcageGeometry,
    LumpedPortDrive,
    SphereInBoxGeometry,
    Study,
)
from tests.complex_mode import complex_only
from tests.package.test_wf16_study_load_result import _synthetic_s, _write
from tests.package.test_wf16_study_tune_capacitors import _study as _four_port_study
from tests.package.test_wf16_study_tune_f_human_copper import (
    OMEGA,
    RECORD_RANK_WIDTH,
    WF12_STEP3_C_MODE_CU_F,
    _study,
    _synthetic,
)
from tests.validation.test_birdcage_sar_mass_averaged import EXACT_IDENTITY_RTOL
from tests.validation.test_port_birdcage_four_port import LEG_D0_REPRODUCTION_BAND
from tests.validation.test_port_lumped_rlc_termination import REDUCTION_BAND

S15_ENV = "FEM_EM_WF16_S15"
OUT_DIR = Path("/workspace/logs/wf16")
N_SLOTS = 16


def _p0(*args):
    if MPI.COMM_WORLD.rank == 0:
        print(*args, flush=True)


def _knob() -> bool:
    raw = os.environ.get(S15_ENV, "").strip()
    if raw == "":
        return False
    if raw == "1":
        return True
    raise ValueError(f"{S15_ENV}={raw!r}: only unset/empty or '1' are accepted")


def _write_named(path, fhz, s, names, comm=MPI.COMM_WORLD):
    """Rank 0 writes one frequency of ``s`` under port ``names``; Barrier."""
    path = Path(path)
    if comm.rank == 0:
        path.parent.mkdir(parents=True, exist_ok=True)
        defs = [
            PortDefinition(port_id=n, positive_tag=100 + i, negative_tag=1, orientation="ring_gap")
            for i, n in enumerate(names)
        ]
        sweep = SimpleNamespace(
            frequency_hz=float(fhz), s_matrix=np.asarray(s, np.complex128), port_ids=tuple(names)
        )
        export_touchstone([sweep], defs, output_path=path, write_csv_companion=False)
    comm.Barrier()
    return path


def _swap12_n(s):
    """U94's ``_swap12`` (``test_wf16_study_load_result.py:66``) on an n x n."""
    n = np.asarray(s).shape[0]
    p = [1, 0] + list(range(2, n))
    return np.asarray(s)[np.ix_(p, p)]


def _expect_ref(st, tag):
    for name, call in (
        ("sar", lambda: st.sar(centres=[(0.0, 0.0, 0.0)])),
        ("b1_plus", lambda: st.b1_plus()),
        ("temperature", lambda: st.temperature(None)),
    ):
        with pytest.raises(TypeError) as ei:
            call()
        _p0(f"[WF-16 s15 {tag} Ref] {name} on the loaded F-human: TypeError: {ei.value}")


@complex_only
def test_wf16_s15_w0_knob_and_synthetic(monkeypatch):
    """W0, no mesh: knob parse; a synthetic 32-port loads (32 ports) and tunes to
    its closed form; a 31-port file is refused; the four-port file still loads."""
    comm = MPI.COMM_WORLD
    for bad in ("2", "yes"):
        monkeypatch.setenv(S15_ENV, bad)
        with pytest.raises(ValueError, match=S15_ENV):
            _knob()
    monkeypatch.delenv(S15_ENV, raising=False)
    assert _knob() is False

    calls = {"n": 0}

    def stub(self, comm=MPI.COMM_WORLD):
        calls["n"] += 1
        raise AssertionError("W0 must not build a mesh")

    for cls in (SphereInBoxGeometry, FHumanBirdcageGeometry, FourPortBirdcageGeometry):
        monkeypatch.setattr(cls, "build", stub)

    s, names, d = _synthetic()
    st = _study()
    fhz = st.frequency_hz
    tag = f"r{comm.size}"
    ok = _write_named(OUT_DIR / f"s15_synthetic_{tag}.s32p", fhz, s, names)
    res = st.load_result(ok)
    assert isinstance(res, study_mod.LoadedPortResult) and st.result is res
    assert len(res.port_ids) == 32 and res.s_matrix.shape == (32, 32)
    assert list(res.port_ids) == names
    rt = float(np.max(np.abs(res.s_matrix - s)))
    t1 = st.tune(1, n_slots=N_SLOTS, comm=comm)
    ref = 1.0 / (OMEGA * d[1].imag)
    e1 = abs(t1.c_mode_farad / ref - 1.0)
    _p0(
        f"\n[WF-16 s15 W0] synthetic 32-port: ports {len(res.port_ids)}, max|dS| {rt:.3e}; "
        f"C_1 {t1.c_mode_farad:.12e} F vs closed {ref:.12e} F, rel {e1:.3e} "
        f"(EXACT_IDENTITY_RTOL {EXACT_IDENTITY_RTOL:.1e})"
    )
    assert e1 <= EXACT_IDENTITY_RTOL
    _expect_ref(st, "W0")

    st31 = _study()
    p31 = _write_named(OUT_DIR / f"s15_synthetic31_{tag}.s31p", fhz, s[:31, :31], names[:31])
    with pytest.raises(ValueError, match="ports") as ei:
        st31.load_result(p31)
    assert st31.result is None
    _p0(f"[WF-16 s15 W0] 31-port file refused: ValueError: {ei.value}")

    fp = _four_port_study(LumpedPortDrive())
    s4 = _synthetic_s()
    p4 = _write(OUT_DIR / f"s15_four_port_{tag}.s4p", [(fp.frequency_hz, s4)])
    r4 = fp.load_result(p4)
    assert len(r4.port_ids) == 4 and float(np.max(np.abs(r4.s_matrix - s4))) <= EXACT_IDENTITY_RTOL
    _p0(f"[WF-16 s15 W0] four-port synthetic still loads: ports {r4.port_ids}")
    assert calls["n"] == 0


@pytest.fixture(scope="module")
def s15():
    if not _knob():
        pytest.skip(f"{S15_ENV} unset: the copper F-human solve is a heavy -n 8 window")
    comm = MPI.COMM_WORLD
    assert comm.size == RECORD_RANK_WIDTH, (comm.size, RECORD_RANK_WIDTH)
    st = _study()
    t0 = time.perf_counter()
    res = st.solve(comm)
    t_solve = time.perf_counter() - t0
    tm = st.tune(1, n_slots=N_SLOTS, comm=comm)
    path = res.write_touchstone(OUT_DIR / f"s15_f_human_cu_d1_n{comm.size}.s32p")
    _p0(f"\n[WF-16 s15] cells {res.cells}, solve {t_solve:.1f} s, -n {comm.size}; wrote {path}")
    return {"path": path, "tune_mem": tm, "fhz": st.frequency_hz}


@complex_only
def test_wf16_s15_tune_on_loaded_f_human(s15):
    comm = MPI.COMM_WORLD
    fresh = _study()
    calls = {"solve": 0, "_solve_ports": 0}

    def _spy(name):
        def _boom(self, *a, **k):
            calls[name] += 1
            raise AssertionError(f"(S0) the loaded Study called {name}")

        return _boom

    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(Study, "solve", _spy("solve"))
        mp.setattr(Study, "_solve_ports", _spy("_solve_ports"))
        loaded = fresh.load_result(s15["path"])
        assert len(loaded.port_ids) == 32
        tl = fresh.tune(1, n_slots=N_SLOTS, comm=comm)
        _expect_ref(fresh, "W1")
    _p0(
        f"[WF-16 s15 S0] class-level spies: solve {calls['solve']}, "
        f"_solve_ports {calls['_solve_ports']}"
    )
    assert calls == {"solve": 0, "_solve_ports": 0}

    cl = tl.c_mode_farad
    cm = s15["tune_mem"].c_mode_farad
    g_cm = abs(cl / WF12_STEP3_C_MODE_CU_F - 1.0)
    g_t = abs(cl / cm - 1.0)
    g_cm_g = comm.allreduce(g_cm, op=MPI.MAX)
    g_t_g = comm.allreduce(g_t, op=MPI.MAX)
    _p0(
        f"[WF-16 s15 CM] loaded C_mode {cl:.12e} F vs record {WF12_STEP3_C_MODE_CU_F:.12e} F: "
        f"rel {g_cm_g:.3e} (REDUCTION_BAND {REDUCTION_BAND:.1e})"
    )
    _p0(
        f"[WF-16 s15 T=] loaded {cl:.12e} F vs in-memory {cm:.12e} F: rel {g_t_g:.3e} "
        f"(REDUCTION_BAND {REDUCTION_BAND:.1e}); off-block loaded "
        f"{tl.off_diagonal_residual:.6e} vs in-memory {s15['tune_mem'].off_diagonal_residual:.6e}"
    )
    _p0(f"[WF-16 s15] port order (file): {list(loaded.port_ids)}")
    assert g_cm_g <= REDUCTION_BAND
    assert g_t_g <= REDUCTION_BAND

    # Negative control: ports 1 <-> 2 swapped in the data, names kept.
    sw = _write_named(
        OUT_DIR / f"s15_swap12_n{comm.size}.s32p",
        s15["fhz"],
        _swap12_n(loaded.s_matrix),
        list(loaded.port_ids),
    )
    st_sw = _study()
    st_sw.load_result(sw)
    try:
        cw = st_sw.tune(1, n_slots=N_SLOTS, comm=comm).c_mode_farad
        miss = abs(cw / WF12_STEP3_C_MODE_CU_F - 1.0)
        miss_g = comm.allreduce(miss, op=MPI.MAX)
        _p0(
            f"[WF-16 s15 Sw] P1<->P2 swapped: C {cw:.12e} F, miss {miss_g:.6e} vs "
            f"LEG_D0_REPRODUCTION_BAND {LEG_D0_REPRODUCTION_BAND:.1e}"
        )
    except ValueError as exc:  # a non-capacitive read is a miss too; printed
        miss_g = float("inf")
        _p0(f"[WF-16 s15 Sw] P1<->P2 swapped: tune raised {exc} (miss = inf)")
    assert miss_g > LEG_D0_REPRODUCTION_BAND
