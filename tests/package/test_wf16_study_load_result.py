"""`WF-16` step 14 (§9 item 661, §10 U94): tune a saved four-port result, no solve.

``Study.load_result(path)`` reads a single-frequency Touchstone file written by
``PortStudyResult.write_touchstone`` into a ``LoadedPortResult`` (S matrix only);
``Study.tune(k=1)`` then runs on it without solving.  Fixture: step 6(b)'s
degree-1 exact-rectangle PEC hole F-small Study (``_study`` of
``test_wf16_study_tune_capacitors.py``), 64 MHz, ``-n 2``, κ off.

Anchors: (RT) the loaded S equals the solved S to ``ANCHOR_TOL``; (CM) the fresh
Study's ``C_mode`` within ``LEG_D0_REPRODUCTION_BAND`` of the in-memory tune and
``.9e``-equal to item 429's ``C_mode``; (S0) the fresh Study's solve entry is
monkeypatched to raise and the tune still returns; (Ref) refusals before
storing.  Negative controls: (Sw) P1 <-> P2 swapped file (asserted only if the
tune is measured to move; printed either way) and (Sw') S scaled by 0.9
(asserted).  Route A on a loaded matrix is not claimed (returned to the 10-28
weekly).  Every band is imported; none is set here.
"""

from __future__ import annotations

import time
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from mpi4py import MPI

from fem_em_solver import LumpedPortDrive, Study
from fem_em_solver import study as study_mod
from fem_em_solver.ports.definitions import PortDefinition
from fem_em_solver.ports.touchstone import export_touchstone

from tests.complex_mode import complex_only
from tests.io.test_touchstone_skrf_crosscheck import ANCHOR_TOL
from tests.package.test_wf16_study_tune_capacitors import ITEM429_C_MODE_F, _study
from tests.validation.test_port_birdcage_four_port import LEG_D0_REPRODUCTION_BAND

OUT_DIR = Path("/workspace/logs/wf16_s14")


def _p0(*args):
    if MPI.COMM_WORLD.rank == 0:
        print(*args, flush=True)


def _defs(n=4):
    return [PortDefinition(port_id=f"P{i + 1}", positive_tag=211 + i, negative_tag=1,
                           orientation="leg_gap_axial_plus_z") for i in range(n)]


def _write(path, s_by_freq, comm=MPI.COMM_WORLD):
    """Rank 0 writes ``[(f, S), ...]`` through ``export_touchstone``; Barrier."""
    path = Path(path)
    n = np.asarray(s_by_freq[0][1]).shape[0]
    if comm.rank == 0:
        path.parent.mkdir(parents=True, exist_ok=True)
        sweeps = [SimpleNamespace(frequency_hz=float(f), s_matrix=np.asarray(s, np.complex128),
                                  port_ids=tuple(f"P{i + 1}" for i in range(n)))
                  for f, s in s_by_freq]
        export_touchstone(sweeps, _defs(n), output_path=path, write_csv_companion=False)
    comm.Barrier()
    return path


def _swap12(s):
    p = [1, 0, 2, 3]
    return np.asarray(s)[np.ix_(p, p)]


def _synthetic_s():
    # A C4-circulant lossy reciprocal 4x4 with capacitive k = 1 mode (coding only).
    a, b, c = 0.2 + 0.5j, -0.3 + 0.1j, 0.05 - 0.2j
    row = np.array([a, b, c, b])
    return np.array([np.roll(row, i) for i in range(4)])


@complex_only
def test_wf16_s14_load_result_coding():
    """(Ref) + coding: refusals before storing, the field-method guards, the
    swap check on a synthetic circulant.  No mesh."""
    comm = MPI.COMM_WORLD
    st = _study(LumpedPortDrive())
    fhz = st.frequency_hz
    s = _synthetic_s()
    tag = f"r{comm.size}"
    two = _write(OUT_DIR / f"coding_two_{tag}.s4p", [(fhz, s), (2 * fhz, s)])
    wrong = _write(OUT_DIR / f"coding_wrong_{tag}.s4p", [(fhz * (1 + 1e-9), s)])
    ok = _write(OUT_DIR / f"coding_ok_{tag}.s4p", [(fhz, s)])
    two_port = _write(OUT_DIR / f"coding_two_port_{tag}.s2p", [(fhz, s[:2, :2])])
    for bad, word in ((wrong, "frequency_hz"), (two_port, "ports")):
        with pytest.raises(ValueError, match=word):
            st.load_result(bad)
        assert st.result is None, f"{bad.name}: refused but stored"
    # `WF-16` step 17 (§9 item 722): the refusal moved from load_result to tune.
    st_two = _study(LumpedPortDrive())
    st_two.load_result(two)
    with pytest.raises(ValueError, match="at_hz"):
        st_two.tune(k=1, comm=comm)
    res = st.load_result(ok)
    assert isinstance(res, study_mod.LoadedPortResult) and st.result is res
    rt = float(np.max(np.abs(res.s_matrix - s)))
    _p0(f"\n[WF-16 s14 W0] synthetic round trip max|dS| = {rt:.3e} (ANCHOR_TOL {ANCHOR_TOL:.0e})")
    assert rt <= ANCHOR_TOL
    for name, call in (
        ("sar", lambda: st.sar(centres=[(0.0, 0.0, 0.0)])),
        ("b1_plus", lambda: st.b1_plus()),
        ("b1_plus_at", lambda: st.b1_plus_at([(0.0, 0.0, 0.0)])),
    ):
        with pytest.raises(TypeError, match="load_result"):
            call()
        _p0(f"[WF-16 s14 W0] {name} on a loaded result: TypeError naming load_result")
    # mean_sar / temperature refuse a port Study before reaching the guard.
    with pytest.raises(TypeError):
        st.mean_sar()
    # Synthetic swap check (printed): a circulant's k = 1 read under P1 <-> P2.
    try:
        cs = st.tune(k=1, comm=comm).c_mode_farad
        st.result = study_mod.LoadedPortResult(_swap12(s), fhz, res.port_ids, "swap")
        cw = st.tune(k=1, comm=comm).c_mode_farad
        _p0(f"[WF-16 s14 W0] synthetic circulant: C_1 {cs:.9e}, P1<->P2 swapped {cw:.9e} "
            f"(rel {cw / cs - 1:+.3e}; printed)")
    except ValueError as exc:
        _p0(f"[WF-16 s14 W0] synthetic swap read not capacitive: {exc} (printed)")


@pytest.fixture(scope="module")
def s14():
    comm = MPI.COMM_WORLD
    out = {}
    s0 = _study(LumpedPortDrive())
    t0 = time.perf_counter()
    res0 = s0.solve(comm)
    out["t_solve"] = time.perf_counter() - t0
    out["S"] = np.asarray(res0.s_matrix, dtype=np.complex128).copy()
    out["tune_mem"] = s0.tune(k=1, comm=comm)
    path = res0.write_touchstone(OUT_DIR / f"s6b_pec_d1_n{comm.size}.s4p")
    out["path"] = path
    out["defs"] = res0.port_definitions
    _p0(f"\n[WF-16 s14] solve {out['t_solve']:.1f} s, cells {res0.cells}; wrote {path}")
    return out


@complex_only
def test_wf16_s14_tune_on_loaded_result(s14):
    comm = MPI.COMM_WORLD
    fresh = _study(LumpedPortDrive())
    calls = {"n": 0}

    def _boom(*a, **k):
        calls["n"] += 1
        raise AssertionError("(S0) the loaded Study solved")

    with pytest.MonkeyPatch.context() as mp:
        # (S0): every solve entry of the fresh Study raises.
        mp.setattr(fresh, "solve", _boom)
        mp.setattr(fresh, "_solve_ports", _boom)
        loaded = fresh.load_result(s14["path"])
        # (RT)
        rt = float(np.max(np.abs(loaded.s_matrix - s14["S"])))
        rt_g = comm.allreduce(rt, op=MPI.MAX)
        _p0(f"[WF-16 s14 RT] max|S_loaded - S_solved| = {rt_g:.3e} (ANCHOR_TOL {ANCHOR_TOL:.0e})")
        assert rt_g <= ANCHOR_TOL
        tl = fresh.tune(k=1, comm=comm)
    assert calls["n"] == 0
    _p0("[WF-16 s14 S0] fresh Study solve entries patched to raise; tune returned, 0 calls")
    # (CM)
    cm = s14["tune_mem"].c_mode_farad
    cl = tl.c_mode_farad
    d = abs(cl / cm - 1.0)
    d_g = comm.allreduce(d, op=MPI.MAX)
    _p0(f"[WF-16 s14 CM] C_mode loaded {cl:.12e} vs in-memory {cm:.12e}: rel {d_g:.3e} "
        f"(LEG_D0_REPRODUCTION_BAND {LEG_D0_REPRODUCTION_BAND:.0e}); "
        f"'{cl:.9e}' vs item 429 '{ITEM429_C_MODE_F:.9e}'")
    assert d_g <= LEG_D0_REPRODUCTION_BAND
    assert f"{cl:.9e}" == f"{ITEM429_C_MODE_F:.9e}"

    # Negative controls on files written back from the loaded matrix.
    fhz = fresh.frequency_hz
    s_l = np.asarray(loaded.s_matrix)

    def _tune_file(name, s):
        p = _write(OUT_DIR / f"{name}_n{comm.size}.s4p", [(fhz, s)])
        st = _study(LumpedPortDrive())
        st.load_result(p)
        return st.tune(k=1, comm=comm).c_mode_farad

    cw = _tune_file("swap12", _swap12(s_l))
    dw = abs(cw / cl - 1.0)
    moved = dw > LEG_D0_REPRODUCTION_BAND
    _p0(f"[WF-16 s14 Sw] P1<->P2 swapped: C {cw:.12e}, miss {dw:.3e} vs band "
        f"{LEG_D0_REPRODUCTION_BAND:.0e} -> {'MOVES (Sw asserted)' if moved else 'INVARIANT (Sw printed; Sw prime asserted)'}")
    if moved:
        assert dw > LEG_D0_REPRODUCTION_BAND
    c9 = _tune_file("scale0p9", 0.9 * s_l)
    d9 = abs(c9 / cl - 1.0)
    _p0(f"[WF-16 s14 Sw'] S x 0.9: C {c9:.12e}, miss {d9:.3e} vs band {LEG_D0_REPRODUCTION_BAND:.0e} "
        "(asserted)")
    assert d9 > LEG_D0_REPRODUCTION_BAND
    # (Ref) on the solved file's fixture: sar on the loaded Study names load_result.
    with pytest.raises(TypeError, match="load_result"):
        fresh.sar(centres=[(0.0, 0.0, 0.0)])
    _p0("[WF-16 s14 Ref] sar on the loaded Study: TypeError naming load_result")
