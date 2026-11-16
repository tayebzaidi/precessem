"""`IMP-8` step 15 (§9 item 856, §10 I73): the self-heated rod on I55's box.

Brief: docs/planning/chunks/IMP-8.md, "§9 item 856", under the "Review ruling
2026-11-06 03:00" paragraph.  I55's graded gel box (``box`` fixture, imported
unchanged), real build, degree 2, CN to t* = 900 s, I55's uniform rod source
(``_source_and_s0`` → ``_rod_source``: DG0 indicator of midpoint r ≤
``ROD_RADIUS`` ∧ |z| ≤ L/2).  The body cells are exactly the source cells —
taken as the positive DG0 dofs of the returned q (owned + ghost, q is
scatter_forward'ed) — and are given (``ROD_K``, ``ROD_RHO_C``, both nominal)
through DG0 ``k`` / ``rho_c``.  Nothing is subtracted (no ball source).

Arm knob ``FEM_EM_IMP8_I73_ARM`` = ``ti`` (default) / ``gel``.  Arm ``gel`` uses
the same DG0 objects with gel values everywhere: the anchor (I55's five θ to one
unit of the sixth significant digit) and I55's (I) gate at ``INTERIOR_BAND``.
Arm ``ti``'s negative control: θ(2, 0) mm moves off I55's 1.772315e+01 by more
than 1e-9 relative.

Asserted (bands imported): (E) ``max_balance_residual`` ≤ ``E_BAND``; (S0)
|∫q − P|/P ≤ ``B_BAND``; cells = ``W0C_CELLS``; body cells = 6 495 (I55's W0c).
Everything else is a 🟡 record, printed, no band.
"""

from __future__ import annotations

import math
import os
import time

import numpy as np
from dolfinx import fem
from mpi4py import MPI

import pytest

from fem_em_solver.thermal import solve_transient_conduction

from tests.mesh.test_birdcage_implant_rod import ROD_RADIUS
from tests.validation.test_imp8_composite_cylinder_closed_form import (
    ROD_K,  # nominal
    ROD_RHO_C,  # nominal
)
from tests.validation.test_imp8_rod_line_source_closed_form import L, MM, T_PRINT, _phi
from tests.validation.test_imp8_rod_line_source_fem import (  # noqa: F401  (box: fixture)
    DEGREE,
    GATE_PTS,
    P,
    PRINT_PTS,
    QP,
    T_INF,
    W0C_CELLS,
    _need,
    _rod_source,  # noqa: F401  (called through _source_and_s0)
    _source_and_s0,
    _theta_pts,
    box,
)
from tests.validation.test_mat10_sar_to_heat_sphere import B_BAND, ROBIN_TAG, _say
from tests.validation.test_mat12_transient_sar_heating_sphere import E_BAND, RHO_C
from tests.validation.test_mat13_transient_perfusion_sphere import INTERIOR_BAND
from tests.validation.test_mat15_coil_driven_heating_phantom import H_GEL, K_GEL, N_MAIN

ARM_ENV = "FEM_EM_IMP8_I73_ARM"
ARMS = ("ti", "gel")
#: I55's selected (source) cell count, W0c.
BODY_CELLS = 6_495
#: I55's printed θ, K/W, t* = 900 s (`20261104T030134Z_IMP-8-I55-W2-item803.log:142–146`).
I55_THETA = {
    (10.0, 0.0): 4.651377e00,
    (20.0, 0.0): 1.105475e00,
    (10.0, 15.0): 2.732755e00,
    (20.0, 15.0): 7.045591e-01,
    (2.0, 0.0): 1.772315e01,
}
CONTROL_FLOOR = 1e-9  # the round-off floor (I55's (E) is 2.098e-14)
READ_PTS = GATE_PTS + PRINT_PTS + ((2.0, 15.0),)


def _arm():
    raw = os.environ.get(ARM_ENV, "").strip()
    if raw == "":
        return "ti"
    if raw in ARMS:
        return raw
    raise ValueError(f"{ARM_ENV} must be unset or one of {ARMS}, got {raw!r}")


def test_k0_arm_knob(monkeypatch):
    monkeypatch.delenv(ARM_ENV, raising=False)
    assert _arm() == "ti"
    for a in ARMS:
        monkeypatch.setenv(ARM_ENV, a)
        assert _arm() == a
    for bad in ("Ti", "hole", "1"):
        monkeypatch.setenv(ARM_ENV, bad)
        with pytest.raises(ValueError):
            _arm()


def _coefficients(msh, q, arm):
    """DG0 k, rho_c on q's space; the body is q's positive dofs (owned + ghost)."""
    W = q.function_space
    body = np.flatnonzero(np.real(q.x.array) > 0.0)
    n_own = W.dofmap.index_map.size_local * W.dofmap.index_map_bs
    n_body = msh.comm.allreduce(int(np.count_nonzero(body < n_own)), op=MPI.SUM)
    k = fem.Function(W)
    rho_c = fem.Function(W)
    k.x.array[:] = K_GEL
    rho_c.x.array[:] = RHO_C
    if arm == "ti":
        k.x.array[body] = ROD_K
        rho_c.x.array[body] = ROD_RHO_C
    k.x.scatter_forward()
    rho_c.x.scatter_forward()
    bad = int(np.count_nonzero(~np.isfinite(k.x.array)) + np.count_nonzero(~np.isfinite(rho_c.x.array)))
    return k, rho_c, n_body, msh.comm.allreduce(bad, op=MPI.SUM)


@pytest.fixture(scope="module")
def body_run(box):
    _need("transient")
    arm = _arm()
    msh, ct, ft, n = box
    q, out = _source_and_s0(msh)
    k, rho_c, n_body, n_bad = _coefficients(msh, q, arm)
    out.update({"arm": arm, "n": n, "n_body": n_body, "n_bad": n_bad})
    v_ratio = (out["a_eff"] / ROD_RADIUS) ** 2  # V_body / (pi a^2 L)
    _say(f"[imp8-i73] arm {arm}; ROD_K {ROD_K:g} W/(m K), ROD_RHO_C {ROD_RHO_C:.6g} J/(m3 K) "
         f"(nominal); K_GEL {K_GEL:g}, RHO_C {RHO_C:.6g}")
    _say(f"[imp8-i73] body cells (positive q dofs, owned, allreduced) {n_body} (I55 W0c "
         f"{BODY_CELLS}); source-selected cells {out['n_sel']}; V_body/(pi a^2 L) {v_ratio:.6f}")
    dt = T_PRINT / N_MAIN
    t1 = time.perf_counter()
    res = solve_transient_conduction(
        msh, ft, k=k, rho_c=rho_c, q=q, h=H_GEL, t_inf=T_INF, t_init=T_INF,
        robin_tag=ROBIN_TAG, dt=dt, n_steps=N_MAIN, theta=0.5, startup_half_steps=2,
        degree=DEGREE)
    out["t_call"] = time.perf_counter() - t1
    out["E"] = float(res.max_balance_residual)
    _say(f"[imp8-i73] CN N={N_MAIN} dt={dt:g} s: t_call {out['t_call']:.1f} s; "
         f"n_factorizations {res.n_factorizations}; ranks {msh.comm.size}")
    th = _theta_pts(res.temperature, READ_PTS)  # one batch
    out["theta"] = dict(zip(READ_PTS, th))
    out["I"] = []
    for (r, z), v in zip(READ_PTS, th):
        cf = QP * _phi(r * MM, z * MM, T_PRINT)
        rel = abs(v / cf - 1.0)
        if (r, z) in GATE_PTS:
            out["I"].append(rel)
        ref = I55_THETA.get((r, z))
        if ref is not None:
            tol = 10.0 ** (math.floor(math.log10(abs(ref))) - 5)
            cmp = (f"I55 {ref:.6e}; theta/theta_55 {v / ref:.7f} (|diff| {abs(v - ref):.3e}, "
                   f"6th-digit unit {tol:.0e})")
        else:
            cmp = "I55 not read here (new point)"
        lab = ("(I) gated on arm gel" if arm == "gel" and (r, z) in GATE_PTS
               else "not a gate: the comparand has no rod" if arm == "ti" else "printed only")
        _say(f"[imp8-i73] ({r:g}, {z:g}) mm t*={T_PRINT:g} s arm {arm}: theta {v:.6e} K/W; {cmp}; "
             f"vs q'/(4 pi k) Phi {cf:.6e}: rel {rel:.3e} ({lab}, band {INTERIOR_BAND:.0e})")
    ax = out["theta"][(2.0, 15.0)] / out["theta"][(2.0, 0.0)]
    move = out["theta"][(2.0, 0.0)] / I55_THETA[(2.0, 0.0)] - 1.0
    out["move"] = move
    _say(f"[imp8-i73] axial ratio theta(2,15)/theta(2,0) arm {arm}: {ax:.6f}")
    _say(f"[imp8-i73] (2, 0) mm move theta/theta_55 - 1 = {move:.4e} (control floor "
         f"{CONTROL_FLOOR:.0e}, asserted on arm ti; predicted sign (not asserted): theta(2,0) "
         f"falls and the axial ratio moves toward 1 on arm ti)")
    _say(f"[imp8-i73] (E) {out['E']:.3e} (band {E_BAND:.0e}); (S0) {out['S0']:.3e} (band "
         f"{B_BAND:.0e}); cells {n} (W0c {W0C_CELLS}); non-finite coefficient dofs {n_bad}")
    return out


def test_w_preconditions(body_run):
    assert body_run["n"] == W0C_CELLS, body_run["n"]
    assert body_run["n_body"] == BODY_CELLS, body_run["n_body"]
    assert body_run["n_bad"] == 0


def test_w_e_balance(body_run):
    assert body_run["E"] <= E_BAND, body_run["E"]


def test_w_s0_source_power(body_run):
    assert body_run["S0"] <= B_BAND, body_run["S0"]


def test_w2_gel_reproduces_i55(body_run):
    if body_run["arm"] != "gel":
        pytest.skip("W2 only (arm gel): the anchor")
    for pt, ref in I55_THETA.items():
        v = body_run["theta"][pt]
        tol = 10.0 ** (math.floor(math.log10(abs(ref))) - 5)
        assert abs(v - ref) <= tol, (pt, v, ref, tol)


def test_w2_gel_interior_against_line_source(body_run):
    if body_run["arm"] != "gel":
        pytest.skip("W2 only (arm gel): I55's (I) re-asserted")
    assert max(body_run["I"]) <= INTERIOR_BAND, body_run["I"]


def test_w1_ti_control_moves(body_run):
    if body_run["arm"] != "ti":
        pytest.skip("W1 only (arm ti): the negative control")
    assert abs(body_run["move"]) > CONTROL_FLOOR, body_run["move"]
