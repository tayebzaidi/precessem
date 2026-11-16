"""`MAT-21` step 12, W2 — `MAT-20`'s dose and six-minute window on the gel box's
FEM point histories under P36's transverse quadratic source (P2, m = 4, CN
N = 40 to t*), against W1's series histories (§9 item 565, §10 chain P39′;
brief `docs/planning/chunks/MAT-21.md` "§9 item 565", review ruling
2026-10-20 03:00).

Sets no band and moves none; every band is imported by name.  Asserted:
preconditions ((Pdeg) 116 025, point sets, validity, T_INF, (R36) P36's
printed digits, snapshot n = 40 equal to the final field, max T < 43);
(Q-h), (Win-h) at ``DOSE_BAR`` / ``AVG_BAR``; (Bd)(i) m_p ≤ ln 4 · M_p and
(Bd)(ii) M_p ≤ ``INTERIOR_BAND`` · max θ_s at each ``DOSE_PTS`` point (the
weekly's first-order form; the exact e^(ln 4·M) − 1 printed beside, never
swapped).  Controls (rule (e)): (Rc′-h) asserted; (c-u-dose) PREDICTED,
measured value printed beside its range, only > ``INTERIOR_BAND`` asserted.
"""

from __future__ import annotations

import math
import time

import numpy as np
import pytest
import ufl
from mpi4py import MPI

from fem_em_solver.thermal import solve_transient_conduction
from fem_em_solver.thermal.dose import cem43, running_average

from tests.complex_mode import complex_only
from tests.validation.test_mat10_sar_to_heat_sphere import ROBIN_TAG, _say
from tests.validation.test_mat12_transient_sar_heating_sphere import RHO_C
from tests.validation.test_mat13_transient_perfusion_sphere import INTERIOR_BAND
from tests.validation.test_mat15_coil_driven_heating_phantom import (
    H_GEL,
    K_GEL,
    N_MAIN,
    T_INF,
)
from tests.validation.test_mat15_finite_cylinder_fem import _max_imag
from tests.validation.test_mat17_duty_cycle_coil_driven import D_BAND
from tests.validation.test_mat21_box_closed_form import B_BOX, C_BOX, Q
from tests.validation.test_mat21_box_closed_form import T_INF as T_INF_0
from tests.validation.test_mat21_box_dose_closed_form import (
    AVG_BAR,
    DOSE_BAR,
    DOSE_PTS,
    HIST_LABELS,
    HIST_PTS,
    T_GRID,
    T_STAR,
    W_B2,
    _q_check,
    _series_hist,
    _win_check,
)
from tests.validation.test_mat21_box_fem import _trans_ex
from tests.validation.test_mat21_box_fem_p2 import DEGREE, _dofs_pred, _rung_p
from tests.validation.test_mat21_box_probe_points_fem import ADJ_PTS, PROBE_PTS, _probe_vals
from tests.validation.test_mat21_box_separable_fem import M4

LN4 = math.log(4.0)
R36 = ("0.063442", "2.893351")  # P36, 20261016T171904Z_MAT-21-step9-W2-quadratic-fem.log:83
# Predictions written first (brief), printed only.
PRED_M = {"wall": (1e-5, 2e-3), "centre": 2.4e-5}
PRED_MP = {"wall": (1e-6, 3e-3), "centre": (1e-7, 1e-3)}
PRED_RATIO = (0.05, 1.0)
PRED_CU = (0.2, 1.0)


class _Clock12:
    """Per-phase wall clock, printed on rank 0 with flush, prefix MAT21-12."""

    def __init__(self, tag):
        self.tag = tag
        self.t = time.perf_counter()

    def lap(self, phase):
        now = time.perf_counter()
        _say(f"MAT21-12 [{self.tag}] clock {phase}: {now - self.t:.1f} s")
        self.t = now


@pytest.fixture(scope="module")
def dose_m4():
    comm = MPI.COMM_WORLD
    t0 = time.perf_counter()
    clk = _Clock12("dose")
    out = {"imag": 0.0}
    _say(f"\nMAT21-12 DEGREE = {DEGREE}; m = {M4}; N_MAIN = {N_MAIN}; ranks {comm.size}; "
         f"HIST_PTS {HIST_PTS}")

    # 1. Build + uniform solve (serves (Pdeg) only).
    msh, ft, ns, n, v_h, st = _rung_p(M4, clk, DEGREE)
    out["ndof_u"] = int(st.temperature.function_space.dofmap.index_map.size_global)

    # 2. P36's quadratic CN (step 9 :180-183) with every step kept.
    xs = ufl.SpatialCoordinate(msh)
    q_quad = Q * 3.0 * (xs[1] * xs[1] + xs[2] * xs[2]) / (B_BOX**2 + C_BOX**2)
    res = solve_transient_conduction(
        msh, ft, k=K_GEL, rho_c=RHO_C, q=q_quad, h=H_GEL, t_inf=T_INF, t_init=T_INF,
        robin_tag=ROBIN_TAG, dt=T_STAR / N_MAIN, n_steps=N_MAIN, theta=0.5,
        startup_half_steps=2, degree=DEGREE, snapshot_steps=tuple(range(N_MAIN + 1)))
    clk.lap(f"CN N={N_MAIN} (quadratic source, {N_MAIN + 1} snapshots)")
    out["ndof_t"] = int(res.temperature.function_space.dofmap.index_map.size_global)
    out["imag"] = max(out["imag"], _max_imag(res.temperature))

    # 3. One batch, one order, on every snapshot (collective).
    pts = np.asarray(HIST_PTS, dtype=float)
    th_h = np.empty((N_MAIN + 1, len(HIST_PTS)))
    ok_all = np.ones(len(HIST_PTS), dtype=bool)
    for k in range(N_MAIN + 1):
        th, ok = _probe_vals(res.snapshots[k], pts)
        th_h[k] = th
        ok_all &= ok
    th_fin, ok_fin = _probe_vals(res.temperature, pts)
    out["ok"] = ok_all & ok_fin
    out["th_fin"] = th_fin
    res.snapshots.clear()
    clk.lap(f"{N_MAIN + 1} + 1 owner-route batches")
    out["th_h"] = th_h

    # 4. Series histories (quadratic and uniform source).
    out["th_s"] = _series_hist(HIST_PTS)
    th_u = np.zeros_like(th_h)
    for k in range(1, T_GRID.size):
        f = _trans_ex(T_GRID[k])
        th_u[k] = np.asarray(f(pts[:, 0], pts[:, 1], pts[:, 2]), dtype=float)
    out["th_u"] = th_u
    clk.lap("series + uniform histories")
    _say(f"MAT21-12 fixture elapsed {time.perf_counter() - t0:.1f} s; max|Im T| {out['imag']:.3e}")
    return out


def _T(th):
    return T_INF + np.asarray(th, dtype=float)


@complex_only
def test_preconditions(dose_m4):
    pred = _dofs_pred(M4)
    _say(f"MAT21-12 (Pdeg) P2 dofs uniform {dose_m4['ndof_u']} / CN {dose_m4['ndof_t']} "
         f"(predicted {pred} = 116025)")
    assert pred == 116025
    assert dose_m4["ndof_u"] == pred and dose_m4["ndof_t"] == pred
    assert set(DOSE_PTS) <= set(PROBE_PTS) and HIST_PTS[3] == ADJ_PTS[0]
    assert bool(np.all(dose_m4["ok"])), dose_m4["ok"]
    assert T_INF == T_INF_0
    th_h = dose_m4["th_h"]
    r36 = (f"{th_h[40, 2]:.6f}", f"{th_h[40, 3]:.6f}")
    same = bool(np.array_equal(th_h[40], dose_m4["th_fin"]))
    _say(f"MAT21-12 (R36) theta_h(t*) centre / (0,0.9B,0) {r36} (P36 {R36}); snapshot n=40 "
         f"array_equal final {same}")
    assert r36 == R36, r36
    assert same
    for j in range(len(HIST_PTS)):
        for arr in (th_h, dose_m4["th_s"]):
            assert float(np.max(_T(arr[:, j]))) < 43.0


@complex_only
def test_q_h_and_rc_h(dose_m4):
    th_h = dose_m4["th_h"]
    for j, pt in enumerate(DOSE_PTS):
        dose, ref, rel, swap, miss = _q_check(_T(th_h[:, j]))
        _say(f"MAT21-12 (Q-h) {HIST_LABELS[j]}: FEM CEM43 {dose:.15e} s vs per-segment quad "
             f"{ref:.15e} s: rel {rel:.3e} (bar {DOSE_BAR:g}, ASSERTED)")
        _say(f"MAT21-12 (Rc'-h) {HIST_LABELS[j]}: R swap miss {miss:.4e} (asserted > {DOSE_BAR:g}); "
             f"factor {swap / ref:.4e}")
        assert rel <= DOSE_BAR, (pt, rel)
        assert miss > DOSE_BAR, (pt, miss)


@complex_only
def test_win_h(dose_m4):
    th_h = dose_m4["th_h"]
    for j, pt in enumerate(DOSE_PTS):
        mx, t_mx, lo, hi, ref, rel = _win_check(th_h[:, j])
        _say(f"MAT21-12 (Win-h) {HIST_LABELS[j]}: max six-minute mean {mx:.15e} K at t_max "
             f"{t_mx:.6f} s; (ii) [{lo:.12e}, {hi:.12e}]; (iii) quad {ref:.15e}, rel {rel:.3e} "
             f"(bar {AVG_BAR:g}, ASSERTED)")
        assert T_GRID[0] + W_B2 <= t_mx <= T_GRID[-1], t_mx
        assert lo <= mx <= hi, (lo, mx, hi)
        assert rel <= AVG_BAR, (pt, rel)


def _bd(th_h, th_s, j):
    d_h = cem43(T_GRID, _T(th_h[:, j]))
    d_s = cem43(T_GRID, _T(th_s[:, j]))
    m = abs(d_h / d_s - 1.0)
    M = float(np.max(np.abs(th_h[:, j] - th_s[:, j])))
    cap = INTERIOR_BAND * float(np.max(th_s[:, j]))
    return d_h, d_s, m, M, cap


@complex_only
def test_bd_dose_sensitivity_bound(dose_m4):
    th_h, th_s = dose_m4["th_h"], dose_m4["th_s"]
    rows = []
    for j in range(len(HIST_PTS)):
        d_h, d_s, m, M, cap = _bd(th_h, th_s, j)
        kind = "centre" if j == 2 else "wall"
        tag = "ASSERTED" if j < 3 else "printed (interior point)"
        _say(f"MAT21-12 (Bd) {HIST_LABELS[j]} [{tag}]: FEM DOSE {d_h:.10e} s; SERIES DOSE "
             f"{d_s:.15e} s; m_p {m:.4e} (predicted {PRED_MP.get(kind)}); M_p {M:.4e} K "
             f"(predicted {PRED_M.get(kind)}); ln4*M_p {LN4 * M:.4e}; exact e^(ln4 M)-1 "
             f"{math.expm1(LN4 * M):.4e}; ratio m_p/(ln4 M_p) {m / (LN4 * M):.4f} (predicted "
             f"{PRED_RATIO}); INTERIOR_BAND*max theta_s {cap:.4e} K; headroom M_p/cap "
             f"{M / cap:.4e} ({cap / M:.1f}x)")
        rows.append((j, m, M, cap))
    for j, m, M, cap in rows[:3]:
        assert m <= LN4 * M, ("(Bd)(i)", HIST_LABELS[j], m, LN4 * M)
        assert M <= cap, ("(Bd)(ii)", HIST_LABELS[j], M, cap)


@complex_only
def test_c_u_dose_control(dose_m4):
    th_h, th_u = dose_m4["th_h"], dose_m4["th_u"]
    for j, pt in enumerate(DOSE_PTS):
        Mu = float(np.max(np.abs(th_h[:, j] - th_u[:, j])))
        r = Mu / float(np.max(th_u[:, j]))
        _say(f"MAT21-12 control (c-u-dose) PREDICTED {HIST_LABELS[j]}: M_u/max theta_u {r:.4e} "
             f"(predicted range {PRED_CU}; asserted only > INTERIOR_BAND; ratio to band "
             f"{r / INTERIOR_BAND:.1f})")
        assert r > INTERIOR_BAND, (pt, r)


@complex_only
def test_printed_records(dose_m4):
    """Printed records (no band); asserts only that they were formed."""
    th_h, th_s = dose_m4["th_h"], dose_m4["th_s"]
    W = float(W_B2)
    dt = float(T_GRID[1] - T_GRID[0])
    d = [cem43(T_GRID, _T(th_h[:, j])) for j in range(len(HIST_PTS))]
    _say(f"MAT21-12 +y / -y FEM dose asymmetry d(+y)/d(-y) - 1 = {d[0] / d[1] - 1.0:+.4e}; "
         f"max|theta_h(+y) - theta_h(-y)| {np.max(np.abs(th_h[:, 0] - th_h[:, 1])):.3e} K")
    for j in range(len(HIST_PTS)):
        a_h, mx_h, t_h = running_average(T_GRID, th_h[:, j], window=W)
        a_s, mx_s, t_s = running_average(T_GRID, th_s[:, j], window=W)
        _m, mx_sh, t_sh = running_average(T_GRID, th_h[:, j], window=W - dt)
        good = ~np.isnan(a_h)
        da = float(np.max(np.abs(a_h[good] - a_s[good])))
        dth = float(np.max(np.abs(th_h[:, j] - th_s[:, j])))
        rhs = dth * (1.0 + 1e-12) + D_BAND * float(np.max(np.abs(a_s[good])))
        _say(f"MAT21-12 {HIST_LABELS[j]}: six-minute peak FEM {mx_h:.9f} K at {t_h:.3f} s vs "
             f"series {mx_s:.9f} K at {t_s:.3f} s; (Sh) W - dt: {mx_sh:.9f} K at {t_sh:.3f} s; "
             f"(FEM-c-win) max|a_h - a_s| {da:.4e} vs max|dtheta|(1+1e-12) + D_BAND max|a_s| "
             f"{rhs:.4e} (printed only; holds {da <= rhs})")
    assert all(np.isfinite(x) for x in d)
