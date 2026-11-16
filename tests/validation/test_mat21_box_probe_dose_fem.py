"""`MAT-21` step 15, W2 — P39′'s dose construction on the gel box's FEM point
histories at P37's remaining F2182 probe points P3–P5 under P36's continuous
transverse quadratic source (P2, m = 4, CN N = 40 to t*), against W1's series
histories (§9 item 571, §10 chain P43; brief `docs/planning/chunks/MAT-21.md`
"§9 item 571", review ruling 2026-10-20 10:30).

Sets no band and moves none; every band is imported by name.  One CN run, read
on two fixed batches per snapshot — item 565's ``HIST_PTS`` first, then
``NEW_PTS`` — never one merged batch.  Asserted: preconditions ((Pdeg)
116 025 on both solves, ``NEW_PTS == PROBE_PTS[3:6]``, validity on both
batches, T_INF, (R36), snapshot n = 40 equal to the final field per batch,
max T < 43); (R39′) the ``HIST_PTS`` batch's dose strings equal
``P39P_FEM_DOSE``; (Rc′-h) at P0–P2; (Q-h), (Win-h), (Bd)(i) m_p ≤ ln 4·M_p
and (Bd)(ii) M_p ≤ ``INTERIOR_BAND``·max θ_s at each ``NEW_PTS`` point (the
weekly's first-order form; exact e^(ln 4·M) − 1 printed, never swapped).
Controls (rule (e)): (Rc′-h) at P3–P5 and (c-u-dose) PREDICTED; the latter
asserted only > ``INTERIOR_BAND``.
"""

from __future__ import annotations

import math
import time

import numpy as np
import pytest
import ufl
from mpi4py import MPI

from fem_em_solver.thermal.dose import cem43, running_average

from tests.complex_mode import complex_only
from tests.validation.test_mat10_sar_to_heat_sphere import _say
from tests.validation.test_mat13_transient_perfusion_sphere import INTERIOR_BAND
from tests.validation.test_mat15_coil_driven_heating_phantom import N_MAIN, T_INF
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
    W_B2,
    _q_check,
    _series_hist,
    _win_check,
)
from tests.validation.test_mat21_box_dose_fem import R36, _bd
from tests.validation.test_mat21_box_duty_dose_fem import P39P_FEM_DOSE, _cn_snap
from tests.validation.test_mat21_box_fem import _trans_ex
from tests.validation.test_mat21_box_fem_p2 import DEGREE, _dofs_pred, _rung_p
from tests.validation.test_mat21_box_probe_dose_closed_form import NEW_LABELS, NEW_PTS
from tests.validation.test_mat21_box_probe_points_fem import PROBE_PTS, _probe_vals
from tests.validation.test_mat21_box_separable_fem import M4

LN4 = math.log(4.0)
# Predictions written first (brief), printed only.
PRED_M = {0: (5e-6, 2e-4), 1: (5e-6, 2e-4), 2: (1e-5, 5e-4)}
PRED_HEADROOM = 3.0
PRED_RATIO = (0.01, 1.0)
PRED_CU = {0: (0.8, 1.0), 1: (0.8, 1.0), 2: (0.6, 1.0)}
PRED_RC = (1e6, 1e7)
PRED_ASYM = 1e-12


class _Clock15:
    """Per-phase wall clock, printed on rank 0 with flush, prefix MAT21-15."""

    def __init__(self, tag):
        self.tag = tag
        self.t = time.perf_counter()

    def lap(self, phase):
        now = time.perf_counter()
        _say(f"MAT21-15 [{self.tag}] clock {phase}: {now - self.t:.1f} s")
        self.t = now


@pytest.fixture(scope="module")
def pdose_m4():
    comm = MPI.COMM_WORLD
    t0 = time.perf_counter()
    clk = _Clock15("pdose")
    out = {"imag": 0.0}
    _say(f"\nMAT21-15 DEGREE = {DEGREE}; m = {M4}; N_MAIN = {N_MAIN}; ranks {comm.size}; "
         f"HIST_PTS {HIST_PTS}; NEW_PTS {NEW_PTS} (two fixed batches)")

    # 1. Build + uniform solve (serves (Pdeg) only).
    msh, ft, ns, n, v_h, st = _rung_p(M4, clk, DEGREE)
    out["ndof_u"] = int(st.temperature.function_space.dofmap.index_map.size_global)

    # 2. P36's quadratic CN with 41 snapshots (item 565's fixture).
    xs = ufl.SpatialCoordinate(msh)
    q_quad = Q * 3.0 * (xs[1] * xs[1] + xs[2] * xs[2]) / (B_BOX**2 + C_BOX**2)
    res = _cn_snap(msh, ft, q_quad, None)
    clk.lap(f"CN N={N_MAIN} (quadratic source, {N_MAIN + 1} snapshots)")
    out["ndof_t"] = int(res.temperature.function_space.dofmap.index_map.size_global)
    out["imag"] = max(out["imag"], _max_imag(res.temperature))

    # 3. Two fixed batches, HIST_PTS then NEW_PTS, on every snapshot (collective).
    p_h = np.asarray(HIST_PTS, dtype=float)
    p_n = np.asarray(NEW_PTS, dtype=float)
    th_h = np.empty((N_MAIN + 1, len(HIST_PTS)))
    th_n = np.empty((N_MAIN + 1, len(NEW_PTS)))
    ok_h = np.ones(len(HIST_PTS), dtype=bool)
    ok_n = np.ones(len(NEW_PTS), dtype=bool)
    for k in range(N_MAIN + 1):
        a, oa = _probe_vals(res.snapshots[k], p_h)
        b, ob = _probe_vals(res.snapshots[k], p_n)
        th_h[k], th_n[k] = a, b
        ok_h &= oa
        ok_n &= ob
    fin_h, okf_h = _probe_vals(res.temperature, p_h)
    fin_n, okf_n = _probe_vals(res.temperature, p_n)
    out["ok_h"], out["ok_n"] = ok_h & okf_h, ok_n & okf_n
    out["fin_h"], out["fin_n"] = fin_h, fin_n
    res.snapshots.clear()
    clk.lap(f"2 x ({N_MAIN + 1} + 1) owner-route batches")
    out["th_h"], out["th_n"] = th_h, th_n

    # 4. Series histories and the uniform-source control at NEW_PTS.
    out["ths_h"] = _series_hist(HIST_PTS)
    out["ths_n"] = _series_hist(NEW_PTS)
    th_u = np.zeros_like(th_n)
    for k in range(1, T_GRID.size):
        f = _trans_ex(T_GRID[k])
        th_u[k] = np.asarray(f(p_n[:, 0], p_n[:, 1], p_n[:, 2]), dtype=float)
    out["th_u"] = th_u
    clk.lap("series + uniform histories")
    _say(f"MAT21-15 fixture elapsed {time.perf_counter() - t0:.1f} s; max|Im T| {out['imag']:.3e}")
    return out


def _T(th):
    return T_INF + np.asarray(th, dtype=float)


@complex_only
def test_preconditions(pdose_m4):
    pred = _dofs_pred(M4)
    _say(f"MAT21-15 (Pdeg) P2 dofs uniform {pdose_m4['ndof_u']} / CN {pdose_m4['ndof_t']} "
         f"(predicted {pred} = 116025)")
    assert pred == 116025
    assert pdose_m4["ndof_u"] == pred and pdose_m4["ndof_t"] == pred
    assert tuple(NEW_PTS) == tuple(PROBE_PTS[3:6])
    _say(f"MAT21-15 validity HIST_PTS {pdose_m4['ok_h'].tolist()} NEW_PTS "
         f"{pdose_m4['ok_n'].tolist()}")
    assert bool(np.all(pdose_m4["ok_h"])), pdose_m4["ok_h"]
    assert bool(np.all(pdose_m4["ok_n"])), pdose_m4["ok_n"]
    assert T_INF == T_INF_0
    th_h, th_n = pdose_m4["th_h"], pdose_m4["th_n"]
    r36 = (f"{th_h[40, 2]:.6f}", f"{th_h[40, 3]:.6f}")
    same_h = bool(np.array_equal(th_h[40], pdose_m4["fin_h"]))
    same_n = bool(np.array_equal(th_n[40], pdose_m4["fin_n"]))
    _say(f"MAT21-15 (R36) theta_h(t*) centre / (0,0.9B,0) {r36} (P36 {R36}); snapshot n=40 "
         f"array_equal final: HIST_PTS {same_h}, NEW_PTS {same_n}")
    assert r36 == R36, r36
    assert same_h and same_n
    for arr in (th_h, th_n, pdose_m4["ths_h"], pdose_m4["ths_n"]):
        for j in range(arr.shape[1]):
            assert float(np.max(_T(arr[:, j]))) < 43.0


@complex_only
def test_r39p_and_rc_h_hist_batch(pdose_m4):
    th_h = pdose_m4["th_h"]
    got = tuple(f"{cem43(T_GRID, _T(th_h[:, j])):.10e}" for j in range(len(DOSE_PTS)))
    _say(f"MAT21-15 (R39') HIST_PTS batch FEM doses {got} (P39P_FEM_DOSE {P39P_FEM_DOSE}; "
         f"strings ASSERTED)")
    assert got == tuple(P39P_FEM_DOSE), got
    for j, pt in enumerate(DOSE_PTS):
        _dose, _ref, _rel, swap, miss = _q_check(_T(th_h[:, j]))
        _say(f"MAT21-15 (Rc'-h) {HIST_LABELS[j]}: R swap miss {miss:.4e} (asserted > "
             f"{DOSE_BAR:g}); factor {swap / _ref:.4e}")
        assert miss > DOSE_BAR, (pt, miss)


@complex_only
def test_q_h_new_points(pdose_m4):
    th_n = pdose_m4["th_n"]
    for j, pt in enumerate(NEW_PTS):
        dose, ref, rel, swap, miss = _q_check(_T(th_n[:, j]))
        _say(f"MAT21-15 (Q-h) {NEW_LABELS[j]}: FEM CEM43 {dose:.10e} s vs per-segment quad "
             f"{ref:.15e} s: rel {rel:.3e} (bar {DOSE_BAR:g}, ASSERTED)")
        _say(f"MAT21-15 control (Rc'-h) PREDICTED {NEW_LABELS[j]}: R swap miss {miss:.4e}; "
             f"factor {swap / ref:.4e} (predicted range {PRED_RC}; printed only, ruling 10:30 (2))")
        assert rel <= DOSE_BAR, (pt, rel)


@complex_only
def test_win_h_new_points(pdose_m4):
    th_n = pdose_m4["th_n"]
    for j, pt in enumerate(NEW_PTS):
        mx, t_mx, lo, hi, ref, rel = _win_check(th_n[:, j])
        _say(f"MAT21-15 (Win-h) {NEW_LABELS[j]}: max six-minute mean {mx:.15e} K at t_max "
             f"{t_mx:.6f} s; (ii) [{lo:.12e}, {hi:.12e}]; (iii) quad {ref:.15e}, rel {rel:.3e} "
             f"(bar {AVG_BAR:g}, ASSERTED)")
        assert T_GRID[0] + W_B2 <= t_mx <= T_GRID[-1], t_mx
        assert lo <= mx <= hi, (lo, mx, hi)
        assert rel <= AVG_BAR, (pt, rel)


@complex_only
def test_bd_dose_sensitivity_bound_new_points(pdose_m4):
    th_n, ths_n = pdose_m4["th_n"], pdose_m4["ths_n"]
    rows = []
    for j in range(len(NEW_PTS)):
        d_h, d_s, m, M, cap = _bd(th_n, ths_n, j)
        _say(f"MAT21-15 (Bd) {NEW_LABELS[j]} [ASSERTED]: FEM DOSE {d_h:.10e} s; SERIES DOSE "
             f"{d_s:.15e} s; m_p {m:.4e}; M_p {M:.4e} K (predicted {PRED_M[j]}); ln4*M_p "
             f"{LN4 * M:.4e}; exact e^(ln4 M)-1 {math.expm1(LN4 * M):.4e}; ratio m_p/(ln4 M_p) "
             f"{m / (LN4 * M):.4f} (predicted {PRED_RATIO}); INTERIOR_BAND*max theta_s "
             f"{cap:.4e} K; headroom {cap / M:.2f}x (predicted >= {PRED_HEADROOM:g}x)")
        rows.append((j, m, M, cap))
    for j, m, M, cap in rows:
        assert m <= LN4 * M, ("(Bd)(i)", NEW_LABELS[j], m, LN4 * M)
        assert M <= cap, ("(Bd)(ii)", NEW_LABELS[j], M, cap)


@complex_only
def test_c_u_dose_control_new_points(pdose_m4):
    th_n, th_u = pdose_m4["th_n"], pdose_m4["th_u"]
    for j, pt in enumerate(NEW_PTS):
        Mu = float(np.max(np.abs(th_n[:, j] - th_u[:, j])))
        r = Mu / float(np.max(th_u[:, j]))
        _say(f"MAT21-15 control (c-u-dose) PREDICTED {NEW_LABELS[j]}: M_u/max theta_u {r:.4e} "
             f"(predicted range {PRED_CU[j]}; asserted only > INTERIOR_BAND; ratio to band "
             f"{r / INTERIOR_BAND:.1f})")
        assert r > INTERIOR_BAND, (pt, r)


@complex_only
def test_printed_records(pdose_m4):
    """Printed records (no band); asserts only that they were formed."""
    th_h, th_n, ths_h, ths_n = (pdose_m4[k] for k in ("th_h", "th_n", "ths_h", "ths_n"))
    W = float(W_B2)
    dt = float(T_GRID[1] - T_GRID[0])
    d_n = [cem43(T_GRID, _T(th_n[:, j])) for j in range(len(NEW_PTS))]
    _say(f"MAT21-15 +x / -x FEM dose asymmetry d(+x)/d(-x) - 1 = {d_n[0] / d_n[1] - 1.0:+.4e} "
         f"(predicted |.| <= {PRED_ASYM:g}); max|theta_h(+x) - theta_h(-x)| "
         f"{np.max(np.abs(th_n[:, 0] - th_n[:, 1])):.3e} K")
    labs = (("P1 +y wall", th_h, ths_h, 0), ("P2 -y wall", th_h, ths_h, 1),
            ("P0 centre", th_h, ths_h, 2)) + tuple(
                (NEW_LABELS[j], th_n, ths_n, j) for j in range(len(NEW_PTS)))
    for lab, a, s, j in labs:
        _say(f"MAT21-15 SIX-POINT DOSE TABLE {lab}: FEM {cem43(T_GRID, _T(a[:, j])):.10e} s; "
             f"series {cem43(T_GRID, _T(s[:, j])):.15e} s")
    for j in range(len(NEW_PTS)):
        a_h, mx_h, t_h = running_average(T_GRID, th_n[:, j], window=W)
        a_s, mx_s, t_s = running_average(T_GRID, ths_n[:, j], window=W)
        _m, mx_sh, t_sh = running_average(T_GRID, th_n[:, j], window=W - dt)
        good = ~np.isnan(a_h)
        da = float(np.max(np.abs(a_h[good] - a_s[good])))
        dth = float(np.max(np.abs(th_n[:, j] - ths_n[:, j])))
        rhs = dth * (1.0 + 1e-12) + D_BAND * float(np.max(np.abs(a_s[good])))
        _say(f"MAT21-15 {NEW_LABELS[j]}: six-minute peak FEM {mx_h:.9f} K at {t_h:.3f} s vs "
             f"series {mx_s:.9f} K at {t_s:.3f} s; (Sh) W - dt: {mx_sh:.9f} K at {t_sh:.3f} s; "
             f"(FEM-c-win) max|a_h - a_s| {da:.4e} vs max|dtheta|(1+1e-12) + D_BAND max|a_s| "
             f"{rhs:.4e} (printed only; holds {da <= rhs})")
    assert all(np.isfinite(x) for x in d_n)
