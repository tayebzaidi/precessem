"""`MAT-21` step 14, W2 — `MAT-20`'s dose and six-minute window on the gel box's
FEM point histories at a second duty fraction, D = 0.25 at P38's period (P36's
transverse quadratic source; P2, m = 4, CN N = 40 to t*) (§9 item 570, §10
chain P42; brief `docs/planning/chunks/MAT-21.md` "§9 item 570", review ruling
2026-10-20 10:30).

Sets no band and moves none; every band is imported by name.  Two CN runs
(default, D = 0.25), each read on item 565's one batch ``HIST_PTS`` on all 41
snapshots.  Asserted: preconditions ((Pdeg) on three solves, (S0) census +
schedule string, validity, (E) ≤ ``E_BAND``, ``T_INF``, max T < 43); (R39′)
the default run's dose strings; (Q-h25) / (Win-h25) at peak and equal mean
power; (Bd25)(i) m_p ≤ ln 4·M_p and (ii) M_p ≤ ``INTERIOR_BAND``·max θ_s
(first-order form, ruling 2026-10-20 03:00 (1)); (Rc′-h25) the R swap misses;
(c-50-dose) the D = 0.5 comparand misses by > ``INTERIOR_BAND`` (predicted).
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
from tests.validation.test_mat12_transient_sar_heating_sphere import E_BAND
from tests.validation.test_mat13_transient_perfusion_sphere import INTERIOR_BAND
from tests.validation.test_mat15_coil_driven_heating_phantom import N_MAIN, T_INF
from tests.validation.test_mat15_finite_cylinder_fem import _max_imag
from tests.validation.test_mat17_duty_cycle_fem import _census
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
from tests.validation.test_mat21_box_dose_fem import _bd
from tests.validation.test_mat21_box_duty25_dose_closed_form import (
    DUTY25,
    P41_EQ_RATIO,
    SCHEDULE_25,
    T_SW25,
    _chi25,
)
from tests.validation.test_mat21_box_duty_closed_form import DT_BOX
from tests.validation.test_mat21_box_duty_dose_closed_form import _duty_hist
from tests.validation.test_mat21_box_duty_dose_fem import LN4, P39P_FEM_DOSE, _cn_snap
from tests.validation.test_mat21_box_fem_p2 import DEGREE, _dofs_pred, _rung_p
from tests.validation.test_mat21_box_probe_points_fem import _probe_vals
from tests.validation.test_mat21_box_separable_fem import M4

# P41's W2 FEM D = 0.5 equal-mean ratios (printed beside), from
# docs/testing/logs/20261020T141240Z_MAT-21-step13-W2-duty-dose-fem.log:104-107.
P41_FEM_EQ = ("1.185888", "1.185888", "1.004496", "1.212846")
PRED_M_WALL = (3e-4, 3e-3)
PRED_M_CENTRE = (5e-6, 5e-5)
PRED_HEAD = {"wall": (3, 15), "centre": (5, 30)}
PRED_RATIO = (0.02, 1.0)
PRED_C50 = (0.3, 2.0)


class _Clock14:
    """Per-phase wall clock, printed on rank 0 with flush, prefix MAT21-14."""

    def __init__(self, tag):
        self.tag = tag
        self.t = time.perf_counter()

    def lap(self, phase):
        now = time.perf_counter()
        _say(f"MAT21-14 [{self.tag}] clock {phase}: {now - self.t:.1f} s")
        self.t = now


@pytest.fixture(scope="module")
def d25_m4():
    comm = MPI.COMM_WORLD
    t0 = time.perf_counter()
    clk = _Clock14("d25")
    out = {"imag": 0.0, "E": {}, "ndof": {}, "th": {}, "ok": {}}
    _say(f"\nMAT21-14 DEGREE = {DEGREE}; m = {M4}; N_MAIN = {N_MAIN}; ranks {comm.size}; "
         f"HIST_PTS {HIST_PTS}")

    msh, ft, ns, n, v_h, st = _rung_p(M4, clk, DEGREE)
    out["ndof"]["uniform"] = int(st.temperature.function_space.dofmap.index_map.size_global)

    xs = ufl.SpatialCoordinate(msh)
    q_quad = Q * 3.0 * (xs[1] * xs[1] + xs[2] * xs[2]) / (B_BOX**2 + C_BOX**2)
    pts = np.asarray(HIST_PTS, dtype=float)

    for label, schedule in (("default", None), ("d25", _chi25)):
        r = _cn_snap(msh, ft, q_quad, schedule)
        clk.lap(f"CN N={N_MAIN} {label} ({N_MAIN + 1} snapshots)")
        out["ndof"][label] = int(r.temperature.function_space.dofmap.index_map.size_global)
        out["imag"] = max(out["imag"], _max_imag(r.temperature))
        out["E"][label] = float(r.max_balance_residual)
        if label == "d25":
            out["sched"] = "".join(str(int(x)) for x in r.source_multipliers)
            _census(r, DT_BOX, 2, DUTY25 * T_STAR, "quad duty 25")
            out["census_ok"] = True
        th_h = np.empty((N_MAIN + 1, len(HIST_PTS)))
        ok_all = np.ones(len(HIST_PTS), dtype=bool)
        for k in range(N_MAIN + 1):
            th, ok = _probe_vals(r.snapshots[k], pts)
            th_h[k] = th
            ok_all &= ok
        r.snapshots.clear()
        del r
        out["th"][label] = th_h
        out["ok"][label] = ok_all
        clk.lap(f"{label}: {N_MAIN + 1} owner-route batches")

    out["th_s"] = _series_hist(HIST_PTS)
    out["th_s25"] = _duty_hist(HIST_PTS, switch_times=T_SW25)
    out["th_s50"] = _duty_hist(HIST_PTS)
    clk.lap("series histories (continuous + D = 0.25 + D = 0.5)")
    _say(f"MAT21-14 fixture elapsed {time.perf_counter() - t0:.1f} s; max|Im T| {out['imag']:.3e}")
    return out


def _T(th):
    return T_INF + np.asarray(th, dtype=float)


@complex_only
def test_preconditions(d25_m4):
    pred = _dofs_pred(M4)
    _say(f"MAT21-14 (Pdeg) P2 dofs {d25_m4['ndof']} (predicted {pred} = 116025)")
    assert pred == 116025
    assert len(d25_m4["ndof"]) == 3
    for key, v in d25_m4["ndof"].items():
        assert v == pred, f"(Pdeg) {key}: {v}"
    _say(f"MAT21-14 (S0) schedule {d25_m4['sched']} (must equal {SCHEDULE_25})")
    assert d25_m4["census_ok"] and d25_m4["sched"] == SCHEDULE_25
    for key, ok in d25_m4["ok"].items():
        assert bool(np.all(ok)), (key, ok)
    _say("MAT21-14 (E) max balance residual: " + ", ".join(
        f"{k_} {v:.3e}" for k_, v in d25_m4["E"].items()) + f" (bound {E_BAND:.0e}; predicted 2.3e-13)")
    assert len(d25_m4["E"]) == 2
    for key, v in d25_m4["E"].items():
        assert v <= E_BAND, (key, v)
    assert T_INF == T_INF_0
    for key, arr in d25_m4["th"].items():
        for j in range(len(HIST_PTS)):
            assert float(np.max(_T(arr[:, j]))) < 43.0, (key, j)
            assert float(np.max(_T(arr[:, j] / DUTY25))) < 43.0, (key, j)


@complex_only
def test_r39p_default_reproduces_p39p_fem_dose(d25_m4):
    th = d25_m4["th"]["default"]
    got = tuple(f"{cem43(T_GRID, _T(th[:, j])):.10e}" for j in range(len(DOSE_PTS)))
    _say(f"MAT21-14 (R39') default FEM doses {got} (P39' {P39P_FEM_DOSE})")
    assert got == P39P_FEM_DOSE, got


@complex_only
def test_q_h25_win_h25_and_rc_h25(d25_m4):
    th_d = d25_m4["th"]["d25"]
    for scale, lab in ((1.0, "peak"), (1.0 / DUTY25, "equal-mean")):
        for j, pt in enumerate(DOSE_PTS):
            th = th_d[:, j] * scale
            dose, ref, rel, swap, miss = _q_check(_T(th))
            mx, t_mx, lo, hi, wref, wrel = _win_check(th)
            _say(f"MAT21-14 (Q-h25) {lab} {HIST_LABELS[j]}: FEM CEM43 {dose:.15e} s vs quad "
                 f"{ref:.15e} s rel {rel:.3e} (bar {DOSE_BAR:g}, ASSERTED); (Rc'-h25) swap miss "
                 f"{miss:.4e} (asserted > {DOSE_BAR:g})")
            _say(f"MAT21-14 (Win-h25) {lab} {HIST_LABELS[j]}: max six-minute mean {mx:.15e} K at "
                 f"{t_mx:.4f} s; (ii) [{lo:.12e}, {hi:.12e}]; (iii) quad {wref:.15e} rel {wrel:.3e} "
                 f"(bar {AVG_BAR:g}, ASSERTED)")
            assert rel <= DOSE_BAR, (lab, pt, rel)
            assert miss > DOSE_BAR, (lab, pt, miss)
            assert T_GRID[0] + W_B2 <= t_mx <= T_GRID[-1], (lab, t_mx)
            assert lo <= mx <= hi, (lab, lo, mx, hi)
            assert wrel <= AVG_BAR, (lab, pt, wrel)


@complex_only
def test_bd25_dose_sensitivity_bound(d25_m4):
    th_h, th_s = d25_m4["th"]["d25"], d25_m4["th_s25"]
    for scale, lab in ((1.0, "peak"), (1.0 / DUTY25, "equal-mean")):
        h, s = th_h * scale, th_s * scale
        for j, pt in enumerate(DOSE_PTS):
            d_h, d_s, m, M, cap = _bd(h, s, j)
            kind = "centre" if j == 2 else "wall"
            pm = PRED_M_CENTRE if kind == "centre" else PRED_M_WALL
            _say(f"MAT21-14 (Bd25) {lab} {HIST_LABELS[j]}: D_h {d_h:.10e} s, D_s {d_s:.15e} s; "
                 f"(i) m_p {m:.4e} <= ln4*M_p {LN4 * M:.4e} (exact e^(ln4 M)-1 "
                 f"{math.expm1(LN4 * M):.4e}, printed only); ratio {m / (LN4 * M):.4f} (predicted "
                 f"{PRED_RATIO}); (ii) M_p {M:.4e} K <= INTERIOR_BAND*max theta_s {cap:.4e} K, "
                 f"headroom {cap / M:.2f}x (predicted M {pm} K at peak, headroom {PRED_HEAD[kind]})")
            assert m <= LN4 * M, ("(Bd25)(i)", lab, pt, m, LN4 * M)
            assert M <= cap, ("(Bd25)(ii)", lab, pt, M, cap)


@complex_only
def test_c50_dose_wrong_schedule_control(d25_m4):
    th_h, th_s50 = d25_m4["th"]["d25"], d25_m4["th_s50"]
    th_s25 = d25_m4["th_s25"]
    for j, pt in enumerate(DOSE_PTS):
        Mw = float(np.max(np.abs(th_h[:, j] - th_s50[:, j])))
        r = Mw / float(np.max(th_s25[:, j]))
        _say(f"MAT21-14 (c-50-dose) {HIST_LABELS[j]}: M_w / max theta_s25 {r:.4f} (asserted > "
             f"INTERIOR_BAND {INTERIOR_BAND:g}; PREDICTED {PRED_C50})")
        assert r > INTERIOR_BAND, (pt, r)


@complex_only
def test_printed_records(d25_m4):
    """Printed records (no band); asserts only that they were formed."""
    th_c, th_d = d25_m4["th"]["default"], d25_m4["th"]["d25"]
    th_s, th_s25 = d25_m4["th_s"], d25_m4["th_s25"]
    W = float(W_B2)
    for j in range(len(HIST_PTS)):
        dc_h = cem43(T_GRID, _T(th_c[:, j]))
        dd_h = cem43(T_GRID, _T(th_d[:, j]))
        de_h = cem43(T_GRID, _T(th_d[:, j] / DUTY25))
        dc_s = cem43(T_GRID, _T(th_s[:, j]))
        dd_s = cem43(T_GRID, _T(th_s25[:, j]))
        de_s = cem43(T_GRID, _T(th_s25[:, j] / DUTY25))
        M = float(np.max(np.abs(th_d[:, j] - th_s25[:, j])))
        cap = INTERIOR_BAND * float(np.max(th_s25[:, j]))
        _a, mx_h, t_h = running_average(T_GRID, th_d[:, j], window=W)
        _a, mx_s, t_s = running_average(T_GRID, th_s25[:, j], window=W)
        tag = "" if j < 3 else " [printed only]"
        _say(f"MAT21-14 DUTY25 DOSE {HIST_LABELS[j]}{tag}: FEM ratio equal-mean {de_h / dc_h:.6f} / "
             f"peak {dd_h / dc_h:.6f} (W1 closed form {de_s / dc_s:.6f} / {dd_s / dc_s:.6f}; P41 FEM "
             f"D = 0.5 equal-mean {P41_FEM_EQ[j]}, series {P41_EQ_RATIO[j]}); FEM D25 dose "
             f"{dd_h:.10e} s vs series {dd_s:.15e} s: m {abs(dd_h / dd_s - 1.0):.4e}, ln4*M "
             f"{LN4 * M:.4e} (M {M:.4e} K), headroom {cap / M:.1f}x; six-minute peak FEM "
             f"{mx_h:.9f} K at {t_h:.2f} s ({t_h / DT_BOX:.2f} DT_BOX) vs series {mx_s:.9f} K at "
             f"{t_s:.2f} s")
        assert np.isfinite(de_h / dc_h)
