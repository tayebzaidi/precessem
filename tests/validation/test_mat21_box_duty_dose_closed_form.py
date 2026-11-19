"""`MAT-21` step 13, W1 — `MAT-20`'s dose and six-minute window on the gel box's
duty-cycled series point histories (P38's schedule, P36's transverse quadratic
source) against the continuous source (§9 item 567, §10 chain P41; brief
`docs/planning/chunks/MAT-21.md` "§9 item 567", review ruling 2026-10-20 03:00).

Numpy, real build.  Sets no band and moves none; every band is imported.
Asserted: (C100-s) the 100 % schedule (``switch_times=(0.0,)``) history is
``np.array_equal`` to item 565's ``_series_hist`` over all 41 rows, its dose
``==`` item 565's dose and its ``:.15e`` string equal to ``P39P_SERIES_DOSE``;
(Q-sd) / (Win-sd) on the duty history at peak power and at equal mean power
(θ_duty / DUTY — a scale, not a solve) at ``DOSE_BAR`` / ``AVG_BAR``, max
T < 43 °C; (Rc′-sd) the R swap misses > ``DOSE_BAR``.  Printed: the duty ratio
of dose beside P38's temperature ratios.
"""

from __future__ import annotations

import numpy as np

from fem_em_solver.thermal.dose import cem43, running_average
from tests.validation.test_mat21_box_closed_form import T_INF as T_INF_0
from tests.validation.test_mat21_box_dose_closed_form import (
    AVG_BAR,
    DOSE_BAR,
    DOSE_PTS,
    HIST_LABELS,
    HIST_PTS,
    T_GRID,
    T_REF,
    W_B2,
    _q_check,
    _series_hist,
    _win_check,
)
from tests.validation.test_mat21_box_duty_closed_form import DUTY, T_STAR, T_SW
from tests.validation.test_mat21_box_quadratic_closed_form import QUAD_TERMS
from tests.validation.test_mat21_box_quadratic_duty_closed_form import _quad_duty

# Item 565's W1 series doses at :.15e, DOSE_PTS order, copied from
# docs/testing/logs/20261020T140156Z_MAT-21-step12-W1-dose-closed-form.log:108-110.
P39P_SERIES_DOSE = ("4.699458195189890e-10", "4.699458195189890e-10", "4.136468998916232e-11")
# P38's W1 closed-form temperature ratios theta_duty(t*)/theta_cont(t*) (brief, printed only).
P38_RATIO = {"(0,+B,0)": 0.48560, "(0,-B,0)": 0.48560, "centre": 0.54466, "(0,0.9B,0)": 0.49468}
PRED_Q = 1e-14
PRED_EQ = {"centre": (1.00, 1.01), "wall": (0.9, 1.3)}
PRED_PEAK = (0.3, 1.0)


def _p(msg):
    print(msg, flush=True)


def _duty_hist(pts, switch_times=T_SW):
    """(41, n_pts) duty series histories; the arrays ``_series_hist`` passes, row 0 = 0.0."""
    P = np.asarray(pts, dtype=float)
    out = np.zeros((T_GRID.size, P.shape[0]))
    for n in range(1, T_GRID.size):
        out[n] = _quad_duty(
            P[:, 0], P[:, 1], P[:, 2], T_GRID[n], switch_times=switch_times, n_terms=QUAD_TERMS
        )
    return out


def _trap(y):
    return float(np.sum(0.5 * (y[1:] + y[:-1]) * np.diff(T_GRID)))


_CACHE = {}


def _hists():
    if not _CACHE:
        _CACHE["s"] = _series_hist(HIST_PTS)
        _CACHE["s1"] = _duty_hist(HIST_PTS, switch_times=(0.0,))
        _CACHE["d"] = _duty_hist(HIST_PTS)
    return _CACHE["s"], _CACHE["s1"], _CACHE["d"]


def test_c100_s_full_schedule_bitwise_p39p():
    th_s, th_s1, _d = _hists()
    eq = bool(np.array_equal(th_s1, th_s))
    _p(
        f"MAT21-13s (C100-s) theta_s1 (switch_times=(0.0,)) array_equal _series_hist over "
        f"{th_s.shape}: {eq} (max |d| {np.max(np.abs(th_s1 - th_s)):.3e} K)"
    )
    assert eq
    for j, pt in enumerate(DOSE_PTS):
        d1 = cem43(T_GRID, T_INF_0 + th_s1[:, j])
        d0 = cem43(T_GRID, T_INF_0 + th_s[:, j])
        _p(
            f"MAT21-13s (C100-s) {HIST_LABELS[j]}: dose(100 %) {d1:.15e} s, dose(P39') {d0:.15e} s, "
            f"== {d1 == d0}; P39' string {P39P_SERIES_DOSE[j]}"
        )
        assert d1 == d0, (pt, d1, d0)
        assert f"{d1:.15e}" == P39P_SERIES_DOSE[j], (pt, f"{d1:.15e}")


def test_q_sd_win_sd_and_rc_sd():
    _s, _s1, th_d = _hists()
    for scale, lab in ((1.0, "peak"), (1.0 / DUTY, "equal-mean")):
        for j, pt in enumerate(DOSE_PTS):
            th = th_d[:, j] * scale
            T = T_INF_0 + th
            assert float(np.max(T)) < T_REF, (lab, pt, float(np.max(T)))
            dose, ref, rel, swap, miss = _q_check(T)
            mx, t_mx, lo, hi, wref, wrel = _win_check(th)
            _p(
                f"MAT21-13s (Q-sd) {lab} {HIST_LABELS[j]}: duty CEM43 {dose:.15e} s vs quad "
                f"{ref:.15e} s rel {rel:.3e} (bar {DOSE_BAR:g}, ASSERTED; predicted <= {PRED_Q:g}); "
                f"max T {np.max(T):.9f} C < 43 ASSERTED"
            )
            _p(
                f"MAT21-13s (Win-sd) {lab} {HIST_LABELS[j]}: max six-minute mean {mx:.15e} K at "
                f"t_max {t_mx:.4f} s (predicted near 2539.35 s at the walls, later at the centre); "
                f"(ii) [{lo:.12e}, {hi:.12e}]; (iii) quad {wref:.15e} rel {wrel:.3e} "
                f"(bar {AVG_BAR:g}, ASSERTED)"
            )
            _p(
                f"MAT21-13s (Rc'-sd) {lab} {HIST_LABELS[j]}: R swap miss {miss:.4e} (asserted > "
                f"{DOSE_BAR:g}); factor {swap / ref:.4e}"
            )
            assert rel <= DOSE_BAR, (lab, pt, rel)
            assert T_GRID[0] + W_B2 <= t_mx <= T_GRID[-1], (lab, t_mx)
            assert lo <= mx <= hi, (lab, lo, mx, hi)
            assert wrel <= AVG_BAR, (lab, pt, wrel)
            assert miss > DOSE_BAR, (lab, pt, miss)


def test_printed_duty_ratio_of_dose():
    """The record the weekly asked for; asserts only that it was formed."""
    th_s, _s1, th_d = _hists()
    W = float(W_B2)
    for j in range(len(HIST_PTS)):
        d_c = cem43(T_GRID, T_INF_0 + th_s[:, j])
        d_pk = cem43(T_GRID, T_INF_0 + th_d[:, j])
        d_eq = cem43(T_GRID, T_INF_0 + th_d[:, j] / DUTY)
        r_t = th_d[-1, j] / th_s[-1, j]
        r_mean = _trap(th_d[:, j]) / _trap(th_s[:, j])
        _a, mx_c, t_c = running_average(T_GRID, th_s[:, j], window=W)
        _a, mx_d, t_d = running_average(T_GRID, th_d[:, j], window=W)
        kind = "centre" if j == 2 else "wall"
        _p(
            f"MAT21-13s DUTY DOSE {HIST_LABELS[j]}: D_cont {d_c:.15e} s; D_duty(peak) {d_pk:.15e} s; "
            f"D_duty(equal mean) {d_eq:.15e} s; ratio equal-mean {d_eq / d_c:.6f} (predicted "
            f"{PRED_EQ[kind]}); ratio peak {d_pk / d_c:.6f} (predicted {PRED_PEAK}); theta ratio at "
            f"t* {r_t:.5f} (P38 {P38_RATIO[HIST_LABELS[j]]:.5f}; x2 {2 * r_t:.5f}); time-mean theta "
            f"ratio {r_mean:.5f} (x2 {r_mean / DUTY:.5f}); six-minute peak cont {mx_c:.6f} K at "
            f"{t_c:.2f} s, duty {mx_d:.6f} K at {t_d:.2f} s (equal mean {mx_d / DUTY:.6f} K); "
            f"t* {T_STAR:.2f} s"
        )
        assert np.isfinite(d_eq / d_c)
