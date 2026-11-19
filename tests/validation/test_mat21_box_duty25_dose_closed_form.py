"""`MAT-21` step 14, W1 — `MAT-20`'s dose and six-minute window on the gel box's
series point histories at a second duty fraction, D = 0.25 at P38's period
(§9 item 570, §10 chain P42; brief `docs/planning/chunks/MAT-21.md`
"§9 item 570", review ruling 2026-10-20 10:30).

Numpy, real build.  Sets no band and moves none; every band is imported.
Asserted: (G25) schedule consistency (midpoint string, closed-form state,
census ≤ 1e-12, switches on the grid); (C100-s) the 100 % schedule bitwise
against item 565; (R41-s) P41's D = 0.5 equal-mean ratio strings; (Q-s25) /
(Win-s25) at peak and equal mean power (θ_25 / DUTY25, a scale); max T < 43;
(Rc′-s25) the R swap misses; (G-ctl) a planted three-step on-phase fails (G25).
Printed: the D = 0.25 vs D = 0.5 dose-ratio table.
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
from tests.validation.test_mat21_box_duty_closed_form import DT_BOX, DUTY, HALF_T, T_STAR
from tests.validation.test_mat21_box_duty_dose_closed_form import (
    P39P_SERIES_DOSE,
    _duty_hist,
    _trap,
)

DUTY25 = 0.25
# Built as k * DT_BOX (never j * P + P / 4): each entry is the same float product as T_GRID[k].
T_SW25 = tuple(k * DT_BOX for k in (0, 2, 8, 10, 16, 18, 24, 26, 32, 34))
SCHEDULE_25 = "11100000011000000110000001100000011000000"
# The solver's 41 (half-)step midpoints: two start-up half-steps, then k = 2 ... 40.
MID_T = (0.25 * DT_BOX, 0.75 * DT_BOX) + tuple((k - 0.5) * DT_BOX for k in range(2, 41))
# P41's W1 equal-mean dose ratios at :.6f, HIST_PTS order, copied from
# docs/testing/logs/20261020T140919Z_MAT-21-step13-W1-duty-dose-closed-form.log:93-96.
P41_EQ_RATIO = ("1.186312", "1.186312", "1.004475", "1.213031")
PRED_EQ = {"wall": (1.2, 2.0), "centre": (1.000, 1.03)}
PRED_PEAK = {"wall": (0.08, 0.25), "centre": (0.96, 0.99)}
PRED_CONVEX = (1.0, 1.7)
PRED_TMAX = (34, 40)


def _p(msg):
    print(msg, flush=True)


def _chi25(t):
    """1 on [jP, jP + P/4), 0 otherwise (P = 8 DT_BOX); t is always a (half-)step midpoint."""
    return 1.0 if int(np.floor(t / DT_BOX)) % 8 < 2 else 0.0


def _chi25_planted(t):
    """(G-ctl): the on-phase one step too long."""
    return 1.0 if int(np.floor(t / DT_BOX)) % 8 < 3 else 0.0


def _sched_string(chi):
    return "".join(str(int(chi(t))) for t in MID_T)


_CACHE = {}


def _hists():
    if not _CACHE:
        _CACHE["s"] = _series_hist(HIST_PTS)
        _CACHE["s1"] = _duty_hist(HIST_PTS, switch_times=(0.0,))
        _CACHE["d50"] = _duty_hist(HIST_PTS)
        _CACHE["d25"] = _duty_hist(HIST_PTS, switch_times=T_SW25)
    return _CACHE["s"], _CACHE["s1"], _CACHE["d50"], _CACHE["d25"]


def test_g25_schedule_consistency():
    s = _sched_string(_chi25)
    _p(f"MAT21-14s (G25) schedule {s} (must equal {SCHEDULE_25}); len {len(s)}")
    assert len(MID_T) == 41
    assert s == SCHEDULE_25
    for t in MID_T:
        n_before = sum(1 for tm in T_SW25 if tm < t)
        assert _chi25(t) == float(n_before % 2), (t, n_before)
    on = sum(T_SW25[2 * j + 1] - T_SW25[2 * j] for j in range(len(T_SW25) // 2))
    rel = abs(on - DUTY25 * T_STAR) / (DUTY25 * T_STAR)
    _p(
        f"MAT21-14s (G25) on-time {on:.12f} s vs DUTY25*T_STAR {DUTY25 * T_STAR:.12f} s rel "
        f"{rel:.3e} (bar 1e-12)"
    )
    assert rel <= 1e-12
    for k, tm in zip((0, 2, 8, 10, 16, 18, 24, 26, 32, 34), T_SW25):
        assert tm == T_GRID[k], (k, tm, T_GRID[k])
    _p(f"MAT21-14s (G25) every T_SW25 entry == its T_GRID sample: True; P = {2 * HALF_T:.4f} s")


def test_g_ctl_planted_schedule_fails():
    s = _sched_string(_chi25_planted)
    flagged = s != SCHEDULE_25
    _p(
        f"MAT21-14s (G-ctl) planted (% 8 < 3) schedule {s}: comparator flags it {flagged} "
        f"(ASSERTED by construction)"
    )
    assert flagged


def test_c100_s_full_schedule_bitwise_p39p():
    th_s, th_s1, _d50, _d25 = _hists()
    eq = bool(np.array_equal(th_s1, th_s))
    _p(f"MAT21-14s (C100-s) switch_times=(0.0,) array_equal _series_hist over {th_s.shape}: {eq}")
    assert eq
    for j, pt in enumerate(DOSE_PTS):
        d1 = cem43(T_GRID, T_INF_0 + th_s1[:, j])
        _p(
            f"MAT21-14s (C100-s) {HIST_LABELS[j]}: dose(100 %) {d1:.15e} s; P39' {P39P_SERIES_DOSE[j]}"
        )
        assert f"{d1:.15e}" == P39P_SERIES_DOSE[j], (pt, f"{d1:.15e}")


def test_r41_s_d50_baseline_strings():
    th_s, _s1, th_d50, _d25 = _hists()
    got = []
    for j in range(len(HIST_PTS)):
        d_c = cem43(T_GRID, T_INF_0 + th_s[:, j])
        d_eq = cem43(T_GRID, T_INF_0 + th_d50[:, j] / DUTY)
        got.append(f"{d_eq / d_c:.6f}")
    _p(f"MAT21-14s (R41-s) D = 0.5 equal-mean ratios {tuple(got)} (P41 {P41_EQ_RATIO})")
    assert tuple(got) == P41_EQ_RATIO


def test_q_s25_win_s25_and_rc_s25():
    _s, _s1, _d50, th_d = _hists()
    for scale, lab in ((1.0, "peak"), (1.0 / DUTY25, "equal-mean")):
        for j, pt in enumerate(DOSE_PTS):
            th = th_d[:, j] * scale
            T = T_INF_0 + th
            assert float(np.max(T)) < T_REF, (lab, pt, float(np.max(T)))
            dose, ref, rel, swap, miss = _q_check(T)
            mx, t_mx, lo, hi, wref, wrel = _win_check(th)
            _p(
                f"MAT21-14s (Q-s25) {lab} {HIST_LABELS[j]}: CEM43 {dose:.15e} s vs quad {ref:.15e} s "
                f"rel {rel:.3e} (bar {DOSE_BAR:g}, ASSERTED; predicted <= 1e-14); max T "
                f"{np.max(T):.9f} C < 43 ASSERTED"
            )
            _p(
                f"MAT21-14s (Win-s25) {lab} {HIST_LABELS[j]}: max six-minute mean {mx:.15e} K at "
                f"t_max {t_mx:.4f} s ({t_mx / DT_BOX:.2f} DT_BOX); (ii) [{lo:.12e}, {hi:.12e}]; "
                f"(iii) quad {wref:.15e} rel {wrel:.3e} (bar {AVG_BAR:g}, ASSERTED; predicted <= 1e-15)"
            )
            _p(
                f"MAT21-14s (Rc'-s25) {lab} {HIST_LABELS[j]}: R swap miss {miss:.4e} (asserted > "
                f"{DOSE_BAR:g}); factor {swap / ref:.4e}"
            )
            assert rel <= DOSE_BAR, (lab, pt, rel)
            assert T_GRID[0] + W_B2 <= t_mx <= T_GRID[-1], (lab, t_mx)
            assert lo <= mx <= hi, (lab, lo, mx, hi)
            assert wrel <= AVG_BAR, (lab, pt, wrel)
            assert miss > DOSE_BAR, (lab, pt, miss)


def test_printed_duty25_ratio_table():
    """The record (D = 0.25 vs D = 0.5); asserts only that it was formed."""
    th_s, _s1, th_d50, th_d = _hists()
    W = float(W_B2)
    for j in range(len(HIST_PTS)):
        kind = "centre" if j == 2 else "wall"
        d_c = cem43(T_GRID, T_INF_0 + th_s[:, j])
        d_pk = cem43(T_GRID, T_INF_0 + th_d[:, j])
        d_eq = cem43(T_GRID, T_INF_0 + th_d[:, j] / DUTY25)
        d50_eq = cem43(T_GRID, T_INF_0 + th_d50[:, j] / DUTY)
        d50_pk = cem43(T_GRID, T_INF_0 + th_d50[:, j])
        r_t = th_d[-1, j] / th_s[-1, j]
        r_mean = _trap(th_d[:, j]) / _trap(th_s[:, j])
        _a, mx_c, t_c = running_average(T_GRID, th_s[:, j], window=W)
        _a, mx_d, t_d = running_average(T_GRID, th_d[:, j], window=W)
        _p(
            f"MAT21-14s DUTY25 DOSE {HIST_LABELS[j]}: D_cont {d_c:.15e} s; D25(peak) {d_pk:.15e} s; "
            f"D25(equal mean) {d_eq:.15e} s; ratio equal-mean {d_eq / d_c:.6f} (predicted "
            f"{PRED_EQ[kind]}; D = 0.5 P41 {P41_EQ_RATIO[j]}); ratio peak {d_pk / d_c:.6f} (predicted "
            f"{PRED_PEAK[kind]}; D = 0.5 {d50_pk / d_c:.6f}); convexity D25,eq/D50,eq "
            f"{d_eq / d50_eq:.6f} (predicted {PRED_CONVEX} at the walls); theta ratio at t* {r_t:.5f} "
            f"(x4 {r_t / DUTY25:.5f}); time-mean theta ratio {r_mean:.5f} (x4 {r_mean / DUTY25:.5f}); "
            f"six-minute peak cont {mx_c:.6f} K at {t_c:.2f} s, D25 {mx_d:.6f} K at {t_d:.2f} s "
            f"({t_d / DT_BOX:.2f} DT_BOX; predicted {PRED_TMAX} at the walls, no weight; equal mean "
            f"{mx_d / DUTY25:.6f} K)"
        )
        assert np.isfinite(d_eq / d_c)
