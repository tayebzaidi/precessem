"""`MAT-21` step 16, W1 — `MAT-20`'s dose and six-minute window on the gel box's
series point histories at P37's probe points P3–P5 under P41's D = 0.5 and
P42's D = 0.25 schedules at P38's period (P36's source) (§9 item 577, §10
chain P44; brief `docs/planning/chunks/MAT-21.md` "§9 item 577", review ruling
2026-10-20 18:00).

Numpy, real build.  Sets no band and moves none; every band is imported.
``HIST_PTS`` and ``NEW_PTS`` are two fixed batches, never merged.  Asserted:
(C100-sP) the 100 % schedule bitwise at ``NEW_PTS`` and its dose strings;
(R43-s) P43's series dose strings; (R42-s) P42's D = 0.25 equal-mean ratio
strings; (Q-sP) / (Win-sP) at P3–P5 on both duty histories at peak and equal
mean power; max T < 43; (Rc′-sP) at P0–P2 on both duty histories.  (Rc′-sP)
at P3–P5 is PREDICTED and printed only (ruling 18:00 (3)).
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
from tests.validation.test_mat21_box_duty25_dose_closed_form import DUTY25, P41_EQ_RATIO, T_SW25
from tests.validation.test_mat21_box_duty_closed_form import DT_BOX, DUTY, T_SW
from tests.validation.test_mat21_box_duty_dose_closed_form import _duty_hist, _trap
from tests.validation.test_mat21_box_probe_dose_closed_form import NEW_LABELS, NEW_PTS, SIX_LABELS

SCHEDS = (("d50", DUTY, T_SW), ("d25", DUTY25, T_SW25))
# P42's W1 D = 0.25 equal-mean ratios at :.6f, HIST_PTS order, copied from
# docs/testing/logs/20261020T185346Z_MAT-21-step14-W1-duty25-dose-closed-form.log:110-113.
P42_EQ_RATIO = ("1.276608", "1.276608", "1.006916", "1.317022")
# P43's W1 series doses at :.15e, NEW_PTS order, copied from
# docs/testing/logs/20261020T190033Z_MAT-21-step15-W1-probe-dose-closed-form.log:86-88.
P43_SERIES_DOSE = ("4.109185736623872e-11", "4.109185736623872e-11", "4.318074242697235e-11")
# Predictions written first (brief), printed only.
PRED_EQ = {"d50": ((1.001, 1.008), (1.004, 1.02)), "d25": ((1.002, 1.012), (1.006, 1.03))}
PRED_PEAK = {"d50": ((0.982, 0.995), (0.955, 0.98)), "d25": ((0.975, 0.99), (0.935, 0.965))}
PRED_CONVEX = (1.000, 1.012)
PRED_RC = (1e6, 1e7)
REC_P02 = {"d50": ("1.186312", "1.004475"), "d25": ("1.276608", "1.006916")}


def _p(msg):
    print(msg, flush=True)


_CACHE = {}


def _hists():
    if not _CACHE:
        _CACHE["s_h"] = _series_hist(HIST_PTS)
        _CACHE["s_n"] = _series_hist(NEW_PTS)
        for key, _d, sw in SCHEDS:
            _CACHE[f"{key}_h"] = _duty_hist(HIST_PTS, switch_times=sw)
            _CACHE[f"{key}_n"] = _duty_hist(NEW_PTS, switch_times=sw)
        _CACHE["s1_n"] = _duty_hist(NEW_PTS, switch_times=(0.0,))
    return _CACHE


def test_c100_sp_full_schedule_bitwise_new_points():
    c = _hists()
    eq = bool(np.array_equal(c["s1_n"], c["s_n"]))
    _p(
        f"MAT21-16s (C100-sP) NEW_PTS switch_times=(0.0,) array_equal _series_hist over "
        f"{c['s_n'].shape}: {eq}"
    )
    assert eq
    got = tuple(f"{cem43(T_GRID, T_INF_0 + c['s1_n'][:, j]):.15e}" for j in range(len(NEW_PTS)))
    _p(f"MAT21-16s (C100-sP) doses(100 %) {got} (P43 {P43_SERIES_DOSE})")
    assert got == P43_SERIES_DOSE


def test_r43_s_series_dose_strings():
    c = _hists()
    got = tuple(f"{cem43(T_GRID, T_INF_0 + c['s_n'][:, j]):.15e}" for j in range(len(NEW_PTS)))
    _p(f"MAT21-16s (R43-s) NEW_PTS series doses {got} (P43 {P43_SERIES_DOSE})")
    assert got == P43_SERIES_DOSE


def test_r42_s_d25_equal_mean_strings():
    c = _hists()
    th_s, th_d = c["s_h"], c["d25_h"]
    got = []
    for j in range(len(HIST_PTS)):
        d_c = cem43(T_GRID, T_INF_0 + th_s[:, j])
        d_eq = cem43(T_GRID, T_INF_0 + th_d[:, j] / DUTY25)
        got.append(f"{d_eq / d_c:.6f}")
    _p(f"MAT21-16s (R42-s) D = 0.25 equal-mean ratios {tuple(got)} (P42 {P42_EQ_RATIO})")
    assert tuple(got) == P42_EQ_RATIO


def test_q_sp_win_sp_new_points_both_schedules():
    c = _hists()
    for key, D, _sw in SCHEDS:
        th_d = c[f"{key}_n"]
        for scale, lab in ((1.0, "peak"), (1.0 / D, "equal-mean")):
            for j, pt in enumerate(NEW_PTS):
                th = th_d[:, j] * scale
                T = T_INF_0 + th
                assert float(np.max(T)) < T_REF, (key, lab, pt, float(np.max(T)))
                dose, ref, rel, swap, miss = _q_check(T)
                mx, t_mx, lo, hi, wref, wrel = _win_check(th)
                _p(
                    f"MAT21-16s (Q-sP) {key} {lab} {NEW_LABELS[j]}: CEM43 {dose:.15e} s vs quad "
                    f"{ref:.15e} s rel {rel:.3e} (bar {DOSE_BAR:g}, ASSERTED; predicted <= 1e-14); "
                    f"max T {np.max(T):.9f} C < 43 ASSERTED"
                )
                _p(
                    f"MAT21-16s (Win-sP) {key} {lab} {NEW_LABELS[j]}: max six-minute mean "
                    f"{mx:.15e} K at t_max {t_mx:.4f} s ({t_mx / DT_BOX:.2f} DT_BOX); (ii) "
                    f"[{lo:.12e}, {hi:.12e}]; (iii) quad {wref:.15e} rel {wrel:.3e} (bar "
                    f"{AVG_BAR:g}, ASSERTED; predicted <= 1e-15)"
                )
                _p(
                    f"MAT21-16s control (Rc'-sP) PREDICTED {key} {lab} {NEW_LABELS[j]}: R swap "
                    f"miss {miss:.4e}; factor {swap / ref:.4e} (predicted range "
                    f"[{PRED_RC[0]:.0e}, {PRED_RC[1]:.0e}]; printed only, ruling 18:00 (3))"
                )
                assert rel <= DOSE_BAR, (key, lab, pt, rel)
                assert T_GRID[0] + W_B2 <= t_mx <= T_GRID[-1], (key, lab, t_mx)
                assert lo <= mx <= hi, (key, lab, lo, mx, hi)
                assert wrel <= AVG_BAR, (key, lab, pt, wrel)


def test_rc_sp_hist_points_both_schedules():
    c = _hists()
    for key, D, _sw in SCHEDS:
        th_d = c[f"{key}_h"]
        for scale, lab in ((1.0, "peak"), (1.0 / D, "equal-mean")):
            for j, pt in enumerate(DOSE_PTS):
                T = T_INF_0 + th_d[:, j] * scale
                assert float(np.max(T)) < T_REF, (key, lab, pt)
                _dose, ref, _rel, swap, miss = _q_check(T)
                _p(
                    f"MAT21-16s (Rc'-sP) {key} {lab} {HIST_LABELS[j]}: R swap miss {miss:.4e} "
                    f"(asserted > {DOSE_BAR:g}); factor {swap / ref:.4e}"
                )
                assert miss > DOSE_BAR, (key, lab, pt, miss)


def test_printed_tables():
    """The record; asserts only that it was formed."""
    c = _hists()
    W = float(W_B2)
    th_s = c["s_n"]
    eqs = {}
    for key, D, _sw in SCHEDS:
        th_d = c[f"{key}_n"]
        for j in range(len(NEW_PTS)):
            kind = 1 if j == 2 else 0
            d_c = cem43(T_GRID, T_INF_0 + th_s[:, j])
            d_pk = cem43(T_GRID, T_INF_0 + th_d[:, j])
            d_eq = cem43(T_GRID, T_INF_0 + th_d[:, j] / D)
            eqs[(key, j)] = d_eq
            r_t = th_d[-1, j] / th_s[-1, j]
            r_mean = _trap(th_d[:, j]) / _trap(th_s[:, j])
            _a, mx_c, t_c = running_average(T_GRID, th_s[:, j], window=W)
            _a, mx_d, t_d = running_average(T_GRID, th_d[:, j], window=W)
            _p(
                f"MAT21-16s DUTY DOSE {key} {NEW_LABELS[j]}: D_cont {d_c:.15e} s; D(peak) "
                f"{d_pk:.15e} s; D(equal mean) {d_eq:.15e} s; ratio equal-mean {d_eq / d_c:.6f} "
                f"(predicted {PRED_EQ[key][kind]}); ratio peak {d_pk / d_c:.6f} (predicted "
                f"{PRED_PEAK[key][kind]}); P0-P2 recorded eq wall/centre {REC_P02[key]}; theta "
                f"ratio at t* {r_t:.5f} (x1/D {r_t / D:.5f}); time-mean theta ratio {r_mean:.5f} "
                f"(x1/D {r_mean / D:.5f}); six-minute peak cont {mx_c:.6f} K at {t_c:.2f} s, duty "
                f"{mx_d:.6f} K at {t_d:.2f} s ({t_d / DT_BOX:.2f} DT_BOX)"
            )
            assert np.isfinite(d_eq / d_c)
    for j in range(len(NEW_PTS)):
        _p(
            f"MAT21-16s convexity D25,eq/D50,eq {NEW_LABELS[j]}: "
            f"{eqs[('d25', j)] / eqs[('d50', j)]:.6f} (predicted {PRED_CONVEX})"
        )
    for key, _d, _sw in SCHEDS:
        _p(
            f"MAT21-16s P3 / P4 {key} series histories np.array_equal: "
            f"{bool(np.array_equal(c[f'{key}_n'][:, 0], c[f'{key}_n'][:, 1]))} (predicted True)"
        )
    for norm in ("peak", "equal-mean"):
        for i, lab in enumerate(SIX_LABELS):
            arr_key, jj = ("h", i) if i < 3 else ("n", i - 3)
            row = [f"cont {cem43(T_GRID, T_INF_0 + c[f's_{arr_key}'][:, jj]):.15e}"]
            for key, D, _sw in SCHEDS:
                sc = 1.0 if norm == "peak" else 1.0 / D
                row.append(
                    f"{key} {cem43(T_GRID, T_INF_0 + c[f'{key}_{arr_key}'][:, jj] * sc):.15e}"
                )
            _p(
                f"MAT21-16s SIX-POINT x THREE-SCHEDULE SERIES DOSE TABLE ({norm}) {lab}: "
                + "; ".join(row)
                + " s"
            )
    _p(f"MAT21-16s (R41-s) reference P41_EQ_RATIO {P41_EQ_RATIO} (asserted by P42's module)")
