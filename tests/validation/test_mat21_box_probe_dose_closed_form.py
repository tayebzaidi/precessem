"""`MAT-21` step 15, W1 — P39′'s dose construction on the gel box's series point
histories at P37's remaining F2182 probe points P3–P5 under P36's continuous
transverse quadratic source (§9 item 571, §10 chain P43; brief
`docs/planning/chunks/MAT-21.md` "§9 item 571", review ruling 2026-10-20 10:30).

Numpy, real build.  Sets no band and moves none; every band is imported.
Asserted at each ``NEW_PTS`` point on the 41-sample CN grid: (Q-s) ``cem43``
against step 3's per-segment ``quad`` at ``DOSE_BAR``, max T < 43 °C; (Win-s)
step 3's six-minute window construction at ``AVG_BAR``.  (R39′-s) the
``DOSE_PTS`` doses of ``_series_hist(HIST_PTS)`` print ``P39P_SERIES_DOSE``.
(Rc′-s) at P3–P5 is PREDICTED and printed only (review ruling 2026-10-20 10:30
(2)); it stays asserted at P0–P2 in item 565's module, which runs in the same
window.  ``NEW_PTS`` is a second fixed batch, never merged with ``HIST_PTS``.
"""

from __future__ import annotations

import numpy as np

from fem_em_solver.thermal.dose import cem43
from tests.validation.test_mat21_box_closed_form import A_BOX, C_BOX
from tests.validation.test_mat21_box_closed_form import T_INF as T_INF_0
from tests.validation.test_mat21_box_dose_closed_form import (
    AVG_BAR,
    DOSE_BAR,
    DOSE_PTS,
    HIST_PTS,
    T_GRID,
    T_REF,
    W_B2,
    _q_check,
    _series_hist,
    _true_dose,
    _win_check,
)
from tests.validation.test_mat21_box_duty_dose_closed_form import P39P_SERIES_DOSE

NEW_PTS = ((A_BOX, 0.0, 0.0), (-A_BOX, 0.0, 0.0), (0.0, 0.0, C_BOX))
NEW_LABELS = ("P3 +x wall", "P4 -x wall", "P5 top face")
SIX_LABELS = ("P1 +y wall", "P2 -y wall", "P0 centre") + NEW_LABELS

# Predictions written first (brief), printed only.
PRED_Q = 1e-14
PRED_WIN_III = 1e-15
PRED_RC = (1e6, 1e7)
PRED_TH = {0: (0.02, 0.065), 1: (0.02, 0.065), 2: (0.05, 0.2)}
PRED_DOSE_X_BASE = (1.00, 1.15)
PRED_FLOOR = 1e-4


def _p(msg):
    print(msg, flush=True)


_CACHE = {}


def _new():
    if "new" not in _CACHE:
        _CACHE["new"] = _series_hist(NEW_PTS)
    return _CACHE["new"]


def test_precondition_points():
    assert len(NEW_PTS) == 3 and len(set(NEW_PTS)) == 3
    assert not set(NEW_PTS) & set(HIST_PTS)
    assert T_GRID.size == 41 and T_GRID[0] == 0.0


def test_q_s_dose_against_quadrature_new_points():
    th = _new()
    for j, pt in enumerate(NEW_PTS):
        T = T_INF_0 + th[:, j]
        assert float(np.max(T)) < T_REF, (pt, float(np.max(T)))
        dose, ref, rel, swap, miss = _q_check(T)
        _p(
            f"MAT21-15s (Q-s) {NEW_LABELS[j]}: series CEM43 {dose:.15e} s vs per-segment quad "
            f"{ref:.15e} s: rel {rel:.3e} (bar {DOSE_BAR:g}, ASSERTED; predicted <= {PRED_Q:g}); "
            f"theta_s(t*) {th[-1, j]:.6f} K (predicted {PRED_TH[j]}); max T {np.max(T):.9f} C "
            f"< 43 ASSERTED"
        )
        _p(
            f"MAT21-15s control (Rc'-s) PREDICTED {NEW_LABELS[j]}: R swap {swap:.6e} s, miss "
            f"{miss:.4e} relative (printed only, ruling 10:30 (2)); factor {swap / ref:.4e} "
            f"(predicted range [{PRED_RC[0]:.0e}, {PRED_RC[1]:.0e}])"
        )
        assert rel <= DOSE_BAR, (pt, rel)


def test_win_s_six_minute_window_new_points():
    th = _new()
    for j, pt in enumerate(NEW_PTS):
        mx, t_mx, lo, hi, ref, rel = _win_check(th[:, j])
        _p(
            f"MAT21-15s (Win-s) {NEW_LABELS[j]}: W {float(W_B2):g} s: max six-minute mean "
            f"{mx:.15e} K at t_max {t_mx:.6f} s; (i) ASSERTED; (ii) [{lo:.12e}, {hi:.12e}] "
            f"ASSERTED; (iii) quad {ref:.15e}, rel {rel:.3e} (bar {AVG_BAR:g}, ASSERTED; "
            f"predicted <= {PRED_WIN_III:g})"
        )
        assert T_GRID[0] + W_B2 <= t_mx <= T_GRID[-1], t_mx
        assert lo <= mx <= hi, (lo, mx, hi)
        assert rel <= AVG_BAR, (pt, rel)


def test_r39p_s_hist_batch_reproduces_series_dose():
    th = _series_hist(HIST_PTS)
    got = tuple(f"{cem43(T_GRID, T_INF_0 + th[:, j]):.15e}" for j in range(len(DOSE_PTS)))
    _p(
        f"MAT21-15s (R39'-s) HIST_PTS batch DOSE_PTS doses {got} (P39P_SERIES_DOSE "
        f"{P39P_SERIES_DOSE}; strings ASSERTED)"
    )
    assert got == tuple(P39P_SERIES_DOSE), got


def test_printed_records_new_points():
    """Printed records; asserts only that the readings were formed."""
    th = _new()
    base = float(0.25 ** (T_REF - T_INF_0) * T_GRID[-1])
    _p(f"MAT21-15s baseline R^(43-20) t* {base:.15e} s")
    _p(
        f"MAT21-15s P3 / P4 series histories np.array_equal: "
        f"{bool(np.array_equal(th[:, 0], th[:, 1]))} (predicted True, printed)"
    )
    fine = _series_hist(NEW_PTS, n_terms=(120, 120, 60))
    _p(
        f"MAT21-15s comparand truncation max_n|hist(QUAD_TERMS) - hist(120,120,60)| per point "
        f"{[f'{float(np.max(np.abs(th[:, j] - fine[:, j]))):.3e}' for j in range(3)]} K "
        f"(printed only)"
    )
    doses = []
    for j, pt in enumerate(NEW_PTS):
        T = T_INF_0 + th[:, j]
        op = cem43(T_GRID, T)
        true, err = _true_dose(pt)
        mx, t_mx, *_r = _win_check(th[:, j])
        doses.append(op)
        _p(
            f"MAT21-15s SERIES DOSE {NEW_LABELS[j]}: 41-sample cem43 {op:.15e} s; true series "
            f"dose {true:.15e} s (quad err est {err:.1e}); sampling floor op/true - 1 "
            f"{op / true - 1.0:+.4e} (predicted |.| <= {PRED_FLOOR:g}); dose/baseline "
            f"{op / base:.6f} (predicted {PRED_DOSE_X_BASE}); six-minute peak {mx:.6f} K at "
            f"{t_mx:.3f} s; theta_s(t*) {th[-1, j]:.6f} K"
        )
        assert np.isfinite(op) and np.isfinite(true)
    six = tuple(P39P_SERIES_DOSE) + tuple(f"{d:.15e}" for d in doses)
    for lab, d in zip(SIX_LABELS, six):
        _p(f"MAT21-15s SIX-POINT SERIES DOSE TABLE {lab}: {d} s")
