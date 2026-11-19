"""`MAT-21` step 12, W1 — `MAT-20`'s dose and six-minute window on the gel box's
series point histories under P36's transverse quadratic source (§9 item 565,
§10 chain P39′; brief `docs/planning/chunks/MAT-21.md` "§9 item 565", review
ruling 2026-10-20 03:00).

Numpy, real build.  Sets no band and moves none; every band is imported.
Asserted at each ``DOSE_PTS`` point on the 41-sample CN grid t_n = n·t*/40:
(Q-s) ``cem43`` against step 3's per-segment ``quad`` at ``DOSE_BAR``, max
T < 43 °C; (Win-s) step 3's six-minute window construction at ``AVG_BAR``;
(A2) `MAT-20` step 1's sphere control reproducing `mat:17`'s record; (Rc′-s)
the R swap misses (Q-s)'s reference by > ``DOSE_BAR``.  Printed only: the true
series dose (adaptive ``quad`` on the closed form) and the 41-sample sampling
floor against it, every series dose at ``:.15e``, the (0, 0.9B, 0) row.

Constants ``DOSE_PTS``, ``HIST_PTS`` (the one evaluation batch, in this order),
``T_GRID``, ``MAT17_REC_A2`` and the helpers ``_series_hist``, ``_q_check``,
``_win_check`` are imported by W2 and by item 567.
"""

from __future__ import annotations

import math

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

from fem_em_solver.thermal.dose import cem43, running_average
from fem_em_solver.utils.analytical import robin_box_quadratic_transient_temperature
from tests.validation.test_mat20_cem43_closed_form import BAR as DOSE_BAR
from tests.validation.test_mat20_cem43_closed_form import T_END, _rel, _T_centre
from tests.validation.test_mat20_running_average import BAR as AVG_BAR
from tests.validation.test_mat20_running_average import W_B2
from tests.validation.test_mat20_tuned_heating_dose import R_BELOW, T_REF, _quad_dose
from tests.validation.test_mat21_box_closed_form import B_BOX, HALF
from tests.validation.test_mat21_box_closed_form import RHO_C as RHO_C_0
from tests.validation.test_mat21_box_closed_form import T_INF as T_INF_0
from tests.validation.test_mat21_box_closed_form import H, K, Q
from tests.validation.test_mat21_box_duty_closed_form import DT_BOX, N_STEPS, T_STAR
from tests.validation.test_mat21_box_quadratic_closed_form import QUAD_TERMS, WEIGHTS

DOSE_PTS = ((0.0, B_BOX, 0.0), (0.0, -B_BOX, 0.0), (0.0, 0.0, 0.0))
HIST_PTS = DOSE_PTS + ((0.0, 0.9 * B_BOX, 0.0),)
HIST_LABELS = ("(0,+B,0)", "(0,-B,0)", "centre", "(0,0.9B,0)")
T_GRID = DT_BOX * np.arange(N_STEPS + 1)
# `mat:17`'s record: examples/materials/17_thermal_dose_duty_cycle.py:106-109 and
# the gate log docs/testing/logs/20261008T171434Z_MAT-20-step1.log.gz (strings).
MAT17_REC_A2 = {
    "t43": "14087.0956163533",
    "quad": "1.708516120835631e+05",
    "op262144": "1.708516120820288e+05",
    "rel": "8.980e-12",
}
LN_R = math.log(R_BELOW)

# Predictions written first (brief), printed only.
PRED_Q = 1e-14
PRED_WIN_III = 1e-15
PRED_FLOOR = {"wall": -7.24e-5, "centre": 8.15e-6}
PRED_RC = (1e5, 1e7)


def _p(msg):
    print(msg, flush=True)


def _series_hist(pts, weights=WEIGHTS, n_terms=QUAD_TERMS):
    """(41, n_pts) series θ histories on ``T_GRID``; row 0 is θ(0) = 0 exactly."""
    P = np.asarray(pts, dtype=float)
    out = np.zeros((T_GRID.size, P.shape[0]))
    for n in range(1, T_GRID.size):
        out[n] = robin_box_quadratic_transient_temperature(
            P[:, 0],
            P[:, 1],
            P[:, 2],
            T_GRID[n],
            half_widths=HALF,
            k=K,
            q=Q,
            h=H,
            rho_c=RHO_C_0,
            weights=weights,
            t_inf=0.0,
            n_terms=n_terms,
        )
    return out


def _q_check(T):
    """Step 3's (Q) construction: (dose, quad reference, rel, swap, swap miss)."""
    dose = cem43(T_GRID, T)
    ref = _quad_dose(T_GRID, T, LN_R)
    swap = cem43(T_GRID, T, r_above=R_BELOW, r_below=0.5)
    return dose, ref, _rel(dose, ref), swap, _rel(swap, ref)


def _win_check(th, window=float(W_B2)):
    """Step 3's (Win) construction on one θ history: (mx, t_mx, lo, hi, ref, rel)."""
    t = T_GRID
    _means, mx, t_mx = running_average(t, th, window=window)
    s = t_mx - window
    inside = t[(t > s) & (t <= t_mx)]
    pts = np.concatenate([[np.interp(s, t, th)], np.interp(inside, t, th)])
    lo, hi = float(np.min(pts)), float(np.max(pts))
    brk = [x for x in t if s < x < t_mx]
    ref, _err = quad(
        lambda x: np.interp(x, t, th),
        s,
        t_mx,
        points=brk or None,
        epsabs=0.0,
        epsrel=1e-13,
        limit=200,
    )
    ref /= window
    return float(mx), float(t_mx), lo, hi, ref, _rel(mx, ref)


def _true_dose(pt):
    """Adaptive quad of R^(43 − T_s(t)) on the closed form over [0, t*] (printed)."""
    x, y, z = pt

    def f(t):
        th = float(
            robin_box_quadratic_transient_temperature(
                x,
                y,
                z,
                t,
                half_widths=HALF,
                k=K,
                q=Q,
                h=H,
                rho_c=RHO_C_0,
                weights=WEIGHTS,
                t_inf=0.0,
                n_terms=QUAD_TERMS,
            )
        )
        return math.exp(LN_R * (T_REF - (T_INF_0 + th)))

    val, err = quad(f, 0.0, T_STAR, epsabs=0.0, epsrel=1e-13, limit=400)
    return val, err


_CACHE = {}


def _hist():
    if "th" not in _CACHE:
        _CACHE["th"] = _series_hist(HIST_PTS)
    return _CACHE["th"]


def test_precondition_grid():
    assert T_GRID.size == 41 and T_GRID[0] == 0.0
    assert abs(T_GRID[-1] / T_STAR - 1.0) <= 1e-15
    assert HIST_PTS[:3] == DOSE_PTS


def test_q_s_dose_against_quadrature():
    th = _hist()
    for j, pt in enumerate(DOSE_PTS):
        T = T_INF_0 + th[:, j]
        assert float(np.max(T)) < T_REF, (pt, float(np.max(T)))
        dose, ref, rel, swap, miss = _q_check(T)
        _p(
            f"MAT21-12s (Q-s) {HIST_LABELS[j]}: series CEM43 {dose:.15e} s vs per-segment quad "
            f"{ref:.15e} s: rel {rel:.3e} (bar {DOSE_BAR:g}, ASSERTED; predicted <= {PRED_Q:g}); "
            f"theta_s(t*) {th[-1, j]:.6f} K; max T {np.max(T):.9f} C < 43 ASSERTED"
        )
        _p(
            f"MAT21-12s (Rc'-s) {HIST_LABELS[j]}: R swap {swap:.6e} s, miss {miss:.4e} relative "
            f"(asserted > {DOSE_BAR:g}); factor {swap / ref:.4e} (predicted range "
            f"[{PRED_RC[0]:.0e}, {PRED_RC[1]:.0e}], printed)"
        )
        assert rel <= DOSE_BAR, (pt, rel)
        assert miss > DOSE_BAR, (pt, miss)


def test_win_s_six_minute_window():
    th = _hist()
    for j, pt in enumerate(DOSE_PTS):
        mx, t_mx, lo, hi, ref, rel = _win_check(th[:, j])
        _p(
            f"MAT21-12s (Win-s) {HIST_LABELS[j]}: W {float(W_B2):g} s: max six-minute mean "
            f"{mx:.15e} K at t_max {t_mx:.6f} s (predicted t* {T_STAR:.3f}); (i) "
            f"[{T_GRID[0] + W_B2:.3f}, {T_GRID[-1]:.3f}] ASSERTED; (ii) [{lo:.12e}, {hi:.12e}] "
            f"ASSERTED; (iii) quad {ref:.15e}, rel {rel:.3e} (bar {AVG_BAR:g}, ASSERTED; "
            f"predicted <= {PRED_WIN_III:g})"
        )
        assert T_GRID[0] + W_B2 <= t_mx <= T_GRID[-1], t_mx
        assert lo <= mx <= hi, (lo, mx, hi)
        assert rel <= AVG_BAR, (pt, rel)


def test_a2_sphere_control_reproduces_mat17_record():
    T0, T1 = _T_centre([0.0, T_END])
    assert T0 < 43.0 < T1, "brentq bracket"
    t43 = brentq(
        lambda s: _T_centre(s)[0] - 43.0, 0.0, T_END, xtol=1e-12, rtol=4 * np.finfo(float).eps
    )

    def f(s):
        T = _T_centre(s)[0]
        return (0.5 if T >= 43.0 else 0.25) ** (43.0 - T)

    i1, _e1 = quad(f, 0.0, t43, epsabs=1e-13, epsrel=1e-13, limit=2000)
    i2, _e2 = quad(f, t43, T_END, epsabs=1e-13, epsrel=1e-13, limit=2000)
    ref = i1 + i2
    t_fine = np.linspace(0.0, T_END, 2**19 + 1)
    d18 = cem43(t_fine[::2], _T_centre(t_fine)[::2])
    err18 = _rel(d18, ref)
    got = {
        "t43": f"{t43:.10f}",
        "quad": f"{ref:.15e}",
        "op262144": f"{d18:.15e}",
        "rel": f"{err18:.3e}",
    }
    _p(
        f"MAT21-12s (A2) sphere control: t43 {got['t43']} quad {got['quad']} op(N=2^18) "
        f"{got['op262144']} rel {got['rel']} (bar {DOSE_BAR:g}, ASSERTED); mat:17 record "
        f"{MAT17_REC_A2} (strings ASSERTED)"
    )
    assert err18 <= DOSE_BAR, err18
    assert got == MAT17_REC_A2, got


def test_printed_true_series_dose_and_records():
    """Printed records; asserts only that the readings were formed."""
    th = _hist()
    _p(
        f"MAT21-12s +y / -y series histories np.array_equal: "
        f"{bool(np.array_equal(th[:, 0], th[:, 1]))} (predicted True, printed)"
    )
    for j, pt in enumerate(HIST_PTS):
        T = T_INF_0 + th[:, j]
        op = cem43(T_GRID, T)
        true, err = _true_dose(pt)
        kind = "centre" if j == 2 else "wall"
        pred = PRED_FLOOR[kind] if j < 3 else float("nan")
        mx, t_mx, *_r = _win_check(th[:, j])
        _p(
            f"MAT21-12s SERIES DOSE {HIST_LABELS[j]}: 41-sample cem43 {op:.15e} s; true series "
            f"dose {true:.15e} s (quad err est {err:.1e}); sampling floor op/true - 1 "
            f"{op / true - 1.0:+.4e} (predicted {pred:+.2e}); six-minute peak {mx:.6f} K at "
            f"{t_mx:.3f} s; theta_s(t*) {th[-1, j]:.6f} K"
        )
        assert np.isfinite(op) and np.isfinite(true)
