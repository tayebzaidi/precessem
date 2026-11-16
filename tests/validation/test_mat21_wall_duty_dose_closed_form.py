"""`MAT-21` step 22, W1 — `MAT-20`'s dose and six-minute window on the duty-cycled
two-layer (gel + wall) slab's closed-form histories (§9 item 698, §10 chain P's
P50a; brief `docs/planning/chunks/MAT-21.md` "§9 item 698").

Numpy (the imported fixture modules bring dolfinx in).  Sets no band and moves
none; every band is imported.  Histories: item 689's
``two_layer_slab_duty_cycle_temperature`` at ``Z_PTS`` on the 321-row grid to
t₂ = 8 t\\*, every switch on a grid point (``PERIOD`` = 20 ``DT``), item 688's
schedule and constants.

Asserted: (C100-w) ``duty = 1.0`` equals ``two_layer_slab_transient_temperature``
to ``DUH_BAND`` at every row and point, and its dose equals the continuous one to
``DOSE_BAR``; (Q-wd) on the duty history at peak power and at equal mean power
(θ_duty / DUTY — a scale, not a solve) ``cem43`` against per-segment ``quad`` at
``DOSE_BAR``, max T < 43 °C; (Win-wd) the six-minute window on the same two scales
at ``AVG_BAR``.  Negative controls: (Rc-wd) the R swap misses > ``DOSE_BAR``;
(Mw-wd) the massless-wall history at the interface differs in dose by > ``DOSE_BAR``.
Printed: each point's D_duty(peak), D_duty(equal mean), D_cont, the equal-mean ratio
beside the bare box's, the walled / massless ratio, the six-minute peak and t_max.

Not claimed: the FEM history, a safety figure (a 20 °C gel at the module's q).
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import quad

from fem_em_solver.thermal.dose import cem43, running_average
from fem_em_solver.utils.analytical import (
    two_layer_slab_duty_cycle_temperature,
    two_layer_slab_transient_temperature,
)

from tests.validation.test_imp7_duty_cycle_fem import DUTY
from tests.validation.test_mat15_coil_driven_heating_phantom import N_MAIN
from tests.validation.test_mat15_finite_cylinder_closed_form import H, Q, T_INF
from tests.validation.test_mat20_cem43_closed_form import BAR as DOSE_BAR
from tests.validation.test_mat20_cem43_closed_form import _rel
from tests.validation.test_mat20_running_average import BAR as AVG_BAR
from tests.validation.test_mat20_running_average import W_B2
from tests.validation.test_mat20_tuned_heating_dose import R_BELOW, T_REF, _quad_dose
from tests.validation.test_mat21_box_closed_form import C_BOX
from tests.validation.test_mat21_box_dose_closed_form import LN_R
from tests.validation.test_mat21_wall_closed_form import RHO_C_W
from tests.validation.test_mat21_wall_duty_closed_form import DUH_BAND, PERIOD
from tests.validation.test_mat21_wall_transient_closed_form import N_SERIES, T_STAR
from tests.validation.test_mat21_wall_transient_fem import MASSLESS, N2, SLAB, _lam

DT = T_STAR / N_MAIN
T_GRID = DT * np.arange(N2 + 1)
Z_PTS = (0.0, 0.9 * C_BOX, C_BOX)
LABELS = ("centre", "near-wall", "interface")
INTERFACE = 2
# The bare box's equal-mean dose ratios (P38's schedule, printed only), item 567's W1:
# docs/testing/logs/20261020T140919Z_MAT-21-step13-W1-duty-dose-closed-form.log:93-96.
BARE_EQ_RATIO = {"centre": 1.004475, "(0,0.9B,0)": 1.213031, "(0,+B,0)": 1.186312}
SCALES = ((1.0, "peak"), (1.0 / DUTY, "equal-mean"))


def _p(msg):
    print(msg, flush=True)


def _duty_hist(rho_c_w, duty=DUTY):
    """(321, 3) θ histories; row 0 = 0.0 exactly (`_duty_series`,
    `test_mat21_wall_duty_fem.py:72`, with ``duty`` and ``rho_c_w`` as arguments)."""
    z = np.asarray(Z_PTS, dtype=float)
    lam = _lam(rho_c_w)
    out = np.zeros((T_GRID.size, z.size))
    for n in range(1, T_GRID.size):
        out[n] = two_layer_slab_duty_cycle_temperature(
            z, T_GRID[n], q=Q, h=H, t_inf=0.0, rho_c_w=rho_c_w, period=PERIOD, duty=duty,
            n_terms=N_SERIES, eigenvalues=lam, **SLAB)
    return out


def _cont_hist():
    z = np.asarray(Z_PTS, dtype=float)
    lam = _lam(RHO_C_W)
    out = np.zeros((T_GRID.size, z.size))
    for n in range(1, T_GRID.size):
        out[n] = two_layer_slab_transient_temperature(
            z, T_GRID[n], q=Q, h=H, t_inf=0.0, rho_c_w=RHO_C_W, n_terms=N_SERIES, eigenvalues=lam,
            **SLAB)
    return out


def _q_check(t, T):
    """`_q_check` of `test_mat21_box_dose_closed_form.py:78`, ``t`` an argument
    (that one closes over the bare box's ``T_GRID``)."""
    dose = cem43(t, T)
    ref = _quad_dose(t, T, LN_R)
    swap = cem43(t, T, r_above=R_BELOW, r_below=0.5)
    return dose, ref, _rel(dose, ref), swap, _rel(swap, ref)


def _win_check(t, th, window=float(W_B2)):
    """`_win_check` of `test_mat21_box_dose_closed_form.py:86`, ``t`` an argument."""
    _means, mx, t_mx = running_average(t, th, window=window)
    s = t_mx - window
    inside = t[(t > s) & (t <= t_mx)]
    pts = np.concatenate([[np.interp(s, t, th)], np.interp(inside, t, th)])
    lo, hi = float(np.min(pts)), float(np.max(pts))
    brk = [x for x in t if s < x < t_mx]
    ref, _err = quad(lambda x: np.interp(x, t, th), s, t_mx, points=brk or None,
                     epsabs=0.0, epsrel=1e-13, limit=200)
    ref /= window
    return float(mx), float(t_mx), lo, hi, ref, _rel(mx, ref)


_CACHE = {}


def _hists():
    if not _CACHE:
        _CACHE["cont"] = _cont_hist()
        _CACHE["d1"] = _duty_hist(RHO_C_W, duty=1.0)
        _CACHE["d"] = _duty_hist(RHO_C_W)
        _CACHE["m"] = _duty_hist(MASSLESS * RHO_C_W)
    return _CACHE


def test_precondition_grid():
    assert T_GRID.size == 321 and T_GRID[0] == 0.0
    assert abs(T_GRID[-1] / (8.0 * T_STAR) - 1.0) <= 1e-14
    # Every switch on a grid point: a period is 20 DT, an on-phase 10.
    assert abs(PERIOD / DT - 20.0) < 1e-12 and abs(DUTY * PERIOD / DT - 10.0) < 1e-12
    assert np.all(_hists()["d"][0] == 0.0)


def test_c100_w_full_duty_equals_transient():
    h = _hists()
    worst = float(np.max(np.abs(h["d1"][1:] - h["cont"][1:]) / np.abs(h["cont"][1:])))
    _p(f"MAT21-22 (C100-w) duty = 1.0 vs transient over {h['cont'].shape}: worst row/point rel "
       f"{worst:.3e} (band {DUH_BAND:g}, ASSERTED)")
    assert worst <= DUH_BAND, worst
    for j, lab in enumerate(LABELS):
        d1 = cem43(T_GRID, T_INF + h["d1"][:, j])
        d0 = cem43(T_GRID, T_INF + h["cont"][:, j])
        _p(f"MAT21-22 (C100-w) {lab}: dose(100 %) {d1:.15e} s, dose(continuous) {d0:.15e} s, rel "
           f"{_rel(d1, d0):.3e} (bar {DOSE_BAR:g}, ASSERTED)")
        assert _rel(d1, d0) <= DOSE_BAR, (lab, d1, d0)


def test_q_wd_win_wd_and_rc_wd():
    th_d = _hists()["d"]
    for scale, lab in SCALES:
        for j, name in enumerate(LABELS):
            th = th_d[:, j] * scale
            T = T_INF + th
            dose, ref, rel, swap, miss = _q_check(T_GRID, T)
            mx, t_mx, lo, hi, wref, wrel = _win_check(T_GRID, th)
            _p(f"MAT21-22 (Q-wd) {lab} {name}: duty CEM43 {dose:.15e} s vs quad {ref:.15e} s rel "
               f"{rel:.3e} (bar {DOSE_BAR:g}, ASSERTED); max T {np.max(T):.9f} C < {T_REF} ASSERTED")
            _p(f"MAT21-22 (Win-wd) {lab} {name}: max six-minute mean {mx:.15e} K at t_max {t_mx:.4f} s; "
               f"(i) [{T_GRID[0] + W_B2:.3f}, {T_GRID[-1]:.3f}]; (ii) [{lo:.12e}, {hi:.12e}]; (iii) "
               f"quad {wref:.15e} rel {wrel:.3e} (bar {AVG_BAR:g}, ASSERTED)")
            _p(f"MAT21-22 (Rc-wd) {lab} {name}: R swap miss {miss:.4e} (asserted > {DOSE_BAR:g}); "
               f"factor {swap / ref:.4e}")
            assert float(np.max(T)) < T_REF, (lab, name, float(np.max(T)))
            assert rel <= DOSE_BAR, (lab, name, rel)
            assert T_GRID[0] + W_B2 <= t_mx <= T_GRID[-1], (lab, t_mx)
            assert lo <= mx <= hi, (lab, lo, mx, hi)
            assert wrel <= AVG_BAR, (lab, name, wrel)
            assert miss > DOSE_BAR, (lab, name, miss)


def test_mw_wd_massless_wall_control():
    h = _hists()
    d_w = cem43(T_GRID, T_INF + h["d"][:, INTERFACE])
    d_m = cem43(T_GRID, T_INF + h["m"][:, INTERFACE])
    miss = _rel(d_m, d_w)
    _p(f"MAT21-22 (Mw-wd) interface: walled {d_w:.15e} s, massless wall {d_m:.15e} s, ratio "
       f"{d_w / d_m:.9f}, miss {miss:.4e} (asserted > {DOSE_BAR:g}; separation "
       f"{miss / DOSE_BAR:.3e} x the bar)")
    assert miss > DOSE_BAR, miss


def test_printed_duty_dose_records():
    """The records (no band); asserts only that they were formed."""
    h = _hists()
    W = float(W_B2)
    for j, lab in enumerate(LABELS):
        d_c = cem43(T_GRID, T_INF + h["cont"][:, j])
        d_pk = cem43(T_GRID, T_INF + h["d"][:, j])
        d_eq = cem43(T_GRID, T_INF + h["d"][:, j] / DUTY)
        d_m = cem43(T_GRID, T_INF + h["m"][:, j])
        _a, mx_c, t_c = running_average(T_GRID, h["cont"][:, j], window=W)
        _a, mx_d, t_d = running_average(T_GRID, h["d"][:, j], window=W)
        _p(f"MAT21-22 DUTY DOSE {lab}: D_cont {d_c:.15e} s; D_duty(peak) {d_pk:.15e} s; D_duty(equal "
           f"mean) {d_eq:.15e} s; equal-mean ratio {d_eq / d_c:.6f} (bare box {BARE_EQ_RATIO}, its own "
           f"schedule, printed); peak ratio {d_pk / d_c:.6f}; walled / massless {d_pk / d_m:.6f}; "
           f"six-minute peak cont {mx_c:.6f} K at {t_c:.2f} s, duty {mx_d:.6f} K at {t_d:.2f} s")
        assert np.isfinite(d_eq / d_c)
