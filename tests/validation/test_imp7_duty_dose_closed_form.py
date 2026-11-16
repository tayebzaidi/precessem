"""`IMP-7` step 7, W1 — `MAT-20`'s dose and six-minute window around the adiabatic
hole under the duty-cycled source (§9 item 699, §10 chain I, I48a; brief
`docs/planning/chunks/IMP-7.md` "§9 item 699").

Numpy (the imported fixture modules bring dolfinx in).  Sets no band and moves
none; every band is imported.  Histories: item 675's
``pec_sphere_implant_duty_cycle_from_modes`` on ``_modes()`` at ``PTS`` on a
161-row grid to 4 t\\*: ``T_STAR`` here is 0.2 (b − a)² / κ ≈ 313.5 s, shorter than
the 360 s window, so a t\\* history could not hold it.  Every switch lies on a
grid point (``PERIOD`` = 20 ``DT``).

Asserted: (C100-i) ``duty = 1.0`` equals the continuous history at every row and
point to ``B1_BAND`` and its dose to ``DOSE_BAR``; (Q-id) ``cem43`` against
per-segment ``quad`` at peak power and at equal mean power (θ / DUTY — a scale, not
a solve) at ``DOSE_BAR``, max T < 43 °C; (Win-id) the six-minute window on the same
two scales at ``AVG_BAR``.  Negative controls: (Rc-id) the R swap misses
> ``DOSE_BAR``; (Pole-id) the isotropic source (β = 0) misses the β = 1 dose at the
pole by > ``DOSE_BAR``.  Printed: each point's D_duty(peak), D_duty(equal mean),
D_cont, the equal-mean ratio, the pole / equator ratio, the six-minute peak and
t_max, beside the bare gel box's ratio.

Not claimed: a metal core, a coil source, a thermal-damage statement, a rate.
The numbers are at `MAT-9`'s q_far and T∞ = 20 °C.
"""

from __future__ import annotations

import numpy as np

from fem_em_solver.thermal.dose import cem43, running_average
from fem_em_solver.utils.analytical import (
    pec_sphere_implant_duty_cycle_from_modes,
    pec_sphere_implant_transient_modes,
)

from tests.validation.test_imp4_pec_cylinder_closed_form import B1_BAND
from tests.validation.test_imp7_duty_cycle_fem import A, B, DUTY, PERIOD
from tests.validation.test_imp7_implant_heating_closed_form import FIX, T_INF
from tests.validation.test_imp7_transient_closed_form import N_DEFAULT, T_STAR, _modes, _theta
from tests.validation.test_mat12_transient_sar_heating_sphere import N_MAIN, RHO_C
from tests.validation.test_mat20_cem43_closed_form import BAR as DOSE_BAR
from tests.validation.test_mat20_cem43_closed_form import _rel
from tests.validation.test_mat20_running_average import BAR as AVG_BAR
from tests.validation.test_mat20_running_average import W_B2
from tests.validation.test_mat20_tuned_heating_dose import T_REF
from tests.validation.test_mat21_wall_duty_dose_closed_form import BARE_EQ_RATIO, _q_check, _win_check

DT = T_STAR / N_MAIN
T_END = 4 * T_STAR
T_GRID = DT * np.arange(4 * N_MAIN + 1)
#: Interior, so I48b can read them on the curved P2 mesh.
PTS = np.array([(0.0, 0.0, 1.05 * A), (1.05 * A, 0.0, 0.0), (0.0, 0.0, 0.95 * B)])
LABELS = ("pole", "equator", "outer pole")
POLE, EQUATOR = 0, 1
SCALES = ((1.0, "peak"), (1.0 / DUTY, "equal-mean"))


def _p(msg):
    print(msg, flush=True)


def _duty_hist(modes, duty=DUTY):
    """(161, 3) θ histories; row 0 = 0.0 exactly."""
    out = np.zeros((T_GRID.size, PTS.shape[0]))
    for n in range(1, T_GRID.size):
        out[n] = pec_sphere_implant_duty_cycle_from_modes(
            PTS, T_GRID[n], modes, rho_c=RHO_C, t_inf=T_INF, period=PERIOD, duty=duty) - T_INF
    return out


def _cont_hist(modes):
    out = np.zeros((T_GRID.size, PTS.shape[0]))
    for n in range(1, T_GRID.size):
        out[n] = _theta(PTS, T_GRID[n], modes)
    return out


_CACHE = {}


def _hists():
    if not _CACHE:
        m = _modes()
        _CACHE["cont"] = _cont_hist(m)
        _CACHE["d1"] = _duty_hist(m, duty=1.0)
        _CACHE["d"] = _duty_hist(m)
        # The control builds its own modes: `_modes()` fixes beta (brief).
        _CACHE["iso"] = _duty_hist(pec_sphere_implant_transient_modes(
            beta=0.0, n_terms=N_DEFAULT, inner="neumann", **FIX))
    return _CACHE


def test_precondition_grid():
    assert T_GRID.size == 161 and T_GRID[0] == 0.0
    assert abs(T_GRID[-1] / T_END - 1.0) <= 1e-14
    assert T_END >= W_B2, (T_END, W_B2)
    assert abs(PERIOD / DT - 20.0) < 1e-12 and abs(DUTY * PERIOD / DT - 10.0) < 1e-12
    assert np.all(_hists()["d"][0] == 0.0)
    _p(f"IMP7-7 grid: t* {T_STAR:.4f} s, T_END {T_END:.4f} s >= W {float(W_B2):g} s, dt {DT:.4f} s, "
       f"period {PERIOD:.4f} s; PTS {PTS.tolist()} ({LABELS})")


def test_c100_i_full_duty_equals_continuous():
    h = _hists()
    worst = float(np.max(np.abs(h["d1"][1:] - h["cont"][1:]) / np.abs(h["cont"][1:])))
    _p(f"IMP7-7 (C100-i) duty = 1.0 vs continuous over {h['cont'].shape}: worst row/point rel "
       f"{worst:.3e} (band {B1_BAND:g}, ASSERTED)")
    assert worst <= B1_BAND, worst
    for j, lab in enumerate(LABELS):
        d1 = cem43(T_GRID, T_INF + h["d1"][:, j])
        d0 = cem43(T_GRID, T_INF + h["cont"][:, j])
        _p(f"IMP7-7 (C100-i) {lab}: dose(100 %) {d1:.15e} s, dose(continuous) {d0:.15e} s, rel "
           f"{_rel(d1, d0):.3e} (bar {DOSE_BAR:g}, ASSERTED)")
        assert _rel(d1, d0) <= DOSE_BAR, (lab, d1, d0)


def test_q_id_win_id_and_rc_id():
    th_d = _hists()["d"]
    for scale, lab in SCALES:
        for j, name in enumerate(LABELS):
            th = th_d[:, j] * scale
            T = T_INF + th
            dose, ref, rel, swap, miss = _q_check(T_GRID, T)
            mx, t_mx, lo, hi, wref, wrel = _win_check(T_GRID, th)
            _p(f"IMP7-7 (Q-id) {lab} {name}: duty CEM43 {dose:.15e} s vs quad {ref:.15e} s rel "
               f"{rel:.3e} (bar {DOSE_BAR:g}, ASSERTED); max T {np.max(T):.9f} C < {T_REF} ASSERTED")
            _p(f"IMP7-7 (Win-id) {lab} {name}: max six-minute mean {mx:.15e} K at t_max {t_mx:.4f} s; "
               f"(i) [{T_GRID[0] + W_B2:.3f}, {T_GRID[-1]:.3f}]; (ii) [{lo:.12e}, {hi:.12e}]; (iii) "
               f"quad {wref:.15e} rel {wrel:.3e} (bar {AVG_BAR:g}, ASSERTED)")
            _p(f"IMP7-7 (Rc-id) {lab} {name}: R swap miss {miss:.4e} (asserted > {DOSE_BAR:g}); "
               f"factor {swap / ref:.4e}")
            assert float(np.max(T)) < T_REF, (lab, name, float(np.max(T)))
            assert rel <= DOSE_BAR, (lab, name, rel)
            assert T_GRID[0] + W_B2 <= t_mx <= T_GRID[-1], (lab, t_mx)
            assert lo <= mx <= hi, (lab, lo, mx, hi)
            assert wrel <= AVG_BAR, (lab, name, wrel)
            assert miss > DOSE_BAR, (lab, name, miss)


def test_pole_id_isotropic_source_control():
    h = _hists()
    d_b = cem43(T_GRID, T_INF + h["d"][:, POLE])
    d_0 = cem43(T_GRID, T_INF + h["iso"][:, POLE])
    miss = _rel(d_0, d_b)
    # 7 digits: the first log's print at .15e / .9f tripped the private-leak guard (held).
    _p(f"IMP7-7 (Pole-id) pole: beta = 1 {d_b:.6e} s, beta = 0 {d_0:.6e} s, ratio {d_b / d_0:.6f}, "
       f"miss {miss:.4e} (asserted > {DOSE_BAR:g}; separation {miss / DOSE_BAR:.3e} x the bar)")
    assert miss > DOSE_BAR, miss


def test_printed_duty_dose_records():
    """The records (no band); asserts only that they were formed."""
    h = _hists()
    W = float(W_B2)
    doses = {}
    for j, lab in enumerate(LABELS):
        d_c = cem43(T_GRID, T_INF + h["cont"][:, j])
        d_pk = cem43(T_GRID, T_INF + h["d"][:, j])
        d_eq = cem43(T_GRID, T_INF + h["d"][:, j] / DUTY)
        doses[j] = (d_c, d_pk, d_eq)
        _a, mx_c, t_c = running_average(T_GRID, h["cont"][:, j], window=W)
        _a, mx_d, t_d = running_average(T_GRID, h["d"][:, j], window=W)
        _p(f"IMP7-7 DUTY DOSE {lab}: D_cont {d_c:.15e} s; D_duty(peak) {d_pk:.15e} s; D_duty(equal "
           f"mean) {d_eq:.15e} s; equal-mean ratio {d_eq / d_c:.6f} (bare gel box {BARE_EQ_RATIO}, "
           f"its own schedule, printed); peak ratio {d_pk / d_c:.6f}; six-minute peak cont "
           f"{mx_c:.6f} K at {t_c:.2f} s, duty {mx_d:.6f} K at {t_d:.2f} s")
        assert np.isfinite(d_eq / d_c)
    for k, name in enumerate(("D_cont", "D_duty(peak)", "D_duty(equal mean)")):
        _p(f"IMP7-7 pole / equator {name}: {doses[POLE][k] / doses[EQUATOR][k]:.6f}")
