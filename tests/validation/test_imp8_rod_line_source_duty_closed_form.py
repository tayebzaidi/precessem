"""`IMP-8` step 4 (§9 item 787, §10 I59): the finite line source under a duty
cycle, closed form, numpy / scipy only.

Brief: docs/planning/chunks/IMP-8.md, "§9 item 787 — 2026-11-03".  Step 0's
kernel Φ (`test_imp8_rod_line_source_closed_form._phi`, imported, never
restated) gives the step response per watt of total power,
θ_step(r, z, t) = q′ / (4πk) · Φ(r, z, t) with q′ = P / L and P = 1 W.  The duty
history is `MAT-17`'s Duhamel rule (`test_mat17_duty_cycle_closed_form.py:92`)
re-implemented on θ_step: on at jP, off at jP + DP, P = 0.5 t* = 450 s,
D ∈ {1, 0.5, 0.25}, read at (2, 15) and (10, 0) mm on the 41-point grid to
t* = 900 s.  Every switch time lies on the grid (asserted).

Asserted:
- (D1) the D = 1 sum (every off edge cancelling the next on edge) equals
  θ_step within ``D1_BAND``, relative to θ_ss at the point.
- (Duh) `MAT-21` step 21a's identity (`test_mat21_wall_duty_closed_form.py:119`,
  "from rest = Duhamel step superposition") applied to Φ.  Φ has no modal duty
  form, so the from-rest route is the convolution of the schedule with the
  impulse kernel, θ(t) = ∫₀ᵗ s(τ) ∂θ_step/∂t(t − τ) dτ, where the z′ integral
  of ∂Φ/∂t is closed: ∂Φ/∂t = e^{−r²/4κt} / (2t) · [erf((z + L/2)/2√(κt))
  − erf((z − L/2)/2√(κt))].  One `quad` in τ per on-interval, against the
  step sum (one `quad` in z′ per edge).  Within ``DUH_BAND`` of θ_ss.
- (Q-d) ``cem43`` against its per-segment quadrature within ``DOSE_BAR``, and
  (Win-d) ``running_average`` (window 360 s) against its quadrature within
  ``AVG_BAR`` (`_q_check_g` / `_win_check_g`, imported).  The temperature
  scale is `MAT-21` step 12's (P39′): T = T_INF + θ with T_INF = 20 °C
  (MAT-15's, the box's), θ here at 1 W.  max T < 43 °C asserted, as there
  (the quadrature reference integrates R_below only).
- Control: the continuous history (D = 1) substituted for D = 0.5 misses the
  D = 0.5 value at t* by > 0.1 relative (the brief's).

Printed, as a record (🟡): θ(t*) per watt, the six-minute-window maximum and
the dose at each duty and point.

Scope: closed form only.  Not claimed: the FEM (I60), a recommended duty.
"""

from __future__ import annotations

import time

import numpy as np
from mpi4py import MPI
from scipy import special
from scipy.integrate import quad

from tests.validation.test_imp8_rod_line_source_closed_form import (
    KAPPA,
    MM,
    QUAD,
    L,
    _phi,
    _phi_steady,
)
from tests.validation.test_mat15_coil_driven_heating_phantom import K_GEL
from tests.validation.test_mat21_box_closed_form import T_INF as T_INF_0
from tests.validation.test_mat21_box_dose_closed_form import AVG_BAR, DOSE_BAR, W_B2
from tests.validation.test_mat21_box_dt_duty125_dose_closed_form import _q_check_g, _win_check_g
from tests.validation.test_mat21_wall_duty_closed_form import D1_BAND, DUH_BAND

P_W = 1.0  # W, total power (per watt)
T_STAR = 900.0  # s
PERIOD = 0.5 * T_STAR  # 450 s
DUTIES = (1.0, 0.5, 0.25)
POINTS = ((2.0, 15.0), (10.0, 0.0))  # (r, z), mm
T_GRID = np.linspace(0.0, T_STAR, 41)  # Δt = 22.5 s
CTRL_MIN = 0.1  # the brief's control threshold
SCALE = P_W / (4.0 * np.pi * K_GEL * L)  # q′ / (4πk), K per unit Φ


def _p(*args):
    if MPI.COMM_WORLD.rank == 0:
        print(*args, flush=True)


def _edges(duty, t_end=T_STAR):
    """Switch instants with signs: on at jP (+1), off at jP + DP (−1), ≤ t_end."""
    out, j = [], 0
    while j * PERIOD <= t_end:
        out.append((j * PERIOD, 1.0))
        if j * PERIOD + duty * PERIOD <= t_end:
            out.append((j * PERIOD + duty * PERIOD, -1.0))
        j += 1
    return out


def _theta_step(r, z, tau):
    return SCALE * _phi(r, z, tau) if tau > 0.0 else 0.0


def _duhamel(r, z, t, duty):
    """Σⱼ (−1)ʲ θ_step(t − tⱼ) over the switch instants with tⱼ < t."""
    return sum(s * _theta_step(r, z, t - te) for te, s in _edges(duty) if t > te)


def _phi_dot(r, z, u):
    """∂Φ/∂t at age u, the z′ integral done in closed form."""
    if u <= 0.0:
        return 0.0
    w = 2.0 * np.sqrt(KAPPA * u)
    lo, hi = (z - 0.5 * L) / w, (z + 0.5 * L) / w
    return np.exp(-r * r / (4.0 * KAPPA * u)) / (2.0 * u) * (special.erf(hi) - special.erf(lo))


def _convolution(r, z, t, duty):
    """From rest: ∫₀ᵗ s(τ) θ̇_step(t − τ) dτ, one `quad` per on-interval."""
    total, j = 0.0, 0
    while j * PERIOD < t:
        a, b = j * PERIOD, min(j * PERIOD + duty * PERIOD, t)
        if b > a:
            total += quad(lambda u: _phi_dot(r, z, u), t - b, t - a, **QUAD)[0]
        j += 1
    return SCALE * total


_CACHE = {}


def _hist(duty, pt):
    key = (duty, pt)
    if key not in _CACHE:
        r, z = pt[0] * MM, pt[1] * MM
        _CACHE[key] = np.array([_duhamel(r, z, t, duty) for t in T_GRID])
    return _CACHE[key]


def _theta_ss(pt):
    return SCALE * _phi_steady(pt[0] * MM, pt[1] * MM)


def test_schedule_on_grid():
    dt = T_GRID[1] - T_GRID[0]
    for duty in DUTIES:
        for te, _s in _edges(duty):
            assert np.min(np.abs(T_GRID - te)) <= 1e-9 * dt, (duty, te)
    _p(
        f"[imp8-4] grid {T_GRID.size} pts, dt {dt:g} s, P {PERIOD:g} s; switches "
        + "; ".join(
            f"D={d:g}: " + ",".join(f"{te:g}{'+' if s > 0 else '-'}" for te, s in _edges(d))
            for d in DUTIES
        )
        + f"; scale q'/(4πk) {SCALE:.9e} K per unit Φ at {P_W:g} W; bands D1 {D1_BAND:g} "
        f"DUH {DUH_BAND:g} DOSE {DOSE_BAR:g} AVG {AVG_BAR:g}"
    )
    assert T_GRID[-1] == T_STAR and W_B2 == 360.0


def test_d1_duty_one_equals_step_response():
    worst = 0.0
    for pt in POINTS:
        r, z = pt[0] * MM, pt[1] * MM
        step = np.array([_theta_step(r, z, t) for t in T_GRID])
        e = float(np.max(np.abs(_hist(1.0, pt) - step))) / _theta_ss(pt)
        _p(f"[imp8-4] (D1) {pt} mm: max |D=1 sum - θ_step| / θ_ss = {e:.3e} (band {D1_BAND:g})")
        worst = max(worst, e)
    assert worst <= D1_BAND


def test_duh_step_sum_equals_convolution():
    t0 = time.perf_counter()
    worst = 0.0
    for duty in DUTIES:
        for pt in POINTS:
            r, z = pt[0] * MM, pt[1] * MM
            conv = np.array([_convolution(r, z, t, duty) for t in T_GRID])
            e = float(np.max(np.abs(_hist(duty, pt) - conv))) / _theta_ss(pt)
            _p(
                f"[imp8-4] (Duh) D={duty:g} {pt} mm: max |step sum - convolution| / θ_ss = "
                f"{e:.3e} (band {DUH_BAND:g}); θ(t*) {_hist(duty, pt)[-1]:.9e} vs {conv[-1]:.9e}"
            )
            worst = max(worst, e)
    _p(f"[imp8-4] (Duh) worst {worst:.3e}; {time.perf_counter() - t0:.1f} s")
    assert worst <= DUH_BAND


def test_q_d_win_d_and_printed_record():
    for duty in DUTIES:
        for pt in POINTS:
            th = _hist(duty, pt)
            T = T_INF_0 + th
            assert float(np.max(T)) < 43.0, (duty, pt, float(np.max(T)))
            dose, ref, rel, _swap, _miss = _q_check_g(T, T_GRID)
            mx, t_mx, lo, hi, wref, wrel = _win_check_g(th, T_GRID)
            _p(
                f"[imp8-4] (Q-d) D={duty:g} {pt} mm: CEM43 {dose:.9e} s vs quad {ref:.9e} s rel "
                f"{rel:.3e} (bar {DOSE_BAR:g}); max T {np.max(T):.6f} C < 43"
            )
            _p(
                f"[imp8-4] (Win-d) D={duty:g} {pt} mm: six-minute max {mx:.9e} K at {t_mx:g} s, "
                f"in [{lo:.6e}, {hi:.6e}]; quad {wref:.9e} rel {wrel:.3e} (bar {AVG_BAR:g})"
            )
            _p(
                f"[imp8-4] (record) D={duty:g} {pt} mm: θ(t*) {th[-1]:.7g} K/W; θ_ss "
                f"{_theta_ss(pt):.7g} K/W; six-minute max {mx:.7g} K/W; CEM43 {dose:.7g} s "
                f"at T = {T_INF_0:g} C + θ(1 W)"
            )
            assert rel <= DOSE_BAR, (duty, pt, rel)
            assert lo <= mx <= hi and T_GRID[0] + W_B2 <= t_mx <= T_GRID[-1], (duty, pt, t_mx)
            assert wrel <= AVG_BAR, (duty, pt, wrel)


def test_control_continuous_for_half_duty():
    least = np.inf
    for pt in POINTS:
        d1, d5 = _hist(1.0, pt)[-1], _hist(0.5, pt)[-1]
        miss = abs(d1 / d5 - 1.0)
        _p(
            f"[imp8-4] (c-cont) {pt} mm: θ_D1(t*) {d1:.7g} vs θ_D0.5(t*) {d5:.7g} K/W, miss "
            f"{miss:.4e} (asserted > {CTRL_MIN:g})"
        )
        least = min(least, miss)
    assert least > CTRL_MIN
