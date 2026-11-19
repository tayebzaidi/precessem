"""`IMP-8` step 7 (§9 item 808, §10 I62, provisional): the two-tip source under a
duty cycle, closed form, numpy / scipy only.

Brief: docs/planning/chunks/IMP-8.md, "§9 item 808 — 2026-11-03".  Item 787's
duty history (`test_imp8_rod_line_source_duty_closed_form.py`, imported, never
edited) with the kernel swapped: step 1's Φ₂ (`_phi_tips`, the tip pair of the
same total power, balls at z = ±L/2), θ_step = q′ / (4πk) · Φ₂ with q′ = P / L and
P = 1 W.  Item 787's step sum and convolution are bound to Φ, so they are
wrapped here with the kernel as an argument; its schedule (`_edges`), grid,
duties, points and scale are imported.

Asserted (bands imported, unmoved):
- (D1) the D = 1 sum equals θ_step within ``D1_BAND`` of θ_ss.
- (Duh) the step sum against the impulse-kernel convolution within ``DUH_BAND``
  of θ_ss.  The time derivative of Φ₂ is closed: each tip term is
  (L/2)·erfc(R/2√(κt))/R, whose t-derivative is (L/2)·e^{−R²/4κt}/(2√(πκ) t^{3/2}).
- (Q-d) ``cem43`` against its quadrature within ``DOSE_BAR``; (Win-d)
  ``running_average`` (360 s) within ``AVG_BAR`` (`_q_check_g` / `_win_check_g`).
  Temperature scale (2026-11-04 03:00 review ruling, item 811): the dose is
  evaluated at P_DOSE_W = 0.5 W, T = 20 °C + 0.5 · θ(1 W), max T < 43 °C asserted
  unchanged (it guards `_q_check_g`'s R_below-only quadrature).  Item 787's line
  doses are at 1 W, so the doses are not comparable in absolute terms; θ, the
  six-minute window and tips / line stay per watt.
- Control (asserted at factor 1, rule (e)): the continuous Φ₂ form misses the
  D = 0.5 history at t* by > ``INTERIOR_BAND`` at both points.  Item 787's
  measured ≥ 1.072 is on Φ, not Φ₂, so it backs no larger factor; the brief's
  0.1 is printed as *predicted*, not asserted.

Printed (🟡 record, no band): θ(t*) per watt at both points and three duties,
beside item 787's line-source values, and tips / line per duty; (a / R₊)²
beside (2, 15) mm (1 mm from the top ball's surface).

Scope: closed form only.  Not claimed: any FEM (I56 / I60), entry 48's rod.
"""

from __future__ import annotations

import time

import numpy as np
from mpi4py import MPI
from scipy.integrate import quad

from tests.mesh.test_birdcage_implant_rod import ROD_RADIUS
from tests.validation.test_imp8_rod_line_source_closed_form import KAPPA, MM, QUAD, L, _phi_tips
from tests.validation.test_imp8_rod_line_source_duty_closed_form import (
    DUTIES,
    PERIOD,
    POINTS,
    SCALE,
    T_GRID,
    T_STAR,
    _edges,
)
from tests.validation.test_imp8_rod_line_source_duty_closed_form import _hist as _hist_line
from tests.validation.test_mat13_transient_perfusion_sphere import INTERIOR_BAND
from tests.validation.test_mat21_box_closed_form import T_INF as T_INF_0
from tests.validation.test_mat21_box_dose_closed_form import AVG_BAR, DOSE_BAR, W_B2
from tests.validation.test_mat21_box_dt_duty125_dose_closed_form import _q_check_g, _win_check_g
from tests.validation.test_mat21_wall_duty_closed_form import D1_BAND, DUH_BAND

#: The power that sets the dose's temperature scale (review ruling 2026-11-04 03:00,
#: docs/planning/chunks/IMP-8.md, item 808 BLOCKED): at 1 W the measured peak
#: θ_tips(t*) = 30.00947 K/W at D = 1, (2, 15) mm gives T = 50.009 °C, past the
#: 43 °C precondition; at 0.5 W the predicted max T is 35.005 °C.
P_DOSE_W = 0.5
MAX_T_PREDICTED = 35.005  # °C, the ruling's predicted max T at P_DOSE_W, printed only
CTRL_PREDICTED = 0.1  # the brief's predicted factor, printed only (rule (e))
#: Item 787's printed θ(t*) per watt for the line, (2, 15) / (10, 0) mm, by duty
#: (`20261103T211919Z_IMP-8-item787-W1.log`), printed beside this record only.
LINE_ITEM787 = {1.0: (9.438477, 4.645356), 0.5: (2.096145, 2.241722), 0.25: (0.8986021, 1.00437)}


def _p(*args):
    if MPI.COMM_WORLD.rank == 0:
        print(*args, flush=True)


def _step(kernel, r, z, tau):
    return SCALE * kernel(r, z, tau) if tau > 0.0 else 0.0


def _duhamel(kernel, r, z, t, duty):
    """Σⱼ (−1)ʲ θ_step(t − tⱼ) over the switch instants with tⱼ < t (item 787's rule)."""
    return sum(s * _step(kernel, r, z, t - te) for te, s in _edges(duty) if t > te)


def _phi_tips_dot(r, z, u, length=L, signs=(1.0, -1.0)):
    """∂Φ₂/∂t at age u, closed form."""
    if u <= 0.0:
        return 0.0
    out = 0.0
    for s in signs:
        big_r = np.hypot(r, z - s * 0.5 * length)
        out += np.exp(-big_r * big_r / (4.0 * KAPPA * u)) / (2.0 * np.sqrt(np.pi * KAPPA) * u**1.5)
    return 0.5 * length * out


def _convolution(kernel_dot, r, z, t, duty):
    """From rest: ∫₀ᵗ s(τ) θ̇_step(t − τ) dτ, one `quad` per on-interval (item 787's rule)."""
    total, j = 0.0, 0
    while j * PERIOD < t:
        a, b = j * PERIOD, min(j * PERIOD + duty * PERIOD, t)
        if b > a:
            total += quad(lambda u: kernel_dot(r, z, u), t - b, t - a, **QUAD)[0]
        j += 1
    return SCALE * total


_CACHE = {}


def _hist(duty, pt):
    key = (duty, pt)
    if key not in _CACHE:
        r, z = pt[0] * MM, pt[1] * MM
        _CACHE[key] = np.array([_duhamel(_phi_tips, r, z, t, duty) for t in T_GRID])
    return _CACHE[key]


def _theta_ss(pt):
    return SCALE * _phi_tips(pt[0] * MM, pt[1] * MM, np.inf)


def test_schedule_and_kernel_derivative():
    dt = T_GRID[1] - T_GRID[0]
    for duty in DUTIES:
        for te, _s in _edges(duty):
            assert np.min(np.abs(T_GRID - te)) <= 1e-9 * dt, (duty, te)
    assert T_GRID[-1] == T_STAR and W_B2 == 360.0
    # The closed ∂Φ₂/∂t against a central difference of Φ₂ (a check on the derivation).
    worst = 0.0
    for pt in POINTS:
        r, z = pt[0] * MM, pt[1] * MM
        for u in (30.0, 225.0, 900.0):
            h = 1e-3 * u
            fd = (_phi_tips(r, z, u + h) - _phi_tips(r, z, u - h)) / (2.0 * h)
            worst = max(worst, abs(fd / _phi_tips_dot(r, z, u) - 1.0))
    r_p = np.hypot(POINTS[0][0] * MM, POINTS[0][1] * MM - 0.5 * L)
    _p(
        f"\n[imp8-7] grid {T_GRID.size} pts, dt {dt:g} s, P {PERIOD:g} s; scale {SCALE:.9e} K per "
        f"unit Φ₂ at 1 W; dΦ₂/dt closed vs central difference worst rel {worst:.3e} (printed); "
        f"(a/R+)^2 at {POINTS[0]} mm = {(ROD_RADIUS / r_p) ** 2:.3e}; bands D1 {D1_BAND:g} "
        f"DUH {DUH_BAND:g} DOSE {DOSE_BAR:g} AVG {AVG_BAR:g}"
    )


def test_d1_duty_one_equals_step_response():
    worst = 0.0
    for pt in POINTS:
        r, z = pt[0] * MM, pt[1] * MM
        step = np.array([_step(_phi_tips, r, z, t) for t in T_GRID])
        e = float(np.max(np.abs(_hist(1.0, pt) - step))) / _theta_ss(pt)
        _p(f"[imp8-7] (D1) {pt} mm: max |D=1 sum - θ_step| / θ_ss = {e:.3e} (band {D1_BAND:g})")
        worst = max(worst, e)
    assert worst <= D1_BAND


def test_duh_step_sum_equals_convolution():
    t0 = time.perf_counter()
    worst = 0.0
    for duty in DUTIES:
        for pt in POINTS:
            r, z = pt[0] * MM, pt[1] * MM
            conv = np.array([_convolution(_phi_tips_dot, r, z, t, duty) for t in T_GRID])
            e = float(np.max(np.abs(_hist(duty, pt) - conv))) / _theta_ss(pt)
            _p(
                f"[imp8-7] (Duh) D={duty:g} {pt} mm: max |step sum - convolution| / θ_ss = "
                f"{e:.3e} (band {DUH_BAND:g}); θ(t*) {_hist(duty, pt)[-1]:.9e} vs {conv[-1]:.9e}"
            )
            worst = max(worst, e)
    _p(f"[imp8-7] (Duh) worst {worst:.3e}; {time.perf_counter() - t0:.1f} s")
    assert worst <= DUH_BAND


def test_q_d_win_d_and_printed_record():
    for duty in DUTIES:
        for i, pt in enumerate(POINTS):
            th = _hist(duty, pt)
            T = T_INF_0 + P_DOSE_W * th
            assert float(np.max(T)) < 43.0, (duty, pt, float(np.max(T)))
            dose, ref, rel, _swap, _miss = _q_check_g(T, T_GRID)
            mx, t_mx, lo, hi, wref, wrel = _win_check_g(th, T_GRID)
            line = float(_hist_line(duty, pt)[-1])
            _p(
                f"[imp8-7] (Q-d) D={duty:g} {pt} mm: CEM43 {dose:.9e} s vs quad {ref:.9e} s rel "
                f"{rel:.3e} (bar {DOSE_BAR:g}); dose at P = {P_DOSE_W:g} W (item 787's line: 1 W); "
                f"max T {np.max(T):.6f} C < 43 (predicted max over all cases {MAX_T_PREDICTED:.3f} C)"
            )
            _p(
                f"[imp8-7] (Win-d) D={duty:g} {pt} mm: six-minute max {mx:.9e} K at {t_mx:g} s, "
                f"in [{lo:.6e}, {hi:.6e}]; quad {wref:.9e} rel {wrel:.3e} (bar {AVG_BAR:g})"
            )
            _p(
                f"[imp8-7] (record) D={duty:g} {pt} mm: θ_tips(t*) {th[-1]:.7g} K/W; θ_ss "
                f"{_theta_ss(pt):.7g} K/W; six-minute max {mx:.7g} K/W; CEM43 {dose:.7g} s; "
                f"θ_line(t*) {line:.7g} K/W (item 787 printed {LINE_ITEM787[duty][i]:.7g}); "
                f"tips / line {th[-1] / line:.7g}"
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
            f"[imp8-7] (c-cont) {pt} mm: θ_D1(t*) {d1:.7g} vs θ_D0.5(t*) {d5:.7g} K/W, miss "
            f"{miss:.4e} (asserted > INTERIOR_BAND {INTERIOR_BAND:g}; predicted {CTRL_PREDICTED:g}, "
            f"not asserted)"
        )
        least = min(least, miss)
    assert least > INTERIOR_BAND
