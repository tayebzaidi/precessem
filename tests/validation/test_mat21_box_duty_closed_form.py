"""`MAT-21` step 2, closed form — the F2182 gel box under `MAT-17` step 3's duty
cycle re-timed to the box (§9 item 424, §10 P26; brief
`docs/planning/chunks/MAT-21.md`).

Schedule: t* = 0.2 t_C, N = 40 CN steps, DT_BOX = t*/40, period P = t*/5,
D = 0.5, HALF_T = P/2 = 4 DT_BOX; switches t_m = m·HALF_T, m = 0 … 9 (even on,
odd off).  The closed form is the ten-run Duhamel sum of step 0's triple
series step response, θ_duty(t) = Σ_{t > t_m} (−1)^m θ_step(t − t_m).

Asserted (bars copy item 406's (R4) literal 1e-10·θ_ss(0),
``test_mat21_box_closed_form.py:247``): (C1) the ten-run sum against the
Duhamel quadrature with χ over the on-intervals; (C2) the exact volume mean
of the truncated series against its Duhamel quadrature at Fo ∈ FOS.
Control (rule (e), predicted, only the miss asserted): (C1-ctl) m = 9 dropped.
Printed records only: the duty history, the six-minute maximum and its
position (``running_average``), the dose per W (``cem43``).

Dolfinx-free by design (smoke window, real build): no ``test_mat17_*`` or
step 1 import.  The shape pins below copy
``test_mat17_duty_cycle_coil_driven.py:68–71`` (STEPS_PER_HALF, N_SWITCH,
DUTY) and ``test_mat15_coil_driven_heating_phantom.py:88`` (N_MAIN = 40);
step 2's FEM module asserts they equal the originals ((P0′)).
"""

from __future__ import annotations

import numpy as np

from fem_em_solver.thermal.dose import cem43, running_average
from fem_em_solver.utils.analytical import robin_box_transient_temperature
from tests.validation.test_mat21_box_closed_form import (
    A_BOX,
    ALPHA,
    B_BOX,
    C_BOX,
    FOS,
    HALF,
    KW,
    R4_PTS,
    RHO_C,
    SERIES_TERMS,
    T_C,
    T_INF,
    Q,
    K,
    _axis_data,
    _say,
    _theta_ss0,
    _uprod,
)

# Pins (see module docstring); asserted equal to the originals in W2 (P0′).
N_STEPS, STEPS_PER_HALF, N_SWITCH, DUTY = 40, 4, 10, 0.5
T_STAR = 0.2 * T_C  # 2 821.5 s, step 1's T_STAR
DT_BOX = T_STAR / N_STEPS  # 70.5375 s
HALF_T = STEPS_PER_HALF * DT_BOX  # 282.15 s
T_SW = tuple(m * HALF_T for m in range(N_SWITCH))
V_BOX = 8.0 * A_BOX * B_BOX * C_BOX
QV = Q * V_BOX  # 49.14 W, the on-phase deposited power
W_AVG = 360.0
R_BAR = 1e-10  # item 406's (R4) literal, test_mat21_box_closed_form.py:247
SAMPLES_PER_HALF = 64


def _theta_step(x, y, z, t):
    """Step 0's series from rest, ``t_inf = 0`` (step 1's ``_trans_ex`` form)."""
    return robin_box_transient_temperature(
        x, y, z, t, half_widths=HALF, rho_c=RHO_C, t_inf=0.0, n_terms=SERIES_TERMS, **KW)


def _duty_cf(x, y, z, t, drop=()):
    """Σ_{m: t > t_m} (−1)^m θ_step(t − t_m); ``drop`` removes runs (controls)."""
    out = 0.0
    for m, tm in enumerate(T_SW):
        if t > tm and m not in drop:
            out = out + (-1.0) ** m * _theta_step(x, y, z, t - tm)
    return out


def _quad_nodes(lo, hi):
    """Step 0's rule (``test_mat21_box_closed_form.py:204–209``): 16-point
    Gauss–Legendre on 80 intervals, log-spaced from 1e-10·hi when lo = 0,
    geometric between lo and hi otherwise."""
    if lo == 0.0:
        edges = np.concatenate([[0.0], np.geomspace(1e-10 * hi, hi, 80)])
    else:
        edges = np.geomspace(lo, hi, 81)
    xg, wg = np.polynomial.legendre.leggauss(16)
    s = (0.5 * (edges[1:, None] - edges[:-1, None]) * (xg + 1.0) + edges[:-1, None]).ravel()
    ws = (0.5 * (edges[1:, None] - edges[:-1, None]) * wg).ravel()
    return s, ws


def _on_lags(t):
    """On-intervals [t_{2j}, t_{2j+1}] ∩ [0, t] mapped to lags τ = t − s."""
    out = []
    for j in range(N_SWITCH // 2):
        a, b = T_SW[2 * j], T_SW[2 * j + 1]
        if t > a:
            out.append((max(0.0, t - b), t - a))
    return out


def _duhamel_duty(x, y, z, t):
    tot = 0.0
    for lo, hi in _on_lags(t):
        s, ws = _quad_nodes(lo, hi)
        tot += np.sum(ws * _uprod(x, y, z, s, SERIES_TERMS))
    return Q / RHO_C * tot


def _mean_axes():
    out = []
    for L, n in zip(HALF, SERIES_TERMS):
        gam, c = _axis_data(L, n)
        out.append((gam / L, c * np.sin(gam) / gam))
    return out


def _box_mean_step(t):
    """Exact volume mean of the truncated step series."""
    (kx, wx), (ky, wy), (kz, wz) = _mean_axes()
    lam = kx[:, None, None] ** 2 + ky[None, :, None] ** 2 + kz[None, None, :] ** 2
    w = wx[:, None, None] * wy[None, :, None] * wz[None, None, :]
    return (Q / K) * float(np.sum(w * (-np.expm1(-ALPHA * lam * t)) / lam))


def _box_mean_duty(t):
    return sum((-1.0) ** m * _box_mean_step(t - tm) for m, tm in enumerate(T_SW) if t > tm)


def _box_mean_duhamel(t):
    s, ws = _quad_nodes(0.0, t)
    prod = np.ones_like(s)
    for kd, wd in _mean_axes():
        prod = prod * (np.exp(-ALPHA * np.outer(s, kd**2)) * wd).sum(1)
    return Q / RHO_C * float(np.sum(ws * prod))


C1_TIMES = (8.5 * HALF_T, 9 * HALF_T, T_STAR)


def test_c1_duhamel_with_chi():
    th0 = _theta_ss0()
    worst = 0.0
    for t in C1_TIMES:
        for p in R4_PTS:
            a = float(_duty_cf(*p, t))
            b = float(_duhamel_duty(*p, t))
            worst = max(worst, abs(a - b) / th0)
            _say(f"MAT21-2 (C1) t = {t:.4f} s {p}: ten-run sum {a:.12f}  Duhamel(chi) {b:.12f}  "
                 f"|d|/theta_ss(0) {abs(a - b) / th0:.3e}")
    _say(f"MAT21-2 (C1) max |sum - Duhamel|/theta_ss(0) = {worst:.3e} (bar {R_BAR:.0e}, "
         f"predicted <= 1e-13)")
    assert worst <= R_BAR


def test_c2_mean_route():
    th0 = _theta_ss0()
    worst = 0.0
    for fo in FOS:
        t = fo * T_C
        a = _box_mean_step(t)
        b = _box_mean_duhamel(t)
        worst = max(worst, abs(a - b) / th0)
        _say(f"MAT21-2 (C2) Fo {fo}: mean series {a:.12f}  Duhamel {b:.12f}  "
             f"|d|/theta_ss(0) {abs(a - b) / th0:.3e}")
    _say(f"MAT21-2 (C2) max |mean - Duhamel|/theta_ss(0) = {worst:.3e} (bar {R_BAR:.0e}, "
         f"predicted <= 1e-13)")
    assert worst <= R_BAR


def test_control_c1_m9_dropped_misses():
    th0 = _theta_ss0()
    misses = []
    for p in R4_PTS:
        a = float(_duty_cf(*p, T_STAR, drop=(9,)))
        b = float(_duhamel_duty(*p, T_STAR))
        misses.append(abs(a - b) / th0)
    _say(f"MAT21-2 control (C1-ctl) m = 9 dropped at t*: miss/theta_ss(0) per point "
         + " / ".join(f"{v:.4e}" for v in misses)
         + f" (predicted ≈ 1.12e-2 at the centre, smaller at (0.9A, 0.8B, 0.8C); "
         f"asserted > {R_BAR:.0e})")
    assert min(misses) > R_BAR


def _histories():
    """Duty histories on 64 samples per half-period (641 samples), the step
    response sampled once and shifted by 64·m samples."""
    n = SAMPLES_PER_HALF * N_SWITCH
    t = np.arange(n + 1) * (HALF_T / SAMPLES_PER_HALF)
    t[-1] = T_STAR
    step_c = np.array([float(_theta_step(0.0, 0.0, 0.0, s)) if s > 0 else 0.0 for s in t])
    step_m = np.array([_box_mean_step(s) if s > 0 else 0.0 for s in t])
    duty_c = np.zeros_like(t)
    duty_m = np.zeros_like(t)
    for m in range(N_SWITCH):
        k = SAMPLES_PER_HALF * m
        duty_c[k:] += (-1.0) ** m * step_c[: n + 1 - k]
        duty_m[k:] += (-1.0) ** m * step_m[: n + 1 - k]
    return t, duty_m, duty_c, step_m, step_c


def test_printed_records():
    for n in (32, 36, 40):
        t = n * DT_BOX
        c = float(_duty_cf(0.0, 0.0, 0.0, t))
        mn = _box_mean_duty(t)
        _say(f"MAT21-2 record duty n = {n} (t = {t:.4f} s): centre {c:.5f} K, volume mean {mn:.5f} K")
    _say(f"MAT21-2 record predicted centre 0.53622 / 0.66914 / 0.66654 K, mean 0.49931 / "
         f"0.62470 / 0.61419 K (±2e-5 K, 30-term stdlib arithmetic)")
    sc = float(_theta_step(0.0, 0.0, 0.0, T_STAR))
    sm = _box_mean_step(T_STAR)
    _say(f"MAT21-2 record step response at t*: centre {sc:.5f} K, mean {sm:.5f} K "
         f"(predicted 1.33568 / 1.23889 K)")

    t, duty_m, duty_c, _, _ = _histories()
    # Cross-check the shifted-sample history against the direct sum at the grid's end.
    _say(f"MAT21-2 record history check at t*: shifted {duty_m[-1]:.12f} vs direct "
         f"{_box_mean_duty(T_STAR):.12f} K (mean)")
    am, am_mx, am_t = running_average(t, duty_m, window=W_AVG)
    ac, ac_mx, ac_t = running_average(t, duty_c, window=W_AVG)
    ok = ~np.isnan(am)
    dm = np.diff(am[ok])
    _say(f"MAT21-2 record six-minute maximum (W = {W_AVG} s, {t.size} samples, "
         f"{int(ok.sum())} windows): mean {am_mx:.5f} K at t_max = {am_t:.4f} s, centre "
         f"{ac_mx:.5f} K at t_max = {ac_t:.4f} s (predicted 0.61684 / 0.66420 K at t* = "
         f"{T_STAR:.4f} s); window mean monotone non-decreasing on [W, t*]: "
         f"{bool(np.all(dm >= 0.0))} (min increment {dm.min():.3e} K)")
    per_w_m, per_w_c = am_mx / QV, ac_mx / QV
    _say(f"MAT21-2 record six-minute peak per W of qV = {QV:.4f} W: mean {per_w_m:.4e} K/W, "
         f"centre {per_w_c:.4e} K/W (predicted 1.2553e-2 / 1.3516e-2); per W of D qV: "
         f"{per_w_m / DUTY:.4e} / {per_w_c / DUTY:.4e} K/W")
    dm_w = cem43(t, T_INF + duty_m / QV)
    dc_w = cem43(t, T_INF + duty_c / QV)
    dm_q = cem43(t, T_INF + duty_m)
    dc_q = cem43(t, T_INF + duty_c)
    _say(f"MAT21-2 record cem43 of T_inf + theta/(qV) (per W of qV): mean {dm_w:.4e} s, centre "
         f"{dc_w:.4e} s (predicted ≈ 4.05e-11 both); at the fixture's q: mean {dm_q:.4e} s, centre "
         f"{dc_q:.4e} s (predicted 6.7e-11 / 6.9e-11); per W of D qV the history is theta/(D qV): "
         f"mean {cem43(t, T_INF + duty_m / (DUTY * QV)):.4e} s, centre "
         f"{cem43(t, T_INF + duty_c / (DUTY * QV)):.4e} s. Records only, no safety statement.")
