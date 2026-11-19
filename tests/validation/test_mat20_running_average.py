"""`MAT-20` step 2 — the trailing running average, closed-form half (§9 item 319).

Numpy only.  Input: `MAT-17` case A's periodic closed form at r = 0
(`robin_sphere_duty_cycle_temperature(..., periodic=True)`, P = 1200 s,
D = 0.5), sampled at M points per period so that the switch times
(0, DP, P, …) sit on the grid.

Asserted: (B1) for W = kP, k = 1, 2, every window mean equals D·θ_ss(0) at
1e-12 relative (the §7 row's bar; no range pre-registered, ruling 319-2);
(B2) for W = 360 s the window maximum lies within one sample of the brentq
root of θ(t) = θ(t − W) on the falling side; (Wc, rule (e), asserted) the
window W = kP + P/2 misses (B1) by more than the bar.  (Pc) printed only:
t_max under a window one sample short.  A linear ramp (exact for any W, also
a fractional one) checks the operator's end-segment arithmetic.
(B3), the FEM history at `D_BAND`, is not executed (ruling 319-3).
"""

from __future__ import annotations

import numpy as np
from mpi4py import MPI
from scipy.optimize import brentq

from fem_em_solver.thermal.dose import running_average
from fem_em_solver.utils.analytical import (
    robin_sphere_duty_cycle_temperature,
    robin_sphere_steady_temperature,
)

# `test_mat17_duty_cycle_closed_form.py:42–48` (case A).
R, K, RHO_C, Q, H = 0.1, 0.6, 4.18e6, 2000.0, 10.0
P, D = 1200.0, 0.5
TH_SS0 = float(robin_sphere_steady_temperature(0.0, radius=R, k=K, q=Q, h=H))
BAR = 1e-12
M_SEQ = (2**12, 2**13)  # samples per period: the printed M / 2M sequence
W_B2 = 360.0


def _p(*a):
    if MPI.COMM_WORLD.rank == 0:
        print(*a, flush=True)


def _theta(t):
    t = np.asarray(t, dtype=float)
    return robin_sphere_duty_cycle_temperature(
        np.zeros_like(t),
        t,
        radius=R,
        k=K,
        rho_c=RHO_C,
        q=Q,
        h=H,
        t_inf=0.0,
        period=P,
        duty=D,
        periodic=True,
    )


def _series(M, n_periods):
    t = np.arange(n_periods * M + 1) * (P / M)
    return t, _theta(t)


def test_linear_ramp_exact_fractional_window():
    t = np.linspace(0.0, 1000.0, 87)  # Δt = 11.63 s; W = 360 is 30.96 samples
    v = 3.0 + 0.01 * t
    for W in (360.0, 360.0 * np.pi / 3.0, 100.0):
        means, mx, t_mx = running_average(t, v, window=W)
        ok = t - W >= 0.0
        ref = 3.0 + 0.01 * (t[ok] - 0.5 * W)
        err = float(np.max(np.abs(means[ok] - ref) / ref))
        _p(f"[MAT-20 ramp] W={W:.4f}: max rel err {err:.3e}; t_max {t_mx} (band 1e-14)")
        assert err <= 1e-14
        assert np.all(np.isnan(means[~ok]))
        assert t_mx == t[-1]
    # 2-D (DG0-shaped) call matches the column-wise 1-D calls exactly.
    V = np.stack([v, 2.0 * v, -v], axis=1)
    m2, mx2, tm2 = running_average(t, V, window=360.0)
    for c in range(3):
        m1, x1, tm1 = running_average(t, V[:, c], window=360.0)
        assert np.array_equal(m2[:, c], m1, equal_nan=True)
        assert mx2[c] == x1 and tm2[c] == tm1


def test_b1_integer_period_mean_and_wc_control():
    target = D * TH_SS0
    for k in (1, 2):
        for M in M_SEQ:
            t, v = _series(M, k + 1)
            means, _, _ = running_average(t, v, window=k * P)
            ok = ~np.isnan(means)
            res = float(np.max(np.abs(means[ok] - target))) / target
            _p(
                f"[MAT-20 B1] k={k} M={M}: {ok.sum()} windows, "
                f"max |mean − D·θ_ss|/(D·θ_ss) = {res:.3e} (bar {BAR:g})"
            )
            assert res <= BAR
            # (Wc) W = kP + P/2: a non-integer number of periods.
            Wc = k * P + 0.5 * P
            t2, v2 = _series(M, k + 2)
            mc, _, _ = running_average(t2, v2, window=Wc)
            okc = ~np.isnan(mc)
            miss = float(np.max(np.abs(mc[okc] - target))) / target
            # Predicted O(ripple): the extra half period's mean departs from
            # D·θ_ss by at most (θ_max − θ_min)/2, diluted by P/2 / Wc.
            pred = 0.5 * (v.max() - v.min()) * (0.5 * P / Wc) / target
            _p(
                f"[MAT-20 Wc] k={k} M={M}: W={Wc:g} miss {miss:.3e} "
                f"(predicted ≤ {pred:.3e}; must exceed {BAR:g})"
            )
            assert miss > BAR


def test_b2_window_max_phase_and_pc():
    # Root of θ(t) = θ(t − W) on the falling side: t ∈ (DP, DP + W) puts
    # t in the off phase and t − W in the on phase.
    g = lambda s: float(_theta(np.array([s]))[0] - _theta(np.array([s - W_B2]))[0])
    a, b = D * P + 1.0, D * P + W_B2 - 1.0
    assert g(a) > 0.0 > g(b), (g(a), g(b))
    root = brentq(g, a, b, xtol=1e-12, rtol=4 * np.finfo(float).eps)
    for M in M_SEQ:
        dt = P / M
        t, v = _series(M, 2)
        # restrict to the second period so the max is unique mod P
        means, _, _ = running_average(t, v, window=W_B2)
        sel = t >= P
        i = np.nanargmax(np.where(sel, means, -np.inf))
        t_max = t[i] - P
        off = t_max - root
        _p(
            f"[MAT-20 B2] M={M}: brentq root {root:.6f} s, t_max {t_max:.6f} s, "
            f"offset {off:+.6f} s = {off / dt:+.3f} samples (Δt {dt:.6f} s, "
            f"W/Δt = {W_B2 / dt:.3f})"
        )
        assert abs(off) <= dt
        # (Pc) printed: window one sample short.
        mp, _, _ = running_average(t, v, window=W_B2 - dt)
        ip = np.nanargmax(np.where(sel, mp, -np.inf))
        _p(
            f"[MAT-20 Pc] M={M}: W−Δt → t_max {t[ip] - P:.6f} s "
            f"(shift {(t[ip] - P) - t_max:+.6f} s; predicted ≈ −Δt/2 = {-dt / 2:+.6f})"
        )
