"""`MAT-21` step 11, W1 — P36's transverse quadratic source under `MAT-17`'s
duty cycle re-timed to the box (step 2's schedule, unchanged), closed form in
numpy (§9 item 524, §10 chain P38; brief `docs/planning/chunks/MAT-21.md`,
review ruling 2026-10-18 03:00).  Numpy only (real build, no dolfinx).

θ_duty(t) = Σ_{t > t_m} (−1)^m θ_q(t − t_m), θ_q =
``robin_box_quadratic_transient_temperature``.  The time factor is linear in
the coefficient tensor, so the ten runs collapse onto one tensor:
F = Σ_m (−1)^m (−expm1(−αΛ(t − t_m))), W = T·((Q/(KΛ))·F), evaluated once.
The three private ``analytical.py`` helpers are imported (no ``src/`` edit;
ruling 2026-10-18 03:00).

Asserted: (D0) at ``switch_times=(0.0,)`` bitwise equal to the public
transient at t* on P36's (Q2) point set, at ``SERIES_TERMS`` and
``QUAD_TERMS``; (D1) one-tensor sum vs Σ(−1)^m of ten public calls ≤
1e-12·θ_ref; (D2) vs the Duhamel quadrature with χ ≤ ``R_BAR``·θ_ref.
Control (D-ctl), predicted (only the miss > ``R_BAR`` asserted): m = 9
dropped at t*.  Printed: the duty/continuous ratio table beside the drafting
values and P26's ratios.
"""

from __future__ import annotations

import numpy as np

from fem_em_solver.utils.analytical import (
    _robin_box_eval,
    _robin_box_lambda,
    _robin_box_quadratic_tensor,
    robin_box_quadratic_transient_temperature,
)
from tests.validation.test_mat21_box_closed_form import (
    A_BOX,
    ALPHA,
    B_BOX,
    C_BOX,
    HALF,
    H,
    K,
    Q,
    R2_PTS,
    R4_PTS,
    RHO_C,
    SERIES_TERMS,
    T_C,
    _say,
)
from tests.validation.test_mat21_box_duty_closed_form import (
    C1_TIMES,
    HALF_T,
    R_BAR,
    T_STAR,
    T_SW,
    _on_lags,
    _quad_nodes,
)
from tests.validation.test_mat21_box_quadratic_closed_form import (
    QUAD_TERMS,
    WEIGHTS,
    _uprod_quad,
    _y_reduced,
)

D1_BAR = 1e-12  # step 0's R0 bar (the 10-16 ruling's class)
THETA_REF_DRAFT = 8.805072  # P36's (Q6) reference


def _public(x, y, z, t, n_terms=QUAD_TERMS):
    return robin_box_quadratic_transient_temperature(
        x, y, z, t, half_widths=HALF, k=K, q=Q, h=H, rho_c=RHO_C, weights=WEIGHTS,
        t_inf=0.0, n_terms=n_terms)


def _duty_tensor(t, switch_times, n_terms, drop):
    axes, T = _robin_box_quadratic_tensor(HALF, K, H, n_terms, WEIGHTS)
    lam = _robin_box_lambda(axes)
    F = None
    for m, tm in enumerate(switch_times):
        if t > tm and m not in drop:
            term = (-1.0) ** m * (-np.expm1(-(ALPHA * lam * float(t - tm))))
            F = term if F is None else F + term
    if F is None:
        F = np.zeros_like(lam)
    W = T * ((Q / (K * lam)) * F)
    return W, axes


def _quad_duty(x, y, z, t, switch_times=T_SW, n_terms=QUAD_TERMS, drop=()):
    """Σ_{t > t_m, m ∉ drop} (−1)^m θ_q(t − t_m) through one coefficient tensor."""
    W, axes = _duty_tensor(t, switch_times, n_terms, drop)
    return _robin_box_eval(W, x, y, z, axes, 0.0, False)


def _quad_duty_mean(t, switch_times=T_SW, n_terms=QUAD_TERMS, drop=()):
    """Exact volume mean of the truncated duty series: W contracted with sin γ/γ."""
    W, axes = _duty_tensor(t, switch_times, n_terms, drop)
    sx, sy, sz = (np.sin(g) / g for (_, g, _) in axes)
    return float(np.einsum("lmn,l,m,n->", W, sx, sy, sz))


def _quad_cont_mean(t, n_terms=QUAD_TERMS):
    return _quad_duty_mean(t, switch_times=(0.0,), n_terms=n_terms)


def _theta_ref():
    return max(abs(_y_reduced(*p, 40, 40)) for p in R2_PTS)


def _duhamel_quad_duty(x, y, z, t):
    tot = 0.0
    for lo, hi in _on_lags(t):
        s, ws = _quad_nodes(lo, hi)
        tot += np.sum(ws * _uprod_quad(x, y, z, s))
    return Q / RHO_C * tot


def _q2_points():
    pts = np.array(R2_PTS)
    g = np.meshgrid(np.linspace(-0.9, 0.9, 7) * A_BOX, np.linspace(-0.9, 0.9, 7) * B_BOX,
                    np.linspace(-0.9, 0.9, 3) * C_BOX, indexing="ij")
    return (np.concatenate([pts[:, 0], g[0].ravel()]), np.concatenate([pts[:, 1], g[1].ravel()]),
            np.concatenate([pts[:, 2], g[2].ravel()]))


# ---------------------------------------------------------------- (D0)
def test_d0_single_switch_bitwise_public():
    x, y, z = _q2_points()
    for nt in (SERIES_TERMS, QUAD_TERMS):
        a = _quad_duty(x, y, z, T_STAR, switch_times=(0.0,), n_terms=nt)
        b = _public(x, y, z, T_STAR, n_terms=nt)
        eq = bool(np.array_equal(a, b))
        _say(f"MAT21-11 (D0) switch_times=(0.0,) vs public transient at t* = {T_STAR:.4f} s, "
             f"n_terms {nt}: array_equal {eq} over {x.size} points (max |d| "
             f"{np.max(np.abs(a - b)):.3e} K)")
        assert eq


# ---------------------------------------------------------------- (D1)
def test_d1_one_tensor_vs_ten_public_calls():
    th_ref = _theta_ref()
    worst = 0.0
    for t in C1_TIMES:
        for p in R4_PTS:
            a = float(_quad_duty(*p, t))
            b = 0.0
            for m, tm in enumerate(T_SW):
                if t > tm:
                    b = b + (-1.0) ** m * float(_public(*p, t - tm))
            worst = max(worst, abs(a - b) / th_ref)
            _say(f"MAT21-11 (D1) t = {t:.4f} s {p}: one-tensor {a:.12f}  ten calls {b:.12f}")
    _say(f"MAT21-11 (D1) theta_ref {th_ref:.6f} K (P36 {THETA_REF_DRAFT}); max |d|/theta_ref "
         f"{worst:.3e} (bar {D1_BAR:.0e}, predicted <= 1e-15)")
    assert worst <= D1_BAR


# ---------------------------------------------------------------- (D2)
def test_d2_duhamel_with_chi():
    th_ref = _theta_ref()
    worst = 0.0
    for t in C1_TIMES:
        for p in R4_PTS:
            a = float(_quad_duty(*p, t, n_terms=SERIES_TERMS))
            b = float(_duhamel_quad_duty(*p, t))
            worst = max(worst, abs(a - b) / th_ref)
            _say(f"MAT21-11 (D2) t = {t:.4f} s {p}: ten-run sum {a:.12f}  Duhamel(chi) {b:.12f}  "
                 f"|d|/theta_ref {abs(a - b) / th_ref:.3e}")
    _say(f"MAT21-11 (D2) max |sum - Duhamel|/theta_ref {worst:.3e} (bar {R_BAR:.0e}, "
         f"predicted <= 1e-13)")
    assert worst <= R_BAR


# ---------------------------------------------------------------- control (D-ctl)
def test_control_d_ctl_m9_dropped_misses():
    th_ref = _theta_ref()
    misses = []
    for p in R4_PTS:
        a = float(_quad_duty(*p, T_STAR, n_terms=SERIES_TERMS, drop=(9,)))
        b = float(_duhamel_quad_duty(*p, T_STAR))
        misses.append(abs(a - b) / th_ref)
    hq = float(_public(0.0, 0.0, 0.0, HALF_T, n_terms=SERIES_TERMS))
    _say("MAT21-11 control (D-ctl) PREDICTED: m = 9 dropped at t*: miss/theta_ref per point "
         + " / ".join(f"{v:.4e}" for v in misses)
         + f" (predicted 8.0e-5 / 1.16e-2 / 2.92e-2; theta_q(HALF_T) at the centre {hq:.4e} K, "
         f"predicted 7.01e-4; asserted only > {R_BAR:.0e}; min/bar {min(misses) / R_BAR:.2e})")
    assert min(misses) > R_BAR


# ---------------------------------------------------------------- printed records
def test_printed_duty_ratio_table():
    rows = (("centre", (0.0, 0.0, 0.0), 0.063466, 0.034567, 0.54466),
            ("(0, B, 0)", (0.0, B_BOX, 0.0), 2.647811, 1.285771, 0.48560),
            ("(0, 0.9B, 0)", (0.0, 0.9 * B_BOX, 0.0), 2.892661, 1.430931, 0.49468),
            ("(A, 0, 0)", (A_BOX, 0.0, 0.0), 0.048327, 0.026173, 0.54159),
            ("(0, 0, C)", (0.0, 0.0, C_BOX), 0.100217, 0.048678, 0.48572))
    for lab, p, dc, dd, dr in rows:
        c = float(_quad_duty(*p, T_STAR, switch_times=(0.0,)))
        d = float(_quad_duty(*p, T_STAR))
        _say(f"MAT21-11 record {lab}: theta_cont(t*) {c:.6f} K (drafting {dc:.6f}), theta_duty(t*) "
             f"{d:.6f} K (drafting {dd:.6f}), ratio {d / c:.5f} (drafting {dr:.5f})")
    mc = _quad_cont_mean(T_STAR)
    md = _quad_duty_mean(T_STAR)
    _say(f"MAT21-11 record volume mean: theta_cont(t*) {mc:.6f} K (drafting 1.209979), theta_duty(t*) "
         f"{md:.6f} K (drafting 0.598518), ratio {md / mc:.5f} (drafting 0.49465)")
    _say("MAT21-11 record P26 uniform ratios (prediction, printed only): volume mean 0.49576 "
         "(0.61419 / 1.23889), centre 0.49902 (0.666540 / 1.335705)")
    _say(f"MAT21-11 record t* = 0.2 t_C = {0.2 * T_C:.4f} s; switches {len(T_SW)}")
