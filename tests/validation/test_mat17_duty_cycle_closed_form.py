"""`MAT-17` step 1 — the duty-cycled Robin gel ball in closed form (§9 item 189).

Numpy only.  Fixture: `MAT-9`'s ball (R 0.1, k 0.6, ρc 4.18e6, q 2000, h 10,
Bi 5/3, T∞ = T(0) = 20), source on while 0 ≤ (t mod P) < DP.
Case A: P = 1200 s, D = 0.5.  Case B: P = 0.2 t_R, D = 0.25.

Anchor: `utils.analytical.robin_sphere_duty_cycle_temperature` (modal form,
derivation in its docstring).  Independent routes: the Fourier series on
`robin_sphere_periodic_harmonic` (no λₙ, no Bₙ) and the Duhamel superposition
of `robin_sphere_transient_temperature` over the switch instants.

Asserted: (H0) complex-w steady state at real w = `MAT-11`; (H1) Robin
residuals; (P1) cycle mean = D·θ_ss; (P2) modal periodic = Fourier at
M = 2¹⁶; (P3) from rest = Duhamel, cycles 0–5; (P4) duty 1 = step response,
duty ∉ (0, 1] raises; (P5) from rest at j* = ⌈30τ₁/P⌉ = periodic.
Negative controls (rule (e), predicted, asserted only as exceeding the
anchor's band): (c1) Dirichlet roots; (c2) geometric sum dropped; (c3) the
e^{−iωt} convention.  All errors are on θ = T − T∞.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from mpi4py import MPI

from fem_em_solver.utils.analytical import (
    _duty_cycle_modal,
    _robin_sphere_complex_w,
    _robin_sphere_complex_w_dR,
    robin_sphere_duty_cycle_temperature,
    robin_sphere_eigenvalues,
    robin_sphere_perfusion_temperature,
    robin_sphere_periodic_harmonic,
    robin_sphere_steady_temperature,
    robin_sphere_transient_coefficients,
    robin_sphere_transient_temperature,
)

R, K, RHO_C, Q, H, T_INF = 0.1, 0.6, 4.18e6, 2000.0, 10.0, 20.0
T_R = R**2 * RHO_C / K
CASES = {"A": (1200.0, 0.5), "B": (0.2 * T_R, 0.25)}
R_EVAL = np.array([0.0, R / 2, R])
FRAC = np.linspace(0.05, 0.95, 8)
PHYS = dict(radius=R, k=K, rho_c=RHO_C, q=Q, h=H)
TH_SS0 = float(robin_sphere_steady_temperature(0.0, radius=R, k=K, q=Q, h=H))
LAM1 = float(robin_sphere_eigenvalues(H * R / K, 1)[0])
TAU1 = T_R / LAM1**2


def _p(*a):
    if MPI.COMM_WORLD.rank == 0:
        print(*a, flush=True)


def _cycle_times(P, D, j=0):
    """Eight times per phase in cycle j, fractions 0.05 … 0.95 of the phase."""
    on = j * P + FRAC * D * P
    off = j * P + D * P + FRAC * (1.0 - D) * P if D < 1.0 else np.empty(0)
    return np.concatenate([on, off])


def _grid(P, D, cycles=(0,)):
    t = np.concatenate([_cycle_times(P, D, j) for j in cycles])
    rr, tt = np.meshgrid(R_EVAL, t, indexing="ij")
    return rr, tt


def _modal(rr, tt, P, D, **kw):
    return _duty_cycle_modal(rr, tt, period=P, duty=D, **PHYS, **kw)


def _fourier(r_vals, t, P, D, M, sign=1.0, block=8192):
    """θ_p = D θ_ss + 2 Re Σ_{m=1}^{M} ûₘ Θ(r; sign·ωₘ) e^{iωₘt}; shape (nr, nt)."""
    t = np.asarray(t, dtype=float)
    out = np.empty((len(r_vals), len(t)))
    for i, r in enumerate(r_vals):
        acc = np.zeros(len(t), dtype=complex)
        for m0 in range(1, M + 1, block):
            m = np.arange(m0, min(m0 + block, M + 1), dtype=float)
            om = 2.0 * np.pi * m / P
            u = (1.0 - np.exp(-1j * om * D * P)) / (1j * om * P)
            th = robin_sphere_periodic_harmonic(r, sign * om, **PHYS)
            acc += np.exp(1j * np.outer(t, om)) @ (u * th)
        ss = float(robin_sphere_steady_temperature(r, radius=R, k=K, q=Q, h=H))
        out[i] = D * ss + 2.0 * acc.real
    return out


def _duhamel(rr, tt, P, D):
    """Σⱼ (−1)ʲ θ_step(t − tⱼ) over switch instants on at jP, off at jP + DP."""
    tmax = float(np.max(tt))
    theta = np.zeros_like(tt)
    j = 0
    while j * P <= tmax:
        for ts, sgn in ((j * P, 1.0), (j * P + D * P, -1.0)):
            if D == 1.0 and sgn < 0:
                continue
            mask = tt > ts
            step = (
                robin_sphere_transient_temperature(
                    rr,
                    np.maximum(tt - ts, 0.0),
                    radius=R,
                    k=K,
                    rho_c=RHO_C,
                    q=Q,
                    h=H,
                    t_inf=T_INF,
                    t_init=T_INF,
                )
                - T_INF
            )
            theta += sgn * np.where(mask, step, 0.0)
        j += 1
    return theta


def _rel(a, b):
    return float(np.max(np.abs(a - b))) / TH_SS0


# ---------------------------------------------------------------------------


def test_h0_complex_w_matches_mat11_at_real_w():
    for w in (240.0, 2400.0, 2.4e7):
        z = _robin_sphere_complex_w(R_EVAL, w, radius=R, k=K, q=Q, h=H)
        ref = robin_sphere_perfusion_temperature(R_EVAL, radius=R, k=K, q=Q, w=w, h=H)
        err = float(np.max(np.abs(z - ref)) / np.max(np.abs(ref)))
        _p(f"[MAT-17 H0] w={w:g}: max rel |Θ(w) − MAT-11| = {err:.3e} (band 1e-13)")
        assert err <= 1e-13


@pytest.mark.parametrize("case", ["A", "B"])
def test_h1_robin_residuals(case):
    P, D = CASES[case]
    for m in (1, 10, 1000):
        w = 1j * 2.0 * np.pi * m / P * RHO_C
        th_R = complex(_robin_sphere_complex_w(R, w, radius=R, k=K, q=Q, h=H))
        dth = complex(_robin_sphere_complex_w_dR(radius=R, k=K, q=Q, h=H, w=w))
        res = abs(-K * dth - H * th_R) / (H * abs(th_R))
        _p(f"[MAT-17 H1 case {case}] m={m}: |−kΘ′−hΘ|/(h|Θ|) = {res:.3e} (band 1e-12)")
        assert res <= 1e-12
    t = _cycle_times(P, D)
    th = _modal(np.full_like(t, R), t, P, D, periodic=True)
    dth = _modal(np.full_like(t, R), t, P, D, periodic=True, want="dR")
    res = float(np.max(np.abs(-K * dth - H * th))) / (Q * R / H)
    _p(f"[MAT-17 H1 case {case}] periodic modal Robin residual / (qR/h) = {res:.3e} (band 1e-10)")
    assert res <= 1e-10


@pytest.mark.parametrize("case", ["A", "B"])
def test_p1_cycle_mean(case):
    P, D = CASES[case]
    lam = robin_sphere_eigenvalues(H * R / K, 400)
    _, B = robin_sphere_transient_coefficients(lam, radius=R, k=K, q=Q, h=H)
    s = lam**2 * K / (RHO_C * R**2)
    G = np.exp(-s * (1 - D) * P) * (-np.expm1(-s * D * P)) / (-np.expm1(-s * P))
    g_end = 1.0 - (1.0 - G) * np.exp(-s * D * P)
    # ∫ on = DP θ_ss − Σ BX (1−G)(1−e^{−sDP})/s ; ∫ off = Σ BX g_end (1−e^{−s(1−D)P})/s
    coef = (g_end * (-np.expm1(-s * (1 - D) * P)) - (1 - G) * (-np.expm1(-s * D * P))) / s
    x = lam * R_EVAL[:, None] / R
    X = np.where(x == 0, 1.0, np.sin(x) / np.where(x == 0, 1.0, x))
    ss = robin_sphere_steady_temperature(R_EVAL, radius=R, k=K, q=Q, h=H)
    mean = D * ss + np.sum(B * X * coef, axis=-1) / P
    err = _rel(mean, D * ss)
    _p(f"[MAT-17 P1 case {case}] cycle mean {mean} vs D·θ_ss {D * ss}: rel {err:.3e} (band 1e-12)")
    assert err <= 1e-12


@pytest.mark.parametrize("case", ["A", "B"])
def test_p2_modal_vs_fourier_and_controls(case):
    P, D = CASES[case]
    t = _cycle_times(P, D)
    rr, tt = np.meshgrid(R_EVAL, t, indexing="ij")
    modal = _modal(rr, tt, P, D, periodic=True)
    ref = _fourier(R_EVAL, t, P, D, 2**16)
    err = _rel(modal, ref)
    _p(
        f"[MAT-17 P2 case {case}] modal vs Fourier M=2^16: {err:.3e}·θ_ss(0) "
        f"(band 1e-9, predicted ≲ 5e-11)"
    )
    # Fourier ladder (printed only)
    sw = np.array([0.0, D * P])
    ref_sw = _fourier(R_EVAL, sw, P, D, 2**16)
    for e in range(10, 17):
        f_e = _fourier(R_EVAL, t, P, D, 2**e)
        f_sw = _fourier(R_EVAL, sw, P, D, 2**e)
        _p(
            f"[MAT-17 ladder case {case}] M=2^{e}: eval-set {_rel(f_e, modal):.3e}, "
            f"switch instants vs M=2^16 {_rel(f_sw, ref_sw):.3e}"
        )
    # (c1) Dirichlet roots, B re-projected
    lam_d = np.pi * np.arange(1, 401)
    c1 = _modal(rr, tt, P, D, periodic=True, lam=lam_d)
    e1 = _rel(c1, ref)
    e1c = float(np.max(np.abs(c1[0] - ref[0]))) / TH_SS0
    _p(
        f"[MAT-17 c1 case {case}] Dirichlet roots: {e1:.3e} (centre {e1c:.3e}; "
        f"predicted ≈ 0.27 centre in A; asserted > 1e-3)"
    )
    # (c3) e^{−iωt} convention
    c3 = _fourier(R_EVAL, t, P, D, 2**16, sign=-1.0)
    e3 = _rel(c3, modal)
    _p(
        f"[MAT-17 c3 case {case}] e^(−iωt) convention: {e3:.3e} "
        f"(predicted ≈ {'1e-2' if case == 'A' else '1e-1'}; asserted > 1e-6)"
    )
    assert err <= 1e-9
    assert e1 > 1e-3
    assert e3 > 1e-6


@pytest.mark.parametrize("case", ["A", "B"])
def test_p3_from_rest_vs_duhamel(case):
    P, D = CASES[case]
    rr, tt = _grid(P, D, cycles=range(6))
    modal = _modal(rr, tt, P, D)
    duh = _duhamel(rr, tt, P, D)
    err = _rel(modal, duh)
    _p(f"[MAT-17 P3 case {case}] from rest vs Duhamel, cycles 0–5: {err:.3e}·θ_ss(0) (band 1e-10)")
    assert err <= 1e-10


def test_p4_duty_one_and_validation():
    P = CASES["A"][0]
    t = np.concatenate([j * P + FRAC * P for j in range(6)])
    rr, tt = np.meshgrid(R_EVAL, t, indexing="ij")
    a = robin_sphere_duty_cycle_temperature(rr, tt, **PHYS, t_inf=T_INF, period=P, duty=1.0) - T_INF
    b = robin_sphere_transient_temperature(rr, tt, **PHYS, t_inf=T_INF, t_init=T_INF) - T_INF
    err = _rel(a, b)
    _p(f"[MAT-17 P4] duty 1 vs step response: {err:.3e}·θ_ss(0) (band 1e-13)")
    assert err <= 1e-13
    for bad in (0.0, -0.1, 1.5):
        with pytest.raises(ValueError):
            robin_sphere_duty_cycle_temperature(0.0, 1.0, **PHYS, t_inf=T_INF, period=P, duty=bad)


@pytest.mark.parametrize("case", ["A", "B"])
def test_p5_from_rest_reaches_periodic_and_records(case):
    P, D = CASES[case]
    jstar = math.ceil(30.0 * TAU1 / P)
    rr, tt = _grid(P, D, cycles=(jstar,))
    rest = _modal(rr, tt, P, D)
    per = _modal(rr, tt, P, D, periodic=True)
    err = _rel(rest, per)
    c2 = _modal(rr, tt, P, D, geometric=False)
    e2 = _rel(c2, per)
    ss0 = TH_SS0
    peak = float(_modal(0.0, D * P * (1 - 1e-12), P, D, periodic=True))
    trough = float(_modal(0.0, P * (1 - 1e-12), P, D, periodic=True))
    _p(
        f"[MAT-17 records] θ_ss(0) = {ss0:.6f} K (pred 12.2222), λ₁ = {LAM1:.6f} "
        f"(pred ≈ 1.91), τ₁ = {TAU1:.6e} s (pred ≈ 1.91e4)"
    )
    _p(
        f"[MAT-17 P5 case {case}] j* = {jstar}: from rest vs periodic {err:.3e}·θ_ss(0) (band 1e-12)"
    )
    _p(
        f"[MAT-17 c2 case {case}] geometric sum dropped: {e2:.3e} (predicted ≈ 0.49; asserted > 1e-3)"
    )
    _p(
        f"[MAT-17 case {case}] centre periodic: mean {D * ss0:.6f} K, peak {peak:.6f} K, "
        f"trough {trough:.6f} K, ripple p-p {peak - trough:.6f} K, "
        f"peak/(D θ_ss) − 1 = {100 * (peak / (D * ss0) - 1):+.3f} %"
        + (" (pred ripple 0.1435 K, +1.2 %)" if case == "A" else "")
    )
    if case == "A":
        rr0, tt0 = _grid(P, D)
        jc = None
        for j in range(0, 400):
            d = _modal(rr0, tt0 + j * P, P, D) - _modal(rr0, tt0 + j * P, P, D, periodic=True)
            if float(np.max(np.abs(d))) <= 1e-3 * D * ss0:
                jc = j
                break
        _p(f"[MAT-17 case A] cycles to within 1e-3·D·θ_ss(0) of periodic: {jc} (pred ≈ 110)")
    assert err <= 1e-12
    assert e2 > 1e-3
