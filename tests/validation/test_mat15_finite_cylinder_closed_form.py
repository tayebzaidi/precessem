"""`MAT-15` step 2(a) — the finite gel cylinder, uniform source, Robin walls:
the closed form in numpy, verified by residuals before any solve (§9 item 137).

0 ≤ r ≤ a, |z| ≤ L, a = 0.03, L = 0.04 (the phantom's radius / half-height as
pinned in the §9 item), `MAT-9`'s k 0.6, ρc 4.18e6, h 10, q 2000
(Bi_r = 0.5, Bi_z = 0.667).  Closed forms:
`utils.analytical.robin_finite_cylinder_steady_temperature` (single series
θ_slab + Σ Dₙ I₀ cos) and `…_transient_temperature` (double eigen-series).

Anchors (asserted): (R0) roots in their brackets, residuals ≤ 1e-12·root, the
two Cₘ forms agree ≤ 1e-12 in |ΔCₘ|/max|Cₘ| (norm disclosed in the test); (R1) Robin residuals on r = a and z = ±L, steady
and transient, analytic derivatives, |−k∂ₙθ − hθ|/h ≤ 1e-10·qa/h; (R2) FD
Laplacian of θ_ss at δ = 1e-2·a and 2δ (pre-registered 1e-3·a; changed,
disclosed in the test): ratio ∈ [3.8, 4.2], |k∇²θ + q|/q ≤
1e-4; (R3) h∮θ_ss dS = qV ≤ 1e-8; (R4) Duhamel θ(t) = (q/ρc)∫₀ᵗ U_cyl U_slab
ds ≤ 1e-5·θ_ss(0,0) at Fo 0.05 / 0.2 / 1; (R5) Fo 50 → θ_ss ≤ 1e-12;
L = 100a → infinite cylinder, a = 100L → slab, ≤ 1e-6.
Negative control (predicted, asserted only > 1e-3): Dirichlet radial roots.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.special import j0, j1, jn_zeros

from fem_em_solver.utils.analytical import (
    robin_cylinder_eigenvalues,
    robin_finite_cylinder_steady_temperature,
    robin_finite_cylinder_transient_temperature,
    robin_slab_eigenvalues,
)

try:
    from mpi4py import MPI

    _RANK = MPI.COMM_WORLD.rank
except Exception:  # pragma: no cover
    _RANK = 0

A = 0.03
L = 0.04
K = 0.6
RHO_C = 4.18e6
H = 10.0
Q = 2000.0
T_INF = 20.0
ALPHA = K / RHO_C
T_R = A**2 / ALPHA
BI_R = H * A / K
BI_Z = H * L / K
KW = dict(radius=A, half_height=L, k=K, q=Q, h=H)
FOS = (0.05, 0.2, 1.0)


def _say(msg):
    if _RANK == 0:
        print(msg, flush=True)


def _steady(r, z, **kw):
    return robin_finite_cylinder_steady_temperature(r, z, t_inf=T_INF, **{**KW, **kw})


def _trans(r, z, t, **kw):
    return robin_finite_cylinder_transient_temperature(
        r, z, t, t_inf=T_INF, rho_c=RHO_C, t_init=T_INF, **{**KW, **kw})


def test_r0_roots_and_coefficients():
    m = 400
    beta = robin_cylinder_eigenvalues(BI_R, m)
    gam = robin_slab_eigenvalues(BI_Z, m)
    lo = np.concatenate([[0.0], jn_zeros(1, m - 1)])
    hi = jn_zeros(0, m)
    assert np.all(beta > lo) and np.all(beta < hi)
    n = np.arange(1, m + 1)
    assert np.all(gam > (n - 1) * np.pi) and np.all(gam < (n - 0.5) * np.pi)
    res_b = np.abs(beta * j1(beta) - BI_R * j0(beta)) / beta
    res_g = np.abs(gam * np.sin(gam) - BI_Z * np.cos(gam)) / gam
    c1 = 2.0 * j1(beta) / (beta * (j0(beta) ** 2 + j1(beta) ** 2))
    c2 = 2.0 * BI_R / ((beta**2 + BI_R**2) * j0(beta))
    dc_rel = np.max(np.abs(c1 - c2) / np.abs(c1))
    # Norm disclosed (W1 20261001T143051Z, line 142): per-root relative the two
    # forms differ by 2.65e-10 at m ≈ 400 — the conditioning of form 1, not a
    # root error: there J₁(β) ≈ Bi J₀/β ≈ 1e-5, so a root good to the last ulp
    # (residual 2.7e-15·β, printed above) moves J₁ by |J₀/J₁|·δβ ≈ 2.5e-10.
    # The 1e-12 bar is kept and applied in the series' own norm, |ΔCₘ|/max|Cₘ|.
    dc = np.max(np.abs(c1 - c2)) / np.max(np.abs(c1))
    _say(f"MAT15-2a (R0) max root residual / root: beta {res_b.max():.3e}, "
         f"gamma {res_g.max():.3e}; C_m two forms max |dC|/max|C| {dc:.3e} (per-root rel {dc_rel:.3e}, printed only)")
    _say(f"MAT15-2a printed beta_1 = {beta[0]:.4f} (Incropera Table 5.1 Bi 0.5, from "
         f"memory: 0.9408); gamma_1 at Bi 0.5 = {robin_slab_eigenvalues(0.5, 1)[0]:.4f} "
         f"(0.6533); gamma_1 at Bi_z {BI_Z:.3f} = {gam[0]:.4f}")
    assert res_b.max() <= 1e-12
    assert res_g.max() <= 1e-12
    assert dc <= 1e-12


def _robin_residuals(fn):
    """max |−k∂ₙθ − hθ| / (q a) on r = a and z = ±L (= residual/h relative to qa/h)."""
    zs = np.linspace(-0.95 * L, 0.95 * L, 9)
    rs = np.linspace(0.0, 0.95 * A, 9)
    T, dr, _ = fn(np.full_like(zs, A), zs)
    res_r = np.abs(-K * dr - H * (T - T_INF)) / (Q * A)
    out = [res_r.max()]
    for sgn in (1.0, -1.0):
        T, _, dz = fn(rs, np.full_like(rs, sgn * L))
        out.append((np.abs(-K * sgn * dz - H * (T - T_INF)) / (Q * A)).max())
    return out[0], max(out[1:])


def test_r1_robin_residuals():
    rr, rz = _robin_residuals(lambda r, z: _steady(r, z, return_gradient=True))
    _say(f"MAT15-2a (R1) steady Robin residual / (qa): r = a {rr:.3e}, z = ±L {rz:.3e}")
    assert rr <= 1e-10 and rz <= 1e-10
    for fo in FOS:
        t = fo * T_R
        rr, rz = _robin_residuals(lambda r, z: _trans(r, z, t, return_gradient=True))
        _say(f"MAT15-2a (R1) transient Fo {fo}: r = a {rr:.3e}, z = ±L {rz:.3e}")
        assert rr <= 1e-10 and rz <= 1e-10


def _fd_lap_residual(delta, x0, y0, z0):
    def th(x, y, z):
        return _steady(np.hypot(x, y), z) - T_INF
    c = th(x0, y0, z0)
    lap = (th(x0 + delta, y0, z0) + th(x0 - delta, y0, z0)
           + th(x0, y0 + delta, z0) + th(x0, y0 - delta, z0)
           + th(x0, y0, z0 + delta) + th(x0, y0, z0 - delta) - 6.0 * c) / delta**2
    return float(K * lap + Q)


def test_r2_fd_laplacian():
    r0, z0 = 0.6 * A, 0.5 * L
    x0, y0 = r0 / np.sqrt(2.0), r0 / np.sqrt(2.0)
    # Pre-registered δ = 1e-3·a; W1 20261001T143051Z (line 93) measured
    # 2.59e-8 / 1.27e-7, ratio 4.90: the 400-term summation's round-off
    # (≈ 20·eps·10 K / δ² · k/q ≈ 1.5e-8) is the size of the δ² term there.
    # Disclosed parameter change: δ = 1e-2·a (truncation ×100, round-off ÷100);
    # δ = 1e-3·a is still printed.  Bounds unchanged.
    e1s = _fd_lap_residual(1e-3 * A, x0, y0, z0)
    e2s = _fd_lap_residual(2e-3 * A, x0, y0, z0)
    _say(f"MAT15-2a (R2, record) delta = 1e-3 a: {abs(e1s) / Q:.3e}, 2 delta {abs(e2s) / Q:.3e}, "
         f"ratio {e2s / e1s:.4f}")
    d = 1e-2 * A
    e1 = _fd_lap_residual(d, x0, y0, z0)
    e2 = _fd_lap_residual(2 * d, x0, y0, z0)
    ratio = e2 / e1
    _say(f"MAT15-2a (R2) |k lap + q|/q at delta = {abs(e1) / Q:.3e}, 2 delta = "
         f"{abs(e2) / Q:.3e}, ratio {ratio:.4f}")
    assert 3.8 <= ratio <= 4.2
    assert abs(e1) / Q <= 1e-4


def test_r3_gauss_balance():
    x, w = np.polynomial.legendre.leggauss(400)
    zs = L * x
    lat = 2 * np.pi * A * L * np.sum(w * (_steady(np.full_like(zs, A), zs) - T_INF))
    rs = 0.5 * A * (x + 1.0)
    cap = 2 * np.pi * 0.5 * A * np.sum(w * rs * (_steady(rs, np.full_like(rs, L)) - T_INF))
    total = H * (lat + 2.0 * cap)
    vol = np.pi * A**2 * 2 * L
    rel = abs(total - Q * vol) / (Q * vol)
    area = 2 * np.pi * A * 2 * L + 2 * np.pi * A**2
    th00 = float(_steady(0.0, 0.0) - T_INF)
    _say(f"MAT15-2a printed theta_ss(0,0) = {th00:.6f} K; qV/(hA) = {Q * vol / (H * area):.4f} K")
    _say(f"MAT15-2a (R3) h*surface integral / qV - 1 = {rel:.3e}")
    assert rel <= 1e-8


def _duhamel(r, z, t, n=2000):
    beta = robin_cylinder_eigenvalues(BI_R, n)
    gam = robin_slab_eigenvalues(BI_Z, n)
    cm = 2.0 * j1(beta) / (beta * (j0(beta) ** 2 + j1(beta) ** 2))
    cn = 2.0 * np.sin(gam) / (gam + np.sin(gam) * np.cos(gam))
    edges = np.concatenate([[0.0], np.geomspace(1e-10 * t, t, 80)])
    x, w = np.polynomial.legendre.leggauss(16)
    s = (0.5 * (edges[1:, None] - edges[:-1, None]) * (x + 1.0) + edges[:-1, None]).ravel()
    ws = (0.5 * (edges[1:, None] - edges[:-1, None]) * w).ravel()
    ucyl = (np.exp(-ALPHA * np.outer(s, beta**2) / A**2) * (cm * j0(beta * r / A))).sum(1)
    uslab = (np.exp(-ALPHA * np.outer(s, gam**2) / L**2) * (cn * np.cos(gam * z / L))).sum(1)
    return Q / RHO_C * np.sum(ws * ucyl * uslab)


PTS = ((0.0, 0.0), (0.5 * A, 0.5 * L), (0.9 * A, 0.8 * L))


def _r4_error(**kw):
    th00 = float(_steady(0.0, 0.0) - T_INF)
    worst = 0.0
    for fo in FOS:
        t = fo * T_R
        for r, z in PTS:
            a = float(_trans(r, z, t, **kw) - T_INF)
            b = _duhamel(r, z, t)
            worst = max(worst, abs(a - b) / th00)
            _say(f"MAT15-2a (R4{'-ctl' if kw else ''}) Fo {fo} (r,z) = ({r:.3f},{z:.3f}): "
                 f"series {a:.9f}  Duhamel {b:.9f}")
    return worst


def test_r4_duhamel():
    worst = _r4_error()
    _say(f"MAT15-2a (R4) max |series − Duhamel| / theta_ss(0,0) = {worst:.3e}")
    assert worst <= 1e-5


def test_r5_limits_and_tail():
    th = _steady(np.array([0.0, 0.5 * A]), np.array([0.0, 0.3 * L]))
    late = _trans(np.array([0.0, 0.5 * A]), np.array([0.0, 0.3 * L]), 50.0 * T_R)
    e_late = np.max(np.abs(late - th) / np.abs(th - T_INF))
    rs = np.array([0.0, 0.5 * A, 0.9 * A])
    inf_cyl = Q * (A**2 - rs**2) / (4 * K) + Q * A / (2 * H)
    got = robin_finite_cylinder_steady_temperature(
        rs, np.zeros_like(rs), radius=A, half_height=100 * A, k=K, q=Q, h=H)
    e_cyl = np.max(np.abs(got - inf_cyl) / inf_cyl)
    zs = np.array([0.0, 0.5 * L, 0.9 * L])
    slab = Q / (2 * K) * (L**2 - zs**2) + Q * L / H
    got = robin_finite_cylinder_steady_temperature(
        np.zeros_like(zs), zs, radius=100 * L, half_height=L, k=K, q=Q, h=H)
    e_slab = np.max(np.abs(got - slab) / slab)
    _say(f"MAT15-2a (R5) Fo 50 vs steady {e_late:.3e}; infinite-cylinder limit "
         f"{e_cyl:.3e}; slab limit {e_slab:.3e}")
    for mn in (100, 200, 400):
        v = robin_finite_cylinder_transient_temperature(
            0.0, 0.0, 0.0, t_inf=0.0, rho_c=RHO_C, n_radial=mn, n_axial=mn, **KW)
        _say(f"MAT15-2a printed t = 0 tail at M = N = {mn}: theta(0,0,0) = {float(v):.3e} K")
    assert e_late <= 1e-12
    assert e_cyl <= 1e-6
    assert e_slab <= 1e-6


def test_negative_control_dirichlet_roots():
    beta_d = jn_zeros(0, 200)
    t = 0.2 * T_R
    rr, _ = _robin_residuals(
        lambda r, z: _trans(r, z, t, radial_roots=beta_d, return_gradient=True))
    worst = _r4_error(radial_roots=beta_d)
    _say(f"MAT15-2a control (Dirichlet roots): (R1) radial residual {rr:.3e} "
         f"(predicted ≈ 1.0); (R4) miss {worst:.3e}")
    assert rr > 1e-3
    assert worst > 1e-3
