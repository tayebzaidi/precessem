"""IMP-6 step 0 (§9 item 403): PEC oblate spheroid in a uniform field, axial
and in-plane drives — the quasi-static closed form, numpy only.

Asserted identities (derivations in ``utils/analytical.py``, IMP-6 block):
(L) the two L_z forms; (Sum) 2L_x + L_z = 1 with L_x computed independently;
(Int) each L_i against the ellipsoid integral by quadrature; (Br) closed vs
series branches of h, g; (Bd) φ = 0 on the body; (Tg) tangential E = 0 on the
body; (Rim) E_x(±a,0,0)·L_x = 1 and E(0,±a,0) = 0; (Pole) E_z(0,0,±c)·L_z = 1;
(F) field vs FD gradient (Richardson / ratio form); (Lap) FD Laplacian
(ratio form); (Far) l = 1 projections equal p_i = E₀a²c/(3L_i); (Sph) the
e → 0 limit; (Disc) the ξ₀ → 0 limit.  Negative controls (rule (e), only the
separation asserted): (c1) the prolate formula in place of the oblate one;
(c2) the regular (growing) solution in place of the decaying one.
"""

import numpy as np
import pytest
from scipy.integrate import quad

from fem_em_solver.utils.analytical import (
    oblate_spheroid_d_functions,
    oblate_spheroid_depolarisation,
    oblate_spheroidal_xi,
    pec_oblate_spheroid_field,
    pec_oblate_spheroid_potential,
    prolate_spheroid_depolarisation,
)

A = 1.0  # equatorial semi-axis (semi_major)
E0 = 1.0
ASPECTS = [1.01, 2.0, 5.0, 20.0, 100.0]
DRIVES = ["axial", "inplane"]


def _ecc(aspect):
    return np.sqrt(1.0 - 1.0 / aspect**2)


def _geom(aspect):
    c = A / aspect
    f = np.sqrt((A - c) * (A + c))
    return c, f, c / f


def _surface_points(a, c, n=200):
    """Fibonacci points on (x²+y²)/a² + z²/c² = 1."""
    i = np.arange(n) + 0.5
    mu = 1.0 - 2.0 * i / n
    ang = np.pi * (1.0 + np.sqrt(5.0)) * i
    s = np.sqrt(1.0 - mu**2)
    return np.column_stack([a * s * np.cos(ang), a * s * np.sin(ang), c * mu])


def _oblate_to_cart(xi, eta, ang, f):
    rho = f * np.sqrt((xi**2 + 1.0) * (1.0 - eta**2))
    return np.column_stack([rho * np.cos(ang), rho * np.sin(ang), f * xi * eta])


def _axis(drive):
    return 2 if drive == "axial" else 0


def _l_of(aspect, drive):
    lx, lz = oblate_spheroid_depolarisation(_ecc(aspect))
    return lz if drive == "axial" else lx


def _l1_projection(phi_fn, radius, drive, n=64):
    """l = 1 projection of φ on the sphere of radius R onto x_i/r:
    axial R²(3/2)∫φ μ dμ (φ-angle 0.3); in-plane R²(3/4π)∫φ sinθ cosφ dΩ
    (n-node Gauss–Legendre in μ × n-point trapezoid in φ)."""
    mu, w = np.polynomial.legendre.leggauss(n)
    s = np.sqrt(1.0 - mu**2)
    if drive == "axial":
        pts = radius * np.column_stack([s * np.cos(0.3), s * np.sin(0.3), mu])
        return radius**2 * 1.5 * np.sum(w * mu * phi_fn(pts))
    ang = 2.0 * np.pi * np.arange(n) / n
    M, P = np.meshgrid(mu, ang, indexing="ij")
    W = np.repeat(w[:, None], n, axis=1) * (2.0 * np.pi / n)
    S = np.sqrt(1.0 - M**2)
    pts = radius * np.column_stack(
        [(S * np.cos(P)).ravel(), (S * np.sin(P)).ravel(), M.ravel()])
    vals = phi_fn(pts).reshape(M.shape)
    return radius**2 * 3.0 / (4.0 * np.pi) * np.sum(W * vals * S * np.cos(P))


# ------------------------------------------------------------- (L) / (Sum)
@pytest.mark.parametrize("aspect", ASPECTS)
def test_depolarisation_forms_and_sum(aspect):
    e = _ecc(aspect)
    lx, lz = oblate_spheroid_depolarisation(e)
    lz_asin = (1.0 / e**2) * (1.0 - np.sqrt(1.0 - e**2) * np.arcsin(e) / e)
    d_l = abs(lz - lz_asin)
    d_sum = abs(2.0 * lx + lz - 1.0)
    print(f"[IMP-6 L] aspect={aspect}: L_x={lx:.6f} L_z={lz:.6f} "
          f"1/L_x={1/lx:.4f} 1/L_z={1/lz:.4f} rim (1/L_x)^2={1/lx**2:.2f}; "
          f"|L_z - asin form|={d_l:.3e} |2L_x+L_z-1|={d_sum:.3e}")
    assert d_l <= 1e-13
    assert d_sum <= 1e-13


def test_printed_values_disc_side():
    for aspect, lz_p, lx_p in [(2.0, 0.527200, 0.236400), (5.0, 0.750484, 0.124758)]:
        lx, lz = oblate_spheroid_depolarisation(_ecc(aspect))
        print(f"[IMP-6 P] aspect={aspect}: L_z={lz:.6f} (predicted {lz_p}) "
              f"L_x={lx:.6f} (predicted {lx_p})")
    lx, _ = oblate_spheroid_depolarisation(_ecc(100.0))
    c = A / 100.0
    print(f"[IMP-6 P] aspect=100: L_x={lx:.6f} vs pi c/(4a)={np.pi*c/(4*A):.6f}; "
          f"p_x/a^3={c/(3*lx):.6f} vs 4/(3pi)={4/(3*np.pi):.6f} "
          f"(predicted 0.007755 / 0.429822)")


# -------------------------------------------------------------------- (Int)
@pytest.mark.parametrize("aspect", [2.0, 5.0, 20.0])
def test_ellipsoid_integral(aspect):
    c = A / aspect
    lx, lz = oblate_spheroid_depolarisation(_ecc(aspect))
    pre = 0.5 * A * A * c

    def integrand(s, ai2):
        return 1.0 / ((ai2 + s) * (A * A + s) * np.sqrt(c * c + s))

    out = {}
    for name, ai2, ref in (("L_x", A * A, lx), ("L_z", c * c, lz)):
        val, err = quad(integrand, 0.0, np.inf, args=(ai2,), epsabs=0.0,
                        epsrel=1e-13, limit=500)
        out[name] = (pre * val, abs(pre * val - ref) / ref, pre * err)
    print(f"[IMP-6 Int] aspect={aspect}: " + "; ".join(
        f"{k} quad={v[0]:.15f} rel={v[1]:.3e} (quad err est {v[2]:.1e})"
        for k, v in out.items()))
    for v in out.values():
        assert v[1] <= 1e-10


# --------------------------------------------------------------------- (Br)
def test_branches():
    xi = np.linspace(5.0, 20.0, 31)
    hc, dhc, gc, dgc = oblate_spheroid_d_functions(xi, branch="closed")
    hs, dhs, gs, dgs = oblate_spheroid_d_functions(xi, branch="series")
    rh = np.max(np.abs(hc - hs) / np.abs(hs))
    rg = np.max(np.abs(gc - gs) / np.abs(gs))
    # Derivatives against a central difference (independent of the formulas).
    xx = np.geomspace(0.05, 50.0, 40)
    dd = 1e-5 * xx
    hp, _, gp, _ = oblate_spheroid_d_functions(xx + dd)
    hm, _, gm, _ = oblate_spheroid_d_functions(xx - dd)
    _, dh, _, dg = oblate_spheroid_d_functions(xx)
    fh = np.max(np.abs((hp - hm) / (2 * dd) - dh) / np.abs(dh))
    fg = np.max(np.abs((gp - gm) / (2 * dd) - dg) / np.abs(dg))
    print(f"[IMP-6 Br] xi in [5, 20]: closed vs series h {rh:.3e}, g {rg:.3e}; "
          f"h' vs FD {fh:.3e}, g' vs FD {fg:.3e} (printed)")
    assert rh <= 1e-13
    assert rg <= 1e-13


# --------------------------------------------------------------- (Bd) (Tg)
@pytest.mark.parametrize("aspect", [2.0, 5.0])
@pytest.mark.parametrize("drive", DRIVES)
def test_boundary_and_tangential(aspect, drive):
    c, f, xi0 = _geom(aspect)
    pts = _surface_points(A, c)
    ax = _axis(drive)
    # (Bd) the exterior expression on the surface (the public potential masks
    # points within 1e-12 of the surface as interior, which would pass
    # trivially): ξ map + D-functions directly.
    xi = oblate_spheroidal_xi(pts, f)
    h, _, g, _ = oblate_spheroid_d_functions(xi)
    h0, _, g0, _ = oblate_spheroid_d_functions(xi0)
    ratio = h / h0[0] if ax == 2 else g / g0[0]
    phi = -E0 * pts[:, ax] * (1.0 - ratio)
    bd = np.max(np.abs(phi)) / (E0 * A)
    # (Tg) tangential E from the public field on the surface.
    E = pec_oblate_spheroid_field(pts, semi_major=A, semi_minor=c, drive=drive,
                                  interior="continue")
    n = pts / np.array([A**2, A**2, c**2])
    n /= np.linalg.norm(n, axis=1)[:, None]
    en = np.sum(E * n, axis=1)
    tg = np.max(np.linalg.norm(E - en[:, None] * n, axis=1)) / E0
    print(f"[IMP-6 Bd/Tg] aspect={aspect} {drive}: max|phi|/(E0 a)={bd:.3e} "
          f"max|E_t|/E0={tg:.3e} max|E_n|/E0={np.max(np.abs(en)):.4f}")
    assert bd <= 1e-12
    assert tg <= 1e-12


# ------------------------------------------------------------ (Rim) (Pole)
@pytest.mark.parametrize("aspect", ASPECTS)
def test_rim_and_pole(aspect):
    c, _, _ = _geom(aspect)
    lx, lz = oblate_spheroid_depolarisation(_ecc(aspect))
    rim = np.array([[A, 0, 0], [-A, 0, 0]], dtype=float)
    side = np.array([[0, A, 0], [0, -A, 0]], dtype=float)
    pole = np.array([[0, 0, c], [0, 0, -c]], dtype=float)
    Er = pec_oblate_spheroid_field(rim, semi_major=A, semi_minor=c, drive="inplane")
    Es = pec_oblate_spheroid_field(side, semi_major=A, semi_minor=c, drive="inplane")
    Ep = pec_oblate_spheroid_field(pole, semi_major=A, semi_minor=c, drive="axial")
    r_rim = np.max(np.abs(Er[:, 0] * lx - 1.0))
    r_side = np.max(np.linalg.norm(Es, axis=1)) / E0
    r_pole = np.max(np.abs(Ep[:, 2] * lz - 1.0))
    print(f"[IMP-6 Rim/Pole] aspect={aspect}: E_x(rim)={Er[0, 0]:.6f} (1/L_x="
          f"{1/lx:.6f}) |E_x L_x-1|={r_rim:.3e}; |E(0,±a,0)|={r_side:.3e}; "
          f"E_z(pole)={Ep[0, 2]:.6f} (1/L_z={1/lz:.6f}) |E_z L_z-1|={r_pole:.3e}")
    assert r_rim <= 1e-12
    assert r_side <= 1e-12
    assert r_pole <= 1e-12


# ---------------------------------------------------------------- (F) (Lap)
def _volume_points(aspect, seed):
    c, f, xi0 = _geom(aspect)
    rng = np.random.default_rng(seed)
    xi = rng.uniform(1.2 * xi0 + 0.05, 3.0, 40)
    eta = rng.uniform(-0.95, 0.95, 40)
    ang = rng.uniform(0.0, 2 * np.pi, 40)
    return c, _oblate_to_cart(xi, eta, ang, f)


@pytest.mark.parametrize("aspect", [2.0, 5.0])
@pytest.mark.parametrize("drive", DRIVES)
def test_field_vs_fd_gradient(aspect, drive):
    c, pts = _volume_points(aspect, 403)

    def phi(p):
        return pec_oblate_spheroid_potential(p, semi_major=A, semi_minor=c, drive=drive)

    def grad_fd(delta):
        out = np.zeros_like(pts)
        for ax in range(3):
            d = np.zeros(3)
            d[ax] = delta
            out[:, ax] = (phi(pts + d) - phi(pts - d)) / (2 * delta)
        return -out

    E = pec_oblate_spheroid_field(pts, semi_major=A, semi_minor=c, drive=drive)
    d0, d1, d2 = grad_fd(1e-3 * A), grad_fd(2e-3 * A), grad_fd(4e-3 * A)
    e1 = np.max(np.linalg.norm(E - d1, axis=1))
    e2 = np.max(np.linalg.norm(E - d2, axis=1))
    rich = np.max(np.linalg.norm(E - (4 * d1 - d2) / 3, axis=1)) / E0
    rich_half = np.max(np.linalg.norm(E - (4 * d0 - d1) / 3, axis=1)) / E0
    ratio = e2 / e1
    r4 = rich / rich_half
    print(f"[IMP-6 F] aspect={aspect} {drive}: |E-FD| {e1:.3e} (2e-3a) "
          f"{e2:.3e} (4e-3a) ratio {ratio:.3f}; Richardson residual/E0 {rich:.3e} "
          f"(2e-3a/4e-3a pair), {rich_half:.3e} (1e-3a/2e-3a pair), ratio {r4:.2f}")
    assert 3.8 <= ratio <= 4.2
    if aspect == 2.0:
        assert rich <= 1e-7
    else:
        # Bound moved, with measurement (W1 20261012T183401Z, this file's
        # first run): at aspect 5 the pre-registered pair's Richardson
        # residual is 3.55e-7 (axial) / 5.91e-7 (in-plane) against the
        # brief's 1e-7, while the second-order ratio reads 3.999 — the
        # residual is the O(δ⁴) truncation term at the sample points nearest
        # the rim (ξ = 1.2ξ₀ + 0.05 sits ≈ 0.02a off a rim whose curvature
        # radius is c²/a = 0.04a), not a gradient defect (a defect reads
        # ratio ≈ 1).  The 1e-7 bound is kept, one δ-halving finer, and the
        # fourth-order signature of the residual itself is asserted.
        assert rich_half <= 1e-7
        assert 14.0 <= r4 <= 18.0


@pytest.mark.parametrize("aspect", [2.0, 5.0])
@pytest.mark.parametrize("drive", DRIVES)
def test_fd_laplacian(aspect, drive):
    c, pts = _volume_points(aspect, 4030)

    def phi(p):
        return pec_oblate_spheroid_potential(p, semi_major=A, semi_minor=c, drive=drive)

    def lap(delta):
        out = -6.0 * phi(pts)
        for ax in range(3):
            d = np.zeros(3)
            d[ax] = delta
            out = out + phi(pts + d) + phi(pts - d)
        return out / delta**2

    l1, l2 = lap(2e-3 * A), lap(4e-3 * A)
    ratio = np.max(np.abs(l2)) / np.max(np.abs(l1))
    rich = np.max(np.abs((4 * l1 - l2) / 3)) * A / E0
    print(f"[IMP-6 Lap] aspect={aspect} {drive}: max|lap| a/E0 {np.max(np.abs(l1)):.3e} "
          f"(2e-3a) {np.max(np.abs(l2)):.3e} (4e-3a) ratio {ratio:.3f}; "
          f"Richardson residual a/E0 {rich:.3e}")
    assert 3.8 <= ratio <= 4.2
    assert rich <= 1e-7


# ------------------------------------------------------------- (Far) (c2)
@pytest.mark.parametrize("aspect", [2.0, 5.0])
@pytest.mark.parametrize("drive", DRIVES)
def test_far_dipole_and_control_c2(aspect, drive):
    c, f, xi0 = _geom(aspect)
    ax = _axis(drive)
    li = _l_of(aspect, drive)
    p = E0 * A * A * c / (3.0 * li)

    def phi_s(pts):
        return pec_oblate_spheroid_potential(pts, semi_major=A, semi_minor=c,
                                             drive=drive, scattered_only=True)

    vals = [_l1_projection(phi_s, R * A, drive) for R in (2.0, 4.0)]
    rel = max(abs(v - p) / p for v in vals)
    print(f"[IMP-6 Far] aspect={aspect} {drive}: p={p:.12f}; projections "
          f"{vals[0]:.12f} {vals[1]:.12f}; max rel {rel:.3e}")
    assert rel <= 1e-12

    # (c2), 10:30 ruling (Q2): the regular solution's ratio in place of
    # D(ξ)/D(ξ₀) — ξ/ξ₀ (axial), √(1+ξ²)/√(1+ξ₀²) (in-plane).  φ_s then grows
    # as x_i·ξ ∼ r², so the projection grows ∝ R⁴ (leading order: ratio ≈ 16).
    def phi_reg(pts):
        xi = oblate_spheroidal_xi(pts, f)
        if ax == 2:
            r = xi / xi0
        else:
            r = np.sqrt(1.0 + xi**2) / np.sqrt(1.0 + xi0**2)
        return E0 * pts[:, ax] * r

    # (c2′), the brief's original wording: D ≡ 1, i.e. φ_s = E₀x_i (R³, ratio 8).
    def phi_one(pts):
        return E0 * pts[:, ax]

    w = [_l1_projection(phi_reg, R * A, drive) for R in (2.0, 4.0)]
    w1 = [_l1_projection(phi_one, R * A, drive) for R in (2.0, 4.0)]
    ratio, ratio1 = w[1] / w[0], w1[1] / w1[0]
    print(f"[IMP-6 c2] aspect={aspect} {drive}: regular ratio in place of D/D0: "
          f"4a/2a projection ratio {ratio:.4f} (predicted ~16, R^4 leading order; "
          f"10:30 ruling text says ~8); D=1 variant {ratio1:.4f} (derived 8, R^3)")
    assert abs(ratio - 1.0) > 1e-3
    assert abs(ratio1 - 1.0) > 1e-3


# ------------------------------------------------------------- (Sph) (Disc)
def test_sphere_limit():
    e = 1e-2
    lx, lz = oblate_spheroid_depolarisation(e)
    cz = (lz - 1.0 / 3.0) / e**2
    cx = (lx - 1.0 / 3.0) / e**2
    print(f"[IMP-6 Sph] e=1e-2: (L_z-1/3)/e^2={cz:.8f} (2/15={2/15:.8f}); "
          f"(L_x-1/3)/e^2={cx:.8f} (-1/15={-1/15:.8f})")
    assert abs(cz - 2.0 / 15.0) <= 1e-3
    assert abs(cx + 1.0 / 15.0) <= 1e-3

    c = A * np.sqrt(1.0 - e**2)
    rng = np.random.default_rng(4031)
    n = 200
    r = rng.uniform(1.2 * A, 3.0 * A, n)
    mu = rng.uniform(-1.0, 1.0, n)
    ang = rng.uniform(0.0, 2 * np.pi, n)
    s = np.sqrt(1.0 - mu**2)
    pts = r[:, None] * np.column_stack([s * np.cos(ang), s * np.sin(ang), mu])
    for drive in DRIVES:
        ax = _axis(drive)
        phi = pec_oblate_spheroid_potential(pts, semi_major=A, semi_minor=c, drive=drive)
        sph = -E0 * pts[:, ax] * (1.0 - A**3 / r**3)
        dev = np.max(np.abs(phi - sph)) / (E0 * A)
        print(f"[IMP-6 Sph] {drive}: max|phi - phi_sphere|/(E0 a), r in [1.2a, 3a]: {dev:.3e}")
        assert dev <= 1e-3


def test_disc_limit():
    c = 1e-3 * A
    e = np.sqrt((1.0 - c) * (1.0 + c))
    xi0 = c / np.sqrt((A - c) * (A + c))
    lx, lz = oblate_spheroid_depolarisation(e)
    coef = (lx - np.pi * xi0 / 4.0) / xi0**2 + 1.0
    print(f"[IMP-6 Disc] c/a=1e-3: xi0={xi0:.6e} L_x={lx:.6e} L_z={lz:.9f}; "
          f"(L_x - pi xi0/4)/xi0^2 + 1 = {coef:.4e} (next term (pi/4)xi0 = "
          f"{np.pi*xi0/4:.4e}); p_x/a^3={c/(3*lx):.6f} (4/(3pi)={4/(3*np.pi):.6f})")
    assert abs(coef) <= 2e-3


# --------------------------------------------------------------------- (c1)
@pytest.mark.parametrize("aspect,pred_x,pred_z", [(2.0, 1.748, 0.329), (5.0, 3.784, 0.0744)])
def test_control_c1_prolate_formula(aspect, pred_x, pred_z):
    e = _ecc(aspect)
    c, _, _ = _geom(aspect)
    lx_obl, lz_obl = oblate_spheroid_depolarisation(e)
    lz_pro = prolate_spheroid_depolarisation(e)
    lx_pro = (1.0 - lz_pro) / 2.0
    Er = pec_oblate_spheroid_field(np.array([[A, 0.0, 0.0]]), semi_major=A,
                                   semi_minor=c, drive="inplane")
    rim = Er[0, 0] * lx_pro
    print(f"[IMP-6 c1] aspect={aspect}: E_x(rim)*L_x^pro = {rim:.4f} "
          f"(predicted {pred_x}); L_z^pro/L_z^obl = {lz_pro/lz_obl:.4f} "
          f"(predicted {pred_z}, printed)")
    assert abs(rim - 1.0) > 0.1
