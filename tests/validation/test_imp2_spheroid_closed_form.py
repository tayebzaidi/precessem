"""IMP-2 step 0 (§9 item 136, chain I3): PEC prolate spheroid in a uniform
axial field — the quasi-static closed form, numpy only.

Asserted identities (derivations in ``utils/analytical.py``, IMP-2 block):
(T) tip identity 1 − ξ₀Q₁′/Q₁ = 1/L_z; (Q) Legendre ODE residual and branch
agreement; (Bd) φ = 0 on the body; (Lap) FD Laplacian with a measured
second-order ratio; (Far) P₁ projection equals the dipole a²c/(3L_z);
(Sph) the e → 0 sphere limit.  Negative controls (c1) sphere potential on the
spheroid surface, (c2) Q₁ → P₁ (the growing solution).
"""

import numpy as np
import pytest

from fem_em_solver.utils.analytical import (
    pec_prolate_spheroid_potential,
    prolate_legendre_q1,
    prolate_spheroid_depolarisation,
    prolate_spheroidal_coordinates,
)

C = 1.0
E0 = 1.0


def _ecc(aspect):
    return np.sqrt(1.0 - 1.0 / aspect**2)


def _surface_points(c, a, n=200):
    """Fibonacci points on the spheroid surface (x²+y²)/a² + z²/c² = 1."""
    i = np.arange(n) + 0.5
    mu = 1.0 - 2.0 * i / n
    phi = np.pi * (1.0 + np.sqrt(5.0)) * i
    s = np.sqrt(1.0 - mu**2)
    return np.column_stack([a * s * np.cos(phi), a * s * np.sin(phi), c * mu])


def _spheroidal_to_cart(xi, eta, phi, f):
    rho = f * np.sqrt((xi**2 - 1.0) * (1.0 - eta**2))
    return np.column_stack([rho * np.cos(phi), rho * np.sin(phi), f * xi * eta])


def _p1_projection(phi_fn, radius, n=200):
    """R²·(3/2)∫ φ(R, μ) μ dμ by Gauss–Legendre (axisymmetric, φ-angle 0.3)."""
    mu, w = np.polynomial.legendre.leggauss(n)
    s = np.sqrt(1.0 - mu**2)
    pts = radius * np.column_stack([s * np.cos(0.3), s * np.sin(0.3), mu])
    return radius**2 * 1.5 * np.sum(w * mu * phi_fn(pts))


# --------------------------------------------------------------------- (T)
@pytest.mark.parametrize("aspect", [1.01, 2.0, 5.0, 20.0])
def test_tip_identity(aspect):
    e = _ecc(aspect)
    xi0 = 1.0 / e
    lz = prolate_spheroid_depolarisation(e)
    q, dq, _ = prolate_legendre_q1(xi0)
    tip = 1.0 - xi0 * dq[0] / q[0]
    rel = abs(tip * lz - 1.0)
    print(f"[IMP-2 T] aspect={aspect}: L_z={lz:.6f} 1/L_z={1/lz:.6f} "
          f"E_tip/E0={tip:.12f} |E_tip*L_z-1|={rel:.3e}")
    assert rel <= 1e-12
    # Same tip field by FD on the potential just outside the tip (printed).
    a = C / aspect
    h = 1e-6 * C
    z = np.array([[0.0, 0.0, C + h], [0.0, 0.0, C + 2 * h]])
    ph = pec_prolate_spheroid_potential(z, semi_major=C, semi_minor=a)
    print(f"[IMP-2 T] aspect={aspect}: FD tip field {-(ph[1]-ph[0])/h:.6f} "
          f"(one-sided, O(h) biased)")


def test_printed_depolarisation_vs_osborn():
    for aspect, osborn in [(2.0, 0.1736), (5.0, 0.0558)]:
        lz = prolate_spheroid_depolarisation(_ecc(aspect))
        print(f"[IMP-2 L] aspect={aspect}: L_z={lz:.6f} (Osborn {osborn}) "
              f"1/L_z={1/lz:.4f} L_x={(1-lz)/2:.6f}")
        assert abs(lz - osborn) < 1e-4  # Osborn's 4-figure table entry


# --------------------------------------------------------------------- (Q)
def test_legendre_ode_and_branches():
    xi = np.geomspace(1.001, 1e3, 50)
    q, dq, d2q = prolate_legendre_q1(xi)
    res = (1.0 - xi**2) * d2q - 2.0 * xi * dq + 2.0 * q
    scale = np.abs((1.0 - xi**2) * d2q) + np.abs(2.0 * xi * dq) + np.abs(2.0 * q)
    rel = np.max(np.abs(res) / scale)
    # Q₁′ against a central difference of Q₁ (independent of the dq formula).
    hh = 1e-5 * (xi - 1.0)
    fd = (prolate_legendre_q1(xi + hh)[0] - prolate_legendre_q1(xi - hh)[0]) / (2 * hh)
    fd_rel = np.max(np.abs(fd - dq) / np.abs(dq))
    print(f"[IMP-2 Q] ODE residual max rel {rel:.3e}; Q1' vs FD max rel {fd_rel:.3e}")
    assert rel <= 1e-10
    assert fd_rel <= 1e-6

    # Branch pairs at the switch points.
    xs = np.array([10.0, 10.5, 12.0])
    qc, dqc, _ = prolate_legendre_q1(xs, branch="closed")
    qs, dqs, _ = prolate_legendre_q1(xs, branch="series")
    dq_rel = np.max(np.abs(qc - qs) / np.abs(qs))
    ddq_rel = np.max(np.abs(dqc - dqs) / np.abs(dqs))
    es = [0.05, 0.049]
    l_rel = max(
        abs(prolate_spheroid_depolarisation(e, branch="closed")
            - prolate_spheroid_depolarisation(e, branch="series"))
        / prolate_spheroid_depolarisation(e, branch="series")
        for e in es
    )
    print(f"[IMP-2 Q] branch agreement: Q1 {dq_rel:.3e}, Q1' {ddq_rel:.3e}, "
          f"L_z {l_rel:.3e}")
    assert dq_rel <= 1e-13
    assert ddq_rel <= 1e-13
    assert l_rel <= 1e-13


# -------------------------------------------------------------------- (Bd)
@pytest.mark.parametrize("aspect", [2.0, 5.0])
def test_boundary_and_control_c1(aspect):
    a = C / aspect
    pts = _surface_points(C, a)
    # Evaluate the exterior expression on the surface: nudge outward by 1 ulp
    # scale is unnecessary — the mask is strict (< 1 − 1e-12), so use the
    # formula directly through scattered_only − (−E0 z).
    e = _ecc(aspect)
    f = C * e
    xi, eta = prolate_spheroidal_coordinates(pts, f)
    q0 = prolate_legendre_q1(1.0 / e)[0][0]
    q = prolate_legendre_q1(xi)[0]
    phi = -E0 * f * eta * (xi - (1.0 / e) * q / q0)
    bd = np.max(np.abs(phi)) / (E0 * C)
    # And through the public function on points 1e-9 c outside the surface.
    n = pts / np.array([a**2, a**2, C**2])
    n /= np.linalg.norm(n, axis=1)[:, None]
    phi_pub = pec_prolate_spheroid_potential(pts + 1e-9 * C * n, semi_major=C, semi_minor=a)
    print(f"[IMP-2 Bd] aspect={aspect}: max|phi|/(E0 c) = {bd:.3e}; "
          f"public fn at +1e-9c: {np.max(np.abs(phi_pub)):.3e}")
    assert bd <= 1e-12
    assert np.max(np.abs(phi_pub)) <= 1e-7  # O(1e-9 c · E_tip)
    if aspect == 2.0:
        r = np.linalg.norm(pts, axis=1)
        cos_t = pts[:, 2] / r
        sph = -E0 * r * cos_t * (1.0 - C**3 / r**3)
        c1 = np.max(np.abs(sph)) / (E0 * C)
        print(f"[IMP-2 c1] sphere(radius c) potential on aspect-2 surface: {c1:.4f} "
              f"(predicted 0.5-1.5)")
        assert c1 > 1e-3


# ------------------------------------------------------------------- (Lap)
@pytest.mark.parametrize("aspect", [2.0, 5.0])
def test_fd_laplacian(aspect):
    a = C / aspect
    e = _ecc(aspect)
    f = C * e
    xi0 = 1.0 / e
    rng = np.random.default_rng(136)
    xi = rng.uniform(1.2 * xi0, 4.0 * xi0, 40)
    eta = rng.uniform(-0.95, 0.95, 40)
    ang = rng.uniform(0.0, 2 * np.pi, 40)
    pts = _spheroidal_to_cart(xi, eta, ang, f)

    def phi(p):
        return pec_prolate_spheroid_potential(p, semi_major=C, semi_minor=a)

    def lap(delta):
        out = -6.0 * phi(pts)
        for ax in range(3):
            d = np.zeros(3)
            d[ax] = delta
            out = out + phi(pts + d) + phi(pts - d)
        return np.abs(out) / delta**2

    # Residual bound at the item's δ = 1e-4·c (pre-registered, unchanged).
    d = 1e-4 * C
    l1 = lap(d)
    l2 = lap(2 * d)
    res = np.max(l1) * C / E0
    print(f"[IMP-2 Lap] aspect={aspect}: max|lap| c/E0 = {res:.3e} (delta=1e-4c), "
          f"{np.max(l2)*C/E0:.3e} (2 delta); ratio {np.max(l2)/np.max(l1):.3f}")
    assert res <= 1e-5
    # Ratio pair moved, with measurement (W1 20261001T142458Z, this file's first
    # run): at δ = 1e-4·c the residual is round-off, not truncation — ratio
    # 0.276 (aspect 2) / 0.281 (aspect 5) ≈ 1/4, the ∝ 1/δ² signature of
    # ε·|φ|/δ² (|φ| ≈ E₀c, ε ≈ 2e-16 → ≈ 1e-6, as measured 8.4e-7 / 1.4e-6).
    # A non-harmonic φ would read ratio ≈ 1 at any δ.  The second-order ratio
    # is therefore asserted where truncation dominates; the δ ladder is printed.
    ladder = {}
    for k in (1e-4, 1e-3, 2e-3, 4e-3, 8e-3):
        ladder[k] = np.max(lap(k * C)) * C / E0
    print(f"[IMP-2 Lap] aspect={aspect}: ladder " +
          ", ".join(f"{k:g}c:{v:.3e}" for k, v in ladder.items()))
    ratio = ladder[4e-3] / ladder[2e-3]
    print(f"[IMP-2 Lap] aspect={aspect}: truncation ratio (4e-3c / 2e-3c) {ratio:.3f}")
    assert 3.8 <= ratio <= 4.2


# ------------------------------------------------------------------- (Far)
@pytest.mark.parametrize("aspect", [2.0, 5.0])
def test_far_dipole_and_control_c2(aspect):
    a = C / aspect
    e = _ecc(aspect)
    f = C * e
    xi0 = 1.0 / e
    lz = prolate_spheroid_depolarisation(e)
    p = E0 * a**2 * C / (3.0 * lz)

    def phi_s(pts):
        return pec_prolate_spheroid_potential(
            pts, semi_major=C, semi_minor=a, scattered_only=True)

    vals = [_p1_projection(phi_s, R * C) for R in (2.0, 4.0)]
    rel = max(abs(v - p) / p for v in vals)
    print(f"[IMP-2 Far] aspect={aspect}: p={p:.12f}; projections {vals[0]:.12f} "
          f"{vals[1]:.12f}; max rel {rel:.3e}")
    assert rel <= 1e-12

    q0 = prolate_legendre_q1(xi0)[0][0]

    def phi_wrong(pts):  # Q₁ → P₁ = ξ: the growing solution
        xi, eta = prolate_spheroidal_coordinates(pts, f)
        return E0 * f * eta * xi0 * xi / q0

    w = [_p1_projection(phi_wrong, R * C) for R in (2.0, 4.0)]
    ratio = w[1] / w[0]
    print(f"[IMP-2 c2] P1 (growing) in place of Q1: 4c/2c ratio {ratio:.4f} (derived 8)")
    assert abs(ratio - 1.0) > 1e-3


# ------------------------------------------------------------------- (Sph)
def test_sphere_limit():
    e = 1e-2
    lz = prolate_spheroid_depolarisation(e)
    coef = (1.0 / lz - 3.0) / e**2
    print(f"[IMP-2 Sph] e=1e-2: (1/L_z-3)/e^2 = {coef:.8f} (6/5)")
    assert abs(coef - 1.2) <= 1e-3

    a = C * np.sqrt(1.0 - e**2)
    rng = np.random.default_rng(1360)
    n = 200
    r = rng.uniform(1.2 * C, 3.0 * C, n)
    mu = rng.uniform(-1.0, 1.0, n)
    ang = rng.uniform(0.0, 2 * np.pi, n)
    s = np.sqrt(1.0 - mu**2)
    pts = r[:, None] * np.column_stack([s * np.cos(ang), s * np.sin(ang), mu])
    phi = pec_prolate_spheroid_potential(pts, semi_major=C, semi_minor=a)
    sph = -E0 * r * mu * (1.0 - C**3 / r**3)
    dev = np.max(np.abs(phi - sph)) / (E0 * C)
    print(f"[IMP-2 Sph] max|phi - phi_sphere|/(E0 c), r in [1.2c, 3c]: {dev:.3e}")
    assert dev <= 1e-3
