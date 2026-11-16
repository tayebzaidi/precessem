"""`MAT-3` step 1, numpy anchors N1–N6 for Debye / Cole–Cole ε̂(ω).

Closed-form identities only; no FEM.  Convention ``e^{+jωt}``: a passive pole
has ``Im ε̂ < 0``.  N4 is the anchor that distinguishes a conjugated pole
(N3's circle is symmetric under conjugation, N4's depressed-arc centre is not).
"""

from __future__ import annotations

import numpy as np

from fem_em_solver.materials import ColeColeMaterial
from fem_em_solver.utils.analytical import complex_permittivity
from fem_em_solver.utils.constants import EPSILON_0

FREQS = np.logspace(3, 11, 400)

EPS_INF = 50.0
DELTA_EPS = 30.0
TAU = 1.0 / (2.0 * np.pi * 100e6)
SIGMA = 0.5


def _rel(a, b):
    return np.max(np.abs(a - b) / np.abs(b))


def test_n1_cole_cole_alpha_zero_equals_debye():
    cc = ColeColeMaterial(EPS_INF, poles=((DELTA_EPS, TAU, 0.0),), sigma_ionic=SIGMA)
    eps = cc.complex_permittivity(FREQS)
    omega = 2.0 * np.pi * FREQS
    wt = omega * TAU
    debye = (EPS_INF + DELTA_EPS * (1.0 - 1j * wt) / (1.0 + wt**2)
             - 1j * SIGMA / (omega * EPSILON_0))
    err = _rel(eps, debye)
    print(f"\n[MAT-3 N1] max rel |CC(alpha=0) - Debye| = {err:.3e}")
    assert err <= 1e-14
    # The classmethod is the same object.
    deb = ColeColeMaterial.debye(EPS_INF, poles=((DELTA_EPS, TAU),), sigma_ionic=SIGMA)
    assert deb == cc


def test_n2_no_poles_equals_complex_permittivity():
    for eps_r, sigma in ((50.0, 0.5), (78.0, 0.7), (1.0, 0.0)):
        mat = ColeColeMaterial.debye(eps_r, poles=((0.0, TAU),), sigma_ionic=sigma)
        eps = mat.complex_permittivity(FREQS)
        ref = np.array([complex_permittivity(eps_r, sigma, f) for f in FREQS])
        err = _rel(eps, ref)
        print(f"\n[MAT-3 N2] eps_r={eps_r}, sigma={sigma}: max rel = {err:.3e}")
        assert err <= 1e-15


def test_n3_debye_semicircle():
    mat = ColeColeMaterial.debye(EPS_INF, poles=((DELTA_EPS, TAU),), sigma_ionic=0.0)
    eps = mat.complex_permittivity(FREQS)
    eps_s = EPS_INF + DELTA_EPS
    dev = np.max(np.abs(np.abs(eps - 0.5 * (eps_s + EPS_INF)) - 0.5 * DELTA_EPS))
    peak = mat.complex_permittivity(1.0 / (2.0 * np.pi * TAU))
    print(f"\n[MAT-3 N3] max |radius - de/2| = {dev:.3e} (bound {1e-12 * DELTA_EPS:.1e}); "
          f"eps'' at wt=1 = {-peak.imag:.15f} (de/2 = {DELTA_EPS / 2})")
    assert dev <= 1e-12 * DELTA_EPS
    assert abs(-peak.imag - 0.5 * DELTA_EPS) <= 1e-12 * DELTA_EPS


def _arc_deviation(eps: np.ndarray, alpha: float) -> float:
    """Max |dist(ε̂, centre) − R| with ε″ = −Im ε̂ signed (not abs)."""
    eps_s = EPS_INF + DELTA_EPS
    centre_re = 0.5 * (eps_s + EPS_INF)
    centre_eps2 = -0.5 * DELTA_EPS * np.tan(0.5 * alpha * np.pi)
    radius = 0.5 * DELTA_EPS / np.cos(0.5 * alpha * np.pi)
    eps2 = -eps.imag  # signed ε″, e^{+jωt}: ε̂ = ε′ − jε″
    dist = np.hypot(eps.real - centre_re, eps2 - centre_eps2)
    return float(np.max(np.abs(dist - radius)))


def test_n4_cole_cole_depressed_arc():
    alpha = 0.1
    mat = ColeColeMaterial(EPS_INF, poles=((DELTA_EPS, TAU, alpha),), sigma_ionic=0.0)
    eps = mat.complex_permittivity(FREQS)
    dev = _arc_deviation(eps, alpha)
    # The conjugated pole (the e^{-jωt} sign) must fail the same identity:
    # this is what N4 exists to catch.
    dev_conj = _arc_deviation(np.conj(eps), alpha)
    print(f"\n[MAT-3 N4] alpha=0.1 arc deviation = {dev:.3e} "
          f"(bound {1e-12 * DELTA_EPS:.1e}); conjugated pole: {dev_conj:.3e}")
    assert dev <= 1e-12 * DELTA_EPS
    assert dev_conj > 1e-3 * DELTA_EPS


def test_n5_passivity():
    mats = (
        ColeColeMaterial(EPS_INF, poles=((DELTA_EPS, TAU, 0.0),), sigma_ionic=SIGMA),
        ColeColeMaterial(EPS_INF, poles=((DELTA_EPS, TAU, 0.1),), sigma_ionic=0.0),
        ColeColeMaterial(4.0, poles=((45.0, 7.23e-12, 0.1), (400.0, 3.18e-9, 0.1),
                                     (2.0e5, 159e-6, 0.2)), sigma_ionic=0.2),
    )
    for mat in mats:
        eps = mat.complex_permittivity(FREQS)
        min_loss = float(np.min(-eps.imag))
        print(f"\n[MAT-3 N5] min(-Im eps) = {min_loss:.6e}")
        assert np.all(-eps.imag > 0.0)


def test_n6_sigma_epsilon_round_trip():
    mat = ColeColeMaterial(EPS_INF, poles=((DELTA_EPS, TAU, 0.1),), sigma_ionic=SIGMA)
    eps = mat.complex_permittivity(FREQS)
    sigma_eff, eps_r = mat.sigma_epsilon_at(FREQS)
    back = np.array([complex_permittivity(e, s, f)
                     for e, s, f in zip(eps_r, sigma_eff, FREQS)])
    err = _rel(back, eps)
    # Scalar path too, via the literal ``*sigma_epsilon_at(f)[::-1]`` spelling.
    err_scalar = max(
        abs(complex_permittivity(*mat.sigma_epsilon_at(f)[::-1], f)
            - mat.complex_permittivity(f)) / abs(mat.complex_permittivity(f))
        for f in (63.87e6, 127.74e6, 1e3, 1e11)
    )
    print(f"\n[MAT-3 N6] round-trip max rel = {err:.3e} (scalar {err_scalar:.3e})")
    assert err <= 1e-14
    assert err_scalar <= 1e-14
    hm = mat.as_homogeneous(127.74e6)
    s, e = mat.sigma_epsilon_at(127.74e6)
    assert (hm.sigma, hm.epsilon_r, hm.mu_r) == (s, e, 1.0)
