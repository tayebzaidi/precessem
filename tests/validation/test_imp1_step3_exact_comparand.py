"""IMP-1 step 3(a) — the exact finite-k l = 1 TM comparand (§9 item 135).

A PEC sphere (radius a) in a homogeneous lossy host (k = k₀√ε_c, principal
root, Im k < 0, e^{+jωt}) under a uniform ẑ drive.  The exact field is the
l = 1 TM pair

    E_reg = (3/2)[2cosθ j₁(x)/x r̂ − sinθ (x j₁)′/x θ̂]          (→ ẑ, k → 0)
    E_out = −j(ka)³[the same with h₁⁽²⁾]                  (→ the a³ dipole)
    E = E₀(E_reg + β_k E_out),  β_k = −(3/2)(x j₁)′/[−j x³ (x h₁⁽²⁾)′] at x = ka

(Jin 3e App. E.2, E.24–E.31), implemented additively in
``fem_em_solver.utils.analytical`` (``pec_sphere_l1_tm_beta`` /
``pec_sphere_l1_tm_field``).  Numpy only: no FEM, no rung, no band.

Anchors (asserted): (X) the real-k identity with `TH-5`'s verified dipole;
(P) n × E = 0 on r = a; (H) the curl-curl residual at |k_H|a = 1 scales as δ²
and is ≤ 1e-4; (Q) the k → 0 limit is the quasi-static field at second order;
(Φ0) the static shell mean 1 + 2(a/c)³.  Negative controls (c1)–(c3) asserted
only > 0.5.  Printed: β_k − 1, `IMP-1`'s Re-E LSQ β on the exact field,
max|Im E|/|Re E|, the box-face departure from the quasi-static field and
ρ_model — I2b's pre-registered prediction beside the banked ρ.
"""

from __future__ import annotations

import numpy as np
from mpi4py import MPI

from fem_em_solver.utils.analytical import (
    complex_permittivity,
    pec_sphere_l1_tm_beta,
    pec_sphere_l1_tm_field,
    pec_sphere_l1_tm_parts,
)
from fem_em_solver.utils.constants import C_0
from tests.validation.test_imp1_pec_sphere_lossy_host import (
    BOX_HALF_WIDTH,
    E0,
    EPSILON_R_HOST,
    FREQUENCY_HZ,
    SHELL_OUTER_FACTOR,
    SIGMA_HIGH,
    SPHERE_RADIUS,
    _dipole_basis,
    _probe_shells,
)
from tests.validation.test_th5_dipole_sphere import _dipole_e, _fd_curl

A = SPHERE_RADIUS
EPS_C = complex_permittivity(EPSILON_R_HOST, SIGMA_HIGH, FREQUENCY_HZ)
K0 = 2.0 * np.pi * FREQUENCY_HZ / C_0
K_FIXTURE = K0 * np.sqrt(EPS_C)  # principal root: Re > 0, Im < 0
K_H = K_FIXTURE / (abs(K_FIXTURE) * A)  # arg k kept, |k_H| a = 1
DELTA_H = 1e-3 * A
PHI_STATIC = 1.0 + 2.0 / SHELL_OUTER_FACTOR**3  # 1.592592593
BANKED_RHO_PERCENT = (0.4473, 0.5694)  # `20261001T033100Z_IMP-1-land-W1.log:545`, a/10 / a/12


def _say(msg):
    if MPI.COMM_WORLD.rank == 0:
        print(msg, flush=True)


def _fibonacci_sphere(n, radius):
    idx = np.arange(n) + 0.5
    polar = np.arccos(1.0 - 2.0 * idx / n)
    azim = np.pi * (1.0 + 5.0**0.5) * idx
    return radius * np.column_stack(
        [np.cos(azim) * np.sin(polar), np.sin(azim) * np.sin(polar), np.cos(polar)]
    )


def _off_axis_points(n, r_lo, r_hi, seed):
    rng = np.random.default_rng(seed)
    r = rng.uniform(r_lo, r_hi, n)
    theta = rng.uniform(np.arcsin(0.2), np.pi - np.arcsin(0.2), n)
    phi = rng.uniform(0.0, 2.0 * np.pi, n)
    return np.column_stack(
        [r * np.sin(theta) * np.cos(phi), r * np.sin(theta) * np.sin(phi), r * np.cos(theta)]
    )


def _qs_field(points, beta=1.0):
    e = beta * _dipole_basis(points).astype(complex)
    e[:, 2] += E0
    return e


def _curlcurl_residual(field, k, points, delta):
    cc = _fd_curl(lambda p: _fd_curl(field, p, delta), points, delta)
    k2e = k**2 * field(points)
    return float(np.max(np.linalg.norm(cc - k2e, axis=1) / np.linalg.norm(k2e, axis=1)))


def _tangential_max(field_values, points):
    r_hat = points / np.linalg.norm(points, axis=1)[:, None]
    return float(np.max(np.linalg.norm(np.cross(r_hat, field_values), axis=1)))


def _shell_mean(field, n_r=64, n_mu=32):
    """⟨|E|²⟩/E₀² over a < r < c (c = 1.5a); the field is axisymmetric in φ."""
    c = SHELL_OUTER_FACTOR * A
    xr, wr = np.polynomial.legendre.leggauss(n_r)
    xm, wm = np.polynomial.legendre.leggauss(n_mu)
    r = 0.5 * (c - A) * xr + 0.5 * (c + A)
    wr = 0.5 * (c - A) * wr
    rr, mm = np.meshgrid(r, xm, indexing="ij")
    ww = np.outer(wr, wm) * rr**2
    pts = np.column_stack(
        [(rr * np.sqrt(1.0 - mm**2)).ravel(), np.zeros(rr.size), (rr * mm).ravel()]
    )
    e2 = np.sum(np.abs(field(pts)) ** 2, axis=1).reshape(rr.shape)
    num = 2.0 * np.pi * np.sum(ww * e2)
    vol = 4.0 * np.pi / 3.0 * (c**3 - A**3)
    return float(num / (E0**2 * vol))


def _exact(k, beta=None):
    return lambda p: pec_sphere_l1_tm_field(p, k=k, radius=A, e0=E0, beta=beta)


def _box_face_points(n=11):
    b = BOX_HALF_WIDTH
    s = np.linspace(-b, b, n)
    u, v = np.meshgrid(s, s, indexing="ij")
    u, v = u.ravel(), v.ravel()
    faces = []
    for axis in range(3):
        for sign in (-1.0, 1.0):
            p = np.zeros((u.size, 3))
            others = [i for i in range(3) if i != axis]
            p[:, axis] = sign * b
            p[:, others[0]] = u
            p[:, others[1]] = v
            faces.append(p)
    return np.vstack(faces)


def test_fixture_wavenumber():
    assert EPSILON_R_HOST == 78.0 and SIGMA_HIGH == 0.57 and FREQUENCY_HZ == 64e6
    assert K_FIXTURE.real > 0.0 and K_FIXTURE.imag < 0.0
    assert abs(abs(K_H) * A - 1.0) < 1e-14
    _say(
        f"\n[IMP-1 3a] eps_c {EPS_C:.6f}; k {K_FIXTURE:.6f} 1/m; |k|a {abs(K_FIXTURE) * A:.6f}; "
        f"arg k {np.degrees(np.angle(K_FIXTURE)):.4f} deg"
    )


def test_x_real_k_identity_with_th5_dipole():
    k, a = 2.0 * np.pi, 0.1
    pts = _off_axis_points(40, 0.12, 0.6, seed=135)
    e_reg, e_out = pec_sphere_l1_tm_parts(pts, k=k, radius=a)
    d = _dipole_e(pts, k=k)
    ref_out = 1j * (k * a) ** 3 * d
    ref_reg = -1.5 * np.real(d)
    err_out = float(
        np.max(np.linalg.norm(e_out - ref_out, axis=1) / np.linalg.norm(ref_out, axis=1))
    )
    err_reg = float(
        np.max(np.linalg.norm(e_reg - ref_reg, axis=1) / np.linalg.norm(ref_reg, axis=1))
    )
    _say(
        f"\n[IMP-1 3a] (X) k = 2pi, a = 0.1, 40 points: |E_out - j(ka)^3 dipole|/|.| {err_out:.3e}; "
        f"|E_reg + 1.5 Re dipole|/|.| {err_reg:.3e} (bound 1e-12)"
    )
    assert err_out <= 1e-12, f"(X) outgoing part {err_out:.3e}"
    assert err_reg <= 1e-12, f"(X) regular part {err_reg:.3e}"


def test_p_tangential_e_vanishes_on_sphere():
    pts = _fibonacci_sphere(200, A)
    results = {}
    for label, k in (("fixture", K_FIXTURE), ("k_H", K_H)):
        e = _exact(k)(pts)
        # the points must hit the exterior evaluation (the normal part ≈ 3cosθ survives)
        assert float(np.max(np.abs(e))) > 1.0, f"(P) {label}: points masked as interior"
        results[label] = _tangential_max(e, pts) / E0
    c2 = _tangential_max(_exact(K_FIXTURE, beta=0.0)(pts), pts) / E0
    _say(
        f"\n[IMP-1 3a] (P) max |r x E|/E0 on 200 points at r = a: fixture {results['fixture']:.3e}; "
        f"k_H {results['k_H']:.3e} (bound 1e-12); (c2) beta = 0: {c2:.4f} (asserted > 0.5)"
    )
    for label, val in results.items():
        assert val <= 1e-12, f"(P) {label}: {val:.3e}"
    assert c2 > 0.5, f"(c2) control {c2:.4f}"


def test_h_curl_curl_residual():
    pts = _off_axis_points(20, 1.1 * A, 3.4 * A, seed=1135)
    res = {d: _curlcurl_residual(_exact(K_H), K_H, pts, d) for d in (DELTA_H, 2.0 * DELTA_H)}
    ratio = res[2.0 * DELTA_H] / res[DELTA_H]
    fixture = _curlcurl_residual(_exact(K_FIXTURE), K_FIXTURE, pts, DELTA_H)
    c1 = _curlcurl_residual(_qs_field, K_H, pts, DELTA_H)
    k_wrong_h = np.conj(K_H)  # the wrong-sign ε_c field (k* on the principal root)
    c3 = _curlcurl_residual(_exact(k_wrong_h), K_H, pts, DELTA_H)
    _say(
        f"\n[IMP-1 3a] (H) max |curl curl E - k^2 E|/|k^2 E| at k_H, 1.1a..3.4a, 20 points: "
        f"delta 1e-3a {res[DELTA_H]:.3e}; 2e-3a {res[2 * DELTA_H]:.3e}; ratio {ratio:.3f} "
        f"(predicted 1-2e-5, bound 1e-4, ratio in [3.8, 4.2])"
    )
    _say(f"[IMP-1 3a] (H) at the fixture's k (printed only): {fixture:.3e}")
    _say(
        f"[IMP-1 3a] (c1) quasi-static field in (H): {c1:.4f} (predicted 1.000); "
        f"(c3) wrong-sign eps_c field: {c3:.4f} (predicted 1.798); asserted > 0.5"
    )
    assert 3.8 <= ratio <= 4.2, f"(H) ratio {ratio:.3f}"
    assert res[DELTA_H] <= 1e-4, f"(H) residual {res[DELTA_H]:.3e}"
    assert c1 > 0.5, f"(c1) {c1:.4f}"
    assert c3 > 0.5, f"(c3) {c3:.4f}"


def test_q_quasi_static_limit_second_order():
    pts = _probe_shells()
    qs = _qs_field(pts)
    err = {}
    for s in (1e-2, 5e-3):
        err[s] = float(np.max(np.linalg.norm(_exact(s * K_FIXTURE)(pts) - qs, axis=1)))
    ratio = err[1e-2] / err[5e-3]
    _say(
        f"\n[IMP-1 3a] (Q) max |E - (z + dipole_basis)| on the probe shells: s = 1e-2 "
        f"{err[1e-2]:.3e}; s = 5e-3 {err[5e-3]:.3e}; ratio {ratio:.3f} "
        f"(bound 1e-5, ratio in [3.8, 4.2])"
    )
    assert err[1e-2] <= 1e-5 and err[5e-3] <= 1e-5, f"(Q) {err}"
    assert 3.8 <= ratio <= 4.2, f"(Q) ratio {ratio:.3f}"


def test_phi0_static_shell_mean():
    phi = _shell_mean(_exact(1e-6 * K_FIXTURE))
    _say(
        f"\n[IMP-1 3a] (Phi0) shell mean at s = 1e-6: {phi:.12f} vs {PHI_STATIC:.12f}, "
        f"|diff| {abs(phi - PHI_STATIC):.3e} (bound 1e-10)"
    )
    assert abs(phi - PHI_STATIC) <= 1e-10


def test_printed_model_error_readout():
    beta_k = pec_sphere_l1_tm_beta(K_FIXTURE, A)
    pts = _probe_shells()
    e = _exact(K_FIXTURE)(pts)
    basis = _dipole_basis(pts)
    uniform = np.zeros_like(pts)
    uniform[:, 2] = E0
    beta_fit = float(np.sum((np.real(e) - uniform) * basis) / np.sum(basis * basis))
    imag_ratio = float(
        np.max(np.linalg.norm(np.imag(e), axis=1)) / np.max(np.linalg.norm(np.real(e), axis=1))
    )
    face = _box_face_points()
    face_dev = float(np.max(np.linalg.norm(_exact(K_FIXTURE)(face) - _qs_field(face), axis=1)))
    phi_exact = _shell_mean(_exact(K_FIXTURE))
    phi_qs_fit = 1.0 + 2.0 * beta_fit**2 / SHELL_OUTER_FACTOR**3
    rho_model = phi_exact / phi_qs_fit - 1.0
    _say(
        f"\n[IMP-1 3a] beta_k - 1 = {beta_k - 1:.4e} (|.| {abs(beta_k - 1) * 100:.4f} %; "
        f"predicted +5.7e-4 - 2.1e-3j, 0.22 %)"
    )
    _say(f"[IMP-1 3a] IMP-1's Re-E LSQ beta on the exact field: {beta_fit:.6f}")
    _say(
        f"[IMP-1 3a] max|Im E|/|Re E| on the probe shells: {imag_ratio:.3e} "
        f"(step 1 FEM record 1.21e-2; predicted 3e-3)"
    )
    _say(
        f"[IMP-1 3a] max|E_exact - E_qs| on the box faces ({face.shape[0]} points): {face_dev:.3e} "
        f"(predicted 3e-3)"
    )
    _say(f"[IMP-1 3a] Phi_exact {phi_exact:.9f}; Phi_qs(beta_fit) {phi_qs_fit:.9f}")
    _say(
        f"[IMP-1 3a] rho_model = {rho_model * 100:+.4f} % beside banked rho "
        f"{BANKED_RHO_PERCENT[0]} / {BANKED_RHO_PERCENT[1]} % (a/10 / a/12); predicted <= 0.35 %"
    )
    assert np.isfinite(rho_model)  # printed readout; the anchors above carry the gate
