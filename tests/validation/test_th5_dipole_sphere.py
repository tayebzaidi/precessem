"""`TH-5` step 1b-0: the Hertzian-dipole field and the first-order ABC's l = 1
reflection on a sphere, verified in numpy before any solve (PROJECT_PLAN §9
item 107; derivation in `docs/planning/chunks/TH-5.md`).

A z-directed Hertzian dipole, e^{+jωt}, ``x = k₀r``, amplitude ``C = 1``::

    E_r    = 2C cosθ (1/x² − j/x³) e^{−jx}
    E_θ    = C sinθ (j f/x) e^{−jx},        f = 1 − j/x − 1/x²
    η₀H_φ  = C sinθ (j g/x) e^{−jx},        g = 1 − j/x

The field is the standard small-dipole result (Balanis, *Antenna Theory*,
eqs. 4-8 / 4-10a–c — cited from memory, **not verified in-repo**; the
converted Jin text does not carry it), so this module verifies it
independently: (n1) it satisfies ∇×∇×E = k₀²E, and (n2) it radiates
4πC²/(3η₀k₀²).

On a sphere r = R the first-order ABC ``n × ∇ × E + jk₀ n × n × E = 0``
(Jin 3e §9.2 (9.44) at s = 1; §9.3 (9.54)) reads ``E_θ = η₀H_φ``.  With
``E = A E_out + B conj(E_out)`` (the incoming wave's magnetic field is
``−conj(H_out)``) the l = 1 reflection is::

    Γ_sph(x) = B/A = −e^{−2jx} / (2x² − 1 + 2jx),     |Γ_sph| = 1/√(1 + 4x⁴)

and, for a general surface impedance ``E_θ = z η₀ H_φ``,
``Γ(z) = (f − z g) e^{−2jx} / (f̄ + z ḡ)`` — zero at ``z = f/g``, the PEC value
``(f/f̄) e^{−2jx}`` at z = 0, ``√(1 + 4x⁴)`` in modulus at z = −1 (γ → −γ).
(n3) checks the condition on the assembled standing wave with H taken from a
central-difference curl of E, (n4) pins the numbers, and (n5) records why the
planar R(θ) of Jin (9.48) cannot band a box (its dipole-weighted rms over a
cube is ≈ 0.113 at any size).

Tolerances here are finite-difference truncation bounds, not physics bands.
"""

from __future__ import annotations

import numpy as np
from mpi4py import MPI

from fem_em_solver.core.boundary_terms import abc_first_order_reflection_te
from fem_em_solver.utils.constants import ETA_0

WAVELENGTH = 1.0
K0 = 2.0 * np.pi / WAVELENGTH
DELTA_N1 = 1e-3  # item 107's step for the curl-curl residual
# Measured 2026-09-30 (`20260930T171523Z_TH-5-step1b0-W1.log:88`): at item 107's
# δ = 1e-3 the nested central-difference residual is 3.614e-4, at 2e-3 it is
# 1.445e-3 (ratio 4.00, pure δ² truncation), i.e. over the 1e-4 bound by
# 3.6× — the item's "step-size finding", not a wrong formula.  The bound is
# unchanged; it is asserted at δ = 2.5e-4 (predicted 3.614e-4/16 = 2.3e-5), and
# the δ² ratio is asserted beside it.
DELTA_N1_FINE = 2.5e-4
DELTA_H = 1e-4  # step for H from a central-difference curl of E, (n3) and (H)


def _say(msg):
    if MPI.COMM_WORLD.rank == 0:
        print(msg, flush=True)


def _f(x):
    return 1.0 - 1j / x - 1.0 / x**2


def _g(x):
    return 1.0 - 1j / x


def _frame(points):
    """r, r̂, θ̂, φ̂, cosθ, sinθ at off-axis points (N, 3)."""
    px, py, pz = points[:, 0], points[:, 1], points[:, 2]
    r = np.linalg.norm(points, axis=1)
    rho = np.hypot(px, py)
    cos_t, sin_t = pz / r, rho / r
    cos_p, sin_p = px / rho, py / rho
    r_hat = points / r[:, None]
    t_hat = np.column_stack([cos_t * cos_p, cos_t * sin_p, -sin_t])
    p_hat = np.column_stack([-sin_p, cos_p, np.zeros_like(r)])
    return r, r_hat, t_hat, p_hat, cos_t, sin_t


def _dipole_e(points, k=K0):
    """Outgoing dipole E (Cartesian, complex), C = 1."""
    r, r_hat, t_hat, _, cos_t, sin_t = _frame(points)
    x = k * r
    ph = np.exp(-1j * x)
    e_r = 2.0 * cos_t * (1.0 / x**2 - 1j / x**3) * ph
    e_t = sin_t * (1j * _f(x) / x) * ph
    return e_r[:, None] * r_hat + e_t[:, None] * t_hat


def _dipole_eta_h(points, k=K0):
    """Outgoing dipole η₀H (Cartesian, complex), C = 1."""
    r, _, _, p_hat, _, sin_t = _frame(points)
    x = k * r
    h_p = sin_t * (1j * _g(x) / x) * np.exp(-1j * x)
    return h_p[:, None] * p_hat


def _gamma_sphere(x):
    """First-order ABC's l = 1 reflection on r = R, x = k₀R."""
    return -np.exp(-2j * x) / (2.0 * x**2 - 1.0 + 2j * x)


def _gamma_impedance(x, z):
    """l = 1 reflection for ``E_θ = z η₀ H_φ`` on r = R."""
    f, g = _f(x), _g(x)
    return (f - z * g) * np.exp(-2j * x) / (np.conj(f) + z * np.conj(g))


def _fd_curl(field, points, delta):
    """Central-difference curl of a vector callable at (N, 3) points."""
    d = []
    for i in range(3):
        step = np.zeros(3)
        step[i] = delta
        d.append((field(points + step) - field(points - step)) / (2.0 * delta))
    # d[i][:, j] = ∂F_j/∂x_i
    return np.column_stack(
        [
            d[1][:, 2] - d[2][:, 1],
            d[2][:, 0] - d[0][:, 2],
            d[0][:, 1] - d[1][:, 0],
        ]
    )


def _off_axis_points(n, r_lo, r_hi, seed):
    rng = np.random.default_rng(seed)
    r = rng.uniform(r_lo, r_hi, n)
    # polar angles with sinθ ≥ 0.2 (central differences lose digits near the axis)
    theta = rng.uniform(np.arcsin(0.2), np.pi - np.arcsin(0.2), n)
    phi = rng.uniform(0.0, 2.0 * np.pi, n)
    return np.column_stack(
        [r * np.sin(theta) * np.cos(phi), r * np.sin(theta) * np.sin(phi), r * np.cos(theta)]
    )


def _sphere_ring(radius, n=20, azimuth=0.3):
    theta = np.linspace(np.arcsin(0.2), np.pi - np.arcsin(0.2), n)
    return np.column_stack(
        [
            radius * np.sin(theta) * np.cos(azimuth),
            radius * np.sin(theta) * np.sin(azimuth),
            radius * np.cos(theta),
        ]
    )


def _abc_residual(gamma, x, delta):
    """max over 20 polar angles of |E_θ − η₀H_φ|/|E_θ| on r = x/k₀ for
    E = E_out + Γ conj(E_out), η₀H = (j/k₀) ∇×E by central differences."""
    pts = _sphere_ring(x / K0)

    def e_total(p):
        e = _dipole_e(p)
        return e + gamma * np.conj(e)

    eta_h = (1j / K0) * _fd_curl(e_total, pts, delta)
    _, _, t_hat, p_hat, _, _ = _frame(pts)
    e_t = np.sum(e_total(pts) * t_hat, axis=1)
    h_p = np.sum(eta_h * p_hat, axis=1)
    return float(np.max(np.abs(e_t - h_p) / np.abs(e_t)))


def test_n1_dipole_satisfies_vector_helmholtz():
    pts = _off_axis_points(20, 0.2, 0.5, seed=107)
    residuals = {}
    for delta in (DELTA_N1, 2.0 * DELTA_N1):
        cc = _fd_curl(lambda p: _fd_curl(_dipole_e, p, delta), pts, delta)
        k2e = K0**2 * _dipole_e(pts)
        residuals[delta] = float(
            np.max(np.linalg.norm(cc - k2e, axis=1) / np.linalg.norm(k2e, axis=1))
        )
    ratio = residuals[2.0 * DELTA_N1] / residuals[DELTA_N1]
    fine = float(
        np.max(
            np.linalg.norm(
                _fd_curl(lambda p: _fd_curl(_dipole_e, p, DELTA_N1_FINE), pts, DELTA_N1_FINE)
                - K0**2 * _dipole_e(pts),
                axis=1,
            )
            / np.linalg.norm(K0**2 * _dipole_e(pts), axis=1)
        )
    )
    _say(
        f"\n[TH-5 1b-0] (n1) max |curl curl E - k0^2 E|/|k0^2 E| over 20 points, "
        f"0.2 <= r <= 0.5: delta {DELTA_N1:g}: {residuals[DELTA_N1]:.3e}; delta "
        f"{2 * DELTA_N1:g}: {residuals[2 * DELTA_N1]:.3e} (ratio {ratio:.3f}, "
        f"delta^2 predicts 4); delta {DELTA_N1_FINE:g}: {fine:.3e}"
    )
    # Pure truncation: the residual scales as δ² (a wrong formula would leave a
    # δ-independent floor), and at a step under it the bound holds.
    assert 3.8 <= ratio <= 4.2, f"(n1) residual does not scale as delta^2: ratio {ratio:.3f}"
    assert fine <= 1e-4, f"(n1) curl-curl residual {fine:.3e} > 1e-4 at delta {DELTA_N1_FINE:g}"

    # (H) η₀H is the curl of E: ∇×E = −jk₀ η₀H — the H that (n2) integrates.
    fd = _fd_curl(_dipole_e, pts, DELTA_H)
    ref = -1j * K0 * _dipole_eta_h(pts)
    h_res = float(np.max(np.linalg.norm(fd - ref, axis=1) / np.linalg.norm(ref, axis=1)))
    _say(
        f"[TH-5 1b-0] (H) max |curl E + j k0 eta0 H|/|k0 eta0 H| at delta "
        f"{DELTA_H:g}: {h_res:.3e}"
    )
    assert h_res <= 1e-5, f"(H) analytic eta0 H is not (j/k0) curl E: {h_res:.3e}"


def test_n2_radiated_power():
    radius = 0.3
    u, w = np.polynomial.legendre.leggauss(64)  # u = cosθ; integrand ∝ (1 − u²)
    sin_t = np.sqrt(1.0 - u**2)
    pts = np.column_stack([radius * sin_t, np.zeros_like(u), radius * u])
    e = _dipole_e(pts)
    eta_h = _dipole_eta_h(pts)
    r_hat = pts / radius
    s_r = 0.5 * np.real(np.sum(np.cross(e, np.conj(eta_h)) * r_hat, axis=1)) / ETA_0
    power = 2.0 * np.pi * radius**2 * float(np.sum(w * s_r))
    closed = 4.0 * np.pi / (3.0 * ETA_0 * K0**2)
    rel = abs(power / closed - 1.0)
    _say(
        f"[TH-5 1b-0] (n2) radiated power at r = {radius}: {power:.12e} vs "
        f"4 pi C^2/(3 eta0 k0^2) = {closed:.12e}: {rel:.2e} relative"
    )
    assert rel <= 1e-10, f"(n2) power {power:.12e} vs {closed:.12e}: {rel:.2e}"


def test_n3_first_order_condition_on_standing_wave():
    for x in (np.pi, 2.0 * np.pi):
        gamma = _gamma_sphere(x)
        res = _abc_residual(gamma, x, DELTA_H)
        res_coarse = _abc_residual(gamma, x, 10.0 * DELTA_H)
        _say(
            f"[TH-5 1b-0] (n3) x = {x:.6f}: Gamma_sph = {gamma.real:+.10f} "
            f"{gamma.imag:+.10f}j, max |E_t - eta0 H_p|/|E_t| = {res:.3e} at delta "
            f"{DELTA_H:g} ({res_coarse:.3e} at {10 * DELTA_H:g})"
        )
        assert res <= 1e-5, f"(n3) x = {x}: first-order condition residual {res:.3e}"

    # Negative control: Γ = 0 misses by |f − g|/|f| (the ceiling, 0.106 at x = π).
    x = np.pi
    res0 = _abc_residual(0.0, x, DELTA_H)
    closed0 = abs(_f(x) - _g(x)) / abs(_f(x))
    _say(
        f"[TH-5 1b-0] (n3) control Gamma = 0, x = pi: residual {res0:.6f} vs "
        f"|f - g|/|f| = {closed0:.6f} (predicted 0.106)"
    )
    assert res0 > 0.05, f"(n3) control: Gamma = 0 satisfies the ABC to {res0:.3e}"


def test_n4_reflection_numbers():
    g_pi, g_2pi = _gamma_sphere(np.pi), _gamma_sphere(2.0 * np.pi)
    m_pi, m_2pi = abs(g_pi), abs(g_2pi)
    small = _gamma_sphere(1e-3)
    x = np.pi
    zero = _gamma_impedance(x, _f(x) / _g(x))
    same = abs(_gamma_impedance(x, 1.0) - g_pi)
    pec = _gamma_impedance(x, 0.0)
    pec_cf = _f(x) / np.conj(_f(x)) * np.exp(-2j * x)
    flipped = abs(_gamma_impedance(x, -1.0))
    _say(
        f"[TH-5 1b-0] (n4) |Gamma_sph(pi)| = {m_pi:.10f} (1/sqrt(1+4x^4) "
        f"{1 / np.sqrt(1 + 4 * np.pi**4):.10f}), |Gamma_sph(2 pi)| = {m_2pi:.10f}, "
        f"ratio {m_pi / m_2pi:.4f}"
    )
    _say(
        f"  |Gamma_sph(1e-3) - 1| = {abs(small - 1):.3e}; Gamma(z = f/g) = {abs(zero):.2e}; "
        f"|Gamma(z = 1) - Gamma_sph| = {same:.2e}; |Gamma_PEC| = {abs(pec):.15f}, "
        f"vs (f/fbar)e^(-2jx) {abs(pec - pec_cf):.2e}; |Gamma(gamma -> -gamma)| = "
        f"{flipped:.4f}"
    )
    assert abs(m_pi - 0.0505957066) <= 1e-9, f"(n4) |Gamma_sph(pi)| = {m_pi:.10f}"
    assert abs(m_2pi - 0.0126641323) <= 1e-9, f"(n4) |Gamma_sph(2pi)| = {m_2pi:.10f}"
    assert abs(m_pi / m_2pi - 3.9952) <= 5e-5, f"(n4) ratio {m_pi / m_2pi:.6f}"
    assert abs(small - 1.0) <= 1e-5, f"(n4) Gamma_sph(1e-3) = {small}"
    assert abs(zero) <= 1e-15, f"(n4) Gamma at Z_s = eta0 f/g is {abs(zero):.2e}"
    assert same <= 1e-15, f"(n4) Gamma(z = 1) differs from Gamma_sph by {same:.2e}"
    assert abs(pec - pec_cf) <= 1e-15, f"(n4) PEC reflection off by {abs(pec - pec_cf):.2e}"
    assert abs(flipped - 19.7645) <= 5e-5, f"(n4) |Gamma| at gamma -> -gamma = {flipped:.6f}"
    assert abs(flipped - np.sqrt(1.0 + 4.0 * x**4)) <= 1e-12


def test_n5_planar_reflection_over_a_cube():
    n_u, n_p = 2000, 4000
    u = -1.0 + (np.arange(n_u) + 0.5) * (2.0 / n_u)  # midpoint in cosθ_z
    p = (np.arange(n_p) + 0.5) * (2.0 * np.pi / n_p)
    uu, pp = np.meshgrid(u, p, indexing="ij")
    ss = np.sqrt(1.0 - uu**2)
    omega = np.stack([ss * np.cos(pp), ss * np.sin(pp), uu])
    cos_inc = np.max(np.abs(omega), axis=0)
    refl = abc_first_order_reflection_te(cos_inc)
    weight = 1.0 - uu**2  # sin²θ_z, the dipole's power pattern
    rho = float(np.sqrt(np.sum(refl**2 * weight) / np.sum(weight)))
    corner = float(abs(abc_first_order_reflection_te(1.0 / np.sqrt(3.0))))
    grid_max = float(np.max(np.abs(refl)))
    _say(
        f"[TH-5 1b-0] (n5) dipole-weighted rms planar R over a cube: rho = {rho:.5f} "
        f"(predicted 0.113); max |R| = {corner:.5f} at the corner direction "
        f"(grid max {grid_max:.5f})"
    )
    assert abs(rho - 0.113) <= 1e-3, f"(n5) rho = {rho:.5f}"
    assert abs(corner - 0.2679) <= 1e-4, f"(n5) corner |R| = {corner:.5f}"
    assert grid_max <= corner + 1e-12


# =====================================================================
# `TH-5` step 1b (§9 item 114): the dipole gate on the sphere, FEM.
#
# Ball r <= R_BALL = 0.5 m (λ₀ = 1 m, so k₀R = π), degree 1.  The drive pins
# every dof of the cells with r ≤ R_DRIVE to E_out (r clamped ≥ R_CLAMP); the
# shell outside carries A E_out + B conj(E_out) and the termination on facet
# tag 1 sets Γ = B/A, read by a volume least-squares projection over
# 0.22 < r < 0.42 (χ DG0 by midpoints, quadrature degree 8) — never by point
# sampling (step 1a's finding).
#
#   (E) Z_s = η₀ f/g: the exact l = 1 impedance, Γ = 0 (premise);
#   (A) the first-order ABC, Γ = Γ_sph(π) = −0.047971 + 0.016085j.
#
# W1 runs with ``BAND = None`` and prints the candidate B = 3|Γ_E(λ₀/40)|
# rounded up to one figure; B is frozen only if B ≤ |Γ_sph(π)|/3 = 1.69e-2.
# W2 asserts |Γ_E| ≤ B and |Γ_A − Γ_sph| ≤ B at λ₀/40, and that PEC misses
# Γ_sph by > B.  The constant is the only edit between the windows.
# =====================================================================

import time  # noqa: E402

import pytest  # noqa: E402
import ufl  # noqa: E402
from dolfinx import fem  # noqa: E402
from dolfinx import mesh as dmesh  # noqa: E402
from dolfinx.fem.petsc import LinearProblem  # noqa: E402
from petsc4py import PETSc  # noqa: E402

from fem_em_solver.core.boundary_terms import (  # noqa: E402
    first_order_abc_term,
    impedance_boundary_term,
)
from fem_em_solver.io.mesh import MeshGenerator  # noqa: E402
from fem_em_solver.utils.constants import C_0  # noqa: E402
from tests.validation.test_th5_abc_reflection import (  # noqa: E402
    _assemble_complex,
    complex_only,
)

R_BALL = 0.5
R_DRIVE = 0.15
R_CLAMP = 0.05
PROJ_R = (0.22, 0.42)
QDEG_1B = {"quadrature_degree": 8}
PLANTED_1B = 0.01 - 0.02j
BAND_CEILING = abs(_gamma_sphere(np.pi)) / 3.0  # 1.69e-2: a band above this is no gate
BAND = 6e-3  # frozen from W1 (`20260930T184708Z_TH-5-step1b-W1.log:125`): 3|Γ_E(λ₀/40)| = 5.31e-3, rounded up


def _e_out_np(k):
    """Numpy drive for ``interpolate``: E_out at x (3, N), r clamped ≥ R_CLAMP.

    Uses sinθ θ̂ = cosθ r̂ − ẑ, so nothing divides by ρ.
    """

    def ev(x):
        p = np.asarray(x[:3], dtype=np.float64).T.copy()
        r = np.linalg.norm(p, axis=1)
        zero = r == 0.0
        p[zero] = np.array([0.0, 0.0, R_CLAMP])
        r[zero] = R_CLAMP
        small = r < R_CLAMP
        p[small] *= (R_CLAMP / r[small])[:, None]
        r[small] = R_CLAMP
        r_hat = p / r[:, None]
        cos_t = r_hat[:, 2]
        s_theta_hat = cos_t[:, None] * r_hat - np.array([0.0, 0.0, 1.0])
        xx = k * r
        ph = np.exp(-1j * xx)
        f = 1.0 - 1j / xx - 1.0 / xx**2
        e = (2.0 * cos_t * (1.0 / xx**2 - 1j / xx**3) * ph)[:, None] * r_hat + (
            1j * f / xx * ph
        )[:, None] * s_theta_hat
        return e.T.astype(PETSc.ScalarType)

    return ev


def _e_out_ufl(msh, k):
    """E_out in UFL (complex via a Constant j), the same sinθ θ̂ identity."""
    x = ufl.SpatialCoordinate(msh)
    j = fem.Constant(msh, PETSc.ScalarType(1j))
    r = ufl.sqrt(ufl.inner(x, x))
    xx = k * r
    ph = ufl.cos(xx) - j * ufl.sin(xx)
    r_hat = x / r
    cos_t = x[2] / r
    ez = ufl.as_vector((0.0, 0.0, 1.0))
    f = 1.0 - j / xx - 1.0 / xx**2
    return ph * (
        2.0 * cos_t * (1.0 / xx**2 - j / xx**3) * r_hat + (j * f / xx) * (cos_t * r_hat - ez)
    )


class _DipoleProjection:
    """Volume LSQ of a vector field onto {E_out, conj(E_out)} over PROJ_R.

    As ``_Projection`` (`test_th5_abc_reflection.py:190-240`), vector-valued:
    M[p][q] = ∫χ inner(φ_q, φ_p), b[p] = ∫χ inner(E, φ_p); Γ = B/A.
    """

    def __init__(self, e, msh, k):
        self.comm = msh.comm
        tdim = msh.topology.dim
        imap = msh.topology.index_map(tdim)
        n = imap.size_local + imap.num_ghosts
        mids = dmesh.compute_midpoints(msh, tdim, np.arange(n, dtype=np.int32))
        rr = np.linalg.norm(mids, axis=1)
        chi = fem.Function(fem.functionspace(msh, ("DG", 0)))
        chi.x.array[:] = ((rr > PROJ_R[0]) & (rr < PROJ_R[1])).astype(PETSc.ScalarType)
        chi.x.scatter_forward()
        phi0 = _e_out_ufl(msh, k)
        phi = (phi0, ufl.conj(phi0))
        dx = ufl.Measure("dx", domain=msh, metadata=QDEG_1B)
        self._m = [
            [fem.form(chi * ufl.inner(phi[q], phi[p]) * dx) for q in range(2)] for p in range(2)
        ]
        self._b = [fem.form(chi * ufl.inner(e, phi[p]) * dx) for p in range(2)]

    def gamma(self):
        M = np.array(
            [[_assemble_complex(self._m[p][q], self.comm) for q in range(2)] for p in range(2)]
        )
        b = np.array([_assemble_complex(self._b[p], self.comm) for p in range(2)])
        a_, b_ = np.linalg.solve(M, b)
        return complex(b_ / a_)


def _ball(comm, res):
    msh, _cells, ft = MeshGenerator.ball_domain(radius=R_BALL, resolution=res, comm=comm)
    n = comm.allreduce(msh.topology.index_map(3).size_local, op=MPI.SUM)
    return msh, ft, n


def _solve_ball(msh, ft, k, termination, degree=1):
    """termination: 'exact' | 'abc' | 'pec' | 'conj'."""
    omega = k * C_0
    V = fem.functionspace(msh, ("N1curl", degree))
    u, v = ufl.TrialFunction(V), ufl.TestFunction(V)
    a = ufl.inner(ufl.curl(u), ufl.curl(v)) * ufl.dx - k * k * ufl.inner(u, v) * ufl.dx
    x_r = k * R_BALL
    if termination == "exact":
        z_s = ETA_0 * _f(x_r) / _g(x_r)
        a = a + impedance_boundary_term(msh, ft, 1, z_s, omega)(u, v)
    elif termination == "abc":
        a = a + first_order_abc_term(msh, ft, 1, omega)(u, v)
    elif termination == "conj":
        a = a + first_order_abc_term(msh, ft, 1, omega, cos_theta=-1.0)(u, v)
    zero_vec = fem.Constant(msh, np.zeros(3, dtype=PETSc.ScalarType))
    L = ufl.inner(zero_vec, v) * ufl.dx

    tdim = msh.topology.dim
    msh.topology.create_connectivity(tdim, tdim)
    drive_cells = dmesh.locate_entities(
        msh, tdim, lambda x: x[0] ** 2 + x[1] ** 2 + x[2] ** 2 <= R_DRIVE**2
    )
    g = fem.Function(V)
    g.interpolate(_e_out_np(k))
    bcs = [fem.dirichletbc(g, fem.locate_dofs_topological(V, tdim, drive_cells))]
    if termination == "pec":
        zero = fem.Function(V)
        zero.x.array[:] = 0.0
        bcs.append(fem.dirichletbc(zero, fem.locate_dofs_topological(V, tdim - 1, ft.find(1))))
    problem = LinearProblem(
        a,
        L,
        bcs=bcs,
        petsc_options={
            "ksp_type": "preonly",
            "pc_type": "lu",
            "pc_factor_mat_solver_type": "mumps",
        },
        petsc_options_prefix=f"fem_em_th5_1b_{termination}_{degree}_",
    )
    e_h = problem.solve()
    e_h.x.scatter_forward()
    return e_h


def _round_up_one_figure(v):
    e = np.floor(np.log10(v))
    return float(np.ceil(v / 10.0**e) * 10.0**e)


@complex_only
@pytest.mark.integration
def test_1b_projection_premise_planted():
    """(P), no solve: a planted φ₀ + Γ_p φ₁ is recovered to 1e-10 at res 0.05."""
    comm = MPI.COMM_WORLD
    msh, _ft, n = _ball(comm, 0.05)
    phi0 = _e_out_ufl(msh, K0)
    planted = phi0 + fem.Constant(msh, PETSc.ScalarType(PLANTED_1B)) * ufl.conj(phi0)
    g = _DipoleProjection(planted, msh, K0).gamma()
    _say(
        f"\n[TH-5 1b (P)] res 0.05 ({n} cells): planted {PLANTED_1B}, recovered "
        f"{g:.12e}, |diff| = {abs(g - PLANTED_1B):.3e}"
    )
    assert abs(g - PLANTED_1B) <= 1e-10, g


@complex_only
@pytest.mark.integration
def test_1b_dipole_sphere_gate():
    """(E) Γ at Z_s = η₀f/g and (A) Γ − Γ_sph(π) under the first-order ABC."""
    comm = MPI.COMM_WORLD
    t0 = time.perf_counter()
    g_sph = _gamma_sphere(np.pi)
    out = {}
    for res in (0.05, 0.025):
        t = time.perf_counter()
        msh, ft, n = _ball(comm, res)
        t_mesh = time.perf_counter() - t
        terms = ("exact", "abc", "pec", "conj") if res == 0.05 else ("exact", "abc")
        for term in terms:
            t = time.perf_counter()
            e_h = _solve_ball(msh, ft, K0, term)
            out[(res, term)] = _DipoleProjection(e_h, msh, K0).gamma()
            _say(
                f"\n[TH-5 1b] res {res} (lambda0/{round(1 / res)}, {n} cells, mesh {t_mesh:.1f} s) "
                f"{term}: Gamma = {out[(res, term)]:.6e} ({time.perf_counter() - t:.1f} s)"
            )
    for res in (0.05, 0.025):
        e_ = abs(out[(res, "exact")])
        a_ = abs(out[(res, "abc")] - g_sph)
        _say(
            f"[TH-5 1b] res {res}: |Gamma_E| = {e_:.4e}, |Gamma_A - Gamma_sph| = {a_:.4e} "
            f"(predicted 3.8e-3 / 9.5e-4 at lambda0/20 / /40, x/ 3; faceting adds 4.6e-4 / 1.2e-4)"
        )
    cand = _round_up_one_figure(3.0 * abs(out[(0.025, "exact")]))
    _say(
        f"[TH-5 1b] band candidate B = 3|Gamma_E(lambda0/40)| rounded up = {cand:.0e} "
        f"(ceiling |Gamma_sph|/3 = {BAND_CEILING:.4e}; predicted 3e-3); "
        f"{'FREEZABLE' if cand <= BAND_CEILING else 'ABOVE THE CEILING'}; BAND = {BAND}"
    )

    pec, conj, abc05 = out[(0.05, "pec")], out[(0.05, "conj")], out[(0.05, "abc")]
    _say(
        f"[TH-5 1b control] PEC @ res 0.05: Gamma = {pec:.6f} (predicted 0.777 - 0.629j), "
        f"|Gamma_PEC|/|Gamma_A| = {abs(pec) / abs(abc05):.2f} (predicted 19.8, 18-22, ceiling 19.76), "
        f"miss ratio |Gamma_PEC - Gamma_sph|/|Gamma_A - Gamma_sph| = "
        f"{abs(pec - g_sph) / abs(abc05 - g_sph):.1f} (predicted 275, 90-830)"
    )
    _say(
        f"[TH-5 1b control] conjugated gamma @ res 0.05: |Gamma| = {abs(conj):.4f} "
        f"(predicted 19.76; a source, never asserted)"
    )

    # Printed only: k₀R = 2π on res 0.025 at λ₀ = 0.5 m; degree 2 at res 0.1.
    msh, ft, n = _ball(comm, 0.025)
    k2 = 2.0 * K0
    g2 = _DipoleProjection(_solve_ball(msh, ft, k2, "abc"), msh, k2).gamma()
    _say(
        f"[TH-5 1b printed] k0R = 2 pi (lambda0 = 0.5 m, res 0.025 = lambda/20, {n} cells): "
        f"|Gamma_A| = {abs(g2):.6f} vs |Gamma_sph(2 pi)| = {abs(_gamma_sphere(2 * np.pi)):.6f}; "
        f"|Gamma_A - Gamma_sph| = {abs(g2 - _gamma_sphere(2 * np.pi)):.4e}"
    )
    msh, ft, n = _ball(comm, 0.1)
    for term in ("exact", "abc"):
        gd2 = _DipoleProjection(_solve_ball(msh, ft, K0, term, degree=2), msh, K0).gamma()
        ref = 0.0 if term == "exact" else g_sph
        _say(
            f"[TH-5 1b printed] degree 2, res 0.1 ({n} cells) {term}: |Gamma - ref| = "
            f"{abs(gd2 - ref):.4e}"
        )
    _say(
        f"[TH-5 1b] box leg not run (printed-only; unpriced in-slot); elapsed "
        f"{time.perf_counter() - t0:.1f} s"
    )

    if BAND is not None:
        assert BAND <= BAND_CEILING, BAND
        assert abs(out[(0.025, "exact")]) <= BAND, out[(0.025, "exact")]
        assert abs(out[(0.025, "abc")] - g_sph) <= BAND, out[(0.025, "abc")]
        assert abs(pec - g_sph) > BAND, pec
