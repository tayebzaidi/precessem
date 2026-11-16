"""`MAT-5` step 1, numpy anchors N1–N6 (§9 item 85).

σ(T) = σ₀[1 + α(T − T₀)] and the fixed-|E| Joule-heated Robin ball on
`MAT-9` step 1's fixture (imported by name): with Q = q₀[1 + αθ],
θ(r) = (1/α)[Bi·R·sin(mr)/(r·D) − 1], D = mR cos mR + (Bi − 1) sin mR,
m² = αq₀/k.  α = 0.02 /K is of the order of saline's (not a material claim).

Predicted (printed): D = 1.04495; θ(0) = 15.114 K, θ(R) = 8.117 K,
θ̄ = 10.875 K, P = 10.200 W (`MAT-9`: 12.222 / 6.667 / 8.889 K, 8.378 W);
λ₁ ≈ 1.9064, α_crit ≈ 0.1090 /K.
"""

from __future__ import annotations

import numpy as np
import pytest

from fem_em_solver.materials import LinearTemperatureConductivity
from fem_em_solver.utils.analytical import (
    robin_sphere_eigenvalues,
    robin_sphere_linear_source_coefficients,
    robin_sphere_linear_source_temperature,
    robin_sphere_steady_temperature,
)
from tests.validation.test_mat9_steady_conduction_sphere import (  # MAT-9 fixture
    H,
    K,
    Q,
    R,
    T_INF,
)

Q0 = Q
ALPHA = 0.02
T0 = T_INF
BI = H * R / K
V = 4.0 / 3.0 * np.pi * R**3


def _theta(r, alpha=ALPHA):
    return robin_sphere_linear_source_temperature(
        r, radius=R, k=K, q0=Q0, alpha=alpha, h=H
    )


def _coeffs(alpha=ALPHA):
    return robin_sphere_linear_source_coefficients(
        radius=R, k=K, q0=Q0, alpha=alpha, h=H
    )


def _dtheta(r, alpha=ALPHA):
    """Analytic θ′ = A (mr cos mr − sin mr)/r²."""
    m, _, A, _ = _coeffs(alpha)
    r = np.asarray(r, dtype=float)
    return A * (m * r * np.cos(m * r) - np.sin(m * r)) / r**2


def _d2theta(r, alpha=ALPHA):
    """Analytic θ″ = A (−m²r² sin mr − 2mr cos mr + 2 sin mr)/r³."""
    m, _, A, _ = _coeffs(alpha)
    r = np.asarray(r, dtype=float)
    return A * (
        -(m * r) ** 2 * np.sin(m * r) - 2.0 * m * r * np.cos(m * r)
        + 2.0 * np.sin(m * r)
    ) / r**3


def _theta_bar(alpha=ALPHA):
    """Closed-form volume mean: (3/R³)∫θ r² dr, ∫₀ᴿ r sin(mr) dr = (sin mR − mR cos mR)/m²."""
    m, D, _, _ = _coeffs(alpha)
    mR = m * R
    return (1.0 / alpha) * (
        3.0 * BI * (np.sin(mR) - mR * np.cos(mR)) / (mR**2 * D) - 1.0
    )


def test_n1_linear_conductivity():
    s0 = 0.57
    c = LinearTemperatureConductivity(sigma_ref=s0, alpha_per_k=ALPHA, t_ref=T0)
    assert c.sigma(T0) == s0  # exactly
    xs = np.array([-0.9, -0.5, 0.1, 0.25, 1.0, 3.0])
    got = c.sigma(T0 + xs / ALPHA)
    rel = np.max(np.abs(got / (s0 * (1.0 + xs)) - 1.0))
    print(f"MAT5 N1 max rel |σ(T0 + x/α) − σ0(1 + x)| = {rel:.3e}")
    assert rel <= 1e-15
    with pytest.raises(ValueError):
        c.sigma(T0 - 1.0 / ALPHA)  # factor exactly 0
    with pytest.raises(ValueError):
        c.sigma(np.array([T0, T0 - 2.0 / ALPHA]))  # one non-positive entry
    with pytest.raises(ValueError):
        LinearTemperatureConductivity(sigma_ref=0.0, alpha_per_k=ALPHA, t_ref=T0)
    # UFL factor is the bare polynomial; on floats it equals the numpy factor.
    assert c.factor_ufl(T0 + 5.0) == pytest.approx(1.0 + 5.0 * ALPHA, rel=1e-15)


def test_n2_small_alpha_reduces_to_mat9():
    rs = np.array([0.0, R / 2, R])
    th = _theta(rs, alpha=1e-6)
    th9 = robin_sphere_steady_temperature(rs, radius=R, k=K, q=Q0, h=H)
    dev = np.max(np.abs(th / th9 - 1.0))
    print(f"MAT5 N2 α=1e-6 max |θ_α/θ_MAT9 − 1| = {dev:.3e} (bound 1e-4)")
    assert dev <= 1e-4


def test_n3_robin_condition():
    lhs = -K * _dtheta(R)
    rhs = H * _theta(R)
    rel = abs(lhs - rhs) / abs(rhs)
    print(f"MAT5 N3 |−kθ′(R) − hθ(R)|/hθ(R) = {rel:.3e}")
    assert rel <= 1e-12


def test_n4_global_balance():
    tb = _theta_bar()
    conv = H * 4.0 * np.pi * R**2 * float(_theta(R))
    src = Q0 * V * (1.0 + ALPHA * tb)
    rel = abs(conv - src) / src
    print(f"MAT5 N4 |h4πR²θ(R) − q0V(1+αθ̄)|/P = {rel:.3e}; P = {src:.4f} W")
    assert rel <= 1e-12


def test_n5_eigen_equation():
    m, D, _, lam1 = _coeffs()
    lam = robin_sphere_eigenvalues(BI, 1)[0]
    d_at = lam * np.cos(lam) + (BI - 1.0) * np.sin(lam)
    print(f"MAT5 N5 |D(λ1/R)| = {abs(d_at):.3e}, λ1 = {lam:.6f}")
    assert abs(d_at) <= 1e-12
    # The guard fires at/after runaway.
    a_crit = lam**2 * K / (Q0 * R**2)
    with pytest.raises(ValueError):
        _theta(0.0, alpha=a_crit * 1.0001)


def test_n6_radial_ode_residual():
    m, _, _, _ = _coeffs()
    rs = np.array([0.2, 0.5, 0.8]) * R
    res = _d2theta(rs) + 2.0 / rs * _dtheta(rs) + m**2 * _theta(rs) + Q0 / K
    worst = float(np.max(np.abs(res)) / (Q0 / K))
    print(f"MAT5 N6 max |θ″ + 2θ′/r + m²θ + q0/k|/(q0/k) = {worst:.3e}")
    assert worst <= 1e-10


def test_predicted_values_printed():
    m, D, _, lam1 = _coeffs()
    th0, thR = float(_theta(0.0)), float(_theta(R))
    tb = _theta_bar()
    P = Q0 * V * (1.0 + ALPHA * tb)
    a_crit = lam1**2 * K / (Q0 * R**2)
    kappa = (m * R / lam1) ** 2
    print(
        f"MAT5 predicted: m²={m**2:.4f} mR={m * R:.5f} D={D:.5f} "
        f"θ(0)={th0:.4f} K θ(R)={thR:.4f} K θ̄={tb:.4f} K P={P:.4f} W "
        f"λ1={lam1:.5f} α_crit={a_crit:.5f} /K κ=(mR/λ1)²={kappa:.5f}"
    )
    print(
        "MAT5 predicted (item 85): D=1.04495 θ(0)=15.114 θ(R)=8.117 "
        "θ̄=10.875 P=10.200 λ1≈1.9064 α_crit≈0.1090 κ≈0.1834"
    )
    # Consistency with MAT-9's quoted fixture (the α → 0 figures).
    assert float(robin_sphere_steady_temperature(0.0, radius=R, k=K, q=Q0, h=H)) \
        == pytest.approx(12.2222222, rel=1e-6)
