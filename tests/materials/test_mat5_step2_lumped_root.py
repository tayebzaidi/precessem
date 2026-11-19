"""`MAT-5` step 2 — the lumped root of the EM-coupled σ(T) fixed point (§9 item 105).

Pure numpy/scipy.  As Bi → 0 the ball is isothermal, σ = σ₀[1 + αθ] is
uniform, and the steady surface rise θ_L solves

    θ = (R/3h) · G · ρ · SAR_cf(σ₀[1 + αθ]),

``SAR_cf`` being `MAT-4`'s quasi-static lossy-sphere closed form
(``post/sar.py: uniform_sphere_sar_closed_form``).  Fixture: σ₀ = 0.285 S/m
(``SIGMA_HIGH``/2), α = 0.02 /K, h = 3 W/m²/K (Bi = hR/k = 0.05), and G chosen
so θ_L = θ* = 50 K, i.e. σ* = 0.57 S/m — `MAT-4`'s gated σ.

Model bound: to first order the absorbed power sees the volume-mean σ, and
θ̄/θ_s = 1 + Bi/5, θ(0)/θ(R) = 1 + Bi/2, so the true θ_s solves the same
equation with α(1 + c), c = Bi/5 ≤ Bi/2.

Anchors (asserted): (N1) the brentq root is 50 K to 1e-12; (N2) the roots at
c = Bi/5 and Bi/2 read −0.2302 % and −0.5727 % off θ_L (to one unit in the
last printed digit — the tolerance is this file's, the readings are the
item's), the latter inside 0.6 %; (N3) a *fixed*-|E| source G·q_cf(σ₀)[1 + αθ]
has no steady state on this fixture: ``robin_sphere_linear_source_coefficients``
raises (mR 0.4331 > λ₁ 0.3854).  Λ and the Picard factor L are printed records.

The fixture constants and the two root functions are imported by name by
``tests/validation/test_mat5_em_coupled_sphere.py``.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.optimize import brentq

from fem_em_solver.post.sar import uniform_sphere_sar_closed_form
from fem_em_solver.utils.analytical import (
    robin_sphere_eigenvalues,
    robin_sphere_linear_source_coefficients,
)
from fem_em_solver.utils.constants import EPSILON_0
from tests.validation.test_lossy_sphere_sar import (
    E0,
    EPSILON_R_SPHERE,
    OMEGA,
    RHO_KG_M3,
    SIGMA_HIGH,
    SPHERE_RADIUS,
)

R = SPHERE_RADIUS
K = 0.6
H_STEP2 = 3.0  # Bi = hR/k = 0.05
BIOT = H_STEP2 * R / K
ALPHA = 0.02  # /K
SIGMA_0 = SIGMA_HIGH / 2.0  # 0.285 S/m
SIGMA_STAR = SIGMA_HIGH  # 0.57 S/m, MAT-4's gated sigma
THETA_STAR = 50.0  # K
T_REF = 20.0  # T0 = T_inf


def q_cf_unit(sigma: float) -> float:
    """ρ·SAR_cf(σ) at E₀ = 1 V/m [W/m³ per (V/m)²] — multiply by G."""
    cf = uniform_sphere_sar_closed_form(
        e0=E0,
        epsilon_r=EPSILON_R_SPHERE,
        sigma=float(sigma),
        omega=OMEGA,
        rho=RHO_KG_M3,
        epsilon_0=EPSILON_0,
    )
    return RHO_KG_M3 * cf["sar_w_per_kg"]


G_STEP2 = 3.0 * H_STEP2 * THETA_STAR / (R * q_cf_unit(SIGMA_STAR))  # 5.6191e8
Q_STAR = G_STEP2 * q_cf_unit(SIGMA_STAR)  # 45 000 W/m^3


def lumped_root(c: float = 0.0) -> float:
    """brentq root of θ = (R/3h)·G·ρ·SAR_cf(σ₀[1 + α(1 + c)θ])."""

    def residual(theta):
        sigma = SIGMA_0 * (1.0 + ALPHA * (1.0 + c) * theta)
        return theta - R / (3.0 * H_STEP2) * G_STEP2 * q_cf_unit(sigma)

    return float(brentq(residual, 0.0, 200.0, xtol=1e-14, rtol=4.0 * np.finfo(float).eps))


def feedback_lambda(sigma: float) -> float:
    """Λ = d ln SAR / d ln σ = 1 − 2t²/((εᵣ + 2)² + t²), t = σ/(ωε₀)."""
    t = sigma / (OMEGA * EPSILON_0)
    return 1.0 - 2.0 * t * t / ((EPSILON_R_SPHERE + 2.0) ** 2 + t * t)


def test_n1_lumped_root_is_50_k():
    theta_l = lumped_root(0.0)
    lam = feedback_lambda(SIGMA_STAR)
    a_th = ALPHA * theta_l
    picard_l = lam * a_th / (1.0 + a_th)
    sigma_zero = (EPSILON_R_SPHERE + 2.0) * OMEGA * EPSILON_0
    print(
        f"\nMAT5s2 N1 G = {G_STEP2:.6e} (predicted 5.6191e8)  q* = {Q_STAR:.6f} W/m^3 "
        f"(predicted 45000)  Bi = {BIOT:.6f}"
    )
    print(
        f"MAT5s2 N1 theta_L = {theta_l:.15f} K  rel = {abs(theta_l - THETA_STAR) / THETA_STAR:.3e} "
        f"(bound 1e-12)  sigma* = {SIGMA_0 * (1.0 + a_th):.15f}"
    )
    print(
        f"MAT5s2 record Lambda(0.57) = {lam:.6f} (predicted -0.6004)  "
        f"Lambda zero at sigma = {sigma_zero:.6f} S/m (predicted 0.2848)  "
        f"Picard L = {picard_l:.6f} (predicted -0.3002)"
    )
    assert abs(theta_l - THETA_STAR) / THETA_STAR <= 1e-12
    assert abs(Q_STAR - 3.0 * H_STEP2 * THETA_STAR / R) <= 1e-9 * Q_STAR


def test_n2_model_bracket_roots():
    theta_l = lumped_root(0.0)
    shift = {}
    for label, c in (("Bi/5", BIOT / 5.0), ("Bi/2", BIOT / 2.0)):
        root = lumped_root(c)
        shift[label] = 100.0 * (root / theta_l - 1.0)
        lam = feedback_lambda(SIGMA_STAR)
        picard_l = lam * ALPHA * theta_l / (1.0 + ALPHA * theta_l)
        print(
            f"MAT5s2 N2 c = {label} = {c:.4f}: root = {root:.6f} K  shift = {shift[label]:+.4f} %  "
            f"(first order L/(1-L)*c = {100.0 * picard_l / (1.0 - picard_l) * c:+.4f} %)"
        )
    print("MAT5s2 N2 predicted -0.2302 % / -0.5727 %; c = Bi/5 root predicted 49.885 K")
    # One unit in the last digit the item prints.
    assert abs(shift["Bi/5"] - (-0.2302)) <= 1.0e-4, shift
    assert abs(shift["Bi/2"] - (-0.5727)) <= 1.0e-4, shift
    assert abs(shift["Bi/2"]) <= 0.6, shift


def test_n3_fixed_field_has_no_steady_state():
    q0 = G_STEP2 * q_cf_unit(SIGMA_0)
    m_r = float(np.sqrt(ALPHA * q0 / K) * R)
    lam1 = float(robin_sphere_eigenvalues(BIOT, 1)[0])
    print(
        f"MAT5s2 N3 q0 = G*q_cf(sigma0) = {q0:.4f} W/m^3  mR = {m_r:.4f} (predicted 0.4331)  "
        f"lambda_1 = {lam1:.4f} (predicted 0.3854)"
    )
    assert m_r > lam1
    with pytest.raises(ValueError, match="thermal runaway"):
        robin_sphere_linear_source_coefficients(radius=R, k=K, q0=q0, alpha=ALPHA, h=H_STEP2)
