"""`MAT-3` step 1b: `TH-6`'s lossy plane wave in a two-pole Cole–Cole medium.

Step 1 (tests/validation/test_mat3_debye_plane_wave.py) ran one Debye pole
(α = 0).  This module runs a broadened, two-pole medium — the shape Gabriel
tissues use — so the evaluator's (jωτ)^(1−α) branch is exercised by a solve:

    ε̂(f) = ε∞ + Σₙ Δεₙ/(1 + (j2πfτₙ)^(1−αₙ)) − jσ/(2πfε₀),
    ε∞ = 4, (Δε, τ, α) = (50, 1/(2π·100 MHz), 0.30), (4000, 1/(2π·1 MHz), 0.35),
    σ = 0.2 S/m

(synthetic, muscle-like: a γ-like pole near 100 MHz, a large β-like pole near
1 MHz).  The fixture, the frozen-frequency route (``as_homogeneous(f)``) and
the comparand ``k = k₀√ε̂(f)`` are step 1's.

Anchors (`TH-6`'s imported 5 % band): (Ev) the evaluator equals an
independent ``cmath`` sum to 1e-12; (i) 127.74 MHz, 12³ → 24³, rel L2 < 5 %,
rate > 0.9, α and β within 5 % at 24³; (ii) 63.87 MHz, 24³, α and β within
5 %.  Negative controls at 16³: (a′) the same poles with every αₙ = 0
(Debye-reduced, closed-form ceiling ≈ 19.9 %), (b′) ε̂(63.87 MHz) solved at
127.74 MHz (ceiling ≈ 29.8 %): each must miss the Cole–Cole α at 127.74 MHz
by more than the band.
"""

from __future__ import annotations

import cmath
import math

import numpy as np
from mpi4py import MPI

from fem_em_solver.materials import ColeColeMaterial
from fem_em_solver.utils.constants import EPSILON_0, MU_0
from tests.complex_mode import complex_only
from tests.validation.test_lossy_plane_wave import (
    MU_R,
    _analytic_alpha_beta,
    _solve_plane_wave,
    _wavenumber,
)
from tests.validation.test_mat3_debye_plane_wave import BAND, F_HIGH, F_LOW

EPS_INF = 4.0
POLES = (
    (50.0, 1.0 / (2.0 * math.pi * 100e6), 0.30),
    (4000.0, 1.0 / (2.0 * math.pi * 1e6), 0.35),
)
SIGMA_IONIC = 0.2

COLE_COLE = ColeColeMaterial(epsilon_inf=EPS_INF, poles=POLES, sigma_ionic=SIGMA_IONIC)
DEBYE_REDUCED = ColeColeMaterial.debye(
    EPS_INF, poles=tuple((d, t) for d, t, _ in POLES), sigma_ionic=SIGMA_IONIC
)


def _eps_cmath(f: float, poles=POLES) -> complex:
    """Independent standard-library ε̂(f): principal-branch ``(1j*w*tau)**(1-alpha)``."""
    w = 2.0 * math.pi * f
    eps = complex(EPS_INF)
    for delta, tau, alpha in poles:
        eps += delta / (1.0 + (1j * w * tau) ** (1.0 - alpha))
    return eps - 1j * SIGMA_IONIC / (w * EPSILON_0)


def _truth(f: float, material: ColeColeMaterial) -> tuple[float, float, complex]:
    """(α, β, ε̂) straight from k₀√ε̂(f), cross-checked against TH-6's helpers."""
    eps = material.complex_permittivity(f)
    k0 = 2.0 * np.pi * f * np.sqrt(MU_0 * EPSILON_0)
    k = k0 * np.sqrt(eps)
    assert k.imag < 0.0
    hm = material.as_homogeneous(f)
    k_th6 = _wavenumber(hm.sigma, frequency_hz=f, epsilon_r=hm.epsilon_r)
    a_th6, b_th6 = _analytic_alpha_beta(hm.sigma, frequency_hz=f, epsilon_r=hm.epsilon_r)
    assert abs(k_th6 - k) <= 1e-12 * abs(k)
    assert abs(a_th6 + k.imag) <= 1e-10 * abs(k)
    assert abs(b_th6 - k.real) <= 1e-10 * abs(k)
    return float(-k.imag), float(k.real), complex(eps)


def _wrong_f_truth() -> tuple[float, float]:
    """Closed-form (α, β) of ε̂(63.87 MHz) frozen and solved at 127.74 MHz."""
    hm = COLE_COLE.as_homogeneous(F_LOW)
    a, b = _analytic_alpha_beta(hm.sigma, frequency_hz=F_HIGH, epsilon_r=hm.epsilon_r)
    return float(a), float(b)


def _solve(n: int, f: float, material: ColeColeMaterial, material_f: float | None = None):
    """Solve at ``f`` with ``material`` frozen at ``material_f`` (default f)."""
    hm = material.as_homogeneous(f if material_f is None else material_f)
    assert hm.mu_r == MU_R
    return _solve_plane_wave(n, hm.sigma, frequency_hz=f, epsilon_r=hm.epsilon_r)


@complex_only
def test_cole_cole_evaluator_matches_cmath_and_controls_can_miss():
    """(Ev) + the table and the two closed-form control ceilings; no mesh."""
    out = [
        "\n[MAT-3 step 1b] Cole-Cole eps_inf=4, (50, 1/(2pi 100 MHz), 0.30), "
        "(4000, 1/(2pi 1 MHz), 0.35), sigma=0.2"
    ]
    preds = {
        F_LOW: "181.505 - 283.983j / 11.8044 / 21.5541",
        F_HIGH: "118.420 - 182.618j / 18.8581 / 34.7047",
    }
    for f in (F_LOW, F_HIGH):
        ev = COLE_COLE.complex_permittivity(f)
        ref = _eps_cmath(f)
        rel = abs(ev - ref) / abs(ref)
        a, b, _ = _truth(f, COLE_COLE)
        out.append(
            f"  {f / 1e6:.2f} MHz: eps_hat={ev.real:.3f} {ev.imag:+.3f}j "
            f"alpha={a:.4f} beta={b:.4f} (pred {preds[f]}); |ev-cmath|/|cmath|={rel:.3e}"
        )
        assert rel <= 1e-12, f"evaluator vs cmath at {f:.4e} Hz: {rel:.3e}"
        # Debye-reduced evaluator against the same cmath sum with alpha = 0.
        ref_d = _eps_cmath(f, poles=tuple((d, t, 0.0) for d, t, _ in POLES))
        rel_d = abs(DEBYE_REDUCED.complex_permittivity(f) - ref_d) / abs(ref_d)
        assert rel_d <= 1e-12, f"Debye-reduced evaluator vs cmath: {rel_d:.3e}"
    # Branch sanity: (jωτ)^(1−α) has arg π(1−α)/2 on the principal branch.
    for _, tau, alpha in POLES:
        z = (1j * 2.0 * math.pi * F_HIGH * tau) ** (1.0 - alpha)
        assert abs(cmath.phase(z) - 0.5 * math.pi * (1.0 - alpha)) <= 1e-12

    alpha_hi, _, _ = _truth(F_HIGH, COLE_COLE)
    alpha_lo, _, _ = _truth(F_LOW, COLE_COLE)
    a_dr, b_dr, _ = _truth(F_HIGH, DEBYE_REDUCED)
    a_dr_lo, _, _ = _truth(F_LOW, DEBYE_REDUCED)
    a_wf, b_wf = _wrong_f_truth()
    ceil_dr = abs(a_dr - alpha_hi) / alpha_hi
    ceil_dr_lo = abs(a_dr_lo - alpha_lo) / alpha_lo
    ceil_wf = abs(a_wf - alpha_hi) / alpha_hi
    out.append(
        f"  ceiling (a') Debye-reduced alpha miss @127.74: {ceil_dr:.4%} (pred 19.91%); "
        f"@63.87: {ceil_dr_lo:.4%} (pred 17.15%)"
    )
    out.append(f"  ceiling (b') wrong-f alpha miss @127.74: {ceil_wf:.4%} (pred 29.77%)")
    if MPI.COMM_WORLD.rank == 0:
        print("\n".join(out))
    assert ceil_dr > 2.0 * BAND, f"Debye-reduced ceiling {ceil_dr:.3%} cannot miss"
    assert ceil_wf > 2.0 * BAND, f"wrong-frequency ceiling {ceil_wf:.3%} cannot miss"


@complex_only
def test_cole_cole_plane_wave_matches_k0_sqrt_eps_hat():
    comm = MPI.COMM_WORLD
    out = ["\n[MAT-3 step 1b] FEM, two-pole Cole-Cole medium"]

    alpha_hi, beta_hi, _ = _truth(F_HIGH, COLE_COLE)
    alpha_lo, beta_lo, _ = _truth(F_LOW, COLE_COLE)

    # (i) 127.74 MHz ladder
    rel_c, a_c, b_c, n_c = _solve(12, F_HIGH, COLE_COLE)
    rel_f, a_f, b_f, n_f = _solve(24, F_HIGH, COLE_COLE)
    rate = float(np.log(rel_c / rel_f) / np.log(2.0))
    ea_f = abs(a_f - alpha_hi) / alpha_hi
    eb_f = abs(b_f - beta_hi) / beta_hi
    out.append(f"  (i) 12^3 ({n_c} cells): rel L2 {rel_c:.5e}, alpha {a_c:.5f}, beta {b_c:.5f}")
    out.append(
        f"  (i) 24^3 ({n_f} cells): rel L2 {rel_f:.5e}, alpha {a_f:.5f} ({ea_f:.4%}), "
        f"beta {b_f:.5f} ({eb_f:.4%}); rate {rate:.4f}"
    )

    # (ii) 63.87 MHz at 24^3
    rel_l, a_l, b_l, _ = _solve(24, F_LOW, COLE_COLE)
    ea_l = abs(a_l - alpha_lo) / alpha_lo
    eb_l = abs(b_l - beta_lo) / beta_lo
    out.append(
        f"  (ii) 24^3 @63.87: rel L2 {rel_l:.5e}, alpha {a_l:.5f} ({ea_l:.4%}), "
        f"beta {b_l:.5f} ({eb_l:.4%})"
    )

    # (a') Debye-reduced: every alpha_n = 0, solved at 127.74 (self-consistent solve)
    _, a_dr, b_dr, _ = _solve(16, F_HIGH, DEBYE_REDUCED)
    miss_a_dr = abs(a_dr - alpha_hi) / alpha_hi
    miss_b_dr = abs(b_dr - beta_hi) / beta_hi
    out.append(
        f"  (a') Debye-reduced 16^3: alpha {a_dr:.5f} miss {miss_a_dr:.4%} (ceiling 19.91%), "
        f"beta miss {miss_b_dr:.4%} (printed only)"
    )

    # (b') wrong frequency: eps_hat(63.87) solved at 127.74
    _, a_wf, b_wf, _ = _solve(16, F_HIGH, COLE_COLE, material_f=F_LOW)
    miss_a_wf = abs(a_wf - alpha_hi) / alpha_hi
    miss_b_wf = abs(b_wf - beta_hi) / beta_hi
    out.append(
        f"  (b') wrong-f 16^3: alpha {a_wf:.5f} miss {miss_a_wf:.4%} (ceiling 29.77%), "
        f"beta miss {miss_b_wf:.4%} (printed only)"
    )

    if comm.rank == 0:
        print("\n".join(out))

    assert rel_f < rel_c
    assert rate > 0.9, f"L2 rate {rate:.3f}"
    assert rel_f < BAND, f"24^3 rel L2 {rel_f:.4e}"
    assert ea_f < BAND, f"alpha @127.74 off by {ea_f:.3%}"
    assert eb_f < BAND, f"beta @127.74 off by {eb_f:.3%}"
    assert ea_l < BAND, f"alpha @63.87 off by {ea_l:.3%}"
    assert eb_l < BAND, f"beta @63.87 off by {eb_l:.3%}"
    assert miss_a_dr > BAND, f"Debye-reduced control inside the band: {miss_a_dr:.3%}"
    assert miss_a_wf > BAND, f"wrong-frequency control inside the band: {miss_a_wf:.3%}"
