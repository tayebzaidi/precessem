"""`MAT-3` step 1: `TH-6`'s lossy plane wave solved in a Debye material.

The fixture is `TH-6`'s (tests/validation/test_lossy_plane_wave.py): an
``L = 0.1 m`` box, Dirichlet data = the exact plane wave ``ẑ·exp(−jkx)`` on
the whole boundary, interior decay/phase slopes fitted along an off-lattice
line.  The material is a single-pole Debye dielectric with ionic conductivity,

    ε̂(f) = ε∞ + Δε/(1 + j2πfτ) − jσ/(2πfε₀),
    ε∞ = 50, Δε = 30, τ = 1/(2π·100 MHz), σ = 0.5 S/m,

evaluated in numpy (``ColeColeMaterial``) and frozen at the solve frequency
into the solver's ``HomogeneousMaterial`` via ``as_homogeneous(f)``.  The
Dirichlet data and the comparand both come from ``k = k₀√ε̂(f)``.  Nothing in
the boundary data sets the interior decay; that is ε̂ acting through the mass
term, so a wrong σ_eff mapping or a conjugated pole moves α.

Anchors (at `TH-6`'s imported 5 % band): (i) 127.74 MHz, 12³ → 24³, rel L2
< 5 %, rate > 0.9, α and β within 5 % at 24³; (ii) 63.87 MHz, 24³, α and β
within 5 %.  Negative controls at 16³: (a) ε̂(63.87 MHz) solved at
127.74 MHz, (b) dispersion dropped (εᵣ = εs = 80, σ = 0.5): each must miss
the 127.74 MHz Debye α by more than the band.
"""

from __future__ import annotations

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

F_HIGH = 127.74e6
F_LOW = 63.87e6
BAND = 0.05  # TH-6's imported 5 % (PROJECT_PLAN §10 MVP criterion)

EPS_INF = 50.0
DELTA_EPS = 30.0
TAU = 1.0 / (2.0 * np.pi * 100e6)
SIGMA_IONIC = 0.5

DEBYE = ColeColeMaterial.debye(EPS_INF, poles=((DELTA_EPS, TAU),), sigma_ionic=SIGMA_IONIC)


def _truth(f: float) -> tuple[float, float, complex]:
    """(α, β, ε̂) straight from k₀√ε̂(f), cross-checked against TH-6's helpers."""
    eps = DEBYE.complex_permittivity(f)
    k0 = 2.0 * np.pi * f * np.sqrt(MU_0 * EPSILON_0)
    k = k0 * np.sqrt(eps)
    assert k.imag < 0.0
    hm = DEBYE.as_homogeneous(f)
    k_th6 = _wavenumber(hm.sigma, frequency_hz=f, epsilon_r=hm.epsilon_r)
    a_th6, b_th6 = _analytic_alpha_beta(hm.sigma, frequency_hz=f, epsilon_r=hm.epsilon_r)
    # The TH-6 Dirichlet data (built from (σ_eff, ε′)) is k₀√ε̂(f).
    assert abs(k_th6 - k) <= 1e-12 * abs(k)
    assert abs(a_th6 + k.imag) <= 1e-10 * abs(k)
    assert abs(b_th6 - k.real) <= 1e-10 * abs(k)
    return float(-k.imag), float(k.real), complex(eps)


def _solve(n: int, f: float, material_f: float | None = None, hm=None):
    """Solve at ``f`` with the material frozen at ``material_f`` (default f)."""
    if hm is None:
        hm = DEBYE.as_homogeneous(f if material_f is None else material_f)
    assert hm.mu_r == MU_R
    return _solve_plane_wave(n, hm.sigma, frequency_hz=f, epsilon_r=hm.epsilon_r)


@complex_only
def test_debye_plane_wave_matches_k0_sqrt_eps_hat():
    comm = MPI.COMM_WORLD
    out = []

    alpha_hi, beta_hi, eps_hi = _truth(F_HIGH)
    alpha_lo, beta_lo, eps_lo = _truth(F_LOW)
    sig_hi, _ = DEBYE.sigma_epsilon_at(F_HIGH)
    sig_lo, _ = DEBYE.sigma_epsilon_at(F_LOW)
    out.append("\n[MAT-3 step 1] Debye eps_inf=50, de=30, tau=1/(2pi 100 MHz), sigma=0.5")
    out.append(f"  127.74 MHz: eps'={eps_hi.real:.4f} eps''={-eps_hi.imag:.4f} "
               f"sigma_eff={sig_hi:.5f} alpha={alpha_hi:.4f} beta={beta_hi:.4f} "
               "(pred 61.399 / 84.919 / 0.6035 / 12.470 / 24.405)")
    out.append(f"  63.87 MHz:  eps'={eps_lo.real:.4f} eps''={-eps_lo.imag:.4f} "
               f"sigma_eff={sig_lo:.5f} alpha={alpha_lo:.4f} beta={beta_lo:.4f} "
               "(pred 71.308 / 154.33 / 0.5484 / 9.404 / 14.704)")
    out.append(f"  alpha ratio 128/64 = {alpha_hi / alpha_lo:.4f} (pred 1.3261)")

    # (i) 127.74 MHz ladder
    rel_c, a_c, b_c, n_c = _solve(12, F_HIGH)
    rel_f, a_f, b_f, n_f = _solve(24, F_HIGH)
    rate = float(np.log(rel_c / rel_f) / np.log(2.0))
    ea_f = abs(a_f - alpha_hi) / alpha_hi
    eb_f = abs(b_f - beta_hi) / beta_hi
    out.append(f"  (i) 12^3 ({n_c} cells): rel L2 {rel_c:.5e}, alpha {a_c:.5f}, beta {b_c:.5f}")
    out.append(f"  (i) 24^3 ({n_f} cells): rel L2 {rel_f:.5e}, alpha {a_f:.5f} ({ea_f:.4%}), "
               f"beta {b_f:.5f} ({eb_f:.4%}); rate {rate:.4f}")

    # (ii) 63.87 MHz at 24^3
    rel_l, a_l, b_l, n_l = _solve(24, F_LOW)
    ea_l = abs(a_l - alpha_lo) / alpha_lo
    eb_l = abs(b_l - beta_lo) / beta_lo
    out.append(f"  (ii) 24^3 @63.87: rel L2 {rel_l:.5e}, alpha {a_l:.5f} ({ea_l:.4%}), "
               f"beta {b_l:.5f} ({eb_l:.4%})")
    out.append(f"  measured alpha ratio 128/64 = {a_f / a_l:.4f}")

    # (a) wrong frequency: eps_hat(63.87) solved at 127.74 (self-consistent solve)
    _, a_wf, b_wf, _ = _solve(16, F_HIGH, material_f=F_LOW)
    miss_a_wf = abs(a_wf - alpha_hi) / alpha_hi
    miss_b_wf = abs(b_wf - beta_hi) / beta_hi
    out.append(f"  (a) wrong-f 16^3: alpha {a_wf:.5f} miss {miss_a_wf:.4%} (pred 11.80%), "
               f"beta miss {miss_b_wf:.4%} (pred 3.02%, printed only)")

    # (b) dispersion dropped: eps_r = eps_s = 80, sigma = 0.5
    static = ColeColeMaterial.debye(EPS_INF + DELTA_EPS, sigma_ionic=SIGMA_IONIC)
    _, a_nd, b_nd, _ = _solve(16, F_HIGH, hm=static.as_homogeneous(F_HIGH))
    miss_a_nd = abs(a_nd - alpha_hi) / alpha_hi
    miss_b_nd = abs(b_nd - beta_hi) / beta_hi
    out.append(f"  (b) no-dispersion 16^3: alpha {a_nd:.5f} miss {miss_a_nd:.4%} (pred 21.80%), "
               f"beta miss {miss_b_nd:.4%} (pred 5.94%, printed only)")

    if comm.rank == 0:
        print("\n".join(out))

    assert rel_f < rel_c
    assert rate > 0.9, f"L2 rate {rate:.3f}"
    assert rel_f < BAND, f"24^3 rel L2 {rel_f:.4e}"
    assert ea_f < BAND, f"alpha @127.74 off by {ea_f:.3%}"
    assert eb_f < BAND, f"beta @127.74 off by {eb_f:.3%}"
    assert ea_l < BAND, f"alpha @63.87 off by {ea_l:.3%}"
    assert eb_l < BAND, f"beta @63.87 off by {eb_l:.3%}"
    assert miss_a_wf > BAND, f"wrong-frequency control inside the band: {miss_a_wf:.3%}"
    assert miss_a_nd > BAND, f"dispersion-dropped control inside the band: {miss_a_nd:.3%}"
