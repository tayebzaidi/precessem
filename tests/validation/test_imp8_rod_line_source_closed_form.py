"""`IMP-8` step 0 (§9 item 759): the finite line source in an infinite
gel medium, closed form, numpy and `scipy.integrate.quad` only.

Brief: docs/planning/chunks/IMP-8.md, "§9 item 759 — 2026-11-02".  A uniform
line of power q′ per length on |z′| ≤ L/2, switched on at t = 0, gives
θ = q′ / (4πk) · Φ with Φ(r, z, t) = ∫ erfc(R / (2√(κt))) / R dz′ and
R² = r² + (z − z′)²: the time integral of the instantaneous point source,
integrated along the line (Carslaw & Jaeger 2nd ed. ch. X, the citation from
memory, the text not being in the repo; the derivation is the brief's).

Asserted: (Inf) the infinite-line limit E₁(r² / 4κt), (St) the steady limit
and its first residual term, (En) the two energy identities, (Far) the far
steady limit L / r, (Mono).  Negative controls (c-κ), (c-fin), (c-L) are
predicted (rule (e)) and asserted only as missing their anchor's band,
measured printed beside predicted.  Scope: a line of zero radius, uniform
power, no walls, no rod heat capacity; no F2182 figure.

Step 1 (§9 item 760, the same chunk file): two tip sources of the same total
power, Φ₂ = (L / 2) Σ_± erfc(R_± / (2√(κt))) / R_±, R_±² = r² + (z ∓ L/2)².
Asserted: (Rate) the midpoint sum of point sources against Φ, (Sym), (Far2).
Control (c-one) is predicted and asserted only as missing (Far2)'s bound.
The (Rate) band is the review's, provisional, for the weekly to ratify.  Two
shapes of a source: which shape a real rod has is not claimed.
"""

from __future__ import annotations

import numpy as np
from mpi4py import MPI
from scipy import special
from scipy.integrate import quad

from tests.mesh.test_birdcage_implant_rod import ROD_LENGTH
from tests.validation.test_imp4_pec_cylinder_closed_form import B1_BAND
from tests.validation.test_mat12_transient_sar_heating_sphere import RHO_C
from tests.validation.test_mat15_coil_driven_heating_phantom import K_GEL

L = ROD_LENGTH  # 0.030, test_birdcage_implant_rod.py:94
KAPPA = K_GEL / RHO_C

# Fixture (rule (f) parameters, not bands).
MM = 1e-3
R_NEAR = 2.0 * MM  # 1 mm off the 1 mm rod's surface; the line is its far-field comparand
T_INF_LIMIT = 10.0  # s, the tips still unfelt at R_NEAR, z = 0
T_FINITE = 120.0  # s, (c-fin): the tips felt
T_STEADY = 1e9  # s
T_LADDER = (10.0, 30.0, 100.0, 300.0, 900.0, 3600.0)  # s
T_PRINT = 900.0  # s
ST_POINTS = ((2.0, 0.0), (2.0, 16.0), (5.0, 30.0), (20.0, 0.0))  # (r, z), mm
QUAD = dict(limit=800, epsabs=1e-13, epsrel=1e-13)
RATE_POINT = (20.0, 0.0)  # (r, z), mm
RATE_LADDER = (4, 8, 16, 32)
RATE_LO, RATE_HI = 1.95, 2.05  # item 760, provisional
SYM_POINT = (2.0, 10.0)  # (r, |z|), mm
FAR2_A = 0.1  # L / r


def _p(*args):
    if MPI.COMM_WORLD.rank == 0:
        print(*args, flush=True)


def _phi(r, z, t, length=L, kappa=KAPPA):
    """Φ(r, z, t), split at z′ = z when the probe's foot lies on the line."""

    def integrand(zp):
        big_r = np.hypot(r, z - zp)
        return special.erfc(big_r / (2.0 * np.sqrt(kappa * t))) / big_r

    lo, hi = -0.5 * length, 0.5 * length
    pieces = [(lo, z), (z, hi)] if lo < z < hi else [(lo, hi)]
    return sum(quad(integrand, a, b, **QUAD)[0] for a, b in pieces)


def _phi_steady(r, z, length=L):
    return np.arcsinh((z + 0.5 * length) / r) - np.arcsinh((z - 0.5 * length) / r)


def _point(big_r, t, kappa=KAPPA):
    return special.erfc(big_r / (2.0 * np.sqrt(kappa * t))) / big_r


def _phi_midpoint(r, z, t, n, length=L):
    """The N-point midpoint sum of Φ's integrand."""
    zp = length * ((np.arange(n) + 0.5) / n - 0.5)
    return length / n * np.sum(_point(np.hypot(r, z - zp), t))


def _phi_tips(r, z, t, length=L, signs=(1.0, -1.0)):
    """Φ₂, the tip pair of the same total power; t = inf is the steady limit."""
    return 0.5 * length * sum(_point(np.hypot(r, z - s * 0.5 * length), t) for s in signs)


def _e1_limit(r, t, kappa=KAPPA):
    return special.exp1(r * r / (4.0 * kappa * t))


def test_inf_infinite_line_limit():
    phi = _phi(R_NEAR, 0.0, T_INF_LIMIT)
    rel = abs(phi / _e1_limit(R_NEAR, T_INF_LIMIT) - 1.0)
    _p(f"[imp8] (Inf) r=2 mm t=10 s Φ={phi:.7g} rel={rel:.3e} band={B1_BAND:.0e}")
    assert rel <= B1_BAND


def test_st_steady_limit_residual_first_term():
    first = L / np.sqrt(np.pi * KAPPA * T_STEADY)
    for r_mm, z_mm in ST_POINTS:
        r, z = r_mm * MM, z_mm * MM
        residual = _phi_steady(r, z) - _phi(r, z, T_STEADY)
        r_max = np.hypot(r, abs(z) + 0.5 * L)
        x = r_max / (2.0 * np.sqrt(KAPPA * T_STEADY))
        bound = x * x / 3.0  # the erf series' own next term over its first
        rel = abs(residual / first - 1.0)
        _p(f"[imp8] (St) ({r_mm:g}, {z_mm:g}) mm rel={rel:.3e} bound={bound:.3e}")
        assert rel <= bound


def test_en_energy_identities():
    e1_mass = sum(
        quad(special.exp1, a, b, limit=800, epsabs=1e-13, epsrel=1e-13)[0]
        for a, b in ((0.0, 1.0), (1.0, np.inf))
    )
    point_mass = quad(lambda u: u * special.erfc(u), 0.0, np.inf, **QUAD)[0]
    line_miss, point_miss = abs(e1_mass - 1.0), abs(point_mass - 0.25)
    _p(f"[imp8] (En) ∫E₁ miss={line_miss:.3e} ∫u erfc miss={point_miss:.3e}")
    assert line_miss <= B1_BAND
    assert point_miss <= B1_BAND


def test_far_steady_limit():
    r = 10.0 * L
    a = L / r
    rel = abs(_phi_steady(r, 0.0) / a - 1.0)
    bound = a * a / 24.0
    _p(f"[imp8] (Far) r=10 L rel={rel:.3e} bound={bound:.3e}")
    assert rel <= bound


def test_mono_rises_below_steady():
    values = [_phi(R_NEAR, 0.0, t) for t in T_LADDER]
    steady = _phi_steady(R_NEAR, 0.0)
    _p("[imp8] (Mono) " + " ".join(f"{v:.7g}" for v in values) + f" steady={steady:.7g}")
    assert all(b > a for a, b in zip(values, values[1:], strict=False))
    assert values[-1] < steady


def test_controls_fire():
    limit = _e1_limit(R_NEAR, T_INF_LIMIT)
    c_kappa = abs(_phi(R_NEAR, 0.0, T_INF_LIMIT, kappa=0.25 * KAPPA) / limit - 1.0)
    c_fin = abs(_phi(R_NEAR, 0.0, T_FINITE) / _e1_limit(R_NEAR, T_FINITE) - 1.0)
    c_l = abs(_phi_steady(R_NEAR, 0.0, length=2.0 * L) / _phi(R_NEAR, 0.0, T_STEADY) - 1.0)
    _p(f"[imp8] (c-κ) measured={c_kappa:.3e} predicted=0.954 (ceiling)")
    _p(f"[imp8] (c-fin) measured={c_fin:.3e} predicted=9.5e-04")
    _p(f"[imp8] (c-L) measured={c_l:.3e} predicted=0.254")
    assert c_kappa > B1_BAND
    assert c_fin > B1_BAND
    assert c_l > 0.1


def test_printed_levels():
    """Printed, no band: Φ and Φ∞ at 900 s, and as K per watt of total power."""
    for z_mm in (0.0, 15.0):
        r, z = R_NEAR, z_mm * MM
        phi, steady = _phi(r, z, T_PRINT), _phi_steady(r, z)
        per_watt = 1.0 / (4.0 * np.pi * K_GEL * L)
        _p(
            f"[imp8] (print) (2, {z_mm:g}) mm t=900 s Φ={phi:.7g} Φ∞={steady:.7g}"
            f" K/W: {phi * per_watt:.7g} / {steady * per_watt:.7g}"
        )
        assert steady > phi > 0.0


def test_rate_midpoint_sum_against_line():
    r, z = RATE_POINT[0] * MM, RATE_POINT[1] * MM
    phi = _phi(r, z, T_PRINT)
    errors = [abs(_phi_midpoint(r, z, T_PRINT, n) / phi - 1.0) for n in RATE_LADDER]
    slope = -np.polyfit(np.log(RATE_LADDER), np.log(errors), 1)[0]
    _p("[imp8] (Rate) rel " + " ".join(f"{e:.3e}" for e in errors) + f" slope={slope:.4f}")
    assert RATE_LO <= slope <= RATE_HI


def test_sym_even_in_z():
    r, z = SYM_POINT[0] * MM, SYM_POINT[1] * MM
    line = abs(_phi(r, z, T_PRINT) / _phi(r, -z, T_PRINT) - 1.0)
    tips = abs(_phi_tips(r, z, T_PRINT) / _phi_tips(r, -z, T_PRINT) - 1.0)
    _p(f"[imp8] (Sym) Φ miss={line:.3e} Φ₂ miss={tips:.3e} band={B1_BAND:.0e}")
    assert line <= B1_BAND
    assert tips <= B1_BAND


def _far2_miss(signs=(1.0, -1.0)):
    a = FAR2_A
    r = L / a
    line, tips = _phi_steady(r, 0.0), _phi_tips(r, 0.0, np.inf, signs=signs)
    return abs(abs(tips - line) / a - a * a / 12.0), line - tips


def test_far2_tip_pair_against_line():
    miss, gap = _far2_miss()
    bound = FAR2_A**4
    _p(f"[imp8] (Far2) a=0.1 miss={miss:.3e} bound={bound:.3e} Φ∞−Φ₂∞={gap:.4e}")
    assert miss <= bound
    assert gap > 0.0


def test_control_one_source_dropped():
    miss = _far2_miss(signs=(1.0,))[0]
    _p(f"[imp8] (c-one) measured={miss:.3e} predicted=0.5 bound={FAR2_A**4:.3e}")
    assert miss > FAR2_A**4


def test_printed_tip_pair_against_line():
    """Printed, no band: Φ₂, Φ and their ratio at 900 s."""
    for z_mm in (15.0, 0.0):
        r, z = R_NEAR, z_mm * MM
        tips, line = _phi_tips(r, z, T_PRINT), _phi(r, z, T_PRINT)
        _p(f"[imp8] (print2) (2, {z_mm:g}) mm Φ₂={tips:.7g} Φ={line:.7g} Φ₂/Φ={tips / line:.7g}")
        assert tips > 0.0
