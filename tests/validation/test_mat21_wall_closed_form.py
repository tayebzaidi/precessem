"""`MAT-21` step 20 — an acrylic container wall on the gel slab, steady, in
closed form (§9 item 660, §10 P48a).

The symmetric 1D slab |z| ≤ C of gel (k, uniform q) inside a source-free wall
C < |z| ≤ C + d (k_w), Robin h outside:
`utils.analytical.two_layer_slab_steady_temperature`; the effective film
h_eff = (1/h + d/k_w)⁻¹: `wall_effective_heat_transfer`.  Box half-thickness C
imported from step 0 (`test_mat21_box_closed_form.py`), gel k, q, h, T∞ from
`MAT-15` (imported, not restated).

Wall constants.  d = 6 mm is the fixture's choice.  k_w = 0.19 W/m/K and
ρ_w c_w = 1190 kg/m³ × 1470 J/kg/K for PMMA (acrylic) are **nominal handbook
values, citation owed**: no PMMA thermal property exists in `docs/references/`
(searched 2026-10-27; the only Plexiglas entry, Jin ch. 11, is a microwave
permittivity).  The identities (W1), (W2), (K∞) are algebraic and do not depend
on the values being right; ε_w and the centre rises are printed only.

Asserted: (W1) on 41 gel points the two-layer gel temperature equals the bare
Robin slab eigen-series at Bi = h_eff C/k (`robin_slab_eigenvalues` +
`_robin_slab_coefficients`, 4000 terms) to ≤ 1e-12 of the centre rise; (W2)
temperature and flux continuity at |z| = C and the Robin residual at
|z| = C + d, ≤ 1e-12 (rise / qC); (K∞) k_w = 1e9 returns the bare Robin series
at h to ≤ 1e-9 of the rise.  Control (Wneg): the bare series at h misses the
two-layer centre rise by > 1e-3 (relative).  Printed: (W3) the 3D box series
at h_eff vs at h (consistency, no claim), ε_w = ρ_w c_w d/(ρc C).
"""

from __future__ import annotations

import numpy as np

from fem_em_solver.utils.analytical import (
    _robin_slab_coefficients,
    robin_box_steady_temperature,
    robin_slab_eigenvalues,
    two_layer_slab_steady_temperature,
    wall_effective_heat_transfer,
)
from tests.validation.test_mat15_finite_cylinder_closed_form import RHO_C, T_INF, H, K, Q
from tests.validation.test_mat21_box_closed_form import C_BOX, HALF, SERIES_TERMS

D_WALL = 0.006  # m — fixture's choice
K_W = 0.19  # W/m/K — PMMA, nominal handbook value, citation owed
RHO_C_W = 1190.0 * 1470.0  # J/m³/K — PMMA, nominal handbook value, citation owed
SLAB_TERMS = 4000  # measured 2026-10-27: 400 → 2.7e-11, 4000 → 2.8e-14, 40000 → 7.0e-16
W_BAND = 1.0e-12
KINF_BAND = 1.0e-9
WNEG_FLOOR = 1.0e-3
Z_GEL = np.linspace(-C_BOX, C_BOX, 41)


def _say(msg):
    print(msg, flush=True)


def _two_layer(z, k_w=K_W, **kw):
    return two_layer_slab_steady_temperature(
        z, half_thickness=C_BOX, wall=D_WALL, k=K, k_w=k_w, q=Q, h=H, t_inf=T_INF, **kw
    )


def _bare_series(z, h):
    """Bare Robin slab at film h: T∞ + (qC²/k) Σ Cₙ cos(γₙ z/C)/γₙ²."""
    gam = robin_slab_eigenvalues(h * C_BOX / K, SLAB_TERMS)
    c = _robin_slab_coefficients(gam)
    z = np.atleast_1d(np.asarray(z, dtype=float))
    return T_INF + (Q * C_BOX**2 / K) * (np.cos(np.outer(z, gam / C_BOX)) * (c / gam**2)).sum(1)


def _rise0():
    return float(_two_layer(0.0) - T_INF)


# ---------------------------------------------------------------- (W1)
def test_w1_gel_side_equals_bare_slab_at_h_eff():
    h_eff = wall_effective_heat_transfer(H, D_WALL, K_W)
    rise = _rise0()
    err = np.abs(_two_layer(Z_GEL) - _bare_series(Z_GEL, h_eff)).max() / rise
    _say(
        f"MAT21-20 (W1) h {H} -> h_eff {h_eff:.12f} W/m2/K (d {D_WALL} m, k_w {K_W}); "
        f"Bi_eff {h_eff * C_BOX / K:.6f}; centre rise {rise:.12f} K; "
        f"max |two-layer - bare series(h_eff)|/rise over {Z_GEL.size} gel points {err:.3e}"
    )
    assert err <= W_BAND


# ---------------------------------------------------------------- (W2)
def test_w2_interface_continuity_and_outer_robin():
    rise = _rise0()
    qc = Q * C_BOX
    worst_t = worst_f = worst_r = 0.0
    for sgn in (1.0, -1.0):
        zi = sgn * C_BOX
        zo = np.nextafter(zi, sgn * np.inf)  # the wall branch, one ulp outside
        tg, dg = _two_layer(zi, return_gradient=True)
        tw, dw = _two_layer(zo, return_gradient=True)
        dt = abs(float(tg - tw)) / rise
        df = abs(float(-K * dg - (-K_W * dw))) / qc
        ze = sgn * (C_BOX + D_WALL)
        te, de = _two_layer(ze, return_gradient=True)
        res = abs(float(-K_W * de * sgn - H * (te - T_INF))) / qc
        out_flux = float(-K_W * de * sgn)
        _say(
            f"MAT21-20 (W2) face {sgn:+.0f}: |dT| at |z|=C /rise {dt:.3e}; "
            f"|flux jump|/(qC) {df:.3e}; Robin residual at |z|=C+d /(qC) {res:.3e}; "
            f"outward flux {out_flux:.12f} W/m2 vs qC {qc:.12f}; "
            f"T(C) {float(tg):.9f} K, T(C+d) {float(te):.9f} K"
        )
        worst_t, worst_f, worst_r = max(worst_t, dt), max(worst_f, df), max(worst_r, res)
    assert worst_t <= W_BAND
    assert worst_f <= W_BAND
    assert worst_r <= W_BAND


# ---------------------------------------------------------------- (W3) + (Print)
def test_w3_box_at_h_eff_and_wall_share_printed():
    h_eff = wall_effective_heat_transfer(H, D_WALL, K_W)
    kw = dict(half_widths=HALF, k=K, q=Q, t_inf=T_INF, n_terms=SERIES_TERMS)
    box_eff = float(robin_box_steady_temperature(0.0, 0.0, 0.0, h=h_eff, **kw) - T_INF)
    box_bare = float(robin_box_steady_temperature(0.0, 0.0, 0.0, h=H, **kw) - T_INF)
    slab_eff = _rise0()
    slab_bare = float(_bare_series(0.0, H)[0] - T_INF)
    eps_w = RHO_C_W * D_WALL / (RHO_C * C_BOX)
    _say(
        f"MAT21-20 (W3, print) box centre rise at h_eff {box_eff:.9f} K vs at h {box_bare:.9f} K "
        f"(ratio {box_eff / box_bare:.6f}); 1D slab with wall {slab_eff:.9f} K vs bare {slab_bare:.9f} K "
        f"(ratio {slab_eff / slab_bare:.6f}); box/slab at h_eff {box_eff / slab_eff:.6f}"
    )
    _say(
        f"MAT21-20 (Print) eps_w = rho_w c_w d/(rho c C) = {eps_w:.6f} "
        f"(rho_w c_w {RHO_C_W:.4e} J/m3/K, nominal, citation owed)"
    )
    assert np.isfinite(box_eff) and np.isfinite(eps_w)  # print-only test; no claim


# ---------------------------------------------------------------- (K∞)
def test_kinf_reduces_to_bare_robin():
    rise_bare = float(_bare_series(0.0, H)[0] - T_INF)
    err = np.abs(_two_layer(Z_GEL, k_w=1.0e9) - _bare_series(Z_GEL, H)).max() / rise_bare
    _say(
        f"MAT21-20 (Kinf) k_w = 1e9: max |two-layer - bare series(h)|/rise {err:.3e} "
        f"(h_eff {wall_effective_heat_transfer(H, D_WALL, 1.0e9):.15f})"
    )
    assert err <= KINF_BAND


# ---------------------------------------------------------------- (Wneg)
def test_negative_control_bare_slab_misses_wall():
    rise = _rise0()
    bare = float(_bare_series(0.0, H)[0] - T_INF)
    miss = abs(rise - bare) / rise
    _say(
        f"MAT21-20 (Wneg) centre rise with wall {rise:.9f} K, bare at h {bare:.9f} K; "
        f"relative miss {miss:.6f} (floor {WNEG_FLOOR})"
    )
    assert miss > WNEG_FLOOR
