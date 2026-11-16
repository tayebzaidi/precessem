"""Unit tests: layered surface impedance closed form (`TH-18` step 1, numpy anchors n1–n6).

The transmission-line recursion of Jin §3.4.2 (3.87)–(3.89) on the
good-conductor line constants of `TH-14` (PROJECT_PLAN §9 item 82). σ = 1e4 S/m
at f = 291.35 MHz (the `TH-14` cavity's TE₁₀₁ frequency).
"""

import numpy as np
import pytest

from fem_em_solver.core.layered_impedance import (
    ConductorLayer,
    DielectricLayer,
    conductor_skin_depth_m,
    conductor_surface_resistance_ohm,
    layered_surface_impedance,
)
from fem_em_solver.utils.constants import C_0, ETA_0, MU_0

SIGMA = 1.0e4
OMEGA = 2.0 * np.pi * 291.35e6
DELTA = conductor_skin_depth_m(OMEGA, SIGMA)
Z_C = (1.0 + 1.0j) * conductor_surface_resistance_ohm(OMEGA, SIGMA)


@pytest.mark.parametrize("z_b", [0.0, ETA_0, 1.0e6])
def test_n1_thick_layer_is_leontovich(z_b):
    z = layered_surface_impedance(OMEGA, [ConductorLayer(SIGMA, 40.0 * DELTA)], z_b)
    rel = abs(z - Z_C) / abs(Z_C)
    print(f"[TH-18] (n1) Z_b={z_b}: |Z_s-Z_c|/|Z_c|={rel:.3e} (bound 1e-12)")
    assert rel <= 1e-12


@pytest.mark.parametrize("z_b", [0.0, ETA_0, 1.0e6])
def test_n2_zero_thickness_is_backing(z_b):
    z = layered_surface_impedance(OMEGA, [ConductorLayer(SIGMA, 0.0)], z_b)
    err = abs(z - z_b) / max(abs(z_b), 1.0)
    print(f"[TH-18] (n2) Z_b={z_b}: |Z_s-Z_b|={err:.3e} (bound 1e-15)")
    assert err <= 1e-15


def test_n3_electrically_thin_sheet_on_pec_is_inductance():
    t = 1.0e-4 * DELTA
    z = layered_surface_impedance(OMEGA, [ConductorLayer(SIGMA, t)], 0.0)
    rel = abs(z / (1.0j * OMEGA * MU_0 * t) - 1.0)
    print(f"[TH-18] (n3) |Z_s/(j w mu0 t) - 1|={rel:.3e} (bound 1e-7)")
    assert rel <= 1e-7


@pytest.mark.parametrize("z_b", [0.0, ETA_0])
def test_n4_cascade_equals_single_layer(z_b):
    two = layered_surface_impedance(
        OMEGA, [ConductorLayer(SIGMA, 0.3 * DELTA), ConductorLayer(SIGMA, 0.2 * DELTA)], z_b
    )
    one = layered_surface_impedance(OMEGA, [ConductorLayer(SIGMA, 0.5 * DELTA)], z_b)
    rel = abs(two - one) / abs(one)
    print(f"[TH-18] (n4) Z_b={z_b}: cascade vs single rel={rel:.3e} (bound 1e-12)")
    assert rel <= 1e-12


def test_n5_dielectric_quarter_and_half_wave():
    eps_r = 4.0
    lam_d = 2.0 * np.pi * C_0 / (OMEGA * np.sqrt(eps_r))
    quarter = layered_surface_impedance(OMEGA, [DielectricLayer(eps_r, lam_d / 4.0)], ETA_0)
    half = layered_surface_impedance(OMEGA, [DielectricLayer(eps_r, lam_d / 2.0)], ETA_0)
    rq = abs(quarter - ETA_0 / 4.0) / (ETA_0 / 4.0)
    rh = abs(half - ETA_0) / ETA_0
    print(f"[TH-18] (n5) lambda/4: {quarter} rel={rq:.3e}; lambda/2: {half} rel={rh:.3e} (bound 1e-12)")
    assert rq <= 1e-12
    assert rh <= 1e-12


def test_n6_passive_real_part():
    worst = np.inf
    for z_b in (0.0, ETA_0 / 2.0, ETA_0):
        for ratio in np.logspace(-3, np.log10(50.0), 200):
            z = layered_surface_impedance(OMEGA, [ConductorLayer(SIGMA, ratio * DELTA)], z_b)
            worst = min(worst, z.real)
            assert z.real >= 0.0, (z_b, ratio, z)
    print(f"[TH-18] (n6) min Re Z_s over t/delta in [1e-3, 50] x 3 backings = {worst:.6e}")
