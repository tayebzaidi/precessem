"""`GEO-27` step 2 (§9 item 377): the 1 g rung's price at degree 2 on exact rectangles.

The exact-rectangle F-small fixture (``leg_sheet_strip_fraction =
GATED_WIDTH_FRACTION``, κ off, solid σ = 800 S/m coil, 50 Ω terminations) built at
``h_p`` = 0.005 m (knob ``h005``) or 0.0025 m (knob ``h0025``, `IMP-3` I14's off
arm), swept at degree 2 with ONE factorisation (``run_n_port_sparameter_sweep``,
the `PORT-19` premise) at ``build_four_port_sweep``'s default 10 MHz.  One env
knob ``FEM_EM_GEO27_STEP2`` in {h005, h0025}; unset => skip; anything else =>
``ValueError``.

Asserted (every band imported, never restated):
    both: (D) d2 dofs = 2E + 2F; (G2) four N_i strip areas = NOMINAL_STRIP_AREA
        within STRIP_BAND; (i) ||S - S^T||/||S|| <= RECIPROCITY_BAND; (P)
        P_src = P_vol + P_sheet per drive <= DISCRETE_IDENTITY_RTOL; (ii) the
        whole-phantom ball vs the tag-3 ``mean_sar`` (power and mass) <=
        EXACT_IDENTITY_RTOL.
    h0025 only: (R) cells within CELL_COUNT_BAND of I14_OFF_CELLS and dofs ==
        I14_OFF_DOFS; (C4-1g) the four cyclic 1 g C4 pairs <= C4_COVARIANCE_BAND;
        (N1g) the mis-paired 1 g control > C4_COVARIANCE_BAND.
Printed, never gated: the price line, the 10 g pairs, the h005 1 g pairs, the
stop-rule line, MAT-4 step 5b's degree-1 ragged comparands.

The centres-and-pairing block is copied from
``tests/validation/test_birdcage_sar_mass_averaged.py:239-331`` at base commit
``cabfed7e`` (``_build_mass_averaged``), with the four fresh factorisations
replaced by the one-factorisation sweep's fields; pairing, centres, radii,
density map, operator and ``QUADRATURE_DEGREE`` are MAT-4's.
"""

from __future__ import annotations

import os
import time

import numpy as np
import pytest
from mpi4py import MPI

from fem_em_solver.ports.sparameters import run_n_port_sparameter_sweep
from fem_em_solver.post.sar import (
    averaging_ball_radius,
    build_density_field,
    mass_averaged_sar,
    mean_sar,
)

from tests.complex_mode import complex_only
from tests.mesh.helpers import global_cell_tag_set
from tests.mesh.test_birdcage_port_sheet_prerequisite import CELL_COUNT_BAND
from tests.mesh.test_birdcage_port_sheets import SHEET_IFACE, _build
from tests.validation.test_birdcage_b1_plus_map import C4_COVARIANCE_BAND
from tests.validation.test_birdcage_b1_quadrature import QUADRATURE_STEP_DEG, _port_index
from tests.validation.test_birdcage_power_identity import (
    DISCRETE_IDENTITY_RTOL,
    _exact_shares,
)
from tests.validation.test_birdcage_sar_1g_rung import PHANTOM_RESOLUTION_1G_RUNG
from tests.validation.test_birdcage_sar_mass_averaged import (
    CENTRE_RADIUS_M,
    EXACT_IDENTITY_RTOL,
    PHANTOM_RHO_KG_PER_M3,
    WHOLE_PHANTOM_BALL_RADIUS_M,
    _global_tag_cell_count,
)
from tests.validation.test_imp3_palace_twin import _sheet_areas
from tests.validation.test_imp3_rod_birdcage import _matrix_quantities, _peak_rss_gib
from tests.validation.test_imp3_tip_ladder_1g import I14_OFF_CELLS, I14_OFF_DOFS
from tests.validation.test_imp3_tip_ratio_1g import NOMINAL_STRIP_AREA
from tests.validation.test_mass_averaged_sar_standard_masses import (
    ONE_GRAM_KG,
    QUADRATURE_DEGREE,
    TEN_GRAM_KG,
)
from tests.validation.test_ops61_exact_strip_ports import STRIP_BAND
from tests.validation.test_port_birdcage_four_port import build_four_port_sweep
from tests.validation.test_port_birdcage_lumped_column import PHANTOM_CELL_TAG
from tests.validation.test_port_lumped_narrowed_sheet import GATED_WIDTH_FRACTION
from tests.validation.test_port_lumped_sheet_sweep import RECIPROCITY_BAND

TAG = "[GEO-27 step2]"
KNOB_NAME = "FEM_EM_GEO27_STEP2"
DEGREE = 2

# Step 1's 0.005 m rung (`20260908T033218Z_GEO-27.log:7152-7155`).
PHANTOM_RESOLUTION_005 = 0.005
RUNGS = {"h005": PHANTOM_RESOLUTION_005, "h0025": PHANTOM_RESOLUTION_1G_RUNG}

# MAT-4 step 5b comparands (degree 1, ragged sheets, same h_p = 0.0025 m, 10 MHz;
# `20260908T200629Z_MAT-4-step5b.log.gz:1990-1993`).  Printed, never bands.
MAT4_STEP5B_ONE_GRAM_PAIRS = {0: 0.000957, 1: 0.001199, 2: 0.001305, 3: 0.001065}
MAT4_STEP5B_ONE_GRAM_CONTROL_K0 = 0.870143

# Predictions (written before W0; printed, never asserted — rule (e)).
PRED = {
    "h005": {
        "cells": (128_000, 136_000), "tag3": 8_500, "dofs": (780_000, 900_000),
        "sweep_s": (130.0, 150.0), "rss_gib": (22.0, 30.0), "window_s": (190.0, 320.0),
        "c4_1g": (5.0e-4, 3.0e-2),
    },
    "h0025": {
        "cells": (201_072, 201_072), "tag3": None, "dofs": (1_284_856, 1_284_856),
        "sweep_s": (230.0, 330.0), "rss_gib": (42.0, 50.0), "window_s": (330.0, 520.0),
        "c4_1g": (1.0e-4, 5.0e-3),
    },
}
PRED_N1G = 0.87
PRED_N1G_CEILING = 1.00
PRED_SEPARATION = 100.0
PRED_II_MAX = 1.0e-13

# Stop rule (printed by W1, applied by the executor before W2).
I14_WINDOW_S, I14_RSS_GIB = 397.0, 44.17
STOP_RSS_GIB, STOP_T_S = 100.0, 540.0
STOP_EXPONENT = 1.37


def _knob():
    raw = os.environ.get(KNOB_NAME)
    if raw is None or raw == "":
        return None
    if raw not in RUNGS:
        raise ValueError(f"{KNOB_NAME}={raw!r}: expected one of {sorted(RUNGS)}")
    return raw


def _say(msg):
    if MPI.COMM_WORLD.rank == 0:
        print(msg, flush=True)


def test_knob_parser_rejects_bad_value(monkeypatch):
    monkeypatch.setenv(KNOB_NAME, "h001")
    with pytest.raises(ValueError):
        _knob()
    monkeypatch.delenv(KNOB_NAME, raising=False)
    assert _knob() is None


@complex_only
def test_geo27_one_gram_rung_degree2():
    mode = _knob()
    if mode is None:
        pytest.skip(f"{KNOB_NAME} unset (h005 | h0025)")
    comm = MPI.COMM_WORLD
    t_body = time.perf_counter()
    h_p = RUNGS[mode]

    # 1. Build at the mode's h_p (sizing to _build, never to the sweep builder).
    t0 = time.perf_counter()
    built = _build(True, phantom_resolution=h_p, leg_sheet_strip_fraction=GATED_WIDTH_FRACTION)
    case = build_four_port_sweep(built=built, build_only=True, degree=DEGREE)
    comm.Barrier()
    t_build = time.perf_counter() - t0
    assert int(case["degree"]) == DEGREE, case["degree"]
    msh = case["mesh"]
    cell_tags = case["cell_tags"]
    cells = int(case["cells"])
    tag3 = _global_tag_cell_count(msh, cell_tags, PHANTOM_CELL_TAG, comm)
    msh.topology.create_entities(1)
    msh.topology.create_entities(2)
    n_e = int(msh.topology.index_map(1).size_global)
    n_f = int(msh.topology.index_map(2).size_global)
    areas = _sheet_areas(msh, case["facet_tags"], comm)
    freq = float(case["problem"].frequency_hz)

    # 2. Before the sweep: cells, (D) via the space after the sweep, (G2), (R).
    _say(f"\n{TAG} mode {mode}: h_p = {h_p} m, f = {freq:.3e} Hz, degree {DEGREE}, "
         f"{cells} cells ({tag3} tag-3), E = {n_e}, F = {n_f}, 2E+2F = {2 * (n_e + n_f)}, "
         f"build {t_build:.2f} s, ranks {comm.size}")
    _say(f"{TAG} (G2) N_i areas " + ", ".join(f"{k} {a:.12e}" for k, a in areas.items())
         + f" vs NOMINAL_STRIP_AREA {NOMINAL_STRIP_AREA:.6e} at STRIP_BAND {STRIP_BAND:.0e}")
    assert len(areas) == 4, areas
    for k, a in areas.items():
        assert abs(a / NOMINAL_STRIP_AREA - 1.0) <= STRIP_BAND, (k, a)
    if mode == "h0025":
        _say(f"{TAG} (R) cells {cells} vs I14_OFF_CELLS {I14_OFF_CELLS} "
             f"(rel {abs(cells / I14_OFF_CELLS - 1.0):.3e}, band {CELL_COUNT_BAND}); "
             f"2E+2F {2 * (n_e + n_f)} vs I14_OFF_DOFS {I14_OFF_DOFS}")
        assert abs(cells / I14_OFF_CELLS - 1.0) <= CELL_COUNT_BAND, (cells, I14_OFF_CELLS)
        assert 2 * (n_e + n_f) == I14_OFF_DOFS, (2 * (n_e + n_f), I14_OFF_DOFS)

    # 3. One-factorisation sweep.
    comm.Barrier()
    t0 = time.perf_counter()
    result = run_n_port_sparameter_sweep(
        case["problem"],
        case["port_defs"],
        lumped_sheet_ports=case["specs"],
        lumped_sheet_facet_tags=case["facet_tags"],
        keep_fields=True,
        degree=DEGREE,
    )
    comm.Barrier()
    t_sweep = time.perf_counter() - t0
    port_ids = tuple(result.port_ids)
    v_space = result.fields[port_ids[0]].e_complex.function_space
    ndofs = int(v_space.dofmap.index_map.size_global * v_space.dofmap.index_map_bs)
    _say(f"{TAG} (D) d2 dofs {ndofs} vs 2E+2F {2 * (n_e + n_f)}; sweep {t_sweep:.2f} s")
    assert ndofs == 2 * (n_e + n_f), (ndofs, n_e, n_f)
    if mode == "h0025":
        assert ndofs == I14_OFF_DOFS, (ndofs, I14_OFF_DOFS)

    mq = _matrix_quantities(result)
    _say(f"{TAG} (i) ||S - S^T||/||S|| = {mq['recip']:.3e} (band {RECIPROCITY_BAND:.0e})")
    assert mq["recip"] <= RECIPROCITY_BAND, mq["recip"]

    omega = 2.0 * np.pi * freq
    p_rel = {}
    for pid in port_ids:
        ex = _exact_shares(case, {"driven": pid, "omega": omega, "fields": result.fields[pid]})
        p_rel[pid] = abs(ex["p_src"] - ex["p_vol"] - ex["sheet_field_total"]) / abs(ex["p_src"])
    _say(f"{TAG} (P) |P_src - P_vol - P_sheet|/|P_src| "
         + ", ".join(f"{p} {r:.3e}" for p, r in p_rel.items())
         + f" (band {DISCRETE_IDENTITY_RTOL:.0e})")
    for pid, r in p_rel.items():
        assert r <= DISCRETE_IDENTITY_RTOL, (pid, r)

    # 4. The copied MAT-4 block (test_birdcage_sar_mass_averaged.py:239-331 @ cabfed7e).
    azimuths = {
        f"P{s['tag'] - SHEET_IFACE}": float(s["azimuth_deg"]) for s in case["sheets"]
    }
    order = sorted(azimuths)
    indices = {
        pid: _port_index(azimuths[pid], azimuths["P1"], QUADRATURE_STEP_DEG)
        for pid in order
    }
    by_k = {indices[pid]: pid for pid in order}
    assert sorted(by_k) == [0, 1, 2, 3], by_k

    all_tags = global_cell_tag_set(msh, cell_tags)
    rho_field = build_density_field(
        msh,
        PHANTOM_RHO_KG_PER_M3,
        cell_tags=cell_tags,
        density_map={int(t): 0.0 for t in sorted(all_tags) if int(t) != PHANTOM_CELL_TAG},
    )
    radii = {
        "1 g": averaging_ball_radius(mass_kg=ONE_GRAM_KG, rho=PHANTOM_RHO_KG_PER_M3),
        "10 g": averaging_ball_radius(mass_kg=TEN_GRAM_KG, rho=PHANTOM_RHO_KG_PER_M3),
    }
    centres = {
        k: (
            CENTRE_RADIUS_M * float(np.cos(np.radians(azimuths[by_k[k]]))),
            CENTRE_RADIUS_M * float(np.sin(np.radians(azimuths[by_k[k]]))),
            0.0,
        )
        for k in range(4)
    }

    def _ball(pid, k, label):
        f = result.fields[pid]
        o = mass_averaged_sar(
            f.e_complex,
            sigma=f.sigma_field,
            rho=rho_field,
            center=centres[k],
            radius=radii[label],
            comm=comm,
            quadrature_degree=QUADRATURE_DEGREE,
        )
        return float(np.real(o["averaged_sar_w_per_kg"]))

    t0 = time.perf_counter()
    one = {(k, k): _ball(by_k[k], k, "1 g") for k in range(4)}
    one.update({(k, (k + 2) % 4): _ball(by_k[k], (k + 2) % 4, "1 g") for k in range(4)})
    ten = {(k, k): _ball(by_k[k], k, "10 g") for k in range(4)}
    whole = mass_averaged_sar(
        result.fields["P1"].e_complex,
        sigma=result.fields["P1"].sigma_field,
        rho=rho_field,
        center=(0.0, 0.0, 0.0),
        radius=WHOLE_PHANTOM_BALL_RADIUS_M,
        comm=comm,
        quadrature_degree=QUADRATURE_DEGREE,
    )
    tagged = mean_sar(
        result.fields["P1"].e_complex,
        sigma=result.fields["P1"].sigma_field,
        rho=PHANTOM_RHO_KG_PER_M3,
        cell_tags=cell_tags,
        comm=comm,
        subdomain_ids=PHANTOM_CELL_TAG,
    )
    t_balls = time.perf_counter() - t0

    def _pair(d, a, b):
        return abs(d[a] - d[b]) / d[b]

    one_pairs = {k: _pair(one, ((k + 1) % 4, (k + 1) % 4), (k, k)) for k in range(4)}
    ten_pairs = {k: _pair(ten, ((k + 1) % 4, (k + 1) % 4), (k, k)) for k in range(4)}
    n1g = {k: _pair(one, (k, (k + 2) % 4), (k, k)) for k in range(4)}
    ii_power = abs(float(np.real(whole["dissipated_power_w"]))
                   / float(np.real(tagged["dissipated_power_w"])) - 1.0)
    ii_mass = abs(float(np.real(whole["mass_kg"]))
                  / (PHANTOM_RHO_KG_PER_M3 * float(np.real(tagged["volume_m3"]))) - 1.0)
    separation = min(n1g.values()) / max(one_pairs.values())

    t_total = time.perf_counter() - t_body
    rss = _peak_rss_gib()
    pr = PRED[mode]

    # 5-7. Print everything, then assert.
    _say(f"{TAG} (ii) whole-phantom ball vs tag-3 mean_sar, P1: power {ii_power:.3e}, "
         f"mass vs rho*V {ii_mass:.3e} (band {EXACT_IDENTITY_RTOL:.0e}; predicted <= {PRED_II_MAX:.0e})")
    _say(f"{TAG} balls a_1g = {radii['1 g'] * 1e3:.4f} mm, a_10g = {radii['10 g'] * 1e3:.4f} mm, "
         f"centres at r0 = {CENTRE_RADIUS_M * 1e3:.1f} mm; by_k = {by_k}; "
         f"2a_1g/h = {2 * radii['1 g'] / h_p:.2f}")
    for k in range(4):
        _say(f"{TAG}   k={k}->{(k + 1) % 4}  1 g C4 {one_pairs[k] * 100:9.4f}%   "
             f"10 g C4 {ten_pairs[k] * 100:9.4f}% (printed)   N1g (c_{(k + 2) % 4}; {k}) "
             f"{n1g[k] * 100:9.4f}%   MAT-4 5b d1 ragged 1 g {MAT4_STEP5B_ONE_GRAM_PAIRS[k] * 100:.4f}%")
    _say(f"{TAG}   SAR_1g diag " + ", ".join(f"{one[(k, k)]:.9e}" for k in range(4))
         + " W/kg; SAR_10g diag " + ", ".join(f"{ten[(k, k)]:.9e}" for k in range(4)) + " W/kg")
    _say(f"{TAG}   worst 1 g C4 {max(one_pairs.values()) * 100:.4f}% "
         f"({'GATED' if mode == 'h0025' else 'PRINTED'} at band {C4_COVARIANCE_BAND * 100:.1f}%; "
         f"predicted [{pr['c4_1g'][0] * 100:.2f}, {pr['c4_1g'][1] * 100:.2f}]%); "
         f"N1g min {min(n1g.values()) * 100:.4f}% (predicted ~{PRED_N1G * 100:.0f}%, ceiling "
         f"{PRED_N1G_CEILING * 100:.0f}%, MAT-4 5b k=0 {MAT4_STEP5B_ONE_GRAM_CONTROL_K0 * 100:.4f}%); "
         f"N-sep {separation:.1f}x (predicted >= {PRED_SEPARATION:.0f}x)")
    _say(f"{TAG} PRICE h_p={h_p} cells={cells} tag3={tag3} E={n_e} F={n_f} d2dofs={ndofs} "
         f"build_s={t_build:.2f} sweep_s={t_sweep:.2f} balls_s={t_balls:.2f} "
         f"body_s={t_total:.2f} peak_rss_sum_gib={rss:.2f} ranks={comm.size}")
    _say(f"{TAG} predicted: cells {pr['cells']}, tag3 {pr['tag3']}, dofs {pr['dofs']}, "
         f"sweep {pr['sweep_s']} s, rss {pr['rss_gib']} GiB, window {pr['window_s']} s")
    if mode == "h005":
        ratio = (I14_OFF_DOFS / ndofs) ** STOP_EXPONENT
        t_hat, r_hat = t_total * ratio, rss * ratio
        fires = r_hat > STOP_RSS_GIB or t_hat > STOP_T_S
        _say(f"{TAG} STOP RULE: T2_hat = {t_hat:.1f} s (limit {STOP_T_S:.0f}), R2_hat = "
             f"{r_hat:.2f} GiB (limit {STOP_RSS_GIB:.0f}); I14 measured {I14_WINDOW_S:.0f} s / "
             f"{I14_RSS_GIB:.2f} GiB -> {'FIRES (do not start W2)' if fires else 'does not fire'}")

    assert ii_power <= EXACT_IDENTITY_RTOL, ii_power
    assert ii_mass <= EXACT_IDENTITY_RTOL, ii_mass
    if mode == "h0025":
        for k, r in one_pairs.items():
            assert r <= C4_COVARIANCE_BAND, (k, r)
        for k, r in n1g.items():
            assert r > C4_COVARIANCE_BAND, (k, r)
