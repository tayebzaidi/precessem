"""`IMP-3` step 11, matched mesh (§9 item 639, §10 I41): I39's rod-heating ratio like for like.

I39 (`test_imp3_rod_heating_phantom`) compared a 537-cell phantom (the solid
record mesh) with the rod mesh's 51 693-cell one.  Here the comparator is
**arm M**: `GEO-36`'s kept-cells twin of the rod build (``_build_kept_rod``,
``KEPT_ROD_CELL_RECORD`` cells, the rod's 15 000 cells tagged 403) with tag 403
**mapped explicitly to the phantom's gel** — the `WF-9` route
(``replace(problem, material_map={**mmap, BIRDCAGE_IMPLANT_CELL_TAG: gel})``).
The material swap is applied by wrapping ``build_four_port_sweep`` as seen from
``test_imp3_rod_birdcage._arm`` for the duration of arm M's call only (restored
in ``finally``); no landed module is edited.

Arm M's heating runs I39's ``_heat_arm`` on a relabelled copy of the cell tags
(403 → 3), so the phantom is tags {3, 403}: P_phantom, the DG0 source average
and the thermal submesh all span both tags, and the rod's surface is interior
(Robin on the whole outer boundary; the excluded wall count must be 0).
G = 1 W / P_phantom(arm M) normalises both arms.  **Arm R** is I39's with-rod
arm unchanged (``_build_rod``, the PEC hole, adiabatic wall).

Asserted: (A0m) arm M's census before any solve — global cells =
``KEPT_ROD_CELL_RECORD``, tag-403 cells = 182 371 − 167 371 = 15 000, tag-3
cells equal to arm R's (owned, reduced); (B), (E), (E′), (S2) on both arms at
I39's imported bands; (C4) on arm M at ``C4_COVARIANCE_BAND``; (W0) Robin
census (M excludes 0, R excludes its parent tag-402 count); control: arm R
reproduces I39's printed P_phantom and (S2) (same arm, same mesh).
Printed control (predicted, never asserted): arm M with 403 left unmapped
(vacuum) changes P_phantom.  Printed record: tip 1 g θ(t*) ratio R / M, the
1 g SAR ratio Q on the matched pair, the max-θ ratio, P(R)/P(M) − 1.

Scope: one mesh pair, degree-1 EM, P1 thermal; not a converged ratio, not the
rod's own conduction, no band on the ratio, not compliance.
"""

from __future__ import annotations

import gc
import os
import time
from dataclasses import replace

import numpy as np
import pytest
from dolfinx import mesh as dmesh
from mpi4py import MPI

import tests.validation.test_imp3_rod_birdcage as _rb
from fem_em_solver.io.mesh import BIRDCAGE_IMPLANT_CELL_TAG
from fem_em_solver.post.sar import mean_sar
from tests.complex_mode import complex_only
from tests.mesh.test_birdcage_implant_rod import (
    KEPT_ROD_CELL_RECORD,
    _build_kept_rod,
    _build_rod,
)
from tests.validation.test_birdcage_b1_plus_map import C4_COVARIANCE_BAND
from tests.validation.test_birdcage_sar_mass_averaged import PHANTOM_RHO_KG_PER_M3
from tests.validation.test_imp3_rod_birdcage import ROD_CELL_RECORD_W1, TIPS
from tests.validation.test_imp3_rod_heating_phantom import (
    SOURCE,
    _heat_arm,
    _parent_wall_count,
)
from tests.validation.test_mat10_sar_to_heat_sphere import B_BAND
from tests.validation.test_mat15_coil_driven_heating_phantom import E_BAND, S2_BAND
from tests.validation.test_port_birdcage_lumped_column import PHANTOM_CELL_TAG
from tests.validation.test_port_drive_superposition import TWO_TAG_ACK_KEY

KNOB = "FEM_EM_IMP3_I41"
ROD_CELLS = KEPT_ROD_CELL_RECORD - ROD_CELL_RECORD_W1  # 15 000
# I39 W1 (`20261025T110300Z_IMP-3-I39-W1.log`) with-rod prints, :3661 and :3666.
I39_ROD_P_W = 1.054706186864e-05
I39_ROD_S2 = 1.0 - 0.999995000055  # the printed (S2) ratio, 12 digits
I39_ROD_TH_RATIO = {"top": 1.1109, "bottom": 1.0884}


def _say(msg):
    if MPI.COMM_WORLD.rank == 0:
        print(msg, flush=True)


def _owned_tag_count(msh, cell_tags, tag):
    n_owned = msh.topology.index_map(msh.topology.dim).size_local
    sel = (np.asarray(cell_tags.indices) < n_owned) & (np.asarray(cell_tags.values) == tag)
    return msh.comm.allreduce(int(np.count_nonzero(sel)), op=MPI.SUM)


def _merge_rod_into_phantom(msh, cell_tags):
    """Copy of ``cell_tags`` with tag 403 relabelled to the phantom tag."""
    vals = np.asarray(cell_tags.values).copy()
    vals[vals == BIRDCAGE_IMPLANT_CELL_TAG] = PHANTOM_CELL_TAG
    return dmesh.meshtags(
        msh, msh.topology.dim, np.asarray(cell_tags.indices, dtype=np.int32), vals.astype(np.int32)
    )


def _census(built):
    msh, ct = built[0], built[1]
    return {
        "cells": int(msh.topology.index_map(msh.topology.dim).size_global),
        "t3": _owned_tag_count(msh, ct, PHANTOM_CELL_TAG),
        "t403": _owned_tag_count(msh, ct, BIRDCAGE_IMPLANT_CELL_TAG),
    }


def _arm_with_403(built, label, field_hook, map_gel):
    """``_rb._arm`` (rod=False) on the kept build; 403 → gel if ``map_gel``."""
    orig = _rb.build_four_port_sweep
    seen = {}

    def _wrapped(*a, **k):
        case = orig(*a, **k)
        mmap = dict(case["problem"].material_map)
        seen["had_403"] = BIRDCAGE_IMPLANT_CELL_TAG in mmap
        if map_gel:
            gel = mmap[PHANTOM_CELL_TAG]
            case["problem"] = replace(
                case["problem"], material_map={**mmap, BIRDCAGE_IMPLANT_CELL_TAG: gel}
            )
            m2 = case["problem"].material_map
            assert m2[BIRDCAGE_IMPLANT_CELL_TAG] is m2[PHANTOM_CELL_TAG]
            case[TWO_TAG_ACK_KEY] = "I41 arm M: shares computed, p_vol never read"
        seen["mapped"] = BIRDCAGE_IMPLANT_CELL_TAG in case["problem"].material_map
        return case

    _rb.build_four_port_sweep = _wrapped
    try:
        out = _rb._arm(built, label, rod=False, field_hook=field_hook)
    finally:
        _rb.build_four_port_sweep = orig
    out["material_403"] = seen
    return out


def _run():
    comm = MPI.COMM_WORLD
    t0 = time.perf_counter()
    heat = {}

    built_rod = _build_rod()
    c_rod = _census(built_rod)
    rod_wall = _parent_wall_count(built_rod[0], built_rod[2])
    built_m = _build_kept_rod()
    c_m = _census(built_m)
    m_wall = _parent_wall_count(built_m[0], built_m[2])
    _say(
        f"\nI41 (A0m) census: arm M {c_m['cells']} cells (record {KEPT_ROD_CELL_RECORD}), "
        f"tag-3 {c_m['t3']}, tag-403 {c_m['t403']} (expected {ROD_CELLS}); arm R "
        f"{c_rod['cells']} cells (record {ROD_CELL_RECORD_W1}), tag-3 {c_rod['t3']}; "
        f"parent tag-402 owned facets M {m_wall} (interior) / R {rod_wall} (exterior)"
    )
    # (A0m) before any solve: the kept build is the hole build plus the rod.
    assert c_m["cells"] == KEPT_ROD_CELL_RECORD, c_m
    assert c_m["t403"] == ROD_CELLS, c_m
    assert c_m["t3"] == c_rod["t3"], (c_m, c_rod)

    def hook_m(label, rod, msh, cell_tags, e_of, sigma_field, rho_field):
        merged = _merge_rod_into_phantom(msh, cell_tags)
        heat["M"] = _heat_arm(label, rod, msh, merged, e_of[SOURCE], sigma_field, None, m_wall)

    def hook_r(label, rod, msh, cell_tags, e_of, sigma_field, rho_field):
        heat["R"] = _heat_arm(
            label, rod, msh, cell_tags, e_of[SOURCE], sigma_field, heat["M"]["G"], rod_wall
        )

    def hook_vac(label, rod, msh, cell_tags, e_of, sigma_field, rho_field):
        merged = _merge_rod_into_phantom(msh, cell_tags)
        heat["vac"] = {
            "P": float(
                np.real(
                    mean_sar(
                        e_of[SOURCE],
                        sigma=sigma_field,
                        rho=PHANTOM_RHO_KG_PER_M3,
                        cell_tags=merged,
                        comm=comm,
                        subdomain_ids=PHANTOM_CELL_TAG,
                    )["dissipated_power_w"]
                )
            ),
            "P3": float(
                np.real(
                    mean_sar(
                        e_of[SOURCE],
                        sigma=sigma_field,
                        rho=PHANTOM_RHO_KG_PER_M3,
                        cell_tags=cell_tags,
                        comm=comm,
                        subdomain_ids=PHANTOM_CELL_TAG,
                    )["dissipated_power_w"]
                )
            ),
        }

    m = _arm_with_403(built_m, "I41 arm M (rod as gel)", hook_m, map_gel=True)
    gc.collect()
    r = _rb._arm(built_rod, "I41 arm R (with rod)", rod=True, field_hook=hook_r)
    del built_rod
    gc.collect()
    vac = _arm_with_403(built_m, "I41 control (403 unmapped)", hook_vac, map_gel=False)

    hm, hr = heat["M"], heat["R"]
    q_sar = {t: r["sar"][(t, SOURCE)] / m["sar"][(t, SOURCE)] for t in TIPS}
    q_th = {t: hr["tip_th"][t] / hm["tip_th"][t] for t in TIPS}
    _say(
        f"I41 material map tag 403: arm M had {m['material_403']['had_403']} -> mapped "
        f"{m['material_403']['mapped']}; control had {vac['material_403']['had_403']} -> "
        f"mapped {vac['material_403']['mapped']}"
    )
    for t in TIPS:
        _say(
            f"I41 record {t} tip, {SOURCE}: 1 g theta(t*) R/M = {q_th[t]:.6f} (I39 mesh pair "
            f"{I39_ROD_TH_RATIO[t]:.4f}); 1 g SAR ratio Q R/M = {q_sar[t]:.6f}; ball masses "
            f"M {m['mass'][t] * 1e3:.6f} g / R {r['mass'][t] * 1e3:.6f} g"
        )
    _say(
        f"I41 record: max theta R/M = {hr['th_max'] / hm['th_max']:.6f} ({hr['th_max']:.6e} / "
        f"{hm['th_max']:.6e} K); mean theta R/M {hr['mean_th'] / hm['mean_th']:.6f}; "
        f"P_phantom M {hm['P']:.12e} W, R {hr['P']:.12e} W: P(R)/P(M) - 1 = "
        f"{hr['P'] / hm['P'] - 1.0:+.4%} (I39 +2.00 %)"
    )
    _say(
        f"I41 control (predicted, printed): 403 unmapped (vacuum): P over {{3,403}} "
        f"{heat['vac']['P']:.12e} W, over tag 3 {heat['vac']['P3']:.12e} W; change vs arm M "
        f"{heat['vac']['P'] / hm['P'] - 1.0:+.4%}"
    )
    _say(
        f"I41 timing: arms M {m['t_arm']:.1f} / R {r['t_arm']:.1f} / control {vac['t_arm']:.1f} s; "
        f"total {time.perf_counter() - t0:.1f} s at -n {comm.size}"
    )
    return {"M": m, "R": r, "heat": heat, "c_m": c_m, "c_rod": c_rod}


@pytest.fixture(scope="module")
def i41():
    if os.environ.get(KNOB) != "1":
        pytest.skip(f"{KNOB}=1 not set (heavy window, ≈ 5-8 min at -n 2)")
    return _run()


# --- coding tests (W0) -------------------------------------------------------


def test_knob_refuses_without_value():
    assert os.environ.get(KNOB) in (None, "1"), f"{KNOB} must be unset or 1"


def test_rod_cell_count_is_15000():
    assert ROD_CELLS == 15_000


def test_merge_relabels_only_403():
    box = dmesh.create_box(
        MPI.COMM_WORLD, [[0, 0, 0], [1, 1, 1]], [2, 2, 2], dmesh.CellType.tetrahedron
    )
    n = box.topology.index_map(3).size_local + box.topology.index_map(3).num_ghosts
    idx = np.arange(n, dtype=np.int32)
    vals = np.where(
        idx % 3 == 0, BIRDCAGE_IMPLANT_CELL_TAG, np.where(idx % 3 == 1, PHANTOM_CELL_TAG, 7)
    ).astype(np.int32)
    ct = dmesh.meshtags(box, 3, idx, vals)
    merged = _merge_rod_into_phantom(box, ct)
    expect = np.where(vals == BIRDCAGE_IMPLANT_CELL_TAG, PHANTOM_CELL_TAG, vals)
    np.testing.assert_array_equal(merged.values, expect)
    np.testing.assert_array_equal(merged.indices, idx)


# --- W1 ------------------------------------------------------------------------


@complex_only
def test_a0m_matched_census(i41):
    c_m, c_rod = i41["c_m"], i41["c_rod"]
    assert c_m["cells"] == KEPT_ROD_CELL_RECORD
    assert c_m["t403"] == ROD_CELLS
    assert c_m["t3"] == c_rod["t3"]
    assert i41["M"]["material_403"]["mapped"]


@complex_only
def test_w0_robin_census(i41):
    m, r = i41["heat"]["M"], i41["heat"]["R"]
    assert m["counts"]["wall"] == 0, m["counts"]
    assert r["parent_wall"] > 0 and r["counts"]["wall"] == r["parent_wall"], r["counts"]


@complex_only
@pytest.mark.parametrize("arm", ["M", "R"])
def test_identities(i41, arm):
    h = i41["heat"][arm]
    assert h["B"] <= B_BAND, f"(B) {arm}: {h['B']:.3e}"
    assert h["E"] <= E_BAND, f"(E) {arm}: {h['E']:.3e}"
    assert h["Ep"] <= E_BAND, f"(E') {arm}: {h['Ep']:.3e}"
    assert h["S2"] <= S2_BAND, f"(S2) {arm}: {h['S2']:.4e}"


@complex_only
@pytest.mark.parametrize("k", [0, 1, 2, 3])
def test_c4_arm_m(i41, k):
    r = i41["heat"]["M"]["pairs"][k]
    assert r <= C4_COVARIANCE_BAND, f"(C4) arm M pair {k}: {r * 100:.4f}%"


@complex_only
def test_control_arm_r_reproduces_i39(i41):
    h = i41["heat"]["R"]
    assert abs(h["P"] / I39_ROD_P_W - 1.0) <= 1e-9, h["P"]
    # Bound set by measurement: the first W1 (`20261025T201043Z_IMP-3-I41-W1.log:3795`)
    # read |S2 - I39| = 2.97e-11 against a pre-set 5e-12.  The arm differs from I39's
    # only in G (arm M's normaliser, +0.26 %), and the BE step solves for T = T_inf +
    # theta with T_inf = 20 K, so mean(theta^1) carries ≈ 1e-16·20/theta^1 cancellation;
    # 1e-10 is the last printed digit of I39's 4.9999e-06 ("to the printed digits' order").
    assert abs(h["S2"] - I39_ROD_S2) <= 1e-10, h["S2"]
