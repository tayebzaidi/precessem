"""`IMP-3` §10 I17 (§9 item 464): the rod cap-ball radius knob.

``implant_rod["tip_ball_radius"]`` [m] is the cap-``Ball`` ``Radius`` itself
(absent / ``None`` keeps today's ``2 * radius``, 2 mm for ``ROD``).  Fixture:
the 1 g-rung exact-rectangle *hole* build at VIn r/4 (``tip_size`` unset),
the `GEO-36` step 4 call.  No solve.

Env knob ``FEM_EM_IMP3_I17`` = ``run`` selects the builds (off, then 2 mm,
then 4 mm, coarse to fine, `GEO-23`); the refusals run ungated.

Asserted (bands imported, never restated):
    (R-b)  bad ``tip_ball_radius`` raises ``ValueError`` matching
           ``implant_rod`` before gmsh is initialised;
    (B0)   the explicit 2 mm build: ``tip_ball_radius_m == 2 * ROD_RADIUS``
           bitwise and cells == ``ROD_1G_RUNG_CELL_RECORD`` (252 118) exactly;
    (A0)   CAD masses, OCC wall area and centroid bitwise equal 2 mm vs 4 mm;
    (A1-A6) per rod build against the window's own off build;
    (D)    N1curl d2 = 2E + 2F, d1 = E, exactly; guard < 1 000 000 cells;
    (Kb)   4 mm: ``tip_ball_radius_m`` bitwise, cells > 252 118, N_ball grows;
    (T1in) 4 mm: every vertex of every near-T1 cell and of the T1 owner cell
           lies within the ball (distance from the top cap centre);
    (T1out) 2 mm: no near-T1 cell lies wholly inside the ball (convexity:
           every near-T1 centroid is >= 2.031 mm from the cap centre).

Printed only (rule (e)): predicted 4 mm cells / N_ball / near-T1 edges, the
near-T1 edge ratio, the guard, and the I30 price + slot-feasibility line.
"""

from __future__ import annotations

import os
import time

import gmsh
import numpy as np
import pytest
from mpi4py import MPI

from fem_em_solver.post.evaluation import point_owner_cells

# Bands, fixture and helpers imported, never restated (`ANS-1`).
from tests.mesh.test_birdcage_conductor_hole import CAD_ANALYTIC_BAND, EXACT
from tests.mesh.test_birdcage_implant_rod import (
    ROD,
    ROD_CENTRE,
    ROD_LENGTH,
    ROD_RADIUS,
    _build_rod,
    _phantom_minus_rod,
    _rod_area,
    _rod_volume,
)
from tests.mesh.test_birdcage_implant_rod_tip import (
    AIR_CELL_TAG,
    FIXTURE_KW,
    PHANTOM_CELL_TAG,
    ROD_1G_RUNG_CELL_RECORD,
    RUNAWAY_CELL_GUARD,
    T1,
    T1_RADIUS,
    _cell_geometry,
    _mean_edge,
    _summarise,
)
from tests.mesh.test_birdcage_leg_gaps import TERMINAL_AREA_BAND
from tests.mesh.test_birdcage_port_sheets import PORT_LOWER, PORT_UPPER, _build
from tests.mesh.test_birdcage_port_tags import LEG_COUNT

KNOB = os.environ.get("FEM_EM_IMP3_I17")
TAG = "[IMP-3 I17]"
TIP_BALL_2MM = 2.0 * ROD_RADIUS
TIP_BALL_4MM = 4.0 * ROD_RADIUS
#: Top cap centre (the ball T1 sits near): ROD_CENTRE + ½·ROD_LENGTH·axis.
TOP_CAP = np.asarray(ROD["centre"], dtype=float) + 0.5 * ROD_LENGTH * np.asarray(
    ROD["axis"], dtype=float
)
#: `ROD_CELL_GUARD_W2` (`tests/validation/test_imp3_rod_birdcage.py:115`),
#: restated for a print only (the I30 decision line).
ROD_CELL_GUARD_W2_PRINTED = 320_000
# GEO-36 W2's unset line (`20261007T170700Z_GEO-36-step4-W2-rungs.log:3752-3753`),
# printed beside (B0) with True / False, never asserted.
GEO36_W2_UNSET_N402 = 4958
GEO36_W2_UNSET_N_BALL = 11938
GEO36_W2_UNSET_N_T1 = 11
GEO36_W2_UNSET_T1_EDGE = 6.983224e-04
# The r/8 window's near-T1 edge ratio (`20261014T070008Z_IMP-3-step7-I15-r8.log:2640`).
R8_EDGE_RATIO_PRINTED = 1.2638
# I30 price basis: item 305 W1 (`20261008T033458Z_IMP-3-item305-W1-d2rod-r4s2.log:
# 2078, 2167, 2726`), rod arm, degree 2, `-n 8`.
PRICE_CELLS = 276_088
PRICE_DOFS = 1_775_954
PRICE_BUILD_S = 74.92
PRICE_SWEEP_S = 267.03
PRICE_WINDOW_S = 459.0
PRICE_RSS_GIB = 53.85
PRICE_EXPONENTS = (1.23, 1.37)
SLOT_WINDOW_CAP_S = 545.0

needs_run = pytest.mark.skipif(KNOB != "run", reason="FEM_EM_IMP3_I17 != run")


def _say(msg):
    if MPI.COMM_WORLD.rank == 0:
        print(msg, flush=True)


def _t1_geometry(msh, comm):
    """Near-T1 vertex geometry relative to TOP_CAP (collective, every rank).

    Returns (n_t1, max vertex distance of near-T1 cells, number of near-T1
    cells whose four vertices all lie within 2 mm and 4 mm of TOP_CAP, T1 owner
    original index, T1 owner mean edge, T1 owner max vertex distance).
    """
    xyz = _cell_geometry(msh)
    centroids = xyz.mean(axis=1)
    near = np.linalg.norm(centroids - T1, axis=1) <= T1_RADIUS
    vdist = np.linalg.norm(xyz - TOP_CAP, axis=2).max(axis=1)  # per cell
    n_t1 = comm.allreduce(int(np.count_nonzero(near)), op=MPI.SUM)
    vmax = comm.allreduce(float(vdist[near].max()) if near.any() else 0.0, op=MPI.MAX)
    cmin = comm.allreduce(
        float(np.linalg.norm(centroids[near] - TOP_CAP, axis=1).min()) if near.any() else np.inf,
        op=MPI.MIN,
    )
    inside = {
        r: comm.allreduce(int(np.count_nonzero(vdist[near] <= r)), op=MPI.SUM)
        for r in (TIP_BALL_2MM, TIP_BALL_4MM)
    }
    owner_orig, owner_local = point_owner_cells(msh, T1[None, :], comm)
    edge, dist = 0.0, 0.0
    if owner_local[0] >= 0:
        c = int(owner_local[0])
        edge = float(_mean_edge(xyz[[c]])[0])
        dist = float(vdist[c])
    return {
        "n_t1": n_t1,
        "t1_vmax": vmax,
        "t1_cmin": cmin,
        "t1_inside": inside,
        "owner": int(owner_orig[0]),
        "owner_edge": comm.allreduce(edge, op=MPI.MAX),
        "owner_vmax": comm.allreduce(dist, op=MPI.MAX),
    }


def _rod_summary(label, request, comm):
    started = time.perf_counter()
    out = _build_rod(implant_rod=dict(ROD, tip_ball_radius=request), **FIXTURE_KW)
    geo = _t1_geometry(out[0], comm)
    s = _summarise(label, request, out, comm)
    del out
    s.update(geo)
    s["wall_t"] = time.perf_counter() - started
    _say(
        f"{TAG} build {label}: near-T1 {s['n_t1']} cells, largest vertex distance from "
        f"top cap {s['t1_vmax']:.6e} m, nearest centroid {s['t1_cmin']:.6e} m, wholly "
        f"inside 2 mm {s['t1_inside'][TIP_BALL_2MM]} / 4 mm {s['t1_inside'][TIP_BALL_4MM]}; "
        f"T1 cell original {s['owner']}, mean edge {s['owner_edge']:.6e} m, largest "
        f"vertex distance {s['owner_vmax']:.6e} m; cells {s['n_cells']} vs "
        f"ROD_CELL_GUARD_W2 {ROD_CELL_GUARD_W2_PRINTED} "
        f"({'over' if s['n_cells'] >= ROD_CELL_GUARD_W2_PRINTED else 'under'})"
    )
    return s


@pytest.fixture(scope="module")
def builds():
    comm = MPI.COMM_WORLD
    off = _summarise("off", None, _build(True, **FIXTURE_KW), comm)
    b2 = _rod_summary("2mm", TIP_BALL_2MM, comm)
    b4 = _rod_summary("4mm", TIP_BALL_4MM, comm)
    return off, b2, b4


# ---------------------------------------------------------------- (R-b) ungated
@pytest.mark.parametrize(
    "value",
    [0.0, -1e-3, float("nan"), float("inf"), ROD_RADIUS, 4.0001e-3, "big", True],
    ids=["zero", "neg", "nan", "inf", "radius", "over_4r", "str", "bool"],
)
def test_rb_refusals_raise_without_meshing(value):
    started = time.perf_counter()
    with pytest.raises(ValueError, match="implant_rod"):
        _build_rod(implant_rod=dict(ROD, tip_ball_radius=value), **FIXTURE_KW)
    assert not gmsh.isInitialized()
    _say(
        f"{TAG} (R-b) tip_ball_radius={value!r} refused in "
        f"{time.perf_counter() - started:.3f} s"
    )


# ---------------------------------------------------------------- (B0)
@needs_run
def test_b0_default_unchanged(builds):
    _, b2, _ = builds
    cad = b2["diag"]["implant_rod_cad"]
    edge_same = f"{b2['t1_edge']:.6e}" == f"{GEO36_W2_UNSET_T1_EDGE:.6e}"
    _say(
        f"{TAG} (B0) 2mm: tip_ball_radius_m {cad['tip_ball_radius_m']!r} vs 2*ROD_RADIUS "
        f"{2.0 * ROD_RADIUS!r}; cells {b2['n_cells']} vs ROD_1G_RUNG_CELL_RECORD "
        f"{ROD_1G_RUNG_CELL_RECORD}; census {b2['census']}; printed beside GEO-36 W2 "
        f"unset: tag 402 {b2['n402']} ({b2['n402'] == GEO36_W2_UNSET_N402}), N_ball "
        f"{b2['n_ball']} ({b2['n_ball'] == GEO36_W2_UNSET_N_BALL}), near-T1 "
        f"{b2['n_t1']} ({b2['n_t1'] == GEO36_W2_UNSET_N_T1}), near-T1 edge "
        f"{b2['t1_edge']:.6e} ({edge_same})"
    )
    assert cad["tip_ball_radius_m"] == 2.0 * ROD_RADIUS
    assert b2["n_cells"] == ROD_1G_RUNG_CELL_RECORD


# ---------------------------------------------------------------- (A0)
@needs_run
def test_a0_cad_is_sizing_invariant(builds):
    _, b2, b4 = builds
    c2, c4 = b2["diag"]["implant_rod_cad"], b4["diag"]["implant_rod_cad"]
    _say(
        f"{TAG} (A0) CAD masses equal "
        f"{b2['diag']['cad_mass_by_group'] == b4['diag']['cad_mass_by_group']}; wall "
        f"area {c2['wall_area_m2']!r} / {c4['wall_area_m2']!r}; centroid "
        f"{c2['wall_centroid_m']!r} / {c4['wall_centroid_m']!r}"
    )
    assert b2["diag"]["cad_mass_by_group"] == b4["diag"]["cad_mass_by_group"]
    assert c2["wall_area_m2"] == c4["wall_area_m2"]
    assert c2["wall_centroid_m"] == c4["wall_centroid_m"]


# ---------------------------------------------------------------- (A1-A6)
@needs_run
def test_a1_to_a6_per_build(builds):
    off, b2, b4 = builds
    ports = range(1, LEG_COUNT + 1)
    expected = {1, AIR_CELL_TAG, PHANTOM_CELL_TAG}
    expected |= {PORT_LOWER + i for i in ports} | {PORT_UPPER + i for i in ports}
    target = _phantom_minus_rod()
    s_off = sum(off["diag"]["cad_mass_by_group"].values())
    low, high = TERMINAL_AREA_BAND
    bound = EXACT * ROD_LENGTH
    _say(f"{TAG} (A5) off {off['n_cells']} cells, tag 402 {off['n402']}")
    assert set(off["census"]) == expected
    assert sum(off["census"].values()) == off["n_cells"]
    assert off["n402"] == 0
    for s in (b2, b4):
        cad = s["diag"]["implant_rod_cad"]
        a1 = abs(s["diag"]["cad_mass_by_group"]["phantom"] / target - 1.0)
        a2 = abs((sum(s["diag"]["cad_mass_by_group"].values()) + cad["rod_mass_m3"]) / s_off - 1.0)
        a3 = abs(cad["wall_area_m2"] / _rod_area() - 1.0)
        a4 = s["area402"] / _rod_area()
        xb, yb, zb = cad["wall_centroid_m"]
        _say(
            f"{TAG} {s['label']}: (A1) {a1:.3e} (band {CAD_ANALYTIC_BAND:.0e}); (A2) "
            f"{a2:.3e} (band {EXACT:.0e}); V_rod,occ/closed "
            f"{cad['rod_mass_m3'] / _rod_volume() - 1.0:+.3e}; (A3) {a3:.3e}; (A4) "
            f"meshed wall ratio {a4:.6f} (band {TERMINAL_AREA_BAND}); (A5) cells "
            f"{s['n_cells']}, census sum {sum(s['census'].values())}, tag 402 "
            f"{s['n402']} ({s['n402_ext']} exterior, {s['n402_phantom']} on tag 3, "
            f"{s['n402_int']} interior); (A6) |x-0.015| {abs(xb - ROD_CENTRE[0]):.3e} "
            f"y {yb:+.3e} z {zb:+.3e} (bound {bound:.1e})"
        )
        assert a1 < CAD_ANALYTIC_BAND
        assert a2 < EXACT
        assert a3 < CAD_ANALYTIC_BAND
        assert low < a4 < high
        assert set(s["census"]) == expected
        assert sum(s["census"].values()) == s["n_cells"]
        assert s["n402"] > 0
        assert s["n402_int"] == 0
        assert s["n402_ext"] == s["n402"]
        assert s["n402_phantom"] == s["n402"]
        assert abs(xb - ROD_CENTRE[0]) <= bound
        assert abs(yb) <= bound
        assert abs(zb) <= bound


# ---------------------------------------------------------------- (D) + guard
@needs_run
def test_d_n1curl_dof_identities(builds):
    for s in builds:
        d = s["dofs"]
        _say(
            f"{TAG} (D) {s['label']}: cells {s['n_cells']} (guard {RUNAWAY_CELL_GUARD}); "
            f"d2 {d[2]} vs 2E+2F {2 * d['E'] + 2 * d['F']}; d1 {d[1]} vs E {d['E']}"
        )
        assert s["n_cells"] < RUNAWAY_CELL_GUARD
        assert d[2] == 2 * d["E"] + 2 * d["F"]
        assert d[1] == d["E"]


# ---------------------------------------------------------------- (Kb)
@needs_run
def test_kb_knob_reached_gmsh(builds):
    _, b2, b4 = builds
    cad = b4["diag"]["implant_rod_cad"]
    _say(
        f"{TAG} (Kb) 4mm: tip_ball_radius_m {cad['tip_ball_radius_m']!r} vs "
        f"{TIP_BALL_4MM!r}; cells {b4['n_cells']} vs {ROD_1G_RUNG_CELL_RECORD} "
        f"(predicted ~350 k, [320 k, 390 k]); N_ball {b4['n_ball']} vs 2mm "
        f"{b2['n_ball']} (predicted ~118 k, [95 k, 135 k]); mean ball edge "
        f"{b4['ball_edge']:.6e} m"
    )
    assert cad["tip_ball_radius_m"] == TIP_BALL_4MM
    assert b4["n_cells"] > ROD_1G_RUNG_CELL_RECORD
    assert b4["n_ball"] > b2["n_ball"]


# ---------------------------------------------------------------- (T1in) / (T1out)
@needs_run
def test_t1_inside_4mm_and_outside_2mm(builds):
    _, b2, b4 = builds
    r4 = b4["diag"]["implant_rod_cad"]["tip_ball_radius_m"]
    ratio = b2["t1_edge"] / b4["t1_edge"]
    _say(
        f"{TAG} (T1in) 4mm: near-T1 {b4['n_t1']} cells (predicted [60, 250]), mean edge "
        f"{b4['t1_edge']:.6e} m (predicted ~3.4e-4, [2.9e-4, 4.0e-4]); largest vertex "
        f"distance {b4['t1_vmax']:.6e} m (predicted [3.1, 3.6] mm) <= {r4!r}; T1 cell "
        f"largest vertex distance {b4['owner_vmax']:.6e} m, mean edge "
        f"{b4['owner_edge']:.6e} m"
    )
    _say(
        f"{TAG} (T1out) 2mm: near-T1 {b2['n_t1']} cells, wholly inside 2 mm "
        f"{b2['t1_inside'][TIP_BALL_2MM]} (asserted 0), nearest centroid "
        f"{b2['t1_cmin']:.6e} m (geometry >= 2.031e-03)"
    )
    _say(
        f"{TAG} near-T1 edge ratio h(2mm)/h(4mm) {b2['t1_edge']:.6e} / "
        f"{b4['t1_edge']:.6e} = {ratio:.4f} (predicted ~2.05, [1.7, 2.4]; the r/8 "
        f"window's h(r/(4rt2))/h(r/8) {R8_EDGE_RATIO_PRINTED})"
    )
    assert b4["n_t1"] > 0
    assert b4["t1_vmax"] <= r4
    assert b4["owner"] >= 0
    assert b4["owner_vmax"] <= r4
    assert b2["n_t1"] > 0
    assert b2["t1_inside"][TIP_BALL_2MM] == 0


# ---------------------------------------------------------------- I30 price
@needs_run
def test_i30_price_printed(builds):
    _, _, b4 = builds
    ratio = b4["dofs"][2] / PRICE_DOFS
    _say(
        f"{TAG} I30 price (rod arm, degree 2, `-n 8`, basis item 305 W1 "
        f"{PRICE_CELLS} cells / d2 {PRICE_DOFS}):"
    )
    _say(f"{TAG}   cells | d2 | build s | e | sweep s | window s | RSS GiB")
    windows = {}
    for e in PRICE_EXPONENTS:
        sweep = PRICE_SWEEP_S * ratio**e
        window = PRICE_WINDOW_S + (sweep - PRICE_SWEEP_S) + (b4["t"] - PRICE_BUILD_S)
        windows[e] = window
        _say(
            f"{TAG}   {b4['n_cells']} | {b4['dofs'][2]} | {b4['t']:.2f} | {e} | "
            f"{sweep:.1f} | {window:.1f} | {PRICE_RSS_GIB * ratio**e:.2f}"
        )
    ok = b4["n_cells"] < ROD_CELL_GUARD_W2_PRINTED and windows[1.37] <= SLOT_WINDOW_CAP_S
    _say(
        f"{TAG} I30 slot-feasible at `-n 8`: cells < {ROD_CELL_GUARD_W2_PRINTED} and "
        f"window <= {SLOT_WINDOW_CAP_S:.0f} s at e = 1.37: {'yes' if ok else 'no'} "
        f"(cells {b4['n_cells']}, window {windows[1.37]:.1f} s)"
    )
    assert b4["dofs"][2] > 0
