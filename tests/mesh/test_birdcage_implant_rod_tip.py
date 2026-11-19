"""`GEO-36` step 4 (§9 item 298, chain I15 part 1): the rod tip-refinement knob.

``implant_rod["tip_size"]`` [m] is the cap-``Ball`` ``VIn`` itself (absent /
``None`` keeps today's ``radius / 4``).  I15's fixture is the *hole* build at
the 1 g rung (``phantom_resolution = PHANTOM_RESOLUTION_1G_RUNG``,
``leg_sheet_strip_fraction = GATED_WIDTH_FRACTION``; the
`tests/validation/test_imp3_tip_ratio_1g.py:121-128` call), so its identities
are step 2's A1-A6.  No solve.

Env knob ``FEM_EM_GEO36_STEP4`` in {``noop``, ``rungs``} selects the builds
(off first, then rod builds coarse to fine, `GEO-23`); the refusals run
ungated.

Asserted (bands imported, never restated):
    (R)  bad ``tip_size`` values raise ``ValueError`` matching ``implant_rod``
         before gmsh is initialised;
    (N0) noop: unset and explicit ``radius / 4`` builds are integer-equal
         (cells, owned census, tag-402 facets); unset vs
         ``ROD_1G_RUNG_CELL_RECORD`` within ``CELL_COUNT_BAND``;
    (A0) CAD masses, OCC wall area and centroid bitwise equal to the unset
         build's; ``tip_size_m`` equals the request bitwise;
    (A1-A6) per rod build against the window's own off build;
    (D)  N1curl dofs: degree 2 = 2E + 2F, degree 1 = E, exactly;
    (G)  rungs: N_ball grows strictly unset -> r/(4 sqrt 2) -> r/8;
    guard: every build < 1 000 000 cells.

Printed only (rule (e), no prior measurement): predicted N_ball growth, mean
ball / near-T1 edge, A4 ratio, cells vs predictions and the 320 000 guard, and
the I15 price per rung.
"""

from __future__ import annotations

import importlib.util
import math
import os
import time
from pathlib import Path

import dolfinx
import gmsh
import numpy as np
import pytest
from mpi4py import MPI

from fem_em_solver.io.mesh import BIRDCAGE_IMPLANT_SURFACE_TAG
from tests.mesh.probe_two_torus_conductor_hole import _exterior_facet_area

# Bands, fixture and helpers imported, never restated (`ANS-1`).
from tests.mesh.test_birdcage_conductor_hole import (
    CAD_ANALYTIC_BAND,
    EXACT,
    _owned_census,
)
from tests.mesh.test_birdcage_implant_rod import (
    ROD,
    ROD_CENTRE,
    ROD_LENGTH,
    ROD_RADIUS,
    _build_rod,
    _facet_check,
    _phantom_minus_rod,
    _rod_area,
    _rod_volume,
)
from tests.mesh.test_birdcage_leg_gaps import TERMINAL_AREA_BAND
from tests.mesh.test_birdcage_port_sheet_prerequisite import CELL_COUNT_BAND
from tests.mesh.test_birdcage_port_sheets import PORT_LOWER, PORT_UPPER, _build
from tests.mesh.test_birdcage_port_tags import LEG_COUNT
from tests.validation.test_birdcage_sar_1g_rung import PHANTOM_RESOLUTION_1G_RUNG
from tests.validation.test_port_lumped_narrowed_sheet import GATED_WIDTH_FRACTION

REPO = Path(__file__).resolve().parents[2]


def _load(name, rel):
    spec = importlib.util.spec_from_file_location(name, REPO / rel)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


_FP = _load("imp3_step4_palace_fields", "scripts/probes/imp3_step4_palace_fields.py")
T1 = np.asarray(_FP.PROBE_POINTS[_FP.PROBE_NAMES.index("T1")], dtype=float)

#: The 1 g-rung exact-rectangle rod (hole) build's global cell count.
#: dolfinx 0.11.0.post0 / gmsh 4.15.2-git-657c8e9: 252 118 at `-n 1`
#: (`20261006T111102Z_IMP-3-item263-W1-build.log`) and at `-n 8`
#: (`20261006T172542Z_IMP-3-item284-W3-d2rod.log`).
ROD_1G_RUNG_CELL_RECORD = 252_118
#: The 1 g-rung off build's count, printed beside A5's off count only
#: (`20261006T111102Z_IMP-3-item263-W1-build.log`).
OFF_1G_RUNG_CELL_PRINTED = 201_072
#: `ROD_CELL_GUARD_W2` (`tests/validation/test_imp3_rod_birdcage.py:115`),
#: restated for a print only: `test_imp3_*` is not imported in the real build.
ROD_CELL_GUARD_W2_PRINTED = 320_000
#: Runaway guard (brief): every rung below this, checked before (D).
RUNAWAY_CELL_GUARD = 1_000_000
#: Near-T1 sampling radius [m] (brief: 0.5 mm).
T1_RADIUS = 0.5e-3

# I15 price inputs (`docs/planning/chunks/IMP-3.md:924,948`, item 284 W3).
PRICE_SWEEP_S = 284.88
PRICE_DOFS = 1_620_748
PRICE_RSS_GIB = 50.97
PRICE_EXPONENTS = (1.23, 1.37)
PRICE_WINDOW_BASE_S = 450.0
PRICE_BUILD_BASE_S = 64.0
HEAVY_CEILING_S = 590.0

KNOB = os.environ.get("FEM_EM_GEO36_STEP4")
TIP_R4 = ROD_RADIUS / 4.0
TIP_R4S2 = ROD_RADIUS / (4.0 * math.sqrt(2.0))
TIP_R8 = ROD_RADIUS / 8.0
#: label -> requested tip_size (None = key absent)
PLANS = {
    "noop": (("unset", None), ("r/4", TIP_R4)),
    "rungs": (("unset", None), ("r/(4rt2)", TIP_R4S2), ("r/8", TIP_R8)),
}
PREDICTED_CELLS = {"r/(4rt2)": "~298 k (<= 345 k)", "r/8": "~427 k (<= 609 k)"}
TAG = "[GEO-36 step 4]"
AIR_CELL_TAG = 2
PHANTOM_CELL_TAG = 3

needs_knob = pytest.mark.skipif(KNOB not in PLANS, reason="FEM_EM_GEO36_STEP4 not in {noop, rungs}")
needs_rungs = pytest.mark.skipif(KNOB != "rungs", reason="rungs window only")
needs_noop = pytest.mark.skipif(KNOB != "noop", reason="noop window only")

FIXTURE_KW = {
    "phantom_resolution": PHANTOM_RESOLUTION_1G_RUNG,
    "leg_sheet_strip_fraction": GATED_WIDTH_FRACTION,
}


def _say(msg):
    if MPI.COMM_WORLD.rank == 0:
        print(msg, flush=True)


def _cell_geometry(msh):
    """(owned cell vertex coordinates [n, 4, 3]) — P1 tets, rank-local."""
    tdim = msh.topology.dim
    n_owned = msh.topology.index_map(tdim).size_local
    dofmap = msh.geometry.dofmap[:n_owned, :4]
    return msh.geometry.x[dofmap]


def _mean_edge(xyz):
    """Per-cell mean of the six edge lengths."""
    pairs = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
    return np.mean([np.linalg.norm(xyz[:, a] - xyz[:, b], axis=1) for a, b in pairs], axis=0)


def _ball_and_t1(msh, cell_tags, ball_radius, comm):
    """Global (N_ball, sum of ball mean-edge, N_t1, sum of t1 mean-edge)."""
    tdim = msh.topology.dim
    n_owned = msh.topology.index_map(tdim).size_local
    xyz = _cell_geometry(msh)
    centroids = xyz.mean(axis=1)
    edges = _mean_edge(xyz)
    tag_of = np.zeros(n_owned, dtype=np.int64)
    idx = np.asarray(cell_tags.indices)
    vals = np.asarray(cell_tags.values)
    keep = idx < n_owned
    tag_of[idx[keep]] = vals[keep]
    rod_c = np.asarray(ROD["centre"], dtype=float)
    rod_a = np.asarray(ROD["axis"], dtype=float)
    in_ball = np.zeros(n_owned, dtype=bool)
    for sign in (-1.0, 1.0):
        cap = rod_c + sign * 0.5 * ROD_LENGTH * rod_a
        in_ball |= np.linalg.norm(centroids - cap, axis=1) <= ball_radius
    in_ball &= tag_of == PHANTOM_CELL_TAG
    near_t1 = np.linalg.norm(centroids - T1, axis=1) <= T1_RADIUS
    return (
        comm.allreduce(int(np.count_nonzero(in_ball)), op=MPI.SUM),
        comm.allreduce(float(edges[in_ball].sum()), op=MPI.SUM),
        comm.allreduce(int(np.count_nonzero(near_t1)), op=MPI.SUM),
        comm.allreduce(float(edges[near_t1].sum()), op=MPI.SUM),
    )


def _dofs(msh):
    msh.topology.create_entities(1)
    msh.topology.create_entities(2)
    e = int(msh.topology.index_map(1).size_global)
    f = int(msh.topology.index_map(2).size_global)
    out = {"E": e, "F": f}
    for k in (1, 2):
        space = dolfinx.fem.functionspace(msh, ("N1curl", k))
        out[k] = int(space.dofmap.index_map.size_global * space.dofmap.index_map_bs)
        del space
    return out


def _summarise(label, request, out, comm):
    msh, cells, facets, diag, elapsed = out
    tdim = msh.topology.dim
    s = {
        "label": label,
        "request": request,
        "n_cells": comm.allreduce(msh.topology.index_map(tdim).size_local, op=MPI.SUM),
        "census": _owned_census(msh, cells, tdim, comm),
        "diag": diag,
        "t": elapsed,
    }
    if label == "off":
        s["n402"] = comm.allreduce(
            int(np.count_nonzero(np.asarray(facets.values) == BIRDCAGE_IMPLANT_SURFACE_TAG)),
            op=MPI.SUM,
        )
    else:
        cad = diag["implant_rod_cad"]
        s["area402"] = _exterior_facet_area(msh, facets, BIRDCAGE_IMPLANT_SURFACE_TAG, comm)
        s["n402"], s["n402_ext"], s["n402_phantom"], s["n402_int"], _ = _facet_check(
            msh, cells, facets, BIRDCAGE_IMPLANT_SURFACE_TAG, comm
        )
        s["n_ball"], ball_edge_sum, s["n_t1"], t1_edge_sum = _ball_and_t1(
            msh, cells, cad["tip_ball_radius_m"], comm
        )
        s["ball_edge"] = ball_edge_sum / s["n_ball"] if s["n_ball"] else float("nan")
        s["t1_edge"] = t1_edge_sum / s["n_t1"] if s["n_t1"] else float("nan")
    if s["n_cells"] >= RUNAWAY_CELL_GUARD:
        raise RuntimeError(
            f"{TAG} {label}: {s['n_cells']} cells >= {RUNAWAY_CELL_GUARD} runaway guard"
        )
    s["dofs"] = _dofs(msh)
    _say(
        f"{TAG} build {label}: tip_size request {request!r}; {s['n_cells']} cells; "
        f"census {s['census']}; mesh_wall_time_s {diag['mesh_wall_time_s']:.2f}; "
        f"elapsed {elapsed:.2f} s; E {s['dofs']['E']} F {s['dofs']['F']}; "
        f"N1curl d1 {s['dofs'][1]} d2 {s['dofs'][2]}"
    )
    if label != "off":
        _say(
            f"{TAG} build {label}: implant_rod_cad tip_size_m "
            f"{diag['implant_rod_cad']['tip_size_m']!r} tip_ball_radius_m "
            f"{diag['implant_rod_cad']['tip_ball_radius_m']!r}; tag 402 {s['n402']} "
            f"facets; N_ball {s['n_ball']}, mean ball edge {s['ball_edge']:.6e} m "
            f"(VIn {diag['implant_rod_cad']['tip_size_m']:.6e} m, ratio "
            f"{s['ball_edge'] / diag['implant_rod_cad']['tip_size_m']:.4f}); near-T1 "
            f"({T1_RADIUS:.1e} m) cells {s['n_t1']}, mean edge {s['t1_edge']:.6e} m"
        )
    del msh, cells, facets
    return s


@pytest.fixture(scope="module")
def builds():
    comm = MPI.COMM_WORLD
    off = _summarise("off", None, _build(True, **FIXTURE_KW), comm)
    rods = []
    for label, tip in PLANS[KNOB]:
        rod = ROD if tip is None else dict(ROD, tip_size=tip)
        rods.append(_summarise(label, tip, _build_rod(implant_rod=rod, **FIXTURE_KW), comm))
    return off, rods


# ---------------------------------------------------------------- (R) ungated
@pytest.mark.parametrize(
    "value",
    [0.0, -1e-4, float("nan"), float("inf"), ROD_RADIUS / 3.0, "fine", True],
    ids=["zero", "neg", "nan", "inf", "r_over_3", "str", "bool"],
)
def test_r_refusals_raise_without_meshing(value):
    started = time.perf_counter()
    with pytest.raises(ValueError, match="implant_rod"):
        _build_rod(implant_rod=dict(ROD, tip_size=value), **FIXTURE_KW)
    assert not gmsh.isInitialized()
    _say(f"{TAG} (R) tip_size={value!r} refused in {time.perf_counter() - started:.3f} s")


# ---------------------------------------------------------------- (N0) noop
@needs_noop
def test_n0_unset_equals_explicit_r4_and_record(builds):
    _, (unset, r4) = builds
    ratio = unset["n_cells"] / ROD_1G_RUNG_CELL_RECORD
    _say(
        f"{TAG} (N0) unset {unset['n_cells']} cells / explicit r/4 {r4['n_cells']}; "
        f"census equal {unset['census'] == r4['census']}; tag 402 {unset['n402']} / "
        f"{r4['n402']}; unset / ROD_1G_RUNG_CELL_RECORD {unset['n_cells']} / "
        f"{ROD_1G_RUNG_CELL_RECORD} = {ratio:.6f} (band {CELL_COUNT_BAND}; predicted "
        f"1.000000)"
    )
    assert unset["n_cells"] == r4["n_cells"]
    assert unset["census"] == r4["census"]
    assert unset["n402"] == r4["n402"]
    assert abs(ratio - 1.0) < CELL_COUNT_BAND


# ---------------------------------------------------------------- (A0)
@needs_knob
def test_a0_cad_is_sizing_invariant(builds):
    _, rods = builds
    ref = rods[0]
    for s in rods:
        cad, rcad = s["diag"]["implant_rod_cad"], ref["diag"]["implant_rod_cad"]
        want = TIP_R4 if s["request"] is None else s["request"]
        _say(
            f"{TAG} (A0) {s['label']}: CAD masses equal "
            f"{s['diag']['cad_mass_by_group'] == ref['diag']['cad_mass_by_group']}, wall "
            f"area {cad['wall_area_m2']!r}, centroid {cad['wall_centroid_m']!r}; "
            f"tip_size_m {cad['tip_size_m']!r} vs {want!r}"
        )
        assert s["diag"]["cad_mass_by_group"] == ref["diag"]["cad_mass_by_group"]
        assert cad["wall_area_m2"] == rcad["wall_area_m2"]
        assert cad["wall_centroid_m"] == rcad["wall_centroid_m"]
        assert cad["tip_size_m"] == want


# ---------------------------------------------------------------- (A1-A6)
@needs_knob
def test_a1_to_a6_per_rung(builds):
    off, rods = builds
    ports = range(1, LEG_COUNT + 1)
    expected = {1, AIR_CELL_TAG, PHANTOM_CELL_TAG}
    expected |= {PORT_LOWER + i for i in ports} | {PORT_UPPER + i for i in ports}
    target = _phantom_minus_rod()
    s_off = sum(off["diag"]["cad_mass_by_group"].values())
    low, high = TERMINAL_AREA_BAND
    bound = EXACT * ROD_LENGTH
    _say(
        f"{TAG} (A5) off {off['n_cells']} cells (1 g rung off printed "
        f"{OFF_1G_RUNG_CELL_PRINTED}; ratio {off['n_cells'] / OFF_1G_RUNG_CELL_PRINTED:.6f}),"
        f" census sum {sum(off['census'].values())}, tag 402 {off['n402']}"
    )
    assert set(off["census"]) == expected
    assert sum(off["census"].values()) == off["n_cells"]
    assert off["n402"] == 0
    for s in rods:
        cad = s["diag"]["implant_rod_cad"]
        a1 = abs(s["diag"]["cad_mass_by_group"]["phantom"] / target - 1.0)
        a2 = abs((sum(s["diag"]["cad_mass_by_group"].values()) + cad["rod_mass_m3"]) / s_off - 1.0)
        a3 = abs(cad["wall_area_m2"] / _rod_area() - 1.0)
        a4 = s["area402"] / _rod_area()
        xb, yb, zb = cad["wall_centroid_m"]
        _say(
            f"{TAG} {s['label']}: (A1) {a1:.3e} (band {CAD_ANALYTIC_BAND:.0e}); (A2) "
            f"{a2:.3e} (band {EXACT:.0e}); V_rod,occ/closed {cad['rod_mass_m3'] / _rod_volume() - 1.0:+.3e};"
            f" (A3) {a3:.3e}; (A4) meshed wall ratio {a4:.6f} (band {TERMINAL_AREA_BAND};"
            f" predicted to rise toward 1 from 0.996254); (A5) cells {s['n_cells']} "
            f"(guard {ROD_CELL_GUARD_W2_PRINTED}; predicted {PREDICTED_CELLS.get(s['label'], '-')}),"
            f" census sum {sum(s['census'].values())}, tag 402 {s['n402']} ({s['n402_ext']} "
            f"exterior, {s['n402_phantom']} on tag 3, {s['n402_int']} interior); (A6) "
            f"|x-0.015| {abs(xb - ROD_CENTRE[0]):.3e} y {yb:+.3e} z {zb:+.3e} (bound {bound:.1e})"
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


# ---------------------------------------------------------------- (D)
@needs_knob
def test_d_n1curl_dof_identities(builds):
    off, rods = builds
    for s in (off, *rods):
        d = s["dofs"]
        _say(
            f"{TAG} (D) {s['label']}: d2 {d[2]} vs 2E+2F {2 * d['E'] + 2 * d['F']}; "
            f"d1 {d[1]} vs E {d['E']}"
        )
        assert d[2] == 2 * d["E"] + 2 * d["F"]
        assert d[1] == d["E"]


# ---------------------------------------------------------------- (G) + price
@needs_rungs
def test_g_ball_grows_and_i15_price(builds):
    _, rods = builds
    unset, mid, fine = rods
    for s in (mid, fine):
        pred = (TIP_R4 / s["request"]) ** 3
        _say(
            f"{TAG} (G) {s['label']}: N_ball {s['n_ball']} / unset {unset['n_ball']} = "
            f"{s['n_ball'] / unset['n_ball']:.4f} (predicted (r/4 / h)^3 = {pred:.3f}); "
            f"cells {s['n_cells']} / unset {unset['n_cells']} = "
            f"{s['n_cells'] / unset['n_cells']:.4f}; near-T1 mean edge "
            f"{s['t1_edge']:.6e} vs unset {unset['t1_edge']:.6e} (predicted unchanged)"
        )
    _say(f"{TAG} I15 price (rod arm, degree 2, `-n 8` basis item 284 W3):")
    _say(f"{TAG}   rung | cells | d2 dofs | build s | e | sweep s | window s | RSS GiB | flag")
    for s in rods:
        ratio = s["dofs"][2] / PRICE_DOFS
        for e in PRICE_EXPONENTS:
            sweep = PRICE_SWEEP_S * ratio**e
            window = PRICE_WINDOW_BASE_S + (sweep - PRICE_SWEEP_S) + (s["t"] - PRICE_BUILD_BASE_S)
            rss = PRICE_RSS_GIB * ratio**e
            flag = "over the heavy ceiling at `-n 8`" if window > HEAVY_CEILING_S else "-"
            _say(
                f"{TAG}   {s['label']} | {s['n_cells']} | {s['dofs'][2]} | {s['t']:.2f} | "
                f"{e} | {sweep:.1f} | {window:.1f} | {rss:.2f} | {flag}"
            )
    assert mid["n_ball"] > unset["n_ball"]
    assert fine["n_ball"] > mid["n_ball"]
