"""`TH-20` step 0a (§9 item 176, chain H H13-0a) — the quadrant cutter on the
F-human hole layout (LH): ring gaps, no leg gaps, 16 legs, ``as_hole=True``,
longitudinal ring sheets.  Mesh only, no solve.

Brief: ``docs/planning/chunks/TH-20.md`` "§9 item 176".  The heavy tests are
selected by environment only (``FEM_EM_TH20_STEP0=1``); the refusal test (Q7)
does no meshing and always runs.  The flag-off full hole mesh is built FIRST,
then the quarter (`GEO-23` contamination).

Anchors (asserted, bands imported):
  (Q0) flag-off full hole mesh: 369 387 cells exactly
       (`20261002T020254Z_WF-12-step1-W1b.log:6176`), tags {2, 3, 101–148, 217–248};
  (Q1) quarter tags = the 22 predicted; 0 owned vertices at x or y < -1e-12;
       owned census sums to size_global;
  (Q2) per cut plane 501 / 502: owned exterior facets on the plane = the tagged
       set, flat to 1e-12 m, area r_e·2z_e - 2 r_leg·L at EXACT_IDENTITY_RTOL;
       cells adjacent to 501 ∪ 502 carry exactly {2, 3}; outer tag 1 area
       4 r_e z_e + 2 r_e² at EXACT_IDENTITY_RTOL, outer ∩ cut = 0;
  (Q3) 0 tag-401 facets on either plane; every owned exterior facet in exactly
       one of {1, 501, 502, 401}; |4 A401_q / A401_f - 1| <=
       TH20_QUARTER_WALL_AREA_RTOL (rule (f) band, §9 item 197);
  (Q3-b) the quarter's tag-401 wall bordering the kept ring-port cells
       (KEPT_RING_PORT_TAGS) = the full mesh's tag-401 wall bordering the SAME
       tags, |A_R,q/A_R,f - 1| <= RING_PORT_WALL_RTOL, equal facet counts; guard
       min(x, y) > 0 over the full mesh's masked facets (asserted);
  (Q4) kept port volumes = CAD mass = full-mesh same-tag volume at EXACT; the 8
       kept sheets' areas = chord·w at SHEET_AREA_BAND; no unkept ordinal;
  (Q5) the D2 port map on the sheet facet centroids (closed form per ring);
  (Q6) phantom 4·CAD_q/CAD_f - 1 within CAD_ANALYTIC_BAND;
  (Q7) the refusals, raised without meshing.
Controls: (C1) the quarter checker on the full mesh, (C2) the mirror misplaced
by π/32 — printed in W1, asserted under ``FEM_EM_TH20_STEP0_CONTROL=1``;
(C3) the cavity rule without the plane walls, asserted c3 >
TH20_QUARTER_WALL_AREA_RTOL under the CONTROL flag (item 197; +3.780 measured
in item 176's W1 / W2); (C4) a planted missing face — the quarter's facets
bordering tag 117 removed from Q3-b — asserted under the CONTROL flag as
ΔN_R = -56 and |ΔA_R/A_R| >= 1e6 · RING_PORT_WALL_RTOL; the same planted face
under Q3 alone is printed, never asserted (it sits inside the Q3 band, which is
why Q3-b exists).

`_plane_readings` / `_facet_nodes` are imported from
`tests/mesh/test_birdcage_symmetry_quadrant.py` unchanged (that module is not
edited).
"""

from __future__ import annotations

import math
import os
import resource
import time

import numpy as np
import pytest
from mpi4py import MPI

STEP0_ENV = "FEM_EM_TH20_STEP0"
CONTROL = os.environ.get("FEM_EM_TH20_STEP0_CONTROL") == "1"

# (Q0) `20261002T020254Z_WF-12-step1-W1b.log:6176`.
FULL_HOLE_CELL_RECORD = 369387
# (Q4) CAD masses of the full hole mesh (`…W1b.log:6174`).
CAD_LEG_BOX = 7.871720117e-6
CAD_RING_LOWER = 8.437757382e-7
CAD_RING_UPPER = 8.723783056e-7
# Closed forms written in the brief.
A_CUT_PREDICTED = 9.8155102e-2
A_OUTER_PREDICTED = 3.16254857e-1
A_CAVITY_QUARTER_CAD_PREDICTED = 5.2369004e-2
C3_PREDICTED = 3.75
PREDICTED_CELLS = (92_000, 140_000)  # central ≈ 115 k
PREDICTED_MESH_S = (25.0, 50.0)

KEPT_LEG = (1, 2, 3, 4)
KEPT_RING = (17, 18, 19, 20, 33, 34, 35, 36)
QUARTER_TAGS = (
    {2, 3}
    | {100 + o for o in KEPT_LEG}
    | {100 + o for o in KEPT_RING}
    | {200 + o for o in KEPT_RING}
)
FULL_TAGS = {2, 3} | set(range(101, 149)) | set(range(217, 249))
# (Q3-b) the kept ring-port cell tags, {117–120, 133–136, 217–220, 233–236}.
KEPT_RING_PORT_TAGS = {100 + o for o in KEPT_RING} | {200 + o for o in KEPT_RING}
CAVITY_TAG = 401

# (Q3) band for THIS comparison (quarter vs full F-human hole mesh, 1.6 mm
# conductor_resolution), registered by the 2026-10-03 03:00 daily review under
# §9 rule (f) (docs/planning/chunks/TH-20.md "Ruling, 2026-10-03 03:00"; the
# 10-03 weekly ratifies or reverts).  Measured by item 182's three -n 1 probe
# logs: 20261003T004412Z_TH-20-probe-wall-h1p6.log,
# 20261003T004644Z_TH-20-probe-wall-h2p4.log,
# 20261003T004837Z_TH-20-probe-wall-h3p2.log: 4q/f - 1 = 1.498e-4 / 3.477e-4 /
# 8.916e-4 at 1.6 / 2.4 / 3.2 mm, decision rule TRIANGULATION (two independent
# inscribed triangulations of the curved wall; the band R ring-port wall exact,
# equal facet counts).  1.69e-3 = 1.9 x the worst rung (8.916e-4, h3p2).
# The imported F_HUMAN_WALL_AREA_RTOL (5e-5) was registered on a different
# comparison (hole vs solid on one full sizing) and is NOT edited here.  The
# band alone catches a missing face only above ≈ 9.6e-5 m^2, or an extra face
# above ≈ 8.0e-5 m^2 — a planted tag-117 ring-port face (5.27e-5 m^2) sits
# inside it — hence the exact sub-identity Q3-b.
TH20_QUARTER_WALL_AREA_RTOL = 1.69e-3
# (C4) the planted tag-117 face: facet count, measured by the h1p6 probe.
C4_PLANTED_FACETS = 56
PLANE_TOL = 1.0e-12


def _on() -> bool:
    return os.environ.get(STEP0_ENV, "").strip() == "1"


def _skip_unless_on():
    if not _on():
        pytest.skip(f"{STEP0_ENV} unset: TH-20 step 0a is selected by environment only")


def _kwargs():
    from tests.validation.test_birdcage_f_human_rung import F_HUMAN_RING_RADIUS, _params

    kwargs = _params(F_HUMAN_RING_RADIUS, scale_sizing=False)
    kwargs["ring_sheet_orientation"] = "longitudinal"
    kwargs["as_hole"] = True
    return kwargs


def _rank0_print(comm, msg):
    if comm.rank == 0:
        print(msg, flush=True)


def _sheet_centroids(msh, cells, ordinals, comm):
    """Area-weighted centroid and area of each ring sheet, owned facets, reduced."""
    import dolfinx

    from fem_em_solver.io.mesh import _interface_facet_tags
    from tests.mesh.test_birdcage_port_sheets import PORT_LOWER, PORT_UPPER, SHEET_IFACE

    tags = _interface_facet_tags(
        msh, cells, {SHEET_IFACE + o: (PORT_LOWER + o, PORT_UPPER + o) for o in ordinals}
    )
    fdim = msh.topology.dim - 1
    n_owned = msh.topology.index_map(fdim).size_local
    idx = np.asarray(tags.indices)
    val = np.asarray(tags.values)
    out = {}
    for o in ordinals:
        f = idx[(val == SHEET_IFACE + o) & (idx < n_owned)]
        if len(f):
            g = dolfinx.cpp.mesh.entities_to_geometry(
                msh._cpp_object, fdim, f.astype(np.int32), False
            )
            p = msh.geometry.x[g][:, :3, :]
            a = 0.5 * np.linalg.norm(np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]), axis=1)
            c = p.mean(axis=1)
            s = np.array([a.sum(), *(a[:, None] * c).sum(axis=0)])
        else:
            s = np.zeros(4)
        s = comm.allreduce(s, op=MPI.SUM)
        out[o] = (float(s[0]), s[1:] / s[0] if s[0] > 0 else np.full(3, np.nan))
    return out


def _owned_ext_tag_check(msh, facets, comm):
    """(Q3) owned exterior facets not tagged in {1, 501, 502, 401}; 401 on planes."""
    import dolfinx

    from fem_em_solver.io.mesh import BIRDCAGE_CUT_X0_TAG, BIRDCAGE_CUT_Y0_TAG
    from tests.mesh.test_birdcage_symmetry_quadrant import _facet_nodes

    fdim = msh.topology.dim - 1
    msh.topology.create_connectivity(fdim, fdim + 1)
    n_owned = msh.topology.index_map(fdim).size_local
    ext = np.asarray(dolfinx.mesh.exterior_facet_indices(msh.topology), dtype=np.int32)
    ext = ext[ext < n_owned]
    idx = np.asarray(facets.indices)
    val = np.asarray(facets.values)
    allowed = {1, BIRDCAGE_CUT_X0_TAG, BIRDCAGE_CUT_Y0_TAG, CAVITY_TAG}
    tagged = set(idx[np.isin(val, list(allowed)) & (idx < n_owned)].tolist())
    untagged = len(set(ext.tolist()) - tagged)
    cav = idx[(val == CAVITY_TAG) & (idx < n_owned)]
    pts = _facet_nodes(msh, cav)
    on_plane = 0
    if len(cav):
        on_plane = int(
            np.count_nonzero(np.all(np.abs(pts[:, :, 0]) < PLANE_TOL, axis=1))
            + np.count_nonzero(np.all(np.abs(pts[:, :, 1]) < PLANE_TOL, axis=1))
        )
    n_cav = comm.allreduce(len(cav), op=MPI.SUM)
    return (
        comm.allreduce(untagged, op=MPI.SUM),
        comm.allreduce(on_plane, op=MPI.SUM),
        n_cav,
    )


def _cavity_neighbour_tags(msh, cells, facets):
    """Owned tag-401 facets and the tag of their single adjacent cell."""
    from tests.validation.test_wf12_f_human_hole import _cell_value

    tdim = msh.topology.dim
    f2c = msh.topology.connectivity(tdim - 1, tdim)
    offsets, links = f2c.offsets, f2c.array
    cv = _cell_value(msh, cells)
    n_owned = msh.topology.index_map(tdim - 1).size_local
    idx = np.asarray(facets.indices)[np.asarray(facets.values) == CAVITY_TAG]
    idx = idx[idx < n_owned].astype(np.int32)
    return idx, cv[links[offsets[idx]]]


def _masked_wall(msh, idx, mask, comm):
    """(area, count, min over vertices of min(x, y)) of owned facets idx[mask],
    reduced over ranks; a rank owning none contributes (0, 0, +inf)."""
    import dolfinx

    fdim = msh.topology.dim - 1
    f = idx[mask]
    area, count, mxy = 0.0, 0, math.inf
    if len(f):
        g = dolfinx.cpp.mesh.entities_to_geometry(msh._cpp_object, fdim, f, False)
        p = msh.geometry.x[g][:, :3, :]
        area = float(
            (0.5 * np.linalg.norm(np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]), axis=1)).sum()
        )
        count = int(len(f))
        mxy = float(np.min(np.minimum(p[:, :, 0], p[:, :, 1])))
    return (
        comm.allreduce(area, op=MPI.SUM),
        comm.allreduce(count, op=MPI.SUM),
        comm.allreduce(mxy, op=MPI.MIN),
    )


_CACHE: dict = {}


def _meshes():
    """Full flag-off hole mesh first, then the quarter (`GEO-23`)."""
    if _CACHE:
        return _CACHE
    from fem_em_solver.io.mesh import MeshGenerator
    from tests.validation.test_wf12_f_human_hole import _owned_census

    comm = MPI.COMM_WORLD
    for name, quadrant in (("full", False), ("quarter", True)):
        kwargs = _kwargs()
        if quadrant:
            kwargs["symmetry_quadrant"] = True
        t0 = time.perf_counter()
        msh, cells, facets, diag = MeshGenerator.birdcage_port_domain(
            comm=comm, return_diagnostics=True, **kwargs
        )
        t_build = time.perf_counter() - t0
        tdim = msh.topology.dim
        msh.topology.create_connectivity(tdim - 1, tdim)
        n_cells = int(msh.topology.index_map(tdim).size_global)
        census = _owned_census(msh, cells, comm)
        _CACHE[name] = {
            "msh": msh,
            "cells": cells,
            "facets": facets,
            "diag": diag,
            "kwargs": kwargs,
            "t_build": t_build,
            "n_cells": n_cells,
            "census": census,
        }
        _rank0_print(
            comm,
            f"\n[TH-20 0a] {name}: {n_cells} cells, build {t_build:.1f} s, gmsh mesh "
            f"{diag.get('mesh_wall_time_s', float('nan')):.1f} s, -n {comm.size}; "
            f"census {census}",
        )
    return _CACHE


def _extents(kwargs, port_dy):
    r_e = (
        kwargs["ring_radius"]
        + max(0.5 * kwargs["leg_width"], kwargs["ring_minor_radius"])
        + port_dy
        + kwargs["air_padding"]
    )
    z_e = (
        max(
            0.5 * kwargs["coil_length"],
            0.5 * kwargs["leg_spacing"] + kwargs["ring_minor_radius"],
            0.5 * kwargs["phantom_height"],
        )
        + kwargs["air_padding"]
    )
    return r_e, z_e


# ---------------------------------------------------------------------------
# (Q7) refusals — no meshing, always runs
# ---------------------------------------------------------------------------
@pytest.mark.parametrize(
    "case",
    [
        "solid_ring_gaps",
        "leg_count_14",
        "transverse",
        "ring_and_leg_gaps",
        "azimuth_offsets",
        "c4_congruent",
    ],
)
def test_th20_q7_lh_refusals(case, monkeypatch):
    from fem_em_solver.io import mesh as mesh_mod

    def _no_gmsh(*_a, **_k):  # any build would reach the model builder
        raise AssertionError("refusal reached the gmsh builder")

    monkeypatch.setattr(mesh_mod.MeshGenerator, "_build_birdcage_port_model", _no_gmsh)
    kwargs = _kwargs()
    kwargs["symmetry_quadrant"] = True
    if case == "solid_ring_gaps":
        kwargs["as_hole"] = False
    elif case == "leg_count_14":
        kwargs["leg_count"] = 14
    elif case == "transverse":
        kwargs["ring_sheet_orientation"] = "transverse"
    elif case == "ring_and_leg_gaps":
        kwargs["leg_gap_length"] = 0.01
    elif case == "azimuth_offsets":
        kwargs["leg_azimuth_offsets_rad"] = [0.0] * kwargs["leg_count"]
    elif case == "c4_congruent":
        kwargs["c4_congruent_sheets"] = True
    with pytest.raises(ValueError):
        mesh_mod.MeshGenerator.birdcage_port_domain(comm=MPI.COMM_WORLD, **kwargs)


# ---------------------------------------------------------------------------
# (Q0) the flag-off record
# ---------------------------------------------------------------------------
def test_th20_q0_full_hole_record():
    _skip_unless_on()
    from tests.mesh.helpers import global_cell_tag_set

    full = _meshes()["full"]
    tags = global_cell_tag_set(full["msh"], full["cells"])
    _rank0_print(
        MPI.COMM_WORLD,
        f"[TH-20 0a] (Q0) full cells {full['n_cells']} "
        f"(record {FULL_HOLE_CELL_RECORD}); tags == FULL_TAGS: {tags == FULL_TAGS}",
    )
    assert full["n_cells"] == FULL_HOLE_CELL_RECORD
    assert tags == FULL_TAGS


# ---------------------------------------------------------------------------
# (Q1)–(Q6) the quarter, plus the controls
# ---------------------------------------------------------------------------
def test_th20_q1_to_q6_quarter_identities():
    _skip_unless_on()
    from fem_em_solver.io.mesh import BIRDCAGE_CUT_X0_TAG, BIRDCAGE_CUT_Y0_TAG
    from tests.mesh.helpers import global_cell_tag_set
    from tests.mesh.probe_two_torus_conductor_hole import _exterior_facet_area
    from tests.mesh.test_birdcage_conductor_hole import CAD_ANALYTIC_BAND, EXACT
    from tests.mesh.test_birdcage_leg_offset import SHEET_AREA_BAND
    from tests.mesh.test_birdcage_symmetry_quadrant import (
        EXACT_IDENTITY_RTOL,
        _plane_readings,
    )
    from tests.mesh.test_coil_phantom_conforming import _tag_volume
    from tests.validation.test_birdcage_f_human_rung import _ring_ports
    from tests.validation.test_wf12_f_human_hole import RING_PORT_WALL_RTOL

    comm = MPI.COMM_WORLD
    m = _meshes()
    full, q = m["full"], m["quarter"]
    msh, cells, facets, diag = q["msh"], q["cells"], q["facets"], q["diag"]
    kw = q["kwargs"]
    fails = []

    # ---- price, printed
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0**2
    rss_sum = comm.allreduce(rss, op=MPI.SUM)
    layout = diag.get("symmetry_quadrant_layout")
    _rank0_print(
        comm,
        f"[TH-20 0a] quarter cells {q['n_cells']} (predicted {PREDICTED_CELLS}); gmsh "
        f"mesh {diag['mesh_wall_time_s']:.1f} s (predicted {PREDICTED_MESH_S}); peak RSS "
        f"rank0 {rss:.2f} GiB, summed {rss_sum:.2f} GiB; layout {layout}",
    )

    # ---- (Q1)
    tags = global_cell_tag_set(msh, cells)
    readings = _plane_readings(msh, cells, facets)
    census_sum = sum(q["census"].values())
    _rank0_print(
        comm,
        f"[TH-20 0a] (Q1) tags {sorted(tags)}; outside vertices "
        f"{readings['n_outside']} (min xy {readings['min_xy']:.3e}); census sum "
        f"{census_sum} vs {q['n_cells']}",
    )
    if tags != QUARTER_TAGS:
        fails.append(f"Q1 tags {sorted(tags)} != {sorted(QUARTER_TAGS)}")
    if readings["n_outside"] != 0:
        fails.append(f"Q1 {readings['n_outside']} vertices outside the quadrant")
    if census_sum != q["n_cells"]:
        fails.append(f"Q1 census {census_sum} != {q['n_cells']}")
    if (
        layout is None
        or layout["layout"] != "LH"
        or tuple(layout["kept_ordinals"]) != KEPT_LEG + KEPT_RING
    ):
        fails.append(f"Q1/Q4 layout diagnostics {layout}")

    # ---- (Q2)
    r_e, z_e = _extents(kw, float(diag["port_box_size_m"][1]))
    a_cut = r_e * 2.0 * z_e - 2.0 * (0.5 * kw["leg_width"]) * kw["coil_length"]
    a_outer = 4.0 * r_e * z_e + 2.0 * r_e * r_e
    _rank0_print(
        comm,
        f"[TH-20 0a] r_e {r_e:.9e} z_e {z_e:.9e}; A_cut closed form "
        f"{a_cut:.9e} (brief {A_CUT_PREDICTED:.8e}); A_outer {a_outer:.9e} "
        f"(brief {A_OUTER_PREDICTED:.9e})",
    )
    for cut in (BIRDCAGE_CUT_X0_TAG, BIRDCAGE_CUT_Y0_TAG):
        r = readings[cut]
        dev = abs(r["area"] / a_cut - 1.0)
        _rank0_print(
            comm,
            f"[TH-20 0a] (Q2) cut {cut}: facets geo {r['n_geo']} tagged "
            f"{r['n_tagged']} mismatch {r['mismatch']} flat {r['flat']:.3e} area "
            f"{r['area']:.12e} rel dev {dev:.3e} (band {EXACT_IDENTITY_RTOL:g})",
        )
        if r["mismatch"] != 0 or r["n_tagged"] == 0:
            fails.append(f"Q2 cut {cut} mismatch {r['mismatch']} tagged {r['n_tagged']}")
        if r["flat"] > PLANE_TOL:
            fails.append(f"Q2 cut {cut} flat {r['flat']:.3e}")
        if not dev <= EXACT_IDENTITY_RTOL:
            fails.append(f"Q2 cut {cut} area dev {dev:.3e}")
    mask = (
        readings[BIRDCAGE_CUT_X0_TAG]["on_plane_mask"]
        | readings[BIRDCAGE_CUT_Y0_TAG]["on_plane_mask"]
    )
    adj = set(int(t) for t in np.unique(readings["_ext_cell_tag"][mask]).tolist())
    adj = set().union(*comm.allgather(adj))
    outer_dev = abs(readings["outer_area"] / a_outer - 1.0)
    _rank0_print(
        comm,
        f"[TH-20 0a] (Q2) cells on the planes {sorted(adj)}; outer area "
        f"{readings['outer_area']:.12e} rel dev {outer_dev:.3e}; outer∩cut "
        f"{readings['outer_cut_overlap']}",
    )
    if adj != {2, 3}:
        fails.append(f"Q2 cells on planes {sorted(adj)}")
    if not outer_dev <= EXACT_IDENTITY_RTOL:
        fails.append(f"Q2 outer area dev {outer_dev:.3e}")
    if readings["outer_cut_overlap"] != 0:
        fails.append(f"Q2 outer∩cut {readings['outer_cut_overlap']}")

    # ---- (Q3)
    untagged, cav_on_plane, n_cav = _owned_ext_tag_check(msh, facets, comm)
    a_q = _exterior_facet_area(msh, facets, CAVITY_TAG, comm)
    a_f = _exterior_facet_area(full["msh"], full["facets"], CAVITY_TAG, comm)
    q3 = abs(4.0 * a_q / a_f - 1.0)
    cad_cav_q = diag["conductor_cavity_cad"]["area_m2"]
    _rank0_print(
        comm,
        f"[TH-20 0a] (Q3) facets 401 {n_cav}, 501 "
        f"{readings[BIRDCAGE_CUT_X0_TAG]['n_tagged']}, 502 "
        f"{readings[BIRDCAGE_CUT_Y0_TAG]['n_tagged']}; 401 on planes {cav_on_plane}; "
        f"untagged exterior {untagged}; A401 q {a_q:.9e} f {a_f:.9e} |4q/f-1| "
        f"{q3:.3e} (band {TH20_QUARTER_WALL_AREA_RTOL:g}, predicted 1.498e-4); "
        f"CAD cavity quarter {cad_cav_q:.9e} (predicted "
        f"{A_CAVITY_QUARTER_CAD_PREDICTED:.7e})",
    )
    if cav_on_plane != 0:
        fails.append(f"Q3 {cav_on_plane} tag-401 facets on a plane")
    if untagged != 0:
        fails.append(f"Q3 {untagged} owned exterior facets outside {{1,501,502,401}}")
    if not q3 <= TH20_QUARTER_WALL_AREA_RTOL:
        fails.append(f"Q3 cavity |4q/f-1| {q3:.3e}")
    # (C3) the cut faces counted as cavity wall; asserted under CONTROL.
    c3 = (
        4.0
        * (a_q + readings[BIRDCAGE_CUT_X0_TAG]["area"] + readings[BIRDCAGE_CUT_Y0_TAG]["area"])
        / a_f
        - 1.0
    )
    _rank0_print(
        comm,
        f"[TH-20 0a] (C3) cavity rule without the plane walls: "
        f"4 A_q'/A_f - 1 = {c3:+.3f} (predicted {C3_PREDICTED:+.2f}; asserted "
        f"> {TH20_QUARTER_WALL_AREA_RTOL:g}: {CONTROL})",
    )
    if CONTROL and not c3 > TH20_QUARTER_WALL_AREA_RTOL:
        fails.append(f"C3 {c3:+.3e} not above the Q3 band")

    # ---- (Q3-b) the ring-port wall restricted to the kept tags, both meshes
    kept = np.array(sorted(KEPT_RING_PORT_TAGS))
    iq, nq = _cavity_neighbour_tags(msh, cells, facets)
    ifl, nfl = _cavity_neighbour_tags(full["msh"], full["cells"], full["facets"])
    a_rq, n_rq, _ = _masked_wall(msh, iq, np.isin(nq, kept), comm)
    a_rf, n_rf, guard = _masked_wall(full["msh"], ifl, np.isin(nfl, kept), comm)
    q3b = abs(a_rq / a_rf - 1.0) if a_rf > 0 else math.inf
    ring_all = ((nfl >= 117) & (nfl <= 148)) | ((nfl >= 217) & (nfl <= 248))
    a_rf_all, n_rf_all, _ = _masked_wall(full["msh"], ifl, ring_all, comm)
    _rank0_print(
        comm,
        f"[TH-20 0a] (Q3-b) kept ring-port wall: A_R,q {a_rq:.10e} "
        f"({n_rq} facets) A_R,f {a_rf:.10e} ({n_rf} facets) |q/f-1| {q3b:.3e} "
        f"(band {RING_PORT_WALL_RTOL:g}; predicted 7.8336145294e-4, 772/772); "
        f"guard min(x,y) over full masked vertices {guard:.4e} m (predicted "
        f"≈ 2.03e-2, asserted > 0)",
    )
    _rank0_print(
        comm,
        f"[TH-20 0a] (Q3-b trap, printed) whole-class 4 A_R,q / "
        f"A_R,f(117-148,217-248) - 1 = {4.0 * a_rq / a_rf_all - 1.0:+.3e} "
        f"(4x{n_rq} vs {n_rf_all} facets; predicted -4.65e-5)",
    )
    if n_rq != n_rf or n_rq == 0:
        fails.append(f"Q3-b facet counts {n_rq} != {n_rf}")
    if not q3b <= RING_PORT_WALL_RTOL:
        fails.append(f"Q3-b |A_R,q/A_R,f - 1| {q3b:.3e}")
    if not guard > 0.0:
        fails.append(f"Q3-b guard min(x,y) {guard:.3e} not > 0")

    # ---- (C4) planted missing face: the quarter's facets bordering tag 117
    a_117, n_117, _ = _masked_wall(msh, iq, nq == 117, comm)
    d_n = -n_117
    d_a = a_117 / a_rq if a_rq > 0 else 0.0
    planted_q3 = 4.0 * (a_q - a_117) / a_f - 1.0
    _rank0_print(
        comm,
        f"[TH-20 0a] (C4) planted missing tag-117 face: dN_R {d_n} "
        f"(predicted -{C4_PLANTED_FACETS}), |dA_R/A_R| {d_a:.4e} (bar "
        f"{1e6 * RING_PORT_WALL_RTOL:g}, predicted 6.731e-2), A_117 {a_117:.7e} m^2; "
        f"asserted={CONTROL}",
    )
    _rank0_print(
        comm,
        f"[TH-20 0a] (C4 under Q3 alone, printed, never asserted) "
        f"4 q'/f - 1 = {planted_q3:+.4e} (predicted -8.653e-4; band "
        f"{TH20_QUARTER_WALL_AREA_RTOL:g})",
    )
    if CONTROL and not (d_n == -C4_PLANTED_FACETS and d_a >= 1e6 * RING_PORT_WALL_RTOL):
        fails.append(f"C4 dN {d_n} dA/A {d_a:.3e}")

    # ---- (Q4) volumes
    cad_q = diag["cad_mass_by_group"]
    cad_f = full["diag"]["cad_mass_by_group"]
    names = {}
    for o in KEPT_LEG:
        names[100 + o] = (f"port_P{o}", CAD_LEG_BOX)
    for o in KEPT_RING:
        names[100 + o] = (f"port_P{o}_lower", CAD_RING_LOWER)
        names[200 + o] = (f"port_P{o}_upper", CAD_RING_UPPER)
    for tag, (name, cad_imported) in sorted(names.items()):
        vq = _tag_volume(msh, cells, tag, comm)
        vf = _tag_volume(full["msh"], full["cells"], tag, comm)
        cq = float(cad_q.get(name, float("nan")))
        d_cad = abs(vq / cq - 1.0)
        d_full = abs(vq / vf - 1.0)
        d_imp = abs(cq / cad_imported - 1.0)
        _rank0_print(
            comm,
            f"[TH-20 0a] (Q4) {name}: V_q {vq:.12e} CAD_q {cq:.12e} V_f "
            f"{vf:.12e} CAD_f {float(cad_f.get(name, float('nan'))):.12e}; dev cad "
            f"{d_cad:.3e} full {d_full:.3e} imported {d_imp:.3e} (band {EXACT:g})",
        )
        if not (d_cad <= EXACT and d_full <= EXACT and d_imp <= EXACT):
            fails.append(f"Q4 {name} dev cad {d_cad:.3e} full {d_full:.3e} imp {d_imp:.3e}")

    # ---- (Q4) sheets + (Q5) the D2 map
    rl = diag["ring_port_layout"]
    sheet_analytic = float(rl["ring_port_sheet_longitudinal_area_m2"])
    sq = _sheet_centroids(msh, cells, KEPT_RING, comm)
    ring_ports = _ring_ports()
    sf = _sheet_centroids(full["msh"], full["cells"], ring_ports, comm)
    for o in KEPT_RING:
        dev = abs(sq[o][0] / sheet_analytic - 1.0)
        _rank0_print(
            comm,
            f"[TH-20 0a] (Q4) sheet P{o}: area {sq[o][0]:.12e} vs chord·w "
            f"{sheet_analytic:.12e} dev {dev:.3e} (band {SHEET_AREA_BAND:g}); "
            f"centroid {sq[o][1]}",
        )
        if not dev <= SHEET_AREA_BAND:
            fails.append(f"Q4 sheet P{o} dev {dev:.3e}")
    full_ord = np.array(ring_ports)
    full_c = np.array([sf[o][1] for o in ring_ports])
    group_ops = {
        "e": np.array([1.0, 1.0, 1.0]),
        "sy": np.array([1.0, -1.0, 1.0]),
        "sx": np.array([-1.0, 1.0, 1.0]),
        "sxsy": np.array([-1.0, -1.0, 1.0]),
    }
    closed = {
        "e": lambda j: j,
        "sy": lambda j: 15 - j,
        "sx": lambda j: (7 - j) % 16,
        "sxsy": lambda j: (j + 8) % 16,
    }
    radius = float(kw["ring_radius"])
    tol5 = EXACT * radius
    images = []
    worst5 = 0.0
    for g, op in group_ops.items():
        for o in KEPT_RING:
            img = sq[o][1] * op
            d = np.linalg.norm(full_c - img, axis=1)
            k = int(np.argmin(d))
            worst5 = max(worst5, float(d[k]))
            target = int(full_ord[k])
            base, j = (17, o - 17) if o < 33 else (33, o - 33)
            expect = base + closed[g](j)
            images.append(target)
            if d[k] > tol5:
                fails.append(f"Q5 {g} P{o}: miss {d[k]:.3e} > {tol5:.3e}")
            if target != expect:
                fails.append(f"Q5 {g} P{o} -> P{target}, closed form P{expect}")
    _rank0_print(
        comm,
        f"[TH-20 0a] (Q5) worst image miss {worst5:.3e} m (tol {tol5:.3e}); "
        f"images cover {len(set(images))} of {len(ring_ports)} ring ordinals",
    )
    if sorted(images) != sorted(ring_ports):
        fails.append(f"Q5 not a bijection: {sorted(images)}")

    # ---- (C2) the mirror misplaced by π/32: φ -> 2α - φ with α = π/32
    alpha = math.pi / 32.0
    ceiling = 2.0 * radius * math.sin(math.pi / 32.0)
    c2_worst = 0.0
    for o in KEPT_RING:
        c = sq[o][1]
        rot = np.array(
            [
                [math.cos(2 * alpha), math.sin(2 * alpha)],
                [math.sin(2 * alpha), -math.cos(2 * alpha)],
            ]
        )
        img = np.array([*(rot @ c[:2]), c[2]])
        miss = float(np.min(np.linalg.norm(full_c - img, axis=1)))
        c2_worst = max(c2_worst, abs(miss / ceiling - 1.0))
    _rank0_print(
        comm,
        f"[TH-20 0a] (C2) misplaced mirror: worst |miss/(2R sin(π/32)) - 1| "
        f"{c2_worst:.3e} (ceiling {ceiling:.7e} m, predicted 2.9405142e-2; "
        f"asserted={CONTROL})",
    )
    if CONTROL and not c2_worst <= EXACT:
        fails.append(f"C2 misplaced-mirror miss dev {c2_worst:.3e}")

    # ---- (Q6) phantom, air printed
    ph = 4.0 * cad_q["phantom"] / cad_f["phantom"] - 1.0
    air = 4.0 * cad_q.get("air", float("nan"))
    _rank0_print(
        comm,
        f"[TH-20 0a] (Q6) phantom 4 CAD_q/CAD_f - 1 = {ph:.3e} (band "
        f"{CAD_ANALYTIC_BAND:g}); air CAD x4 {air:.9e} vs full "
        f"{cad_f.get('air', float('nan')):.9e}",
    )
    if not abs(ph) <= CAD_ANALYTIC_BAND:
        fails.append(f"Q6 phantom {ph:.3e}")

    # ---- (C1) the quarter checker on the flag-off full mesh
    rf = _plane_readings(full["msh"], full["cells"], full["facets"])
    c1 = (
        rf["n_outside"],
        rf[BIRDCAGE_CUT_X0_TAG]["n_geo"],
        rf[BIRDCAGE_CUT_Y0_TAG]["n_geo"],
        rf[BIRDCAGE_CUT_X0_TAG]["area"],
        rf[BIRDCAGE_CUT_Y0_TAG]["area"],
    )
    _rank0_print(
        comm,
        f"[TH-20 0a] (C1) full mesh: vertices outside {c1[0]} (predicted > 0); "
        f"plane facets {c1[1]} / {c1[2]} (predicted 0); cut area {c1[3]:.3e} / "
        f"{c1[4]:.3e} (predicted 0); asserted={CONTROL}",
    )
    if CONTROL and not (c1[0] > 0 and c1[1] == 0 and c1[2] == 0 and c1[3] == 0.0 and c1[4] == 0.0):
        fails.append(f"C1 {c1}")

    _rank0_print(comm, f"[TH-20 0a] failures: {fails}")
    assert not fails, "; ".join(fails)
