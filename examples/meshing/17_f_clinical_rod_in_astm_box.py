"""Example (`EX-137`): the kept 1 mm x 30 mm rod inside F-clinical-1.5T's ASTM F2182 gel box.

`GEO-36` step 8 (`tests/validation/test_geo36_rod_f_clinical.py`, §9 item 530,
commit `c07031df`) builds the clinical-scale fixture of `mesh:15` (32 legs,
64 ring-gap sheets, closed PEC can, conductors as holes, the 0.42 x 0.09 x
0.65 m gel box) with the kept, uncoated insulated-geometry rod
(``implant_rod``, r = 1 mm, l = 30 mm, tags 403 / 402) at x = Lx/4 = 0.105 m
in C2's coil frame, and gates (G0)-(G10). This example builds that same
fixture from the gate's own ``_kwargs()`` (never retyped, no centre typed),
puts it in ParaView, recomputes the gate's identities on its own mesh through
imported helpers and bands (owned facets / cells only, reduced; the only UFL
``dS`` is the gate's ``_interior_facet_check``), and adds one reading the gate
does not make: the rod's placement **read off the mesh vertices** (Cap).

Example-side anchors: (F) / (R3) as `mesh:15` with the rod inside; (Cap) every
owned tag-403 vertex lies within the capsule at ``_rod_centre()``.
Negative controls, asserted: (N1r), (N2r), (cB) a mis-placed rod is never
built, (cR) the capsule test against step 6's centre puts no vertex inside,
(cF) one box face missing.

Scope: a mesh, not a solve. No PEC or sigma condition, tip field, SAR,
heating; no ASTM F2182 placement or compliance claim (the +x placement is a
modelling choice, review ruling 7); the ``-n 2`` width is not claimed.

Run it through the example runner::

    ./scripts/run_examples.sh -e mesh:17 -n 4 -t 600
"""

from __future__ import annotations

import re
import resource
import sys
import time
from pathlib import Path

import numpy as np
import pytest
from mpi4py import MPI

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from fem_em_solver.io.mesh import (  # noqa: E402
    BIRDCAGE_IMPLANT_CELL_TAG,
    BIRDCAGE_IMPLANT_SURFACE_TAG,
    MeshGenerator,
    _interface_facet_tags,
)
from fem_em_solver.io.paraview_utils import (  # noqa: E402
    adopt_host_ownership,
    write_xdmf_with_tags,
)
from fem_em_solver.post.setup_figure import write_setup_figure  # noqa: E402

import tests.validation.test_geo36_rod_f_clinical as g36  # noqa: E402
from tests.mesh.test_birdcage_conductor_hole import EXACT  # noqa: E402
from tests.mesh.test_birdcage_leg_gaps import TERMINAL_AREA_BAND  # noqa: E402
from tests.mesh.test_birdcage_port_sheets import (  # noqa: E402
    PORT_LOWER,
    PORT_UPPER,
    SHEET_IFACE,
)
from tests.validation.gmsh_stop_rule import (  # noqa: E402
    arm_gmsh_stop_rule,
    gmsh_stop_armed,
)
from tests.validation.test_geo38_f_clinical_build import (  # noqa: E402
    LEG_COUNT,
    _area,
    _facet_coords,
    _owned,
    _ring_ports,
)
from tests.validation.test_ops61_exact_strip_ports import STRIP_BAND  # noqa: E402
from tests.validation.test_wf12_f_human_hole import (  # noqa: E402
    _census_checks,
    _owned_census,
)

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"  # committed, EX-57
BASENAME = "meshing_17_f_clinical_rod_in_astm_box"

# The box interface facet tag: unused by the 64 sheet tags (SHEET_IFACE + i)
# and by the rod wall (BIRDCAGE_IMPLANT_SURFACE_TAG = 402).
BOX_IFACE = 400

# Records, example-local, each with its log line (`GEO-36` step 8 W1, -n 4).
W1_LOG = "20261019T200702Z_GEO-36-rodbox-W1-n4.log"
REC = {
    "tag2_cells": 1_063_136,  # W1_LOG:12004
    "tag3_cells": 85_122,  # W1_LOG:12004
    "tag403_cells": 14_834,  # W1_LOG:12004
    "n402": 4_954,  # owned 402 facets = interior = 3|403 pairs, W1_LOG:12011
    "n_gmsh_total": 2_882,  # W1_LOG:12014
}
CELLS_W1 = 1_217_025  # W1_LOG:12002, PRINTED beside this run's count, never asserted
U2_W1 = 8_327_394  # W1_LOG:12017, cited, never recomputed (no N1curl space here)


def _cell_diameter_function(mesh):
    """Per-cell max vertex-to-vertex distance as a DG0 ``cell_diameter`` field."""
    import dolfinx.cpp.mesh as cpp_mesh
    from dolfinx import fem

    tdim = mesh.topology.dim
    imap = mesh.topology.index_map(tdim)
    n_local = imap.size_local + imap.num_ghosts
    cells = np.arange(n_local, dtype=np.int32)
    h = cpp_mesh.h(mesh._cpp_object, tdim, cells)

    v0 = fem.functionspace(mesh, ("DG", 0))
    diam = fem.Function(v0, name="cell_diameter")
    cell_dofs = v0.dofmap.list.reshape(-1)
    diam.x.array[cell_dofs[cells]] = h
    diam.x.scatter_forward()
    return diam


def _merge_facet_tags(msh, *tag_sets):
    """One sorted, unique facet MeshTags; no facet may carry two tags."""
    import dolfinx

    idx = np.concatenate([np.asarray(t.indices, dtype=np.int32) for t in tag_sets])
    val = np.concatenate([np.asarray(t.values, dtype=np.int32) for t in tag_sets])
    uniq, first = np.unique(idx, return_index=True)
    assert len(uniq) == len(idx), "a facet carries two tags (sheet / box interface / rod wall)"
    return dolfinx.mesh.meshtags(msh, msh.topology.dim - 1, uniq, val[first])


def _capsule_distance(P, centre, axis, half):
    """Distance of each row of ``P`` to the segment ``centre`` +- half * axis."""
    d = P - np.asarray(centre, dtype=float)
    t = d @ axis
    perp = d - t[:, None] * axis
    ax = np.maximum(np.abs(t) - half, 0.0)
    return np.sqrt((perp**2).sum(axis=1) + ax**2)


def main() -> None:
    import dolfinx

    comm = MPI.COMM_WORLD
    started = time.perf_counter()

    kw = g36._kwargs()
    ring_ports = _ring_ports()
    box = kw["phantom_box"]
    Lx, Ly, Lz = box["size_m"]
    off = box["table_offset_m"]
    L3 = Lx * Ly * Lz
    V = g36._rod_volume()
    A = g36._rod_area()
    r_rod = g36.ROD_RADIUS
    l_rod = g36.ROD_LENGTH
    centre = g36._rod_centre()
    axis = np.asarray(g36.ROD["axis"], dtype=float)
    axis = axis / np.linalg.norm(axis)
    half = 0.5 * l_rod

    # ---- structural asserts at start-up -----------------------------------
    assert comm.size == 4, f"EX-137 is a -n 4 example, comm.size = {comm.size}"
    assert np.issubdtype(dolfinx.default_scalar_type, np.floating), (
        "the meshing group runs the real build"
    )
    assert g36.ROD_OFFSET_M is None, g36.ROD_OFFSET_M
    assert kw["implant_rod"]["keep_cells"] is True
    assert "coat_radius" not in kw["implant_rod"]
    assert kw["implant_rod"]["centre"] == centre
    assert kw["as_hole"] is True
    assert box["table_offset_m"] == 0.0
    assert len(ring_ports) == 2 * LEG_COUNT, len(ring_ports)

    if comm.rank == 0:
        print("=" * 72)
        print("EX-137 - F-clinical-1.5T: ASTM F2182 gel box with the kept 1 x 30 mm rod")
        print("=" * 72)
        print(
            f"\n[geometry] legs {LEG_COUNT}  ring ports {len(ring_ports)}  "
            f"ring_radius {kw['ring_radius']} m  shield_radius {kw['shield_radius']} m"
            f"\n[box]      size_m (Lx, Ly, Lz) = ({Lx}, {Ly}, {Lz}) m  table_offset {off} m"
            f"\n[rod]      centre {centre} (Lx/4 toward +x in C2's coil frame), axis "
            f"{g36.ROD['axis']}, r {r_rod} m, l {l_rod} m, kept, uncoated"
            f"\n[rod]      V_rod {V:.12e} m^3, A_rod {A:.12e} m^2, L^3 {L3:.15e} m^3"
            "\n[scope]    mesh only: no solve, no PEC / sigma condition, no tip field, no SAR; "
            "no ASTM placement claim",
            flush=True,
        )

    # ---- (cB): a mis-placed rod is never built (pure, every rank) ----------
    R_cyl = float(np.hypot(0.5 * Lx, 0.5 * Ly + abs(off)))
    breach = abs(g36.Y_MISREAD_CENTRE[1] - off) + r_rod - 0.5 * Ly
    cb_msg = None
    try:
        MeshGenerator._normalise_implant_rod(
            dict(g36.ROD, keep_cells=True, centre=g36.Y_MISREAD_CENTRE),
            phantom_radius=R_cyl,
            phantom_height=Lz,
            as_hole=True,
            phantom_box=dict(box),
        )
    except ValueError as exc:
        cb_msg = str(exc)
    cb_refused = cb_msg is not None and re.search(
        "implant_rod is not strictly inside the phantom_box", cb_msg
    ) is not None
    cb_ok = MeshGenerator._normalise_implant_rod(
        dict(g36.ROD, keep_cells=True, centre=centre),
        phantom_radius=R_cyl,
        phantom_height=Lz,
        as_hole=True,
        phantom_box=dict(box),
    )
    cb_accepts = cb_ok["keep_cells"] is True and cb_ok["centre"] == centre
    if comm.rank == 0:
        print(
            f"\n[EX-137] (cB) y-misread centre {g36.Y_MISREAD_CENTRE}: refused={cb_refused}; "
            f"breach |c_y - off| + r - Ly/2 = {breach:.6e} m (W1 6.100000e-02 m, {W1_LOG}:68); "
            f"the rod's own centre accepted={cb_accepts}"
            f"\n[EX-137] (cB) message: {cb_msg}",
            flush=True,
        )
    assert cb_refused, cb_msg
    assert cb_accepts, cb_ok
    assert breach > 0.0, breach

    # ---- build (timer captured before the figure and XDMF) -----------------
    armed = gmsh_stop_armed()
    if comm.rank == 0:
        print(f"\n[EX-137] gmsh stop rule armed: {armed}", flush=True)
    t0 = time.perf_counter()
    if armed:
        with pytest.MonkeyPatch.context() as mp:
            arm_gmsh_stop_rule(mp, comm, "EX-137")
            msh, cells, facets, diag = MeshGenerator.birdcage_port_domain(
                comm=comm, return_diagnostics=True, **kw
            )
    else:
        msh, cells, facets, diag = MeshGenerator.birdcage_port_domain(
            comm=comm, return_diagnostics=True, **kw
        )
    t_build = time.perf_counter() - t0
    if comm.rank == 0:
        print(f"\n[EX-137] build done {t_build:.1f} s at -n {comm.size}", flush=True)

    tdim = msh.topology.dim
    fdim = tdim - 1
    n_cells = int(msh.topology.index_map(tdim).size_global)
    # (G0), before any facet reading
    if comm.rank == 0:
        print(
            f"[EX-137] (G0) {n_cells} cells vs guard {g36.F_CLINICAL_ROD_CELL_GUARD}",
            flush=True,
        )
    assert n_cells < g36.F_CLINICAL_ROD_CELL_GUARD, n_cells

    # ---- setup figure (flag-gated; after the build timer) ------------------
    t_fig0 = time.perf_counter()
    write_setup_figure(
        msh,
        cells,
        FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names={
            2: "air",
            3: "phantom (ASTM F2182 gel box)",
            BIRDCAGE_IMPLANT_CELL_TAG: "implant rod (kept, uncoated)",
            **{PORT_LOWER + i: f"leg L{i} port box" for i in range(1, LEG_COUNT + 1)},
            **{PORT_LOWER + i: f"ring port P{i} lower" for i in ring_ports},
            **{PORT_UPPER + i: f"ring port P{i} upper" for i in ring_ports},
        },
        hide_tags=(2,),
        translucent_tags=(3,),
        slice_normal=(0.0, 1.0, 0.0),
        slice_origin=(0.0, 0.0, 0.0),
        show_edges_on_slice=False,
        title=(
            "mesh:17 -- F-clinical coil, ASTM F2182 gel box, kept 1 x 30 mm rod "
            "at x = Lx/4 (GEO-36 step 8)"
        ),
        comm=comm,
    )
    t_fig = time.perf_counter() - t_fig0

    # ---- common mesh bookkeeping -------------------------------------------
    msh.topology.create_entities(fdim)
    msh.topology.create_connectivity(fdim, tdim)
    msh.topology.create_entity_permutations()  # every rank, before any dS
    census = _owned_census(msh, cells, comm)

    R = kw["ring_radius"]
    r_leg = 0.5 * kw["leg_width"]
    layout = diag["ring_port_layout"]
    chord = float(layout["ring_port_gap_chord_m"])
    width = float(layout["ring_port_box_width_m"])
    sheet_ref = chord * width
    n_gmsh = {int(k): int(v) for k, v in diag["ring_sheet_gmsh_triangles"].items()}

    # ---- (G1), (G2): CAD ----------------------------------------------------
    rod = diag["implant_rod_cad"]
    cad = dict(diag["cad_mass_by_group"])
    rod_mass = float(rod["rod_mass_m3"])
    kept = float(rod["kept_mass_m3"])
    g1_mass = abs(rod_mass / V - 1.0)
    g1_kept = abs(kept / V - 1.0)
    g1_area = abs(float(rod["wall_area_m2"]) / A - 1.0)
    g2_403 = abs(cad["implant"] / V - 1.0)
    g2_3 = abs(cad["phantom"] / (L3 - V) - 1.0)
    g2_sum = abs((cad["phantom"] + cad["implant"]) / L3 - 1.0)
    n1_miss = abs(cad["phantom"] / L3 - 1.0)
    n1_pred = V / (L3 - V)

    # ---- (G3), (R3), (Cap): owned tag-3 / tag-403 tets, reduced -------------
    cv = np.asarray(cells.values)
    ci = np.asarray(cells.indices)
    n_owned_cells = msh.topology.index_map(tdim).size_local

    def _tets(tag):
        sel = ci[(cv == tag) & (ci < n_owned_cells)]
        if len(sel) == 0:
            return np.zeros((0, 4, 3))
        g = dolfinx.cpp.mesh.entities_to_geometry(
            msh._cpp_object, tdim, np.asarray(sel, dtype=np.int32), False
        )
        return msh.geometry.x[g[:, :4]]

    def _vol(X):
        if len(X) == 0:
            return 0.0
        return float(
            np.abs(
                np.einsum(
                    "ij,ij->i",
                    X[:, 1] - X[:, 0],
                    np.cross(X[:, 2] - X[:, 0], X[:, 3] - X[:, 0]),
                )
            ).sum()
            / 6.0
        )

    X3 = _tets(g36.PHANTOM_CELL_TAG)
    X403 = _tets(BIRDCAGE_IMPLANT_CELL_TAG)
    vol3 = comm.allreduce(_vol(X3), op=MPI.SUM)
    vol403 = comm.allreduce(_vol(X403), op=MPI.SUM)
    if len(X3):
        P3 = X3.reshape(-1, 3)
        lo_local, hi_local = P3.min(axis=0), P3.max(axis=0)
        r3_local = float(np.hypot(P3[:, 0], P3[:, 1]).max(initial=-np.inf))
    else:
        lo_local, hi_local = np.full(3, np.inf), np.full(3, -np.inf)
        r3_local = -np.inf
    lo = np.array([comm.allreduce(float(x), op=MPI.MIN) for x in lo_local])
    hi = np.array([comm.allreduce(float(x), op=MPI.MAX) for x in hi_local])
    r_box_mesh = comm.allreduce(r3_local, op=MPI.MAX)
    g3_vol = abs((vol3 + vol403) / L3 - 1.0)
    lo_ref = np.array([-0.5 * Lx, -0.5 * Ly + off, -0.5 * Lz])
    hi_ref = np.array([0.5 * Lx, 0.5 * Ly + off, 0.5 * Lz])
    g3_ext = float(
        max(np.max(np.abs(lo - lo_ref) / np.abs(lo_ref)), np.max(np.abs(hi - hi_ref) / np.abs(hi_ref)))
    )
    r_box_cad = float(diag["phantom_box_circumradius_m"])
    r3_dev = abs(r_box_mesh / r_box_cad - 1.0)

    # (Cap): every owned tag-403 vertex against the capsule at _rod_centre().
    P403 = X403.reshape(-1, 3)
    if len(P403):
        dcap = _capsule_distance(P403, centre, axis, half) / r_rod - 1.0
        k = int(np.argmax(dcap))
        cap_local = (float(dcap[k]), P403[k].tolist())
        dwrong = _capsule_distance(P403, g36.ROD["centre"], axis, half) / r_rod - 1.0
        cr_local = float(dwrong.min())
        x403_lo = P403.min(axis=0)
        x403_hi = P403.max(axis=0)
    else:
        cap_local = (-np.inf, [float("nan")] * 3)
        cr_local = np.inf
        x403_lo, x403_hi = np.full(3, np.inf), np.full(3, -np.inf)
    cap_all = comm.allgather(cap_local)
    cap_max, cap_vertex = max(cap_all, key=lambda t: t[0])
    cr_min = comm.allreduce(cr_local, op=MPI.MIN)
    v403_lo = np.array([comm.allreduce(float(x), op=MPI.MIN) for x in x403_lo])
    v403_hi = np.array([comm.allreduce(float(x), op=MPI.MAX) for x in x403_hi])
    wall_clearance = 0.5 * Lx - float(v403_hi[0])

    # Rod wall (402) vertices, owned facets only.
    fi = np.asarray(facets.indices)
    fv = np.asarray(facets.values)
    w_owned = _owned(msh, fi[fv == BIRDCAGE_IMPLANT_SURFACE_TAG])
    wxyz = _facet_coords(msh, w_owned).reshape(-1, 3)
    if len(wxyz):
        dw = _capsule_distance(wxyz, centre, axis, half) / r_rod - 1.0
        tw = (wxyz - np.asarray(centre)) @ axis
        lateral = np.abs(tw) < half * (1.0 - 1e-6)
        wall_all_local = float(np.abs(dw).max(initial=-np.inf))
        wall_lat_local = float(np.abs(dw[lateral]).max(initial=-np.inf))
    else:
        wall_all_local = wall_lat_local = -np.inf
    wall_all = comm.allreduce(wall_all_local, op=MPI.MAX)
    wall_lat = comm.allreduce(wall_lat_local, op=MPI.MAX)

    # ---- (G4), (N2r) --------------------------------------------------------
    cad_sum = float(sum(cad.values()))
    m_cond = float(diag["conductor_cad_mass_m3"])
    m_air = float(diag["air_domain_cad_mass_m3"])
    g4 = abs((cad_sum + m_cond) / m_air - 1.0)
    n2_miss = abs((cad_sum + m_cond + rod_mass) / m_air - 1.0)
    n4_cond = abs(m_cond / g36.C2_CONDUCTOR_CAD_MASS_M3 - 1.0)

    # ---- (G5): the 402 interface (the gate's helper, dS licensed) ----------
    n402, n402_ext, n402_int, n402_pair, area402, moments = g36._interior_facet_check(
        msh, cells, facets, BIRDCAGE_IMPLANT_SURFACE_TAG, comm
    )
    g5_ratio = area402 / A

    # ---- (G6), (G7) ---------------------------------------------------------
    c = rod["wall_centroid_m"]
    g6 = max(abs(c[i] - centre[i]) for i in range(3))
    census_fails = _census_checks(census, n_cells, ring_ports)

    # ---- (G8): ring sheets, owned-facet numpy sums only --------------------
    sheet_tags = _interface_facet_tags(
        msh, cells, {SHEET_IFACE + i: (PORT_LOWER + i, PORT_UPPER + i) for i in ring_ports}
    )
    sv = np.asarray(sheet_tags.values)
    si = np.asarray(sheet_tags.indices)
    local = {i: _facet_coords(msh, _owned(msh, si[sv == SHEET_IFACE + i])) for i in ring_ports}
    gathered = comm.allgather(local)
    sheet = {i: np.concatenate([g[i] for g in gathered]) for i in ring_ports}
    n_owned = {i: len(sheet[i]) for i in ring_ports}
    g8_bad = {i: (n_owned[i], n_gmsh.get(i)) for i in ring_ports if n_owned[i] != n_gmsh.get(i)}
    g8_area = {i: abs(float(_area(sheet[i]).sum()) / sheet_ref - 1.0) for i in ring_ports}

    # ---- (F): the box surface, rod inside ----------------------------------
    assert BOX_IFACE not in {SHEET_IFACE + i for i in ring_ports}
    assert BOX_IFACE != BIRDCAGE_IMPLANT_SURFACE_TAG
    box_tags = _interface_facet_tags(msh, cells, {BOX_IFACE: (3, 2)})
    bi = np.asarray(box_tags.indices)
    bf = _owned(msh, bi[np.asarray(box_tags.values) == BOX_IFACE])
    bxyz = _facet_coords(msh, bf)
    barea = _area(bxyz)
    H = (0.5 * Lx, 0.5 * Ly, 0.5 * Lz)
    C = (0.0, off, 0.0)
    face_area = {}
    n_claimed = np.zeros(len(bf), dtype=np.int64)
    for kk in range(3):
        for s in (-1, 1):
            on = np.all(np.abs((bxyz[:, :, kk] - C[kk]) - s * H[kk]) <= STRIP_BAND * H[kk], axis=1)
            n_claimed += on
            face_area[(kk, s)] = float(barea[on].sum())
    keys = sorted(face_area)
    loc = np.array(
        [face_area[kk] for kk in keys]
        + [float(barea.sum()), float(len(bf)), float(np.sum(n_claimed == 0)),
           float(np.sum(n_claimed > 1))]
    )
    tot = np.zeros_like(loc)
    comm.Allreduce(loc, tot, op=MPI.SUM)
    face_g = {kk: float(tot[j]) for j, kk in enumerate(keys)}
    box_area_owned, n_box_owned, n_unclass, n_multi = (
        float(tot[6]), int(tot[7]), int(tot[8]), int(tot[9]))
    surf_ref = 2.0 * (Lx * Ly + Ly * Lz + Lx * Lz)
    f_dev = abs(box_area_owned / surf_ref - 1.0)
    pair_ref = {0: 2.0 * Ly * Lz, 1: 2.0 * Lx * Lz, 2: 2.0 * Lx * Ly}
    pair_dev = {kk: abs((face_g[(kk, -1)] + face_g[(kk, 1)]) / pair_ref[kk] - 1.0) for kk in range(3)}
    cf_area = box_area_owned - face_g[(1, 1)]
    cf_dev = abs(cf_area / surf_ref - 1.0)

    rss = comm.allreduce(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0**2,
                         op=MPI.SUM)

    # ---- ParaView output ----------------------------------------------------
    t_x0 = time.perf_counter()
    wall_idx = fi[fv == BIRDCAGE_IMPLANT_SURFACE_TAG]
    wall_tags = dolfinx.mesh.meshtags(
        msh, fdim, np.asarray(wall_idx, dtype=np.int32),
        np.full(len(wall_idx), BIRDCAGE_IMPLANT_SURFACE_TAG, dtype=np.int32))
    merged = _merge_facet_tags(msh, sheet_tags, box_tags, wall_tags)
    OUTPUT_DIR.mkdir(exist_ok=True)
    diam = _cell_diameter_function(msh)
    path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_combined",
        msh,
        cells,
        {"cell_diameter": diam},
        comm=comm,
        facet_tags=merged,
    )
    adopt_host_ownership(OUTPUT_DIR, comm=comm)
    t_xdmf = time.perf_counter() - t_x0

    # ---- print everything, then assert (all values are global) --------------
    if comm.rank == 0:
        wa = max(g8_area, key=g8_area.get)
        low, high = TERMINAL_AREA_BAND
        print(
            f"\n[EX-137] -n {comm.size}: {n_cells} cells (PRINTED, not asserted; W1 {CELLS_W1}, "
            f"{W1_LOG}:12002; +{n_cells - CELLS_W1}); build {t_build:.1f} s (gmsh "
            f"{diag.get('mesh_wall_time_s', float('nan')):.2f} s); figure {t_fig:.1f} s; XDMF "
            f"{t_xdmf:.1f} s; summed peak RSS {rss:.3f} GiB (W1 2.688)"
            f"\n[EX-137] census {census}"
            f"\n[EX-137] (G1) rod_mass {rod_mass:.12e} / kept {kept:.12e} vs V_rod {V:.12e}: "
            f"rel {g1_mass:.3e} / {g1_kept:.3e} (W1 4.441e-16, {W1_LOG}:12005); OCC wall "
            f"{float(rod['wall_area_m2']):.12e} vs A_rod {A:.12e}: rel {g1_area:.3e} "
            f"(band {STRIP_BAND:.0e})"
            f"\n[EX-137] (G2) CAD(403)/V {g2_403:.3e}; CAD(3) {cad['phantom']:.15e} vs L^3 - V "
            f"{L3 - V:.15e}: rel {g2_3:.3e}; CAD(3)+CAD(403) vs L^3 {L3:.15e}: rel {g2_sum:.3e} "
            f"({W1_LOG}:12006)"
            f"\n[EX-137] (N1r) CAD(3) vs C2's L^3 misses by {n1_miss:.6e} (predicted "
            f"{n1_pred:.6e}; W1 4.006372e-06, {W1_LOG}:12007; {n1_miss / STRIP_BAND:.3e}x "
            f"STRIP_BAND)"
            f"\n[EX-137] (G3) meshed tag-3 {vol3:.15e} + tag-403 {vol403:.15e} vs L^3: rel "
            f"{g3_vol:.3e}; tag-3 extents lo {lo.tolist()} hi {hi.tolist()}: worst rel "
            f"{g3_ext:.3e} (W1 2.220e-16 / 0, {W1_LOG}:12008)"
            f"\n[EX-137] (G4) CAD groups {cad_sum:.15e} + conductor {m_cond:.15e} vs can "
            f"{m_air:.15e}: rel {g4:.3e} (W1 6.457e-11, {W1_LOG}:12009)"
            f"\n[EX-137] (N2r) pre-repair reading (+ V_rod) misses by {n2_miss:.6e} (predicted "
            f"2.366e-07; W1 2.365686e-07, {W1_LOG}:12010; {n2_miss / STRIP_BAND:.3e}x STRIP_BAND)"
            f"\n[EX-137] (G5) 402 owned {n402}, exterior {n402_ext}, interior {n402_int}, "
            f"3|403 pairs {n402_pair}; meshed dS(402)/A_rod {g5_ratio:.6f} (band "
            f"{TERMINAL_AREA_BAND}; W1 0.996243, {W1_LOG}:12011)"
            f"\n[EX-137] (G6) OCC wall centroid {c} vs {centre}: max dev {g6:.3e} "
            f"(bound {EXACT * l_rod:.1e}); meshed "
            f"{tuple(m / area402 for m in moments) if area402 else None}"
            f"\n[EX-137] (G7) census failures {census_fails or 'none'}; implant_spheroid_cad "
            f"{diag.get('implant_spheroid_cad')}"
            f"\n[EX-137] (G8) owned sheet facets == N_gmsh on {len(ring_ports) - len(g8_bad)}"
            f"/{len(ring_ports)}; mismatches {g8_bad or 'none'}; N_gmsh "
            f"{min(n_gmsh.values())}..{max(n_gmsh.values())}, total {sum(n_gmsh.values())} "
            f"(record {REC['n_gmsh_total']}, {W1_LOG}:12014); worst area P{wa} "
            f"{g8_area[wa]:.3e}"
            f"\n[EX-137] REC: tag-2 {census.get(2)} (record {REC['tag2_cells']}), tag-3 "
            f"{census.get(3)} (record {REC['tag3_cells']}), tag-403 "
            f"{census.get(BIRDCAGE_IMPLANT_CELL_TAG)} (record {REC['tag403_cells']}), 402 "
            f"{n402}/{n402_int}/{n402_pair} (record {REC['n402']}); {W1_LOG}:12004, :12011"
            f"\n[EX-137] (F) box interface: {n_box_owned} owned facets, area {box_area_owned:.12e} "
            f"vs 2(LxLy+LyLz+LxLz) = {surf_ref:.12e}: rel {f_dev:.3e}; unclassified {n_unclass}, "
            f"multi-claimed {n_multi}"
            f"\n[EX-137] (F) face pairs: |x| {face_g[(0, -1)] + face_g[(0, 1)]:.9e} (ref "
            f"{pair_ref[0]:.9e}, rel {pair_dev[0]:.3e}); |y| {face_g[(1, -1)] + face_g[(1, 1)]:.9e} "
            f"(ref {pair_ref[1]:.9e}, rel {pair_dev[1]:.3e}); |z| "
            f"{face_g[(2, -1)] + face_g[(2, 1)]:.9e} (ref {pair_ref[2]:.9e}, rel {pair_dev[2]:.3e})"
            f"\n[EX-137] (R3) max owned tag-3 vertex hypot {r_box_mesh:.12f} vs "
            f"phantom_box_circumradius_m {r_box_cad:.12f}: rel {r3_dev:.3e}"
            f"\n[EX-137] (Cap) max over owned tag-403 vertices of d/r - 1 = {cap_max:.3e} "
            f"(worst vertex {cap_vertex}; bound {STRIP_BAND:.0e}); (cR) min d/r - 1 against "
            f"ROD['centre'] {g36.ROD['centre']} = {cr_min:.6f} (predicted ~88)"
            f"\n[EX-137] (Cap, printed) max |d/r - 1| over owned 402-facet vertices: all "
            f"{wall_all:.3e}, lateral only {wall_lat:.3e}; tag-403 vertex extents lo "
            f"{v403_lo.tolist()} hi {v403_hi.tolist()}; wall clearance Lx/2 - max x(403) = "
            f"{wall_clearance * 1e3:.4f} mm (Lx/2 - (x_rod + r) = "
            f"{(0.5 * Lx - centre[0] - r_rod) * 1e3:.4f} mm)"
            f"\n[EX-137] (cB) refused={cb_refused}, own centre accepted={cb_accepts}, breach "
            f"{breach:.6e} m"
            f"\n[EX-137] (cF) box without its y=+Ly/2 face: area rel miss {cf_dev:.3e} "
            f"(LxLz/total = {Lx * Lz / surf_ref:.3e}) > {STRIP_BAND:.0e}"
            f"\n[EX-137] (N4r) tag-403 cells {census.get(BIRDCAGE_IMPLANT_CELL_TAG, 0)}; meshed "
            f"403 / V {vol403 / V:.6f} (W1 0.986659); meshed 402 / A {g5_ratio:.6f}; conductor "
            f"vs C2 {g36.C2_CONDUCTOR_CAD_MASS_M3:.15e}: rel {n4_cond:.3e} (printed only)"
            f"\n[EX-137] unknowns U2 = {U2_W1} cited from {W1_LOG}:12017 (no N1curl space built here)"
            f"\n[paraview] wrote {path}",
            flush=True,
        )
        _ = (low, high, R, r_leg)

    # Negative controls first.
    assert n1_miss > STRIP_BAND, n1_miss
    assert n2_miss > STRIP_BAND, n2_miss
    assert cf_dev > STRIP_BAND, cf_dev
    assert cr_min > STRIP_BAND, cr_min
    # Records (integers exact; the cell count is printed, not asserted).
    assert census.get(2) == REC["tag2_cells"], census.get(2)
    assert census.get(3) == REC["tag3_cells"], census.get(3)
    assert census.get(BIRDCAGE_IMPLANT_CELL_TAG) == REC["tag403_cells"], census
    assert n402 == REC["n402"] and n402_int == REC["n402"] and n402_pair == REC["n402"], (
        n402, n402_int, n402_pair)
    assert sum(n_gmsh.values()) == REC["n_gmsh_total"], sum(n_gmsh.values())
    # (G1)
    assert g1_mass < STRIP_BAND and g1_kept < STRIP_BAND, (rod_mass, kept, V)
    assert g1_area < STRIP_BAND, (rod["wall_area_m2"], A)
    # (G2)
    assert g2_403 < STRIP_BAND, g2_403
    assert g2_3 < STRIP_BAND, g2_3
    assert g2_sum < STRIP_BAND, g2_sum
    # (G3)
    assert g3_vol <= STRIP_BAND, (vol3, vol403, L3)
    assert g3_ext <= STRIP_BAND, (lo, hi)
    # (G4)
    assert g4 <= STRIP_BAND, (cad_sum, m_cond, m_air)
    # (G5)
    low, high = TERMINAL_AREA_BAND
    assert n402 > 0 and n402_ext == 0, (n402, n402_ext)
    assert n402_int == n402 and n402_pair == n402, (n402, n402_int, n402_pair)
    assert low < g5_ratio < high, g5_ratio
    # (G6)
    assert g6 <= EXACT * l_rod, c
    # (G7)
    assert not census_fails, census_fails
    assert g36.PHANTOM_CELL_TAG in census and BIRDCAGE_IMPLANT_CELL_TAG in census, census
    assert 1 not in census and 404 not in census, census
    assert diag.get("implant_spheroid_cad") is None
    # (G8)
    assert len(n_gmsh) == 2 * LEG_COUNT, sorted(n_gmsh)
    assert not g8_bad, g8_bad
    for i, dev in g8_area.items():
        assert dev <= STRIP_BAND, (i, dev)
    # (F), (R3)
    assert n_unclass == 0 and n_multi == 0, (n_unclass, n_multi)
    assert f_dev <= STRIP_BAND, (box_area_owned, surf_ref)
    for kk, dev in pair_dev.items():
        assert dev <= STRIP_BAND, (kk, dev)
    assert r3_dev <= STRIP_BAND, (r_box_mesh, r_box_cad)
    # (Cap)
    assert cap_max <= STRIP_BAND, (cap_max, cap_vertex)

    if comm.rank == 0:
        print(
            "\nAll assertions hold: the rod-in-box fixture matches the GEO-36 step 8 gate on "
            "this run's mesh; every tag-403 vertex lies inside the capsule at the rod's own "
            "centre and none inside it at step 6's; the mis-placed rod is refused and the "
            "five other negative controls separate. "
            f"Total elapsed {time.perf_counter() - started:.1f} s.",
            flush=True,
        )


if __name__ == "__main__":
    main()
