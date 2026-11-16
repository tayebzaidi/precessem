"""Example (`EX-117`): the F-clinical-1.5T body coil with the ASTM F2182 gel box.

`GEO-38` step 1a (`tests/validation/test_geo38_f_clinical_build.py`) and step 2
(`tests/validation/test_geo38_f_clinical_box.py`) gate a clinical-scale
fixture: 32 legs, 64 ring-gap sheets, a closed PEC can (``shield_radius``),
conductors as holes, and the 0.42 x 0.09 x 0.65 m gel box (``phantom_box``)
as the load. This example builds it once, at the gate's own
``_box_kwargs()`` (never retyped), puts it in ParaView with the 64 sheets and
the box interface as facet tags, and reads two things the gate does not make:
the box **surface** (its six planar faces recovered from the tag-3 / tag-2
interface) and its **reach** in the bore (the farthest box vertex from the
axis against the CAD circumradius the generator exported).

Every anchor is the gate's, recomputed on this run's mesh through the
imported helpers (owned facets / cells only, reduced; no UFL ``dS``):
(T1) owned facets of each sheet = ``ring_sheet_gmsh_triangles``; (T2) per-sheet
area = chord*w; (B1) tag-3 volume = Lx*Ly*Lz; (B2) tag-3 extents; (V1) CAD
partition; (V2) can vs pi R_s^2 2 z_e; (K) census. Example-side: (F) the box
surface closed form and its six faces; (R3) the circumradius identity.

Negative controls, asserted: (c1) one P33 facet dropped; (cF) one box face
missing; (cG) local copies (owned + ghost) are not facets, asserted only if
this partition holds a ghost copy of a box facet.

Scope: a mesh, not a solve. No port gate, S-matrix, tuning, resonance, B1 or
SAR; provisional modelling choices of `GEO-38` stay provisional.

Run it through the example runner::

    ./scripts/run_examples.sh -e mesh:15 -n 4 -t 600
"""

from __future__ import annotations

import math
import resource
import sys
import time
from pathlib import Path

import numpy as np
from mpi4py import MPI

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from fem_em_solver.io.mesh import MeshGenerator, _interface_facet_tags  # noqa: E402
from fem_em_solver.io.paraview_utils import (  # noqa: E402
    adopt_host_ownership,
    write_xdmf_with_tags,
)
from fem_em_solver.post.setup_figure import write_setup_figure  # noqa: E402

from tests.mesh.test_birdcage_conductor_hole import CAD_ANALYTIC_BAND  # noqa: E402
from tests.mesh.test_birdcage_port_sheet_prerequisite import CELL_COUNT_BAND  # noqa: E402
from tests.mesh.test_birdcage_port_sheets import (  # noqa: E402
    PORT_LOWER,
    PORT_UPPER,
    SHEET_IFACE,
)
from tests.validation.test_geo38_f_clinical_box import _box_kwargs  # noqa: E402
from tests.validation.test_geo38_f_clinical_build import (  # noqa: E402
    DESIGN_CORNER_CLEARANCE_M,
    F_CLINICAL_R1_CELL_RECORD,
    F_CLINICAL_R1_UNKNOWNS_RECORD,
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
BASENAME = "meshing_15_f_clinical_coil_astm_box"

# The box interface facet tag: unused by the 64 sheet tags (SHEET_IFACE + i).
BOX_IFACE = 400

# Records, example-local, each with its log line (`GEO-38` step 2 W3, -n 4).
W3_LOG = "20261015T185056Z_GEO-38-s2-W3-n4.log"
REC = {
    "cells": 1_149_150,  # W3_LOG:11898 (PRINTED here, never asserted: 18:00 ruling)
    "tag2_cells": 1_061_834,  # W3_LOG:11898
    "tag3_cells": 33_372,  # W3_LOG:11898
    "n_gmsh_total": 2_882,  # W3_LOG:11899
}


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
    assert len(uniq) == len(idx), "a facet carries two tags (sheet vs box interface)"
    return dolfinx.mesh.meshtags(msh, msh.topology.dim - 1, uniq, val[first])


def main() -> None:
    import dolfinx

    comm = MPI.COMM_WORLD
    started = time.perf_counter()

    kw = _box_kwargs()
    ring_ports = _ring_ports()
    Lx, Ly, Lz = kw["phantom_box"]["size_m"]
    off = kw["phantom_box"]["table_offset_m"]

    # ---- structural asserts at start-up -----------------------------------
    assert comm.size == 4, f"EX-117 is a -n 4 example, comm.size = {comm.size}"
    assert np.issubdtype(dolfinx.default_scalar_type, np.floating), (
        "the meshing group runs the real build"
    )
    assert len(ring_ports) == 2 * LEG_COUNT, len(ring_ports)
    assert off == 0.0, off

    if comm.rank == 0:
        print("=" * 72)
        print("EX-117 - F-clinical-1.5T: 32-leg body coil, closed PEC can, ASTM F2182 box")
        print("=" * 72)
        print(
            f"\n[geometry] legs {LEG_COUNT}  ring ports {len(ring_ports)} (= 2*N_legs)  "
            f"ring_radius {kw['ring_radius']} m  shield_radius {kw['shield_radius']} m"
            f"\n[box]      size_m (Lx, Ly, Lz) = ({Lx}, {Ly}, {Lz}) m  table_offset {off} m"
            "\n[scope]    mesh only: no solve, no port gate, no S, no tuning, no SAR; "
            "modelling choices provisional",
            flush=True,
        )

    # ---- build (timer captured before the figure and XDMF) -----------------
    t0 = time.perf_counter()
    msh, cells, facets, diag = MeshGenerator.birdcage_port_domain(
        comm=comm, return_diagnostics=True, **kw
    )
    t_build = time.perf_counter() - t0
    if comm.rank == 0:
        print(f"\n[EX-117] build done {t_build:.1f} s at -n {comm.size}", flush=True)

    # ---- setup figure (flag-gated; after the build timer) ------------------
    t_fig0 = time.perf_counter()
    write_setup_figure(
        msh,
        cells,
        FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names={
            2: "air",
            3: "phantom (ASTM F2182 gel box)",
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
            "mesh:15 -- F-clinical-1.5T: 32-leg body coil, closed PEC can, "
            "64 ring-gap sheets, ASTM F2182 gel box (GEO-38 C1/C2)"
        ),
        comm=comm,
    )
    t_fig = time.perf_counter() - t_fig0

    # ---- common mesh bookkeeping -------------------------------------------
    tdim = msh.topology.dim
    fdim = tdim - 1
    msh.topology.create_entities(fdim)
    msh.topology.create_connectivity(fdim, tdim)
    n_cells = int(msh.topology.index_map(tdim).size_global)
    census = _owned_census(msh, cells, comm)

    R = kw["ring_radius"]
    r_leg = 0.5 * kw["leg_width"]
    R_s = float(diag["shield_radius_m"])
    z_e = 0.5 * float(diag["shield_length_m"])
    layout = diag["ring_port_layout"]
    chord = float(layout["ring_port_gap_chord_m"])
    width = float(layout["ring_port_box_width_m"])
    sheet_ref = chord * width
    n_gmsh = {int(k): int(v) for k, v in diag["ring_sheet_gmsh_triangles"].items()}

    # ---- (T1), (T2): owned sheet facets -------------------------------------
    sheet_tags = _interface_facet_tags(
        msh, cells, {SHEET_IFACE + i: (PORT_LOWER + i, PORT_UPPER + i) for i in ring_ports}
    )
    sv = np.asarray(sheet_tags.values)
    si = np.asarray(sheet_tags.indices)
    local = {}
    local_copy = {}
    for i in ring_ports:
        f_all = si[sv == SHEET_IFACE + i]
        f = _owned(msh, f_all)
        a = _area(_facet_coords(msh, f))
        local[i] = (int(len(f)), float(a.sum()), float(a[0]) if len(a) else 0.0)
        local_copy[i] = int(len(f_all))
    gathered = comm.allgather(local)
    n_owned = {i: sum(g[i][0] for g in gathered) for i in ring_ports}
    area = {i: sum(g[i][1] for g in gathered) for i in ring_ports}
    n_local_copy = {i: comm.allreduce(local_copy[i], op=MPI.SUM) for i in ring_ports}
    t1_bad = {i: (n_owned[i], n_gmsh.get(i)) for i in ring_ports if n_owned[i] != n_gmsh.get(i)}
    area_dev = {i: abs(area[i] / sheet_ref - 1.0) for i in ring_ports}

    # ---- (c1): one P33 facet dropped, by construction ----------------------
    p33 = ring_ports[0]
    owner = next((r for r, g in enumerate(gathered) if g[p33][0] > 0), -1)
    c1_count = n_owned[p33] - 1
    c1_area = area[p33] - (gathered[owner][p33][2] if owner >= 0 else 0.0)
    c1_dev = abs(c1_area / sheet_ref - 1.0)

    # ---- (B1), (B2), (R3): owned tag-3 tets ---------------------------------
    cv = np.asarray(cells.values)
    ci = np.asarray(cells.indices)
    n_owned_cells = msh.topology.index_map(tdim).size_local
    ph = ci[(cv == 3) & (ci < n_owned_cells)]
    if len(ph):
        g = dolfinx.cpp.mesh.entities_to_geometry(
            msh._cpp_object, tdim, np.asarray(ph, dtype=np.int32), False
        )
        X = msh.geometry.x[g[:, :4]]
        vol_local = float(
            np.abs(
                np.einsum(
                    "ij,ij->i",
                    X[:, 1] - X[:, 0],
                    np.cross(X[:, 2] - X[:, 0], X[:, 3] - X[:, 0]),
                )
            ).sum()
            / 6.0
        )
        P = X.reshape(-1, 3)
        lo_local, hi_local = P.min(axis=0), P.max(axis=0)
        r_local = float(np.hypot(P[:, 0], P[:, 1]).max(initial=-np.inf))
    else:
        vol_local = 0.0
        lo_local, hi_local = np.full(3, np.inf), np.full(3, -np.inf)
        r_local = -np.inf
    box_vol = comm.allreduce(vol_local, op=MPI.SUM)
    lo = np.array([comm.allreduce(float(x), op=MPI.MIN) for x in lo_local])
    hi = np.array([comm.allreduce(float(x), op=MPI.MAX) for x in hi_local])
    r_box_mesh = comm.allreduce(r_local, op=MPI.MAX)
    box_ref = Lx * Ly * Lz
    b1_dev = abs(box_vol / box_ref - 1.0)
    lo_ref = np.array([-0.5 * Lx, -0.5 * Ly + off, -0.5 * Lz])
    hi_ref = np.array([0.5 * Lx, 0.5 * Ly + off, 0.5 * Lz])
    b2_dev = float(
        max(np.max(np.abs(lo - lo_ref) / np.abs(lo_ref)), np.max(np.abs(hi - hi_ref) / np.abs(hi_ref)))
    )
    r_box_cad = float(diag["phantom_box_circumradius_m"])
    r3_dev = abs(r_box_mesh / r_box_cad - 1.0)
    corner_radius = float(diag["leg_box_corner_radius_m"])
    corner_clear = R_s - corner_radius
    rung_clear = R - r_box_cad - r_leg

    # ---- (V1), (V2) ---------------------------------------------------------
    cad_groups = dict(diag["cad_mass_by_group"])
    cad_sum = float(sum(cad_groups.values()))
    m_cond = float(diag["conductor_cad_mass_m3"])
    m_air = float(diag["air_domain_cad_mass_m3"])
    v1_dev = abs((cad_sum + m_cond) / m_air - 1.0)
    can_ref = math.pi * R_s**2 * 2.0 * z_e
    v2_dev = abs(m_air / can_ref - 1.0)

    # ---- (K) ----------------------------------------------------------------
    census_fails = _census_checks(census, n_cells, ring_ports)

    # ---- (F): the box surface ----------------------------------------------
    assert BOX_IFACE not in {SHEET_IFACE + i for i in ring_ports}
    box_tags = _interface_facet_tags(msh, cells, {BOX_IFACE: (3, 2)})
    bi = np.asarray(box_tags.indices)
    bf_all = bi[np.asarray(box_tags.values) == BOX_IFACE]
    bf = _owned(msh, bf_all)
    bxyz = _facet_coords(msh, bf)
    barea = _area(bxyz)
    H = (0.5 * Lx, 0.5 * Ly, 0.5 * Lz)
    C = (0.0, off, 0.0)
    face_area = {}  # (axis, sign) -> local area
    n_claimed = np.zeros(len(bf), dtype=np.int64)
    for k in range(3):
        for s in (-1, 1):
            on = np.all(np.abs((bxyz[:, :, k] - C[k]) - s * H[k]) <= STRIP_BAND * H[k], axis=1)
            n_claimed += on
            face_area[(k, s)] = float(barea[on].sum())
    keys = sorted(face_area)
    loc = np.array([face_area[k] for k in keys] + [float(barea.sum()), float(len(bf)),
                                                   float(np.sum(n_claimed == 0)),
                                                   float(np.sum(n_claimed > 1)),
                                                   float(len(bf_all)),
                                                   float(_area(_facet_coords(msh, bf_all)).sum())])
    tot = np.zeros_like(loc)
    comm.Allreduce(loc, tot, op=MPI.SUM)
    face_g = {k: float(tot[j]) for j, k in enumerate(keys)}
    box_area_owned, n_box_owned, n_unclass, n_multi, n_box_local, box_area_local = (
        float(tot[6]), int(tot[7]), int(tot[8]), int(tot[9]), int(tot[10]), float(tot[11]))
    surf_ref = 2.0 * (Lx * Ly + Ly * Lz + Lx * Lz)
    f_dev = abs(box_area_owned / surf_ref - 1.0)
    pair_ref = {0: 2.0 * Ly * Lz, 1: 2.0 * Lx * Lz, 2: 2.0 * Lx * Ly}
    pair_dev = {k: abs((face_g[(k, -1)] + face_g[(k, 1)]) / pair_ref[k] - 1.0) for k in range(3)}
    # (cF): the five faces other than y = +(Ly/2 + off).
    cf_area = box_area_owned - face_g[(1, 1)]
    cf_dev = abs(cf_area / surf_ref - 1.0)
    # (cG): local copies vs owned.
    ghost_box = n_box_local - n_box_owned
    cg_dev = abs(box_area_local / box_area_owned - 1.0)

    rss = comm.allreduce(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0**2,
                         op=MPI.SUM)

    # ---- ParaView output ----------------------------------------------------
    t_x0 = time.perf_counter()
    merged = _merge_facet_tags(msh, sheet_tags, box_tags)
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
    a0_rel = n_cells / F_CLINICAL_R1_CELL_RECORD - 1.0
    if comm.rank == 0:
        wa = max(area_dev, key=area_dev.get)
        dup = sum(n_local_copy.values()) - sum(n_owned.values())
        t3 = census.get(3)
        print(
            f"\n[EX-117] -n {comm.size}: {n_cells} cells; build {t_build:.1f} s (gmsh "
            f"{diag.get('mesh_wall_time_s', float('nan')):.1f} s); figure {t_fig:.1f} s; XDMF "
            f"{t_xdmf:.1f} s; summed peak RSS {rss:.3f} GiB"
            f"\n[EX-117] census {census}"
            f"\n[EX-117] (A0) cells {n_cells} (PRINTED, not asserted) vs W3 record {REC['cells']} "
            f"({W3_LOG}:11898) and F_CLINICAL_R1_CELL_RECORD {F_CLINICAL_R1_CELL_RECORD}: rel "
            f"{a0_rel:.4e} (CELL_COUNT_BAND {CELL_COUNT_BAND}, not applied)"
            f"\n[EX-117] (T1) owned sheet facets == N_gmsh on {len(ring_ports) - len(t1_bad)}"
            f"/{len(ring_ports)} sheets; mismatches {t1_bad or 'none'}; N_gmsh "
            f"{min(n_gmsh.values())}..{max(n_gmsh.values())}, total {sum(n_gmsh.values())} "
            f"(record {REC['n_gmsh_total']}, {W3_LOG}:11899)"
            f"\n[EX-117] (T2) chord*w = {sheet_ref:.12e}; worst area P{wa} {area_dev[wa]:.3e} "
            f"(band {STRIP_BAND:.0e})"
            f"\n[EX-117] (B1) box meshed volume {box_vol:.15e} vs {box_ref:.15e} m^3: rel {b1_dev:.3e}"
            f"\n[EX-117] (B2) tag-3 extents lo {lo.tolist()} hi {hi.tolist()}: worst rel {b2_dev:.3e}"
            f"\n[EX-117] (V1) CAD partition rel {v1_dev:.3e}; (V2) air domain {m_air:.15e} vs "
            f"pi R_s^2 2z_e {can_ref:.15e}: rel {v2_dev:.3e} (band {CAD_ANALYTIC_BAND:.0e})"
            f"\n[EX-117] (K) census failures {census_fails or 'none'}; tag-3 cells {t3} "
            f"(record {REC['tag3_cells']}), tag-2 cells {census.get(2)} (record {REC['tag2_cells']})"
            f"\n[EX-117] (F) box interface: {n_box_owned} owned facets, area {box_area_owned:.12e} vs "
            f"2(LxLy+LyLz+LxLz) = {surf_ref:.12e}: rel {f_dev:.3e}; unclassified {n_unclass}, "
            f"multi-claimed {n_multi}"
            f"\n[EX-117] (F) face pairs: |x| {face_g[(0, -1)] + face_g[(0, 1)]:.9e} (ref "
            f"{pair_ref[0]:.9e}, rel {pair_dev[0]:.3e}); |y| {face_g[(1, -1)] + face_g[(1, 1)]:.9e} "
            f"(ref {pair_ref[1]:.9e}, rel {pair_dev[1]:.3e}); |z| "
            f"{face_g[(2, -1)] + face_g[(2, 1)]:.9e} (ref {pair_ref[2]:.9e}, rel {pair_dev[2]:.3e})"
            f"\n[EX-117] (R3) max owned tag-3 vertex hypot {r_box_mesh:.12f} vs "
            f"phantom_box_circumradius_m {r_box_cad:.12f}: rel {r3_dev:.3e}"
            f"\n[EX-117] clearance: R - r_box - r_leg = {rung_clear * 1e3:.3f} mm "
            f"(R {R}, r_leg {r_leg}); shield R_s - leg_box_corner_radius = {corner_clear * 1e3:.4f} mm "
            f"(design >= {DESIGN_CORNER_CLEARANCE_M * 1e3:.1f} mm)"
            f"\n[EX-117] volume fraction tag-3 / can = {box_vol / m_air:.4f} (B1 / V2 air domain)"
            f"\n[EX-117] unknowns record (printed, no N1curl space built here): "
            f"F_CLINICAL_R1_UNKNOWNS_RECORD {F_CLINICAL_R1_UNKNOWNS_RECORD}"
            f"\n[EX-117] (c1) P{p33} one facet dropped on rank {owner}: count {c1_count} vs "
            f"N_gmsh-1 {n_gmsh[p33] - 1}; area rel {c1_dev:.3e} > {STRIP_BAND:.0e}"
            f"\n[EX-117] (cF) box without its y=+Ly/2 face: area rel miss {cf_dev:.3e} "
            f"(LxLz/total = {Lx * Lz / surf_ref:.3e}) > {STRIP_BAND:.0e}"
            f"\n[EX-117] (cG) box facet local copies {n_box_local} vs owned {n_box_owned} "
            f"(ghost copies {ghost_box}); area local-copy sum {box_area_local:.9e} vs owned "
            f"{box_area_owned:.9e}: rel {cg_dev:.3e}"
            + ("" if ghost_box > 0 else "\n[EX-117] (cG) no ghost copy on this partition")
            + f"\n[EX-117] (c2 sheets, printed) {sum(n_local_copy.values())} local copies vs "
            f"{sum(n_owned.values())} owned (excess {dup})"
            f"\n[paraview] wrote {path}",
            flush=True,
        )

    # Negative controls first.
    assert owner >= 0
    assert c1_count == n_gmsh[p33] - 1, (c1_count, n_gmsh[p33])
    assert c1_dev > STRIP_BAND, c1_dev
    assert cf_dev > STRIP_BAND, cf_dev
    if ghost_box > 0:
        assert cg_dev > STRIP_BAND, (box_area_local, box_area_owned)
    # Records (integers exact; the cell count is printed, not asserted).
    assert census.get(2) == REC["tag2_cells"], census.get(2)
    assert census.get(3) == REC["tag3_cells"], census.get(3)
    assert sum(n_gmsh.values()) == REC["n_gmsh_total"], sum(n_gmsh.values())
    # (T1), (T2)
    assert len(n_gmsh) == 2 * LEG_COUNT, sorted(n_gmsh)
    assert not t1_bad, t1_bad
    for i, dev in area_dev.items():
        assert dev <= STRIP_BAND, (i, dev)
    # (B1), (B2)
    assert b1_dev <= STRIP_BAND, (box_vol, box_ref)
    assert b2_dev <= STRIP_BAND, (lo, hi)
    # (V1), (V2)
    assert v1_dev <= STRIP_BAND, (cad_sum, m_cond, m_air)
    assert v2_dev <= CAD_ANALYTIC_BAND, (m_air, can_ref)
    # (K)
    assert 3 in census and 1 not in census, census
    assert not census_fails, census_fails
    # (F)
    assert n_unclass == 0 and n_multi == 0, (n_unclass, n_multi)
    assert f_dev <= STRIP_BAND, (box_area_owned, surf_ref)
    for k, dev in pair_dev.items():
        assert dev <= STRIP_BAND, (k, dev)
    # (R3)
    assert r3_dev <= STRIP_BAND, (r_box_mesh, r_box_cad)

    if comm.rank == 0:
        print(
            "\nAll assertions hold: the F-clinical fixture's sheets, box and can match the "
            "GEO-38 gates on this run's mesh; the box surface and reach read exactly; the "
            f"three negative controls separate. Total elapsed {time.perf_counter() - started:.1f} s.",
            flush=True,
        )


if __name__ == "__main__":
    main()
