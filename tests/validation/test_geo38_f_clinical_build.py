"""`GEO-38` step 1a (§9 item 468, chain C1 (a)) — F-clinical-1.5T built:
CAD identities at ``-n 2``, cells and unknowns at two sizing rungs.

Selected by environment only (``FEM_EM_GEO38_STEP1A`` in {``r1``, ``r2``}),
never by ``-k``; skipped when unset. Brief: docs/planning/chunks/GEO-38.md,
"§9 item 468" and "Review ruling 2026-10-15 03:00".

Anchors (asserted, owned facets / vertices, allreduced):
  (R) every vertex of the 64 longitudinal ring sheets lies on u = R
      (u = x cos φ_c + y sin φ_c) within ``STRIP_BAND`` 1e-9;
  (L) per sheet (z_max + z_min)/2 = ±leg_spacing/2 (rung length 0.45 m) and
      the cavity wall's max z = coil_length/2, at 1e-9;
  (S) outer vertices off the caps sit on hypot(x, y) = R_s, max|z| over the
      outer facets = z_extent, at 1e-9; no cavity-wall vertex reaches R_s and
      the outer tag (1) is exactly exterior ∖ cavity (the trap check);
  (N) 32 legs from the cavity-wall end caps at z = +coil_length/2;
  (A) 64 sheet areas = chord·w within ``SHEET_AREA_BAND`` 1e-9;
  (C) ring-port census N_ports = 2·N_legs (64; the review ruling 2026-10-15
      03:00 returns the §7 "census = 32" to the weekly, printed beside it);
  (U) U₁ = E_global, U₂ = 2·E_global + 2·F_global.
Negative control (c1), asserted by construction: one P33 facet dropped moves
P33's area off chord·w by more than ``SHEET_AREA_BAND``.
Scope: geometry and price only — no solve, no port gate, no physics.
"""

from __future__ import annotations

import math
import os
import resource
import time

import numpy as np
import pytest
from mpi4py import MPI

STEP1A_ENV = "FEM_EM_GEO38_STEP1A"
LEG_COUNT = 32
IDENTITY_BAND = 1.0e-9
WEEKLY_CENSUS_LETTER = 32  # printed beside the measured census (rule (e)), never asserted

# The parameter set (brief "Shared: the F-clinical-1.5T parameter set"); the
# single source of truth for items 468–470.
F_CLINICAL_15T_KWARGS = dict(
    leg_count=LEG_COUNT,
    ring_radius=0.355,
    leg_spacing=0.45,
    coil_length=0.45 + 2.0 * 0.0321428571,
    leg_width=0.012,
    ring_minor_radius=0.004,
    ring_gap_length=0.0171428571,
    # dx 21.43 -> 4 mm (review ruling 2026-10-15 10:30, item 478): the
    # axis-aligned leg boxes' outer corners must clear the can; at dx = 4 mm
    # the corner max is 0.3668149 m (5.685 mm clearance, boxes at 39.375°).
    port_box_size=(0.004, 0.004, 0.0214285714),
    phantom_radius=0.10969,
    phantom_height=0.65,
    air_padding=(0.03 / 0.07) * 0.355,
    shield_radius=0.3725,
    as_hole=True,
    emit_port_sheets=True,
    ring_sheet_orientation="longitudinal",
    resolution=0.015,
)
# Attempt 1's table (item 468): dx = 21.43 mm put the leg-box corners at
# 0.37486 m > R_s (`20261015T114056Z_GEO-38-s1a-W2-r2.log`); the corner
# refusal must reject it before gmsh (negative control for the refusal).
ATTEMPT1_PORT_BOX_SIZE = (0.0214285714, 0.004, 0.0214285714)
# Provisional modelling choice (review 2026-10-15 03:00 / 10:30), not a band.
DESIGN_CORNER_CLEARANCE_M = 5.0e-3
CONDUCTOR_RESOLUTION = {"r1": 1.6e-3, "r2": 2.4e-3}

# Printed predictions (brief, "Printed predictions, written first").
PREDICTED_CELLS = {"r1": (1.0e6, (0.8e6, 1.5e6)), "r2": (0.76e6, (0.62e6, 1.15e6))}
PREDICTED_U1 = {"r1": 1.34e6, "r2": 1.03e6}
PREDICTED_U2 = {"r1": 6.9e6, "r2": 5.3e6}
PREDICTED_BUILD_S = {"r1": (340.0, (280.0, 500.0)), "r2": (260.0, (210.0, 380.0))}
F_HUMAN_N1 = 369_387

# Version-tagged cell records ((1*) licence), registered from this step's
# measured values; items 469 / 470 assert them at CELL_COUNT_BAND. Measured
# 2026-10-15 on the default image at -n 2, dx = 4 mm (item 478):
# r1 `20261015T170721Z_GEO-38-s1a-W1-r1-a2.log`, r2
# `20261015T171259Z_GEO-38-s1a-W2-r2-a2.log`. Read only by items 469 / 470;
# this module's windowed path never reads them (set after its last window).
F_CLINICAL_R1_CELL_RECORD = 1_147_366
F_CLINICAL_R2_CELL_RECORD = 912_804
# Unknown-count records at the same windows: (U1 = E, U2 = 2E + 2F).
F_CLINICAL_R1_UNKNOWNS_RECORD = (1_525_193, 7_885_674)
F_CLINICAL_R2_UNKNOWNS_RECORD = (1_184_518, 6_179_214)


def _rung():
    return os.environ.get(STEP1A_ENV, "").strip()


def _f_clinical_kwargs(rung):
    kw = dict(F_CLINICAL_15T_KWARGS)
    kw["conductor_resolution"] = CONDUCTOR_RESOLUTION[rung]
    return kw


def _ring_ports():
    return list(range(LEG_COUNT + 1, 3 * LEG_COUNT + 1))


def _facet_coords(msh, facets):
    """(n, 3, 3) vertex coordinates of the given local facets."""
    import dolfinx

    fdim = msh.topology.dim - 1
    if len(facets) == 0:
        return np.zeros((0, 3, 3))
    geom = dolfinx.cpp.mesh.entities_to_geometry(
        msh._cpp_object, fdim, np.asarray(facets, dtype=np.int32), False)
    return msh.geometry.x[geom[:, :3]]


def _owned(msh, facets):
    n_owned = msh.topology.index_map(msh.topology.dim - 1).size_local
    facets = np.asarray(facets)
    return facets[facets < n_owned]


def _area(xyz):
    if len(xyz) == 0:
        return np.zeros(0)
    return 0.5 * np.linalg.norm(np.cross(xyz[:, 1] - xyz[:, 0], xyz[:, 2] - xyz[:, 0]), axis=1)


def test_geo38_step1a_f_clinical_identities():
    rung = _rung()
    if rung not in CONDUCTOR_RESOLUTION:
        pytest.skip(f"{STEP1A_ENV} unset: GEO-38 step 1a is selected by environment only")
    import dolfinx
    from dolfinx import fem

    from fem_em_solver.io.mesh import (
        BIRDCAGE_CONDUCTOR_SURFACE_TAG,
        MeshGenerator,
        _interface_facet_tags,
    )
    from tests.mesh.test_birdcage_leg_offset import SHEET_AREA_BAND
    from tests.mesh.test_birdcage_port_sheets import PORT_LOWER, PORT_UPPER, SHEET_IFACE
    from tests.validation.test_ops61_exact_strip_ports import STRIP_BAND
    from tests.validation.test_wf12_f_human_hole import _census_checks, _owned_census

    comm = MPI.COMM_WORLD
    kw = _f_clinical_kwargs(rung)
    ring_ports = _ring_ports()

    t0 = time.perf_counter()
    msh, cells, facets, diag = MeshGenerator.birdcage_port_domain(
        comm=comm, return_diagnostics=True, **kw)
    t_build = time.perf_counter() - t0
    tdim = msh.topology.dim
    fdim = tdim - 1
    msh.topology.create_entities(1)
    msh.topology.create_entities(fdim)
    msh.topology.create_connectivity(fdim, tdim)
    n_cells = int(msh.topology.index_map(tdim).size_global)
    census = _owned_census(msh, cells, comm)

    R = kw["ring_radius"]
    r_leg = 0.5 * kw["leg_width"]
    R_s = float(diag["shield_radius_m"])
    z_e = 0.5 * float(diag["shield_length_m"])
    layout = diag["ring_port_layout"]
    sheet_analytic = float(layout["ring_port_sheet_longitudinal_area_m2"])
    chord = float(layout["ring_port_gap_chord_m"])
    width = float(layout["ring_port_box_width_m"])

    # --- sheets: (R), (L), (A)
    sheet_tags = _interface_facet_tags(
        msh, cells, {SHEET_IFACE + i: (PORT_LOWER + i, PORT_UPPER + i) for i in ring_ports})
    sv = np.asarray(sheet_tags.values)
    si = np.asarray(sheet_tags.indices)
    local_sheets = {}
    for i in ring_ports:
        f = _owned(msh, si[sv == SHEET_IFACE + i])
        local_sheets[i] = (f, _facet_coords(msh, f))
    gathered = comm.allgather({i: v[1] for i, v in local_sheets.items()})
    sheet = {i: np.concatenate([g[i] for g in gathered]) for i in ring_ports}
    sheet_area = {i: float(_area(sheet[i]).sum()) for i in ring_ports}
    sheet_nf = {i: len(sheet[i]) for i in ring_ports}
    sheet_dev = {i: abs(a / sheet_analytic - 1.0) for i, a in sheet_area.items()}
    r_dev, z_mid = {}, {}
    for i, xyz in sheet.items():
        p = xyz.reshape(-1, 3)
        if len(p) == 0:
            r_dev[i], z_mid[i] = float("inf"), float("nan")
            continue
        # φ_c is the gap's mid-azimuth (k + 1/2)·2π/N, not the vertex mean (W1
        # r1: the vertex-mean estimate read 1.96e-5 from an uneven node spread).
        dphi = 2.0 * math.pi / LEG_COUNT
        phi = (math.floor(math.atan2(p[:, 1].mean(), p[:, 0].mean()) / dphi) + 0.5) * dphi
        u = p[:, 0] * math.cos(phi) + p[:, 1] * math.sin(phi)
        r_dev[i] = float(np.max(np.abs(u / R - 1.0)))
        z_mid[i] = 0.5 * (p[:, 2].max() + p[:, 2].min())
    half_sp = 0.5 * kw["leg_spacing"]
    l_dev = {i: abs(abs(z) / half_sp - 1.0) for i, z in z_mid.items()}
    z_signs = sorted(int(np.sign(z)) for z in z_mid.values())
    separation = (max(z_mid.values()) - min(z_mid.values()))

    # --- cavity wall and outer facets: (L) leg span, (S), (N)
    fv = np.asarray(facets.values)
    fi = np.asarray(facets.indices)
    cav = _owned(msh, fi[fv == BIRDCAGE_CONDUCTOR_SURFACE_TAG])
    cav_xyz = _facet_coords(msh, cav)
    ext = _owned(msh, dolfinx.mesh.exterior_facet_indices(msh.topology))
    outer = np.setdiff1d(ext, cav)
    tag1 = _owned(msh, fi[fv == 1])
    n_outer_mismatch = comm.allreduce(
        int(len(np.setxor1d(outer, tag1))), op=MPI.SUM)
    n_tag_collision = comm.allreduce(int(len(np.intersect1d(tag1, cav))), op=MPI.SUM)
    n_ext = comm.allreduce(int(len(ext)), op=MPI.SUM)
    n_outer = comm.allreduce(int(len(outer)), op=MPI.SUM)
    n_cav = comm.allreduce(int(len(cav)), op=MPI.SUM)
    out_p = _facet_coords(msh, outer).reshape(-1, 3)
    cav_p = cav_xyz.reshape(-1, 3)
    cav_zmax = comm.allreduce(float(cav_p[:, 2].max()) if len(cav_p) else -np.inf, op=MPI.MAX)
    cav_hmax = comm.allreduce(
        float(np.hypot(cav_p[:, 0], cav_p[:, 1]).max()) if len(cav_p) else 0.0, op=MPI.MAX)
    if len(cav_p):
        k = int(np.argmax(np.hypot(cav_p[:, 0], cav_p[:, 1])))
        print(f"[GEO-38 s1a] rank {comm.rank} cavity max-hypot vertex {cav_p[k]}", flush=True)
    lat = out_p[np.abs(out_p[:, 2]) < z_e * (1.0 - 1.0e-9)]
    shield_dev = comm.allreduce(
        float(np.max(np.abs(np.hypot(lat[:, 0], lat[:, 1]) / R_s - 1.0))) if len(lat) else 0.0,
        op=MPI.MAX)
    n_lat = comm.allreduce(int(len(lat)), op=MPI.SUM)
    out_zmax = comm.allreduce(float(np.abs(out_p[:, 2]).max()) if len(out_p) else 0.0, op=MPI.MAX)
    half_cl = 0.5 * kw["coil_length"]
    leg_dev = abs(cav_zmax / half_cl - 1.0)
    # (N): cavity facets with all three vertices on z = +coil_length/2.
    if len(cav_xyz):
        on_cap = np.all(np.abs(cav_xyz[:, :, 2] / half_cl - 1.0) <= IDENTITY_BAND, axis=1)
        cap_p = cav_xyz[on_cap].reshape(-1, 3)
    else:
        cap_p = np.zeros((0, 3))
    dth = 2.0 * math.pi / LEG_COUNT
    bins = np.mod(np.rint(np.arctan2(cap_p[:, 1], cap_p[:, 0]) / dth).astype(int), LEG_COUNT)
    th = bins * dth
    axis_d = np.hypot(cap_p[:, 0] - R * np.cos(th), cap_p[:, 1] - R * np.sin(th))
    leg_bins = set()
    for s in comm.allgather(set(int(b) for b in bins)):
        leg_bins |= s
    axis_rel = comm.allreduce(float((axis_d / r_leg).max()) if len(axis_d) else 0.0, op=MPI.MAX)

    # --- (S') leg-box corner identity: max vertex hypot over cells 101..132
    # equals the refusal's own `leg_box_corner_radius_m` (planar boxes, so the
    # corner is a mesh vertex). Local cells incl. ghosts; MAX is idempotent.
    cv = np.asarray(cells.values)
    ci = np.asarray(cells.indices)
    box_cells = ci[(cv >= PORT_LOWER + 1) & (cv <= PORT_LOWER + LEG_COUNT)]
    if len(box_cells):
        bgeom = dolfinx.cpp.mesh.entities_to_geometry(
            msh._cpp_object, tdim, np.asarray(box_cells, dtype=np.int32), False)
        bp = msh.geometry.x[bgeom.reshape(-1)]
        box_hmax_local = float(np.hypot(bp[:, 0], bp[:, 1]).max())
    else:
        box_hmax_local = 0.0
    box_hmax = comm.allreduce(box_hmax_local, op=MPI.MAX)
    corner_ref = float(diag["leg_box_corner_radius_m"])
    corner_dev = abs(box_hmax - corner_ref)
    corner_clear = R_s - corner_ref

    # --- (C) census
    census_fails = _census_checks(census, n_cells, ring_ports)
    leg_boxes = [PORT_LOWER + k for k in range(1, LEG_COUNT + 1)]
    missing_leg_boxes = [t for t in leg_boxes if t not in census]
    empty_sheets = [i for i in ring_ports if sheet_nf[i] == 0]
    n_ports = len([i for i in ring_ports if sheet_nf[i] > 0])

    # --- (U) unknowns
    E = int(msh.topology.index_map(1).size_global)
    F = int(msh.topology.index_map(2).size_global)
    V1 = fem.functionspace(msh, ("N1curl", 1))
    U1 = int(V1.dofmap.index_map.size_global * V1.dofmap.index_map_bs)
    del V1
    V2 = fem.functionspace(msh, ("N1curl", 2))
    U2 = int(V2.dofmap.index_map.size_global * V2.dofmap.index_map_bs)
    del V2

    # --- (c1) planted dropped facet on P33, by construction.
    p33 = ring_ports[0]
    have = comm.allgather(len(local_sheets[p33][0]) > 0)
    owner = have.index(True) if any(have) else -1
    xyz = local_sheets[p33][1]
    local_a = float(_area(xyz[1:] if comm.rank == owner else xyz).sum())
    area_drop = comm.allreduce(local_a, op=MPI.SUM)
    c1_dev = abs(area_drop / (chord * width) - 1.0)

    rss = comm.allreduce(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0**2,
                         op=MPI.SUM)

    if comm.rank == 0:
        rho = U1 / F_HUMAN_N1
        pc, prange = PREDICTED_CELLS[rung]
        pb, brange = PREDICTED_BUILD_S[rung]
        worst_a = max(sheet_dev, key=sheet_dev.get)
        worst_r = max(r_dev, key=r_dev.get)
        worst_l = max(l_dev, key=l_dev.get)
        print(
            f"\n[GEO-38 s1a] rung {rung} (conductor_resolution {kw['conductor_resolution']}), "
            f"-n {comm.size}: {n_cells} cells (PREDICTED ≈{pc:.3g}, range {prange}); build "
            f"{t_build:.1f} s (PREDICTED ≈{pb} s, range {brange}; gmsh "
            f"{diag.get('mesh_wall_time_s', float('nan')):.1f} s); summed peak RSS {rss:.3f} GiB\n"
            f"[GEO-38 s1a] census {census}\n"
            f"[GEO-38 s1a] (R) rung diameter 2R = {2 * R:.6f} m: worst P{worst_r} max|u/R-1| "
            f"{r_dev[worst_r]:.3e} (band {STRIP_BAND:.0e})\n"
            f"[GEO-38 s1a] (L) sheet z-mid worst P{worst_l} ||z|/(s/2)-1| {l_dev[worst_l]:.3e}; "
            f"signs {z_signs.count(-1)} lower / {z_signs.count(1)} upper; separation "
            f"{separation:.12f} m (leg_spacing {kw['leg_spacing']}); cavity max z {cav_zmax:.12f} "
            f"vs coil_length/2 {half_cl:.12f}: rel {leg_dev:.3e}; leg span {2 * cav_zmax:.9f} m\n"
            f"[GEO-38 s1a] (S) shield R_s {R_s} (Ø {2 * R_s:.6f} m), length {2 * z_e:.9f} m, CAD air "
            f"mass {diag['air_domain_cad_mass_m3']:.9e} m^3 (π R_s² 2z_e = "
            f"{math.pi * R_s ** 2 * 2 * z_e:.9e}); lateral vertices {n_lat}: max|hypot/R_s-1| "
            f"{shield_dev:.3e}; outer max|z| {out_zmax:.12f} vs z_e {z_e:.12f} rel "
            f"{abs(out_zmax / z_e - 1.0):.3e}; cavity max hypot {cav_hmax:.9f} (< R_s(1-1e-9) "
            f"{R_s * (1 - 1e-9):.9f}); exterior {n_ext} = outer {n_outer} + cavity {n_cav}; "
            f"tag-1 vs outer mismatches {n_outer_mismatch}, tag-1∩cavity {n_tag_collision}\n"
            f"[GEO-38 s1a] (S') leg-box max vertex hypot {box_hmax:.12f} vs leg_box_corner_radius_m "
            f"{corner_ref:.12f}: |diff| {corner_dev:.3e} m (band {IDENTITY_BAND:.0e}); corner "
            f"clearance R_s - corner {corner_clear * 1e3:.4f} mm (design >= "
            f"{DESIGN_CORNER_CLEARANCE_M * 1e3:.1f} mm); ring-box corner "
            f"{float(diag['ring_box_corner_radius_m']):.9f} m\n"
            f"[GEO-38 s1a] (N) leg bins {len(leg_bins)} (of {LEG_COUNT}); max axis distance / "
            f"r_leg {axis_rel:.12f}\n"
            f"[GEO-38 s1a] (A) 64 sheets vs chord*w = {chord:.9e}*{width:.9e} = {chord * width:.9e} "
            f"(layout {sheet_analytic:.9e}); worst P{worst_a} {sheet_area[worst_a]:.9e} rel "
            f"{sheet_dev[worst_a]:.3e} (band {SHEET_AREA_BAND:.0e}); facets per sheet "
            f"{min(sheet_nf.values())}..{max(sheet_nf.values())}\n"
            f"[GEO-38 s1a] (C) N_ports {n_ports} = 2*N_legs {2 * LEG_COUNT} (weekly's letter "
            f"{WEEKLY_CENSUS_LETTER}, PREDICTED != measured, printed only); census failures "
            f"{census_fails or 'none'}; missing leg boxes {missing_leg_boxes or 'none'}\n"
            f"[GEO-38 s1a] (U) E {E}, F {F}; U1 {U1} (= E: {U1 == E}; PREDICTED ≈{PREDICTED_U1[rung]:.3g}); "
            f"U2 {U2} = 2E+2F {2 * E + 2 * F}: {U2 == 2 * E + 2 * F} (PREDICTED ≈{PREDICTED_U2[rung]:.3g})\n"
            f"[GEO-38 s1a] (c1) P{p33} with one facet dropped on rank {owner}: area {area_drop:.9e}, "
            f"rel {c1_dev:.3e} > band {SHEET_AREA_BAND:.0e} (ceiling ≈ 1/{sheet_nf[p33]} = "
            f"{1.0 / max(sheet_nf[p33], 1):.3e})\n"
            f"[GEO-38 s1a] PRICE (prediction for the weekly) degree-2 xl at -n 16: rho = U1/369387 = "
            f"{rho:.4f}; t ≈ t_b + {287.1 * rho ** 0.84:.1f} s; memory {89.31 * rho ** 0.84:.1f} GiB "
            f"(x rho^(4/3): {89.31 * rho ** (4 / 3):.1f} GiB)",
            flush=True,
        )

    # Negative control first (by construction).
    assert owner >= 0
    assert c1_dev > SHEET_AREA_BAND, (area_drop, chord * width, c1_dev)
    # (A) and (C)
    assert not empty_sheets, empty_sheets
    assert n_ports == 2 * LEG_COUNT, n_ports
    for i, dev in sheet_dev.items():
        assert dev < SHEET_AREA_BAND, (i, sheet_area[i], sheet_analytic, dev)
    assert not census_fails, census_fails
    assert not missing_leg_boxes, missing_leg_boxes
    # (R)
    for i, dev in r_dev.items():
        assert dev <= STRIP_BAND, (i, dev)
    # (L)
    for i, dev in l_dev.items():
        assert dev <= IDENTITY_BAND, (i, z_mid[i], dev)
    assert z_signs.count(-1) == LEG_COUNT and z_signs.count(1) == LEG_COUNT, z_signs
    assert abs(separation / kw["leg_spacing"] - 1.0) <= IDENTITY_BAND, separation
    assert leg_dev <= IDENTITY_BAND, (cav_zmax, half_cl)
    # (S)
    assert n_lat > 0
    assert shield_dev <= IDENTITY_BAND, shield_dev
    assert abs(out_zmax / z_e - 1.0) <= IDENTITY_BAND, (out_zmax, z_e)
    assert cav_hmax < R_s * (1.0 - 1.0e-9), cav_hmax
    assert n_outer_mismatch == 0 and n_tag_collision == 0, (n_outer_mismatch, n_tag_collision)
    # (S')
    assert corner_dev <= IDENTITY_BAND, (box_hmax, corner_ref, corner_dev)
    assert corner_clear >= DESIGN_CORNER_CLEARANCE_M, corner_clear
    # (N)
    assert len(leg_bins) == LEG_COUNT, sorted(leg_bins)
    assert axis_rel <= 1.0 + IDENTITY_BAND, axis_rel
    # (U)
    assert U1 == E, (U1, E)
    assert U2 == 2 * E + 2 * F, (U2, E, F)


def _refuse_before_gmsh(monkeypatch, **kw):
    """Call `birdcage_port_domain` with the builder replaced by a tripwire, so a
    passing refusal never meshes and a missing one fails fast (no gmsh)."""
    from fem_em_solver.io.mesh import MeshGenerator

    def _tripwire(*a, **k):
        raise AssertionError("reached the gmsh builder: the shield_radius refusal did not fire")

    monkeypatch.setattr(MeshGenerator, "_build_birdcage_port_model", staticmethod(_tripwire))
    return MeshGenerator.birdcage_port_domain(comm=MPI.COMM_WORLD, return_diagnostics=True, **kw)


def test_geo38_shield_refuses_attempt1_table(monkeypatch):
    """Negative control for the corner refusal: attempt 1's dx = 21.43 mm table
    passed the old face-centre check and leaked the leg boxes past R_s (item
    468 W2); the corner formula must refuse it before gmsh."""
    kw = _f_clinical_kwargs("r2")
    kw["port_box_size"] = ATTEMPT1_PORT_BOX_SIZE
    with pytest.raises(ValueError, match="leg-box corner") as exc:
        _refuse_before_gmsh(monkeypatch, **kw)
    if MPI.COMM_WORLD.rank == 0:
        print(f"\n[GEO-38 s1a] refusal (attempt-1 table, dx 21.43 mm): {exc.value}", flush=True)


def test_geo38_shield_refuses_leg_gap_length(monkeypatch):
    """`shield_radius` is high-pass only (review 2026-10-15 10:30, item 3)."""
    kw = _f_clinical_kwargs("r2")
    # Leg gaps with leg sheets, solid coil: `as_hole` and the ring-gap keywords
    # are dropped so no earlier refusal fires first (smoke run
    # `20261015T170229Z_GEO-38-s1a-a2-refusal.log` hit `as_hole needs ...`).
    for key in ("ring_gap_length", "ring_sheet_orientation", "as_hole"):
        kw.pop(key)
    kw["leg_gap_length"] = 0.01
    with pytest.raises(ValueError, match="shield_radius cannot be combined with leg_gap_length"):
        _refuse_before_gmsh(monkeypatch, **kw)
