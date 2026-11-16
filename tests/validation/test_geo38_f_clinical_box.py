"""`GEO-38` step 2 (§9 item 470, chain C2) — the ASTM F2182 box as a tagged
load in the F-clinical-1.5T coil: ring-sheet tags complete and unduplicated
at ``-n 1 / 2 / 4``, congruence, box volume and placement, CAD partition.

Selected by environment only (``FEM_EM_GEO38_STEP2=1``), never by ``-k``;
skipped when unset. Brief: docs/planning/chunks/GEO-38.md, "§9 item 470".

Anchors (asserted at the window's width; owned facets / cells, reduced):
  (T1) per ring sheet, owned facets of tag 210+i (reduced) = N_gmsh[i]
       (``ring_sheet_gmsh_triangles``, `GEO-37`'s (A2));
  (T2) per sheet: area = chord·w, vertices on u = R, z-extent = w, chord
       extent = chord, all at ``STRIP_BAND``;
  (B1) Σ owned tag-3 tet volumes = 0.65 × 0.42 × 0.09 m³ at ``STRIP_BAND``;
  (B2) tag-3 vertex extents ±0.21 (x), ±0.045 + offset (y), ±0.325 (z);
  (V1) Σ cad_mass_by_group + conductor_cad_mass_m3 = air_domain_cad_mass_m3
       at ``STRIP_BAND``;
  (V2) air_domain_cad_mass_m3 vs π R_s² · 2z_e at ``CAD_ANALYTIC_BAND``;
  (K)  tag 3 present, no tag 1, the census sums to the global cell count.
(A0) cells vs ``F_CLINICAL_R1_CELL_RECORD`` printed (the record is asserted
only if the move is under ``CELL_COUNT_BAND``; else a box record is
registered from the measurement).
Negative control (c1), asserted by construction: one P33 facet dropped
reads N_gmsh − 1 and misses chord·w by more than ``STRIP_BAND``; (c2) the
local-copy double count (owned + ghost facets) printed only.
No UFL ``dS`` over the sheets: owned-facet numpy sums only.
Scope: load geometry and mesh gates only — no solve, no SAR.
"""

from __future__ import annotations

import math
import os
import resource
import time

import numpy as np
import pytest
from mpi4py import MPI

from tests.validation.test_geo38_f_clinical_build import (
    F_CLINICAL_R1_CELL_RECORD,
    LEG_COUNT,
    _area,
    _f_clinical_kwargs,
    _facet_coords,
    _owned,
    _ring_ports,
)

STEP2_ENV = "FEM_EM_GEO38_STEP2"
TABLE_OFFSET_M = 0.0  # centred; a nonzero value only with a cited source


def _box_kwargs():
    from tests.validation.test_mat21_box_closed_form import A_BOX, B_BOX, C_BOX

    kw = _f_clinical_kwargs("r1")
    kw.pop("phantom_radius")
    kw.pop("phantom_height")
    # Lz = 2A on the coil axis, Lx = 2B horizontal, Ly = 2C vertical gel depth.
    kw["phantom_box"] = {
        "size_m": (2.0 * B_BOX, 2.0 * C_BOX, 2.0 * A_BOX),
        "table_offset_m": TABLE_OFFSET_M,
    }
    return kw


def test_geo38_step2_astm_box_mesh_gates():
    if os.environ.get(STEP2_ENV, "").strip() != "1":
        pytest.skip(f"{STEP2_ENV} unset: GEO-38 step 2 is selected by environment only")
    import dolfinx

    from fem_em_solver.io.mesh import MeshGenerator, _interface_facet_tags
    from tests.mesh.test_birdcage_conductor_hole import CAD_ANALYTIC_BAND
    from tests.mesh.test_birdcage_port_sheet_prerequisite import CELL_COUNT_BAND
    from tests.mesh.test_birdcage_port_sheets import PORT_LOWER, PORT_UPPER, SHEET_IFACE
    from tests.validation.test_ops61_exact_strip_ports import STRIP_BAND
    from tests.validation.test_wf12_f_human_hole import _census_checks, _owned_census

    comm = MPI.COMM_WORLD
    kw = _box_kwargs()
    ring_ports = _ring_ports()
    Lx, Ly, Lz = kw["phantom_box"]["size_m"]
    off = kw["phantom_box"]["table_offset_m"]

    t0 = time.perf_counter()
    msh, cells, facets, diag = MeshGenerator.birdcage_port_domain(
        comm=comm, return_diagnostics=True, **kw)
    t_build = time.perf_counter() - t0
    if comm.rank == 0:
        print(f"\n[GEO-38 s2] build done {t_build:.1f} s at -n {comm.size}", flush=True)
    tdim = msh.topology.dim
    fdim = tdim - 1
    msh.topology.create_entities(fdim)
    msh.topology.create_connectivity(fdim, tdim)
    n_cells = int(msh.topology.index_map(tdim).size_global)
    census = _owned_census(msh, cells, comm)

    R = kw["ring_radius"]
    R_s = float(diag["shield_radius_m"])
    z_e = 0.5 * float(diag["shield_length_m"])
    layout = diag["ring_port_layout"]
    chord = float(layout["ring_port_gap_chord_m"])
    width = float(layout["ring_port_box_width_m"])
    sheet_ref = chord * width
    n_gmsh = {int(k): int(v) for k, v in diag["ring_sheet_gmsh_triangles"].items()}

    # --- (T1), (T2): owned sheet facets, numpy sums, gathered once.
    sheet_tags = _interface_facet_tags(
        msh, cells, {SHEET_IFACE + i: (PORT_LOWER + i, PORT_UPPER + i) for i in ring_ports})
    sv = np.asarray(sheet_tags.values)
    si = np.asarray(sheet_tags.indices)
    local_sheets, local_all = {}, {}
    for i in ring_ports:
        f_all = si[sv == SHEET_IFACE + i]
        local_all[i] = int(len(f_all))
        f = _owned(msh, f_all)
        local_sheets[i] = (f, _facet_coords(msh, f))
    gathered = comm.allgather({i: v[1] for i, v in local_sheets.items()})
    sheet = {i: np.concatenate([g[i] for g in gathered]) for i in ring_ports}
    n_owned = {i: len(sheet[i]) for i in ring_ports}
    n_local_copy = {i: comm.allreduce(local_all[i], op=MPI.SUM) for i in ring_ports}
    t1_bad = {i: (n_owned[i], n_gmsh.get(i)) for i in ring_ports if n_owned[i] != n_gmsh.get(i)}
    dphi = 2.0 * math.pi / LEG_COUNT
    area_dev, u_dev, zext_dev, chord_dev = {}, {}, {}, {}
    for i, xyz in sheet.items():
        area_dev[i] = abs(float(_area(xyz).sum()) / sheet_ref - 1.0)
        p = xyz.reshape(-1, 3)
        if len(p) == 0:
            u_dev[i] = zext_dev[i] = chord_dev[i] = float("inf")
            continue
        phi = (math.floor(math.atan2(p[:, 1].mean(), p[:, 0].mean()) / dphi) + 0.5) * dphi
        u = p[:, 0] * math.cos(phi) + p[:, 1] * math.sin(phi)
        v = -p[:, 0] * math.sin(phi) + p[:, 1] * math.cos(phi)
        u_dev[i] = float(np.max(np.abs(u / R - 1.0)))
        zext_dev[i] = abs((p[:, 2].max() - p[:, 2].min()) / width - 1.0)
        chord_dev[i] = abs((v.max() - v.min()) / chord - 1.0)

    # --- (B1), (B2): owned tag-3 tets, reduced.
    cv = np.asarray(cells.values)
    ci = np.asarray(cells.indices)
    n_owned_cells = msh.topology.index_map(tdim).size_local
    ph = ci[(cv == 3) & (ci < n_owned_cells)]
    if len(ph):
        g = dolfinx.cpp.mesh.entities_to_geometry(
            msh._cpp_object, tdim, np.asarray(ph, dtype=np.int32), False)
        X = msh.geometry.x[g[:, :4]]
        vol_local = float(np.abs(np.einsum(
            "ij,ij->i", X[:, 1] - X[:, 0],
            np.cross(X[:, 2] - X[:, 0], X[:, 3] - X[:, 0]))).sum() / 6.0)
        P = X.reshape(-1, 3)
        lo_local, hi_local = P.min(axis=0), P.max(axis=0)
    else:
        vol_local = 0.0
        lo_local, hi_local = np.full(3, np.inf), np.full(3, -np.inf)
    box_vol = comm.allreduce(vol_local, op=MPI.SUM)
    lo = np.array([comm.allreduce(float(x), op=MPI.MIN) for x in lo_local])
    hi = np.array([comm.allreduce(float(x), op=MPI.MAX) for x in hi_local])
    box_ref = Lx * Ly * Lz
    b1_dev = abs(box_vol / box_ref - 1.0)
    lo_ref = np.array([-0.5 * Lx, -0.5 * Ly + off, -0.5 * Lz])
    hi_ref = np.array([0.5 * Lx, 0.5 * Ly + off, 0.5 * Lz])
    b2_dev = float(max(np.max(np.abs(lo - lo_ref) / np.abs(lo_ref)),
                       np.max(np.abs(hi - hi_ref) / np.abs(hi_ref))))

    # --- (V1), (V2): CAD-internal partition and the can.
    cad_groups = dict(diag["cad_mass_by_group"])
    cad_sum = float(sum(cad_groups.values()))
    m_cond = float(diag["conductor_cad_mass_m3"])
    m_air = float(diag["air_domain_cad_mass_m3"])
    v1_dev = abs((cad_sum + m_cond) / m_air - 1.0)
    can_ref = math.pi * R_s**2 * 2.0 * z_e
    v2_dev = abs(m_air / can_ref - 1.0)

    # --- (K) census.
    census_fails = _census_checks(census, n_cells, ring_ports)

    # --- (c1) planted dropped facet on P33, by construction.
    p33 = ring_ports[0]
    have = comm.allgather(len(local_sheets[p33][0]) > 0)
    owner = have.index(True) if any(have) else -1
    xyz33 = local_sheets[p33][1]
    drop = comm.rank == owner
    c1_count = comm.allreduce(len(xyz33) - (1 if drop else 0), op=MPI.SUM)
    c1_area = comm.allreduce(float(_area(xyz33[1:] if drop else xyz33).sum()), op=MPI.SUM)
    c1_dev = abs(c1_area / sheet_ref - 1.0)

    rss = comm.allreduce(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0**2,
                         op=MPI.SUM)

    if comm.rank == 0:
        a0_rel = abs(n_cells / F_CLINICAL_R1_CELL_RECORD - 1.0)
        wa = max(area_dev, key=area_dev.get)
        wu = max(u_dev, key=u_dev.get)
        wz = max(zext_dev, key=zext_dev.get)
        wc = max(chord_dev, key=chord_dev.get)
        dup = sum(n_local_copy[i] - n_owned[i] for i in ring_ports)
        print(
            f"[GEO-38 s2] -n {comm.size}: {n_cells} cells; build {t_build:.1f} s (gmsh "
            f"{diag.get('mesh_wall_time_s', float('nan')):.1f} s); summed peak RSS {rss:.3f} GiB\n"
            f"[GEO-38 s2] census {census}\n"
            f"[GEO-38 s2] (T1) owned sheet facets == N_gmsh on {len(ring_ports) - len(t1_bad)}"
            f"/{len(ring_ports)} sheets; mismatches {t1_bad or 'none'}; N_gmsh "
            f"{min(n_gmsh.values())}..{max(n_gmsh.values())}, total {sum(n_gmsh.values())}; "
            f"owned total {sum(n_owned.values())}\n"
            f"[GEO-38 s2] (T2) chord*w = {chord:.12e}*{width:.12e} = {sheet_ref:.12e}; worst area "
            f"P{wa} {area_dev[wa]:.3e}; worst |u/R-1| P{wu} {u_dev[wu]:.3e}; worst z-extent/w "
            f"P{wz} {zext_dev[wz]:.3e}; worst chord extent P{wc} {chord_dev[wc]:.3e} "
            f"(band {STRIP_BAND:.0e})\n"
            f"[GEO-38 s2] (B1) box meshed volume {box_vol:.15e} vs {box_ref:.15e} m^3: rel "
            f"{b1_dev:.3e}\n"
            f"[GEO-38 s2] (B2) tag-3 extents lo {lo.tolist()} hi {hi.tolist()} vs lo "
            f"{lo_ref.tolist()} hi {hi_ref.tolist()}: worst rel {b2_dev:.3e}\n"
            f"[GEO-38 s2] (V1) CAD groups {cad_sum:.15e} + conductor {m_cond:.15e} vs air domain "
            f"{m_air:.15e} m^3: rel {v1_dev:.3e} (groups {cad_groups})\n"
            f"[GEO-38 s2] (V2) air domain {m_air:.15e} vs pi R_s^2 2z_e {can_ref:.15e}: rel "
            f"{v2_dev:.3e} (band {CAD_ANALYTIC_BAND:.0e})\n"
            f"[GEO-38 s2] (K) census failures {census_fails or 'none'}\n"
            f"[GEO-38 s2] (A0) cells {n_cells} vs F_CLINICAL_R1_CELL_RECORD "
            f"{F_CLINICAL_R1_CELL_RECORD}: rel {a0_rel:.4e} (CELL_COUNT_BAND {CELL_COUNT_BAND}; "
            f"{'within' if a0_rel < CELL_COUNT_BAND else 'outside: register a box record'}; "
            f"printed)\n"
            f"[GEO-38 s2] (c1) P{p33} one facet dropped on rank {owner}: count {c1_count} vs "
            f"N_gmsh-1 {n_gmsh[p33] - 1}; area rel {c1_dev:.3e} > {STRIP_BAND:.0e} (ceiling ≈ "
            f"1/{n_gmsh[p33]} = {1.0 / max(n_gmsh[p33], 1):.3e})\n"
            f"[GEO-38 s2] (c2) local-copy (owned + ghost) double count, printed only: "
            f"{sum(n_local_copy.values())} local copies vs {sum(n_owned.values())} owned "
            f"(excess {dup})",
            flush=True,
        )

    # Negative control first (by construction).
    assert owner >= 0
    assert c1_count == n_gmsh[p33] - 1, (c1_count, n_gmsh[p33])
    assert c1_dev > STRIP_BAND, c1_dev
    # (T1)
    assert len(n_gmsh) == 2 * LEG_COUNT, sorted(n_gmsh)
    assert not t1_bad, t1_bad
    # (T2)
    for d in (area_dev, u_dev, zext_dev, chord_dev):
        for i, dev in d.items():
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
