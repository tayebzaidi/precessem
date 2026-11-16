"""`WF-12` step 1 (§10 chain H12a, §9 item 129) — F-human as a copper hole
mesh: CAD identities and a one-drive cost probe.

Selected by environment only (``FEM_EM_WF12_STEP1=1``), never by ``-k``.

``test_wf12_step1_hole_mesh_identities`` (real build, ``-n 1`` — a failed gmsh
build deadlocks at ranks > 1, `GEO-23`) meshes the F-human rung (`GEO-25`
branch B, longitudinal ring sheets, `_params(F_HUMAN_RING_RADIUS,
scale_sizing=False)`) twice: ``as_hole=True`` and the solid control.

Anchors (asserted, every band imported):
  (i)   the 32 longitudinal ring sheets mesh ``chord · w`` from
        ``ring_port_layout`` (``ring_port_sheet_longitudinal_area_m2``) within
        ``SHEET_AREA_BAND`` 1e-9 (`test_birdcage_leg_offset.py:104`);
  (ii)  the CAD partition box − coil: hole groups + the solid's conductor CAD
        mass = the solid's groups within ``EXACT`` 1e-9, and hole groups =
        analytic box − coil within ``CAD_ANALYTIC_BAND`` 1e-7
        (`test_birdcage_conductor_hole.py:117,123`); the analytic box is the
        generator's own extents (`io/mesh.py:4232–4233`) evaluated on the
        call's kwargs and the realised ``port_box_size_m[1]``;
  (iii) the meshed cavity wall (tag ``BIRDCAGE_CONDUCTOR_SURFACE_TAG``) equals
        the solid route's whole conductor interface (tag 1 against every other
        tag) within the fixture band ``F_HUMAN_WALL_AREA_RTOL`` 5e-5 (ruled
        2026-10-01 18:00 on item 147's probe; the two-torus
        ``INTERFACE_AREA_RTOL`` 1e-5 is printed beside it);
  (iii-b) the ring-port-class part of that wall (tags 117–148 / 217–248) is
        equal on the two routes within 1e-8 relative, with equal facet counts;
  (iv)  the hole census has no tag 1, has air, phantom and all 64 ring-port
        halves, and sums (owned, allreduced) to the global cell count.
Negative control (asserted, binary): the solid route reproduces
``F_HUMAN_BRANCH_B_CELL_RECORD`` within ``CELL_COUNT_BAND`` with tag 1
present, and the hole mesh has strictly fewer cells.

``test_wf12_step1_one_drive_price`` (complex, ``-n 8``): the hole built by
``_build_f_human(as_hole=True)``, one P17 drive with all 32 ring ports at
50 Ω, the copper (σ = 5.8e7) Leontovich wall on the cavity
(`test_th14_birdcage_copper._leontovich_term`), the outer box pinned only.
Re-asserts (iv) under reduction and asserts
  (P)   P_src = P_phantom + P_(Ω∖phantom) + P_surf + Σ P_sheets within
        ``DISCRETE_IDENTITY_RTOL`` 1e-6 (`test_birdcage_power_identity.py:227`).
Predicted, printed only: dropping P_surf from (P) misses by the copper share
(`ANS-6` F-small 64 MHz read 0.448).  Price (cells, mesh time, unknowns,
factorisation, RSS per rank and summed) printed against the item's written
predictions.

Scope: one drive, degree 1, 64 MHz. Not the 32×32, a PEC column, C16
spreads or |B1+| (H12b).
"""

from __future__ import annotations

import os
import resource
import time

import numpy as np
import pytest
from mpi4py import MPI

from tests.complex_mode import complex_only

STEP1_ENV = "FEM_EM_WF12_STEP1"
DRIVEN = "P17"
COPPER_SIGMA = 5.8e7

# Predictions written before the first run (§9 item 129, "Printed against
# predictions written first").
PREDICTED_HOLE_CELLS = (300_000, 420_000)  # central ≈ 350 k
PREDICTED_MESH_S = 114.6  # ≈ the solid's
PREDICTED_UNKNOWNS = 420_000
PREDICTED_FACTOR_S = (20.0, 35.0)
PREDICTED_RSS_GIB = (7.0, 9.0)
# The H12b price line: 2 factorisations + 62 held drives + post-processing.
H12B_TIME_CEILING_S = 540.0
H12B_RSS_CEILING_GIB = 64.0
# `ANS-6` F-small 64 MHz copper share — the predicted miss of the P_surf-blind
# reading (printed, never asserted).
PREDICTED_COPPER_SHARE = 0.448

AIR_CELL_TAG = 2
SOLID_CONDUCTOR_IFACE = 900

# Anchor (iii) band for THIS fixture (F-human, 1.6 mm conductor_resolution),
# re-registered by the 2026-10-01 18:00 review (§9 rule (f); the 10-03 weekly
# ratifies or reverts; brief docs/planning/chunks/WF-12.md "Ruled 2026-10-01
# 18:00").  Measured by item 147's probe (six -n 1 logs:
# 20261001T214637Z_WF-12-probe-wall-h1p6-hole.log,
# 20261001T214837Z_WF-12-probe-wall-h1p6-solid.log,
# 20261001T215102Z_WF-12-probe-wall-h2p4-hole.log,
# 20261001T215223Z_WF-12-probe-wall-h2p4-solid.log,
# 20261001T215402Z_WF-12-probe-wall-h3p2-hole.log,
# 20261001T215509Z_WF-12-probe-wall-h3p2-solid.log): hole-minus-solid wall
# area D = +2.40e-6 / +4.84e-6 / -5.51e-6 m^2 at 1.6 / 2.4 / 3.2 mm changes
# sign, sits 100 % on the air-facing wall, follows the air-class facet-count
# difference (+32 / +68 / -50), and the ring-port wall is identical on both
# routes -- two independent inscribed triangulations of a 0.208 m^2 curved
# wall, not a missing face.  |D|/A_solid = 1.155e-5 / 2.341e-5 / 2.684e-5, so
# the two-torus INTERFACE_AREA_RTOL 1e-5 (3.7e-7 measured there) cannot be met
# by construction here.  5e-5 = 1.9x the worst rung, 4.3x the 1.6 mm rung used;
# it still catches any missing face > 5e-5 * 0.2078 m^2 ≈ 1.0e-5 m^2.  The
# two-torus / F-small band is unchanged.
F_HUMAN_WALL_AREA_RTOL = 5.0e-5
# (iii-b), same ruling: the ring-port-class wall area of the two routes is
# equal (measured 5.4e-10 relative at 1.6 mm, identical to 11 digits and in
# facet count at 2.4 / 3.2 mm), with equal facet counts.
RING_PORT_WALL_RTOL = 1.0e-8


# The neighbour-class split below is lifted from the measurement probe
# tests/mesh/probe_wf12_cavity_wall_area.py (§9 item 147: `_class_of`,
# `_cell_value`, `_owned_facet_split`), reduced to the two classes the (iii-b)
# anchor reads.  Not imported: the probe is imported by nothing.
def _wall_class(tag: int) -> str:
    if tag == AIR_CELL_TAG:
        return "air"
    if 117 <= tag <= 148 or 217 <= tag <= 248:
        return "ring_port"
    return f"other({tag})"


def _cell_value(msh, cell_tags):
    """Per local cell (owned + ghost) tag value, ghosts filled by scatter."""
    from dolfinx import fem

    space = fem.functionspace(msh, ("DG", 0))
    f = fem.Function(space)
    f.x.array[:] = 0.0
    c2d = space.dofmap.list.reshape(-1)
    f.x.array[c2d[cell_tags.indices]] = cell_tags.values
    f.x.scatter_forward()
    return np.rint(np.real(f.x.array[c2d])).astype(np.int32)


def _wall_split(msh, facet_indices, neighbour_tag, comm):
    """Per-class (area, facet count) over OWNED facets, allgathered."""
    import dolfinx

    fdim = msh.topology.dim - 1
    n_owned = msh.topology.index_map(fdim).size_local
    keep = facet_indices < n_owned
    facets = facet_indices[keep]
    nbr = neighbour_tag[keep]
    geom_ids = dolfinx.cpp.mesh.entities_to_geometry(
        msh._cpp_object, fdim, facets.astype(np.int32), False)
    x = msh.geometry.x
    p0, p1, p2 = x[geom_ids[:, 0]], x[geom_ids[:, 1]], x[geom_ids[:, 2]]
    area = 0.5 * np.linalg.norm(np.cross(p1 - p0, p2 - p0), axis=1)
    local = {}
    for t, a in zip(nbr, area):
        s = local.setdefault(_wall_class(int(t)), [0.0, 0])
        s[0] += float(a)
        s[1] += 1
    out = {}
    for d in comm.allgather(local):
        for c, (a, n) in d.items():
            s = out.setdefault(c, [0.0, 0])
            s[0] += a
            s[1] += n
    return out


def _hole_wall_split(msh, cells, facets, cavity_tag, comm):
    tdim = msh.topology.dim
    f2c = msh.topology.connectivity(tdim - 1, tdim)
    offsets, links = f2c.offsets, f2c.array
    cv = _cell_value(msh, cells)
    idx = np.asarray(facets.indices)[np.asarray(facets.values) == cavity_tag]
    return _wall_split(msh, idx, cv[links[offsets[idx]]], comm)


def _solid_wall_split(msh, cells, iface, iface_tag, conductor_tag, comm):
    tdim = msh.topology.dim
    f2c = msh.topology.connectivity(tdim - 1, tdim)
    offsets, links = f2c.offsets, f2c.array
    cv = _cell_value(msh, cells)
    idx = np.asarray(iface.indices)[np.asarray(iface.values) == iface_tag]
    a = cv[links[offsets[idx]]]
    b = cv[links[offsets[idx] + 1]]
    return _wall_split(msh, idx, np.where(a == conductor_tag, b, a), comm)


def _on() -> bool:
    return os.environ.get(STEP1_ENV, "").strip() == "1"


def _skip_unless_on():
    if not _on():
        pytest.skip(f"{STEP1_ENV} unset: WF-12 step 1 is selected by environment only")


def _rss_gib(comm):
    """Peak RSS of this rank (GiB) and the summed peak across ranks."""
    mine = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0**2
    per_rank = comm.gather(mine, root=0)
    return mine, per_rank, comm.allreduce(mine, op=MPI.SUM)


def _analytic_box(kwargs, port_dy):
    """`io/mesh.py:4232–4233`'s own extents on this call's kwargs."""
    leg_radius_eff = 0.5 * kwargs["leg_width"]
    radial = (kwargs["ring_radius"] + max(leg_radius_eff, kwargs["ring_minor_radius"])
              + port_dy + kwargs["air_padding"])
    z = (max(0.5 * kwargs["coil_length"],
             0.5 * kwargs["leg_spacing"] + kwargs["ring_minor_radius"],
             0.5 * kwargs["phantom_height"]) + kwargs["air_padding"])
    return 8.0 * radial * radial * z


def _owned_census(msh, tags, comm):
    from tests.mesh.probe_two_torus_conductor_hole import _census

    tdim = msh.topology.dim
    n_owned = msh.topology.index_map(tdim).size_local
    owned = np.asarray(tags.indices) < n_owned
    return _census(np.asarray(tags.values)[owned], comm)


def _census_checks(census, n_cells, ring_ports):
    """(iv): returns a list of failure strings (empty ⇒ green)."""
    from tests.mesh.test_birdcage_port_sheets import PORT_LOWER, PORT_UPPER
    from tests.validation.test_port_birdcage_lumped_column import (
        CONDUCTOR_CELL_TAG,
        PHANTOM_CELL_TAG,
    )

    fails = []
    if CONDUCTOR_CELL_TAG in census:
        fails.append(f"tag {CONDUCTOR_CELL_TAG} present ({census[CONDUCTOR_CELL_TAG]} cells)")
    need = {AIR_CELL_TAG, PHANTOM_CELL_TAG}
    need |= {PORT_LOWER + i for i in ring_ports} | {PORT_UPPER + i for i in ring_ports}
    missing = sorted(need - set(census))
    if missing:
        fails.append(f"missing tags {missing}")
    if sum(census.values()) != n_cells:
        fails.append(f"census sums to {sum(census.values())} != {n_cells}")
    return fails


def _ring_sheet_areas(msh, cell_tags, ring_ports, comm):
    from fem_em_solver.io.mesh import _interface_facet_tags
    from tests.mesh.probe_two_torus_conductor_hole import _facet_area
    from tests.mesh.test_birdcage_port_sheets import PORT_LOWER, PORT_UPPER, SHEET_IFACE

    tags = _interface_facet_tags(
        msh, cell_tags,
        {SHEET_IFACE + i: (PORT_LOWER + i, PORT_UPPER + i) for i in ring_ports},
    )
    return {i: _facet_area(msh, tags, SHEET_IFACE + i, comm) for i in ring_ports}


def test_wf12_step1_hole_mesh_identities():
    """(i)–(iv) and the solid negative control, real build, ``-n 1``."""
    _skip_unless_on()
    from fem_em_solver.io.mesh import (
        BIRDCAGE_CONDUCTOR_SURFACE_TAG,
        MeshGenerator,
        _interface_facet_tags,
    )
    from tests.mesh.probe_two_torus_conductor_hole import _exterior_facet_area, _facet_area
    from tests.mesh.test_birdcage_conductor_hole import (
        CAD_ANALYTIC_BAND,
        EXACT,
        INTERFACE_AREA_RTOL,
    )
    from tests.mesh.test_birdcage_leg_offset import SHEET_AREA_BAND
    from tests.mesh.test_birdcage_port_sheet_prerequisite import CELL_COUNT_BAND
    from tests.validation.test_birdcage_f_human_rung import (
        F_HUMAN_BRANCH_B_CELL_RECORD,
        F_HUMAN_RING_RADIUS,
        LEG_COUNT,
        _params,
        _ring_ports,
    )
    from tests.validation.test_port_birdcage_lumped_column import CONDUCTOR_CELL_TAG

    comm = MPI.COMM_WORLD
    ring_ports = _ring_ports()
    assert len(ring_ports) == 2 * LEG_COUNT

    def build(as_hole):
        kwargs = _params(F_HUMAN_RING_RADIUS, scale_sizing=False)
        kwargs["ring_sheet_orientation"] = "longitudinal"
        if as_hole:
            kwargs["as_hole"] = True
        t0 = time.perf_counter()
        msh, cells, facets, diag = MeshGenerator.birdcage_port_domain(
            comm=comm, return_diagnostics=True, **kwargs)
        t_build = time.perf_counter() - t0
        tdim = msh.topology.dim
        msh.topology.create_connectivity(tdim - 1, tdim)
        msh.topology.create_entity_permutations()
        n_cells = int(msh.topology.index_map(tdim).size_global)
        return {"msh": msh, "cells": cells, "facets": facets, "diag": diag, "kwargs": kwargs,
                "t_build": t_build, "n_cells": n_cells,
                "census": _owned_census(msh, cells, comm)}

    hole = build(True)
    if comm.rank == 0:
        print(f"\n[WF-12 step1] hole: {hole['n_cells']} cells, build {hole['t_build']:.1f} s "
              f"(gmsh mesh_wall_time_s {hole['diag'].get('mesh_wall_time_s', float('nan')):.1f}), "
              f"-n {comm.size}; census {hole['census']}", flush=True)
    hole["sheet_area"] = _ring_sheet_areas(hole["msh"], hole["cells"], ring_ports, comm)
    hole["cavity_area"] = _exterior_facet_area(
        hole["msh"], hole["facets"], BIRDCAGE_CONDUCTOR_SURFACE_TAG, comm)
    hole["split"] = _hole_wall_split(
        hole["msh"], hole["cells"], hole["facets"], BIRDCAGE_CONDUCTOR_SURFACE_TAG, comm)
    _mine, _per, hole_rss = _rss_gib(comm)

    solid = build(False)
    solid_tags = sorted(t for t in solid["census"] if t != CONDUCTOR_CELL_TAG)
    iface = _interface_facet_tags(
        solid["msh"], solid["cells"],
        {SOLID_CONDUCTOR_IFACE: tuple((CONDUCTOR_CELL_TAG, t) for t in solid_tags)})
    solid["interface_area"] = _facet_area(solid["msh"], iface, SOLID_CONDUCTOR_IFACE, comm)
    solid["split"] = _solid_wall_split(
        solid["msh"], solid["cells"], iface, SOLID_CONDUCTOR_IFACE, CONDUCTOR_CELL_TAG, comm)
    _mine, _per, both_rss = _rss_gib(comm)

    layout = hole["diag"]["ring_port_layout"]
    sheet_analytic = float(layout["ring_port_sheet_longitudinal_area_m2"])
    chord = float(layout["ring_port_gap_chord_m"])
    width = float(layout["ring_port_box_width_m"])
    sheet_dev = {i: abs(a / sheet_analytic - 1.0) for i, a in hole["sheet_area"].items()}

    coil = float(solid["diag"]["cad_mass_by_group"]["conductor"])
    hole_sum = float(sum(hole["diag"]["cad_mass_by_group"].values()))
    solid_sum = float(sum(solid["diag"]["cad_mass_by_group"].values()))
    box = _analytic_box(hole["kwargs"], float(hole["diag"]["port_box_size_m"][1]))
    cut_ratio = (hole_sum + coil) / solid_sum
    hole_ratio = hole_sum / (box - coil)
    solid_ratio = solid_sum / box
    cavity_dev = abs(hole["cavity_area"] / solid["interface_area"] - 1.0)
    empty = [0.0, 0]
    h_ring = hole["split"].get("ring_port", empty)
    s_ring = solid["split"].get("ring_port", empty)
    h_air = hole["split"].get("air", empty)
    s_air = solid["split"].get("air", empty)
    ring_dev = (abs(h_ring[0] / s_ring[0] - 1.0) if s_ring[0] > 0.0 else float("inf"))
    a_cad = float(hole["diag"].get("conductor_cavity_cad", {}).get("area_m2", float("nan")))
    census_fails = _census_checks(hole["census"], hole["n_cells"], ring_ports)
    record_ratio = solid["n_cells"] / F_HUMAN_BRANCH_B_CELL_RECORD

    if comm.rank == 0:
        worst = max(sheet_dev, key=sheet_dev.get)
        print(
            f"[WF-12 step1] solid: {solid['n_cells']} cells (record "
            f"{F_HUMAN_BRANCH_B_CELL_RECORD}, ratio {record_ratio:.6f}, band {CELL_COUNT_BAND}), "
            f"build {solid['t_build']:.1f} s; census {solid['census']}\n"
            f"[WF-12 step1] (i) 32 ring sheets vs chord*w = {chord:.9e}*{width:.9e} = "
            f"{chord * width:.9e} (layout {sheet_analytic:.9e}); worst P{worst} "
            f"{hole['sheet_area'][worst]:.9e} rel {sheet_dev[worst]:.3e} (band {SHEET_AREA_BAND:.0e})\n"
            f"[WF-12 step1] (ii) hole CAD {hole_sum:.12e} + coil {coil:.12e} vs solid CAD "
            f"{solid_sum:.12e}: ratio-1 {cut_ratio - 1.0:+.3e} (band {EXACT:.0e}); box {box:.12e}: "
            f"hole/(box-coil)-1 {hole_ratio - 1.0:+.3e}, solid/box-1 {solid_ratio - 1.0:+.3e} "
            f"(band {CAD_ANALYTIC_BAND:.0e})\n"
            f"[WF-12 step1] (iii) cavity {hole['cavity_area']:.9e} vs solid conductor interface "
            f"{solid['interface_area']:.9e} (tag 1 vs {solid_tags[:3]}... {len(solid_tags)} tags): "
            f"rel {cavity_dev:.3e} (F-human band {F_HUMAN_WALL_AREA_RTOL:.0e}; two-torus band "
            f"{INTERFACE_AREA_RTOL:.0e})\n"
            f"[WF-12 step1] (iii-b) ring-port wall hole {h_ring[0]:.10e} ({h_ring[1]} facets) vs "
            f"solid {s_ring[0]:.10e} ({s_ring[1]} facets): rel {ring_dev:.3e} (band "
            f"{RING_PORT_WALL_RTOL:.0e}), facet counts equal {h_ring[1] == s_ring[1]}\n"
            f"[WF-12 step1] air wall hole {h_air[0]:.10e} ({h_air[1]} facets) vs solid "
            f"{s_air[0]:.10e} ({s_air[1]} facets): air-class facet difference hole-solid "
            f"{h_air[1] - s_air[1]:+d}, area difference {h_air[0] - s_air[0]:+.6e} m^2; "
            f"classes hole {sorted(hole['split'])} solid {sorted(solid['split'])}\n"
            f"[WF-12 step1] A_CAD {a_cad:.10e} m^2: hole A/A_CAD-1 "
            f"{hole['cavity_area'] / a_cad - 1.0:+.6e}, solid A/A_CAD-1 "
            f"{solid['interface_area'] / a_cad - 1.0:+.6e}\n"
            f"[WF-12 step1] (iv) census failures: {census_fails or 'none'}\n"
            f"[WF-12 step1] PRICE hole cells {hole['n_cells']} (PREDICTED "
            f"{PREDICTED_HOLE_CELLS[0]}-{PREDICTED_HOLE_CELLS[1]}, ≈350 k), hole/solid "
            f"{hole['n_cells'] / solid['n_cells']:.4f}; hole build {hole['t_build']:.1f} s "
            f"(PREDICTED ≈{PREDICTED_MESH_S} s); peak RSS after hole {hole_rss:.3f} GiB, after "
            f"both {both_rss:.3f} GiB (summed over {comm.size} rank(s))",
            flush=True,
        )

    # Negative control first (binary).
    assert abs(record_ratio - 1.0) < CELL_COUNT_BAND, (solid["n_cells"], record_ratio)
    assert CONDUCTOR_CELL_TAG in solid["census"], solid["census"]
    assert hole["n_cells"] < solid["n_cells"], (hole["n_cells"], solid["n_cells"])
    # (i)
    assert len(hole["sheet_area"]) == 2 * LEG_COUNT
    for i, dev in sheet_dev.items():
        assert dev < SHEET_AREA_BAND, (i, hole["sheet_area"][i], sheet_analytic, dev)
    # (ii)
    assert abs(cut_ratio - 1.0) < EXACT, (hole_sum, coil, solid_sum, cut_ratio)
    assert abs(hole_ratio - 1.0) < CAD_ANALYTIC_BAND, (hole_sum, box, coil, hole_ratio)
    assert abs(solid_ratio - 1.0) < CAD_ANALYTIC_BAND, (solid_sum, box, solid_ratio)
    # (iii)
    assert cavity_dev < F_HUMAN_WALL_AREA_RTOL, (hole["cavity_area"], solid["interface_area"])
    # (iii-b)
    assert ring_dev <= RING_PORT_WALL_RTOL, (h_ring, s_ring, ring_dev)
    assert h_ring[1] == s_ring[1] and h_ring[1] > 0, (h_ring, s_ring)
    # (iv)
    assert not census_fails, census_fails


@complex_only
def test_wf12_step1_one_drive_price():
    """(iv) under reduction and (P) on one copper P17 drive, ``-n 8``."""
    _skip_unless_on()
    import ufl
    from dolfinx import fem

    from fem_em_solver.core.cavity import surface_resistance_ohm
    from fem_em_solver.ports.lumped import lumped_port_bilinear_term, lumped_port_linear_term
    from tests.validation.test_birdcage_f_human_rung import (
        F_HUMAN_BRANCH_B_CELL_RECORD,
        _ring_ports,
    )
    from tests.validation.test_birdcage_power_identity import (
        DISCRETE_IDENTITY_RTOL,
        _sheet_field_dissipation_w,
        _source_power_w,
    )
    from tests.validation.test_port_birdcage_lumped_column import PHANTOM_CELL_TAG
    from tests.validation.test_th14_birdcage_copper import _leontovich_term, _surface_loss_w
    from tests.validation.test_wf7_f_human_b1_quadrature import DEGREE, _build_f_human

    comm = MPI.COMM_WORLD
    t_window = time.perf_counter()
    t0 = time.perf_counter()
    built = _build_f_human(comm, as_hole=True)
    t_build = time.perf_counter() - t0
    ctx = built["ctx"]
    msh, tags_f, omega = ctx["msh"], ctx["tags_f"], ctx["omega"]
    census = _owned_census(msh, ctx["cell_tags"], comm)
    census_fails = _census_checks(census, built["cells"], _ring_ports())
    oc = built["outer_census"]
    _v = fem.functionspace(msh, ("N1curl", DEGREE))
    unknowns = int(_v.dofmap.index_map.size_global * _v.dofmap.index_map_bs)
    del _v
    if comm.rank == 0:
        print(f"\n[WF-12 step1 W2] hole: {built['cells']} cells (solid record "
              f"{F_HUMAN_BRANCH_B_CELL_RECORD}), build {t_build:.1f} s at -n {comm.size}; "
              f"unknowns {unknowns} (PREDICTED ≈{PREDICTED_UNKNOWNS}); census {census}; "
              f"(iv) failures {census_fails or 'none'}; outer-box census {oc}", flush=True)
    assert not census_fails, census_fails
    assert oc["tag_collision"] == 0 and oc["cavity"] > 0, oc
    assert oc["cavity_exterior"] == oc["cavity"] and oc["outer"] + oc["cavity"] == oc["exterior"], oc

    z_s = (1.0 + 1.0j) * surface_resistance_ohm(omega, COPPER_SIGMA)
    wall = _leontovich_term(msh, built["facets"], z_s, omega)
    port_sheets = [sp.sheet(driven=(sp.port_id == DRIVEN)) for sp in ctx["specs"]]
    driven_sheet = next(s for s in port_sheets if s.port_id == DRIVEN)
    bilinear = [wall] + [
        (lambda trial, test, _s=s: lumped_port_bilinear_term(
            msh, tags_f, _s, trial, test, omega_rad_per_s=omega))
        for s in port_sheets
    ]
    linear = [lambda test, _s=driven_sheet: lumped_port_linear_term(
        msh, tags_f, _s, test, omega_rad_per_s=omega)]
    comm.Barrier()
    t0 = time.perf_counter()
    fields = ctx["solver"].solve(
        current_density=None, project_source=False,
        extra_bilinear_terms=bilinear, extra_linear_terms=linear)
    comm.Barrier()
    t_solve = time.perf_counter() - t0
    rss_mine, rss_per, rss_sum = _rss_gib(comm)

    e = fields.e_complex
    sigma_f = fields.sigma_field
    dx = ufl.Measure("dx", domain=msh, subdomain_data=ctx["cell_tags"])

    def _vol(measure):
        form = fem.form(0.5 * sigma_f * ufl.inner(e, e) * measure)
        return float(np.real(comm.allreduce(fem.assemble_scalar(form), op=MPI.SUM)))

    acct = {"mesh": msh, "facet_tags": tags_f}
    p_sheets = float(sum(_sheet_field_dissipation_w(acct, s, e, omega) for s in port_sheets))
    p_src = float(_source_power_w(acct, driven_sheet, e, omega))
    p_total = _vol(dx)
    p_phantom = _vol(dx(PHANTOM_CELL_TAG))
    p_rest = p_total - p_phantom
    p_surf = _surface_loss_w(msh, built["facets"], z_s, e, comm)
    residual = abs(p_src - p_sheets - p_phantom - p_rest - p_surf) / abs(p_src)
    blind = abs(p_src - p_sheets - p_phantom - p_rest) / abs(p_src)
    p_in = p_src - p_sheets
    t_total = time.perf_counter() - t_window

    if comm.rank == 0:
        n = comm.size
        ok_t = t_solve * 2 + 62 * 3.0  # printed estimate only, see below
        print(
            f"[WF-12 step1 W2] Z_s = {z_s.real:.6e}(1+j) ohm; P17 solve (assembly + "
            f"factorisation + solve) {t_solve:.1f} s (PREDICTED {PREDICTED_FACTOR_S[0]:.0f}-"
            f"{PREDICTED_FACTOR_S[1]:.0f} s)\n"
            f"[WF-12 step1 W2] (P) P_src {p_src:.9e} W, sum P_sheet {p_sheets:.9e} W, P_phantom "
            f"{p_phantom:.9e} W, P_(Omega\\phantom) {p_rest:.9e} W, P_surf {p_surf:.9e} W; "
            f"residual/P_src {residual:.3e} (ASSERTED <= {DISCRETE_IDENTITY_RTOL:g})\n"
            f"[WF-12 step1 W2] NEG (PREDICTED, printed): P_surf-blind miss {blind:.4e} vs predicted "
            f"≈ copper share {PREDICTED_COPPER_SHARE} (ANS-6 F-small); measured P_surf/P_in "
            f"{p_surf / p_in:.6e}, P_phantom/P_in {p_phantom / p_in:.6e}\n"
            f"[WF-12 step1 W2] RSS per rank (GiB) "
            f"{', '.join(f'{r:.3f}' for r in rss_per)}; summed {rss_sum:.3f} GiB (PREDICTED "
            f"{PREDICTED_RSS_GIB[0]:.0f}-{PREDICTED_RSS_GIB[1]:.0f} GiB) at -n {n}\n"
            f"[WF-12 step1 W2] window {t_total:.1f} s; H12b price line (2 factorisations + 62 "
            f"held drives + post) vs <= {H12B_TIME_CEILING_S:.0f} s and <= "
            f"{H12B_RSS_CEILING_GIB:.0f} GiB at -n 8: one factorised drive here {t_solve:.1f} s, "
            f"summed RSS {rss_sum:.3f} GiB (held-drive cost not measured in this window; "
            f"naive 2x factor + 62x3 s estimate {ok_t:.0f} s)",
            flush=True,
        )
    assert p_surf > 0.0, p_surf
    assert residual <= DISCRETE_IDENTITY_RTOL, residual
