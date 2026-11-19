"""`GEO-36` step 8 (§9 item 530, §10 I33): the kept 1 mm x 30 mm rod in
F-clinical's ASTM box — CAD and mesh identities only, no solve.

`_box_kwargs()` (item 470's r1 kwargs with the ASTM box, `GEO-38` step 2) plus
``implant_rod = dict(ROD, keep_cells=True, centre=(x_rod, off, 0))``, axis ẑ —
the box's 650 mm axis (C2's coil frame: Lz = 0.65 on the coil axis, Lx = 0.42
horizontal, Ly = 0.09 the gel depth). ``x_rod`` defaults to Lx/4 = 0.105 m
toward the +x wall (review ruling 7, 2026-10-18 03:00): a modelling choice, no
ASTM-placement claim. Brief: docs/planning/chunks/GEO-36.md, "§9 item 530".

The build test is selected by environment only (``FEM_EM_GEO36_ROD_BOX=1``),
never by ``-k``; it skips when unset. The refusal tests are pure and always run.

Anchors (asserted at the window's width; owned cells / facets, reduced):
  (G0) cells < F_CLINICAL_ROD_CELL_GUARD before any facet reading;
  (G1) rod_mass = kept_mass = V_rod and OCC wall area = A_rod at STRIP_BAND;
  (G2) CAD(403) = V_rod, CAD(3) = L³ − V_rod, CAD(3)+CAD(403) = L³ at STRIP_BAND;
  (G3) owned tag-3 + tag-403 tet volumes = L³, tag-3 extents = box at STRIP_BAND;
  (G4) Σ cad_mass_by_group + conductor_cad_mass_m3 = can at STRIP_BAND;
  (G5) 402 interior, between tag 3 and 403, meshed dS(402)/A_rod in
       TERMINAL_AREA_BAND;
  (G6) OCC wall centroid at the rod centre to EXACT·ℓ;
  (G7) census clean, tags 3 and 403 present, no tag 1, implant_spheroid_cad None;
  (G8) per ring sheet owned facets = N_gmsh, area = chord·w at STRIP_BAND;
  (G9) width invariance: the printed ``(G9) ints`` line equal at -n 4 / -n 2;
  (G10) refusals (pure).
Negative controls: (N1r) C2's comparand L³ misses CAD(3) by V/(L³ − V), asserted
> STRIP_BAND, predicted 4.006e-6; (N2r) the pre-repair conductor reading misses
(G4) by V/can, asserted > STRIP_BAND, predicted 2.366e-7; (N3r) the y-frame
misreading (0, 0.105, 0) accepted by the cylinder test, refused by the box test
(both asserted, W0); (N4r) tag-403 cells, meshed 403 / V, 402 / A and the
conductor mass against C2's are printed only.
"""

from __future__ import annotations

import os
import resource
import time
from typing import Optional

import numpy as np
import pytest
from mpi4py import MPI

from tests.mesh.test_birdcage_implant_rod import (
    ROD,
    ROD_LENGTH,
    ROD_RADIUS,
    _interior_facet_check,
    _rod_area,
    _rod_volume,
)
from tests.validation.gmsh_stop_rule import arm_gmsh_stop_rule, gmsh_stop_armed
from tests.validation.test_geo38_f_clinical_box import _box_kwargs
from tests.validation.test_geo38_f_clinical_build import (
    LEG_COUNT,
    _area,
    _facet_coords,
    _owned,
    _refuse_before_gmsh,
    _ring_ports,
)

ROD_BOX_ENV = "FEM_EM_GEO36_ROD_BOX"

#: Guard, checked before any facet reading: predicted 1 149 150 + [50 k, 100 k]
#: plus headroom (the disc's guard, item 482).
F_CLINICAL_ROD_CELL_GUARD = 1_400_000
#: C2's conductor CAD mass, a printed comparand only
#: (`20261015T185056Z_GEO-38-s2-W3-n4.log:11903`).
C2_CONDUCTOR_CAD_MASS_M3 = 1.994040484843362e-03
#: Rod x offset [m]; ``None`` ⇒ 0.25 * Lx from `_box_kwargs()`. Set only with a
#: cited source (the `TABLE_OFFSET_M` precedent).
ROD_OFFSET_M: Optional[float] = None

PHANTOM_CELL_TAG = 3

#: (N3r): the literal "±y wall by a quarter of the width" read in C2's frame.
Y_MISREAD_CENTRE = (0.0, 0.105, 0.0)


def _say(msg):
    if MPI.COMM_WORLD.rank == 0:
        print(msg, flush=True)


def _rod_centre():
    box = _box_kwargs()["phantom_box"]
    x = 0.25 * float(box["size_m"][0]) if ROD_OFFSET_M is None else float(ROD_OFFSET_M)
    return (x, float(box["table_offset_m"]), 0.0)


def _kwargs():
    kw = _box_kwargs()
    kw["implant_rod"] = dict(ROD, keep_cells=True, centre=_rod_centre())
    return kw


# --- (G10) / (N3r): pure refusals, always run. ------------------------------


def test_n3r_y_misreading_cylinder_accepts_box_refuses(monkeypatch):
    from fem_em_solver.io.mesh import MeshGenerator

    box = _box_kwargs()["phantom_box"]
    Lx, Ly, Lz = box["size_m"]
    off = box["table_offset_m"]
    R = float(np.hypot(0.5 * Lx, 0.5 * Ly + abs(off)))
    req = dict(ROD, keep_cells=True, centre=Y_MISREAD_CENTRE)
    # Cylinder test (phantom_box=None): accepted.
    out = MeshGenerator._normalise_implant_rod(
        dict(req), phantom_radius=R, phantom_height=Lz, as_hole=True
    )
    margin = R - (float(np.hypot(*Y_MISREAD_CENTRE[:2])) + ROD_RADIUS)
    breach = abs(Y_MISREAD_CENTRE[1] - off) + ROD_RADIUS - 0.5 * Ly
    _say(
        f"\n[GEO-36 s8] (N3r) R {R:.6f} m, H {Lz} m: cylinder test accepts with margin "
        f"{margin:.6f} m; box breach |c_y - off| + r - Ly/2 = {breach:.6e} m"
    )
    assert out["keep_cells"] is True
    assert margin > 0.0
    assert breach > 0.06  # 0.106 - 0.045 = 0.061 m
    # Box test: refused before gmsh, through the public path.
    kw = _box_kwargs()
    kw["implant_rod"] = dict(req)
    with pytest.raises(
        ValueError, match="implant_rod is not strictly inside the phantom_box"
    ) as exc:
        _refuse_before_gmsh(monkeypatch, **kw)
    _say(f"[GEO-36 s8] (N3r) box refusal: {exc.value}")


def test_g10_coated_rod_with_box_under_as_hole_refused(monkeypatch):
    kw = _box_kwargs()
    assert kw.get("as_hole") is True
    kw["implant_rod"] = dict(
        ROD, keep_cells=True, centre=_rod_centre(), coat_radius=2.0 * ROD_RADIUS
    )
    with pytest.raises(
        ValueError, match="with as_hole unless keep_cells is True and no coat_radius"
    ):
        _refuse_before_gmsh(monkeypatch, **kw)


def test_g10_hole_rod_with_box_under_as_hole_refused(monkeypatch):
    kw = _box_kwargs()
    assert kw.get("as_hole") is True
    kw["implant_rod"] = dict(ROD, keep_cells=False, centre=_rod_centre())
    with pytest.raises(
        ValueError, match="with as_hole unless keep_cells is True and no coat_radius"
    ):
        _refuse_before_gmsh(monkeypatch, **kw)


def test_g10_kept_rod_with_box_accepted():
    from fem_em_solver.io.mesh import MeshGenerator

    box = _box_kwargs()["phantom_box"]
    Lx, Ly, Lz = box["size_m"]
    R = float(np.hypot(0.5 * Lx, 0.5 * Ly + abs(box["table_offset_m"])))
    out = MeshGenerator._normalise_implant_rod(
        dict(ROD, keep_cells=True, centre=_rod_centre()),
        phantom_radius=R,
        phantom_height=Lz,
        as_hole=True,
        phantom_box=dict(box),
    )
    assert out["keep_cells"] is True
    assert out["centre"] == _rod_centre()
    assert abs(out["centre"][0] - 0.105) < 1e-15  # Lx/4 for Lx = 0.42


# --- The build (environment-selected). ---------------------------------------


def test_geo36_step8_rod_in_astm_box(monkeypatch):
    if os.environ.get(ROD_BOX_ENV, "").strip() != "1":
        pytest.skip(f"{ROD_BOX_ENV} unset: GEO-36 step 8 is selected by environment only")
    import dolfinx

    from fem_em_solver.io.mesh import (
        BIRDCAGE_IMPLANT_CELL_TAG,
        BIRDCAGE_IMPLANT_SURFACE_TAG,
        MeshGenerator,
        _interface_facet_tags,
    )
    from tests.mesh.test_birdcage_conductor_hole import EXACT
    from tests.mesh.test_birdcage_leg_gaps import TERMINAL_AREA_BAND
    from tests.mesh.test_birdcage_port_sheets import PORT_LOWER, PORT_UPPER, SHEET_IFACE
    from tests.validation.test_ops61_exact_strip_ports import STRIP_BAND
    from tests.validation.test_wf12_f_human_hole import _census_checks, _owned_census

    comm = MPI.COMM_WORLD
    kw = _kwargs()
    ring_ports = _ring_ports()
    Lx, Ly, Lz = kw["phantom_box"]["size_m"]
    off = kw["phantom_box"]["table_offset_m"]
    L3 = Lx * Ly * Lz
    V = _rod_volume()
    A = _rod_area()
    centre = _rod_centre()
    _say(
        f"\n[GEO-36 s8] rod centre {centre} (Lx/4 toward +x), axis {ROD['axis']}; "
        f"V_rod {V:.12e}, A_rod {A:.12e}, L^3 {L3:.15e}"
    )

    armed = gmsh_stop_armed()
    _say(f"[GEO-36 s8] gmsh stop rule armed: {armed}")
    if armed:
        arm_gmsh_stop_rule(monkeypatch, comm, "GEO-36 s8")
    t0 = time.perf_counter()
    msh, cells, facets, diag = MeshGenerator.birdcage_port_domain(
        comm=comm, return_diagnostics=True, **kw
    )
    t_build = time.perf_counter() - t0
    _say(f"\n[GEO-36 s8] build done {t_build:.1f} s at -n {comm.size}")
    tdim = msh.topology.dim
    fdim = tdim - 1
    n_cells = int(msh.topology.index_map(tdim).size_global)
    _say(f"[GEO-36 s8] (G0) {n_cells} cells vs guard {F_CLINICAL_ROD_CELL_GUARD}")
    assert n_cells < F_CLINICAL_ROD_CELL_GUARD, n_cells  # (G0), before any facet reading

    msh.topology.create_entities(fdim)
    msh.topology.create_connectivity(fdim, tdim)
    msh.topology.create_entity_permutations()  # every rank, before any dS
    census = _owned_census(msh, cells, comm)

    # --- (G1), (G2): CAD.
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

    # --- (G3): owned tag-3 + tag-403 tets, reduced.
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
                    "ij,ij->i", X[:, 1] - X[:, 0], np.cross(X[:, 2] - X[:, 0], X[:, 3] - X[:, 0])
                )
            ).sum()
            / 6.0
        )

    X3 = _tets(PHANTOM_CELL_TAG)
    X403 = _tets(BIRDCAGE_IMPLANT_CELL_TAG)
    vol3 = comm.allreduce(_vol(X3), op=MPI.SUM)
    vol403 = comm.allreduce(_vol(X403), op=MPI.SUM)
    if len(X3):
        P = X3.reshape(-1, 3)
        lo_local, hi_local = P.min(axis=0), P.max(axis=0)
    else:
        lo_local, hi_local = np.full(3, np.inf), np.full(3, -np.inf)
    lo = np.array([comm.allreduce(float(x), op=MPI.MIN) for x in lo_local])
    hi = np.array([comm.allreduce(float(x), op=MPI.MAX) for x in hi_local])
    g3_vol = abs((vol3 + vol403) / L3 - 1.0)
    lo_ref = np.array([-0.5 * Lx, -0.5 * Ly + off, -0.5 * Lz])
    hi_ref = np.array([0.5 * Lx, 0.5 * Ly + off, 0.5 * Lz])
    g3_ext = float(
        max(
            np.max(np.abs(lo - lo_ref) / np.abs(lo_ref)),
            np.max(np.abs(hi - hi_ref) / np.abs(hi_ref)),
        )
    )

    # --- (G4), (N2r): the can partition with the rod.
    cad_sum = float(sum(cad.values()))
    m_cond = float(diag["conductor_cad_mass_m3"])
    m_air = float(diag["air_domain_cad_mass_m3"])
    g4 = abs((cad_sum + m_cond) / m_air - 1.0)
    n2_miss = abs((cad_sum + m_cond + rod_mass) / m_air - 1.0)
    n4_cond = abs(m_cond / C2_CONDUCTOR_CAD_MASS_M3 - 1.0)

    # --- (G5): the 402 interface (ghosted mesh, dS licensed).
    n402, n402_ext, n402_int, n402_pair, area402, moments = _interior_facet_check(
        msh, cells, facets, BIRDCAGE_IMPLANT_SURFACE_TAG, comm
    )
    g5_ratio = area402 / A

    # --- (G6) centroid.
    c = rod["wall_centroid_m"]
    g6 = max(abs(c[i] - centre[i]) for i in range(3))

    # --- (G7) census.
    census_fails = _census_checks(census, n_cells, ring_ports)

    # --- (G8): ring sheets, owned-facet numpy sums only.
    layout = diag["ring_port_layout"]
    chord = float(layout["ring_port_gap_chord_m"])
    width = float(layout["ring_port_box_width_m"])
    sheet_ref = chord * width
    n_gmsh = {int(k): int(v) for k, v in diag["ring_sheet_gmsh_triangles"].items()}
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

    # --- Unknowns (price input for a rod solve), memory.
    msh.topology.create_entities(1)
    msh.topology.create_entities(2)
    E = comm.allreduce(msh.topology.index_map(1).size_local, op=MPI.SUM)
    F = comm.allreduce(msh.topology.index_map(2).size_local, op=MPI.SUM)
    rss = comm.allreduce(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0**2, op=MPI.SUM)

    if comm.rank == 0:
        wa = max(g8_area, key=g8_area.get)
        g9_ints = (
            n_cells,
            n402,
            n402_ext,
            n402_int,
            n402_pair,
            int(census.get(PHANTOM_CELL_TAG, 0)),
            int(census.get(BIRDCAGE_IMPLANT_CELL_TAG, 0)),
            sum(n_owned.values()),
            E,
            F,
        )
        print(
            f"[GEO-36 s8] -n {comm.size}: {n_cells} cells (predicted 1 149 150 + [50 k, 100 k], "
            f"point 1 216 304; +{n_cells - 1149150} over C2); build {t_build:.1f} s (gmsh "
            f"{diag.get('mesh_wall_time_s', float('nan')):.2f} s; predicted T + [-40, +40] s, "
            f"T = 441.11 s); summed peak RSS {rss:.3f} GiB (predicted ≈ 2.7)\n"
            f"[GEO-36 s8] census {census}\n"
            f"[GEO-36 s8] (G1) rod_mass {rod_mass:.12e} / kept {kept:.12e} vs V_rod {V:.12e}: "
            f"rel {g1_mass:.3e} / {g1_kept:.3e}; OCC wall {float(rod['wall_area_m2']):.12e} vs "
            f"A_rod {A:.12e}: rel {g1_area:.3e} (band {STRIP_BAND:.0e})\n"
            f"[GEO-36 s8] (G2) CAD(403)/V {g2_403:.3e}; CAD(3) {cad['phantom']:.15e} vs L^3 - V "
            f"{L3 - V:.15e}: rel {g2_3:.3e}; CAD(3)+CAD(403) vs L^3 {L3:.15e}: rel {g2_sum:.3e}\n"
            f"[GEO-36 s8] (N1r) CAD(3) vs C2's L^3 misses by {n1_miss:.6e} (predicted "
            f"{n1_pred:.6e}; {n1_miss / STRIP_BAND:.3e}x STRIP_BAND)\n"
            f"[GEO-36 s8] (G3) meshed tag-3 {vol3:.15e} + tag-403 {vol403:.15e} vs L^3: rel "
            f"{g3_vol:.3e}; tag-3 extents lo {lo.tolist()} hi {hi.tolist()}: worst rel "
            f"{g3_ext:.3e}\n"
            f"[GEO-36 s8] (G4) CAD groups {cad_sum:.15e} + conductor {m_cond:.15e} vs can "
            f"{m_air:.15e}: rel {g4:.3e} (groups {cad})\n"
            f"[GEO-36 s8] (N2r) pre-repair reading (+ V_rod) misses by {n2_miss:.6e} (predicted "
            f"2.366e-07; {n2_miss / STRIP_BAND:.3e}x STRIP_BAND)\n"
            f"[GEO-36 s8] (G5) 402 owned {n402}, exterior {n402_ext}, interior {n402_int}, "
            f"3|403 pairs {n402_pair}; meshed dS(402)/A_rod {g5_ratio:.6f} (band "
            f"{TERMINAL_AREA_BAND})\n"
            f"[GEO-36 s8] (G6) OCC wall centroid {c} vs {centre}: max dev {g6:.3e} "
            f"(bound {EXACT * ROD_LENGTH:.1e}); meshed "
            f"{tuple(m / area402 for m in moments) if area402 else None}\n"
            f"[GEO-36 s8] (G7) census failures {census_fails or 'none'}; implant_spheroid_cad "
            f"{diag.get('implant_spheroid_cad')}\n"
            f"[GEO-36 s8] (G8) owned sheet facets == N_gmsh on {len(ring_ports) - len(g8_bad)}"
            f"/{len(ring_ports)}; mismatches {g8_bad or 'none'}; N_gmsh "
            f"{min(n_gmsh.values())}..{max(n_gmsh.values())}, total {sum(n_gmsh.values())}; "
            f"worst area P{wa} {g8_area[wa]:.3e}\n"
            f"[GEO-36 s8] (G9) ints {g9_ints}; CAD floats rod {rod_mass!r} wall "
            f"{float(rod['wall_area_m2'])!r} phantom {cad['phantom']!r} implant "
            f"{cad['implant']!r} conductor {m_cond!r} can {m_air!r}\n"
            f"[GEO-36 s8] (N4r) tag-403 cells {census.get(BIRDCAGE_IMPLANT_CELL_TAG, 0)} "
            f"(step 6: 15 025); meshed 403 / V {vol403 / V:.6f} (step 6: 0.986763); meshed "
            f"402 / A {g5_ratio:.6f} (step 6: 0.996276); conductor vs C2 "
            f"{C2_CONDUCTOR_CAD_MASS_M3:.15e}: rel {n4_cond:.3e} (predicted <= 1e-9; printed "
            f"only)\n"
            f"[GEO-36 s8] unknowns E {E}, F {F}; U1 = E = {E}, U2 = 2E + 2F = {2 * E + 2 * F}",
            flush=True,
        )

    # Negative controls (N1r), (N2r).
    assert n1_miss > STRIP_BAND, n1_miss
    assert n2_miss > STRIP_BAND, n2_miss
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
    assert g6 <= EXACT * ROD_LENGTH, c
    # (G7)
    assert not census_fails, census_fails
    assert PHANTOM_CELL_TAG in census and BIRDCAGE_IMPLANT_CELL_TAG in census, census
    assert 1 not in census, census
    assert diag.get("implant_spheroid_cad") is None
    # (G8)
    assert len(n_gmsh) == 2 * LEG_COUNT, sorted(n_gmsh)
    assert not g8_bad, g8_bad
    for i, dev in g8_area.items():
        assert dev <= STRIP_BAND, (i, dev)


# --- §9 item 560 (§10 I35): the column knob on step 8's module. -------------
# Brief: docs/planning/chunks/GEO-36.md, "§9 item 560", and the 2026-10-20
# 03:00 review ruling beneath it (rod σ = 800 S/m with the gel's ε_r, IMP-3
# step 8's precedent; "(D-rod)" read as I33's (G1)/(G2)/(G4)/(G7)).
#
# FEM_EM_GEO36_S8_COLUMN: unset/empty ⇒ skip; ``build`` ⇒ build-only proof
# (this module's W1); ``1`` ⇒ entry 45 only: entry 40's column body run
# unchanged on the rod fixture through `_route_column` (four monkeypatches).
# The gmsh stop rule is never armed here: a priced rebuild (I33).

from tests.complex_mode import complex_only  # noqa: E402

S8_COLUMN_ENV = "FEM_EM_GEO36_S8_COLUMN"
#: I33's measured cells (`20261019T201444Z_GEO-36-rodbox-W2-n2.log`; equal at
#: -n 4 in `20261019T200702Z_GEO-36-rodbox-W1-n4.log`), registered as the
#: column's (A0) comparand under the (1*) licence, not a new target.
F_CLINICAL_ROD_CELL_RECORD = 1_217_025
#: I33's (U1 = E, U2 = 2E + 2F), same log: E 1 606 744, F 2 556 953.
F_CLINICAL_ROD_UNKNOWNS_RECORD = (1_606_744, 8_327_394)
#: entry 40's 1.5 T r1 record, read back in (K1) / (N-A0).
_R1_CELL_RECORD_EXPECTED = 1_147_366
_STEP1B_ENV = "FEM_EM_GEO38_STEP1B"


def _s8_column_mode():
    raw = os.environ.get(S8_COLUMN_ENV, "").strip()
    if raw == "":
        return None
    if raw in ("build", "1"):
        return raw
    raise ValueError(f"{S8_COLUMN_ENV} must be unset, 'build' or '1', got {raw!r}")


def _rod_material_minmax(msh, cell_tags, problem, comm, tag):
    """Global (min, max) of σ and ε_r over owned cells of ``tag`` (IMP-3 (M0))."""
    from fem_em_solver.core.time_harmonic import build_material_fields
    from tests.mesh.test_birdcage_implant_rod import _cell_values_with_ghosts

    n_owned = msh.topology.index_map(msh.topology.dim).size_local
    tagv = _cell_values_with_ghosts(msh, cell_tags)
    owned = np.flatnonzero(tagv[:n_owned] == tag)
    sigma_f, eps_f = build_material_fields(
        msh,
        problem.material,
        cell_tags=problem.cell_tags,
        material_map=problem.material_map,
        phantom_material=problem.phantom_material,
        phantom_tag=problem.phantom_tag,
    )
    res = {"n": comm.allreduce(int(len(owned)), op=MPI.SUM)}
    for name, f in (("sigma", sigma_f), ("eps", eps_f)):
        dofs = f.function_space.dofmap.list[:, 0]
        vals = np.real(f.x.array[dofs[owned]])
        lo = float(vals.min()) if vals.size else np.inf
        hi = float(vals.max()) if vals.size else -np.inf
        res[name] = (comm.allreduce(lo, op=MPI.MIN), comm.allreduce(hi, op=MPI.MAX))
    return res


def _build_rod_f_clinical(comm):
    """`col._build_f_clinical`'s dict shape on the step-8 rod fixture, with the
    rod's material wired and (D-rod) / (M) / (N-M) / (C) / (U1) asserted."""
    import dataclasses

    from fem_em_solver.core import HomogeneousMaterial, TimeHarmonicSolver
    from fem_em_solver.io.mesh import BIRDCAGE_IMPLANT_CELL_TAG, MeshGenerator
    from tests.validation import test_geo38_f_clinical_column as col
    from tests.validation import test_geo38_f_clinical_matrix as c3
    from tests.validation.test_lossy_sphere_fullwave import SALINE_SIGMA
    from tests.validation.test_ops61_exact_strip_ports import STRIP_BAND
    from tests.validation.test_port_gap_voltage_impedance import SIGMA_WIRE_S_PER_M
    from tests.validation.test_wf7_f_human_b1_quadrature import DEGREE
    from tests.validation.test_wf12_f_human_hole import _census_checks, _owned_census

    kw = _kwargs()
    t0 = time.perf_counter()
    msh, cell_tags, facets, diag = MeshGenerator.birdcage_port_domain(
        comm=comm, return_diagnostics=True, **kw
    )
    t_build = time.perf_counter() - t0
    n_cells = int(msh.topology.index_map(3).size_global)
    _say(
        f"\n[GEO-36 I35] build {t_build:.1f} s at -n {comm.size} (gmsh "
        f"{diag.get('mesh_wall_time_s', float('nan')):.2f} s; PREDICTED ≈ 445 s [410, 560]); "
        f"cells {n_cells}"
    )
    t1 = time.perf_counter()
    chord = float(diag["ring_port_layout"]["ring_port_gap_chord_m"])
    ctx, outer_census, _w = col._context_from_mesh(msh, cell_tags, facets, chord, comm)

    # --- the rod's material (ruling: σ_wire with the gel's ε_r, μ_r); one solver.
    problem_gel_only = ctx["solver"].problem
    mmap = dict(problem_gel_only.material_map)
    gel = mmap[PHANTOM_CELL_TAG]
    problem_rod = dataclasses.replace(
        problem_gel_only,
        material_map={
            **mmap,
            BIRDCAGE_IMPLANT_CELL_TAG: HomogeneousMaterial(
                sigma=SIGMA_WIRE_S_PER_M, epsilon_r=gel.epsilon_r, mu_r=gel.mu_r
            ),
        },
    )
    del ctx["solver"]
    ctx["solver"] = TimeHarmonicSolver(problem_rod, degree=DEGREE)

    # --- (D-rod): I33's (G1) / (G2) / (G4) / (G7) at their imported bands.
    Lx, Ly, Lz = kw["phantom_box"]["size_m"]
    L3 = Lx * Ly * Lz
    V = _rod_volume()
    A = _rod_area()
    rod = diag["implant_rod_cad"]
    cad = dict(diag["cad_mass_by_group"])
    g1 = max(
        abs(float(rod["rod_mass_m3"]) / V - 1.0),
        abs(float(rod["kept_mass_m3"]) / V - 1.0),
        abs(float(rod["wall_area_m2"]) / A - 1.0),
    )
    g2 = max(
        abs(cad["implant"] / V - 1.0),
        abs(cad["phantom"] / (L3 - V) - 1.0),
        abs((cad["phantom"] + cad["implant"]) / L3 - 1.0),
    )
    g4 = abs(
        (float(sum(cad.values())) + float(diag["conductor_cad_mass_m3"]))
        / float(diag["air_domain_cad_mass_m3"])
        - 1.0
    )
    census = _owned_census(msh, cell_tags, comm)
    g7_fails = _census_checks(census, n_cells, _ring_ports())

    # --- (M) / (N-M): owned tag-403 and tag-3 cells, reduced.
    m_rod = _rod_material_minmax(msh, cell_tags, problem_rod, comm, BIRDCAGE_IMPLANT_CELL_TAG)
    m_gel = _rod_material_minmax(msh, cell_tags, problem_rod, comm, PHANTOM_CELL_TAG)
    n_m = _rod_material_minmax(msh, cell_tags, problem_gel_only, comm, BIRDCAGE_IMPLANT_CELL_TAG)

    # --- (C): 64 specs P33..P96, each with owned facets.
    geo = c3._sheet_geometry(ctx, comm)
    ids = [sp.port_id for sp in ctx["specs"]]
    c_fails = c3._census_failures(ids, {pid: g[0] for pid, g in geo.items()})

    # --- (U1): degree-1 global dofs = E.
    V_h = ctx["solver"].function_space()
    u1 = int(V_h.dofmap.index_map.size_global * V_h.dofmap.index_map_bs)
    t_ctx = time.perf_counter() - t1

    _say(
        f"[GEO-36 I35] (D-rod) (G1) {g1:.3e}, (G2) {g2:.3e}, (G4) {g4:.3e} (band "
        f"{STRIP_BAND:.0e}, ASSERTED; PREDICTED ≤ 4.4e-16 / ≤ 4.4e-16 / 6.457e-11); (G7) "
        f"census failures {g7_fails or 'none'}, census {census}, implant_spheroid_cad "
        f"{diag.get('implant_spheroid_cad')}\n"
        f"[GEO-36 I35] (M) tag 403 ({m_rod['n']} cells): sigma {m_rod['sigma']} eps_r "
        f"{m_rod['eps']} vs ({SIGMA_WIRE_S_PER_M}, gel eps_r {gel.epsilon_r}); tag 3 "
        f"({m_gel['n']} cells): sigma {m_gel['sigma']} vs SALINE_SIGMA {SALINE_SIGMA} (ASSERTED "
        f"exact)\n"
        f"[GEO-36 I35] (N-M) unpatched problem, tag 403: sigma {n_m['sigma']} (air default "
        f"0.0, ASSERTED)\n"
        f"[GEO-36 I35] (C) {len(ids)} specs, owned facets "
        f"{min(g[0] for g in geo.values())}..{max(g[0] for g in geo.values())}, failures "
        f"{c_fails or 'none'}; (U1) dofs {u1} vs E {F_CLINICAL_ROD_UNKNOWNS_RECORD[0]} "
        f"(ASSERTED exact); context + checks {t_ctx:.1f} s (PREDICTED ≈ 10–20 s)"
    )
    assert g1 < STRIP_BAND, g1
    assert g2 < STRIP_BAND, g2
    assert g4 <= STRIP_BAND, g4
    assert not g7_fails, g7_fails
    assert PHANTOM_CELL_TAG in census and BIRDCAGE_IMPLANT_CELL_TAG in census, census
    assert 1 not in census, census
    assert diag.get("implant_spheroid_cad") is None
    assert m_rod["n"] > 0 and m_gel["n"] > 0, (m_rod, m_gel)
    assert m_rod["sigma"] == (SIGMA_WIRE_S_PER_M, SIGMA_WIRE_S_PER_M), m_rod
    assert m_rod["eps"] == (gel.epsilon_r, gel.epsilon_r), m_rod
    assert m_gel["sigma"] == (SALINE_SIGMA, SALINE_SIGMA), m_gel
    assert n_m["sigma"] == (0.0, 0.0), n_m
    assert not c_fails, c_fails
    assert u1 == F_CLINICAL_ROD_UNKNOWNS_RECORD[0], u1
    return {
        "ctx": ctx,
        "cells": n_cells,
        "facets": facets,
        "shield_radius_m": float(diag["shield_radius_m"]),
        "outer_census": outer_census,
        "problem_gel_only": problem_gel_only,
        "t_build": t_build,
    }


def _route_column(monkeypatch):
    """Entry 40's column body, unchanged, on the rod fixture (C8's no-edit route)."""
    from tests.validation import test_geo38_f_clinical_build as s1a
    from tests.validation import test_geo38_f_clinical_column as col

    monkeypatch.setenv(_STEP1B_ENV, "1")
    monkeypatch.setattr(col, "_build_f_clinical", _build_rod_f_clinical)
    monkeypatch.setattr(s1a, "F_CLINICAL_R1_CELL_RECORD", F_CLINICAL_ROD_CELL_RECORD)
    monkeypatch.setattr(s1a, "F_CLINICAL_R1_UNKNOWNS_RECORD", F_CLINICAL_ROD_UNKNOWNS_RECORD)


@complex_only
def test_geo36_i35_rod_column(monkeypatch):
    mode = _s8_column_mode()
    if mode is None:
        pytest.skip(f"{S8_COLUMN_ENV} unset: GEO-36 I35 is selected by environment only")
    from tests.validation import test_geo38_f_clinical_column as col

    if mode == "1":
        _route_column(monkeypatch)
        col.test_geo38_step1b_f_clinical_copper_column()
        return

    # mode == "build": the column body's pre-column asserts, then return.
    from fem_em_solver.io.mesh import BIRDCAGE_CONDUCTOR_SURFACE_TAG
    from tests.mesh.test_birdcage_port_sheet_prerequisite import CELL_COUNT_BAND
    from tests.validation.test_wf12_f_human_hole import _census_checks, _owned_census, _rss_gib

    comm = MPI.COMM_WORLD
    t_window = time.perf_counter()
    built = _build_rod_f_clinical(comm)
    ctx = built["ctx"]
    msh, cell_tags = ctx["msh"], ctx["cell_tags"]
    facets = built["facets"]
    n_cells = built["cells"]
    census = _owned_census(msh, cell_tags, comm)
    census_fails = _census_checks(census, n_cells, _ring_ports())
    oc = built["outer_census"]
    cell_rel = n_cells / F_CLINICAL_ROD_CELL_RECORD - 1.0
    R_s = float(built["shield_radius_m"])
    fv = np.asarray(facets.values)
    fi = np.asarray(facets.indices)
    cav_p = _facet_coords(msh, _owned(msh, fi[fv == BIRDCAGE_CONDUCTOR_SURFACE_TAG])).reshape(-1, 3)
    cav_hmax = comm.allreduce(
        float(np.hypot(cav_p[:, 0], cav_p[:, 1]).max()) if len(cav_p) else 0.0, op=MPI.MAX
    )
    _mine, per, tot = _rss_gib(comm)
    _say(
        f"[GEO-36 I35] (A0) cells {n_cells} vs record {F_CLINICAL_ROD_CELL_RECORD}, rel "
        f"{cell_rel:.3e} (band {CELL_COUNT_BAND}, ASSERTED; PREDICTED 0); (iv) census failures "
        f"{census_fails or 'none'}; outer census {oc}; (S) cavity max hypot {cav_hmax:.9f} m < "
        f"R_s {R_s} (ASSERTED)"
    )
    _say(
        f"[GEO-36 I35] BUILD-ONLY return: t_build {built['t_build']:.1f} s; RSS per rank "
        f"{', '.join(f'{r:.2f}' for r in (per or []))} GiB, summed {tot:.3f} GiB (PREDICTED ≈ "
        f"2.2); window {time.perf_counter() - t_window:.1f} s (PREDICTED ≈ 465 s)"
    )
    _say(
        f"[GEO-36 I35] for entry 45 (printed): rho = {n_cells}/369 387 = {n_cells / 369_387:.4f}; "
        f"conservative column 39.2*(E/500 568)^2 = "
        f"{39.2 * (F_CLINICAL_ROD_UNKNOWNS_RECORD[0] / 500_568) ** 2:.0f} s (entry 40: 364 s)"
    )
    assert abs(cell_rel) <= CELL_COUNT_BAND, cell_rel
    assert not census_fails, census_fails
    assert oc["tag_collision"] == 0 and oc["cavity"] > 0, oc
    assert (
        oc["cavity_exterior"] == oc["cavity"] and oc["outer"] + oc["cavity"] == oc["exterior"]
    ), oc
    assert cav_hmax < R_s, (cav_hmax, R_s)


# --- I35 coding tests (no gmsh); skip whenever the knob is set. -------------
_s8_coding = pytest.mark.skipif(
    os.environ.get(S8_COLUMN_ENV, "").strip() != "",
    reason=f"{S8_COLUMN_ENV} set: coding tests skip in build windows",
)


class _Routed(Exception):
    def __init__(self, kwargs):
        super().__init__("routed")
        self.kwargs = kwargs


@_s8_coding
def test_i35_k0_knob_parsing(monkeypatch):
    monkeypatch.delenv(S8_COLUMN_ENV, raising=False)
    assert _s8_column_mode() is None
    for ok in ("build", "1"):
        monkeypatch.setenv(S8_COLUMN_ENV, ok)
        assert _s8_column_mode() == ok
    for bad in ("2", "x"):
        monkeypatch.setenv(S8_COLUMN_ENV, bad)
        with pytest.raises(ValueError):
            _s8_column_mode()


def _call_column_into_recorder(monkeypatch):
    from fem_em_solver.io.mesh import MeshGenerator
    from tests.validation import test_geo38_f_clinical_column as col

    def _rec(*args, **kwargs):
        raise _Routed(kwargs)

    monkeypatch.setattr(MeshGenerator, "birdcage_port_domain", staticmethod(_rec))
    with pytest.raises(_Routed) as exc:
        col.test_geo38_step1b_f_clinical_copper_column()
    return exc.value.kwargs


@_s8_coding
def test_i35_k1_routing_tripwire(monkeypatch):
    from tests.validation import test_geo38_f_clinical_build as s1a
    from tests.validation import test_geo38_f_clinical_column as col

    assert os.environ.get(_STEP1B_ENV) is None
    original = col._build_f_clinical
    with monkeypatch.context() as m:
        _route_column(m)
        kw = _call_column_into_recorder(m)
        assert s1a.F_CLINICAL_R1_CELL_RECORD == F_CLINICAL_ROD_CELL_RECORD
    rod = kw.get("implant_rod")
    _say(
        f"\n[GEO-36 I35] (K1) routed kwargs: implant_rod {rod}; phantom_box "
        f"{kw.get('phantom_box')}"
    )
    assert rod is not None and rod["keep_cells"] is True, rod
    assert tuple(rod["centre"]) == _rod_centre(), rod
    assert kw.get("phantom_box") is not None, sorted(kw)
    assert col._build_f_clinical is original
    assert s1a.F_CLINICAL_R1_CELL_RECORD == _R1_CELL_RECORD_EXPECTED
    assert os.environ.get(_STEP1B_ENV) is None


@_s8_coding
def test_i35_n_route_unpatched_reaches_15t_fixture(monkeypatch):
    """(N-route): STEP1B=1 alone reaches the 1.5 T r1 kwargs, with no rod."""
    assert os.environ.get(_STEP1B_ENV) is None
    with monkeypatch.context() as m:
        m.setenv(_STEP1B_ENV, "1")
        kw = _call_column_into_recorder(m)
    _say(f"\n[GEO-36 I35] (N-route) unpatched kwargs implant_rod: {kw.get('implant_rod')}")
    assert kw.get("implant_rod") is None, kw.get("implant_rod")
    assert os.environ.get(_STEP1B_ENV) is None


@_s8_coding
def test_i35_n_a0_15t_record_misses_rod_mesh():
    """(N-A0): the 1.5 T record misses the rod mesh by a fixed ratio of records."""
    from tests.mesh.test_birdcage_port_sheet_prerequisite import CELL_COUNT_BAND

    miss = abs(F_CLINICAL_ROD_CELL_RECORD / _R1_CELL_RECORD_EXPECTED - 1.0)
    _say(
        f"\n[GEO-36 I35] (N-A0) |{F_CLINICAL_ROD_CELL_RECORD}/{_R1_CELL_RECORD_EXPECTED} - 1| = "
        f"{miss:.3e} = {miss / CELL_COUNT_BAND:.2f}x CELL_COUNT_BAND (PREDICTED 6.07e-2, 6.1x)"
    )
    assert miss > CELL_COUNT_BAND, miss
