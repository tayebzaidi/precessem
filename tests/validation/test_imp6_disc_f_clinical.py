"""`IMP-6` step 4 (§9 item 482, §10 I31): the kept aspect-5 disc in F-clinical's
ASTM box — CAD and mesh identities only, no solve.

`_box_kwargs()` (item 470's r1 kwargs with the ASTM box, `GEO-38` step 2) plus
``implant_spheroid = dict(SPHEROID, keep_cells=True)`` (`GEO-36` step 7's body
and placement: a = 6 mm, c = 1.2 mm, centre (0.015, 0, 0), axis x̂ — a
modelling choice, no ASTM-site claim). Brief: docs/planning/chunks/IMP-6.md,
"§9 item 482".

The build test is selected by environment only (``FEM_EM_IMP6_STEP4=1``), never
by ``-k``; it skips when unset. The refusal tests are pure and always run.

Anchors (asserted at the window's width; owned cells / facets, reduced):
  (D0) cells < F_CLINICAL_DISC_CELL_GUARD before any facet reading;
  (D1) body mass = 4πa²c/3 and OCC wall area = A at CAD_ANALYTIC_BAND;
  (D2) CAD(403) = kept = body, CAD(3) = L³ − body, CAD(3)+CAD(403) = L³ at
       STRIP_BAND;
  (D3) owned tag-3 + tag-403 tet volumes = L³, tag-3 extents = box at STRIP_BAND;
  (D4) Σ cad_mass_by_group + conductor_cad_mass_m3 = can at STRIP_BAND;
  (D5) 402 interior, between tag 3 and 403, meshed dS(402)/A in TERMINAL_AREA_BAND;
  (D6) OCC wall centroid at SPHEROID_CENTRE to EXACT·a;
  (D7) rim fidelity (n_mask > 0, n_inside 0, n403 > 0, n403_outside 0);
  (D8) census clean, tags 3 and 403 present, no tag 1, implant_rod_cad None;
  (D9) per ring sheet owned facets = N_gmsh, area = chord·w at STRIP_BAND;
  (D10) refusals (pure).
Negative controls: (N1) C2's comparand L³ misses CAD(3) by body/(L³ − body),
asserted > STRIP_BAND, factor predicted 7.36e3×; (N2) item 470's pre-repair
conductor reading misses (D4) by V/can, asserted > STRIP_BAND, factor
predicted 435×; (N3) a 1 mm poke-out accepted by the cylinder test, refused
by the box test (both asserted, W0); (N4) tag-403 count, meshed 403 / V and
the conductor mass against C2's are printed only.
"""

from __future__ import annotations

import math
import os
import resource
import time

import numpy as np
import pytest
from mpi4py import MPI

from tests.mesh.test_birdcage_implant_spheroid import (
    SPHEROID,
    SPHEROID_CENTRE,
    SPHEROID_SEMI_MAJOR,
    _rim_census,
    _spheroid_area,
    _spheroid_volume,
)
from tests.validation.gmsh_stop_rule import arm_gmsh_stop_rule, gmsh_stop_armed
from tests.validation.test_geo38_f_clinical_box import _box_kwargs
from tests.validation.test_geo38_f_clinical_build import (
    LEG_COUNT,
    _area,
    _facet_coords,
    _owned,
    _ring_ports,
)

STEP4_ENV = "FEM_EM_IMP6_STEP4"

#: Guard, checked before any facet reading: predicted 1 149 150 + [35 k, 100 k]
#: plus ≈ 17 % headroom (F-human's 600 000 sat over 417 653).
F_CLINICAL_DISC_CELL_GUARD = 1_400_000
#: C2's conductor CAD mass, a printed comparand only
#: (`20261015T185056Z_GEO-38-s2-W3-n4.log:11903`).
C2_CONDUCTOR_CAD_MASS_M3 = 1.994040484843362e-03

PHANTOM_CELL_TAG = 3

#: (N3): the planted poke-out — |c_y| + a = 0.046 > Ly/2 = 0.045.
POKE_OUT_CENTRE = (0.015, 0.040, 0.0)


def _say(msg):
    if MPI.COMM_WORLD.rank == 0:
        print(msg, flush=True)


def _kwargs():
    kw = _box_kwargs()
    kw["implant_spheroid"] = dict(SPHEROID, keep_cells=True)
    return kw


def _refuse_before_gmsh(monkeypatch, **kw):
    """`birdcage_port_domain` with the builder replaced by a tripwire."""
    from fem_em_solver.io.mesh import MeshGenerator

    def _tripwire(*a, **k):
        raise AssertionError("reached the gmsh builder: the refusal did not fire")

    monkeypatch.setattr(MeshGenerator, "_build_birdcage_port_model", staticmethod(_tripwire))
    return MeshGenerator.birdcage_port_domain(comm=MPI.COMM_WORLD, return_diagnostics=True, **kw)


# --- (D10) / (N3): pure refusals, always run. -------------------------------


def test_n3_poke_out_cylinder_accepts_box_refuses(monkeypatch):
    from fem_em_solver.io.mesh import MeshGenerator

    box = _box_kwargs()["phantom_box"]
    Lx, Ly, Lz = box["size_m"]
    off = box["table_offset_m"]
    R = float(np.hypot(0.5 * Lx, 0.5 * Ly + abs(off)))
    req = dict(SPHEROID, centre=POKE_OUT_CENTRE, keep_cells=True)
    # Cylinder test (phantom_box=None): accepted, margin ≈ 0.166 m.
    out = MeshGenerator._normalise_implant_spheroid(
        dict(req), phantom_radius=R, phantom_height=Lz, as_hole=True
    )
    margin = R - (float(np.hypot(*POKE_OUT_CENTRE[:2])) + SPHEROID_SEMI_MAJOR)
    breach = abs(POKE_OUT_CENTRE[1] - off) + SPHEROID_SEMI_MAJOR - 0.5 * Ly
    _say(f"\n[IMP-6 s4] (N3) R {R:.6f} m, H {Lz} m: cylinder test accepts with margin "
         f"{margin:.6f} m; box breach |c_y - off| + a - Ly/2 = {breach:.6e} m")
    assert out["keep_cells"] is True
    assert margin > 0.16
    assert breach > 0.0
    # Box test: refused before gmsh.
    kw = _box_kwargs()
    kw["implant_spheroid"] = dict(req)
    with pytest.raises(ValueError, match="not strictly inside the phantom_box") as exc:
        _refuse_before_gmsh(monkeypatch, **kw)
    _say(f"[IMP-6 s4] (N3) box refusal: {exc.value}")


def test_d10_rod_and_disc_with_box_refused(monkeypatch):
    # `GEO-36` step 8 (§9 item 530, review ruling 6): the rod now composes
    # with the box, so this re-points to the one-implant-per-build refusal.
    from tests.mesh.test_birdcage_implant_rod import ROD

    kw = _kwargs()
    kw["implant_rod"] = dict(ROD, keep_cells=True)
    with pytest.raises(ValueError, match="implant_spheroid cannot be combined with implant_rod"):
        _refuse_before_gmsh(monkeypatch, **kw)


def test_d10_hole_spheroid_under_as_hole_still_refused(monkeypatch):
    kw = _box_kwargs()
    assert kw.get("as_hole") is True
    kw["implant_spheroid"] = dict(SPHEROID, keep_cells=False)
    with pytest.raises(ValueError, match="under as_hole needs 'keep_cells': True"):
        _refuse_before_gmsh(monkeypatch, **kw)


def test_d10_kept_with_box_accepted():
    from fem_em_solver.io.mesh import MeshGenerator

    box = _box_kwargs()["phantom_box"]
    Lx, Ly, Lz = box["size_m"]
    R = float(np.hypot(0.5 * Lx, 0.5 * Ly + abs(box["table_offset_m"])))
    out = MeshGenerator._normalise_implant_spheroid(
        dict(SPHEROID, keep_cells=True), phantom_radius=R, phantom_height=Lz,
        as_hole=True, phantom_box=dict(box),
    )
    assert out["keep_cells"] is True
    assert out["centre"] == tuple(SPHEROID_CENTRE)


# --- The build (environment-selected). ---------------------------------------


def test_imp6_step4_disc_in_astm_box(monkeypatch):
    if os.environ.get(STEP4_ENV, "").strip() != "1":
        pytest.skip(f"{STEP4_ENV} unset: IMP-6 step 4 is selected by environment only")
    import dolfinx

    from fem_em_solver.io.mesh import (
        BIRDCAGE_IMPLANT_CELL_TAG,
        BIRDCAGE_IMPLANT_SURFACE_TAG,
        MeshGenerator,
        _interface_facet_tags,
    )
    from tests.mesh.test_birdcage_conductor_hole import CAD_ANALYTIC_BAND, EXACT
    from tests.mesh.test_birdcage_implant_rod import _interior_facet_check
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
    V = _spheroid_volume()
    A = _spheroid_area()

    armed = gmsh_stop_armed()
    _say(f"[IMP-6 s4] gmsh stop rule armed: {armed}")
    if armed:
        arm_gmsh_stop_rule(monkeypatch, comm, "IMP-6 s4")
    t0 = time.perf_counter()
    msh, cells, facets, diag = MeshGenerator.birdcage_port_domain(
        comm=comm, return_diagnostics=True, **kw)
    t_build = time.perf_counter() - t0
    _say(f"\n[IMP-6 s4] build done {t_build:.1f} s at -n {comm.size}")
    tdim = msh.topology.dim
    fdim = tdim - 1
    n_cells = int(msh.topology.index_map(tdim).size_global)
    _say(f"[IMP-6 s4] (D0) {n_cells} cells vs guard {F_CLINICAL_DISC_CELL_GUARD}")
    assert n_cells < F_CLINICAL_DISC_CELL_GUARD, n_cells  # (D0), before any facet reading

    msh.topology.create_entities(fdim)
    msh.topology.create_connectivity(fdim, tdim)
    msh.topology.create_entity_permutations()  # every rank, before any dS (3b-iv)
    census = _owned_census(msh, cells, comm)

    # --- (D1), (D2): CAD.
    sph = diag["implant_spheroid_cad"]
    cad = dict(diag["cad_mass_by_group"])
    body = float(sph["body_mass_m3"])
    d1_mass = abs(body / V - 1.0)
    d1_area = abs(float(sph["wall_area_m2"]) / A - 1.0)
    d2_403 = abs(cad["implant"] / body - 1.0)
    d2_kept = abs(float(sph["kept_mass_m3"]) / body - 1.0)
    d2_3 = abs(cad["phantom"] / (L3 - body) - 1.0)
    d2_sum = abs((cad["phantom"] + cad["implant"]) / L3 - 1.0)
    n1_miss = abs(cad["phantom"] / L3 - 1.0)
    n1_pred = body / (L3 - body)

    # --- (D3): owned tag-3 + tag-403 tets, reduced.
    cv = np.asarray(cells.values)
    ci = np.asarray(cells.indices)
    n_owned_cells = msh.topology.index_map(tdim).size_local

    def _tets(tag):
        sel = ci[(cv == tag) & (ci < n_owned_cells)]
        if len(sel) == 0:
            return np.zeros((0, 4, 3))
        g = dolfinx.cpp.mesh.entities_to_geometry(
            msh._cpp_object, tdim, np.asarray(sel, dtype=np.int32), False)
        return msh.geometry.x[g[:, :4]]

    def _vol(X):
        if len(X) == 0:
            return 0.0
        return float(np.abs(np.einsum(
            "ij,ij->i", X[:, 1] - X[:, 0],
            np.cross(X[:, 2] - X[:, 0], X[:, 3] - X[:, 0]))).sum() / 6.0)

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
    d3_vol = abs((vol3 + vol403) / L3 - 1.0)
    lo_ref = np.array([-0.5 * Lx, -0.5 * Ly + off, -0.5 * Lz])
    hi_ref = np.array([0.5 * Lx, 0.5 * Ly + off, 0.5 * Lz])
    d3_ext = float(max(np.max(np.abs(lo - lo_ref) / np.abs(lo_ref)),
                       np.max(np.abs(hi - hi_ref) / np.abs(hi_ref))))

    # --- (D4), (N2): the can partition with the body.
    cad_sum = float(sum(cad.values()))
    m_cond = float(diag["conductor_cad_mass_m3"])
    m_air = float(diag["air_domain_cad_mass_m3"])
    d4 = abs((cad_sum + m_cond) / m_air - 1.0)
    n2_miss = abs((cad_sum + m_cond + body) / m_air - 1.0)
    n4_cond = abs(m_cond / C2_CONDUCTOR_CAD_MASS_M3 - 1.0)

    # --- (D5): the 402 interface (ghosted mesh, dS licensed).
    n402, n402_ext, n402_int, n402_pair, area402, moments = _interior_facet_check(
        msh, cells, facets, BIRDCAGE_IMPLANT_SURFACE_TAG, comm)
    d5_ratio = area402 / A

    # --- (D6) centroid.
    c = sph["wall_centroid_m"]
    d6 = max(abs(c[i] - SPHEROID_CENTRE[i]) for i in range(3))

    # --- (D7) rim fidelity.
    n_mask, n_in, s_min, n403_mid, n_out, s_max = _rim_census(msh, cells, comm)

    # --- (D8) census.
    census_fails = _census_checks(census, n_cells, ring_ports)

    # --- (D9): ring sheets, owned-facet numpy sums only.
    layout = diag["ring_port_layout"]
    chord = float(layout["ring_port_gap_chord_m"])
    width = float(layout["ring_port_box_width_m"])
    sheet_ref = chord * width
    n_gmsh = {int(k): int(v) for k, v in diag["ring_sheet_gmsh_triangles"].items()}
    sheet_tags = _interface_facet_tags(
        msh, cells, {SHEET_IFACE + i: (PORT_LOWER + i, PORT_UPPER + i) for i in ring_ports})
    sv = np.asarray(sheet_tags.values)
    si = np.asarray(sheet_tags.indices)
    local = {i: _facet_coords(msh, _owned(msh, si[sv == SHEET_IFACE + i])) for i in ring_ports}
    gathered = comm.allgather(local)
    sheet = {i: np.concatenate([g[i] for g in gathered]) for i in ring_ports}
    n_owned = {i: len(sheet[i]) for i in ring_ports}
    d9_bad = {i: (n_owned[i], n_gmsh.get(i)) for i in ring_ports if n_owned[i] != n_gmsh.get(i)}
    d9_area = {i: abs(float(_area(sheet[i]).sum()) / sheet_ref - 1.0) for i in ring_ports}

    # --- Unknowns (price input for the solve step), memory.
    msh.topology.create_entities(1)
    msh.topology.create_entities(2)
    E = comm.allreduce(msh.topology.index_map(1).size_local, op=MPI.SUM)
    F = comm.allreduce(msh.topology.index_map(2).size_local, op=MPI.SUM)
    rss = comm.allreduce(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0**2,
                         op=MPI.SUM)

    if comm.rank == 0:
        wa = max(d9_area, key=d9_area.get)
        print(
            f"[IMP-6 s4] -n {comm.size}: {n_cells} cells (predicted 1 149 150 + [35 k, 100 k]; "
            f"+{n_cells - 1149150} over C2); build {t_build:.1f} s (gmsh "
            f"{diag.get('mesh_wall_time_s', float('nan')):.1f} s; predicted ≈ 400 s [380, 510]); "
            f"summed peak RSS {rss:.3f} GiB\n"
            f"[IMP-6 s4] census {census}\n"
            f"[IMP-6 s4] (D1) body {body:.12e} vs V {V:.12e}: rel {d1_mass:.3e}; OCC wall "
            f"{float(sph['wall_area_m2']):.12e} vs A {A:.12e}: rel {d1_area:.3e} "
            f"(band {CAD_ANALYTIC_BAND:.0e})\n"
            f"[IMP-6 s4] (D2) CAD(403)/body {d2_403:.3e}; kept/body {d2_kept:.3e}; CAD(3) "
            f"{cad['phantom']:.15e} vs L^3 - body {L3 - body:.15e}: rel {d2_3:.3e}; "
            f"CAD(3)+CAD(403) vs L^3 {L3:.15e}: rel {d2_sum:.3e} (band {STRIP_BAND:.0e})\n"
            f"[IMP-6 s4] (N1) CAD(3) vs C2's L^3 misses by {n1_miss:.6e} (predicted "
            f"{n1_pred:.6e}; {n1_miss / STRIP_BAND:.3e}x STRIP_BAND)\n"
            f"[IMP-6 s4] (D3) meshed tag-3 {vol3:.15e} + tag-403 {vol403:.15e} vs L^3: rel "
            f"{d3_vol:.3e}; tag-3 extents lo {lo.tolist()} hi {hi.tolist()}: worst rel "
            f"{d3_ext:.3e}\n"
            f"[IMP-6 s4] (D4) CAD groups {cad_sum:.15e} + conductor {m_cond:.15e} vs can "
            f"{m_air:.15e}: rel {d4:.3e} (groups {cad})\n"
            f"[IMP-6 s4] (N2) pre-repair reading (+ body) misses by {n2_miss:.6e} (predicted "
            f"4.350023e-07; {n2_miss / STRIP_BAND:.3e}x STRIP_BAND)\n"
            f"[IMP-6 s4] (D5) 402 owned {n402}, exterior {n402_ext}, interior {n402_int}, "
            f"3|403 pairs {n402_pair}; meshed dS(402)/A {d5_ratio:.6f} (band "
            f"{TERMINAL_AREA_BAND})\n"
            f"[IMP-6 s4] (D6) OCC wall centroid {c} vs {SPHEROID_CENTRE}: max dev {d6:.3e} "
            f"(bound {EXACT * SPHEROID_SEMI_MAJOR:.1e}); meshed "
            f"{tuple(m / area402 for m in moments) if area402 else None}\n"
            f"[IMP-6 s4] (D7) rim census n_mask {n_mask}, inside {n_in} (min s {s_min:.6f}); "
            f"n403 {n403_mid}, outside {n_out} (max s {s_max:.9f})\n"
            f"[IMP-6 s4] (D8) census failures {census_fails or 'none'}; implant_rod_cad "
            f"{diag.get('implant_rod_cad')}\n"
            f"[IMP-6 s4] (D9) owned sheet facets == N_gmsh on {len(ring_ports) - len(d9_bad)}"
            f"/{len(ring_ports)}; mismatches {d9_bad or 'none'}; N_gmsh "
            f"{min(n_gmsh.values())}..{max(n_gmsh.values())}, total {sum(n_gmsh.values())}; "
            f"worst area P{wa} {d9_area[wa]:.3e}\n"
            f"[IMP-6 s4] (N4) tag-403 cells {census.get(BIRDCAGE_IMPLANT_CELL_TAG, 0)} "
            f"(F-human 14 047); meshed 403 / V {vol403 / V:.6f} (F-human 0.995984); "
            f"conductor vs C2 {C2_CONDUCTOR_CAD_MASS_M3:.15e}: rel {n4_cond:.3e} "
            f"(predicted <= 1e-9; printed only)\n"
            f"[IMP-6 s4] unknowns E {E}, F {F}; U1 = E = {E}, U2 = 2E + 2F = {2 * E + 2 * F}",
            flush=True,
        )

    # Negative controls (N1), (N2).
    assert n1_miss > STRIP_BAND, n1_miss
    assert n2_miss > STRIP_BAND, n2_miss
    # (D1)
    assert d1_mass < CAD_ANALYTIC_BAND, (body, V)
    assert d1_area < CAD_ANALYTIC_BAND, (sph["wall_area_m2"], A)
    # (D2)
    assert d2_403 < STRIP_BAND and d2_kept < STRIP_BAND, (d2_403, d2_kept)
    assert d2_3 < STRIP_BAND, d2_3
    assert d2_sum < STRIP_BAND, d2_sum
    # (D3)
    assert d3_vol <= STRIP_BAND, (vol3, vol403, L3)
    assert d3_ext <= STRIP_BAND, (lo, hi)
    # (D4)
    assert d4 <= STRIP_BAND, (cad_sum, m_cond, m_air)
    # (D5)
    low, high = TERMINAL_AREA_BAND
    assert n402 > 0 and n402_ext == 0, (n402, n402_ext)
    assert n402_int == n402 and n402_pair == n402, (n402, n402_int, n402_pair)
    assert low < d5_ratio < high, d5_ratio
    # (D6)
    assert d6 <= EXACT * SPHEROID_SEMI_MAJOR, c
    # (D7)
    assert n_mask > 0 and n_in == 0, (n_mask, n_in)
    assert n403_mid > 0 and n_out == 0, (n403_mid, n_out)
    # (D8)
    assert not census_fails, census_fails
    assert PHANTOM_CELL_TAG in census and BIRDCAGE_IMPLANT_CELL_TAG in census, census
    assert 1 not in census, census
    assert diag.get("implant_rod_cad") is None
    # (D9)
    assert len(n_gmsh) == 2 * LEG_COUNT, sorted(n_gmsh)
    assert not d9_bad, d9_bad
    for i, dev in d9_area.items():
        assert dev <= STRIP_BAND, (i, dev)
