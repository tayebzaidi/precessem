"""`GEO-38` C20 (§9 item 716, §10 chain C): the ASTM box and the kept 1 mm x
30 mm rod inside F-clinical-3T — CAD and mesh identities only, no solve.

C8's box kwargs (``_box_kwargs()`` on the 3 T r1 table, the route of
`test_geo38_f_clinical_3t.py:284–295`) plus ``implant_rod = dict(ROD,
keep_cells=True, centre=(x, off, 0))`` with ``x = 0.25 * Lx`` = 0.105 m and
``off = table_offset_m``, as `test_geo36_rod_f_clinical._rod_centre` / `_kwargs`
(`:88–97`) do on the 1.5 T box: I33's placement, a modelling choice, no
ASTM-placement claim. Those two read the 1.5 T ``_box_kwargs`` and are not
called here. Brief: docs/planning/chunks/GEO-38.md, "§9 item 716".

The build test is selected by environment only (``FEM_EM_GEO38_C20=1``), never
by ``-k``; it skips when unset. The coding tests are pure and always run.

Anchors (imported bands, unmoved). The (G1), (G2), (G4) and (G7) readings are
copied from `test_geo36_rod_f_clinical.py:215–224`, `:266–269`, `:283`, `:364–369`,
`:374` and `:383–386`, which are inline there:
  (G0) cells < F_CLINICAL_ROD_CELL_GUARD before any facet reading;
  (G1) rod_mass = kept_mass = V_rod and OCC wall area = A_rod at STRIP_BAND;
  (G2) CAD(403) = V_rod, CAD(3) = L³ − V_rod, their sum = L³ at STRIP_BAND;
  (G4) Σ cad_mass_by_group + conductor_cad_mass_m3 = the 3 T can at STRIP_BAND;
  (G7) census clean with tags 3 and 403, no tag 1, implant_spheroid_cad None;
  (64) 64 ring sheets in the census and in ring_sheet_gmsh_triangles;
  (T3) shield Ø 0.743, length 0.70, rung Ø 0.704 at IDENTITY_BAND, as C11 reads it;
  (Bore) every box corner's radius from the coil axis < the rung radius 0.352 m.
Negative control (asserted): `F_CLINICAL_ROD_CELL_RECORD`, the 1.5 T rod record,
misses the 3 T count by > CELL_COUNT_BAND.
Scope: the mesh and its identities. Not claimed: a solve, a tuned coil, a price,
a 3 T gel value, an ASTM placement.
"""

from __future__ import annotations

import os
import resource
import time

import numpy as np
import pytest
from mpi4py import MPI

import tests.validation.test_geo38_f_clinical_3t as t3
import tests.validation.test_geo38_f_clinical_box as s2
from tests.mesh.test_birdcage_implant_rod import ROD, _rod_area, _rod_volume
from tests.validation.test_geo36_rod_f_clinical import (
    F_CLINICAL_ROD_CELL_GUARD,
    F_CLINICAL_ROD_CELL_RECORD,
)
from tests.validation.test_geo36_rod_f_clinical import _kwargs as _kwargs_15t
from tests.validation.test_geo38_f_clinical_build import (
    IDENTITY_BAND,
    LEG_COUNT,
    _refuse_before_gmsh,
    _ring_ports,
)

C20_ENV = "FEM_EM_GEO38_C20"
#: Variables a C20 window must never read (the 1.5 T rod and neighbouring knobs).
FOREIGN_ENVS = (
    "FEM_EM_GEO36_ROD_BOX",
    "FEM_EM_GEO36_S8_COLUMN",
    "FEM_EM_GEO36_S8_DRIVE",
    "FEM_EM_GEO36_S8_DRIVE_DEGREE",
    "FEM_EM_GEO36_S8_DRIVE_C_PF",
    "FEM_EM_GEO38_C8",
    "FEM_EM_GEO38_C11",
    "FEM_EM_GEO38_STEP1A",
    "FEM_EM_GEO38_STEP1B",
    "FEM_EM_GEO38_STEP2",
    "FEM_EM_GMSH_STOP_RULE",
)

PHANTOM_CELL_TAG = 3
RUNG_RADIUS_3T_M = 0.5 * t3.RUNG_DIAMETER_3T_M

#: Rod and box on the 3 T table, version-tagged ((1*) licence), default image,
#: 2026-11-01: `20261101T180220Z_GEO-38-C20-W1-build.log:11920`, `:11929` (-n 2).
F_CLINICAL_3T_ROD_CELL_RECORD = 1_128_430
# (U1 = E, U2 = 2E + 2F): E 1 508 256, F 2 383 258.
F_CLINICAL_3T_ROD_UNKNOWNS_RECORD = (1_508_256, 7_783_028)

#: Printed, never asserted (brief).
PREDICTED_BUILD_S = (430.0, 520.0)


def _say(msg):
    if MPI.COMM_WORLD.rank == 0:
        print(f"[GEO-38 C20] {msg}", flush=True)


def _guard_foreign():
    bad = [k for k in FOREIGN_ENVS if os.environ.get(k) is not None]
    if bad:
        raise ValueError(f"C20 never reads these variables; unset them: {bad}")


def _mode():
    raw = os.environ.get(C20_ENV, "").strip()
    if raw == "":
        return None
    if raw == "1":
        return raw
    raise ValueError(f"{C20_ENV} must be unset or '1', got {raw!r}")


def _box_kwargs_3t():
    """C8's route: `_box_kwargs()` with the 3 T table in its namespace."""
    with pytest.MonkeyPatch.context() as m:
        m.setattr(s2, "_f_clinical_kwargs", t3._f_clinical_3t_kwargs)
        return s2._box_kwargs()


def _rod_centre_3t():
    box = _box_kwargs_3t()["phantom_box"]
    return (0.25 * float(box["size_m"][0]), float(box["table_offset_m"]), 0.0)


def _kwargs_3t():
    kw = _box_kwargs_3t()
    kw["implant_rod"] = dict(ROD, keep_cells=True, centre=_rod_centre_3t())
    return kw


def _bore_gap(kw):
    """Rung radius less the largest box-corner radius from the coil axis."""
    Lx, Ly, _Lz = kw["phantom_box"]["size_m"]
    off = kw["phantom_box"]["table_offset_m"]
    corners = [
        np.hypot(sx * 0.5 * Lx, sy * 0.5 * Ly + off) for sx in (-1.0, 1.0) for sy in (-1.0, 1.0)
    ]
    return float(kw["ring_radius"] - max(corners)), float(max(corners))


# --- Coding tests: pure, no gmsh. ---------------------------------------------


def test_c20_k0_knob_and_foreign_guard(monkeypatch):
    for k in FOREIGN_ENVS + (C20_ENV,):
        monkeypatch.delenv(k, raising=False)
    assert _mode() is None
    _guard_foreign()
    monkeypatch.setenv(C20_ENV, "1")
    assert _mode() == "1"
    for bad in ("2", "build"):
        monkeypatch.setenv(C20_ENV, bad)
        with pytest.raises(ValueError):
            _mode()
    for k in FOREIGN_ENVS:
        with monkeypatch.context() as m:
            m.setenv(k, "1")
            with pytest.raises(ValueError):
                _guard_foreign()


def test_c20_k1_kwargs_are_the_3t_table_with_box_and_rod():
    kw = _kwargs_3t()
    kw15 = _kwargs_15t()
    for key, val in t3.F_CLINICAL_3T_OVERRIDES.items():
        assert kw[key] == val, key
    assert kw["ring_radius"] != kw15["ring_radius"]
    assert kw["coil_length"] != kw15["coil_length"]
    assert kw["shield_radius"] != kw15["shield_radius"]
    assert kw["phantom_box"] == kw15["phantom_box"]
    assert tuple(kw["phantom_box"]["size_m"]) == (0.42, 0.09, 0.65)
    assert "phantom_radius" not in kw and "phantom_height" not in kw
    rod = kw["implant_rod"]
    assert rod == kw15["implant_rod"]
    assert rod["keep_cells"] is True
    assert abs(rod["centre"][0] - 0.105) < 1e-15 and rod["centre"][1:] == (0.0, 0.0)
    # The swap is scoped: `_box_kwargs()` reads the 1.5 T table again afterwards.
    assert s2._box_kwargs()["ring_radius"] == kw15["ring_radius"]


def test_c20_k2_bore():
    gap, corner = _bore_gap(_kwargs_3t())
    _say(
        f"(Bore) largest box-corner radius {corner:.9f} m, rung radius {RUNG_RADIUS_3T_M} m, "
        f"gap {gap:.9f} m (ASSERTED > 0)"
    )
    assert gap > 0.0, (gap, corner)
    # The control: a box as wide as the rung circle does not fit.
    wide = dict(
        ring_radius=RUNG_RADIUS_3T_M,
        phantom_box=dict(size_m=(2.0 * RUNG_RADIUS_3T_M, 0.09, 0.65), table_offset_m=0.0),
    )
    assert _bore_gap(wide)[0] < 0.0


def test_c20_k3_kwargs_reach_the_builder(monkeypatch):
    # The tripwire fires on rank 0; `birdcage_port_domain` reports it elsewhere.
    if MPI.COMM_WORLD.rank == 0:
        expected = dict(expected_exception=AssertionError, match="reached the gmsh builder")
    else:
        expected = dict(expected_exception=RuntimeError, match="generation failed on rank 0")
    with pytest.raises(**expected):
        _refuse_before_gmsh(monkeypatch, **_kwargs_3t())


# --- The build (environment-selected). ----------------------------------------


def test_geo38_c20_rod_in_3t_box(monkeypatch):
    if _mode() is None:
        pytest.skip(f"{C20_ENV} unset: GEO-38 C20 is selected by environment only")
    _guard_foreign()
    from dolfinx import fem

    from fem_em_solver.io.mesh import BIRDCAGE_IMPLANT_CELL_TAG, MeshGenerator
    from tests.mesh.test_birdcage_port_sheet_prerequisite import CELL_COUNT_BAND
    from tests.validation.test_ops61_exact_strip_ports import STRIP_BAND
    from tests.validation.test_wf12_f_human_hole import _census_checks, _owned_census

    comm = MPI.COMM_WORLD
    kw = _kwargs_3t()
    ring_ports = _ring_ports()
    Lx, Ly, Lz = kw["phantom_box"]["size_m"]
    L3 = Lx * Ly * Lz
    V = _rod_volume()
    A = _rod_area()
    gap, corner = _bore_gap(kw)
    _say(
        f"rod centre {kw['implant_rod']['centre']}, axis {ROD['axis']}; V_rod {V:.12e}, A_rod "
        f"{A:.12e}, L^3 {L3:.15e}; (Bore) corner {corner:.9f} m vs rung radius "
        f"{kw['ring_radius']} m, gap {gap:.9f} m; PREDICTED build {PREDICTED_BUILD_S} s; "
        f"/proc/loadavg {t3._loadavg()}"
    )

    t0 = time.perf_counter()
    with t3._GmshStopRule(monkeypatch, t3.GMSH_STOP_S, "C20"):
        msh, cells, _facets, diag = MeshGenerator.birdcage_port_domain(
            comm=comm, return_diagnostics=True, **kw
        )
    t_build = time.perf_counter() - t0
    tdim = msh.topology.dim
    fdim = tdim - 1
    n_cells = int(msh.topology.index_map(tdim).size_global)
    _say(
        f"build done {t_build:.1f} s at -n {comm.size}; gmsh "
        f"{diag.get('mesh_wall_time_s', float('nan')):.2f} s; (G0) {n_cells} cells vs guard "
        f"{F_CLINICAL_ROD_CELL_GUARD}; /proc/loadavg {t3._loadavg()}"
    )
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

    # --- (G4): the 3 T can's partition with the rod.
    cad_sum = float(sum(cad.values()))
    m_cond = float(diag["conductor_cad_mass_m3"])
    m_air = float(diag["air_domain_cad_mass_m3"])
    g4 = abs((cad_sum + m_cond) / m_air - 1.0)

    # --- (G7), (64): census.
    census_fails = _census_checks(census, n_cells, ring_ports)
    n_gmsh = {int(k): int(v) for k, v in diag["ring_sheet_gmsh_triangles"].items()}

    # --- (T3), as C11 reads it.
    shield_d = 2.0 * float(diag["shield_radius_m"])
    shield_l = float(diag["shield_length_m"])
    rung_d = 2.0 * kw["ring_radius"]

    # --- Unknowns: U1 from a function-space build, U2 = 2E + 2F.
    msh.topology.create_entities(1)
    msh.topology.create_entities(2)
    E = comm.allreduce(msh.topology.index_map(1).size_local, op=MPI.SUM)
    F = comm.allreduce(msh.topology.index_map(2).size_local, op=MPI.SUM)
    V1 = fem.functionspace(msh, ("N1curl", 1))
    U1 = int(V1.dofmap.index_map.size_global * V1.dofmap.index_map_bs)
    U2 = 2 * E + 2 * F
    rss = comm.allreduce(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0**2, op=MPI.SUM)

    # --- The negative control, and the 3 T box beside it (printed).
    miss_15t = abs(n_cells / F_CLINICAL_ROD_CELL_RECORD - 1.0)
    added = n_cells - t3.F_CLINICAL_3T_BOX_CELL_RECORD
    _say(
        f"-n {comm.size}: {n_cells} cells; beside F_CLINICAL_3T_BOX_CELL_RECORD "
        f"{t3.F_CLINICAL_3T_BOX_CELL_RECORD}: rod adds {added} (I33 added ≈ 70 000 on 1.5 T; "
        f"printed); summed peak RSS {rss:.3f} GiB\n"
        f"[GEO-38 C20] census {census}\n"
        f"[GEO-38 C20] (G1) rod_mass {rod_mass:.12e} / kept {kept:.12e} vs V_rod {V:.12e}: rel "
        f"{g1_mass:.3e} / {g1_kept:.3e}; OCC wall {float(rod['wall_area_m2']):.12e} vs A_rod "
        f"{A:.12e}: rel {g1_area:.3e} (band {STRIP_BAND:.0e})\n"
        f"[GEO-38 C20] (G2) CAD(403)/V {g2_403:.3e}; CAD(3) {cad['phantom']:.15e} vs L^3 - V "
        f"{L3 - V:.15e}: rel {g2_3:.3e}; CAD(3)+CAD(403) vs L^3: rel {g2_sum:.3e}\n"
        f"[GEO-38 C20] (G4) CAD groups {cad_sum:.15e} + conductor {m_cond:.15e} vs can "
        f"{m_air:.15e}: rel {g4:.3e} (groups {cad})\n"
        f"[GEO-38 C20] (G7) census failures {census_fails or 'none'}; implant_spheroid_cad "
        f"{diag.get('implant_spheroid_cad')}; tag 403 cells "
        f"{census.get(BIRDCAGE_IMPLANT_CELL_TAG, 0)}\n"
        f"[GEO-38 C20] (64) ring ports {len(ring_ports)}, gmsh sheets {len(n_gmsh)}, "
        f"N_gmsh {min(n_gmsh.values())}..{max(n_gmsh.values())}, total {sum(n_gmsh.values())}\n"
        f"[GEO-38 C20] (T3) shield Ø {shield_d:.12f} (vs {t3.SHIELD_DIAMETER_3T_M}), length "
        f"{shield_l:.12f} (vs {t3.SHIELD_LENGTH_3T_M}), rung Ø {rung_d:.12f} (vs "
        f"{t3.RUNG_DIAMETER_3T_M}) (IDENTITY_BAND {IDENTITY_BAND:.0e})\n"
        f"[GEO-38 C20] (N) |{n_cells}/{F_CLINICAL_ROD_CELL_RECORD} - 1| = {miss_15t:.3e} = "
        f"{miss_15t / CELL_COUNT_BAND:.2f}x CELL_COUNT_BAND {CELL_COUNT_BAND}\n"
        f"[GEO-38 C20] unknowns E {E}, F {F}; U1 {U1} (function space), U2 = 2E + 2F = {U2}; "
        f"records {F_CLINICAL_3T_ROD_CELL_RECORD}, {F_CLINICAL_3T_ROD_UNKNOWNS_RECORD}"
    )

    # Negative control: the 1.5 T rod record misses, so the 3 T table reached gmsh.
    assert miss_15t > CELL_COUNT_BAND, miss_15t
    # (Bore)
    assert gap > 0.0, (gap, corner)
    # (T3)
    assert abs(shield_d / t3.SHIELD_DIAMETER_3T_M - 1.0) <= IDENTITY_BAND, shield_d
    assert abs(shield_l / t3.SHIELD_LENGTH_3T_M - 1.0) <= IDENTITY_BAND, shield_l
    assert abs(rung_d / t3.RUNG_DIAMETER_3T_M - 1.0) <= IDENTITY_BAND, rung_d
    # (G1)
    assert g1_mass < STRIP_BAND and g1_kept < STRIP_BAND, (rod_mass, kept, V)
    assert g1_area < STRIP_BAND, (rod["wall_area_m2"], A)
    # (G2)
    assert g2_403 < STRIP_BAND, g2_403
    assert g2_3 < STRIP_BAND, g2_3
    assert g2_sum < STRIP_BAND, g2_sum
    # (G4)
    assert g4 <= STRIP_BAND, (cad_sum, m_cond, m_air)
    # (G7), (64)
    assert not census_fails, census_fails
    assert PHANTOM_CELL_TAG in census and BIRDCAGE_IMPLANT_CELL_TAG in census, census
    assert 1 not in census, census
    assert diag.get("implant_spheroid_cad") is None
    assert len(ring_ports) == 2 * LEG_COUNT == 64
    assert sorted(n_gmsh) == sorted(ring_ports), sorted(n_gmsh)
    # Unknowns; the records bind once registered.
    assert U1 == E, (U1, E)
    if F_CLINICAL_3T_ROD_CELL_RECORD is not None:
        assert abs(n_cells / F_CLINICAL_3T_ROD_CELL_RECORD - 1.0) <= CELL_COUNT_BAND, n_cells
        assert (U1, U2) == F_CLINICAL_3T_ROD_UNKNOWNS_RECORD, (U1, U2)
