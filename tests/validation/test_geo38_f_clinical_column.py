"""GEO-38 step 1b (§9 items 469 / 485, chain C1 (b)): one degree-1 copper column
on F-clinical-1.5T at ``-n 8``, 64 MHz.

Selected by ``FEM_EM_GEO38_STEP1B`` only:

* ``=1`` — attempt 1's in-memory build + column route, unchanged (the route an
  ``xl`` window would lift).
* ``=cache`` — item 485 window A: build r1 once, write mesh + cell tags + facet
  tags to ``output/geo38_s1b_cache/f_clinical_r1.xdmf`` (+ ``.h5``, JSON
  sidecar) with `TH-11` step 5a's idiom, read it back in the same window; the
  round trip is an **exact integer identity** (global cells = the record
  1 147 366, per-tag owned cell and facet counts equal; review ruling
  2026-10-15 18:00).
* ``=column`` — item 485 window B: read the cache with
  ``GhostMode.shared_facet`` (the 64 sheets are interior facets), check it
  against the sidecar (same exact identity, and the sidecar's versions equal
  the running image), rebuild the context with the same helper the in-memory
  route uses, and run attempt 1's body from "(A0)" on.  Never builds: a
  missing cache or sidecar fails.

The fixture is item 478's r1 build (``F_CLINICAL_15T_KWARGS`` /
``_f_clinical_kwargs("r1")``) on the copper hole route: 64 ring-gap sheets at
``TERMINATED_PORT_IMPEDANCE_OHM``, the Leontovich wall on tag 401 with
Z_s = (1 + j)·R_s(ω, 5.8e7) (first in the factorising call's bilinear list,
`WF-12.md:109`), the shield can PEC through ``_outer_box_tags``.  P33
factorises, P34 is held on the same factor.

Context builder: a local F-clinical copy of ``_build_f_human``'s ``as_hole``
branch (`test_wf7_f_human_b1_quadrature.py:155–248`) — a copy, not an edit.

Price rule (item 469 step 0, evaluated before code from 478's record):
N₁ = 1 147 366, ρ = N₁/369 387 = 3.1061, t_b = 323.3 s (478 W1) ⇒
T = 323.3 + 101.6 + 13.1 + 92.1 = 530.0 s ≤ 540 and M = 46.2 GiB ≤ 64 ⇒ GO.
Conservative column (printed, not binding) 39.2·(U₁/500 568)² = 363.9 s.
"""

from __future__ import annotations

import json
import math
import os
import time
from pathlib import Path

import numpy as np
import pytest
from mpi4py import MPI

from tests.complex_mode import complex_only

STEP1B_ENV = "FEM_EM_GEO38_STEP1B"
COPPER_SIGMA = 5.8e7
DRIVEN = "P33"
HELD = "P34"
F_HUMAN_N1 = 369_387
F_HUMAN_COL_S = 39.2
F_HUMAN_HELD_S = 2.88
COL_EXPONENT = 0.84
PREDICTED_T_S = 530.0
PREDICTED_M_GIB = 46.2
PREDICTED_BLIND_FACTOR = 100.0  # rule (e), printed: F-human read 123.3x (`WF-12.md:142`)
PREDICTED_COIL_SHARE = 0.370  # F-human, printed
PREDICTED_RECIPROCITY = 1e-14  # F-human 9.167e-15, printed
ATTEMPT1_T_BUILD_N8_S = 361.6  # item 469 W1 build at -n 8, printed beside the sidecar's t_b


def _say(comm, msg):
    if comm.rank == 0:
        print(f"[GEO-38 s1b] {msg}", flush=True)


def _selector() -> str:
    return os.environ.get(STEP1B_ENV, "").strip()


def _sheet_tags_and_widths(msh, cell_tags, chord, comm):
    """The 64 interior sheet tags and widths area / chord (the fixture's route)."""
    from fem_em_solver.io.mesh import _interface_facet_tags
    from tests.mesh.test_birdcage_port_sheets import PORT_LOWER, PORT_UPPER, SHEET_IFACE
    from tests.mesh.test_birdcage_port_terminals import _interface_area_or_zero
    from tests.validation.test_geo38_f_clinical_build import LEG_COUNT, _ring_ports

    ring_ports = _ring_ports()
    assert len(ring_ports) == 2 * LEG_COUNT
    tags_f = _interface_facet_tags(
        msh, cell_tags, {SHEET_IFACE + i: (PORT_LOWER + i, PORT_UPPER + i) for i in ring_ports}
    )
    tdim = msh.topology.dim
    msh.topology.create_connectivity(tdim - 1, tdim)
    msh.topology.create_entity_permutations()
    widths = {
        int(i): _interface_area_or_zero(msh, tags_f, SHEET_IFACE + i, comm) / chord
        for i in ring_ports
    }
    return tags_f, widths


def _context_from_mesh(msh, cell_tags, facets, chord, comm):
    """`_build_f_clinical`'s lines after the `birdcage_port_domain` call — the
    one helper both the in-memory (=1) and the cached (=column) route call."""
    from fem_em_solver.core import HomogeneousMaterial, TimeHarmonicProblem, TimeHarmonicSolver
    from fem_em_solver.ports.lumped import LumpedSheetPortSpec
    from tests.mesh.test_birdcage_port_sheets import SHEET_IFACE
    from tests.mesh.test_birdcage_ring_gaps_scaleup import _ring_gap_frame
    from tests.validation.test_geo38_f_clinical_build import LEG_COUNT
    from tests.validation.test_lossy_sphere_fullwave import SALINE_EPSILON_R, SALINE_SIGMA
    from tests.validation.test_port_birdcage_four_port import TERMINATED_PORT_IMPEDANCE_OHM
    from tests.validation.test_port_birdcage_lumped_column import PHANTOM_CELL_TAG
    from tests.validation.test_th14_birdcage_copper import OUTER_BOX_TAG, _outer_box_tags
    from tests.validation.test_wf7_f_human_b1_quadrature import DEGREE, FREQUENCY_HZ

    tags_f, widths = _sheet_tags_and_widths(msh, cell_tags, chord, comm)
    specs = []
    for i, width in widths.items():
        phi_hat, _centre = _ring_gap_frame(i, LEG_COUNT)
        specs.append(
            LumpedSheetPortSpec(
                port_id=f"P{i}",
                facet_tag=int(SHEET_IFACE + i),
                port_impedance_ohm=TERMINATED_PORT_IMPEDANCE_OHM,
                gap_height_m=chord,
                sheet_width_m=width,
                drive_direction=tuple(float(c) for c in phi_hat),
                drive_voltage_v=1.0 + 0.0j,
                interior=True,
            )
        )
    outer_tags, outer_census = _outer_box_tags(msh, facets, comm)
    problem = TimeHarmonicProblem(
        mesh=msh,
        frequency_hz=FREQUENCY_HZ,
        material=HomogeneousMaterial(sigma=0.0, epsilon_r=1.0, mu_r=1.0),
        cell_tags=cell_tags,
        material_map={
            PHANTOM_CELL_TAG: HomogeneousMaterial(
                sigma=SALINE_SIGMA, epsilon_r=SALINE_EPSILON_R, mu_r=1.0
            )
        },
        boundary_condition="pec_zero_tangential_a",
        facet_tags=outer_tags,
        pec_facet_tags=(OUTER_BOX_TAG,),
    )
    ctx = {
        "msh": msh,
        "tags_f": tags_f,
        "cell_tags": cell_tags,
        "omega": 2.0 * np.pi * FREQUENCY_HZ,
        "specs": specs,
        "solver": TimeHarmonicSolver(problem, degree=DEGREE),
    }
    return ctx, outer_census, widths


def _build_f_clinical(comm):
    """F-clinical copy of `_build_f_human`'s ``as_hole`` branch (leg_count 32,
    ring ordinals 33..96, ``_ring_gap_frame(i, 32)``) — the in-memory route."""
    from fem_em_solver.io.mesh import MeshGenerator
    from tests.validation.test_geo38_f_clinical_build import _f_clinical_kwargs

    kwargs = _f_clinical_kwargs("r1")
    msh, cell_tags, facets, diag = MeshGenerator.birdcage_port_domain(
        comm=comm, return_diagnostics=True, **kwargs
    )
    n_cells = int(msh.topology.index_map(3).size_global)
    chord = float(diag["ring_port_layout"]["ring_port_gap_chord_m"])
    ctx, outer_census, _widths = _context_from_mesh(msh, cell_tags, facets, chord, comm)
    return {
        "ctx": ctx,
        "cells": n_cells,
        "facets": facets,
        "shield_radius_m": float(diag["shield_radius_m"]),
        "outer_census": outer_census,
    }


# ---- item 485: the disk cache (window A writes, window B reads) -------------
# `TH-11` step 5a's idiom (`test_coil_loading_larmor_mesh_cache.py:196–216`),
# except that the read is `GhostMode.shared_facet`: the 64 sheets are interior
# facets and need the cell on both sides (`io/mesh.py:4905–4917`).  Named
# grids with an explicit geometry xpath; `create_connectivity(2, 3)` before the
# facet write and before the facet read.
CACHE_DIR_ENV = "FEM_EM_GEO38_CACHE_DIR"
MESH_NAME = "mesh"
CELL_TAGS_NAME = "cell_tags"
FACET_TAGS_NAME = "facet_tags"
_GEOMETRY_XPATH = f"/Xdmf/Domain/Grid[@Name='{MESH_NAME}']/Geometry"


def _cache_path() -> Path:
    root = Path(os.environ.get(CACHE_DIR_ENV, "output/geo38_s1b_cache"))
    return root / "f_clinical_r1.xdmf"


def _owned_tag_counts(msh, tags, comm):
    """Global per-tag counts over **owned** entities only (5a's rule: ghosts
    would be double-counted on a ghosted mesh)."""
    size_local = msh.topology.index_map(tags.dim).size_local
    owned = np.asarray(tags.indices) < size_local
    values = np.asarray(tags.values, dtype=np.int64)[owned]
    local = {int(v): int(c) for v, c in zip(*np.unique(values, return_counts=True))}
    merged = {}
    for part in comm.allgather(local):
        for v, c in part.items():
            merged[v] = merged.get(v, 0) + c
    return merged


def _cache_signature(msh, cell_tags, facets, comm):
    return {
        "cells": int(comm.allreduce(msh.topology.index_map(3).size_local, op=MPI.SUM)),
        "cell_counts": _owned_tag_counts(msh, cell_tags, comm),
        "facet_counts": _owned_tag_counts(msh, facets, comm),
    }


def _sig_from_json(d):
    return {
        "cells": int(d["cells"]),
        "cell_counts": {int(k): int(v) for k, v in d["cell_counts"].items()},
        "facet_counts": {int(k): int(v) for k, v in d["facet_counts"].items()},
    }


def _versions():
    import dolfinx
    import gmsh

    return {"dolfinx": str(dolfinx.__version__), "gmsh": str(gmsh.__version__)}


def _read_cache(comm):
    """Read the cached mesh with `shared_facet` ghosts -> (msh, cell_tags, facets, t_read)."""
    from dolfinx import io
    from dolfinx import mesh as dmesh

    path = _cache_path()
    comm.Barrier()
    t0 = time.perf_counter()
    with io.XDMFFile(comm, str(path), "r") as xdmf:
        msh = xdmf.read_mesh(ghost_mode=dmesh.GhostMode.shared_facet, name=MESH_NAME)
        msh.topology.create_connectivity(msh.topology.dim - 1, msh.topology.dim)
        cell_tags = xdmf.read_meshtags(msh, name=CELL_TAGS_NAME)
        facets = xdmf.read_meshtags(msh, name=FACET_TAGS_NAME)
    comm.Barrier()
    return msh, cell_tags, facets, time.perf_counter() - t0


def _cache_size_mib(path: Path) -> float:
    return sum(p.stat().st_size for p in (path, path.with_suffix(".h5")) if p.exists()) / 1024.0**2


def test_geo38_step1b_cache_r1():
    """Window A (``=cache``): build r1 once, write it, read it back in the same
    window.  Exact integer identities only (no tolerance)."""
    if _selector() != "cache":
        pytest.skip(f"{STEP1B_ENV} != cache: window A is selected by environment only")
    from dolfinx import io

    from fem_em_solver.io.mesh import MeshGenerator
    from tests.validation.test_geo38_f_clinical_build import (
        F_CLINICAL_R1_CELL_RECORD,
        _f_clinical_kwargs,
    )

    comm = MPI.COMM_WORLD
    path = _cache_path()
    if comm.rank == 0:
        path.parent.mkdir(parents=True, exist_ok=True)
    comm.Barrier()
    t0 = time.perf_counter()
    msh, cell_tags, facets, diag = MeshGenerator.birdcage_port_domain(
        comm=comm, return_diagnostics=True, **_f_clinical_kwargs("r1")
    )
    comm.Barrier()
    t_build = time.perf_counter() - t0
    chord = float(diag["ring_port_layout"]["ring_port_gap_chord_m"])
    written = _cache_signature(msh, cell_tags, facets, comm)
    _say(
        comm,
        f"cache build: {t_build:.1f} s at -n {comm.size}; {written['cells']} cells; "
        f"cell tags {written['cell_counts']}; facet tags {written['facet_counts']}",
    )
    _tags_f, widths = _sheet_tags_and_widths(msh, cell_tags, chord, comm)

    msh.name = MESH_NAME
    cell_tags.name = CELL_TAGS_NAME
    facets.name = FACET_TAGS_NAME
    msh.topology.create_connectivity(msh.topology.dim - 1, msh.topology.dim)
    t1 = time.perf_counter()
    with io.XDMFFile(comm, str(path), "w") as xdmf:
        xdmf.write_mesh(msh)
        xdmf.write_meshtags(cell_tags, msh.geometry, geometry_xpath=_GEOMETRY_XPATH)
        xdmf.write_meshtags(facets, msh.geometry, geometry_xpath=_GEOMETRY_XPATH)
    comm.Barrier()
    t_write = time.perf_counter() - t1
    if comm.rank == 0:
        side = {
            **written,
            "shield_radius_m": float(diag["shield_radius_m"]),
            "ring_port_gap_chord_m": chord,
            "t_build_s": t_build,
            "n_ranks": comm.size,
            "versions": _versions(),
            "sheet_widths_m": widths,
        }
        path.with_suffix(".json").write_text(json.dumps(side, indent=1, sort_keys=True))
    comm.Barrier()

    msh_r, ct_r, ft_r, t_read = _read_cache(comm)
    read = _cache_signature(msh_r, ct_r, ft_r, comm)
    _tf_r, widths_r = _sheet_tags_and_widths(msh_r, ct_r, chord, comm)
    w_rel = max(abs(widths_r[i] / widths[i] - 1.0) for i in widths)
    _say(
        comm,
        f"read-back sheet widths vs built: max rel {w_rel:.3e} over {len(widths)} "
        f"sheets (printed; catches a ghost-mode slip before window B)",
    )
    size_mib = _cache_size_mib(path) if comm.rank == 0 else 0.0
    _say(
        comm,
        f"cache {path} ({size_mib:.1f} MiB): write {t_write:.1f} s, read-back "
        f"{t_read:.1f} s (shared_facet, -n {comm.size}); read {read['cells']} cells; "
        f"cell tags equal {read['cell_counts'] == written['cell_counts']}; facet tags "
        f"equal {read['facet_counts'] == written['facet_counts']} (EXACT, ASSERTED); "
        f"cells vs record {F_CLINICAL_R1_CELL_RECORD} (EXACT, ASSERTED); versions "
        f"{_versions()}",
    )
    assert written["cells"] == F_CLINICAL_R1_CELL_RECORD, written["cells"]
    assert read == written, (read, written)


def _built_from_cache(comm):
    """Window B (``=column``): the cache, never a build.  Refuses on a missing
    cache/sidecar, a version mismatch, or any count mismatch."""
    from tests.validation.test_geo38_f_clinical_build import F_CLINICAL_R1_CELL_RECORD

    path = _cache_path()
    side_path = path.with_suffix(".json")
    present = path.exists() and path.with_suffix(".h5").exists() and side_path.exists()
    assert (
        present
    ), f"cache {path} or its .h5 / .json sidecar is missing: refusing (never build here)"
    side = json.loads(side_path.read_text())
    assert side["versions"] == _versions(), (side["versions"], _versions())
    msh, cell_tags, facets, t_read = _read_cache(comm)
    read = _cache_signature(msh, cell_tags, facets, comm)
    expected = _sig_from_json(side)
    _say(
        comm,
        f"cache read: {t_read:.1f} s at -n {comm.size} (shared_facet), "
        f"{_cache_size_mib(path):.1f} MiB; {read['cells']} cells; vs sidecar: cells "
        f"{read['cells'] == expected['cells']}, cell tags "
        f"{read['cell_counts'] == expected['cell_counts']}, facet tags "
        f"{read['facet_counts'] == expected['facet_counts']} (EXACT, ASSERTED); sidecar "
        f"t_b {side['t_build_s']:.1f} s at -n {side['n_ranks']}, versions {side['versions']}",
    )
    assert read == expected, (read, expected)
    assert read["cells"] == F_CLINICAL_R1_CELL_RECORD, read["cells"]
    chord = float(side["ring_port_gap_chord_m"])
    t0 = time.perf_counter()
    ctx, outer_census, widths = _context_from_mesh(msh, cell_tags, facets, chord, comm)
    w_a = {int(k): float(v) for k, v in side["sheet_widths_m"].items()}
    w_rel = max(abs(widths[i] / w_a[i] - 1.0) for i in w_a)
    _say(
        comm,
        f"context from cache: {time.perf_counter() - t0:.1f} s; sheet widths vs window A "
        f"max rel {w_rel:.3e} over {len(w_a)} sheets (printed)",
    )
    return {
        "ctx": ctx,
        "cells": int(msh.topology.index_map(3).size_global),
        "facets": facets,
        "shield_radius_m": float(side["shield_radius_m"]),
        "outer_census": outer_census,
        "t_read": t_read,
        "t_build_sidecar": float(side["t_build_s"]),
        "n_build_sidecar": int(side["n_ranks"]),
    }


@complex_only
def test_geo38_step1b_f_clinical_copper_column():
    mode = _selector()
    if mode not in ("1", "column"):
        pytest.skip(
            f"{STEP1B_ENV} not in (1, column): GEO-38 step 1b is selected by environment only"
        )
    import ufl
    from dolfinx import fem

    from fem_em_solver.core.cavity import surface_resistance_ohm
    from fem_em_solver.io.mesh import BIRDCAGE_CONDUCTOR_SURFACE_TAG
    from fem_em_solver.ports.lumped import (
        lumped_port_bilinear_term,
        lumped_port_linear_term,
        sheet_terminal_current,
    )
    from tests.mesh.test_birdcage_port_sheet_prerequisite import CELL_COUNT_BAND
    from tests.validation.test_birdcage_power_identity import (
        DISCRETE_IDENTITY_RTOL,
        _sheet_field_dissipation_w,
        _source_power_w,
    )
    from tests.validation.test_geo38_f_clinical_build import (
        F_CLINICAL_R1_CELL_RECORD,
        F_CLINICAL_R1_UNKNOWNS_RECORD,
        _facet_coords,
        _owned,
        _ring_ports,
    )
    from tests.validation.test_port_birdcage_four_port import TERMINATED_PORT_IMPEDANCE_OHM
    from tests.validation.test_port_birdcage_lumped_column import PHANTOM_CELL_TAG
    from tests.validation.test_port_lumped_sheet_sweep import RECIPROCITY_BAND
    from tests.validation.test_th14_birdcage_copper import _leontovich_term, _surface_loss_w
    from tests.validation.test_wf12_f_human_hole import _census_checks, _owned_census, _rss_gib

    comm = MPI.COMM_WORLD
    t_window = time.perf_counter()
    _say(
        comm,
        f"start ({STEP1B_ENV}={mode}): -n {comm.size}, degree 1, 64 MHz, drives "
        f"{DRIVEN} (factorise) + {HELD} (held); PREDICTED T {PREDICTED_T_S} s, "
        f"M {PREDICTED_M_GIB} GiB",
    )

    # ---- fixture, (A0), (iv) ---------------------------------------------
    t0 = time.perf_counter()
    built = _built_from_cache(comm) if mode == "column" else _build_f_clinical(comm)
    t_here = time.perf_counter() - t0
    if mode == "column":
        t_build = built["t_build_sidecar"]
        t_b_note = (
            f"t_b {t_build:.1f} s from the sidecar (window A, -n "
            f"{built['n_build_sidecar']}); item 469's in-memory build "
            f"{ATTEMPT1_T_BUILD_N8_S} s (-n 8); t_read {built['t_read']:.1f} s"
        )
    else:
        t_build = t_here
        t_b_note = f"t_b {t_build:.1f} s (-n {comm.size}, in-memory)"
    ctx = built["ctx"]
    msh, tags_f, omega, cell_tags = ctx["msh"], ctx["tags_f"], ctx["omega"], ctx["cell_tags"]
    facets = built["facets"]
    n_cells = built["cells"]
    rho = n_cells / F_HUMAN_N1
    _say(
        comm,
        f"fixture done: {t_here:.1f} s, cells {n_cells}, rho = N1/{F_HUMAN_N1} = "
        f"{rho:.4f}, t={time.perf_counter() - t_window:.1f} s; {t_b_note}",
    )
    census = _owned_census(msh, cell_tags, comm)
    census_fails = _census_checks(census, n_cells, _ring_ports())
    oc = built["outer_census"]
    cell_rel = n_cells / F_CLINICAL_R1_CELL_RECORD - 1.0
    R_s = float(built["shield_radius_m"])
    fv = np.asarray(facets.values)
    fi = np.asarray(facets.indices)
    cav_p = _facet_coords(msh, _owned(msh, fi[fv == BIRDCAGE_CONDUCTOR_SURFACE_TAG])).reshape(-1, 3)
    cav_hmax = comm.allreduce(
        float(np.hypot(cav_p[:, 0], cav_p[:, 1]).max()) if len(cav_p) else 0.0, op=MPI.MAX
    )
    _say(
        comm,
        f"(A0) cells {n_cells} vs record {F_CLINICAL_R1_CELL_RECORD}, rel {cell_rel:.3e} "
        f"(band {CELL_COUNT_BAND}, ASSERTED); (iv) census failures "
        f"{census_fails or 'none'}; outer census {oc}; (S) cavity max hypot "
        f"{cav_hmax:.9f} m < R_s {R_s} (ASSERTED)",
    )
    assert abs(cell_rel) <= CELL_COUNT_BAND, cell_rel
    assert not census_fails, census_fails
    assert oc["tag_collision"] == 0 and oc["cavity"] > 0, oc
    assert (
        oc["cavity_exterior"] == oc["cavity"] and oc["outer"] + oc["cavity"] == oc["exterior"]
    ), oc
    assert cav_hmax < R_s, (cav_hmax, R_s)

    z_s = (1.0 + 1.0j) * surface_resistance_ohm(omega, COPPER_SIGMA)
    wall = _leontovich_term(msh, facets, z_s, omega)
    z_p = complex(TERMINATED_PORT_IMPEDANCE_OHM)
    dx = ufl.Measure("dx", domain=msh, subdomain_data=cell_tags)
    acct = {"mesh": msh, "facet_tags": tags_f}
    solver = ctx["solver"]

    def _drive(driven_id, *, held):
        """`test_wf12_f_human_copper_matrix.py`'s `_drive`, wall term first."""
        sheets = [sp.sheet(driven=(sp.port_id == driven_id)) for sp in ctx["specs"]]
        dsheet = next(s for s in sheets if s.port_id == driven_id)
        linear = [
            lambda test, _s=dsheet: lumped_port_linear_term(
                msh, tags_f, _s, test, omega_rad_per_s=omega
            )
        ]
        comm.Barrier()
        t0 = time.perf_counter()
        if held:
            fields = solver.solve_with_held_factorization(extra_linear_terms=linear)
        else:
            bil = [wall] + [
                (
                    lambda trial, test, _s=s: lumped_port_bilinear_term(
                        msh, tags_f, _s, trial, test, omega_rad_per_s=omega
                    )
                )
                for s in sheets
            ]
            _say(comm, f"column start t={time.perf_counter() - t_window:.1f} s")
            fields = solver.solve(
                current_density=None,
                project_source=False,
                extra_bilinear_terms=bil,
                extra_linear_terms=linear,
            )
        comm.Barrier()
        t_solve = time.perf_counter() - t0
        _say(
            comm,
            f"{'held drive' if held else 'column (assembly + factorisation + solve)'} done "
            f"({driven_id}): {t_solve:.1f} s, t={time.perf_counter() - t_window:.1f} s",
        )
        if not held:
            _m, per, tot = _rss_gib(comm)
            _say(
                comm,
                f"RSS after column: per rank "
                f"{', '.join(f'{r:.2f}' for r in (per or []))} GiB, summed {tot:.2f} GiB",
            )
        t1 = time.perf_counter()
        e = fields.e_complex
        currents = {s.port_id: sheet_terminal_current(msh, tags_f, s, e, comm) for s in sheets}
        v_src = complex(dsheet.source_voltage_v)
        s_col = {}
        for s in sheets:
            v = complex(s.source_voltage_v) - currents[s.port_id] * z_p
            s_col[s.port_id] = (v - z_p * currents[s.port_id]) / v_src
        sigma_f = fields.sigma_field

        def _vol(measure):
            form = fem.form(0.5 * sigma_f * ufl.inner(e, e) * measure)
            return float(np.real(comm.allreduce(fem.assemble_scalar(form), op=MPI.SUM)))

        p_sheets = float(sum(_sheet_field_dissipation_w(acct, s, e, omega) for s in sheets))
        p_src = float(_source_power_w(acct, dsheet, e, omega))
        p_total = _vol(dx)
        p_phantom = _vol(dx(PHANTOM_CELL_TAG))
        p_surf = _surface_loss_w(msh, facets, z_s, e, comm)
        del fields
        t_post = time.perf_counter() - t1
        p_vol = p_total
        d = {
            "s": s_col,
            "t_solve": t_solve,
            "t_post": t_post,
            "p_src": p_src,
            "p_sheets": p_sheets,
            "p_phantom": p_phantom,
            "p_rest": p_total - p_phantom,
            "p_surf": p_surf,
        }
        d["rel"] = abs(p_src - p_sheets - p_vol - p_surf) / abs(p_src)
        d["blind"] = abs(p_src - p_sheets - p_vol) / abs(p_src)
        d["share"] = p_surf / (p_src - p_sheets)
        _say(
            comm,
            f"{driven_id}: post {t_post:.1f} s; S_driven {s_col[driven_id]:.6f}; sum|S|^2 "
            f"{sum(abs(v) ** 2 for v in s_col.values()):.6f} | P_src {p_src:.9e} sumP_sheet "
            f"{p_sheets:.9e} P_phantom {p_phantom:.9e} P_rest {d['p_rest']:.9e} P_surf "
            f"{p_surf:.9e} W; (P) rel {d['rel']:.3e}; NEG blind {d['blind']:.4e} = "
            f"{d['blind'] / DISCRETE_IDENTITY_RTOL:.1f}x band (PREDICTED ≈"
            f"{PREDICTED_BLIND_FACTOR:.0f}x); coil share {d['share']:.4f} (F-human "
            f"{PREDICTED_COIL_SHARE})",
        )
        return d

    col = _drive(DRIVEN, held=False)
    held = _drive(HELD, held=True)
    rss_mine, rss_per, rss_sum = _rss_gib(comm)

    # ---- (R1) one reciprocity pair ---------------------------------------
    s_34_33 = col["s"][HELD]
    s_33_34 = held["s"][DRIVEN]
    s_max = max(max(abs(v) for v in col["s"].values()), max(abs(v) for v in held["s"].values()))
    recip = abs(s_33_34 - s_34_33) / s_max

    # ---- printed price -----------------------------------------------------
    u1 = F_CLINICAL_R1_UNKNOWNS_RECORD[0]
    t_col = col["t_solve"]
    t_held = held["t_solve"]
    exp_meas = math.log(t_col / F_HUMAN_COL_S) / math.log(rho)
    c3 = t_build + t_col + 63 * t_held
    deg2_t = t_build + 287.1 * rho**COL_EXPONENT
    _say(
        comm,
        f"PRICE: t_b {t_build:.1f} s ({t_b_note}); t_col {t_col:.1f} s "
        f"(PREDICTED {F_HUMAN_COL_S * rho ** COL_EXPONENT:.1f}; conservative "
        f"{F_HUMAN_COL_S * (u1 / 500_568) ** 2:.1f}); t_held {t_held:.2f} s (PREDICTED "
        f"{F_HUMAN_HELD_S * rho ** (4 / 3):.2f}); U1 {u1}; measured exponent "
        f"ln(t_col/39.2)/ln rho = {exp_meas:.3f} (weekly {COL_EXPONENT}); RSS per rank "
        f"{', '.join(f'{r:.2f}' for r in (rss_per or []))} GiB, summed {rss_sum:.2f} GiB "
        f"(PREDICTED {PREDICTED_M_GIB})",
    )
    _say(
        comm,
        f"PROJECTIONS: C3 = t_b + t_col + 63 t_held = {c3:.0f} s -> "
        f"{'heavy' if c3 <= 590 else 'xl'}; degree-2 column at -n 16 t_b + 287.1 rho^0.84 = "
        f"{deg2_t:.0f} s, memory 89.31 rho^0.84 = {89.31 * rho ** COL_EXPONENT:.1f} GiB "
        f"(rho^(4/3): {89.31 * rho ** (4 / 3):.1f} GiB)",
    )
    if mode == "column":
        c3_cached = built["t_read"] + t_col + 63 * t_held
        _say(
            comm,
            f"PROJECTIONS (cached): t_read + t_col + 63 t_held = {c3_cached:.0f} s; "
            f"with item 469's -n 8 build {ATTEMPT1_T_BUILD_N8_S + t_col + 63 * t_held:.0f} s",
        )
    _say(
        comm,
        f"(P) {DRIVEN} rel {col['rel']:.3e}, {HELD} rel {held['rel']:.3e} (band "
        f"{DISCRETE_IDENTITY_RTOL:g}, ASSERTED); P_surf {col['p_surf']:.3e} W > 0 (ASSERTED); "
        f"NEG blind miss {col['blind']:.4e} > band (ASSERTED); (R1) |S33,34 - S34,33|/max|S| "
        f"{recip:.3e} (band {RECIPROCITY_BAND:.0e}, ASSERTED; PREDICTED ≈"
        f"{PREDICTED_RECIPROCITY:.0e})",
    )
    _say(comm, f"window {time.perf_counter() - t_window:.1f} s")
    assert col["rel"] <= DISCRETE_IDENTITY_RTOL, col["rel"]
    assert held["rel"] <= DISCRETE_IDENTITY_RTOL, held["rel"]
    assert col["p_surf"] > 0.0 and held["p_surf"] > 0.0, (col["p_surf"], held["p_surf"])
    assert col["blind"] > DISCRETE_IDENTITY_RTOL, col["blind"]
    assert recip <= RECIPROCITY_BAND, recip
