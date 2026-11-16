"""`GEO-34` step 1a — the F-small birdcage through a STEP round trip.

The four-port rung is built exactly as
`tests/mesh/test_birdcage_port_sheets.py::_build(True)` builds it, with the
generator's additive ``step_export_path`` set, so gmsh's own OCC writer emits
the finished CAD model as STEP before meshing.  `io/step_import.import_step`
re-imports it, claims every volume for one group of the JSON config beside the
test (`tests/mesh/data/geo34_fsmall_birdcage_step.json`), sizes and meshes it.
The Gmsh OCC writer stands in for an independent CAD kernel (CadQuery is not in
the image): this exercises the STEP reader, entity renumbering and the group
map, not kernel independence (that is step 1c).

Gates (every band imported, none restated): the cell-tag census equals the
generator mesh's; the four port sheets rebuilt dolfinx-side have `GEO-18`'s
closed-form area ``SHEET_AREA_M2`` within ``SHEET_AREA_TOL``; meshed / CAD
conductor volume ≥ ``CAD_MASS_GATE``.  Negative controls: the config with P1's
lower port-box group deleted, and a group whose selector matches nothing, both
raise ``ValueError`` naming the offender.  Printed, never gated: cell count vs
the primitive rung's 116 085, per-group CAD volume generator vs import, wall
times.  Mesh-level only — no port model, no solve.
"""

from __future__ import annotations

import copy
import os
import tempfile
import time
from pathlib import Path

import numpy as np
import pytest
import ufl
from dolfinx import default_scalar_type, fem
from mpi4py import MPI

from fem_em_solver.io.mesh import _interface_facet_tags
from fem_em_solver.io.step_import import import_step, load_step_config
from tests.mesh.helpers import global_cell_tag_set
from tests.mesh.test_birdcage_conductor_sizing import CAD_MASS_GATE
from tests.mesh.test_birdcage_phantom_resolution import (
    DEFAULT_CELL_COUNT,
    SHEET_AREA_M2,
    SHEET_AREA_TOL,
)
from tests.mesh.test_birdcage_port_sheet_prerequisite import CONDUCTOR_RESOLUTION
from tests.mesh.test_birdcage_port_sheets import (
    PORT_LOWER,
    PORT_UPPER,
    SHEET_IFACE,
    _build,
)
from tests.mesh.test_birdcage_port_tags import LEG_COUNT, RESOLUTION
from tests.mesh.test_birdcage_port_terminals import _interface_area_or_zero
from tests.mesh.test_coil_phantom_conforming import _tag_volume

CONFIG_PATH = Path(__file__).parent / "data" / "geo34_fsmall_birdcage_step.json"

# The F-small air box is 0.24 x 0.24 x 0.20 m (the import's volume 2 bbox,
# `20260930T201605Z_GEO-34-step2-W1.log:1850`), so its six walls have area
# 2(ab + bc + ca).  Known-issues 2026-09-30 fix, PROJECT_PLAN §9 item 119.
AIR_BOX_ABC_M = (0.24, 0.24, 0.20)
WALL_AREA_M2 = 2.0 * (
    AIR_BOX_ABC_M[0] * AIR_BOX_ABC_M[1]
    + AIR_BOX_ABC_M[1] * AIR_BOX_ABC_M[2]
    + AIR_BOX_ABC_M[2] * AIR_BOX_ABC_M[0]
)
WALL_AREA_REL_TOL = 1.0e-12


def wall_census(mesh, facet_tags, tag, comm):
    """``(import's owned tag facets, owned exterior facets, tag facets not on
    exactly one cell)``, each allreduced.

    "Exterior" is counted independently of the import: owned facets with one
    attached cell, which is the exterior only under ``GhostMode.shared_facet``
    (``import_step`` partitions with it), where a partition-boundary facet
    sees both of its cells.
    """
    tdim = mesh.topology.dim
    mesh.topology.create_connectivity(tdim - 1, tdim)
    conn = mesh.topology.connectivity(tdim - 1, tdim)
    ncell = np.diff(conn.offsets)
    n_owned = mesh.topology.index_map(tdim - 1).size_local
    n_ext = int(np.count_nonzero(ncell[:n_owned] == 1))
    if facet_tags is None:
        n_tag, n_bad = 0, 0
    else:
        idx = np.asarray(facet_tags.indices)
        sel = idx[(np.asarray(facet_tags.values) == tag) & (idx < n_owned)]
        n_tag = int(sel.size)
        n_bad = int(np.count_nonzero(ncell[sel] != 1))
    return tuple(int(comm.allreduce(v, op=MPI.SUM)) for v in (n_tag, n_ext, n_bad))


def old_bbox_wall_predicate(b, model_bbox):
    """The wall predicate as it stood at ``4857675d``
    (``src/fem_em_solver/io/step_import.py:198–208``), copied verbatim — the
    negative control's pinned source, never ``HEAD``."""
    x0, y0, z0, x1, y1, z1 = model_bbox
    tol = 1.0e-9 * max(x1 - x0, y1 - y0, z1 - z0)
    return (
        abs(b[0] - b[3]) < tol and (abs(b[0] - x0) < tol or abs(b[0] - x1) < tol)
        or abs(b[1] - b[4]) < tol and (abs(b[1] - y0) < tol or abs(b[1] - y1) < tol)
        or abs(b[2] - b[5]) < tol and (abs(b[2] - z0) < tol or abs(b[2] - z1) < tol)
    )


@pytest.fixture(scope="module")
def round_trip():
    """Generator mesh (+ STEP export) and the imported mesh, built once."""
    comm = MPI.COMM_WORLD
    step_dir = tempfile.mkdtemp(prefix="geo34_") if comm.rank == 0 else None
    step_dir = comm.bcast(step_dir, root=0)
    step_path = os.path.join(step_dir, "fsmall_birdcage.step")

    gen = _build(True, step_export_path=step_path)
    started = time.perf_counter()
    imported = import_step(step_path, CONFIG_PATH, comm=comm)
    import_elapsed = time.perf_counter() - started
    return {"step_path": step_path, "gen": gen, "imp": imported,
            "import_elapsed_s": import_elapsed}


def test_config_sizes_are_the_gates():
    """The JSON's sizes are the gate's constants, not restatements that drift."""
    config = load_step_config(CONFIG_PATH)
    assert config["resolution"] == RESOLUTION
    conductor = [g for g in config["volume_groups"] if g["name"] == "conductor"][0]
    assert conductor["grading"]["size"] == CONDUCTOR_RESOLUTION


def test_step_round_trip_reproduces_the_fixture(round_trip):
    comm = MPI.COMM_WORLD
    ports = list(range(1, LEG_COUNT + 1))
    g_mesh, g_cells, _, g_diag, g_elapsed = round_trip["gen"]
    mesh, cells, _, diag = round_trip["imp"]

    g_tags = global_cell_tag_set(g_mesh, g_cells)
    tags = global_cell_tag_set(mesh, cells)
    n_cells = mesh.topology.index_map(3).size_global
    n_cells_gen = g_mesh.topology.index_map(3).size_global

    halves = {i: (PORT_LOWER + i, PORT_UPPER + i) for i in ports}
    sheet_tags = _interface_facet_tags(
        mesh, cells, {SHEET_IFACE + i: halves[i] for i in ports}
    )
    sheet_area = {
        i: _interface_area_or_zero(mesh, sheet_tags, SHEET_IFACE + i, comm)
        for i in ports
    }
    conductor_meshed = _tag_volume(mesh, cells, 1, comm)
    conductor_cad = diag["cad_mass_by_group"]["conductor"]

    if comm.rank == 0:
        print(
            f"\n[GEO-34 1a] generator: cells={n_cells_gen} (record "
            f"{DEFAULT_CELL_COUNT}) rung={g_elapsed:.2f} s "
            f"mesh={g_diag['mesh_wall_time_s']:.2f} s"
            f"\n[GEO-34 1a] import: {diag['n_volumes']} volumes, cells={n_cells} "
            f"(primitive rung {DEFAULT_CELL_COUNT}, ratio "
            f"{n_cells / DEFAULT_CELL_COUNT:.6f}; printed, not gated) "
            f"import={diag['import_wall_time_s']:.2f} s "
            f"mesh={diag['mesh_wall_time_s']:.2f} s "
            f"total={round_trip['import_elapsed_s']:.2f} s"
            f"\n[GEO-34 1a] model bbox={diag['model_bbox']}"
            f"\n[GEO-34 1a] cell tags generator={sorted(g_tags)} import={sorted(tags)}",
            flush=True,
        )
        for name, mass in diag["cad_mass_by_group"].items():
            ref = g_diag["cad_mass_by_group"].get(name)
            rel = float("nan") if ref is None else mass / ref - 1.0
            print(
                f"[GEO-34 1a] CAD volume {name}: import={mass:.12e} "
                f"generator={ref if ref is None else f'{ref:.12e}'} rel diff={rel:.3e} "
                f"({diag['volumes_by_group'][name]} volume(s))",
                flush=True,
            )
        for i in ports:
            print(
                f"[GEO-34 1a] P{i} sheet area={sheet_area[i]:.12e} m^2 "
                f"/ SHEET_AREA_M2 - 1 = {sheet_area[i] / SHEET_AREA_M2 - 1.0:.3e} "
                f"(tol {SHEET_AREA_TOL})",
                flush=True,
            )
        print(
            f"[GEO-34 1a] conductor meshed={conductor_meshed:.9e} "
            f"CAD={conductor_cad:.9e} ratio={conductor_meshed / conductor_cad:.6f} "
            f"(gate >= {CAD_MASS_GATE})",
            flush=True,
        )

    # (i) tag census: the generator's set, and every group claimed >= 1 volume.
    assert tags == g_tags, (
        f"imported cell tags {sorted(tags)} != generator's {sorted(g_tags)}"
    )
    assert all(n >= 1 for n in diag["volumes_by_group"].values())
    # (ii) the four sheets at GEO-18's closed form.
    for i in ports:
        assert abs(sheet_area[i] / SHEET_AREA_M2 - 1.0) < SHEET_AREA_TOL, (
            f"P{i}: imported sheet area {sheet_area[i]:.12e} m^2 vs the closed "
            f"form {SHEET_AREA_M2:.12e} m^2"
        )
    # (iii) the CAD-mass gate on the imported conductor.
    assert conductor_meshed / conductor_cad >= CAD_MASS_GATE, (
        f"imported conductor meshed/CAD={conductor_meshed / conductor_cad:.6f} "
        f"< {CAD_MASS_GATE}"
    )


def test_unmapped_port_group_is_refused(round_trip):
    """Negative control: P1's lower port-box group deleted ⇒ ValueError."""
    config = load_step_config(CONFIG_PATH)
    config["volume_groups"] = [
        g for g in config["volume_groups"] if g["name"] != "port_P1_lower"
    ]
    with pytest.raises(ValueError, match=r"claimed by no config group.*centroid=\(") as exc:
        import_step(round_trip["step_path"], config, comm=MPI.COMM_WORLD)
    if MPI.COMM_WORLD.rank == 0:
        print(f"\n[GEO-34 1a] control (P1 lower deleted): {exc.value}", flush=True)


def test_selector_matching_nothing_is_refused(round_trip):
    """Negative control: a group whose selector matches nothing ⇒ ValueError."""
    config = load_step_config(CONFIG_PATH)
    config["volume_groups"] = copy.deepcopy(config["volume_groups"]) + [
        {"name": "port_P9_nowhere", "tag": 109,
         "selectors": [{"centroid_in_box": [5.0, 5.0, 5.0, 6.0, 6.0, 6.0]}]}
    ]
    with pytest.raises(ValueError, match="matched no volume: port_P9_nowhere") as exc:
        import_step(round_trip["step_path"], config, comm=MPI.COMM_WORLD)
    if MPI.COMM_WORLD.rank == 0:
        print(f"\n[GEO-34 1a] control (selector matches nothing): {exc.value}", flush=True)


def test_outer_boundary_is_the_exterior(round_trip):
    """The import's wall group is the mesh's exterior (known-issues 2026-09-30).

    (A1) owned tag-1 facets = owned exterior facets, each on exactly one cell;
    (A2) ∫ds(1) = 2(ab + bc + ca) of the air box to ``WALL_AREA_REL_TOL``;
    control (binary): the pinned old bbox-plane predicate accepts 0 walls.
    """
    comm = MPI.COMM_WORLD
    config = load_step_config(CONFIG_PATH)
    wall = int(config["outer_boundary_tag"])
    mesh, _cells, facets, diag = round_trip["imp"]
    assert facets is not None, "import_step returned facet_tags=None"

    n_tag, n_ext, n_bad = wall_census(mesh, facets, wall, comm)
    one = fem.Constant(mesh, default_scalar_type(1.0))
    ds = ufl.Measure("ds", domain=mesh, subdomain_data=facets)
    area_local = fem.assemble_scalar(fem.form(one * ds(wall)))
    area = float(np.real(comm.allreduce(area_local, op=MPI.SUM)))
    rel = area / WALL_AREA_M2 - 1.0

    surfaces = diag["outer_boundary_surfaces"]
    bboxes = diag["outer_boundary_bboxes"]
    model_bbox = diag["model_bbox"]
    old_accepts = [s for s, b in zip(surfaces, bboxes)
                   if old_bbox_wall_predicate(b, model_bbox)]
    x0, y0, z0, x1, y1, z1 = model_bbox
    tol = 1.0e-9 * max(x1 - x0, y1 - y0, z1 - z0)
    if comm.rank == 0:
        cad = diag["outer_boundary_cad_area"]
        print(
            f"\n[GEO-34 wall] cells={mesh.topology.index_map(3).size_global}; "
            f"wall surfaces {len(surfaces)} {surfaces}"
            f"\n[GEO-34 wall] CAD wall area {cad:.12e} m^2 vs 2(ab+bc+ca) "
            f"{WALL_AREA_M2:.12e} (rel {cad / WALL_AREA_M2 - 1.0:.3e}; printed)"
            f"\n[GEO-34 wall] (A1) owned tag-{wall} facets {n_tag}, owned exterior "
            f"facets {n_ext}, tag-{wall} facets not on exactly one cell {n_bad}"
            f"\n[GEO-34 wall] (A2) int ds({wall}) = {area:.15e} m^2 vs "
            f"{WALL_AREA_M2:.15e} (rel {rel:.3e}, tol {WALL_AREA_REL_TOL:.0e})"
            f"\n[GEO-34 wall] model bbox {model_bbox}; old predicate tol {tol:.3e}",
            flush=True,
        )
        lo, hi = model_bbox[:3], model_bbox[3:]
        for s, b in zip(surfaces, bboxes):
            k = int(np.argmin([abs(b[i] - b[i + 3]) for i in range(3)]))
            thick = abs(b[k] - b[k + 3])
            offset = min(abs(b[k] - lo[k]), abs(b[k] - hi[k]))
            print(
                f"[GEO-34 wall] surface {s}: axis {'xyz'[k]}, bbox thickness "
                f"{thick:.3e}, offset from model plane {offset:.3e} (old tol "
                f"{tol:.3e}); old predicate accepts: "
                f"{old_bbox_wall_predicate(b, model_bbox)}",
                flush=True,
            )
        print(f"[GEO-34 wall] control: old predicate accepts {len(old_accepts)} "
              f"wall(s) {old_accepts}; new finds {len(surfaces)}", flush=True)

    assert n_tag == n_ext, (n_tag, n_ext)
    assert n_bad == 0, n_bad
    assert abs(rel) <= WALL_AREA_REL_TOL, (area, WALL_AREA_M2, rel)
    assert len(old_accepts) == 0, old_accepts
    assert len(surfaces) >= 1
