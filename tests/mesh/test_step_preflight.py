"""`GEO-34` step 2 — the no-solve STEP preflight carries the solve path's census.

Fixture: step 1a's export (the F-small four-port rung through gmsh's OCC STEP
writer, `tests/mesh/test_step_import.py`'s ``round_trip``) into a broadcast
tempdir; config = `geo34_fsmall_birdcage_step_congruent.json` plus
``port_sheets`` (the four gap sheets 211–214 between each port's lower/upper
box), in a temp JSON.  ``scripts/run_step_case.py STEP CONFIG --preflight
--out DIR`` is called in-process (``main([...])``, collective).

Census C = owned cells and volume per cell tag, facets per facet tag, area per
sheet.  **C_file** is read on rank 0 from the written ``.xdmf`` / ``.h5`` with
``xml``, ``h5py`` and numpy only (tet determinants, triangle cross products;
every ``DataItem`` followed from the light data, no hard-coded HDF5 path) and
broadcast.  **C_solve** is in memory: a second ``import_step`` on the same STEP
and config plus the solve path's sheet rebuild
(`tests/validation/test_port_birdcage_four_port.py:349–352`), reduced by the
gates' helpers.

Asserted: (i) C_file = C_solve — counts exactly on the 11 cell tags, the wall
tag and 211–214; volumes, their sum and sheet areas ≤ ``CENSUS_REL_TOL``;
(ii) each sheet's file area = ``SHEET_AREA_M2`` within ``SHEET_AREA_TOL``;
(iii) PNG, ``.xdmf`` and ``.h5`` exist.  Negative controls: (n1) the preflight
from 1a's (unbound) config is unequal to the congruent C_solve; (n2, binary)
written without ``facet_tags=``, the file has no facet grid.

Standard tier, real build, ``-n 2``; windows in PROJECT_PLAN §9 item 118.
"""

from __future__ import annotations

import importlib.util
import json
import os
import tempfile
import time
import xml.etree.ElementTree as ET
from pathlib import Path

import h5py
import numpy as np
import pytest
from mpi4py import MPI

from fem_em_solver.io.mesh import _interface_facet_tags
from fem_em_solver.io.step_import import import_step, load_step_config
from fem_em_solver.io.step_preflight import (
    FACET_TAGS_NAME,
    port_sheet_interfaces,
    write_preflight_files,
)
from tests.mesh.helpers import tag_volume
from tests.mesh.test_birdcage_phantom_resolution import (
    DEFAULT_CELL_COUNT,
    SHEET_AREA_M2,
    SHEET_AREA_TOL,
    _global_tag_cell_count,
)
from tests.mesh.test_birdcage_port_sheets import _build
from tests.mesh.test_birdcage_port_terminals import (
    _global_facet_count,
    _interface_area_or_zero,
)
from tests.mesh.test_step_import import wall_census

REPO = Path(__file__).resolve().parents[2]
DATA = Path(__file__).parent / "data"
CONGRUENT_CONFIG = DATA / "geo34_fsmall_birdcage_step_congruent.json"
UNBOUND_CONFIG = DATA / "geo34_fsmall_birdcage_step.json"

# The import re-meshes identically (116 662 cells in five logs), but no
# reading has been compared past ten digits, hence 1e-9 (predicted ≤ 1e-13;
# printed).  Pre-registered in §9 item 118.
CENSUS_REL_TOL = 1.0e-9
# Printed beside the file's cell count, never gated: step 1b′'s congruent
# import (`20260926T093715Z_GEO-34-step1bprime.log`).
CONGRUENT_CELL_RECORD = 116685

PORT_SHEETS = [
    {"name": f"P{i}", "tag": 210 + i, "between": [f"port_P{i}_lower", f"port_P{i}_upper"]}
    for i in range(1, 5)
]


def _load_runner():
    spec = importlib.util.spec_from_file_location(
        "run_step_case", REPO / "scripts" / "run_step_case.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# --- C_file: xml + h5py + numpy, rank 0 --------------------------------------


def _data(item, base):
    """Array of one ``DataItem`` (HDF or XML), shaped by ``Dimensions``."""
    dims = [int(d) for d in item.get("Dimensions").split()]
    fmt = item.get("Format", "XML")
    if fmt == "HDF":
        h5_name, h5_path = item.text.strip().split(":")
        with h5py.File(base / h5_name, "r") as h5:
            arr = h5[h5_path][...]
    elif fmt == "XML":
        arr = np.array(item.text.split(), dtype=float)
    else:
        raise ValueError(f"unhandled DataItem Format {fmt!r}")
    return np.asarray(arr).reshape(dims)


def read_file_census(xdmf_path):
    """Census from the written file: rank-0 only, no dolfinx."""
    xdmf_path = Path(xdmf_path)
    base = xdmf_path.parent
    domain = ET.parse(xdmf_path).getroot().find("Domain")
    grids = domain.findall("Grid")
    tet = [g for g in grids if g.find("Topology") is not None
           and g.find("Topology").get("TopologyType") == "Tetrahedron"]
    tri = [g for g in grids if g.find("Topology") is not None
           and g.find("Topology").get("TopologyType") == "Triangle"]
    assert len(tet) == 1, f"{len(tet)} tetrahedral grids"
    mesh_grid = tet[0]
    x = _data(mesh_grid.find("Geometry").find("DataItem"), base)
    topo = _data(mesh_grid.find("Topology").find("DataItem"), base).astype(np.int64)
    attr = [a for a in mesh_grid.findall("Attribute") if a.get("Name") == "CellTags"]
    assert len(attr) == 1 and attr[0].get("Center") == "Cell"
    ctag = np.rint(_data(attr[0].find("DataItem"), base).reshape(-1)).astype(np.int64)
    assert ctag.size == topo.shape[0], (ctag.size, topo.shape)
    p = x[topo]
    vol = np.abs(np.einsum(
        "ij,ij->i", np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]), p[:, 3] - p[:, 0]
    )) / 6.0
    census = {
        "n_grids": len(grids),
        "n_facet_grids": len(tri),
        "cells_total": int(topo.shape[0]),
        "cells": {int(t): int(np.count_nonzero(ctag == t)) for t in np.unique(ctag)},
        "volume": {int(t): float(vol[ctag == t].sum()) for t in np.unique(ctag)},
        "facets": {},
        "area": {},
    }
    if len(tri) == 1:
        fgrid = tri[0]
        geom = fgrid.find("Geometry")
        fx = x if geom is None else _data(geom.find("DataItem"), base)
        ftopo = _data(fgrid.find("Topology").find("DataItem"), base).astype(np.int64)
        fattr = fgrid.findall("Attribute")
        assert len(fattr) == 1, [a.get("Name") for a in fattr]
        census["facet_attr_name"] = fattr[0].get("Name")
        fval = np.rint(_data(fattr[0].find("DataItem"), base).reshape(-1)).astype(np.int64)
        q = fx[ftopo]
        area = 0.5 * np.linalg.norm(np.cross(q[:, 1] - q[:, 0], q[:, 2] - q[:, 0]), axis=1)
        for t in np.unique(fval):
            census["facets"][int(t)] = int(np.count_nonzero(fval == t))
            census["area"][int(t)] = float(area[fval == t].sum())
    return census


def _bcast_file_census(xdmf, comm):
    payload = None
    if comm.rank == 0:
        try:
            payload = ("ok", read_file_census(xdmf))
        except Exception as exc:  # noqa: BLE001 — re-raised on every rank
            payload = ("err", repr(exc))
    status, value = comm.bcast(payload, root=0)
    if status != "ok":
        raise RuntimeError(f"file census failed on rank 0: {value}")
    return value


# --- fixtures ----------------------------------------------------------------


def _write_config(src, path):
    config = load_step_config(src)
    config["port_sheets"] = PORT_SHEETS
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(config, handle, indent=1)


@pytest.fixture(scope="module")
def case():
    """1a's STEP export and both temp configs, in one broadcast tempdir."""
    comm = MPI.COMM_WORLD
    root = tempfile.mkdtemp(prefix="geo34_pf_") if comm.rank == 0 else None
    root = Path(comm.bcast(root, root=0))
    step = root / "fsmall_birdcage.step"
    t0 = time.perf_counter()
    gen = _build(True, step_export_path=str(step))
    t_export = time.perf_counter() - t0
    del gen
    if comm.rank == 0:
        _write_config(CONGRUENT_CONFIG, root / "congruent_ports.json")
        _write_config(UNBOUND_CONFIG, root / "unbound_ports.json")
        print(f"\n[GEO-34 2] generator rung + STEP export {t_export:.2f} s", flush=True)
    comm.barrier()
    return {"root": root, "step": step,
            "congruent": root / "congruent_ports.json",
            "unbound": root / "unbound_ports.json"}


def _run_preflight(case, config_key, out_name):
    comm = MPI.COMM_WORLD
    out = case["root"] / out_name
    t0 = time.perf_counter()
    rc = _load_runner().main(
        [str(case["step"]), str(case[config_key]), "--preflight", "--out", str(out)]
    )
    elapsed = time.perf_counter() - t0
    assert rc == 0
    census = _bcast_file_census(out / "preflight.xdmf", comm)
    if comm.rank == 0:
        print(f"[GEO-34 2] preflight ({config_key}) {elapsed:.2f} s; file cells "
              f"{census['cells_total']} (congruent record {CONGRUENT_CELL_RECORD}, "
              f"primitive {DEFAULT_CELL_COUNT}; printed, not gated)", flush=True)
    return {"out": out, "census": census}


@pytest.fixture(scope="module")
def preflight(case):
    return _run_preflight(case, "congruent", "congruent")


@pytest.fixture(scope="module")
def solve_census(case):
    """C_solve: a second import + the solve path's sheet rebuild, in memory."""
    comm = MPI.COMM_WORLD
    config = load_step_config(case["congruent"])
    t0 = time.perf_counter()
    msh, cell_tags, facet_tags, _diag = import_step(case["step"], config, comm=comm)
    t_import = time.perf_counter() - t0
    tdim = msh.topology.dim
    msh.topology.create_connectivity(tdim - 1, tdim)
    msh.topology.create_entity_permutations()
    sheets = port_sheet_interfaces(config)
    tags_f = _interface_facet_tags(msh, cell_tags, sheets)
    cell_list = [int(g["tag"]) for g in config["volume_groups"]]
    wall = int(config["outer_boundary_tag"])
    # The wall is counted independently of the import — owned facets with one
    # attached cell (the import partitions with GhostMode.shared_facet, so a
    # partition-boundary facet sees both cells) — and the import's own
    # outer_boundary_tag count is recorded beside it (asserted equal in
    # test_census; known-issues 2026-09-30 fix, §9 item 119).
    import_wall, n_wall, wall_not_one_cell = wall_census(msh, facet_tags, wall, comm)
    census = {
        "cells_total": int(msh.topology.index_map(tdim).size_global),
        "cells": {t: _global_tag_cell_count(msh, cell_tags, t, comm) for t in cell_list},
        "volume": {t: tag_volume(msh, cell_tags, t, comm) for t in cell_list},
        "facets": {wall: n_wall},
        "area": {},
        "import_wall_facets": import_wall,
        "wall_not_one_cell": wall_not_one_cell,
        "facet_tags_is_none": facet_tags is None,
        "bbox": _diag["model_bbox"],
    }
    for t in sorted(sheets):
        census["facets"][t] = _global_facet_count(msh, tags_f, t, comm)
        census["area"][t] = _interface_area_or_zero(msh, tags_f, t, comm)
    t_total = time.perf_counter() - t0
    if comm.rank == 0:
        print(f"[GEO-34 2] C_solve: import {t_import:.2f} s, census "
              f"{t_total - t_import:.2f} s; cells {census['cells_total']}"
              f"\n[GEO-34 2] C_solve cells per tag {census['cells']}"
              f"\n[GEO-34 2] C_solve facets per tag {census['facets']} (wall = owned "
              f"exterior facets; import's own wall-tag facets {import_wall}, of which "
              f"not on exactly one cell {wall_not_one_cell}; facet_tags is None: "
              f"{facet_tags is None})",
              flush=True)
    return {"census": census, "mesh": msh, "cell_tags": cell_tags, "config": config}


def _compare(c_file, c_solve):
    """(count mismatches, worst volume/area relative difference, detail)."""
    mism = []
    if c_file["cells_total"] != c_solve["cells_total"]:
        mism.append(("cells_total", c_file["cells_total"], c_solve["cells_total"]))
    for t, n in c_solve["cells"].items():
        if c_file["cells"].get(t, 0) != n:
            mism.append((f"cell {t}", c_file["cells"].get(t, 0), n))
    for t, n in c_solve["facets"].items():
        if c_file["facets"].get(t, 0) != n:
            mism.append((f"facet {t}", c_file["facets"].get(t, 0), n))
    rel = {}
    for t, v in c_solve["volume"].items():
        rel[f"vol {t}"] = abs(c_file["volume"].get(t, 0.0) / v - 1.0)
    vs = sum(c_solve["volume"].values())
    vf = sum(c_file["volume"].get(t, 0.0) for t in c_solve["volume"])
    rel["vol sum"] = abs(vf / vs - 1.0)
    for t, a in c_solve["area"].items():
        rel[f"area {t}"] = abs(c_file["area"].get(t, 0.0) / a - 1.0)
    return mism, max(rel.values()), rel


# --- anchors -----------------------------------------------------------------


def test_census(preflight, solve_census):
    """(i) C_file = C_solve: counts exactly, volumes/areas ≤ CENSUS_REL_TOL."""
    comm = MPI.COMM_WORLD
    c_file, c_solve = preflight["census"], solve_census["census"]
    mism, worst, rel = _compare(c_file, c_solve)
    if comm.rank == 0:
        print(f"\n[GEO-34 2] (i) file cells per tag {c_file['cells']}"
              f"\n[GEO-34 2] (i) file facets per tag {c_file['facets']} "
              f"(facet array {c_file.get('facet_attr_name')!r})"
              f"\n[GEO-34 2] (i) count mismatches {mism}"
              f"\n[GEO-34 2] (i) worst volume/area rel diff {worst:.3e} (tol "
              f"{CENSUS_REL_TOL:.0e}); " + ", ".join(f"{k} {v:.2e}" for k, v in rel.items()),
              flush=True)
    # (A1) the import's own wall group is the exterior, facet for facet.
    wall = int(solve_census["config"]["outer_boundary_tag"])
    assert not c_solve["facet_tags_is_none"]
    assert c_solve["import_wall_facets"] == c_solve["facets"][wall], (
        c_solve["import_wall_facets"], c_solve["facets"][wall])
    assert c_solve["wall_not_one_cell"] == 0, c_solve["wall_not_one_cell"]
    assert c_file["n_facet_grids"] == 1
    assert c_file.get("facet_attr_name") == FACET_TAGS_NAME
    assert set(c_file["cells"]) == set(c_solve["cells"]), (c_file["cells"], c_solve["cells"])
    assert set(c_file["facets"]) == set(c_solve["facets"]), (c_file["facets"], c_solve["facets"])
    assert not mism, mism
    assert worst <= CENSUS_REL_TOL, rel


def test_sheet_areas(preflight, solve_census):
    """(ii) each sheet's area read from the file is `GEO-18`'s closed form."""
    comm = MPI.COMM_WORLD
    area = preflight["census"]["area"]
    rel = {t: area.get(t, 0.0) / SHEET_AREA_M2 - 1.0 for t in (211, 212, 213, 214)}
    if comm.rank == 0:
        for t, r in rel.items():
            print(f"[GEO-34 2] (ii) sheet {t}: file area {area.get(t, 0.0):.12e} m^2, "
                  f"/ SHEET_AREA_M2 - 1 = {r:.3e}; facets "
                  f"{preflight['census']['facets'].get(t)}", flush=True)
        wall = preflight["census"]["area"].get(1)
        x0, y0, z0, x1, y1, z1 = solve_census["census"]["bbox"]
        a, b, c = x1 - x0, y1 - y0, z1 - z0
        box = 2.0 * (a * b + b * c + c * a)
        print(f"[GEO-34 2] wall (tag 1) file area {wall} m^2 vs model-bbox surface "
              f"{box:.12e} m^2 (rel {(wall or 0.0) / box - 1.0:.3e}; printed, not gated)",
              flush=True)
    for t, r in rel.items():
        assert abs(r) <= SHEET_AREA_TOL, (t, r)


def test_artifacts(preflight):
    """(iii) the PNG and both data files exist and are non-empty."""
    out = preflight["out"]
    for name in ("preflight.xdmf", "preflight.h5", "preflight_setup.png"):
        path = out / name
        assert path.is_file() and path.stat().st_size > 0, path


def test_no_facet_grid(case, solve_census):
    """(n2, binary) written without ``facet_tags=``, the file has no facet grid."""
    comm = MPI.COMM_WORLD
    out = case["root"] / "no_facets"
    write_preflight_files(solve_census["mesh"], solve_census["cell_tags"], None, out,
                          comm, figure=False)
    census = _bcast_file_census(out / "preflight.xdmf", comm)
    if comm.rank == 0:
        print(f"[GEO-34 2] (n2) grids {census['n_grids']}, facet grids "
              f"{census['n_facet_grids']}, facets {census['facets']}", flush=True)
    assert census["n_facet_grids"] == 0
    assert census["facets"] == {}
    assert census["cells_total"] == solve_census["census"]["cells_total"]


def test_mismatched_config(case, solve_census):
    """(n1, asserted) the preflight from 1a's unbound config ≠ congruent C_solve."""
    comm = MPI.COMM_WORLD
    unbound = _run_preflight(case, "unbound", "unbound")
    mism, worst, _rel = _compare(unbound["census"], solve_census["census"])
    if comm.rank == 0:
        c = unbound["census"]
        print(f"\n[GEO-34 2] (n1) unbound file cells {c['cells_total']} vs congruent "
              f"C_solve {solve_census['census']['cells_total']}; sheet facets "
              f"{[c['facets'].get(t) for t in (211, 212, 213, 214)]} vs "
              f"{[solve_census['census']['facets'][t] for t in (211, 212, 213, 214)]}"
              f"\n[GEO-34 2] (n1) count mismatches {len(mism)}: {mism}; worst "
              f"volume/area rel diff {worst:.3e}", flush=True)
    assert mism, "the unbound config's preflight matched the congruent census"
