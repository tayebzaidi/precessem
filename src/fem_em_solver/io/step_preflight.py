"""`GEO-34` step 2 — the no-solve STEP preflight.

``write_step_preflight(step, config, out_dir, comm)`` imports a STEP file
through :func:`~fem_em_solver.io.step_import.import_step`, rebuilds the port
gap sheets declared by the config's ``port_sheets`` key dolfinx-side (the
solve path's own construction, :func:`~fem_em_solver.io.mesh._interface_facet_tags`,
merged into the import's facet tags so the outer wall and the sheets share one
facet grid), writes the combined XDMF / HDF5 (cell tags as the ``CellTags``
DG0 array, facet tags as their own grid) and the two-panel setup figure, and
prints its tag census.  Nothing is solved.

``port_sheets`` is a list of ``{"name": "P1", "tag": 211, "between":
["port_P1_lower", "port_P1_upper"]}``: the sheet is the interface between the
two named volume groups.  ``import_step`` ignores the key.  The file carries
the **full** gap sheet; any narrowing of a drive sheet is a solve-side choice.

Collective over ``comm``: the config is read on rank 0 and broadcast, and
every failure — including a setup figure that did not render
(``write_setup_figure`` returns ``None`` on a failed render) — is raised on
every rank.
"""

from __future__ import annotations

import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
from mpi4py import MPI

from fem_em_solver.io.mesh import _interface_facet_tags
from fem_em_solver.io.paraview_utils import write_xdmf_with_tags
from fem_em_solver.io.step_import import import_step, load_step_config
from fem_em_solver.post.setup_figure import MAX_TITLE_CHARS, write_setup_figure

FACET_TAGS_NAME = "facet_tags"
STEM = "preflight"


def port_sheet_interfaces(config) -> Dict[int, Tuple[int, int]]:
    """``{sheet tag: (lower cell tag, upper cell tag)}`` from ``port_sheets``."""
    groups = {g["name"]: int(g["tag"]) for g in config["volume_groups"]}
    out: Dict[int, Tuple[int, int]] = {}
    for sheet in config.get("port_sheets", []):
        lower, upper = sheet["between"]
        for name in (lower, upper):
            if name not in groups:
                raise ValueError(
                    f"step_preflight: port sheet {sheet['name']!r} names group "
                    f"{name!r}, which is not a config volume group"
                )
        out[int(sheet["tag"])] = (groups[lower], groups[upper])
    return out


def build_preflight_mesh(step, config, comm: MPI.Intracomm = MPI.COMM_WORLD):
    """``(mesh, cell_tags, facet_tags_with_sheets, diag, config)``; collective."""
    cfg = comm.bcast(load_step_config(config) if comm.rank == 0 else None, root=0)
    msh, cell_tags, facet_tags, diag = import_step(step, cfg, comm=comm)
    # The wall is import_step's own outer_boundary_tag group (the model's
    # exterior skin).  The dolfinx-side re-tagging fallback that stood here
    # (2026-09-30, while the import produced no wall group) is removed: a
    # silent re-tagging branch would hide the next regression.
    msh.topology.create_connectivity(msh.topology.dim - 1, msh.topology.dim)
    merged = _interface_facet_tags(
        msh, cell_tags, port_sheet_interfaces(cfg), existing=facet_tags
    )
    merged.name = FACET_TAGS_NAME
    return msh, cell_tags, merged, diag, cfg


def write_preflight_files(
    msh, cell_tags, facet_tags, out_dir, comm, *, config=None, title="", figure=True
):
    """Write ``preflight.xdmf`` / ``.h5`` (+ the setup PNG); collective.

    Returns ``{"xdmf", "h5", "png", "t_write_s", "t_figure_s"}`` (paths as
    strings, identical on every rank).
    """
    out_dir = Path(out_dir)
    if comm.rank == 0:
        out_dir.mkdir(parents=True, exist_ok=True)
    comm.barrier()
    base = out_dir / STEM
    t0 = time.perf_counter()
    write_xdmf_with_tags(base, msh, cell_tags, {}, comm=comm, facet_tags=facet_tags)
    t_write = time.perf_counter() - t0

    png = None
    t_fig = 0.0
    if figure:
        names = {int(g["tag"]): g["name"] for g in config["volume_groups"]}
        hide = tuple(t for t, n in names.items() if "air" in n.lower())
        translucent = tuple(t for t, n in names.items() if "phantom" in n.lower())
        t0 = time.perf_counter()
        written = write_setup_figure(
            msh, cell_tags, out_dir / f"{STEM}_setup.png",
            region_names=names, title=title, hide_tags=hide,
            translucent_tags=translucent, enabled=True, comm=comm,
        )
        t_fig = time.perf_counter() - t0
        ok = comm.bcast(written is not None if comm.rank == 0 else None, root=0)
        if not ok:
            raise RuntimeError(
                "step_preflight: write_setup_figure returned None (the render "
                "failed; see the [setup-figure] line above)"
            )
        png = str(out_dir / f"{STEM}_setup.png")
    return {
        "xdmf": str(base.with_suffix(".xdmf")),
        "h5": str(base.with_suffix(".h5")),
        "png": png,
        "t_write_s": t_write,
        "t_figure_s": t_fig,
    }


def _owned_count(values, indices, n_owned, tag) -> int:
    return int(np.count_nonzero((np.asarray(values) == tag) & (np.asarray(indices) < n_owned)))


def memory_counts(msh, cell_tags, facet_tags, cell_tag_list, facet_tag_list, comm):
    """Owned cells per cell tag and owned facets per facet tag, reduced."""
    tdim = msh.topology.dim
    n_cells = msh.topology.index_map(tdim).size_local
    n_facets = msh.topology.index_map(tdim - 1).size_local
    cells = {
        t: comm.allreduce(_owned_count(cell_tags.values, cell_tags.indices, n_cells, t))
        for t in cell_tag_list
    }
    facets = {
        t: comm.allreduce(_owned_count(facet_tags.values, facet_tags.indices, n_facets, t))
        for t in facet_tag_list
    }
    return cells, facets


def write_step_preflight(
    step, config, out_dir, comm: MPI.Intracomm = MPI.COMM_WORLD, *, title: Optional[str] = None
):
    """Import, tag, write and render; no solve.  Collective; returns a record."""
    t0 = time.perf_counter()
    msh, cell_tags, facets, diag, cfg = build_preflight_mesh(step, config, comm)
    t_import = time.perf_counter() - t0
    n_cells = int(msh.topology.index_map(msh.topology.dim).size_global)
    sheets = port_sheet_interfaces(cfg)
    if title is None:
        title = (
            f"GEO-34 preflight: {Path(step).name}, {n_cells} cells, "
            f"{len(cfg['volume_groups'])} volume groups, {len(sheets)} port sheets"
        )
    if len(title) > MAX_TITLE_CHARS:  # the default must never trip the guard
        title = f"GEO-34 preflight: {n_cells} cells, {len(sheets)} port sheets"
    files = write_preflight_files(
        msh, cell_tags, facets, out_dir, comm, config=cfg, title=title
    )
    cell_list = [int(g["tag"]) for g in cfg["volume_groups"]]
    facet_list = ([int(cfg["outer_boundary_tag"])] if cfg.get("outer_boundary_tag") is not None else []) + sorted(sheets)
    cells, fcounts = memory_counts(msh, cell_tags, facets, cell_list, facet_list, comm)
    if comm.rank == 0:
        print(
            f"[preflight] {Path(step).name}: cells {n_cells}; import+tag "
            f"{t_import:.2f} s (import {diag['import_wall_time_s']:.2f} s, mesh "
            f"{diag['mesh_wall_time_s']:.2f} s); write {files['t_write_s']:.2f} s; "
            f"figure {files['t_figure_s']:.2f} s at -n {comm.size}"
            f"\n[preflight] cells per tag {cells}"
            f"\n[preflight] facets per tag {fcounts}"
            f"\n[preflight] wrote {files['xdmf']}, {files['h5']}, {files['png']}",
            flush=True,
        )
    record = dict(files)
    record.update(
        {"cells": n_cells, "cells_per_tag": cells, "facets_per_tag": fcounts,
         "t_import_s": t_import}
    )
    return record
