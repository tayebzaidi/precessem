"""`GEO-37` step 2 (§9 item 267, chain U11): UFL ``dS`` over the native
interior facet tag on the 0.11 image.

Probe script (never collected: ``python_files = test_*.py``).  Item 222's
fixture -- unit box, embedded rectangle 0.8 x 0.6 at z = 0.47 tagged 201, box
faces tagged 101, h = 0.08 -- is built by importing item 222's
``_build_model`` verbatim from ``tests/mesh/probe_geo37_interior_facet.py``
and read through the project shim ``_model_to_mesh`` with an explicit cell
partitioner of ghost mode ``--ghost {shared_facet,none}``.  The integrals are
then *assembled* (the operation the 2026-08 `PORT-1` 3b-iv hang was localised
to, which item 222 deliberately avoided):

(I1) ``1 dS(201)`` = 0.48; (I2) ``g('+') dS(201)``, g = 2 + 3x + 5y, = 2.88;
(I3) ``(u('+') + u('-')) dS(201)`` with DG0 u = 1 below / 3 above = 1.92
(``u('+') dS(201)`` printed only); (I4) ``1 ds(101)`` = 6; (c1)
``1 dS(999)`` (absent tag) ; (c2) I1 on a re-tagged 201 with one owned facet
dropped.  Every ``assemble_scalar`` is ``allreduce``d.  Rank 0 writes a JSON
record whose ``completed: true`` is written last; the asserting reader is
``tests/mesh/test_geo37_interior_dS.py``.
"""

import argparse
import importlib.util
import json
import os
import time

import dolfinx
import gmsh
import numpy as np
import ufl
from dolfinx import fem
from mpi4py import MPI

from fem_em_solver.io.mesh import (
    _model_to_mesh,
    _raise_geometry_failure_on_every_rank,
)

_spec = importlib.util.spec_from_file_location(
    "probe_geo37_interior_facet",
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "probe_geo37_interior_facet.py"),
)
_p222 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_p222)

H = _p222.H
RECT = _p222.RECT
TAG_RECT, TAG_BOX = _p222.TAG_RECT, _p222.TAG_BOX
TAG_ABSENT = 999
QD = {"quadrature_degree": 2}  # POST-5 step 1: pin, or FFCx can run > 9 min


def _assemble(comm, form_expr):
    """Collective: allreduced scalar, or the exception text (recorded)."""
    try:
        local = fem.assemble_scalar(fem.form(form_expr))
        return float(comm.allreduce(float(np.real(local)), op=MPI.SUM)), None
    except BaseException as e:  # noqa: BLE001 -- recorded, not swallowed silently
        return None, f"{type(e).__name__}: {e}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--ghost", choices=("shared_facet", "none"), required=True)
    args = ap.parse_args()

    comm = MPI.COMM_WORLD
    rank = comm.rank
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)

    t0 = time.perf_counter()
    err, n_gmsh = None, None
    if rank == 0:
        try:
            n_gmsh = _p222._build_model()
        except BaseException as e:  # noqa: BLE001 -- GEO-23 pattern
            err = e
            gmsh.finalize()
    _raise_geometry_failure_on_every_rank(comm, 0, "GEO-37 step 2 probe", err, H)
    n_gmsh = comm.bcast(n_gmsh, root=0)
    t1 = time.perf_counter()

    mode = (
        dolfinx.mesh.GhostMode.shared_facet
        if args.ghost == "shared_facet"
        else dolfinx.mesh.GhostMode.none
    )
    # 0.11: max_facet_to_cell_links is required; 2 for a tet manifold mesh
    # (io/mesh.py OPS-18 step 3 comment).
    partitioner = dolfinx.mesh.create_cell_partitioner(mode, 2)
    mesh, cell_tags, facet_tags = _model_to_mesh(
        gmsh.model, comm, 0, gdim=3, partitioner=partitioner
    )
    gmsh.finalize()
    comm.Barrier()
    t2 = time.perf_counter()

    tdim = mesh.topology.dim
    fdim = tdim - 1
    imap_c = mesh.topology.index_map(tdim)
    n_cells = comm.allreduce(imap_c.size_local, op=MPI.SUM)
    n_ghost_cells = comm.allreduce(imap_c.num_ghosts, op=MPI.SUM)
    mesh.topology.create_connectivity(fdim, tdim)
    n_owned_f = mesh.topology.index_map(fdim).size_local

    # Globally reduced facet-tag set (for the (D5) validator call in the test).
    tag_set = sorted(
        {
            int(t)
            for part in comm.allgather(sorted({int(v) for v in facet_tags.values}))
            for t in part
        }
    )

    dS = ufl.Measure("dS", domain=mesh, subdomain_data=facet_tags, metadata=QD)
    ds = ufl.Measure("ds", domain=mesh, subdomain_data=facet_tags, metadata=QD)
    x = ufl.SpatialCoordinate(mesh)
    g = 2 + 3 * x[0] + 5 * x[1]

    V0 = fem.functionspace(mesh, ("DG", 0))
    u = fem.Function(V0)
    n_all_cells = imap_c.size_local + imap_c.num_ghosts
    mids = dolfinx.mesh.compute_midpoints(mesh, tdim, np.arange(n_all_cells, dtype=np.int32))
    u.x.array[:n_all_cells] = np.where(mids[:, 2] > RECT[2], 3.0, 1.0)
    u.x.scatter_forward()

    one = fem.Constant(mesh, dolfinx.default_scalar_type(1.0))
    rec = {
        "comm_size": comm.size,
        "ghost": args.ghost,
        "N_gmsh": n_gmsh,
        "cells": n_cells,
        "cells_ghost": n_ghost_cells,
        "h": H,
        "global_facet_tag_set": tag_set,
    }

    # Owned-facet count / local-copy count of 201 (item 222's arithmetic, context).
    local201 = facet_tags.find(TAG_RECT)
    owned201 = local201[local201 < n_owned_f]
    rec["count_201"] = int(comm.allreduce(len(owned201), op=MPI.SUM))
    rec["per_rank_owned_201"] = comm.gather(int(len(owned201)), root=0)
    rec["per_rank_local_201"] = comm.gather(int(len(local201)), root=0)
    # Facets of 201 that have only one local cell (the unghosted hazard).
    f2c = mesh.topology.connectivity(fdim, tdim)
    one_sided = sum(1 for f in local201 if len(f2c.links(int(f))) < 2)
    rec["per_rank_one_sided_201"] = comm.gather(int(one_sided), root=0)

    t3 = time.perf_counter()
    for key, expr in (
        ("I1", one * dS(TAG_RECT)),
        ("I2", g("+") * dS(TAG_RECT)),
        ("I3", (u("+") + u("-")) * dS(TAG_RECT)),
        ("I3_plus_only", u("+") * dS(TAG_RECT)),
        ("I4", one * ds(TAG_BOX)),
        ("c1_absent", one * dS(TAG_ABSENT)),
    ):
        val, e = _assemble(comm, expr)
        rec[key] = val
        rec[f"{key}_error"] = e
        if rank == 0:
            print(
                f"[GEO-37 s2] ghost={args.ghost} n={comm.size} {key} = {val!r}"
                + (f" ERROR {e}" if e else ""),
                flush=True,
            )

    # (c2) planted dropped facet: lowest rank owning any 201 facet drops its first.
    owners = comm.allgather(len(owned201) > 0)
    drop_rank = owners.index(True) if True in owners else -1
    keep = np.ones(len(facet_tags.indices), dtype=bool)
    if rank == drop_rank:
        first = int(owned201[0])
        keep &= ~((facet_tags.indices == first) & (facet_tags.values == TAG_RECT))
    ft_drop = dolfinx.mesh.meshtags(
        mesh,
        fdim,
        facet_tags.indices[keep].astype(np.int32),
        facet_tags.values[keep].astype(np.int32),
    )
    dS_drop = ufl.Measure("dS", domain=mesh, subdomain_data=ft_drop, metadata=QD)
    val, e = _assemble(comm, one * dS_drop(TAG_RECT))
    rec["drop_rank"] = drop_rank
    rec["c2_drop"] = val
    rec["c2_drop_error"] = e
    t4 = time.perf_counter()

    rec["seconds_mesh"] = t1 - t0
    rec["seconds_import"] = t2 - t1
    rec["seconds_assemble"] = t4 - t3
    if rank == 0:
        rec["completed"] = True
        with open(args.out, "w") as f:
            json.dump(rec, f, indent=1)
        print(f"[GEO-37 s2] ghost={args.ghost} n={comm.size} record: {json.dumps(rec)}", flush=True)


if __name__ == "__main__":
    main()
