"""`GEO-37` step 1 (§9 item 222): does an interior dim-2 physical group reach
DolfinX's ``facet_tags`` on the 0.11 image?

Probe script (never collected: ``python_files = test_*.py``).  A unit box with
one embedded interior rectangle (0.8 x 0.6 in z = 0.47, area 0.48) tagged 201,
the six box faces tagged 101, imported through the project shim
``_model_to_mesh`` with the default partitioner.  Areas are measured as numpy
owned-facet sums reduced with ``allreduce`` -- never a UFL ``dS`` integral (the
2026-08 hang, `PORT-1` 3b-iv, was that assembly on an unghosted mesh).

Rank 0 writes a JSON record (``--out``) whose ``completed: true`` is written
last; the asserting test is ``tests/mesh/test_geo37_interior_facet.py``.
"""

import argparse
import json
import time

import gmsh
import numpy as np
from mpi4py import MPI

import dolfinx

from fem_em_solver.io.mesh import (
    _model_to_mesh,
    _raise_geometry_failure_on_every_rank,
)

RECT = (0.1, 0.2, 0.47, 0.8, 0.6)  # x, y, z, w, h
RECT_AREA = 0.8 * 0.6
H = 0.08
TAG_VOL, TAG_RECT, TAG_BOX = 1, 201, 101


def _build_model():
    """Rank 0 only.  Returns N_gmsh (rectangle triangle count)."""
    occ = gmsh.model.occ
    box = occ.addBox(0, 0, 0, 1, 1, 1)
    rect = occ.addRectangle(*RECT)
    occ.fragment([(3, box)], [(2, rect)])
    occ.synchronize()

    vols = gmsh.model.getEntities(3)
    if len(vols) != 1:
        raise RuntimeError(f"GEO-37 fixture: expected 1 volume, got {vols}")
    vol = vols[0][1]

    rect_tags, box_tags = [], []
    for _, s in gmsh.model.getEntities(2):
        x0, y0, z0, x1, y1, z1 = gmsh.model.getBoundingBox(2, s)
        # OCC bounding boxes are padded by the shape tolerance (~1e-7), so the
        # brief's literal 1e-12 z-span classifier found no face in window W1
        # (20261005T033147Z_GEO-37-W1.log); classify at 1e-6 and keep the
        # exact area check (getMass, 1e-12) below as the guard.
        if abs(z1 - z0) < 1e-6 and abs(0.5 * (z0 + z1) - RECT[2]) < 1e-6:
            rect_tags.append(s)
        else:
            box_tags.append(s)
    if len(rect_tags) != 1:
        raise RuntimeError(f"GEO-37 fixture: expected 1 interior face, got {rect_tags}")
    rect_s = rect_tags[0]
    mass = occ.getMass(2, rect_s)
    if abs(mass - RECT_AREA) > 1e-12:
        raise RuntimeError(f"GEO-37 fixture: interior face area {mass!r} != 0.48")
    if len(box_tags) != 6:
        raise RuntimeError(f"GEO-37 fixture: expected 6 box faces, got {box_tags}")
    gmsh.model.addPhysicalGroup(3, [vol], TAG_VOL)
    gmsh.model.addPhysicalGroup(2, [rect_s], TAG_RECT)
    gmsh.model.addPhysicalGroup(2, box_tags, TAG_BOX)

    gmsh.option.setNumber("Mesh.MeshSizeMax", H)
    gmsh.option.setNumber("Mesh.MeshSizeMin", H)
    gmsh.model.mesh.generate(3)
    _, elem_tags, _ = gmsh.model.mesh.getElements(2, rect_s)
    n_gmsh = int(sum(len(t) for t in elem_tags))
    return n_gmsh


def _tri_areas(mesh, facets):
    if len(facets) == 0:
        return np.zeros(0)
    g = dolfinx.mesh.entities_to_geometry(mesh, 2, np.asarray(facets, dtype=np.int32))
    x = mesh.geometry.x
    a, b, c = x[g[:, 0]], x[g[:, 1]], x[g[:, 2]]
    return 0.5 * np.linalg.norm(np.cross(b - a, c - a), axis=1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    comm = MPI.COMM_WORLD
    rank = comm.rank
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)

    t0 = time.perf_counter()
    err, n_gmsh = None, None
    if rank == 0:
        try:
            n_gmsh = _build_model()
        except BaseException as e:  # noqa: BLE001 -- GEO-23 pattern
            err = e
            gmsh.finalize()
    _raise_geometry_failure_on_every_rank(comm, 0, "GEO-37 probe", err, H)
    n_gmsh = comm.bcast(n_gmsh, root=0)
    t1 = time.perf_counter()

    mesh, cell_tags, facet_tags = _model_to_mesh(gmsh.model, comm, 0, gdim=3)
    gmsh.finalize()
    comm.Barrier()
    t2 = time.perf_counter()

    tdim = mesh.topology.dim
    n_cells = comm.allreduce(mesh.topology.index_map(tdim).size_local, op=MPI.SUM)
    n_ghost_cells = comm.allreduce(mesh.topology.index_map(tdim).num_ghosts, op=MPI.SUM)
    n_owned = mesh.topology.index_map(2).size_local

    rec = {
        "comm_size": comm.size,
        "N_gmsh": n_gmsh,
        "cells": n_cells,
        "cells_ghost": n_ghost_cells,
        "h": H,
        "rect_area_closed_form": RECT_AREA,
    }
    owned_201_local = None
    for t in (TAG_RECT, TAG_BOX):
        local = facet_tags.find(t)
        owned = local[local < n_owned]
        a_owned = _tri_areas(mesh, owned)
        a_local = _tri_areas(mesh, local)
        rec[f"count_{t}"] = int(comm.allreduce(len(owned), op=MPI.SUM))
        rec[f"area_{t}"] = float(comm.allreduce(float(a_owned.sum()), op=MPI.SUM))
        rec[f"count_local_{t}"] = int(comm.allreduce(len(local), op=MPI.SUM))
        rec[f"area_local_{t}"] = float(comm.allreduce(float(a_local.sum()), op=MPI.SUM))
        rec[f"per_rank_owned_{t}"] = comm.gather(int(len(owned)), root=0)
        rec[f"per_rank_local_{t}"] = comm.gather(int(len(local)), root=0)
        if t == TAG_RECT:
            owned_201_local = a_owned

    # (c1) planted control: lowest rank owning any 201 facet drops its first.
    has = len(owned_201_local) > 0
    owners = comm.allgather(has)
    drop_rank = owners.index(True) if True in owners else -1
    a_kept = owned_201_local[1:] if rank == drop_rank else owned_201_local
    rec["drop_rank"] = drop_rank
    rec["count_drop"] = int(comm.allreduce(len(a_kept), op=MPI.SUM))
    rec["area_drop"] = float(comm.allreduce(float(a_kept.sum()), op=MPI.SUM))
    t3 = time.perf_counter()

    rec["seconds_mesh"] = t1 - t0
    rec["seconds_import"] = t2 - t1
    rec["seconds_measure"] = t3 - t2
    if rank == 0:
        rec["completed"] = True
        with open(args.out, "w") as f:
            json.dump(rec, f, indent=1)
        print(f"[GEO-37] n={comm.size} record: {json.dumps(rec)}", flush=True)


if __name__ == "__main__":
    main()
