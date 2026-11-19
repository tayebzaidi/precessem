"""Measurement-only probe (§9 item 147, `WF-12` step 1 unblock): hole vs solid
conductor-wall area at F-human scale, swept over ``conductor_resolution``.

Asserts nothing, imported by nothing.  One fresh process per (rung, route):

    python3 -u tests/mesh/probe_wf12_cavity_wall_area.py --route hole|solid \
        --conductor-resolution 1.6e-3

Builder = item 129 W1's ``build()``: ``_params(F_HUMAN_RING_RADIUS,
scale_sizing=False)``, ``ring_sheet_orientation="longitudinal"``,
``as_hole=True`` on the hole route only, ``conductor_resolution`` overridden
in the kwargs dict, called on ``MeshGenerator.birdcage_port_domain``.

Prints: cells, build time, peak RSS, total wall area A (hole: exterior tag
``BIRDCAGE_CONDUCTOR_SURFACE_TAG``; solid: interface tag 1 vs every other
tag), A split by neighbour class with facet counts (owned facets only),
|Σ_class − A|, and on the hole route the CAD cavity area and A/A_CAD − 1.
"""

from __future__ import annotations

import argparse
import resource
import time

import numpy as np
from mpi4py import MPI

AIR = 2
PHANTOM = 3
SOLID_CONDUCTOR_IFACE = 900


def _class_of(tag: int) -> str:
    if tag == AIR:
        return "air(2)"
    if tag == PHANTOM:
        return "phantom(3)"
    if 101 <= tag <= 116:
        return "leg-port(101-116)"
    if 117 <= tag <= 148 or 217 <= tag <= 248:
        return "ring-port(117-148/217-248)"
    return f"other({tag})"


def _owned_facet_split(msh, facet_indices, neighbour_tag, comm):
    """Per-class (area, count) over owned facets, triangle areas from geometry."""
    fdim = msh.topology.dim - 1
    n_owned = msh.topology.index_map(fdim).size_local
    keep = facet_indices < n_owned
    facets = facet_indices[keep]
    nbr = neighbour_tag[keep]
    import dolfinx

    geom_ids = dolfinx.cpp.mesh.entities_to_geometry(
        msh._cpp_object, fdim, facets.astype(np.int32), False
    )
    x = msh.geometry.x
    p0, p1, p2 = x[geom_ids[:, 0]], x[geom_ids[:, 1]], x[geom_ids[:, 2]]
    area = 0.5 * np.linalg.norm(np.cross(p1 - p0, p2 - p0), axis=1)
    local = {}
    for t, a in zip(nbr, area):
        c = _class_of(int(t))
        s = local.setdefault(c, [0.0, 0])
        s[0] += float(a)
        s[1] += 1
    gathered = comm.allgather(local)
    out = {}
    for d in gathered:
        for c, (a, n) in d.items():
            s = out.setdefault(c, [0.0, 0])
            s[0] += a
            s[1] += n
    return out


def _cell_value(msh, cell_tags):
    import dolfinx

    space = dolfinx.fem.functionspace(msh, ("DG", 0))
    f = dolfinx.fem.Function(space)
    f.x.array[:] = 0.0
    c2d = space.dofmap.list.reshape(-1)
    f.x.array[c2d[cell_tags.indices]] = cell_tags.values
    f.x.scatter_forward()
    return np.rint(np.real(f.x.array[c2d])).astype(np.int32)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--route", choices=("hole", "solid"), required=True)
    ap.add_argument("--conductor-resolution", type=float, required=True)
    args = ap.parse_args()

    from fem_em_solver.io.mesh import (
        BIRDCAGE_CONDUCTOR_SURFACE_TAG,
        MeshGenerator,
        _interface_facet_tags,
    )
    from tests.mesh.probe_two_torus_conductor_hole import (
        _census,
        _exterior_facet_area,
        _facet_area,
    )
    from tests.validation.test_birdcage_f_human_rung import (
        F_HUMAN_RING_RADIUS,
        _params,
    )

    comm = MPI.COMM_WORLD
    tag = f"[WF-12 probe {args.route} h={args.conductor_resolution:.2e}]"

    kwargs = _params(F_HUMAN_RING_RADIUS, scale_sizing=False)
    kwargs["ring_sheet_orientation"] = "longitudinal"
    kwargs["conductor_resolution"] = args.conductor_resolution
    if args.route == "hole":
        kwargs["as_hole"] = True
    t0 = time.perf_counter()
    msh, cells, facets, diag = MeshGenerator.birdcage_port_domain(
        comm=comm, return_diagnostics=True, **kwargs
    )
    t_build = time.perf_counter() - t0
    tdim = msh.topology.dim
    fdim = tdim - 1
    msh.topology.create_entities(fdim)
    msh.topology.create_connectivity(fdim, tdim)
    n_cells = int(msh.topology.index_map(tdim).size_global)
    n_owned_cells = msh.topology.index_map(tdim).size_local
    census = _census(cells.values, comm, indices=cells.indices, n_owned=n_owned_cells)
    rss_gib = comm.allreduce(
        resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0**2, op=MPI.SUM
    )
    if comm.rank == 0:
        print(
            f"{tag} cells {n_cells}; build {t_build:.1f} s (gmsh mesh_wall_time_s "
            f"{diag.get('mesh_wall_time_s', float('nan'))}); peak RSS summed "
            f"{rss_gib:.3f} GiB; -n {comm.size}; census {census}",
            flush=True,
        )

    cell_value = _cell_value(msh, cells)
    f2c = msh.topology.connectivity(fdim, tdim)
    offsets, links = f2c.offsets, f2c.array

    if args.route == "hole":
        A = _exterior_facet_area(msh, facets, BIRDCAGE_CONDUCTOR_SURFACE_TAG, comm)
        idx = np.asarray(facets.indices)[
            np.asarray(facets.values) == BIRDCAGE_CONDUCTOR_SURFACE_TAG
        ]
        nbr = cell_value[links[offsets[idx]]]
    else:
        others = sorted(t for t in census if t != 1)
        iface = _interface_facet_tags(
            msh, cells, {SOLID_CONDUCTOR_IFACE: tuple((1, t) for t in others)}
        )
        A = _facet_area(msh, iface, SOLID_CONDUCTOR_IFACE, comm)
        idx = np.asarray(iface.indices)[np.asarray(iface.values) == SOLID_CONDUCTOR_IFACE]
        a = cell_value[links[offsets[idx]]]
        b = cell_value[links[offsets[idx] + 1]]
        nbr = np.where(a == 1, b, a)

    split = _owned_facet_split(msh, np.asarray(idx), nbr, comm)
    total = sum(v[0] for v in split.values())
    if comm.rank == 0:
        print(f"{tag} wall area A = {A:.10e} m^2", flush=True)
        for c in sorted(split):
            print(
                f"{tag}   class {c}: area {split[c][0]:.10e} m^2, facets {split[c][1]}", flush=True
            )
        print(
            f"{tag} sum_class {total:.10e}; |sum_class - A| = {abs(total - A):.3e}; "
            f"facets {sum(v[1] for v in split.values())}",
            flush=True,
        )
        if args.route == "hole":
            cad = diag.get("conductor_cavity_cad", {})
            a_cad = float(cad.get("area_m2", float("nan")))
            print(
                f"{tag} A_CAD {a_cad:.10e} m^2 ({cad.get('surfaces')} surfaces); "
                f"A/A_CAD - 1 = {A / a_cad - 1.0:+.6e}; A_CAD - A = {a_cad - A:.6e}",
                flush=True,
            )


if __name__ == "__main__":
    main()
