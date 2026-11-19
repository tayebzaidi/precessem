"""GEO-33 step 0a: the TE101 box with its y = b/2 facet plane duplicated.

The cavity [0, 1.0] x [0, 0.8] x [0, 0.6] is built as the lower half
[0, 1.0] x [0, 0.4] x [0, 0.6] (``create_box`` on COMM_SELF, (9, 4, 6)
tetrahedra) mirrored across y = 0.4 on TH-16's pattern
(``tests/validation/test_th16_symmetry_plane_driven.py:116-154``, cells
reordered ``[1, 0, 2, 3]`` for the reflection) **without merging the plane
vertices**: the upper half's vertices are offset by n_v, so the plane is a
crack with two lips, each an exterior facet set. ``merge=True`` merges the
plane vertices (the single-valued sheet; step 0b's 1.000 control).

Anchors (asserted, exact construction identities; §9 item 149):
(L1) |lip lower| = |lip upper| = a*d = 0.6 within 1e-9 relative;
(L2) 108 owned facets per lip (9*6*2), sorted midpoints pair up between
     the lips within 1e-12 m (conforming twins);
(L3) every owned facet on y = 0.4 has exactly one local adjacent cell, no
     cell touches both lips, and the two lips share no global vertex;
(L4) total boundary = 2(ab + bd + ad) + 2ad = 4.96 and volume 0.48, 1e-9;
(L5) 2592 cells, 700 vertices (global);
(L6) the same readings at -n 1 and -n 2 (two harness windows).
Negative control (asserted binary): merge=True gives 0 lip facets,
boundary 3.76, 630 vertices.

The builder lives in ``fem_em_solver.core.cracked_box`` since step 0b.
"""

import numpy as np
import pytest
import ufl
from dolfinx import default_scalar_type, fem
from dolfinx import mesh as dmesh
from mpi4py import MPI

from fem_em_solver.core.cracked_box import (
    LIP_LOWER,
    LIP_UPPER,
    X0,
    X1,
    Y0,
    Y1,
    Z0,
    Z1,
    cracked_box_mesh,
)

A, B, D = 1.0, 0.8, 0.6  # EDGES of test_cavity_leontovich_q.py:42
Y_PLANE = B / 2.0
PLANE_TOL = 1e-12  # exact plane coordinate; never widened
TWIN_TOL = 1e-12  # m
REL_TOL = 1e-9
FACE_TOL = 1e-9


EXPECTED_LIP_AREA = A * D  # 0.6
EXPECTED_LIP_FACETS = 9 * 6 * 2  # 108
EXPECTED_BOUNDARY_CRACKED = 2 * (A * B + B * D + A * D) + 2 * A * D  # 4.96
EXPECTED_BOUNDARY_MERGED = 2 * (A * B + B * D + A * D)  # 3.76
EXPECTED_VOLUME = A * B * D  # 0.48
EXPECTED_CELLS = 2 * 6 * 9 * 4 * 6  # 2592
EXPECTED_VERTICES_CRACKED = 2 * 10 * 5 * 7  # 700
EXPECTED_VERTICES_MERGED = 2 * 10 * 5 * 7 - 10 * 7  # 630


def _allsum(comm, v):
    return comm.allreduce(v, op=MPI.SUM)


def _integral(msh, measure):
    one = fem.Constant(msh, default_scalar_type(1.0))
    local = fem.assemble_scalar(fem.form(one * measure))
    return float(np.real(_allsum(msh.comm, local)))


def _owned(msh, facets):
    n_local = msh.topology.index_map(2).size_local
    return facets[facets < n_local]


def _readings(msh, tags):
    comm = msh.comm
    ds = ufl.Measure("ds", domain=msh, subdomain_data=tags)
    r = {
        "cells": msh.topology.index_map(3).size_global,
        "vertices": msh.topology.index_map(0).size_global,
        "boundary": _integral(msh, ufl.ds(domain=msh)),
        "volume": _integral(msh, ufl.dx(domain=msh)),
        "lip_lower_area": _integral(msh, ds(LIP_LOWER)),
        "lip_upper_area": _integral(msh, ds(LIP_UPPER)),
        "faces_area": sum(_integral(msh, ds(t)) for t in (X0, X1, Y0, Y1, Z0, Z1)),
    }
    lo = _owned(msh, tags.find(LIP_LOWER))
    up = _owned(msh, tags.find(LIP_UPPER))
    r["lip_lower_facets"] = _allsum(comm, int(lo.size))
    r["lip_upper_facets"] = _allsum(comm, int(up.size))
    return r, lo, up


@pytest.fixture(scope="module")
def cracked():
    msh, tags, spread = cracked_box_mesh(merge=False)
    if MPI.COMM_WORLD.rank == 0:
        print(f"[GEO-33 0a] half-box plane y-spread max|y-0.4| = {spread!r}")
    if msh is None:
        pytest.fail(f"plane coordinate not exact: y-spread {spread!r} (stop; never widen)")
    return msh, tags


def test_lips_decoupled_and_conforming(cracked):
    msh, tags = cracked
    comm = msh.comm
    r, lo, up = _readings(msh, tags)
    if comm.rank == 0:
        print(
            f"[GEO-33 0a] cracked -n {comm.size}: " + ", ".join(f"{k}={v!r}" for k, v in r.items())
        )

    # (L1) each lip a*d.
    for key in ("lip_lower_area", "lip_upper_area"):
        rel = abs(r[key] / EXPECTED_LIP_AREA - 1.0)
        if comm.rank == 0:
            print(f"[GEO-33 0a] L1 {key} = {r[key]!r}, rel {rel:.3e}")
        assert rel < REL_TOL, (key, r[key])

    # (L2) facet counts and conforming twins.
    assert r["lip_lower_facets"] == EXPECTED_LIP_FACETS, r["lip_lower_facets"]
    assert r["lip_upper_facets"] == EXPECTED_LIP_FACETS, r["lip_upper_facets"]
    mlo = dmesh.compute_midpoints(msh, 2, lo.astype(np.int32)) if lo.size else np.zeros((0, 3))
    mup = dmesh.compute_midpoints(msh, 2, up.astype(np.int32)) if up.size else np.zeros((0, 3))
    glo = comm.gather(mlo, root=0)
    gup = comm.gather(mup, root=0)
    twin = None
    if comm.rank == 0:
        glo, gup = np.vstack(glo), np.vstack(gup)

        def _sorted(p):
            k = np.round(p, 9)
            return p[np.lexsort((k[:, 2], k[:, 1], k[:, 0]))]

        twin = float(np.max(np.linalg.norm(_sorted(glo) - _sorted(gup), axis=1)))
        print(
            f"[GEO-33 0a] L2 lip facets {r['lip_lower_facets']}/{r['lip_upper_facets']}, "
            f"max twin midpoint distance {twin:.3e} m"
        )
    twin = comm.bcast(twin, root=0)
    assert twin < TWIN_TOL, twin

    # (L3) one adjacent cell per plane facet; no cell on both lips; no shared vertex.
    tdim, fdim = 3, 2
    msh.topology.create_connectivity(fdim, tdim)
    msh.topology.create_connectivity(tdim, fdim)
    msh.topology.create_connectivity(fdim, 0)
    f2c = msh.topology.connectivity(fdim, tdim)
    plane = _owned(
        msh, dmesh.locate_entities(msh, fdim, lambda p: np.abs(p[1] - Y_PLANE) < PLANE_TOL)
    )
    links = np.array([f2c.links(f).size for f in plane], dtype=np.int64)
    n_plane = _allsum(comm, int(plane.size))
    n_bad = _allsum(comm, int(np.sum(links != 1)))
    lo_all, up_all = tags.find(LIP_LOWER), tags.find(LIP_UPPER)
    c2f = msh.topology.connectivity(tdim, fdim)
    n_cells_local = msh.topology.index_map(tdim).size_local
    lo_set, up_set = set(lo_all.tolist()), set(up_all.tolist())
    both = 0
    for c in range(n_cells_local):
        fs = c2f.links(c)
        if any(int(f) in lo_set for f in fs) and any(int(f) in up_set for f in fs):
            both += 1
    both = _allsum(comm, both)
    f2v = msh.topology.connectivity(fdim, 0)
    vmap = msh.topology.index_map(0)

    def _gverts(fs):
        if fs.size == 0:
            return np.zeros(0, dtype=np.int64)
        loc = np.unique(np.concatenate([f2v.links(f) for f in fs]))
        return vmap.local_to_global(loc.astype(np.int32))

    vlo = comm.gather(_gverts(lo), root=0)
    vup = comm.gather(_gverts(up), root=0)
    shared = None
    if comm.rank == 0:
        vlo, vup = np.unique(np.concatenate(vlo)), np.unique(np.concatenate(vup))
        shared = int(np.intersect1d(vlo, vup).size)
        print(
            f"[GEO-33 0a] L3 owned facets on y=0.4: {n_plane}, with != 1 local cell: "
            f"{n_bad}; cells touching both lips: {both}; lip vertices "
            f"{vlo.size}/{vup.size}, shared: {shared}"
        )
    shared = comm.bcast(shared, root=0)
    assert n_plane == 2 * EXPECTED_LIP_FACETS, n_plane
    assert n_bad == 0, n_bad
    assert both == 0, both
    assert shared == 0, shared

    # (L4) boundary and volume.
    rb = abs(r["boundary"] / EXPECTED_BOUNDARY_CRACKED - 1.0)
    rv = abs(r["volume"] / EXPECTED_VOLUME - 1.0)
    rf = abs(r["faces_area"] / EXPECTED_BOUNDARY_MERGED - 1.0)
    if comm.rank == 0:
        print(
            f"[GEO-33 0a] L4 boundary {r['boundary']!r} (rel {rb:.3e} vs 4.96), "
            f"volume {r['volume']!r} (rel {rv:.3e} vs 0.48), six faces "
            f"{r['faces_area']!r} (rel {rf:.3e} vs 3.76)"
        )
    assert rb < REL_TOL, r["boundary"]
    assert rv < REL_TOL, r["volume"]
    assert rf < REL_TOL, r["faces_area"]

    # (L5) global counts.
    if comm.rank == 0:
        print(f"[GEO-33 0a] L5 cells {r['cells']}, vertices {r['vertices']}")
    assert r["cells"] == EXPECTED_CELLS, r["cells"]
    assert r["vertices"] == EXPECTED_VERTICES_CRACKED, r["vertices"]


def test_merged_control_is_single_valued():
    msh, tags, spread = cracked_box_mesh(merge=True)
    assert msh is not None, spread
    comm = msh.comm
    r, _, _ = _readings(msh, tags)
    short = 1.0 - r["boundary"] / EXPECTED_BOUNDARY_CRACKED
    if comm.rank == 0:
        print(
            f"[GEO-33 0a] control merge=True -n {comm.size}: lip facets "
            f"{r['lip_lower_facets']}/{r['lip_upper_facets']}, boundary "
            f"{r['boundary']!r} ({100 * short:.2f} % short of 4.96), vertices "
            f"{r['vertices']}, cells {r['cells']}, volume {r['volume']!r}"
        )
    assert r["lip_lower_facets"] == 0 and r["lip_upper_facets"] == 0
    assert abs(r["boundary"] / EXPECTED_BOUNDARY_MERGED - 1.0) < REL_TOL, r["boundary"]
    assert r["vertices"] == EXPECTED_VERTICES_MERGED, r["vertices"]
    assert r["cells"] == EXPECTED_CELLS, r["cells"]
