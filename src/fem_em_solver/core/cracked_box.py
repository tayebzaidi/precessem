"""GEO-33: the TE101 box with its y = b/2 facet plane duplicated (a crack).

Moved from ``tests/mesh/test_geo33_cracked_cavity_lips.py`` (step 0a) by step
0b. The cavity [0, 1.0] x [0, 0.8] x [0, 0.6] is built as the lower half
[0, 1.0] x [0, 0.4] x [0, 0.6] (``create_box`` on COMM_SELF, (9, 4, 6)
tetrahedra) mirrored across y = 0.4 on TH-16's pattern (cells reordered
``[1, 0, 2, 3]`` for the reflection) **without merging the plane vertices**:
the upper half's vertices are offset by n_v, so the plane is a crack with two
lips, each an exterior facet set. ``merge=True`` merges the plane vertices
(the single-valued sheet; step 0b's 1.000 control).

Importable in the real and the complex build (no complex-only imports).
"""
import basix.ufl
import numpy as np
import ufl
from dolfinx import mesh as dmesh
from mpi4py import MPI

A, B, D = 1.0, 0.8, 0.6
Y_PLANE = B / 2.0
HALF_DIVISIONS = (9, 4, 6)
PLANE_TOL = 1e-12  # exact plane coordinate; never widened
FACE_TOL = 1e-9

X0, X1, Y0, Y1, Z0, Z1 = 1, 2, 3, 4, 5, 6
LIP_LOWER, LIP_UPPER = 11, 12


def _distribute(x, cells, comm):
    """Rank-0 arrays -> distributed mesh (TH-16's ``_distribute``)."""
    if comm.rank != 0:
        x = np.empty((0, 3), dtype=np.float64)
        cells = np.empty((0, 4), dtype=np.int64)
    domain = ufl.Mesh(basix.ufl.element("Lagrange", "tetrahedron", 1, shape=(3,)))
    partitioner = dmesh.create_cell_partitioner(dmesh.GhostMode.shared_facet, 2)
    return dmesh.create_mesh(
        comm, cells=np.ascontiguousarray(cells), x=np.ascontiguousarray(x), e=domain,
        partitioner=partitioner,
    )


def _half_arrays():
    h = dmesh.create_box(
        MPI.COMM_SELF,
        [np.zeros(3), np.array([A, Y_PLANE, D])],
        HALF_DIVISIONS,
        dmesh.CellType.tetrahedron,
    )
    x = np.asarray(h.geometry.x, dtype=np.float64)[:, :3].copy()
    # The accessor the DeprecationWarning names ("Use dofmaps[0] instead").
    dm = h.geometry.dofmaps
    dm = dm(0) if callable(dm) else dm[0]
    cells = np.asarray(dm, dtype=np.int64).reshape(-1, 4).copy()
    return x, cells


def _tag(msh):
    tdim, fdim = 3, 2
    msh.topology.create_connectivity(fdim, tdim)
    ext = dmesh.exterior_facet_indices(msh.topology).astype(np.int32)
    mid = dmesh.compute_midpoints(msh, fdim, ext) if ext.size else np.zeros((0, 3))
    vals = np.zeros(ext.size, dtype=np.int32)
    for tag, axis, c in (
        (X0, 0, 0.0), (X1, 0, A), (Y0, 1, 0.0), (Y1, 1, B), (Z0, 2, 0.0), (Z1, 2, D),
    ):
        vals[np.abs(mid[:, axis] - c) < FACE_TOL] = tag
    on_plane = np.abs(mid[:, 1] - Y_PLANE) < PLANE_TOL
    if on_plane.any():
        f2c = msh.topology.connectivity(fdim, tdim)
        adj = np.array([f2c.links(f)[0] for f in ext[on_plane]], dtype=np.int32)
        cmid = dmesh.compute_midpoints(msh, tdim, adj)
        vals[np.flatnonzero(on_plane)] = np.where(
            cmid[:, 1] < Y_PLANE, LIP_LOWER, LIP_UPPER
        ).astype(np.int32)
    keep = vals != 0
    idx, vals = ext[keep], vals[keep]
    order = np.argsort(idx, kind="stable")
    return dmesh.meshtags(msh, fdim, idx[order], vals[order])


def _septum_tags(msh, septum_tag):
    """Interior facets on y = 0.4 (the merged mesh's single-valued sheet)."""
    fdim = 2
    msh.topology.create_connectivity(fdim, 3)
    idx = dmesh.locate_entities(msh, fdim, lambda p: np.abs(p[1] - Y_PLANE) < PLANE_TOL)
    idx = np.unique(np.asarray(idx, dtype=np.int32))
    vals = np.full(idx.size, int(septum_tag), dtype=np.int32)
    return dmesh.meshtags(msh, fdim, idx, vals)


def cracked_box_mesh(merge: bool, comm=MPI.COMM_WORLD, septum_tag=None):
    """Build the mirrored box; ``merge=False`` duplicates the y = 0.4 plane.

    Returns ``(mesh, facet_tags, plane_spread)``, plus ``septum_tags`` as a
    fourth element when ``septum_tag`` is given. ``plane_spread`` is the
    largest |y - 0.4| among the half-box's plane vertices (must be exactly 0;
    otherwise ``mesh`` is ``None``).
    """
    x = cells = None
    spread = None
    if comm.rank == 0:
        xh, ch = _half_arrays()
        near = np.abs(xh[:, 1] - Y_PLANE) < 1e-6
        spread = float(np.max(np.abs(xh[near, 1] - Y_PLANE))) if near.any() else np.inf
        xm = xh.copy()
        xm[:, 1] = 2.0 * Y_PLANE - xm[:, 1]
        x = np.vstack([xh, xm])
        cells = np.vstack([ch, ch[:, [1, 0, 2, 3]] + xh.shape[0]])
        if merge:
            key = np.round(x, 9)
            _, first, inverse = np.unique(
                key, axis=0, return_index=True, return_inverse=True
            )
            x = x[first]
            cells = inverse.reshape(-1)[cells]
    spread = comm.bcast(spread, root=0)
    if spread != 0.0:
        return (None, None, spread) if septum_tag is None else (None, None, spread, None)
    msh = _distribute(x, cells, comm)
    tags = _tag(msh)
    if septum_tag is None:
        return msh, tags, spread
    return msh, tags, spread, _septum_tags(msh, septum_tag)
