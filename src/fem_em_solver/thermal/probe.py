"""Virtual temperature probe: a ball-averaged temperature time series (`MAT-19`).

`MAT-19` step 1 (PROJECT_PLAN.md §9 item 243, brief
``docs/planning/chunks/MAT-19.md``).  A probe is a ball ``B(c, a)`` of stated
volume; its reading is the mean of the temperature over that ball,

    P(t) = ∫_{B ∩ Ω} θ dV / ∫_{B ∩ Ω} dV,

computed as an **integral over owned cells**, never as a point sample.  The
point-owner rule of ``post.evaluation.evaluate_vector_field_parallel`` raises
on ``create_submesh`` meshes unless ``parent=`` is given (`POST-10`), and the
thermal fields of `MAT-10` / `MAT-12` live on such a submesh; a ball average
needs no point owner.  Each
rank sums over its owned cells only (``index_map(3).size_local``), then one
``allreduce(SUM)`` of ``(∫θ, ∫1)`` on every rank — each cell counted exactly
once at every width.

Algorithm (pure numpy): each candidate cell is sampled on the centroids of the
``8**levels`` equal-volume sub-tetrahedra of ``levels``-fold red refinement of
the reference tetrahedron (weight ``|K|/8**levels``); the indicator is a
float64 squared-distance comparison; ``θ = λ·θ_K`` is exact for P1.

Restrictions (raise): scalar Lagrange degree 1 on an affine tetrahedral mesh.
The ball is given in mesh coordinates.  Not claimed: a real probe's thermal
mass, self-heating or lead artefact.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Sequence

import numpy as np
from mpi4py import MPI

__all__ = ["ProbeBall", "ProbeSeries", "ball_average", "probe_series"]


@dataclass(frozen=True)
class ProbeBall:
    """A probe volume: the ball of ``radius`` about ``centre`` (mesh coordinates)."""

    centre: tuple
    radius: float

    def __post_init__(self):
        c = tuple(float(np.real(v)) for v in self.centre)
        if len(c) != 3:
            raise ValueError("centre must have three coordinates")
        object.__setattr__(self, "centre", c)
        r = float(self.radius)
        if not r > 0.0:
            raise ValueError(f"radius must be positive, got {self.radius!r}")
        object.__setattr__(self, "radius", r)


@dataclass
class ProbeSeries:
    """Ball-averaged series: ``values[i, j]`` is probe ``j`` at ``times[i]``."""

    steps: np.ndarray
    times: np.ndarray
    values: np.ndarray
    volumes: np.ndarray
    nominal_volumes: np.ndarray
    max_imag_ratio: float


def _red_children(t):
    """The 8 children of red refinement of tetrahedron ``t`` (4 × dim rows)."""
    x0, x1, x2, x3 = t
    m01, m02, m03 = 0.5 * (x0 + x1), 0.5 * (x0 + x2), 0.5 * (x0 + x3)
    m12, m13, m23 = 0.5 * (x1 + x2), 0.5 * (x1 + x3), 0.5 * (x2 + x3)
    return [
        np.array([x0, m01, m02, m03]),
        np.array([m01, x1, m12, m13]),
        np.array([m02, m12, x2, m23]),
        np.array([m03, m13, m23, x3]),
        np.array([m01, m02, m03, m13]),
        np.array([m01, m02, m12, m13]),
        np.array([m02, m03, m13, m23]),
        np.array([m02, m12, m13, m23]),
    ]


def _lattice(levels):
    """Barycentric centroids of the 8**levels equal-volume red sub-tetrahedra."""
    tets = [np.eye(4)]
    for _ in range(int(levels)):
        tets = [c for t in tets for c in _red_children(t)]
    return np.array([t.mean(axis=0) for t in tets])  # (8**L, 4)


def _check_space(V):
    msh = V.mesh
    from dolfinx.mesh import CellType

    if msh.topology.dim != 3 or msh.topology.cell_type != CellType.tetrahedron:
        raise ValueError("probe: tetrahedral 3-D mesh required")
    if msh.geometry.cmap.degree != 1:
        raise ValueError("probe: affine (degree-1) geometry required")
    if V.dofmap.list.shape[1] != 4 or V.dofmap.index_map_bs != 1 or V.dofmap.bs != 1:
        raise ValueError("probe: scalar Lagrange degree-1 space required")


def _geometry(V, cells, probe, lam):
    """Per-cell dof weights for one ball: returns (dofs (m,4), W (m,4), n_inside)."""
    c = np.asarray(probe.centre, dtype=np.float64)
    a = probe.radius
    xdof = np.real(V.tabulate_dof_coordinates()).astype(np.float64)
    cells = np.asarray(cells, dtype=np.int32)
    if cells.size == 0:
        return np.zeros((0, 4), np.int32), np.zeros((0, 4)), 0
    dofs = V.dofmap.list[cells]  # (n, 4)
    X = xdof[dofs]  # (n, 4, 3): the P1 vertices, ordered like the dofs
    cen = X.mean(axis=1)
    rho = np.sqrt(np.max(np.sum((X - cen[:, None, :]) ** 2, axis=2), axis=1))
    keep = np.sqrt(np.sum((cen - c) ** 2, axis=1)) - rho <= a
    dofs, X = dofs[keep], X[keep]
    if dofs.shape[0] == 0:
        return np.zeros((0, 4), np.int32), np.zeros((0, 4)), 0
    vol = np.abs(np.linalg.det(X[:, 1:, :] - X[:, :1, :])) / 6.0  # (m,)
    P = np.einsum("pk,mkd->mpd", lam, X)  # (m, npts, 3)
    inside = np.sum((P - c) ** 2, axis=2) <= a * a  # (m, npts)
    w = vol[:, None] / lam.shape[0] * inside  # (m, npts)
    W = w @ lam  # (m, 4)
    hit = np.any(inside, axis=1)
    return dofs[hit], W[hit], int(np.count_nonzero(hit))


def _lattice_points(V, cells, probe, levels=3):
    """Physical lattice points inside the ball and their weights (for `MAT-19` (G))."""
    lam = _lattice(levels)
    c = np.asarray(probe.centre, dtype=np.float64)
    xdof = np.real(V.tabulate_dof_coordinates()).astype(np.float64)
    cells = np.asarray(cells, dtype=np.int32)
    if cells.size == 0:
        return np.zeros((0, 3)), np.zeros(0)
    X = xdof[V.dofmap.list[cells]]
    vol = np.abs(np.linalg.det(X[:, 1:, :] - X[:, :1, :])) / 6.0
    P = np.einsum("pk,mkd->mpd", lam, X)
    inside = np.sum((P - c) ** 2, axis=2) <= probe.radius ** 2
    w = np.broadcast_to(vol[:, None] / lam.shape[0], inside.shape)
    return P[inside], w[inside]


def _owned_cells(V):
    return np.arange(V.mesh.topology.index_map(3).size_local, dtype=np.int32)


def _ball_average_on_cells(temperature, probe, cells, *, levels=3):
    """Ball average over the given local ``cells`` (collective; internal).

    The public path passes owned cells only; passing ghosts double-counts
    (the planted control of `MAT-19` (U3c)).
    """
    V = temperature.function_space
    _check_space(V)
    dofs, W, _ = _geometry(V, cells, probe, _lattice(levels))
    th = temperature.x.array
    num = complex(np.sum(W * th[dofs])) if dofs.shape[0] else 0.0 + 0.0j
    vol = float(np.sum(W))
    comm = V.mesh.comm
    tot = comm.allreduce(np.array([num.real, num.imag, vol]), op=MPI.SUM)
    return float(tot[0] / tot[2]), float(tot[2])


def ball_average(temperature, probe, *, levels=3):
    """``(value, volume)``: mean of ``temperature`` over ``probe`` ∩ mesh (collective)."""
    V = temperature.function_space
    return _ball_average_on_cells(temperature, probe, _owned_cells(V), levels=levels)


def probe_series(snapshots: Mapping[int, object], dt: float, probes: Sequence[ProbeBall],
                 *, levels=3) -> ProbeSeries:
    """Probe readings at every snapshot step ``n`` (time ``n·dt``), collective.

    Geometry is precomputed once per probe; each snapshot is one weighted sum
    per probe and one ``allreduce`` on every rank.
    """
    steps = np.array(sorted(int(s) for s in snapshots), dtype=int)
    if steps.size == 0:
        raise ValueError("no snapshots")
    first = snapshots[int(steps[0])]
    V = first.function_space
    _check_space(V)
    comm = V.mesh.comm
    lam = _lattice(levels)
    owned = _owned_cells(V)
    geo = [_geometry(V, owned, p, lam) for p in probes]
    vols_loc = np.array([float(np.sum(g[1])) for g in geo])
    volumes = comm.allreduce(vols_loc, op=MPI.SUM)
    used = np.unique(np.concatenate([g[0].ravel() for g in geo])) if geo else np.zeros(0, int)

    vals = np.empty((steps.size, len(probes)))
    im_max, re_max = 0.0, 0.0
    for i, s in enumerate(steps):
        th = snapshots[int(s)].x.array
        num = np.array([np.real(np.sum(g[1] * th[g[0]])) if g[0].shape[0] else 0.0
                        for g in geo])
        vals[i] = comm.allreduce(num, op=MPI.SUM) / volumes
        if used.size:
            im_max = max(im_max, float(np.max(np.abs(np.imag(th[used])))))
            re_max = max(re_max, float(np.max(np.abs(np.real(th[used])))))
    im_max = comm.allreduce(im_max, op=MPI.MAX)
    re_max = comm.allreduce(re_max, op=MPI.MAX)
    ratio = im_max / re_max if re_max > 0.0 else 0.0
    nominal = np.array([4.0 * np.pi * p.radius ** 3 / 3.0 for p in probes])
    return ProbeSeries(steps=steps, times=steps * float(dt), values=vals,
                       volumes=np.asarray(volumes), nominal_volumes=nominal,
                       max_imag_ratio=float(ratio))
