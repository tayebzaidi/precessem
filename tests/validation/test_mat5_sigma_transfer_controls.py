"""`MAT-5` (D′) against planted failures — the σ-transfer identity, no EM solve
(§9 item 109, the 2026-09-30 10:30 review's audit follow-up).

(D′) (`test_mat5_em_coupled_sphere.py:330–336`) compares ∫σ x·x over the
parent's tag-1 cells with the same integral on the submesh.  It read exactly
0 in `MAT-5`'s closing window; whether it can *see* a wrong transfer was
argued, not measured.  Here the sub σ is the closed-form profile's σ(T) per
cell — θ(r) = a − b r² from `robin_sphere_steady_temperature` at q*, h = 3 —
and the parent starts at the sub field's volume mean σ̄ (what
``material_map[1]`` carries before the overwrite).

(p0) the committed `_sub_to_parent` reads ≤ `D_BAND` — asserted.  Planted
failures, each asserted only as a miss of `D_BAND` (sizes printed beside
their predictions, rule (e)): (n1) no transfer (predicted 1.42e-3), (n2) a
per-rank shifted transfer, (n3) a per-rank reversed transfer.  (n4) a
ghost-blind transfer (owned cells only, no ``scatter_forward``) is predicted
invisible at ``-n 2`` — printed, not asserted either way.
"""

from __future__ import annotations

import numpy as np
import ufl
from dolfinx import fem
from dolfinx import mesh as dmesh
from mpi4py import MPI

from fem_em_solver.materials.temperature import LinearTemperatureConductivity
from fem_em_solver.utils.analytical import robin_sphere_steady_temperature
from tests.complex_mode import complex_only
from tests.materials.test_mat5_step2_lumped_root import (
    ALPHA,
    H_STEP2,
    Q_STAR,
    SIGMA_0,
    T_REF,
)
from tests.materials.test_phantom_material_model import build_preset_sphere_mesh
from tests.validation.test_lossy_sphere_sar import SPHERE_RADIUS, SPHERE_TAG
from tests.validation.test_mat5_em_coupled_sphere import D_BAND, _real, _sub_to_parent
from tests.validation.test_mat10_sar_to_heat_sphere import (
    T_INF,
    K,
    _allsum,
    _say,
    _sphere_submesh,
)

R = SPHERE_RADIUS
N_PARENT = 74_020  # the row's counts (MAT-10 / MAT-5 closing log)
N_SUB = 19_496
N1_PREDICTED = 1.42e-3  # (Δσ/σ̄)(5/7 − 3/5), continuum quadratic profile


def _moment_rel(msh, cell_tags, sub, sig_par, sig_sub):
    comm = msh.comm
    xp = ufl.SpatialCoordinate(msh)
    xs = ufl.SpatialCoordinate(sub)
    dx1 = ufl.Measure("dx", domain=msh, subdomain_data=cell_tags)(SPHERE_TAG)
    m_par = _real(_allsum(comm, sig_par * ufl.inner(xp, xp) * dx1))
    m_sub = _real(_allsum(comm, sig_sub * ufl.inner(xs, xs) * ufl.dx(domain=sub)))
    return abs(m_par - m_sub) / abs(m_sub), m_par, m_sub


@complex_only
def test_dprime_sees_planted_transfer_failures():
    comm = MPI.COMM_WORLD
    msh, cell_tags, _ = build_preset_sphere_mesh(comm)
    sub, parent_ids = _sphere_submesh(msh, cell_tags)
    n_parent = msh.topology.index_map(3).size_global
    n_sub = sub.topology.index_map(3).size_global
    _say(
        f"\nMAT5-D' cells: parent = {n_parent}  sub = {n_sub}  (row: {N_PARENT} / {N_SUB}); "
        f"-n {comm.size}"
    )

    # Sub σ from the closed-form profile, per cell midpoint (no solve).
    th_c = (
        float(robin_sphere_steady_temperature(0.0, radius=R, k=K, q=Q_STAR, h=H_STEP2, t_inf=T_INF))
        - T_INF
    )
    th_R = (
        float(robin_sphere_steady_temperature(R, radius=R, k=K, q=Q_STAR, h=H_STEP2, t_inf=T_INF))
        - T_INF
    )
    law = LinearTemperatureConductivity(sigma_ref=SIGMA_0, alpha_per_k=ALPHA, t_ref=T_REF)
    Ws = fem.functionspace(sub, ("DG", 0))
    sig_sub = fem.Function(Ws)
    simap = sub.topology.index_map(3)
    n_sub_local = simap.size_local + simap.num_ghosts
    cells = np.arange(n_sub_local, dtype=np.int32)
    mids = dmesh.compute_midpoints(sub, 3, cells)
    r2 = np.sum(mids**2, axis=1)
    theta = th_c - (th_c - th_R) / R**2 * r2
    sdofs_all = Ws.dofmap.list[cells, 0]
    sig_sub.x.array[sdofs_all] = law.sigma(T_INF + theta)
    sig_sub.x.scatter_forward()

    vol = _real(_allsum(comm, 1.0 * ufl.dx(domain=sub)))
    sigma_bar = _real(_allsum(comm, sig_sub * ufl.dx(domain=sub))) / vol
    d_sigma = SIGMA_0 * ALPHA * (th_c - th_R)
    _say(
        f"MAT5-D' profile: theta(0) = {th_c:.6f} K  theta(R) = {th_R:.6f} K; sigma_bar = "
        f"{sigma_bar:.9f}  delta sigma = {d_sigma:.6f} (predicted 0.5729 / 0.007125)"
    )

    Wp = fem.functionspace(msh, ("DG", 0))

    def parent_at_mean():
        f = fem.Function(Wp)
        f.x.array[:] = sigma_bar
        f.x.scatter_forward()
        return f

    pdofs = Wp.dofmap.list[parent_ids, 0]
    sdofs = Ws.dofmap.list[np.arange(parent_ids.size), 0]

    # (p0) the committed transfer.
    p0 = parent_at_mean()
    _sub_to_parent(p0, sig_sub, parent_ids)
    rel_p0, m_par, m_sub = _moment_rel(msh, cell_tags, sub, p0, sig_sub)
    _say(
        f"MAT5-D' (p0) committed _sub_to_parent: parent = {m_par:.15e}  sub = {m_sub:.15e}  "
        f"rel = {rel_p0:.3e}  (bound {D_BAND:.0e}; predicted 0.0)"
    )

    # (n1) no transfer.
    n1 = parent_at_mean()
    rel_n1, _, _ = _moment_rel(msh, cell_tags, sub, n1, sig_sub)
    _say(
        f"MAT5-D' (n1) no transfer: rel = {rel_n1:.4e}  (predicted {N1_PREDICTED:.2e}; "
        f"asserted > {D_BAND:.0e})"
    )

    # (n2) shifted by one sub dof, per rank.
    n2 = parent_at_mean()
    n2.x.array[pdofs] = sig_sub.x.array[np.roll(sdofs, 1)]
    n2.x.scatter_forward()
    rel_n2, _, _ = _moment_rel(msh, cell_tags, sub, n2, sig_sub)
    _say(
        f"MAT5-D' (n2) shifted transfer (np.roll by 1 on each rank's local+ghost list): "
        f"rel = {rel_n2:.4e}  (size not predicted; asserted > {D_BAND:.0e})"
    )

    # (n3) reversed, per rank.
    n3 = parent_at_mean()
    n3.x.array[pdofs] = sig_sub.x.array[sdofs[::-1]]
    n3.x.scatter_forward()
    rel_n3, _, _ = _moment_rel(msh, cell_tags, sub, n3, sig_sub)
    _say(
        f"MAT5-D' (n3) reversed transfer (sdofs[::-1] on each rank's local+ghost list): "
        f"rel = {rel_n3:.4e}  (size not predicted, expected above (n2); asserted > "
        f"{D_BAND:.0e})"
    )

    # (n4) ghost-blind: owned sub cells only, no scatter_forward — printed only.
    n4 = parent_at_mean()
    n_owned = simap.size_local
    n4.x.array[pdofs[:n_owned]] = sig_sub.x.array[sdofs[:n_owned]]
    rel_n4, _, _ = _moment_rel(msh, cell_tags, sub, n4, sig_sub)
    _say(
        f"MAT5-D' (n4) ghost-blind transfer (owned cells only, no scatter_forward): "
        f"rel = {rel_n4:.3e}  (predicted 0 at -n 2 — not detected; printed, not asserted)"
    )

    assert (
        n_parent == N_PARENT and n_sub == N_SUB
    ), f"fixture moved: {n_parent} / {n_sub} cells, row {N_PARENT} / {N_SUB}"
    assert rel_p0 <= D_BAND, f"(p0) committed transfer {rel_p0:.3e} > {D_BAND:.0e}"
    assert rel_n1 > D_BAND, f"(n1) no transfer not seen: {rel_n1:.3e}"
    assert rel_n2 > D_BAND, f"(n2) shifted transfer not seen: {rel_n2:.3e}"
    assert rel_n3 > D_BAND, f"(n3) reversed transfer not seen: {rel_n3:.3e}"
