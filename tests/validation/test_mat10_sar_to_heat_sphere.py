"""`MAT-10` step 1 — `MAT-4`'s FEM SAR heats a Robin-cooled ball (§9 item 88).

The lossy sphere of `MAT-4` step 1 (R = 0.01 m, σ 0.57 S/m, εᵣ 78, 64 MHz,
h_sphere = R/10, h_far = R/5, complex build) is solved once; the heat source
q_h is the DG0 cell average of G·½σ|E|² (G = (5 kV/m / E₀)² = 2.5e7, applied
by linearity to q only), transferred onto a ``create_submesh`` of cell tag 1,
and drives item 85's ``solve_steady_conduction`` with h = 100 W/m²/K
(Bi = hR/k = 5/3, `MAT-9`'s) on every exterior facet of the submesh.

Closed form: θ(r) = `robin_sphere_steady_temperature(r, radius=R, k=0.6,
q=q_cf, h=100, t_inf=20)` − 20 with q_cf = G·ρ·SAR_cf (2002.1 W/m³):
θ(0) = 0.12235 K, θ(R) = 0.06674 K.  Errors are taken on θ, never on T.

Anchors (asserted): (P) uniform q_cf on the same submesh, interior
(midpoint |x|² < 0.81R²) relative L2 on θ ≤ 1 % (`MAT-9` (i)); (A) θ_h[q_h]
against θ ≤ 11 % — composed: `MAT-4`'s 10 % mean-SAR band plus (P)'s 1 %,
never tightened after reading; (B) ∫_sub q_h dV = G·``dissipated_power_w``
over tag 1 on the parent ≤ 1e-10; (C) ∮h(T_h − T∞)dS = ∫q_h dV ≤ 1e-9.
Negative controls (asserted only as "misses (A)'s 11 %"): (a) q from Re E
only (predicted ≈ 79.1 % miss); (b) missing ½ (predicted ≈ 107.1 %).
"""

from __future__ import annotations

import time

import dolfinx
import numpy as np
import ufl
from dolfinx import default_scalar_type, fem
from dolfinx import mesh as dmesh
from mpi4py import MPI

from fem_em_solver.post.sar import uniform_sphere_sar_closed_form
from fem_em_solver.thermal.steady import solve_steady_conduction
from fem_em_solver.utils.analytical import robin_sphere_steady_temperature
from fem_em_solver.utils.constants import EPSILON_0
from tests.complex_mode import complex_only
from tests.validation.test_lossy_sphere_sar import (
    E0,
    EPSILON_R_SPHERE,
    OMEGA,
    RHO_KG_M3,
    SIGMA_HIGH,
    SPHERE_RADIUS,
    SPHERE_TAG,
    _solve_lossy_sphere,
)

R = SPHERE_RADIUS
K = 0.6
H = 100.0  # Bi = hR/k = 5/3, MAT-9's
T_INF = 20.0
G = (5.0e3 / E0) ** 2  # 2.5e7, linear scale on q only
ROBIN_TAG = 1
QDEG = {"quadrature_degree": 4}

P_BAND = 1.0e-2  # (P): MAT-9 (i)
A_BAND = 0.11  # (A): MAT-4's 10 % + (P)'s 1 %, composed; never tightened
B_BAND = 1.0e-10
C_BAND = 1.0e-9


def _say(msg):
    if MPI.COMM_WORLD.rank == 0:
        print(msg, flush=True)


def _allsum(comm, form):
    return comm.allreduce(fem.assemble_scalar(fem.form(form)), op=MPI.SUM)


def _sphere_submesh(msh, cell_tags):
    """Cell tag 1 of the parent as its own mesh (`TH-14` pattern)."""
    tdim = msh.topology.dim
    imap = msh.topology.index_map(tdim)
    n_parent = imap.size_local + imap.num_ghosts
    parent_values = np.full(n_parent, -1, dtype=np.int32)
    parent_values[cell_tags.indices] = cell_tags.values
    keep = np.flatnonzero(parent_values == SPHERE_TAG).astype(np.int32)
    sub, emap, _vmap, _gmap = dmesh.create_submesh(msh, tdim, keep)
    simap = sub.topology.index_map(tdim)
    n_sub = simap.size_local + simap.num_ghosts
    parent_ids = np.asarray(
        emap.sub_topology_to_topology(np.arange(n_sub, dtype=np.int32), False),
        dtype=np.int32,
    )
    return sub, parent_ids


def _robin_tags(sub):
    """Every exterior facet of the submesh (the sphere surface) -> ROBIN_TAG."""
    tdim = sub.topology.dim
    fdim = tdim - 1
    sub.topology.create_connectivity(fdim, tdim)
    ext = np.sort(np.asarray(dmesh.exterior_facet_indices(sub.topology), dtype=np.int32))
    return dmesh.meshtags(sub, fdim, ext, np.full(ext.size, ROBIN_TAG, dtype=np.int32))


def _cell_average_dg0(msh, cell_tags, integrand):
    """Parent DG0 cell average of ``integrand`` over tag 1: ∫_K f·conj(w) / |K|."""
    W = fem.functionspace(msh, ("DG", 0))
    w = ufl.TestFunction(W)
    dx1 = ufl.Measure("dx", domain=msh, subdomain_data=cell_tags)(SPHERE_TAG)
    one = fem.Constant(msh, default_scalar_type(1.0))

    def _vec(form):
        b = fem.assemble_vector(fem.form(form))
        b.scatter_reverse(dolfinx.la.InsertMode.add)
        b.scatter_forward()
        return b.array.copy()

    num = _vec(ufl.inner(integrand, w) * dx1)
    vol = _vec(ufl.inner(one, w) * dx1)
    avg = np.zeros_like(num)
    mask = np.abs(vol) > 0.0
    avg[mask] = num[mask] / vol[mask]
    return W, avg


def _to_sub(W_parent, avg_parent, sub, parent_ids):
    """Map a parent DG0 array to a sub DG0 Function through both dofmaps."""
    Ws = fem.functionspace(sub, ("DG", 0))
    q = fem.Function(Ws)
    pdofs = W_parent.dofmap.list[parent_ids, 0]
    sdofs = Ws.dofmap.list[np.arange(parent_ids.size), 0]
    q.x.array[sdofs] = avg_parent[pdofs]
    q.x.scatter_forward()
    return q


def _interior_indicator(sub):
    tdim = sub.topology.dim
    imap = sub.topology.index_map(tdim)
    n = imap.size_local + imap.num_ghosts
    mids = dmesh.compute_midpoints(sub, tdim, np.arange(n, dtype=np.int32))
    W = fem.functionspace(sub, ("DG", 0))
    chi = fem.Function(W)
    chi.x.array[:] = (np.sum(mids**2, axis=1) < 0.81 * R**2).astype(float)
    chi.x.scatter_forward()
    return chi


def _theta_cf(sub, q_cf):
    x = ufl.SpatialCoordinate(sub)
    r2 = ufl.inner(x, x)
    return q_cf * R / (3.0 * H) + q_cf * (R**2 - r2) / (6.0 * K)


def _interior_rel(sub, chi, diff, ref):
    """sqrt(∫χ|diff|² / ∫χ|ref|²), complex-safe (inner, np.real)."""
    comm = sub.comm
    dx = ufl.dx(domain=sub, metadata=QDEG)
    e = _allsum(comm, chi * ufl.inner(diff, diff) * dx)
    n = _allsum(comm, chi * ufl.inner(ref, ref) * dx)
    return float(np.sqrt(np.real(e) / np.real(n)))


@complex_only
def test_mat10_fem_sar_heats_robin_ball():
    comm = MPI.COMM_WORLD
    t0 = time.perf_counter()

    run = _solve_lossy_sphere(SIGMA_HIGH, R / 10.0, R / 5.0, return_fields=True)
    msh, cell_tags = run["mesh"], run["cell_tags"]
    E, sigma = run["e_complex"], run["sigma_field"]
    t_em = time.perf_counter() - t0

    cf = uniform_sphere_sar_closed_form(
        e0=E0,
        epsilon_r=EPSILON_R_SPHERE,
        sigma=SIGMA_HIGH,
        omega=OMEGA,
        rho=RHO_KG_M3,
        epsilon_0=EPSILON_0,
    )
    q_cf = G * RHO_KG_M3 * cf["sar_w_per_kg"]
    th0 = robin_sphere_steady_temperature(0.0, radius=R, k=K, q=q_cf, h=H, t_inf=T_INF) - T_INF
    thR = robin_sphere_steady_temperature(R, radius=R, k=K, q=q_cf, h=H, t_inf=T_INF) - T_INF
    _say(
        f"MAT10 q_cf = {q_cf:.6f} W/m^3  theta(0) = {th0:.6f} K  theta(R) = {thR:.6f} K "
        f"(predicted 2002.1 / 0.12235 / 0.06674)"
    )

    sub, parent_ids = _sphere_submesh(msh, cell_tags)
    ft = _robin_tags(sub)
    n_parent = msh.topology.index_map(3).size_global
    n_sub = sub.topology.index_map(3).size_global
    _say(f"MAT10 cells: parent = {n_parent}  sub = {n_sub}  (predicted 74020 / ≈ 2e4)")

    # q_h: DG0 cell average of G·½σ|E|² on the parent, transferred to the sub.
    q_expr = G * 0.5 * sigma * ufl.inner(E, E)
    Wp, q_par = _cell_average_dg0(msh, cell_tags, q_expr)
    q_h = _to_sub(Wp, q_par, sub, parent_ids)
    n_own = sub.topology.index_map(3).size_local
    qmax = comm.allreduce(float(np.max(np.real(q_h.x.array[:n_own]))), op=MPI.MAX)
    qmin = comm.allreduce(float(np.min(np.real(q_h.x.array[:n_own]))), op=MPI.MIN)
    _say(f"MAT10 q_h max = {qmax:.6f}  min = {qmin:.6f} W/m^3  (q_cf {q_cf:.6f})")

    # (B) transfer: ∫_sub q_h dV against G·dissipated_power_w on the parent.
    int_q_sub = float(np.real(_allsum(comm, q_h * ufl.dx(domain=sub))))
    p_parent = G * run["dissipated_power_w"]
    err_B = abs(int_q_sub - p_parent) / abs(p_parent)
    _say(
        f"MAT10 (B) int_sub q_h = {int_q_sub:.15e}  G*P_diss = {p_parent:.15e}  "
        f"rel = {err_B:.3e}  (bound {B_BAND:.0e}, predicted ≈ 1e-15)"
    )

    chi = _interior_indicator(sub)
    theta = _theta_cf(sub, q_cf)

    # (P) uniform q_cf on the same submesh — first complex-build thermal solve.
    res_P = solve_steady_conduction(
        sub, ft, k=K, q=q_cf, h=H, t_inf=T_INF, robin_tag=ROBIN_TAG, degree=1
    )
    err_P = _interior_rel(sub, chi, (res_P.temperature - T_INF) - theta, theta)
    bal_P = abs(res_P.convected_power - res_P.source_power) / abs(res_P.source_power)
    _say(
        f"MAT10 (P) uniform q_cf: interior rel L2 = {err_P:.4e}  balance = {bal_P:.3e}  "
        f"(bound {P_BAND:.0e}, predicted ≈ 3.2e-3)"
    )

    # (A) FEM SAR as the source.
    res_A = solve_steady_conduction(
        sub, ft, k=K, q=q_h, h=H, t_inf=T_INF, robin_tag=ROBIN_TAG, degree=1
    )
    err_A = _interior_rel(sub, chi, (res_A.temperature - T_INF) - theta, theta)
    leg_offset = _interior_rel(sub, chi, res_A.temperature - res_P.temperature, theta)
    _say(
        f"MAT10 (A) theta_h[q_h] interior rel L2 = {err_A:.4e}  (bound {A_BAND}, "
        f"predicted ≈ 3.5–3.9e-2)"
    )
    _say(
        f"MAT10 (A) legs: offset theta_h[q_h - q_cf]/theta = {leg_offset:.4e}  "
        f"(mean-SAR offset {run['mean_sar'] / cf['sar_w_per_kg'] - 1.0:+.4e}); "
        f"(P) leg = {err_P:.4e}"
    )

    # (C) discrete balance.
    err_C = abs(res_A.convected_power - res_A.source_power) / abs(res_A.source_power)
    _say(
        f"MAT10 (C) conv = {res_A.convected_power:.15e}  src = {res_A.source_power:.15e}  "
        f"rel = {err_C:.3e}  (bound {C_BAND:.0e})"
    )

    # Surface-mean ratio (record).
    ds = ufl.Measure("ds", domain=sub, subdomain_data=ft)(ROBIN_TAG)
    area = float(np.real(_allsum(comm, fem.Constant(sub, default_scalar_type(1.0)) * ds)))
    th_s = float(np.real(_allsum(comm, (res_A.temperature - T_INF) * ds))) / area
    _say(f"MAT10 surface-mean theta_h[q_h] / theta(R) = {th_s / thR:.6f}  (record)")

    # Negative controls: thermal-only re-solves.
    q_re_expr = G * 0.5 * sigma * ufl.inner(ufl.real(E), ufl.real(E))
    _, q_re_par = _cell_average_dg0(msh, cell_tags, q_re_expr)
    q_re = _to_sub(Wp, q_re_par, sub, parent_ids)
    res_a = solve_steady_conduction(
        sub, ft, k=K, q=q_re, h=H, t_inf=T_INF, robin_tag=ROBIN_TAG, degree=1
    )
    err_a = _interior_rel(sub, chi, (res_a.temperature - T_INF) - theta, theta)
    q_two = fem.Function(q_h.function_space)
    q_two.x.array[:] = 2.0 * q_h.x.array
    q_two.x.scatter_forward()
    res_b = solve_steady_conduction(
        sub, ft, k=K, q=q_two, h=H, t_inf=T_INF, robin_tag=ROBIN_TAG, degree=1
    )
    err_b = _interior_rel(sub, chi, (res_b.temperature - T_INF) - theta, theta)
    _say(f"MAT10 control (a) Re-E-only: interior rel L2 = {err_a:.4e}  (predicted ≈ 0.791)")
    _say(f"MAT10 control (b) missing 1/2: interior rel L2 = {err_b:.4e}  (predicted ≈ 1.071)")
    _say(f"MAT10 elapsed: EM {t_em:.1f} s, total {time.perf_counter() - t0:.1f} s")

    assert err_B <= B_BAND, f"(B) transfer {err_B:.3e} > {B_BAND:.0e}"
    assert err_P <= P_BAND, f"(P) interior {err_P:.4e} > 1 %"
    assert err_C <= C_BAND, f"(C) balance {err_C:.3e} > {C_BAND:.0e}"
    assert err_A <= A_BAND, f"(A) {err_A:.4e} > 11 % (composed band)"
    assert err_a > A_BAND, f"control (a) Re-E-only {err_a:.4e} did not miss 11 %"
    assert err_b > A_BAND, f"control (b) missing 1/2 {err_b:.4e} did not miss 11 %"
