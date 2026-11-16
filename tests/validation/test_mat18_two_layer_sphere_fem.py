"""`MAT-18` step 2(a): conforming two-region ball + piecewise k, steady FEM.

Gel core (r < a = 0.1, k₁ 0.6, q₁ 2000) inside a container-wall shell
(a < r < b = 0.11, k₂ 0.19, q₂ 0), Robin h 10 / T∞ 20 on r = b, P1.
Comparand: step 1's ``two_layer_sphere_steady_temperature`` (closed form,
no series), interpolated into P2.  Bands are `MAT-9` step 1's literals
(``test_mat9_steady_conduction_sphere.py:139, 154, 156``).

Anchors (asserted): (G) conformity / census on every rung; (S1) whole-ball
P1 L2 rate ≥ 1.7; (S2) interior rel L2 ≤ 1 % at 0.003536; (S3) balance
≤ 1e-9 on every rung; (D) DG0 k ≡ K1 vs scalar K1 ≤ 1e-13; (E) equal
properties vs `robin_sphere_steady_temperature(radius=b)` interior ≤ 1 %.
Negative controls (rule (e), predicted, asserted only as missing band):
(c1) uniform k, predicted interior miss 0.197, asserted > 1e-2; (c2) source
in both regions, predicted balance miss 0.331, asserted > 1e-3.
"""

from __future__ import annotations

import time

import numpy as np
import pytest
import ufl
from dolfinx import default_scalar_type, fem, mesh as dmesh
from mpi4py import MPI

from fem_em_solver.io.mesh import MeshGenerator
from fem_em_solver.post.evaluation import evaluate_vector_field_parallel
from fem_em_solver.thermal.steady import solve_steady_conduction
from fem_em_solver.utils.analytical import (
    robin_sphere_steady_temperature,
    two_layer_sphere_steady_temperature,
)
from tests.validation.test_mat18_two_layer_sphere_closed_form import (
    A_,
    B_,
    H,
    K1,
    K2,
    Q1,
    Q2,
    T_INF,
)

ROBIN_TAG = 1
LADDER = (0.01, 0.007071, 0.005, 0.003536)
CELL_GUARD = 800_000
QDEG = {"quadrature_degree": 4}
V_CORE = 4.0 / 3.0 * np.pi * A_**3
AREA_B = 4.0 * np.pi * B_**2

_MESHES: dict = {}


def _say(msg):
    if MPI.COMM_WORLD.rank == 0:
        print(msg, flush=True)


def _ball(res):
    if res not in _MESHES:
        t0 = time.perf_counter()
        out = MeshGenerator.two_layer_ball_domain(
            inner_radius=A_, outer_radius=B_, resolution=res
        )
        _MESHES[res] = (*out, time.perf_counter() - t0)
    return _MESHES[res]


def _allsum(comm, form):
    return float(
        np.real(comm.allreduce(fem.assemble_scalar(fem.form(form)), op=MPI.SUM))
    )


def _dg0(msh, ct, v_core, v_shell):
    """DG0 coefficient set from the cell tags, then ghost-updated."""
    W = fem.functionspace(msh, ("DG", 0))
    f = fem.Function(W)
    f.x.array[:] = np.nan
    f.x.array[ct.find(1)] = v_core
    f.x.array[ct.find(2)] = v_shell
    f.x.scatter_forward()
    return f


def _indicator(msh, mask_fn):
    tdim = msh.topology.dim
    imap = msh.topology.index_map(tdim)
    n = imap.size_local + imap.num_ghosts
    mids = dmesh.compute_midpoints(msh, tdim, np.arange(n, dtype=np.int32))
    W = fem.functionspace(msh, ("DG", 0))
    chi = fem.Function(W)
    chi.x.array[:] = mask_fn(mids).astype(float)
    chi.x.scatter_forward()
    return chi


def _comparand(msh, fn):
    """θ_exact = T_exact − T∞ interpolated into P2."""
    V2 = fem.functionspace(msh, ("Lagrange", 2))
    th = fem.Function(V2)
    th.interpolate(lambda x: fn(np.sqrt(np.sum(x**2, axis=0))) - T_INF)
    th.x.scatter_forward()
    return th


def _two_layer(r):
    return two_layer_sphere_steady_temperature(
        r, a=A_, b=B_, k1=K1, k2=K2, q1=Q1, q2=Q2, h=H, t_inf=T_INF
    )


def _errors(msh, ct, T_h, theta):
    """(whole, interior, shell-only) relative L2 errors on θ."""
    comm = msh.comm
    dx = ufl.dx(domain=msh, metadata=QDEG)
    d = (T_h - T_INF) - theta
    chi_in = _indicator(msh, lambda m: np.sum(m**2, axis=1) < 0.81 * B_**2)
    chi_sh = _dg0(msh, ct, 0.0, 1.0)
    out = []
    for w in (None, chi_in, chi_sh):
        if w is None:
            out.append(np.sqrt(_allsum(comm, d**2 * dx) / _allsum(comm, theta**2 * dx)))
        else:
            out.append(
                np.sqrt(_allsum(comm, w * d**2 * dx) / _allsum(comm, w * theta**2 * dx))
            )
    return tuple(float(v) for v in out)


def _census(msh, ct, ft):
    """(G): returns (n_core, n_shell, n_global, core_rmax, shell_rmin,
    shell_rmax, n_nonexterior_tag1_facets, area)."""
    comm = msh.comm
    tdim = msh.topology.dim
    fdim = tdim - 1
    imap = msh.topology.index_map(tdim)
    nloc = imap.size_local
    c1 = ct.find(1)
    c2 = ct.find(2)
    c1 = c1[c1 < nloc]
    c2 = c2[c2 < nloc]
    n1 = comm.allreduce(int(c1.size), op=MPI.SUM)
    n2 = comm.allreduce(int(c2.size), op=MPI.SUM)
    nglob = imap.size_global

    gdm = msh.geometry.dofmap
    gdm = np.asarray(gdm).reshape(imap.size_local + imap.num_ghosts, -1)
    X = msh.geometry.x

    def rr(cells):
        if cells.size == 0:
            return np.empty(0)
        return np.linalg.norm(X[gdm[cells]], axis=2).ravel()

    r1 = rr(c1)
    r2 = rr(c2)
    core_rmax = comm.allreduce(float(r1.max()) if r1.size else 0.0, op=MPI.MAX)
    shell_rmin = comm.allreduce(float(r2.min()) if r2.size else np.inf, op=MPI.MIN)
    shell_rmax = comm.allreduce(float(r2.max()) if r2.size else 0.0, op=MPI.MAX)

    msh.topology.create_connectivity(fdim, tdim)
    ext = set(dmesh.exterior_facet_indices(msh.topology).tolist())
    f1 = ft.find(ROBIN_TAG)
    bad = sum(1 for f in f1.tolist() if f not in ext)
    bad = comm.allreduce(bad, op=MPI.SUM)

    ds = ufl.Measure("ds", domain=msh, subdomain_data=ft)
    area = _allsum(comm, fem.Constant(msh, default_scalar_type(1.0)) * ds(ROBIN_TAG))
    return n1, n2, nglob, core_rmax, shell_rmin, shell_rmax, bad, area


def _vcore(msh, ct):
    chi = _dg0(msh, ct, 1.0, 0.0)
    return _allsum(msh.comm, chi * ufl.dx(domain=msh))


def test_ladder_conformity_rate_interior_balance():
    comm = MPI.COMM_WORLD
    rows = []
    for res in LADDER:
        msh, ct, ft, t_build = _ball(res)
        n1, n2, nglob, crmax, srmin, srmax, bad, area = _census(msh, ct, ft)
        _say(
            f"MAT18 G res={res} cells={nglob} core={n1} shell={n2} "
            f"core_rmax/a-1={crmax / A_ - 1:.3e} shell_rmin/a-1={srmin / A_ - 1:.3e} "
            f"shell_rmax/b-1={srmax / B_ - 1:.3e} nonexterior_tag1={bad} "
            f"area/4πb²-1={area / AREA_B - 1:.4e}"
        )
        # (G) conformity and census.
        assert n1 + n2 == nglob, f"(G) census {n1}+{n2} != {nglob}"
        assert crmax <= A_ * (1 + 1e-9), f"(G) core vertex at r={crmax!r} > a"
        assert srmin >= A_ * (1 - 1e-9), f"(G) shell vertex at r={srmin!r} < a"
        assert srmax <= B_ * (1 + 1e-9), f"(G) shell vertex at r={srmax!r} > b"
        assert bad == 0, f"(G) {bad} tag-1 facets are not exterior"
        assert abs(area / AREA_B - 1) <= 1e-2, f"(G) area deficit {area / AREA_B - 1:.3e}"

        if res == LADDER[-1]:
            _say(f"MAT18 cell guard: {nglob} cells at res={res} (guard {CELL_GUARD})")
            assert nglob <= CELL_GUARD, f"cell guard: {nglob} > {CELL_GUARD}; stop"

        k_f = _dg0(msh, ct, K1, K2)
        q_f = _dg0(msh, ct, Q1, Q2)
        t0 = time.perf_counter()
        sol = solve_steady_conduction(
            msh, ft, k=k_f, q=q_f, h=H, t_inf=T_INF, robin_tag=ROBIN_TAG, degree=1
        )
        t_solve = time.perf_counter() - t0
        theta = _comparand(msh, _two_layer)
        e_all, e_in, e_sh = _errors(msh, ct, sol.temperature, theta)
        vol = _allsum(comm, fem.Constant(msh, default_scalar_type(1.0)) * ufl.dx(domain=msh))
        vc = _vcore(msh, ct)
        h_eff = (vol / nglob) ** (1.0 / 3.0)
        bal = abs(sol.convected_power - Q1 * vc) / (Q1 * vc)
        rows.append((res, nglob, h_eff, e_all, e_in))
        _say(
            f"MAT18 P1 res={res} cells={nglob} h_eff={h_eff:.5e} build={t_build:.1f}s "
            f"solve={t_solve:.1f}s whole_rel_L2={e_all:.4e} interior_rel_L2={e_in:.4e} "
            f"shell_rel_L2={e_sh:.4e} Vcore_h/Vcore-1={vc / V_CORE - 1:.4e} "
            f"area_deficit={area / AREA_B - 1:.4e} balance={bal:.3e} "
            f"src={sol.source_power:.6f}W conv={sol.convected_power:.6f}W"
        )
        # (S3) balance on every rung.
        assert bal <= 1e-9, f"(S3) balance {bal:.3e} > 1e-9 at res={res}"

        if res == LADDER[-1]:
            pts = np.array([[0.0, 0.0, 0.0], [A_, 0.0, 0.0], [0.999 * B_, 0.0, 0.0]])
            vals, ok = evaluate_vector_field_parallel(sol.temperature, pts)
            if comm.rank == 0:
                for lab, v, o, ref in zip(
                    ("0", "a", "0.999b"), np.ravel(vals[:, 0]) if vals.ndim == 2 else vals,
                    ok, (14.2550, 8.6994, 5.5096),
                ):
                    th = float(np.real(v)) - T_INF if o else float("nan")
                    print(
                        f"MAT18 theta_h(r={lab}) = {th:.4f} K (valid={bool(o)}); "
                        f"step-1 record at r={'b' if lab == '0.999b' else lab}: {ref} K",
                        flush=True,
                    )

    h = np.array([r[2] for r in rows])
    e = np.array([r[3] for r in rows])
    for i in range(len(rows) - 1):
        _say(
            f"MAT18 pairwise rate {LADDER[i]}->{LADDER[i + 1]}: "
            f"{np.log(e[i] / e[i + 1]) / np.log(h[i] / h[i + 1]):.3f}"
        )
    rate = float(np.polyfit(np.log(h), np.log(e), 1)[0])
    _say(f"MAT18 (S1) LSQ P1 whole-ball L2 rate = {rate:.3f} (predicted ≈ 2.0, band ≥ 1.7)")
    _say(f"MAT18 (S2) finest interior rel L2 = {rows[-1][4]:.4e} (predicted 4e-4..1.5e-3, band ≤ 1e-2)")
    assert rate >= 1.7, f"(S1) P1 rate {rate:.3f} < 1.7"
    assert rows[-1][4] <= 1e-2, f"(S2) interior {rows[-1][4]:.4e} > 1 %"


def test_d_coefficient_path():
    """(D) DG0 k ≡ K1 against scalar k = K1, same piecewise q, at 0.01."""
    msh, ct, ft, _ = _ball(0.01)
    comm = msh.comm
    q_f = _dg0(msh, ct, Q1, Q2)
    k_f = _dg0(msh, ct, K1, K1)
    s_f = solve_steady_conduction(
        msh, ft, k=k_f, q=q_f, h=H, t_inf=T_INF, robin_tag=ROBIN_TAG, degree=1
    )
    s_c = solve_steady_conduction(
        msh, ft, k=K1, q=q_f, h=H, t_inf=T_INF, robin_tag=ROBIN_TAG, degree=1
    )
    dx = ufl.dx(domain=msh, metadata=QDEG)
    d = np.sqrt(
        _allsum(comm, (s_f.temperature - s_c.temperature) ** 2 * dx)
        / _allsum(comm, (s_c.temperature - T_INF) ** 2 * dx)
    )
    _say(f"MAT18 (D) DG0-vs-scalar whole rel diff = {d:.3e} (predicted ≤ 1e-15, band ≤ 1e-13)")
    assert d <= 1e-13, f"(D) {d:.3e} > 1e-13"


def test_e_equal_properties():
    """(E) DG0 k ≡ K1, q ≡ Q1 in both regions vs the homogeneous Robin ball."""
    msh, ct, ft, _ = _ball(0.005)
    k_f = _dg0(msh, ct, K1, K1)
    q_f = _dg0(msh, ct, Q1, Q1)
    sol = solve_steady_conduction(
        msh, ft, k=k_f, q=q_f, h=H, t_inf=T_INF, robin_tag=ROBIN_TAG, degree=1
    )
    theta = _comparand(
        msh,
        lambda r: robin_sphere_steady_temperature(
            r, radius=B_, k=K1, q=Q1, h=H, t_inf=T_INF
        ),
    )
    e_all, e_in, e_sh = _errors(msh, ct, sol.temperature, theta)
    _say(
        f"MAT18 (E) equal properties res=0.005: interior={e_in:.4e} whole={e_all:.4e} "
        f"shell={e_sh:.4e} (predicted interior ≈ 8e-4, band ≤ 1e-2)"
    )
    assert e_in <= 1e-2, f"(E) interior {e_in:.4e} > 1 %"


def test_c1_uniform_k_control():
    """(c1) scalar k = K1 with piecewise q against the two-layer comparand."""
    msh, ct, ft, _ = _ball(0.005)
    q_f = _dg0(msh, ct, Q1, Q2)
    sol = solve_steady_conduction(
        msh, ft, k=K1, q=q_f, h=H, t_inf=T_INF, robin_tag=ROBIN_TAG, degree=1
    )
    theta = _comparand(msh, _two_layer)
    e_all, e_in, e_sh = _errors(msh, ct, sol.temperature, theta)
    _say(
        f"MAT18 (c1) uniform k: interior miss measured {e_in:.4f} vs predicted 0.197 "
        f"(asserted > 1e-2); whole {e_all:.4f}"
    )
    assert e_in > 1e-2, f"(c1) control did not fire: {e_in:.4e}"


def test_c2_source_both_regions_control():
    """(c2) q₁ in both regions, balance taken against q₁V_core,h."""
    msh, ct, ft, _ = _ball(0.005)
    k_f = _dg0(msh, ct, K1, K2)
    sol = solve_steady_conduction(
        msh, ft, k=k_f, q=Q1, h=H, t_inf=T_INF, robin_tag=ROBIN_TAG, degree=1
    )
    vc = _vcore(msh, ct)
    miss = abs(sol.convected_power - Q1 * vc) / (Q1 * vc)
    _say(
        f"MAT18 (c2) source in both regions: balance miss measured {miss:.4f} vs "
        f"predicted 0.331 (asserted > 1e-3)"
    )
    assert miss > 1e-3, f"(c2) control did not fire: {miss:.4e}"
