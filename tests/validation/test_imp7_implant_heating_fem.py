"""`IMP-7` step 1 (§9 item 594, §10 I38): P1 FEM of the implant-heated gel
ball against step 0's closed form, on the metal core and on the adiabatic hole.

Brief: docs/planning/chunks/IMP-7.md, "§9 item 594 — 2026-10-22" and the
2026-10-22 03:00 review ruling.  Fixture (step 0's, rule (f)): a = 5 mm,
b = 4a, gel k₂ 0.6, q_far 2000, h 100, T∞ 20, core k₁ 20; source = the β = 1
quasi-static dipole SAR of the PEC sphere, shell only.

Asserted (both fixtures): (G) census / conformity, cell guard, (S3) balance
≤ 1e-9 every rung, (S1) LSQ P1 whole-domain L2 rate ≥ 1.7 over three rungs,
(S2) finest interior rel L2 ≤ 1e-2; (Src′) / (Th′) ≤ `B1_BAND`.  Negative
controls (rule (e), predicted): isotropic source and β = 0 source, asserted
only as interior miss > 1e-2.  Printed only: the l = 2 readout at r₀ = 1.1a.

Step 2 (§9 item 630, §10 I40, prefix `IMP7-2`): the same ladder at P2 through
an additive `degree` keyword (default 1, so step 1 is unchanged); DG3 source,
P3 comparand, quadrature degree 6.  Asserted: S2_BAND / S3_BAND and the β = 0
control; printed only: the P2 rate and the l = 2 readout.
"""

from __future__ import annotations

import time

import numpy as np
import ufl
from dolfinx import default_scalar_type, fem, mesh as dmesh
from mpi4py import MPI

from fem_em_solver.io.mesh import MeshGenerator
from fem_em_solver.post.evaluation import evaluate_vector_field_parallel
from fem_em_solver.thermal.steady import solve_steady_conduction
from fem_em_solver.utils.analytical import (
    pec_sphere_implant_coefficients,
    pec_sphere_implant_radial,
    pec_sphere_implant_source,
    pec_sphere_implant_steady_temperature,
)
from tests.validation.test_imp1_pec_sphere_lossy_host import SPHERE_RADIUS
from tests.validation.test_imp4_pec_cylinder_closed_form import B1_BAND
from tests.validation.test_imp7_implant_heating_closed_form import (
    H,
    K_CORE,
    K_GEL,
    Q_FAR,
    T_INF,
)
from tests.validation.test_mat18_two_layer_sphere_fem import (
    CELL_GUARD,
    QDEG,
    _allsum,
    _dg0,
    _indicator,
)

# MAT-18 step 2(a) bands — literals in `assert` lines, restated by value, not widened.
S3_BAND = 1e-9  # balance: test_mat18_two_layer_sphere_fem.py:226 (MAT-9 :139)
S1_RATE = 1.7  # LSQ P1 L2 rate: test_mat18_two_layer_sphere_fem.py:253 (MAT-9 :154)
S2_BAND = 1e-2  # interior rel L2: test_mat18_two_layer_sphere_fem.py:254 (MAT-9 :156)
G_AREA_BAND = 1e-2  # (G) area: test_mat18_two_layer_sphere_fem.py:198

A = SPHERE_RADIUS  # 0.005
B = 4.0 * A
BETA = 1.0
ROBIN_TAG = 1
LADDER = (1.414e-3, 1.0e-3, 0.7071e-3)  # rule (f), ratio √2
CONTROL_RES = 1.0e-3
AREA_B = 4.0 * np.pi * B**2
AREA_A = 4.0 * np.pi * A**2
FIX = dict(radius=A, outer_radius=B, k_gel=K_GEL, q_far=Q_FAR, h=H, beta=BETA)
R0_READ = 1.1 * A

# Step 2 (§9 item 630, §10 I40): additive degree route.  The step-1 tests keep
# degree 1 through every default below.  At P2 the source is DG3 (interpolation
# error one order above the solution's), the comparand is P3 and the error
# forms use quadrature degree 6 (`MAT-21` step 3's rule), so the comparand
# floor stays below the P2 error.
DEG_ROUTE = {
    1: dict(q_dg=2, comparand=2, qmeta=QDEG),
    2: dict(q_dg=3, comparand=3, qmeta={"quadrature_degree": 6}),
}

_MESHES: dict = {}


def _say(msg):
    if MPI.COMM_WORLD.rank == 0:
        print(msg, flush=True)


def _ball(res):
    if res not in _MESHES:
        t0 = time.perf_counter()
        out = MeshGenerator.two_layer_ball_domain(
            inner_radius=A, outer_radius=B, resolution=res
        )
        _MESHES[res] = (*out, time.perf_counter() - t0)
    return _MESHES[res]


def _coeffs(k_core):
    return pec_sphere_implant_coefficients(k_core=k_core, **FIX)


# ---- unmasked shell continuations, x of shape (3, N) -------------------------


def _q_shell(x, beta=BETA, drop_l2=False):
    r = np.sqrt(np.sum(x * x, axis=0))
    u = (A / r) ** 3
    mu = x[2] / r
    p2 = 0.5 * (3.0 * mu * mu - 1.0)
    l2 = 0.0 if drop_l2 else (4.0 * beta * u + 2.0 * beta**2 * u**2) * p2
    return Q_FAR * (1.0 + 2.0 * beta**2 * u**2 + l2)


def _theta_shell(x, coeffs):
    r = np.sqrt(np.sum(x * x, axis=0))
    mu = x[2] / r
    p2 = 0.5 * (3.0 * mu * mu - 1.0)
    return (pec_sphere_implant_radial(r, coeffs, 0)
            + pec_sphere_implant_radial(r, coeffs, 2) * p2)


def _theta_metal(x):
    return pec_sphere_implant_steady_temperature(
        x.T, k_core=K_CORE, t_inf=T_INF, **FIX
    ) - T_INF


# ---- helpers (MAT-18's forms, parametrised by a and b) -----------------------


def _source_dg2(msh, cells, fn, dg=2):
    W = fem.functionspace(msh, ("DG", dg))
    q = fem.Function(W)
    q.x.array[:] = 0.0
    if cells is None:
        q.interpolate(fn)
    else:
        q.interpolate(fn, cells0=np.asarray(cells, dtype=np.int32))
    q.x.scatter_forward()
    return q


def _p2(msh, fn, p=2):
    V2 = fem.functionspace(msh, ("Lagrange", p))
    th = fem.Function(V2)
    th.interpolate(fn)
    th.x.scatter_forward()
    return th


def _errors(msh, chi_sh, T_h, theta, b, qmeta=QDEG):
    """(whole, interior |x|² < 0.81b², shell-only or None) relative L2 on θ."""
    comm = msh.comm
    dx = ufl.dx(domain=msh, metadata=qmeta)
    d = (T_h - T_INF) - theta
    chi_in = _indicator(msh, lambda m: np.sum(m**2, axis=1) < 0.81 * b**2)
    out = [np.sqrt(_allsum(comm, d**2 * dx) / _allsum(comm, theta**2 * dx))]
    for w in (chi_in, chi_sh):
        if w is None:
            out.append(float("nan"))
        else:
            out.append(np.sqrt(_allsum(comm, w * d**2 * dx)
                               / _allsum(comm, w * theta**2 * dx)))
    return tuple(float(v) for v in out)


def _vertex_radii(msh, cells):
    tdim = msh.topology.dim
    imap = msh.topology.index_map(tdim)
    gdm = np.asarray(msh.geometry.dofmap).reshape(imap.size_local + imap.num_ghosts, -1)
    if cells.size == 0:
        return np.empty(0)
    return np.linalg.norm(msh.geometry.x[gdm[cells]], axis=2).ravel()


def _rmin_rmax(comm, r):
    lo = comm.allreduce(float(r.min()) if r.size else np.inf, op=MPI.MIN)
    hi = comm.allreduce(float(r.max()) if r.size else 0.0, op=MPI.MAX)
    return lo, hi


def _census(msh, ct, ft, a, b):
    """(G): (n_core, n_shell, n_global, core_rmax, shell_rmin, shell_rmax,
    n_nonexterior_tag1_facets, area)."""
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
    _, core_rmax = _rmin_rmax(comm, _vertex_radii(msh, c1))
    shell_rmin, shell_rmax = _rmin_rmax(comm, _vertex_radii(msh, c2))
    msh.topology.create_connectivity(fdim, tdim)
    ext = set(dmesh.exterior_facet_indices(msh.topology).tolist())
    f1 = ft.find(ROBIN_TAG)
    bad = comm.allreduce(sum(1 for f in f1.tolist() if f not in ext), op=MPI.SUM)
    ds = ufl.Measure("ds", domain=msh, subdomain_data=ft)
    area = _allsum(comm, fem.Constant(msh, default_scalar_type(1.0)) * ds(ROBIN_TAG))
    return n1, n2, nglob, core_rmax, shell_rmin, shell_rmax, bad, area


def _check_g(msh, ct, ft, res):
    n1, n2, nglob, crmax, srmin, srmax, bad, area = _census(msh, ct, ft, A, B)
    _say(
        f"IMP7-1 G res={res} cells={nglob} core={n1} shell={n2} "
        f"core_rmax/a-1={crmax / A - 1:.3e} shell_rmin/a-1={srmin / A - 1:.3e} "
        f"shell_rmax/b-1={srmax / B - 1:.3e} nonexterior_tag1={bad} "
        f"area/4πb²-1={area / AREA_B - 1:.4e}"
    )
    assert n1 + n2 == nglob, f"(G) census {n1}+{n2} != {nglob}"
    assert crmax <= A * (1 + 1e-9), f"(G) core vertex at r={crmax!r} > a"
    assert srmin >= A * (1 - 1e-9), f"(G) shell vertex at r={srmin!r} < a"
    assert srmax <= B * (1 + 1e-9), f"(G) shell vertex at r={srmax!r} > b"
    assert bad == 0, f"(G) {bad} tag-1 facets are not exterior"
    assert abs(area / AREA_B - 1) <= G_AREA_BAND, f"(G) area {area / AREA_B - 1:.3e}"
    return n1, n2, nglob


def _shell_submesh(msh, ct):
    """Cell tag 2 of the parent as its own mesh (`MAT-10` `_sphere_submesh`)."""
    tdim = msh.topology.dim
    imap = msh.topology.index_map(tdim)
    n_parent = imap.size_local + imap.num_ghosts
    vals = np.full(n_parent, -1, dtype=np.int32)
    vals[ct.indices] = ct.values
    keep = np.flatnonzero(vals == 2).astype(np.int32)
    sub = dmesh.create_submesh(msh, tdim, keep)[0]
    return sub


def _robin_outer(sub, a, b):
    """Exterior facets with midpoint |x| > (a + b)/2 -> ROBIN_TAG; the r = a
    facets -> 2 (natural, zero flux; tagged only to measure their area)."""
    tdim = sub.topology.dim
    fdim = tdim - 1
    sub.topology.create_connectivity(fdim, tdim)
    ext = np.sort(np.asarray(dmesh.exterior_facet_indices(sub.topology), dtype=np.int32))
    mids = dmesh.compute_midpoints(sub, fdim, ext)
    outer = np.linalg.norm(mids, axis=1) > 0.5 * (a + b)
    tags = np.where(outer, ROBIN_TAG, 2).astype(np.int32)
    return dmesh.meshtags(sub, fdim, ext, tags)


def _solve(msh, ft, k, q, degree=1):
    t0 = time.perf_counter()
    sol = solve_steady_conduction(
        msh, ft, k=k, q=q, h=H, t_inf=T_INF, robin_tag=ROBIN_TAG, degree=degree
    )
    return sol, time.perf_counter() - t0


def _volume(msh):
    return _allsum(msh.comm, fem.Constant(msh, default_scalar_type(1.0)) * ufl.dx(domain=msh))


def _l2_readout(T_h, label, pred):
    """R₂,h(r₀) by Legendre projection, plus θ_h at pole / equator; printed only."""
    comm = MPI.COMM_WORLD
    mu, w = np.polynomial.legendre.leggauss(16)
    phi = 2.0 * np.pi * np.arange(8) / 8.0
    M, P = np.meshgrid(mu, phi, indexing="ij")
    s = np.sqrt(1.0 - M**2)
    pts = np.column_stack([
        (R0_READ * s * np.cos(P)).ravel(), (R0_READ * s * np.sin(P)).ravel(),
        (R0_READ * M).ravel(),
    ])
    pts = np.vstack([pts, [[0.0, 0.0, R0_READ], [R0_READ, 0.0, 0.0]]])
    vals, ok = evaluate_vector_field_parallel(T_h, pts)
    if comm.rank == 0:
        v = np.real(np.ravel(vals[:, 0]) if np.ndim(vals) == 2 else np.ravel(vals)) - T_INF
        ok = np.asarray(ok, dtype=bool)
        ring = v[:-2].reshape(16, 8)
        ok_ring = ok[:-2].reshape(16, 8)
        p2 = 0.5 * (3.0 * mu**2 - 1.0)
        r2h = 2.5 * float(np.sum(w * p2 * ring.mean(axis=1)))
        print(
            f"IMP7-1 l=2 readout [{label}] r0=1.1a: R2_h={r2h:.4e} (pred {pred[0]:.4e}) "
            f"pole={v[-2]:.6f} (pred {pred[1]:.6f}) equator={v[-1]:.6f} "
            f"(pred {pred[2]:.6f}) K; all_valid={bool(ok_ring.all() and ok[-2:].all())} "
            f"— printed, no band",
            flush=True,
        )


def _rate_gates(rows, label):
    h = np.array([r[2] for r in rows])
    e = np.array([r[3] for r in rows])
    for i in range(len(rows) - 1):
        _say(
            f"IMP7-1 [{label}] pairwise rate {LADDER[i]}->{LADDER[i + 1]}: "
            f"{np.log(e[i] / e[i + 1]) / np.log(h[i] / h[i + 1]):.3f}"
        )
    rate = float(np.polyfit(np.log(h), np.log(e), 1)[0])
    _say(f"IMP7-1 [{label}] (S1) LSQ P1 whole L2 rate = {rate:.3f} "
         f"(predicted 1.8–2.1, band ≥ {S1_RATE})")
    _say(f"IMP7-1 [{label}] (S2) finest interior rel L2 = {rows[-1][4]:.4e} "
         f"(predicted 5e-4..5e-3, band ≤ {S2_BAND})")
    assert rate >= S1_RATE, f"(S1) [{label}] P1 rate {rate:.3f} < {S1_RATE}"
    assert rows[-1][4] <= S2_BAND, f"(S2) [{label}] interior {rows[-1][4]:.4e} > {S2_BAND}"


# ---- tests --------------------------------------------------------------------


def test_imp7_1_source_guard():
    """(Src′) / (Th′): the unmasked continuations equal step 0's gated functions."""
    rng = np.random.default_rng(594)
    n = 512
    r = rng.uniform(1.001 * A, B, n)
    mu = rng.uniform(-1.0, 1.0, n)
    phi = rng.uniform(0.0, 2.0 * np.pi, n)
    s = np.sqrt(1.0 - mu * mu)
    pts = np.column_stack([r * s * np.cos(phi), r * s * np.sin(phi), r * mu])
    q_ref = pec_sphere_implant_source(pts, radius=A, q_far=Q_FAR, beta=BETA)
    q_new = _q_shell(pts.T)
    e_src = float(np.max(np.abs(q_new - q_ref)) / np.max(np.abs(q_ref)))
    c0 = _coeffs(0.0)
    th_ref = pec_sphere_implant_steady_temperature(pts, k_core=0.0, t_inf=T_INF, **FIX) - T_INF
    th_new = _theta_shell(pts.T, c0)
    e_th = float(np.max(np.abs(th_new - th_ref)) / np.max(np.abs(th_ref)))
    _say(f"IMP7-1 (Src') {e_src:.3e}  (Th') {e_th:.3e}  (band B1_BAND {B1_BAND:.0e})")
    assert e_src <= B1_BAND, f"(Src') {e_src:.3e}"
    assert e_th <= B1_BAND, f"(Th') {e_th:.3e}"


def _metal_rung(res, q_fn=_q_shell, label="metal", degree=1):
    comm = MPI.COMM_WORLD
    rt = DEG_ROUTE[degree]
    msh, ct, ft, t_build = _ball(res)
    k_f = _dg0(msh, ct, K_CORE, K_GEL)
    q_f = _source_dg2(msh, ct.find(2), q_fn, dg=rt["q_dg"])
    sol, t_solve = _solve(msh, ft, k_f, q_f, degree=degree)
    theta = _p2(msh, _theta_metal, p=rt["comparand"])
    chi_sh = _dg0(msh, ct, 0.0, 1.0)
    e_all, e_in, e_sh = _errors(msh, chi_sh, sol.temperature, theta, B,
                                qmeta=rt["qmeta"])
    nglob = msh.topology.index_map(3).size_global
    vol = _volume(msh)
    vcore = _allsum(comm, _dg0(msh, ct, 1.0, 0.0) * ufl.dx(domain=msh))
    h_eff = (vol / nglob) ** (1.0 / 3.0)
    bal = abs(sol.convected_power - sol.source_power) / sol.source_power
    return msh, sol, t_build, t_solve, (res, nglob, h_eff, e_all, e_in, e_sh), bal, vcore


def test_imp7_1_metal_core_ladder():
    P = _coeffs(K_CORE)["power"]
    rows = []
    for res in LADDER:
        msh, ct, ft, _ = _ball(res)
        n1, n2, nglob = _check_g(msh, ct, ft, res)
        if res == LADDER[-1]:
            _say(f"IMP7-1 cell guard: {nglob} cells at res={res} (guard {CELL_GUARD})")
            assert nglob <= CELL_GUARD, f"cell guard: {nglob} > {CELL_GUARD}; stop"
        msh, sol, t_b, t_s, row, bal, vcore = _metal_rung(res)
        rows.append(row)
        _say(
            f"IMP7-1 [metal] P1 res={res} cells={row[1]} h_eff={row[2]:.5e} "
            f"build={t_b:.1f}s solve={t_s:.1f}s whole_rel_L2={row[3]:.4e} "
            f"interior_rel_L2={row[4]:.4e} shell_rel_L2={row[5]:.4e} "
            f"Vcore_h/Vcore-1={vcore / (4 / 3 * np.pi * A**3) - 1:.4e} "
            f"int_q_h/P-1={sol.source_power / P - 1:.4e} balance={bal:.3e} "
            f"src={sol.source_power:.6e}W conv={sol.convected_power:.6e}W"
        )
        assert bal <= S3_BAND, f"(S3) balance {bal:.3e} > {S3_BAND} at res={res}"
        if res == LADDER[-1]:
            _l2_readout(sol.temperature, "metal", (1.1966e-2, 0.361610, 0.343661))
    _rate_gates(rows, "metal")


def test_imp7_1_metal_core_controls():
    for lab, fn, pred in (
        ("c-iso", lambda x: _q_shell(x, drop_l2=True), 2.730e-2),
        ("c-beta0", lambda x: _q_shell(x, beta=0.0), 5.671e-2),
    ):
        _, _, _, t_s, row, bal, _ = _metal_rung(CONTROL_RES, q_fn=fn)
        _say(
            f"IMP7-1 [metal] ({lab}) res={CONTROL_RES}: interior miss measured "
            f"{row[4]:.4e} vs predicted {pred:.3e} (asserted > {S2_BAND}); "
            f"whole {row[3]:.4e}; solve={t_s:.1f}s"
        )
        assert row[4] > S2_BAND, f"({lab}) control did not fire: {row[4]:.4e}"


def _hole_rung(res, q_fn=_q_shell, check=False, degree=1):
    comm = MPI.COMM_WORLD
    rt = DEG_ROUTE[degree]
    msh, ct, ft, t_build = _ball(res)
    c2 = ct.find(2)
    n2 = comm.allreduce(int(np.sum(c2 < msh.topology.index_map(3).size_local)), op=MPI.SUM)
    sub = _shell_submesh(msh, ct)
    sft = _robin_outer(sub, A, B)
    nsub = sub.topology.index_map(3).size_global
    ds = ufl.Measure("ds", domain=sub, subdomain_data=sft)
    one = fem.Constant(sub, default_scalar_type(1.0))
    area_b = _allsum(comm, one * ds(ROBIN_TAG))
    area_a = _allsum(comm, one * ds(2))
    if check:
        srmin, srmax = _rmin_rmax(
            comm, _vertex_radii(sub, np.arange(sub.topology.index_map(3).size_local,
                                               dtype=np.int32)))
        _say(
            f"IMP7-1 G-hole res={res} sub_cells={nsub} parent_tag2={n2} "
            f"shell_rmin/a-1={srmin / A - 1:.3e} shell_rmax/b-1={srmax / B - 1:.3e} "
            f"robin_area/4πb²-1={area_b / AREA_B - 1:.4e} "
            f"r=a_area/4πa²-1={area_a / AREA_A - 1:.4e} (printed)"
        )
        assert nsub == n2, f"(G) hole cells {nsub} != parent tag-2 {n2}"
        assert srmin >= A * (1 - 1e-9), f"(G) hole vertex at r={srmin!r} < a"
        assert srmax <= B * (1 + 1e-9), f"(G) hole vertex at r={srmax!r} > b"
        assert abs(area_b / AREA_B - 1) <= G_AREA_BAND, (
            f"(G-hole) Robin area {area_b / AREA_B - 1:.3e}")
    q_f = _source_dg2(sub, None, q_fn, dg=rt["q_dg"])
    sol, t_solve = _solve(sub, sft, K_GEL, q_f, degree=degree)
    c0 = _coeffs(0.0)
    theta = _p2(sub, lambda x: _theta_shell(x, c0), p=rt["comparand"])
    nan_loc = int(np.count_nonzero(~np.isfinite(theta.x.array)))
    nan = comm.allreduce(nan_loc, op=MPI.SUM)
    assert nan == 0, f"hole comparand holds {nan} NaN dofs"
    e_all, e_in, _ = _errors(sub, None, sol.temperature, theta, B, qmeta=rt["qmeta"])
    vol = _volume(sub)
    h_eff = (vol / nsub) ** (1.0 / 3.0)
    bal = abs(sol.convected_power - sol.source_power) / sol.source_power
    return sol, t_build, t_solve, (res, nsub, h_eff, e_all, e_in), bal, vol


def test_imp7_1_hole_ladder():
    P = _coeffs(0.0)["power"]
    v_shell = 4.0 / 3.0 * np.pi * (B**3 - A**3)
    rows = []
    for res in LADDER:
        msh, ct, ft, _ = _ball(res)
        nglob = msh.topology.index_map(3).size_global
        if res == LADDER[-1]:
            _say(f"IMP7-1 cell guard (parent): {nglob} cells at res={res} (guard {CELL_GUARD})")
            assert nglob <= CELL_GUARD, f"cell guard: {nglob} > {CELL_GUARD}; stop"
        sol, t_b, t_s, row, bal, vol = _hole_rung(res, check=True)
        rows.append(row)
        _say(
            f"IMP7-1 [hole] P1 res={res} cells={row[1]} h_eff={row[2]:.5e} "
            f"build={t_b:.1f}s solve={t_s:.1f}s whole_rel_L2={row[3]:.4e} "
            f"interior_rel_L2={row[4]:.4e} Vshell_h/Vshell-1={vol / v_shell - 1:.4e} "
            f"int_q_h/P-1={sol.source_power / P - 1:.4e} balance={bal:.3e} "
            f"src={sol.source_power:.6e}W conv={sol.convected_power:.6e}W"
        )
        assert bal <= S3_BAND, f"(S3) hole balance {bal:.3e} > {S3_BAND} at res={res}"
        if res == LADDER[-1]:
            _l2_readout(sol.temperature, "hole", (4.4793e-2, 0.394437, 0.327248))
    _rate_gates(rows, "hole")


def test_imp7_1_hole_controls():
    for lab, fn, pred in (
        ("c-iso", lambda x: _q_shell(x, drop_l2=True), 3.603e-2),
        ("c-beta0", lambda x: _q_shell(x, beta=0.0), 6.195e-2),
    ):
        _, _, t_s, row, bal, _ = _hole_rung(CONTROL_RES, q_fn=fn)
        _say(
            f"IMP7-1 [hole] ({lab}) res={CONTROL_RES}: interior miss measured "
            f"{row[4]:.4e} vs predicted {pred:.3e} (asserted > {S2_BAND}); "
            f"whole {row[3]:.4e}; solve={t_s:.1f}s"
        )
        assert row[4] > S2_BAND, f"({lab}) control did not fire: {row[4]:.4e}"


# ---- step 2 (§9 item 630, §10 I40): P2 on both cases ---------------------------
# Asserted (order-free; S2_BAND / S3_BAND are this module's names): finest
# interior rel L2 ≤ S2_BAND, balance ≤ S3_BAND on every rung, and the β = 0
# control (at CONTROL_RES) missing S2_BAND.  Printed only, no band: the P2 L2
# rate (MAT-18 has no degree-2 rate gate) and the l = 2 readout at 1.1a beside
# step 1's P1 values and the series.

P1_READ = {"metal": (1.0174e-2, 0.359626, 0.344320),
           "hole": (4.3742e-2, 0.393152, 0.327785)}  # IMP-7.md:429 (step-1 W1:71 / W2:68)
SERIES_READ = {"metal": (1.1966e-2, 0.361610, 0.343661),
               "hole": (4.4793e-2, 0.394437, 0.327248)}  # IMP-7.md:367-368
P1_BETA0 = {"metal": 5.965e-2, "hole": 6.243e-2}  # IMP-7.md:429


def _p2_rate_print(rows, label):
    h = np.array([r[2] for r in rows])
    e = np.array([r[3] for r in rows])
    for i in range(len(rows) - 1):
        _say(f"IMP7-2 [{label}] P2 pairwise rate {LADDER[i]}->{LADDER[i + 1]}: "
             f"{np.log(e[i] / e[i + 1]) / np.log(h[i] / h[i + 1]):.3f}")
    rate = float(np.polyfit(np.log(h), np.log(e), 1)[0])
    _say(f"IMP7-2 [{label}] P2 LSQ whole L2 rate = {rate:.3f} (printed, no band)")


def _p2_case(label, rung):
    rows = []
    for res in LADDER:
        sol, t_b, t_s, row, bal = rung(res)
        rows.append(row)
        _say(f"IMP7-2 [{label}] P2 res={res} cells={row[1]} h_eff={row[2]:.5e} "
             f"build={t_b:.1f}s solve={t_s:.1f}s whole_rel_L2={row[3]:.4e} "
             f"interior_rel_L2={row[4]:.4e} balance={bal:.3e}")
        assert bal <= S3_BAND, f"(S3) [{label}] P2 balance {bal:.3e} > {S3_BAND} at res={res}"
        if res == LADDER[-1]:
            p1 = P1_READ[label]
            _say(f"IMP7-2 [{label}] P1 readout (step 1): R2_h={p1[0]:.4e} "
                 f"pole={p1[1]:.6f} equator={p1[2]:.6f} K")
            _l2_readout(sol.temperature, f"{label} P2", SERIES_READ[label])
    _p2_rate_print(rows, label)
    _say(f"IMP7-2 [{label}] P2 finest interior rel L2 = {rows[-1][4]:.4e} "
         f"(band <= S2_BAND {S2_BAND})")
    assert rows[-1][4] <= S2_BAND, f"[{label}] P2 interior {rows[-1][4]:.4e} > {S2_BAND}"
    _, _, t_s, row, bal = rung(CONTROL_RES, q_fn=lambda x: _q_shell(x, beta=0.0))
    _say(f"IMP7-2 [{label}] P2 (c-beta0) res={CONTROL_RES}: interior miss {row[4]:.4e} "
         f"(P1 read {P1_BETA0[label]:.3e}; asserted > {S2_BAND}); "
         f"whole {row[3]:.4e}; solve={t_s:.1f}s")
    assert row[4] > S2_BAND, f"[{label}] P2 (c-beta0) control did not fire: {row[4]:.4e}"


def test_imp7_2_metal_core_p2():
    def rung(res, q_fn=_q_shell):
        _, sol, t_b, t_s, row, bal, _ = _metal_rung(res, q_fn=q_fn, degree=2)
        return sol, t_b, t_s, row, bal
    _p2_case("metal", rung)


def test_imp7_2_hole_p2():
    def rung(res, q_fn=_q_shell):
        sol, t_b, t_s, row, bal, _ = _hole_rung(res, q_fn=q_fn, degree=2)
        return sol, t_b, t_s, row, bal
    _p2_case("hole", rung)
