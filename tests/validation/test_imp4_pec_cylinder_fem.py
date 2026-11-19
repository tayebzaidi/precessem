"""`IMP-4` step 1 (§9 item 206, lane 3): FEM for the infinite PEC circular
cylinder in `IMP-1`'s lossy host at finite k, against step 0's closed form.

Fixture: `TH-15`'s hole route — the a = 1 mm cylinder cut from a coaxial
c = 10 mm host annulus of length L = 4 mm by
``MeshGenerator.cylindrical_domain(as_hole=True)`` (facet tags 1 outer wall,
2 cavity, 3 end annuli; h(ρ) = ρ/6).  The outer wall carries step 0's exact
finite-k tangential trace; the end planes carry the condition the exact
z-invariant field already satisfies (Jin §1.4.2 eqs (1.44)/(1.45), §11.2.1
eq (11.53)): **TM_z n = 0** has n̂ × E = 0 on z = const, so the ends are PEC
(pec tags 1, 2, 3); **TE_z n = 1** has n̂ × H = 0 there, so the ends are
natural / PMC (pec tags 1, 2).  Truncation and end effects are therefore
exact; only discretisation and faceting remain.

Gated (degree 2, each drive): (i) |β − 1| ≤ `TH15_BETA_BAND`, β the complex
LSQ amplitude of the closed-form scattered field on 24-point circles at
ρ = 1.2a, 1.5a, 2a (z = 0); (ii) the annulus (ρ < 1.5a) mean |E|² against
the closed form integrated by numpy quadrature over the same owned cells
within `MAT4_BAND`; (iii) cavity dofs zero.  W0: cell guard, tag areas,
(M′) numpy-vs-UFL quadrature identity, (S) the incident/scattered split.
Controls (W2): natural cavity — TM (pec 1, 3) misses (i) and (ii); TE
(pec 1) misses (i) only (Φ cannot discriminate in TE, see the brief).
"""

from __future__ import annotations

import time

import numpy as np
import pytest
import ufl
from dolfinx import fem
from mpi4py import MPI
from petsc4py import PETSc
from scipy import special

from fem_em_solver.core import (
    HomogeneousMaterial,
    TimeHarmonicProblem,
    TimeHarmonicSolver,
)
from fem_em_solver.io.mesh import MeshGenerator
from fem_em_solver.post.evaluation import evaluate_vector_field_parallel
from fem_em_solver.utils.analytical import (
    pec_cylinder_te1_field,
    pec_cylinder_tm0_current,
    pec_cylinder_tm0_field,
)
from tests.complex_mode import complex_only
from tests.validation.test_imp1_pec_sphere_lossy_host import (
    CAVITY_TAG,
    DEGREE,
    E0,
    EPSILON_R_HOST,
    FREQUENCY_HZ,
    MAT4_BAND,
    OUTER_BOUNDARY_TAG,
    SIGMA_HIGH,
    TH15_BETA_BAND,
    _assemble_real,
    _cavity_dofs,
)
from tests.validation.test_imp1_step3_exact_comparand import K_FIXTURE
from tests.validation.test_imp2_spheroid_fem import (
    QDEG,
    _numpy_quadrature_mean,
    _owned_masked_cells,
    _rss_gib,
)

A = 1.0e-3
C = 10.0e-3
LENGTH = 4.0e-3
END_TAG = 3
RESOLUTION_WALL = A / 6.0
RESOLUTION_FAR = C / 6.0
MASK_FACTOR = 1.5
PROBE_FACTORS = (1.2, 1.5, 2.0)
N_PROBE = 24
CELL_GUARD = 60_000
TAG_BAND = 0.02
M_BAND = 1e-12
S_BAND = 1e-12

OMEGA = 2.0 * np.pi * FREQUENCY_HZ
MU_0 = 4e-7 * np.pi
K = complex(K_FIXTURE)
EPS_C = K**2 / (OMEGA**2 * MU_0)

PEC_TAGS = {"tm": (OUTER_BOUNDARY_TAG, CAVITY_TAG, END_TAG), "te": (OUTER_BOUNDARY_TAG, CAVITY_TAG)}

_MESH_CACHE: dict = {}


def _say(msg):
    if MPI.COMM_WORLD.rank == 0:
        print(msg, flush=True)


def _mesh():
    if "m" not in _MESH_CACHE:
        t0 = time.perf_counter()
        out = MeshGenerator.cylindrical_domain(
            inner_radius=A,
            outer_radius=C,
            length=LENGTH,
            comm=MPI.COMM_WORLD,
            as_hole=True,
            resolution_wall=RESOLUTION_WALL,
            resolution_far=RESOLUTION_FAR,
        )
        _MESH_CACHE["m"] = (out, time.perf_counter() - t0)
    return _MESH_CACHE["m"][0]


def _ncells(msh):
    return int(msh.comm.allreduce(msh.topology.index_map(msh.topology.dim).size_local, op=MPI.SUM))


# ------------------------------------------------------------------ fields
def _cf(drive, interior="zero"):
    if drive == "tm":
        return lambda p: pec_cylinder_tm0_field(p, k=K, radius=A, e0=E0, interior=interior)
    return lambda p: pec_cylinder_te1_field(p, k=K, radius=A, omega=OMEGA, e0=E0, interior=interior)


def _dirichlet(drive):
    f = _cf(drive, "zero")

    def trace(x):
        if x.shape[1] == 0:
            return np.zeros((3, 0), dtype=PETSc.ScalarType)
        return np.ascontiguousarray(f(x.T).T).astype(PETSc.ScalarType)

    return trace


def _te_mode(points, kind):
    """TE n = 1 E with radial function J1 (kind 'J') or H1^(2) (kind 'H'),
    H0 = −2jωε_c E0/k: E_ρ = −H0 f sinφ/(jωε_c ρ), E_φ = −H0 f′ cosφ/(jωε_c)."""
    pts = np.asarray(points, dtype=float)
    out = np.zeros((pts.shape[0], 3), dtype=complex)
    if pts.shape[0] == 0:
        return out
    rho = np.hypot(pts[:, 0], pts[:, 1])
    if np.any(rho == 0.0):
        raise ValueError("TE mode sampled on the axis")
    cp, sp = pts[:, 0] / rho, pts[:, 1] / rho
    z = K * rho
    if kind == "J":
        f, fp = special.jv(1, z), K * special.jvp(1, z, 1)
    else:
        f, fp = special.hankel2(1, z), K * special.h2vp(1, z, 1)
    h0 = -2j * OMEGA * EPS_C * E0 / K
    pref = -h0 / (1j * OMEGA * EPS_C)
    e_rho = pref * f * sp / rho
    e_phi = pref * fp * cp
    out[:, 0] = e_rho * cp - e_phi * sp
    out[:, 1] = e_rho * sp + e_phi * cp
    return out


def _incident(drive):
    if drive == "tm":

        def inc(p):
            p = np.asarray(p, dtype=float)
            out = np.zeros((p.shape[0], 3), dtype=complex)
            out[:, 2] = E0 * special.jv(0, K * np.hypot(p[:, 0], p[:, 1]))
            return out

        return inc
    return lambda p: _te_mode(p, "J")


def _static_te_numpy(points):
    p = np.asarray(points, dtype=float)
    x, y = p[:, 0], p[:, 1]
    r4 = (x * x + y * y) ** 2
    out = np.zeros_like(p)
    out[:, 0] = 2.0 * E0 * A * A * x * y / r4
    out[:, 1] = E0 + E0 * A * A * (y * y - x * x) / r4
    return out


# ------------------------------------------------------------------ probes / masks
def _circle(rho, z=0.0, n=N_PROBE):
    ang = (np.arange(n) + 0.5) * 2.0 * np.pi / n
    return np.column_stack([rho * np.cos(ang), rho * np.sin(ang), np.full(n, z)])


def _probe_points():
    return np.vstack([_circle(s * A) for s in PROBE_FACTORS])


def _annulus_mask(msh):
    dg0 = fem.functionspace(msh, ("DG", 0))
    chi = fem.Function(dg0)
    r2 = (MASK_FACTOR * A) ** 2
    chi.interpolate(lambda x: ((x[0] ** 2 + x[1] ** 2) < r2).astype(np.float64))
    return chi


def _continuous_annulus_mean(field_fn, n_r=64, n_phi=64):
    """⟨|E|²⟩/E₀² over a < ρ < 1.5a (z = 0), Gauss–Legendre in ρ × midpoint in φ."""
    xg, wg = np.polynomial.legendre.leggauss(n_r)
    r0, r1 = A, MASK_FACTOR * A
    rho = 0.5 * (r1 - r0) * xg + 0.5 * (r1 + r0)
    wr = 0.5 * (r1 - r0) * wg * rho
    ang = (np.arange(n_phi) + 0.5) * 2.0 * np.pi / n_phi
    rr, aa = np.meshgrid(rho, ang, indexing="ij")
    pts = np.column_stack([(rr * np.cos(aa)).ravel(), (rr * np.sin(aa)).ravel(), np.zeros(rr.size)])
    e2 = np.sum(np.abs(field_fn(pts)) ** 2, axis=1).reshape(rr.shape)
    return float(np.sum(wr[:, None] * e2) / (np.sum(wr) * n_phi) / E0**2)


def _eval_complex(fields, points, comm):
    re_vals, ok_re = evaluate_vector_field_parallel(fields.e_real, points, comm)
    im_vals, ok_im = evaluate_vector_field_parallel(fields.e_imag, points, comm)
    if not np.all(ok_re & ok_im):
        raise RuntimeError("probe points were not evaluated")
    return np.real(re_vals[:, :3]) + 1j * np.real(im_vals[:, :3])


def _h_wall_over_a(msh):
    dmap = msh.geometry.dofmap
    n_owned = msh.topology.index_map(msh.topology.dim).size_local
    v = msh.geometry.x[dmap[:n_owned]]
    if n_owned:
        mid = v.mean(axis=1)
        near = np.hypot(mid[:, 0], mid[:, 1]) < 1.1 * A
        edges = np.stack(
            [np.linalg.norm(v[:, i] - v[:, j], axis=1) for i in range(4) for j in range(i + 1, 4)],
            axis=1,
        )
        local = float(np.max(edges[near])) if np.any(near) else 0.0
    else:
        local = 0.0
    return msh.comm.allreduce(local, op=MPI.MAX) / A


# ------------------------------------------------------------------ solve
def _solve_cylinder(drive, pec_tags, *, return_fields=False, trace=None):
    comm = MPI.COMM_WORLD
    t0 = time.perf_counter()
    msh, cell_tags, facet_tags = _mesh()
    problem = TimeHarmonicProblem(
        mesh=msh,
        frequency_hz=FREQUENCY_HZ,
        material=HomogeneousMaterial(sigma=SIGMA_HIGH, epsilon_r=EPSILON_R_HOST),
        cell_tags=cell_tags,
        facet_tags=facet_tags,
        boundary_condition="pec_zero_tangential_a",
        dirichlet_e_field=_dirichlet(drive) if trace is None else trace,
        pec_facet_tags=tuple(pec_tags),
    )
    solver = TimeHarmonicSolver(problem, degree=DEGREE)
    fields = solver.solve()
    t_solve = time.perf_counter() - t0

    pts = _probe_points()
    e_h = _eval_complex(fields, pts, comm)
    e_inc = _incident(drive)(pts)
    e_s = _cf(drive)(pts) - e_inc
    beta = complex(np.sum(np.conj(e_s) * (e_h - e_inc)) / np.sum(np.abs(e_s) ** 2))

    zc = [_eval_complex(fields, _circle(1.5 * A, z), comm) for z in (0.0, LENGTH / 4, -LENGTH / 4)]
    z_inv = float(
        max(
            np.max(np.linalg.norm(zc[1] - zc[0], axis=1)),
            np.max(np.linalg.norm(zc[2] - zc[0], axis=1)),
        )
        / np.max(np.linalg.norm(zc[0], axis=1))
    )

    chi = _annulus_mask(msh)
    e = fields.e_complex
    dx = ufl.Measure("dx", domain=msh, metadata={"quadrature_degree": QDEG})
    num = _assemble_real(chi * ufl.inner(e, e) * dx, comm)
    vol = _assemble_real(chi * dx, comm)
    phi_h = num / (E0**2 * vol)
    cells = _owned_masked_cells(msh, chi)
    n_mask = int(comm.allreduce(int(cells.size), op=MPI.SUM))
    phi_cf_h = _numpy_quadrature_mean(msh, cells, _cf(drive, "continue"))
    phi_cont = _continuous_annulus_mean(_cf(drive, "continue"))

    v_space = solver.function_space()
    cavity = _cavity_dofs(msh, facet_tags, v_space)
    n_cavity = int(comm.allreduce(int(cavity.size), op=MPI.SUM))
    sol = np.asarray(e.x.array)
    local_max = float(np.max(np.abs(sol[cavity]))) if cavity.size else 0.0
    cavity_max = float(comm.allreduce(local_max, op=MPI.MAX))
    ndofs = int(v_space.dofmap.index_map.size_global * v_space.dofmap.index_map_bs)
    h_wall = _h_wall_over_a(msh)
    rss = _rss_gib(comm)  # collective
    out = dict(
        beta=beta,
        phi_h=phi_h,
        phi_cf_h=phi_cf_h,
        phi_cont=phi_cont,
        n_mask=n_mask,
        n_cavity=n_cavity,
        cavity_max=cavity_max,
        ncells=_ncells(msh),
        ndofs=ndofs,
        h_wall=h_wall,
        z_inv=z_inv,
        seconds=time.perf_counter() - t0,
        t_solve=t_solve,
        rss=rss,
    )
    if return_fields:
        out["fields"] = fields
    return out


def _report(label, drive, pec_tags, d):
    miss = d["phi_h"] / d["phi_cf_h"] - 1.0
    _say(
        f"\n[IMP-4 step 1] {label}: drive {drive}, pec tags {pec_tags}; cells "
        f"{d['ncells']}, degree-{DEGREE} dofs {d['ndofs']}, {d['seconds']:.1f} s "
        f"(solve {d['t_solve']:.1f} s), RSS(sum) {d['rss']:.2f} GiB, h_wall/a "
        f"{d['h_wall']:.3f}, annulus cells {d['n_mask']}"
    )
    _say(
        f"  beta = {d['beta'].real:.6f}{d['beta'].imag:+.6f}j  |beta - 1| = "
        f"{abs(d['beta'] - 1):.4%} (band {TH15_BETA_BAND:.4%})"
    )
    _say(
        f"  Phi_h = {d['phi_h']:.6e}, Phi_cf,h = {d['phi_cf_h']:.6e}, continuous "
        f"Phi = {d['phi_cont']:.6e} (Phi_cf,h/cont - 1 = "
        f"{d['phi_cf_h'] / d['phi_cont'] - 1:+.4%}); Phi_h/Phi_cf,h - 1 = "
        f"{miss:+.4%} (band {MAT4_BAND:.0%})"
    )
    _say(
        f"  cavity dofs {d['n_cavity']}, max |E| {d['cavity_max']:.3e}; "
        f"z-invariance {d['z_inv']:.3e}"
    )
    return miss


# --------------------------------------------------------------------- W0
def test_census_and_numpy():
    comm = MPI.COMM_WORLD
    _say(
        f"\n[IMP-4 step 1] k = {K:.6g} 1/m, |k|a = {abs(K) * A:.4f}, |k|c = "
        f"{abs(K) * C:.4f}; a {A * 1e3} mm, c {C * 1e3} mm, L {LENGTH * 1e3} mm"
    )
    msh, _, facet_tags = _mesh()
    ncells = _ncells(msh)
    one = fem.Constant(msh, PETSc.ScalarType(1.0))
    ds = ufl.Measure("ds", domain=msh, subdomain_data=facet_tags, metadata={"quadrature_degree": 2})
    refs = {
        OUTER_BOUNDARY_TAG: 2 * np.pi * C * LENGTH,
        CAVITY_TAG: 2 * np.pi * A * LENGTH,
        END_TAG: 2 * np.pi * (C**2 - A**2),
    }
    areas = {t: _assemble_real(one * ds(t), comm) for t in refs}
    nfacets = {t: int(comm.allreduce(int(facet_tags.find(t).size), op=MPI.SUM)) for t in refs}

    # (M') numpy quadrature vs UFL cell integral of the static TE field.
    chi = _annulus_mask(msh)
    cells = _owned_masked_cells(msh, chi)
    n_mask = int(comm.allreduce(int(cells.size), op=MPI.SUM))
    mq = _numpy_quadrature_mean(msh, cells, _static_te_numpy)
    x = ufl.SpatialCoordinate(msh)
    r4 = (x[0] ** 2 + x[1] ** 2) ** 2
    e_st = ufl.as_vector(
        (2.0 * E0 * A * A * x[0] * x[1] / r4, E0 + E0 * A * A * (x[1] ** 2 - x[0] ** 2) / r4, 0.0)
    )
    dx = ufl.Measure("dx", domain=msh, metadata={"quadrature_degree": QDEG})
    mu_ufl = _assemble_real(chi * ufl.dot(e_st, e_st) * dx, comm) / (
        E0**2 * _assemble_real(chi * dx, comm)
    )
    m_err = abs(mq / mu_ufl - 1.0)
    h_wall = _h_wall_over_a(msh)
    rss = _rss_gib(comm)

    _say(
        f"[IMP-4 census] cells {ncells} (guard {CELL_GUARD}), mesh "
        f"{_MESH_CACHE['m'][1]:.1f} s, annulus-mask cells {n_mask}, h_wall/a "
        f"{h_wall:.3f}, RSS(sum) {rss:.2f} GiB"
    )
    for t, ref in refs.items():
        _say(
            f"[IMP-4 Tag] tag {t}: facets {nfacets[t]}, area {areas[t]:.9e} vs "
            f"{ref:.9e}, rel {areas[t] / ref - 1:+.4%} (band {TAG_BAND:.0%})"
        )
    _say(
        f"[IMP-4 M'] numpy {mq:.15f} vs UFL {mu_ufl:.15f} rel {m_err:.3e}; "
        f"continuous 1 + <a^4/rho^4> = {1 + 4.0 / 9.0:.4f}"
    )

    # (S) the incident/scattered split, 64 points in [1.2a, 2a].
    rng = np.random.default_rng(206)
    rho = np.linspace(1.2 * A, 2.0 * A, 64)
    ang = rng.uniform(0.0, 2.0 * np.pi, 64)
    pts = np.column_stack([rho * np.cos(ang), rho * np.sin(ang), np.zeros(64)])
    tm = _cf("tm")(pts)[:, 2]
    ratio = (tm - E0 * special.jv(0, K * rho)) / special.hankel2(0, K * rho)
    target = -E0 * special.jv(0, K * A) / special.hankel2(0, K * A)
    s_tm = float(np.max(np.abs(ratio - target)) / abs(target))
    coef = special.jvp(1, K * A, 1) / special.h2vp(1, K * A, 1)
    scat = _cf("te")(pts) - _te_mode(pts, "J")
    pred = -coef * _te_mode(pts, "H")
    s_te = float(np.max(np.linalg.norm(scat - pred, axis=1)) / np.max(np.linalg.norm(pred, axis=1)))
    _say(
        f"[IMP-4 S] TM const ratio vs -J0(ka)/H0(ka) = {target:.6g}: rel {s_tm:.3e}; "
        f"TE scattered vs -J1'(ka)/H1'(ka) H-mode (coef {coef:.6g}): rel {s_te:.3e}"
    )

    assert ncells <= CELL_GUARD, f"cell guard: {ncells} > {CELL_GUARD}"
    for t, ref in refs.items():
        assert abs(areas[t] / ref - 1.0) <= TAG_BAND, (t, areas[t], ref)
    assert n_mask > 0
    assert m_err <= M_BAND
    assert s_tm <= S_BAND
    assert s_te <= S_BAND


# --------------------------------------------------------------------- W1
@complex_only
@pytest.mark.integration
def test_both_drives():
    results = {}
    for drive in ("tm", "te"):
        d = _solve_cylinder(drive, PEC_TAGS[drive])
        _report("gate", drive, PEC_TAGS[drive], d)
        if drive == "tm":
            i_cf = pec_cylinder_tm0_current(k=K, radius=A, omega=OMEGA, e0=E0)
            i_fit = d["beta"] * i_cf
            _say(
                f"  TM current beta*I = {i_fit:.6g} A per V/m (|I| {abs(i_fit):.6g}, "
                f"closed form {i_cf:.6g}, |I_cf| {abs(i_cf):.6g})"
            )
        results[drive] = d
    for drive, d in results.items():
        assert d["n_cavity"] > 0, drive
        assert d["cavity_max"] < 1e-12, (drive, d["cavity_max"])
        assert abs(d["beta"] - 1.0) <= TH15_BETA_BAND, (drive, d["beta"], d["h_wall"])
        assert abs(d["phi_h"] / d["phi_cf_h"] - 1.0) <= MAT4_BAND, (
            drive,
            d["phi_h"],
            d["phi_cf_h"],
            d["phi_cont"],
            d["z_inv"],
        )


# --------------------------------------------------------------------- W2
@complex_only
@pytest.mark.integration
def test_natural_cavity_controls():
    # (c1) TM, natural cavity (PMC wire): predicted |beta-1| ~0.49 [0.3, 0.7],
    # Phi miss ~ +7900 % [+3000, +15000] %.
    d1 = _solve_cylinder("tm", (OUTER_BOUNDARY_TAG, END_TAG))
    m1 = _report("control c1", "tm", (OUTER_BOUNDARY_TAG, END_TAG), d1)
    _say(
        f"  c1 predicted: beta ~ 0.53-0.13j, |beta-1| ~ 0.49 [0.3, 0.7]; Phi miss "
        f"~ +7900 % [+3000, +15000] %; measured |beta-1| {abs(d1['beta'] - 1):.4f}, "
        f"miss {m1:+.2%}"
    )
    # (c2) TE, natural cavity (void cylinder): predicted beta ~ -0.980,
    # |beta-1| ~ 1.98; Phi miss ~ 1 - A^2 ~ 3.9 % (printed only).
    d2 = _solve_cylinder("te", (OUTER_BOUNDARY_TAG,))
    m2 = _report("control c2", "te", (OUTER_BOUNDARY_TAG,), d2)
    _say(
        f"  c2 predicted: beta ~ -0.980, |beta-1| ~ 1.98; Phi miss ~ 3.9 % (not "
        f"asserted); measured |beta-1| {abs(d2['beta'] - 1):.4f}, miss {m2:+.2%}"
    )
    assert abs(d1["beta"] - 1.0) > TH15_BETA_BAND
    assert abs(m1) > MAT4_BAND
    assert abs(d2["beta"] - 1.0) > TH15_BETA_BAND
