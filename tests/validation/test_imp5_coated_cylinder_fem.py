"""`IMP-5` step 1 (§9 item 211, lane 3): FEM for the dielectric-coated
infinite PEC circular cylinder in `IMP-1`'s lossy host at finite k, against
step 0's closed form (``coated_pec_cylinder_*``).

Fixture: `IMP-4` step 1's hole route with a conforming εᵣ 3 lossless sleeve
a < ρ < b (b = 1.5a) meshed as its own cell tag 4 by
``MeshGenerator.cylindrical_domain(as_hole=True, coat_radius=b)`` and assigned
through ``material_map``.  The outer wall carries step 0's exact coated trace;
the ends carry `IMP-4`'s polarisation-specific conditions (TM PEC: pec tags
1, 2, 3; TE natural: pec tags 1, 2).  The ρ = b interface needs no condition:
the curl-curl weak form carries tangential E and tangential μ⁻¹ curl E
continuity naturally, and Nédélec elements leave the normal component free to
jump (≈ 59× in TE).

Gated (degree 2, each drive): (i) |β − 1| ≤ `TH15_BETA_BAND` (β the complex
LSQ amplitude of the closed-form scattered field on 24-point circles at
ρ = 1.2b, 1.5b, 2b); (ii) sleeve (cell tag 4) and (iii) near-lead host
annulus (cell tag 2, midpoint ρ < 1.5b) mean |E|² against the closed form
integrated by numpy quadrature over the same owned cells with the region
**forced by cell tag** within `MAT4_BAND`; (iv) cavity dofs zero.  W0: cell
guard, (R) the ``coat_radius=None`` bit-identity on `IMP-4` W0's 25 568
cells, tag areas, sleeve/host volumes, (M′) on both masks, (S) the
incident/scattered split.  Controls (W2): TE coating omitted; TM natural
cavity.
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
from fem_em_solver.utils.analytical import (
    _coated_cylinder_coefficients,
    coated_pec_cylinder_te1_field,
    coated_pec_cylinder_tm0_current,
    coated_pec_cylinder_tm0_field,
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
from tests.validation.test_imp1_step3_exact_comparand import K0, K_FIXTURE
from tests.validation.test_imp2_spheroid_fem import (
    QDEG,
    _numpy_quadrature_mean,
    _owned_masked_cells,
    _rss_gib,
)
from tests.validation.test_imp4_pec_cylinder_fem import (
    CELL_GUARD,
    END_TAG,
    LENGTH,
    M_BAND,
    PEC_TAGS,
    RESOLUTION_FAR,
    RESOLUTION_WALL,
    S_BAND,
    TAG_BAND,
    A,
    C,
    _circle,
    _eval_complex,
    _h_wall_over_a,
    _incident,
    _ncells,
    _say,
    _static_te_numpy,
    _te_mode,
)

B = 1.5 * A
COAT_TAG = 4
HOST_TAG = 2
EPS_R_COAT = 3.0
K_COAT = complex(K0 * np.sqrt(3.0))
K = complex(K_FIXTURE)
OMEGA = 2.0 * np.pi * FREQUENCY_HZ
HOST_MASK_FACTOR = 1.5
PROBE_FACTORS = (1.2, 1.5, 2.0)
IMP4_W0_CELLS = 25_568  # 20261003T200309Z_IMP-4-step1-W0.log:465

_MESH_CACHE: dict = {}


def _mesh(coated=True):
    key = "coat" if coated else "bare"
    if key not in _MESH_CACHE:
        t0 = time.perf_counter()
        out = MeshGenerator.cylindrical_domain(
            inner_radius=A,
            outer_radius=C,
            length=LENGTH,
            comm=MPI.COMM_WORLD,
            as_hole=True,
            resolution_wall=RESOLUTION_WALL,
            resolution_far=RESOLUTION_FAR,
            coat_radius=B if coated else None,
        )
        _MESH_CACHE[key] = (out, time.perf_counter() - t0)
    return _MESH_CACHE[key][0]


# ------------------------------------------------------------------ fields
def _cf(drive, region=None, interior="zero", condition="exact"):
    kw = dict(
        k_host=K, k_coat=K_COAT, radius=A, coat_radius=B, e0=E0, interior=interior, region=region
    )
    if drive == "tm":
        return lambda p: coated_pec_cylinder_tm0_field(p, **kw)
    return lambda p: coated_pec_cylinder_te1_field(p, omega=OMEGA, **kw)


def _dirichlet(drive):
    f = _cf(drive, region=None, interior="zero")

    def trace(x):
        if x.shape[1] == 0:
            return np.zeros((3, 0), dtype=PETSc.ScalarType)
        return np.ascontiguousarray(f(x.T).T).astype(PETSc.ScalarType)

    return trace


# ------------------------------------------------------------------ masks
def _all_cell_midpoints(msh):
    imap = msh.topology.index_map(msh.topology.dim)
    n = imap.size_local + imap.num_ghosts
    v = msh.geometry.x[msh.geometry.dofmap[:n]]
    return v.mean(axis=1)


def _tag_mask(msh, cell_tags, tag, rho_max=None):
    """DG0 indicator from cell tags (numpy; never a radius test at an
    interface).  ``rho_max`` additionally restricts by cell midpoint."""
    dg0 = fem.functionspace(msh, ("DG", 0))
    chi = fem.Function(dg0)
    cells = cell_tags.find(tag)
    if rho_max is not None and cells.size:
        mid = _all_cell_midpoints(msh)[cells]
        cells = cells[(mid[:, 0] ** 2 + mid[:, 1] ** 2) < rho_max**2]
    dofs = dg0.dofmap.list[cells, 0] if cells.size else np.zeros(0, dtype=np.int32)
    chi.x.array[:] = 0.0
    chi.x.array[dofs] = 1.0
    chi.x.scatter_forward()
    return chi


def _masks(msh, cell_tags):
    return {
        "coat": _tag_mask(msh, cell_tags, COAT_TAG),
        "host": _tag_mask(msh, cell_tags, HOST_TAG, HOST_MASK_FACTOR * B),
    }


def _continuous_mean(field_fn, r0, r1, n_r=64, n_phi=64):
    """⟨|E|²⟩/E₀² over r0 < ρ < r1 (z = 0), Gauss–Legendre in ρ × midpoint in φ."""
    xg, wg = np.polynomial.legendre.leggauss(n_r)
    rho = 0.5 * (r1 - r0) * xg + 0.5 * (r1 + r0)
    wr = 0.5 * (r1 - r0) * wg * rho
    ang = (np.arange(n_phi) + 0.5) * 2.0 * np.pi / n_phi
    rr, aa = np.meshgrid(rho, ang, indexing="ij")
    pts = np.column_stack([(rr * np.cos(aa)).ravel(), (rr * np.sin(aa)).ravel(), np.zeros(rr.size)])
    e2 = np.sum(np.abs(field_fn(pts)) ** 2, axis=1).reshape(rr.shape)
    return float(np.sum(wr[:, None] * e2) / (np.sum(wr) * n_phi) / E0**2)


MASK_RANGE = {"coat": (A, B), "host": (B, HOST_MASK_FACTOR * B)}


def _sleeve_edge_over_t(msh, cell_tags):
    cells = cell_tags.find(COAT_TAG)
    n_owned = msh.topology.index_map(msh.topology.dim).size_local
    cells = cells[cells < n_owned]
    if cells.size:
        v = msh.geometry.x[msh.geometry.dofmap[cells]]
        local = float(
            max(
                np.max(np.linalg.norm(v[:, i] - v[:, j], axis=1))
                for i in range(4)
                for j in range(i + 1, 4)
            )
        )
    else:
        local = 0.0
    return msh.comm.allreduce(local, op=MPI.MAX) / (B - A)


def _probe_points():
    return np.vstack([_circle(s * B) for s in PROBE_FACTORS])


# ------------------------------------------------------------------ solve
def _solve_coated(drive, pec_tags, coated=True, *, return_fields=False, trace=None):
    comm = MPI.COMM_WORLD
    t0 = time.perf_counter()
    msh, cell_tags, facet_tags = _mesh(True)
    mmap = {COAT_TAG: HomogeneousMaterial(sigma=0.0, epsilon_r=EPS_R_COAT)} if coated else None
    problem = TimeHarmonicProblem(
        mesh=msh,
        frequency_hz=FREQUENCY_HZ,
        material=HomogeneousMaterial(sigma=SIGMA_HIGH, epsilon_r=EPSILON_R_HOST),
        material_map=mmap,
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
    e_s = _cf(drive, region="host")(pts) - e_inc
    beta = complex(np.sum(np.conj(e_s) * (e_h - e_inc)) / np.sum(np.abs(e_s) ** 2))

    zc = [_eval_complex(fields, _circle(1.5 * B, z), comm) for z in (0.0, LENGTH / 4, -LENGTH / 4)]
    z_inv = float(
        max(
            np.max(np.linalg.norm(zc[1] - zc[0], axis=1)),
            np.max(np.linalg.norm(zc[2] - zc[0], axis=1)),
        )
        / np.max(np.linalg.norm(zc[0], axis=1))
    )

    e = fields.e_complex
    dx = ufl.Measure("dx", domain=msh, metadata={"quadrature_degree": QDEG})
    phi = {}
    for name, chi in _masks(msh, cell_tags).items():
        num = _assemble_real(chi * ufl.inner(e, e) * dx, comm)
        vol = _assemble_real(chi * dx, comm)
        cells = _owned_masked_cells(msh, chi)
        n_mask = int(comm.allreduce(int(cells.size), op=MPI.SUM))
        phi_cf = _numpy_quadrature_mean(msh, cells, _cf(drive, region=name))
        phi_cont = _continuous_mean(_cf(drive, region=name), *MASK_RANGE[name])
        phi[name] = dict(h=num / (E0**2 * vol), cf=phi_cf, cont=phi_cont, n=n_mask)

    v_space = solver.function_space()
    cavity = _cavity_dofs(msh, facet_tags, v_space)
    n_cavity = int(comm.allreduce(int(cavity.size), op=MPI.SUM))
    sol = np.asarray(e.x.array)
    local_max = float(np.max(np.abs(sol[cavity]))) if cavity.size else 0.0
    cavity_max = float(comm.allreduce(local_max, op=MPI.MAX))
    ndofs = int(v_space.dofmap.index_map.size_global * v_space.dofmap.index_map_bs)
    h_wall = _h_wall_over_a(msh)
    sleeve_edge = _sleeve_edge_over_t(msh, cell_tags)
    ntag = {
        t: int(
            comm.allreduce(
                int(np.sum(cell_tags.find(t) < msh.topology.index_map(3).size_local)), op=MPI.SUM
            )
        )
        for t in (HOST_TAG, COAT_TAG)
    }
    rss = _rss_gib(comm)  # collective, outside any rank-0 branch
    out = dict(
        beta=beta,
        phi=phi,
        n_cavity=n_cavity,
        cavity_max=cavity_max,
        ncells=_ncells(msh),
        ntag=ntag,
        ndofs=ndofs,
        h_wall=h_wall,
        sleeve_edge=sleeve_edge,
        z_inv=z_inv,
        seconds=time.perf_counter() - t0,
        t_solve=t_solve,
        rss=rss,
    )
    if return_fields:
        out["fields"] = fields
    return out


def _report(label, drive, pec_tags, d, coated=True):
    _say(
        f"\n[IMP-5 step 1] {label}: drive {drive}, pec tags {pec_tags}, coated "
        f"{coated}; cells {d['ncells']} (per tag {d['ntag']}), degree-{DEGREE} dofs "
        f"{d['ndofs']}, {d['seconds']:.1f} s (solve {d['t_solve']:.1f} s), RSS(sum) "
        f"{d['rss']:.2f} GiB, h_wall/a {d['h_wall']:.3f}, sleeve max edge/(b-a) "
        f"{d['sleeve_edge']:.3f}"
    )
    _say(
        f"  beta = {d['beta'].real:.6f}{d['beta'].imag:+.6f}j  |beta - 1| = "
        f"{abs(d['beta'] - 1):.4%} (band {TH15_BETA_BAND:.4%})"
    )
    misses = {}
    for name, p in d["phi"].items():
        misses[name] = p["h"] / p["cf"] - 1.0
        _say(
            f"  [{name}] cells {p['n']}: Phi_h = {p['h']:.6e}, Phi_cf,h = "
            f"{p['cf']:.6e}, continuous Phi = {p['cont']:.6e} (Phi_cf,h/cont - 1 = "
            f"{p['cf'] / p['cont'] - 1:+.4%}); Phi_h/Phi_cf,h - 1 = "
            f"{misses[name]:+.4%} (band {MAT4_BAND:.0%})"
        )
    _say(
        f"  cavity dofs {d['n_cavity']}, max |E| {d['cavity_max']:.3e}; "
        f"z-invariance {d['z_inv']:.3e}"
    )
    return misses


def _assert_gates(drive, d, misses):
    assert d["n_cavity"] > 0, drive
    assert d["cavity_max"] < 1e-12, (drive, d["cavity_max"])
    assert abs(d["beta"] - 1.0) <= TH15_BETA_BAND, (drive, d["beta"], d["h_wall"])
    for name in ("coat", "host"):
        p = d["phi"][name]
        assert abs(misses[name]) <= MAT4_BAND, (drive, name, p, d["z_inv"])


# --------------------------------------------------------------------- W0
def test_census_and_numpy():
    comm = MPI.COMM_WORLD
    _say(
        f"\n[IMP-5 step 1] k_host = {K:.6g} 1/m, k_coat = {K_COAT:.6g} 1/m, "
        f"r = k1^2/k2^2 = {K_COAT**2 / K**2:.4g}; a {A * 1e3} mm, b {B * 1e3} mm, "
        f"c {C * 1e3} mm, L {LENGTH * 1e3} mm"
    )
    # (R) the coat_radius=None branch reproduces IMP-4 W0's cell count.
    bare, _, _ = _mesh(False)
    n_bare = _ncells(bare)
    msh, cell_tags, facet_tags = _mesh(True)
    ncells = _ncells(msh)
    n_owned = msh.topology.index_map(3).size_local
    ntag = {
        t: int(comm.allreduce(int(np.sum(cell_tags.find(t) < n_owned)), op=MPI.SUM))
        for t in (HOST_TAG, COAT_TAG)
    }
    one = fem.Constant(msh, PETSc.ScalarType(1.0))
    ds = ufl.Measure("ds", domain=msh, subdomain_data=facet_tags, metadata={"quadrature_degree": 2})
    refs = {
        OUTER_BOUNDARY_TAG: 2 * np.pi * C * LENGTH,
        CAVITY_TAG: 2 * np.pi * A * LENGTH,
        END_TAG: 2 * np.pi * (C**2 - A**2),
    }
    areas = {t: _assemble_real(one * ds(t), comm) for t in refs}
    nfacets = {t: int(comm.allreduce(int(facet_tags.find(t).size), op=MPI.SUM)) for t in refs}
    dxt = ufl.Measure("dx", domain=msh, subdomain_data=cell_tags, metadata={"quadrature_degree": 2})
    vrefs = {COAT_TAG: np.pi * (B**2 - A**2) * LENGTH, HOST_TAG: np.pi * (C**2 - B**2) * LENGTH}
    vols = {t: _assemble_real(one * dxt(t), comm) for t in vrefs}

    # (M') numpy quadrature vs UFL cell integral of IMP-4's static TE field.
    x = ufl.SpatialCoordinate(msh)
    r4 = (x[0] ** 2 + x[1] ** 2) ** 2
    e_st = ufl.as_vector(
        (2.0 * E0 * A * A * x[0] * x[1] / r4, E0 + E0 * A * A * (x[1] ** 2 - x[0] ** 2) / r4, 0.0)
    )
    dx = ufl.Measure("dx", domain=msh, metadata={"quadrature_degree": QDEG})
    m_err, n_mask = {}, {}
    for name, chi in _masks(msh, cell_tags).items():
        cells = _owned_masked_cells(msh, chi)
        n_mask[name] = int(comm.allreduce(int(cells.size), op=MPI.SUM))
        mq = _numpy_quadrature_mean(msh, cells, _static_te_numpy)
        mu = _assemble_real(chi * ufl.dot(e_st, e_st) * dx, comm) / (
            E0**2 * _assemble_real(chi * dx, comm)
        )
        m_err[name] = abs(mq / mu - 1.0)
        _say(
            f"[IMP-5 M'] {name} mask cells {n_mask[name]}: numpy {mq:.15f} vs UFL "
            f"{mu:.15f} rel {m_err[name]:.3e}"
        )
    h_wall = _h_wall_over_a(msh)
    sleeve_edge = _sleeve_edge_over_t(msh, cell_tags)
    rss = _rss_gib(comm)

    _say(f"[IMP-5 R] coat_radius=None cells {n_bare} vs IMP-4 W0 {IMP4_W0_CELLS}")
    _say(
        f"[IMP-5 census] cells {ncells} (guard {CELL_GUARD}), per tag {ntag}, mesh "
        f"{_MESH_CACHE['coat'][1]:.1f} s, h_wall/a {h_wall:.3f}, sleeve max "
        f"edge/(b-a) {sleeve_edge:.3f}, RSS(sum) {rss:.2f} GiB"
    )
    for t, ref in refs.items():
        _say(
            f"[IMP-5 Tag] tag {t}: facets {nfacets[t]}, area {areas[t]:.9e} vs "
            f"{ref:.9e}, rel {areas[t] / ref - 1:+.4%} (band {TAG_BAND:.0%})"
        )
    for t, ref in vrefs.items():
        _say(
            f"[IMP-5 Vol] tag {t}: volume {vols[t]:.9e} vs {ref:.9e}, rel "
            f"{vols[t] / ref - 1:+.4%} (band {TAG_BAND:.0%})"
        )

    # (S) the incident/scattered split, 64 points in [1.2b, 2b].
    rng = np.random.default_rng(211)
    rho = np.linspace(1.2 * B, 2.0 * B, 64)
    ang = rng.uniform(0.0, 2.0 * np.pi, 64)
    pts = np.column_stack([rho * np.cos(ang), rho * np.sin(ang), np.zeros(64)])
    co = dict(k_host=K, k_coat=K_COAT, radius=A, coat_radius=B)
    _, s_coef = _coated_cylinder_coefficients("tm0", **co)
    _, t_coef = _coated_cylinder_coefficients("te1", **co)
    tm = _cf("tm", region="host")(pts)[:, 2]
    ratio = (tm - E0 * special.jv(0, K * rho)) / special.hankel2(0, K * rho)
    s_tm = float(np.max(np.abs(ratio - E0 * s_coef)) / abs(E0 * s_coef))
    scat = _cf("te", region="host")(pts) - _te_mode(pts, "J")
    pred = t_coef * _te_mode(pts, "H")
    s_te = float(np.max(np.linalg.norm(scat - pred, axis=1)) / np.max(np.linalg.norm(pred, axis=1)))
    t0_bare = -special.jvp(1, K * A, 1) / special.h2vp(1, K * A, 1)
    _say(
        f"[IMP-5 S] TM ratio vs s = {s_coef:.6g}: rel {s_tm:.3e}; TE scattered vs "
        f"t = {t_coef:.6g} H-mode: rel {s_te:.3e} (t/t0 = {t_coef / t0_bare:.5g})"
    )

    assert ncells <= CELL_GUARD, f"cell guard: {ncells} > {CELL_GUARD}"
    assert n_bare == IMP4_W0_CELLS, (n_bare, IMP4_W0_CELLS)
    for t, ref in refs.items():
        assert abs(areas[t] / ref - 1.0) <= TAG_BAND, (t, areas[t], ref)
    for t, ref in vrefs.items():
        assert abs(vols[t] / ref - 1.0) <= TAG_BAND, (t, vols[t], ref)
    for name in ("coat", "host"):
        assert n_mask[name] > 0, name
        assert m_err[name] <= M_BAND, (name, m_err[name])
    assert s_tm <= S_BAND
    assert s_te <= S_BAND


# --------------------------------------------------------------------- W1
@complex_only
@pytest.mark.integration
def test_tm_drive():
    d = _solve_coated("tm", PEC_TAGS["tm"])
    misses = _report("gate", "tm", PEC_TAGS["tm"], d)
    i_cf = coated_pec_cylinder_tm0_current(
        k_host=K, k_coat=K_COAT, radius=A, coat_radius=B, omega=OMEGA, e0=E0
    )
    i_fit = d["beta"] * i_cf
    _say(
        f"  TM current beta*I = {i_fit:.6g} A per V/m (closed form {i_cf:.6g}); "
        f"predicted |beta-1| <~ 0.2 %, coat miss +1..+3 %, host miss +0.3..+2 %, "
        f"Phi_cont coat 3.640e-3 host 2.296e-2"
    )
    _assert_gates("tm", d, misses)


@complex_only
@pytest.mark.integration
def test_te_drive():
    d = _solve_coated("te", PEC_TAGS["te"])
    misses = _report("gate", "te", PEC_TAGS["te"], d)
    _say(
        "  predicted |beta-1| <~ 2 %, |coat miss| <~ 2 %, |host miss| <~ 1.5 %, "
        "Phi_cont coat 17.99 host 1.411"
    )
    _assert_gates("te", d, misses)


# --------------------------------------------------------------------- W2
@complex_only
@pytest.mark.integration
def test_control_te_coating_omitted():
    d = _solve_coated("te", PEC_TAGS["te"], coated=False)
    m = _report("control c1", "te", PEC_TAGS["te"], d, coated=False)
    _say(
        f"  c1 predicted: beta ~ -0.476-0.036j, |beta-1| ~ 1.476 [1.3, 1.65]; coat "
        f"miss ~ -91.4 % [-95, -85] %; host miss ~ -17.8 % [-22, -14] %; measured "
        f"|beta-1| {abs(d['beta'] - 1):.4f}, coat {m['coat']:+.2%}, host "
        f"{m['host']:+.2%}"
    )
    assert abs(d["beta"] - 1.0) > TH15_BETA_BAND
    assert abs(m["coat"]) > MAT4_BAND
    assert abs(m["host"]) > MAT4_BAND


@complex_only
@pytest.mark.integration
def test_control_tm_natural_cavity():
    tags = (OUTER_BOUNDARY_TAG, END_TAG)
    d = _solve_coated("tm", tags)
    m = _report("control c2", "tm", tags, d)
    _say(
        f"  c2 predicted: beta ~ 0.590-0.123j, |beta-1| ~ 0.428 [0.25, 0.6]; coat "
        f"miss ~ +7940 % [+3000, +15000] %; host miss ~ +1175 % [+500, +2500] %; "
        f"measured |beta-1| {abs(d['beta'] - 1):.4f}, coat {m['coat']:+.2%}, host "
        f"{m['host']:+.2%}"
    )
    assert abs(d["beta"] - 1.0) > TH15_BETA_BAND
    assert abs(m["coat"]) > MAT4_BAND
    assert abs(m["host"]) > MAT4_BAND
