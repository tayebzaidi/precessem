"""`IMP-2` step 1 (§9 item 144, chain I4): FEM for a PEC prolate spheroid in
`IMP-1`'s lossy host under the quasi-static axial trace, aspect 2 and 5.

Fixture: `TH-15`'s hole route — the spheroid (c = 6 mm along z, a = c/aspect)
cut from a box of half-extents (c, c, 2c) by the new
``MeshGenerator.spheroid_in_box_domain(as_hole=True)``; the outer walls carry
the closed-form tangential trace ``pec_prolate_spheroid_field`` (step 0's
potential differentiated), the cavity wall zero.  ``(|k| c)² = 1.154 %`` is
under a third of the tighter band (4.886 % / 3), so the quasi-static closed
form is the reference at the stated accuracy.

Gated (degree 2, each aspect): (i) β — the least-squares amplitude of the
closed-form scattered field on 24 Fibonacci points on each of the 1.2× and
1.5× scaled spheroids — within `TH-15`'s band; (ii) the tip-ball mean
``Φ_tip,h = ⟨|E|²⟩/E₀²`` (cells with midpoint within r_t = 2a²/c of either
tip) against the closed form integrated by numpy quadrature over the same
cells within `MAT-4`'s 10 %; (iii) the cavity dofs are zero.  (F) and (M′)
are numpy / quadrature identities.  Control: the aspect-5 mesh driven by the
transverse trace E₀x̂ must miss both gates.
"""

from __future__ import annotations

import resource
import time

import basix
import numpy as np
import pytest
import ufl
from dolfinx import fem
from mpi4py import MPI
from petsc4py import PETSc

from fem_em_solver.core import (
    HomogeneousMaterial,
    TimeHarmonicProblem,
    TimeHarmonicSolver,
)
from fem_em_solver.io.mesh import MeshGenerator
from fem_em_solver.post.evaluation import evaluate_vector_field_parallel
from fem_em_solver.utils.analytical import (
    pec_prolate_spheroid_field,
    pec_prolate_spheroid_potential,
    prolate_spheroid_depolarisation,
)
from tests.complex_mode import complex_only
from tests.validation.test_imp1_pec_sphere_lossy_host import (
    CAVITY_TAG,
    DEGREE,
    E0,
    EPSILON_R_HOST,
    MAT4_BAND,
    OUTER_BOUNDARY_TAG,
    SPHERE_RADIUS,
    TH15_BETA_BAND,
    _assemble_real,
    _cavity_dofs,
    _phi_closed_form_on_mask,
)
from tests.validation.test_lossy_sphere_sar import FREQUENCY_HZ, SIGMA_HIGH

SEMI_MAJOR = 0.006
ASPECTS = (2.0, 5.0)
BOX_HALF_EXTENTS = (SEMI_MAJOR, SEMI_MAJOR, 2.0 * SEMI_MAJOR)
CELL_GUARD = 95_000
QDEG = 8

_MESH_CACHE: dict = {}


def _semi_minor(aspect):
    return SEMI_MAJOR / aspect


def _mesh(aspect):
    if aspect not in _MESH_CACHE:
        a = _semi_minor(aspect)
        rho_t = a * a / SEMI_MAJOR
        t0 = time.perf_counter()
        out = MeshGenerator.spheroid_in_box_domain(
            semi_major=SEMI_MAJOR,
            semi_minor=a,
            box_half_extents=BOX_HALF_EXTENTS,
            resolution_tip=rho_t / 3.0,
            resolution_surface=a / 6.0,
            resolution_far=SEMI_MAJOR / 4.0,
            comm=MPI.COMM_WORLD,
        )
        _MESH_CACHE[aspect] = (out, time.perf_counter() - t0)
    return _MESH_CACHE[aspect][0]


def _ncells(msh):
    return int(msh.comm.allreduce(msh.topology.index_map(msh.topology.dim).size_local, op=MPI.SUM))


def _rss_gib(comm):
    kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return comm.allreduce(kb, op=MPI.SUM) / 1024.0**2


def _axial_field(aspect, interior="zero"):
    a = _semi_minor(aspect)

    def field(points):
        return pec_prolate_spheroid_field(
            points, semi_major=SEMI_MAJOR, semi_minor=a, e0=E0, interior=interior
        )

    return field


def _dirichlet_axial(aspect):
    f = _axial_field(aspect, "zero")
    return lambda x: np.ascontiguousarray(f(x.T).T).astype(PETSc.ScalarType)


def _dirichlet_transverse(aspect):
    a = _semi_minor(aspect)
    c = SEMI_MAJOR

    def field(x):
        inside = (x[0] ** 2 + x[1] ** 2) / a**2 + x[2] ** 2 / c**2 < 1.0 - 1e-12
        out = np.zeros((3, x.shape[1]), dtype=PETSc.ScalarType)
        out[0] = np.where(inside, 0.0, E0)
        return out

    return field


def _probe_points(aspect):
    """24 Fibonacci points on each of the 1.2× and 1.5× scaled spheroids."""
    a = _semi_minor(aspect)
    n = 24
    i = np.arange(n) + 0.5
    pol = np.arccos(1.0 - 2.0 * i / n)
    az = np.pi * (1.0 + 5.0**0.5) * i
    out = []
    for s in (1.2, 1.5):
        out.append(
            np.column_stack(
                [
                    s * a * np.cos(az) * np.sin(pol),
                    s * a * np.sin(az) * np.sin(pol),
                    s * SEMI_MAJOR * np.cos(pol),
                ]
            )
        )
    return np.vstack(out)


def _tip_mask(msh, aspect):
    a = _semi_minor(aspect)
    rt = 2.0 * a * a / SEMI_MAJOR
    c = SEMI_MAJOR
    dg0 = fem.functionspace(msh, ("DG", 0))
    chi = fem.Function(dg0)
    chi.interpolate(
        lambda x: (
            (x[0] ** 2 + x[1] ** 2 + (x[2] - c) ** 2 < rt * rt)
            | (x[0] ** 2 + x[1] ** 2 + (x[2] + c) ** 2 < rt * rt)
        ).astype(np.float64)
    )
    return chi


def _owned_masked_cells(msh, chi):
    n_owned = msh.topology.index_map(msh.topology.dim).size_local
    vals = np.real(np.asarray(chi.x.array))
    dofmap = chi.function_space.dofmap.list
    cells = np.arange(n_owned, dtype=np.int32)
    flags = vals[dofmap[cells, 0]] > 0.5
    return cells[flags]


def _numpy_quadrature_mean(msh, cells, field_fn):
    """Σ w|E|²|J| / Σ w|J| over ``cells`` (owned), basix degree-8 rule; reduced."""
    dmap = msh.geometry.dofmap
    if dmap.shape[1] != 4:
        raise RuntimeError(f"expected affine tets, geometry dofmap {dmap.shape}")
    pts_ref, wts = basix.make_quadrature(basix.CellType.tetrahedron, QDEG)
    xg = msh.geometry.x
    v = xg[dmap[cells]]  # (n, 4, 3)
    jac = np.stack([v[:, 1] - v[:, 0], v[:, 2] - v[:, 0], v[:, 3] - v[:, 0]], axis=2)
    det = np.abs(np.linalg.det(jac)) if len(cells) else np.zeros(0)
    phys = v[:, 0][:, None, :] + np.einsum("cij,qj->cqi", jac, pts_ref)
    if len(cells):
        e = field_fn(phys.reshape(-1, 3)).reshape(len(cells), len(wts), 3)
        e2 = np.sum(np.abs(e) ** 2, axis=2)
        num = float(np.sum(det[:, None] * wts[None, :] * e2))
        vol = float(np.sum(det) * np.sum(wts))
    else:
        num = vol = 0.0
    num = msh.comm.allreduce(num, op=MPI.SUM)
    vol = msh.comm.allreduce(vol, op=MPI.SUM)
    return num / (E0**2 * vol)


def _imp1_dipole_numpy(points):
    """`IMP-1`'s β = 1 field E₀ẑ + E₀R³(3z r/r⁵ − ẑ/r³), R = 5 mm (for (M′))."""
    r2 = np.sum(points**2, axis=1)
    r = np.sqrt(r2)
    out = SPHERE_RADIUS**3 * E0 * (3.0 * points[:, 2:3] * points / r[:, None] ** 5)
    out[:, 2] += E0 - E0 * SPHERE_RADIUS**3 / r**3
    return out


def _solve_spheroid(aspect, drive, *, return_fields=False):
    comm = MPI.COMM_WORLD
    t0 = time.perf_counter()
    msh, cell_tags, facet_tags = _mesh(aspect)
    dirichlet = _dirichlet_axial(aspect) if drive == "axial" else _dirichlet_transverse(aspect)
    problem = TimeHarmonicProblem(
        mesh=msh,
        frequency_hz=FREQUENCY_HZ,
        material=HomogeneousMaterial(sigma=SIGMA_HIGH, epsilon_r=EPSILON_R_HOST),
        cell_tags=cell_tags,
        facet_tags=facet_tags,
        boundary_condition="pec_zero_tangential_a",
        dirichlet_e_field=dirichlet,
        pec_facet_tags=(OUTER_BOUNDARY_TAG, CAVITY_TAG),
    )
    solver = TimeHarmonicSolver(problem, degree=DEGREE)
    fields = solver.solve()

    points = _probe_points(aspect)
    re_vals, ok_re = evaluate_vector_field_parallel(fields.e_real, points, comm)
    im_vals, ok_im = evaluate_vector_field_parallel(fields.e_imag, points, comm)
    if not np.all(ok_re & ok_im):
        raise RuntimeError("probe points were not evaluated")
    e_re = np.real(re_vals[:, :3])
    e_im = np.real(im_vals[:, :3])
    basis = _axial_field(aspect)(points)
    basis[:, 2] -= E0
    uniform = np.zeros_like(points)
    uniform[:, 0 if drive == "transverse" else 2] = E0
    beta = float(np.sum((e_re - uniform) * basis) / np.sum(basis * basis))
    imag_ratio = float(
        np.max(np.linalg.norm(e_im, axis=1))
        / max(float(np.max(np.linalg.norm(e_re, axis=1))), 1e-300)
    )

    chi = _tip_mask(msh, aspect)
    e = fields.e_complex
    dx = ufl.Measure("dx", domain=msh, metadata={"quadrature_degree": QDEG})
    num = _assemble_real(chi * ufl.inner(e, e) * dx, comm)
    vol = _assemble_real(chi * dx, comm)
    phi_h = num / (E0**2 * vol)
    cells = _owned_masked_cells(msh, chi)
    n_mask = int(comm.allreduce(int(cells.size), op=MPI.SUM))
    phi_cf_h = _numpy_quadrature_mean(msh, cells, _axial_field(aspect, "continue"))

    cavity = _cavity_dofs(msh, facet_tags, solver.function_space())
    n_cavity = int(comm.allreduce(int(cavity.size), op=MPI.SUM))
    sol = np.asarray(e.x.array)
    local_max = float(np.max(np.abs(sol[cavity]))) if cavity.size else 0.0
    cavity_max = float(comm.allreduce(local_max, op=MPI.MAX))

    # h_tip / ρ_t: longest edge of the cells touching the tip ball (printed).
    a = _semi_minor(aspect)
    rho_t = a * a / SEMI_MAJOR
    dmap = msh.geometry.dofmap
    v = msh.geometry.x[dmap[cells]] if cells.size else np.zeros((0, 4, 3))
    if cells.size:
        edges = [
            np.linalg.norm(v[:, i] - v[:, j], axis=1) for i in range(4) for j in range(i + 1, 4)
        ]
        mid = v.mean(axis=1)
        near = (
            np.minimum(
                np.linalg.norm(mid - [0, 0, SEMI_MAJOR], axis=1),
                np.linalg.norm(mid + [0, 0, SEMI_MAJOR], axis=1),
            )
            < rho_t
        )
        h_tip_local = float(np.max(np.max(edges, axis=0)[near])) if np.any(near) else 0.0
    else:
        h_tip_local = 0.0
    h_tip = comm.allreduce(h_tip_local, op=MPI.MAX)

    out = dict(
        beta=beta,
        imag_ratio=imag_ratio,
        phi_h=phi_h,
        phi_cf_h=phi_cf_h,
        n_mask=n_mask,
        n_cavity=n_cavity,
        cavity_max=cavity_max,
        ncells=_ncells(msh),
        h_tip_over_rho=h_tip / rho_t,
        seconds=time.perf_counter() - t0,
        rss=_rss_gib(comm),
    )
    if return_fields:  # additive (EX-77, rule (a)): the example reads the solved field
        out["fields"] = fields
    return out


# --------------------------------------------------------------------- W0
@complex_only
def test_census_and_numpy():
    comm = MPI.COMM_WORLD
    c = SEMI_MAJOR
    omega = 2.0 * np.pi * FREQUENCY_HZ
    k = omega * np.sqrt(
        4e-7
        * np.pi
        * 8.8541878128e-12
        * (EPSILON_R_HOST - 1j * SIGMA_HIGH / (omega * 8.8541878128e-12))
    )
    corner = float(np.linalg.norm(BOX_HALF_EXTENTS))
    if comm.rank == 0:
        print(
            f"\n[IMP-2 step 1] |k| = {abs(k):.3f} 1/m (k = {k:.3f}); (|k|c)^2 = "
            f"{(abs(k) * c) ** 2:.4%} (third of beta band {TH15_BETA_BAND / 3:.4%}); "
            f"corner {corner * 1e3:.2f} mm, (|k|R)^2 = {(abs(k) * corner) ** 2:.4%} "
            f"(quasi-static box trace; I2a face departure 1.56e-2)"
        )

    for aspect in ASPECTS:
        a = _semi_minor(aspect)
        e = np.sqrt(1.0 - 1.0 / aspect**2)
        lz = prolate_spheroid_depolarisation(e)
        field = _axial_field(aspect, "continue")

        # (F) field vs central-difference gradient of step 0's potential.
        rng = np.random.default_rng(144)
        f = c * e
        xi0 = 1.0 / e
        xi = rng.uniform(1.2 * xi0, 4.0 * xi0, 40)
        eta = rng.uniform(-0.95, 0.95, 40)
        ang = rng.uniform(0.0, 2.0 * np.pi, 40)
        rho = f * np.sqrt((xi**2 - 1.0) * (1.0 - eta**2))
        pts = np.column_stack([rho * np.cos(ang), rho * np.sin(ang), f * xi * eta])

        def fd_grad(delta):
            g = np.zeros_like(pts)
            for ax in range(3):
                d = np.zeros(3)
                d[ax] = delta
                g[:, ax] = -(
                    pec_prolate_spheroid_potential(pts + d, semi_major=c, semi_minor=a)
                    - pec_prolate_spheroid_potential(pts - d, semi_major=c, semi_minor=a)
                ) / (2.0 * delta)
            return g

        exact = field(pts)
        d1 = 2e-3 * c
        g1 = fd_grad(d1)
        g2 = fd_grad(2.0 * d1)
        r1 = float(np.max(np.abs(exact - g1))) / E0
        r2 = float(np.max(np.abs(exact - g2))) / E0
        rich = float(np.max(np.abs(exact - (4.0 * g1 - g2) / 3.0))) / E0
        # Bound moved, with measurement (probe before W0, 2026-10-01 15:10): the
        # raw central difference at δ = 2e-3·c is truncation-bound, not a
        # gradient error — residual 1.1e-6 (aspect 2) / 6.2e-7 (aspect 5),
        # ratio 4.000 at δ/2δ (the O(δ²) signature; a wrong gradient reads
        # ratio ≈ 1).  The pre-registered 1e-7 is asserted on the Richardson
        # combination (4 D(δ) − D(2δ))/3 at the same δ, and the ratio is asserted.
        if comm.rank == 0:
            print(
                f"[IMP-2 F] aspect {aspect}: |E - FD| at d=2e-3c {r1:.3e}, 4e-3c "
                f"{r2:.3e}, ratio {r2 / r1:.3f}; Richardson {rich:.3e}"
            )
        assert 3.8 <= r2 / r1 <= 4.2
        assert rich <= 1e-7

        # Tangential E on the surface.
        n = 200
        i = np.arange(n) + 0.5
        mu = 1.0 - 2.0 * i / n
        az = np.pi * (1.0 + np.sqrt(5.0)) * i
        s = np.sqrt(1.0 - mu**2)
        surf = np.column_stack([a * s * np.cos(az), a * s * np.sin(az), c * mu])
        nrm = surf / np.array([a**2, a**2, c**2])
        nrm /= np.linalg.norm(nrm, axis=1)[:, None]
        es = field(surf)
        tang = es - np.sum(es * nrm, axis=1)[:, None] * nrm
        tmax = float(np.max(np.linalg.norm(tang, axis=1))) / E0
        tip = field(np.array([[0.0, 0.0, c], [0.0, 0.0, -c]]))
        tip_err = float(np.max(np.abs(tip[:, 2] * lz - E0))) / E0
        if comm.rank == 0:
            print(
                f"[IMP-2 F] aspect {aspect}: max |E_t| on surface {tmax:.3e}; "
                f"|E_z(tip) L_z - 1| {tip_err:.3e}; pointwise tip (1/L_z)^2 = "
                f"{1.0 / lz**2:.1f}"
            )
        assert tmax <= 1e-12
        assert tip_err <= 1e-12

        # Census + (M').
        msh, _, facet_tags = _mesh(aspect)
        ncells = _ncells(msh)
        n_cav_facets = int(comm.allreduce(int(facet_tags.find(CAVITY_TAG).size), op=MPI.SUM))
        n_wall_facets = int(
            comm.allreduce(int(facet_tags.find(OUTER_BOUNDARY_TAG).size), op=MPI.SUM)
        )
        chi = _tip_mask(msh, aspect)
        cells = _owned_masked_cells(msh, chi)
        n_mask = int(comm.allreduce(int(cells.size), op=MPI.SUM))
        mq = _numpy_quadrature_mean(msh, cells, _imp1_dipole_numpy)
        mu_ufl = _phi_closed_form_on_mask(msh, chi, 1.0)
        m_err = abs(mq / mu_ufl - 1.0)
        rss = _rss_gib(comm)  # collective: outside the rank-0 branch
        if comm.rank == 0:
            print(
                f"[IMP-2 census] aspect {aspect}: a = {a * 1e3:.3f} mm, cells {ncells}, "
                f"mesh {_MESH_CACHE[aspect][1]:.1f} s, cavity facets {n_cav_facets}, "
                f"wall facets {n_wall_facets}, tip-mask cells {n_mask}, RSS(sum) "
                f"{rss:.2f} GiB"
            )
            print(
                f"[IMP-2 M'] aspect {aspect}: numpy {mq:.15f} vs UFL {mu_ufl:.15f} "
                f"rel {m_err:.3e}"
            )
        assert ncells <= CELL_GUARD, f"cell guard: {ncells} > {CELL_GUARD}"
        assert n_cav_facets > 0 and n_wall_facets > 0 and n_mask > 0
        assert m_err <= 1e-12


# --------------------------------------------------------------------- W1
@complex_only
@pytest.mark.integration
def test_both_aspects():
    comm = MPI.COMM_WORLD
    results = {}
    for aspect in ASPECTS:
        d = _solve_spheroid(aspect, "axial")
        results[aspect] = d
        lz = prolate_spheroid_depolarisation(np.sqrt(1.0 - 1.0 / aspect**2))
        if comm.rank == 0:
            print(
                f"\n[IMP-2 step 1] aspect {aspect}: cells {d['ncells']}, {d['seconds']:.1f} s, "
                f"RSS(sum) {d['rss']:.2f} GiB, tip-mask cells {d['n_mask']}, "
                f"h_tip/rho_t {d['h_tip_over_rho']:.3f}"
            )
            print(
                f"  beta = {d['beta']:.6f}  |beta - 1| = {abs(d['beta'] - 1):.4%} "
                f"(band {TH15_BETA_BAND:.4%}); |Im E|/|Re E| = {d['imag_ratio']:.3e}"
            )
            print(
                f"  Phi_tip,h = {d['phi_h']:.6f}, Phi_tip,cf,h = {d['phi_cf_h']:.6f}, "
                f"ratio {d['phi_h'] / d['phi_cf_h']:.6f}, miss "
                f"{abs(d['phi_h'] / d['phi_cf_h'] - 1):.4%} (band {MAT4_BAND:.0%}); "
                f"pointwise tip (1/L_z)^2 = {1 / lz**2:.1f}"
            )
            print(f"  cavity dofs {d['n_cavity']}, max |E| {d['cavity_max']:.3e}")
    for aspect, d in results.items():
        assert d["n_cavity"] > 0
        assert d["cavity_max"] < 1e-12
        assert abs(d["beta"] - 1.0) <= TH15_BETA_BAND, (aspect, d["beta"], d["h_tip_over_rho"])
        assert abs(d["phi_h"] / d["phi_cf_h"] - 1.0) <= MAT4_BAND, (
            aspect,
            d["phi_h"],
            d["phi_cf_h"],
            d["h_tip_over_rho"],
        )


# --------------------------------------------------------------------- W2
@complex_only
@pytest.mark.integration
def test_transverse_control():
    comm = MPI.COMM_WORLD
    d = _solve_spheroid(5.0, "transverse")
    miss = abs(d["phi_h"] / d["phi_cf_h"] - 1.0)
    if comm.rank == 0:
        print(
            f"\n[IMP-2 control] aspect 5, transverse drive: cells {d['ncells']}, "
            f"{d['seconds']:.1f} s, RSS(sum) {d['rss']:.2f} GiB"
        )
        print(
            f"  beta = {d['beta']:.6f} (predicted ~0), |beta - 1| = {abs(d['beta'] - 1):.4f} "
            f"(band {TH15_BETA_BAND:.4%})"
        )
        print(
            f"  Phi_tip,h = {d['phi_h']:.6f} vs axial Phi_tip,cf,h = {d['phi_cf_h']:.6f}: "
            f"miss {miss:.4%} (predicted 60-100 %)"
        )
        print(f"  cavity dofs {d['n_cavity']}, max |E| {d['cavity_max']:.3e}")
    assert abs(d["beta"] - 1.0) > TH15_BETA_BAND
    assert miss > 0.10
