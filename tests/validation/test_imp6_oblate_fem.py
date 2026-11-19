"""`IMP-6` step 1 (§9 item 404): FEM for a PEC **oblate** spheroid in `IMP-1`'s
lossy host under the quasi-static **in-plane** trace, aspect 2 and 5.

Fixture: `IMP-2`'s hole route — the oblate spheroid (a = 6 mm equatorial, in
x and y; c = a/aspect along z) cut from a box of half-extents (1.75a, 1.75a, a)
by ``MeshGenerator.spheroid_in_box_domain(oblate=True)``; the outer walls carry
the closed-form tangential trace ``pec_oblate_spheroid_field(drive="inplane")``
(item 403), the cavity wall zero.  ``(|k| a)² = 1.154 %`` is under a third of
the tighter band (4.886 % / 3).

Gated (degree 2, each aspect, imported bands): (i) β — the least-squares
amplitude of the closed-form in-plane scattered field on 24 Fibonacci points
on each of the 1.2× and 1.5× scaled spheroids — within `TH-15`'s band; (ii)
the rim-ball mean ``Φ_rim,h = ⟨|E|²⟩/E₀²`` (cells with midpoint within
r_r = 2c²/a of (±a, 0, 0)) against the closed form integrated by numpy
quadrature over the same cells within `MAT-4`'s 10 %; (iii) the cavity dofs
are zero.  (M′) is a quadrature identity.  Control (W2): the aspect-5 mesh
driven by the axial trace must miss both gates.
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
    oblate_spheroid_depolarisation,
    pec_oblate_spheroid_field,
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

SEMI_MAJOR = 0.006  # a, equatorial (x, y)
ASPECTS = (2.0, 5.0)
BOX_HALF_EXTENTS = (1.75 * SEMI_MAJOR, 1.75 * SEMI_MAJOR, SEMI_MAJOR)
CELL_GUARD = 95_000
QDEG = 8

_MESH_CACHE: dict = {}


def _semi_minor(aspect):
    return SEMI_MAJOR / aspect


def _rho_r(aspect):
    c = _semi_minor(aspect)
    return c * c / SEMI_MAJOR


def _mesh(aspect):
    if aspect not in _MESH_CACHE:
        c = _semi_minor(aspect)
        rho_r = _rho_r(aspect)
        t0 = time.perf_counter()
        out = MeshGenerator.spheroid_in_box_domain(
            semi_major=SEMI_MAJOR,
            semi_minor=c,
            box_half_extents=BOX_HALF_EXTENTS,
            resolution_tip=rho_r / 3.0,
            resolution_surface=SEMI_MAJOR / 12.0,
            resolution_far=SEMI_MAJOR / 4.0,
            comm=MPI.COMM_WORLD,
            oblate=True,
        )
        _MESH_CACHE[aspect] = (out, time.perf_counter() - t0)
    return _MESH_CACHE[aspect][0]


def _ncells(msh):
    return int(msh.comm.allreduce(msh.topology.index_map(msh.topology.dim).size_local, op=MPI.SUM))


def _rss_gib(comm):
    kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return comm.allreduce(kb, op=MPI.SUM) / 1024.0**2


def _field(aspect, drive, interior="zero"):
    c = _semi_minor(aspect)

    def field(points):
        return pec_oblate_spheroid_field(
            points,
            semi_major=SEMI_MAJOR,
            semi_minor=c,
            drive=drive,
            e0=E0,
            interior=interior,
        )

    return field


def _dirichlet(aspect, drive):
    f = _field(aspect, drive, "zero")
    return lambda x: np.ascontiguousarray(f(x.T).T).astype(PETSc.ScalarType)


def _probe_points(aspect):
    """24 Fibonacci points on each of the 1.2× and 1.5× scaled oblate spheroids."""
    c = _semi_minor(aspect)
    n = 24
    i = np.arange(n) + 0.5
    pol = np.arccos(1.0 - 2.0 * i / n)
    az = np.pi * (1.0 + 5.0**0.5) * i
    out = []
    for s in (1.2, 1.5):
        out.append(
            np.column_stack(
                [
                    s * SEMI_MAJOR * np.cos(az) * np.sin(pol),
                    s * SEMI_MAJOR * np.sin(az) * np.sin(pol),
                    s * c * np.cos(pol),
                ]
            )
        )
    return np.vstack(out)


def _rim_mask(msh, aspect):
    c = _semi_minor(aspect)
    rr = 2.0 * c * c / SEMI_MAJOR
    a = SEMI_MAJOR
    dg0 = fem.functionspace(msh, ("DG", 0))
    chi = fem.Function(dg0)
    chi.interpolate(
        lambda x: (
            ((x[0] - a) ** 2 + x[1] ** 2 + x[2] ** 2 < rr * rr)
            | ((x[0] + a) ** 2 + x[1] ** 2 + x[2] ** 2 < rr * rr)
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


def _body_volume_deficit(msh, aspect):
    """(8 bx by bz − Σ air-cell volumes)/(4πa²c/3) − 1 (owned cells, reduced)."""
    n_owned = msh.topology.index_map(msh.topology.dim).size_local
    v = msh.geometry.x[msh.geometry.dofmap[np.arange(n_owned)]]
    jac = np.stack([v[:, 1] - v[:, 0], v[:, 2] - v[:, 0], v[:, 3] - v[:, 0]], axis=2)
    vol = msh.comm.allreduce(float(np.sum(np.abs(np.linalg.det(jac))) / 6.0), op=MPI.SUM)
    bx, by, bz = BOX_HALF_EXTENTS
    c = _semi_minor(aspect)
    body = 4.0 * np.pi * SEMI_MAJOR**2 * c / 3.0
    return (8.0 * bx * by * bz - vol) / body - 1.0


def _imp1_dipole_numpy(points):
    """`IMP-1`'s β = 1 field E₀ẑ + E₀R³(3z r/r⁵ − ẑ/r³), R = 5 mm (for (M′))."""
    r2 = np.sum(points**2, axis=1)
    r = np.sqrt(r2)
    out = SPHERE_RADIUS**3 * E0 * (3.0 * points[:, 2:3] * points / r[:, None] ** 5)
    out[:, 2] += E0 - E0 * SPHERE_RADIUS**3 / r**3
    return out


def _beta(e_re, points, aspect, basis_drive, uniform_axis):
    basis = _field(aspect, basis_drive)(points)
    basis[:, 0 if basis_drive == "inplane" else 2] -= E0
    uniform = np.zeros_like(points)
    uniform[:, uniform_axis] = E0
    return float(np.sum((e_re - uniform) * basis) / np.sum(basis * basis))


def _solve_oblate(aspect, drive, *, return_fields=False):
    comm = MPI.COMM_WORLD
    t0 = time.perf_counter()
    msh, cell_tags, facet_tags = _mesh(aspect)
    problem = TimeHarmonicProblem(
        mesh=msh,
        frequency_hz=FREQUENCY_HZ,
        material=HomogeneousMaterial(sigma=SIGMA_HIGH, epsilon_r=EPSILON_R_HOST),
        cell_tags=cell_tags,
        facet_tags=facet_tags,
        boundary_condition="pec_zero_tangential_a",
        dirichlet_e_field=_dirichlet(aspect, drive),
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
    # Trap (brief): the uniform part follows the *drive*; the basis is always
    # the in-plane scattered field minus E₀x̂.
    uniform_axis = 0 if drive == "inplane" else 2
    beta = _beta(e_re, points, aspect, "inplane", uniform_axis)
    beta_own = _beta(e_re, points, aspect, drive, uniform_axis)
    imag_ratio = float(
        np.max(np.linalg.norm(e_im, axis=1))
        / max(float(np.max(np.linalg.norm(e_re, axis=1))), 1e-300)
    )

    chi = _rim_mask(msh, aspect)
    e = fields.e_complex
    dx = ufl.Measure("dx", domain=msh, metadata={"quadrature_degree": QDEG})
    num = _assemble_real(chi * ufl.inner(e, e) * dx, comm)
    vol = _assemble_real(chi * dx, comm)
    phi_h = num / (E0**2 * vol)
    cells = _owned_masked_cells(msh, chi)
    n_mask = int(comm.allreduce(int(cells.size), op=MPI.SUM))
    phi_cf_h = _numpy_quadrature_mean(msh, cells, _field(aspect, "inplane", "continue"))

    cavity = _cavity_dofs(msh, facet_tags, solver.function_space())
    n_cavity = int(comm.allreduce(int(cavity.size), op=MPI.SUM))
    sol = np.asarray(e.x.array)
    local_max = float(np.max(np.abs(sol[cavity]))) if cavity.size else 0.0
    cavity_max = float(comm.allreduce(local_max, op=MPI.MAX))

    # h_rim / ρ_r: longest edge of mask cells with midpoint within ρ_r of a rim point.
    rho_r = _rho_r(aspect)
    v = msh.geometry.x[msh.geometry.dofmap[cells]] if cells.size else np.zeros((0, 4, 3))
    if cells.size:
        edges = [
            np.linalg.norm(v[:, i] - v[:, j], axis=1) for i in range(4) for j in range(i + 1, 4)
        ]
        mid = v.mean(axis=1)
        near = (
            np.minimum(
                np.linalg.norm(mid - [SEMI_MAJOR, 0, 0], axis=1),
                np.linalg.norm(mid + [SEMI_MAJOR, 0, 0], axis=1),
            )
            < rho_r
        )
        h_local = float(np.max(np.max(edges, axis=0)[near])) if np.any(near) else 0.0
    else:
        h_local = 0.0
    h_rim = comm.allreduce(h_local, op=MPI.MAX)

    out = dict(
        beta=beta,
        beta_own=beta_own,
        imag_ratio=imag_ratio,
        phi_h=phi_h,
        phi_cf_h=phi_cf_h,
        n_mask=n_mask,
        n_cavity=n_cavity,
        cavity_max=cavity_max,
        ncells=_ncells(msh),
        h_rim_over_rho=h_rim / rho_r,
        deficit=_body_volume_deficit(msh, aspect),
        seconds=time.perf_counter() - t0,
        rss=_rss_gib(comm),
    )
    if return_fields:  # additive (EX-111, rule (a)): the example reads the solved field
        out["fields"] = fields
    return out


# --------------------------------------------------------------------- W0
@complex_only
def test_census_and_numpy():
    comm = MPI.COMM_WORLD
    a = SEMI_MAJOR
    omega = 2.0 * np.pi * FREQUENCY_HZ
    eps0 = 8.8541878128e-12
    k = omega * np.sqrt(4e-7 * np.pi * eps0 * (EPSILON_R_HOST - 1j * SIGMA_HIGH / (omega * eps0)))
    corner = float(np.linalg.norm(BOX_HALF_EXTENTS))
    if comm.rank == 0:
        print(
            f"\n[IMP-6 step 1] |k| = {abs(k):.3f} 1/m; (|k|a)^2 = {(abs(k) * a) ** 2:.4%} "
            f"(third of beta band {TH15_BETA_BAND / 3:.4%}); corner {corner * 1e3:.2f} mm, "
            f"(|k|R)^2 = {(abs(k) * corner) ** 2:.4%} (I2a face departure 1.56e-2)"
        )
    for aspect in ASPECTS:
        c = _semi_minor(aspect)
        e = np.sqrt(1.0 - 1.0 / aspect**2)
        l_x, l_z = oblate_spheroid_depolarisation(e)
        msh, _, facet_tags = _mesh(aspect)
        ncells = _ncells(msh)
        n_cav = int(comm.allreduce(int(facet_tags.find(CAVITY_TAG).size), op=MPI.SUM))
        n_wall = int(comm.allreduce(int(facet_tags.find(OUTER_BOUNDARY_TAG).size), op=MPI.SUM))
        chi = _rim_mask(msh, aspect)
        cells = _owned_masked_cells(msh, chi)
        n_mask = int(comm.allreduce(int(cells.size), op=MPI.SUM))
        mq = _numpy_quadrature_mean(msh, cells, _imp1_dipole_numpy)
        mu_ufl = _phi_closed_form_on_mask(msh, chi, 1.0)
        m_err = abs(mq / mu_ufl - 1.0)
        deficit = _body_volume_deficit(msh, aspect)
        rss = _rss_gib(comm)  # collective: outside the rank-0 branch
        if comm.rank == 0:
            print(
                f"[IMP-6 census] aspect {aspect}: c = {c * 1e3:.3f} mm, rho_r = "
                f"{_rho_r(aspect) * 1e3:.3f} mm, cells {ncells}, mesh "
                f"{_MESH_CACHE[aspect][1]:.1f} s, cavity facets {n_cav}, wall facets "
                f"{n_wall}, rim-mask cells {n_mask}, body-volume deficit {deficit:.3e}, "
                f"RSS(sum) {rss:.2f} GiB; L_x {l_x:.6f}, pointwise rim (1/L_x)^2 = "
                f"{1.0 / l_x**2:.2f}"
            )
            print(
                f"[IMP-6 M'] aspect {aspect}: numpy {mq:.15f} vs UFL {mu_ufl:.15f} "
                f"rel {m_err:.3e}"
            )
        assert ncells <= CELL_GUARD, f"cell guard: {ncells} > {CELL_GUARD}"
        assert n_cav > 0 and n_wall > 0 and n_mask > 0
        assert m_err <= 1e-12


def _report_and_gate(aspect, d):
    comm = MPI.COMM_WORLD
    l_x, _ = oblate_spheroid_depolarisation(np.sqrt(1.0 - 1.0 / aspect**2))
    if comm.rank == 0:
        print(
            f"\n[IMP-6 step 1] aspect {aspect}, in-plane: cells {d['ncells']}, "
            f"{d['seconds']:.1f} s, RSS(sum) {d['rss']:.2f} GiB, rim-mask cells "
            f"{d['n_mask']}, h_rim/rho_r {d['h_rim_over_rho']:.3f}, body-volume "
            f"deficit {d['deficit']:.3e}"
        )
        print(
            f"  beta = {d['beta']:.6f}  |beta - 1| = {abs(d['beta'] - 1):.4%} "
            f"(band {TH15_BETA_BAND:.4%}, predicted <~1.5 %); |Im E|/|Re E| = "
            f"{d['imag_ratio']:.3e}"
        )
        print(
            f"  Phi_rim,h = {d['phi_h']:.6f}, Phi_rim,cf,h = {d['phi_cf_h']:.6f}, "
            f"ratio {d['phi_h'] / d['phi_cf_h']:.6f}, miss "
            f"{abs(d['phi_h'] / d['phi_cf_h'] - 1):.4%} (band {MAT4_BAND:.0%}, "
            f"predicted <~3 %); pointwise rim (1/L_x)^2 = {1 / l_x**2:.2f}"
        )
        print(f"  cavity dofs {d['n_cavity']}, max |E| {d['cavity_max']:.3e}")
    assert d["n_cavity"] > 0
    assert d["cavity_max"] < 1e-12
    assert abs(d["beta"] - 1.0) <= TH15_BETA_BAND, (aspect, d["beta"], d["h_rim_over_rho"])
    assert abs(d["phi_h"] / d["phi_cf_h"] - 1.0) <= MAT4_BAND, (
        aspect,
        d["phi_h"],
        d["phi_cf_h"],
        d["h_rim_over_rho"],
    )


# --------------------------------------------------------------------- W1a
@complex_only
@pytest.mark.integration
def test_in_plane_aspect2():
    _report_and_gate(2.0, _solve_oblate(2.0, "inplane"))


# --------------------------------------------------------------------- W1b
@complex_only
@pytest.mark.integration
def test_in_plane_aspect5():
    _report_and_gate(5.0, _solve_oblate(5.0, "inplane"))


# --------------------------------------------------------------------- W2
@complex_only
@pytest.mark.integration
def test_axial_control():
    comm = MPI.COMM_WORLD
    d = _solve_oblate(5.0, "axial")
    miss = abs(d["phi_h"] / d["phi_cf_h"] - 1.0)
    if comm.rank == 0:
        print(
            f"\n[IMP-6 control] aspect 5, axial drive: cells {d['ncells']}, "
            f"{d['seconds']:.1f} s, RSS(sum) {d['rss']:.2f} GiB"
        )
        print(
            f"  beta (in-plane basis) = {d['beta']:.6f} (predicted ~0), |beta - 1| = "
            f"{abs(d['beta'] - 1):.4f} (band {TH15_BETA_BAND:.4%}); beta_axial "
            f"(axial basis, record) = {d['beta_own']:.6f}"
        )
        print(
            f"  Phi_rim,h = {d['phi_h']:.6f} vs in-plane Phi_rim,cf,h = "
            f"{d['phi_cf_h']:.6f}: miss {miss:.4%} (predicted 90-100 %)"
        )
        print(f"  cavity dofs {d['n_cavity']}, max |E| {d['cavity_max']:.3e}")
    assert abs(d["beta"] - 1.0) > TH15_BETA_BAND
    assert miss > 0.10
