"""`IMP-1` step 1: a PEC-sphere implant in a lossy host under an imposed
quasi-static field, degree 2 (PROJECT_PLAN §9 item 91).

Outside a perfectly conducting sphere of radius ``a`` in a uniform field
``E₀ẑ`` the quasi-static field is::

    E = E₀ẑ + β E₀ a³ (3 (ẑ·r̂) r̂ − ẑ)/r³,     β = 1

(`test_pec_sphere_hole.py:8`).  Over a shell ``a < r < c`` the cross term
between the uniform part and the dipole integrates to zero, and
``∫|E_dip|² = (8π/3) β² E₀² a³ (1 − (a/c)³)``, so the shell-mean enhancement is::

    Φ(c) = ⟨|E|²⟩/E₀² = 1 + 2 β² (a/c)³ = 1.5926   at c = 1.5 a.

Since local SAR is ``σ|E|²/(2ρ)`` in a homogeneous host, Φ is the shell-mean SAR
enhancement the implant produces over the far-field SAR.

Fixture: `TH-15`'s hole route (``sphere_in_box_domain(as_hole=True)``) scaled
by 1/10 — a = 5 mm in a 10 mm half-width box — in `MAT-4`'s host (σ = 0.57 S/m,
εᵣ = 78, 64 MHz): ``|k_host| b = 0.179`` (`MAT-4`'s gated reach),
``|k| a = 0.090``.  The cavity wall (facet tag 2) and the outer wall (tag 1)
carry the closed-form tangential trace; the cavity's trace is zero.

Gated (degree 2): (i) ``|β − 1|`` within `TH-15`'s band (imported);
(ii) ``Φ_h`` against the closed form integrated over the *same* discrete mask
within `MAT-4`'s 10 %; (iii) the cavity dofs are zero.  The natural-cavity
negative control (β → −½, Φ_void = 1.1481) is asserted only to miss.
"""

from __future__ import annotations

import numpy as np
import pytest
import ufl
from mpi4py import MPI
from petsc4py import PETSc

from dolfinx import fem

from fem_em_solver.core import (
    HomogeneousMaterial,
    TimeHarmonicProblem,
    TimeHarmonicSolver,
)
from fem_em_solver.io.mesh import MeshGenerator
from fem_em_solver.post.evaluation import evaluate_vector_field_parallel
from fem_em_solver.post.sar import point_sar

from tests.complex_mode import complex_only
from tests.validation.test_lossy_sphere_sar import (
    FREQUENCY_HZ,
    RHO_KG_M3,
    SIGMA_HIGH,
)
from tests.validation.test_pec_sphere_hole import BAND as TH15_BETA_BAND

# `MAT-4`'s band.  It is a literal `0.10` at `test_lossy_sphere_sar.py:321`
# (there is no module constant to import); restated here by value, not widened.
MAT4_BAND = 0.10

SPHERE_RADIUS = 0.005
BOX_HALF_WIDTH = 0.010
RESOLUTION_SPHERE = SPHERE_RADIUS / 6.0
RESOLUTION_FAR = SPHERE_RADIUS / 3.0
EPSILON_R_HOST = 78.0
E0 = 1.0

OUTER_BOUNDARY_TAG = 1
CAVITY_TAG = 2

BETA_PEC = 1.0
BETA_VOID = -0.5
SHELL_OUTER_FACTOR = 1.5
PHI_CLOSED_FORM_PEC = 1.0 + 2.0 * BETA_PEC**2 / SHELL_OUTER_FACTOR**3  # 1.5926
PHI_CLOSED_FORM_VOID = 1.0 + 2.0 * BETA_VOID**2 / SHELL_OUTER_FACTOR**3  # 1.1481

DEGREE = 2


def _pec_exterior_numpy():
    """β = 1 closed form at radius ``a``, zero for r < a (degree-2 interpolation
    points on the cavity chords lie strictly inside the sphere)."""
    a = SPHERE_RADIUS

    def field(x: np.ndarray) -> np.ndarray:
        r2 = x[0] ** 2 + x[1] ** 2 + x[2] ** 2
        r2_safe = np.where(r2 > 0.0, r2, 1.0)
        r5 = r2_safe**2.5
        z = x[2]
        common = BETA_PEC * E0 * a**3 / r5
        ex_out = common * 3.0 * z * x[0]
        ey_out = common * 3.0 * z * x[1]
        ez_out = E0 + common * (3.0 * z * z - r2_safe)
        inside = r2 < a * a
        return np.array(
            [
                np.where(inside, 0.0, ex_out),
                np.where(inside, 0.0, ey_out),
                np.where(inside, 0.0, ez_out),
            ],
            dtype=PETSc.ScalarType,
        )

    return field


def _dipole_basis(points: np.ndarray) -> np.ndarray:
    r = np.linalg.norm(points, axis=1)
    r_hat = points / r[:, None]
    z_hat = np.array([0.0, 0.0, 1.0])
    cos_theta = r_hat[:, 2]
    return (
        E0
        * SPHERE_RADIUS**3
        / r[:, None] ** 3
        * (3.0 * cos_theta[:, None] * r_hat - z_hat[None, :])
    )


def _probe_shells() -> np.ndarray:
    """24 Fibonacci points on each of r = 1.2a and 1.5a (`test_pec_sphere_hole.py:145`)."""
    n_per_shell = 24
    indices = np.arange(n_per_shell) + 0.5
    phi = np.arccos(1.0 - 2.0 * indices / n_per_shell)
    theta = np.pi * (1.0 + 5.0**0.5) * indices
    shells = []
    for factor in (1.2, 1.5):
        rr = factor * SPHERE_RADIUS
        shells.append(
            np.column_stack(
                [
                    rr * np.cos(theta) * np.sin(phi),
                    rr * np.sin(theta) * np.sin(phi),
                    rr * np.cos(phi),
                ]
            )
        )
    return np.vstack(shells)


_MESH_CACHE: dict = {}


def _hole_mesh(resolution_sphere=RESOLUTION_SPHERE, resolution_far=RESOLUTION_FAR):
    """Step 1's mesh by default; step 2 passes its own rung (cache keyed by the pair)."""
    key = (resolution_sphere, resolution_far)
    if key not in _MESH_CACHE:
        _MESH_CACHE[key] = MeshGenerator.sphere_in_box_domain(
            sphere_radius=SPHERE_RADIUS,
            box_half_width=BOX_HALF_WIDTH,
            resolution_sphere=resolution_sphere,
            resolution_far=resolution_far,
            comm=MPI.COMM_WORLD,
            as_hole=True,
        )
    return _MESH_CACHE[key]


def _shell_mask(msh):
    """DG0 indicator of cells whose midpoint has r² < c² (c = 1.5a); numpy only."""
    c = SHELL_OUTER_FACTOR * SPHERE_RADIUS
    dg0 = fem.functionspace(msh, ("DG", 0))
    chi = fem.Function(dg0)
    chi.interpolate(
        lambda x: ((x[0] ** 2 + x[1] ** 2 + x[2] ** 2) < c * c).astype(np.float64)
    )
    return chi


def _assemble_real(form_expr, comm) -> float:
    local = fem.assemble_scalar(fem.form(form_expr))
    return float(comm.allreduce(float(np.real(local)), op=MPI.SUM))


def _phi_closed_form_on_mask(msh, chi, beta: float) -> float:
    """The closed-form ⟨|E|²⟩/E₀² over the same discrete mask (UFL, qdeg 8)."""
    comm = msh.comm
    x = ufl.SpatialCoordinate(msh)
    r2 = ufl.dot(x, x)
    r = ufl.sqrt(r2)
    z_hat = ufl.as_vector((0.0, 0.0, 1.0))
    e_cf = E0 * z_hat + beta * E0 * SPHERE_RADIUS**3 * (
        3.0 * x[2] * x / r**5 - z_hat / r**3
    )
    dx = ufl.Measure("dx", domain=msh, metadata={"quadrature_degree": 8})
    num = _assemble_real(chi * ufl.dot(e_cf, e_cf) * dx, comm)
    vol = _assemble_real(chi * dx, comm)
    return num / (E0**2 * vol)


def _cavity_dofs(msh, facet_tags, v_space) -> np.ndarray:
    fdim = msh.topology.dim - 1
    msh.topology.create_connectivity(fdim, msh.topology.dim)
    facets = np.asarray(facet_tags.find(CAVITY_TAG), dtype=np.int32)
    return fem.locate_dofs_topological(v_space, fdim, facets)


def _solve(pec_facet_tags, degree, resolution_sphere=RESOLUTION_SPHERE,
           resolution_far=RESOLUTION_FAR, dirichlet_e_field=None):
    """``dirichlet_e_field=None`` is the quasi-static trace ``_pec_exterior_numpy()``
    (byte-identical to the pre-step-3(b) behaviour); step 3(b) passes the exact
    finite-k trace (§9 item 146)."""
    comm = MPI.COMM_WORLD
    if dirichlet_e_field is None:
        dirichlet_e_field = _pec_exterior_numpy()
    msh, cell_tags, facet_tags = _hole_mesh(resolution_sphere, resolution_far)
    problem = TimeHarmonicProblem(
        mesh=msh,
        frequency_hz=FREQUENCY_HZ,
        material=HomogeneousMaterial(sigma=SIGMA_HIGH, epsilon_r=EPSILON_R_HOST),
        cell_tags=cell_tags,
        facet_tags=facet_tags,
        boundary_condition="pec_zero_tangential_a",
        dirichlet_e_field=dirichlet_e_field,
        pec_facet_tags=pec_facet_tags,
    )
    solver = TimeHarmonicSolver(problem, degree=degree)
    fields = solver.solve()

    points = _probe_shells()
    real_values, valid_real = evaluate_vector_field_parallel(fields.e_real, points, comm)
    imag_values, valid_imag = evaluate_vector_field_parallel(fields.e_imag, points, comm)
    if not np.all(valid_real & valid_imag):
        raise RuntimeError("exterior probe points were not evaluated")
    e_re = np.real(real_values[:, :3])
    e_im = np.real(imag_values[:, :3])
    basis = _dipole_basis(points)
    uniform = np.zeros_like(points)
    uniform[:, 2] = E0
    beta = float(np.sum((e_re - uniform) * basis) / np.sum(basis * basis))
    imag_ratio = float(
        np.max(np.linalg.norm(e_im, axis=1))
        / max(float(np.max(np.linalg.norm(e_re, axis=1))), 1e-300)
    )

    chi = _shell_mask(msh)
    e = fields.e_complex
    dx = ufl.Measure("dx", domain=msh, metadata={"quadrature_degree": 8})
    num = _assemble_real(chi * ufl.inner(e, e) * dx, comm)
    vol = _assemble_real(chi * dx, comm)
    phi_h = num / (E0**2 * vol)
    phi_cf_h = _phi_closed_form_on_mask(msh, chi, BETA_PEC)

    cavity_dofs = _cavity_dofs(msh, facet_tags, solver.function_space())
    n_cavity = int(comm.allreduce(int(cavity_dofs.size), op=MPI.SUM))
    solution = np.asarray(e.x.array)
    local_max = float(np.max(np.abs(solution[cavity_dofs]))) if cavity_dofs.size else 0.0
    cavity_max = float(comm.allreduce(local_max, op=MPI.MAX))

    ncells = int(
        comm.allreduce(msh.topology.index_map(msh.topology.dim).size_local, op=MPI.SUM)
    )
    return dict(
        beta=beta,
        imag_ratio=imag_ratio,
        phi_h=phi_h,
        phi_cf_h=phi_cf_h,
        n_cavity=n_cavity,
        cavity_max=cavity_max,
        ncells=ncells,
        fields=fields,
    )


def _sar_ratio_records(fields):
    """SAR/SAR₀ at r = 1.2a, polar 3° and 87° (printed only), with closed forms."""
    rr = 1.2 * SPHERE_RADIUS
    azimuth = 0.3
    out = []
    for theta_deg in (3.0, 87.0):
        th = np.deg2rad(theta_deg)
        p = np.array(
            [[rr * np.sin(th) * np.cos(azimuth), rr * np.sin(th) * np.sin(azimuth), rr * np.cos(th)]]
        )
        sar = point_sar(
            fields.e_real, fields.e_imag, p, sigma=SIGMA_HIGH, rho=RHO_KG_M3
        )
        sar0 = SIGMA_HIGH * E0**2 / (2.0 * RHO_KG_M3)
        q = (SPHERE_RADIUS / rr) ** 3
        closed = (np.cos(th) * (1.0 + 2.0 * q)) ** 2 + (np.sin(th) * (1.0 - q)) ** 2
        out.append((theta_deg, float(np.real(np.ravel(sar)[0])) / sar0, float(closed)))
    return out


@complex_only
@pytest.mark.integration
def test_pec_implant_in_lossy_host_beta_and_shell_sar_enhancement():
    comm = MPI.COMM_WORLD
    d2 = _solve((OUTER_BOUNDARY_TAG, CAVITY_TAG), DEGREE)
    sar_records = _sar_ratio_records(d2["fields"])
    d1 = _solve((OUTER_BOUNDARY_TAG, CAVITY_TAG), 1)

    beta_miss = abs(d2["beta"] - BETA_PEC)
    phi_miss = abs(d2["phi_h"] / d2["phi_cf_h"] - 1.0)
    phi_miss_d1 = abs(d1["phi_h"] / d1["phi_cf_h"] - 1.0)

    if comm.rank == 0:
        omega = 2.0 * np.pi * FREQUENCY_HZ
        print(
            f"\n[IMP-1 step 1] PEC sphere a = {SPHERE_RADIUS} m in box b = "
            f"{BOX_HALF_WIDTH} m, host sigma = {SIGMA_HIGH} S/m, eps_r = "
            f"{EPSILON_R_HOST}, f = {FREQUENCY_HZ / 1e6:.3f} MHz "
            f"(omega = {omega:.4e}), {d2['ncells']} cells"
        )
        print(
            f"  degree 2: beta = {d2['beta']:.6f}  |beta - 1| = {beta_miss:.4%} "
            f"(band {TH15_BETA_BAND:.4%})"
        )
        print(
            f"  degree 2: Phi_h = {d2['phi_h']:.6f}, Phi_cf,h = {d2['phi_cf_h']:.6f}, "
            f"|Phi_h/Phi_cf,h - 1| = {phi_miss:.4%} (band {MAT4_BAND:.0%})"
        )
        print(
            f"  [record] Phi_cf,h = {d2['phi_cf_h']:.6f} vs continuous closed form "
            f"{PHI_CLOSED_FORM_PEC:.6f} ({d2['phi_cf_h'] / PHI_CLOSED_FORM_PEC - 1:+.4%})"
        )
        print(f"  [record] degree 2 |Im E|/|Re E| = {d2['imag_ratio']:.4e}")
        print(
            f"  [record] degree 1 same mesh: beta = {d1['beta']:.6f} "
            f"(|beta-1| = {abs(d1['beta'] - 1):.4%}), Phi_h = {d1['phi_h']:.6f}, "
            f"miss {phi_miss_d1:.4%}, |Im E|/|Re E| = {d1['imag_ratio']:.4e}"
        )
        for theta_deg, ratio, closed in sar_records:
            print(
                f"  [record] SAR/SAR0 at r = 1.2a, polar {theta_deg:.0f} deg: "
                f"{ratio:.4f} (quasi-static closed form {closed:.4f})"
            )
        print(
            f"  cavity dofs (degree 2, reduced) = {d2['n_cavity']}, max |E| = "
            f"{d2['cavity_max']:.3e}"
        )

    assert d2["n_cavity"] > 0, "no dofs located on cavity facet tag 2"
    assert d2["cavity_max"] < 1e-12, (
        f"degree-2 cavity dofs carry max |E| = {d2['cavity_max']:.3e}: the hole "
        "route leaks at degree 2"
    )
    assert beta_miss <= TH15_BETA_BAND, (
        f"degree-2 beta = {d2['beta']:.6f} misses 1 by {beta_miss:.4%} > "
        f"{TH15_BETA_BAND:.4%} (degree 1: {d1['beta']:.6f}; |Im E|/|Re E| = "
        f"{d2['imag_ratio']:.3e})"
    )
    assert phi_miss <= MAT4_BAND, (
        f"degree-2 Phi_h = {d2['phi_h']:.6f} vs Phi_cf,h = {d2['phi_cf_h']:.6f}: "
        f"miss {phi_miss:.4%} > {MAT4_BAND:.0%} (degree 1 miss {phi_miss_d1:.4%})"
    )


@complex_only
@pytest.mark.integration
def test_natural_cavity_negative_control():
    """Drop tag 2 from the Dirichlet set: the cavity is a void, β → −½,
    Φ_void = 1.1481 (27.9 % below the PEC closed form, also the ceiling)."""
    comm = MPI.COMM_WORLD
    ctl = _solve((OUTER_BOUNDARY_TAG,), DEGREE)
    delta_beta = abs(ctl["beta"] - BETA_PEC)
    phi_miss = abs(ctl["phi_h"] / ctl["phi_cf_h"] - 1.0)
    if comm.rank == 0:
        print(
            f"\n[IMP-1 step 1] negative control, pec_facet_tags=(1,), degree 2: "
            f"beta = {ctl['beta']:.6f} (void {BETA_VOID:+.2f}), |dbeta| = "
            f"{delta_beta:.4f} (predicted 1.5)"
        )
        print(
            f"  Phi_h = {ctl['phi_h']:.6f} (void closed form {PHI_CLOSED_FORM_VOID:.4f}) "
            f"vs Phi_cf,h(PEC) = {ctl['phi_cf_h']:.6f}: miss {phi_miss:.4%} "
            f"(predicted 27.9 %)"
        )
    assert delta_beta > TH15_BETA_BAND, (
        f"natural cavity beta = {ctl['beta']:.6f} within the band of 1"
    )
    assert phi_miss > MAT4_BAND, (
        f"natural cavity Phi miss {phi_miss:.4%} not over {MAT4_BAND:.0%}: the Phi "
        "gate does not resolve the PEC condition"
    )
