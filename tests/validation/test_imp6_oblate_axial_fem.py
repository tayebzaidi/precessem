"""`IMP-6` step 3 (§9 item 409): the **axial** drive gated at aspect 2 and 5 on
step 1's meshes, against the pole enhancement E_z = E₀/L_z at (0, 0, ±c).

Fixture: step 1's (`test_imp6_oblate_fem.py`) oblate spheroid in `IMP-1`'s lossy
host, a = 6 mm, aspect 2 / 5, same meshes (cached by aspect, unrefined at the
pole).  The outer walls carry the closed-form tangential trace
``pec_oblate_spheroid_field(drive="axial")``, the cavity wall zero.

Gated (degree 2, each aspect, imported bands): (i) β_axial — least-squares
amplitude of the closed-form *axial* scattered field, uniform part on z — within
`TH-15`'s band; (ii) the pole-ball mean ``Φ_pole,h = ⟨|E|²⟩/E₀²`` (cells with
midpoint within r_p = a/5 of (0, 0, ±c)) against the axial closed form
integrated by numpy quadrature over the same cells, within `MAT-4`'s band;
(iii) cavity dofs zero.  (M″) is the pole-mask quadrature identity.  Control
(W2): aspect-5 mesh, in-plane drive, compared on the axial basis and against
the axial Φ_cf — must miss both gates (all factors predicted).

Only private names are imported from step 1's module (rule: never import a
``test_*`` name; pytest would collect and re-run step 1's solves).
"""

from __future__ import annotations

import time

import numpy as np
import pytest
import ufl
from dolfinx import fem
from mpi4py import MPI

from fem_em_solver.core import (
    HomogeneousMaterial,
    TimeHarmonicProblem,
    TimeHarmonicSolver,
)
from fem_em_solver.post.evaluation import evaluate_vector_field_parallel
from fem_em_solver.utils.analytical import oblate_spheroid_depolarisation
from tests.complex_mode import complex_only
from tests.validation.test_imp1_pec_sphere_lossy_host import (
    CAVITY_TAG,
    DEGREE,
    E0,
    EPSILON_R_HOST,
    MAT4_BAND,
    OUTER_BOUNDARY_TAG,
    TH15_BETA_BAND,
    _assemble_real,
    _cavity_dofs,
    _phi_closed_form_on_mask,
)
from tests.validation.test_imp6_oblate_fem import (
    CELL_GUARD,
    QDEG,
    SEMI_MAJOR,
    _beta,
    _body_volume_deficit,
    _dirichlet,
    _field,
    _imp1_dipole_numpy,
    _mesh,
    _ncells,
    _numpy_quadrature_mean,
    _owned_masked_cells,
    _probe_points,
    _rss_gib,
    _semi_minor,
)
from tests.validation.test_lossy_sphere_sar import FREQUENCY_HZ, SIGMA_HIGH

ASPECTS = (2.0, 5.0)
POLE_RADIUS = SEMI_MAJOR / 5.0  # r_p = a/5, fixed by the brief before any solve
STEP1_CELLS = {2.0: 39_142, 5.0: 32_676}  # printed beside, not asserted
# Rule (e): stdlib Monte Carlo of item 403's closed form over the continuous
# ball (brief, item 409); printed beside the measured value, never asserted.
PREDICTED_PHI_CF = {2.0: 3.11, 5.0: 1.67}


def _pole_mask(msh, aspect):
    c = _semi_minor(aspect)
    rp2 = POLE_RADIUS * POLE_RADIUS
    dg0 = fem.functionspace(msh, ("DG", 0))
    chi = fem.Function(dg0)
    chi.interpolate(
        lambda x: (
            (x[0] ** 2 + x[1] ** 2 + (x[2] - c) ** 2 < rp2)
            | (x[0] ** 2 + x[1] ** 2 + (x[2] + c) ** 2 < rp2)
        ).astype(np.float64)
    )
    return chi


def _h_pole(msh, cells, comm):
    """Longest edge over the owned pole-mask cells (reduced MAX)."""
    if cells.size:
        v = msh.geometry.x[msh.geometry.dofmap[cells]]
        edges = [
            np.linalg.norm(v[:, i] - v[:, j], axis=1) for i in range(4) for j in range(i + 1, 4)
        ]
        h_local = float(np.max(edges))
    else:
        h_local = 0.0
    return comm.allreduce(h_local, op=MPI.MAX)


def _solve_pole(aspect, drive):
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
    # Trap (brief): the basis follows the comparison (axial, always); the
    # uniform part follows the drive.
    uniform_axis = 2 if drive == "axial" else 0
    beta = _beta(e_re, points, aspect, "axial", uniform_axis)
    imag_ratio = float(
        np.max(np.linalg.norm(e_im, axis=1))
        / max(float(np.max(np.linalg.norm(e_re, axis=1))), 1e-300)
    )

    chi = _pole_mask(msh, aspect)
    e = fields.e_complex
    dx = ufl.Measure("dx", domain=msh, metadata={"quadrature_degree": QDEG})
    num = _assemble_real(chi * ufl.inner(e, e) * dx, comm)
    vol = _assemble_real(chi * dx, comm)
    phi_h = num / (E0**2 * vol)
    cells = _owned_masked_cells(msh, chi)
    n_mask = int(comm.allreduce(int(cells.size), op=MPI.SUM))
    # Comparand: the *axial* closed form for both drives.
    phi_cf_h = _numpy_quadrature_mean(msh, cells, _field(aspect, "axial", "continue"))

    cavity = _cavity_dofs(msh, facet_tags, solver.function_space())
    n_cavity = int(comm.allreduce(int(cavity.size), op=MPI.SUM))
    sol = np.asarray(e.x.array)
    local_max = float(np.max(np.abs(sol[cavity]))) if cavity.size else 0.0
    cavity_max = float(comm.allreduce(local_max, op=MPI.MAX))

    return dict(
        beta=beta,
        imag_ratio=imag_ratio,
        phi_h=phi_h,
        phi_cf_h=phi_cf_h,
        n_mask=n_mask,
        mask_vol=vol,
        n_cavity=n_cavity,
        cavity_max=cavity_max,
        ncells=_ncells(msh),
        h_pole=_h_pole(msh, cells, comm),
        deficit=_body_volume_deficit(msh, aspect),
        seconds=time.perf_counter() - t0,
        rss=_rss_gib(comm),
    )


# --------------------------------------------------------------------- W0
@complex_only
def test_pole_census():
    comm = MPI.COMM_WORLD
    for aspect in ASPECTS:
        c = _semi_minor(aspect)
        _, l_z = oblate_spheroid_depolarisation(np.sqrt(1.0 - 1.0 / aspect**2))
        t0 = time.perf_counter()
        msh, _, facet_tags = _mesh(aspect)
        t_mesh = time.perf_counter() - t0
        ncells = _ncells(msh)
        chi = _pole_mask(msh, aspect)
        cells = _owned_masked_cells(msh, chi)
        n_mask = int(comm.allreduce(int(cells.size), op=MPI.SUM))
        mq = _numpy_quadrature_mean(msh, cells, _imp1_dipole_numpy)
        mu_ufl = _phi_closed_form_on_mask(msh, chi, 1.0)
        m_err = abs(mq / mu_ufl - 1.0)
        phi_cf = _numpy_quadrature_mean(msh, cells, _field(aspect, "axial", "continue"))
        dx = ufl.Measure("dx", domain=msh, metadata={"quadrature_degree": QDEG})
        mask_vol = _assemble_real(chi * dx, comm)
        h_pole = _h_pole(msh, cells, comm)
        deficit = _body_volume_deficit(msh, aspect)
        rss = _rss_gib(comm)  # collective: outside the rank-0 branch
        if comm.rank == 0:
            print(
                f"\n[IMP-6 step3 census] aspect {aspect}: c = {c * 1e3:.3f} mm, r_p = "
                f"{POLE_RADIUS * 1e3:.3f} mm, cells {ncells} (step 1: "
                f"{STEP1_CELLS[aspect]}, {'match' if ncells == STEP1_CELLS[aspect] else 'MISMATCH'}), "
                f"mesh {t_mesh:.1f} s, pole-mask cells {n_mask}, mask volume "
                f"{mask_vol * 1e9:.3f} mm^3 (both poles), h_pole {h_pole * 1e3:.3f} mm, "
                f"body-volume deficit {deficit:.3e}, RSS(sum) {rss:.2f} GiB"
            )
            print(
                f"[IMP-6 step3 M''] aspect {aspect}: numpy {mq:.15f} vs UFL {mu_ufl:.15f} "
                f"rel {m_err:.3e}"
            )
            print(
                f"[IMP-6 step3 cf] aspect {aspect}: L_z {l_z:.6f}, Phi_pole,cf,h = {phi_cf:.6f} "
                f"(predicted MC ~{PREDICTED_PHI_CF[aspect]:.2f}); pointwise (1/L_z)^2 = "
                f"{1.0 / l_z**2:.3f}"
            )
        assert ncells <= CELL_GUARD, f"cell guard: {ncells} > {CELL_GUARD}"
        assert n_mask > 0
        assert m_err <= 1e-12


def _report_and_gate(aspect, d):
    comm = MPI.COMM_WORLD
    _, l_z = oblate_spheroid_depolarisation(np.sqrt(1.0 - 1.0 / aspect**2))
    miss = d["phi_h"] / d["phi_cf_h"] - 1.0
    if comm.rank == 0:
        print(
            f"\n[IMP-6 step 3] aspect {aspect}, axial: cells {d['ncells']}, "
            f"{d['seconds']:.1f} s, RSS(sum) {d['rss']:.2f} GiB, pole-mask cells "
            f"{d['n_mask']}, h_pole {d['h_pole'] * 1e3:.3f} mm, body-volume deficit "
            f"{d['deficit']:.3e}"
        )
        print(
            f"  beta_axial = {d['beta']:.6f}  |beta - 1| = {abs(d['beta'] - 1):.4%} "
            f"(band {TH15_BETA_BAND:.4%}, predicted <~1.5 %); |Im E|/|Re E| = "
            f"{d['imag_ratio']:.3e}"
        )
        print(
            f"  Phi_pole,h = {d['phi_h']:.6f}, Phi_pole,cf,h = {d['phi_cf_h']:.6f} "
            f"(predicted MC ~{PREDICTED_PHI_CF[aspect]:.2f}), ratio "
            f"{d['phi_h'] / d['phi_cf_h']:.6f}, signed miss {miss:+.4%} (band "
            f"{MAT4_BAND:.0%}, predicted <~2 %); pointwise (1/L_z)^2 = {1 / l_z**2:.3f}"
        )
        print(f"  cavity dofs {d['n_cavity']}, max |E| {d['cavity_max']:.3e}")
    assert d["n_cavity"] > 0
    assert d["cavity_max"] < 1e-12
    assert abs(d["beta"] - 1.0) <= TH15_BETA_BAND, (aspect, d["beta"], d["h_pole"])
    assert abs(miss) <= MAT4_BAND, (aspect, d["phi_h"], d["phi_cf_h"], d["h_pole"])


# --------------------------------------------------------------------- W1a
@complex_only
@pytest.mark.integration
def test_axial_aspect2():
    _report_and_gate(2.0, _solve_pole(2.0, "axial"))


# --------------------------------------------------------------------- W1b
@complex_only
@pytest.mark.integration
def test_axial_aspect5():
    _report_and_gate(5.0, _solve_pole(5.0, "axial"))


# --------------------------------------------------------------------- W2
@complex_only
@pytest.mark.integration
def test_inplane_control():
    comm = MPI.COMM_WORLD
    d = _solve_pole(5.0, "inplane")
    miss = abs(d["phi_h"] / d["phi_cf_h"] - 1.0)
    if comm.rank == 0:
        print(
            f"\n[IMP-6 step3 control] aspect 5, in-plane drive: cells {d['ncells']}, "
            f"{d['seconds']:.1f} s, RSS(sum) {d['rss']:.2f} GiB, pole-mask cells {d['n_mask']}"
        )
        print(
            f"  beta (axial basis, uniform on x) = {d['beta']:.6f} (predicted ~0; step 1's "
            f"mirror read -0.001796), |beta - 1| = {abs(d['beta'] - 1):.4f} "
            f"(band {TH15_BETA_BAND:.4%})"
        )
        print(
            f"  Phi_pole,h = {d['phi_h']:.6f} (predicted MC ~0.038) vs axial "
            f"Phi_pole,cf,h = {d['phi_cf_h']:.6f}: miss {miss:.4%} (predicted 90-100 %, "
            f"MC ~97.7 %)"
        )
        print(f"  cavity dofs {d['n_cavity']}, max |E| {d['cavity_max']:.3e}")
    assert abs(d["beta"] - 1.0) > TH15_BETA_BAND
    assert miss > MAT4_BAND
