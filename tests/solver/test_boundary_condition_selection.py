"""Boundary-condition selection tests for time-harmonic solver API."""

from __future__ import annotations

import numpy as np
import pytest
import ufl
from dolfinx import fem
from dolfinx import mesh as dmesh
from mpi4py import MPI

from fem_em_solver.core import (
    HomogeneousMaterial,
    TimeHarmonicBoundaryCondition,
    TimeHarmonicProblem,
    TimeHarmonicSolver,
    normalize_boundary_condition,
)
from fem_em_solver.io.mesh import MeshGenerator

from fem_em_solver.utils.constants import EPSILON_0, MU_0

from tests.complex_mode import complex_only
from tests.validation.test_time_harmonic_mms import (
    BOX_L,
    EPSILON_R,
    FREQUENCY_HZ,
    MU_R,
    SIGMA,
)

# GEO-23 step 2c (2026-08-28): the global cell count step 2b only printed,
# now asserted.  Reference is step 1's -n 1 ladder on this generator, measured
# on the dolfinx 0.11 image; the 0.7.2-era sizing (resolution=0.04) no longer
# meshes at all, so there is no older count to compare against.  Step 2b
# reproduced 1213 exactly (0.00%) at -n 1 and -n 2.
N_CELLS_REF = 1213


def _make_problem(boundary_condition: TimeHarmonicBoundaryCondition | str) -> TimeHarmonicProblem:
    # GEO-23 step 2b (2026-08-28): 0.04 does not mesh on dolfinx 0.11 — gmsh
    # throws "Invalid boundary mesh (overlapping facets)".  Step 1's -n 1
    # ladder on this generator is monotone and its coarsest meshing rung is
    # 0.032 (1 213 cells); the finer rungs are 0.0256 (1 769) / 0.02048
    # (2 478) / 0.016384 (3 834).  Moved to the coarsest meshing rung.
    mesh, cell_tags, facet_tags = MeshGenerator.cylindrical_domain(
        inner_radius=0.01,
        outer_radius=0.08,
        length=0.12,
        resolution=0.032,
        comm=MPI.COMM_WORLD,
    )
    # size_global is already a global quantity on every rank — do not sum
    # size_local on top of it.
    n_global = mesh.topology.index_map(mesh.topology.dim).size_global
    if MPI.COMM_WORLD.rank == 0:
        print(f"GEO-23 step 2b: cylindrical_domain(resolution=0.032) -> {n_global} cells")
    assert abs(n_global / N_CELLS_REF - 1) <= 0.01, (
        f"cylindrical_domain(resolution=0.032) meshed {n_global} cells, "
        f"outside +/-1% of the GEO-23 step 1 ladder reference {N_CELLS_REF}"
    )
    return TimeHarmonicProblem(
        mesh=mesh,
        frequency_hz=127.74e6,
        material=HomogeneousMaterial(sigma=0.7, epsilon_r=78.0, mu_r=1.0),
        cell_tags=cell_tags,
        facet_tags=facet_tags,
        boundary_condition=boundary_condition,
    )


def test_normalize_boundary_condition_accepts_enum_and_string_values():
    assert (
        normalize_boundary_condition(TimeHarmonicBoundaryCondition.NATURAL)
        == TimeHarmonicBoundaryCondition.NATURAL
    )
    assert normalize_boundary_condition("pec_zero_tangential_a") == TimeHarmonicBoundaryCondition.PEC_ZERO_TANGENTIAL_A


def test_normalize_boundary_condition_rejects_unknown_value():
    with pytest.raises(ValueError, match="boundary_condition"):
        normalize_boundary_condition("perfectly_matched_layer")


def test_time_harmonic_solver_boundary_natural_selects_empty_dirichlet_set():
    solver = TimeHarmonicSolver(_make_problem(TimeHarmonicBoundaryCondition.NATURAL), degree=1)

    bcs, selected_bc, dirichlet_dof_count = solver.build_boundary_conditions()

    assert selected_bc == TimeHarmonicBoundaryCondition.NATURAL
    assert bcs == []
    assert dirichlet_dof_count == 0


# ---------------------------------------------------------------------------
# TH-3 (§9 item 89): each boundary-condition option solves to its own
# manufactured field, and zero-data PEC constrains exactly the exterior edges.
# Replaces the placeholder `..._pec_is_applied_to_solve_path` (len(bcs) == 1,
# dirichlet_dof_count > 0, unchecked solve).  Box, frequency and material are
# TH-1's (imported, never restated); bands are TH-1's
# (`test_time_harmonic_mms.py:187-197`): relative L2 at 24^3 < 0.10, 24^3 below
# 12^3, rate > 0.9.
#
# E_pec = (sin ky sin kz, sin kz sin kx, sin kx sin ky): zero tangential trace
#   on every face of [0, L]^3, divergence-free, curl curl E = 2k^2 E.
# E_nat = (sin kx cos ky, -cos kx sin ky, 0) = curl(psi z): (curl E)_z =
#   2k sin kx sin ky, so n x curl E = 0 on every face (PMC, the NATURAL
#   builder's docstring), divergence-free, curl curl E = 2k^2 E.
# Source: -j w mu0 J = (2k^2/mu_r - k0^2 eps_c) E  (TH-1's with k^2 -> 2k^2).
# ---------------------------------------------------------------------------

TH3_SIZES = (12, 24)
#: TH-1's bands, named (EX-64, rule (a)); values unchanged from the literals.
TH3_REL_L2_BAND = 0.10
TH3_RATE_FLOOR = 0.9


def _th3_field_ufl(field: str, x):
    k = np.pi / BOX_L
    if field == "pec":
        return ufl.as_vector(
            [
                ufl.sin(k * x[1]) * ufl.sin(k * x[2]),
                ufl.sin(k * x[2]) * ufl.sin(k * x[0]),
                ufl.sin(k * x[0]) * ufl.sin(k * x[1]),
            ]
        )
    if field == "nat":
        return ufl.as_vector(
            [
                ufl.sin(k * x[0]) * ufl.cos(k * x[1]),
                -ufl.cos(k * x[0]) * ufl.sin(k * x[1]),
                0.0 * x[2],
            ]
        )
    raise ValueError(field)


def _th3_rhs_scale() -> complex:
    k = np.pi / BOX_L
    omega = 2.0 * np.pi * FREQUENCY_HZ
    k0_squared = omega * omega * MU_0 * EPSILON_0
    epsilon_c = EPSILON_R - 1j * SIGMA / (omega * EPSILON_0)
    return (2.0 * k * k / MU_R) - k0_squared * epsilon_c


def _th3_current(field: str):
    omega = 2.0 * np.pi * FREQUENCY_HZ
    scale = _th3_rhs_scale() / (-1j * omega * MU_0)

    def current_density(x):
        return scale * _th3_field_ufl(field, x)

    return current_density


def _solve_mms(mode: str, field: str, n: int, return_fields: bool = False):
    """Solve ``field``'s manufactured source under boundary ``mode`` on an n^3 box.

    Returns ``(relative_l2, owned_constrained_dofs_global, ncells_global)``;
    with ``return_fields=True`` (EX-64, additive) also ``msh, fields``.
    """
    comm = MPI.COMM_WORLD
    msh = dmesh.create_box(
        comm,
        [np.array([0.0, 0.0, 0.0]), np.array([BOX_L, BOX_L, BOX_L])],
        [n, n, n],
        cell_type=dmesh.CellType.tetrahedron,
    )
    problem = TimeHarmonicProblem(
        mesh=msh,
        frequency_hz=FREQUENCY_HZ,
        material=HomogeneousMaterial(sigma=SIGMA, epsilon_r=EPSILON_R, mu_r=MU_R),
        boundary_condition=mode,
    )
    solver = TimeHarmonicSolver(problem, degree=1)

    # The API's dirichlet_dof_count is rank-local and includes ghosts: count
    # owned dofs (index < size_local) and allreduce those.
    bcs, _, _ = solver.build_boundary_conditions()
    size_local = solver.function_space().dofmap.index_map.size_local
    owned_local = 0
    for bc in bcs:
        dofs = np.asarray(bc.dof_indices()[0])
        owned_local += int(np.count_nonzero(dofs < size_local))
    owned = comm.allreduce(owned_local, op=MPI.SUM)

    fields = solver.solve(current_density=_th3_current(field), project_source=False)

    x = ufl.SpatialCoordinate(msh)
    exact = _th3_field_ufl(field, x)
    diff = fields.e_complex - exact
    # assemble_scalar is rank-local: allreduce before the square root.
    err_sq = comm.allreduce(fem.assemble_scalar(fem.form(ufl.inner(diff, diff) * ufl.dx)), op=MPI.SUM)
    ref_sq = comm.allreduce(fem.assemble_scalar(fem.form(ufl.inner(exact, exact) * ufl.dx)), op=MPI.SUM)
    ncells = msh.topology.index_map(msh.topology.dim).size_global
    result = (float(np.sqrt(abs(err_sq)) / np.sqrt(abs(ref_sq))), int(owned), int(ncells))
    if return_fields:
        return (*result, msh, fields)
    return result


def _th3_exterior_edges(n: int) -> int:
    # Surface of the n^3 tet box: V = 6n^2 + 2, F = 12n^2 triangles, E = 18n^2.
    return 18 * n * n


@complex_only
@pytest.mark.parametrize(
    "mode, field, expected_constrained",
    [
        (TimeHarmonicBoundaryCondition.PEC_ZERO_TANGENTIAL_A, "pec", _th3_exterior_edges),
        (TimeHarmonicBoundaryCondition.NATURAL, "nat", lambda n: 0),
    ],
    ids=["pec", "natural"],
)
def test_th3_boundary_mode_solves_to_its_manufactured_field(mode, field, expected_constrained):
    comm = MPI.COMM_WORLD
    rel = {}
    for n in TH3_SIZES:
        rel[n], owned, ncells = _solve_mms(mode.value, field, n)
        if comm.rank == 0:
            print(
                f"\n[TH-3] mode={mode.value} field=E_{field} n={n}^3 ({ncells} cells): "
                f"rel L2 = {rel[n]:.6e}, owned constrained dofs = {owned} "
                f"(expected {expected_constrained(n)})"
            )
        # (0) exact count, allreduced owned dofs.
        assert owned == expected_constrained(n), (
            f"{mode.value} at {n}^3 constrained {owned} owned dofs, expected {expected_constrained(n)}"
        )
    coarse, fine = rel[TH3_SIZES[0]], rel[TH3_SIZES[1]]
    rate = float(np.log(coarse / fine) / np.log(TH3_SIZES[1] / TH3_SIZES[0]))
    if comm.rank == 0:
        print(f"[TH-3] mode={mode.value}: rate 12^3 -> 24^3 = {rate:.4f}")
    # (i) TH-1's bands.
    assert fine < TH3_REL_L2_BAND, f"{mode.value}: relative L2 at 24^3 {fine:.4e} >= {TH3_REL_L2_BAND}"
    assert fine < coarse, f"{mode.value}: refinement did not reduce the error ({coarse:.4e} -> {fine:.4e})"
    assert rate > TH3_RATE_FLOOR, f"{mode.value}: L2 rate {rate:.4f} <= {TH3_RATE_FLOOR}"


@complex_only
@pytest.mark.parametrize(
    "mode, field",
    [
        (TimeHarmonicBoundaryCondition.PEC_ZERO_TANGENTIAL_A, "nat"),
        (TimeHarmonicBoundaryCondition.NATURAL, "pec"),
    ],
    ids=["nat-source-under-pec", "pec-source-under-natural"],
)
def test_th3_swapped_boundary_mode_misses_the_field(mode, field):
    """Negative control: the other mode's source under this mode's boundary.

    Predicted miss O(1) (approximately 30-100 %); asserted only as > 0.10.
    """
    comm = MPI.COMM_WORLD
    n = TH3_SIZES[0]
    rel, owned, ncells = _solve_mms(mode.value, field, n)
    if comm.rank == 0:
        print(
            f"\n[TH-3 control] mode={mode.value} field=E_{field} n={n}^3 ({ncells} cells): "
            f"rel L2 = {rel:.6e}, owned constrained dofs = {owned}"
        )
    assert rel > TH3_REL_L2_BAND, f"swap control {mode.value}/E_{field} reads {rel:.4e}, not a miss"
