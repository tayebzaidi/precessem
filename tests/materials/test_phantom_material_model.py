"""Tests for gelled saline phantom material integration (chunk C2)."""

from __future__ import annotations

import numpy as np
import ufl
from dolfinx import fem
from mpi4py import MPI

from fem_em_solver.core import (
    HomogeneousMaterial,
    TimeHarmonicProblem,
    TimeHarmonicSolver,
    build_material_fields,
)
from fem_em_solver.io.mesh import MeshGenerator
from fem_em_solver.materials import GelledSalinePhantomMaterial
from tests.complex_mode import complex_only

# GEO-23 step 2c (2026-08-28): the global cell count step 2b only printed,
# now asserted.  Reference is step 1's -n 1 ladder on this generator, measured
# on the dolfinx 0.11 image; the 0.7.2-era sizing (resolution=0.03) no longer
# meshes at all, so there is no older count to compare against.  Step 2b
# reproduced 5464 exactly (0.00%) at -n 1 and -n 2.
N_CELLS_REF = 5464


def _global_bounds_for_tagged_cells(field: fem.Function, cell_tags, tag: int, comm: MPI.Intracomm):
    """Return global (min, max, sample_count) over DG0 dofs in cells with given tag."""
    q0 = field.function_space
    cells = cell_tags.indices[cell_tags.values == tag]

    local_vals = []
    for cell in cells:
        dofs = q0.dofmap.cell_dofs(int(cell))
        local_vals.extend(field.x.array[dofs])

    if local_vals:
        # sigma / epsilon_r are real material coefficients, but the complex build
        # stores every dof array as complex128. Drop the (identically zero)
        # imaginary part explicitly instead of letting float() do it and emit a
        # ComplexWarning; the max |Im| is asserted below on the reduced values.
        local_array = np.asarray(local_vals).real
        local_min = float(np.min(local_array))
        local_max = float(np.max(local_array))
        local_imag = float(np.max(np.abs(np.imag(np.asarray(local_vals)))))
        local_count = int(len(local_vals))
    else:
        local_min = float("inf")
        local_max = float("-inf")
        local_imag = 0.0
        local_count = 0

    global_min = comm.allreduce(local_min, op=MPI.MIN)
    global_max = comm.allreduce(local_max, op=MPI.MAX)
    global_imag = comm.allreduce(local_imag, op=MPI.MAX)
    global_count = comm.allreduce(local_count, op=MPI.SUM)
    # Collective: every rank asserts the same reduced value, so a violation fails
    # everywhere rather than deadlocking the surviving ranks at the next barrier.
    assert global_imag == 0.0, f"material coefficient has a non-zero imaginary part: {global_imag}"
    return global_min, global_max, global_count


def test_gelled_saline_preset_table_supports_low_mid_high_conductivity_variants():
    low = GelledSalinePhantomMaterial.from_preset("low", frequency_hz=127.74e6)
    mid = GelledSalinePhantomMaterial.from_preset("mid", frequency_hz=127.74e6)
    high = GelledSalinePhantomMaterial.from_preset("high", frequency_hz=127.74e6)

    assert GelledSalinePhantomMaterial.preset_names() == ("high", "low", "mid")
    assert low.sigma < mid.sigma < high.sigma

    for material in (low, mid, high):
        material.validate()
        assert material.epsilon_r > 0.0
        assert np.isfinite(material.conduction_plus_displacement)


def test_gelled_saline_frequency_adjustment_hook_produces_finite_terms():
    base = GelledSalinePhantomMaterial.from_preset("mid", frequency_hz=64.0e6)
    adjusted = base.with_frequency_adjustment(
        128.0e6,
        conductivity_exponent=0.15,
        permittivity_exponent=-0.05,
    )

    assert adjusted.frequency_hz == 128.0e6
    assert adjusted.sigma > 0.0
    assert adjusted.epsilon_r > 0.0
    assert np.isfinite(adjusted.displacement_term)
    assert np.isfinite(adjusted.conduction_plus_displacement)
    assert np.isfinite(adjusted.complex_admittivity.real)
    assert np.isfinite(adjusted.complex_admittivity.imag)


def test_gelled_saline_material_container_frequency_term_is_finite():
    material = GelledSalinePhantomMaterial(sigma=0.7, epsilon_r=78.0, frequency_hz=127.74e6)
    material.validate()

    assert material.omega > 0.0
    assert np.isfinite(material.displacement_term)
    assert material.conduction_plus_displacement > material.sigma


@complex_only
def test_phantom_material_assignment_and_time_harmonic_pipeline_wiring():
    """Assign phantom material to phantom-tagged cells and verify solver returns these fields."""
    comm = MPI.COMM_WORLD

    # GEO-23 step 2b (2026-08-28): 0.03 does not mesh on dolfinx 0.11 — gmsh
    # throws "Invalid boundary mesh (overlapping facets)".  Step 1's -n 1
    # ladder on this generator is monotone and its coarsest meshing rung is
    # 0.024 (5 464 cells); the finer rungs are 0.0192 (9 330) / 0.01536
    # (16 177) / 0.012288 (28 485).  Moved to the coarsest meshing rung.
    mesh, cell_tags, facet_tags = MeshGenerator.coil_phantom_domain(
        coil_major_radius=0.07,
        coil_minor_radius=0.010,
        coil_separation=0.08,
        phantom_radius=0.03,
        phantom_height=0.08,
        air_padding=0.04,
        resolution=0.024,
        comm=comm,
    )
    # size_global is already a global quantity on every rank — do not sum
    # size_local on top of it.
    n_global = mesh.topology.index_map(mesh.topology.dim).size_global
    if comm.rank == 0:
        print(f"GEO-23 step 2b: coil_phantom_domain(resolution=0.024) -> {n_global} cells")
    assert abs(n_global / N_CELLS_REF - 1) <= 0.01, (
        f"coil_phantom_domain(resolution=0.024) meshed {n_global} cells, "
        f"outside +/-1% of the GEO-23 step 1 ladder reference {N_CELLS_REF}"
    )

    background = HomogeneousMaterial(sigma=0.0, epsilon_r=1.0, mu_r=1.0)
    phantom = GelledSalinePhantomMaterial.from_preset("mid", frequency_hz=127.74e6)

    sigma_field, epsilon_r_field = build_material_fields(
        mesh,
        background,
        cell_tags=cell_tags,
        phantom_material=phantom,
        phantom_tag=3,
    )

    sigma_min, sigma_max, sigma_count = _global_bounds_for_tagged_cells(
        sigma_field, cell_tags, tag=3, comm=comm
    )
    eps_min, eps_max, eps_count = _global_bounds_for_tagged_cells(
        epsilon_r_field, cell_tags, tag=3, comm=comm
    )

    assert sigma_count > 0 and eps_count > 0, "phantom tag must own at least one DG0 dof"
    assert np.isclose(sigma_min, phantom.sigma) and np.isclose(sigma_max, phantom.sigma)
    assert np.isclose(eps_min, phantom.epsilon_r) and np.isclose(eps_max, phantom.epsilon_r)

    # Verify the wiring is active in the solve pipeline: the returned fields
    # carry the assigned material maps.  `MAT-1` (2026-09-29): the
    # `DummyMagnetostaticSolver` monkeypatch that used to wrap this solve was
    # dead — `TimeHarmonicSolver.solve` never builds a MagnetostaticSolver — so
    # it is removed; this is the same real curl-curl solve it always was, and
    # the asserts below are unchanged.  The quantitative gate on the preset
    # path is `test_gelled_saline_presets_match_lossy_sphere_closed_form`.
    problem = TimeHarmonicProblem(
        mesh=mesh,
        frequency_hz=127.74e6,
        material=background,
        cell_tags=cell_tags,
        facet_tags=facet_tags,
        phantom_material=phantom,
        phantom_tag=3,
    )
    solver = TimeHarmonicSolver(problem, degree=1)

    def current_density(x):
        return ufl.as_vector([0.0, 0.0, 1.0])

    fields = solver.solve(current_density=current_density, subdomain_ids=[1, 2], gauge_penalty=1e-3)

    assert fields.sigma_field is not None
    assert fields.epsilon_r_field is not None

    fs_min, fs_max, _ = _global_bounds_for_tagged_cells(
        fields.sigma_field, cell_tags, tag=3, comm=comm
    )
    fe_min, fe_max, _ = _global_bounds_for_tagged_cells(
        fields.epsilon_r_field, cell_tags, tag=3, comm=comm
    )

    assert np.isclose(fs_min, phantom.sigma) and np.isclose(fs_max, phantom.sigma)
    assert np.isclose(fe_min, phantom.epsilon_r) and np.isclose(fe_max, phantom.epsilon_r)


# ---------------------------------------------------------------------------
# `MAT-1` (2026-09-29): the gelled-saline presets through the real solve,
# against `MAT-4` step 1's lossy-sphere closed form.
#
# Fixture, constants and probe points are `MAT-4` step 1's, imported (names
# only).  `_solve_lossy_sphere` / `_exact_exterior_numpy` there hard-code
# eps_r = 78, so the wall data below is `test_lossy_sphere_sar.py:98–124`
# re-parameterised on the case's own (sigma, eps_r).
#
# Bands: `MAT-4`'s mean-SAR and Im/Re bands are the literal 0.10 at
# `test_lossy_sphere_sar.py:321` and `:347`; there is no named constant to
# import, so the same value is bound here once and cited, never re-chosen.
# ---------------------------------------------------------------------------

from dolfinx.fem.petsc import assemble_matrix  # noqa: E402
from petsc4py import PETSc  # noqa: E402

from fem_em_solver.core import build_mu_r_field  # noqa: E402
from fem_em_solver.post.evaluation import evaluate_vector_field_parallel  # noqa: E402
from fem_em_solver.post.sar import mean_sar, uniform_sphere_sar_closed_form  # noqa: E402
from fem_em_solver.utils.constants import EPSILON_0  # noqa: E402
from tests.validation.test_lossy_sphere_sar import (  # noqa: E402
    BOX_HALF_WIDTH,
    E0,
    EPSILON_R_SPHERE,
    FREQUENCY_HZ,
    K0,
    OMEGA,
    RHO_KG_M3,
    SIGMA_HIGH,
    SPHERE_RADIUS,
    SPHERE_TAG,
    _interior_probe_points,
)

MAT4_BAND = 0.10  # = test_lossy_sphere_sar.py:321 (mean SAR) and :347 (Im/Re)

# `MAT-1` premise, re-anchored 2026-09-29 18:00 review (§9 item 83; known-issues
# `MAT-1` entry).  Attempt 1 gated the *solve outputs* of the two paths to
# 1e-12 and went red on an identical re-run: 9.047e-13 / 5.962e-13 then
# 2.544e-11 / 1.200e-11 (SAR / <E_z>, `20260929T183740Z_MAT-1-W1.log:485`,
# `20260929T184225Z_MAT-1-W1-final.log:485`), because the `-n 2` complex MUMPS
# solve on this fixture does not reproduce itself to 1e-12 (map-vs-map
# 9.797e-13 … 4.45e-11 threaded, 9.047e-13 … 4.24e-11 at threads = 1,
# `20260929T183957Z_MAT-1-probe-noise.log:105–106`,
# `20260929T184102Z_MAT-1-probe-omp1.log:106–107`).  The 1e-12 bound is not
# widened; it is withdrawn as a gate on a solve output.  The path's actual
# claim — the phantom branch *assembles the same material* as material_map —
# is asserted upstream of the solve, exactly: bitwise DG0 sigma / eps_r / mu_r
# on every owned dof, and the curl-curl operator built from them.  The solve
# difference is printed beside the measured floor below; above
# SOLVE_RECORD_FLAG it is flagged in the journal and gates nothing.
SOLVE_FLOOR_RECORD = (5.962e-13, 4.450e-11)  # measured -n 2 floor, logs cited above
SOLVE_RECORD_FLAG = 1e-9
# Control (c2)'s relative sigma perturbation on the operator.  A 1-ULP sigma
# is not used there: it can round away inside the assembled entries.
SIGMA_REL_PERTURBATION = 1e-9


def _dg0_owned_diff(f_a: fem.Function, f_b: fem.Function, comm: MPI.Intracomm):
    """Compare two DG0 fields on the **owned** dofs only.

    ``.x.array`` is rank-local and its ghost tail is a copy of another rank's
    owned value, so the owned slice ``[:size_local]`` is the whole field
    exactly once.  Both fields must share the per-rank layout (asserted on the
    reduced verdict, so a mismatch fails on every rank).  Arrays are compared
    whole as complex128 — never ``.real`` first.

    Returns (max |delta| over ranks, differing count summed, owned count summed).
    """
    im_a = f_a.function_space.dofmap.index_map
    im_b = f_b.function_space.dofmap.index_map
    bs_a = f_a.function_space.dofmap.index_map_bs
    bs_b = f_b.function_space.dofmap.index_map_bs
    local_same = (
        im_a.size_local == im_b.size_local
        and im_a.num_ghosts == im_b.num_ghosts
        and bs_a == bs_b == 1
    )
    same_layout = bool(comm.allreduce(bool(local_same), op=MPI.LAND))
    assert same_layout, "DG0 fields do not share the per-rank owned/ghost layout"
    n_owned = im_a.size_local
    a = np.asarray(f_a.x.array[:n_owned], dtype=np.complex128)
    b = np.asarray(f_b.x.array[:n_owned], dtype=np.complex128)
    delta = np.abs(a - b)
    local_max = float(np.max(delta)) if n_owned else 0.0
    local_diff = int(np.count_nonzero(a != b))
    return (
        float(comm.allreduce(local_max, op=MPI.MAX)),
        int(comm.allreduce(local_diff, op=MPI.SUM)),
        int(comm.allreduce(int(n_owned), op=MPI.SUM)),
    )


def _owned_count_equal(f: fem.Function, value: float, comm: MPI.Intracomm) -> int:
    """Owned DG0 dofs whose complex value is exactly ``value + 0j`` (reduced)."""
    n_owned = f.function_space.dofmap.index_map.size_local
    arr = np.asarray(f.x.array[:n_owned], dtype=np.complex128)
    return int(comm.allreduce(int(np.count_nonzero(arr == complex(value, 0.0))), op=MPI.SUM))


def _volume_operator(solver, sigma_f, epsilon_r_f, mu_r_f):
    """The curl-curl operator ``solve`` factors, assembled with no bcs."""
    mat = assemble_matrix(fem.form(solver.bilinear_form(sigma_f, epsilon_r_f, mu_r_field=mu_r_f)))
    mat.assemble()
    return mat


# `_same_pattern` / `_difference`: copied verbatim from
# `tests/validation/test_port19_matrix_premise.py:112–137` (`PORT-19` step 1).
# Copied, not imported: importing that module pulls in
# `test_port_birdcage_four_port` and its ten-module birdcage chain.
# `_difference`'s census allreduce is guarded by a rank-local shape test, so
# it is only called here after `_same_pattern` has asserted equal CSR shapes
# on every rank.
def _same_pattern(mat_a, mat_b, comm):
    """Rank-local CSR structure compared per rank; verdict reduced with AND."""
    ia_a, ja_a, _ = mat_a.getValuesCSR()
    ia_b, ja_b, _ = mat_b.getValuesCSR()
    local = bool(np.array_equal(ia_a, ia_b) and np.array_equal(ja_a, ja_b))
    return bool(comm.allreduce(local, op=MPI.LAND))


def _difference(mat_k, mat_1, comm):
    """``||A_k - A_1||_inf`` (global) plus a rank-reduced differing-entry census."""
    diff = mat_k.copy()
    diff.axpy(-1.0, mat_1, structure=PETSc.Mat.Structure.DIFFERENT_NONZERO_PATTERN)
    norm_inf = float(diff.norm(PETSc.NormType.NORM_INFINITY))
    _, _, vals_k = mat_k.getValuesCSR()
    _, _, vals_1 = mat_1.getValuesCSR()
    if vals_k.shape == vals_1.shape:
        changed = int(np.count_nonzero(vals_k != vals_1))
        max_abs = float(np.max(np.abs(vals_k - vals_1))) if vals_k.size else 0.0
    else:
        changed, max_abs = -1, float("nan")
    diff.destroy()
    return (
        norm_inf,
        int(comm.allreduce(changed, op=MPI.SUM)) if changed >= 0 else -1,
        float(comm.allreduce(max_abs, op=MPI.MAX)),
    )


def _epsilon_c_case(sigma: float, epsilon_r: float) -> complex:
    return complex(epsilon_r, -sigma / (OMEGA * EPSILON_0))


def _exterior_case(sigma: float, epsilon_r: float):
    """`test_lossy_sphere_sar.py:98–124` on the case's own (sigma, eps_r)."""
    epsilon_c = _epsilon_c_case(sigma, epsilon_r)
    beta = (epsilon_c - 1.0) / (epsilon_c + 2.0)
    e_in = 3.0 * E0 / (epsilon_c + 2.0)
    R = SPHERE_RADIUS

    def field(x: np.ndarray) -> np.ndarray:
        r = np.sqrt(x[0] ** 2 + x[1] ** 2 + x[2] ** 2)
        r_safe = np.where(r > 0.0, r, 1.0)
        z = x[2]
        common = beta * E0 * R**3 / r_safe**5
        ex_out = common * 3.0 * z * x[0]
        ey_out = common * 3.0 * z * x[1]
        ez_out = E0 + common * (3.0 * z * z - r_safe**2)
        inside = r < R
        ex = np.where(inside, 0.0, ex_out)
        ey = np.where(inside, 0.0, ey_out)
        ez = np.where(inside, e_in, ez_out)
        return np.array([ex, ey, ez], dtype=PETSc.ScalarType)

    return field


def _solve_case(msh, cell_tags, sigma, epsilon_r, *, phantom=None):
    """One lossy-sphere solve.  The sphere material goes in by ``material_map``
    unless ``phantom`` is given; then it goes in as ``phantom_material`` on
    ``phantom_tag = SPHERE_TAG`` (the default 3 raises on this mesh).
    """
    comm = MPI.COMM_WORLD
    kwargs = {}
    if phantom is None:
        kwargs["material_map"] = {SPHERE_TAG: HomogeneousMaterial(sigma=sigma, epsilon_r=epsilon_r)}
    else:
        kwargs["phantom_material"] = phantom
        kwargs["phantom_tag"] = SPHERE_TAG
    problem = TimeHarmonicProblem(
        mesh=msh,
        frequency_hz=FREQUENCY_HZ,
        material=HomogeneousMaterial(sigma=0.0, epsilon_r=1.0),
        cell_tags=cell_tags,
        boundary_condition="pec_zero_tangential_a",
        dirichlet_e_field=_exterior_case(sigma, epsilon_r),
        **kwargs,
    )
    solver = TimeHarmonicSolver(problem, degree=1)
    fields = solver.solve()
    sar = mean_sar(
        fields.e_complex,
        sigma=fields.sigma_field,
        rho=RHO_KG_M3,
        cell_tags=cell_tags,
        subdomain_ids=SPHERE_TAG,
        comm=comm,
    )
    points = _interior_probe_points()
    real_values, valid_real = evaluate_vector_field_parallel(fields.e_real, points, comm)
    imag_values, valid_imag = evaluate_vector_field_parallel(fields.e_imag, points, comm)
    valid = valid_real & valid_imag
    if not np.all(valid):
        raise RuntimeError(f"{int((~valid).sum())} interior probe points were not evaluated")
    ez = np.real(real_values[:, 2]) + 1j * np.real(imag_values[:, 2])
    # `fields.sigma_field` / `epsilon_r_field` / `mu_r_field` are the objects
    # `bilinear_form` consumed (`core/time_harmonic.py:636–642`, `:749–751`).
    return {
        "mean_sar": float(sar["mean_sar_w_per_kg"]),
        "ez_mean": complex(np.mean(ez)),
        "fields": fields,
        "solver": solver,
    }


def _closed_case(sigma: float, epsilon_r: float) -> dict:
    c = uniform_sphere_sar_closed_form(
        e0=E0,
        epsilon_r=epsilon_r,
        sigma=sigma,
        omega=OMEGA,
        rho=RHO_KG_M3,
        epsilon_0=EPSILON_0,
    )
    e_in = 3.0 * E0 / (_epsilon_c_case(sigma, epsilon_r) + 2.0)
    c["im_re"] = e_in.imag / e_in.real
    c["k_in_r"] = float(K0 * np.sqrt(c["epsilon_c_magnitude"]) * SPHERE_RADIUS)
    return c


# PROJECT_PLAN §9 item 81 predictions, printed beside the measurement:
# (closed-form SAR [W/kg], closed Im/Re, |k_in|R, sigma-blind inflation ceiling)
_PREDICTED = {
    "low": (8.92e-8, 1.541, 0.164, 3.38),
    "mid": (6.65e-8, 2.633, 0.201, 7.93),
    "high": (5.08e-8, 3.781, 0.234, 15.3),
}


def build_preset_sphere_mesh(comm):
    """The `MAT-1` sphere-in-box mesh, lifted verbatim from the gate body so
    `examples/materials/03_gelled_saline_presets.py` (`EX-65`) builds the
    same mesh by import.  Returns the generator's (mesh, cell_tags, facet_tags)."""
    return MeshGenerator.sphere_in_box_domain(
        sphere_radius=SPHERE_RADIUS,
        box_half_width=BOX_HALF_WIDTH,
        resolution_sphere=SPHERE_RADIUS / 10.0,
        resolution_far=SPHERE_RADIUS / 5.0,
        comm=comm,
    )


@complex_only
def test_gelled_saline_presets_match_lossy_sphere_closed_form():
    """Each preset, passed as ``phantom_material``, gives mean SAR and interior
    Im/Re E_z within `MAT-4`'s 10 % of ``uniform_sphere_sar_closed_form``.

    Premise first (§9 item 83): at sigma = SIGMA_HIGH, eps_r = 78 the phantom
    path assembles bitwise the same DG0 sigma / eps_r / mu_r and the same
    curl-curl operator as the ``material_map`` path; a 1-ULP sigma control on
    the arrays and a 1e-9 relative sigma control on the operator read > 0.
    The solve-level difference is printed beside the measured floor only."""
    comm = MPI.COMM_WORLD
    msh, cell_tags, _ = build_preset_sphere_mesh(comm)
    n_global = msh.topology.index_map(msh.topology.dim).size_global

    # --- Premise: the phantom branch assembles the same field as material_map.
    via_map = _solve_case(msh, cell_tags, SIGMA_HIGH, EPSILON_R_SPHERE)
    via_phantom = _solve_case(
        msh,
        cell_tags,
        SIGMA_HIGH,
        EPSILON_R_SPHERE,
        phantom=GelledSalinePhantomMaterial(SIGMA_HIGH, EPSILON_R_SPHERE, FREQUENCY_HZ),
    )
    premise_sar = abs(via_phantom["mean_sar"] - via_map["mean_sar"]) / abs(via_map["mean_sar"])
    premise_ez = abs(via_phantom["ez_mean"] - via_map["ez_mean"]) / abs(via_map["ez_mean"])
    c_premise = _closed_case(SIGMA_HIGH, EPSILON_R_SPHERE)
    premise_map_err = (
        abs(via_map["mean_sar"] - c_premise["sar_w_per_kg"]) / c_premise["sar_w_per_kg"]
    )

    f_map = via_map["fields"]
    f_ph = via_phantom["fields"]

    # (P1) the three DG0 arrays the two solves' operators consumed.
    p1 = {
        "sigma": _dg0_owned_diff(f_ph.sigma_field, f_map.sigma_field, comm),
        "eps_r": _dg0_owned_diff(f_ph.epsilon_r_field, f_map.epsilon_r_field, comm),
        "mu_r": _dg0_owned_diff(f_ph.mu_r_field, f_map.mu_r_field, comm),
    }
    # (P1') non-vacuity: the sphere dofs carry sigma = SIGMA_HIGH on both
    # paths, and the two paths really built two distinct field objects.
    n_sigma_high_map = _owned_count_equal(f_map.sigma_field, SIGMA_HIGH, comm)
    n_sigma_high_ph = _owned_count_equal(f_ph.sigma_field, SIGMA_HIGH, comm)
    distinct_objects = f_map.sigma_field is not f_ph.sigma_field

    # (P2) the curl-curl operator each path's fields assemble, each through
    # its own solver (as `solve` built it).  Pattern before norm.
    a_map = _volume_operator(
        via_map["solver"], f_map.sigma_field, f_map.epsilon_r_field, f_map.mu_r_field
    )
    a_ph = _volume_operator(
        via_phantom["solver"], f_ph.sigma_field, f_ph.epsilon_r_field, f_ph.mu_r_field
    )
    a_map_norm = float(a_map.norm(PETSc.NormType.NORM_INFINITY))
    p2_pattern = _same_pattern(a_ph, a_map, comm)

    if comm.rank == 0:
        print(
            f"\n[MAT-1] lossy sphere R = {SPHERE_RADIUS} m, h = R/10 in / R/5 out, "
            f"{n_global} cells, f = {FREQUENCY_HZ / 1e6:.3f} MHz, rho = {RHO_KG_M3}"
        )
        for key, (dmax, ndiff, nown) in p1.items():
            print(
                f"  (P1) {key:5s}: max |delta| over owned dofs = {dmax!r}, differing = "
                f"{ndiff}, owned = {nown} (predicted 0.0 / 0 / 74020; global cells {n_global})"
            )
        print(
            f"  (P1') owned dofs with sigma == {SIGMA_HIGH}: map {n_sigma_high_map}, "
            f"phantom {n_sigma_high_ph}; distinct field objects = {distinct_objects}"
        )
        print(
            f"  (P2) same CSR pattern (per rank, LAND) = {p2_pattern}; ||A_map||_inf = {a_map_norm:.6e}"
        )

    for key, (dmax, ndiff, nown) in p1.items():
        assert (
            dmax == 0.0 and ndiff == 0
        ), f"(P1) phantom path writes a different {key}: max |delta| {dmax!r}, {ndiff} dofs differ"
        assert nown == n_global, f"(P1) {key}: owned dofs sum to {nown}, not the {n_global} cells"
    assert (
        n_sigma_high_map > 0 and n_sigma_high_ph == n_sigma_high_map
    ), f"(P1') sigma == {SIGMA_HIGH} on {n_sigma_high_map} (map) / {n_sigma_high_ph} (phantom) owned dofs"
    assert distinct_objects, "(P1') both paths returned the same sigma Function object"
    assert p2_pattern, "(P2) operators assembled from the two paths differ in nonzero pattern"

    p2_norm, p2_changed, p2_max = _difference(a_ph, a_map, comm)
    if comm.rank == 0:
        print(
            f"  (P2) ||A_phantom - A_map||_inf = {p2_norm!r}, differing entries = {p2_changed}, "
            f"max |entry delta| = {p2_max!r} (predicted 0.0)"
        )
    assert a_map_norm > 0.0, "(P2) A_map is the zero operator"
    assert p2_norm == 0.0, f"(P2) ||A_phantom - A_map||_inf = {p2_norm!r}, predicted 0.0"

    # --- Negative controls (no solve).
    # (c1) a 1-ULP sigma through the phantom path: the array compare must see it.
    sigma_ulp = float(np.nextafter(SIGMA_HIGH, 1.0))
    background = HomogeneousMaterial(sigma=0.0, epsilon_r=1.0)
    sigma_c1, eps_c1 = build_material_fields(
        msh,
        background,
        cell_tags=cell_tags,
        phantom_material=GelledSalinePhantomMaterial(sigma_ulp, EPSILON_R_SPHERE, FREQUENCY_HZ),
        phantom_tag=SPHERE_TAG,
    )
    mu_c1 = build_mu_r_field(msh, background, cell_tags=cell_tags)
    c1_sigma = _dg0_owned_diff(sigma_c1, f_map.sigma_field, comm)
    c1_eps = _dg0_owned_diff(eps_c1, f_map.epsilon_r_field, comm)
    c1_mu = _dg0_owned_diff(mu_c1, f_map.mu_r_field, comm)

    # (c2) sigma * (1 + 1e-9) on the operator: the matrix compare must see it.
    sigma_c2 = fem.Function(f_map.sigma_field.function_space, name="sigma_c2")
    sigma_c2.x.array[:] = f_map.sigma_field.x.array * (1.0 + SIGMA_REL_PERTURBATION)
    sigma_c2.x.scatter_forward()
    a_c2 = _volume_operator(via_map["solver"], sigma_c2, f_map.epsilon_r_field, f_map.mu_r_field)
    c2_pattern = _same_pattern(a_c2, a_map, comm)

    if comm.rank == 0:
        print(
            f"  (c1) sigma = nextafter({SIGMA_HIGH}, 1) = {sigma_ulp!r}: max |d sigma| = "
            f"{c1_sigma[0]!r} (predicted 2^-53 = {2.0**-53!r}), differing = {c1_sigma[1]} "
            f"(predicted {n_sigma_high_map}); d eps_r = {c1_eps[0]!r} ({c1_eps[1]}), "
            f"d mu_r = {c1_mu[0]!r} ({c1_mu[1]}) (predicted 0)"
        )
        print(f"  (c2) sigma * (1 + {SIGMA_REL_PERTURBATION:g}): same CSR pattern = {c2_pattern}")
    assert (
        c1_sigma[0] > 0.0 and c1_sigma[1] > 0
    ), f"(c1) a 1-ULP sigma reads {c1_sigma[0]!r}: the array compare is blind"
    assert c2_pattern, "(c2) perturbed-sigma operator changed its nonzero pattern"
    c2_norm, c2_changed, c2_max = _difference(a_c2, a_map, comm)
    if comm.rank == 0:
        print(
            f"  (c2) ||A'' - A_map||_inf = {c2_norm!r}, relative {c2_norm / a_map_norm:.3e} "
            f"(predicted <~ 1e-13), differing entries = {c2_changed}, max |entry delta| = {c2_max!r}"
        )
    assert c2_norm > 0.0, "(c2) a 1e-9 relative sigma reads 0: the operator compare is blind"
    for mat in (a_map, a_ph, a_c2):
        mat.destroy()

    # (P3) and the solve-level record.
    if comm.rank == 0:
        flag = (
            f"  ** above {SOLVE_RECORD_FLAG:g}: flag in the journal **"
            if max(premise_sar, premise_ez) > SOLVE_RECORD_FLAG
            else ""
        )
        print(
            f"  solve-level record (printed, gates nothing): material_map SAR = "
            f"{via_map['mean_sar']:.9e}, phantom SAR = {via_phantom['mean_sar']:.9e}, "
            f"rel diff SAR = {premise_sar:.3e}, rel diff <E_z> = {premise_ez:.3e} "
            f"(measured -n 2 floor {SOLVE_FLOOR_RECORD[0]:.3e} … {SOLVE_FLOOR_RECORD[1]:.3e}){flag}"
        )
        print(
            f"  (P3) material_map vs closed {c_premise['sar_w_per_kg']:.6e}: {premise_map_err:.3%} "
            f"(predicted 3.536%, band {MAT4_BAND:.0%})"
        )
    assert (
        premise_map_err < MAT4_BAND
    ), f"material_map premise run is {premise_map_err:.2%} off the closed form"

    # --- Presets through the phantom path.
    names = ("low", "mid", "high")
    rows = {}
    for name in names:
        mat = GelledSalinePhantomMaterial.from_preset(name, frequency_hz=FREQUENCY_HZ)
        run = _solve_case(msh, cell_tags, mat.sigma, mat.epsilon_r, phantom=mat)
        c = _closed_case(mat.sigma, mat.epsilon_r)
        im_re = run["ez_mean"].imag / run["ez_mean"].real
        t = c["loss_tangent_numerator"]
        rows[name] = {
            "mat": mat,
            "closed": c,
            "run": run,
            "sar_err": abs(run["mean_sar"] - c["sar_w_per_kg"]) / c["sar_w_per_kg"],
            "im_re": im_re,
            "im_re_err": abs(im_re - c["im_re"]) / abs(c["im_re"]),
            # sigma-blind inflation ceiling ((eps_r+2)^2 + t^2) / (eps_r+2)^2
            "inflation": ((mat.epsilon_r + 2.0) ** 2 + t**2) / (mat.epsilon_r + 2.0) ** 2,
            # what a sigma-blind solve returns: E_in = 3 E0 / (eps_r + 2)
            "sar_blind": mat.sigma * (3.0 * E0 / (mat.epsilon_r + 2.0)) ** 2 / (2.0 * RHO_KG_M3),
        }

    ratio_fem = rows["high"]["run"]["mean_sar"] / rows["low"]["run"]["mean_sar"]
    ratio_closed = rows["high"]["closed"]["sar_w_per_kg"] / rows["low"]["closed"]["sar_w_per_kg"]
    ratio_blind = rows["high"]["sar_blind"] / rows["low"]["sar_blind"]
    separation = ratio_blind / ratio_fem
    separation_ceiling = ratio_blind / ratio_closed

    if comm.rank == 0:
        for name in names:
            r = rows[name]
            p = _PREDICTED[name]
            print(
                f"  preset {name:4s} sigma = {r['mat'].sigma:.2f}, eps_r = {r['mat'].epsilon_r:.0f}: "
                f"|k_in|R = {r['closed']['k_in_r']:.4f} (pred {p[2]}), "
                f"closed SAR = {r['closed']['sar_w_per_kg']:.4e} (pred {p[0]:.3g}), "
                f"FEM SAR = {r['run']['mean_sar']:.4e} err {r['sar_err']:.3%}; "
                f"Im/Re closed {r['closed']['im_re']:.4f} (pred {p[1]}), FEM {r['im_re']:.4f} "
                f"err {r['im_re_err']:.3%}; sigma-blind inflation ceiling "
                f"{r['inflation']:.2f}x (pred {p[3]})"
            )
        print(
            f"  sigma-blind control (printed, not asserted, rule (e)): FEM SAR_high/SAR_low = "
            f"{ratio_fem:.4f}, closed {ratio_closed:.4f}, sigma-blind {ratio_blind:.4f} "
            f"=> separation {separation:.3f} (ceiling {separation_ceiling:.3f}; "
            f"predicted blind 2.5788, closed 0.5692, separation 4.531)"
        )

    for name in names:
        r = rows[name]
        assert r["sar_err"] < MAT4_BAND, (
            f"preset {name}: mean SAR {r['run']['mean_sar']:.4e} vs closed "
            f"{r['closed']['sar_w_per_kg']:.4e} ({r['sar_err']:.2%}, |k_in|R = "
            f"{r['closed']['k_in_r']:.3f})"
        )
        assert r["im_re_err"] < MAT4_BAND, (
            f"preset {name}: Im/Re E_z {r['im_re']:.4f} vs closed "
            f"{r['closed']['im_re']:.4f} ({r['im_re_err']:.2%}, |k_in|R = "
            f"{r['closed']['k_in_r']:.3f})"
        )
