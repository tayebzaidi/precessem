"""`WF-15` step 1 (§9 item 240, chain U3): `Study` v0 on `MAT-4` step 1's sphere.

A ``Study`` composes a geometry generator, a `MAT-1` material, an imposed
field, a solve (degree 2 by default) and a SAR result.  Asserted (brief:
docs/planning/chunks/WF-15.md#§9-item-240):

(F)   the imposed field is `MAT-4`'s wall data, bitwise;
(Deg) the default degree is 2, and the solved space is degree 2;
(M)   the geometry generator equals the direct generator (global cell count);
(U)   Study vs direct call upstream of the solve, bitwise: DG0 sigma / eps_r /
      mu_r on owned dofs and the curl-curl operator (pattern + ||dA||_inf == 0);
(P)   SAR post-processing bitwise against ``mean_sar`` called directly;
(D1)  the direct copy reproduces `MAT-1`'s map record within ``L0_BAND``;
(I)   Study vs direct end to end within ``L0_BAND`` (SAR and <E_z>);
(R)   `MAT-4`'s degree-1 figure, ``L0_BAND`` and the ``.6e`` string;
(G)   inside `MAT-4`'s band (mean SAR, Im/Re E_z; spread < 0.02) at both degrees;
(Pre) the "mid" preset by name reproduces `MAT-1`'s ``.4e`` record.

Negative controls: (c1) 1-ULP sigma on the arrays and (c2) sigma*(1+1e-9) on
the operator at degree 1, (c3) wrong material through (I)'s comparator —
asserted; (c2) at degree 2 and (c4) degree mismatch — printed (predicted).

No solve output is compared bitwise: two `-n 2` MUMPS solves are not
(known-issues `MAT-1` entry, ``SOLVE_FLOOR_RECORD``).
"""

from __future__ import annotations

import datetime
import json
import os
import resource
import time
from types import SimpleNamespace

import numpy as np
import pytest
from dolfinx import fem
from mpi4py import MPI
from petsc4py import PETSc

from fem_em_solver import (
    GelledSalinePhantomMaterial,
    QuasiStaticSphereField,
    SphereInBoxGeometry,
    Study,
)
from fem_em_solver.core import (
    HomogeneousMaterial,
    TimeHarmonicProblem,
    TimeHarmonicSolver,
    build_material_fields,
    build_mu_r_field,
)
from fem_em_solver.post.evaluation import evaluate_vector_field_parallel
from fem_em_solver.post.sar import mean_sar, uniform_sphere_sar_closed_form
from fem_em_solver.utils.constants import EPSILON_0

from tests.complex_mode import complex_only
from tests.materials.test_phantom_material_model import (
    MAT4_BAND,
    SOLVE_FLOOR_RECORD,
    _difference,
    _dg0_owned_diff,
    _same_pattern,
    _volume_operator,
    build_preset_sphere_mesh,
)
from tests.validation.test_lossy_sphere_sar import (
    BOX_HALF_WIDTH,
    E0,
    EPSILON_R_SPHERE,
    FREQUENCY_HZ,
    OMEGA,
    RHO_KG_M3,
    SIGMA_HIGH,
    SPHERE_RADIUS,
    SPHERE_TAG,
    _exact_exterior_numpy,
    _interior_probe_points,
)
from tests.validation.test_mat5_em_coupled_sphere import L0_BAND

# `MAT-1`'s map-path record at sigma = 0.57 on this fixture, degree 1
# (`20260930T003229Z_MAT-1-reanchor-W1.log:494`) = `MAT-4` step 1's 0.11-image
# figure 8.291558e-08 (`20260921T111802Z_OPS-55-callers-complex.log.gz:1137`).
MAT4_DEG1_RECORD = 8.291558128e-08
MAT4_DEG1_PRINT = "8.291558e-08"  # `MAT-4`'s own print format, test_lossy_sphere_sar.py:303
# `MAT-1`'s "mid" preset record, `.4e` (`20260930T003229Z_MAT-1-reanchor-W1.log:497`).
MAT1_MID_PRINT = "6.8832e-08"
# Interior spread bound, by value with citation: test_lossy_sphere_sar.py:331.
SPREAD_BAND = 0.02
SIGMA_REL_PERTURBATION = 1e-9  # = MAT-1's (c2), test_phantom_material_model.py:243
# §9 item 251: the degree-2 Study leg (W2a) hands its solve outputs to the
# direct leg (W2b) across processes; /logs/ is gitignored (.gitignore:82).
DEG2_HANDOFF = "/workspace/logs/wf15/deg2_study.json"

GEOMETRY = SphereInBoxGeometry(
    sphere_radius=SPHERE_RADIUS,
    box_half_width=BOX_HALF_WIDTH,
    resolution_sphere=SPHERE_RADIUS / 10.0,
    resolution_far=SPHERE_RADIUS / 5.0,
    sphere_tag=SPHERE_TAG,
)


def _material_057():
    return GelledSalinePhantomMaterial(SIGMA_HIGH, EPSILON_R_SPHERE, FREQUENCY_HZ)


def _study(material=None, **kwargs):
    return Study(
        GEOMETRY,
        _material_057() if material is None else material,
        QuasiStaticSphereField(e0=E0),
        frequency_hz=FREQUENCY_HZ,
        density_kg_m3=RHO_KG_M3,
        **kwargs,
    )


def _ez_probe(fields, comm):
    points = _interior_probe_points()
    real_values, valid_real = evaluate_vector_field_parallel(fields.e_real, points, comm)
    imag_values, valid_imag = evaluate_vector_field_parallel(fields.e_imag, points, comm)
    valid = valid_real & valid_imag
    if not np.all(valid):
        raise RuntimeError(f"{int((~valid).sum())} interior probe points were not evaluated")
    return np.real(real_values[:, 2]) + 1j * np.real(imag_values[:, 2])


def _exterior_case(sigma, epsilon_r):
    """`test_phantom_material_model.py:336–357` (`MAT-1`'s copy of `MAT-4`'s wall data)."""
    epsilon_c = complex(epsilon_r, -sigma / (OMEGA * EPSILON_0))
    beta = (epsilon_c - 1.0) / (epsilon_c + 2.0)
    e_in = 3.0 * E0 / (epsilon_c + 2.0)
    R = SPHERE_RADIUS

    def field(x):
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


def _direct_call(msh, cell_tags, sigma, epsilon_r, degree, solve=True):
    """`MAT-1`'s ``_solve_case`` material_map branch
    (`test_phantom_material_model.py:360–407`) with a ``degree`` parameter.

    ``solve=False`` (§9 item 251, W2a) builds the same problem and solver but
    stops upstream of the factorisation: it returns the DG0 material fields
    built by the same calls ``TimeHarmonicSolver.solve`` makes
    (`core/time_harmonic.py:476–491`, ``build_material_fields`` then
    ``build_mu_r_field`` with the problem's own arguments) and the solver,
    which is all ``_volume_operator`` needs."""
    comm = MPI.COMM_WORLD
    problem = TimeHarmonicProblem(
        mesh=msh,
        frequency_hz=FREQUENCY_HZ,
        material=HomogeneousMaterial(sigma=0.0, epsilon_r=1.0),
        cell_tags=cell_tags,
        boundary_condition="pec_zero_tangential_a",
        dirichlet_e_field=_exterior_case(sigma, epsilon_r),
        material_map={SPHERE_TAG: HomogeneousMaterial(sigma=sigma, epsilon_r=epsilon_r)},
    )
    solver = TimeHarmonicSolver(problem, degree=degree)
    if not solve:
        sigma_field, epsilon_r_field = build_material_fields(
            solver.mesh,
            problem.material,
            cell_tags=problem.cell_tags,
            material_map=problem.material_map,
            phantom_material=problem.phantom_material,
            phantom_tag=problem.phantom_tag,
        )
        mu_r_field = build_mu_r_field(
            solver.mesh,
            problem.material,
            cell_tags=problem.cell_tags,
            material_map=problem.material_map,
        )
        fields = SimpleNamespace(
            sigma_field=sigma_field, epsilon_r_field=epsilon_r_field, mu_r_field=mu_r_field
        )
        return {"fields": fields, "solver": solver}
    t0 = time.perf_counter()
    fields = solver.solve()
    seconds = time.perf_counter() - t0
    sar = mean_sar(
        fields.e_complex,
        sigma=fields.sigma_field,
        rho=RHO_KG_M3,
        cell_tags=cell_tags,
        subdomain_ids=SPHERE_TAG,
        comm=comm,
    )
    ez = _ez_probe(fields, comm)
    return {
        "mean_sar": float(sar["mean_sar_w_per_kg"]),
        "ez_mean": complex(np.mean(ez)),
        "fields": fields,
        "solver": solver,
        "seconds": seconds,
    }


def _summed_rss_gib(comm):
    kib = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return comm.allreduce(int(kib), op=MPI.SUM) / 1024.0**2


def _closed(sigma, epsilon_r):
    c = uniform_sphere_sar_closed_form(
        e0=E0, epsilon_r=epsilon_r, sigma=sigma, omega=OMEGA, rho=RHO_KG_M3, epsilon_0=EPSILON_0
    )
    e_in = 3.0 * E0 / (complex(epsilon_r, -sigma / (OMEGA * EPSILON_0)) + 2.0)
    c["im_re"] = e_in.imag / e_in.real
    return c


def _run_leg(degree):
    """Study (sigma 0.57) solved, then the direct call on the Study's own mesh."""
    comm = MPI.COMM_WORLD
    study = _study(degree=degree) if degree is not None else _study()
    t0 = time.perf_counter()
    res = study.solve(comm)
    study_seconds = time.perf_counter() - t0
    sar_dict = study.mean_sar()
    ez = _ez_probe(res.fields, comm)
    direct = _direct_call(res.mesh, res.cell_tags, SIGMA_HIGH, EPSILON_R_SPHERE, study.degree)
    return {
        "study": study,
        "res": res,
        "sar_dict": sar_dict,
        "sar": float(sar_dict["mean_sar_w_per_kg"]),
        "ez": ez,
        "ez_mean": complex(np.mean(ez)),
        "study_seconds": study_seconds,
        "direct": direct,
        "rss_gib": _summed_rss_gib(comm),
    }


@pytest.fixture(scope="module")
def deg1():
    return _run_leg(1)


def _upstream_identity(leg, tag, comm):
    """(U) + (P): bitwise upstream of the solve, and the SAR post-processing."""
    res = leg["res"]
    f_s = res.fields
    f_d = leg["direct"]["fields"]
    n_global = res.mesh.topology.index_map(res.mesh.topology.dim).size_global
    u1 = {
        "sigma": _dg0_owned_diff(f_s.sigma_field, f_d.sigma_field, comm),
        "eps_r": _dg0_owned_diff(f_s.epsilon_r_field, f_d.epsilon_r_field, comm),
        "mu_r": _dg0_owned_diff(f_s.mu_r_field, f_d.mu_r_field, comm),
    }
    a_s = _volume_operator(res.solver, f_s.sigma_field, f_s.epsilon_r_field, f_s.mu_r_field)
    a_d = _volume_operator(
        leg["direct"]["solver"], f_d.sigma_field, f_d.epsilon_r_field, f_d.mu_r_field
    )
    a_d_norm = float(a_d.norm(PETSc.NormType.NORM_INFINITY))
    pattern = _same_pattern(a_s, a_d, comm)
    sar_direct_pp = mean_sar(
        f_s.e_complex,
        sigma=f_s.sigma_field,
        rho=RHO_KG_M3,
        cell_tags=res.cell_tags,
        subdomain_ids=SPHERE_TAG,
        comm=comm,
    )
    if comm.rank == 0:
        for key, (dmax, ndiff, nown) in u1.items():
            print(
                f"  (U) {tag} {key:5s}: max |delta| owned = {dmax!r}, differing = {ndiff}, "
                f"owned = {nown} (global cells {n_global})"
            )
        print(f"  (U) {tag} same CSR pattern = {pattern}; ||A_direct||_inf = {a_d_norm:.6e}")
        print(
            f"  (P) {tag} Study mean_sar = {leg['sar']!r}, direct mean_sar on Study fields = "
            f"{sar_direct_pp['mean_sar_w_per_kg']!r}"
        )
    for key, (dmax, ndiff, nown) in u1.items():
        assert dmax == 0.0 and ndiff == 0, f"(U) {tag} {key}: max |delta| {dmax!r}, {ndiff} differ"
        assert nown == n_global, f"(U) {tag} {key}: owned {nown} != {n_global} cells"
    assert a_d_norm > 0.0
    assert pattern, f"(U) {tag}: operators differ in nonzero pattern"
    norm, changed, max_abs = _difference(a_s, a_d, comm)
    if comm.rank == 0:
        print(
            f"  (U) {tag} ||A_study - A_direct||_inf = {norm!r}, differing entries = {changed}, "
            f"max |entry delta| = {max_abs!r} (predicted 0.0)"
        )
    assert norm == 0.0, f"(U) {tag}: ||A_study - A_direct||_inf = {norm!r}"
    assert leg["sar_dict"] == sar_direct_pp, f"(P) {tag}: {leg['sar_dict']} != {sar_direct_pp}"
    a_s.destroy()
    return a_d, a_d_norm


def _end_to_end(leg, tag, comm):
    """(I): Study vs direct within L0_BAND; returns (d_sar, d_ez)."""
    d = leg["direct"]
    d_sar = abs(leg["sar"] - d["mean_sar"]) / abs(d["mean_sar"])
    d_ez = abs(leg["ez_mean"] - d["ez_mean"]) / abs(d["ez_mean"])
    if comm.rank == 0:
        print(
            f"  (I) {tag} SAR study = {leg['sar']:.10e}, direct = {d['mean_sar']:.10e}, "
            f"rel = {d_sar:.3e}; <E_z> rel = {d_ez:.3e} (band L0_BAND {L0_BAND:.0e}; "
            f"deg-1 floor record {SOLVE_FLOOR_RECORD[0]:.3e} … {SOLVE_FLOOR_RECORD[1]:.3e})"
        )
    assert d_sar <= L0_BAND, f"(I) {tag} SAR rel {d_sar:.3e} > {L0_BAND:.0e}"
    assert d_ez <= L0_BAND, f"(I) {tag} <E_z> rel {d_ez:.3e} > {L0_BAND:.0e}"
    return d_sar, d_ez


def _closed_form_band(leg, tag, comm):
    """(G): inside MAT-4's imported band."""
    c = _closed(SIGMA_HIGH, EPSILON_R_SPHERE)
    err = abs(leg["sar"] - c["sar_w_per_kg"]) / c["sar_w_per_kg"]
    ez_mean = leg["ez_mean"]
    im_re = ez_mean.imag / ez_mean.real
    im_re_err = abs(im_re - c["im_re"]) / abs(c["im_re"])
    spread = float(np.max(np.abs(leg["ez"] - ez_mean)) / abs(ez_mean))
    vol = leg["sar_dict"]["volume_m3"] / (4.0 / 3.0 * np.pi * SPHERE_RADIUS**3)
    if comm.rank == 0:
        print(
            f"  (G) {tag} closed SAR = {c['sar_w_per_kg']:.6e}, FEM = {leg['sar']:.6e} "
            f"err {err:.3%}; Im/Re closed {c['im_re']:.4f}, FEM {im_re:.4f} err {im_re_err:.3%}; "
            f"spread {spread:.3%}; V_mesh/V_exact = {vol:.6f} (bands {MAT4_BAND:.0%} / "
            f"{MAT4_BAND:.0%} / {SPREAD_BAND:.0%})"
        )
    assert err <= MAT4_BAND, f"(G) {tag} mean SAR {err:.2%} off the closed form"
    assert im_re_err <= MAT4_BAND, f"(G) {tag} Im/Re {im_re_err:.2%} off the closed form"
    assert spread < SPREAD_BAND, f"(G) {tag} interior spread {spread:.2%}"
    return err, im_re_err


def _sigma_perturbed_operator(leg, a_d, a_d_norm, tag, comm):
    """(c2): sigma * (1 + 1e-9) on the direct operator; returns ||A'' - A||_inf."""
    f_d = leg["direct"]["fields"]
    sigma_c2 = fem.Function(f_d.sigma_field.function_space, name="sigma_c2")
    sigma_c2.x.array[:] = f_d.sigma_field.x.array * (1.0 + SIGMA_REL_PERTURBATION)
    sigma_c2.x.scatter_forward()
    a_c2 = _volume_operator(leg["direct"]["solver"], sigma_c2, f_d.epsilon_r_field, f_d.mu_r_field)
    pattern = _same_pattern(a_c2, a_d, comm)
    assert pattern, f"(c2) {tag}: perturbed operator changed its nonzero pattern"
    norm, changed, _ = _difference(a_c2, a_d, comm)
    if comm.rank == 0:
        print(
            f"  (c2) {tag} sigma*(1+{SIGMA_REL_PERTURBATION:g}): ||A'' - A_direct||_inf = {norm!r}, "
            f"relative {norm / a_d_norm:.3e}, differing entries = {changed}"
        )
    a_c2.destroy()
    return norm


@complex_only
def test_wf15_imposed_field_is_mat4_wall_data_bitwise():
    """(F): the Study's Dirichlet callable equals `MAT-4`'s, bitwise."""
    comm = MPI.COMM_WORLD
    study = _study()
    mine = study.dirichlet_field()
    theirs = _exact_exterior_numpy(SIGMA_HIGH)
    msh, _, _ = GEOMETRY.build(comm)
    n_owned_nodes = msh.geometry.index_map().size_local
    nodes = np.ascontiguousarray(msh.geometry.x[:n_owned_nodes].T)
    h = BOX_HALF_WIDTH
    corners = np.array(
        [[sx * h, sy * h, sz * h] for sx in (-1, 1) for sy in (-1, 1) for sz in (-1, 1)]
    ).T
    sets = {
        "interior probes": np.ascontiguousarray(_interior_probe_points().T),
        "box corners": corners,
        "owned geometry nodes": nodes,
    }
    verdicts = {}
    for name, pts in sets.items():
        a = mine(pts)
        b = theirs(pts)
        local = bool(a.dtype == b.dtype and a.shape == b.shape and np.array_equal(a, b))
        verdicts[name] = (
            bool(comm.allreduce(local, op=MPI.LAND)),
            int(comm.allreduce(pts.shape[1], op=MPI.SUM)),
        )
    if comm.rank == 0:
        print(f"\n[WF-15] (F) epsilon_c = {study.epsilon_c()!r}")
        for name, (ok, n) in verdicts.items():
            print(f"  (F) {name}: array_equal = {ok} over {n} points")
    for name, (ok, n) in verdicts.items():
        assert ok, f"(F) the Study's wall data differs from MAT-4's on {name}"
        assert n > 0


def test_wf15_study_default_degree_is_two():
    """(Deg), structural half: the default is 2; the explicit argument is kept."""
    assert _study().degree == 2
    assert _study(degree=1).degree == 1
    assert Study(GEOMETRY, "mid", QuasiStaticSphereField(), frequency_hz=FREQUENCY_HZ).degree == 2


@complex_only
def test_wf15_degree1_reproduces_mat4_record(deg1):
    """Record leg (degree=1): (M), (U), (P), (D1), (I), (R), (G), (c1), (c2)."""
    comm = MPI.COMM_WORLD
    leg = deg1
    res = leg["res"]
    n_study = res.mesh.topology.index_map(res.mesh.topology.dim).size_global
    msh_d, _, _ = build_preset_sphere_mesh(comm)
    n_direct = msh_d.topology.index_map(msh_d.topology.dim).size_global
    space_degree = res.fields.e_complex.function_space.ufl_element().degree
    if comm.rank == 0:
        print(
            f"\n[WF-15] degree-1 record leg: (M) Study cells {n_study}, direct generator "
            f"{n_direct} (version-tagged record 74020); space degree {space_degree}; "
            f"Study solve+mesh {leg['study_seconds']:.1f} s, direct solve "
            f"{leg['direct']['seconds']:.1f} s, summed ru_maxrss {leg['rss_gib']:.2f} GiB"
        )
    assert n_study == n_direct, f"(M) {n_study} != {n_direct}"
    assert space_degree == 1

    a_d, a_d_norm = _upstream_identity(leg, "deg1", comm)

    # (D1) the direct copy is faithful to MAT-1's map record.
    d1 = abs(leg["direct"]["mean_sar"] - MAT4_DEG1_RECORD) / MAT4_DEG1_RECORD
    # (R) the Study reproduces MAT-4 step 1's figure.
    r1 = abs(leg["sar"] - MAT4_DEG1_RECORD) / MAT4_DEG1_RECORD
    r2 = f"{leg['sar']:.6e}"
    if comm.rank == 0:
        print(
            f"  (D1) direct SAR = {leg['direct']['mean_sar']:.10e} vs record "
            f"{MAT4_DEG1_RECORD:.9e}: rel {d1:.3e}"
        )
        print(f"  (R1) Study SAR rel to record {r1:.3e}; (R2) '{r2}' vs '{MAT4_DEG1_PRINT}'")
    assert d1 <= L0_BAND, f"(D1) direct copy drifted: rel {d1:.3e}"
    _end_to_end(leg, "deg1", comm)
    assert r1 <= L0_BAND, f"(R1) rel {r1:.3e}"
    assert r2 == MAT4_DEG1_PRINT, f"(R2) {r2} != {MAT4_DEG1_PRINT}"
    _closed_form_band(leg, "deg1", comm)

    # (c1) a 1-ULP sigma through the Study's material path, no solve.
    sigma_ulp = float(np.nextafter(SIGMA_HIGH, 1.0))
    background = HomogeneousMaterial(sigma=0.0, epsilon_r=1.0)
    sigma_c1, _ = build_material_fields(
        res.mesh,
        background,
        cell_tags=res.cell_tags,
        phantom_material=GelledSalinePhantomMaterial(sigma_ulp, EPSILON_R_SPHERE, FREQUENCY_HZ),
        phantom_tag=SPHERE_TAG,
    )
    c1 = _dg0_owned_diff(sigma_c1, leg["direct"]["fields"].sigma_field, comm)
    if comm.rank == 0:
        print(
            f"  (c1) sigma = {sigma_ulp!r}: max |d sigma| = {c1[0]!r} (predicted 2^-53 = "
            f"{2.0**-53!r}), differing = {c1[1]}"
        )
    assert c1[0] > 0.0 and c1[1] > 0, f"(c1) blind: {c1}"

    c2 = _sigma_perturbed_operator(leg, a_d, a_d_norm, "deg1", comm)
    assert c2 > 0.0, "(c2) deg1: a 1e-9 relative sigma reads 0 on the operator"
    a_d.destroy()


@complex_only
def test_wf15_preset_leg_reproduces_mat1_record(deg1):
    """(Pre): the "mid" preset by name reproduces MAT-1's record; (c3) asserted."""
    comm = MPI.COMM_WORLD
    study = _study(material="mid", degree=1)
    t0 = time.perf_counter()
    study.solve(comm)
    seconds = time.perf_counter() - t0
    sar = float(study.mean_sar()["mean_sar_w_per_kg"])
    mat = study.material
    c = _closed(mat.sigma, mat.epsilon_r)
    err = abs(sar - c["sar_w_per_kg"]) / c["sar_w_per_kg"]
    sar_direct = deg1["direct"]["mean_sar"]
    c3 = abs(sar - sar_direct) / sar_direct
    if comm.rank == 0:
        print(
            f"\n[WF-15] (Pre) preset 'mid' sigma = {mat.sigma}, eps_r = {mat.epsilon_r}: "
            f"SAR = {sar:.10e} ('{sar:.4e}' vs MAT-1 '{MAT1_MID_PRINT}'), closed "
            f"{c['sar_w_per_kg']:.6e}, err {err:.3%} (recorded 3.548%); solve+mesh {seconds:.1f} s"
        )
        print(
            f"  (c3) |SAR_mid - SAR_direct(0.57)| / SAR_direct = {c3:.4%} "
            f"(ceiling 16.99 %, {c3 / L0_BAND:.2e} x L0_BAND)"
        )
    assert (mat.sigma, mat.epsilon_r) == (0.75, 78.0)
    assert f"{sar:.4e}" == MAT1_MID_PRINT, f"(Pre) {sar:.4e} != {MAT1_MID_PRINT}"
    assert err <= MAT4_BAND, f"(Pre) {err:.2%} off the closed form"
    assert c3 > L0_BAND, f"(c3) the end-to-end comparator is blind: {c3:.3e}"


def _global_dofs(V):
    return int(V.dofmap.index_map.size_global * V.dofmap.index_map_bs)


@complex_only
def test_wf15_default_degree2_study_leg():
    """W2a (§9 item 251): one degree-2 Study solve. (Deg), (U), (P), (G) asserted;
    the direct path is built without a solve; (c2), (c4) printed; hand-off written."""
    comm = MPI.COMM_WORLD
    if comm.rank == 0 and os.path.exists(DEG2_HANDOFF):
        os.remove(DEG2_HANDOFF)
    comm.Barrier()

    study = _study()  # the default degree
    t0 = time.perf_counter()
    res = study.solve(comm)
    study_seconds = time.perf_counter() - t0
    sar_dict = study.mean_sar()
    ez = _ez_probe(res.fields, comm)
    leg = {
        "study": study,
        "res": res,
        "sar_dict": sar_dict,
        "sar": float(sar_dict["mean_sar_w_per_kg"]),
        "ez": ez,
        "ez_mean": complex(np.mean(ez)),
        "direct": _direct_call(
            res.mesh, res.cell_tags, SIGMA_HIGH, EPSILON_R_SPHERE, study.degree, solve=False
        ),
    }
    rss_gib = _summed_rss_gib(comm)
    V = res.fields.e_complex.function_space
    space_degree = V.ufl_element().degree
    ndofs = _global_dofs(V)
    n_cells = int(res.mesh.topology.index_map(res.mesh.topology.dim).size_global)
    if comm.rank == 0:
        print(
            f"\n[WF-15] degree-2 Study leg: Study.degree = {study.degree}, "
            f"res.degree = {res.degree}, space degree {space_degree}; {n_cells} cells, "
            f"{ndofs} N1curl dofs ({ndofs / n_cells:.2f} per cell); "
            f"Study solve+mesh {study_seconds:.1f} s (gate 450 s), summed ru_maxrss "
            f"{rss_gib:.2f} GiB"
        )
    assert study.degree == 2 and res.degree == 2
    assert space_degree == 2, f"(Deg) solved space degree {space_degree}"

    a_d, a_d_norm = _upstream_identity(leg, "deg2", comm)
    err, im_re_err = _closed_form_band(leg, "deg2", comm)

    # Printed, not asserted (predicted).
    c2 = _sigma_perturbed_operator(leg, a_d, a_d_norm, "deg2", comm)
    a_d.destroy()
    c4 = leg["sar"] / MAT4_DEG1_RECORD - 1.0
    rss_end = _summed_rss_gib(comm)

    # Hand-off: after every reduction, rank 0 only, no comm call inside.
    payload = {
        "n_cells": n_cells,
        "n_dofs": ndofs,
        "mean_sar": repr(leg["sar"]),
        "ez_mean_real": repr(leg["ez_mean"].real),
        "ez_mean_imag": repr(leg["ez_mean"].imag),
        "study_seconds": study_seconds,
        "written_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
    }
    if comm.rank == 0:
        print(
            f"  (c2) deg2 ||A'' - A_direct||_inf > 0: {c2 > 0.0} (predicted O(1e-9), printed only)\n"
            f"  (c4) re-posed: SAR_d2 / MAT4_DEG1_RECORD - 1 = {c4:+.4%} "
            f"(predicted [-3.5, 0] %, printed only)\n"
            f"  degree-2 closed-form error {err:.4%} (predicted [0.2, 4] %), Im/Re error "
            f"{im_re_err:.4%} (predicted <= 0.6 %)\n"
            f"  WF-15 record: degree-2 mean SAR = {leg['sar']:.10e} W/kg; summed ru_maxrss "
            f"{rss_end:.2f} GiB"
        )
        os.makedirs(os.path.dirname(DEG2_HANDOFF), exist_ok=True)
        with open(DEG2_HANDOFF, "w") as fh:
            json.dump(payload, fh, indent=1)
        print(f"  hand-off written to {DEG2_HANDOFF}: {payload}")
    comm.Barrier()


@complex_only
def test_wf15_default_degree2_direct_leg():
    """W2b (§9 item 251): one direct degree-2 solve on the Study's own generator;
    precondition (hand-off file, equal cell / dof counts) and (I) at degree 2."""
    comm = MPI.COMM_WORLD
    exists = os.path.exists(DEG2_HANDOFF)
    assert exists, f"precondition: {DEG2_HANDOFF} missing (run W2a first)"
    with open(DEG2_HANDOFF) as fh:
        handoff = json.load(fh)
    msh, cell_tags, _ = GEOMETRY.build(comm)
    n_cells = int(msh.topology.index_map(msh.topology.dim).size_global)
    direct = _direct_call(msh, cell_tags, SIGMA_HIGH, EPSILON_R_SPHERE, 2)
    ndofs = _global_dofs(direct["fields"].e_complex.function_space)
    rss_gib = _summed_rss_gib(comm)
    sar_study = float(handoff["mean_sar"])
    ez_study = complex(float(handoff["ez_mean_real"]), float(handoff["ez_mean_imag"]))
    if comm.rank == 0:
        print(
            f"\n[WF-15] degree-2 direct leg: hand-off {handoff}\n"
            f"  precondition: cells {n_cells} vs {handoff['n_cells']}, dofs {ndofs} vs "
            f"{handoff['n_dofs']}; direct solve {direct['seconds']:.1f} s, summed ru_maxrss "
            f"{rss_gib:.2f} GiB"
        )
    assert n_cells == handoff["n_cells"], f"precondition: cells {n_cells} != {handoff['n_cells']}"
    assert ndofs == handoff["n_dofs"], f"precondition: dofs {ndofs} != {handoff['n_dofs']}"
    leg = {"sar": sar_study, "ez_mean": ez_study, "direct": direct}
    _end_to_end(leg, "deg2 (cross-process)", comm)
