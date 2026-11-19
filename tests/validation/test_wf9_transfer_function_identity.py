"""`WF-9` step 0 (§9 item 183): the ISO 10974 Tier-3 identity in exact discrete form.

One structured box in `IMP-1`'s host (σ 0.57 S/m, εᵣ 78, 64 MHz), N1curl
degree 1, with and without a σ = 800 S/m square lead (2 × 2 × 24 mm, cell tag
2).  Since the only difference between the two operators is
ΔA = jωμ₀Δσ·M_lead, ``q_with − q_without`` is the with-lead response to the
impressed current ``J_eq = Δσ·E_inc`` on the lead, and since A = Aᵀ (real
basis, scalar ε_c) the 24 stepped-source transfer functions equal one
tip-source adjoint solve read along the lead (Lorentz reciprocity, Jin §6.4
eqs (6.78)–(6.85); the Hermitian route (6.86)–(6.88) is control C1).

Asserted (bands imported, see the brief ``docs/planning/chunks/WF-9.md``):
(A0) census, (A1) A = Aᵀ, (T1) superposition, (T2) reciprocity, (T3) the
row's identity −jωμ₀Δσ∫χ_lead w·E_inc = q_with − q_without, (H) held vs
fresh, and controls C1 / C2 wired.  Recorded only: the ISO path fold Q_N for
N = 3, 6, 12, 24 (its band is the weekly's).

Traps honoured: PEC mode on every problem; ``project_source=False``; zero-trace
sources only on the zero-trace factor; q/TF/folds via components (no
``ufl.inner`` conjugation); complex assembly reduced over ranks; DG0 masks
from numpy midpoints over owned + ghost cells, counts over owned cells only.
"""

from __future__ import annotations

import gc
import os
import resource
import time

import dolfinx
import numpy as np
import ufl
from dolfinx import fem
from dolfinx.fem.petsc import assemble_matrix
from mpi4py import MPI
from petsc4py import PETSc

from fem_em_solver.core import (
    HomogeneousMaterial,
    TimeHarmonicProblem,
    TimeHarmonicSolver,
)
from fem_em_solver.core.time_harmonic import build_material_fields
from tests.complex_mode import complex_only
from tests.mesh.test_birdcage_symmetry_quadrant import EXACT_IDENTITY_RTOL
from tests.validation.test_birdcage_power_identity import DISCRETE_IDENTITY_RTOL
from tests.validation.test_imp1_pec_sphere_lossy_host import EPSILON_R_HOST
from tests.validation.test_lossy_sphere_sar import E0, FREQUENCY_HZ, SIGMA_HIGH
from tests.validation.test_port_gap_voltage_impedance import SIGMA_WIRE_S_PER_M

MU_0 = 4.0e-7 * np.pi
OMEGA = 2.0 * np.pi * FREQUENCY_HZ

# A1's bound is `test_time_harmonic_mms.py:255`'s literal (no module constant).
SYMMETRY_RTOL = 1.0e-10

HALF_XY = 0.008
HALF_Z = 0.020
N_CELLS = (16, 16, 40)
CELL_COUNT = 61_440
LEAD_HALF_XY = 0.001
LEAD_HALF_Z = 0.012
N_LAYERS = 24
LEAD_CELL_COUNT = 576
BLOCK_CELL_COUNT = 192
LEAD_VOLUME_M3 = 9.6e-8
BALL_CENTRE = np.array([0.0, 0.0, 0.012])
BALL_RADIUS = 0.003

HOST_TAG = 1
LEAD_TAG = 2

DELTA_SIGMA = SIGMA_WIRE_S_PER_M - SIGMA_HIGH
J_BLOCK = 1.0  # A/m²; every recorded quantity is a ratio
LOAD = -1j * OMEGA * MU_0

TF_EXCHANGE_DIR = "/workspace/logs"


def _build():
    comm = MPI.COMM_WORLD
    msh = dolfinx.mesh.create_box(
        comm,
        [np.array([-HALF_XY, -HALF_XY, -HALF_Z]), np.array([HALF_XY, HALF_XY, HALF_Z])],
        list(N_CELLS),
        cell_type=dolfinx.mesh.CellType.tetrahedron,
    )
    tdim = msh.topology.dim
    imap = msh.topology.index_map(tdim)
    n_owned = imap.size_local
    n_local = n_owned + imap.num_ghosts
    cells = np.arange(n_local, dtype=np.int32)
    mid = dolfinx.mesh.compute_midpoints(msh, tdim, cells)
    x, y, z = mid[:, 0], mid[:, 1], mid[:, 2]

    lead = (np.abs(x) < LEAD_HALF_XY) & (np.abs(y) < LEAD_HALF_XY) & (np.abs(z) < LEAD_HALF_Z)
    block = (x > 0.004) & (x < 0.006) & (np.abs(y) < 0.001) & (z > 0.002) & (z < 0.010)
    ball = (~lead) & (np.linalg.norm(mid - BALL_CENTRE, axis=1) < BALL_RADIUS)
    layer_index = np.floor((z + LEAD_HALF_Z) / 1.0e-3).astype(int)
    layers = [lead & (layer_index == j) for j in range(N_LAYERS)]

    q0 = fem.functionspace(msh, ("DG", 0))
    dof_of_cell = q0.dofmap.list[:, 0]

    def mask(flags):
        f = fem.Function(q0)
        f.x.array[:] = 0.0
        f.x.array[dof_of_cell[cells]] = flags.astype(float)
        f.x.scatter_forward()
        return f

    def owned_count(flags):
        return comm.allreduce(int(np.count_nonzero(flags[:n_owned])), op=MPI.SUM)

    tag_values = np.where(lead, LEAD_TAG, HOST_TAG).astype(np.int32)
    cell_tags = dolfinx.mesh.meshtags(msh, tdim, cells, tag_values)

    return {
        "comm": comm,
        "mesh": msh,
        "n_owned": n_owned,
        "flags": {"lead": lead, "block": block, "ball": ball, "layers": layers},
        "owned_count": owned_count,
        "chi_lead": mask(lead),
        "chi_src": mask(block),
        "chi_b": mask(ball),
        "chi_layers": [mask(f) for f in layers],
        "cell_tags": cell_tags,
    }


def _scalar(comm, expr):
    return comm.allreduce(fem.assemble_scalar(fem.form(expr)), op=MPI.SUM)


def _host():
    return HomogeneousMaterial(sigma=SIGMA_HIGH, epsilon_r=EPSILON_R_HOST)


def _lead_map():
    return {LEAD_TAG: HomogeneousMaterial(sigma=SIGMA_WIRE_S_PER_M, epsilon_r=EPSILON_R_HOST)}


def _problem(fx, *, with_lead, dirichlet):
    kwargs = {}
    if dirichlet:
        kwargs["dirichlet_e_field"] = lambda xx: np.vstack(
            (np.zeros(xx.shape[1]), np.zeros(xx.shape[1]), E0 * np.ones(xx.shape[1]))
        )
    return TimeHarmonicProblem(
        mesh=fx["mesh"],
        frequency_hz=FREQUENCY_HZ,
        material=_host(),
        cell_tags=fx["cell_tags"],
        material_map=_lead_map() if with_lead else None,
        boundary_condition="pec_zero_tangential_a",
        **kwargs,
    )


def _free(solver):
    lp = solver._linear_problem
    solver._linear_problem = None
    if lp is not None:
        try:
            lp.solver.destroy()
        except Exception:
            pass
    del lp
    gc.collect()


def _rss_gib(comm):
    kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return comm.allreduce(kb, op=MPI.SUM) / 1024.0**2


def test_wf9_census():
    fx = _build()
    comm = fx["comm"]
    fl = fx["flags"]
    cnt = fx["owned_count"]
    n_cells = fx["mesh"].topology.index_map(3).size_global
    n_lead = cnt(fl["lead"])
    n_block = cnt(fl["block"])
    n_ball = cnt(fl["ball"])
    layer_counts = [cnt(f) for f in fl["layers"]]

    stacked = np.sum(np.vstack(fl["layers"]).astype(int), axis=0)
    partition_ok = comm.allreduce(
        bool(np.array_equal(stacked, fl["lead"].astype(int))), op=MPI.LAND
    )
    disjoint_local = not (
        np.any(fl["ball"] & fl["lead"])
        or np.any(fl["ball"] & fl["block"])
        or np.any(fl["lead"] & fl["block"])
    )
    disjoint = comm.allreduce(bool(disjoint_local), op=MPI.LAND)
    lead_volume = _scalar(comm, fx["chi_lead"] * ufl.dx(domain=fx["mesh"])).real
    ball_volume = _scalar(comm, fx["chi_b"] * ufl.dx(domain=fx["mesh"])).real

    if comm.rank == 0:
        print(f"\n[WF-9 step0 census] ranks={comm.size}")
        print(
            f"  cells={n_cells} lead={n_lead} block={n_block} ball={n_ball} (ball predicted [400, 700])"
        )
        print(f"  layer counts={layer_counts}")
        print(
            f"  lead volume={lead_volume:.12e} m^3 (expect {LEAD_VOLUME_M3:.3e}); ball volume={ball_volume:.6e} m^3"
        )
        print(
            f"  partition sum_j chi_j == chi_lead: {partition_ok}; masks pairwise disjoint: {disjoint}"
        )

    assert n_cells == CELL_COUNT
    assert n_lead == LEAD_CELL_COUNT
    assert n_block == BLOCK_CELL_COUNT
    assert n_ball > 0
    assert partition_ok
    assert disjoint
    assert abs(lead_volume - LEAD_VOLUME_M3) <= EXACT_IDENTITY_RTOL * LEAD_VOLUME_M3


def _run_identities():
    t_start = time.perf_counter()
    fx = _build()
    comm = fx["comm"]
    msh = fx["mesh"]
    chi_b, chi_lead, chi_src = fx["chi_b"], fx["chi_lead"], fx["chi_src"]
    chi_layers = fx["chi_layers"]
    dx = ufl.dx(domain=msh)
    zhat = ufl.as_vector((0.0, 0.0, 1.0))
    times = {}

    v_ball = _scalar(comm, chi_b * dx).real

    def q(e):
        return _scalar(comm, chi_b * e[2] * dx) / v_ball

    # ---- (A1) symmetry of the with-lead operator -------------------------------
    sym_solver = TimeHarmonicSolver(_problem(fx, with_lead=True, dirichlet=False), degree=1)
    sigma_f, eps_f = build_material_fields(
        msh, _host(), cell_tags=fx["cell_tags"], material_map=_lead_map()
    )
    mat = assemble_matrix(fem.form(sym_solver.bilinear_form(sigma_f, eps_f)))
    mat.assemble()
    a_norm = mat.norm(PETSc.NormType.FROBENIUS)
    mat_t = mat.copy()
    mat_t.transpose()
    mat_t.axpy(-1.0, mat, structure=PETSc.Mat.Structure.DIFFERENT_NONZERO_PATTERN)
    sym_norm = mat_t.norm(PETSc.NormType.FROBENIUS)
    ndofs = sym_solver.function_space().dofmap.index_map.size_global
    mat.destroy()
    mat_t.destroy()
    del sym_solver

    # ---- D1: Dirichlet E0 ẑ, lead out / in ------------------------------------
    t0 = time.perf_counter()
    s_od = TimeHarmonicSolver(_problem(fx, with_lead=False, dirichlet=True), degree=1)
    e_inc_d1 = s_od.solve().e_complex
    _free(s_od)
    times["S_oD"] = time.perf_counter() - t0

    t0 = time.perf_counter()
    s_wd = TimeHarmonicSolver(_problem(fx, with_lead=True, dirichlet=True), degree=1)
    e_with_d1 = s_wd.solve().e_complex
    _free(s_wd)
    times["S_wD"] = time.perf_counter() - t0

    # ---- D2: impressed ẑ block, zero trace, lead out ---------------------------
    t0 = time.perf_counter()
    s_o0 = TimeHarmonicSolver(_problem(fx, with_lead=False, dirichlet=False), degree=1)
    e_inc_d2 = s_o0.solve(
        current_density=lambda xx: J_BLOCK * chi_src * zhat, project_source=False
    ).e_complex
    _free(s_o0)
    times["S_o0"] = time.perf_counter() - t0

    # ---- S_w0: zero trace, lead in; adjoint solve on ℓ = ẑχ_b/V_b ---------------
    t0 = time.perf_counter()
    s_w0 = TimeHarmonicSolver(_problem(fx, with_lead=True, dirichlet=False), degree=1)
    ell = lambda v: ufl.inner((chi_b / v_ball) * zhat, v) * ufl.dx  # noqa: E731
    w = s_w0.solve(extra_linear_terms=[ell]).e_complex
    times["S_w0_factor+adjoint"] = time.perf_counter() - t0

    def held(load_vec_expr):
        """Held back-substitution on S_w0 for J (a UFL vector); load −jωμ₀∫J·v̄."""
        return s_w0.solve_with_held_factorization(
            extra_linear_terms=[lambda v: LOAD * ufl.inner(load_vec_expr, v) * ufl.dx]
        ).e_complex

    t0 = time.perf_counter()
    e_with_d2 = held(J_BLOCK * chi_src * zhat)
    tf_step = np.array([q(held(DELTA_SIGMA * chi_layers[j] * zhat)) for j in range(N_LAYERS)])
    part_d1 = np.array([q(held(DELTA_SIGMA * chi_layers[j] * e_inc_d1)) for j in range(N_LAYERS)])
    part_d2 = np.array([q(held(DELTA_SIGMA * chi_layers[j] * e_inc_d2)) for j in range(N_LAYERS)])
    times["held_73"] = time.perf_counter() - t0
    _free(s_w0)

    # ---- (H) fresh TF_23 --------------------------------------------------------
    t0 = time.perf_counter()
    s_fresh = TimeHarmonicSolver(_problem(fx, with_lead=True, dirichlet=False), degree=1)
    tf23_fresh = q(
        s_fresh.solve(
            current_density=lambda xx: DELTA_SIGMA * chi_layers[N_LAYERS - 1] * zhat,
            project_source=False,
        ).e_complex
    )
    _free(s_fresh)
    times["fresh_H"] = time.perf_counter() - t0

    # ---- read-backs (components only: no conjugation) ---------------------------
    tf_adj = np.array(
        [LOAD * DELTA_SIGMA * _scalar(comm, chi_layers[j] * w[2] * dx) for j in range(N_LAYERS)]
    )
    tf_herm = np.array(
        [
            LOAD * DELTA_SIGMA * _scalar(comm, chi_layers[j] * ufl.conj(w[2]) * dx)
            for j in range(N_LAYERS)
        ]
    )

    drives = {
        "D1": {"inc": e_inc_d1, "with": e_with_d1, "part": part_d1},
        "D2": {"inc": e_inc_d2, "with": e_with_d2, "part": part_d2},
    }
    results = {}
    for name, d in drives.items():
        q_with = q(d["with"])
        q_without = q(d["inc"])
        dq = q_with - q_without
        t1_sum = np.sum(d["part"])
        t3_fold = LOAD * DELTA_SIGMA * _scalar(comm, chi_lead * ufl.dot(w, d["inc"]) * dx)
        c2_fold = LOAD * DELTA_SIGMA * _scalar(comm, chi_lead * ufl.dot(w, d["with"]) * dx)
        transverse = (
            LOAD
            * DELTA_SIGMA
            * _scalar(comm, chi_lead * (w[0] * d["inc"][0] + w[1] * d["inc"][1]) * dx)
        )
        ez_int = np.array(
            [_scalar(comm, chi_layers[j] * d["inc"][2] * dx) for j in range(N_LAYERS)]
        )
        along_layer = np.array(
            [
                LOAD * DELTA_SIGMA * _scalar(comm, chi_layers[j] * w[2] * d["inc"][2] * dx)
                for j in range(N_LAYERS)
            ]
        )
        layer_vol = np.array([_scalar(comm, chi_layers[j] * dx).real for j in range(N_LAYERS)])
        folds = {}
        for n_seg in (3, 6, 12, 24):
            per = N_LAYERS // n_seg
            q_n = 0.0 + 0.0j
            along = 0.0 + 0.0j
            for k in range(n_seg):
                sl = slice(k * per, (k + 1) * per)
                tf_k = np.sum(tf_step[sl])
                e_bar = np.sum(ez_int[sl]) / np.sum(layer_vol[sl])
                q_n += tf_k * e_bar
                along += np.sum(along_layer[sl])
            folds[n_seg] = (q_n, along - q_n)
        results[name] = dict(
            q_with=q_with,
            q_without=q_without,
            dq=dq,
            t1_sum=t1_sum,
            t3_fold=t3_fold,
            c2_fold=c2_fold,
            transverse=transverse,
            folds=folds,
        )

    rss = _rss_gib(comm)
    elapsed = time.perf_counter() - t_start

    # W1/W2 TF agreement through a gitignored exchange file.
    tf_agreement = None
    if comm.rank == 0:
        os.makedirs(TF_EXCHANGE_DIR, exist_ok=True)
        np.save(os.path.join(TF_EXCHANGE_DIR, f"wf9_step0_tf_n{comm.size}.npy"), tf_step)
        for other in (2, 4):
            if other == comm.size:
                continue
            path = os.path.join(TF_EXCHANGE_DIR, f"wf9_step0_tf_n{other}.npy")
            if os.path.exists(path):
                ref = np.load(path)
                tf_agreement = (other, float(np.max(np.abs(tf_step - ref)) / np.max(np.abs(ref))))

    t2_dev = float(np.max(np.abs(tf_step - tf_adj)) / np.max(np.abs(tf_adj)))
    c1_dev = float(np.max(np.abs(tf_herm - tf_adj)) / np.max(np.abs(tf_adj)))
    h_dev = float(abs(tf23_fresh - tf_step[-1]) / abs(tf_step[-1]))

    if comm.rank == 0:
        print(f"\n[WF-9 step0 identities] ranks={comm.size} dofs={ndofs}")
        print(f"  (A1) ||A-A^T||_F/||A||_F = {sym_norm / a_norm:.3e}  (||A||_F={a_norm:.6e})")
        print(f"  (T2) max|TF_step-TF_adj|/max|TF_adj| = {t2_dev:.3e}")
        print(f"  (H)  |TF23_fresh-TF23_held|/|TF23| = {h_dev:.3e}")
        print(f"  (C1) Hermitian read-back deviation = {c1_dev:.6%}  (predicted [50%, 200%])")
        for name, r in results.items():
            dq = r["dq"]
            print(f"  --- {name} ---")
            print(f"    q_without={r['q_without']:.10e} q_with={r['q_with']:.10e}")
            print(
                f"    |q_with|/|q_without| = {abs(r['q_with']) / abs(r['q_without']):.6f}"
                + ("  (D1 predicted [1.2, 10])" if name == "D1" else "")
            )
            print(f"    (T1) |sum_j - dq|/|dq| = {abs(r['t1_sum'] - dq) / abs(dq):.3e}")
            print(f"    (T3) |fold(w,E_inc) - dq|/|dq| = {abs(r['t3_fold'] - dq) / abs(dq):.3e}")
            print(
                f"    (C2) |fold(w,E_with) - dq|/|dq| = {abs(r['c2_fold'] - dq) / abs(dq):.6%}  (predicted [80%, 105%])"
            )
            print(f"    transverse term / dq = {abs(r['transverse']) / abs(dq):.6e}")
            print(f"    Q_exact = {dq:.10e}")
            for n_seg, (q_n, cov) in r["folds"].items():
                print(
                    f"    N={n_seg:2d} Q_N={q_n:.10e} |Q_N-Q|/|Q|={abs(q_n - dq) / abs(dq):.6e} "
                    f"cov/|Q|={abs(cov) / abs(dq):.6e} |Q_N|^2/|Q|^2={abs(q_n) ** 2 / abs(dq) ** 2:.6f}"
                )
        print("  TF_j profile (j, z_mid mm, |TF|, arg deg):")
        for j in range(N_LAYERS):
            print(
                f"    {j:2d} {-11.5 + j:6.1f} {abs(tf_step[j]):.6e} {np.degrees(np.angle(tf_step[j])):9.4f}"
            )
        print("  times: " + ", ".join(f"{k}={v:.2f}s" for k, v in times.items()))
        print(f"  summed max RSS over ranks = {rss:.3f} GiB; elapsed = {elapsed:.2f} s")
        if tf_agreement is None:
            print("  W1/W2 TF agreement: no other-width TF file yet")
        else:
            print(
                f"  W1/W2 TF agreement vs n={tf_agreement[0]}: max rel = {tf_agreement[1]:.3e} (predicted <= 1e-10)"
            )

    assert sym_norm <= SYMMETRY_RTOL * a_norm
    assert t2_dev <= DISCRETE_IDENTITY_RTOL
    assert h_dev <= DISCRETE_IDENTITY_RTOL
    assert (
        c1_dev > DISCRETE_IDENTITY_RTOL
    ), "C1 not wired: Hermitian read-back equals the adjoint TF"
    for name, r in results.items():
        dq = r["dq"]
        assert abs(r["t1_sum"] - dq) <= DISCRETE_IDENTITY_RTOL * abs(dq), f"T1 {name}"
        assert abs(r["t3_fold"] - dq) <= DISCRETE_IDENTITY_RTOL * abs(dq), f"T3 {name}"
        assert abs(r["c2_fold"] - dq) > DISCRETE_IDENTITY_RTOL * abs(dq), f"C2 {name} not wired"

    return {
        "fx": fx,
        "q": q,
        "ell": ell,
        "w": w,
        "drives": drives,
        "results": results,
        "tf_step": tf_step,
        "tf_adj": tf_adj,
        "tf_herm": tf_herm,
        "ndofs": ndofs,
        "sym_ratio": float(sym_norm / a_norm),
        "t2_dev": t2_dev,
        "h_dev": h_dev,
        "c1_dev": c1_dev,
        "v_ball": v_ball,
        "times": times,
        "rss": rss,
        "elapsed": elapsed,
    }


@complex_only
def test_wf9_identities():
    _run_identities()
