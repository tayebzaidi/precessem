"""`WF-9` step 1 Part B (§9 item 224, chain I8): the lead in the coil-driven phantom.

The F-small untuned 4-leg birdcage at 64 MHz (50 Ohm lumped sheets, degree 1)
with `GEO-36` step 3's capped rod (r = 1 mm, tip to tip 32 mm, centre
(0.015, 0, 0)) **kept** as cell tag `BIRDCAGE_IMPLANT_CELL_TAG` and switched by
material only: σ = 800 S/m (with) against the gel (without), εᵣ and μᵣ the
gel's in both.  One mesh, so ΔA = jωμ₀Δσ·M_lead exactly and, since A = Aᵀ
(real basis, symmetric sheet terms), the adjoint tip solve w folded with the
lead-free field equals the direct with-lead tip response (Jin §6.4 eqs
(6.78)–(6.85); the Hermitian route (6.86)–(6.88) is control C1).

Asserted (bands imported; brief ``docs/planning/chunks/WF-9.md`` § item 224):
(A0) census, (M0) material read-back on the implant cells, (A1) A = Aᵀ,
(i′) loaded reciprocity both arms, (T3) the Done-when identity per drive,
(H) held vs fresh, (PF) the 1-D path fold at N = 32 within `PATH_FOLD_BAND`,
(BF) the Born fold misses (> `PATH_FOLD_BAND`), controls C1 / C2 wired.

Not covered: a measured / ISO 10974-calibrated TF, heating or SAR, a PEC or
insulated or real lead, h-convergence, quadrature drives, F-human, Palace.
"""

from __future__ import annotations

import gc
import os
import resource
import time
from dataclasses import replace

import numpy as np
import pytest
import ufl
from dolfinx import fem
from dolfinx.fem.petsc import assemble_matrix
from mpi4py import MPI
from petsc4py import PETSc

from fem_em_solver.core import HomogeneousMaterial, TimeHarmonicSolver
from fem_em_solver.core.time_harmonic import build_material_fields
from fem_em_solver.io.mesh import BIRDCAGE_IMPLANT_CELL_TAG
from fem_em_solver.ports.lumped import (
    lumped_port_bilinear_term,
    lumped_port_linear_term,
    sheet_terminal_current,
)
from fem_em_solver.utils.constants import MU_0
from tests.complex_mode import complex_only
from tests.mesh.test_birdcage_implant_rod import (
    KEPT_ROD_CELL_GUARD,
    KEPT_ROD_CELL_RECORD,
    _build_kept_rod,
    _cell_values_with_ghosts,
)
from tests.mesh.test_birdcage_port_sheet_prerequisite import CELL_COUNT_BAND
from tests.validation.test_birdcage_b1_plus_map import _solve_driven
from tests.validation.test_birdcage_power_identity import (
    DISCRETE_IDENTITY_RTOL,
    STEP2D_FREQUENCY_HZ,
)
from tests.validation.test_imp3_rod_birdcage import TIPS
from tests.validation.test_port_birdcage_four_port import build_four_port_sweep
from tests.validation.test_port_birdcage_lumped_column import PHANTOM_CELL_TAG
from tests.validation.test_port_gap_voltage_impedance import SIGMA_WIRE_S_PER_M
from tests.validation.test_port_lumped_sheet_sweep import RECIPROCITY_BAND
from tests.validation.test_wf9_transfer_function_identity import (
    BALL_RADIUS,
    SYMMETRY_RTOL,
    _free,
)

#: The 1-D path-fold band: PROJECT_PLAN §10, the 2026-10-03 weekly's `WF-9`
#: Done-when ruling ("the 1-D path fold ∫TF(z)·E_tan(z) dz reproduces it
#: within 5 %"); the Born fold must miss by more than the same band.
PATH_FOLD_BAND = 0.05

OMEGA = 2.0 * np.pi * float(STEP2D_FREQUENCY_HZ)
LOAD = -1j * OMEGA * MU_0
PORT_IDS = ("P1", "P2", "P3", "P4")

N_SLABS = 32
SLAB_Z0 = -0.016
SLAB_DZ = 0.001
FOLD_NS = (4, 8, 16, 32)

TF_EXCHANGE_DIR = "/workspace/logs"

# --- `WF-9` step 2 (§9 item 276, chain I16): the same fixture at degree 2 ----
#: Opt-in Nédélec degree for every solve of this module.  Unset ⇒ 1, so step
#: 1's nodes run byte for byte as before; 2 runs `test_wf9_step2_degree2`.
DEGREE_ENV = "FEM_EM_WF9_DEGREE"

#: Step 1's degree-1 global dofs on this fixture
#: (`20261004T201655Z_WF-9-step1-W1.log:1895`); (Deg) requires > 4× this.
DEG1_DOFS_RECORD = 216_731

#: Step 1's degree-1 (q_without, q_with) per drive, copied at `.10e` from
#: `20261004T201655Z_WF-9-step1-W1.log:1900, 1913, 1926, 1939`.  Printed
#: order move only; never asserted.
STEP1_Q_RECORD = {
    "P1": (-1.2753864198e-01 - 1.8713433801e-01j, -3.1505058232e-01 - 4.2222761708e-01j),
    "P2": (-1.9768599288e-02 - 2.2756234126e-02j, -4.8687195593e-02 - 5.2392004656e-02j),
    "P3": (4.1543136576e-02 + 4.8934145739e-02j, 1.0469981634e-01 + 1.1224765880e-01j),
    "P4": (-2.1013373374e-02 - 2.4680695145e-02j, -5.0248941990e-02 - 5.4629859795e-02j),
}

#: Printed-only predictions for step 2 (brief § item 276).
PRED_ORDER_MOVE = (0.005, 0.30)
PRED_RSS_P2_N8 = (20.0, 32.0)


def _degree():
    raw = os.environ.get(DEGREE_ENV, "").strip()
    if not raw:
        return 1
    d = int(raw)
    if d not in (1, 2):
        raise ValueError(f"{DEGREE_ENV}={raw!r}: only 1 or 2")
    return d


def _space_degree(space):
    return int(space.ufl_element().degree)


# --- predictions (printed, never asserted; brief § item 224) -----------------
PRED_PF = {"P1": (1e-4, 2e-2), "P3": (1e-4, 2e-2), "P2": (3e-3, 4e-2), "P4": (3e-3, 4e-2)}
PRED_BF = (1.0, 30.0)
PRED_C1 = (0.10, 2.00)
PRED_C2 = (0.80, 1.05)
PRED_RATIO = {"P1": (1.5, 6.0), "P3": (1.3, 5.0), "P2": (0.8, 3.0), "P4": (0.8, 3.0)}


def _say(msg):
    if MPI.COMM_WORLD.rank == 0:
        print(msg, flush=True)


def _in(x, band):
    return "IN" if band[0] <= x <= band[1] else "OUT"


def _scalar(comm, expr):
    return comm.allreduce(fem.assemble_scalar(fem.form(expr)), op=MPI.SUM)


def _rss_gib(comm):
    kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return comm.allreduce(float(kb), op=MPI.SUM) / 1024.0**2


def _build_case():
    comm = MPI.COMM_WORLD
    t0 = time.perf_counter()
    built = _build_kept_rod()
    msh = built[0]
    tdim = msh.topology.dim
    n_cells = int(msh.topology.index_map(tdim).size_global)
    assert n_cells <= KEPT_ROD_CELL_GUARD, (n_cells, KEPT_ROD_CELL_GUARD)
    t_mesh = time.perf_counter() - t0

    case = build_four_port_sweep(
        frequency_hz=STEP2D_FREQUENCY_HZ, built=built, build_only=True, degree=_degree()
    )
    problem = case["problem"]
    mmap = dict(problem.material_map)
    gel = mmap[PHANTOM_CELL_TAG]
    delta_sigma = float(SIGMA_WIRE_S_PER_M - gel.sigma)
    lead_mat = HomogeneousMaterial(sigma=SIGMA_WIRE_S_PER_M, epsilon_r=gel.epsilon_r, mu_r=gel.mu_r)
    p_with = replace(problem, material_map={**mmap, BIRDCAGE_IMPLANT_CELL_TAG: lead_mat})
    p_without = replace(problem, material_map={**mmap, BIRDCAGE_IMPLANT_CELL_TAG: gel})

    imap = msh.topology.index_map(tdim)
    n_owned = imap.size_local
    n_local = n_owned + imap.num_ghosts
    cells = np.arange(n_local, dtype=np.int32)
    mid = msh.geometry.x[msh.geometry.dofmap[cells]].mean(axis=1)
    tag = _cell_values_with_ghosts(msh, built[1])
    lead = tag == BIRDCAGE_IMPLANT_CELL_TAG
    ball = (tag == PHANTOM_CELL_TAG) & (
        np.linalg.norm(mid - np.asarray(TIPS["top"]), axis=1) < BALL_RADIUS
    )
    slab_idx = np.floor((mid[:, 2] - SLAB_Z0) / SLAB_DZ).astype(int)
    slabs = [lead & (slab_idx == k) for k in range(N_SLABS)]

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

    return {
        "comm": comm,
        "mesh": msh,
        "built": built,
        "sweep": case,
        "n_cells": n_cells,
        "n_owned": n_owned,
        "cells": cells,
        "dof_of_cell": dof_of_cell,
        "flags": {"lead": lead, "ball": ball, "slabs": slabs},
        "owned_count": owned_count,
        "chi_lead": mask(lead),
        "chi_b": mask(ball),
        "chi_k": [mask(f) for f in slabs],
        "gel": gel,
        "delta_sigma": delta_sigma,
        "p_with": p_with,
        "p_without": p_without,
        "t_mesh": t_mesh,
    }


@pytest.fixture(scope="module")
def _case():
    fx = _build_case()
    yield fx
    gc.collect()


def _material_minmax(fx, problem):
    """(M0): rank-reduced min/max DG0 σ and εᵣ on owned implant cells."""
    comm = fx["comm"]
    sigma_f, eps_f = build_material_fields(
        fx["mesh"],
        problem.material,
        cell_tags=problem.cell_tags,
        material_map=problem.material_map,
        phantom_material=problem.phantom_material,
        phantom_tag=problem.phantom_tag,
    )
    owned_lead = np.flatnonzero(fx["flags"]["lead"][: fx["n_owned"]])
    out = {}
    for name, f in (("sigma", sigma_f), ("eps", eps_f)):
        dofs = f.function_space.dofmap.list[:, 0]
        vals = np.real(f.x.array[dofs[owned_lead]])
        lo = float(vals.min()) if vals.size else np.inf
        hi = float(vals.max()) if vals.size else -np.inf
        out[name] = (comm.allreduce(lo, op=MPI.MIN), comm.allreduce(hi, op=MPI.MAX))
    return out


@complex_only
def test_wf9_step1_census(_case):
    fx = _case
    comm = fx["comm"]
    fl = fx["flags"]
    cnt = fx["owned_count"]
    n_lead = cnt(fl["lead"])
    n_ball = cnt(fl["ball"])
    slab_counts = [cnt(f) for f in fl["slabs"]]
    stacked = np.sum(np.vstack(fl["slabs"]).astype(int), axis=0)
    partition_ok = comm.allreduce(
        bool(np.array_equal(stacked, fl["lead"].astype(int))), op=MPI.LAND
    )
    disjoint = comm.allreduce(bool(not np.any(fl["ball"] & fl["lead"])), op=MPI.LAND)
    m_with = _material_minmax(fx, fx["p_with"])
    m_without = _material_minmax(fx, fx["p_without"])
    gel = fx["gel"]
    rel_cells = abs(fx["n_cells"] / KEPT_ROD_CELL_RECORD - 1.0)

    _say(f"\n[WF-9 step1 census] ranks={comm.size}")
    _say(
        f"  (A0) cells={fx['n_cells']} vs KEPT_ROD_CELL_RECORD {KEPT_ROD_CELL_RECORD} "
        f"rel {rel_cells:.3e} (band {CELL_COUNT_BAND:g}); mesh+guard {fx['t_mesh']:.1f} s"
    )
    _say(f"  implant={n_lead} ball={n_ball} slabs={slab_counts}")
    _say(f"  sum_k chi_k == chi_lead cell by cell: {partition_ok}; ball/lead disjoint: {disjoint}")
    _say(
        f"  (M0) with: sigma [{m_with['sigma'][0]:.6g}, {m_with['sigma'][1]:.6g}] eps_r "
        f"[{m_with['eps'][0]:.6g}, {m_with['eps'][1]:.6g}]; without: sigma "
        f"[{m_without['sigma'][0]:.6g}, {m_without['sigma'][1]:.6g}] eps_r "
        f"[{m_without['eps'][0]:.6g}, {m_without['eps'][1]:.6g}]; gel sigma {gel.sigma} "
        f"eps_r {gel.epsilon_r}; delta_sigma {fx['delta_sigma']:.6f}"
    )

    assert rel_cells <= CELL_COUNT_BAND
    assert n_lead > 0 and n_ball > 0
    assert all(c > 0 for c in slab_counts), slab_counts
    assert partition_ok
    assert disjoint
    assert m_with["sigma"] == (SIGMA_WIRE_S_PER_M, SIGMA_WIRE_S_PER_M), m_with
    assert m_with["eps"] == (gel.epsilon_r, gel.epsilon_r), m_with
    assert m_without["sigma"] == (gel.sigma, gel.sigma), m_without
    assert m_without["eps"] == (gel.epsilon_r, gel.epsilon_r), m_without


def _sheets(specs, driven=None):
    return [spec.sheet(driven=(spec.port_id == driven)) for spec in specs]


def _arm(fx, problem, *, first, extra_loads):
    """Factor once with ``first`` (a drive id or the tip load ``"ell"``); held
    back-substitutions for every other entry of ``extra_loads``.

    Returns {key: e_complex} and per-drive terminal currents.
    """
    sweep = fx["sweep"]
    msh = fx["mesh"]
    comm = fx["comm"]
    tags_f = sweep["facet_tags"]
    specs = sweep["specs"]
    held_sheets = _sheets(specs)
    bil = [
        (
            lambda trial, test, _s=s: lumped_port_bilinear_term(
                msh, tags_f, _s, trial, test, omega_rad_per_s=OMEGA
            )
        )
        for s in held_sheets
    ]

    def load_of(key):
        if key == "ell":
            return [extra_loads["ell"]]
        drv = next(s for s in _sheets(specs, key) if s.port_id == key)
        return [
            lambda v, _s=drv: lumped_port_linear_term(msh, tags_f, _s, v, omega_rad_per_s=OMEGA)
        ]

    solver = TimeHarmonicSolver(problem, degree=int(sweep["degree"]))
    fx.setdefault("space_degrees", []).append(_space_degree(solver.function_space()))
    fields = {}
    t0 = time.perf_counter()
    fields[first] = solver.solve(
        current_density=None,
        project_source=False,
        extra_bilinear_terms=bil,
        extra_linear_terms=load_of(first),
    ).e_complex
    t_factor = time.perf_counter() - t0
    t0 = time.perf_counter()
    for key in extra_loads["order"]:
        if key == first:
            continue
        fields[key] = solver.solve_with_held_factorization(
            extra_linear_terms=load_of(key)
        ).e_complex
    t_held = time.perf_counter() - t0
    _free(solver)

    currents = {}
    v_src = {}
    for d in PORT_IDS:
        if d not in fields:
            continue
        sh = _sheets(specs, d)
        currents[d] = {
            s.port_id: sheet_terminal_current(msh, tags_f, s, fields[d], comm) for s in sh
        }
        v_src[d] = complex(next(s for s in sh if s.port_id == d).source_voltage_v)
    y = np.array(
        [[currents[j][i] / v_src[j] for j in PORT_IDS] for i in PORT_IDS], dtype=np.complex128
    )
    return fields, y, t_factor, t_held


def _symmetry(fx):
    sweep = fx["sweep"]
    msh = fx["mesh"]
    tags_f = sweep["facet_tags"]
    problem = fx["p_with"]
    solver = TimeHarmonicSolver(problem, degree=int(sweep["degree"]))
    v_space = solver.function_space()
    fx.setdefault("space_degrees", []).append(_space_degree(v_space))
    sigma_f, eps_f = build_material_fields(
        msh,
        problem.material,
        cell_tags=problem.cell_tags,
        material_map=problem.material_map,
        phantom_material=problem.phantom_material,
        phantom_tag=problem.phantom_tag,
    )
    a = solver.bilinear_form(sigma_f, eps_f)
    trial, test = ufl.TrialFunction(v_space), ufl.TestFunction(v_space)
    for s in _sheets(sweep["specs"]):
        a = a + lumped_port_bilinear_term(msh, tags_f, s, trial, test, omega_rad_per_s=OMEGA)
    mat = assemble_matrix(fem.form(a))
    mat.assemble()
    a_norm = mat.norm(PETSc.NormType.FROBENIUS)
    mat_t = mat.copy()
    mat_t.transpose()
    mat_t.axpy(-1.0, mat, structure=PETSc.Mat.Structure.DIFFERENT_NONZERO_PATTERN)
    sym_norm = mat_t.norm(PETSc.NormType.FROBENIUS)
    ndofs = v_space.dofmap.index_map.size_global
    mat.destroy()
    mat_t.destroy()
    del solver
    gc.collect()
    return sym_norm, a_norm, ndofs


def _run_step1_identities(fx, degree=1):
    """Step 1's identities body, shared by step 2 (`degree=2`).

    ``degree`` must equal the fixture's sweep degree (``FEM_EM_WF9_DEGREE``);
    every solver in the body takes ``fx["sweep"]["degree"]``.
    """
    comm = fx["comm"]
    msh = fx["mesh"]
    assert int(fx["sweep"]["degree"]) == int(degree), (fx["sweep"]["degree"], degree)
    fx["space_degrees"] = []
    t_start = time.perf_counter()
    dx = ufl.dx(domain=msh)
    zhat = ufl.as_vector((0.0, 0.0, 1.0))
    chi_b, chi_lead, chi_k = fx["chi_b"], fx["chi_lead"], fx["chi_k"]
    ds_ = fx["delta_sigma"]
    times = {}

    v_ball = _scalar(comm, chi_b * dx).real

    def q(e):
        return _scalar(comm, chi_b * e[2] * dx) / v_ball

    ell = lambda v: ufl.inner((chi_b / v_ball) * zhat, v) * ufl.dx  # noqa: E731

    # ---- without arm: factor on P1, held P2..P4 and the tip load (w0) ---------
    t0 = time.perf_counter()
    e_o, y_o, tf_o, th_o = _arm(
        fx,
        fx["p_without"],
        first="P1",
        extra_loads={"ell": ell, "order": PORT_IDS + ("ell",)},
    )
    w0 = e_o.pop("ell")
    times["without_factor"] = tf_o
    times["without_held4"] = th_o

    # ---- with arm: factor on the tip load (w), held P1..P4 --------------------
    e_w, y_w, tf_w, th_w = _arm(
        fx,
        fx["p_with"],
        first="ell",
        extra_loads={"ell": ell, "order": ("ell",) + PORT_IDS},
    )
    w = e_w.pop("ell")
    times["with_factor"] = tf_w
    times["with_held4"] = th_w

    # ---- (H) fresh P1 on its own solver ----------------------------------------
    t0 = time.perf_counter()
    fresh = _solve_driven(dict(fx["sweep"], problem=fx["p_with"]), "P1")
    q_fresh = q(fresh["fields"].e_complex)
    fx["space_degrees"].append(_space_degree(fresh["fields"].e_complex.function_space))
    del fresh
    gc.collect()
    times["fresh_H"] = time.perf_counter() - t0

    # ---- (A1) ------------------------------------------------------------------
    t0 = time.perf_counter()
    sym_norm, a_norm, ndofs = _symmetry(fx)
    times["A1"] = time.perf_counter() - t0

    # ---- read-backs (components / ufl.dot only: no conjugation) ---------------
    tf_k = np.array([LOAD * ds_ * _scalar(comm, chi_k[k] * w[2] * dx) for k in range(N_SLABS)])
    v_k = np.array([_scalar(comm, chi_k[k] * dx).real for k in range(N_SLABS)])

    res = {}
    for d in PORT_IDS:
        eo, ew = e_o[d], e_w[d]
        q_with = q(ew)
        q_without = q(eo)
        dq = q_with - q_without
        t3 = LOAD * ds_ * _scalar(comm, chi_lead * ufl.dot(w, eo) * dx)
        born = LOAD * ds_ * _scalar(comm, chi_lead * ufl.dot(w0, eo) * dx)
        c1 = LOAD * ds_ * _scalar(comm, chi_lead * ufl.dot(ufl.conj(w), eo) * dx)
        c2 = LOAD * ds_ * _scalar(comm, chi_lead * ufl.dot(w, ew) * dx)
        transverse = LOAD * ds_ * _scalar(comm, chi_lead * (w[0] * eo[0] + w[1] * eo[1]) * dx)
        ez_int = np.array([_scalar(comm, chi_k[k] * eo[2] * dx) for k in range(N_SLABS)])
        along = np.array(
            [LOAD * ds_ * _scalar(comm, chi_k[k] * w[2] * eo[2] * dx) for k in range(N_SLABS)]
        )
        folds = {}
        for n_seg in FOLD_NS:
            per = N_SLABS // n_seg
            qn = 0.0 + 0.0j
            al = 0.0 + 0.0j
            for j in range(n_seg):
                sl = slice(j * per, (j + 1) * per)
                qn += np.sum(tf_k[sl]) * (np.sum(ez_int[sl]) / np.sum(v_k[sl]))
                al += np.sum(along[sl])
            folds[n_seg] = (qn, al - qn)
        res[d] = dict(
            q_with=q_with,
            q_without=q_without,
            dq=dq,
            t3=t3,
            born=born,
            c1=c1,
            c2=c2,
            transverse=transverse,
            folds=folds,
            q_held=q_with,
        )

    rss = _rss_gib(comm)
    elapsed = time.perf_counter() - t_start

    # Trap 2 (item 276): degree-2 exchange files and partner widths are
    # distinct, so a degree-2 run never reads step 1's degree-1 files.
    if int(degree) == 1:
        tf_stem, partners = "wf9_step1_tf_n", (2, 4)
    else:
        tf_stem, partners = f"wf9_step2_tf_p{int(degree)}_n", (8, 12)
    tf_agreement = None
    if comm.rank == 0:
        os.makedirs(TF_EXCHANGE_DIR, exist_ok=True)
        np.save(os.path.join(TF_EXCHANGE_DIR, f"{tf_stem}{comm.size}.npy"), tf_k)
        for other in partners:
            if other == comm.size:
                continue
            path = os.path.join(TF_EXCHANGE_DIR, f"{tf_stem}{other}.npy")
            if os.path.exists(path):
                ref = np.load(path)
                tf_agreement = (other, float(np.max(np.abs(tf_k - ref)) / np.max(np.abs(ref))))

    recip = {
        "with": float(np.linalg.norm(y_w - y_w.T) / np.linalg.norm(y_w)),
        "without": float(np.linalg.norm(y_o - y_o.T) / np.linalg.norm(y_o)),
    }
    h_dev = float(abs(q_fresh - res["P1"]["q_held"]) / abs(res["P1"]["q_held"]))

    _say(
        f"\n[WF-9 step1 identities] ranks={comm.size} degree={int(degree)} dofs={ndofs} "
        f"V_ball={v_ball:.6e} m^3; space degrees per factor {fx['space_degrees']}"
    )
    _say(
        f"  (A1) ||A-A^T||_F/||A||_F = {sym_norm / a_norm:.3e} (bound {SYMMETRY_RTOL:g}; ||A||_F={a_norm:.6e})"
    )
    _say(
        f"  (i') ||Y-Y^T||/||Y||: with {recip['with']:.3e}, without {recip['without']:.3e} "
        f"(band {RECIPROCITY_BAND:g})"
    )
    _say(
        f"  (H) |q_fresh - q_held|/|q_held| (P1, with) = {h_dev:.3e} (band {DISCRETE_IDENTITY_RTOL:g})"
    )
    for d in PORT_IDS:
        r = res[d]
        dq = r["dq"]
        ratio = abs(r["q_with"]) / abs(r["q_without"])
        pf = abs(r["folds"][32][0] - dq) / abs(dq)
        bf = abs(r["born"] - dq) / abs(dq)
        _say(f"  --- {d} ---")
        _say(f"    q_without={r['q_without']:.10e} q_with={r['q_with']:.10e}")
        _say(
            f"    |q_with/q_without| = {ratio:.6f} (predicted [{PRED_RATIO[d][0]:g}, "
            f"{PRED_RATIO[d][1]:g}], {_in(ratio, PRED_RATIO[d])})"
        )
        _say(
            f"    (T3) |fold(w,E_without) - dq|/|dq| = {abs(r['t3'] - dq) / abs(dq):.3e}; "
            f"q_without + fold = {r['q_without'] + r['t3']:.10e} vs q_with"
        )
        _say(
            f"    (PF) |Q32 - dq|/|dq| = {pf:.6e} (band {PATH_FOLD_BAND:g}; predicted "
            f"[{PRED_PF[d][0]:g}, {PRED_PF[d][1]:g}], {_in(pf, PRED_PF[d])})"
        )
        _say(
            f"    (BF) |Born - dq|/|dq| = {bf:.6%} (must exceed {PATH_FOLD_BAND:.0%}; predicted "
            f"[{PRED_BF[0]:.0%}, {PRED_BF[1]:.0%}], {_in(bf, PRED_BF)})"
        )
        c1v = abs(r["c1"] - dq) / abs(dq)
        c2v = abs(r["c2"] - dq) / abs(dq)
        _say(
            f"    (C1) Hermitian read-back: {c1v:.6%} (predicted [{PRED_C1[0]:.0%}, {PRED_C1[1]:.0%}], {_in(c1v, PRED_C1)})"
        )
        _say(
            f"    (C2) total-field fold (E_with): {c2v:.6%} (predicted [{PRED_C2[0]:.0%}, {PRED_C2[1]:.0%}], {_in(c2v, PRED_C2)})"
        )
        _say(f"    error split: transverse/|dq| = {abs(r['transverse']) / abs(dq):.6e}")
        for n_seg, (qn, cov) in r["folds"].items():
            _say(
                f"    N={n_seg:2d} Q_N={qn:.10e} |Q_N-dq|/|dq|={abs(qn - dq) / abs(dq):.6e} "
                f"within-slab cov/|dq|={abs(cov) / abs(dq):.6e}"
            )
    _say("  TF_k profile (k, z_mid mm, V_k m^3, |TF|, arg deg):")
    for k in range(N_SLABS):
        _say(
            f"    {k:2d} {(SLAB_Z0 + (k + 0.5) * SLAB_DZ) * 1e3:6.1f} {v_k[k]:.4e} "
            f"{abs(tf_k[k]):.6e} {np.degrees(np.angle(tf_k[k])):9.4f}"
        )
    _say("  times: " + ", ".join(f"{k}={v:.2f}s" for k, v in times.items()))
    rss_pred = (
        "predicted <= 3.5 at -n 2"
        if int(degree) == 1
        else f"predicted [{PRED_RSS_P2_N8[0]:g}, {PRED_RSS_P2_N8[1]:g}] at -n 8"
    )
    _say(f"  summed max RSS over ranks = {rss:.3f} GiB ({rss_pred}); elapsed = {elapsed:.2f} s")
    if tf_agreement is None:
        _say("  W1/W2 TF agreement: no other-width TF file yet")
    else:
        _say(
            f"  W1/W2 TF agreement vs n={tf_agreement[0]}: max rel = {tf_agreement[1]:.3e} (predicted <= 1e-10)"
        )

    assert sym_norm <= SYMMETRY_RTOL * a_norm
    assert recip["with"] <= RECIPROCITY_BAND, recip
    assert recip["without"] <= RECIPROCITY_BAND, recip
    assert h_dev <= DISCRETE_IDENTITY_RTOL, h_dev
    for d in PORT_IDS:
        r = res[d]
        dq = r["dq"]
        assert abs(r["t3"] - dq) <= DISCRETE_IDENTITY_RTOL * abs(dq), f"T3 {d}"
        assert abs(r["folds"][32][0] - dq) <= PATH_FOLD_BAND * abs(dq), f"PF {d}"
        assert abs(r["born"] - dq) > PATH_FOLD_BAND * abs(dq), f"BF {d} does not miss"
        assert abs(r["c1"] - dq) > DISCRETE_IDENTITY_RTOL * abs(dq), f"C1 {d} not wired"
        assert abs(r["c2"] - dq) > DISCRETE_IDENTITY_RTOL * abs(dq), f"C2 {d} not wired"

    return {
        "q": q,
        "v_ball": v_ball,
        "e_without": e_o,
        "e_with": e_w,
        "w": w,
        "w0": w0,
        "res": res,
        "tf_k": tf_k,
        "v_k": v_k,
        "y_with": y_w,
        "y_without": y_o,
        "recip": recip,
        "sym_ratio": sym_norm / a_norm,
        "ndofs": ndofs,
        "h_dev": h_dev,
        "times": times,
        "rss": rss,
        "elapsed": elapsed,
        "tf_agreement": tf_agreement,
        "degree": int(degree),
        "space_degrees": list(fx["space_degrees"]),
    }


@complex_only
@pytest.mark.skipif(_degree() != 1, reason=f"{DEGREE_ENV} != 1 (step 1 is degree 1)")
def test_wf9_step1_identities(_case):
    _run_step1_identities(_case, degree=1)


@complex_only
@pytest.mark.skipif(_degree() != 2, reason=f"{DEGREE_ENV} != 2 (step 2 opt-in)")
def test_wf9_step2_degree2(_case):
    """`WF-9` step 2 (§9 item 276): step 1's identities at degree 2, plus
    (Deg) and the printed order move against `STEP1_Q_RECORD`."""
    fx = _case
    comm = fx["comm"]
    out = _run_step1_identities(fx, degree=2)
    res = out["res"]

    _say(f"\n[WF-9 step2 degree 2] ranks={comm.size}")
    _say(
        f"  (Deg) space degrees per factor {out['space_degrees']}; global dofs {out['ndofs']} "
        f"vs 4 x {DEG1_DOFS_RECORD} = {4 * DEG1_DOFS_RECORD} (ratio {out['ndofs'] / DEG1_DOFS_RECORD:.3f})"
    )
    _say("  order move (printed only; vs step 1 degree-1 record):")
    for d in PORT_IDS:
        qo1, qw1 = STEP1_Q_RECORD[d]
        dq1 = qw1 - qo1
        dq2 = res[d]["dq"]
        r = dq2 / dq1
        mv = abs(r - 1.0)
        _say(
            f"    {d}: |dq2/dq1 - 1| = {mv:.6e} (predicted [{PRED_ORDER_MOVE[0]:g}, "
            f"{PRED_ORDER_MOVE[1]:g}], {_in(mv, PRED_ORDER_MOVE)}); arg(dq2/dq1) = "
            f"{np.degrees(np.angle(r)):.6f} deg; |q_with2/q_with1 - 1| = "
            f"{abs(res[d]['q_with'] / qw1 - 1.0):.6e}"
        )
    if comm.rank == 0:
        path = os.path.join(TF_EXCHANGE_DIR, "wf9_step1_tf_n2.npy")
        if os.path.exists(path):
            tf1 = np.load(path)
            move = float(np.max(np.abs(out["tf_k"] - tf1)) / np.max(np.abs(tf1)))
            print(f"  TF order move max_k |TF2-TF1|/max|TF1| = {move:.6e}", flush=True)
        else:
            print("  TF order move: wf9_step1_tf_n2.npy absent", flush=True)

    assert all(p == 2 for p in out["space_degrees"]), out["space_degrees"]
    assert len(out["space_degrees"]) == 4, out["space_degrees"]
    assert out["ndofs"] > 4 * DEG1_DOFS_RECORD, out["ndofs"]
