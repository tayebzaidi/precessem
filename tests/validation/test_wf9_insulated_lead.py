"""`WF-9` step 3 (§9 item 317, chain I18b): the insulated lead's transfer function.

`GEO-36` step 5's coated build (item 316): the kept rod (tag
`BIRDCAGE_IMPLANT_CELL_TAG` 403, r = 1 mm) inside a shaft sleeve (tag
`BIRDCAGE_IMPLANT_COAT_CELL_TAG` 404, b = 1.5 mm) in the gel (tag 3) of the
coil-driven F-small birdcage at 64 MHz, degree 1.  Three material states on
**one** mesh, switched by material only:

* without: 403 and 404 gel (the lead-free ISO 10974 background);
* bare:    403 σ = 800 S/m with the gel's εᵣ, 404 gel (step 1's lead);
* coated:  403 as bare, 404 εᵣ 3, σ 0, μᵣ 1 (`IMP-5.md:29`).

The fold weight is the per-cell complex contrast read off the solver's own
operator: the mass term of ``TimeHarmonicSolver.bilinear_form`` is
``−k₀² ε_c E·v`` with ``ε_c = εᵣ − jσ/(ωε₀)``, so A_state − A_without =
−k₀² Δε_c M and the adjoint identity Δq = −wᵀ ΔA e_without reads
Δq = Σ c·w·E_without with c = k₀² Δε_c (= LOAD·Δσ + ω²μ₀ε₀Δεᵣ).  (R1)
asserts that on 403 (εᵣ ≡ gel) c reduces to step 1's LOAD·Δσ cell by cell.

Asserted (bands imported; brief ``docs/planning/chunks/WF-9.md`` § item 317):
(A0) census, (M0) read-back, (R1), (A1) A = Aᵀ per state, (T3c) the adjoint
identity per drive on both arms, (H) held vs fresh (coated P1), (PF) N = 32
path fold within `PATH_FOLD_BAND`, (BF) Born fold misses; controls C3c / C1 /
C2 asserted only as wired.  Printed only: TF profiles, coated/bare ratios,
tip ratios, coat share of the fold, times, RSS.

Not covered: a measured or ISO-calibrated TF, heating, SAR, a real lead,
degree 2, h-convergence, the `IMP-5` series.
"""

from __future__ import annotations

import gc
import os
import time
from dataclasses import replace

import numpy as np
import pytest
import ufl
from dolfinx import fem
from mpi4py import MPI

from fem_em_solver.core import HomogeneousMaterial
from fem_em_solver.core.time_harmonic import build_material_fields
from fem_em_solver.io.mesh import (
    BIRDCAGE_IMPLANT_CELL_TAG,
    BIRDCAGE_IMPLANT_COAT_CELL_TAG,
)
from fem_em_solver.utils.constants import EPSILON_0, MU_0
from tests.complex_mode import complex_only
from tests.mesh.test_birdcage_implant_rod import (
    KEPT_ROD_CELL_GUARD,
    _build_rod,
    _cell_values_with_ghosts,
)
from tests.mesh.test_birdcage_implant_rod_coat import COATED_ROD
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
from tests.validation.test_wf9_lead_in_coil_phantom import (
    FOLD_NS,
    LOAD,
    N_SLABS,
    OMEGA,
    PATH_FOLD_BAND,
    PORT_IDS,
    PRED_RATIO,
    SLAB_DZ,
    SLAB_Z0,
    STEP1_Q_RECORD,
    _arm,
    _material_minmax,
    _scalar,
    _symmetry,
)
from tests.validation.test_wf9_transfer_function_identity import (
    BALL_RADIUS,
    SYMMETRY_RTOL,
)

#: Item 316's printed coated-build record
#: (`20261008T141656Z_GEO-36-step5-W2.log:5494`, `-n 4`: 190 747 cells,
#: tags 3 / 403 / 404 = 41 602 / 15 035 / 18 495).
COATED_CELL_RECORD = 190_747

#: `IMP-5.md:29` order-of-magnitude silicone sleeve.
COAT_MAT = HomogeneousMaterial(sigma=0.0, epsilon_r=3.0, mu_r=1.0)

#: (R1) "exactly": k₀²·(−jΔσ/(ωε₀)) vs −jωμ₀Δσ differ only by the rounding
#: of k₀²/(ωε₀) against ωμ₀ (a few ulp); 1e-13 relative is that floor.
R1_RTOL = 1e-13

K0_SQ = OMEGA * OMEGA * MU_0 * EPSILON_0
TF_EXCHANGE_DIR = "/workspace/logs"

# Printed-only predictions (brief § item 317; never asserted).
PRED_C1 = (0.10, 2.00)
PRED_C2 = (0.80, 1.05)


def _say(msg):
    if MPI.COMM_WORLD.rank == 0:
        print(msg, flush=True)


def _in(x, band):
    return "IN" if band[0] <= x <= band[1] else "OUT"


def _build_case():
    comm = MPI.COMM_WORLD
    t0 = time.perf_counter()
    built = _build_rod(implant_rod=dict(COATED_ROD))
    msh = built[0]
    tdim = msh.topology.dim
    n_cells = int(msh.topology.index_map(tdim).size_global)
    assert n_cells <= KEPT_ROD_CELL_GUARD, (n_cells, KEPT_ROD_CELL_GUARD)
    t_mesh = time.perf_counter() - t0

    case = build_four_port_sweep(
        frequency_hz=STEP2D_FREQUENCY_HZ, built=built, build_only=True, degree=1
    )
    problem = case["problem"]
    mmap = dict(problem.material_map)
    gel = mmap[PHANTOM_CELL_TAG]
    rod_mat = HomogeneousMaterial(sigma=SIGMA_WIRE_S_PER_M, epsilon_r=gel.epsilon_r, mu_r=gel.mu_r)
    states = {
        "without": {BIRDCAGE_IMPLANT_CELL_TAG: gel, BIRDCAGE_IMPLANT_COAT_CELL_TAG: gel},
        "bare": {BIRDCAGE_IMPLANT_CELL_TAG: rod_mat, BIRDCAGE_IMPLANT_COAT_CELL_TAG: gel},
        "coated": {BIRDCAGE_IMPLANT_CELL_TAG: rod_mat, BIRDCAGE_IMPLANT_COAT_CELL_TAG: COAT_MAT},
    }
    problems = {k: replace(problem, material_map={**mmap, **v}) for k, v in states.items()}
    for k, p in problems.items():
        # Every implant tag set explicitly in every state (never a default).
        assert BIRDCAGE_IMPLANT_CELL_TAG in p.material_map, k
        assert BIRDCAGE_IMPLANT_COAT_CELL_TAG in p.material_map, k

    imap = msh.topology.index_map(tdim)
    n_owned = imap.size_local
    n_local = n_owned + imap.num_ghosts
    cells = np.arange(n_local, dtype=np.int32)
    mid = msh.geometry.x[msh.geometry.dofmap[cells]].mean(axis=1)
    tag = _cell_values_with_ghosts(msh, built[1])
    rod = tag == BIRDCAGE_IMPLANT_CELL_TAG
    coat = tag == BIRDCAGE_IMPLANT_COAT_CELL_TAG
    lead = rod | coat
    ball = (tag == PHANTOM_CELL_TAG) & (
        np.linalg.norm(mid - np.asarray(TIPS["top"]), axis=1) < BALL_RADIUS
    )
    slab_idx = np.floor((mid[:, 2] - SLAB_Z0) / SLAB_DZ).astype(int)
    slabs_lead = [lead & (slab_idx == k) for k in range(N_SLABS)]
    slabs_rod = [rod & (slab_idx == k) for k in range(N_SLABS)]

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
        "q0": q0,
        "flags": {
            "rod": rod,
            "coat": coat,
            "lead": lead,
            "ball": ball,
            "slabs_lead": slabs_lead,
            "slabs_rod": slabs_rod,
        },
        "owned_count": owned_count,
        "chi_rod": mask(rod),
        "chi_coat": mask(coat),
        "chi_lead": mask(lead),
        "chi_b": mask(ball),
        "chi_k_lead": [mask(f) for f in slabs_lead],
        "chi_k_rod": [mask(f) for f in slabs_rod],
        "gel": gel,
        "problems": problems,
        "t_mesh": t_mesh,
        # `_arm` / `_symmetry` read these keys.
        "p_with": problems["coated"],
        "p_without": problems["without"],
    }


@pytest.fixture(scope="module")
def _case():
    fx = _build_case()
    yield fx
    gc.collect()


def _minmax_on(fx, problem, flag):
    view = {
        "comm": fx["comm"],
        "mesh": fx["mesh"],
        "n_owned": fx["n_owned"],
        "flags": {"lead": flag},
    }
    return _material_minmax(view, problem)


def _eps_c_cells(fx, problem):
    """ε_c = εᵣ − jσ/(ωε₀) per local cell, as ``bilinear_form`` builds it."""
    sigma_f, eps_f = build_material_fields(
        fx["mesh"],
        problem.material,
        cell_tags=problem.cell_tags,
        material_map=problem.material_map,
        phantom_material=problem.phantom_material,
        phantom_tag=problem.phantom_tag,
    )
    cells = fx["cells"]
    sig = sigma_f.x.array[sigma_f.function_space.dofmap.list[:, 0][cells]]
    eps = eps_f.x.array[eps_f.function_space.dofmap.list[:, 0][cells]]
    return eps - 1j * sig / (OMEGA * EPSILON_0)


def _contrast(fx, state):
    """DG0 c = k₀² (ε_c,state − ε_c,without), from the solver's own operator."""
    vals = K0_SQ * (
        _eps_c_cells(fx, fx["problems"][state]) - _eps_c_cells(fx, fx["problems"]["without"])
    )
    c = fem.Function(fx["q0"])
    c.x.array[:] = 0.0
    n_owned = fx["n_owned"]
    c.x.array[fx["dof_of_cell"][fx["cells"][:n_owned]]] = vals[:n_owned]
    c.x.scatter_forward()
    return c, vals


@complex_only
def test_wf9_step3_census(_case):
    fx = _case
    comm = fx["comm"]
    fl = fx["flags"]
    cnt = fx["owned_count"]
    n_rod, n_coat, n_ball = cnt(fl["rod"]), cnt(fl["coat"]), cnt(fl["ball"])
    slab_counts = [cnt(f) for f in fl["slabs_lead"]]
    stacked = np.sum(np.vstack(fl["slabs_lead"]).astype(int), axis=0)
    partition_ok = comm.allreduce(
        bool(np.array_equal(stacked, fl["lead"].astype(int))), op=MPI.LAND
    )
    disjoint = comm.allreduce(bool(not np.any(fl["ball"] & fl["lead"])), op=MPI.LAND)
    rel_cells = abs(fx["n_cells"] / COATED_CELL_RECORD - 1.0)
    gel = fx["gel"]
    mm = {
        (s, t): _minmax_on(fx, fx["problems"][s], fl[t])
        for s in ("without", "bare", "coated")
        for t in ("rod", "coat")
    }
    expect = {
        ("without", "rod"): (gel.sigma, gel.epsilon_r),
        ("without", "coat"): (gel.sigma, gel.epsilon_r),
        ("bare", "rod"): (SIGMA_WIRE_S_PER_M, gel.epsilon_r),
        ("bare", "coat"): (gel.sigma, gel.epsilon_r),
        ("coated", "rod"): (SIGMA_WIRE_S_PER_M, gel.epsilon_r),
        ("coated", "coat"): (COAT_MAT.sigma, COAT_MAT.epsilon_r),
    }

    _say(f"\n[WF-9 step3 census] ranks={comm.size}")
    _say(
        f"  (A0) cells={fx['n_cells']} vs COATED_CELL_RECORD {COATED_CELL_RECORD} "
        f"rel {rel_cells:.3e} (band {CELL_COUNT_BAND:g}); mesh+guard {fx['t_mesh']:.1f} s"
    )
    _say(f"  403={n_rod} 404={n_coat} ball={n_ball} lead slabs={slab_counts}")
    _say(f"  slab partition of rod∪coat: {partition_ok}; ball/lead disjoint: {disjoint}")
    for key, v in mm.items():
        _say(
            f"  (M0) {key[0]:7s} on {key[1]}: sigma [{v['sigma'][0]:.6g}, {v['sigma'][1]:.6g}] "
            f"eps_r [{v['eps'][0]:.6g}, {v['eps'][1]:.6g}]; expected {expect[key]}"
        )

    assert rel_cells <= CELL_COUNT_BAND
    assert n_rod > 0 and n_coat > 0 and n_ball > 0
    assert all(c > 0 for c in slab_counts), slab_counts
    assert partition_ok
    assert disjoint
    for key, (s_exp, e_exp) in expect.items():
        assert mm[key]["sigma"] == (s_exp, s_exp), (key, mm[key])
        assert mm[key]["eps"] == (e_exp, e_exp), (key, mm[key])


def _run_step3_identities(fx):
    comm = fx["comm"]
    msh = fx["mesh"]
    fl = fx["flags"]
    n_owned = fx["n_owned"]
    fx["space_degrees"] = []
    t_start = time.perf_counter()
    dx = ufl.dx(domain=msh)
    zhat = ufl.as_vector((0.0, 0.0, 1.0))
    chi_b = fx["chi_b"]
    times = {}

    # ---- contrasts from the operator; (R1) ------------------------------------
    c_fn, c_vals, r1 = {}, {}, {}
    ds_ = float(SIGMA_WIRE_S_PER_M - fx["gel"].sigma)
    ref = LOAD * ds_
    for arm in ("bare", "coated"):
        c_fn[arm], c_vals[arm] = _contrast(fx, arm)
        own = c_vals[arm][:n_owned]
        rod_o = fl["rod"][:n_owned]
        off_o = ~fl["lead"][:n_owned]
        dev = float(np.max(np.abs(own[rod_o] - ref))) if np.any(rod_o) else 0.0
        off = float(np.max(np.abs(own[off_o]))) if np.any(off_o) else 0.0
        r1[arm] = (comm.allreduce(dev, op=MPI.MAX) / abs(ref), comm.allreduce(off, op=MPI.MAX))

    def tag_mean(arm, flag):
        own = c_vals[arm][:n_owned][flag[:n_owned]]
        s = comm.allreduce(complex(np.sum(own)), op=MPI.SUM)
        n = comm.allreduce(int(own.size), op=MPI.SUM)
        return s / max(n, 1)

    v_ball = _scalar(comm, chi_b * dx).real

    def q(e):
        return _scalar(comm, chi_b * e[2] * dx) / v_ball

    ell = lambda v: ufl.inner((chi_b / v_ball) * zhat, v) * ufl.dx  # noqa: E731

    # ---- three arms; one factorised solver alive at a time (`_arm` frees) ------
    e_o, y_o, tf_, th_ = _arm(
        fx,
        fx["problems"]["without"],
        first="P1",
        extra_loads={"ell": ell, "order": PORT_IDS + ("ell",)},
    )
    w0 = e_o.pop("ell")
    times["without_factor"], times["without_held5"] = tf_, th_
    arms = {}
    for arm in ("bare", "coated"):
        e_w, y_w, tf_, th_ = _arm(
            fx,
            fx["problems"][arm],
            first="ell",
            extra_loads={"ell": ell, "order": ("ell",) + PORT_IDS},
        )
        w = e_w.pop("ell")
        arms[arm] = {"e": e_w, "y": y_w, "w": w}
        times[f"{arm}_factor"], times[f"{arm}_held4"] = tf_, th_

    # ---- (H) fresh coated P1 ---------------------------------------------------
    t0 = time.perf_counter()
    fresh = _solve_driven(dict(fx["sweep"], problem=fx["problems"]["coated"]), "P1")
    q_fresh = q(fresh["fields"].e_complex)
    del fresh
    gc.collect()
    times["fresh_H"] = time.perf_counter() - t0

    # ---- (A1) per state --------------------------------------------------------
    t0 = time.perf_counter()
    sym = {}
    for s, p in fx["problems"].items():
        sn, an, ndofs = _symmetry(dict(fx, p_with=p))
        sym[s] = sn / an
    times["A1x3"] = time.perf_counter() - t0

    # ---- read-backs (components / ufl.dot only) --------------------------------
    res = {}
    tf_prof = {}
    for arm in ("bare", "coated"):
        c = c_fn[arm]
        w = arms[arm]["w"]
        chi_k = fx["chi_k_lead"] if arm == "coated" else fx["chi_k_rod"]
        chi_l = fx["chi_lead"] if arm == "coated" else fx["chi_rod"]
        tf_k = np.array([_scalar(comm, chi_k[k] * c * w[2] * dx) for k in range(N_SLABS)])
        v_k = np.array([_scalar(comm, chi_k[k] * dx).real for k in range(N_SLABS)])
        tf_prof[arm] = (tf_k, v_k)
        res[arm] = {}
        for d in PORT_IDS:
            eo, ew = e_o[d], arms[arm]["e"][d]
            q_with, q_without = q(ew), q(eo)
            dq = q_with - q_without
            t3 = _scalar(comm, chi_l * c * ufl.dot(w, eo) * dx)
            born = _scalar(comm, chi_l * c * ufl.dot(w0, eo) * dx)
            c1 = _scalar(comm, chi_l * c * ufl.dot(ufl.conj(w), eo) * dx)
            c2 = _scalar(comm, chi_l * c * ufl.dot(w, ew) * dx)
            rod_only = _scalar(comm, fx["chi_rod"] * c * ufl.dot(w, eo) * dx)
            coat_part = _scalar(comm, fx["chi_coat"] * c * ufl.dot(w, eo) * dx)
            ez_int = np.array([_scalar(comm, chi_k[k] * eo[2] * dx) for k in range(N_SLABS)])
            folds = {}
            for n_seg in FOLD_NS:
                per = N_SLABS // n_seg
                qn = 0.0 + 0.0j
                for j in range(n_seg):
                    sl = slice(j * per, (j + 1) * per)
                    if np.sum(v_k[sl]) > 0.0:
                        qn += np.sum(tf_k[sl]) * (np.sum(ez_int[sl]) / np.sum(v_k[sl]))
                folds[n_seg] = qn
            res[arm][d] = dict(
                q_with=q_with,
                q_without=q_without,
                dq=dq,
                t3=t3,
                born=born,
                c1=c1,
                c2=c2,
                rod_only=rod_only,
                coat_part=coat_part,
                folds=folds,
            )

    rss_kb = float(__import__("resource").getrusage(__import__("resource").RUSAGE_SELF).ru_maxrss)
    rss = comm.allreduce(rss_kb, op=MPI.SUM) / 1024.0**2
    elapsed = time.perf_counter() - t_start

    tf_agreement = None
    if comm.rank == 0:
        os.makedirs(TF_EXCHANGE_DIR, exist_ok=True)
        stack = np.vstack([tf_prof["bare"][0], tf_prof["coated"][0]])
        np.save(os.path.join(TF_EXCHANGE_DIR, f"wf9_step3_tf_n{comm.size}.npy"), stack)
        for other in (2, 4):
            if other == comm.size:
                continue
            path = os.path.join(TF_EXCHANGE_DIR, f"wf9_step3_tf_n{other}.npy")
            if os.path.exists(path):
                ref_tf = np.load(path)
                tf_agreement = (
                    other,
                    float(np.max(np.abs(stack - ref_tf)) / np.max(np.abs(ref_tf))),
                )

    q_held = res["coated"]["P1"]["q_with"]
    h_dev = float(abs(q_fresh - q_held) / abs(q_held))
    recip = {
        a: float(np.linalg.norm(arms[a]["y"] - arms[a]["y"].T) / np.linalg.norm(arms[a]["y"]))
        for a in arms
    }
    recip["without"] = float(np.linalg.norm(y_o - y_o.T) / np.linalg.norm(y_o))

    _say(
        f"\n[WF-9 step3 identities] ranks={comm.size} dofs={ndofs} V_ball={v_ball:.6e} m^3; "
        f"space degrees per factor {fx['space_degrees']}"
    )
    for arm in ("bare", "coated"):
        _say(
            f"  (R1) {arm}: max|c-LOAD*dsigma|/|LOAD*dsigma| on 403 = {r1[arm][0]:.3e} "
            f"(rtol {R1_RTOL:g}); max|c| off lead = {r1[arm][1]:.3e}; "
            f"mean c on 403 = {tag_mean(arm, fl['rod']):.6e}, on 404 = "
            f"{tag_mean(arm, fl['coat']):.6e}"
        )
    _say(
        "  (A1) ||A-A^T||_F/||A||_F: "
        + ", ".join(f"{s} {v:.3e}" for s, v in sym.items())
        + f" (bound {SYMMETRY_RTOL:g})"
    )
    _say("  (i') ||Y-Y^T||/||Y|| (printed): " + ", ".join(f"{s} {v:.3e}" for s, v in recip.items()))
    _say(
        f"  (H) |q_fresh - q_held|/|q_held| (coated P1) = {h_dev:.3e} (band {DISCRETE_IDENTITY_RTOL:g})"
    )
    for arm in ("bare", "coated"):
        for d in PORT_IDS:
            r = res[arm][d]
            dq = r["dq"]
            ratio = abs(r["q_with"]) / abs(r["q_without"])
            t3v = abs(r["t3"] - dq) / abs(dq)
            pf = abs(r["folds"][32] - dq) / abs(dq)
            bf = abs(r["born"] - dq) / abs(dq)
            c3 = abs(r["rod_only"] - dq) / abs(dq)
            c1v = abs(r["c1"] - dq) / abs(dq)
            c2v = abs(r["c2"] - dq) / abs(dq)
            _say(f"  --- {arm} {d} ---")
            _say(f"    q_without={r['q_without']:.10e} q_with={r['q_with']:.10e} dq={dq:.10e}")
            _say(
                f"    |q_with/q_without| = {ratio:.6f} (step 1 PRED_RATIO "
                f"[{PRED_RATIO[d][0]:g}, {PRED_RATIO[d][1]:g}], {_in(ratio, PRED_RATIO[d])})"
            )
            if arm == "bare":
                qo1, qw1 = STEP1_Q_RECORD[d]
                _say(
                    f"    step-1 dq record (other mesh, printed only) = {qw1 - qo1:.10e}; "
                    f"|dq/dq_step1| = {abs(dq) / abs(qw1 - qo1):.6f}"
                )
            _say(
                f"    (T3c) |fold(c w, E_without) - dq|/|dq| = {t3v:.3e} (rtol {DISCRETE_IDENTITY_RTOL:g})"
            )
            _say(f"    (PF) |Q32 - dq|/|dq| = {pf:.6e} (band {PATH_FOLD_BAND:g})")
            _say(f"    (BF) |Born - dq|/|dq| = {bf:.6%} (must exceed {PATH_FOLD_BAND:.0%})")
            _say(
                f"    (C3c) rod-only contrast: {c3:.6%} (wired: > {DISCRETE_IDENTITY_RTOL:g}); "
                f"coat share |fold_404/fold| = {abs(r['coat_part']) / abs(r['t3']):.6f}"
            )
            _say(
                f"    (C1) Hermitian: {c1v:.6%} (predicted [{PRED_C1[0]:.0%}, {PRED_C1[1]:.0%}], "
                f"{_in(c1v, PRED_C1)}); (C2) total field: {c2v:.6%} (predicted "
                f"[{PRED_C2[0]:.0%}, {PRED_C2[1]:.0%}], {_in(c2v, PRED_C2)})"
            )
            for n_seg, qn in r["folds"].items():
                _say(f"    N={n_seg:2d} |Q_N-dq|/|dq|={abs(qn - dq) / abs(dq):.6e}")
    tfb, vb = tf_prof["bare"]
    tfc, vc = tf_prof["coated"]
    _say(
        "  TF profile (k, z_mid mm, V_bare, V_coated, |TF_bare|, |TF_coated|, |TF_c/TF_b|, "
        "arg_b deg, arg_c deg):"
    )
    for k in range(N_SLABS):
        rat = abs(tfc[k]) / abs(tfb[k]) if abs(tfb[k]) > 0 else float("nan")
        _say(
            f"    {k:2d} {(SLAB_Z0 + (k + 0.5) * SLAB_DZ) * 1e3:6.1f} {vb[k]:.4e} {vc[k]:.4e} "
            f"{abs(tfb[k]):.6e} {abs(tfc[k]):.6e} {rat:.6f} "
            f"{np.degrees(np.angle(tfb[k])):9.4f} {np.degrees(np.angle(tfc[k])):9.4f}"
        )
    shaft_tip = {}
    for arm, tfa in (("bare", tfb), ("coated", tfc)):
        shaft = np.mean(np.abs(tfa[4:28]))
        tips = 0.5 * (np.abs(tfa[0]) + np.abs(tfa[-1]))
        shaft_tip[arm] = shaft / tips
        _say(f"  shaft(|z|<12 mm mean)/tip(end slabs mean) |TF| {arm} = {shaft / tips:.6f}")
    _say("  times: " + ", ".join(f"{k}={v:.2f}s" for k, v in times.items()))
    _say(f"  summed max RSS over ranks = {rss:.3f} GiB; elapsed = {elapsed:.2f} s")
    if tf_agreement is None:
        _say("  W1/W2 TF agreement: no other-width TF file yet")
    else:
        _say(f"  W1/W2 TF agreement vs n={tf_agreement[0]}: max rel = {tf_agreement[1]:.3e}")

    for arm in ("bare", "coated"):
        assert r1[arm][0] <= R1_RTOL, (arm, r1[arm])
        assert r1[arm][1] == 0.0, (arm, r1[arm])
    for s, v in sym.items():
        assert v <= SYMMETRY_RTOL, (s, v)
    assert h_dev <= DISCRETE_IDENTITY_RTOL, h_dev
    for arm in ("bare", "coated"):
        for d in PORT_IDS:
            r = res[arm][d]
            dq = r["dq"]
            assert abs(r["t3"] - dq) <= DISCRETE_IDENTITY_RTOL * abs(dq), f"T3c {arm} {d}"
            assert abs(r["folds"][32] - dq) <= PATH_FOLD_BAND * abs(dq), f"PF {arm} {d}"
            assert abs(r["born"] - dq) > PATH_FOLD_BAND * abs(dq), f"BF {arm} {d} does not miss"
            assert abs(r["c1"] - dq) > DISCRETE_IDENTITY_RTOL * abs(dq), f"C1 {arm} {d}"
            assert abs(r["c2"] - dq) > DISCRETE_IDENTITY_RTOL * abs(dq), f"C2 {arm} {d}"
            if arm == "coated":
                assert abs(r["rod_only"] - dq) > DISCRETE_IDENTITY_RTOL * abs(dq), f"C3c {d}"

    return dict(
        res=res,
        tf_prof=tf_prof,
        r1=r1,
        sym=sym,
        recip=recip,
        h_dev=h_dev,
        ndofs=ndofs,
        v_ball=v_ball,
        shaft_tip=shaft_tip,
        times=times,
        rss=rss,
        elapsed=elapsed,
        tf_agreement=tf_agreement,
    )


@complex_only
def test_wf9_step3_identities(_case):
    _run_step3_identities(_case)
