"""Example (`EX-121`): point probes read straight off the gel-ball submesh.

`mat:13` reads the heated gel ball through ball-averaged probes and sets them
beside a closed-form point value; `mat:16` reads the tag-3 submesh only as
owned-cell integrals.  This one is the first example that reads a FEM field
**at points on a ``create_submesh`` mesh**: `MAT-10`'s Robin-cooled gel
sphere (19 496 cells) cut from the 74 020-cell sphere-in-box EM mesh and
heated by the FEM SAR.  A ``create_submesh`` mesh carries no original cell
index, so the owner rule takes ``parent=(parent_mesh, parent_cells)`` (`POST-10`)
and keys each submesh cell by its parent's original cell index.

Part A is `POST-10`'s own cube on the gate's fixture and band.  Part B is
`MAT-10`'s ball: the gate's four anchors (B) / (P) / (C) / (A) and controls
(a) / (b), then six probe points B0-B5 (centre, 3.5 mm, 6.5 mm off-axis,
9.5 mm, a mesh vertex, outside) read through ``owner=RULE, parent=...``.

**Import, never restate** (the `ANS-1` rule).  Every fixture helper and band
comes from ``tests/post/test_submesh_owner_rule.py`` or
``tests/validation/test_mat10_sar_to_heat_sphere.py`` by import; no pytest
fixture or ``test_*`` function is called.  No new band is set.

Asserted: Part A (S0), (K), (I), (N), (W), (C1)-(C4) at the POST-10 bands;
Part B (B), (P), (C), (A) at the MAT-10 bands; on the ball (I-q), (I-o), (O),
(N-ball), (W-ball) at band 0 / ``CLOSED_FORM_TOL``; (C2-ball) the missing
``parent`` raises on every rank; MAT-10 controls (a) / (b) miss ``A_BAND``.
Printed, not asserted: every pointwise theta against the closed form, the
2 mm ball probe beside the point value, the owner=None comparison, the B4
incident-cell spread (cV) and the ball-vs-point offset (cB).

Run it::

    ./scripts/run_examples.sh -e mat:21 -n 2 -t 300

Outputs (``examples/materials/paraview_output/``):
``materials_21_point_probes_on_gel_ball_submesh_combined.xdmf`` and
``materials_21_point_probes_on_gel_ball_submesh_probes.csv``.
"""

from __future__ import annotations

import csv
import sys
import time
from pathlib import Path

import numpy as np
import ufl
from dolfinx import default_scalar_type, fem, geometry
from mpi4py import MPI

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from fem_em_solver.io.paraview_utils import (  # noqa: E402
    adopt_host_ownership,
    write_xdmf_time_series,
)
from fem_em_solver.post.evaluation import (  # noqa: E402
    evaluate_vector_field_parallel,
    point_owner_cells,
)
from fem_em_solver.post.sar import uniform_sphere_sar_closed_form  # noqa: E402
from fem_em_solver.post.setup_figure import write_setup_figure  # noqa: E402
from fem_em_solver.thermal.probe import ProbeBall, ball_average  # noqa: E402
from fem_em_solver.thermal.steady import solve_steady_conduction  # noqa: E402
from fem_em_solver.utils.analytical import robin_sphere_steady_temperature  # noqa: E402
from fem_em_solver.utils.constants import EPSILON_0  # noqa: E402

from tests.post.test_submesh_owner_rule import (  # noqa: E402
    A_VEC,
    B_VEC,
    CLOSED_FORM_TOL,
    INSIDE,
    LABELS,
    OUTSIDE,
    POINTS,
    RULE,
    SUB_CELLS_GLOBAL,
    _bitwise_equal_across_ranks,
    _dg0,
    _fixture,
    _oracle,
    _original,
)
from tests.validation.test_mat10_sar_to_heat_sphere import (  # noqa: E402
    A_BAND,
    B_BAND,
    C_BAND,
    E0,
    EPSILON_R_SPHERE,
    G,
    H,
    K,
    OMEGA,
    P_BAND,
    R,
    RHO_KG_M3,
    ROBIN_TAG,
    SIGMA_HIGH,
    SPHERE_TAG,
    T_INF,
    _allsum,
    _cell_average_dg0,
    _interior_indicator,
    _interior_rel,
    _robin_tags,
    _solve_lossy_sphere,
    _sphere_submesh,
    _theta_cf,
    _to_sub,
)

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"  # committed, EX-57
BASENAME = "materials_21_point_probes_on_gel_ball_submesh"
BOX_HALF_WIDTH = 0.02  # test_lossy_sphere_sar.py:69 (the air box the figure shows)
A_BALL = 2.0e-3  # the mat:13 probe radius, a stated fixture parameter

#: Probe points (example parameters, not records; identical on every rank).
B_LABELS = ("B0 centre", "B1 3.5 mm on x", "B2 6.5 mm off-axis", "B3 near wall 9.5 mm",
            "B4 mesh vertex", "B5 outside (air)")
B_PTS = np.array(
    [
        [0.0, 0.0, 0.0],
        [3.5e-3, 0.0, 0.0],
        list(6.5e-3 * np.ones(3) / np.sqrt(3.0)),  # mat:13's c2
        [0.0, 0.0, 9.5e-3],
        [0.0, 0.0, 5.0e-3],  # replaced by the nearest mesh vertex below
        [0.0, 0.0, 1.5e-2],
    ],
    dtype=np.float64,
)
B_IN = (0, 1, 2, 3, 4)
B_OUT = 5

#: Records: (value, log:line).  Printed beside the reading, never asserted.
POST10_OWNERS = [217, 223, 126, 229, 132, -1]
POST10_OWNER_LOGS = "20261005T200223Z_POST-10-step1-W1-real.log:98 / ...W2.log:146 / ...W3.log:77"
POST10_N1 = {2: (1.554e-15, "20261005T200223Z_POST-10-step1-W1-real.log:99"),
             4: (2.220e-15, "20261005T200334Z_POST-10-step1-W2.log:147"),
             1: (1.776e-15, "20261005T200500Z_POST-10-step1-W3.log:78")}
MAT10_LOG = "20260930T021802Z_MAT-10-step1-W1.log.gz (zcat line numbers)"
REC = {
    "q_cf": (2002.088514, ":157"), "th0": (0.122350, ":157"), "thR": (0.066736, ":157"),
    "cells": ((74020, 19496), ":158"), "B": (2.807e-15, ":160"), "P": (3.1923e-03, ":161"),
    "A": (3.2464e-02, ":162"), "C": (1.233e-13, ":164"), "a": (7.9323e-01, ":166"),
    "b": (1.0649e+00, ":167"), "t_em": (15.5, ":168"), "t_tot": (23.6, ":168"),
}

_T0 = time.perf_counter()


def _say(comm, msg):
    if comm.rank == 0:
        print(msg, flush=True)


def _clock(comm, what):
    _say(comm, f"[EX-121 clock] {what}: {time.perf_counter() - _T0:.1f} s since start")


def _ev(f, pts, comm, **kw):
    v, ok = evaluate_vector_field_parallel(f, pts, comm, **kw)
    return np.asarray(v).reshape(len(pts), -1), np.asarray(ok, dtype=bool)


def _raises(comm, exc_type, fn):
    """Call ``fn`` on every rank, then allgather the raise flags (after the call)."""
    raised, msg = False, ""
    try:
        fn()
    except exc_type as exc:
        raised, msg = True, str(exc)
    return comm.allgather(raised), msg


def _claims(comm, local_cells):
    return comm.allreduce((np.asarray(local_cells) >= 0).astype(np.int64), op=MPI.SUM)


def _closed_form_theta(r, q_cf):
    return robin_sphere_steady_temperature(
        r, radius=R, k=K, q=q_cf, h=H, t_inf=T_INF
    ) - T_INF


# --------------------------------------------------------------------------------
# Part A: POST-10's own cube, on the gate's fixture and band
# --------------------------------------------------------------------------------
def part_a(comm):
    msh, keep, sub, parent_cells = _fixture(comm)
    sets = _oracle(msh, keep, comm)
    parent = (msh, parent_cells)
    orig = _original(msh)

    # (S0)
    simap = sub.topology.index_map(3)
    n_sub_owned = simap.size_local
    gathered = np.concatenate(comm.allgather(_original(msh)[parent_cells[:n_sub_owned]]))
    s0 = (simap.size_global == SUB_CELLS_GLOBAL
          and len(parent_cells) == simap.size_local + simap.num_ghosts
          and gathered.size == SUB_CELLS_GLOBAL
          and np.unique(gathered).size == SUB_CELLS_GLOBAL)
    _say(comm, f"[A S0] sub size_global {simap.size_global} (gate {SUB_CELLS_GLOBAL}); kept originals "
               f"{gathered.size}, unique {np.unique(gathered).size}; ASSERTED")
    assert s0

    sub_oo, sub_local = point_owner_cells(sub, POINTS, comm, parent=parent)
    par_oo, _ = point_owner_cells(msh, POINTS, comm, candidate_cells=keep)
    f_par = _dg0(msh, orig)
    f_sub = _dg0(sub, orig[parent_cells])
    sub_vals, sub_valid = _ev(f_sub, POINTS, comm, owner=RULE, parent=parent)
    par_vals, par_valid = _ev(f_par, POINTS, comm, owner=RULE, candidate_cells=keep)
    sub_vals, par_vals = sub_vals[:, 0], par_vals[:, 0]

    v_n1 = fem.functionspace(sub, ("N1curl", 1))
    f_n1 = fem.Function(v_n1)
    f_n1.interpolate(lambda x: A_VEC[:, None] + np.cross(B_VEC, x.T).T)
    n1_vals, n1_valid = _ev(f_n1, POINTS, comm, owner=RULE, parent=parent)
    n1_exact = A_VEC[None, :] + np.cross(B_VEC[None, :], POINTS)
    err_n1 = float(np.max(np.abs(n1_vals[list(INSIDE)] - n1_exact[list(INSIDE)])))
    claims = _claims(comm, sub_local)

    owners = [int(x) for x in sub_oo]
    _say(comm, f"[A owner] P0-P5: {owners}")
    _say(comm, f"[A] owner integers beside POST-10 record {POST10_OWNERS} ({POST10_OWNER_LOGS}): "
               f"{'identical' if owners == POST10_OWNERS else 'DIFFER'}")
    rec_n1 = POST10_N1[comm.size]
    _say(comm, f"[A] N1curl max |err| P0-P4 {err_n1:.3e} (record {rec_n1[0]:.3e}, {rec_n1[1]}); "
               f"ASSERTED <= {CLOSED_FORM_TOL:.0e}")
    for i, lab in enumerate(LABELS):
        s = sets[i]
        _say(comm, f"[A] {lab:15s} kept set size {len(s)}, min {s[0] if s else None}, "
                   f"max {s[-1] if s else None}, sep {(s[-1] - s[0]) if s else None}; "
                   f"owner {owners[i]}; DG0 {sub_vals[i].real:.1f}; claims {int(claims[i])}")

    # (C1)
    for i in (1, 2, 4):
        assert len(sets[i]) >= 2 and sets[i][-1] - sets[i][0] >= 1, (i, sets[i])
    _say(comm, "[A C1] oracle keep-set spreads at P1 / P2 / P4: "
               + " / ".join(str(sets[i][-1] - sets[i][0]) for i in (1, 2, 4))
               + " (POST-10 outcome 9 / 131 / 121); ASSERTED >= 1 against band 0")
    for i in INSIDE:
        assert len(sets[i]) >= 1 and len(set(sets[i])) == len(sets[i])
    assert sets[OUTSIDE] == []
    # (K)
    for i in INSIDE:
        assert int(sub_oo[i]) == sets[i][0], (i, sub_oo, sets[i])
        assert bool(sub_valid[i])
        assert np.imag(sub_vals[i]) == 0.0
        assert sub_vals[i].real == float(sets[i][0]), (i, sub_vals[i], sets[i])
        assert int(claims[i]) == 1, (i, claims)
    assert int(sub_oo[OUTSIDE]) == -1 and not bool(sub_valid[OUTSIDE]) and int(claims[OUTSIDE]) == 0
    # (I)
    assert np.array_equal(sub_oo, par_oo), (sub_oo, par_oo)
    assert np.array_equal(sub_vals, par_vals), (sub_vals, par_vals)
    assert np.array_equal(sub_valid, par_valid), (sub_valid, par_valid)
    # (N)
    assert all(bool(n1_valid[i]) for i in INSIDE) and not bool(n1_valid[OUTSIDE])
    assert err_n1 <= CLOSED_FORM_TOL, err_n1
    # (W)
    for arr in (sub_oo, sub_vals, n1_vals):
        assert _bitwise_equal_across_ranks(comm, np.ascontiguousarray(arr))
    _say(comm, "[A] (K) (I) (N) (W) hold; ASSERTED at the POST-10 bands (K, I, W band 0)")

    # (C2) cube
    flags, message = _raises(comm, RuntimeError,
                             lambda: evaluate_vector_field_parallel(f_sub, POINTS, comm, owner=RULE))
    _say(comm, f"[A C2] no parent: flags {flags}; message {message!r}")
    assert flags == [True] * comm.size and "parent=" in message and "create_submesh" in message
    # (C3)
    flags3, message3 = _raises(
        comm, ValueError,
        lambda: evaluate_vector_field_parallel(f_sub, POINTS, comm, parent=parent))
    _say(comm, f"[A C3] parent without the rule: flags {flags3}; message {message3!r}")
    assert flags3 == [True] * comm.size and "parent is valid only" in message3
    # (C4)
    flags4, message4 = _raises(
        comm, ValueError,
        lambda: evaluate_vector_field_parallel(
            f_sub, POINTS, comm, owner=RULE, parent=(msh, parent_cells[:-1])))
    _say(comm, f"[A C4] short parent_cells: flags {flags4}; message {message4!r}")
    assert flags4 == [True] * comm.size and "parent" in message4
    return owners


# --------------------------------------------------------------------------------
# Part B: MAT-10's gel ball
# --------------------------------------------------------------------------------
def part_b(comm):
    t0 = time.perf_counter()
    run = _solve_lossy_sphere(SIGMA_HIGH, R / 10.0, R / 5.0, return_fields=True)
    msh, cell_tags = run["mesh"], run["cell_tags"]
    E, sigma = run["e_complex"], run["sigma_field"]
    t_em = time.perf_counter() - t0

    write_setup_figure(
        msh,
        cell_tags,
        FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names={1: "gel sphere (create_submesh, 19 496 cells)", 2: "air (parent only)"},
        translucent_tags=(2,),
        slice_normal=(0.0, 1.0, 0.0),
        slice_origin=(0.0, 0.0, 0.0),
        window_size=(1400, 630),
        title="mat:21 - probes B0 centre, B1 3.5 mm, B2 6.5 mm, B3 9.5 mm, B4 vertex, B5 outside",
        comm=comm,
    )

    cf = uniform_sphere_sar_closed_form(
        e0=E0, epsilon_r=EPSILON_R_SPHERE, sigma=SIGMA_HIGH, omega=OMEGA,
        rho=RHO_KG_M3, epsilon_0=EPSILON_0,
    )
    q_cf = G * RHO_KG_M3 * cf["sar_w_per_kg"]
    th0 = float(_closed_form_theta(0.0, q_cf))
    thR = float(_closed_form_theta(R, q_cf))
    _say(comm, f"[B] q_cf = {q_cf:.6f} W/m^3  theta(0) = {th0:.6f} K  theta(R) = {thR:.6f} K "
               f"(MAT-10 record {REC['q_cf'][0]} / {REC['th0'][0]} / {REC['thR'][0]}, "
               f"{MAT10_LOG}{REC['q_cf'][1]})")

    # the gate's submesh and the parent keep set (the _sphere_submesh idiom, four lines)
    sub, parent_ids = _sphere_submesh(msh, cell_tags)
    ft = _robin_tags(sub)
    imap = msh.topology.index_map(3)
    n_parent_loc = imap.size_local + imap.num_ghosts
    parent_values = np.full(n_parent_loc, -1, dtype=np.int32)
    parent_values[cell_tags.indices] = cell_tags.values
    keep = np.flatnonzero(parent_values == SPHERE_TAG).astype(np.int32)
    parent = (msh, parent_ids)
    n_parent = imap.size_global
    n_sub = sub.topology.index_map(3).size_global
    n_own = sub.topology.index_map(3).size_local
    _say(comm, f"[B] cells: parent = {n_parent}  sub = {n_sub} (record {REC['cells'][0]}, "
               f"{MAT10_LOG}{REC['cells'][1]})")

    # q_h as the gate builds it, plus the parent DG0 as a Function
    q_expr = G * 0.5 * sigma * ufl.inner(E, E)
    Wp, q_par_arr = _cell_average_dg0(msh, cell_tags, q_expr)
    q_h = _to_sub(Wp, q_par_arr, sub, parent_ids)
    q_par = fem.Function(Wp)
    q_par.x.array[:] = q_par_arr
    q_par.x.scatter_forward()
    qmax = comm.allreduce(float(np.max(np.real(q_h.x.array[:n_own]))), op=MPI.MAX)
    qmin = comm.allreduce(float(np.min(np.real(q_h.x.array[:n_own]))), op=MPI.MIN)
    _say(comm, f"[B] q_h max = {qmax:.6f}  min = {qmin:.6f} W/m^3 (brief: 2582.4 / 1701.1)")

    # (B)
    int_q_sub = float(np.real(_allsum(comm, q_h * ufl.dx(domain=sub))))
    p_parent = G * run["dissipated_power_w"]
    err_B = abs(int_q_sub - p_parent) / abs(p_parent)
    _say(comm, f"[B] (B) rel = {err_B:.3e} (record {REC['B'][0]:.3e}, {REC['B'][1]}); "
               f"ASSERTED <= {B_BAND:.0e}")

    chi = _interior_indicator(sub)
    theta = _theta_cf(sub, q_cf)
    res_P = solve_steady_conduction(sub, ft, k=K, q=q_cf, h=H, t_inf=T_INF,
                                    robin_tag=ROBIN_TAG, degree=1)
    err_P = _interior_rel(sub, chi, (res_P.temperature - T_INF) - theta, theta)
    _say(comm, f"[B] (P) uniform q_cf interior rel L2 = {err_P:.4e} (record {REC['P'][0]:.4e}, "
               f"{REC['P'][1]}); ASSERTED <= {P_BAND:.0e}")
    res_A = solve_steady_conduction(sub, ft, k=K, q=q_h, h=H, t_inf=T_INF,
                                    robin_tag=ROBIN_TAG, degree=1)
    err_A = _interior_rel(sub, chi, (res_A.temperature - T_INF) - theta, theta)
    _say(comm, f"[B] (A) theta_h[q_h] interior rel L2 = {err_A:.4e} (record {REC['A'][0]:.4e}, "
               f"{REC['A'][1]}); ASSERTED <= {A_BAND}")
    err_C = abs(res_A.convected_power - res_A.source_power) / abs(res_A.source_power)
    _say(comm, f"[B] (C) balance rel = {err_C:.3e} (record {REC['C'][0]:.3e}, {REC['C'][1]}); "
               f"ASSERTED <= {C_BAND:.0e}")

    q_re_expr = G * 0.5 * sigma * ufl.inner(ufl.real(E), ufl.real(E))
    _, q_re_par = _cell_average_dg0(msh, cell_tags, q_re_expr)
    q_re = _to_sub(Wp, q_re_par, sub, parent_ids)
    res_a = solve_steady_conduction(sub, ft, k=K, q=q_re, h=H, t_inf=T_INF,
                                    robin_tag=ROBIN_TAG, degree=1)
    err_a = _interior_rel(sub, chi, (res_a.temperature - T_INF) - theta, theta)
    q_two = fem.Function(q_h.function_space)
    q_two.x.array[:] = 2.0 * q_h.x.array
    q_two.x.scatter_forward()
    res_b = solve_steady_conduction(sub, ft, k=K, q=q_two, h=H, t_inf=T_INF,
                                    robin_tag=ROBIN_TAG, degree=1)
    err_b = _interior_rel(sub, chi, (res_b.temperature - T_INF) - theta, theta)
    _say(comm, f"[B] control (a) Re-E-only miss {err_a:.4e} (record {REC['a'][0]:.4e}, {REC['a'][1]}); "
               f"control (b) missing 1/2 miss {err_b:.4e} (record {REC['b'][0]:.4e}, {REC['b'][1]}); "
               f"ASSERTED > {A_BAND}")
    _say(comm, f"[B] elapsed: EM {t_em:.1f} s (record {REC['t_em'][0]}), fixture+solves "
               f"{time.perf_counter() - t0:.1f} s (record total {REC['t_tot'][0]}, {REC['t_tot'][1]})")
    assert err_B <= B_BAND, err_B
    assert err_P <= P_BAND, err_P
    assert err_C <= C_BAND, err_C
    assert err_A <= A_BAND, err_A
    assert err_a > A_BAND and err_b > A_BAND, (err_a, err_b)
    _clock(comm, "MAT-10 anchors done")

    # ---- B4: the submesh vertex nearest (0, 0, 5e-3) -----------------------------
    target = B_PTS[4].copy()
    xs = sub.geometry.x
    d = np.linalg.norm(xs - target[None, :], axis=1)
    idx = int(np.lexsort((xs[:, 2], xs[:, 1], xs[:, 0], d))[0])
    best = (float(d[idx]), float(xs[idx, 0]), float(xs[idx, 1]), float(xs[idx, 2]))
    allbest = comm.allgather(best)
    b4 = min(allbest)
    pts = B_PTS.copy()
    pts[4] = np.array(b4[1:], dtype=np.float64)
    _say(comm, f"[B4] xyz: ({pts[4, 0]:.17e}, {pts[4, 1]:.17e}, {pts[4, 2]:.17e}) distance to "
               f"(0, 0, 5e-3): {b4[0]:.3e} m")

    # ---- owner integers and the identities ------------------------------------------
    sub_oo, sub_local = point_owner_cells(sub, pts, comm, parent=parent)
    par_oo, par_local = point_owner_cells(msh, pts, comm, candidate_cells=keep)
    claims_sub = _claims(comm, sub_local)
    claims_par = _claims(comm, par_local)
    owners = [int(x) for x in sub_oo]
    _say(comm, f"[B owner] B0-B5: {owners}")

    # (I-q)
    qs, qs_ok = _ev(q_h, pts, comm, owner=RULE, parent=parent)
    qp, qp_ok = _ev(q_par, pts, comm, owner=RULE, candidate_cells=keep)
    qs, qp = qs[:, 0], qp[:, 0]
    _say(comm, f"[B I-q] owner ints equal {np.array_equal(sub_oo, par_oo)}; q_h bitwise equal "
               f"{np.array_equal(qs, qp)}; validity equal {np.array_equal(qs_ok, qp_ok)}; "
               f"q_h at B0-B4: " + " ".join(f"{qs[i].real:.6f}" for i in B_IN) + "; ASSERTED band 0")
    assert np.array_equal(sub_oo, par_oo), (sub_oo, par_oo)
    assert np.array_equal(qs, qp), (qs, qp)
    assert np.array_equal(qs_ok, qp_ok)

    # (I-o)
    orig = _original(msh)
    o_par = _dg0(msh, orig)
    o_sub = _dg0(sub, orig[parent_ids])
    ov_s, ov_s_ok = _ev(o_sub, pts, comm, owner=RULE, parent=parent)
    ov_p, ov_p_ok = _ev(o_par, pts, comm, owner=RULE, candidate_cells=keep)
    ov_s, ov_p = ov_s[:, 0], ov_p[:, 0]
    assert np.array_equal(ov_s, ov_p), (ov_s, ov_p)
    assert np.array_equal(ov_s_ok, ov_p_ok)
    for i in B_IN:
        assert bool(ov_s_ok[i]) and np.imag(ov_s[i]) == 0.0
        assert ov_s[i].real == float(sub_oo[i]), (i, ov_s[i], sub_oo[i])
        assert int(claims_sub[i]) == 1 and int(claims_par[i]) == 1, (i, claims_sub, claims_par)
    _say(comm, "[B I-o] planted original-index DG0 bitwise equal sub vs parent; value == owner_original "
               "at B0-B4; one claim each; ASSERTED band 0")

    # (O)
    assert int(sub_oo[B_OUT]) == -1 and not bool(qs_ok[B_OUT]) and int(claims_sub[B_OUT]) == 0
    assert int(par_oo[B_OUT]) == -1 and not bool(qp_ok[B_OUT]) and int(claims_par[B_OUT]) == 0
    _say(comm, f"[B O] B5: owner {int(sub_oo[B_OUT])} / {int(par_oo[B_OUT])}, valid "
               f"{bool(qs_ok[B_OUT])} / {bool(qp_ok[B_OUT])}, claims {int(claims_sub[B_OUT])} / "
               f"{int(claims_par[B_OUT])} (submesh / parent keep); ASSERTED")

    # (N-ball)
    v_n1 = fem.functionspace(sub, ("N1curl", 1))
    f_n1 = fem.Function(v_n1)
    f_n1.interpolate(lambda x: A_VEC[:, None] + np.cross(B_VEC, x.T).T)
    n1_vals, n1_ok = _ev(f_n1, pts, comm, owner=RULE, parent=parent)
    n1_exact = A_VEC[None, :] + np.cross(B_VEC[None, :], pts)
    err_n1 = float(np.max(np.abs(n1_vals[list(B_IN)] - n1_exact[list(B_IN)])))
    _say(comm, f"[B N-ball] N1curl a + b x x max |err| B0-B4 {err_n1:.3e}; ASSERTED <= {CLOSED_FORM_TOL:.0e}")
    assert all(bool(n1_ok[i]) for i in B_IN) and not bool(n1_ok[B_OUT])
    assert err_n1 <= CLOSED_FORM_TOL, err_n1

    # (C2-ball)
    flags, message = _raises(comm, RuntimeError,
                             lambda: evaluate_vector_field_parallel(q_h, pts, comm, owner=RULE))
    _say(comm, f"[B C2-ball] no parent on the ball submesh: flags {flags}; message {message!r}")
    assert flags == [True] * comm.size and "parent=" in message and "create_submesh" in message

    # ---- theta readings (printed only) -------------------------------------------------
    thP, thP_ok = _ev(res_P.temperature, pts, comm, owner=RULE, parent=parent)
    thA, thA_ok = _ev(res_A.temperature, pts, comm, owner=RULE, parent=parent)
    thP = np.real(thP[:, 0]) - T_INF
    thA = np.real(thA[:, 0]) - T_INF
    r_pts = np.linalg.norm(pts, axis=1)
    th_ex = np.array([float(_closed_form_theta(min(r, R), q_cf)) for r in r_pts])
    _say(comm, "[B theta] printed only, no pointwise theta band is gated. "
               "Predicted |theta_h[q_cf]/theta - 1| <~ 1e-2; theta_h[q_h]/theta - 1 ~ +3 to +4 %")
    for i in B_IN:
        _say(comm, f"[B theta] {B_LABELS[i]:22s} r={r_pts[i]:.4e} m theta_cf={th_ex[i]:.6f} K  "
                   f"theta_h[q_cf]={thP[i]:.6f} ({thP[i] / th_ex[i] - 1:+.4e})  "
                   f"theta_h[q_h]={thA[i]:.6f} ({thA[i] / th_ex[i] - 1:+.4e})")

    # owner=None against owner=RULE (continuous P1: predicted <= 1e-15 relative)
    try:
        legacy, _ = _ev(res_A.temperature, pts, comm)
        legacy_err = None
    except Exception as exc:  # printed only
        legacy, legacy_err = None, repr(exc)
    if legacy is not None:
        dmax = max(abs(np.real(legacy[i, 0]) - T_INF - thA[i]) for i in B_IN)
        _say(comm, f"[B owner=None] theta_h[q_h] owner=None vs owner=RULE parent=: max |diff| {dmax:.3e} K "
                   f"(relative {dmax / np.max(np.abs(thA[list(B_IN)])):.3e}; predicted <= 1e-15)")
    else:
        _say(comm, f"[B owner=None] raised: {legacy_err}")

    # 2 mm ball probe beside the point value (cB)
    ball = {}
    for i in (0, 1, 2):
        val, vol = ball_average(res_P.temperature, ProbeBall(tuple(pts[i]), A_BALL))
        ball[i] = float(val) - T_INF
        _say(comm, f"[B ball] {B_LABELS[i]:22s} ball_2mm={ball[i]:.6f} K point={thP[i]:.6f} K "
                   f"ball/point - 1 = {ball[i] / thP[i] - 1:+.4e} (V {vol:.4e} m^3)")
    cB = -q_cf * A_BALL**2 / (10.0 * K) / th0
    _say(comm, f"[B cB] predicted centre offset -q a^2/(10k)/theta(0) = {cB:+.4e} (-1.09 %); "
               f"|offset| / P_BAND = {abs(cB) / P_BAND:.2f}x; printed, not asserted")

    # B4 incident-cell spread (cV)
    tree = geometry.bb_tree(sub, 3)
    b4p = pts[4:5]
    hits = geometry.compute_collisions_points(tree, b4p)
    coll = geometry.compute_colliding_cells(sub, hits, b4p)
    links = np.asarray(coll.links(0))
    owned = links[links < n_own]
    qdm = q_h.function_space.dofmap
    local = [(float(np.real(q_h.x.array[qdm.cell_dofs(int(c))[0]])), int(orig[parent_ids[int(c)]]))
             for c in owned]
    inc = [x for part in comm.allgather(local) for x in part]
    inc.sort(key=lambda t: t[1])
    vals = [v for v, _ in inc]
    spread = (max(vals) - min(vals)) if vals else float("nan")
    q_rule = float(qs[4].real)
    lowest = inc[0][0] if inc else float("nan")
    _say(comm, f"[B cV] B4 incident sub cells {len(inc)} (predicted >= 4); q_h spread {spread:.6f} W/m^3 "
               f"beside the rule's value {q_rule:.6f}; factor spread/rule {spread / q_rule:.4e} "
               f"(ceiling: whole-ball range {(qmax - qmin) / q_cf:.3f}); "
               f"rule value == q_h of the lowest-original incident cell: {q_rule == lowest}; printed, not asserted")
    _say(comm, "[B cV] incident (q_h, owner_original): " + str(inc))

    # (W-ball)
    for name, arr in (("owner", sub_oo), ("q_h", qs), ("theta_q_cf", thP), ("theta_q_h", thA),
                      ("n1", n1_vals)):
        ok = _bitwise_equal_across_ranks(comm, np.ascontiguousarray(arr))
        assert ok, name
    _say(comm, "[B W-ball] owner ints, q_h, theta (both), N1 values bitwise equal on every rank; ASSERTED")

    # ---- outputs -------------------------------------------------------------------------
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    V1 = res_P.temperature.function_space

    def _p1(name, arr):
        f = fem.Function(V1, name=name)
        f.x.array[:] = arr
        f.x.scatter_forward()
        return f

    th_cf_f = _p1("theta_q_cf", np.real(res_P.temperature.x.array) - T_INF)
    th_h_f = _p1("theta_q_h", np.real(res_A.temperature.x.array) - T_INF)
    th_ex_f = fem.Function(V1, name="theta_exact")
    th_ex_f.interpolate(lambda x: np.asarray(
        _closed_form_theta(np.minimum(np.sqrt(np.sum(x**2, axis=0)), R), q_cf), dtype=np.float64))
    th_ex_f.x.scatter_forward()
    th_err_f = _p1("theta_error", np.real(th_cf_f.x.array) - np.real(th_ex_f.x.array))
    q_out = fem.Function(q_h.function_space, name="q_h")
    q_out.x.array[:] = np.real(q_h.x.array)
    q_out.x.scatter_forward()
    o_out = fem.Function(o_sub.function_space, name="owner_original_dg0")
    o_out.x.array[:] = np.real(o_sub.x.array)
    o_out.x.scatter_forward()
    xdmf_path, _ = write_xdmf_time_series(
        OUTPUT_DIR / f"{BASENAME}_combined", sub, None,
        [(0.0, {"theta_q_cf": th_cf_f, "theta_q_h": th_h_f, "theta_exact": th_ex_f,
                "theta_error": th_err_f, "q_h": q_out, "owner_original_dg0": o_out})],
        comm=comm,
    )
    csv_path = OUTPUT_DIR / f"{BASENAME}_probes.csv"
    if comm.rank == 0:
        with csv_path.open("w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["label", "x_m", "y_m", "z_m", "owner_original", "valid", "q_h",
                        "theta_q_cf", "theta_q_h", "theta_exact", "ball_2mm", "incident_spread_q_h"])
            for i in range(6):
                inside = i in B_IN
                w.writerow([
                    B_LABELS[i], *(f"{v:.9e}" for v in pts[i]), int(sub_oo[i]), int(bool(qs_ok[i])),
                    f"{qs[i].real:.9e}" if inside else "",
                    f"{thP[i]:.9e}" if inside else "", f"{thA[i]:.9e}" if inside else "",
                    f"{th_ex[i]:.9e}" if inside else "",
                    f"{ball[i]:.9e}" if i in ball else "",
                    f"{spread:.9e}" if i == 4 else "",
                ])
    adopt_host_ownership(OUTPUT_DIR, comm=comm)
    _say(comm, f"[output] {xdmf_path}")
    _say(comm, f"[output] {csv_path}")
    return owners, pts


def main() -> None:
    comm = MPI.COMM_WORLD
    if not np.issubdtype(np.dtype(default_scalar_type), np.complexfloating):
        raise RuntimeError(
            "this example needs the complex DolfinX build "
            "(source /usr/local/bin/dolfinx-complex-mode); the runner's `mat:` "
            "group sources it automatically"
        )
    assert comm.size in (2, 4), f"run at -n 2 (record) or -n 4 (control), got {comm.size}"
    _say(comm, "=" * 78)
    _say(comm, f"EX-121 -- point probes on the gel-ball submesh (-n {comm.size})")
    _say(comm, "=" * 78)
    _say(comm, f"[bands] imported unchanged: CLOSED_FORM_TOL {CLOSED_FORM_TOL:.0e}, B_BAND {B_BAND:.0e}, "
               f"P_BAND {P_BAND:.0e}, A_BAND {A_BAND}, C_BAND {C_BAND:.0e}; box half width {BOX_HALF_WIDTH}")
    part_a(comm)
    _clock(comm, "Part A done")
    part_b(comm)
    _say(comm, f"[timing] total {time.perf_counter() - _T0:.1f} s at -n {comm.size}")
    _say(comm, "\nAll assertions hold.")


if __name__ == "__main__":
    main()
