"""Example (`EX-123`): the `Study` transient API on thermal-only gel balls.

`mri:24` and earlier drive `Study` with an EM solve.  This one uses no EM
solve at all: a `Study` with ``material=None`` and a `UniformHeatSource`, and
``Study.temperature(RobinCooling, transient=TransientSpec, duty=DutyCycle)``
over a `ThermalBallGeometry` with or without a `ContainerWall`.  It reproduces
`MAT-17` step 2 (cases A and B, the duty-cycled ball) and `MAT-18` step 2(b)
(the walled ball at t*) **through the five public dataclasses**, and composes
the two once (printed only, identities asserted).

**It asserts, and it does not re-implement** (the `ANS-1` rule).  Fixture,
helpers, comparands, bands and records are imported from
``tests/package/test_wf15_study_transient.py`` (`WF-15` step 3, item 330) as
``g``.  The one substitution (the `mat:11` precedent): the duty legs use the
complex-safe ``m13._interior_error`` where the gate uses `MAT-9`'s, which calls
``float()`` on a complex value.  The gate ran in the real build; this example
runs in the complex build, the Study path's first complex-build run.

Run it::

    ./scripts/run_examples.sh -e mri:25 -n 2 -t 590

Output (``examples/mri/paraview_output/``):
``mri_25_study_transient_duty_and_wall_duty_combined.xdmf`` (bare ball) and
``mri_25_study_transient_duty_and_wall_wall_combined.xdmf`` (walled ball).
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
from dolfinx import default_scalar_type, fem
from mpi4py import MPI

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from fem_em_solver import (  # noqa: E402
    ContainerWall,
    DutyCycle,
    RobinCooling,
    Study,
    ThermalBallGeometry,
    TransientSpec,
    UniformHeatSource,
)
from fem_em_solver.io.paraview_utils import (  # noqa: E402
    adopt_host_ownership,
    write_xdmf_with_tags,
)
from fem_em_solver.post.evaluation import evaluate_vector_field_parallel  # noqa: E402
from fem_em_solver.post.setup_figure import write_setup_figure  # noqa: E402

import tests.package.test_wf15_study_transient as g  # noqa: E402
from tests.validation import test_mat9_transient_conduction_sphere as mat9  # noqa: E402
from tests.validation import test_mat13_transient_perfusion_sphere as m13  # noqa: E402

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"  # committed, EX-57
BASENAME = "mri_25_study_transient_duty_and_wall"

#: Gate records (log:line), printed beside each reading.  The .4e strings are
#: asserted at 1e-3 relative through the gate module's own constants; these
#: 10-digit values are printed only.
REC = {
    "A_cells": 151548,            # 20261016 W1a log :88
    "A_s0": (49, 46),             # W1a :88-89
    "A_it": 8.146e-15,            # W1a :90
    "A_g": 4.3869840018e-04,      # W1a :91
    "A_bal": 2.759e-14,           # W1a :91
    "B_s0": (9, 12),              # W1a :93
    "B_it": 2.080e-14,            # W1a :95
    "B_g": 1.0682780964e-03,      # W1a :96
    "B_bal": 1.061e-14,           # W1a :96
    "c1": 9.7087e-01,             # W1a :98
    "wall_cells": 215542,         # W1b :88
    "c3_sp": 2.7755575615628914e-17,  # W1b :89
    "wall_it": 1.433e-14,         # W1b :90
    "wall_g": 8.2066325346e-04,   # W1b :91
    "wall_whole": 1.0874e-03,     # W1b :91
    "wall_shell": 2.2248e-03,     # W1b :91
    "wall_bal": 1.146e-14,        # W1b :91
    "c2": 5.3140e-02,             # W1b :93
}


def _say(comm, msg):
    if comm.rank == 0:
        print(msg, flush=True)


def _raises(exc, fn):
    try:
        fn()
    except exc:
        return True
    return False


def _max_im(comm, fs):
    m = 0.0
    for f in fs:
        n = f.function_space.dofmap.index_map.size_local * f.function_space.dofmap.index_map_bs
        if n:
            m = max(m, float(np.max(np.abs(np.imag(f.x.array[:n])))))
    return comm.allreduce(m, op=MPI.MAX)


def _theta_at(T, pts):
    vals, ok = evaluate_vector_field_parallel(T, pts)
    out = []
    for i in range(pts.shape[0]):
        if not bool(ok[i]):
            raise RuntimeError(f"point {pts[i]} not found")
        out.append(float(np.real(np.ravel(vals[i])[0])) - g.T_INF)
    return np.array(out)


def _named(V, arr_fn, name):
    f = fem.Function(V, name=name)
    arr_fn(f)
    f.x.scatter_forward()
    return f


def _p1_minus_tinf(T, name):
    f = fem.Function(T.function_space, name=name)
    f.x.array[:] = np.real(T.x.array) - g.T_INF
    f.x.scatter_forward()
    return f


def _p1_of(u, V1, name):
    f = fem.Function(V1, name=name)
    f.interpolate(u)
    f.x.scatter_forward()
    return f


def _diff(a, b, V1, name):
    f = fem.Function(V1, name=name)
    f.x.array[:] = a.x.array - b.x.array
    f.x.scatter_forward()
    return f


def main() -> None:
    comm = MPI.COMM_WORLD
    started = time.perf_counter()

    if not np.issubdtype(np.dtype(default_scalar_type), np.complexfloating):
        raise RuntimeError(
            "this example needs the complex DolfinX build "
            "(source /usr/local/bin/dolfinx-complex-mode); the runner's `mri:` "
            "group sources it automatically"
        )
    assert comm.size == 2, f"EX-123 is a -n 2 example (comm.size {comm.size})"

    # Structural: the walled core is MAT-9's ball; the complex-safe twin shares
    # its constants; the case table is the gate's.
    assert (g.K1, g.RC1, g.Q1, g.A_) == (g.K, g.RHO_C, g.Q, g.R)
    assert g.CASES["A"] == (g.P_A, g.D_A)
    assert (m13.R, m13.K, m13.Q, m13.H, m13.T_INF, m13.RHO_C, m13.ROBIN_TAG, m13.QDEG) == (
        g.R, g.K, g.Q, g.H, g.T_INF, g.RHO_C, mat9.ROBIN_TAG, mat9.QDEG)

    _say(comm, "=" * 78)
    _say(comm, "EX-123 -- the Study transient API: duty cycle + container wall, no EM solve")
    _say(comm, "=" * 78)
    _say(comm, f"[gate] bands imported: L0 {g.L0_BAND:.0e}, interior {g.INTERIOR_BAND:.0e}, "
               f"balance {g.BALANCE_BAND:.0e}; build complex, -n {comm.size}")

    all_T = []  # every theta the example produces, for (Im)
    asserts = []  # (label, ok) -- everything asserted at the end of its leg

    # ------------------------------------------------------------ (K) no solve
    # Copied from test_wf15_study_transient.py:186-219 (provenance), pytest.raises
    # replaced by a local helper.
    spec = TransientSpec(150.0, 94, g.RHO_C)
    assert (spec.theta, spec.startup_half_steps, spec.t_init_c) == (0.5, 2, None)
    geo = ThermalBallGeometry(g.A_, g.RES, wall=ContainerWall(g.B_ - g.A_, g.K2, g.RC2))
    assert geo.outer_radius == g.B_
    assert ThermalBallGeometry(g.R, g.RES).outer_radius == g.R
    n_chi = 0
    for P, D in (g.CASES["A"], g.CASES["B"]):
        d, chi = DutyCycle(P, D), g._chi(P, D)
        ts = np.concatenate([np.linspace(0.0, 13 * P, 4001), np.arange(14) * P,
                             np.arange(14) * P + D * P])
        assert all(d(float(t)) == chi(float(t)) for t in ts), "(chi)"
        n_chi += ts.size
    bads = (lambda: TransientSpec(0.0, 10, g.RHO_C), lambda: TransientSpec(1.0, 0, g.RHO_C),
            lambda: TransientSpec(1.0, 10, g.RHO_C, theta=1.5),
            lambda: TransientSpec(1.0, 10, g.RHO_C, startup_half_steps=1),
            lambda: DutyCycle(100.0, 0.0), lambda: DutyCycle(100.0, 1.5),
            lambda: DutyCycle(-1.0, 0.5), lambda: ContainerWall(0.0, g.K2, g.RC2),
            lambda: UniformHeatSource(0.0), lambda: ThermalBallGeometry(0.0, g.RES))
    for i, bad in enumerate(bads):
        assert _raises(ValueError, bad), f"(K) ValueError construction {i}"
    assert _raises(TypeError, lambda: Study(ThermalBallGeometry(g.R, g.RES), "gel",
                                            UniformHeatSource(g.Q)))
    assert _raises(TypeError, lambda: Study(ThermalBallGeometry(g.R, g.RES), None, None))
    st = Study(ThermalBallGeometry(g.R, g.RES), None, UniformHeatSource(g.Q))
    rc = RobinCooling(g.K, g.H, g.T_INF)
    assert _raises(TypeError, st.solve)
    assert _raises(TypeError, lambda: st.temperature(rc))
    assert _raises(TypeError, lambda: st.temperature(rc, transient=spec, duty=(g.P_A, g.D_A)))
    assert _raises(ValueError, lambda: st.temperature(
        RobinCooling(g.K, g.H, g.T_INF, source_scale=2.0), transient=spec))
    assert st.result is None
    _say(comm, f"[K] spec defaults (0.5, 2, None), outer_radius, {len(bads)} ValueError "
               f"constructions, 5 refusals on a thermal-only Study, result None: green (no solve)")
    _say(comm, f"[chi] DutyCycle(P, D)(t) == MAT-17 chi(t) exactly on {n_chi} points "
               f"(cases A and B), band 0: green")

    # ------------------------------------------------------ Part A, MAT-17 step 2
    t_a0 = time.perf_counter()
    msh, ctb, ft = g._ball(g.RES)
    n_glob = msh.topology.index_map(3).size_global
    rows = {}
    exacts = {}
    apis = {}
    for lab, P, D, dt, t_end, on_time, rec_print, rec_g, rec_it, rec_bal, rec_s0 in (
        ("A", g.P_A, g.D_A, 150.0, g.T_END_A, g.ON_TIME_A, g.MAT17_A_PRINT,
         REC["A_g"], REC["A_it"], REC["A_bal"], REC["A_s0"]),
        ("B", g.P_B, g.D_B, g.DT_B, g.T_END_B, g.ON_TIME_B, g.MAT17_B_PRINT,
         REC["B_g"], REC["B_it"], REC["B_bal"], REC["B_s0"]),
    ):
        n = g._n(t_end, dt)
        t0 = time.perf_counter()
        api = g._duty_study(dt, n, DutyCycle(P, D))
        t_api = time.perf_counter() - t0
        t0 = time.perf_counter()
        direct = g._run17(msh, ft, dt, n, 0.5, 2, schedule=g._chi(P, D))
        t_dir = time.perf_counter() - t0
        ca, cb, geo_eq = g._same_mesh(api.mesh, msh, comm)
        mult_eq = np.array_equal(np.asarray(api.source_multipliers),
                                 np.asarray(direct.source_multipliers))
        n_on, n_off = g._census(api, dt, 2, on_time, f"EX-123 (S0) case {lab} API")
        T_api = g._on(msh, api.temperature)
        it = g._theta_owned_rel(T_api, direct.temperature, comm)
        t0 = time.perf_counter()
        ex = g._exact_duty(msh, t_end, P, D)
        t_ex = time.perf_counter() - t0
        e = m13._interior_error(msh, T_api, ex)
        e_dir = m13._interior_error(msh, direct.temperature, ex)
        rows[lab] = (e, n, dt)
        exacts[lab] = ex
        apis[lab] = (api, T_api)
        all_T += [api.temperature, direct.temperature]
        _say(comm, f"[A] (UT) case {lab}: cells {ca} / {cb} (direct {n_glob}); geometry equal "
                   f"{geo_eq}; multipliers equal {mult_eq} ({len(api.source_multipliers)} (half-)steps)")
        _say(comm, f"[A] (S0) case {lab}: on/off {n_on}/{n_off} (record {rec_s0[0]}/{rec_s0[1]})")
        _say(comm, f"[A] (IT) case {lab}: {it:.3e} (band {g.L0_BAND:.0e}; record "
                   f"{REC['A_it'] if lab == 'A' else REC['B_it']:.3e})")
        _say(comm, f"[A] (G) case {lab} N={n} dt={dt:.4f}: interior_rel_L2 {e:.10e} (band "
                   f"{g.INTERIOR_BAND:.0e}); direct {e_dir:.10e}; balance API "
                   f"{api.max_balance_residual:.3e} direct {direct.max_balance_residual:.3e} "
                   f"(record {rec_bal:.3e})")
        _say(comm, f"[A] (R) case {lab}: {e:.4e} vs record {rec_print} (string equal "
                   f"{f'{e:.4e}' == rec_print}); 10-digit record {rec_g:.10e}, rel diff "
                   f"{abs(e - rec_g) / rec_g:.2e} (printed only, predicted <= 1e-9)")
        _say(comm, f"[A] wall times case {lab}: API solve+mesh {t_api:.1f}s, direct {t_dir:.1f}s, "
                   f"comparand {t_ex:.1f}s")
        assert ca == cb == n_glob, f"(UT) cells case {lab}"
        if lab == "A":
            assert ca == REC["A_cells"], "(UT) cells vs record"
        assert geo_eq, f"(UT) geometry case {lab}"
        assert mult_eq, f"(UT) multipliers case {lab}"
        assert it <= g.L0_BAND, f"(IT) case {lab} {it:.3e}"
        assert api.max_balance_residual <= g.BALANCE_BAND, f"(G) balance case {lab}"
        assert e <= g.INTERIOR_BAND, f"(G) case {lab} {e:.4e}"
        assert abs(e - float(rec_print)) / float(rec_print) <= 1e-3, f"(R) case {lab} {e:.4e}"

    # (c1) the schedule dropped, case A
    e_a, n_a, dt_a = rows["A"]
    c1 = g._duty_study(dt_a, n_a, None)
    all_T.append(c1.temperature)
    T_c1 = g._on(msh, c1.temperature)
    e_c1 = m13._interior_error(msh, T_c1, exacts["A"])
    _say(comm, f"[control] (c1) no DutyCycle, case A: interior miss {e_c1:.4e} (record "
               f"{REC['c1']:.4e}; {e_c1 / g.INTERIOR_BAND:.0f}x the band); source_multipliers "
               f"{c1.source_multipliers!r}  ASSERTED > {g.INTERIOR_BAND:.0e}")
    assert c1.source_multipliers == []
    assert e_c1 > g.INTERIOR_BAND, f"(c1) did not miss: {e_c1:.4e}"
    t_partA = time.perf_counter() - t_a0
    _say(comm, f"[A] Part A leg {t_partA:.1f}s")

    # ------------------------------------------------------ Part B, MAT-18 2(b)
    t_b0 = time.perf_counter()
    mshw, ct, ftw, nglob, k_f, rc_f, q_f, theta = g._fine()
    assert nglob <= g.CELL_GUARD
    t0 = time.perf_counter()
    apiw = g._wall_study(g.RC2)
    t_apiw = time.perf_counter() - t0
    directw, t_dirw = g._run18(mshw, ftw, k=k_f, rho_c=rc_f, q=q_f, n=g.N_CN, theta=0.5, startup=2)
    caw, cbw, geow = g._same_mesh(apiw.mesh, mshw, comm)
    dk = g._dg0_owned_diff(apiw.k, k_f, comm)
    drc = g._dg0_owned_diff(apiw.rho_c, rc_f, comm)
    dq = g._dg0_owned_diff(apiw.q, q_f, comm)
    _say(comm, f"[B] (UT) wall: cells {caw} / {cbw} (record {g.MAT18_T1_CELLS}); geometry equal "
               f"{geow}; DG0 k max|d| {dk[0]!r} ({dk[1]} differ), rho_c {drc[0]!r} ({drc[1]}), "
               f"q {dq[0]!r} ({dq[1]}) of {dk[2]} owned")
    # (c3) comparator planted: one owned DG0 entry on rank 0 (gate :311-319)
    k_p = fem.Function(apiw.k.function_space)
    k_p.x.array[:] = apiw.k.x.array
    sp = float(np.spacing(np.real(k_p.x.array[0]))) if comm.rank == 0 else 0.0
    if comm.rank == 0:
        k_p.x.array[0] = k_p.x.array[0] + sp
    sp = comm.bcast(sp, root=0)
    dpl = g._dg0_owned_diff(k_p, k_f, comm)
    _say(comm, f"[control] (c3) planted spacing {sp!r} (record {REC['c3_sp']!r}): read back "
               f"max|d| {dpl[0]!r}, differing {dpl[1]}  ASSERTED equal and 1")
    T_apiw = g._on(mshw, apiw.temperature)
    itw = g._theta_owned_rel(T_apiw, directw.temperature, comm)
    e_all, e_in, e_sh = g._errors(mshw, ct, T_apiw, theta)
    _, e_dirw, _ = g._errors(mshw, ct, directw.temperature, theta)
    all_T += [apiw.temperature, directw.temperature]
    _say(comm, f"[B] (IT) wall: {itw:.3e} (band {g.L0_BAND:.0e}; record {REC['wall_it']:.3e})")
    _say(comm, f"[B] (G) (T1) CN+Rannacher N={g.N_CN}: interior {e_in:.10e} (band "
               f"{g.INTERIOR_BAND:.0e}); whole {e_all:.4e} (record {REC['wall_whole']:.4e}), shell "
               f"{e_sh:.4e} (record {REC['wall_shell']:.4e}); direct {e_dirw:.10e}; balance API "
               f"{apiw.max_balance_residual:.3e} direct {directw.max_balance_residual:.3e} "
               f"(record {REC['wall_bal']:.3e})")
    _say(comm, f"[B] (R) (T1): {e_in:.4e} vs record {g.MAT18_T1_PRINT} (string equal "
               f"{f'{e_in:.4e}' == g.MAT18_T1_PRINT}); 10-digit record {REC['wall_g']:.10e}, rel "
               f"diff {abs(e_in - REC['wall_g']) / REC['wall_g']:.2e} (printed only, predicted <= 1e-9)")
    c2 = g._wall_study(g.RC1)
    all_T.append(c2.temperature)
    _, e_c2, _ = g._errors(mshw, ct, g._on(mshw, c2.temperature), theta)
    _say(comm, f"[control] (c2) wall rho_c = RC1: interior miss {e_c2:.4e} (record "
               f"{g.MAT18_C1_UNIFORM_RHOC:.4e}; "
               f"{e_c2 / g.INTERIOR_BAND:.1f}x the band)  ASSERTED > {g.INTERIOR_BAND:.0e}; "
               f"balance {c2.max_balance_residual:.3e}")
    _say(comm, f"[B] wall times: API solve+mesh {t_apiw:.1f}s, direct {t_dirw:.1f}s; Part B leg "
               f"{time.perf_counter() - t_b0:.1f}s")
    assert caw == cbw == g.MAT18_T1_CELLS == REC["wall_cells"], "(UT) cells"
    assert geow, "(UT) geometry"
    assert dk[0] == 0.0 and drc[0] == 0.0 and dq[0] == 0.0, "(UT) DG0 coefficients"
    assert dpl[0] == sp and dpl[1] == 1, "(c3) planted spacing"
    assert itw <= g.L0_BAND, f"(IT) wall {itw:.3e}"
    assert apiw.max_balance_residual <= g.BALANCE_BAND, "(G) wall balance"
    assert c2.max_balance_residual <= g.BALANCE_BAND, "(G) balance on (c2)"
    assert e_in <= g.INTERIOR_BAND, f"(G) (T1) {e_in:.4e}"
    assert abs(e_in - float(g.MAT18_T1_PRINT)) / float(g.MAT18_T1_PRINT) <= 1e-3, "(R) (T1)"
    assert e_c2 > g.INTERIOR_BAND, f"(c2) did not miss: {e_c2:.4e}"

    # ------------------------------------------- (X) duty on the walled ball
    t_x0 = time.perf_counter()
    geo_x = ThermalBallGeometry(g.A_, g.RES, wall=ContainerWall(g.B_ - g.A_, g.K2, g.RC2))
    comp = Study(geo_x, None, UniformHeatSource(g.Q1)).temperature(
        RobinCooling(g.K1, g.H, g.T_INF), transient=TransientSpec(150.0, 94, g.RC1),
        duty=DutyCycle(g.P_A, g.D_A))
    all_T.append(comp.temperature)
    cx = comp.mesh.topology.index_map(3).size_global
    n_on_x, n_off_x = g._census(comp, 150.0, 2, g.ON_TIME_A, "EX-123 (X) walled duty")
    mult_x = np.array_equal(np.asarray(comp.source_multipliers),
                            np.asarray(apis["A"][0].source_multipliers))
    _say(comm, f"[X] composition: cells {cx} (record {g.MAT18_T1_CELLS}); multipliers equal to "
               f"case A's {mult_x}; balance {comp.max_balance_residual:.3e} (band "
               f"{g.BALANCE_BAND:.0e}; predicted ~1e-14); wall time {time.perf_counter() - t_x0:.1f}s")
    assert cx == g.MAT18_T1_CELLS, "(X) cells"
    assert mult_x, "(X) multipliers"
    assert comp.max_balance_residual <= g.BALANCE_BAND, "(X) balance"

    pts = np.array([[0.0, 0.0, 0.0], [0.05, 0.0, 0.0], [0.09, 0.0, 0.0]])
    th_w = _theta_at(comp.temperature, pts)
    th_b = _theta_at(apis["A"][0].temperature, pts)
    ratio = th_w / th_b
    _say(comm, f"[X] theta_h(t_end = {g.T_END_A:.0f} s) at x = 0, 0.05, 0.09 m: walled "
               f"{np.array2string(th_w, precision=4)}, bare {np.array2string(th_b, precision=4)}, "
               f"ratio {np.array2string(ratio, precision=4)}")
    _say(comm, f"[X] predicted centre ratio in [0.99, 1.03]: {0.99 <= ratio[0] <= 1.03}; r = 0.09 "
               f"ratio above centre ratio and in [1.03, 1.5]: "
               f"{bool(ratio[2] > ratio[0] and 1.03 <= ratio[2] <= 1.5)} (printed only; no closed "
               f"form, no accuracy claimed)")

    # ------------------------------------------------------------------ (Im)
    im = _max_im(comm, all_T)
    _say(comm, f"[complex] max |Im T_h| over {len(all_T)} solves = {im:.3e}  ASSERTED == 0")
    assert im == 0.0, f"max |Im T_h| = {im:.3e}"

    # ----------------------------------------------------------------- outputs
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    V1 = apis["A"][1].function_space
    fa = {}
    for lab in ("A", "B"):
        T_on = apis[lab][1]
        th = _p1_minus_tinf(T_on, f"theta_api_{lab}")
        ex1 = _p1_of(exacts[lab], V1, f"theta_exact_{lab}")
        fa[f"theta_api_{lab}"] = th
        fa[f"theta_exact_{lab}"] = ex1
        fa[f"theta_err_{lab}"] = _diff(th, ex1, V1, f"theta_err_{lab}")
    fa["theta_no_duty_A"] = _p1_minus_tinf(T_c1, "theta_no_duty_A")
    xd_duty = OUTPUT_DIR / f"{BASENAME}_duty_combined"
    write_xdmf_with_tags(xd_duty, msh, ctb, fa, comm=comm)

    V1w = T_apiw.function_space
    fw = {}
    th_t1 = _p1_minus_tinf(T_apiw, "theta_wall_T1")
    th_ser = _p1_of(theta, V1w, "theta_wall_series")
    fw["theta_wall_T1"] = th_t1
    fw["theta_wall_series"] = th_ser
    fw["theta_wall_err"] = _diff(th_t1, th_ser, V1w, "theta_wall_err")
    cc, cd, geoc = g._same_mesh(comp.mesh, mshw, comm)
    if cc == cd and geoc:
        fw["theta_wall_duty_tend"] = _p1_minus_tinf(g._on(mshw, comp.temperature),
                                                    "theta_wall_duty_tend")
    for nm, src in (("k", apiw.k), ("rho_c", apiw.rho_c), ("q", apiw.q)):
        f = fem.Function(src.function_space, name=nm)
        f.x.array[:] = np.real(src.x.array)
        f.x.scatter_forward()
        fw[nm] = f
    xd_wall = OUTPUT_DIR / f"{BASENAME}_wall_combined"
    write_xdmf_with_tags(xd_wall, mshw, ct, fw, comm=comm)
    adopt_host_ownership(OUTPUT_DIR, comm=comm)
    _say(comm, f"[output] {xd_duty}.xdmf -- P1 fields on the bare ball ({len(fa)} fields)")
    _say(comm, f"[output] {xd_wall}.xdmf -- P1 and DG0 fields on the walled ball ({len(fw)} fields)")

    write_setup_figure(
        apiw.mesh,
        apiw.cell_tags,
        FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names={
            1: f"gel core (k {g.K1} W/m/K, rho c {g.RC1:.3g})",
            2: f"container wall (k {g.K2} W/m/K, rho c {g.RC2:.3g})",
        },
        translucent_tags=(2,),
        slice_normal=(0.0, 1.0, 0.0),
        slice_origin=(0.0, 0.0, 0.0),
        title="mri:25 -- Study transient API: walled gel ball (MAT-18 2(b)) + duty cycle (MAT-17 step 2), no EM solve",
        window_size=(1400, 630),
        comm=comm,
    )

    _say(comm, f"\n[timing] total {time.perf_counter() - started:.1f} s at -n {comm.size}; "
               f"summed ru_maxrss {g._rss_gib(comm):.2f} GiB")
    _say(comm, "\nAll assertions hold.")


if __name__ == "__main__":
    main()
