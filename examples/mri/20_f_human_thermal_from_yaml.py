#!/usr/bin/env python3
"""mri:20 (EX-114) -- `femem run` thermal on the tuned copper F-human from YAML.

One EM solve (``run_study_file`` on the tracked ``f_human_tuned_copper_thermal.yaml``,
read in place) answers the file's thermal knobs -- ``source_scale``, ``h``, ``balls`` --
and the temperature field the API does not return is rebuilt from the API's own
thermal set-up (copied verbatim) and written for ParaView.

Every band and helper is imported from
``tests/package/test_wf1_yaml_f_human_thermal_solve.py`` (``g``); the record digits
live only in the example-local ``REC`` dict, each beside its gate-log line.  New
values are printed at <= 6 significant digits (the private-leak guard).

Run (heavy, complex build, 8 ranks):  ./run_examples.sh -e mri:20 -n 8 -t 590
Setup figure: add -e FEM_EM_SETUP_FIGURES=1 after ``exec -T``.
"""

from __future__ import annotations

import dataclasses
import json
import sys
import time
from pathlib import Path

import numpy as np
import ufl
from dolfinx import default_scalar_type, fem
from dolfinx.mesh import compute_midpoints
from mpi4py import MPI

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from fem_em_solver.io.paraview_utils import adopt_host_ownership, write_xdmf_with_tags  # noqa: E402
from fem_em_solver.post.setup_figure import write_setup_figure  # noqa: E402
from fem_em_solver.thermal.transient import solve_transient_conduction  # noqa: E402

import tests.package.test_wf1_yaml_f_human_thermal_solve as g  # noqa: E402
from tests.mesh.test_birdcage_port_sheets import PORT_LOWER, PORT_UPPER  # noqa: E402

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"  # committed, EX-57
BASENAME = "mri_20_f_human_thermal_from_yaml"
GATE_LOG = "20261014T020254Z_WF-1-s10b-W1.log"
H_INPUTS = (1.0, 100.0)  # example-local inputs beside the file's h (cited, not records)

# ---- the gate's records, each beside its log:line (compared at printed digits) ----
REC = {
    "cells": 369387,                      # GATE_LOG:6283
    "v_sub": "2.210693836095e-03",        # :6285
    "n_sub": 3306,                        # :6285
    "source_power": "2.821343985102460e-01",  # :6286
    "eta": "0.282134",                    # :6288
    "theta_tstar": "4.926110849651e-02",  # :6297
    "six_min": "4.277005459539e-02",      # :6308
    "vol_six_min": "3.070252443e-02",     # :6310
    "n1": "4.275342672e-02",              # :6319
    "n3": "9.852221699e-02",              # :6321
}


def _digits(s: str) -> int:
    return len(s.split("e")[0].replace(".", "").replace("-", "").replace("+", ""))


def main() -> None:
    comm = MPI.COMM_WORLD
    rank0 = comm.rank == 0
    started = time.perf_counter()
    if not np.issubdtype(np.dtype(default_scalar_type), np.complexfloating):
        raise RuntimeError("this example needs the complex DolfinX build "
                           "(source /usr/local/bin/dolfinx-complex-mode)")

    def p0(msg):
        if rank0:
            print(msg, flush=True)

    def rec_check(label, got, want):
        """A covered record is compared at max(10^-(d-2), REUSE_REPRODUCTION_RTOL),
        d its printed significant digits (review ruling 2026-10-14 10:30); exact
        (integer) records stay exact."""
        w = float(want) if isinstance(want, str) else float(want)
        tol = (max(10.0 ** (-(_digits(want) - 2)), g.REUSE_REPRODUCTION_RTOL)
               if isinstance(want, str) else 0.0)
        rel = abs(got - w) / abs(w)
        p0(f"[record] {label}: matches {GATE_LOG} to rel {rel:.3e} (limit {tol:.0e})")
        assert rel <= tol, f"{label}: {got!r} vs {want!r} (rel {rel:.3e})"

    p0("=" * 78)
    p0("EX-114 -- femem run thermal on the tuned copper F-human from YAML (mri:20)")
    p0("=" * 78)
    p0(f"[records] every covered row is beside {GATE_LOG}:6282-6321")

    # ---- inputs and structural asserts (the gate's pattern, :105-138) ---------------
    series = None
    if rank0:
        assert g.SERIES_JSON and (g.REPO / g.SERIES_JSON).is_file(), \
            f"{g.SERIES_JSON} absent: BLOCKED (rule (d))"
        series = json.loads((g.REPO / g.SERIES_JSON).read_text())
    series = comm.bcast(series, root=0)
    assert comm.size == g.RECORD_RANK_WIDTH, f"records are -n {g.RECORD_RANK_WIDTH}; got {comm.size}"
    _spec, kind, _outputs, fh = g._load(g.FHT_FILE)
    assert kind == "f_human_birdcage", kind
    assert fh.cooling.source_scale == float(series["G_per_w"])
    assert fh.window_s == g.W_B2
    assert series["ranks"] == g.RECORD_RANK_WIDTH and series["cells"] == g.HOLE_CELL_RECORD
    G = float(series["G_per_w"])

    # ---- the one EM solve, through the file's own entry point -----------------------
    calls = {"n": 0}
    real = g.study_mod._study_solve_tuned_quadrature

    def counting(*args, **kwargs):
        calls["n"] += 1
        return real(*args, **kwargs)

    g.study_mod._study_solve_tuned_quadrature = counting
    try:
        comm.Barrier()
        t0 = time.perf_counter()
        run = g.run_study_file(g.FHT_FILE, output_dir=OUTPUT_DIR / "mri_20_run", comm=comm)
        comm.Barrier()
        elapsed = time.perf_counter() - t0
    finally:
        g.study_mod._study_solve_tuned_quadrature = real
    res, summary, th = run.result, run.summary, run.thermal
    p0(f"[run] run_study_file {elapsed:.1f} s at -n {comm.size} (build "
       f"{res.times.get('build', 0.0):.1f} s, node pair {res.times.get('node_pair', 0.0):.1f} s); "
       f"solve calls {calls['n']}")

    def once(where):
        assert calls["n"] == 1, f"(Once) {where}: {calls['n']} EM solves"

    once("after run_study_file")

    # EX-57 setup figure -- right after the built mesh is in hand.
    c_star = [float(x) for x in series["c_star_m"]]
    names = {2: "air", g.PHANTOM_CELL_TAG: "gel phantom"}
    for n, label in ((17, "feed P17 (0 deg)"), (21, "feed P21 (90 deg)")):
        names[PORT_LOWER + n] = label
        names[PORT_UPPER + n] = label
    write_setup_figure(
        res.mesh, res.cell_tags, FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names=names, hide_tags=(2,), translucent_tags=(g.PHANTOM_CELL_TAG,),
        slice_normal=(0.0, 0.0, 1.0), slice_origin=(0.0, 0.0, c_star[2]),
        title="mri:20 -- femem run thermal on the tuned copper F-human from YAML: gel phantom, "
              "10 g ball at c*_theta (WF-1 step 10b)",
        window_size=(1400, 630), comm=comm,
    )

    def _red(form):
        return float(np.real(comm.allreduce(fem.assemble_scalar(form), op=MPI.SUM)))

    # ---- (A0), (V), (B), (E) --------------------------------------------------------
    tp = summary["temperature"]
    cell_rel = summary["cells"] / g.HOLE_CELL_RECORD - 1.0
    sub, _ids = g.study_mod._tag_submesh(res.mesh, res.cell_tags, run.study.geometry.phantom_tag)
    one = fem.Constant(sub, default_scalar_type(1.0))
    v_sub = _red(fem.form(one * ufl.dx(domain=sub)))
    r_v = abs(v_sub / g.V_PHANTOM_H9P_RECORD_M3 - 1.0)
    f = res.fields["ccw"]
    tagged = g.mean_sar(f.e_complex, sigma=f.sigma_field, rho=g.PHANTOM_RHO_KG_PER_M3,
                        cell_tags=res.cell_tags, subdomain_ids=g.PHANTOM_CELL_TAG, comm=comm)
    p_ph = float(np.real(tagged["dissipated_power_w"]))
    r_b = abs(tp["source_power_w"] - G * p_ph) / (G * p_ph)
    eta = p_ph * G
    p0(f"[anchor] (A0) cells {summary['cells']} vs record {g.HOLE_CELL_RECORD}, rel {cell_rel:.3e} "
       f"(band {g.CELL_COUNT_BAND}, ASSERTED); {GATE_LOG}:6283")
    p0(f"[anchor] (V) V_sub {v_sub:.12e} m^3: rel {r_v:.3e} (band {g.V_BAND:g}, ASSERTED); "
       f"n_sub {th.n_sub} vs JSON {series['n_sub']}; {GATE_LOG}:6285")
    p0(f"[anchor] (B) source_power_w {tp['source_power_w']:.15e} W: rel {r_b:.3e} "
       f"(band {g.B_BAND:g}, ASSERTED); {GATE_LOG}:6286")
    p0(f"[anchor] (E) max balance residual {tp['max_balance_residual']:.3e} "
       f"(band {g.E_BAND:g}, ASSERTED); {GATE_LOG}:6287")
    p0(f"[anchor] eta = P_phantom / P_acc,ccw {eta:.6f}; {GATE_LOG}:6288")
    assert abs(cell_rel) <= g.CELL_COUNT_BAND, cell_rel
    assert summary["cells"] == int(res.cells) == th.cells
    assert run.result is run.study.result
    assert r_v <= g.V_BAND, r_v
    assert th.n_sub == series["n_sub"]
    assert r_b <= g.B_BAND, r_b
    assert tp["max_balance_residual"] <= g.E_BAND, tp["max_balance_residual"]
    rec_check("cells", float(summary["cells"]), float(REC["cells"]))
    rec_check("V_sub", v_sub, REC["v_sub"])
    assert th.n_sub == REC["n_sub"]
    rec_check("source power", tp["source_power_w"], REC["source_power"])
    rec_check("eta", eta, REC["eta"])
    once("after (A0)/(V)/(B)/(E)")

    # ---- (T*), (Ser), (6min) --------------------------------------------------------
    assert tp["sense"] == "ccw" and len(tp["balls"]) == 1
    ball = tp["balls"][0]
    th_star = ball["theta_k"][-1]
    th_json = float(series["theta_ball10g_k"][-1])
    r_t = abs(th_star / th_json - 1.0)
    t_sum = np.array(tp["times_s"])
    t_json = np.array(series["times_s"])
    assert t_sum.shape == t_json.shape
    d_t = float(np.max(np.abs(t_sum - t_json)))
    b_sum = np.array(ball["theta_k"])
    b_json = np.array(series["theta_ball10g_k"])
    m_sum = np.array(tp["volume_mean_k"])
    m_json = np.array(series["theta_mean_k"])
    r_ball = float(np.max(np.abs(b_sum - b_json)) / np.max(np.abs(b_json)))
    r_mean = float(np.max(np.abs(m_sum - m_json)) / np.max(np.abs(m_json)))
    _m, peak_json, tmax_json = g.running_average(t_json, b_json, window=g.W_B2)
    six = tp["six_minute"]
    peak, tmax = six["balls"][0]["peak_k"], six["balls"][0]["t_max_s"]
    r_6 = abs(peak / float(peak_json) - 1.0)
    _m2, peak_own, _t_own = g.running_average(t_sum, b_sum, window=g.W_B2)
    r_own = abs(peak - float(peak_own)) / abs(float(peak_own))
    vol_peak = six["volume_mean"]["peak_k"]
    p0(f"[anchor] (T*) theta-bar(c*_theta, t* = {tp['times_s'][-1]:.1f} s) {th_star:.12e} K: rel "
       f"{r_t:.3e} (band {g.REUSE_REPRODUCTION_RTOL:g}, ASSERTED); {GATE_LOG}:6297")
    p0(f"[anchor] (Ser) {t_sum.size} samples; max|dt| {d_t:.3e} s (<= 1e-12 ASSERTED); ball "
       f"{r_ball:.3e}, volume-mean {r_mean:.3e} (band {g.REUSE_REPRODUCTION_RTOL:g}, ASSERTED); "
       f"{GATE_LOG}:6307")
    p0(f"[anchor] (6min) window {six['window_s']} s: peak {peak:.12e} K at {tmax} s: rel {r_6:.3e} "
       f"(band {g.REUSE_REPRODUCTION_RTOL:g}, ASSERTED; t_max exact ASSERTED); API vs operator "
       f"{r_own:.3e} (band {g.AVG_BAR:g}); {GATE_LOG}:6308-6309")
    p0(f"[anchor] volume-mean six-minute peak {vol_peak:.9e} K at "
       f"{six['volume_mean']['t_max_s']} s (printed); {GATE_LOG}:6310")
    assert r_t <= g.REUSE_REPRODUCTION_RTOL, r_t
    assert d_t <= 1e-12, d_t
    assert r_ball <= g.REUSE_REPRODUCTION_RTOL, r_ball
    assert r_mean <= g.REUSE_REPRODUCTION_RTOL, r_mean
    assert six["window_s"] == g.W_B2
    assert r_6 <= g.REUSE_REPRODUCTION_RTOL, r_6
    assert tmax == float(tmax_json), (tmax, tmax_json)
    assert r_own <= g.AVG_BAR, r_own
    rec_check("theta-bar(t*)", th_star, REC["theta_tstar"])
    rec_check("six-minute peak", peak, REC["six_min"])
    rec_check("volume-mean six-minute peak", vol_peak, REC["vol_six_min"])
    once("after (T*)/(Ser)/(6min)")

    # ---- the file's knobs, answered without a second EM solve ------------------------
    def temp(cooling, centres=None):
        balls = fh.balls if centres is None else {"mass_g": fh.balls["mass_g"], "centres": centres}
        return run.study.temperature(cooling, transient=fh.transient, balls=balls,
                                     window_s=fh.window_s, sense="ccw", comm=comm)

    # (Lin) 2G doubles every sample.
    t0 = time.perf_counter()
    n3 = temp(dataclasses.replace(fh.cooling, source_scale=2.0 * fh.cooling.source_scale))
    t_n3 = time.perf_counter() - t0
    lin_ball = abs(n3.balls[0]["theta_k"][-1] / (2.0 * th_star) - 1.0)
    lin_peak = abs(n3.six_minute["balls"][0]["peak_k"] / (2.0 * peak) - 1.0)
    lin_mean = abs(n3.volume_mean_k[-1] / (2.0 * m_sum[-1]) - 1.0)
    r_n3 = abs(n3.balls[0]["theta_k"][-1] / th_json - 1.0)
    p0(f"[Lin] 2G / (2 x G): ball {lin_ball:.3e}, six-minute peak {lin_peak:.3e}, volume mean "
       f"{lin_mean:.3e} (band {g.REUSE_REPRODUCTION_RTOL:g}, ASSERTED); {t_n3:.1f} s")
    assert max(lin_ball, lin_peak, lin_mean) <= g.REUSE_REPRODUCTION_RTOL
    once("after (Lin)")

    # (Balls) two balls in one call; ball 1 is the mirror image.
    m_c = [float(x) for x in g._mirror_xy(np.array([c_star], dtype=float),
                                          float(res.az_mirror))[0]]
    t0 = time.perf_counter()
    nb = temp(fh.cooling, centres=[c_star, m_c])
    t_nb = time.perf_counter() - t0
    b0 = np.array(nb.balls[0]["theta_k"])
    r_b0 = float(np.max(np.abs(b0 - b_sum)) / np.max(np.abs(b_sum)))
    n1_val = nb.balls[1]["theta_k"][-1]
    r_n1 = abs(n1_val / th_json - 1.0)
    p0(f"[Balls] ball 0 vs the one-ball series: max|d|/max {r_b0:.3e} (band "
       f"{g.REUSE_REPRODUCTION_RTOL:g}, ASSERTED); c*_theta = ({c_star[0]:.6g}, {c_star[1]:.6g}, "
       f"{c_star[2]:.6g}) m, M(c*_theta) = ({m_c[0]:.6g}, {m_c[1]:.6g}, {m_c[2]:.6g}) m; {t_nb:.1f} s")
    p0(f"[Balls] (N1) mirror ball theta-bar(t*) {n1_val:.9e} K, miss {r_n1:.6e} "
       f"(> {g.REUSE_REPRODUCTION_RTOL:g} ASSERTED); {GATE_LOG}:6319")
    assert r_b0 <= g.REUSE_REPRODUCTION_RTOL
    rec_check("(N1) mirror ball", n1_val, REC["n1"])
    assert r_n1 > g.REUSE_REPRODUCTION_RTOL, f"(N1) blind comparator {r_n1:.3e}"
    once("after (Balls)")

    # (N2) the window one sample short; (N3) from the 2G call.
    dt = float(series["dt_s"])
    _m3, peak_short, _t3 = g.running_average(tp["times_s"], tp["balls"][0]["theta_k"],
                                             window=g.W_B2 - dt)
    r_n2 = abs(float(peak_short) / float(peak_json) - 1.0)
    p0(f"[control] (N2) window {g.W_B2 - dt} s peak misses by {r_n2:.6e} "
       f"(> {g.REUSE_REPRODUCTION_RTOL:g} ASSERTED); {GATE_LOG}:6320")
    p0(f"[control] (N3) 2G theta-bar(t*) {n3.balls[0]['theta_k'][-1]:.9e} K, miss {r_n3:.6f} "
       f"(> {g.D_CONTROL_FLOOR:g} ASSERTED); {GATE_LOG}:6321")
    assert r_n2 > g.REUSE_REPRODUCTION_RTOL, f"(N2) blind comparator {r_n2:.3e}"
    assert r_n3 > g.D_CONTROL_FLOOR, f"(N3) blind comparator {r_n3:.3e}"
    rec_check("(N3) 2G ball", n3.balls[0]["theta_k"][-1], REC["n3"])

    # h table (printed only).
    rho_c = float(fh.transient.rho_c_j_per_m3_k)
    t_star = float(tp["times_s"][-1])
    times_arr = np.array(tp["times_s"])
    p_src = float(tp["source_power_w"])
    adiabatic = p_src * times_arr / (rho_c * v_sub)
    h_runs = {}
    for h_val in H_INPUTS:
        h_runs[h_val] = temp(dataclasses.replace(fh.cooling, h_w_per_m2_k=h_val))
    h_runs[float(fh.cooling.h_w_per_m2_k)] = None  # the file's own run (summary)
    p0(f"[h table] adiabatic line P_src t*/(rho c V_sub) at t* = {adiabatic[-1]:.6g} K "
       f"(printed; predicted 0.03829)")
    for h_val in sorted(h_runs):
        if h_runs[h_val] is None:
            bt, pk, tm, vm = th_star, peak, tmax, m_sum[-1]
        else:
            r_h = h_runs[h_val]
            bt = r_h.balls[0]["theta_k"][-1]
            pk = r_h.six_minute["balls"][0]["peak_k"]
            tm = r_h.six_minute["balls"][0]["t_max_s"]
            vm = r_h.volume_mean_k[-1]
        p0(f"[h table] h = {h_val:g} W/m2K: ball theta-bar(t*) {bt:.6g} K, six-minute peak "
           f"{pk:.6g} K at {tm:.1f} s, volume mean(t*) {vm:.6g} K, volume mean / adiabatic "
           f"{vm / adiabatic[-1]:.6g}")
    once("after the h table")

    # ---- the field: the API's thermal set-up, copied verbatim ------------------------
    sm = g.study_mod
    cooling, transient = fh.cooling, fh.transient
    res_c = sm._tuned_result(run.study, comm)
    rcomm = res_c.comm
    tag = run.study.geometry.phantom_tag
    sub, parent_ids = sm._tag_submesh(res_c.mesh, res_c.cell_tags, tag)
    ft = sm._exterior_robin_tags(sub)
    fc = res_c.fields["ccw"]
    Wp, avg = sm._dg0_cell_average(
        res_c.mesh, res_c.cell_tags, 0.5 * fc.sigma_field * ufl.inner(fc.e_complex, fc.e_complex), tag
    )
    avg *= cooling.source_scale  # after the average (ruling 361-4)
    q = sm._dg0_to_submesh(Wp, avg, sub, parent_ids)

    def _red_c(form):
        return float(np.real(rcomm.allreduce(fem.assemble_scalar(form), op=MPI.SUM)))

    one_c = fem.Constant(sub, default_scalar_type(1.0))
    v_sub_c = _red_c(fem.form(one_c * ufl.dx(domain=sub)))
    source_power_c = _red_c(fem.form(q * ufl.dx(domain=sub)))
    n_steps = int(transient.n_steps)
    dt_c = float(transient.dt_s)
    t_inf = float(cooling.t_inf_c)
    t_init = t_inf if transient.t_init_c is None else float(transient.t_init_c)
    steps = tuple(range(n_steps + 1))
    trun = solve_transient_conduction(
        sub, ft, k=cooling.k_w_per_m_k, rho_c=transient.rho_c_j_per_m3_k, q=q,
        h=cooling.h_w_per_m2_k, t_inf=t_inf, t_init=t_init, robin_tag=sm._ROBIN_TAG,
        dt=dt_c, n_steps=n_steps, theta=float(transient.theta),
        startup_half_steps=int(transient.startup_half_steps), degree=int(cooling.degree),
        snapshot_steps=steps,
    )
    snaps = trun.snapshots
    thf = fem.Function(snaps[0].function_space)
    f_mean = fem.form(thf * ufl.dx(domain=sub))
    mean_c = []
    for n in steps:
        thf.x.array[:] = snaps[n].x.array - t_inf
        thf.x.scatter_forward()
        mean_c.append(_red_c(f_mean) / v_sub_c)
    mean_c = np.array(mean_c)
    r_copy = float(np.max(np.abs(mean_c - m_sum)) / np.max(np.abs(m_sum)))
    r_copy_p = abs(source_power_c - tp["source_power_w"]) / abs(tp["source_power_w"])
    p0(f"[Copy] copied thermal set-up volume-mean series vs the summary: max|d|/max {r_copy:.3e} "
       f"(band {g.REUSE_REPRODUCTION_RTOL:g}, ASSERTED); source power rel {r_copy_p:.3e} "
       f"(band {g.B_BAND:g}, ASSERTED)")
    assert r_copy <= g.REUSE_REPRODUCTION_RTOL, f"(Copy) {r_copy:.3e}"
    assert r_copy_p <= g.B_BAND, f"(Copy) power {r_copy_p:.3e}"
    # (Copy-c) the one-sample-shifted copy must miss.
    r_shift = float(np.max(np.abs(mean_c[1:] - m_sum[:-1])) / np.max(np.abs(m_sum)))
    p0(f"[control] (Copy-c) copied series shifted one sample misses by {r_shift:.6g} "
       f"(> {g.REUSE_REPRODUCTION_RTOL:g} ASSERTED; predicted 0.025-0.05)")
    assert r_shift > g.REUSE_REPRODUCTION_RTOL, f"(Copy-c) blind comparator {r_shift:.3e}"
    once("after the copied thermal solve")

    # ---- per-vertex maps on the phantom submesh --------------------------------------
    V = snaps[0].function_space
    bs = V.dofmap.index_map_bs
    n_own = V.dofmap.index_map.size_local * bs
    times_c = np.array([n * dt_c for n in steps])
    a_mat = np.empty((len(steps), n_own))
    for n in steps:
        a_mat[n, :] = np.real(snaps[n].x.array[:n_own]) - t_inf
    theta_t = a_mat[-1, :].copy()
    if n_own > 0:
        _mm, pk_vert, _tm = g.running_average(times_c, a_mat, window=fh.window_s)
    else:
        pk_vert = np.zeros(0)
    th_map = fem.Function(V, name="theta_tstar_K")
    th_map.x.array[:n_own] = theta_t
    th_map.x.scatter_forward()
    pk_map = fem.Function(V, name="theta_six_minute_peak_K")
    pk_map.x.array[:n_own] = pk_vert
    pk_map.x.scatter_forward()
    coords = V.tabulate_dof_coordinates()[:n_own]
    loc_max = float(np.max(pk_vert, initial=-np.inf))
    loc_pos = coords[int(np.argmax(pk_vert))].tolist() if n_own > 0 else [0.0, 0.0, 0.0]
    everyone = comm.allgather((loc_max, loc_pos))
    hot_val, hot_pos = max(everyone, key=lambda t: t[0])
    thf.x.array[:] = np.real(snaps[steps[-1]].x.array) - t_inf
    thf.x.scatter_forward()
    vm_vertex = _red_c(f_mean) / v_sub_c
    p0(f"[field] hottest vertex six-minute peak {hot_val:.6g} K at ({hot_pos[0]:.6g}, "
       f"{hot_pos[1]:.6g}, {hot_pos[2]:.6g}) m beside the 10 g ball's {peak:.6g} K "
       f"(ratio {hot_val / peak:.6g}; pointwise, no compliance claim; predicted above the ball)")
    p0(f"[field] volume mean of the per-vertex theta(t*) {vm_vertex:.6g} K vs the API's last "
       f"volume-mean sample {m_sum[-1]:.6g} K")

    tdim = sub.topology.dim
    dg0 = q.function_space
    imap = sub.topology.index_map(tdim)
    n_cells = imap.size_local + imap.num_ghosts
    mids = compute_midpoints(sub, tdim, np.arange(n_cells, dtype=np.int32))
    cell_dof = np.asarray(dg0.dofmap.list).reshape(-1)[:n_cells]
    radius = float(ball["radius_m"])

    def indicator(name, centre):
        fn = fem.Function(dg0, name=name)
        inside = np.linalg.norm(mids - np.array(centre)[None, :], axis=1) < radius
        fn.x.array[cell_dof[inside]] = 1.0
        fn.x.scatter_forward()
        return fn

    q_map = fem.Function(dg0, name="q_source_W_per_m3")
    q_map.x.array[:] = np.real(q.x.array)
    q_map.x.scatter_forward()
    fields = {
        "theta_tstar_K": th_map,
        "theta_six_minute_peak_K": pk_map,
        "q_source_W_per_m3": q_map,
        "ball10g_cstar": indicator("ball10g_cstar", c_star),
        "ball10g_mirror": indicator("ball10g_mirror", m_c),
    }

    # ---- outputs: CSV, then XDMF, then ownership ---------------------------------------
    OUTPUT_DIR.mkdir(exist_ok=True)
    csv_path = OUTPUT_DIR / f"{BASENAME}_series.csv"
    if rank0:
        with open(csv_path, "w") as fcsv:
            fcsv.write("t_s,ball_cstar_K,ball_mirror_K,volume_mean_K,volume_mean_h1_K,"
                       "volume_mean_h100_K,ball_h1_K,ball_h100_K,adiabatic_K\n")
            h1, h100 = h_runs[H_INPUTS[0]], h_runs[H_INPUTS[1]]
            for i, t in enumerate(tp["times_s"]):
                row = [t, b_sum[i], nb.balls[1]["theta_k"][i], m_sum[i],
                       h1.volume_mean_k[i], h100.volume_mean_k[i],
                       h1.balls[0]["theta_k"][i], h100.balls[0]["theta_k"][i], adiabatic[i]]
                fcsv.write(",".join(f"{v:.6g}" for v in row) + "\n")
    xdmf_path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_combined", sub, None, fields, comm=comm
    )
    adopt_host_ownership(OUTPUT_DIR, comm=comm)
    p0(f"[output] {xdmf_path} and {csv_path}")
    once("at the end")

    p0(f"[done] elapsed {time.perf_counter() - started:.1f} s at -n {comm.size}")
    p0("All assertions hold (EX-114): one EM solve, the file's thermal knobs, the field.")


if __name__ == "__main__":
    main()
