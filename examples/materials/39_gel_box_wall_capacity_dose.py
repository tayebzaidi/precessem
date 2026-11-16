"""Example (`EX-186`): the D = 0.25 thermal dose under the capacitive wall -- the closed form
(P70) and the FEM (P68).

`mat:38` shows the drift through a wall with heat capacity and states that the D = 0.25 dose
is not its claim.  This page's output is the CEM43 dose (``cem43``) under the duty cycle
D = 0.25 on item 682's 1D-reduced walled gel box (m = 4, 119 808 cells), with the wall's
heat capacity in the drift.

Two halves.

- **Closed form (numpy)**, P70's helpers imported by name
  (``test_mat21_wall_capacity_duty25_dose_closed_form.py`` and, through it, P64's ``Z``,
  ``WAITS``, ``_dose``, ``_duty_series``, ``_hists`` and item 823's ``_drift_cap_z``):
  the composite-drift dose, the massless-wall (``_drift_1d``) dose and the six-face h_eff
  dose at P5 / ``ADJ_PT`` / P0 for t_w in {0, T_C}.
- **FEM**, item 835's run body written out here (fixtures cannot be imported), its plain
  helpers imported by name: the reference run (theta_0 = 0) and the t_w = 0 run
  (theta_0 = g z), both under ``_chi25`` to t*, every snapshot read at P5 / ``ADJ_PT``.

No ``test_*`` name is imported.  No gate edit, no ``src/`` edit.

Asserted (all collected first, then asserted): (Q-d), (Win-d) with their controls and
(D50-cf) with the drop-run-0 control, at the imported ``DOSE_BAR`` / ``AVG_BAR`` (P70's
anchors); the vectorised ``_drift_cap_z`` against the scalar ``_drift_cap``; the control
that the massless-wall dose misses the composite dose at P5, t_w = 0, by more than
``DOSE_BAR``; FEM (I) on the t_w = 0 run against the composite drift <= ``INTERIOR_BAND``,
(I) on the reference run, (E) <= ``E_BAND`` on both runs, (S0), (Q-wf) at ``DOSE_BAR`` on all
four histories, the t_w = 0 / reference dose control at P5; max|Im theta| = 0.
Printed only: P70's table, the FEM dose at P5 beside the composite, the geometry versus
capacity split, I0h (the h_eff comparand).

The wall's k and rho_c are nominal PMMA, citation owed (W3).  One rung, one time step.
Not claimed: the six-wall box, a dose band, a recommended wait, a cited PMMA rho_c, a tissue
dose.

Run it::

    ./run_examples.sh -e mat:39 -n 2 -t 480

Output: ``examples/materials/paraview_output/
materials_39_gel_box_wall_capacity_dose_combined.xdmf`` (theta at t* on the t_w = 0 run and
its composite comparand).
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

from fem_em_solver.io.paraview_utils import (  # noqa: E402
    adopt_host_ownership,
    write_xdmf_with_tags,
)
from fem_em_solver.post.setup_figure import write_setup_figure  # noqa: E402
from fem_em_solver.thermal import solve_transient_conduction  # noqa: E402

from tests.validation.test_imp7_duty_cycle_fem import S0_BAND  # noqa: E402
from tests.validation.test_mat12_transient_sar_heating_sphere import E_BAND  # noqa: E402
from tests.validation.test_mat13_transient_perfusion_sphere import INTERIOR_BAND  # noqa: E402
from tests.validation.test_mat15_coil_driven_heating_phantom import (  # noqa: E402
    H_GEL,
    K_GEL,
    N_MAIN,
)
from tests.validation.test_mat15_finite_cylinder_closed_form import H, K, Q, T_INF  # noqa: E402
from tests.validation.test_mat15_finite_cylinder_fem import _max_imag, _rel  # noqa: E402
from tests.validation.test_mat17_duty_cycle_fem import _census  # noqa: E402
from tests.validation.test_mat21_box_closed_form import A_BOX, B_BOX, C_BOX, T_C  # noqa: E402
from tests.validation.test_mat21_box_dose_closed_form import (  # noqa: E402
    AVG_BAR,
    DOSE_BAR,
    T_GRID,
    W_B2,
)
from tests.validation.test_mat21_box_dt_duty125_dose_closed_form import (  # noqa: E402
    _q_check_g,
    _win_check_g,
)
from tests.validation.test_mat21_box_duty25_dose_closed_form import (  # noqa: E402
    DUTY25,
    SCHEDULE_25,
    T_SW25,
    _chi25,
    _sched_string,
)
from tests.validation.test_mat21_box_duty_closed_form import DUTY, T_SW  # noqa: E402
from tests.validation.test_mat21_initial_gradient_closed_form import G_GRAD  # noqa: E402
from tests.validation.test_mat21_pre_equilibration_closed_form import ADJ_PT, P5  # noqa: E402
from tests.validation.test_mat21_pre_equilibration_duty_dose_closed_form import LABS  # noqa: E402
from tests.validation.test_mat21_wall_capacity_drift_closed_form import (  # noqa: E402
    _drift_cap,
    _drift_cap_z,
)
from tests.validation.test_mat21_pre_equilibration_duty_dose_closed_form import PTS  # noqa: E402
from tests.validation.test_mat21_wall_capacity_duty25_dose_closed_form import (  # noqa: E402
    P65_FEM_D0_P5,
    _cap_hists,
    _z1d_hists,
)
from tests.validation.test_mat21_wall_closed_form import D_WALL, RHO_C_W  # noqa: E402
from tests.validation.test_mat21_wall_fem import (  # noqa: E402
    GEL,
    ROBIN,
    WALL,
    _materials,
    _robin_facets,
    _walled_box,
)
from tests.validation.test_mat21_wall_pre_equilibration_duty25_dose_closed_form import (  # noqa: E402
    WAITS,
    Z,
    _dose,
    _duty_series,
    _hists,
)
from tests.validation.test_mat21_wall_pre_equilibration_fem import (  # noqa: E402
    DT,
    PRE_BAR,
    _drift_z,
    _read,
    _series_z,
    _ZField,
)
from tests.validation.test_mat21_wall_transient_closed_form import T_STAR  # noqa: E402
from tests.validation.test_mat21_wall_transient_fem import RUNG, _rho_c_field  # noqa: E402

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"  # committed, EX-57
BASENAME = "materials_39_gel_box_wall_capacity_dose"
HIST_LABELS = ("P5", "ADJ_PT")
HIST_PTS = np.array((P5, ADJ_PT), dtype=float)


def _say(comm, msg):
    if comm.rank == 0:
        print(msg, flush=True)


def closed_form(comm, gates):
    """P70's anchors through its helpers (numpy, rank-identical)."""
    h = _hists()
    cap = _cap_hists()
    z1 = _z1d_hists()
    six = h["drift"]

    pre = 0.0
    # vectorised _drift_cap_z against the scalar _drift_cap (precondition), on PTS rows
    for wl, w in WAITS:
        for n in (0, 7, T_GRID.size - 1):
            for j, p in enumerate(PTS):
                s = _drift_cap(p, T_GRID[n] + w)
                pre = max(pre, abs(cap[wl][n, j] - s) / max(abs(s), 1e-300) if s
                          else abs(cap[wl][n, j]))
    _say(comm, f"[cf] vectorised _drift_cap_z vs scalar _drift_cap max rel {pre:.3e} "
               f"(asserted <= 1e-12)")
    gates["cf_precondition_cap_z"] = pre <= 1e-12

    d25 = h["25"]
    worst = {"q": 0.0, "w": 0.0, "rc": np.inf}
    for j, lab in enumerate(LABS):
        for wl, _w in WAITS:
            th = d25[:, j] + cap[wl][:, j]
            _dose_v, _quad, rel, _swap, miss = _q_check_g(T_INF + th, T_GRID)
            mx, t_mx, lo, hi, _ref, wrel = _win_check_g(th, T_GRID)
            gates[f"Qd_{lab}_{wl}"] = rel <= DOSE_BAR
            gates[f"Qd_swap_{lab}_{wl}"] = miss > DOSE_BAR
            gates[f"Windloc_{lab}_{wl}"] = lo <= mx <= hi and T_GRID[0] + W_B2 <= t_mx <= T_GRID[-1]
            gates[f"Wind_{lab}_{wl}"] = wrel <= AVG_BAR
            worst["q"] = max(worst["q"], rel)
            worst["w"] = max(worst["w"], wrel)
            worst["rc"] = min(worst["rc"], miss)
    _say(comm, f"[cf] (Q-d) worst rel {worst['q']:.3e} (bar {DOSE_BAR:g}, asserted); (Win-d) worst "
               f"rel {worst['w']:.3e} (bar {AVG_BAR:g}, asserted); R swap least miss "
               f"{worst['rc']:.4e} (asserted > {DOSE_BAR:g})")

    cf50 = _duty_series(DUTY)
    w50, l50 = 0.0, np.inf
    for j, lab in enumerate(LABS):
        d_cf = _dose(cf50[:, j])
        rel = abs(_dose(h["50"][:, j]) / d_cf - 1.0)
        miss = abs(_dose(h["50drop"][:, j]) / d_cf - 1.0)
        gates[f"D50cf_{lab}"] = rel <= DOSE_BAR
        gates[f"D50cf_drop_{lab}"] = miss > DOSE_BAR
        w50, l50 = max(w50, rel), min(l50, miss)
    _say(comm, f"[cf] (D50-cf) worst {w50:.3e} (bar {DOSE_BAR:g}, asserted); drop-run-0 least "
               f"miss {l50:.4e} (asserted > {DOSE_BAR:g})")

    j = 0  # P5
    ctrl = {}
    for wl, _w in WAITS:
        d_c = _dose(d25[:, j] + cap[wl][:, j])
        d_m = _dose(d25[:, j] + z1[wl][:, j])
        ctrl[wl] = (d_c, d_m)
        _say(comm, f"[ctrl] {LABS[j]} t_w {wl}: composite-drift dose {d_c:.9e} s, massless-wall "
                   f"(_drift_1d) dose {d_m:.9e} s, |massless/composite - 1| {abs(d_m / d_c - 1):.4e}, "
                   f"composite/massless - 1 {d_c / d_m - 1:+.4e}"
                   + (f" (ASSERTED > {DOSE_BAR:g})" if wl == "0" else " (printed)"))
    gates["ctrl_massless_misses_composite_P5_tw0"] = (
        abs(ctrl["0"][1] / ctrl["0"][0] - 1.0) > DOSE_BAR)

    _say(comm, "[TABLE-cf] point | D_ref s | D(0) s | D(T_C) s | D(T_C)/D(0) | D(0)/D_ref-1 | "
               "D(T_C)/D_ref-1 || six-face h_eff D(0)/D_ref-1 | D(T_C)/D_ref-1 || z-only "
               "_drift_1d D(0)/D_ref-1 | D(T_C)/D_ref-1   (composite first)")
    for j, lab in enumerate(LABS):
        ref = _dose(d25[:, j])
        d0 = _dose(d25[:, j] + cap["0"][:, j])
        dc = _dose(d25[:, j] + cap["T_C"][:, j])
        s0 = _dose(d25[:, j] + six[(j, "0")])
        sc = _dose(d25[:, j] + six[(j, "T_C")])
        m0 = _dose(d25[:, j] + z1["0"][:, j])
        mc = _dose(d25[:, j] + z1["T_C"][:, j])
        _say(comm, f"[TABLE-cf] {lab} | {ref:.6e} | {d0:.6e} | {dc:.6e} | {dc / d0:.6f} | "
                   f"{d0 / ref - 1:+.6e} | {dc / ref - 1:+.6e} || {s0 / ref - 1:+.6e} | "
                   f"{sc / ref - 1:+.6e} || {m0 / ref - 1:+.6e} | {mc / ref - 1:+.6e}")
    for wl, _w in WAITS:
        gap6 = max(float(np.max(np.abs(six[(j, wl)] - z1[wl][:, j]))) for j in (0, 1))
        gapc = max(float(np.max(np.abs(cap[wl][:, j] - z1[wl][:, j]))) for j in (0, 1))
        _say(comm, f"[SIDE] t_w {wl}: max |six-face - _drift_1d| on T_GRID at P5/ADJ_PT "
                   f"{gap6:.3e} K (geometry); max |composite - _drift_1d| {gapc:.3e} K (capacity)")
    return ctrl


def main() -> None:
    comm = MPI.COMM_WORLD
    started = time.perf_counter()

    if not np.issubdtype(np.dtype(default_scalar_type), np.complexfloating):
        raise RuntimeError(
            "this example expects the complex DolfinX build the runner's `mat:` group "
            "sources (source /usr/local/bin/dolfinx-complex-mode)"
        )

    _say(comm, "=" * 78)
    _say(comm, "EX-186 -- the D = 0.25 dose under the capacitive wall (MAT-21 P70 / P68)")
    _say(comm, "=" * 78)
    _say(comm, f"[gate] imported: DOSE_BAR {DOSE_BAR:g}, AVG_BAR {AVG_BAR:g}, INTERIOR_BAND "
               f"{INTERIOR_BAND:.0e}, E_BAND {E_BAND:.0e}, S0_BAND {S0_BAND:g}, PRE_BAR "
               f"{PRE_BAR:.0e}; D {DUTY25}, schedule {SCHEDULE_25}")
    assert K == K_GEL and H == H_GEL, (K, K_GEL, H, H_GEL)

    gates: dict = {}
    ctrl = closed_form(comm, gates)

    msh, ns, nw = _walled_box(RUNG, (False, False, True))
    chi, kf, inter, ct = _materials(msh)
    rcf = _rho_c_field(msh)
    ft = _robin_facets(msh, (A_BOX, B_BOX, C_BOX + D_WALL), (2,))
    n = msh.topology.index_map(3).size_global
    write_setup_figure(
        msh,
        ct,
        FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names={GEL: "gel, q pulsed (D = 0.25), theta_0 = g z", WALL: "acrylic wall (nominal)"},
        title="mat:39 -- gel box with acrylic walls on z = +/-C, Robin on the outer wall "
              "faces, D = 0.25 source, T(0) = T_inf + g z: y = 0 slice",
        slice_normal=(0.0, 1.0, 0.0),
        slice_origin=(0.0, 0.0, 0.0),
        comm=comm,
    )
    _say(comm, f"[mesh] walled box m = {RUNG} cells {tuple(ns)} N_w = {nw}: n = {n}; "
               f"t* {T_STAR:.3f} s, dt {DT:.4f} s, N {N_MAIN}; ranks {comm.size}")

    grad = lambda x: T_INF + G_GRAD * x[2]  # noqa: E731  (every cell, gel and wall)
    kw = dict(k=kf, q=Q * chi, h=H_GEL, t_inf=T_INF, t_init=T_INF, robin_tag=ROBIN, theta=0.5,
              startup_half_steps=2, dt=DT, n_steps=N_MAIN,
              snapshot_steps=tuple(range(N_MAIN + 1)), source_schedule=_chi25)
    starts = {"ref": None, "tw0": grad}
    res, E, hist, hist_ok, sched, s0ok, src = {}, {}, {}, {}, {}, {}, {}
    for key, expr in starts.items():
        tq = time.perf_counter()
        res[key] = solve_transient_conduction(msh, ft, t_init_expr=expr, rho_c=rcf, **kw)
        E[key] = float(np.real(res[key].max_balance_residual))
        src[key] = float(np.real(res[key].source_power))
        sched[key] = "".join(str(int(x)) for x in res[key].source_multipliers)
        try:
            _census(res[key], DT, 2, DUTY25 * T_STAR, f"EX-186 {key}")
            s0ok[key] = True
        except AssertionError:
            s0ok[key] = False
        hs = np.empty((N_MAIN + 1, HIST_PTS.shape[0]))
        ok_all = [True] * HIST_PTS.shape[0]
        for i in range(N_MAIN + 1):
            hs[i], ok = _read(res[key].snapshots[i], HIST_PTS)
            ok_all = [a and b for a, b in zip(ok_all, ok)]
        res[key].snapshots.clear()
        hist[key], hist_ok[key] = hs, ok_all
        _say(comm, f"[run {key}] {time.perf_counter() - tq:.1f} s; (E) {E[key]:.3e} "
                   f"(E_BAND {E_BAND:.0e}); source_power {src[key]:.9e} W; schedule "
                   f"{sched[key]}")
    imag = max(_max_imag(r.temperature) for r in res.values())

    tq = time.perf_counter()
    zf = _ZField()
    s_c = zf(msh, _series_z(RHO_C_W), T_STAR)
    s_d = fem.Function(s_c.function_space)
    s_d.x.array[:] = s_c.x.array
    for i, tm in enumerate(T_SW25):
        if i:
            s_d.x.array[:] += (-1.0) ** i * zf(msh, _series_z(RHO_C_W), T_STAR - tm).x.array
    d_0 = zf(msh, _drift_cap_z, T_STAR)
    d_0h = zf(msh, _drift_z, T_STAR)
    pre = comm.allreduce(zf.worst, op=MPI.MAX)
    _say(comm, f"[series] z-only comparands: {time.perf_counter() - tq:.1f} s; precondition "
               f"worst rel {pre:.3e} (bar {PRE_BAR:.0e})")

    th = {key: r.temperature - T_INF for key, r in res.items()}
    Ir = _rel(msh, th["ref"], s_d, inter)
    I0 = _rel(msh, th["tw0"], s_d + d_0, inter)
    I0h = _rel(msh, th["tw0"], s_d + d_0h, inter)
    Ic = _rel(msh, th["tw0"], s_c + d_0, inter)
    _say(comm, f"[I] t_w = 0 vs composite drift {I0:.6e} (bound {INTERIOR_BAND:.0e}, asserted; "
               f"P68 records 5.513765e-03); reference {Ir:.6e} (P68 records 1.849598706e-03); "
               f"I0h vs h_eff drift {I0h:.9e} (printed; P68 records 2.942439114e-02)")
    _say(comm, f"[control] continuous comparand on t_w = 0: {Ic:.4e} (asserted > "
               f"{INTERIOR_BAND:.0e})")

    for key in res:
        gates[f"E_{key}"] = E[key] <= E_BAND
        gates[f"S0_census_{key}"] = s0ok[key]
        gates[f"S0_schedule_{key}"] = sched[key] == SCHEDULE_25
    gates["schedule_string"] = _sched_string(_chi25) == SCHEDULE_25
    gates["E_source_power_identical"] = src["ref"] > 0.0 and src["ref"] == src["tw0"]
    for key in res:
        for j, lab in enumerate(HIST_LABELS):
            dose, quad, rel, _swap, miss = _q_check_g(T_INF + hist[key][:, j], T_GRID)
            _say(comm, f"[Q-wf] {key} {lab}: D {dose:.6e} s, quad {quad:.6e} s, rel {rel:.3e} "
                       f"(bar {DOSE_BAR:g}, asserted); R swap miss {miss:.4e} (printed)")
            gates[f"Qwf_{key}_{lab}"] = rel <= DOSE_BAR
    d_fem0 = _dose(hist["tw0"][:, 0])
    d_femr = _dose(hist["ref"][:, 0])
    miss = abs(d_fem0 / d_femr - 1.0)
    _say(comm, f"[control] P5: FEM t_w = 0 dose {d_fem0:.6e} s vs reference {d_femr:.6e} s, "
               f"miss {miss:.4e} (asserted > {DOSE_BAR:g})")
    d_comp = ctrl["0"][0]
    d_mless = ctrl["0"][1]
    _say(comm, f"[FEM] P5 t_w = 0: FEM D(0) {d_fem0:.6e} s beside composite {d_comp:.6e} s "
               f"(FEM/composite - 1 {d_fem0 / d_comp - 1:+.4e}), massless {d_mless:.6e} s "
               f"(FEM/massless - 1 {d_fem0 / d_mless - 1:+.4e}); P65 record {P65_FEM_D0_P5:.6e} s "
               f"(printed)")

    gates.update({
        "valid_points": all(all(o) for o in hist_ok.values()),
        "precondition_unique_z": pre <= PRE_BAR,
        "I_tw0_composite": I0 <= INTERIOR_BAND,
        "I_ref": Ir <= INTERIOR_BAND,
        "control_continuous_misses_I": Ic > INTERIOR_BAND,
        "control_tw0_dose_misses_ref_P5": miss > DOSE_BAR,
        "max_imag_zero": imag == 0.0,
    })
    _say(comm, f"max|Im theta| = {imag:.3e}")

    V = res["ref"].temperature.function_space
    outs = {}
    f = fem.Function(V, name="theta_tw0_t_star")
    f.x.array[:] = np.real(res["tw0"].temperature.x.array) - T_INF
    f.x.scatter_forward()
    outs[f.name] = f
    comp = fem.Function(V, name="theta_comparand_tw0")
    d_p1 = fem.Function(V)
    comp.interpolate(s_d)
    d_p1.interpolate(d_0)
    comp.x.array[:] = comp.x.array + d_p1.x.array
    comp.x.scatter_forward()
    outs[comp.name] = comp
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_combined", msh, ct, outs, comm=comm)
    adopt_host_ownership(OUTPUT_DIR, comm=comm)
    _say(comm, f"\n[output] {path} -- theta(t*) on the t_w = 0 run and its composite comparand")

    bad = {k: v for k, v in gates.items() if not v}
    _say(comm, f"[gates] {len(gates)} collected, {len(bad)} red: {bad}")
    _say(comm, f"\n[timing] total {time.perf_counter() - started:.1f} s at -n {comm.size}")
    assert not bad, bad
    _say(comm, "\nAll imported gates hold.")


if __name__ == "__main__":
    main()
