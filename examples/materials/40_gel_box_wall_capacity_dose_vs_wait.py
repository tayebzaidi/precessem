"""Example (`EX-189`): the D = 0.25 thermal dose against the wait, with and without the wall's
heat capacity in the drift -- the closed form (P73) and one FEM point at t_w = 1/lambda_1 (P74).

`mat:39` shows the D = 0.25 dose under the capacitive wall at t_w in {0, T_C}, FEM at t_w = 0
only.  This page's output is ``cem43`` at D = 0.25 as a function of the wait
t_w in {0, 0.5, 1, 2}/lambda_1 (``W41``), composite drift against the massless-wall
(``_drift_1d``) drift, on item 682's 1D-reduced walled gel box (m = 4, 119 808 cells), plus an
FEM dose at t_w = 1/lambda_1 (``TW1``), which no earlier example shows.

Two halves.

- **Closed form (numpy)**, P73's helpers imported by name
  (``test_mat21_wall_capacity_dose_sweep_closed_form.py``: ``W41``, ``_cap_hist_w``,
  ``_d_comp``, ``_d_mass``, ``_side_row``, and through P64's module ``Z``, ``_dose``,
  ``_hists``): the composite-drift dose, the massless-wall dose and their ratio at P5 /
  ``ADJ_PT`` / P0 at the four waits.
- **FEM**, item 885's knob run (P68's module with ``FEM_EM_MAT21_P74_WAIT=1``) written out
  here (fixtures cannot be imported), its plain helpers imported by name: the reference run
  (theta_0 = 0) and ``tw1`` -- (a) a source-off wait of ``TW1`` from theta_0 = g z, then (b) the
  D = 0.25 call from (a)'s temperature -- every snapshot of the reference run and of (b) read at
  P5 / ``ADJ_PT``.

No ``test_*`` name is imported.  No gate edit, no ``src/`` edit.

Asserted (all collected first, then asserted): P73's (Lim) at t_w = 0 against P70's
construction at ``B1_BAND`` and the reproduction of P70's printed digits, (Mono) at P5 /
``ADJ_PT`` (one sign, each |dD|/D > ``DOSE_BAR``) and constancy at P0 within ``DOSE_BAR``,
(Q-d) with the R-swap control and (Win-d) per point and wait; the control that the massless-wall
dose misses the composite dose at P5 at every wait by more than ``DOSE_BAR``; P74's (I) on tw1
against the composite drift at t* + ``TW1`` <= ``INTERIOR_BAND``, (E) <= ``E_BAND`` on the
reference run and on both tw1 calls, the hand-off precondition (<= ``HANDOFF_BAR``), (S0) on
(b) and (a)'s multipliers all 0, (Q-wf) at ``DOSE_BAR``, the controls (c-tw) and (c-D) (tw1's
FEM dose at P5 misses P68's t_w = 0 record by more than ``DOSE_BAR``); max|Im theta| = 0.
Printed only: the dose-by-wait table, the thickness side table (the heading's reading), the
FEM dose at t_w = 1/lambda_1 beside P73's ``_d_comp``.

The wall's k and rho_c are nominal PMMA, citation owed (W3).  One rung, one time step.
Not claimed: the six-wall box, a dose band, a recommended wait or thickness, a cited PMMA
rho_c, a tissue dose.

Run it::

    ./run_examples.sh -e mat:40 -n 2 -t 480

Output: ``examples/materials/paraview_output/
materials_40_gel_box_wall_capacity_dose_vs_wait_combined.xdmf`` (theta at t* on the tw1 run
and its composite comparand).
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

from tests.validation.test_imp4_pec_cylinder_closed_form import B1_BAND  # noqa: E402
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
from tests.validation.test_mat21_box_closed_form import A_BOX, B_BOX, C_BOX  # noqa: E402
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
)
from tests.validation.test_mat21_initial_gradient_closed_form import G_GRAD  # noqa: E402
from tests.validation.test_mat21_pre_equilibration_closed_form import ADJ_PT, P5  # noqa: E402
from tests.validation.test_mat21_pre_equilibration_duty_dose_closed_form import LABS  # noqa: E402
from tests.validation.test_mat21_wall_capacity_dose_sweep_closed_form import (  # noqa: E402
    D_FACTORS,
    P70_D0,
    P70_RATIO_P5,
    W41,
    W41_LABELS,
    _cap_hist_w,
    _d_comp,
    _d_mass,
    _side_row,
)
from tests.validation.test_mat21_wall_capacity_drift_closed_form import (  # noqa: E402
    LAM1_W,
    _drift_cap_z,
)
from tests.validation.test_mat21_wall_capacity_duty25_dose_closed_form import (  # noqa: E402
    P65_FEM_D0_P5,
    _cap_hists,
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
    Z,
    _dose,
    _hists,
)
from tests.validation.test_mat21_wall_pre_equilibration_duty25_dose_fem import (  # noqa: E402
    DT_W,
    HANDOFF_BAR,
    HIST_PTS,
    N_W,
    PREDICTED_I_TW1,
    REPRO,
    TW0_FEM_SERIES_P5,
    TW1,
)
from tests.validation.test_mat21_wall_pre_equilibration_fem import (  # noqa: E402
    DT,
    PRE_BAR,
    _read,
    _series_z,
    _ZField,
)
from tests.validation.test_mat21_wall_transient_closed_form import T_STAR  # noqa: E402
from tests.validation.test_mat21_wall_transient_fem import RUNG, _rho_c_field  # noqa: E402

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"  # committed, EX-57
BASENAME = "materials_40_gel_box_wall_capacity_dose_vs_wait"
HIST_LABELS = ("P5", "ADJ_PT")


def _say(comm, msg):
    if comm.rank == 0:
        print(msg, flush=True)


def closed_form(comm, gates):
    """P73's anchors through its helpers (numpy, rank-identical)."""
    d25 = _hists()["25"]
    cap70 = _cap_hists()["0"]

    # (Lim) + reproduction of P70's printed digits
    for j, lab in enumerate(LABS):
        d70 = _dose(d25[:, j] + cap70[:, j])
        d0 = _d_comp(j, 0.0)
        rel = abs(d0 - d70) / d70
        val, tol = P70_D0[j]
        _say(comm, f"[cf] (Lim) {lab}: D(t_w=0) {d0:.9e} s vs P70 construction {d70:.9e} s, rel "
                   f"{rel:.3e} (band {B1_BAND:g}, asserted); vs P70 printed {val:.9e} s |diff| "
                   f"{abs(d0 - val):.3e} s (bound {tol:g}, asserted)")
        gates[f"Lim_{lab}"] = rel <= B1_BAND
        gates[f"Lim_repro_{lab}"] = abs(d0 - val) <= tol
    r_p5 = _d_comp(0, 0.0) / _d_mass(0, 0.0) - 1.0
    gates["Lim_repro_ratio_P5"] = abs(r_p5 - P70_RATIO_P5[0]) <= P70_RATIO_P5[1]
    _say(comm, f"[cf] (Lim) P5 composite/massless - 1 at t_w = 0 {r_p5:+.7e} (P70 {P70_RATIO_P5[0]:+.4e},"
               f" bound {P70_RATIO_P5[1]:g}, asserted)")

    # (Mono)
    doses = {}
    for j, lab in enumerate(LABS):
        ds = np.array([_d_comp(j, w) for w in W41])
        doses[lab] = ds
        if j < 2:
            diff = np.diff(ds)
            sgn = np.sign(diff)
            sep = np.abs(diff) / ds[:-1]
            gates[f"Mono_{lab}"] = bool(np.all(sgn == sgn[0]) and sgn[0] != 0
                                        and np.all(sep > DOSE_BAR))
            _say(comm, f"[cf] (Mono) {lab}: dD signs {sgn.tolist()}, |dD|/D "
                       f"{', '.join(f'{s:.4e}' for s in sep)} (asserted one sign, > {DOSE_BAR:g}); "
                       f"sign {'decreasing' if sgn[0] < 0 else 'increasing'}")
        else:
            dev = float(np.max(np.abs(ds / ds[0] - 1.0)))
            gates[f"Const_{lab}"] = dev <= DOSE_BAR
            _say(comm, f"[cf] (Mono) {lab}: max |D/D(0) - 1| {dev:.3e} (asserted <= {DOSE_BAR:g})")

    # (Q-d), (Win-d)
    worst = {"q": 0.0, "w": 0.0, "rc": np.inf}
    for j, lab in enumerate(LABS):
        for w, wl in zip(W41, W41_LABELS):
            th = d25[:, j] + _cap_hist_w(w)[:, j]
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

    # control: massless wall misses the composite at P5, every wait
    ctrl = {}
    for w, wl in zip(W41, W41_LABELS):
        d_c, d_m = _d_comp(0, w), _d_mass(0, w)
        ctrl[wl] = (d_c, d_m)
        miss = abs(d_m / d_c - 1.0)
        gates[f"ctrl_massless_P5_{wl}"] = miss > DOSE_BAR
        _say(comm, f"[ctrl] P5 t_w {wl}: composite {d_c:.9e} s, massless {d_m:.9e} s, |miss| "
                   f"{miss:.4e} (ASSERTED > {DOSE_BAR:g}); composite/massless - 1 "
                   f"{d_c / d_m - 1.0:+.4e}")

    # printed: the dose-by-wait table
    _say(comm, "[TABLE-w] point | t_w | t_w s | D(t_w) s | D(t_w)/D_ref - 1 | massless D s | "
               "composite/massless - 1")
    for j, lab in enumerate(LABS):
        ref = _dose(d25[:, j])
        for w, wl in zip(W41, W41_LABELS):
            d = _d_comp(j, w)
            dm = _d_mass(j, w)
            _say(comm, f"[TABLE-w] {lab} | {wl} | {w:.2f} | {d:.6e} | {d / ref - 1.0:+.6e} | "
                       f"{dm:.6e} | {d / dm - 1.0:+.4e}")
        _say(comm, f"[TABLE-w] {lab} D_ref {ref:.6e} s")

    # printed: the thickness side table, the heading's reading
    _say(comm, "[SIDE] thickness side table at t_w = 0, the heading's reading, printed only: "
               "d/D_WALL | d m | lambda1_odd /s | D P5 | D ADJ_PT | D P0 | D/D(D_WALL) - 1 at P5")
    base, _l = _side_row(D_WALL)
    for f in D_FACTORS:
        d = f * D_WALL
        try:
            row, l1 = _side_row(d)
        except Exception as exc:  # the table never gates
            _say(comm, f"[SIDE] {f:g} | {d:.4f} | REFUSED: {exc!r} (row skipped)")
            continue
        _say(comm, f"[SIDE] {f:g} | {d:.4f} | {l1:.6e} | {row[0]:.6e} | {row[1]:.6e} | "
                   f"{row[2]:.6e} | {row[0] / base[0] - 1.0:+.6e}")
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
    _say(comm, "EX-189 -- the D = 0.25 dose against the wait, wall capacity (MAT-21 P73 / P74)")
    _say(comm, "=" * 78)
    _say(comm, f"[gate] imported: B1_BAND {B1_BAND:g}, DOSE_BAR {DOSE_BAR:g}, AVG_BAR "
               f"{AVG_BAR:g}, INTERIOR_BAND {INTERIOR_BAND:.0e}, E_BAND {E_BAND:.0e}, S0_BAND "
               f"{S0_BAND:g}, PRE_BAR {PRE_BAR:.0e}, HANDOFF_BAR {HANDOFF_BAR:.0e}; D {DUTY25}, "
               f"schedule {SCHEDULE_25}; TW1 {TW1:.6f} s (|TW1 LAM1_W - 1| "
               f"{abs(TW1 * LAM1_W - 1.0):.3e}), N_W {N_W}, DT_W {DT_W:.6f} s")
    assert K == K_GEL and H == H_GEL, (K, K_GEL, H, H_GEL)

    gates: dict = {}
    ctrl = closed_form(comm, gates)
    gates["tw1_lam"] = abs(TW1 * LAM1_W - 1.0) <= 1e-12
    gates["n_w"] = N_W == 62

    msh, ns, nw = _walled_box(RUNG, (False, False, True))
    chi, kf, inter, ct = _materials(msh)
    rcf = _rho_c_field(msh)
    ft = _robin_facets(msh, (A_BOX, B_BOX, C_BOX + D_WALL), (2,))
    n = msh.topology.index_map(3).size_global
    write_setup_figure(
        msh,
        ct,
        FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names={GEL: "gel, q pulsed (D = 0.25) after the wait", WALL: "acrylic wall (nominal)"},
        title="mat:40 -- walled gel box, wait 1/lambda_1 from T(0) = T_inf + g z, "
              "then D = 0.25: y = 0 slice",
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

    res, E, hist, hist_ok, sched, s0ok, src = {}, {}, {}, {}, {}, {}, {}

    def _record(key, tq):
        E[key] = float(np.real(res[key].max_balance_residual))
        src[key] = float(np.real(res[key].source_power))
        sched[key] = "".join(str(int(x)) for x in res[key].source_multipliers)
        try:
            _census(res[key], DT, 2, DUTY25 * T_STAR, f"EX-189 {key}")
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

    tq = time.perf_counter()
    res["ref"] = solve_transient_conduction(msh, ft, t_init_expr=None, rho_c=rcf, **kw)
    _record("ref", tq)

    # tw1 (a): the source-off wait from theta_0 = g z; (b): the D = 0.25 call from its end state
    tq = time.perf_counter()
    kw_a = dict(kw, dt=DT_W, n_steps=N_W, snapshot_steps=(), source_schedule=lambda t: 0.0)
    res_a = solve_transient_conduction(msh, ft, t_init_expr=grad, rho_c=rcf, **kw_a)
    e_wait = float(np.real(res_a.max_balance_residual))
    mult_a = [float(np.real(x)) for x in res_a.source_multipliers]
    wait_mult0 = len(mult_a) > 0 and all(x == 0.0 for x in mult_a)
    _say(comm, f"[tw1 (a)] wait: {N_W} steps of {DT_W:.6f} s ({N_W * DT_W:.4f} s), "
               f"{time.perf_counter() - tq:.1f} s; (E) {e_wait:.3e} (E_BAND {E_BAND:.0e}); "
               f"multipliers all 0 {wait_mult0} (n = {len(mult_a)})")
    tq = time.perf_counter()
    res["tw1"] = solve_transient_conduction(msh, ft, t_init_expr=res_a.temperature, rho_c=rcf,
                                            **kw)
    loc = float(np.max(np.abs(res["tw1"].snapshots[0].x.array - res_a.temperature.x.array)))
    handoff = comm.allreduce(loc, op=MPI.MAX)
    handoff_grid = abs(N_W * DT_W - TW1)
    _say(comm, f"[hand-off] max |snap0(b) - T(a)| {handoff:.3e} K (bar {HANDOFF_BAR:.0e}, "
               f"asserted); |N_W DT_W - TW1| {handoff_grid:.3e} s (bar {1e-9 * DT:.3e})")
    _record("tw1", tq)
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
    d_1 = zf(msh, _drift_cap_z, T_STAR + TW1)
    pre = comm.allreduce(zf.worst, op=MPI.MAX)
    _say(comm, f"[series] z-only comparands: {time.perf_counter() - tq:.1f} s; precondition "
               f"worst rel {pre:.3e} (bar {PRE_BAR:.0e})")

    th = {key: r.temperature - T_INF for key, r in res.items()}
    Ir = _rel(msh, th["ref"], s_d, inter)
    I1 = _rel(msh, th["tw1"], s_d + d_1, inter)
    I1c0 = _rel(msh, th["tw1"], s_d + d_0, inter)
    lo1, hi1 = PREDICTED_I_TW1
    _say(comm, f"[I] tw1 vs composite drift at t* + 1/LAM1_W {I1:.6e} (bound {INTERIOR_BAND:.0e}, "
               f"asserted; item 885 predicted [{lo1:.1e}, {hi1:.1e}], printed, inside "
               f"{lo1 <= I1 <= hi1}); reference {Ir:.9e} (P68 records {REPRO['Ir'][0]:.6e}, bound "
               f"{REPRO['Ir'][1]:g}, asserted)")
    _say(comm, f"[control] (c-tw) tw1 against the t_w = 0 comparand: {I1c0:.4e} (asserted > "
               f"{INTERIOR_BAND:.0e})")

    for key in ("ref", "tw1"):
        gates[f"E_{key}"] = E[key] <= E_BAND
        gates[f"S0_census_{key}"] = s0ok[key]
        gates[f"S0_schedule_{key}"] = sched[key] == SCHEDULE_25
    gates["E_tw1_wait"] = e_wait <= E_BAND
    gates["S0_tw1_wait_mult0"] = wait_mult0
    gates["handoff"] = handoff <= HANDOFF_BAR and handoff_grid <= 1e-9 * DT
    gates["E_source_power_identical"] = src["ref"] > 0.0 and src["ref"] == src["tw1"]
    for key in ("ref", "tw1"):
        for j, lab in enumerate(HIST_LABELS):
            dose, quad, rel, _swap, miss = _q_check_g(T_INF + hist[key][:, j], T_GRID)
            _say(comm, f"[Q-wf] {key} {lab}: D {dose:.6e} s, quad {quad:.6e} s, rel {rel:.3e} "
                       f"(bar {DOSE_BAR:g}, asserted); R swap miss {miss:.4e} (printed)")
            gates[f"Qwf_{key}_{lab}"] = rel <= DOSE_BAR

    d_f1 = _dose(hist["tw1"][:, 0])
    d_ser = _d_comp(0, TW1)
    d_m1 = _d_mass(0, TW1)
    d_ref = _dose(hist["ref"][:, 0])
    miss_cd = abs(d_f1 / P65_FEM_D0_P5 - 1.0)
    _say(comm, f"[control] (c-D) P5: tw1 FEM D {d_f1:.6e} s vs P68 tw0 record {P65_FEM_D0_P5:.6e} s, "
               f"miss {miss_cd:.4e} (asserted > {DOSE_BAR:g})")
    _say(comm, f"[FEM] P5 t_w = 1/lambda_1: FEM D {d_f1:.6e} s beside P73 _d_comp(0, TW1) "
               f"{d_ser:.6e} s (FEM/series - 1 {d_f1 / d_ser - 1.0:+.4e}; P68 tw0 "
               f"{TW0_FEM_SERIES_P5:+.3e}, printed); massless {d_m1:.6e} s (FEM/massless - 1 "
               f"{d_f1 / d_m1 - 1.0:+.4e}); FEM reference-run D {d_ref:.6e} s (tw1/ref - 1 "
               f"{d_f1 / d_ref - 1.0:+.6e}); composite/massless - 1 "
               f"{ctrl[W41_LABELS[2]][0] / ctrl[W41_LABELS[2]][1] - 1.0:+.4e} (printed)")

    gates.update({
        "valid_points": all(all(o) for o in hist_ok.values()),
        "precondition_unique_z": pre <= PRE_BAR,
        "I_tw1_composite": I1 <= INTERIOR_BAND,
        "I_ref_repro": abs(Ir - REPRO["Ir"][0]) <= REPRO["Ir"][1],
        "control_c_tw": I1c0 > INTERIOR_BAND,
        "control_c_D": miss_cd > DOSE_BAR,
        "max_imag_zero": imag == 0.0,
    })
    _say(comm, f"max|Im theta| = {imag:.3e}")

    V = res["ref"].temperature.function_space
    outs = {}
    f = fem.Function(V, name="theta_tw1_t_star")
    f.x.array[:] = np.real(res["tw1"].temperature.x.array) - T_INF
    f.x.scatter_forward()
    outs[f.name] = f
    comp = fem.Function(V, name="theta_comparand_tw1")
    d_p1 = fem.Function(V)
    comp.interpolate(s_d)
    d_p1.interpolate(d_1)
    comp.x.array[:] = comp.x.array + d_p1.x.array
    comp.x.scatter_forward()
    outs[comp.name] = comp
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_combined", msh, ct, outs, comm=comm)
    adopt_host_ownership(OUTPUT_DIR, comm=comm)
    _say(comm, f"\n[output] {path} -- theta(t*) on the tw1 run and its composite comparand")

    bad = {k: v for k, v in gates.items() if not v}
    _say(comm, f"[gates] {len(gates)} collected, {len(bad)} red: {bad}")
    _say(comm, f"\n[timing] total {time.perf_counter() - started:.1f} s at -n {comm.size}")
    assert not bad, bad
    _say(comm, "\nAll imported gates hold.")


if __name__ == "__main__":
    main()
