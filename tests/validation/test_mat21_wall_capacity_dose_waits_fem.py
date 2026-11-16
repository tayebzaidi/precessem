"""`MAT-21` step 41c — P73's wait sweep in FEM at t_w = 0.5/λ₁ and 2/λ₁ (§9 item 897, §10 P78,
chain P; brief `docs/planning/chunks/MAT-21.md` "§9 item 897").

Sets no band and moves none; every band is imported.  P74's construction (item 885, the
``FEM_EM_MAT21_P74_WAIT`` branch of `test_mat21_wall_pre_equilibration_duty25_dose_fem.py`'s
``runs``), copied — P68's knob refuses ``2`` in a landed W0 assertion, so it is not extended.
Runs: ref (θ₀ = 0) and tw0 (θ₀ = g·z on every cell), P68's calls verbatim; for each wait w in
``WAIT_KEYS`` (P73's ``W41[1:4]``) P74's two calls: (a) a source-off wait from g·z, ``N_w`` =
ceil(w / DT) steps of ``DT_w`` = w / ``N_w``, then (b) the tw0 call verbatim from (a)'s
temperature.  tw0_mw and its massless sums are dropped (not this step's).

Comparand on the unique dof z: P68's walled uniform sum ``s_d`` plus item 823's composite drift
``_drift_cap_z`` at t* + w (``d_w``); ``Ih`` with the h_eff ``_drift_z`` and the control ``Ic0``
with the t_w = 0 drift ``d_0``.

Asserted (collected, item 668): reproductions of P74's (I) ref / tw0 / tw1 and FEM D at P5 for
tw0 / tw1 (half the last printed digit, the I41 ruling); (I) tw05 / tw2 ≤ ``INTERIOR_BAND``;
(E) on every call; (S0) census + schedule on every main call, every wait's multipliers 0; the
hand-off per wait; (Q-wf) on every history; the unique-z precondition; (Mono) on the FEM doses
tw0, tw05, tw1, tw2 at P5 and ``ADJ_PT`` (one sign, |ΔD|/D > ``DOSE_BAR``).  Controls (asserted,
factor 1 at the gate's own band): (c-tw) ``Ic0`` > ``INTERIOR_BAND``; (c-D) tw05 / tw2 P5 dose
misses tw0's by > ``DOSE_BAR``.  Predicted factors are printed only (rule (e)).

Not claimed: the T_C wait on this schedule, the thickness reading, a dose band, a rate, a
recommended wait or thickness.
"""

from __future__ import annotations

import time

import numpy as np
import pytest
from dolfinx import fem
from mpi4py import MPI

from fem_em_solver.thermal import solve_transient_conduction

from tests.validation.test_imp7_duty_cycle_fem import S0_BAND
from tests.validation.test_mat10_sar_to_heat_sphere import _say
from tests.validation.test_mat12_transient_sar_heating_sphere import E_BAND
from tests.validation.test_mat13_transient_perfusion_sphere import INTERIOR_BAND
from tests.validation.test_mat15_coil_driven_heating_phantom import H_GEL, N_MAIN
from tests.validation.test_mat15_finite_cylinder_closed_form import Q, T_INF
from tests.validation.test_mat15_finite_cylinder_fem import _rel
from tests.validation.test_mat17_duty_cycle_fem import _census
from tests.validation.test_mat21_box_closed_form import A_BOX, B_BOX, C_BOX
from tests.validation.test_mat21_box_dose_closed_form import DOSE_BAR, T_GRID
from tests.validation.test_mat21_box_dt_duty125_dose_closed_form import _q_check_g
from tests.validation.test_mat21_box_duty25_dose_closed_form import (
    DUTY25,
    SCHEDULE_25,
    T_SW25,
    _chi25,
)
from tests.validation.test_mat21_initial_gradient_closed_form import G_GRAD
from tests.validation.test_mat21_pre_equilibration_duty_dose_closed_form import _dose
from tests.validation.test_mat21_wall_capacity_dose_sweep_closed_form import (
    P70_CEIL_P5,
    W41,
    _d_comp,
)
from tests.validation.test_mat21_wall_capacity_drift_closed_form import (
    LAM1_W,
    WAIT_C,
    _drift_cap,
    _drift_cap_z,
)
from tests.validation.test_mat21_wall_closed_form import D_WALL, RHO_C_W
from tests.validation.test_mat21_wall_fem import ROBIN, _materials, _robin_facets, _walled_box
from tests.validation.test_mat21_wall_pre_equilibration_duty25_dose_fem import (
    GRID_BAR,
    HANDOFF_BAR,
    HIST_LABELS,
    HIST_PTS,
    REPRO,
    TW1,
    _grid_offsets,
    _hist_series,
)
from tests.validation.test_mat21_wall_pre_equilibration_fem import (
    DT,
    PRE_BAR,
    _drift_z,
    _read,
    _series_z,
    _ZField,
)
from tests.validation.test_mat21_wall_transient_closed_form import T_STAR
from tests.validation.test_mat21_wall_transient_fem import RUNG, _rho_c_field

P = "MAT21-41c"
WAIT_KEYS = {"tw05": W41[1], "tw1": W41[2], "tw2": W41[3]}
WAIT_CS = {"tw05": WAIT_C[1], "tw1": WAIT_C[2], "tw2": WAIT_C[3]}
N_WAIT = {k: int(np.ceil(w / DT)) for k, w in WAIT_KEYS.items()}
DT_WAIT = {k: w / N_WAIT[k] for k, w in WAIT_KEYS.items()}
N_WAIT_PRED = {"tw05": 31, "tw1": 62, "tw2": 124}
NEW_KEYS = ("tw05", "tw2")
MONO_KEYS = ("tw0", "tw05", "tw1", "tw2")
#: Reproductions of P74's record (`20261108T225830Z_MAT-21-P74-W2.log`), bound = half the last
#: printed digit (the I41 ruling): (I) ref (:81, REPRO["Ir"]), tw0 (:81), tw1 (:75); FEM D at P5
#: tw0 (:99, :117) and tw1 (:117).
REPRO_I = {"ref": REPRO["Ir"], "tw0": (5.513765e-03, 5e-10), "tw1": (4.712673e-03, 5e-10)}
REPRO_D_P5 = {"tw0": (7.253626e-11, 5e-18), "tw1": (5.642660e-11, 5e-18)}
#: Predicted, printed only (item 897).
PRED_I = {"tw05": (4.7e-03, 5.6e-03), "tw2": (1.8e-03, 4.8e-03)}
PRED_CTW = {"tw05": 0.16, "tw2": 0.35}  # one-mode scaling of P74's I1c0 2.5314e-01
CTW_CEIL = 0.40  # ||d_0|| / ||s_d|| on the interior, same scaling (w -> inf)
PRED_CD = {"tw05": 0.150, "tw2": 0.279}  # 1 - P73's D(w)/D(0), P73-W1.log:66, :68
P73_DW_D0 = {"tw05": 0.849690, "tw2": 0.720512}
P73_D_DREF_P5 = {"tw05": 2.359837e-01, "tw2": 4.807859e-02}  # P73-W1.log:66, :68
P67_DRIFT_P5 = {"tw05": 1.098177e-01, "tw1": 6.658629e-02, "tw2": 2.449520e-02}  # P67-W1.log
TW1_IH_RATIO = 4.077  # P74-W2.log:75
FEM_COMP_P5 = {"tw0": 4.756e-03, "tw1": 4.900e-04}  # P74-W2.log:117 (tw0: P68-W1.log:91)


def test_w0_waits():
    """W0 (no mesh): W41[2] == TW1, N_w, the grid, w·LAM1_W, and P68's grid premise."""
    same = W41[2] == TW1
    grid = {k: abs(N_WAIT[k] * DT_WAIT[k] - w) for k, w in WAIT_KEYS.items()}
    lam = {k: abs(w * LAM1_W - WAIT_CS[k]) for k, w in WAIT_KEYS.items()}
    sw, gr = _grid_offsets()
    _say(f"\n{P} W0 W41[2] == TW1 bitwise {same} ({W41[2]!r} / {TW1!r}); N_w {N_WAIT} (must be "
         f"{N_WAIT_PRED}); DT_w " + ", ".join(f"{k} {v:.6f} s" for k, v in DT_WAIT.items()))
    _say(f"{P} W0 |N_w DT_w - w| {grid} (bar {1e-9 * DT:.3e}); |w LAM1_W - WAIT_C| {lam} (bar "
         f"1e-12); P68 grid premise T_SW25 {sw:.3e}, T_GRID {gr:.3e} steps (bar {GRID_BAR:.0e})")
    gates = {"tw1_bitwise": same, "n_w": N_WAIT == N_WAIT_PRED,
             "grid": all(v <= 1e-9 * DT for v in grid.values()),
             "lam": all(v <= 1e-12 for v in lam.values()),
             "p68_grid": sw <= GRID_BAR and gr <= GRID_BAR}
    _say(f"{P} W0 gates: {gates}")
    assert all(gates.values()), gates


@pytest.fixture(scope="module")
def runs():
    comm = MPI.COMM_WORLD
    t0 = time.perf_counter()
    msh, ns, nw = _walled_box(RUNG, (False, False, True))
    chi, kf, inter, _ct = _materials(msh)
    rcf = _rho_c_field(msh)
    ft = _robin_facets(msh, (A_BOX, B_BOX, C_BOX + D_WALL), (2,))
    _say(f"\n{P} walled box m = {RUNG} cells {tuple(ns)} N_w = {nw}: n = "
         f"{msh.topology.index_map(3).size_global}; D = {DUTY25}; CN N = {N_MAIN} + 2 BE "
         f"half-steps; waits {WAIT_KEYS}; batch {HIST_LABELS}; ranks {comm.size}")
    grad = lambda x: T_INF + G_GRAD * x[2]  # noqa: E731  (every cell, gel and wall)
    kw = dict(k=kf, q=Q * chi, h=H_GEL, t_inf=T_INF, t_init=T_INF, robin_tag=ROBIN,
              theta=0.5, startup_half_steps=2, dt=DT, n_steps=N_MAIN,
              snapshot_steps=tuple(range(N_MAIN + 1)), source_schedule=_chi25)
    starts = {"ref": None, "tw0": grad}
    starts.update({k: None for k in WAIT_KEYS})  # expr is call (a)'s temperature, set below
    out = {"E": {}, "sched": {}, "s0": {}, "hist": {}, "hist_ok": {}, "wall": {},
           "mult0": {}, "handoff": {}, "handoff_grid": {}, "keys": tuple(starts)}
    res = {}
    for key, expr in starts.items():
        tq = time.perf_counter()
        if key in WAIT_KEYS:  # P74 (a): the source-off wait, P68's kw and rcf
            w, n_w, dt_w = WAIT_KEYS[key], N_WAIT[key], DT_WAIT[key]
            kw_a = dict(kw, dt=dt_w, n_steps=n_w, snapshot_steps=(),
                        source_schedule=lambda t: 0.0)
            res_a = solve_transient_conduction(msh, ft, t_init_expr=grad, rho_c=rcf, **kw_a)
            out["E"][key + "_wait"] = float(res_a.max_balance_residual)
            mult_a = [float(x) for x in res_a.source_multipliers]
            out["mult0"][key] = len(mult_a) > 0 and all(x == 0.0 for x in mult_a)
            out["wall"][key + "_wait"] = time.perf_counter() - tq
            _say(f"{P} {key} (a) wait: {n_w} steps of {dt_w:.6f} s ({n_w * dt_w:.4f} s), "
                 f"{out['wall'][key + '_wait']:.1f} s; (E) {out['E'][key + '_wait']:.3e} (E_BAND "
                 f"{E_BAND:.0e}); multipliers all 0 {out['mult0'][key]} (n = {len(mult_a)})")
            expr = res_a.temperature
            tq = time.perf_counter()
        res[key] = solve_transient_conduction(msh, ft, t_init_expr=expr, rho_c=rcf, **kw)
        ts = time.perf_counter()
        if key in WAIT_KEYS:  # hand-off: (b)'s snapshot 0 == (a)'s final temperature
            loc = float(np.max(np.abs(res[key].snapshots[0].x.array
                                      - res_a.temperature.x.array)))
            out["handoff"][key] = comm.allreduce(loc, op=MPI.MAX)
            out["handoff_grid"][key] = abs(N_WAIT[key] * DT_WAIT[key] - WAIT_KEYS[key])
            _say(f"{P} {key} hand-off: max |snap0(b) - T(a)| {out['handoff'][key]:.3e} K (bar "
                 f"{HANDOFF_BAR:.0e}); |N_w DT_w - w| {out['handoff_grid'][key]:.3e} s (bar "
                 f"{1e-9 * DT:.3e})")
        out["E"][key] = float(res[key].max_balance_residual)
        out["sched"][key] = "".join(str(int(x)) for x in res[key].source_multipliers)
        try:  # item 747's census, run by run
            _census(res[key], DT, 2, DUTY25 * T_STAR, f"wall capacity waits d25 {key}")
            out["s0"][key] = True
        except AssertionError:
            out["s0"][key] = False
        hist = np.empty((N_MAIN + 1, HIST_PTS.shape[0]))
        ok_all = [True] * HIST_PTS.shape[0]
        for n in range(N_MAIN + 1):
            hist[n], ok = _read(res[key].snapshots[n], HIST_PTS)
            ok_all = [a and b for a, b in zip(ok_all, ok)]
        res[key].snapshots.clear()
        out["hist"][key], out["hist_ok"][key] = hist, ok_all
        out["wall"][key] = time.perf_counter() - tq
        _say(f"{P} run {key}: solve {ts - tq:.1f} s, reads {time.perf_counter() - ts:.1f} s; "
             f"(E) {out['E'][key]:.3e} (E_BAND {E_BAND:.0e}); source_power "
             f"{float(res[key].source_power):.9e} W")

    tq = time.perf_counter()
    zf = _ZField()
    s_c = zf(msh, _series_z(RHO_C_W), T_STAR)  # m = 0 (P68 :303-308)
    s_d = fem.Function(s_c.function_space)
    s_d.x.array[:] = s_c.x.array
    for i, tm in enumerate(T_SW25):
        if i:
            s_d.x.array[:] += (-1.0) ** i * zf(msh, _series_z(RHO_C_W), T_STAR - tm).x.array
    d_0 = zf(msh, _drift_cap_z, T_STAR)
    th = {key: r.temperature - T_INF for key, r in res.items()}
    out["I"] = {"ref": _rel(msh, th["ref"], s_d, inter),
                "tw0": _rel(msh, th["tw0"], s_d + d_0, inter)}
    out["Ih"], out["Ic0"], out["drift_p5"] = {}, {}, {}
    for key, w in WAIT_KEYS.items():
        d_w = zf(msh, _drift_cap_z, T_STAR + w)
        d_wh = zf(msh, _drift_z, T_STAR + w)
        out["I"][key] = _rel(msh, th[key], s_d + d_w, inter)
        out["Ih"][key] = _rel(msh, th[key], s_d + d_wh, inter)
        out["Ic0"][key] = _rel(msh, th[key], s_d + d_0, inter)
        out["drift_p5"][key] = _drift_cap(tuple(HIST_PTS[0]), T_STAR + w)
    out["pre"] = comm.allreduce(zf.worst, op=MPI.MAX)
    _say(f"{P} z-only comparands on the unique dof z ({zf.n_unique} values per rank max): "
         f"{time.perf_counter() - tq:.1f} s; precondition worst rel {out['pre']:.3e} (bar "
         f"{PRE_BAR:.0e})")
    for key in WAIT_KEYS:
        tag = ""
        if key in PRED_I:
            lo, hi = PRED_I[key]
            tag = (f"; predicted [{lo:.1e}, {hi:.1e}], printed, inside "
                   f"{lo <= out['I'][key] <= hi}; (c-tw) predicted ≈ {PRED_CTW[key]:.2f}, "
                   f"ceiling ≈ {CTW_CEIL:.2f} (printed)")
        _say(f"{P} (I) {key} vs composite drift at t* + w {out['I'][key]:.6e} (bound "
             f"{INTERIOR_BAND:.0e}){tag}; Ih vs h_eff drift {out['Ih'][key]:.6e}, Ih / I "
             f"{out['Ih'][key] / out['I'][key]:.3f}x (tw1 {TW1_IH_RATIO}x, printed); control "
             f"(c-tw) vs the t_w = 0 comparand {out['Ic0'][key]:.4e} (asserted > "
             f"{INTERIOR_BAND:.0e})")
        _say(f"{P} composite drift at P5 at t* + w ({key}): {out['drift_p5'][key]:.6e} K beside "
             f"P67's {P67_DRIFT_P5[key]:.6e} K (printed)")
    _say(f"{P} (I) ref {out['I']['ref']:.6e}, tw0 {out['I']['tw0']:.6e} (reproductions)")
    _say(f"{P} elapsed {time.perf_counter() - t0:.1f} s")
    return out


def test_all_gates(runs):
    """Collect every gate before asserting (item 668), so one red hides no other."""
    r = runs
    gates = {"valid_points": all(all(o) for o in r["hist_ok"].values()),
             "precondition_unique_z": r["pre"] <= PRE_BAR}
    for key, (val, bnd) in REPRO_I.items():
        _say(f"{P} reproduction (I) {key}: {r['I'][key]:.9e} vs {val:.6e}, |diff| "
             f"{abs(r['I'][key] - val):.3e} (bound {bnd:.0e}, ASSERTED)")
        gates[f"repro_I_{key}"] = abs(r["I"][key] - val) <= bnd
    j5 = HIST_LABELS.index("P5")
    dose = {key: [_dose(r["hist"][key][:, j]) for j in range(len(HIST_LABELS))]
            for key in r["keys"]}
    for key, (val, bnd) in REPRO_D_P5.items():
        _say(f"{P} reproduction FEM D P5 {key}: {dose[key][j5]:.9e} s vs {val:.6e} s, |diff| "
             f"{abs(dose[key][j5] - val):.3e} (bound {bnd:.0e}, ASSERTED)")
        gates[f"repro_D_P5_{key}"] = abs(dose[key][j5] - val) <= bnd
    for key in NEW_KEYS:
        gates[f"I_{key}"] = r["I"][key] <= INTERIOR_BAND
        gates[f"control_c_tw_{key}"] = r["Ic0"][key] > INTERIOR_BAND
    for key in r["keys"]:
        gates[f"E_{key}"] = r["E"][key] <= E_BAND
        gates[f"S0_census_{key}"] = r["s0"][key]
        gates[f"S0_schedule_{key}"] = r["sched"][key] == SCHEDULE_25
        _say(f"{P} (S0) {key}: schedule {r['sched'][key]} (must equal {SCHEDULE_25})")
    for key in WAIT_KEYS:
        gates[f"E_{key}_wait"] = r["E"][key + "_wait"] <= E_BAND
        gates[f"S0_{key}_wait_mult0"] = r["mult0"][key]
        gates[f"handoff_{key}"] = (r["handoff"][key] <= HANDOFF_BAR
                                   and r["handoff_grid"][key] <= 1e-9 * DT)
    for key in r["keys"]:
        for j, lab in enumerate(HIST_LABELS):
            dq, quad, rel, _swap, _miss = _q_check_g(T_INF + r["hist"][key][:, j], T_GRID)
            _say(f"{P} (Q-wf) {key} {lab}: D {dq:.6e} s, quad {quad:.6e} s, rel {rel:.3e} "
                 f"(bar {DOSE_BAR:g}, ASSERTED)")
            gates[f"Q_wf_{key}_{lab}"] = rel <= DOSE_BAR
    for key in NEW_KEYS:  # (c-D)
        miss = abs(dose[key][j5] / dose["tw0"][j5] - 1.0)
        _say(f"{P} control (c-D) P5: {key} FEM D {dose[key][j5]:.6e} s vs tw0 "
             f"{dose['tw0'][j5]:.6e} s, miss {miss:.4e} (asserted > {DOSE_BAR:g}; predicted "
             f"{PRED_CD[key]:.3f}, printed; ceiling D(0)/D_ref - 1 {P70_CEIL_P5:+.6e})")
        gates[f"control_c_D_{key}"] = miss > DOSE_BAR
    for j, lab in enumerate(HIST_LABELS):  # (Mono) on the FEM
        ds = [dose[k][j] for k in MONO_KEYS]
        diffs = [b - a for a, b in zip(ds[:-1], ds[1:])]
        signs = {int(np.sign(d)) for d in diffs}
        rels = [abs(d) / a for d, a in zip(diffs, ds[:-1])]
        sign = ("decreasing" if signs == {-1} else "increasing" if signs == {1}
                else f"mixed {signs}")
        _say(f"{P} (Mono) {lab}: D {dict(zip(MONO_KEYS, (f'{d:.6e}' for d in ds)))} s; "
             f"sign {sign} (P73 predicts decreasing, printed); smallest |dD|/D {min(rels):.4e} "
             f"(asserted > {DOSE_BAR:g})")
        gates[f"Mono_{lab}"] = len(signs) == 1 and 0 not in signs and min(rels) > DOSE_BAR
    _say(f"{P} gates: {gates}")
    assert all(gates.values()), gates


def test_printed_dose_table(runs):
    """The record (no band); asserts only that it was formed."""
    uni = _hist_series()[0]
    _say(f"{P} TABLE point | run | D FEM s | composite series D s | P73 _d_comp s | "
         "|history/_d_comp - 1| | h_eff series D s | FEM/composite-1 | FEM/h_eff-1 | "
         "FEM D/D_ref-1 | composite D/D_ref-1")
    for j, lab in enumerate(HIST_LABELS):
        d_rf = _dose(runs["hist"]["ref"][:, j])
        d_rs = _dose(uni[:, j])
        for key, w in WAIT_KEYS.items():
            _u, wc, wh, _d = _hist_series(w)
            d_f = _dose(runs["hist"][key][:, j])
            d_c = _dose(wc[:, j])
            d_p = _d_comp(j, w)
            d_hh = _dose(wh[:, j])
            _say(f"{P} TABLE {lab} | {key} | {d_f:.6e} | {d_c:.6e} | {d_p:.6e} | "
                 f"{abs(d_c / d_p - 1.0):.3e} | {d_hh:.6e} | {d_f / d_c - 1.0:+.3e} | "
                 f"{d_f / d_hh - 1.0:+.3e} | {d_f / d_rf - 1.0:+.6e} | {d_c / d_rs - 1.0:+.6e}")
            assert np.isfinite(d_f / d_c)
    _say(f"{P} beside (printed): FEM/composite-1 at P5 tw0 {FEM_COMP_P5['tw0']:+.3e}, tw1 "
         f"{FEM_COMP_P5['tw1']:+.3e}; P73 D/D_ref-1 at P5 tw05 {P73_D_DREF_P5['tw05']:+.6e}, "
         f"tw2 {P73_D_DREF_P5['tw2']:+.6e}; P73 D(w)/D(0) {P73_DW_D0}")
    _say(f"{P} wall time per call (solve + reads): "
         + ", ".join(f"{k} {v:.1f} s" for k, v in runs["wall"].items()))
