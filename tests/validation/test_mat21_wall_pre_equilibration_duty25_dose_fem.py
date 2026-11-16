"""`MAT-21` step 37 — the walled box under D = 0.25 from the drifted start, FEM (§9 item 809,
§10 chain P, P65; brief `docs/planning/chunks/MAT-21.md` "§9 item 809").

Sets no band and moves none; every band is imported.  Item 747's route (the wall-free twin,
`test_mat21_pre_equilibration_duty_dose_fem.py`) on item 682's 1D-reduced walled box, built as
item 805 builds it (`test_mat21_wall_pre_equilibration_fem.py`, imported by name: the build, the
unique-z comparand route ``_ZField``, the point read).  Two runs, reference (θ₀ = 0) and t_w = 0
(θ₀ = g·z on every cell, gel and wall), ``source_schedule=_chi25`` (P64's schedule, ``T_SW25``,
P = 2 ``HALF_T``), to t\\*, ``snapshot_steps=tuple(range(N_MAIN + 1))``; each run's snapshots are
read on one fixed batch (``HIST_PTS`` = P5, ``ADJ_PT``) and cleared.

Comparand: the walled uniform sum Σ(−1)^m ``_series(z, t* − t_m, RHO_C_W, lam)`` over
``T_SW25`` on the unique dof z, plus ``_drift_1d(p, t*)`` on the t_w = 0 run.  Histories:
item 789's ``_uniform_w`` with ``_drift_1d`` in place of the six-face ``_drift`` (the
difference printed).

Asserted (collected, item 668): (I) interior rel L2 ≤ ``INTERIOR_BAND`` at t\\* on both runs;
(E) ≤ ``E_BAND`` on both; (S0) item 747's ``_census`` at ``S0_BAND`` on each run's recorded
multipliers, and the schedule string; (Q-wf) ``cem43`` against ``_q_check_g``'s quadrature at
``DOSE_BAR`` on all four histories.  Controls (asserted, item 747's): the t_w = 0 FEM dose at P5
misses the reference run's by > ``DOSE_BAR``; the continuous comparand misses (I) by
> ``INTERIOR_BAND``.  Printed (records, no band): D at P5 and ``ADJ_PT``, FEM beside item 789's
walled series, the signed gap.  W0, no solve: ``DT_BOX`` == item 682's Δt, every ``T_SW25``
entry on the step grid, the schedule string and census arithmetic.

Attempt 2 (§9 item 812, the "Review ruling 2026-11-04 03:00"): a third run ``tw0_mw``, the
t_w = 0 run with ρ_wc_w = ``MASSLESS`` · ``RHO_C_W`` (item 805's R3 form), against
Σ(−1)^m ``_series_z(MASSLESS·RHO_C_W)``(z, t* − t_m) over ``T_SW25`` + ``_drift_1d`` — exact up to
discretisation on that twin.  Gated (I-mw25), (E), (S0), (Q-wf) on it; printed (no band) the
FEM-vs-FEM wall-capacity term (tw0 − tw0_mw) at P5, ``ADJ_PT``, P0 at t\\* and (I-mw25) at
snapshot ``N_MAIN − 1`` beside ``N_MAIN`` (the CN-ringing signature).  At attempt 2 the t_w = 0
run's (I) against the h_eff drift was pre-registered red in [2.6e-02, 3.3e-02] (attempt 1
measured 2.942439e-02, `20261104T061055Z_MAT-21-P65-W1.log`); attempt 3 retires that prediction.

Attempt 3 (§9 item 835, §10 P68, the "Review ruling 2026-11-05 03:00"): the t_w = 0 run's
comparand drift is item 823's composite (capacitive-wall) odd-mode drift ``_drift_cap_z``
(`test_mat21_wall_capacity_drift_closed_form.py`) in place of the h_eff ``_drift_z``; (I) on tw0
is asserted against it (no prediction line).  (I-mw25) keeps ``_drift_z`` — exact on the massless
twin ((Lim-c)) — and ``d_prev`` with it.  ``I0h`` (tw0 against the h_eff comparand) is printed and
reproduction-controlled at attempts 1 / 2's 2.942439e-02 (bound 5e-9), with (I) ref 1.849599e-03
and (I-mw25) 5.708485e-03 (bounds 5e-10; half the print's last digit, the I41 ruling).  The tw0
history comparand uses ``_drift_cap``; the ``_drift_1d`` history is printed beside it.

Step 41b (§9 item 885, §10 P74; brief `docs/planning/chunks/MAT-21.md` "§9 item 885"): the knob
``FEM_EM_MAT21_P74_WAIT`` (unset / empty: the landed three runs; ``1``: adds run ``tw1``; anything
else raises ``ValueError``, read inside the fixture, never at import).  ``tw1`` is two calls: (a) a
source-off wait from θ₀ = g·z for ``TW1`` = 1/``LAM1_W`` (``N_W`` steps of ``DT_W``), then (b) the
tw0 call verbatim from (a)'s temperature.  Gated with the knob on: (I) tw1 against the composite
drift at t* + ``TW1``, (E) on both calls, (S0) on (b) and (a)'s multipliers all 0, (Q-wf) on tw1's
histories, the hand-off precondition; controls (c-tw) (tw1 vs the t_w = 0 comparand misses (I))
and (c-D) (tw1's FEM dose at P5 misses tw0's by > ``DOSE_BAR``).  Printed: I1h, TABLE rows for
tw1 beside P73's ``_d_comp``.

Not claimed: the T_C wait in FEM, other waits, a dose band, a rate, a recommended wait.  ``K_W`` / ``RHO_C_W``
nominal PMMA, citation owed (W3).
"""

from __future__ import annotations

import os
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
from tests.validation.test_mat15_coil_driven_heating_phantom import H_GEL, K_GEL, N_MAIN
from tests.validation.test_mat15_finite_cylinder_closed_form import H, K, Q, T_INF
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
    _sched_string,
)
from tests.validation.test_mat21_box_duty_closed_form import DT_BOX, HALF_T
from tests.validation.test_mat21_initial_gradient_closed_form import G_GRAD, GC
from tests.validation.test_mat21_pre_equilibration_closed_form import ADJ_PT, P5
from tests.validation.test_mat21_pre_equilibration_closed_form import T_INF as T_INF_P59
from tests.validation.test_mat21_pre_equilibration_duty_dose_closed_form import _dose
from tests.validation.test_mat21_wall_capacity_dose_sweep_closed_form import _d_comp
from tests.validation.test_mat21_wall_capacity_drift_closed_form import (
    LAM1_W,
    WAITS,
    _drift_cap,
    _drift_cap_z,
)
from tests.validation.test_mat21_wall_closed_form import D_WALL, RHO_C_W
from tests.validation.test_mat21_wall_fem import ROBIN, _materials, _robin_facets, _walled_box
from tests.validation.test_mat21_wall_pre_equilibration_closed_form import _drift, _drift_1d
from tests.validation.test_mat21_wall_pre_equilibration_duty25_dose_closed_form import _uniform_w
from tests.validation.test_mat21_wall_pre_equilibration_fem import (
    DT,
    I_ADJ,
    I_P0,
    I_P5,
    PRE_BAR,
    READ_PTS,
    _drift_z,
    _mids,
    _read,
    _series_z,
    _ZField,
)
from tests.validation.test_mat21_wall_transient_closed_form import T_STAR
from tests.validation.test_mat21_wall_transient_fem import MASSLESS, RUNG, _rho_c_field

HIST_LABELS = ("P5", "ADJ_PT")
HIST_PTS = np.array((P5, ADJ_PT), dtype=float)
KEYS = ("ref", "tw0", "tw0_mw")
TABLE_KEYS = ("ref", "tw0")  # item 789's walled series exists for these two only
#: Reproduction controls (§9 item 835): attempts 1 / 2's printed values, bound = half the last
#: printed digit (`20261104T061055Z_MAT-21-P65-W1.log:64`, `20261104T103557Z_MAT-21-P65-W1m.log:66-67`).
REPRO = {"I0h": (2.942439e-02, 5e-9), "Ir": (1.849599e-03, 5e-10), "Imw": (5.708485e-03, 5e-10)}
#: Predicted, printed only (item 835): (I) tw0 against the composite, the drifted-start FEM floor.
PREDICTED_I_TW0_CAP = (5e-3, 8e-3)
#: Item 805's R1 − R3 at t* under the full source (P5, ADJ_PT, P0), K — item 812's brief.
ITEM805_R1_R3 = {"P5": -6.931e-02, "ADJ_PT": -5.664e-02, "P0": -1.042e-02}
WC_PTS = (("P5", I_P5), ("ADJ_PT", I_ADJ), ("P0", I_P0))
#: Item 789's walled D(0) at P5 (t_w = 0, peak scale, P59''s ``_dose``: cem43 of T_INF + θ on
#: ``T_GRID``, no equal-mean normalisation), MAT21-35 TABLE.
ITEM789_D0_P5 = 6.959718e-11
GRID_BAR = 1e-12
#: Item 885 (P74): the wait 1/LAM1_W (P67's module) and its source-off step grid.
TW1 = WAITS[2]
N_W = int(np.ceil(TW1 / DT))
DT_W = TW1 / N_W
P74_KNOB = "FEM_EM_MAT21_P74_WAIT"
HANDOFF_BAR = 1e-12  # K: (b)'s snapshot 0 vs (a)'s final temperature, max |diff| allreduced
#: Predicted, printed only (item 885): (I) tw1 between P68's reference floor and its tw0.
PREDICTED_I_TW1 = (1.8e-03, 5.6e-03)
#: Predicted, printed (P67-W1.log:68, :70): composite drift at P5 at t* and t* + 1/LAM1_W, K.
PREDICTED_DRIFT_P5 = (1.819866e-01, 6.658629e-02)
TW0_FEM_SERIES_P5 = 4.756e-03  # P68-W1.log:91, tw0's FEM/series - 1 at P5, printed beside


def _p74_wait():
    """The P74 knob: unset / '' -> False, '1' -> True, anything else raises (fixture-time)."""
    v = os.environ.get(P74_KNOB, "")
    if v == "":
        return False
    if v == "1":
        return True
    raise ValueError(f"{P74_KNOB}={v!r}: must be unset, empty or '1'")


def _grid_offsets():
    """max |t_m / DT − round(t_m / DT)| over ``T_SW25`` and over ``T_GRID``."""
    sw = max(abs(t / DT - round(t / DT)) for t in T_SW25)
    gr = float(np.max(np.abs(T_GRID / DT - np.arange(T_GRID.size))))
    return sw, gr


def _hist_series(tw=0.0):
    """Walled histories at P5 / ADJ_PT: (ref, tw0 with item 823's ``_drift_cap``, tw0 with
    ``_drift_1d`` (item 789's h_eff form, printed beside), six-face gap).  ``tw`` > 0 (item 885)
    evaluates both drifts at t + tw (the tw1 rows)."""
    uni = _uniform_w(T_SW25)[:, :2]
    dc = np.array([[_drift_cap(tuple(p), t + tw) for p in HIST_PTS] for t in T_GRID])
    d1 = np.array([[_drift_1d(tuple(p), t + tw) for p in HIST_PTS] for t in T_GRID])
    d6 = np.array([[_drift(tuple(p), t + tw) for p in HIST_PTS] for t in T_GRID])
    return uni, uni + dc, uni + d1, float(np.max(np.abs(d6 - d1)))


def test_w0_p74_knob(monkeypatch):
    """W0 (item 885, no mesh): the knob's accepted / refused values, TW1, N_W, DT_W."""
    seen = {}
    for v, want in (("", False), ("1", True)):
        monkeypatch.setenv(P74_KNOB, v)
        seen[repr(v)] = _p74_wait() is want
    monkeypatch.delenv(P74_KNOB, raising=False)
    seen["unset"] = _p74_wait() is False
    for v in ("2", "yes"):
        monkeypatch.setenv(P74_KNOB, v)
        try:
            _p74_wait()
            seen[repr(v)] = False
        except ValueError:
            seen[repr(v)] = True
    monkeypatch.delenv(P74_KNOB, raising=False)
    tl = abs(TW1 * LAM1_W - 1.0)
    grid = abs(N_W * DT_W - TW1)
    _say(f"\nMAT21-41b W0 knob accepted/refused as specified: {seen}; TW1 {TW1:.6f} s, "
         f"|TW1 LAM1_W - 1| {tl:.3e} (bar 1e-12); N_W {N_W} (must be 62), DT_W {DT_W:.6f} s, "
         f"|N_W DT_W - TW1| {grid:.3e} s (bar {1e-9 * DT:.3e})")
    gates = {"knob": all(seen.values()), "tw1_lam": tl <= 1e-12, "n_w": N_W == 62,
             "grid": grid <= 1e-9 * DT}
    _say(f"MAT21-41b W0 gates: {gates}")
    assert all(gates.values()), gates


def test_w0_premises():
    """W0 (no solve): ``DT_BOX`` vs item 682's Δt, the grid, the schedule, the census."""
    t0 = time.perf_counter()
    rel_dt = abs(DT_BOX / DT - 1.0)
    sw, gr = _grid_offsets()
    s = _sched_string(_chi25)
    mult = [_chi25(t) for t in _mids(N_MAIN)]
    dts = np.array([0.5 * DT] * 2 + [DT] * (len(mult) - 2))
    rel_on = abs(float(np.dot(mult, dts)) - DUTY25 * T_STAR) / (DUTY25 * T_STAR)
    _say(f"\nMAT21-37 W0 DT_BOX {DT_BOX:.12f} s vs item 682 dt {DT:.12f} s rel {rel_dt:.3e} "
         f"(bar {GRID_BAR:.0e}); T_SW25 max off-grid {sw:.3e} steps, T_GRID max off-grid {gr:.3e} "
         f"steps; P = 2 HALF_T {2 * HALF_T:.4f} s; T_STAR {T_STAR:.4f} s, N {N_MAIN}")
    _say(f"MAT21-37 W0 schedule {s} (must equal {SCHEDULE_25}); census sum chi dt vs D t* rel "
         f"{rel_on:.3e} (band {S0_BAND:g}); on {int(sum(mult))} of {len(mult)}")
    _say(f"MAT21-37 W0 T_INF {T_INF} (FEM read) vs P59' dose baseline {T_INF_P59}; K == K_GEL "
         f"{K == K_GEL}, H == H_GEL {H == H_GEL}; elapsed {time.perf_counter() - t0:.1f} s")
    assert rel_dt <= GRID_BAR, (DT_BOX, DT)
    assert sw <= GRID_BAR and gr <= GRID_BAR, (sw, gr)
    assert s == SCHEDULE_25, s
    assert rel_on <= S0_BAND, rel_on
    assert K == K_GEL and H == H_GEL and T_INF == T_INF_P59


@pytest.fixture(scope="module")
def runs():
    comm = MPI.COMM_WORLD
    t0 = time.perf_counter()
    wait_on = _p74_wait()  # item 885: refused here, never at import
    msh, ns, nw = _walled_box(RUNG, (False, False, True))
    chi, kf, inter, _ct = _materials(msh)
    rcf = _rho_c_field(msh)
    rc0 = fem.Function(rcf.function_space)  # item 805's R3 form (`..._wall_pre_equilibration_fem.py:206`)
    rc0.x.array[:] = np.where(rcf.x.array == RHO_C_W, MASSLESS * RHO_C_W, rcf.x.array)
    ft = _robin_facets(msh, (A_BOX, B_BOX, C_BOX + D_WALL), (2,))
    _say(f"\nMAT21-37 walled box m = {RUNG} cells {tuple(ns)} N_w = {nw}: n = "
         f"{msh.topology.index_map(3).size_global}; D = {DUTY25}; CN N = {N_MAIN} + 2 BE "
         f"half-steps, {N_MAIN + 1} snapshots per run; batch {HIST_LABELS}; ranks {comm.size}")
    grad = lambda x: T_INF + G_GRAD * x[2]  # noqa: E731  (every cell, gel and wall)
    kw = dict(k=kf, q=Q * chi, h=H_GEL, t_inf=T_INF, t_init=T_INF, robin_tag=ROBIN,
              theta=0.5, startup_half_steps=2, dt=DT, n_steps=N_MAIN,
              snapshot_steps=tuple(range(N_MAIN + 1)), source_schedule=_chi25)
    starts = {"ref": (None, rcf), "tw0": (grad, rcf), "tw0_mw": (grad, rc0)}
    if wait_on:
        starts["tw1"] = (None, rcf)  # expr is call (a)'s temperature, set below
    out = {"E": {}, "sched": {}, "s0": {}, "hist": {}, "hist_ok": {}, "wall": {}, "wc": {},
           "keys": tuple(starts), "wait_on": wait_on}
    res = {}
    prev_mw = None
    for key, (expr, rho_c) in starts.items():
        tq = time.perf_counter()
        if key == "tw1":  # item 885 (a): the source-off wait, P68's kw and rcf
            kw_a = dict(kw, dt=DT_W, n_steps=N_W, snapshot_steps=(),
                        source_schedule=lambda t: 0.0)
            res_a = solve_transient_conduction(msh, ft, t_init_expr=grad, rho_c=rho_c, **kw_a)
            out["E"]["tw1_wait"] = float(res_a.max_balance_residual)
            mult_a = [float(x) for x in res_a.source_multipliers]
            out["wait_mult0"] = len(mult_a) > 0 and all(x == 0.0 for x in mult_a)
            out["wall"]["tw1_wait"] = time.perf_counter() - tq
            _say(f"MAT21-41b tw1 (a) wait: {N_W} steps of {DT_W:.6f} s ({N_W * DT_W:.4f} s), "
                 f"{out['wall']['tw1_wait']:.1f} s; (E) {out['E']['tw1_wait']:.3e} (E_BAND "
                 f"{E_BAND:.0e}); multipliers all 0 {out['wait_mult0']} (n = {len(mult_a)})")
            expr = res_a.temperature
            tq = time.perf_counter()
        res[key] = solve_transient_conduction(msh, ft, t_init_expr=expr, rho_c=rho_c, **kw)
        ts = time.perf_counter()
        if key == "tw1":  # hand-off precondition: (b)'s snapshot 0 == (a)'s final temperature
            loc = float(np.max(np.abs(res[key].snapshots[0].x.array
                                      - res_a.temperature.x.array)))
            out["handoff"] = comm.allreduce(loc, op=MPI.MAX)
            out["handoff_grid"] = abs(N_W * DT_W - TW1)
            _say(f"MAT21-41b hand-off: max |snap0(b) - T(a)| {out['handoff']:.3e} K (bar "
                 f"{HANDOFF_BAR:.0e}); |N_W DT_W - TW1| {out['handoff_grid']:.3e} s (bar "
                 f"{1e-9 * DT:.3e})")
        out["E"][key] = float(res[key].max_balance_residual)
        out["sched"][key] = "".join(str(int(x)) for x in res[key].source_multipliers)
        try:  # item 747's census, run by run (a stale closure would read the wrong duty)
            _census(res[key], DT, 2, DUTY25 * T_STAR, f"wall pre-equilibration d25 {key}")
            out["s0"][key] = True
        except AssertionError:
            out["s0"][key] = False
        hist = np.empty((N_MAIN + 1, HIST_PTS.shape[0]))
        ok_all = [True] * HIST_PTS.shape[0]
        for n in range(N_MAIN + 1):
            hist[n], ok = _read(res[key].snapshots[n], HIST_PTS)
            ok_all = [a and b for a, b in zip(ok_all, ok)]
        if key == "tw0_mw":  # kept for the (I-mw25) ringing print at N_MAIN - 1
            prev_mw = res[key].snapshots[N_MAIN - 1]
        res[key].snapshots.clear()
        if key in ("tw0", "tw0_mw", "tw1"):
            out["wc"][key], ok = _read(res[key].temperature, READ_PTS)
            out["hist_ok"][key + "_wc"] = ok
        out["hist"][key], out["hist_ok"][key] = hist, ok_all
        out["wall"][key] = time.perf_counter() - tq
        _say(f"MAT21-37 run {key}: solve {ts - tq:.1f} s, reads {time.perf_counter() - ts:.1f} s; "
             f"(E) {out['E'][key]:.3e} (E_BAND {E_BAND:.0e}); source_power "
             f"{float(res[key].source_power):.9e} W")

    tq = time.perf_counter()
    zf = _ZField()
    s_c = zf(msh, _series_z(RHO_C_W), T_STAR)  # m = 0, also the continuous comparand
    s_d = fem.Function(s_c.function_space)
    s_d.x.array[:] = s_c.x.array
    for i, tm in enumerate(T_SW25):
        if i:
            s_d.x.array[:] += (-1.0) ** i * zf(msh, _series_z(RHO_C_W), T_STAR - tm).x.array
    d_0 = zf(msh, _drift_cap_z, T_STAR)  # item 835: the composite (capacitive-wall) drift
    d_0h = zf(msh, _drift_z, T_STAR)  # the h_eff drift: exact on the massless twin ((Lim-c))

    def _mw_sum(t_ev):
        """Σ(−1)^m ``_series_z(MASSLESS·RHO_C_W)``(z, t_ev − t_m) over ``T_SW25`` (t_m < t_ev)."""
        acc = None
        for i, tm in enumerate(T_SW25):
            if tm < t_ev:
                f = zf(msh, _series_z(MASSLESS * RHO_C_W), t_ev - tm)
                if acc is None:
                    acc = f
                else:
                    acc.x.array[:] += (-1.0) ** i * f.x.array
        return acc

    s_mw = _mw_sum(T_STAR)
    t_prev = (N_MAIN - 1) * DT
    s_mw_prev = _mw_sum(t_prev)
    d_prev = zf(msh, _drift_z, t_prev)
    out["pre"] = comm.allreduce(zf.worst, op=MPI.MAX)
    _say(f"MAT21-37 {len(T_SW25) + 1} z-only comparands on the unique dof z ({zf.n_unique} values "
         f"per rank max): {time.perf_counter() - tq:.1f} s; precondition worst rel "
         f"{out['pre']:.3e} (bar {PRE_BAR:.0e})")

    th = {key: r.temperature - T_INF for key, r in res.items()}
    out["Ir"] = _rel(msh, th["ref"], s_d, inter)
    out["I0"] = _rel(msh, th["tw0"], s_d + d_0, inter)
    out["Ic"] = _rel(msh, th["tw0"], s_c + d_0, inter)
    out["Imw"] = _rel(msh, th["tw0_mw"], s_mw + d_0h, inter)
    out["I0h"] = _rel(msh, th["tw0"], s_d + d_0h, inter)
    out["Imw_prev"] = _rel(msh, prev_mw - T_INF, s_mw_prev + d_prev, inter)
    if wait_on:  # item 885's comparands at t* + TW1
        d_1 = zf(msh, _drift_cap_z, T_STAR + TW1)
        d_1h = zf(msh, _drift_z, T_STAR + TW1)
        out["pre"] = comm.allreduce(zf.worst, op=MPI.MAX)
        out["I1"] = _rel(msh, th["tw1"], s_d + d_1, inter)
        out["I1h"] = _rel(msh, th["tw1"], s_d + d_1h, inter)
        out["I1c0"] = _rel(msh, th["tw1"], s_d + d_0, inter)
        lo1, hi1 = PREDICTED_I_TW1
        _say(f"MAT21-41b (I) tw1 vs composite drift at t* + 1/LAM1_W {out['I1']:.6e} (bound "
             f"{INTERIOR_BAND:.0e}, ASSERTED; predicted [{lo1:.1e}, {hi1:.1e}], printed, inside "
             f"{lo1 <= out['I1'] <= hi1}); I1h vs h_eff drift {out['I1h']:.6e} (printed); "
             f"ceiling I1h / I1 {out['I1h'] / out['I1']:.3f}x; control (c-tw) vs the t_w = 0 "
             f"comparand {out['I1c0']:.4e} (asserted > {INTERIOR_BAND:.0e})")
        dp = (_drift_cap(tuple(P5), T_STAR), _drift_cap(tuple(P5), T_STAR + TW1))
        _say(f"MAT21-41b composite drift at P5: t* {dp[0]:.6e} K, t* + 1/LAM1_W {dp[1]:.6e} K "
             f"(predicted {PREDICTED_DRIFT_P5[0]:.6e} / {PREDICTED_DRIFT_P5[1]:.6e}, printed)")
    lo, hi = PREDICTED_I_TW0_CAP
    _say(f"MAT21-37 (I) t_w = 0 vs composite drift {out['I0']:.6e} (bound {INTERIOR_BAND:.0e}, "
         f"ASSERTED; predicted [{lo:.0e}, {hi:.0e}], printed, inside {lo <= out['I0'] <= hi}); "
         f"I0h vs h_eff drift {out['I0h']:.6e} (printed; reproduction control)")
    _say(f"MAT21-37 ceiling (printed): I0h / (I-mw25) {out['I0h'] / out['Imw']:.3f}x; "
         f"I0 / (I-mw25) {out['I0'] / out['Imw']:.3f}x")
    _say(f"MAT21-37 (I-mw25) tw0_mw interior rel L2 vs its exact massless-wall comparand "
         f"{out['Imw']:.6e} (bound {INTERIOR_BAND:.0e}); whole "
         f"{_rel(msh, th['tw0_mw'], s_mw + d_0h):.6e} (printed)")
    _say(f"MAT21-37 (I-mw25) ringing print: snapshot N_MAIN - 1 (t {t_prev:.4f} s) "
         f"{out['Imw_prev']:.6e} vs N_MAIN {out['Imw']:.6e}, ratio "
         f"{out['Imw_prev'] / out['Imw']:.4f} (printed, no band)")
    _say(f"MAT21-37 (I) reference: interior rel L2 {out['Ir']:.6e}; t_w = 0: {out['I0']:.6e}; "
         f"bound {INTERIOR_BAND:.0e}; whole {_rel(msh, th['ref'], s_d):.6e} / "
         f"{_rel(msh, th['tw0'], s_d + d_0):.6e} (printed); control, continuous comparand on "
         f"t_w = 0: {out['Ic']:.4e} (asserted > {INTERIOR_BAND:.0e})")
    _say(f"MAT21-37 elapsed {time.perf_counter() - t0:.1f} s")
    return out


def test_all_gates(runs):
    """Collect every gate before asserting (item 668), so one red hides no other."""
    r = runs
    gates = {"valid_points": all(all(o) for o in r["hist_ok"].values()),
             "precondition_unique_z": r["pre"] <= PRE_BAR,
             "I_ref": r["Ir"] <= INTERIOR_BAND,
             "I_tw0": r["I0"] <= INTERIOR_BAND,
             "I_mw25_tw0_mw": r["Imw"] <= INTERIOR_BAND,
             "control_continuous_misses_I": r["Ic"] > INTERIOR_BAND}
    for name, (val, bnd) in REPRO.items():  # item 835's reproduction controls (asserted)
        _say(f"MAT21-37 reproduction {name}: {r[name]:.9e} vs {val:.6e}, |diff| "
             f"{abs(r[name] - val):.3e} (bound {bnd:.0e}, ASSERTED)")
        gates[f"repro_{name}"] = abs(r[name] - val) <= bnd
    for key in r["keys"]:
        gates[f"E_{key}"] = r["E"][key] <= E_BAND
        gates[f"S0_census_{key}"] = r["s0"][key]
        gates[f"S0_schedule_{key}"] = r["sched"][key] == SCHEDULE_25
        _say(f"MAT21-37 (S0) {key}: schedule {r['sched'][key]} (must equal {SCHEDULE_25})")
    for key in r["keys"]:
        for j, lab in enumerate(HIST_LABELS):
            dose, quad, rel, _swap, miss = _q_check_g(T_INF + r["hist"][key][:, j], T_GRID)
            _say(f"MAT21-37 (Q-wf) {key} {lab}: D {dose:.6e} s, quad {quad:.6e} s, rel {rel:.3e} "
                 f"(bar {DOSE_BAR:g}, ASSERTED); R swap miss {miss:.4e} (printed)")
            gates[f"Q_wf_{key}_{lab}"] = rel <= DOSE_BAR
    j = HIST_LABELS.index("P5")
    d_0 = _dose(r["hist"]["tw0"][:, j])
    d_r = _dose(r["hist"]["ref"][:, j])
    miss = abs(d_0 / d_r - 1.0)
    _say(f"MAT21-37 control P5: t_w = 0 FEM D {d_0:.6e} s vs reference {d_r:.6e} s, miss "
         f"{miss:.4e} (asserted > {DOSE_BAR:g})")
    gates["control_tw0_dose_misses_ref_P5"] = miss > DOSE_BAR
    if r["wait_on"]:  # item 885 (P74)
        d_1 = _dose(r["hist"]["tw1"][:, j])
        miss1 = abs(d_1 / d_0 - 1.0)
        _say(f"MAT21-41b control (c-D) P5: tw1 FEM D {d_1:.6e} s vs tw0 {d_0:.6e} s, miss "
             f"{miss1:.4e} (asserted > {DOSE_BAR:g})")
        gates["I_tw1"] = r["I1"] <= INTERIOR_BAND
        gates["E_tw1_wait"] = r["E"]["tw1_wait"] <= E_BAND
        gates["S0_tw1_wait_mult0"] = r["wait_mult0"]
        gates["handoff"] = r["handoff"] <= HANDOFF_BAR and r["handoff_grid"] <= 1e-9 * DT
        gates["control_c_tw"] = r["I1c0"] > INTERIOR_BAND
        gates["control_c_D"] = miss1 > DOSE_BAR
    _say(f"MAT21-37 gates: {gates}")
    assert all(gates.values()), (gates, {k: v for k, v in r.items() if k not in ("hist",)})


def test_printed_dose_table(runs):
    """The record (no band); asserts only that it was formed."""
    uni, tw0, tw0h, d16 = _hist_series()
    ser = {"ref": uni, "tw0": tw0}
    _say(f"MAT21-37 SIDE |six-face _drift - _drift_1d| on T_GRID at P5/ADJ_PT: max {d16:.3e} K "
         f"({d16 / GC:.3e} g C) (printed)")
    _say("MAT21-37 TABLE point | run | D FEM s | D walled series (tw0: item 823 _drift_cap) s | "
         "FEM/series-1 | FEM - series s | FEM D/D_ref-1 | series D/D_ref-1 | max |hist FEM - "
         "series| K   [peak scale: P59' _dose, cem43 of T_INF + theta on T_GRID]")
    for j, lab in enumerate(HIST_LABELS):
        d_h = {key: _dose(runs["hist"][key][:, j]) for key in TABLE_KEYS}
        d_s = {key: _dose(ser[key][:, j]) for key in TABLE_KEYS}
        for key in TABLE_KEYS:
            hg = float(np.max(np.abs(runs["hist"][key][:, j] - ser[key][:, j])))
            _say(f"MAT21-37 TABLE {lab} | {key} | {d_h[key]:.6e} | {d_s[key]:.6e} | "
                 f"{d_h[key] / d_s[key] - 1.0:+.3e} | {d_h[key] - d_s[key]:+.3e} | "
                 f"{d_h[key] / d_h['ref'] - 1.0:+.6e} | {d_s[key] / d_s['ref'] - 1.0:+.6e} | "
                 f"{hg:.3e}")
            assert np.isfinite(d_h[key] / d_s[key])
        d_hh = _dose(tw0h[:, j])
        _say(f"MAT21-37 TABLE {lab} | tw0 h_eff (_drift_1d) series | {d_hh:.6e} | FEM/h_eff-1 "
             f"{d_h['tw0'] / d_hh - 1.0:+.3e} | composite/h_eff-1 {d_s['tw0'] / d_hh - 1.0:+.3e} "
             f"| max |hist FEM - h_eff series| "
             f"{float(np.max(np.abs(runs['hist']['tw0'][:, j] - tw0h[:, j]))):.3e} K")
    d0 = _dose(runs["hist"]["tw0"][:, 0])
    _say(f"MAT21-37 P5 t_w = 0: FEM D {d0:.6e} s beside item 789's walled D(0) "
         f"{ITEM789_D0_P5:.6e} s (six-face drift), signed gap {d0 - ITEM789_D0_P5:+.3e} s, "
         f"rel {d0 / ITEM789_D0_P5 - 1.0:+.3e}; one rung, one dt, no rate")
    _say("MAT21-37 WALLCAP point | tw0 - tw0_mw at t* K (FEM vs FEM, D = 0.25) | item 805 R1 - R3 "
         "K (full source) | ratio   [the wall's heat capacity, printed, no band]")
    for lab, i in WC_PTS:
        dwc = float(runs["wc"]["tw0"][i] - runs["wc"]["tw0_mw"][i])
        ref805 = ITEM805_R1_R3[lab]
        _say(f"MAT21-37 WALLCAP {lab} | {dwc:+.6e} | {ref805:+.3e} | {dwc / ref805:.4f}")
        assert np.isfinite(dwc)
    if runs["wait_on"]:  # item 885 (P74): the tw1 rows
        _u1, tw1c, tw1h, _d = _hist_series(TW1)
        _say("MAT21-41b TABLE point | tw1 D FEM s | composite series D s | P73 _d_comp s | "
             "|history/_d_comp - 1| | h_eff series D s | FEM/composite-1 (tw0 P5 "
             f"{TW0_FEM_SERIES_P5:+.3e}) | FEM/h_eff-1 | FEM D/D_ref-1 | composite D/D_ref-1")
        for j, lab in enumerate(HIST_LABELS):
            d_f = _dose(runs["hist"]["tw1"][:, j])
            d_c = _dose(tw1c[:, j])
            d_p = _d_comp(j, TW1)
            d_hh = _dose(tw1h[:, j])
            d_rf = _dose(runs["hist"]["ref"][:, j])
            d_rs = _dose(uni[:, j])
            _say(f"MAT21-41b TABLE {lab} | {d_f:.6e} | {d_c:.6e} | {d_p:.6e} | "
                 f"{abs(d_c / d_p - 1.0):.3e} | {d_hh:.6e} | {d_f / d_c - 1.0:+.3e} | "
                 f"{d_f / d_hh - 1.0:+.3e} | {d_f / d_rf - 1.0:+.6e} | {d_c / d_rs - 1.0:+.6e}")
            assert np.isfinite(d_f / d_c)
    _say("MAT21-37 wall time per run (solve + reads): "
         + ", ".join(f"{k} {v:.1f} s" for k, v in runs["wall"].items()))
