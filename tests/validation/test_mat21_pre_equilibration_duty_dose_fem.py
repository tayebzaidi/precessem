"""`MAT-21` step 33 — the wait in a dose under D = 0.25 on the FEM history, the FEM half of P59
(§9 item 747, §10 chain P, P60, provisional under rule (6); brief
`docs/planning/chunks/MAT-21.md` "§9 item 747", Review ruling 2026-11-02 03:00).

Sets no band and moves none; every band is imported.  Two runs on ``RUNGS[-1]``, reference (from
uniform T∞) and t_w = 0 (from T∞ + g·z), P52's call with ``source_schedule=_chi25`` and
``snapshot_steps=tuple(range(N_MAIN + 1))``.  Each run's snapshots are read on one fixed batch,
``HIST_PTS`` = P5, ``ADJ_PT`` (item 537), then cleared, so one run's snapshots are held at a
time.  The comparand is item 753's uniform sum at ``T_SW25`` plus the drift, as the ruling amends
"item 746's superposition".  Real build.

Asserted: (I) ≤ ``INTERIOR_BAND`` at t* on both runs, the reference's comparand being the sum
without the drift; (E) ≤ ``E_BAND`` on both; (S0) ``_census`` at ``S0_BAND`` on the recorded
multipliers, and the schedule string; (Q-wf) ``cem43`` against ``_q_check_g``'s quadrature at
``DOSE_BAR`` on all four histories.  W0, no solve: the schedule, the grid, and item 753's
t_w = 0 rows reproduced within one unit in the print's last digit.  Negative controls: the
t_w = 0 FEM dose at P5 misses the reference run's by > ``DOSE_BAR``; the continuous superposition
misses (I) by > ``INTERIOR_BAND``.  Printed (**no band**): D at P5 and ``ADJ_PT``, FEM vs item
753, and the signed gap, one rung, one Δt, no rate.

Not claimed: a dose band, a rate, a recommended wait.
"""

from __future__ import annotations

import time

import numpy as np
import pytest
from dolfinx import fem
from dolfinx import mesh as dmesh
from mpi4py import MPI

from fem_em_solver.thermal import solve_transient_conduction
from tests.validation.test_imp7_duty_cycle_fem import S0_BAND
from tests.validation.test_mat10_sar_to_heat_sphere import ROBIN_TAG, _robin_tags, _say
from tests.validation.test_mat12_transient_sar_heating_sphere import E_BAND, RHO_C
from tests.validation.test_mat13_transient_perfusion_sphere import INTERIOR_BAND
from tests.validation.test_mat15_coil_driven_heating_phantom import H_GEL, K_GEL, N_MAIN, T_INF
from tests.validation.test_mat15_finite_cylinder_fem import _rel
from tests.validation.test_mat17_duty_cycle_fem import _census
from tests.validation.test_mat21_box_closed_form import A_BOX, B_BOX, C_BOX, Q
from tests.validation.test_mat21_box_dose_closed_form import DOSE_BAR, T_GRID
from tests.validation.test_mat21_box_dt_duty125_dose_closed_form import _q_check_g
from tests.validation.test_mat21_box_duty25_dose_closed_form import (
    DUTY25,
    SCHEDULE_25,
    T_SW25,
    _chi25,
    _sched_string,
)
from tests.validation.test_mat21_box_duty_closed_form import DT_BOX
from tests.validation.test_mat21_box_fem import (
    BASE,
    RUNGS,
    T_STAR,
    _exact,
    _interior_chi,
    _trans_ex,
)
from tests.validation.test_mat21_initial_gradient_closed_form import G_GRAD
from tests.validation.test_mat21_initial_gradient_fem import _grad_ex
from tests.validation.test_mat21_pre_equilibration_closed_form import ADJ_PT, P5
from tests.validation.test_mat21_pre_equilibration_duty_dose_closed_form import (
    _dose,
    _drifts,
    _uniform,
)
from tests.validation.test_mat21_pre_equilibration_fem import _read

HIST_LABELS = ("P5", "ADJ_PT")
HIST_PTS = np.array((P5, ADJ_PT), dtype=float)
KEYS = ("ref", "tw0")
#: Item 753's D / D_ref − 1, D = 0.25, t_w = 0, peak, `20261102T103526Z_MAT-21-step32-W1.log`.
ITEM753_RATIO = {"P5": 3.772126e-01, "ADJ_PT": 3.957192e-01}
REPRO_BOUND = 1e-6  # one unit in the print's last digit
_CACHE = {}


def _series(j, key):
    """Item 753's history at ``HIST_PTS[j]``: the uniform sum, plus the drift on "tw0"."""
    if not _CACHE:
        for i, p in enumerate((P5, ADJ_PT)):
            _CACHE[i] = (_uniform(p, T_SW25), _drifts(p, 0.0))
    sum_h, drift = _CACHE[j]
    return sum_h if key == "ref" else sum_h + drift


@pytest.fixture(scope="module")
def runs():
    comm = MPI.COMM_WORLD
    t0 = time.perf_counter()
    m = RUNGS[-1]
    ns = [b * m for b in BASE]
    msh = dmesh.create_box(
        comm,
        [np.array([-A_BOX, -B_BOX, -C_BOX]), np.array([A_BOX, B_BOX, C_BOX])],
        ns,
        cell_type=dmesh.CellType.tetrahedron,
    )
    ft = _robin_tags(msh)
    chi = _interior_chi(msh)
    _say(
        f"\nMAT21-33 rung m = {m} {tuple(ns)}: cells {msh.topology.index_map(3).size_global}; "
        f"D = {DUTY25}; CN N = {N_MAIN} + 2 BE half-steps, {N_MAIN + 1} snapshots per run; "
        f"batch {HIST_LABELS}; ranks {comm.size}"
    )

    kw = dict(
        k=K_GEL,
        rho_c=RHO_C,
        q=Q,
        h=H_GEL,
        t_inf=T_INF,
        t_init=T_INF,
        robin_tag=ROBIN_TAG,
        dt=T_STAR / N_MAIN,
        n_steps=N_MAIN,
        theta=0.5,
        startup_half_steps=2,
        degree=1,
        snapshot_steps=tuple(range(N_MAIN + 1)),
        source_schedule=_chi25,
    )
    starts = {"ref": None, "tw0": lambda x: T_INF + G_GRAD * x[2]}
    out = {"E": {}, "sched": {}, "hist": {}, "hist_ok": {}}
    res = {}
    for key, expr in starts.items():
        tq = time.perf_counter()
        res[key] = solve_transient_conduction(msh, ft, t_init_expr=expr, **kw)
        ts = time.perf_counter()
        out["E"][key] = float(res[key].max_balance_residual)
        out["sched"][key] = "".join(str(int(x)) for x in res[key].source_multipliers)
        _census(res[key], DT_BOX, 2, DUTY25 * T_STAR, f"pre-equilibration d25 {key}")

        # One batch, one order, on every snapshot; read and clear.
        hist = np.empty((N_MAIN + 1, HIST_PTS.shape[0]))
        ok_all = [True] * HIST_PTS.shape[0]
        for n in range(N_MAIN + 1):
            hist[n], ok = _read(res[key].snapshots[n], HIST_PTS)
            ok_all = [a and b for a, b in zip(ok_all, ok)]
        res[key].snapshots.clear()
        out["hist"][key], out["hist_ok"][key] = hist, ok_all
        _say(
            f"MAT21-33 run {key}: solve {ts - tq:.1f} s, reads {time.perf_counter() - ts:.1f} s; "
            f"(E) {out['E'][key]:.3e} (E_BAND {E_BAND:.0e})"
        )

    tq = time.perf_counter()
    ex_c, npts = _exact(msh, _trans_ex(T_STAR))  # m = 0, also the continuous series
    ex_d = fem.Function(ex_c.function_space)
    ex_d.x.array[:] = ex_c.x.array
    for i, tm in enumerate(T_SW25):
        if i:
            g_m, _ = _exact(msh, _trans_ex(T_STAR - tm))
            ex_d.x.array[:] += (-1.0) ** i * g_m.x.array
    ex_0, _ = _exact(msh, _grad_ex(T_STAR))
    _say(
        f"MAT21-33 {len(T_SW25) + 1} series interpolations ({npts} P2 points): "
        f"{time.perf_counter() - tq:.1f} s"
    )

    th = {key: r.temperature - T_INF for key, r in res.items()}
    out["Ir"] = _rel(msh, th["ref"], ex_d, chi)
    out["I0"] = _rel(msh, th["tw0"], ex_d + ex_0, chi)
    out["Ic"] = _rel(msh, th["tw0"], ex_c + ex_0, chi)
    _say(
        f"MAT21-33 (I) reference: interior rel L2 {out['Ir']:.6e}; t_w = 0: {out['I0']:.6e}; "
        f"bound {INTERIOR_BAND:.0e}; control, continuous superposition on t_w = 0: "
        f"{out['Ic']:.4e} (asserted > {INTERIOR_BAND:.0e})"
    )
    _say(f"MAT21-33 elapsed {time.perf_counter() - t0:.1f} s")
    return out


def test_w0_schedule_grid_and_item_753():
    """W0 (numpy): the schedule, the grid, and item 753's rows, with no solve."""
    s = _sched_string(_chi25)
    assert s == SCHEDULE_25, s
    assert S0_BAND == 1e-12  # `_census` carries the same bar as a literal
    assert N_MAIN == 40 and RUNGS[-1] == 6 and T_GRID.size == N_MAIN + 1
    assert abs(T_STAR / N_MAIN / DT_BOX - 1.0) <= 1e-14
    assert abs(T_GRID[-1] / T_STAR - 1.0) <= 1e-14
    for j, lab in enumerate(HIST_LABELS):
        ratio = _dose(_series(j, "tw0")) / _dose(_series(j, "ref")) - 1.0
        diff = abs(ratio - ITEM753_RATIO[lab])
        _say(
            f"MAT21-33 W0 {lab}: series D / D_ref - 1 {ratio:+.6e} vs item 753 "
            f"{ITEM753_RATIO[lab]:+.6e}, |delta| {diff:.3e} (bound {REPRO_BOUND:g})"
        )
        assert diff <= REPRO_BOUND, (lab, ratio)
    _say(f"MAT21-33 W0 schedule {s}; {T_GRID.size} rows, DT {DT_BOX:.4f} s")


def test_preconditions_every_history_read(runs):
    for key, ok in runs["hist_ok"].items():
        assert all(ok), (key, ok)


def test_s0_schedule_on_the_recorded_multipliers(runs):
    # `_census` asserted the sum at S0_BAND in the fixture, run by run.
    for key, s in runs["sched"].items():
        _say(f"MAT21-33 (S0) {key}: schedule {s} (must equal {SCHEDULE_25})")
        assert s == SCHEDULE_25, (key, s)


def test_i_both_runs(runs):
    assert runs["Ir"] <= INTERIOR_BAND, runs["Ir"]
    assert runs["I0"] <= INTERIOR_BAND, runs["I0"]


def test_e_both_runs(runs):
    for key, e in runs["E"].items():
        assert e <= E_BAND, (key, e)


def test_q_wf_on_all_four_histories(runs):
    rows = []
    for key, hist in runs["hist"].items():
        for j, lab in enumerate(HIST_LABELS):
            dose, quad, rel, _swap, miss = _q_check_g(T_INF + hist[:, j], T_GRID)
            _say(
                f"MAT21-33 (Q-wf) {key} {lab}: D {dose:.6e} s, quad {quad:.6e} s, rel {rel:.3e} "
                f"(bar {DOSE_BAR:g}, ASSERTED); R swap miss {miss:.4e} (printed)"
            )
            rows.append((key, lab, rel))
    for key, lab, rel in rows:
        assert rel <= DOSE_BAR, (key, lab, rel)


def test_control_gradient_dose_misses_the_reference_at_p5(runs):
    j = HIST_LABELS.index("P5")
    d_0 = _dose(runs["hist"]["tw0"][:, j])
    d_r = _dose(runs["hist"]["ref"][:, j])
    miss = abs(d_0 / d_r - 1.0)
    _say(
        f"MAT21-33 control P5: t_w = 0 FEM D {d_0:.6e} s vs reference {d_r:.6e} s, miss "
        f"{miss:.4e} (asserted > {DOSE_BAR:g})"
    )
    assert miss > DOSE_BAR, miss


def test_control_continuous_superposition_misses_i(runs):
    assert runs["Ic"] > INTERIOR_BAND, runs["Ic"]


def test_printed_dose_table(runs):
    """The record (no band); asserts only that it was formed."""
    _say(
        "MAT21-33 TABLE point | run | D FEM s | D item 753 s | FEM/series-1 | FEM D/D_ref-1 | "
        "series D/D_ref-1"
    )
    for j, lab in enumerate(HIST_LABELS):
        d_h = {key: _dose(runs["hist"][key][:, j]) for key in KEYS}
        d_s = {key: _dose(_series(j, key)) for key in KEYS}
        for key in KEYS:
            _say(
                f"MAT21-33 TABLE {lab} | {key} | {d_h[key]:.6e} | {d_s[key]:.6e} | "
                f"{d_h[key] / d_s[key] - 1.0:+.3e} | {d_h[key] / d_h['ref'] - 1.0:+.6e} | "
                f"{d_s[key] / d_s['ref'] - 1.0:+.6e}"
            )
            assert np.isfinite(d_h[key] / d_s[key])
