"""`IMP-7` step 9 — the FEM-vs-series dose gap around the curved adiabatic hole, separating
h from Δt (§9 item 710, §10 chain I, I48c; brief `docs/planning/chunks/IMP-7.md` "§9 item
710").

Sets no band and moves none; every band is imported.  Item 700's fixture and read route,
imported from ``test_imp7_duty_dose_fem.py``, untouched, on two new runs to 4 t\\*, selected
by ``FEM_EM_IMP7_S9`` because together they risk the 600 s ceiling:

- ``h``: ``LADDER[1]``, Δt = ``DT``, 160 steps, 161 reads;
- ``dt``: ``LADDER[0]``, Δt = ``DT`` / 2, 320 steps, read at the even steps only, the same
  161 rows of ``T_GRID``.

Asserted, on each run: (T) interior rel L2 at 4 t\\* ≤ ``INTERIOR_BAND``; (E) ≤ ``E_BAND``;
(S0) census ≤ ``S0_BAND``; every point valid; max T < 43 °C; (Grid) the run's read times
equal ``T_GRID`` to 1e-12 relative; the continuous series at ``T_END`` misses >
``INTERIOR_BAND``.  Printed (records, **no band**): the gaps of the run beside ``ITEM_700``.
Real build.

Not claimed: a rate, a dose band, a converged dose, a Δt recommendation, a metal core.
"""

from __future__ import annotations

import os
import time

import numpy as np
import pytest
from mpi4py import MPI

from fem_em_solver.thermal.dose import cem43, running_average

from tests.validation.test_imp7_duty_cycle_fem import DUTY, PERIOD, S0_BAND, _chi, _duty_fn
from tests.validation.test_imp7_duty_dose_closed_form import (
    DT,
    LABELS,
    PTS,
    SCALES,
    T_END,
    T_GRID,
    _duty_hist,
)
from tests.validation.test_imp7_duty_dose_fem import (
    N_RUN,
    ON_TIME_END,
    STARTUP,
    _census,
    _err,
    _probe,
    _shell_fixture,
)
from tests.validation.test_imp7_implant_heating_closed_form import T_INF
from tests.validation.test_imp7_transient_closed_form import _modes
from tests.validation.test_imp7_transient_fem import LADDER, _series_fn, _transient
from tests.validation.test_mat12_transient_sar_heating_sphere import E_BAND, N_MAIN
from tests.validation.test_mat13_transient_perfusion_sphere import INTERIOR_BAND
from tests.validation.test_mat20_running_average import W_B2
from tests.validation.test_mat20_tuned_heating_dose import T_REF

#: run -> (res, dt, n_steps, the step stride that lands on T_GRID)
RUNS = {"h": (LADDER[1], DT, N_RUN, 1), "dt": (LADDER[0], 0.5 * DT, 2 * N_RUN, 2)}
#: Row 1, copied from `20261031T020326Z_IMP-7-step8-W1-r2.log:123–125`, never re-run: point ->
#: (peak-scale dose rel, equal-mean dose rel, six-minute peak Δ / K, max_n|Δθ| / K).
ITEM_700 = {"pole": (1.245e-04, 2.481e-04, -8.242e-05, 1.285e-03),
            "equator": (1.795e-05, 3.556e-05, -7.155e-06, 2.181e-04),
            "outer pole": (2.605e-06, 4.630e-06, +4.165e-06, 6.484e-05)}
GRID_BAR = 1e-12  # the brief's (Grid)


def _say(msg):
    if MPI.COMM_WORLD.rank == 0:
        print(msg, flush=True)


def _census_dt(mult, dt):
    """`_census` at any Δt: the same two half-steps, of that Δt."""
    dts = [0.5 * dt] * STARTUP + [dt] * (len(mult) - STARTUP)
    tot = float(np.dot(np.asarray(mult), np.asarray(dts)))
    return abs(tot - ON_TIME_END) / ON_TIME_END, tot


def _mids(dt, n_steps):
    """The schedule's evaluation times (`transient.py:352–357`)."""
    return ([(i + 0.5) * (0.5 * dt) for i in range(STARTUP)]
            + [(STARTUP // 2 + i + 0.5) * dt for i in range(n_steps - STARTUP // 2)])


def _read_times(dt, n_steps, stride):
    return dt * np.arange(0, n_steps + 1, stride)


def _gaps(th_h, th_s):
    """The four records of one run at each of `LABELS`."""
    W = float(W_B2)
    rows = {}
    for j, name in enumerate(LABELS):
        d_h = cem43(T_GRID, T_INF + th_h[:, j])
        d_s = cem43(T_GRID, T_INF + th_s[:, j])
        de_h = cem43(T_GRID, T_INF + th_h[:, j] / DUTY)
        de_s = cem43(T_GRID, T_INF + th_s[:, j] / DUTY)
        _a, mx_h, _t = running_average(T_GRID, th_h[:, j], window=W)
        _a, mx_s, _t = running_average(T_GRID, th_s[:, j], window=W)
        rows[name] = (abs(d_h / d_s - 1.0), abs(de_h / de_s - 1.0), mx_h - mx_s,
                      float(np.max(np.abs(th_h[:, j] - th_s[:, j]))))
    return rows


@pytest.fixture(scope="module")
def run():
    key = os.environ.get("FEM_EM_IMP7_S9")
    if key not in RUNS:
        pytest.skip("FEM_EM_IMP7_S9 selects the run: h or dt")
    res, dt, n_steps, stride = RUNS[key]
    t0 = time.perf_counter()
    out = {"key": key}
    sub, sft, q_f, t_build, parent = _shell_fixture(res)
    _say(f"\nIMP7-9 run {key}: res={res} sub_cells={sub.topology.index_map(3).size_global} "
         f"build={t_build:.1f}s; ranks {sub.comm.size}; dt {dt:.5f} s, {n_steps} steps to "
         f"{T_END:.4f} s, reads every {stride} step(s)")
    steps = tuple(range(0, n_steps + 1, stride))
    cn, ts = _transient(sub, sft, q_f, dt=dt, n_steps=n_steps, theta=0.5,
                        startup_half_steps=STARTUP, source_schedule=_chi(PERIOD, DUTY),
                        snapshot_steps=steps)
    out["E"] = float(cn.max_balance_residual)
    out["S0"], tot = _census_dt(cn.source_multipliers, dt)
    _say(f"IMP7-9 CN solve {ts:.1f}s; (E) {out['E']:.3e} (E_BAND {E_BAND:.0e}); (S0) sum chi*dt "
         f"{tot:.9f} s vs {ON_TIME_END:.9f} s rel {out['S0']:.3e} (bar {S0_BAND:.0e}); "
         f"{len(cn.source_multipliers)} (half-)steps")

    tq = time.perf_counter()
    th_h = np.empty((len(steps), PTS.shape[0]))
    ok_all = np.ones(PTS.shape[0], dtype=bool)
    for row, k in enumerate(steps):
        th, ok = _probe(cn.snapshots[k], PTS, parent)
        th_h[row] = th
        ok_all &= ok
    out["th"], out["ok"] = th_h, ok_all
    out["grid"] = float(np.max(np.abs(_read_times(dt, n_steps, stride)[1:] / T_GRID[1:] - 1.0)))
    _say(f"IMP7-9 reads: {len(steps)} batches of {PTS.shape[0]} in "
         f"{time.perf_counter() - tq:.1f}s; valid {ok_all.tolist()}; (Grid) {out['grid']:.3e} "
         f"(bar {GRID_BAR:.0e})")

    tq = time.perf_counter()
    e_all, out["T"] = _err(sub, cn.snapshots[n_steps], _duty_fn(T_END))
    _say(f"IMP7-9 (T) at 4 t*: interior rel L2 {out['T']:.4e} (INTERIOR_BAND {INTERIOR_BAND:.0e}); "
         f"whole {e_all:.4e} (printed); {time.perf_counter() - tq:.1f}s")
    tq = time.perf_counter()
    out["cC"] = _err(sub, cn.snapshots[n_steps], _series_fn(T_END, "neumann"))[1]
    _say(f"IMP7-9 control (c-cont) continuous series at T_END: interior miss {out['cC']:.4e} "
         f"(asserted > {INTERIOR_BAND:.0e}); {time.perf_counter() - tq:.1f}s")
    cn.snapshots.clear()

    tq = time.perf_counter()
    out["th_s"] = _duty_hist(_modes())
    _say(f"IMP7-9 series history {out['th_s'].shape} in {time.perf_counter() - tq:.1f}s; "
         f"fixture elapsed {time.perf_counter() - t0:.1f}s")
    return out


def test_imp7_9_w0_grids_and_both_censuses():
    """W0 (numpy): the two grids and both census arithmetics, with no solve."""
    assert T_GRID.size == N_RUN + 1 and N_MAIN == N_RUN // 4
    chi = _chi(PERIOD, DUTY)
    for key, (_res, dt, n_steps, stride) in RUNS.items():
        grid = float(np.max(np.abs(_read_times(dt, n_steps, stride)[1:] / T_GRID[1:] - 1.0)))
        mult = [chi(t) for t in _mids(dt, n_steps)]
        rel, tot = _census_dt(mult, dt)
        _say(f"IMP7-9 W0 run {key}: {len(mult)} (half-)steps, on {int(sum(mult))}; sum chi*dt "
             f"{tot:.9f} s vs {ON_TIME_END:.9f} s rel {rel:.3e} (bar {S0_BAND:.0e}); (Grid) "
             f"{grid:.3e} (bar {GRID_BAR:.0e}); period {PERIOD / dt:.9f} steps")
        assert grid <= GRID_BAR, (key, grid)
        assert rel <= S0_BAND, (key, rel)
    mult = [chi(t) for t in _mids(DT, N_RUN)]
    assert _census_dt(mult, DT) == _census(mult)


def test_preconditions(run):
    assert run["E"] <= E_BAND, run["E"]
    assert run["S0"] <= S0_BAND, run["S0"]
    assert run["grid"] <= GRID_BAR, run["grid"]
    assert bool(np.all(run["ok"])), run["ok"]
    assert np.all(np.isfinite(run["th"]))
    assert run["th"].shape == T_GRID.shape + (PTS.shape[0],)
    for scale, lab in SCALES:
        mx = T_INF + np.max(run["th"], axis=0) * scale
        _say(f"IMP7-9 max T {lab}: {np.array2string(mx, precision=6)} C (< {T_REF} ASSERTED)")
        assert np.all(mx < T_REF), (lab, mx)


def test_t_at_end_and_continuous_control(run):
    e, c = run["T"], run["cC"]
    _say(f"IMP7-9 gates run {run['key']}: (T) {e:.4e} <= {INTERIOR_BAND:.0e} "
         f"{'PASS' if e <= INTERIOR_BAND else 'FAIL'}; (c-cont) {c:.4e} > {INTERIOR_BAND:.0e} "
         f"{'PASS' if c > INTERIOR_BAND else 'FAIL'}")
    assert e <= INTERIOR_BAND, e
    assert c > INTERIOR_BAND, c


def test_printed_gap_table(run):
    """The record (no band, this item sets none); asserts only that it was formed."""
    rows = _gaps(run["th"], run["th_s"])
    _say("IMP7-9 GAP point | run | peak-scale dose rel | equal-mean dose rel | six-minute peak "
         "delta (K) | max_n|dtheta| (K)")
    for name in LABELS:
        for tag, r in (("item 700", ITEM_700[name]), (run["key"], rows[name])):
            _say(f"IMP7-9 GAP {name} | {tag} | {r[0]:.3e} | {r[1]:.3e} | {r[2]:+.3e} | {r[3]:.3e}")
        assert np.all(np.isfinite(rows[name]))
