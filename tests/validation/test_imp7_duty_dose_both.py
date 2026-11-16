"""`IMP-7` step 12 — the fourth corner of the FEM-vs-series dose gap: ``LADDER[1]`` and half
Δt together (§9 item 725, chain I, I51, provisional; brief `docs/planning/chunks/IMP-7.md`
"§9 item 725").

Sets no band and moves none; every band is imported.  Item 710's route, imported from
``test_imp7_duty_dose_ladder.py``, untouched, with ``run``'s body written out for this one
run: Δt = ``DT`` / 2, 320 steps, read at the even steps, the same 161 rows of ``T_GRID``.

Asserted: (T) interior rel L2 at 4 t\\* ≤ ``INTERIOR_BAND``; (E) ≤ ``E_BAND``; (S0) census ≤
``S0_BAND``; (Grid) ≤ ``GRID_BAR``; every point valid; max T < 43 °C; the continuous series
at ``T_END`` misses > ``INTERIOR_BAND``.  Printed (a record, **no band**): the four-row table
at pole / equator / outer pole.  Real build.

Not claimed: a rate, a dose band, a converged dose, a metal core.
"""

from __future__ import annotations

import time

import numpy as np
import pytest
from mpi4py import MPI

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
    _err,
    _probe,
    _shell_fixture,
)
from tests.validation.test_imp7_duty_dose_ladder import (
    GRID_BAR,
    ITEM_700,
    _census_dt,
    _gaps,
    _mids,
    _read_times,
)
from tests.validation.test_imp7_implant_heating_closed_form import T_INF
from tests.validation.test_imp7_transient_closed_form import _modes
from tests.validation.test_imp7_transient_fem import LADDER, _series_fn, _transient
from tests.validation.test_mat12_transient_sar_heating_sphere import E_BAND
from tests.validation.test_mat13_transient_perfusion_sphere import INTERIOR_BAND
from tests.validation.test_mat20_tuned_heating_dose import T_REF

RUN = (LADDER[1], 0.5 * DT, 2 * N_RUN, 2)
#: Rows 2 and 3, copied from `20261101T120136Z_IMP-7-step9-W1a.log:111–115` and
#: `20261101T120345Z_IMP-7-step9-W1b.log:100–104`, never re-run: as `ITEM_700`.
ITEM_710_DT = {"pole": (9.462e-05, 1.903e-04, -6.774e-05, 3.413e-04),
               "equator": (1.732e-05, 3.412e-05, -7.082e-06, 3.889e-05),
               "outer pole": (2.484e-06, 5.219e-06, +4.917e-06, 1.697e-05)}
ITEM_710_H = {"pole": (1.310e-04, 2.608e-04, -8.709e-05, 1.292e-03),
              "equator": (1.317e-05, 2.686e-05, -9.836e-06, 2.321e-04),
              "outer pole": (4.920e-06, 9.417e-06, +2.462e-07, 5.992e-05)}


def _say(msg):
    if MPI.COMM_WORLD.rank == 0:
        print(msg, flush=True)


@pytest.fixture(scope="module")
def run():
    res, dt, n_steps, stride = RUN
    t0 = time.perf_counter()
    out = {}
    sub, sft, q_f, t_build, parent = _shell_fixture(res)
    _say(f"\nIMP7-12 run: res={res} sub_cells={sub.topology.index_map(3).size_global} "
         f"build={t_build:.1f}s; ranks {sub.comm.size}; dt {dt:.5f} s, {n_steps} steps to "
         f"{T_END:.4f} s, reads every {stride} step(s)")
    steps = tuple(range(0, n_steps + 1, stride))
    cn, ts = _transient(sub, sft, q_f, dt=dt, n_steps=n_steps, theta=0.5,
                        startup_half_steps=STARTUP, source_schedule=_chi(PERIOD, DUTY),
                        snapshot_steps=steps)
    out["E"] = float(cn.max_balance_residual)
    out["S0"], tot = _census_dt(cn.source_multipliers, dt)
    _say(f"IMP7-12 CN solve {ts:.1f}s; (E) {out['E']:.3e} (E_BAND {E_BAND:.0e}); (S0) sum chi*dt "
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
    _say(f"IMP7-12 reads: {len(steps)} batches of {PTS.shape[0]} in "
         f"{time.perf_counter() - tq:.1f}s; valid {ok_all.tolist()}; (Grid) {out['grid']:.3e} "
         f"(bar {GRID_BAR:.0e})")

    tq = time.perf_counter()
    e_all, out["T"] = _err(sub, cn.snapshots[n_steps], _duty_fn(T_END))
    _say(f"IMP7-12 (T) at 4 t*: interior rel L2 {out['T']:.4e} (INTERIOR_BAND "
         f"{INTERIOR_BAND:.0e}); whole {e_all:.4e} (printed); {time.perf_counter() - tq:.1f}s")
    tq = time.perf_counter()
    out["cC"] = _err(sub, cn.snapshots[n_steps], _series_fn(T_END, "neumann"))[1]
    _say(f"IMP7-12 control (c-cont) continuous series at T_END: interior miss {out['cC']:.4e} "
         f"(asserted > {INTERIOR_BAND:.0e}); {time.perf_counter() - tq:.1f}s")
    cn.snapshots.clear()

    tq = time.perf_counter()
    out["th_s"] = _duty_hist(_modes())
    _say(f"IMP7-12 series history {out['th_s'].shape} in {time.perf_counter() - tq:.1f}s; "
         f"fixture elapsed {time.perf_counter() - t0:.1f}s")
    return out


def test_imp7_12_w0_grid_and_census():
    """W0 (numpy): this run's grid and census arithmetic, with no solve."""
    _res, dt, n_steps, stride = RUN
    grid = float(np.max(np.abs(_read_times(dt, n_steps, stride)[1:] / T_GRID[1:] - 1.0)))
    chi = _chi(PERIOD, DUTY)
    mult = [chi(t) for t in _mids(dt, n_steps)]
    rel, tot = _census_dt(mult, dt)
    _say(f"IMP7-12 W0: {len(mult)} (half-)steps, on {int(sum(mult))}; sum chi*dt {tot:.9f} s vs "
         f"{ON_TIME_END:.9f} s rel {rel:.3e} (bar {S0_BAND:.0e}); (Grid) {grid:.3e} (bar "
         f"{GRID_BAR:.0e}); res {_res}")
    assert grid <= GRID_BAR, grid
    assert rel <= S0_BAND, rel
    assert _res == LADDER[1] and n_steps == 2 * N_RUN


def test_preconditions(run):
    assert run["E"] <= E_BAND, run["E"]
    assert run["S0"] <= S0_BAND, run["S0"]
    assert run["grid"] <= GRID_BAR, run["grid"]
    assert bool(np.all(run["ok"])), run["ok"]
    assert np.all(np.isfinite(run["th"]))
    assert run["th"].shape == T_GRID.shape + (PTS.shape[0],)
    for scale, lab in SCALES:
        mx = T_INF + np.max(run["th"], axis=0) * scale
        _say(f"IMP7-12 max T {lab}: {np.array2string(mx, precision=6)} C (< {T_REF} ASSERTED)")
        assert np.all(mx < T_REF), (lab, mx)


def test_t_at_end_and_continuous_control(run):
    e, c = run["T"], run["cC"]
    _say(f"IMP7-12 gates: (T) {e:.4e} <= {INTERIOR_BAND:.0e} "
         f"{'PASS' if e <= INTERIOR_BAND else 'FAIL'}; (c-cont) {c:.4e} > {INTERIOR_BAND:.0e} "
         f"{'PASS' if c > INTERIOR_BAND else 'FAIL'}")
    assert e <= INTERIOR_BAND, e
    assert c > INTERIOR_BAND, c


def test_printed_four_corner_table(run):
    """The record (no band, this item sets none); asserts only that it was formed."""
    rows = _gaps(run["th"], run["th_s"])
    _say("IMP7-12 GAP point | run | peak-scale dose rel | equal-mean dose rel | six-minute peak "
         "delta (K) | max_n|dtheta| (K)")
    for name in LABELS:
        for tag, r in (("item 700", ITEM_700[name]), ("dt", ITEM_710_DT[name]),
                       ("h", ITEM_710_H[name]), ("both", rows[name])):
            _say(f"IMP7-12 GAP {name} | {tag} | {r[0]:.3e} | {r[1]:.3e} | {r[2]:+.3e} | {r[3]:.3e}")
        assert np.all(np.isfinite(rows[name]))
