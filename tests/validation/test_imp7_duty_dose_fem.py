"""`IMP-7` step 8 — `MAT-20`'s dose and six-minute window on the FEM point histories
around the curved adiabatic hole, under the duty-cycled source, to 4 t\\* (§9 item
700, §10 chain I, I48b; brief `docs/planning/chunks/IMP-7.md` "§9 item 700").

Sets no band and moves none; every band is imported.  Item 675's fixture at
``LADDER[0]`` (P2 on the geometry-order-2 sub-mesh), CN with two BE start-up
half-steps, ``DT`` = t\\*/40, 160 steps, 161 snapshots, each read on I48a's ``PTS``
in one fixed batch on every rank.  The read is `_probe_vals`' owner route with
``parent=`` added: this mesh is a ``create_submesh`` mesh, which carries no original
cell index (`POST-10`), so `_shell_fixture` is `_fixture` keeping the parent map.
Real build.

Asserted: preconditions ((E) balance ≤ ``E_BAND``, (S0) census to ``T_END`` ≤
``S0_BAND``, every point valid, max T < 43 °C); (Q-hd) / (Win-hd) on the FEM history
at peak and equal mean power at ``DOSE_BAR`` / ``AVG_BAR``; (T) the interior rel L2
at ``T_END`` against I48a's series ≤ ``INTERIOR_BAND``.  Negative controls: the R
swap misses > ``DOSE_BAR``; the continuous series at ``T_END`` misses >
``INTERIOR_BAND``.  Printed (records, **no band**): the error at t\\* and 2 t\\*, the
FEM dose against I48a's series dose at each point, the six-minute peak of each.

Not claimed: dose convergence in h or Δt (I47's interior error does not shrink with
h at fixed Δt), a rate, a metal core, a coil source, a thermal-damage statement.
"""

from __future__ import annotations

import time

import numpy as np
import pytest
from dolfinx import mesh as dmesh
from mpi4py import MPI

from fem_em_solver.post.evaluation import evaluate_vector_field_parallel
from fem_em_solver.thermal.dose import cem43, running_average
from tests.validation import test_imp7_implant_heating_fem as m7
from tests.validation.test_imp7_curved_p2 import _ball_curved, _cmap_degree
from tests.validation.test_imp7_duty_cycle_fem import (
    DUTY,
    ON_TIME,
    PERIOD,
    S0_BAND,
    A,
    B,
    _chi,
    _duty_fn,
)
from tests.validation.test_imp7_duty_dose_closed_form import (
    AVG_BAR,
    DOSE_BAR,
    DT,
    LABELS,
    PTS,
    SCALES,
    T_END,
    T_GRID,
    _duty_hist,
)
from tests.validation.test_imp7_implant_heating_closed_form import T_INF
from tests.validation.test_imp7_transient_closed_form import T_STAR, _modes
from tests.validation.test_imp7_transient_fem import (
    LADDER,
    RT,
    _series_fn,
    _transient,
)
from tests.validation.test_mat12_transient_sar_heating_sphere import E_BAND, N_MAIN
from tests.validation.test_mat13_transient_perfusion_sphere import INTERIOR_BAND
from tests.validation.test_mat20_running_average import W_B2
from tests.validation.test_mat20_tuned_heating_dose import T_REF
from tests.validation.test_mat21_wall_duty_dose_closed_form import _q_check, _win_check

N_RUN = 4 * N_MAIN
ON_TIME_END = 4 * ON_TIME  # eight periods
STARTUP = 2


def _say(msg):
    if MPI.COMM_WORLD.rank == 0:
        print(msg, flush=True)


def _shell_fixture(res):
    """`_fixture` (`test_imp7_transient_fem.py:95`) on `_shell_submesh`'s cells, keeping
    the parent map `create_submesh` returns, which `_shell_submesh` drops."""
    msh, ct, _, t_build = _ball_curved(res)
    tdim = msh.topology.dim
    imap = msh.topology.index_map(tdim)
    vals = np.full(imap.size_local + imap.num_ghosts, -1, dtype=np.int32)
    vals[ct.indices] = ct.values
    keep = np.flatnonzero(vals == 2).astype(np.int32)
    sub, emap, _v, _g = dmesh.create_submesh(msh, tdim, keep)
    assert _cmap_degree(sub) == 2
    sim = sub.topology.index_map(tdim)
    n_sub = sim.size_local + sim.num_ghosts
    parent_cells = emap.sub_topology_to_topology(np.arange(n_sub, dtype=np.int32), False)
    sft = m7._robin_outer(sub, A, B)
    q_f = m7._source_dg2(sub, None, m7._q_shell, dg=RT["q_dg"])
    return sub, sft, q_f, t_build, (msh, parent_cells)


def _probe(T, pts, parent):
    """`_probe_vals` (`test_mat21_box_probe_points_fem.py:106`) with ``parent=``."""
    vals, ok = evaluate_vector_field_parallel(T, pts, owner="lowest_original_cell", parent=parent)
    ok = np.asarray(ok, dtype=bool)
    th = np.array(
        [float(np.real(np.ravel(vals[i])[0])) - T_INF if ok[i] else np.nan for i in range(len(pts))]
    )
    return th, ok


def _census(mult):
    dts = [0.5 * DT] * STARTUP + [DT] * (len(mult) - STARTUP)
    tot = float(np.dot(np.asarray(mult), np.asarray(dts)))
    return abs(tot - ON_TIME_END) / ON_TIME_END, tot


def _err(sub, T_h, fn):
    th = m7._p2(sub, fn, p=RT["comparand"])
    nan = sub.comm.allreduce(int(np.count_nonzero(~np.isfinite(th.x.array))), op=MPI.SUM)
    assert nan == 0, f"comparand holds {nan} NaN dofs"
    return m7._errors(sub, None, T_h, th, B, qmeta=RT["qmeta"])[:2]


@pytest.fixture(scope="module")
def fem_dose():
    t0 = time.perf_counter()
    out = {}
    res = LADDER[0]
    sub, sft, q_f, t_build, parent = _shell_fixture(res)
    _say(
        f"\nIMP7-8 res={res} sub_cells={sub.topology.index_map(3).size_global} "
        f"build={t_build:.1f}s; ranks {sub.comm.size}; dt {DT:.4f} s, {N_RUN} steps to "
        f"{T_END:.4f} s, period {PERIOD:.4f} s, duty {DUTY}"
    )
    cn, ts = _transient(
        sub,
        sft,
        q_f,
        dt=DT,
        n_steps=N_RUN,
        theta=0.5,
        startup_half_steps=STARTUP,
        source_schedule=_chi(PERIOD, DUTY),
        snapshot_steps=tuple(range(N_RUN + 1)),
    )
    out["E"] = float(cn.max_balance_residual)
    out["S0"], tot = _census(cn.source_multipliers)
    _say(
        f"IMP7-8 CN solve {ts:.1f}s; (E) {out['E']:.3e} (E_BAND {E_BAND:.0e}); (S0) sum chi*dt "
        f"{tot:.9f} s vs {ON_TIME_END:.9f} s rel {out['S0']:.3e} (bar {S0_BAND:.0e}); "
        f"{len(cn.source_multipliers)} (half-)steps"
    )

    tq = time.perf_counter()
    th_h = np.empty((N_RUN + 1, PTS.shape[0]))
    ok_all = np.ones(PTS.shape[0], dtype=bool)
    for k in range(N_RUN + 1):
        th, ok = _probe(cn.snapshots[k], PTS, parent)
        th_h[k] = th
        ok_all &= ok
    out["th"], out["ok"] = th_h, ok_all
    _say(
        f"IMP7-8 reads: {N_RUN + 1} batches of {PTS.shape[0]} in "
        f"{time.perf_counter() - tq:.1f}s; valid {ok_all.tolist()}"
    )

    # The error at t*, 2 t*, 4 t* (the brief's negative-result print), taken up front.
    out["T"] = {}
    for n, s in ((N_MAIN, "t*"), (2 * N_MAIN, "2 t*"), (N_RUN, "4 t*")):
        tq = time.perf_counter()
        t = DT * n
        e_all, e_in = _err(sub, cn.snapshots[n], _duty_fn(t))
        out["T"][s] = e_in
        _say(
            f"IMP7-8 (T) at {s} ({t:.4f} s): interior rel L2 {e_in:.4e} (INTERIOR_BAND "
            f"{INTERIOR_BAND:.0e}); whole {e_all:.4e} (printed); {time.perf_counter() - tq:.1f}s"
        )
    tq = time.perf_counter()
    out["cC"] = _err(sub, cn.snapshots[N_RUN], _series_fn(T_END, "neumann"))[1]
    _say(
        f"IMP7-8 control (c-cont) continuous series at T_END: interior miss {out['cC']:.4e} "
        f"(asserted > {INTERIOR_BAND:.0e}; item 675 read it at t*); "
        f"{time.perf_counter() - tq:.1f}s"
    )
    cn.snapshots.clear()

    tq = time.perf_counter()
    out["th_s"] = _duty_hist(_modes())
    _say(
        f"IMP7-8 series history {out['th_s'].shape} in {time.perf_counter() - tq:.1f}s; "
        f"fixture elapsed {time.perf_counter() - t0:.1f}s"
    )
    return out


def test_imp7_8_w0_imports_grid_census():
    """W0 (numpy): the imports, the grid, the census arithmetic with no solve."""
    assert T_GRID.size == N_RUN + 1 and abs(T_GRID[-1] / T_END - 1.0) <= 1e-14
    assert T_END == 4 * T_STAR
    chi = _chi(PERIOD, DUTY)
    mids = [(i + 0.5) * (0.5 * DT) for i in range(STARTUP)]
    mids += [(STARTUP // 2 + i + 0.5) * DT for i in range(N_RUN - STARTUP // 2)]
    mult = [chi(t) for t in mids]
    rel, tot = _census(mult)
    _say(
        f"IMP7-8 W0: {len(mult)} (half-)steps, on {int(sum(mult))}; sum chi*dt {tot:.9f} s vs "
        f"{ON_TIME_END:.9f} s rel {rel:.3e} (bar {S0_BAND:.0e})"
    )
    assert rel <= S0_BAND, rel


def test_preconditions(fem_dose):
    assert fem_dose["E"] <= E_BAND, fem_dose["E"]
    assert fem_dose["S0"] <= S0_BAND, fem_dose["S0"]
    assert bool(np.all(fem_dose["ok"])), fem_dose["ok"]
    assert np.all(np.isfinite(fem_dose["th"]))
    for scale, lab in SCALES:
        mx = T_INF + np.max(fem_dose["th"], axis=0) * scale
        _say(f"IMP7-8 max T {lab}: {np.array2string(mx, precision=6)} C (< {T_REF} ASSERTED)")
        assert np.all(mx < T_REF), (lab, mx)


def test_t_at_end_and_continuous_control(fem_dose):
    e = fem_dose["T"]["4 t*"]
    _say(
        f"IMP7-8 gates: (T) {e:.4e} <= {INTERIOR_BAND:.0e} {'PASS' if e <= INTERIOR_BAND else 'FAIL'}; "
        f"(c-cont) {fem_dose['cC']:.4e} > {INTERIOR_BAND:.0e} "
        f"{'PASS' if fem_dose['cC'] > INTERIOR_BAND else 'FAIL'}"
    )
    assert e <= INTERIOR_BAND, fem_dose["T"]
    assert fem_dose["cC"] > INTERIOR_BAND, fem_dose["cC"]


def test_q_hd_win_hd_and_rc_hd(fem_dose):
    th_d = fem_dose["th"]
    for scale, lab in SCALES:
        for j, name in enumerate(LABELS):
            th = th_d[:, j] * scale
            dose, ref, rel, swap, miss = _q_check(T_GRID, T_INF + th)
            mx, t_mx, lo, hi, wref, wrel = _win_check(T_GRID, th)
            _say(
                f"IMP7-8 (Q-hd) {lab} {name}: FEM CEM43 {dose:.6e} s vs quad {ref:.6e} s rel "
                f"{rel:.3e} (bar {DOSE_BAR:g}, ASSERTED); (Rc-hd) swap miss {miss:.4e} (asserted "
                f"> {DOSE_BAR:g}), factor {swap / ref:.3e}"
            )
            _say(
                f"IMP7-8 (Win-hd) {lab} {name}: max six-minute mean {mx:.6e} K at {t_mx:.2f} s; "
                f"[{lo:.6e}, {hi:.6e}]; quad {wref:.6e} rel {wrel:.3e} (bar {AVG_BAR:g}, ASSERTED)"
            )
            assert rel <= DOSE_BAR, (lab, name, rel)
            assert miss > DOSE_BAR, (lab, name, miss)
            assert T_GRID[0] + W_B2 <= t_mx <= T_GRID[-1], (lab, name, t_mx)
            assert lo <= mx <= hi, (lab, name, lo, mx, hi)
            assert wrel <= AVG_BAR, (lab, name, wrel)


def test_printed_fem_vs_series_records(fem_dose):
    """The records (no band, this item sets none); asserts only that they were formed."""
    th_h, th_s = fem_dose["th"], fem_dose["th_s"]
    W = float(W_B2)
    for j, name in enumerate(LABELS):
        d_h = cem43(T_GRID, T_INF + th_h[:, j])
        d_s = cem43(T_GRID, T_INF + th_s[:, j])
        de_h = cem43(T_GRID, T_INF + th_h[:, j] / DUTY)
        de_s = cem43(T_GRID, T_INF + th_s[:, j] / DUTY)
        _a, mx_h, t_h = running_average(T_GRID, th_h[:, j], window=W)
        _a, mx_s, t_s = running_average(T_GRID, th_s[:, j], window=W)
        dm = float(np.max(np.abs(th_h[:, j] - th_s[:, j])))
        _say(
            f"IMP7-8 FEM vs SERIES {name}: peak-scale dose {d_h:.6e} s vs {d_s:.6e} s, rel "
            f"{abs(d_h / d_s - 1.0):.3e}; equal-mean {de_h:.6e} s vs {de_s:.6e} s, rel "
            f"{abs(de_h / de_s - 1.0):.3e}; six-minute peak {mx_h:.6f} K at {t_h:.2f} s vs "
            f"{mx_s:.6f} K at {t_s:.2f} s, delta {mx_h - mx_s:+.3e} K; max_n|dtheta| {dm:.3e} K"
        )
        assert np.isfinite(d_h / d_s)
