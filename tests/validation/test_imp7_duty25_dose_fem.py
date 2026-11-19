"""`IMP-7` step 11 — the FEM dose around the curved adiabatic hole under D = 0.25 (§9 item
719, §10 chain I, I50; brief `docs/planning/chunks/IMP-7.md` "§9 item 719").

Sets no band and moves none; every band is imported.  Item 700's fixture and read route,
imported from ``test_imp7_duty_dose_fem.py``, untouched, at ``LADDER[0]``: CN with two BE
start-up half-steps, ``DT``, 160 steps to 4 t\\*, the schedule ``_chi(PERIOD, 0.25)`` (an
on-phase of 5 ``DT``), 161 reads on I49's ``P_POLE`` and ``P_MID`` in one fixed batch.
``_duty_fn`` and ``ON_TIME_END`` are bound to D = 0.5, so ``_duty_fn_at`` and
``ON_END_25`` carry the duty here.  Real build.

Asserted: (E) ≤ ``E_BAND``, (S0) census ≤ ``S0_BAND``, every point valid, max T < 43 °C;
(T) interior rel L2 at 4 t\\* against the D = 0.25 series ≤ ``INTERIOR_BAND``; (Q-hd) /
(Win-hd) on the FEM history at peak and equal mean power at ``DOSE_BAR`` / ``AVG_BAR``.
Negative controls: the continuous and the D = 0.5 series at ``T_END`` miss >
``INTERIOR_BAND``; the R swap misses > ``DOSE_BAR``.  Printed (records, **no band**): the
dose gap beside item 700's, the six-minute peak of each, (Bd)(ii)'s M_p against the cap.

Not claimed: a rate, a dose band, convergence, a metal core, a thermal-damage statement.
"""

from __future__ import annotations

import time

import numpy as np
import pytest
from mpi4py import MPI

from fem_em_solver.thermal.dose import cem43, running_average
from fem_em_solver.utils.analytical import pec_sphere_implant_duty_cycle_from_modes
from tests.validation.test_imp7_duty_cycle_fem import PERIOD, S0_BAND, A, B, _chi
from tests.validation.test_imp7_duty_dose_closed_form import (
    AVG_BAR,
    DOSE_BAR,
    DT,
    T_END,
    T_GRID,
)
from tests.validation.test_imp7_duty_dose_fem import (
    N_RUN,
    ON_TIME_END,
    STARTUP,
    _err,
    _probe,
    _shell_fixture,
)
from tests.validation.test_imp7_duty_matrix_closed_form import (
    M_LABELS,
    M_PTS,
    _duty_hist_at,
    _scales,
)
from tests.validation.test_imp7_implant_heating_closed_form import T_INF
from tests.validation.test_imp7_transient_closed_form import _modes
from tests.validation.test_imp7_transient_fem import INTERP_CHUNK, LADDER, _series_fn, _transient
from tests.validation.test_mat12_transient_sar_heating_sphere import E_BAND, RHO_C
from tests.validation.test_mat13_transient_perfusion_sphere import INTERIOR_BAND
from tests.validation.test_mat20_running_average import W_B2
from tests.validation.test_mat20_tuned_heating_dose import T_REF
from tests.validation.test_mat21_wall_duty_dose_closed_form import _q_check, _win_check

D25 = 0.25
ON_END_25 = ON_TIME_END * (D25 / 0.5)  # the brief's census target
#: D = 0.5 at the pole, copied from `20261031T020326Z_IMP-7-step8-W1-r2.log:123`, never
#: re-run: (peak-scale dose rel, equal-mean dose rel, six-minute peak FEM / K, six-minute
#: peak series / K, max_n|Δθ| / K).  Item 700 read no mid point, so mid has no row.
ITEM_700_POLE = (1.245e-04, 2.481e-04, 0.179697, 0.179779, 1.285e-03)
#: The D = 0.5 cap at the pole, `INTERIOR_BAND`·max θ_s: I49's printed table.
CAP_50_POLE = 2.1490e-03


def _say(msg):
    if MPI.COMM_WORLD.rank == 0:
        print(msg, flush=True)


def _mids():
    """The schedule's evaluation times (`transient.py:352–357`)."""
    return [(i + 0.5) * (0.5 * DT) for i in range(STARTUP)] + [
        (STARTUP // 2 + i + 0.5) * DT for i in range(N_RUN - STARTUP // 2)
    ]


def _census_25(mult):
    dts = [0.5 * DT] * STARTUP + [DT] * (len(mult) - STARTUP)
    tot = float(np.dot(np.asarray(mult), np.asarray(dts)))
    return abs(tot - ON_END_25) / ON_END_25, tot


def _duty_fn_at(t, duty):
    """`_duty_fn` (bound to ``DUTY``) at any duty."""
    modes = _modes()

    def ev(x):
        n = x.shape[1]
        out = np.empty(n)
        for s in range(0, n, INTERP_CHUNK):
            out[s : s + INTERP_CHUNK] = (
                pec_sphere_implant_duty_cycle_from_modes(
                    x[:, s : s + INTERP_CHUNK].T,
                    t,
                    modes,
                    rho_c=RHO_C,
                    t_inf=T_INF,
                    period=PERIOD,
                    duty=duty,
                )
                - T_INF
            )
        return out

    return ev


@pytest.fixture(scope="module")
def fem_dose():
    t0 = time.perf_counter()
    out = {}
    res = LADDER[0]
    sub, sft, q_f, t_build, parent = _shell_fixture(res)
    _say(
        f"\nIMP7-11 res={res} sub_cells={sub.topology.index_map(3).size_global} "
        f"build={t_build:.1f}s; ranks {sub.comm.size}; dt {DT:.4f} s, {N_RUN} steps to "
        f"{T_END:.4f} s, period {PERIOD:.4f} s, duty {D25}"
    )
    cn, ts = _transient(
        sub,
        sft,
        q_f,
        dt=DT,
        n_steps=N_RUN,
        theta=0.5,
        startup_half_steps=STARTUP,
        source_schedule=_chi(PERIOD, D25),
        snapshot_steps=tuple(range(N_RUN + 1)),
    )
    out["E"] = float(cn.max_balance_residual)
    out["S0"], tot = _census_25(cn.source_multipliers)
    _say(
        f"IMP7-11 CN solve {ts:.1f}s; (E) {out['E']:.3e} (E_BAND {E_BAND:.0e}); (S0) sum chi*dt "
        f"{tot:.9f} s vs {ON_END_25:.9f} s rel {out['S0']:.3e} (bar {S0_BAND:.0e}); "
        f"{len(cn.source_multipliers)} (half-)steps; first two {list(cn.source_multipliers[:2])}"
    )

    tq = time.perf_counter()
    th_h = np.empty((N_RUN + 1, M_PTS.shape[0]))
    ok_all = np.ones(M_PTS.shape[0], dtype=bool)
    for k in range(N_RUN + 1):
        th, ok = _probe(cn.snapshots[k], M_PTS, parent)
        th_h[k] = th
        ok_all &= ok
    out["th"], out["ok"] = th_h, ok_all
    _say(
        f"IMP7-11 reads: {N_RUN + 1} batches of {M_PTS.shape[0]} in "
        f"{time.perf_counter() - tq:.1f}s; valid {ok_all.tolist()}"
    )

    tq = time.perf_counter()
    e_all, out["T"] = _err(sub, cn.snapshots[N_RUN], _duty_fn_at(T_END, D25))
    _say(
        f"IMP7-11 (T) at 4 t*: interior rel L2 {out['T']:.4e} (INTERIOR_BAND "
        f"{INTERIOR_BAND:.0e}); whole {e_all:.4e} (printed); {time.perf_counter() - tq:.1f}s"
    )
    tq = time.perf_counter()
    out["cC"] = _err(sub, cn.snapshots[N_RUN], _series_fn(T_END, "neumann"))[1]
    out["c50"] = _err(sub, cn.snapshots[N_RUN], _duty_fn_at(T_END, 0.5))[1]
    _say(
        f"IMP7-11 controls at T_END: continuous series miss {out['cC']:.4e}, D = 0.5 series miss "
        f"{out['c50']:.4e} (each asserted > {INTERIOR_BAND:.0e}); "
        f"{time.perf_counter() - tq:.1f}s"
    )
    cn.snapshots.clear()

    tq = time.perf_counter()
    out["th_s"] = _duty_hist_at(_modes(), M_PTS, D25)
    _say(
        f"IMP7-11 series history {out['th_s'].shape} in {time.perf_counter() - tq:.1f}s; "
        f"fixture elapsed {time.perf_counter() - t0:.1f}s"
    )
    return out


def test_imp7_11_w0_grid_and_census():
    """W0 (numpy): the imports, the grid, the census arithmetic with no solve."""
    assert T_GRID.size == N_RUN + 1 and abs(T_GRID[-1] / T_END - 1.0) <= 1e-14
    assert abs(D25 * PERIOD / DT - 5.0) < 1e-12
    assert A < M_PTS[0, 2] < M_PTS[1, 2] < B
    assert list(M_LABELS) == ["pole", "mid"]
    chi = _chi(PERIOD, D25)
    mult = [chi(t) for t in _mids()]
    rel, tot = _census_25(mult)
    _say(
        f"IMP7-11 W0: {len(mult)} (half-)steps, on {int(sum(mult))}; first two {mult[:2]}; sum "
        f"chi*dt {tot:.9f} s vs {ON_END_25:.9f} s rel {rel:.3e} (bar {S0_BAND:.0e})"
    )
    assert mult[:2] == [1.0, 1.0]
    assert rel <= S0_BAND, rel


def test_preconditions(fem_dose):
    assert fem_dose["E"] <= E_BAND, fem_dose["E"]
    assert fem_dose["S0"] <= S0_BAND, fem_dose["S0"]
    assert bool(np.all(fem_dose["ok"])), fem_dose["ok"]
    assert np.all(np.isfinite(fem_dose["th"]))
    for scale, lab in _scales(D25):
        mx = T_INF + np.max(fem_dose["th"], axis=0) * scale
        _say(f"IMP7-11 max T {lab}: {np.array2string(mx, precision=6)} C (< {T_REF} ASSERTED)")
        assert np.all(mx < T_REF), (lab, mx)


def test_t_at_end_and_series_controls(fem_dose):
    e, c, c50 = fem_dose["T"], fem_dose["cC"], fem_dose["c50"]
    _say(
        f"IMP7-11 gates: (T) {e:.4e} <= {INTERIOR_BAND:.0e} "
        f"{'PASS' if e <= INTERIOR_BAND else 'FAIL'}; (c-cont) {c:.4e}, (c-50) {c50:.4e}, each > "
        f"{INTERIOR_BAND:.0e} {'PASS' if min(c, c50) > INTERIOR_BAND else 'FAIL'}"
    )
    assert e <= INTERIOR_BAND, e
    assert c > INTERIOR_BAND, c
    assert c50 > INTERIOR_BAND, c50


def test_q_hd_win_hd_and_rc_hd(fem_dose):
    th_d = fem_dose["th"]
    for scale, lab in _scales(D25):
        for j, name in enumerate(M_LABELS):
            th = th_d[:, j] * scale
            dose, ref, rel, swap, miss = _q_check(T_GRID, T_INF + th)
            mx, t_mx, lo, hi, wref, wrel = _win_check(T_GRID, th)
            _say(
                f"IMP7-11 (Q-hd) {lab} {name}: FEM CEM43 {dose:.6e} s vs quad {ref:.6e} s rel "
                f"{rel:.3e} (bar {DOSE_BAR:g}, ASSERTED); (Rc-hd) swap miss {miss:.4e} (asserted "
                f"> {DOSE_BAR:g}), factor {swap / ref:.3e}"
            )
            _say(
                f"IMP7-11 (Win-hd) {lab} {name}: max six-minute mean {mx:.6e} K at {t_mx:.2f} s; "
                f"[{lo:.6e}, {hi:.6e}]; quad {wref:.6e} rel {wrel:.3e} (bar {AVG_BAR:g}, ASSERTED)"
            )
            assert rel <= DOSE_BAR, (lab, name, rel)
            assert miss > DOSE_BAR, (lab, name, miss)
            assert T_GRID[0] + W_B2 <= t_mx <= T_GRID[-1], (lab, name, t_mx)
            assert lo <= mx <= hi, (lab, name, lo, mx, hi)
            assert wrel <= AVG_BAR, (lab, name, wrel)


def test_printed_gap_and_headroom_records(fem_dose):
    """The records (no band, this item sets none); asserts only that they were formed.
    M_p above the cap is a reading for the weekly, not a red (brief)."""
    th_h, th_s = fem_dose["th"], fem_dose["th_s"]
    W = float(W_B2)
    p = ITEM_700_POLE
    _say(
        f"IMP7-11 GAP item 700 D=0.5 pole: peak-scale dose rel {p[0]:.3e}; equal-mean rel "
        f"{p[1]:.3e}; six-minute peak {p[2]:.6f} K vs {p[3]:.6f} K; max_n|dtheta| {p[4]:.3e} K; "
        f"cap {CAP_50_POLE:.4e} K; headroom cap/M_p {CAP_50_POLE / p[4]:.3f}"
    )
    for j, name in enumerate(M_LABELS):
        d_h = cem43(T_GRID, T_INF + th_h[:, j])
        d_s = cem43(T_GRID, T_INF + th_s[:, j])
        de_h = cem43(T_GRID, T_INF + th_h[:, j] / D25)
        de_s = cem43(T_GRID, T_INF + th_s[:, j] / D25)
        _a, mx_h, t_h = running_average(T_GRID, th_h[:, j], window=W)
        _a, mx_s, t_s = running_average(T_GRID, th_s[:, j], window=W)
        M = float(np.max(np.abs(th_h[:, j] - th_s[:, j])))
        cap = INTERIOR_BAND * float(np.max(th_s[:, j]))
        _say(
            f"IMP7-11 GAP D=0.25 {name}: peak-scale dose {d_h:.6e} s vs {d_s:.6e} s, rel "
            f"{abs(d_h / d_s - 1.0):.3e}; equal-mean {de_h:.6e} s vs {de_s:.6e} s, rel "
            f"{abs(de_h / de_s - 1.0):.3e}; six-minute peak {mx_h:.6f} K at {t_h:.2f} s vs "
            f"{mx_s:.6f} K at {t_s:.2f} s, delta {mx_h - mx_s:+.3e} K; max_n|dtheta| {M:.3e} K"
        )
        _say(
            f"IMP7-11 (Bd-ii) D=0.25 {name}: M_p {M:.4e} K vs cap {cap:.4e} K, headroom cap/M_p "
            f"{cap / M:.3f} (a reading, no gate); t(M_p) "
            f"{T_GRID[int(np.argmax(np.abs(th_h[:, j] - th_s[:, j])))]:.2f} s"
        )
        assert np.isfinite(d_h / d_s)
