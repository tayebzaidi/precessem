"""`IMP-7` step 13 — the FEM dose around the curved adiabatic hole at D = 1, the third M_p
point (§9 item 726, §10 chain I, I52, provisional under rule (6); brief
`docs/planning/chunks/IMP-7.md` "§9 item 726").

Sets no band and moves none; every band is imported.  Item 719's module body at D = 1 with
no ``source_schedule`` (the constant source), at ``LADDER[0]``: CN with two BE start-up
half-steps, ``DT``, 160 steps to 4 t\\*, 161 reads on I49's ``P_POLE`` and ``P_MID`` in one
fixed batch.  Real build.

Asserted: (E) ≤ ``E_BAND``, every point valid, max T < 43 °C; (T) interior rel L2 at 4 t\\*
against the continuous series ≤ ``INTERIOR_BAND``; (Q-hd) / (Win-hd) on the FEM history on
the peak scale at ``DOSE_BAR`` / ``AVG_BAR``, at pole and mid (equal mean is the same scale
at D = 1, and ``SCALES``' factor is bound to D = 0.5, so it is not read); (C100) the D = 1
history equals the continuous series at ``B1_BAND``.  (S0) has no content at D = 1: the
constant source records no multipliers, and ``_chi(PERIOD, 1.0)`` reads 1 at every
midpoint.  Negative controls: the D = 0.5 series at ``T_END`` misses > ``INTERIOR_BAND``;
the R swap misses > ``DOSE_BAR``.  Printed (records, **no band**): the M_p table at the pole
at D = 1, 0.5, 0.25, against the prediction M_p ≈ 1.28e-03 K at t ≈ 7.84 s.

Not claimed: a dose accuracy, a rate, any change to (Bd), which stays a reading.
"""

from __future__ import annotations

import time

import numpy as np
import pytest
from mpi4py import MPI

from fem_em_solver.thermal.dose import cem43, running_average
from tests.validation.test_imp4_pec_cylinder_closed_form import B1_BAND
from tests.validation.test_imp7_duty25_dose_fem import (
    CAP_50_POLE,
    ITEM_700_POLE,
    _duty_fn_at,
    _mids,
)
from tests.validation.test_imp7_duty_cycle_fem import PERIOD, S0_BAND, A, B, _chi
from tests.validation.test_imp7_duty_dose_closed_form import (
    AVG_BAR,
    DOSE_BAR,
    DT,
    T_END,
    T_GRID,
)
from tests.validation.test_imp7_duty_dose_fem import N_RUN, STARTUP, _err, _probe, _shell_fixture
from tests.validation.test_imp7_duty_matrix_closed_form import M_LABELS, M_PTS, _duty_hist_at
from tests.validation.test_imp7_implant_heating_closed_form import T_INF
from tests.validation.test_imp7_transient_closed_form import _modes, _theta
from tests.validation.test_imp7_transient_fem import LADDER, _series_fn, _transient
from tests.validation.test_mat12_transient_sar_heating_sphere import E_BAND
from tests.validation.test_mat13_transient_perfusion_sphere import INTERIOR_BAND
from tests.validation.test_mat20_running_average import W_B2
from tests.validation.test_mat20_tuned_heating_dose import T_REF
from tests.validation.test_mat21_wall_duty_dose_closed_form import _q_check, _win_check

D100 = 1.0
#: D = 0.25 at the pole, copied from `20261101T150144Z_IMP-7-step11-W1.log:131`, never
#: re-run: (M_p / K, t(M_p) / s, cap / K).
ITEM_719_POLE = (1.2845e-03, 7.84, 1.2410e-03)
#: The D = 1 prediction, the brief's; printed beside the reading, never asserted.
PREDICTED_M_P = 1.28e-03
PREDICTED_T_M_P = 7.84


def _say(msg):
    if MPI.COMM_WORLD.rank == 0:
        print(msg, flush=True)


def _census_100():
    """The multipliers the schedule would record at D = 1, and their on-time."""
    mult = [_chi(PERIOD, D100)(t) for t in _mids()]
    dts = [0.5 * DT] * STARTUP + [DT] * (len(mult) - STARTUP)
    tot = float(np.dot(np.asarray(mult), np.asarray(dts)))
    return mult, abs(tot - T_END) / T_END, tot


def _continuous_hist(pts, modes):
    """The continuous series at ``pts``; row 0 = 0.0 exactly, as `_duty_hist_at`.  I49's
    route: ``_series_fn`` rebuilds the modes per call, which cost W0 its 120 s."""
    out = np.zeros((T_GRID.size, pts.shape[0]))
    for n in range(1, T_GRID.size):
        out[n] = _theta(pts, T_GRID[n], modes)
    return out


@pytest.fixture(scope="module")
def fem_dose():
    t0 = time.perf_counter()
    out = {}
    res = LADDER[0]
    sub, sft, q_f, t_build, parent = _shell_fixture(res)
    _say(
        f"\nIMP7-13 res={res} sub_cells={sub.topology.index_map(3).size_global} "
        f"build={t_build:.1f}s; ranks {sub.comm.size}; dt {DT:.4f} s, {N_RUN} steps to "
        f"{T_END:.4f} s, constant source (D = {D100})"
    )
    cn, ts = _transient(
        sub,
        sft,
        q_f,
        dt=DT,
        n_steps=N_RUN,
        theta=0.5,
        startup_half_steps=STARTUP,
        snapshot_steps=tuple(range(N_RUN + 1)),
    )
    out["E"] = float(cn.max_balance_residual)
    out["mult"] = list(cn.source_multipliers)
    _say(
        f"IMP7-13 CN solve {ts:.1f}s; (E) {out['E']:.3e} (E_BAND {E_BAND:.0e}); recorded "
        f"multipliers {len(out['mult'])} (constant source, (S0) has no content)"
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
        f"IMP7-13 reads: {N_RUN + 1} batches of {M_PTS.shape[0]} in "
        f"{time.perf_counter() - tq:.1f}s; valid {ok_all.tolist()}"
    )

    tq = time.perf_counter()
    e_all, out["T"] = _err(sub, cn.snapshots[N_RUN], _series_fn(T_END, "neumann"))
    _say(
        f"IMP7-13 (T) at 4 t*: interior rel L2 {out['T']:.4e} (INTERIOR_BAND "
        f"{INTERIOR_BAND:.0e}); whole {e_all:.4e} (printed); {time.perf_counter() - tq:.1f}s"
    )
    tq = time.perf_counter()
    out["c50"] = _err(sub, cn.snapshots[N_RUN], _duty_fn_at(T_END, 0.5))[1]
    _say(
        f"IMP7-13 control at T_END: D = 0.5 series miss {out['c50']:.4e} (asserted > "
        f"{INTERIOR_BAND:.0e}); {time.perf_counter() - tq:.1f}s"
    )
    cn.snapshots.clear()

    tq = time.perf_counter()
    out["th_s"] = _duty_hist_at(_modes(), M_PTS, D100)
    _say(
        f"IMP7-13 series history {out['th_s'].shape} in {time.perf_counter() - tq:.1f}s; "
        f"fixture elapsed {time.perf_counter() - t0:.1f}s"
    )
    return out


def test_imp7_13_w0_grid_census_and_c100():
    """W0 (numpy): the imports, the grid, (S0)'s content, (C100), with no solve."""
    assert T_GRID.size == N_RUN + 1 and abs(T_GRID[-1] / T_END - 1.0) <= 1e-14
    assert A < M_PTS[0, 2] < M_PTS[1, 2] < B
    assert list(M_LABELS) == ["pole", "mid"]
    mult, rel, tot = _census_100()
    _say(
        f"IMP7-13 W0: {len(mult)} (half-)steps, all 1: {all(m == 1.0 for m in mult)}; sum "
        f"chi*dt {tot:.6f} s vs {T_END:.6f} s rel {rel:.3e} (bar {S0_BAND:.0e})"
    )
    assert all(m == 1.0 for m in mult)
    assert rel <= S0_BAND, rel

    tq = time.perf_counter()
    modes = _modes()
    d1 = _duty_hist_at(modes, M_PTS, D100)
    cont = _continuous_hist(M_PTS, modes)
    worst = float(np.max(np.abs(d1[1:] - cont[1:]) / np.abs(cont[1:])))
    _say(
        f"IMP7-13 (C100) D = 1 vs continuous over {cont.shape}: worst rel {worst:.3e} (band "
        f"{B1_BAND:g}, ASSERTED); {time.perf_counter() - tq:.1f}s"
    )
    assert worst <= B1_BAND, worst


def test_preconditions(fem_dose):
    assert fem_dose["E"] <= E_BAND, fem_dose["E"]
    assert fem_dose["mult"] == [], fem_dose["mult"]
    assert bool(np.all(fem_dose["ok"])), fem_dose["ok"]
    assert np.all(np.isfinite(fem_dose["th"]))
    mx = T_INF + np.max(fem_dose["th"], axis=0)
    _say(f"IMP7-13 max T peak scale: {np.array2string(mx, precision=6)} C (< {T_REF} ASSERTED)")
    assert np.all(mx < T_REF), mx


def test_t_at_end_and_series_control(fem_dose):
    e, c50 = fem_dose["T"], fem_dose["c50"]
    _say(
        f"IMP7-13 gates: (T) {e:.4e} <= {INTERIOR_BAND:.0e} "
        f"{'PASS' if e <= INTERIOR_BAND else 'FAIL'}; (c-50) {c50:.4e} > {INTERIOR_BAND:.0e} "
        f"{'PASS' if c50 > INTERIOR_BAND else 'FAIL'}"
    )
    assert e <= INTERIOR_BAND, e
    assert c50 > INTERIOR_BAND, c50


def test_q_hd_win_hd_and_rc_hd(fem_dose):
    th_d = fem_dose["th"]
    for j, name in enumerate(M_LABELS):
        th = th_d[:, j]
        dose, ref, rel, swap, miss = _q_check(T_GRID, T_INF + th)
        mx, t_mx, lo, hi, wref, wrel = _win_check(T_GRID, th)
        _say(
            f"IMP7-13 (Q-hd) peak {name}: FEM CEM43 {dose:.6e} s vs quad {ref:.6e} s rel "
            f"{rel:.3e} (bar {DOSE_BAR:g}, ASSERTED); (Rc-hd) swap miss {miss:.4e} (asserted "
            f"> {DOSE_BAR:g}), factor {swap / ref:.3e}"
        )
        _say(
            f"IMP7-13 (Win-hd) peak {name}: max six-minute mean {mx:.6e} K at {t_mx:.2f} s; "
            f"[{lo:.6e}, {hi:.6e}]; quad {wref:.6e} rel {wrel:.3e} (bar {AVG_BAR:g}, ASSERTED)"
        )
        assert rel <= DOSE_BAR, (name, rel)
        assert miss > DOSE_BAR, (name, miss)
        assert T_GRID[0] + W_B2 <= t_mx <= T_GRID[-1], (name, t_mx)
        assert lo <= mx <= hi, (name, lo, mx, hi)
        assert wrel <= AVG_BAR, (name, wrel)


def test_printed_gap_and_m_p_table(fem_dose):
    """The records (no band, this item sets none); asserts only that they were formed.
    An M_p far from the prediction refutes I50's reading, which is the outcome's to say."""
    th_h, th_s = fem_dose["th"], fem_dose["th_s"]
    W = float(W_B2)
    rows = {}
    for j, name in enumerate(M_LABELS):
        d_h = cem43(T_GRID, T_INF + th_h[:, j])
        d_s = cem43(T_GRID, T_INF + th_s[:, j])
        _a, mx_h, t_h = running_average(T_GRID, th_h[:, j], window=W)
        _a, mx_s, t_s = running_average(T_GRID, th_s[:, j], window=W)
        dth = np.abs(th_h[:, j] - th_s[:, j])
        M, t_M = float(np.max(dth)), float(T_GRID[int(np.argmax(dth))])
        cap = INTERIOR_BAND * float(np.max(th_s[:, j]))
        rows[name] = (M, t_M, cap)
        _say(
            f"IMP7-13 GAP D=1 {name}: peak-scale dose {d_h:.6e} s vs {d_s:.6e} s, rel "
            f"{abs(d_h / d_s - 1.0):.3e}; six-minute peak {mx_h:.6f} K at {t_h:.2f} s vs "
            f"{mx_s:.6f} K at {t_s:.2f} s, delta {mx_h - mx_s:+.3e} K"
        )
        _say(
            f"IMP7-13 (Bd-ii) D=1 {name}: M_p {M:.4e} K vs cap {cap:.4e} K, headroom cap/M_p "
            f"{cap / M:.3f} (a reading, no gate); t(M_p) {t_M:.2f} s"
        )
        assert np.isfinite(d_h / d_s)
    M, t_M, cap = rows["pole"]
    m50, m25, t25, cap25 = ITEM_700_POLE[4], *ITEM_719_POLE
    _say("IMP7-13 M_p TABLE pole (D, M_p K, t(M_p) s, cap K, headroom):")
    _say(f"IMP7-13   1.00 {M:.4e} {t_M:.2f} {cap:.4e} {cap / M:.3f} (this run)")
    _say(f"IMP7-13   0.50 {m50:.4e} -- {CAP_50_POLE:.4e} {CAP_50_POLE / m50:.3f} (item 700)")
    _say(f"IMP7-13   0.25 {m25:.4e} {t25:.2f} {cap25:.4e} {cap25 / m25:.3f} (item 719)")
    _say(
        f"IMP7-13 prediction (no gate): M_p {PREDICTED_M_P:.2e} K at t {PREDICTED_T_M_P:.2f} s; "
        f"measured / predicted M_p {M / PREDICTED_M_P:.3f}"
    )
