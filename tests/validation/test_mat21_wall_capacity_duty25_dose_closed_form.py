"""`MAT-21` step 39 — the dose under the capacitive wall, D = 0.25, closed form
(§9 item 834, §10 P70; brief `docs/planning/chunks/MAT-21.md` "§9 item 834", bound by the
"Review ruling 2026-11-05 03:00" above item 823).

Numpy.  Sets no band and moves none.  P64's uniform-sum dose
(`test_mat21_wall_pre_equilibration_duty25_dose_closed_form.py`, imported unmoved) with item
823's composite (capacitive-wall) drift in place of the h_eff drift:
θ(t; t_w) = Σ_{t > t_m} (−1)^m θ_wall(t − t_m) + `_drift_cap_z`(Z, t + t_w).  Schedule
``T_SW25`` (D = 0.25) on the box's ``T_GRID``; waits {0, ``T_C``}; points P5 (z = C),
``ADJ_PT`` (0.9 C), P0 (z = 0).

Asserted: (Q-d) ``cem43`` against its quadrature at ``DOSE_BAR`` with the R-swap control and
(Win-d) the running average at ``AVG_BAR``, per point and wait, with the composite drift;
(D50-cf) the walled uniform sum at ``T_SW`` against ``two_layer_slab_duty_cycle_temperature``
(``side="left"``) at ``DOSE_BAR`` with the drop-run-0 control (drift-independent, re-asserted);
the vectorised ``_drift_cap_z`` against item 823's scalar ``_drift_cap`` (precondition).
Control (asserted, factor 1 at the gate's own band): the massless-wall (``_drift_1d``) dose
misses the composite-drift dose at P5, t_w = 0, by > ``DOSE_BAR``; the ≈ 4 % predicted from
P65's FEM tw0 against the h_eff series (+4.222e-02, `20261104T103557Z_MAT-21-P65-W1m.log:87`)
is printed beside it, never asserted.
Printed (🟡 record, no band): (Sym-d) at P0 (0 by construction); D_ref, D(0), D(T_C),
D(T_C)/D(0) and the biases at three heights, beside P64's six-face h_eff form and the z-only
``_drift_1d`` form; the composite D(0) at P5 beside P65's FEM 7.253626e-11 s; N vs 2N at
``T_GRID`` row 0.

Not claimed: a dose band, any FEM, a recommended wait, a cited PMMA ρc.
"""

from __future__ import annotations

import numpy as np

from tests.validation.test_mat21_box_closed_form import T_C
from tests.validation.test_mat21_box_dose_closed_form import AVG_BAR, DOSE_BAR, T_GRID, W_B2
from tests.validation.test_mat21_box_dt_duty125_dose_closed_form import _q_check_g, _win_check_g
from tests.validation.test_mat21_box_duty25_dose_closed_form import DUTY25, T_SW25
from tests.validation.test_mat21_box_duty_closed_form import DT_BOX, DUTY, T_SW
from tests.validation.test_mat21_initial_gradient_closed_form import N_TERMS
from tests.validation.test_mat21_pre_equilibration_duty_dose_closed_form import LABS, PTS, T_INF
from tests.validation.test_mat21_wall_capacity_drift_closed_form import _drift_cap, _drift_cap_z
from tests.validation.test_mat21_wall_pre_equilibration_closed_form import _drift_1d
from tests.validation.test_mat21_wall_pre_equilibration_duty25_dose_closed_form import (
    PERIOD_BOX,
    WAITS,
    Z,
    _dose,
    _duty_series,
    _hists,
)

# P65's FEM record (`20261104T103557Z_MAT-21-P65-W1m.log:87`): P5 tw0 FEM dose, FEM/series − 1.
P65_FEM_D0_P5 = 7.253626e-11
P65_FEM_VS_HEFF = 4.222e-02
_CAP = {}


def _say(msg):
    print(msg, flush=True)


def _cap_hists(n=N_TERMS):
    """{wait label: (41, 3) composite drift on T_GRID at Z}, vectorised (never per point)."""
    if n not in _CAP:
        _CAP[n] = {wl: np.array([_drift_cap_z(Z, t + w, n=n) for t in T_GRID])
                   for wl, w in WAITS}
    return _CAP[n]


def _z1d_hists():
    return {wl: np.array([[_drift_1d(p, t + w) for p in PTS] for t in T_GRID])
            for wl, w in WAITS}


def test_preconditions():
    assert all(min(abs(T_GRID - t)) <= 1e-9 * DT_BOX for t in T_SW25)
    assert all(min(abs(T_GRID - t)) <= 1e-9 * DT_BOX for t in T_SW)
    assert abs(PERIOD_BOX / DT_BOX - 8.0) < 1e-12 and DUTY == 0.5 and DUTY25 == 0.25
    assert [wl for wl, _w in WAITS] == ["0", "T_C"] and WAITS[1][1] == T_C
    cap = _cap_hists()
    worst = 0.0
    for wl, w in WAITS:
        for n in (0, 7, T_GRID.size - 1):
            for j, p in enumerate(PTS):
                s = _drift_cap(p, T_GRID[n] + w)
                worst = max(worst, abs(cap[wl][n, j] - s) / max(abs(s), 1e-300) if s else
                            abs(cap[wl][n, j]))
    _say(f"MAT21-39 {T_GRID.size} rows, DT {DT_BOX:.4f} s, P {PERIOD_BOX:.4f} s; "
         f"{len(T_SW25)} + {len(T_SW)} switches on the grid; vectorised _drift_cap_z vs item 823 "
         f"scalar _drift_cap max rel {worst:.3e} (asserted <= 1e-12); bars DOSE {DOSE_BAR:g}, "
         f"AVG {AVG_BAR:g}; N_TERMS {N_TERMS}")
    assert worst <= 1e-12, worst
    assert np.all(_hists()["25"][0] == 0.0)


def test_q_d_win_d_and_rc_with_composite_drift():
    h = _hists()
    cap = _cap_hists()
    d25 = h["25"]
    worst = {"q": 0.0, "w": 0.0, "rc": np.inf, "sym": 0.0}
    gates = []
    for j, lab in enumerate(LABS):
        d_free = _dose(d25[:, j])
        for wl, _w in WAITS:
            th = d25[:, j] + cap[wl][:, j]
            dose, _quad, rel, _swap, miss = _q_check_g(T_INF + th, T_GRID)
            mx, t_mx, lo, hi, _ref, wrel = _win_check_g(th, T_GRID)
            gates.append(((lab, wl, "Q"), rel <= DOSE_BAR, rel))
            gates.append(((lab, wl, "swap"), miss > DOSE_BAR, miss))
            gates.append(((lab, wl, "win-loc"),
                          lo <= mx <= hi and T_GRID[0] + W_B2 <= t_mx <= T_GRID[-1], t_mx))
            gates.append(((lab, wl, "Win"), wrel <= AVG_BAR, wrel))
            worst["q"] = max(worst["q"], rel)
            worst["w"] = max(worst["w"], wrel)
            worst["rc"] = min(worst["rc"], miss)
            if j == 2:
                worst["sym"] = max(worst["sym"], abs(dose / d_free - 1.0))
    _say(f"MAT21-39 (Q-d) worst rel {worst['q']:.3e} (bar {DOSE_BAR:g}, ASSERTED); (Win-d) worst "
         f"rel {worst['w']:.3e} (bar {AVG_BAR:g}, ASSERTED); R swap least miss "
         f"{worst['rc']:.4e} (asserted > {DOSE_BAR:g}); (Sym-d) P0 worst rel "
         f"{worst['sym']:.3e} (printed, 0 by construction)")
    bad = [g for g in gates if not g[1]]
    assert not bad, bad


def test_d50_cf_against_the_two_layer_duty_series():
    h = _hists()
    cf50 = _duty_series(DUTY)
    worst, least = 0.0, np.inf
    gates = []
    for j, lab in enumerate(LABS):
        d_sum = _dose(h["50"][:, j])
        d_cf = _dose(cf50[:, j])
        d_drop = _dose(h["50drop"][:, j])
        rel = abs(d_sum / d_cf - 1.0)
        miss = abs(d_drop / d_cf - 1.0)
        worst, least = max(worst, rel), min(least, miss)
        _say(f"MAT21-39 (D50-cf) {lab}: walled sum dose {d_sum:.6e} s, two-layer duty series dose "
             f"{d_cf:.6e} s, rel {rel:.3e} (bar {DOSE_BAR:g}, ASSERTED); drop=(0,) miss {miss:.4e} "
             f"(asserted > {DOSE_BAR:g})")
        gates.append(((lab, "rel"), rel <= DOSE_BAR, rel))
        gates.append(((lab, "drop"), miss > DOSE_BAR, miss))
    _say(f"MAT21-39 (D50-cf) worst {worst:.3e}; drop=(0,) least miss {least:.4e}")
    bad = [g for g in gates if not g[1]]
    assert not bad, bad


def test_control_massless_wall_drift_misses_composite_dose_at_p5():
    h = _hists()
    cap = _cap_hists()
    z1 = _z1d_hists()
    j = 0  # P5, z = C
    out = {}
    for wl, _w in WAITS:
        d_c = _dose(h["25"][:, j] + cap[wl][:, j])
        d_m = _dose(h["25"][:, j] + z1[wl][:, j])
        out[wl] = (d_c, d_m, abs(d_m / d_c - 1.0), d_c / d_m - 1.0)
        _say(f"MAT21-39 (Ctrl) {LABS[j]} t_w {wl}: composite-drift dose {d_c:.9e} s, massless-wall "
             f"(_drift_1d) dose {d_m:.9e} s, rel miss {out[wl][2]:.4e}, composite/massless-1 "
             f"{out[wl][3]:+.4e}" + (f" (ASSERTED > {DOSE_BAR:g}; predicted ≈ "
                                      f"{P65_FEM_VS_HEFF:+.3e}, P65-W1m.log:87, printed)"
                                      if wl == "0" else " (printed)"))
    assert out["0"][2] > DOSE_BAR, out["0"]


def test_table_printed():
    """The record (no band); asserts only that it was formed."""
    h = _hists()
    cap = _cap_hists()
    cap2 = _cap_hists(2 * N_TERMS)
    z1 = _z1d_hists()
    six = h["drift"]
    _say("MAT21-39 TABLE point | D_ref (no drift) s | D(t_w=0) s | D(t_w=T_C) s | D(T_C)/D(0) | "
         "D(0)/D_ref-1 | D(T_C)/D_ref-1 || six-face h_eff (P64) D(0)/D_ref-1 | D(T_C)/D_ref-1 || "
         "z-only _drift_1d D(0)/D_ref-1 | D(T_C)/D_ref-1")
    for j, lab in enumerate(LABS):
        ref = _dose(h["25"][:, j])
        d0 = _dose(h["25"][:, j] + cap["0"][:, j])
        dc = _dose(h["25"][:, j] + cap["T_C"][:, j])
        s0 = _dose(h["25"][:, j] + six[(j, "0")])
        sc = _dose(h["25"][:, j] + six[(j, "T_C")])
        m0 = _dose(h["25"][:, j] + z1["0"][:, j])
        mc = _dose(h["25"][:, j] + z1["T_C"][:, j])
        _say(f"MAT21-39 TABLE {lab} | {ref:.6e} | {d0:.6e} | {dc:.6e} | {dc / d0:.6f} | "
             f"{d0 / ref - 1.0:+.6e} | {dc / ref - 1.0:+.6e} || {s0 / ref - 1.0:+.6e} | "
             f"{sc / ref - 1.0:+.6e} || {m0 / ref - 1.0:+.6e} | {mc / ref - 1.0:+.6e}")
        assert ref > 0.0
    for wl, _w in WAITS:
        gap6 = max(float(np.max(np.abs(six[(j, wl)] - z1[wl][:, j]))) for j in (0, 1))
        gapc = max(float(np.max(np.abs(cap[wl][:, j] - z1[wl][:, j]))) for j in (0, 1))
        _say(f"MAT21-39 SIDE t_w {wl}: max |six-face - _drift_1d| on T_GRID at P5/ADJ_PT "
             f"{gap6:.3e} K (geometry; P65-W1m.log:84 read 1.961e-04); max |composite - "
             f"_drift_1d| {gapc:.3e} K (capacity)")
    d0 = _dose(h["25"][:, 0] + cap["0"][:, 0])
    _say(f"MAT21-39 P65 P5 t_w=0: composite D(0) {d0:.6e} s vs P65 FEM {P65_FEM_D0_P5:.6e} s "
         f"(P65-W1m.log:87), FEM/composite-1 {P65_FEM_D0_P5 / d0 - 1.0:+.4e} (printed)")
    for wl, _w in WAITS:
        r0 = np.abs(cap[wl][0] - cap2[wl][0])
        dn = max(abs(_dose(h["25"][:, j] + cap[wl][:, j])
                     / _dose(h["25"][:, j] + cap2[wl][:, j]) - 1.0) for j in range(Z.size))
        _say(f"MAT21-39 (N vs 2N) t_w {wl}: T_GRID row 0 |drift_N - drift_2N| at P5/ADJ_PT/P0 "
             f"{r0[0]:.3e} / {r0[1]:.3e} / {r0[2]:.3e} K; dose worst rel {dn:.3e} (printed)")
