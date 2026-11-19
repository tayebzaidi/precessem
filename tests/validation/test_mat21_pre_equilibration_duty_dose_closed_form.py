"""`MAT-21` step 32 — the wait in a dose under D = 0.25, closed form (§9 item 753, re-posing item
746, §10 chain P, P59, provisional under rule (6); brief `docs/planning/chunks/MAT-21.md` "§9
item 753", Review ruling 2026-11-02 03:00).

Numpy.  Sets no band and moves none.  θ(t; t_w, D) = ``_uniform(p, T_SW25)``, the uniform-source
duty sum Σ (−1)^m θ(t − t_m), plus the drift at t + t_w, at P5, ``ADJ_PT`` and P0 on the box's
``T_GRID``, for t_w = 0 and ``T_C``.  ``T_INF`` (MAT-15's) and ``T_INF_0`` (the box's) are
asserted equal; the histories add to ``T_INF``.

Asserted: (Q-d) ``cem43`` against its quadrature at ``DOSE_BAR`` and (Win-d) the running average
at ``AVG_BAR``, per point and wait; (Sym-d) at P0 the dose equals the drift-free duty dose at
``DOSE_BAR``; (D50-cf) at each point the dose of ``_uniform(p, T_SW)`` against the dose of
``T_INF`` + ``_duty_cf(*p, t)`` at ``DOSE_BAR``.  Negative controls: D = 0.5 misses D = 0.25 by
> ``DOSE_BAR``; the R swap misses by > ``DOSE_BAR``; ``_duty_cf(..., drop=(0,))`` misses (D50-cf)
by > ``DOSE_BAR``.  Printed (**no band**): D / D_ref − 1 at P5 and ``ADJ_PT`` for D = 1, 0.5,
0.25, peak and equal-mean, on the quadratic ``_duty_hist`` and on the uniform sum; the one-switch
``_uniform`` against ``_combined``, true by construction, so it closes nothing (rule (h)).

Not claimed: a recommended wait, a clinical duty.
"""

from __future__ import annotations

import numpy as np

from tests.validation.test_mat21_box_closed_form import T_C
from tests.validation.test_mat21_box_closed_form import T_INF as T_INF_0
from tests.validation.test_mat21_box_dose_closed_form import AVG_BAR, DOSE_BAR, T_GRID, W_B2
from tests.validation.test_mat21_box_dt_duty125_dose_closed_form import _q_check_g, _win_check_g
from tests.validation.test_mat21_box_duty25_dose_closed_form import DUTY25, T_SW25
from tests.validation.test_mat21_box_duty_closed_form import DT_BOX, DUTY, T_SW, _duty_cf
from tests.validation.test_mat21_box_duty_dose_closed_form import _duty_hist
from tests.validation.test_mat21_box_probe_points_fem import PROBE_LABELS, PROBE_PTS
from tests.validation.test_mat21_initial_gradient_closed_form import T_STAR
from tests.validation.test_mat21_pre_equilibration_closed_form import (
    ADJ_PT,
    P5,
    T_INF,
    _combined,
    _drift,
    _sourced,
)

P0 = PROBE_PTS[0]
PTS = (P5, ADJ_PT, P0)
LABS = (PROBE_LABELS[5], "ADJ_PT", PROBE_LABELS[0])
WAITS = (("0", 0.0), ("T_C", T_C))
SCHEDS = (("1", 1.0, (0.0,)), ("0.5", DUTY, T_SW), ("0.25", DUTY25, T_SW25))
_CACHE = {}


def _say(msg):
    print(msg, flush=True)


def _dose(th):
    return _q_check_g(T_INF + th, T_GRID)[0]


def _drifts(p, t_w):
    return np.array([_drift(p, t + t_w) for t in T_GRID])


def _uniform(p, sw):
    """The uniform-source duty sum, Σ_{t > t_m} (−1)^m θ(t − t_m); row 0 = 0.0."""
    return np.array(
        [
            sum((-1.0) ** m * _sourced(p, t - tm) for m, tm in enumerate(sw) if t > tm)
            for t in T_GRID
        ],
        dtype=float,
    )


def _hists():
    if not _CACHE:
        _CACHE["quad"] = {lab: _duty_hist(PTS, switch_times=sw) for lab, _d, sw in SCHEDS}
        _CACHE["uniform"] = {
            lab: np.column_stack([_uniform(p, sw) for p in PTS]) for lab, _d, sw in SCHEDS
        }
        _CACHE["drift"] = {(j, wl): _drifts(p, w) for j, p in enumerate(PTS) for wl, w in WAITS}
    return _CACHE["quad"], _CACHE["uniform"], _CACHE["drift"]


def test_preconditions():
    assert T_INF == T_INF_0, (T_INF, T_INF_0)
    assert abs(T_GRID[-1] / T_STAR - 1.0) <= 1e-14 and W_B2 <= T_STAR
    assert all(min(abs(T_GRID - t)) <= 1e-9 * DT_BOX for t in T_SW25)
    assert all(min(abs(T_GRID - t)) <= 1e-9 * DT_BOX for t in T_SW)
    _say(
        f"MAT21-32 T_INF {T_INF} == T_INF_0 {T_INF_0}; {T_GRID.size} rows, DT {DT_BOX:.4f} s, "
        f"T_C {T_C:.1f} s; {len(T_SW25)} + {len(T_SW)} switches on the grid; bars DOSE "
        f"{DOSE_BAR:g}, AVG {AVG_BAR:g}"
    )


def test_q_d_win_d_sym_d_and_rc():
    _quad_h, uni, drift = _hists()
    d25 = uni["0.25"]
    d50 = uni["0.5"]
    worst = {"q": 0.0, "w": 0.0, "rc": np.inf, "sym": 0.0, "d50": np.inf}
    for j, lab in enumerate(LABS):
        d_free = _dose(d25[:, j])
        for wl, _w in WAITS:
            th = d25[:, j] + drift[(j, wl)]
            dose, _quad, rel, _swap, miss = _q_check_g(T_INF + th, T_GRID)
            mx, t_mx, lo, hi, _ref, wrel = _win_check_g(th, T_GRID)
            assert rel <= DOSE_BAR, (lab, wl, rel)
            assert miss > DOSE_BAR, (lab, wl, miss)
            assert lo <= mx <= hi and T_GRID[0] + W_B2 <= t_mx <= T_GRID[-1], (lab, wl, t_mx)
            assert wrel <= AVG_BAR, (lab, wl, wrel)
            worst["q"] = max(worst["q"], rel)
            worst["w"] = max(worst["w"], wrel)
            worst["rc"] = min(worst["rc"], miss)
            d5 = abs(_dose(d50[:, j] + drift[(j, wl)]) / dose - 1.0)
            worst["d50"] = min(worst["d50"], d5)
            assert d5 > DOSE_BAR, (lab, wl, d5)
            if j == 2:
                sym = abs(dose / d_free - 1.0)
                worst["sym"] = max(worst["sym"], sym)
                assert sym <= DOSE_BAR, (wl, sym)
    _say(
        f"MAT21-32 (Q-d) worst rel {worst['q']:.3e} (bar {DOSE_BAR:g}, ASSERTED); (Win-d) worst "
        f"rel {worst['w']:.3e} (bar {AVG_BAR:g}, ASSERTED); (Sym-d) P0 worst rel {worst['sym']:.3e} "
        f"(ASSERTED); R swap least miss {worst['rc']:.4e}, D = 0.5 least miss {worst['d50']:.4e} "
        f"(asserted > {DOSE_BAR:g})"
    )


def test_table_printed():
    quad, uni, drift = _hists()
    _say(
        "MAT21-32 TABLE source | D | point | t_w | D_ref peak s | peak D/D_ref-1 | equal-mean "
        "D/D_ref-1"
    )
    for src in ("quadratic", "uniform"):
        for slab, duty, _sw in SCHEDS:
            for j, lab in enumerate(LABS[:2]):
                hist = quad[slab][:, j] if src == "quadratic" else uni[slab][:, j]
                for wl, _w in WAITS:
                    dr = drift[(j, wl)]
                    d_pk = _dose(hist)
                    d_eq = _dose(hist / duty)
                    _say(
                        f"MAT21-32 TABLE {src} | {slab} | {lab} | {wl} | {d_pk:.6e} | "
                        f"{_dose(hist + dr) / d_pk - 1.0:+.6e} | "
                        f"{_dose(hist / duty + dr) / d_eq - 1.0:+.6e}"
                    )


def test_d50_cf_against_the_box_duty_series():
    _quad_h, uni, drift = _hists()
    worst = 0.0
    least = np.inf
    for j, lab in enumerate(LABS):
        p = PTS[j]
        d_sum = _dose(uni["0.5"][:, j])
        d_cf = _dose(np.array([float(_duty_cf(*p, t)) for t in T_GRID]))
        d_drop = _dose(np.array([float(_duty_cf(*p, t, drop=(0,))) for t in T_GRID]))
        rel = abs(d_sum / d_cf - 1.0)
        miss = abs(d_drop / d_cf - 1.0)
        worst = max(worst, rel)
        least = min(least, miss)
        _say(
            f"MAT21-32 (D50-cf) {lab}: sum dose {d_sum:.6e} s, _duty_cf dose {d_cf:.6e} s, rel "
            f"{rel:.3e} (bar {DOSE_BAR:g}, ASSERTED); drop=(0,) miss {miss:.4e} (asserted > "
            f"{DOSE_BAR:g})"
        )
        assert rel <= DOSE_BAR, (lab, rel)
        assert miss > DOSE_BAR, (lab, miss)
    _say(f"MAT21-32 (D50-cf) worst {worst:.3e}; drop=(0,) least miss {least:.4e}")
    for j, lab in enumerate(LABS):
        for wl, w in WAITS:
            ref = _dose(np.array([_combined(PTS[j], t, w) for t in T_GRID]))
            rel_u = abs(_dose(uni["1"][:, j] + drift[(j, wl)]) / ref - 1.0)
            _say(
                f"MAT21-32 one-switch sum vs _combined {lab} t_w {wl}: rel {rel_u:.4e} (printed, "
                "by construction, closes nothing)"
            )
