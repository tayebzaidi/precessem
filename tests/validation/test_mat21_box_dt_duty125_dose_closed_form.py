"""`MAT-21` step 17, W1 — `MAT-20`'s dose and six-minute window on the gel box's
series point histories at all six F2182 points under a D = 0.125 schedule at
P38's period (P36's source), at Δt = t*/40 and t*/80 (§9 item 592, §10 chain
P45; brief `docs/planning/chunks/MAT-21.md` "§9 item 592", review ruling
2026-10-22 03:00).

Numpy, real build.  Sets no band and moves none; every band is imported.
``HIST_PTS`` and ``NEW_PTS`` are two fixed batches, never merged.  Every
inherited helper is bound to the 41-point grid; the ``*_g`` local copies take
the grid as an argument and are proven bitwise against their originals on the
41-point grid by (G-copy) before any t*/80 use.  Asserted: (G-copy); (G125) /
(G-ctl) the schedule strings; (C100-s40) bitwise and P44's ``cont`` strings;
(R44-d25s) P44's ``d25`` strings; (Q-s) / (Win-s) on the D = 0.125 history at
both Δt, six points, peak and equal mean; max T < 43; (Rc′-s) at t*/40.
Printed only: (C100-s80), (X-grid), (Rc′-s) at t*/80, ratio / convexity /
sampling / window tables.
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import quad

from fem_em_solver.thermal.dose import cem43, running_average
from fem_em_solver.utils.analytical import robin_box_quadratic_transient_temperature
from tests.validation.test_mat13_transient_perfusion_sphere import INTERIOR_BAND
from tests.validation.test_mat20_cem43_closed_form import _rel
from tests.validation.test_mat20_tuned_heating_dose import R_BELOW, _quad_dose
from tests.validation.test_mat21_box_closed_form import HALF
from tests.validation.test_mat21_box_closed_form import RHO_C as RHO_C_0
from tests.validation.test_mat21_box_closed_form import T_INF as T_INF_0
from tests.validation.test_mat21_box_closed_form import H, K, Q
from tests.validation.test_mat21_box_dose_closed_form import (
    AVG_BAR,
    DOSE_BAR,
    DOSE_PTS,
    HIST_LABELS,
    HIST_PTS,
    LN_R,
    T_GRID,
    T_REF,
    W_B2,
    _q_check,
    _series_hist,
    _win_check,
)
from tests.validation.test_mat21_box_duty25_dose_closed_form import (
    MID_T,
    SCHEDULE_25,
    T_SW25,
    _chi25,
    _sched_string,
)
from tests.validation.test_mat21_box_duty_closed_form import DT_BOX, T_STAR
from tests.validation.test_mat21_box_duty_dose_closed_form import (
    P39P_SERIES_DOSE,
    _duty_hist,
    _trap,
)
from tests.validation.test_mat21_box_probe_dose_closed_form import NEW_LABELS, NEW_PTS, SIX_LABELS
from tests.validation.test_mat21_box_probe_duty_dose_closed_form import P43_SERIES_DOSE
from tests.validation.test_mat21_box_quadratic_closed_form import QUAD_TERMS, WEIGHTS
from tests.validation.test_mat21_box_quadratic_duty_closed_form import _quad_duty

DUTY125 = 0.125
# Physical switch times built as k * DT_BOX (P42's grid-ULP trap); one tuple for both grids.
T_SW125 = tuple(k * DT_BOX for k in (0, 1, 8, 9, 16, 17, 24, 25, 32, 33))
SCHEDULE_125 = "11" + "0" * 7 + ("1" + "0" * 7) * 4
SCHEDULE_125_80 = "111" + "0" * 14 + ("11" + "0" * 14) * 4
N80 = 80
DT80 = T_STAR / N80
T_GRID80 = DT80 * np.arange(N80 + 1)
# The solver's 81 (half-)step midpoints at N = 80: two start-up half-steps, then k = 2 ... 80.
MID_T80 = (0.25 * DT80, 0.75 * DT80) + tuple((k - 0.5) * DT80 for k in range(2, 81))
# P44's W1 D = 0.25 series doses (peak) at :.15e, SIX_LABELS order, copied from
# docs/testing/logs/20261021T005053Z_MAT-21-step16-W1-probe-duty-dose-closed-form.log:181-186.
P44_D25_SERIES_DOSE = (
    "7.275169262715064e-11",
    "7.275169262715064e-11",
    "4.047569053358871e-11",
    "4.039333276168868e-11",
    "4.039333276168868e-11",
    "4.094517681689437e-11",
)
# P41 / P42 / P44 recorded equal-mean ratios (printed beside), SIX_LABELS order (wall, wall,
# centre from P41 / P42 HIST_PTS logs; P3-P5 from P44 W1).
REC_EQ = {"d50": ("1.186312", "1.186312", "1.004475"), "d25": ("1.276608", "1.276608", "1.006916")}
# Predictions written first (brief), printed only.  Index: 0 walls, 1 centre, 2 +-x, 3 P5.
PRED_EQ = ((1.28, 1.5), (1.006, 1.015), (1.005, 1.012), (1.010, 1.022))
PRED_PEAK = ((0.06, 0.11), (0.970, 0.980), (0.968, 0.985), (0.915, 0.945))
PRED_RC = (1e6, 1e7)
PRED_XGRID = 1e-14
SIX = tuple(("h", j, HIST_LABELS[j]) for j in range(len(DOSE_PTS))) + tuple(
    ("n", j, NEW_LABELS[j]) for j in range(len(NEW_PTS))
)
SIX_KIND = (0, 0, 1, 2, 2, 3)
BATCHES = (("h", HIST_PTS), ("n", NEW_PTS))


def _p(msg):
    print(msg, flush=True)


def _chi125(t):
    """1 on [jP, jP + P/8), 0 otherwise (P = 8 DT_BOX); t is always a (half-)step midpoint."""
    return 1.0 if int(np.floor(t / DT_BOX)) % 8 < 1 else 0.0


def _chi125_planted(t):
    """(G-ctl): the on-phase one DT_BOX too long."""
    return 1.0 if int(np.floor(t / DT_BOX)) % 8 < 2 else 0.0


# --- grid-parameterised local copies (bodies verbatim; grid as an argument) ---------------


def _series_hist_g(pts, t, weights=WEIGHTS, n_terms=QUAD_TERMS):
    """``_series_hist`` on grid ``t``; row 0 is θ(0) = 0 exactly."""
    P = np.asarray(pts, dtype=float)
    out = np.zeros((t.size, P.shape[0]))
    for n in range(1, t.size):
        out[n] = robin_box_quadratic_transient_temperature(
            P[:, 0],
            P[:, 1],
            P[:, 2],
            t[n],
            half_widths=HALF,
            k=K,
            q=Q,
            h=H,
            rho_c=RHO_C_0,
            weights=weights,
            t_inf=0.0,
            n_terms=n_terms,
        )
    return out


def _duty_hist_g(pts, t, switch_times):
    """``_duty_hist`` on grid ``t``; row 0 = 0.0."""
    P = np.asarray(pts, dtype=float)
    out = np.zeros((t.size, P.shape[0]))
    for n in range(1, t.size):
        out[n] = _quad_duty(
            P[:, 0], P[:, 1], P[:, 2], t[n], switch_times=switch_times, n_terms=QUAD_TERMS
        )
    return out


def _q_check_g(T, t):
    """``_q_check`` on grid ``t``: (dose, quad reference, rel, swap, swap miss)."""
    dose = cem43(t, T)
    ref = _quad_dose(t, T, LN_R)
    swap = cem43(t, T, r_above=R_BELOW, r_below=0.5)
    return dose, ref, _rel(dose, ref), swap, _rel(swap, ref)


def _win_check_g(th, t, window=float(W_B2)):
    """``_win_check`` on grid ``t``: (mx, t_mx, lo, hi, ref, rel)."""
    _means, mx, t_mx = running_average(t, th, window=window)
    s = t_mx - window
    inside = t[(t > s) & (t <= t_mx)]
    pts = np.concatenate([[np.interp(s, t, th)], np.interp(inside, t, th)])
    lo, hi = float(np.min(pts)), float(np.max(pts))
    brk = [x for x in t if s < x < t_mx]
    ref, _err = quad(
        lambda x: np.interp(x, t, th),
        s,
        t_mx,
        points=brk or None,
        epsabs=0.0,
        epsrel=1e-13,
        limit=200,
    )
    ref /= window
    return float(mx), float(t_mx), lo, hi, ref, _rel(mx, ref)


def _trap_g(y, t):
    return float(np.sum(0.5 * (y[1:] + y[:-1]) * np.diff(t)))


def _sched_string_g(chi, mids):
    return "".join(str(int(chi(t))) for t in mids)


def _T_g(th):
    return T_INF_0 + np.asarray(th, dtype=float)


def _bd_g(th_h, th_s, j, t):
    """``_bd`` (test_mat21_box_dose_fem.py:187-193) on grid ``t``."""
    d_h = cem43(t, _T_g(th_h[:, j]))
    d_s = cem43(t, _T_g(th_s[:, j]))
    m = abs(d_h / d_s - 1.0)
    M = float(np.max(np.abs(th_h[:, j] - th_s[:, j])))
    cap = INTERIOR_BAND * float(np.max(th_s[:, j]))
    return d_h, d_s, m, M, cap


_CACHE = {}


def _hists():
    """Module cache: {dt: {sched: {batch: array}}}; dt in (40, 80), sched in cont/d125/s1/d25."""
    if not _CACHE:
        c40, c80 = {}, {}
        for b, P in BATCHES:
            c40.setdefault("cont", {})[b] = _series_hist(P)
            c40.setdefault("d125", {})[b] = _duty_hist(P, switch_times=T_SW125)
            c40.setdefault("s1", {})[b] = _duty_hist(P, switch_times=(0.0,))
            c40.setdefault("d25", {})[b] = _duty_hist(P, switch_times=T_SW25)
            c80.setdefault("cont", {})[b] = _series_hist_g(P, T_GRID80)
            c80.setdefault("d125", {})[b] = _duty_hist_g(P, T_GRID80, T_SW125)
            c80.setdefault("s1", {})[b] = _duty_hist_g(P, T_GRID80, (0.0,))
            c80.setdefault("d25", {})[b] = _duty_hist_g(P, T_GRID80, T_SW25)
        _CACHE[40] = c40
        _CACHE[80] = c80
    return _CACHE


GRIDS = ((40, T_GRID), (80, T_GRID80))


def test_g_copy_local_copies_bitwise_on_41_point_grid():
    c = _hists()
    res = {}
    for b, P in BATCHES:
        res[f"_series_hist_g {b}"] = bool(
            np.array_equal(_series_hist_g(P, T_GRID), c[40]["cont"][b])
        )
        res[f"_duty_hist_g {b}"] = bool(
            np.array_equal(_duty_hist_g(P, T_GRID, T_SW125), c[40]["d125"][b])
        )
    th = c[40]["d125"]["n"][:, 2]
    res["_q_check_g"] = _q_check_g(T_INF_0 + th, T_GRID) == _q_check(T_INF_0 + th)
    res["_win_check_g"] = _win_check_g(th, T_GRID) == _win_check(th)
    res["_trap_g"] = _trap_g(th, T_GRID) == _trap(th)
    res["_sched_string_g"] = _sched_string_g(_chi25, MID_T) == _sched_string(_chi25)
    from tests.validation.test_mat21_box_dose_fem import _bd

    th_s = c[40]["cont"]["h"]
    th_h = c[40]["d125"]["h"]
    res["_bd_g"] = all(_bd_g(th_h, th_s, j, T_GRID) == _bd(th_h, th_s, j) for j in range(3))
    _p(f"MAT21-17s (G-copy) local copies vs imported originals on the 41-point grid: {res}")
    assert all(res.values()), res


def test_g125_schedule_strings_and_planted_control():
    s40 = _sched_string_g(_chi125, MID_T)
    s80 = _sched_string_g(_chi125, MID_T80)
    _p(
        f"MAT21-17s (G125) N=40 {s40} (len {len(s40)}, on {s40.count('1')}); expected "
        f"{SCHEDULE_125}"
    )
    _p(
        f"MAT21-17s (G125) N=80 {s80} (len {len(s80)}, on {s80.count('1')}); expected "
        f"{SCHEDULE_125_80}"
    )
    assert s40 == SCHEDULE_125 and s80 == SCHEDULE_125_80
    assert (len(s40), s40.count("1"), len(s80), s80.count("1")) == (41, 6, 81, 11)
    p40 = _sched_string_g(_chi125_planted, MID_T)
    p80 = _sched_string_g(_chi125_planted, MID_T80)
    _p(
        f"MAT21-17s (G-ctl) planted (% 8 < 2) N=40 {p40} (== SCHEDULE_25 {p40 == SCHEDULE_25}); "
        f"N=80 {p80}"
    )
    assert p40 != SCHEDULE_125 and p80 != SCHEDULE_125_80
    assert p40 == SCHEDULE_25


def test_c100_s40_bitwise_and_p44_cont_strings():
    c = _hists()
    eqs = {b: bool(np.array_equal(c[40]["s1"][b], c[40]["cont"][b])) for b, _P in BATCHES}
    _p(f"MAT21-17s (C100-s40) switch_times=(0.0,) array_equal _series_hist over 41 rows: {eqs}")
    assert all(eqs.values()), eqs
    got = tuple(f"{cem43(T_GRID, T_INF_0 + c[40]['s1'][b][:, j]):.15e}" for b, j, _l in SIX)
    ref = tuple(P39P_SERIES_DOSE) + tuple(P43_SERIES_DOSE)
    _p(f"MAT21-17s (C100-s40) doses(100 %) {got} (P44 cont {ref})")
    assert got == ref
    e80 = {b: bool(np.array_equal(c[80]["s1"][b], c[80]["cont"][b])) for b, _P in BATCHES}
    _p(f"MAT21-17s (C100-s80) t*/80 array_equal (printed only; predicted True): {e80}")


def test_r44_d25s_series_dose_strings():
    c = _hists()
    got = tuple(f"{cem43(T_GRID, T_INF_0 + c[40]['d25'][b][:, j]):.15e}" for b, j, _l in SIX)
    _p(f"MAT21-17s (R44-d25s) t*/40 D = 0.25 series doses {got} (P44 {P44_D25_SERIES_DOSE})")
    assert got == P44_D25_SERIES_DOSE


def test_q_s_win_s_d125_both_dt_six_points():
    c = _hists()
    bad = []
    for n, t in GRIDS:
        for scale, lab in ((1.0, "peak"), (1.0 / DUTY125, "equal-mean")):
            for b, j, name in SIX:
                th = c[n]["d125"][b][:, j] * scale
                T = T_INF_0 + th
                for key in ("cont", "d25", "d125"):
                    sc = 1.0 if key == "cont" else scale
                    assert float(np.max(T_INF_0 + c[n][key][b][:, j] * sc)) < T_REF, (n, key, name)
                dose, ref, rel, swap, miss = _q_check_g(T, t)
                mx, t_mx, lo, hi, wref, wrel = _win_check_g(th, t)
                _p(
                    f"MAT21-17s (Q-s) t*/{n} d125 {lab} {name}: CEM43 {dose:.15e} s vs quad "
                    f"{ref:.15e} s rel {rel:.3e} (bar {DOSE_BAR:g}, ASSERTED; predicted <= 1e-14); "
                    f"max T {np.max(T):.9f} C < 43 ASSERTED"
                )
                _p(
                    f"MAT21-17s (Win-s) t*/{n} d125 {lab} {name}: max six-minute mean {mx:.15e} K "
                    f"at t_max {t_mx:.4f} s; (ii) [{lo:.12e}, {hi:.12e}]; (iii) quad {wref:.15e} "
                    f"rel {wrel:.3e} (bar {AVG_BAR:g}, ASSERTED; predicted <= 1e-15)"
                )
                if not (
                    rel <= DOSE_BAR
                    and t[0] + W_B2 <= t_mx <= t[-1]
                    and lo <= mx <= hi
                    and wrel <= AVG_BAR
                ):
                    bad.append((n, lab, name, rel, t_mx, lo, mx, hi, wrel))
    assert not bad, bad


def test_rc_s_controls():
    c = _hists()
    bad = []
    for n, t in GRIDS:
        for scale, lab in ((1.0, "peak"), (1.0 / DUTY125, "equal-mean")):
            for b, j, name in SIX:
                _d, ref, _r, swap, miss = _q_check_g(T_INF_0 + c[n]["d125"][b][:, j] * scale, t)
                tag = (
                    f"asserted > {DOSE_BAR:g}"
                    if n == 40
                    else f"PREDICTED factor {PRED_RC}; printed only, ruling 03:00 (6)"
                )
                _p(
                    f"MAT21-17s (Rc'-s) t*/{n} d125 {lab} {name}: R swap miss {miss:.4e} ({tag}); "
                    f"factor {swap / ref:.4e}"
                )
                if n == 40 and not miss > DOSE_BAR:
                    bad.append((lab, name, miss))
    assert not bad, bad


def test_printed_tables():
    """The record; asserts only that it was formed."""
    c = _hists()
    W = float(W_B2)
    for key in ("cont", "d125"):
        for b, _P in BATCHES:
            a40 = c[40][key][b]
            a80 = c[80][key][b]
            x = float(np.max(np.abs(a80[::2] - a40)) / np.max(a40))
            _p(
                f"MAT21-17s (X-grid) {key} {b}: max_n |theta80[2n] - theta40[n]| / max theta40 "
                f"{x:.3e} (predicted <= {PRED_XGRID:g}); T_GRID80[::2] == T_GRID "
                f"{bool(np.array_equal(T_GRID80[::2], T_GRID))}"
            )
    eqs = {}
    for n, t in GRIDS:
        for b, j, name in SIX:
            k = SIX_KIND[SIX.index((b, j, name))]
            th_s = c[n]["cont"][b][:, j]
            d_c = cem43(t, T_INF_0 + th_s)
            for key, D in (("d25", 0.25), ("d125", DUTY125)):
                th_d = c[n][key][b][:, j]
                d_pk = cem43(t, T_INF_0 + th_d)
                d_eq = cem43(t, T_INF_0 + th_d / D)
                eqs[(n, key, name)] = d_eq
                if key == "d25":
                    continue
                r_t = th_d[-1] / th_s[-1]
                r_mean = _trap_g(th_d, t) / _trap_g(th_s, t)
                _a, mx_c, t_c = running_average(t, th_s, window=W)
                _a, mx_d, t_d = running_average(t, th_d, window=W)
                rec = REC_EQ["d50"][j] + " / " + REC_EQ["d25"][j] if b == "h" else "P44 W1 log"
                _p(
                    f"MAT21-17s DUTY DOSE t*/{n} d125 {name}: D_cont {d_c:.15e} s; D(peak) "
                    f"{d_pk:.15e} s; D(equal mean) {d_eq:.15e} s; ratio equal-mean "
                    f"{d_eq / d_c:.6f} (predicted {PRED_EQ[k]}; d50 / d25 recorded {rec}); ratio "
                    f"peak {d_pk / d_c:.6f} (predicted {PRED_PEAK[k]}); theta ratio at t* "
                    f"{r_t:.5f} (x8 {r_t / DUTY125:.5f}); time-mean theta ratio {r_mean:.5f} "
                    f"(x8 {r_mean / DUTY125:.5f}); six-minute peak cont {mx_c:.6f} K at {t_c:.2f} "
                    f"s, d125 {mx_d:.6f} K at {t_d:.2f} s ({t_d / DT_BOX:.2f} DT_BOX)"
                )
                assert np.isfinite(d_eq / d_c)
    for n, _t in GRIDS:
        for _b, _j, name in SIX:
            _p(
                f"MAT21-17s convexity t*/{n} D12.5,eq/D25,eq {name}: "
                f"{eqs[(n, 'd125', name)] / eqs[(n, 'd25', name)]:.6f}"
            )
    for key in ("cont", "d25", "d125"):
        for b, j, name in SIX:
            D = {"cont": 1.0, "d25": 0.25, "d125": DUTY125}[key]
            d41 = cem43(T_GRID, T_INF_0 + c[40][key][b][:, j])
            d81 = cem43(T_GRID80, T_INF_0 + c[80][key][b][:, j])
            _p(
                f"MAT21-17s SAMPLING (41 vs 81 samples) {key} {name}: D41 {d41:.15e} s, D81 "
                f"{d81:.15e} s, D41/D81 - 1 {d41 / d81 - 1.0:+.4e} (D = {D})"
            )
    for n, t in GRIDS:
        for norm in ("peak", "equal-mean"):
            for b, j, name in SIX:
                row = []
                for key, D in (("cont", 1.0), ("d25", 0.25), ("d125", DUTY125)):
                    sc = 1.0 if norm == "peak" else 1.0 / D
                    row.append(f"{key} {cem43(t, T_INF_0 + c[n][key][b][:, j] * sc):.15e}")
                _p(
                    f"MAT21-17s SIX-POINT x THREE-SCHEDULE SERIES DOSE TABLE t*/{n} ({norm}) "
                    f"{name}: " + "; ".join(row) + " s"
                )
    for n, _t in GRIDS:
        a = c[n]["d125"]["n"]
        _p(
            f"MAT21-17s P3 / P4 d125 t*/{n} series histories np.array_equal: "
            f"{bool(np.array_equal(a[:, 0], a[:, 1]))} (predicted True)"
        )
    _p(f"MAT21-17s labels {SIX_LABELS}; SIX order {[x[2] for x in SIX]}")
