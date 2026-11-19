"""`MAT-21` step 45 — step 6's end-weighted source per deposited watt at entry 59's six probe
points and on C6's 18-point lattice, beside P75's uniform U3 and P80's eddy source, closed form
(§9 item 921, §10 P81; brief `docs/planning/chunks/MAT-21.md` "§9 item 921").

Numpy (complex build only because P80's module imports P79 -> P61 / C6).  Sets no band and
moves none; every band is imported.  Step 6's ratified q = Q[1 - cos(pi x / A)] (``PROFILES``)
is read per deposited watt (/ ``p75.P_DEP``; (P-eq) asserts all three pictures deposit the same
watt) at P75's ``PROBE_PTS`` as one batch and at C6's lattice (``p80._data()["lat_pts"]``) as one
batch.  Asserted (collected, then asserted — item 668): preconditions incl. (P-eq); (Q2-s) the
uniform profile is P75's U3 bitwise on P75's batch; (Rep-s) step 6's strings on step 6's own
scalar calls; (Lim); (pw); (Lin) superposition of PROFILES + FLIPPED = 2 UNIFORM; MAT-20's
(Q-d) / (Win-d) on the per-watt end histories at the six probes and at c*; (Cem) at 24
points; (Ball) the closed-form 10 g-ball mean against a 64-node cubature; controls (c-u),
(c-e), (c-ball).  Entry 52's SAR10g shape is printed beside the end-source ball-mean shape,
never gated.

Not claimed: that the clinical SAR is end-weighted, eddy-like or any mixture, which picture
entry 59 follows, any entry-59 / entry-60 figure or ratio to them, a fitted source, a dose
band, a compliance figure, a 10 g-averaged theta.  Per deposited gel watt throughout; entry 59
is per accepted watt (theta_acc = theta * P_gel / P_acc).
"""

from __future__ import annotations

import numpy as np

from fem_em_solver.thermal.dose import running_average
from fem_em_solver.utils.analytical import (
    robin_box_separable_steady_temperature,
    robin_box_separable_transient_temperature,
)
from tests.validation import test_geo38_f_clinical_box_loaded as c6
from tests.validation import test_mat21_eddy_source_probe_points_closed_form as p80
from tests.validation import test_mat21_wall_uniform_source_probe_points_closed_form as p75
from tests.validation.test_imp4_pec_cylinder_closed_form import B1_BAND
from tests.validation.test_mat13_transient_perfusion_sphere import INTERIOR_BAND
from tests.validation.test_mat21_box_closed_form import (
    A_BOX,
    B_BOX,
    C_BOX,
    HALF,
)
from tests.validation.test_mat21_box_closed_form import RHO_C as RHO_C_0
from tests.validation.test_mat21_box_closed_form import SERIES_TERMS, T_C, T_INF, H, K, Q
from tests.validation.test_mat21_box_dose_closed_form import AVG_BAR, DOSE_BAR, T_GRID, W_B2
from tests.validation.test_mat21_box_dt_duty125_dose_closed_form import _q_check_g, _win_check_g
from tests.validation.test_mat21_box_separable_closed_form import (
    FLIPPED,
    PROFILES,
    SEP_TERMS,
    UNIFORM,
    _sep,
    _sep_t,
)

EPS = float(np.finfo(float).eps)
# Step 6's printed records (`20261015T130111Z_MAT-21-step6-W1-closed-form.log.gz:198-199`).
STEP6_STEADY = ("2.741147", "14.919189", "8.541895")
STEP6_TSTAR = ("0.024797", "2.537123")
# U3's (40,40,40) vs (80,80,80) truncation at t* (K), printed only
# (`20261108T134258Z_MAT-21-P75-W1.log:93`).
P75_BQ3D = (7.164e-06, 5.391e-05, 5.391e-05, 1.764e-04, 1.764e-04, 9.008e-06)
BQ_TERMS = (160, 80, 80)
# Predictions written first (brief), printed only; (value, lo, hi) with None = no range.
PRED_CU = (
    (0.9814, 0.978, 0.985),
    (0.98, 0.978, 0.985),
    (0.98, 0.978, 0.985),
    (0.98, 0.90, 1.00),
    (0.98, 0.90, 1.00),
    (0.98, 0.978, 0.985),
)
PRED_CE = (
    (0.609, None, None),
    (0.9925, 0.991, 0.994),
    (0.9925, 0.991, 0.994),
    (42.5, 40.7, 42.9),
    (42.5, 40.7, 42.9),
    (0.80, 0.76, 0.84),
)
PRED_CBALL = {50: 8.482e-03, 150: 1.523e-04, 250: 4.761e-04}
PRED_LAT_EU = {50: (0.131, 0.12, 0.14), 150: (0.882, 0.86, 0.90), 250: (1.735, 1.70, 1.76)}
PRED_SHAPE = {50: 1.0, 150: 7.582, 250: 15.06}
PRED_COL = (0.0664, 0.5035)
PRED_TH_P0 = 5.0462e-04
PRED_EU360_P0 = 2.4e-03
PRED_CSTAR_EE = (0.171, 0.157, 0.183)

_CACHE = {}


def _say(msg):
    print(msg, flush=True)


def _probe_arr():
    return np.asarray(p75.PROBE_PTS, dtype=float)


def _end_at(P, t, q=Q, profiles=PROFILES, n_terms=SEP_TERMS):
    P = np.asarray(P, dtype=float)
    return np.asarray(
        robin_box_separable_transient_temperature(
            P[:, 0],
            P[:, 1],
            P[:, 2],
            t,
            half_widths=HALF,
            k=K,
            q=q,
            h=H,
            rho_c=RHO_C_0,
            profiles=profiles,
            t_inf=0.0,
            n_terms=n_terms,
        ),
        dtype=float,
    )


def _hist(P, profiles=PROFILES):
    """(41, n) absolute history on T_GRID, one batch per row; row 0 is 0 exactly."""
    out = np.zeros((T_GRID.size, len(P)))
    for n in range(1, T_GRID.size):
        out[n] = _end_at(P, T_GRID[n], profiles=profiles)
    return out


def _data():
    if not _CACHE:
        d80 = p80._data()
        P = _probe_arr()
        lp = d80["lat_pts"]
        s_abs, sl_abs = _hist(P), _hist(lp)
        u_abs, ul_abs = _hist(P, UNIFORM), _hist(lp, UNIFORM)
        _CACHE.update(
            S_abs=s_abs,
            SL_abs=sl_abs,
            S=s_abs / p75.P_DEP,
            SL=sl_abs / p75.P_DEP,
            UU=u_abs / p75.P_DEP,
            UUL=ul_abs / p75.P_DEP,
            S360=_end_at(P, float(W_B2)) / p75.P_DEP,
            d80=d80,
        )
    return _CACHE


def _q_end(x):
    return Q * (1.0 - np.cos(np.pi * np.asarray(x, dtype=float) / A_BOX))


def _ball_mean_end(p, r, shell=False):
    """Closed-form mean of q_end over the ball (p, r): Q[1 - cos(pi p_a / A) * f(u)],
    u = pi r / A, f = 3 j1(u) / u (ball) or sin u / u (``shell``, the (c-ball) control)."""
    u = np.pi * r / A_BOX
    f = np.sin(u) / u if shell else 3.0 * (np.sin(u) - u * np.cos(u)) / u**3
    return float(Q * (1.0 - np.cos(np.pi * float(p[0]) / A_BOX) * f))


def _ball_cub_end(p, r):
    """Polar axis along a: x_a = p_a + r s mu; 8-node GL in s (weight s^2), 8-node GL in mu."""
    xs, ws = np.polynomial.legendre.leggauss(8)
    s, w_s = 0.5 * (xs + 1.0), 0.5 * ws * (0.5 * (xs + 1.0)) ** 2
    mu, w_mu = np.polynomial.legendre.leggauss(8)
    S, MU = np.meshgrid(s, mu, indexing="ij")
    W = np.outer(w_s, w_mu)
    return float(np.sum(W * _q_end(float(p[0]) + r * S * MU)) / np.sum(W))


def _dz(lab):
    return abs(int(lab.split("dz")[1].rstrip("mm")))


def _in(v, lo, hi):
    return lo is None or (lo <= v <= hi)


# ------------------------------------------------------------------ tests
def test_preconditions_p_eq():
    gates = []
    d = _data()
    d80 = d["d80"]
    tg = (p75.T_GRID is T_GRID) or bool(np.array_equal(p75.T_GRID, T_GRID))
    gates.append((("T_GRID",), tg, 0))
    dts = abs(p75.T_STAR - 0.2 * T_C)
    gates.append((("T_STAR", p75.T_STAR, 0.2 * T_C), dts <= 1e-9, dts))
    pr_ok = PROFILES == ((1.0, -1.0, 1), (1.0, 0.0, 0), (1.0, 0.0, 0))
    gates.append((("PROFILES", PROFILES), pr_ok, 0))
    m = 1.0
    for a, b, j in PROFILES:
        m *= a if b == 0.0 else a + b * np.sin(j * np.pi) / (j * np.pi)
    p_end = Q * 8.0 * A_BOX * B_BOX * C_BOX * m
    peq = abs(p_end / p75.P_DEP - 1.0)
    gates.append((("P-eq", p_end, p75.P_DEP), peq <= 4 * EPS, peq))
    lat52, _star = p80._parse_c52()
    gates.append((("c52 n", len(lat52)), len(lat52) == 18, len(lat52)))
    ci = p80._cstar_idx()
    gates.append((("c*", ci), d80["lat"][ci][0] == p80.C52_STAR, ci))
    r = float(c6._a10())
    worst = -np.inf
    for (lab, pt, why), p in zip(d80["lat"], d80["lat_pts"]):
        if pt is None:
            gates.append((("lattice", lab, why), False, None))
            continue
        clear = max(abs(p[k]) + r - HALF[k] for k in range(3))
        worst = max(worst, clear)
        gates.append((("ball inside", lab, clear), clear <= 0.0, clear))
    gates.append((("lattice n",), len(d80["lat"]) == 18, len(d80["lat"])))
    gates.append(
        (("row0",), bool(np.all(d["S_abs"][0] == 0.0) and np.all(d["SL_abs"][0] == 0.0)), 0)
    )
    _say(
        f"MAT21-45 T_GRID is p75.T_GRID {p75.T_GRID is T_GRID} (ASSERTED); p75.T_STAR "
        f"{p75.T_STAR:.6f} s vs 0.2 T_C {0.2 * T_C:.6f} s, |d| {dts:.3e} s (bound 1e-9, ASSERTED)"
    )
    _say(f"MAT21-45 PROFILES {PROFILES} == ((1,-1,1),(1,0,0),(1,0,0)) {pr_ok} (ASSERTED)")
    _say(
        f"MAT21-45 (P-eq) end P_dep {p_end:.6f} W vs p75.P_DEP {p75.P_DEP:.6f} W, |ratio-1| "
        f"{peq:.3e} (bound 4 eps {4 * EPS:.3e}, ASSERTED)"
    )
    _say(
        f"MAT21-45 entry 52 parse: {len(lat52)} lattice matches (18 ASSERTED); c* idx {ci} -> "
        f"{d80['lat'][ci][0]!r} == {p80.C52_STAR!r} (ASSERTED)"
    )
    _say(
        f"MAT21-45 lattice {len(d80['lat'])} points; 10 g ball r {r:.6f} m, worst max_d(|p_d| + r "
        f"- L_d) {worst:.6f} m (<= 0 ASSERTED); row 0 of both histories 0 exactly (ASSERTED)"
    )
    bad = [g for g in gates if not g[1]]
    assert not bad, bad


def test_q2s_rep_s_lim_pw_lin():
    gates = []
    P = _probe_arr()
    ts = p75.T_STAR
    # (Q2-s) on P75's own PROBE_PTS batch
    qu = _end_at(P, ts, profiles=UNIFORM, n_terms=SERIES_TERMS)
    u3 = p75._u3(ts)
    q2 = bool(np.array_equal(qu, u3))
    dq = float(np.max(np.abs(qu - u3)))
    gates.append((("Q2-s", dq), q2, dq))
    _say(
        f"MAT21-45 (Q2-s) separable(UNIFORM, SERIES_TERMS) vs p75._u3 at t* on PROBE_PTS: "
        f"np.array_equal {q2}, max |diff| {dq:.3e} K (ASSERTED bitwise)"
    )
    # (Rep-s) on step 6's own scalar calls
    pts = ((0.0, 0.0, 0.0), (0.9 * A_BOX, 0.0, 0.0), (0.9 * A_BOX, 0.9 * B_BOX, 0.9 * C_BOX))
    st = tuple(f"{float(_sep(*p, n_terms=SEP_TERMS) - T_INF):.6f}" for p in pts)
    tt = tuple(
        f"{float(_sep_t(x, 0.0, 0.0, 0.2 * T_C, n_terms=SEP_TERMS) - T_INF):.6f}"
        for x in (0.0, 0.9 * A_BOX)
    )
    gates.append((("Rep-s steady", st), st == STEP6_STEADY, 0))
    gates.append((("Rep-s t*", tt), tt == STEP6_TSTAR, 0))
    _say(
        f"MAT21-45 (Rep-s) step 6 steady (80,40,40) {st} vs {STEP6_STEADY}: {st == STEP6_STEADY}; "
        f"transient t* {tt} vs {STEP6_TSTAR}: {tt == STEP6_TSTAR} (strings ASSERTED, step 6's "
        f"own scalar _sep / _sep_t calls)"
    )
    # (Lim)
    ss = np.asarray(
        robin_box_separable_steady_temperature(
            P[:, 0],
            P[:, 1],
            P[:, 2],
            half_widths=HALF,
            k=K,
            q=Q,
            h=H,
            profiles=PROFILES,
            t_inf=0.0,
            n_terms=SEP_TERMS,
        ),
        dtype=float,
    )
    tl = _end_at(P, 1000.0 * T_C)
    lim = float(np.max(np.abs(tl - ss) / np.abs(ss)))
    gates.append((("Lim",), lim <= B1_BAND, lim))
    _say(
        f"MAT21-45 (Lim) end transient(1000 T_C) vs separable steady (SEP_TERMS), six probes: max "
        f"rel {lim:.3e} (band {B1_BAND:g}, ASSERTED; step 6 (Q6) read 1.296e-15 at 100 t_C)"
    )
    # (pw)
    t1 = _end_at(P, ts) / p75.P_DEP
    t2 = _end_at(P, ts, q=2.0 * Q) / (2.0 * p75.P_DEP)
    pw = float(np.max(np.abs(t2 - t1) / np.abs(t1)))
    gates.append((("pw",), pw <= B1_BAND, pw))
    _say(
        f"MAT21-45 (pw) theta(2Q)/(2P) vs theta(Q)/P at t*, six probes: max rel {pw:.3e} (band "
        f"{B1_BAND:g}, ASSERTED)"
    )
    # (Lin) one batch for all three profiles
    tp = _end_at(P, ts, profiles=PROFILES)
    tf = _end_at(P, ts, profiles=FLIPPED)
    tu = _end_at(P, ts, profiles=UNIFORM)
    lin = float(np.max(np.abs(tp + tf - 2.0 * tu)) / np.max(np.abs(tu)))
    gates.append((("Lin",), lin <= B1_BAND, lin))
    _say(
        f"MAT21-45 (Lin) max|theta(PROFILES) + theta(FLIPPED) - 2 theta(UNIFORM)| / max "
        f"theta(UNIFORM) at t*, SEP_TERMS, PROBE_PTS batch: {lin:.3e} (band {B1_BAND:g}, ASSERTED)"
    )
    bad = [g for g in gates if not g[1]]
    assert not bad, bad


def test_mat20_forms_and_cem_bracket():
    gates = []
    d = _data()
    d80 = d["d80"]
    ci = p80._cstar_idx()
    rows = [(lab, d["S"][:, i]) for i, lab in enumerate(p75.PROBE_LABELS)]
    rows.append((f"c* {p80.C52_STAR}", d["SL"][:, ci]))
    worst = {"q": 0.0, "w": 0.0, "rc": np.inf}
    for lab, th in rows:
        _dose, _ref, rel, _swap, miss = _q_check_g(T_INF + th, T_GRID)
        mx, t_mx, lo, hi, _r, wrel = _win_check_g(th, T_GRID)
        gates.append(((lab, "Q"), rel <= DOSE_BAR, rel))
        gates.append(((lab, "swap"), miss > DOSE_BAR, miss))
        gates.append(
            ((lab, "win-loc"), lo <= mx <= hi and T_GRID[0] + W_B2 <= t_mx <= T_GRID[-1], t_mx)
        )
        gates.append(((lab, "Win"), wrel <= AVG_BAR, wrel))
        worst["q"], worst["w"] = max(worst["q"], rel), max(worst["w"], wrel)
        worst["rc"] = min(worst["rc"], miss)
    _say(
        f"MAT21-45 (Q-d) worst rel {worst['q']:.3e} (bar {DOSE_BAR:g}, ASSERTED); (Win-d) worst "
        f"rel {worst['w']:.3e} (bar {AVG_BAR:g}, ASSERTED), window located (ASSERTED); R swap least "
        f"miss {worst['rc']:.4e} (asserted > {DOSE_BAR:g}); per-watt end source at 6 probes + c*"
    )
    allrows = [(lab, d["S"][:, i]) for i, lab in enumerate(p75.PROBE_LABELS)]
    allrows += [(lab, d["SL"][:, j]) for j, (lab, _p, _w) in enumerate(d80["lat"])]
    worst_m, n_ok = np.inf, 0
    for lab, th in allrows:
        cem, lo, hi = p80._cem_bracket(th)
        ok = lo <= cem <= hi
        n_ok += int(ok)
        gates.append(((lab, "Cem", lo, cem, hi), ok, cem))
        worst_m = min(worst_m, min(cem - lo, hi - cem) / cem)
    _say(
        f"MAT21-45 (Cem) p80._cem_bracket per watt at {len(allrows)} points (6 probe + 18 "
        f"lattice): inside at {n_ok} / {len(allrows)} (all ASSERTED); least relative margin "
        f"{worst_m:.3e}"
    )
    bad = [g for g in gates if not g[1]]
    assert not bad, bad


def test_ball():
    gates = []
    d80 = _data()["d80"]
    r = float(c6._a10())
    worst = 0.0
    for (lab, _pt, _w), p in zip(d80["lat"], d80["lat_pts"]):
        bm, bc = _ball_mean_end(p, r), _ball_cub_end(p, r)
        e = abs(bc - bm) / bm
        worst = max(worst, e)
        gates.append(((lab, "Ball", bm, bc), e <= B1_BAND, e))
    _say(
        f"MAT21-45 (Ball) |cub - closed| / closed max over 18 lattice balls {worst:.3e} (band "
        f"{B1_BAND:g}, ASSERTED; brief expects near 1e-14); r {r:.6f} m; 64-node rule, u "
        f"{np.pi * r / A_BOX:.4f}"
    )
    bad = [g for g in gates if not g[1]]
    assert not bad, bad


def test_controls_c_u_c_e_c_ball():
    d = _data()
    d80 = d["d80"]
    gates = []
    end = d["S"][-1]
    for i, lab in enumerate(p75.PROBE_LABELS):
        rt = float(end[i] / d80["U3"][i])
        cu = abs(rt - 1.0)
        gates.append(((lab, "c-u", cu), cu > INTERIOR_BAND, cu))
        pv, lo, hi = PRED_CU[i]
        _say(
            f"MAT21-45 (c-u) {lab}: end/U3 at t* {rt:.6f}, |ratio-1| {cu:.4f} (ASSERTED > "
            f"{INTERIOR_BAND:g}); predicted {pv} range [{lo}, {hi}] inside {_in(cu, lo, hi)}; "
            f"separation {cu / INTERIOR_BAND:.1f}x"
        )
    for i, lab in enumerate(p75.PROBE_LABELS):
        re_ = float(end[i] / d80["E"][-1, i])
        ce = abs(re_ - 1.0)
        gates.append(((lab, "c-e", ce), ce > INTERIOR_BAND, ce))
        pv, lo, hi = PRED_CE[i]
        _say(
            f"MAT21-45 (c-e) {lab}: end/eddy at t* {re_:.6f}, |ratio-1| {ce:.4f} (ASSERTED > "
            f"{INTERIOR_BAND:g}); predicted {pv} range [{lo}, {hi}] inside {_in(ce, lo, hi)}"
        )
    r = float(c6._a10())
    by_dz = {}
    for (lab, _pt, _w), p in zip(d80["lat"], d80["lat_pts"]):
        bc = _ball_cub_end(p, r)
        miss = abs(_ball_mean_end(p, r, shell=True) / bc - 1.0)
        gates.append(((lab, "c-ball", miss), miss > B1_BAND, miss))
        by_dz.setdefault(_dz(lab), []).append(miss)
    for dz in sorted(by_dz):
        lo_, hi_ = min(by_dz[dz]), max(by_dz[dz])
        pv = PRED_CBALL[dz]
        ok = 0.95 * pv <= lo_ and hi_ <= 1.05 * pv
        _say(
            f"MAT21-45 (c-ball) shell sin u/u vs cub miss at dz+-{dz}mm (6 points): [{lo_:.4e}, "
            f"{hi_:.4e}] (all > {B1_BAND:g} ASSERTED); predicted {pv:.3e} +-5 % inside {ok}"
        )
    bad = [g for g in gates if not g[1]]
    assert not bad, bad


def test_record_printed():
    """The record (no band); asserts only that it was formed."""
    d = _data()
    d80 = d["d80"]
    P = _probe_arr()
    ts = p75.T_STAR
    u360 = p75._u3(float(W_B2)) / p75.P_DEP
    _say(
        f"MAT21-45 per deposited gel watt: P_dep {p75.P_DEP:.6f} W (all three sources); entry 59 "
        f"is per accepted watt, theta_acc = theta * P_gel / P_acc"
    )
    _say(
        "MAT21-45 PROBE point | end/P t* K/W | end/P 360 s | U3/P t* (P80) | eddy/P t* (P80) | "
        "end/U3 | end/eddy | cem43/W s | Cem lo | Cem hi | max 360 s mean/W K"
    )
    for i, lab in enumerate(p75.PROBE_LABELS):
        th = d["S"][:, i]
        cem, lo, hi = p80._cem_bracket(th)
        _m, mx, _t = running_average(T_GRID, th, window=float(W_B2))
        _say(
            f"MAT21-45 PROBE {lab} | {th[-1]:.7g} | {d['S360'][i]:.7g} | {d80['U3'][i]:.7g} | "
            f"{d80['E'][-1, i]:.7g} | {th[-1] / d80['U3'][i]:.7g} | {th[-1] / d80['E'][-1, i]:.7g}"
            f" | {cem:.7g} | {lo:.7g} | {hi:.7g} | {float(mx):.7g}"
        )
    _say(f"MAT21-45 predicted end/P at P0, t* {PRED_TH_P0:.4e} K/W; measured {d['S'][-1, 0]:.4e}")
    eu360 = float(d["S360"][0] / u360[0])
    _say(
        f"MAT21-45 end/U3 at 360 s at P0 {eu360:.4e} (predicted approx {PRED_EU360_P0:.1e}, "
        f"wall-free 1 - (1 - e^-kt)/(kt))"
    )
    # (Bq) truncation, printed only
    bq = np.abs(_end_at(P, ts, n_terms=BQ_TERMS) - _end_at(P, ts)) / p75.P_DEP
    for i, lab in enumerate(p75.PROBE_LABELS):
        _say(
            f"MAT21-45 (Bq) {lab}: |end(SEP_TERMS) - end{BQ_TERMS}| at t* {bq[i] * p75.P_DEP:.3e} K "
            f"({bq[i]:.3e} K/W); P75 U3 (40^3 vs 80^3) {P75_BQ3D[i]:.3e} K (printed only)"
        )
    # lattice
    labs = [lab for lab, _p, _w in d80["lat"]]
    ci = p80._cstar_idx()
    r = float(c6._a10())
    bms = np.array([_ball_mean_end(p, r) for p in d80["lat_pts"]])
    lat52, _star = p80._parse_c52()
    sar = {m[0]: float(m[2]) for m in lat52}
    eL, uL, edL = d["SL"][-1], d80["U3_lat"], d80["L"][-1]
    _say(
        "MAT21-45 LATTICE label | box (a, b, c) m | end/P t* K/W | end/U3 | end/eddy | "
        "ball-mean/ball-mean(c*) | entry 52 SAR10g/SAR10g(c*) ccw"
    )
    for j, lab in enumerate(labs):
        p = d80["lat_pts"][j]
        _say(
            f"MAT21-45 LATTICE {lab} | ({p[0]:+.3f}, {p[1]:+.3f}, {p[2]:+.3f}) | {eL[j]:.7g} | "
            f"{eL[j] / uL[j]:.7g} | {eL[j] / edL[j]:.7g} | {bms[j] / bms[ci]:.4g} | "
            f"{sar[lab] / sar[p80.C52_STAR]:.4g}"
        )
    for dz in (50, 150, 250):
        js = [j for j, lab in enumerate(labs) if _dz(lab) == dz]
        eu = [eL[j] / uL[j] for j in js]
        sh = [bms[j] / bms[ci] for j in js]
        spread = (max(eL[j] for j in js) - min(eL[j] for j in js)) / max(eL[j] for j in js)
        pv, lo, hi = PRED_LAT_EU[dz]
        _say(
            f"MAT21-45 dz+-{dz}mm: end/U3 [{min(eu):.4g}, {max(eu):.4g}] (predicted {pv} range "
            f"[{lo}, {hi}] inside {all(lo <= v <= hi for v in eu)}); end/P relative spread over "
            f"6 points {spread:.3e} (predicted approx 1e-5); ball-mean shape/c* [{min(sh):.4g}, "
            f"{max(sh):.4g}] (predicted {PRED_SHAPE[dz]})"
        )
    col = {lab: bms[j] for j, lab in enumerate(labs) if lab.startswith("x+0mm")}
    c50 = col["x+0mm dz-50mm"] / col["x+0mm dz-250mm"]
    c150 = col["x+0mm dz-150mm"] / col["x+0mm dz-250mm"]
    s_m = (
        sar["x+0mm dz-50mm"] / sar["x+0mm dz-250mm"],
        sar["x+0mm dz-150mm"] / sar["x+0mm dz-250mm"],
    )
    s_p = (
        sar["x+0mm dz+50mm"] / sar["x+0mm dz+250mm"],
        sar["x+0mm dz+150mm"] / sar["x+0mm dz+250mm"],
    )
    _say(
        f"MAT21-45 x+0mm column ball-mean ratios dz50/dz250 {c50:.4f}, dz150/dz250 {c150:.4f} "
        f"(predicted {PRED_COL}) beside entry 52 SAR10g - side {s_m[0]:.4f} / {s_m[1]:.4f}, + side "
        f"{s_p[0]:.4f} / {s_p[1]:.4f} (printed only)"
    )
    ee = float(eL[ci] / edL[ci])
    pv, lo, hi = PRED_CSTAR_EE
    _say(
        f"MAT21-45 end/eddy at c* {p80.C52_STAR!r} {ee:.4f} (predicted {pv}, range [{lo}, {hi}] "
        f"inside {lo <= ee <= hi}; printed only)"
    )
    near = [labs[j] for j in range(len(labs)) if eL[j] >= np.max(eL) * (1.0 - 1e-4)]
    pred_near = sorted(lab for lab in labs if _dz(lab) == 250)
    _say(
        f"MAT21-45 end tie set within 1e-4 rel of max: {near} (predicted the six dz+-250mm: "
        f"{sorted(near) == pred_near}); no exact argmax is read (P80 disclosure (1))"
    )
    # comparison bound 2 theta(UNIFORM) - theta_end over 24 points x 40 rows
    cmp_p = 2.0 * d["UU"][1:] - d["S"][1:]
    cmp_l = 2.0 * d["UUL"][1:] - d["SL"][1:]
    mn = float(min(np.min(cmp_p), np.min(cmp_l)))
    _say(
        f"MAT21-45 comparison bound min over 24 points x 40 rows of 2 theta(UNIFORM) - theta_end "
        f"per watt {mn:.3e} K/W, sign {'+' if mn >= 0 else '-'} (printed only; a negative value "
        f"would be truncation)"
    )
    assert np.all(np.isfinite(eL)) and bms[ci] > 0.0
