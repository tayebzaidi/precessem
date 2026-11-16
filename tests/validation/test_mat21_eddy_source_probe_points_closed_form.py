"""`MAT-21` step 44 — P36's eddy source per deposited watt at entry 59's six probe points
and on C6's 18-point lattice, beside P75's uniform U3, closed form (§9 item 913, §10 P80;
brief `docs/planning/chunks/MAT-21.md` "§9 item 913").

Numpy (complex build only because P79's module imports P61 / C6).  Sets no band and moves
none; every band is imported.  P36's long-prism CP eddy shape q = 3Q(y² + z²)/(B² + C²)
(``WEIGHTS``) is read per deposited watt (÷ ``p75.P_DEP``; (P-eq) asserts both sources
deposit the same watt) at P37's six points — P39′'s ``HIST_PTS`` batch and P43's
``NEW_PTS`` batch, each its own ``_series_hist`` call, permuted by ``MAP`` — and at C6's
lattice (mid-depth, c = 0), with P75's U3 beside.  Asserted (collected, then asserted —
item 668): preconditions incl. (P-eq); (Q2-p) U3 is the unit-weight member bitwise;
(Rep-s) P39′'s / P43's dose strings and the six θ_s(t*) strings on their own batches;
(Lim); (pw); MAT-20's (Q-d) / (Win-d) on the per-watt eddy histories at the six points and
at c*; (Cem) the min / max dose bracket at every probe and lattice point; (Ball) the
closed-form 10 g-ball mean against a tensor cubature; (K-lat) custody of entry 52's lines;
controls (c-u) and (c-ball).  Entry 52's SAR10g shape is printed beside the eddy ball-mean
shape, never gated.

Not claimed: that the clinical SAR is eddy-like, which picture entry 59 follows, any
entry-59 / entry-60 figure or ratio to them, a dose band, a compliance figure, a
non-quasi-static model, a 10 g-averaged θ.  Per deposited gel watt throughout; entry 59 is
per accepted watt (θ_acc = θ·P_gel/P_acc) and is never set beside these without that.
"""

from __future__ import annotations

import io
import re
from pathlib import Path

import numpy as np

from fem_em_solver.thermal.dose import cem43, running_average
from tests.log_records import read_log_text
from fem_em_solver.utils.analytical import (
    robin_box_quadratic_steady_temperature,
    robin_box_quadratic_transient_temperature,
    robin_box_transient_temperature,
)

from tests.validation import test_geo38_f_clinical_box_loaded as c6
from tests.validation import test_mat21_clinical_probe_reader as p79
from tests.validation import test_mat21_wall_uniform_source_probe_points_closed_form as p75
from tests.validation.test_imp4_pec_cylinder_closed_form import B1_BAND
from tests.validation.test_mat13_transient_perfusion_sphere import INTERIOR_BAND
from tests.validation.test_mat21_box_closed_form import (
    A_BOX,
    B_BOX,
    C_BOX,
    HALF,
    SERIES_TERMS,
    T_C,
    T_INF,
    H,
    K,
    Q,
)
from tests.validation.test_mat21_box_closed_form import RHO_C as RHO_C_0
from tests.validation.test_mat21_box_dose_closed_form import (
    AVG_BAR,
    DOSE_BAR,
    HIST_LABELS,
    HIST_PTS,
    T_GRID,
    W_B2,
    _series_hist,
)
from tests.validation.test_mat21_box_dt_duty125_dose_closed_form import _q_check_g, _win_check_g
from tests.validation.test_mat21_box_duty_dose_closed_form import P39P_SERIES_DOSE
from tests.validation.test_mat21_box_probe_dose_closed_form import NEW_LABELS, NEW_PTS
from tests.validation.test_mat21_box_probe_duty_dose_closed_form import P43_SERIES_DOSE
from tests.validation.test_mat21_box_quadratic_closed_form import QUAD_TERMS, UNIT, WEIGHTS

EPS = float(np.finfo(float).eps)
# Probe order from the two landed batches: P0 <- HIST_PTS[2], P1 <- HIST_PTS[0],
# P2 <- HIST_PTS[1], P3-P5 <- NEW_PTS[0:3].  ("H", i) / ("N", i) = batch, column.
MAP = (("H", 2), ("H", 0), ("H", 1), ("N", 0), ("N", 1), ("N", 2))
# Step 10's theta_N (quadratic steady, QUAD_TERMS), PROBE order, K:
# `20261018T111612Z_MAT-21-step10-W2a-points-n2.log.gz:76-81`.
STEP10_THETA_N = (4.276700, 10.890636, 10.890636, 1.822844, 1.822844, 3.127848)
STEP10_TOL = 5e-7
# θ_s(t*) at Q = 2000 W/m^3, PROBE order, `:.6f`:
# `20261020T140156Z_MAT-21-step12-W1-dose-closed-form.log.gz:94,96,98` (P0, P1, P2) and
# `20261020T190033Z_MAT-21-step15-W1-probe-dose-closed-form.log.gz:70,72,74` (P3-P5).
THETA_TSTAR = ("0.063466", "2.647811", "2.647811", "0.048327", "0.048327", "0.100217")
C52_LOG = Path(__file__).resolve().parents[2] / "docs/testing/logs/20261103T080009Z_GEO-38-C6.log"
C52_STAR = "x-105mm dz-50mm"  # `:12307`
# P2 m = 4 FEM's max |θ_h − θ_s| (K) on this source, cited and printed only:
# `20261020T140522Z_MAT-21-step12-W2-dose-fem.log.gz:85-87` (±y, centre) and
# `20261020T190411Z_MAT-21-step15-W2-probe-dose-fem.log.gz:97-99` (±x, P5).
FEM_FLOOR_MP = (2.3933e-05, 1.0716e-03, 1.0716e-03, 6.9684e-05, 6.9684e-05, 1.9878e-04)
# P43's series truncation |hist(QUAD_TERMS) − hist(120,120,60)| (K), printed only:
# `…step15-W1-probe-dose-closed-form.log.gz:85`.
P43_TRUNC = {3: 4.832e-07, 4: 4.832e-07, 5: 7.243e-07}
# Predictions written first (brief), printed only.
PRED_CU = (0.9525, 1.494, 1.494, 0.9545, 0.9545, 0.9065)
PRED_TH_PW = (1.2915e-03, 5.3883e-02, 5.3883e-02, 9.835e-04, 9.835e-04, 2.0394e-03)
PRED_FLOOR_RATIO = (3.8e-4, 4.0e-4, 4.0e-4, 1.44e-3, 1.44e-3, 1.98e-3)
PRED_CEM_P12 = (4.0096e-11, 4.321e-11)
PRED_CSTAR_RATIO = (0.5, 1.0)

RE_LAT = re.compile(r"^\[GEO-38 C6a\] lattice (.+?): \|B1\+\|_ccw (\S+) T, SAR10g (\S+) W/kg, "
                    r"mass (\S+) kg$")
RE_STAR = re.compile(r"^\[GEO-38 C6a\] c\* (.+?): SAR10g (\S+) W/kg \(ccw\)")

_CACHE = {}


def _say(msg):
    print(msg, flush=True)


def _probe_arr():
    return np.asarray(p75.PROBE_PTS, dtype=float)


def _eddy_at(P, t, q=Q, weights=WEIGHTS, n_terms=QUAD_TERMS):
    P = np.asarray(P, dtype=float)
    return np.asarray(robin_box_quadratic_transient_temperature(
        P[:, 0], P[:, 1], P[:, 2], t, half_widths=HALF, k=K, q=q, h=H, rho_c=RHO_C_0,
        weights=weights, t_inf=0.0, n_terms=n_terms), dtype=float)


def _lattice():
    """[(label, box point or None, reason)] in C6's lattice order."""
    off = float(c6._box()[0][1])
    out = []
    for lab, X in c6._lattice():
        pt, why = p79._to_box(X, off)
        out.append((lab, pt, why))
    return out, off


def _data():
    if not _CACHE:
        hist = _series_hist(HIST_PTS)
        new = _series_hist(NEW_PTS)
        src = {"H": hist, "N": new}
        e_abs = np.column_stack([src[b][:, i] for b, i in MAP])
        lat, off = _lattice()
        lat_pts = np.asarray([pt for _lab, pt, _w in lat], dtype=float)
        l_abs = _series_hist(lat_pts)
        u3_lat = np.asarray(robin_box_transient_temperature(
            lat_pts[:, 0], lat_pts[:, 1], lat_pts[:, 2], p75.T_STAR, half_widths=HALF, k=K, q=Q,
            h=H, rho_c=RHO_C_0, t_inf=0.0, n_terms=SERIES_TERMS), dtype=float)
        _CACHE.update(hist=hist, new=new, E_abs=e_abs, E=e_abs / p75.P_DEP,
                      lat=lat, off=off, lat_pts=lat_pts, L_abs=l_abs, L=l_abs / p75.P_DEP,
                      U3=p75._u3(p75.T_STAR) / p75.P_DEP, U3_lat=u3_lat / p75.P_DEP)
    return _CACHE


def _q_src(p):
    """P36's source at points ``p`` (..., 3), W/m^3."""
    p = np.asarray(p, dtype=float)
    w0, wx, wy, wz = WEIGHTS
    return Q * (w0 + wx * (p[..., 0] / HALF[0]) ** 2 + wy * (p[..., 1] / HALF[1]) ** 2
                + wz * (p[..., 2] / HALF[2]) ** 2)


def _ball_mean(p, r, k5=5.0):
    """Closed-form mean of P36's source over the ball (p, r): per axis p_d² + r²/5.
    ``k5 = 3`` is the (c-ball) control (the shell mean)."""
    w0, wx, wy, wz = WEIGHTS
    return Q * (w0 + sum(w * (float(p[d]) ** 2 + r * r / k5) / HALF[d] ** 2
                         for d, w in enumerate((wx, wy, wz))))


def _ball_cub(p, r):
    """Tensor rule on the ball: 4-node GL in s (weight s²), 4-node GL in μ, 8 φ."""
    xs, ws = np.polynomial.legendre.leggauss(4)
    s, w_s = 0.5 * (xs + 1.0), 0.5 * ws * (0.5 * (xs + 1.0)) ** 2
    mu, w_mu = np.polynomial.legendre.leggauss(4)
    phi = 2.0 * np.pi * np.arange(8) / 8.0
    S, MU, PHI = np.meshgrid(s, mu, phi, indexing="ij")
    W = np.einsum("i,j,k->ijk", w_s, w_mu, np.ones(8))
    st = np.sqrt(1.0 - MU**2)
    d = np.stack([st * np.cos(PHI), st * np.sin(PHI), MU], axis=-1)
    pts = np.asarray(p, dtype=float) + r * S[..., None] * d
    return float(np.sum(W * _q_src(pts)) / np.sum(W))


def _parse_c52():
    lat, star = [], []
    # `.log` or `.log.gz`, whichever retention left (§9 item 997); StringIO
    # iteration keeps the file-iteration line set (no str.splitlines()).
    with io.StringIO(read_log_text(C52_LOG, errors="strict")) as fh:
        for line in fh:
            line = line.rstrip("\n")
            m = RE_LAT.match(line)
            if m:
                lat.append((m.group(1), m.group(2), m.group(3), m.group(4)))
            m = RE_STAR.match(line)
            if m:
                star.append((m.group(1), m.group(2)))
    return lat, star


def _cstar_idx():
    return [lab for lab, _pt, _w in _data()["lat"]].index(C52_STAR)


# ------------------------------------------------------------------ tests
def test_preconditions_p_eq():
    gates = []
    d = _data()
    mapped = tuple((HIST_PTS if b == "H" else NEW_PTS)[i] for b, i in MAP)
    same = all(float(a) == float(b) for pm, pp in zip(mapped, p75.PROBE_PTS)
               for a, b in zip(pm, pp)) and len(mapped) == len(p75.PROBE_PTS)
    gates.append((("MAP", mapped, p75.PROBE_PTS), same, 0))
    tg = (p75.T_GRID is T_GRID) or bool(np.array_equal(p75.T_GRID, T_GRID))
    gates.append((("T_GRID",), tg, p75.T_GRID is T_GRID))
    dts = abs(p75.T_STAR - float(T_GRID[-1]))
    gates.append((("T_STAR", p75.T_STAR, float(T_GRID[-1])), dts <= 1e-9, dts))
    wsum = abs(WEIGHTS[2] + WEIGHTS[3] - 3.0)
    gates.append((("weights", WEIGHTS), WEIGHTS[0] == 0.0 and WEIGHTS[1] == 0.0
                   and wsum <= 4 * EPS, wsum))
    p_eddy = Q * 8.0 * A_BOX * B_BOX * C_BOX * (WEIGHTS[0] + (WEIGHTS[1] + WEIGHTS[2]
                                                                + WEIGHTS[3]) / 3.0)
    peq = abs(p_eddy / p75.P_DEP - 1.0)
    gates.append((("P-eq", p_eddy, p75.P_DEP), peq <= 4 * EPS, peq))
    off = d["off"]
    bound = 4.0 * EPS * (abs(off) + C_BOX)
    worst_c, r = 0.0, float(c6._a10())
    worst_wall = -np.inf
    for lab, pt, why in d["lat"]:
        if pt is None:
            gates.append((("lattice", lab, why), False, None))
            continue
        worst_c = max(worst_c, abs(pt[2]))
        clear = max(abs(pt[k]) + r - HALF[k] for k in range(3))
        worst_wall = max(worst_wall, clear)
        gates.append((("ball inside", lab, clear), clear <= 0.0, clear))
    gates.append((("lattice c=0",), worst_c <= bound, worst_c))
    gates.append((("lattice n",), len(d["lat"]) == 18, len(d["lat"])))
    gates.append((("row0",), bool(np.all(d["E_abs"][0] == 0.0) and np.all(d["L_abs"][0] == 0.0)),
                  0))
    _say(f"MAT21-44 MAP -> {mapped} vs p75.PROBE_PTS {p75.PROBE_PTS}: equal {same} (ASSERTED)")
    _say(f"MAT21-44 T_GRID is p75.T_GRID {p75.T_GRID is T_GRID}; p75.T_STAR {p75.T_STAR:.9f} s vs "
         f"T_GRID[-1] {float(T_GRID[-1]):.9f} s, |d| {dts:.3e} s (bound 1e-9, ASSERTED)")
    _say(f"MAT21-44 WEIGHTS {WEIGHTS}; |w_y + w_z - 3| {wsum:.3e} (bound 4 eps {4 * EPS:.3e}, "
         f"ASSERTED)")
    _say(f"MAT21-44 (P-eq) eddy P_dep {p_eddy:.9f} W vs p75.P_DEP {p75.P_DEP:.9f} W, |ratio-1| "
         f"{peq:.3e} (bound 4 eps, ASSERTED)")
    _say(f"MAT21-44 lattice: {len(d['lat'])} points, off {off!r} m, max |c| {worst_c:.3e} m (bound "
         f"{bound:.3e}, ASSERTED); 10 g ball r {r:.6f} m, worst max_d(|p_d| + r - L_d) "
         f"{worst_wall:.6f} m (<= 0 ASSERTED)")
    bad = [g for g in gates if not g[1]]
    assert not bad, bad


def test_q2p_rep_s_lim_pw():
    gates = []
    d = _data()
    P = _probe_arr()
    # (Q2-p) on P75's own PROBE_PTS batch
    qu = _eddy_at(P, p75.T_STAR, weights=UNIT, n_terms=SERIES_TERMS)
    u3 = p75._u3(p75.T_STAR)
    q2p = bool(np.array_equal(qu, u3))
    gates.append((("Q2-p", float(np.max(np.abs(qu - u3)))), q2p, 0))
    _say(f"MAT21-44 (Q2-p) quadratic(UNIT, SERIES_TERMS) vs p75._u3 at t*: np.array_equal {q2p}, "
         f"max |diff| {float(np.max(np.abs(qu - u3))):.3e} K (ASSERTED bitwise)")
    # (Rep-s) on the two landed batches, absolute (Q = 2000 W/m^3)
    hist, new = d["hist"], d["new"]
    got39 = tuple(f"{cem43(T_GRID, T_INF + hist[:, j]):.15e}" for j in range(3))
    got43 = tuple(f"{cem43(T_GRID, T_INF + new[:, j]):.15e}" for j in range(3))
    gates.append((("Rep-s P39'", got39), got39 == tuple(P39P_SERIES_DOSE), 0))
    gates.append((("Rep-s P43", got43), got43 == tuple(P43_SERIES_DOSE), 0))
    th = tuple(f"{v:.6f}" for v in d["E_abs"][-1])
    gates.append((("Rep-s theta(t*)", th), th == THETA_TSTAR, 0))
    _say(f"MAT21-44 (Rep-s) HIST_PTS batch {HIST_LABELS[:3]} doses {got39} vs P39P_SERIES_DOSE: "
         f"{got39 == tuple(P39P_SERIES_DOSE)} (strings ASSERTED)")
    _say(f"MAT21-44 (Rep-s) NEW_PTS batch {NEW_LABELS} doses {got43} vs P43_SERIES_DOSE: "
         f"{got43 == tuple(P43_SERIES_DOSE)} (strings ASSERTED)")
    _say(f"MAT21-44 (Rep-s) theta_s(t*) PROBE order {th} vs {THETA_TSTAR}: {th == THETA_TSTAR} "
         f"(strings ASSERTED)")
    # (Lim)
    ss = np.asarray(robin_box_quadratic_steady_temperature(
        P[:, 0], P[:, 1], P[:, 2], half_widths=HALF, k=K, q=Q, h=H, weights=WEIGHTS, t_inf=0.0,
        n_terms=QUAD_TERMS), dtype=float)
    tl = _eddy_at(P, 1000.0 * T_C)
    lim = float(np.max(np.abs(tl - ss) / np.abs(ss)))
    gates.append((("Lim",), lim <= B1_BAND, lim))
    _say(f"MAT21-44 (Lim) eddy transient(1000 T_C) vs quadratic steady (QUAD_TERMS) max rel "
         f"{lim:.3e} (band {B1_BAND:g}, ASSERTED; P36 read 8.307e-16)")
    for i, lab in enumerate(p75.PROBE_LABELS):
        dv = abs(ss[i] - STEP10_THETA_N[i])
        gates.append(((lab, "step10"), dv <= STEP10_TOL, dv))
        _say(f"MAT21-44 (Lim) {lab}: steady {ss[i]:.9f} K, transient(1000 T_C) {tl[i]:.9f} K, step "
             f"10 theta_N {STEP10_THETA_N[i]:.6f} K, |diff| {dv:.3e} (bound {STEP10_TOL:g}, ASSERTED)")
    # (pw)
    t1 = _eddy_at(P, p75.T_STAR) / p75.P_DEP
    t2 = _eddy_at(P, p75.T_STAR, q=2.0 * Q) / (2.0 * p75.P_DEP)
    pw = float(np.max(np.abs(t2 - t1) / np.abs(t1)))
    gates.append((("pw",), pw <= B1_BAND, pw))
    _say(f"MAT21-44 (pw) theta(2Q)/(2P) vs theta(Q)/P at t*, six points: max rel {pw:.3e} (band "
         f"{B1_BAND:g}, ASSERTED)")
    bad = [g for g in gates if not g[1]]
    assert not bad, bad


def _cem_bracket(th):
    tst = float(T_GRID[-1])
    lo = tst * 0.25 ** (43.0 - T_INF - float(np.min(th)))
    hi = tst * 0.25 ** (43.0 - T_INF - float(np.max(th)))
    return cem43(T_GRID, T_INF + th), lo, hi


def test_mat20_forms_and_cem_bracket():
    gates = []
    d = _data()
    ci = _cstar_idx()
    rows = [(lab, d["E"][:, i]) for i, lab in enumerate(p75.PROBE_LABELS)]
    rows.append((f"c* {C52_STAR}", d["L"][:, ci]))
    worst = {"q": 0.0, "w": 0.0, "rc": np.inf}
    for lab, th in rows:
        _dose, _ref, rel, _swap, miss = _q_check_g(T_INF + th, T_GRID)
        mx, t_mx, lo, hi, _r, wrel = _win_check_g(th, T_GRID)
        gates.append(((lab, "Q"), rel <= DOSE_BAR, rel))
        gates.append(((lab, "swap"), miss > DOSE_BAR, miss))
        gates.append(((lab, "win-loc"),
                      lo <= mx <= hi and T_GRID[0] + W_B2 <= t_mx <= T_GRID[-1], t_mx))
        gates.append(((lab, "Win"), wrel <= AVG_BAR, wrel))
        worst["q"], worst["w"] = max(worst["q"], rel), max(worst["w"], wrel)
        worst["rc"] = min(worst["rc"], miss)
    _say(f"MAT21-44 (Q-d) worst rel {worst['q']:.3e} (bar {DOSE_BAR:g}, ASSERTED); (Win-d) worst "
         f"rel {worst['w']:.3e} (bar {AVG_BAR:g}, ASSERTED); R swap least miss {worst['rc']:.4e} "
         f"(asserted > {DOSE_BAR:g}); per-watt eddy at 6 probe points + c*")
    allrows = [(lab, d["E"][:, i]) for i, lab in enumerate(p75.PROBE_LABELS)]
    allrows += [(lab, d["L"][:, j]) for j, (lab, _p, _w) in enumerate(d["lat"])]
    worst_m = np.inf
    n_tstar = 0
    for lab, th in allrows:
        cem, lo, hi = _cem_bracket(th)
        ok = lo <= cem <= hi
        gates.append(((lab, "Cem", lo, cem, hi), ok, cem))
        worst_m = min(worst_m, min(cem - lo, hi - cem) / cem)
        n_tstar += int(float(np.max(th)) == float(th[-1]))
    _say(f"MAT21-44 (Cem) t* 0.25^(43-T_INF-min) <= cem43 <= t* 0.25^(43-T_INF-max) per watt at "
         f"{len(allrows)} points (6 probe + 18 lattice): all inside "
         f"{all(g[1] for g in gates if g[0][1] == 'Cem')} (ASSERTED); least relative margin "
         f"{worst_m:.3e}; max theta = theta(t*) at {n_tstar} / {len(allrows)} (printed)")
    bad = [g for g in gates if not g[1]]
    assert not bad, bad


def test_ball_and_k_lat():
    gates = []
    d = _data()
    r = float(c6._a10())
    worst = 0.0
    for (lab, _pt, _w), p in zip(d["lat"], d["lat_pts"]):
        bm, bc = _ball_mean(p, r), _ball_cub(p, r)
        e = abs(bc - bm) / bm
        worst = max(worst, e)
        gates.append(((lab, "Ball", bm, bc), e <= B1_BAND, e))
    _say(f"MAT21-44 (Ball) |cub - closed| / closed max over 18 lattice balls {worst:.3e} (band "
         f"{B1_BAND:g}, ASSERTED); r {r:.6f} m; 128-node rule")
    lat, star = _parse_c52()
    labs = {lab for lab, _X in c6._lattice()}
    gates.append((("K-lat n", len(lat)), len(lat) == 18, len(lat)))
    gates.append((("K-lat labels",), {m[0] for m in lat} == labs and len({m[0] for m in lat}) == 18,
                  0))
    fmt_bad = [m for m in lat if f"{float(m[2]):.9e}" != m[2]]
    gates.append((("K-lat fmt", fmt_bad), not fmt_bad, len(fmt_bad)))
    tok = {m[0]: m[2] for m in lat}
    star_ok = (len(star) == 1 and star[0][0] == C52_STAR and star[0][1] == tok.get(C52_STAR))
    gates.append((("K-lat c*", star), star_ok, 0))
    _say(f"MAT21-44 (K-lat) {C52_LOG.name}: {len(lat)} lattice matches (18 ASSERTED), label set == "
         f"c6._lattice() {({m[0] for m in lat} == labs)} (ASSERTED), .9e re-format mismatches "
         f"{len(fmt_bad)} (0 ASSERTED); c* line {star} label {C52_STAR!r}, token == lattice token "
         f"{star_ok} (ASSERTED)")
    bad = [g for g in gates if not g[1]]
    assert not bad, bad


def test_controls_c_u_c_ball():
    d = _data()
    gates = []
    for i, lab in enumerate(p75.PROBE_LABELS):
        rt = float(d["E"][-1, i] / d["U3"][i])
        cu = abs(rt - 1.0)
        gates.append(((lab, "c-u", cu), cu > INTERIOR_BAND, cu))
        _say(f"MAT21-44 (c-u) {lab}: eddy/U3 at t* {rt:.6f}, |ratio-1| {cu:.4f} (ASSERTED > "
             f"{INTERIOR_BAND:g}; predicted {PRED_CU[i]}); separation {cu / INTERIOR_BAND:.1f}x")
    r = float(c6._a10())
    mins = {}
    for (lab, _pt, _w), p in zip(d["lat"], d["lat_pts"]):
        bc = _ball_cub(p, r)
        miss = abs(_ball_mean(p, r, k5=3.0) / bc - 1.0)
        gates.append(((lab, "c-ball", miss), miss > B1_BAND, miss))
        key = "x+0mm" if lab.startswith("x+0mm") else "x+-105mm"
        mins.setdefault(key, []).append(miss)
    pred_105 = (4 * r * r / 15) / (0.105**2 + 2 * r * r / 5)
    _say(f"MAT21-44 (c-ball) shell r^2/3 vs cub miss: x+0mm [{min(mins['x+0mm']):.9f}, "
         f"{max(mins['x+0mm']):.9f}] (predicted exactly 2/3); x+-105mm "
         f"[{min(mins['x+-105mm']):.4e}, {max(mins['x+-105mm']):.4e}] (predicted {pred_105:.4e}, "
         f"range [3e-3, 5e-3]); all > {B1_BAND:g} ASSERTED")
    bad = [g for g in gates if not g[1]]
    assert not bad, bad


def test_record_printed():
    """The record (no band); asserts only that it was formed."""
    d = _data()
    P = _probe_arr()
    e360 = _eddy_at(P, float(W_B2)) / p75.P_DEP
    _say(f"MAT21-44 per deposited gel watt: P_dep {p75.P_DEP:.6f} W (both sources)")
    _say("MAT21-44 PROBE point | eddy/P t* K/W | eddy/P 360 s | U3/P t* | eddy/U3 | cem43/W s | "
         "Cem lo | Cem hi | max 360 s mean/W K | FEM floor/max theta_s | P43 trunc K")
    for i, lab in enumerate(p75.PROBE_LABELS):
        th = d["E"][:, i]
        cem, lo, hi = _cem_bracket(th)
        _m, mx, _t = running_average(T_GRID, th, window=float(W_B2))
        fr = FEM_FLOOR_MP[i] / float(np.max(d["E_abs"][:, i]))
        tr = P43_TRUNC.get(i, float("nan"))
        _say(f"MAT21-44 PROBE {lab} | {th[-1]:.7g} | {e360[i]:.7g} | {d['U3'][i]:.7g} | "
             f"{th[-1] / d['U3'][i]:.7g} | {cem:.7g} | {lo:.7g} | {hi:.7g} | {float(mx):.7g} | "
             f"{fr:.3e} (predicted {PRED_FLOOR_RATIO[i]:.2e}) | {tr:.3e}")
        _say(f"MAT21-44 PROBE {lab} predicted eddy/P t* {PRED_TH_PW[i]:.4e} K/W")
    _say(f"MAT21-44 predicted cem43/W at P1 / P2 inside {PRED_CEM_P12} s")
    labs = [lab for lab, _p, _w in d["lat"]]
    ci = _cstar_idx()
    r = float(c6._a10())
    bms = np.array([_ball_mean(p, r) for p in d["lat_pts"]])
    lat52, _star = _parse_c52()
    sar = {m[0]: float(m[2]) for m in lat52}
    _say("MAT21-44 LATTICE label | box (a, b, c) m | eddy/P t* K/W | U3/P t* | eddy/U3 | "
         "ball-mean/ball-mean(c*) | entry 52 SAR10g/SAR10g(c*) ccw")
    eL, uL = d["L"][-1], d["U3_lat"]
    for j, lab in enumerate(labs):
        p = d["lat_pts"][j]
        _say(f"MAT21-44 LATTICE {lab} | ({p[0]:+.3f}, {p[1]:+.3f}, {p[2]:+.3f}) | {eL[j]:.7g} | "
             f"{uL[j]:.7g} | {eL[j] / uL[j]:.7g} | {bms[j] / bms[ci]:.4g} | "
             f"{sar[lab] / sar[C52_STAR]:.4g}")
    e_max = [labs[j] for j in range(len(labs)) if eL[j] == np.max(eL)]
    e_near = [labs[j] for j in range(len(labs)) if eL[j] >= np.max(eL) * (1.0 - 1e-12)]
    u_max = [labs[j] for j in range(len(labs)) if uL[j] == np.max(uL)]
    u_near = [labs[j] for j in range(len(labs)) if uL[j] >= np.max(uL) * (1.0 - 1e-12)]
    _say(f"MAT21-44 eddy argmax (exact) {e_max}; within 1e-12 rel {e_near} (predicted the four "
         f"x+-105mm dz+-50mm, {C52_STAR!r} among them: {C52_STAR in e_near})")
    _say(f"MAT21-44 U3 argmax (exact) {u_max}; within 1e-12 rel {u_near} (predicted x+0mm "
         f"dz+-50mm)")
    _say(f"MAT21-44 eddy/U3 at c* {eL[ci] / uL[ci]:.6f} (predicted in {PRED_CSTAR_RATIO}; "
         f"point-source basis w_y/4 {WEIGHTS[2] / 4:.6f})")
    pred_shape = (2 * r * r / 5) / (0.105**2 + 2 * r * r / 5)
    for lab in ("x+0mm dz-50mm", "x+0mm dz+50mm", "x-105mm dz-250mm"):
        j = labs.index(lab)
        _say(f"MAT21-44 shape {lab} / c*: eddy ball mean {bms[j] / bms[ci]:.4g} (predicted "
             f"{pred_shape:.4g} at x+0mm, 1 at x-105mm dz-250mm) vs entry 52 SAR10g "
             f"{sar[lab] / sar[C52_STAR]:.4g} (printed only)")
    _say("MAT21-44 quasi-static: shape per deposited watt is frequency-independent, so the same "
         "table serves entry 60 (stated, not asserted)")
    assert np.all(np.isfinite(eL)) and bms[ci] > 0.0
