"""`MAT-21` step 46 — the field-strength heating prediction, closed form (§9 item 935, §10 P83;
brief `docs/planning/chunks/MAT-21.md` "§9 item 935").

Numpy (complex build only because P79 / P80 import P61 / C6, and C33's reader).  Sets no band
and moves none; every band is imported.  Three continuous-source pictures per deposited gel
watt on the unwalled six-face box — uniform U3 (P75 / P79), P36's eddy (P80), step 6's end
source (P81) — plus P75's 1D-reduced walled uniform W_src at z = 0, each read at each side's
c* (C6a's SAR10g argmax on its 18-point lattice: entry 52 at 1.5 T, entry 58 at 3 T) and
scaled per accepted watt by s_side = P_box / P_acc from the side's ``.9e`` normaliser tokens.
The SAR10g columns are read from the two **tracked** C6a logs through the ``ans:11`` page's
``parse_c6_log`` (``c33._load_page()``); C33's strings from its tracked closing log through
``EX-193``'s ``read_c33_rows`` (loaded by path).  No entry-59 / entry-60 log is opened.

Asserted (collected, then asserted — item 668): preconditions incl. (M) the box mass string;
(K) C33's re-read by string; (Sym) the mirror identity; (Sup); (Lim); (Bq); MAT-20's (Q-d) /
(Win-d); (Cem); (Ball); controls (c-u), (c-e), (c-sym).  Printed (🟡 record, no band): s and
the power share, θ(c*, t*) per accepted watt, the dose with its unheated baseline, the 10 g
ball model's shape part, the tie sets, (c-pow).

Not claimed: any entry-59 / entry-60 figure or ratio to them (C39); which picture the clinical
heating follows; an eddy / end-weighted walled box; a dose lever, a dose band, a compliance
figure, a 10 g-averaged θ; that 3.6935 is the clinical heating lever.

One disclosed local read: the (i) line's ``mass … kg`` token, which the page's parser does
not carry, is taken by a plain string split of the one ``[GEO-38 C6a] (i) ball radius`` line
(no regex of the SAR10g columns is re-implemented).
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

from fem_em_solver.thermal.dose import cem43, running_average
from fem_em_solver.utils.analytical import (
    robin_box_quadratic_steady_temperature,
    robin_box_separable_steady_temperature,
    robin_box_steady_temperature,
    robin_box_transient_temperature,
)

from tests.validation import test_geo38_c33_field_strength_reader as c33
from tests.validation import test_geo38_f_clinical_box_loaded as c6
from tests.validation import test_mat21_clinical_probe_reader as p79
from tests.validation import test_mat21_eddy_source_probe_points_closed_form as p80
from tests.validation import test_mat21_end_source_probe_points_closed_form as p81
from tests.validation import test_mat21_wall_uniform_source_probe_points_closed_form as p75
from tests.validation.test_birdcage_sar_mass_averaged import PHANTOM_RHO_KG_PER_M3
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
from tests.validation.test_mat21_box_dose_closed_form import AVG_BAR, DOSE_BAR, T_GRID, W_B2
from tests.validation.test_mat21_box_dt_duty125_dose_closed_form import _q_check_g, _win_check_g
from tests.validation.test_mat21_box_quadratic_closed_form import QUAD_TERMS, WEIGHTS
from tests.validation.test_mat21_box_separable_closed_form import PROFILES, SEP_TERMS

REPO = Path(__file__).resolve().parents[2]
EX193_PATH = REPO / "examples" / "mri" / "41_f_clinical_field_strength_lever.py"
EX193_MOD = "ex193_field_strength_lever"
C33_LOG = REPO / "docs" / "testing" / "logs" / "20261109T103601Z_GEO-38-C33-W1-n1-r2.log"
SIDES = ("1.5T", "3T")
PICS = ("u", "e", "n", "w")
PIC_NAME = {"u": "uniform U3", "e": "eddy (P36)", "n": "end (step 6)", "w": "walled u (z = 0)"}
MASS_TOKEN = "2.457000000e+01"  # `…C6.log:12308`, `…C6-3T.log:12257` (string-compared at run time)
C33_RATIO_TOKENS = ("1.410682e+00", "1.344035e+00", "4.101440e+00", "3.693500e+00")  # :127–130
NON_MIRROR = "x+105mm dz+150mm"
# Predictions written first (brief), printed only; never a band.
PRED_S = {"1.5T": 5.285093e-02, "3T": 1.952049e-01}
PRED_SHARE = 3.6935
PRED_TH_DEP = {"u": 0.0271814, "e": 0.02078249, "n": 0.00356031, "w": 0.02701082}  # K/W per dep. W
PRED_TH_TOL = 1e-4
PRED_D0 = 4.009592657894245e-11
PRED_EXC_U15 = (9.9e-04, 9.0e-04, 1.1e-03)
PRED_EXC_RATIO = (3.698, 3.694, 3.705)
PRED_CEM_RATIO = 1.003
PRED_BALL = {"u": 4.070004e-02, "e": 2.93740e-02, "n": 4.7221e-03}
PRED_MP = {"1.5T": {"u": 1.0573, "e": 1.4649, "n": 9.113},
           "3T": {"u": 1.1741, "e": 1.6268, "n": 10.12}}
PRED_SHAPE = 1.110448
PRED_TIES = {"1.5T": {"x-105mm dz-50mm", "x+105mm dz-50mm", "x-105mm dz+50mm"},
             "3T": {"x+105mm dz+50mm", "x+105mm dz-50mm"}}
PRED_CU = {"n/u": 0.8690, "e/u": 0.2354, "n/e": 0.8287}
PRED_CSYM = 5.73
PRED_CPOW = 0.1104

_CACHE = {}


def _say(msg):
    print(msg, flush=True)


def _ex193():
    if EX193_MOD in sys.modules:
        return sys.modules[EX193_MOD]
    spec = importlib.util.spec_from_file_location(EX193_MOD, EX193_PATH)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[EX193_MOD] = mod
    spec.loader.exec_module(mod)
    return mod


def _mass_token(page, path):
    hits = [ln for ln in page._read_log(path).splitlines() if "[GEO-38 C6a] (i) ball radius" in ln]
    if len(hits) != 1:
        return None, len(hits)
    return hits[0].split("; mass ", 1)[1].split(" kg", 1)[0], 1


def _data():
    if _CACHE:
        return _CACHE
    page = c33._load_page()
    p52 = page.parse_c6_log(page.SOURCE_LOG)
    hit = c33.find_log(c33.ENTRY58_GLOB)
    p58 = page.parse_c6_log(hit[0])
    parsed = {"1.5T": p52, "3T": p58}
    d80, d81 = p80._data(), p81._data()
    labels = [lab for lab, _p, _w in d80["lat"]]
    wall = p75._hists()["Wsrc_raw"][:, 2] / p75.P_DEP
    sd = {}
    for sname, p in parsed.items():
        nm = p["normalisers"]
        s = float(nm["p_box"]) / float(nm["p_acc"])
        lab = p["c_star"]["label"]
        idx = labels.index(lab)
        pt = tuple(float(v) for v in d80["lat_pts"][idx])
        dep = {"u": p79._u3_hist_at(pt), "e": np.array(d80["L"][:, idx]),
               "n": np.array(d81["SL"][:, idx]), "w": np.array(wall)}
        sd[sname] = {"s": s, "label": lab, "idx": idx, "pt": pt, "dep": dep,
                     "acc": {k: v * s for k, v in dep.items()}}
    _CACHE.update(page=page, parsed=parsed, hit=hit, d80=d80, d81=d81, labels=labels, sd=sd,
                  m_box=PHANTOM_RHO_KG_PER_M3 * 8.0 * A_BOX * B_BOX * C_BOX)
    return _CACHE


def _fail(gates):
    return [g for g in gates if not g[1]]


# ------------------------------------------------------------------ tests
def test_preconditions_and_m():
    gates = []
    d = _data()
    dts = abs(p75.T_STAR - float(T_GRID[-1]))
    gates.append((("T_STAR", p75.T_STAR), dts <= 1e-9, dts))
    lab6 = [lab for lab, _X in c6._lattice()]
    gates.append((("lattice labels",), d["labels"] == lab6, 0))
    gates.append((("entry 58 kind", d["hit"]), d["hit"] is not None and d["hit"][1] == "tracked", 0))
    worst_c = 0.0
    for lab, pt, why in d["d80"]["lat"]:
        gates.append((("inside", lab, why), pt is not None, why))
        if pt is not None:
            worst_c = max(worst_c, abs(pt[2]))
    bound = 4.0 * float(np.finfo(float).eps) * (abs(d["d80"]["off"]) + C_BOX)
    gates.append((("c=0", worst_c), worst_c <= bound, worst_c))
    for sname in SIDES:
        gates.append((("c* on lattice", sname), d["sd"][sname]["label"] in d["labels"], 0))
    mstr = f"{d['m_box']:.9e}"
    toks = {}
    for sname in SIDES:
        tok, n = _mass_token(d["page"], d["page"].SOURCE_LOG if sname == "1.5T" else d["hit"][0])
        toks[sname] = tok
        gates.append((("M", sname, tok, n), n == 1 and tok == mstr == MASS_TOKEN, n))
    _say(f"MAT21-46 p75.T_STAR {p75.T_STAR:.6f} s vs T_GRID[-1] {float(T_GRID[-1]):.6f} s, |d| "
         f"{dts:.3e} (bound 1e-9, ASSERTED); lattice labels == c6._lattice() in order "
         f"{d['labels'] == lab6} (ASSERTED); entry 58 {d['hit'][0].name} kind {d['hit'][1]} "
         f"(tracked ASSERTED)")
    _say(f"MAT21-46 lattice {len(d['labels'])} points all inside, max |c| {worst_c:.3e} m (bound "
         f"{bound:.3e}, ASSERTED); c* 1.5 T {d['sd']['1.5T']['label']!r} idx {d['sd']['1.5T']['idx']},"
         f" 3 T {d['sd']['3T']['label']!r} idx {d['sd']['3T']['idx']} (on lattice ASSERTED)")
    _say(f"MAT21-46 (M) M_box = rho 8ABC {d['m_box']:.6e} kg; .9e string == both logs' (i) mass "
         f"token: 1.5 T {toks['1.5T'] == mstr}, 3 T {toks['3T'] == mstr} (ASSERTED)")
    bad = _fail(gates)
    assert not bad, bad


def test_k_c33_reread():
    gates = []
    d = _data()
    ex = _ex193()
    p52, p58 = d["parsed"]["1.5T"], d["parsed"]["3T"]
    for sname, p, node in (("1.5T", p52, c33.ENTRY52_NODE_PF), ("3T", p58, c33.ENTRY58_NODE_PF)):
        res = c33.check_side(p, node, sname)
        for n, ok, _det in res:
            gates.append((("check_side", n), ok, 0))
        _say(f"MAT21-46 (K) c33.check_side {sname} node {node}: {sum(ok for _n, ok, _ in res)} / "
             f"{len(res)} pass (all ASSERTED)")
    rows, (cs15, cs3) = ex.read_c33_rows(d["page"]._read_log(C33_LOG), p52["log"], p58["log"])
    s15, s3 = c33.side(p52), c33.side(p58)
    r = c33.ratios(s15, s3)
    keys = list(r)
    mine = tuple(f"{r[k]:.6e}" for k in keys)
    theirs = tuple(rows[k][2] for k in keys)
    gates.append((("ratios", mine, theirs), mine == theirs == C33_RATIO_TOKENS, 0))
    n_side = 0
    for k in keys:
        fmt = (lambda v: f"{v:.6e}") if k == "lattice_mean_per_sqrt_p_acc" else str
        for j, s_ in ((0, s15), (1, s3)):
            ok = fmt(s_[k]) == rows[k][j]
            n_side += int(ok)
            gates.append((("side string", k, j), ok, 0))
    labs_ok = (s15["c_star"]["label"], s3["c_star"]["label"]) == (cs15, cs3) == \
        (d["sd"]["1.5T"]["label"], d["sd"]["3T"]["label"])
    gates.append((("c* labels",), labs_ok, 0))
    n_fmt = 0
    for p in (p52, p58):
        for row in p["lattice"]:
            n_fmt += int(f"{float(row['sar10g']):.9e}" == row["sar10g"])
    gates.append((("K-lat fmt", n_fmt), n_fmt == 36, n_fmt))
    nm15, nm3 = p52["normalisers"], p58["normalisers"]
    share = d["sd"]["3T"]["s"] / d["sd"]["1.5T"]["s"]
    toks = (nm15["p_box"], nm15["p_acc"], nm3["p_box"], nm3["p_acc"],
            nm15["p_box_over_p_acc"], nm3["p_box_over_p_acc"])
    bnd = ex.rel_bound(toks)
    rel = abs(share - float(C33_RATIO_TOKENS[3])) / float(C33_RATIO_TOKENS[3])
    gates.append((("share", share, rel, bnd), rel <= bnd, rel))
    _say(f"MAT21-46 (K) c33.ratios at .6e {mine} vs read_c33_rows {theirs} vs brief "
         f"{C33_RATIO_TOKENS}: equal {mine == theirs == C33_RATIO_TOKENS} (strings ASSERTED)")
    _say(f"MAT21-46 (K) per-side strings equal {n_side} / 8 (ASSERTED); c* labels {cs15!r} | {cs3!r}"
         f" == side() == parse {labs_ok} (ASSERTED); 36 SAR10g tokens .9e re-format equal {n_fmt} / "
         f"36 (ASSERTED)")
    _say(f"MAT21-46 (K) full-token power share s3/s1.5 {share:.7f} vs C33 {C33_RATIO_TOKENS[3]}: "
         f"rel {rel:.3e} (<= sum of 6 half-units {bnd:.3e}, ASSERTED)")
    jp = c33.HOST_JSON_DIR / c33.C33_JSON_NAME
    if jp.exists():
        j = json.loads(jp.read_text())
        jr = tuple(j["ratios"][k] for k in keys)
        gates.append((("host json ratios", jr), jr == mine, 0))
        _say(f"MAT21-46 (K) host JSON {jp.name} present (custody only): ratios {jr} == rebuilt "
             f"{jr == mine} (ASSERTED)")
    else:
        _say(f"MAT21-46 (K) host JSON {jp} absent (custody only, never a red)")
    bad = _fail(gates)
    assert not bad, bad


def test_sym_sup_lim_bq_ball():
    gates = []
    d = _data()
    sd = d["sd"]
    ts = float(T_GRID[-1])
    worst_sym = {}
    for k in PICS:
        a, b = sd["1.5T"]["dep"][k], sd["3T"]["dep"][k]
        q_t = abs(b[-1] / a[-1] - 1.0)
        q_h = float(np.max(np.abs(b[1:] / a[1:] - 1.0)))
        worst_sym[k] = (q_t, q_h, float(np.max(np.abs(b[1:] - a[1:]))))
        gates.append((("Sym", k, q_t, q_h), q_t <= B1_BAND and q_h <= B1_BAND, max(q_t, q_h)))
        gates.append((("row0", k), a[0] == 0.0 and b[0] == 0.0, 0))
    for k in PICS:
        q_t, q_h, dmax = worst_sym[k]
        _say(f"MAT21-46 (Sym) {PIC_NAME[k]}: |theta(c*3)/theta(c*1.5) - 1| t* {q_t:.3e}, rows 1-40 "
             f"{q_h:.3e}, max |diff| {dmax:.3e} K/W (band {B1_BAND:g}, ASSERTED)")
    # (Sup)
    w_sup = 0.0
    w_pw = 0.0
    for sname in SIDES:
        pt = sd[sname]["pt"]
        P = np.asarray([pt], dtype=float)
        u1 = p79._u3_at(pt, ts)
        u2 = float(np.asarray(robin_box_transient_temperature(
            P[:, 0], P[:, 1], P[:, 2], ts, half_widths=HALF, k=K, q=2.0 * Q, h=H, rho_c=RHO_C_0,
            t_inf=0.0, n_terms=SERIES_TERMS), dtype=float)[0])
        e1, e2 = float(p80._eddy_at(P, ts)[0]), float(p80._eddy_at(P, ts, q=2.0 * Q)[0])
        n1, n2 = float(p81._end_at(P, ts)[0]), float(p81._end_at(P, ts, q=2.0 * Q)[0])
        for one, two in ((u1, u2), (e1, e2), (n1, n2)):
            w_sup = max(w_sup, abs(two - 2.0 * one) / abs(2.0 * one))
            w_pw = max(w_pw, abs(two / (2.0 * p75.P_DEP) - one / p75.P_DEP) / abs(one / p75.P_DEP))
    for t in (ts, float(W_B2)):
        s1, s2 = float(p75._src_q(t, Q)[2]), float(p75._src_q(t, 2.0 * Q)[2])
        w_sup = max(w_sup, abs(s2 - 2.0 * s1) / abs(2.0 * s1))
    gates.append((("Sup",), w_sup <= B1_BAND, w_sup))
    gates.append((("pw",), w_pw <= B1_BAND, w_pw))
    _say(f"MAT21-46 (Sup) theta(2Q) vs 2 theta(Q), u / e / n at both c* (t*) and walled z = 0 (t*, "
         f"360 s): max rel {w_sup:.3e}; per-watt invariance max rel {w_pw:.3e} (band {B1_BAND:g}, "
         f"ASSERTED)")
    # (Lim)
    P2 = np.asarray([sd[s]["pt"] for s in SIDES], dtype=float)
    tl = 1000.0 * T_C
    ss_u = np.asarray(robin_box_steady_temperature(
        P2[:, 0], P2[:, 1], P2[:, 2], half_widths=HALF, k=K, q=Q, h=H, t_inf=0.0,
        n_terms=SERIES_TERMS), dtype=float)
    tr_u = np.array([p79._u3_at(tuple(p), tl) for p in P2])
    ss_e = np.asarray(robin_box_quadratic_steady_temperature(
        P2[:, 0], P2[:, 1], P2[:, 2], half_widths=HALF, k=K, q=Q, h=H, weights=WEIGHTS, t_inf=0.0,
        n_terms=QUAD_TERMS), dtype=float)
    tr_e = p80._eddy_at(P2, tl)
    ss_n = np.asarray(robin_box_separable_steady_temperature(
        P2[:, 0], P2[:, 1], P2[:, 2], half_widths=HALF, k=K, q=Q, h=H, profiles=PROFILES,
        t_inf=0.0, n_terms=SEP_TERMS), dtype=float)
    tr_n = p81._end_at(P2, tl)
    lims = {k: float(np.max(np.abs(tr - ss) / np.abs(ss)))
            for k, tr, ss in (("u", tr_u, ss_u), ("e", tr_e, ss_e), ("n", tr_n, ss_n))}
    for k, v in lims.items():
        gates.append((("Lim", k), v <= B1_BAND, v))
    _say(f"MAT21-46 (Lim) transient(1000 T_C) vs steady at both c*: u {lims['u']:.3e}, e "
         f"{lims['e']:.3e}, n {lims['n']:.3e} (band {B1_BAND:g}, ASSERTED; P80 / P81 read 4.872e-16 "
         f"/ 6.772e-16)")
    # (Bq)
    bq = max(abs(float(p75._src_q(t, Q)[2]) - float(p75._src_q(t, Q, n=2 * p75.N_SERIES)[2]))
             / abs(float(p75._src_q(t, Q, n=2 * p75.N_SERIES)[2])) for t in (ts, float(W_B2)))
    gates.append((("Bq",), bq <= B1_BAND, bq))
    _say(f"MAT21-46 (Bq) walled W_src z = 0, N_SERIES {p75.N_SERIES} vs 2N at t* and 360 s: max rel "
         f"{bq:.3e} (band {B1_BAND:g}, ASSERTED; P75 read 0.000e+00)")
    # (Ball)
    r = float(c6._a10())
    wb = {"e": 0.0, "n": 0.0}
    for sname in SIDES:
        pt = sd[sname]["pt"]
        bm, bc = p80._ball_mean(pt, r), p80._ball_cub(pt, r)
        wb["e"] = max(wb["e"], abs(bc - bm) / bm)
        bm, bc = p81._ball_mean_end(pt, r), p81._ball_cub_end(pt, r)
        wb["n"] = max(wb["n"], abs(bc - bm) / bm)
    for k, v in wb.items():
        gates.append((("Ball", k), v <= B1_BAND, v))
    _say(f"MAT21-46 (Ball) |cub - closed| / closed at both c*: eddy {wb['e']:.3e}, end {wb['n']:.3e} "
         f"(band {B1_BAND:g}, ASSERTED; P80 3.150e-16, P81 3.944e-14); r {r:.6f} m")
    bad = _fail(gates)
    assert not bad, bad


def test_mat20_forms_and_cem():
    gates = []
    sd = _data()["sd"]
    worst = {"q": 0.0, "w": 0.0, "rc": np.inf, "m": np.inf}
    n = 0
    for sname in SIDES:
        for k in PICS:
            th = sd[sname]["acc"][k]
            _dose, _ref, rel, _swap, miss = _q_check_g(T_INF + th, T_GRID)
            mx, t_mx, lo, hi, _r, wrel = _win_check_g(th, T_GRID)
            gates.append(((sname, k, "Q"), rel <= DOSE_BAR, rel))
            gates.append(((sname, k, "swap"), miss > DOSE_BAR, miss))
            gates.append(((sname, k, "win-loc"),
                          lo <= mx <= hi and T_GRID[0] + W_B2 <= t_mx <= T_GRID[-1], t_mx))
            gates.append(((sname, k, "Win"), wrel <= AVG_BAR, wrel))
            cem, clo, chi = p80._cem_bracket(th)
            gates.append(((sname, k, "Cem", clo, cem, chi), clo <= cem <= chi, cem))
            worst["q"], worst["w"] = max(worst["q"], rel), max(worst["w"], wrel)
            worst["rc"] = min(worst["rc"], miss)
            worst["m"] = min(worst["m"], min(cem - clo, chi - cem) / cem)
            n += 1
    _say(f"MAT21-46 (Q-d) worst rel {worst['q']:.3e} (bar {DOSE_BAR:g}, ASSERTED); R swap least miss "
         f"{worst['rc']:.4e} (asserted > {DOSE_BAR:g}); (Win-d) worst rel {worst['w']:.3e} (bar "
         f"{AVG_BAR:g}, ASSERTED), window located (ASSERTED); {n} per-accepted-watt histories")
    _say(f"MAT21-46 (Cem) p80._cem_bracket on {n} histories: all inside "
         f"{all(g[1] for g in gates if g[0][2] == 'Cem')} (ASSERTED); least relative margin "
         f"{worst['m']:.3e}")
    bad = _fail(gates)
    assert not bad, bad


def test_controls():
    gates = []
    d = _data()
    sd = d["sd"]
    for sname in SIDES:
        dep = sd[sname]["dep"]
        u, e, n = dep["u"][-1], dep["e"][-1], dep["n"][-1]
        for tag, v in (("n/u", n / u), ("e/u", e / u), ("n/e", n / e)):
            c = abs(v - 1.0)
            gates.append(((sname, tag, c), c > INTERIOR_BAND, c))
            _say(f"MAT21-46 (c-{'e' if tag == 'n/e' else 'u'}) {sname} c* {sd[sname]['label']}: "
                 f"{tag} at t* {v:.6f}, |ratio-1| {c:.4f} (ASSERTED > {INTERIOR_BAND:g}; predicted "
                 f"{PRED_CU[tag]}); separation {c / INTERIOR_BAND:.1f}x")
    j = d["labels"].index(NON_MIRROR)
    i15 = sd["1.5T"]["idx"]
    sl = d["d81"]["SL"]
    cs = abs(sl[-1, j] / sl[-1, i15] - 1.0)
    gates.append((("c-sym", cs), cs > B1_BAND, cs))
    _say(f"MAT21-46 (c-sym) end picture, c*3 replaced by {NON_MIRROR!r}: |quotient - 1| {cs:.4f} "
         f"(ASSERTED > {B1_BAND:g}; predicted approx {PRED_CSYM})")
    bad = _fail(gates)
    assert not bad, bad


def _in(v, lo, hi):
    return lo <= v <= hi


def test_record_printed():
    """The record (no band); asserts only that it was formed."""
    d = _data()
    sd = d["sd"]
    nm = {s: d["parsed"][s]["normalisers"] for s in SIDES}
    share = sd["3T"]["s"] / sd["1.5T"]["s"]
    _say(f"MAT21-46 normaliser: per accepted watt = per deposited gel watt x P_box/P_acc; P_dep "
         f"{p75.P_DEP:.2f} W")
    for s in SIDES:
        _say(f"MAT21-46 s {s} = P_box/P_acc from the .9e tokens {sd[s]['s']:.7g} (predicted "
             f"{PRED_S[s]:.6e}); C33's 7-digit token {nm[s]['p_box_over_p_acc']}")
    _say(f"MAT21-46 power share s3/s1.5 {share:.7g} (predicted approx {PRED_SHARE})")
    _say("MAT21-46 THETA picture | 1.5 T K/W per acc. W | 3 T | 3T/1.5T | per dep. W (both) | "
         "predicted per dep. W, within 1e-4")
    for k in PICS:
        a, b = sd["1.5T"]["acc"][k][-1], sd["3T"]["acc"][k][-1]
        dep = sd["1.5T"]["dep"][k][-1]
        ok = abs(dep / PRED_TH_DEP[k] - 1.0) <= PRED_TH_TOL
        _say(f"MAT21-46 THETA {PIC_NAME[k]} | {a:.6e} | {b:.6e} | {b / a:.7g} | {dep:.7g} | "
             f"{PRED_TH_DEP[k]} inside {ok}")
    _say(f"MAT21-46 the pictures' 3T/1.5T theta lever = power share x (Sym) 1 = {share:.7g} for "
         f"every picture (the power-share part)")
    # dose
    d0 = cem43(T_GRID, np.full(T_GRID.size, T_INF))
    _say(f"MAT21-46 D0 unheated 20 C dose {d0:.7g} s (predicted {PRED_D0:.7g}, |rel| "
         f"{abs(d0 / PRED_D0 - 1):.2e})")
    _say("MAT21-46 DOSE side picture | cem43 per acc. W s | D0 s | excess D-D0 s | excess/D0 | "
         "max 360 s mean K | at t s | at t*")
    exc = {}
    cemv = {}
    for s in SIDES:
        for k in PICS:
            th = sd[s]["acc"][k]
            cem = cem43(T_GRID, T_INF + th)
            _m, mx, tmx = running_average(T_GRID, th, window=float(W_B2))
            exc[(s, k)], cemv[(s, k)] = cem - d0, cem
            _say(f"MAT21-46 DOSE {s} {PIC_NAME[k]} | {cem:.7g} | {d0:.7g} | {cem - d0:.6e} | "
                 f"{(cem - d0) / d0:.6e} | {float(mx):.6e} | {float(tmx):.6g} | "
                 f"{abs(float(tmx) - float(T_GRID[-1])) <= 1e-9}")
    e15 = exc[("1.5T", "u")] / d0
    pv, lo, hi = PRED_EXC_U15
    _say(f"MAT21-46 excess/D0 under u at 1.5 T {e15:.4e} (predicted {pv}, range [{lo}, {hi}] inside "
         f"{_in(e15, lo, hi)})")
    pv, lo, hi = PRED_EXC_RATIO
    for k in PICS:
        er = exc[("3T", k)] / exc[("1.5T", k)]
        cr = cemv[("3T", k)] / cemv[("1.5T", k)]
        _say(f"MAT21-46 DOSE ratio 3T/1.5T {PIC_NAME[k]}: excess {er:.5g} (predicted {pv}, range "
             f"[{lo}, {hi}] inside {_in(er, lo, hi)}); cem43 {cr:.5g} (predicted approx "
             f"{PRED_CEM_RATIO}; printed only -- a cem43 per accepted watt is the 1 W history's dose, "
             f"approx D0 plus an excess, not a scalable per-watt dose and never a lever)")
    # ball model
    r = float(c6._a10())
    m_box = d["m_box"]
    mp = {}
    for s in SIDES:
        pt = sd[s]["pt"]
        pic = {"u": 1.0 / m_box, "e": p80._ball_mean(pt, r) / (Q * m_box),
               "n": p81._ball_mean_end(pt, r) / (Q * m_box)}
        meas = float(nm[s]["sar_per_box_w"])
        meas_acc = float(nm[s]["sar_per_p_acc"])
        _say(f"MAT21-46 BALL {s} c* {sd[s]['label']}: measured SAR10g(c*) per box W {meas:.6e}, per "
             f"acc. W {meas_acc:.6e} W/kg/W")
        mp[s] = {}
        for k in ("u", "e", "n"):
            mp[s][k] = meas / pic[k]
            _say(f"MAT21-46 BALL {s} {PIC_NAME[k]}: picture per box W {pic[k]:.6e} (predicted "
                 f"{PRED_BALL[k]}), per acc. W {pic[k] * sd[s]['s']:.6e}; measured/picture per box W "
                 f"{mp[s][k]:.5g} (predicted {PRED_MP[s][k]}), per acc. W "
                 f"{meas_acc / (pic[k] * sd[s]['s']):.5g}")
    for k in ("u", "e", "n"):
        _say(f"MAT21-46 BALL shape part 3T/1.5T {PIC_NAME[k]}: {mp['3T'][k] / mp['1.5T'][k]:.7g} "
             f"(predicted {PRED_SHAPE})")
    # tie sets
    d80, d81 = d["d80"], d["d81"]
    for s in SIDES:
        lat = d["parsed"][s]["lattice"]
        sar = np.array([float(x["sar10g"]) for x in lat])
        ties = [x["label"] for x, v in zip(lat, sar) if v >= 0.99 * sar.max()]
        _say(f"MAT21-46 TIES {s} within 1 % of max SAR10g: {ties} (predicted {sorted(PRED_TIES[s])}"
             f": {set(ties) == PRED_TIES[s]})")
        vals = {k: [] for k in ("u", "e", "n")}
        for lab in ties:
            j = d["labels"].index(lab)
            th = {"u": d80["U3_lat"][j], "e": d80["L"][-1, j], "n": d81["SL"][-1, j]}
            for k in vals:
                vals[k].append(th[k] * sd[s]["s"])
            _say(f"MAT21-46 TIES {s} {lab}: theta(t*) per acc. W u {th['u'] * sd[s]['s']:.6e}, e "
                 f"{th['e'] * sd[s]['s']:.6e}, n {th['n'] * sd[s]['s']:.6e} K/W")
        _say(f"MAT21-46 TIES {s} relative spread over the set: "
             + ", ".join(f"{k} {(max(v) - min(v)) / max(v):.2e}" for k, v in vals.items()))
    cpow = abs(float(C33_RATIO_TOKENS[2]) / float(C33_RATIO_TOKENS[3]) - 1.0)
    _say(f"MAT21-46 (c-pow) |{C33_RATIO_TOKENS[2]} / {C33_RATIO_TOKENS[3]} - 1| {cpow:.4f} "
         f"(predicted {PRED_CPOW}): the SAR lever is not the power share alone (printed)")
    assert np.isfinite(share) and d0 > 0.0
