"""`MAT-21` step 47 — P86: P76's record on entry 60 (3 T).  Entry 60's six probe θ and
hotspot θ / cem43 / six-minute window per accepted watt beside the three closed-form
pictures per deposited watt (§9 item 955, §10 P86; brief `docs/planning/chunks/MAT-21.md`
"§9 item 955").

Complex build (P79 imports P61 / C6), no solve.  Entry 60 is located by glob
(``c46.find_log(c46.ENTRY60_GLOB)``, no knob) and read through **P79's own**
``_parse_p61_text`` on the lines containing ``[GEO-38 P61] `` only: entry 60's cw block is
cut by interleaved pytest node-id lines, which P79's consecutive-run block rule does not
bridge (control (c-knob) shows it).  The filter is the only new code on the read path;
(K-filter) shows it is the identity on entry 59 and (K-C46) ties every token and line
number to C46's independent reader on entry 60.

The pictures (P75's U3, P80's eddy, P81's end-weighted) are per deposited gel watt and
depend only on the box and the gel's thermal constants; entry 60 heated with the same
(1.5 T) gel constants (C31, entry 60 `:139`).  The only 3 T input is s = P_gel/P_acc, read
from entry 60's own tokens.  No 3 T gel claim.

Asserted (collected, then asserted; print-precision custody, C33's (D-*) pattern, no new
band): (Real), (Gel), (K-60), (K-C46), (K-filter), (K-s), (K-norm), (R58), (K-twin),
(P83-tie), (Id-3); controls (c-knob), (c-mirror).  Printed, no band (a 🟡 record): the
per-sense table, the hotspot dose against D0, the six-minute window.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
from mpi4py import MPI

from fem_em_solver.thermal.dose import cem43, running_average

from tests.validation import test_geo38_c33_field_strength_reader as c33
from tests.validation import test_geo38_c39_thermal_reader as c46
from tests.validation import test_geo38_f_clinical_box_loaded as c6
from tests.validation import test_mat21_clinical_probe_reader as p79
from tests.validation import test_mat21_p76_three_pictures as p76
from tests.validation.test_geo36_rod_f_clinical_heating import H_GEL, K_GEL
from tests.validation.test_geo36_rod_f_clinical_heating import RHO_C as RHO_C_GEL
from tests.validation.test_imp4_pec_cylinder_closed_form import B1_BAND
from tests.validation.test_mat21_box_closed_form import RHO_C, T_INF, H, K
from tests.validation.test_mat21_box_dose_closed_form import AVG_BAR, DOSE_BAR, T_GRID, W_B2

PFX = "MAT21-47"
PROBE_LABELS, P_DEP = p76.PROBE_LABELS, p76.P_DEP
JSON_PATH = p79.REPO / "logs" / "mat21" / "p86-three-pictures-3t.json"
FILTER = "[GEO-38 P61] "
WANT_BASE = {"e60": "20261111T080008Z_GEO-38-P33-3T.log", "e59": "20261110T080008Z_GEO-38-P33.log",
             "e58": "20261109T080008Z_GEO-38-C6-3T.log"}
WANT_LINES60 = {"ccw": [list(range(12275, 12283))],
                "cw": [[12293, 12294] + list(range(12298, 12304))]}
WANT_LINES59 = {"ccw": [list(range(12325, 12333))], "cw": [list(range(12343, 12351))]}
WANT_STAR = "x+105mm dz+50mm"
WANT_STAR_IDX = 15  # P83's `…P83-W1.log:45`
MIRROR_STAR = "x-105mm dz-50mm"  # entry 59's c*
N_VALUES = 28
# (Gel) the I32 twin's geometry tokens and C31's gel line (entry 60 `:139`), by string.
V_SUB_TOK, ROBIN_TOK = "2.457000000000e-02", "7.386000000e-01"
C31_LINE_NO = 139
C31_LINE = ("[GEO-38 C31] the gel sigma / eps_r and the thermal constants are the 1.5 T values "
            "(provisional until W3; no 3 T gel claim)")
RE_TWIN_GEO = re.compile(r"\[GEO-36 I32\] \[gel twin (cw|ccw)\] sense .*?; V_sub (\S+) m\^3; "
                         r"Robin area (\S+) m\^2;")
# (P83-tie) P83's printed 3 T theta(c*, t*) per accepted watt (`20261110T103343Z_MAT-21-P83-W1.log:79–81`).
P83_3T = {"U3": "5.305942e-03", "eddy": "4.056844e-03", "end": "6.949899e-04"}
# (Sym) band of P83 (`…P83-W1.log:54–56`, read 0) reused for (c-mirror); not a new band.
SYM_BAND = 1e-12
# Brief's expected tables (printed beside; a miss > 1e-4 rel is disclosed, never red).
EXPECT_CCW = {
    "P0 centre": (0.006751649, 0.2484, 5.228, 13.38),
    "P1 +y wall": (0.06519120, 3.017, 1.210, 170.1),
    "P2 -y wall": (0.06216465, 2.877, 1.154, 162.2),
    "P3 +x wall": (0.01505990, 0.6971, 15.31, 0.3519),
    "P4 -x wall": (0.01581514, 0.7320, 16.08, 0.3695),
    "P5 top face": (0.01211746, 0.5553, 5.942, 31.11),
    "hotspot": (0.03338110, 1.228, 1.606, 9.376),
}
EXPECT_CW = {"P0 centre": (0.3106, 6.536, 16.73), "P1 +y wall": (2.701, 1.083, 152.3),
             "P5 top face": (0.7153, 7.654, 40.08), "hotspot": (1.184, 1.549, 9.043)}
EXPECT_EXCESS = {"ccw": (1.79681e-13, 4.4813e-03), "cw": (1.56517e-13, 3.9036e-03)}
EXPECT_PIC_EXCESS_CCW = {"U3": 1.4878e-13, "eddy": 1.1154e-13, "end": 1.8641e-14}
EXPECT_PIC_WINDOW_CCW = {"U3": 4.9759e-03, "eddy": 3.7923e-03, "end": 6.4679e-04}

_CACHE: dict = {}


def _say(msg):
    if MPI.COMM_WORLD.rank == 0:
        print(msg, flush=True)


def _real_path(base: Path) -> Path:
    """The file on disk (the base name, or its ``.gz`` once housekeeping gzipped it)."""
    return base if base.is_file() else Path(str(base) + ".gz")


def _logs():
    if "logs" not in _CACHE:
        _CACHE["logs"] = {"e60": c46.find_log(c46.ENTRY60_GLOB),
                          "e59": c46.find_log(c46.ENTRY59_GLOB),
                          "e58": c33.find_log(c33.ENTRY58_GLOB)}
    return _CACHE["logs"]


def _p61_only(path):
    """Joined text of the lines containing ``[GEO-38 P61] `` and their 1-based log line
    numbers (filtered line i ↔ ``nums[i - 1]``)."""
    keep, nums = [], []
    for n, raw in enumerate(c46._read_log(path).splitlines(), start=1):
        if FILTER in raw:
            keep.append(raw)
            nums.append(n)
    return "\n".join(keep), nums


def _parse_filtered(path):
    """P79's own parser on the filtered text, block lines mapped back to the log's."""
    text, nums = _p61_only(path)
    parsed = p79._parse_p61_text(text)
    for rec in parsed.values():
        rec["_lines"] = [[nums[i - 1] for i in blk] for blk in rec["_lines"]]
    if not parsed:
        raise ValueError(f"{path}: no P61 sense block found")
    for s, rec in parsed.items():
        missing = [lab for lab in PROBE_LABELS if lab not in rec["probes"]]
        if missing:
            raise ValueError(f"{path}: sense {s} lacks probes {missing}")
    return parsed


def _parse60():
    if "p60" not in _CACHE:
        path = _logs()["e60"][0]
        _CACHE["p60"] = (str(path), c46._on_root(_parse_filtered, path))
    return _CACHE["p60"]


def _gates_assert(gates):
    bad = [g for g in gates if not g[1]]
    assert not bad, bad


# ------------------------------------------------------------------ tests
def test_real_and_gel():
    gates = []
    logs = _logs()
    for key, want in WANT_BASE.items():
        hit = logs[key]
        ok = hit is not None and hit[1] == "tracked" and hit[0].name == want
        gates.append((("Real", key, hit), ok, 0))
        _say(f"{PFX} (Real) {key}: {hit[0].name if hit else None} kind {hit[1] if hit else None} "
             f"(want {want}, tracked): {ok} (ASSERTED)")
    # (Gel) constants, exactly
    for nm, a, b in (("H", H, H_GEL), ("K", K, K_GEL), ("RHO_C", RHO_C, RHO_C_GEL)):
        gates.append((("Gel const", nm, a, b), a == b, 0))
        _say(f"{PFX} (Gel) box {nm} {a!r} == I32 gel twin {b!r}: {a == b} (exact, ASSERTED)")
    # (Gel) entry 60's I32 twin lines vs entry 59's, by string
    h60 = c46._on_root(c46.parse_heat_log, logs["e60"][0])
    h59 = c46._on_root(c46.parse_heat_log, logs["e59"][0])
    txt60 = c46._on_root(c46._read_log, logs["e60"][0]).splitlines()
    txt59 = c46._on_root(c46._read_log, logs["e59"][0]).splitlines()
    for s in c46.SENSES:
        for tag, h, txt in (("e60", h60, txt60), ("e59", h59, txt59)):
            ln = h["senses"][s]["i32"]["src"]["line"]
            m = RE_TWIN_GEO.search(txt[ln - 1])
            got = m.groups() if m else None
            ok = got is not None and got == (s, V_SUB_TOK, ROBIN_TOK)
            gates.append((("Gel twin geo", tag, s, ln, got), ok, 0))
            _say(f"{PFX} (Gel) {tag} {s} I32 twin :{ln}: V_sub {got[1] if got else None} m^3 (want "
                 f"{V_SUB_TOK}), Robin area {got[2] if got else None} m^2 (want {ROBIN_TOK}): {ok} "
                 f"(strings, ASSERTED)")
        p60 = {lab: v["pt"] for lab, v in h60["senses"][s]["i32"]["probes"].items()}
        p59 = {lab: v["pt"] for lab, v in h59["senses"][s]["i32"]["probes"].items()}
        ok = p60 == p59 and sorted(p60) == sorted(PROBE_LABELS)
        gates.append((("Gel probe pts", s, p60, p59), ok, 0))
        _say(f"{PFX} (Gel) {s} I32 twin probe coordinate strings e60 == e59 over "
             f"{len(p60)} labels (P1 +y wall {p60.get('P1 +y wall')!r}): {ok} (strings, ASSERTED)")
    c31 = txt60[C31_LINE_NO - 1] if len(txt60) >= C31_LINE_NO else None
    ok = c31 == C31_LINE
    gates.append((("Gel C31 line", c31), ok, 0))
    _say(f"{PFX} (Gel) entry 60 :{C31_LINE_NO} == C31's 1.5 T gel line: {ok} (string, ASSERTED)")
    _gates_assert(gates)


def test_k60_kc46_kfilter():
    gates = []
    path, parsed = _parse60()
    _say(f"{PFX} entry 60 {path!r}: parsed through p79._parse_p61_text on the {FILTER!r} lines")
    # (K-60)
    bad, n = p79._reformat_misses(parsed)
    gates.append((("K-60 reformat", bad, n), not bad and n == N_VALUES, n))
    lines = {s: r["_lines"] for s, r in parsed.items()}
    gates.append((("K-60 lines", lines), lines == WANT_LINES60, 0))
    hr = p79._hot_row(parsed)  # raises if a printed c* is > 5e-7 m off its lattice point
    lat = c6._lattice()
    labs = [lab for lab, _x in lat]
    lab_c, lab_w = parsed["ccw"]["hot"]["label"], parsed["cw"]["hot"]["label"]
    dpt = max(abs(a - b) for a, b in zip(dict(lat).get(lab_c, (np.inf,) * 3),
                                         parsed["ccw"]["hot"]["point"]))
    ok_star = (lab_c == WANT_STAR and lab_w == WANT_STAR and labs.index(WANT_STAR) == WANT_STAR_IDX
               and hr["ccw"]["row"] is not None and hr["cw"]["row"] is not None and dpt <= 5e-7)
    gates.append((("K-60 c*", lab_c, lab_w, dpt), ok_star, dpt))
    _say(f"{PFX} (K-60) {n} values re-formatted (want {N_VALUES}), misses {bad or 'none'}; block "
         f"lines {lines} (want {WANT_LINES60}): {lines == WANT_LINES60}; c* ccw {lab_c!r} cw {lab_w!r}"
         f" (want {WANT_STAR!r}, lattice index {labs.index(lab_c) if lab_c in labs else None}, want "
         f"{WANT_STAR_IDX}), printed vs lattice {dpt:.3e} m (<= 5e-7) (ASSERTED)")
    # (K-C46)
    h60 = c46._on_root(c46.parse_heat_log, path)
    for s, rec in parsed.items():
        c = h60["senses"][s]["p61"]
        blk = rec["_lines"][0]
        pairs = [("P_gel", rec["P_gel"]["tok"], c["src"]["p_gel"], blk[0], c["src"]["line"]),
                 ("P_acc", rec["P_acc"]["tok"], c["src"]["p_acc"], blk[0], c["src"]["line"]),
                 ("factor", rec["factor"]["tok"], c["src"]["factor"], blk[0], c["src"]["line"])]
        h = rec["hot"]
        hl = c["hot"]["line"]
        pairs += [("hot label", "c* " + h["label"], c["hot"]["label"], blk[1], hl),
                  ("hot point", "(" + h["point_tok"] + ")", c["hot"]["point"], blk[1], hl),
                  ("theta", h["th"]["tok"], c["hot"]["theta_point"], blk[1], hl),
                  ("final", h["final"]["tok"], c["hot"]["theta_final"], blk[1], hl),
                  ("cem43", h["cem43"]["tok"], c["hot"]["cem43"], blk[1], hl),
                  ("ra_max", h["ra_max"]["tok"], c["hot"]["ra_max"], blk[1], hl),
                  ("ra_t", h["ra_tmax"]["tok"], c["hot"]["ra_t"], blk[1], hl)]
        cp = c["probes"]
        for (lab, v), ln in zip(rec["probes"].items(), blk[2:]):
            cl = cp.get(lab, {})
            pairs.append((f"probe {lab}", v["tok"], cl.get("theta"), ln, cl.get("line")))
        miss = [p for p in pairs if p[1] != p[2] or p[3] != p[4]]
        ok = not miss and len(pairs) == 16 and sorted(cp) == sorted(PROBE_LABELS)
        gates.append((("K-C46", s, miss), ok, len(miss)))
        _say(f"{PFX} (K-C46) {s}: {len(pairs)} tokens (P_gel, P_acc, factor; hotspot label, point, "
             f"theta, final, cem43, ra max, ra t; six probes) and their line numbers vs "
             f"c46.parse_heat_log({h60['log']}): mismatches {miss or 'none'} (strings, ASSERTED)")
    # (K-filter) on entry 59
    p59 = _logs()["e59"][0]
    filt = c46._on_root(_parse_filtered, p59)
    unf = c46._on_root(p79._parse_p61, _real_path(p59))
    diff = sorted(set(filt) ^ set(unf)) + [s for s in filt if s in unf and filt[s] != unf[s]]
    l59 = {s: r["_lines"] for s, r in filt.items()}
    ok = not diff and l59 == WANT_LINES59 and filt == unf
    gates.append((("K-filter", diff, l59), ok, len(diff)))
    _say(f"{PFX} (K-filter) entry 59 {p59.name}: filtered parse == p79._parse_p61 on every token, "
         f"raw block text and block lines; differing senses {diff or 'none'}; lines {l59} (want "
         f"{WANT_LINES59}): {ok} (ASSERTED)")
    _gates_assert(gates)


def test_ks_knorm_r58_ktwin():
    gates = []
    path, parsed = _parse60()
    # (K-s)
    for s, rec in parsed.items():
        pg, pa, f = rec["P_gel"], rec["P_acc"], rec["factor"]
        miss = abs(pg["val"] / pa["val"] / f["val"] - 1.0)
        bnd = sum(c33._half_ulp_rel(t["tok"]) for t in (pg, pa, f))
        gates.append((("K-s", s, miss, bnd), miss <= bnd, miss))
        _say(f"{PFX} (K-s) {s}: P_gel {pg['tok']} / P_acc {pa['tok']} / factor {f['tok']} - 1 = "
             f"{miss:.3e} (half-ulp sum {bnd:.3e}, ASSERTED)")
    # (K-norm) and (R58)
    page = c33._load_page()
    e60 = c33._on_root(page.parse_c6_log, path)
    e58 = c33._on_root(page.parse_c6_log, _logs()["e58"][0])
    nm = e60["normalisers"]
    for key, tok_key in (("P_gel", "p_box"), ("P_acc", "p_acc")):
        mine = format(parsed["ccw"][key]["val"], ".9e")
        ok = mine == nm[tok_key]
        gates.append((("K-norm", key, mine, nm[tok_key]), ok, 0))
        _say(f"{PFX} (K-norm) ccw {key} {parsed['ccw'][key]['tok']} as .9e {mine} vs entry 60's C6a "
             f"normalisers {tok_key} {nm[tok_key]} (:{nm['line']}): {ok} (ASSERTED)")
    f60, f58 = p76._r52_fields(e60), p76._r52_fields(e58)
    diff = [k for k in f58 if f60[k] != f58[k]]
    ok = not diff and len(f60["lattice"]) == 18 and f60["c_star"]["label"] == WANT_STAR
    gates.append((("R58", diff), ok, diff))
    _say(f"{PFX} (R58) parse_c6_log({e60['log']}) vs entry 58 {e58['log']}: lattice "
         f"{len(f60['lattice'])} rows (label, b1p_ccw, sar10g, mass_kg), c* {f60['c_star']['label']!r}, "
         f"normalisers, S2: differing fields {diff or 'none'} (strings, ASSERTED)")
    # (K-twin)
    twin, gp = c33._on_root(p76._parse_twin, path)
    worst = {}
    for s, rec in parsed.items():
        g = gp.get(s, [])
        gp_ok = len(g) == 1 and g[0][2] == "1.000000000000e+00" and g[0][1] == rec["P_gel"]["tok"]
        gates.append((("K-twin G*P", s, g), gp_ok, 0))
        _say(f"{PFX} (K-twin) {s}: twin G*P token {g[0][2] if g else None} (want 1.000000000000e+00),"
             f" twin P_gel {g[0][1] if g else None} vs P61 {rec['P_gel']['tok']}: {gp_ok} "
             f"(strings, ASSERTED)")
        f = rec["factor"]
        tw = twin.get(s, {})
        worst[s] = (0.0, None)
        for lab in PROBE_LABELS:
            if lab not in tw:
                gates.append((("K-twin missing", s, lab), False, None))
                continue
            th = rec["probes"][lab]
            miss = abs(th["val"] / (f["val"] * float(tw[lab])) - 1.0)
            bnd = sum(c33._half_ulp_rel(t) for t in (th["tok"], f["tok"], tw[lab]))
            if miss / bnd > worst[s][0]:
                worst[s] = (miss / bnd, lab)
            gates.append((("K-twin", s, lab, miss, bnd), miss <= bnd, miss))
            _say(f"{PFX} (K-twin) {s} {lab}: theta_acc {th['tok']} / ({f['tok']} x twin {tw[lab]}) "
                 f"- 1 = {miss:.3e} (half-ulp sum {bnd:.3e}, ASSERTED)")
    _say(f"{PFX} (K-twin) worst miss / bound: ccw {worst['ccw'][0]:.3f} ({worst['ccw'][1]}), cw "
         f"{worst['cw'][0]:.3f} ({worst['cw'][1]}) (predicted 0.750 ccw P0 / 0.580 cw)")
    _gates_assert(gates)


def test_p83_tie_and_id3():
    gates = []
    path, parsed = _parse60()
    pic = p76._pictures()
    t = p76._rows(parsed, pic)
    page = c33._load_page()
    nm = c33._on_root(page.parse_c6_log, path)["normalisers"]
    h = t["ccw"]["hot"]
    fac_tok = parsed["ccw"]["factor"]["tok"]
    norm_b = c33._half_ulp_rel(nm["p_box"]) + c33._half_ulp_rel(nm["p_acc"])
    for k, tok in P83_3T.items():
        got = h["s"] * h["pictures"][k]["per_W"]
        miss = abs(got / float(tok) - 1.0)
        bnd = c33._half_ulp_rel(tok) + c33._half_ulp_rel(fac_tok) + norm_b
        gates.append((("P83-tie", k, miss, bnd), miss <= bnd, miss))
        _say(f"{PFX} (P83-tie) ccw c* {h['label']}: s {fac_tok} x {k}/P {h['pictures'][k]['per_W']:.7g}"
             f" = {got:.7g} vs P83's 3 T {tok} K per acc. W: rel {miss:.3e} (bound {bnd:.3e} = "
             f"half-ulps of the P83 token, the factor and the two .9e normalisers {nm['p_box']} / "
             f"{nm['p_acc']}, ASSERTED)")
    # (Id-3)
    for col, key, k in (("eddy", "E", 1), ("end", "S", 2)):
        syn = p79._synthetic()
        s = syn["cw"]["factor"]["val"]
        for i, lab in enumerate(PROBE_LABELS):
            v = s * float(pic[key][-1, i])
            syn["cw"]["probes"][lab] = {"tok": repr(v), "val": v}
        r = p76._rows(syn, pic)["cw"]["rows"]
        w = max(abs(r[lab]["ratio"][k] - 1.0) for lab in PROBE_LABELS)
        gates.append((("Id-3", col, w), w <= B1_BAND, w))
        _say(f"{PFX} (Id-3) synthetic s {s:.9e}, theta_acc planted s x {col}/P at the six points: "
             f"max |FEM/{col} - 1| {w:.3e} (band {B1_BAND:g}, ASSERTED; P76 read 0 / 2.220e-16)")
    _gates_assert(gates)


def test_controls_cknob_cmirror(monkeypatch):
    gates = []
    path, parsed = _parse60()
    real60 = _real_path(Path(path))
    # (c-knob): P79's unfiltered route on entry 60 (the only place the knob names entry 60)
    monkeypatch.setenv(p79.KNOB, str(real60))
    msg = None
    try:
        c46._on_root(p79._knob)
    except ValueError as exc:
        msg = str(exc)
    monkeypatch.delenv(p79.KNOB)
    raised = msg is not None and "sense cw lacks probes" in msg
    unf = c46._on_root(p79._parse_p61, real60)
    cw_n = len(unf.get("cw", {}).get("probes", {"x": 1}))
    ccw_eq = unf.get("ccw") == parsed["ccw"]
    ok = raised and cw_n == 0 and ccw_eq
    gates.append((("c-knob", msg, cw_n, ccw_eq), ok, cw_n))
    _say(f"{PFX} (c-knob) {p79.KNOB}={real60.name}: p79._knob() raised ValueError naming 'sense cw "
         f"lacks probes' {raised} ({msg!r}); unfiltered cw probes {cw_n} / 6 (filtered 6 / 6); "
         f"unfiltered ccw block == filtered {ccw_eq} (ASSERTED; ceiling: boolean, separation the "
         f"whole cw probe set)")
    # (c-mirror): the pictures at 3 T's c* equal those at entry 59's c*
    pic = p76._pictures()
    labs = [lab for lab, _pt, _w in pic["lat"]]
    j3, j15 = labs.index(WANT_STAR), labs.index(MIRROR_STAR)
    cols = {"U3": (p79._u3_hist_at(pic["lat"][j3][1]), p79._u3_hist_at(pic["lat"][j15][1])),
            "eddy": (np.asarray(pic["L"][:, j3]), np.asarray(pic["L"][:, j15])),
            "end": (np.asarray(pic["SL"][:, j3]), np.asarray(pic["SL"][:, j15]))}
    for k, (a, b) in cols.items():
        rel = abs(a[-1] / b[-1] - 1.0)
        dmax = float(np.max(np.abs(a - b)))
        gates.append((("c-mirror", k, rel, dmax), rel <= SYM_BAND and dmax <= SYM_BAND, rel))
        _say(f"{PFX} (c-mirror) {k}: theta({WANT_STAR})/theta({MIRROR_STAR}) - 1 at t* {rel:.3e}, "
             f"max |diff| over T_GRID {dmax:.3e} K/W (P83's (Sym) band {SYM_BAND:g}, ASSERTED)")
    _gates_assert(gates)


def _disc(got, want):
    return abs(got / want - 1.0)


def test_three_pictures_record_and_json():
    """The record (no band): asserts only that it was formed and written."""
    path, parsed = _parse60()
    pic = p76._pictures()
    t = p76._rows(parsed, pic)
    _say(f"{PFX} per deposited gel watt P_dep {P_DEP:.6f} W (all three pictures); entry 60 per "
         f"accepted watt: P_gel ccw {parsed['ccw']['P_gel']['tok']} W, cw {parsed['cw']['P_gel']['tok']}"
         f" W; theta per deposited watt = theta_acc / s, s = P_gel/P_acc of the same sense")
    _say(f"{PFX} TABLE point | sense | s | theta_acc K | theta_acc/s K/W | U3/P | eddy/P | end/P | "
         f"FEM/U3 | FEM/eddy | FEM/end")
    misses = []
    for s in ("ccw", "cw"):
        exp = EXPECT_CCW if s == "ccw" else EXPECT_CW
        for lab, r in t[s]["rows"].items():
            ru, re_, rn = r["ratio"]
            _say(f"{PFX} TABLE {lab} | {s} | {r['s']:.9e} | {r['theta_acc']:+.7g} | "
                 f"{r['fem_per_W']:.7g} | {r['U3_per_W']:.7g} | {r['eddy_per_W']:.7g} | "
                 f"{r['end_per_W']:.7g} | {ru:.7g} | {re_:.7g} | {rn:.7g}")
            if lab in exp:
                got = ((r["fem_per_W"],) if s == "ccw" else ()) + tuple(r["ratio"])
                names = (("fem",) if s == "ccw" else ()) + ("U3", "eddy", "end")
                misses += [(s, lab, nm, g, w, _disc(g, w)) for g, w, nm in zip(got, exp[lab], names)]
        h = t[s]["hot"]
        pu, pe, pn = (h["pictures"][k] for k in ("U3", "eddy", "end"))
        _say(f"{PFX} TABLE hotspot c* {h['label']} | {s} | {h['s']:.9e} | {h['theta_acc']:+.7g} | "
             f"{h['fem_per_W']:.7g} | {pu['per_W']:.7g} (u by _u3_hist_at) | {pe['per_W']:.7g} "
             f"(L[-1, {h['idx']}]) | {pn['per_W']:.7g} (SL[-1, {h['idx']}]) | {pu['ratio']:.7g} | "
             f"{pe['ratio']:.7g} | {pn['ratio']:.7g}")
        got = ((h["fem_per_W"],) if s == "ccw" else ()) + (pu["ratio"], pe["ratio"], pn["ratio"])
        names = (("fem",) if s == "ccw" else ()) + ("U3", "eddy", "end")
        misses += [(s, "hotspot", nm, g, w, _disc(g, w)) for g, w, nm in zip(got, exp["hotspot"], names)]
    over = [(s, lab, nm, f"{g:.7g}", w, f"{rel:.2e}") for s, lab, nm, g, w, rel in misses if rel > 1e-4]
    _say(f"{PFX} expected tables (brief, ccw 7 rows x 4, cw 4 rows x 3): {len(misses)} cells, max rel "
         f"miss {max(m[5] for m in misses):.3e}; cells above 1e-4 (disclosed, the brief's 4 s.f. "
         f"rounding reaches 5e-4): {over or 'none'}")
    # The hotspot dose
    for s in ("ccw", "cw"):
        h = t[s]["hot"]
        ex, exr = EXPECT_EXCESS[s]
        _say(f"{PFX} DOSE {s} c* {h['label']}: entry 60 cem43 per P_acc watt {h['cem43_tok']} s; "
             f"D0 = cem43(T_GRID, T_INF) {h['D0']:.15e} s (cited {p76.D0_CITED:.15e}, rel "
             f"{abs(h['D0'] / p76.D0_CITED - 1.0):.2e}); excess {h['excess']:.6g} s (expected "
             f"{ex:g}, rel {_disc(h['excess'], ex):.1e}), excess/D0 {h['excess'] / h['D0']:.5g} "
             f"(expected {exr:g}); at s {h['s']:.9e}: " + "; ".join(
                 f"{k} cem43 {p['cem43']:.9e} excess {p['excess']:.6g} s ({p['excess'] / h['D0']:.5g})"
                 for k, p in h["pictures"].items())
             + f" (DOSE_BAR {DOSE_BAR:g} unused; printed only)")
    hc = t["ccw"]["hot"]
    _say(f"{PFX} DOSE ccw pictures vs the brief's (P83's 3 T DOSE rows `…P83-W1.log:90–92`): " + "; ".join(
        f"{k} excess {hc['pictures'][k]['excess']:.5g} vs {w:g} (rel {_disc(hc['pictures'][k]['excess'], w):.1e})"
        for k, w in EXPECT_PIC_EXCESS_CCW.items()))
    _say(f"{PFX} DOSE note: a cem43 per accepted watt is approximately the unheated baseline D0 and "
         f"is not a scalable per-watt dose (cem43 is not linear in power); the excess is printed as "
         f"a record, never as a dose lever")
    # The six-minute window
    for s in ("ccw", "cw"):
        h = t[s]["hot"]
        still = h["ra_tmax"] == float(T_GRID[-1])
        _say(f"{PFX} WINDOW {s} c*: entry 60 running mean (W {W_B2:g} s) max {h['ra_max']:.9e} K at "
             f"t {h['ra_tmax']:.1f} s (t* {float(T_GRID[-1]):.1f} s; still rising {still}); at s: "
             + "; ".join(f"{k} max {p['ra_max']:.7g} K at t {p['ra_tmax']:.1f} s"
                         + (f" (expected {EXPECT_PIC_WINDOW_CCW[k]:g}, rel "
                            f"{_disc(p['ra_max'], EXPECT_PIC_WINDOW_CCW[k]):.1e})" if s == "ccw" else "")
                         for k, p in h["pictures"].items())
             + f" (AVG_BAR {AVG_BAR:g} unused; printed only)")
    # Notes
    rc = t["ccw"]["rows"]
    _say(f"{PFX} NOTE c* is C6a's ccw SAR10g lattice argmax at 3 T, not the thermal maximum: per "
         f"accepted watt, ccw, P1 +y / P2 -y read {rc['P1 +y wall']['theta_acc'] / hc['theta_acc']:.4g}x"
         f" / {rc['P2 -y wall']['theta_acc'] / hc['theta_acc']:.4g}x c*'s theta (expected 1.953x / "
         f"1.862x); the cw arm is heated at ccw's c*")
    _say(f"{PFX} NOTE probe labels are P37's box frame; the printed coordinates are the fixture's "
         f"(P1 +y wall prints at (0.21, 0.0, 0.0)); rows are matched by label, never by coordinate")
    _say(f"{PFX} NOTE the thermal constants are the 1.5 T gel's (C31, entry 60 :{C31_LINE_NO}) on "
         f"both sides; the only 3 T input is s; no 3 T gel claim")
    doc = {"entry60_log": path, "P_dep_W": P_DEP, "D0_s": hc["D0"],
           "normaliser": "pictures per deposited gel watt; entry 60 per accepted watt, s = P_gel/P_acc",
           "senses": t}
    if MPI.COMM_WORLD.rank == 0:
        JSON_PATH.parent.mkdir(parents=True, exist_ok=True)
        JSON_PATH.write_text(json.dumps(doc, indent=1, default=list))
    MPI.COMM_WORLD.barrier()
    back = json.loads(JSON_PATH.read_text())
    _say(f"{PFX} JSON {JSON_PATH} keys {sorted(back)}; senses {sorted(back['senses'])}")
    assert back["entry60_log"] == path and sorted(back["senses"]) == ["ccw", "cw"]
    assert all(np.isfinite(r["ratio"]).all() for s in t for r in t[s]["rows"].values())
