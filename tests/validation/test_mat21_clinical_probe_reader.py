"""`MAT-21` step 43 — P76's reader with entry 59's columns blank (§9 item 898, §10 P79;
brief `docs/planning/chunks/MAT-21.md` "§9 item 898").

Complex build (imports P61's module), no solve.  Three parts:

- ``_parse_p61(path)`` — the lines P61's `_run_heat_box`
  (`test_geo38_f_clinical_box_heating.py:240–255`) prints: per sense key (``cw`` /
  ``ccw``) the ``sources, sense`` line (P_gel, P_acc, factor), the ``hotspot c*`` line
  (label, rounded c*, θ_acc(t*), final-field read, cem43, running-mean max and its
  time) and every ``probe <label>, sense`` line.  Every value keeps its printed token.
  A sense block printed twice must match byte for byte, else ``ValueError``.
- ``_to_box(X, off)`` — the inverse of I32's `_map_box_point`
  (`test_geo36_rod_f_clinical_heating.py:145–152`): (X, Y, Z) → (Z, X, Y − off),
  refusing (``None`` + reason) a point outside the box.
- ``_p76_rows(parsed)`` — P76's six rows (P75's uniform-source figures, P52's FEM
  sourced floor, entry 59's columns per sense, ``None`` when ``parsed`` is ``None``)
  plus a hotspot row, written to ``logs/mat21/p76-table.json`` by rank 0.

Normalisers: P75 is per deposited gel watt (P_dep = 49.14 W); entry 59 is per accepted
watt, θ_acc = θ·P_gel/P_acc.  The factor s = P_gel/P_acc is stated on every row.

**Placement disclosed:** P61 prints a cem43 only at the hotspot (no probe history is
printed), so the entry-59 cem43 column exists only on the hotspot row; each probe row
carries the uniform comparand's cem43 per accepted watt (s·U3_hist/P_dep) with the
entry-59 cem43 cell ``None`` ("not printed by P61").

Knob ``FEM_EM_MAT21_P76_LOG`` (read inside the test): unset/empty ⇒ entry 59's columns
blank (this step); a path must exist and parse with all six ``PROBE_LABELS`` on every
sense present, else ``ValueError``.

Not claimed: any entry-59 figure or ratio, a non-uniformity reading, entry 60, a dose
band, a compliance figure.
"""

from __future__ import annotations

import gzip
import json
import os
import re
from pathlib import Path

import numpy as np
from mpi4py import MPI

from fem_em_solver.thermal.dose import cem43, running_average
from fem_em_solver.utils.analytical import robin_box_transient_temperature

from tests.validation import test_geo38_f_clinical_box_heating as p61
from tests.validation import test_geo38_f_clinical_box_loaded as c6
from tests.validation import test_mat21_wall_uniform_source_probe_points_closed_form as p75
from tests.validation.test_mat21_box_closed_form import (
    A_BOX,
    B_BOX,
    C_BOX,
    HALF,
    SERIES_TERMS,
    H,
    K,
    Q,
)
from tests.validation.test_mat21_box_closed_form import RHO_C as RHO_C_0

_hists, _u3, P_DEP = p75._hists, p75._u3, p75.P_DEP
PROBE_PTS, PROBE_LABELS = p75.PROBE_PTS, p75.PROBE_LABELS
T_GRID, T_INF, T_STAR = p75.T_GRID, p75.T_INF, p75.T_STAR
B1_BAND, DOSE_BAR, W_B2 = p75.B1_BAND, p75.DOSE_BAR, p75.W_B2

KNOB = "FEM_EM_MAT21_P76_LOG"
REPO = Path(__file__).resolve().parents[2]
JSON_PATH = REPO / "logs" / "mat21" / "p76-table.json"
LOGS = REPO / "docs" / "testing" / "logs"
STUB_LOG = LOGS / "20261105T103249Z_GEO-38-P61b-W0.log"
EPS = np.finfo(float).eps
SENSE_KEY = {v: k for k, v in p61.SENSE_TUPLE.items()}

# P52's FEM sourced value at the six points, K at P_dep (unwalled six-face box, m = 6,
# P1, CN N = 40; `20261101T150646Z_MAT-21-P52-W1.log:76–81`, the t_w = 0 rows).
P52_FEM = (1.337056, 1.065488, 1.065488, 1.065258, 1.065258, 1.074071)
# The same rows' "P51 sourced" column (K), (Rep-52)'s reference.
P51_SRC = (1.335697, 1.061827, 1.061827, 1.061678, 1.061678, 1.072345)
REP52_TOL = 5e-7
# (Rep): P75's printed digits (`20261108T134258Z_MAT-21-P75-W1.log:71–76`, `:86`, `:91`).
REP_U3 = {0: 0.02718146, 1: 0.02160819, 3: 0.02160517, 5: 0.02182225}  # K/W, within 5e-9
REP_WS = {0: (0.02701082, 5e-9), 5: (0.0211361, 5e-8)}                 # K/W
REP_CEM = {0: 4.086591e-11, 5: 4.073623e-11}                           # s/W, within 5e-18
# (Id): the stub's cw factor (`…P61b-W0.log:129`).
S_STUB = 2.858205069e-05
PLANT = 1e-6

_PFX = r"\[GEO-38 P61\] "
_TUP = r"(\(1,[+-]1j\))"
RE_SRC = re.compile(_PFX + r"sources, sense " + _TUP + r" \((cw|ccw), heated\): P_gel (\S+) W, "
                    r"P_acc (\S+) W \(lumped line cut, ruling 2026-10-12\); per-accepted-watt "
                    r"factor P_gel/P_acc \((cw|ccw)\) (\S+)$")
RE_HOT = re.compile(_PFX + r"hotspot c\* (.+?) \((.+?)\), sense " + _TUP + r": theta\(t\*\) per "
                    r"P_acc watt (\S+) K \(point; final-field read (\S+) K\); cem43 per P_acc "
                    r"watt (\S+) s .*running mean \(W \d+ s\) max (\S+) K at t (\S+) s \(printed\)$")
RE_PRB = re.compile(_PFX + r"probe (.+?), sense " + _TUP + r": theta\(t\*\) per P_acc watt (\S+) K "
                    r"\(printed\)$")
FMT = {"P_gel": ".12e", "P_acc": ".12e", "factor": ".9e", "th": ".9e", "final": ".9e",
       "cem43": ".9e", "ra_max": ".9e", "ra_tmax": ".1f", "probe": "+.9e"}


def _say(msg):
    if MPI.COMM_WORLD.rank == 0:
        print(msg, flush=True)


# ------------------------------------------------------------------ the map
def _to_box(X, off):
    """(X, Y, Z) in the fixture frame → (a, b, c) = (Z, X, Y − off); ``(None, reason)``
    outside the box (never evaluated there)."""
    x, y, z = (float(v) for v in X)
    p = (z, x, y - off)
    if abs(p[0]) <= A_BOX and abs(p[1]) <= B_BOX and abs(p[2]) <= C_BOX:
        return p, None
    return None, f"{tuple(X)} maps to {p}, outside |a|<={A_BOX}, |b|<={B_BOX}, |c|<={C_BOX}"


def _u3_at(pt, t):
    a, b, c = pt
    return float(np.asarray(robin_box_transient_temperature(
        np.array([a]), np.array([b]), np.array([c]), t, half_widths=HALF, k=K, q=Q, h=H,
        rho_c=RHO_C_0, t_inf=0.0, n_terms=SERIES_TERMS), dtype=float)[0])


def _u3_hist_at(pt):
    """Per-P_dep-watt history on T_GRID, row 0 = 0 (never evaluated)."""
    out = np.zeros(T_GRID.size)
    for n in range(1, T_GRID.size):
        out[n] = _u3_at(pt, float(T_GRID[n]))
    return out / P_DEP


# ------------------------------------------------------------------ the parser
def _tok(s):
    return {"tok": s, "val": float(s)}


def _parse_p61_text(text):
    """Per sense key: tokens/values of the sources, hotspot and probe lines, plus the
    1-based line numbers of every occurrence of its block."""
    lines = text.splitlines()
    blocks = []  # [sense, [line numbers], [raw lines], record]
    cur = None
    for n, raw in enumerate(lines, start=1):
        m = RE_SRC.search(raw)
        if m:
            tup, s, pg, pa, s2, fac = m.groups()
            if SENSE_KEY[tup] != s or s2 != s:
                raise ValueError(f"line {n}: sense keys disagree {tup!r} {s!r} {s2!r}")
            cur = [s, [n], [raw], {"tuple": tup, "P_gel": _tok(pg), "P_acc": _tok(pa),
                                   "factor": _tok(fac), "hot": None, "probes": {}}]
            blocks.append(cur)
            continue
        m = RE_HOT.search(raw) or RE_PRB.search(raw)
        if m and cur is not None and m.re is RE_HOT:
            lab, pt, tup, th, fin, dose, ram, rat = m.groups()
            if tup != cur[3]["tuple"]:
                raise ValueError(f"line {n}: hotspot sense {tup} inside {cur[3]['tuple']} block")
            cur[1].append(n)
            cur[2].append(raw)
            cur[3]["hot"] = {"label": lab, "point_tok": pt,
                             "point": tuple(float(v) for v in pt.split(",")),
                             "th": _tok(th), "final": _tok(fin), "cem43": _tok(dose),
                             "ra_max": _tok(ram), "ra_tmax": _tok(rat)}
            continue
        if m and cur is not None and m.re is RE_PRB:
            lab, tup, th = m.groups()
            if tup != cur[3]["tuple"]:
                raise ValueError(f"line {n}: probe sense {tup} inside {cur[3]['tuple']} block")
            cur[1].append(n)
            cur[2].append(raw)
            cur[3]["probes"][lab] = _tok(th)
            continue
        cur = None  # a block is the consecutive run from its sources line
    out = {}
    for s, nums, raws, rec in blocks:
        text_b = "\n".join(raws)
        if s in out:
            if text_b != out[s]["_text"]:
                raise ValueError(f"sense {s}: block at lines {nums} differs from the block at "
                                 f"lines {out[s]['_lines'][0]}")
            out[s]["_lines"].append(nums)
            continue
        if rec["hot"] is None:
            raise ValueError(f"sense {s}: block at lines {nums} has no hotspot line")
        rec["_text"], rec["_lines"] = text_b, [nums]
        out[s] = rec
    return out


def _parse_p61(path):
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as f:
        return _parse_p61_text(f.read())


def _reformat_misses(parsed):
    """Every parsed value re-formatted with its print format vs its token (strings)."""
    bad, n = [], 0
    for s, rec in parsed.items():
        items = [(k, rec[k]) for k in ("P_gel", "P_acc", "factor")]
        items += [(k, rec["hot"][k]) for k in ("th", "final", "cem43", "ra_max", "ra_tmax")]
        items += [("probe", v) for v in rec["probes"].values()]
        for k, v in items:
            n += 1
            if format(v["val"], FMT[k]) != v["tok"]:
                bad.append((s, k, v["tok"], format(v["val"], FMT[k])))
    return bad, n


def _knob():
    raw = os.environ.get(KNOB, "").strip()
    if raw == "":
        return None, None
    path = Path(raw)
    if not path.is_file():
        raise ValueError(f"{KNOB}={raw!r}: no such file")
    parsed = _parse_p61(path)
    if not parsed:
        raise ValueError(f"{KNOB}={raw!r}: no P61 sense block found")
    for s, rec in parsed.items():
        missing = [lab for lab in PROBE_LABELS if lab not in rec["probes"]]
        if missing:
            raise ValueError(f"{KNOB}={raw!r}: sense {s} lacks probes {missing}")
    return str(path), parsed


# ------------------------------------------------------------------ the table
def _hot_row(parsed):
    """Hotspot row per sense: c*'s lattice point by label, U3 there, entry 59's hotspot
    columns.  ``None`` + reason when the label is not on C6's lattice."""
    if parsed is None:
        return None
    lat = dict(c6._lattice())
    off = float(c6._box()[0][1])
    h = _hists()
    out = {}
    for s, rec in parsed.items():
        hot = rec["hot"]
        if hot["label"] not in lat:
            out[s] = {"row": None, "reason": f"c* label {hot['label']!r} is not on C6's lattice"}
            continue
        X = lat[hot["label"]]
        dpt = max(abs(a - b) for a, b in zip(X, hot["point"]))
        if dpt > 5e-7:
            raise ValueError(f"c* {hot['label']}: printed {hot['point']} vs lattice {X}, {dpt:.3e} m")
        pt, why = _to_box(X, off)
        if pt is None:
            out[s] = {"row": None, "reason": why}
            continue
        uh = _u3_hist_at(pt)
        fac = rec["factor"]["val"]
        comp = fac * uh[-1]
        out[s] = {"row": {
            "label": hot["label"], "fixture_pt": X, "box_pt": pt, "print_vs_lattice_m": dpt,
            "U3_per_W": uh[-1],
            "Wsrc_per_W_z0": float(h["Wsrc_raw"][-1, 2]) / P_DEP, "Wdrift_K": 0.0,
            "s": fac, "theta_acc": hot["th"]["val"], "comparand": comp,
            "ratio": hot["th"]["val"] / comp,
            "cem43_e59": hot["cem43"]["val"],
            "cem43_comp": float(cem43(T_GRID, T_INF + fac * uh)),
        }, "reason": None}
    return out


def _p76_rows(parsed=None):
    h = _hists()
    u3t = _u3(T_STAR) / P_DEP
    rows = {}
    for i, lab in enumerate(PROBE_LABELS):
        uh = h["U3"][:, i]
        row = {
            "box_pt": tuple(float(v) for v in PROBE_PTS[i]),
            "tag": "1D: z = 0" if i in (1, 2, 3, 4) else "",
            "U3_per_W": float(u3t[i]),
            "Wsrc_per_W": float(h["Wsrc"][-1, i]),
            "Wsd_per_W": float(h["Wsd"][-1, i]),
            "cem43_per_W_U3": float(cem43(T_GRID, T_INF + uh)),
            "ra_max_per_W_U3": float(running_average(T_GRID, uh, window=W_B2)[1]),
            "P52_FEM_per_W": P52_FEM[i] / P_DEP,
            "e59": None,
        }
        if parsed is not None:
            row["e59"] = {}
            for s, rec in parsed.items():
                fac = rec["factor"]["val"]
                th = rec["probes"][lab]["val"]
                comp = fac * float(u3t[i])
                row["e59"][s] = {"s": fac, "theta_acc": th, "comparand": comp,
                                 "ratio": th / comp, "cem43_e59": None,
                                 "cem43_e59_note": "not printed by P61 at the probes",
                                 "cem43_comp": float(cem43(T_GRID, T_INF + fac * uh))}
        rows[lab] = row
    return {"rows": rows, "hot": _hot_row(parsed)}


def _synthetic(s=S_STUB, plant=None):
    """A parsed record (cw) whose probe θ_acc are s·U3(t*)/P_dep, its hotspot at a
    lattice point with θ_acc and cem43 of s·U3_hist/P_dep there.  Exact floats."""
    u3t = _u3(T_STAR) / P_DEP
    probes = {}
    for i, lab in enumerate(PROBE_LABELS):
        v = s * float(u3t[i]) * (1.0 + (plant[1] if plant and plant[0] == i else 0.0))
        probes[lab] = {"tok": repr(v), "val": v}
    lab, X = c6._lattice()[9]
    pt, _why = _to_box(X, float(c6._box()[0][1]))
    uh = _u3_hist_at(pt)
    d = float(cem43(T_GRID, T_INF + s * uh))
    f = lambda v: {"tok": repr(v), "val": v}  # noqa: E731
    hot = {"label": lab, "point_tok": "", "point": tuple(round(v, 6) for v in X),
           "th": f(s * uh[-1]), "final": f(s * uh[-1]), "cem43": f(d), "ra_max": f(0.0),
           "ra_tmax": f(0.0)}
    return {"cw": {"tuple": p61.SENSE_TUPLE["cw"], "P_gel": f(s), "P_acc": f(1.0), "factor": f(s),
                   "hot": hot, "probes": probes}}


# ------------------------------------------------------------------ tests
def test_k_map_and_preconditions():
    gates = []
    off = float(c6._box()[0][1])
    bound = 4.0 * EPS * (abs(off) + C_BOX)
    fx = p61._probes()
    for i, lab in enumerate(PROBE_LABELS):
        pt, why = _to_box(fx[lab], off)
        d = np.inf if pt is None else max(abs(a - b) for a, b in zip(pt, PROBE_PTS[i]))
        _say(f"MAT21-43 (K-map) {lab}: fixture {fx[lab]} -> box {pt} vs {PROBE_PTS[i]}: |d| {d:.3e} m "
             f"(bound {bound:.3e}, ASSERTED){'' if why is None else ' ' + why}")
        gates.append((("K-map", lab), d <= bound, d))
    lat = c6._lattice()
    worst_c = 0.0
    for lab, X in lat:
        pt, why = _to_box(X, off)
        if pt is None:
            gates.append((("lattice", lab, why), False, None))
            continue
        worst_c = max(worst_c, abs(pt[2]))
    gates.append((("lattice c=0",), worst_c <= bound, worst_c))
    _say(f"MAT21-43 (K-map) off {off!r} m; {len(lat)} lattice points all inside the box, max |c| "
         f"{worst_c:.3e} m (bound {bound:.3e}, ASSERTED)")
    ref = np.arange(41) * p61.T_STAR_BOX / 40.0
    tg = float(np.max(np.abs(T_GRID - ref)) / p61.T_STAR_BOX)
    gates.append((("T_GRID",), T_GRID.size == 41 and tg <= 1e-12, tg))
    rw = p61.RUNNING_WINDOW_S == 360.0 == float(W_B2)
    gates.append((("window",), rw, (p61.RUNNING_WINDOW_S, W_B2)))
    _say(f"MAT21-43 (pre) T_GRID {T_GRID.size} rows vs arange(41) T_STAR_BOX/40: max rel {tg:.3e} "
         f"(bound 1e-12); RUNNING_WINDOW_S {p61.RUNNING_WINDOW_S} == 360.0 == W_B2 {W_B2}: {rw}")
    # (c-perm): the un-permuted inverse misses P1 by sqrt(2) B.
    X1 = fx[PROBE_LABELS[1]]
    up = (X1[0], X1[1] - off, X1[2])
    dp = float(np.linalg.norm(np.subtract(up, PROBE_PTS[1])))
    _say(f"MAT21-43 (c-perm) un-permuted inverse at P1: {up} vs {PROBE_PTS[1]}, miss {dp:.6f} m "
         f"(sqrt(2) B {np.sqrt(2) * B_BOX:.6f}; ASSERTED > 1e-6)")
    gates.append((("c-perm",), dp > 1e-6, dp))
    bad = [g for g in gates if not g[1]]
    assert not bad, bad


def test_k_parse_on_the_stub_log():
    gates = []
    parsed = _parse_p61(STUB_LOG)
    bad, n = _reformat_misses(parsed)
    gates.append((("reformat",), not bad and n > 0, bad))
    lines = {s: r["_lines"] for s, r in parsed.items()}
    want = {"cw": [[129, 130, 131, 132], [188, 189, 190, 191]], "ccw": [[178, 179, 180, 181]]}
    gates.append((("lines",), lines == want, lines))
    pg, dose = parsed["cw"]["P_gel"]["val"], parsed["cw"]["hot"]["cem43"]["val"]
    gates.append((("P_gel cw",), parsed["cw"]["P_gel"]["tok"] == "9.765625000000e-07", pg))
    gates.append((("cem43 cw",), parsed["cw"]["hot"]["cem43"]["tok"] == "4.009604022e-11", dose))
    hr = _hot_row(parsed)
    reason = hr["cw"]["reason"]
    gates.append((("hot None",), hr["cw"]["row"] is None and "k8-1" in (reason or ""), reason))
    _say(f"MAT21-43 (K-parse) {STUB_LOG.name}: senses {sorted(parsed)}; block lines {lines} "
         f"(want {want}); {n} values re-formatted, misses {bad or 'none'}; P_gel(cw) "
         f"{parsed['cw']['P_gel']['tok']} W, cem43(cw) {parsed['cw']['hot']['cem43']['tok']} s; "
         f"probes cw {list(parsed['cw']['probes'])}; hotspot row: {reason} (ASSERTED)")
    # (c-dup): one digit changed in the second cw block.
    text = STUB_LOG.read_text().splitlines()
    old = text[187]
    text[187] = old.replace("9.765625000000e-07", "9.765625000001e-07", 1)
    raised = None
    try:
        _parse_p61_text("\n".join(text))
    except ValueError as e:
        raised = str(e)
    _say(f"MAT21-43 (c-dup) line 188 edited ({old != text[187]}): ValueError {raised!r} (ASSERTED)")
    gates.append((("c-dup",), old != text[187] and raised is not None, raised))
    bad_g = [g for g in gates if not g[1]]
    assert not bad_g, bad_g


def test_rep_and_rep52():
    gates = []
    h = _hists()
    u = _u3(T_STAR)
    for i, v in REP_U3.items():
        d = abs(float(u[i]) / P_DEP - v)
        gates.append((("U3/P", i), d <= 5e-9, d))
        _say(f"MAT21-43 (Rep) U3/P {PROBE_LABELS[i]} t*: {float(u[i]) / P_DEP:.12e} vs {v} |d| {d:.3e} "
             f"(bound 5e-9)")
    for i, (v, tol) in REP_WS.items():
        d = abs(float(h["Wsrc"][-1, i]) - v)
        gates.append((("Wsrc/P", i), d <= tol, d))
        _say(f"MAT21-43 (Rep) W_src/P {PROBE_LABELS[i]} t*: {float(h['Wsrc'][-1, i]):.12e} vs {v} "
             f"|d| {d:.3e} (bound {tol:g})")
    for i, v in REP_CEM.items():
        c = float(cem43(T_GRID, T_INF + h["U3"][:, i]))
        d = abs(c - v)
        gates.append((("cem43/W", i), d <= 5e-18, d))
        _say(f"MAT21-43 (Rep) cem43/W U3 {PROBE_LABELS[i]}: {c:.12e} s vs {v} |d| {d:.3e} "
             f"(bound 5e-18)")
    for i, lab in enumerate(PROBE_LABELS):
        d = abs(float(u[i]) - P51_SRC[i])
        gates.append((("Rep-52", lab), d <= REP52_TOL, d))
        _say(f"MAT21-43 (Rep-52) {lab}: _u3(t*) {float(u[i]):.9f} K vs P52's P51 sourced "
             f"{P51_SRC[i]:.6f} |d| {d:.3e} (bound {REP52_TOL:g})")
    bad = [g for g in gates if not g[1]]
    assert not bad, bad


def test_id_and_c_plant_on_a_synthetic_entry59():
    gates = []
    t = _p76_rows(_synthetic())
    miss = {lab: abs(r["e59"]["cw"]["ratio"] - 1.0) for lab, r in t["rows"].items()}
    w = max(miss.values())
    gates.append((("Id ratio",), w <= B1_BAND, w))
    hr = t["hot"]["cw"]["row"]
    dm = abs(hr["cem43_e59"] / hr["cem43_comp"] - 1.0)
    hm = abs(hr["ratio"] - 1.0)
    gates.append((("Id cem43",), dm <= DOSE_BAR, dm))
    gates.append((("Id hot ratio",), hm <= B1_BAND, hm))
    _say(f"MAT21-43 (Id) s {S_STUB:.9e}: max |ratio-1| over six points {w:.3e} (band {B1_BAND:g}, "
         f"ASSERTED; predicted <= 1e-15); hotspot {hr['label']} box {hr['box_pt']}: |ratio-1| "
         f"{hm:.3e}, cem43 pair {hr['cem43_e59']:.9e} / {hr['cem43_comp']:.9e} s rel {dm:.3e} "
         f"(bar {DOSE_BAR:g}, ASSERTED)")
    tp = _p76_rows(_synthetic(plant=(3, PLANT)))
    mp = {lab: r["e59"]["cw"]["ratio"] - 1.0 for lab, r in tp["rows"].items()}
    m3 = mp[PROBE_LABELS[3]]
    others = max(abs(v) for k, v in mp.items() if k != PROBE_LABELS[3])
    gates.append((("c-plant ceiling",), abs(m3 - PLANT) <= 1e-12, m3))
    gates.append((("c-plant > band",), m3 > B1_BAND, m3))
    gates.append((("c-plant others",), others <= B1_BAND, others))
    _say(f"MAT21-43 (c-plant) P3 x (1 + {PLANT:g}): miss {m3:.12e} (|miss - plant| "
         f"{abs(m3 - PLANT):.3e} <= 1e-12, > {B1_BAND:g}, ASSERTED); others max {others:.3e} "
         f"(<= {B1_BAND:g}, ASSERTED)")
    bad = [g for g in gates if not g[1]]
    assert not bad, bad


def _fmt(v, f=".7g"):
    return "—" if v is None else format(v, f)


def test_table_and_json():
    """The record (no band): entry 59's columns blank unless the knob names a log."""
    path, parsed = _knob()
    t = _p76_rows(parsed)
    _say(f"MAT21-43 knob {KNOB}: {path!r}; P_dep {P_DEP:.6f} W; entry 59 normaliser: per accepted "
         f"watt, factor s = P_gel/P_acc stated per row")
    _say("MAT21-43 TABLE point | box pt (m) | U3/P K/W | W_src/P | W_src/P+W_drift | cem43/W U3 s "
         "| max 360 s mean/W U3 K | P52 FEM/P | FEM/U3-1 | e59 sense | s | theta_acc | s U3/P | "
         "ratio | cem43 e59 | cem43 s U3/P")
    for lab, r in t["rows"].items():
        tag = f" ({r['tag']})" if r["tag"] else ""
        head = (f"MAT21-43 TABLE {lab}{tag} | {r['box_pt']} | {_fmt(r['U3_per_W'])} | "
                f"{_fmt(r['Wsrc_per_W'])} | {_fmt(r['Wsd_per_W'])} | {_fmt(r['cem43_per_W_U3'])} | "
                f"{_fmt(r['ra_max_per_W_U3'])} | {_fmt(r['P52_FEM_per_W'])} | "
                f"{r['P52_FEM_per_W'] / r['U3_per_W'] - 1.0:+.3e}")
        if r["e59"] is None:
            _say(head + " | — | — | — | — | — | — | —")
            continue
        for s, e in r["e59"].items():
            _say(head + f" | {s} | {e['s']:.9e} | {e['theta_acc']:+.7g} | {e['comparand']:.7g} | "
                        f"{e['ratio']:.7g} | — ({e['cem43_e59_note']}) | {e['cem43_comp']:.7g}")
    if t["hot"] is None:
        _say("MAT21-43 TABLE hotspot c* | — (entry 59's log not given) | — | — | — | — | — | — | —")
    else:
        for s, hh in t["hot"].items():
            r = hh["row"]
            if r is None:
                _say(f"MAT21-43 TABLE hotspot c* ({s}) | None: {hh['reason']}")
                continue
            _say(f"MAT21-43 TABLE hotspot c* {r['label']} ({s}) | {r['box_pt']} | {r['U3_per_W']:.7g} | "
                 f"{r['Wsrc_per_W_z0']:.7g} (z = 0) | drift 0 | s {r['s']:.9e} | {r['theta_acc']:+.7g} "
                 f"| {r['comparand']:.7g} | {r['ratio']:.7g} | {r['cem43_e59']:.7g} | "
                 f"{r['cem43_comp']:.7g}")
    rows = t["rows"]
    _say("MAT21-43 U3/P vs P52 FEM/P, FEM/U3-1: " + ", ".join(
        f"{lab.split()[0]} {r['P52_FEM_per_W']:.7g}/{r['U3_per_W']:.7g} "
        f"{r['P52_FEM_per_W'] / r['U3_per_W'] - 1.0:+.3e}" for lab, r in rows.items()))
    doc = {"entry59_log": path, "P_dep_W": P_DEP, "T_STAR_s": T_STAR,
           "normaliser": "P75 per deposited gel watt; entry 59 per accepted watt, s = P_gel/P_acc",
           "rows": rows, "hot": t["hot"]}
    if MPI.COMM_WORLD.rank == 0:
        JSON_PATH.parent.mkdir(parents=True, exist_ok=True)
        JSON_PATH.write_text(json.dumps(doc, indent=1, default=list))
    MPI.COMM_WORLD.barrier()
    back = json.loads(JSON_PATH.read_text())
    _say(f"MAT21-43 JSON {JSON_PATH} keys {sorted(back)}; row keys {sorted(back['rows'][PROBE_LABELS[0]])}")
    assert back["entry59_log"] == path and sorted(back["rows"]) == sorted(PROBE_LABELS)
