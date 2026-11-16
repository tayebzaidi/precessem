"""`MAT-21` step 43b — P76: entry 59's six probe θ and hotspot θ / cem43 / six-minute window
per accepted watt beside the three closed-form pictures per deposited watt (§9 item 936,
§10 P76; brief `docs/planning/chunks/MAT-21.md` "§9 item 936").

Complex build (P79 imports P61 / C6), no solve.  Entry 59 is read **only** through
``p79._knob()`` (knob ``FEM_EM_MAT21_P76_LOG``, read inside a test, never at import).  The
pictures: P75's uniform U3 (``p79._p76_rows`` / ``p79._u3_hist_at``), P80's eddy source
(``p80._data()`` ``E`` / ``L``) and P81's end-weighted source (``p81._data()`` ``S`` /
``SL``), all per deposited gel watt (P_dep = 49.14 W).  Entry 59 is per accepted watt,
θ_acc = θ·P_gel/P_acc, so θ per deposited watt = θ_acc / s with s = P_gel/P_acc stated on
every row.

Asserted (collected, then asserted; print-precision custody, C33's (D-*) pattern, no new
band): (K-59) P79's reformat 0 misses / 28 values, block lines, ccw c* on C6's lattice;
(K-s) the factor against its two powers; (K-norm) the ccw powers against entry 59's own C6a
normalisers; (K-twin) θ_acc against factor × the I32 gel twin's per-watt θ; (R52) entry 59's
C6a lattice / c* / normalisers / S2 strings equal entry 52's; (Id-3) P79's synthetic record
planted on the eddy and end columns reads ratio 1 within ``B1_BAND``.

Printed, no band (a 🟡 record): the per-sense table, the hotspot dose against the unheated
baseline D0, the six-minute window.  Not claimed: the hotspot calorimetric slope, entry 59's
decision rule, a per-probe dose, cem43 per accepted watt as a dose lever, which picture the
heating follows, a fitted source, entry 60, a compliance or tissue claim.
"""

from __future__ import annotations

import json
import re

import numpy as np
from mpi4py import MPI

from fem_em_solver.thermal.dose import cem43, running_average

from tests.validation import test_geo38_c33_field_strength_reader as c33
from tests.validation import test_geo38_f_clinical_box_loaded as c6
from tests.validation import test_mat21_clinical_probe_reader as p79
from tests.validation import test_mat21_eddy_source_probe_points_closed_form as p80
from tests.validation import test_mat21_end_source_probe_points_closed_form as p81
from tests.validation.test_imp4_pec_cylinder_closed_form import B1_BAND
from tests.validation.test_mat21_box_closed_form import T_INF
from tests.validation.test_mat21_box_dose_closed_form import AVG_BAR, DOSE_BAR, T_GRID, W_B2

PROBE_LABELS, P_DEP = p79.PROBE_LABELS, p79.P_DEP
JSON_PATH = p79.REPO / "logs" / "mat21" / "p76-three-pictures.json"
WANT_LINES = {"ccw": [list(range(12325, 12333))], "cw": [list(range(12343, 12351))]}
WANT_STAR = "x-105mm dz-50mm"
N_VALUES = 28  # 2 senses x (3 + 5 + 6)
# The baseline D0 = cem43(T_GRID, T_INF) = t* 0.25^23, cited
# (`20261020T190033Z_MAT-21-step15-W1-probe-dose-closed-form.log.gz:83`), printed beside.
D0_CITED = 4.009592657894245e-11
# Brief's expected ccw FEM/picture ratios (printed beside; a miss > 1e-4 rel is disclosed).
EXPECT_CCW = {
    "P0 centre": (0.1108, 2.331, 5.967), "P1 +y wall": (3.892, 1.561, 219.5),
    "P2 -y wall": (3.599, 1.443, 203.0), "P3 +x wall": (0.6548, 14.38, 0.3305),
    "P4 -x wall": (0.7066, 15.52, 0.3567), "P5 top face": (0.2166, 2.318, 12.14),
    "hotspot": (1.122, 1.468, 8.566),
}
RE_TWIN = re.compile(r"^\[GEO-36 I32\] \[gel twin (cw|ccw)\] probe (.+?) \(([^)]*)\), sense "
                     r"(\(1,[+-]1j\)): theta\(t\*\) (\S+) K \(valid True; printed\)$")
RE_TWIN_GP = re.compile(r"^\[GEO-36 I32\] \[gel twin (cw|ccw)\] sense (\(1,[+-]1j\)): P_gel (\S+) W; "
                        r"G (\S+) 1/W; G\*P (\S+) W;")
R52_KEYS = ("c_star", "normalisers", "s2")


def _say(msg):
    if MPI.COMM_WORLD.rank == 0:
        print(msg, flush=True)


def _read_text(path):
    page = c33._load_page()
    return page._read_log(path)


def _parse_twin(path):
    """Per sense: {label: token} of the I32 gel-twin probe lines, and the G*P tokens."""
    probes, gp = {}, {}
    for n, raw in enumerate(_read_text(path).splitlines(), start=1):
        m = RE_TWIN.search(raw)
        if m:
            s, lab, _pt, tup, th = m.groups()
            if p79.SENSE_KEY[tup] != s:
                raise ValueError(f"line {n}: twin sense {s} vs {tup}")
            if lab in probes.setdefault(s, {}):
                raise ValueError(f"line {n}: twin {s} {lab} printed twice")
            probes[s][lab] = th
            continue
        m = RE_TWIN_GP.search(raw)
        if m:
            gp.setdefault(m.group(1), []).append((n, m.group(3), m.group(5)))
    return probes, gp


def _strip(d):
    return {k: v for k, v in d.items() if k != "line"}


def _r52_fields(parsed):
    lat = [(r["label"], r["b1p_ccw"], r["sar10g"], r["mass_kg"]) for r in parsed["lattice"]]
    return {"lattice": lat, **{k: _strip(parsed[k]) for k in R52_KEYS}}


def _pictures():
    d80, d81 = p80._data(), p81._data()
    return {"E": d80["E"], "L": d80["L"], "lat": d80["lat"], "S": d81["S"], "SL": d81["SL"]}


def _rows(parsed, pic):
    """Per sense: six probe rows + the hotspot row (FEM per deposited watt vs three pictures),
    the hotspot dose and the six-minute window.  Entry 59's figures are per accepted watt."""
    t3 = p79._p76_rows(None)["rows"]
    labs = [lab for lab, _pt, _w in pic["lat"]]
    out = {}
    for s, rec in parsed.items():
        fac = rec["factor"]["val"]
        rows = {}
        for i, lab in enumerate(PROBE_LABELS):
            th = rec["probes"][lab]["val"]
            fem = th / fac
            u, e, en = t3[lab]["U3_per_W"], float(pic["E"][-1, i]), float(pic["S"][-1, i])
            rows[lab] = {"s": fac, "theta_acc": th, "fem_per_W": fem, "U3_per_W": u,
                         "eddy_per_W": e, "end_per_W": en,
                         "ratio": (fem / u, fem / e, fem / en)}
        hot = rec["hot"]
        hr = None
        if hot["label"] in labs:
            j = labs.index(hot["label"])
            pt = pic["lat"][j][1]
            uh = p79._u3_hist_at(pt)
            hists = {"U3": uh, "eddy": np.asarray(pic["L"][:, j]), "end": np.asarray(pic["SL"][:, j])}
            fem = hot["th"]["val"] / fac
            d0 = float(cem43(T_GRID, T_INF * np.ones_like(T_GRID)))
            pics = {}
            for k, hh in hists.items():
                ra = running_average(T_GRID, fac * hh, window=float(W_B2))
                dk = float(cem43(T_GRID, T_INF + fac * hh))
                pics[k] = {"per_W": float(hh[-1]), "ratio": fem / float(hh[-1]), "cem43": dk,
                           "excess": dk - d0, "ra_max": float(ra[1]), "ra_tmax": float(ra[2])}
            hr = {"label": hot["label"], "box_pt": pt, "idx": j, "s": fac,
                  "theta_acc": hot["th"]["val"], "fem_per_W": fem,
                  "cem43_tok": hot["cem43"]["tok"], "cem43": hot["cem43"]["val"], "D0": d0,
                  "excess": hot["cem43"]["val"] - d0,
                  "ra_max": hot["ra_max"]["val"], "ra_tmax": hot["ra_tmax"]["val"],
                  "pictures": pics}
        out[s] = {"rows": rows, "hot": hr}
    return out


def _knob_parsed():
    path, parsed = c33._on_root(p79._knob)
    assert parsed is not None, f"{p79.KNOB} unset: P76 needs entry 59's log"
    return path, parsed


# ------------------------------------------------------------------ tests
def test_custody_k59_ks_knorm_ktwin_r52():
    gates = []
    path, parsed = _knob_parsed()
    _say(f"MAT21-43b knob {p79.KNOB}: {path!r}")
    # (K-59)
    bad, n = p79._reformat_misses(parsed)
    gates.append((("K-59 reformat", bad, n), not bad and n == N_VALUES, n))
    lines = {s: r["_lines"] for s, r in parsed.items()}
    gates.append((("K-59 lines", lines), lines == WANT_LINES, 0))
    hr = p79._hot_row(parsed)  # raises if the printed c* is > 5e-7 m off its lattice point
    lat = dict(c6._lattice())
    lab = parsed["ccw"]["hot"]["label"]
    dpt = max(abs(a - b) for a, b in zip(lat.get(lab, (np.inf,) * 3), parsed["ccw"]["hot"]["point"]))
    ok_star = lab == WANT_STAR and lab in lat and hr["ccw"]["row"] is not None and dpt <= 5e-7
    gates.append((("K-59 c*", lab, dpt), ok_star, dpt))
    _say(f"MAT21-43b (K-59) {n} values re-formatted (want {N_VALUES}), misses {bad or 'none'}; "
         f"block lines {lines} (want {WANT_LINES}): {lines == WANT_LINES}; ccw c* {lab!r} on C6's "
         f"lattice {lab in lat}, printed vs lattice {dpt:.3e} m (<= 5e-7) (ASSERTED)")
    # (K-s)
    for s, rec in parsed.items():
        pg, pa, f = rec["P_gel"], rec["P_acc"], rec["factor"]
        miss = abs(pg["val"] / pa["val"] / f["val"] - 1.0)
        bnd = sum(c33._half_ulp_rel(t["tok"]) for t in (pg, pa, f))
        gates.append((("K-s", s, miss, bnd), miss <= bnd, miss))
        _say(f"MAT21-43b (K-s) {s}: P_gel {pg['tok']} / P_acc {pa['tok']} / factor {f['tok']} - 1 "
             f"= {miss:.3e} (half-ulp sum {bnd:.3e}, ASSERTED)")
    # (K-norm) and (R52)
    page = c33._load_page()
    e59 = c33._on_root(page.parse_c6_log, path)
    e52 = c33._on_root(page.parse_c6_log, page.SOURCE_LOG)
    nm = e59["normalisers"]
    for key, tok_key in (("P_gel", "p_box"), ("P_acc", "p_acc")):
        mine = format(parsed["ccw"][key]["val"], ".9e")
        ok = mine == nm[tok_key]
        gates.append((("K-norm", key, mine, nm[tok_key]), ok, 0))
        _say(f"MAT21-43b (K-norm) ccw {key} {parsed['ccw'][key]['tok']} as .9e {mine} vs entry 59's "
             f"C6a normalisers {tok_key} {nm[tok_key]} (:{nm['line']}): {ok} (ASSERTED)")
    f59, f52 = _r52_fields(e59), _r52_fields(e52)
    diff = [k for k in f52 if f59[k] != f52[k]]
    gates.append((("R52", diff), not diff and len(f59["lattice"]) == 18, diff))
    _say(f"MAT21-43b (R52) parse_c6_log({e59['log']}) vs SOURCE_LOG {e52['log']}: lattice "
         f"{len(f59['lattice'])} rows (label, b1p_ccw, sar10g, mass_kg), c* {f59['c_star']['label']!r}, "
         f"normalisers, S2: differing fields {diff or 'none'} (strings, ASSERTED)")
    # (K-twin)
    twin, gp = c33._on_root(_parse_twin, path)
    worst = 0.0
    for s, rec in parsed.items():
        g = gp.get(s, [])
        gp_ok = len(g) == 1 and g[0][2] == "1.000000000000e+00" and g[0][1] == rec["P_gel"]["tok"]
        gates.append((("K-twin G*P", s, g), gp_ok, 0))
        _say(f"MAT21-43b (K-twin) {s}: twin G*P token {g[0][2] if g else None} (want "
             f"1.000000000000e+00), twin P_gel {g[0][1] if g else None} vs P61 {rec['P_gel']['tok']}: "
             f"{gp_ok} (strings, ASSERTED)")
        f = rec["factor"]
        tw = twin.get(s, {})
        for lab in PROBE_LABELS:
            if lab not in tw:
                gates.append((("K-twin missing", s, lab), False, None))
                continue
            th = rec["probes"][lab]
            miss = abs(th["val"] / (f["val"] * float(tw[lab])) - 1.0)
            bnd = sum(c33._half_ulp_rel(t) for t in (th["tok"], f["tok"], tw[lab]))
            worst = max(worst, miss / bnd)
            gates.append((("K-twin", s, lab, miss, bnd), miss <= bnd, miss))
            _say(f"MAT21-43b (K-twin) {s} {lab}: theta_acc {th['tok']} / ({f['tok']} x twin "
                 f"{tw[lab]}) - 1 = {miss:.3e} (half-ulp sum {bnd:.3e}, ASSERTED)")
    _say(f"MAT21-43b (K-twin) worst miss / bound {worst:.3f} over {2 * len(PROBE_LABELS)} rows")
    bad_g = [g for g in gates if not g[1]]
    assert not bad_g, bad_g


def test_id3_synthetic_on_eddy_and_end():
    gates = []
    pic = _pictures()
    for col, key, k in (("eddy", "E", 1), ("end", "S", 2)):
        syn = p79._synthetic()
        s = syn["cw"]["factor"]["val"]
        for i, lab in enumerate(PROBE_LABELS):
            v = s * float(pic[key][-1, i])
            syn["cw"]["probes"][lab] = {"tok": repr(v), "val": v}
        r = _rows(syn, pic)["cw"]["rows"]
        w = max(abs(r[lab]["ratio"][k] - 1.0) for lab in PROBE_LABELS)
        gates.append((("Id-3", col, w), w <= B1_BAND, w))
        _say(f"MAT21-43b (Id-3) synthetic s {s:.9e}, theta_acc planted s x {col}/P at the six points: "
             f"max |FEM/{col} - 1| {w:.3e} (band {B1_BAND:g}, ASSERTED)")
    bad = [g for g in gates if not g[1]]
    assert not bad, bad


def test_three_pictures_record_and_json():
    """The record (no band): asserts only that it was formed and written."""
    path, parsed = _knob_parsed()
    pic = _pictures()
    t = _rows(parsed, pic)
    _say(f"MAT21-43b per deposited gel watt P_dep {P_DEP:.6f} W (all three pictures); entry 59 per "
         f"accepted watt: P_gel ccw {parsed['ccw']['P_gel']['tok']} W, cw {parsed['cw']['P_gel']['tok']}"
         f" W; theta per deposited watt = theta_acc / s, s = P_gel/P_acc of the same sense")
    _say("MAT21-43b TABLE point | sense | s | theta_acc K | theta_acc/s K/W | U3/P | eddy/P | end/P | "
         "FEM/U3 | FEM/eddy | FEM/end")
    misses = []
    for s in ("ccw", "cw"):
        for lab, r in t[s]["rows"].items():
            ru, re_, rn = r["ratio"]
            _say(f"MAT21-43b TABLE {lab} | {s} | {r['s']:.9e} | {r['theta_acc']:+.7g} | "
                 f"{r['fem_per_W']:.7g} | {r['U3_per_W']:.7g} | {r['eddy_per_W']:.7g} | "
                 f"{r['end_per_W']:.7g} | {ru:.7g} | {re_:.7g} | {rn:.7g}")
            if s == "ccw":
                for got, want, nm in zip(r["ratio"], EXPECT_CCW[lab], ("U3", "eddy", "end")):
                    rel = abs(got / want - 1.0)
                    misses.append((lab, nm, got, want, rel))
        h = t[s]["hot"]
        pu, pe, pn = (h["pictures"][k] for k in ("U3", "eddy", "end"))
        _say(f"MAT21-43b TABLE hotspot c* {h['label']} | {s} | {h['s']:.9e} | {h['theta_acc']:+.7g} | "
             f"{h['fem_per_W']:.7g} | {pu['per_W']:.7g} (u by _u3_hist_at) | {pe['per_W']:.7g} "
             f"(L[-1, {h['idx']}]) | {pn['per_W']:.7g} (SL[-1, {h['idx']}]) | {pu['ratio']:.7g} | "
             f"{pe['ratio']:.7g} | {pn['ratio']:.7g}")
        if s == "ccw":
            for k, want, nm in zip(("U3", "eddy", "end"), EXPECT_CCW["hotspot"], ("U3", "eddy", "end")):
                got = h["pictures"][k]["ratio"]
                misses.append(("hotspot", nm, got, want, abs(got / want - 1.0)))
    # Reproduction of the brief's expected ccw table (4 significant figures, so its own
    # rounding is up to 5e-4 relative; a miss > 1e-4 is disclosed, never red).
    over = [(lab, nm, f"{got:.7g}", want, f"{rel:.2e}") for lab, nm, got, want, rel in misses
            if rel > 1e-4]
    _say(f"MAT21-43b expected ccw table (brief, 4 s.f.): {len(misses)} cells, max rel miss "
         f"{max(m[4] for m in misses):.3e}; cells above 1e-4 (disclosed, rounding of the brief's "
         f"4 s.f. print reaches 5e-4): {over or 'none'}")
    # The hotspot dose
    for s in ("ccw", "cw"):
        h = t[s]["hot"]
        _say(f"MAT21-43b DOSE {s} c* {h['label']}: entry 59 cem43 per P_acc watt {h['cem43_tok']} s; "
             f"D0 = cem43(T_GRID, T_INF) {h['D0']:.15e} s (cited {D0_CITED:.15e}, rel "
             f"{abs(h['D0'] / D0_CITED - 1.0):.2e}); excess {h['excess']:.6g} s (+-5e-21 at the "
             f".9e print), excess/D0 {h['excess'] / h['D0']:.5g}; at s {h['s']:.9e}: " + "; ".join(
                 f"{k} cem43 {p['cem43']:.9e} excess {p['excess']:.6g} s ({p['excess'] / h['D0']:.5g})"
                 for k, p in h["pictures"].items())
             + f" (DOSE_BAR {DOSE_BAR:g} unused; printed only)")
    _say("MAT21-43b DOSE note: a cem43 per accepted watt is approximately the unheated baseline "
         "D0 and is not a scalable per-watt dose (cem43 is not linear in power); the excess is "
         "printed as a record, never as a dose lever")
    # The six-minute window
    for s in ("ccw", "cw"):
        h = t[s]["hot"]
        still = h["ra_tmax"] == float(T_GRID[-1])
        _say(f"MAT21-43b WINDOW {s} c*: entry 59 running mean (W {W_B2:g} s) max {h['ra_max']:.9e} K "
             f"at t {h['ra_tmax']:.1f} s (t* {float(T_GRID[-1]):.1f} s; still rising {still}); at s: "
             + "; ".join(f"{k} max {p['ra_max']:.7g} K at t {p['ra_tmax']:.1f} s"
                         for k, p in h["pictures"].items())
             + f" (AVG_BAR {AVG_BAR:g} unused; printed only)")
    # Notes
    rc = t["ccw"]["rows"]
    hc = t["ccw"]["hot"]
    _say(f"MAT21-43b NOTE c* is C6's ccw SAR10g lattice argmax, not the thermal maximum: per "
         f"accepted watt, ccw, P1 +y / P2 -y read {rc['P1 +y wall']['theta_acc'] / hc['theta_acc']:.4g}x"
         f" / {rc['P2 -y wall']['theta_acc'] / hc['theta_acc']:.4g}x c*'s theta "
         f"({rc['P1 +y wall']['theta_acc']:.9e} and {rc['P2 -y wall']['theta_acc']:.9e} over "
         f"{hc['theta_acc']:.9e} K); the cw arm is heated at ccw's c*")
    _say("MAT21-43b NOTE probe labels are P37's box frame; the printed coordinates are the "
         "fixture's (P1 +y wall prints at (0.21, 0.0, 0.0)); rows are matched by label, never by "
         "coordinate")
    doc = {"entry59_log": path, "P_dep_W": P_DEP, "D0_s": t["ccw"]["hot"]["D0"],
           "normaliser": "pictures per deposited gel watt; entry 59 per accepted watt, s = P_gel/P_acc",
           "senses": t}
    if MPI.COMM_WORLD.rank == 0:
        JSON_PATH.parent.mkdir(parents=True, exist_ok=True)
        JSON_PATH.write_text(json.dumps(doc, indent=1, default=list))
    MPI.COMM_WORLD.barrier()
    back = json.loads(JSON_PATH.read_text())
    _say(f"MAT21-43b JSON {JSON_PATH} keys {sorted(back)}; senses {sorted(back['senses'])}")
    assert back["entry59_log"] == path and sorted(back["senses"]) == ["ccw", "cw"]
    assert all(np.isfinite(r["ratio"]).all() for s in t for r in t[s]["rows"].values())
