"""GEO-38 C18 (§9 item 866, §10 chain C18): the degree-2 ``C_mode`` (k = 1) of
F-clinical-1.5T's copper 64×64 at 64 MHz and C13's class anatomy at degree 2
beside degree 1 — read from entry 54's saved npz
(``logs/geo38/c3-cu-s64-d2.npz``, gitignored scratch, read only).  Item 655's
(C15's) route on the 1.5 T degree-2 file.

Knob ``FEM_EM_GEO38_C18`` (C10's ``_knob``, imported): unset / empty ⇒ the
gated test skips; ``1`` ⇒ it runs and an absent input fails; anything else ⇒
``ValueError``.

Asserted:

* **(RR2)** reciprocity / σ_max of the loaded degree-2 64×64 string-equal entry
  54's prints (``20261107T080008Z_GEO-38-C5d.log:14813–14814``) and inside
  ``RECIPROCITY_BAND`` / ``MATRIX_PASSIVITY_CEILING``;
* **(CT2)** the 34 class lines equal ``C5d.log:14815–14848`` under C13's
  ``_norm``, 34 / 34.  Those lines are **not pasted here**: the C5d log is held
  by the private-leak guard (gitignored ``logs/leak-held/``) and a ``.9e`` mean
  can collide with the private reference set by chance, so the module reads
  them at run time from the held log, after a sha256 check (the brief's
  fallback);
* **(Im)** Im λ₁ > 0 (C10's assertion);
* **(Map)** the geometry map is unique 64 / 64, the name route agrees 64 / 64,
  P33 / P41 land on C10's slots;
* **(CD)** negative control: the degree-1 file through the same code gives C13's
  (RR) strings (≠ the degree-2 strings) and C10's ``C_mode`` 195.0265 pF;
* **(CN)** negative control: P33's z flipped ⇒ a class changes n and (CT2)
  no longer matches 34 / 34.

Printed only: ``C_mode`` (d2) beside C10's, the degree-2 root and root −
``C_mode``; the per-class order reading d1 / d2 and its SUMMARY; C13's abs
statistic at each order; the crossing spreads.  Nothing asserted on a ratio.
Writes ``/workspace/logs/geo38/c18-c-mode-d2.json`` (C10's keys) after every
anchor holds.  Not claimed: a root or tuning statement, mesh-vs-physics, the
known-issues 2026-10-27 retirement, any AED value.
"""

from __future__ import annotations

import hashlib
import json
import re
import time

import numpy as np
import pytest

from tests.package.test_wf16_load_result_f_clinical import C10_C_MODE_PF_RECORD
from tests.validation.test_geo38_c13_class_anatomy import (
    C3_TABLE_LINES,
    RECIPROCITY_RECORD_STR,
    SIGMA_MAX_RECORD_STR,
    _abs_spreads,
    _norm,
    _table_lines,
)
from tests.validation.test_geo38_f_clinical_c_mode import (
    C3_KEYS,
    FREQUENCY_HZ,
    PRED_C3_P33_P41_SLOTS,
    REFERENCE_IMPEDANCE_OHM,
    WORKSPACE,
    _knob,
    map_from_geometry,
    map_from_names_64,
    ring_c_modes,
)
from tests.validation.test_geo38_f_clinical_matrix import (
    N_AZIMUTH_CLASSES_32,
    N_SYMMETRY_CLASSES_32,
    SAVE_PATH,
    SAVE_PATH_D2,
    _azimuth_class_32,
    _class_spreads,
    _class_table,
    _expected_class_sizes,
)
from tests.validation.test_port_birdcage_lumped_column import ADJACENT_SPREAD_BAND
from tests.validation.test_port_birdcage_ring_column import _reciprocity_ratio
from tests.validation.test_port_birdcage_ring_matrix import MATRIX_PASSIVITY_CEILING
from tests.validation.test_port_lumped_sheet_sweep import RECIPROCITY_BAND

KNOB_C18 = "FEM_EM_GEO38_C18"
C18_JSON = WORKSPACE / "logs" / "geo38" / "c18-c-mode-d2.json"

# ---- records (not bands), each with its log line ------------------------------
# Entry 54's own prints: 20261107T080008Z_GEO-38-C5d.log:14813–14814 (log held).
RECIPROCITY_D2_RECORD_STR = "1.488791e-13"
SIGMA_MAX_D2_RECORD_STR = "0.999932072111"

# (CT2)'s record, read at run time (never copied into a tracked file).
C5D_HELD_LOG = WORKSPACE / "logs" / "leak-held" / "20261107T080008Z_GEO-38-C5d.log"
C5D_HELD_SHA256 = "260ab8980279c09f048fbf4120729e454293be71ec3416b6216dbadbdab2ec84"
C5D_TABLE_SPAN = (14815, 14848)  # 1-indexed, inclusive
C3_PREFIX = "[GEO-38 C3] "

# Degree-2 root, entry 56: P33 (1, +1j) 185.762904 pF (xl-ledger.md 2026-11-05 row,
# 20261105T080008Z_GEO-38-C5b-prime.log:12814); P41 (1, -1j) 185.820370 pF (:12817);
# P41 (1, +1j) and P33 (1, -1j) "none" in the grid (:12815–12816).
ROOT_D2_PF = 185.762904
CROSSINGS_D2_PF = {"P33 (1,+1j)": 185.762904, "P41 (1,-1j)": 185.820370,
                   "P41 (1,+1j)": None, "P33 (1,-1j)": None}
# Degree-1 at 0.25 pF, entry 51 (xl-ledger.md 2026-10-28 row,
# 20261028T070008Z_GEO-38-C4r.log:13260, :13263; :13261–13262 "do not cross").
CROSSINGS_D1_PF = {"P33 (1,+1j)": 193.599624, "P41 (1,-1j)": 193.716245,
                   "P41 (1,+1j)": None, "P33 (1,-1j)": None}
# Root − C_mode at degree 1 (C26, 20261103T103936Z_GEO-38-C26-W1.log).
ROOT_MINUS_CMODE_D1_PCT = -0.73
# C13's statistic at degree 1, 20261025T093429Z_GEO-38-C13-W1.log:237.
C13_OTHER1_OVER_SAME1_D1 = 0.8883

_SPREAD_RE = re.compile(r"spread\s+([0-9.]+)%")


def c18_enabled():
    return _knob(KNOB_C18)


def _say(msg):
    print(f"[GEO-38 C18] {msg}", flush=True)


def _rel_spread(line):
    m = _SPREAD_RE.search(line)
    if m is None:
        raise ValueError(f"no spread in {line!r}")
    return float(m.group(1)) / 100.0


def _held_table_lines(path=C5D_HELD_LOG):
    raw = path.read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    lines = raw.decode("utf-8").splitlines()
    lo, hi = C5D_TABLE_SPAN
    span = lines[lo - 1:hi]
    out = []
    for ln in span:
        if not ln.startswith(C3_PREFIX):
            raise ValueError(f"held line lacks the C3 prefix: {ln[:40]!r}")
        out.append(ln[len(C3_PREFIX):])
    return digest, out


def _ct_misses(lines, ref):
    if len(lines) != len(ref):
        return [("<length>", f"{len(lines)} != {len(ref)}")]
    return [(a, b) for a, b in zip(lines, ref) if _norm(a) != _norm(b)]


def _spread_of(crossings):
    vals = [v for v in crossings.values() if v is not None]
    return (max(vals) - min(vals), len(vals)) if len(vals) >= 2 else (None, len(vals))


def _load(path):
    assert path.exists(), f"{path} absent; never re-solve"
    with np.load(path, allow_pickle=False) as z:
        sv = {k: z[k] for k in z.files}
    missing = [k for k in C3_KEYS if k not in sv]
    assert not missing, f"{path.name} lacks {missing}"
    return sv


def _read(path):
    """C10's reader body (`test_geo38_f_clinical_c_mode.py:440–474`), returned."""
    sv = _load(path)
    s = np.array(sv["s64"], dtype=np.complex128)
    order = [str(x) for x in sv["order"]]
    az = np.asarray(sv["azimuth_deg"], dtype=float)
    zc = np.asarray(sv["z_m"], dtype=float)
    geo = map_from_geometry(az, zc)
    agree = sum(g == m for g, m in zip(geo, map_from_names_64(order)))
    ratio = _reciprocity_ratio(s)
    sigma_max = float(np.linalg.norm(s, 2))
    r = ring_c_modes(s, geo, freq_hz=FREQUENCY_HZ, n_slots=32, z0=REFERENCE_IMPEDANCE_OHM)
    return dict(sv=sv, s=s, order=order, az=az, zc=zc, geo=geo, agree=agree,
                ratio=ratio, sigma_max=sigma_max, r=r,
                rr=(f"{ratio:.6e}", f"{sigma_max:.12f}"))


# --------------------------------------------------------------------------- W0
def test_c18_knob_and_records(monkeypatch):
    """Ungated: the knob, the degree-1 record parse, the crossing-spread helper."""
    monkeypatch.delenv(KNOB_C18, raising=False)
    assert c18_enabled() is False
    monkeypatch.setenv(KNOB_C18, "1")
    assert c18_enabled() is True
    for bad in ("0", "yes", "build"):
        monkeypatch.setenv(KNOB_C18, bad)
        with pytest.raises(ValueError):
            c18_enabled()
    monkeypatch.delenv(KNOB_C18, raising=False)
    rel = [_rel_spread(ln) for ln in C3_TABLE_LINES]
    assert len(rel) == N_SYMMETRY_CLASSES_32
    assert max(rel) == pytest.approx(0.009262, abs=1e-12)  # C3.log:14768's worst
    assert (RECIPROCITY_D2_RECORD_STR, SIGMA_MAX_D2_RECORD_STR) != (
        RECIPROCITY_RECORD_STR, SIGMA_MAX_RECORD_STR)
    sp, n = _spread_of(CROSSINGS_D2_PF)
    assert n == 2 and sp == pytest.approx(185.820370 - 185.762904, abs=1e-9)
    with pytest.raises(ValueError):
        _rel_spread("no spread token")


# --------------------------------------------------------------------------- W1
def test_geo38_c18_degree2_c_mode():
    if not c18_enabled():
        pytest.skip(f"{KNOB_C18} unset (reads gitignored logs/geo38 npz + held log)")
    t0 = time.perf_counter()
    path = WORKSPACE / SAVE_PATH_D2
    d2 = _read(path)
    s, order, az, zc, geo, r = d2["s"], d2["order"], d2["az"], d2["zc"], d2["geo"], d2["r"]
    sv = d2["sv"]
    _say(f"{path} ({path.stat().st_size} B); s64 {s.shape} {s.dtype}; order {order[0]}…{order[-1]}; "
         f"cells {int(sv['cells'])}, n_ranks {int(sv['n_ranks'])}")
    assert s.shape == (64, 64), s.shape
    assert order == [f"P{i}" for i in range(33, 97)], order[:4]

    # ---- (RR2) ----
    rr = d2["rr"]
    _say(f"(RR2) reciprocity {rr[0]} vs entry 54 {RECIPROCITY_D2_RECORD_STR} (band {RECIPROCITY_BAND:.0e}); "
         f"sigma_max {rr[1]} vs entry 54 {SIGMA_MAX_D2_RECORD_STR} (ceiling {MATRIX_PASSIVITY_CEILING}) "
         f"(string equality, ASSERTED)")

    # ---- (CT2) against the held log ----
    digest, c5d_lines = _held_table_lines()
    _say(f"(CT2) record read at run time from {C5D_HELD_LOG.name} lines {C5D_TABLE_SPAN[0]}–"
         f"{C5D_TABLE_SPAN[1]}; sha256 {digest} (expected {C5D_HELD_SHA256}, ASSERTED)")
    classes, res = _class_table(az, zc, classify=_azimuth_class_32)
    sizes = {k: len(v) for k, v in classes.items()}
    spreads = _class_spreads(s, classes)
    mine = _table_lines(spreads, N_AZIMUTH_CLASSES_32)
    ct_miss = _ct_misses(mine, c5d_lines)
    _say(f"(CT2) classes {len(classes)} (expected {N_SYMMETRY_CLASSES_32}), sizes equal expected "
         f"{sizes == _expected_class_sizes()}, worst pair residual {res:.3e} deg; lines matching "
         f"C5d.log:14815–14848 {len(c5d_lines) - len(ct_miss)}/{len(c5d_lines)} (C13's _norm; ASSERTED)")
    for a, b in ct_miss:
        _say(f"(CT2) MISS mine '{_norm(a)}'")

    # ---- (Map) ----
    p33_p41 = (geo[0], geo[8])
    _say(f"(Map) geometry map unique {len(set(geo))}/64; name-route agreement {d2['agree']}/64; "
         f"P33 -> {geo[0]}, P41 -> {geo[8]} (C10's {PRED_C3_P33_P41_SLOTS}; ASSERTED)")

    # ---- C_mode at degree 2 ----
    lam1 = r["modes"][1]["lambda_anti"]
    for k in range(17):
        m = r["modes"][k]
        lam = m["lambda_anti"]
        _say(f"  k={k:2d} lambda_anti {lam.real:+.9e} {lam.imag:+.9e}j Ohm anti {m['anti']:.6f} "
             f"C {m['c_anti_f'] * 1e12:+.4f} pF")
    c_pf = r["c_mode_f"] * 1e12
    _say(f"C_mode(d2) {r['c_mode_f']:.12e} F = {c_pf:.4f} pF at {FREQUENCY_HZ / 1e6:g} MHz; C_2 "
         f"{r['c_k2_f']:.12e} F = {r['c_k2_f'] * 1e12:.4f} pF; anti(k=1) {r['modes'][1]['anti']:.6f}; "
         f"off-block {r['off_block']:.6e}; (Im) Im lambda_1 {lam1.imag:+.9e} Ohm (ASSERTED > 0)")

    # ---- (CD) the degree-1 file through the same route ----
    d1 = _read(WORKSPACE / SAVE_PATH)
    c1_str = f"{d1['r']['c_mode_f'] * 1e12:.4f}"
    _say(f"(CD) {SAVE_PATH.name}: reciprocity {d1['rr'][0]} sigma_max {d1['rr'][1]} vs C13 "
         f"{RECIPROCITY_RECORD_STR} / {SIGMA_MAX_RECORD_STR}; C_mode {c1_str} pF vs C10 "
         f"{C10_C_MODE_PF_RECORD} (string equality, ASSERTED)")

    # ---- (CN) P33's z flipped ----
    zf = zc.copy()
    zf[0] = -zf[0]
    cls_f, _r = _class_table(az, zf, classify=_azimuth_class_32)
    n_changed = sum(1 for k in classes if len(cls_f.get(k, [])) != len(classes[k]))
    lines_f = _table_lines(_class_spreads(s, cls_f), N_AZIMUTH_CLASSES_32) if len(cls_f) == len(classes) else []
    cn_miss = len(_ct_misses(lines_f, c5d_lines))
    _say(f"(CN) {order[0]} z sign flipped: classes {len(cls_f)}, classes with changed n {n_changed} "
         f"(ASSERTED >= 1), (CT2) lines differing {cn_miss} (ASSERTED > 0)")

    assert digest == C5D_HELD_SHA256, digest
    assert rr == (RECIPROCITY_D2_RECORD_STR, SIGMA_MAX_D2_RECORD_STR), rr
    assert d2["ratio"] <= RECIPROCITY_BAND, d2["ratio"]
    assert d2["sigma_max"] <= MATRIX_PASSIVITY_CEILING, d2["sigma_max"]
    assert len(classes) == N_SYMMETRY_CLASSES_32 and sizes == _expected_class_sizes(), sizes
    assert len(c5d_lines) == N_SYMMETRY_CLASSES_32 and not ct_miss, ct_miss[:3]
    assert len(geo) == 64 and len(set(geo)) == 64, sorted(geo)
    assert d2["agree"] == 64, d2["agree"]
    assert p33_p41 == PRED_C3_P33_P41_SLOTS, p33_p41
    assert lam1.imag > 0.0, lam1
    assert d1["rr"] == (RECIPROCITY_RECORD_STR, SIGMA_MAX_RECORD_STR), d1["rr"]
    assert c1_str == C10_C_MODE_PF_RECORD, c1_str
    assert n_changed >= 1 and cn_miss > 0, (n_changed, cn_miss)

    # ---- JSON (only after every anchor holds) ----
    C18_JSON.parent.mkdir(parents=True, exist_ok=True)
    C18_JSON.write_text(json.dumps({
        "source": str(path), "freq_hz": FREQUENCY_HZ, "z0_ohm": REFERENCE_IMPEDANCE_OHM,
        "n_slots": 32, "c_mode_f": r["c_mode_f"], "c_k2_f": r["c_k2_f"],
        "off_block": r["off_block"], "reciprocity": d2["ratio"], "sigma_max": d2["sigma_max"],
        "name_route_agreement": d2["agree"], "cells": int(sv["cells"]), "n_ranks": int(sv["n_ranks"]),
        "modes": {str(k): {"lambda_anti": [m["lambda_anti"].real, m["lambda_anti"].imag],
                           "anti": m["anti"], "c_anti_f": m["c_anti_f"]}
                  for k, m in r["modes"].items()},
    }, indent=2))
    _say(f"wrote {C18_JSON} ({C18_JSON.stat().st_size} B)")

    # ---- the record (printed only) ----
    c10 = float(C10_C_MODE_PF_RECORD)
    _say(f"C_mode d2 {c_pf:.4f} pF beside C10 d1 {C10_C_MODE_PF_RECORD} pF: move "
         f"{(c_pf / c10 - 1.0) * 100:+.4f}%")
    _say(f"root d2 {ROOT_D2_PF:.6f} pF (entry 56, P33 (1,+1j)); root - C_mode d2 "
         f"{(ROOT_D2_PF / c_pf - 1.0) * 100:+.4f}% beside d1 {ROOT_MINUS_CMODE_D1_PCT:+.2f}% (C26)")

    rel_d1 = [_rel_spread(ln) for ln in C3_TABLE_LINES]
    keys = [(same, st) for same in (True, False) for st in range(N_AZIMUTH_CLASSES_32)]
    ratios = {}
    for key, d1rel in zip(keys, rel_d1):
        d2rel = spreads[key][0]
        ratios[key] = d1rel / d2rel
        _say(f"order {'same ' if key[0] else 'other'} ring, {key[1]:2d} steps  rel d1 {d1rel * 100:8.4f}%  "
             f"rel d2 {d2rel * 100:8.4f}%  d1/d2 {ratios[key]:.4f}")
    w1 = keys[int(np.argmax(rel_d1))]
    w2 = max(keys, key=lambda k: spreads[k][0])
    vals = np.array(list(ratios.values()))

    def lbl(k):
        return f"{'same' if k[0] else 'other'} ring, {k[1]} step(s)"

    _say(f"SUMMARY order reading: worst d1 [{lbl(w1)}] {max(rel_d1) * 100:.4f}% ratio {ratios[w1]:.4f}; "
         f"worst d2 [{lbl(w2)}] {spreads[w2][0] * 100:.4f}% ratio {ratios[w2]:.4f}; ratio d1/d2 min "
         f"{vals.min():.4f} median {float(np.median(vals)):.4f} max {vals.max():.4f}; classes with ratio "
         f">= 2: {int((vals >= 2.0).sum())}/34; classes above ADJACENT_SPREAD_BAND {ADJACENT_SPREAD_BAND} "
         f"d1 {sum(v > ADJACENT_SPREAD_BAND for v in rel_d1)} d2 "
         f"{sum(spreads[k][0] > ADJACENT_SPREAD_BAND for k in keys)} (printed, not asserted)")

    abs2 = _abs_spreads(s, classes)
    abs1 = _abs_spreads(d1["s"], _class_table(d1["az"], d1["zc"], classify=_azimuth_class_32)[0])
    st1 = abs1[(False, 1)] / abs1[(True, 1)]
    st2 = abs2[(False, 1)] / abs2[(True, 1)]
    _say(f"C13 statistic (other ring 1 step abs / same ring 1 step abs): d1 {abs1[(False, 1)]:.6e} / "
         f"{abs1[(True, 1)]:.6e} = {st1:.4f} (C13 W1:237 {C13_OTHER1_OVER_SAME1_D1}); d2 "
         f"{abs2[(False, 1)]:.6e} / {abs2[(True, 1)]:.6e} = {st2:.4f}")

    for tag, cr in (("d2 (entry 56)", CROSSINGS_D2_PF), ("d1 at 0.25 pF (entry 51)", CROSSINGS_D1_PF)):
        sp, n = _spread_of(cr)
        txt = ", ".join(f"{k} {'none in grid' if v is None else f'{v:.6f} pF'}" for k, v in cr.items())
        _say(f"crossings {tag}: {txt}; spread over the {n} crossing(s) in grid "
             f"{'n/a' if sp is None else f'{sp:.6f} pF'} (the four-series spread is not available: "
             f"two series do not cross inside the grid)")
    _say(f"body {time.perf_counter() - t0:.2f} s")
