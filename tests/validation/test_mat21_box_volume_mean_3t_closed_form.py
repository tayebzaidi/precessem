r"""`MAT-21` step 52 — the box volume mean ⟨θ⟩(t*) of the six-face Robin box on entry 58's
3 T degree-1 lattice, closed form (§9 item 1008, §10 P92; brief `docs/planning/chunks/MAT-21.md`
"§9 item 1008").

P91′'s series and identities (``test_mat21_box_volume_mean_closed_form.py``, imported as
``p91``: ``lattice_split`` / ``mean_lattice`` / ``mean_uniform_split``), re-run on the 3 T
weights ``p89._data()["w"]["3T"]`` (entry 58's 18 ``sar10g`` strings, normalised).  Not a knob on
P91′: that module and its node IDs are untouched.  No ``src/`` edit, no landed module edited,
C6 import-only, never ``p89._gmat``, no fixture or ``test_*`` name imported.  Sets no band and
moves none (``B1_BAND`` imported).  Numpy / scipy; complex build because P89′ imports C33 / C6.

Asserted: (K) entry 60's 18 labels and ``sar10g`` strings, parsed here through
``page.parse_c6_log``, equal entry 58's by string (P86's (R58) re-asserted); (K-91) P91′'s
records re-computed here equal W1b ``:78`` / ``:79``, read via ``read_log_text(errors="strict")``
with exactly one distinct match each; (Lim-∞), (Bq-N) split N₀ → 2N₀, (Bq-split) split(N₀) vs
direct(2N₀), (Sym) a → −a and b → −b per-ball invariance, (N-eq-m) the 1.5 T weights at
mirrored labels, all at ``B1_BAND``; (N-h) at h = 1e-11 on N_h = (40, 40, 40) (below the bracket
flip; never ``m_short`` there) ≤ 1e-12; negative control (N-lat) |⟨θ⟩₃/⟨θ⟩₁.₅ − 1| > ``B1_BAND``;
(Env) min ⟨θ_i⟩ ≤ ⟨θ⟩₃ ≤ max ⟨θ_i⟩.  Printed (🟡 record, ≤ 7 digits): ⟨θ⟩_lat,3T(t*), 3 T / 1.5 T,
lattice / uniform, surface fraction, uniform, the per-accepted-watt lever (context), per-ball
table.

Not claimed: entry 60's or 70's FEM ⟨θ⟩; a 3 T degree-2 lattice; a 3 T gel claim; anything about
the near-face SAR.
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np

from tests.log_records import read_log_text
from tests.validation import test_geo38_c39_thermal_reader as c46
from tests.validation import test_mat21_box_volume_mean_closed_form as p91
from tests.validation import test_mat21_lattice_ball_superposition_closed_form as p89
from tests.validation.test_imp4_pec_cylinder_closed_form import B1_BAND
from tests.validation.test_mat21_box_closed_form import RHO_C, T_C
from tests.validation.test_mat21_box_dose_closed_form import T_GRID

PFX = "MAT21-52"
REPO = Path(__file__).resolve().parents[2]
W1B_LOG = REPO / "docs" / "testing" / "logs" / "20261116T030255Z_MAT-21-P91-W1b.log"
RE_W1B = {
    "uniform": (re.compile(r"MAT21-51 RECORD <theta>_uniform\(t\*\) (\S+) K/W"), 78, "0.02521149"),
    "lattice": (re.compile(r"MAT21-51 RECORD <theta>_lattice 1\.5 T deg 1\(t\*\) (\S+) K/W"), 79,
                "0.02713084"),
}
N_L = p91.N_L
N2 = tuple(2 * v for v in N_L)
N_H = p91.N_H
H_SMALL = p91.H_SMALL
N_HCTRL = p91.N_HCTRL
V_BOX = p91.V_BOX
ENV_PRED = 2.252e-4
# Predictions (brief arithmetic), printed only.
PRED_L3 = (0.02713068, 0.0271267, 0.0271328)
PRED_Q = (0.9999941, 0.99977, 1.00023)
PRED_LU = 1.076124
PRED_SURF = 0.012443
PRED_NLAT = 5.95e-6
PRED_BQS = 4.8e-13
PRED_NH = (0.0, 1.2e-13)
PRED_LEVER = 3.69348
C39_LEVER = "4.042515"
P89_LEVER = "4.078321"


def _say(msg):
    print(msg, flush=True)


def _fail(gates):
    return [g for g in gates if not g[1]]


def _in(v, lo, hi):
    return lo <= v <= hi


def _ts():
    return float(T_GRID[-1])


def _w(side):
    return np.asarray(p89._data()["w"][side], dtype=float)


def _pts():
    return np.asarray(p89._data()["pts"], dtype=float)


def _per_ball(w, pts, t):
    """Per-ball box means ⟨θ_i⟩(t) per watt of ball i, through P91′'s split routes."""
    _, mod_per = p91.mean_lattice(w, pts, t, n=N_L, s0=p89.S0)
    return p89.S0 / (RHO_C * V_BOX) + mod_per


# ------------------------------------------------------------------ tests
def test_preconditions_and_k():
    gates = []
    d = p89._data()
    w3, pts = _w("3T"), _pts()
    gates.append((("n", len(w3)), len(w3) == 18 and len(pts) == 18, len(w3)))
    s = float(w3.sum())
    gates.append((("sum w3", s), abs(s - 1.0) <= 4 * p91.EPS, s - 1.0))
    mir = d["mirror"]
    gates.append((("mirror complete", mir), all(m is not None for m in mir), 0))
    _say(f"{PFX} lattice {len(w3)} balls (3T, entry 58 weights), |sum w - 1| {abs(s - 1.0):.3e} "
         f"(<= 4 eps ASSERTED); mirror map complete {all(m is not None for m in mir)} (ASSERTED)")
    # (K) entry 60 parsed here vs entry 58 (p89's parse), by string, looked up by label
    hit = c46.find_log(c46.ENTRY60_GLOB)
    gates.append((("entry 60 found", hit), hit is not None, 0))
    assert hit is not None, gates
    e60 = d["page"].parse_c6_log(hit[0])
    e58 = d["parsed"]["3T"]
    by60 = {r["label"]: r["sar10g"] for r in e60["lattice"]}
    by58 = {r["label"]: r["sar10g"] for r in e58["lattice"]}
    lab_ok = (len(e60["lattice"]) == 18 and len(e58["lattice"]) == 18 and len(by60) == 18
              and set(by60) == set(by58))
    diff = sorted(lab for lab in by58 if by60.get(lab) != by58[lab])
    gates.append((("K labels", sorted(by60)), lab_ok, 0))
    gates.append((("K sar10g", diff), not diff, diff))
    _say(f"{PFX} (K) parse_c6_log({Path(str(hit[0])).name} [{hit[1]}]) vs entry 58 "
         f"{Path(str(d['hit'][0])).name}: {len(by60)} / {len(by58)} labels, label sets equal "
         f"{set(by60) == set(by58)}; sar10g strings differing {diff or 'none'} (strings, ASSERTED)")
    bad = _fail(gates)
    assert not bad, bad


def test_k91_custody():
    gates = []
    ts = _ts()
    u, _, _ = p91.mean_uniform_split(ts)
    l15, _, _ = p91.lattice_split(_w("1.5T"), _pts(), ts)
    text = read_log_text(W1B_LOG, errors="strict")
    lines = text.splitlines()
    for key, mine in (("uniform", u), ("lattice", l15)):
        rx, want_line, want = RE_W1B[key]
        hits = [(i + 1, m.group(1)) for i, ln in enumerate(lines) for m in [rx.search(ln)] if m]
        distinct = sorted({v for _i, v in hits})
        tok = f"{mine:.7g}"
        ok = (len(distinct) == 1 and distinct[0] == tok == want
              and want_line in [i for i, _v in hits])
        gates.append((("K-91", key, tok, distinct, hits), ok, 0))
        _say(f"{PFX} (K-91) {key}: re-computed {tok} vs W1b {W1B_LOG.name} matches "
             f"{len(hits)} (lines {[i for i, _v in hits]}), distinct {distinct} (brief "
             f"{want} at :{want_line}) (ASSERTED)")
    bad = _fail(gates)
    assert not bad, bad


def test_identities():
    gates = []
    d = p89._data()
    w3, w15, pts = _w("3T"), _w("1.5T"), _pts()
    ts = _ts()
    # (Lim-inf) on the direct form
    tl = 1000.0 * T_C
    a, _ = p91.mean_lattice(w3, pts, tl)
    b, _ = p91.mean_lattice(w3, pts, None, steady=True)
    rel = abs(a - b) / abs(b)
    gates.append((("Lim-inf", a, b), rel <= B1_BAND, rel))
    _say(f"{PFX} (Lim-inf) 3T <theta>(1000 T_C) {a:.7g} vs steady {b:.7g}: rel {rel:.3e} "
         f"(band {B1_BAND:g}, ASSERTED)")
    # (Bq-N) split N0 -> 2N0
    l0, free, lmod = p91.lattice_split(w3, pts, ts)
    l1, _, _ = p91.lattice_split(w3, pts, ts, n=N2)
    rel = abs(l1 - l0) / abs(l0)
    gates.append((("Bq-N", l0, l1), rel <= B1_BAND, rel))
    _say(f"{PFX} (Bq-N) 3T split N {N_L} -> 2N {N2}: {l0:.7g} -> {l1:.7g} K/W, rel {rel:.3e} "
         f"(band {B1_BAND:g}, ASSERTED; predicted <= 1e-15)")
    # (Bq-split) split(N0) vs direct(2N0) — never direct(N0)
    ld, _ = p91.mean_lattice(w3, pts, ts, n=N2)
    rel = abs(l0 - ld) / abs(ld)
    gates.append((("Bq-split", l0, ld), rel <= B1_BAND, rel))
    _say(f"{PFX} (Bq-split) s0 {p89.S0:g} s: free {free:.7g} + modal [s0, t*] {lmod:.7g} = "
         f"{l0:.7g} (N0) vs direct {ld:.7g} (2N0 {N2}): rel {rel:.3e} (band {B1_BAND:g}, "
         f"ASSERTED; predicted approx {PRED_BQS:g})")
    # (Sym) per-ball means under a -> -a and b -> -b separately
    per = _per_ball(w3, pts, ts)
    for ax, tag in ((0, "a"), (1, "b")):
        fl = pts.copy()
        fl[:, ax] = -fl[:, ax]
        perf = _per_ball(w3, fl, ts)
        res = float(np.max(np.abs(perf - per) / np.abs(per)))
        gates.append((("Sym", tag, res), res <= B1_BAND, res))
        _say(f"{PFX} (Sym) {tag} -> -{tag}: max_i |<theta_i>' - <theta_i>| / <theta_i> {res:.3e} "
             f"(band {B1_BAND:g}, ASSERTED; predicted 0)")
    # (N-eq-m) 1.5 T weights at mirrored labels reproduce the 1.5 T value
    mir = np.asarray(d["mirror"], dtype=int)
    l15, _, _ = p91.lattice_split(w15, pts, ts)
    l15m, _, _ = p91.lattice_split(w15[mir], pts, ts)
    rel = abs(l15m - l15) / abs(l15)
    gates.append((("N-eq-m", l15, l15m), rel <= B1_BAND, rel))
    _say(f"{PFX} (N-eq-m) 1.5T weights at mirrored labels {l15m:.7g} vs 1.5T {l15:.7g}: rel "
         f"{rel:.3e} (band {B1_BAND:g}, ASSERTED)")
    bad = _fail(gates)
    assert not bad, bad


def test_negative_controls():
    gates = []
    w3, w15, pts = _w("3T"), _w("1.5T"), _pts()
    ts = _ts()
    stored = ts / (RHO_C * V_BOX)
    # (N-h) at h = 1e-11 on N_h, direct route only (never m_short at this h)
    lt, _ = p91.mean_lattice(w3, pts, ts, n=N_H, h=H_SMALL)
    rel = abs(lt / stored - 1.0)
    gates.append((("N-h", lt, stored), rel <= N_HCTRL, rel))
    _say(f"{PFX} (N-h) 3T h {H_SMALL:g}, N_h {N_H}: <theta>(t*) {lt:.7g} vs P_dep t*/(rho c V) "
         f"{stored:.7g}: |ratio - 1| {rel:.3e} (ceiling {N_HCTRL:g}, ASSERTED; predicted in "
         f"[{PRED_NH[0]:g}, {PRED_NH[1]:g}] inside {_in(rel, *PRED_NH)})")
    # (N-lat) 1.5 T lattice in place of 3 T
    l3, _, _ = p91.lattice_split(w3, pts, ts)
    l15, _, _ = p91.lattice_split(w15, pts, ts)
    dev = abs(l3 / l15 - 1.0)
    gates.append((("N-lat", dev), dev > B1_BAND, dev))
    _say(f"{PFX} (N-lat) |<theta>_3T / <theta>_1.5T - 1| {dev:.3e} (> {B1_BAND:g} ASSERTED; "
         f"predicted approx {PRED_NLAT:g}, {dev / B1_BAND:.3g} x the band)")
    # (Env) min_i <theta_i> <= <theta>_3 <= max_i <theta_i>
    per = _per_ball(w3, pts, ts)
    lo, hi = float(per.min()), float(per.max())
    ok = lo <= l3 <= hi
    gates.append((("Env", lo, l3, hi), ok, 0))
    _say(f"{PFX} (Env) min <theta_i> {lo:.7g} <= <theta>_3T {l3:.7g} <= max <theta_i> {hi:.7g}: "
         f"{ok} (ASSERTED); envelope (max - min) / <theta>_3T {(hi - lo) / l3:.4e} (brief "
         f"{ENV_PRED:g}); dot(w, per) vs split rel {abs(float(np.dot(w3, per)) - l3) / l3:.3e}")
    bad = _fail(gates)
    assert not bad, bad


def test_record_printed():
    d = p89._data()
    w3, w15, pts = _w("3T"), _w("1.5T"), _pts()
    ts = _ts()
    stored = ts / (RHO_C * V_BOX)
    u, _, _ = p91.mean_uniform_split(ts)
    l3, _, _ = p91.lattice_split(w3, pts, ts)
    l15, _, _ = p91.lattice_split(w15, pts, ts)
    q = l3 / l15
    lu = l3 / u
    surf = 1.0 - l3 / stored
    lever = float(p89.P83_S[2]) * q
    _say(f"{PFX} records from the s0-split routes (N0 {N_L}, s0 {p89.S0:g} s), per deposited watt")
    pv, lo, hi = PRED_L3
    _say(f"{PFX} RECORD <theta>_lat,3T deg 1(t*) {l3:.7g} K/W (predicted approx {pv}, range "
         f"[{lo}, {hi}] inside {_in(l3, lo, hi)})")
    pv, lo, hi = PRED_Q
    _say(f"{PFX} RECORD 3T / 1.5T {q:.7g} (predicted approx {pv}, range [{lo}, {hi}] inside "
         f"{_in(q, lo, hi)}); 1.5T {l15:.7g} K/W")
    _say(f"{PFX} RECORD lattice / uniform {lu:.7g} (predicted approx {PRED_LU})")
    _say(f"{PFX} RECORD surface fraction {surf:.6f} (predicted approx {PRED_SURF}); stored "
         f"fraction {l3 / stored:.6f}; P_dep t*/(rho c V) {stored:.6f} K/W")
    _say(f"{PFX} RECORD uniform {u:.7g} K/W (re-print, frequency-independent per deposited watt)")
    _say(f"{PFX} context: per-accepted-watt lever P83_S share {p89.P83_S[2]} x ratio = "
         f"{lever:.6g} (predicted approx {PRED_LEVER}), beside C39's {C39_LEVER} and P89's "
         f"{P89_LEVER} (context only, not asserted)")
    per = _per_ball(w3, pts, ts)
    labels = d["labels"]
    _say(f"{PFX} per-ball table (3T): label | w | <theta_i>(t*) K/W per ball watt")
    for lab, wi, pi in zip(labels, w3, per):
        _say(f"{PFX}   {lab:>18s} | {wi:.7g} | {pi:.7g}")
    assert l3 > 0.0 and u > 0.0
