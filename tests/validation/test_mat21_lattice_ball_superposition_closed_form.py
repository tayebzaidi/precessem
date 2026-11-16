r"""`MAT-21` step 49 — the 18-point lattice as 10 g-ball sources, closed form (§9 item 963, §10
P89; brief `docs/planning/chunks/MAT-21.md` "§9 item 963", review note 2026-11-12 03:00).

Numpy / scipy (complex build only because P80 imports P79 → P61 / C6, and C33's reader uses
``mpi4py``).  Sets no band and moves none; every band is imported (``B1_BAND``).

The source.  Each of C6a's 18 lattice points i (``p80._lattice()``: labels in C6 order, box
points p_i = (a, b, c) with c = 0) carries a uniform heat source on the ball B_i of radius
r = ``c6._a10()`` about p_i, of strength q_i = w_i / V_B, w_i = SAR_i / Σ_j SAR_j — one
deposited watt in total.  SAR_i are each side's ccw SAR10g strings from the **tracked** C6a logs
(entry 52 at 1.5 T, entry 58 at 3 T) read through the ``ans:11`` page's ``parse_c6_log``
(``c33._load_page()``); C33's own JSON carries no lattice.  ρ, P_acc and the drive amplitude
cancel exactly from θ per deposited watt.

The closed form (derived here; analytical.py has no ball source in the Robin box and its box
series are even-only).  The gel box |a| ≤ A, |b| ≤ B, |c| ≤ C, ρc θ_t = k∇²θ + Σ q_i 1_{B_i},
−k∂_nθ = hθ on all six faces, θ(0) = 0, has

    θ(x, t) = (1/ρc) ∫_0^t ds ∫ G(x, y; s) q(y) dy,
    G(x, y; s) = Σ_modes Φ(x) Φ(y) / 𝒩 · e^{−αΛs},        α = k/ρc.

Split the time integral at s₀ ≤ t (Ewald-style):

  * s ∈ [0, s₀] — the walls are neglected and G is the free-space Gaussian of variance
    σ² = 2αs per axis.  Its integral over the ball B(p, r) is the probability P∞(d; s) that a
    3D Gaussian about x lies in the ball, d = |x − p|:
        d > 0:  P∞ = ½[erf((r−d)/(σ√2)) + erf((r+d)/(σ√2))]
                     + (σ/(d√(2π)))·[e^{−(r+d)²/(2σ²)} − e^{−(r−d)²/(2σ²)}],
        d = 0:  P∞ = erf(w) − (2/√π)·w·e^{−w²},  w = r/(σ√2)   (the limit form).
    The d > 0 bracket is evaluated as ½[erfc((d−r)/(σ√2)) − erfc((d+r)/(σ√2))] — the same
    function (erf(−z) = −erf z, erf = 1 − erfc), written so that a distant ball's ≈ 1e-50
    is not the difference of two unit numbers.  Every ball is ≥ δ = 0.031635 m from a face, so
    the first image carries ≤ exp(−(2δ)²/(4αs₀)) ≈ e^{−69.7} at s₀ = 100 s (printed).
    The piece is (q_i/ρc) ∫_0^{s₀} P∞(d_i; s) ds, by ``scipy.integrate.quad`` at
    ``epsabs=0, epsrel=1e-13, limit=200`` (``_quad_dose``'s settings).
  * s ∈ [s₀, t] — the mode sum.  A product mode Φ = φ_a(a) φ_b(b) φ_c(c) with ∇²Φ = −ΛΦ is a
    sum of plane waves e^{ik·y}, |k|² = Λ, and the ball mean of e^{ik·y} is e^{ik·p}·F(|k|r),
    F(u) = 3(sin u − u cos u)/u³ = 3 j₁(u)/u (P81's ``_ball_mean_end`` identity).  So
    ∫_{B_i} Φ(y) dy = V_B F(√Λ r) Φ(p_i) exactly, and
        θ₂ = (1/k) Σ_modes F(√Λ r) Φ(p_i) Φ(x)/𝒩 · w_i · [e^{−αΛs₀} − e^{−αΛt}]/Λ,
    the bracket written e^{−αΛs₀}·(−expm1(−αΛ(t − s₀))).  F is evaluated as
    3·``scipy.special.spherical_jn(1, u)``/u (the same F; scipy takes the Bessel route below
    u = 1, where the sin − u cos form loses 1–2 digits at u_min ≈ 0.24 — the brief's trap).

Modes.  Per axis a and b both families: even cos(γs/L), γ tan γ = Bi
(``robin_slab_eigenvalues``), 𝒩 = L(γ + sin γ cos γ)/γ; odd sin(γs/L), γ cot γ = −Bi
(``robin_slab_odd_eigenvalues``), 𝒩 = L(γ − sin γ cos γ)/γ; Bi = hL/k.  Axis c: even only —
every ball centre and both c* sit at c = 0, so every odd-c term carries sin 0 = 0, and the
c-factor is 1/𝒩_c for every ball.  Λ = Σ_d (γ_d/L_d)².  Truncation per axis at
α(γ/L)²s₀ ≥ 40: N = (180, 120, 30) per family at s₀ = 100 s; (250, 160, 40) at s₀ = 50 s.
Contraction: G_lm = Σ_n T_lmn/𝒩_c,n, T = F·e^{−αΛs₀}·(−expm1(−αΛ(t−s₀)))/Λ; ball i at x
contributes (1/k)·Σ_lm [φ_a,l(x_a)φ_a,l(p_ia)/𝒩_a,l] G_lm [φ_b,m(x_b)φ_b,m(p_ib)/𝒩_b,m].
The steady form is the same split with T_ss = F·e^{−αΛs₀}/Λ.

Asserted (collected, then asserted — item 668): preconditions (18 points, labels, c = 0, balls
inside, c* labels, the exact mirror map, t* ≥ s₀); (K) c33.check_side, the 36 ``.9e`` tokens,
host C6a node JSON custody — per the review ruling 2026-11-12 10:30 (item 972),
``c33.custody_json``'s mismatches are split: a string-only one (both strings parse to the same
float, e.g. its ``.6f`` dz against the log's signed ``+.6f``) is printed, and only value
mismatches are asserted absent; (K-s) P83's s strings; (Lim); (Bq)(i) N0 → 2N0 and (ii) the split
identity s₀ = 100 → 50 s; controls (N-eq-m) mirrored column and (N-eq-u) flat column, both
SP = 1 at ``B1_BAND``.  Printed (🟡 record, no band): weights, θ(c*, t*) per deposited watt
with its self / neighbour and short / modal splits, SP, the lever beside C39 / P83 / C33, the
lattice-local reading, (c-lit), the verdict line.  C39's three lines are read per the same
ruling as one distinct matched string each (C39 ran at -n 2 and prints its table once per rank),
the match count ≥ 1 printed.

Not claimed: that the clinical SAR is 18 balls, a fitted source, a field between lattice points;
any entry-59 / 60 figure beyond C39's printed ratios; a θ lever band, a dose, a 10 g-averaged θ,
a 3 T gel claim; that the verdict line locates the excess.
"""

from __future__ import annotations

import math
import re
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.special import erf, erfc, spherical_jn

from fem_em_solver.utils.analytical import robin_slab_eigenvalues, robin_slab_odd_eigenvalues
from tests.log_records import read_log_text

from tests.validation import test_geo38_c33_field_strength_reader as c33
from tests.validation import test_geo38_f_clinical_box_loaded as c6
from tests.validation import test_mat21_eddy_source_probe_points_closed_form as p80
from tests.validation.test_birdcage_sar_mass_averaged import PHANTOM_RHO_KG_PER_M3
from tests.validation.test_imp4_pec_cylinder_closed_form import B1_BAND
from tests.validation.test_mat13_transient_perfusion_sphere import INTERIOR_BAND
from tests.validation.test_mat21_box_closed_form import A_BOX, B_BOX, C_BOX, HALF, T_C, H, K
from tests.validation.test_mat21_box_closed_form import RHO_C as RHO_C_0
from tests.validation.test_mat21_box_dose_closed_form import T_GRID

REPO = Path(__file__).resolve().parents[2]
C39_LOG = REPO / "docs" / "testing" / "logs" / "20261111T103024Z_GEO-38-C39-W1.log"
RE_C39 = {
    "ccw point": re.compile(r"^\[GEO-38 C46\]\s+ccw hotspot point theta per P_acc: "
                            r"(\S+) \| (\S+) \| (\S+)$"),
    "ccw ball": re.compile(r"^\[GEO-38 C46\]\s+ccw ball theta per P_acc \(derived = I32 ball x "
                           r"factor\): (\S+) \| (\S+) \| (\S+)$"),
    "cw point": re.compile(r"^\[GEO-38 C46\]\s+cw hotspot point theta per P_acc: "
                           r"(\S+) \| (\S+) \| (\S+)$"),
}
SIDES = ("1.5T", "3T")
CSTAR = {"1.5T": "x-105mm dz-50mm", "3T": "x+105mm dz+50mm"}
P83_S = ("0.05285093", "0.1952049", "3.6935")  # `…P83-W1.log:75–77`
C33_SAR_LEVER = "4.101440e+00"  # `…C33-W1-n1-r2.log:129`
NEEDED_SP = "1.094495"  # `20261111T104011Z_EX-194.log:66`
S0 = 100.0
N0 = (180, 120, 30)
S0_HALF = 50.0
N1 = (250, 160, 40)
EPS = float(np.finfo(float).eps)
ALPHA = K / RHO_C_0
QUAD_KW = dict(epsabs=0.0, epsrel=1e-13, limit=200)
# Predictions written first (brief), printed only; never a band.
PRED_W = {"1.5T": 0.08182573, "3T": 0.09035190}
PRED_W_TOL = 1e-6
PRED_TH = {"1.5T": (0.920, 0.89, 0.94), "3T": (1.016, 0.98, 1.04)}
PRED_NB = 1e-3
PRED_SP = (1.1042, 1.103, 1.106)
PRED_SELF_SP = 1.104199
PRED_LEVER = (4.0784, 4.074, 4.085)
PRED_FRAC = (1.10, 1.08, 1.13)
PRED_LOCAL = (4.1018, 4.100, 4.104)
PRED_CLIT = (0.9789, 0.977, 0.981)

_CACHE = {}
_AXES = {}
_GMAT = {}
_SHORT = {}


def _say(msg):
    print(msg, flush=True)


def _fail(gates):
    return [g for g in gates if not g[1]]


def _in(v, lo, hi):
    return lo <= v <= hi


# ------------------------------------------------------------------ the closed form
def _F(u):
    """Ball factor 3 j1(u)/u (= 3(sin u − u cos u)/u³)."""
    u = np.asarray(u, dtype=float)
    return 3.0 * spherical_jn(1, u) / u


def _axis(L, n, odd_too=True):
    """(k = γ/L, 𝒩, odd mask) for one axis: even family then (optionally) odd family."""
    key = (L, n, odd_too)
    if key not in _AXES:
        bi = H * L / K
        ge = robin_slab_eigenvalues(bi, n)
        ks, ns, od = [ge / L], [L * (ge + np.sin(ge) * np.cos(ge)) / ge], [np.zeros(n, bool)]
        if odd_too:
            go = robin_slab_odd_eigenvalues(bi, n)
            ks.append(go / L)
            ns.append(L * (go - np.sin(go) * np.cos(go)) / go)
            od.append(np.ones(n, bool))
        _AXES[key] = (np.concatenate(ks), np.concatenate(ns), np.concatenate(od))
    return _AXES[key]


def _phi(kk, odd, s):
    arg = np.outer(np.asarray(s, dtype=float), kk)
    return np.where(odd[None, :], np.sin(arg), np.cos(arg))


def _gmat(t, s0, n, steady=False):
    """G_lm = Σ_n T_lmn / 𝒩_c,n (cached)."""
    key = (None if steady else float(t), float(s0), tuple(n))
    if key not in _GMAT:
        r = float(c6._a10())
        ka, _na, _oa = _axis(A_BOX, n[0])
        kb, _nb, _ob = _axis(B_BOX, n[1])
        kc, nc, _oc = _axis(C_BOX, n[2], odd_too=False)
        lam = (ka**2)[:, None, None] + (kb**2)[None, :, None] + (kc**2)[None, None, :]
        T = _F(np.sqrt(lam) * r)
        T *= np.exp(-ALPHA * lam * s0)
        T /= lam
        if not steady:
            T *= -np.expm1(-ALPHA * lam * (t - s0))
        del lam
        _GMAT[key] = np.tensordot(T, 1.0 / nc, axes=([2], [0]))
        del T
    return _GMAT[key]


def _modal(x, pts, t, s0, n, steady=False):
    """Per-ball modal kernels g_i(x) [K per W of ball power] (vector over balls)."""
    G = _gmat(t, s0, n, steady)
    ka, na, oa = _axis(A_BOX, n[0])
    kb, nb, ob = _axis(B_BOX, n[1])
    ax = _phi(ka, oa, [x[0]])[0] / na
    bx = _phi(kb, ob, [x[1]])[0] / nb
    pa = _phi(ka, oa, pts[:, 0]) * ax[None, :]
    pb = _phi(kb, ob, pts[:, 1]) * bx[None, :]
    return np.einsum("il,lm,im->i", pa, G, pb) / K


def _pinf(d, s, r):
    if s <= 0.0:
        return 1.0 if d == 0.0 else 0.0
    sig = math.sqrt(2.0 * ALPHA * s)
    if d == 0.0:
        w = r / (sig * math.sqrt(2.0))
        return float(erf(w) - 2.0 / math.sqrt(math.pi) * w * math.exp(-w * w))
    z1, z2 = (d - r) / (sig * math.sqrt(2.0)), (d + r) / (sig * math.sqrt(2.0))
    return float(0.5 * (erfc(z1) - erfc(z2))
                 + sig / (d * math.sqrt(2.0 * math.pi)) * (math.exp(-z2 * z2) - math.exp(-z1 * z1)))


def _short(x, pts, s0):
    """Per-ball short-time kernels [K per W of ball power]: (1/(ρc V_B)) ∫_0^{s0} P∞ ds."""
    key = (tuple(x), float(s0))
    if key not in _SHORT:
        r = float(c6._a10())
        vb = 4.0 / 3.0 * math.pi * r**3
        out = np.empty(len(pts))
        for i, p in enumerate(pts):
            d = float(np.sqrt(np.sum((np.asarray(x) - p) ** 2)))
            val, _err = quad(lambda s, d=d: _pinf(d, s, r), 0.0, s0, **QUAD_KW)
            out[i] = val / (RHO_C_0 * vb)
        _SHORT[key] = out
    return _SHORT[key]


def _kern(x, pts, t, s0=S0, n=N0, steady=False):
    """(short, modal) per-ball kernels at x."""
    return _short(x, pts, s0), _modal(x, pts, t, s0, n, steady)


def _theta(side, w, t, s0=S0, n=N0, steady=False):
    d = _data()
    x = d["pt"][side]
    sh, mo = _kern(x, d["pts"], t, s0, n, steady)
    return float(np.dot(w, sh + mo))


# ------------------------------------------------------------------ data
def _data():
    if _CACHE:
        return _CACHE
    page = c33._load_page()
    hit = c33.find_log(c33.ENTRY58_GLOB)
    parsed = {"1.5T": page.parse_c6_log(page.SOURCE_LOG), "3T": page.parse_c6_log(hit[0])}
    lat, off = p80._lattice()
    labels = [lab for lab, _p, _w in lat]
    pts = np.asarray([pt if pt is not None else (np.nan,) * 3 for _l, pt, _w in lat], dtype=float)
    sar, w, s, idx, pt = {}, {}, {}, {}, {}
    for sd, p in parsed.items():
        by = {row["label"]: row["sar10g"] for row in p["lattice"]}
        sar[sd] = np.array([float(by[lab]) for lab in labels])
        w[sd] = sar[sd] / sar[sd].sum()
        nm = p["normalisers"]
        s[sd] = float(nm["p_box"]) / float(nm["p_acc"])
        idx[sd] = labels.index(p["c_star"]["label"])
        pt[sd] = tuple(float(v) for v in pts[idx[sd]])
    mirror = []
    for i in range(len(labels)):
        m = [j for j in range(len(labels))
             if pts[j, 0] == -pts[i, 0] and pts[j, 1] == -pts[i, 1] and pts[j, 2] == 0.0]
        mirror.append(m[0] if len(m) == 1 else None)
    _CACHE.update(page=page, hit=hit, parsed=parsed, lat=lat, off=off, labels=labels, pts=pts,
                  sar=sar, w=w, s=s, idx=idx, pt=pt, mirror=mirror, r=float(c6._a10()))
    return _CACHE


def _sp(w15, w3, t=None, s0=S0, n=N0, steady=False):
    t = float(T_GRID[-1]) if t is None else t
    return _theta("3T", w3, t, s0, n, steady) / _theta("1.5T", w15, t, s0, n, steady)


# ------------------------------------------------------------------ tests
def test_preconditions():
    gates = []
    d = _data()
    lab6 = [lab for lab, _X in c6._lattice()]
    gates.append((("n", len(d["labels"])), len(d["labels"]) == 18, len(d["labels"])))
    gates.append((("labels == c6._lattice()",), d["labels"] == lab6, 0))
    for sd in SIDES:
        pl = {row["label"] for row in d["parsed"][sd]["lattice"]}
        gates.append((("parsed labels", sd), pl == set(lab6) and len(d["parsed"][sd]["lattice"]) == 18,
                      0))
    r = d["r"]
    bound = 4.0 * EPS * (abs(d["off"]) + C_BOX)
    worst_c, worst_wall = 0.0, -np.inf
    for lab, p, why in d["lat"]:
        gates.append((("inside", lab, why), p is not None, why))
        if p is None:
            continue
        worst_c = max(worst_c, abs(p[2]))
        clear = max(abs(p[k]) + r - HALF[k] for k in range(3))
        worst_wall = max(worst_wall, clear)
        gates.append((("ball inside", lab, clear), clear <= 0.0, clear))
    gates.append((("c=0", worst_c), worst_c <= bound, worst_c))
    for sd in SIDES:
        lab = d["parsed"][sd]["c_star"]["label"]
        gates.append((("c* label", sd, lab), lab == CSTAR[sd], lab))
    mok = all(m is not None for m in d["mirror"])
    gates.append((("mirror map",), mok, d["mirror"]))
    mc = mok and d["mirror"][d["idx"]["1.5T"]] == d["idx"]["3T"]
    gates.append((("mirror c*",), mc, 0))
    ts = float(T_GRID[-1])
    gates.append((("t* >= S0", ts), ts >= S0, ts))
    _say(f"MAT21-49 lattice {len(d['labels'])} points, labels == c6._lattice() in order "
         f"{d['labels'] == lab6} (ASSERTED); both parses carry the same 18 labels (ASSERTED)")
    _say(f"MAT21-49 max |c| {worst_c:.3e} m (bound 4 eps (|off| + C) {bound:.3e}, ASSERTED); 10 g "
         f"ball r {r:.6f} m, worst max_d(|p_d| + r - L_d) {worst_wall:.6f} m (<= 0 ASSERTED; P80 "
         f"read -0.031635)")
    _say(f"MAT21-49 c* 1.5 T {d['parsed']['1.5T']['c_star']['label']!r} idx {d['idx']['1.5T']}, "
         f"3 T {d['parsed']['3T']['c_star']['label']!r} idx {d['idx']['3T']} (labels ASSERTED); "
         f"mirror (a, b, 0) -> (-a, -b, 0) exact on all 18 {mok}, maps c*1.5 to c*3 {mc} "
         f"(ASSERTED)")
    _say(f"MAT21-49 t* = T_GRID[-1] {ts:.6f} s >= S0 {S0:g} s (ASSERTED); T_C {T_C:.4f} s, "
         f"alpha {ALPHA:.6e} m^2/s")
    bad = _fail(gates)
    assert not bad, bad


def test_k_and_ks():
    gates = []
    d = _data()
    gates.append((("entry 58 kind", d["hit"]), d["hit"] is not None and d["hit"][1] == "tracked",
                  0))
    for sd, node in (("1.5T", c33.ENTRY52_NODE_PF), ("3T", c33.ENTRY58_NODE_PF)):
        res = c33.check_side(d["parsed"][sd], node, sd)
        for nme, ok, det in res:
            gates.append((("check_side", nme, det), ok, 0))
        _say(f"MAT21-49 (K) c33.check_side {sd} node {node}: {sum(ok for _n, ok, _ in res)} / "
             f"{len(res)} pass (all ASSERTED)")
    n_fmt = 0
    for sd in SIDES:
        for row in d["parsed"][sd]["lattice"]:
            n_fmt += int(f"{float(row['sar10g']):.9e}" == row["sar10g"])
    gates.append((("K-lat fmt", n_fmt), n_fmt == 36, n_fmt))
    _say(f"MAT21-49 (K) 36 lattice SAR10g tokens re-format .9e equal {n_fmt} / 36 (ASSERTED); "
         f"entry 52 {d['parsed']['1.5T']['log']}, entry 58 {d['hit'][0].name} kind {d['hit'][1]} "
         f"(tracked ASSERTED)")
    for sd, node in (("1.5T", c33.ENTRY52_NODE_PF), ("3T", c33.ENTRY58_NODE_PF)):
        jp = c33.HOST_JSON_DIR / f"c6-node-{node}.json"
        present, nf, bad = c33.custody_json(d["parsed"][sd], jp)
        if present:
            # Review ruling 2026-11-12 10:30 (item 972): a string-only mismatch (both strings
            # parse to the same float -- c33.custody_json's `.6f` dz vs the C6a log's signed
            # `{z:+.6f}`) is printed, not a red; a value mismatch is asserted absent.
            str_only, value_bad = [], []
            for nm_, a_, b_ in bad:
                try:
                    same = float(a_) == float(b_)
                except (TypeError, ValueError):
                    same = False
                (str_only if same else value_bad).append((nm_, a_, b_))
            for nm_, a_, b_ in str_only:
                _say(f"MAT21-49 (K) string-only: {nm_} {a_} | {b_}")
            gates.append((("custody value", sd, value_bad[:1]), not value_bad, len(value_bad)))
            _say(f"MAT21-49 (K) host C6a node JSON {jp.name} present (custody; the host JSON was "
                 f"last rewritten by entry 60, so this compares the log against entry 60's "
                 f"rewrite): {nf} fields, {len(bad)} string mismatches = {len(str_only)} "
                 f"string-only + {len(value_bad)} value (0 value ASSERTED)")
        else:
            _say(f"MAT21-49 (K) host C6a node JSON {jp} absent (custody only, never a red)")
    _say(f"MAT21-49 (K) C33's own JSON {c33.C33_JSON_NAME} is lattice-free (side() + ratios only); "
         f"not read -- the lattice comes from the tracked C6a logs")
    s15, s3 = d["s"]["1.5T"], d["s"]["3T"]
    ks = (f"{s15:.7g}", f"{s3:.7g}", f"{s3 / s15:.5g}")
    gates.append((("K-s", ks, P83_S), ks == P83_S, 0))
    _say(f"MAT21-49 (K-s) s = P_box/P_acc from the .9e tokens: 1.5 T {ks[0]}, 3 T {ks[1]}, share "
         f"{ks[2]} vs P83 {P83_S}: equal {ks == P83_S} (strings ASSERTED)")
    bad = _fail(gates)
    assert not bad, bad


def test_lim_bq():
    gates = []
    d = _data()
    ts = float(T_GRID[-1])
    tl = 1000.0 * T_C
    lim = {}
    for sd in SIDES:
        tr = _theta(sd, d["w"][sd], tl)
        ss = _theta(sd, d["w"][sd], None, steady=True)
        lim[sd] = abs(tr - ss) / abs(ss)
        gates.append((("Lim", sd, tr, ss), lim[sd] <= B1_BAND, lim[sd]))
    _say(f"MAT21-49 (Lim) transient(1000 T_C) vs steady split form at each c*: 1.5 T "
         f"{lim['1.5T']:.3e}, 3 T {lim['3T']:.3e} (band {B1_BAND:g}, ASSERTED; predicted <= 1e-15)")
    base = {sd: _theta(sd, d["w"][sd], ts) for sd in SIDES}
    base["SP"] = base["3T"] / base["1.5T"]
    two_n = tuple(2 * v for v in N0)
    alt_i = {sd: _theta(sd, d["w"][sd], ts, S0, two_n) for sd in SIDES}
    alt_i["SP"] = alt_i["3T"] / alt_i["1.5T"]
    alt_ii = {sd: _theta(sd, d["w"][sd], ts, S0_HALF, N1) for sd in SIDES}
    alt_ii["SP"] = alt_ii["3T"] / alt_ii["1.5T"]
    for tag, alt in (("(i) N0 -> 2 N0", alt_i), ("(ii) s0 100 -> 50 s at N1", alt_ii)):
        rels = {k: abs(alt[k] - base[k]) / abs(base[k]) for k in ("1.5T", "3T", "SP")}
        for k, v in rels.items():
            gates.append((("Bq", tag, k, base[k], alt[k]), v <= B1_BAND, v))
        _say(f"MAT21-49 (Bq) {tag}: rel theta1.5(c*1.5, t*) {rels['1.5T']:.3e}, theta3(c*3, t*) "
             f"{rels['3T']:.3e}, SP {rels['SP']:.3e} (band {B1_BAND:g}, ASSERTED; predicted "
             f"{'<= 1e-15' if tag.startswith('(i)') else '<= 1e-13'})")
    for sd in SIDES:
        x = d["pt"][sd]
        sh100, mo100 = _kern(x, d["pts"], ts)
        sh50, mo50 = _kern(x, d["pts"], ts, S0_HALF, N1)
        w = d["w"][sd]
        _say(f"MAT21-49 (Bq) {sd} split: s0 100 s short {np.dot(w, sh100):.9e} + modal "
             f"{np.dot(w, mo100):.9e} = {base[sd]:.9e}; s0 50 s short {np.dot(w, sh50):.9e} + modal "
             f"{np.dot(w, mo50):.9e} = {alt_ii[sd]:.9e} K/W")
    bad = _fail(gates)
    assert not bad, bad


def test_controls():
    gates = []
    d = _data()
    w15 = d["w"]["1.5T"]
    w3m = np.array([w15[d["mirror"][i]] for i in range(len(w15))])
    spm = _sp(w15, w3m)
    flat = np.full(len(w15), 1.0 / len(w15))
    spu = _sp(flat, flat)
    gates.append((("N-eq-m", spm), abs(spm - 1.0) <= B1_BAND, abs(spm - 1.0)))
    gates.append((("N-eq-u", spu), abs(spu - 1.0) <= B1_BAND, abs(spu - 1.0)))
    _say(f"MAT21-49 (N-eq-m) 3 T weights := mirrored 1.5 T weights: |SP - 1| {abs(spm - 1.0):.3e} "
         f"(band {B1_BAND:g}, ASSERTED; predicted <= 1e-15)")
    _say(f"MAT21-49 (N-eq-u) flat weights 1/18 both sides: |SP - 1| {abs(spu - 1.0):.3e} (band "
         f"{B1_BAND:g}, ASSERTED; P83's even-source case)")
    spl = _sp(w15, w15)
    pv, lo, hi = PRED_CLIT
    _say(f"MAT21-49 (c-lit) the 1.5 T column on both sides at the two c*: SP {spl:.7f} (predicted "
         f"approx {pv}, range [{lo}, {hi}] inside {_in(spl, lo, hi)}; printed only; |SP - 1| "
         f"{abs(spl - 1.0):.4f} vs INTERIOR_BAND {INTERIOR_BAND:g})")
    bad = _fail(gates)
    assert not bad, bad


def _c39():
    text = read_log_text(C39_LOG, errors="strict").splitlines()  # .log or .log.gz (item 997)
    out = {}
    for k, rx in RE_C39.items():
        hits = [m.groups() for m in (rx.match(ln) for ln in text) if m]
        distinct = set(hits)
        out[k] = (len(hits), distinct, next(iter(distinct)) if len(distinct) == 1 else None)
    return out


def test_record_printed():
    """The record (no band); asserts only that it was formed and C39's lines read one distinct
    string each (count >= 1; C39 ran at -n 2 and prints its table per rank)."""
    d = _data()
    ts = float(T_GRID[-1])
    r = d["r"]
    vb = 4.0 / 3.0 * math.pi * r**3
    _say("MAT21-49 normaliser: theta per deposited watt of the 18-ball source (w_i = SAR_i / sum "
         "SAR); P_dep, P_box, P_acc never mixed; s from the .9e tokens enters only the lever")
    delta = min(HALF[k] - abs(p[k]) - r for _l, p, _w in d["lat"] for k in range(3))
    expo = (2.0 * delta) ** 2 / (4.0 * ALPHA * S0)
    _say(f"MAT21-49 image factor at s0 {S0:g} s: exp(-(2 x {delta:.6f})^2 / (4 alpha s0)) = "
         f"exp(-{expo:.4g}) = {math.exp(-expo):.3e}")
    th, self_t = {}, {}
    for sd in SIDES:
        sar, w, i = d["sar"][sd], d["w"][sd], d["idx"][sd]
        sh, mo = _kern(d["pt"][sd], d["pts"], ts)
        tot = float(np.dot(w, sh + mo))
        slf = float(w[i] * (sh[i] + mo[i]))
        th[sd], self_t[sd] = tot, slf
        nb = tot - slf
        _say(f"MAT21-49 {sd} sum SAR10g {sar.sum():.7g} W/kg; w(c*) {w[i]:.8f} (predicted "
             f"{PRED_W[sd]:.8f}, within {PRED_W_TOL:g}: {abs(w[i] - PRED_W[sd]) <= PRED_W_TOL})")
        pv, lo, hi = PRED_TH[sd]
        _say(f"MAT21-49 {sd} theta(c*, t*) per dep. W {tot:.7g} K/W (predicted approx {pv}, range "
             f"[{lo}, {hi}] inside {_in(tot, lo, hi)}); self ball {slf:.7g}, 17-ball neighbours "
             f"{nb:.4e}, neighbours / total {nb / tot:.3e} (predicted <= {PRED_NB:g}: "
             f"{nb / tot <= PRED_NB}); short {float(np.dot(w, sh)):.7g} + modal "
             f"{float(np.dot(w, mo)):.7g}")
    sp = th["3T"] / th["1.5T"]
    self_q = d["w"]["3T"][d["idx"]["3T"]] / d["w"]["1.5T"][d["idx"]["1.5T"]]
    pv, lo, hi = PRED_SP
    _say(f"MAT21-49 SP = theta3(c*3)/theta1.5(c*1.5) per dep. W {sp:.7f} (predicted approx {pv}, "
         f"range [{lo}, {hi}] inside {_in(sp, lo, hi)}); self-only w3(c*)/w1.5(c*) {self_q:.7f} "
         f"(predicted {PRED_SELF_SP}); neighbours' share SP - self-only {sp - self_q:.3e}; inside "
         f"[1.09, 1.12] {_in(sp, 1.09, 1.12)}")
    share = d["s"]["3T"] / d["s"]["1.5T"]
    lever = share * sp
    c39 = _c39()
    pv, lo, hi = PRED_LEVER
    _say(f"MAT21-49 LEVER predicted 3 T / 1.5 T = share {share:.7g} x SP {sp:.7f} = {lever:.6f} "
         f"(predicted approx {pv}, range [{lo}, {hi}] inside {_in(lever, lo, hi)})")
    for k, (n, ds, g) in c39.items():
        _say(f"MAT21-49 LEVER beside C39 {k} (`{C39_LOG.name}`, {n} matches, {len(ds)} distinct; "
             f"distinct == 1 and count >= 1 ASSERTED): "
             f"{g[2] if g else None}" + (" (cw heated at ccw's c*)" if k == "cw point" else ""))
    _say(f"MAT21-49 LEVER beside P83 {P83_S[2]} and C33's SAR lever {C33_SAR_LEVER}")
    needed = float(NEEDED_SP) - 1.0
    frac = (sp - 1.0) / needed
    pv, lo, hi = PRED_FRAC
    _say(f"MAT21-49 SP - 1 {sp - 1.0:.6f} vs needed {needed:.6f} (EX-194): fraction {frac:.4f} "
         f"(predicted approx {pv}, range [{lo}, {hi}] inside {_in(frac, lo, hi)})")
    loc = {}
    for sd in SIDES:
        nm = d["parsed"][sd]["normalisers"]
        sh, mo = _kern(d["pt"][sd], d["pts"], ts)
        loc[sd] = float(np.sum(d["sar"][sd] / float(nm["p_acc"]) * PHANTOM_RHO_KG_PER_M3 * vb
                               * (sh + mo)))
        _say(f"MAT21-49 LOCAL {sd} sum_i (SAR_i/P_acc) rho theta_i(c*) {loc[sd]:.7g} K per acc. W "
             f"(NOT per deposited watt: the balls do not deposit the box's power)")
    lr = loc["3T"] / loc["1.5T"]
    pv, lo, hi = PRED_LOCAL
    _say(f"MAT21-49 LOCAL 3 T / 1.5 T {lr:.6f} (predicted approx {pv}, range [{lo}, {hi}] inside "
         f"{_in(lr, lo, hi)}); C33's SAR lever {C33_SAR_LEVER}")
    if abs(sp - 1.0 - needed) <= 0.25 * needed:
        verdict = "lattice resolves the excess"
    elif sp - 1.0 > 1.25 * needed:
        verdict = "the excess lives above the lattice's 100 mm spacing"
    else:
        verdict = "the excess lives below the lattice's 100 mm spacing"
    _say(f"MAT21-49 VERDICT (printed, the weekly reads it; no gate): {verdict} "
         f"(SP - 1 {sp - 1.0:.6f} vs 0.094495 +/- 25 %)")
    assert all(len(ds) == 1 and n >= 1 for n, ds, _g in c39.values()), c39
    assert np.isfinite(lever)
