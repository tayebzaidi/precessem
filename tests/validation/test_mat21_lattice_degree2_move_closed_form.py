r"""`MAT-21` step 50 — P89′'s ball-source closed form on entry 62's degree-2 lattice: the
predicted degree-1 → degree-2 move of θ(c*, t*) per deposited watt (§9 item 990, §10 P90;
brief `docs/planning/chunks/MAT-21.md` "§9 item 990" and the review ruling 2026-11-15 10:30).

Numpy / scipy; complex build only because P89 imports C33 (``mpi4py``) and C6.  Sets no band
and moves none; ``B1_BAND`` and ``INTERIOR_BAND`` are imported.

The closed form is P89′'s, imported unchanged: ``p89._theta("1.5T", w, t)`` — the split
(free-space Gaussian ball integral on [0, s₀], Robin-box odd/even mode sum with the exact ball
factor F = 3 j₁(u)/u on [s₀, t]) at the 1.5 T c* point.  Degree 1: ``p89._data()["w"]["1.5T"]``
(entry 52).  Degree 2: entry 62's 18 SAR10g strings, read through the ``ans:11`` page's
``degree2_block()`` (C37's route: ``parse_c6_log(C19_LOG)``, node ``C19_NODE_PF``) and, for
(K), the page's ``parse_c6_log(page.C19_LOG)`` called here; looked up by label in
``p89._data()["labels"]`` order, w_i = SAR_i / Σ SAR.  c* is ``x-105mm dz-50mm`` on both
orders, so both θ are evaluated at the same point; the move is θ₂/θ₁ − 1.

Asserted (collected, then asserted — item 668): preconditions (18 distinct labels equal to P89's
as a set, node 185.7629, c* label on both sides, t* ≥ S0); (K) labels, sar10g and c* label /
SAR by string against the C37 block; (Lim) transient at 1000·T_C vs steady; (Bq)(i) N0 → 2 N0
and (ii) s₀ 100 → 50 s at N1, on θ₂ and on the move; the negative control (factor 1, the
weekly's, asserted on the raw SAR10g strings per the 10:30 ruling): |SAR₂(c*)/SAR₁(c*) − 1| >
INTERIOR_BAND.  Printed (🟡 record): Σ SAR10g, w(c*), θ per side with self / neighbour split,
the move beside entry 69's 10 % clause, C37's per-P_acc move string (context only), the
normalised-weight move (unasserted).

Not claimed: entry 69's FEM move or its probes P0 … P5; that the clinical SAR is 18 balls; any
per-accepted-watt figure; a 3 T degree-2 statement.
"""

from __future__ import annotations

import numpy as np

from tests.validation import test_mat21_lattice_ball_superposition_closed_form as p89
from tests.validation.test_imp4_pec_cylinder_closed_form import B1_BAND
from tests.validation.test_mat13_transient_perfusion_sphere import INTERIOR_BAND
from tests.validation.test_mat21_box_closed_form import T_C
from tests.validation.test_mat21_box_dose_closed_form import T_GRID

SIDE = "1.5T"
CSTAR = "x-105mm dz-50mm"
NODE = "185.7629"
CLAUSE = 0.10  # entry 69's decision rule, `docs/testing/xl-pending.md:2601` (printed only)
# Predictions (brief), printed only; never a band.
PRED_SUM2 = 1.639563e-03
PRED_W2 = 0.08212568
PRED_W_TOL = 1e-6
PRED_TH2 = (0.92418, 0.9235, 0.9250)
PRED_MOVE = (3.67e-03, 3.0e-03, 4.4e-03)
PRED_CTRL = -2.349e-02

_C = {}


def _say(msg):
    print(msg, flush=True)


def _fail(gates):
    return [g for g in gates if not g[1]]


def _d2():
    if _C:
        return _C
    d = p89._data()
    page = d["page"]
    block = page.degree2_block()
    own = page.parse_c6_log(page.C19_LOG)
    by = {}
    dup = []
    for row in block["parsed"]["lattice"]:
        if row["label"] in by:
            dup.append(row["label"])
        by[row["label"]] = row["sar10g"]
    sar2 = (
        np.array([float(by[lab]) for lab in d["labels"]])
        if not dup and set(by) == set(d["labels"])
        else None
    )
    _C.update(
        page=page,
        block=block,
        own=own,
        by=by,
        dup=dup,
        sar2=sar2,
        w2=None if sar2 is None else sar2 / sar2.sum(),
    )
    return _C


def test_preconditions_and_k():
    gates = []
    d = p89._data()
    e = _d2()
    blk, own = e["block"], e["own"]
    p2 = blk["parsed"]
    labs2 = [row["label"] for row in p2["lattice"]]
    gates.append((("dup labels", e["dup"]), not e["dup"], e["dup"]))
    gates.append(
        (("18 labels == p89 set",), len(labs2) == 18 and set(labs2) == set(d["labels"]), len(labs2))
    )
    gates.append(
        (("node", p2["node"]["c_pf"]), p2["node"]["c_pf"] == NODE == e["page"].C19_NODE_PF, 0)
    )
    gates.append((("c* deg2", p2["c_star"]["label"]), p2["c_star"]["label"] == CSTAR, 0))
    gates.append(
        (
            ("c* deg1", d["parsed"][SIDE]["c_star"]["label"]),
            d["parsed"][SIDE]["c_star"]["label"] == CSTAR == p89.CSTAR[SIDE],
            0,
        )
    )
    ts = float(T_GRID[-1])
    gates.append((("t* >= S0", ts), ts >= p89.S0, ts))
    # (K) by string, the module's own parse against the C37 block's.
    k_lab = sum(a["label"] == b["label"] for a, b in zip(own["lattice"], p2["lattice"]))
    k_sar = sum(a["sar10g"] == b["sar10g"] for a, b in zip(own["lattice"], p2["lattice"]))
    k_cs = int(own["c_star"]["label"] == p2["c_star"]["label"]) + int(
        own["c_star"]["sar10g"] == p2["c_star"]["sar10g"]
    )
    n_ok = len(own["lattice"]) == len(p2["lattice"]) == 18
    gates.append((("K labels", k_lab), n_ok and k_lab == 18, k_lab))
    gates.append((("K sar10g", k_sar), n_ok and k_sar == 18, k_sar))
    gates.append((("K c*", k_cs), k_cs == 2, k_cs))
    cs_in_lat = e["by"].get(p2["c_star"]["label"]) == p2["c_star"]["sar10g"]
    _say(
        f"MAT21-50 entry 62 {blk['source_log']}: {len(labs2)} lattice labels, distinct "
        f"{not e['dup']}, == P89's label set {set(labs2) == set(d['labels'])} (ASSERTED); node "
        f"{p2['node']['c_pf']} pF (== C19_NODE_PF ASSERTED)"
    )
    _say(
        f"MAT21-50 c* degree 1 {d['parsed'][SIDE]['c_star']['label']!r}, degree 2 "
        f"{p2['c_star']['label']!r} (both == {CSTAR!r} ASSERTED); C37 c_star_move "
        f"{blk['c_star_move']!r} (printed); c* SAR string in lattice row equal {cs_in_lat}"
    )
    _say(
        f"MAT21-50 (K) own parse_c6_log(C19_LOG) vs degree2_block()['parsed']: labels "
        f"{k_lab} / 18, sar10g {k_sar} / 18, c* label + SAR {k_cs} / 2 (strings ASSERTED)"
    )
    _say(f"MAT21-50 t* = T_GRID[-1] {ts:.6f} s >= S0 {p89.S0:g} s (ASSERTED)")
    bad = _fail(gates)
    assert not bad, bad


def test_lim_bq():
    gates = []
    d = p89._data()
    e = _d2()
    assert e["w2"] is not None
    w1, w2 = d["w"][SIDE], e["w2"]
    ts = float(T_GRID[-1])
    tl = 1000.0 * T_C
    tr = p89._theta(SIDE, w2, tl)
    ss = p89._theta(SIDE, w2, None, steady=True)
    lim = abs(tr - ss) / abs(ss)
    gates.append((("Lim", tr, ss), lim <= B1_BAND, lim))
    _say(
        f"MAT21-50 (Lim) degree-2 theta(c*) transient(1000 T_C) vs steady: {lim:.3e} (band "
        f"{B1_BAND:g}, ASSERTED; P89 read 0.000e+00 on degree 1)"
    )
    two_n = tuple(2 * v for v in p89.N0)
    base = {"th2": p89._theta(SIDE, w2, ts), "th1": p89._theta(SIDE, w1, ts)}
    base["move"] = base["th2"] / base["th1"] - 1.0
    for tag, (s0, n) in (
        ("(i) N0 -> 2 N0", (p89.S0, two_n)),
        ("(ii) s0 100 -> 50 s at N1", (p89.S0_HALF, p89.N1)),
    ):
        a2 = p89._theta(SIDE, w2, ts, s0, n)
        a1 = p89._theta(SIDE, w1, ts, s0, n)
        r_th = abs(a2 - base["th2"]) / abs(base["th2"])
        mv = a2 / a1 - 1.0
        r_mv = abs(mv - base["move"]) / abs(base["move"])
        gates.append((("Bq th2", tag, a2), r_th <= B1_BAND, r_th))
        gates.append((("Bq move", tag, mv), r_mv <= B1_BAND, r_mv))
        _say(
            f"MAT21-50 (Bq) {tag}: rel theta2(c*, t*) {r_th:.3e}, rel move {r_mv:.3e} (band "
            f"{B1_BAND:g}, ASSERTED; P89 read 2.170e-15 / 3.738e-15 on theta1)"
        )
    bad = _fail(gates)
    assert not bad, bad


def test_control_raw_strings():
    """Negative control (factor 1, asserted, the weekly's; 10:30 ruling: on the raw strings)."""
    d = p89._data()
    e = _d2()
    s1 = d["parsed"][SIDE]["c_star"]["sar10g"]
    s2 = e["block"]["parsed"]["c_star"]["sar10g"]
    rel = float(s2) / float(s1) - 1.0
    w1c = d["w"][SIDE][d["idx"][SIDE]]
    w2c = e["w2"][d["labels"].index(CSTAR)]
    wmv = w2c / w1c - 1.0
    _say(
        f"MAT21-50 (N-ctrl) raw SAR10g(c*) strings degree 1 {s1} -> degree 2 {s2}: rel "
        f"{rel:+.4e} (|.| > INTERIOR_BAND {INTERIOR_BAND:g} ASSERTED; predicted {PRED_CTRL:+.3e})"
    )
    _say(
        f"MAT21-50 (N-ctrl, printed, unasserted) normalised weight w(c*) move {wmv:+.4e} "
        f"(|.| > INTERIOR_BAND {abs(wmv) > INTERIOR_BAND}; the control's per-watt re-posing is "
        f"the weekly's)"
    )
    assert abs(rel) > INTERIOR_BAND, (s1, s2, rel)


def test_record_printed():
    d = p89._data()
    e = _d2()
    ts = float(T_GRID[-1])
    i = d["labels"].index(CSTAR)
    assert i == d["idx"][SIDE]
    _say(
        "MAT21-50 normaliser: theta per deposited watt of the 18-ball source (w_i = SAR_i / sum "
        "SAR, each order per its own ccw drive); P_acc and P_box never mixed in"
    )
    sh, mo = p89._kern(d["pt"][SIDE], d["pts"], ts)
    th = {}
    for tag, sar, w in (
        ("degree 1", d["sar"][SIDE], d["w"][SIDE]),
        ("degree 2", e["sar2"], e["w2"]),
    ):
        tot = float(np.dot(w, sh + mo))
        slf = float(w[i] * (sh[i] + mo[i]))
        nb = tot - slf
        th[tag] = tot
        _say(
            f"MAT21-50 {tag} sum SAR10g {sar.sum():.7g} W/kg; w(c*) {w[i]:.8f}; theta(c*, t*) "
            f"per dep. W {tot:.7g} K/W; self {slf:.7g}, neighbours {nb:.4e}, neighbours / total "
            f"{nb / tot:.3e}"
        )
    _say(
        f"MAT21-50 predictions: degree-2 sum {PRED_SUM2:.6e}; w2(c*) {PRED_W2:.8f} within "
        f"{PRED_W_TOL:g}: {abs(e['w2'][i] - PRED_W2) <= PRED_W_TOL}; theta2 approx {PRED_TH2[0]} "
        f"in [{PRED_TH2[1]}, {PRED_TH2[2]}]: {PRED_TH2[1] <= th['degree 2'] <= PRED_TH2[2]}"
    )
    move = th["degree 2"] / th["degree 1"] - 1.0
    pv, lo, hi = PRED_MOVE
    _say(
        f"MAT21-50 MOVE theta2/theta1 - 1 at c* {move:+.3e} (predicted approx {pv:+.2e}, range "
        f"[{lo:+.1e}, {hi:+.1e}] inside {lo <= move <= hi}); beside entry 69's clause {CLAUSE:.2f}"
    )
    fires = abs(move) > CLAUSE
    _say(
        f"MAT21-50 closed form predicts the 10 % clause {'does' if fires else 'does not'} fire "
        f"at c* (printed, never asserted; P0 ... P5 not posed)"
    )
    c37 = [m for m in e["block"]["moves"] if m["row"] == "SAR10g(c*) per P_acc"]
    _say(
        f"MAT21-50 context only: C37's 'SAR10g(c*) per P_acc' move "
        f"{c37[0]['move'] if c37 else None} (per ACCEPTED watt, not per deposited watt)"
    )
    assert np.isfinite(move) and len(c37) == 1
