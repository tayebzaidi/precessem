"""`WF-10` H27 (§9 item 491, §10 chain H27): the 18-point 10 g lattice at degree 2,
read from `xl` entry 35's window, its mirror gate re-derived, and the per-point
d1 -> d2 move recorded.

Post-processing only (brief re-tier, review ruling 2026-10-15 18:00): numpy and
regex over two harness logs, no mesh, no solve, never ``build_tuned_quadrature``.

* (K1) coding control on known data — the parser + recomputation on item 359's
  degree-1 table (``20261010T171705Z_WF-10-step2s-W1.log``) reproduces its
  printed (iii-a) / (iii-b) / per-row rel / z_c spreads to <= 6e-7 absolute
  (the window prints ``.4f`` %, i.e. +-5e-7).  ASSERTED.
* (R1) the same reproduction on entry 35's degree-2 table.  ASSERTED.
* Label identity — the 18 degree-2 labels equal the degree-1 labels in order,
  and every ccw ball centre (c* included) equals ``_lattice(az17, z_c)`` for
  its label to 1e-3 mm (the printed ``.3f``).  ASSERTED.
* (M2) recomputed (iii-a), (iii-b) <= ``C4_COVARIANCE_BAND``.  ASSERTED.
* No-flip negative control — ccw at a z_c point against ccw at its
  feed-bisector image (rows az17+0 <-> +90, +180 <-> +270), relative L2 over
  the 12 z_c points, > ``C4_COVARIANCE_BAND``.  ASSERTED (predicted [10, 30] %;
  17.93 % at degree 1).

Printed only (rule (e), a record, not a gate): the per-point d1 -> d2 move
(ccw and cw(M), ``d2/d1 - 1``, max |move| and RMS); the same move per watt
reaching the phantom (review ruling 2026-10-16 03:00); the d2 spread beside
item 359's, H18's and H24's; c* and SAR10g(c*) d1 -> d2.

Scope: degree 2 on the 0.015 m phantom at the degree-1 root C_in = 118.2437 pF,
two-port drive.  No absolute-SAR, tuning or compliance claim; the mirror gate is
a symmetry identity, the moves are an order record (two orders, one mesh), not a
convergence rate.
"""

from __future__ import annotations

import os
import re

import numpy as np
import pytest

from tests.log_records import glob_log_paths, read_log_text

H27_ENV = "FEM_EM_WF10_H27"
LOG_DIR = "docs/testing/logs"
# Item 359's degree-1 table, 0.015 m phantom — the only table differenced (trap 1).
D1_TABLE_LOG = "20261010T171705Z_WF-10-step2s-W1.log"
# Degree-1 ccw P_phantom for the per-watt move (03:00 ruling): step-4b re-gate log,
# the ``[WF-10 step2] [ccw] w = ... P_phantom`` line (``:6325``).
D1_PPH_LOG = "20261013T111006Z_WF-10-step4b-W1-regate.log"
# Entry 35's window (exactly one, condition 2).
D2_LOG_GLOB = "20261016T*_WF-10-step4-H20.log"
# H18's spread (d1, 0.0075 m phantom) — printed beside, never differenced (trap 1).
H18_LOG = "20261011T125304Z_WF-10-step3b-W3.log"
# H24's four-port spread, r = 15 / 30 / 45 mm (`20261013T141801Z_WF-10-step5-W1.log`,
# item 422 outcome) — printed beside only.
H24_FOUR_PORT_SPREAD_PCT = (3.2122, 4.2668, 0.7809)

PRINT_ROUNDING_ABS = 6.0e-7  # `.4f` % prints are +-5e-7 absolute
COORD_TOL_MM = 1.0e-3  # `.3f` mm prints
EXPECTED_POINTS = 18
REDUCED_POINTS = 12
# (K1) known values at degree 1, from the brief (drafter's stdlib re-computation).
K1_IIIB = 0.0127886
K1_SPREADS = (0.167267, 0.209986, 0.186041)
K1_NOFLIP_PREDICTED = 0.1793
PREDICTED_NOFLIP = (0.10, 0.30)
PREDICTED_MIRROR = (0.003, 0.03)

_NUM = r"([-+]?\d+\.\d+(?:e[-+]?\d+)?)"
ROW_RE = re.compile(
    r"^\[WF-10 step2s\]   (r=\d+ dz=[-+]\d+ az17[-+]\d+)\s+ccw "
    + _NUM
    + r" cw\(M\) "
    + _NUM
    + r" W/kg, rel "
    + _NUM
    + r"%(\s+<- c\*)?\s*$"
)
HEAD_RE = re.compile(r"^\[WF-10 step2s\] ccw 10 g table \((\d+) points(, REDUCED)?\)")
CSTAR_RE = re.compile(
    r"^\[WF-10 step2s\] c\* = (r=\d+ dz=[-+]\d+ az17[-+]\d+) \("
    + _NUM
    + r", "
    + _NUM
    + r", "
    + _NUM
    + r"\) mm, SAR10g\(c\*\) "
    + _NUM
)
IIIA_RE = re.compile(r"^\[WF-10 step2s\] \(iii-a\) mirror at c\*: " + _NUM + r"%")
IIIB_RE = re.compile(
    r"^\[WF-10 step2s\] \(iii-b\) mirror lattice rel L2 over (\d+) points: " + _NUM + r"%"
)
SPREAD_RE = re.compile(
    r"^\[WF-10 step2s\] azimuthal spread of the four z_c rows .*?: "
    r"r=15 mm " + _NUM + r"%, r=30 mm " + _NUM + r"%, r=45 mm " + _NUM + r"%"
)
H18_SPREAD_RE = re.compile(
    r"^\[WF-10 step3b\] azimuthal spread of the four z_c rows .*?: "
    r"r=15 mm " + _NUM + r"%.*?r=30 mm " + _NUM + r"%.*?r=45 mm " + _NUM + r"%"
)
GEOM_RE = re.compile(
    r"^\[WF-10 step2\] phantom bbox .*centre \("
    + _NUM
    + r", "
    + _NUM
    + r", "
    + _NUM
    + r"\) m,.*az17 "
    + _NUM
    + r" deg, mirror plane "
    + _NUM
    + r" deg"
)
BALL_RE = re.compile(
    r"^\[WF-10 step2\]   ball ccw (r=\d+ dz=[-+]\d+ az17[-+]\d+)\s+at \(\s*"
    + _NUM
    + r",\s*"
    + _NUM
    + r",\s*"
    + _NUM
    + r"\) mm: SAR10g "
    + _NUM
)
PPH_CCW_RE = re.compile(r"^\[WF-10 step2\] \[ccw\] w = .*P_phantom " + _NUM + r" ")
SIGMA_P_RE = re.compile(r"^\[WF-10 step2s\] \(Sigma_P\) P_phantom " + _NUM)


def _path(name):
    return os.path.join(LOG_DIR, name)


def _d2_log():
    found = glob_log_paths(_path(D2_LOG_GLOB))
    assert len(found) == 1, f"condition 2: expected exactly one {D2_LOG_GLOB}, found {found}"
    return found[0]


def parse_table(path):
    """Parse one 10 g lattice window: the rank-0 table, c*, (iii-a/b), spreads, the
    ccw ball centres, phantom geometry and the sar10g module's ccw P_phantom.

    Only the ``[WF-10 step2s]   r=`` rows are read for the table (trap 2); the
    count must match the header (18, or 12 if REDUCED).
    """
    lines = read_log_text(path).splitlines()
    head_idx = [i for i, ln in enumerate(lines) if HEAD_RE.match(ln)]
    assert len(head_idx) == 1, f"{path}: {len(head_idx)} table headers"
    h = head_idx[0]
    m = HEAD_RE.match(lines[h])
    n_head, reduced = int(m.group(1)), bool(m.group(2))
    rows = []
    for ln in lines[h + 1 :]:
        mm = ROW_RE.match(ln)
        if not mm:
            break
        rows.append(
            (
                mm.group(1),
                float(mm.group(2)),
                float(mm.group(3)),
                float(mm.group(4)) / 100.0,
                bool(mm.group(5)),
            )
        )
    assert len(rows) == n_head, (path, len(rows), n_head)
    assert n_head == (REDUCED_POINTS if reduced else EXPECTED_POINTS), (n_head, reduced)

    def one(rx, start=h):
        hits = [rx.match(ln) for ln in lines[start:]]
        hits = [x for x in hits if x]
        assert hits, f"{path}: no line matching {rx.pattern[:50]}"
        return hits[0]

    cs = one(CSTAR_RE)
    iiia = float(one(IIIA_RE).group(1)) / 100.0
    ib = one(IIIB_RE)
    sp = one(SPREAD_RE)
    sigma_p = float(one(SIGMA_P_RE).group(1))
    # Geometry / ccw balls / ccw P_phantom: the sar10g module's block, i.e. the last
    # occurrence BEFORE the table (entry 35's log holds the b1_sar module first).
    geo = [GEOM_RE.match(ln) for ln in lines[:h]]
    geo = [g for g in geo if g][-1]
    pph = [PPH_CCW_RE.match(ln) for ln in lines[:h]]
    pph = [p for p in pph if p][-1]
    balls = {}
    for ln in lines[:h]:
        b = BALL_RE.match(ln)
        if b:
            balls[b.group(1)] = (
                np.array([float(b.group(i)) for i in (2, 3, 4)]),
                float(b.group(5)),
            )
    return {
        "path": path,
        "reduced": reduced,
        "labels": [r[0] for r in rows],
        "ccw": np.array([r[1] for r in rows]),
        "cw": np.array([r[2] for r in rows]),
        "rel_printed": np.array([r[3] for r in rows]),
        "star_flags": [r[4] for r in rows],
        "cstar_label": cs.group(1),
        "cstar_mm": np.array([float(cs.group(i)) for i in (2, 3, 4)]),
        "cstar_sar": float(cs.group(5)),
        "iiia_printed": iiia,
        "iiib_n": int(ib.group(1)),
        "iiib_printed": float(ib.group(2)) / 100.0,
        "spreads_printed": np.array([float(sp.group(i)) for i in (1, 2, 3)]) / 100.0,
        "sigma_p": sigma_p,
        "p_phantom_ccw": float(pph.group(1)),
        "z_c": float(geo.group(3)),
        "az17": float(geo.group(4)),
        "az_m": float(geo.group(5)),
        "balls": balls,
    }


def recompute(t):
    """(iii-a), (iii-b), per-row rel, z_c spreads and the no-flip control from the
    parsed numbers, with the window's own formulas (``_relative_l2``, ``:365–373``)."""
    from tests.validation.test_birdcage_b1_plus_map import _relative_l2
    from tests.validation.test_wf7_f_human_sar10g_hotspot import LATTICE_R_M

    labels, ccw, cw = t["labels"], t["ccw"], t["cw"]
    k_star = int(np.argmax(ccw))
    mask = np.ones(len(labels), dtype=bool)
    out = {
        "cstar_label": labels[k_star],
        "iiia": abs(cw[k_star] - ccw[k_star]) / ccw[k_star],
        "iiib": _relative_l2(cw, ccw, mask),
        "rel": np.abs(cw / ccw - 1.0),
    }
    spreads = []
    for r in LATTICE_R_M:
        tag = f"r={r * 1e3:.0f} dz=+0 "
        vals = [ccw[i] for i, lab in enumerate(labels) if lab.startswith(tag)]
        assert len(vals) == 4, (tag, vals)
        spreads.append((max(vals) - min(vals)) / np.mean(vals))
    out["spreads"] = np.array(spreads)
    # z_c rows: index d*3 + i over d in (0, 90, 180, 270); image rows 0<->1, 2<->3.
    rows = ccw[:12].reshape(4, 3)
    img = rows[[1, 0, 3, 2]]
    out["noflip"] = float(np.linalg.norm(img - rows) / np.linalg.norm(rows))
    out["flip12"] = float(np.linalg.norm(cw[:12] - ccw[:12]) / np.linalg.norm(ccw[:12]))
    return out


def _check_reproduction(t, rc, tag):
    """(K1) / (R1): recomputed vs printed, |delta| <= PRINT_ROUNDING_ABS."""
    d_a = abs(rc["iiia"] - t["iiia_printed"])
    d_b = abs(rc["iiib"] - t["iiib_printed"])
    d_rel = float(np.max(np.abs(rc["rel"] - t["rel_printed"])))
    d_sp = float(np.max(np.abs(rc["spreads"] - t["spreads_printed"])))
    print(
        f"[WF-10 H27] ({tag}) {os.path.basename(t['path'])}: {len(t['labels'])} rows; "
        f"c* {rc['cstar_label']} (logged {t['cstar_label']}); (iii-a) {rc['iiia']:.7e} vs "
        f"printed {t['iiia_printed']:.6e} |d| {d_a:.2e}; (iii-b) {rc['iiib']:.7e} vs "
        f"{t['iiib_printed']:.6e} |d| {d_b:.2e}; per-row rel max |d| {d_rel:.2e}; spreads "
        f"{np.array2string(rc['spreads'], precision=7)} max |d| {d_sp:.2e} (all <= "
        f"{PRINT_ROUNDING_ABS:g} ASSERTED)",
        flush=True,
    )
    assert rc["cstar_label"] == t["cstar_label"], (rc["cstar_label"], t["cstar_label"])
    assert [i for i, f in enumerate(t["star_flags"]) if f] == [t["labels"].index(t["cstar_label"])]
    assert abs(t["cstar_sar"] - t["ccw"][t["labels"].index(t["cstar_label"])]) == 0.0
    assert t["iiib_n"] == len(t["labels"])
    assert d_a <= PRINT_ROUNDING_ABS, f"({tag}) (iii-a) {d_a:.3e}"
    assert d_b <= PRINT_ROUNDING_ABS, f"({tag}) (iii-b) {d_b:.3e}"
    assert d_rel <= PRINT_ROUNDING_ABS, f"({tag}) per-row rel {d_rel:.3e}"
    assert d_sp <= PRINT_ROUNDING_ABS, f"({tag}) spreads {d_sp:.3e}"


def _check_label_geometry(t, tag):
    """Every ccw ball centre (c* included) is ``_lattice(az17, z_c)`` for its label."""
    from tests.validation.test_wf10_f_human_tuned_sar10g import _lattice

    lat = {lab: np.array(c) * 1e3 for lab, *_x, c in _lattice(t["az17"], t["z_c"])}
    worst = 0.0
    for lab in t["labels"]:
        assert lab in t["balls"], f"({tag}) no ccw ball line for {lab}"
        worst = max(worst, float(np.max(np.abs(t["balls"][lab][0] - lat[lab]))))
        assert t["balls"][lab][1] == t["ccw"][t["labels"].index(lab)], lab
    d_star = float(np.max(np.abs(t["cstar_mm"] - lat[t["cstar_label"]])))
    print(
        f"[WF-10 H27] ({tag}) label geometry: az17 {t['az17']:.9f} deg, z_c {t['z_c']:.6e} m; "
        f"worst ccw ball centre vs _lattice {worst:.2e} mm, c* {d_star:.2e} mm (<= "
        f"{COORD_TOL_MM:g} mm ASSERTED)",
        flush=True,
    )
    assert worst <= COORD_TOL_MM, worst
    assert d_star <= COORD_TOL_MM, d_star


def test_wf10_h27_parser_coding():
    """(L0)-style: the row regex reads a synthetic table, rejects a truncated one,
    and the no-flip permutation pairs rows 0<->1, 2<->3."""
    good = (
        "[WF-10 step2s]   r=45 dz=+0 az17+0          ccw 2.873208988e-03 cw(M) "
        "2.844402845e-03 W/kg, rel 1.0026%   <- c*"
    )
    m = ROW_RE.match(good)
    assert m and m.group(1) == "r=45 dz=+0 az17+0" and m.group(5)
    assert float(m.group(2)) == 2.873208988e-03 and float(m.group(4)) == 1.0026
    plain = good.replace("   <- c*", "")
    assert ROW_RE.match(plain) and not ROW_RE.match(plain).group(5)
    assert not ROW_RE.match("[WF-10 step2]   ball ccw r=15 dz=+0 az17+0  at (...)")
    # Permutation: a C4-perfect row set with rows 0 == 1 and 2 == 3 gives no-flip 0.
    ccw = np.array([1.0, 2.0, 3.0] * 2 + [4.0, 5.0, 6.0] * 2 + [7.0] * 6)
    rc_rows = ccw[:12].reshape(4, 3)
    assert np.array_equal(rc_rows[[1, 0, 3, 2]], rc_rows)
    print(
        "[WF-10 H27 L0] row regex, c* flag, ball-line rejection, image permutation "
        "(all ASSERTED)",
        flush=True,
    )


def test_wf10_h27_k1_degree1_control():
    """(K1): the parser + recomputation on item 359's degree-1 table."""
    t = parse_table(_path(D1_TABLE_LOG))
    assert not t["reduced"] and len(t["labels"]) == EXPECTED_POINTS
    rc = recompute(t)
    _check_reproduction(t, rc, "K1 d1")
    _check_label_geometry(t, "K1 d1")
    d_known = max(
        abs(rc["iiib"] - K1_IIIB), float(np.max(np.abs(rc["spreads"] - np.array(K1_SPREADS))))
    )
    print(
        f"[WF-10 H27] (K1 d1) vs the brief's known values (iii-b) {K1_IIIB} / spreads "
        f"{K1_SPREADS}: max |d| {d_known:.2e} (<= {PRINT_ROUNDING_ABS:g} ASSERTED); no-flip "
        f"{rc['noflip'] * 100:.4f}% (brief 17.93%), flip over 12 z_c {rc['flip12'] * 100:.4f}% "
        f"(brief 1.47%)",
        flush=True,
    )
    assert d_known <= PRINT_ROUNDING_ABS, d_known
    assert abs(rc["noflip"] - K1_NOFLIP_PREDICTED) <= 5e-5, rc["noflip"]


def test_wf10_h27_degree2_lattice():
    if os.environ.get(H27_ENV, "").strip() != "1":
        pytest.skip(f"{H27_ENV} unset: W1 (entry 35's degree-2 log), environment-selected")
    from tests.validation.test_birdcage_b1_plus_map import C4_COVARIANCE_BAND

    d2_path = _d2_log()
    tail = read_log_text(d2_path).splitlines()[-30:]
    assert "- Status: 0" in tail and "[capture] rc=0" in tail, "condition 3: entry 35 not green"

    d1 = parse_table(_path(D1_TABLE_LOG))
    d2 = parse_table(d2_path)
    rc1, rc2 = recompute(d1), recompute(d2)

    # ---- (R1) + label identity ------------------------------------------------------
    _check_reproduction(d2, rc2, "R1 d2")
    _check_label_geometry(d2, "R1 d2")
    n = len(d2["labels"])
    assert d1["labels"][:n] == d2["labels"], (d1["labels"], d2["labels"])
    assert abs(d2["az17"] - d1["az17"]) == 0.0 and abs(d2["z_c"] - d1["z_c"]) == 0.0
    if d2["reduced"]:
        print(
            f"[WF-10 H27] REDUCED ({n} of {EXPECTED_POINTS}): dz = +-50 mm move unrecorded",
            flush=True,
        )

    # ---- (M2) the gate + no-flip control --------------------------------------------
    print(
        f"[WF-10 H27] (M2) d2 (iii-a) mirror at c*: {rc2['iiia'] * 100:.4f}% ; (iii-b) rel L2 "
        f"over {n}: {rc2['iiib'] * 100:.4f}% (band {C4_COVARIANCE_BAND * 100:.1f}%, ASSERTED; "
        f"PREDICTED {PREDICTED_MIRROR[0] * 100:g}-{PREDICTED_MIRROR[1] * 100:g}%); d1 "
        f"{rc1['iiia'] * 100:.4f}% / {rc1['iiib'] * 100:.4f}%",
        flush=True,
    )
    print(
        f"[WF-10 H27] (NC) no-flip rel L2 over 12 z_c: d2 {rc2['noflip'] * 100:.4f}% vs flip "
        f"{rc2['flip12'] * 100:.4f}% (> band {C4_COVARIANCE_BAND * 100:.1f}% ASSERTED; "
        f"PREDICTED {PREDICTED_NOFLIP[0] * 100:g}-{PREDICTED_NOFLIP[1] * 100:g}%); d1 "
        f"{rc1['noflip'] * 100:.4f}% vs {rc1['flip12'] * 100:.4f}%; separation "
        f"{rc2['noflip'] / C4_COVARIANCE_BAND:.2f}x band",
        flush=True,
    )

    # ---- readout (printed, rule (e)) -------------------------------------------------
    mv_ccw = d2["ccw"] / d1["ccw"][:n] - 1.0
    mv_cw = d2["cw"] / d1["cw"][:n] - 1.0
    p1, p2 = d1["p_phantom_ccw"], d2["p_phantom_ccw"]
    # 03:00 ruling: the d1 ccw P_phantom from the step-4b re-gate log.
    rg = [PPH_CCW_RE.match(ln) for ln in read_log_text(_path(D1_PPH_LOG)).splitlines()]
    rg = [x for x in rg if x]
    p1_regate = float(rg[-1].group(1)) if rg else None
    pw_ratio = (p1_regate if p1_regate is not None else np.nan) / p2
    mv_ccw_w = (1.0 + mv_ccw) * pw_ratio - 1.0
    mv_cw_w = (1.0 + mv_cw) * pw_ratio - 1.0
    print(
        f"[WF-10 H27] per-point d1 -> d2 move (d2/d1 - 1), raw and per P_phantom (d1 ccw "
        f"P_phantom {p1_regate!r} W, {D1_PPH_LOG}; d2 {p2:.9e} W; 359's table log reads "
        f"{p1:.9e} W):",
        flush=True,
    )
    for i, lab in enumerate(d2["labels"]):
        print(
            f"[WF-10 H27]   {lab:<22s} ccw {d1['ccw'][i]:.6e} -> {d2['ccw'][i]:.6e} "
            f"{mv_ccw[i] * 100:+8.3f}% | cw(M) {d1['cw'][i]:.6e} -> {d2['cw'][i]:.6e} "
            f"{mv_cw[i] * 100:+8.3f}% | per-W ccw {mv_ccw_w[i] * 100:+7.3f}% cw(M) "
            f"{mv_cw_w[i] * 100:+7.3f}%",
            flush=True,
        )
    allmv = np.concatenate([mv_ccw, mv_cw])
    allmv_w = np.concatenate([mv_ccw_w, mv_cw_w])
    print(
        f"[WF-10 H27] move raw: max |move| {np.max(np.abs(allmv)) * 100:.3f}%, RMS "
        f"{np.sqrt(np.mean(allmv ** 2)) * 100:.3f}%, mean {np.mean(allmv) * 100:+.3f}% over "
        f"{allmv.size} values (PREDICTED 0.5-8%)",
        flush=True,
    )
    print(
        f"[WF-10 H27] move per P_phantom: max |move| {np.max(np.abs(allmv_w)) * 100:.3f}%, RMS "
        f"{np.sqrt(np.mean(allmv_w ** 2)) * 100:.3f}%, mean {np.mean(allmv_w) * 100:+.3f}% "
        f"(P_phantom d1/d2 {pw_ratio:.6f}; printed, never asserted)",
        flush=True,
    )
    h18 = [H18_SPREAD_RE.match(ln) for ln in read_log_text(_path(H18_LOG)).splitlines()]
    h18 = [x for x in h18 if x][0]
    h18_sp = [float(h18.group(i)) for i in (1, 2, 3)]
    print(
        "[WF-10 H27] z_c spread per radius r = 15 / 30 / 45 mm: d2 "
        + " / ".join(f"{v * 100:.4f}" for v in rc2["spreads"])
        + " %; 359 (d1, 0.015 m) "
        + " / ".join(f"{v * 100:.4f}" for v in rc1["spreads"])
        + " %; H18 (d1, 0.0075 m) "
        + " / ".join(f"{v:.4f}" for v in h18_sp)
        + " %; H24 four-port "
        + " / ".join(f"{v:.4f}" for v in H24_FOUR_PORT_SPREAD_PCT)
        + " % (PREDICTED d2 12-22%)",
        flush=True,
    )
    s1, s2 = d1["cstar_sar"], d2["cstar_sar"]
    print(
        f"[WF-10 H27] c*: d1 {d1['cstar_label']} -> d2 {d2['cstar_label']} "
        f"({'unchanged' if d1['cstar_label'] == d2['cstar_label'] else 'MOVED'}); SAR10g(c*) "
        f"{s1:.9e} -> {s2:.9e} W/kg ({s2 / s1 - 1:+.4%}); per P_phantom "
        f"{s1 / p1_regate if p1_regate else np.nan:.6e} -> {s2 / p2:.6e} W/kg/W "
        f"({s2 / s1 * pw_ratio - 1:+.4%})",
        flush=True,
    )
    if p1_regate is None:
        print(
            f"[WF-10 H27] d1 ccw P_phantom NOT printed in {D1_PPH_LOG}: raw move only", flush=True
        )

    assert rc2["iiia"] <= C4_COVARIANCE_BAND, f"(iii-a) {rc2['iiia']:.4e} — never re-band"
    assert rc2["iiib"] <= C4_COVARIANCE_BAND, f"(iii-b) {rc2['iiib']:.4e} — never re-band"
    assert rc2["noflip"] > C4_COVARIANCE_BAND, f"no-flip {rc2['noflip']:.4e} at or below band"
