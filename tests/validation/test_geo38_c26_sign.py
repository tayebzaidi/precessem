"""GEO-38 C26 (§9 item 776, §10 chain C, lane 1 → chain step, no solve): the
circuit-vs-root sign, read on the circuit layer of the two saved degree-1
64×64 matrices.

Premise (corrected by the item): ``C_mode`` is ESR-invariant by construction —
a uniform r adds r·I to Z, which in the k domain adds r to every block's
eigenvalues, moving Re λ and leaving Im λ (so ``C_mode = 1/(ω Im λ₁)``)
alone.  (b) makes that a gate on both matrices.  What can carry the sign is
the circuit-layer root: Im Y_in at P33 under the (1, +1j) drive, every
non-feed gap terminated in Z_C (r = 0) or Z_RC (r = 0.2 Ω) and each feed
driven across its own (r + C) branch — C4's "circuit-only Y_in"
(``test_geo38_f_clinical_esr_root.py:626–629``), computed here by the same
helper, ``h17._current_route_y_in``, on the npz alone (no per-node field
data).

Measured:

* (a) ``C_mode``(k = 1) by C10's ``ring_c_modes`` at an explicit frequency;
  control: 195.0265 pF (C10 W1 ``:106``) and 43.9719 pF (C15 W1 ``:76``)
  reproduced at ``.4f`` (ASSERTED).
* (b) ``ring_mode_rows`` at ``esr_ohm`` 0.2 against 0, ``n_slots = 32`` — the
  two-ring 64-port reading (the brief's "n_slots = 64" names the port count;
  ``ring_mode_rows`` takes slots per ring and raises on 64 for a 64×64):
  |ΔIm λ₁|/|Im λ₁(0)| and |ΔRe λ₁ − r|/r ≤ ``REDUCTION_BAND`` (ASSERTED).
* (R) reciprocity / σ_max on both matrices (C10's bands, ASSERTED).
* (C4) the helper reproduces C4's printed circuit-only Y_in at the 195.0265 pF
  node (``20261026T070008Z_GEO-38-C4.log:13033–13034``) to 1e-8 rel
  (ASSERTED; the log prints 10 significant digits).
* (c) Im Y_in(P33, (1,+1j)) on 11 nodes over ±1 % about each ``C_mode``, per
  twin and per r; − → + crossing by the bracketing pair (10-31 ruling); none
  ⇒ widened once to ±3 % (31 nodes).  Printed beside the in-model roots
  (193.5996 / 44.0615 pF) with signed gaps; RECORD, never asserted.
* (c-ctl) the 1.5 T npz read at 128 MHz (wrong ω): no − → + crossing within
  1 % of the 64 MHz reference root (ASSERTED, per the brief's gate list).

Leak guard: node capacitances print at ``.4f`` only (≤ 7 significant
digits); no ``.12e`` node form is printed.  JSON → ``logs/geo38/c26-sign.json``
(gitignored).  Knob ``FEM_EM_GEO38_C26``: unset ⇒ the gated test skips; ``1``
⇒ runs, an absent npz fails (never regenerated); anything else ⇒
``ValueError``.  Not claimed: a root, a tuning statement, anything at
degree 2.
"""

from __future__ import annotations

import json
import time

import numpy as np
import pytest

import tests.validation.test_wf10_f_human_copper_esr_root as h17
from fem_em_solver.ports.circuit import s_to_z
from fem_em_solver.ports.ring_modes import ring_mode_rows
from tests.validation.test_geo38_f_clinical_3t import FREQUENCY_3T_HZ
from tests.validation.test_geo38_f_clinical_3t_matrix import C11_SAVE_PATH
from tests.validation.test_geo38_f_clinical_c_mode import (
    FREQUENCY_HZ,
    WORKSPACE,
    _knob,
    map_from_geometry,
    ring_c_modes,
)
from tests.validation.test_geo38_f_clinical_matrix import SAVE_PATH
from tests.validation.test_port_lumped_rlc_termination import REDUCTION_BAND
from tests.validation.test_port_package_sparameters import REFERENCE_IMPEDANCE_OHM

KNOB_C26 = "FEM_EM_GEO38_C26"
C26_JSON = WORKSPACE / "logs" / "geo38" / "c26-sign.json"
ESR_OHM = 0.2
FEEDS = ("P33", "P41")
FEEDS_SLOTS = [(0, "bottom"), (8, "bottom")]
SENSE = (1.0, 1.0j)
N_GRID_1 = 11  # ±1 %, ≥ 9 nodes
N_GRID_3 = 31  # ±3 %, the one widening
CTL_MOVE_MIN_REL = 0.01  # the brief's "> 1 %"
C4_REPRO_RTOL = 1.0e-8

# ---- record comparands (transcribed, never bands) ------------------------------
C_MODE_RECORD_STR = {"1.5T": "195.0265", "3T": "43.9719"}  # C10 W1 :106; C15 W1 :76
IN_MODEL_ROOT_PF = {"1.5T": 193.5996, "3T": 44.0615}  # entry 52's node; entry 49 (C16)
# C4 (entry 50) :13033 / :13034, node 195.0265 pF, P33 (1,+1j), circuit-only:
C4_Y_RECORD = {0.2: 4.684470286e-02 + 5.484324358e-02j, 0.0: 4.126229463e-02 + 2.607149565e-02j}
C4_NODE_PF = 195.0265

TWINS = {
    "1.5T": {"path": SAVE_PATH, "freq_hz": FREQUENCY_HZ},
    "3T": {"path": C11_SAVE_PATH, "freq_hz": FREQUENCY_3T_HZ},
}


def c26_enabled():
    return _knob(KNOB_C26)


def _say(msg):
    print(f"[GEO-38 C26] {msg}", flush=True)


def _load(path):
    full = WORKSPACE / path
    assert full.exists(), f"{full} absent: BLOCKED, never regenerated"
    with np.load(full, allow_pickle=False) as z:
        sv = {k: z[k] for k in z.files}
    s = np.array(sv["s64"], dtype=np.complex128)
    order = [str(x) for x in sv["order"]]
    assert s.shape == (64, 64), s.shape
    assert order == [f"P{i}" for i in range(33, 97)], order[:4]
    geo = map_from_geometry(sv["azimuth_deg"], sv["z_m"])
    assert len(set(geo)) == 64
    assert geo[0] == FEEDS_SLOTS[0] and geo[8] == FEEDS_SLOTS[1], (geo[0], geo[8])
    return full, s, order, geo


def _z_feeds_first(s, order):
    """Z of the npz in C4's order: [P33, P41] + others by number."""
    others = sorted((p for p in order if p not in FEEDS), key=lambda p: int(p[1:]))
    new = list(FEEDS) + others
    ix = [order.index(p) for p in new]
    return s_to_z(s[np.ix_(ix, ix)], REFERENCE_IMPEDANCE_OHM)


def _y_in_p33(z64, c_pf, freq_hz, r_ohm):
    z_load = h17._z_rc(c_pf * 1e-12, freq_hz, r_ohm=r_ohm)
    return h17._current_route_y_in(z64, [0, 1], z_load, z_load, SENSE, REFERENCE_IMPEDANCE_OHM)[0]


def _grid(center_pf, half_rel, n):
    return np.round(np.linspace(center_pf * (1 - half_rel), center_pf * (1 + half_rel), n), 4)


def _scan(z64, center_pf, freq_hz, r_ohm, half_rel, n):
    nodes = _grid(center_pf, half_rel, n)
    b = np.array([_y_in_p33(z64, c, freq_hz, r_ohm).imag for c in nodes])
    all_x = h17._sign_change_roots(nodes, b)
    up = [
        (a, bb, r) for a, bb, r in all_x if np.interp(a, nodes, b) < 0.0 < np.interp(bb, nodes, b)
    ]
    return nodes, b, all_x, up


def _root_reading(z64, center_pf, freq_hz, r_ohm, label):
    out = {}
    for half, n in ((0.01, N_GRID_1), (0.03, N_GRID_3)):
        nodes, b, all_x, up = _scan(z64, center_pf, freq_hz, r_ohm, half, n)
        _say(
            f"{label} r={r_ohm:g} grid +-{half * 100:g}% ({n} nodes): Im Y_in(P33,(1,+1j)) "
            + ", ".join(f"{c:.4f}:{v:+.6e}" for c, v in zip(nodes, b))
            + " S"
        )
        _say(
            f"{label} r={r_ohm:g} +-{half * 100:g}%: sign changes "
            + (", ".join(f"[{a:.4f}, {bb:.4f}] -> {r:.4f} pF" for a, bb, r in all_x) or "none")
            + "; - -> + "
            + (", ".join(f"{r:.4f}" for _a, _b, r in up) or "none")
        )
        out = {
            "half_rel": half,
            "nodes_pf": nodes.tolist(),
            "im_y_s": b.tolist(),
            "sign_changes_pf": [list(x) for x in all_x],
            "up_roots_pf": [r for *_x, r in up],
        }
        if up:
            break
    roots = out["up_roots_pf"]
    if len(roots) > 1:
        roots = [min(roots, key=lambda r: abs(r - center_pf))]
        _say(
            f"{label} r={r_ohm:g}: {len(out['up_roots_pf'])} - -> + crossings; nearest C_mode "
            f"read: {roots[0]:.4f} pF"
        )
    out["root_pf"] = roots[0] if roots else None
    return out


def _gap_pct(x, ref):
    return None if x is None else 100.0 * (x - ref) / ref


def test_c26_knob(monkeypatch):
    """(K) knob parsing: unset / empty ⇒ off, ``1`` ⇒ on, anything else raises."""
    monkeypatch.delenv(KNOB_C26, raising=False)
    assert c26_enabled() is False
    monkeypatch.setenv(KNOB_C26, " ")
    assert c26_enabled() is False
    monkeypatch.setenv(KNOB_C26, "1")
    assert c26_enabled() is True
    monkeypatch.setenv(KNOB_C26, "yes")
    with pytest.raises(ValueError):
        c26_enabled()


def test_geo38_c26_sign():
    """C26: (a), (b), (R), (C4), (c) and (c-ctl) on both saved 64×64s."""
    if not c26_enabled():
        pytest.skip(f"{KNOB_C26} unset")
    from tests.validation.test_port_birdcage_ring_column import _reciprocity_ratio
    from tests.validation.test_port_birdcage_ring_matrix import MATRIX_PASSIVITY_CEILING
    from tests.validation.test_port_lumped_sheet_sweep import RECIPROCITY_BAND

    t0 = time.perf_counter()
    z0 = REFERENCE_IMPEDANCE_OHM
    rec = {"esr_ohm": ESR_OHM, "z0_ohm": z0, "sense": "(1,+1j)", "feed": "P33", "twins": {}}
    fails = []
    z_by = {}
    for name, tw in TWINS.items():
        full, s, order, geo = _load(tw["path"])
        f = tw["freq_hz"]
        ratio = _reciprocity_ratio(s)
        sig = float(np.linalg.norm(s, 2))
        _say(
            f"{name}: {full} ({full.stat().st_size} B) at {f / 1e6:g} MHz; reciprocity "
            f"{ratio:.6e} (band {RECIPROCITY_BAND:.0e}, ASSERTED); sigma_max {sig:.12f} "
            f"(ceiling {MATRIX_PASSIVITY_CEILING}, ASSERTED)"
        )
        if not ratio <= RECIPROCITY_BAND:
            fails.append(f"(R) {name} reciprocity {ratio:.3e}")
        if not sig <= MATRIX_PASSIVITY_CEILING:
            fails.append(f"(R) {name} sigma_max {sig:.6f}")
        # (a)
        r = ring_c_modes(s, geo, freq_hz=f, n_slots=32, z0=z0)
        c_mode_pf = r["c_mode_f"] * 1e12
        ok_a = f"{c_mode_pf:.4f}" == C_MODE_RECORD_STR[name]
        _say(
            f"(a) {name} C_mode {r['c_mode_f']:.12e} F = {c_mode_pf:.4f} pF; record "
            f"{C_MODE_RECORD_STR[name]} ({'reproduced' if ok_a else 'NOT reproduced'} at .4f, ASSERTED)"
        )
        if not ok_a:
            fails.append(f"(a) {name} C_mode {c_mode_pf:.4f} != {C_MODE_RECORD_STR[name]}")
        # (b)
        kw = dict(
            port_map=geo,
            freq_hz=f,
            c_f=r["c_mode_f"],
            feeds=FEEDS_SLOTS,
            drives={"(1,+1j)": list(SENSE)},
            z0=z0,
            n_slots=32,
        )
        rows0 = ring_mode_rows(s, esr_ohm=0.0, **kw)
        rows_r = ring_mode_rows(s, esr_ohm=ESR_OHM, **kw)
        l0, lr = rows0["modes"]["1"]["lambda_anti"], rows_r["modes"]["1"]["lambda_anti"]
        d_im = abs(lr[1] - l0[1]) / abs(l0[1])
        d_re = abs((lr[0] - l0[0]) - ESR_OHM) / ESR_OHM
        c_rows0 = rows0["c_mode_f"] * 1e12
        c_rows_r = rows_r["c_mode_f"] * 1e12
        _say(
            f"(b) {name} lambda_1 r=0 {l0[0]:+.9e} {l0[1]:+.9e}j, r={ESR_OHM:g} {lr[0]:+.9e} "
            f"{lr[1]:+.9e}j Ohm; |dIm|/|Im| {d_im:.3e}, |dRe - r|/r {d_re:.3e} (band "
            f"{REDUCTION_BAND:g}, ASSERTED); C_mode rows r=0 {c_rows0:.4f}, r={ESR_OHM:g} "
            f"{c_rows_r:.4f} pF; rows vs ring_c_modes {abs(c_rows0 - c_mode_pf) / c_mode_pf:.3e}"
        )
        if not (d_im <= REDUCTION_BAND and d_re <= REDUCTION_BAND):
            fails.append(f"(b) {name} d_im {d_im:.3e} d_re {d_re:.3e}")
        z64 = _z_feeds_first(s, order)
        z_by[name] = z64
        # (C4) reproduction, 1.5 T only (C4 is the 1.5 T fixture)
        c4 = {}
        if name == "1.5T":
            for r_ohm, y_rec in C4_Y_RECORD.items():
                y = _y_in_p33(z64, C4_NODE_PF, f, r_ohm)
                rel = abs(y - y_rec) / abs(y_rec)
                _say(
                    f"(C4) node {C4_NODE_PF:.4f} pF r={r_ohm:g}: Y_in P33 {y:.9e} S vs C4 log "
                    f"{y_rec:.9e} S; rel {rel:.3e} (<= {C4_REPRO_RTOL:g}, ASSERTED)"
                )
                c4[str(r_ohm)] = rel
                if not rel <= C4_REPRO_RTOL:
                    fails.append(f"(C4) r={r_ohm:g} rel {rel:.3e}")
        # (c)
        readings = {}
        for r_ohm in (0.0, ESR_OHM):
            readings[str(r_ohm)] = _root_reading(z64, c_mode_pf, f, r_ohm, name)
        root_im = IN_MODEL_ROOT_PF[name]
        gap_im = _gap_pct(root_im, c_mode_pf)
        line = [f"in-model {root_im:.4f} pF ({gap_im:+.3f} %)"]
        for k, v in readings.items():
            g = _gap_pct(v["root_pf"], c_mode_pf)
            v["gap_pct"] = g
            line.append(
                f"circuit r={float(k):g} "
                + (
                    "no root within +-3 %"
                    if v["root_pf"] is None
                    else f"{v['root_pf']:.4f} pF ({g:+.3f} %)"
                )
            )
        _say(
            f"(c) {name} C_mode {c_mode_pf:.4f} pF; "
            + "; ".join(line)
            + " (signed gap (root - C_mode)/C_mode; RECORD, never asserted)"
        )
        rec["twins"][name] = {
            "npz": str(full),
            "freq_hz": f,
            "reciprocity": ratio,
            "sigma_max": sig,
            "c_mode_pf": c_mode_pf,
            "c_mode_record_reproduced": ok_a,
            "lambda1_r0": l0,
            "lambda1_r": lr,
            "esr_d_im": d_im,
            "esr_d_re": d_re,
            "c4_repro_rel": c4,
            "in_model_root_pf": root_im,
            "in_model_gap_pct": gap_im,
            "circuit": readings,
        }

    # (c-ctl): the 1.5 T npz at 128 MHz (wrong ω)
    t15 = rec["twins"]["1.5T"]
    ref = (
        t15["circuit"][str(ESR_OHM)]["root_pf"]
        or t15["circuit"]["0.0"]["root_pf"]
        or t15["c_mode_pf"]
    )
    ctl = {}
    moved_all = True
    for r_ohm in (0.0, ESR_OHM):
        nodes, b, all_x, up = _scan(
            z_by["1.5T"], t15["c_mode_pf"], FREQUENCY_3T_HZ, r_ohm, 0.03, N_GRID_3
        )
        near = [x for *_y, x in up if abs(x - ref) / ref <= CTL_MOVE_MIN_REL]
        moved_all &= not near
        _say(
            f"(c-ctl) 1.5T npz at {FREQUENCY_3T_HZ / 1e6:g} MHz r={r_ohm:g}, +-3 % grid: - -> + "
            + (", ".join(f"{x:.4f}" for *_y, x in up) or "none")
            + f"; within 1 % of the 64 MHz "
            f"reference {ref:.4f} pF: {len(near)} (ASSERTED 0 ⇒ moved > 1 %); Im Y_in range "
            f"{b.min():+.6e} .. {b.max():+.6e} S"
        )
        ctl[str(r_ohm)] = {"up_roots_pf": [x for *_y, x in up], "near": near}
    if not moved_all:
        fails.append("(c-ctl) the wrong-ω read did not move the root by > 1 %")
    rec["control_wrong_omega"] = {"reference_pf": ref, **ctl}

    # The sign reading (RECORD)
    def _reproduces(r_key):
        out = []
        for name, t in rec["twins"].items():
            g = t["circuit"][r_key]["gap_pct"]
            out.append(g is not None and np.sign(g) == np.sign(t["in_model_gap_pct"]))
        return all(out)

    any_root = any(
        t["circuit"][k]["root_pf"] is not None
        for t in rec["twins"].values()
        for k in ("0.0", str(ESR_OHM))
    )
    if _reproduces("0.0"):
        reading = "termination's (r = 0 already reproduces both in-model signs)"
    elif _reproduces(str(ESR_OHM)):
        reading = "circuit's (r = 0.2 reproduces both in-model signs; r = 0 does not)"
    elif not any_root:
        reading = "no circuit root within +-3 % (on either twin, either r)"
    else:
        reading = "fixture's (neither r reproduces both in-model signs)"
    rec["reading"] = reading
    _say(f"SIGN READING (RECORD, never asserted): {reading}")
    C26_JSON.parent.mkdir(parents=True, exist_ok=True)
    C26_JSON.write_text(json.dumps(rec, indent=1, default=str))
    _say(f"wrote {C26_JSON} ({C26_JSON.stat().st_size} B); body {time.perf_counter() - t0:.2f} s")
    assert not fails, fails
