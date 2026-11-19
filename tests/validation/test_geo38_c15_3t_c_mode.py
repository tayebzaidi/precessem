"""GEO-38 C15 (§9 item 655, §10 chain C15): the 3 T ``C_mode`` (k = 1) of
F-clinical-3T's degree-1 copper 64×64 at 128 MHz, and C13's class anatomy on
the same matrix — read from entry 47's saved npz
(``logs/geo38/c11-cu-s64-d1-128mhz.npz``, gitignored scratch, read only).

Knob ``FEM_EM_GEO38_C15`` (C10's ``_knob`` pattern, imported): unset / empty ⇒
the gated test skips; ``1`` ⇒ it runs and an absent input fails; anything
else ⇒ ``ValueError``.

Asserted (the brief's anchors and controls):

* **(RR3)** reciprocity (``_reciprocity_ratio``) and σ_max of the loaded 64×64
  string-equal entry 47's prints (``.6e`` / ``.12f``) and inside
  ``RECIPROCITY_BAND`` / ``MATRIX_PASSIVITY_CEILING``;
* **(CT3)** the 34 class lines rebuilt from the npz (C3's ``_class_table`` /
  ``_class_spreads``) equal ``20261027T070008Z_GEO-38-C11.log:14750–14783``
  under C13's ``_norm``, 34 / 34;
* **(Im)** Im λ₁ > 0 (C10's assertion);
* **(Map)** the geometry map is unique 64 / 64 and P33 / P41 land on C10's
  slots;
* **(CF)** negative control: the same matrix read at 64 MHz gives a ``C_mode``
  more than 10 % away from the 128 MHz reading (expected ratio 2: C = 1/(ω Im λ));
* **(CN)** negative control: P33's z sign flipped ⇒ at least one class changes
  n and (CT3) no longer matches 34 / 34.

Printed only: every class above ``ADJACENT_SPREAD_BAND`` with n, mean |S|,
absolute and relative spread, and its absolute spread over "same ring, 1
step"'s; a SUMMARY line with the largest ratio.  Nothing is asserted on the
ratio.  Writes ``/workspace/logs/geo38/c15-c-mode-d1-128mhz.json`` with C10's
keys (entry 49's grid centre) after the anchors hold.

No solve, no band, no glyph.  Not claimed: a root (entry 49), any tuning
statement, whether the classes come from physics or the mesh, any AED value.
"""

from __future__ import annotations

import json
import time

import numpy as np
import pytest

from fem_em_solver.ports.circuit import z_to_s
from tests.validation.test_geo38_c13_class_anatomy import _abs_spreads, _norm, _table_lines
from tests.validation.test_geo38_f_clinical_3t import FREQUENCY_3T_HZ
from tests.validation.test_geo38_f_clinical_3t_matrix import C11_SAVE_PATH
from tests.validation.test_geo38_f_clinical_c_mode import (
    C3_KEYS,
    FREQUENCY_HZ,
    PRED_C3_P33_P41_SLOTS,
    REFERENCE_IMPEDANCE_OHM,
    WORKSPACE,
    _block_circulant,
    _canonical_map,
    _knob,
    _shuffled,
    _sym_blocks,
    map_from_geometry,
    map_from_names_64,
    ring_c_modes,
)
from tests.validation.test_geo38_f_clinical_matrix import (
    N_AZIMUTH_CLASSES_32,
    N_SYMMETRY_CLASSES_32,
    _azimuth_class_32,
    _class_spreads,
    _class_table,
    _expected_class_sizes,
)
from tests.validation.test_port_birdcage_lumped_column import ADJACENT_SPREAD_BAND
from tests.validation.test_port_birdcage_ring_column import _reciprocity_ratio
from tests.validation.test_port_birdcage_ring_matrix import MATRIX_PASSIVITY_CEILING
from tests.validation.test_port_lumped_sheet_sweep import RECIPROCITY_BAND

KNOB_C15 = "FEM_EM_GEO38_C15"
C15_JSON = WORKSPACE / "logs" / "geo38" / "c15-c-mode-d1-128mhz.json"
CF_MIN_REL = 0.10  # (CF): the brief's "> 10 %" (expected ratio 2)

# ---- records (not bands), each with its log line ------------------------------
# Entry 47's own prints: 20261027T070008Z_GEO-38-C11.log:14748–14749.
RECIPROCITY_3T_RECORD_STR = "8.858135e-15"
SIGMA_MAX_3T_RECORD_STR = "0.999805773843"

# Entry 47's printed 34-class table, verbatim payloads after the "[GEO-38 C3] "
# prefix: 20261027T070008Z_GEO-38-C11.log:14750–14783.  A record, not a band.
C11_TABLE_LINES = (
    "     same  ring,  0 steps  n =   64  mean |S| = 5.312882206e-01  spread   0.1444%",
    "     same  ring,  1 steps  n =  128  mean |S| = 3.736487718e-01  spread   0.1383%",
    "     same  ring,  2 steps  n =  128  mean |S| = 2.293004228e-01  spread   0.1806%",
    "     same  ring,  3 steps  n =  128  mean |S| = 1.170521811e-01  spread   0.3130%",
    "     same  ring,  4 steps  n =  128  mean |S| = 5.224440149e-02  spread   0.5223%",
    "     same  ring,  5 steps  n =  128  mean |S| = 3.010408181e-02  spread   0.8987%   <-- worst",
    "     same  ring,  6 steps  n =  128  mean |S| = 2.634068527e-02  spread   0.6096%",
    "     same  ring,  7 steps  n =  128  mean |S| = 2.187722578e-02  spread   0.4335%",
    "     same  ring,  8 steps  n =  128  mean |S| = 1.585856530e-02  spread   0.4064%",
    "     same  ring,  9 steps  n =  128  mean |S| = 1.036206401e-02  spread   0.4031%",
    "     same  ring, 10 steps  n =  128  mean |S| = 6.407567659e-03  spread   0.4743%",
    "     same  ring, 11 steps  n =  128  mean |S| = 3.980183241e-03  spread   0.5644%",
    "     same  ring, 12 steps  n =  128  mean |S| = 2.542537359e-03  spread   0.6209%",
    "     same  ring, 13 steps  n =  128  mean |S| = 1.615401678e-03  spread   0.6547%",
    "     same  ring, 14 steps  n =  128  mean |S| = 1.118375634e-03  spread   0.3962%",
    "     same  ring, 15 steps  n =  128  mean |S| = 1.076503237e-03  spread   0.2190%",
    "     same  ring, 16 steps  n =   64  mean |S| = 1.134943168e-03  spread   0.2012%",
    "     other ring,  0 steps  n =   64  mean |S| = 4.584144285e-01  spread   0.1279%",
    "     other ring,  1 steps  n =  128  mean |S| = 9.942989475e-02  spread   0.5840%",
    "     other ring,  2 steps  n =  128  mean |S| = 7.978403455e-02  spread   0.5603%",
    "     other ring,  3 steps  n =  128  mean |S| = 1.036100206e-01  spread   0.2771%",
    "     other ring,  4 steps  n =  128  mean |S| = 8.766615949e-02  spread   0.2362%",
    "     other ring,  5 steps  n =  128  mean |S| = 6.003921031e-02  spread   0.3147%",
    "     other ring,  6 steps  n =  128  mean |S| = 3.574433450e-02  spread   0.4306%",
    "     other ring,  7 steps  n =  128  mean |S| = 1.967399298e-02  spread   0.5050%",
    "     other ring,  8 steps  n =  128  mean |S| = 1.140637107e-02  spread   0.5745%",
    "     other ring,  9 steps  n =  128  mean |S| = 7.895836683e-03  spread   0.4356%",
    "     other ring, 10 steps  n =  128  mean |S| = 5.969981602e-03  spread   0.4231%",
    "     other ring, 11 steps  n =  128  mean |S| = 4.349922215e-03  spread   0.4436%",
    "     other ring, 12 steps  n =  128  mean |S| = 2.902100637e-03  spread   0.4715%",
    "     other ring, 13 steps  n =  128  mean |S| = 1.736975937e-03  spread   0.5668%",
    "     other ring, 14 steps  n =  128  mean |S| = 1.060061571e-03  spread   0.5498%",
    "     other ring, 15 steps  n =  128  mean |S| = 9.961472550e-04  spread   0.3370%",
    "     other ring, 16 steps  n =   64  mean |S| = 1.071518656e-03  spread   0.2693%",
)
WORST_MARKER = "<-- worst"


def c15_enabled():
    return _knob(KNOB_C15)


def _say(msg):
    print(f"[GEO-38 C15] {msg}", flush=True)


def _ct_misses(lines):
    if len(lines) != len(C11_TABLE_LINES):
        return [("<length>", f"{len(lines)} != {len(C11_TABLE_LINES)}")]
    return [(a, b) for a, b in zip(lines, C11_TABLE_LINES) if _norm(a) != _norm(b)]


# --------------------------------------------------------------------------- W0
def test_c15_knob_and_records(monkeypatch):
    """Ungated: the knob, the record table's shape, and the frequency constants."""
    monkeypatch.delenv(KNOB_C15, raising=False)
    assert c15_enabled() is False
    monkeypatch.setenv(KNOB_C15, "1")
    assert c15_enabled() is True
    for bad in ("0", "yes", "build"):
        monkeypatch.setenv(KNOB_C15, bad)
        with pytest.raises(ValueError):
            c15_enabled()
    monkeypatch.delenv(KNOB_C15, raising=False)
    assert len(C11_TABLE_LINES) == N_SYMMETRY_CLASSES_32
    assert sum(WORST_MARKER in ln for ln in C11_TABLE_LINES) == 1
    assert FREQUENCY_3T_HZ == 2.0 * FREQUENCY_HZ, (FREQUENCY_3T_HZ, FREQUENCY_HZ)


def test_c15_reader_frequency_closed_form():
    """Ungated (CF)'s coding twin: a synthetic N = 32 closed form built for
    ``c_mode`` at 128 MHz reads back ``c_mode`` at 128 MHz and 2·``c_mode`` at
    64 MHz, ports shuffled."""
    n = 32
    c_true = 1.0e-10
    omega = 2.0 * np.pi * FREQUENCY_3T_HZ
    z = _block_circulant(_sym_blocks(n, omega, c_true), n)
    zs, pm = _shuffled(z, n, seed=15)
    s = z_to_s(zs, REFERENCE_IMPEDANCE_OHM)
    r128 = ring_c_modes(s, pm, freq_hz=FREQUENCY_3T_HZ, n_slots=n, z0=REFERENCE_IMPEDANCE_OHM)
    r64 = ring_c_modes(s, pm, freq_hz=FREQUENCY_HZ, n_slots=n, z0=REFERENCE_IMPEDANCE_OHM)
    rel128 = abs(r128["c_mode_f"] / c_true - 1.0)
    ratio = r64["c_mode_f"] / r128["c_mode_f"]
    print(
        f"\n[GEO-38 C15] (K) closed form N=32: C_mode(128 MHz) rel {rel128:.3e}; "
        f"C_mode(64)/C_mode(128) {ratio:.12f}",
        flush=True,
    )
    assert rel128 <= 1e-9, rel128
    assert abs(ratio - 2.0) <= 1e-12, ratio
    assert pm != _canonical_map(n)


# --------------------------------------------------------------------------- W1
def test_geo38_c15_3t_c_mode():
    if not c15_enabled():
        pytest.skip(f"{KNOB_C15} unset (reads gitignored logs/geo38 npz)")
    t0 = time.perf_counter()
    path = WORKSPACE / C11_SAVE_PATH
    assert path.exists(), f"{path} absent: entry 47's matrix; never re-solve"
    with np.load(path, allow_pickle=False) as z:
        sv = {k: z[k] for k in z.files}
    missing = [k for k in C3_KEYS if k not in sv]
    assert not missing, f"C11 npz lacks {missing}"
    s = np.array(sv["s64"], dtype=np.complex128)
    order = [str(x) for x in sv["order"]]
    az = np.asarray(sv["azimuth_deg"], dtype=float)
    zc = np.asarray(sv["z_m"], dtype=float)
    _say(
        f"{path} ({path.stat().st_size} B); s64 {s.shape} {s.dtype}; order {order[0]}…{order[-1]}; "
        f"cells {int(sv['cells'])}, n_ranks {int(sv['n_ranks'])}"
    )
    assert s.shape == (64, 64), s.shape
    assert order == [f"P{i}" for i in range(33, 97)], order[:4]

    # ---- (RR3) ----
    ratio = _reciprocity_ratio(s)
    sigma_max = float(np.linalg.norm(s, 2))
    rr = (f"{ratio:.6e}", f"{sigma_max:.12f}")
    _say(
        f"(RR3) reciprocity {rr[0]} vs entry 47 {RECIPROCITY_3T_RECORD_STR} (band {RECIPROCITY_BAND:.0e}); "
        f"sigma_max {rr[1]} vs entry 47 {SIGMA_MAX_3T_RECORD_STR} (ceiling {MATRIX_PASSIVITY_CEILING}) "
        f"(string equality, ASSERTED)"
    )

    # ---- (CT3) ----
    classes, res = _class_table(az, zc, classify=_azimuth_class_32)
    sizes = {k: len(v) for k, v in classes.items()}
    spreads = _class_spreads(s, classes)
    mine = _table_lines(spreads, N_AZIMUTH_CLASSES_32)
    ct_miss = _ct_misses(mine)
    _say(
        f"(CT3) classes {len(classes)} (expected {N_SYMMETRY_CLASSES_32}), sizes equal expected "
        f"{sizes == _expected_class_sizes()}, worst pair residual {res:.3e} deg; lines matching "
        f"C11.log:14750–14783 {len(C11_TABLE_LINES) - len(ct_miss)}/{len(C11_TABLE_LINES)} "
        f"(normalisation: payload after the '[GEO-38 C3] ' prefix, '{WORST_MARKER}' dropped, "
        f"whitespace runs collapsed; ASSERTED)"
    )
    for a, b in ct_miss:
        _say(f"(CT3) MISS mine '{_norm(a)}' vs C11 '{_norm(b)}'")

    # ---- (Map) ----
    geo = map_from_geometry(az, zc)
    agree = sum(g == m for g, m in zip(geo, map_from_names_64(order)))
    p33_p41 = (geo[0], geo[8])
    _say(
        f"(Map) geometry map unique {len(set(geo))}/64; name-route agreement {agree}/64; P33 -> {geo[0]}, "
        f"P41 -> {geo[8]} (C10's {PRED_C3_P33_P41_SLOTS}; ASSERTED)"
    )

    # ---- C_mode at 128 MHz, and (CF) at 64 MHz ----
    r = ring_c_modes(s, geo, freq_hz=FREQUENCY_3T_HZ, n_slots=32, z0=REFERENCE_IMPEDANCE_OHM)
    lam1 = r["modes"][1]["lambda_anti"]
    for k in range(17):
        m = r["modes"][k]
        lam = m["lambda_anti"]
        _say(
            f"  k={k:2d} lambda_anti {lam.real:+.9e} {lam.imag:+.9e}j Ohm anti {m['anti']:.6f} "
            f"C {m['c_anti_f'] * 1e12:+.4f} pF"
        )
    _say(
        f"C_mode {r['c_mode_f']:.12e} F = {r['c_mode_f'] * 1e12:.4f} pF at {FREQUENCY_3T_HZ / 1e6:g} MHz; "
        f"C_2 {r['c_k2_f']:.12e} F = {r['c_k2_f'] * 1e12:.4f} pF; anti(k=1) {r['modes'][1]['anti']:.6f}; "
        f"off-block {r['off_block']:.6e}; (Im) Im lambda_1 {lam1.imag:+.9e} Ohm (ASSERTED > 0)"
    )
    r64 = ring_c_modes(s, geo, freq_hz=FREQUENCY_HZ, n_slots=32, z0=REFERENCE_IMPEDANCE_OHM)
    cf_ratio = r64["c_mode_f"] / r["c_mode_f"]
    _say(
        f"(CF) C_mode read at {FREQUENCY_HZ / 1e6:g} MHz {r64['c_mode_f'] * 1e12:.4f} pF; ratio 64/128 "
        f"{cf_ratio:.6f} (expected 2; ASSERTED |ratio - 1| > {CF_MIN_REL})"
    )

    # ---- (CN) negative control: P33's z sign flipped ----
    zf = zc.copy()
    zf[0] = -zf[0]
    cls_f, _r = _class_table(az, zf, classify=_azimuth_class_32)
    n_changed = sum(1 for k in classes if len(cls_f.get(k, [])) != len(classes[k]))
    lines_f = (
        _table_lines(_class_spreads(s, cls_f), N_AZIMUTH_CLASSES_32)
        if len(cls_f) == len(classes)
        else []
    )
    cn_miss = len(_ct_misses(lines_f))
    _say(
        f"(CN) {order[0]} z sign flipped: classes {len(cls_f)}, classes with changed n {n_changed} "
        f"(ASSERTED >= 1), (CT3) lines differing {cn_miss} (ASSERTED > 0)"
    )

    assert rr == (RECIPROCITY_3T_RECORD_STR, SIGMA_MAX_3T_RECORD_STR), rr
    assert ratio <= RECIPROCITY_BAND, ratio
    assert sigma_max <= MATRIX_PASSIVITY_CEILING, sigma_max
    assert len(classes) == N_SYMMETRY_CLASSES_32 and sizes == _expected_class_sizes(), sizes
    assert not ct_miss, ct_miss[:3]
    assert len(geo) == 64 and len(set(geo)) == 64, sorted(geo)
    assert p33_p41 == PRED_C3_P33_P41_SLOTS, p33_p41
    assert lam1.imag > 0.0, lam1
    assert abs(cf_ratio - 1.0) > CF_MIN_REL, cf_ratio
    assert n_changed >= 1 and cn_miss > 0, (n_changed, cn_miss)

    # ---- entry 49's grid centre (written only after the anchors hold) ----
    C15_JSON.parent.mkdir(parents=True, exist_ok=True)
    C15_JSON.write_text(
        json.dumps(
            {
                "source": str(path),
                "freq_hz": FREQUENCY_3T_HZ,
                "z0_ohm": REFERENCE_IMPEDANCE_OHM,
                "n_slots": 32,
                "c_mode_f": r["c_mode_f"],
                "c_k2_f": r["c_k2_f"],
                "off_block": r["off_block"],
                "reciprocity": ratio,
                "sigma_max": sigma_max,
                "name_route_agreement": agree,
                "cells": int(sv["cells"]),
                "n_ranks": int(sv["n_ranks"]),
                "modes": {
                    str(k): {
                        "lambda_anti": [m["lambda_anti"].real, m["lambda_anti"].imag],
                        "anti": m["anti"],
                        "c_anti_f": m["c_anti_f"],
                    }
                    for k, m in r["modes"].items()
                },
            },
            indent=2,
        )
    )
    _say(f"wrote {C15_JSON} ({C15_JSON.stat().st_size} B)")

    # ---- the anatomy (printed only; entry 47's decision-rule action) ----
    abs_sp = _abs_spreads(s, classes)
    ref_key = (True, 1)
    ref_abs = abs_sp[ref_key]
    _say(
        f"reference class 'same ring, 1 step': rel {spreads[ref_key][0] * 100:.4f}%, abs {ref_abs:.6e}"
    )
    over = [
        k
        for same in (True, False)
        for k in [(same, st) for st in range(N_AZIMUTH_CLASSES_32)]
        if spreads[k][0] > ADJACENT_SPREAD_BAND
    ]
    ratios = {}
    for key in over:
        sp, cnt, mean = spreads[key]
        ratios[key] = abs_sp[key] / ref_abs
        _say(
            f"anatomy {'same ' if key[0] else 'other'} ring, {key[1]:2d} steps  n = {cnt:4d}  "
            f"mean |S| {mean:.9e}  abs spread {abs_sp[key]:.6e}  rel {sp * 100:8.4f}%  "
            f"abs / (same ring, 1 step) {ratios[key]:.4f}"
        )
    worst = max(ratios, key=ratios.get)
    n_gt1 = sum(v > 1.0 for v in ratios.values())
    _say(
        f"SUMMARY classes above ADJACENT_SPREAD_BAND {ADJACENT_SPREAD_BAND}: {len(over)}/34; largest "
        f"abs-spread ratio {ratios[worst]:.4f} on {'same' if worst[0] else 'other'} ring, {worst[1]} steps; "
        f"classes with ratio > 1: {n_gt1} (printed, not asserted)"
    )
    _say(f"body {time.perf_counter() - t0:.2f} s")
