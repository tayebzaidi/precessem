"""`ANS-8` H25 (§9 item 516, §10 chain H): SPEC v2 rows on the saved degree-2 copper 32x32.

No solve.  Entry 36's saved matrix ``logs/wf12/h16b-cu-s32-d2.npz`` (degree 2,
copper via Leontovich, kappa off, the 369 387-cell hole mesh, 64 MHz) is run
through the one shared ``ring_mode_rows_v2`` (via ``ans:8``'s ``rows_v2``),
beside degree-1 copper (``h12c-cu-s32-a2.npz``) and the two kappa-off
sigma = 800 matrices (``logs/wf7/``), each tied to its window by reproduction.
The ESR uniform-shift identity (`PORT-24`'s (S)) is asserted on the new matrix,
which is what licenses quoting the lossless ``C_mode`` beside the r = 0.2 ohm
root.

This module prints our numbers only: no band, tolerance or target is set; every
band is imported by name.  The record constants below are comparands (the
361-3 precedent), each with its log line.  Every 2b share and leak is the
**circuit route** (Norton feeds), never the in-model (field-route) row 2.

Knob ``FEM_EM_ANS8_H25``: unset => the gated test skips; ``"1"`` => runs, and an
absent degree-2 npz **fails**; anything else => ``ValueError``.
"""

from __future__ import annotations

import importlib.util
import os
import sys
import time
from pathlib import Path

import numpy as np
import pytest

KNOB = "FEM_EM_ANS8_H25"
REPO = Path(__file__).resolve().parents[2]
CASE = REPO / "examples" / "ansys_benchmarks" / "ans8_f_human_quadrature_mode_content_64MHz"
SCRIPT = CASE / "08_f_human_quadrature_mode_content.py"
FREQ_HZ = 64.0e6
N_PORTS = 32
N_PER_RING = 16
H12C_A2_NPZ = "/workspace/logs/wf12/h12c-cu-s32-a2.npz"
CIRCUIT_LABEL = (
    "circuit route (Norton feeds); on copper the circuit <-> field identity at "
    "C_mode is the pre-registered H12c record (1.581167); not the in-model row"
)

# ---- record comparands (not bands), each with its log line ------------------
# Entry 36's window: docs/testing/logs/20261018T070008Z_WF-12-step5-H16c.log.
SIGMA_MAX_CU_D2_RECORD = 0.999983227022  # :8036
S1717_CU_D2_RECORD = 0.402723 + 0.391780j  # :8057 (printed beside)
RECIPROCITY_CU_D2_RECORD = 9.964981e-14  # :8035 (printed beside)
CLASS_MEAN_CU_D2_RECORD = {
    (True, 0): 5.614468065e-01,  # :8038
    (True, 1): 2.075482324e-01,  # :8039
    (True, 2): 1.028625914e-01,  # :8040
    (True, 3): 6.953867595e-02,  # :8041
    (True, 4): 6.146862913e-02,  # :8042
    (True, 5): 6.110272147e-02,  # :8043
    (True, 6): 6.237623938e-02,  # :8044
    (True, 7): 6.349044468e-02,  # :8045
    (True, 8): 6.390400628e-02,  # :8046
    (False, 0): 6.791200966e-01,  # :8047
    (False, 1): 6.952027509e-02,  # :8048
    (False, 2): 4.169587581e-02,  # :8049
    (False, 3): 5.801530821e-02,  # :8050
    (False, 4): 6.191764103e-02,  # :8051
    (False, 5): 6.317682134e-02,  # :8052
    (False, 6): 6.407670910e-02,  # :8053
    (False, 7): 6.476908015e-02,  # :8054
    (False, 8): 6.503497059e-02,  # :8055
}
# The same lines' d1 column: h12c-cu-s32.npz read by name.
CLASS_MEAN_CU_D1_RECORD = {
    (True, 0): 5.475316379e-01,
    (True, 1): 2.078336513e-01,
    (True, 2): 1.033277230e-01,
    (True, 3): 7.022870785e-02,
    (True, 4): 6.204664148e-02,
    (True, 5): 6.156144157e-02,
    (True, 6): 6.278578084e-02,
    (True, 7): 6.388862094e-02,
    (True, 8): 6.430136797e-02,
    (False, 0): 6.957296292e-01,
    (False, 1): 6.632270829e-02,
    (False, 2): 4.298222110e-02,
    (False, 3): 5.884048498e-02,
    (False, 4): 6.252567795e-02,
    (False, 5): 6.369516838e-02,
    (False, 6): 6.455612010e-02,
    (False, 7): 6.523454123e-02,
    (False, 8): 6.549744183e-02,
}
# H12c attempt 2: docs/testing/logs/20261006T170222Z_WF-12-step3-a2-W1.log.gz (decompressed lines).
C_MODE_CU_D1_RECORD_F = 1.181185799070e-10  # :6366
OFF_BLOCK_CU_D1_RECORD = 4.195925e-03  # :6366
P17_2B_CU_D1_RECORD = (0.993012, 0.005612, 0.998482)  # :6382 feed share, feed leak, top share
CW_2B_CU_D1_RECORD = (0.993217, 0.998703, 0.005618, 0.000013)  # :6387 share b/t, leak b/t
# Item 366: docs/testing/logs/20261011T123443Z_ANS-8-step1c-W1.log:116.
C_MODE_S800_RECORD_F = {2: 1.082287732e-10, 1: 1.131081186e-10}
ROOT_D1_PF = 118.2437  # the degree-1 ESR root (SPEC v3), for the printed offset only

# ---- predictions (written first; printed, never asserted — rule (e)) --------
PRED_CMODE_CU_D2_PF_WINDOW = (114.37, 116.68)  # within 1 % of the d2 root
PRED_CMODE_POINT_OFFSET_PF = 115.399  # d1 C_mode/root offset carried to the d2 root
PRED_CMODE_POINT_S800_PF = 113.023  # sigma = 800's -4.314 % order move on 118.1186
PRED_ANTI_FLOOR = 0.99
PRED_CCW_SHARE_FLOOR = 0.99
PRED_CCW_LEAK = (0.002, 0.02)
PRED_D2_COLUMN_CONTROL = (2.0e-4, 4.0e-4)
PRED_PLANT_MISS = 3.0e-2


def _enabled() -> bool:
    raw = os.environ.get(KNOB, "").strip()
    if not raw:
        return False
    if raw != "1":
        raise ValueError(f"{KNOB} must be unset or '1', got {raw!r}")
    return True


def _need():
    if not _enabled():
        pytest.skip(f"{KNOB} unset (H25 reads gitignored logs/wf12, logs/wf7 npz)")


def _load_ans8():
    spec = importlib.util.spec_from_file_location("ans8_h25", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["ans8_h25"] = mod
    spec.loader.exec_module(mod)
    return mod


def _names_az(order):
    """``(az, zc)`` from SPEC names alone (the d2 file has no slot keys)."""
    from fem_em_solver.ports.ring_modes import ring_port_map_from_names
    from tests.validation.test_port_birdcage_ring_matrix import AZIMUTH_STEP_DEG

    pm = ring_port_map_from_names(order)
    az = np.array([j * AZIMUTH_STEP_DEG for j, _r in pm], dtype=float)
    zc = np.array([1.0 if r == "top" else -1.0 for _j, r in pm], dtype=float)
    return az, zc


def _slot_az(saved):
    from tests.validation.test_port_birdcage_ring_matrix import AZIMUTH_STEP_DEG

    return (
        np.asarray(saved["slot_k"], dtype=float) * AZIMUTH_STEP_DEG,
        np.asarray(saved["slot_sign"], dtype=float),
    )


def _rels(spreads, record):
    return {k: abs(spreads[k][2] / record[k] - 1.0) for k in record}


def _esr_shift(s, port_map, c_f, feeds, esr):
    """``(worst |dRe - r|/r, worst |dIm|/|Im|, |dC|/C, per-k dRe)`` for ``esr`` vs 0."""
    from fem_em_solver.ports.ring_modes import ring_mode_rows
    from tests.validation.test_wf10_f_human_copper_esr_root import ESR_OHM

    kw = dict(
        port_map=port_map, freq_hz=FREQ_HZ, c_f=c_f, feeds=feeds, drives={"p17_alone": (1.0, 0.0)}
    )
    r0 = ring_mode_rows(s, **kw)
    rr = ring_mode_rows(s, esr_ohm=esr, **kw)
    d_re, d_im = [], []
    for k in range(N_PER_RING // 2 + 1):
        l0 = complex(*r0["modes"][str(k)]["lambda_anti"])
        lr = complex(*rr["modes"][str(k)]["lambda_anti"])
        d_re.append(abs((lr.real - l0.real) - ESR_OHM) / ESR_OHM)
        d_im.append(abs(lr.imag - l0.imag) / abs(l0.imag))
    dc = abs(rr["c_mode_f"] / r0["c_mode_f"] - 1.0)
    return max(d_re), max(d_im), dc, d_re


# --------------------------------------------------------------------------- W0
def test_ans8_h25_coding(monkeypatch):
    """Ungated: knob parser; names = slot synthesis (shuffled D16); (Sigma k); (S)."""
    from fem_em_solver.ports.ring_modes import k_spectrum, ring_port_map_from_names
    from tests.unit.test_port24_esr_circuit import SHIFT_RTOL
    from tests.validation.test_wf7_f_human_ring_matrix import _class_means

    monkeypatch.delenv(KNOB, raising=False)
    assert _enabled() is False
    monkeypatch.setenv(KNOB, "")
    assert _enabled() is False
    monkeypatch.setenv(KNOB, "1")
    assert _enabled() is True
    for bad in ("0", "yes", "2"):
        monkeypatch.setenv(KNOB, bad)
        with pytest.raises(ValueError):
            _enabled()

    # Synthetic D16 32x32 under a shuffled port order.
    rng = np.random.default_rng(516)
    order = [f"P{n}" for n in rng.permutation(np.arange(17, 49))]
    pm = ring_port_map_from_names(order)
    slot_k = np.array([j for j, _r in pm])
    slot_sign = np.array([1.0 if r == "top" else -1.0 for _j, r in pm])
    val = {
        (same, d): complex(rng.uniform(0.01, 1.0), rng.uniform(-1.0, 1.0))
        for same in (True, False)
        for d in range(N_PER_RING // 2 + 1)
    }
    s = np.zeros((N_PORTS, N_PORTS), dtype=np.complex128)
    for i in range(N_PORTS):
        for j in range(N_PORTS):
            d = abs(int(slot_k[i]) - int(slot_k[j]))
            d = min(d, N_PER_RING - d)
            s[i, j] = val[(bool(slot_sign[i] * slot_sign[j] > 0), d)]
    c_n, sp_n, res_n = _class_means(s, *_names_az(order))
    c_s, sp_s, res_s = _class_means(s, *_slot_az({"slot_k": slot_k, "slot_sign": slot_sign}))
    print(
        f"[ANS-8 H25 W0] shuffled D16: names {len(sp_n)} classes, slot {len(sp_s)}; "
        f"bitwise {sp_n == sp_s}; max spread {max(v[0] for v in sp_n.values()):.3e}",
        flush=True,
    )
    assert len(sp_n) == 18 and sp_n == sp_s
    assert all(v[0] == 0.0 for v in sp_n.values())

    # (Sigma k) on a synthetic current vector.
    ks = k_spectrum(rng.normal(size=16) + 1j * rng.normal(size=16))
    fr = ks["fractions"]
    sk = abs(sum(fr) - 1.0)
    print(f"[ANS-8 H25 W0] (Sigma k) |sum - 1| {sk:.3e}", flush=True)
    assert sk <= 1e-12
    assert ks["share_k1"] == fr[1] + fr[15] and ks["leak_k0"] == fr[0]

    # (S) on a synthetic passive 32x32 (canonical order, P17..P48).
    names = [f"P{n}" for n in range(17, 49)]
    a = rng.normal(size=(N_PORTS, N_PORTS)) + 1j * rng.normal(size=(N_PORTS, N_PORTS))
    a = a + a.T  # reciprocal
    s_syn = 0.9 * a / np.linalg.norm(a, 2)
    pm_c = ring_port_map_from_names(names)
    w_re, w_im, w_c, _ = _esr_shift(s_syn, pm_c, 1.0e-10, [(0, "bottom"), (4, "bottom")], 0.2)
    print(
        f"[ANS-8 H25 W0] (S) synthetic: dRe {w_re:.3e}, dIm {w_im:.3e}, dC {w_c:.3e} "
        f"(band {SHIFT_RTOL:.0e})",
        flush=True,
    )
    assert max(w_re, w_im, w_c) <= SHIFT_RTOL


# --------------------------------------------------------------------------- W1
def test_ans8_h25_copper_degree2_rows():
    _need()
    t0 = time.perf_counter()
    from fem_em_solver.ports.ring_modes import ring_mode_rows_v2, ring_port_map_from_names
    from tests.unit.test_port24_esr_circuit import SHIFT_RTOL
    from tests.validation.test_ans8_degree2_rows import SIGMA_MAX_D2_RECORD_BOXED  # noqa: F401
    from tests.validation.test_ans8_degree2_rows import (  # noqa: F401  (item 366's comparand, imported per brief)
        CLASS_MEAN_D2_RECORD_BOXED,
        _summary,
    )
    from tests.validation.test_port_birdcage_ring_column import (
        CONTROL_COLUMN_SCALE,
        _reciprocity_ratio,
    )
    from tests.validation.test_port_birdcage_ring_matrix import (
        MATRIX_PASSIVITY_CEILING,
        REUSE_REPRODUCTION_RTOL,
    )
    from tests.validation.test_port_lumped_sheet_sweep import RECIPROCITY_BAND
    from tests.validation.test_wf7_f_human_ring_matrix import (
        BOXED_NPZ_BASENAMES,
        CLASS_MEAN_D1_RECORD_BOXED,
        SAVE_CONTROL_FACTOR,
        SAVE_ROOT,
        _class_means,
    )
    from tests.validation.test_wf10_f_human_copper_esr_root import C_MODE_CU_PF, ESR_OHM
    from tests.validation.test_wf10_f_human_tuned_b1_sar import C_IN_F_HUMAN_D2_PF
    from tests.validation.test_wf12_f_human_copper_matrix import HOLE_CELL_RECORD
    from tests.validation.test_wf12_f_human_copper_s32_degree2 import (
        ENTRY36_SAVE,
        H12C_NPZ,
        UNKNOWNS_RECORD,
        _load_s32,
    )

    ans8 = _load_ans8()
    fails = []

    def say(msg):
        print(f"[ANS-8 H25] {msg}", flush=True)

    def check(label, ok, msg):
        say(f"{label} {msg} — {'green' if ok else 'RED'}")
        if not ok:
            fails.append(label)

    assert os.path.exists(ENTRY36_SAVE), f"{ENTRY36_SAVE} absent: entry 36's matrix; never re-solve"
    for p in (H12C_NPZ, H12C_A2_NPZ):
        assert os.path.exists(p), f"{p} absent"

    # ---- (F) d2 copper file identity ----
    d2 = _load_s32(ENTRY36_SAVE)
    a2 = _load_s32(H12C_A2_NPZ)
    d1 = _load_s32(H12C_NPZ)
    say(f"d2 keys {sorted(d2)}")
    s2 = np.array(d2["s32"], dtype=np.complex128)
    o2 = d2["order"]
    pm2 = ring_port_map_from_names(o2)
    bij = (
        len(set(pm2)) == N_PORTS
        and {r for _j, r in pm2} == {"bottom", "top"}
        and all(
            sorted(j for j, r in pm2 if r == ring) == list(range(N_PER_RING))
            for ring in ("bottom", "top")
        )
    )
    f_ok = (
        int(d2["degree"]) == 2
        and bool(d2["kappa_applied"]) is False
        and int(d2["cells"]) == HOLE_CELL_RECORD
        and int(d2["n_dofs"]) == UNKNOWNS_RECORD
        and float(d2["freq_hz"]) == FREQ_HZ
        and int(d2["ranks"]) == 16
        and s2.shape == (N_PORTS, N_PORTS)
        and o2[0] == "P17"
        and bij
        and set(o2) == set(a2["order"])
    )
    check(
        "(F) d2 copper",
        f_ok,
        f"{Path(ENTRY36_SAVE).name} ({Path(ENTRY36_SAVE).stat().st_size} B): degree "
        f"{int(d2['degree'])}, kappa_applied {bool(d2['kappa_applied'])}, cells {int(d2['cells'])} "
        f"(record {HOLE_CELL_RECORD}), n_dofs {int(d2['n_dofs'])} (record {UNKNOWNS_RECORD}), "
        f"f {float(d2['freq_hz']):.6e} Hz, ranks {int(d2['ranks'])}, shape {s2.shape}, "
        f"order[0] {o2[0]}, names bijection {bij}, set(order) == a2 {set(o2) == set(a2['order'])}",
    )

    # ---- (F1) d1 copper files ----
    f1 = []
    for tag, z in (("h12c", d1), ("h12c-a2", a2)):
        s = np.asarray(z["s32"])
        mism = [
            (p, m, n)
            for p, m, n in zip(
                z["order"], ans8._slot_port_map(z), ring_port_map_from_names(z["order"])
            )
            if m != n
        ]
        cells = int(z["cells"]) if "cells" in z else None
        ok = (
            s.shape == (N_PORTS, N_PORTS)
            and not mism
            and (cells is None or cells == HOLE_CELL_RECORD)
        )
        f1.append(ok)
        say(
            f"(F1) {tag}: shape {s.shape}, slot map vs names {len(mism)} mismatches {mism[:3]}, "
            f"cells {cells if cells is not None else 'key absent (printed)'} (record {HOLE_CELL_RECORD})"
        )
    check("(F1) d1 copper", all(f1), f"both h12c files {f1}")
    if fails:
        pytest.fail(f"(F)/(F1) red {fails}: the file is not what entry 36 wrote — stop")

    # ---- (Syn) names synthesis = slot synthesis on h12c-cu-s32.npz ----
    s1 = np.array(d1["s32"], dtype=np.complex128)
    _c, sp_names1, _r = _class_means(s1, *_names_az(d1["order"]))
    _c, sp_slot1, _r = _class_means(s1, *_slot_az(d1))
    check(
        "(Syn) names = slot synthesis",
        len(sp_names1) == 18 and sp_names1 == sp_slot1,
        f"{len(sp_names1)} / {len(sp_slot1)} classes, bitwise {sp_names1 == sp_slot1}",
    )
    if fails:
        pytest.fail(f"(Syn) red {fails}: the names synthesis is wrong — no d2 class mean read")

    def class_table(tag, spreads, rels, record):
        for k in record:
            say(
                f"   {tag} {'same ' if k[0] else 'other'} ring, {k[1]} steps  n = {spreads[k][1]:3d}  "
                f"mean |S| {spreads[k][2]:.10e} vs record {record[k]:.10e}, rel {rels[k]:.3e}"
            )

    # ---- (L-d1) degree-1 reproduction ----
    rels1 = _rels(sp_names1, CLASS_MEAN_CU_D1_RECORD)
    check(
        "(L-d1) d1 class means",
        max(rels1.values()) <= REUSE_REPRODUCTION_RTOL,
        f"worst rel {max(rels1.values()):.3e} (band {REUSE_REPRODUCTION_RTOL:g}; predicted <= 5e-10)",
    )
    class_table("Cu d1", sp_names1, rels1, CLASS_MEAN_CU_D1_RECORD)
    s_a2 = np.array(a2["s32"], dtype=np.complex128)
    pm_a2 = ring_port_map_from_names(a2["order"])
    v_a2 = ans8.rows_v2(s_a2, pm_a2, a2["order"], FREQ_HZ)
    c1 = v_a2["c_own_f"]
    r_key = abs(c1 / float(a2["c_mode"]) - 1.0)
    r_rec = abs(c1 / C_MODE_CU_D1_RECORD_F - 1.0)
    check(
        "(L-d1) C_mode",
        max(r_key, r_rec) <= REUSE_REPRODUCTION_RTOL and abs(c1 * 1e12 - C_MODE_CU_PF) <= 5e-5,
        f"c_own_f {c1:.12e} vs npz c_mode rel {r_key:.3e}, vs record rel {r_rec:.3e} (band "
        f"{REUSE_REPRODUCTION_RTOL:g}; :6367 read 1.776e-15); {c1 * 1e12:.6f} pF vs C_MODE_CU_PF "
        f"{C_MODE_CU_PF} (|d| <= 5e-5 pF)",
    )
    ob1 = v_a2["fixed_c"]["off_block"]
    check(
        "(L-d1) off-block",
        abs(ob1 / OFF_BLOCK_CU_D1_RECORD - 1.0) <= ans8.OFF_BLOCK_RTOL,
        f"{ob1:.6e} vs {OFF_BLOCK_CU_D1_RECORD:.6e} (band {ans8.OFF_BLOCK_RTOL:g})",
    )
    sm1 = _summary(v_a2)
    p17, cw = sm1["drives"]["p17_alone"], sm1["drives"]["cw"]
    got_p17 = (p17["feed_share"], p17["feed_leak"], p17["other_share"])
    got_cw = (cw["feed_share"], cw["other_share"], cw["feed_leak"], cw["other_leak"])
    d_sh = max(
        abs(a - b) for a, b in zip(got_p17 + got_cw, P17_2B_CU_D1_RECORD + CW_2B_CU_D1_RECORD)
    )
    check(
        "(L-d1) 2b shares/leaks",
        d_sh <= ans8.SHARE_ATOL,
        f"p17_alone {tuple(round(x, 6) for x in got_p17)}, cw {tuple(round(x, 6) for x in got_cw)}; "
        f"worst |d| {d_sh:.3e} (band {ans8.SHARE_ATOL:g})",
    )
    if any(f.startswith("(L-d1)") for f in fails):
        pytest.fail(f"(L-d1) red {fails}: d1 records do not reproduce — park")

    # ---- (G1') and (L-d2) ----
    rec2 = _reciprocity_ratio(s2)
    sig2 = float(np.linalg.norm(s2, 2))
    check(
        "(G1') d2 reciprocity",
        rec2 <= RECIPROCITY_BAND,
        f"{rec2:.6e} (band {RECIPROCITY_BAND:.0e}; window {RECIPROCITY_CU_D2_RECORD:.6e})",
    )
    check(
        "(G1') d2 sigma_max",
        sig2 <= MATRIX_PASSIVITY_CEILING,
        f"{sig2:.12f} (ceiling {MATRIX_PASSIVITY_CEILING})",
    )
    r_sig = abs(sig2 / SIGMA_MAX_CU_D2_RECORD - 1.0)
    check(
        "(L-d2) sigma_max",
        r_sig <= REUSE_REPRODUCTION_RTOL,
        f"recomputed {sig2:.12f} vs record {SIGMA_MAX_CU_D2_RECORD:.12f}: rel {r_sig:.3e} "
        f"(band {REUSE_REPRODUCTION_RTOL:g})",
    )
    _c, sp2, res2 = _class_means(s2, *_names_az(o2))
    rels2 = _rels(sp2, CLASS_MEAN_CU_D2_RECORD)
    check(
        "(L-d2) class means",
        len(sp2) == 18 and max(rels2.values()) <= REUSE_REPRODUCTION_RTOL,
        f"{len(sp2)} classes, azimuth residual {res2:.3e} deg; worst rel {max(rels2.values()):.3e} "
        f"(band {REUSE_REPRODUCTION_RTOL:g}; predicted <= 1e-9)",
    )
    class_table("Cu d2", sp2, rels2, CLASS_MEAN_CU_D2_RECORD)
    i17 = o2.index("P17")
    say(
        f"S17,17 d2 {s2[i17, i17].real:.6f}{s2[i17, i17].imag:+.6f}j beside record "
        f"{S1717_CU_D2_RECORD.real:.6f}{S1717_CU_D2_RECORD.imag:+.6f}j (printed)"
    )

    # ---- (L-sigma) sigma = 800 reproduction ----
    wf7 = {}
    for deg, name in zip((1, 2), BOXED_NPZ_BASENAMES):
        p = Path(SAVE_ROOT) / "wf7" / name
        assert p.exists(), f"{p} absent"
        with np.load(p, allow_pickle=False) as z:
            sv = {k: z[k] for k in z.files}
        names = [str(x) for x in sv["order"]]
        pm = ans8._slot_port_map(sv)
        wf7[deg] = {
            "saved": sv,
            "s": np.array(sv["s32"], dtype=np.complex128),
            "order": names,
            "pm": pm,
            "v2": ans8.rows_v2(np.array(sv["s32"], dtype=np.complex128), pm, names, FREQ_HZ),
        }
        rr = abs(wf7[deg]["v2"]["c_own_f"] / C_MODE_S800_RECORD_F[deg] - 1.0)
        check(
            f"(L-sigma) sigma800 d{deg} C_mode",
            rr <= REUSE_REPRODUCTION_RTOL,
            f"{wf7[deg]['v2']['c_own_f']:.9e} vs {C_MODE_S800_RECORD_F[deg]:.9e}: rel {rr:.3e}",
        )

    # ---- negative controls (asserted against the floor) ----
    floor = SAVE_CONTROL_FACTOR * REUSE_REPRODUCTION_RTOL
    sv1 = wf7[1]
    ctl = sv1["s"].copy()
    ctl[:, sv1["order"].index("P17")] *= CONTROL_COLUMN_SCALE
    _c, spc, _r = _class_means(ctl, *_slot_az(sv1["saved"]))
    nc = max(_rels(spc, CLASS_MEAN_D1_RECORD_BOXED).values())
    check(
        "(NC-s800 d1) P17 x 1.01",
        nc >= floor,
        f"worst rel {nc:.3e} = {nc / REUSE_REPRODUCTION_RTOL:.1f}x band (bar {SAVE_CONTROL_FACTOR:.0f}x; "
        f"item 366 read 3.147e-04)",
    )
    for tag, s, order, rec in (
        ("Cu d2", s2, o2, CLASS_MEAN_CU_D2_RECORD),
        ("Cu d1", s1, d1["order"], CLASS_MEAN_CU_D1_RECORD),
    ):
        ctl = s.copy()
        ctl[:, order.index("P17")] *= CONTROL_COLUMN_SCALE
        _c, spc, _r = _class_means(ctl, *_names_az(order))
        nc = max(_rels(spc, rec).values())
        check(
            f"(NC-{tag}) P17 x 1.01",
            nc >= floor,
            f"worst rel {nc:.3e} = {nc / REUSE_REPRODUCTION_RTOL:.1f}x band (bar "
            f"{SAVE_CONTROL_FACTOR:.0f}x); predicted {PRED_D2_COLUMN_CONTROL}: "
            f"{PRED_D2_COLUMN_CONTROL[0] <= nc <= PRED_D2_COLUMN_CONTROL[1]} (printed)",
        )
    if fails:
        pytest.fail(f"anchors red {fails} — rows not read")

    # ---- rows on d2 copper; (B0), (Sigma k), (S) ----
    t_rows = time.perf_counter()
    v_d2 = ans8.rows_v2(s2, pm2, o2, FREQ_HZ)
    b0 = (
        v_d2["own_c_mode"]["c_f"] == v_d2["fixed_c"]["c_mode_f"]
        and v_d2["c_own_f"] == v_d2["fixed_c"]["c_mode_f"]
    )
    check(
        "(B0) own-C wiring",
        b0,
        f"own c_f {v_d2['own_c_mode']['c_f']!r} vs fixed c_mode_f {v_d2['fixed_c']['c_mode_f']!r} "
        f"vs c_own_f {v_d2['c_own_f']!r} (bitwise)",
    )
    sk_worst, sk_bit = 0.0, True
    blocks = list(v_d2["own_c_mode"]["drives"].values()) + [
        c["drives"][lbl] for lbl, c in v_d2["controls"]["own_c_mode"].items()
    ]
    for blk in blocks:
        for ring in ("bottom", "top"):
            rr = blk["rings"][ring]
            fr = rr["fractions"]
            sk_worst = max(sk_worst, abs(sum(fr) - 1.0))
            sk_bit &= rr["share_k1"] == fr[1] + fr[15] and rr["leak_k0"] == fr[0]
    check(
        "(Sigma k)",
        sk_worst <= 1e-12 and sk_bit,
        f"{len(blocks)} drives/controls x 2 rings: worst |sum - 1| {sk_worst:.3e} (<= 1e-12); "
        f"share/leak bitwise {sk_bit}",
    )
    feeds2 = [pm2[o2.index(n)] for n in ans8.SPEC_FEED_NAMES]
    c2 = v_d2["c_own_f"]
    w_re, w_im, w_c, per_k = _esr_shift(s2, pm2, c2, feeds2, ESR_OHM)
    check(
        "(S) ESR uniform shift, Cu d2",
        max(w_re, w_im, w_c) <= SHIFT_RTOL,
        f"max |dRe - r|/r {w_re:.3e}, |dIm|/|Im| {w_im:.3e}, |dC_mode|/C_mode {w_c:.3e} "
        f"(band {SHIFT_RTOL:.0e}; PORT-24 sigma800 2.889e-15 / 4.571e-16)",
    )
    plant = np.full(N_PORTS, ESR_OHM)
    plant[o2.index("P17")] = 0.4
    _wr, _wi, _wc, per_k_p = _esr_shift(s2, pm2, c2, feeds2, plant)
    miss = max(per_k_p)
    check(
        "(P) planted 0.4 ohm on P17",
        miss >= SAVE_CONTROL_FACTOR * SHIFT_RTOL,
        f"max_k |dRe - r|/r {miss:.3e} (floor {SAVE_CONTROL_FACTOR * SHIFT_RTOL:.0e}); predicted "
        f"≈ {PRED_PLANT_MISS:.0e} (printed); per k " + " ".join(f"{x:.2e}" for x in per_k_p),
    )
    if fails:
        pytest.fail(f"identity red {fails}")

    # ---- the record: four columns ----
    cols = {
        "Cu d2": _summary(v_d2),
        "Cu d1 (a2)": sm1,
        "s800 d2 k-off": _summary(wf7[2]["v2"]),
        "s800 d1 k-off": _summary(wf7[1]["v2"]),
    }
    keys = list(cols)
    say(
        f"rows_v2 + identities: {time.perf_counter() - t_rows:.3f} s; columns {keys}; every 2b "
        f"share/leak: {CIRCUIT_LABEL}"
    )

    def line(label, fn, fmt, cs=cols):
        say(f"{label:<38}" + " | ".join(f"{k}: {format(fn(cs[k]), fmt)}" for k in cs))

    line("row 1  C_mode (F)", lambda c: c["c_mode_f"], ".9e")
    line("row 1  C_k2 (F)", lambda c: c["c_k2_f"], ".9e")
    line("row 1  off_block", lambda c: c["off_block"], ".6e")
    line("row 1  anti(k = 1)", lambda c: c["anti_k1"], ".6f")
    line("row 4a f0 one-point vs SPEC C (Hz)", lambda c: c["f0_hz"], ".2f")
    line(
        "row 4a 64 MHz*sqrt(C_mode/C_in(d2))",
        lambda c: FREQ_HZ * np.sqrt(c["c_mode_f"] / (C_IN_F_HUMAN_D2_PF * 1e-12)),
        ".2f",
    )
    for drv in ans8.SPEC_DRIVES:
        line(f"row 2b {drv} feed-ring share", lambda c, d=drv: c["drives"][d]["feed_share"], ".6f")
        line(f"row 2b {drv} feed-ring leak k0", lambda c, d=drv: c["drives"][d]["feed_leak"], ".6f")
        line(
            f"row 2b {drv} other-ring share", lambda c, d=drv: c["drives"][d]["other_share"], ".6f"
        )
        line(
            f"row 2b {drv} other-ring leak k0", lambda c, d=drv: c["drives"][d]["other_leak"], ".6f"
        )
        for ring in ("feed", "other"):
            say(
                f"   2b {drv} {ring}-ring |k| folds 0..8: "
                + " | ".join(
                    f"{k}: " + " ".join(f"{x:.4f}" for x in cols[k]["drives"][drv][f"{ring}_folds"])
                    for k in keys
                )
            )
    for lbl in ans8.SPEC_CONTROLS:
        for ring in ("bottom", "top"):
            line(
                f"row 2c {lbl} {ring} share", lambda c, a=lbl, r=ring: c["controls"][a][r][0], ".6f"
            )
            line(
                f"row 2c {lbl} {ring} leak k0",
                lambda c, a=lbl, r=ring: c["controls"][a][r][1],
                ".3e",
            )

    # ESR rows (production convention), d2 copper.
    v_esr = ring_mode_rows_v2(
        s2,
        port_map=pm2,
        freq_hz=FREQ_HZ,
        c_spec_f=ans8.SPEC_C_F,
        feeds=feeds2,
        drives=ans8.SPEC_DRIVES,
        controls=ans8.spec_controls(pm2, o2),
        esr_ohm=ESR_OHM,
    )
    esr_cols = {
        "Cu d2 lossless": cols["Cu d2"],
        "Cu d2 ESR (production convention)": _summary(v_esr),
    }
    for drv in ans8.SPEC_DRIVES:
        line(
            f"row 2b {drv} feed share [ESR]",
            lambda c, d=drv: c["drives"][d]["feed_share"],
            ".6f",
            esr_cols,
        )
        line(
            f"row 2b {drv} feed leak [ESR]",
            lambda c, d=drv: c["drives"][d]["feed_leak"],
            ".6f",
            esr_cols,
        )
        line(
            f"row 2b {drv} other share [ESR]",
            lambda c, d=drv: c["drives"][d]["other_share"],
            ".6f",
            esr_cols,
        )
    line("row 1  C_mode (F) [ESR]", lambda c: c["c_mode_f"], ".12e", esr_cols)

    # Printed controls.
    cu2 = cols["Cu d2"]
    v_c2 = ans8.rows_v2(s2, pm2, o2, FREQ_HZ, c_spec_f=cu2["c_k2_f"])
    nc_share = v_c2["fixed_c"]["drives"]["p17_alone"]["rings"]["bottom"]["share_k1"]
    say(
        f"(printed) C2 control, Cu d2: P17-alone feed-ring share at C_k2 {cu2['c_k2_f']:.6e} F = "
        f"{nc_share:.6f} (< {ans8.NEGATIVE_CONTROL_CEILING}: {nc_share < ans8.NEGATIVE_CONTROL_CEILING}; "
        f"H12c d1 Thevenin route 0.105707, not the same comparison)"
    )
    ctl_leak = max(v for c in cu2["controls"].values() for _sh, v in c.values())
    say(
        f"(printed) 2c controls' k = 0 leak at own C, Cu d2: max {ctl_leak:.3e} "
        f"(<= {ans8.CONTROL_LEAK_BAR:.0e}: {ctl_leak <= ans8.CONTROL_LEAK_BAR})"
    )

    # Order moves d2 - d1, copper and sigma = 800.
    for tag, a, b in (
        ("Cu", "Cu d2", "Cu d1 (a2)"),
        ("s800 k-off", "s800 d2 k-off", "s800 d1 k-off"),
    ):
        ca, cb = cols[a], cols[b]
        say(
            f"order move {tag} (d2 - d1)/d1: C_mode {ca['c_mode_f'] / cb['c_mode_f'] - 1:+.4e}, C_k2 "
            f"{ca['c_k2_f'] / cb['c_k2_f'] - 1:+.4e}, f0 {ca['f0_hz'] / cb['f0_hz'] - 1:+.4e}; 2b ccw d feed "
            f"share {ca['drives']['ccw']['feed_share'] - cb['drives']['ccw']['feed_share']:+.6f}, d feed "
            f"leak {ca['drives']['ccw']['feed_leak'] - cb['drives']['ccw']['feed_leak']:+.6f}"
        )

    # The H25 lines.
    cpf = cu2["c_mode_f"] * 1e12
    off = (cpf - C_IN_F_HUMAN_D2_PF) / C_IN_F_HUMAN_D2_PF
    hit = PRED_CMODE_CU_D2_PF_WINDOW[0] <= cpf <= PRED_CMODE_CU_D2_PF_WINDOW[1]
    say(
        f"H25 C_mode(d2, Cu) = {cpf:.4f} pF beside C_MODE_CU_PF {C_MODE_CU_PF} (d1) and "
        f"C_IN_F_HUMAN_D2_PF {C_IN_F_HUMAN_D2_PF} (d2 ESR root); (C_mode - root)/root {off:+.4e}; "
        f"order move on copper {cpf / C_MODE_CU_PF - 1:+.4e}"
    )
    say(
        f"H25 (predicted, printed) window {PRED_CMODE_CU_D2_PF_WINDOW} pF: {'HIT' if hit else 'MISS'}; "
        f"point predictions {PRED_CMODE_POINT_OFFSET_PF} (d {cpf - PRED_CMODE_POINT_OFFSET_PF:+.4f}) / "
        f"{PRED_CMODE_POINT_S800_PF} (d {cpf - PRED_CMODE_POINT_S800_PF:+.4f}) pF; degree-1 offset "
        f"(118.1186 - {ROOT_D1_PF}) / {ROOT_D1_PF} = {(C_MODE_CU_PF - ROOT_D1_PF) / ROOT_D1_PF:+.4e}"
    )
    ccw = cu2["drives"]["ccw"]
    say(
        f"H25 (predicted, printed) anti {cu2['anti_k1']:.6f} >= {PRED_ANTI_FLOOR}: "
        f"{cu2['anti_k1'] >= PRED_ANTI_FLOOR}; 2b ccw feed share {ccw['feed_share']:.6f} >= "
        f"{PRED_CCW_SHARE_FLOOR}: {ccw['feed_share'] >= PRED_CCW_SHARE_FLOOR}; leak {ccw['feed_leak']:.6f} "
        f"in {PRED_CCW_LEAK}: {PRED_CCW_LEAK[0] <= ccw['feed_leak'] <= PRED_CCW_LEAK[1]} [{CIRCUIT_LABEL}]"
    )
    say(f"elapsed {time.perf_counter() - t0:.2f} s")
    assert not fails, fails
