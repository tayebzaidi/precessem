"""`ANS-8` step 1c (§9 item 366, §10 H15c): SPEC v2 rows on the saved degree-2 32x32.

No solve.  Two gitignored matrices saved by
`tests/validation/test_wf7_f_human_ring_matrix.py` (H15a / H15b, ``:917–931``)
— raw S at 50 ohm, **kappa off**, the boxed F-human sigma = 800 mesh, 64 MHz —
are run through the one shared ``ring_mode_rows_v2`` (via ``ans:8``'s
``rows_v2``), beside the tracked degree-1 **kappa-on** ``metrics.json``
``rows_v2`` (loaded, never retyped).  The order move is read at matched kappa
(degree 2 kappa off vs degree 1 kappa off); the kappa move at degree 1 is printed
beside it, so kappa and order never mix in one difference.

This module prints our numbers only: no band, tolerance or target is set; every
band is imported by name.  The two record constants below are comparands (the
361-3 precedent), each with its log line.

Knob ``FEM_EM_ANS8_STEP1C``: unset => the gated test skips; ``"1"`` => runs, and an
absent degree-2 npz **fails**; anything else => ``ValueError``.
"""

from __future__ import annotations

import importlib.util
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import pytest

KNOB = "FEM_EM_ANS8_STEP1C"
REPO = Path(__file__).resolve().parents[2]
CASE = REPO / "examples" / "ansys_benchmarks" / "ans8_f_human_quadrature_mode_content_64MHz"
SCRIPT = CASE / "08_f_human_quadrature_mode_content.py"
METRICS = CASE / "metrics.json"
FREQ_HZ = 64.0e6
N_PORTS = 32
N_PER_RING = 16

# Record comparands for (L-d2): the H15b window's own reading of the file it
# wrote (`20261009T070009Z_WF-7-step0e-H15b.log`), at its printed 10 digits.
SIGMA_MAX_D2_RECORD_BOXED = 0.999913975  # :12320, :12343
CLASS_MEAN_D2_RECORD_BOXED = {
    (True, 0): 5.573876225e-01,   # :12321
    (True, 1): 1.993723855e-01,   # :12322
    (True, 2): 1.022202605e-01,   # :12323
    (True, 3): 6.892532731e-02,   # :12324
    (True, 4): 5.998949132e-02,   # :12325
    (True, 5): 5.930025714e-02,   # :12326
    (True, 6): 6.052613677e-02,   # :12327
    (True, 7): 6.165646758e-02,   # :12328
    (True, 8): 6.208039048e-02,   # :12329
    (False, 0): 6.255404023e-01,  # :12330
    (False, 1): 6.686527944e-02,  # :12331
    (False, 2): 3.560407807e-02,  # :12332
    (False, 3): 5.477784226e-02,  # :12333
    (False, 4): 5.964939474e-02,  # :12334
    (False, 5): 6.121184642e-02,  # :12335
    (False, 6): 6.220602265e-02,  # :12336
    (False, 7): 6.293280415e-02,  # :12337
    (False, 8): 6.320841099e-02,  # :12338
}

# Predictions (written first, printed, never asserted — rule (e)).
PRED_CMODE_ORDER_MOVE = (0.005, 0.06)      # |d2 - d1| / d1, matched kappa
PRED_SHARE_ATOL = 0.03                     # 2b ccw feed-ring share vs degree 1 (kappa on)
PRED_LEAK_ATOL = 0.03                      # 2b ccw leak vs degree 1 (kappa on)
PRED_ANTI_FLOOR = 0.99
PRED_D2_CONTROL = (2.0e-4, 4.0e-4)         # P17 x 1.01 vs CLASS_MEAN_D2_RECORD_BOXED


def _enabled() -> bool:
    raw = os.environ.get(KNOB, "").strip()
    if not raw:
        return False
    if raw != "1":
        raise ValueError(f"{KNOB} must be unset or '1', got {raw!r}")
    return True


def _need():
    if not _enabled():
        pytest.skip(f"{KNOB} unset (H15c reads gitignored logs/wf7/*.npz)")


def _load_ans8():
    spec = importlib.util.spec_from_file_location("ans8_step1c", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["ans8_step1c"] = mod
    spec.loader.exec_module(mod)
    return mod


def _synth_az(saved):
    """``az = slot_k * AZIMUTH_STEP_DEG``, ``zc = slot_sign`` (the npz has no azimuths)."""
    from tests.validation.test_port_birdcage_ring_matrix import AZIMUTH_STEP_DEG
    az = np.asarray(saved["slot_k"], dtype=float) * AZIMUTH_STEP_DEG
    zc = np.asarray(saved["slot_sign"], dtype=float)
    return az, zc


def _class_mean_rels(s, saved, record):
    from tests.validation.test_wf7_f_human_ring_matrix import _class_means
    az, zc = _synth_az(saved)
    classes, spreads, resid = _class_means(s, az, zc)
    rels = {k: abs(spreads[k][2] / record[k] - 1.0) for k in record}
    return classes, spreads, resid, rels


# --------------------------------------------------------------------------- W0
def test_ans8_step1c_coding(monkeypatch):
    """Ungated: the knob parser, and ``_class_means`` on a synthetic D16 matrix."""
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

    # Synthetic D16-symmetric 32x32: |S_ij| depends only on (same ring, folded steps).
    slot_k = np.array(list(range(N_PER_RING)) * 2)
    slot_sign = np.array([-1.0] * N_PER_RING + [1.0] * N_PER_RING)
    rng = np.random.default_rng(366)
    val = {(same, d): complex(rng.uniform(0.01, 1.0), rng.uniform(-1.0, 1.0))
           for same in (True, False) for d in range(N_PER_RING // 2 + 1)}
    s = np.zeros((N_PORTS, N_PORTS), dtype=np.complex128)
    for i in range(N_PORTS):
        for j in range(N_PORTS):
            d = abs(int(slot_k[i]) - int(slot_k[j]))
            d = min(d, N_PER_RING - d)
            s[i, j] = val[(bool(slot_sign[i] * slot_sign[j] > 0), d)]
    record = {k: abs(v) for k, v in val.items()}
    classes, spreads, resid, rels = _class_mean_rels(s, {"slot_k": slot_k, "slot_sign": slot_sign},
                                                     record)
    print(f"[ANS-8 1c W0] synthetic D16: {len(classes)} classes, max spread "
          f"{max(v[0] for v in spreads.values()):.3e}, worst azimuth residual {resid:.3e} deg, "
          f"max class-mean rel {max(rels.values()):.3e}", flush=True)
    assert len(classes) == 18
    assert all(v[0] == 0.0 for v in spreads.values())
    assert max(rels.values()) <= 1e-12
    assert resid <= 1e-9


# --------------------------------------------------------------------------- W1
def _summary(v2):
    """The printed rows of one ``rows_v2`` dict (ours or the tracked JSON)."""
    own, fixed = v2["own_c_mode"], v2["fixed_c"]
    out = {"c_mode_f": float(v2["c_own_f"]), "c_k2_f": float(fixed["c_k2_f"]),
           "off_block": float(fixed["off_block"]), "anti_k1": float(v2["anti_k1"]),
           "f0_hz": float(v2["f0_one_point_hz"]), "drives": {}, "controls": {}}
    for drv, blk in own["drives"].items():
        rb, rt = blk["rings"]["bottom"], blk["rings"]["top"]
        out["drives"][drv] = {"feed_share": rb["share_k1"], "feed_leak": rb["leak_k0"],
                              "other_share": rt["share_k1"], "other_leak": rt["leak_k0"],
                              "feed_folds": rb["abs_k_folds"], "other_folds": rt["abs_k_folds"]}
    for lbl, ctl in v2["controls"]["own_c_mode"].items():
        rings = ctl["drives"][lbl]["rings"]
        out["controls"][lbl] = {r: (rings[r]["share_k1"], rings[r]["leak_k0"])
                                for r in ("bottom", "top")}
    return out


def test_ans8_step1c_degree2_rows():
    _need()
    t0 = time.perf_counter()
    from fem_em_solver.ports.ring_modes import ring_port_map_from_names
    from tests.validation.test_port_birdcage_ring_column import _reciprocity_ratio
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
        SIGMA_MAX_D1_RECORD_BOXED,
        WF7_BOXED_CELL_RECORD,
    )
    from tests.validation.test_port_birdcage_ring_column import CONTROL_COLUMN_SCALE

    ans8 = _load_ans8()
    fails = []

    def say(msg):
        print(f"[ANS-8 1c] {msg}", flush=True)

    def check(label, ok, msg):
        say(f"{label} {msg} — {'green' if ok else 'RED'}")
        if not ok:
            fails.append(label)

    paths = {deg: Path(SAVE_ROOT) / "wf7" / name
             for deg, name in zip((1, 2), BOXED_NPZ_BASENAMES)}
    assert paths[2].exists(), f"{paths[2]} absent: entry 29 is re-posed by the weekly; never re-solve"
    have_d1 = paths[1].exists()
    if not have_d1:
        say(f"{paths[1]} absent — degree 2 alone (printed); (L-d1) not run")

    data = {}
    for deg in ((1, 2) if have_d1 else (2,)):
        with np.load(paths[deg], allow_pickle=False) as z:
            saved = {k: z[k] for k in z.files}
        s = np.array(saved["s32"], dtype=np.complex128)
        order = [str(p) for p in saved["order"]]
        # (F) file identity.
        slots_ok = all(sorted(saved["slot_k"][saved["slot_sign"] == sg].tolist())
                       == list(range(N_PER_RING)) for sg in (1.0, -1.0))
        f_ok = (int(saved["degree"]) == deg and bool(saved["kappa_applied"]) is False
                and int(saved["cells"]) == WF7_BOXED_CELL_RECORD
                and float(saved["freq_hz"]) == FREQ_HZ and s.shape == (N_PORTS, N_PORTS)
                and len(order) == N_PORTS and slots_ok)
        check(f"(F) degree {deg}", f_ok,
              f"{paths[deg].name} ({paths[deg].stat().st_size} B): degree {int(saved['degree'])}, "
              f"kappa_applied {bool(saved['kappa_applied'])}, cells {int(saved['cells'])} "
              f"(record {WF7_BOXED_CELL_RECORD}), f {float(saved['freq_hz']):.6e} Hz, shape {s.shape}, "
              f"16 slots per ring {slots_ok}, ranks {int(saved['ranks'])}, kappa key "
              f"{float(saved['kappa']):.6e} (printed)")
        # (G2') slot map.
        measured = ans8._slot_port_map(saved)
        named = ring_port_map_from_names(order)
        mism = [(p, m, n) for p, m, n in zip(order, measured, named) if m != n]
        check(f"(G2') degree {deg} slot map", not mism,
              f"_slot_port_map vs ring_port_map_from_names(order): {len(mism)} mismatches {mism[:3]}")
        data[deg] = {"saved": saved, "s": s, "order": order, "port_map": named}
    if fails:
        pytest.fail(f"(F)/(G2') red {fails}: the file is not what entry 29 wrote — stop")

    # (G1') degree 2.
    s2 = data[2]["s"]
    rec2 = _reciprocity_ratio(s2)
    sig2 = float(np.linalg.norm(s2, 2))
    check("(G1') degree 2 reciprocity", rec2 <= RECIPROCITY_BAND,
          f"{rec2:.6e} (band {RECIPROCITY_BAND:.0e}; window 9.581e-14)")
    check("(G1') degree 2 sigma_max", sig2 <= MATRIX_PASSIVITY_CEILING,
          f"{sig2:.9f} (ceiling {MATRIX_PASSIVITY_CEILING}); margin {MATRIX_PASSIVITY_CEILING - sig2:+.3e}")

    def class_table(tag, spreads, rels, record):
        for k in record:
            say(f"   {tag} {'same ' if k[0] else 'other'} ring, {k[1]} steps  n = {spreads[k][1]:3d}  "
                f"mean |S| {spreads[k][2]:.10e} vs record {record[k]:.10e}, rel {rels[k]:.3e}")

    # (L-d1) degree-1 reproduction — the az-synthesis check.
    if have_d1:
        s1 = data[1]["s"]
        sig1 = float(np.linalg.norm(s1, 2))
        sig1_rel = abs(sig1 / SIGMA_MAX_D1_RECORD_BOXED - 1.0)
        classes1, spreads1, resid1, rels1 = _class_mean_rels(s1, data[1]["saved"], CLASS_MEAN_D1_RECORD_BOXED)
        w1 = max(rels1.values())
        check("(L-d1) sigma_max", sig1_rel <= REUSE_REPRODUCTION_RTOL,
              f"{sig1:.9f} vs {SIGMA_MAX_D1_RECORD_BOXED:.9f}: rel {sig1_rel:.3e} "
              f"(band {REUSE_REPRODUCTION_RTOL:g}; H15a read 2.596e-10)")
        check("(L-d1) class means", len(classes1) == 18 and w1 <= REUSE_REPRODUCTION_RTOL,
              f"{len(classes1)} classes, worst azimuth residual {resid1:.3e} deg; worst rel {w1:.3e} "
              f"(band {REUSE_REPRODUCTION_RTOL:g}; H15a read 4.279e-10)")
        class_table("d1", spreads1, rels1, CLASS_MEAN_D1_RECORD_BOXED)
        # Asserted negative control (backed): P17 x 1.01 on the degree-1 array.
        ctl1 = s1.copy()
        ctl1[:, data[1]["order"].index("P17")] *= CONTROL_COLUMN_SCALE
        _c, _s, _r, rels1c = _class_mean_rels(ctl1, data[1]["saved"], CLASS_MEAN_D1_RECORD_BOXED)
        nc1 = max(rels1c.values())
        check("(NC-d1) P17 x 1.01", nc1 >= SAVE_CONTROL_FACTOR * REUSE_REPRODUCTION_RTOL,
              f"worst class-mean rel {nc1:.3e} = {nc1 / REUSE_REPRODUCTION_RTOL:.1f}x the band "
              f"(bar {SAVE_CONTROL_FACTOR:.0f}x; H15a read 3.147e-04)")
        if any(f.startswith("(L-d1)") for f in fails):
            pytest.fail(f"(L-d1) red {fails}: the az synthesis is wrong — degree 2 not read")

    # (L-d2) the H15c anchor.
    sig2_key = float(data[2]["saved"]["sigma_max"])
    r_key = abs(sig2 / sig2_key - 1.0)
    r_rec = abs(sig2 / SIGMA_MAX_D2_RECORD_BOXED - 1.0)
    check("(L-d2) sigma_max", r_key <= REUSE_REPRODUCTION_RTOL and r_rec <= REUSE_REPRODUCTION_RTOL,
          f"recomputed {sig2:.12f} vs npz key {sig2_key:.12f} rel {r_key:.3e}; vs record "
          f"{SIGMA_MAX_D2_RECORD_BOXED:.9f} rel {r_rec:.3e} (band {REUSE_REPRODUCTION_RTOL:g})")
    classes2, spreads2, resid2, rels2 = _class_mean_rels(s2, data[2]["saved"], CLASS_MEAN_D2_RECORD_BOXED)
    w2 = max(rels2.values())
    check("(L-d2) class means", len(classes2) == 18 and w2 <= REUSE_REPRODUCTION_RTOL,
          f"{len(classes2)} classes, worst azimuth residual {resid2:.3e} deg; worst rel {w2:.3e} "
          f"(band {REUSE_REPRODUCTION_RTOL:g}; predicted <= 1e-9)")
    class_table("d2", spreads2, rels2, CLASS_MEAN_D2_RECORD_BOXED)
    ctl2 = s2.copy()
    ctl2[:, data[2]["order"].index("P17")] *= CONTROL_COLUMN_SCALE
    _c, _s, _r, rels2c = _class_mean_rels(ctl2, data[2]["saved"], CLASS_MEAN_D2_RECORD_BOXED)
    nc2 = max(rels2c.values())
    check("(NC-d2) P17 x 1.01", nc2 >= SAVE_CONTROL_FACTOR * REUSE_REPRODUCTION_RTOL,
          f"worst class-mean rel {nc2:.3e} = {nc2 / REUSE_REPRODUCTION_RTOL:.1f}x the band "
          f"(bar {SAVE_CONTROL_FACTOR:.0f}x); predicted {PRED_D2_CONTROL}: "
          f"{PRED_D2_CONTROL[0] <= nc2 <= PRED_D2_CONTROL[1]} (printed)")
    if fails:
        pytest.fail(f"anchors red {fails} — rows not read")

    # ---- rows (records, printed) ----
    t_rows = time.perf_counter()
    v2 = {deg: ans8.rows_v2(d["s"], d["port_map"], d["order"], FREQ_HZ) for deg, d in data.items()}
    tracked = json.loads(METRICS.read_text())["rows_v2"]
    cols = {"d2 k-off": _summary(v2[2])}
    if have_d1:
        cols["d1 k-off"] = _summary(v2[1])
    cols["d1 k-on"] = _summary(tracked)
    keys = list(cols)
    say(f"rows_v2 on {len(v2)} matrices: {time.perf_counter() - t_rows:.3f} s; columns {keys} "
        f"(d1 k-on = tracked metrics.json rows_v2, loaded)")

    def line(label, fn, fmt):
        say(f"{label:<38}" + " | ".join(f"{k}: {format(fn(cols[k]), fmt)}" for k in keys))

    line("row 1  C_mode (F)", lambda c: c["c_mode_f"], ".9e")
    line("row 1  C_k2 (F)", lambda c: c["c_k2_f"], ".9e")
    line("row 1  off_block", lambda c: c["off_block"], ".6e")
    line("row 1  anti(k = 1)", lambda c: c["anti_k1"], ".6f")
    line("row 4a f0 one-point vs SPEC C (Hz)", lambda c: c["f0_hz"], ".2f")
    for drv in ans8.SPEC_DRIVES:
        line(f"row 2b {drv} feed-ring share", lambda c, d=drv: c["drives"][d]["feed_share"], ".6f")
        line(f"row 2b {drv} feed-ring leak k0", lambda c, d=drv: c["drives"][d]["feed_leak"], ".6f")
        line(f"row 2b {drv} other-ring share", lambda c, d=drv: c["drives"][d]["other_share"], ".6f")
        for ring in ("feed", "other"):
            say(f"   2b {drv} {ring}-ring |k| folds 0..8: " + " | ".join(
                f"{k}: " + " ".join(f"{x:.4f}" for x in cols[k]["drives"][drv][f"{ring}_folds"])
                for k in keys))
    for lbl in ans8.SPEC_CONTROLS:
        for ring in ("bottom", "top"):
            line(f"row 2c {lbl} {ring} share", lambda c, a=lbl, r=ring: c["controls"][a][r][0], ".6f")
            line(f"row 2c {lbl} {ring} leak k0", lambda c, a=lbl, r=ring: c["controls"][a][r][1], ".3e")

    # The C2 control (printed only on degree 2): P17-alone feed-ring share at C_k2.
    d2 = cols["d2 k-off"]
    v2_c2 = ans8.rows_v2(s2, data[2]["port_map"], data[2]["order"], FREQ_HZ, c_spec_f=d2["c_k2_f"])
    nc_share = v2_c2["fixed_c"]["drives"]["p17_alone"]["rings"]["bottom"]["share_k1"]
    say(f"(printed) C2 control, degree 2: P17-alone feed-ring share at C_k2 {d2['c_k2_f']:.6e} F = "
        f"{nc_share:.6f} (< {ans8.NEGATIVE_CONTROL_CEILING}: {nc_share < ans8.NEGATIVE_CONTROL_CEILING}; "
        f"degree 1 k-on 0.099704)")
    ctl_leak = max(v for c in d2["controls"].values() for _sh, v in c.values())
    say(f"(printed) 2c controls' k = 0 leak at own C, degree 2: max {ctl_leak:.3e} "
        f"(<= {ans8.CONTROL_LEAK_BAR:.0e}: {ctl_leak <= ans8.CONTROL_LEAK_BAR})")

    # The order move at matched kappa, and the kappa move at degree 1.
    def moves(a, b, tag):
        rel = {q: a[q] / b[q] - 1.0 for q in ("c_mode_f", "c_k2_f", "f0_hz")}
        say(f"{tag}: (a - b)/b C_mode {rel['c_mode_f']:+.4e}, C_k2 {rel['c_k2_f']:+.4e}, "
            f"f0 {rel['f0_hz']:+.4e}")
        for drv in ans8.SPEC_DRIVES:
            da, db = a["drives"][drv], b["drives"][drv]
            say(f"   {tag} 2b {drv}: d feed share {da['feed_share'] - db['feed_share']:+.6f}, "
                f"d feed leak {da['feed_leak'] - db['feed_leak']:+.6f}, "
                f"d other share {da['other_share'] - db['other_share']:+.6f}")
        return rel

    if have_d1:
        om = moves(cols["d2 k-off"], cols["d1 k-off"], "order move (d2 - d1, kappa off)")
        moves(cols["d1 k-on"], cols["d1 k-off"], "kappa move (d1 k-on - d1 k-off)")
        in_c = PRED_CMODE_ORDER_MOVE[0] <= abs(om["c_mode_f"]) <= PRED_CMODE_ORDER_MOVE[1]
        say(f"(predicted, printed) |C_mode order move| {abs(om['c_mode_f']):.4e} in "
            f"{PRED_CMODE_ORDER_MOVE}: {in_c}; f0 move / C_mode move "
            f"{om['f0_hz'] / om['c_mode_f']:+.4f} (predicted ≈ +0.5)")
    ref = cols["d1 k-on"]["drives"]["ccw"]
    dsh = d2["drives"]["ccw"]["feed_share"] - ref["feed_share"]
    dlk = d2["drives"]["ccw"]["feed_leak"] - ref["feed_leak"]
    say(f"(predicted, printed) 2b ccw feed-ring share d2 {d2['drives']['ccw']['feed_share']:.6f} vs "
        f"d1 k-on {ref['feed_share']:.6f}: d {dsh:+.6f} within {PRED_SHARE_ATOL}: {abs(dsh) <= PRED_SHARE_ATOL}; "
        f"leak d2 {d2['drives']['ccw']['feed_leak']:.6f} vs {ref['feed_leak']:.6f}: d {dlk:+.6f} within "
        f"{PRED_LEAK_ATOL}: {abs(dlk) <= PRED_LEAK_ATOL}; anti d2 {d2['anti_k1']:.6f} >= {PRED_ANTI_FLOOR}: "
        f"{d2['anti_k1'] >= PRED_ANTI_FLOOR}")
    say(f"elapsed {time.perf_counter() - t0:.2f} s")
    assert not fails, fails
