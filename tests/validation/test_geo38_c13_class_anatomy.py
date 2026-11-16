"""GEO-38 C13 (§9 item 622, §10 chain C13): the anatomy of the one C32 × mirror
class that spreads 0.9262 % on F-clinical-1.5T's degree-1 copper 64×64
("other ring, 1 step", `20261024T070008Z_GEO-38-C3.log:14768`), read from C3's
saved matrix (``logs/geo38/c3-cu-s64-d1.npz``, gitignored scratch, read only).

Knob ``FEM_EM_GEO38_C13`` (C10's ``_knob`` pattern, imported): unset / empty ⇒
the gated test skips; ``1`` ⇒ it runs and an absent input fails; anything
else ⇒ ``ValueError``.

Asserted (the weekly's anchor, nothing else on the fixture):

* **(RR)** reciprocity (``_reciprocity_ratio``) and σ_max re-read from the npz
  equal C10's prints at ``.6e`` / ``.12f`` string equality;
* **(CT)** the 34-class table rebuilt from the npz alone (``azimuth_deg`` /
  ``z_m`` → C3's own ``_class_table`` / ``_class_spreads``) reproduces C3's
  printed table line for line — an identity on the saved bytes;
* **(CN)** negative control: one port's ``z_m`` sign flipped ⇒ at least one
  class line differs from C3's (the classifier is not blind to the geometry).

Printed only: per-class n / mean |S| / absolute spread / relative spread; the
"other ring, 1 step" class's 128 members sorted by |S|; the same statistics on
F-human's degree-1 copper 32×32 (``h12c-cu-s32-a2.npz``, classes from
``map_from_slot_keys``) beside `WF-12` step 2's 3.2748 %; the abs-spread ratio
(worst class / "same ring, 1 step") for both fixtures; the transpose control.

No solve, no band, no glyph.  Not claimed: the known-issues retirement
condition (the 10-26 weekly's), any degree-2 reading, any AED value.
"""

from __future__ import annotations

import time

import numpy as np
import pytest

from tests.validation.test_geo38_f_clinical_c_mode import (
    C3_KEYS,
    H12C_A2_NPZ,
    WORKSPACE,
    _knob,
    map_from_geometry,
    map_from_slot_keys,
)
from tests.validation.test_geo38_f_clinical_matrix import (
    LEG_COUNT_32,
    N_AZIMUTH_CLASSES_32,
    N_SYMMETRY_CLASSES_32,
    PREDICTED_SPREAD_F_HUMAN,
    SAVE_PATH,
    _azimuth_class_32,
    _class_spreads,
    _class_table,
    _expected_class_sizes,
)
from tests.validation.test_port_birdcage_ring_column import _reciprocity_ratio

KNOB_C13 = "FEM_EM_GEO38_C13"

# ---- records (not bands), each with its log line ------------------------------
# C10 W1 re-read of the same npz: 20261024T093056Z_GEO-38-C10-W1.log:88
# (= C3's own prints, 20261024T070008Z_GEO-38-C3.log:14748–14749).
RECIPROCITY_RECORD_STR = "1.887168e-14"
SIGMA_MAX_RECORD_STR = "0.999970196163"

# C3's printed 34-class table, verbatim payloads after the "[GEO-38 C3] " prefix:
# 20261024T070008Z_GEO-38-C3.log:14750–14783.  A record, not a band.
C3_TABLE_LINES = (
    "     same  ring,  0 steps  n =   64  mean |S| = 4.524475151e-01  spread   0.1275%",
    "     same  ring,  1 steps  n =  128  mean |S| = 3.564592759e-01  spread   0.1680%",
    "     same  ring,  2 steps  n =  128  mean |S| = 1.786953961e-01  spread   0.1903%",
    "     same  ring,  3 steps  n =  128  mean |S| = 9.072620988e-02  spread   0.2416%",
    "     same  ring,  4 steps  n =  128  mean |S| = 6.197950666e-02  spread   0.2331%",
    "     same  ring,  5 steps  n =  128  mean |S| = 5.054011744e-02  spread   0.2160%",
    "     same  ring,  6 steps  n =  128  mean |S| = 4.028888776e-02  spread   0.2326%",
    "     same  ring,  7 steps  n =  128  mean |S| = 3.090778913e-02  spread   0.2832%",
    "     same  ring,  8 steps  n =  128  mean |S| = 2.330721425e-02  spread   0.2540%",
    "     same  ring,  9 steps  n =  128  mean |S| = 1.743371621e-02  spread   0.2432%",
    "     same  ring, 10 steps  n =  128  mean |S| = 1.288740175e-02  spread   0.3169%",
    "     same  ring, 11 steps  n =  128  mean |S| = 9.384239818e-03  spread   0.2857%",
    "     same  ring, 12 steps  n =  128  mean |S| = 6.859199453e-03  spread   0.2491%",
    "     same  ring, 13 steps  n =  128  mean |S| = 5.384955757e-03  spread   0.1843%",
    "     same  ring, 14 steps  n =  128  mean |S| = 4.883488436e-03  spread   0.1264%",
    "     same  ring, 15 steps  n =  128  mean |S| = 4.904352113e-03  spread   0.0783%",
    "     same  ring, 16 steps  n =   64  mean |S| = 4.969564071e-03  spread   0.0640%",
    "     other ring,  0 steps  n =   64  mean |S| = 6.049332318e-01  spread   0.1249%",
    "     other ring,  1 steps  n =  128  mean |S| = 5.741663544e-02  spread   0.9262%   <-- worst",
    "     other ring,  2 steps  n =  128  mean |S| = 8.744101876e-02  spread   0.2129%",
    "     other ring,  3 steps  n =  128  mean |S| = 1.004498174e-01  spread   0.1491%",
    "     other ring,  4 steps  n =  128  mean |S| = 7.830309896e-02  spread   0.1868%",
    "     other ring,  5 steps  n =  128  mean |S| = 5.585925938e-02  spread   0.2426%",
    "     other ring,  6 steps  n =  128  mean |S| = 4.053226710e-02  spread   0.2610%",
    "     other ring,  7 steps  n =  128  mean |S| = 3.046261845e-02  spread   0.2523%",
    "     other ring,  8 steps  n =  128  mean |S| = 2.316699936e-02  spread   0.2554%",
    "     other ring,  9 steps  n =  128  mean |S| = 1.747134854e-02  spread   0.2364%",
    "     other ring, 10 steps  n =  128  mean |S| = 1.294368897e-02  spread   0.3217%",
    "     other ring, 11 steps  n =  128  mean |S| = 9.415944671e-03  spread   0.2439%",
    "     other ring, 12 steps  n =  128  mean |S| = 6.869743489e-03  spread   0.2400%",
    "     other ring, 13 steps  n =  128  mean |S| = 5.381813602e-03  spread   0.1776%",
    "     other ring, 14 steps  n =  128  mean |S| = 4.871808264e-03  spread   0.1119%",
    "     other ring, 15 steps  n =  128  mean |S| = 4.888582837e-03  spread   0.0933%",
    "     other ring, 16 steps  n =   64  mean |S| = 4.952693246e-03  spread   0.0615%",
)
WORST_MARKER = "<-- worst"

# `WF-12` step 2's F-human copper record (`docs/planning/chunks/WF-12.md:142`,
# `20261002T094751Z_WF-12-step2-W2.log:6391–6410`): worst class (other ring, 2 steps).
WF12_WORST_CLASS_RECORD = (False, 2)
WF12_WORST_SPREAD_RECORD = 0.032748


def c13_enabled():
    return _knob(KNOB_C13)


def _norm(line):
    """(CT) normalisation: drop the worst marker, collapse whitespace runs."""
    return " ".join(line.replace(WORST_MARKER, "").split())


def _table_lines(spreads, n_steps):
    """C3's line format (`test_geo38_f_clinical_matrix.py:455–456`), unmarked."""
    out = []
    for same_ring in (True, False):
        for steps in range(n_steps):
            sp, cnt, mean = spreads[(same_ring, steps)]
            out.append(f"     {'same ' if same_ring else 'other'} ring, {steps:2d} steps  n = {cnt:4d}"
                       f"  mean |S| = {mean:.9e}  spread {sp * 100:8.4f}%")
    return out


def _abs_spreads(s, classes):
    """{key: max|S| − min|S|} — new arithmetic (``_class_spreads`` returns rel, n, mean)."""
    out = {}
    for key, pairs in classes.items():
        m = np.array([abs(s[i, j]) for i, j in pairs])
        out[key] = float(m.max() - m.min())
    return out


def _slot_class(n_slots):
    def classify(k_i, k_j):
        d = abs(int(k_i) - int(k_j)) % n_slots
        return min(d, n_slots - d), 0.0
    return classify


def _say(msg):
    print(f"[GEO-38 C13] {msg}", flush=True)


def _print_table(tag, spreads, abs_sp, n_steps):
    for same_ring in (True, False):
        for steps in range(n_steps):
            sp, cnt, mean = spreads[(same_ring, steps)]
            _say(f"{tag} {'same ' if same_ring else 'other'} ring, {steps:2d} steps  n = {cnt:4d}  "
                 f"mean |S| {mean:.9e}  abs spread {abs_sp[(same_ring, steps)]:.6e}  "
                 f"rel {sp * 100:8.4f}%")


def _anatomy(tag, s, pairs, label):
    """Members of one class sorted by |S|; sub-means by the signed top − bottom slot offset."""
    rows = sorted(((abs(s[i, j]), i, j) for i, j in pairs))
    for m, i, j in rows:
        _say(f"{tag} member ({i:2d},{j:2d}) {label(i)} -> {label(j)}  |S| {m:.9e}  "
             f"arg S {np.degrees(np.angle(s[i, j])):+9.4f} deg")
    lo, hi = rows[0], rows[-1]
    _say(f"{tag} MIN ({lo[1]},{lo[2]}) {label(lo[1])} -> {label(lo[2])} |S| {lo[0]:.9e}; "
         f"MAX ({hi[1]},{hi[2]}) {label(hi[1])} -> {label(hi[2])} |S| {hi[0]:.9e}")
    return rows


# --------------------------------------------------------------------------- W0
def test_c13_knob(monkeypatch):
    """Ungated: the knob (unset/empty ⇒ off, '1' ⇒ on, else ValueError) and the
    (CT) normaliser on C3's own record lines."""
    monkeypatch.delenv(KNOB_C13, raising=False)
    assert c13_enabled() is False
    monkeypatch.setenv(KNOB_C13, "1")
    assert c13_enabled() is True
    for bad in ("0", "yes", "build"):
        monkeypatch.setenv(KNOB_C13, bad)
        with pytest.raises(ValueError):
            c13_enabled()
    monkeypatch.delenv(KNOB_C13, raising=False)
    assert len(C3_TABLE_LINES) == N_SYMMETRY_CLASSES_32
    assert sum(WORST_MARKER in ln for ln in C3_TABLE_LINES) == 1


# --------------------------------------------------------------------------- W1
def test_geo38_c13_class_anatomy():
    if not c13_enabled():
        pytest.skip(f"{KNOB_C13} unset (reads gitignored logs/geo38 npz)")
    t0 = time.perf_counter()
    path = WORKSPACE / SAVE_PATH
    assert path.exists(), f"{path} absent: entry 44's matrix; never re-solve"
    with np.load(path, allow_pickle=False) as z:
        sv = {k: z[k] for k in z.files}
    missing = [k for k in C3_KEYS if k not in sv]
    assert not missing, f"C3 npz lacks {missing}"
    s = np.array(sv["s64"], dtype=np.complex128)
    order = [str(x) for x in sv["order"]]
    az = np.asarray(sv["azimuth_deg"], dtype=float)
    zc = np.asarray(sv["z_m"], dtype=float)
    _say(f"{path} ({path.stat().st_size} B); s64 {s.shape} {s.dtype}; order {order[0]}…{order[-1]}")

    # ---- (RR) ----
    ratio = _reciprocity_ratio(s)
    sigma_max = float(np.linalg.norm(s, 2))
    rr = (f"{ratio:.6e}", f"{sigma_max:.12f}")
    _say(f"(RR) reciprocity {rr[0]} vs C10 {RECIPROCITY_RECORD_STR}; sigma_max {rr[1]} vs C10 "
         f"{SIGMA_MAX_RECORD_STR} (string equality, ASSERTED)")

    # ---- (CT) ----
    classes, res = _class_table(az, zc, classify=_azimuth_class_32)
    sizes = {k: len(v) for k, v in classes.items()}
    spreads = _class_spreads(s, classes)
    mine = _table_lines(spreads, N_AZIMUTH_CLASSES_32)
    ct_miss = [(a, b) for a, b in zip(mine, C3_TABLE_LINES) if _norm(a) != _norm(b)]
    _say(f"(CT) classes {len(classes)} (expected {N_SYMMETRY_CLASSES_32}), sizes equal expected "
         f"{sizes == _expected_class_sizes()}, worst pair residual {res:.3e} deg; lines matching "
         f"C3.log:14750–14783 {len(mine) - len(ct_miss)}/{len(C3_TABLE_LINES)} (normalisation: payload "
         f"after the '[GEO-38 C3] ' prefix, '{WORST_MARKER}' dropped, whitespace runs collapsed; ASSERTED)")
    for a, b in ct_miss:
        _say(f"(CT) MISS mine '{_norm(a)}' vs C3 '{_norm(b)}'")

    # ---- (CN) negative control: port 0's z sign flipped ----
    zf = zc.copy()
    zf[0] = -zf[0]
    cls_f, _r = _class_table(az, zf, classify=_azimuth_class_32)
    sp_f = _class_spreads(s, cls_f)
    n_changed = sum(1 for k in classes if len(cls_f.get(k, [])) != len(classes[k]))
    lines_f = _table_lines(sp_f, N_AZIMUTH_CLASSES_32) if len(cls_f) == len(classes) else []
    cn_diff = (len(lines_f) != len(C3_TABLE_LINES)
               or sum(_norm(a) != _norm(b) for a, b in zip(lines_f, C3_TABLE_LINES)))
    _say(f"(CN) {order[0]} z sign flipped: classes {len(cls_f)}, classes with changed n {n_changed} "
         f"(predicted >= 2), differing lines {cn_diff} (ASSERTED > 0)")

    # ---- transpose control (printed only) ----
    sp_t = _class_spreads(s.T, classes)
    t_match = sum(_norm(a) == _norm(b)
                  for a, b in zip(_table_lines(sp_t, N_AZIMUTH_CLASSES_32), C3_TABLE_LINES))
    _say(f"transpose control (printed): s64.T reproduces {t_match}/{len(C3_TABLE_LINES)} class lines "
         f"(predicted {len(C3_TABLE_LINES)})")

    assert rr == (RECIPROCITY_RECORD_STR, SIGMA_MAX_RECORD_STR), rr
    assert len(classes) == N_SYMMETRY_CLASSES_32 and sizes == _expected_class_sizes(), sizes
    assert not ct_miss, ct_miss[:3]
    assert cn_diff, "flip control reproduced C3's table: classifier blind to z"

    # ---- the record (printed only, after the anchor holds) ----
    abs_sp = _abs_spreads(s, classes)
    _print_table("F-clinical", spreads, abs_sp, N_AZIMUTH_CLASSES_32)
    worst = max(spreads, key=lambda k: spreads[k][0])
    geo = map_from_geometry(az, zc)

    def label64(i):
        return f"{order[i]}(slot {geo[i][0]:2d}, {geo[i][1]:6s})"

    _say(f"F-clinical worst class {worst}: rel {spreads[worst][0] * 100:.4f}%, abs {abs_sp[worst]:.6e}")
    rows = _anatomy("F-clinical other1", s, classes[(False, 1)], label64)
    by_off = {}
    for m, i, j in rows:
        bi, ti = (i, j) if geo[i][1] == "bottom" else (j, i)
        half = LEG_COUNT_32 // 2
        d = (geo[ti][0] - geo[bi][0] + half) % LEG_COUNT_32 - half
        by_off.setdefault(d, []).append(m)
    for d, ms in sorted(by_off.items()):
        ms = np.array(ms)
        _say(f"F-clinical other1 top slot − bottom slot = {d:+d}: n {len(ms)}, mean |S| {ms.mean():.9e}, "
             f"min {ms.min():.9e}, max {ms.max():.9e}, abs spread {ms.max() - ms.min():.6e}")
    ratio64 = abs_sp[worst] / abs_sp[(True, 1)]
    _say(f"F-clinical abs-spread ratio worst {worst} / (same ring, 1 step): {abs_sp[worst]:.6e} / "
         f"{abs_sp[(True, 1)]:.6e} = {ratio64:.4f}")

    # ---- F-human degree-1 copper 32×32 (printed only) ----
    if not H12C_A2_NPZ.exists():
        _say("F-human comparand not available (file absent)")
    else:
        with np.load(H12C_A2_NPZ, allow_pickle=False) as z:
            fh = {k: z[k] for k in z.files}
        if "slot_k" not in fh or "slot_sign" not in fh:
            _say("F-human comparand not available (no slot_k / slot_sign)")
        else:
            s32 = np.array(fh["s32"], dtype=np.complex128)
            names = [str(x) for x in fh["order"]]
            pm = map_from_slot_keys(fh["slot_k"], fh["slot_sign"])
            n_slots = len({k for k, _r in pm})
            k32 = np.array([k for k, _r in pm], dtype=float)
            z32 = np.array([1.0 if r == "top" else -1.0 for _k, r in pm])
            cls32, _r = _class_table(k32, z32, classify=_slot_class(n_slots))
            sp32 = _class_spreads(s32, cls32)
            ab32 = _abs_spreads(s32, cls32)
            n_steps32 = n_slots // 2 + 1
            _say(f"F-human {H12C_A2_NPZ.name}: s32 {s32.shape}; slots {n_slots}, map unique "
                 f"{len(set(pm))}/{len(pm)}, classes {len(cls32)}")
            _print_table("F-human", sp32, ab32, n_steps32)
            w32 = max(sp32, key=lambda k: sp32[k][0])
            _say(f"F-human worst class {w32}: rel {sp32[w32][0] * 100:.4f}% beside WF-12 step 2 "
                 f"{WF12_WORST_SPREAD_RECORD * 100:.4f}% on {WF12_WORST_CLASS_RECORD} (WF-12.md:142; "
                 f"PREDICTED_SPREAD_F_HUMAN {PREDICTED_SPREAD_F_HUMAN}); printed, not asserted "
                 f"(classification route differs)")

            def label32(i):
                return f"{names[i]}(slot {pm[i][0]:2d}, {pm[i][1]:6s})"

            lo = min(cls32[w32], key=lambda p: abs(s32[p]))
            hi = max(cls32[w32], key=lambda p: abs(s32[p]))
            _say(f"F-human worst-class MIN {lo} {label32(lo[0])} -> {label32(lo[1])} |S| "
                 f"{abs(s32[lo]):.9e}; MAX {hi} {label32(hi[0])} -> {label32(hi[1])} |S| {abs(s32[hi]):.9e}")
            if (False, 1) in sp32:
                _say(f"F-human (other ring, 1 step): rel {sp32[(False, 1)][0] * 100:.4f}%, abs "
                     f"{ab32[(False, 1)]:.6e}")
            ratio32 = ab32[w32] / ab32[(True, 1)]
            _say(f"F-human abs-spread ratio worst {w32} / (same ring, 1 step): {ab32[w32]:.6e} / "
                 f"{ab32[(True, 1)]:.6e} = {ratio32:.4f}")
            _say(f"SUMMARY abs-spread ratio (worst / same ring 1 step): F-clinical {ratio64:.4f} "
                 f"[{worst}], F-human {ratio32:.4f} [{w32}]")
    _say(f"body {time.perf_counter() - t0:.2f} s")
