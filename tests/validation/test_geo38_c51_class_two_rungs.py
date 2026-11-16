"""GEO-38 C51 (§9 item 1003, §10 chain C51): F-clinical-1.5T's degree-1 copper
64×64 class table on two h-rungs — r1 (entry 44, ``logs/geo38/c3-cu-s64-d1.npz``,
`20261024T070008Z_GEO-38-C3.log`) and r2 (entry 65, the coarser rung,
``logs/geo38/c3-cu-s64-d1-r2.npz``, `20261116T080008Z_GEO-38-C35w.log`).

Both logs are parsed (``read_log_text(..., errors="strict")``, plain or ``.gz``)
by :func:`parse_c3_log`: every key needs **exactly one distinct match**, count
printed; any other count raises ``ValueError`` naming it.  Both npz files are
gitignored scratch, read only (``SAVE_PATH_R2`` / ``ENTRY44_D1_PATH`` literals —
never C3's rebindable ``SAVE_PATH``).

Knob ``FEM_EM_GEO38_C51`` (C10's ``_knob``): unset / empty ⇒ the gated test
skips; ``1`` ⇒ it runs and an absent input fails; anything else ⇒ ``ValueError``.

Asserted (the anchor): (P1) entry 44's parse reproduces C13's landed record;
(RR), (CT) per rung, npz against its own log by value; (Cus) the r2 npz's
custody; (C13r) the r1 C13 abs-spread statistic; (CN) the r1 npz read as r2
fails (CT) and (RR) against entry 65.  Printed only: the per-class r2/r1
relative-spread record (r2 is the COARSER rung: "falls ≥ 2× with finer h"
means r2/r1 ≥ 2), absolute spreads and the C13 statistic on both rungs.  No
solve, band, verdict or branch.
"""

from __future__ import annotations

import re
import time

import numpy as np
import pytest

from tests.log_records import read_log_text
from tests.validation.test_geo38_c13_class_anatomy import (
    C3_TABLE_LINES,
    RECIPROCITY_RECORD_STR,
    SIGMA_MAX_RECORD_STR,
    _abs_spreads,
    _norm,
    _table_lines,
)
from tests.validation.test_geo38_f_clinical_build import F_CLINICAL_R2_CELL_RECORD
from tests.validation.test_geo38_f_clinical_c_mode import WORKSPACE, _knob
from tests.validation.test_geo38_f_clinical_matrix import (
    ENTRY44_D1_PATH,
    N_AZIMUTH_CLASSES_32,
    N_SYMMETRY_CLASSES_32,
    SAVE_PATH_R2,
    _azimuth_class_32,
    _class_spreads,
    _class_table,
    _expected_class_sizes,
)
from tests.validation.test_port_birdcage_ring_column import _reciprocity_ratio

KNOB_C51 = "FEM_EM_GEO38_C51"
LOG_DIR = WORKSPACE / "docs" / "testing" / "logs"
ENTRY44_LOG = LOG_DIR / "20261024T070008Z_GEO-38-C3.log"
ENTRY65_LOG = LOG_DIR / "20261116T080008Z_GEO-38-C35w.log"
C13_W1_LOG = LOG_DIR / "20261025T093429Z_GEO-38-C13-W1.log"

# ---- records (not bands), each with its log line ------------------------------
# C13 W1: 20261025T093429Z_GEO-38-C13-W1.log:237 (abs-spread ratio worst / same ring 1 step).
C13_STAT_R1_RECORD_STR = "0.8883"
EXPECTED_N_RANKS = 16  # entry 65's window (`-n 16`, log :12508)
EXPECTED_ORDER = [f"P{i}" for i in range(33, 97)]
# The drafter's arithmetic at the 03:00 11-16 review (printed beside, never asserted).
PRED_ASYM_RATIO = 1.127
PRED_MEDIAN = 1.501
PRED_COUNT_GE_2X = 3

_P = r"\[GEO-38 C3\] "
_RE_CLASS = re.compile(_P + r"(?P<payload>\s+(?P<ring>same |other) ring,\s*(?P<steps>\d+) steps\s+n =\s*"
                       r"(?P<n>\d+)\s+mean \|S\| = (?P<mean>\S+)\s+spread\s+(?P<spread>\S+)%[^\n]*)")
_RE_RECIP = re.compile(_P + r"\(ii\) reciprocity (\S+) ")
_RE_SIGMA = re.compile(_P + r"\(iii\) sigma_max (\S+) ")
_RE_A0 = re.compile(_P + r"fixture: [^\n]*?\(A0\) cells (\d+) vs record (\d+)")
_RE_SAVED = re.compile(_P + r"saved (\S+); window")
_RE_C13 = re.compile(r"\[GEO-38 C13\] F-clinical abs-spread ratio worst [^\n]*= (\S+)[ \t]*$", re.M)


def c51_enabled():
    return _knob(KNOB_C51)


def _say(msg):
    print(f"[GEO-38 C51] {msg}", flush=True)


def _one(name, values, tag):
    """Exactly one distinct value, total and distinct counts printed; else ValueError."""
    distinct = sorted(set(values))
    _say(f"{tag} parse key {name}: {len(values)} match(es), {len(distinct)} distinct")
    if len(distinct) != 1:
        raise ValueError(f"{tag} key {name}: {len(distinct)} distinct matches (count {len(distinct)}), "
                         f"expected exactly 1")
    return distinct[0]


def parse_c3_log(text, tag="log"):
    """C3's printed (ii), (iii), (A0), 34 class lines and saved path.

    Regexes anchored on the ``[GEO-38 C3] `` prefix, not on line start.  Returns
    ``{"recip", "sigma", "cells", "cells_record", "saved", "table"}`` where
    ``table[(same_ring, steps)] = (payload, n, mean_str, spread_str)``."""
    out = {
        "recip": _one("(ii) reciprocity", _RE_RECIP.findall(text), tag),
        "sigma": _one("(iii) sigma_max", _RE_SIGMA.findall(text), tag),
        "saved": _one("saved", _RE_SAVED.findall(text), tag),
    }
    cells, cells_rec = _one("(A0) cells", _RE_A0.findall(text), tag)
    out["cells"], out["cells_record"] = int(cells), int(cells_rec)
    found = {}
    for m in _RE_CLASS.finditer(text):
        key = (m["ring"] == "same ", int(m["steps"]))
        found.setdefault(key, []).append((m["payload"].rstrip(), int(m["n"]), m["mean"], m["spread"]))
    expected = [(same, k) for same in (True, False) for k in range(N_AZIMUTH_CLASSES_32)]
    extra = sorted(set(found) - set(expected))
    if extra:
        raise ValueError(f"{tag}: unexpected class keys {extra}")
    table = {}
    counts = []
    for key in expected:
        vals = found.get(key, [])
        distinct = sorted(set(vals))
        counts.append((len(vals), len(distinct)))
        if len(distinct) != 1:
            raise ValueError(f"{tag} class {key}: {len(distinct)} distinct matches (count {len(distinct)}), "
                             f"expected exactly 1")
        table[key] = distinct[0]
    _say(f"{tag} parse key class lines: {len(table)} classes, matches per class "
         f"{sorted(set(c[0] for c in counts))}, distinct per class {sorted(set(c[1] for c in counts))}")
    out["table"] = table
    return out


def _ordered_keys():
    return [(same, k) for same in (True, False) for k in range(N_AZIMUTH_CLASSES_32)]


def _p1(text44):
    """(P1): entry 44's parse equals C13's landed record."""
    p = parse_c3_log(text44, tag="entry44")
    mine = [p["table"][k][0] for k in _ordered_keys()]
    hits = sum(_norm(a) == _norm(b) for a, b in zip(mine, C3_TABLE_LINES))
    rr = (p["recip"], p["sigma"])
    _say(f"(P1) entry 44 parse vs C13 C3_TABLE_LINES under _norm: {hits}/{len(C3_TABLE_LINES)}; (ii)/(iii) "
         f"{rr} vs C13 records {(RECIPROCITY_RECORD_STR, SIGMA_MAX_RECORD_STR)} (string equality, ASSERTED)")
    assert len(mine) == len(C3_TABLE_LINES) == N_SYMMETRY_CLASSES_32
    assert hits == len(C3_TABLE_LINES), [(a, b) for a, b in zip(mine, C3_TABLE_LINES) if _norm(a) != _norm(b)][:3]
    assert rr == (RECIPROCITY_RECORD_STR, SIGMA_MAX_RECORD_STR), rr
    return p


def _load(path):
    assert path.exists(), f"{path} absent (gitignored scratch; never re-solve)"
    with np.load(path, allow_pickle=False) as z:
        sv = {k: z[k] for k in z.files}
    return {
        "s": np.array(sv["s64"], dtype=np.complex128),
        "order": [str(x) for x in sv["order"]],
        "az": np.asarray(sv["azimuth_deg"], dtype=float),
        "zc": np.asarray(sv["z_m"], dtype=float),
        "cells": int(sv["cells"]),
        "n_ranks": int(sv["n_ranks"]),
    }


def _rr(d, parsed, tag):
    """(RR) strings and by-value equality against a parsed log."""
    mine = (f"{_reciprocity_ratio(d['s']):.6e}", f"{float(np.linalg.norm(d['s'], 2)):.12f}")
    log = (parsed["recip"], parsed["sigma"])
    by_value = all(float(a) == float(b) for a, b in zip(mine, log))
    _say(f"(RR) {tag}: npz {mine} vs log {log}; by value {by_value}; string identity {mine == log}")
    return mine, by_value


def _ct(d, parsed, tag):
    """(CT): classes from the npz, sizes, and each line against the log (n exact,
    mean .9e and spread .4f by value).  Returns (classes, spreads, n_match, misses)."""
    classes, res = _class_table(d["az"], d["zc"], classify=_azimuth_class_32)
    sizes_ok = len(classes) == N_SYMMETRY_CLASSES_32 and {k: len(v) for k, v in classes.items()} == \
        _expected_class_sizes()
    spreads = _class_spreads(d["s"], classes)
    misses = []
    for key in _ordered_keys():
        sp, cnt, mean = spreads[key]
        _payload, n_log, mean_log, spread_log = parsed["table"][key]
        ok = (cnt == n_log and float(f"{mean:.9e}") == float(mean_log)
              and float(f"{sp * 100:.4f}") == float(spread_log))
        if not ok:
            misses.append((key, (cnt, f"{mean:.9e}", f"{sp * 100:.4f}"), (n_log, mean_log, spread_log)))
    n_match = N_SYMMETRY_CLASSES_32 - len(misses)
    strings = _table_lines(spreads, N_AZIMUTH_CLASSES_32)
    n_str = sum(_norm(a) == _norm(parsed["table"][k][0]) for a, k in zip(strings, _ordered_keys()))
    _say(f"(CT) {tag}: classes {len(classes)}, sizes equal expected {sizes_ok}, worst pair residual "
         f"{res:.3e} deg; lines matching log by value {n_match}/{N_SYMMETRY_CLASSES_32} (string "
         f"identity under _norm {n_str}/{N_SYMMETRY_CLASSES_32}, printed)")
    return classes, spreads, sizes_ok, misses


# --------------------------------------------------------------------------- W0
def test_c51_knob_and_parser(monkeypatch):
    """Ungated: the knob, (P1) on entry 44's tracked log, entry 65's parse, and the
    doubled-table control (two windows' tables in one text ⇒ ValueError, count 2)."""
    monkeypatch.delenv(KNOB_C51, raising=False)
    assert c51_enabled() is False
    monkeypatch.setenv(KNOB_C51, "1")
    assert c51_enabled() is True
    for bad in ("0", "yes", "r2"):
        monkeypatch.setenv(KNOB_C51, bad)
        with pytest.raises(ValueError):
            c51_enabled()
    monkeypatch.delenv(KNOB_C51, raising=False)
    assert str(SAVE_PATH_R2) == "logs/geo38/c3-cu-s64-d1-r2.npz"
    assert str(ENTRY44_D1_PATH) == "logs/geo38/c3-cu-s64-d1.npz"

    text44 = read_log_text(ENTRY44_LOG, errors="strict")
    text65 = read_log_text(ENTRY65_LOG, errors="strict")
    _p1(text44)
    p65 = parse_c3_log(text65, tag="entry65")
    assert len(p65["table"]) == N_SYMMETRY_CLASSES_32
    with pytest.raises(ValueError, match=r"count 2") as exc:
        parse_c3_log(text65 + "\n" + text44, tag="doubled")
    _say(f"(CN-W0) doubled log (entry 65 + entry 44 text) raises: {exc.value}")
    table44 = "\n".join(ln for ln in text44.splitlines() if _RE_CLASS.search(ln))
    with pytest.raises(ValueError, match=r"class .*count 2") as exc:
        parse_c3_log(text65 + "\n" + table44, tag="doubled-table")
    _say(f"(CN-W0) doubled table (entry 65 + entry 44's {len(table44.splitlines())} class lines) raises: "
         f"{exc.value}")


# --------------------------------------------------------------------------- W1
def test_geo38_c51_class_two_rungs():
    if not c51_enabled():
        pytest.skip(f"{KNOB_C51} unset (reads gitignored logs/geo38 npz)")
    t0 = time.perf_counter()
    text44 = read_log_text(ENTRY44_LOG, errors="strict")
    text65 = read_log_text(ENTRY65_LOG, errors="strict")
    text13 = read_log_text(C13_W1_LOG, errors="strict")

    # ---- (P1) ----
    p44 = _p1(text44)
    p65 = parse_c3_log(text65, tag="entry65")

    path1 = WORKSPACE / ENTRY44_D1_PATH
    path2 = WORKSPACE / SAVE_PATH_R2
    d1 = _load(path1)
    d2 = _load(path2)
    for tag, path in (("r1", path1), ("r2", path2)):
        st = path.stat()
        _say(f"{tag} {path}: {st.st_size} B, mtime {time.strftime('%Y-%m-%dT%H:%M:%S%z', time.localtime(st.st_mtime))}")

    # ---- (RR) ----
    rr1, rr1_ok = _rr(d1, p44, "r1 npz vs entry 44")
    rr2, rr2_ok = _rr(d2, p65, "r2 npz vs entry 65")

    # ---- (CT) ----
    cls1, sp1, sz1, miss1 = _ct(d1, p44, "r1 npz vs entry 44")
    cls2, sp2, sz2, miss2 = _ct(d2, p65, "r2 npz vs entry 65")
    for m in miss1 + miss2:
        _say(f"(CT) MISS {m}")

    # ---- (Cus) ----
    order_ok = d2["order"] == EXPECTED_ORDER
    saved_ok = p65["saved"] == str(SAVE_PATH_R2)
    _say(f"(Cus) r2 npz cells {d2['cells']} vs F_CLINICAL_R2_CELL_RECORD {F_CLINICAL_R2_CELL_RECORD} vs log (A0) "
         f"{p65['cells']} (record {p65['cells_record']}); n_ranks {d2['n_ranks']} (expected {EXPECTED_N_RANKS}); "
         f"order {d2['order'][0]}…{d2['order'][-1]} equal P33…P96 {order_ok}; saved line '{p65['saved']}' "
         f"== SAVE_PATH_R2 {saved_ok} (ASSERTED)")
    _say(f"(Cus, printed) r1 npz cells {d1['cells']} vs entry 44 (A0) {p44['cells']}; n_ranks {d1['n_ranks']}; "
         f"saved line '{p44['saved']}'")

    assert rr1_ok and rr2_ok, (rr1, rr2)
    assert sz1 and sz2, (sz1, sz2)
    assert not miss1 and not miss2, (miss1[:3], miss2[:3])
    assert d2["cells"] == F_CLINICAL_R2_CELL_RECORD == p65["cells"] == p65["cells_record"], d2["cells"]
    assert d2["n_ranks"] == EXPECTED_N_RANKS, d2["n_ranks"]
    assert order_ok, d2["order"][:3]
    assert saved_ok, p65["saved"]

    # ---- (C13r) ----
    abs1 = _abs_spreads(d1["s"], cls1)
    abs2 = _abs_spreads(d2["s"], cls2)
    w1 = max(sp1, key=lambda k: sp1[k][0])
    w2 = max(sp2, key=lambda k: sp2[k][0])
    stat1 = f"{abs1[w1] / abs1[(True, 1)]:.4f}"
    stat2 = f"{abs2[w2] / abs2[(True, 1)]:.4f}"
    c13_log = _one("C13 F-clinical abs-spread ratio", _RE_C13.findall(text13), "C13-W1")
    _say(f"(C13r) r1 worst {w1} abs {abs1[w1]:.6e} / (same ring, 1 step) {abs1[(True, 1)]:.6e} = {stat1} vs "
         f"record {C13_STAT_R1_RECORD_STR} (C13-W1.log:237 parse '{c13_log}') (string equality, ASSERTED)")
    assert stat1 == C13_STAT_R1_RECORD_STR == c13_log, (stat1, c13_log)

    # ---- (CN) the r1 npz read as r2 against entry 65 ----
    rr_cn, rr_cn_same = _rr(d1, p65, "CN: r1 npz vs entry 65")
    _c, _s, _sz, miss_cn = _ct(d1, p65, "CN: r1 npz vs entry 65")
    rr_differ = all(float(a) != float(b) for a, b in zip(rr_cn, (p65["recip"], p65["sigma"])))
    _say(f"(CN) r1 npz read as r2: (CT) differing lines {len(miss_cn)}/{N_SYMMETRY_CLASSES_32} (PREDICTED 34; "
         f"ASSERTED >= 1); (RR) both differ {rr_differ} (ASSERTED)")
    assert len(miss_cn) >= 1, "r1 npz reproduces entry 65's table: (CT) blind to the rung"
    assert rr_differ, rr_cn

    # ---- the record (printed only, after the anchor holds; r2 is the COARSER rung) ----
    _say("RECORD (printed, from the log strings): class | r1 rel % | r2 rel % | r2/r1")
    ratios = {}
    for key in _ordered_keys():
        r1s, r2s = p44["table"][key][3], p65["table"][key][3]
        ratios[key] = float(r2s) / float(r1s)
        same, k = key
        _say(f"  {'same ' if same else 'other'} ring, {k:2d} steps | {r1s} | {r2s} | {ratios[key]:.4f}")
    asym = (False, 1)
    _say(f"asymmetric class (other ring, 1 step): r2 {p65['table'][asym][3]}% / r1 {p44['table'][asym][3]}% "
         f"= {ratios[asym]:.4f} (drafter's arithmetic {PRED_ASYM_RATIO}, printed)")
    vals = np.array([ratios[k] for k in _ordered_keys()])
    kmin = min(ratios, key=ratios.get)
    kmax = max(ratios, key=ratios.get)
    ge2 = [k for k in _ordered_keys() if ratios[k] >= 2.0]
    _say(f"r2/r1 median {np.median(vals):.4f} (drafter {PRED_MEDIAN}); min {ratios[kmin]:.4f} at {kmin}; max "
         f"{ratios[kmax]:.4f} at {kmax}; count >= 2x {len(ge2)} {ge2} (drafter {PRED_COUNT_GE_2X})")
    _say("RECORD (printed, from the npz after custody): class | abs spread r1 | abs spread r2 | r2/r1")
    for key in _ordered_keys():
        same, k = key
        _say(f"  {'same ' if same else 'other'} ring, {k:2d} steps | {abs1[key]:.6e} | {abs2[key]:.6e} | "
             f"{abs2[key] / abs1[key]:.4f}")
    _say(f"C13 statistic (abs spread worst / same ring 1 step): r1 {stat1} [{w1}], r2 {stat2} [{w2}]; "
         f"(other ring, 1 step) / (same ring, 1 step): r1 {abs1[asym] / abs1[(True, 1)]:.4f}, "
         f"r2 {abs2[asym] / abs2[(True, 1)]:.4f}")
    _say(f"body {time.perf_counter() - t0:.2f} s")
