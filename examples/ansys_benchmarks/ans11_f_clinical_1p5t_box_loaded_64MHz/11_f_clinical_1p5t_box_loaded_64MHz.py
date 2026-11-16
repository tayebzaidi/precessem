"""`ANS-11` — the `ans:11` record page: our degree-1 column of the
F-clinical-1.5T ASTM-box-loaded SPEC, read from entry 52's tracked log (no solve).

§9 item 833 (docs/planning/chunks/ANS-11.md, "Review ruling 2026-11-05 03:00"):
a **record reader** on the `ans:4` item 417 precedent.  It parses the tracked C6
log ``docs/testing/logs/20261103T080009Z_GEO-38-C6.log`` (``xl`` entry 52, `-n 8`,
868 s) and builds no mesh.  The gitignored run JSON (log `:12311`) is not
reproducible from a fresh clone, so the tracked log is the source.

**A degree-1 record beside C19's degree-2 record** (C37, §9 item 983): the degree-2
column is parsed from ``xxl`` entry 62's tracked log
``docs/testing/logs/20261114T080009Z_GEO-38-C19.log`` (`-n 16`, degree 2, SPEC v2's
node 185.7629 pF) by the same ``parse_c6_log``; entry 62's host JSONs are a custody
check (J62) only and the committed render never reads them.  The degree-1 → degree-2
move per row and the 10-24 exit's 2 % order rule are printed records.  Neither
column is adjudicated (Done-when (2)), and the AED column is blank (AED values are
private; none is read or written here).

Anchors (asserted here and in ``tests/unit/test_ans11_record_page.py``):

* (X) the 18 lattice labels equal C6's `_lattice()` in order; every SPEC-quoted
  value (``SPEC.md:120–126``, parsed from the SPEC, never retyped) equals its log
  token; the SPEC's 36 six-digit lattice entries (``:132–134``) equal
  ``f"{float(tok):.6e}"``;
* (D) derived quantities from the parsed tokens, each at half a unit in its last
  printed digit: `_p_acc_spec` against P_acc, SAR10g(c*)/P_acc, P_box/P_acc, the
  18-point lattice mean (and its per-√P_acc value), argmax(SAR) = c*, the c* SAR
  string = its lattice line's string, S₁₂ string = S₂₁ string.

Printed, not asserted (rule (e)): (N-sense) the cw weights (1, +1j) through the
same `_p_acc_spec`.

Run::

    ./scripts/run_examples.sh -e ans:11 -n 1 -t 180
"""

from __future__ import annotations

import gzip
import json
import math
import os
import re
import sys
from pathlib import Path

import numpy as np

# The runner puts only ``src`` on PYTHONPATH; the repo root goes on sys.path so
# C6's module imports (the ans:4 / ans:9 pattern).
_REPO_ROOT = Path(__file__).resolve().parents[3]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from tests.validation.test_geo38_f_clinical_box_loaded import (  # noqa: E402
    BRIEF_LATTICE_N,
    C6_JSON_DIR,
    C6A_LATTICE_DZ_M,
    K8_JSON_KEYS,
    _lattice,
    _p_acc_spec,
)
from tests.validation.test_geo36_rod_f_clinical_drive import (  # noqa: E402
    I36_JSON_DIR,  # noqa: F401  (the root `_i36_json_path` writes under; entry 62's call 2)
    _i36_json_path,
)
from tests.validation.test_geo38_c18_degree2_c_mode import ROOT_D2_PF  # noqa: E402
from tests.validation.test_geo38_f_clinical_esr_root import FREQUENCY_HZ_C4  # noqa: E402
from tests.validation.test_wf10_f_human_copper_esr_root import (  # noqa: E402
    ESR_OHM,
    _thevenin_feed_impedance,
    _z_rc,
)

CASE_DIR = Path(__file__).resolve().parent
SPEC_PATH = CASE_DIR / "SPEC.md"
METRICS_JSON = CASE_DIR / "metrics.json"
COMPARISON_MD = CASE_DIR / "COMPARISON.md"
#: Entry 52's C6 log (tracked; retention may gzip it — read either suffix, cite the base name).
SOURCE_LOG_BASENAME = "20261103T080009Z_GEO-38-C6.log"
SOURCE_LOG = _REPO_ROOT / "docs" / "testing" / "logs" / SOURCE_LOG_BASENAME
SKIP_PRIVATE_ENV = "FEM_EM_ANS_SKIP_PRIVATE"
Z0_OHM = 50.0
CCW = (1.0, -1.0j)
CW = (1.0, 1.0j)

_TAG = r"\[GEO-38 C6a\] "
_C = r"([+-]?\d\.\d+e[+-]\d+[+-]\d\.\d+e[+-]\d+j)"  # a printed complex token
_F = r"([+-]?\d+\.\d+(?:e[+-]\d+)?)"  # a printed real token
_LABEL = r"(x[+-]\d+mm dz[+-]\d+mm)"

#: Singleton lines: key -> (regex, field names).  ``re.search``, never ``startswith``
#: (interleaved pytest output sits on log lines).
_SINGLETONS = {
    "node": (_TAG + r"\(Z\) node (\d+\.\d+) pF: Z_RC " + _C + r" Ohm, Z_th " + _C
             + r" Ohm, V_s " + _C + r" V",
             ("c_pf", "z_rc", "z_th", "v_s")),
    "s2": (_TAG + r"kinds \[[^\]]*\]; S2 " + _C + r", " + _C + r"; " + _C + r", " + _C
           + r"; \(R1\)",
           ("s11", "s12", "s21", "s22")),
    "y_in_ccw": (_TAG + r"Y_in ccw \([^:]*\): P33 " + _F + r" " + _F + r"j S, P41 "
                 + _F + r" " + _F + r"j S",
                 ("p33_re", "p33_im", "p41_re", "p41_im")),
    "y_in_cw": (_TAG + r"Y_in cw \([^:]*\): P33 " + _F + r" " + _F + r"j S, P41 "
                + _F + r" " + _F + r"j S",
                ("p33_re", "p33_im", "p41_re", "p41_im")),
    "means": (_TAG + r"mean \|B1\+\|_ccw " + _F + r" T over (\d+) / (\d+) tag-3 centroids, "
              r"mean \|B1-\| " + _F + r" T, \(Sense\) " + _F + r" \(printed\); CV " + _F,
              ("b1p_ccw", "n_used", "n_total", "b1m_ccw", "sense", "cv")),
    "c_star": (_TAG + r"c\* " + _LABEL + r": SAR10g " + _F + r" W/kg \(ccw\), \(1,\+1j\) ball "
               + _F + r" W/kg; r " + _F + r" m, dz " + _F + r" m, az " + _F + r" deg",
               ("label", "sar10g", "sar10g_cw", "r_m", "dz_m", "az_deg")),
    "normalisers": (_TAG + r"normalisers: P_box " + _F + r" W, P_src \(field partition, ccw\) "
                    + _F + r" W, P_acc \(lumped line cut, ruling 2026-10-12\) " + _F
                    + r" W; P_box/P_src " + _F + r", P_box/P_acc " + _F
                    + r"; SAR10g\(c\*\) per box watt " + _F + r", per P_acc " + _F
                    + r", per P_src " + _F + r" W/kg/W; mean \|B1\+\|_ccw / sqrt\(W\): box "
                    + _F + r", P_acc " + _F + r" T",
                    ("p_box", "p_src", "p_acc", "p_box_over_p_src", "p_box_over_p_acc",
                     "sar_per_box_w", "sar_per_p_acc", "sar_per_p_src",
                     "b1_per_sqrt_box_w", "b1_per_sqrt_p_acc")),
    "wrote": (_TAG + r"wrote (\S+)", ("path",)),
    "price": (_TAG + r"price: build " + _F + r" s; columns " + _F + r", " + _F
              + r" s; readouts " + _F + r" s; .*window " + _F + r" s",
              ("build_s", "column_p33_s", "column_p41_s", "readouts_s", "window_s")),
}
_LATTICE_RE = re.compile(_TAG + r"lattice " + _LABEL + r": \|B1\+\|_ccw " + _F
                         + r" T, SAR10g " + _F + r" W/kg, mass " + _F + r" kg")


# ---------------------------------------------------------------------------
# Readers (pure).
# ---------------------------------------------------------------------------
def _read_log(path) -> str:
    """A harness log's text, from ``path`` or, once retention has gzipped it,
    ``path`` + ``.gz`` (`ans:4` `_read_log`, `OPS-75`)."""
    path = str(path)
    if not os.path.exists(path) and os.path.exists(f"{path}.gz"):
        path = f"{path}.gz"
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt", errors="replace") as fh:
        return fh.read()


def parse_c6_log(path) -> dict:
    """Pure: entry 52's C6a readout lines, every token kept as its printed
    **string** with its 1-based line number.  Raises ``ValueError`` if a
    singleton line is missing or repeated, or the lattice count is not
    ``BRIEF_LATTICE_N``."""
    lines = _read_log(path).splitlines()
    compiled = {k: re.compile(rx) for k, (rx, _) in _SINGLETONS.items()}
    out: dict = {"log": Path(str(path)).name.removesuffix(".gz"), "lattice": []}
    for lineno, line in enumerate(lines, start=1):
        for key, rx in compiled.items():
            m = rx.search(line)
            if m:
                if key in out:
                    raise ValueError(f"{key}: second match at :{lineno} "
                                     f"(first :{out[key]['line']})")
                out[key] = {"line": lineno, **dict(zip(_SINGLETONS[key][1], m.groups()))}
        m = _LATTICE_RE.search(line)
        if m:
            out["lattice"].append({"line": lineno, "label": m.group(1), "b1p_ccw": m.group(2),
                                   "sar10g": m.group(3), "mass_kg": m.group(4)})
    missing = [k for k in _SINGLETONS if k not in out]
    if missing:
        raise ValueError(f"log lines not found: {missing}")
    if len(out["lattice"]) != BRIEF_LATTICE_N:
        raise ValueError(f"lattice lines: {len(out['lattice'])} != BRIEF_LATTICE_N "
                         f"{BRIEF_LATTICE_N}")
    return out


def _spec_num(s: str) -> str:
    return s.replace("−", "-").replace(" ", "")


def parse_spec(path=SPEC_PATH) -> dict:
    """Pure: the SPEC's quoted entry-52 values (``SPEC.md:120–126``) and its
    six-digit lattice table (``:132–134``), as strings with line numbers."""
    lines = Path(path).read_text().splitlines()
    cnum = r"(−?\d\.\d+e[+-]\d+ [+−] \d\.\d+e[+-]\d+j)"
    num = r"(−?\d+\.\d+(?:e[+-]\d+)?)"
    rows = {
        "s2": (r"^\| S₁₁, S₁₂ = S₂₁, S₂₂ \| " + cnum + "; " + cnum + "; " + cnum,
               ("s11", "s12", "s22")),
        "y_in_ccw": (r"^\| Y_in ccw, P33 / P41 \[S\] \| " + cnum + " / " + cnum,
                     ("p33", "p41")),
        "means": (r"^\| mean `\\\|B₁⁺\\\|`_ccw over 400 gel-cell centroids / mean `\\\|B₁⁻\\\|` \| "
                  + num + " T / " + num + " T; ratio " + num + "; CV " + num,
                  ("b1p_ccw", "b1m_ccw", "sense", "cv")),
        "lattice_mean": (r"^\| 18-point lattice mean `\\\|B₁⁺\\\|`_ccw \| " + num + " T",
                         ("b1p_ccw",)),
        "c_star": (r"^\| c\\\* \| x " + num + ", z " + num + r" \(r " + num + ", az "
                   + r"(\d+)°\): SAR10g " + num + r" W/kg \(ccw\); cw ball " + num,
                   ("x_m", "z_m", "r_m", "az_deg", "sar10g", "sar10g_cw")),
        "powers": (r"^\| P_box / P_src / P_acc \[W\] \| " + num + " / " + num + " / " + num,
                   ("p_box", "p_src", "p_acc")),
        "per_p_acc": (r"^\| SAR10g\(c\\\*\) per P_acc; mean `\\\|B₁⁺\\\|`/√P_acc \| " + num
                      + " W/kg/W; " + num + r" T/√W \(lattice mean: " + num + r", derived\)",
                      ("sar_per_p_acc", "b1_per_sqrt_p_acc", "lattice_mean_per_sqrt_p_acc")),
    }
    out: dict = {"lattice": []}
    for lineno, line in enumerate(lines, start=1):
        for key, (rx, names) in rows.items():
            m = re.search(rx, line)
            if m:
                if key in out:
                    raise ValueError(f"SPEC {key}: second match at :{lineno}")
                out[key] = {"line": lineno,
                            **{n: _spec_num(g) for n, g in zip(names, m.groups())}}
        m = re.match(r"^\| ([−+]?0(?:\.105)?) \| (.*) \|$", line)
        if m:
            cells = [c.strip() for c in m.group(2).split("|")]
            pairs = [tuple(_spec_num(v) for v in c.split(" / ")) for c in cells]
            if len(pairs) != len(C6A_LATTICE_DZ_M) or any(len(p) != 2 for p in pairs):
                raise ValueError(f"SPEC lattice row :{lineno} has {len(pairs)} cells")
            out["lattice"].append({"line": lineno, "x_m": _spec_num(m.group(1)), "cells": pairs})
    missing = [k for k in rows if k not in out]
    if missing:
        raise ValueError(f"SPEC rows not found: {missing}")
    if len(out["lattice"]) != 3:
        raise ValueError(f"SPEC lattice rows: {len(out['lattice'])} != 3")
    return out


# ---------------------------------------------------------------------------
# Checks (pure).  Each returns ``[(name, ok, detail)]``.
# ---------------------------------------------------------------------------
def check_extraction(parsed: dict, spec: dict) -> list:
    res = []
    labels = [r["label"] for r in parsed["lattice"]]
    want = [lab for lab, _ in _lattice()]
    res.append(("(X) lattice labels == _lattice()", labels == want,
                f"{len(labels)} labels, first mismatch "
                f"{next((i for i, (a, b) in enumerate(zip(labels, want)) if a != b), None)}"))

    def eq(name, spec_s, log_val, kind=float):
        ok = kind(spec_s) == log_val
        res.append((f"(X) {name}", ok, f"SPEC {spec_s} vs log {log_val!r}"))

    s2, ss = parsed["s2"], spec["s2"]
    for k in ("s11", "s12", "s22"):
        eq(f"S {k}", ss[k], complex(s2[k]), complex)
    for p in ("p33", "p41"):
        y = parsed["y_in_ccw"]
        eq(f"Y_in ccw {p}", spec["y_in_ccw"][p],
           complex(float(y[f"{p}_re"]), float(y[f"{p}_im"])), complex)
    for k in ("b1p_ccw", "b1m_ccw", "sense", "cv"):
        eq(f"means {k}", spec["means"][k], float(parsed["means"][k]))
    cs, sc = parsed["c_star"], spec["c_star"]
    for k in ("sar10g", "sar10g_cw", "r_m", "az_deg"):
        eq(f"c* {k}", sc[k], float(cs[k]))
    eq("c* z (= dz, box centred)", sc["z_m"], float(cs["dz_m"]))
    lab_pt = dict(_lattice())[cs["label"]]
    res.append(("(X) c* x (from the c* label's lattice point)",
                math.isclose(float(sc["x_m"]), lab_pt[0], abs_tol=1e-12),
                f"SPEC {sc['x_m']} vs lattice point x {lab_pt[0]!r}"))
    nm = parsed["normalisers"]
    for k in ("p_box", "p_src", "p_acc"):
        eq(f"P {k}", spec["powers"][k], float(nm[k]))
    eq("SAR10g(c*) per P_acc", spec["per_p_acc"]["sar_per_p_acc"], float(nm["sar_per_p_acc"]))
    eq("mean |B1+|/sqrt(P_acc)", spec["per_p_acc"]["b1_per_sqrt_p_acc"],
       float(nm["b1_per_sqrt_p_acc"]))
    # The 36 six-digit lattice entries.
    by_label = {r["label"]: r for r in parsed["lattice"]}
    pts = _lattice()
    bad = []
    for row in spec["lattice"]:
        x = float(row["x_m"])
        for iz, (b_s, s_s) in enumerate(row["cells"]):
            dz = C6A_LATTICE_DZ_M[iz]
            hits = [lab for lab, p in pts if math.isclose(p[0], x, abs_tol=1e-12)
                    and math.isclose(p[2], dz, abs_tol=1e-12)]
            if len(hits) != 1:
                bad.append((row["line"], iz, "no unique lattice point"))
                continue
            r = by_label[hits[0]]
            for spec_s, tok in ((b_s, r["b1p_ccw"]), (s_s, r["sar10g"])):
                if spec_s != f"{float(tok):.6e}":
                    bad.append((row["line"], hits[0], spec_s, tok))
    n = sum(len(r["cells"]) * 2 for r in spec["lattice"])
    res.append(("(X) SPEC six-digit lattice entries == f'{tok:.6e}'", n == 36 and not bad,
                f"{n} entries, mismatches {bad[:3]}"))
    return res


def derived(parsed: dict) -> dict:
    """Floats from the parsed tokens (never from the SPEC)."""
    s = parsed["s2"]
    s2 = np.array([[complex(s["s11"]), complex(s["s12"])],
                   [complex(s["s21"]), complex(s["s22"])]])
    v_s = complex(parsed["node"]["v_s"])
    nm = parsed["normalisers"]
    p_acc_log, p_box = float(nm["p_acc"]), float(nm["p_box"])
    b1 = np.array([float(r["b1p_ccw"]) for r in parsed["lattice"]])
    sar = np.array([float(r["sar10g"]) for r in parsed["lattice"]])
    lat_mean = float(b1.mean())
    return {
        "p_acc_ccw": _p_acc_spec(s2, CCW, v_s, Z0_OHM),
        "p_acc_cw": _p_acc_spec(s2, CW, v_s, Z0_OHM),
        "sar_cstar_per_p_acc": float(parsed["c_star"]["sar10g"]) / p_acc_log,
        "p_box_per_p_acc": p_box / p_acc_log,
        "lattice_mean_b1p": lat_mean,
        "lattice_mean_b1p_per_sqrt_p_acc": lat_mean / math.sqrt(p_acc_log),
        "argmax_label": parsed["lattice"][int(np.argmax(sar))]["label"],
        "b1p_per_sqrt_p_acc": [float(v) / math.sqrt(p_acc_log) for v in b1],
        "sar_per_p_acc": [float(v) / p_acc_log for v in sar],
    }


def check_derived(parsed: dict, spec: dict, d: dict) -> list:
    res = []

    def near(name, val, target_s, tol, rel=False):
        t = float(target_s)
        err = abs(val - t) / abs(t) if rel else abs(val - t)
        res.append((f"(D) {name}", err <= tol,
                    f"{val:.6e} vs printed {target_s}: {'rel ' if rel else ''}{err:.3e} "
                    f"(<= {tol:.0e})"))

    near("_p_acc_spec(S2, (1,-1j), V_s, 50) vs SPEC P_acc", d["p_acc_ccw"],
         spec["powers"]["p_acc"], 1e-8, rel=True)
    near("SAR10g(c*)/P_acc vs SPEC", d["sar_cstar_per_p_acc"],
         spec["per_p_acc"]["sar_per_p_acc"], 5e-10)
    near("P_box/P_acc vs log's printed ratio", d["p_box_per_p_acc"],
         parsed["normalisers"]["p_box_over_p_acc"], 5e-9)
    near("18-point lattice mean |B1+| vs SPEC", d["lattice_mean_b1p"],
         spec["lattice_mean"]["b1p_ccw"], 5e-14)
    near("lattice mean / sqrt(P_acc) vs SPEC", d["lattice_mean_b1p_per_sqrt_p_acc"],
         spec["per_p_acc"]["lattice_mean_per_sqrt_p_acc"], 5e-13)
    cs = parsed["c_star"]
    res.append(("(D) argmax(lattice SAR10g) == c* label", d["argmax_label"] == cs["label"],
                f"argmax {d['argmax_label']} vs c* {cs['label']}"))
    lat = {r["label"]: r for r in parsed["lattice"]}.get(cs["label"])
    res.append(("(D) c* SAR string == its lattice line's string",
                lat is not None and lat["sar10g"] == cs["sar10g"],
                f"c* {cs['sar10g']} vs lattice {lat and lat['sar10g']}"))
    s = parsed["s2"]
    res.append(("(D) S12 string == S21 string", s["s12"] == s["s21"],
                f"{s['s12']} vs {s['s21']}"))
    return res


# ---------------------------------------------------------------------------
# C36 (§9 item 871): the entry-62 reader, the degree-2 mapping, (K) and (D).
# ---------------------------------------------------------------------------
#: COMPARISON.md's Quantity labels, keyed with their row number (rows 2 and 4 share
#: the "per √W" label).
Q_S11, Q_S12, Q_S22 = (1, "S₁₁"), (1, "S₁₂ = S₂₁"), (1, "S₂₂")
Q_Y33, Q_Y41 = (1, "Y_in ccw P33 [S]"), (1, "Y_in ccw P41 [S]")
Q_LAT = (2, "18-point lattice mean `\\|B₁⁺\\|`_ccw [T] (the comparison row)")
Q_LAT_P = (2, "the same per √W of P_acc [T/√W]")
Q_CSTAR = (3, "c\\* (largest of the 18 SAR10g)")
Q_SAR = (3, "SAR10g(c\\*) [W/kg]; per W of P_acc")
Q_CW = (3, "cw (1, +1j) ball at c\\* [W/kg] (record)")
Q_GRID = (3, "box-wide 10 g maximum and location (record)")
Q_MEAN = (4, "400-centroid gel mean `\\|B₁⁺\\|`_ccw [T] (second record)")
Q_MEAN_P = (4, "the same per √W of P_acc [T/√W]")
Q_P = (5, "P_box / P_src / P_acc [W]")
Q_GEL = (5, "gel fraction P_box / P_acc")
Q_CU = (5, "copper / capacitor-ESR fractions under ccw")
Q_SENSE = (6, "mean `\\|B₁⁺\\|` / mean `\\|B₁⁻\\|` over the gel (sense)")
Q_ROWS = (Q_S11, Q_S12, Q_S22, Q_Y33, Q_Y41, Q_LAT, Q_LAT_P, Q_CSTAR, Q_SAR, Q_CW, Q_GRID,
          Q_MEAN, Q_MEAN_P, Q_P, Q_GEL, Q_CU, Q_SENSE)
#: The node C19 runs at (SPEC v2, ``SPEC.md:14``).
C19_NODE_PF = "185.7629"
#: SPEC tokens COMPARISON.md does not carry in the SPEC's form ((K) prints them).
K_NOT_CARRIED = ("c_star.az_deg", "c_star.r_m", "c_star.x_m", "c_star.z_m", "means.b1m_ccw",
                 "means.cv")


def read_c19(c_pf=C19_NODE_PF, c6_dir=None, i36_dir=None) -> dict:
    """Pure reader over entry 62's two JSONs (absent gives ``None``, no error).

    C6 at degree 2: ``c6-node-<c_pf>-d2.json`` (``test_geo38_f_clinical_box_loaded.py:739–741``);
    I36b at degree 2: ``_i36_json_path(c_pf, 2)``.  A present file is validated; anything
    else raises ``ValueError`` naming the file and the key."""
    c6_path = os.path.join(c6_dir or C6_JSON_DIR, f"c6-node-{c_pf}-d2.json")
    i36_path = _i36_json_path(float(c_pf), 2)
    if i36_dir is not None:
        i36_path = os.path.join(i36_dir, os.path.basename(i36_path))
    out = {"c6": None, "i36": None, "paths": [c6_path, i36_path]}
    if os.path.exists(c6_path):
        with open(c6_path) as fh:
            rec = json.load(fh)
        keys, k8 = set(rec), set(K8_JSON_KEYS)
        if keys not in (k8, k8 | {"grid"}):
            raise ValueError(f"{c6_path}: key set != K8_JSON_KEYS (+grid); missing "
                             f"{sorted(k8 - keys)}, extra {sorted(keys - k8 - {'grid'})}")
        if f"{rec['c_pf']:.4f}" != c_pf:
            raise ValueError(f"{c6_path}: key c_pf {rec['c_pf']:.4f} != {c_pf}")
        out["c6"] = rec
    if os.path.exists(i36_path):
        with open(i36_path) as fh:
            rec = json.load(fh)
        if rec.get("degree") != 2:
            raise ValueError(f"{i36_path}: key degree {rec.get('degree')!r} != 2")
        if not isinstance(rec.get("c_pf"), (int, float)) or f"{rec['c_pf']:.4f}" != c_pf:
            raise ValueError(f"{i36_path}: key c_pf {rec.get('c_pf')!r} != {c_pf}")
        out["i36"] = rec
    return out


def degree2_cells(c19) -> dict:
    """``{(row, Quantity): cell}`` for the degree-2 column; every cell ``""`` without a
    C6 record.  Raw values at the C6 module's print formats, derived ones at the page's
    (<= 7 significant digits).  The I36 record maps to no row (the SPEC has no rod row)."""
    cells = {q: "" for q in Q_ROWS}
    rec = (c19 or {}).get("c6")
    if rec is None:
        return cells
    s2 = [[complex(*v) for v in row] for row in rec["s2"]]
    yc = [complex(*v) for v in rec["y_in"]["ccw"]]
    p_acc = float(rec["p_acc"])
    lat = float(np.mean([r["b_plus_ccw"] for r in rec["lattice"]]))
    mbp = float(rec["mean_b_plus"]["ccw"])

    def ycell(y):
        return f"{y.real:+.9e} {y.imag:+.9e}j (Im {'−' if y.imag < 0 else '+'})"

    cells.update({
        Q_S11: f"{s2[0][0]:.10e}", Q_S12: f"{s2[0][1]:.10e}", Q_S22: f"{s2[1][1]:.10e}",
        Q_Y33: ycell(yc[0]), Q_Y41: ycell(yc[1]),
        Q_LAT: f"{lat:.5e}", Q_LAT_P: f"{lat / math.sqrt(p_acc):.4e}",
        Q_CSTAR: str(rec["star_label"]),
        Q_SAR: f"{rec['sar_star']:.9e}; {rec['sar_star'] / p_acc:.6e}",
        Q_CW: f"{rec['sar_star_cw']:.9e}",
        Q_MEAN: f"{mbp:.9e}", Q_MEAN_P: f"{mbp / math.sqrt(p_acc):.6e}",
        Q_P: f"{rec['p_box']:.9e} / {rec['p_src']:.9e} / {p_acc:.9e}",
        Q_GEL: f"{rec['p_box'] / p_acc:.6e}",
        Q_SENSE: f"{mbp / float(rec['mean_b_minus']['ccw']):.6f}",
    })
    if "grid" in rec:
        cells[Q_GRID] = f"{rec['grid']['star_label']}: {rec['grid']['sar_star']:.9e}"
    return cells


def _k_norm(t: str) -> str:
    """(K)'s N: drop every ``(…)`` parenthetical, delete spaces, − to -, strip one leading +."""
    prev = None
    while prev != t:
        prev, t = t, re.sub(r"\([^()]*\)", "", t)
    t = t.replace(" ", "").replace("−", "-")
    return t[1:] if t.startswith("+") else t


def comparison_cells(text: str) -> dict:
    """Pure: COMPARISON's main table ``{(row, Quantity): (degree-1, degree-2)}`` and the
    lattice table ``{label: (|B1+| cell, SAR cell)}``."""
    main, lattice = {}, {}
    for line in text.splitlines():
        if not line.startswith("| ") or line.startswith("|---"):
            continue
        cells = [c.strip() for c in line[1:-1].split(" |")]
        if len(cells) == 5 and cells[0].isdigit():
            main[(int(cells[0]), cells[1])] = (cells[2], cells[3])
        elif len(cells) == 6 and re.fullmatch(r"x[+-]\d+mm dz[+-]\d+mm", cells[0]):
            lattice[cells[0]] = (cells[2], cells[4])
    return {"main": main, "lattice": lattice}


def check_k(spec: dict, comparison_text: str) -> list:
    """(K): every SPEC token occurs, under N, in its COMPARISON degree-1 cell; the SPEC's
    36 lattice entries equal ``f"{float(cell):.6e}"``; the not-carried list is the named 6."""
    c = comparison_cells(comparison_text)
    m = c["main"]
    sp = spec
    want = [
        ("s2.s11", sp["s2"]["s11"], Q_S11), ("s2.s12", sp["s2"]["s12"], Q_S12),
        ("s2.s22", sp["s2"]["s22"], Q_S22),
        ("y_in_ccw.p33", sp["y_in_ccw"]["p33"], Q_Y33),
        ("y_in_ccw.p41", sp["y_in_ccw"]["p41"], Q_Y41),
        ("means.b1p_ccw", sp["means"]["b1p_ccw"], Q_MEAN),
        ("means.sense", sp["means"]["sense"], Q_SENSE),
        ("lattice_mean.b1p_ccw", sp["lattice_mean"]["b1p_ccw"], Q_LAT),
        ("c_star.sar10g", sp["c_star"]["sar10g"], Q_SAR),
        ("c_star.sar10g_cw", sp["c_star"]["sar10g_cw"], Q_CW),
        ("powers.p_box", sp["powers"]["p_box"], Q_P),
        ("powers.p_src", sp["powers"]["p_src"], Q_P),
        ("powers.p_acc", sp["powers"]["p_acc"], Q_P),
        ("per_p_acc.sar_per_p_acc", sp["per_p_acc"]["sar_per_p_acc"], Q_SAR),
        ("per_p_acc.b1_per_sqrt_p_acc", sp["per_p_acc"]["b1_per_sqrt_p_acc"], Q_MEAN_P),
        ("per_p_acc.lattice_mean_per_sqrt_p_acc", sp["per_p_acc"]["lattice_mean_per_sqrt_p_acc"],
         Q_LAT_P),
    ]
    res = []
    for name, tok, key in want:
        cell = m.get(key, (None,))[0]
        ok = cell is not None and _k_norm(tok) in _k_norm(cell)
        res.append((f"(K) {name}", ok, f"N(SPEC) {_k_norm(tok)} in N(row {key[0]} degree-1 "
                                      f"cell) {None if cell is None else _k_norm(cell)}"))
    bad, n = [], 0
    for row in sp["lattice"]:
        x = float(row["x_m"])
        for iz, pair in enumerate(row["cells"]):
            lab = f"x{round(x * 1000):+d}mm dz{round(C6A_LATTICE_DZ_M[iz] * 1000):+d}mm"
            got = c["lattice"].get(lab)
            for j, spec_s in enumerate(pair):
                n += 1
                if got is None or spec_s != f"{float(got[j]):.6e}":
                    bad.append((lab, j, spec_s, got and got[j]))
    res.append(("(K) SPEC 36 lattice entries == f'{float(cell):.6e}'", n == 36 and not bad,
                f"{n} entries, mismatches {bad[:3]}"))
    carried = {name for name, _, _ in want}
    tokens = {f"{k}.{f}" for k in ("means", "c_star") for f in sp[k] if f != "line"}
    not_carried = sorted(tokens - carried)
    res.append(("(K) not-carried SPEC tokens == the named 6",
                not_carried == sorted(K_NOT_CARRIED), f"not carried: {not_carried}"))
    return res


def node_figures(c_f, frequency_hz=FREQUENCY_HZ_C4, r_ohm=ESR_OHM, z0=Z0_OHM) -> dict:
    """The SPEC's R_p ‖ C figures at C = ``c_f`` [F]; Z_RC / Z_th from WF-10's module
    (``_z_rc`` ``:188–191``, ``_thevenin_feed_impedance`` ``:144``), never re-implemented."""
    w = 2.0 * math.pi * float(frequency_hz)
    q = 1.0 / (w * c_f * r_ohm)
    z_rc = complex(_z_rc(c_f, frequency_hz, r_ohm))
    return {"q_c": q, "r_p": 1.0 / (w * w * c_f * c_f * r_ohm), "r_p_exact": r_ohm * (1.0 + q * q),
            "c_p_pf": c_f / (1.0 + q ** -2) * 1e12, "z_rc": z_rc,
            "z_th": complex(_thevenin_feed_impedance(z_rc, z0))}


def parse_spec_node(path=SPEC_PATH) -> dict:
    """Pure: SPEC v2's node strings (``:14``, ``:75``) and the R_p ‖ C / Z_th figures at
    ``:75`` / ``:77`` (v2, and the v1 parenthetical beside them), by regex, as strings."""
    text = Path(path).read_text()
    f = r"(\d+\.\d+)"
    pats = {
        "node_14": r"ruled 2026-11-07 \(weekly, SPEC v2\): C = " + f + " pF",
        "node_75": r"C = " + f + r" pF exactly\*\*",
        "v2_z_rc": r"`Z_RC = 0\.2 − " + f + r"j Ω` at 64 MHz",
        "v2_rp": r"R_p = 1/\(ω²C²r\) = " + f + r" Ω ‖ " + f + r" pF\*\*, Q_C = " + f
                 + r" \(exact C_p = " + f + " pF",
        "v1": r"R_p = " + f + r" Ω ‖ " + f + r" pF, Q_C = " + f + r", C_p = " + f
              + r" pF, `Z_RC = 0\.2 − " + f + r"j Ω`",
        "v2_z_th": r"Z_th = 50 ‖ Z_RC = " + f + " − " + f + r"j Ω` at the v2 node",
        "v1_z_th": r"record column's value was `" + f + " − " + f + r"j Ω`",
    }
    out = {}
    for k, rx in pats.items():
        hits = re.findall(rx, text)
        if len(hits) != 1:
            raise ValueError(f"SPEC node pattern {k}: {len(hits)} matches")
        out[k] = hits[0] if isinstance(hits[0], tuple) else (hits[0],)
    return out


def _at(val: float, s: str) -> bool:
    """``val`` printed at ``s``'s own decimals equals ``s``."""
    return f"{val:.{len(s.split('.')[1])}f}" == s


def _ulps(val: float, s: str) -> int:
    """``val - s`` in units of ``s``'s last printed digit (an integer; leak-safe)."""
    return round((val - float(s)) * 10 ** len(s.split(".")[1]))


def check_d(sn: dict, parsed: dict) -> tuple:
    """(D-node), (D-v1), (D-v2) as ``[(name, ok, detail)]``, and the printed-only
    node-C lines (rule (e)).  Details print derived values at <= 7 digits."""
    res = [("(D-node) SPEC node :14 == :75 == f'{ROOT_D2_PF:.4f}'",
            sn["node_14"][0] == sn["node_75"][0] == f"{ROOT_D2_PF:.4f}",
            f":14 {sn['node_14'][0]}, :75 {sn['node_75'][0]}, ROOT_D2_PF {ROOT_D2_PF:.4f}")]

    def figs(name, d, rp, q, cp, zrc, zth):
        for lab, v, st in (("R_p", d["r_p"], rp), ("Q_C", d["q_c"], q), ("C_p", d["c_p_pf"], cp),
                           ("-Im Z_RC", -d["z_rc"].imag, zrc), ("Re Z_th", d["z_th"].real, zth[0]),
                           ("-Im Z_th", -d["z_th"].imag, zth[1])):
            res.append((f"{name} {lab} at the SPEC's printed digits", _at(v, st),
                        f"{v:.7g} vs SPEC {st} ({_ulps(v, st):+d} in its last digit)"))

    tok = parsed["node"]
    d1 = node_figures(float(tok["c_pf"]) * 1e-12)
    zrc_s, zth_s = f"{d1['z_rc']:.9e}", f"{d1['z_th']:.9e}"
    res.append(("(D-v1) f'{z_rc:.9e}' == the C6 log's node Z_RC token", zrc_s == tok["z_rc"],
                f"at C {tok['c_pf']} pF, log `:{tok['line']}`"))
    res.append(("(D-v1) f'{z_th:.9e}' == the C6 log's node Z_th token", zth_s == tok["z_th"],
                f"at C {tok['c_pf']} pF, log `:{tok['line']}`"))
    v1 = sn["v1"]
    res.append(("(D-v1) SPEC v1 parenthetical C == the log's node", v1[1] == tok["c_pf"],
                f"SPEC {v1[1]} vs log {tok['c_pf']}"))
    figs("(D-v1)", d1, v1[0], v1[2], v1[3], v1[4], sn["v1_z_th"])
    d2 = node_figures(ROOT_D2_PF * 1e-12)
    v2 = sn["v2_rp"]
    figs("(D-v2)", d2, v2[0], v2[2], v2[3], sn["v2_z_rc"][0], sn["v2_z_th"])
    dn = node_figures(float(sn["node_75"][0]) * 1e-12)
    printed = [
        f"node C {sn['node_75'][0]} pF: R_p {dn['r_p']:.1f} (SPEC v2 {v2[0]}), Q_C "
        f"{dn['q_c']:.1f} (SPEC v2 {v2[2]}), C_p {dn['c_p_pf']:.4f} (SPEC v2 {v2[3]})",
        f"node C, against the SPEC v2 strings in units of their last digit: -Im Z_RC "
        f"{_ulps(-dn['z_rc'].imag, sn['v2_z_rc'][0]):+d}, Re Z_th "
        f"{_ulps(dn['z_th'].real, sn['v2_z_th'][0]):+d}, -Im Z_th "
        f"{_ulps(-dn['z_th'].imag, sn['v2_z_th'][1]):+d}",
        f"exact R_p = r(1 + Q^2): {dn['r_p_exact']:.1f} at the node C, {d2['r_p_exact']:.1f} at "
        f"ROOT_D2_PF (the SPEC prints the high-Q form, {v2[0]})",
    ]
    return res, printed


def node_section_lines(spec_node: dict) -> list:
    """C43 (§9 item 880): SPEC v2's node section, built only from ``node_figures`` /
    ``parse_spec_node`` output ``spec_node``.  The v2 figures are derived at the
    full-digit root the SPEC's node names (``ROOT_D2_PF`` when the node string equals
    ``f'{ROOT_D2_PF:.4f}'``; otherwise at the node as printed, so a wrong node fails
    by arithmetic, control (N-node)).  Derived values print at the SPEC's own
    decimals (a PASS reprints the SPEC's string) or <= 7 digits."""
    node = spec_node["node_75"][0]
    named = node == spec_node["node_14"][0] == f"{ROOT_D2_PF:.4f}"
    c_pf = ROOT_D2_PF if named else float(node)
    d2 = node_figures(c_pf * 1e-12)
    dn = node_figures(float(node) * 1e-12)
    v2 = spec_node["v2_rp"]
    zrc_s = spec_node["v2_z_rc"][0]
    zth_s = spec_node["v2_z_th"]
    L = [f"[{'PASS' if named else 'FAIL'}] (C43 node) C = {node} pF on all 64 gaps, r = 0.2 Ω, "
         f"64 MHz: SPEC :14 {spec_node['node_14'][0]} == :75 {node} == f'{{ROOT_D2_PF:.4f}}' "
         f"{ROOT_D2_PF:.4f} (entry 56's crossing {ROOT_D2_PF:.6f} pF, SPEC.md:17–19)"]
    for lab, v, st, unit in (("R_p", d2["r_p"], v2[0], " Ω (high-Q form 1/(ω²C²r))"),
                             ("Q_C", d2["q_c"], v2[2], ""),
                             ("C_p", d2["c_p_pf"], v2[3], " pF"),
                             ("-Im Z_RC", -d2["z_rc"].imag, zrc_s, " Ω (Z_RC = 0.2 − " + zrc_s + "j)"),
                             ("Re Z_th", d2["z_th"].real, zth_s[0], " Ω"),
                             ("-Im Z_th", -d2["z_th"].imag, zth_s[1],
                              f" Ω (Z_th = {zth_s[0]} − {zth_s[1]}j)")):
        ok = _at(v, st)
        shown = f"{v:.{len(st.split('.')[1])}f}" if ok else f"{v:.7g}"
        L.append(f"[{'PASS' if ok else 'FAIL'}] (C43 v2) {lab} derived at the root = {shown}{unit} "
                 f"vs SPEC v2 {st} ({_ulps(v, st):+d} in its last digit)")
    L.append(f"(C43 node C, printed, not asserted; the weekly's to rule) at C = {node} pF: "
             f"R_p {dn['r_p']:.1f} ({_ulps(dn['r_p'], v2[0]):+d}), Q_C {dn['q_c']:.1f} "
             f"({_ulps(dn['q_c'], v2[2]):+d}), C_p {dn['c_p_pf']:.4f} "
             f"({_ulps(dn['c_p_pf'], v2[3]):+d}), -Im Z_RC ({_ulps(-dn['z_rc'].imag, zrc_s):+d}), "
             f"Re Z_th ({_ulps(dn['z_th'].real, zth_s[0]):+d}), -Im Z_th "
             f"({_ulps(-dn['z_th'].imag, zth_s[1]):+d}), moves in units of the SPEC v2 string's last digit")
    L.append(f"(C43) R_p by r(1 + Q²) printed beside the SPEC's high-Q form: {dn['r_p_exact']:.1f} Ω "
             f"at the node C, {d2['r_p_exact']:.1f} Ω at the root; high-Q form {d2['r_p']:.1f} Ω "
             f"(SPEC v2 {v2[0]})")
    return L


# ---------------------------------------------------------------------------
# metrics.json and COMPARISON.md (pure).  Derived prints stay at <= 7
# significant digits (the private-leak hook's 8-digit floor).
# ---------------------------------------------------------------------------
#: C49 (§9 item 931): ANS-12's tracked C6 log (``xl`` entry 58), in the shape of
#: ``SOURCE_LOG_BASENAME``.  C33's module is not imported here (C33 loads this page).
SOURCE_LOG_3T_BASENAME = "20261109T080008Z_GEO-38-C6-3T.log"
SOURCE_LOG_3T = _REPO_ROOT / "docs" / "testing" / "logs" / SOURCE_LOG_3T_BASENAME
COMPARISON_MD_3T = (_REPO_ROOT / "examples" / "ansys_benchmarks"
                    / "ans12_f_clinical_3t_box_loaded_128MHz" / "COMPARISON.md")
#: The 11-09 weekly's tie clause (`ANS-11` v2.1; ANS-12 SPEC ``:150``): SAR10g within
#: 1 % of the lattice maximum.  The weekly's value, not a band set here.
TIE_REL = 0.01

#: C49: the literal strings of each case's COMPARISON.md.  ``"ans11"`` is today's text.
CASES = {
    "ans11": {
        "case": "ans:11",
        "column": "ours, degree 1 (xl entry 52) - a record, never adjudicated",
        "source_log": SOURCE_LOG_BASENAME,
        "degree2": "pending (xxl entry 62)",
        "title": "# ANS-11 — comparison (our degree-1 and degree-2 records; AED column blank)",
        "preamble_head": "F-clinical-1.5T loaded with the ASTM F2182 gel box, two-port ccw drive "
                         "(V_33, V_41) = (1, −1j) at 64 MHz, C = ",
        "node_note": "(entry 52's node)",
        "adjudication": "**Degree 1 is a record, never adjudicated. The degree-2 column is "
                        "C19's record** (`xxl` entry 62, `-n 16`, degree 2, at SPEC v2's node "
                        "185.7629 pF + 0.2 Ω, parsed from `20261114T080009Z_GEO-38-C19.log` "
                        "at the cited lines), **still un-adjudicated** (Done-when (2) is the "
                        "weekly's). ",
        "script": "11_f_clinical_1p5t_box_loaded_64MHz.py",
        "h1": "ours, degree 1 (entry 52, record)",
        "h2": "ours, degree 2 (entry 62, record)",
        "grid": "not read in entry 52 (the C6b reader was off)",
        "cu": "not printed by entry 52 (per-column partitions only, `:12282–12284`)",
        "entry": "`xl` entry 52, `-n 8`",
        "not_claimed": "- Not claimed: any adjudication or AED number (Done-when (2)), the rod "
                       "rows (entry 62's call 2), a setup figure, and any reading of the "
                       "order rule's verdict for other clinical figures (the weekly's). The "
                       "moves compare entry 52 (degree 1, C = 193.5996 pF) with entry 62 "
                       "(degree 2, C = 185.7629 pF): each move carries the node change as "
                       "well as the order change.",
        "tie": False,
    },
    "ans12": {
        "case": "ans:12",
        "column": "ours, degree 1 (xl entry 58) - a record, never adjudicated",
        "source_log": SOURCE_LOG_3T_BASENAME,
        "degree2": "pending (entry 66)",
        "title": "# ANS-12 — comparison (our degree-1 record; degree-2 and AED columns blank)",
        "preamble_head": "F-clinical-3T loaded with the ASTM F2182 gel box, two-port ccw drive "
                         "(V_33, V_41) = (1, −1j) at 128 MHz, C = ",
        "node_note": "(entry 57's degree-1 root; not the node to enter, SPEC v2)",
        "adjudication": "**Degree 1 is a record, never adjudicated; the "
                        "adjudication column is degree 2** (entry 66, pending; its node is "
                        "SPEC v2's). ",
        "script": "11_f_clinical_1p5t_box_loaded_64MHz.py --case ans12",
        "h1": "ours, degree 1 (entry 58, record)",
        "h2": "ours, degree 2 (entry 66, pending)",
        "grid": "not read in entry 58 (the C6b reader was off)",
        "cu": "not printed by entry 58 (per-column partitions only, `:12229–12231`)",
        "entry": "`xl` entry 58, `-n 8`",
        "not_claimed": "- Not claimed: the degree-2 rows, any adjudication or AED number, "
                       "the node to enter (SPEC v2).",
        "tie": True,
    },
}


def tie_set(parsed, rel=TIE_REL) -> dict:
    """Lattice balls with SAR10g within ``rel`` of the lattice maximum, and the nearest
    outsider: ``(label, SAR string, f"{100 * (s / smax - 1):.2f}")``."""
    rows = sorted(parsed["lattice"], key=lambda r: -float(r["sar10g"]))
    smax = float(rows[0]["sar10g"])
    pct = [(r["label"], r["sar10g"], f"{100 * (float(r['sar10g']) / smax - 1):.2f}") for r in rows]
    inside = [t for t, r in zip(pct, rows) if float(r["sar10g"]) >= (1 - rel) * smax]
    return {"members": inside, "outsider": pct[len(inside)]}


def _derived_strings(d: dict) -> dict:
    """``derived(...)``'s floats at the page's print formats (<= 7 significant digits)."""
    return {
        "p_acc_ccw_w": f"{d['p_acc_ccw']:.6e}",
        "p_acc_cw_w (N-sense, printed)": f"{d['p_acc_cw']:.4e}",
        "sar10g_cstar_per_p_acc_w_per_kg_per_w": f"{d['sar_cstar_per_p_acc']:.6e}",
        "p_box_per_p_acc": f"{d['p_box_per_p_acc']:.6e}",
        "lattice_mean_b1p_ccw_t": f"{d['lattice_mean_b1p']:.5e}",
        "lattice_mean_b1p_ccw_per_sqrt_p_acc_t": f"{d['lattice_mean_b1p_per_sqrt_p_acc']:.4e}",
        "lattice_b1p_ccw_per_sqrt_p_acc_t": [f"{v:.5e}" for v in d["b1p_per_sqrt_p_acc"]],
        "lattice_sar10g_per_p_acc_w_per_kg_per_w": [f"{v:.5e}" for v in d["sar_per_p_acc"]],
    }


# ---------------------------------------------------------------------------
# C37 (§9 item 983): the degree-2 column from entry 62's tracked log, the moves,
# the order rule, (D2), (X2) and the (J62) custody check against the host JSON.
# ---------------------------------------------------------------------------
#: ``xxl`` entry 62's tracked log (C19, degree 2, `-n 16`; either suffix reads).
C19_LOG_BASENAME = "20261114T080009Z_GEO-38-C19.log"
C19_LOG = _REPO_ROOT / "docs" / "testing" / "logs" / C19_LOG_BASENAME
#: The 10-24 exit's pre-registered order-rule threshold (|move| > 2 % on a per-watt
#: row ⇒ clinical figures quoted at degree 2 only).  The exit's number, not a band set here.
ORDER_RULE_REL = 0.02
ORDER_RULE_FIRES = "clinical figures quoted at degree 2 only"
ORDER_RULE_QUIET = "the 2 % order rule does not fire"
_GRID_RE = re.compile(
    _TAG + r"\(C6b\) grid n (\d+) \([^)]*\); g\* (x[+-]\d+\.\d+mm y[+-]\d+\.\d+mm dz[+-]\d+\.\d+mm) "
    r"\(([+-]?\d+\.\d+), ([+-]?\d+\.\d+), ([+-]?\d+\.\d+)\) m: SAR10g " + _F
    + r" W/kg \(ccw\), \(1,\+1j\) ball " + _F + r" W/kg")


def parse_c19_grid(path) -> dict:
    """Pure: the box-wide grid's one result line (entry 62 ``:12407``) as printed
    strings; ``ValueError`` on zero or more than one match."""
    hits = []
    for lineno, line in enumerate(_read_log(path).splitlines(), start=1):
        m = _GRID_RE.search(line)
        if m:
            hits.append({"line": lineno, **dict(zip(
                ("n", "label", "x_m", "y_m", "z_m", "sar10g", "sar10g_cw"), m.groups()))})
    if len(hits) != 1:
        raise ValueError(f"(C6b) grid result lines: {len(hits)} matches "
                         f"(at {[h['line'] for h in hits]}), want exactly 1")
    return hits[0]


def _move_rows(p: dict, d: dict) -> list:
    """``[(row, per_watt, value)]`` from one column's parsed tokens and ``derived``."""
    s, y, mn, cs, nm = p["s2"], p["y_in_ccw"], p["means"], p["c_star"], p["normalisers"]
    p_acc = float(nm["p_acc"])
    return [
        ("|S11|", False, abs(complex(s["s11"]))),
        ("|S22|", False, abs(complex(s["s22"]))),
        ("|Y_in ccw P33|", False, abs(complex(float(y["p33_re"]), float(y["p33_im"])))),
        ("|Y_in ccw P41|", False, abs(complex(float(y["p41_re"]), float(y["p41_im"])))),
        ("18-point lattice mean |B1+|_ccw", False, d["lattice_mean_b1p"]),
        ("18-point lattice mean |B1+|_ccw per sqrt(P_acc)", True,
         d["lattice_mean_b1p_per_sqrt_p_acc"]),
        ("SAR10g(c*)", False, float(cs["sar10g"])),
        ("SAR10g(c*) per P_acc", True, d["sar_cstar_per_p_acc"]),
        ("cw (1,+1j) ball at c*", False, float(cs["sar10g_cw"])),
        ("400-centroid mean |B1+|_ccw", False, float(mn["b1p_ccw"])),
        ("400-centroid mean |B1+|_ccw per sqrt(P_acc)", True,
         float(mn["b1p_ccw"]) / math.sqrt(p_acc)),
        ("P_box", False, float(nm["p_box"])),
        ("P_src", False, float(nm["p_src"])),
        ("P_acc", False, p_acc),
        ("gel fraction P_box / P_acc", True, d["p_box_per_p_acc"]),
        ("sense mean |B1+| / mean |B1-|", False, float(mn["sense"])),
    ]


def moves_and_rule(p1: dict, d1: dict, p2: dict, d2: dict) -> tuple:
    """The degree-1 → degree-2 move per row, ``f"{x2 / x1 - 1:+.3e}"`` (a ratio, never an
    AED figure), and the order rule's printed verdict over the per-watt rows."""
    moves, fired = [], []
    for (row, pw, x1), (row2, _, x2) in zip(_move_rows(p1, d1), _move_rows(p2, d2)):
        assert row == row2
        rel = x2 / x1 - 1.0
        moves.append({"row": row, "per_watt": pw, "move": f"{rel:+.3e}"})
        if pw and abs(rel) > ORDER_RULE_REL:
            fired.append(row)
    sentence = ORDER_RULE_FIRES if fired else ORDER_RULE_QUIET
    return moves, {"threshold": f"{ORDER_RULE_REL:.2f}", "rows_over": fired,
                   "per_watt_rows": [m["row"] for m in moves if m["per_watt"]],
                   "sentence": sentence}


def degree2_block(path=C19_LOG, source=None) -> dict:
    """C37: the degree-2 column from entry 62's tracked log (``parse_c6_log`` unedited,
    ``derived``, ``parse_c19_grid``), with the moves against entry 52 (``source``,
    default ``SOURCE_LOG``) and the order rule.  ``ValueError`` if the log's node is
    not ``C19_NODE_PF``."""
    p2 = parse_c6_log(path)
    if p2["node"]["c_pf"] != C19_NODE_PF:
        raise ValueError(f"entry 62 node {p2['node']['c_pf']} pF != C19_NODE_PF {C19_NODE_PF}")
    d2 = derived(p2)
    g = parse_c19_grid(path)
    p1 = parse_c6_log(SOURCE_LOG if source is None else source)
    moves, rule = moves_and_rule(p1, derived(p1), p2, d2)
    return {"source_log": Path(str(path)).name.removesuffix(".gz"), "parsed": p2,
            "derived": _derived_strings(d2), "grid": g, "moves": moves, "order_rule": rule,
            "c_star_move": f"{p1['c_star']['label']} -> {p2['c_star']['label']}"}


def check_d2(block: dict) -> list:
    """(D2): the log-internal identities on entry 62, at ``check_derived``'s tolerance for
    the same quantity (S12 / S21 strings printed, never asserted)."""
    p = block["parsed"]
    s, nm, cs = p["s2"], p["normalisers"], p["c_star"]
    s2 = np.array([[complex(s["s11"]), complex(s["s12"])], [complex(s["s21"]), complex(s["s22"])]])
    p_acc = float(nm["p_acc"])
    lat_b = np.array([float(r["b1p_ccw"]) for r in p["lattice"]])
    lat_s = np.array([float(r["sar10g"]) for r in p["lattice"]])
    res = []

    def near(name, val, tok, tol, rel=False):
        t = float(tok)
        err = abs(val - t) / abs(t) if rel else abs(val - t)
        res.append((f"(D2) {name}", err <= tol, f"{val:.6e} vs printed {tok}: "
                                               f"{'rel ' if rel else ''}{err:.3e} (<= {tol:.0e})"))

    near("_p_acc_spec(S2, (1,-1j), V_s, 50) vs :%d P_acc" % nm["line"],
         _p_acc_spec(s2, CCW, complex(p["node"]["v_s"]), Z0_OHM), nm["p_acc"], 1e-8, rel=True)
    near("SAR10g(c*)/P_acc vs the log's per-P_acc token", float(cs["sar10g"]) / p_acc,
         nm["sar_per_p_acc"], 5e-10)
    near("P_box/P_acc vs the log's printed ratio", float(nm["p_box"]) / p_acc,
         nm["p_box_over_p_acc"], 5e-9)
    tok = nm["b1_per_sqrt_p_acc"]
    half = 0.5 * 10.0 ** (math.floor(math.log10(abs(float(tok)))) - 6)
    near("centroid mean / sqrt(P_acc) vs the log's token at its printed digits",
         float(p["means"]["b1p_ccw"]) / math.sqrt(p_acc), tok, half)
    am = p["lattice"][int(np.argmax(lat_s))]["label"]
    res.append(("(D2) argmax(lattice SAR10g) == c* label", am == cs["label"],
                f"argmax {am} vs c* {cs['label']}"))
    lat = {r["label"]: r for r in p["lattice"]}.get(cs["label"])
    res.append(("(D2) c* SAR string == its lattice line's string",
                lat is not None and lat["sar10g"] == cs["sar10g"],
                f"c* {cs['sar10g']} vs lattice {lat and lat['sar10g']}"))
    res.append(("(D2) lattice mean from 18 lines == the block's derived string",
                f"{float(lat_b.mean()):.5e}" == block["derived"]["lattice_mean_b1p_ccw_t"],
                f"{float(lat_b.mean()):.5e} (arithmetic on :{p['lattice'][0]['line']}–"
                f"{p['lattice'][-1]['line']})"))
    return res


def x2_misses(text: str, block: dict) -> list:
    """(X2): every degree-2 cell of ``text`` carries entry 62's token string(s) for its row,
    with its ``:line``; the copper / ESR cell is empty; 18 degree-2 lattice rows."""
    c = comparison_cells(text)["main"]
    p, dv, g = block["parsed"], block["derived"], block["grid"]
    s, y, mn, cs, nm = (p[k] for k in ("s2", "y_in_ccw", "means", "c_star", "normalisers"))
    want = {
        Q_S11: [s["s11"], f":{s['line']}"], Q_S12: [s["s12"]], Q_S22: [s["s22"]],
        Q_Y33: [y["p33_re"], y["p33_im"], f":{y['line']}"], Q_Y41: [y["p41_re"], y["p41_im"]],
        Q_LAT: [dv["lattice_mean_b1p_ccw_t"], f":{p['lattice'][0]['line']}"],
        Q_LAT_P: [dv["lattice_mean_b1p_ccw_per_sqrt_p_acc_t"]],
        Q_CSTAR: [cs["label"], cs["r_m"], cs["dz_m"], cs["az_deg"], f":{cs['line']}"],
        Q_SAR: [cs["sar10g"], nm["sar_per_p_acc"], f":{nm['line']}"], Q_CW: [cs["sar10g_cw"]],
        Q_GRID: [g["label"], g["sar10g"], f":{g['line']}"],
        Q_MEAN: [mn["b1p_ccw"], f"{mn['n_used']} / {mn['n_total']}", f":{mn['line']}"],
        Q_MEAN_P: [nm["b1_per_sqrt_p_acc"]],
        Q_P: [nm["p_box"], nm["p_src"], nm["p_acc"]], Q_GEL: [nm["p_box_over_p_acc"]],
        Q_SENSE: [mn["sense"], f":{mn['line']}"],
    }
    miss = [(q, t) for q, toks in want.items() for t in toks if t not in c.get(q, ("", ""))[1]]
    if c.get(Q_CU, ("", None))[1] != "":
        miss.append((Q_CU, "non-empty"))
    rows = re.findall(r"^\| (x[+-]\d+mm dz[+-]\d+mm) \| 2 \| `:(\d+)` \| (\S+) \| \S+ \| (\S+) \| \S+ \|$",
                      text, re.M)
    if rows != [(r["label"], str(r["line"]), r["b1p_ccw"], r["sar10g"]) for r in p["lattice"]]:
        miss.append(("degree-2 lattice rows", len(rows)))
    return miss


def custody_c19(block: dict, c19) -> tuple | None:
    """(J62): each raw field the C6 JSON carries and the log prints, formatted in the log
    token's own print format and compared as strings.  ``(n_fields, mismatches,
    print_form_only)``; ``None`` without the JSON.  A value-equal, print-form-only
    difference is listed, not red (10:30 11-12 ruling)."""
    rec = (c19 or {}).get("c6")
    if rec is None:
        return None
    p, g = block["parsed"], block["grid"]
    pairs = []
    for (i, j), k in zip(((0, 0), (0, 1), (1, 0), (1, 1)), ("s11", "s12", "s21", "s22")):
        pairs.append((f"s2.{k}", f"{complex(*rec['s2'][i][j]):.10e}", p["s2"][k], complex))
    for sense in ("ccw", "cw"):
        for idx, port in enumerate(("p33", "p41")):
            v = complex(*rec["y_in"][sense][idx])
            tok = p[f"y_in_{sense}"]
            pairs.append((f"y_in.{sense}.{port}", f"{v.real:+.9e} {v.imag:+.9e}",
                          f"{tok[f'{port}_re']} {tok[f'{port}_im']}", None))
    lat = {r["label"]: r for r in rec["lattice"]}
    for r in p["lattice"]:
        jr = lat.get(r["label"], {})
        for key, jkey in (("b1p_ccw", "b_plus_ccw"), ("sar10g", "sar10g")):
            jv = jr.get(jkey)
            pairs.append((f"lattice.{r['label']}.{key}",
                          None if jv is None else f"{jv:.9e}", r[key], float))
    cs = p["c_star"]
    pairs += [("c_star.label", str(rec["star_label"]), cs["label"], None),
              ("c_star.sar10g", f"{rec['sar_star']:.9e}", cs["sar10g"], float),
              ("c_star.sar10g_cw", f"{rec['sar_star_cw']:.9e}", cs["sar10g_cw"], float),
              ("mean_b_plus.ccw", f"{rec['mean_b_plus']['ccw']:.9e}", p["means"]["b1p_ccw"], float)]
    for k in ("p_box", "p_src", "p_acc"):
        pairs.append((k, f"{rec[k]:.9e}", p["normalisers"][k], float))
    jg = rec.get("grid") or {}
    pairs += [("grid.star_label", jg.get("star_label"), g["label"], None),
              ("grid.sar10g", None if "sar_star" not in jg else f"{jg['sar_star']:.9e}",
               g["sar10g"], float)]
    mism, pform = [], []
    for name, js, tok, kind in pairs:
        if js == tok:
            continue
        if js is not None and kind is not None and kind(js) == kind(tok):
            pform.append((name, js, tok))
        else:
            mism.append((name, js, tok))
    return len(pairs), mism, pform


def build_metrics(parsed: dict, d: dict, case: str = "ans11") -> dict:
    cs = CASES[case]
    return {
        "case": cs["case"],
        "column": cs["column"],
        "source_log": cs["source_log"],
        "parsed": parsed,
        "derived": _derived_strings(d),
        "degree2_c19": degree2_block() if case == "ans11" else cs["degree2"],
        "aed": "never written here (CLAUDE.md hard rule)",
    }


def _imsign(s: str) -> str:
    return "−" if s.startswith("-") else "+"


def _column_cells(p: dict, dv: dict, grid_cell: str, cu_cell: str) -> list:
    """One column's ``[(Quantity, cell)]`` from parsed tokens and derived strings (the
    degree-1 cell formats; the degree-2 column reuses them row by row)."""
    s, y, mn, cs, nm = p["s2"], p["y_in_ccw"], p["means"], p["c_star"], p["normalisers"]
    return [
        (Q_S11, f"{s['s11']} (`:{s['line']}`)"),
        (Q_S12, f"{s['s12']}"),
        (Q_S22, f"{s['s22']}"),
        (Q_Y33, f"{y['p33_re']} {y['p33_im']}j (Im {_imsign(y['p33_im'])}; `:{y['line']}`)"),
        (Q_Y41, f"{y['p41_re']} {y['p41_im']}j (Im {_imsign(y['p41_im'])})"),
        (Q_LAT, f"{dv['lattice_mean_b1p_ccw_t']} (arithmetic on `:{p['lattice'][0]['line']}–"
                f"{p['lattice'][-1]['line']}`)"),
        (Q_LAT_P, f"{dv['lattice_mean_b1p_ccw_per_sqrt_p_acc_t']}"),
        (Q_CSTAR, f"{cs['label']} (r {cs['r_m']} m, dz {cs['dz_m']} m, az {cs['az_deg']}°; "
                  f"`:{cs['line']}`)"),
        (Q_SAR, f"{cs['sar10g']}; {nm['sar_per_p_acc']} (`:{nm['line']}`)"),
        (Q_CW, f"{cs['sar10g_cw']}"),
        (Q_GRID, grid_cell),
        (Q_MEAN, f"{mn['b1p_ccw']} over {mn['n_used']} / {mn['n_total']} tag-3 centroids "
                 f"(`:{mn['line']}`)"),
        (Q_MEAN_P, f"{nm['b1_per_sqrt_p_acc']}"),
        (Q_P, f"{nm['p_box']} / {nm['p_src']} / {nm['p_acc']}"),
        (Q_GEL, f"{nm['p_box_over_p_acc']}"),
        (Q_CU, cu_cell),
        (Q_SENSE, f"{mn['sense']} (`:{mn['line']}`)"),
    ]


def render_comparison(m: dict, c19=None, case: str = "ans11") -> str:
    k = CASES[case]
    p, dv = m["parsed"], m["derived"]
    log = m["source_log"]
    blk = m.get("degree2_c19") if isinstance(m.get("degree2_c19"), dict) else None
    L = [
        k["title"],
        "",
        k["preamble_head"]
        + f"{p['node']['c_pf']} pF + 0.2 Ω on all 64 gaps {k['node_note']}, raw values "
        "at V_th = 1 V per feed. " + k["adjudication"] + "Every "
        f"value is the printed token of `{log}` at the cited line (parsed by "
        f"`{k['script']}`; no solve). AED values are never "
        "written here (PROJECT_PLAN §5.4). Rows follow [SPEC.md](SPEC.md).",
        "",
        f"| Row | Quantity | {k['h1']} | {k['h2']} | AED First Order |",
        "|---|---|---|---|---|",
    ]
    d1 = _column_cells(p, dv, k["grid"], k["cu"])
    if c19 is not None or blk is None:
        d2 = degree2_cells(c19)  # C36's stub path (or no degree-2 record)
    else:
        g = blk["grid"]
        d2 = dict(_column_cells(blk["parsed"], blk["derived"],
                                f"{g['label']}: {g['sar10g']} (`:{g['line']}`)", ""))
    for q, c1 in d1:
        c2 = d2.get(q, "")
        L.append(f"| {q[0]} | {q[1]} | {c1} |{f' {c2} ' if c2 else ' '}| |")
    L += [
        "",
        "## Row 2 / 3 lattice (ours, degree 1; raw at V_th = 1 V, and per P_acc)",
        "",
        "| Point | log line | `\\|B₁⁺\\|`_ccw [T] | per √W P_acc [T/√W] | SAR10g [W/kg] "
        "| per W P_acc [W/kg/W] |",
        "|---|---|---|---|---|---|",
    ]
    for r, b, q in zip(p["lattice"], dv["lattice_b1p_ccw_per_sqrt_p_acc_t"],
                       dv["lattice_sar10g_per_p_acc_w_per_kg_per_w"]):
        L.append(f"| {r['label']} | `:{r['line']}` | {r['b1p_ccw']} | {b} | {r['sar10g']} | {q} |")
    if k["tie"]:
        t = tie_set(p)
        L += ["", "## Row 3 tie set", "",
              f"Lattice balls with SAR10g within {100 * TIE_REL:.0f} % of the lattice maximum "
              "(the 11-09 weekly's tie clause, `ANS-11` v2.1; SPEC `:150`), and the nearest "
              "outsider:", "",
              "| Point | SAR10g [W/kg] | vs the maximum [%] | role |", "|---|---|---|---|"]
        for i, (lab, sar, pc) in enumerate(t["members"]):
            L.append(f"| {lab} | {sar} | {pc} | {'maximum' if i == 0 else 'member'} |")
        o = t["outsider"]
        L.append(f"| {o[0]} | {o[1]} | {o[2]} | nearest outsider |")
    if blk is not None:
        p2, dv2 = blk["parsed"], blk["derived"]
        L += [
            "",
            "## Row 2 / 3 lattice (ours, degree 2; raw at V_th = 1 V, and per P_acc)",
            "",
            "| Point | degree | log line | `\\|B₁⁺\\|`_ccw [T] | per √W P_acc [T/√W] "
            "| SAR10g [W/kg] | per W P_acc [W/kg/W] |",
            "|---|---|---|---|---|---|---|",
        ]
        for r, b, q in zip(p2["lattice"], dv2["lattice_b1p_ccw_per_sqrt_p_acc_t"],
                           dv2["lattice_sar10g_per_p_acc_w_per_kg_per_w"]):
            L.append(f"| {r['label']} | 2 | `:{r['line']}` | {r['b1p_ccw']} | {b} | "
                     f"{r['sar10g']} | {q} |")
        s2 = p2["s2"]
        rule = blk["order_rule"]
        L += [
            "",
            "## Degree 1 → degree 2 (records)",
            "",
            "Move = degree-2 value / degree-1 value − 1, from the two tracked logs' tokens "
            "(a ratio of our own numbers; no AED figure). Per-watt rows are the order "
            "rule's.",
            "",
            "| Quantity | per-watt row | move |",
            "|---|---|---|",
        ]
        for mv in blk["moves"]:
            L.append(f"| {mv['row']} | {'yes' if mv['per_watt'] else 'no'} | {mv['move']} |")
        over = ", ".join(rule["rows_over"]) or "none"
        L += [
            "",
            f"- Order rule (the 10-24 exit's, |move| > {100 * float(rule['threshold']):.0f} % on "
            f"a per-watt row; a printed verdict, not an assertion): rows over: {over}. "
            f"**{rule['sentence'][0].upper() + rule['sentence'][1:]}.**",
            f"- c\\*: {blk['c_star_move'].replace('->', '→')} (degree 1 → degree 2). The "
            f"box-wide g\\* is degree-2 only (entry 52 ran no grid): "
            f"{blk['grid']['label']} ({blk['grid']['x_m']}, {blk['grid']['y_m']}, "
            f"{blk['grid']['z_m']}) m, `:{blk['grid']['line']}`; no claim about a maximum "
            "outside the grid.",
            f"- Degree-2 S₁₂ and S₂₁ strings (`:{s2['line']}`): {s2['s12']} / {s2['s21']} — "
            "they differ in the last printed digit; reciprocity is the window's asserted "
            "(R1) on that line, not a string identity.",
        ]
    L += [
        "",
        "## Provenance",
        "",
        f"- Source: `docs/testing/logs/{log}` ({k['entry']}; price line "
        f"`:{p['price']['line']}`: build {p['price']['build_s']} s, window "
        f"{p['price']['window_s']} s).",
    ]
    if blk is not None:
        pr = blk["parsed"]["price"]
        L.append(f"- Degree-2 source: `docs/testing/logs/{blk['source_log']}` (`xxl` entry 62, "
                 f"`-n 16`; price line `:{pr['line']}`: build {pr['build_s']} s, window "
                 f"{pr['window_s']} s). Entry 62's host JSON is a custody check only (the "
                 "page's (J62)); this file never reads it.")
    L += [
        f"- P_acc re-derived from the printed S and V_s through C6's `_p_acc_spec`: "
        f"{dv['p_acc_ccw_w']} W (ccw); the cw weights give {dv['p_acc_cw_w (N-sense, printed)']}"
        " W (N-sense, printed only).",
        k["not_claimed"],
        "",
    ]
    return "\n".join(L)


def main_ans12(rank: int) -> int:
    """C49 (§9 item 931): parse entry 58, render the ``"ans12"`` case, write ANS-12's
    ``COMPARISON.md`` only (no ANS-12 ``metrics.json``) and print.  Degree 2 blank
    (``c19=None``: the 3 T node is SPEC v2's, entry 66)."""
    parsed = parse_c6_log(SOURCE_LOG_3T)
    m = build_metrics(parsed, derived(parsed), case="ans12")
    text = render_comparison(m, None, case="ans12")
    if rank == 0:
        t = tie_set(parsed)
        print(f"[ANS-12 C49] source {SOURCE_LOG_3T_BASENAME}; lattice {len(parsed['lattice'])} "
              f"lines :{parsed['lattice'][0]['line']}–{parsed['lattice'][-1]['line']}")
        print(f"[ANS-12 C49] tie set (TIE_REL {TIE_REL}): "
              + "; ".join(f"{a} {b} ({c})" for a, b, c in t["members"])
              + f"; nearest outsider {t['outsider'][0]} {t['outsider'][1]} ({t['outsider'][2]})")
        COMPARISON_MD_3T.write_text(text)
        print(f"[ANS-12 C49] wrote {COMPARISON_MD_3T.relative_to(_REPO_ROOT)} "
              f"({len(text.encode())} B)")
    return 0


# ---------------------------------------------------------------------------
# C54 (§9 item 959): ANS-13's thermal comparison.  A separate table beside ``CASES``
# (``render_comparison`` / ``build_metrics`` read EM fields); the ``"ans11"`` /
# ``"ans12"`` paths are untouched.  C46's reader, P61's probes and RHO_C are imported
# lazily inside the functions below, so the import closure C33 / ``mri:41`` load is
# unchanged.  Never imports the ``ans:13`` page (it loads this one) or C33's module.
# ---------------------------------------------------------------------------
COMPARISON_MD_13 = (_REPO_ROOT / "examples" / "ansys_benchmarks"
                    / "ans13_f_clinical_box_heating_64_128MHz" / "COMPARISON.md")
THERMAL_CASES = {
    "ans13": {
        "case": "ans:13",
        "title": "# ANS-13 — comparison (our degree-1 records; degree-2 and AED columns blank)",
        "h59": "ours, degree 1, 1.5 T (entry 59, record)",
        "h60": "ours, degree 1, 3 T (entry 60, record)",
        "h2": "ours, degree 2 (entries 69 / 70, pending)",
        "aed": "AED First Order + IcepakFEA",
        "script": "11_f_clinical_1p5t_box_loaded_64MHz.py --case ans13",
        "row5": "P87 prints it",
    },
}
_GT13 = r"\[GEO-36 I32\] \[gel twin (\w+)\] "
_RX13_VSUB = re.compile(_GT13 + r"sense (\(1,[+-]1j\)): .*?V_sub (\S+) m\^3")
_RX13_B = re.compile(_GT13 + r"\(B\) int_sub q (\S+) vs G\*P")
_RX13_CN = re.compile(_GT13 + r"CN N=(\d+) to t\* .*?\(E'\) stored (\S+) \+ surface (\S+) W: "
                      r".*?\(S2\) (\S+) \(band")
_RX13_ELAPSED = re.compile(r"^- Elapsed \(s\): (\d+)$")


def thermal_extras(path, block) -> dict:
    """Pure: inside one C46 block, the lines C46 does not parse — per sense the
    ``V_sub``, (B) and CN / (E′) / (S2) tokens — plus the footer's elapsed seconds."""
    from tests.validation.test_geo38_c39_thermal_reader import SENSES
    from tests.validation.test_geo38_c39_thermal_reader import _read_log as read13
    lines = read13(path).splitlines()
    per = {s: {} for s in SENSES}
    for i in range(block["head"]["line"], block["end"] + 1):
        ln = lines[i - 1]
        for key, rx in (("vsub", _RX13_VSUB), ("b", _RX13_B), ("cn", _RX13_CN)):
            m = rx.search(ln)
            if not m:
                continue
            sn = m.group(1)
            if sn not in per or key in per[sn]:
                raise ValueError(f"{key} line for {sn!r} at :{i} unknown or repeated")
            if key == "vsub":
                per[sn][key] = {"v_sub": m.group(3), "line": i}
            elif key == "b":
                per[sn][key] = {"int_q": m.group(2), "line": i}
            else:
                per[sn][key] = {"n": m.group(2), "stored": m.group(3), "surface": m.group(4),
                                "s2": m.group(5), "line": i}
    for sn in SENSES:
        miss = [k for k in ("vsub", "b", "cn") if k not in per[sn]]
        if miss:
            raise ValueError(f"block :{block['head']['line']}, sense {sn}: lacks {miss}")
    el = [m.group(1) for ln in lines if (m := _RX13_ELAPSED.search(ln))]
    if len(el) != 1:
        raise ValueError(f"Elapsed footer: {len(el)} matches")
    return {"senses": per, "elapsed_s": el[0]}


def read_ans13(log59=None, log60=None) -> tuple:
    """Parse entries 59 / 60 through C46's ``parse_heat_log`` (tracked logs only).
    Returns ``(blocks, extra)`` keyed ``"59"`` / ``"60"``."""
    from tests.validation.test_geo38_c39_thermal_reader import (
        ENTRY59_GLOB, ENTRY60_GLOB, find_log, parse_heat_log)
    paths = {}
    for k, glob, given in (("59", ENTRY59_GLOB, log59), ("60", ENTRY60_GLOB, log60)):
        if given is None:
            hit = find_log(glob)
            if hit is None or hit[1] != "tracked":
                raise ValueError(f"entry {k}'s log: {hit!r} (must be tracked)")
            given = hit[0]
        paths[k] = Path(given)
    blocks = {k: parse_heat_log(p) for k, p in paths.items()}
    extra = {k: thermal_extras(paths[k], blocks[k]) for k in paths}
    return blocks, extra


def thermal_cells(blocks, extra) -> list:
    """``[(row, quantity, {"59": cell, "60": cell})]`` — the degree-1 cells, each the
    log's raw token(s) with their ``:line``; derived cells (c\\* quotients at ``.6e``,
    the row-3 closed form at ``.4e``) recomputed from the parsed tokens."""
    from tests.validation.test_mat21_box_closed_form import RHO_C
    row5 = THERMAL_CASES["ans13"]["row5"]

    def per(k):
        s, x = blocks[k]["senses"], extra[k]["senses"]
        ccw, cw = s["ccw"], s["cw"]
        hc, fc = ccw["p61"]["hot"], ccw["p61"]["src"]
        hw, fw = cw["p61"]["hot"], cw["p61"]["src"]
        bl = ccw["i32"]["ball"]
        out = [(1, "Total power deposited in the box by the mapped source [W] (ccw (B) "
                   "int_sub q, per deposited watt)",
                f"{x['ccw']['b']['int_q']} (`:{x['ccw']['b']['line']}`)")]
        for lab, p in ccw["i32"]["probes"].items():
            out.append((2, f"θ(t\\*) at {lab} per deposited watt [K] (ccw)",
                        f"{p['theta']} (`:{p['line']}`)"))
        q = float(hc["theta_point"]) / float(fc["factor"])
        out.append((2, "θ(t\\*) at c\\* per deposited watt [K] (ccw; per P_acc / factor; "
                       "1 g ball beside)",
                    f"{q:.6e} (`:{hc['line']}` / `:{fc['line']}`); 1 g ball {bl['theta']} "
                    f"(`:{bl['line']}`)"))
        qw = float(hw["theta_point"]) / float(fw["factor"])
        out.append((2, "cw (1, +1j): θ(t\\*) at c\\* per deposited watt [K] (per P_acc / "
                       "factor); cw factor",
                    f"{qw:.6e} (`:{hw['line']}` / `:{fw['line']}`); factor {fw['factor']}"))
        v = x["ccw"]["vsub"]
        cn = x["ccw"]["cn"]
        out.append((3, "Initial volume-mean slope [K/s]: closed form 1 W / (ρc·V_sub) (the logs "
                       "print (S2) = abs(slope / closed form − 1), not the slope)",
                    f"{1.0 / (RHO_C * float(v['v_sub'])):.4e} (V_sub {v['v_sub']} m³, "
                    f"`:{v['line']}`); (S2) {cn['s2']} (`:{cn['line']}`)"))
        out.append((4, "Surface heat at t\\* per deposited watt [W] (ccw (E′))",
                    f"{cn['surface']} (`:{cn['line']}`)"))
        out.append((5, "Volume-mean rise ⟨θ⟩(t\\*) [K]", row5))
        out.append((6, "P_box / P_acc (the factor), ccw",
                    f"{fc['factor']} (`:{fc['line']}`)"))
        return out

    a, b = per("59"), per("60")
    if [r[:2] for r in a] != [r[:2] for r in b]:
        raise ValueError("entries 59 / 60 give different row sets")
    return [(r, qn, {"59": c59, "60": c60}) for (r, qn, c59), (_, _, c60) in zip(a, b)]


def render_comparison_thermal(blocks, extra) -> str:
    """ANS-13's ``COMPARISON.md`` text: two degree-1 columns, degree 2 and AED blank."""
    k = THERMAL_CASES["ans13"]
    L = [
        k["title"],
        "",
        "F-clinical box heating at 1.5 T (64 MHz) and 3 T (128 MHz): the gel box driven by "
        "the coil's ccw (1, −1j) SAR, every value **per deposited watt** (SPEC item 4) "
        "unless the row says otherwise. **Degree 1 is a record, never adjudicated; the "
        "adjudication columns are degree 2** (entries 69 / 70, PENDING on P87). Every value "
        "is the printed token of entry 59's or entry 60's log at the cited line (parsed by "
        f"`{k['script']}` through C46's reader; no solve). AED values are never written "
        "here (PROJECT_PLAN §5.4). Rows follow [SPEC.md](SPEC.md).",
        "",
        f"| Row | Quantity | {k['h59']} | {k['h60']} | {k['h2']} | {k['aed']} |",
        "|---|---|---|---|---|---|",
    ]
    for r, qn, c in thermal_cells(blocks, extra):
        L.append(f"| {r} | {qn} | {c['59']} | {c['60']} | | |")
    L += ["", "## Provenance", ""]
    for key, lab in (("59", "1.5 T"), ("60", "3 T")):
        bk, x = blocks[key], extra[key]
        L.append(f"- Entry {key} ({lab}): `docs/testing/logs/{bk['log']}` (node "
                 f"{bk['head']['node']} pF, `-n {bk['head']['ranks']}`, "
                 f"{x['elapsed_s']} s; CN N = {x['senses']['ccw']['cn']['n']}; P61 block "
                 f"`:{bk['head']['line']}–{bk['end']}`).")
    L += [
        "- Not claimed: any degree-2 value or reader (entries 69 / 70, P87); row 3's slope "
        "value (only (S2) is printed); rows 5 and 7 as rows; an adjudication or any AED "
        "number; a SPEC edit.",
        "",
    ]
    return "\n".join(L)


def comparison_cells_thermal(text: str) -> dict:
    """Pure: the thermal table ``{(row, Quantity): (1.5 T, 3 T, degree 2, AED)}`` (a
    thermal twin of ``comparison_cells``)."""
    out = {}
    for line in text.splitlines():
        if not line.startswith("| ") or line.startswith("|---"):
            continue
        cells = [c.strip() for c in line[1:-1].split(" |")]
        if len(cells) == 6 and cells[0].isdigit():
            key = (int(cells[0]), cells[1])
            if key in out:
                raise ValueError(f"thermal row {key} repeated")
            out[key] = tuple(cells[2:])
    return out


def render_ans13(log59=None, log60=None) -> str:
    """(R-cmp) hook for the ``ans:13`` page: the ``"ans13"`` render, no write."""
    blocks, extra = read_ans13(log59, log60)
    return render_comparison_thermal(blocks, extra)


def main_ans13(rank: int) -> int:
    """C54 (§9 item 959): parse entries 59 / 60, render the ``"ans13"`` case, write
    ANS-13's ``COMPARISON.md`` only (no ANS-13 ``metrics.json``) and print."""
    import tests.validation.test_geo38_f_clinical_box_heating as p61
    blocks, extra = read_ans13()
    text = render_comparison_thermal(blocks, extra)
    labels = list(blocks["59"]["senses"]["ccw"]["i32"]["probes"])
    probes = list(p61._probes())
    ok = labels == probes and list(blocks["60"]["senses"]["ccw"]["i32"]["probes"]) == probes
    if rank == 0:
        for key in ("59", "60"):
            print(f"[ANS-13 C54] entry {key}: {blocks[key]['log']} block "
                  f":{blocks[key]['head']['line']}–{blocks[key]['end']}, node "
                  f"{blocks[key]['head']['node']} pF")
        print(f"[ANS-13 C54] probes {labels} vs _probes() {probes}: {'ok' if ok else 'MISMATCH'}")
        if ok:
            COMPARISON_MD_13.write_text(text)
            print(f"[ANS-13 C54] wrote {COMPARISON_MD_13.relative_to(_REPO_ROOT)} "
                  f"({len(text.encode())} B)")
    return 0 if ok else 1


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--case", choices=sorted(CASES) + sorted(THERMAL_CASES), default="ans11")
    args = ap.parse_args(argv)
    try:
        from mpi4py import MPI
        rank = MPI.COMM_WORLD.rank
    except ImportError:  # pragma: no cover
        rank = 0
    if args.case == "ans12":
        return main_ans12(rank)
    if args.case == "ans13":
        return main_ans13(rank)
    parsed = parse_c6_log(SOURCE_LOG)
    spec = parse_spec(SPEC_PATH)
    d = derived(parsed)
    k_res = check_k(spec, COMPARISON_MD.read_text())
    spec_node = parse_spec_node(SPEC_PATH)
    d_res, node_printed = check_d(spec_node, parsed)
    node_lines = node_section_lines(spec_node)
    c19 = read_c19()
    results = check_extraction(parsed, spec) + check_derived(parsed, spec, d) + k_res + d_res
    metrics = build_metrics(parsed, d)
    blk = metrics["degree2_c19"]
    text = render_comparison(metrics, None)  # the tracked-log render; never read_c19()'s
    results += check_d2(blk)
    xm = x2_misses(text, blk)
    results.append(("(X2) every degree-2 cell carries its entry-62 token(s) and :line",
                    not xm, f"misses {xm}"))
    cust = custody_c19(blk, c19)
    if cust is not None:
        results.append(("(J62) host C6 JSON formatted to the log's print == the log token",
                        not cust[1], f"n_fields {cust[0]}, mismatches {cust[1]}, "
                                     f"print_form_only {cust[2]}"))
    if rank == 0:
        print(f"[ANS-11 page] source {SOURCE_LOG_BASENAME}; lattice {len(parsed['lattice'])} lines "
              f":{parsed['lattice'][0]['line']}–{parsed['lattice'][-1]['line']}")
        for name, ok, detail in results:
            print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}")
        print("[ANS-11 page] SPEC v2's node section (C43):")
        for line in node_lines:
            print(line)
        for line in node_printed:
            print(f"[ANS-11 page] (printed, not asserted; rule (e)) {line}")
        print("[ANS-11 page] read_c19: " + "; ".join(
            f"{path} {'present' if c19[k] is not None else 'absent'}"
            for k, path in zip(("c6", "i36"), c19["paths"])))
        print(f"[ANS-11 page] (N-sense, printed) P_acc cw {d['p_acc_cw']:.4e} W vs ccw "
              f"{d['p_acc_ccw']:.4e} W, ratio {d['p_acc_cw'] / d['p_acc_ccw']:.4f}")
        skip = os.environ.get(SKIP_PRIVATE_ENV) == "1"
        print(f"[ANS-11 page] private side: none on this page (no AED read or write); "
              f"{SKIP_PRIVATE_ENV}={'1' if skip else 'unset'}")
        p2 = blk["parsed"]
        print(f"[ANS-11 C37] degree-2 source {blk['source_log']}; node {p2['node']['c_pf']} pF "
              f"`:{p2['node']['line']}`; lattice :{p2['lattice'][0]['line']}–"
              f"{p2['lattice'][-1]['line']}; grid `:{blk['grid']['line']}`")
        for q, c2 in _column_cells(p2, blk["derived"], f"{blk['grid']['label']}: "
                                   f"{blk['grid']['sar10g']} (`:{blk['grid']['line']}`)", ""):
            print(f"[ANS-11 C37] degree-2 row {q[0]} | {q[1]} | {c2}")
        for mv in blk["moves"]:
            print(f"[ANS-11 C37] move {mv['row']}: {mv['move']}"
                  f"{' (per-watt)' if mv['per_watt'] else ''}")
        rule = blk["order_rule"]
        print(f"[ANS-11 C37] order rule (|move| > {rule['threshold']} on {rule['per_watt_rows']}): "
              f"rows over {rule['rows_over']} => {rule['sentence']} (printed verdict, not asserted)")
        print(f"[ANS-11 C37] c* {blk['c_star_move']}; g* {blk['grid']['label']} "
              f"({blk['grid']['x_m']}, {blk['grid']['y_m']}, {blk['grid']['z_m']}) m beside c*")
        print(f"[ANS-11 C37] (printed) degree-2 S12 {p2['s2']['s12']} vs S21 {p2['s2']['s21']} "
              f"(string equality not asserted; (R1) in-window)")
        if cust is None:
            print("[ANS-11 C37] (J62) host C6 JSON absent; custody not run (the tracked-log "
                  "route is asserted regardless)")
        else:
            print(f"[ANS-11 C37] (J62) n_fields {cust[0]}, mismatches {len(cust[1])}, "
                  f"print_form_only {cust[2]}")
            pts = np.asarray((c19["c6"].get("grid") or {}).get("points") or [[np.nan] * 3])
            print(f"[ANS-11 C37] (custody print, not rendered) host grid x extent "
                  f"[{1e3 * np.nanmin(pts[:, 0]):+.1f}, {1e3 * np.nanmax(pts[:, 0]):+.1f}] mm, "
                  f"z extent [{1e3 * np.nanmin(pts[:, 2]):+.1f}, {1e3 * np.nanmax(pts[:, 2]):+.1f}] "
                  f"mm; g* x {1e3 * float(blk['grid']['x_m']):+.1f} mm")
        METRICS_JSON.write_text(json.dumps(metrics, indent=1, sort_keys=True) + "\n")
        COMPARISON_MD.write_text(text)
        print(f"[ANS-11 page] wrote {METRICS_JSON.name}, {COMPARISON_MD.name}")
    failed = [n for n, ok, _ in results if not ok]
    failed += [ln for ln in node_lines if ln.startswith("[FAIL]")]
    if failed:
        if rank == 0:
            print(f"[ANS-11 page] FAILED: {failed}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
