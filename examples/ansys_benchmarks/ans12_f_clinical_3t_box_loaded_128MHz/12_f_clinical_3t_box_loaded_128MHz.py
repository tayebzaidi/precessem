"""`ANS-12` — the `ans:12` record page: our degree-1 column of the
F-clinical-3T ASTM-box-loaded SPEC, read from entry 58's tracked log (no solve).

§9 item 930 (docs/planning/chunks/ANS-12.md, "§9 item 930"; §10 C48): a **record
reader** on the `ans:11` page's precedent.  It parses the tracked C6a log
``docs/testing/logs/20261109T080008Z_GEO-38-C6-3T.log`` (``xl`` entry 58, `-n 8`,
780 s) through C33's reader and the `ans:11` page's parser, and builds no mesh.

**A degree-1 record at the degree-1 root 44.0525 pF, not the node to enter; the
adjudication column is entry 66** (degree 2, ``xxl``, not yet run).  No AED value is
read or written here.

Anchors (asserted here and in ``tests/unit/test_ans12_record_page.py``):

* (X) the 18 lattice labels equal C6a's `_lattice()` in order; the node string is
  ``ENTRY58_NODE_PF``; every SPEC-quoted token of the quantity table (parsed from the
  SPEC by regex, never retyped) equals its parsed log field; the SPEC's 36 six-digit
  lattice entries equal ``f"{float(tok):.6e}"``; the SPEC header's gate tokens occur
  on their C6a lines; every SPEC token is matched to exactly one field;
* (D) C33's ``check_side`` on entry 58, (D-lat) and (D-latP);
* (R) the lattice table rendered from the parsed tokens equals ``SPEC.md``'s byte for
  byte; the row-3 tie set (SAR10g within 1 % of the lattice maximum) and its
  percentages equal the SPEC's strings;
* (R-cmp) conditional on item 931 (C49): absent until then.

Printed, not asserted (rule (e)): (N-sense), the cw / ccw P_acc ratio; the SPEC's
"Log line" column beside the parsed line numbers (a SPEC transcription question for
the weekly, never asserted, never edited).

Run::

    ./scripts/run_examples.sh -e ans:12 -n 1 -t 180
"""

from __future__ import annotations

import math
import re
import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[3]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

import tests.validation.test_geo38_c33_field_strength_reader as c33  # noqa: E402
import tests.validation.test_geo38_f_clinical_3t_box_loaded as c21  # noqa: E402

_M = c21._modules()
_box = _M["box"]
_lattice = _box._lattice
BRIEF_LATTICE_N = _box.BRIEF_LATTICE_N
C6A_LATTICE_DZ_M = _box.C6A_LATTICE_DZ_M
FREQUENCY_3T_HZ = _M["t3"].FREQUENCY_3T_HZ

CASE_DIR = Path(__file__).resolve().parent
SPEC_PATH = CASE_DIR / "SPEC.md"
COMPARISON_MD = CASE_DIR / "COMPARISON.md"
#: The 11-09 weekly's tie clause (`ANS-11` v2.1; this SPEC `:150`): SAR10g within 1 %
#: of our lattice maximum.  The weekly's value, not a band set here.
TIE_REL = 0.01

_TAG = r"\[GEO-38 C6a\] "
_F = r"([+-]?\d+\.\d+(?:e[+-]\d+)?)"
#: The one line the `ans:11` parser does not cover (entry 58's log `:12235`).
_C21_RE = re.compile(_TAG + r"\(C21\) mean \|B1\+\|_ccw: centroid mean " + _F
                     + r" T over (\d+) tag-3 centroids; lattice mean " + _F
                     + r" T over the 18-point lattice")


def _page():
    """The `ans:11` page module, through C33's loader (collective: call on every rank)."""
    return c33._load_page()


# ---------------------------------------------------------------------------
# Readers (pure).
# ---------------------------------------------------------------------------
def entry58_log():
    hit = c33.find_log(c33.ENTRY58_GLOB)
    if hit is None or hit[1] != "tracked":
        raise ValueError(f"entry 58's log: {hit!r} (must be tracked)")
    return hit[0]


def parse_c21_line(path) -> dict:
    lines = _page()._read_log(path).splitlines()
    hits = [(i, m) for i, ln in enumerate(lines, start=1) if (m := _C21_RE.search(ln))]
    if len(hits) != 1:
        raise ValueError(f"(C21) ccw mean line: {len(hits)} matches")
    i, m = hits[0]
    return {"line": i, "centroid": m.group(1), "n_centroids": m.group(2),
            "lattice_mean": m.group(3)}


_CNUM = r"([+−]?\d\.\d+e[+-]\d+ [+−] \d\.\d+e[+-]\d+j)"
_NUM = r"([+−]?\d+\.\d+(?:e[+-]\d+)?)"
_CITE = r" \| `:(\d+)` \|"


def parse_spec_3t(path=SPEC_PATH) -> dict:
    """Pure: the SPEC's quantity rows (``:165–173``), lattice rows (``:180–182``),
    tie-set sentence (``:184–186``) and header gate tokens (``:158–161``), by regex,
    as strings with their SPEC line numbers (and the SPEC's cited log line)."""
    text = Path(path).read_text()
    lines = text.splitlines()
    rows = {
        "s2": (r"^\| S₁₁, S₁₂ = S₂₁, S₂₂ \| " + _CNUM + "; " + _CNUM + "; " + _CNUM + _CITE,
               ("s11", "s12", "s22", "cited")),
        "y_in_ccw": (r"^\| Y_in ccw, P33 / P41 \[S\] \| " + _CNUM + " / " + _CNUM
                     + r" \(Im of opposite signs\)" + _CITE, ("p33", "p41", "cited")),
        "y_in_cw": (r"^\| Y_in cw, P33 / P41 \[S\] \| " + _CNUM + " / " + _CNUM + _CITE,
                    ("p33", "p41", "cited")),
        "lattice_mean": (r"^\| 18-point lattice mean `\\\|B₁⁺\\\|`_ccw \(the comparison row\) \| "
                         + _NUM + " T" + _CITE, ("b1p_ccw", "cited")),
        "means": (r"^\| 400-centroid gel mean `\\\|B₁⁺\\\|`_ccw / mean `\\\|B₁⁻\\\|` \| " + _NUM
                  + " T / " + _NUM + r" T; \(Sense\) " + _NUM + "; CV " + _NUM + _CITE,
                  ("b1p_ccw", "b1m_ccw", "sense", "cv", "cited")),
        "c_star": (r"^\| c\\\* \| x " + _NUM + ", z " + _NUM + r" \(r " + _NUM + r", az (\d+)°\): "
                   r"SAR10g " + _NUM + r" W/kg \(ccw\); cw ball " + _NUM + _CITE,
                   ("x_m", "z_m", "r_m", "az_deg", "sar10g", "sar10g_cw", "cited")),
        "powers": (r"^\| P_box / P_src / P_acc \[W\] \| " + _NUM + " / " + _NUM + " / " + _NUM
                   + _CITE, ("p_box", "p_src", "p_acc", "cited")),
        "per_p_acc": (r"^\| SAR10g\(c\\\*\) per P_acc; centroid mean `\\\|B₁⁺\\\|`/√P_acc \| " + _NUM
                      + " W/kg/W; " + _NUM + r" T/√W \(lattice mean: " + _NUM
                      + r" T/√W, C33\)" + _CITE,
                      ("sar_per_p_acc", "b1_per_sqrt_p_acc", "lattice_mean_per_sqrt_p_acc",
                       "cited")),
        "ratios": (r"^\| 3 T / 1\.5 T per accepted watt \(C33 reader, record\) \| lattice "
                   r"`\\\|B₁⁺\\\|`/√P_acc " + _NUM + "; centroid " + _NUM + r"; SAR10g\(c\\\*\) "
                   + _NUM + "; P_box / P_acc " + _NUM + r" \| `([^`]+)` \|",
                   ("lattice_mean_per_sqrt_p_acc", "centroid_mean_per_sqrt_p_acc",
                    "sar_per_p_acc", "p_box_over_p_acc", "cited")),
    }
    out: dict = {"lattice": [], "lattice_block": None}
    for lineno, line in enumerate(lines, start=1):
        for key, (rx, names) in rows.items():
            m = re.search(rx, line)
            if m:
                if key in out:
                    raise ValueError(f"SPEC {key}: second match at :{lineno}")
                out[key] = {"line": lineno, **{n: (g if n == "cited" else _spec_num(g))
                                               for n, g in zip(names, m.groups())}}
        if re.match(r"^\| x \| dz ", line):
            if out["lattice_block"] is not None:
                raise ValueError(f"SPEC lattice header: second match at :{lineno}")
            out["lattice_block"] = {"line": lineno, "text": "\n".join(lines[lineno - 1:lineno + 4])}
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
    if len(out["lattice"]) != 3 or out["lattice_block"] is None:
        raise ValueError(f"SPEC lattice rows: {len(out['lattice'])} != 3")
    flat = re.sub(r"\s+", " ", text)
    lab = r"(x[+−]\d+ mm dz[+−]\d+ mm)"
    m = re.findall(r"Row 3's tie set on this column: " + lab + r" \(" + _NUM + r"\), " + lab
                   + r" \(" + _NUM + ", " + _NUM + r" %\), " + lab + r" \(" + _NUM + ", " + _NUM
                   + r" %, outside 1 %\)", flat)
    if len(m) != 1:
        raise ValueError(f"SPEC tie sentence: {len(m)} matches")
    g = m[0]
    nl = lambda s: _spec_num(s)  # noqa: E731  ("x+105 mm dz−50 mm" -> "x+105mmdz-50mm")
    out["tie"] = {"members": [(nl(g[0]), g[1], None), (nl(g[2]), _spec_num(g[3]), _spec_num(g[4]))],
                  "outsider": (nl(g[5]), _spec_num(g[6]), _spec_num(g[7]))}
    h = re.findall(r"`-n 8`, (\d+) s, `memory\.peak` (\d+\.\d+) GiB; every gate green \(\(R1\) "
                   + _NUM + r", \(P\) " + _NUM + " / " + _NUM + " with NEG " + _NUM + " / " + _NUM
                   + r", \(P_RC\) " + _NUM + r", coverage \(i\) " + _NUM + " / " + _NUM
                   + "; build " + _NUM + " s, columns " + _NUM + r" \+ " + _NUM + " s, readouts "
                   + _NUM + r" s, `:(\d+)`\)", flat)
    if len(h) != 1:
        raise ValueError(f"SPEC header gate sentence: {len(h)} matches")
    names = ("elapsed_s", "peak_gib", "r1", "p_p33", "p_p41", "neg_p33", "neg_p41", "p_rc",
             "cov_power", "cov_mass", "build_s", "column_p33_s", "column_p41_s", "readouts_s",
             "cited")
    out["header"] = dict(zip(names, h[0]))
    return out


def _spec_num(s: str) -> str:
    return s.replace("−", "-").replace(" ", "")


def gate_tokens(path) -> list:
    """``[(name, log_line, log_token)]``: each header gate token as the log prints it,
    each found by its context regex on exactly one line."""
    lines = _page()._read_log(path).splitlines()
    ctx = [
        (("r1",), _TAG + r"kinds .*\(R1\) \|\|S-S\^T\|\|_F/\|\|S\|\|_F (\S+) \(band"),
        (("p_p33", "neg_p33"), _TAG + r"\(P\) column P33: .* rel (\S+) \(band.*NEG blind (\S+) \("),
        (("p_p41", "neg_p41"), _TAG + r"\(P\) column P41: .* rel (\S+) \(band.*NEG blind (\S+) \("),
        (("p_rc",), _TAG + r"\(P_RC\) column P33 at the node: rel (\S+) \(band"),
        (("cov_power", "cov_mass"), _TAG + r"\(i\) ball radius .* miss (\S+); mass .* miss (\S+) \(band"),
        (("build_s", "column_p33_s", "column_p41_s", "readouts_s"),
         _TAG + r"price: build (\S+) s; columns (\S+), (\S+) s; readouts (\S+) s;"),
        (("elapsed_s",), r"^- Elapsed \(s\): (\d+)$"),
    ]
    out = []
    for names, rx in ctx:
        hits = [(i, m) for i, ln in enumerate(lines, start=1) if (m := re.search(rx, ln))]
        if len(hits) != 1:
            raise ValueError(f"gate line {names}: {len(hits)} matches")
        i, m = hits[0]
        out += [(n, i, g) for n, g in zip(names, m.groups())]
    peak = [i for i, ln in enumerate(lines) if ln.strip() == "[XL] memory.peak bytes:"]
    if len(peak) != 1:
        raise ValueError(f"memory.peak header: {len(peak)} matches")
    raw = lines[peak[0] + 1].strip()
    out.append(("peak_gib", peak[0] + 2, f"{int(raw) / 2 ** 30:.2f}"))
    return out


def read_all(log_path=None, spec_path=SPEC_PATH) -> dict:
    page = _page()
    log_path = entry58_log() if log_path is None else log_path
    p58 = page.parse_c6_log(log_path)
    p52 = page.parse_c6_log(page.SOURCE_LOG)
    return {"page": page, "log_path": log_path, "p58": p58, "p52": p52,
            "c21": parse_c21_line(log_path), "spec": parse_spec_3t(spec_path),
            "gates": gate_tokens(log_path), "d": page.derived(p58)}


# ---------------------------------------------------------------------------
# Checks (pure): each returns ``[(name, ok, detail)]``.
# ---------------------------------------------------------------------------
def _cstr(re_s, im_s):
    return f"{re_s}{im_s}j"


def x_pairs(ctx) -> list:
    """``[(spec_name, spec_token, parsed_token, mode)]`` — one per SPEC token."""
    p, sp, c = ctx["p58"], ctx["spec"], ctx["c21"]
    s, yc, yw, mn, cs, nm = (p["s2"], p["y_in_ccw"], p["y_in_cw"], p["means"], p["c_star"],
                             p["normalisers"])
    pt = dict(_lattice())[cs["label"]]
    sd = c33.side(p)
    r = c33.ratios(c33.side(ctx["p52"]), sd)
    e6 = "{:.6e}".format
    pairs = [("s2.s11", sp["s2"]["s11"], s["s11"], "str"),
             ("s2.s12", sp["s2"]["s12"], s["s12"], "str"),
             ("s2.s22", sp["s2"]["s22"], s["s22"], "str")]
    for key, y in (("y_in_ccw", yc), ("y_in_cw", yw)):
        for port in ("p33", "p41"):
            pairs.append((f"{key}.{port}", sp[key][port], _cstr(y[f"{port}_re"], y[f"{port}_im"]),
                          "str"))
    pairs += [
        ("lattice_mean.b1p_ccw", sp["lattice_mean"]["b1p_ccw"], c["lattice_mean"], "str"),
        ("means.b1p_ccw", sp["means"]["b1p_ccw"], mn["b1p_ccw"], "str"),
        ("means.b1m_ccw", sp["means"]["b1m_ccw"], mn["b1m_ccw"], "str"),
        ("means.sense", sp["means"]["sense"], mn["sense"], "str"),
        ("means.cv", sp["means"]["cv"], mn["cv"], "str"),
        ("c_star.sar10g", sp["c_star"]["sar10g"], cs["sar10g"], "str"),
        ("c_star.sar10g_cw", sp["c_star"]["sar10g_cw"], cs["sar10g_cw"], "str"),
        # The SPEC prints r / az at fewer decimals than the log: compared as floats.
        ("c_star.r_m", sp["c_star"]["r_m"], cs["r_m"], "float"),
        ("c_star.az_deg", sp["c_star"]["az_deg"], cs["az_deg"], "float"),
        ("c_star.x_m", sp["c_star"]["x_m"], repr(pt[0]), "point"),
        ("c_star.z_m", sp["c_star"]["z_m"], repr(pt[2]), "point"),
        ("powers.p_box", sp["powers"]["p_box"], nm["p_box"], "str"),
        ("powers.p_src", sp["powers"]["p_src"], nm["p_src"], "str"),
        ("powers.p_acc", sp["powers"]["p_acc"], nm["p_acc"], "str"),
        ("per_p_acc.sar_per_p_acc", sp["per_p_acc"]["sar_per_p_acc"], nm["sar_per_p_acc"], "str"),
        ("per_p_acc.b1_per_sqrt_p_acc", sp["per_p_acc"]["b1_per_sqrt_p_acc"],
         nm["b1_per_sqrt_p_acc"], "str"),
        ("per_p_acc.lattice_mean_per_sqrt_p_acc", sp["per_p_acc"]["lattice_mean_per_sqrt_p_acc"],
         e6(sd["lattice_mean_per_sqrt_p_acc"]), "str"),
    ]
    for k in ("lattice_mean_per_sqrt_p_acc", "centroid_mean_per_sqrt_p_acc", "sar_per_p_acc",
              "p_box_over_p_acc"):
        pairs.append((f"ratios.{k}", e6(float(sp["ratios"][k])), e6(r[k]), "str"))
    return pairs


def _pair_ok(spec_s, tok, mode):
    if mode == "str":
        return spec_s == tok
    if mode == "float":
        return float(spec_s) == float(tok)
    return math.isclose(float(spec_s), float(tok), abs_tol=1e-12)


def check_x(ctx) -> list:
    p, sp = ctx["p58"], ctx["spec"]
    res = []
    labels = [r["label"] for r in p["lattice"]]
    want = [lab for lab, _ in _lattice()]
    res.append(("(X) 18 lattice labels == _lattice() in order", labels == want,
                f"{len(labels)} labels"))
    res.append(("(X) node string == ENTRY58_NODE_PF", p["node"]["c_pf"] == c33.ENTRY58_NODE_PF,
                f"log {p['node']['c_pf']} (`:{p['node']['line']}`) vs {c33.ENTRY58_NODE_PF}"))
    pairs = x_pairs(ctx)
    for name, spec_s, tok, mode in pairs:
        res.append((f"(X) {name}", _pair_ok(spec_s, tok, mode),
                    f"SPEC {spec_s} vs parsed {tok} ({mode})"))
    by_label = {r["label"]: r for r in p["lattice"]}
    bad, n = [], 0
    for row in sp["lattice"]:
        x = float(row["x_m"])
        for iz, pair in enumerate(row["cells"]):
            lab = f"x{x * 1e3:+.0f}mm dz{C6A_LATTICE_DZ_M[iz] * 1e3:+.0f}mm"
            r = by_label.get(lab)
            for spec_s, key in zip(pair, ("b1p_ccw", "sar10g")):
                n += 1
                if r is None or spec_s != f"{float(r[key]):.6e}":
                    bad.append((row["line"], lab, spec_s, r and r[key]))
    res.append(("(X) SPEC 36 six-digit lattice entries == f'{float(tok):.6e}'",
                n == 36 and not bad, f"{n} entries, mismatches {bad[:3]}"))
    hdr = sp["header"]
    for name, line, tok in ctx["gates"]:
        res.append((f"(X) header gate {name}", hdr[name] == tok,
                    f"SPEC {hdr[name]} vs log `:{line}` {tok}"))
    # Every SPEC token matched to exactly one field, none unmatched.
    spec_fields = {f"{k}.{f}" for k in ("s2", "y_in_ccw", "y_in_cw", "lattice_mean", "means",
                                        "c_star", "powers", "per_p_acc", "ratios")
                   for f in sp[k] if f not in ("line", "cited")}
    names = [nm for nm, *_ in pairs]
    gate_names = [n for n, _, _ in ctx["gates"]]
    hdr_fields = {f for f in hdr if f != "cited"}
    ok = (sorted(names) == sorted(spec_fields) and len(set(names)) == len(names)
          and sorted(gate_names) == sorted(hdr_fields) and len(set(gate_names)) == len(gate_names))
    res.append(("(X) every SPEC token matched to exactly one field, none unmatched", ok,
                f"{len(names)} quantity tokens of {len(spec_fields)}, {len(gate_names)} gate "
                f"tokens of {len(hdr_fields)}, 36 lattice entries"))
    return res


def check_d(ctx) -> list:
    p, c, sp, d = ctx["p58"], ctx["c21"], ctx["spec"], ctx["d"]
    res = list(c33.check_side(p, c33.ENTRY58_NODE_PF, "3 T (entry 58)"))
    toks = [r["b1p_ccw"] for r in p["lattice"]]
    tgt = c["lattice_mean"]
    bound = c33._half_ulp_rel(tgt) + max(c33._half_ulp_rel(t) for t in toks)
    err = abs(d["lattice_mean_b1p"] - float(tgt)) / float(tgt)
    res.append(("(D-lat) mean of the 18 tokens vs the (C21) line's lattice mean", err <= bound,
                f"{d['lattice_mean_b1p']:.6e} vs printed {tgt} (`:{c['line']}`): rel {err:.3e} "
                f"(<= {bound:.3e})"))
    res.append(("(D-lat) the (C21) line's centroid mean == means.b1p_ccw (string)",
                c["centroid"] == p["means"]["b1p_ccw"],
                f"{c['centroid']} vs {p['means']['b1p_ccw']}"))
    tgt2 = sp["per_p_acc"]["lattice_mean_per_sqrt_p_acc"]
    bound2 = (c33._half_ulp_rel(tgt2) + max(c33._half_ulp_rel(t) for t in toks)
              + 0.5 * c33._half_ulp_rel(p["normalisers"]["p_acc"]))
    err2 = abs(d["lattice_mean_b1p_per_sqrt_p_acc"] - float(tgt2)) / float(tgt2)
    res.append(("(D-latP) lattice mean / sqrt(P_acc) vs SPEC", err2 <= bound2,
                f"{d['lattice_mean_b1p_per_sqrt_p_acc']:.6e} vs SPEC {tgt2}: rel {err2:.3e} "
                f"(<= {bound2:.3e})"))
    return res


def _signed(v, fmt):
    s = format(v, fmt)
    return s.replace("-", "−")


def render_lattice_table(parsed) -> str:
    """The SPEC's lattice table (``SPEC.md:178–182``) from the parsed tokens."""
    by_label = {r["label"]: r for r in parsed["lattice"]}
    dz = C6A_LATTICE_DZ_M
    head = ["dz " + _signed(dz[0], "+.2f")] + [_signed(z, "+.2f") for z in dz[1:]]
    L = ["| x | " + " | ".join(head) + " |", "|" + "---|" * (len(dz) + 1)]
    xs = []
    for _, pt in _lattice():
        if pt[0] not in xs:
            xs.append(pt[0])
    for x in xs:
        xl = "0" if x == 0.0 else _signed(x, "+.3f")
        cells = []
        for z in dz:
            r = by_label[f"x{x * 1e3:+.0f}mm dz{z * 1e3:+.0f}mm"]
            cells.append(f"{float(r['b1p_ccw']):.6e} / {float(r['sar10g']):.6e}")
        L.append(f"| {xl} | " + " | ".join(cells) + " |")
    return "\n".join(L)


def tie_set(parsed, rel=TIE_REL) -> dict:
    """Lattice balls with SAR10g within ``rel`` of the lattice maximum, and the nearest
    outsider: ``(label, SAR string, f"{100 * (s / smax - 1):.2f}")``."""
    rows = sorted(parsed["lattice"], key=lambda r: -float(r["sar10g"]))
    smax = float(rows[0]["sar10g"])
    pct = [(r["label"], r["sar10g"], f"{100 * (float(r['sar10g']) / smax - 1):.2f}") for r in rows]
    inside = [t for t, r in zip(pct, rows) if float(r["sar10g"]) >= (1 - rel) * smax]
    return {"members": inside, "outsider": pct[len(inside)]}


def check_r(ctx) -> list:
    p, sp = ctx["p58"], ctx["spec"]
    render = render_lattice_table(p)
    blk = sp["lattice_block"]
    res = [("(R) lattice table rendered from the log == SPEC bytes", render == blk["text"],
            f"{len(render.encode())} B vs SPEC :{blk['line']}–{blk['line'] + 4} "
            f"{len(blk['text'].encode())} B")]
    t = tie_set(p)
    lab = lambda s: s.replace(" ", "")  # noqa: E731
    got_m = [(lab(a), f"{float(b):.6e}") for a, b, _ in t["members"]]
    want_m = [(a, b) for a, b, _ in sp["tie"]["members"]]
    res.append(("(R) tie set (SAR10g within 1 % of the lattice max) == SPEC's members",
                got_m == want_m, f"{[a for a, _, _ in t['members']]} vs SPEC {[a for a, _ in want_m]}"))
    res.append(("(R) tie member percentage string == SPEC's",
                [c for _, _, c in t["members"][1:]] == [c for _, _, c in sp["tie"]["members"][1:]],
                f"{[c for _, _, c in t['members']]} vs SPEC {sp['tie']['members'][1][2]}"))
    o, so = t["outsider"], sp["tie"]["outsider"]
    res.append(("(R) nearest outsider == SPEC's (label, SAR, percentage)",
                (lab(o[0]), f"{float(o[1]):.6e}", o[2]) == so, f"{o[0]} {o[2]} vs SPEC {so}"))
    return res


def check_r_cmp(ctx):
    """(R-cmp): live only once item 931 (C49) has landed; ``(results, message)``."""
    page = ctx["page"]
    has_case = "ans12" in getattr(page, "CASES", {})
    if not COMPARISON_MD.exists() and not has_case:
        return [], "absent — item 931 (C49)"
    if not (COMPARISON_MD.exists() and has_case):
        return [("(R-cmp) COMPARISON.md and the ans:11 page's 'ans12' case both present", False,
                 f"COMPARISON.md {COMPARISON_MD.exists()}, case {has_case}")], "inconsistent"
    m = page.build_metrics(ctx["p58"], ctx["d"], case="ans12")
    render = page.render_comparison(m, None, case="ans12")
    text = COMPARISON_MD.read_text()
    return [("(R-cmp) re-render through the 'ans12' case == COMPARISON.md bytes", render == text,
             f"{len(render.encode())} B vs {len(text.encode())} B")], "live"


# ---------------------------------------------------------------------------
# C58 (§9 item 989): the entry-66 degree-2 reader, written ahead of its window.
# It finds ``xxl`` entry 66's tracked log through the ledger and parses it with the
# `ans:11` page's parsers, unedited; (K) runs the same reader on entry 62's log and
# compares C37's committed parse (`ans:11`'s tracked ``metrics.json``) by string.
# ---------------------------------------------------------------------------
XXL_LEDGER = _REPO_ROOT / "docs" / "testing" / "xxl-ledger.md"
TRACKED_LOG_DIR = _REPO_ROOT / "docs" / "testing" / "logs"
ENTRY66_CHUNK = "GEO-38-C19-3T"
ENTRY62_LOG_BASENAME = "20261114T080009Z_GEO-38-C19.log"
#: C37's (J62) custody count on entry 62 (``custody_c19``'s pair names) — a record.
D2_FIELDS_N = 53
#: The chunk column matched exactly between pipes: the 1.5 T row (``GEO-38-C19``) cannot match.
_ROW66_RE = re.compile(r"^\| \d{4}-\d{2}-\d{2} \| GEO-38-C19-3T \| `(\d{8}T\d{6}Z_GEO-38-C19-3T\.log)`")


def entry66_scan(ledger=XXL_LEDGER) -> tuple:
    """``(rows, distinct base names)`` of the ledger's ``GEO-38-C19-3T`` rows."""
    names = [m.group(1) for ln in Path(ledger).read_text().splitlines()
             if (m := _ROW66_RE.search(ln))]
    return len(names), sorted(set(names))


def entry66_log(ledger=XXL_LEDGER, log_dir=TRACKED_LOG_DIR, verbose=True):
    """Entry 66's tracked log path, or ``None`` while the ledger names no row.
    ``ValueError`` on more than one distinct base name, or a named log that is not
    under ``log_dir`` with either suffix (a held log is the review's)."""
    rows, names = entry66_scan(ledger)
    if verbose:
        print(f"[ANS-12 entry66] ledger {Path(ledger).name}: rows {rows}, distinct {len(names)}")
    if not names:
        return None
    if len(names) > 1:
        raise ValueError(f"entry 66 ledger rows name {len(names)} distinct logs: {names}")
    path = Path(log_dir) / names[0]
    if not (path.exists() or Path(f"{path}.gz").exists()):
        raise ValueError(f"entry 66 log {names[0]} not under {log_dir} (either suffix): "
                         "the page reads tracked logs only")
    return path


def _d2_triples(p: dict, g) -> list:
    """``[(name, token, line)]`` in ``custody_c19``'s pair order and names; the grid pair
    only when ``g`` is not ``None``."""
    out = [(f"s2.{k}", p["s2"][k], p["s2"]["line"]) for k in ("s11", "s12", "s21", "s22")]
    for sense in ("ccw", "cw"):
        tok = p[f"y_in_{sense}"]
        for port in ("p33", "p41"):
            out.append((f"y_in.{sense}.{port}", f"{tok[f'{port}_re']} {tok[f'{port}_im']}",
                        tok["line"]))
    for r in p["lattice"]:
        for key in ("b1p_ccw", "sar10g"):
            out.append((f"lattice.{r['label']}.{key}", r[key], r["line"]))
    cs, nm = p["c_star"], p["normalisers"]
    out += [("c_star.label", cs["label"], cs["line"]), ("c_star.sar10g", cs["sar10g"], cs["line"]),
            ("c_star.sar10g_cw", cs["sar10g_cw"], cs["line"]),
            ("mean_b_plus.ccw", p["means"]["b1p_ccw"], p["means"]["line"])]
    out += [(k, nm[k], nm["line"]) for k in ("p_box", "p_src", "p_acc")]
    if g is not None:
        out += [("grid.star_label", g["label"], g["line"]), ("grid.sar10g", g["sar10g"], g["line"])]
    return out


def degree2_fields(path) -> dict:
    """Pure: a C19-shaped log's degree-2 fields through the `ans:11` page's parsers."""
    page = _page()
    p = page.parse_c6_log(path)
    n = sum(1 for ln in page._read_log(path).splitlines() if page._GRID_RE.search(ln))
    if n > 1:
        raise ValueError(f"(C6b) grid result lines: {n}, want 0 or 1")
    g = page.parse_c19_grid(path) if n == 1 else None
    tr = _d2_triples(p, g)
    return {"log": p["log"], "node": (p["node"]["c_pf"], p["node"]["line"]),
            "fields": [(a, b) for a, b, _ in tr], "lines": [(a, c) for a, _, c in tr],
            "grid_lines": n}


def fields_from_metrics(path=None) -> tuple:
    """``(fields, block)`` from `ans:11`'s tracked ``metrics.json["degree2_c19"]``."""
    import json

    path = _page().METRICS_JSON if path is None else path
    blk = json.loads(Path(path).read_text())["degree2_c19"]
    return [(a, b) for a, b, _ in _d2_triples(blk["parsed"], blk["grid"])], blk


def check_k_d2(res: dict, metrics_fields: list, block: dict) -> list:
    """(K-n), (K) exact, (K-src): ``[(name, ok, detail)]``."""
    page = _page()
    f = res["fields"]
    out = [("(K-n) degree-2 fields on entry 62 == D2_FIELDS_N", len(f) == D2_FIELDS_N,
            f"{len(f)} fields vs {D2_FIELDS_N}")]
    want = dict(metrics_fields)
    mism = [(n, t, want.get(n)) for n, t in f if want.get(n) != t]
    extra = sorted(set(want) - {n for n, _ in f})
    out.append(("(K) every field token == metrics.json degree2_c19 string",
                not mism and not extra and len(want) == len(f),
                f"{len(f)} fields, mismatches {len(mism)} {mism[:3]}, unmatched metrics {extra}"))
    src = block.get("source_log")
    ok = (src == ENTRY62_LOG_BASENAME == page.C19_LOG_BASENAME and res["node"][0] == page.C19_NODE_PF)
    out.append(("(K-src) source_log == ENTRY62_LOG_BASENAME == C19_LOG_BASENAME, node == C19_NODE_PF",
                ok, f"source_log {src}; node {res['node'][0]} (`:{res['node'][1]}`) vs "
                    f"{page.C19_NODE_PF}"))
    return out


def render_d2(res: dict) -> list:
    """Printed lines for one parsed degree-2 log."""
    L = [f"[ANS-12 page] degree 2: {res['log']} node {res['node'][0]} pF (`:{res['node'][1]}`); "
         f"grid_lines {res['grid_lines']}; {len(res['fields'])} fields"]
    if res["grid_lines"] == 0:
        L.append("[ANS-12 page] degree 2: grid: not run")
    for (n, t), (_, ln) in zip(res["fields"], res["lines"]):
        L.append(f"[ANS-12 page]   {n} | :{ln} | {t}")
    return L


def spec_line_table(ctx) -> list:
    """``[(row, SPEC-cited line, parsed line)]`` — printed, never asserted."""
    p, sp, c = ctx["p58"], ctx["spec"], ctx["c21"]
    rows = [("S", sp["s2"]["cited"], p["s2"]["line"]),
            ("Y_in ccw", sp["y_in_ccw"]["cited"], p["y_in_ccw"]["line"]),
            ("Y_in cw", sp["y_in_cw"]["cited"], p["y_in_cw"]["line"]),
            ("lattice mean", sp["lattice_mean"]["cited"], c["line"]),
            ("centroid means", sp["means"]["cited"], p["means"]["line"]),
            ("c*", sp["c_star"]["cited"], p["c_star"]["line"]),
            ("P_box / P_src / P_acc", sp["powers"]["cited"], p["normalisers"]["line"]),
            ("per P_acc", sp["per_p_acc"]["cited"], p["normalisers"]["line"]),
            ("header price", sp["header"]["cited"], p["price"]["line"])]
    return [(r, int(a), b) for r, a, b in rows]


def all_results(ctx):
    rc, msg = check_r_cmp(ctx)
    return check_x(ctx) + check_d(ctx) + check_r(ctx) + rc, msg


def main() -> int:
    from mpi4py import MPI

    comm = MPI.COMM_WORLD
    _page()  # collective import on every rank (EX-193, deviation 2)
    ctx = read_all()
    results, cmp_msg = all_results(ctx)
    k62 = degree2_fields(TRACKED_LOG_DIR / ENTRY62_LOG_BASENAME)
    mf, mblk = fields_from_metrics()
    results = results + check_k_d2(k62, mf, mblk)
    e66 = entry66_log(verbose=comm.rank == 0)
    d66 = None if e66 is None else degree2_fields(e66)
    if comm.rank == 0:
        p, d = ctx["p58"], ctx["d"]
        print(f"[ANS-12 page] source {p['log']} ({c33.find_log(c33.ENTRY58_GLOB)[1]}); "
              f"node {p['node']['c_pf']} pF (`:{p['node']['line']}`), f = "
              f"{FREQUENCY_3T_HZ / 1e6:.0f} MHz (t3.FREQUENCY_3T_HZ); lattice "
              f"{len(p['lattice'])} lines :{p['lattice'][0]['line']}–{p['lattice'][-1]['line']}")
        print("[ANS-12 page] a degree-1 record at the degree-1 root; not the node to enter "
              "(SPEC v2); the adjudication column is entry 66")
        for name, ok, detail in results:
            print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}")
        print(f"[ANS-12 page] (R-cmp) {cmp_msg}")
        print("[ANS-12 page] lattice (ours, degree 1; raw at V_th = 1 V, and per P_acc):")
        print("[ANS-12 page]   point | log line | |B1+|_ccw [T] | per sqrt(W) P_acc | SAR10g [W/kg] "
              "| per W P_acc")
        for r, b, s in zip(p["lattice"], d["b1p_per_sqrt_p_acc"], d["sar_per_p_acc"]):
            print(f"[ANS-12 page]   {r['label']} | :{r['line']} | {r['b1p_ccw']} | {b:.5e} | "
                  f"{r['sar10g']} | {s:.5e}")
        t = tie_set(p)
        print("[ANS-12 page] row-3 tie set (SAR10g within 1 % of the lattice maximum): "
              + "; ".join(f"{a} {b} ({c} %)" for a, b, c in t["members"])
              + f"; nearest outsider {t['outsider'][0]} {t['outsider'][1]} ({t['outsider'][2]} %)")
        print("[ANS-12 page] SPEC 'Log line' vs parsed line (printed, never asserted; a SPEC "
              "transcription question for the weekly):")
        for row, a, b in spec_line_table(ctx):
            print(f"[ANS-12 page]   {row}: SPEC :{a} | parsed :{b} | "
                  f"{'same' if a == b else 'MISMATCH'}")
        print(f"[ANS-12 page] (N-sense, printed) P_acc cw {d['p_acc_cw']:.4e} W vs ccw "
              f"{d['p_acc_ccw']:.4e} W, ratio {d['p_acc_cw'] / d['p_acc_ccw']:.4f}")
        (n0, l0), (n1, l1) = k62["lines"][0], k62["lines"][-1]
        print(f"[ANS-12 page] (K) entry 62 through the degree-2 reader: {len(k62['fields'])} "
              f"fields, first {n0} (`:{l0}`), last {n1} (`:{l1}`)")
        if d66 is None:
            rows, names = entry66_scan()
            print(f"[ANS-12 page] degree 2: (Pend) entry 66 pending — xxl-ledger.md names no "
                  f"{ENTRY66_CHUNK} row (rows {rows}, distinct {len(names)})")
        else:
            for ln in render_d2(d66):
                print(ln)
        print("[ANS-12 page] writes nothing; no AED read or write")
    failed = [n for n, ok, _ in results if not ok]
    if failed:
        if comm.rank == 0:
            print(f"[ANS-12 page] FAILED: {failed}")
        return 1
    return 0


if __name__ == "__main__":
    try:
        rc = main()
    except BaseException:
        import traceback

        traceback.print_exc()
        from mpi4py import MPI

        MPI.COMM_WORLD.Abort(1)
        raise
    sys.exit(rc)
