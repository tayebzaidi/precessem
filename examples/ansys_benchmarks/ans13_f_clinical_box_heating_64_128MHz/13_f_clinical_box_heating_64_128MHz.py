"""`ANS-13` — the `ans:13` record page: our degree-1 columns of the F-clinical
box-heating SPEC (1.5 T and 3 T), read from ``xl`` entries 59 / 60's tracked logs
(no solve, no mesh).

§9 item 958 (docs/planning/chunks/ANS-13.md, "§9 item 958"; §10 C53): a **record
reader**.  It parses ``docs/testing/logs/20261110T080008Z_GEO-38-P33.log`` (entry 59,
1.5 T, node 193.5996 pF) and ``20261111T080008Z_GEO-38-P33-3T.log`` (entry 60, 3 T,
node 44.0525 pF) through C46's reader (``tests/validation/test_geo38_c39_thermal_reader.py``),
imports P61's / C31's modules for the probe set and the knob names, and never
re-implements them.

**Degree-1 records at the degree-1 roots; the adjudication columns are entries 69 / 70**
(degree 2, PENDING on §10 P87).  No AED value is read or written here.

Anchors (asserted here and in ``tests/unit/test_ans13_record_page.py``):

* (X) C46's ``check_block`` on both blocks; SPEC row 2's seven coordinates; custody
  (each log's command line carries its knob, ``both`` and its node); every token of the
  SPEC's record table (``SPEC.md:140–148``) matched to exactly one parsed field, with
  its line cite; the thermal table (``:78–86``) against the imported constants;
* (D) (D-dep) θ per P_acc / factor = I32's per-deposited-watt token (24 checks);
  (D-c*) the SPEC's per-deposited-watt c* figures; (D-ratio) C46's ``ratios`` = C39's
  22 strings and the SPEC's 4.0425 / 4.5397; (D-lever) 3.6935 / 4.1014; (D-cal)
  9.737e-06 = 1 / (ρc · V_sub);
* (R) the record table rendered from the parsed tokens equals ``SPEC.md:140–148``
  byte for byte;
* (R-cmp) conditional on item 959 (C54): absent until ANS-13's ``COMPARISON.md``
  exists.  Once it does, the `ans:11` page (loaded through C33's ``_load_page``) must
  expose ``render_ans13()`` returning the file's text; the page asserts byte identity.

Writes nothing.

Run::

    ./scripts/run_examples.sh -e ans:13 -n 1 -t 180
"""

from __future__ import annotations

import math
import re
import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[3]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

import tests.validation.test_geo38_c39_thermal_reader as c46  # noqa: E402
import tests.validation.test_geo38_f_clinical_3t_box_heating as c31  # noqa: E402
import tests.validation.test_geo38_f_clinical_box_heating as p61  # noqa: E402
from tests.validation.test_geo38_c39_thermal_reader import (  # noqa: E402
    ENTRY59_GLOB,
    ENTRY59_NODE_PF,
    ENTRY60_GLOB,
    ENTRY60_NODE_PF,
    SENSES,
    _half_ulp_rel,
    _read_log,
    check_block,
    find_log,
    parse_heat_log,
    ratios,
    side,
)
from tests.validation.test_geo38_f_clinical_box_heating import (  # noqa: E402
    HEAT_ENV,
    HEAT_SENSE_ENV,
    RUNNING_WINDOW_S,
    SENSE_TUPLE,
    T_STAR_BOX,
    _probes,
)
from tests.validation.test_mat15_coil_driven_heating_phantom import (  # noqa: E402
    H_GEL,
    K_GEL,
    N_MAIN,
    T_INF,
)
from tests.validation.test_mat21_box_closed_form import A_BOX, B_BOX, C_BOX, RHO_C  # noqa: E402
from tests.validation.test_mat21_box_fem import T_STAR  # noqa: E402

CASE_DIR = Path(__file__).resolve().parent
SPEC_PATH = CASE_DIR / "SPEC.md"
COMPARISON_MD = CASE_DIR / "COMPARISON.md"
LOG_DIR = _REPO_ROOT / "docs" / "testing" / "logs"
C39_LOG = LOG_DIR / "20261111T103024Z_GEO-38-C39-W1.log"
C33_LOG = LOG_DIR / "20261109T103557Z_GEO-38-C33-W1-r2.log"
TAG = "[ANS-13 page]"
#: (key, glob, node, label, knob env) — entry 59 is P61's own window, entry 60 C31's.
ENTRIES = (("59", ENTRY59_GLOB, ENTRY59_NODE_PF, "1.5 T", HEAT_ENV),
           ("60", ENTRY60_GLOB, ENTRY60_NODE_PF, "3 T", c31.C31_ENV))

# ---------------------------------------------------------------------------
# Local regexes for the lines C46 does not parse (applied inside each block).
# ---------------------------------------------------------------------------
_GT = r"\[GEO-36 I32\] \[gel twin (\w+)\] "
_RX_VSUB = re.compile(_GT + r"sense (\(1,[+-]1j\)): .*?V_sub (\S+) m\^3; Robin area (\S+) m\^2")
_RX_B = re.compile(_GT + r"\(B\) int_sub q (\S+) vs G\*P (\S+):")
_RX_CN = re.compile(_GT + r"CN N=(\d+) to t\* (\S+) s, dt (\S+) s .*?\(E'\) stored (\S+) \+ "
                    r"surface (\S+) W: .*?\(S2\) (\S+) \(band")
_RX_NORM = re.compile(r"\[GEO-38 C6a\] normalisers: .*?P_box/P_acc (\S+); SAR10g\(c\*\) per box "
                      r"watt (\S+), per P_acc (\S+), per P_src")
_RX_ELAPSED = re.compile(r"^- Elapsed \(s\): (\d+)$")
_RX_NODE_KNOB = re.compile(r"FEM_EM_GEO36_S8_DRIVE_C_PF=([0-9.]+)")
_RX_C39_HEAD = re.compile(r"\[GEO-38 C46\] C39 table: ")
_RX_C39_ROW = re.compile(r"\[GEO-38 C46\]   (ccw|cw) (.+?): (\S+) \| (\S+) \| (\S+)$")
_RX_C33_LEVER = re.compile(r"\[GEO-38 C44\]   (SAR10g\(c\*\) per P_acc \[W/kg/W\]|P_box / P_acc): "
                           r"(\S+) \| (\S+) \| (\S+)$")


def _norm(tok: str) -> str:
    """A token as compared: Unicode minus to '-', a leading '+' stripped."""
    t = tok.strip().replace("−", "-")
    return t[1:] if t.startswith("+") else t


def _signed_spec(tok: str) -> str:
    """A log token as the SPEC prints it (leading '+' dropped)."""
    return tok[1:] if tok.startswith("+") else tok


# ---------------------------------------------------------------------------
# Readers (pure).
# ---------------------------------------------------------------------------
def entry_log(key):
    glob = dict((k, g) for k, g, *_ in ENTRIES)[key]
    hit = find_log(glob)
    if hit is None or hit[1] != "tracked":
        raise ValueError(f"entry {key}'s log: {hit!r} (must be tracked)")
    return hit[0]


def _one(hits, what):
    if len(hits) != 1:
        raise ValueError(f"{what}: {len(hits)} matches")
    return hits[0]


def parse_extras(path, block) -> dict:
    """The lines inside ``block`` that C46 does not parse, per sense, plus the footer,
    the ``memory.peak`` bytes and the command line (outside the block)."""
    lines = _read_log(path).splitlines()
    start, end = block["head"]["line"], block["end"]
    per = {s: {} for s in SENSES}
    norm = []
    for i in range(start, end + 1):
        ln = lines[i - 1]
        for key, rx in (("vsub", _RX_VSUB), ("b", _RX_B), ("cn", _RX_CN)):
            m = rx.search(ln)
            if m:
                sn = m.group(1)
                if sn not in per or key in per[sn]:
                    raise ValueError(f"{key} line for {sn!r} at :{i} unknown or repeated")
                if key == "vsub":
                    if SENSE_TUPLE[sn] != m.group(2):
                        raise ValueError(f"V_sub line :{i}: sense {sn} tuple {m.group(2)}")
                    per[sn][key] = {"v_sub": m.group(3), "robin": m.group(4), "line": i}
                elif key == "b":
                    per[sn][key] = {"int_q": m.group(2), "gp": m.group(3), "line": i}
                else:
                    per[sn][key] = {"n": m.group(2), "t_star": m.group(3), "dt": m.group(4),
                                    "stored": m.group(5), "surface": m.group(6),
                                    "s2": m.group(7), "line": i}
        m = _RX_NORM.search(ln)
        if m:
            norm.append({"p_box_over_p_acc": m.group(1), "sar_per_box_w": m.group(2),
                         "sar_per_p_acc": m.group(3), "line": i})
    for sn in SENSES:
        miss = [k for k in ("vsub", "b", "cn") if k not in per[sn]]
        if miss:
            raise ValueError(f"block :{start}, sense {sn}: lacks {miss}")
    el = _one([(i, m.group(1)) for i, ln in enumerate(lines, start=1)
               if (m := _RX_ELAPSED.search(ln))], "Elapsed footer")
    pk = _one([i for i, ln in enumerate(lines, start=1) if ln.strip() == "[XL] memory.peak bytes:"],
              "memory.peak header")
    cmd = next((i for i, ln in enumerate(lines, start=1) if ln.startswith("- Command:")), None)
    if cmd is None:
        raise ValueError("no '- Command:' line")
    return {"senses": per, "norm": _one(norm, "normalisers line"),
            "elapsed": {"s": el[1], "line": el[0]},
            "peak": {"bytes": lines[pk].strip(), "line": pk + 1},
            "cmd": {"text": lines[cmd - 1], "line": cmd}}


def parse_c39(path=C39_LOG) -> dict:
    """C39's table(s) under each ``C39 table:`` header: ``[{(sense, key): (v15, v3, r, line)}]``."""
    lines = _read_log(path).splitlines()
    tables = []
    cur = None
    for i, ln in enumerate(lines, start=1):
        if _RX_C39_HEAD.search(ln):
            cur = {"head": i, "rows": {}}
            tables.append(cur)
            continue
        m = _RX_C39_ROW.search(ln)
        if m and cur is not None:
            k = (m.group(1), m.group(2))
            if k in cur["rows"]:
                raise ValueError(f"C39 row {k} repeated at :{i}")
            cur["rows"][k] = (m.group(3), m.group(4), m.group(5), i)
        elif cur is not None and cur["rows"]:
            cur = None
    return {"log": Path(str(path)).name, "tables": tables}


def parse_c33_levers(path=C33_LOG) -> dict:
    """C33's lever lines: the first occurrence of each (``(v15, v3, ratio, line, count)``)."""
    out, count = {}, {}
    for i, ln in enumerate(_read_log(path).splitlines(), start=1):
        m = _RX_C33_LEVER.search(ln)
        if m:
            count[m.group(1)] = count.get(m.group(1), 0) + 1
            out.setdefault(m.group(1), (m.group(2), m.group(3), m.group(4), i))
    return {k: v + (count[k],) for k, v in out.items()}


_E = r"(\d\.\d+e[+-]\d+)"
_CITE1 = r" \(`:(\d+)`\)"


def parse_spec13(path=SPEC_PATH) -> dict:
    """Pure: the SPEC's record table (``:140–148``), its per-deposited-watt and ratio
    sentences, row 2's coordinates, row 3's calorimetry figure and the thermal table,
    by regex (never by line number), as strings with their SPEC line numbers."""
    text = Path(path).read_text()
    lines = text.splitlines()
    rows = {
        "header": (r"^\| Quantity, ccw \(1, −1j\) \| 1\.5 T — entry 59, `(\S+)` \((\d+) s, "
                   r"(\d+\.\d+) GiB, node ([0-9.]+) pF\) \| 3 T — entry 60, `(\S+)` \((\d+) s, "
                   r"(\d+\.\d+) GiB, node ([0-9.]+) pF\) \|$",
                   ("log59", "el59", "gib59", "node59", "log60", "el60", "gib60", "node60")),
        "factor": (r"^\| P_box / P_acc \(the factor; log \"per-accepted-watt factor\"\) \| " + _E
                   + _CITE1 + r" \| " + _E + _CITE1 + r" \|$",
                   ("f59", "c_f59", "f60", "c_f60")),
        "cstar": (r"^\| θ\(t\\\*\) at c\\\* per P_acc watt \[K/W\] \| " + _E + r" at (\([^)]*\))"
                  + _CITE1 + r" \| " + _E + r" at (\([^)]*\))" + _CITE1 + r" \|$",
                  ("t59", "pt59", "c_t59", "t60", "pt60", "c_t60")),
        "probes": (r"^\| θ\(t\\\*\) at P0 … P5 per P_acc watt \[K/W\] \| ((?:\S+, ){5}\S+) "
                   r"\(`:(\d+)–(\d+)`\) \| ((?:\S+, ){5}\S+) \(`:(\d+)–(\d+)`\) \|$",
                   ("p59", "c_p59a", "c_p59b", "p60", "c_p60a", "c_p60b")),
        "cw": (r"^\| cw \(1, \+1j\): θ\(t\\\*\) at c\\\* per P_acc watt; factor \| " + _E + "; " + _E
               + r" \(`:(\d+)`, `:(\d+)`\) \| " + _E + "; " + _E + r" \(`:(\d+)`, `:(\d+)`\) \|$",
               ("t59", "f59", "c_t59", "c_f59", "t60", "f60", "c_t60", "c_f60")),
        "ball": (r"^\| 1 g ball mean θ\(t\\\*\) at c\\\* at 1 W deposited \[K\] \| " + _E + _CITE1
                 + r" \| " + _E + _CITE1 + r" \|$", ("b59", "c_b59", "b60", "c_b60")),
        "surface": (r"^\| Surface heat at t\\\* per deposited watt \(row 4\) \| " + _E
                    + r" W \(`:(\d+)`, \(E′\)\) \| " + _E + r" W \(`:(\d+)`, \(E′\)\) \|$",
                    ("s59", "c_s59", "s60", "c_s60")),
        "ra": (r"^\| 360 s running-mean maximum at c\\\* per P_acc watt \| " + _E
               + r" K \(still rising at t\\\*\) \| " + _E + r" K \(still rising at t\\\*\) \|$",
               ("r59", "r60")),
    }
    out: dict = {}
    for lineno, line in enumerate(lines, start=1):
        for key, (rx, names) in rows.items():
            m = re.search(rx, line)
            if m:
                if key in out:
                    raise ValueError(f"SPEC {key}: second match at :{lineno}")
                out[key] = {"line": lineno, **dict(zip(names, m.groups()))}
    missing = [k for k in rows if k not in out]
    if missing:
        raise ValueError(f"SPEC record rows not found: {missing}")
    h = out["header"]["line"]
    if lines[h] != "|---|---|---|":
        raise ValueError(f"SPEC record table separator at :{h + 1}: {lines[h]!r}")
    want = [h + 2 + j for j in range(7)]
    got = [out[k]["line"] for k in ("factor", "cstar", "probes", "cw", "ball", "surface", "ra")]
    if got != want:
        raise ValueError(f"SPEC record rows at {got}, expected {want}")
    out["table"] = {"line": h, "text": "\n".join(lines[h - 1:h + 8])}
    flat = re.sub(r"\s+", " ", text)
    m = re.findall(r"review arithmetic: c\\\* " + _E + r" K at 1\.5 T, " + _E + r" K at 3 T", flat)
    out["per_dep"] = dict(zip(("c59", "c60"), _one(m, "SPEC per-deposited-watt sentence")))
    m = re.findall(r"is (\d\.\d+) ccw / (\d\.\d+) cw\*\* \(C39, `(\S+?):(\d+)`, `:(\d+)`\), against "
                   r"the SAR10g lever (\d\.\d+) \(C33\) and the closed-form pictures' (\d\.\d+) \(P83",
                   flat)
    out["ratio"] = dict(zip(("ccw", "cw", "log", "c_ccw", "c_cw", "sar_lever", "fac_lever"),
                            _one(m, "SPEC ratio sentence")))
    row2 = _one([(i, ln) for i, ln in enumerate(lines, start=1) if ln.startswith("| 2 | **Point rise")],
                "SPEC row 2")
    m = re.search(r"\*\*c\\\*\*\* \(1\.5 T: (\([^)]*\)); 3 T: (\([^)]*\))", row2[1])
    if not m:
        raise ValueError("SPEC row 2: c* coordinates not found")
    pts = re.findall(r"\*\*(P\d)\*\*[^(]*?(\([^)]*\))", row2[1])
    out["row2"] = {"line": row2[0], "cstar59": m.group(1), "cstar60": m.group(2), "probes": pts}
    m = re.findall(r"equals 1 W / \(ρc·V\) = " + _E + " K/s", flat)
    out["row3"] = {"slope": _one(m, "SPEC row 3 calorimetry figure")}
    th = {}
    for key, rx in (
            ("box", r"(\d\.\d+) \(x\) × (\d\.\d+) \(y\) × (\d\.\d+) \(z\) m, centred on the origin"),
            ("vr", r"V = (\S+) m³, Robin area (\S+) m²; entry 59 `:(\d+)`"),
            ("k", r"^\| Thermal conductivity k \| \*\*(\S+) W/\(m·K\)\*\*"),
            ("rhoc", r"^\| Volumetric heat capacity ρc \| \*\*(\S+) J/\(m³·K\)\*\*"),
            ("h", r"\*\*h = (\S+) W/\(m²·K\)\*\* to ambient \*\*T∞ = (\S+) °C\*\*"),
            ("tstar", r"^\| Duration t\\\* \| \*\*(\S+) s\*\* \(= 0\.2"),
            ("dt", r"^\| Time step \| ≤ (\S+) s \((\d+) steps or finer; ours is Crank–Nicolson, "
                   r"N = (\d+)\)")):
        hits = [(i, mm.groups()) for i, ln in enumerate(lines, start=1) if (mm := re.search(rx, ln))]
        i, g = _one(hits, f"SPEC thermal table {key}")
        th[key] = {"line": i, "g": g}
    out["thermal"] = th
    return out


def read_all(log59=None, log60=None, spec_path=SPEC_PATH, c39_path=C39_LOG, swap=False) -> dict:
    """Parse both entries (each through ``parse_heat_log`` once), C39's log and the SPEC.
    ``swap`` reads entry 59's block as the 3 T column (the (N-swap) control)."""
    paths = {"59": entry_log("59") if log59 is None else Path(log59),
             "60": entry_log("60") if log60 is None else Path(log60)}
    if swap:
        paths = {"59": paths["60"], "60": paths["59"]}
    blocks = {k: parse_heat_log(p) for k, p in paths.items()}
    extras = {k: parse_extras(paths[k], blocks[k]) for k in paths}
    return {"paths": paths, "blocks": blocks, "extras": extras, "c39": parse_c39(c39_path),
            "c33": parse_c33_levers(), "spec": parse_spec13(spec_path)}


# ---------------------------------------------------------------------------
# Rendering.
# ---------------------------------------------------------------------------
def _coord(v: float) -> str:
    return "0" if v == 0.0 else format(v, "+.3f").replace("-", "−")


def _point_str(pt_tok: str) -> str:
    vals = [float(v) for v in pt_tok.strip("()").split(",")]
    return "(" + ", ".join(_coord(v) for v in vals) + ")"


def _gib(b: str) -> str:
    return f"{int(b) / 2 ** 30:.2f}"


def render_table(ctx) -> str:
    """``SPEC.md:140–148`` rendered from the parsed tokens."""
    b, x = ctx["blocks"], ctx["extras"]

    def col(k):
        return (f"entry {k}, `{b[k]['log']}` ({x[k]['elapsed']['s']} s, "
                f"{_gib(x[k]['peak']['bytes'])} GiB, node {b[k]['head']['node']} pF)")

    def cell(k, row):
        s = b[k]["senses"]
        ccw, cw = s["ccw"], s["cw"]
        if row == "factor":
            src = ccw["p61"]["src"]
            return f"{src['factor']} (`:{src['line']}`)"
        if row == "cstar":
            h = ccw["p61"]["hot"]
            return f"{h['theta_point']} at {_point_str(h['point'])} (`:{h['line']}`)"
        if row == "probes":
            pr = list(ccw["p61"]["probes"].values())
            return (", ".join(_signed_spec(p["theta"]) for p in pr)
                    + f" (`:{pr[0]['line']}–{pr[-1]['line']}`)")
        if row == "cw":
            h, src = cw["p61"]["hot"], cw["p61"]["src"]
            return f"{h['theta_point']}; {src['factor']} (`:{h['line']}`, `:{src['line']}`)"
        if row == "ball":
            bl = ccw["i32"]["ball"]
            return f"{bl['theta']} (`:{bl['line']}`)"
        if row == "surface":
            cn = x[k]["senses"]["ccw"]["cn"]
            return f"{float(cn['surface']):.9e} W (`:{cn['line']}`, (E′))"
        if row == "ra":
            return f"{ccw['p61']['hot']['ra_max']} K (still rising at t\\*)"
        raise KeyError(row)

    names = [("factor", 'P_box / P_acc (the factor; log "per-accepted-watt factor")'),
             ("cstar", "θ(t\\*) at c\\* per P_acc watt [K/W]"),
             ("probes", "θ(t\\*) at P0 … P5 per P_acc watt [K/W]"),
             ("cw", "cw (1, +1j): θ(t\\*) at c\\* per P_acc watt; factor"),
             ("ball", "1 g ball mean θ(t\\*) at c\\* at 1 W deposited [K]"),
             ("surface", "Surface heat at t\\* per deposited watt (row 4)"),
             ("ra", f"{RUNNING_WINDOW_S:.0f} s running-mean maximum at c\\* per P_acc watt")]
    L = [f"| Quantity, ccw (1, −1j) | 1.5 T — {col('59')} | 3 T — {col('60')} |", "|---|---|---|"]
    for row, label in names:
        L.append(f"| {label} | {cell('59', row)} | {cell('60', row)} |")
    return "\n".join(L)


# ---------------------------------------------------------------------------
# Checks (pure): each returns ``[(name, ok, detail)]``.
# ---------------------------------------------------------------------------
def x_pairs(ctx) -> list:
    """``[(spec_field, spec_token, parsed_token, kind)]`` — one per SPEC record-table
    token (kind ``tok`` a value, ``cite`` a line number)."""
    sp, b, x = ctx["spec"], ctx["blocks"], ctx["extras"]
    pairs = []
    for k in ("59", "60"):
        s = b[k]["senses"]
        ccw, cw = s["ccw"], s["cw"]
        h = sp["header"]
        pairs += [(f"header.log{k}", h[f"log{k}"], b[k]["log"], "tok"),
                  (f"header.el{k}", h[f"el{k}"], x[k]["elapsed"]["s"], "tok"),
                  (f"header.gib{k}", h[f"gib{k}"], _gib(x[k]["peak"]["bytes"]), "tok"),
                  (f"header.node{k}", h[f"node{k}"], b[k]["head"]["node"], "tok"),
                  (f"factor.f{k}", sp["factor"][f"f{k}"], ccw["p61"]["src"]["factor"], "tok"),
                  (f"factor.c_f{k}", sp["factor"][f"c_f{k}"], ccw["p61"]["src"]["line"], "cite"),
                  (f"cstar.t{k}", sp["cstar"][f"t{k}"], ccw["p61"]["hot"]["theta_point"], "tok"),
                  (f"cstar.pt{k}", sp["cstar"][f"pt{k}"], _point_str(ccw["p61"]["hot"]["point"]),
                   "tok"),
                  (f"cstar.c_t{k}", sp["cstar"][f"c_t{k}"], ccw["p61"]["hot"]["line"], "cite")]
        pr = list(ccw["p61"]["probes"].values())
        spec_pr = [t.strip() for t in sp["probes"][f"p{k}"].split(",")]
        for j, (st, pv) in enumerate(zip(spec_pr, pr)):
            pairs.append((f"probes.p{k}[{j}]", st, pv["theta"], "tok"))
        pairs += [(f"probes.c_p{k}a", sp["probes"][f"c_p{k}a"], pr[0]["line"], "cite"),
                  (f"probes.c_p{k}b", sp["probes"][f"c_p{k}b"], pr[-1]["line"], "cite"),
                  (f"cw.t{k}", sp["cw"][f"t{k}"], cw["p61"]["hot"]["theta_point"], "tok"),
                  (f"cw.f{k}", sp["cw"][f"f{k}"], cw["p61"]["src"]["factor"], "tok"),
                  (f"cw.c_t{k}", sp["cw"][f"c_t{k}"], cw["p61"]["hot"]["line"], "cite"),
                  (f"cw.c_f{k}", sp["cw"][f"c_f{k}"], cw["p61"]["src"]["line"], "cite"),
                  (f"ball.b{k}", sp["ball"][f"b{k}"], ccw["i32"]["ball"]["theta"], "tok"),
                  (f"ball.c_b{k}", sp["ball"][f"c_b{k}"], ccw["i32"]["ball"]["line"], "cite"),
                  (f"surface.s{k}", sp["surface"][f"s{k}"],
                   f"{float(x[k]['senses']['ccw']['cn']['surface']):.9e}", "tok"),
                  (f"surface.c_s{k}", sp["surface"][f"c_s{k}"], x[k]["senses"]["ccw"]["cn"]["line"],
                   "cite"),
                  (f"ra.r{k}", sp["ra"][f"r{k}"], ccw["p61"]["hot"]["ra_max"], "tok")]
        if len(spec_pr) != len(pr):
            pairs.append((f"probes.p{k}.count", str(len(spec_pr)), str(len(pr)), "tok"))
    return pairs


def _pair_ok(spec_s, parsed, kind):
    if kind == "cite":
        return int(spec_s) == int(parsed)
    return _norm(str(spec_s)) == _norm(str(parsed))


def _spec_record_fields(sp) -> set:
    """Every token the SPEC's record table carries (the probe lists split into six)."""
    scalar = {f"{k}.{f}" for k in ("header", "factor", "cstar", "probes", "cw", "ball", "surface", "ra")
              for f in sp[k] if f not in ("line", "p59", "p60")}
    probes = {f"probes.p{k}[{j}]" for k in ("59", "60")
              for j in range(len(sp["probes"][f"p{k}"].split(",")))}
    return scalar | probes


def _close(a, b):
    return math.isclose(float(a), float(b), rel_tol=0.0, abs_tol=1e-12)


def check_x(ctx) -> list:
    b, x, sp = ctx["blocks"], ctx["extras"], ctx["spec"]
    names = list(_probes())
    res = [("(X) c31.p61 is p61 (the 3 T column is P61's body under C21's patches)",
            c31.p61 is p61, f"{c31.p61.__name__}")]
    for key, _glob, node, label, knob in ENTRIES:
        bb = b[key]
        res += check_block(bb, node, names, f"{label} (entry {key})")
        cmd = x[key]["cmd"]["text"]
        nk = _RX_NODE_KNOB.findall(cmd)
        ok = (f"{knob}=1" in cmd and f"{HEAT_SENSE_ENV}=both" in cmd
              and nk == [bb["head"]["node"]] and bb["head"]["node"] == node)
        res.append((f"(X) custody {label}: command carries {knob}=1, {HEAT_SENSE_ENV}=both, node",
                    ok, f":{x[key]['cmd']['line']} knob {f'{knob}=1' in cmd}, sense "
                        f"{f'{HEAT_SENSE_ENV}=both' in cmd}, node {nk} vs block {bb['head']['node']}"))
    # SPEC row 2's seven coordinates.
    pr = _probes()
    r2 = sp["row2"]
    spec_pts = {lab: tuple(float(_norm(v)) for v in pt.strip("()").split(",")) for lab, pt in r2["probes"]}
    bad = []
    for nm in names:
        lab = nm.split()[0]
        want = tuple(float(v) for v in pr[nm])
        got = spec_pts.get(lab)
        if got is None or len(got) != 3 or not all(_close(a, c) for a, c in zip(got, want)):
            bad.append((lab, got, want))
    for k, ck in (("59", "cstar59"), ("60", "cstar60")):
        hp = tuple(float(v) for v in b[k]["senses"]["ccw"]["p61"]["hot"]["point"].strip("()").split(","))
        got = tuple(float(_norm(v)) for v in r2[ck].strip("()").split(","))
        if not all(_close(a, c) for a, c in zip(got, hp)):
            bad.append((ck, got, hp))
    res.append(("(X) SPEC row 2's seven coordinates == _probes() and the parsed hotspot points",
                not bad and len(spec_pts) == 6,
                f"SPEC :{r2['line']}: {len(spec_pts)} probes + 2 c*; mismatches {bad}"))
    # Every record-table token.
    pairs = x_pairs(ctx)
    for name, spec_s, parsed, kind in pairs:
        res.append((f"(X) {name}", _pair_ok(spec_s, parsed, kind),
                    f"SPEC {spec_s} vs parsed {parsed} ({kind})"))
    fields = _spec_record_fields(sp)
    pn = [n for n, *_ in pairs]
    ok = sorted(pn) == sorted(fields) and len(set(pn)) == len(pn)
    res.append(("(X) every SPEC record-table token matched to exactly one parsed field",
                ok, f"{len(set(pn) & fields)} of {len(fields)}, unmatched "
                    f"{sorted(fields - set(pn))}, extra {sorted(set(pn) - fields)}"))
    # The thermal table.
    th = sp["thermal"]
    v59 = x["59"]["senses"]["ccw"]["vsub"]
    cn59 = x["59"]["senses"]["ccw"]["cn"]
    box = th["box"]["g"]
    checks = [
        ("k == K_GEL", float(th["k"]["g"][0]) == K_GEL, f"{th['k']['g'][0]} vs {K_GEL}"),
        ("ρc == RHO_C", float(th["rhoc"]["g"][0]) == RHO_C, f"{th['rhoc']['g'][0]} vs {RHO_C}"),
        ("h == H_GEL", float(th["h"]["g"][0]) == H_GEL, f"{th['h']['g'][0]} vs {H_GEL}"),
        ("T∞ == T_INF", float(th["h"]["g"][1]) == T_INF, f"{th['h']['g'][1]} vs {T_INF}"),
        ("t* == f'{T_STAR:.1f}' and T_STAR_BOX", th["tstar"]["g"][0] == f"{T_STAR:.1f}"
         and T_STAR == T_STAR_BOX, f"{th['tstar']['g'][0]} vs {T_STAR:.1f}"),
        ("Δt == T_STAR / N_MAIN and the log's dt", math.isclose(float(th["dt"]["g"][0]),
                                                               T_STAR / N_MAIN, rel_tol=1e-12)
         and cn59["dt"] == f"{T_STAR / N_MAIN:.6f}",
         f"{th['dt']['g'][0]} vs {T_STAR / N_MAIN:.6f}; log dt {cn59['dt']} (`:{cn59['line']}`)"),
        ("N == N_MAIN (both SPEC figures and both CN lines)",
         int(th["dt"]["g"][1]) == N_MAIN == int(th["dt"]["g"][2])
         and all(x[k]["senses"][s]["cn"]["n"] == str(N_MAIN) for k in ("59", "60") for s in SENSES),
         f"{th['dt']['g'][1]}, {th['dt']['g'][2]} vs {N_MAIN}"),
        ("box == 2 × (B_BOX, C_BOX, A_BOX)",
         list(box) == [f"{2 * B_BOX:.2f}", f"{2 * C_BOX:.2f}", f"{2 * A_BOX:.2f}"],
         f"{box} vs {(2 * B_BOX, 2 * C_BOX, 2 * A_BOX)}"),
        ("V == the V_sub token", math.isclose(float(th["vr"]["g"][0]), float(v59["v_sub"]),
                                              rel_tol=1e-12), f"{th['vr']['g'][0]} vs {v59['v_sub']}"),
        ("Robin area == the token", math.isclose(float(th["vr"]["g"][1]), float(v59["robin"]),
                                                 rel_tol=1e-12),
         f"{th['vr']['g'][1]} vs {v59['robin']}"),
        ("V_sub / Robin area tokens identical across logs and senses",
         len({(x[k]["senses"][s]["vsub"]["v_sub"], x[k]["senses"][s]["vsub"]["robin"])
              for k in ("59", "60") for s in SENSES}) == 1, "4 lines"),
    ]
    for n, ok, d in checks:
        res.append((f"(X) thermal table {n}", ok, d))
    return res


def check_d(ctx) -> list:
    b, x, sp, c39 = ctx["blocks"], ctx["extras"], ctx["spec"], ctx["c39"]
    res = []

    def rel(tag, val, tok, inputs):
        bound = _half_ulp_rel(tok) + sum(_half_ulp_rel(t) for t in inputs)
        err = abs(val - float(tok)) / abs(float(tok))
        res.append((tag, err <= bound, f"{val:.6e} vs {tok}: rel {err:.3e} (<= {bound:.3e})"))

    for key, _g, _n, label, _k in ENTRIES:
        for sn in SENSES:
            s = b[key]["senses"][sn]
            fac = s["p61"]["src"]["factor"]
            for pn, pv in s["p61"]["probes"].items():
                itok = s["i32"]["probes"][pn]["theta"]
                rel(f"(D-dep) {label} {sn} {pn}: θ per P_acc / factor == I32 per-deposited-watt",
                    float(pv["theta"]) / float(fac), itok, (pv["theta"], fac))
    for key, ck in (("59", "c59"), ("60", "c60")):
        s = b[key]["senses"]["ccw"]["p61"]
        rel(f"(D-c*) entry {key} ccw θ(c*) per P_acc / factor == SPEC per-deposited-watt",
            float(s["hot"]["theta_point"]) / float(s["src"]["factor"]), sp["per_dep"][ck],
            (s["hot"]["theta_point"], s["src"]["factor"]))
    # (D-ratio)
    r = ratios(side(b["59"]), side(b["60"]))
    tabs = c39["tables"]
    first = tabs[0]["rows"] if tabs else {}
    bad = []
    n = 0
    for sn in SENSES:
        for k, v in r[sn].items():
            n += 1
            got = first.get((sn, k))
            if got is None or got[2] != f"{v:.6e}":
                bad.append((sn, k, f"{v:.6e}", got and got[2]))
    res.append(("(D-ratio) C46 ratios(side(59), side(60)) .6e == C39's strings",
                n == 22 and len(first) == 22 and not bad,
                f"{n} ratios vs {len(first)} C39 rows (`{c39['log']}:{tabs[0]['head'] if tabs else '-'}`); "
                f"mismatches {bad}"))
    same = len(tabs) >= 1 and all({k: v[:3] for k, v in t["rows"].items()}
                                  == {k: v[:3] for k, v in first.items()} for t in tabs)
    res.append(("(D-ratio) every printed copy of C39's table identical", same,
                f"{len(tabs)} copies at {[t['head'] for t in tabs]}"))
    hk = "hotspot point theta per P_acc"
    rt = sp["ratio"]
    for sn in SENSES:
        row = first.get((sn, hk))
        ok = (row is not None and rt[sn] == f"{r[sn][hk]:.4f}"
              and int(rt[f"c_{sn}"]) == row[3] and rt["log"] == c39["log"])
        res.append((f"(D-ratio) SPEC {sn} {rt[sn]} == f'{{r:.4f}}' of the C39 line",
                    ok, f"{r[sn][hk]:.6e} → {r[sn][hk]:.4f}; SPEC cite `{rt['log']}:{rt[f'c_{sn}']}` "
                        f"vs parsed :{row[3] if row else '-'}"))
    f59 = float(b["59"]["senses"]["ccw"]["p61"]["src"]["factor"])
    f60 = float(b["60"]["senses"]["ccw"]["p61"]["src"]["factor"])
    res.append(("(D-lever) SPEC fac lever == f'{fac60 / fac59:.4f}'", rt["fac_lever"] == f"{f60 / f59:.4f}",
                f"{f60 / f59:.6e} vs SPEC {rt['fac_lever']}"))
    s59 = float(x["59"]["norm"]["sar_per_p_acc"])
    s60 = float(x["60"]["norm"]["sar_per_p_acc"])
    res.append(("(D-lever) SPEC SAR10g lever == f'{s60 / s59:.4f}'", rt["sar_lever"] == f"{s60 / s59:.4f}",
                f"{s60 / s59:.6e} vs SPEC {rt['sar_lever']} (normalisers `:{x['59']['norm']['line']}` / "
                f"`:{x['60']['norm']['line']}`)"))
    for k in ("59", "60"):
        v = x[k]["senses"]["ccw"]["vsub"]["v_sub"]
        res.append((f"(D-cal) entry {k}: SPEC row 3 == f'{{1 / (RHO_C * V):.3e}}'",
                    sp["row3"]["slope"] == f"{1.0 / (RHO_C * float(v)):.3e}",
                    f"{1.0 / (RHO_C * float(v)):.6e} vs SPEC {sp['row3']['slope']} (V_sub {v})"))
    return res


def check_r(ctx) -> list:
    render = render_table(ctx)
    t = ctx["spec"]["table"]
    return [("(R) record table rendered from the logs == SPEC bytes", render == t["text"],
             f"{len(render.encode())} B vs SPEC :{t['line']}–{t['line'] + 8} {len(t['text'].encode())} B")]


def check_r_cmp():
    """(R-cmp): live only once item 959 (C54) has written ANS-13's COMPARISON.md."""
    if not COMPARISON_MD.exists():
        return [], "absent — item 959 (C54)"
    from tests.validation import test_geo38_c33_field_strength_reader as c33

    page = c33._load_page()  # collective (EX-193): every rank reaches this line
    fn = getattr(page, "render_ans13", None)
    if fn is None:
        return [("(R-cmp) the ans:11 page exposes render_ans13()", False,
                 "COMPARISON.md exists but the ans:11 page has no render_ans13()")], "inconsistent"
    render = fn()
    text = COMPARISON_MD.read_text()
    return [("(R-cmp) re-render through the 'ans13' case == COMPARISON.md bytes", render == text,
             f"{len(render.encode())} B vs {len(text.encode())} B")], "live"


def spec_line_table(ctx) -> list:
    """``[(item, SPEC-cited line, parsed line)]`` for the cites outside the record table
    (printed, never asserted: reading 2)."""
    th = ctx["spec"]["thermal"]
    v = ctx["extras"]["59"]["senses"]["ccw"]["vsub"]
    return [("thermal table V / Robin area (entry 59)", int(th["vr"]["g"][2]), v["line"])]


def all_results(ctx):
    rc, msg = check_r_cmp()
    return check_x(ctx) + check_d(ctx) + check_r(ctx) + rc, msg


def main() -> int:
    from mpi4py import MPI

    comm = MPI.COMM_WORLD
    ctx = read_all()
    results, cmp_msg = all_results(ctx)
    if comm.rank == 0:
        b, x = ctx["blocks"], ctx["extras"]
        for key, glob, node, label, _k in ENTRIES:
            print(f"{TAG} entry {key} ({label}): {b[key]['log']} ({find_log(glob)[1]}); block "
                  f":{b[key]['head']['line']} → gates :{b[key]['end']}; node {b[key]['head']['node']} pF, "
                  f"degree {b[key]['head']['degree']}, -n {b[key]['head']['ranks']}; "
                  f"{x[key]['elapsed']['s']} s, {_gib(x[key]['peak']['bytes'])} GiB")
        print(f"{TAG} degree-1 records at the degree-1 roots; the adjudication columns are "
              "entries 69 / 70 (degree 2, PENDING on P87)")
        for name, ok, detail in results:
            print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}")
        print(f"{TAG} (R-cmp) {cmp_msg}")
        print(f"{TAG} per probe: θ(t*) per P_acc watt | factor | per deposited watt (quotient) | "
              "I32 gel-twin token (per deposited watt)")
        for key, _g, _n, label, _k in ENTRIES:
            for sn in SENSES:
                s = b[key]["senses"][sn]
                fac = s["p61"]["src"]["factor"]
                h = s["p61"]["hot"]
                q = float(h["theta_point"]) / float(fac)
                print(f"{TAG}   {label} {sn} c* {h['label']}: {h['theta_point']} | {fac} | {q:.6e} | "
                      f"1 g ball {s['i32']['ball']['theta']}"
                      + ("" if sn == "ccw" else " (no printed point comparator: printed only)"))
                for pn, pv in s["p61"]["probes"].items():
                    print(f"{TAG}   {label} {sn} {pn}: {pv['theta']} | {fac} | "
                          f"{float(pv['theta']) / float(fac):.6e} | {s['i32']['probes'][pn]['theta']}")
        for key, _g, _n, label, _k in ENTRIES:
            n = x[key]["norm"]
            print(f"{TAG} {label} factor identity: P61 factor P_gel/P_acc "
                  f"{b[key]['senses']['ccw']['p61']['src']['factor']} (tag 403 counted as gel) is the "
                  f"SPEC's P_box / P_acc; the normaliser's 6-digit P_box/P_acc {n['p_box_over_p_acc']} "
                  f"(`:{n['line']}`) is printed, never used")
            for sn in SENSES:
                e = x[key]["senses"][sn]
                print(f"{TAG} {label} {sn}: (B) int_sub q {e['b']['int_q']} (`:{e['b']['line']}`); (E') "
                      f"stored {e['cn']['stored']} + surface {e['cn']['surface']} W; (S2) {e['cn']['s2']} "
                      f"(`:{e['cn']['line']}`)")
        lev = ctx["c33"]
        f59 = float(b["59"]["senses"]["ccw"]["p61"]["src"]["factor"])
        f60 = float(b["60"]["senses"]["ccw"]["p61"]["src"]["factor"])
        s59 = float(x["59"]["norm"]["sar_per_p_acc"])
        s60 = float(x["60"]["norm"]["sar_per_p_acc"])
        for k, mine in (("SAR10g(c*) per P_acc [W/kg/W]", s60 / s59), ("P_box / P_acc", f60 / f59)):
            got = lev.get(k)
            print(f"{TAG} (D-lever, printed) C33 {k}: {got[2] if got else 'absent'} "
                  f"(`{C33_LOG.name}:{got[3] if got else '-'}`, first of {got[4] if got else 0} "
                  f"copies) | page {mine:.6e}")
        tabs = ctx["c39"]["tables"]
        print(f"{TAG} C39 table multiplicity in `{ctx['c39']['log']}`: {len(tabs)} copies at "
              f"{[t['head'] for t in tabs]}, {[len(t['rows']) for t in tabs]} rows")
        print(f"{TAG} SPEC record table (re-rendered):")
        for ln in render_table(ctx).splitlines():
            print(f"{TAG}   {ln}")
        print(f"{TAG} SPEC line cites outside the record table (printed, never asserted):")
        for item, a, p in spec_line_table(ctx):
            print(f"{TAG}   {item}: SPEC :{a} | parsed :{p} | {'same' if a == p else 'MISMATCH'}")
        print(f"{TAG} writes nothing; no AED read or write")
    failed = [n for n, ok, _ in results if not ok]
    if failed:
        if comm.rank == 0:
            print(f"{TAG} FAILED: {failed}")
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
