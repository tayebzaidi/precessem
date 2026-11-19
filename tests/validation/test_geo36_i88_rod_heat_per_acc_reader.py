"""`GEO-36` / chain I I88 — the rod's tip 1 g θ per accepted watt, 1.5 T and 3 T
(§9 item 946; docs/planning/chunks/GEO-36.md "§9 item 946").  No solve, no mesh.

Entry 67's readout is "the 3 T / 1.5 T ratio of the tip θ per accepted watt (C40's
twin, a record)".  The normaliser is P61's (``test_geo38_f_clinical_box_heating.py``
docstring): θ_acc = θ · P_gel / P_acc(cw), P_acc the lumped line cut (ruling
2026-10-12, ``_p_acc_spec``).  The I32 heating log prints θ at G = 1 / P_gel(twin)
with no S-matrix and no P_acc, so the rod route reaches P_acc by **string custody**
between logs: entry 48 (I36a, same mesh / drive / node) prints S2 on both arms,
entry 59 (P61) prints the twin's cw factor.

* twin θ_acc = θ_twin × factor(cw)  (entry 59's printed factor);
* rod θ_acc = θ_rod × P_gel_twin / P_acc_rod, P_acc_rod = ``_p_acc_spec(S2_M, (1, 1j),
  V_s, 50)``; the twin-P_acc version and P_acc_rod / P_acc_twin − 1 printed beside.

The sense is (1,+1j) — C6's **cw**; the ccw P_acc of entries 52 / 58 is never used.

Gates (no band; string equality or half a unit of the printed digit): (X48) (X55)
(X59), (Cust-Z), (Cust-P), (D-Pacc) at C33's 1e-8, (D-f), (Stub-3T), (Form-63),
(Self), (Pend).  Controls: (N-sense), (N-drive), (N-digit), (N-blocks).

Knob ``FEM_EM_GEO36_I88``: unset ⇒ the gated read skips; ``1`` ⇒ entries 48 / 55 / 59
required, 60 / 61 / 63 / 67 optional (absent ⇒ "pending", nothing written); ``both``
⇒ all of 48 / 59 / 60 / 61 / 63 / 67 required, 63 replaces 55 as the 1.5 T θ source
(55 printed beside), rank 0 writes ``/workspace/logs/geo36/i88-rod-heat-acc-1p5t-3t.json``;
anything else ⇒ ``ValueError``.  No log token of >= 8 significant digits is a literal
here: every number checked is parsed at run time.
"""

from __future__ import annotations

import copy
import json
import os
import re
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

# Module objects only — never a ``test_*`` name (I68').
import tests.validation.test_geo36_c40_rod_frequency_reader as c45  # noqa: E402
import tests.validation.test_geo38_c39_thermal_reader as c46  # noqa: E402
import tests.validation.test_geo38_f_clinical_box_loaded as c6  # noqa: E402

I88_ENV = "FEM_EM_GEO36_I88"
TAG = "[GEO-36 I88]"
CW = "(1,+1j)"
CW_W = (1.0, 1.0j)
Z0 = 50.0
#: C33's (D-Pacc) bound (``test_geo38_c33_field_strength_reader.py:161–164``).
DPACC_REL = 1e-8
NODE_1P5 = "193.5996"
NODE_3T = c45.ENTRY61_NODE_PF
BALLS_1P5 = ("top", "bottom")
#: Entry 55's gates string: custody of a red window (its figure is a pre-figure).
ENTRY55_GATES = "['rod (S2-tip) 6.420e-03']"
TRACKED_LOG_DIR = REPO / "docs" / "testing" / "logs"
HELD_LOG_DIR = REPO / "logs" / "leak-held"
#: Landed windows, pinned by name.
ENTRY48_LOG = c45.ENTRY48_LOG
ENTRY55_LOG = TRACKED_LOG_DIR / "20261106T080008Z_GEO-36-step10-I32.log"
ENTRY59_LOG = TRACKED_LOG_DIR / "20261110T080008Z_GEO-38-P33.log"
#: The brief's cited custody lines (1-based).
ENTRY48_Z_LINE, ENTRY55_Z_LINE, ENTRY59_Z_LINE = 12062, 12062, 12065
ENTRY55_SHARE_LINE, ENTRY55_TWIN_PGEL_LINE = 12497, 12499
ENTRY48_ARM_M_YIN_LINE, ENTRY48_ARM_R_YIN_LINE = 12494, 12498
ENTRY48_ARM_M_S2_LINE, ENTRY48_ARM_R_S2_LINE = 12491, 12495
ENTRY59_CW_LINE, ENTRY59_CCW_LINE = 12343, 12325
#: Pending windows, found by glob.
ENTRY60_GLOB = c46.ENTRY60_GLOB
ENTRY61_GLOB = c45.ENTRY61_GLOB
ENTRY63_GLOB = "*_GEO-36-step10-I32-prime.log"
#: Provisional: entry 67's READY review confirms it against its XL_CHUNK or edits it.
ENTRY67_GLOB = "*_GEO-36-step12-I32-3T.log"
#: 3 T stubs.
C32_STUB_LOG, C32_STUB_LINE = c45.C32_STUB_LOG, c45.C32_STUB_LINE
C31_STUB_LOG, C31_STUB_LINE = c46.C31_STUB_LOG, c46.C31_STUB_LINE
C31_Z_LINE, C31_S2_LINE, C31_CW_LINE = 87, 91, 117
I85_STUB_LOG = TRACKED_LOG_DIR / "20261110T105744Z_GEO-36-I85-W0b.log"
I85_STUB_LINE = 80
I85_Z_LINE, I85_SHARE_LINE, I85_TWIN_PGEL_LINE = 81, 88, 91
C32_ARM_M_YIN_LINE, C32_ARM_R_YIN_LINE = 85, 89
#: Entry 63's form (I32b L-W0b).
I32B_FORM_LOG = TRACKED_LOG_DIR / "20261108T181211Z_GEO-36-I32b-L-W0b.log"
I32B_FORM_LINE = 95
I32B_ROD_THETA_LINES, I32B_TWIN_THETA_LINES = (121, 122), (110, 111)
HOST_JSON_DIR = Path("/workspace/logs/geo36")
I88_JSON_NAME = "i88-rod-heat-acc-1p5t-3t.json"

_on_root = c45._on_root
_read_log = c45._read_log
_half_ulp_rel = c46._half_ulp_rel


def _half_ulp_abs(tok):
    return _half_ulp_rel(tok) * abs(float(tok))


def _find(glob):
    """``c45.find_log`` over this module's (monkeypatchable) roots."""
    return c45.find_log(glob, ((TRACKED_LOG_DIR, "tracked"), (HELD_LOG_DIR, "held")))


# ---------------------------------------------------------------------------
# Pure readers.
# ---------------------------------------------------------------------------
_T = r"(\(1,[+-]1j\))"
_RX_HEAD = re.compile(r"start: -n (\d+), I32 body; node C ([0-9.]+) pF, degree (\d), sense " + _T)
_RX_GATES = re.compile(r"\[GEO-36 I32\] gates: failures (.*?) \(ASSERTED\)")
_RX_I32_Z = re.compile(
    r"\[GEO-36 I32\] \(M2\) / \(Z\) node ([0-9.]+) pF: Z_RC (\S+) Ohm, "
    r"Z_th (\S+) Ohm, V_s (\S+) V"
)
_RX_SHARE = re.compile(
    r"\[GEO-36 I32\] rod share of P_vol, sense " + _T + r": (\S+) \(P_vol (\S+) W, "
    r"gel ([^,\s]+), rod ([^;\s]+); printed\)"
)
_RX_SRC = re.compile(
    r"\[GEO-36 I32\] \[(gel twin|rod)\] sense " + _T + r": P_gel (\S+) W; G (\S+) 1/W; "
    r"G\*P (\S+) W"
)
#: One θ regex for the I32a and I32b forms (the tail after "K" differs).
_RX_THETA = re.compile(
    r"\[GEO-36 I32\] \[(gel twin|rod)\] 1 g ball \(a (\S+) mm\) at (\w+), sense "
    + _T
    + r": theta\(t\*\) (\S+) K"
)
#: Any of the three (Z) forms: I32's (M2) / (Z), I36a's, C6a's.
_RX_Z_ANY = re.compile(
    r"\[GEO-3[68] (I32|I36a|C6a)\] .*?\(Z\) node ([0-9.]+) pF: Z_RC (\S+) Ohm, "
    r"Z_th (\S+) Ohm, V_s (\S+) V"
)
_CT = r"([^,;\s]+)"
_RX_C6A_S2 = re.compile(
    r"\[GEO-38 C6a\] kinds .*?; S2 " + _CT + ", " + _CT + "; " + _CT + ", " + _CT + ";"
)
ARMS = {"gel twin": "twin", "rod": "rod"}


def parse_i32_block(path, at_line=None) -> dict:
    """Pure: one I32 block (a ``start: -n N, I32 body`` header to the next ``[GEO-36 I32]
    gates: failures`` line), tokens kept as strings.  A window log holds exactly one
    block; a stub selects its block by ``at_line`` (the header's 1-based line).
    ``ValueError`` on any missing / repeated line."""
    lines = _read_log(path).splitlines()
    name = Path(str(path)).name.removesuffix(".gz")
    heads = [i for i, ln in enumerate(lines, start=1) if _RX_HEAD.search(ln)]
    if at_line is None:
        if len(heads) != 1:
            raise ValueError(
                f"{name}: {len(heads)} I32 blocks (lines {heads[:8]}); a window "
                f"log holds exactly one"
            )
        start = heads[0]
    else:
        if at_line not in heads:
            raise ValueError(f"{name}: no I32 header at :{at_line}")
        start = at_line
    end = next((i for i in range(start, len(lines) + 1) if _RX_GATES.search(lines[i - 1])), None)
    if end is None:
        raise ValueError(f"{name}: block at :{start} has no gates line")
    inner = [n for n in heads if start < n < end]
    if inner:
        raise ValueError(f"{name}: block at :{start} holds another header at :{inner[0]}")
    h = _RX_HEAD.search(lines[start - 1])
    out = {
        "log": name,
        "head": {
            "line": start,
            "ranks": h.group(1),
            "node": h.group(2),
            "degree": h.group(3),
            "sense": h.group(4),
        },
        "gates": {"failures": _RX_GATES.search(lines[end - 1]).group(1), "line": end},
        "arms": {"twin": {"balls": {}}, "rod": {"balls": {}}},
    }

    def _put(d, key, val, i):
        if key in d:
            raise ValueError(f"{name}: {key} repeated at :{i}")
        d[key] = val

    for i in range(start, end + 1):
        ln = lines[i - 1]
        m = _RX_I32_Z.search(ln)
        if m:
            _put(
                out,
                "z",
                {
                    "node": m.group(1),
                    "z_rc": m.group(2),
                    "z_th": m.group(3),
                    "v_s": m.group(4),
                    "line": i,
                },
                i,
            )
            continue
        m = _RX_SHARE.search(ln)
        if m:
            _put(
                out,
                "share",
                {
                    "sense": m.group(1),
                    "share": m.group(2),
                    "p_vol": m.group(3),
                    "gel": m.group(4),
                    "rod": m.group(5),
                    "line": i,
                },
                i,
            )
            continue
        m = _RX_SRC.search(ln)
        if m:
            _put(
                out["arms"][ARMS[m.group(1)]],
                "src",
                {
                    "sense": m.group(2),
                    "p_gel": m.group(3),
                    "g": m.group(4),
                    "gp": m.group(5),
                    "line": i,
                },
                i,
            )
            continue
        m = _RX_THETA.search(ln)
        if m:
            _put(
                out["arms"][ARMS[m.group(1)]]["balls"],
                m.group(3),
                {"a_mm": m.group(2), "sense": m.group(4), "theta": m.group(5), "line": i},
                i,
            )
    for key in ("z", "share"):
        if key not in out:
            raise ValueError(f"{name}: block at :{start} lacks its {key} line")
    for arm, a in out["arms"].items():
        if "src" not in a or not a["balls"]:
            raise ValueError(
                f"{name}: block at :{start}, arm {arm}, lacks its source or ball lines"
            )
    if list(out["arms"]["twin"]["balls"]) != list(out["arms"]["rod"]["balls"]):
        raise ValueError(f"{name}: ball names differ between arms")
    return out


def parse_z_line(path, at_line=None, tag=None) -> dict:
    """Pure: one (Z) line (any of I32 / I36a / C6a), tokens as strings.  Without
    ``at_line`` the log must hold exactly one (of ``tag`` when given)."""
    lines = _read_log(path).splitlines()
    name = Path(str(path)).name.removesuffix(".gz")
    hits = [
        (i, m)
        for i, ln in enumerate(lines, start=1)
        if (m := _RX_Z_ANY.search(ln)) and (tag is None or m.group(1) == tag)
    ]
    if at_line is None:
        if len(hits) != 1:
            raise ValueError(f"{name}: {len(hits)} (Z) lines; a window log holds exactly one")
        i, m = hits[0]
    else:
        sel = [(i, m) for i, m in hits if i == at_line]
        if not sel:
            raise ValueError(f"{name}: no (Z) line at :{at_line}")
        i, m = sel[0]
    return {
        "log": name,
        "tag": m.group(1),
        "node": m.group(2),
        "z_rc": m.group(3),
        "z_th": m.group(4),
        "v_s": m.group(5),
        "line": i,
    }


def parse_c6a_s2(path, at_line) -> dict:
    """Pure: the C6a ``kinds …; S2 …`` line at ``at_line``, tokens as strings."""
    lines = _read_log(path).splitlines()
    m = _RX_C6A_S2.search(lines[at_line - 1])
    if m is None:
        raise ValueError(f"{Path(str(path)).name}: no C6a S2 line at :{at_line}")
    return {
        "s11": m.group(1),
        "s12": m.group(2),
        "s21": m.group(3),
        "s22": m.group(4),
        "line": at_line,
    }


def parse_p61_cw_sources(path, at_line) -> dict:
    """Pure: the ``[GEO-38 P61] sources, sense (1,+1j) (cw, heated)`` line at ``at_line``
    (C46's regex), shaped as ``parse_heat_log``'s ``senses.cw.p61.src``.  C31's :117 sits
    in its cw-only block (header :83), which ``parse_heat_log`` refuses by design."""
    lines = _read_log(path).splitlines()
    m = c46._RX_P61_SRC.search(lines[at_line - 1])
    if m is None or m.group(1) != CW or m.group(2) != "cw" or m.group(5) != "cw":
        raise ValueError(f"{Path(str(path)).name}: no cw P61 sources line at :{at_line}")
    return {
        "log": Path(str(path)).name,
        "senses": {
            "cw": {
                "p61": {
                    "src": {
                        "p_gel": m.group(3),
                        "p_acc": m.group(4),
                        "factor": m.group(6),
                        "line": at_line,
                    }
                }
            }
        },
    }


def _s2_matrix(s):
    return np.array(
        [[complex(s["s11"]), complex(s["s12"])], [complex(s["s21"]), complex(s["s22"])]]
    )


def p_acc_of(s2_tokens, v_s_tok):
    return c6._p_acc_spec(_s2_matrix(s2_tokens), CW_W, complex(v_s_tok), Z0)


def _ztok(z):
    return (z["node"], z["z_rc"], z["z_th"], z["v_s"])


# ---------------------------------------------------------------------------
# Checks: ``[(name, ok, detail)]``.
# ---------------------------------------------------------------------------
def check_x(i32, b36, heat, node, balls, name, gates=None):
    """(X..): one block each, node token, balls, every θ line (1,+1j)."""
    res = [
        (f"{name} (X-I32) header node == {node}", i32["head"]["node"] == node, i32["head"]["node"]),
        (f"{name} (X-I32) header sense {CW}", i32["head"]["sense"] == CW, i32["head"]["sense"]),
        (f"{name} (X-I32) share sense {CW}", i32["share"]["sense"] == CW, i32["share"]["sense"]),
        (f"{name} (X-I36) node == {node}", b36["node"]["c_pf"] == node, b36["node"]["c_pf"]),
        (f"{name} (X-P61) node == {node}", heat["head"]["node"] == node, heat["head"]["node"]),
    ]
    for arm, a in i32["arms"].items():
        res.append(
            (
                f"{name} (X-I32) {arm} balls == {list(balls)}",
                tuple(a["balls"]) == tuple(balls),
                list(a["balls"]),
            )
        )
        res.append(
            (f"{name} (X-I32) {arm} source sense {CW}", a["src"]["sense"] == CW, a["src"]["sense"])
        )
        for b, v in a["balls"].items():
            res.append((f"{name} (X-I32) {arm} {b} theta sense {CW}", v["sense"] == CW, v["sense"]))
    g_t, g_r = i32["arms"]["twin"]["src"]["g"], i32["arms"]["rod"]["src"]["g"]
    res.append(
        (
            f"{name} (X-G) rod G == twin G (rod theta is per twin P_gel)",
            g_t == g_r,
            f"{g_r} vs {g_t}",
        )
    )
    if gates is not None:
        res.append(
            (
                f"{name} (X-gates) failures == {gates}",
                i32["gates"]["failures"] == gates,
                i32["gates"]["failures"],
            )
        )
    return res


def check_cust_z(zs, name):
    toks = [_ztok(z) for z in zs]
    return [
        (
            f"{name} (Cust-Z) node / Z_RC / Z_th / V_s string-equal across "
            f"{[z['log'] + ':' + str(z['line']) for z in zs]}",
            all(t == toks[0] for t in toks),
            repr(toks),
        )
    ]


def check_cust_p(i32, b36, heat, name):
    t_src = i32["arms"]["twin"]["src"]
    cw = heat["senses"]["cw"]["p61"]["src"]
    sh, am, ar = i32["share"], b36["arms"]["M"], b36["arms"]["R"]
    got = (sh["share"], sh["p_vol"], sh["gel"], sh["rod"])
    want = (am["share"], am["p_vol"], am["gel"], am["rod"])
    d = abs(float(ar["p_vol"]) - float(t_src["p_gel"]))
    bound = _half_ulp_abs(ar["p_vol"])
    return [
        (
            f"{name} (Cust-P) I32 twin P_gel == P61 cw P_gel (string)",
            t_src["p_gel"] == cw["p_gel"],
            f"{i32['log']}:{t_src['line']} {t_src['p_gel']} vs {heat['log']}:{cw['line']} {cw['p_gel']}",
        ),
        (
            f"{name} (Cust-P) I32 rod share line == I36a arm M (share, P_vol, gel, rod)",
            got == want,
            f"{i32['log']}:{sh['line']} {got} vs {b36['log']}:{am['y_in']['line']} {want}",
        ),
        (
            f"{name} (Cust-P) I36a arm R P_vol vs I32 twin P_gel within half a unit",
            d <= bound,
            f"{ar['p_vol']} vs {t_src['p_gel']}: |d| {d:.3e} (<= {bound:.3e})",
        ),
    ]


def check_d(b36, heat, v_s, name):
    cw = heat["senses"]["cw"]["p61"]["src"]
    p = p_acc_of(b36["arms"]["R"]["s2"], v_s)
    err = abs(p - float(cw["p_acc"])) / abs(float(cw["p_acc"]))
    f = float(cw["p_gel"]) / float(cw["p_acc"])
    fb = _half_ulp_rel(cw["factor"]) + _half_ulp_rel(cw["p_gel"]) + _half_ulp_rel(cw["p_acc"])
    ferr = abs(f - float(cw["factor"])) / abs(float(cw["factor"]))
    return [
        (
            f"{name} (D-Pacc) _p_acc_spec(S2_R, {CW}, V_s, {Z0}) vs P61 cw P_acc",
            err <= DPACC_REL,
            f"{p:.12e} vs printed {cw['p_acc']}: rel {err:.3e} (<= {DPACC_REL:.0e})",
        ),
        (
            f"{name} (D-f) cw factor == P_gel / P_acc",
            ferr <= fb,
            f"{f:.12e} vs printed {cw['factor']}: rel {ferr:.3e} (<= {fb:.3e})",
        ),
    ]


def report(res, what):
    for n, ok, d in res:
        print(f"{TAG} {'PASS' if ok else 'FAIL'} {n}: {d}")
    bad = [(n, d) for n, ok, d in res if not ok]
    if bad:
        raise AssertionError(f"{what}: {bad}")
    return res


# ---------------------------------------------------------------------------
# Figures.
# ---------------------------------------------------------------------------
def figures(i32, b36, heat):
    """θ per accepted watt on both arms (floats), (1,+1j) = cw throughout."""
    v_s = i32["z"]["v_s"]
    factor = float(heat["senses"]["cw"]["p61"]["src"]["factor"])
    p_gel_twin = float(i32["arms"]["twin"]["src"]["p_gel"])
    p_acc_rod = p_acc_of(b36["arms"]["M"]["s2"], v_s)
    p_acc_twin = p_acc_of(b36["arms"]["R"]["s2"], v_s)
    out = {
        "factor_cw": factor,
        "p_gel_twin": p_gel_twin,
        "p_acc_rod": p_acc_rod,
        "p_acc_twin": p_acc_twin,
        "p_acc_rod_over_twin_m1": p_acc_rod / p_acc_twin - 1.0,
        "twin": {},
        "rod": {},
        "rod_twin_pacc": {},
    }
    for b, v in i32["arms"]["twin"]["balls"].items():
        out["twin"][b] = float(v["theta"]) * factor
    for b, v in i32["arms"]["rod"]["balls"].items():
        out["rod"][b] = float(v["theta"]) * p_gel_twin / p_acc_rod
        out["rod_twin_pacc"][b] = float(v["theta"]) * p_gel_twin / p_acc_twin
    return out


def ratios(f15, f3):
    """3 T / 1.5 T per tip ball and arm.  Ball names must match, else ``ValueError``."""
    out = {}
    for arm in ("twin", "rod"):
        if list(f15[arm]) != list(f3[arm]):
            raise ValueError(f"ball names differ: {list(f15[arm])} vs {list(f3[arm])}")
        for b in f15[arm]:
            out[f"{arm}.{b}"] = f3[arm][b] / f15[arm][b]
    return out


def print_figures(f, label):
    print(
        f"{TAG} {label}: {CW} = cw; factor P_gel/P_acc (cw) {f['factor_cw']:.6e}; P_acc rod "
        f"{f['p_acc_rod']:.6e} W, twin {f['p_acc_twin']:.6e} W; P_acc_rod / P_acc_twin - 1 "
        f"{f['p_acc_rod_over_twin_m1']:+.6e}"
    )
    for b in f["twin"]:
        print(
            f"{TAG} {label}: {CW} = cw; tip 1 g {b}: twin theta_acc {f['twin'][b]:.6e} K/W; rod "
            f"theta_acc {f['rod'][b]:.6e} K/W (own P_acc; twin-P_acc version "
            f"{f['rod_twin_pacc'][b]:.6e})"
        )


def print_ratios(r, label):
    for k, v in r.items():
        print(f"{TAG} {label}: {CW} = cw; ratio 3T/1.5T theta_acc {k}: {v:.6e}")


# ---------------------------------------------------------------------------
# Readers on the landed / pending logs.
# ---------------------------------------------------------------------------
def _parse_i32(path, at_line=None):
    return _on_root(parse_i32_block, path, at_line)


def _parse_i36(path, at_line=None):
    return _on_root(c45.parse_i36_log, path, at_line)


def _parse_heat(path, at_line=None):
    return _on_root(c46.parse_heat_log, path, at_line)


def _parse_z(path, at_line=None, tag=None):
    return _on_root(parse_z_line, path, at_line, tag)


def read_1p5t(i32_path=ENTRY55_LOG, gates=ENTRY55_GATES, label="1.5 T (entry 55)"):
    """X48 / X55 / X59, Cust-Z, Cust-P, D-Pacc, D-f and the figures; all asserted."""
    i32 = _parse_i32(i32_path)
    b48 = _parse_i36(ENTRY48_LOG)
    h59 = _parse_heat(ENTRY59_LOG)
    z48 = _parse_z(ENTRY48_LOG, tag="I36a")
    z59 = _parse_z(ENTRY59_LOG, tag="C6a")
    res = check_x(i32, b48, h59, NODE_1P5, BALLS_1P5, label, gates)
    res += [
        (
            f"{label} (X48) I36a balls / degree / eta",
            tuple(b48["tip"]["balls"]) == c45.ENTRY48_BALLS
            and b48["tip"]["degree"] == "1"
            and b48["tip"]["eta_valid"] == "True",
            f"{list(b48['tip']['balls'])}, {b48['tip']['degree']}, {b48['tip']['eta_valid']}",
        )
    ]
    res += c45.check_side(b48, NODE_1P5, f"{label} (X48) entry 48")
    res += check_cust_z([i32["z"] | {"log": i32["log"]}, z48, z59], label)
    res += check_cust_p(i32, b48, h59, label)
    res += check_d(b48, h59, z48["v_s"], label)
    if i32_path == ENTRY55_LOG:
        lines = {
            "(Z) 55": (i32["z"]["line"], ENTRY55_Z_LINE),
            "(Z) 48": (z48["line"], ENTRY48_Z_LINE),
            "(Z) 59": (z59["line"], ENTRY59_Z_LINE),
            "share 55": (i32["share"]["line"], ENTRY55_SHARE_LINE),
            "twin P_gel 55": (i32["arms"]["twin"]["src"]["line"], ENTRY55_TWIN_PGEL_LINE),
            "arm M Y_in 48": (b48["arms"]["M"]["y_in"]["line"], ENTRY48_ARM_M_YIN_LINE),
            "arm R Y_in 48": (b48["arms"]["R"]["y_in"]["line"], ENTRY48_ARM_R_YIN_LINE),
            "arm M S2 48": (b48["arms"]["M"]["line"], ENTRY48_ARM_M_S2_LINE),
            "arm R S2 48": (b48["arms"]["R"]["line"], ENTRY48_ARM_R_S2_LINE),
            "cw 59": (h59["senses"]["cw"]["p61"]["src"]["line"], ENTRY59_CW_LINE),
            "ccw 59": (h59["senses"]["ccw"]["p61"]["src"]["line"], ENTRY59_CCW_LINE),
        }
        res.append(
            (
                f"{label} (X-lines) the brief's cited lines",
                all(a == b for a, b in lines.values()),
                repr(lines),
            )
        )
    report(res, label)
    f = figures(i32, b48, h59)
    print_figures(
        f,
        label + (", pre-figure (red window, entry 63 decides)" if i32_path == ENTRY55_LOG else ""),
    )
    return {"i32": i32, "i36": b48, "heat": h59, "z36": z48, "zheat": z59, "fig": f}


def read_3t(p60, p61, p67, label="3 T (entries 60 / 61 / 67)"):
    i32 = _parse_i32(p67)
    b61 = _parse_i36(p61)
    h60 = _parse_heat(p60)
    z61 = _parse_z(p61, tag="I36a")
    z60 = _parse_z(p60, tag="C6a")
    i32_balls = tuple(i32["arms"]["twin"]["balls"])
    res = check_x(i32, b61, h60, NODE_3T, i32_balls, label)
    res += c45.check_side(b61, NODE_3T, f"{label} entry 61")
    res += check_cust_z([i32["z"] | {"log": i32["log"]}, z61, z60], label)
    res += check_cust_p(i32, b61, h60, label)
    res += check_d(b61, h60, z61["v_s"], label)
    report(res, label)
    f = figures(i32, b61, h60)
    print_figures(f, label)
    return {"i32": i32, "i36": b61, "heat": h60, "fig": f}


def _i88_mode():
    raw = os.environ.get(I88_ENV, "").strip()
    if raw == "":
        return None
    if raw in ("1", "both"):
        return raw
    raise ValueError(f"{I88_ENV} must be unset, '1' or 'both', got {raw!r}")


def _i88_body(mode):
    """The gated body; returns the summary dict.  Writes only at ``both``."""
    s55 = read_1p5t()
    hits = {
        k: _on_root(_find, g)
        for k, g in (
            ("60", ENTRY60_GLOB),
            ("61", ENTRY61_GLOB),
            ("63", ENTRY63_GLOB),
            ("67", ENTRY67_GLOB),
        )
    }
    for k, h in hits.items():
        print(f"{TAG} entry {k}: " + (f"{h[0].name} ({h[1]})" if h else "pending"))
    missing = [k for k, h in hits.items() if h is None]
    if mode == "both" and missing:
        raise AssertionError(
            f"entries {missing} absent: I88 at 'both' needs 48 / 59 / 60 / 61 / 63 / 67"
        )
    s15 = s55
    if hits["63"] is not None:
        s63 = read_1p5t(hits["63"][0], gates="none", label="1.5 T (entry 63)")
        if mode == "both":
            s15 = s63
            print(f"{TAG} 1.5 T theta source: entry 63; entry 55 printed beside")
    s3 = None
    r = None
    if all(hits[k] is not None for k in ("60", "61", "67")):
        s3 = read_3t(hits["60"][0], hits["61"][0], hits["67"][0])
        r = ratios(s15["fig"], s3["fig"])
        print_ratios(r, f"3 T / 1.5 T ({s15['i32']['log']})")
    else:
        print(
            f"{TAG} 3 T / 1.5 T ratio: pending (entries {[k for k in ('60', '61', '67') if hits[k] is None]})"
        )
    out = {
        "sense": f"{CW} = cw",
        "fig_1p5t": s15["fig"],
        "fig_entry55": s55["fig"],
        "fig_3t": s3["fig"] if s3 else None,
        "ratios": {k: f"{v:.6e}" for k, v in r.items()} if r else None,
        "logs": {k: (h[0].name if h else None) for k, h in hits.items()},
    }
    if mode == "both":

        def _write():
            HOST_JSON_DIR.mkdir(parents=True, exist_ok=True)
            p = HOST_JSON_DIR / I88_JSON_NAME
            p.write_text(json.dumps(out, indent=1))
            return str(p)

        print(f"{TAG} wrote {_on_root(_write)}")
    return out


# ---------------------------------------------------------------------------
# Coding tests (no knob).
# ---------------------------------------------------------------------------
def test_i88_custody_1p5t():
    """(X48) (X55) (X59), (Cust-Z), (Cust-P), (D-Pacc), (D-f) and the 1.5 T pre-figure."""
    s = read_1p5t()
    f = s["fig"]
    assert list(f["twin"]) == list(BALLS_1P5) and list(f["rod"]) == list(BALLS_1P5)
    print(
        f"{TAG} predicted (drafter, unlogged): P_acc_rod / P_acc_twin - 1 ~ -5.8e-6; twin ~ "
        f"1.2839e-03 / 1.2874e-03 K/W; rod ~ 1.3981e-03 / 1.4030e-03 K/W (top / bottom); printed"
    )


def test_i88_stub_3t():
    """(Stub-3T): the 3 T side's custody on C32 / C31 / I85 stubs — stub, form proof, not a record."""
    i85 = _parse_i32(I85_STUB_LOG, I85_STUB_LINE)
    b32 = _parse_i36(C32_STUB_LOG, C32_STUB_LINE)
    h31 = _on_root(parse_p61_cw_sources, C31_STUB_LOG, C31_CW_LINE)
    z32 = _parse_z(C32_STUB_LOG, C32_STUB_LINE)
    z31 = _parse_z(C31_STUB_LOG, C31_Z_LINE)
    s2_31 = _on_root(parse_c6a_s2, C31_STUB_LOG, C31_S2_LINE)
    name = "3 T stub (C32 / C31 / I85)"
    res = check_cust_z([z32, i85["z"] | {"log": i85["log"]}, z31], name)
    res.append(
        (
            f"{name} (Z) line numbers C32 :81, I85 :81, C31 :87",
            (z32["line"], i85["z"]["line"], z31["line"]) == (C32_STUB_LINE, I85_Z_LINE, C31_Z_LINE),
            f"{(z32['line'], i85['z']['line'], z31['line'])}",
        )
    )
    s2r = {k: b32["arms"]["R"]["s2"][k] for k in ("s11", "s12", "s21", "s22")}
    s2c = {k: s2_31[k] for k in ("s11", "s12", "s21", "s22")}
    res.append((f"{name} C32 arm R S2 == C31 :91 S2 (string)", s2r == s2c, f"{s2r} vs {s2c}"))
    cw = h31["senses"]["cw"]["p61"]["src"]
    p = p_acc_of(b32["arms"]["R"]["s2"], z32["v_s"])
    err = abs(p - float(cw["p_acc"])) / abs(float(cw["p_acc"]))
    res.append(
        (
            f"{name} (D-Pacc) _p_acc_spec vs C31 :{cw['line']} cw P_acc",
            err <= DPACC_REL and cw["line"] == C31_CW_LINE,
            f"{p:.12e} vs {cw['p_acc']}: rel {err:.3e} (<= 1e-08)",
        )
    )
    f = figures(i85, b32, h31)
    res.append(
        (
            f"{name} P_acc_rod / P_acc_twin - 1 == 0.0 exactly",
            f["p_acc_rod_over_twin_m1"] == 0.0,
            f"{f['p_acc_rod_over_twin_m1']!r}",
        )
    )
    sh, am = i85["share"], b32["arms"]["M"]
    got = (sh["share"], sh["p_vol"], sh["gel"], sh["rod"])
    want = (am["share"], am["p_vol"], am["gel"], am["rod"])
    res.append(
        (
            f"{name} I85 :88 == C32 :85",
            got == want
            and sh["line"] == I85_SHARE_LINE
            and am["y_in"]["line"] == C32_ARM_M_YIN_LINE,
            f"{got} vs {want}",
        )
    )
    tp, ar = i85["arms"]["twin"]["src"], b32["arms"]["R"]
    d = abs(float(tp["p_gel"]) - float(ar["p_vol"]))
    res.append(
        (
            f"{name} I85 :91 twin P_gel vs C32 :89 arm R P_vol within half a unit",
            d <= _half_ulp_abs(ar["p_vol"])
            and tp["line"] == I85_TWIN_PGEL_LINE
            and ar["y_in"]["line"] == C32_ARM_R_YIN_LINE,
            f"{tp['p_gel']} vs {ar['p_vol']}: |d| {d:.3e} (<= {_half_ulp_abs(ar['p_vol']):.3e})",
        )
    )
    res.append(
        (
            f"{name} I85 gates failures none",
            i85["gates"]["failures"] == "none",
            i85["gates"]["failures"],
        )
    )
    report(res, name)
    print_figures(f, "3 T stub, form proof, not a record")
    s15 = read_1p5t()
    print_ratios(
        ratios(
            {
                "twin": {"top": f["twin"]["a"], "bottom": f["twin"]["b"]},
                "rod": {"top": f["rod"]["a"], "bottom": f["rod"]["b"]},
            },
            s15["fig"],
        ),
        "stub, form proof, not a record (stub a / b against 1.5 T top / bottom)",
    )
    print(
        f"{TAG} (Cust-P) at 3 T is not stub-provable (C31's P61 stub feeds a different box "
        f"field); its first real check is entry 60 against entry 67 at 'both'"
    )


def test_i88_form63():
    """(Form-63): ``parse_i32_block`` reads entry 63's form (I32b L-W0b)."""
    b = _parse_i32(I32B_FORM_LOG, I32B_FORM_LINE)
    b2 = _parse_i32(I32B_FORM_LOG)
    assert b == b2
    lines = _on_root(lambda: _read_log(I32B_FORM_LOG).splitlines())
    for arm, want_lines in (("rod", I32B_ROD_THETA_LINES), ("twin", I32B_TWIN_THETA_LINES)):
        balls = b["arms"][arm]["balls"]
        got = [(v["line"], v["theta"]) for v in balls.values()]
        want = [(ln, _RX_THETA.search(lines[ln - 1]).group(5)) for ln in want_lines]
        print(f"{TAG} (Form-63) {arm}: parsed {got} vs lines {want}")
        assert got == want, (arm, got, want)
    assert b["gates"]["failures"] == "none"


def test_i88_self():
    """(Self): the 1.5 T side against itself gives every ratio exactly 1."""
    s = read_1p5t()
    r = ratios(s["fig"], copy.deepcopy(s["fig"]))
    print_ratios(r, "(Self) 1.5 T / 1.5 T")
    assert len(r) == 4 and all(v == 1.0 for v in r.values()), r
    assert all(f"{v:.6e}" == "1.000000e+00" for v in r.values())


def _empty_log_roots(monkeypatch, tmp_path):
    """Point the pending-entry roots at two empty tmp dirs (§9 item 918)."""

    def _make():
        t, h = Path(tmp_path) / "tracked-empty", Path(tmp_path) / "held-empty"
        t.mkdir(parents=True, exist_ok=True)
        h.mkdir(parents=True, exist_ok=True)
        return str(t), str(h)

    t, h = _on_root(_make)
    monkeypatch.setitem(globals(), "TRACKED_LOG_DIR", Path(t))
    monkeypatch.setitem(globals(), "HELD_LOG_DIR", Path(h))


def _json_stamp(p):
    return _on_root(lambda: (p.stat().st_mtime_ns, p.read_bytes()) if p.exists() else None)


def test_i88_pend(monkeypatch, tmp_path, capsys):
    """(Pend): entries 60 / 61 / 63 / 67 absent ⇒ "pending", nothing written."""
    _empty_log_roots(monkeypatch, tmp_path)
    before = _json_stamp(HOST_JSON_DIR / I88_JSON_NAME)
    out = _i88_body("1")
    assert out["ratios"] is None and out["fig_3t"] is None
    assert all(v is None for v in out["logs"].values())
    assert _json_stamp(HOST_JSON_DIR / I88_JSON_NAME) == before
    txt = capsys.readouterr().out
    sys.stdout.write(txt)
    for k in ("60", "61", "63", "67"):
        assert f"{TAG} entry {k}: pending" in txt


def test_i88_n_both(monkeypatch, tmp_path):
    """At ``both`` on empty roots the body fails naming the entries; nothing written."""
    _empty_log_roots(monkeypatch, tmp_path)
    before = _json_stamp(HOST_JSON_DIR / I88_JSON_NAME)
    with pytest.raises(AssertionError, match=r"entries \['60', '61', '63', '67'\] absent"):
        _i88_body("both")
    assert _json_stamp(HOST_JSON_DIR / I88_JSON_NAME) == before


def test_i88_knob(monkeypatch):
    monkeypatch.setenv(I88_ENV, "bogus")
    with pytest.raises(ValueError):
        _i88_mode()
    monkeypatch.setenv(I88_ENV, "both")
    assert _i88_mode() == "both"
    monkeypatch.delenv(I88_ENV)
    assert _i88_mode() is None


def test_i88_n_sense():
    """(N-sense): entry 55's twin θ with entry 59's **ccw** factor; separation fixed by
    two logged factors, asserted above their combined half-unit bound."""
    h59 = _parse_heat(ENTRY59_LOG)
    i32 = _parse_i32(ENTRY55_LOG)
    fc, fcc = (
        h59["senses"]["cw"]["p61"]["src"]["factor"],
        h59["senses"]["ccw"]["p61"]["src"]["factor"],
    )
    sep = float(fcc) / float(fc) - 1.0
    bound = _half_ulp_rel(fc) + _half_ulp_rel(fcc)
    for b, v in i32["arms"]["twin"]["balls"].items():
        print(
            f"{TAG} (N-sense) twin {b}: theta x ccw factor {float(v['theta']) * float(fcc):.6e} vs "
            f"x cw factor {float(v['theta']) * float(fc):.6e} K/W"
        )
    print(f"{TAG} (N-sense) separation ccw/cw - 1 = {sep:.6f} (> bound {bound:.3e})")
    assert sep > bound
    assert f"{sep:.6f}" == "0.268866"


def test_i88_n_drive():
    """(N-drive): (Cust-Z) on C32 :81's 3 T tokens against entry 48's fails."""
    z32 = _parse_z(C32_STUB_LOG, C32_STUB_LINE)
    z48 = _parse_z(ENTRY48_LOG, tag="I36a")
    res = check_cust_z([z32, z48], "(N-drive)")
    print(f"{TAG} (N-drive) {res}")
    assert not res[0][1]
    with pytest.raises(AssertionError, match="Cust-Z"):
        report(res, "(N-drive)")


def test_i88_n_digit():
    """(N-digit): a one-digit in-memory change to entry 55's :12499 token fails (Cust-P)."""
    i32 = _parse_i32(ENTRY55_LOG)
    b48 = _parse_i36(ENTRY48_LOG)
    h59 = _parse_heat(ENTRY59_LOG)
    assert all(ok for _, ok, _ in check_cust_p(i32, b48, h59, "clean"))
    bad = copy.deepcopy(i32)
    tok = bad["arms"]["twin"]["src"]["p_gel"]
    mant, e = tok.split("e")
    mant = mant[:-1] + str((int(mant[-1]) + 1) % 10)
    bad["arms"]["twin"]["src"]["p_gel"] = mant + "e" + e
    res = check_cust_p(bad, b48, h59, "(N-digit)")
    failed = [n for n, ok, _ in res if not ok]
    print(f"{TAG} (N-digit) {tok} -> {bad['arms']['twin']['src']['p_gel']}: failed {failed}")
    assert failed == ["(N-digit) (Cust-P) I32 twin P_gel == P61 cw P_gel (string)"], failed


def test_i88_n_blocks():
    """(N-blocks): I85's W0b (five I32 headers) without ``at_line`` raises."""
    with pytest.raises(ValueError, match="I32 blocks"):
        _parse_i32(I85_STUB_LOG)


# ---------------------------------------------------------------------------
# The gated test.
# ---------------------------------------------------------------------------
def test_geo36_i88_reader():
    mode = _i88_mode()
    if mode is None:
        pytest.skip(f"{I88_ENV} unset")
    before = _on_root(
        lambda: sorted(p.name for p in HOST_JSON_DIR.glob("*")) if HOST_JSON_DIR.is_dir() else []
    )
    _i88_body(mode)
    if mode == "1":
        after = _on_root(
            lambda: (
                sorted(p.name for p in HOST_JSON_DIR.glob("*")) if HOST_JSON_DIR.is_dir() else []
            )
        )
        assert after == before, "wrote under logs/geo36 at knob 1"
        print(f"{TAG} knob 1: nothing written under logs/geo36")
