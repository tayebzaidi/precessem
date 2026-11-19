"""`GEO-38` C46 — C39's 1.5 T / 3 T thermal reader, pre-built on P61b's and C31's
stub logs (§9 item 907; docs/planning/chunks/GEO-38.md "§9 item 907").  No solve,
no mesh.

C39 reads entries 59 / 60: θ at the hotspot ball and P37's six points at t* per
accepted watt under both senses, ``cem43`` and the six-minute running mean at the
hotspot, the 3 T / 1.5 T ratio of each.  **Correction of spelling (for the weekly):**
neither P61's module (entry 59's) nor C31's (entry 60's) writes a JSON — their results
are print lines only — so this reader parses the window logs' ``[GEO-38 P61]`` /
``[GEO-36 I32] [gel twin <s>]`` lines.

Anchors: (Stub-59) / (Stub-60) ``check_block`` green on P61b W0 :146 and C31 W0r :130
under both senses — (D-f) factor = P_gel / P_acc, (D-scale) P61's probe θ = I32's θ ×
factor, (Cem) t*·R^(43−T_INF) ≤ cem43 ≤ t*·R^(43−T_INF−θ_point), (RA) running-mean
max ≤ θ_point at t = t*, (T*), (Node); (Stub-dry) C39's table on the two stubs; (Pend)
entries 59 / 60 absent.  Controls: (N-cem), (N-scale), (N-blocks), (N-sense).

Knob ``FEM_EM_GEO38_C39``: unset skip; ``1`` each present entry read and checked, an
absent one prints "pending", nothing written; ``both`` both required, rank 0 writes
``/workspace/logs/geo38/c39-thermal-1p5t-3t.json`` after every check holds.  No log
token of >= 8 significant digits is a literal here: every number checked is parsed.
"""

from __future__ import annotations

import gzip
import inspect
import json
import math
import os
import re
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from fem_em_solver.thermal.dose import cem43  # noqa: E402
from tests.validation.test_geo38_f_clinical_box_heating import (  # noqa: E402
    RUNNING_WINDOW_S,
    SENSE_TUPLE,
    T_STAR_BOX,
    _probes,
)
from tests.validation.test_mat15_coil_driven_heating_phantom import T_INF  # noqa: E402

C39_ENV = "FEM_EM_GEO38_C39"
#: Entry 59's XL_CHUNK (cannot match ``…-P33-3T.log``).
ENTRY59_GLOB = "*_GEO-38-P33.log"
#: Entry 60's XL_CHUNK as pre-registered; its READY review writes the final one.
ENTRY60_GLOB = "*_GEO-38-P33-3T.log"
ENTRY59_NODE_PF = "193.5996"
ENTRY60_NODE_PF = "44.0525"
TRACKED_LOG_DIR = REPO / "docs" / "testing" / "logs"
HELD_LOG_DIR = REPO / "logs" / "leak-held"
P61B_STUB_LOG = TRACKED_LOG_DIR / "20261105T103249Z_GEO-38-P61b-W0.log"
P61B_STUB_LINE = 146
P61B_STUB_NODE_PF = "193.5996"
P61B_CW_ONLY_LINE = 97
C31_STUB_LOG = TRACKED_LOG_DIR / "20261105T181050Z_GEO-38-C31-W0r.log"
C31_STUB_LINE = 130
C31_STUB_NODE_PF = "44.0615"
HOST_JSON_DIR = Path("/workspace/logs/geo38")
C39_JSON_NAME = "c39-thermal-1p5t-3t.json"
SENSES = ("ccw", "cw")
TAG = "[GEO-38 C46]"
#: cem43's own R below t_ref and t_ref, read from its signature (dose.py:45).
_SIG = inspect.signature(cem43).parameters
R_BELOW = float(_SIG["r_below"].default)
T_REF = float(_SIG["t_ref"].default)
_TUPLE_SENSE = {v: k for k, v in SENSE_TUPLE.items()}


def _comm():
    from mpi4py import MPI

    return MPI.COMM_WORLD


def _on_root(fn, *args, **kwargs):
    """Run ``fn`` on rank 0 only, bcast its result or its exception to every rank."""
    comm = _comm()
    out = None
    if comm.rank == 0:
        try:
            out = ("ok", fn(*args, **kwargs))
        except BaseException as exc:  # re-raised on every rank below
            out = ("err", exc)
    out = comm.bcast(out, root=0)
    if out[0] == "err":
        raise out[1]
    return out[1]


# ---------------------------------------------------------------------------
# Pure readers.
# ---------------------------------------------------------------------------
def _read_log(path) -> str:
    """A harness log's text, from ``path`` or, once retention has gzipped it,
    ``path`` + ``.gz``."""
    path = str(path)
    if not os.path.exists(path) and os.path.exists(f"{path}.gz"):
        path = f"{path}.gz"
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt", errors="replace") as fh:
        return fh.read()


def find_log(glob, dirs=None):
    """``(path, "tracked" | "held")`` of the one log matching ``glob`` (either suffix)
    under ``docs/testing/logs/`` then ``logs/leak-held/`` (read only), or ``None``.
    Two matches raise ``ValueError``.  (Item 905's spelling, own copy.)"""
    if dirs is None:
        dirs = ((TRACKED_LOG_DIR, "tracked"), (HELD_LOG_DIR, "held"))
    hits = []
    for d, kind in dirs:
        d = Path(d)
        if not d.is_dir():
            continue
        for p in sorted(set(d.glob(glob)) | set(d.glob(glob + ".gz"))):
            hits.append((p, kind))
    if len(hits) > 1:
        raise ValueError(f"{glob}: {len(hits)} matches {[str(p) for p, _ in hits]}")
    if not hits:
        return None
    p, kind = hits[0]
    return (Path(str(p).removesuffix(".gz")), kind)


_T = r"(\(1,[+-]1j\))"
_RX_HEAD = re.compile(
    r"^start: -n (\d+), P61 body; node C ([0-9.]+) pF, degree (\d); "
    r".*?thermal senses \[([^\]]*)\]"
)
_RX_GATES = re.compile(r"\[GEO-38 P61\] gates: failures")
_RX_I32_SRC = re.compile(
    r"\[GEO-36 I32\] \[gel twin (\w+)\] sense " + _T + r": P_gel (\S+) W; G (\S+) 1/W; G\*P (\S+) W"
)
_RX_I32_BALL = re.compile(
    r"\[GEO-36 I32\] \[gel twin (\w+)\] 1 g ball \(a (\S+) mm\) at (.+?), "
    r"sense " + _T + r": theta\(t\*\) (\S+) K"
)
_RX_I32_PROBE = re.compile(
    r"\[GEO-36 I32\] \[gel twin (\w+)\] probe (.+?) \(([^)]*)\), sense "
    + _T
    + r": theta\(t\*\) (\S+) K \(valid (\w+)"
)
_RX_I32_CN = re.compile(r"\[GEO-36 I32\] \[gel twin (\w+)\] CN N=(\d+) to t\* (\S+) s")
_RX_P61_SRC = re.compile(
    r"\[GEO-38 P61\] sources, sense " + _T + r" \((\w+), heated\): "
    r"P_gel (\S+) W, P_acc (\S+) W .*?factor P_gel/P_acc \((\w+)\) (\S+)\s*$"
)
_RX_P61_HOT = re.compile(
    r"\[GEO-38 P61\] hotspot (.+?) (\([^)]*\)), sense "
    + _T
    + r": theta\(t\*\) per P_acc watt (\S+) K \(point; final-field read (\S+) K\); "
    r"cem43 per P_acc watt (\S+) s .*?running mean \(W (\S+) s\) max (\S+) K "
    r"at t (\S+) s"
)
_RX_P61_PROBE = re.compile(
    r"\[GEO-38 P61\] probe (.+?), sense " + _T + r": theta\(t\*\) per P_acc watt (\S+) K"
)


def _senses_of(raw):
    return [s.strip().strip("'\"") for s in raw.split(",") if s.strip()]


def parse_heat_log(path, at_line=None) -> dict:
    """Pure: one P61 block (its ``start: … P61 body`` header to the next ``gates:``
    line), every token kept as its printed **string**.  A window log must hold exactly
    one block; a stub selects its block by ``at_line`` (the header's 1-based line).
    The block must carry both senses ``['ccw', 'cw']``.  ``ValueError`` otherwise, and
    on any missing / repeated / sense-inconsistent line."""
    lines = _read_log(path).splitlines()
    name = Path(str(path)).name.removesuffix(".gz")
    heads = [i for i, ln in enumerate(lines, start=1) if _RX_HEAD.search(ln)]
    if at_line is None:
        if len(heads) != 1:
            raise ValueError(
                f"{name}: {len(heads)} P61 blocks (lines {heads[:8]}); a window "
                f"log holds exactly one"
            )
        start = heads[0]
    else:
        if at_line not in heads:
            raise ValueError(f"{name}: no P61 header at :{at_line}")
        start = at_line
    h = _RX_HEAD.search(lines[start - 1])
    senses = _senses_of(h.group(4))
    if senses != list(SENSES):
        raise ValueError(
            f"{name}: block at :{start} has senses {senses}; both {list(SENSES)} " f"required"
        )
    end = next((i for i in range(start, len(lines) + 1) if _RX_GATES.search(lines[i - 1])), None)
    if end is None:
        raise ValueError(f"{name}: block at :{start} has no gates line")
    inner = [n for n in heads if start < n < end]
    if inner:
        raise ValueError(f"{name}: block at :{start} holds another header at :{inner[0]}")
    per = {s: {"i32": {"probes": {}}, "p61": {"probes": {}}} for s in SENSES}

    def _sense(word, tup, i):
        if word not in per:
            raise ValueError(f"{name}: unknown sense {word!r} at :{i}")
        if SENSE_TUPLE[word] != tup:
            raise ValueError(f"{name}: sense {word} with tuple {tup} at :{i}")
        return per[word]

    def _put(d, key, val, i):
        if key in d:
            raise ValueError(f"{name}: {key} repeated at :{i}")
        d[key] = val

    for i in range(start, end + 1):
        ln = lines[i - 1]
        m = _RX_I32_SRC.search(ln)
        if m:
            s = _sense(m.group(1), m.group(2), i)
            _put(
                s["i32"],
                "src",
                {"p_gel": m.group(3), "g": m.group(4), "gp": m.group(5), "line": i},
                i,
            )
            continue
        m = _RX_I32_BALL.search(ln)
        if m:
            s = _sense(m.group(1), m.group(4), i)
            _put(
                s["i32"],
                "ball",
                {"a_mm": m.group(2), "label": m.group(3), "theta": m.group(5), "line": i},
                i,
            )
            continue
        m = _RX_I32_PROBE.search(ln)
        if m:
            s = _sense(m.group(1), m.group(4), i)
            _put(
                s["i32"]["probes"],
                m.group(2),
                {"pt": m.group(3), "theta": m.group(5), "valid": m.group(6), "line": i},
                i,
            )
            continue
        m = _RX_I32_CN.search(ln)
        if m:
            if m.group(1) not in per:
                raise ValueError(f"{name}: unknown sense {m.group(1)!r} at :{i}")
            _put(
                per[m.group(1)]["i32"], "cn", {"n": m.group(2), "t_star": m.group(3), "line": i}, i
            )
            continue
        m = _RX_P61_SRC.search(ln)
        if m:
            if m.group(2) != m.group(5):
                raise ValueError(
                    f"{name}: sources line at :{i} names {m.group(2)} and {m.group(5)}"
                )
            s = _sense(m.group(2), m.group(1), i)
            _put(
                s["p61"],
                "src",
                {"p_gel": m.group(3), "p_acc": m.group(4), "factor": m.group(6), "line": i},
                i,
            )
            continue
        m = _RX_P61_HOT.search(ln)
        if m:
            if m.group(3) not in _TUPLE_SENSE:
                raise ValueError(f"{name}: hotspot tuple {m.group(3)} at :{i}")
            s = per[_TUPLE_SENSE[m.group(3)]]
            _put(
                s["p61"],
                "hot",
                {
                    "label": m.group(1),
                    "point": m.group(2),
                    "theta_point": m.group(4),
                    "theta_final": m.group(5),
                    "cem43": m.group(6),
                    "ra_w": m.group(7),
                    "ra_max": m.group(8),
                    "ra_t": m.group(9),
                    "line": i,
                },
                i,
            )
            continue
        m = _RX_P61_PROBE.search(ln)
        if m:
            if m.group(2) not in _TUPLE_SENSE:
                raise ValueError(f"{name}: probe tuple {m.group(2)} at :{i}")
            s = per[_TUPLE_SENSE[m.group(2)]]
            _put(s["p61"]["probes"], m.group(1), {"theta": m.group(3), "line": i}, i)
    for sn, s in per.items():
        for part, keys in (("i32", ("src", "ball", "cn")), ("p61", ("src", "hot"))):
            miss = [k for k in keys if k not in s[part]]
            if miss:
                raise ValueError(f"{name}: block at :{start}, sense {sn}, lacks {part} {miss}")
        if not s["i32"]["probes"] or not s["p61"]["probes"]:
            raise ValueError(f"{name}: block at :{start}, sense {sn}, lacks probe lines")
    return {
        "log": name,
        "head": {
            "line": start,
            "ranks": h.group(1),
            "node": h.group(2),
            "degree": h.group(3),
            "senses": senses,
        },
        "end": end,
        "senses": per,
    }


def _half_ulp_rel(tok):
    """Half a unit in the printed token's last digit, relative to its value."""
    t = tok.strip()
    mant, _, exp = t.lower().partition("e")
    dec = mant.split(".")[1] if "." in mant else ""
    abs_ulp = 0.5 * 10.0 ** ((int(exp) if exp else 0) - len(dec))
    return abs_ulp / abs(float(t))


def cem_bracket(t_star_tok, theta_tok, t_inf=T_INF):
    """``(lower, upper)``: t*·R^(t_ref − t_inf) and t*·R^(t_ref − t_inf − θ)."""
    t, th = float(t_star_tok), float(theta_tok)
    return (t * R_BELOW ** (T_REF - t_inf), t * R_BELOW ** (T_REF - t_inf - th))


def check_block(block, node_pf, probe_names=None, name=""):
    """Derived checks on one block, both senses: ``[(check, ok, detail)]``.  Each
    holds at its half-unit bound (Σ half-ulp rel over the tokens involved)."""
    res = []

    def rel_check(tag, val, tok, inputs):
        bound = _half_ulp_rel(tok) + sum(_half_ulp_rel(t) for t in inputs)
        err = abs(val - float(tok)) / abs(float(tok))
        res.append(
            (
                f"{name} {tag}",
                err <= bound,
                f"{val:.12e} vs printed {tok}: rel {err:.3e} (<= {bound:.3e})",
            )
        )

    res.append(
        (
            f"{name} (Node) node token == {node_pf}",
            block["head"]["node"] == node_pf,
            f"log node {block['head']['node']} vs {node_pf}",
        )
    )
    for sn in SENSES:
        s = block["senses"][sn]
        i32, p61 = s["i32"], s["p61"]
        src, hot = p61["src"], p61["hot"]
        fac = src["factor"]
        rel_check(
            f"{sn} (D-f) factor = P_gel/P_acc",
            float(src["p_gel"]) / float(src["p_acc"]),
            fac,
            (src["p_gel"], src["p_acc"]),
        )
        res.append(
            (
                f"{name} {sn} (P_gel) P61 == I32 string",
                src["p_gel"] == i32["src"]["p_gel"],
                f"{src['p_gel']} vs {i32['src']['p_gel']}",
            )
        )
        same = list(p61["probes"]) == list(i32["probes"])
        res.append(
            (
                f"{name} {sn} (D-scale) probe names P61 == I32",
                same,
                f"{list(p61['probes'])} vs {list(i32['probes'])}",
            )
        )
        for pn, pv in p61["probes"].items():
            if pn not in i32["probes"]:
                continue
            itok = i32["probes"][pn]["theta"]
            rel_check(
                f"{sn} (D-scale) probe {pn}", float(itok) * float(fac), pv["theta"], (itok, fac)
            )
        t_tok = i32["cn"]["t_star"]
        th = float(hot["theta_point"])
        if th < T_REF - T_INF:
            lo, up = cem_bracket(t_tok, hot["theta_point"])
            c = float(hot["cem43"])
            h = _half_ulp_rel(hot["cem43"]) + _half_ulp_rel(t_tok)
            ok = lo * (1.0 - h) <= c <= up * (1.0 + h)
            res.append(
                (
                    f"{name} {sn} (Cem) bracket",
                    ok,
                    f"lower {lo:.12e} <= cem43 {hot['cem43']} <= upper {up:.12e} (slack rel "
                    f"{h:.2e}); cem/lower - 1 {c / lo - 1.0:.3e}, width up/lower - 1 "
                    f"{up / lo - 1.0:.3e}; R {R_BELOW}, t_ref {T_REF}, T_INF {T_INF}, t* {t_tok}, "
                    f"theta_point {hot['theta_point']}",
                )
            )
        else:
            print(
                f"{TAG} {name} {sn} (Cem) bracket not applicable: theta_point {th} >= "
                f"{T_REF - T_INF}"
            )
        hr = _half_ulp_rel(hot["ra_max"]) + _half_ulp_rel(hot["theta_point"])
        ra_ok = float(hot["ra_max"]) <= th * (1.0 + hr)
        res.append(
            (
                f"{name} {sn} (RA) running max <= theta_point",
                ra_ok,
                f"{hot['ra_max']} <= {hot['theta_point']} (slack rel {hr:.2e})",
            )
        )
        res.append(
            (
                f"{name} {sn} (RA) its t == t* at 0.1 s",
                hot["ra_t"] == f"{T_STAR_BOX:.1f}",
                f"t {hot['ra_t']} vs {T_STAR_BOX:.1f}",
            )
        )
        res.append(
            (
                f"{name} {sn} (RA) window == {RUNNING_WINDOW_S:.0f} s",
                hot["ra_w"] == f"{RUNNING_WINDOW_S:.0f}",
                f"W {hot['ra_w']}",
            )
        )
        res.append(
            (
                f"{name} {sn} (T*) t* == T_STAR_BOX .3f",
                t_tok == f"{T_STAR_BOX:.3f}",
                f"{t_tok} vs {T_STAR_BOX:.3f}",
            )
        )
        if probe_names is not None:
            res.append(
                (
                    f"{name} {sn} (Probes) names == p61._probes()",
                    list(p61["probes"]) == list(probe_names),
                    f"{list(p61['probes'])} vs {list(probe_names)}",
                )
            )
        else:
            print(f"{TAG} {name} {sn} (Probes) stub block's own probes: {list(p61['probes'])}")
        rel = abs(float(hot["theta_final"]) / th - 1.0)
        print(
            f"{TAG} {name} {sn} point read {hot['theta_point']} | final-field read "
            f"{hot['theta_final']} K: rel {rel:.3e} (printed, never asserted)"
        )
    return res


def require_block(block, node_pf, probe_names=None, name=""):
    res = check_block(block, node_pf, probe_names, name)
    for n, ok, d in res:
        print(f"{TAG} {'PASS' if ok else 'FAIL'} {n}: {d}")
    bad = [(n, d) for n, ok, d in res if not ok]
    if bad:
        raise AssertionError(f"{name}: {bad}")
    return res


def side(block):
    """C39's readout per sense, as floats (derived quantities labelled)."""
    out = {}
    for sn in SENSES:
        s = block["senses"][sn]
        fac = float(s["p61"]["src"]["factor"])
        hot = s["p61"]["hot"]
        lo, _up = cem_bracket(s["i32"]["cn"]["t_star"], hot["theta_point"])
        d = {
            "ball theta per P_acc (derived = I32 ball x factor)": float(s["i32"]["ball"]["theta"])
            * fac,
            "hotspot point theta per P_acc": float(hot["theta_point"]),
        }
        for pn, pv in s["p61"]["probes"].items():
            d[f"probe {pn}"] = float(pv["theta"])
        d["cem43 per P_acc"] = float(hot["cem43"])
        d["cem43 excess over t*.R^(43-T_INF)"] = float(hot["cem43"]) - lo
        d[f"running mean ({RUNNING_WINDOW_S:.0f} s) max"] = float(hot["ra_max"])
        out[sn] = d
    return out


def ratios(s15, s3):
    """3 T / 1.5 T ratio of each readout, per sense (floats; print at ``.6e``).  The
    keys must match, otherwise ``ValueError``."""
    out = {}
    for sn in SENSES:
        if list(s15[sn]) != list(s3[sn]):
            raise ValueError(f"{sn}: readout keys differ: {list(s15[sn])} vs {list(s3[sn])}")
        out[sn] = {k: s3[sn][k] / s15[sn][k] for k in s15[sn]}
    return out


def print_table(b15, b3=None, stub=False):
    """C39's table; returns ``(rows, ratios)``.  ``rows``: the printed lines."""
    pend = "pending (xl entry 60)"
    mark = "STUB " if stub else ""
    s15 = side(b15)
    s3 = side(b3) if b3 is not None else None
    r = ratios(s15, s3) if s3 is not None else None
    rows = [
        f"{TAG} {mark}C39 table: 1.5 T ({b15['log']}:{b15['head']['line']}, node "
        f"{b15['head']['node']} pF) | 3 T "
        f"({b3['log'] + ':' + str(b3['head']['line']) + ', node ' + b3['head']['node'] + ' pF' if b3 else pend})"
        f" | ratio 3T/1.5T (.6e)"
    ]
    for sn in SENSES:
        for k, v in s15[sn].items():
            v3 = format(s3[sn][k], ".9e") if s3 else pend
            rr = format(r[sn][k], ".6e") if r else pend
            rows.append(f"{TAG} {mark}  {sn} {k}: {v:.9e} | {v3} | {rr}")
    for ln in rows:
        print(ln)
    return rows, r


def _c39_mode():
    raw = os.environ.get(C39_ENV, "").strip()
    if raw == "":
        return None
    if raw in ("1", "both"):
        return raw
    raise ValueError(f"{C39_ENV} must be unset, '1' or 'both', got {raw!r}")


def _parse(path, at_line=None):
    return _on_root(parse_heat_log, path, at_line)


def _c39_body(mode):
    """The gated body; returns the summary dict.  Writes only at ``both``."""
    names = list(_probes())
    got = {}
    for key, glob, node, label in (
        ("59", ENTRY59_GLOB, ENTRY59_NODE_PF, "1.5 T"),
        ("60", ENTRY60_GLOB, ENTRY60_NODE_PF, "3 T"),
    ):
        hit = _on_root(find_log, glob)
        if hit is None:
            print(f"{TAG} entry {key} ({label}): pending (no {glob} tracked or held)")
            got[key] = None
            continue
        path, kind = hit
        print(f"{TAG} entry {key} ({label}): {path.name} ({kind})")
        b = _parse(path)
        require_block(b, node, names, f"{label} (entry {key})")
        got[key] = b
    if mode == "both":
        missing = [k for k, v in got.items() if v is None]
        if missing:
            raise AssertionError(f"entries {missing} absent: C39 needs both entry 59 and entry 60")
    r = None
    if got["59"] is not None:
        _rows, r = print_table(got["59"], got["60"])
    out = {
        "entry59": got["59"],
        "entry60": got["60"],
        "ratios": (
            {sn: {k: f"{v:.6e}" for k, v in d.items()} for sn, d in r.items()} if r else None
        ),
    }
    if mode == "both":

        def _write():
            HOST_JSON_DIR.mkdir(parents=True, exist_ok=True)
            p = HOST_JSON_DIR / C39_JSON_NAME
            p.write_text(json.dumps(out, indent=1))
            return str(p)

        print(f"{TAG} wrote {_on_root(_write)}")
    return out


# ---------------------------------------------------------------------------
# Coding tests (no knob).
# ---------------------------------------------------------------------------
def _tmp_dir(tmp_path):
    return _on_root(lambda: str(tmp_path))


def _stubs():
    return (_parse(P61B_STUB_LOG, P61B_STUB_LINE), _parse(C31_STUB_LOG, C31_STUB_LINE))


def _print_worst(res, label):
    for kind in ("(D-f)", "(D-scale)"):
        rows = [d for n, _ok, d in res if kind in n and "names" not in n]
        print(f"{TAG} {label} {kind}: {len(rows)} checks; " + "; ".join(rows))


def test_c46_stub59():
    """(Stub-59): P61b W0 :146 (node 193.5996), both senses."""
    b, _ = _stubs()
    assert b["head"]["senses"] == list(SENSES) and b["head"]["degree"] == "1"
    res = require_block(b, P61B_STUB_NODE_PF, None, "1.5 T stub (P61b W0 :146)")
    _print_worst(res, "Stub-59")
    for sn in SENSES:
        assert len(b["senses"][sn]["p61"]["probes"]) >= 1
        assert any(f"{sn} (Cem)" in n for n, _o, _d in res), "(Cem) not applied"


def test_c46_stub60():
    """(Stub-60): C31 W0r :130 (node 44.0615), both senses."""
    _, b = _stubs()
    assert b["head"]["senses"] == list(SENSES) and b["head"]["degree"] == "1"
    res = require_block(b, C31_STUB_NODE_PF, None, "3 T stub (C31 W0r :130)")
    _print_worst(res, "Stub-60")
    for sn in SENSES:
        assert any(f"{sn} (Cem)" in n for n, _o, _d in res), "(Cem) not applied"


def test_c46_stub_dry():
    """(Stub-dry): C39's table on the two stub blocks, "STUB" on every line, every cell
    non-empty, each ratio the quotient of its parsed values."""
    b15, b3 = _stubs()
    rows, r = print_table(b15, b3, stub=True)
    assert all("STUB" in ln for ln in rows), rows
    s15, s3 = side(b15), side(b3)
    n_cells = 0
    for ln in rows[1:]:
        cells = ln.split(": ", 1)[1].split(" | ")
        assert len(cells) == 3 and all(c.strip() for c in cells), ln
        n_cells += 1
    assert n_cells == sum(len(s15[sn]) for sn in SENSES)
    for sn in SENSES:
        assert s15[sn] and list(s15[sn]) == list(s3[sn])
        for k in s15[sn]:
            assert math.isfinite(s15[sn][k]) and s15[sn][k] != 0.0, (sn, k)
            assert r[sn][k] == s3[sn][k] / s15[sn][k], (sn, k)
        hot15, hot3 = b15["senses"][sn]["p61"]["hot"], b3["senses"][sn]["p61"]["hot"]
        assert r[sn]["hotspot point theta per P_acc"] == (
            float(hot3["theta_point"]) / float(hot15["theta_point"])
        )
        k_ball = "ball theta per P_acc (derived = I32 ball x factor)"
        q = (
            float(b3["senses"][sn]["i32"]["ball"]["theta"])
            * float(b3["senses"][sn]["p61"]["src"]["factor"])
        ) / (
            float(b15["senses"][sn]["i32"]["ball"]["theta"])
            * float(b15["senses"][sn]["p61"]["src"]["factor"])
        )
        assert r[sn][k_ball] == q


def test_c46_n_cem():
    """(N-cem): the bracket at T_INF + 1 K excludes each stub's cem43; the lower bound
    sits a factor <= 1/R = 4 above it."""
    for b in _stubs():
        for sn in SENSES:
            s = b["senses"][sn]
            hot = s["p61"]["hot"]
            c = float(hot["cem43"])
            lo1, _ = cem_bracket(s["i32"]["cn"]["t_star"], hot["theta_point"], T_INF + 1.0)
            h = _half_ulp_rel(hot["cem43"]) + _half_ulp_rel(s["i32"]["cn"]["t_star"])
            ratio = lo1 / c
            print(
                f"{TAG} (N-cem) {b['log']}:{b['head']['line']} {sn}: lower(T_INF+1) {lo1:.9e} vs "
                f"cem43 {hot['cem43']}: ratio {ratio:.6f} (excluded > 1, ceiling 1/R "
                f"{1.0 / R_BELOW:.1f})"
            )
            assert lo1 * (1.0 - h) > c, (lo1, c)
            assert 1.0 < ratio <= (1.0 / R_BELOW) * (1.0 + h), ratio


def test_c46_n_scale(tmp_path):
    """(N-scale): one digit of a P61 probe token changed in a tmp copy makes (D-scale)
    red on that probe, and only there."""
    d = Path(_tmp_dir(tmp_path)) / "nscale"
    b = _parse(C31_STUB_LOG, C31_STUB_LINE)
    pn = next(iter(b["senses"]["ccw"]["p61"]["probes"]))
    pv = b["senses"]["ccw"]["p61"]["probes"][pn]

    def _make():
        lines = _read_log(C31_STUB_LOG).split("\n")
        i = pv["line"] - 1
        tok = pv["theta"]
        j = tok.index(".") + 4
        new = tok[:j] + str((int(tok[j]) + 1) % 10) + tok[j + 1 :]
        assert f"per P_acc watt {tok} K" in lines[i]
        lines[i] = lines[i].replace(f"per P_acc watt {tok} K", f"per P_acc watt {new} K", 1)
        d.mkdir(parents=True, exist_ok=True)
        out = d / C31_STUB_LOG.name
        out.write_text("\n".join(lines))
        return str(out)

    path = _on_root(_make)
    bb = _parse(path, C31_STUB_LINE)
    res = check_block(bb, C31_STUB_NODE_PF, None, "N-scale")
    failed = {n: dd for n, ok, dd in res if not ok}
    print(f"{TAG} (N-scale) probe {pn!r} ccw token changed; failed: {failed}")
    assert list(failed) == [f"N-scale ccw (D-scale) probe {pn}"], failed


def test_c46_n_blocks(tmp_path):
    """(N-blocks): P61b W0 without ``at_line`` raises (a ['cw'] block at :97 and a
    both-sense block at :146); a tmp log holding two both-sense blocks also raises."""
    with pytest.raises(ValueError, match="P61 blocks"):
        _parse(P61B_STUB_LOG)
    d = Path(_tmp_dir(tmp_path)) / "nblocks"
    b = _parse(P61B_STUB_LOG, P61B_STUB_LINE)

    def _make():
        lines = _read_log(P61B_STUB_LOG).split("\n")
        blk = lines[b["head"]["line"] - 1 : b["end"]]
        d.mkdir(parents=True, exist_ok=True)
        out = d / "20990101T000000Z_GEO-38-P33.log"
        out.write_text("\n".join(blk + ["", "between"] + blk))
        return str(out)

    path = _on_root(_make)
    with pytest.raises(ValueError, match="2 P61 blocks"):
        _parse(path)

    # one both-sense block alone is a valid window log (the window rule's positive arm)
    def _one():
        lines = _read_log(P61B_STUB_LOG).split("\n")
        out = d / "one" / "20990101T000000Z_GEO-38-P33.log"
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text("\n".join(lines[b["head"]["line"] - 1 : b["end"]]))
        return str(out)

    one = _parse(_on_root(_one))
    assert _no_lines(one["senses"]) == _no_lines(b["senses"])


def _no_lines(d):
    """``d`` with every ``line`` key dropped (a block's tokens, position-free)."""
    if isinstance(d, dict):
        return {k: _no_lines(v) for k, v in d.items() if k != "line"}
    return d


def test_c46_n_sense():
    """(N-sense): ``at_line=97`` (senses ['cw']) raises: both senses required."""
    with pytest.raises(ValueError, match="both"):
        _parse(P61B_STUB_LOG, P61B_CW_ONLY_LINE)


def _empty_log_roots(monkeypatch, tmp_path):
    """Point ``find_log``'s default roots at two empty tmp dirs (rank 0 makes them,
    bcast), so the premise test reads entries 59 / 60 absent by construction
    (§9 item 918)."""

    def _make():
        t, h = Path(tmp_path) / "tracked-empty", Path(tmp_path) / "held-empty"
        t.mkdir(parents=True, exist_ok=True)
        h.mkdir(parents=True, exist_ok=True)
        return str(t), str(h)

    t, h = _on_root(_make)
    monkeypatch.setitem(globals(), "TRACKED_LOG_DIR", Path(t))
    monkeypatch.setitem(globals(), "HELD_LOG_DIR", Path(h))


def test_c46_pend_and_knob(monkeypatch, tmp_path):
    """(Pend): entries 59 / 60 absent on empty log roots; at ``1`` two "pending",
    nothing written; at ``both`` the body fails naming the entries; the knob's shape."""
    _empty_log_roots(monkeypatch, tmp_path)
    assert _on_root(find_log, ENTRY59_GLOB) is None
    assert _on_root(find_log, ENTRY60_GLOB) is None
    before = _on_root(
        lambda: HOST_JSON_DIR.is_dir() and sorted(p.name for p in HOST_JSON_DIR.iterdir())
    )
    # Rank 0's (mtime_ns, bytes) of C39's JSON, or None: "nothing written" holds whether
    # or not C39's own window has already written it (§9 item 917).
    stamp = lambda p: (p.stat().st_mtime_ns, p.read_bytes()) if p.exists() else None  # noqa: E731
    json_before = _on_root(stamp, HOST_JSON_DIR / C39_JSON_NAME)
    out = _c39_body("1")
    assert out["entry59"] is None and out["entry60"] is None and out["ratios"] is None
    with pytest.raises(AssertionError, match="entry 59 and entry 60"):
        _c39_body("both")
    after = _on_root(
        lambda: HOST_JSON_DIR.is_dir() and sorted(p.name for p in HOST_JSON_DIR.iterdir())
    )
    assert after == before
    assert _on_root(stamp, HOST_JSON_DIR / C39_JSON_NAME) == json_before
    monkeypatch.setenv(C39_ENV, "bogus")
    with pytest.raises(ValueError):
        _c39_mode()


def test_c46_probe_names():
    """P61's ``_probes()`` (P37's six) — the names a window block must carry."""
    names = list(_probes())
    print(f"{TAG} p61._probes(): {names}")
    assert len(names) == 6 and len(set(names)) == 6


# ---------------------------------------------------------------------------
# The gated test.
# ---------------------------------------------------------------------------
def test_geo38_c39_reader():
    mode = _c39_mode()
    if mode is None:
        pytest.skip(f"{C39_ENV} unset")
    before = _on_root(
        lambda: sorted(p.name for p in HOST_JSON_DIR.glob("*")) if HOST_JSON_DIR.is_dir() else []
    )
    _c39_body(mode)
    if mode == "1":
        after = _on_root(
            lambda: (
                sorted(p.name for p in HOST_JSON_DIR.glob("*")) if HOST_JSON_DIR.is_dir() else []
            )
        )
        assert after == before, "wrote under logs/geo38 at knob 1"
        print(f"{TAG} knob 1: nothing written under logs/geo38")
    return None
