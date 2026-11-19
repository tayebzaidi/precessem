"""`GEO-38` / chain C C45 — C40's rod 1.5 T / 3 T reader (§9 item 906;
docs/planning/chunks/GEO-38.md "§9 item 906").  No solve, no mesh.

C40 reads entries 48 and 61: tip 1 g Q on the four balls, eta_rod(T1), the rod's
power share, Im Y_in at each root — the frequency lever on one rod model and one
definition.  The host JSONs C40 names (``logs/geo36/i36-node-<pF>.json``) live under
the gitignored ``/logs/`` and are not reproducible from a clone, so this reader reads
the **tracked** ``[GEO-36 I36a]`` lines of the I36a logs (item 905's route; the
11-05 03:00 ruling "the tracked log is the source"); the host JSON is a custody
cross-check only.  No existing module parses these lines, so the regexes are new.

Anchors: (X48) one block on entry 48's log, node 193.5996, the four balls, degree 1,
eta valid; (D48) Q = SAR_M / SAR_R per ball and share = rod / P_vol per arm at the
half-unit bounds; (K-C32) C32's print-only ``ENTRY48_Q_RANGE`` / ``ENTRY48_ETA_T1``
equal the parsed min / max Q at ``.6f`` and eta at ``.3f``; (J48) the host JSON's
fields format to their log strings; (Stub61) C32's tracked stub block (3 T body under
C21's patches) parses; (Self) a renamed 3 T copy of entry 48's log reads with ratios
1; (Pend) entry 61 absent.  Controls: (N-digit), (N-drop), (N-node), (N-blocks),
(N-balls).

Knob ``FEM_EM_GEO36_C40``: unset skip; ``1`` entry 48 required, entry 61 optional
(nothing written); ``both`` both required, rank 0 writes ``c40-rod-1p5t-3t.json``
after every check holds; anything else ``ValueError``.  No log token of >= 8
significant digits is a literal here: every number checked is parsed at run time.
"""

from __future__ import annotations

import gzip
import json
import os
import re
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

C40_ENV = "FEM_EM_GEO36_C40"
#: Entry 61's XL_CHUNK (xl-pending.md entry 61).
ENTRY61_GLOB = "*_GEO-36-step11-I36-3T.log"
#: Entry 61's command node, at four decimals.
ENTRY61_NODE_PF = "44.0525"
ENTRY48_NODE_PF = "193.5996"
TRACKED_LOG_DIR = REPO / "docs" / "testing" / "logs"
HELD_LOG_DIR = REPO / "logs" / "leak-held"
#: Entry 48's tracked window log (I36a at the 1.5 T node).
ENTRY48_LOG = TRACKED_LOG_DIR / "20261102T080008Z_GEO-36-step9-I36.log"
#: C32's tracked stub log (I36a's lines on the 3 T body under C21's patches).
C32_STUB_LOG = TRACKED_LOG_DIR / "20261105T133319Z_GEO-36-C32-W0.log"
C32_STUB_LINE = 81
ENTRY48_BALLS = ("top", "bottom", "top_apex", "bottom_apex")
HOST_JSON_DIR = Path("/workspace/logs/geo36")
C40_JSON_NAME = "c40-rod-1p5t-3t.json"
TAG = "[GEO-36 C45]"


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
    ``path`` + ``.gz`` (the ``ans:11`` page's ``_read_log``)."""
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


_RX_NODE = re.compile(r"\[GEO-36 I36a\].*?\(Z\) node ([0-9.]+) pF")
_RX_WROTE = re.compile(r"\[GEO-36 I36a\] wrote (\S+)")
_C = r"([^,;\s]+)"
_RX_ARM = re.compile(
    r"\[GEO-36 I36a\] arm (M|R) \(([^)]*)\): kinds .*?; S2 "
    + _C
    + ", "
    + _C
    + "; "
    + _C
    + ", "
    + _C
    + ";"
)
_RX_YIN = re.compile(
    r"\[GEO-36 I36a\]\s+Y_in \(1,\+1j\): (P\d+) (\S+) (\S+)j S, "
    r"(P\d+) (\S+) (\S+)j S .*?P_vol (\S+) W, gel ([^,\s]+), "
    r"rod (\S+) \(rod share (\S+), printed\)"
)
_RX_TIP = re.compile(r"\[GEO-36 I36a\] tip 1 g ratio Q = .*?degree (\d):")
_RX_BALL = re.compile(r"(\w+) ([0-9.]+) \(M (\S+), R (\S+) W/kg, mass (\S+) kg\)")
_RX_ETA = re.compile(r"eta_rod\(T1 = \((.*?)\)\) = \|E_M\|/\|E_R\| (\S+) \(valid (True|False)")


def parse_i36_log(path, at_line=None) -> dict:
    """Pure: one I36a block (a ``(Z) node`` line to the next ``wrote`` line), every
    token kept as its printed **string** with its 1-based line number.  A window log
    must hold exactly one block; a stub log selects its block by ``at_line``, the line
    number of its ``(Z) node`` line.  ``ValueError`` on any missing / repeated line."""
    lines = _read_log(path).splitlines()
    name = Path(str(path)).name.removesuffix(".gz")
    nodes = [i for i, ln in enumerate(lines, start=1) if _RX_NODE.search(ln)]
    if at_line is None:
        if len(nodes) != 1:
            raise ValueError(
                f"{name}: {len(nodes)} I36a blocks (lines {nodes[:6]}); "
                f"a window log holds exactly one"
            )
        start = nodes[0]
    else:
        if at_line not in nodes:
            raise ValueError(f"{name}: no (Z) node line at :{at_line}")
        start = at_line
    end = None
    for i in range(start, len(lines) + 1):
        if _RX_WROTE.search(lines[i - 1]):
            end = i
            break
    if end is None:
        raise ValueError(f"{name}: block at :{start} has no wrote line")
    nxt = [n for n in nodes if start < n < end]
    if nxt:
        raise ValueError(f"{name}: block at :{start} holds another node line at :{nxt[0]}")
    out = {
        "log": name,
        "node": {"c_pf": _RX_NODE.search(lines[start - 1]).group(1), "line": start},
        "wrote": {"path": _RX_WROTE.search(lines[end - 1]).group(1), "line": end},
        "arms": {},
    }
    cur = None
    tip = None
    for i in range(start, end + 1):
        ln = lines[i - 1]
        m = _RX_ARM.search(ln)
        if m:
            arm = m.group(1)
            if arm in out["arms"]:
                raise ValueError(f"{name}: arm {arm} repeated at :{i}")
            out["arms"][arm] = {
                "label": m.group(2),
                "line": i,
                "s2": {"s11": m.group(3), "s12": m.group(4), "s21": m.group(5), "s22": m.group(6)},
            }
            cur = arm
            continue
        m = _RX_YIN.search(ln)
        if m:
            if cur is None or "y_in" in out["arms"][cur]:
                raise ValueError(f"{name}: Y_in line at :{i} without its arm header")
            out["arms"][cur]["y_in"] = {
                "line": i,
                "feeds": (m.group(1), m.group(4)),
                m.group(1): (m.group(2), m.group(3)),
                m.group(4): (m.group(5), m.group(6)),
            }
            out["arms"][cur]["p_vol"] = m.group(7)
            out["arms"][cur]["gel"] = m.group(8)
            out["arms"][cur]["rod"] = m.group(9)
            out["arms"][cur]["share"] = m.group(10)
            continue
        m = _RX_TIP.search(ln)
        if m:
            if tip is not None:
                raise ValueError(f"{name}: tip line repeated at :{i}")
            seg_end = ln.find("eta_rod(")
            seg = ln[m.end() : seg_end if seg_end >= 0 else len(ln)]
            balls = {}
            for b in _RX_BALL.finditer(seg):
                if b.group(1) in balls:
                    raise ValueError(f"{name}: ball {b.group(1)} repeated at :{i}")
                balls[b.group(1)] = {
                    "q": b.group(2),
                    "sar_m": b.group(3),
                    "sar_r": b.group(4),
                    "mass": b.group(5),
                }
            e = _RX_ETA.search(ln)
            if not balls or e is None:
                raise ValueError(f"{name}: tip line at :{i} without balls or eta")
            tip = {
                "line": i,
                "degree": m.group(1),
                "balls": balls,
                "t1": e.group(1),
                "eta": e.group(2),
                "eta_valid": e.group(3),
            }
    for arm in ("M", "R"):
        if arm not in out["arms"] or "y_in" not in out["arms"][arm]:
            raise ValueError(f"{name}: block at :{start} lacks arm {arm} or its Y_in line")
    if tip is None:
        raise ValueError(f"{name}: block at :{start} lacks the tip line")
    out["tip"] = tip
    return out


def _half_ulp_rel(tok):
    """Half a unit in the printed token's last digit, relative to its value."""
    t = tok.strip()
    mant, _, exp = t.lower().partition("e")
    dec = mant.split(".")[1] if "." in mant else ""
    abs_ulp = 0.5 * 10.0 ** ((int(exp) if exp else 0) - len(dec))
    return abs_ulp / abs(float(t))


def check_side(block, node_pf, name):
    """Spec-free derived checks on one block: ``[(check, ok, detail)]``, every check
    named with ``name``."""
    res = []

    def rel_check(tag, val, tok, inputs):
        bound = _half_ulp_rel(tok) + sum(_half_ulp_rel(t) for t in inputs)
        err = abs(val - float(tok)) / abs(float(tok))
        res.append(
            (
                f"{name} {tag}",
                err <= bound,
                f"{val:.6e} vs printed {tok}: rel {err:.3e} (<= {bound:.3e})",
            )
        )

    for ball, b in block["tip"]["balls"].items():
        rel_check(
            f"(D-Q) {ball}", float(b["sar_m"]) / float(b["sar_r"]), b["q"], (b["sar_m"], b["sar_r"])
        )
    for arm in ("M", "R"):
        a = block["arms"][arm]
        rel_check(
            f"(D-share) arm {arm}",
            float(a["rod"]) / float(a["p_vol"]),
            a["share"],
            (a["rod"], a["p_vol"]),
        )
    res.append(
        (
            f"{name} (Node) node token == {node_pf}",
            block["node"]["c_pf"] == node_pf,
            f"log node {block['node']['c_pf']} vs {node_pf}",
        )
    )
    res.append(
        (
            f"{name} (Deg) degree 1",
            block["tip"]["degree"] == "1",
            f"degree {block['tip']['degree']}",
        )
    )
    res.append(
        (
            f"{name} (eta-valid) eta valid True",
            block["tip"]["eta_valid"] == "True",
            f"valid {block['tip']['eta_valid']}",
        )
    )
    return res


def require_side(block, node_pf, name):
    res = check_side(block, node_pf, name)
    for n, ok, d in res:
        print(f"{TAG} {'PASS' if ok else 'FAIL'} {n}: {d}")
    bad = [(n, d) for n, ok, d in res if not ok]
    if bad:
        raise AssertionError(f"{name}: {bad}")
    return res


def ratios(s15, s3):
    """3 T / 1.5 T ratios as floats (print at ``.6e``): Q per ball, eta, arm M's share.
    The ball names must match, otherwise ``ValueError``."""
    b15, b3 = s15["tip"]["balls"], s3["tip"]["balls"]
    if list(b15) != list(b3):
        raise ValueError(f"ball names differ: {list(b15)} vs {list(b3)}")
    out = {f"q.{k}": float(b3[k]["q"]) / float(b15[k]["q"]) for k in b15}
    out["eta_t1"] = float(s3["tip"]["eta"]) / float(s15["tip"]["eta"])
    out["share.M"] = float(s3["arms"]["M"]["share"]) / float(s15["arms"]["M"]["share"])
    return out


def print_table(s15, s3=None):
    pend = "pending (xl entry 61)"
    r = ratios(s15, s3) if s3 is not None else None
    print(
        f"{TAG} C40 table: 1.5 T ({s15['log']}, node {s15['node']['c_pf']} pF) | "
        f"3 T ({s3['log'] + ', node ' + s3['node']['c_pf'] + ' pF' if s3 else pend}) "
        f"| ratio 3T/1.5T"
    )

    def rr(key):
        return format(r[key], ".6e") if r else pend

    for ball, b in s15["tip"]["balls"].items():
        q3 = s3["tip"]["balls"][ball]["q"] if s3 else pend
        print(f"{TAG}   tip 1 g Q {ball}: {b['q']} | {q3} | {rr('q.' + ball)}")
    print(
        f"{TAG}   eta_rod(T1): {s15['tip']['eta']} | {s3['tip']['eta'] if s3 else pend} | "
        f"{rr('eta_t1')}"
    )
    print(
        f"{TAG}   rod share arm M (rod): {s15['arms']['M']['share']} | "
        f"{s3['arms']['M']['share'] if s3 else pend} | {rr('share.M')}"
    )
    print(
        f"{TAG}   rod share arm R (gel twin, numerical floor): {s15['arms']['R']['share']} | "
        f"{s3['arms']['R']['share'] if s3 else pend} (printed beside, never a ratio)"
    )
    for arm in ("M", "R"):
        f15 = s15["arms"][arm]["y_in"]["feeds"]
        f3 = s3["arms"][arm]["y_in"]["feeds"] if s3 else (None, None)
        for i in range(2):
            v15 = s15["arms"][arm]["y_in"][f15[i]][1]
            v3 = s3["arms"][arm]["y_in"][f3[i]][1] if s3 else pend
            print(
                f"{TAG}   Im Y_in (1,+1j) arm {arm} feed {i + 1} ({f15[i]} | "
                f"{f3[i] if s3 else '-'}) [S]: {v15} | {v3} (side by side, never a ratio)"
            )
    return r


def custody_json(block, path):
    """Rank 0 reads the host JSON when present; each C40 field formatted as its log
    token is printed must equal that token's string.  ``(present, n_fields, mismatches)``."""

    def _read():
        p = Path(path)
        return json.loads(p.read_text()) if p.exists() else None

    j = _on_root(_read)
    if j is None:
        print(f"{TAG} (J) host JSON absent: {path}")
        return (False, 0, [])
    f6, f6e, f9, fs = "{:.6f}".format, "{:.6e}".format, "{:.9e}".format, "{:+.9e}".format
    pairs = []
    tip = block["tip"]
    if list(j["q"]) != list(tip["balls"]):
        pairs.append(("q.names", repr(list(j["q"])), repr(list(tip["balls"]))))
    for ball, b in tip["balls"].items():
        pairs.append((f"q.{ball}", f6(j["q"][ball]) if ball in j["q"] else "absent", b["q"]))
    pairs.append(("eta_t1", f6(j["eta_t1"]), tip["eta"]))
    for arm in ("M", "R"):
        a, ja = block["arms"][arm], j[arm]
        for ip, feed in enumerate(a["y_in"]["feeds"]):
            for ipart in range(2):
                pairs.append(
                    (
                        f"{arm}.y_in.{feed}.{('re', 'im')[ipart]}",
                        fs(ja["y_in"][ip][ipart]),
                        a["y_in"][feed][ipart],
                    )
                )
        ps = ja["p_sup"]
        pairs.append((f"{arm}.share (rod / total)", f6e(ps["rod"] / ps["total"]), a["share"]))
        pairs.append((f"{arm}.p_sup.total", f9(ps["total"]), a["p_vol"]))
        pairs.append((f"{arm}.p_sup.rod", f9(ps["rod"]), a["rod"]))
        for ball, b in tip["balls"].items():
            key = "sar_m" if arm == "M" else "sar_r"
            pairs.append((f"{arm}.sar.{ball}", f6e(ja["sar"][ball]), b[key]))
    pairs.append(("feeds", repr(list(j["feeds"])), repr(list(block["arms"]["M"]["y_in"]["feeds"]))))
    pairs.append(("degree", str(j["degree"]), tip["degree"]))
    bad = [(n, a, b) for n, a, b in pairs if a != b]
    print(
        f"{TAG} (J) host JSON {Path(path).name}: {len(pairs)} fields, "
        f"{len(bad)} mismatches{'; first ' + repr(bad[0]) if bad else ''}"
    )
    return (True, len(pairs), bad)


def _c40_mode():
    raw = os.environ.get(C40_ENV, "").strip()
    if raw == "":
        return None
    if raw in ("1", "both"):
        return raw
    raise ValueError(f"{C40_ENV} must be unset, '1' or 'both', got {raw!r}")


def _parse(path, at_line=None):
    return _on_root(parse_i36_log, path, at_line)


def _k_c32(b48):
    """(K-C32): C32's print-only constants against the parsed entry 48 block."""
    from tests.validation.test_geo36_rod_f_clinical_3t_drive import ENTRY48_ETA_T1, ENTRY48_Q_RANGE

    qs = [float(b["q"]) for b in b48["tip"]["balls"].values()]
    got = (f"{min(qs):.6f}", f"{max(qs):.6f}", f"{float(b48['tip']['eta']):.3f}")
    want = (f"{ENTRY48_Q_RANGE[0]:.6f}", f"{ENTRY48_Q_RANGE[1]:.6f}", f"{ENTRY48_ETA_T1:.3f}")
    print(
        f"{TAG} (K-C32) parsed min / max Q, eta {got} vs C32's ENTRY48_Q_RANGE / "
        f"ENTRY48_ETA_T1 {want}: {'PASS' if got == want else 'FAIL'}"
    )
    return got, want


def _read_entry48():
    b48 = _parse(ENTRY48_LOG)
    print(
        f"{TAG} (X48) parsed {b48['log']}: node {b48['node']['c_pf']} pF (line "
        f":{b48['node']['line']}), block to :{b48['wrote']['line']}; balls "
        f"{list(b48['tip']['balls'])}, degree {b48['tip']['degree']}, eta valid "
        f"{b48['tip']['eta_valid']}; feeds {b48['arms']['M']['y_in']['feeds']}"
    )
    assert b48["node"]["c_pf"] == ENTRY48_NODE_PF, b48["node"]
    assert tuple(b48["tip"]["balls"]) == ENTRY48_BALLS, list(b48["tip"]["balls"])
    assert b48["tip"]["degree"] == "1" and b48["tip"]["eta_valid"] == "True"
    require_side(b48, ENTRY48_NODE_PF, "1.5 T (entry 48)")
    got, want = _k_c32(b48)
    assert got == want, f"(K-C32) {got} vs {want}"
    return b48


def _c40_body(mode):
    """The gated body; returns the summary dict.  Writes only at ``both``."""
    b48 = _read_entry48()
    present, n, bad = custody_json(b48, HOST_JSON_DIR / f"i36-node-{ENTRY48_NODE_PF}.json")
    assert present, "(J48) host JSON absent on this host"
    assert not bad, f"(J48) {len(bad)} of {n} mismatches, first {bad[0]}"
    hit = _on_root(find_log, ENTRY61_GLOB)
    b61 = None
    if hit is None:
        print(f"{TAG} entry 61: pending (no {ENTRY61_GLOB} tracked or held)")
        if mode == "both":
            raise AssertionError(f"entry 61 (xl, {ENTRY61_GLOB}) absent: C40 needs both")
    else:
        path61, kind = hit
        print(f"{TAG} entry 61: {path61.name} ({kind})")
        b61 = _parse(path61)
        require_side(b61, ENTRY61_NODE_PF, "3 T (entry 61)")
    r = print_table(b48, b61)
    out = {
        "entry48": b48,
        "entry61": b61,
        "ratios": {k: f"{v:.6e}" for k, v in r.items()} if r else None,
    }
    if mode == "both":

        def _write():
            HOST_JSON_DIR.mkdir(parents=True, exist_ok=True)
            p = HOST_JSON_DIR / C40_JSON_NAME
            p.write_text(json.dumps(out, indent=1))
            return str(p)

        print(f"{TAG} wrote {_on_root(_write)}")
    return out


# ---------------------------------------------------------------------------
# Coding tests (no knob).
# ---------------------------------------------------------------------------
def _renamed_copy(tmp_dir, node_to=ENTRY61_NODE_PF, q_digit=False, drop_tip=False):
    """Rank 0: entry 48's log copied as a 3 T name with only the node token rewritten;
    optionally ``top``'s Q changed in its 5th decimal, or the tip line deleted."""
    b48 = parse_i36_log(ENTRY48_LOG)
    lines = _read_log(ENTRY48_LOG).split("\n")
    i_node, i_tip = b48["node"]["line"] - 1, b48["tip"]["line"] - 1
    old = b48["node"]["c_pf"]
    assert f"(Z) node {old} pF" in lines[i_node]
    lines[i_node] = lines[i_node].replace(f"(Z) node {old} pF", f"(Z) node {node_to} pF", 1)
    if q_digit:
        tok = b48["tip"]["balls"]["top"]["q"]
        i = tok.index(".") + 4  # the 4th decimal: 4.460568 -> 4.460668
        new = tok[:i] + str((int(tok[i]) + 1) % 10) + tok[i + 1 :]
        assert f"top {tok} (M " in lines[i_tip]
        lines[i_tip] = lines[i_tip].replace(f"top {tok} (M ", f"top {new} (M ", 1)
    if drop_tip:
        del lines[i_tip]
    d = Path(tmp_dir)
    d.mkdir(parents=True, exist_ok=True)
    out = d / "20990101T000000Z_GEO-36-step11-I36-3T.log"
    out.write_text("\n".join(lines))
    return str(out)


def _tmp_dir(tmp_path):
    return _on_root(lambda: str(tmp_path))


def test_c45_x48_d48_k_c32():
    """(X48), (D48), (K-C32) on entry 48's tracked log (no knob)."""
    b48 = _read_entry48()
    assert len(b48["tip"]["balls"]) == 4


def test_c45_stub61():
    """(Stub61): C32's stub block at :81 (the 3 T body under C21's patches) parses."""
    b = _parse(C32_STUB_LOG, C32_STUB_LINE)
    print(
        f"{TAG} (Stub61) {b['log']}:{b['node']['line']}: node {b['node']['c_pf']}, balls "
        f"{ {k: v['q'] for k, v in b['tip']['balls'].items()} }, eta {b['tip']['eta']}, "
        f"share M {b['arms']['M']['share']}, feeds {b['arms']['M']['y_in']['feeds']}"
    )
    assert list(b["tip"]["balls"]) == ["a", "b"]
    for v in b["tip"]["balls"].values():
        assert v["q"] == "4.000000"
        assert float(v["sar_m"]) / float(v["sar_r"]) == float(v["q"])
    assert b["tip"]["balls"]["a"]["sar_m"] == "4.500000e-03"
    assert b["tip"]["balls"]["a"]["sar_r"] == "1.125000e-03"
    assert b["tip"]["eta"] == "2.000000"
    assert b["arms"]["M"]["share"] == "1.041667e-02"
    assert f"{float(b['arms']['M']['rod']) / float(b['arms']['M']['p_vol']):.6e}" == "1.041667e-02"
    assert b["arms"]["M"]["y_in"]["feeds"] == ("P1", "P2")
    assert b["arms"]["R"]["y_in"]["feeds"] == ("P1", "P2")
    res = check_side(b, b["node"]["c_pf"], "3 T stub (C32 :81)")
    for n, ok, d in res:
        print(f"{TAG} {'PASS' if ok else 'FAIL'} {n}: {d}")
    dq = [ok for n, ok, _ in res if "(D-Q)" in n or "(D-share)" in n]
    assert len(dq) == 4 and all(dq), res
    assert all(ok for _, ok, _ in res), res


def test_c45_self_renamed_copy(tmp_path):
    """(Self): the 3 T path on a renamed copy; every ratio is 1."""
    b48 = _parse(ENTRY48_LOG)
    d = _tmp_dir(tmp_path)
    path = _on_root(_renamed_copy, Path(d) / "self")
    hit = _on_root(find_log, ENTRY61_GLOB, ((Path(d) / "self", "tracked"),))
    assert hit is not None and hit[1] == "tracked" and hit[0].name == Path(path).name, hit
    b3 = _parse(hit[0])
    assert b3["node"]["c_pf"] == ENTRY61_NODE_PF
    require_side(b3, ENTRY61_NODE_PF, "3 T (self copy)")
    r = print_table(b48, b3)
    assert all(f"{v:.6e}" == "1.000000e+00" for v in r.values()), r
    assert b3["arms"] == b48["arms"] and b3["tip"] == b48["tip"]


def test_c45_n_digit(tmp_path):
    """(N-digit): ``top``'s Q changed in one digit makes (D-Q) red on ``top``."""
    d = _tmp_dir(tmp_path)
    path = _on_root(_renamed_copy, Path(d) / "ndigit", ENTRY61_NODE_PF, True)
    b3 = _parse(path)
    res = check_side(b3, ENTRY61_NODE_PF, "3 T (N-digit)")
    failed = {n: dd for n, ok, dd in res if not ok}
    print(f"{TAG} (N-digit) failed: {failed}")
    assert list(failed) == ["3 T (N-digit) (D-Q) top"], failed
    with pytest.raises(AssertionError, match=r"\(D-Q\) top"):
        require_side(b3, ENTRY61_NODE_PF, "3 T (N-digit)")


def test_c45_n_drop(tmp_path):
    """(N-drop): the tip line deleted raises ``ValueError``."""
    d = _tmp_dir(tmp_path)
    path = _on_root(_renamed_copy, Path(d) / "ndrop", ENTRY61_NODE_PF, False, True)
    with pytest.raises(ValueError, match="tip line"):
        _parse(path)


def test_c45_n_node(tmp_path):
    """(N-node): the 3 T name with the node left at entry 48's makes (Node) raise."""
    d = _tmp_dir(tmp_path)
    old = _parse(ENTRY48_LOG)["node"]["c_pf"]
    path = _on_root(_renamed_copy, Path(d) / "nnode", old)
    b3 = _parse(path)
    with pytest.raises(AssertionError, match=r"\(Node\)"):
        require_side(b3, ENTRY61_NODE_PF, "3 T (N-node)")


def test_c45_n_blocks():
    """(N-blocks): C32's stub log without ``at_line`` (several blocks) raises."""
    with pytest.raises(ValueError, match="I36a blocks"):
        _parse(C32_STUB_LOG)


def test_c45_n_balls():
    """(N-balls): entry 48's block paired with the stub block (a / b) raises."""
    b48 = _parse(ENTRY48_LOG)
    stub = _parse(C32_STUB_LOG, C32_STUB_LINE)
    with pytest.raises(ValueError, match="ball names differ"):
        ratios(b48, stub)


def test_c45_n_two(tmp_path):
    """``find_log``: two files matching the glob raise ``ValueError``."""
    d = Path(_tmp_dir(tmp_path)) / "ntwo"

    def _make():
        d.mkdir(parents=True, exist_ok=True)
        (d / "20990101T000000Z_GEO-36-step11-I36-3T.log").write_text("a\n")
        (d / "20990102T000000Z_GEO-36-step11-I36-3T.log").write_text("b\n")

    _on_root(_make)
    with pytest.raises(ValueError, match="2 matches"):
        _on_root(find_log, ENTRY61_GLOB, ((d, "tracked"),))


def _empty_log_roots(monkeypatch, tmp_path):
    """Point ``find_log``'s default roots at two empty tmp dirs (rank 0 makes them,
    bcast), so a premise test reads entry 61 absent by construction (§9 item 918)."""

    def _make():
        t, h = Path(tmp_path) / "tracked-empty", Path(tmp_path) / "held-empty"
        t.mkdir(parents=True, exist_ok=True)
        h.mkdir(parents=True, exist_ok=True)
        return str(t), str(h)

    t, h = _on_root(_make)
    monkeypatch.setitem(globals(), "TRACKED_LOG_DIR", Path(t))
    monkeypatch.setitem(globals(), "HELD_LOG_DIR", Path(h))


def _json_stamp(p):
    """Rank 0's ``(mtime_ns, bytes)`` of ``p``, or ``None`` if absent: "nothing written"
    holds whether or not C40's own window has already written the file (§9 item 917)."""
    return _on_root(lambda: (p.stat().st_mtime_ns, p.read_bytes()) if p.exists() else None)


def test_c45_n_both(monkeypatch, tmp_path):
    """At ``both`` on empty log roots the gated body fails, naming entry 61."""
    _empty_log_roots(monkeypatch, tmp_path)
    before = _json_stamp(HOST_JSON_DIR / C40_JSON_NAME)
    monkeypatch.setenv(C40_ENV, "both")
    mode = _c40_mode()
    assert mode == "both"
    with pytest.raises(AssertionError, match="entry 61"):
        _c40_body(mode)
    assert _json_stamp(HOST_JSON_DIR / C40_JSON_NAME) == before


def test_c45_pend_and_knob(monkeypatch, tmp_path):
    """(Pend) entry 61 absent on empty log roots; the knob's shape."""
    _empty_log_roots(monkeypatch, tmp_path)
    assert _on_root(find_log, ENTRY61_GLOB) is None
    monkeypatch.setenv(C40_ENV, "bogus")
    with pytest.raises(ValueError):
        _c40_mode()


# ---------------------------------------------------------------------------
# The gated test.
# ---------------------------------------------------------------------------
def test_geo36_c40_reader():
    mode = _c40_mode()
    if mode is None:
        pytest.skip(f"{C40_ENV} unset")
    before = _on_root(
        lambda: sorted(p.name for p in HOST_JSON_DIR.glob("*")) if HOST_JSON_DIR.is_dir() else []
    )
    out = _c40_body(mode)
    if mode == "1":
        after = _on_root(
            lambda: (
                sorted(p.name for p in HOST_JSON_DIR.glob("*")) if HOST_JSON_DIR.is_dir() else []
            )
        )
        assert after == before, "wrote under logs/geo36 at knob 1"
        assert out["entry61"] is None or out["ratios"] is not None
        print(f"{TAG} knob 1: nothing written under logs/geo36")
