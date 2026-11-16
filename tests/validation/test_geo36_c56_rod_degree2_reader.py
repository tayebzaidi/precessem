"""`GEO-36` / chain I C56 — the 1.5 T clinical rod at degree 2 (§9 item 988;
docs/planning/chunks/GEO-38.md "§9 item 988").  No solve, no mesh.

C45's landed ``parse_i36_log`` reads entry 62's call-2 ``[GEO-36 I36a]`` block
(degree 2, node 185.7629 pF) beside entry 48's (degree 1, node 193.5996 pF).  Each
figure's degree-1 -> degree-2 move is printed, then the 10-24 exit's order rule's
verdict.  The verdict is printed, never asserted.  Every move mixes the order change
with the node change; C56 cannot separate the two and says so in its verdict line.

C45's module is imported, never edited: its ``check_side`` hard-codes ``(Deg) degree
1``, so ``check_degree2`` requires that to be the only red and asserts its own (Deg2).
No log token of >= 8 significant digits is a literal here: every number checked is
parsed at run time.
"""

from __future__ import annotations

import ast
import json
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

import tests.validation.test_geo36_c40_rod_frequency_reader as c45  # noqa: E402

ENTRY62_LOG = c45.TRACKED_LOG_DIR / "20261114T080009Z_GEO-38-C19.log"
ENTRY62_NODE_PF = "185.7629"
ENTRY62_JSON = c45.HOST_JSON_DIR / "i36-node-185.7629-d2.json"
ENTRY48_JSON = c45.HOST_JSON_DIR / f"i36-node-{c45.ENTRY48_NODE_PF}.json"
ENTRY62_FEEDS = ("P33", "P41")
#: The 10-24 exit's pre-registered threshold, as the ``ans:11`` page's
#: ``ORDER_RULE_REL``; not a band set here.
ORDER_RULE_REL = 0.02
ORDER_RULE_FIRES = ("clinical rod figures quoted at degree 2 only, degree-1 values as "
                    "records beside them")
ORDER_RULE_QUIET = "the 2 % order rule does not fire"
TAG = "[GEO-36 C56]"
#: The 11-14 weekly's quoted degree-2 figures (PROJECT_PLAN.md, entry-62 ruling:
#: "tip 1 g Q 4.145 / 4.133 / 4.116 / 4.089, eta_rod(T1) 3.443491, share 3.1332e-05").
WEEKLY_QUOTED = {
    "q.top": ("4.145", ".3f"),
    "q.bottom": ("4.133", ".3f"),
    "q.top_apex": ("4.116", ".3f"),
    "q.bottom_apex": ("4.089", ".3f"),
    "eta_t1": ("3.443491", None),
    "share.M": ("3.1332e-05", ".4e"),
}
N_RULE_ROWS = 10


def _parse(path, at_line=None):
    return c45._on_root(c45.parse_i36_log, path, at_line)


def _tmp_dir(tmp_path):
    return c45._on_root(lambda: str(tmp_path))


# ---------------------------------------------------------------------------
# Pure helpers.
# ---------------------------------------------------------------------------
def check_degree2(block, node_pf, name):
    """C45's ``check_side`` with its ``(Deg) degree 1`` red required as the **only**
    red (detail ``degree 2``), plus this module's (Deg2).  ``AssertionError`` naming
    every other failure; returns the full result list."""
    res = c45.check_side(block, node_pf, name)
    for n, ok, d in res:
        print(f"{TAG} {'PASS' if ok else 'FAIL'} {n}: {d}")
    deg_name = f"{name} (Deg) degree 1"
    failed = [(n, d) for n, ok, d in res if not ok]
    other = [(n, d) for n, d in failed if n != deg_name]
    deg = [(n, d) for n, d in failed if n == deg_name]
    if deg != [(deg_name, "degree 2")]:
        other.append((f"{name} (Deg) expected red with detail 'degree 2'", repr(deg)))
    ok2 = block["tip"]["degree"] == "2"
    res.append((f"{name} (Deg2) degree 2", ok2, f"degree {block['tip']['degree']}"))
    print(f"{TAG} {'PASS' if ok2 else 'FAIL'} {name} (Deg2) degree 2: "
          f"degree {block['tip']['degree']}")
    if not ok2:
        other.append((f"{name} (Deg2) degree 2", f"degree {block['tip']['degree']}"))
    if other:
        raise AssertionError(f"{name}: {other}")
    return res


def rod_rows(block):
    """``[(row, cls, value, shown)]``: the rule rows (``ratio`` / ``per-watt``) and the
    printed-only rows (``printed``, never in the rule).  ``shown`` is the log's own
    token for a parsed figure, ``.4e`` for a derived one."""
    tip, M, R = block["tip"], block["arms"]["M"], block["arms"]["R"]
    rows = []
    for ball, b in tip["balls"].items():
        rows.append((f"q.{ball}", "ratio", float(b["q"]), b["q"]))
    rows.append(("eta_t1", "ratio", float(tip["eta"]), tip["eta"]))
    rows.append(("share.M", "ratio", float(M["share"]), M["share"]))
    for ball, b in tip["balls"].items():
        v = float(b["sar_m"]) / float(M["p_vol"])
        rows.append((f"sar_m_per_p_vol.{ball}", "per-watt", v, f"{v:.4e}"))
    rows.append(("share.R", "printed", float(R["share"]), R["share"]))
    for ball, b in tip["balls"].items():
        v = float(b["sar_r"]) / float(R["p_vol"])
        rows.append((f"sar_r_per_p_vol.{ball}", "printed", v, f"{v:.4e}"))
    rows.append(("p_vol.M", "printed", float(M["p_vol"]), M["p_vol"]))
    rows.append(("rod.M", "printed", float(M["rod"]), M["rod"]))
    return rows


def moves_from_rows(rows1, rows2):
    """Per row ``x2 / x1 - 1`` at ``+.3e`` (C37's form); the rule on the non-printed
    rows.  Row names / classes must match, otherwise ``ValueError``."""
    if [(r[0], r[1]) for r in rows1] != [(r[0], r[1]) for r in rows2]:
        raise ValueError("row names / classes differ")
    table = []
    rows_over, rule_rows = [], []
    for r1, r2 in zip(rows1, rows2):
        row, cls = r1[0], r1[1]
        mv = r2[2] / r1[2] - 1
        over = cls != "printed" and abs(mv) > ORDER_RULE_REL
        if cls != "printed":
            rule_rows.append(row)
        if over:
            rows_over.append(row)
        table.append({"row": row, "cls": cls, "d1": r1[3], "d2": r2[3],
                      "move": f"{mv:+.3e}", "over": over})
    verdict = {"threshold": ORDER_RULE_REL, "rows_over": rows_over, "rule_rows": rule_rows,
               "sentence": ORDER_RULE_FIRES if rows_over else ORDER_RULE_QUIET}
    return table, verdict


def moves_and_rule(b1, b2):
    if list(b1["tip"]["balls"]) != list(b2["tip"]["balls"]):
        raise ValueError("ball names differ")
    return moves_from_rows(rod_rows(b1), rod_rows(b2))


def print_c56(b1, b2):
    table, verdict = moves_and_rule(b1, b2)
    print(f"{TAG} C56 table: degree 1 ({b1['log']}, node {b1['node']['c_pf']} pF) | "
          f"degree 2 ({b2['log']}, node {b2['node']['c_pf']} pF) | move | class | over")
    for t in table:
        over = "-" if t["cls"] == "printed" else ("over" if t["over"] else "under")
        print(f"{TAG}   {t['row']}: {t['d1']} | {t['d2']} | {t['move']} | {t['cls']} | {over}")
    for arm in ("M", "R"):
        f1 = b1["arms"][arm]["y_in"]["feeds"]
        f2 = b2["arms"][arm]["y_in"]["feeds"]
        for i in range(2):
            print(f"{TAG}   Im Y_in (1,+1j) arm {arm} feed {i + 1} ({f1[i]} | {f2[i]}) [S]: "
                  f"{b1['arms'][arm]['y_in'][f1[i]][1]} | {b2['arms'][arm]['y_in'][f2[i]][1]} "
                  f"(side by side, never a ratio)")
    print(f"{TAG} order rule (threshold {ORDER_RULE_REL}; the node moved "
          f"{b1['node']['c_pf']} -> {b2['node']['c_pf']} pF with the order): rows over "
          f"{verdict['rows_over']} -> \"{verdict['sentence']}\"")
    return table, verdict


def _value_equal(a, b):
    for conv in (float, complex):
        try:
            return conv(a) == conv(b)
        except (TypeError, ValueError):
            pass
    try:
        return ast.literal_eval(a) == ast.literal_eval(b)
    except (ValueError, SyntaxError):
        return False


def custody_values(block, path):
    """C45's ``custody_json`` with the 10:30 11-12 ruling: a mismatch equal as values
    is ``print_form_only`` (printed), any other is a value mismatch.  ``None`` when the
    host JSON is absent."""
    present, n, bad = c45.custody_json(block, path)
    if not present:
        print(f"{TAG} (J) {Path(path).name} absent; custody not run")
        return None
    form = [b for b in bad if _value_equal(b[1], b[2])]
    value = [b for b in bad if not _value_equal(b[1], b[2])]
    print(f"{TAG} (J) {Path(path).name}: n_fields {n}, value mismatches {value}, "
          f"print_form_only {form}")
    return {"n_fields": n, "value": value, "print_form_only": form}


# ---------------------------------------------------------------------------
# Tests.
# ---------------------------------------------------------------------------
def _read_entry62():
    n_node = c45._on_root(
        lambda: sum(1 for ln in c45._read_log(ENTRY62_LOG).splitlines()
                    if c45._RX_NODE.search(ln)))
    print(f"{TAG} {ENTRY62_LOG.name}: [GEO-36 I36a] (Z) node lines: {n_node}")
    b62 = _parse(ENTRY62_LOG)
    print(f"{TAG} (X62) parsed {b62['log']}: node {b62['node']['c_pf']} pF, block "
          f":{b62['node']['line']}-:{b62['wrote']['line']}; tip :{b62['tip']['line']}; balls "
          f"{list(b62['tip']['balls'])}, degree {b62['tip']['degree']}, eta valid "
          f"{b62['tip']['eta_valid']}; feeds M {b62['arms']['M']['y_in']['feeds']}, "
          f"R {b62['arms']['R']['y_in']['feeds']}")
    assert n_node == 1
    return b62


def _read_entry48():
    b48 = _parse(c45.ENTRY48_LOG)
    print(f"{TAG} (X48) parsed {b48['log']}: node {b48['node']['c_pf']} pF, block "
          f":{b48['node']['line']}-:{b48['wrote']['line']}; degree {b48['tip']['degree']}")
    return b48


def test_c56_x62_d62():
    """(X62) one block on entry 62's log; (D62) ``check_degree2`` green."""
    b62 = _read_entry62()
    assert b62["node"]["c_pf"] == ENTRY62_NODE_PF
    assert tuple(b62["tip"]["balls"]) == c45.ENTRY48_BALLS, list(b62["tip"]["balls"])
    for arm in ("M", "R"):
        assert b62["arms"][arm]["y_in"]["feeds"] == ENTRY62_FEEDS
    assert b62["tip"]["degree"] == "2"
    assert b62["tip"]["eta_valid"] == "True"
    check_degree2(b62, ENTRY62_NODE_PF, "1.5 T d2 (entry 62)")


def test_c56_x48_d48():
    """(X48) / (D48) as C45's, without ``_k_c32``."""
    b48 = _read_entry48()
    assert b48["node"]["c_pf"] == c45.ENTRY48_NODE_PF
    assert tuple(b48["tip"]["balls"]) == c45.ENTRY48_BALLS
    assert b48["tip"]["degree"] == "1" and b48["tip"]["eta_valid"] == "True"
    c45.require_side(b48, c45.ENTRY48_NODE_PF, "1.5 T d1 (entry 48)")


def test_c56_k_weekly():
    """(K-wk): the parsed degree-2 tokens at the weekly's digits."""
    b62 = _parse(ENTRY62_LOG)
    parsed = {f"q.{k}": v["q"] for k, v in b62["tip"]["balls"].items()}
    parsed["eta_t1"] = b62["tip"]["eta"]
    parsed["share.M"] = b62["arms"]["M"]["share"]
    got = {k: (parsed[k] if fmt is None else format(float(parsed[k]), fmt))
           for k, (_, fmt) in WEEKLY_QUOTED.items()}
    want = {k: s for k, (s, _) in WEEKLY_QUOTED.items()}
    print(f"{TAG} (K-wk) parsed {got} vs weekly {want}: {'PASS' if got == want else 'FAIL'}")
    assert got == want


def test_c56_moves_rule():
    """The record: the C56 table and the order rule's verdict (printed, not asserted)."""
    b48, b62 = _read_entry48(), _parse(ENTRY62_LOG)
    table, verdict = print_c56(b48, b62)
    rule = [t for t in table if t["cls"] != "printed"]
    assert len(rule) == N_RULE_ROWS
    assert [t["row"] for t in rule] == verdict["rule_rows"]
    assert {t["cls"] for t in rule} == {"ratio", "per-watt"}
    assert [t["row"] for t in table if t["cls"] == "ratio"] == (
        [f"q.{b}" for b in c45.ENTRY48_BALLS] + ["eta_t1", "share.M"])
    assert [t["row"] for t in table if t["cls"] == "per-watt"] == [
        f"sar_m_per_p_vol.{b}" for b in c45.ENTRY48_BALLS]
    assert [t["row"] for t in table if t["cls"] == "printed"] == (
        ["share.R"] + [f"sar_r_per_p_vol.{b}" for b in c45.ENTRY48_BALLS]
        + ["p_vol.M", "rod.M"])

    # Every move from the parsed strings, independently of rod_rows.
    def direct(b, row):
        t, M, R = b["tip"], b["arms"]["M"], b["arms"]["R"]
        kind, _, ball = row.partition(".")
        if kind == "q":
            return float(t["balls"][ball]["q"])
        if row == "eta_t1":
            return float(t["eta"])
        if kind == "share":
            return float(b["arms"][ball]["share"])
        if kind == "sar_m_per_p_vol":
            return float(t["balls"][ball]["sar_m"]) / float(M["p_vol"])
        if kind == "sar_r_per_p_vol":
            return float(t["balls"][ball]["sar_r"]) / float(R["p_vol"])
        if row == "p_vol.M":
            return float(M["p_vol"])
        if row == "rod.M":
            return float(M["rod"])
        raise KeyError(row)

    for t in table:
        assert t["move"] == f"{direct(b62, t['row']) / direct(b48, t['row']) - 1:+.3e}", t
    assert verdict["sentence"] == (ORDER_RULE_FIRES if verdict["rows_over"] else ORDER_RULE_QUIET)


def test_c56_j_custody():
    """(J62r) / (J48): value custody against the host JSONs where present."""
    for path, block in ((ENTRY62_JSON, _parse(ENTRY62_LOG)),
                        (ENTRY48_JSON, _parse(c45.ENTRY48_LOG))):
        out = custody_values(block, path)
        if out is None:
            continue
        assert out["value"] == [], out["value"]


def test_c56_n_deg():
    """(N-deg): entry 48's degree-1 block is refused, naming (Deg2)."""
    b48 = _parse(c45.ENTRY48_LOG)
    with pytest.raises(AssertionError, match=r"\(Deg2\)"):
        check_degree2(b48, c45.ENTRY48_NODE_PF, "N-deg (entry 48)")


def _copy62(tmp_dir, q_digit=False, dup_node=False):
    """Rank 0: entry 62's log copied to ``tmp_dir``; optionally top's Q changed in its
    4th decimal, or the (Z) node line duplicated in place."""
    b = c45.parse_i36_log(ENTRY62_LOG)
    lines = c45._read_log(ENTRY62_LOG).split("\n")
    i_node, i_tip = b["node"]["line"] - 1, b["tip"]["line"] - 1
    if q_digit:
        tok = b["tip"]["balls"]["top"]["q"]
        i = tok.index(".") + 4  # the 4th decimal: 4.145285 -> 4.145385
        new = tok[:i] + str((int(tok[i]) + 1) % 10) + tok[i + 1:]
        assert f"top {tok} (M " in lines[i_tip]
        lines[i_tip] = lines[i_tip].replace(f"top {tok} (M ", f"top {new} (M ", 1)
    if dup_node:
        lines.insert(i_node, lines[i_node])
    d = Path(tmp_dir)
    d.mkdir(parents=True, exist_ok=True)
    out = d / ENTRY62_LOG.name
    out.write_text("\n".join(lines))
    return str(out)


def test_c56_n_digit62(tmp_path):
    """(N-digit62): top's Q one digit off makes (D-Q) top red beside (Deg)."""
    d = _tmp_dir(tmp_path)
    path = c45._on_root(_copy62, Path(d) / "ndigit", True)
    b = _parse(path)
    name = "N-digit62"
    res = c45.check_side(b, ENTRY62_NODE_PF, name)
    failed = sorted(n for n, ok, _ in res if not ok)
    print(f"{TAG} (N-digit62) failed: {failed}")
    assert failed == sorted([f"{name} (D-Q) top", f"{name} (Deg) degree 1"]), failed
    with pytest.raises(AssertionError, match=r"\(D-Q\) top"):
        check_degree2(b, ENTRY62_NODE_PF, name)


def test_c56_n_blocks62(tmp_path):
    """(N-blocks62): the (Z) node line duplicated makes ``parse_i36_log`` raise."""
    d = _tmp_dir(tmp_path)
    path = c45._on_root(_copy62, Path(d) / "nblocks", False, True)
    with pytest.raises(ValueError, match="2 I36a blocks"):
        _parse(path)


def test_c56_n_rule():
    """(N-rule): x 1.019 QUIET, x 1.021 FIRES naming the row, a printed row x 1.5 QUIET."""
    base = [("q.top", "ratio", 1.0, "1"), ("share.R", "printed", 1.0, "1")]

    def scaled(fq, fp):
        return [("q.top", "ratio", fq, "-"), ("share.R", "printed", fp, "-")]

    _, v = moves_from_rows(base, scaled(1.019, 1.0))
    assert v["rows_over"] == [] and v["sentence"] == ORDER_RULE_QUIET, v
    _, v = moves_from_rows(base, scaled(1.021, 1.0))
    assert v["rows_over"] == ["q.top"] and v["sentence"] == ORDER_RULE_FIRES, v
    _, v = moves_from_rows(base, scaled(1.0, 1.5))
    assert v["rows_over"] == [] and v["sentence"] == ORDER_RULE_QUIET, v
    assert v["rule_rows"] == ["q.top"]


def test_c56_n_custody(tmp_path):
    """(N-custody): the d2 JSON with q.top scaled by 1 + 1e-5 gives one value mismatch."""
    present = c45._on_root(lambda: ENTRY62_JSON.exists())
    if not present:
        print(f"{TAG} (N-custody) {ENTRY62_JSON} absent; skipped (a clone)")
        return
    d = _tmp_dir(tmp_path)

    def _make():
        j = json.loads(ENTRY62_JSON.read_text())
        j["q"]["top"] = j["q"]["top"] * (1 + 1e-5)
        p = Path(d) / "ncustody" / ENTRY62_JSON.name
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(json.dumps(j))
        return str(p)

    path = c45._on_root(_make)
    out = custody_values(_parse(ENTRY62_LOG), path)
    assert out is not None
    assert [b[0] for b in out["value"]] == ["q.top"], out["value"]
