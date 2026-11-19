"""`WF-17` step 49 / U159 (§9 item 971): the clinical page's rod-in-box section is the logs.

``docs/guide/clinical_coil_fixtures.md`` section "Rod-in-box records (entries 48 /
61)" quotes `xl` entries 48 (1.5 T) and 61 (3 T) of `GEO-36` I36 (Q per tip ball,
eta_rod(T1), arm M's rod share) and `GEO-36` C40's 3 T / 1.5 T ratios.  Every quoted
digit is compared, as a string, with the token C45's ``parse_i36_log`` reads from the
tracked logs.  Knobless, no mesh, no solve.

Gates: (Chk) C45's checks on both tracked blocks; (Q-tab) table cells == parsed
tokens, node strings in the provenance sentence; (Q-ratio) ratio cells == C45's
``ratios`` at ``.6e`` == C40's record-log tuple (custody: first two tokens == the
parsed strings); (Q-all) no other decimal in the section but ``1.5``; (Disc)
disclaimer tokens, no "AED" / "HFSS"; (D) docrefs dead 0.  Controls: (c-digit),
(c-node), (c-sec).

C40's record log prints the C45 table twice: once from C45's self-copy control
(entry 48 against itself, ratios 1.000000e+00) and once from the real 1.5 T / 3 T
read.  (Q-ratio) therefore selects, per key, the tuples whose first two tokens are
the parsed entry 48 / entry 61 strings (custody), and requires exactly one distinct
such tuple.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from tests.package.test_geo38_clinical_fixture_page import PAGE  # noqa: E402
from tests.package.test_wf17_first_study_tutorial import (  # noqa: E402
    LINK_PASS_RE,
    _checker_stdout,
)
from tests.unit.test_doc_reference_exit_codes import run_checker  # noqa: E402
from tests.validation.test_geo36_c40_rod_frequency_reader import (  # noqa: E402
    ENTRY48_BALLS,
    ENTRY48_LOG,
    ENTRY48_NODE_PF,
    ENTRY61_GLOB,
    ENTRY61_NODE_PF,
    TRACKED_LOG_DIR,
    _on_root,
    _read_log,
    check_side,
    find_log,
    parse_i36_log,
    ratios,
    require_side,
)

TAG = "[WF-17 s49]"
HEADING = "## Rod-in-box records (entries 48 / 61)"
C40_GLOB = "*_GEO-36-C40-W1-n1.log"
ROW_KEYS = tuple((f"Q {b}", f"q.{b}") for b in ENTRY48_BALLS) + (
    ("η_rod(T1)", "eta_t1"),
    ("rod share (arm M)", "share.M"),
)
C40_LABEL = {
    **{f"q.{b}": f"tip 1 g Q {b}" for b in ENTRY48_BALLS},
    "eta_t1": "eta_rod(T1)",
    "share.M": "rod share arm M (rod)",
}
C40_RE = re.compile(r"^\[GEO-36 C45\]   (.+?): (\S+) \| (\S+) \| (\S+)$")
DISCLAIMERS = (
    "degree 1",
    "1.5 T gel constants",
    "Box containment is not asserted",
    "No published comparand",
    "A record",
)
DECIMAL_RE = re.compile(r"\d+\.\d+")


# ---------------------------------------------------------------------------
# Pure readers.
# ---------------------------------------------------------------------------
def section(text: str) -> str:
    """The heading line to the next ``## `` heading; ``LookupError`` unless the
    heading occurs exactly once."""
    heads = [m.start() for m in re.finditer(r"^" + re.escape(HEADING) + r"[ \t]*$", text, re.M)]
    if len(heads) != 1:
        raise LookupError(f"{HEADING!r}: {len(heads)} occurrences (exactly one required)")
    start = heads[0]
    nxt = re.search(r"^## ", text[start + len(HEADING) :], re.M)
    return text[start : start + len(HEADING) + nxt.start()] if nxt else text[start:]


def table(sec: str) -> dict[str, list[str]]:
    """``{row label: [1.5 T, 3 T, ratio]}`` from the section's table rows."""
    rows = {}
    for line in sec.splitlines():
        line = line.strip()
        if not line.startswith("|") or set(line) <= set("|-: "):
            continue
        cells = [c.strip() for c in line.strip("|").split("|")]
        if cells[0] == "Quantity":
            continue
        if cells[0] in rows:
            raise ValueError(f"row {cells[0]!r} repeated")
        rows[cells[0]] = cells[1:]
    return rows


def _val(block, key):
    if key.startswith("q."):
        return block["tip"]["balls"][key[2:]]["q"]
    if key == "eta_t1":
        return block["tip"]["eta"]
    return block["arms"]["M"]["share"]


def blocks() -> dict[str, dict]:
    """``{"1.5 T": block48, "3 T": block61}`` from the tracked logs."""
    hit = _on_root(find_log, ENTRY61_GLOB)
    assert hit is not None, f"entry 61: no {ENTRY61_GLOB}"
    path61, kind = hit
    print(f"{TAG} entry 48: {ENTRY48_LOG.name}; entry 61: {path61.name} ({kind})")
    assert kind == "tracked", kind
    return {"1.5 T": _on_root(parse_i36_log, ENTRY48_LOG), "3 T": _on_root(parse_i36_log, path61)}


def expected_rows(bl) -> dict[str, list[str]]:
    r = ratios(bl["1.5 T"], bl["3 T"])
    return {
        label: [_val(bl["1.5 T"], key), _val(bl["3 T"], key), format(r[key], ".6e")]
        for label, key in ROW_KEYS
    }


def tab_mismatches(sec: str, exp) -> list[tuple]:
    got = table(sec)
    bad = [(k, None, None) for k in got if k not in exp]
    for k, want in exp.items():
        have = got.get(k)
        if have is None:
            bad.append((k, None, want))
            continue
        for j, (h, w) in enumerate(zip(have, want)):
            if h != w:
                bad.append((k, j, h, w))
        if len(have) != len(want):
            bad.append((k, "len", len(have), len(want)))
    return bad


def c40_log() -> Path:
    hit = _on_root(find_log, C40_GLOB, ((TRACKED_LOG_DIR, "tracked"),))
    assert hit is not None, f"no tracked {C40_GLOB}"
    return hit[0]


def c40_tuples(path) -> dict[str, set[tuple[str, str, str]]]:
    by_label = {}
    for ln in _on_root(_read_log, path).splitlines():
        m = C40_RE.match(ln)
        if m:
            by_label.setdefault(m.group(1), set()).add(m.groups()[1:])
    return {key: by_label.get(lab, set()) for key, lab in C40_LABEL.items()}


def _bl_cache():
    if not hasattr(_bl_cache, "v"):
        _bl_cache.v = blocks()
    return _bl_cache.v


def _sec():
    return section(PAGE.read_text(encoding="utf-8"))


# ---------------------------------------------------------------------------
# Gates.
# ---------------------------------------------------------------------------
def test_s49_chk_c45_checks_on_tracked_logs():
    bl = _bl_cache()
    require_side(bl["1.5 T"], ENTRY48_NODE_PF, "1.5 T (entry 48)")
    require_side(bl["3 T"], ENTRY61_NODE_PF, "3 T (entry 61)")


def test_s49_q_tab_cells_are_the_parsed_tokens():
    bl = _bl_cache()
    sec = _sec()
    exp = expected_rows(bl)
    for k, v in exp.items():
        print(f"{TAG} (Q-tab) {k}: parsed {v}")
    bad = tab_mismatches(sec, exp)
    print(f"{TAG} (Q-tab) {len(exp)} rows x 3 cells compared; mismatches {bad}")
    assert bad == []
    for label in ("1.5 T", "3 T"):
        node = bl[label]["node"]["c_pf"]
        print(f"{TAG} (Q-tab) {label} node {node} pF on the page: {f'{node} pF' in sec}")
        assert f"{node} pF" in sec, node


def test_s49_q_ratio_c45_and_c40_custody():
    bl = _bl_cache()
    r = ratios(bl["1.5 T"], bl["3 T"])
    got = table(_sec())
    path = c40_log()
    tup = c40_tuples(path)
    for label, key in ROW_KEYS:
        t15, t3 = _val(bl["1.5 T"], key), _val(bl["3 T"], key)
        mine = format(r[key], ".6e")
        cust = {t for t in tup[key] if t[0] == t15 and t[1] == t3}
        print(
            f"{TAG} (Q-ratio) {key}: C40 tuples {sorted(tup[key])}; custody {sorted(cust)}; "
            f"C45 ratios {mine}; page {got[label][2]}"
        )
        assert len(cust) == 1, (key, tup[key])
        # (Q-self), 10:30 ruling: every other tuple is C45's self-copy control.
        others = set(tup[key]) - cust
        print(f"{TAG} (Q-self) {key}: others {sorted(others)}")
        assert others <= {(t15, t15, "1.000000e+00")}, (key, sorted(others))
        ((_a, _b, rr),) = cust
        assert rr == mine == got[label][2]
    assert path.name in _sec()


def test_s49_q_all_no_other_decimal():
    bl = _bl_cache()
    sec = _sec()
    exp = set()
    for v in expected_rows(bl).values():
        exp.update(v)
    exp.update(bl[k]["node"]["c_pf"] for k in ("1.5 T", "3 T"))
    rest = sec
    for s in sorted(exp, key=len, reverse=True):
        cnt = sec.count(s)
        print(f"{TAG} (Q-all) {s!r}: {cnt} occurrence(s)")
        assert cnt >= 1, s
        rest = rest.replace(s, "")
    stray = [t for t in DECIMAL_RE.findall(rest) if t != "1.5"]
    print(f"{TAG} (Q-all) {len(exp)} expected strings; stray decimals {stray}")
    assert stray == []


def test_s49_disc_disclaimers():
    page = PAGE.read_text(encoding="utf-8")
    sec = section(page)
    for d in DISCLAIMERS:
        print(f"{TAG} (Disc) {d!r}: {sec.count(d)}")
        assert d in sec, d
    assert "AED" not in page and "HFSS" not in page
    assert "~" not in sec


def test_s49_d_docrefs_guide():
    rc, counts = _on_root(run_checker, "--docs-root", "docs/guide", "--check-links")
    link = LINK_PASS_RE.search(
        _on_root(_checker_stdout, "--docs-root", "docs/guide", "--check-links")
    )
    assert link is not None, "(D) no 'Link pass:' line"
    print(
        f"{TAG} (D) dead={counts['dead']} links={link.group(1)} pages={link.group(2)} "
        f"link-dead={link.group(3)} rc={rc}"
    )
    assert counts["dead"] == 0
    assert int(link.group(3)) == 0


# ---------------------------------------------------------------------------
# Negative controls.
# ---------------------------------------------------------------------------
def test_s49_c_digit_one_mismatch():
    bl = _bl_cache()
    text = PAGE.read_text(encoding="utf-8")
    tok = bl["3 T"]["tip"]["balls"]["top"]["q"]
    new = tok[:-1] + str((int(tok[-1]) + 1) % 10)
    assert text.count(tok) == 1
    planted = text.replace(tok, new, 1)
    bad = tab_mismatches(section(planted), expected_rows(bl))
    print(f"{TAG} (c-digit) {tok} -> {new}: {bad}")
    assert bad == [("Q top", 1, new, tok)]


def test_s49_c_node_only_node_fails():
    bl = _bl_cache()
    res = check_side(bl["3 T"], ENTRY48_NODE_PF, "c-node")
    failed = [n for n, ok, _d in res if not ok]
    print(f"{TAG} (c-node) entry 61 against {ENTRY48_NODE_PF}: failed {failed} of {len(res)}")
    assert failed == [f"c-node (Node) node token == {ENTRY48_NODE_PF}"]


def test_s49_c_sec_missing_or_doubled_heading():
    text = PAGE.read_text(encoding="utf-8")
    assert text.count(HEADING) == 1
    with pytest.raises(LookupError):
        section(text.replace(HEADING, "## Something else"))
    with pytest.raises(LookupError):
        section(text.replace(HEADING, HEADING + "\n\n" + HEADING))
    print(f"{TAG} (c-sec) removed and doubled heading both raise LookupError")
