"""`WF-17` step 46 / U156 (§9 item 960): the clinical page's heating section is the logs.

``docs/guide/clinical_coil_fixtures.md`` section "Clinical heating records (entries
59 / 60)" quotes `xl` entries 59 (1.5 T, P33) and 60 (3 T, P33-3T) per accepted
watt, and `GEO-38` C39's 3 T / 1.5 T hotspot ratio.  Every quoted digit is compared,
as a string, with the token C46's parser reads from the tracked logs (no number of
>= 8 significant digits is a literal here).  Knobless, no mesh, no solve.

Gates: (Chk) C46's checks on both tracked blocks; (Q-tab) table cells == parsed
tokens; (Q-ratio) ratio cells == C46's ``ratios`` at ``.6e`` == C39's log tuple
(custody: its first two tokens == the parsed ``theta_point``); (Q-all) no other
decimal in the section but ``1.5``; (Disc) disclaimer tokens, no "AED" / "HFSS";
(D) docrefs dead 0.  Controls: (c-digit), (c-node), (c-sec).
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
from tests.validation.test_geo38_c39_thermal_reader import (  # noqa: E402
    ENTRY59_GLOB,
    ENTRY59_NODE_PF,
    ENTRY60_GLOB,
    ENTRY60_NODE_PF,
    SENSES,
    TRACKED_LOG_DIR,
    _on_root,
    _probes,
    _read_log,
    check_block,
    find_log,
    parse_heat_log,
    ratios,
    require_block,
    side,
)

TAG = "[WF-17 s46]"
HEADING = "## Clinical heating records (entries 59 / 60)"
C39_LOG = TRACKED_LOG_DIR / "20261111T103024Z_GEO-38-C39-W1.log"
RATIO_KEY = "hotspot point theta per P_acc"
C39_RE = re.compile(r"^\[GEO-38 C46\]   (ccw|cw) hotspot point theta per P_acc: "
                    r"(\S+) \| (\S+) \| (\S+)$")
PROBE_COLS = ("P1 +y wall", "P2 -y wall")
FIELDS = (("1.5 T", "59", ENTRY59_GLOB, ENTRY59_NODE_PF),
          ("3 T", "60", ENTRY60_GLOB, ENTRY60_NODE_PF))
DISCLAIMERS = ("degree 1", "1.5 T gel constants", "20 °C baseline", "No container wall",
               "No compliance claim")
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
    nxt = re.search(r"^## ", text[start + len(HEADING):], re.M)
    return text[start: start + len(HEADING) + nxt.start()] if nxt else text[start:]


def table(sec: str) -> dict[tuple[str, str], list[str]]:
    """``{(field, sense): [c*, θ, P1, P2]}`` from the section's table rows."""
    rows = {}
    for line in sec.splitlines():
        line = line.strip()
        if not line.startswith("|") or set(line) <= set("|-: "):
            continue
        cells = [c.strip() for c in line.strip("|").split("|")]
        if cells[0].startswith("Field"):
            continue
        field, _, sense = cells[0].rpartition(", ")
        if (field, sense) in rows:
            raise ValueError(f"row {cells[0]!r} repeated")
        rows[(field, sense)] = cells[1:]
    return rows


def blocks() -> dict[str, dict]:
    """``{"1.5 T": block, "3 T": block}`` from the tracked logs (both must be tracked)."""
    out = {}
    for label, key, glob, _node in FIELDS:
        hit = _on_root(find_log, glob)
        assert hit is not None, f"entry {key}: no {glob}"
        path, kind = hit
        print(f"{TAG} entry {key} ({label}): {path.name} ({kind})")
        assert kind == "tracked", (key, kind)
        out[label] = _on_root(parse_heat_log, path)
    return out


def expected_rows(bl: dict[str, dict]) -> dict[tuple[str, str], list[str]]:
    """The page table as the parser reads the logs (strings)."""
    exp = {}
    for label, *_ in FIELDS:
        for sn in SENSES:
            p61 = bl[label]["senses"][sn]["p61"]
            exp[(label, sn)] = [p61["hot"]["point"], p61["hot"]["theta_point"],
                                *(p61["probes"][c]["theta"] for c in PROBE_COLS)]
    r = ratios(side(bl["1.5 T"]), side(bl["3 T"]))
    for sn in SENSES:
        exp[("3 T / 1.5 T", sn)] = ["each its own", format(r[sn][RATIO_KEY], ".6e"),
                                    "not read", "not read"]
    return exp


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


def c39_tuples() -> dict[str, set[tuple[str, str, str]]]:
    out = {sn: set() for sn in SENSES}
    n = {sn: 0 for sn in SENSES}
    for ln in _on_root(_read_log, C39_LOG).splitlines():
        m = C39_RE.match(ln)
        if m:
            out[m.group(1)].add(m.groups()[1:])
            n[m.group(1)] += 1
    print(f"{TAG} C39 log {C39_LOG.name}: matching lines per sense {n}")
    return out


def _bl_cache():
    if not hasattr(_bl_cache, "v"):
        _bl_cache.v = blocks()
    return _bl_cache.v


def _sec():
    return section(PAGE.read_text(encoding="utf-8"))


# ---------------------------------------------------------------------------
# Gates.
# ---------------------------------------------------------------------------
def test_s46_chk_c46_checks_on_tracked_logs():
    bl = _bl_cache()
    names = list(_probes())
    for label, key, _glob, node in FIELDS:
        b = bl[label]
        assert b["head"]["degree"] == "1" and b["head"]["senses"] == list(SENSES)
        require_block(b, node, names, f"{label} (entry {key})")


def test_s46_q_tab_cells_are_the_parsed_tokens():
    bl = _bl_cache()
    sec = _sec()
    exp = expected_rows(bl)
    for k, v in exp.items():
        print(f"{TAG} (Q-tab) {k[0]}, {k[1]}: parsed {v}")
    bad = tab_mismatches(sec, exp)
    print(f"{TAG} (Q-tab) {len(exp)} rows x 4 cells compared; mismatches {bad}")
    assert bad == []
    # provenance sentence: node and t* strings
    for label, *_ in FIELDS:
        node = bl[label]["head"]["node"]
        assert f"{node} pF" in sec, node
        for sn in SENSES:
            t_star = bl[label]["senses"][sn]["i32"]["cn"]["t_star"]
            assert f"t\\* = {t_star} s" in sec, t_star
            print(f"{TAG} (Q-tab) {label} {sn}: node {node} pF, t* {t_star} s on the page")


def test_s46_q_ratio_c46_and_c39_custody():
    bl = _bl_cache()
    r = ratios(side(bl["1.5 T"]), side(bl["3 T"]))
    got = table(_sec())
    tup = c39_tuples()
    for sn in SENSES:
        mine = format(r[sn][RATIO_KEY], ".6e")
        assert len(tup[sn]) == 1, (sn, tup[sn])
        (t15, t3, rr), = tup[sn]
        print(f"{TAG} (Q-ratio) {sn}: C39 tuple ({t15} | {t3} | {rr}); C46 ratios {mine}; "
              f"page {got[('3 T / 1.5 T', sn)][1]}")
        assert t15 == bl["1.5 T"]["senses"][sn]["p61"]["hot"]["theta_point"]
        assert t3 == bl["3 T"]["senses"][sn]["p61"]["hot"]["theta_point"]
        assert rr == mine == got[("3 T / 1.5 T", sn)][1]
    assert C39_LOG.name in _sec()


def test_s46_q_all_no_other_decimal():
    bl = _bl_cache()
    sec = _sec()
    exp = set()
    for v in expected_rows(bl).values():
        exp.update(c for c in v if DECIMAL_RE.search(c))
    for label, *_ in FIELDS:
        exp.add(bl[label]["head"]["node"])
        exp.update(bl[label]["senses"][sn]["i32"]["cn"]["t_star"] for sn in SENSES)
    rest = sec
    for s in sorted(exp, key=len, reverse=True):
        cnt = sec.count(s)
        print(f"{TAG} (Q-all) {s!r}: {cnt} occurrence(s)")
        assert cnt >= 1, s
        rest = rest.replace(s, "")
    stray = [t for t in DECIMAL_RE.findall(rest) if t != "1.5"]
    print(f"{TAG} (Q-all) {len(exp)} expected strings; stray decimals {stray}")
    assert stray == []


def test_s46_disc_disclaimers():
    page = PAGE.read_text(encoding="utf-8")
    sec = section(page)
    for d in DISCLAIMERS:
        print(f"{TAG} (Disc) {d!r}: {sec.count(d)}")
        assert d in sec, d
    assert "AED" not in page and "HFSS" not in page
    assert "~" not in sec


def test_s46_d_docrefs_guide():
    rc, counts = _on_root(run_checker, "--docs-root", "docs/guide", "--check-links")
    link = LINK_PASS_RE.search(_on_root(_checker_stdout, "--docs-root", "docs/guide",
                                        "--check-links"))
    assert link is not None, "(D) no 'Link pass:' line"
    print(f"{TAG} (D) dead={counts['dead']} links={link.group(1)} pages={link.group(2)} "
          f"link-dead={link.group(3)} rc={rc}")
    assert counts["dead"] == 0
    assert int(link.group(3)) == 0


# ---------------------------------------------------------------------------
# Negative controls.
# ---------------------------------------------------------------------------
def test_s46_c_digit_one_mismatch():
    bl = _bl_cache()
    text = PAGE.read_text(encoding="utf-8")
    tok = bl["3 T"]["senses"]["ccw"]["p61"]["hot"]["theta_point"]
    j = tok.index("e") - 1
    new = tok[:j] + str((int(tok[j]) + 1) % 10) + tok[j + 1:]
    assert text.count(tok) == 1
    planted = text.replace(tok, new, 1)
    bad = tab_mismatches(section(planted), expected_rows(bl))
    print(f"{TAG} (c-digit) {tok} -> {new}: {bad}")
    assert bad == [(("3 T", "ccw"), 1, new, tok)]


def test_s46_c_node_only_node_fails():
    bl = _bl_cache()
    res = check_block(bl["3 T"], ENTRY59_NODE_PF, list(_probes()), "c-node")
    failed = [n for n, ok, _d in res if not ok]
    print(f"{TAG} (c-node) entry 60 against {ENTRY59_NODE_PF}: failed {failed} of {len(res)}")
    assert failed == [f"c-node (Node) node token == {ENTRY59_NODE_PF}"]


def test_s46_c_sec_missing_or_doubled_heading():
    text = PAGE.read_text(encoding="utf-8")
    assert text.count(HEADING) == 1
    with pytest.raises(LookupError):
        section(text.replace(HEADING, "## Something else"))
    with pytest.raises(LookupError):
        section(text.replace(HEADING, HEADING + "\n\n" + HEADING))
    print(f"{TAG} (c-sec) removed and doubled heading both raise LookupError")
