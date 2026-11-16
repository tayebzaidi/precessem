"""`WF-17` step 50 / U160 (§9 item 975): the clinical page's closed-form section is the logs.

``docs/guide/clinical_coil_fixtures.md`` section "Closed-form gel-box records (MAT-21
steps 48 / 49)" quotes `MAT-21` step 48 (P88: the thermal-dose excess over the
unheated baseline, per deposited watt, against its linear term) and step 49 (P89':
the 18-ball lattice's SP, the predicted 3 T / 1.5 T lever, the verdict) beside
`GEO-38` C39's three printed ratios.  Every quoted digit is compared, as a string,
with the group of an anchored regex over the tracked log; each regex must yield
exactly one distinct group-tuple (count printed; 10:30 ruling of 2026-11-12, never
"the first match").  Knobless, no mesh, no solve; imports no MAT-21 module and
neither sibling page module.

Gates: (Log) four tracked logs, footer ``- Status: 0``; (Q-tab) Table A / B cells,
the per-watt sentence and the verdict == the regex groups; (Q-cust) SP and share
custody (P89' :75 vs :76, P83 :77 vs :76), C39 tokens == C39's anchored third group,
ccw point == the heating section's ``3 T / 1.5 T, ccw`` cell; (Q-all) no other
decimal in the section but ``1.5``; (Disc) disclaimer tokens, no "AED" / "HFSS",
no ``~``; (D) docrefs dead 0, link-dead 0, rc != 1.  Controls: (c-digit), (c-dup),
(c-stub), (c-sec).
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
from tests.log_records import read_log_text, resolve_log_path  # noqa: E402

TAG = "[WF-17 s50]"
HEADING = "## Closed-form gel-box records (MAT-21 steps 48 / 49)"
HEAT_HEADING = "## Clinical heating records (entries 59 / 60)"
LOGS = REPO / "docs" / "testing" / "logs"
P88_LOG = LOGS / "20261112T120219Z_MAT-21-P88-W1.log"
P89_LOG = LOGS / "20261112T182020Z_MAT-21-P89-W1b.log"
P83_LOG = LOGS / "20261110T103343Z_MAT-21-P83-W1.log"
C39_LOG = LOGS / "20261111T103024Z_GEO-38-C39-W1.log"
ALL_LOGS = (P88_LOG, P89_LOG, P83_LOG, C39_LOG)

STATUS_RE = r"^- Status: (\S+)$"
P_DEP = ("0.1", "1", "10", "100")


def p88_row_re(p: str) -> str:
    return (r"^MAT21-48 TABLE " + re.escape(p)
            + r" \| (?:\S+ \| ){5}(\S+) \| \S+ \| \S+ \| (\S+) \| \(")


PERWATT_RE = r"^MAT21-48 per-watt excess at 100 W / at 0\.1 W (\S+) \("
SP_RE = r"^MAT21-49 SP = theta3\(c\*3\)/theta1\.5\(c\*1\.5\) per dep\. W (\S+) \("
LEVER_RE = r"^MAT21-49 LEVER predicted 3 T / 1\.5 T = share (\S+) x SP (\S+) = (\S+) \("
VERDICT_RE = r"^MAT21-49 VERDICT \(printed, the weekly reads it; no gate\): (.+?) \(SP - 1 "
P83_SHARE_RE = r"^MAT21-46 power share s3/s1\.5 (\S+) \("


def p89_beside_re(label: str) -> str:
    return (r"^MAT21-49 LEVER beside C39 " + re.escape(label)
            + r" \(`[^`]+`, \d+ matches, \d+ distinct; [^)]*\): (\S+)")


def c39_re(label: str) -> str:
    return (r"^\[GEO-38 C46\]   (" + re.escape(label)
            + r") theta per P_acc[^:]*: (\S+) \| (\S+) \| (\S+)$")


C39_UNANCHORED = r"(ccw hotspot point) theta per P_acc[^:]*: (\S+) \| (\S+) \| (\S+)$"
# Table B: (row label, P89' beside label or None, C39 label or None)
C39_ROWS = (("C39 ccw point", "ccw point", "ccw hotspot point"),
            ("C39 ccw ball", "ccw ball", "ccw ball"),
            ("C39 cw point (cw heated at ccw's c\\*)", "cw point", "cw hotspot point"))
DISCLAIMERS = ("Closed form", "Not the clinical SAR", "Per deposited watt",
               "1.5 T gel constants", "No dose per watt", "A record")
DECIMAL_RE = re.compile(r"\d+\.\d+")


# ---------------------------------------------------------------------------
# Pure readers.
# ---------------------------------------------------------------------------
def section(text: str, heading: str) -> str:
    """The heading line to the next ``## `` heading; ``LookupError`` unless the
    heading occurs exactly once."""
    heads = [m.start() for m in re.finditer(r"^" + re.escape(heading) + r"[ \t]*$", text, re.M)]
    if len(heads) != 1:
        raise LookupError(f"{heading!r}: {len(heads)} occurrences (exactly one required)")
    start = heads[0]
    nxt = re.search(r"^## ", text[start + len(heading):], re.M)
    return text[start: start + len(heading) + nxt.start()] if nxt else text[start:]


def read_one_text(text: str, regex: str, name: str = "") -> tuple[tuple[str, ...], int]:
    """(groups, count): raises unless the distinct group-tuples number exactly 1
    and the count >= 1.  Prints the count."""
    hits = [m.groups() for m in re.finditer(regex, text, re.M)]
    distinct = set(hits)
    print(f"{TAG} read_one {name}: {len(hits)} match(es), {len(distinct)} distinct "
          f"{sorted(distinct)}")
    if len(hits) < 1 or len(distinct) != 1:
        raise ValueError(f"{name} {regex!r}: {len(hits)} matches, {len(distinct)} distinct")
    return hits[0], len(hits)


def read_one(path: Path, regex: str) -> tuple[tuple[str, ...], int]:
    return read_one_text(read_log_text(path, errors="strict"), regex, path.name)


def table_rows(sec: str) -> dict[str, list[str]]:
    """``{first cell: other cells}`` over every table row of the section."""
    rows = {}
    for line in sec.splitlines():
        line = line.strip()
        if not line.startswith("|") or set(line) <= set("|-: "):
            continue
        cells = [c.strip() for c in line.strip("|").split("|")]
        if cells[0] in ("P_dep (W)", "Quantity"):
            continue
        if cells[0] in rows:
            raise ValueError(f"row {cells[0]!r} repeated")
        rows[cells[0]] = cells[1:]
    return rows


def expected() -> dict:
    """Every quoted value, read from the logs by its named regex."""
    exp = {"A": {}, "B": {}}
    for p in P_DEP:
        g, _n = read_one(P88_LOG, p88_row_re(p))
        exp["A"][p] = list(g)
    exp["perwatt"] = read_one(P88_LOG, PERWATT_RE)[0][0]
    sp = read_one(P89_LOG, SP_RE)[0][0]
    share, lever_sp, lever = read_one(P89_LOG, LEVER_RE)[0]
    exp["lever_sp"] = lever_sp
    exp["B"]["SP"] = [sp]
    exp["B"]["power share (3 T / 1.5 T)"] = [share]
    exp["B"]["predicted lever = share × SP"] = [lever]
    for row, beside, _c39 in C39_ROWS:
        exp["B"][row] = [read_one(P89_LOG, p89_beside_re(beside))[0][0]]
    exp["verdict"] = read_one(P89_LOG, VERDICT_RE)[0][0]
    return exp


def q_tab_mismatches(sec: str, exp: dict) -> list[tuple]:
    got = table_rows(sec)
    want = {**exp["A"], **exp["B"]}
    bad = [(k, None, got[k]) for k in got if k not in want]
    for k, w in want.items():
        h = got.get(k)
        if h is None:
            bad.append((k, None, w))
            continue
        if h != w:
            bad.append((k, h, w))
    sentence = f"at 100 W is {exp['perwatt']} times that at 0.1 W"
    if sentence not in " ".join(sec.split()):
        bad.append(("perwatt", exp["perwatt"]))
    if f'"{exp["verdict"]}"' not in sec:
        bad.append(("verdict", exp["verdict"]))
    return bad


def _exp_cache():
    if not hasattr(_exp_cache, "v"):
        _exp_cache.v = expected()
    return _exp_cache.v


def _page():
    return PAGE.read_text(encoding="utf-8")


def _sec():
    return section(_page(), HEADING)


# ---------------------------------------------------------------------------
# Gates.
# ---------------------------------------------------------------------------
def test_s50_log_status_zero():
    for p in ALL_LOGS:
        assert resolve_log_path(p).is_file(), p
        (st,), n = read_one(p, STATUS_RE)
        print(f"{TAG} (Log) {p.name}: Status {st} (count {n})")
        assert st == "0", (p.name, st)


def test_s50_q_tab_cells_are_the_log_groups():
    exp = _exp_cache()
    sec = _sec()
    print(f"{TAG} (Q-tab) Table A parsed {exp['A']}")
    print(f"{TAG} (Q-tab) Table B parsed {exp['B']}")
    print(f"{TAG} (Q-tab) per-watt {exp['perwatt']}; verdict {exp['verdict']!r}")
    got = table_rows(sec)
    print(f"{TAG} (Q-tab) page rows {got}")
    bad = q_tab_mismatches(sec, exp)
    print(f"{TAG} (Q-tab) {len(exp['A'])} + {len(exp['B'])} rows compared; mismatches {bad}")
    assert bad == []


def test_s50_q_cust_custody():
    exp = _exp_cache()
    sp = exp["B"]["SP"][0]
    print(f"{TAG} (Q-cust) SP :75 {sp} vs lever line {exp['lever_sp']}")
    assert sp == exp["lever_sp"]
    p83 = read_one(P83_LOG, P83_SHARE_RE)[0][0]
    share = exp["B"]["power share (3 T / 1.5 T)"][0]
    print(f"{TAG} (Q-cust) P83 share {p83} vs lever line share {share}")
    assert p83 == share
    for row, _beside, lab in C39_ROWS:
        g, n = read_one(C39_LOG, c39_re(lab))
        print(f"{TAG} (Q-cust) C39 {lab}: {g} count {n}; page {exp['B'][row][0]}")
        assert n == 2, (lab, n)
        assert g[3] == exp["B"][row][0], (lab, g, exp["B"][row])
    heat = table_rows_heating(_page())
    cell = heat["3 T / 1.5 T, ccw"]
    print(f"{TAG} (Q-cust) heating section '3 T / 1.5 T, ccw' cell {cell}")
    assert cell == exp["B"]["C39 ccw point"][0]


def table_rows_heating(page: str) -> dict[str, str]:
    """``{row label: θ(t*) at c* cell}`` from the heating section (own reader)."""
    out = {}
    for line in section(page, HEAT_HEADING).splitlines():
        line = line.strip()
        if not line.startswith("|") or set(line) <= set("|-: "):
            continue
        cells = [c.strip() for c in line.strip("|").split("|")]
        out[cells[0]] = cells[2]
    return out


def _expected_strings(exp) -> set[str]:
    """The decimal-bearing quoted strings (P_dep ``0.1`` included) and the verdict;
    the integer P_dep labels are not removed, so they cannot mask a stray decimal."""
    s = set(P_DEP) | {exp["perwatt"], exp["verdict"]}
    for v in exp["A"].values():
        s.update(v)
    for v in exp["B"].values():
        s.update(v)
    return {t for t in s if DECIMAL_RE.search(t) or t == exp["verdict"]}


def test_s50_q_all_no_other_decimal():
    exp = _exp_cache()
    sec = _sec()
    rest = sec
    strings = _expected_strings(exp)
    for s in sorted(strings, key=len, reverse=True):
        cnt = sec.count(s)
        print(f"{TAG} (Q-all) {s!r}: {cnt} occurrence(s)")
        assert cnt >= 1, s
        rest = rest.replace(s, "")
    stray = [t for t in DECIMAL_RE.findall(rest) if t != "1.5"]
    print(f"{TAG} (Q-all) {len(strings)} expected strings; stray decimals {stray}")
    assert stray == []


def test_s50_disc_disclaimers():
    page = _page()
    sec = section(page, HEADING)
    for d in DISCLAIMERS:
        print(f"{TAG} (Disc) {d!r}: {sec.count(d)}")
        assert d in sec, d
    assert "AED" not in page and "HFSS" not in page
    assert "~" not in sec
    assert "cem43" not in sec


def test_s50_d_docrefs_guide():
    rc, counts = run_checker("--docs-root", "docs/guide", "--check-links")
    link = LINK_PASS_RE.search(_checker_stdout("--docs-root", "docs/guide", "--check-links"))
    assert link is not None, "(D) no 'Link pass:' line"
    print(f"{TAG} (D) dead={counts['dead']} links={link.group(1)} pages={link.group(2)} "
          f"link-dead={link.group(3)} rc={rc}")
    assert counts["dead"] == 0
    assert int(link.group(3)) == 0
    assert rc != 1


# ---------------------------------------------------------------------------
# Negative controls.
# ---------------------------------------------------------------------------
def test_s50_c_digit_one_mismatch():
    exp = _exp_cache()
    text = _page()
    tok = "4.078321"
    assert text.count(tok) == 1
    planted = text.replace(tok, "4.078322", 1)
    bad = q_tab_mismatches(section(planted, HEADING), exp)
    print(f"{TAG} (c-digit) {tok} -> 4.078322: {bad}")
    assert bad == [("predicted lever = share × SP", ["4.078322"], [tok])]


def test_s50_c_dup_doubled_and_changed():
    lines = read_log_text(P89_LOG, errors="strict").splitlines()
    lever = [ln for ln in lines if re.match(LEVER_RE, ln)]
    assert len(lever) == 1
    base = "\n".join(lines)
    g1, n1 = read_one_text(base, LEVER_RE, "c-dup base")
    g2, n2 = read_one_text(base + "\n" + lever[0], LEVER_RE, "c-dup doubled")
    print(f"{TAG} (c-dup) doubled: {g2} count {n2} (base {g1} count {n1})")
    assert g2 == g1 and n1 == 1 and n2 == 2
    changed = lever[0].replace("4.078321", "4.078331")
    assert changed != lever[0]
    with pytest.raises(ValueError):
        read_one_text(base + "\n" + changed, LEVER_RE, "c-dup changed")
    print(f"{TAG} (c-dup) one-digit-changed copy raises (two distinct)")


def test_s50_c_stub_unanchored_sees_two():
    lines = read_log_text(C39_LOG, errors="strict").splitlines()
    unanch = {m.groups() for ln in lines if (m := re.search(C39_UNANCHORED, ln))}
    print(f"{TAG} (c-stub) unanchored: {len(unanch)} distinct {sorted(unanch)}")
    assert len(unanch) == 2
    g, n = read_one(C39_LOG, c39_re("ccw hotspot point"))
    print(f"{TAG} (c-stub) anchored: 1 distinct {g}, count {n}")


def test_s50_c_sec_missing_or_doubled_heading():
    text = _page()
    assert text.count(HEADING) == 1
    with pytest.raises(LookupError):
        section(text.replace(HEADING, "## Something else"), HEADING)
    with pytest.raises(LookupError):
        section(text.replace(HEADING, HEADING + "\n\n" + HEADING), HEADING)
    print(f"{TAG} (c-sec) removed and doubled heading both raise LookupError")
