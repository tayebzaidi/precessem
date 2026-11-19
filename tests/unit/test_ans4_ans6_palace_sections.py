"""`ANS-4` / `ANS-6` public "Palace at matched order" sections (§9 item 462, §10 F38).

Render-only, on scratch copies of each case's committed ``metrics.json`` (neither case
is re-run; ``_write_comparison(private=True)`` is never called — it reads
``aed_results/``):

(X)       every Palace / ours figure in either section is traced: the cell string occurs
          verbatim on the cited line of the tracked log, and in the committed section;
(R4)/(R6) the live writer, public mode, reproduces the committed ``COMPARISON.md`` byte
          for byte, the Palace heading exactly once, as the file's tail;
(P)       ``PALACE_SECTION`` is referenced in ``_write_comparison`` only under the
          ``if not private`` branch (AST check), so the private bytes cannot change;
(B)       no bridge: neither section names AED except the closing "No AED value.", nor
          ``aed_results`` / ``COMPARISON_private`` / ``docs/private``;
(R0)      the repo files are unchanged by the test.

Negative controls (asserted, by construction): (N-cite) a cell cited on the wrong line
fails (X); (N-digit) one Palace digit changed in a copy of the committed file differs from
the render; (N-old) the pinned pre-change writers (``4989ebdd`` / ``c7b1b682``, copied by
host ``git show`` to ``logs/f38_pin/<repo-relative path>``) render Palace-heading count 0
and differ from the committed files.
"""

from __future__ import annotations

import ast
import gzip
import inspect
import json
import sys
import textwrap
import types
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
LOGS = REPO / "docs" / "testing" / "logs"
PIN = REPO / "logs" / "f38_pin"
BENCH = REPO / "examples" / "ansys_benchmarks"

CASES = {
    "4": dict(
        case=BENCH / "ans4_birdcage_four_port_10_64_128MHz",
        script="04_birdcage_four_port_10_64_128MHz.py",
        pin_sha="4989ebdd",
        heading="## Palace at matched order (public; `OPS-61` F21 / F34 / F36)",
        digit=("1.917e-10", "1.918e-10"),
    ),
    "6": dict(
        case=BENCH / "ans6_copper_birdcage_four_port_10_64_128MHz",
        script="06_copper_birdcage_four_port_10_64_128MHz.py",
        pin_sha="c7b1b682",
        heading="## Palace at matched order (public; `OPS-61` steps 2a–2c)",
        digit=("2.197e-10", "2.198e-10"),
    ),
}

#: (tag, cell, log, 1-indexed line) — every figure of the two sections.
TRACE = [
    # ans:6
    ("6", "1.612e-10", "20261005T093623Z_OPS-61-item241-W6-read.log", 56),
    ("6", "1.586e-10", "20261005T093623Z_OPS-61-item241-W6-read.log", 56),
    ("6", "2.033e-10", "20261005T093623Z_OPS-61-item241-W6-read.log", 56),
    ("6", "1.031e-07", "20261005T093623Z_OPS-61-item241-W6-read.log", 57),
    ("6", "2.263e-07", "20261005T093623Z_OPS-61-item241-W6-read.log", 57),
    ("6", "1.860e-07", "20261005T093623Z_OPS-61-item241-W6-read.log", 57),
    ("6", "1.340e-10", "20261005T022223Z_OPS-61-item234-W5-read.log", 49),
    ("6", "1.120e-10", "20261005T022223Z_OPS-61-item234-W5-read.log", 49),
    ("6", "2.197e-10", "20261005T022223Z_OPS-61-item234-W5-read.log", 49),
    ("6", "1.274e-10", "20261015T004434Z_OPS-61-item454-W7-read.log", 70),
    ("6", "4.521e-11", "20261015T004434Z_OPS-61-item454-W7-read.log", 70),
    ("6", "6.621e-11", "20261015T004434Z_OPS-61-item454-W7-read.log", 70),
    ("6", "0.1206", "20261015T004434Z_OPS-61-item454-W7-read.log", 73),
    ("6", "0.0506", "20261015T004434Z_OPS-61-item454-W7-read.log", 73),
    ("6", "0.1091", "20261015T004434Z_OPS-61-item454-W7-read.log", 73),
    ("6", "0.1193", "20261012T070009Z_ANS-6-step5b.log", 4594),
    ("6", "0.0501", "20261012T070009Z_ANS-6-step5b.log", 4594),
    ("6", "0.1081", "20261012T070009Z_ANS-6-step5b.log", 4594),
    # item 457 (F37 PEC)
    ("6", "1.378e-10", "20261015T095049Z_OPS-61-item457-W6-read.log", 70),
    ("6", "1.153e-10", "20261015T095049Z_OPS-61-item457-W6-read.log", 70),
    ("6", "1.484e-10", "20261015T095049Z_OPS-61-item457-W6-read.log", 70),
    (
        "6",
        "self 0.1193%  adjacent 0.0501%  opposite 0.1081%",
        "20261015T095049Z_OPS-61-item457-W6-read.log",
        73,
    ),
    ("6", "3.581e-08", "20261015T095049Z_OPS-61-item457-W6-read.log", 89),
    ("6", "2.515e-07", "20261015T095049Z_OPS-61-item457-W6-read.log", 89),
    ("6", "1.790e-07", "20261015T095049Z_OPS-61-item457-W6-read.log", 89),
    ("6", "v0.18.1", "20261005T093309Z_OPS-61-item241-W5-palace-hole-cu-o2.log", 45),
    ("6", "80798", "20261005T093309Z_OPS-61-item241-W5-palace-hole-cu-o2.log", 68),
    ("6", "570332", "20261005T093309Z_OPS-61-item241-W5-palace-hole-cu-o2.log", 113),
    ("6", "v0.18.1", "20261015T003418Z_OPS-61-item454-W3-palace-hole-h0075-cu-e1.log", 45),
    ("6", "226904", "20261015T003418Z_OPS-61-item454-W3-palace-hole-h0075-cu-e1.log", 68),
    ("6", "1538030", "20261015T003418Z_OPS-61-item454-W3-palace-hole-h0075-cu-e1.log", 104),
    ("6", "226904", "20261015T094114Z_OPS-61-item457-W2-palace-hole-h0075-pec-e1.log", 68),
    ("6", "1538030", "20261015T094114Z_OPS-61-item457-W2-palace-hole-h0075-pec-e1.log", 101),
    # ans:4
    (
        "4",
        "self 2.464e-10  adjacent 9.668e-11  opposite 9.883e-11",
        "20261004T214342Z_OPS-61-item229-W4-read.log",
        41,
    ),
    (
        "4",
        "self 2.542e-10  adjacent 9.997e-11  opposite 6.236e-11",
        "20261004T214342Z_OPS-61-item229-W4-read.log",
        48,
    ),
    (
        "4",
        "self 9.930e-11  adjacent 8.685e-11  opposite 1.495e-10",
        "20261004T214342Z_OPS-61-item229-W4-read.log",
        60,
    ),
    ("4", "1.917e-10", "20261013T124505Z_OPS-61-item418-W6-read.log", 65),
    ("4", "1.861e-11", "20261013T124505Z_OPS-61-item418-W6-read.log", 65),
    ("4", "7.373e-11", "20261013T124505Z_OPS-61-item418-W6-read.log", 65),
    ("4", "0.2527", "20261013T124505Z_OPS-61-item418-W6-read.log", 67),
    ("4", "0.0447", "20261013T124505Z_OPS-61-item418-W6-read.log", 67),
    ("4", "0.1375", "20261013T124505Z_OPS-61-item418-W6-read.log", 67),
    ("4", "4.757e-11", "20261013T171447Z_OPS-61-item434-W6-read.log", 65),
    ("4", "8.022e-11", "20261013T171447Z_OPS-61-item434-W6-read.log", 65),
    ("4", "3.630e-11", "20261013T171447Z_OPS-61-item434-W6-read.log", 65),
    ("4", "0.2345", "20261013T171447Z_OPS-61-item434-W6-read.log", 67),
    ("4", "0.0950", "20261013T171447Z_OPS-61-item434-W6-read.log", 67),
    ("4", "0.2098", "20261013T171447Z_OPS-61-item434-W6-read.log", 67),
    ("4", "1.275e-10", "20261013T173812Z_OPS-61-item435-W6-read.log", 65),
    ("4", "1.189e-10", "20261013T173812Z_OPS-61-item435-W6-read.log", 65),
    ("4", "7.497e-11", "20261013T173812Z_OPS-61-item435-W6-read.log", 65),
    ("4", "0.2204", "20261013T173812Z_OPS-61-item435-W6-read.log", 67),
    ("4", "0.1590", "20261013T173812Z_OPS-61-item435-W6-read.log", 67),
    ("4", "0.2863", "20261013T173812Z_OPS-61-item435-W6-read.log", 67),
    ("4", "6.308e-11", "20261013T221851Z_OPS-61-item437-W6-read.log", 65),
    ("4", "5.760e-11", "20261013T221851Z_OPS-61-item437-W6-read.log", 65),
    ("4", "1.679e-11", "20261013T221851Z_OPS-61-item437-W6-read.log", 65),
    ("4", "0.0144", "20261013T221851Z_OPS-61-item437-W6-read.log", 84),
    ("4", "0.0088", "20261013T221851Z_OPS-61-item437-W6-read.log", 84),
    ("4", "0.0161", "20261013T221851Z_OPS-61-item437-W6-read.log", 84),
    # F39 (items 458 / 458b)
    ("4", "9.151e-11", "20261015T111439Z_OPS-61-item458b-W3-read.log", 67),
    ("4", "8.137e-11", "20261015T111439Z_OPS-61-item458b-W3-read.log", 67),
    ("4", "1.143e-10", "20261015T111439Z_OPS-61-item458b-W3-read.log", 67),
    (
        "4",
        "self 0.0569% (5.689e-04)  adjacent 0.0235% (2.353e-04)  opposite 0.0611%",
        "20261015T111439Z_OPS-61-item458b-W3-read.log",
        87,
    ),
    ("4", "v0.18.1", "20261004T213356Z_OPS-61-item229-W2-palace-o2-f10.log", 45),
    ("4", "117206", "20261004T213356Z_OPS-61-item229-W2-palace-o2-f10.log", 68),
    ("4", "753062", "20261004T213356Z_OPS-61-item229-W2-palace-o2-f10.log", 110),
    ("4", "v0.18.1", "20261013T123438Z_OPS-61-item418-W2-palace-h0075-e1.log", 45),
    ("4", "282198", "20261013T123438Z_OPS-61-item418-W2-palace-h0075-e1.log", 68),
    ("4", "1822094", "20261013T123438Z_OPS-61-item418-W2-palace-h0075-e1.log", 101),
    ("4", "v0.18.1", "20261013T215659Z_OPS-61-item437-W2-palace-h0050-e1.log", 45),
    ("4", "593478", "20261013T215659Z_OPS-61-item437-W2-palace-h0050-e1.log", 68),
    ("4", "3833208", "20261013T215659Z_OPS-61-item437-W2-palace-h0050-e1.log", 101),
]

#: section strings for cells whose section form differs from the log form
SECTION_FORM = {
    "80798": "80 798",
    "570332": "570 332",
    "226904": "226 904",
    "1538030": "1 538 030",
    "117206": "117 206",
    "753062": "753 062",
    "282198": "282 198",
    "1822094": "1 822 094",
    "593478": "593 478",
    "3833208": "3 833 208",
    "self 0.1193%  adjacent 0.0501%  opposite 0.1081%": "0.1193 / 0.0501 / 0.1081",
    "self 2.464e-10  adjacent 9.668e-11  opposite 9.883e-11": "2.464e-10 / 9.668e-11 / 9.883e-11",
    "self 2.542e-10  adjacent 9.997e-11  opposite 6.236e-11": "2.542e-10 / 9.997e-11 / 6.236e-11",
    "self 9.930e-11  adjacent 8.685e-11  opposite 1.495e-10": "9.930e-11 / 8.685e-11 / 1.495e-10",
    "self 0.0569% (5.689e-04)  adjacent 0.0235% (2.353e-04)  opposite 0.0611%": "0.0569 / 0.0235 / 0.0611 — ours 0.0569 / 0.0235 / 0.0611",
}


def _log_line(name: str, line: int) -> str:
    p = LOGS / name
    if p.exists():
        text = p.read_text(errors="replace")
    else:
        text = gzip.open(LOGS / (name + ".gz"), "rt", errors="replace").read()
    return text.splitlines()[line - 1]


def _exec_as(name: str, source: Path, as_path: Path):
    mod = types.ModuleType(name)
    mod.__file__ = str(as_path)
    sys.modules[name] = mod
    exec(compile(source.read_text(), str(source), "exec"), mod.__dict__)
    return mod


def _render(mod, case: Path, tmp: Path) -> bytes:
    (tmp / "metrics.json").write_bytes((case / "metrics.json").read_bytes())
    old = mod.CASE_DIR
    mod.CASE_DIR = tmp
    try:
        mod._write_comparison(json.loads((case / "metrics.json").read_text()))
    finally:
        mod.CASE_DIR = old
    return (tmp / "COMPARISON.md").read_bytes()


def _section(text: str, heading: str) -> str:
    return text[text.index(heading) :]


def test_traced():
    """(X) and (N-cite)."""
    sections = {
        t: _section((c["case"] / "COMPARISON.md").read_text(), c["heading"])
        for t, c in CASES.items()
    }
    bad = []
    for tag, cell, log, line in TRACE:
        if cell not in _log_line(log, line):
            bad.append(f"log {log}:{line} lacks {cell!r}")
        if SECTION_FORM.get(cell, cell) not in sections[tag]:
            bad.append(f"ans:{tag} section lacks {SECTION_FORM.get(cell, cell)!r}")
    print(f"\n(X) {len(TRACE)} traced cells, {len(bad)} failure(s): {bad!r}")
    assert not bad
    # (N-cite) a cell cited on the wrong line fails the trace
    wrong = _log_line("20261013T124505Z_OPS-61-item418-W6-read.log", 64)
    print(f"(N-cite) 1.917e-10 at 418 W6:64 present = {'1.917e-10' in wrong}")
    assert "1.917e-10" not in wrong


@pytest.mark.parametrize("tag", ["4", "6"])
def test_palace_section(tag, tmp_path):
    c = CASES[tag]
    case, script = c["case"], c["script"]
    rel = (case / script).relative_to(REPO)
    pin_script = PIN / rel
    if not pin_script.exists():
        pytest.fail(f"{pin_script} missing: pin it by host `git show {c['pin_sha']}:{rel}` first")
    committed = (case / "COMPARISON.md").read_bytes()
    metrics_before = (case / "metrics.json").read_bytes()
    h_json = case / "h_table_s.json"
    h_before = h_json.read_bytes() if h_json.exists() else None

    # (R)
    mod = _exec_as(f"ans{tag}_live_palace_ut", case / script, case / script)
    (tmp_path / "new").mkdir()
    new = _render(mod, case, tmp_path / "new")
    text = new.decode()
    n_head = text.count(c["heading"])
    print(
        f"\n(R{tag}) render {len(new)} B vs committed {len(committed)} B: "
        f"identical={new == committed}; Palace heading count {n_head}"
    )
    assert n_head == 1
    assert new == committed
    section = _section(text, c["heading"])
    assert text.endswith(section) and "\n## " not in section[len(c["heading"]) :]
    assert mod.PALACE_SECTION.strip() == section.strip()

    # (P) PALACE_SECTION referenced only under `if not private`
    src = textwrap.dedent(inspect.getsource(mod._write_comparison))
    fn = ast.parse(src).body[0]
    guarded, total = 0, 0
    for node in ast.walk(fn):
        if isinstance(node, ast.Name) and node.id == "PALACE_SECTION":
            total += 1
    for node in ast.walk(fn):
        if (
            isinstance(node, ast.If)
            and isinstance(node.test, ast.UnaryOp)
            and isinstance(node.test.op, ast.Not)
            and isinstance(node.test.operand, ast.Name)
            and node.test.operand.id == "private"
        ):
            guarded += sum(
                1 for n in ast.walk(node) if isinstance(n, ast.Name) and n.id == "PALACE_SECTION"
            )
    print(f"(P{tag}) PALACE_SECTION references {total}, under `if not private` {guarded}")
    assert total >= 1 and guarded == total

    # (B) no bridge
    aed = section.count("AED")
    bad = [s for s in ("aed_results", "COMPARISON_private", "docs/private") if s in section]
    print(f"(B{tag}) 'AED' count {aed} (only 'No AED value.'), forbidden tokens {bad!r}")
    assert aed == 1 and section.rstrip().endswith("No AED value.")
    assert not bad
    for ext in (".json", ".msh", ".csv"):
        assert ext not in section

    # (R0)
    assert (case / "COMPARISON.md").read_bytes() == committed
    assert (case / "metrics.json").read_bytes() == metrics_before
    if h_before is not None:
        assert h_json.read_bytes() == h_before

    # (N-digit)
    a, b = c["digit"]
    planted = text.replace(a, b, 1).encode()
    assert planted != new
    print(f"(N-digit{tag}) {a!r} -> {b!r}: planted file differs from render = {planted != new}")

    # (N-old)
    old_mod = _exec_as(f"ans{tag}_pin_{c['pin_sha']}_palace_ut", pin_script, case / script)
    (tmp_path / "old").mkdir()
    old = _render(old_mod, case, tmp_path / "old")
    n_old = old.decode().count(c["heading"])
    print(
        f"(N-old{tag}) pinned {c['pin_sha']} render Palace-heading count {n_old}; "
        f"differs from committed = {old != committed}; pre-section prefix equal = "
        f"{committed.startswith(old)}"
    )
    assert n_old == 0
    assert old != committed
    assert committed == old + mod.PALACE_SECTION.encode()
