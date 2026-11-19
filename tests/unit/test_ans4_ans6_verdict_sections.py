"""`ANS-4` / `ANS-6` public verdict sections (§9 item 416, §10 F33) — render-only.

On scratch copies of each case's committed ``metrics.json`` (neither case is re-run):

(Pre4)/(Pre6) the pinned pre-change writer (``5d966d96`` / ``69fa0777``, copied by a host
              ``git show`` to ``logs/f33_pin/<repo-relative path>``) reproduces the pre-edit
              ``COMPARISON.md`` (copied from ``0cdd1cb9`` to the same tree) byte for byte —
              the premise that the committed file is the writer's output;
(R4)/(R6)     the edited writer reproduces the committed ``COMPARISON.md`` byte for byte, the
              verdict heading exactly once;
(Δ)           a line-level difflib from the pre-edit to the committed file touches only the
              scope-paragraph lines and the verdict section;
(V)           each verdict carries its levels, each exactly once inside the section;
(N-digit)     one verdict digit changed in a copy of the committed file differs from the render;
(N-old)       the pinned render has verdict-heading count 0 and differs from the committed file.

The pinned source is executed with ``__file__`` set to the live script's path, so its
``parents[3]`` repo-root resolution and sibling-example loads are the live tree's.
"""

from __future__ import annotations

import difflib
import json
import re
import sys
import types
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
PIN = REPO / "logs" / "f33_pin"
BENCH = REPO / "examples" / "ansys_benchmarks"
HEADING = "## Verdict (public; adjudicated by the weekly review)"

CASES = {
    "4": dict(
        case=BENCH / "ans4_birdcage_four_port_10_64_128MHz",
        script="04_birdcage_four_port_10_64_128MHz.py",
        pin_sha="5d966d96",
        levels=("0.73 %", "1 %", "0.1 %", "0.014 / 0.057 / 0.066 %"),
        digit=("0.73 %", "0.74 %"),
        # the stale clause's first and last words, pre-edit (scope paragraph)
        scope_old=(
            "absolute-accuracy claim is made from our side; adjudication is the next weekly",
            "review's, after the AED numbers land.",
        ),
        scope_new=(
            "absolute-accuracy claim is made from our side; the adjudication is recorded in",
            "the verdict section below.",
        ),
    ),
    "6": dict(
        case=BENCH / "ans6_copper_birdcage_four_port_10_64_128MHz",
        script="06_copper_birdcage_four_port_10_64_128MHz.py",
        pin_sha="69fa0777",
        levels=("0.33 %", "0.21 %"),
        digit=("0.33 %", "0.34 %"),
        scope_old=(
            "absolute `S₁₁` / `Z_in` claim at 10 or 64 MHz (`PORT-21`). Adjudication is a",
            "future weekly review's, after the AED numbers land.",
        ),
        scope_new=(
            "absolute `S₁₁` / `Z_in` claim at 10 or 64 MHz (`PORT-21`). The adjudication",
            "is recorded in the verdict section below.",
        ),
    ),
}


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


def _count_level(text: str, level: str) -> int:
    # "1 %" must not match inside "0.1 %": no digit or dot before the level.
    return len(re.findall(r"(?<![\d.])" + re.escape(level), text))


def _verdict(text: str) -> str:
    start = text.index(HEADING)
    nxt = text.find("\n## ", start + len(HEADING))
    return text[start:nxt]


@pytest.mark.parametrize("tag", ["4", "6"])
def test_verdict_section(tag, tmp_path):
    c = CASES[tag]
    case, script = c["case"], c["script"]
    rel = (case / script).relative_to(REPO)
    pin_script = PIN / rel
    pin_comparison = PIN / rel.parent / "COMPARISON.md"
    for p in (pin_script, pin_comparison):
        if not p.exists():
            pytest.fail(f"{p} missing: pin it by host `git show` first (§9 item 416 window 1)")
    committed = (case / "COMPARISON.md").read_bytes()
    pre_edit = pin_comparison.read_bytes()
    metrics_before = (case / "metrics.json").read_bytes()

    # (Pre)
    old_mod = _exec_as(f"ans{tag}_pin_{c['pin_sha']}_ut", pin_script, case / script)
    (tmp_path / "old").mkdir()
    old = _render(old_mod, case, tmp_path / "old")
    print(
        f"\n(Pre{tag}) pinned {c['pin_sha']} render {len(old)} B vs pre-edit {len(pre_edit)} B: "
        f"identical={old == pre_edit}"
    )
    if old != pre_edit:
        for i, (a, b) in enumerate(zip(old.decode().splitlines(), pre_edit.decode().splitlines())):
            if a != b:
                print(f"  first differing line {i + 1}: render {a!r} / pre-edit {b!r}")
                break
    assert old == pre_edit

    # (R)
    new_mod = _exec_as(f"ans{tag}_live_ut", case / script, case / script)
    (tmp_path / "new").mkdir()
    new = _render(new_mod, case, tmp_path / "new")
    n_head = new.decode().count(HEADING)
    print(
        f"(R{tag}) render {len(new)} B vs committed {len(committed)} B: identical={new == committed}; "
        f"verdict heading count {n_head}"
    )
    assert n_head == 1
    assert new == committed

    # (R0) the repo files are untouched by the test
    assert (case / "COMPARISON.md").read_bytes() == committed
    assert (case / "metrics.json").read_bytes() == metrics_before

    # (Δ) the diff touches only the scope clause and the verdict section
    old_lines = pre_edit.decode().splitlines()
    committed_text = committed.decode()
    # §9 item 417 (ruling (2) of the 2026-10-13 03:00 review): `ans:4`'s writer now
    # appends the h-sequence section after Provenance. It is item 417's diff, not
    # item 416's — gated byte for byte by tests/unit/test_ans4_h_table.py — so it is
    # cut here, at its heading (present exactly once, as the file's tail), before the
    # item-416 diff is classified. Nothing else is exempted.
    if tag == "4":
        h_head = (
            "\n## h-sequence of the adjudicated pairing (ours; degree 2, exact rectangles, κ off)\n"
        )
        assert committed_text.count(h_head) == 1
        committed_text = committed_text[: committed_text.index(h_head) + 1]
    # §9 item 462 (licensed, 2026-10-15 03:00 review): `ans:6`'s writer now appends
    # the public Palace section after Provenance — item 462's diff, gated byte for
    # byte by tests/unit/test_ans4_ans6_palace_sections.py — so it is cut here at its
    # heading (present exactly once), exactly as tag "4" above. Nothing else exempted.
    if tag == "6":
        p_head = "\n## Palace at matched order (public; `OPS-61` steps 2a–2c)\n"
        assert committed_text.count(p_head) == 1
        committed_text = committed_text[: committed_text.index(p_head) + 1]
    new_lines = committed_text.splitlines()
    text = committed.decode()
    verdict = _verdict(text)
    v_lines = set(verdict.splitlines())
    sm = difflib.SequenceMatcher(a=old_lines, b=new_lines, autojunk=False)
    removed, added = [], []
    for op, i1, i2, j1, j2 in sm.get_opcodes():
        if op != "equal":
            removed += old_lines[i1:i2]
            added += new_lines[j1:j2]
    scope_first, scope_last = c["scope_old"]
    i0 = old_lines.index(scope_first)
    i1 = old_lines.index(scope_last)
    allowed_removed = set(old_lines[i0 : i1 + 1])
    # added lines: the re-pointed scope clause (two lines) and verdict-section lines
    stray_added = [
        ln for ln in added if ln not in v_lines and ln != "" and ln not in c["scope_new"]
    ]
    stray_removed = [ln for ln in removed if ln not in allowed_removed]
    print(
        f"(Δ{tag}) removed {len(removed)} line(s), added {len(added)} line(s); "
        f"stray removed {stray_removed!r}, stray added {stray_added!r}"
    )
    assert removed and set(removed) <= allowed_removed
    assert not stray_added
    assert all(ln in added for ln in c["scope_new"])

    # (V) the levels, each exactly once inside the section
    counts = {lvl: _count_level(verdict, lvl) for lvl in c["levels"]}
    print(f"(V{tag}) level counts inside the verdict: {counts}")
    assert all(n == 1 for n in counts.values())
    assert "AED (First Order)" not in verdict  # the section quotes no AED column

    # (N-digit)
    a, b = c["digit"]
    planted = text.replace(a, b, 1).encode()
    assert planted != text.encode()
    print(f"(N-digit{tag}) {a!r} -> {b!r}: planted file differs from render = {planted != new}")
    assert planted != new

    # (N-old)
    n_old = old.decode().count(HEADING)
    print(
        f"(N-old{tag}) pinned render verdict-heading count {n_old}; differs from committed = "
        f"{old != committed}"
    )
    assert n_old == 0
    assert old != committed
