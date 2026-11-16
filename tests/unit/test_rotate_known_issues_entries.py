"""OPS-66: the known-issues rotation finds undated RETIRED headings.

Before OPS-66, ``rotate_known_issues`` took an entry to be a ``###`` heading or a
dated ``## YYYY-MM-DD`` one, so four retired entries headed ``## RETIRED <date>`` /
``## ✅ RETIRED <date>`` never rotated. The entry pattern is now the module
constant ``KI_ENTRY_RE``.

Anchors (§9 item 139):
  A1  the four headings, verbatim, are entries and are closed;
  A2  over every non-fenced heading of the live file and the archive, old ⊆ new
      and new − old ⊆ {closed headings}; new − old is exactly the four;
  A3  no section header becomes an entry; the struck-through "🟡 OPEN" clauses
      inside two of the four stay closed;
  NC  (negative control, binary) the pinned pre-change tool's dry-run on a copy
      of the live file reports none of the four, the new tool all four.

Durability: the four are moved to the archive by the 2026-10-03 weekly. The
"live file" for A2's count and NC is therefore the working known-issues.md while
it still carries all four, else known-issues.md as pinned at KI_PINNED_SHA (its
last commit before OPS-66). No live line numbers are hard-coded.
"""
from __future__ import annotations

import importlib.util
import re
import subprocess
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
TOOL = REPO_ROOT / "scripts" / "maintenance" / "rotate_plan_archive.py"
KI = REPO_ROOT / "docs" / "testing" / "known-issues.md"
KI_ARCHIVE = REPO_ROOT / "docs" / "testing" / "known-issues-archive.md"
SCRATCH = REPO_ROOT / "logs" / "ops66-negative-control"

# The tool's last commit before OPS-66 — never HEAD: (brief, "Traps").
TOOL_PINNED_SHA = "e5af2cae111b35c95a5b17022ae299b7e5d1139e"
# known-issues.md's last commit before OPS-66 (fallback once the four are rotated).
KI_PINNED_SHA = "32fabbec07b000f03074a161e1096379e9ddcdb4"

OLD_ENTRY_RE = re.compile(r"^(### |## \d{4}-\d{2}-\d{2})")

FOUR = (
    "## RETIRED 2026-09-20 (`OPS-59` (b), 07:30 slot) — `tests/unit/test_setup_figure_title.py::` "
    "the call-site count pin (`EXPECTED_EXAMPLE_CALL_SITES = 18`) was stale on `main`",
    "## ✅ RETIRED 2026-09-21 (`OPS-59` (a), 06:00 implementer slot) — ~~2026-09-19 — "
    "`tests/unit/test_doc_reference_exit_codes.py::test_the_in_tree_exemption_cannot_silently_widen` "
    "is red on `main`: the `OPS-44` artifact pin was never widened for the `EX-57` setup figures~~",
    "## ✅ RETIRED 2026-09-20 (`PORT-20` steps 1–3, 04:30 + 09:00 implementer slots) — ~~🟡 OPEN "
    "2026-09-19 (`ANS-3` adjudication, weekly review) — on the gap-voltage / impressed-current route "
    "`run_n_port_sparameter_sweep` reports an `S` that is **not the 50 Ω S-matrix**; its `z_matrix` "
    "is right~~",
    "## ✅ RETIRED 2026-09-27 (03:00 daily review; fixed by `4719e82`, 2026-09-26 09:00 slot) — was: "
    "2026-09-20 — 🟡 OPEN (`PORT-20` step 3, 09:00 implementer slot) — `ans:3`'s tracked "
    "`metrics.json` / `COMPARISON.md` still carry the **pre-`PORT-20`** S table; the example cannot "
    "be re-run in a slot without rewriting the gitignored private comparison",
)

SECTION_HEADERS = (
    "## How to check a failure against this baseline",
    "## Failing tests",
    "## Non-test issues",
    "## Recording a new entry",
)


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


tool = _load("rotate_plan_archive_ops66", TOOL)


def _git(*args: str) -> subprocess.CompletedProcess:
    return subprocess.run(
        ["git", "-c", f"safe.directory={REPO_ROOT}", "-C", str(REPO_ROOT), *args],
        capture_output=True, text=True, check=False,
    )


def _git_show(sha: str, rel: str) -> str:
    proc = _git("show", f"{sha}:{rel}")
    assert proc.returncode == 0, f"git show {sha}:{rel} failed: {proc.stderr.strip()}"
    return proc.stdout


def _headings(text: str) -> list[str]:
    """Every heading line outside a fenced block, the tool's own scan."""
    out, in_fence = [], False
    for l in text.split("\n"):
        if l.lstrip().startswith("```"):
            in_fence = not in_fence
        elif not in_fence and re.match(r"^#{1,3} ", l):
            out.append(l)
    return out


def _live_text() -> tuple[str, str]:
    text = KI.read_text(encoding="utf-8")
    lines = set(text.split("\n"))
    if all(h in lines for h in FOUR):
        return text, "working known-issues.md"
    return _git_show(KI_PINNED_SHA, "docs/testing/known-issues.md"), f"known-issues.md@{KI_PINNED_SHA[:8]}"


# --- A1 ---------------------------------------------------------------------

@pytest.mark.parametrize("heading", FOUR, ids=["OPS-59b", "OPS-59a", "PORT-20", "PORT-20-ans3"])
def test_a1_the_four_headings_are_entries_and_closed(heading):
    assert tool.KI_ENTRY_RE.match(heading), heading[:60]
    assert tool._ki_closed(heading), heading[:60]
    assert not OLD_ENTRY_RE.match(heading), "the old pattern already matched — not the defect"


# --- A2 ---------------------------------------------------------------------

def _new_minus_old(text: str) -> tuple[list[str], set[str], set[str], set[str]]:
    heads = _headings(text)
    old = {h for h in heads if OLD_ENTRY_RE.match(h)}
    new = {h for h in heads if tool.KI_ENTRY_RE.match(h)}
    return heads, old, new, new - old


def test_a2_old_subset_of_new_and_difference_is_exactly_the_four():
    live, src = _live_text()
    archive = KI_ARCHIVE.read_text(encoding="utf-8")
    diffs = {}
    for label, text in (("live", live), ("archive", archive)):
        heads, old, new, extra = _new_minus_old(text)
        assert old <= new, f"{label}: old ⊄ new: {sorted(old - new)[:3]}"
        unclosed = [h for h in extra if not tool._ki_closed(h)]
        assert not unclosed, f"{label}: new − old carries non-closed headings: {unclosed[:3]}"
        diffs[label] = extra
        print(f"[OPS-66 A2] {label} ({src if label == 'live' else 'known-issues-archive.md'}): "
              f"headings={len(heads)} old={len(old)} new={len(new)} new-old={len(extra)}")
    assert diffs["live"] == set(FOUR), (
        f"live new − old is {len(diffs['live'])}, expected exactly the four: "
        f"{[h[:60] for h in diffs['live']]}")
    # OPS-66 step 3 (§9 item 393): the former clause `archive ⊆ FOUR` asserted a
    # snapshot of the archive, which later rotations broke by moving further
    # closed undated headings there (known-issues, 2026-10-12). The docstring's
    # A2 is "new − old ⊆ closed" (asserted above for both files) and "the four
    # are found"; the control shows the snapshot clause still fails on today's
    # archive, so the replacement is what turns this green.
    beyond = diffs["archive"] - set(FOUR)
    print(f"[OPS-66 A2] archive new-old beyond the four: {len(beyond)}")
    for h in sorted(beyond):
        print(f"[OPS-66 A2]   {h[:60]}")
    assert beyond, "control: the old snapshot clause would pass on today's archive"
    assert set(FOUR) <= diffs["archive"] | diffs["live"], (
        [h[:60] for h in set(FOUR) - diffs["archive"] - diffs["live"]])
    # Second control: a planted admitted-but-unclosed heading must be caught by
    # the closed check. The brief's form first; the first admitted-and-unclosed
    # form is used.
    planted = None
    for form in ("## RETIRED-NOT 2026-10-12 planted", "## RETIRED 2026-10-12 PARTIAL planted"):
        admitted, closed = bool(tool.KI_ENTRY_RE.match(form)), tool._ki_closed(form)
        print(f"[OPS-66 A2] planted {form!r}: admitted={admitted} closed={closed}")
        if admitted and not closed and planted is None:
            planted = form
    assert planted is not None, "no planted form is admitted-but-unclosed"
    _, _, _, extra = _new_minus_old(archive + "\n" + planted + "\n")
    assert [h for h in extra if not tool._ki_closed(h)] == [planted]


# --- A3 ---------------------------------------------------------------------

def test_a3_section_headers_are_not_entries():
    for h in SECTION_HEADERS + (tool.KI_INDEX_HEADING,):
        assert not tool.KI_ENTRY_RE.match(h), f"section header became an entry: {h}"
    # And over the live file: every `## ` heading the new pattern admits is
    # either dated or one of the closed undated-retired forms.
    live, _ = _live_text()
    for h in _headings(live):
        if h.startswith("## ") and tool.KI_ENTRY_RE.match(h):
            assert re.match(r"^## \d{4}-\d{2}-\d{2}", h) or tool._ki_closed(h), h[:80]


def test_a3_struck_open_clauses_stay_closed():
    struck = [h for h in FOUR if "🟡 OPEN" in h]
    assert len(struck) == 2, f"expected two headings with a struck 🟡 OPEN clause, got {len(struck)}"
    for h in struck:
        assert not h.startswith("## 🟡") and tool._ki_closed(h), h[:60]


# --- NC: pinned pre-change tool vs new tool ---------------------------------

def _dry_run_count(module, ki_path: Path, capsys, monkeypatch) -> tuple[int, list[str]]:
    monkeypatch.setattr(module, "KNOWN_ISSUES", ki_path)
    capsys.readouterr()
    module.rotate_known_issues("2026-10-01", dry_run=True)
    out = capsys.readouterr().out
    m = re.search(r"^(\d+) RETIRED entries", out, re.M)
    n = int(m.group(1)) if m else 0
    if not m:
        assert "no RETIRED entries" in out, out
    listed = [l.strip() for l in out.split("\n") if l.startswith("  - ")]
    return n, listed


def test_negative_control_pinned_tool_reports_none_new_reports_four(capsys, monkeypatch):
    SCRATCH.mkdir(parents=True, exist_ok=True)
    old_path = SCRATCH / "rotate_plan_archive_e5af2cae.py"
    old_path.write_text(_git_show(TOOL_PINNED_SHA, "scripts/maintenance/rotate_plan_archive.py"),
                        encoding="utf-8")
    old_tool = _load("rotate_plan_archive_pinned_e5af2cae", old_path)
    assert not hasattr(old_tool, "KI_ENTRY_RE"), "pinned tool already carries KI_ENTRY_RE"
    live, src = _live_text()
    ki_copy = SCRATCH / "known-issues.copy.md"
    ki_copy.write_text(live, encoding="utf-8")
    before = ki_copy.read_bytes()

    n_old, listed_old = _dry_run_count(old_tool, ki_copy, capsys, monkeypatch)
    n_new, listed_new = _dry_run_count(tool, ki_copy, capsys, monkeypatch)
    with capsys.disabled():
        print(f"\n[OPS-66 NC] source={src} pinned={TOOL_PINNED_SHA[:8]} old={n_old} new={n_new}")
        for l in listed_new:
            print(f"[OPS-66 NC]   {l}")
    assert ki_copy.read_bytes() == before, "a dry run wrote the scratch copy"
    assert not any(any(h[3:60] in l for h in FOUR) for l in listed_old)
    assert (n_old, n_new) == (0, 4)
    for h in FOUR:
        assert any(re.sub(r"^#+\s*", "", h)[:100] in l for l in listed_new), h[:60]
