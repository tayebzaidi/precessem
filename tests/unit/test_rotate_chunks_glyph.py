"""OPS-66 step 2 (§9 item 384): ``rotate_plan_archive.py chunks`` no longer
doubles a state line's leading glyph.

The replacement cell is ``<row glyph> <state line>``; a state line written from
the row's last ruling ("✅ 2026-10-10 …") made it "✅ ✅ 2026-10-10 …" (the
2026-10-10 weekly's twelve rows, de-duplicated by hand in 397b69d0).

Anchors (asserted):
  G1  glyph-led state → exactly one glyph in the new status cell;
  G2  no-glyph state → plan bytes, chunk-file bytes, --dry-run and --census
      stdout identical to the pinned pre-change tool;
  G3  a state opening with a different glyph is refused, nothing written;
  G4  the variation-selector glyph ⚠️ (U+26A0 U+FE0F) is stripped whole;
  NC  (negative control, binary) the pinned pre-change tool doubles on the G1
      fixture (2 glyphs), the new tool does not (1).
Printed only: the real-data replay of the weekly's rotation on nine rows.
Every fixture is synthetic under tmp_path; the live plan is never touched.
"""
from __future__ import annotations

import importlib.util
import re
import subprocess
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
TOOL = REPO_ROOT / "scripts" / "maintenance" / "rotate_plan_archive.py"
SCRATCH = REPO_ROOT / "logs" / "ops66-step2-control"

# The tool's last commit before the fix — never HEAD: (brief, "Traps").
TOOL_PINNED_SHA = "15bdce61ca16199c6b1f3e947b6f13ad020e6ed8"
PRE_ROTATION_SHA = "38aa4afd"   # 397b69d0^, the plan before the weekly's rotation
ROTATED_SHA = "397b69d0"        # the weekly's rotation, glyphs de-duplicated by hand
REPLAY_IDS = ("OPS-70", "GEO-36", "MAT-15", "WF-1", "WF-10", "WF-12", "WF-15", "WF-16", "WF-17")


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _git(*args: str) -> subprocess.CompletedProcess:
    return subprocess.run(
        ["git", "-c", f"safe.directory={REPO_ROOT}", "-C", str(REPO_ROOT), *args],
        capture_output=True, text=True, check=False,
    )


def _git_show(sha: str, rel: str) -> str:
    proc = _git("show", f"{sha}:{rel}")
    assert proc.returncode == 0, f"git show {sha}:{rel} failed: {proc.stderr.strip()}"
    return proc.stdout


new = _load("rotate_plan_archive_ops66s2_new", TOOL)


@pytest.fixture(scope="module")
def old():
    SCRATCH.mkdir(parents=True, exist_ok=True)
    path = SCRATCH / "rotate_plan_archive_pinned.py"
    path.write_text(_git_show(TOOL_PINNED_SHA, "scripts/maintenance/rotate_plan_archive.py"),
                    encoding="utf-8")
    return _load("rotate_plan_archive_ops66s2_old", path)


PLAN = """# Synthetic plan

## 7. Backlog

| ID | Title | Status | Tier |
|---|---|---|---|
| `ZZ-1` | Synthetic one — long history 2026 10 11 | ✅ 2026-10-10 closed after a long history. More text. | smoke |
| `ZZ-2` | Synthetic two — history 2026 | 🟡 open 2026 partial. Done when the gate passes. | smoke |
| `ZZ-3` | Synthetic three — history 2026 10 | ⚠️ 2026-10-10 placeholder backed. | smoke |

## 8. Next

tail
"""


def _fixture(tmp_path: Path, tag: str, spec_text: str) -> tuple[Path, Path, Path]:
    d = tmp_path / tag
    d.mkdir()
    plan = d / "PROJECT_PLAN.md"
    plan.write_text(PLAN, encoding="utf-8")
    spec = d / "spec.txt"
    spec.write_text(spec_text, encoding="utf-8")
    return spec, plan, d / "chunks"


def _row(plan: Path, cid: str) -> str:
    return next(l for l in plan.read_text(encoding="utf-8").split("\n") if l.startswith(f"| `{cid}` |"))


def _status_cell(row: str) -> str:
    """Text from the title/status boundary to ``*History:``."""
    head = row.split("*History:", 1)[0]
    pipes = new._split_unescaped(head)
    return head[pipes[2] + 1:]


def _count(row: str, glyph: str) -> int:
    return _status_cell(row).count(glyph)


def test_g1_glyph_led_state_and_negative_control(tmp_path, old):
    spec_text = "=== ZZ-1\n✅ 2026-10-10 state.\n"
    res = {}
    for tag, mod in (("new", new), ("old", old)):
        spec, plan, chunks = _fixture(tmp_path, tag, spec_text)
        rc = mod.rotate_chunks(spec, plan, chunks, dry_run=False, date="2026-10-11")
        assert rc == 0, tag
        row = _row(plan, "ZZ-1")
        res[tag] = (_count(row, "✅"), _status_cell(row))
        print(f"G1/NC {tag}: glyphs={res[tag][0]} cell={res[tag][1]!r}")
    assert res["new"][1].startswith(" ✅ 2026-10-10 state.")
    # NC, the binary: old doubles, new does not.
    assert (res["old"][0], res["new"][0]) == (2, 1)


@pytest.mark.parametrize("with_dw", [False, True])
def test_g2_no_glyph_state_byte_identical(tmp_path, old, capsys, with_dw):
    spec_text = "=== ZZ-2\nstate without a glyph 2026.\n"
    if with_dw:
        spec_text += "Done-when: Done when the gate passes.\n"
    out = {}
    for tag, mod in (("new", new), ("old", old)):
        spec, plan, chunks = _fixture(tmp_path, tag, spec_text)
        capsys.readouterr()
        assert mod.rotate_chunks(spec, plan, chunks, dry_run=True, date="2026-10-11") == 0
        dry = capsys.readouterr().out
        assert mod.census_chunks(plan, 0) == 0
        census = capsys.readouterr().out
        assert mod.rotate_chunks(spec, plan, chunks, dry_run=False, date="2026-10-11") == 0
        capsys.readouterr()
        out[tag] = (plan.read_bytes(), (chunks / "ZZ-2.md").read_bytes(), dry, census)
    for k, name in enumerate(("plan", "chunk", "dry-run stdout", "census stdout")):
        a, b = out["new"][k], out["old"][k]
        if a != b:
            first = next((i for i, (x, y) in enumerate(zip(a, b)) if x != y), min(len(a), len(b)))
            pytest.fail(f"G2 {name} differs at byte {first}")
    with capsys.disabled():
        print(f"\nG2 (with_dw={with_dw}): plan/chunk/dry-run/census byte-identical "
              f"({len(out['new'][0])} / {len(out['new'][1])} / {len(out['new'][2])} / {len(out['new'][3])} B)")


def test_g3_mismatched_glyph_refused(tmp_path, capsys):
    spec, plan, chunks = _fixture(tmp_path, "g3", "=== ZZ-1\n🟡 2026-10-10 state.\n")
    before = plan.read_bytes()
    rc = new.rotate_chunks(spec, plan, chunks, dry_run=False, date="2026-10-11")
    printed = capsys.readouterr().out
    with capsys.disabled():
        print(f"\nG3: rc={rc} {printed.strip()!r}")
    assert rc == 1
    assert "REFUSED: ZZ-1: state line opens with 🟡" in printed
    assert plan.read_bytes() == before
    assert not (chunks / "ZZ-1.md").exists()


def test_g4_variation_selector_glyph(tmp_path, old):
    counts = {}
    for tag, mod in (("new", new), ("old", old)):
        spec, plan, chunks = _fixture(tmp_path, tag, "=== ZZ-3\n⚠️ 2026-10-10 x.\n")
        assert mod.rotate_chunks(spec, plan, chunks, dry_run=False, date="2026-10-11") == 0
        row = _row(plan, "ZZ-3")
        counts[tag] = _count(row, "⚠️")
        if tag == "new":
            assert _status_cell(row).startswith(" ⚠️ 2026-10-10 x.")
            assert "️️" not in row and _status_cell(row).count("⚠") == 1
    print(f"G4: old={counts['old']} new={counts['new']}")
    assert counts == {"old": 2, "new": 1}


def test_replay_real_rotation_printed_only(tmp_path, old):
    """Printed only: the weekly's rotation re-run on nine rows; no assertion."""
    try:
        pre = _git_show(PRE_ROTATION_SHA, "PROJECT_PLAN.md")
        post = _git_show(ROTATED_SHA, "PROJECT_PLAN.md").split("\n")
    except AssertionError as e:
        print(f"REPLAY skipped: {e}")
        return
    equal = doubled = 0
    for cid in REPLAY_IDS:
        try:
            target = next(l for l in post if l.startswith(f"| `{cid}` |"))
            i = post.index(target)
            p = new._parse_row(target, cid, new._ncols(post, i))
            span, glyph = p["span"], p["glyph"]
            m = re.search(r" \| (?:\*\*)?(?=" + re.escape(glyph) + ")", span)
            status = span[m.end():]
            state = re.split(r" \*\*Done-when:\*\* | \*History:", status, 1)[0].strip()
            dwm = re.search(r"\*\*Done-when:\*\* (.*?) \*History:", status)
            spec_text = f"=== {cid}\n{state}\n" + (f"Done-when: {dwm.group(1)}\n" if dwm else "")
            res = {}
            for tag, mod in (("new", new), ("old", old)):
                d = tmp_path / f"replay-{cid}-{tag}"
                d.mkdir()
                plan = d / "PROJECT_PLAN.md"
                plan.write_text(pre, encoding="utf-8")
                spec = d / "spec.txt"
                spec.write_text(spec_text, encoding="utf-8")
                rc = mod.rotate_chunks(spec, plan, d / "chunks", dry_run=False, date="2026-10-11")
                res[tag] = (rc, _row(plan, cid))
            eq = res["new"][0] == 0 and res["new"][1] == target
            dbl = res["old"][0] == 0 and _count(res["old"][1], glyph) == 2
            equal += eq
            doubled += dbl
            print(f"REPLAY {cid:7} new==397b69d0: {eq}  old doubled: {dbl}  (rc new {res['new'][0]}, old {res['old'][0]})")
        except Exception as e:  # printed only, never a red
            print(f"REPLAY {cid:7} not reconstructed: {type(e).__name__}: {e}")
    print(f"REPLAY summary: equal {equal}/{len(REPLAY_IDS)}, old doubled {doubled}/{len(REPLAY_IDS)}")
