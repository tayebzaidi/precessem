"""OPS-65: retention's citing set covers the chunk history files and the XL ledger.

`housekeeping.cited_basenames()` decides which logs are gating (never deleted).
Before OPS-65 it read PROJECT_PLAN.md, plan-archive.md, known-issues.md and
three globs only, so logs cited solely from `docs/planning/chunks/*.md` or
`docs/testing/xl-ledger.md` were swept (214 deleted by 2026-09-27).

Negative control: set HOUSEKEEPING_PATH to a copy of housekeeping.py as of the
pre-change commit; both citation tests must fail.
"""
from __future__ import annotations

import importlib.util
import os
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
DEFAULT = ROOT / "scripts" / "maintenance" / "housekeeping.py"
CITATION = re.compile(r"\d{8}T\d{6}Z_[A-Za-z0-9._-]+\.log")

# The citing locations as they stood before OPS-65 (housekeeping.py:42-47 at e1964dff).
LEGACY_FILES = [
    ROOT / "PROJECT_PLAN.md",
    ROOT / "docs" / "planning" / "plan-archive.md",
    ROOT / "docs" / "testing" / "known-issues.md",
]
LEGACY_GLOBS = ["examples/**/*.md", "docs/validation/**/*.md", "docs/ports/**/*.md"]
KEEP_FILE = ROOT / "docs" / "testing" / "retention-keep.txt"


def _load():
    path = Path(os.environ.get("HOUSEKEEPING_PATH", DEFAULT))
    spec = importlib.util.spec_from_file_location("housekeeping_under_test", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    print(f"housekeeping under test: {path}")
    return mod


def _cites(path: Path) -> set[str]:
    return set(CITATION.findall(path.read_text(encoding="utf-8", errors="replace")))


def _legacy_cited() -> set[str]:
    files = list(LEGACY_FILES)
    for g in LEGACY_GLOBS:
        files.extend(ROOT.glob(g))
    out: set[str] = set()
    for f in files:
        if f.is_file():
            out |= _cites(f)
    if KEEP_FILE.is_file():
        for line in KEEP_FILE.read_text().splitlines():
            line = line.split("#", 1)[0].strip()
            if line:
                out.add(line.removesuffix(".gz"))
    return out


@pytest.fixture(scope="module")
def cited():
    return _load().cited_basenames()


@pytest.mark.parametrize(
    "source",
    ["docs/planning/chunks/WF-6.md", "docs/testing/xl-ledger.md"],
)
def test_basename_cited_only_from_new_location_is_gating(cited, source):
    only_here = sorted(_cites(ROOT / source) - _legacy_cited())
    print(f"{source}: {len(only_here)} basenames cited only here (vs legacy set)")
    assert only_here, f"fixture premise: {source} must cite a basename no legacy source cites"
    print(f"  example: {only_here[0]}")
    missing = [b for b in only_here if b not in cited]
    print(f"  absent from cited_basenames(): {len(missing)}")
    assert not missing, f"{len(missing)} of {len(only_here)} not gating, e.g. {missing[:3]}"
