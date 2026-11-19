"""WF-15 step 9 (§9 item 662, U95): ``ColeColeMaterial`` joins the public surface.

(Δ1)  ``public_surface.txt`` differs from its parent commit's copy by exactly
      one added name, ``ColeColeMaterial``, and nothing removed. The parent is
      pinned by sha (``PARENT_SHA``), never ``HEAD:`` — a control that reads
      ``HEAD:`` of the file it changes goes stale on landing (``OPS-43``).
      The parent file is read with ``git show`` when git and the object are
      available; otherwise the test FAILS (never skips).
(Imp) ``from fem_em_solver import ColeColeMaterial`` is the same object as
      ``fem_em_solver.materials.ColeColeMaterial``, and the name is in
      ``fem_em_solver.__all__`` and in the snapshot.
(c)   Negative control: the diff helper returns two added names on a synthetic
      list carrying a second, planted extra name, so (Δ1) can fail.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

import fem_em_solver

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
SNAPSHOT_PATH = HERE / "public_surface.txt"
SNAPSHOT_REL = "tests/package/public_surface.txt"
# The commit item 662 lands on (main before U95), pinned — never HEAD:.
PARENT_SHA = "44159a26"
ADDED_NAME = "ColeColeMaterial"


def _names(text: str) -> list[str]:
    return [ln.strip() for ln in text.splitlines() if ln.strip()]


def surface_diff(parent: list[str], child: list[str]) -> tuple[list[str], list[str]]:
    """Return (added, removed) names, each sorted."""
    p, c = set(parent), set(child)
    return sorted(c - p), sorted(p - c)


def _parent_snapshot() -> list[str]:
    out = subprocess.run(
        # safe.directory: the container user does not own the bind mount
        # (first W1, 2026-10-27: "detected dubious ownership").
        [
            "git",
            "-c",
            f"safe.directory={REPO}",
            "-C",
            str(REPO),
            "show",
            f"{PARENT_SHA}:{SNAPSHOT_REL}",
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert (
        out.returncode == 0
    ), f"git show {PARENT_SHA}:{SNAPSHOT_REL} failed ({out.returncode}): {out.stderr.strip()}"
    return _names(out.stdout)


def test_wf15_s9_surface_differs_from_parent_by_one_name():
    parent = _parent_snapshot()
    child = _names(SNAPSHOT_PATH.read_text())
    added, removed = surface_diff(parent, child)
    print(
        f"\n[WF-15 s9] parent {PARENT_SHA}: {len(parent)} names; now {len(child)}; "
        f"added {added}; removed {removed}",
        flush=True,
    )
    assert added == [ADDED_NAME]
    assert removed == []
    assert len(child) == len(parent) + 1
    assert child == sorted(child), "snapshot is not sorted"


def test_wf15_s9_import_identity_and_all():
    from fem_em_solver import ColeColeMaterial as top
    from fem_em_solver.materials import ColeColeMaterial as sub

    print(
        f"\n[WF-15 s9] top is sub: {top is sub}; len(__all__) {len(fem_em_solver.__all__)}",
        flush=True,
    )
    assert top is sub
    assert ADDED_NAME in fem_em_solver.__all__
    assert sorted(fem_em_solver.__all__) == _names(SNAPSHOT_PATH.read_text())


def test_wf15_s9_negative_control_planted_extra_name():
    parent = ["CapacitorSheets", "Study"]
    planted = sorted(parent + [ADDED_NAME, "PlantedExtraName"])
    added, removed = surface_diff(parent, planted)
    print(f"\n[WF-15 s9] control added {added}; removed {removed}", flush=True)
    assert added == [ADDED_NAME, "PlantedExtraName"]
    assert len(added) == 2
    assert added != [ADDED_NAME]
    assert removed == []
