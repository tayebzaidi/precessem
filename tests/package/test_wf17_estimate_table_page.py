"""`WF-17` step 30 (§9 item 827, §10 U130): the estimate table as a generated page.

Brief: docs/planning/chunks/WF-17.md, "§9 item 827".  Generator:
``scripts/docs/generate_estimate_page.py``; page: ``docs/guide/estimates.md``.
Run under ``mpiexec -n 2`` (the page states its launch width).

W0 (no mesh):
  discovery -- globbing ``examples/`` and ``docs/guide/`` for study mappings (a
  ``femem:`` key; gitignored ``paraview_output/`` copies excluded) yields exactly
  the pinned 8 files + 6 fences; the tier rule on three synthetic estimates;
  (c-ref) the thermal-box file through the same row builder returns rc 2 with
  ``_refuse_thermal_only``'s message as its refusal cell, CLI stderr = API
  message; (c-stale) the tracked page with one cell mutated differs from the
  tracked text and the diff names the file; (c-width) ``_predict`` seconds at
  1 and 2 ranks have ratio exactly 2.0 on one fitted row.

W1 (knob ``FEM_EM_WF17_S30=1``, complex build, -n 2, five real meshes):
  (CLI=API) for each of the 14 entries every printed field equals
  ``repr(getattr(est, f))`` (``_parse`` / ``_diff`` of ``test_wf13_cli_estimate``),
  the field set complete, rc 0; (Page) fresh render == tracked page byte for
  byte; (B1) real builds per rank == distinct geometries (5) and memo hits ==
  2 x 14 - 5.
"""

from __future__ import annotations

import dataclasses
import importlib.util
import os
import time
from pathlib import Path

import pytest
from mpi4py import MPI

from tests.package.test_wf13_cli_estimate import _diff, _parse

ROOT = Path(__file__).resolve().parent.parent.parent
GENERATOR = ROOT / "scripts" / "docs" / "generate_estimate_page.py"
PAGE = ROOT / "docs" / "guide" / "estimates.md"
KNOB = "FEM_EM_WF17_S30"
THERMAL_BOX = ROOT / "tests" / "package" / "wf1" / "thermal_box_steady_m4.yaml"
THERMAL_BOX_MESSAGE = (
    "estimate is not defined on a thermal-only Study "
    "(ThermalBoxGeometry has no EM solve; use temperature())"
)
N_ENTRIES = 14
N_GEOMETRIES = 5  # sphere, F-small four-port, F-small hole, rod, F-human hole (brief)


def _load_generator():
    spec = importlib.util.spec_from_file_location("wf17_generate_estimate_page", GENERATOR)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


gen = _load_generator()


def _p0(msg):
    if MPI.COMM_WORLD.rank == 0:
        print(msg, flush=True)


# ---------------------------------------------------------------- W0
def test_wf17_s30_discovery_equals_pinned_list():
    files, fences = gen.discover()
    _p0(f"\n[WF-17 s30] discovered files {files}\n[WF-17 s30] discovered fences {fences}")
    assert files == sorted(gen.FILES)
    assert fences == sorted(gen.FENCES)
    assert len(gen.FILES) + len(gen.FENCES) == N_ENTRIES
    assert not any("clinical" in f for f in files + fences)


def test_wf17_s30_tier_rule_synthetic():
    from fem_em_solver.study import ResourceEstimate

    def est(ranks, seconds):
        return ResourceEstimate(
            cells=1,
            unknowns=1,
            degree=2,
            ranks=ranks,
            predicted_peak_gib=1.0,
            predicted_seconds=seconds,
            rows_used=(),
        )

    got = [
        gen.tier_for(est(2, 30.0)),
        gen.tier_for(est(2, 180.5)),
        gen.tier_for(est(16, 10.0)),
        gen.tier_for(est(12, 1200.0)),
        gen.tier_for(est(2, 1200.1)),
    ]
    _p0(f"\n[WF-17 s30] tiers (2,30) (2,180.5) (16,10) (12,1200) (2,1200.1): {got}")
    assert got == ["smoke", "heavy", "xl", "heavy", "xl"]


def test_wf17_s30_c_ref_thermal_box_refusal():
    comm = MPI.COMM_WORLD
    row = gen.row_for("c-ref", THERMAL_BOX, comm)
    _p0(
        f"\n[WF-17 s30] (c-ref) rc={row['rc']} api_error={row['api_error']!r} "
        f"stderr={row['stderr']!r}"
    )
    assert row["rc"] == 2
    assert row["est"] is None
    assert row["api_error"] == THERMAL_BOX_MESSAGE
    if comm.rank == 0:
        assert row["stderr"].strip() == f"femem: error: {THERMAL_BOX_MESSAGE}"
    page = gen.render([row])
    cell = f"rc 2: {THERMAL_BOX_MESSAGE}"
    assert cell in page


def test_wf17_s30_c_stale_mutated_cell_names_file(tmp_path):
    text = PAGE.read_text(encoding="utf-8")
    row = next(
        ln for ln in text.splitlines() if ln.startswith("| `examples/mri/09_study_from_yaml_sphere")
    )
    mutated = text.replace(row, row.replace("`74020`", "`74021`", 1), 1)
    assert mutated != text
    planted = tmp_path / "estimates.md"
    planted.write_text(mutated, encoding="utf-8")
    diff = gen.diff_against(planted, text)
    _p0("\n[WF-17 s30] (c-stale) diff:\n" + "\n".join(diff))
    assert planted.read_bytes() != text.encode("utf-8")
    assert diff and str(planted) in diff[0]
    assert sum(ln.startswith("-|") for ln in diff) == 1
    assert gen.diff_against(PAGE, text) == []


def test_wf17_s30_c_width_seconds_ratio():
    from fem_em_solver.study import PRICE_TABLE, _predict

    row = next(r for r in PRICE_TABLE if r["unknowns"] and r["seconds"] and r["degree"] == 2)
    s1 = _predict(row["unknowns"], 2, 1)[1]
    s2 = _predict(row["unknowns"], 2, 2)[1]
    _p0(f"\n[WF-17 s30] (c-width) row {row['row']}: s(1)={s1!r} s(2)={s2!r} ratio={s1 / s2!r}")
    assert s1 / s2 == 2.0


# ---------------------------------------------------------------- W1
def test_wf17_s30_cli_equals_api_and_page_fresh(monkeypatch):
    if os.environ.get(KNOB) != "1":
        pytest.skip(f"{KNOB} unset: the five-mesh estimate window is environment-selected")
    import fem_em_solver.study_yaml as study_yaml
    from fem_em_solver.study import ResourceEstimate, Study

    def _boom(*a, **k):
        raise AssertionError("solve reached")

    monkeypatch.setattr(Study, "solve", _boom)
    monkeypatch.setattr(study_yaml, "run_study_file", _boom)
    comm = MPI.COMM_WORLD
    assert comm.size == gen.LAUNCH_RANKS, comm.size
    t0 = time.perf_counter()
    text, rows, count = gen.build(comm)
    elapsed = time.perf_counter() - t0
    names = [f.name for f in dataclasses.fields(ResourceEstimate)]
    bad = None
    if comm.rank == 0:
        bad = {}
        for r in rows:
            parsed = _parse(r["stdout"], names)
            complete = set(parsed) == set(names) and all(len(v) == 1 for v in parsed.values())
            differ = _diff(parsed, r["est"]) if (complete and r["est"] is not None) else None
            print(
                f"[WF-17 s30] (CLI=API) {r['label']}: rc={r['rc']} complete={complete} "
                f"differ={differ}",
                flush=True,
            )
            if r["rc"] != 0 or not complete or differ != set() or r["stderr"]:
                bad[r["label"]] = (r["rc"], complete, differ, r["stderr"])
    bad = comm.bcast(bad, root=0)
    diff = gen.diff_against(PAGE, text)
    _p0(
        f"[WF-17 s30] (B1) builds={count['builds']} hits={count['hits']} "
        f"geometries={count['keys']}; (Page) diff lines {len(diff)}; elapsed {elapsed:.1f} s"
    )
    assert len(rows) == N_ENTRIES
    assert bad == {}
    assert diff == []
    assert PAGE.read_bytes() == text.encode("utf-8")
    assert count["builds"] == N_GEOMETRIES
    assert count["hits"] == 2 * N_ENTRIES - N_GEOMETRIES
