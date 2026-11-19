"""`WF-13` step 12 (§9 item 663, chain U96): F-clinical-1.5T through ``Study.estimate`` only.

Brief: docs/planning/chunks/WF-13.md, "§9 item 663".  No mesh is built anywhere
in this module: every ``MeshGenerator`` entry and ``gmsh.initialize`` are patched
to raise in every node that reaches the package.

Asserted:
  (FC1) ``Study.estimate()`` on the ``f_clinical_birdcage`` kind returns cells
        1 147 366 and unknowns 1 525 193, equal to `GEO-38`'s test constants
        (``tests/validation/test_geo38_f_clinical_build.py:86``, ``:89``);
        ``rows_used == ("FC1",)``; seconds / GiB / ranks equal to the measured
        row (entry 46's Z50 pair 203.3 s, memory.peak 49.07 GiB, -n 16).
  (NoB) no builder is called (the patches above raise; ``estimate`` returns).
  (Ref) ``solve()`` raises ``ValueError`` containing "xl tier only" before any
        builder; ``femem run`` without ``--estimate`` returns 2 on every rank and
        writes no output directory; with ``--estimate`` it returns 0.
  (CLI) the CLI's printed estimate equals the API's field by field (repr).
  (Iso) ``PRICE_TABLE`` carries no clinical row (the per-degree fits are untouched).
Negative controls (asserted): degree 3, which has no measured row (degree 2 has
FC2 since step 13, §9 item 715), returns ``not_priced`` with NaN (never a
fitted number, although ``_predict`` would fit one); a planted second clinical
row with a different ``peak_gib`` changes the estimate (the row is what is read);
a non-default geometry is refused.
"""

from __future__ import annotations

import dataclasses
import math
from pathlib import Path

import pytest
from mpi4py import MPI

FC = Path(__file__).resolve().parent / "wf13" / "f_clinical_15t.yaml"


def _say(msg):
    if MPI.COMM_WORLD.rank == 0:
        print(msg, flush=True)


@pytest.fixture
def no_build(monkeypatch):
    """Every MeshGenerator entry and gmsh.initialize raise (NoB)."""
    import gmsh

    from fem_em_solver.io.mesh import MeshGenerator

    calls = []

    def _boom(*a, **k):
        calls.append(1)
        raise AssertionError("a mesh builder was called")

    patched = []
    for name, val in list(vars(MeshGenerator).items()):
        if name.startswith("__"):
            continue
        if isinstance(val, (staticmethod, classmethod)) or callable(val):
            monkeypatch.setattr(MeshGenerator, name, staticmethod(_boom))
            patched.append(name)
    monkeypatch.setattr(gmsh, "initialize", _boom)
    assert "birdcage_port_domain" in patched
    return calls


def _study(degree=1):
    from fem_em_solver.study_yaml import load_study

    study = load_study(FC)
    if degree != 1:
        study.degree = degree
    return study


def test_wf13_s12_fc1_estimate_reads_records(no_build):
    from fem_em_solver.study import (
        F_CLINICAL_R1_CELL_RECORD,
        F_CLINICAL_R1_UNKNOWNS_RECORD,
        PRICE_TABLE,
        PRICE_TABLE_CLINICAL,
        FClinicalBirdcageGeometry,
    )
    from tests.validation import test_geo38_f_clinical_build as geo38

    assert F_CLINICAL_R1_CELL_RECORD == geo38.F_CLINICAL_R1_CELL_RECORD == 1_147_366
    assert F_CLINICAL_R1_UNKNOWNS_RECORD == geo38.F_CLINICAL_R1_UNKNOWNS_RECORD
    study = _study()
    assert isinstance(study.geometry, FClinicalBirdcageGeometry)
    g = study.geometry
    assert (g.leg_count, g.ring_radius, g.shield_radius) == (
        geo38.LEG_COUNT,
        geo38.F_CLINICAL_15T_KWARGS["ring_radius"],
        geo38.F_CLINICAL_15T_KWARGS["shield_radius"],
    )
    assert g.resolution == geo38.F_CLINICAL_15T_KWARGS["resolution"]
    assert g.conductor_resolution == geo38.CONDUCTOR_RESOLUTION["r1"]

    est = study.estimate(MPI.COMM_WORLD)
    row, _fc2 = PRICE_TABLE_CLINICAL  # FC2: `WF-13` step 13 (§9 item 715)
    _say(f"\n[WF-13 s12] FC1: {est}")
    assert est.cells == 1_147_366
    assert est.unknowns == 1_525_193 == geo38.F_CLINICAL_R1_UNKNOWNS_RECORD[0]
    assert est.degree == 1
    assert est.rows_used == ("FC1",)
    assert est.ranks == row["ranks"] == 16
    assert est.predicted_seconds == row["seconds"] == 203.3
    assert est.predicted_peak_gib == row["peak_gib"] == 49.07
    assert row["build_seconds"] == 357.7 and row["cells"] == est.cells
    assert row["unknowns"] == est.unknowns
    # memory.peak bytes at 20261026T070008Z_GEO-38-C4.log:14211, in GiB to 2 d.p.
    assert round(52_689_350_656 / 2**30, 2) == row["peak_gib"]
    assert study.result is None
    assert no_build == []
    # (Iso) the shared fit is untouched.
    assert all(r["row"] != "FC1" and r["unknowns"] != 1_525_193 for r in PRICE_TABLE)
    assert len(PRICE_TABLE) == 14


def test_wf13_s12_negative_controls(no_build, monkeypatch):
    import fem_em_solver.study as study_mod

    # Degree 3 (§9 item 715): no row, so not priced, never fitted (the degree-2
    # fit would give a finite number).
    est2 = _study(degree=3).estimate(MPI.COMM_WORLD)
    fitted = study_mod._predict(7_885_674, 2, 16)
    _say(
        f"\n[WF-13 s12] degree 3: {est2}; the per-degree fit would have said "
        f"{fitted[0]:.3f} GiB / {fitted[1]:.3f} s"
    )
    assert est2.rows_used == ()
    assert est2.unknowns == 22_523_541  # U_3 from E, F, T, not degree 2's record (§9 item 731)
    assert any("degree 3" in s for s in est2.not_priced)
    assert math.isnan(est2.predicted_peak_gib) and math.isnan(est2.predicted_seconds)
    assert math.isfinite(fitted[0]) and math.isfinite(fitted[1])

    # Planted second clinical row: the estimate moves to it.
    base = _study().estimate(MPI.COMM_WORLD)
    plant = dict(study_mod.PRICE_TABLE_CLINICAL[0], row="PLANT", peak_gib=99.0)
    monkeypatch.setattr(
        study_mod, "PRICE_TABLE_CLINICAL", study_mod.PRICE_TABLE_CLINICAL + (plant,)
    )
    moved = _study().estimate(MPI.COMM_WORLD)
    _say(
        f"[WF-13 s12] planted row: peak {base.predicted_peak_gib} -> {moved.predicted_peak_gib}, "
        f"rows {base.rows_used} -> {moved.rows_used}"
    )
    assert moved.predicted_peak_gib == 99.0 != base.predicted_peak_gib
    assert moved.rows_used == ("PLANT",)

    # A non-default geometry is refused (the records price the r1 fixture only).
    study = _study()
    study.geometry = dataclasses.replace(study.geometry, conductor_resolution=2.4e-3)
    with pytest.raises(ValueError, match="r1 fixture is priced"):
        study.estimate(MPI.COMM_WORLD)
    assert no_build == []


def test_wf13_s12_refusal_and_cli(no_build, capsys, tmp_path):
    from fem_em_solver.cli import main
    from tests.package.test_wf13_cli_estimate import _diff, _parse

    comm = MPI.COMM_WORLD
    study = _study()
    with pytest.raises(ValueError, match="xl tier only") as exc:
        study.solve(comm)
    assert "f_clinical_birdcage" in str(exc.value)
    assert no_build == []

    est = study.estimate(comm)
    out_dir = tmp_path / "out"
    capsys.readouterr()
    rc_run = main(["run", str(FC), "--output-dir", str(out_dir)])
    err_run = capsys.readouterr().err
    rc_est = main(["run", "--estimate", str(FC), "--output-dir", str(out_dir)])
    out = capsys.readouterr().out
    names = {f.name for f in dataclasses.fields(est)}
    diff = _diff(_parse(out, names), est) if comm.rank == 0 else set()
    diff = comm.bcast(diff, root=0)
    rcs_run = comm.allgather(rc_run)
    rcs_est = comm.allgather(rc_est)
    _say(f"\n[WF-13 s12] femem run: rc per rank {rcs_run}; stderr {err_run.strip()!r}")
    _say(f"[WF-13 s12] femem run --estimate: rc per rank {rcs_est}; CLI vs API diff {diff}")
    assert rcs_run == [2] * comm.size
    if comm.rank == 0:
        assert "xl tier only" in err_run
    assert rcs_est == [0] * comm.size
    assert diff == set()
    assert not out_dir.exists()
    assert no_build == []
