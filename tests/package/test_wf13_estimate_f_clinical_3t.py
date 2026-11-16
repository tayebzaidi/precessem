"""`WF-13` step 15 (§9 item 782, §10 U118): the F-clinical-3T twin through ``Study.estimate``.

Brief: docs/planning/chunks/WF-13.md, "§9 item 782".  No mesh anywhere: the
``no_build`` fixture of step 12's module is imported.

Asserted:
  (FC3) the 3 T geometry (``field_strength_t=3.0``) at degree 1 returns cells
        1 057 895 and unknowns 1 425 765, equal to `GEO-38` C8's 3 T r1 records
        (``tests/validation/test_geo38_f_clinical_3t.py:105``, ``:108``);
        ``rows_used == ("FC3",)``; seconds / GiB / ranks equal to the measured
        row, which is checked against entry 49's log lines (``:14125`` Z50 pair
        135.1 s, FC1's convention; ``:14183`` memory.peak 39.39 GiB).
  (FC4) at degree 2 it returns ``not_priced`` with NaN and U_2 = 7 335 988 exact.
  (NoB) no builder is called.
  (Ref) ``solve()`` raises "xl tier only" before any builder; ``femem run`` on
        the 3 T YAML returns 2 on every rank, no output directory;
        ``--estimate`` returns 0, and the CLI equals the API by ``repr``.
  (Key) YAML ``field_strength_t`` other than 1.5 / 3.0 is a ``StudyConfigError``;
        an unknown key beside it is still refused; 1.5 gives the default geometry.
Negative controls (asserted): the 1.5 T geometry still returns FC1 (the 3 T row
does not shadow it under the "last row of a degree wins" lookup); a planted
second 3 T row with a different ``peak_gib`` moves the 3 T estimate (the row is
what is read) and leaves FC1 unmoved.
"""

from __future__ import annotations

import dataclasses
import gzip
import math
from pathlib import Path

import pytest
from mpi4py import MPI

from tests.package.test_wf13_estimate_f_clinical import no_build  # noqa: F401  (fixture)

HERE = Path(__file__).resolve().parent
FC3T = HERE / "wf13" / "f_clinical_3t.yaml"
FC15 = HERE / "wf13" / "f_clinical_15t.yaml"
LOG = HERE.parents[1] / "docs" / "testing" / "logs" / "20261101T080008Z_GEO-38-C16.log"


def _say(msg):
    if MPI.COMM_WORLD.rank == 0:
        print(msg, flush=True)


def _log_lines():
    if LOG.exists():
        return LOG.read_text().splitlines()
    gz = LOG.with_name(LOG.name + ".gz")
    with gzip.open(gz, "rt") as fh:
        return fh.read().splitlines()


def _study(path=FC3T, degree=1):
    from fem_em_solver.study_yaml import load_study

    study = load_study(path)
    if degree != 1:
        study.degree = degree
    return study


def test_wf13_s15_fc3_and_fc4(no_build):
    from fem_em_solver.study import (
        F_CLINICAL_3T_R1_CELL_RECORD,
        F_CLINICAL_3T_R1_UNKNOWNS_RECORD,
        PRICE_TABLE,
        PRICE_TABLE_CLINICAL_3T,
        FClinicalBirdcageGeometry,
    )
    from tests.validation import test_geo38_f_clinical_3t as c8

    assert F_CLINICAL_3T_R1_CELL_RECORD == c8.F_CLINICAL_3T_R1_CELL_RECORD == 1_057_895
    assert F_CLINICAL_3T_R1_UNKNOWNS_RECORD == c8.F_CLINICAL_3T_R1_UNKNOWNS_RECORD
    study = _study()
    assert study.geometry == FClinicalBirdcageGeometry(field_strength_t=3.0)

    est = study.estimate(MPI.COMM_WORLD)
    (row,) = PRICE_TABLE_CLINICAL_3T
    _say(f"\n[WF-13 s15] FC3: {est}")
    assert est.cells == 1_057_895 == row["cells"]
    assert est.unknowns == 1_425_765 == row["unknowns"] == c8.F_CLINICAL_3T_R1_UNKNOWNS_RECORD[0]
    assert est.degree == 1
    assert est.rows_used == ("FC3",)
    assert est.ranks == row["ranks"] == 16
    assert est.predicted_seconds == row["seconds"] == 135.1
    assert est.predicted_peak_gib == row["peak_gib"] == 39.39
    assert row["build_seconds"] == 327.6
    assert est.not_priced == ("thermal transient",)

    # The row against entry 49's log, line by line (1-based line numbers).
    lines = _log_lines()
    price = lines[14125 - 1]
    assert "build 327.6 s, Z50 pair 135.1 s" in price and "at -n 16" in price
    assert lines[14182 - 1].strip() == "[XL] memory.peak bytes:"
    peak_bytes = int(lines[14183 - 1].strip())
    _say(f"[WF-13 s15] log :14125 {price[:80]!r}...; :14183 memory.peak {peak_bytes} B "
         f"= {peak_bytes / 2**30:.4f} GiB")
    assert round(peak_bytes / 2**30, 2) == row["peak_gib"]

    # FC4: degree 2 is not priced; U_2 is the record, exactly.
    est2 = _study(degree=2).estimate(MPI.COMM_WORLD)
    _say(f"[WF-13 s15] FC4 (degree 2): {est2}")
    assert est2.rows_used == ()
    assert est2.cells == 1_057_895
    assert est2.unknowns == 7_335_988 == c8.F_CLINICAL_3T_R1_UNKNOWNS_RECORD[1]
    assert any("degree 2" in s for s in est2.not_priced)
    assert math.isnan(est2.predicted_peak_gib) and math.isnan(est2.predicted_seconds)

    # (Iso) the shared fit carries no 3 T row.
    assert all(r["row"] != "FC3" and r["unknowns"] != 1_425_765 for r in PRICE_TABLE)
    assert study.result is None
    assert no_build == []


def test_wf13_s15_controls(no_build, monkeypatch):
    import fem_em_solver.study as study_mod

    # The 1.5 T geometry still returns FC1 (no shadowing by the 3 T row).
    fc1 = _study(FC15).estimate(MPI.COMM_WORLD)
    fc3 = _study().estimate(MPI.COMM_WORLD)
    _say(f"\n[WF-13 s15] control 1.5 T: rows {fc1.rows_used} {fc1.predicted_seconds} s "
         f"{fc1.predicted_peak_gib} GiB; 3 T: rows {fc3.rows_used}")
    assert fc1.rows_used == ("FC1",) and fc1.unknowns == 1_525_193
    assert fc1.predicted_seconds == 203.3 and fc1.predicted_peak_gib == 49.07
    assert fc1.cells != fc3.cells

    # Planted second 3 T row: the 3 T estimate moves, FC1 does not.
    plant = dict(study_mod.PRICE_TABLE_CLINICAL_3T[0], row="PLANT", peak_gib=99.0)
    monkeypatch.setattr(study_mod, "PRICE_TABLE_CLINICAL_3T",
                        study_mod.PRICE_TABLE_CLINICAL_3T + (plant,))
    moved = _study().estimate(MPI.COMM_WORLD)
    fc1_after = _study(FC15).estimate(MPI.COMM_WORLD)
    _say(f"[WF-13 s15] planted 3 T row: peak {fc3.predicted_peak_gib} -> "
         f"{moved.predicted_peak_gib}; FC1 peak {fc1_after.predicted_peak_gib}")
    assert moved.predicted_peak_gib == 99.0 and moved.rows_used == ("PLANT",)
    assert fc1_after == fc1
    assert no_build == []


def test_wf13_s15_yaml_key(tmp_path):
    from fem_em_solver.study import FClinicalBirdcageGeometry
    from fem_em_solver.study_yaml import StudyConfigError, load_study

    text = FC3T.read_text()
    p15 = tmp_path / "b15.yaml"
    p15.write_text(text.replace("field_strength_t: 3.0", "field_strength_t: 1.5"))
    assert load_study(p15).geometry == FClinicalBirdcageGeometry()

    bad = tmp_path / "b7.yaml"
    bad.write_text(text.replace("field_strength_t: 3.0", "field_strength_t: 7.0"))
    with pytest.raises(StudyConfigError, match="field_strength_t") as exc:
        load_study(bad)
    _say(f"\n[WF-13 s15] 7.0 refused: {exc.value}")

    unk = tmp_path / "unk.yaml"
    unk.write_text(text.replace("field_strength_t: 3.0",
                                "field_strength_t: 3.0\n  ring_radius: 0.352"))
    with pytest.raises(StudyConfigError, match="ring_radius") as exc2:
        load_study(unk)
    _say(f"[WF-13 s15] unknown key still refused: {exc2.value}")


def test_wf13_s15_refusal_and_cli(no_build, capsys, tmp_path):
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
    rc_run = main(["run", str(FC3T), "--output-dir", str(out_dir)])
    err_run = capsys.readouterr().err
    rc_est = main(["run", "--estimate", str(FC3T), "--output-dir", str(out_dir)])
    out = capsys.readouterr().out
    names = {f.name for f in dataclasses.fields(est)}
    diff = _diff(_parse(out, names), est) if comm.rank == 0 else set()
    diff = comm.bcast(diff, root=0)
    rcs_run = comm.allgather(rc_run)
    rcs_est = comm.allgather(rc_est)
    _say(f"\n[WF-13 s15] femem run (3 T): rc per rank {rcs_run}; stderr {err_run.strip()!r}")
    _say(f"[WF-13 s15] femem run --estimate (3 T): rc per rank {rcs_est}; CLI vs API diff {diff}")
    assert rcs_run == [2] * comm.size
    if comm.rank == 0:
        assert "xl tier only" in err_run
    assert rcs_est == [0] * comm.size
    assert diff == set()
    assert est.rows_used == ("FC3",)
    assert not out_dir.exists()
    assert no_build == []
