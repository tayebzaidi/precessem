"""`WF-13` step 16 (§9 item 783, §10 U119): entry 48's rod in the clinical box through ``Study.estimate``.

Brief: docs/planning/chunks/WF-13.md, "§9 item 783".  No mesh anywhere: the
``no_build`` fixture of step 12's module is imported.

Asserted:
  (Rod)  ``F_CLINICAL_ROD_ENTRY48`` equals `GEO-36`'s ``ROD`` moved to
         ``_rod_centre()``; the cell / unknowns records equal
         ``test_geo36_rod_f_clinical.py:410`` / ``:412``.
  (FCR)  ``Study.estimate()`` returns FC-ROD1: cells 1 217 025, U 1 606 744,
         ranks 8, 263.5 s (arm M, FC1's one-solve convention), 42.39 GiB, each
         checked against entry 48's log lines (``:12503``, ``:12626`` / ``:12627``).
  (NoB)  no builder is called.
  (Ref)  ``solve``, ``solve_sweep`` and the field route raise "xl tier only";
         ``femem run`` returns 2 on every rank; ``--estimate`` returns 0 and the
         CLI equals the API by ``repr``.
  (Only) a second implant (the disc), a moved rod, the four-port rod, and the
         rod at 3 T are refused; YAML keeps the unchanged message.
Negative controls (asserted): the rod-free FC1 is unchanged; a planted second
rod row moves the rod estimate and leaves FC1 equal.
"""

from __future__ import annotations

import dataclasses
import gzip
from pathlib import Path

import pytest
import yaml
from mpi4py import MPI

from tests.package.test_wf13_estimate_f_clinical import no_build  # noqa: F401  (fixture)

HERE = Path(__file__).resolve().parent
ROD_YAML = HERE / "wf13" / "f_clinical_rod.yaml"
FC15 = HERE / "wf13" / "f_clinical_15t.yaml"
LOG = HERE.parents[1] / "docs" / "testing" / "logs" / "20261102T080008Z_GEO-36-step9-I36.log"
UNCHANGED = "an implant is defined only on a four_port_birdcage, not 'f_clinical_birdcage'"


def _say(msg):
    if MPI.COMM_WORLD.rank == 0:
        print(msg, flush=True)


def _log_lines():
    if LOG.exists():
        return LOG.read_text().splitlines()
    with gzip.open(LOG.with_name(LOG.name + ".gz"), "rt") as fh:
        return fh.read().splitlines()


def _study(path=ROD_YAML):
    from fem_em_solver.study_yaml import load_study

    return load_study(path)


def _write(tmp_path, name, mapping):
    path = tmp_path / f"{name}-{MPI.COMM_WORLD.rank}.yaml"
    path.write_text(yaml.safe_dump(mapping, sort_keys=False))
    return path


def test_wf13_s16_fc_rod1(no_build):
    from fem_em_solver.study import (
        F_CLINICAL_ROD_CELL_RECORD,
        F_CLINICAL_ROD_ENTRY48,
        F_CLINICAL_ROD_UNKNOWNS_RECORD,
        PRICE_TABLE,
        PRICE_TABLE_CLINICAL_ROD,
        FClinicalBirdcageGeometry,
    )
    from tests.validation import test_geo36_rod_f_clinical as g36

    assert F_CLINICAL_ROD_ENTRY48 == dict(g36.ROD, centre=g36._rod_centre())
    assert F_CLINICAL_ROD_CELL_RECORD == g36.F_CLINICAL_ROD_CELL_RECORD == 1_217_025
    assert F_CLINICAL_ROD_UNKNOWNS_RECORD == g36.F_CLINICAL_ROD_UNKNOWNS_RECORD
    study = _study()
    assert study.geometry == FClinicalBirdcageGeometry(implant_rod=F_CLINICAL_ROD_ENTRY48)

    est = study.estimate(MPI.COMM_WORLD)
    (row,) = PRICE_TABLE_CLINICAL_ROD
    _say(f"\n[WF-13 s16] FC-ROD1: {est}")
    assert est.cells == 1_217_025 == row["cells"]
    assert est.unknowns == 1_606_744 == row["unknowns"] == F_CLINICAL_ROD_UNKNOWNS_RECORD[0]
    assert est.degree == 1
    assert est.rows_used == ("FC-ROD1",)
    assert est.ranks == row["ranks"] == 8
    assert est.predicted_seconds == row["seconds"] == 263.5
    assert est.predicted_peak_gib == row["peak_gib"] == 42.39
    assert row["build_seconds"] == 461.7
    assert est.not_priced == ("thermal transient",)

    lines = _log_lines()
    price = lines[12503 - 1]
    assert "price: build 461.7 s; arm M 263.5 s, arm R 262.3 s" in price
    assert lines[12626 - 1].strip() == "[XL] memory.peak bytes:"
    peak_bytes = int(lines[12627 - 1].strip())
    _say(f"[WF-13 s16] log :12503 {price!r}; :12627 memory.peak {peak_bytes} B "
         f"= {peak_bytes / 2**30:.4f} GiB")
    assert round(peak_bytes / 2**30, 2) == row["peak_gib"]

    assert all(r["row"] != "FC-ROD1" and r["unknowns"] != 1_606_744 for r in PRICE_TABLE)
    assert study.result is None
    assert no_build == []


def test_wf13_s16_refusals_and_cli(no_build, capsys, tmp_path):
    from fem_em_solver.cli import main
    from fem_em_solver.study import _study_solve_port_fields
    from tests.package.test_wf13_cli_estimate import _diff, _parse

    comm = MPI.COMM_WORLD
    study = _study()
    refused = []
    for name, call in (("solve", lambda: study.solve(comm)),
                       ("solve_sweep", lambda: study.solve_sweep([64.0e6], comm)),
                       ("field route", lambda: _study_solve_port_fields(study, comm))):
        with pytest.raises(ValueError, match="xl tier only") as exc:
            call()
        refused.append(f"[WF-13 s16] (Ref) {name}: {exc.value}")
    assert no_build == []

    est = study.estimate(comm)
    out_dir = tmp_path / "out"
    capsys.readouterr()
    rc_run = main(["run", str(ROD_YAML), "--output-dir", str(out_dir)])
    err_run = capsys.readouterr().err
    rc_est = main(["run", "--estimate", str(ROD_YAML), "--output-dir", str(out_dir)])
    out = capsys.readouterr().out
    names = {f.name for f in dataclasses.fields(est)}
    diff = _diff(_parse(out, names), est) if comm.rank == 0 else set()
    diff = comm.bcast(diff, root=0)
    rcs_run, rcs_est = comm.allgather(rc_run), comm.allgather(rc_est)
    _say("\n" + "\n".join(refused))
    _say(f"[WF-13 s16] femem run (rod): rc per rank {rcs_run}; stderr {err_run.strip()!r}")
    _say(f"[WF-13 s16] femem run --estimate (rod): rc per rank {rcs_est}; CLI vs API diff {diff}")
    assert rcs_run == [2] * comm.size
    if comm.rank == 0:
        assert "xl tier only" in err_run
    assert rcs_est == [0] * comm.size
    assert diff == set()
    assert est.rows_used == ("FC-ROD1",)
    assert not out_dir.exists()
    assert no_build == []


def test_wf13_s16_only_this_rod(tmp_path):
    from fem_em_solver.study import F_CLINICAL_ROD_ENTRY48, FClinicalBirdcageGeometry
    from fem_em_solver.study_yaml import StudyConfigError, load_study
    from tests.mesh.test_birdcage_implant_rod import ROD as FOUR_PORT_ROD

    good = yaml.safe_load(ROD_YAML.read_text())
    rod = good["implant"]["rod"]
    disc = {"radius_m": 0.01, "thickness_m": 0.002, "centre_m": [0.105, 0.0, 0.0]}
    four_port = {"radius_m": FOUR_PORT_ROD["radius"], "length_m": FOUR_PORT_ROD["length"],
                 "centre_m": list(FOUR_PORT_ROD["centre"]), "axis": list(FOUR_PORT_ROD["axis"])}
    cases = {
        "disc_beside": dict(good, implant={"rod": rod, "disc": disc}),
        "disc_alone": dict(good, implant={"disc": disc}),
        "moved": dict(good, implant={"rod": dict(rod, centre_m=[0.106, 0.0, 0.0])}),
        "longer": dict(good, implant={"rod": dict(rod, length_m=0.031)}),
        "four_port_rod": dict(good, implant={"rod": four_port}),
        "at_3t": dict(good, geometry={"kind": "f_clinical_birdcage", "field_strength_t": 3.0}),
    }
    for name, mapping in cases.items():
        with pytest.raises(StudyConfigError) as err:
            load_study(_write(tmp_path, name, mapping))
        _say(f"[WF-13 s16] (Only) {name}: {err.value.field}: {err.value.reason}")
        assert err.value.field == "implant", name
        assert err.value.reason == UNCHANGED, name

    for name, kw in (("disc", dict(implant_rod={"radius": 0.01, "thickness": 0.002})),
                     ("four_port_rod", dict(implant_rod=dict(FOUR_PORT_ROD))),
                     ("at_3t", dict(implant_rod=F_CLINICAL_ROD_ENTRY48, field_strength_t=3.0))):
        with pytest.raises(ValueError, match="entry 48's rod"):
            FClinicalBirdcageGeometry(**kw)
        _say(f"[WF-13 s16] (Only) Python {name}: refused")


def test_wf13_s16_controls(no_build, monkeypatch):
    import fem_em_solver.study as study_mod

    fc1 = _study(FC15).estimate(MPI.COMM_WORLD)
    rod = _study().estimate(MPI.COMM_WORLD)
    _say(f"\n[WF-13 s16] control FC1: rows {fc1.rows_used} {fc1.predicted_seconds} s "
         f"{fc1.predicted_peak_gib} GiB, U {fc1.unknowns}; rod: rows {rod.rows_used}")
    assert fc1.rows_used == ("FC1",) and fc1.unknowns == 1_525_193 and fc1.cells == 1_147_366
    assert fc1.predicted_seconds == 203.3 and fc1.predicted_peak_gib == 49.07

    plant = dict(study_mod.PRICE_TABLE_CLINICAL_ROD[0], row="PLANT", peak_gib=99.0)
    monkeypatch.setattr(study_mod, "PRICE_TABLE_CLINICAL_ROD",
                        study_mod.PRICE_TABLE_CLINICAL_ROD + (plant,))
    moved = _study().estimate(MPI.COMM_WORLD)
    fc1_after = _study(FC15).estimate(MPI.COMM_WORLD)
    _say(f"[WF-13 s16] planted rod row: peak {rod.predicted_peak_gib} -> "
         f"{moved.predicted_peak_gib}; FC1 peak {fc1_after.predicted_peak_gib}")
    assert moved.predicted_peak_gib == 99.0 and moved.rows_used == ("PLANT",)
    assert fc1_after == fc1
    assert no_build == []
