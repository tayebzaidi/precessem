"""`WF-13` step 13 (§9 item 715, chain U106): the F-clinical-1.5T degree-2 price row FC2.

Brief: docs/planning/chunks/WF-13.md, "§9 item 715".  No mesh is built anywhere
in this module (``no_build``, imported from step 12's module).

Asserted:
  (FC2) ``Study.estimate()`` at degree 2 returns cells 1 147 366, unknowns
        7 885 674 equal to ``F_CLINICAL_R1_UNKNOWNS_RECORD[1]``,
        ``rows_used == ("FC2",)``; seconds / GiB / ranks equal to the row, and
        the row equal to `xxl` entry 50's log (build, Z50 pair, memory.peak).
  (NoB) no builder is called.
  (Ref) ``solve()`` at degree 2 raises before any builder with both tiers named;
        ``femem run`` returns 2 on every rank, no output directory;
        ``--estimate`` returns 0.
  (CLI) the CLI's printed estimate equals the API's field by field (repr).
Negative controls (asserted): a planted second FC2 row with a different
``peak_gib`` changes the estimate (the row is what is read); the degree-1 row
misses on unknowns by 5.17×.
"""

from __future__ import annotations

import dataclasses
import gzip
import re
from pathlib import Path

import pytest
import yaml
from mpi4py import MPI

from tests.package.test_wf13_estimate_f_clinical import FC, _study, no_build  # noqa: F401

LOG = (Path(__file__).resolve().parents[2]
       / "docs/testing/logs/20261031T070008Z_GEO-38-C5a.log")
FC2_LITERALS = dict(row="FC2", degree=2, ranks=16, cells=1147366, unknowns=7885674,
                    seconds=2892.1, build_seconds=341.8, peak_gib=472.42)
# 7 885 674 / 1 525 193 to two decimals, the brief's control.
UNKNOWN_MISS = 5.17


def _say(msg):
    if MPI.COMM_WORLD.rank == 0:
        print(msg, flush=True)


def _log_lines():
    """The log's lines; the weekly sweep gzips a cited log at 7 days (retention §1)."""
    if LOG.exists():
        return LOG.read_text(errors="replace").splitlines()
    with gzip.open(f"{LOG}.gz", "rt", errors="replace") as fh:
        return fh.read().splitlines()


def test_wf13_s13_fc2_estimate_reads_row(no_build):  # noqa: F811
    from fem_em_solver.study import (
        F_CLINICAL_R1_CELL_RECORD,
        F_CLINICAL_R1_UNKNOWNS_RECORD,
        PRICE_TABLE,
        PRICE_TABLE_CLINICAL,
    )
    from tests.validation import test_geo38_f_clinical_build as geo38

    (fc1, row) = PRICE_TABLE_CLINICAL
    assert row == FC2_LITERALS
    est = _study(degree=2).estimate(MPI.COMM_WORLD)
    _say(f"\n[WF-13 s13] FC2: {est}")
    assert est.cells == 1_147_366 == F_CLINICAL_R1_CELL_RECORD == row["cells"]
    assert est.unknowns == 7_885_674 == F_CLINICAL_R1_UNKNOWNS_RECORD[1] == row["unknowns"]
    assert F_CLINICAL_R1_UNKNOWNS_RECORD == geo38.F_CLINICAL_R1_UNKNOWNS_RECORD
    assert est.degree == 2
    assert est.rows_used == ("FC2",)
    assert est.not_priced == ("thermal transient",)
    assert est.ranks == row["ranks"] == 16
    assert est.predicted_seconds == row["seconds"] == 2892.1
    assert est.predicted_peak_gib == row["peak_gib"] == 472.42

    # The row is its log, build excluded from the seconds (FC1's convention).
    lines = _log_lines()
    price = lines[12595 - 1]
    m = re.search(r"price: build ([\d.]+) s, Z50 pair ([\d.]+) s", price)
    assert m is not None, price
    assert (float(m[1]), float(m[2])) == (row["build_seconds"], row["seconds"])
    peak_bytes = int(lines[12650 - 1])
    assert lines[12649 - 1].startswith("[XL] memory.peak bytes")
    assert peak_bytes == 507_256_434_688
    assert round(peak_bytes / 2**30, 2) == row["peak_gib"]
    assert "- Status: 0" in lines[12654 - 1]
    _say(f"[WF-13 s13] log: {m[0]}; memory.peak {peak_bytes} B")

    # The fit never sees a clinical row: FC1's (Iso) at degree 2.
    assert all(r["row"] != "FC2" and r["unknowns"] != 7_885_674 for r in PRICE_TABLE)
    assert fc1["row"] == "FC1"
    assert no_build == []


def test_wf13_s13_negative_controls(no_build, monkeypatch):  # noqa: F811
    import fem_em_solver.study as study_mod

    # Planted second FC2 row: the estimate moves to it.
    base = _study(degree=2).estimate(MPI.COMM_WORLD)
    plant = dict(study_mod.PRICE_TABLE_CLINICAL[1], row="PLANT", peak_gib=99.0)
    with monkeypatch.context() as patch:  # `no_build`'s patches stay in force
        patch.setattr(study_mod, "PRICE_TABLE_CLINICAL",
                      study_mod.PRICE_TABLE_CLINICAL + (plant,))
        moved = _study(degree=2).estimate(MPI.COMM_WORLD)
    _say(f"\n[WF-13 s13] planted row: peak {base.predicted_peak_gib} -> "
         f"{moved.predicted_peak_gib}, rows {base.rows_used} -> {moved.rows_used}")
    assert moved.predicted_peak_gib == 99.0 != base.predicted_peak_gib
    assert moved.rows_used == ("PLANT",)

    # The degree-1 row is not the degree-2 reading.
    d1 = _study().estimate(MPI.COMM_WORLD)
    ratio = base.unknowns / d1.unknowns
    _say(f"[WF-13 s13] degree-1 row misses unknowns by {ratio:.4f}x; "
         f"seconds {d1.predicted_seconds} vs {base.predicted_seconds}")
    assert round(ratio, 2) == UNKNOWN_MISS
    assert d1.rows_used == ("FC1",) and base.rows_used == ("FC2",)


def test_wf13_s13_refusal_and_cli(no_build, capsys, tmp_path):  # noqa: F811
    from fem_em_solver.cli import main
    from tests.package.test_wf13_cli_estimate import _diff, _parse

    comm = MPI.COMM_WORLD
    study = _study(degree=2)
    with pytest.raises(ValueError, match="xl tier only") as exc:
        study.solve(comm)
    msg = str(exc.value)
    _say(f"\n[WF-13 s13] solve refused: {msg}")
    assert "f_clinical_birdcage" in msg and "<= 3 root-grid nodes" in msg
    assert "xxl tier" in msg and "column set" in msg
    assert no_build == []

    # A rank-local temporary directory would differ per rank; rank 0's is shared.
    base = comm.bcast(str(tmp_path) if comm.rank == 0 else None, root=0)
    file_d2 = Path(base) / "f_clinical_15t_d2.yaml"
    if comm.rank == 0:
        doc = yaml.safe_load(FC.read_text())
        doc["study"]["degree"] = 2
        file_d2.write_text(yaml.safe_dump(doc))
    comm.Barrier()

    est = study.estimate(comm)
    out_dir = Path(base) / "out"
    capsys.readouterr()
    rc_run = main(["run", str(file_d2), "--output-dir", str(out_dir)])
    err_run = capsys.readouterr().err
    rc_est = main(["run", "--estimate", str(file_d2), "--output-dir", str(out_dir)])
    out = capsys.readouterr().out
    names = {f.name for f in dataclasses.fields(est)}
    diff = _diff(_parse(out, names), est) if comm.rank == 0 else set()
    diff = comm.bcast(diff, root=0)
    rcs_run = comm.allgather(rc_run)
    rcs_est = comm.allgather(rc_est)
    _say(f"[WF-13 s13] femem run: rc per rank {rcs_run}; stderr {err_run.strip()!r}")
    _say(f"[WF-13 s13] femem run --estimate: rc per rank {rcs_est}; CLI vs API diff {diff}")
    assert rcs_run == [2] * comm.size
    if comm.rank == 0:
        assert "xl tier only" in err_run and "xxl tier" in err_run
    assert rcs_est == [0] * comm.size
    assert diff == set()
    assert not out_dir.exists()
    assert no_build == []
