"""`WF-13` step 8 (§9 item 480): ``femem run --estimate`` on the F-human YAML.

Brief: docs/planning/chunks/WF-13.md, "§9 item 480".  `WF-16` step 8's refusal
(``_refuse_f_human`` raises ``TypeError`` before any mesh) is not touched; only
its CLI exit path.  W0 (``20261015T214205Z_WF-13-s8-W0-observe.log``, before the
fix) observed an uncaught ``TypeError`` traceback on both ranks, status 1.

Asserted, in the degree-1 refusal node's pattern
(``test_wf13_cli_estimate.py::test_wf13_cli_estimate_degree1_passes_api_error``):
  ``Study.estimate`` on the loaded F-human Study raises ``TypeError`` (the API
  message); with solves and the geometry's ``build`` patched to raise,
  ``femem run --estimate`` returns 2 on every rank; rank 0's stderr carries
  exactly ``femem: error: <API message>``, other ranks none; the output
  directory does not exist.

Negative control (asserted, backed by W0's log): the pre-fix handler (a local
copy of ``cli._estimate`` catching ``ValueError`` only) re-raises ``TypeError``.

Re-pointed 2026-10-20 by `WF-13` step 10 (§9 item 562), which lifted the
estimate refusal: the node now gates the F-human CLI = API path on a stand-in
mesh (``build`` patched to a 2x2x2 unit cube), and its control asserts that the
``_refuse_f_human`` helper, called directly, still raises item 480's message.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest
from mpi4py import MPI

from tests.complex_mode import complex_only

FH = Path(__file__).resolve().parent / "wf1" / "f_human_tuned_copper.yaml"


def _say(msg):
    if MPI.COMM_WORLD.rank == 0:
        print(msg, flush=True)


def _pre_fix_estimate(path, comm) -> int:
    """Local copy of ``cli._estimate`` as of 24088335 (catches ``ValueError`` only)."""
    import dataclasses

    from fem_em_solver.study_yaml import load_study

    study = load_study(path)
    try:
        est = study.estimate(comm)
    except ValueError as exc:
        if comm.rank == 0:
            print(f"femem: error: {exc}", file=sys.stderr, flush=True)
        return 2
    if comm.rank == 0:
        for f in dataclasses.fields(est):
            print(f"{f.name}: {getattr(est, f.name)!r}")
        sys.stdout.flush()
    return 0


# FH2: §9 item 585 — the helper's message updated ((R), review ruling
# 2026-10-22 03:00); the control's meaning is unchanged.
REFUSAL_MESSAGE = (
    "estimate is not defined on an FHumanBirdcageGeometry Study "
    "(WF-16 step 8: S-parameters only; estimate() is priced, WF-13 steps 10-11)"
)


@complex_only
def test_wf13_cli_estimate_f_human_refusal(monkeypatch, capsys, tmp_path):
    # FH1: §9 item 562 — re-pointed to the lifted contract on a stand-in mesh.
    import dataclasses

    import dolfinx.mesh
    from dolfinx import fem

    import fem_em_solver.study_yaml as study_yaml
    from fem_em_solver.cli import main
    from fem_em_solver.study import FHumanBirdcageGeometry, Study, _refuse_f_human
    from fem_em_solver.study_yaml import load_study
    from tests.package.test_wf13_cli_estimate import _diff, _parse

    comm = MPI.COMM_WORLD

    def _boom(*a, **k):
        raise AssertionError("solve reached")

    def _stand_in(self, c=MPI.COMM_WORLD, *a, **k):
        return {"mesh": dolfinx.mesh.create_unit_cube(c, 2, 2, 2)}

    monkeypatch.setattr(Study, "solve", _boom)
    monkeypatch.setattr(study_yaml, "run_study_file", _boom)
    monkeypatch.setattr(FHumanBirdcageGeometry, "build", _stand_in)

    study = load_study(FH)
    assert isinstance(study.geometry, FHumanBirdcageGeometry)
    est = study.estimate(comm)
    cube = dolfinx.mesh.create_unit_cube(comm, 2, 2, 2)
    V = fem.functionspace(cube, ("N1curl", 1))
    n_ref = int(V.dofmap.index_map.size_global * V.dofmap.index_map_bs)
    assert est.cells == 48
    assert est.unknowns == n_ref
    assert "FH1" in est.rows_used
    assert est.not_priced == ("thermal transient",)
    assert study.result is None

    out_dir = tmp_path / "out"
    capsys.readouterr()
    rc = main(["run", "--estimate", str(FH), "--output-dir", str(out_dir)])
    out = capsys.readouterr().out
    rcs = comm.allgather(rc)
    names = {f.name for f in dataclasses.fields(est)}
    diff = _diff(_parse(out, names), est) if comm.rank == 0 else set()
    diff = comm.bcast(diff, root=0)
    _say(f"\n[WF-13 s10] stand-in API: cells {est.cells} unknowns {est.unknowns} "
         f"(ref {n_ref}) rows {est.rows_used} not_priced {est.not_priced}")
    _say(f"[WF-13 s10] femem run --estimate F-human (stand-in): rc per rank {rcs}; diff {diff}")
    assert rcs == [0] * comm.size
    assert diff == set()
    assert not out_dir.exists()

    # Control (asserted, by construction): the helper still refuses when called.
    with pytest.raises(TypeError) as ctl:
        _refuse_f_human(study, "estimate")
    assert str(ctl.value) == REFUSAL_MESSAGE
    _say(f"[WF-13 s10] control: _refuse_f_human(study, 'estimate') -> TypeError {ctl.value}")
