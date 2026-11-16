"""`WF-13` step 4 (§9 item 325, chain U28): ``femem run --estimate <file>``.

Brief: docs/planning/chunks/WF-13.md, "§9 item 325".  Run in-process under
``mpiexec -n 2`` (rank-local bugs show only there).  Asserted:

(E) field-by-field equality, CLI vs API, same process, same comm: rank 0 parses
    the printed ``name: repr(value)`` lines; for every field of
    ``ResourceEstimate`` ``repr(printed) == repr(api)``; the printed field set is
    the dataclass's field set (no missing / duplicate line).  Verdict bcast.
(R) the gated record through the CLI: ``unknowns == WF15_DEG2_UNKNOWNS_RECORD``,
    ``cells == WF15_CELLS_RECORD`` (step 2's (U), imported), ``degree == 2``,
    ``ranks == comm.size``.
(S) no solve: ``Study.solve`` and ``study_yaml.run_study_file`` patched to raise;
    ``--estimate`` returns 0 through them; ``--output-dir`` is not created.
(X) refusals: ``--ranks 2 --estimate`` returns 2 naming both flags, with a
    marker-writing stub launcher that is never reached; a degree-1 file returns
    2 on every rank with the API's ``ValueError`` message unchanged (review
    ruling (1), 2026-10-08).  A thermal-only ball file (`WF-13` step 9, §9 item
    507) returns 2 on every rank with the API's pre-mesh ``TypeError`` message
    and no mesh built; its control (the refusal patched to a no-op) reproduces
    the uncaught ``AttributeError`` (``sphere_tag``) after one mesh build per rank
    that the item's W0 logged before the fix.

Negative controls (asserted): (c1) the solve detector fires on a plain
``femem run``; (c2) the (E) comparator reports exactly ``ranks`` against a
one-field-perturbed estimate; (c3) the degree-1 count (``PRICE_TABLE`` P1)
differs from the printed degree-2 count.
"""

from __future__ import annotations

import ast
import dataclasses
import os
from pathlib import Path

import pytest
from mpi4py import MPI

from tests.complex_mode import complex_only
from tests.package.test_wf13_study_estimate import WF15_CELLS_RECORD, WF15_DEG2_UNKNOWNS_RECORD
from tests.package.test_wf1_study_yaml import SPHERE_FILE

DEG2_FILE = Path(__file__).resolve().parent / "wf13" / "sphere_mat4_deg2.yaml"
DEG1_MESSAGE = "fewer than two priced rows at degree 1"
THERMAL_FILE = Path(__file__).resolve().parent / "wf1" / "thermal_ball_duty_a.yaml"
THERMAL_MESSAGE = (
    "estimate is not defined on a thermal-only Study "
    "(ThermalBallGeometry has no EM solve; use temperature())"
)


def _say(msg):
    if MPI.COMM_WORLD.rank == 0:
        print(msg, flush=True)


def _patch_solves(monkeypatch):
    import fem_em_solver.study_yaml as study_yaml
    from fem_em_solver.study import Study

    def _boom(*a, **k):
        raise AssertionError("solve reached")

    monkeypatch.setattr(Study, "solve", _boom)
    monkeypatch.setattr(study_yaml, "run_study_file", _boom)


def _parse(out: str, names) -> dict:
    """Printed ``name: repr`` lines whose key is a field name -> list of raw values."""
    got = {}
    for ln in out.splitlines():
        key, sep, val = ln.partition(": ")
        if sep and key in names:
            got.setdefault(key, []).append(val)
    return got


def _diff(parsed: dict, est) -> set:
    """Fields whose printed repr differs from ``repr`` of ``est``'s value."""
    return {
        f.name
        for f in dataclasses.fields(est)
        if repr(ast.literal_eval(parsed[f.name][0])) != repr(getattr(est, f.name))
    }


@complex_only
def test_wf13_cli_estimate_equals_api_and_record(monkeypatch, capsys, tmp_path):
    """(E), (R), (S), (c1), (c2), (c3)."""
    from fem_em_solver.cli import main
    from fem_em_solver.study import PRICE_TABLE, ResourceEstimate
    from fem_em_solver.study_yaml import load_study

    comm = MPI.COMM_WORLD
    _patch_solves(monkeypatch)
    out_dir = tmp_path / "out"
    capsys.readouterr()
    rc = main(["run", "--estimate", str(DEG2_FILE), "--output-dir", str(out_dir)])
    out = capsys.readouterr().out
    api = load_study(DEG2_FILE).estimate(comm)
    names = [f.name for f in dataclasses.fields(ResourceEstimate)]

    verdict = None
    if comm.rank == 0:
        parsed = _parse(out, names)
        print("\n[WF-13 s4] printed by `femem run --estimate`:\n" + out.rstrip(), flush=True)
        print(f"[WF-13 s4] API estimate: {api!r}", flush=True)
        complete = set(parsed) == set(names) and all(len(v) == 1 for v in parsed.values())
        differ = _diff(parsed, api) if complete else None
        c2 = _diff(parsed, dataclasses.replace(api, ranks=api.ranks + 1)) if complete else None
        cli = {k: ast.literal_eval(v[0]) for k, v in parsed.items()} if complete else {}
        verdict = dict(complete=complete, differ=differ, c2=c2, cli=cli)
        print(
            f"[WF-13 s4] (E) field set complete {complete}; differing fields {differ}; "
            f"(c2) vs ranks+1: {c2}\n"
            f"  (R) unknowns {cli.get('unknowns')} vs {WF15_DEG2_UNKNOWNS_RECORD}; "
            f"cells {cli.get('cells')} vs {WF15_CELLS_RECORD}; degree {cli.get('degree')}; "
            f"ranks {cli.get('ranks')} vs comm.size {comm.size}",
            flush=True,
        )
    else:
        # rank 1 prints no estimate line
        assert _parse(out, names) == {}, out
    verdict = comm.bcast(verdict, root=0)

    # (S)
    assert rc == 0
    assert not out_dir.exists()
    # (E)
    assert verdict["complete"]
    assert verdict["differ"] == set()
    # (c2) the comparator can fail: exactly the perturbed field
    assert verdict["c2"] == {"ranks"}
    # (R)
    cli = verdict["cli"]
    assert cli["unknowns"] == WF15_DEG2_UNKNOWNS_RECORD
    assert cli["cells"] == WF15_CELLS_RECORD
    assert cli["degree"] == 2
    assert cli["ranks"] == comm.size
    # (c3) degree confusion: the degree-1 count is not what was printed
    p1 = [r for r in PRICE_TABLE if r["row"] == "P1"][0]["unknowns"]
    _say(f"  (c3) PRICE_TABLE P1 degree-1 unknowns {p1} vs printed {cli['unknowns']}")
    assert p1 != cli["unknowns"]

    # (c1) the solve detector is live: plain `femem run` reaches it
    with pytest.raises(AssertionError, match="solve reached"):
        main(["run", str(DEG2_FILE), "--output-dir", str(tmp_path / "c1")])
    _say("  (c1) plain `femem run` hit 'solve reached' (detector live)")


@complex_only
def test_wf13_cli_estimate_degree1_passes_api_error(monkeypatch, capsys, tmp_path):
    """(X) degree 1: rc 2 on every rank, the API's message unchanged, no solve.

    Since `WF-13` step 6 (§9 item 471a) completed ``PRICE_TABLE`` P3, a degree-1
    file is priced; the refusal path is kept by restricting the table in-test to
    P1 only (``_predict`` patched — its default ``table=`` is bound at definition).
    """
    import fem_em_solver.study as study_mod
    from fem_em_solver.cli import main
    from fem_em_solver.study_yaml import load_study

    orig = study_mod._predict
    p1_only = tuple(row for row in study_mod.PRICE_TABLE if row["row"] not in ("P3", "T49", "FH1"))  # T49: item 479; FH1: §9 item 562
    monkeypatch.setattr(
        study_mod, "_predict", lambda u, d, r, n_rhs=1: orig(u, d, r, n_rhs, table=p1_only)
    )
    comm = MPI.COMM_WORLD
    with pytest.raises(ValueError) as api_exc:
        load_study(SPHERE_FILE).estimate(comm)
    assert str(api_exc.value) == DEG1_MESSAGE

    _patch_solves(monkeypatch)
    out_dir = tmp_path / "out"
    capsys.readouterr()
    rc = main(["run", "--estimate", str(SPHERE_FILE), "--output-dir", str(out_dir)])
    err = capsys.readouterr().err
    lines = [ln for ln in err.splitlines() if ln.startswith("femem: error:")]
    _say(f"\n[WF-13 s4] (X) degree 1: rc {rc}; rank-0 stderr {lines!r}")
    rcs = comm.allgather(rc)
    assert rcs == [2] * comm.size
    if comm.rank == 0:
        assert lines == [f"femem: error: {api_exc.value}"]
    else:
        assert lines == []
    assert not out_dir.exists()


def test_wf13_cli_estimate_thermal_ball_refusal(monkeypatch, capsys, tmp_path):
    """(X) thermal-only ball (`WF-13` step 9, §9 item 507): rc 2 on every rank,
    the API's pre-mesh ``TypeError`` message unchanged, ``ThermalBallGeometry.build``
    called 0 times.  Control: the refusal patched to a no-op reproduces W0's
    uncaught ``AttributeError`` (``sphere_tag``) after one real build per rank
    (`20261017T093035Z_WF-13-s9-W0-observe.log`, status 1)."""
    import fem_em_solver.study as study_mod
    from fem_em_solver.cli import main
    from fem_em_solver.study import ThermalBallGeometry
    from fem_em_solver.study_yaml import load_study

    comm = MPI.COMM_WORLD

    # (1) control: refusal seam a no-op, real build counted
    ctrl_calls = []
    real_build = ThermalBallGeometry.build
    with monkeypatch.context() as m:
        m.setattr(study_mod, "_refuse_thermal_only", lambda study, method: None)

        def _counted(self, *a, **k):
            ctrl_calls.append(1)
            return real_build(self, *a, **k)

        m.setattr(ThermalBallGeometry, "build", _counted)
        with pytest.raises(AttributeError, match="sphere_tag") as ctrl_exc:
            load_study(THERMAL_FILE).estimate(comm)
    ctrl_counts = comm.allgather(len(ctrl_calls))
    _say(f"\n[WF-13 s9] control (refusal no-op): {type(ctrl_exc.value).__name__}: "
         f"{ctrl_exc.value}; build calls per rank {ctrl_counts}")
    assert ctrl_counts == [1] * comm.size

    # (2) the API refuses before any mesh
    _patch_solves(monkeypatch)
    calls = []

    def _no_build(self, *a, **k):
        calls.append(1)
        raise AssertionError("ThermalBallGeometry.build reached")

    monkeypatch.setattr(ThermalBallGeometry, "build", _no_build)
    study = load_study(THERMAL_FILE)
    assert isinstance(study.geometry, ThermalBallGeometry)
    with pytest.raises(TypeError) as api_exc:
        study.estimate(comm)
    _say(f"[WF-13 s9] API: TypeError: {api_exc.value}")
    assert str(api_exc.value) == THERMAL_MESSAGE

    # (3) the CLI
    out_dir = tmp_path / "out"
    capsys.readouterr()
    rc = main(["run", "--estimate", str(THERMAL_FILE), "--output-dir", str(out_dir)])
    err = capsys.readouterr().err
    lines = [ln for ln in err.splitlines() if ln.startswith("femem: error:")]

    # (4) anchor
    rcs = comm.allgather(rc)
    build_counts = comm.allgather(len(calls))
    _say(f"[WF-13 s9] CLI rc per rank {rcs}; rank-0 stderr {lines!r}; "
         f"build calls per rank {build_counts}")
    assert rcs == [2] * comm.size
    if comm.rank == 0:
        assert lines == [f"femem: error: {api_exc.value}"]
    else:
        assert lines == []
    assert not out_dir.exists()
    assert build_counts == [0] * comm.size


def test_wf13_cli_estimate_refuses_ranks(capsys, tmp_path, monkeypatch):
    """(X) ``--ranks 2 --estimate``: rc 2 naming both flags; the stub launcher
    (writes a marker) is never reached.  No MPI import needed, no solve."""
    from fem_em_solver.cli import main

    marker = tmp_path / "launched"
    stub = tmp_path / "stub_mpiexec"
    stub.write_text(f"#!/bin/sh\ntouch {marker}\n")
    stub.chmod(0o755)
    monkeypatch.setenv("FEMEM_MPIEXEC", str(stub))
    out_dir = tmp_path / "out"
    capsys.readouterr()
    rc = main(["run", "--ranks", "2", "--estimate", str(DEG2_FILE), "--output-dir", str(out_dir)])
    cap = capsys.readouterr()
    print(f"\n[WF-13 s4] (X) --ranks 2 --estimate: rc {rc}; stderr {cap.err.strip()!r}; "
          f"marker exists {marker.exists()}", flush=True)
    assert rc == 2
    assert "--ranks" in cap.err and "--estimate" in cap.err
    assert "femem: launching" not in cap.out
    assert not marker.exists()
    assert not out_dir.exists()
    assert os.environ["FEMEM_MPIEXEC"] == str(stub)
