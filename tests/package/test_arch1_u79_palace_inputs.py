"""`ARCH-1` step 1e (§10 U79): ``femem run --backend palace`` — the `WF-1` YAML
route — refuses a hole-wall birdcage (``geometry.as_hole``) and a launch wider
than one rank by field ``--backend`` with exit code 2, before any mesh and
before anything is written; a degree other than 1 or 2 is shown refused at
parse (field ``study.degree``), unchanged.

Brief: docs/planning/chunks/ARCH-1.md, "§9 item 601".  No node builds a mesh
or runs Palace; the module imports in both builds.

Anchors (asserted): (H) the hole file is refused by field, rc 2, at both
widths; (W) a 2-rank launch is refused by field, rc 2, and an ``outputs.sar``
file is still refused by ``outputs.sar`` (the allow-list precedes the width
check); (D) degree 3 is refused at parse with 0 ``_run_palace`` calls; (A')
the admissible file reaches the solve once.  Controls (rule (e)):
(c-helper), asserted, by construction — the helper no-op'd lets the hole
file reach ``Study.solve``; (c-old), printed, never asserted — the
``ValueError`` today's route raises in ``solve_palace``.
"""

from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from mpi4py import MPI

from fem_em_solver import FourPortBirdcageGeometry, cli
from fem_em_solver import study_yaml
from fem_em_solver.study_yaml import _load, _run_study_file_palace

from tests.package.test_arch1_u76_palace_outputs import _all, _field, _stub_builds, _write
from tests.package.test_wf1_study_yaml import WF1, _fresh_dir, _mapping, _rc_all

REPO = Path(__file__).resolve().parents[2]
SCRATCH = REPO / "logs" / "arch1" / "u79"  # gitignored (logs/)
HOLE = WF1 / "f_small_hole_pec_control.yaml"
FOUR_PORT = WF1 / "f_small_four_port.yaml"
HOLE_MSG = "palace: geometry.as_hole is not supported"
EM_OUTPUT_KEYS = {"s_matrix", "z_matrix", "touchstone", "sar", "b1_plus", "probes",
                  "tune", "temperature", "s_matrices"}


def _say(capsys, msg):
    if MPI.COMM_WORLD.rank == 0:
        with capsys.disabled():
            print(msg, flush=True)


def _refusal_row(f, tag, capsys, comm):
    """(field, msg, rc (min, max), rank-0 stderr, nothing written on any rank)."""
    out_py = SCRATCH / "out" / f"{tag}_py"
    out_cli = SCRATCH / "out" / f"{tag}_cli"
    _fresh_dir(out_py, comm)
    _fresh_dir(out_cli, comm)
    g, msg = _field(lambda: _run_study_file_palace(f, output_dir=out_py, comm=comm))
    capsys.readouterr()
    rc = cli.main(["run", str(f), "--output-dir", str(out_cli), "--backend", "palace"])
    err = capsys.readouterr().err
    rc_all = _rc_all(rc, comm)
    unwritten = _all(comm, not out_py.exists() and not out_cli.exists())
    return g, msg, rc_all, err, unwritten


def _fake_solve(spy):
    """U76's (A) stand-in result (test_arch1_u76_palace_outputs.py:235-240)."""

    def fake_solve(self, comm=MPI.COMM_WORLD, *, backend="dolfinx", **kw):
        if backend != "palace":
            raise AssertionError(f"spy: backend {backend!r}")
        spy["n"] += 1
        return SimpleNamespace(s_matrix=np.eye(4), cells=0, order=1,
                               case_dir=kw.get("case_dir"), sha256={}, palace_total_s=None)

    return fake_solve


# ------------------------------------------------------------------ (H)


def test_u79_hole_refused_by_field(monkeypatch, capsys):
    """(H): any width; the as_hole check precedes the width check."""
    comm = MPI.COMM_WORLD
    calls = _stub_builds(monkeypatch)
    g, msg, rc_all, err, unwritten = _refusal_row(HOLE, "hole", capsys, comm)
    builds = comm.allreduce(calls["n"], op=MPI.SUM)
    _say(capsys, f"[U79 H] n={comm.size} {HOLE.name}: field {g!r}; msg {msg!r}; rc {rc_all}; "
         f"nothing written {unwritten}; builds {builds}")
    assert g == "--backend", (g, msg)
    assert msg is not None and msg.endswith(HOLE_MSG), msg
    assert rc_all == (2, 2), rc_all
    if comm.rank == 0:
        assert "femem: error:" in err and HOLE_MSG in err, err
    assert unwritten
    assert builds == 0


# ------------------------------------------------------------------ (W)


def test_u79_width_refused_by_field(monkeypatch, capsys):
    """(W): a 2-rank launch, plus the order row (allow-list first)."""
    comm = MPI.COMM_WORLD
    if comm.size == 1:
        assert comm.size == 1
        _say(capsys, "[U79 W] skipped: comm.size == 1 (the width row needs mpiexec -n 2)")
        pytest.skip("width row needs comm.size > 1 (W0b at -n 2)")
    calls = _stub_builds(monkeypatch)
    want = f"palace: a {comm.size}-rank launch is not supported"
    g, msg, rc_all, err, unwritten = _refusal_row(FOUR_PORT, "width", capsys, comm)
    _say(capsys, f"[U79 W] n={comm.size} {FOUR_PORT.name}: field {g!r}; msg {msg!r}; "
         f"rc {rc_all}; nothing written {unwritten}")
    assert g == "--backend", (g, msg)
    assert msg is not None and f": {want}" in msg, msg
    assert rc_all == (2, 2), rc_all
    if comm.rank == 0:
        assert "femem: error:" in err and want in err, err
    assert unwritten

    base = _mapping(FOUR_PORT)
    sar = _mapping(WF1 / "f_small_sar_b1.yaml")["outputs"]["sar"]
    f_sar = _write("u79_width_sar.yaml", {**base, "outputs": {"sar": sar}}, comm)
    g2, msg2, rc2, err2, unwritten2 = _refusal_row(f_sar, "width_sar", capsys, comm)
    builds = comm.allreduce(calls["n"], op=MPI.SUM)
    _say(capsys, f"[U79 W order] n={comm.size} +outputs.sar: field {g2!r}; msg {msg2!r}; "
         f"rc {rc2}; nothing written {unwritten2}; builds {builds}")
    assert g2 == "--backend", (g2, msg2)
    assert msg2 is not None and "palace: outputs.sar is not supported" in msg2, msg2
    assert "rank launch" not in msg2, msg2
    assert rc2 == (2, 2), rc2
    assert unwritten2
    assert builds == 0


# ------------------------------------------------------------------ (D)


def test_u79_degree_refused_at_parse(monkeypatch, capsys):
    """(D): the loader's existing refusal (study_yaml.py, ``study.degree``)."""
    comm = MPI.COMM_WORLD
    calls = _stub_builds(monkeypatch)
    spy = {"n": 0}

    def spy_run_palace(*a, **kw):
        spy["n"] += 1
        raise AssertionError("_run_palace must not be reached on a parse refusal")

    monkeypatch.setattr(study_yaml, "_run_palace", spy_run_palace)
    base = _mapping(FOUR_PORT)
    f = _write("u79_degree3.yaml", {**base, "study": {**base["study"], "degree": 3}}, comm)
    g, msg = _field(lambda: _load(f))
    out = SCRATCH / "out" / "degree3_cli"
    _fresh_dir(out, comm)
    capsys.readouterr()
    rc = cli.main(["run", str(f), "--output-dir", str(out), "--backend", "palace"])
    err = capsys.readouterr().err
    rc_all = _rc_all(rc, comm)
    unwritten = _all(comm, not out.exists())
    runs = comm.allreduce(spy["n"], op=MPI.SUM)
    _say(capsys, f"[U79 D] degree 3: _load field {g!r}; msg {msg!r}; cli rc {rc_all}; "
         f"_run_palace calls {runs}; nothing written {unwritten}")
    assert g == "study.degree", (g, msg)
    assert rc_all == (2, 2), rc_all
    if comm.rank == 0:
        assert "study.degree" in err, err
    assert runs == 0
    assert unwritten
    assert calls["n"] == 0


# ------------------------------------------------------------------ (A'), controls


def test_u79_admitted_unchanged(monkeypatch, capsys):
    """(A'), (c-helper), (c-old): serial (``_run_palace`` writes on every rank)."""
    comm = MPI.COMM_WORLD
    if comm.size != 1:
        pytest.skip("serial node: the Palace route refuses a parallel comm by design")
    calls = _stub_builds(monkeypatch)
    spy = {"n": 0}
    monkeypatch.setattr("fem_em_solver.study.Study.solve", _fake_solve(spy))

    out = SCRATCH / "out" / "admit_plain"
    _fresh_dir(out, comm)
    run = _run_study_file_palace(FOUR_PORT, output_dir=out, comm=comm)
    summary = json.loads((out / "summary.json").read_text())
    keys = set(summary) & EM_OUTPUT_KEYS
    _say(capsys, f"[U79 A'] {FOUR_PORT.name}: solves {spy['n']}; EM keys {sorted(keys)}; "
         f"run.summary == file {run.summary == summary}")
    assert spy["n"] == 1
    assert keys == {"s_matrix"}

    # (c-helper), asserted: the helper no-op'd -> the hole file reaches the solve.
    monkeypatch.setattr(study_yaml, "_palace_input_refusal", lambda *a, **kw: None)
    out_c = SCRATCH / "out" / "c_helper"
    _fresh_dir(out_c, comm)
    before = spy["n"]
    _run_study_file_palace(HOLE, output_dir=out_c, comm=comm)
    _say(capsys, f"[U79 c-helper] helper no-op'd, {HOLE.name}: solves {spy['n'] - before}")
    assert spy["n"] - before == 1
    assert calls["n"] == 0

    # (c-old), printed, never asserted: the real Study.solve with build stubbed.
    monkeypatch.undo()
    monkeypatch.setattr(study_yaml, "_palace_input_refusal", lambda *a, **kw: None)
    old = {"n": 0}

    def raising_build(self, comm=MPI.COMM_WORLD):
        old["n"] += 1
        raise RuntimeError("c-old: build reached")

    monkeypatch.setattr(FourPortBirdcageGeometry, "build", raising_build)
    out_o = SCRATCH / "out" / "c_old"
    _fresh_dir(out_o, comm)
    try:
        _run_study_file_palace(HOLE, output_dir=out_o, comm=comm)
        got = "no exception"
    except Exception as exc:  # noqa: BLE001 - printed control
        got = f"{type(exc).__name__}({str(exc)!r})"
    _say(capsys, f"[U79 c-old] helper no-op'd, real solve, {HOLE.name}: {got}; "
         f"build calls {old['n']} (predicted ValueError(\"backend='palace': as_hole=True "
         f"is not supported\"), build calls 0)")
