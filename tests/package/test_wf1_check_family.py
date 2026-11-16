"""`WF-1` step 9 (§9 item 401): ``femem check`` prints the parsed geometry kind.

No mesh, no solve; run at ``-n 2``.

(F)   every tracked non-malformed YAML under ``tests/package/wf1/``: ``cli.main(["check", f])``
      returns 0 on every rank and its ``family:`` line equals ``study_yaml._load(f)[1]``, the
      kind ``run_study_file`` records;
(F-h) the known-issues letter: F-human prints ``f_human_birdcage``, ``sphere_mat4.yaml``
      ``sphere_in_box``, ``f_small_four_port.yaml`` ``four_port_birdcage``;
(U)   on those last two files the other five printed lines equal the pinned pre-change
      ``cli.py`` (``b01dc790``, copied to ``logs/cli_b01dc790.py`` by a host ``git show``,
      loaded by path inside the ``fem_em_solver`` package);
(N)   negative control: the pre-change ``check`` on the F-human file printed
      ``family: four_port_birdcage`` (``20261012T123618Z_WF-1-s7a-W0.log:68``), which
      differs from what (F) now prints.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest
from mpi4py import MPI

from fem_em_solver import cli
from fem_em_solver.study_yaml import _load

REPO = Path(__file__).resolve().parents[2]
WF1 = REPO / "tests" / "package" / "wf1"
OLD_CLI = REPO / "logs" / "cli_b01dc790.py"
GOOD = sorted(f for f in WF1.glob("*.yaml") if not f.name.startswith("malformed_"))
PRE_CHANGE_FH_FAMILY_LINE = "family: four_port_birdcage"  # 20261012T123618Z_WF-1-s7a-W0.log:68


def _rc_all(rc, comm):
    # copied from tests/package/test_wf1_study_yaml.py:169 (a test module's helper)
    return comm.allreduce(int(rc), op=MPI.MIN), comm.allreduce(int(rc), op=MPI.MAX)


def _p0(msg, capsys):
    # past capsys, which would otherwise swallow the record into the next readouterr()
    if MPI.COMM_WORLD.rank == 0:
        with capsys.disabled():
            print(msg, flush=True)


def _check_lines(main, f, capsys):
    capsys.readouterr()
    rc = main(["check", str(f)])
    out = capsys.readouterr().out.splitlines()
    return rc, out


@pytest.fixture(scope="module")
def old_cli():
    if not OLD_CLI.exists():
        pytest.fail(f"{OLD_CLI} missing: run `git show b01dc790:src/fem_em_solver/cli.py > {OLD_CLI}`")
    spec = importlib.util.spec_from_file_location("fem_em_solver._cli_b01dc790", OLD_CLI)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_wf1_check_prints_the_parsed_kind(capsys):
    comm = MPI.COMM_WORLD
    _p0(f"\n[WF-1 s9] (F) {len(GOOD)} files: {[f.name for f in GOOD]}", capsys)
    assert len(GOOD) >= 10
    printed = {}
    for f in GOOD:
        rc, out = _check_lines(cli.main, f, capsys)
        assert _rc_all(rc, comm) == (0, 0), f.name
        kind = _load(f)[1]
        if comm.rank == 0:
            fam = [ln for ln in out if ln.startswith("family: ")]
            assert len(fam) == 1, (f.name, out)
            printed[f.name] = fam[0]
            _p0(f"  (F) {f.name}: printed {fam[0]!r}, _load kind {kind!r}", capsys)
            assert fam[0] == f"family: {kind}"
    if comm.rank == 0:
        # (F-h)
        assert printed["f_human_tuned_copper.yaml"] == "family: f_human_birdcage"
        assert printed["sphere_mat4.yaml"] == "family: sphere_in_box"
        assert printed["f_small_four_port.yaml"] == "family: four_port_birdcage"
        # (N)
        _p0(f"  (N) pre-change F-human line {PRE_CHANGE_FH_FAMILY_LINE!r} vs now "
            f"{printed['f_human_tuned_copper.yaml']!r}", capsys)
        assert printed["f_human_tuned_copper.yaml"] != PRE_CHANGE_FH_FAMILY_LINE


def test_wf1_check_other_lines_unchanged(old_cli, capsys):
    comm = MPI.COMM_WORLD
    for name in ("sphere_mat4.yaml", "f_small_four_port.yaml"):
        rc_new, new = _check_lines(cli.main, WF1 / name, capsys)
        rc_old, old = _check_lines(old_cli.main, WF1 / name, capsys)
        assert _rc_all(rc_new, comm) == (0, 0) and _rc_all(rc_old, comm) == (0, 0)
        if comm.rank == 0:
            _p0(f"\n  (U) {name}: new {new}\n      old {old}", capsys)
            assert len(new) == len(old) == 6
            assert new[0] == old[0]  # both families agree on these two files
            assert new[1:] == old[1:]
