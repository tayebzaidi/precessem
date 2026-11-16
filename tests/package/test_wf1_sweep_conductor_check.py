"""`WF-1` step 23: ``Study.solve_sweep`` refuses a Study with a conductor by name
(the YAML refusal's twin, same message), and ``femem check`` lists ``sweep_hz`` on
a sweep file — the `WF-16` sweep gains a refusal, nothing else.

Brief: docs/planning/chunks/WF-1.md, "§9 item 589" and the review rulings
2026-10-22 03:00 / 10:30.  No knob: every node builds nothing (every ``build`` is
stubbed to raise), runs in either build, at ``-n 2``.

Anchors: (Py) both copper Studies (hole from ``f_small_copper.yaml``; solid from
the F-small parts with ``LeontovichWall(COPPER_SIGMA)``) raise ``TypeError`` with
``"conductor: " + _SWEEP_CONDUCTOR_REASON``, before the list check, 0 builds;
(Y) the YAML twin: field ``study.sweep_hz``, ``reason == _SWEEP_CONDUCTOR_REASON``
(string equality), ``femem check`` / ``run`` rc (2, 2), no output directory;
(C) ``femem check`` on the sweep file prints one ``sweep_hz:`` line right after
``frequency_hz:``, its key a key of the file's ``study`` mapping, its value equal
to the file's list and to ``FREQUENCIES`` by ``float.hex``; (U) every other line
equals the pinned pre-change ``cli.py`` (``acf99b1b``; its ``check`` branch is
the same as ``6f0c8453``'s, which edited one docstring only).  Controls (rule
(e)): (c-reach) asserted — the conductor-free Study reaches ``build``;
(c-check) asserted — the pinned ``check`` prints no ``sweep_hz:`` line;
(c-old) printed, not asserted (ruling 03:00 (4)) — the pinned pre-change
``study.py`` (``6f0c8453``) lets a copper Study reach ``build`` silently.
"""

from __future__ import annotations

import ast
import subprocess
import sys
import types
from pathlib import Path

import pytest
import yaml
from mpi4py import MPI

from fem_em_solver import (
    FHumanBirdcageGeometry,
    FourPortBirdcageGeometry,
    LeontovichWall,
    SphereInBoxGeometry,
    StudyConfigError,
    ThermalBallGeometry,
    ThermalBoxGeometry,
    Study,
    cli,
    load_study,
)
from fem_em_solver.study import COPPER_SIGMA, _SWEEP_CONDUCTOR_REASON
from fem_em_solver.study_yaml import _load, _study_from_mapping

from tests.package.test_wf16_study_sweep import _port_study
from tests.package.test_wf1_study_yaml import LOGS, WF1, _fresh_dir, _mapping, _rc_all
from tests.package.test_wf1_yaml_sweep import FOUR_PORT, GOOD_SWEEP, SWEEP, _field, _with_sweep
from tests.validation.test_ans4_exact_strip_records import FREQUENCIES

REPO = Path(__file__).resolve().parents[2]
CLI_PATH = "src/fem_em_solver/cli.py"
STUDY_PATH = "src/fem_em_solver/study.py"
#: the pre-change ``cli.py`` (rule: pin a commit, never ``HEAD:``)
PRE_CHANGE_CLI_SHA = "acf99b1b"
#: the pre-change ``study.py`` (its last edit before this step; c-old, printed only)
PRE_CHANGE_STUDY_SHA = "6f0c8453"
COPPER_FILE = WF1 / "f_small_copper.yaml"
BUILD_CLASSES = (SphereInBoxGeometry, FHumanBirdcageGeometry, FourPortBirdcageGeometry,
                 ThermalBallGeometry, ThermalBoxGeometry)


class _Sentinel(Exception):
    pass


def _say(capsys, msg):
    if MPI.COMM_WORLD.rank == 0:
        with capsys.disabled():
            print(msg, flush=True)


def _git_show(sha, path):
    return subprocess.run(
        ["git", "-c", f"safe.directory={REPO}", "-C", str(REPO), "show", f"{sha}:{path}"],
        check=True, capture_output=True, text=True).stdout


def _pinned(sha, path, name, absent):
    """PORT-24's in-container loader (``tests/unit/test_port24_esr_circuit.py``), copied."""
    src = _git_show(sha, path)
    assert absent not in src, f"pinned {sha}:{path} already carries {absent!r}"
    mod = types.ModuleType(name)
    mod.__package__ = "fem_em_solver"
    sys.modules[name] = mod  # dataclasses look their module up here
    try:
        exec(compile(src, f"{sha}:{path}", "exec"), mod.__dict__)
    finally:
        sys.modules.pop(name, None)
    return mod


@pytest.fixture
def stub_build(monkeypatch):
    calls = {"n": 0}

    def stub(self, comm=MPI.COMM_WORLD):
        calls["n"] += 1
        raise _Sentinel("build reached")

    for cls in BUILD_CLASSES:
        monkeypatch.setattr(cls, "build", stub)
    return calls


def _copper_studies():
    hole = load_study(COPPER_FILE)
    base, _ = _port_study(1)
    solid = Study(base.geometry, base.material, base.excitation,
                  frequency_hz=base.frequency_hz, degree=1,
                  conductor=LeontovichWall(COPPER_SIGMA))
    return {"hole": hole, "solid": solid}


# ------------------------------------------------------------------ (Py)


def test_wf1_s23_python_sweep_refuses_conductor(stub_build, capsys):
    """(Py), (c-reach); (c-old) printed."""
    freqs = tuple(FREQUENCIES.values())
    want = "conductor: " + _SWEEP_CONDUCTOR_REASON
    for label, st in _copper_studies().items():
        assert st.conductor is not None, label
        assert isinstance(st.geometry, FourPortBirdcageGeometry), label
        n0 = stub_build["n"]
        with pytest.raises(TypeError) as exc:
            st.solve_sweep(freqs)
        got = str(exc.value)
        with pytest.raises(TypeError) as exc_bad:
            st.solve_sweep(())
        _say(capsys, f"\n[WF-1 s23] (Py) {label}: as_hole {st.geometry.as_hole}; "
             f"refusal {got!r}; equal {got == want}; builds {stub_build['n'] - n0}; "
             f"bad list () -> {type(exc_bad.value).__name__}: {str(exc_bad.value)!r}")
        assert got == want, (label, got)
        assert str(exc_bad.value) == want, (label, str(exc_bad.value))
        assert stub_build["n"] == n0, label
        assert not hasattr(st, "sweep_result"), label
        assert st.result is None, label

    # (c-reach) asserted: the same F-small Study without a conductor reaches build.
    plain, _ = _port_study(1)
    assert plain.conductor is None
    n0 = stub_build["n"]
    with pytest.raises(_Sentinel):
        plain.solve_sweep(freqs)
    reached = stub_build["n"] - n0
    _say(capsys, f"  (c-reach) conductor-free F-small: build calls {reached}")
    assert reached == 1

    # (c-old) printed, not asserted (ruling 2026-10-22 03:00 (4)): the pinned
    # pre-change solve_sweep on a copper Study of its own classes.
    old = _pinned(PRE_CHANGE_STUDY_SHA, STUDY_PATH, "fem_em_solver._study_6f0c8453",
                  "_SWEEP_CONDUCTOR_REASON")
    old_calls = {"n": 0}

    def old_stub(self, comm=MPI.COMM_WORLD):
        old_calls["n"] += 1
        raise _Sentinel("old build reached")

    old.FourPortBirdcageGeometry.build = old_stub
    ost = old.Study(old.FourPortBirdcageGeometry.f_small(), plain.material,
                    old.LumpedPortDrive(), frequency_hz=plain.frequency_hz, degree=1,
                    conductor=old.LeontovichWall(old.COPPER_SIGMA))
    try:
        ost.solve_sweep(freqs)
        outcome = "returned"
    except _Sentinel:
        outcome = "reached build (conductor ignored)"
    except Exception as e:  # printed only
        outcome = f"{type(e).__name__}: {e}"
    _say(capsys, f"  (c-old) {PRE_CHANGE_STUDY_SHA} solve_sweep on a solid copper Study: "
         f"{outcome}; build calls {old_calls['n']} (printed, not asserted)")


# ------------------------------------------------------------------ (Y)


def test_wf1_s23_yaml_twin_and_cli(stub_build, capsys):
    """(Y): the YAML refusal and the Python refusal carry the same reason."""
    comm = MPI.COMM_WORLD
    m = _with_sweep(_mapping(COPPER_FILE), GOOD_SWEEP)
    field = _field(lambda: _study_from_mapping(m))
    try:
        _study_from_mapping(m)
        reason = None
    except StudyConfigError as exc:
        reason = exc.reason
    _say(capsys, f"\n[WF-1 s23] (Y) field {field!r}; reason {reason!r}; "
         f"== _SWEEP_CONDUCTOR_REASON {reason == _SWEEP_CONDUCTOR_REASON}")
    assert field == "study.sweep_hz"
    assert reason == _SWEEP_CONDUCTOR_REASON

    d = LOGS / "s23"
    _fresh_dir(d, comm)
    f = d / "copper_sweep.yaml"
    out = d / "copper_sweep_out"
    if comm.rank == 0:
        d.mkdir(parents=True, exist_ok=True)
        f.write_text(yaml.safe_dump(m, sort_keys=False))
    comm.Barrier()
    rc_check = _rc_all(cli.main(["check", str(f)]), comm)
    rc_run = _rc_all(cli.main(["run", str(f), "--output-dir", str(out)]), comm)
    absent = comm.allreduce(int(not out.exists()), op=MPI.MIN)
    _say(capsys, f"  (Y) femem check rc {rc_check}; femem run rc {rc_run}; "
         f"output dir absent on every rank {bool(absent)}; builds {stub_build['n']}")
    assert rc_check == (2, 2)
    assert rc_run == (2, 2)
    assert absent == 1
    assert stub_build["n"] == 0


# ------------------------------------------------------------------ (C), (U)


def _check_lines(main, f, capsys):
    capsys.readouterr()
    rc = main(["check", str(f)])
    return rc, capsys.readouterr().out.splitlines()


def test_wf1_s23_check_lists_sweep_hz(stub_build, capsys):
    """(C), (U), (c-check)."""
    comm = MPI.COMM_WORLD
    old = _pinned(PRE_CHANGE_CLI_SHA, CLI_PATH, "fem_em_solver._cli_acf99b1b", "sweep_hz")

    # (C)
    rc, lines = _check_lines(cli.main, SWEEP, capsys)
    rc_c = _rc_all(rc, comm)
    file_study = yaml.safe_load(SWEEP.read_text())["study"]
    c_table = []
    if comm.rank == 0:
        sw = [i for i, ln in enumerate(lines) if ln.startswith("sweep_hz: ")]
        fq = [i for i, ln in enumerate(lines) if ln.startswith("frequency_hz: ")]
        assert len(sw) == 1, lines
        assert len(fq) == 1, lines
        assert sw[0] == fq[0] + 1, lines
        key, _, val = lines[sw[0]].partition(": ")
        parsed = ast.literal_eval(val)
        c_table = [key, val, key in file_study,
                   [float(x).hex() for x in parsed] == [float(x).hex() for x in file_study[key]],
                   [float(x).hex() for x in parsed]
                   == [float(x).hex() for x in FREQUENCIES.values()]]
        assert key == "sweep_hz"
        assert key in file_study
        assert len(parsed) == len(file_study[key])
        assert c_table[3], (parsed, file_study[key])
        assert c_table[4], (parsed, FREQUENCIES)
    assert rc_c == (0, 0)

    # (U) + (c-check)
    files = sorted(f for f in WF1.glob("*.yaml") if not f.name.startswith("malformed_"))
    files += sorted((REPO / "examples" / "mri").glob("*.yaml"))
    rows = []
    n_sweep = 0
    c_check = None
    for f in files:
        rc_new, new = _check_lines(cli.main, f, capsys)
        rc_old, prev = _check_lines(old.main, f, capsys)
        rcs = (_rc_all(rc_new, comm), _rc_all(rc_old, comm))
        assert rcs == ((0, 0), (0, 0)), (f.name, rcs)
        is_sweep = "sweep_hz" in (yaml.safe_load(f.read_text()).get("study") or {})
        n_sweep += int(is_sweep)
        if comm.rank == 0:
            if is_sweep:
                extra = [ln for ln in new if ln.startswith("sweep_hz: ")]
                assert len(extra) == 1, (f.name, new)
                i = new.index(extra[0])
                assert new[:i] + new[i + 1:] == prev, (f.name, new, prev)
                old_has = any(ln.startswith("sweep_hz") for ln in prev)
                assert not old_has, (f.name, prev)  # (c-check)
                if f == SWEEP:
                    c_check = old_has
            else:
                assert new == prev, (f.name, new, prev)
            rows.append((f.relative_to(REPO).as_posix(), len(new), len(prev), is_sweep))

    _say(capsys, f"\n[WF-1 s23] (C) femem check {SWEEP.name}: rc {rc_c}; "
         f"key/value/key-in-file/hex=file/hex=FREQUENCIES {c_table}")
    _say(capsys, f"  (U) files {len(files)}; sweep files {n_sweep} (predicted 1); "
         f"pinned cli {PRE_CHANGE_CLI_SHA}")
    for r in rows:
        _say(capsys, f"    {r[0]}: new {r[1]} lines, pinned {r[2]} lines, sweep {r[3]}")
    _say(capsys, f"  (c-check) pinned check on {SWEEP.name} prints sweep_hz: {c_check}")
    assert n_sweep == 1
    if comm.rank == 0:
        assert c_check is False
    assert FOUR_PORT in files
