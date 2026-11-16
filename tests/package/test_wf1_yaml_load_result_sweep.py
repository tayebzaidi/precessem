"""`WF-1` step 27 (§9 item 730, §10 U111): YAML ``tune.at_hz`` on a multi-frequency
``load_result:``.

The file is U71's sweep output, ``logs/wf16/wf16_sweep_f_small_p1_n2.s4p`` (gitignored;
absent ⇒ fail, never re-solve).  For each frequency of the file: (I=) the summary and its
JSON are ``float.hex``-equal to ``Study.load_result`` + ``tune(k=1, at_hz=f)``; (C=)
``C_mode`` at ``REDUCTION_BAND`` against item 722's record column; (S0) 0 solves; (CLI =
API) rc (0, 0), summary diff [].  Negative control (asserted): the 128 MHz summary read
against the 64 MHz ``C_mode`` misses by more than ``LEG_D0_REPRODUCTION_BAND``, and the
ceiling from the file's own columns equals it.

(Ref): ``tune.at_hz`` outside the span, NaN, a non-number, on a one-frequency file and
with no ``load_result:`` each refuse by that field, rc (2, 2), no output directory.  The
multi-frequency F-clinical refusal, unexercised by name at item 722, gets one assertion.

Not claimed: interpolation between frequencies, a multi-frequency F-human path, a
resonance.  ``summary["frequency_hz"]`` stays the Study's; the selected frequency is
``summary["tune"]["frequency_hz"]``.
"""

from __future__ import annotations

import copy
import json
from pathlib import Path

import numpy as np
import pytest
from mpi4py import MPI

from fem_em_solver import LumpedPortDrive, Study, StudyConfigError, cli, load_study, run_study_file
from fem_em_solver.ports.touchstone import load_touchstone

from tests.complex_mode import complex_only
from tests.package.test_wf13_estimate_f_clinical import FC
from tests.package.test_wf16_load_result_sweep import FREQS, SWEEP, _zk
from tests.package.test_wf16_study_load_result import _synthetic_s, _write
from tests.package.test_wf1_study_yaml import _fresh_dir, _rc_all
from tests.package.test_wf1_yaml_load_result import (
    _hex,
    _loaded_mapping,
    _no_mesh,
    _summary_diff,
    _tune_block_hex_equal,
    _write_yaml,
)
from tests.validation.test_port_birdcage_four_port import LEG_D0_REPRODUCTION_BAND
from tests.validation.test_port_lumped_rlc_termination import REDUCTION_BAND

OUT_DIR = Path("/workspace/logs/wf1_s27")  # gitignored
#: Item 722's record column on U71's file, F, in ``FREQS`` order.
ITEM722_C_MODE_F = (2.472675544637e-09, 6.166035421775e-11, 1.513380716932e-11)


def _p0(*args):
    if MPI.COMM_WORLD.rank == 0:
        print(*args, flush=True)


def _doc(target, tune):
    doc = copy.deepcopy(_loaded_mapping(str(target)))
    doc["tune"] = tune
    return doc


@complex_only
def test_wf1_s27_yaml_at_hz_equals_api(capsys):
    """(I=), (C=), (S0), (CLI = API) and the control, on U71's file."""
    comm = MPI.COMM_WORLD
    assert SWEEP.is_file(), f"{SWEEP} is absent: this test never re-solves"
    freqs, _stack, _ids = load_touchstone(SWEEP)
    assert tuple(float(f) for f in freqs) == FREQS, freqs
    root = OUT_DIR / f"run_n{comm.size}"
    _fresh_dir(root, comm)
    calls = {"n": 0}
    summaries = {}
    apis = {}
    results = {}
    with pytest.MonkeyPatch.context() as mp:
        _no_mesh(mp, calls)
        for i, f in enumerate(FREQS):
            yml = _write_yaml(root / f"f{i}.yaml", _doc(SWEEP, {"k": 1, "at_hz": f}))
            run = run_study_file(yml, output_dir=root / f"api{i}", comm=comm)
            st = load_study(yml)
            st.load_result(SWEEP)
            apis[f] = st.tune(k=1, comm=comm, at_hz=f)
            capsys.readouterr()
            rc = cli.main(["run", str(yml), "--output-dir", str(root / f"cli{i}")])
            capsys.readouterr()
            summaries[f] = run.summary
            results[f] = (i, rc)
    n_calls = comm.allreduce(calls["n"], op=MPI.SUM)
    _p0(f"\n[WF-1 s27 S0] solve, _solve_ports, build and gmsh.initialize patched: {n_calls} calls")

    bad = {}
    for f in FREQS:
        i, rc = results[f]
        api, block = apis[f], summaries[f]["tune"]
        with open(root / f"api{i}" / "summary.json") as fh:
            js_api = json.load(fh)
        with open(root / f"cli{i}" / "summary.json") as fh:
            js_cli = json.load(fh)
        rel = abs(block["c_mode_farad"] / ITEM722_C_MODE_F[i] - 1.0)
        rcs = _rc_all(rc, comm)
        diff = _summary_diff(js_api, js_cli)
        hex_bad = _tune_block_hex_equal(block, api) + _tune_block_hex_equal(js_api["tune"], api)
        _p0(f"[WF-1 s27 I=] at_hz {f:.3e}: selected {block['frequency_hz']:.3e}, C "
            f"{block['c_mode_farad'].hex()} vs api {_hex(api.c_mode_farad)}, differing {hex_bad}")
        _p0(f"[WF-1 s27 C=] {f:.3e} Hz: {block['c_mode_farad']:.12e} vs record "
            f"{ITEM722_C_MODE_F[i]:.12e}, rel {rel:.3e} (REDUCTION_BAND {REDUCTION_BAND:.0e})")
        _p0(f"[WF-1 s27 CLI=API] {f:.3e} Hz: rc {rcs}, fields differing {diff}; summary "
            f"frequency_hz {js_api['frequency_hz']:.3e} (the Study's)")
        bad[f] = (hex_bad, rel <= REDUCTION_BAND, rcs, diff, block["frequency_hz"] == f)
    assert n_calls == 0
    for f, (hex_bad, in_band, rcs, diff, selected) in bad.items():
        assert hex_bad == [] and in_band and rcs == (0, 0) and diff == [] and selected, f

    # Control: the 128 MHz summary against the 64 MHz C_mode, and the ceiling from the file.
    c64 = summaries[FREQS[1]]["tune"]["c_mode_farad"]
    c128 = summaries[FREQS[2]]["tune"]["c_mode_farad"]
    miss = abs(c128 / c64 - 1.0)
    ceiling = abs(1.0 - (FREQS[1] * _zk(_stack[1])[1].imag) / (FREQS[2] * _zk(_stack[2])[1].imag))
    _p0(f"[WF-1 s27 ctl] 128 MHz against 64 MHz: miss {miss:.6f}, ceiling {ceiling:.6f} "
        f"(asserted > {LEG_D0_REPRODUCTION_BAND:.0e})")
    assert miss > LEG_D0_REPRODUCTION_BAND
    assert abs(miss - ceiling) <= REDUCTION_BAND


@complex_only
def test_wf1_s27_refusals(capsys):
    """(Ref) by field ``tune.at_hz``, rc (2, 2), no output directory, no mesh."""
    comm = MPI.COMM_WORLD
    root = OUT_DIR / f"ref_n{comm.size}"
    _fresh_dir(root, comm)
    s = _synthetic_s()
    multi = _write(root / "multi.s4p", [(FREQS[0], s), (FREQS[1], 0.9 * s), (FREQS[2], 0.8 * s)])
    one = _write(root / "one.s4p", [(FREQS[1], s)])
    no_load = _doc(multi, {"k": 1, "at_hz": FREQS[1]})
    del no_load["load_result"]
    cases = {
        "below the span": _doc(multi, {"k": 1, "at_hz": FREQS[0] * 0.5}),
        "above the span": _doc(multi, {"k": 1, "at_hz": FREQS[2] * 1.5}),
        "NaN": _doc(multi, {"k": 1, "at_hz": float("nan")}),
        "non-number": _doc(multi, {"k": 1, "at_hz": "64 MHz"}),
        "one-frequency file": _doc(one, {"k": 1, "at_hz": FREQS[1]}),
        "no load_result": no_load,
    }
    calls = {"n": 0}
    readings = {}
    reasons = {}
    with pytest.MonkeyPatch.context() as mp:
        _no_mesh(mp, calls)
        for i, (name, doc) in enumerate(cases.items()):
            yml = _write_yaml(root / f"case{i}.yaml", doc)
            out = root / f"out{i}"
            with pytest.raises(StudyConfigError) as err:
                run_study_file(yml, output_dir=out, comm=comm)
            capsys.readouterr()
            rc = cli.main(["run", str(yml), "--output-dir", str(out)])
            capsys.readouterr()
            readings[name] = (err.value.field, _rc_all(rc, comm), out.exists())
            reasons[name] = err.value.reason
    assert comm.allreduce(calls["n"], op=MPI.SUM) == 0
    for name, reading in readings.items():  # printed here, past capsys
        _p0(f"\n[WF-1 s27 Ref] {name}: {reading[0]}: {reasons[name]}; rc {reading[1]}, "
            f"output dir {reading[2]}")
    for name, reading in readings.items():
        assert reading == ("tune.at_hz", (2, 2), False), name


@complex_only
def test_wf1_s27_f_clinical_multi_frequency_refused():
    """Item 722's disclosure: a multi-frequency ``.s64p`` stays refused, nothing stored."""
    comm = MPI.COMM_WORLD
    st = load_study(FC)
    assert isinstance(st, Study) and isinstance(st.excitation, LumpedPortDrive)
    root = OUT_DIR / f"fc_n{comm.size}"
    _fresh_dir(root, comm)
    s64 = 0.1 * np.eye(64, dtype=np.complex128)
    path = _write(root / "ring.s64p", [(st.frequency_hz, s64), (st.frequency_hz * 2.0, s64)])
    with pytest.raises(ValueError, match="exactly one frequency") as err:
        st.load_result(path)
    assert st.result is None
    _p0(f"\n[WF-1 s27 FC] multi-frequency .s64p: {err.value}")


@complex_only
def test_ops76_no_at_hz_refused_by_field(capsys):
    """`OPS-76` (§9 item 754): (Rc), (S0), (O), (API=); a one-frequency file still tunes."""
    comm = MPI.COMM_WORLD
    root = OUT_DIR / f"ops76_n{comm.size}"
    _fresh_dir(root, comm)
    s = _synthetic_s()
    multi = _write(root / "multi.s4p", [(FREQS[0], s), (FREQS[1], 0.9 * s), (FREQS[2], 0.8 * s)])
    yml = _write_yaml(root / "no_at_hz.yaml", _doc(multi, {"k": 1}))
    one = _write(root / "one.s4p", [(load_study(yml).frequency_hz, s)])
    one_yml = _write_yaml(root / "one.yaml", _doc(one, {"k": 1}))
    out = root / "out"
    calls = {"n": 0}
    with pytest.MonkeyPatch.context() as mp:
        _no_mesh(mp, calls)
        capsys.readouterr()
        rc = cli.main(["run", str(yml), "--output-dir", str(out)])
        stderr = capsys.readouterr().err
        st = load_study(yml)
        st.load_result(multi)
        with pytest.raises(ValueError, match="at_hz") as err:
            st.tune(1, comm=comm)
        run = run_study_file(one_yml, output_dir=root / "one_out", comm=comm)
        st_one = load_study(one_yml)
        st_one.load_result(one)
        api = st_one.tune(k=1, comm=comm)
    n_calls = comm.allreduce(calls["n"], op=MPI.SUM)
    rcs = _rc_all(rc, comm)
    named = comm.bcast("tune.at_hz" in stderr, root=0)
    wrote = comm.allreduce(out.exists(), op=MPI.LOR)
    hex_bad = _tune_block_hex_equal(run.summary["tune"], api)
    with capsys.disabled():  # printed past capsys
        _p0(f"\n[OPS-76 Rc] rc {rcs}; rank 0 stderr names tune.at_hz: {named}: {stderr.strip()}")
        _p0(f"[OPS-76 S0] solve, _solve_ports, build and gmsh.initialize patched: {n_calls} calls")
        _p0(f"[OPS-76 O] output directory written: {wrote}")
        _p0(f"[OPS-76 API=] tune(1): {type(err.value).__name__}: {err.value}")
        _p0(f"[OPS-76 ctl] one-frequency file, no at_hz: C {run.summary['tune']['c_mode_farad'].hex()}"
            f" vs api {_hex(api.c_mode_farad)}, differing {hex_bad}")
    assert rcs == (2, 2)
    assert named
    assert n_calls == 0
    assert not wrote
    assert type(err.value) is ValueError
    assert hex_bad == []
