"""`WF-1` step 26 (§9 item 714, §10 U102): ``load_result:`` on an ``f_human_birdcage``,
the YAML twin of `WF-16` step 15.

The 10-31 weekly ruled the design question: ``conductor`` is admitted beside
``load_result:`` on this kind, because it documents what the saved column was solved
with and nothing is solved; every solve-time key is refused as U99 refuses it.  No key
is new: ``tune: {k: 1}`` reaches ``Study.tune`` with ``n_slots`` unset, which is the
geometry's ``leg_count``, 16.

Anchors: (I=) the run's ``tune`` block equals the API's ``load_result`` +
``tune(1, n_slots=16)`` by ``float.hex``; (CM) ``C_mode`` against
``WF12_STEP3_C_MODE_CU_F`` at ``REDUCTION_BAND``; (S0) 0 solves and 0 ``build`` calls;
(CLI = API) ``femem run`` writes a summary equal field by field (``elapsed_s``
excepted) to the API run's; (Ref) a sweep, ``thermal``, a field output or a Palace
backend beside a present ``conductor`` is refused by field before any mesh.  Negative
control (asserted): the swapped file misses by more than ``LEG_D0_REPRODUCTION_BAND``.

Opt-in: ``FEM_EM_WF1_S26=1`` for the record node (the file absent under the knob fails,
never skips, never re-solves; any other value raises).  The no-mesh node always runs.
Scope: the YAML twin of item 685.  Not claimed: that the saved file was solved with the
stated conductor, a sweep, the two-port drive, a schema-stability statement.
"""

from __future__ import annotations

import copy
import json
import os
import shutil
from pathlib import Path

import numpy as np
import pytest
from mpi4py import MPI

from fem_em_solver import StudyConfigError, cli, load_study, run_study_file
from fem_em_solver import study as study_mod
from fem_em_solver import study_yaml
from tests.complex_mode import complex_only
from tests.package.test_wf1_study_yaml import WF1, _fresh_dir, _mapping, _rc_all
from tests.package.test_wf1_yaml_load_result import (
    _hex,
    _no_mesh,
    _summary_diff,
    _tune_block_hex_equal,
    _write_yaml,
)
from tests.package.test_wf16_load_result_f_human import N_SLOTS, _write_named
from tests.package.test_wf16_study_tune_f_human_copper import WF12_STEP3_C_MODE_CU_F, _synthetic
from tests.validation.test_port_birdcage_four_port import LEG_D0_REPRODUCTION_BAND
from tests.validation.test_port_lumped_rlc_termination import REDUCTION_BAND

S26_ENV = "FEM_EM_WF1_S26"
FILE = WF1 / "f_human_hole_copper_tune.yaml"
SAVED = Path("/workspace/logs/wf16/s15_f_human_cu_d1_n8.s32p")  # U103's; gitignored
SWAPPED = Path("/workspace/logs/wf16/s15_swap12_n8.s32p")
OUT_DIR = Path("/workspace/logs/wf1_s26")  # gitignored; tmp_path would differ per rank


def _p0(*args):
    if MPI.COMM_WORLD.rank == 0:
        print(*args, flush=True)


def _knob() -> bool:
    raw = os.environ.get(S26_ENV, "").strip()
    if raw == "":
        return False
    if raw == "1":
        return True
    raise ValueError(f"{S26_ENV}={raw!r}: only unset/empty or '1' are accepted")


def _loaded_mapping(target_name):
    m = copy.deepcopy(_mapping(FILE))
    assert "conductor" in m and m["tune"] == {"k": 1}
    m["load_result"] = target_name
    return m


@complex_only
def test_wf1_s26_coding_and_refusals(monkeypatch, capsys):
    """Knob, parse, (I=) on a synthetic 32-port, (Ref).  No mesh."""
    comm = MPI.COMM_WORLD
    lines = []
    say = lines.append  # capsys holds stdout: printed past it at the end
    for bad in ("2", "yes"):
        monkeypatch.setenv(S26_ENV, bad)
        with pytest.raises(ValueError, match=S26_ENV):
            _knob()
    monkeypatch.delenv(S26_ENV, raising=False)
    assert _knob() is False

    d = OUT_DIR / f"coding_r{comm.size}"
    _fresh_dir(d, comm)
    st0 = load_study(FILE)
    s, names, _d = _synthetic()
    target = _write_named(d / "synthetic.s32p", st0.frequency_hz, s, names)
    yml = _write_yaml(d / "loaded.yaml", _loaded_mapping("synthetic.s32p"))
    calls = {"n": 0}
    with pytest.MonkeyPatch.context() as mp:
        _no_mesh(mp, calls)
        st, kind, outputs, _cooling = study_yaml._load(yml)
        say(
            f"\n[WF-1 s26 W0] parse: kind {kind}, outputs {outputs!r}; conductor "
            f"{type(st.conductor).__name__}"
        )
        assert kind == "f_human_birdcage"
        assert outputs == {"tune": {"k": 1}, "load_result": target.resolve()}

        run = run_study_file(yml, output_dir=d / "out", comm=comm)
        api_st = load_study(yml)
        api_st.load_result(target)
        api = api_st.tune(k=1, n_slots=N_SLOTS, comm=comm)
        bad = _tune_block_hex_equal(run.summary["tune"], api)
        with open(d / "out" / "summary.json") as fh:
            js = json.load(fh)
        say(
            f"[WF-1 s26 W0] synthetic (I=): C {run.summary['tune']['c_mode_farad'].hex()} vs api "
            f"{_hex(api.c_mode_farad)}; differing {bad}"
        )
        assert bad == []
        assert isinstance(run.study.result, study_mod.LoadedPortResult)
        assert len(js["port_ids"]) == 32 and "cells" not in js and "touchstone" not in js

        # (Ref) beside load_result: and a present conductor, by field.
        base = _loaded_mapping("synthetic.s32p")
        variants = {}
        m = copy.deepcopy(base)
        m["study"]["sweep_hz"] = [6.3e7, 6.4e7]
        variants["study.sweep_hz"] = m
        m = copy.deepcopy(base)
        m["thermal"] = {}
        variants["thermal"] = m
        for key, val in (
            ("sar", {"mass_g": 10.0, "centres": [[0.0, 0.0, 0.0]]}),
            ("b1_plus", {}),
            ("probes", {"points": [[0.0, 0.0, 0.0]]}),
        ):
            m = copy.deepcopy(base)
            m["outputs"] = {key: val}
            variants[f"outputs.{key}"] = m
        for label, mm in variants.items():
            f = _write_yaml(d / f"refused_{label.replace('.', '_')}.yaml", mm)
            with pytest.raises(StudyConfigError) as ei:
                load_study(f)
            o = d / f"out_{f.stem}"
            rc_run = _rc_all(cli.main(["run", str(f), "--output-dir", str(o)]), comm)
            exists = comm.allreduce(int(o.exists()), op=MPI.MAX)
            say(
                f"  (Ref) {label}: field {ei.value.field!r}, run rc {rc_run}, out exists {exists}; "
                f"{ei.value}"
            )
            assert ei.value.field == label
            assert "load_result" in str(ei.value)
            assert rc_run == (2, 2) and exists == 0

        o = d / "out_palace"
        rc_p = _rc_all(
            cli.main(["run", str(yml), "--backend", "palace", "--output-dir", str(o)]), comm
        )
        with pytest.raises(StudyConfigError) as ei:
            study_yaml._run_study_file_palace(yml, output_dir=o, comm=comm)
        say(f"  (Ref) --backend palace: rc {rc_p}; field {ei.value.field!r}; {ei.value}")
        assert rc_p == (2, 2)

        # The two_port_quadrature default is refused by the LumpedPortDrive check, by field.
        m = copy.deepcopy(base)
        m["excitation"] = {"kind": "two_port_quadrature", "c_in_f": 1.18e-10}
        f = _write_yaml(d / "refused_quadrature.yaml", m)
        with pytest.raises(StudyConfigError) as ei:
            load_study(f)
        say(f"  (Ref) excitation {m['excitation']}: field {ei.value.field!r}; {ei.value}")
        assert ei.value.field == "tune"  # the tune parser's refusal comes first
    builds = comm.allreduce(calls["n"], op=MPI.SUM)
    say(f"[WF-1 s26 W0] mesh builds / solves reached: {builds}")
    with capsys.disabled():
        _p0("\n".join(lines))
    assert builds == 0


@complex_only
def test_wf1_s26_yaml_load_result_equals_api(capsys):
    """(I=), (CM), (S0), (CLI = API) and the swapped-file control."""
    if not _knob():
        pytest.skip(f"{S26_ENV} unset: the saved file lives in the gitignored logs/")
    comm = MPI.COMM_WORLD
    lines = []
    say = lines.append  # capsys holds stdout: printed past it at the end
    for src in (SAVED, SWAPPED):
        assert src.exists(), f"{src} is absent: restore it, never re-solve"
    d = OUT_DIR / f"record_n{comm.size}"
    _fresh_dir(d, comm)
    if comm.rank == 0:
        d.mkdir(parents=True, exist_ok=True)
        shutil.copy(SAVED, d / SAVED.name)
        shutil.copy(SWAPPED, d / SWAPPED.name)
    comm.Barrier()
    yml = _write_yaml(d / "loaded.yaml", _loaded_mapping(SAVED.name))
    yml_sw = _write_yaml(d / "loaded_swap.yaml", _loaded_mapping(SWAPPED.name))

    calls = {"n": 0}

    def _c_or_inf(y, out):
        try:
            c = run_study_file(y, output_dir=out, comm=comm).summary["tune"]["c_mode_farad"]
            return c, abs(c / WF12_STEP3_C_MODE_CU_F - 1.0)
        except ValueError as exc:  # a non-capacitive read is a miss too; printed
            say(f"[WF-1 s26] {y.name}: tune raised {exc} (miss = inf)")
            return float("nan"), float("inf")

    with pytest.MonkeyPatch.context() as mp:
        _no_mesh(mp, calls)
        run = run_study_file(yml, output_dir=d / "api", comm=comm)
        st = load_study(yml)
        st.load_result(d / SAVED.name)
        api = st.tune(k=1, n_slots=N_SLOTS, comm=comm)
        capsys.readouterr()
        rc = _rc_all(cli.main(["run", str(yml), "--output-dir", str(d / "cli")]), comm)
        capsys.readouterr()
        c_sw, miss_sw = _c_or_inf(yml_sw, d / "swap")
    n_calls = comm.allreduce(calls["n"], op=MPI.SUM)

    bad = _tune_block_hex_equal(run.summary["tune"], api)
    with open(d / "api" / "summary.json") as fh:
        js_api = json.load(fh)
    bad_js = _tune_block_hex_equal(js_api["tune"], api)
    with open(d / "cli" / "summary.json") as fh:
        js_cli = json.load(fh)
    diff = _summary_diff(js_api, js_cli)
    c = run.summary["tune"]["c_mode_farad"]
    rel = abs(c / WF12_STEP3_C_MODE_CU_F - 1.0)
    say(
        f"\n[WF-1 s26 S0] Study.solve / _solve_ports / build / gmsh.initialize patched to raise: "
        f"{n_calls} calls"
    )
    say(
        f"[WF-1 s26 I=] summary C {c.hex()} vs api {_hex(api.c_mode_farad)}; z_k {len(api.z_k)}; "
        f"differing {bad} (json {bad_js}); ports {len(js_api['port_ids'])}"
    )
    say(
        f"[WF-1 s26 CM] C_mode vs WF12_STEP3_C_MODE_CU_F: rel {rel:.3e} (REDUCTION_BAND "
        f"{REDUCTION_BAND:.0e}; item 685 read 3.421e-13)"
    )
    say(f"[WF-1 s26 CLI] femem run rc {rc}; summary diff against the API run {diff}")
    say(
        f"[WF-1 s26 Sw] ports 1 and 2 swapped: C {c_sw:.12e} F, miss {miss_sw:.6e} vs "
        f"LEG_D0_REPRODUCTION_BAND {LEG_D0_REPRODUCTION_BAND:.1e}"
    )
    with capsys.disabled():
        _p0("\n".join(lines))
    assert n_calls == 0
    assert isinstance(run.study.result, study_mod.LoadedPortResult)
    assert bad == [] and bad_js == []
    assert np.isfinite(rel) and rel <= REDUCTION_BAND
    assert rc == (0, 0)
    assert diff == []
    assert miss_sw > LEG_D0_REPRODUCTION_BAND
