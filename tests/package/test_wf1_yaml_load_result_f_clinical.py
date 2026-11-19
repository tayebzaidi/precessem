"""`WF-1` step 25 (§9 item 697, §10 U105): ``load_result:`` on an
``f_clinical_birdcage``, the YAML twin of `WF-16` step 16.

``f_clinical_birdcage`` takes no key but ``kind``, so ``load_result`` needs no
``conductor`` (U102, the F-human's design question, stays the weekly's).
``run_study_file`` reaches ``_run_load_result`` with no mesh and no solve; a file
without ``load_result:`` still returns rc 2, "xl tier only".

Anchors: (I=) the run's ``tune`` block equals the API's ``load_result`` + ``tune(1)``
by ``float.hex``; (CM-r) ``.4f`` = ``C10_C_MODE_PF_RECORD``, imported; (S0) 0 solves
and 0 ``build`` calls; (CLI = API) ``femem run`` writes a summary equal field by field
(``elapsed_s`` excepted) to the API run's; (K2) ``tune: {k: 2}`` misses the record by
more than ``LEG_D0_REPRODUCTION_BAND``; (Ref) a sweep, a field output or a Palace
backend beside ``load_result:`` is refused before any mesh.  Negative control
(asserted): the swapped file misses by more than ``LEG_D0_REPRODUCTION_BAND``.

Opt-in: ``FEM_EM_WF1_S25=1`` for the record node (the npz absent under the knob
fails, never skips, never re-solves; any other value raises).  The no-mesh node
always runs.  Scope: the YAML twin of item 696.  Not claimed: the F-human YAML
route, a page, a schema-stability statement.
"""

from __future__ import annotations

import copy
import json
import os
from pathlib import Path

import numpy as np
import pytest
from mpi4py import MPI

from fem_em_solver import Study, StudyConfigError, cli, load_study, run_study_file
from fem_em_solver import study as study_mod
from fem_em_solver import study_yaml
from fem_em_solver.ports.circuit import z_to_s
from tests.complex_mode import complex_only
from tests.package.test_wf1_study_yaml import _fresh_dir, _mapping, _rc_all
from tests.package.test_wf1_yaml_load_result import (
    _hex,
    _no_mesh,
    _summary_diff,
    _tune_block_hex_equal,
    _write_yaml,
)
from tests.package.test_wf13_estimate_f_clinical import FC
from tests.package.test_wf16_load_result_f_clinical import C10_C_MODE_PF_RECORD, N_SLOTS
from tests.package.test_wf16_load_result_f_human import _swap12_n, _write_named
from tests.validation.test_birdcage_sar_mass_averaged import EXACT_IDENTITY_RTOL
from tests.validation.test_geo38_f_clinical_c_mode import (
    C10_JSON,
    WORKSPACE,
    _block_circulant,
    _load,
    _sym_blocks,
)
from tests.validation.test_geo38_f_clinical_matrix import SAVE_PATH
from tests.validation.test_port_birdcage_four_port import LEG_D0_REPRODUCTION_BAND

S25_ENV = "FEM_EM_WF1_S25"
OUT_DIR = Path("/workspace/logs/wf1_s25")  # gitignored; tmp_path would differ per rank


def _p0(*args):
    if MPI.COMM_WORLD.rank == 0:
        print(*args, flush=True)


def _knob() -> bool:
    raw = os.environ.get(S25_ENV, "").strip()
    if raw == "":
        return False
    if raw == "1":
        return True
    raise ValueError(f"{S25_ENV}={raw!r}: only unset/empty or '1' are accepted")


def _loaded_mapping(target_name, k=1):
    m = copy.deepcopy(_mapping(FC))
    m["tune"] = {"k": k}
    m["load_result"] = target_name
    return m


@complex_only
def test_wf1_s25_coding_and_refusals(monkeypatch, capsys):
    """Knob, parse, (I=) on a synthetic 64-port, (Ref), and the refusal without
    ``load_result:``.  No mesh."""
    comm = MPI.COMM_WORLD
    lines = []
    say = lines.append  # capsys holds stdout: printed past it at the end
    for bad in ("2", "yes"):
        monkeypatch.setenv(S25_ENV, bad)
        with pytest.raises(ValueError, match=S25_ENV):
            _knob()
    monkeypatch.delenv(S25_ENV, raising=False)
    assert _knob() is False

    d = OUT_DIR / f"coding_r{comm.size}"
    _fresh_dir(d, comm)
    st0 = load_study(FC)
    c_ref = 195.0e-12
    z = _block_circulant(_sym_blocks(N_SLOTS, st0.omega, c_ref), N_SLOTS)
    s = z_to_s(z, float(st0.excitation.z0_ohm))
    names = [f"P{i}" for i in range(N_SLOTS + 1, 3 * N_SLOTS + 1)]
    target = _write_named(d / "synthetic.s64p", st0.frequency_hz, s, names)
    yml = _write_yaml(d / "loaded.yaml", _loaded_mapping("synthetic.s64p"))
    calls = {"n": 0}
    with pytest.MonkeyPatch.context() as mp:
        _no_mesh(mp, calls)
        st, kind, outputs, cooling = study_yaml._load(yml)
        say(f"\n[WF-1 s25 W0] parse: kind {kind}, outputs {outputs!r}")
        assert kind == "f_clinical_birdcage" and cooling is None
        assert outputs == {"tune": {"k": 1}, "load_result": target.resolve()}

        run = run_study_file(yml, output_dir=d / "out", comm=comm)
        api_st = load_study(yml)
        api_st.load_result(target)
        api = api_st.tune(k=1, comm=comm)
        bad = _tune_block_hex_equal(run.summary["tune"], api)
        with open(d / "out" / "summary.json") as fh:
            js = json.load(fh)
        bad_js = _tune_block_hex_equal(js["tune"], api)
        e1 = abs(run.summary["tune"]["c_mode_farad"] / c_ref - 1.0)
        say(
            f"[WF-1 s25 W0] synthetic (I=): C {run.summary['tune']['c_mode_farad'].hex()} vs api "
            f"{_hex(api.c_mode_farad)}; differing {bad} (json {bad_js}); closed form rel {e1:.3e}"
        )
        assert bad == [] and bad_js == []
        assert e1 <= EXACT_IDENTITY_RTOL
        assert isinstance(run.study.result, study_mod.LoadedPortResult)
        assert len(js["port_ids"]) == 64 and "cells" not in js and "touchstone" not in js

        # (Ref) beside load_result:, by field, through load_study and `femem run`.
        base = _loaded_mapping("synthetic.s64p")
        variants = {}
        m = copy.deepcopy(base)
        m["study"]["sweep_hz"] = [6.3e7, 6.4e7]
        variants["study.sweep_hz"] = m
        m = copy.deepcopy(base)
        m["outputs"] = {"sar": {"mass_g": 10.0, "centres": [[0.0, 0.0, 0.0]]}}
        variants["outputs.sar"] = m
        m = copy.deepcopy(base)
        m["outputs"] = {"b1_plus": {}}
        variants["outputs.b1_plus"] = m
        m = copy.deepcopy(base)
        del m["tune"]
        variants["<no tune>"] = m
        for label, mm in variants.items():
            f = _write_yaml(
                d / f"refused_{label.replace('.', '_').strip('<>').replace(' ', '_')}.yaml", mm
            )
            with pytest.raises(StudyConfigError) as ei:
                load_study(f)
            o = d / f"out_{f.stem}"
            rc_run = _rc_all(cli.main(["run", str(f), "--output-dir", str(o)]), comm)
            exists = comm.allreduce(int(o.exists()), op=MPI.MAX)
            say(
                f"  (Ref) {label}: field {ei.value.field!r}, run rc {rc_run}, out exists {exists}; "
                f"{ei.value}"
            )
            assert "load_result" in str(ei.value)
            assert rc_run == (2, 2) and exists == 0

        # --backend palace is refused before any mesh (the xl refusal, ahead of Palace's own).
        o = d / "out_palace"
        rc_p = _rc_all(
            cli.main(["run", str(yml), "--backend", "palace", "--output-dir", str(o)]), comm
        )
        with pytest.raises(StudyConfigError) as ei:
            study_yaml._run_study_file_palace(yml, output_dir=o, comm=comm)
        say(f"  (Ref) --backend palace: rc {rc_p}; field {ei.value.field!r}; {ei.value}")
        assert rc_p == (2, 2)

        # Without load_result: the file still returns rc 2, "xl tier only", tune or no tune.
        for label, mm in (("plain", _mapping(FC)), ("tune", {**_mapping(FC), "tune": {"k": 1}})):
            f = _write_yaml(d / f"no_load_{label}.yaml", mm)
            o = d / f"out_no_load_{label}"
            capsys.readouterr()
            rc = _rc_all(cli.main(["run", str(f), "--output-dir", str(o)]), comm)
            err = capsys.readouterr().err
            with pytest.raises(StudyConfigError, match="xl tier only") as ei:
                run_study_file(f, output_dir=o, comm=comm)
            say(f"  (Ref) no load_result ({label}): rc {rc}; {ei.value}")
            assert rc == (2, 2)
            if comm.rank == 0:
                assert "xl tier only" in err
    builds = comm.allreduce(calls["n"], op=MPI.SUM)
    say(f"[WF-1 s25 W0] mesh builds / solves reached: {builds}")
    with capsys.disabled():
        _p0("\n".join(lines))
    assert builds == 0


@complex_only
def test_wf1_s25_yaml_load_result_equals_api(capsys):
    """(I=), (CM-r), (S0), (CLI = API), (K2) and the swapped-file control."""
    if not _knob():
        pytest.skip(f"{S25_ENV} unset: the saved C3 matrix lives in the gitignored logs/")
    comm = MPI.COMM_WORLD
    lines = []
    say = lines.append  # capsys holds stdout: printed past it at the end
    npz = WORKSPACE / SAVE_PATH
    assert Path(npz).exists(), f"{npz} is absent: restore it, never re-solve"
    data = _load(npz)
    c10 = float(json.loads(Path(C10_JSON).read_text())["c_mode_f"])
    s64 = np.asarray(data["s64"], dtype=np.complex128)
    order = data["order"]
    d = OUT_DIR / f"record_n{comm.size}"
    _fresh_dir(d, comm)
    fhz = load_study(FC).frequency_hz
    path = _write_named(d / "c3_cu_d1.s64p", fhz, s64, order)
    swap = _write_named(d / "swap12.s64p", fhz, _swap12_n(s64), order)
    yml = _write_yaml(d / "loaded.yaml", _loaded_mapping(path.name))
    yml2 = _write_yaml(d / "loaded_k2.yaml", _loaded_mapping(path.name, k=2))
    yml_sw = _write_yaml(d / "loaded_swap.yaml", _loaded_mapping(swap.name))

    calls = {"n": 0}

    def _c_or_inf(y, out):
        try:
            c = run_study_file(y, output_dir=out, comm=comm).summary["tune"]["c_mode_farad"]
            return c, abs(c / c10 - 1.0)
        except ValueError as exc:  # a non-capacitive read is a miss too; printed
            say(f"[WF-1 s25] {y.name}: tune raised {exc} (miss = inf)")
            return float("nan"), float("inf")

    with pytest.MonkeyPatch.context() as mp:
        _no_mesh(mp, calls)
        run = run_study_file(yml, output_dir=d / "api", comm=comm)
        st = load_study(yml)
        st.load_result(path)
        api = st.tune(k=1, comm=comm)
        capsys.readouterr()
        rc = _rc_all(cli.main(["run", str(yml), "--output-dir", str(d / "cli")]), comm)
        capsys.readouterr()
        c_k2, miss_k2 = _c_or_inf(yml2, d / "k2")
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
    pf = f"{c * 1e12:.4f}"
    say(
        f"\n[WF-1 s25 S0] Study.solve / _solve_ports / build / gmsh.initialize patched to raise: "
        f"{n_calls} calls"
    )
    say(
        f"[WF-1 s25 I=] summary C {c.hex()} vs api {_hex(api.c_mode_farad)}; z_k {len(api.z_k)}; "
        f"differing {bad} (json {bad_js}); rel to C10 {abs(c / c10 - 1.0):.3e} (printed)"
    )
    say(f"[WF-1 s25 CM-r] {pf} pF vs record {C10_C_MODE_PF_RECORD}")
    say(f"[WF-1 s25 CLI] femem run rc {rc}; summary diff against the API run {diff}")
    say(
        f"[WF-1 s25 K2] tune k = 2: C {c_k2:.12e} F, miss {miss_k2:.6e} vs "
        f"LEG_D0_REPRODUCTION_BAND {LEG_D0_REPRODUCTION_BAND:.1e}"
    )
    say(
        f"[WF-1 s25 Sw] ports 1 and 2 swapped: C {c_sw:.12e} F, miss {miss_sw:.6e} vs "
        f"LEG_D0_REPRODUCTION_BAND {LEG_D0_REPRODUCTION_BAND:.1e}"
    )
    with capsys.disabled():
        _p0("\n".join(lines))
    assert n_calls == 0
    assert isinstance(run.study.result, study_mod.LoadedPortResult)
    assert isinstance(st, Study)
    assert bad == [] and bad_js == []
    assert pf == C10_C_MODE_PF_RECORD
    assert rc == (0, 0)
    assert diff == []
    assert miss_k2 > LEG_D0_REPRODUCTION_BAND
    assert miss_sw > LEG_D0_REPRODUCTION_BAND
