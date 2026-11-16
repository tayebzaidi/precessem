"""`WF-1` step 24 (§9 item 674, §10 U99): ``load_result:`` in the YAML.

``run_study_file`` on a four-port file with ``load_result: <path>`` and a
``tune:`` section calls ``Study.load_result`` + ``Study.tune`` (U94, `WF-16`
step 14) with no mesh and no solve, and writes the same ``tune`` summary block
the solve route writes.  Fixture: step 6(b)'s degree-1 exact-rectangle PEC hole
F-small Study (``_study`` of ``test_wf16_study_tune_capacitors.py``, imported as
U94's test does), 64 MHz, ``-n 2``, κ off; its YAML twin is
``wf1/f_small_hole_pec_tune.yaml`` plus the new key.

Anchors: (I=) the YAML summary's ``tune`` equals the same Study's
``load_result`` + ``tune(1)`` by ``float.hex``; (CM) ``C_mode`` ``.9e``-equal to
``ITEM429_C_MODE_F`` and within ``LEG_D0_REPRODUCTION_BAND``; (S0) ``Study.solve``
/ ``Study._solve_ports`` patched to raise, 0 calls; (CLI = API) ``femem run`` at
``-n 2`` writes a summary equal field by field (``elapsed_s`` excepted) to the
API run's; (Ref) every refused field raises ``StudyConfigError`` naming the field
and ``load_result`` (rc 2 through ``femem run``) with 0 mesh builds.  Negative
control (asserted): the P1 <-> P2 swapped file misses ``ITEM429_C_MODE_F`` by
more than ``LEG_D0_REPRODUCTION_BAND``.  Every band and record is imported.
"""

from __future__ import annotations

import copy
import json
import time
from pathlib import Path

import numpy as np
import pytest
from mpi4py import MPI

from fem_em_solver import LumpedPortDrive, Study, StudyConfigError, cli, load_study, run_study_file
from fem_em_solver import study as study_mod
from fem_em_solver import study_yaml
from fem_em_solver.study import (
    FHumanBirdcageGeometry,
    FourPortBirdcageGeometry,
    SphereInBoxGeometry,
    ThermalBallGeometry,
    ThermalBoxGeometry,
)

from tests.complex_mode import complex_only
from tests.package.test_wf16_study_load_result import _swap12, _synthetic_s, _write
from tests.package.test_wf16_study_tune_capacitors import ITEM429_C_MODE_F, _study
from tests.package.test_wf1_study_yaml import WF1, _fresh_dir, _mapping, _rc_all, _spec_diff
from tests.validation.test_port_birdcage_four_port import LEG_D0_REPRODUCTION_BAND

OUT_DIR = Path("/workspace/logs/wf1_s24")  # gitignored
TUNE_FILE = WF1 / "f_small_hole_pec_tune.yaml"
REFUSED = ("conductor", "thermal", "study.sweep_hz", "outputs.sar", "outputs.b1_plus",
           "outputs.probes")


def _p0(*args):
    if MPI.COMM_WORLD.rank == 0:
        print(*args, flush=True)


def _hex(x):
    return float(x).hex()


def _write_yaml(path, mapping, comm=MPI.COMM_WORLD):
    import yaml

    path = Path(path)
    if comm.rank == 0:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(yaml.safe_dump(mapping, sort_keys=False))
    comm.Barrier()
    return path


def _loaded_mapping(target_name):
    m = copy.deepcopy(_mapping(TUNE_FILE))
    m["load_result"] = target_name
    return m


def _no_mesh(mp, calls):
    """Every mesh entry point and every solve entry raises (and counts)."""
    import gmsh

    def boom(*a, **k):
        calls["n"] += 1
        raise AssertionError("a load_result run reached a mesh build or a solve")

    for cls in (SphereInBoxGeometry, FHumanBirdcageGeometry, FourPortBirdcageGeometry,
                ThermalBallGeometry, ThermalBoxGeometry):
        mp.setattr(cls, "build", boom)
    mp.setattr(gmsh, "initialize", boom)
    mp.setattr(Study, "solve", boom)
    mp.setattr(Study, "_solve_ports", boom)


def _tune_block_hex_equal(block, api):
    """(I=) field by field; returns the list of differing fields."""
    bad = []
    if block["c_mode_farad"].hex() != _hex(api.c_mode_farad):
        bad.append("c_mode_farad")
    if block["off_diagonal_residual"].hex() != _hex(api.off_diagonal_residual):
        bad.append("off_diagonal_residual")
    if block["frequency_hz"].hex() != _hex(api.frequency_hz):
        bad.append("frequency_hz")
    if block["k"] != int(api.k):
        bad.append("k")
    if len(block["z_k"]) != len(api.z_k):
        bad.append("z_k.len")
    for i, ((re, im), z) in enumerate(zip(block["z_k"], api.z_k)):
        if re.hex() != _hex(complex(z).real) or im.hex() != _hex(complex(z).imag):
            bad.append(f"z_k[{i}]")
    return bad


def _summary_diff(a: dict, b: dict) -> list:
    keys = (set(a) | set(b)) - {"elapsed_s"}
    return sorted(k for k in keys if a.get(k) != b.get(k))


@complex_only
def test_wf1_s24_load_result_coding_and_refusals(capsys):
    """Parse, (I=) on a synthetic circulant, ``femem check``, and (Ref).  No mesh."""
    comm = MPI.COMM_WORLD
    tag = f"r{comm.size}"
    d = OUT_DIR / f"coding_{tag}"
    _fresh_dir(d, comm)
    fhz = _study(LumpedPortDrive()).frequency_hz
    s = _synthetic_s()
    target = _write(d / "synthetic.s4p", [(fhz, s)])
    yml = _write_yaml(d / "loaded.yaml", _loaded_mapping("synthetic.s4p"))
    calls = {"n": 0}
    with pytest.MonkeyPatch.context() as mp:
        _no_mesh(mp, calls)
        # Parse: the Study is step 11's file's; the key rides in outputs, resolved.
        st, kind, outputs, cooling = study_yaml._load(yml)
        diff = _spec_diff(st, load_study(TUNE_FILE))
        _p0(f"\n[WF-1 s24 W0] parse: kind {kind}, outputs {outputs!r}, spec diff vs "
            f"{TUNE_FILE.name} {diff}")
        assert kind == "four_port_birdcage" and cooling is None and diff == []
        assert outputs == {"tune": {"k": 1}, "load_result": target.resolve()}
        # Unchanged files carry no key.
        assert "load_result" not in study_yaml._load(TUNE_FILE)[2]

        # (I=) on the synthetic matrix: YAML route vs the API route.
        run = run_study_file(yml, output_dir=d / "out", comm=comm)
        api_st = load_study(yml)
        api_st.load_result(target)
        api = api_st.tune(k=1, comm=comm)
        bad = _tune_block_hex_equal(run.summary["tune"], api)
        with open(d / "out" / "summary.json") as fh:
            js = json.load(fh)
        bad_js = _tune_block_hex_equal(js["tune"], api)
        _p0(f"[WF-1 s24 W0] synthetic (I=): C {run.summary['tune']['c_mode_farad'].hex()} vs api "
            f"{_hex(api.c_mode_farad)}; differing {bad} (json {bad_js})")
        assert bad == [] and bad_js == []
        assert isinstance(run.study.result, study_mod.LoadedPortResult)
        assert "cells" not in js and "touchstone" not in js
        assert js["load_result"] == str(target.resolve())

        # `femem check` prints the key.
        pre = capsys.readouterr().out
        with capsys.disabled():
            if comm.rank == 0:
                print(pre, end="", flush=True)
        rc = cli.main(["check", str(yml)])
        out = capsys.readouterr().out
        if comm.rank == 0:
            _p0("[WF-1 s24 W0] femem check:\n" + out.rstrip())
            assert f"load_result: {str(target.resolve())!r}" in out.splitlines()
        assert _rc_all(rc, comm) == (0, 0)

        # (Ref) by field, through load_study and `femem run`.
        variants = {}
        base = _loaded_mapping("synthetic.s4p")
        m = copy.deepcopy(base)
        m["conductor"] = {"kind": "leontovich", "sigma_s_per_m": 5.8e7}
        variants["conductor"] = m
        m = copy.deepcopy(base)
        m["thermal"] = {"k_w_per_m_k": 0.5, "h_w_per_m2_k": 10.0, "t_inf_c": 20.0}
        variants["thermal"] = m
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
        m["outputs"] = {"probes": {"points": [[0.0, 0.0, 0.0]]}}
        variants["outputs.probes"] = m
        assert tuple(variants) == REFUSED
        m = copy.deepcopy(base)
        del m["tune"]
        variants["<no tune>"] = m
        fields = {}
        for label, mm in variants.items():
            f = _write_yaml(d / f"refused_{label.replace('.', '_').strip('<>').replace(' ', '_')}.yaml", mm)
            with pytest.raises(StudyConfigError) as ei:
                load_study(f)
            o = d / f"out_{f.stem}"
            rc_run = _rc_all(cli.main(["run", str(f), "--output-dir", str(o)]), comm)
            exists = comm.allreduce(int(o.exists()), op=MPI.MAX)
            fields[label] = ei.value.field
            _p0(f"  (Ref) {label}: field {ei.value.field!r}, run rc {rc_run}, out exists {exists}; {ei.value}")
            assert type(ei.value) is StudyConfigError
            assert "load_result" in str(ei.value)
            assert rc_run == (2, 2) and exists == 0
        assert {k: fields[k] for k in REFUSED} == {k: k for k in REFUSED}, fields
        assert fields["<no tune>"] == "load_result"

        # Non-four-port geometries keep their refusal (parse, no mesh).  The F-human left this
        # list with `WF-1` step 26 (§9 item 714): the 10-31 weekly admitted it.
        for name in ("sphere_mat4.yaml", "thermal_ball_wall.yaml"):
            mm = copy.deepcopy(_mapping(WF1 / name))
            mm["load_result"] = "synthetic.s4p"
            f = _write_yaml(d / f"refused_{Path(name).stem}.yaml", mm)
            with pytest.raises(StudyConfigError) as ei:
                load_study(f)
            _p0(f"  (Ref) {name} + load_result: field {ei.value.field!r}; {ei.value}")
            assert type(ei.value) is StudyConfigError

        # --backend palace refuses by name (before any mesh).
        o = d / "out_palace"
        rc_p = _rc_all(cli.main(["run", str(yml), "--backend", "palace", "--output-dir", str(o)]), comm)
        with pytest.raises(StudyConfigError) as ei:
            study_yaml._run_study_file_palace(yml, output_dir=o, comm=comm)
        _p0(f"  (Ref) --backend palace: rc {rc_p}; field {ei.value.field!r}; {ei.value}")
        assert rc_p == (2, 2) and "load_result" in str(ei.value)
    builds = comm.allreduce(calls["n"], op=MPI.SUM)
    _p0(f"[WF-1 s24 W0] mesh builds / solves reached: {builds}")
    assert builds == 0


@pytest.fixture(scope="module")
def s24():
    comm = MPI.COMM_WORLD
    out = {}
    d = OUT_DIR
    if comm.rank == 0:
        d.mkdir(parents=True, exist_ok=True)
    comm.Barrier()
    s0 = _study(LumpedPortDrive())
    t0 = time.perf_counter()
    res0 = s0.solve(comm)
    out["t_solve"] = time.perf_counter() - t0
    out["path"] = Path(res0.write_touchstone(d / f"s6b_pec_d1_n{comm.size}.s4p"))
    out["S"] = np.asarray(res0.s_matrix, dtype=np.complex128).copy()
    out["yaml"] = _write_yaml(d / f"s6b_pec_d1_n{comm.size}_load.yaml",
                              _loaded_mapping(out["path"].name))
    sw = _write(d / f"swap12_n{comm.size}.s4p", [(s0.frequency_hz, _swap12(out["S"]))])
    out["yaml_swap"] = _write_yaml(d / f"swap12_n{comm.size}_load.yaml", _loaded_mapping(sw.name))
    _p0(f"\n[WF-1 s24] solve {out['t_solve']:.1f} s, cells {res0.cells}; wrote {out['path']}")
    return out


@complex_only
def test_wf1_s24_yaml_load_result_equals_api(s24, capsys):
    """(I=), (CM), (S0), (CLI = API) and the swapped-file control."""
    comm = MPI.COMM_WORLD
    calls = {"n": 0}
    d_api = OUT_DIR / f"run_api_n{comm.size}"
    d_cli = OUT_DIR / f"run_cli_n{comm.size}"
    d_sw = OUT_DIR / f"run_swap_n{comm.size}"
    for d in (d_api, d_cli, d_sw):
        _fresh_dir(d, comm)
    with pytest.MonkeyPatch.context() as mp:
        _no_mesh(mp, calls)
        t0 = time.perf_counter()
        run = run_study_file(s24["yaml"], output_dir=d_api, comm=comm)
        t_run = time.perf_counter() - t0
        st = load_study(s24["yaml"])
        st.load_result(s24["path"])
        api = st.tune(k=1, comm=comm)
        capsys.readouterr()
        t0 = time.perf_counter()
        rc = cli.main(["run", str(s24["yaml"]), "--output-dir", str(d_cli)])
        t_cli = time.perf_counter() - t0
        capsys.readouterr()
        run_sw = run_study_file(s24["yaml_swap"], output_dir=d_sw, comm=comm)
    n_calls = comm.allreduce(calls["n"], op=MPI.SUM)

    # (S0)
    _p0(f"\n[WF-1 s24 S0] Study.solve / _solve_ports / build / gmsh.initialize patched to raise: "
        f"{n_calls} calls; run_study_file {t_run:.3f} s, femem run {t_cli:.3f} s")
    assert n_calls == 0
    assert isinstance(run.study.result, study_mod.LoadedPortResult)
    rt = comm.allreduce(float(np.max(np.abs(run.study.result.s_matrix - s24["S"]))), op=MPI.MAX)
    _p0(f"[WF-1 s24 RT] max|S_loaded - S_solved| = {rt:.3e} (printed)")

    # (I=)
    bad = _tune_block_hex_equal(run.summary["tune"], api)
    with open(d_api / "summary.json") as fh:
        js_api = json.load(fh)
    bad_js = _tune_block_hex_equal(js_api["tune"], api)
    _p0(f"[WF-1 s24 I=] summary C {run.summary['tune']['c_mode_farad'].hex()} vs api "
        f"{_hex(api.c_mode_farad)}; residual {run.summary['tune']['off_diagonal_residual'].hex()} vs "
        f"{_hex(api.off_diagonal_residual)}; z_k {len(api.z_k)}; differing {bad} (json {bad_js})")
    assert bad == [] and bad_js == []

    # (CM)
    c = run.summary["tune"]["c_mode_farad"]
    rel = abs(c / ITEM429_C_MODE_F - 1.0)
    _p0(f"[WF-1 s24 CM] C_mode {c:.12e} F; '{c:.9e}' vs item 429 '{ITEM429_C_MODE_F:.9e}'; rel "
        f"{rel:.3e} (LEG_D0_REPRODUCTION_BAND {LEG_D0_REPRODUCTION_BAND:.0e})")
    assert f"{c:.9e}" == f"{ITEM429_C_MODE_F:.9e}"
    assert rel <= LEG_D0_REPRODUCTION_BAND

    # (CLI = API)
    rcs = _rc_all(rc, comm)
    with open(d_cli / "summary.json") as fh:
        js_cli = json.load(fh)
    diff = _summary_diff(js_api, js_cli)
    _p0(f"[WF-1 s24 CLI=API] femem run rc {rcs}; summary keys {sorted(js_cli)}; "
        f"fields differing (elapsed_s excepted) {diff}")
    assert rcs == (0, 0)
    assert diff == []
    assert js_cli["ranks"] == comm.size == 2

    # Negative control (asserted): P1 <-> P2 swapped file.
    cw = run_sw.summary["tune"]["c_mode_farad"]
    miss = abs(cw / ITEM429_C_MODE_F - 1.0)
    _p0(f"[WF-1 s24 Sw] swapped file C {cw:.12e} F, miss {miss:.5f} vs band "
        f"{LEG_D0_REPRODUCTION_BAND:.0e} (U94 measured 0.197)")
    assert miss > LEG_D0_REPRODUCTION_BAND
