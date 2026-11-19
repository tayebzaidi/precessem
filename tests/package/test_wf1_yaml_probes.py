"""`WF-1` step 6 (§9 item 345, chain U34): ``outputs.probes`` in the YAML.

Schema v0's ``outputs`` gains ``probes: [{centre, radius}, …]`` on a
``sphere_in_box`` with a ``thermal:`` section; ``run_study_file`` drives
``Study.probe(probes, cooling)`` after the temperature block and writes the
readings into the JSON summary (brief: docs/planning/chunks/WF-1.md, "§9 item
345").  The probe values come from the thermal solve inside ``Study.probe`` (a
second, cheap thermal solve beside the temperature block's; one EM solve) —
disclosed, and the two θ's are compared (printed).

Asserted, on `MAT-12`'s fixture at EM degree 1, P1 thermal, ``-n 2``:
(K)   the parsed probes equal ``[WHOLE] + PROBES`` by ``float.hex``; the cooling
      equals ``RobinCooling(K, H, T_INF, source_scale=G)`` (diff ``[]``); the
      Study spec diff vs ``_study(degree=1)`` is ``[]``;
(A0)  ``run.probes.thermal.submesh`` global cells == ``N_SUB``;
(A)   P8's whole-ball identity through YAML: ``values_k[0]`` vs ∫θ dx / V_sub
      and the whole-ball volume vs V_sub, each <= ``B_BAND``;
(U′)  each summary value == ``ball_average`` on ``run.probes.theta`` exactly;
      JSON == memory exactly;
(IT)  vs an independently solved Python Study's ``probe``, per probe relative
      <= ``L0_BAND``;
(X)   two malformed files raise ``StudyConfigError`` naming ``outputs.probes`` /
      ``outputs.probes[1].radius``; ``cli.main`` run / check rc 2 on every rank.
Negative control (asserted): (c1) c₂ displaced 0.1 mm outward (``C2_OUT``)
misses the summary's c₂ value by > ``WIDTH_BAND``.

Every band, constant and helper is imported, none restated.
"""

from __future__ import annotations

import json
import time

import numpy as np
import pytest
import ufl
from dolfinx import default_scalar_type, fem
from mpi4py import MPI

from fem_em_solver import RobinCooling, StudyConfigError, cli, load_study, run_study_file
from fem_em_solver.thermal.probe import ball_average
from tests.complex_mode import complex_only
from tests.package.test_wf1_study_yaml import LOGS, WF1, _fresh_dir, _rc_all, _spec_diff
from tests.package.test_wf1_yaml_thermal import _hex_diff
from tests.package.test_wf15_study_lossy_sphere import _study
from tests.package.test_wf15_study_probe import C2_OUT, PROBES, WHOLE
from tests.package.test_wf15_study_temperature import _theta_owned_rel
from tests.validation.test_mat5_em_coupled_sphere import L0_BAND
from tests.validation.test_mat10_sar_to_heat_sphere import B_BAND, T_INF, G, H, K, _allsum
from tests.validation.test_mat19_virtual_temperature_probe import N_SUB, WIDTH_BAND

PROBES_FILE = WF1 / "sphere_mat10_probes.yaml"
MALFORMED_PROBES = {
    "malformed_probes_without_thermal.yaml": "outputs.probes",
    "malformed_probe_radius_zero.yaml": "outputs.probes[1].radius",
}


def _p0(*args):
    if MPI.COMM_WORLD.rank == 0:
        print(*args, flush=True)


def _probe_hex(p):
    return [float(x).hex() for x in p.centre], float(p.radius).hex()


def test_wf1_probes_yaml_parse_and_refusals():
    """(K) and (X): no mesh, no solve."""
    comm = MPI.COMM_WORLD
    from fem_em_solver.study_yaml import _load

    study, kind, outputs, cooling = _load(PROBES_FILE)
    parsed = outputs["probes"]
    ref = [WHOLE] + PROBES
    hex_equal = [_probe_hex(a) == _probe_hex(b) for a, b in zip(parsed, ref)]
    d_cool = _hex_diff(cooling, RobinCooling(K, H, T_INF, source_scale=G))
    d_spec = _spec_diff(study, _study(degree=1))
    _p0(
        f"\n[WF-1 s6] (K) {len(parsed)} probes parsed; hex equal per probe {hex_equal}; "
        f"cooling diff {d_cool}; spec diff vs _study(degree=1) {d_spec}"
    )
    for p in parsed:
        _p0(f"  (K) {_probe_hex(p)}")
    assert kind == "sphere_in_box"
    assert isinstance(parsed, tuple) and len(parsed) == len(ref) == 4
    assert all(hex_equal), hex_equal
    assert d_cool == [], d_cool
    assert d_spec == [], d_spec

    out = LOGS / "malformed_probes"
    fields = []
    for name, expected in MALFORMED_PROBES.items():
        path = WF1 / name
        with pytest.raises(StudyConfigError) as ei:
            load_study(path)
        exc = ei.value
        _p0(f"  (X) {name}: field {exc.field!r} (expected {expected!r}); {exc}")
        assert type(exc) is StudyConfigError
        assert exc.field == expected
        assert exc.source.endswith(name)
        fields.append(exc.field)
        _fresh_dir(out, comm)
        rc_run = cli.main(["run", str(path), "--output-dir", str(out)])
        rc_check = cli.main(["check", str(path)])
        run_all, check_all = _rc_all(rc_run, comm), _rc_all(rc_check, comm)
        exists = comm.allreduce(bool(out.exists()), op=MPI.LOR)
        _p0(f"  (X) {name}: run rc (min, max) {run_all}, check rc {check_all}, out exists {exists}")
        assert run_all == (2, 2)
        assert check_all == (2, 2)
        assert not exists
    assert fields == list(MALFORMED_PROBES.values()), fields


@pytest.fixture(scope="module")
def case():
    comm = MPI.COMM_WORLD
    d = LOGS / PROBES_FILE.stem
    _fresh_dir(d, comm)
    t0 = time.perf_counter()
    run = run_study_file(PROBES_FILE, output_dir=d, comm=comm)
    t_run = time.perf_counter() - t0
    with open(d / "summary.json") as fh:
        summary = json.load(fh)
    _p0(
        f"\n[WF-1 s6] run_study_file {t_run:.1f} s (summary elapsed {summary['elapsed_s']:.1f} s); "
        f"summary['probes'] {summary['probes']}"
    )
    # Disclosure: the temperature block's θ vs the probe's θ (two thermal solves).
    two_theta, _ = _theta_owned_rel(run.thermal.temperature, run.probes.thermal.temperature, comm)
    _p0(
        f"  temperature-block θ vs probe θ owned max rel {two_theta:.3e} (printed; one EM solve: "
        f"{run.study.result is run.result})"
    )
    return {"run": run, "summary": summary, "t_run": t_run}


@complex_only
def test_wf1_probes_yaml_whole_ball_identity(case):
    """(A0), (A): P8's whole-ball identity on the YAML run's own θ."""
    comm = MPI.COMM_WORLD
    run, summary = case["run"], case["summary"]
    pr = run.probes
    assert pr is not None
    assert run.study.result is run.result, "(Once) a second EM solve happened"
    sub, theta = pr.thermal.submesh, pr.theta
    cells = int(sub.topology.index_map(3).size_global)
    v_sub = float(
        np.real(_allsum(comm, fem.Constant(sub, default_scalar_type(1.0)) * ufl.dx(domain=sub)))
    )
    mean = float(np.real(_allsum(comm, theta * ufl.dx(domain=sub)))) / v_sub
    a_val = abs(summary["probes"]["values_k"][0] - mean) / abs(mean)
    a_vol = abs(summary["probes"]["volumes_m3"][0] / v_sub - 1.0)
    _p0(
        f"\n[WF-1 s6] (A0) sub cells {cells} (N_SUB {N_SUB}); (A) |P - int/V|/|int/V| = "
        f"{a_val:.3e}, |V_h/V_sub - 1| = {a_vol:.3e} (band {B_BAND:.0e}); V_sub {v_sub:.10e}; "
        f"whole-ball mean {mean:.10e} K"
    )
    assert cells == N_SUB, f"(A0) {cells}"
    assert a_val <= B_BAND, f"(A) value {a_val:.3e}"
    assert a_vol <= B_BAND, f"(A) volume {a_vol:.3e}"


@complex_only
def test_wf1_probes_yaml_equals_api(case):
    """(U′): summary == ``ball_average`` on the run's θ exactly; JSON == memory."""
    run, summary = case["run"], case["summary"]
    pr = run.probes
    direct = [float(ball_average(pr.theta, p)[0]) for p in pr.probes]
    mem = run.summary["probes"]
    equal = summary["probes"]["values_k"] == direct
    json_equal = summary["probes"] == mem
    _p0(
        f"\n[WF-1 s6] (U′) summary values {summary['probes']['values_k']!r}; direct "
        f"{direct!r}; equal {equal}; JSON == memory {json_equal}"
    )
    assert equal, (summary["probes"]["values_k"], direct)
    assert json_equal
    assert mem["values_k"] == [float(v) for v in pr.values]


@complex_only
def test_wf1_probes_yaml_vs_python_study(case):
    """(IT): an independently solved Python Study's ``probe`` within ``L0_BAND``."""
    comm = MPI.COMM_WORLD
    summary = case["summary"]
    t0 = time.perf_counter()
    py = _study(degree=1)
    py.solve(comm)
    pr_py = py.probe([WHOLE] + PROBES, RobinCooling(K, H, T_INF, source_scale=G))
    t_py = time.perf_counter() - t0
    yv = np.asarray(summary["probes"]["values_k"])
    rel = np.abs(np.asarray(pr_py.values) - yv) / np.abs(yv)
    _p0(
        f"\n[WF-1 s6] (IT) per probe rel {rel.tolist()!r}; max {float(np.max(rel)):.3e} "
        f"(band {L0_BAND:.0e}); Python Study + probe {t_py:.1f} s"
    )
    assert float(np.max(rel)) <= L0_BAND, f"(IT) {float(np.max(rel)):.3e}"


@complex_only
def test_wf1_probes_yaml_negative_control(case):
    """(c1): c₂ displaced 0.1 mm outward misses the summary's c₂ by > ``WIDTH_BAND``."""
    run, summary = case["run"], case["summary"]
    v_out = float(ball_average(run.probes.theta, C2_OUT)[0])
    c2 = summary["probes"]["values_k"][3]
    c1 = abs(v_out - c2) / abs(c2)
    _p0(
        f"\n[WF-1 s6] (c1) c2 +0.1 mm: |d|/theta = {c1:.4e} (item 326 read 7.5333e-03; "
        f"band {WIDTH_BAND:.0e})"
    )
    assert c1 > WIDTH_BAND, f"(c1) {c1:.3e} did not miss {WIDTH_BAND:.0e}"
