"""`WF-1` step 17: ``outputs.probes`` on a ``thermal_box`` file drives ``Study.probe``
through ``femem run`` and reproduces `WF-15` step 6's (§9 item 523, U60) box probes.

Brief: docs/planning/chunks/WF-1.md, "§9 item 539" and the review ruling
2026-10-18 18:00.  The solve legs run in the real build (the build of U60's
records, ``20261018T123416Z_WF-15-s6-W2-box.log:12``) and only with
``FEM_EM_WF1_S17=1`` (then they assert the real build and ``-n 2``, never
skip).  The parse node builds no mesh and runs in either build.

Anchors: (K) the YAML Study = the Python Study and the parsed probes / cooling
/ transient / duty by ``float.hex``; (X) refusals naming their field; (UT) the
file's probe mesh = the Python Study's; (I-Py) θ_file = θ_Python at
``L0_BAND``; (R-box) U60's five readings at ``INTERIOR_BAND`` and to their
printed ``.4e`` digits; (A-box) the whole-box identity at ``B_BAND``; (U′) the
summary is its own result; (I-tr) the duty pass-through, an identity only.
Controls (rule (e)): (c-wall) the h-doubled tracked file, asserted, backed by
U60's 3.5712e-01; (c-ball) the ball's empty-probe-list / no-thermal refusals,
pinned elsewhere (W0 / RG1; retargeted to ``probes: []`` by §9 item 545);
(c0) by construction (RG7a / RG7b).
"""

from __future__ import annotations

import json
import os
import time

import numpy as np
import pytest
from mpi4py import MPI

from fem_em_solver import (
    DutyCycle,
    FHumanBirdcageGeometry,
    FourPortBirdcageGeometry,
    RobinCooling,
    SphereInBoxGeometry,
    Study,
    StudyConfigError,
    ThermalBallGeometry,
    ThermalBoxGeometry,
    TransientSpec,
    UniformHeatSource,
    cli,
    load_study,
    run_study_file,
)
from fem_em_solver.study_yaml import _load, _study_from_mapping

from tests.complex_mode import IS_COMPLEX
from tests.package.test_wf1_study_yaml import (
    LOGS,
    WF1,
    _cmp,
    _fresh_dir,
    _mapping,
    _rc_all,
    _spec_diff,
)
from tests.package.test_wf15_study_probe_thermal import (
    A_BOXPROBE,
    BOX_CENTRES,
    BOX_DIVISIONS,
    BOX_PROBES,
    WHOLE_BOX,
    _domain_mean,
    _series_at,
)
from tests.package.test_wf15_study_transient import _on, _same_mesh, _theta_owned_rel
from tests.validation.test_mat5_em_coupled_sphere import L0_BAND
from tests.validation.test_mat10_sar_to_heat_sphere import B_BAND
from tests.validation.test_mat12_transient_sar_heating_sphere import RHO_C
from tests.validation.test_mat13_transient_perfusion_sphere import INTERIOR_BAND
from tests.validation.test_mat15_coil_driven_heating_phantom import H_GEL, K_GEL, T_INF
from tests.validation.test_mat21_box_closed_form import HALF, Q
from tests.validation.test_mat21_box_duty_closed_form import DT_BOX, HALF_T, N_STEPS

S17_ENV = "FEM_EM_WF1_S17"
STEADY = WF1 / "thermal_box_steady_m4_probes.yaml"
WRONG_WALL = WF1 / "thermal_box_wrong_wall_m4_probes.yaml"
DUTY = WF1 / "thermal_box_duty_m4_probes.yaml"
MALFORMED = WF1 / "malformed_thermal_box_probe_radius_zero.yaml"
BALL = WF1 / "thermal_ball_duty_a.yaml"

# U60's records (real build), transcribed with their log lines.
U60_R_BOX = ("2.6992e-03", "2.6936e-03", "2.7066e-03", "2.6468e-03", "2.6443e-03")
# ^ 20261018T123416Z_WF-15-s6-W2-box.log:64-68
U60_C_WALL = 3.5712e-01  # ibid. :72 (printed beside; the assertion is > INTERIOR_BAND)
U60_A_BOX = 2.495e-13  # ibid. :63 (printed beside, never compared)


def _knob() -> bool:
    raw = os.environ.get(S17_ENV, "").strip()
    if raw == "":
        return False
    if raw == "1":
        return True
    raise ValueError(f"{S17_ENV}={raw!r}: only unset/empty or '1' are accepted")


def _p0(*args):
    if MPI.COMM_WORLD.rank == 0:
        print(*args, flush=True)


def _say(capsys, msg):
    if MPI.COMM_WORLD.rank == 0:
        with capsys.disabled():
            print(msg, flush=True)


def _all(comm, flag):
    return bool(comm.allreduce(bool(flag), op=MPI.LAND))


def _field(fn):
    try:
        fn()
    except StudyConfigError as exc:
        return exc.field
    return None


def _py_study():
    return Study(ThermalBoxGeometry(HALF, BOX_DIVISIONS), None, UniformHeatSource(Q))


def _probe_hex(probes):
    return [([float(x).hex() for x in pb.centre], float(pb.radius).hex()) for pb in probes]


# ------------------------------------------------------------------ no mesh


def test_wf1_s17_parse_and_refusals(monkeypatch, capsys):
    """(K), (X): no mesh, no solve, either build."""
    _knob()
    comm = MPI.COMM_WORLD
    calls = {"n": 0}

    def stub(self, comm=MPI.COMM_WORLD):
        calls["n"] += 1
        raise AssertionError("build must not be called by a parse")

    for cls in (SphereInBoxGeometry, FHumanBirdcageGeometry, FourPortBirdcageGeometry,
                ThermalBallGeometry, ThermalBoxGeometry):
        monkeypatch.setattr(cls, "build", stub)

    good = (STEADY, WRONG_WALL, DUTY)
    diffs = {f.name: _spec_diff(load_study(f), _py_study()) for f in good}
    _say(capsys, f"\n[WF-1 s17] (K) spec diffs {diffs}")
    assert all(d == [] for d in diffs.values()), diffs

    want = {
        STEADY: (RobinCooling(K_GEL, H_GEL, T_INF), None, None, [WHOLE_BOX, *BOX_PROBES]),
        WRONG_WALL: (RobinCooling(K_GEL, 2.0 * H_GEL, T_INF), None, None, BOX_PROBES[:1]),
        DUTY: (RobinCooling(K_GEL, H_GEL, T_INF), TransientSpec(DT_BOX, N_STEPS, RHO_C),
               DutyCycle(2.0 * HALF_T, 0.5), list(BOX_PROBES)),
    }
    for f, (cool, tr, du, probes) in want.items():
        _, kind, outputs, th = _load(f)
        out = []
        _cmp(th.cooling, cool, "cooling", out)
        _cmp(th.transient, tr, "transient", out)
        _cmp(th.duty, du, "duty", out)
        got_p, want_p = _probe_hex(outputs["probes"]), _probe_hex(probes)
        _say(capsys, f"  (K) {f.name}: kind {kind!r}, {len(outputs['probes'])} probes, hex equal "
             f"{got_p == want_p}; thermal diffs {out}")
        assert kind == "thermal_box"
        assert set(outputs) == {"probes"}
        assert got_p == want_p, (f.name, got_p, want_p)
        assert out == [], (f.name, out)

    for f in good:
        capsys.readouterr()
        rc = cli.main(["check", str(f)])
        lines = capsys.readouterr().out.splitlines()
        rc_all = _rc_all(rc, comm)
        _say(capsys, f"  (K) femem check {f.name}: rc {rc_all}")
        assert rc_all == (0, 0)
        if comm.rank == 0:
            assert [ln for ln in lines if ln.startswith("family: ")] == ["family: thermal_box"]

    # (X) the tracked malformed file through load, run and check
    assert _field(lambda: load_study(MALFORMED)) == "outputs.probes[1].radius"
    d = LOGS / "thermal_box_probes_malformed"
    _fresh_dir(d, comm)
    capsys.readouterr()
    rc_run = cli.main(["run", str(MALFORMED), "--output-dir", str(d)])
    rc_check = cli.main(["check", str(MALFORMED)])
    capsys.readouterr()
    run_all, check_all = _rc_all(rc_run, comm), _rc_all(rc_check, comm)
    exists = comm.allreduce(bool(d.exists()), op=MPI.LOR)
    _say(capsys, f"  (X) malformed radius 0: run rc {run_all}, check rc {check_all}, out exists "
         f"{exists}")
    assert run_all == (2, 2) and check_all == (2, 2)
    assert not exists

    s = _mapping(STEADY)
    base = {k: v for k, v in s.items() if k != "outputs"}
    ok = {"centre": [0.0, 0.0, 0.0], "radius": 0.01}
    rows = [
        ({**base, "outputs": {"probes": []}}, "outputs.probes"),
        ({**base, "outputs": {"probes": {"centre": [0.0, 0.0, 0.0]}}}, "outputs.probes"),
        ({**base, "outputs": {"probes": [1]}}, "outputs.probes[0]"),
        ({**base, "outputs": {"probes": [{**ok, "mass": 1.0}]}}, "outputs.probes[0].mass"),
        ({**base, "outputs": {"probes": [{"centre": [0.0, 0.0, 0.0]}]}},
         "outputs.probes[0].radius"),
        ({**base, "outputs": {"probes": [{**ok, "centre": [0.0, 0.0]}]}},
         "outputs.probes[0].centre"),
        ({**base, "outputs": {"probes": [ok], "sar": {}}}, "outputs.sar"),
        ({**base, "outputs": {"probes": [ok]}}, None),
    ]
    got = [_field(lambda m=m: _study_from_mapping(m)) for m, _ in rows]
    for (_, want_f), g in zip(rows, got):
        _say(capsys, f"  (X) want {want_f!r} got {g!r}")
    assert got == [f for _, f in rows], [(f, g) for (_, f), g in zip(rows, got) if f != g]
    # (c-ball) the ball refuses an empty probe list (retargeted by §9 item 545,
    # review ruling 2026-10-19 03:00 Q1: a valid ball probe is admitted since then).
    b = _mapping(BALL)
    ball_f = _field(lambda: _study_from_mapping({**b, "outputs": {"probes": []}}))
    _say(capsys, f"  (c-ball) thermal_ball with an empty probe list: field {ball_f!r}; build "
         f"calls {calls['n']}")
    assert ball_f == "outputs.probes"
    assert calls["n"] == 0


# ------------------------------------------------------------------ solves


def _require():
    if not _knob():
        pytest.skip(f"{S17_ENV} unset: the thermal-box probe solve legs are environment-"
                    "selected (real build, -n 2)")
    comm = MPI.COMM_WORLD
    assert not IS_COMPLEX, "U60's box probe records are real-build records"
    assert comm.size == 2, comm.size
    return comm


def _check_summary(run, comm, label):
    """(U′) the summary probes are the run's own values; JSON on disk = memory."""
    s = run.summary
    js = json.loads(run.paths["summary"].read_text())
    vals_ok = s["probes"]["values_k"] == [float(v) for v in run.probes.values]
    rest = {k: v for k, v in s.items() if k != "elapsed_s"}
    js_rest = {k: v for k, v in js.items() if k != "elapsed_s"}
    same_json = _all(comm, js_rest == rest and (comm.rank != 0 or js == s))
    _p0(f"  (U′) {label}: values == run.probes.values {vals_ok}; JSON == summary {same_json}; "
        f"keys {sorted(s)}")
    assert vals_ok
    assert same_json
    assert s["family"] == "thermal_box"
    assert run.result is None
    assert s["probes"]["radii"] == [float(pb.radius) for pb in run.probes.probes]


def test_wf1_s17_box_probes_reproduce_u60():
    comm = _require()
    t00 = time.perf_counter()
    out = LOGS / "thermal_box_probes_steady"
    _fresh_dir(out, comm)
    run = run_study_file(STEADY, output_dir=out, comm=comm)
    t_file = time.perf_counter() - t00
    pr = run.probes
    py = _py_study().probe([WHOLE_BOX, *BOX_PROBES], RobinCooling(K_GEL, H_GEL, T_INF))
    m_f = pr.theta.function_space.mesh
    m_p = py.theta.function_space.mesh
    ca, cb, geo = _same_mesh(m_f, m_p, comm)
    # (I-Py) on the temperatures (T - T_INF normalises; θ itself is the rise).
    ip = _theta_owned_rel(_on(m_p, pr.thermal.temperature), py.thermal.temperature, comm)
    # run.thermal (the temperature block) vs the probe's re-solve
    two = _theta_owned_rel(_on(run.thermal.temperature.function_space.mesh,
                               pr.thermal.temperature), run.thermal.temperature, comm)
    vals_f = [float(v) for v in pr.values]
    vals_p = [float(v) for v in py.values]
    _p0(f"WF1S17 (UT) steady: cells {ca} / {cb}; geometry equal {geo}; file run {t_file:.1f}s")
    _p0(f"WF1S17 (I-Py) steady: owned max|dT|/max|T-Tinf| = {ip:.3e} (band {L0_BAND:.0e}); "
        f"temperature block vs probe re-solve {two:.3e}")
    _p0(f"WF1S17 probe values file {vals_f!r}")
    _p0(f"WF1S17 probe values py   {vals_p!r} (printed only)")
    mean_h, v_h = _domain_mean(pr.theta, comm)
    a_val = abs(vals_f[0] - mean_h) / abs(mean_h)
    a_vol = abs(float(pr.volumes[0]) / v_h - 1.0)
    _p0(f"WF1S17 (A-box) |P0 - mean| / |mean| = {a_val:.3e}, |V0/V_h - 1| = {a_vol:.3e} (band "
        f"{B_BAND:.0e}; U60 {U60_A_BOX:.3e})")
    mv = Q * A_BOXPROBE**2 / (10.0 * K_GEL)
    ratios, digits = [], []
    for j, c in enumerate(BOX_CENTRES):
        cf = _series_at(c) - mv
        r = abs(vals_f[1 + j] / cf - 1.0)
        ratios.append(r)
        digits.append(f"{r:.4e}")
        _p0(f"WF1S17 (R-box) c{j} = {c}: P = {vals_f[1 + j]:.6e}, ball mean = {cf:.6e}, rel = "
            f"{r:.6e} -> {r:.4e} (U60 {U60_R_BOX[j]}; band {INTERIOR_BAND:.0e})")
    # (c-wall) through the tracked file
    out_w = LOGS / "thermal_box_probes_wrong_wall"
    _fresh_dir(out_w, comm)
    run_w = run_study_file(WRONG_WALL, output_dir=out_w, comm=comm)
    c_wall = abs(float(run_w.probes.values[0]) / (_series_at(BOX_CENTRES[0]) - mv) - 1.0)
    _p0(f"WF1S17 control (c-wall) wrong-wall file at c0: miss {c_wall:.4e} (U60 "
        f"{U60_C_WALL:.4e}; equal digits {f'{c_wall:.4e}' == f'{U60_C_WALL:.4e}'}; asserted > "
        f"{INTERIOR_BAND:.0e})")
    _check_summary(run, comm, "steady")
    _check_summary(run_w, comm, "wrong wall")
    assert ca == cb and geo
    assert ip <= L0_BAND, f"(I-Py) {ip:.3e}"
    assert two <= L0_BAND, f"temperature block vs probe {two:.3e}"
    assert a_val <= B_BAND and a_vol <= B_BAND
    for r in ratios:
        assert r <= INTERIOR_BAND
    assert tuple(digits) == U60_R_BOX, digits
    assert c_wall > INTERIOR_BAND
    _p0(f"WF1S17 W1 steady leg {time.perf_counter() - t00:.1f}s")


def test_wf1_s17_transient_probes_equal_python():
    comm = _require()
    t00 = time.perf_counter()
    out = LOGS / "thermal_box_probes_duty"
    _fresh_dir(out, comm)
    run = run_study_file(DUTY, output_dir=out, comm=comm)
    pr = run.probes
    py = _py_study().probe(
        BOX_PROBES, RobinCooling(K_GEL, H_GEL, T_INF),
        transient=TransientSpec(DT_BOX, N_STEPS, RHO_C), duty=DutyCycle(2.0 * HALF_T, 0.5),
    )
    m_f = pr.theta.function_space.mesh
    m_p = py.theta.function_space.mesh
    ca, cb, geo = _same_mesh(m_f, m_p, comm)
    ip = _theta_owned_rel(_on(m_p, pr.thermal.temperature), py.thermal.temperature, comm)
    mult = np.array_equal(np.asarray(pr.thermal.source_multipliers),
                          np.asarray(py.thermal.source_multipliers))
    mult_s = np.array_equal(np.asarray(run.summary["temperature"]["source_multipliers"]),
                            np.asarray(py.thermal.source_multipliers))
    _p0(f"WF1S17 (I-tr) duty: cells {ca} / {cb}; geometry equal {geo}; owned max|dT|/max|T-Tinf| "
        f"= {ip:.3e} (band {L0_BAND:.0e}); multipliers equal {mult} / summary {mult_s} "
        f"({len(py.thermal.source_multipliers)} entries)")
    _p0(f"WF1S17 (I-tr) probe values file {[float(v) for v in pr.values]!r}")
    _p0(f"WF1S17 (I-tr) probe values py   {[float(v) for v in py.values]!r} (printed; no "
        "comparand is posed)")
    _check_summary(run, comm, "duty")
    assert ca == cb and geo
    assert ip <= L0_BAND, f"(I-tr) {ip:.3e}"
    assert mult and mult_s
    assert run.summary["temperature"]["duty"] == {"period_s": 2.0 * HALF_T, "duty": 0.5}
    _p0(f"WF1S17 W1 duty leg {time.perf_counter() - t00:.1f}s")
