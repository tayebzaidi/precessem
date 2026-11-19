"""`WF-1` step 19: ``outputs.probes`` on a ``thermal_ball`` file drives ``Study.probe``
through ``femem run`` and reproduces `WF-15` step 6's (§9 item 523, U60) ball probes,
`MAT-19`'s virtual-probe gate.

Brief: docs/planning/chunks/WF-1.md, "§9 item 545" and the review ruling
2026-10-19 03:00.  The solve legs run in the real build (the build of U60's
records, ``20261018T123401Z_WF-15-s6-W1-ball.log:50``) and only with
``FEM_EM_WF1_S19=1`` (then they assert the real build and ``-n 2``, never
skip).  The parse node builds no mesh and runs in either build.

Anchors: (K) the YAML Study = the Python Study and the parsed probes / cooling
/ transient / duty by ``float.hex``; (X) refusals naming their field, including
a walled ball with probes (ruling Q3); (UT) the file's probe mesh = the Python
Study's; (I-Py) T_file = T_Python at ``L0_BAND``; (T*-ball) U60's three readings
at ``P_BAND`` and to their printed ``.4e`` digits; (A-ball) the whole-ball
identity at ``B_BAND``; (U′) the summary is its own result; (I-du) the duty
pass-through, an identity only.  Controls (rule (e)): (c-disp) the displaced
probe through the tracked file, asserted, backed by U60's 5.4208e-02; (c-sar)
the ball ``outputs.sar`` / no-thermal sphere refusals, pinned elsewhere (W0 /
RG1); (c0) by construction (RG6a / RG6b).
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
from fem_em_solver.thermal.probe import ProbeBall
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
    BALL_PROBES,
    DISP_PROBE,
    WHOLE_BALL,
    _domain_mean,
)
from tests.package.test_wf15_study_transient import RES, _on, _same_mesh, _theta_owned_rel
from tests.validation.test_mat5_em_coupled_sphere import L0_BAND
from tests.validation.test_mat9_transient_conduction_sphere import RHO_C as RHO_C9
from tests.validation.test_mat9_transient_conduction_sphere import T_INF as T_INF9
from tests.validation.test_mat9_transient_conduction_sphere import H as H9
from tests.validation.test_mat9_transient_conduction_sphere import K as K9
from tests.validation.test_mat9_transient_conduction_sphere import Q as Q9
from tests.validation.test_mat9_transient_conduction_sphere import R as R9
from tests.validation.test_mat10_sar_to_heat_sphere import B_BAND, P_BAND, T_INF, H, K, R
from tests.validation.test_mat12_transient_sar_heating_sphere import N_MAIN, RHO_C, T_STAR
from tests.validation.test_mat17_duty_cycle_closed_form import CASES
from tests.validation.test_mat17_duty_cycle_fem import T_END_A, _n
from tests.validation.test_mat19_virtual_temperature_probe import (
    A_PROBE,
    CENTRES,
    DISPLACED,
    _ball_mean,
    _q_cf,
)

S19_ENV = "FEM_EM_WF1_S19"
M19 = WF1 / "thermal_ball_m19_probes.yaml"
DUTY = WF1 / "thermal_ball_duty_a_probes.yaml"
DUTY_PLAIN = WF1 / "thermal_ball_duty_a.yaml"
WALL = WF1 / "thermal_ball_wall.yaml"
MALFORMED = WF1 / "malformed_thermal_ball_probe_radius_zero.yaml"

DUTY_PROBES = [ProbeBall((0.0, 0.0, 0.0), 2.0 * R9), ProbeBall((0.0, 0.0, 0.0), 0.02)]

# U60's records (real build), transcribed with their log lines.
U60_T_BALL = ("9.7501e-04", "5.0491e-04", "9.4257e-04")
# ^ 20261018T123401Z_WF-15-s6-W1-ball.log:79-81
U60_C_DISP = 5.4208e-02  # ibid. :83 (printed beside; the assertion is > P_BAND)
U60_A_BALL = 3.007e-16  # ibid. :78 (printed beside, never compared)
U60_BALL_CELLS = 19489  # ibid. :77 (printed beside)


def _knob() -> bool:
    raw = os.environ.get(S19_ENV, "").strip()
    if raw == "":
        return False
    if raw == "1":
        return True
    raise ValueError(f"{S19_ENV}={raw!r}: only unset/empty or '1' are accepted")


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


def _py_m19():
    return Study(ThermalBallGeometry(R, 1.0e-3), None, UniformHeatSource(_q_cf()))


def _py_duty():
    return Study(ThermalBallGeometry(R9, RES), None, UniformHeatSource(Q9))


def _probe_hex(probes):
    return [([float(x).hex() for x in pb.centre], float(pb.radius).hex()) for pb in probes]


# ------------------------------------------------------------------ no mesh


def test_wf1_s19_parse_and_refusals(monkeypatch, capsys):
    """(K), (X): no mesh, no solve, either build."""
    _knob()
    comm = MPI.COMM_WORLD
    calls = {"n": 0}

    def stub(self, comm=MPI.COMM_WORLD):
        calls["n"] += 1
        raise AssertionError("build must not be called by a parse")

    for cls in (
        SphereInBoxGeometry,
        FHumanBirdcageGeometry,
        FourPortBirdcageGeometry,
        ThermalBallGeometry,
        ThermalBoxGeometry,
    ):
        monkeypatch.setattr(cls, "build", stub)

    diffs = {
        M19.name: _spec_diff(load_study(M19), _py_m19()),
        DUTY.name: _spec_diff(load_study(DUTY), load_study(DUTY_PLAIN)),
    }
    _say(capsys, f"\n[WF-1 s19] (K) spec diffs {diffs}; q repr {_q_cf()!r}")
    assert all(d == [] for d in diffs.values()), diffs

    want = {
        M19: (
            RobinCooling(K, H, T_INF),
            TransientSpec(T_STAR / N_MAIN, N_MAIN, RHO_C),
            None,
            [WHOLE_BALL, *BALL_PROBES, DISP_PROBE],
        ),
        DUTY: (
            RobinCooling(K9, H9, T_INF9),
            TransientSpec(150.0, _n(T_END_A, 150.0), RHO_C9),
            DutyCycle(*CASES["A"]),
            DUTY_PROBES,
        ),
    }
    for f, (cool, tr, du, probes) in want.items():
        _, kind, outputs, th = _load(f)
        out = []
        _cmp(th.cooling, cool, "cooling", out)
        _cmp(th.transient, tr, "transient", out)
        _cmp(th.duty, du, "duty", out)
        got_p, want_p = _probe_hex(outputs["probes"]), _probe_hex(probes)
        _say(
            capsys,
            f"  (K) {f.name}: kind {kind!r}, {len(outputs['probes'])} probes, hex equal "
            f"{got_p == want_p}; thermal diffs {out}",
        )
        if got_p != want_p:
            _say(
                capsys,
                f"    want reprs {[([repr(float(x)) for x in pb.centre], repr(float(pb.radius))) for pb in probes]}",
            )
        assert kind == "thermal_ball"
        assert set(outputs) == {"probes"}
        assert got_p == want_p, (f.name, got_p, want_p)
        assert out == [], (f.name, out)

    for f in (M19, DUTY):
        capsys.readouterr()
        rc = cli.main(["check", str(f)])
        lines = capsys.readouterr().out.splitlines()
        rc_all = _rc_all(rc, comm)
        _say(capsys, f"  (K) femem check {f.name}: rc {rc_all}")
        assert rc_all == (0, 0)
        if comm.rank == 0:
            assert [ln for ln in lines if ln.startswith("family: ")] == ["family: thermal_ball"]

    # (X) the tracked malformed file through load, run and check
    assert _field(lambda: load_study(MALFORMED)) == "outputs.probes[1].radius"
    d = LOGS / "thermal_ball_probes_malformed"
    _fresh_dir(d, comm)
    capsys.readouterr()
    rc_run = cli.main(["run", str(MALFORMED), "--output-dir", str(d)])
    rc_check = cli.main(["check", str(MALFORMED)])
    capsys.readouterr()
    run_all, check_all = _rc_all(rc_run, comm), _rc_all(rc_check, comm)
    exists = comm.allreduce(bool(d.exists()), op=MPI.LOR)
    _say(
        capsys,
        f"  (X) malformed radius 0: run rc {run_all}, check rc {check_all}, out exists "
        f"{exists}",
    )
    assert run_all == (2, 2) and check_all == (2, 2)
    assert not exists

    s = _mapping(M19)
    base = {k: v for k, v in s.items() if k != "outputs"}
    ok = {"centre": [0.0, 0.0, 0.0], "radius": 0.002}
    w = _mapping(WALL)
    rows = [
        ({**base, "outputs": {"probes": []}}, "outputs.probes"),
        ({**base, "outputs": {"probes": {"centre": [0.0, 0.0, 0.0]}}}, "outputs.probes"),
        ({**base, "outputs": {"probes": [1]}}, "outputs.probes[0]"),
        ({**base, "outputs": {"probes": [{**ok, "mass": 1.0}]}}, "outputs.probes[0].mass"),
        (
            {**base, "outputs": {"probes": [{"centre": [0.0, 0.0, 0.0]}]}},
            "outputs.probes[0].radius",
        ),
        (
            {**base, "outputs": {"probes": [{**ok, "centre": [0.0, 0.0]}]}},
            "outputs.probes[0].centre",
        ),
        ({**base, "outputs": {"probes": [ok], "sar": {}}}, "outputs.sar"),
        # ruling 2026-10-19 03:00 Q3: a walled ball with probes is refused
        ({**w, "outputs": {"probes": [ok]}}, "outputs.probes"),
        ({**base, "outputs": {"probes": [ok]}}, None),
    ]
    got = [_field(lambda m=m: _study_from_mapping(m)) for m, _ in rows]
    for (_, want_f), g in zip(rows, got):
        _say(capsys, f"  (X) want {want_f!r} got {g!r}")
    assert got == [f for _, f in rows], [(f, g) for (_, f), g in zip(rows, got) if f != g]
    _say(capsys, f"  build calls {calls['n']}")
    assert calls["n"] == 0


# ------------------------------------------------------------------ solves


def _require():
    if not _knob():
        pytest.skip(
            f"{S19_ENV} unset: the thermal-ball probe solve legs are environment-"
            "selected (real build, -n 2)"
        )
    comm = MPI.COMM_WORLD
    assert not IS_COMPLEX, "U60's ball probe records are real-build records"
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
    _p0(
        f"  (U′) {label}: values == run.probes.values {vals_ok}; JSON == summary {same_json}; "
        f"keys {sorted(s)}"
    )
    assert vals_ok
    assert same_json
    assert s["family"] == "thermal_ball"
    assert run.result is None
    assert s["probes"]["radii"] == [float(pb.radius) for pb in run.probes.probes]


def test_wf1_s19_ball_probes_reproduce_u60():
    comm = _require()
    t00 = time.perf_counter()
    out = LOGS / "thermal_ball_probes_m19"
    _fresh_dir(out, comm)
    run = run_study_file(M19, output_dir=out, comm=comm)
    t_file = time.perf_counter() - t00
    pr = run.probes
    q = _q_cf()
    py = _py_m19().probe(
        [WHOLE_BALL, *BALL_PROBES, DISP_PROBE],
        RobinCooling(K, H, T_INF),
        transient=TransientSpec(T_STAR / N_MAIN, N_MAIN, RHO_C),
    )
    m_f = pr.theta.function_space.mesh
    m_p = py.theta.function_space.mesh
    ca, cb, geo = _same_mesh(m_f, m_p, comm)
    ip = _theta_owned_rel(_on(m_p, pr.thermal.temperature), py.thermal.temperature, comm)
    two = _theta_owned_rel(
        _on(run.thermal.temperature.function_space.mesh, pr.thermal.temperature),
        run.thermal.temperature,
        comm,
    )
    vals_f = [float(v) for v in pr.values]
    vals_p = [float(v) for v in py.values]
    _p0(
        f"WF1S19 (UT) m19: cells {ca} / {cb} (U60 {U60_BALL_CELLS}); geometry equal {geo}; "
        f"file run {t_file:.1f}s"
    )
    _p0(
        f"WF1S19 (I-Py) m19: owned max|dT|/max|T-Tinf| = {ip:.3e} (band {L0_BAND:.0e}); "
        f"temperature block vs probe re-solve {two:.3e}"
    )
    _p0(f"WF1S19 probe values file {vals_f!r}")
    _p0(f"WF1S19 probe values py   {vals_p!r} (printed only)")
    mean_h, v_h = _domain_mean(pr.theta, comm)
    a_val = abs(vals_f[0] - mean_h) / abs(mean_h)
    a_vol = abs(float(pr.volumes[0]) / v_h - 1.0)
    _p0(
        f"WF1S19 (A-ball) |P0 - mean| / |mean| = {a_val:.3e}, |V0/V_h - 1| = {a_vol:.3e} (band "
        f"{B_BAND:.0e}; U60 {U60_A_BALL:.3e})"
    )
    ratios, digits = [], []
    for j, c in enumerate(CENTRES):
        cf = _ball_mean(float(np.linalg.norm(c)), A_PROBE, T_STAR, q)
        r = abs(vals_f[1 + j] / cf - 1.0)
        ratios.append(r)
        digits.append(f"{r:.4e}")
        _p0(
            f"WF1S19 (T*-ball) c{j}: P = {vals_f[1 + j]:.6e}, ball mean = {cf:.6e}, rel = "
            f"{r:.6e} -> {r:.4e} (U60 {U60_T_BALL[j]}; band {P_BAND:.0e})"
        )
    cf2 = _ball_mean(float(np.linalg.norm(CENTRES[2])), A_PROBE, T_STAR, q)
    c_disp = abs(vals_f[-1] / cf2 - 1.0)
    sep = abs(_ball_mean(float(np.linalg.norm(DISPLACED)), A_PROBE, T_STAR, q) / cf2 - 1.0)
    _p0(
        f"WF1S19 control (c-disp) displaced vs c2's comparand: miss {c_disp:.4e} (U60 "
        f"{U60_C_DISP:.4e}; equal digits {f'{c_disp:.4e}' == f'{U60_C_DISP:.4e}'}; comparand "
        f"separation {sep:.4e}; asserted > {P_BAND:.0e})"
    )
    _check_summary(run, comm, "m19")
    assert ca == cb and geo
    assert ip <= L0_BAND, f"(I-Py) {ip:.3e}"
    assert two <= L0_BAND, f"temperature block vs probe {two:.3e}"
    assert a_val <= B_BAND and a_vol <= B_BAND
    for r in ratios:
        assert r <= P_BAND
    assert tuple(digits) == U60_T_BALL, digits
    assert c_disp > P_BAND
    _p0(f"WF1S19 W1 m19 leg {time.perf_counter() - t00:.1f}s")


def test_wf1_s19_duty_probes_equal_python():
    comm = _require()
    t00 = time.perf_counter()
    out = LOGS / "thermal_ball_probes_duty"
    _fresh_dir(out, comm)
    run = run_study_file(DUTY, output_dir=out, comm=comm)
    pr = run.probes
    py = _py_duty().probe(
        DUTY_PROBES,
        RobinCooling(K9, H9, T_INF9),
        transient=TransientSpec(150.0, _n(T_END_A, 150.0), RHO_C9),
        duty=DutyCycle(*CASES["A"]),
    )
    m_f = pr.theta.function_space.mesh
    m_p = py.theta.function_space.mesh
    ca, cb, geo = _same_mesh(m_f, m_p, comm)
    ip = _theta_owned_rel(_on(m_p, pr.thermal.temperature), py.thermal.temperature, comm)
    mult = np.array_equal(
        np.asarray(pr.thermal.source_multipliers), np.asarray(py.thermal.source_multipliers)
    )
    mult_s = np.array_equal(
        np.asarray(run.summary["temperature"]["source_multipliers"]),
        np.asarray(py.thermal.source_multipliers),
    )
    mean_h, v_h = _domain_mean(pr.theta, comm)
    vals_f = [float(v) for v in pr.values]
    a_val = abs(vals_f[0] - mean_h) / abs(mean_h)
    a_vol = abs(float(pr.volumes[0]) / v_h - 1.0)
    _p0(
        f"WF1S19 (I-du) duty: cells {ca} / {cb}; geometry equal {geo}; owned max|dT|/max|T-Tinf| "
        f"= {ip:.3e} (band {L0_BAND:.0e}); multipliers equal {mult} / summary {mult_s} "
        f"({len(py.thermal.source_multipliers)} entries)"
    )
    _p0(
        f"WF1S19 (I-du) whole ball |P0 - mean| / |mean| = {a_val:.3e}, |V0/V_h - 1| = "
        f"{a_vol:.3e} (band {B_BAND:.0e})"
    )
    _p0(f"WF1S19 (I-du) probe values file {vals_f!r}")
    _p0(
        f"WF1S19 (I-du) probe values py   {[float(v) for v in py.values]!r} (printed; no "
        "comparand is posed)"
    )
    _check_summary(run, comm, "duty")
    assert ca == cb and geo
    assert ip <= L0_BAND, f"(I-du) {ip:.3e}"
    assert mult and mult_s
    assert a_val <= B_BAND and a_vol <= B_BAND
    assert run.summary["temperature"]["duty"] == {"period_s": CASES["A"][0], "duty": CASES["A"][1]}
    _p0(f"WF1S19 W1 duty leg {time.perf_counter() - t00:.1f}s")
