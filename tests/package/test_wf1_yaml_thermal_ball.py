"""`WF-1` step 14: ``femem run`` on thermal-only ball files reproduces `WF-15`
step 3's `MAT-17` (duty cycle) and `MAT-18` (container wall) records.

Brief: docs/planning/chunks/WF-1.md, "§9 item 500" and its review ruling
(2026-10-16 18:00).  The solve legs are gated in the real build, the build of
the records; they run only with ``FEM_EM_WF1_S14=1`` (and then assert the real
build and ``-n 2``, never skip).  The parse node builds no mesh and runs in
either build.

Anchors: (K) the YAML Study = the Python Study (``_spec_diff``) and the parsed
``RobinCooling`` / ``TransientSpec`` / ``DutyCycle`` by ``float.hex``; (X)
refusals naming their field; (UT) bitwise upstream of the solve (cells, mesh
geometry, multipliers, DG0 k / rho_c / q) with (S0); (I-Py) the file's theta =
the Python Study's at ``L0_BAND``; (G) ``INTERIOR_BAND`` / ``BALANCE_BAND``;
(R) the records at ``.4e``; (I=) the summary = its own thermal result.
Controls (rule (e), asserted): (c1) the schedule dropped through a file,
(c2) the wall's rho_c dropped through a file, (c0) dropped keys at parse.
"""

from __future__ import annotations

import copy
import json
import os
import time

import numpy as np
import pytest
from mpi4py import MPI

from fem_em_solver import (
    ContainerWall,
    DutyCycle,
    FHumanBirdcageGeometry,
    FourPortBirdcageGeometry,
    RobinCooling,
    SphereInBoxGeometry,
    Study,
    StudyConfigError,
    ThermalBallGeometry,
    TransientSpec,
    UniformHeatSource,
    cli,
    load_study,
    run_study_file,
)
from fem_em_solver.study_yaml import _load, _study_from_mapping
from tests.complex_mode import IS_COMPLEX
from tests.materials.test_phantom_material_model import _dg0_owned_diff
from tests.package.test_wf1_study_yaml import (
    LOGS,
    WF1,
    _cmp,
    _fresh_dir,
    _mapping,
    _rc_all,
    _spec_diff,
)
from tests.package.test_wf15_study_transient import (
    MAT17_A_PRINT,
    MAT17_B_PRINT,
    MAT18_C1_UNIFORM_RHOC,
    MAT18_T1_CELLS,
    MAT18_T1_PRINT,
    RES,
    _duty_study,
    _exact_duty,
    _on,
    _same_mesh,
    _theta_owned_rel,
    _wall_study,
)
from tests.validation.test_mat5_em_coupled_sphere import L0_BAND
from tests.validation.test_mat9_transient_conduction_sphere import (
    RHO_C,
    T_INF,
    H,
    K,
    Q,
    R,
    _interior_error,
)
from tests.validation.test_mat17_duty_cycle_closed_form import CASES
from tests.validation.test_mat17_duty_cycle_fem import (
    BALANCE_BAND,
    DT_B,
    INTERIOR_BAND,
    ON_TIME_A,
    ON_TIME_B,
    T_END_A,
    T_END_B,
    _census,
    _n,
)
from tests.validation.test_mat18_two_layer_sphere_closed_form import (
    A_,
    B_,
    K1,
    K2,
    Q1,
    RC1,
    RC2,
)
from tests.validation.test_mat18_two_layer_sphere_transient_fem import (
    CELL_GUARD,
    N_CN,
    T_STAR,
    _errors,
    _fine,
)

S14_ENV = "FEM_EM_WF1_S14"
DUTY_A = WF1 / "thermal_ball_duty_a.yaml"
DUTY_B = WF1 / "thermal_ball_duty_b.yaml"
WALL = WF1 / "thermal_ball_wall.yaml"
MALFORMED_RHO_C = WF1 / "malformed_thermal_ball_rho_c_string.yaml"
SPHERE_THERMAL = WF1 / "sphere_mat10_thermal.yaml"
P_A, D_A = CASES["A"]
P_B, D_B = CASES["B"]

# `WF-15` step 3's readings on the same fixtures (printed beside, never asserted).
WF15_A_READ = 4.3869840018e-04  # 20261008T213834Z_WF-15-s3-W1a-duty.log:91
WF15_B_READ = 1.0682780964e-03  # 20261008T213834Z_WF-15-s3-W1a-duty.log:96
WF15_T1_READ = 8.2066325346e-04  # 20261008T214046Z_WF-15-s3-W1b-wall.log:91
WF15_C1_READ = 9.7087e-01  # 20261008T213834Z_WF-15-s3-W1a-duty.log:98
WF15_C2_READ = 5.3140e-02  # 20261008T214046Z_WF-15-s3-W1b-wall.log:93
WF15_CELLS = 151548  # 20261008T213834Z_WF-15-s3-W1a-duty.log:89


def _knob() -> bool:
    raw = os.environ.get(S14_ENV, "").strip()
    if raw == "":
        return False
    if raw == "1":
        return True
    raise ValueError(f"{S14_ENV}={raw!r}: only unset/empty or '1' are accepted")


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


def _without(m, *path):
    m = copy.deepcopy(m)
    d = m
    for p in path[:-1]:
        d = d[p]
    del d[path[-1]]
    return m


def _with(m, value, *path):
    m = copy.deepcopy(m)
    d = m
    for p in path[:-1]:
        d = d[p]
    d[path[-1]] = value
    return m


def _write_yaml(mapping, directory, name, comm):
    """Rank 0 writes ``mapping`` (repr floats via safe_dump), then ``Barrier``."""
    import yaml

    _fresh_dir(directory, comm)
    f = directory / name
    if comm.rank == 0:
        directory.mkdir(parents=True, exist_ok=True)
        f.write_text(yaml.safe_dump(mapping, sort_keys=False))
    comm.Barrier()
    return f


def _hex_tree_equal(a, b):
    if isinstance(a, float) and isinstance(b, float):
        return a.hex() == b.hex()
    if isinstance(a, dict) and isinstance(b, dict):
        return a.keys() == b.keys() and all(_hex_tree_equal(a[k], b[k]) for k in a)
    if isinstance(a, list) and isinstance(b, list):
        return len(a) == len(b) and all(_hex_tree_equal(x, y) for x, y in zip(a, b))
    return type(a) is type(b) and a == b


# ------------------------------------------------------------------ no mesh


def test_wf1_s14_parse_and_refusals(monkeypatch, capsys):
    """(K), (X), (c0): no mesh, no solve, either build."""
    _knob()  # a bad knob value raises here too
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
    ):
        monkeypatch.setattr(cls, "build", stub)

    # (K) the Study spec
    py_ball = Study(ThermalBallGeometry(R, RES), None, UniformHeatSource(Q))
    py_wall = Study(
        ThermalBallGeometry(A_, RES, wall=ContainerWall(B_ - A_, K2, RC2)),
        None,
        UniformHeatSource(Q1),
    )
    studies = {f: load_study(f) for f in (DUTY_A, DUTY_B, WALL)}
    diffs = {
        DUTY_A.name: _spec_diff(studies[DUTY_A], py_ball),
        DUTY_B.name: _spec_diff(studies[DUTY_B], py_ball),
        WALL.name: _spec_diff(studies[WALL], py_wall),
    }
    _say(capsys, f"\n[WF-1 s14] (K) spec diffs {diffs}")
    assert all(d == [] for d in diffs.values()), diffs
    for st in studies.values():
        assert st.frequency_hz is None and st.material is None

    # (K) the fourth element
    want = {
        DUTY_A: (
            RobinCooling(K, H, T_INF),
            TransientSpec(150.0, _n(T_END_A, 150.0), RHO_C),
            DutyCycle(*CASES["A"]),
        ),
        DUTY_B: (
            RobinCooling(K, H, T_INF),
            TransientSpec(DT_B, _n(T_END_B, DT_B), RHO_C),
            DutyCycle(*CASES["B"]),
        ),
        WALL: (RobinCooling(K1, H, T_INF), TransientSpec(T_STAR / N_CN, N_CN, RC1), None),
    }
    for f, (cool, tr, du) in want.items():
        _, kind, outputs, th = _load(f)
        out = []
        _cmp(th.cooling, cool, "cooling", out)
        _cmp(th.transient, tr, "transient", out)
        _cmp(th.duty, du, "duty", out)
        _say(
            capsys,
            f"  (K) {f.name}: kind {kind!r}, outputs {outputs!r}, thermal diffs {out}; "
            f"dt_s {th.transient.dt_s.hex()} n_steps {th.transient.n_steps}",
        )
        assert kind == "thermal_ball"
        assert outputs == {}
        assert out == [], (f.name, out)
    assert _load(WALL)[3].duty is None

    # (K) femem check
    for f in (DUTY_A, DUTY_B, WALL):
        capsys.readouterr()
        rc = cli.main(["check", str(f)])
        lines = capsys.readouterr().out.splitlines()
        rc_all = _rc_all(rc, comm)
        _say(capsys, f"  (K) femem check {f.name}: rc {rc_all}; rank-0 lines {lines}")
        assert rc_all == (0, 0)
        if comm.rank == 0:
            assert [ln for ln in lines if ln.startswith("family: ")] == ["family: thermal_ball"]
            assert "excitation: UniformHeatSource" in lines

    # (X) the tracked string-float file, load + cli run / check
    assert _field(lambda: load_study(MALFORMED_RHO_C)) == "thermal.transient.rho_c_j_per_m3_k"
    d = LOGS / "thermal_ball_malformed"
    _fresh_dir(d, comm)
    capsys.readouterr()
    rc_run = cli.main(["run", str(MALFORMED_RHO_C), "--output-dir", str(d)])
    rc_check = cli.main(["check", str(MALFORMED_RHO_C)])
    capsys.readouterr()
    run_all, check_all = _rc_all(rc_run, comm), _rc_all(rc_check, comm)
    exists = comm.allreduce(bool(d.exists()), op=MPI.LOR)
    _say(
        capsys,
        f"  (X) malformed rho_c string: run rc {run_all}, check rc {check_all}, "
        f"out exists {exists}",
    )
    assert run_all == (2, 2)
    assert check_all == (2, 2)
    assert not exists

    # (X) the variant table
    a = _mapping(DUTY_A)
    w = _mapping(WALL)
    sph = _mapping(SPHERE_THERMAL)
    tr_a = a["thermal"]["transient"]
    rows = [
        ({**a, "study": {"frequency_hz": 64000000.0}}, "study"),
        ({**a, "material": {"sigma": 0.5, "epsilon_r": 78.0}}, "material"),
        ({**a, "conductor": {"kind": "leontovich"}}, "conductor"),
        ({**a, "tune": {"k": 1}}, "tune"),
        (_without(a, "thermal"), "thermal"),
        (_without(a, "excitation"), "excitation"),
        (_with(a, 0.0, "geometry", "radius"), "geometry.radius"),
        (_without(a, "geometry", "resolution"), "geometry.resolution"),
        (_with(a, 0.1, "geometry", "radius_m"), "geometry.radius_m"),
        (_with(a, [0.01], "geometry", "wall"), "geometry.wall"),
        (_with(w, -0.01, "geometry", "wall", "thickness_m"), "geometry.wall.thickness_m"),
        (_with(w, 0.19, "geometry", "wall", "k2"), "geometry.wall.k2"),
        (_with(a, "quasi_static_sphere_field", "excitation", "kind"), "excitation.kind"),
        (_with(a, 0.0, "excitation", "q_w_per_m3"), "excitation.q_w_per_m3"),
        (_with(a, 2.0, "thermal", "source_scale"), "thermal.source_scale"),
        (_with(a, {"mass_g": 10.0}, "thermal", "balls"), "thermal.balls"),
        (_without(a, "thermal", "transient"), "thermal"),
        (_with(a, 0, "thermal", "transient", "n_steps"), "thermal.transient.n_steps"),
        (_with(a, 1.5, "thermal", "transient", "theta"), "thermal.transient.theta"),
        (
            _with(a, 1, "thermal", "transient", "startup_half_steps"),
            "thermal.transient.startup_half_steps",
        ),
        (_with(a, 25.0, "thermal", "transient", "t_init_c"), "thermal.transient.t_init_c"),
        (_with(a, 0.0, "thermal", "duty", "duty"), "thermal.duty.duty"),
        (_with(a, 1.5, "thermal", "duty", "duty"), "thermal.duty.duty"),
        (_with(a, -1.0, "thermal", "duty", "period_s"), "thermal.duty.period_s"),
        (_without(a, "thermal", "duty", "duty"), "thermal.duty.duty"),
        (
            {**a, "outputs": {"probes": []}},  # §9 item 545 (ruling 2026-10-19 03:00 Q1)
            "outputs.probes",
        ),
        ({**a, "outputs": {"sar": {"mass_g": 10.0, "centres": [[0.0, 0.0, 0.0]]}}}, "outputs.sar"),
        (_with(sph, {"period_s": 1200.0, "duty": 0.5}, "thermal", "duty"), "thermal.duty"),
        (_with(sph, copy.deepcopy(tr_a), "thermal", "transient"), "thermal.transient"),
        (_with(sph, "uniform_heat_source", "excitation", "kind"), "excitation.kind"),
    ]
    got = [_field(lambda m=m: _study_from_mapping(m)) for m, _ in rows]
    for (m, want_f), g in zip(rows, got):
        _say(capsys, f"  (X) want {want_f!r} got {g!r}")
    assert got == [f for _, f in rows], [(f, g) for (_, f), g in zip(rows, got) if f != g]

    # (c0) dropped keys, at parse
    c1 = _write_yaml(_without(a, "thermal", "duty"), LOGS / "thermal_ball_c0", "c1.yaml", comm)
    th_c1 = _load(c1)[3]
    d_wall = _spec_diff(load_study(WALL), _study_from_mapping(_without(w, "geometry", "wall")))
    d_ab = _spec_diff(load_study(DUTY_A), load_study(DUTY_B))
    th_a, th_b = _load(DUTY_A)[3], _load(DUTY_B)[3]
    ab = []
    _cmp(th_a.transient, th_b.transient, "transient", ab)
    _cmp(th_a.duty, th_b.duty, "duty", ab)
    _say(
        capsys,
        f"  (c0) c1 duty {th_c1.duty!r}; wall dropped diff {d_wall}; A vs B spec diff "
        f"{d_ab}, thermal diff {ab}; build calls {calls['n']}",
    )
    assert th_c1.duty is None and th_a.duty == DutyCycle(*CASES["A"])
    assert d_wall == ["geometry.wall"]
    assert d_ab == []
    assert ab == ["transient.dt_s", "transient.n_steps", "duty.period_s", "duty.duty"]
    assert calls["n"] == 0


# ------------------------------------------------------------------ solves


def _solve_preamble():
    if not _knob():
        pytest.skip(
            f"{S14_ENV} unset: the thermal-ball solve legs are environment-selected "
            "(real build, -n 2)"
        )
    comm = MPI.COMM_WORLD
    assert not IS_COMPLEX, "the records are real-build records"
    assert comm.size == 2, comm.size
    return comm


def _check_summary(run, comm, label):
    """(I=) the summary is its own thermal result (hex), and the JSON is the summary."""
    th, s = run.thermal, run.summary
    T = th.temperature
    n = T.function_space.dofmap.index_map.size_local * T.function_space.dofmap.index_map_bs
    loc = np.real(np.asarray(T.x.array[:n])) - T_INF
    tmax = float(comm.allreduce(float(np.max(loc)) if loc.size else -np.inf, op=MPI.MAX))
    st = s["temperature"]
    js = json.loads(run.paths["summary"].read_text())
    # ``elapsed_s`` is each rank's own wall clock and only rank 0's summary is
    # written (every route does this), so it is compared on rank 0 only; every
    # other key on every rank.  Measured: the whole-dict compare on every rank
    # read False on rank 1 by elapsed_s alone (20261017T004010Z_WF-1-s14-W1a-duty.log:108).
    rest = {k: v for k, v in s.items() if k != "elapsed_s"}
    js_rest = {k: v for k, v in js.items() if k != "elapsed_s"}
    same_json = _all(
        comm, _hex_tree_equal(js_rest, rest) and (comm.rank != 0 or _hex_tree_equal(js, s))
    )
    _p0(
        f"  (I=) {label}: theta_max {st['theta_max_k']!r} vs recomputed {tmax!r}; balance "
        f"{st['max_balance_residual']!r}; power {st['source_power_w']!r}; JSON == summary "
        f"{same_json}; keys {sorted(s)}"
    )
    assert st["theta_max_k"].hex() == tmax.hex()
    assert st["max_balance_residual"].hex() == float(th.max_balance_residual).hex()
    assert st["source_power_w"].hex() == float(np.real(th.source_power)).hex()
    assert same_json
    assert "elapsed_s" in js
    assert s["family"] == "thermal_ball"
    for key in ("frequency_hz", "degree", "touchstone"):
        assert key not in s
    assert run.result is None and run.study.result is None
    assert st["max_balance_residual"] <= BALANCE_BAND, f"(G) balance {label}"


def _owned_hex_equal(Ta, Tb, comm):
    n = Tb.function_space.dofmap.index_map.size_local
    return _all(comm, np.array_equal(np.asarray(Ta.x.array[:n]), np.asarray(Tb.x.array[:n])))


def test_wf1_s14_duty_reproduces_mat17():
    comm = _solve_preamble()
    t00 = time.perf_counter()
    rows = {}
    for lab, f, P, D, dt, t_end, on_time, rec, prior in (
        ("A", DUTY_A, P_A, D_A, 150.0, T_END_A, ON_TIME_A, MAT17_A_PRINT, WF15_A_READ),
        ("B", DUTY_B, P_B, D_B, DT_B, T_END_B, ON_TIME_B, MAT17_B_PRINT, WF15_B_READ),
    ):
        n = _n(t_end, dt)
        out = LOGS / f"thermal_ball_{lab.lower()}"
        _fresh_dir(out, comm)
        t0 = time.perf_counter()
        run = run_study_file(f, output_dir=out, comm=comm)
        t_file = time.perf_counter() - t0
        py = _duty_study(dt, n, DutyCycle(P, D))
        th = run.thermal
        # (UT) + (S0)
        ca, cb, geo = _same_mesh(th.mesh, py.mesh, comm)
        mult = np.array_equal(np.asarray(th.source_multipliers), np.asarray(py.source_multipliers))
        mult_s = run.summary["temperature"]["source_multipliers"] == [
            float(x) for x in th.source_multipliers
        ]
        _census(th, dt, 2, on_time, f"WF1S14 (S0) case {lab} file")
        # (I-Py)
        T_py = _on(th.mesh, py.temperature)
        ip = _theta_owned_rel(T_py, th.temperature, comm)
        hx = _owned_hex_equal(T_py, th.temperature, comm)
        # (G) + (R)
        ex = _exact_duty(th.mesh, t_end, P, D)
        e = _interior_error(th.mesh, th.temperature, ex)
        rows[lab] = (th.mesh, ex)
        _p0(
            f"WF1S14 (UT) case {lab}: summary cells {run.summary['cells']}, cells {ca} / {cb} "
            f"(WF-15 {WF15_CELLS}); geometry equal {geo}; multipliers equal {mult} "
            f"(summary {mult_s}, {len(th.source_multipliers)} (half-)steps)"
        )
        _p0(
            f"WF1S14 (I-Py) case {lab}: owned max|dθ|/max|θ| = {ip:.3e} (band {L0_BAND:.0e}, "
            f"predicted 0 to 2e-14); owned arrays bitwise equal {hx} (printed, not asserted)"
        )
        _p0(
            f"WF1S14 (G) case {lab} dt={dt!r} N={n}: interior_rel_L2={e:.10e} (band "
            f"{INTERIOR_BAND:.0e}) WF-15 read {prior:.10e}; max_balance="
            f"{run.summary['temperature']['max_balance_residual']:.3e}; file run {t_file:.1f}s"
        )
        _p0(f"WF1S14 (R) case {lab}: {e:.4e} vs record {rec}")
        _check_summary(run, comm, f"case {lab}")
        assert run.summary["cells"] == ca == cb, f"(UT) cells case {lab}"
        assert geo, f"(UT) geometry case {lab}"
        assert mult and mult_s, f"(UT) multipliers case {lab}"
        assert run.summary["temperature"]["duty"] == {"period_s": P, "duty": D}
        assert ip <= L0_BAND, f"(I-Py) case {lab} {ip:.3e}"
        assert e <= INTERIOR_BAND, f"(G) case {lab} {e:.4e}"
        assert f"{e:.4e}" == rec, f"(R) case {lab} {e:.4e} != {rec}"

    # (c1) the schedule dropped through a file, case A.
    c1f = _write_yaml(
        _without(_mapping(DUTY_A), "thermal", "duty"), LOGS / "thermal_ball_c1", "c1.yaml", comm
    )
    out = LOGS / "thermal_ball_c1" / "out"
    c1 = run_study_file(c1f, output_dir=out, comm=comm)
    msh_a, ex_a = rows["A"]
    e_c1 = _interior_error(msh_a, _on(msh_a, c1.thermal.temperature), ex_a)
    _p0(
        f"WF1S14 (c1) duty dropped through a file, case A: interior miss {e_c1:.4e} vs WF-15 "
        f"API {WF15_C1_READ:.4e} (asserted > {INTERIOR_BAND:.0e}); summary duty "
        f"{c1.summary['temperature']['duty']!r}, multipliers "
        f"{c1.summary['temperature']['source_multipliers']!r}"
    )
    assert c1.summary["temperature"]["duty"] is None
    assert c1.summary["temperature"]["source_multipliers"] == []
    assert e_c1 > INTERIOR_BAND, f"(c1) did not miss: {e_c1:.4e}"
    _p0(f"WF1S14 W1a leg {time.perf_counter() - t00:.1f}s")


def test_wf1_s14_wall_reproduces_mat18_2b():
    comm = _solve_preamble()
    t00 = time.perf_counter()
    msh, ct, ft, nglob, k_f, rc_f, q_f, theta = _fine()
    assert nglob <= CELL_GUARD  # before any solve
    out = LOGS / "thermal_ball_wall"
    _fresh_dir(out, comm)
    t0 = time.perf_counter()
    run = run_study_file(WALL, output_dir=out, comm=comm)
    t_file = time.perf_counter() - t0
    py = _wall_study(RC2)
    th = run.thermal
    # (UT)
    ca, cb, geo = _same_mesh(th.mesh, py.mesh, comm)
    _, cf, geo_f = _same_mesh(th.mesh, msh, comm)
    dk = _dg0_owned_diff(th.k, k_f, comm)
    drc = _dg0_owned_diff(th.rho_c, rc_f, comm)
    dq = _dg0_owned_diff(th.q, q_f, comm)
    _p0(
        f"WF1S14 (UT) wall: summary cells {run.summary['cells']}, cells {ca} / {cb} / fine {cf} "
        f"(record {MAT18_T1_CELLS}); geometry equal py {geo} fine {geo_f}; DG0 k max|d| "
        f"{dk[0]!r}, rho_c {drc[0]!r}, q {dq[0]!r} of {dk[2]} owned"
    )
    # (I-Py)
    T_py = _on(th.mesh, py.temperature)
    ip = _theta_owned_rel(T_py, th.temperature, comm)
    hx = _owned_hex_equal(T_py, th.temperature, comm)
    # (G) + (R)
    e_all, e_in, e_sh = _errors(msh, ct, _on(msh, th.temperature), theta)
    _p0(
        f"WF1S14 (I-Py) wall: owned max|dθ|/max|θ| = {ip:.3e} (band {L0_BAND:.0e}); owned arrays "
        f"bitwise equal {hx} (printed, not asserted)"
    )
    _p0(
        f"WF1S14 (G) (T1) N={N_CN} res={RES}: interior_rel_L2={e_in:.10e} (band "
        f"{INTERIOR_BAND:.0e}) whole={e_all:.4e} shell={e_sh:.4e}; WF-15 read {WF15_T1_READ:.10e}; "
        f"max_balance={run.summary['temperature']['max_balance_residual']:.3e}; file run "
        f"{t_file:.1f}s"
    )
    _p0(f"WF1S14 (R) (T1): {e_in:.4e} vs record {MAT18_T1_PRINT}")
    _check_summary(run, comm, "wall")

    # (c2) the wall's rho_c dropped through a file (RC1 in the wall).
    c2f = _write_yaml(
        _with(_mapping(WALL), 4180000.0, "geometry", "wall", "rho_c_j_per_m3_k"),
        LOGS / "thermal_ball_c2",
        "c2.yaml",
        comm,
    )
    c2 = run_study_file(c2f, output_dir=LOGS / "thermal_ball_c2" / "out", comm=comm)
    _, e_c2, _ = _errors(msh, ct, _on(msh, c2.thermal.temperature), theta)
    _p0(
        f"WF1S14 (c2) wall rho_c = RC1 through a file: interior miss {e_c2:.4e} vs WF-15 API "
        f"{WF15_C2_READ:.4e} and MAT-18 prior {MAT18_C1_UNIFORM_RHOC:.4e} (asserted > "
        f"{INTERIOR_BAND:.0e}); max_balance={c2.summary['temperature']['max_balance_residual']:.3e}"
    )
    _p0(f"WF1S14 W1b leg {time.perf_counter() - t00:.1f}s")

    assert run.summary["cells"] == ca == cb == cf == MAT18_T1_CELLS, "(UT) cells"
    assert geo and geo_f, "(UT) geometry"
    assert dk[0] == 0.0 and drc[0] == 0.0 and dq[0] == 0.0, "(UT) DG0 coefficients"
    assert run.summary["temperature"]["wall"] == {
        "thickness_m": B_ - A_,
        "k_w_per_m_k": K2,
        "rho_c_j_per_m3_k": RC2,
    }
    assert run.summary["temperature"]["duty"] is None
    assert ip <= L0_BAND, f"(I-Py) {ip:.3e}"
    assert e_in <= INTERIOR_BAND, f"(G) (T1) {e_in:.4e}"
    assert f"{e_in:.4e}" == MAT18_T1_PRINT, f"(R) {e_in:.4e} != {MAT18_T1_PRINT}"
    assert c2.summary["temperature"]["max_balance_residual"] <= BALANCE_BAND
    assert e_c2 > INTERIOR_BAND, f"(c2) did not miss: {e_c2:.4e}"
