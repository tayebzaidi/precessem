"""`WF-1` step 15: ``femem run`` on thermal-only ``thermal_box`` files reproduces
`WF-15` step 5's (§9 item 521a) `MAT-21` records.

Brief: docs/planning/chunks/WF-1.md, "§9 item 521" / "521b" and the review
ruling 2026-10-18 03:00.  The solve legs are gated in the complex build (the
build of the `MAT-21` records); they run only with ``FEM_EM_WF1_S15=1`` (and
then assert the complex build and ``-n 2``, never skip).  The parse node builds
no mesh and runs in either build.

Anchors: (K) the YAML Study = the Python Study (``_spec_diff``) and the parsed
``RobinCooling`` / ``TransientSpec`` / ``DutyCycle`` by ``float.hex``; (X)
refusals naming their field; (UT) the file's mesh = the Python Study's, and
multipliers equal; (I-Py) the file's theta = the Python Study's at
``L0_BAND``; (R) the records at ``.4e`` with ``INTERIOR_BAND`` / ``E_BAND`` /
``S2_GAUSS_BAND``; (I=) the summary = its own thermal result.
Controls (rule (e)): (c-wall) the h-doubled tracked file, asserted, backed by
521a's (c-wall) 3.8048e-01 (``20261018T112154Z_WF-15-s5-W1-steady.log:70``);
(c3) through the duty file,
asserted, backed by P26's 5.0282e-01; (c0) the wall rides in the cooling.
"""

from __future__ import annotations

import copy
import json
import os
import time

import numpy as np
import pytest
import ufl
from dolfinx import default_scalar_type, fem
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
from tests.package.test_wf15_study_thermal_box import (
    P26_C3,
    P26_T_DUTY,
    STEP1_M4_CELLS,
    STEP1_M4_WHOLE,
    STEP1_M6_CELLS,
)
from tests.package.test_wf15_study_transient import _on, _same_mesh, _theta_owned_rel
from tests.validation.test_mat5_em_coupled_sphere import L0_BAND
from tests.validation.test_mat10_sar_to_heat_sphere import ROBIN_TAG, _allsum
from tests.validation.test_mat12_transient_sar_heating_sphere import E_BAND, RHO_C
from tests.validation.test_mat13_transient_perfusion_sphere import INTERIOR_BAND
from tests.validation.test_mat15_coil_driven_heating_phantom import H_GEL, K_GEL, N_MAIN, T_INF
from tests.validation.test_mat15_finite_cylinder_fem import S2_GAUSS_BAND, _rel
from tests.validation.test_mat21_box_closed_form import HALF, Q
from tests.validation.test_mat21_box_duty_closed_form import DT_BOX, HALF_T, T_STAR, _duty_cf
from tests.validation.test_mat21_box_fem import _exact, _interior_chi, _steady_ex, _trans_ex

S15_ENV = "FEM_EM_WF1_S15"
STEADY = WF1 / "thermal_box_steady_m4.yaml"
DUTY = WF1 / "thermal_box_duty_m6.yaml"
WRONG_WALL = WF1 / "thermal_box_wrong_wall_m4.yaml"
MALFORMED_DUTY = WF1 / "malformed_thermal_box_duty_without_transient.yaml"
BALL = WF1 / "thermal_ball_duty_a.yaml"
DIV4 = (52, 32, 8)
DIV6 = (78, 48, 12)

# 521a's (c-wall) on the same fixture through the API
# (20261018T112154Z_WF-15-s5-W1-steady.log:70; printed beside, the assertion
# is > INTERIOR_BAND).
WF15_S5_C_WALL = 3.8048e-01


def _knob() -> bool:
    raw = os.environ.get(S15_ENV, "").strip()
    if raw == "":
        return False
    if raw == "1":
        return True
    raise ValueError(f"{S15_ENV}={raw!r}: only unset/empty or '1' are accepted")


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


def _hex_tree_equal(a, b):
    if isinstance(a, float) and isinstance(b, float):
        return a.hex() == b.hex()
    if isinstance(a, dict) and isinstance(b, dict):
        return a.keys() == b.keys() and all(_hex_tree_equal(a[k], b[k]) for k in a)
    if isinstance(a, list) and isinstance(b, list):
        return len(a) == len(b) and all(_hex_tree_equal(x, y) for x, y in zip(a, b))
    return type(a) is type(b) and a == b


def _py_study(divs):
    return Study(ThermalBoxGeometry(HALF, divs), None, UniformHeatSource(Q))


# ------------------------------------------------------------------ no mesh


def test_wf1_s15_parse_and_refusals(monkeypatch, capsys):
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
        ThermalBoxGeometry,
    ):
        monkeypatch.setattr(cls, "build", stub)

    # (K) the Study spec
    want_spec = {STEADY: _py_study(DIV4), WRONG_WALL: _py_study(DIV4), DUTY: _py_study(DIV6)}
    diffs = {f.name: _spec_diff(load_study(f), s) for f, s in want_spec.items()}
    _say(capsys, f"\n[WF-1 s15] (K) spec diffs {diffs}")
    assert all(d == [] for d in diffs.values()), diffs
    for f in want_spec:
        st = load_study(f)
        assert st.frequency_hz is None and st.material is None
        assert all(type(n) is int for n in st.geometry.divisions)

    # (K) the fourth element
    want = {
        STEADY: (RobinCooling(K_GEL, H_GEL, T_INF), None, None),
        WRONG_WALL: (RobinCooling(K_GEL, 2.0 * H_GEL, T_INF), None, None),
        DUTY: (
            RobinCooling(K_GEL, H_GEL, T_INF),
            TransientSpec(DT_BOX, N_MAIN, RHO_C),
            DutyCycle(2.0 * HALF_T, 0.5),
        ),
    }
    for f, (cool, tr, du) in want.items():
        _, kind, outputs, th = _load(f)
        out = []
        _cmp(th.cooling, cool, "cooling", out)
        _cmp(th.transient, tr, "transient", out)
        _cmp(th.duty, du, "duty", out)
        _say(capsys, f"  (K) {f.name}: kind {kind!r}, outputs {outputs!r}, thermal diffs {out}")
        assert kind == "thermal_box"
        assert outputs == {}
        assert out == [], (f.name, out)
    assert _load(STEADY)[3].transient is None and _load(STEADY)[3].duty is None
    th_d = _load(DUTY)[3]
    _say(
        capsys,
        f"  (K) duty dt_s {th_d.transient.dt_s.hex()} (DT_BOX {DT_BOX.hex()}); period "
        f"{th_d.duty.period_s.hex()} (2 HALF_T {(2.0 * HALF_T).hex()})",
    )

    # (K) femem check
    for f in (STEADY, WRONG_WALL, DUTY):
        capsys.readouterr()
        rc = cli.main(["check", str(f)])
        lines = capsys.readouterr().out.splitlines()
        rc_all = _rc_all(rc, comm)
        _say(capsys, f"  (K) femem check {f.name}: rc {rc_all}; rank-0 lines {lines}")
        assert rc_all == (0, 0)
        if comm.rank == 0:
            assert [ln for ln in lines if ln.startswith("family: ")] == ["family: thermal_box"]
            assert "excitation: UniformHeatSource" in lines

    # (X) the tracked malformed file, load + cli run / check
    assert _field(lambda: load_study(MALFORMED_DUTY)) == "thermal.duty"
    d = LOGS / "thermal_box_malformed"
    _fresh_dir(d, comm)
    capsys.readouterr()
    rc_run = cli.main(["run", str(MALFORMED_DUTY), "--output-dir", str(d)])
    rc_check = cli.main(["check", str(MALFORMED_DUTY)])
    capsys.readouterr()
    run_all, check_all = _rc_all(rc_run, comm), _rc_all(rc_check, comm)
    exists = comm.allreduce(bool(d.exists()), op=MPI.LOR)
    _say(
        capsys,
        f"  (X) malformed duty without transient: run rc {run_all}, check rc "
        f"{check_all}, out exists {exists}",
    )
    assert run_all == (2, 2)
    assert check_all == (2, 2)
    assert not exists

    # (X) the variant table
    s = _mapping(STEADY)
    du = _mapping(DUTY)
    rows = [
        ({**s, "study": {"frequency_hz": 64000000.0}}, "study"),
        ({**s, "material": {"sigma": 0.5, "epsilon_r": 78.0}}, "material"),
        ({**s, "conductor": {"kind": "leontovich"}}, "conductor"),
        ({**s, "tune": {"k": 1}}, "tune"),
        (_with(s, [0.325, 0.21], "geometry", "half_widths_m"), "geometry.half_widths_m"),
        (_with(s, [0.325, 0.0, 0.045], "geometry", "half_widths_m"), "geometry.half_widths_m"),
        (_with(s, [0.325, "a", 0.045], "geometry", "half_widths_m"), "geometry.half_widths_m"),
        (_with(s, [52, 0, 8], "geometry", "divisions"), "geometry.divisions"),
        (_with(s, [52, 32.0, 8], "geometry", "divisions"), "geometry.divisions"),
        (
            _with(
                s,
                {"thickness_m": 0.005, "k_w_per_m_k": 0.19, "rho_c_j_per_m3_k": 1.0e6},
                "geometry",
                "wall",
            ),
            "geometry.wall",
        ),
        (_with(s, 0.1, "geometry", "radius"), "geometry.radius"),
        (_with(s, 0.0, "excitation", "q_w_per_m3"), "excitation.q_w_per_m3"),
        (_with(s, 2.0, "thermal", "source_scale"), "thermal.source_scale"),
        (_with(du, 0, "thermal", "transient", "n_steps"), "thermal.transient.n_steps"),
        (_with(du, 1.5, "thermal", "duty", "duty"), "thermal.duty.duty"),
        # §9 item 539 (review ruling 2026-10-18 18:00): a valid probe is admitted
        # on a box since step 17; the row's input moved to a still-refused form.
        ({**s, "outputs": {"probes": []}}, "outputs.probes"),
    ]
    got = [_field(lambda m=m: _study_from_mapping(m)) for m, _ in rows]
    for (m, want_f), g in zip(rows, got):
        _say(capsys, f"  (X) want {want_f!r} got {g!r}")
    assert got == [f for _, f in rows], [(f, g) for (_, f), g in zip(rows, got) if f != g]
    # The ball still refuses a missing transient (item 500's pinned row, beside the box).
    assert (
        _field(lambda: _study_from_mapping(_without(_mapping(BALL), "thermal", "transient")))
        == "thermal"
    )
    assert _field(lambda: _study_from_mapping(_without(du, "thermal", "duty"))) is None

    # (c0) by construction: the wall rides in the cooling, not the Study
    d_ww = _spec_diff(load_study(STEADY), load_study(WRONG_WALL))
    cw = []
    _cmp(_load(STEADY)[3].cooling, _load(WRONG_WALL)[3].cooling, "cooling", cw)
    _say(
        capsys,
        f"  (c0) steady vs wrong wall: spec diff {d_ww}; cooling diff {cw}; build "
        f"calls {calls['n']}",
    )
    assert d_ww == []
    assert cw == ["cooling.h_w_per_m2_k"]
    assert calls["n"] == 0


# ------------------------------------------------------------------ solves


def _solve_preamble():
    if not _knob():
        pytest.skip(
            f"{S15_ENV} unset: the thermal-box solve legs are environment-selected "
            "(complex build, -n 2)"
        )
    comm = MPI.COMM_WORLD
    assert IS_COMPLEX, "the MAT-21 records are complex-build records"
    assert comm.size == 2, comm.size
    return comm


def _check_summary(run, comm, label, transient):
    """(I=) the summary is its own thermal result (hex), and the JSON is the summary."""
    th, s = run.thermal, run.summary
    T = th.temperature
    n = T.function_space.dofmap.index_map.size_local * T.function_space.dofmap.index_map_bs
    loc = np.real(np.asarray(T.x.array[:n])) - T_INF
    tmax = float(comm.allreduce(float(np.max(loc)) if loc.size else -np.inf, op=MPI.MAX))
    st = s["temperature"]
    js = json.loads(run.paths["summary"].read_text())
    # elapsed_s compared on rank 0 only (item 500's measured fix,
    # 20261017T004010Z_WF-1-s14-W1a-duty.log:108).
    rest = {k: v for k, v in s.items() if k != "elapsed_s"}
    js_rest = {k: v for k, v in js.items() if k != "elapsed_s"}
    same_json = _all(
        comm, _hex_tree_equal(js_rest, rest) and (comm.rank != 0 or _hex_tree_equal(js, s))
    )
    _p0(
        f"  (I=) {label}: theta_max {st['theta_max_k']!r} vs recomputed {tmax!r}; power "
        f"{st['source_power_w']!r}; JSON == summary {same_json}; keys {sorted(s)}; "
        f"temperature keys {sorted(st)}"
    )
    assert st["theta_max_k"].hex() == tmax.hex()
    assert st["source_power_w"].hex() == float(np.real(th.source_power)).hex()
    assert same_json
    assert "elapsed_s" in js
    assert s["family"] == "thermal_box"
    for key in ("frequency_hz", "degree", "touchstone"):
        assert key not in s
    assert run.result is None and run.study.result is None
    if transient:
        assert st["max_balance_residual"].hex() == float(th.max_balance_residual).hex()
        assert st["source_multipliers"] == [float(m) for m in th.source_multipliers]
        assert "convected_power_w" not in st
    else:
        assert st["transient"] is None and st["duty"] is None
        assert st["convected_power_w"].hex() == float(np.real(th.convected_power)).hex()
        assert "max_balance_residual" not in st and "source_multipliers" not in st


def test_wf1_s15_steady_reproduces_mat21_m4():
    comm = _solve_preamble()
    t00 = time.perf_counter()
    out = LOGS / "thermal_box_steady"
    _fresh_dir(out, comm)
    run = run_study_file(STEADY, output_dir=out, comm=comm)
    t_file = time.perf_counter() - t00
    th = run.thermal
    py = _py_study(DIV4).temperature(RobinCooling(K_GEL, H_GEL, T_INF))
    msh = th.submesh
    # (UT)
    ca, cb, geo = _same_mesh(msh, py.submesh, comm)
    # (I-Py)
    T_py = _on(msh, py.temperature)
    ip = _theta_owned_rel(T_py, th.temperature, comm)
    n_own = th.temperature.function_space.dofmap.index_map.size_local
    hx = _all(
        comm,
        np.array_equal(
            np.asarray(T_py.x.array[:n_own]), np.asarray(th.temperature.x.array[:n_own])
        ),
    )
    # (R)
    th_h = th.temperature - T_INF
    ex, npts = _exact(msh, _steady_ex())
    whole = _rel(msh, th_h, ex)
    chi = _interior_chi(msh)
    interior = _rel(msh, th_h, ex, chi)
    ds = ufl.Measure("ds", domain=msh, subdomain_data=th.facet_tags)(ROBIN_TAG)
    conv = float(np.real(_allsum(comm, H_GEL * th_h * ds)))
    one = fem.Constant(msh, default_scalar_type(1.0))
    vh = float(np.real(_allsum(comm, one * ufl.dx(domain=msh))))
    gauss = abs(conv - Q * vh) / (Q * vh)
    # (c-wall) through the tracked file
    out_w = LOGS / "thermal_box_wrong_wall"
    _fresh_dir(out_w, comm)
    run_w = run_study_file(WRONG_WALL, output_dir=out_w, comm=comm)
    ca_w, _cb_w, geo_w = _same_mesh(run_w.thermal.submesh, msh, comm)
    c_wall = _rel(msh, _on(msh, run_w.thermal.temperature) - T_INF, ex, chi)
    _p0(
        f"WF1S15 (UT) steady: summary cells {run.summary['cells']}, cells {ca} / {cb} (record "
        f"{STEP1_M4_CELLS}); geometry equal {geo}; file run {t_file:.1f}s"
    )
    _p0(
        f"WF1S15 (I-Py) steady: owned max|dθ|/max|θ| = {ip:.3e} (hex {float(ip).hex()}; band "
        f"{L0_BAND:.0e}); owned arrays bitwise equal {hx} (printed, not asserted)"
    )
    _p0(
        f"WF1S15 (R) steady: whole L2 rel = {whole:.6e} -> {whole:.4e} (record {STEP1_M4_WHOLE}); "
        f"interior {interior:.6e} (band {INTERIOR_BAND:.0e}); Gauss {gauss:.3e} (band "
        f"{S2_GAUSS_BAND:.0e}); {npts} P2 points"
    )
    _p0(
        f"WF1S15 control (c-wall) wrong-wall file: interior miss {c_wall:.4e} (521a API "
        f"{WF15_S5_C_WALL:.4e}; asserted > {INTERIOR_BAND:.0e}); mesh equal {geo_w}, cells {ca_w}"
    )
    _check_summary(run, comm, "steady", transient=False)
    assert run.summary["cells"] == ca == cb == STEP1_M4_CELLS
    assert geo and geo_w
    assert ip <= L0_BAND, f"(I-Py) steady {ip:.3e}"
    assert f"{whole:.4e}" == STEP1_M4_WHOLE
    assert interior <= INTERIOR_BAND
    assert gauss <= S2_GAUSS_BAND
    assert run_w.summary["temperature"]["transient"] is None
    assert c_wall > INTERIOR_BAND, f"(c-wall) did not miss: {c_wall:.4e}"
    _p0(f"WF1S15 W1a leg {time.perf_counter() - t00:.1f}s")


def test_wf1_s15_duty_reproduces_p26():
    comm = _solve_preamble()
    t00 = time.perf_counter()
    out = LOGS / "thermal_box_duty"
    _fresh_dir(out, comm)
    run = run_study_file(DUTY, output_dir=out, comm=comm)
    t_file = time.perf_counter() - t00
    th = run.thermal
    py = _py_study(DIV6).temperature(
        RobinCooling(K_GEL, H_GEL, T_INF),
        transient=TransientSpec(DT_BOX, N_MAIN, RHO_C),
        duty=DutyCycle(period_s=2.0 * HALF_T, duty=0.5),
    )
    msh = th.mesh
    ca, cb, geo = _same_mesh(msh, py.mesh, comm)
    mult = np.array_equal(np.asarray(th.source_multipliers), np.asarray(py.source_multipliers))
    T_py = _on(msh, py.temperature)
    ip = _theta_owned_rel(T_py, th.temperature, comm)
    th_h = th.temperature - T_INF
    chi = _interior_chi(msh)
    ex, npts = _exact(msh, lambda x, y, z: _duty_cf(x, y, z, T_STAR))
    t_duty = _rel(msh, th_h, ex, chi)
    ex_s, _ = _exact(msh, _trans_ex(T_STAR))
    c3 = _rel(msh, th_h, ex_s, chi)
    bal = run.summary["temperature"]["max_balance_residual"]
    _p0(
        f"WF1S15 (UT) duty: summary cells {run.summary['cells']}, cells {ca} / {cb} (record "
        f"{STEP1_M6_CELLS}); geometry equal {geo}; multipliers equal {mult} "
        f"({len(th.source_multipliers)} entries); file run {t_file:.1f}s"
    )
    _p0(
        f"WF1S15 (I-Py) duty: owned max|dθ|/max|θ| = {ip:.3e} (hex {float(ip).hex()}; band "
        f"{L0_BAND:.0e})"
    )
    _p0(
        f"WF1S15 (R) duty: (T-duty) interior rel = {t_duty:.6e} -> {t_duty:.4e} (record "
        f"{P26_T_DUTY}; band {INTERIOR_BAND:.0e}); balance {bal:.3e} (band {E_BAND:.0e}); "
        f"{npts} P2 points"
    )
    _p0(
        f"WF1S15 control (c3) step response vs duty file: miss {c3:.4e} (P26 {P26_C3:.4e}; "
        f"asserted > {INTERIOR_BAND:.0e})"
    )
    _check_summary(run, comm, "duty", transient=True)
    assert run.summary["cells"] == ca == cb == STEP1_M6_CELLS
    assert geo
    assert mult
    assert run.summary["temperature"]["duty"] == {"period_s": 2.0 * HALF_T, "duty": 0.5}
    assert ip <= L0_BAND, f"(I-Py) duty {ip:.3e}"
    assert t_duty <= INTERIOR_BAND
    assert f"{t_duty:.4e}" == P26_T_DUTY
    assert bal <= E_BAND
    assert c3 > INTERIOR_BAND
    _p0(f"WF1S15 W1b leg {time.perf_counter() - t00:.1f}s")
