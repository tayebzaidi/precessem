"""`WF-1` step 29 (§9 item 873, §10 U134): the YAML twin of `WF-15` step 15's `MAT-12` transient.

The EM-sphere study file gains ``thermal: {transient: {t_end_s, n_steps,
rho_c_j_per_m3_k, theta, startup_half_steps}}`` (dt = t_end_s / n_steps), and
``femem run`` reproduces U131's (A_t) / (S2) records on the same mesh and Δt
(brief: docs/planning/chunks/WF-1.md, "§9 item 873").

U131's records (`WF-15` step 15, docs/testing/logs/20261107T013321Z_WF-15-s15-W1.log):
``:138`` (A_t) interior rel L2 = 3.4708293773e-02 [3.4708e-02]; ``:139`` (S2)
ratio 0.999950006 (|ratio − 1| = 4.9994e-05); ``:141`` (c-qs) 1.0108e+00;
``:142`` EM 15.8 s.

Gates (bands imported, none set):
(K)    the parsed Study spec equals ``_study(degree=1)``'s; the parsed
       ``RobinCooling`` / ``TransientSpec`` ``repr``s equal U131's; dt by
       ``float.hex``; ``femem check`` rc (0, 0);
(Ref)  ``thermal.duty`` / ``thermal.initial`` / ``outputs.probes`` refused with
       the API's own ``str(exc)``; a four-port file naming a transient refused at
       ``thermal`` with the landed reason; ``dt_s`` refused at
       ``thermal.transient`` with ``_SPHERE_TRANSIENT_DT_REFUSAL`` (item 892); no build, ``result`` None;
(UT), (Once), (A_t) / (R-A_t), (S2) / (R-S2), (E), (IT), (CLI = API).
Control (asserted, backed): (c-qs) the steady θ for θ(t*) misses ``A_BAND``.

Solve node runs with ``FEM_EM_WF1_S29=1`` (complex build, ``-n 2``).
"""

from __future__ import annotations

import json
import os
import time

import numpy as np
import pytest
import ufl
from dolfinx import default_scalar_type, fem
from mpi4py import MPI

from fem_em_solver import (
    FourPortBirdcageGeometry,
    RobinCooling,
    SphereInBoxGeometry,
    StudyConfigError,
    TransientSpec,
    cli,
    load_study,
    run_study_file,
)
from fem_em_solver.post.sar import uniform_sphere_sar_closed_form
from fem_em_solver.study_yaml import (
    _SPHERE_DUTY_REFUSAL,
    _SPHERE_TRANSIENT_DT_REFUSAL,
    _BallThermal,
    _parse_mapping,
)
from fem_em_solver.thermal.probe import ProbeBall
from fem_em_solver.utils.constants import EPSILON_0
from tests.complex_mode import complex_only
from tests.package.test_wf1_study_yaml import (
    LOGS,
    WF1,
    _fresh_dir,
    _mapping,
    _rc_all,
    _spec_diff,
)
from tests.package.test_wf1_yaml_thermal_box_wall import _field, _with, _without
from tests.package.test_wf15_study_lossy_sphere import _study
from tests.package.test_wf15_study_sphere_transient import (
    MAT12_A_T_PRINT,
    MAT12_S2_PRINT,
    _owned_rel,
)
from tests.package.test_wf15_study_temperature import MAT10_SUB_CELLS
from tests.validation.test_lossy_sphere_sar import (
    E0,
    EPSILON_R_SPHERE,
    OMEGA,
    RHO_KG_M3,
    SIGMA_HIGH,
)
from tests.validation.test_mat5_em_coupled_sphere import L0_BAND
from tests.validation.test_mat10_sar_to_heat_sphere import (
    A_BAND,
    G,
    H,
    K,
    _allsum,
    _interior_indicator,
    _interior_rel,
)
from tests.validation.test_mat12_transient_sar_heating_sphere import (
    DT_1,
    E_BAND,
    N_MAIN,
    RHO_C,
    S2_BAND,
    T_STAR,
    _exact_theta,
    _mean,
)

S29_ENV = "FEM_EM_WF1_S29"
needs_s29 = pytest.mark.skipif(
    os.environ.get(S29_ENV) != "1", reason=f"opt-in solve: set {S29_ENV}=1"
)
A_FILE = WF1 / "sphere_mat12_transient_a.yaml"
S_FILE = WF1 / "sphere_mat12_transient_s.yaml"
DUTY_FILE = WF1 / "malformed_sphere_transient_duty.yaml"
PORT_FILE = WF1 / "malformed_transient_on_port.yaml"
PORT_LANDED = WF1 / "malformed_thermal_on_port.yaml"
# U131's printed (A_t) to ten digits (20261107T013321Z_WF-15-s15-W1.log:138), printed beside ours.
U131_A_T_10E = "3.4708293773e-02"
U131_S2_RATIO = "0.999950006"  # :139


def _spec_a():
    return TransientSpec(T_STAR / N_MAIN, N_MAIN, RHO_C)


def _spec_s():
    return TransientSpec(DT_1, 1, RHO_C, theta=1.0, startup_half_steps=0)


def _p0(*args):
    if MPI.COMM_WORLD.rank == 0:
        print(*args, flush=True)


def _thermal(m, source="<m>"):
    return _parse_mapping(m, source)[3]


def _exc(fn):
    try:
        fn()
    except StudyConfigError as exc:
        return exc
    return None


def test_wf1_s29_parse_and_refusals(monkeypatch):
    """W0: (K) and (Ref); no mesh, no solve."""
    comm = MPI.COMM_WORLD
    calls = {"n": 0}

    def stub(self, comm=MPI.COMM_WORLD):
        calls["n"] += 1
        raise AssertionError("build must not be called by a parse")

    monkeypatch.setattr(SphereInBoxGeometry, "build", stub)
    monkeypatch.setattr(FourPortBirdcageGeometry, "build", stub)

    # (K)
    cool_ref = RobinCooling(K, H, 0.0, source_scale=G)
    for path, spec_ref, dt_ref in ((A_FILE, _spec_a(), T_STAR / N_MAIN), (S_FILE, _spec_s(), DT_1)):
        d_spec = _spec_diff(load_study(path), _study(degree=1))
        th = _thermal(_mapping(path), str(path))
        rc_check = _rc_all(cli.main(["check", str(path)]), comm)
        _p0(
            f"\n[WF-1 s29] (K) {path.name}: spec diff {d_spec}; cooling {th.cooling!r}; "
            f"transient {th.transient!r} (U131 {spec_ref!r}); dt {th.transient.dt_s.hex()} "
            f"(U131 {dt_ref.hex()}); check rc {rc_check}"
        )
        assert isinstance(th, _BallThermal) and th.duty is None
        assert d_spec == [], d_spec
        assert repr(th.cooling) == repr(cool_ref), (th.cooling, cool_ref)
        assert repr(th.transient) == repr(spec_ref), (th.transient, spec_ref)
        assert th.transient.dt_s.hex() == dt_ref.hex()
        assert rc_check == (0, 0)

    # (Ref) the API's own messages, from U131's no-solve calls.
    # The duty reason is YAML-owned since §9 item 891 (the API accepts duty on the
    # sphere since U135, WF-15 step 17): compared with the module constant, never
    # with a live API call.
    study = _study(degree=1)
    _p0(
        f"  (Ref) duty split: the API accepts DutyCycle on the EM sphere since U135 "
        f"(WF-15 step 17); YAML refuses thermal.duty (no sphere duty key commissioned): "
        f"{_SPHERE_DUTY_REFUSAL!r}"
    )
    with pytest.raises(ValueError) as e_init:
        study.temperature(
            cool_ref,
            transient=TransientSpec(T_STAR / N_MAIN, N_MAIN, RHO_C, initial_gradient_k_per_m=1.0),
        )
    with pytest.raises(ValueError) as e_probe:
        study.probe([ProbeBall((0.0, 0.0, 0.0), 1e-3)], cool_ref, transient=_spec_a())
    landed = _exc(lambda: load_study(PORT_LANDED))

    mA = _mapping(A_FILE)
    cases = {
        "duty file": (_exc(lambda: load_study(DUTY_FILE)), "thermal.duty", _SPHERE_DUTY_REFUSAL),
        "initial": (
            _exc(lambda: _thermal(_with(mA, {"gradient_k_per_m": 1.0}, "thermal", "initial"))),
            "thermal.initial",
            str(e_init.value),
        ),
        "probes": (
            _exc(
                lambda: _thermal(
                    _with(mA, [{"centre": [0.0, 0.0, 0.0], "radius": 1e-3}], "outputs", "probes")
                )
            ),
            "outputs.probes",
            str(e_probe.value),
        ),
        "port file": (_exc(lambda: load_study(PORT_FILE)), "thermal", landed.reason),
    }
    for label, (exc, field, reason) in cases.items():
        _p0(
            f"  (Ref) {label}: field {getattr(exc, 'field', None)!r} (expected {field!r}); "
            f"reason {getattr(exc, 'reason', None)!r}; API {reason!r}"
        )
        assert type(exc) is StudyConfigError, label
        assert exc.field == field, (label, exc.field)
        assert exc.reason == reason, (label, exc.reason, reason)
    dt_m = _with(
        _without(mA, "thermal", "transient", "t_end_s"),
        T_STAR / N_MAIN,
        "thermal",
        "transient",
        "dt_s",
    )
    exc_dt = _exc(lambda: _thermal(dt_m))
    _p0(f"  (Ref) dt_s for t_end_s: field {exc_dt.field!r}; reason {exc_dt.reason!r}")
    # §9 item 892 (18:00 ruling): refused at the block, with the YAML-owned
    # reason, as item 873's brief fixed it (row 322 of the s14 module unedited).
    assert exc_dt.field == "thermal.transient", exc_dt.field
    assert exc_dt.reason == _SPHERE_TRANSIENT_DT_REFUSAL, exc_dt.reason
    assert _field(lambda: _thermal(dt_m)) == "thermal.transient"

    out = LOGS / "malformed_sphere_transient"
    for path in (DUTY_FILE, PORT_FILE):
        _fresh_dir(out, comm)
        run_all = _rc_all(cli.main(["run", str(path), "--output-dir", str(out)]), comm)
        check_all = _rc_all(cli.main(["check", str(path)]), comm)
        exists = comm.allreduce(bool(out.exists()), op=MPI.LOR)
        _p0(f"  (Ref) {path.name}: run rc {run_all}, check rc {check_all}, out exists {exists}")
        assert run_all == (2, 2) and check_all == (2, 2)
        assert not exists

    _p0(f"  (Ref) build calls {calls['n']}; study.result {study.result}")
    assert calls["n"] == 0
    assert study.result is None


def _a_t(th):
    cf = uniform_sphere_sar_closed_form(
        e0=E0,
        epsilon_r=EPSILON_R_SPHERE,
        sigma=SIGMA_HIGH,
        omega=OMEGA,
        rho=RHO_KG_M3,
        epsilon_0=EPSILON_0,
    )
    q_cf = G * RHO_KG_M3 * cf["sar_w_per_kg"]
    sub = th.mesh
    theta = _exact_theta(sub, T_STAR, q_cf)
    return _interior_rel(sub, _interior_indicator(sub), th.temperature - theta, theta), q_cf


def _s2(th, comm):
    subS = th.mesh
    one = fem.Constant(subS, default_scalar_type(1.0))
    vol = float(np.real(_allsum(comm, one * ufl.dx(domain=subS))))
    qbar = _mean(subS, th.q, vol)
    thbar = _mean(subS, th.temperature, vol)
    ratio = RHO_C * thbar / (DT_1 * qbar)
    return abs(ratio - 1.0), ratio


def _close(a, b, path=""):
    """Key sets equal; numbers within L0_BAND relative; returns the mismatching paths."""
    bad = []
    if isinstance(a, dict):
        if set(a) != set(b):
            return [f"{path}: keys {sorted(set(a) ^ set(b))}"]
        for k in a:
            if k == "elapsed_s":
                continue
            bad += _close(a[k], b[k], f"{path}.{k}")
    elif isinstance(a, list):
        if len(a) != len(b):
            return [f"{path}: len"]
        for i, (x, y) in enumerate(zip(a, b)):
            bad += _close(x, y, f"{path}[{i}]")
    elif isinstance(a, float) and not isinstance(a, bool):
        if path.endswith("max_balance_residual"):  # MUMPS noise floor: both under E_BAND
            if not (a <= E_BAND and b <= E_BAND):
                bad.append(path)
        elif abs(a - b) > L0_BAND * max(abs(a), abs(b)):
            bad.append(f"{path}: {a!r} vs {b!r}")
    elif a != b:
        bad.append(f"{path}: {a!r} vs {b!r}")
    return bad


@complex_only
@needs_s29
def test_wf1_s29_sphere_transient_reproduces_u131():
    """W1: (UT), (Once), (A_t)/(R-A_t), (S2)/(R-S2), (E), (IT), (CLI = API), (c-qs)."""
    comm = MPI.COMM_WORLD
    runs = {}
    for label, path in (("A", A_FILE), ("S", S_FILE)):
        out = LOGS / path.stem
        _fresh_dir(out, comm)
        t0 = time.perf_counter()
        run = run_study_file(path, output_dir=out, comm=comm)
        wall = time.perf_counter() - t0
        with open(out / "summary.json") as fh:
            disk = json.load(fh)
        runs[label] = (run, disk, wall)

    runA, diskA, wallA = runs["A"]
    runS, diskS, wallS = runs["S"]
    thA, thS = runA.thermal, runS.thermal
    n_sub = {k: r.thermal.mesh.topology.index_map(3).size_global for k, (r, _, _) in runs.items()}
    once = {k: r.study.result is r.result for k, (r, _, _) in runs.items()}
    dts = {k: r.thermal.transient.dt_s.hex() for k, (r, _, _) in runs.items()}
    _p0(
        f"\n[WF-1 s29] (UT) sub cells {n_sub} (MAT10_SUB_CELLS {MAT10_SUB_CELLS}); (Once) {once}; "
        f"dt A {dts['A']} (T_STAR/N_MAIN {(T_STAR / N_MAIN).hex()}), S {dts['S']} "
        f"(DT_1 {DT_1.hex()}); wall A {wallA:.1f} s, S {wallS:.1f} s"
    )

    err_A, q_cf = _a_t(thA)
    _p0(
        f"  (A_t) interior rel L2 = {err_A:.10e} [{err_A:.4e}] (band {A_BAND}; U131 "
        f"{U131_A_T_10E}, record {MAT12_A_T_PRINT}); q_cf {q_cf:.6f} W/m^3"
    )
    s2, ratio = _s2(thS, comm)
    _p0(
        f"  (S2) ratio = {ratio:.9f} |ratio - 1| = {s2:.10e} [{s2:.4e}] (band {S2_BAND:.0e}; "
        f"U131 ratio {U131_S2_RATIO}, record {MAT12_S2_PRINT})"
    )
    eA = runA.summary["temperature"]["max_balance_residual"]
    eS = runS.summary["temperature"]["max_balance_residual"]
    _p0(f"  (E) max balance A {eA:.3e}  S {eS:.3e}  (band {E_BAND:.0e})")

    # (c-qs) the steady θ on the solved A Study (no re-solve).
    cool = _thermal(_mapping(A_FILE), str(A_FILE)).cooling
    thQ = runA.study.temperature(cool)
    once_q = runA.study.result is runA.result
    subQ = thQ.submesh
    thetaQ = _exact_theta(subQ, T_STAR, q_cf)
    err_qs = _interior_rel(subQ, _interior_indicator(subQ), thQ.temperature - thetaQ, thetaQ)
    _p0(
        f"  (c-qs) steady theta for theta(t*): {err_qs:.4e} (must miss {A_BAND}; U131 "
        f"1.0108e+00); once {once_q}"
    )

    # (IT-run) §9 item 891: the YAML θ against the same run's Study, U131's spec;
    # no second EM solve (the same `run` object, never a fresh load_study).
    cool_api = RobinCooling(K, H, 0.0, source_scale=G)
    once_run_pre = {"A": runA.study.result is runA.result, "S": runS.study.result is runS.result}
    rA = runA.study.temperature(cool_api, transient=_spec_a())
    rS = runS.study.temperature(cool_api, transient=_spec_s())
    once_run_post = {"A": runA.study.result is runA.result, "S": runS.study.result is runS.result}
    itrA, lrA = _owned_rel(thA.temperature, rA.temperature, comm)
    itrS, lrS = _owned_rel(thS.temperature, rS.temperature, comm)
    _p0(
        f"  (IT-run) A {itrA:.3e}  S {itrS:.3e} (band {L0_BAND:.0e}); layouts {lrA}/{lrS}; "
        f"once before {once_run_pre} after {once_run_post}"
    )
    # (Q-rep) the YAML route's two runs: q_h and integrated source power.
    q_rep, l_rep = _owned_rel(thA.q, thS.q, comm)
    p_rep = abs(float(np.real(thA.source_power)) - float(np.real(thS.source_power)))
    p_rep_rel = p_rep / abs(float(np.real(thA.source_power)))
    _p0(
        f"  (Q-rep) q_h(A run) vs q_h(S run) owned max rel {q_rep:.3e} (layout {l_rep}); "
        f"source power |diff| {p_rep:.3e} W, rel {p_rep_rel:.3e}"
    )

    # (IT-ind) an independently solved Python Study, U131's specs.
    t0 = time.perf_counter()
    study_py = _study(degree=1)
    study_py.solve(comm)
    cool_py = RobinCooling(K, H, 0.0, source_scale=G)
    pyA = study_py.temperature(cool_py, transient=_spec_a())
    pyS = study_py.temperature(cool_py, transient=_spec_s())
    t_py = time.perf_counter() - t0
    itA, layA = _owned_rel(thA.temperature, pyA.temperature, comm)
    itS, layS = _owned_rel(thS.temperature, pyS.temperature, comm)
    qA_ind = _owned_rel(thA.q, pyA.q, comm)[0]
    qS_ind = _owned_rel(thS.q, pyS.q, comm)[0]
    _p0(
        f"  (IT-ind) A {itA:.3e} (band {L0_BAND:.0e}, asserted); layouts {layA}/{layS}; "
        f"Python twin {t_py:.1f} s"
    )
    _p0(
        f"  (IT-ind) S {itS:.3e} printed: pointwise cross-solve spread (known-issues `-n 2` "
        f"MUMPS observation, 2026-11-08 row); a band is the weekly's; q_h owned max rel "
        f"A {qA_ind:.3e}  S {qS_ind:.3e}"
    )
    if q_rep >= 1e-9:
        verdict = "floor of the solve (Q-rep >= 1e-9): recorded"
    elif q_rep < 1e-10 and itS > 1e-9:
        verdict = "route and twin build different q_h (Q-rep < 1e-10, IT-ind S > 1e-9): PARK"
    else:
        verdict = "between 1e-10 and 1e-9: record and stop"
    _p0(f"  (Q-rep) decision: {verdict}")
    # (IT) per-ingredient comparison (the brief's Negative-result print; printed only).
    for k, f_th, p_th in (("A", thA, pyA), ("S", thS, pyS)):
        cells = (
            f_th.mesh.topology.index_map(3).size_global,
            p_th.mesh.topology.index_map(3).size_global,
        )
        nf = tuple(
            comm.allreduce(int(np.sum(t.values == 1)), op=MPI.SUM)
            for t in (f_th.facet_tags, p_th.facet_tags)
        )
        spec = (
            f_th.transient.dt_s.hex(),
            p_th.transient.dt_s.hex(),
            f_th.transient.theta,
            p_th.transient.theta,
            f_th.transient.startup_half_steps,
            p_th.transient.startup_half_steps,
        )
        q_rel = _owned_rel(f_th.q, p_th.q, comm)[0] if hasattr(f_th.q, "function_space") else "n/a"
        sp = (repr(float(np.real(f_th.source_power))), repr(float(np.real(p_th.source_power))))
        _p0(
            f"  (IT-ingr) {k}: cells {cells}; Robin facets (tag 1, owned+ghost sum) {nf}; "
            f"dt/theta/startup {spec}; q_h owned max rel diff {q_rel}; source_power {sp}"
        )

    # (CLI = API by repr)
    json_eq = {}
    for k, (r, disk, _) in runs.items():
        a = {kk: v for kk, v in disk.items() if kk != "elapsed_s"}
        b = {kk: v for kk, v in r.summary.items() if kk != "elapsed_s"}
        json_eq[k] = a == b and all(repr(a[kk]) == repr(b[kk]) for kk in a)
    field_eq = {}
    for k, (r, _, _) in runs.items():
        t, th = r.summary["temperature"], r.thermal
        field_eq[k] = (
            repr(t["max_balance_residual"]) == repr(float(th.max_balance_residual))
            and repr(t["source_power_w"]) == repr(float(np.real(th.source_power)))
            and repr(t["source_multipliers"]) == repr([float(m) for m in th.source_multipliers])
        )
    out_cli = LOGS / f"{A_FILE.stem}_cli"
    _fresh_dir(out_cli, comm)
    t0 = time.perf_counter()
    rc_cli = _rc_all(cli.main(["run", str(A_FILE), "--output-dir", str(out_cli)]), comm)
    t_cli = time.perf_counter() - t0
    with open(out_cli / "summary.json") as fh:
        disk_cli = json.load(fh)
    cli_bad = _close(runA.summary, disk_cli)
    _p0(
        f"  (CLI=API) JSON round trip {json_eq}; summary fields vs result by repr {field_eq}; "
        f"CLI rc {rc_cli}, {t_cli:.1f} s; CLI vs API mismatches {cli_bad}"
    )
    _p0(
        f"  summary['temperature'] A {runA.summary['temperature']['transient']}; theta_max "
        f"{runA.summary['temperature']['theta_max_k']!r}"
    )

    assert all(once.values()) and once_q, "(Once) a second EM solve happened"
    assert all(n == MAT10_SUB_CELLS for n in n_sub.values()), f"(UT) {n_sub}"
    assert dts["A"] == (T_STAR / N_MAIN).hex() and dts["S"] == DT_1.hex(), dts
    assert eA <= E_BAND and eS <= E_BAND, f"(E) A {eA:.3e} S {eS:.3e}"
    assert err_A <= A_BAND, f"(A_t) {err_A:.4e} > {A_BAND}"
    assert f"{err_A:.4e}" == MAT12_A_T_PRINT, f"(R-A_t) {err_A:.4e}"
    assert s2 <= S2_BAND, f"(S2) {s2:.4e} > {S2_BAND:.0e}"
    assert f"{s2:.4e}" == MAT12_S2_PRINT, f"(R-S2) {s2:.4e}"
    assert all(once_run_pre.values()) and all(once_run_post.values()), "(IT-run) once"
    assert lrA and lrS, "(IT-run) owned layouts differ"
    assert itrA <= L0_BAND and itrS <= L0_BAND, f"(IT-run) A {itrA:.3e} S {itrS:.3e}"
    assert layA, "(IT-ind) owned layouts differ"
    assert itA <= L0_BAND, f"(IT-ind) A {itA:.3e}"
    assert all(json_eq.values()), json_eq
    assert all(field_eq.values()), field_eq
    assert rc_cli == (0, 0), rc_cli
    assert cli_bad == [], cli_bad
    assert err_qs > A_BAND, f"(c-qs) {err_qs:.4e} did not miss {A_BAND}"
