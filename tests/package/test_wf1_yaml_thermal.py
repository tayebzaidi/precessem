"""`WF-1` step 3 (§9 item 295, chain U20): a ``thermal:`` section in the YAML.

``femem``'s schema v0 gains one optional top-level section that drives
``Study.temperature(RobinCooling)`` (brief: docs/planning/chunks/WF-1.md,
"§9 item 295").  Asserted:

(K)  the YAML-built ``RobinCooling`` equals ``RobinCooling(K, H, T_INF,
     source_scale=G)`` field by field by ``float.hex``; the Study spec equals
     ``_study(degree=1)``'s; omitted ``degree`` / ``source_scale`` give the
     dataclass defaults;
(X)  ``thermal`` on the F-small file and ``thermal.k_w_per_m_k: 0.0`` raise
     ``StudyConfigError`` naming the field, in-process and through
     ``cli.main`` (rc 2 on every rank, nothing written);
(IT) the ``run_study_file`` θ vs ``Study.temperature`` on an independently
     solved ``_study(degree=1)``, owned max rel ≤ ``L0_BAND``;
(R)  `MAT-10`'s (A) / (P) digits and sub-cell count on the YAML θ; (B), (C)
     bands; ``summary["temperature"]`` equals the in-memory values exactly.
Control (asserted): (c1) the control file's 2G cooling on the solved YAML
Study misses ``A_BAND``.

Every band, constant, record and helper is imported, none restated.
"""

from __future__ import annotations

import copy
import dataclasses
import json
import time
from pathlib import Path

import numpy as np
import pytest
import ufl
from mpi4py import MPI

from fem_em_solver import RobinCooling, StudyConfigError, cli, load_study, run_study_file
from fem_em_solver.post.sar import uniform_sphere_sar_closed_form
from fem_em_solver.study_yaml import _parse_mapping
from fem_em_solver.thermal.steady import solve_steady_conduction
from fem_em_solver.utils.constants import EPSILON_0

from tests.complex_mode import complex_only
from tests.package.test_wf15_study_lossy_sphere import _study
from tests.package.test_wf15_study_temperature import (
    MAT10_A_PRINT,
    MAT10_P_PRINT,
    MAT10_SUB_CELLS,
    _theta_owned_rel,
)
from tests.package.test_wf1_study_yaml import (
    LOGS,
    WF1,
    _fresh_dir,
    _mapping,
    _rc_all,
    _spec_diff,
)
from tests.validation.test_lossy_sphere_sar import (
    E0,
    EPSILON_R_SPHERE,
    OMEGA,
    RHO_KG_M3,
    SIGMA_HIGH,
)
from tests.validation.test_mat10_sar_to_heat_sphere import (
    A_BAND,
    B_BAND,
    C_BAND,
    G,
    H,
    K,
    P_BAND,
    ROBIN_TAG,
    T_INF,
    _allsum,
    _interior_indicator,
    _interior_rel,
    _theta_cf,
)
from tests.validation.test_mat5_em_coupled_sphere import L0_BAND

THERMAL_FILE = WF1 / "sphere_mat10_thermal.yaml"
CONTROL_2G_FILE = WF1 / "sphere_mat10_thermal_control_2g.yaml"
MALFORMED_THERMAL = {
    "malformed_thermal_on_port.yaml": "thermal",
    "malformed_thermal_k_zero.yaml": "thermal.k_w_per_m_k",
}


def _p0(*args):
    if MPI.COMM_WORLD.rank == 0:
        print(*args, flush=True)


def _hex_diff(a: RobinCooling, b: RobinCooling) -> list:
    out = []
    for f in dataclasses.fields(RobinCooling):
        x, y = getattr(a, f.name), getattr(b, f.name)
        if type(x) is not type(y):
            out.append(f.name)
        elif isinstance(x, float):
            if x.hex() != y.hex():
                out.append(f.name)
        elif x != y:
            out.append(f.name)
    return out


def _cooling(path):
    return _parse_mapping(_mapping(path), str(path))[3]


def test_wf1_thermal_yaml_parse_and_refusals():
    """(K) and (X): no mesh, no solve."""
    comm = MPI.COMM_WORLD
    c_yaml = _cooling(THERMAL_FILE)
    ref = RobinCooling(K, H, T_INF, source_scale=G)
    d_cool = _hex_diff(c_yaml, ref)
    d_spec = _spec_diff(load_study(THERMAL_FILE), _study(degree=1))
    c2 = _cooling(CONTROL_2G_FILE)
    d_c2 = _hex_diff(c2, ref)
    _p0(f"\n[WF-1 s3] (K) cooling {c_yaml}; diff vs RobinCooling(K, H, T_INF, G) {d_cool}; "
        f"spec diff vs _study(degree=1) {d_spec}; control diff {d_c2} "
        f"(source_scale {c2.source_scale!r} vs 2G {2.0 * G!r})")
    assert d_cool == [], d_cool
    assert d_spec == [], d_spec
    assert d_c2 == ["source_scale"] and c2.source_scale == 2.0 * G

    # Omitted optional keys: the dataclass defaults, never a loader default.
    m = copy.deepcopy(_mapping(THERMAL_FILE))
    del m["thermal"]["source_scale"]
    c_def = _parse_mapping(m, "<omitted>")[3]
    defaults = RobinCooling(K, H, T_INF)
    _p0(f"  (K) omitted source_scale / degree: {c_def.source_scale!r} / {c_def.degree!r} "
        f"(defaults {defaults.source_scale!r} / {defaults.degree!r})")
    assert _hex_diff(c_def, defaults) == []
    assert c_def.degree == 1 and c_def.source_scale == 1.0
    # No thermal section: no cooling.
    assert _parse_mapping(_mapping(WF1 / "sphere_mat4.yaml"), "<m>")[3] is None

    out = LOGS / "malformed_thermal"
    fields = []
    for name, expected in MALFORMED_THERMAL.items():
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
    assert len(set(fields)) == 2, fields


@complex_only
def test_wf1_thermal_yaml_reproduces_mat10_record():
    """(IT), (R), (c1) through ``run_study_file`` at degree 1, -n 2."""
    comm = MPI.COMM_WORLD
    out = LOGS / THERMAL_FILE.stem
    _fresh_dir(out, comm)
    t0 = time.perf_counter()
    run = run_study_file(THERMAL_FILE, output_dir=out, comm=comm)
    t_yaml = time.perf_counter() - t0
    study, res, th = run.study, run.result, run.thermal
    assert th is not None
    assert study.result is res, "(Once) a second EM solve happened"
    with open(out / "summary.json") as fh:
        summary = json.load(fh)
    _p0(f"\n[WF-1 s3] run_study_file {t_yaml:.1f} s (summary elapsed "
        f"{summary['elapsed_s']:.1f} s); summary['temperature'] {summary['temperature']}")

    # (R) JSON round trip equals the in-memory values exactly.
    T = th.temperature
    imap = T.function_space.dofmap.index_map
    n_own = imap.size_local * T.function_space.dofmap.index_map_bs
    loc = np.real(np.asarray(T.x.array[:n_own])) - T_INF
    theta_max = comm.allreduce(float(np.max(loc)) if loc.size else -np.inf, op=MPI.MAX)
    mem = {
        "source_power_w": float(np.real(th.source_power)),
        "convected_power_w": float(np.real(th.convected_power)),
        "theta_max_k": float(theta_max),
    }
    same_json = summary["temperature"] == mem and run.summary["temperature"] == mem
    _p0(f"  (R) in-memory {mem}; JSON equal {same_json}")

    # (R) MAT-10's gates on the YAML θ.
    cf = uniform_sphere_sar_closed_form(
        e0=E0, epsilon_r=EPSILON_R_SPHERE, sigma=SIGMA_HIGH, omega=OMEGA,
        rho=RHO_KG_M3, epsilon_0=EPSILON_0,
    )
    q_cf = G * RHO_KG_M3 * cf["sar_w_per_kg"]
    sub, ft = th.submesh, th.facet_tags
    n_sub = sub.topology.index_map(3).size_global
    chi = _interior_indicator(sub)
    theta = _theta_cf(sub, q_cf)
    err_A = _interior_rel(sub, chi, (T - T_INF) - theta, theta)
    res_P = solve_steady_conduction(sub, ft, k=K, q=q_cf, h=H, t_inf=T_INF,
                                    robin_tag=ROBIN_TAG, degree=1)
    err_P = _interior_rel(sub, chi, (res_P.temperature - T_INF) - theta, theta)
    int_q = float(np.real(_allsum(comm, th.heat_source * ufl.dx(domain=sub))))
    p_parent = th.cooling.source_scale * study.mean_sar()["dissipated_power_w"]
    err_B = abs(int_q - p_parent) / abs(p_parent)
    err_C = abs(th.convected_power - th.source_power) / abs(th.source_power)
    _p0(f"  (R) sub cells {n_sub} (record {MAT10_SUB_CELLS}); (A) = {err_A:.10e} [{err_A:.4e}] "
        f"vs {MAT10_A_PRINT} (band {A_BAND}); (P) = {err_P:.10e} [{err_P:.4e}] vs "
        f"{MAT10_P_PRINT} (band {P_BAND:.0e}); (B) = {err_B:.3e} (band {B_BAND:.0e}); "
        f"(C) = {err_C:.3e} (band {C_BAND:.0e})")

    # (c1) missing ½: the control file's cooling on the already-solved Study.
    c2 = _cooling(CONTROL_2G_FILE)
    th2 = study.temperature(c2)
    once_c1 = study.result is res
    sub2 = th2.submesh
    theta2 = _theta_cf(sub2, q_cf)
    c1 = _interior_rel(sub2, _interior_indicator(sub2), (th2.temperature - T_INF) - theta2, theta2)
    _p0(f"  (c1) control file source_scale {c2.source_scale!r}: (A) comparator {c1:.4e} "
        f"(band {A_BAND}); no re-solve {once_c1}")

    # (IT) an independently solved Python Study.
    t0 = time.perf_counter()
    study_py = _study(degree=1)
    study_py.solve(comm)
    th_py = study_py.temperature(RobinCooling(K, H, T_INF, source_scale=G))
    t_py = time.perf_counter() - t0
    it, th_max = _theta_owned_rel(T, th_py.temperature, comm)
    _p0(f"  (IT) owned max|dT| / max|T - Tinf| = {it:.3e} (band {L0_BAND:.0e}); "
        f"theta_py max {th_max:.6f} K; Python Study + thermal {t_py:.1f} s")

    assert summary["family"] == "sphere_in_box" and summary["degree"] == 1
    assert same_json, (summary["temperature"], mem)
    assert n_sub == MAT10_SUB_CELLS, f"(R) sub cells {n_sub}"
    assert f"{err_A:.4e}" == MAT10_A_PRINT, f"(R) (A) {err_A:.4e}"
    assert f"{err_P:.4e}" == MAT10_P_PRINT, f"(R) (P) {err_P:.4e}"
    assert err_A <= A_BAND and err_P <= P_BAND
    assert err_B <= B_BAND, f"(B) {err_B:.3e}"
    assert err_C <= C_BAND, f"(C) {err_C:.3e}"
    assert once_c1, "(c1) re-solved"
    assert c1 > A_BAND, f"(c1) {c1:.4e} did not miss {A_BAND}"
    assert it <= L0_BAND, f"(IT) {it:.3e} > {L0_BAND:.0e}"
