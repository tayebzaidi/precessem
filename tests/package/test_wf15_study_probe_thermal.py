"""`WF-15` step 6 (§9 item 523, §10 U60): ``Study.probe`` on a thermal-only
Study, against `MAT-19`'s and `MAT-21`'s comparands.

``Study.probe(..., transient=...)`` passes the ``TransientSpec`` through to
``Study.temperature`` on a ``ThermalBallGeometry``; on a ``ThermalBoxGeometry``
the steady branch is reached unchanged.  Anchors (brief:
docs/planning/chunks/WF-15.md, "§9 item 523"):

(H)       ``omega`` / ``epsilon_c()`` / ``dirichlet_field()`` on a thermal-only
          Study raise ``StudyConfigError`` naming ``frequency_hz`` /
          ``material`` / ``excitation``, no mesh built; unchanged on an EM Study.
(A-ball)  the whole-ball probe reads ∫θ dx / V_h and V_h within ``B_BAND``.
(T*-ball) the three ``CENTRES`` probes at t* read the exact Helmholtz ball
          mean ``_ball_mean`` within ``P_BAND`` (the gate, not the digits).
(M)       every box probe ball lies inside the 0.9-interior box.
(R-box)   five box probes read θ_series(c) − Q a²/(10K) within ``INTERIOR_BAND``.
(A-box)   the whole-box probe identity within ``B_BAND``.
Controls (rule (e)): (c-disp) and (c-wall) *predicted*, only the miss asserted;
(c-mv) printed only; (c-guard) asserted, backed by W0a–c; (c-pass) by
construction.  Opt-in solves: ``FEM_EM_WF15_S6=1`` (real build, ``-n 2``).
"""

from __future__ import annotations

import os
import time

import numpy as np
import pytest
import ufl
from dolfinx import default_scalar_type, fem
from mpi4py import MPI

import fem_em_solver.study as study_mod
from fem_em_solver import (
    RobinCooling,
    Study,
    ThermalBallGeometry,
    ThermalBoxGeometry,
    TransientSpec,
    UniformHeatSource,
)
from fem_em_solver.study_yaml import StudyConfigError
from fem_em_solver.thermal.probe import ProbeBall
from fem_em_solver.utils.constants import EPSILON_0

from tests.complex_mode import IS_COMPLEX
from tests.package.test_wf15_study_lossy_sphere import FREQUENCY_HZ, _study
from tests.validation.test_mat10_sar_to_heat_sphere import (
    B_BAND,
    H,
    K,
    P_BAND,
    R,
    T_INF,
    _allsum,
    _say,
)
from tests.validation.test_mat12_transient_sar_heating_sphere import N_MAIN, RHO_C, T_STAR
from tests.validation.test_mat13_transient_perfusion_sphere import INTERIOR_BAND
from tests.validation.test_mat15_coil_driven_heating_phantom import H_GEL, K_GEL
from tests.validation.test_mat15_coil_driven_heating_phantom import T_INF as T_INF_GEL
from tests.validation.test_mat19_virtual_temperature_probe import (
    A_PROBE,
    CENTRES,
    DISPLACED,
    _ball_mean,
    _q_cf,
)
from tests.validation.test_mat21_box_closed_form import A_BOX, B_BOX, C_BOX, HALF, Q
from tests.validation.test_mat21_box_fem import _steady_ex

S6_ENV = "FEM_EM_WF15_S6"

# Records printed beside, never compared (transcribed with their log lines).
MAT19_T_MAX = 1.7315e-03  # 20261005T095808Z_MAT-19-step1-W1.log.gz:203 (max over 3 x 40)
MAT19_C1 = 5.421e-02  # ibid. :209, the displaced probe at n = 40
STEP1_M4_WHOLE = 3.459853e-03  # 20261013T142646Z_MAT-21-step1r-W1-steady.log:102
S5_C_WALL = 3.8048e-01  # item 521a's (c-wall) interior miss (commit f437514a)

WHOLE_BALL = ProbeBall((0.0, 0.0, 0.0), 2.0 * R)  # every ball vertex lies at r <= R
BALL_PROBES = [ProbeBall(c, A_PROBE) for c in CENTRES]
DISP_PROBE = ProbeBall(DISPLACED, A_PROBE)

A_BOXPROBE = 0.0135  # 2 h_eff at m = 4 (h_eff 6.750514e-03, ...step1r-W1-steady.log:102)
BOX_CENTRES = [
    (0.0, 0.0, 0.0),
    (0.0, 0.5 * B_BOX, 0.0),
    (0.0, -0.5 * B_BOX, 0.0),
    (0.5 * A_BOX, 0.0, 0.0),
    (-0.5 * A_BOX, 0.0, 0.0),
]
BOX_PROBES = [ProbeBall(c, A_BOXPROBE) for c in BOX_CENTRES]
WHOLE_BOX = ProbeBall((0.0, 0.0, 0.0), 0.4)  # beyond the half-diagonal, 0.3896 m
BOX_DIVISIONS = (52, 32, 8)  # m = 4 of MAT-21's BASE (13, 8, 2)

MSG = "is not defined on a thermal-only Study"
# (c-guard): what W0a / W0b / W0c logged before the edit
# (20261018T123111Z / 123115Z / 123117Z_WF-15-s6-W0*.log).
UNGUARDED = {
    "omega": (TypeError, "unsupported operand type"),
    "epsilon_c": (AttributeError, "'NoneType' object has no attribute 'epsilon_r'"),
    "dirichlet_field": (AttributeError, "object has no attribute 'dirichlet_field'"),
}
FIELDS = {"omega": "frequency_hz", "epsilon_c": "material", "dirichlet_field": "excitation"}


def _knob() -> bool:
    raw = os.environ.get(S6_ENV, "").strip()
    if raw == "":
        return False
    if raw == "1":
        return True
    raise ValueError(f"{S6_ENV}={raw!r}: only unset/empty or '1' are accepted")


def _require_solve_env():
    if not _knob():
        pytest.skip(f"{S6_ENV}=1 opts in to the thermal probe solve legs")
    assert not IS_COMPLEX, "the WF-15 step-6 solve legs run in the real build"
    assert MPI.COMM_WORLD.size == 2, "the WF-15 step-6 solve legs run at -n 2"


def _call(st, name):
    return st.omega if name == "omega" else getattr(st, name)()


def _domain_mean(theta, comm):
    msh = theta.function_space.mesh
    one = fem.Constant(msh, default_scalar_type(1.0))
    integral = float(np.real(_allsum(comm, theta * ufl.dx(domain=msh))))
    v_h = float(np.real(_allsum(comm, one * ufl.dx(domain=msh))))
    return integral / v_h, v_h


def _series_at(c, **kw):
    f = _steady_ex(**kw)
    return float(np.asarray(f(np.array([c[0]]), np.array([c[1]]), np.array([c[2]]))).ravel()[0])


def test_wf15_s6_helpers_refuse(monkeypatch):
    calls = []

    def _no_build(self, comm=MPI.COMM_WORLD):
        calls.append(1)
        raise AssertionError("build must not be called")

    monkeypatch.setattr(ThermalBallGeometry, "build", _no_build)
    monkeypatch.setattr(ThermalBoxGeometry, "build", _no_build)
    studies = {
        "ball": Study(ThermalBallGeometry(0.1, 0.005), None, UniformHeatSource(2000.0)),
        "box": Study(ThermalBoxGeometry(HALF, BOX_DIVISIONS), None, UniformHeatSource(Q)),
    }
    # (H)
    for label, st in studies.items():
        for name, field in FIELDS.items():
            with pytest.raises(StudyConfigError) as exc:
                _call(st, name)
            assert exc.value.field == field, (label, name, exc.value.field)
            assert isinstance(exc.value, ValueError)
            assert MSG in str(exc.value) and name in str(exc.value)
            _say(f"WF15-6 (H) {label}.{name}: {type(exc.value).__name__} field={exc.value.field!r}")
    # (c-guard) asserted, backed by W0a-c: the helper removed restores the raw exceptions.
    with monkeypatch.context() as m:
        m.setattr(study_mod, "_refuse_thermal_only_attr", lambda *a, **k: None)
        for label, st in studies.items():
            for name, (etype, frag) in UNGUARDED.items():
                with pytest.raises(etype) as exc:
                    _call(st, name)
                assert not isinstance(exc.value, StudyConfigError)
                assert frag in str(exc.value), (label, name, str(exc.value))
                _say(f"WF15-6 control (c-guard) {label}.{name}: {type(exc.value).__name__}")
    # (c-pass) by construction: no transient on the ball still raises TypeError.
    with pytest.raises(TypeError, match="needs transient=TransientSpec"):
        studies["ball"].probe(BALL_PROBES, RobinCooling(K, H, T_INF))
    assert calls == []
    # Unchanged on an EM Study (no solve), hex within one process.
    em = _study()
    assert em.omega.hex() == (2.0 * np.pi * FREQUENCY_HZ).hex()
    eps = em.epsilon_c()
    ref = complex(em.material.epsilon_r, -em.material.sigma / (em.omega * EPSILON_0))
    assert eps.real.hex() == ref.real.hex() and eps.imag.hex() == ref.imag.hex()
    _say(f"WF15-6 (H) EM omega {em.omega!r}, epsilon_c {eps!r}: hex-equal to the formulas")


def test_wf15_s6_ball_probe_reproduces_mat19_gate():
    _require_solve_env()
    comm = MPI.COMM_WORLD
    q = _q_cf()
    st = Study(ThermalBallGeometry(R, 1.0e-3), None, UniformHeatSource(q))
    t0 = time.perf_counter()
    res = st.probe([WHOLE_BALL, *BALL_PROBES, DISP_PROBE], RobinCooling(K, H, T_INF),
                   transient=TransientSpec(T_STAR / N_MAIN, N_MAIN, RHO_C))
    _say(f"WF15-6 ball Study.probe (CN N={N_MAIN}, t*={T_STAR:.4f} s): "
         f"{time.perf_counter() - t0:.1f} s")
    msh = res.theta.function_space.mesh
    cells = msh.topology.index_map(3).size_global
    mean_h, v_h = _domain_mean(res.theta, comm)
    p0, vol0 = float(res.values[0]), float(res.volumes[0])
    a_val = abs(p0 - mean_h) / abs(mean_h)
    a_vol = abs(vol0 / v_h - 1.0)
    _say(f"WF15-6 ball cells {cells}; V_h {v_h:.6e} (4/3 pi R^3 {4.0 * np.pi * R**3 / 3.0:.6e}); "
         f"max_imag_ratio {res.max_imag_ratio:.1e}")
    _say(f"WF15-6 (A-ball) |P - mean| / |mean| = {a_val:.3e}, |V_probe / V_h - 1| = {a_vol:.3e} "
         f"(bound {B_BAND:.0e})")
    ratios = []
    for j, c in enumerate(CENTRES):
        d = float(np.linalg.norm(c))
        cf = _ball_mean(d, A_PROBE, T_STAR, q)
        p = float(res.values[1 + j])
        r = abs(p / cf - 1.0)
        ratios.append(r)
        _say(f"WF15-6 (T*-ball) c{j} |c| = {d:.4e}: P = {p:.6e}, ball mean = {cf:.6e}, "
             f"rel = {r:.4e} (bound {P_BAND:.0e}; predicted <= 3e-3); volume "
             f"{float(res.volumes[1 + j]):.6e} / nominal {float(res.nominal_volumes[1 + j]):.6e}")
    _say(f"WF15-6 (T*-ball) max rel = {max(ratios):.4e} (MAT-19's own-mesh max over 40 times "
         f"{MAT19_T_MAX:.4e}, printed only)")
    d2 = float(np.linalg.norm(CENTRES[2]))
    cf2 = _ball_mean(d2, A_PROBE, T_STAR, q)
    p_disp = float(res.values[-1])
    c_disp = abs(p_disp / cf2 - 1.0)
    sep = abs(_ball_mean(float(np.linalg.norm(DISPLACED)), A_PROBE, T_STAR, q) / cf2 - 1.0)
    _say(f"WF15-6 control (c-disp) displaced vs c2's comparand: {c_disp:.4e} (predicted 5.4e-2, "
         f"MAT-19 {MAT19_C1:.4e}; comparand separation {sep:.4e}; asserted > {P_BAND:.0e})")
    assert a_val <= B_BAND
    assert a_vol <= B_BAND
    for r in ratios:
        assert r <= P_BAND
    assert c_disp > P_BAND


def test_wf15_s6_box_probe_robin_series():
    _require_solve_env()
    comm = MPI.COMM_WORLD
    # (M) every probe ball inside the 0.9-interior box.
    for c in BOX_CENTRES:
        for x, hw in zip(c, HALF):
            assert abs(x) + A_BOXPROBE <= 0.9 * hw, (c, hw)
    geo = ThermalBoxGeometry(HALF, BOX_DIVISIONS)
    t0 = time.perf_counter()
    res = Study(geo, None, UniformHeatSource(Q)).probe(
        [WHOLE_BOX, *BOX_PROBES], RobinCooling(K_GEL, H_GEL, T_INF_GEL))
    _say(f"WF15-6 box Study.probe steady m=4: {time.perf_counter() - t0:.1f} s")
    msh = res.theta.function_space.mesh
    cells = msh.topology.index_map(3).size_global
    mean_h, v_h = _domain_mean(res.theta, comm)
    p0, vol0 = float(res.values[0]), float(res.volumes[0])
    a_val = abs(p0 - mean_h) / abs(mean_h)
    a_vol = abs(vol0 / v_h - 1.0)
    _say(f"WF15-6 box cells {cells}; V_h {v_h:.6e} (8ABC {8.0 * A_BOX * B_BOX * C_BOX:.6e}); "
         f"max_imag_ratio {res.max_imag_ratio:.1e}")
    _say(f"WF15-6 (A-box) |P - mean| / |mean| = {a_val:.3e}, |V_probe / V_h - 1| = {a_vol:.3e} "
         f"(bound {B_BAND:.0e})")
    mv = Q * A_BOXPROBE**2 / (10.0 * K_GEL)
    ratios, vals = [], []
    for j, c in enumerate(BOX_CENTRES):
        point = _series_at(c)
        cf = point - mv
        p = float(res.values[1 + j])
        r = abs(p / cf - 1.0)
        ratios.append(r)
        vals.append(p)
        _say(f"WF15-6 (R-box) c{j} = {c}: P = {p:.6e}, theta_series = {point:.6e}, "
             f"ball mean = {cf:.6e}, rel = {r:.4e} (bound {INTERIOR_BAND:.0e}; predicted <= 5e-3; "
             f"m=4 whole-box L2 {STEP1_M4_WHOLE:.6e}); volume {float(res.volumes[1 + j]):.6e} / "
             f"nominal {float(res.nominal_volumes[1 + j]):.6e}")
    _say(f"WF15-6 (R-box) max rel = {max(ratios):.4e}")
    _say(f"WF15-6 mirror pairs (printed only): y {vals[1] - vals[2]:.3e} K, "
         f"x {vals[3] - vals[4]:.3e} K")
    point0 = _series_at(BOX_CENTRES[0])
    _say(f"WF15-6 control (c-mv) point theta_series(c0) {point0:.6e} vs ball mean "
         f"{point0 - mv:.6e}: Q a^2/(10K) = {mv:.4e} K ({mv / point0:.3e} relative; printed only)")
    res2 = Study(geo, None, UniformHeatSource(Q)).probe(
        BOX_PROBES[:1], RobinCooling(K_GEL, 2.0 * H_GEL, T_INF_GEL))
    c_wall = abs(float(res2.values[0]) / (point0 - mv) - 1.0)
    _say(f"WF15-6 control (c-wall) h doubled at c0: {c_wall:.4e} (predicted 0.35, range 0.30-0.42; "
         f"item 521a interior {S5_C_WALL:.4e}; asserted > {INTERIOR_BAND:.0e})")
    assert a_val <= B_BAND
    assert a_vol <= B_BAND
    for r in ratios:
        assert r <= INTERIOR_BAND
    assert c_wall > INTERIOR_BAND
