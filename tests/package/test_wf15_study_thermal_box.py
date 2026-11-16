"""`WF-15` step 5 (§9 item 521a, §10 U59 part 1): `MAT-21`'s box through the API.

``ThermalBoxGeometry`` (``create_box``, Robin cooling on all six faces) drives
``Study.temperature`` steady and Crank–Nicolson with a duty cycle.  The solve
legs are gated in the complex build, the build of the step-1 / P26 records.

Asserted: (K0) spec and refusals, no mesh; (UT) the Study's mesh / facet tags
equal ``_rung``'s; (I-dir) θ_Study = θ_direct ≤ ``L0_BAND``; (S1-m4) step 1's
m = 4 whole-box record at ``.4e``, interior ≤ ``INTERIOR_BAND``, Gauss ≤
``S2_GAUSS_BAND``; (S0) duty multipliers equal the direct run's; (T-duty) P26's
record at ``.4e``, interior ≤ ``INTERIOR_BAND``, balance ≤ ``E_BAND``.
Controls (rule (e)): (c3) the unscheduled comparand misses (asserted, backed by
P26's (c3)); (c-wall) Robin h doubled misses (*predicted*, only the miss is
asserted).  Opt-in solves: ``FEM_EM_WF15_S5=1``.
"""

from __future__ import annotations

import os

import numpy as np
import pytest
import ufl
from dolfinx import default_scalar_type, fem
from mpi4py import MPI

import fem_em_solver.study as study_mod
from fem_em_solver import (
    DutyCycle,
    RobinCooling,
    Study,
    ThermalBallGeometry,
    ThermalBoxGeometry,
    TransientSpec,
    UniformHeatSource,
)

from tests.complex_mode import IS_COMPLEX
from tests.package.test_wf13_cli_estimate import THERMAL_MESSAGE
from tests.package.test_wf15_study_transient import _on, _same_mesh, _theta_owned_rel
from tests.validation.test_mat5_em_coupled_sphere import L0_BAND
from tests.validation.test_mat10_sar_to_heat_sphere import ROBIN_TAG, _allsum, _say
from tests.validation.test_mat12_transient_sar_heating_sphere import E_BAND, RHO_C
from tests.validation.test_mat13_transient_perfusion_sphere import INTERIOR_BAND
from tests.validation.test_mat15_coil_driven_heating_phantom import H_GEL, K_GEL, N_MAIN, T_INF
from tests.validation.test_mat15_finite_cylinder_fem import S2_GAUSS_BAND, _max_imag, _rel
from tests.validation.test_mat21_box_closed_form import HALF, Q
from tests.validation.test_mat21_box_duty_closed_form import DT_BOX, HALF_T, T_STAR, _duty_cf
from tests.validation.test_mat21_box_duty_fem import M_RUNG, _chi_box, _cn
from tests.validation.test_mat21_box_fem import (
    BASE,
    _Clock,
    _exact,
    _interior_chi,
    _rung,
    _steady_ex,
    _trans_ex,
)

S5_ENV = "FEM_EM_WF15_S5"

# Records (transcribed with their log lines).
STEP1_M4_WHOLE = "3.4599e-03"  # 20261013T142646Z_MAT-21-step1r-W1-steady.log:112 (3.459853e-03 at :102)
STEP1_M4_CELLS = 79872  # ibid. :102
STEP1_M4_P1_DOFS = 15741  # ibid. :102
STEP1_M6_CELLS = 269568  # ibid. :106
P26_T_DUTY = "9.4014e-04"  # 20261014T003506Z_MAT-21-step2-a2-W3-closed-form.log:92 (9.401399e-04)
P26_C3 = 5.0282e-01  # ibid. :93


def _knob() -> bool:
    raw = os.environ.get(S5_ENV, "").strip()
    if raw == "":
        return False
    if raw == "1":
        return True
    raise ValueError(f"{S5_ENV}={raw!r}: only unset/empty or '1' are accepted")


def _require_solve_env():
    if not _knob():
        pytest.skip(f"{S5_ENV}=1 opts in to the box solve legs")
    assert IS_COMPLEX, "the WF-15 step-5 solve legs are gated in the complex build"
    assert MPI.COMM_WORLD.size == 2, "the WF-15 step-5 solve legs run at -n 2"


def _all(comm, flag):
    return bool(comm.allreduce(int(bool(flag)), op=MPI.MIN))


def _same_tags(ft_a, ft_b, comm):
    ok = np.array_equal(ft_a.indices, ft_b.indices) and np.array_equal(ft_a.values, ft_b.values)
    return _all(comm, ok)


def _box(m):
    return ThermalBoxGeometry(HALF, tuple(b * m for b in BASE))


def test_wf15_box_spec_and_refusals(monkeypatch):
    calls = []

    def _no_build(self, comm=MPI.COMM_WORLD):
        calls.append(1)
        raise AssertionError("build must not be called")

    monkeypatch.setattr(ThermalBoxGeometry, "build", _no_build)
    geo = ThermalBoxGeometry(HALF, (52, 32, 8))
    assert geo.wall is None
    st = Study(geo, None, UniformHeatSource(Q))
    assert st.material is None and st.frequency_hz is None
    with pytest.raises(TypeError):
        Study(geo, "gel", UniformHeatSource(Q))
    with pytest.raises(TypeError):
        Study(geo, None, object())
    with pytest.raises(TypeError):
        st.solve()
    with pytest.raises(TypeError):
        st.mean_sar()
    with pytest.raises(TypeError, match="ThermalBoxGeometry has no EM solve"):
        st.estimate()
    cool = RobinCooling(K_GEL, H_GEL, T_INF)
    with pytest.raises(ValueError, match="duty needs a transient"):
        st.temperature(cool, duty=DutyCycle(period_s=2.0 * HALF_T, duty=0.5))
    with pytest.raises(ValueError, match="source_scale"):
        st.temperature(RobinCooling(K_GEL, H_GEL, T_INF, source_scale=2.0))
    for hw, dv in [((0.1, 0.1), (1, 1, 1)), ((0.1, 0.0, 0.1), (1, 1, 1)),
                   ((0.1, 0.1, 0.1), (1, 0, 1))]:
        with pytest.raises(ValueError):
            ThermalBoxGeometry(hw, dv)
    assert calls == []
    # The ball's estimate message is unchanged.
    ball = Study(ThermalBallGeometry(0.05, 0.01), None, UniformHeatSource(Q))
    with pytest.raises(TypeError) as exc:
        ball.estimate()
    assert str(exc.value) == THERMAL_MESSAGE
    # (c0) by construction: the box branch is taken (no frequency_hz needed).
    assert isinstance(study_mod.Study(geo, None, UniformHeatSource(Q)).excitation, UniformHeatSource)


def test_wf15_box_steady_reproduces_mat21_step1_m4():
    _require_solve_env()
    comm = MPI.COMM_WORLD
    clk = _Clock("wf15-s5-steady")
    msh, ft, _ns, n, v_h, st = _rung(4, clk)
    cool = RobinCooling(K_GEL, H_GEL, T_INF)
    res = Study(_box(4), None, UniformHeatSource(Q)).temperature(cool)
    clk.lap("Study steady m=4")
    smsh = res.submesh
    ca, cb, geo_eq = _same_mesh(smsh, msh, comm)
    tags_eq = _same_tags(res.facet_tags, ft, comm)
    _say(f"WF15-5 (UT) cells Study {ca} / direct {cb} (record {STEP1_M4_CELLS}); "
         f"geometry equal {geo_eq}; facet tags equal {tags_eq}")
    assert ca == cb == STEP1_M4_CELLS and geo_eq and tags_eq
    T = _on(msh, res.temperature)
    d = _theta_owned_rel(T, st.temperature, comm)
    n_own = st.temperature.function_space.dofmap.index_map.size_local
    hex_eq = _all(comm, np.array_equal(T.x.array[:n_own], st.temperature.x.array[:n_own]))
    _say(f"WF15-5 (I-dir) steady theta_Study vs direct: {d:.3e} (bound {L0_BAND:.0e}); "
         f"owned arrays bitwise {hex_eq} (printed only)")
    assert d <= L0_BAND
    th_h = T - T_INF
    ex, npts = _exact(msh, _steady_ex())
    clk.lap(f"steady series interpolation ({npts} P2 points)")
    whole = _rel(msh, th_h, ex)
    chi = _interior_chi(msh)
    interior = _rel(msh, th_h, ex, chi)
    p1 = T.function_space.dofmap.index_map.size_global
    ds = ufl.Measure("ds", domain=msh, subdomain_data=ft)(ROBIN_TAG)
    conv = float(np.real(_allsum(comm, H_GEL * th_h * ds)))
    one = fem.Constant(msh, default_scalar_type(1.0))
    vh = float(np.real(_allsum(comm, one * ufl.dx(domain=msh))))
    gauss = abs(conv - Q * vh) / (Q * vh)
    _say(f"WF15-5 (S1-m4) whole L2 rel = {whole:.6e} -> {whole:.4e} (record {STEP1_M4_WHOLE}); "
         f"interior = {interior:.6e} (bound {INTERIOR_BAND:.0e}; new reading at the step's band); "
         f"P1 dofs {p1} (record {STEP1_M4_P1_DOFS}); Gauss {gauss:.3e} (bound {S2_GAUSS_BAND:.0e}); "
         f"max imag {_max_imag(res.temperature):.1e}")
    # (c-wall) predicted: Robin h doubled through the API.
    res2 = Study(_box(4), None, UniformHeatSource(Q)).temperature(RobinCooling(K_GEL, 2.0 * H_GEL, T_INF))
    T2 = _on(msh, res2.temperature)
    c_wall = _rel(msh, T2 - T_INF, ex, chi)
    _say(f"WF15-5 control (c-wall) h doubled: interior miss = {c_wall:.4e} "
         f"(predicted 0.37, range 0.32-0.42, printed only; asserted > {INTERIOR_BAND:.0e})")
    assert f"{whole:.4e}" == STEP1_M4_WHOLE
    assert interior <= INTERIOR_BAND
    assert p1 == STEP1_M4_P1_DOFS
    assert gauss <= S2_GAUSS_BAND
    assert c_wall > INTERIOR_BAND


def test_wf15_box_duty_reproduces_mat21_p26():
    _require_solve_env()
    comm = MPI.COMM_WORLD
    clk = _Clock("wf15-s5-duty")
    msh, ft, _ns, n, v_h, _st = _rung(M_RUNG, clk)
    direct = _cn(msh, ft, _chi_box)
    direct.snapshots.clear()
    clk.lap(f"direct CN N={N_MAIN}")
    period = 2.0 * HALF_T
    res = Study(_box(M_RUNG), None, UniformHeatSource(Q)).temperature(
        RobinCooling(K_GEL, H_GEL, T_INF),
        transient=TransientSpec(DT_BOX, N_MAIN, RHO_C),
        duty=DutyCycle(period_s=period, duty=0.5),
    )
    clk.lap("Study CN duty")
    ca, cb, geo_eq = _same_mesh(res.mesh, msh, comm)
    tags_eq = _same_tags(res.facet_tags, ft, comm)
    _say(f"WF15-5 (UT) cells Study {ca} / direct {cb} (record {STEP1_M6_CELLS}); "
         f"geometry equal {geo_eq}; facet tags equal {tags_eq}")
    assert ca == cb == STEP1_M6_CELLS and geo_eq and tags_eq
    mult_eq = np.array_equal(np.asarray(res.source_multipliers),
                             np.asarray(direct.source_multipliers))
    _say(f"WF15-5 (S0) multipliers equal {mult_eq} ({len(res.source_multipliers)} entries)")
    assert mult_eq
    T = _on(msh, res.temperature)
    d = _theta_owned_rel(T, direct.temperature, comm)
    _say(f"WF15-5 (I-dir) duty theta_Study vs direct: {d:.3e} (bound {L0_BAND:.0e})")
    assert d <= L0_BAND
    th = T - T_INF
    chi = _interior_chi(msh)
    ex, npts = _exact(msh, lambda x, y, z: _duty_cf(x, y, z, T_STAR))
    clk.lap(f"ten-run series interpolation ({npts} P2 points)")
    t_duty = _rel(msh, th, ex, chi)
    ex_s, _ = _exact(msh, _trans_ex(T_STAR))
    clk.lap("step-response interpolation (c3)")
    c3 = _rel(msh, th, ex_s, chi)
    bal = float(res.max_balance_residual)
    _say(f"WF15-5 (T-duty) interior rel = {t_duty:.6e} -> {t_duty:.4e} (record {P26_T_DUTY}; "
         f"bound {INTERIOR_BAND:.0e}); balance {bal:.3e} (bound {E_BAND:.0e}); "
         f"(c3) miss {c3:.4e} (P26 {P26_C3:.4e}; asserted > {INTERIOR_BAND:.0e})")
    assert t_duty <= INTERIOR_BAND
    assert f"{t_duty:.4e}" == P26_T_DUTY
    assert bal <= E_BAND
    assert c3 > INTERIOR_BAND
