"""`MAT-21` step 11, W2 — P36's transverse quadratic source q = 3Q(y² + z²)/(B² + C²)
under `MAT-17`'s duty cycle re-timed to the box (step 2's schedule), P2 CN N = 40
to t* on step 3's m = 4 rung, against W1's one-tensor ten-run closed form (§9
item 524, §10 chain P38; brief `docs/planning/chunks/MAT-21.md`, review ruling
2026-10-18 03:00).

Sets no band and moves none; every band is imported by name.  Asserted:
(P0)/(P0′) pins; (Pdeg) 116 025 dofs on all four solves; (S0) census and
schedule string; (U) step 3's digits; (R36) the default run reproduces P36's
printed digits; (D′-quad) χ ≡ 1 vs default ≤ ``D_BAND``; (T-qduty) duty CN
interior vs the duty closed form ≤ ``INTERIOR_BAND``; (E-qduty) balance ≤
``E_BAND`` on the three CN runs.  Negative control (c-step), predicted (only
the separation > ``INTERIOR_BAND`` asserted): the duty CN against the
continuous series.
"""

from __future__ import annotations

import time

import numpy as np
import pytest
import ufl
from mpi4py import MPI

from fem_em_solver.thermal import solve_transient_conduction

from tests.complex_mode import complex_only
from tests.validation.test_mat10_sar_to_heat_sphere import _allsum, _say
from tests.validation.test_mat10_sar_to_heat_sphere import ROBIN_TAG
from tests.validation.test_mat12_transient_sar_heating_sphere import E_BAND, RHO_C
from tests.validation.test_mat13_transient_perfusion_sphere import INTERIOR_BAND
from tests.validation.test_mat15_coil_driven_heating_phantom import H_GEL, K_GEL, N_MAIN, T_INF
from tests.validation.test_mat15_finite_cylinder_fem import _max_imag
from tests.validation.test_mat17_duty_cycle_coil_driven import (
    D_BAND,
    DUTY,
    N_SWITCH,
    STEPS_PER_HALF,
)
from tests.validation.test_mat17_duty_cycle_fem import _census
from tests.validation.test_mat21_box_closed_form import B_BOX, C_BOX, Q, T_C
from tests.validation.test_mat21_box_closed_form import T_INF as T_INF_0
from tests.validation.test_mat21_box_duty_closed_form import DT_BOX, HALF_T, N_STEPS
from tests.validation.test_mat21_box_duty_closed_form import DUTY as DUTY_PIN
from tests.validation.test_mat21_box_duty_closed_form import N_SWITCH as N_SWITCH_PIN
from tests.validation.test_mat21_box_duty_closed_form import STEPS_PER_HALF as SPH_PIN
from tests.validation.test_mat21_box_duty_closed_form import T_STAR as T_STAR_2
from tests.validation.test_mat21_box_duty_fem import _chi_box
from tests.validation.test_mat21_box_fem import T_STAR, _interior_chi, _steady_ex
from tests.validation.test_mat21_box_fem_p2 import (
    CMP_DEGREE,
    DEGREE,
    QDEG_P,
    _dofs_pred,
    _exact_p,
    _rel_p,
    _rung_p,
)
from tests.validation.test_mat21_box_quadratic_closed_form import QUAD_TERMS, WEIGHTS
from tests.validation.test_mat21_box_quadratic_duty_closed_form import _quad_duty, _quad_duty_mean
from tests.validation.test_mat21_box_quadratic_fem import _quad_floor, _quad_trans_ex
from tests.validation.test_mat21_box_separable_fem import M4, STEP3_M4, _point_vals

SCHEDULE = "11111000011110000111100001111000011110000"  # step 2's (S0)
R36 = {"interior": "2.694340e-04", "whole": "4.310154e-04"}  # P36 W2 log :83


class _Clock11:
    """Per-phase wall clock, printed on rank 0 with flush, prefix MAT21-11."""

    def __init__(self, tag):
        self.tag = tag
        self.t = time.perf_counter()

    def lap(self, phase):
        now = time.perf_counter()
        _say(f"MAT21-11 [{self.tag}] clock {phase}: {now - self.t:.1f} s")
        self.t = now


def _cn(msh, ft, q_quad, schedule):
    kw = {} if schedule is None else {"source_schedule": schedule}
    return solve_transient_conduction(
        msh, ft, k=K_GEL, rho_c=RHO_C, q=q_quad, h=H_GEL, t_inf=T_INF, t_init=T_INF,
        robin_tag=ROBIN_TAG, dt=T_STAR / N_MAIN, n_steps=N_MAIN, theta=0.5,
        startup_half_steps=2, degree=DEGREE, **kw)


@pytest.fixture(scope="module")
def qduty_m4():
    comm = MPI.COMM_WORLD
    t0 = time.perf_counter()
    clk = _Clock11("qduty")
    out = {"imag": 0.0, "E": {}, "ndof": {}}
    _say(f"\nMAT21-11 DEGREE = {DEGREE}, comparand P{CMP_DEGREE}, quadrature {QDEG_P}; m = {M4}; "
         f"weights {WEIGHTS}; QUAD_TERMS {QUAD_TERMS}; ranks {comm.size}; DT_BOX {DT_BOX:.6f} s, "
         f"HALF_T {HALF_T:.4f} s")

    # 1. Uniform control (U), step 3's exact code path.
    msh, ft, ns, n, v_h, st = _rung_p(M4, clk, DEGREE)
    out["ndof"]["uniform"] = int(st.temperature.function_space.dofmap.index_map.size_global)
    out["imag"] = max(out["imag"], _max_imag(st.temperature))
    th_u = st.temperature - T_INF
    ex_u, npts = _exact_p(msh, _steady_ex())
    chi = _interior_chi(msh)
    clk.lap(f"uniform series interpolation ({npts} P3 points) + chi")
    out["u_whole"] = _rel_p(msh, th_u, ex_u)
    out["u_int"] = _rel_p(msh, th_u, ex_u, chi)
    _say(f"MAT21-11 (U) m = {M4} {tuple(ns)}: cells {n}; P2 dofs {out['ndof']['uniform']}; whole "
         f"{out['u_whole']:.6e} (step 3 {STEP3_M4['whole']}); interior {out['u_int']:.6e} "
         f"(step 3 {STEP3_M4['interior']})")
    del ex_u
    clk.lap("uniform readings")

    xs = ufl.SpatialCoordinate(msh)
    q_quad = Q * 3.0 * (xs[1] * xs[1] + xs[2] * xs[2]) / (B_BOX**2 + C_BOX**2)
    pts = [[0.0, 0.0, 0.0], [0.0, 0.9 * B_BOX, 0.0]]

    def _go(label, schedule):
        r = _cn(msh, ft, q_quad, schedule)
        clk.lap(f"CN N={N_MAIN} {label}")
        out["ndof"][label] = int(r.temperature.function_space.dofmap.index_map.size_global)
        out["imag"] = max(out["imag"], _max_imag(r.temperature))
        out["E"][label] = float(r.max_balance_residual)
        th = r.temperature - T_INF
        mean = float(np.real(_allsum(comm, th * ufl.dx(domain=msh)))) / v_h
        return r, th, mean, _point_vals(r.temperature, pts)

    # 2. Default schedule (P36's run), read against ex_t.
    r_def, th_def, m_def, pv_def = _go("default", None)
    ex_t, _ = _exact_p(msh, _quad_trans_ex(T_STAR))
    clk.lap(f"continuous series interpolation ({QUAD_TERMS})")
    out["r36_int"] = _rel_p(msh, th_def, ex_t, chi)
    out["r36_whole"] = _rel_p(msh, th_def, ex_t)
    _say(f"MAT21-11 (R36) default CN vs continuous quadratic series: interior {out['r36_int']:.6e} "
         f"(P36 {R36['interior']}), whole {out['r36_whole']:.6e} (P36 {R36['whole']})")

    # 3. chi = 1 (the 100 % duty run).
    r_one, th_one, m_one, pv_one = _go("chi1", lambda t: 1.0)
    out["dprime"] = _rel_p(msh, th_one, th_def)
    _say(f"MAT21-11 (D'-quad) chi = 1 vs default at t*: whole rel L2 {out['dprime']:.3e} "
         f"(bound D_BAND {D_BAND:.0e}, predicted <= 1e-14; P29 4.481e-15)")
    del r_def, r_one

    # 4. The duty run and its census.
    r_d, th_d, m_d, pv_d = _go("duty", _chi_box)
    out["sched"] = "".join(str(int(x)) for x in r_d.source_multipliers)
    _census(r_d, DT_BOX, 2, DUTY * T_STAR, "quad duty")
    out["census_ok"] = True

    # 5. The duty closed form.
    ex_d, npts_d = _exact_p(msh, lambda x, y, z: _quad_duty(x, y, z, T_STAR))
    clk.lap(f"duty series interpolation ({npts_d} P3 points, {QUAD_TERMS})")

    # 6. Readings.
    out["T"] = _rel_p(msh, th_d, ex_d, chi)
    out["T_whole"] = _rel_p(msh, th_d, ex_d)
    out["c_step"] = _rel_p(msh, th_d, ex_t, chi)
    clk.lap("duty readings")
    cf_mean = _quad_duty_mean(T_STAR)
    cf_c = float(_quad_duty(0.0, 0.0, 0.0, T_STAR))
    cf_b = float(_quad_duty(0.0, 0.9 * B_BOX, 0.0, T_STAR))
    _say(f"MAT21-11 (S0) schedule (2 half-steps, then 39 steps): {out['sched']} (must equal step 2's "
         f"{SCHEDULE})")
    _say(f"MAT21-11 (T-qduty) duty CN interior vs duty closed form at t* {out['T']:.6e} (bound "
         f"{INTERIOR_BAND:.0e}, predicted 2.2e-4, range 8e-5 to 1e-3); whole {out['T_whole']:.6e}")
    _say(f"MAT21-11 (E-qduty) max balance residual: " + ", ".join(
        f"{k_} {v:.3e}" for k_, v in out["E"].items()) + f" (bound {E_BAND:.0e}, predicted <= 5e-13)")
    _say(f"MAT21-11 control (c-step) PREDICTED: duty CN vs continuous series interior "
         f"{out['c_step']:.4e} (predicted 0.51, range 0.49-0.53; step 4 (c3) 0.50300; asserted only "
         f"> INTERIOR_BAND; ratio to band {out['c_step'] / INTERIOR_BAND:.1f})")
    _say(f"MAT21-11 record volume means at t*: <theta_duty> {m_d:.6f} K (closed form {cf_mean:.6f}), "
         f"<theta_chi1> {m_one:.6f} K, <theta_default> {m_def:.6f} K")
    _say(f"MAT21-11 record duty drop (FEM): <theta_duty>/<theta_chi1> {m_d / m_one:.5f}; centre "
         f"{pv_d[0] / pv_one[0]:.5f}; (0, 0.9B, 0) {pv_d[1] / pv_one[1]:.5f} (W1 closed form 0.49465 / "
         f"0.54466 / 0.49468; P26 uniform 0.49576 mean / 0.49902 centre)")
    _say(f"MAT21-11 record point values at t*: duty theta_h(0) {pv_d[0]:.6f} K (closed form "
         f"{cf_c:.6f}), theta_h(0,0.9B,0) {pv_d[1]:.6f} K (closed form {cf_b:.6f}); chi1 "
         f"{pv_one[0]:.6f} / {pv_one[1]:.6f} K; default {pv_def[0]:.6f} / {pv_def[1]:.6f} K")
    _say(f"MAT21-11 comparand floor _quad_floor() {_quad_floor():.4e} (drafting 1.2484e-6); "
         f"max|Im theta| {out['imag']:.3e}; fixture elapsed {time.perf_counter() - t0:.1f} s")
    return out


def test_precondition_pins():
    assert (N_STEPS, SPH_PIN, N_SWITCH_PIN, DUTY_PIN) == (N_MAIN, STEPS_PER_HALF, N_SWITCH, DUTY)
    assert T_STAR == 0.2 * T_C
    assert T_STAR == T_STAR_2
    assert T_INF == T_INF_0
    assert HALF_T == STEPS_PER_HALF * DT_BOX
    assert DT_BOX == T_STAR / N_MAIN
    assert DEGREE == 2 and CMP_DEGREE == 3


@complex_only
def test_precondition_degree_dofs_m4(qduty_m4):
    pred = _dofs_pred(M4)
    _say(f"MAT21-11 (Pdeg) P2 dofs {qduty_m4['ndof']} (predicted {pred} = 116025)")
    assert pred == 116025
    assert len(qduty_m4["ndof"]) == 4
    for key, v in qduty_m4["ndof"].items():
        assert v == pred, f"(Pdeg) {key}: {v}"


@complex_only
def test_s0_census_and_schedule(qduty_m4):
    assert qduty_m4["census_ok"]
    assert qduty_m4["sched"] == SCHEDULE


@complex_only
def test_u_uniform_control_reproduces_step3(qduty_m4):
    assert f"{qduty_m4['u_whole']:.6e}" == STEP3_M4["whole"]
    assert f"{qduty_m4['u_int']:.6e}" == STEP3_M4["interior"]


@complex_only
def test_r36_default_reproduces_p36(qduty_m4):
    assert f"{qduty_m4['r36_int']:.6e}" == R36["interior"]
    assert f"{qduty_m4['r36_whole']:.6e}" == R36["whole"]


@complex_only
def test_d_prime_quad_chi_one_is_default(qduty_m4):
    assert qduty_m4["dprime"] <= D_BAND


@complex_only
def test_t_qduty_interior_against_closed_form(qduty_m4):
    assert qduty_m4["T"] <= INTERIOR_BAND, f"(T-qduty) {qduty_m4['T']:.4e}"


@complex_only
def test_e_qduty_every_cn_balances(qduty_m4):
    assert len(qduty_m4["E"]) == 3
    for name, r in qduty_m4["E"].items():
        assert r <= E_BAND, f"(E-qduty) {name}: {r:.3e}"


@complex_only
def test_control_c_step_continuous_series_misses(qduty_m4):
    assert qduty_m4["c_step"] > INTERIOR_BAND
