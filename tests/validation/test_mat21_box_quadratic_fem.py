"""`MAT-21` step 9, W2 — the gel box under the transverse quadratic source
q = 3Q(y² + z²)/(B² + C²) (the long-prism CP eddy-SAR shape), P2 on step 3's
m = 4 rung, steady and CN N = 40 to t* (§9 item 497, §10 chain P36; brief
`docs/planning/chunks/MAT-21.md`, review ruling 2026-10-16 10:30).

Sets no band and moves none; every band is imported by name.  Asserted:
(P0) constants; (Pdeg) 116 025 dofs on all three solves; (U) step 3's digits;
(S1-quad) steady interior vs the quadratic series ≤ ``INTERIOR_BAND``;
(S2-quad) Gauss ≤ ``S2_GAUSS_BAND``; (T-quad) CN interior at t* ≤
``INTERIOR_BAND``; (E-quad) CN balance ≤ ``E_BAND``.  Controls, PREDICTED
(10:30 ruling, rule (e)): (c-src), (c-src′), (c-z), (c-src-T) print their
measured value beside the predicted range and assert only > ``INTERIOR_BAND``.
"""

from __future__ import annotations

import time

import numpy as np
import pytest
import ufl
from mpi4py import MPI

from fem_em_solver.thermal import solve_steady_conduction, solve_transient_conduction
from fem_em_solver.utils.analytical import (
    robin_box_quadratic_steady_temperature,
    robin_box_quadratic_transient_temperature,
    robin_slab_eigenvalues,
    robin_slab_quadratic_coefficients,
)
from tests.complex_mode import complex_only
from tests.validation.test_mat10_sar_to_heat_sphere import ROBIN_TAG, _allsum, _say
from tests.validation.test_mat12_transient_sar_heating_sphere import E_BAND, RHO_C
from tests.validation.test_mat13_transient_perfusion_sphere import INTERIOR_BAND
from tests.validation.test_mat15_coil_driven_heating_phantom import (
    H_GEL,
    K_GEL,
    N_MAIN,
    T_INF,
)
from tests.validation.test_mat15_finite_cylinder_fem import S2_GAUSS_BAND, _max_imag
from tests.validation.test_mat21_box_closed_form import (
    A_BOX,
    B_BOX,
    C_BOX,
    HALF,
)
from tests.validation.test_mat21_box_closed_form import RHO_C as RHO_C_0
from tests.validation.test_mat21_box_closed_form import SERIES_TERMS, H, K, Q
from tests.validation.test_mat21_box_fem import T_STAR, V_BOX, _interior_chi, _steady_ex, _trans_ex
from tests.validation.test_mat21_box_fem_p2 import (
    CMP_DEGREE,
    DEGREE,
    QDEG_P,
    _dofs_pred,
    _exact_p,
    _rel_p,
    _rung_p,
)
from tests.validation.test_mat21_box_quadratic_closed_form import (
    QUAD_TERMS,
    WEIGHTS,
    WEIGHTS_NOZ,
)
from tests.validation.test_mat21_box_separable_fem import M4, STEP3_M4, _point_vals


class _Clock9:
    """Per-phase wall clock, printed on rank 0 with flush, prefix MAT21-9."""

    def __init__(self, tag):
        self.tag = tag
        self.t = time.perf_counter()

    def lap(self, phase):
        now = time.perf_counter()
        _say(f"MAT21-9 [{self.tag}] clock {phase}: {now - self.t:.1f} s")
        self.t = now


def _quad_ex(weights=WEIGHTS, n_terms=QUAD_TERMS):
    def f(x, y, z):
        return robin_box_quadratic_steady_temperature(
            x, y, z, half_widths=HALF, k=K, q=Q, h=H, weights=weights, t_inf=0.0, n_terms=n_terms
        )

    return f


def _quad_trans_ex(t, weights=WEIGHTS, n_terms=QUAD_TERMS):
    def f(x, y, z):
        return robin_box_quadratic_transient_temperature(
            x,
            y,
            z,
            t,
            half_widths=HALF,
            k=K,
            q=Q,
            h=H,
            rho_c=RHO_C_0,
            weights=weights,
            t_inf=0.0,
            n_terms=n_terms,
        )

    return f


def _quad_floor(n_terms=QUAD_TERMS, weights=WEIGHTS):
    S, S2 = [], []
    for L, n in zip(HALF, n_terms):
        gam = robin_slab_eigenvalues(H * L / K, n)
        c = 2.0 * np.sin(gam) / (gam + np.sin(gam) * np.cos(gam))
        e = robin_slab_quadratic_coefficients(gam)
        S.append(np.sum(c * np.sin(gam) / gam))
        S2.append(np.sum(e * np.sin(gam) / gam))
    w0, wx, wy, wz = weights
    tgt = (
        w0 * S[0] * S[1] * S[2]
        + wx * S2[0] * S[1] * S[2]
        + wy * S[0] * S2[1] * S[2]
        + wz * S[0] * S[1] * S2[2]
    )
    return 1.0 - tgt


@pytest.fixture(scope="module")
def quad_m4():
    comm = MPI.COMM_WORLD
    t0 = time.perf_counter()
    clk = _Clock9("quad")
    out = {"imag": 0.0}
    _say(
        f"\nMAT21-9 DEGREE = {DEGREE}, comparand P{CMP_DEGREE}, quadrature {QDEG_P}; m = {M4}; "
        f"weights {WEIGHTS}; QUAD_TERMS {QUAD_TERMS}; ranks {comm.size}"
    )

    # 1. Uniform control (U), step 3's exact code path (as steps 7 / 8).
    msh, ft, ns, n, v_h, st = _rung_p(M4, clk, DEGREE)
    ds = ufl.Measure("ds", domain=msh, subdomain_data=ft)(ROBIN_TAG)
    out["ndof_u"] = int(st.temperature.function_space.dofmap.index_map.size_global)
    out["imag"] = max(out["imag"], _max_imag(st.temperature))
    th_u = st.temperature - T_INF
    ex_u, npts = _exact_p(msh, _steady_ex())
    chi = _interior_chi(msh)
    clk.lap(f"uniform series interpolation ({npts} P3 points) + chi")
    out["u_whole"] = _rel_p(msh, th_u, ex_u)
    out["u_int"] = _rel_p(msh, th_u, ex_u, chi)
    _say(
        f"MAT21-9 (U) m = {M4} {tuple(ns)}: cells {n}; P2 dofs {out['ndof_u']}; whole "
        f"{out['u_whole']:.6e} (step 3 {STEP3_M4['whole']}); interior {out['u_int']:.6e} "
        f"(step 3 {STEP3_M4['interior']})"
    )
    clk.lap("uniform readings")

    # 2. Quadratic steady solve.
    xs = ufl.SpatialCoordinate(msh)
    q_quad = Q * 3.0 * (xs[1] * xs[1] + xs[2] * xs[2]) / (B_BOX**2 + C_BOX**2)
    sq = solve_steady_conduction(
        msh, ft, k=K_GEL, q=q_quad, h=H_GEL, t_inf=T_INF, robin_tag=ROBIN_TAG, degree=DEGREE
    )
    clk.lap("quadratic steady solve")
    out["ndof_q"] = int(sq.temperature.function_space.dofmap.index_map.size_global)
    out["imag"] = max(out["imag"], _max_imag(sq.temperature))
    th_q = sq.temperature - T_INF
    conv = float(np.real(_allsum(comm, H_GEL * th_q * ds)))
    sp = float(np.real(sq.source_power))
    out["gauss"] = abs(conv - sp) / sp
    _say(
        f"MAT21-9 (S2-quad) Gauss |h oint theta_h - source_power|/source_power = "
        f"{out['gauss']:.3e} (bound {S2_GAUSS_BAND:.0e}, predicted <= 5e-12); "
        f"source_power/(Q V_h) - 1 = {sp / (Q * v_h) - 1.0:+.3e} (predicted <= 1e-13)"
    )
    ex_q, npts_q = _exact_p(msh, _quad_ex())
    clk.lap(f"quadratic series interpolation ({npts_q} P3 points, {QUAD_TERMS})")
    out["q_int"] = _rel_p(msh, th_q, ex_q, chi)
    q_whole = _rel_p(msh, th_q, ex_q)
    ex_z, _ = _exact_p(msh, _quad_ex(WEIGHTS_NOZ, SERIES_TERMS))
    clk.lap(f"no-z series interpolation ({SERIES_TERMS})")
    out["c_src"] = _rel_p(msh, th_q, ex_u, chi)
    out["c_srcp"] = _rel_p(msh, th_u, ex_q, chi)
    out["c_z"] = _rel_p(msh, th_q, ex_z, chi)
    pv = _point_vals(
        sq.temperature, [[0.0, 0.0, 0.0], [0.0, 0.9 * B_BOX, 0.0], [0.9 * A_BOX, 0.0, 0.0]]
    )
    _say(
        f"MAT21-9 (S1-quad) steady interior vs quadratic series {out['q_int']:.6e} (bound "
        f"{INTERIOR_BAND:.0e}, predicted 1.5e-5 to 4e-5, range 8e-6 to 1.5e-4); whole {q_whole:.6e}; "
        f"theta_h(0) {pv[0]:.6f} K (series 4.276700), theta_h(0,0.9B,0) {pv[1]:.6f} K "
        f"(series 13.337014), theta_h(0.9A,0,0) {pv[2]:.6f} K (series 2.682231)"
    )
    _say(
        f"MAT21-9 control (c-src) PREDICTED: quadratic FEM vs uniform series {out['c_src']:.4e} "
        f"(predicted 0.449, range 0.43-0.47; asserted only > INTERIOR_BAND)"
    )
    _say(
        f"MAT21-9 control (c-src') PREDICTED: uniform FEM vs quadratic series {out['c_srcp']:.4e} "
        f"(predicted 0.514, range 0.50-0.53; asserted only > INTERIOR_BAND)"
    )
    _say(
        f"MAT21-9 control (c-z) PREDICTED: quadratic FEM vs no-z series {out['c_z']:.4e} "
        f"(predicted 0.0477, range 0.040-0.056; asserted only > INTERIOR_BAND; "
        f"ratio to band {out['c_z'] / INTERIOR_BAND:.2f})"
    )
    clk.lap("quadratic steady readings")

    # 3. CN N = N_MAIN to t*.
    res = solve_transient_conduction(
        msh,
        ft,
        k=K_GEL,
        rho_c=RHO_C,
        q=q_quad,
        h=H_GEL,
        t_inf=T_INF,
        t_init=T_INF,
        robin_tag=ROBIN_TAG,
        dt=T_STAR / N_MAIN,
        n_steps=N_MAIN,
        theta=0.5,
        startup_half_steps=2,
        degree=DEGREE,
    )
    clk.lap(f"CN N={N_MAIN} (quadratic source)")
    out["ndof_t"] = int(res.temperature.function_space.dofmap.index_map.size_global)
    out["imag"] = max(out["imag"], _max_imag(res.temperature))
    out["E"] = float(res.max_balance_residual)
    th_t = res.temperature - T_INF
    ex_t, _ = _exact_p(msh, _quad_trans_ex(T_STAR))
    ex_ut, _ = _exact_p(msh, _trans_ex(T_STAR))
    out["t_int"] = _rel_p(msh, th_t, ex_t, chi)
    t_whole = _rel_p(msh, th_t, ex_t)
    out["c_srcT"] = _rel_p(msh, th_t, ex_ut, chi)
    clk.lap("transient interpolation + readings")
    tv = _point_vals(res.temperature, [[0.0, 0.0, 0.0], [0.0, 0.9 * B_BOX, 0.0]])
    sp_t = getattr(res, "source_power", float("nan"))
    _say(
        f"MAT21-9 (T-quad) CN N={N_MAIN} at t*: interior vs quadratic series {out['t_int']:.6e} "
        f"(bound {INTERIOR_BAND:.0e}, predicted 1.2e-4 to 3e-4, range 5e-5 to 1e-3); whole "
        f"{t_whole:.6e}; theta_h(0, t*) {tv[0]:.6f} K (series 0.063466); theta_h(0,0.9B,0, t*) "
        f"{tv[1]:.6f} K (series 2.892661)"
    )
    _say(
        f"MAT21-9 (E-quad) CN max balance residual {out['E']:.3e} (bound {E_BAND:.0e}, predicted "
        f"<= 5e-13); transient source_power/(Q V_h) - 1 = "
        f"{float(np.real(sp_t)) / (Q * v_h) - 1.0:+.3e}"
    )
    _say(
        f"MAT21-9 control (c-src-T) PREDICTED: CN quadratic vs uniform transient series "
        f"{out['c_srcT']:.4e} (predicted 0.698, range 0.68-0.72; asserted only > INTERIOR_BAND)"
    )

    # 4.
    fl = _quad_floor()
    _say(
        f"MAT21-9 comparand floor 1 - target at QUAD_TERMS {fl:.4e} (drafting 1.2484e-6); "
        f"(S1-quad)/floor {out['q_int'] / fl:.1f}; max|Im theta| {out['imag']:.3e}; "
        f"V_BOX/V_h - 1 {V_BOX / v_h - 1.0:+.3e}; fixture elapsed {time.perf_counter() - t0:.1f} s"
    )
    return out


def test_precondition_constants():
    assert K == K_GEL and H == H_GEL and RHO_C_0 == RHO_C


@complex_only
def test_precondition_degree_dofs_m4(quad_m4):
    pred = _dofs_pred(M4)
    _say(
        f"MAT21-9 (Pdeg) P2 dofs uniform {quad_m4['ndof_u']} / quadratic {quad_m4['ndof_q']} / "
        f"CN {quad_m4['ndof_t']} (predicted {pred} = 116025)"
    )
    assert pred == 116025
    for key in ("ndof_u", "ndof_q", "ndof_t"):
        assert quad_m4[key] == pred


@complex_only
def test_u_uniform_control_reproduces_step3(quad_m4):
    assert f"{quad_m4['u_whole']:.6e}" == STEP3_M4["whole"]
    assert f"{quad_m4['u_int']:.6e}" == STEP3_M4["interior"]


@complex_only
def test_s1_quad_interior(quad_m4):
    assert quad_m4["q_int"] <= INTERIOR_BAND


@complex_only
def test_s2_quad_gauss_balance(quad_m4):
    assert quad_m4["gauss"] <= S2_GAUSS_BAND


@complex_only
def test_t_quad_transient_interior(quad_m4):
    assert quad_m4["t_int"] <= INTERIOR_BAND


@complex_only
def test_e_quad_cn_balance(quad_m4):
    assert quad_m4["E"] <= E_BAND


@complex_only
def test_control_c_src_uniform_series_misses(quad_m4):
    assert quad_m4["c_src"] > INTERIOR_BAND


@complex_only
def test_control_c_src_prime_quadratic_series_misses(quad_m4):
    assert quad_m4["c_srcp"] > INTERIOR_BAND


@complex_only
def test_control_c_z_no_z_series_misses(quad_m4):
    assert quad_m4["c_z"] > INTERIOR_BAND


@complex_only
def test_control_c_src_t_uniform_transient_misses(quad_m4):
    assert quad_m4["c_srcT"] > INTERIOR_BAND
