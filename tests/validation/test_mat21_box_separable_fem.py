"""`MAT-21` step 6, W2 — P2 FEM on step 3's m = 4 rung under the separable
source q(x) = Q·[1 − cos(πx/A)] (§9 item 467, §10 chain P32; brief
`docs/planning/chunks/MAT-21.md`), steady and CN N = 40 to t*.

Sets no band and moves none; every band is imported by name as step 3 does.
Asserted: (P0) constants; (Pdeg) global P2 dofs = 116 025; (U) the uniform
control's whole / interior readings to step 3's printed digit; (S1-sep)
steady interior ≤ ``INTERIOR_BAND``; (S2-sep) Gauss ≤ ``S2_GAUSS_BAND``;
(T-sep) CN interior at t* ≤ ``INTERIOR_BAND``; (E-sep) CN balance ≤
``E_BAND``.  Controls (predicted; separations only): (c-src), (c-src′),
(c-src-T) > ``INTERIOR_BAND``.

The cosine is a full period on [−A, A]; on ``create_box``'s x-periodic split
both automatically estimated source quadratures integrate it to zero, so the
balances stay at round-off.  A half-period profile loses this.
"""

from __future__ import annotations

import time

import numpy as np
import pytest
import ufl
from mpi4py import MPI

from fem_em_solver.post.evaluation import evaluate_vector_field_parallel
from fem_em_solver.thermal import solve_steady_conduction, solve_transient_conduction
from fem_em_solver.utils.analytical import (
    robin_box_separable_steady_temperature,
    robin_box_separable_transient_temperature,
    robin_slab_cosine_coefficients,
    robin_slab_eigenvalues,
)

from tests.complex_mode import complex_only
from tests.validation.test_mat10_sar_to_heat_sphere import ROBIN_TAG, _allsum, _say
from tests.validation.test_mat12_transient_sar_heating_sphere import E_BAND, RHO_C
from tests.validation.test_mat13_transient_perfusion_sphere import INTERIOR_BAND
from tests.validation.test_mat15_coil_driven_heating_phantom import H_GEL, K_GEL, N_MAIN, T_INF
from tests.validation.test_mat15_finite_cylinder_fem import S2_GAUSS_BAND, _max_imag
from tests.validation.test_mat21_box_closed_form import A_BOX, HALF, H, K, Q
from tests.validation.test_mat21_box_closed_form import RHO_C as RHO_C_0
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

PROFILES = ((1.0, -1.0, 1), (1.0, 0.0, 0), (1.0, 0.0, 0))
SEP_TERMS = (80, 40, 40)
M4 = 4
STEP3_M4 = {"whole": "2.489307e-05", "interior": "1.669359e-05"}  # step 3 W1 :105-106


class _Clock6:
    """Per-phase wall clock, printed on rank 0 with flush, prefix MAT21-6."""

    def __init__(self, tag):
        self.tag = tag
        self.t = time.perf_counter()

    def lap(self, phase):
        now = time.perf_counter()
        _say(f"MAT21-6 [{self.tag}] clock {phase}: {now - self.t:.1f} s")
        self.t = now


def _sep_ex():
    def f(x, y, z):
        return robin_box_separable_steady_temperature(
            x, y, z, half_widths=HALF, k=K, q=Q, h=H, profiles=PROFILES, t_inf=0.0,
            n_terms=SEP_TERMS)
    return f


def _sep_trans_ex(t):
    def f(x, y, z):
        return robin_box_separable_transient_temperature(
            x, y, z, t, half_widths=HALF, k=K, q=Q, h=H, rho_c=RHO_C_0, profiles=PROFILES,
            t_inf=0.0, n_terms=SEP_TERMS)
    return f


def _sep_floor():
    out = 1.0
    for L, n, (a, b, j) in zip(HALF, SEP_TERMS, PROFILES):
        gam = robin_slab_eigenvalues(H * L / K, n)
        c = 2.0 * np.sin(gam) / (gam + np.sin(gam) * np.cos(gam))
        p = c if b == 0.0 else a * c + b * robin_slab_cosine_coefficients(gam, j)
        out *= np.sum(p * np.sin(gam) / gam)
    return 1.0 - out


def _point_vals(T, pts):
    vals, ok = evaluate_vector_field_parallel(T, np.asarray(pts, dtype=float))
    return [float(np.real(np.ravel(vals[i])[0])) - T_INF if bool(ok[i]) else float("nan")
            for i in range(len(pts))]


@pytest.fixture(scope="module")
def sep_m4():
    comm = MPI.COMM_WORLD
    t0 = time.perf_counter()
    clk = _Clock6("sep")
    out = {"imag": 0.0}
    _say(f"\nMAT21-6 DEGREE = {DEGREE}, comparand P{CMP_DEGREE}, quadrature {QDEG_P}; m = {M4}; "
         f"profiles {PROFILES}; SEP_TERMS {SEP_TERMS}")

    # 1. Uniform control (U), step 3's exact code path.
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
    gauss_u = abs(float(np.real(_allsum(comm, H_GEL * th_u * ds))) - Q * v_h) / (Q * v_h)
    thc_u = _point_vals(st.temperature, np.zeros((1, 3)))[0]
    _say(f"MAT21-6 (U) m = {M4} {tuple(ns)}: cells {n}; P2 dofs {out['ndof_u']}; whole "
         f"{out['u_whole']:.6e} (step 3 {STEP3_M4['whole']}); interior {out['u_int']:.6e} "
         f"(step 3 {STEP3_M4['interior']}); theta_h(0) {thc_u:.6f} K (step 3 12.003960, printed); "
         f"Gauss {gauss_u:.3e} (printed)")
    clk.lap("uniform readings")

    # 2. Separable steady solve.
    xs = ufl.SpatialCoordinate(msh)
    q_sep = Q * (1.0 - ufl.cos(ufl.pi * xs[0] / A_BOX))
    ss = solve_steady_conduction(msh, ft, k=K_GEL, q=q_sep, h=H_GEL, t_inf=T_INF,
                                 robin_tag=ROBIN_TAG, degree=DEGREE)
    clk.lap("separable steady solve")
    out["ndof_s"] = int(ss.temperature.function_space.dofmap.index_map.size_global)
    out["imag"] = max(out["imag"], _max_imag(ss.temperature))
    th_s = ss.temperature - T_INF
    conv = float(np.real(_allsum(comm, H_GEL * th_s * ds)))
    out["gauss_s"] = abs(conv - Q * v_h) / (Q * v_h)
    sp = float(np.real(ss.source_power))
    _say(f"MAT21-6 (S2-sep) Gauss |h oint theta_h - Q V_h|/(Q V_h) = {out['gauss_s']:.3e} "
         f"(bound {S2_GAUSS_BAND:.0e}, predicted <= 5e-12); source_power/(Q V_h) - 1 = "
         f"{sp / (Q * v_h) - 1.0:+.3e} (predicted <= 1e-13); V_h/V - 1 = {v_h / V_BOX - 1.0:+.3e}")
    ex_s, npts_s = _exact_p(msh, _sep_ex())
    clk.lap(f"separable series interpolation ({npts_s} P3 points, {SEP_TERMS})")
    out["s_int"] = _rel_p(msh, th_s, ex_s, chi)
    out["s_whole"] = _rel_p(msh, th_s, ex_s)
    pv = _point_vals(ss.temperature, [[0.0, 0.0, 0.0], [0.9 * A_BOX, 0.0, 0.0]])
    out["c_src"] = _rel_p(msh, th_s, ex_u, chi)
    out["c_src_p"] = _rel_p(msh, th_u, ex_s, chi)
    _say(f"MAT21-6 (S1-sep) steady interior {out['s_int']:.6e} (bound {INTERIOR_BAND:.0e}, "
         f"predicted 1.5e-5 to 4e-5, range 8e-6 to 1.5e-4); whole {out['s_whole']:.6e} "
         f"(predicted 2.5e-5 to 6e-5); theta_h(0) {pv[0]:.6f} K (series 2.741147), "
         f"theta_h(0.9A,0,0) {pv[1]:.6f} K (series 14.919189)")
    _say(f"MAT21-6 control (c-src) sep FEM vs uniform series: {out['c_src']:.4e} (predicted 0.513, "
         f"range 0.50-0.53); (c-src') uniform FEM vs sep series: {out['c_src_p']:.4e} "
         f"(predicted 0.527, range 0.51-0.54)")
    clk.lap("separable steady readings + crosses")

    # 3. CN N = N_MAIN to t*.
    res = solve_transient_conduction(
        msh, ft, k=K_GEL, rho_c=RHO_C, q=q_sep, h=H_GEL, t_inf=T_INF, t_init=T_INF,
        robin_tag=ROBIN_TAG, dt=T_STAR / N_MAIN, n_steps=N_MAIN, theta=0.5,
        startup_half_steps=2, degree=DEGREE)
    clk.lap(f"CN N={N_MAIN}")
    out["imag"] = max(out["imag"], _max_imag(res.temperature))
    out["E"] = float(res.max_balance_residual)
    th_t = res.temperature - T_INF
    ex_t, _ = _exact_p(msh, _sep_trans_ex(T_STAR))
    out["t_int"] = _rel_p(msh, th_t, ex_t, chi)
    t_whole = _rel_p(msh, th_t, ex_t)
    ex_ut, _ = _exact_p(msh, _trans_ex(T_STAR))
    out["c_src_t"] = _rel_p(msh, th_t, ex_ut, chi)
    clk.lap("transient interpolations + readings")
    tc = _point_vals(res.temperature, np.zeros((1, 3)))[0]
    sp_t = getattr(res, "source_power", float("nan"))
    _say(f"MAT21-6 (T-sep) CN N={N_MAIN} at t*: interior {out['t_int']:.6e} (bound "
         f"{INTERIOR_BAND:.0e}, predicted 1.3e-4 to 4e-4; step 3 uniform 1.322904e-4); whole "
         f"{t_whole:.6e}; theta_h(0, t*) {tc:.6f} K (series 0.024797)")
    _say(f"MAT21-6 (E-sep) CN max balance residual {out['E']:.3e} (bound {E_BAND:.0e}, predicted "
         f"<= 5e-13); transient source_power/(Q V_h) - 1 = "
         f"{float(np.real(sp_t)) / (Q * v_h) - 1.0:+.3e}")
    _say(f"MAT21-6 control (c-src-T) CN sep vs uniform transient series: {out['c_src_t']:.4e} "
         f"(predicted 0.657, range 0.64-0.67)")

    # 4. Floors and imaginary parts.
    fl = _sep_floor()
    _say(f"MAT21-6 comparand floor 1 - prod S' at SEP_TERMS {fl:.4e} (predicted 2.2196e-6); "
         f"ratios steady {out['s_int'] / fl:.2f}, transient {out['t_int'] / fl:.2f}; max|Im theta| "
         f"{out['imag']:.3e}; fixture elapsed {time.perf_counter() - t0:.1f} s")
    return out


def test_precondition_constants():
    assert K == K_GEL and H == H_GEL and RHO_C_0 == RHO_C


@complex_only
def test_precondition_degree_dofs_m4(sep_m4):
    pred = _dofs_pred(M4)
    _say(f"MAT21-6 (Pdeg) P2 dofs uniform {sep_m4['ndof_u']} / sep {sep_m4['ndof_s']} "
         f"(predicted {pred} = 116025)")
    assert pred == 116025
    assert sep_m4["ndof_u"] == pred and sep_m4["ndof_s"] == pred


@complex_only
def test_u_uniform_control_reproduces_step3(sep_m4):
    assert f"{sep_m4['u_whole']:.6e}" == STEP3_M4["whole"]
    assert f"{sep_m4['u_int']:.6e}" == STEP3_M4["interior"]


@complex_only
def test_s1_sep_interior(sep_m4):
    assert sep_m4["s_int"] <= INTERIOR_BAND


@complex_only
def test_s2_sep_gauss_balance(sep_m4):
    assert sep_m4["gauss_s"] <= S2_GAUSS_BAND


@complex_only
def test_t_sep_transient_interior(sep_m4):
    assert sep_m4["t_int"] <= INTERIOR_BAND


@complex_only
def test_e_sep_cn_balance(sep_m4):
    assert sep_m4["E"] <= E_BAND


@complex_only
def test_control_csrc_uniform_series_misses(sep_m4):
    assert sep_m4["c_src"] > INTERIOR_BAND


@complex_only
def test_control_csrc_prime_sep_series_misses(sep_m4):
    assert sep_m4["c_src_p"] > INTERIOR_BAND


@complex_only
def test_control_csrc_t_uniform_transient_misses(sep_m4):
    assert sep_m4["c_src_t"] > INTERIOR_BAND
