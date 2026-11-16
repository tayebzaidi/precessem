"""`MAT-21` step 8, W2 — P32's separable source q(x) = Q·[1 − cos(πx/A)]
handed to the thermal solvers three ways (UFL expression, DG0 cell average,
P2 nodal interpolant) on step 3's m = 4 rung at P2; steady for all three, CN
N = 40 to t* for the P2 source (§9 item 483, §10 chain P35; brief
`docs/planning/chunks/MAT-21.md`, review ruling 2026-10-15 18:00).

Sets no band and moves none; every band is imported by name.  Asserted:
(P0) constants; (Pdeg) 116 025 dofs on all four solves; (U) step 3's digits;
(S6) the UFL-source run reproduces step 6 to the printed digit; (S7) the DG0
steady run reproduces step 7 to the printed digit; (S1-p2s) P2-source steady
interior vs P32's series ≤ ``INTERIOR_BAND``; (S2-p2s) Gauss ≤
``S2_GAUSS_BAND``; (T-p2s) CN interior at t* ≤ ``INTERIOR_BAND``; (E-p2s) CN
balance ≤ ``E_BAND``.  Control (predicted; separation only): (c-p2s) P2-source
run vs the uniform series > ``INTERIOR_BAND``.  Printed, NOT asserted (18:00
ruling, rule (e)): (c-ord) (S1-p2s) < (S1-dg0a), labelled predicted; a_h; a_p.
"""

from __future__ import annotations

import time

import numpy as np
import pytest
import ufl
from dolfinx import fem, la
from dolfinx import mesh as dmesh
from mpi4py import MPI

from fem_em_solver.thermal import solve_steady_conduction, solve_transient_conduction

from tests.complex_mode import complex_only
from tests.validation.test_mat10_sar_to_heat_sphere import ROBIN_TAG, _allsum, _say
from tests.validation.test_mat12_transient_sar_heating_sphere import E_BAND, RHO_C
from tests.validation.test_mat13_transient_perfusion_sphere import INTERIOR_BAND
from tests.validation.test_mat15_coil_driven_heating_phantom import (
    H_GEL,
    K_GEL,
    N_MAIN,
    T_INF,
    _cell_average_dg0_tag,
)
from tests.validation.test_mat15_finite_cylinder_fem import S2_GAUSS_BAND, _max_imag
from tests.validation.test_mat21_box_closed_form import A_BOX, H, K, Q
from tests.validation.test_mat21_box_closed_form import RHO_C as RHO_C_0
from tests.validation.test_mat21_box_dg0_source_fem import STEP6_S1, _sep_ex_a
from tests.validation.test_mat21_box_fem import T_STAR, V_BOX, _interior_chi, _steady_ex
from tests.validation.test_mat21_box_fem_p2 import (
    CMP_DEGREE,
    DEGREE,
    QDEG_P,
    _dofs_pred,
    _exact_p,
    _rel_p,
    _rung_p,
)
from tests.validation.test_mat21_box_separable_fem import (
    M4,
    PROFILES,
    SEP_TERMS,
    STEP3_M4,
    _point_vals,
    _sep_ex,
    _sep_floor,
    _sep_trans_ex,
)

# Step 7 W2 printed digits (`20261015T200342Z_MAT-21-step7-W2-dg0-source.log:81–82`).
STEP7_S1 = {"dg0_int": "3.419236e-04", "dg0a_int": "1.179001e-04"}


class _Clock8:
    """Per-phase wall clock, printed on rank 0 with flush, prefix MAT21-8."""

    def __init__(self, tag):
        self.tag = tag
        self.t = time.perf_counter()

    def lap(self, phase):
        now = time.perf_counter()
        _say(f"MAT21-8 [{self.tag}] clock {phase}: {now - self.t:.1f} s")
        self.t = now


@pytest.fixture(scope="module")
def p2s_m4():
    comm = MPI.COMM_WORLD
    t0 = time.perf_counter()
    clk = _Clock8("p2s")
    out = {"imag": 0.0}
    _say(f"\nMAT21-8 DEGREE = {DEGREE}, comparand P{CMP_DEGREE}, quadrature {QDEG_P}; m = {M4}; "
         f"profiles {PROFILES}; SEP_TERMS {SEP_TERMS}; ranks {comm.size}")

    # (i) Uniform control (U), step 3's exact code path (as step 7).
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
    _say(f"MAT21-8 (U) m = {M4} {tuple(ns)}: cells {n}; P2 dofs {out['ndof_u']}; whole "
         f"{out['u_whole']:.6e} (step 3 {STEP3_M4['whole']}); interior {out['u_int']:.6e} "
         f"(step 3 {STEP3_M4['interior']})")
    clk.lap("uniform readings")

    # (ii) UFL-source separable steady solve, exactly as steps 6 / 7.
    xs = ufl.SpatialCoordinate(msh)
    q_sep = Q * (1.0 - ufl.cos(ufl.pi * xs[0] / A_BOX))
    ss = solve_steady_conduction(msh, ft, k=K_GEL, q=q_sep, h=H_GEL, t_inf=T_INF,
                                 robin_tag=ROBIN_TAG, degree=DEGREE)
    clk.lap("UFL-source steady solve")
    out["ndof_s"] = int(ss.temperature.function_space.dofmap.index_map.size_global)
    out["imag"] = max(out["imag"], _max_imag(ss.temperature))
    th_s = ss.temperature - T_INF
    ex_s, npts_s = _exact_p(msh, _sep_ex())
    clk.lap(f"separable series interpolation ({npts_s} P3 points, {SEP_TERMS})")
    out["s6_int"] = _rel_p(msh, th_s, ex_s, chi)
    out["s6_whole"] = _rel_p(msh, th_s, ex_s)
    _say(f"MAT21-8 (S6) UFL-source steady interior {out['s6_int']:.6e} (step 6 "
         f"{STEP6_S1['interior']}); whole {out['s6_whole']:.6e} (step 6 {STEP6_S1['whole']})")
    clk.lap("UFL-source readings")

    # (iii) The DG0 hand-off and its steady solve, exactly as step 7 (no DG0 CN leg).
    tdim = msh.topology.dim
    n_loc = msh.topology.index_map(tdim).size_local
    ct = dmesh.meshtags(msh, tdim, np.arange(n_loc, dtype=np.int32),
                        np.ones(n_loc, dtype=np.int32))
    W, cbar = _cell_average_dg0_tag(msh, ct, 1, ufl.cos(ufl.pi * xs[0] / A_BOX))
    q_h = fem.Function(W)
    q_h.x.array[:] = Q * (1.0 - cbar)
    q_h.x.scatter_forward()
    w = ufl.TestFunction(W)
    one = fem.Constant(msh, np.asarray(1.0, dtype=q_h.x.array.dtype))
    bvol = fem.assemble_vector(fem.form(ufl.inner(one, w) * ufl.dx))
    bvol.scatter_reverse(la.InsertMode.add)
    n_own = W.dofmap.index_map.size_local
    vol = np.real(bvol.array[:n_own])
    cb = np.real(cbar[:n_own])
    a_loc = float(np.sum(vol * cb * cb))
    a_h = comm.allreduce(a_loc, op=MPI.SUM) / (V_BOX / 2.0)
    out["a_h"] = a_h
    _say(f"MAT21-8 a_h = {a_h:.9f} (step 7 0.999391791; printed, not asserted); "
         f"1 - a_h = {1.0 - a_h:.4e}")
    clk.lap("DG0 cell average + a_h")

    sd = solve_steady_conduction(msh, ft, k=K_GEL, q=q_h, h=H_GEL, t_inf=T_INF,
                                 robin_tag=ROBIN_TAG, degree=DEGREE)
    clk.lap("DG0-source steady solve")
    out["ndof_d"] = int(sd.temperature.function_space.dofmap.index_map.size_global)
    out["imag"] = max(out["imag"], _max_imag(sd.temperature))
    th_d = sd.temperature - T_INF
    ex_a, npts_a = _exact_p(msh, _sep_ex_a(a_h))
    clk.lap(f"a_h series interpolation ({npts_a} P3 points)")
    out["d_int"] = _rel_p(msh, th_d, ex_s, chi)
    out["da_int"] = _rel_p(msh, th_d, ex_a, chi)
    _say(f"MAT21-8 (S7) DG0 steady interior vs P32 series {out['d_int']:.6e} (step 7 "
         f"{STEP7_S1['dg0_int']}); vs a_h series {out['da_int']:.6e} (step 7 "
         f"{STEP7_S1['dg0a_int']})")
    clk.lap("DG0 steady readings")

    # (iv) The P2 source: nodal interpolant of q on the solution space.
    V2 = fem.functionspace(msh, ("Lagrange", DEGREE))
    q_p = fem.Function(V2)
    q_p.interpolate(lambda x: Q * (1.0 - np.cos(np.pi * x[0] / A_BOX)))
    q_p.x.scatter_forward()
    dxq = ufl.dx(domain=msh, metadata=QDEG_P)
    a_p = float(np.real(_allsum(comm, (1.0 - q_p / Q) * ufl.cos(ufl.pi * xs[0] / A_BOX) * dxq))) \
        / (V_BOX / 2.0)
    out["a_p"] = a_p
    _say(f"MAT21-8 a_p = {a_p:.9f}; 1 - a_p = {1.0 - a_p:+.4e} (predicted |1 - a_p| <= 1e-5; "
         f"printed, not asserted)")
    clk.lap("P2 interpolation + a_p")

    # (v) P2-source steady solve and readings.
    sp2 = solve_steady_conduction(msh, ft, k=K_GEL, q=q_p, h=H_GEL, t_inf=T_INF,
                                  robin_tag=ROBIN_TAG, degree=DEGREE)
    clk.lap("P2-source steady solve")
    out["ndof_p"] = int(sp2.temperature.function_space.dofmap.index_map.size_global)
    out["imag"] = max(out["imag"], _max_imag(sp2.temperature))
    th_p = sp2.temperature - T_INF
    conv = float(np.real(_allsum(comm, H_GEL * th_p * ds)))
    sp = float(np.real(sp2.source_power))
    out["gauss_p"] = abs(conv - sp) / sp
    _say(f"MAT21-8 (S2-p2s) Gauss |h oint theta_h - source_power|/source_power = "
         f"{out['gauss_p']:.3e} (bound {S2_GAUSS_BAND:.0e}, predicted <= 5e-12); "
         f"source_power/(Q V_h) - 1 = {sp / (Q * v_h) - 1.0:+.3e} (predicted round-off)")
    out["p_int"] = _rel_p(msh, th_p, ex_s, chi)
    p_whole = _rel_p(msh, th_p, ex_s)
    pa_int = _rel_p(msh, th_p, ex_a, chi)
    out["c_p2s"] = _rel_p(msh, th_p, ex_u, chi)
    pv = _point_vals(sp2.temperature, [[0.0, 0.0, 0.0], [0.9 * A_BOX, 0.0, 0.0]])
    _say(f"MAT21-8 (S1-p2s) P2-source steady interior vs P32 series {out['p_int']:.6e} (bound "
         f"{INTERIOR_BAND:.0e}, predicted 3.2e-5 to 4.5e-5, range 2e-5 to 8e-5); whole "
         f"{p_whole:.6e}; vs a_h series {pa_int:.6e} (printed); theta_h(0) {pv[0]:.6f} K "
         f"(series 2.741147), theta_h(0.9A,0,0) {pv[1]:.6f} K (series 14.919189)")
    _say(f"MAT21-8 control (c-p2s) P2-source run vs uniform series {out['c_p2s']:.4e} "
         f"(predicted 0.512, range 0.50-0.53; asserted only > INTERIOR_BAND)")
    ratio = out["da_int"] / out["p_int"]
    _say(f"MAT21-8 (c-ord) PREDICTED, printed not asserted: (S1-p2s) {out['p_int']:.6e} < "
         f"(S1-dg0a) {out['da_int']:.6e}? {out['p_int'] < out['da_int']}; ratio (S1-dg0a)/(S1-p2s) "
         f"{ratio:.3f} (predicted approx 3.4, ceiling approx 3.4)")
    clk.lap("P2-source steady readings")

    # (vi) CN N = N_MAIN to t* with the P2 source.
    res = solve_transient_conduction(
        msh, ft, k=K_GEL, rho_c=RHO_C, q=q_p, h=H_GEL, t_inf=T_INF, t_init=T_INF,
        robin_tag=ROBIN_TAG, dt=T_STAR / N_MAIN, n_steps=N_MAIN, theta=0.5,
        startup_half_steps=2, degree=DEGREE)
    clk.lap(f"CN N={N_MAIN} (P2 source)")
    out["imag"] = max(out["imag"], _max_imag(res.temperature))
    out["E"] = float(res.max_balance_residual)
    th_t = res.temperature - T_INF
    ex_t, _ = _exact_p(msh, _sep_trans_ex(T_STAR))
    out["t_int"] = _rel_p(msh, th_t, ex_t, chi)
    t_whole = _rel_p(msh, th_t, ex_t)
    clk.lap("transient interpolation + readings")
    tc = _point_vals(res.temperature, np.zeros((1, 3)))[0]
    sp_t = getattr(res, "source_power", float("nan"))
    _say(f"MAT21-8 (T-p2s) CN N={N_MAIN} at t*: interior vs P32 series {out['t_int']:.6e} "
         f"(bound {INTERIOR_BAND:.0e}, predicted 1.39e-4 to 1.7e-4, range 1.2e-4 to 3e-4); whole "
         f"{t_whole:.6e}; theta_h(0, t*) {tc:.6f} K (series 0.024797); beside step 7 (T-dg0) "
         f"5.299989e-4 (ratio {5.299989e-4 / out['t_int']:.3f}) and its a_h reading 3.825551e-4 "
         f"(ratio {3.825551e-4 / out['t_int']:.3f})")
    _say(f"MAT21-8 (E-p2s) CN max balance residual {out['E']:.3e} (bound {E_BAND:.0e}, predicted "
         f"<= 5e-13); transient source_power/(Q V_h) - 1 = "
         f"{float(np.real(sp_t)) / (Q * v_h) - 1.0:+.3e}")

    # (vii) Source-path differences, printed.
    dp = _rel_p(msh, th_p, th_s, chi)
    dd = _rel_p(msh, th_d, th_s, chi)
    _say(f"MAT21-8 source-path differences (interior, rel. to UFL run): P2 - UFL {dp:.6e} "
         f"(predicted <= 1e-5); DG0 - UFL {dd:.6e} (predicted approx 3.4e-4)")
    clk.lap("source-path differences")

    # (viii)
    fl = _sep_floor()
    _say(f"MAT21-8 comparand floor 1 - prod S' at SEP_TERMS {fl:.4e} (step 6 2.2196e-6); "
         f"max|Im theta| {out['imag']:.3e}; fixture elapsed {time.perf_counter() - t0:.1f} s")
    return out


def test_precondition_constants():
    assert K == K_GEL and H == H_GEL and RHO_C_0 == RHO_C


@complex_only
def test_precondition_degree_dofs_m4(p2s_m4):
    pred = _dofs_pred(M4)
    _say(f"MAT21-8 (Pdeg) P2 dofs uniform {p2s_m4['ndof_u']} / UFL {p2s_m4['ndof_s']} / DG0 "
         f"{p2s_m4['ndof_d']} / P2-source {p2s_m4['ndof_p']} (predicted {pred} = 116025)")
    assert pred == 116025
    for key in ("ndof_u", "ndof_s", "ndof_d", "ndof_p"):
        assert p2s_m4[key] == pred


@complex_only
def test_u_uniform_control_reproduces_step3(p2s_m4):
    assert f"{p2s_m4['u_whole']:.6e}" == STEP3_M4["whole"]
    assert f"{p2s_m4['u_int']:.6e}" == STEP3_M4["interior"]


@complex_only
def test_s6_ufl_source_reproduces_step6(p2s_m4):
    assert f"{p2s_m4['s6_int']:.6e}" == STEP6_S1["interior"]
    assert f"{p2s_m4['s6_whole']:.6e}" == STEP6_S1["whole"]


@complex_only
def test_s7_dg0_source_reproduces_step7(p2s_m4):
    assert f"{p2s_m4['d_int']:.6e}" == STEP7_S1["dg0_int"]
    assert f"{p2s_m4['da_int']:.6e}" == STEP7_S1["dg0a_int"]


@complex_only
def test_s1_p2s_interior(p2s_m4):
    assert p2s_m4["p_int"] <= INTERIOR_BAND


@complex_only
def test_s2_p2s_gauss_balance(p2s_m4):
    assert p2s_m4["gauss_p"] <= S2_GAUSS_BAND


@complex_only
def test_t_p2s_transient_interior(p2s_m4):
    assert p2s_m4["t_int"] <= INTERIOR_BAND


@complex_only
def test_e_p2s_cn_balance(p2s_m4):
    assert p2s_m4["E"] <= E_BAND


@complex_only
def test_control_cp2s_uniform_series_misses(p2s_m4):
    assert p2s_m4["c_p2s"] > INTERIOR_BAND
