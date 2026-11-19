"""`MAT-21` step 7, W2 — P32's separable source q(x) = Q·[1 − cos(πx/A)]
through the DG0 cell-average hand-off coil SAR uses (`MAT-10` / `MAT-12` /
`MAT-15`), on step 3's m = 4 rung at P2, steady and CN N = 40 to t*
(§9 item 481, §10 chain P34; brief `docs/planning/chunks/MAT-21.md`).

Sets no band and moves none; every band is imported by name.  Asserted:
(P0) constants; (Pdeg) 116 025 dofs; (U) step 3's digits; (S6) the UFL-source
run reproduces step 6 to the printed digit; (S1-dg0) DG0 steady interior vs
P32's series ≤ ``INTERIOR_BAND``; (S1-dg0a) the same vs the series at the DG0
source's cosine amplitude a_h ≤ ``INTERIOR_BAND``; (S2-dg0) Gauss ≤
``S2_GAUSS_BAND``; (T-dg0) CN interior at t* ≤ ``INTERIOR_BAND``; (E-dg0) CN
balance ≤ ``E_BAND``.  Controls (predicted; separation / ordering only):
(c-dg0) DG0 run vs the uniform series > ``INTERIOR_BAND``; (c-a) UFL run vs
the a_h series > (S1-dg0a).  a_h is printed, never asserted.
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
from fem_em_solver.utils.analytical import (
    robin_box_separable_steady_temperature,
    robin_box_separable_transient_temperature,
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
    _cell_average_dg0_tag,
)
from tests.validation.test_mat15_finite_cylinder_fem import S2_GAUSS_BAND, _max_imag
from tests.validation.test_mat21_box_closed_form import A_BOX, HALF
from tests.validation.test_mat21_box_closed_form import RHO_C as RHO_C_0
from tests.validation.test_mat21_box_closed_form import H, K, Q
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

# Step 6 W2 printed digits (`20261015T130420Z_MAT-21-step6-W2-separable-fem.log:99`).
STEP6_S1 = {"interior": "3.486633e-05", "whole": "4.327217e-05"}


class _Clock7:
    """Per-phase wall clock, printed on rank 0 with flush, prefix MAT21-7."""

    def __init__(self, tag):
        self.tag = tag
        self.t = time.perf_counter()

    def lap(self, phase):
        now = time.perf_counter()
        _say(f"MAT21-7 [{self.tag}] clock {phase}: {now - self.t:.1f} s")
        self.t = now


def _profiles_a(a_h):
    return ((1.0, -a_h, PROFILES[0][2]), PROFILES[1], PROFILES[2])


def _sep_ex_a(a_h):
    prof = _profiles_a(a_h)

    def f(x, y, z):
        return robin_box_separable_steady_temperature(
            x, y, z, half_widths=HALF, k=K, q=Q, h=H, profiles=prof, t_inf=0.0, n_terms=SEP_TERMS
        )

    return f


def _sep_trans_ex_a(a_h, t):
    prof = _profiles_a(a_h)

    def f(x, y, z):
        return robin_box_separable_transient_temperature(
            x,
            y,
            z,
            t,
            half_widths=HALF,
            k=K,
            q=Q,
            h=H,
            rho_c=RHO_C_0,
            profiles=prof,
            t_inf=0.0,
            n_terms=SEP_TERMS,
        )

    return f


@pytest.fixture(scope="module")
def dg0_m4():
    comm = MPI.COMM_WORLD
    t0 = time.perf_counter()
    clk = _Clock7("dg0")
    out = {"imag": 0.0}
    _say(
        f"\nMAT21-7 DEGREE = {DEGREE}, comparand P{CMP_DEGREE}, quadrature {QDEG_P}; m = {M4}; "
        f"profiles {PROFILES}; SEP_TERMS {SEP_TERMS}; ranks {comm.size}"
    )

    # (i) Uniform control (U), step 3's exact code path.
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
        f"MAT21-7 (U) m = {M4} {tuple(ns)}: cells {n}; P2 dofs {out['ndof_u']}; whole "
        f"{out['u_whole']:.6e} (step 3 {STEP3_M4['whole']}); interior {out['u_int']:.6e} "
        f"(step 3 {STEP3_M4['interior']})"
    )
    clk.lap("uniform readings")

    # (ii) UFL-source separable steady solve, exactly as step 6.
    xs = ufl.SpatialCoordinate(msh)
    q_sep = Q * (1.0 - ufl.cos(ufl.pi * xs[0] / A_BOX))
    ss = solve_steady_conduction(
        msh, ft, k=K_GEL, q=q_sep, h=H_GEL, t_inf=T_INF, robin_tag=ROBIN_TAG, degree=DEGREE
    )
    clk.lap("UFL-source steady solve")
    out["ndof_s"] = int(ss.temperature.function_space.dofmap.index_map.size_global)
    out["imag"] = max(out["imag"], _max_imag(ss.temperature))
    th_s = ss.temperature - T_INF
    ex_s, npts_s = _exact_p(msh, _sep_ex())
    clk.lap(f"separable series interpolation ({npts_s} P3 points, {SEP_TERMS})")
    out["s6_int"] = _rel_p(msh, th_s, ex_s, chi)
    out["s6_whole"] = _rel_p(msh, th_s, ex_s)
    _say(
        f"MAT21-7 (S6) UFL-source steady interior {out['s6_int']:.6e} (step 6 "
        f"{STEP6_S1['interior']}); whole {out['s6_whole']:.6e} (step 6 {STEP6_S1['whole']})"
    )
    clk.lap("UFL-source readings")

    # (iii) The DG0 hand-off: cell average of cos(pi x / A) over every cell.
    tdim = msh.topology.dim
    n_loc = msh.topology.index_map(tdim).size_local
    ct = dmesh.meshtags(msh, tdim, np.arange(n_loc, dtype=np.int32), np.ones(n_loc, dtype=np.int32))
    W, cbar = _cell_average_dg0_tag(msh, ct, 1, ufl.cos(ufl.pi * xs[0] / A_BOX))
    q_h = fem.Function(W)
    q_h.x.array[:] = Q * (1.0 - cbar)
    q_h.x.scatter_forward()
    # Owned-cell volumes (DG0 array carries ghosts; sum owned only, then allreduce).
    w = ufl.TestFunction(W)
    one = fem.Constant(msh, np.asarray(1.0, dtype=q_h.x.array.dtype))
    bvol = fem.assemble_vector(fem.form(ufl.inner(one, w) * ufl.dx))
    bvol.scatter_reverse(la.InsertMode.add)
    n_own = W.dofmap.index_map.size_local
    vol = np.real(bvol.array[:n_own])
    cb = np.real(cbar[:n_own])
    a_loc = float(np.sum(vol * cb * cb))
    a_h = comm.allreduce(a_loc, op=MPI.SUM) / (V_BOX / 2.0)
    vol_tot = comm.allreduce(float(np.sum(vol)), op=MPI.SUM)
    out["a_h"] = a_h
    _say(
        f"MAT21-7 a_h = {a_h:.9f}; 1 - a_h = {1.0 - a_h:.4e} (predicted 6.08e-4, range "
        f"3e-4 to 1.22e-3; printed, not asserted); sum owned vol_K / V - 1 = "
        f"{vol_tot / V_BOX - 1.0:+.3e}; owned cells this rank {n_own}"
    )
    clk.lap("DG0 cell average + a_h")

    # (iv) DG0-source steady solve.
    sd = solve_steady_conduction(
        msh, ft, k=K_GEL, q=q_h, h=H_GEL, t_inf=T_INF, robin_tag=ROBIN_TAG, degree=DEGREE
    )
    clk.lap("DG0-source steady solve")
    out["ndof_d"] = int(sd.temperature.function_space.dofmap.index_map.size_global)
    out["imag"] = max(out["imag"], _max_imag(sd.temperature))
    th_d = sd.temperature - T_INF
    conv = float(np.real(_allsum(comm, H_GEL * th_d * ds)))
    sp = float(np.real(sd.source_power))
    out["gauss_d"] = abs(conv - sp) / sp
    _say(
        f"MAT21-7 (S2-dg0) Gauss |h oint theta_h - source_power|/source_power = "
        f"{out['gauss_d']:.3e} (bound {S2_GAUSS_BAND:.0e}, predicted <= 5e-12); "
        f"source_power/(Q V_h) - 1 = {sp / (Q * v_h) - 1.0:+.3e}; V_h/V - 1 = "
        f"{v_h / V_BOX - 1.0:+.3e}"
    )

    # (v) Comparands and readings.
    ex_a, npts_a = _exact_p(msh, _sep_ex_a(a_h))
    clk.lap(f"a_h series interpolation ({npts_a} P3 points)")
    out["d_int"] = _rel_p(msh, th_d, ex_s, chi)
    d_whole = _rel_p(msh, th_d, ex_s)
    out["da_int"] = _rel_p(msh, th_d, ex_a, chi)
    out["c_dg0"] = _rel_p(msh, th_d, ex_u, chi)
    out["c_a"] = _rel_p(msh, th_s, ex_a, chi)
    pv = _point_vals(sd.temperature, [[0.0, 0.0, 0.0], [0.9 * A_BOX, 0.0, 0.0]])
    _say(
        f"MAT21-7 (S1-dg0) DG0 steady interior vs P32 series {out['d_int']:.6e} (bound "
        f"{INTERIOR_BAND:.0e}, predicted 2.9e-4 to 3.6e-4, range 1.2e-4 to 7e-4); whole "
        f"{d_whole:.6e}; theta_h(0) {pv[0]:.6f} K (series 2.741147), theta_h(0.9A,0,0) "
        f"{pv[1]:.6f} K (series 14.919189)"
    )
    _say(
        f"MAT21-7 (S1-dg0a) DG0 steady interior vs a_h series {out['da_int']:.6e} (bound "
        f"{INTERIOR_BAND:.0e}, predicted 3e-5 to 6e-5, range 1.5e-5 to 1.5e-4)"
    )
    _say(
        f"MAT21-7 control (c-dg0) DG0 run vs uniform series {out['c_dg0']:.4e} (predicted 0.512, "
        f"range 0.50-0.53); (c-a) UFL run vs a_h series {out['c_a']:.4e} (predicted approx "
        f"3.2e-4; asserted only (c-a) > (S1-dg0a))"
    )
    clk.lap("DG0 steady readings")

    # CN N = N_MAIN to t* with the DG0 source.
    res = solve_transient_conduction(
        msh,
        ft,
        k=K_GEL,
        rho_c=RHO_C,
        q=q_h,
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
    clk.lap(f"CN N={N_MAIN} (DG0 source)")
    out["imag"] = max(out["imag"], _max_imag(res.temperature))
    out["E"] = float(res.max_balance_residual)
    th_t = res.temperature - T_INF
    ex_t, _ = _exact_p(msh, _sep_trans_ex(T_STAR))
    ex_ta, _ = _exact_p(msh, _sep_trans_ex_a(a_h, T_STAR))
    out["t_int"] = _rel_p(msh, th_t, ex_t, chi)
    t_whole = _rel_p(msh, th_t, ex_t)
    out["ta_int"] = _rel_p(msh, th_t, ex_ta, chi)
    clk.lap("transient interpolations + readings")
    tc = _point_vals(res.temperature, np.zeros((1, 3)))[0]
    sp_t = getattr(res, "source_power", float("nan"))
    _say(
        f"MAT21-7 (T-dg0) CN N={N_MAIN} at t*: interior vs P32 series {out['t_int']:.6e} "
        f"(bound {INTERIOR_BAND:.0e}, predicted 3e-4 to 5e-4, range 1.5e-4 to 1e-3); whole "
        f"{t_whole:.6e}; vs a_h series {out['ta_int']:.6e} (printed); theta_h(0, t*) "
        f"{tc:.6f} K (series 0.024797)"
    )
    _say(
        f"MAT21-7 (E-dg0) CN max balance residual {out['E']:.3e} (bound {E_BAND:.0e}, predicted "
        f"<= 5e-13); transient source_power/(Q V_h) - 1 = "
        f"{float(np.real(sp_t)) / (Q * v_h) - 1.0:+.3e}"
    )

    fl = _sep_floor()
    _say(
        f"MAT21-7 comparand floor 1 - prod S' at SEP_TERMS {fl:.4e} (step 6 2.2196e-6); "
        f"max|Im theta| {out['imag']:.3e}; fixture elapsed {time.perf_counter() - t0:.1f} s"
    )
    return out


def test_precondition_constants():
    assert K == K_GEL and H == H_GEL and RHO_C_0 == RHO_C


@complex_only
def test_precondition_degree_dofs_m4(dg0_m4):
    pred = _dofs_pred(M4)
    _say(
        f"MAT21-7 (Pdeg) P2 dofs uniform {dg0_m4['ndof_u']} / UFL {dg0_m4['ndof_s']} / DG0 "
        f"{dg0_m4['ndof_d']} (predicted {pred} = 116025)"
    )
    assert pred == 116025
    assert dg0_m4["ndof_u"] == pred and dg0_m4["ndof_s"] == pred and dg0_m4["ndof_d"] == pred


@complex_only
def test_u_uniform_control_reproduces_step3(dg0_m4):
    assert f"{dg0_m4['u_whole']:.6e}" == STEP3_M4["whole"]
    assert f"{dg0_m4['u_int']:.6e}" == STEP3_M4["interior"]


@complex_only
def test_s6_ufl_source_reproduces_step6(dg0_m4):
    assert f"{dg0_m4['s6_int']:.6e}" == STEP6_S1["interior"]
    assert f"{dg0_m4['s6_whole']:.6e}" == STEP6_S1["whole"]


@complex_only
def test_s1_dg0_interior(dg0_m4):
    assert dg0_m4["d_int"] <= INTERIOR_BAND


@complex_only
def test_s1_dg0a_interior(dg0_m4):
    assert dg0_m4["da_int"] <= INTERIOR_BAND


@complex_only
def test_s2_dg0_gauss_balance(dg0_m4):
    assert dg0_m4["gauss_d"] <= S2_GAUSS_BAND


@complex_only
def test_t_dg0_transient_interior(dg0_m4):
    assert dg0_m4["t_int"] <= INTERIOR_BAND


@complex_only
def test_e_dg0_cn_balance(dg0_m4):
    assert dg0_m4["E"] <= E_BAND


@complex_only
def test_control_cdg0_uniform_series_misses(dg0_m4):
    assert dg0_m4["c_dg0"] > INTERIOR_BAND


@complex_only
def test_control_ca_ordering(dg0_m4):
    assert dg0_m4["c_a"] > dg0_m4["da_int"]
