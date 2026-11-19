"""`WF-15` step 2 (§9 item 260, chain U10): ``Study.temperature(RobinCooling)``.

On `MAT-4` step 1's sphere the Study's solved field drives `MAT-10`'s steady
Robin-cooled conduction.  Anchors (brief: docs/planning/chunks/WF-15.md,
"§9 item 260"):

(K)    spec defaults / refusals, no solve;
(Once) ``temperature`` after ``solve`` keeps ``study.result is res``;
(UT)   the Study builds `MAT-10`'s thermal problem bitwise upstream of the
       thermal solve (submesh / Robin facet counts, parent map, q_h);
(IT)   θ end to end against `MAT-10`'s own helper chain ≤ ``L0_BAND``;
(G)    `MAT-10`'s (A) / (P) / (B) / (C) on the Study's outputs, imported bands;
(R)    `MAT-10`'s recorded (A) / (P) digits at ``.4e`` (degree 1 only);
(Deg)  default EM degree 2, thermal space degree 1 (W2).
Controls: (c1) missing ½ via ``source_scale = 2·G`` misses ``A_BAND``
(asserted at degree 1, printed at degree 2); (c2) the (UT) comparator reads a
planted ``np.spacing`` at full resolution (asserted).

Every band, constant, helper and comparator is imported, none restated.
"""

from __future__ import annotations

import time

import numpy as np
import pytest
import ufl
from dolfinx import default_scalar_type, fem
from mpi4py import MPI

from fem_em_solver import (
    FourPortBirdcageGeometry,
    LumpedPortDrive,
    RobinCooling,
    Study,
)
from fem_em_solver.post.sar import uniform_sphere_sar_closed_form
from fem_em_solver.thermal.steady import solve_steady_conduction
from fem_em_solver.utils.analytical import robin_sphere_steady_temperature
from fem_em_solver.utils.constants import EPSILON_0
from tests.complex_mode import complex_only
from tests.materials.test_phantom_material_model import _dg0_owned_diff
from tests.package.test_wf15_study_lossy_sphere import (
    _material_057,
    _study,
    _summed_rss_gib,
)
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
    B_BAND,
    C_BAND,
    P_BAND,
    ROBIN_TAG,
    T_INF,
    G,
    H,
    K,
    _allsum,
    _cell_average_dg0,
    _interior_indicator,
    _interior_rel,
    _robin_tags,
    _say,
    _sphere_submesh,
    _theta_cf,
    _to_sub,
)

# `MAT-10`'s recorded digits (20260930T021802Z_MAT-10-step1-W1.log:158–167,
# reproduced 2026-10-03 at 20261003T185646Z_MAT-18-step2b-W3a-regate.log).
MAT10_SUB_CELLS = 19496
MAT10_A_PRINT = "3.2464e-02"
MAT10_P_PRINT = "3.1923e-03"
MAT10_B_MISSING_HALF_PRINT = "1.0649e+00"
MAT10_SURFACE_RATIO_RECORD = 1.033687


def _space_degree(V):
    d = V.ufl_element().degree
    return d() if callable(d) else d


def _owned_facet_count(ft, comm):
    topo = ft.topology if hasattr(ft, "topology") else ft.mesh.topology
    fimap = topo.index_map(ft.dim)
    n = int(np.count_nonzero(np.asarray(ft.indices) < fimap.size_local))
    return int(comm.allreduce(n, op=MPI.SUM))


def _theta_owned_rel(T_a, T_b, comm):
    """Owned max |T_a − T_b| / max |T_b − T∞|, reduced (P1, bs 1)."""
    n_a = T_a.function_space.dofmap.index_map.size_local
    n_b = T_b.function_space.dofmap.index_map.size_local
    same = bool(comm.allreduce(bool(n_a == n_b), op=MPI.LAND))
    assert same, "(IT) temperature fields do not share the owned layout"
    a = np.asarray(T_a.x.array[:n_a], dtype=np.complex128)
    b = np.asarray(T_b.x.array[:n_b], dtype=np.complex128)
    dmax = comm.allreduce(float(np.max(np.abs(a - b), initial=0.0)), op=MPI.MAX)
    nmax = comm.allreduce(float(np.max(np.abs(b - T_INF), initial=0.0)), op=MPI.MAX)
    return dmax / nmax, nmax


def _temperature_leg(study, res, comm):
    """Drive the Study's temperature and `MAT-10`'s chain on the same inputs."""
    rd = {}
    t0 = time.perf_counter()
    spec = RobinCooling(K, H, T_INF, source_scale=G)
    th = study.temperature(spec)
    rd["t_bridge_thermal"] = time.perf_counter() - t0
    rd["once"] = study.result is res

    # `MAT-10`'s own helper chain on the Study's E / σ / mesh / tags.
    msh, cell_tags = res.mesh, res.cell_tags
    E, sigma = res.fields.e_complex, res.fields.sigma_field
    sub_h, pid_h = _sphere_submesh(msh, cell_tags)
    ft_h = _robin_tags(sub_h)
    q_expr = G * 0.5 * sigma * ufl.inner(E, E)
    Wp, q_par = _cell_average_dg0(msh, cell_tags, q_expr)
    q_hh = _to_sub(Wp, q_par, sub_h, pid_h)

    # (UT)
    rd["n_sub_study"] = th.submesh.topology.index_map(3).size_global
    rd["n_sub_helper"] = sub_h.topology.index_map(3).size_global
    rd["n_fac_study"] = _owned_facet_count(th.facet_tags, comm)
    rd["n_fac_helper"] = _owned_facet_count(ft_h, comm)
    local_pc = bool(np.array_equal(np.asarray(th.parent_cells), np.asarray(pid_h)))
    local_ft = bool(
        np.array_equal(np.asarray(th.facet_tags.indices), np.asarray(ft_h.indices))
        and np.array_equal(np.asarray(th.facet_tags.values), np.asarray(ft_h.values))
    )
    rd["parent_equal"] = bool(comm.allreduce(local_pc, op=MPI.LAND))
    rd["facets_equal"] = bool(comm.allreduce(local_ft, op=MPI.LAND))
    rd["q_diff"] = _dg0_owned_diff(th.heat_source, q_hh, comm)
    _say(
        f"WF15T (UT) sub cells {rd['n_sub_study']} / {rd['n_sub_helper']}; Robin facets "
        f"{rd['n_fac_study']} / {rd['n_fac_helper']}; parent_cells equal {rd['parent_equal']}; "
        f"facet tags equal {rd['facets_equal']}; q_h max|d| {rd['q_diff'][0]!r}, "
        f"differing {rd['q_diff'][1]} / {rd['q_diff'][2]} owned"
    )

    # (c2) planted control of the (UT) comparator.
    planted = fem.Function(th.heat_source.function_space)
    planted.x.array[:] = th.heat_source.x.array
    n_own = planted.function_space.dofmap.index_map.size_local
    sp = 0.0
    if comm.rank == 0 and n_own > 0:
        v = planted.x.array[0]
        sp = float(np.spacing(np.real(v)))
        planted.x.array[0] = v + sp
    sp = comm.bcast(sp, root=0)
    rd["c2"] = _dg0_owned_diff(planted, th.heat_source, comm)
    rd["c2_spacing"] = sp
    _say(
        f"WF15T (c2) planted spacing {sp!r}: read max|d| {rd['c2'][0]!r}, differing "
        f"{rd['c2'][1]}"
    )

    # (IT) MAT-10's thermal call on the helper-built inputs.
    res_h = solve_steady_conduction(
        sub_h, ft_h, k=K, q=q_hh, h=H, t_inf=T_INF, robin_tag=ROBIN_TAG, degree=1
    )
    rd["it"], th_max = _theta_owned_rel(th.temperature, res_h.temperature, comm)
    _say(
        f"WF15T (IT) owned max|dT| / max|T - Tinf| = {rd['it']:.3e}  (band {L0_BAND:.0e}); "
        f"theta_h max = {th_max:.6f} K"
    )

    # (G) MAT-10's gates on the Study's outputs.
    cf = uniform_sphere_sar_closed_form(
        e0=E0,
        epsilon_r=EPSILON_R_SPHERE,
        sigma=SIGMA_HIGH,
        omega=OMEGA,
        rho=RHO_KG_M3,
        epsilon_0=EPSILON_0,
    )
    q_cf = G * RHO_KG_M3 * cf["sar_w_per_kg"]
    sub, ft = th.submesh, th.facet_tags
    chi = _interior_indicator(sub)
    theta = _theta_cf(sub, q_cf)
    rd["err_A"] = _interior_rel(sub, chi, (th.temperature - T_INF) - theta, theta)
    res_P = solve_steady_conduction(
        sub, ft, k=K, q=q_cf, h=H, t_inf=T_INF, robin_tag=ROBIN_TAG, degree=1
    )
    rd["err_P"] = _interior_rel(sub, chi, (res_P.temperature - T_INF) - theta, theta)
    int_q = float(np.real(_allsum(comm, th.heat_source * ufl.dx(domain=sub))))
    p_parent = spec.source_scale * study.mean_sar()["dissipated_power_w"]
    rd["err_B"] = abs(int_q - p_parent) / abs(p_parent)
    rd["err_C"] = abs(th.convected_power - th.source_power) / abs(th.source_power)
    rd["once_after_sar"] = study.result is res
    ds = ufl.Measure("ds", domain=sub, subdomain_data=ft)(ROBIN_TAG)
    area = float(np.real(_allsum(comm, fem.Constant(sub, default_scalar_type(1.0)) * ds)))
    th_s = float(np.real(_allsum(comm, (th.temperature - T_INF) * ds))) / area
    R = study.geometry.sphere_radius
    thR = robin_sphere_steady_temperature(R, radius=R, k=K, q=q_cf, h=H, t_inf=T_INF) - T_INF
    rd["surface_ratio"] = th_s / thR
    rd["thermal_degree"] = _space_degree(th.temperature.function_space)
    _say(
        f"WF15T (G) (A) = {rd['err_A']:.10e} [{rd['err_A']:.4e}] (band {A_BAND}); "
        f"(P) = {rd['err_P']:.10e} [{rd['err_P']:.4e}] (band {P_BAND:.0e}); "
        f"(B) = {rd['err_B']:.3e} (band {B_BAND:.0e}); (C) = {rd['err_C']:.3e} "
        f"(band {C_BAND:.0e})"
    )
    _say(
        f"WF15T surface-mean theta_h / theta(R) = {rd['surface_ratio']:.6f} "
        f"(degree-1 record {MAT10_SURFACE_RATIO_RECORD}); thermal degree "
        f"{rd['thermal_degree']}; q_cf {q_cf:.6f} W/m^3"
    )

    # (c1) missing ½ through the public knob.
    th2 = study.temperature(RobinCooling(K, H, T_INF, source_scale=2.0 * G))
    sub2 = th2.submesh  # each temperature() call builds its own submesh
    theta2 = _theta_cf(sub2, q_cf)
    rd["c1"] = _interior_rel(
        sub2, _interior_indicator(sub2), (th2.temperature - T_INF) - theta2, theta2
    )
    _say(
        f"WF15T (c1) source_scale = 2G: (A) comparator = {rd['c1']:.4e} "
        f"(MAT-10 record {MAT10_B_MISSING_HALF_PRINT}; band {A_BAND})"
    )
    _say(
        f"WF15T bridge + thermal (Study) {rd['t_bridge_thermal']:.1f} s; summed "
        f"ru_maxrss {_summed_rss_gib(comm):.2f} GiB"
    )
    return rd


def _assert_common(rd):
    assert rd["once"] and rd["once_after_sar"], "(Once) a second EM solve happened"
    assert rd["n_sub_study"] == rd["n_sub_helper"], "(UT) submesh cell counts differ"
    assert rd["n_fac_study"] == rd["n_fac_helper"], "(UT) Robin facet counts differ"
    assert rd["parent_equal"], "(UT) parent_cells differ"
    assert rd["facets_equal"], "(UT) Robin facet tags differ"
    assert rd["q_diff"][0] == 0.0 and rd["q_diff"][1] == 0, f"(UT) q_h {rd['q_diff']}"
    assert rd["c2_spacing"] > 0.0
    assert rd["c2"][0] == rd["c2_spacing"] and rd["c2"][1] == 1, f"(c2) {rd['c2']}"
    assert rd["it"] <= L0_BAND, f"(IT) {rd['it']:.3e} > {L0_BAND:.0e}"
    assert rd["err_A"] <= A_BAND, f"(A) {rd['err_A']:.4e} > {A_BAND}"
    assert rd["err_P"] <= P_BAND, f"(P) {rd['err_P']:.4e} > {P_BAND}"
    assert rd["err_B"] <= B_BAND, f"(B) {rd['err_B']:.3e} > {B_BAND:.0e}"
    assert rd["err_C"] <= C_BAND, f"(C) {rd['err_C']:.3e} > {C_BAND:.0e}"
    assert rd["thermal_degree"] == 1, f"thermal degree {rd['thermal_degree']}"


def test_wf15_temperature_spec_and_refusals():
    """(K): defaults, ValueError on non-positive k / h / scale, TypeError on a port Study."""
    spec = RobinCooling(K, H, T_INF)
    assert spec.degree == 1
    assert spec.source_scale == 1.0
    for bad in (
        dict(k_w_per_m_k=0.0, h_w_per_m2_k=H, t_inf_c=T_INF),
        dict(k_w_per_m_k=-K, h_w_per_m2_k=H, t_inf_c=T_INF),
        dict(k_w_per_m_k=K, h_w_per_m2_k=0.0, t_inf_c=T_INF),
        dict(k_w_per_m_k=K, h_w_per_m2_k=-H, t_inf_c=T_INF),
        dict(k_w_per_m_k=K, h_w_per_m2_k=H, t_inf_c=T_INF, source_scale=0.0),
        dict(k_w_per_m_k=K, h_w_per_m2_k=H, t_inf_c=T_INF, source_scale=-1.0),
    ):
        with pytest.raises(ValueError):
            RobinCooling(**bad)
    port = Study(
        FourPortBirdcageGeometry.f_small(), _material_057(), LumpedPortDrive(), frequency_hz=64e6
    )
    with pytest.raises(TypeError):
        port.temperature(spec)
    assert port.result is None  # no mesh was built
    _say(
        "WF15T (K) defaults degree 1 / scale 1.0; 6 ValueErrors; port Study TypeError, "
        "result None"
    )


@complex_only
def test_wf15_temperature_degree1_reproduces_mat10_record():
    """W1: (Once), (UT), (IT), (G), (R) and (c1) / (c2) at EM degree 1."""
    comm = MPI.COMM_WORLD
    study = _study(degree=1)
    t0 = time.perf_counter()
    res = study.solve(comm)
    t_em = time.perf_counter() - t0
    _say(f"WF15T degree-1 EM solve + mesh {t_em:.1f} s")
    rd = _temperature_leg(study, res, comm)
    _say(
        f"WF15T (R) (A) {rd['err_A']:.4e} vs {MAT10_A_PRINT}; (P) {rd['err_P']:.4e} vs "
        f"{MAT10_P_PRINT}"
    )
    _assert_common(rd)
    assert rd["n_sub_study"] == MAT10_SUB_CELLS, f"(UT) sub cells {rd['n_sub_study']}"
    assert f"{rd['err_A']:.4e}" == MAT10_A_PRINT, f"(R) (A) {rd['err_A']:.4e}"
    assert f"{rd['err_P']:.4e}" == MAT10_P_PRINT, f"(R) (P) {rd['err_P']:.4e}"
    assert rd["c1"] > A_BAND, f"(c1) {rd['c1']:.4e} did not miss {A_BAND}"


@complex_only
def test_wf15_temperature_default_degree2_leg():
    """W2: (Deg), (Once), (UT), (IT), (G) at the default EM degree 2; (c1) printed."""
    comm = MPI.COMM_WORLD
    study = _study()  # the default degree
    t0 = time.perf_counter()
    res = study.solve(comm)
    t_em = time.perf_counter() - t0
    em_degree = _space_degree(res.fields.e_complex.function_space)
    _say(
        f"WF15T degree-2 EM solve + mesh {t_em:.1f} s; study.degree {study.degree}; "
        f"EM space degree {em_degree}; summed ru_maxrss {_summed_rss_gib(comm):.2f} GiB"
    )
    rd = _temperature_leg(study, res, comm)
    _say(
        f"WF15T degree-2 (A) record = {rd['err_A']:.4e} (predicted [0.1, 1.0] %); "
        f"(c1) {rd['c1']:.4e} (predicted ≈ 1.00, printed only)"
    )
    assert study.degree == 2 and res.degree == 2, "(Deg) study degree"
    assert em_degree == 2, f"(Deg) EM space degree {em_degree}"
    _assert_common(rd)
