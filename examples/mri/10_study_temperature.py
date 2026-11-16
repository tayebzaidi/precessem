"""Example (`EX-94`): a temperature from one ``Study`` -- steady conduction on the
SAR-heated gel sphere.

`WF-15`'s ``Study`` solves the lossy sphere (`MAT-4` step 1: imposed 1 V/m field,
64 MHz) and ``study.temperature(RobinCooling(...))`` hands the solved field to
`MAT-10`'s steady Robin-cooled conduction: the DG0 heat source q = G * sigma |E|^2 / 2
on the sphere submesh, conduction k, convection h to T_inf.  This example makes the
two calls at **degree 1** (the width and degree of `MAT-10`'s recorded digits) and
reads `MAT-10`'s (A) / (P) / (B) / (C) off the returned ``ThermalStudyResult``.

Every constant, band and helper is imported from the gate modules (the ANS-1 rule);
k, h, T_inf, G, sigma, rho and R are read from imports or from ``study.geometry``.
The only literals are the control's factor 10 and the 2 of (c1).

Asserted:

* (R)    (A) and (P) at ``.4e`` equal `MAT-10`'s recorded prints;
* (G)    (A) <= ``A_BAND``, (P) <= ``P_BAND``, (B) <= ``B_BAND``, (C) <= ``C_BAND``;
* (Once) ``study.result is res`` after ``temperature()`` and after ``mean_sar()``;
* (Sub)  the sphere submesh has `MAT-10`'s recorded global cell count.

Negative controls:

* (c-h)  asserted: ``h -> 10 h`` -- the (A) comparator (closed form with H baked in),
  on that call's own submesh, exceeds ``A_BAND``;
* (c-h') printed, not asserted: the FEM surface-rise ratio theta_s(10 H) / theta_s(H)
  against the closed form's factor, and their relative difference;
* (c1)   asserted: ``source_scale = 2 G`` (the missing 1/2) -- the (A) comparator
  equals `MAT-10`'s recorded print and exceeds ``A_BAND``.

Run it through the example runner (the ``mri:`` group sources the complex build)::

    ./scripts/run_examples.sh -e mri:10 -n 2 -t 180

Not claimed: transient heating, perfusion, tissue or a coil; degree-2 thermal
accuracy; a gate on the h-scaling ratio; anything `WF-15` / `MAT-10` do not claim.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import ufl
from mpi4py import MPI

from dolfinx import default_scalar_type, fem

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from fem_em_solver import RobinCooling  # noqa: E402
from fem_em_solver.io.paraview_utils import (  # noqa: E402
    adopt_host_ownership,
    write_xdmf_with_tags,
)
from fem_em_solver.post.setup_figure import write_setup_figure  # noqa: E402
from fem_em_solver.post.sar import uniform_sphere_sar_closed_form  # noqa: E402
from fem_em_solver.thermal.steady import solve_steady_conduction  # noqa: E402
from fem_em_solver.utils.analytical import robin_sphere_steady_temperature  # noqa: E402
from fem_em_solver.utils.constants import EPSILON_0  # noqa: E402

from tests.package.test_wf15_study_lossy_sphere import _study  # noqa: E402
from tests.package.test_wf15_study_temperature import (  # noqa: E402
    MAT10_A_PRINT,
    MAT10_B_MISSING_HALF_PRINT,
    MAT10_P_PRINT,
    MAT10_SUB_CELLS,
    MAT10_SURFACE_RATIO_RECORD,
)
from tests.validation.test_lossy_sphere_sar import (  # noqa: E402
    E0,
    EPSILON_R_SPHERE,
    OMEGA,
    RHO_KG_M3,
    SIGMA_HIGH,
)
from tests.validation.test_mat10_sar_to_heat_sphere import (  # noqa: E402
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

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"
BASENAME = "mri_10_study_temperature"
DEGREE = 1  # MAT-10's record degree; the Study default is 2.
AIR_TAG = 2
H_FACTOR = 10.0  # the (c-h) control's only literal


def _a_comparator(th, q_cf):
    """(A) on this call's own submesh: interior rel L2 of (T - T_inf) vs the closed form."""
    sub = th.submesh
    chi = _interior_indicator(sub)
    theta = _theta_cf(sub, q_cf)
    return _interior_rel(sub, chi, (th.temperature - T_INF) - theta, theta)


def _surface_mean(th, comm):
    """Mean of (T - T_inf) over the Robin surface of this call's submesh (collective)."""
    sub = th.submesh
    ds = ufl.Measure("ds", domain=sub, subdomain_data=th.facet_tags)(ROBIN_TAG)
    area = float(np.real(_allsum(comm, fem.Constant(sub, default_scalar_type(1.0)) * ds)))
    return float(np.real(_allsum(comm, (th.temperature - T_INF) * ds))) / area


def main() -> None:
    comm = MPI.COMM_WORLD
    started = time.perf_counter()

    # ---- the Study: two calls, one EM solve --------------------------------
    study = _study(degree=DEGREE)
    t0 = time.perf_counter()
    res = study.solve(comm)
    t_em = time.perf_counter() - t0

    spec = RobinCooling(K, H, T_INF, source_scale=G)
    t0 = time.perf_counter()
    th = study.temperature(spec)
    t_th = time.perf_counter() - t0
    once_after_temperature = study.result is res

    # EX-57 setup figure; no-op unless FEM_EM_SETUP_FIGURES=1.
    geometry = study.geometry
    write_setup_figure(
        res.mesh,
        res.cell_tags,
        FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names={
            geometry.sphere_tag: "gel sphere (Robin-cooled surface)",
            AIR_TAG: "air box (imposed-field walls)",
        },
        hide_tags=(AIR_TAG,),
        translucent_tags=(geometry.sphere_tag,),
        slice_normal=(0.0, 0.0, 1.0),
        slice_origin=(0.0, 0.0, 0.0),
        title="mri:10 -- Study temperature on the SAR-heated sphere (WF-15, degree 1)",
        comm=comm,
    )

    # ---- (A) / (P) / (B) / (C), term for term as WF-15's _temperature_leg ---
    cf = uniform_sphere_sar_closed_form(
        e0=E0, epsilon_r=EPSILON_R_SPHERE, sigma=SIGMA_HIGH, omega=OMEGA,
        rho=RHO_KG_M3, epsilon_0=EPSILON_0,
    )
    q_cf = G * RHO_KG_M3 * cf["sar_w_per_kg"]
    sub, ft = th.submesh, th.facet_tags
    chi = _interior_indicator(sub)
    theta = _theta_cf(sub, q_cf)
    err_a = _interior_rel(sub, chi, (th.temperature - T_INF) - theta, theta)
    res_p = solve_steady_conduction(
        sub, ft, k=K, q=q_cf, h=H, t_inf=T_INF, robin_tag=ROBIN_TAG, degree=1
    )
    err_p = _interior_rel(sub, chi, (res_p.temperature - T_INF) - theta, theta)
    int_q = float(np.real(_allsum(comm, th.heat_source * ufl.dx(domain=sub))))
    p_parent = spec.source_scale * study.mean_sar()["dissipated_power_w"]
    err_b = abs(int_q - p_parent) / abs(p_parent)
    err_c = abs(th.convected_power - th.source_power) / abs(th.source_power)
    once_after_sar = study.result is res
    n_sub = int(sub.topology.index_map(sub.topology.dim).size_global)

    radius = geometry.sphere_radius
    th_s = _surface_mean(th, comm)
    th_r = (
        robin_sphere_steady_temperature(radius, radius=radius, k=K, q=q_cf, h=H, t_inf=T_INF)
        - T_INF
    )
    surface_ratio = th_s / th_r

    # ---- (c-h): h -> 10 h ---------------------------------------------------
    t0 = time.perf_counter()
    th_h = study.temperature(RobinCooling(K, H_FACTOR * H, T_INF, source_scale=G))
    t_th_h = time.perf_counter() - t0
    c_h = _a_comparator(th_h, q_cf)
    th_s_h = _surface_mean(th_h, comm)
    fem_factor = th_s_h / th_s
    cf_factor = (
        robin_sphere_steady_temperature(
            radius, radius=radius, k=K, q=q_cf, h=H_FACTOR * H, t_inf=T_INF
        )
        - T_INF
    ) / th_r
    factor_rel = abs(fem_factor - cf_factor) / abs(cf_factor)

    # ---- (c1): the missing half ---------------------------------------------
    th_2 = study.temperature(RobinCooling(K, H, T_INF, source_scale=2.0 * G))
    c1 = _a_comparator(th_2, q_cf)
    once_final = study.result is res

    # ---- output fields on the sphere submesh (every collective on every rank)
    theta_out = fem.Function(th.temperature.function_space, name="theta")
    theta_out.x.array[:] = np.real(th.temperature.x.array) - T_INF
    theta_out.x.scatter_forward()
    q_out = fem.Function(th.heat_source.function_space, name="q")
    q_out.x.array[:] = np.real(th.heat_source.x.array)
    q_out.x.scatter_forward()
    sar_out = fem.Function(th.heat_source.function_space, name="SAR")
    sar_out.x.array[:] = np.real(th.heat_source.x.array) / (G * RHO_KG_M3)
    sar_out.x.scatter_forward()
    OUTPUT_DIR.mkdir(exist_ok=True)
    xdmf_path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_combined",
        sub,
        None,
        {"theta": theta_out, "SAR": sar_out, "q": q_out},
        comm=comm,
    )
    adopt_host_ownership(OUTPUT_DIR, comm=comm)

    # ---- prints (rank 0, no collective inside) -----------------------------
    if comm.rank == 0:
        print("=" * 72)
        print("EX-94 -- a temperature from one Study: steady conduction on the sphere (WF-15)")
        print("=" * 72)
        print(
            f"[study] degree {study.degree}, -n {comm.size}, k {K:g}, h {H:g}, T_inf {T_INF:g}, "
            f"source_scale G {G:g}, R {radius:g} m; EM solve + mesh {t_em:.1f} s; "
            f"temperature() {t_th:.2f} s"
        )
        print(f"[(Sub)] sphere submesh {n_sub} cells (record {MAT10_SUB_CELLS})")
        print(f"[(Once)] after temperature {once_after_temperature}; after mean_sar {once_after_sar}")
        print(
            f"[(R)] (A) {err_a:.4e} vs {MAT10_A_PRINT}; (P) {err_p:.4e} vs {MAT10_P_PRINT}"
        )
        print(
            f"[(G)] (A) = {err_a:.10e} (band {A_BAND}); (P) = {err_p:.10e} (band {P_BAND:.0e}); "
            f"(B) = {err_b:.3e} (band {B_BAND:.0e}); (C) = {err_c:.3e} (band {C_BAND:.0e})"
        )
        print(
            f"[surface] surface-mean theta_h / theta(R) = {surface_ratio:.6f} "
            f"(degree-1 record {MAT10_SURFACE_RATIO_RECORD}; printed, not asserted)"
        )
        print(
            f"[(c-h)] h -> {H_FACTOR:g} h ({t_th_h:.2f} s): (A) comparator = {c_h:.4e} "
            f"(must exceed A_BAND {A_BAND}; predicted [0.55, 0.70], stop if < 0.3)"
        )
        print(
            f"[(c-h')] FEM surface-rise ratio {fem_factor:.6f}; closed-form factor "
            f"{cf_factor:.6f}; rel diff {factor_rel:.3e} (predicted <= 1e-2; printed, not asserted)"
        )
        print(
            f"[(c1)] source_scale = 2G: (A) comparator = {c1:.4e} (MAT-10 record "
            f"{MAT10_B_MISSING_HALF_PRINT}; must exceed A_BAND {A_BAND})"
        )
        print(f"[paraview] wrote {xdmf_path}")

    assert once_after_temperature and once_after_sar and once_final, (
        "(Once) a second EM solve happened"
    )
    assert n_sub == MAT10_SUB_CELLS, f"(Sub) sub cells {n_sub}"
    assert f"{err_a:.4e}" == MAT10_A_PRINT, f"(R) (A) {err_a:.4e}"
    assert f"{err_p:.4e}" == MAT10_P_PRINT, f"(R) (P) {err_p:.4e}"
    assert err_a <= A_BAND, f"(G) (A) {err_a:.4e} > {A_BAND}"
    assert err_p <= P_BAND, f"(G) (P) {err_p:.4e} > {P_BAND}"
    assert err_b <= B_BAND, f"(G) (B) {err_b:.3e} > {B_BAND:.0e}"
    assert err_c <= C_BAND, f"(G) (C) {err_c:.3e} > {C_BAND:.0e}"
    assert c_h > A_BAND, f"(c-h) blind: {c_h:.4e} <= {A_BAND}"
    assert f"{c1:.4e}" == MAT10_B_MISSING_HALF_PRINT, f"(c1) {c1:.4e}"
    assert c1 > A_BAND, f"(c1) blind: {c1:.4e} <= {A_BAND}"

    if comm.rank == 0:
        print(
            f"\nAll anchors and controls hold. Total elapsed "
            f"{time.perf_counter() - started:.1f} s.",
            flush=True,
        )


if __name__ == "__main__":
    main()
