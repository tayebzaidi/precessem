"""Example (`EX-187`): the SAR-heated gel sphere over time, through one ``Study``.

`WF-15` step 15 (U131) drives `MAT-12`'s theta-scheme transient of the Robin-cooled gel
ball from ``Study.temperature(cool, transient=TransientSpec(...))``.  This example
writes that route as a user would, at EM **degree 1** (the degree of `MAT-12`'s
records): one EM solve of `MAT-4` step 1's imposed-field sphere, then

* A: ``TransientSpec(T_STAR / N_MAIN, N_MAIN, RHO_C)`` -- CN, N = 40 to t*;
* S: ``TransientSpec(DT_1, 1, RHO_C, theta=1.0, startup_half_steps=0)`` -- one BE step;
* Q: ``study.temperature(cool)`` -- the steady solve, for the control and the XDMF.

T_inf = T(0) = 0 as `MAT-12` runs it (theta = T).  Every constant, band, helper and
printed record is imported by name from the gate modules (the ANS-1 rule); there is no
gate-module edit and no ``src/`` edit.

Asserted (all collected first):

* (Once)  ``study.result is res`` after the three calls;  (UT) sub cells = `MAT10_SUB_CELLS`;
* (I-dir) owned max |theta_Study - theta_direct| / max |theta_direct| <= ``L0_BAND`` on A and S,
          theta_direct = `MAT-12`'s own ``_cn`` / ``_be_one`` on the Study's submesh;
* (A_t)   interior rel L2 against the series <= ``A_BAND``; its ``.4e`` = `MAT12_A_T_PRINT`;
* (S2)    calorimetric ratio on S, |ratio - 1| <= ``S2_BAND``; its ``.4e`` = `MAT12_S2_PRINT`;
* (E)     per-step balance <= ``E_BAND`` on A and S;
* (c-qs)  asserted control: the steady theta standing in for theta(t*) misses ``A_BAND``
          (`MAT-12` control (b) record 1.0108e+00);
* max |Im theta| = 0 on A (the group runs complex).

Printed only: theta at the sphere centre, FEM and closed form, at t* and at steady state,
and theta(t*)/theta(inf).

Run it through the example runner (the ``mri:`` group sources the complex build)::

    ./scripts/run_examples.sh -e mri:40 -n 2 -t 180

Not claimed: EM degree 2, duty cycle, perfusion, probes or a YAML key on the sphere, a
coil, tissue, a compliance figure.
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

from fem_em_solver import RobinCooling, TransientSpec  # noqa: E402
from fem_em_solver.io.paraview_utils import (  # noqa: E402
    adopt_host_ownership,
    write_xdmf_with_tags,
)
from fem_em_solver.post.setup_figure import write_setup_figure  # noqa: E402
from fem_em_solver.post.sar import uniform_sphere_sar_closed_form  # noqa: E402
from fem_em_solver.utils.analytical import robin_sphere_steady_temperature  # noqa: E402
from fem_em_solver.utils.constants import EPSILON_0  # noqa: E402

from tests.package.test_wf15_study_lossy_sphere import _study  # noqa: E402
from tests.package.test_wf15_study_sphere_transient import (  # noqa: E402
    MAT12_A_T_PRINT,
    MAT12_S2_PRINT,
    _owned_rel,
)
from tests.package.test_wf15_study_temperature import MAT10_SUB_CELLS  # noqa: E402
from tests.validation.test_lossy_sphere_sar import (  # noqa: E402
    E0,
    EPSILON_R_SPHERE,
    OMEGA,
    RHO_KG_M3,
    SIGMA_HIGH,
)
from tests.validation.test_mat10_sar_to_heat_sphere import (  # noqa: E402
    A_BAND,
    G,
    H,
    K,
    R,
    _allsum,
    _interior_indicator,
    _interior_rel,
    _say,
)
from tests.validation.test_mat12_transient_sar_heating_sphere import (  # noqa: E402
    DT_1,
    E_BAND,
    N_MAIN,
    RHO_C,
    S2_BAND,
    T_STAR,
    _be_one,
    _cn,
    _exact_theta,
    _mean,
    _series,
)
from tests.validation.test_mat21_wall_fem import _centre  # noqa: E402
from tests.validation.test_mat5_em_coupled_sphere import L0_BAND  # noqa: E402

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"
BASENAME = "mri_40_study_sar_driven_transient_sphere"
DEGREE = 1  # MAT-12's record degree; the Study default is 2.
AIR_TAG = 2
T_INF = 0.0  # MAT-12's T_inf = T(0): theta = T


def _max_imag(T):
    n_own = T.function_space.dofmap.index_map.size_local
    loc = float(np.max(np.abs(np.imag(T.x.array[:n_own])))) if n_own else 0.0
    return T.function_space.mesh.comm.allreduce(loc, op=MPI.MAX)


def _same_layout(Va, Vb, comm):
    """True on every rank iff the two spaces have identical owned dof coordinates."""
    na = Va.dofmap.index_map.size_local
    nb = Vb.dofmap.index_map.size_local
    ok = bool(na == nb)
    if ok and na:
        ok = bool(np.allclose(Va.tabulate_dof_coordinates()[:na],
                              Vb.tabulate_dof_coordinates()[:nb], rtol=0.0, atol=1e-12))
    return bool(comm.allreduce(ok, op=MPI.LAND))


def main() -> None:
    comm = MPI.COMM_WORLD
    started = time.perf_counter()

    # ---- the Study: one EM solve ---------------------------------------------
    study = _study(degree=DEGREE)
    t0 = time.perf_counter()
    res = study.solve(comm)
    t_em = time.perf_counter() - t0
    cool = RobinCooling(K, H, T_INF, source_scale=G)

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
        title="mri:40 -- SAR-driven transient of the gel sphere through one Study (degree 1)",
        comm=comm,
    )

    # ---- three temperature calls, no second EM solve --------------------------
    t1 = time.perf_counter()
    thA = study.temperature(cool, transient=TransientSpec(T_STAR / N_MAIN, N_MAIN, RHO_C))
    thS = study.temperature(cool, transient=TransientSpec(DT_1, 1, RHO_C, theta=1.0,
                                                          startup_half_steps=0))
    t_th = time.perf_counter() - t1
    dA = _cn(thA.mesh, thA.facet_tags, thA.q)
    dS = _be_one(thS.mesh, thS.facet_tags, thS.q)
    thQ = study.temperature(cool)
    once = study.result is res

    # ---- gates ----------------------------------------------------------------
    n_sub = int(thA.mesh.topology.index_map(3).size_global)
    idA, layA = _owned_rel(thA.temperature, dA.temperature, comm)
    idS, layS = _owned_rel(thS.temperature, dS.temperature, comm)

    cf = uniform_sphere_sar_closed_form(
        e0=E0, epsilon_r=EPSILON_R_SPHERE, sigma=SIGMA_HIGH, omega=OMEGA,
        rho=RHO_KG_M3, epsilon_0=EPSILON_0,
    )
    q_cf = G * RHO_KG_M3 * cf["sar_w_per_kg"]
    sub = thA.mesh
    chi = _interior_indicator(sub)
    theta_cf = _exact_theta(sub, T_STAR, q_cf)
    err_A = _interior_rel(sub, chi, thA.temperature - theta_cf, theta_cf)

    subS = thS.mesh
    one = fem.Constant(subS, default_scalar_type(1.0))
    vol = float(np.real(_allsum(comm, one * ufl.dx(domain=subS))))
    qbar = _mean(subS, thS.q, vol)
    thbar = _mean(subS, thS.temperature, vol)
    ratio = RHO_C * thbar / (DT_1 * qbar)
    s2 = abs(ratio - 1.0)

    eA, eS = thA.max_balance_residual, thS.max_balance_residual

    subQ = thQ.submesh
    thetaQ = _exact_theta(subQ, T_STAR, q_cf)
    err_qs = _interior_rel(subQ, _interior_indicator(subQ), thQ.temperature - thetaQ, thetaQ)
    imag = _max_imag(thA.temperature)

    # ---- printed only: centre values ------------------------------------------
    c_A = _centre(thA.temperature)
    c_Q = _centre(thQ.temperature)
    c_cf_t = float(np.real(_series(np.zeros(1), T_STAR, q_cf)[0]))
    c_cf_inf = float(
        robin_sphere_steady_temperature(0.0, radius=R, k=K, q=q_cf, h=H, t_inf=T_INF)
    ) - T_INF

    # ---- output fields on the sphere submesh -----------------------------------
    V = thA.temperature.function_space
    theta_A = fem.Function(V, name="theta_Study_tstar")
    theta_A.x.array[:] = np.real(thA.temperature.x.array) - T_INF
    theta_A.x.scatter_forward()
    theta_series = fem.Function(V, name="theta_closed_form_tstar")
    theta_series.interpolate(theta_cf)
    theta_series.x.scatter_forward()
    steady = fem.Function(V, name="theta_steady")
    same = _same_layout(V, thQ.temperature.function_space, comm)
    if same:
        steady.x.array[:] = np.real(thQ.temperature.x.array) - T_INF
    else:  # layouts differ: write the closed-form steady theta instead (labelled in the print)
        steady.interpolate(lambda x: q_cf * R / (3.0 * H)
                           + q_cf * (R**2 - np.sum(x**2, axis=0)) / (6.0 * K))
    steady.x.scatter_forward()
    OUTPUT_DIR.mkdir(exist_ok=True)
    xdmf_path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_combined",
        sub,
        None,
        {"theta_Study_tstar": theta_A, "theta_closed_form_tstar": theta_series,
         "theta_steady": steady},
        comm=comm,
    )
    adopt_host_ownership(OUTPUT_DIR, comm=comm)

    # ---- prints -----------------------------------------------------------------
    _say("=" * 72)
    _say("EX-187 -- SAR-driven transient of the gel sphere through one Study (WF-15 step 15)")
    _say("=" * 72)
    _say(f"[study] degree {study.degree}, -n {comm.size}, k {K:g}, h {H:g}, T_inf {T_INF:g}, "
         f"G {G:g}, R {R:g} m; t* {T_STAR:.4f} s, A: dt {thA.transient.dt_s:.6f} s x "
         f"{thA.transient.n_steps} (theta {thA.transient.theta}, startup "
         f"{thA.transient.startup_half_steps}); S: dt {thS.transient.dt_s:.6e} s x "
         f"{thS.transient.n_steps}; q_cf {q_cf:.6f} W/m^3")
    _say(f"[(UT)] sub cells {n_sub} (MAT10_SUB_CELLS {MAT10_SUB_CELLS}); (Once) {once}")
    _say(f"[(I-dir)] A {idA:.3e}  S {idS:.3e}  (band {L0_BAND:.0e}); layouts {layA} / {layS}")
    _say(f"[(A_t)] interior rel L2 = {err_A:.10e} [{err_A:.4e}] (band {A_BAND}; "
         f"record {MAT12_A_T_PRINT})")
    _say(f"[(S2)] rho c mean(theta^1)/(dt_1 mean q_h) = {ratio:.9f}  |ratio - 1| = {s2:.4e} "
         f"(band {S2_BAND:.0e}; record {MAT12_S2_PRINT})")
    _say(f"[(E)] max balance A {eA:.3e}  S {eS:.3e}  (band {E_BAND:.0e})")
    _say(f"[(c-qs)] steady theta for theta(t*): interior rel L2 = {err_qs:.4e} "
         f"(must miss {A_BAND}; MAT-12 control (b) 1.0108e+00)")
    _say(f"[max|Im theta|] = {imag:.3e}")
    _say(f"[centre] theta(0, t*): FEM {c_A:.7f} K, closed form {c_cf_t:.7f} K; "
         f"steady theta(0): FEM {c_Q:.7f} K, closed form {c_cf_inf:.7f} K")
    _say(f"[centre] theta(t*)/theta(inf): FEM {c_A / c_Q:.7f}, closed form "
         f"{c_cf_t / c_cf_inf:.7f}")
    _say(f"[xdmf] steady field is the {'FEM steady solve' if same else 'closed-form steady theta'}; "
         f"wrote {xdmf_path}")
    _say(f"[time] EM solve + mesh {t_em:.1f} s; two transients {t_th:.2f} s")

    assert once, "(Once) a second EM solve happened"
    assert n_sub == MAT10_SUB_CELLS, f"(UT) sub cells {n_sub}"
    assert layA and layS, "(I-dir) owned layouts differ"
    assert idA <= L0_BAND and idS <= L0_BAND, f"(I-dir) A {idA:.3e} S {idS:.3e}"
    assert err_A <= A_BAND, f"(A_t) {err_A:.4e} > {A_BAND}"
    assert f"{err_A:.4e}" == MAT12_A_T_PRINT, f"(A_t) print {err_A:.4e}"
    assert s2 <= S2_BAND, f"(S2) {s2:.4e} > {S2_BAND:.0e}"
    assert f"{s2:.4e}" == MAT12_S2_PRINT, f"(S2) print {s2:.4e}"
    assert eA <= E_BAND and eS <= E_BAND, f"(E) A {eA:.3e} S {eS:.3e}"
    assert err_qs > A_BAND, f"(c-qs) {err_qs:.4e} did not miss {A_BAND}"
    assert imag == 0.0, f"max|Im theta| {imag:.3e}"

    _say(f"\nAll anchors and controls hold. Total elapsed {time.perf_counter() - started:.1f} s.")


if __name__ == "__main__":
    main()
