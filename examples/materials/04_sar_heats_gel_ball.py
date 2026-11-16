"""Example (`EX-67`): SAR heats the gel ball — the steady temperature, the rise
curve, and the calorimetric slope SAR/c.

Every other example stops at SAR; none imports ``fem_em_solver.thermal``. This
one takes `mat:2`'s lossy sphere (R = 0.01 m, sigma 0.57 S/m, eps_r 78, 64 MHz,
imposed uniform field), turns G * 1/2 sigma |E|^2 into a DG0 heat source q_h on
the sphere's own submesh, and shows three readings:

* the Robin-cooled **steady** field (h = 100 W/m^2/K, Bi = 5/3) — `MAT-10`;
* the **rise** from theta = 0 to t* = 0.2 R^2 rho c / k — `MAT-12`, where the
  steady answer stands in about twice too hot;
* the **first-instant slope**: one backward-Euler step, volume mean, equals
  mean SAR / c — the identity phantom calorimetry rests on — `MAT-12` (S2).

**It asserts, and it does not re-implement** (the `ANS-1` rule). The fixture,
helpers and bands are imported from
``tests/validation/test_mat10_sar_to_heat_sphere.py`` and
``tests/validation/test_mat12_transient_sar_heating_sphere.py``; the EM solve is
``_solve_lossy_sphere`` from ``test_lossy_sphere_sar.py``. The one change to a
gate module is additive and disclosed: ``_cn(sub, ft, q, snapshot_steps=())``
passes the keyword through to ``solve_transient_conduction`` (`EX-67-mat12-regate`).

Asserted, each as its gate writes it: (B) transfer <= ``B_BAND``; steady (P)
<= ``P_BAND``, (A) <= ``A_BAND``, (C) <= ``C_BAND``; transient (P_t) <=
``P_BAND``, (A_t) <= ``A_BAND``; (S2) <= ``S2_BAND``; (E) <= ``E_BAND``; (inf)
<= ``INF_BAND``. Negative controls, asserted only as "misses ``A_BAND``": (a)
the Re-E-only source through the steady solve; (b) the steady theta_h[q_h]
standing in for theta(t*).

Errors are on theta = T - T_inf, never on T: the steady solves use
``t_inf=T_INF`` (`MAT-10`), the transient ones 0 (`MAT-12`).

Run it::

    ./run_examples.sh -e mat:4 -n 2 -t 180

Output (on the submesh): ``examples/materials/paraview_output/
materials_04_sar_heats_gel_ball_steady_combined.xdmf`` (DG0 ``q_h``, P1
``theta_steady``), ``..._transient_combined.xdmf`` (P1 ``theta`` at steps 0,
4, ..., 40) and ``..._rise_curve.csv``.
"""

from __future__ import annotations

import csv
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

from fem_em_solver.io.paraview_utils import (  # noqa: E402
    adopt_host_ownership,
    write_xdmf_time_series,
    write_xdmf_with_tags,
)
from fem_em_solver.post.sar import uniform_sphere_sar_closed_form  # noqa: E402
from fem_em_solver.post.setup_figure import write_setup_figure  # noqa: E402
from fem_em_solver.thermal import (  # noqa: E402
    solve_steady_conduction,
    solve_transient_conduction,
)
from fem_em_solver.utils.constants import EPSILON_0  # noqa: E402

from tests.validation.test_lossy_sphere_sar import (  # noqa: E402
    E0,
    EPSILON_R_SPHERE,
    OMEGA,
    RHO_KG_M3,
    SIGMA_HIGH,
    SPHERE_TAG,
    _solve_lossy_sphere,
)
from tests.validation.test_mat10_sar_to_heat_sphere import (  # noqa: E402
    A_BAND,
    B_BAND,
    C_BAND,
    G,
    H,
    K,
    P_BAND,
    R,
    ROBIN_TAG,
    T_INF,
    _allsum,
    _cell_average_dg0,
    _interior_indicator,
    _interior_rel,
    _robin_tags,
    _sphere_submesh,
    _theta_cf,
    _to_sub,
)
from tests.validation.test_mat12_transient_sar_heating_sphere import (  # noqa: E402
    DT_1,
    E_BAND,
    INF_BAND,
    N_MAIN,
    RHO_C,
    S2_BAND,
    T_R,
    T_STAR,
    _be_one,
    _cn,
    _exact_theta,
    _mean,
    _whole_rel,
)

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"  # committed, EX-57
BASENAME = "materials_04_sar_heats_gel_ball"
FRAMES = tuple(range(0, N_MAIN + 1, 4))

#: Gate records — `20260930T021802Z_MAT-10-step1-W1.log:157-167` and
#: `20260930T110230Z_MAT-12-W1.log:158-170`. Printed beside this run, a
#: difference reported and never asserted.
RECORD = {
    "B": 2.807e-15, "P": 3.1923e-03, "A": 3.2464e-02, "C": 1.233e-13,
    "P_t": 1.6845e-03, "A_t": 3.4708e-02, "S2": 4.9994e-05,
    "slope_mean": 4.958818e-04, "sar_c": 4.959066e-04, "inf": 8.856e-16,
    "ctl_a": 7.9323e-01, "ctl_b": 1.0108e00,
}


def _say(comm, msg):
    if comm.rank == 0:
        print(msg, flush=True)


def _rec(key, value, fmt=".4e"):
    return f"{value:{fmt}} (record {RECORD[key]:{fmt}}, diff {value - RECORD[key]:+.1e})"


def main() -> None:
    comm = MPI.COMM_WORLD
    started = time.perf_counter()

    if not np.issubdtype(np.dtype(default_scalar_type), np.complexfloating):
        raise RuntimeError(
            "this example needs the complex DolfinX build "
            "(source /usr/local/bin/dolfinx-complex-mode); the runner's `mat:` "
            "group sources it automatically"
        )

    _say(comm, "=" * 78)
    _say(comm, "EX-67 -- SAR heats the gel ball: steady field, rise curve, slope SAR/c")
    _say(comm, "=" * 78)
    _say(comm, f"\n[gate] bands imported: B {B_BAND:.0e}, P {P_BAND:.0e}, A {A_BAND}, "
               f"C {C_BAND:.0e}, S2 {S2_BAND:.0e}, E {E_BAND:.0e}, inf {INF_BAND:.0e}")

    # ---- EM: mat:2's lossy sphere, one solve -------------------------------
    t = time.perf_counter()
    run = _solve_lossy_sphere(SIGMA_HIGH, R / 10.0, R / 5.0, return_fields=True)
    msh, cell_tags = run["mesh"], run["cell_tags"]
    E, sigma = run["e_complex"], run["sigma_field"]
    t_em = time.perf_counter() - t

    write_setup_figure(
        msh,
        cell_tags,
        FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names={SPHERE_TAG: "sphere (gel, heated)", 2: "air"},
        hide_tags=(2,),
        translucent_tags=(SPHERE_TAG,),
        slice_normal=(0.0, 0.0, 1.0),
        slice_origin=(0.0, 0.0, 0.0),
        title="mat:4 -- imposed-field SAR heats the lossy gel ball (MAT-10 / MAT-12 fixture), Robin h = 100",
        comm=comm,
    )

    cf = uniform_sphere_sar_closed_form(
        e0=E0, epsilon_r=EPSILON_R_SPHERE, sigma=SIGMA_HIGH, omega=OMEGA,
        rho=RHO_KG_M3, epsilon_0=EPSILON_0,
    )
    q_cf = G * RHO_KG_M3 * cf["sar_w_per_kg"]

    sub, parent_ids = _sphere_submesh(msh, cell_tags)
    ft = _robin_tags(sub)
    n_parent = msh.topology.index_map(3).size_global
    n_sub = sub.topology.index_map(3).size_global
    _say(comm, f"\n[mesh] parent {n_parent} cells, sphere submesh {n_sub} cells; EM solve {t_em:.1f} s")

    q_expr = G * 0.5 * sigma * ufl.inner(E, E)
    Wp, q_par = _cell_average_dg0(msh, cell_tags, q_expr)
    q_h = _to_sub(Wp, q_par, sub, parent_ids)
    q_h.name = "q_h"
    q_re_expr = G * 0.5 * sigma * ufl.inner(ufl.real(E), ufl.real(E))
    _, q_re_par = _cell_average_dg0(msh, cell_tags, q_re_expr)
    q_re = _to_sub(Wp, q_re_par, sub, parent_ids)

    one = fem.Constant(sub, default_scalar_type(1.0))
    vol = float(np.real(_allsum(comm, one * ufl.dx(domain=sub))))
    qbar = _mean(sub, q_h, vol)
    _say(comm, f"[source] q_cf = {q_cf:.6f} W/m^3 (closed form), mean q_h = {qbar:.6f} W/m^3 "
               f"(q_h/q_cf - 1 = {qbar / q_cf - 1.0:+.4e}); mean SAR/c = {qbar / RHO_C:.6e} K/s")

    # (B) transfer
    int_q_sub = float(np.real(_allsum(comm, q_h * ufl.dx(domain=sub))))
    p_parent = G * run["dissipated_power_w"]
    err_B = abs(int_q_sub - p_parent) / abs(p_parent)
    _say(comm, f"\n(B) int_sub q_h = {int_q_sub:.15e} W vs G*P_diss = {p_parent:.15e} W: "
               f"rel {_rec('B', err_B, '.3e')}  ASSERTED <= {B_BAND:.0e}")

    # ---- steady (MAT-10) -----------------------------------------------------
    chi = _interior_indicator(sub)
    theta_ss = _theta_cf(sub, q_cf)
    res_P = solve_steady_conduction(sub, ft, k=K, q=q_cf, h=H, t_inf=T_INF,
                                    robin_tag=ROBIN_TAG, degree=1)
    err_P = _interior_rel(sub, chi, (res_P.temperature - T_INF) - theta_ss, theta_ss)
    res_A = solve_steady_conduction(sub, ft, k=K, q=q_h, h=H, t_inf=T_INF,
                                    robin_tag=ROBIN_TAG, degree=1)
    err_A = _interior_rel(sub, chi, (res_A.temperature - T_INF) - theta_ss, theta_ss)
    err_C = abs(res_A.convected_power - res_A.source_power) / abs(res_A.source_power)
    _say(comm, f"[steady] (P) uniform q_cf: interior rel L2 {_rec('P', err_P)}  ASSERTED <= {P_BAND:.0e}")
    _say(comm, f"[steady] (A) q_h: interior rel L2 {_rec('A', err_A)}  ASSERTED <= {A_BAND}")
    _say(comm, f"[steady] (C) balance {_rec('C', err_C, '.3e')}  ASSERTED <= {C_BAND:.0e}")

    # ---- transient (MAT-12), theta = T, t_inf = T(0) = 0 ------------------------
    theta_star = _exact_theta(sub, T_STAR, q_cf)
    res_Pt = _cn(sub, ft, q_cf)
    err_Pt = _interior_rel(sub, chi, res_Pt.temperature - theta_star, theta_star)
    res_At = _cn(sub, ft, q_h, snapshot_steps=FRAMES)
    err_At = _interior_rel(sub, chi, res_At.temperature - theta_star, theta_star)
    _say(comm, f"\n[transient] t_R = {T_R:.4f} s, t* = {T_STAR:.4f} s, CN+Rannacher N = {N_MAIN}")
    _say(comm, f"[transient] (P_t) uniform q_cf at t*: {_rec('P_t', err_Pt)}  ASSERTED <= {P_BAND:.0e}")
    _say(comm, f"[transient] (A_t) q_h at t*: {_rec('A_t', err_At)}  ASSERTED <= {A_BAND}")

    # (S2) calorimetric slope: one BE step of dt_1 from theta = 0
    res_S2 = _be_one(sub, ft, q_h)
    thbar1 = _mean(sub, res_S2.temperature, vol)
    ratio_S2 = RHO_C * thbar1 / (DT_1 * qbar)
    err_S2 = abs(ratio_S2 - 1.0)
    _say(comm, f"[slope] (S2) rho c mean(theta^1)/(dt_1 mean q_h) = {ratio_S2:.9f}: "
               f"|ratio - 1| {_rec('S2', err_S2)}  ASSERTED <= {S2_BAND:.0e}")
    _say(comm, f"[slope] mean slope {thbar1 / DT_1:.6e} K/s (record {RECORD['slope_mean']:.6e}) "
               f"on mean SAR/c {qbar / RHO_C:.6e} K/s (record {RECORD['sar_c']:.6e})")

    # (inf) long backward Euler against the steady solve with the same q_h
    res_inf = solve_transient_conduction(
        sub, ft, k=K, rho_c=RHO_C, q=q_h, h=H, t_inf=0.0, t_init=0.0,
        robin_tag=ROBIN_TAG, dt=5.0 * T_R, n_steps=20, theta=1.0, degree=1,
    )
    st0 = solve_steady_conduction(sub, ft, k=K, q=q_h, h=H, t_inf=0.0,
                                  robin_tag=ROBIN_TAG, degree=1)
    err_inf = _whole_rel(sub, res_inf.temperature - st0.temperature, st0.temperature)
    _say(comm, f"[transient] (inf) BE 20 x 5 t_R vs steady[q_h]: {_rec('inf', err_inf, '.3e')}  "
               f"ASSERTED <= {INF_BAND:.0e}")

    balances = {
        "P_t": res_Pt.max_balance_residual, "A_t": res_At.max_balance_residual,
        "S2": res_S2.max_balance_residual, "inf": res_inf.max_balance_residual,
    }
    worst = max(balances, key=balances.get)
    _say(comm, f"[transient] (E) worst per-step balance {balances[worst]:.3e} on run {worst}  "
               f"ASSERTED <= {E_BAND:.0e}")

    # ---- negative controls ------------------------------------------------------
    res_a = solve_steady_conduction(sub, ft, k=K, q=q_re, h=H, t_inf=T_INF,
                                    robin_tag=ROBIN_TAG, degree=1)
    err_a = _interior_rel(sub, chi, (res_a.temperature - T_INF) - theta_ss, theta_ss)
    err_b = _interior_rel(sub, chi, st0.temperature - theta_star, theta_star)
    _say(comm, f"\n[control] (a) Re-E-only source, steady: {_rec('ctl_a', err_a)}  "
               f"ASSERTED > {A_BAND}")
    _say(comm, f"[control] (b) steady theta_h[q_h] for theta(t*): {_rec('ctl_b', err_b)}  "
               f"ASSERTED > {A_BAND}")

    # ---- the rise curve (printed and written, never gated) -----------------------
    dt = T_STAR / N_MAIN
    steady_mean = _mean(sub, st0.temperature, vol)
    rows = []
    _say(comm, "\n[rise] step   t [s]    mean theta_h [K]   calorimetric q_h t/(rho c) [K]   "
               "steady mean [K]  (printed, not gated)")
    for step in FRAMES:
        f = res_At.snapshots[step]
        m = _mean(sub, f, vol)
        line = qbar * step * dt / RHO_C
        rows.append((step, step * dt, m, line, steady_mean))
        _say(comm, f"[rise] {step:4d}  {step * dt:8.3f}   {m:.6e}       {line:.6e}"
                   f"                   {steady_mean:.6e}")

    # ---- assertions ---------------------------------------------------------------
    assert err_B <= B_BAND, f"(B) {err_B:.3e}"
    assert err_P <= P_BAND, f"(P) {err_P:.4e}"
    assert err_A <= A_BAND, f"(A) {err_A:.4e}"
    assert err_C <= C_BAND, f"(C) {err_C:.3e}"
    assert err_Pt <= P_BAND, f"(P_t) {err_Pt:.4e}"
    assert err_At <= A_BAND, f"(A_t) {err_At:.4e}"
    assert err_S2 <= S2_BAND, f"(S2) {err_S2:.4e}"
    assert balances[worst] <= E_BAND, f"(E) {worst} {balances[worst]:.3e}"
    assert err_inf <= INF_BAND, f"(inf) {err_inf:.3e}"
    assert err_a > A_BAND, f"control (a) {err_a:.4e} did not miss {A_BAND}"
    assert err_b > A_BAND, f"control (b) {err_b:.4e} did not miss {A_BAND}"

    # ---- outputs (on the submesh) -----------------------------------------------------
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    V1 = res_A.temperature.function_space
    theta_steady = fem.Function(V1, name="theta_steady")
    theta_steady.x.array[:] = res_A.temperature.x.array - T_INF
    theta_steady.x.scatter_forward()
    steady_path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_steady_combined", sub, None,
        {"q_h": q_h, "theta_steady": theta_steady}, comm=comm,
    )
    frames = []
    for step in FRAMES:
        f = res_At.snapshots[step]
        f.name = "theta"
        frames.append((step * dt, {"theta": f}))
    series_path, _ = write_xdmf_time_series(
        OUTPUT_DIR / f"{BASENAME}_transient_combined", sub, None, frames, comm=comm,
    )
    table = OUTPUT_DIR / f"{BASENAME}_rise_curve.csv"
    if comm.rank == 0:
        with table.open("w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["step", "t_s", "mean_theta_h_K", "calorimetric_K", "steady_mean_K"])
            for step, tt, m, line, sm in rows:
                w.writerow([step, f"{tt:.6f}", f"{m:.9e}", f"{line:.9e}", f"{sm:.9e}"])
    adopt_host_ownership(OUTPUT_DIR, comm=comm)
    _say(comm, f"\n[output] {steady_path} -- DG0 q_h, P1 theta_steady (theta = T - T_inf)")
    _say(comm, f"[output] {series_path} -- P1 theta at steps {FRAMES[0]}, {FRAMES[1]}, ..., {FRAMES[-1]}")
    _say(comm, f"[output] {table}")

    _say(comm, f"\n[timing] total {time.perf_counter() - started:.1f} s at -n {comm.size} "
               f"(EM {t_em:.1f} s)")
    _say(comm, "\nAll imported gates hold.")


if __name__ == "__main__":
    main()
