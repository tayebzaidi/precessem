"""Example (`EX-73`): SAR heats a perfused ball — the rise with the Pennes sink
switched on.

`mat:4` heats the lossy gel ball with SAR but has no sink; `mat:6` has the sink
but a uniform source. This is the first example with both. It takes `mat:4`'s
fixture (R = 0.01 m, sigma 0.57 S/m, eps_r 78, 64 MHz, imposed uniform field,
Robin h = 100 W/m^2/K, arterial = ambient = T(0) = 20 C), turns G * 1/2 sigma
|E|^2 into the DG0 source q_h on the sphere submesh, and runs `MAT-14`'s
Crank-Nicolson rise to t* = 0.2 R^2 rho c / k with w = 24 000 W/m^3/K beside
w = 0. w = 24 000 (m R = 2) is `MAT-13`'s dimensionless twin -- a numerical
choice, **not a tissue value**.

**It asserts, and it does not re-implement** (the `ANS-1` rule). Everything is
imported from ``tests/validation/test_mat14_sar_pennes_transient_sphere.py``
(`MAT-14`), ``test_mat10_sar_to_heat_sphere.py``,
``test_mat12_transient_sar_heating_sphere.py`` and ``test_lossy_sphere_sar.py``;
the names `MAT-14` re-imports are asserted ``is`` the originals. No gate edit.

Asserted, each as `MAT-14` writes it: (P_w') <= ``P_BAND``; (A_w) <=
``A_BAND``; (E) every run <= ``E_BAND`` and recomputed from snapshots 39 / 40
<= ``E_BAND``; (inf) <= ``INF_BAND``; (w -> 0) the formatted errors equal
`MAT-14`'s ``MAT12_P_RECORD`` / ``MAT12_A_RECORD``. Negative control, asserted
only as a band miss: (a') w dropped with q_h, against the perfused series,
misses ``A_BAND``.

Errors are on theta = T - T_inf (T_inf = T(0) = 20 in every main run).

Run it::

    ./run_examples.sh -e mat:7 -n 2 -t 180

Output (on the submesh): ``examples/materials/paraview_output/
materials_07_sar_heats_perfused_ball_transient_combined.xdmf`` (P1 theta at
steps 0, 4, ..., 40 for w = 24 000 and w = 0, plus DG0 ``q_h``),
``..._steady_combined.xdmf`` (perfused and unperfused steady theta under q_h)
and ``..._rise_curve.csv``.
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
from fem_em_solver.post.setup_figure import write_setup_figure  # noqa: E402
from fem_em_solver.thermal import (  # noqa: E402
    solve_steady_conduction,
    solve_transient_conduction,
)

from tests.validation import test_lossy_sphere_sar as lss  # noqa: E402
from tests.validation import test_mat10_sar_to_heat_sphere as m10  # noqa: E402
from tests.validation import test_mat12_transient_sar_heating_sphere as m12  # noqa: E402
from tests.validation import test_mat14_sar_pennes_transient_sphere as m14  # noqa: E402
from tests.validation.test_lossy_sphere_sar import (  # noqa: E402
    SIGMA_HIGH,
    SPHERE_TAG,
    _solve_lossy_sphere,
)
from tests.validation.test_mat10_sar_to_heat_sphere import (  # noqa: E402
    A_BAND,
    G,
    H,
    K,
    P_BAND,
    R,
    ROBIN_TAG,
    _allsum,
    _cell_average_dg0,
    _interior_indicator,
    _interior_rel,
    _robin_tags,
    _sphere_submesh,
    _to_sub,
)
from tests.validation.test_mat12_transient_sar_heating_sphere import (  # noqa: E402
    E_BAND,
    INF_BAND,
    N_MAIN,
    RHO_C,
    T_R,
    T_STAR,
    _cn,
    _exact_theta,
    _whole_rel,
)
from tests.validation.test_mat14_sar_pennes_transient_sphere import (  # noqa: E402
    MAT12_A_RECORD,
    MAT12_P_RECORD,
    T_INF,
    W,
    _cn_w,
    _exact_theta_w,
    _q_cf,
)

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"  # committed, EX-57
BASENAME = "materials_07_sar_heats_perfused_ball"
FRAMES = tuple(range(0, N_MAIN + 1, 4))

#: Gate records — `20261001T004928Z_MAT-14-W1.log:162-171`. Printed beside this
#: run, a difference reported and never asserted.
RECORD = {
    "P_w": 1.6089e-03, "A_w": 3.4979e-02, "E": 3.308e-14, "inf": 3.105e-13,
    "ctl_a1": 4.3826e-01, "split": (0.2134, 0.3888, 0.3978),
}


def _say(comm, msg):
    if comm.rank == 0:
        print(msg, flush=True)


def _rec(key, value, fmt=".4e"):
    return f"{value:{fmt}} (record {RECORD[key]:{fmt}}, diff {value - RECORD[key]:+.1e})"


def _theta_of(T_h, name):
    """theta = T - T_inf as a new function in T_h's space."""
    f = fem.Function(T_h.function_space, name=name)
    f.x.array[:] = T_h.x.array - T_INF
    f.x.scatter_forward()
    return f


def _mean(comm, sub, f, vol):
    return float(np.real(_allsum(comm, f * ufl.dx(domain=sub)))) / vol


def main() -> None:
    comm = MPI.COMM_WORLD
    started = time.perf_counter()

    if not np.issubdtype(np.dtype(default_scalar_type), np.complexfloating):
        raise RuntimeError(
            "this example needs the complex DolfinX build "
            "(source /usr/local/bin/dolfinx-complex-mode); the runner's `mat:` "
            "group sources it automatically"
        )

    # Structural check: every name MAT-14 re-imports is the original object.
    assert m14.A_BAND is A_BAND and m14.P_BAND is P_BAND
    assert m14.G is G and m14.H is H and m14.K is K and m14.R is R
    assert m14.ROBIN_TAG is ROBIN_TAG
    assert m14._allsum is _allsum and m14._cell_average_dg0 is _cell_average_dg0
    assert m14._interior_indicator is _interior_indicator
    assert m14._interior_rel is _interior_rel and m14._robin_tags is _robin_tags
    assert m14._sphere_submesh is _sphere_submesh and m14._to_sub is _to_sub
    assert m14.E_BAND is E_BAND and m14.INF_BAND is INF_BAND
    assert m14.N_MAIN is N_MAIN and m14.RHO_C is RHO_C
    assert m14.T_R is T_R and m14.T_STAR is T_STAR
    assert m14._cn is _cn and m14._exact_theta is _exact_theta
    assert m14._whole_rel is _whole_rel
    assert m14.SIGMA_HIGH is SIGMA_HIGH and m14._solve_lossy_sphere is _solve_lossy_sphere
    assert lss.SIGMA_HIGH is SIGMA_HIGH and m10.A_BAND is A_BAND and m12.E_BAND is E_BAND

    _say(comm, "=" * 78)
    _say(comm, "EX-73 -- SAR heats a perfused ball: the rise with the Pennes sink on")
    _say(comm, "=" * 78)
    _say(comm, f"\n[gate] bands imported: P {P_BAND:.0e}, A {A_BAND}, E {E_BAND:.0e}, "
               f"inf {INF_BAND:.0e}; w = {W:g} W/m^3/K (m R = {np.sqrt(W / K) * R:.4f}; "
               f"numerical choice, not a tissue value), t* = {T_STAR:.4f} s, N = {N_MAIN}")

    # ---- EM: mat:2's lossy sphere, one solve reused for every leg ---------------
    t = time.perf_counter()
    run = _solve_lossy_sphere(SIGMA_HIGH, R / 10.0, R / 5.0, return_fields=True)
    msh, cell_tags = run["mesh"], run["cell_tags"]
    E, sigma = run["e_complex"], run["sigma_field"]
    t_em = time.perf_counter() - t

    write_setup_figure(
        msh,
        cell_tags,
        FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names={SPHERE_TAG: "sphere (gel, heated, Pennes sink w = 24 000)", 2: "air"},
        hide_tags=(2,),
        translucent_tags=(SPHERE_TAG,),
        slice_normal=(0.0, 0.0, 1.0),
        slice_origin=(0.0, 0.0, 0.0),
        title="mat:7 -- imposed-field SAR heats the perfused gel ball (MAT-14 fixture), Robin h = 100",
        comm=comm,
    )

    q_cf, cf = _q_cf()
    sub, parent_ids = _sphere_submesh(msh, cell_tags)
    ft = _robin_tags(sub)
    n_parent = msh.topology.index_map(3).size_global
    n_sub = sub.topology.index_map(3).size_global
    _say(comm, f"\n[mesh] parent {n_parent} cells, sphere submesh {n_sub} cells; EM solve {t_em:.1f} s")

    q_expr = G * 0.5 * sigma * ufl.inner(E, E)
    Wp, q_par = _cell_average_dg0(msh, cell_tags, q_expr)
    q_h = _to_sub(Wp, q_par, sub, parent_ids)
    q_h.name = "q_h"

    one = fem.Constant(sub, default_scalar_type(1.0))
    vol = float(np.real(_allsum(comm, one * ufl.dx(domain=sub))))
    qbar = _mean(comm, sub, q_h, vol)
    _say(comm, f"[source] q_cf = {q_cf:.6f} W/m^3, mean q_h = {qbar:.6f} W/m^3 "
               f"(q_h/q_cf - 1 = {qbar / q_cf - 1.0:+.4e})")

    chi = _interior_indicator(sub)
    theta_w = _exact_theta_w(sub, T_STAR, q_cf)

    # ---- (P_w') and (A_w): perfused CN rises ----------------------------------------
    res_P = _cn_w(sub, ft, q_cf)
    err_P = _interior_rel(sub, chi, res_P.temperature - T_INF - theta_w, theta_w)
    res_A = _cn_w(sub, ft, q_h, snapshot_steps=FRAMES + (N_MAIN - 1,))
    err_A = _interior_rel(sub, chi, res_A.temperature - T_INF - theta_w, theta_w)
    _say(comm, f"\n[transient] CN + Rannacher, N = {N_MAIN}, t* = {T_STAR:.4f} s, w = {W:g}")
    _say(comm, f"[transient] (P_w') uniform q_cf: interior rel L2 {_rec('P_w', err_P)}  "
               f"ASSERTED <= {P_BAND:.0e}")
    _say(comm, f"[transient] (A_w) q_h: interior rel L2 {_rec('A_w', err_A)}  "
               f"ASSERTED <= {A_BAND}")

    # ---- (E) recomputed from snapshots 39 / 40 --------------------------------------
    T39, T40 = res_A.snapshots[N_MAIN - 1], res_A.snapshots[N_MAIN]
    ds = ufl.Measure("ds", domain=sub, subdomain_data=ft)
    dx = ufl.dx(domain=sub)
    dt = T_STAR / N_MAIN
    stor = float(np.real(_allsum(comm, RHO_C * (T40 - T39) / dt * dx)))
    conv = float(np.real(_allsum(comm, H * (0.5 * T40 + 0.5 * T39 - T_INF) * ds(ROBIN_TAG))))
    perf = float(np.real(_allsum(comm, W * (0.5 * T40 + 0.5 * T39 - T_INF) * dx)))
    src = res_A.source_power
    bal = abs(stor + conv + perf - src) / src
    _say(comm, f"[transient] (E) recomputed last-step balance {_rec('E', bal, '.3e')}  "
               f"ASSERTED <= {E_BAND:.0e}")
    _say(comm, f"[transient] power split at t* storage / convected / perfused = "
               f"{stor / src:.4f} / {conv / src:.4f} / {perf / src:.4f} "
               f"(record {RECORD['split'][0]:.4f} / {RECORD['split'][1]:.4f} / "
               f"{RECORD['split'][2]:.4f}; printed)")

    # ---- (inf) long BE with q_h and w vs the perfused steady solve ----------------------
    res_inf = solve_transient_conduction(
        sub, ft, k=K, rho_c=RHO_C, q=q_h, h=H, t_inf=T_INF, t_init=T_INF,
        robin_tag=ROBIN_TAG, dt=5.0 * T_R, n_steps=20, theta=1.0, degree=1,
        perfusion_w_per_m3_k=W,
    )
    st_h = solve_steady_conduction(sub, ft, k=K, q=q_h, h=H, t_inf=T_INF,
                                   robin_tag=ROBIN_TAG, degree=1, perfusion_w_per_m3_k=W)
    err_inf = _whole_rel(sub, res_inf.temperature - st_h.temperature,
                         st_h.temperature - T_INF)
    _say(comm, f"[transient] (inf) BE 20 x 5 t_R vs steady[q_h, w]: {_rec('inf', err_inf, '.3e')}  "
               f"ASSERTED <= {INF_BAND:.0e}")

    # ---- (w -> 0): MAT-12's own legs ----------------------------------------------------------
    theta_0 = _exact_theta(sub, T_STAR, q_cf)
    r12P = _cn(sub, ft, q_cf)
    r12A = _cn(sub, ft, q_h)
    e12P = _interior_rel(sub, chi, r12P.temperature - theta_0, theta_0)
    e12A = _interior_rel(sub, chi, r12A.temperature - theta_0, theta_0)
    _say(comm, f"[w->0] MAT-12 _cn legs: q_cf {e12P:.4e} (MAT-14 record {MAT12_P_RECORD}), "
               f"q_h {e12A:.4e} (record {MAT12_A_RECORD})  ASSERTED string-equal")

    # ---- negative control (a'): w dropped with q_h, against the perfused series ---------------
    r_a1 = _cn_w(sub, ft, q_h, w=0.0, snapshot_steps=FRAMES)
    e_a1 = _interior_rel(sub, chi, r_a1.temperature - T_INF - theta_w, theta_w)
    _say(comm, f"[control] (a') w dropped, q_h, against the perfused series: "
               f"{_rec('ctl_a1', e_a1)}  ASSERTED > {A_BAND}")

    balances = {
        "P_w'": res_P.max_balance_residual, "A_w": res_A.max_balance_residual,
        "inf": res_inf.max_balance_residual, "w0_P": r12P.max_balance_residual,
        "w0_A": r12A.max_balance_residual, "ctl_a'": r_a1.max_balance_residual,
    }
    worst = max(balances, key=balances.get)
    _say(comm, f"[transient] (E) worst per-step balance {balances[worst]:.3e} on run {worst}  "
               f"ASSERTED <= {E_BAND:.0e}")

    # ---- the rise curve (printed and written, never gated) -------------------------------------
    st_0 = solve_steady_conduction(sub, ft, k=K, q=q_h, h=H, t_inf=T_INF,
                                   robin_tag=ROBIN_TAG, degree=1)
    steady_w = _mean(comm, sub, _theta_of(st_h.temperature, "t"), vol)
    steady_0 = _mean(comm, sub, _theta_of(st_0.temperature, "t"), vol)
    rows = []
    _say(comm, f"\n[rise] step   t [s]    mean theta_h [K] w={W:g}   w=0        "
               f"steady mean {steady_w:.4f} K perfused / {steady_0:.4f} K unperfused "
               f"(printed, not gated)")
    for step in FRAMES:
        m_w = _mean(comm, sub, _theta_of(res_A.snapshots[step], "t"), vol)
        m_0 = _mean(comm, sub, _theta_of(r_a1.snapshots[step], "t"), vol)
        rows.append((step, step * dt, m_w, m_0))
        _say(comm, f"[rise] {step:4d}  {step * dt:8.3f}   {m_w:.6e}       {m_0:.6e}")

    # ---- assertions ---------------------------------------------------------------------------------
    assert err_P <= P_BAND, f"(P_w') {err_P:.4e}"
    assert err_A <= A_BAND, f"(A_w) {err_A:.4e}"
    assert balances[worst] <= E_BAND, f"(E) {worst} {balances[worst]:.3e}"
    assert bal <= E_BAND, f"(E) recomputed {bal:.3e}"
    assert err_inf <= INF_BAND, f"(inf) {err_inf:.3e}"
    assert f"{e12P:.4e}" == MAT12_P_RECORD, f"(w->0) q_cf {e12P:.4e}"
    assert f"{e12A:.4e}" == MAT12_A_RECORD, f"(w->0) q_h {e12A:.4e}"
    assert e_a1 > A_BAND, f"control (a') {e_a1:.4e} did not miss {A_BAND}"

    # ---- outputs (on the submesh; CG1 exports are interpolants, nothing asserted) ----------------
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    steady_path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_steady_combined", sub, None,
        {"theta_perfused": _theta_of(st_h.temperature, "theta_perfused"),
         "theta_unperfused": _theta_of(st_0.temperature, "theta_unperfused")}, comm=comm,
    )
    frames = [
        (step * dt, {
            "theta_w24000": _theta_of(res_A.snapshots[step], "theta_w24000"),
            "theta_w0": _theta_of(r_a1.snapshots[step], "theta_w0"),
            "q_h": q_h,
        })
        for step in FRAMES
    ]
    series_path, _ = write_xdmf_time_series(
        OUTPUT_DIR / f"{BASENAME}_transient_combined", sub, None, frames, comm=comm,
    )
    table = OUTPUT_DIR / f"{BASENAME}_rise_curve.csv"
    if comm.rank == 0:
        with table.open("w", newline="") as fh:
            wr = csv.writer(fh)
            wr.writerow(["step", "t_s", "mean_theta_w24000_K", "mean_theta_w0_K",
                         "steady_mean_w24000_K", "steady_mean_w0_K"])
            for step, tt, m_w, m_0 in rows:
                wr.writerow([step, f"{tt:.6f}", f"{m_w:.9e}", f"{m_0:.9e}",
                             f"{steady_w:.9e}", f"{steady_0:.9e}"])
    adopt_host_ownership(OUTPUT_DIR, comm=comm)
    _say(comm, f"\n[output] {steady_path} -- P1 theta_perfused, theta_unperfused (theta = T - T_inf)")
    _say(comm, f"[output] {series_path} -- P1 theta_w24000, theta_w0 and DG0 q_h at steps "
               f"{FRAMES[0]}, {FRAMES[1]}, ..., {FRAMES[-1]}")
    _say(comm, f"[output] {table}")

    _say(comm, f"\n[timing] total {time.perf_counter() - started:.1f} s at -n {comm.size} "
               f"(EM {t_em:.1f} s)")
    _say(comm, "\nAll imported gates hold.")


if __name__ == "__main__":
    main()
