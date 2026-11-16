"""`MAT-12` — transient SAR heating of the gel ball (§9 item 100).

`MAT-10`'s fixture (one EM solve of `MAT-4`'s lossy sphere, R = 0.01 m, the
tag-1 submesh, q_h = DG0 cell average of G·½σ|E|², h = 100 W/m²/K, Bi 5/3,
k 0.6) driven through `MAT-9`'s θ-scheme with ρc = 4.18e6 J/m³/K, so
t_R = R²ρc/k = 696.67 s.  **T∞ = T(0) = 0 in every solve** (θ = T).

Closed form: θ(r, t) = `robin_sphere_transient_temperature` with
q_cf = G·ρ·SAR_cf = 2002.09 W/m³ (series kept in numpy, interpolated into
P2); at t* = 0.2 t_R = 139.33 s: θ(0) = 0.059564 K, θ(R) = 0.035599 K.

Anchors (asserted): (P_t) uniform q_cf, CN + Rannacher, N = 40, at t*:
interior relative L2 on θ ≤ 1 % (`MAT-10`'s ``P_BAND``); (A_t) the same run
with q_h ≤ 11 % (``A_BAND``, composed, never tightened); (S) one backward-Euler
step Δt₁ = 1e-5 t_R from θ = 0: (S1) uniform q_cf, interior
‖ρcθ¹/Δt₁ − q_cf‖/‖q_cf‖ ≤ 1e-3; (S2) q_h, volume means
|ρc·θ̄¹/(Δt₁·q̄_h) − 1| ≤ 1e-4 — the calorimetric identity; (E) per-step
balance ≤ 1e-9 on every run; (∞) BE 20 × 5 t_R with q_h equals
``solve_steady_conduction(q=q_h)`` to ≤ 1e-8.
Negative controls (asserted only as "misses (A_t)'s 11 %"): (a) q from Re E
only (predicted ≈ 79.1 %); (b) quasi-static, the steady θ_h[q_h] standing in
for θ(t*) (predicted ≈ 101 %).
"""

from __future__ import annotations

import time

import numpy as np
import ufl
from dolfinx import default_scalar_type, fem
from mpi4py import MPI

from fem_em_solver.post.sar import uniform_sphere_sar_closed_form
from fem_em_solver.thermal import solve_steady_conduction, solve_transient_conduction
from fem_em_solver.utils.analytical import robin_sphere_transient_temperature
from fem_em_solver.utils.constants import EPSILON_0

from tests.complex_mode import complex_only
from tests.validation.test_lossy_sphere_sar import (
    E0,
    EPSILON_R_SPHERE,
    OMEGA,
    RHO_KG_M3,
    SIGMA_HIGH,
    _solve_lossy_sphere,
)
from tests.validation.test_mat10_sar_to_heat_sphere import (
    A_BAND,
    G,
    H,
    K,
    P_BAND,
    QDEG,
    R,
    ROBIN_TAG,
    _allsum,
    _cell_average_dg0,
    _interior_indicator,
    _interior_rel,
    _robin_tags,
    _say,
    _sphere_submesh,
    _to_sub,
)

C_P = 4180.0  # J/kg/K
RHO_C = RHO_KG_M3 * C_P  # 4.18e6 J/m³/K, MAT-9's
T_R = R**2 * RHO_C / K  # 696.67 s
T_STAR = 0.2 * T_R  # Fo 0.2
N_MAIN = 40
DT_1 = 1.0e-5 * T_R  # 6.967 ms, the calorimetric step

S1_BAND = 1.0e-3
S2_BAND = 1.0e-4
E_BAND = 1.0e-9
INF_BAND = 1.0e-8


def _series(r, t, q_cf):
    return robin_sphere_transient_temperature(
        r, t, radius=R, k=K, rho_c=RHO_C, q=q_cf, h=H, t_inf=0.0, t_init=0.0,
    )


def _exact_theta(sub, t, q_cf):
    """θ(t) from the series, interpolated in numpy into P2 (chunked)."""
    V2 = fem.functionspace(sub, ("Lagrange", 2))
    f = fem.Function(V2)

    def ev(x):
        r = np.sqrt(np.sum(x**2, axis=0))
        out = np.empty(r.shape)
        for s in range(0, r.size, 20000):
            out[s:s + 20000] = _series(r[s:s + 20000], t, q_cf)
        return out

    f.interpolate(ev)
    f.x.scatter_forward()
    return f


def _whole_rel(sub, diff, ref):
    """sqrt(∫|diff|² / ∫|ref|²) over the whole ball, complex-safe."""
    dx = ufl.dx(domain=sub, metadata=QDEG)
    e = _allsum(sub.comm, ufl.inner(diff, diff) * dx)
    n = _allsum(sub.comm, ufl.inner(ref, ref) * dx)
    return float(np.sqrt(np.real(e) / np.real(n)))


def _mean(sub, f, vol):
    return float(np.real(_allsum(sub.comm, f * ufl.dx(domain=sub)))) / vol


def _cn(sub, ft, q, snapshot_steps=()):
    """CN, N = 40 to t*, two BE half-steps first (`MAT-9` (i)'s scheme).

    ``snapshot_steps`` passes through to ``solve_transient_conduction``
    (additive, `EX-67`'s frames; the gate never sets it).
    """
    return solve_transient_conduction(
        sub, ft, k=K, rho_c=RHO_C, q=q, h=H, t_inf=0.0, t_init=0.0,
        robin_tag=ROBIN_TAG, dt=T_STAR / N_MAIN, n_steps=N_MAIN, theta=0.5,
        startup_half_steps=2, degree=1, snapshot_steps=snapshot_steps,
    )


def _be_one(sub, ft, q):
    """One backward-Euler step of Δt₁ from θ = 0."""
    return solve_transient_conduction(
        sub, ft, k=K, rho_c=RHO_C, q=q, h=H, t_inf=0.0, t_init=0.0,
        robin_tag=ROBIN_TAG, dt=DT_1, n_steps=1, theta=1.0, degree=1,
    )


@complex_only
def test_mat12_transient_sar_heating_of_the_robin_ball():
    comm = MPI.COMM_WORLD
    t0 = time.perf_counter()

    run = _solve_lossy_sphere(SIGMA_HIGH, R / 10.0, R / 5.0, return_fields=True)
    msh, cell_tags = run["mesh"], run["cell_tags"]
    E, sigma = run["e_complex"], run["sigma_field"]
    t_em = time.perf_counter() - t0

    cf = uniform_sphere_sar_closed_form(
        e0=E0, epsilon_r=EPSILON_R_SPHERE, sigma=SIGMA_HIGH, omega=OMEGA,
        rho=RHO_KG_M3, epsilon_0=EPSILON_0,
    )
    q_cf = G * RHO_KG_M3 * cf["sar_w_per_kg"]
    th0 = float(_series(np.array([0.0]), T_STAR, q_cf)[0])
    thR = float(_series(np.array([R]), T_STAR, q_cf)[0])
    _say(f"MAT12 q_cf = {q_cf:.6f} W/m^3  t_R = {T_R:.4f} s  t* = {T_STAR:.4f} s  "
         f"dt = {T_STAR / N_MAIN:.4f} s  dt_1 = {DT_1 * 1e3:.4f} ms")
    _say(f"MAT12 series: theta(0, t*) = {th0:.6f} K  theta(R, t*) = {thR:.6f} K  "
         f"(predicted 0.059564 / 0.035599); initial slope q_cf/(rho c) = "
         f"{q_cf / RHO_C:.4e} K/s (predicted 4.7897e-4)")

    sub, parent_ids = _sphere_submesh(msh, cell_tags)
    ft = _robin_tags(sub)
    n_sub = sub.topology.index_map(3).size_global
    _say(f"MAT12 cells: sub = {n_sub}  (MAT-10's 19496)")

    q_expr = G * 0.5 * sigma * ufl.inner(E, E)
    Wp, q_par = _cell_average_dg0(msh, cell_tags, q_expr)
    q_h = _to_sub(Wp, q_par, sub, parent_ids)
    q_re_expr = G * 0.5 * sigma * ufl.inner(ufl.real(E), ufl.real(E))
    _, q_re_par = _cell_average_dg0(msh, cell_tags, q_re_expr)
    q_re = _to_sub(Wp, q_re_par, sub, parent_ids)

    one = fem.Constant(sub, default_scalar_type(1.0))
    vol = float(np.real(_allsum(comm, one * ufl.dx(domain=sub))))
    qbar = _mean(sub, q_h, vol)
    _say(f"MAT12 volume = {vol:.6e} m^3  mean q_h = {qbar:.6f} W/m^3  "
         f"(q_h/q_cf - 1 = {qbar / q_cf - 1.0:+.4e}); mean SAR/c = "
         f"{qbar / RHO_C:.6e} K/s")

    chi = _interior_indicator(sub)
    theta = _exact_theta(sub, T_STAR, q_cf)

    # (P_t) uniform q_cf — the first complex-build transient.
    t1 = time.perf_counter()
    res_P = _cn(sub, ft, q_cf)
    t_cn = time.perf_counter() - t1
    err_P = _interior_rel(sub, chi, res_P.temperature - theta, theta)
    _say(f"MAT12 (P_t) uniform q_cf, CN+Rannacher N={N_MAIN}: interior rel L2 = "
         f"{err_P:.4e}  (bound {P_BAND:.0e}, predicted ≈ 1.7e-3)  max_balance = "
         f"{res_P.max_balance_residual:.3e}  factorizations = {res_P.n_factorizations}  "
         f"steps = {len(res_P.balance_residuals)}  solve {t_cn:.1f} s")

    # (A_t) FEM SAR as the source.
    res_A = _cn(sub, ft, q_h)
    err_A = _interior_rel(sub, chi, res_A.temperature - theta, theta)
    leg_offset = _interior_rel(sub, chi, res_A.temperature - res_P.temperature, theta)
    _say(f"MAT12 (A_t) theta_h[q_h](t*) interior rel L2 = {err_A:.4e}  (bound {A_BAND}, "
         f"predicted ≈ 3.3–3.8e-2)  max_balance = {res_A.max_balance_residual:.3e}")
    _say(f"MAT12 (A_t) legs: offset theta_h[q_h - q_cf]/theta = {leg_offset:.4e}  "
         f"(mean-SAR offset {run['mean_sar'] / cf['sar_w_per_kg'] - 1.0:+.4e}); "
         f"(P_t) leg = {err_P:.4e}")

    # (S1) one BE step, uniform q_cf: the interior slope is q_cf/(ρc).
    res_S1 = _be_one(sub, ft, q_cf)
    qcf_c = fem.Constant(sub, default_scalar_type(q_cf))
    slope_S1 = RHO_C * res_S1.temperature / DT_1
    err_S1 = _interior_rel(sub, chi, slope_S1 - qcf_c, qcf_c)
    n_own = res_S1.temperature.function_space.dofmap.index_map.size_local
    th1_max = comm.allreduce(
        float(np.max(np.real(res_S1.temperature.x.array[:n_own]))), op=MPI.MAX
    )
    _say(f"MAT12 (S1) uniform q_cf, one BE step: interior ||rho c theta/dt_1 - q_cf||/||q_cf|| "
         f"= {err_S1:.4e}  (bound {S1_BAND:.0e}, predicted ≈ 1e-4)  max theta^1 = "
         f"{th1_max:.4e} K  balance = {res_S1.max_balance_residual:.3e}")

    # (S2) one BE step, q_h: the volume-mean slope is mean SAR / c.
    res_S2 = _be_one(sub, ft, q_h)
    thbar = _mean(sub, res_S2.temperature, vol)
    ratio_S2 = RHO_C * thbar / (DT_1 * qbar)
    ds = ufl.Measure("ds", domain=sub, subdomain_data=ft)(ROBIN_TAG)
    leak = float(np.real(_allsum(comm, H * res_S2.temperature * ds))) / (qbar * vol)
    _say(f"MAT12 (S2) q_h, one BE step: rho c mean(theta^1)/(dt_1 mean q_h) = "
         f"{ratio_S2:.9f}  |ratio - 1| = {abs(ratio_S2 - 1.0):.4e}  (bound {S2_BAND:.0e}); "
         f"surface deficit oint h theta^1 / int q_h = {leak:.4e}  "
         f"(uniform-source estimate 3 Bi dt_1/t_R = {3.0 * (H * R / K) * DT_1 / T_R:.1e})  "
         f"mean slope = {thbar / DT_1:.6e} K/s  balance = {res_S2.max_balance_residual:.3e}")
    qbar_c = fem.Constant(sub, default_scalar_type(qbar))
    cell_S2 = _whole_rel(sub, RHO_C * res_S2.temperature / DT_1 - q_h, q_h)
    spread = _whole_rel(sub, q_h - qbar_c, q_h)
    _say(f"MAT12 (S2) cell-wise, printed only: ||rho c theta^1/dt_1 - q_h||/||q_h|| = "
         f"{cell_S2:.4e}  beside 2||q_h - mean||/||q_h|| + 1e-3 = {2.0 * spread + 1e-3:.4e}")

    # (∞) long backward Euler against the steady solve with the same q_h.
    res_inf = solve_transient_conduction(
        sub, ft, k=K, rho_c=RHO_C, q=q_h, h=H, t_inf=0.0, t_init=0.0,
        robin_tag=ROBIN_TAG, dt=5.0 * T_R, n_steps=20, theta=1.0, degree=1,
    )
    st = solve_steady_conduction(sub, ft, k=K, q=q_h, h=H, t_inf=0.0,
                                 robin_tag=ROBIN_TAG, degree=1)
    err_inf = _whole_rel(sub, res_inf.temperature - st.temperature, st.temperature)
    _say(f"MAT12 (inf) BE 20 x 5 t_R vs steady[q_h]: rel L2 = {err_inf:.3e}  "
         f"(bound {INF_BAND:.0e})  max_balance = {res_inf.max_balance_residual:.3e}")

    # Negative controls.
    res_a = _cn(sub, ft, q_re)
    err_a = _interior_rel(sub, chi, res_a.temperature - theta, theta)
    res_a1 = _be_one(sub, ft, q_re)
    ratio_a = RHO_C * _mean(sub, res_a1.temperature, vol) / (DT_1 * qbar)
    _say(f"MAT12 control (a) Re-E-only: interior rel L2 = {err_a:.4e}  (predicted ≈ 0.791); "
         f"its (S2) ratio against mean q_h = {ratio_a:.6f}  (printed)")
    err_b = _interior_rel(sub, chi, st.temperature - theta, theta)
    _say(f"MAT12 control (b) quasi-static steady[q_h] for theta(t*): interior rel L2 = "
         f"{err_b:.4e}  (predicted ≈ 1.01)")

    balances = {
        "P_t": res_P.max_balance_residual, "A_t": res_A.max_balance_residual,
        "S1": res_S1.max_balance_residual, "S2": res_S2.max_balance_residual,
        "inf": res_inf.max_balance_residual, "ctl_a": res_a.max_balance_residual,
        "ctl_a_S2": res_a1.max_balance_residual,
    }
    worst = max(balances, key=balances.get)
    n_steps = sum(len(r.balance_residuals)
                  for r in (res_P, res_A, res_S1, res_S2, res_inf, res_a, res_a1))
    _say(f"MAT12 (E) worst per-step balance over {n_steps} steps = {balances[worst]:.3e} "
         f"on run {worst}  (bound {E_BAND:.0e})")
    _say(f"MAT12 elapsed: EM {t_em:.1f} s, total {time.perf_counter() - t0:.1f} s")

    assert balances[worst] <= E_BAND, f"(E) {worst} balance {balances[worst]:.3e} > 1e-9"
    assert err_P <= P_BAND, f"(P_t) interior {err_P:.4e} > 1 %"
    assert err_S1 <= S1_BAND, f"(S1) slope {err_S1:.4e} > 1e-3"
    assert abs(ratio_S2 - 1.0) <= S2_BAND, f"(S2) |ratio - 1| {abs(ratio_S2 - 1.0):.4e} > 1e-4"
    assert err_inf <= INF_BAND, f"(inf) {err_inf:.3e} > 1e-8"
    assert err_A <= A_BAND, f"(A_t) {err_A:.4e} > 11 % (composed band)"
    assert err_a > A_BAND, f"control (a) Re-E-only {err_a:.4e} did not miss 11 %"
    assert err_b > A_BAND, f"control (b) quasi-static {err_b:.4e} did not miss 11 %"
