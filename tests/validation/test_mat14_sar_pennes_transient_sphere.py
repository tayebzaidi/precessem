"""`MAT-14` — SAR-driven transient Pennes on the lossy sphere (§9 item 121).

`MAT-10`'s FEM SAR (q_h = DG0 cell average of G·½σ|E|² on the tag-1
submesh of `MAT-4`'s lossy sphere, R = 0.01 m) through `MAT-13`'s perfused
θ-scheme: k 0.6, h 100 (Bi 5/3), ρc 4.18e6, **w = 24 000 W/m³/K so mR = 2**
— `MAT-13`'s exact dimensionless twin (a numerical choice, not a tissue
claim) — arterial = ambient, **T∞ = T(0) = 20** in every main run (at
T∞ = 0 a missing ∫wT∞v is invisible).

Closed form: `robin_sphere_perfusion_transient_temperature` with
q = q_cf = G·ρ·`uniform_sphere_sar_closed_form` (uniform: the quasi-static
interior field is, `post/sar.py:377–395`).

Anchors: (N-sim, numpy) θ/(qR²/k) equals `MAT-13`'s at ρ ∈ {0, ½, 1},
Fo ∈ {0.05, 0.2, 1} to ≤ 1e-12; (P_w′) uniform q_cf, CN + Rannacher N 40
to t*: interior rel L2 on θ ≤ `P_BAND`; (A_w) q_h ≤ `A_BAND` (composed);
(P_ss) steady perfused q_cf vs `robin_sphere_perfusion_temperature`
≤ `P_BAND`; (E) per-step balance ≤ `E_BAND` on every run, recomputed on
(A_w) from snapshots 39 / 40; (∞) BE 20 × 5 t_R with q_h and w equals
``solve_steady_conduction(q=q_h, perfusion_w_per_m3_k=W)`` ≤ `INF_BAND`;
(S1) one BE step at T∞ = 0 with w ≤ `S1_BAND`; (S2) printed only;
(w→0) `MAT-12`'s own `_cn` legs reproduce its 1.6845e-03 / 3.4708e-02.
Negative controls (asserted only as band misses): (a) w dropped, q_cf;
(a′) w dropped, q_h; (b) Re-E-only with w; (c) perfused term dropped
from (E). Quasi-static steady field for θ(t*): printed.
"""

from __future__ import annotations

import time

import numpy as np
import ufl
from dolfinx import default_scalar_type, fem
from mpi4py import MPI

from fem_em_solver.post.sar import uniform_sphere_sar_closed_form
from fem_em_solver.thermal import solve_steady_conduction, solve_transient_conduction
from fem_em_solver.utils.analytical import (
    robin_sphere_perfusion_temperature,
    robin_sphere_perfusion_transient_temperature,
)
from fem_em_solver.utils.constants import EPSILON_0
from tests.complex_mode import complex_only
from tests.validation import test_mat13_transient_perfusion_sphere as m13
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
    P_BAND,
    ROBIN_TAG,
    G,
    H,
    K,
    R,
    _allsum,
    _cell_average_dg0,
    _interior_indicator,
    _interior_rel,
    _robin_tags,
    _say,
    _sphere_submesh,
    _to_sub,
)
from tests.validation.test_mat12_transient_sar_heating_sphere import (
    DT_1,
    E_BAND,
    INF_BAND,
    N_MAIN,
    RHO_C,
    S1_BAND,
    T_R,
    T_STAR,
    _cn,
    _exact_theta,
    _whole_rel,
)

W = 24000.0  # W/m³/K: mR = 2, MAT-13's twin
T_INF = 20.0
# (w→0) `MAT-12`'s own printed digits (20260930T110230Z_MAT-12-W1.log:161–162).
MAT12_P_RECORD = "1.6845e-03"
MAT12_A_RECORD = "3.4708e-02"


def _series_w(r, t, q, w=W, n_terms=400):
    return robin_sphere_perfusion_transient_temperature(
        r,
        t,
        radius=R,
        k=K,
        rho_c=RHO_C,
        q=q,
        w=w,
        h=H,
        t_inf=T_INF,
        t_init=T_INF,
        n_terms=n_terms,
    )


def _cn_w(sub, ft, q, w=W, snapshot_steps=()):
    """CN, N = 40 to t*, Rannacher start, T∞ = T(0) = 20, perfusion w."""
    return solve_transient_conduction(
        sub,
        ft,
        k=K,
        rho_c=RHO_C,
        q=q,
        h=H,
        t_inf=T_INF,
        t_init=T_INF,
        robin_tag=ROBIN_TAG,
        dt=T_STAR / N_MAIN,
        n_steps=N_MAIN,
        theta=0.5,
        startup_half_steps=2,
        degree=1,
        snapshot_steps=snapshot_steps,
        perfusion_w_per_m3_k=w,
    )


def _exact_theta_w(sub, t, q, w=W):
    """θ = T − T∞ from the perfused series, numpy into P2 (chunked)."""
    V2 = fem.functionspace(sub, ("Lagrange", 2))
    f = fem.Function(V2)

    def ev(x):
        r = np.sqrt(np.sum(x**2, axis=0))
        out = np.empty(r.shape)
        for s in range(0, r.size, 20000):
            out[s : s + 20000] = _series_w(r[s : s + 20000], t, q, w=w) - T_INF
        return out

    f.interpolate(ev)
    f.x.scatter_forward()
    return f


def _q_cf():
    cf = uniform_sphere_sar_closed_form(
        e0=E0,
        epsilon_r=EPSILON_R_SPHERE,
        sigma=SIGMA_HIGH,
        omega=OMEGA,
        rho=RHO_KG_M3,
        epsilon_0=EPSILON_0,
    )
    return G * RHO_KG_M3 * cf["sar_w_per_kg"], cf


# ------------------------------------------------------------------ numpy first


def test_nsim_scaled_series_is_mat13s():
    q_cf, _ = _q_cf()
    alpha = K / RHO_C
    alpha13 = m13.K / m13.RHO_C
    assert abs(H * R / K - m13.H * m13.R / m13.K) < 1e-15, "Bi differs"
    assert abs(W * R**2 / K - m13.W * m13.R**2 / m13.K) < 1e-12, "(mR)^2 differs"
    worst = 0.0
    for fo in (0.05, 0.2, 1.0):
        rho = np.array([0.0, 0.5, 1.0])
        a = (_series_w(rho * R, fo * R**2 / alpha, q_cf) - T_INF) / (q_cf * R**2 / K)
        b = (m13._series_w(rho * m13.R, fo * m13.R**2 / alpha13) - m13.T_INF) / (
            m13.Q * m13.R**2 / m13.K
        )
        worst = max(worst, float(np.max(np.abs(a - b))))
    th = _series_w(np.array([0.0, R]), T_STAR, q_cf) - T_INF
    _say(
        f"MAT14 (N-sim) max |theta/(qR^2/k) - MAT-13's| over rho in (0, 1/2, 1), "
        f"Fo in (0.05, 0.2, 1) = {worst:.3e} (bar 1e-12); theta(0, t*) = {th[0]:.6f} K, "
        f"theta(R, t*) = {th[1]:.6f} K (predicted 0.042031 / 0.026036)"
    )
    assert worst <= 1e-12, f"(N-sim) {worst:.3e}"


# ------------------------------------------------------------------ FEM anchors


@complex_only
def test_mat14_sar_driven_pennes_transient():
    comm = MPI.COMM_WORLD
    t0 = time.perf_counter()

    run = _solve_lossy_sphere(SIGMA_HIGH, R / 10.0, R / 5.0, return_fields=True)
    msh, cell_tags = run["mesh"], run["cell_tags"]
    E, sigma = run["e_complex"], run["sigma_field"]
    t_em = time.perf_counter() - t0

    q_cf, cf = _q_cf()
    _say(
        f"MAT14 q_cf = {q_cf:.6f} W/m^3 (MAT-12's 2002.088514)  w = {W}  "
        f"mR = {np.sqrt(W / K) * R:.6f}  Bi = {H * R / K:.6f}  t* = {T_STAR:.4f} s"
    )

    sub, parent_ids = _sphere_submesh(msh, cell_tags)
    ft = _robin_tags(sub)
    _say(f"MAT14 cells: sub = {sub.topology.index_map(3).size_global}  (MAT-10's 19496)")

    q_expr = G * 0.5 * sigma * ufl.inner(E, E)
    Wp, q_par = _cell_average_dg0(msh, cell_tags, q_expr)
    q_h = _to_sub(Wp, q_par, sub, parent_ids)
    q_re_expr = G * 0.5 * sigma * ufl.inner(ufl.real(E), ufl.real(E))
    _, q_re_par = _cell_average_dg0(msh, cell_tags, q_re_expr)
    q_re = _to_sub(Wp, q_re_par, sub, parent_ids)

    chi = _interior_indicator(sub)
    theta_w = _exact_theta_w(sub, T_STAR, q_cf)

    # (P_w′) uniform q_cf.
    t1 = time.perf_counter()
    res_P = _cn_w(sub, ft, q_cf)
    t_cn = time.perf_counter() - t1
    err_P = _interior_rel(sub, chi, res_P.temperature - T_INF - theta_w, theta_w)
    _say(
        f"MAT14 (P_w') uniform q_cf, CN+Rannacher N={N_MAIN}: interior rel L2 = "
        f"{err_P:.4e}  (bound {P_BAND:.0e}, predicted ≈ 1.64e-3)  max_balance = "
        f"{res_P.max_balance_residual:.3e}  perfused_power = {res_P.perfused_power:.4e} W  "
        f"solve {t_cn:.1f} s"
    )

    # (A_w) FEM SAR as the source, snapshots for (E).
    res_A = _cn_w(sub, ft, q_h, snapshot_steps=(39, 40))
    err_A = _interior_rel(sub, chi, res_A.temperature - T_INF - theta_w, theta_w)
    leg = _interior_rel(sub, chi, res_A.temperature - res_P.temperature, theta_w)
    _say(
        f"MAT14 (A_w) theta_h[q_h](t*) interior rel L2 = {err_A:.4e}  (bound {A_BAND}, "
        f"predicted ≈ 3.4–3.8e-2)  max_balance = {res_A.max_balance_residual:.3e}"
    )
    _say(
        f"MAT14 (A_w) offset leg theta_h[q_h - q_cf]/theta = {leg:.4e}  (mean-SAR offset "
        f"{run['mean_sar'] / cf['sar_w_per_kg'] - 1.0:+.4e}, plan +3.5364e-2)"
    )

    # (E) recomputed from snapshots 39 / 40 (CN step) and control (c).
    T39, T40 = res_A.snapshots[39], res_A.snapshots[40]
    ds = ufl.Measure("ds", domain=sub, subdomain_data=ft)
    dx = ufl.dx(domain=sub)
    dt = T_STAR / N_MAIN
    stor = float(np.real(_allsum(comm, RHO_C * (T40 - T39) / dt * dx)))
    conv = float(np.real(_allsum(comm, H * (0.5 * T40 + 0.5 * T39 - T_INF) * ds(ROBIN_TAG))))
    perf = float(np.real(_allsum(comm, W * (0.5 * T40 + 0.5 * T39 - T_INF) * dx)))
    src = res_A.source_power
    bal = abs(stor + conv + perf - src) / src
    bal_c = abs(stor + conv - src) / src
    _say(
        f"MAT14 (E) (A_w) recomputed last-step balance = {bal:.3e} (bound {E_BAND:.0e}); "
        f"split storage/convected/perfused = {stor / src:.4f}/{conv / src:.4f}/"
        f"{perf / src:.4f} (predicted 0.21/0.39/0.40); perf vs result {perf:.6e}/"
        f"{res_A.perfused_power:.6e}"
    )
    _say(
        f"MAT14 control (c) perfused term dropped from (E): {bal_c:.4e} (asserted > 1e-3; "
        f"predicted ≈ 0.398)"
    )

    # (P_ss) steady perfused, q_cf.
    st_cf = solve_steady_conduction(
        sub,
        ft,
        k=K,
        q=q_cf,
        h=H,
        t_inf=T_INF,
        robin_tag=ROBIN_TAG,
        degree=1,
        perfusion_w_per_m3_k=W,
    )
    V2 = fem.functionspace(sub, ("Lagrange", 2))
    th_ss = fem.Function(V2)
    th_ss.interpolate(
        lambda x: robin_sphere_perfusion_temperature(
            np.sqrt(np.sum(x**2, axis=0)), radius=R, k=K, q=q_cf, w=W, h=H
        )
    )
    th_ss.x.scatter_forward()
    err_ss = _interior_rel(sub, chi, st_cf.temperature - T_INF - th_ss, th_ss)
    _say(
        f"MAT14 (P_ss) steady perfused q_cf vs closed form: interior rel L2 = {err_ss:.4e} "
        f"(bound {P_BAND:.0e}, predicted ≈ 1.8e-3)"
    )

    # (∞) long BE with q_h and w vs the perfused steady solve with q_h.
    res_inf = solve_transient_conduction(
        sub,
        ft,
        k=K,
        rho_c=RHO_C,
        q=q_h,
        h=H,
        t_inf=T_INF,
        t_init=T_INF,
        robin_tag=ROBIN_TAG,
        dt=5.0 * T_R,
        n_steps=20,
        theta=1.0,
        degree=1,
        perfusion_w_per_m3_k=W,
    )
    st_h = solve_steady_conduction(
        sub, ft, k=K, q=q_h, h=H, t_inf=T_INF, robin_tag=ROBIN_TAG, degree=1, perfusion_w_per_m3_k=W
    )
    err_inf = _whole_rel(sub, res_inf.temperature - st_h.temperature, st_h.temperature - T_INF)
    _say(
        f"MAT14 (inf) BE 20 x 5 t_R vs steady[q_h, w]: rel L2 on theta = {err_inf:.3e} "
        f"(bound {INF_BAND:.0e})  max_balance = {res_inf.max_balance_residual:.3e}; "
        f"perfused power transient/steady = {res_inf.perfused_power:.8e}/"
        f"{st_h.perfused_power:.8e}"
    )

    # (S1) one BE step at T∞ = 0 with w, q_cf; (S2) q_h, printed only.
    res_S1 = solve_transient_conduction(
        sub,
        ft,
        k=K,
        rho_c=RHO_C,
        q=q_cf,
        h=H,
        t_inf=0.0,
        t_init=0.0,
        robin_tag=ROBIN_TAG,
        dt=DT_1,
        n_steps=1,
        theta=1.0,
        degree=1,
        perfusion_w_per_m3_k=W,
    )
    qcf_c = fem.Constant(sub, default_scalar_type(q_cf))
    err_S1 = _interior_rel(sub, chi, RHO_C * res_S1.temperature / DT_1 - qcf_c, qcf_c)
    res_S2 = solve_transient_conduction(
        sub,
        ft,
        k=K,
        rho_c=RHO_C,
        q=q_h,
        h=H,
        t_inf=0.0,
        t_init=0.0,
        robin_tag=ROBIN_TAG,
        dt=DT_1,
        n_steps=1,
        theta=1.0,
        degree=1,
        perfusion_w_per_m3_k=W,
    )
    one = fem.Constant(sub, default_scalar_type(1.0))
    vol = float(np.real(_allsum(comm, one * dx)))
    qbar = float(np.real(_allsum(comm, q_h * dx))) / vol
    thbar = float(np.real(_allsum(comm, res_S2.temperature * dx))) / vol
    ratio_S2 = RHO_C * thbar / (DT_1 * qbar)
    _say(
        f"MAT14 (S1) one BE step at T_inf = 0 with w: interior ||rho c theta/dt_1 - q_cf||"
        f"/||q_cf|| = {err_S1:.4e} (bound {S1_BAND:.0e}, predicted ≈ 5e-5)  balance = "
        f"{res_S1.max_balance_residual:.3e}"
    )
    _say(
        f"MAT14 (S2) [printed only] q_h with w: |rho c mean(theta^1)/(dt_1 mean q_h) - 1| = "
        f"{abs(ratio_S2 - 1.0):.4e} (predicted ≈ 9.0e-5, 90 % of 1e-4)  balance = "
        f"{res_S2.max_balance_residual:.3e}"
    )

    # (w→0) MAT-12's own legs against its own series.
    theta_0 = _exact_theta(sub, T_STAR, q_cf)
    r12P = _cn(sub, ft, q_cf)
    r12A = _cn(sub, ft, q_h)
    e12P = _interior_rel(sub, chi, r12P.temperature - theta_0, theta_0)
    e12A = _interior_rel(sub, chi, r12A.temperature - theta_0, theta_0)
    _say(
        f"MAT14 (w->0) MAT-12 _cn legs: q_cf {e12P:.4e} (record {MAT12_P_RECORD}), "
        f"q_h {e12A:.4e} (record {MAT12_A_RECORD})"
    )

    # Negative controls.
    r_a = _cn_w(sub, ft, q_cf, w=0.0)
    e_a = _interior_rel(sub, chi, r_a.temperature - T_INF - theta_w, theta_w)
    r_a1 = _cn_w(sub, ft, q_h, w=0.0)
    e_a1 = _interior_rel(sub, chi, r_a1.temperature - T_INF - theta_w, theta_w)
    r_b = _cn_w(sub, ft, q_re)
    e_b = _interior_rel(sub, chi, r_b.temperature - T_INF - theta_w, theta_w)
    e_qs = _interior_rel(sub, chi, st_h.temperature - T_INF - theta_w, theta_w)
    _say(
        f"MAT14 control (a) w dropped, q_cf: {e_a:.4e} (asserted > {P_BAND}; predicted "
        f"≈ 0.389); (a') w dropped, q_h: {e_a1:.4e} (asserted > {A_BAND}; predicted ≈ "
        f"0.43–0.45); (b) Re-E-only with w: {e_b:.4e} (asserted > {A_BAND}; predicted "
        f"≈ 0.79–0.80)"
    )
    _say(
        f"MAT14 [printed] quasi-static steady[q_h, w] for theta(t*): {e_qs:.4e} "
        f"(predicted ≈ 0.28–0.32)"
    )

    balances = {
        "P_w'": res_P.max_balance_residual,
        "A_w": res_A.max_balance_residual,
        "inf": res_inf.max_balance_residual,
        "S1": res_S1.max_balance_residual,
        "S2": res_S2.max_balance_residual,
        "w0_P": r12P.max_balance_residual,
        "w0_A": r12A.max_balance_residual,
        "ctl_a": r_a.max_balance_residual,
        "ctl_a'": r_a1.max_balance_residual,
        "ctl_b": r_b.max_balance_residual,
    }
    worst = max(balances, key=balances.get)
    _say(
        f"MAT14 (E) worst per-step balance over {len(balances)} runs = "
        f"{balances[worst]:.3e} on run {worst} (bound {E_BAND:.0e})"
    )
    _say(f"MAT14 elapsed: EM {t_em:.1f} s, total {time.perf_counter() - t0:.1f} s")

    assert balances[worst] <= E_BAND, f"(E) {worst} {balances[worst]:.3e}"
    assert bal <= E_BAND, f"(E) recomputed {bal:.3e}"
    assert err_P <= P_BAND, f"(P_w') {err_P:.4e} > 1 %"
    assert err_ss <= P_BAND, f"(P_ss) {err_ss:.4e} > 1 %"
    assert err_inf <= INF_BAND, f"(inf) {err_inf:.3e}"
    assert err_S1 <= S1_BAND, f"(S1) {err_S1:.4e}"
    assert err_A <= A_BAND, f"(A_w) {err_A:.4e} > 11 % (composed band)"
    assert f"{e12P:.4e}" == MAT12_P_RECORD, f"(w->0) q_cf {e12P:.4e}"
    assert f"{e12A:.4e}" == MAT12_A_RECORD, f"(w->0) q_h {e12A:.4e}"
    assert e_a > P_BAND, f"control (a) {e_a:.4e} did not miss 1 %"
    assert e_a1 > A_BAND, f"control (a') {e_a1:.4e} did not miss 11 %"
    assert e_b > A_BAND, f"control (b) {e_b:.4e} did not miss 11 %"
    assert bal_c > 1e-3, f"control (c) {bal_c:.4e}"
