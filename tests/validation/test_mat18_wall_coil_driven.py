"""`MAT-18` step 3 — the container wall on the coil-driven phantom (§9 item 270, chain P11).

Fixture: `MAT-15` step 1's (`_run_heated(return_fields=True)`, imported, not
edited): `ANS-2`'s F-small at 10 MHz, the 0.0025 m phantom rung, the four
single-leg drives, G = 1 W / P₁, P1 thermal on the tag-3 submesh.

Wall model: the submesh has no wall cells, so `MAT-18`'s shell (t_w = b − a,
k₂, outer Robin h) enters as its **series thermal resistance** in the Robin
coefficient — the steady limit of the conforming shell, no wall heat capacity:
h_eff = 1 / (1/h + t_w/k₂) (planar), one scalar for the lateral faces and the
caps.  The lateral cylindrical value is printed beside it and never used.

Asserted: (R) on the sphere the series reduction is exact in steady state —
`two_layer_sphere_steady_temperature` against `robin_sphere_steady_temperature`
with h_sph, ≤ 1e-12; (E′) ρc∫Δθ/Δt + ∮h_eff θ^(n+½) = G·P_k from the EM side
≤ ``E_BAND`` and (E) every wall run's balance ≤ ``E_BAND``; (L0) t_w = 0
reproduces P1's volume-mean θ(t*) per drive ≤ ``B_BAND``; (C4) the 10 g ball
means of θ at t* with the wall, four cyclic pairs ≤ ``C4_COVARIANCE_BAND``.
Negative controls: (d) doubled source misses (E′) by > ``D_CONTROL_FLOOR``
(asserted, prior same-fixture measurement 1.000000); (W-ctl) the wall is felt
— predicted +2 … +6 % on ⟨θ⟩(t*), asserted only as a miss > 1e-6.

Scope: series-resistance wall, P1's untuned 10 MHz σ = 800 S/m fixture, P1
thermal, gel constants, no perfusion.  Not: a conforming wall, wall storage,
contact resistance, a real container, tuned / copper / 64 MHz, absolute
temperature, `ANS-9`.  Probe B is not printed (its centres live on `MAT-19`'s
fixture) — deferred.
"""

from __future__ import annotations

import time

import numpy as np
import pytest
import ufl
from dolfinx import default_scalar_type, fem
from mpi4py import MPI

from fem_em_solver.post.sar import averaging_ball_radius
from fem_em_solver.thermal import solve_transient_conduction
from fem_em_solver.utils.analytical import (
    robin_sphere_steady_temperature,
    two_layer_sphere_steady_temperature,
)

from tests.complex_mode import complex_only
from tests.validation.test_birdcage_b1_plus_map import (
    C4_COVARIANCE_BAND,
    PHANTOM_RHO_KG_PER_M3,
)
from tests.validation.test_mass_averaged_sar_standard_masses import (
    QUADRATURE_DEGREE,
    TEN_GRAM_KG,
)
from tests.validation.test_mat10_sar_to_heat_sphere import (
    B_BAND,
    ROBIN_TAG,
    _allsum,
    _say,
)
from tests.validation.test_mat12_transient_sar_heating_sphere import E_BAND, RHO_C
from tests.validation.test_mat15_coil_driven_heating_phantom import (
    D_CONTROL_FLOOR,
    H_GEL,
    K_GEL,
    N_MAIN,
    T_INF,
    T_STAR,
    _max_imag,
    _run_heated,
)
from tests.validation.test_mat15_coil_driven_heating_phantom import A as A_PHANTOM
from tests.validation.test_mat18_two_layer_sphere_closed_form import A_, B_, H, K1, K2, Q1
from tests.validation.test_mat18_two_layer_sphere_closed_form import T_INF as T_INF_18

T_WALL = B_ - A_  # 0.01 m
R_BAND = 1.0e-12  # (R): numpy round-off; MAT-18 step 1's steady band is 1e-13
W_CONTROL_FLOOR = 1.0e-6  # (W-ctl) predicted +2 … +6 %, asserted only as a miss


def _h_eff(h, t_w, k_w):
    """Planar series resistance of a wall of thickness t_w, conductivity k_w."""
    return 1.0 / (1.0 / h + t_w / k_w)


def _h_lateral(h, t_w, k_w, a):
    """Cylindrical series value referred to r = a (printed, never used)."""
    return 1.0 / (a * np.log(1.0 + t_w / a) / k_w + a / ((a + t_w) * h))


def _transient_h(sub, ft, q, h, *, dt, n_steps, theta, startup_half_steps=0, snapshot_steps=()):
    """`MAT-15` step 1's ``_transient`` with ``h`` as an argument."""
    return solve_transient_conduction(
        sub, ft, k=K_GEL, rho_c=RHO_C, q=q, h=h, t_inf=T_INF, t_init=T_INF,
        robin_tag=ROBIN_TAG, dt=dt, n_steps=n_steps, theta=theta,
        startup_half_steps=startup_half_steps, degree=1,
        snapshot_steps=snapshot_steps,
    )


def _cn_h(sub, ft, q, h):
    """CN, N = 40 to t*, two BE half-steps first; snapshots N−1, N (step 1's ``_cn``)."""
    return _transient_h(sub, ft, q, h, dt=T_STAR / N_MAIN, n_steps=N_MAIN, theta=0.5,
                        startup_half_steps=2, snapshot_steps=(N_MAIN - 1, N_MAIN))


def _e_prime_h(sub, ds, res, deposited, h):
    """ρc∫(θᴺ−θᴺ⁻¹)/Δt + ∮hθ^(n+½) against the EM-side deposited power."""
    comm = sub.comm
    tn, tm = res.snapshots[N_MAIN], res.snapshots[N_MAIN - 1]
    dt = T_STAR / N_MAIN
    stor = float(np.real(_allsum(comm, RHO_C * (tn - tm) / dt * ufl.dx(domain=sub))))
    conv = float(np.real(_allsum(comm, h * (0.5 * (tn + tm) - T_INF) * ds)))
    return abs(stor + conv - deposited) / abs(deposited), stor, conv


def _r_precondition():
    """(R): two-layer sphere steady state = one-layer Robin sphere with h_sph."""
    h_sph = 1.0 / (A_**2 * (1.0 / A_ - 1.0 / B_) / K2 + A_**2 / (B_**2 * H))
    r = np.linspace(0.0, A_, 20)
    two = two_layer_sphere_steady_temperature(r, a=A_, b=B_, k1=K1, k2=K2, q1=Q1, q2=0.0,
                                              h=H, t_inf=T_INF_18)
    one = robin_sphere_steady_temperature(r, radius=A_, k=K1, q=Q1, h=h_sph, t_inf=T_INF_18)
    rel = float(np.max(np.abs(two - one)) / np.max(np.abs(two - T_INF_18)))
    return rel, h_sph


@pytest.fixture(scope="module")
def walled():
    return _run_walled()


def _run_walled(*, return_fields=False):
    """The step-3 fixture body, lifted (EX-98, rule (a)) so a script can call it.

    ``return_fields=True`` adds ``base`` (the `_run_heated(return_fields=True)` dict) and
    ``wall_runs`` ({k: the h_eff TransientResult}); otherwise identical to the fixture.
    """
    comm = MPI.COMM_WORLD
    out = {}
    out["R"], h_sph = _r_precondition()
    _say(f"\nMAT18s3 (R) two-layer sphere steady vs Robin sphere with h_sph = {h_sph:.12f}: "
         f"rel {out['R']:.3e} (bound {R_BAND:.0e})")
    if out["R"] > R_BAND:  # stop before any solve (Negative result, (R) red)
        return out

    h_eff = _h_eff(H_GEL, T_WALL, K2)
    h_lat = _h_lateral(H_GEL, T_WALL, K2, A_PHANTOM)
    h_zero = _h_eff(H_GEL, 0.0, K2)
    out.update(h_eff=h_eff, h_lat=h_lat)
    _say(f"MAT18s3 wall t_w = {T_WALL:.4f} m, k2 = {K2}, h = {H_GEL}: h_eff (planar, used) = "
         f"{h_eff:.9f} W/m^2/K (predicted ≈ 6.552);  lateral cylindrical at a = {A_PHANTOM} m = "
         f"{h_lat:.9f} W/m^2/K (printed, not used);  t_w = 0 limit h = {h_zero!r}")

    base = _run_heated(return_fields=True)
    t0 = time.perf_counter()
    sub, ft, G, P = base["sub"], base["ft"], base["G"], base["P"]
    v_sub = base["v_sub"]
    ds = ufl.Measure("ds", domain=sub, subdomain_data=ft)(ROBIN_TAG)

    # Ball-mean forms (step 1's :231–261, centres from the fixture).
    a10 = averaging_ball_radius(mass_kg=TEN_GRAM_KG, rho=PHANTOM_RHO_KG_PER_M3)
    centre_c = fem.Constant(sub, np.zeros(3, dtype=default_scalar_type))
    xs = ufl.SpatialCoordinate(sub)
    off = xs - centre_c
    inside = ufl.conditional(ufl.lt(ufl.real(ufl.dot(off, off)), a10**2), 1.0, 0.0)
    dxq = ufl.dx(domain=sub, metadata={"quadrature_degree": int(QUADRATURE_DEGREE)})
    V1 = fem.functionspace(sub, ("Lagrange", 1))
    Ws = fem.functionspace(sub, ("DG", 0))
    th_ball = fem.Function(V1)
    q_ball = fem.Function(Ws)
    f_vol = fem.form(inside * dxq)
    f_th = fem.form(th_ball * inside * dxq)
    f_q = fem.form(q_ball * inside * dxq)

    def _red(form):
        return float(np.real(comm.allreduce(fem.assemble_scalar(form), op=MPI.SUM)))

    centres = base["centres"]
    ball_vol = {}
    for j in range(4):
        centre_c.value[:] = centres[j]
        ball_vol[j] = _red(f_vol)

    def _ball_mean(form, j):
        centre_c.value[:] = centres[j]
        return _red(form) / ball_vol[j]

    def _mean_theta(T):
        return float(np.real(_allsum(comm, (T - T_INF) * ufl.dx(domain=sub)))) / v_sub

    out.update(Ep={}, E={}, L0={}, W={}, ratio={}, th_own={}, q_own={}, th_p1={},
               mean_w={}, mean_p1={})
    tmax_imag = 0.0
    wall_runs = {}
    for k in range(4):
        q_h = base["q_h"][k]
        dep = G * P[k]
        mean_p1 = _mean_theta(base["T_star"][k])

        t1 = time.perf_counter()
        res = _cn_h(sub, ft, q_h, h_eff)
        wall_runs[k] = res
        t_cn = time.perf_counter() - t1
        ep, stor, conv = _e_prime_h(sub, ds, res, dep, h_eff)
        out["Ep"][k] = ep
        out["E"][f"wall{k}"] = res.max_balance_residual
        mean_w = _mean_theta(res.temperature)

        res0 = _cn_h(sub, ft, q_h, h_zero)
        mean_0 = _mean_theta(res0.temperature)
        out["L0"][k] = abs(mean_0 / mean_p1 - 1.0)
        out["W"][k] = abs(mean_w / mean_p1 - 1.0)
        out["ratio"][k] = mean_w / mean_p1
        out["mean_w"][k], out["mean_p1"][k] = mean_w, mean_p1

        th_ball.x.array[:] = res.temperature.x.array - T_INF
        th_ball.x.scatter_forward()
        out["th_own"][k] = _ball_mean(f_th, k)
        th_ball.x.array[:] = base["T_star"][k].x.array - T_INF
        th_ball.x.scatter_forward()
        out["th_p1"][k] = _ball_mean(f_th, k)
        q_ball.x.array[:] = q_h.x.array
        q_ball.x.scatter_forward()
        out["q_own"][k] = _ball_mean(f_q, k)
        tmax_imag = max(tmax_imag, _max_imag(res.temperature), _max_imag(res0.temperature))

        _say(f"MAT18s3 drive k={k} ({base['by_k'][k]}): wall CN N={N_MAIN} solve {t_cn:.1f} s;  "
             f"(E) max_balance = {res.max_balance_residual:.3e} (bound {E_BAND:.0e})")
        _say(f"MAT18s3   (E') stored {stor:.12e} W + surface {conv:.12e} W vs deposited G*P_k "
             f"{dep:.12e} W: rel {ep:.3e} (bound {E_BAND:.0e})")
        _say(f"MAT18s3   <theta>(t*): wall {mean_w:.9f} K  P1 {mean_p1:.9f} K  ratio {mean_w / mean_p1:.6f} "
             f"(wall effect {100 * (mean_w / mean_p1 - 1.0):+.4f} %, predicted +2 … +6 %; (W-ctl) "
             f"asserted > {W_CONTROL_FLOOR:.0e});  (L0) t_w=0 {mean_0:.12f} K rel {out['L0'][k]:.3e} "
             f"(bound {B_BAND:.0e})")
        _say(f"MAT18s3   10 g ball theta(c_k; k) at t*: wall {out['th_own'][k]:.9e} K  P1 "
             f"{out['th_p1'][k]:.9e} K  ratio {out['th_own'][k] / out['th_p1'][k]:.6f};  "
             f"q_h(c_k; k) = {out['q_own'][k]:.9e} W/m^3")

    # Negative control (d): doubled source on drive k = 0, with the wall.
    q0 = base["q_h"][0]
    q_two = fem.Function(q0.function_space)
    q_two.x.array[:] = 2.0 * q0.x.array
    q_two.x.scatter_forward()
    res_d = _cn_h(sub, ft, q_two, h_eff)
    out["ctl_d"], _, _ = _e_prime_h(sub, ds, res_d, G * P[0], h_eff)
    _say(f"MAT18s3 control (d) doubled source, k=0, wall: (E') vs G*P_1 = {out['ctl_d']:.6f} "
         f"(prior 1.000000, asserted > {D_CONTROL_FLOOR:.0e}); internal (E) = "
         f"{res_d.max_balance_residual:.3e}")

    out["pairs"] = {k: abs(out["th_own"][(k + 1) % 4] - out["th_own"][k]) / out["th_own"][k]
                    for k in range(4)}
    out["q_pairs"] = {k: abs(out["q_own"][(k + 1) % 4] - out["q_own"][k]) / out["q_own"][k]
                      for k in range(4)}
    for k in range(4):
        _say(f"MAT18s3 (C4) k={k}->{(k + 1) % 4}: wall theta 10 g pair {out['pairs'][k] * 100:.4f}% "
             f"(bound {C4_COVARIANCE_BAND * 100:.1f}%, predicted ≈ step 1's 0.020-0.049%);  "
             f"q_h pair {out['q_pairs'][k] * 100:.4f}% (printed)")
    kpk = max(out["th_own"], key=out["th_own"].get)
    _say(f"MAT18s3 peak 10 g theta with wall = {out['th_own'][kpk]:.9e} K on k={kpk}; P1 "
         f"{out['th_p1'][kpk]:.9e} K; ratio {out['th_own'][kpk] / out['th_p1'][kpk]:.6f}")
    worst = max(out["E"], key=out["E"].get)
    out["max_imag"] = tmax_imag
    _say(f"MAT18s3 (E) worst wall balance over {len(out['E'])} runs = {out['E'][worst]:.3e} on "
         f"{worst} (bound {E_BAND:.0e}); max|Im theta| = {tmax_imag:.3e} (predicted 0)")
    _say(f"MAT18s3 elapsed after P1 fixture: {time.perf_counter() - t0:.1f} s")
    if return_fields:
        out["base"] = base
        out["wall_runs"] = wall_runs
    return out


@complex_only
def test_r_series_resistance_is_exact_on_the_sphere(walled):
    assert walled["R"] <= R_BAND, f"(R) {walled['R']:.3e} > {R_BAND:.0e}"


@complex_only
@pytest.mark.parametrize("k", [0, 1, 2, 3])
def test_e_prime_wall_stored_plus_surface_equals_em_deposited(walled, k):
    assert walled["Ep"][k] <= E_BAND, f"(E') drive {k}: {walled['Ep'][k]:.3e} > {E_BAND:.0e}"


@complex_only
def test_e_every_wall_run_balances(walled):
    for name, r in walled["E"].items():
        assert r <= E_BAND, f"(E) run {name}: {r:.3e} > {E_BAND:.0e}"


@complex_only
@pytest.mark.parametrize("k", [0, 1, 2, 3])
def test_l0_zero_wall_reproduces_p1(walled, k):
    assert walled["L0"][k] <= B_BAND, f"(L0) drive {k}: {walled['L0'][k]:.3e} > {B_BAND:.0e}"


@complex_only
@pytest.mark.parametrize("k", [0, 1, 2, 3])
def test_c4_wall_ten_gram_theta_follows_the_drive(walled, k):
    r = walled["pairs"][k]
    assert r <= C4_COVARIANCE_BAND, f"(C4) pair {k}: {r * 100:.4f}% > 5 %"


@complex_only
def test_control_d_doubled_source_misses_the_balance(walled):
    assert walled["ctl_d"] > D_CONTROL_FLOOR, (
        f"control (d) doubled source reads (E') {walled['ctl_d']:.3e} — did not miss"
    )


@complex_only
@pytest.mark.parametrize("k", [0, 1, 2, 3])
def test_control_w_the_wall_is_felt(walled, k):
    assert walled["W"][k] > W_CONTROL_FLOOR, (
        f"(W-ctl) drive {k}: wall-on <theta>(t*) misses P1 by {walled['W'][k]:.3e} — h_eff not felt"
    )
