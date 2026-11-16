"""`MAT-15` step 4 — the Pennes sink on the coil-driven phantom (§9 item 315, chain P15).

Fixture: `MAT-15` step 1's (`_run_heated(return_fields=True)`, imported, not
edited): `ANS-2`'s F-small at 10 MHz, σ = 800 S/m, the 0.0025 m phantom rung,
the four single-leg drives, G = 1 W / P₁, P1 thermal on the tag-3 submesh with
the gel constants (k 0.6, ρc 4.18e6, h 10, T∞ = T(0) = 20).

The sink is ``solve_transient_conduction(..., perfusion_w_per_m3_k=w)`` (the
`MAT-11` / `MAT-13` term w(T − T∞), arterial = ambient), at w ∈ {24, 240, 2400}
W/m³/K (``W``/10, ``W``, ``W_STRONG`` from `MAT-11`) and the limit value
w_ε = 1e-8, which is non-zero on purpose so it takes the solver's ``perfused``
branch.  **Gel has no perfusion**: these w are a numerical exercise of the sink
on a coil-driven source, not a gel or tissue property.

Asserted: (E′) ρc∫(θᴺ−θᴺ⁻¹)/Δt + ∮hθ^(n+½) + ∫wθ^(n+½) = G·P_k from the EM side
≤ ``E_BAND`` (four drives × three w) and (E) every sink run's balance ≤
``E_BAND``; (L0) ⟨θ⟩(t*) at w_ε against P1's ≤ ``B_BAND`` per drive; (Mo)
⟨θ⟩(t*) and the own-centre 10 g θ̄ strictly decreasing over w ∈ {0 (P1), 24,
240, 2400}, every relative gap > ``MO_GAP_FLOOR``; (C4) the own-centre 10 g θ̄
in cyclic pairs ≤ ``C4_COVARIANCE_BAND`` at each w.  Negative controls: (d)
doubled source at w = ``W`` misses (E′) by > ``D_CONTROL_FLOOR`` (asserted,
prior same-fixture measurement 1.000000); (P-ctl) the sink dropped from the
identity at w = ``W`` misses by the perfused share, predicted [0.04, 0.08],
asserted only > ``D_CONTROL_FLOOR``.

Scope: a uniform numerical sink on P1's untuned 10 MHz σ = 800 S/m fixture, P1
thermal, gel constants.  Not: any gel or tissue perfusion value, a Pennes
closed form on the cylinder, σ(T) or w(T), tuned / copper / 64 MHz, absolute
temperature or compliance, `ANS-9`.
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
from tests.validation.test_mat11_perfusion_sphere import W, W_STRONG
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

W_VALUES = (W / 10.0, W, W_STRONG)  # 24, 240, 2400 W/m^3/K — fixture parameters, not bands
W_EPS = 1.0e-8  # the (L0) limit value; non-zero so the perfused branch is taken
MO_GAP_FLOOR = 1.0e-6  # (Mo) each consecutive relative gap; ≥ 1e3× under the smallest predicted


def _transient_w(sub, ft, q, w, *, dt, n_steps, theta, startup_half_steps=0, snapshot_steps=()):
    """`MAT-15` step 1's ``_transient`` with the Pennes sink ``w`` as an argument."""
    return solve_transient_conduction(
        sub, ft, k=K_GEL, rho_c=RHO_C, q=q, h=H_GEL, t_inf=T_INF, t_init=T_INF,
        robin_tag=ROBIN_TAG, dt=dt, n_steps=n_steps, theta=theta,
        startup_half_steps=startup_half_steps, degree=1,
        snapshot_steps=snapshot_steps, perfusion_w_per_m3_k=w,
    )


def _cn_w(sub, ft, q, w):
    """CN, N = 40 to t*, two BE half-steps first; snapshots N−1, N (step 1's ``_cn``)."""
    return _transient_w(sub, ft, q, w, dt=T_STAR / N_MAIN, n_steps=N_MAIN, theta=0.5,
                        startup_half_steps=2, snapshot_steps=(N_MAIN - 1, N_MAIN))


def _e_prime_w(sub, ds, res, deposited, w):
    """ρc∫(θᴺ−θᴺ⁻¹)/Δt + ∮hθ^(n+½) + ∫wθ^(n+½) against the EM-side deposited power.

    Returns (rel, stor, conv, perf); the (P-ctl) control reads stor + conv alone.
    """
    comm = sub.comm
    tn, tm = res.snapshots[N_MAIN], res.snapshots[N_MAIN - 1]
    dt = T_STAR / N_MAIN
    half = 0.5 * (tn + tm) - T_INF
    stor = float(np.real(_allsum(comm, RHO_C * (tn - tm) / dt * ufl.dx(domain=sub))))
    conv = float(np.real(_allsum(comm, H_GEL * half * ds)))
    perf = float(np.real(_allsum(comm, w * half * ufl.dx(domain=sub))))
    return abs(stor + conv + perf - deposited) / abs(deposited), stor, conv, perf


@pytest.fixture(scope="module")
def perfused():
    return _run_perfused()


def _run_perfused(*, return_fields=False):
    """The gate's whole body.

    ``return_fields=True`` additionally returns ``base`` (step 1's
    ``_run_heated(return_fields=True)`` dict) and ``results`` keyed ``(k, w)``
    for w in ``(W_EPS,) + W_VALUES`` (EX-125, rule (a), additive); prints and
    the other keys are unchanged.
    """
    comm = MPI.COMM_WORLD
    results = {}
    base = _run_heated(return_fields=True)
    t0 = time.perf_counter()
    sub, ft, G, P = base["sub"], base["ft"], base["G"], base["P"]
    v_sub = base["v_sub"]
    ds = ufl.Measure("ds", domain=sub, subdomain_data=ft)(ROBIN_TAG)

    # Ball-mean forms (item 270's block, centres from the fixture).
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

    def _th_own(T, k):
        th_ball.x.array[:] = T.x.array - T_INF
        th_ball.x.scatter_forward()
        return _ball_mean(f_th, k)

    out = {"Ep": {}, "E": {}, "L0": {}, "mean": {}, "th": {}, "q_own": {}, "pctl": {},
           "perf_power": {}}
    tmax_imag = 0.0
    _say(f"\nMAT15s4 w values {W_VALUES} W/m^3/K (numerical exercise; gel has no perfusion), "
         f"w_eps = {W_EPS:.0e}; t* = {T_STAR:.3f} s, N = {N_MAIN}")
    for k in range(4):
        q_h = base["q_h"][k]
        dep = G * P[k]
        mean_p1 = _mean_theta(base["T_star"][k])
        th_p1 = _th_own(base["T_star"][k], k)
        out["mean"][(k, 0.0)] = mean_p1
        out["th"][(k, 0.0)] = th_p1
        q_ball.x.array[:] = q_h.x.array
        q_ball.x.scatter_forward()
        out["q_own"][k] = _ball_mean(f_q, k)

        t1 = time.perf_counter()
        res_e = _cn_w(sub, ft, q_h, W_EPS)
        t_e = time.perf_counter() - t1
        if return_fields:
            results[(k, W_EPS)] = res_e
        out["E"][f"k{k}_weps"] = res_e.max_balance_residual
        mean_e = _mean_theta(res_e.temperature)
        out["L0"][k] = abs(mean_e / mean_p1 - 1.0)
        ep_e, _, _, perf_e = _e_prime_w(sub, ds, res_e, dep, W_EPS)
        tmax_imag = max(tmax_imag, _max_imag(res_e.temperature))
        _say(f"MAT15s4 drive k={k} ({base['by_k'][k]}): P1 <theta>(t*) {mean_p1:.12f} K, 10 g "
             f"theta(c_k) {th_p1:.9e} K;  q_h(c_k) {out['q_own'][k]:.9e} W/m^3")
        _say(f"MAT15s4   w_eps: CN {t_e:.1f} s; <theta>(t*) {mean_e:.12f} K; (L0) rel "
             f"{out['L0'][k]:.3e} (bound {B_BAND:.0e}, estimate w_eps t*/(2 rho c) = "
             f"{W_EPS * T_STAR / (2.0 * RHO_C):.2e}); (E') {ep_e:.3e} (printed); perf {perf_e:.3e} W; "
             f"(E) {res_e.max_balance_residual:.3e}")

        for w in W_VALUES:
            t1 = time.perf_counter()
            res = _cn_w(sub, ft, q_h, w)
            t_cn = time.perf_counter() - t1
            if return_fields:
                results[(k, w)] = res
            ep, stor, conv, perf = _e_prime_w(sub, ds, res, dep, w)
            out["Ep"][(k, w)] = ep
            out["E"][f"k{k}_w{w:g}"] = res.max_balance_residual
            mean_w = _mean_theta(res.temperature)
            th_w = _th_own(res.temperature, k)
            out["mean"][(k, w)] = mean_w
            out["th"][(k, w)] = th_w
            out["perf_power"][(k, w)] = res.perfused_power
            tmax_imag = max(tmax_imag, _max_imag(res.temperature))
            if w == W:
                out["pctl"][k] = abs(stor + conv - dep) / abs(dep)
            _say(f"MAT15s4   w = {w:g}: CN {t_cn:.1f} s; (E) {res.max_balance_residual:.3e}; (E') "
                 f"stor {stor:.12e} + conv {conv:.12e} + perf {perf:.12e} vs G*P_k {dep:.12e} W: "
                 f"rel {ep:.3e} (bound {E_BAND:.0e}); result.perfused_power {res.perfused_power:.12e} W")
            _say(f"MAT15s4     <theta>(t*) {mean_w:.9f} K, ratio to P1 {mean_w / mean_p1:.6f};  "
                 f"10 g theta(c_k) {th_w:.9e} K, ratio {th_w / th_p1:.6f}"
                 + (f";  (P-ctl) sink dropped: {out['pctl'][k]:.6f} (predicted [0.04, 0.08], "
                    f"asserted > {D_CONTROL_FLOOR:.0e})" if w == W else ""))

    # (Mo) monotone in w, per drive: P1, 24, 240, 2400.
    seq = (0.0,) + W_VALUES
    out["mo_gaps"] = {}
    for k in range(4):
        for key in ("mean", "th"):
            vals = [out[key][(k, w)] for w in seq]
            gaps = [(vals[i] - vals[i + 1]) / vals[i] for i in range(len(vals) - 1)]
            out["mo_gaps"][(k, key)] = gaps
            _say(f"MAT15s4 (Mo) k={k} {key}: relative gaps "
                 + "  ".join(f"{g:.6e}" for g in gaps)
                 + f" (each asserted > {MO_GAP_FLOOR:.0e})")

    # (C4) per w, own-centre 10 g theta in cyclic pairs.
    out["pairs"] = {}
    for w in seq:
        for k in range(4):
            a, b = out["th"][(k, w)], out["th"][((k + 1) % 4, w)]
            out["pairs"][(w, k)] = abs(b - a) / a
        _say(f"MAT15s4 (C4) w = {w:g}: theta 10 g pairs "
             + "  ".join(f"{out['pairs'][(w, k)] * 100:.4f}%" for k in range(4))
             + f" (bound {C4_COVARIANCE_BAND * 100:.1f}%)")
    q_pairs = [abs(out["q_own"][(k + 1) % 4] - out["q_own"][k]) / out["q_own"][k] for k in range(4)]
    _say("MAT15s4 q_h 10 g pairs (printed): " + "  ".join(f"{p * 100:.4f}%" for p in q_pairs))
    for w in W_VALUES:
        kpk = max(range(4), key=lambda j: out["th"][(j, w)])
        _say(f"MAT15s4 peak 10 g theta at w = {w:g}: {out['th'][(kpk, w)]:.9e} K on k={kpk}; "
             f"mean ratio to P1 over drives "
             f"{np.mean([out['mean'][(j, w)] / out['mean'][(j, 0.0)] for j in range(4)]):.6f}")

    # Negative control (d): doubled source on drive k = 0 at w = W.
    q0 = base["q_h"][0]
    q_two = fem.Function(q0.function_space)
    q_two.x.array[:] = 2.0 * q0.x.array
    q_two.x.scatter_forward()
    res_d = _cn_w(sub, ft, q_two, W)
    out["ctl_d"], _, _, _ = _e_prime_w(sub, ds, res_d, G * P[0], W)
    _say(f"MAT15s4 control (d) doubled source, k=0, w = {W:g}: (E') vs G*P_1 = {out['ctl_d']:.6f} "
         f"(prior 1.000000, asserted > {D_CONTROL_FLOOR:.0e}); internal (E) = "
         f"{res_d.max_balance_residual:.3e}")

    worst = max(out["E"], key=out["E"].get)
    out["max_imag"] = tmax_imag
    _say(f"MAT15s4 (E) worst balance over {len(out['E'])} runs = {out['E'][worst]:.3e} on {worst} "
         f"(bound {E_BAND:.0e}); max|Im theta| = {tmax_imag:.3e} (predicted 0)")
    _say(f"MAT15s4 elapsed after P1 fixture: {time.perf_counter() - t0:.1f} s")
    if return_fields:
        out.update(base=base, results=results)
    return out


@complex_only
@pytest.mark.parametrize("w", W_VALUES)
@pytest.mark.parametrize("k", [0, 1, 2, 3])
def test_e_prime_with_sink_equals_em_deposited(perfused, k, w):
    r = perfused["Ep"][(k, w)]
    assert r <= E_BAND, f"(E') drive {k}, w = {w:g}: {r:.3e} > {E_BAND:.0e}"


@complex_only
def test_e_every_sink_run_balances(perfused):
    for name, r in perfused["E"].items():
        assert r <= E_BAND, f"(E) run {name}: {r:.3e} > {E_BAND:.0e}"


@complex_only
@pytest.mark.parametrize("k", [0, 1, 2, 3])
def test_l0_vanishing_sink_reproduces_p1(perfused, k):
    assert perfused["L0"][k] <= B_BAND, f"(L0) drive {k}: {perfused['L0'][k]:.3e} > {B_BAND:.0e}"


@complex_only
@pytest.mark.parametrize("key", ["mean", "th"])
@pytest.mark.parametrize("k", [0, 1, 2, 3])
def test_mo_theta_strictly_decreasing_in_w(perfused, k, key):
    gaps = perfused["mo_gaps"][(k, key)]
    for i, g in enumerate(gaps):
        assert g > MO_GAP_FLOOR, f"(Mo) drive {k} {key} gap {i}: {g:.3e} not > {MO_GAP_FLOOR:.0e}"


@complex_only
@pytest.mark.parametrize("w", W_VALUES)
@pytest.mark.parametrize("k", [0, 1, 2, 3])
def test_c4_sink_ten_gram_theta_follows_the_drive(perfused, k, w):
    r = perfused["pairs"][(w, k)]
    assert r <= C4_COVARIANCE_BAND, f"(C4) w = {w:g} pair {k}: {r * 100:.4f}% > 5 %"


@complex_only
def test_control_d_doubled_source_misses_the_balance(perfused):
    assert perfused["ctl_d"] > D_CONTROL_FLOOR, (
        f"control (d) doubled source reads (E') {perfused['ctl_d']:.3e} — did not miss"
    )


@complex_only
@pytest.mark.parametrize("k", [0, 1, 2, 3])
def test_control_pctl_dropping_the_sink_misses(perfused, k):
    assert perfused["pctl"][k] > D_CONTROL_FLOOR, (
        f"(P-ctl) drive {k}: sink-free identity reads {perfused['pctl'][k]:.3e} — sink not felt"
    )
