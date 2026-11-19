"""`MAT-17` step 3 — duty cycle on the coil-driven phantom (§9 item 271, chain P12).

Fixture: `MAT-15` step 1's (`_run_heated(return_fields=True)`, imported, not
edited): `ANS-2`'s F-small at 10 MHz, drive k = 0 (P1), G = 1 W / P₁, P1
thermal on the tag-3 submesh, step 1's CN grid (t* = 0.2 t_R, N = 40, two BE
start-up half-steps).

Schedule (pinned): period P = t*/5 = 8 steps, D = 0.5, χ = 1 on
[jP, jP + P/2), 0 on [jP + P/2, (j + 1)P).  Every switch t_m = m·P/2
(m = 0 … 9) is a grid point, so every (half-)step midpoint χ is exactly 0 / 1.

Asserted: (Δ) Duhamel superposition θ[χ] − T∞ = Σ (−1)^m (θ[H_m] − T∞), each
H_m a full run from t = 0, whole-submesh relative L2 ≤ ``L_BAND`` (`MAT-17`'s
(L) bar); (D′) χ ≡ 1 reproduces step 1's θ(t*) for k = 0 ≤ ``D_BAND`` (`MAT-17`'s
(D) bar); (E) every run's balance ≤ ``E_BAND``; (E′) from the EM side on the
duty run, |ρc∫Δθ/Δt + ∮hθ^(n+½) − χ_n·G·P₀| / (G·P₀) ≤ ``E_BAND`` at the last
step N (off, χ = 0) and at the last on step n = 36 (χ = 1).
Negative controls: (d) doubled source, asserted > ``D_CONTROL_FLOOR`` — on the
last ON step (n = 36): at N the source is off, a doubled source is invisible
to the balance there (θ is linear in q, stor + conv = 0 for both), so the
step-N reading is printed only (predicted ≈ 0); (Δ-ctl) m = 9 dropped,
predicted 0.15–0.3, asserted only as a miss > ``L_BAND``.

Scope: a 50 % duty, P = t*/5 schedule on P1's untuned 10 MHz drive k = 0, P1
thermal, gel constants, Robin h 10, no perfusion.  Not: any pulse sequence,
SAR-averaging or regulatory statement, a closed form for the coil-driven θ,
the tuned copper fixture, 64 MHz, `ANS-9`.
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
from tests.validation.test_birdcage_b1_plus_map import PHANTOM_RHO_KG_PER_M3
from tests.validation.test_mass_averaged_sar_standard_masses import (
    QUADRATURE_DEGREE,
    TEN_GRAM_KG,
)
from tests.validation.test_mat10_sar_to_heat_sphere import ROBIN_TAG, _allsum, _say
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

# `MAT-17`'s bars (`test_mat17_duty_cycle_fem.py:183` (D), `:201` (L)); that
# module asserts them as literals, so they are restated here, not moved.
L_BAND = 1.0e-12
D_BAND = 1.0e-13

DT = T_STAR / N_MAIN
STEPS_PER_HALF = 4  # P/2 = t*/10 = 4 steps
HALF = STEPS_PER_HALF * DT
DUTY = 0.5
N_SWITCH = 10  # t_m = m·P/2, m = 0 … 9
N_ON_LAST = N_MAIN - STEPS_PER_HALF  # 36: end of the last on phase
SNAPS = (N_MAIN - 2 * STEPS_PER_HALF, N_ON_LAST - 1, N_ON_LAST, N_MAIN - 1, N_MAIN)


def _chi_duty(t):
    """1 on [jP, jP + P/2), 0 on [jP + P/2, (j + 1)P); t is always a midpoint."""
    return 1.0 if int(np.floor(t / HALF)) % 2 == 0 else 0.0


def _heaviside(m):
    t_m = m * HALF
    return lambda t: 1.0 if t > t_m else 0.0


def _transient_s(sub, ft, q, schedule, *, snapshot_steps=SNAPS):
    """`MAT-15` step 1's ``_cn`` (``:145–157``) plus ``source_schedule=``."""
    return solve_transient_conduction(
        sub,
        ft,
        k=K_GEL,
        rho_c=RHO_C,
        q=q,
        h=H_GEL,
        t_inf=T_INF,
        t_init=T_INF,
        robin_tag=ROBIN_TAG,
        dt=DT,
        n_steps=N_MAIN,
        theta=0.5,
        startup_half_steps=2,
        degree=1,
        snapshot_steps=snapshot_steps,
        source_schedule=schedule,
    )


def _e_prime_at(sub, ds, res, n, deposited):
    """|ρc∫(θⁿ−θⁿ⁻¹)/Δt + ∮hθ^(n−½) − χ_n·G·P₀| / (G·P₀), χ_n from the run."""
    comm = sub.comm
    tn, tm = res.snapshots[n], res.snapshots[n - 1]
    # multipliers: two start-up half-steps, then full steps 2 … N
    chi_n = res.source_multipliers[n] if n >= 2 else None
    stor = float(np.real(_allsum(comm, RHO_C * (tn - tm) / DT * ufl.dx(domain=sub))))
    conv = float(np.real(_allsum(comm, H_GEL * (0.5 * (tn + tm) - T_INF) * ds)))
    return abs(stor + conv - chi_n * deposited) / abs(deposited), stor, conv, chi_n


def _run_duty(*, return_fields=False):
    """`duty()`'s body, lifted to module level (`EX-99`, rule (a)); ``return_fields=True``
    adds the fields and runs an example needs to the returned dict."""
    comm = MPI.COMM_WORLD
    base = _run_heated(return_fields=True)
    t0 = time.perf_counter()
    sub, ft, G, P = base["sub"], base["ft"], base["G"], base["P"]
    v_sub = base["v_sub"]
    ds = ufl.Measure("ds", domain=sub, subdomain_data=ft)(ROBIN_TAG)
    q0 = base["q_h"][0]
    dep = G * P[0]
    V1 = fem.functionspace(sub, ("Lagrange", 1))
    work = fem.Function(V1)

    def _mean_theta(T):
        return float(np.real(_allsum(comm, (T - T_INF) * ufl.dx(domain=sub)))) / v_sub

    def _l2sq(arr):
        work.x.array[:] = arr
        work.x.scatter_forward()
        return float(np.real(_allsum(comm, ufl.inner(work, work) * ufl.dx(domain=sub))))

    out = {"E": {}}
    tmax_imag = 0.0

    # Duty run.
    t1 = time.perf_counter()
    r_duty = _transient_s(sub, ft, q0, _chi_duty)
    t_duty = time.perf_counter() - t1
    chis = r_duty.source_multipliers
    out["chi_binary"] = all(c in (0.0, 1.0) for c in chis)
    out["E"]["duty"] = r_duty.max_balance_residual
    tmax_imag = max(tmax_imag, _max_imag(r_duty.temperature))
    th_duty = r_duty.temperature.x.array - T_INF
    _say(
        f"\nMAT17s3 schedule: dt = {DT:.4f} s, P/2 = {HALF:.4f} s ({STEPS_PER_HALF} steps), "
        f"D = {DUTY}; chi per (half-)step ({len(chis)}) = "
        + "".join(str(int(c)) if c in (0.0, 1.0) else "?" for c in chis)
        + f"; on-time {sum(c for c in chis[2:]) * DT + 0.5 * DT * sum(chis[:2]):.4f} s "
        f"of t* = {T_STAR:.4f} s; duty run {t_duty:.1f} s, (E) {r_duty.max_balance_residual:.3e}"
    )

    # (E′) on the duty run: last step (off) and last on step.
    out["Ep_N"], stor, conv, chiN = _e_prime_at(sub, ds, r_duty, N_MAIN, dep)
    _say(
        f"MAT17s3 (E') step N={N_MAIN} (chi_N = {chiN}): stored {stor:.12e} W + surface "
        f"{conv:.12e} W vs chi_N*G*P_1 = {chiN * dep:.12e} W, normalised by G*P_1 = {dep:.12e} W: "
        f"rel {out['Ep_N']:.3e} (bound {E_BAND:.0e})"
    )
    out["Ep_on"], stor, conv, chi36 = _e_prime_at(sub, ds, r_duty, N_ON_LAST, dep)
    _say(
        f"MAT17s3 (E') step n={N_ON_LAST} (chi = {chi36}): stored {stor:.12e} W + surface "
        f"{conv:.12e} W vs chi*G*P_1 = {chi36 * dep:.12e} W: rel {out['Ep_on']:.3e} "
        f"(bound {E_BAND:.0e})"
    )

    # (D′) χ ≡ 1 against step 1's T_star[0].
    r_one = _transient_s(sub, ft, q0, lambda t: 1.0)
    out["E"]["chi1"] = r_one.max_balance_residual
    tmax_imag = max(tmax_imag, _max_imag(r_one.temperature))
    ref = base["T_star"][0].x.array - T_INF
    out["D"] = np.sqrt(_l2sq(r_one.temperature.x.array - T_INF - ref) / _l2sq(ref))
    _say(
        f"MAT17s3 (D') chi=1 vs step 1 T_star[k=0]: whole-submesh rel L2 {out['D']:.3e} "
        f"(bound {D_BAND:.0e}); (E) {r_one.max_balance_residual:.3e}"
    )

    # (Δ) Duhamel superposition over ten full-from-zero Heaviside runs.
    terms = []
    for m in range(N_SWITCH):
        r_m = _transient_s(sub, ft, q0, _heaviside(m))
        out["E"][f"H{m}"] = r_m.max_balance_residual
        tmax_imag = max(tmax_imag, _max_imag(r_m.temperature))
        terms.append(r_m.temperature.x.array - T_INF)
        _say(
            f"MAT17s3   H_{m} (t_m = {m * HALF:.4f} s): <theta>(t*) = "
            f"{_mean_theta(r_m.temperature):.12e} K  ||theta||^2 = {_l2sq(terms[-1]):.6e}  "
            f"(E) {r_m.max_balance_residual:.3e}"
        )
    signs = [(-1.0) ** m for m in range(N_SWITCH)]
    sup = sum(s * t for s, t in zip(signs, terms))
    nrm = _l2sq(th_duty)
    out["Delta"] = np.sqrt(_l2sq(th_duty - sup) / nrm)
    sup9 = sup - signs[9] * terms[9]
    out["Delta_ctl"] = np.sqrt(_l2sq(th_duty - sup9) / nrm)
    work.x.array[:] = sup9
    work.x.scatter_forward()
    mean_sup9 = float(np.real(_allsum(comm, work * ufl.dx(domain=sub)))) / v_sub
    _say(
        f"MAT17s3 (Delta) theta[chi] - sum_m (-1)^m theta[H_m], m = 0..9: whole-submesh rel L2 "
        f"{out['Delta']:.3e} (bound {L_BAND:.0e}, predicted <= 1e-13)"
    )
    _say(
        f"MAT17s3 control (Delta-ctl) m = 9 dropped: rel L2 {out['Delta_ctl']:.6f} "
        f"(predicted 0.15-0.3, asserted > {L_BAND:.0e}); mean of the 9-term sum "
        f"{mean_sup9:.9f} K"
    )

    # (d) doubled source on the duty schedule.
    q_two = fem.Function(q0.function_space)
    q_two.x.array[:] = 2.0 * q0.x.array
    q_two.x.scatter_forward()
    r_d = _transient_s(sub, ft, q_two, _chi_duty)
    out["ctl_d"], _, _, _ = _e_prime_at(sub, ds, r_d, N_ON_LAST, dep)
    out["ctl_d_N"], _, _, _ = _e_prime_at(sub, ds, r_d, N_MAIN, dep)
    _say(
        f"MAT17s3 control (d) doubled source, duty: (E') vs chi*G*P_1 at n={N_ON_LAST} (on) = "
        f"{out['ctl_d']:.6f} (prior 1.000000, asserted > {D_CONTROL_FLOOR:.0e}); at N={N_MAIN} "
        f"(off, blind by linearity, printed) = {out['ctl_d_N']:.3e}; internal (E) "
        f"{r_d.max_balance_residual:.3e}"
    )

    # Records.
    mean_duty = _mean_theta(r_duty.temperature)
    mean_p1 = _mean_theta(base["T_star"][0])
    out["mean_ratio"] = mean_duty / (DUTY * mean_p1)
    _say(
        f"MAT17s3 record <theta>(t*): duty {mean_duty:.9f} K  vs D*P1 = {DUTY} * {mean_p1:.9f} K "
        f"= {DUTY * mean_p1:.9f} K  ratio {out['mean_ratio']:.6f} (predicted 0.97-1.00)"
    )
    s_lo, s_hi = SNAPS[0], N_ON_LAST
    m_lo = _mean_theta(r_duty.snapshots[s_lo])
    m_hi = _mean_theta(r_duty.snapshots[s_hi])
    _say(
        f"MAT17s3 record ripple, last on/off pair: <theta> n={s_lo} {m_lo:.9f} K -> n={s_hi} (on end) "
        f"{m_hi:.9f} K -> n={N_MAIN} (off end) {mean_duty:.9f} K;  (peak - end)/end = "
        f"{(m_hi - mean_duty) / mean_duty * 100:.4f} %, (peak - start)/start = "
        f"{(m_hi - m_lo) / m_lo * 100:.4f} %"
    )

    # 10 g ball mean at c_0 (step 1's forms, centre from the fixture).
    a10 = averaging_ball_radius(mass_kg=TEN_GRAM_KG, rho=PHANTOM_RHO_KG_PER_M3)
    centre_c = fem.Constant(sub, np.asarray(base["centres"][0], dtype=default_scalar_type))
    xs = ufl.SpatialCoordinate(sub)
    off = xs - centre_c
    inside = ufl.conditional(ufl.lt(ufl.real(ufl.dot(off, off)), a10**2), 1.0, 0.0)
    dxq = ufl.dx(domain=sub, metadata={"quadrature_degree": int(QUADRATURE_DEGREE)})
    th_ball = fem.Function(V1)
    f_vol = fem.form(inside * dxq)
    f_th = fem.form(th_ball * inside * dxq)

    def _red(form):
        return float(np.real(comm.allreduce(fem.assemble_scalar(form), op=MPI.SUM)))

    bv = _red(f_vol)
    th_ball.x.array[:] = th_duty
    th_ball.x.scatter_forward()
    b_duty = _red(f_th) / bv
    th_ball.x.array[:] = ref
    th_ball.x.scatter_forward()
    b_p1 = _red(f_th) / bv
    _say(
        f"MAT17s3 record 10 g ball theta(c_0) at t*: duty {b_duty:.9e} K  P1 {b_p1:.9e} K  "
        f"duty/(D*P1) = {b_duty / (DUTY * b_p1):.6f}"
    )

    worst = max(out["E"], key=out["E"].get)
    out["max_imag"] = tmax_imag
    _say(
        f"MAT17s3 (E) worst balance over {len(out['E'])} runs = {out['E'][worst]:.3e} on {worst} "
        f"(bound {E_BAND:.0e}); max|Im theta| = {tmax_imag:.3e} (predicted 0); "
        f"elapsed after P1 fixture {time.perf_counter() - t0:.1f} s"
    )
    if return_fields:
        out.update(
            base=base,
            r_duty=r_duty,
            r_one=r_one,
            r_d=r_d,
            q0=q0,
            dep=dep,
            sub=sub,
            ft=ft,
            ds=ds,
            ball10_duty=b_duty,
            ball10_p1=b_p1,
        )
    return out


@pytest.fixture(scope="module")
def duty():
    return _run_duty()


@complex_only
def test_schedule_midpoints_are_binary(duty):
    assert duty["chi_binary"], "a midpoint chi is not 0/1 — schedule-grid bug"


@complex_only
def test_d_prime_chi_one_reproduces_p1(duty):
    assert duty["D"] <= D_BAND, f"(D') {duty['D']:.3e} > {D_BAND:.0e}"


@complex_only
def test_delta_duhamel_superposition(duty):
    assert duty["Delta"] <= L_BAND, f"(Delta) {duty['Delta']:.3e} > {L_BAND:.0e}"


@complex_only
def test_e_every_run_balances(duty):
    for name, r in duty["E"].items():
        assert r <= E_BAND, f"(E) run {name}: {r:.3e} > {E_BAND:.0e}"


@complex_only
def test_e_prime_duty_balance_from_em_side(duty):
    assert duty["Ep_N"] <= E_BAND, f"(E') step N: {duty['Ep_N']:.3e} > {E_BAND:.0e}"
    assert duty["Ep_on"] <= E_BAND, f"(E') last on step: {duty['Ep_on']:.3e} > {E_BAND:.0e}"


@complex_only
def test_control_d_doubled_source_misses_the_balance(duty):
    assert (
        duty["ctl_d"] > D_CONTROL_FLOOR
    ), f"control (d) doubled source reads (E') {duty['ctl_d']:.3e} — did not miss"


@complex_only
def test_control_delta_dropped_term_misses(duty):
    assert (
        duty["Delta_ctl"] > L_BAND
    ), f"(Delta-ctl) m = 9 dropped reads {duty['Delta_ctl']:.3e} — did not miss"
