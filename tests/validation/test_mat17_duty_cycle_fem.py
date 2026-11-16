"""`MAT-17` step 2 — duty-cycled source in `solve_transient_conduction` (§9 item 200).

Fixture: `MAT-9`'s ball (R 0.1, k 0.6, ρc 4.18e6, q 2000, h 10, Bi 5/3,
T∞ = T(0) = 20), P1, source on while 0 ≤ (t mod P) < DP (step 1's `CASES`):
case A (P 1200 s, D 0.5) to t_end = 11P + 900 s = 14 100 s (mid off phase of
cycle 11); case B (P 0.2 t_R, D 0.25) to t_end = 2.5 DP (mid off phase of
cycle 0), Δt = DP/8.  Every Δt divides the phases, so every switch is a step
boundary and χ(t_mid) is exactly 0 or 1.  Comparand:
`robin_sphere_duty_cycle_temperature(…, periodic=False)` evaluated in numpy
at P2 dof coordinates (chunks of 20 000 points).  Errors on θ = T − T∞.

Asserted: (S0) Σ χₙΔtₙ = on-time to ≤ 1e-12 relative, every scheduled run;
(D) χ ≡ 1 equals the default call to ≤ 1e-13 (BE, Δt 300, 0.01 m);
(L) θ[χ] + θ[1 − χ] = θ[1] to ≤ 1e-12 (BE, Δt 150, 0.01 m);
(i) CN + Rannacher on 0.005 m, interior rel L2 ≤ 1 % in each case;
(ii) BE ladder Δt 300 … 18.75 on 0.01 m, case A, whole-ball self-difference
LSQ rate ∈ [0.9, 1.2]; (iii) every run's max balance ≤ 1e-9.
Bands imported by value from `MAT-9` step 2
(`test_mat9_transient_conduction_sphere.py:274–275, 301–302`).
Negative controls (rule (e), predicted, asserted only as missing the band):
(c1) schedule ignored vs case A's comparand, > 1e-2 (pred. ≈ 1.0);
(c2) schedule dropped from the balance on the off step 93→94, > 1e-3
(pred. 1.000); (c3) half-phase shift χ(t − 300 s), printed only
(pred. 2–4 %).
"""

from __future__ import annotations

import time

import numpy as np
import ufl
from dolfinx import fem
from mpi4py import MPI

from fem_em_solver.post.evaluation import evaluate_vector_field_parallel
from fem_em_solver.thermal import solve_transient_conduction
from fem_em_solver.utils.analytical import robin_sphere_duty_cycle_temperature
from tests.validation.test_mat9_transient_conduction_sphere import (
    H,
    K,
    Q,
    QDEG,
    R,
    RHO_C,
    ROBIN_TAG,
    T_INF,
    T_R,
    _allsum,
    _ball,
    _interior_error,
    _rel_diff,
    _say,
)
from tests.validation.test_mat17_duty_cycle_closed_form import CASES

P_A, D_A = CASES["A"]
P_B, D_B = CASES["B"]
T_END_A = 11 * P_A + 900.0  # 14 100 s, middle of cycle 11's off phase
T_END_B = 2.5 * D_B * P_B  # middle of cycle 0's off phase
DT_B = D_B * P_B / 8.0
ON_TIME_A = 12 * D_A * P_A  # 7 200 s
ON_TIME_B = D_B * P_B
BE_DTS = (300.0, 150.0, 75.0, 37.5, 18.75)

# MAT-9 step 2's bands, by value (literals at
# test_mat9_transient_conduction_sphere.py:274–275, 301–302).
INTERIOR_BAND = 1e-2
BALANCE_BAND = 1e-9
BE_RATE_BAND = (0.9, 1.2)

_CACHE: dict = {}
_T0 = time.perf_counter()


def _stamp(label):
    _say(MPI.COMM_WORLD, f"MAT17F t+{time.perf_counter() - _T0:.1f}s {label}")


def _chi(period, duty, shift=0.0):
    def chi(t):
        return 1.0 if ((t - shift) % period) < duty * period else 0.0

    return chi


def _n(t_end, dt):
    n = int(round(t_end / dt))
    assert abs(n * dt - t_end) <= 1e-9 * t_end
    return n


def _run(msh, ft, dt, n, theta, startup, schedule=None, snapshot_steps=()):
    t0 = time.perf_counter()
    s = solve_transient_conduction(
        msh, ft, k=K, rho_c=RHO_C, q=Q, h=H, t_inf=T_INF, t_init=T_INF,
        robin_tag=ROBIN_TAG, dt=dt, n_steps=n, theta=theta,
        startup_half_steps=startup, degree=1, snapshot_steps=snapshot_steps,
        source_schedule=schedule,
    )
    s.wall = time.perf_counter() - t0
    return s


def _census(s, dt, startup, on_time, label):
    """(S0): Σ χₙ Δtₙ against the on-time in [0, t_end]."""
    m = startup
    dts = [0.5 * dt] * m + [dt] * (len(s.source_multipliers) - m)
    tot = float(np.dot(np.asarray(s.source_multipliers), np.asarray(dts)))
    rel = abs(tot - on_time) / on_time
    n_on = int(sum(1 for x in s.source_multipliers if x == 1.0))
    n_off = int(sum(1 for x in s.source_multipliers if x == 0.0))
    _say(MPI.COMM_WORLD, f"MAT17F (S0) {label}: sum chi*dt={tot:.6f}s vs on-time {on_time:.6f}s "
         f"rel={rel:.3e} (bar 1e-12); (half-)steps on={n_on} off={n_off} "
         f"(of {len(s.source_multipliers)}, {m} half-steps)")
    assert rel <= 1e-12, f"(S0) {label}"
    return n_on, n_off


def _exact_duty(msh, t, period, duty):
    """θ(t) from step 1's closed form, interpolated in numpy into P2 (chunked)."""
    V2 = fem.functionspace(msh, ("Lagrange", 2))
    f = fem.Function(V2)

    def ev(x):
        # The comparand depends on r only, and interpolation points repeat
        # across cells: evaluate once per distinct r (identical values; the
        # per-point form cost 37 s per call on 0.005 m, probe
        # 20261003T170842Z_MAT-17-step2-probe.log:63,66).
        r_all = np.sqrt(np.sum(x**2, axis=0))
        r, inv = np.unique(r_all, return_inverse=True)
        out = np.empty(r.shape)
        for a in range(0, r.size, 20000):
            rr = r[a:a + 20000]
            out[a:a + 20000] = robin_sphere_duty_cycle_temperature(
                rr, np.full(rr.shape, t), radius=R, k=K, rho_c=RHO_C, q=Q, h=H,
                t_inf=T_INF, period=period, duty=duty, periodic=False,
            ) - T_INF
        return out[inv.reshape(-1)]

    t0 = time.perf_counter()
    f.interpolate(ev)
    _stamp(f"closed-form P2 interpolation t={t:.1f} ({time.perf_counter() - t0:.1f}s)")
    return f


def _theta_centre(T):
    try:
        vals, ok = evaluate_vector_field_parallel(T, np.zeros((1, 3)))
        return float(np.real(np.ravel(vals[0])[0])) - T_INF if bool(ok[0]) else float("nan")
    except Exception as exc:  # printed quantity only
        _say(MPI.COMM_WORLD, f"MAT17F centre evaluation failed: {exc!r}")
        return float("nan")


def _closed(r, t, period, duty):
    r = np.atleast_1d(np.asarray(r, dtype=float))
    return robin_sphere_duty_cycle_temperature(
        r, np.full(r.shape, t), radius=R, k=K, rho_c=RHO_C, q=Q, h=H, t_inf=T_INF,
        period=period, duty=duty, periodic=False,
    ) - T_INF


# ------------------------------------------------------------- (D), (L), (S0)


def test_d_default_path_and_l_linearity():
    comm = MPI.COMM_WORLD
    _stamp("start (D)/(L)")
    msh, _, ft = _ball(0.01)
    _stamp("0.01 m mesh ready")

    # (D) χ ≡ 1 against the default call: BE, Δt 300, 47 steps.
    n = _n(T_END_A, 300.0)
    s_def = _run(msh, ft, 300.0, n, 1.0, 0)
    s_one = _run(msh, ft, 300.0, n, 1.0, 0, schedule=lambda t: 1.0)
    _census(s_one, 300.0, 0, T_END_A, "(D) chi=1 BE dt=300")
    d_D = _rel_diff(msh, s_one.temperature, s_def.temperature)
    _say(comm, f"MAT17F (D) chi=1 vs default, BE dt=300 N={n}: whole-ball rel={d_D:.3e} "
         f"(bar 1e-13, predicted <= 1e-15); balance default={s_def.max_balance_residual:.3e} "
         f"chi=1={s_one.max_balance_residual:.3e}; wall {s_def.wall:.1f}s / {s_one.wall:.1f}s")
    assert s_def.source_multipliers == []
    assert d_D <= 1e-13, "(D)"

    # (L) θ[χ] + θ[1 − χ] = θ[1]: BE, Δt 150, 94 steps.
    n = _n(T_END_A, 150.0)
    chi = _chi(P_A, D_A)
    s_a = _run(msh, ft, 150.0, n, 1.0, 0, schedule=chi)
    s_b = _run(msh, ft, 150.0, n, 1.0, 0, schedule=lambda t: 1.0 - chi(t))
    s_1 = _run(msh, ft, 150.0, n, 1.0, 0)
    _census(s_a, 150.0, 0, ON_TIME_A, "(L) chi BE dt=150")
    _census(s_b, 150.0, 0, T_END_A - ON_TIME_A, "(L) 1-chi BE dt=150")
    dx = ufl.dx(domain=msh, metadata=QDEG)
    e = _allsum(comm, (s_a.temperature + s_b.temperature - T_INF - s_1.temperature) ** 2 * dx)
    nrm = _allsum(comm, (s_1.temperature - T_INF) ** 2 * dx)
    d_L = float(np.sqrt(np.real(e) / np.real(nrm)))
    worst = max(s.max_balance_residual for s in (s_def, s_one, s_a, s_b, s_1))
    _say(comm, f"MAT17F (L) theta[chi]+theta[1-chi] vs theta[1], BE dt=150 N={n}: whole-ball "
         f"rel={d_L:.3e} (bar 1e-12, predicted ≈ 1e-15); (iii) worst balance={worst:.3e}; "
         f"wall {s_a.wall:.1f}s / {s_b.wall:.1f}s / {s_1.wall:.1f}s")
    assert d_L <= 1e-12, "(L)"
    assert worst <= BALANCE_BAND, "(iii) on (D)/(L) runs"


# ------------------------------------------------------- (i) and controls


def test_i_interior_and_controls_c1_c2_c3():
    comm = MPI.COMM_WORLD
    _stamp("start (i)")
    msh, _, ft = _ball(0.005)
    _stamp("0.005 m mesh ready")
    n_glob = msh.topology.index_map(3).size_global
    chi_a = _chi(P_A, D_A)
    chi_b = _chi(P_B, D_B)

    # (i) case A: CN + Rannacher, Δt 150, N 94.
    n_a = _n(T_END_A, 150.0)
    s_a = _run(msh, ft, 150.0, n_a, 0.5, 2, schedule=chi_a,
               snapshot_steps=(n_a - 1, n_a))
    n_on, n_off = _census(s_a, 150.0, 2, ON_TIME_A, f"(i) case A CN N={n_a}")
    ex_a = _exact_duty(msh, T_END_A, P_A, D_A)
    e_a = _interior_error(msh, s_a.temperature, ex_a)
    th0_a = _theta_centre(s_a.temperature)
    th0_ex_a = float(_closed(0.0, T_END_A, P_A, D_A)[0])
    _say(comm, f"MAT17F (i) case A CN+Rannacher dt=150 N={n_a} res=0.005 cells={n_glob}: "
         f"interior_rel_L2={e_a:.4e} (band 1e-2, predicted 5e-4 – 2e-3) "
         f"max_balance={s_a.max_balance_residual:.3e} wall {s_a.wall:.1f}s; "
         f"theta_h(0)={th0_a:.6f} K vs closed form {th0_ex_a:.6f} K; "
         f"(half-)steps on={n_on} off={n_off} (predicted 48 on / 46 off full-step equivalents)")

    # (i) case B: CN + Rannacher, Δt DP/8, N 20.
    n_b = _n(T_END_B, DT_B)
    s_b = _run(msh, ft, DT_B, n_b, 0.5, 2, schedule=chi_b)
    _census(s_b, DT_B, 2, ON_TIME_B, f"(i) case B CN N={n_b}")
    ex_b = _exact_duty(msh, T_END_B, P_B, D_B)
    e_b = _interior_error(msh, s_b.temperature, ex_b)
    th0_b = _theta_centre(s_b.temperature)
    th0_ex_b = float(_closed(0.0, T_END_B, P_B, D_B)[0])
    _say(comm, f"MAT17F (i) case B CN+Rannacher dt={DT_B:.4f} N={n_b}: interior_rel_L2="
         f"{e_b:.4e} (band 1e-2, predicted 5e-4 – 2e-3) max_balance="
         f"{s_b.max_balance_residual:.3e} wall {s_b.wall:.1f}s; theta_h(0)={th0_b:.6f} K "
         f"vs closed form {th0_ex_b:.6f} K")

    for lab, e, s, t_end, P, D, dt in (("A", e_a, s_a, T_END_A, P_A, D_A, 150.0),
                                       ("B", e_b, s_b, T_END_B, P_B, D_B, DT_B)):
        if e > INTERIOR_BAND:  # negative-result diagnostics (brief)
            _say(comm, f"MAT17F (i) RED case {lab}: closed form at r=0,R/2,R = "
                 f"{_closed([0.0, R / 2, R], t_end, P, D)}; last 12 chi = "
                 f"{s.source_multipliers[-12:]}")

    # Printed: CN without the Rannacher start, case A.
    s_ns = _run(msh, ft, 150.0, n_a, 0.5, 0, schedule=chi_a)
    e_ns = _interior_error(msh, s_ns.temperature, ex_a)
    _say(comm, f"MAT17F printed: case A CN no-startup N={n_a}: interior_rel_L2={e_ns:.4e} "
         f"max_balance={s_ns.max_balance_residual:.3e} wall {s_ns.wall:.1f}s")
    _census(s_ns, 150.0, 0, ON_TIME_A, "case A CN no-startup")

    # (c1) schedule ignored: default call against case A's comparand.
    s_c1 = _run(msh, ft, 150.0, n_a, 0.5, 2)
    e_c1 = _interior_error(msh, s_c1.temperature, ex_a)
    _say(comm, f"MAT17F control (c1) schedule ignored: interior_rel_L2={e_c1:.4e} "
         f"(asserted > 1e-2; predicted ≈ 1.0, range 0.8–1.2) wall {s_c1.wall:.1f}s")

    # (c2) schedule dropped from the balance on the off step 93→94.
    T93, T94 = s_a.snapshots[n_a - 1], s_a.snapshots[n_a]
    ds = ufl.Measure("ds", domain=msh, subdomain_data=ft)
    dx = ufl.dx(domain=msh)
    stor = _allsum(comm, RHO_C * (T94 - T93) / 150.0 * dx)
    conv = _allsum(comm, H * (0.5 * T94 + 0.5 * T93 - T_INF) * ds(ROBIN_TAG))
    bal_c2 = abs(float(np.real(stor + conv)) - s_a.source_power) / s_a.source_power
    _say(comm, f"MAT17F control (c2) schedule dropped from balance, step {n_a - 1}->{n_a} "
         f"(chi={s_a.source_multipliers[-1]}): balance={bal_c2:.4e} "
         f"(asserted > 1e-3; predicted 1.000)")

    # (c3) half-phase shift, printed only.
    s_c3 = _run(msh, ft, 150.0, n_a, 0.5, 2, schedule=_chi(P_A, D_A, shift=300.0))
    e_c3 = _interior_error(msh, s_c3.temperature, ex_a)
    _say(comm, f"MAT17F control (c3) chi(t-300s), printed only: interior_rel_L2={e_c3:.4e} "
         f"(predicted 2–4 %) wall {s_c3.wall:.1f}s")

    worst = max(s.max_balance_residual for s in (s_a, s_b, s_ns, s_c1, s_c3))
    _say(comm, f"MAT17F (iii) worst balance on the 0.005 m runs={worst:.3e} "
         f"(bar 1e-9, predicted <= 5e-14)")
    assert worst <= BALANCE_BAND, "(iii) on the (i) runs"
    assert e_a <= INTERIOR_BAND, f"(i) case A {e_a:.4e} > 1 %"
    assert e_b <= INTERIOR_BAND, f"(i) case B {e_b:.4e} > 1 %"
    assert e_c1 > 1e-2, "(c1)"
    assert bal_c2 > 1e-3, "(c2)"


# --------------------------------------------------------------- (ii) rate


def _ladder(theta, startup, label):
    comm = MPI.COMM_WORLD
    msh, _, ft = _ball(0.01)
    chi = _chi(P_A, D_A)
    sols, worst, last = [], 0.0, None
    t0 = time.perf_counter()
    for dt in BE_DTS:
        n = _n(T_END_A, dt)
        s = _run(msh, ft, dt, n, theta, startup, schedule=chi)
        _census(s, dt, startup, ON_TIME_A, f"(ii) {label} dt={dt}")
        worst = max(worst, s.max_balance_residual)
        _say(comm, f"MAT17F (ii) {label} dt={dt} N={n} max_balance={s.max_balance_residual:.3e} "
             f"wall {s.wall:.1f}s")
        sols.append(s.temperature)
        last = s
        _stamp(f"{label} dt={dt} done")
    d = [_rel_diff(msh, sols[i], sols[i + 1]) for i in range(len(sols) - 1)]
    lsq = float(np.polyfit(np.log(np.asarray(BE_DTS[: len(d)])), np.log(d), 1)[0])
    pair = [float(np.log(d[i] / d[i + 1]) / np.log(2.0)) for i in range(len(d) - 1)]
    _say(comm, f"MAT17F (ii) {label} d (dt=300,150,75,37.5; whole ball, ratio 2) = "
         f"{[f'{x:.4e}' for x in d]} LSQ rate={lsq:.3f} pairwise={[f'{p:.3f}' for p in pair]} "
         f"({time.perf_counter() - t0:.1f}s)")
    return lsq, pair, worst, last, msh


def test_ii_be_rate_and_iii_balance():
    comm = MPI.COMM_WORLD
    lsq_be, pair_be, w_be, last_be, msh = _ladder(1.0, 0, "BE")
    ex = _exact_duty(msh, T_END_A, P_A, D_A)
    e752 = _interior_error(msh, last_be.temperature, ex)
    _say(comm, f"MAT17F printed: BE N=752 interior_rel_L2 vs closed form={e752:.4e} "
         f"(predicted ≈ 2.8e-3, the 0.01 m floor)")
    lsq_cn, _, w_cn, _, _ = _ladder(0.5, 2, "CN+Rannacher (printed only)")
    _say(comm, f"MAT17F printed: CN+Rannacher LSQ rate case A={lsq_cn:.3f} (predicted ≈ 2.0); "
         f"(iii) worst balance BE {w_be:.3e}, CN {w_cn:.3e}")
    assert w_be <= BALANCE_BAND and w_cn <= BALANCE_BAND, "(iii)"
    assert BE_RATE_BAND[0] <= lsq_be <= BE_RATE_BAND[1], (
        f"(ii) BE LSQ rate {lsq_be:.3f}, pairwise {pair_be}")
