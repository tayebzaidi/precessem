"""`MAT-18` step 2(b): piecewise k and ρc in the transient solver (§9 item 205).

2(a)'s conforming two-region ball (gel core r < a = 0.1, k₁ 0.6, ρc₁ 4.18e6,
q₁ 2000; container-wall shell a < r < b = 0.11, k₂ 0.19, ρc₂ 1.7e6, q₂ 0;
Robin h 10 / T∞ 20 on r = b), P1, T(0) = T∞, source on for t > 0, to
t* = 0.2·a²ρc₁/k₁ (Fo 0.2 on the core time; numerically `MAT-9`'s T_END).
Comparand: step 1's ``two_layer_sphere_transient_temperature`` (400 terms),
evaluated in numpy at P2 dof coordinates in 20 000-point chunks (`MAT-9`'s
``_exact_theta``).  Bands are `MAT-9` step 2's literals
(``test_mat9_transient_conduction_sphere.py:274–275, 301, 321–322``) and
2(a)'s (D) literal.

Anchors (asserted): (T1) CN + Rannacher N 40 at 0.005 with DG0 k and ρc,
interior rel L2 ≤ 1 %; (T2) every run's max per-step balance ≤ 1e-9
(normalised by ∫q = q₁V_core,h); (T3) BE 20 × 5 t_R equals the steady solve
≤ 1e-8 at 0.01; (D) DG0 k ≡ K1, ρc ≡ RC1 vs scalar, BE N 10 at 0.01, ≤ 1e-13;
(E) equal properties vs `robin_sphere_transient_temperature(radius=b)`,
interior ≤ 1 %.  Negative controls (rule (e), predicted, asserted only as
missing the anchor's band): (c1) uniform ρc = RC1, predicted interior miss
≈ 0.05 (0.02–0.10), asserted > 1e-2; (c2) ρc₁ storage weight in the balance,
predicted ≈ 0.07 (0.03–0.15), asserted > 1e-3.
"""

from __future__ import annotations

import time

import numpy as np
import ufl
from dolfinx import default_scalar_type, fem
from mpi4py import MPI

from fem_em_solver.post.evaluation import evaluate_vector_field_parallel
from fem_em_solver.thermal import solve_steady_conduction, solve_transient_conduction
from fem_em_solver.utils.analytical import (
    robin_sphere_transient_temperature,
    two_layer_sphere_transient_temperature,
)
from tests.validation.test_mat18_two_layer_sphere_closed_form import (
    A_,
    B_,
    K1,
    K2,
    Q1,
    Q2,
    RC1,
    RC2,
    T_INF,
    T_R,
    H,
)
from tests.validation.test_mat18_two_layer_sphere_fem import (
    _allsum,
    _ball,
    _dg0,
    _indicator,
    _vcore,
)

ROBIN_TAG = 1
T_STAR = 0.2 * T_R
N_CN = 40
CELL_GUARD = 400_000
QDEG = {"quadrature_degree": 4}
V_CORE = 4.0 / 3.0 * np.pi * A_**3
V_SHELL = 4.0 / 3.0 * np.pi * (B_**3 - A_**3)
CHUNK = 20_000

_CACHE: dict = {}


def _say(msg):
    if MPI.COMM_WORLD.rank == 0:
        print(msg, flush=True)


def _two_layer(r, t):
    return two_layer_sphere_transient_temperature(
        r,
        t,
        a=A_,
        b=B_,
        k1=K1,
        k2=K2,
        rho_c1=RC1,
        rho_c2=RC2,
        q1=Q1,
        q2=Q2,
        h=H,
        t_inf=T_INF,
        t_init=T_INF,
    )


def _homog(r, t):
    return robin_sphere_transient_temperature(
        r,
        t,
        radius=B_,
        k=K1,
        rho_c=RC1,
        q=Q1,
        h=H,
        t_inf=T_INF,
        t_init=T_INF,
    )


def _exact_theta(msh, fn, t):
    """θ(t) = T − T∞ from a series, in numpy at the P2 dof coordinates
    (owned + ghost), chunked.  Lagrange dofs are point values, so this is the
    P2 interpolant.  Measured 2026-10-03 (`…184829Z_MAT-18-step2b-probe-
    comparand.log`): one 20 000-point two-layer call costs 2.2 s, so
    ``f.interpolate``'s per-cell points (10 per tet, ≈ 1.1M per rank at 0.005)
    overran W1's 180 s (`…184348Z_MAT-18-step2b-W1.log`, status 124); the
    dof coordinates are ≈ 7× fewer."""
    V2 = fem.functionspace(msh, ("Lagrange", 2))
    f = fem.Function(V2)
    x = V2.tabulate_dof_coordinates()
    r = np.sqrt(np.sum(x**2, axis=1))
    out = np.empty(r.shape)
    times = []
    t0 = time.perf_counter()
    for s in range(0, r.size, CHUNK):
        tc = time.perf_counter()
        out[s : s + CHUNK] = fn(r[s : s + CHUNK], t) - T_INF
        times.append(time.perf_counter() - tc)
    f.x.array[:] = out
    f.x.scatter_forward()
    return f, time.perf_counter() - t0, (max(times) if times else 0.0)


def _errors(msh, ct, T_h, theta):
    """(whole, interior, shell-only) relative L2 errors on θ."""
    comm = msh.comm
    dx = ufl.dx(domain=msh, metadata=QDEG)
    d = (T_h - T_INF) - theta
    chi_in = _indicator(msh, lambda m: np.sum(m**2, axis=1) < 0.81 * B_**2)
    chi_sh = _dg0(msh, ct, 0.0, 1.0)
    out = [np.sqrt(_allsum(comm, d**2 * dx) / _allsum(comm, theta**2 * dx))]
    for w in (chi_in, chi_sh):
        out.append(np.sqrt(_allsum(comm, w * d**2 * dx) / _allsum(comm, w * theta**2 * dx)))
    return tuple(float(v) for v in out)


def _rel_diff(msh, a, b):
    """Whole-ball ‖(a − T∞) − (b − T∞)‖ / ‖b − T∞‖ (`MAT-9`'s form)."""
    comm = msh.comm
    dx = ufl.dx(domain=msh, metadata=QDEG)
    return float(np.sqrt(_allsum(comm, (a - b) ** 2 * dx) / _allsum(comm, (b - T_INF) ** 2 * dx)))


def _census(res):
    """Mesh, tags, cell count and volume deficits (printed)."""
    msh, ct, ft, t_build = _ball(res)
    comm = msh.comm
    nglob = msh.topology.index_map(3).size_global
    vol = _allsum(comm, fem.Constant(msh, default_scalar_type(1.0)) * ufl.dx(domain=msh))
    vc = _vcore(msh, ct)
    vs = vol - vc
    _say(
        f"MAT18T census res={res} cells={nglob} build={t_build:.1f}s "
        f"Vcore_h/Vcore-1={vc / V_CORE - 1:.4e} Vshell_h/Vshell-1={vs / V_SHELL - 1:.4e}"
    )
    return msh, ct, ft, nglob


def _run(msh, ft, *, k, rho_c, q, n, theta, startup, dt=None, snapshot_steps=()):
    t0 = time.perf_counter()
    s = solve_transient_conduction(
        msh,
        ft,
        k=k,
        rho_c=rho_c,
        q=q,
        h=H,
        t_inf=T_INF,
        t_init=T_INF,
        robin_tag=ROBIN_TAG,
        dt=T_STAR / n if dt is None else dt,
        n_steps=n,
        theta=theta,
        startup_half_steps=startup,
        degree=1,
        snapshot_steps=snapshot_steps,
    )
    return s, time.perf_counter() - t0


def _fine():
    """0.005 rung: guard, coefficients, comparand (cached)."""
    if "fine" not in _CACHE:
        msh, ct, ft, nglob = _census(0.005)
        _say(f"MAT18T cell guard: {nglob} cells at res=0.005 (guard {CELL_GUARD})")
        assert nglob <= CELL_GUARD, f"cell guard: {nglob} > {CELL_GUARD}; stop before the solve"
        k_f = _dg0(msh, ct, K1, K2)
        rc_f = _dg0(msh, ct, RC1, RC2)
        q_f = _dg0(msh, ct, Q1, Q2)
        theta, t_cmp, t_chunk = _exact_theta(msh, _two_layer, T_STAR)
        _say(
            f"MAT18T two-layer comparand at t*={T_STAR:.1f}s: {t_cmp:.1f}s "
            f"(max chunk call {t_chunk:.2f}s, predicted ≈ 0.5s)"
        )
        _CACHE["fine"] = (msh, ct, ft, nglob, k_f, rc_f, q_f, theta)
    return _CACHE["fine"]


def test_t1_cn_two_layer_t2_balance_controls():
    comm = MPI.COMM_WORLD
    msh, ct, ft, nglob, k_f, rc_f, q_f, theta = _fine()
    dt = T_STAR / N_CN

    # (T1) CN + Rannacher, DG0 k and ρc.
    s, t_solve = _run(
        msh, ft, k=k_f, rho_c=rc_f, q=q_f, n=N_CN, theta=0.5, startup=2, snapshot_steps=(39, 40)
    )
    e_all, e_in, e_sh = _errors(msh, ct, s.temperature, theta)
    vc = _vcore(msh, ct)
    _say(
        f"MAT18T (T1) CN+Rannacher N={N_CN} res=0.005 cells={nglob}: "
        f"interior_rel_L2={e_in:.4e} (predicted 5e-4..2e-3, band ≤ 1e-2) "
        f"whole={e_all:.4e} shell={e_sh:.4e} solve={t_solve:.1f}s"
    )
    _say(
        f"MAT18T (T2) on (T1): max_balance={s.max_balance_residual:.3e} "
        f"(predicted ≤ 5e-14, band ≤ 1e-9); src={s.source_power:.6f}W "
        f"q1*Vcore_h={Q1 * vc:.6f}W factorizations={s.n_factorizations}"
    )

    pts = np.array([[0.0, 0.0, 0.0], [A_, 0.0, 0.0], [0.999 * B_, 0.0, 0.0]])
    vals, ok = evaluate_vector_field_parallel(s.temperature, pts)
    if comm.rank == 0:
        vv = np.ravel(vals[:, 0]) if vals.ndim == 2 else np.ravel(vals)
        ref = _two_layer(np.array([0.0, A_, 0.999 * B_]), T_STAR) - T_INF
        for lab, v, o, rv in zip(("0", "a", "0.999b"), vv, ok, ref):
            th = float(np.real(v)) - T_INF if o else float("nan")
            print(
                f"MAT18T theta_h(r={lab}, t*) = {th:.4f} K (valid={bool(o)}); "
                f"comparand {rv:.4f} K",
                flush=True,
            )

    # (c2) ρc₁ storage weight everywhere in the balance, CN step 39 → 40.
    T39, T40 = s.snapshots[39], s.snapshots[40]
    ds = ufl.Measure("ds", domain=msh, subdomain_data=ft)
    stor1 = _allsum(comm, RC1 * (T40 - T39) / dt * ufl.dx(domain=msh))
    conv = _allsum(comm, H * (0.5 * T40 + 0.5 * T39 - T_INF) * ds(ROBIN_TAG))
    bal_c2 = abs(stor1 + conv - s.source_power) / s.source_power
    _say(
        f"MAT18T (c2) uniform ρc in the balance: miss measured {bal_c2:.4e} vs "
        f"predicted ≈ 0.07 (0.03..0.15) (asserted > 1e-3)"
    )

    # (c1) uniform ρc = RC1 everywhere, DG0 k, same piecewise q.
    s1, t1 = _run(msh, ft, k=k_f, rho_c=RC1, q=q_f, n=N_CN, theta=0.5, startup=2)
    _, e1_in, _ = _errors(msh, ct, s1.temperature, theta)
    _say(
        f"MAT18T (c1) uniform ρc: interior miss measured {e1_in:.4e} vs predicted "
        f"≈ 0.05 (0.02..0.10) (asserted > 1e-2); max_balance={s1.max_balance_residual:.3e} "
        f"solve={t1:.1f}s"
    )

    # Printed only: CN without the Rannacher start.
    s0, t0 = _run(msh, ft, k=k_f, rho_c=rc_f, q=q_f, n=N_CN, theta=0.5, startup=0)
    _, e0_in, _ = _errors(msh, ct, s0.temperature, theta)
    _say(
        f"MAT18T CN no-startup N={N_CN} (printed only): interior={e0_in:.4e} "
        f"max_balance={s0.max_balance_residual:.3e} solve={t0:.1f}s"
    )

    assert s.max_balance_residual <= 1e-9, "(T2) on the (T1) run"
    assert s1.max_balance_residual <= 1e-9, "(T2) on the (c1) run"
    assert s0.max_balance_residual <= 1e-9, "(T2) on the no-startup run"
    assert e_in <= 1e-2, f"(T1) interior {e_in:.4e} > 1 %"
    assert e1_in > 1e-2, f"(c1) control did not fire: {e1_in:.4e}"
    assert bal_c2 > 1e-3, f"(c2) control did not fire: {bal_c2:.4e}"


def test_e_equal_properties():
    """(E) DG0 k ≡ K1, ρc ≡ RC1, q ≡ Q1 vs the homogeneous Robin ball (radius b)."""
    msh, ct, ft, nglob, *_ = _fine()
    k_f = _dg0(msh, ct, K1, K1)
    rc_f = _dg0(msh, ct, RC1, RC1)
    q_f = _dg0(msh, ct, Q1, Q1)
    s, t_solve = _run(msh, ft, k=k_f, rho_c=rc_f, q=q_f, n=N_CN, theta=0.5, startup=2)
    theta, t_cmp, _ = _exact_theta(msh, _homog, T_STAR)
    e_all, e_in, e_sh = _errors(msh, ct, s.temperature, theta)
    _say(
        f"MAT18T (E) equal properties res=0.005 CN+Rannacher N={N_CN}: interior={e_in:.4e} "
        f"whole={e_all:.4e} shell={e_sh:.4e} (predicted 5e-4..2e-3, band ≤ 1e-2); "
        f"max_balance={s.max_balance_residual:.3e} solve={t_solve:.1f}s comparand={t_cmp:.1f}s"
    )
    assert s.max_balance_residual <= 1e-9, "(T2) on the (E) run"
    assert e_in <= 1e-2, f"(E) interior {e_in:.4e} > 1 %"


def test_d_coefficient_path():
    """(D) DG0 k ≡ K1, ρc ≡ RC1 vs scalar K1 / RC1, piecewise q, BE N 10, 0.01."""
    msh, ct, ft, _ = _census(0.01)
    q_f = _dg0(msh, ct, Q1, Q2)
    k_f = _dg0(msh, ct, K1, K1)
    rc_f = _dg0(msh, ct, RC1, RC1)
    s_f, _ = _run(msh, ft, k=k_f, rho_c=rc_f, q=q_f, n=10, theta=1.0, startup=0)
    s_c, _ = _run(msh, ft, k=K1, rho_c=RC1, q=q_f, n=10, theta=1.0, startup=0)
    d = _rel_diff(msh, s_f.temperature, s_c.temperature)
    _say(
        f"MAT18T (D) DG0-vs-scalar BE N=10 whole rel diff = {d:.3e} "
        f"(predicted ≤ 1e-15, band ≤ 1e-13); balances {s_f.max_balance_residual:.3e} / "
        f"{s_c.max_balance_residual:.3e}"
    )
    assert s_f.max_balance_residual <= 1e-9 and s_c.max_balance_residual <= 1e-9, "(T2) on (D)"
    assert d <= 1e-13, f"(D) {d:.3e} > 1e-13"


def test_t3_long_be_equals_steady():
    """(T3) BE 20 × 5 t_R with piecewise k, ρc, q equals the steady solve, 0.01."""
    msh, ct, ft, _ = _census(0.01)
    k_f = _dg0(msh, ct, K1, K2)
    rc_f = _dg0(msh, ct, RC1, RC2)
    q_f = _dg0(msh, ct, Q1, Q2)
    s, _ = _run(msh, ft, k=k_f, rho_c=rc_f, q=q_f, n=20, theta=1.0, startup=0, dt=5.0 * T_R)
    st = solve_steady_conduction(
        msh, ft, k=k_f, q=q_f, h=H, t_inf=T_INF, robin_tag=ROBIN_TAG, degree=1
    )
    diff = _rel_diff(msh, s.temperature, st.temperature)
    _say(
        f"MAT18T (T3) BE 20 x 5 t_R vs steady: rel diff={diff:.3e} "
        f"(predicted ≈ 1e-14, band ≤ 1e-8) max_balance={s.max_balance_residual:.3e}"
    )
    assert s.max_balance_residual <= 1e-9, "(T2) on (T3)"
    assert diff <= 1e-8, f"(T3) {diff:.3e} > 1e-8"


def test_cn_self_difference_printed():
    """Printed: CN + Rannacher self-difference at 0.01, N 20 / 40 / 80."""
    msh, ct, ft, _ = _census(0.01)
    k_f = _dg0(msh, ct, K1, K2)
    rc_f = _dg0(msh, ct, RC1, RC2)
    q_f = _dg0(msh, ct, Q1, Q2)
    sols, worst = [], 0.0
    for n in (20, 40, 80):
        s, _ = _run(msh, ft, k=k_f, rho_c=rc_f, q=q_f, n=n, theta=0.5, startup=2)
        worst = max(worst, s.max_balance_residual)
        sols.append(s.temperature)
    d1 = _rel_diff(msh, sols[0], sols[1])
    d2 = _rel_diff(msh, sols[1], sols[2])
    _say(
        f"MAT18T CN+Rannacher self-difference res=0.01: d(20,40)={d1:.4e} "
        f"d(40,80)={d2:.4e} rate={np.log(d1 / d2) / np.log(2.0):.3f} "
        f"(printed; predicted ≈ 2.0) worst balance={worst:.3e}"
    )
    assert worst <= 1e-9, "(T2) on the CN ladder"
