"""`MAT-9` step 2 — transient conduction in the Robin-cooled ball (§9 item 72).

Fixture: step 1's ball (R 0.1 m, k 0.6, Q 2000 W/m³, h 10, Bi 5/3, T∞ 20 °C)
plus ρc 4.18e6 J/m³/K, so t_R = R²ρc/k = 69 667 s.  T(0) = T∞, Q switched on
at t = 0, t_end = 0.2 t_R (Fo 0.2).  Errors on θ = T − T∞, interior mask
|x|² < 0.81R² (step 1's).  Closed form:
`utils.analytical.robin_sphere_transient_temperature` (series kept out of
UFL: evaluated in numpy at P2 dof coordinates).

Numpy anchors (asserted): eigenvalues (Bi = 1 → (n − ½)π; Bi = 5/3 residual),
Bₙ/Cₙ vs quad, t = 0 series ≈ 0, 400 vs 800 terms, Fo 50 → steady.
FEM anchors (asserted): (i) CN + Rannacher, N 40, 0.005 m: interior error
≤ 1 %; (ii) 0.01 m, N 10…160 (step ratio 2), whole-ball self-differences:
LSQ rate BE ∈ [0.9, 1.2], CN ≥ 1.8; (iii) per-step balance ≤ 1e-9; (iv) BE
20 × 5 t_R equals the steady solve to ≤ 1e-8.  Negative controls (each
asserts only against its own anchor's bar, magnitude *predicted*):
(a) quasi-static on (i), > 1 % (pred. ≈ 90 %); (b) storage dropped from
(iii), > 1e-3 (pred. ≈ 0.47); (c) BE rate against CN's bar, < 1.8 (pred. 1.0).
"""

from __future__ import annotations

import time

import numpy as np
import ufl
from dolfinx import fem
from dolfinx import mesh as dmesh
from mpi4py import MPI
from scipy.integrate import quad

from fem_em_solver.io.mesh import MeshGenerator
from fem_em_solver.thermal import solve_steady_conduction, solve_transient_conduction
from fem_em_solver.utils.analytical import (
    robin_sphere_eigenvalues,
    robin_sphere_steady_temperature,
    robin_sphere_transient_coefficients,
    robin_sphere_transient_temperature,
)

R = 0.1
K = 0.6
Q = 2000.0
H = 10.0
T_INF = 20.0
RHO_C = 4.18e6
BI = H * R / K
T_R = R**2 * RHO_C / K
T_END = 0.2 * T_R
ROBIN_TAG = 1
QDEG = {"quadrature_degree": 4}
LADDER_N = (10, 20, 40, 80, 160)
THETA_SS0 = Q * R / (3 * H) + Q * R**2 / (6 * K)

_MESHES: dict = {}
_CACHE: dict = {}


def _ball(res):
    if res not in _MESHES:
        _MESHES[res] = MeshGenerator.ball_domain(radius=R, resolution=res)
    return _MESHES[res]


def _say(comm, msg):
    if comm.rank == 0:
        print(msg, flush=True)


def _allsum(comm, form):
    return comm.allreduce(fem.assemble_scalar(fem.form(form)), op=MPI.SUM)


def _series(r, t, n_terms=400):
    return robin_sphere_transient_temperature(
        r,
        t,
        radius=R,
        k=K,
        rho_c=RHO_C,
        q=Q,
        h=H,
        t_inf=T_INF,
        t_init=T_INF,
        n_terms=n_terms,
    )


def _interior_indicator(msh):
    tdim = msh.topology.dim
    imap = msh.topology.index_map(tdim)
    n = imap.size_local + imap.num_ghosts
    mids = dmesh.compute_midpoints(msh, tdim, np.arange(n, dtype=np.int32))
    W = fem.functionspace(msh, ("DG", 0))
    chi = fem.Function(W)
    chi.x.array[:] = (np.sum(mids**2, axis=1) < 0.81 * R**2).astype(float)
    chi.x.scatter_forward()
    return chi


def _exact_theta(msh, t):
    """θ(t) from the series, interpolated in numpy into P2 (chunked)."""
    V2 = fem.functionspace(msh, ("Lagrange", 2))
    f = fem.Function(V2)

    def ev(x):
        r = np.sqrt(np.sum(x**2, axis=0))
        out = np.empty(r.shape)
        for s in range(0, r.size, 20000):
            out[s : s + 20000] = _series(r[s : s + 20000], t) - T_INF
        return out

    f.interpolate(ev)
    return f


def _interior_error(msh, T_h, theta_ex):
    comm = msh.comm
    chi = _interior_indicator(msh)
    dx = ufl.dx(domain=msh, metadata=QDEG)
    e = _allsum(comm, chi * (T_h - T_INF - theta_ex) ** 2 * dx)
    n = _allsum(comm, chi * theta_ex**2 * dx)
    return float(np.sqrt(e / n))


def _rel_diff(msh, a, b):
    """Whole-ball ‖(a − T∞) − (b − T∞)‖ / ‖b − T∞‖."""
    comm = msh.comm
    dx = ufl.dx(domain=msh, metadata=QDEG)
    e = _allsum(comm, (a - b) ** 2 * dx)
    n = _allsum(comm, (b - T_INF) ** 2 * dx)
    return float(np.sqrt(e / n))


def _transient(msh, ft, n, theta, startup, snapshot_steps=()):
    return solve_transient_conduction(
        msh,
        ft,
        k=K,
        rho_c=RHO_C,
        q=Q,
        h=H,
        t_inf=T_INF,
        t_init=T_INF,
        robin_tag=ROBIN_TAG,
        dt=T_END / n,
        n_steps=n,
        theta=theta,
        startup_half_steps=startup,
        degree=1,
        snapshot_steps=snapshot_steps,
    )


def _rates(d):
    """(LSQ, pairwise) rates of d_N against Δt, step ratio 2."""
    d = np.asarray(d)
    dts = T_END / np.asarray(LADDER_N[: len(d)], dtype=float)
    lsq = float(np.polyfit(np.log(dts), np.log(d), 1)[0])
    pair = [float(np.log(d[i] / d[i + 1]) / np.log(2.0)) for i in range(len(d) - 1)]
    return lsq, pair


def _ladder(theta, startup, label):
    """Whole-ball self-differences on the 0.01 m mesh, N 10…160."""
    key = (theta, startup)
    if key in _CACHE:
        return _CACHE[key]
    comm = MPI.COMM_WORLD
    msh, _, ft = _ball(0.01)
    theta_ex = _exact_theta(msh, T_END)
    sols, worst = [], 0.0
    t0 = time.perf_counter()
    for n in LADDER_N:
        s = _transient(msh, ft, n, theta, startup)
        worst = max(worst, s.max_balance_residual)
        e_ser = _interior_error(msh, s.temperature, theta_ex)
        _say(
            comm,
            f"MAT9T {label} N={n} dt={T_END / n:.2f}s series_interior_rel_L2="
            f"{e_ser:.4e} max_balance={s.max_balance_residual:.3e} "
            f"factorizations={s.n_factorizations}",
        )
        sols.append(s.temperature)
    d = [_rel_diff(msh, sols[i], sols[i + 1]) for i in range(len(sols) - 1)]
    lsq, pair = _rates(d)
    _say(
        comm,
        f"MAT9T {label} d_N (N=10,20,40,80; whole ball, ratio 2) = "
        f"{[f'{x:.4e}' for x in d]}  LSQ rate={lsq:.3f} pairwise="
        f"{[f'{p:.3f}' for p in pair]} ({time.perf_counter() - t0:.1f}s)",
    )
    _CACHE[key] = (d, lsq, pair, worst)
    return _CACHE[key]


# ---------------------------------------------------------------- numpy anchors


def test_numpy_eigenvalues():
    lam1 = robin_sphere_eigenvalues(1.0, 50)
    n = np.arange(1, 51)
    err = np.max(np.abs(lam1 - (n - 0.5) * np.pi))
    # Bi = 5/3 residual asserted for n ≤ 20 (the range the coefficient anchor
    # uses).  Measured 2026-09-29 (MAT-9-step2-dev, 11:20Z): over n ≤ 50 the
    # max residual read 1.529e-12, which is the double-precision floor
    # |f'(λ)|·ulp(λ) ≈ λ·ulp(λ) ≈ 157 × 2.8e-14 at n = 50, not a root error;
    # the n ≤ 50 value and that floor are printed, not asserted.
    lam50 = robin_sphere_eigenvalues(BI, 50)
    res50 = np.abs((1 - BI) * np.sin(lam50) - lam50 * np.cos(lam50))
    floor50 = float(np.max(lam50 * np.spacing(lam50)))
    lam = lam50[:20]
    resid = float(np.max(res50[:20]))
    _say(
        MPI.COMM_WORLD,
        f"MAT9T numpy: Bi=1 max|λn-(n-½)π| (n≤50)={err:.3e}; "
        f"Bi=5/3 max root residual n≤20={resid:.3e}, n≤50={np.max(res50):.3e} "
        f"(printed; float floor λ·ulp(λ) at n≤50 = {floor50:.3e}); "
        f"λ1={lam[0]:.6f} (predicted ≈ 1.907)",
    )
    assert err <= 1e-12
    assert resid <= 1e-12

    # Scale-correct anchor over the series' whole range (§9 item 76, 2026-09-29
    # 10:30 review ruling): |f(λ)| has a float floor ∝ λ²·eps, so the gated
    # quantity is the implied relative root error e_n = |f(λn)|/(|f'(λn)|·λn),
    # f(λ) = (1 − Bi) sin λ − λ cos λ, f'(λ) = −Bi cos λ + λ sin λ, over the
    # n ≤ 4000 terms the θ(r, 0) anchor sums.  Predicted ≲ 1e-15 (brentq
    # rtol = 4·eps).
    def rel_root_err(x):
        f = (1 - BI) * np.sin(x) - x * np.cos(x)
        fp = -BI * np.cos(x) + x * np.sin(x)
        return np.abs(f) / (np.abs(fp) * x)

    nn = np.arange(1, 4001)
    lam4k = robin_sphere_eigenvalues(BI, 4000)
    e = rel_root_err(lam4k)
    bracketed = bool(np.all((lam4k > (nn - 1) * np.pi) & (lam4k < nn * np.pi)))
    # Negative control (asserted only as > 1e-8): the Bi = 1 roots (n − ½)π
    # sit ≈ (Bi − 1)/λ off the Bi = 5/3 roots; predicted e_n ≈ (2/3)/λ²
    # (≈ 2.7e-5 at n = 50, ≈ 4e-9 at n = 4000) and max e_n ≥ 1e-2 (n = 1).
    e_ctl = rel_root_err((nn - 0.5) * np.pi)
    _say(
        MPI.COMM_WORLD,
        f"MAT9T numpy: Bi=5/3 max e_n (n≤4000)={np.max(e):.3e} at "
        f"n={int(np.argmax(e)) + 1} (bar 1e-14, predicted ≲ 1e-15); bracketed={bracketed}; "
        f"control max e_n={np.max(e_ctl):.3e} at n={int(np.argmax(e_ctl)) + 1} "
        f"(predicted ≥ 1e-2), e_50={e_ctl[49]:.3e} (pred. ≈ 2.7e-5), "
        f"e_4000={e_ctl[-1]:.3e} (pred. ≈ 4e-9)",
    )
    assert np.max(e) <= 1e-14
    assert bracketed
    assert np.max(e_ctl) > 1e-8


def test_numpy_coefficients_vs_quad():
    lam = robin_sphere_eigenvalues(BI, 20)
    C, B = robin_sphere_transient_coefficients(lam, radius=R, k=K, q=Q, h=H)
    worst = 0.0
    for i, l in enumerate(lam):
        X = lambda p, l=l: np.sin(l * p) / (l * p) if p > 0 else 1.0
        nrm = quad(lambda p: p**2 * X(p) ** 2, 0, 1, epsabs=1e-14, epsrel=1e-13, limit=200)[0]
        c = quad(lambda p: p**2 * X(p), 0, 1, epsabs=1e-14, epsrel=1e-13, limit=200)[0] / nrm
        th = lambda p: robin_sphere_steady_temperature(p * R, radius=R, k=K, q=Q, h=H)
        b = (
            quad(lambda p: p**2 * th(p) * X(p), 0, 1, epsabs=1e-14, epsrel=1e-13, limit=200)[0]
            / nrm
        )
        worst = max(worst, abs(c - C[i]), abs(b - B[i]))
    _say(MPI.COMM_WORLD, f"MAT9T numpy: max |Bn,Cn - quad| (n ≤ 20) = {worst:.3e}")
    assert worst <= 1e-10


def test_numpy_series_limits():
    comm = MPI.COMM_WORLD
    r = np.linspace(0.0, R, 201)
    t0 = np.max(np.abs(_series(r, 0.0, n_terms=4000) - T_INF))
    alpha = K / RHO_C
    fos = np.array([0.05, 0.1, 0.2, 1.0])
    d48 = max(
        np.max(np.abs(_series(r, fo * R**2 / alpha, 400) - _series(r, fo * R**2 / alpha, 800)))
        for fo in fos
    )
    ss = np.max(
        np.abs(
            _series(r, 50 * R**2 / alpha)
            - robin_sphere_steady_temperature(r, radius=R, k=K, q=Q, h=H, t_inf=T_INF)
        )
    )
    _say(
        comm,
        f"MAT9T numpy: max|θ(r,0)| (4000 terms)={t0:.3e} "
        f"(bar {1e-6 * THETA_SS0:.3e}); 400 vs 800 at Fo≥0.05={d48:.3e}; "
        f"Fo 50 vs steady={ss:.3e}",
    )
    assert t0 <= 1e-6 * THETA_SS0
    assert d48 <= 1e-12
    assert ss <= 1e-12


# ------------------------------------------------------------------ FEM anchors


def test_fem_i_cn_accuracy_and_control_a_b():
    comm = MPI.COMM_WORLD
    msh, _, ft = _ball(0.005)
    n_glob = msh.topology.index_map(3).size_global
    t0 = time.perf_counter()
    s = _transient(msh, ft, 40, 0.5, 2, snapshot_steps=(39, 40))
    t_solve = time.perf_counter() - t0
    theta_ex = _exact_theta(msh, T_END)
    e_i = _interior_error(msh, s.temperature, theta_ex)
    _say(
        comm,
        f"MAT9T (i) CN+Rannacher N=40 res=0.005 cells={n_glob}: interior_rel_L2="
        f"{e_i:.4e} (predicted ≈ 1e-3) max_balance={s.max_balance_residual:.3e} "
        f"factorizations={s.n_factorizations} solve {t_solve:.1f}s",
    )
    assert s.max_balance_residual <= 1e-9, "(iii) on the (i) run"
    assert e_i <= 1e-2, f"(i) {e_i:.4e} > 1 %"

    # (a) quasi-static: the steady field stands in for T(t_end).
    st = solve_steady_conduction(msh, ft, k=K, q=Q, h=H, t_inf=T_INF, robin_tag=ROBIN_TAG, degree=1)
    e_a = _interior_error(msh, st.temperature, theta_ex)
    _say(
        comm,
        f"MAT9T control (a) quasi-static: interior_rel_L2={e_a:.4e} "
        "(asserted > 1e-2; predicted ≈ 0.90)",
    )
    assert e_a > 1e-2

    # (b) storage term dropped from (iii) at the end state (CN step 39→40).
    T39, T40 = s.snapshots[39], s.snapshots[40]
    ds = ufl.Measure("ds", domain=msh, subdomain_data=ft)
    conv = _allsum(comm, H * (0.5 * T40 + 0.5 * T39 - T_INF) * ds(ROBIN_TAG))
    bal_b = abs(float(conv) - s.source_power) / s.source_power
    _say(
        comm,
        f"MAT9T control (b) storage dropped: balance={bal_b:.4e} "
        "(asserted > 1e-3; predicted ≈ 0.47)",
    )
    assert bal_b > 1e-3


def test_fem_ii_rates_iii_balance_and_control_c():
    comm = MPI.COMM_WORLD
    d_be, lsq_be, _, w_be = _ladder(1.0, 0, "BE")
    d_cn, lsq_cn, _, w_cn = _ladder(0.5, 2, "CN+Rannacher")
    _ladder(0.5, 0, "CN-no-startup (printed only)")
    _say(comm, f"MAT9T (iii) worst per-step balance: BE {w_be:.3e}, CN {w_cn:.3e}")
    assert w_be <= 1e-9 and w_cn <= 1e-9, "(iii)"
    assert 0.9 <= lsq_be <= 1.2, f"(ii) BE LSQ rate {lsq_be:.3f}"
    assert lsq_cn >= 1.8, f"(ii) CN LSQ rate {lsq_cn:.3f}"
    _say(
        comm,
        f"MAT9T control (c) BE rate vs CN bar: {lsq_be:.3f} " "(asserted < 1.8; predicted 1.0)",
    )
    assert lsq_be < 1.8


def test_fem_iv_long_be_equals_steady():
    comm = MPI.COMM_WORLD
    msh, _, ft = _ball(0.01)
    s = solve_transient_conduction(
        msh,
        ft,
        k=K,
        rho_c=RHO_C,
        q=Q,
        h=H,
        t_inf=T_INF,
        t_init=T_INF,
        robin_tag=ROBIN_TAG,
        dt=5.0 * T_R,
        n_steps=20,
        theta=1.0,
    )
    st = solve_steady_conduction(msh, ft, k=K, q=Q, h=H, t_inf=T_INF, robin_tag=ROBIN_TAG, degree=1)
    diff = _rel_diff(msh, s.temperature, st.temperature)
    _say(
        comm,
        f"MAT9T (iv) BE 20 x 5 t_R vs steady: rel diff={diff:.3e} "
        f"(predicted ≈ 1e-13) max_balance={s.max_balance_residual:.3e}",
    )
    assert s.max_balance_residual <= 1e-9
    assert diff <= 1e-8
