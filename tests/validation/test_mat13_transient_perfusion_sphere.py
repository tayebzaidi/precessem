"""`MAT-13` — the Pennes perfusion sink in the θ-scheme (§9 item 108).

Transient Pennes on `MAT-9` step 2's Robin ball (R 0.1, k 0.6, h 10, q 2000,
ρc 4.18e6, T∞ = T(0) = 20, t_R 69 667 s, t* = 0.2 t_R) with `MAT-11`'s
w = 240 W/m³/K (x = mR = 2), arterial = ambient.  Closed form:
`utils.analytical.robin_sphere_perfusion_transient_temperature` — θ = θ_ss^w
− Σ Bₙ^w e^{−(λₙ² + x²)Fo} Xₙ, Bₙ^w = (qR²/k) Cₙ/(λₙ² + x²) (Green's second
identity on θ_ss^w and Xₙ; no textbook equation pinned).

Numpy anchors (asserted): (N1) Bₙ^w vs quad; (N2) θ(r, 0) ≈ 0, 400 vs 800
terms, Fo 50 → θ_ss^w; (N3) w → 0 against `MAT-9`'s series; (N4) the
time-convolution identity θ^w = e^{−βt}θ⁰ + β∫e^{−βτ}θ⁰dτ, β = w/ρc, which
sees the exponent N3 cannot.
FEM anchors (asserted): (P_w) CN + Rannacher N 40, 0.01 m: interior error
≤ 1 %; (T) CN LSQ rate ≥ 1.8 over N 10…160; (E) per-step balance ≤ 1e-9,
recomputed from snapshots 39, 40; (∞) BE 20 × 5 t_R equals the perfused
steady solve to ≤ 1e-8; (S) one BE step at T∞ = 0 is ρcθ¹/Δt₁ = q to 1e-3;
(w→0) the unperfused run reproduces `MAT-9`'s 1.7085e-03.
Negative controls (asserted only as missing their anchor's bar): (a) w
dropped, (b) sign flipped, (c) perfused term dropped from (E).
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
from fem_em_solver.thermal.transient import _solve_transient
from fem_em_solver.utils.analytical import (
    robin_sphere_eigenvalues,
    robin_sphere_perfusion_temperature,
    robin_sphere_perfusion_transient_coefficients,
    robin_sphere_perfusion_transient_temperature,
    robin_sphere_transient_coefficients,
    robin_sphere_transient_temperature,
)

# `MAT-9` step 2's fixture (test_mat9_transient_conduction_sphere.py:41–50).
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
RES = 0.01

W = 240.0  # `MAT-11`'s w (test_mat11_perfusion_sphere.py:55): x = mR = 2
X2 = W * R**2 / K
# Bands restated by value, not widened: INTERIOR_BAND / BALANCE_BAND
# (test_mat11_perfusion_sphere.py:59,61), S1_BAND (test_mat12…:74).
INTERIOR_BAND = 1e-2
BALANCE_BAND = 1e-9
S1_BAND = 1e-3
# (w→0): `MAT-9`'s own series error on this rung (…MAT-12-W2-ruleC.log:127).
MAT9_RECORD = "1.7085e-03"

_MESHES: dict = {}


def _ball(res=RES):
    if res not in _MESHES:
        _MESHES[res] = MeshGenerator.ball_domain(radius=R, resolution=res)
    return _MESHES[res]


def _say(msg):
    if MPI.COMM_WORLD.rank == 0:
        print(msg, flush=True)


def _allsum(comm, form):
    return comm.allreduce(fem.assemble_scalar(fem.form(form)), op=MPI.SUM)


def _series_w(r, t, w=W, n_terms=400):
    return robin_sphere_perfusion_transient_temperature(
        r,
        t,
        radius=R,
        k=K,
        rho_c=RHO_C,
        q=Q,
        w=w,
        h=H,
        t_inf=T_INF,
        t_init=T_INF,
        n_terms=n_terms,
    )


def _series_0(r, t, n_terms=400):
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


def _theta_ss_w(r):
    return robin_sphere_perfusion_temperature(r, radius=R, k=K, q=Q, w=W, h=H)


def _interior_indicator(msh):
    tdim = msh.topology.dim
    imap = msh.topology.index_map(tdim)
    n = imap.size_local + imap.num_ghosts
    mids = dmesh.compute_midpoints(msh, tdim, np.arange(n, dtype=np.int32))
    Wdg = fem.functionspace(msh, ("DG", 0))
    chi = fem.Function(Wdg)
    chi.x.array[:] = (np.sum(mids**2, axis=1) < 0.81 * R**2).astype(float)
    chi.x.scatter_forward()
    return chi


def _exact_theta(msh, series):
    """θ(t*) from a series, interpolated in numpy into P2 (chunked)."""
    V2 = fem.functionspace(msh, ("Lagrange", 2))
    f = fem.Function(V2)

    def ev(x):
        r = np.sqrt(np.sum(x**2, axis=0))
        out = np.empty(r.shape)
        for s in range(0, r.size, 20000):
            out[s : s + 20000] = series(r[s : s + 20000], T_END) - T_INF
        return out

    f.interpolate(ev)
    return f


def _interior_error(msh, T_h, theta_ex):
    comm = msh.comm
    chi = _interior_indicator(msh)
    dx = ufl.dx(domain=msh, metadata=QDEG)
    e = _allsum(comm, chi * (T_h - T_INF - theta_ex) ** 2 * dx)
    n = _allsum(comm, chi * theta_ex**2 * dx)
    return float(np.sqrt(np.real(e) / np.real(n)))


def _rel_diff(msh, a, b):
    comm = msh.comm
    dx = ufl.dx(domain=msh, metadata=QDEG)
    e = _allsum(comm, (a - b) ** 2 * dx)
    n = _allsum(comm, (b - T_INF) ** 2 * dx)
    return float(np.sqrt(np.real(e) / np.real(n)))


def _transient(msh, ft, n, *, w=W, theta=0.5, startup=2, snapshot_steps=(), solver=None):
    solver = solver or solve_transient_conduction
    return solver(
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
        perfusion_w_per_m3_k=w,
    )


# ---------------------------------------------------------------- numpy anchors


def test_n1_coefficients_vs_quad():
    lam = robin_sphere_eigenvalues(BI, 20)
    C, B = robin_sphere_perfusion_transient_coefficients(lam, radius=R, k=K, q=Q, w=W, h=H)
    _, B0 = robin_sphere_transient_coefficients(lam, radius=R, k=K, q=Q, h=H)
    worst = 0.0
    for i, l in enumerate(lam):
        X = lambda p, l=l: np.sin(l * p) / (l * p) if p > 0 else 1.0
        nrm = quad(lambda p: p**2 * X(p) ** 2, 0, 1, epsabs=1e-14, epsrel=1e-13, limit=200)[0]
        b = (
            quad(
                lambda p: p**2 * float(_theta_ss_w(p * R)) * X(p),
                0,
                1,
                epsabs=1e-14,
                epsrel=1e-13,
                limit=200,
            )[0]
            / nrm
        )
        worst = max(worst, abs(b - B[i]))
    # MAT-12's Bₙ (x = 0) in place of Bₙ^w: the first wrong way to code it.
    wrong = float(np.max(np.abs(B0 - B)))
    _say(
        f"MAT13 (N1) max |Bn^w - quad| (n <= 20) = {worst:.3e} (bar 1e-10); "
        f"|Bn(x=0) - Bn^w| max = {wrong:.3e} (printed)"
    )
    assert worst <= 1e-10
    assert wrong > 1e-3


def test_n2_series_limits():
    r = np.linspace(0.0, R, 201)
    t0 = float(np.max(np.abs(_series_w(r, 0.0, n_terms=4000) - T_INF)))
    alpha = K / RHO_C
    fos = np.array([0.05, 0.1, 0.2, 1.0])
    d48 = max(
        float(
            np.max(
                np.abs(
                    _series_w(r, fo * R**2 / alpha, n_terms=400)
                    - _series_w(r, fo * R**2 / alpha, n_terms=800)
                )
            )
        )
        for fo in fos
    )
    ss = float(np.max(np.abs(_series_w(r, 50 * R**2 / alpha) - T_INF - _theta_ss_w(r))))
    bar0 = 1e-6 * float(_theta_ss_w(0.0))
    _say(
        f"MAT13 (N2) max|theta(r,0)| (4000 terms) = {t0:.3e} (bar {bar0:.3e}); "
        f"400 vs 800 at Fo >= 0.05 = {d48:.3e}; Fo 50 vs theta_ss^w = {ss:.3e}"
    )
    assert t0 <= bar0
    assert d48 <= 1e-12
    assert ss <= 1e-12


def test_n3_small_w_limit():
    rs = np.array([0.0, 0.5 * R, R])
    a = _series_w(rs, T_END, w=6e-3)
    b = _series_0(rs, T_END)
    d = float(np.max(np.abs(a - b)))
    _say(
        f"MAT13 (N3) w = 6e-3 vs MAT-9 series at t*, r in (0, R/2, R): max diff = "
        f"{d:.3e} K (bar 1e-4; hand 9.3e-6)"
    )
    assert d <= 1e-4


def test_n4_time_convolution_identity():
    # Root the 400 eigenvalues once (the public function re-roots per call).
    lam = robin_sphere_eigenvalues(BI, 400)
    _, B0 = robin_sphere_transient_coefficients(lam, radius=R, k=K, q=Q, h=H)
    from fem_em_solver.utils.analytical import robin_sphere_steady_temperature

    alpha = K / RHO_C
    beta = W / RHO_C
    worst = 0.0
    lines = []
    for r in (0.0, 0.5 * R, R):
        rho = r / R
        x = lam * rho
        X = np.where(x == 0.0, 1.0, np.sin(x) / np.where(x == 0.0, 1.0, x))
        th_ss0 = float(robin_sphere_steady_temperature(r, radius=R, k=K, q=Q, h=H))

        def theta0(tau, X=X, th_ss0=th_ss0):
            return th_ss0 - float(np.sum(B0 * np.exp(-(lam**2) * alpha * tau / R**2) * X))

        integral = quad(
            lambda tau: np.exp(-beta * tau) * theta0(tau),
            0.0,
            T_END,
            epsabs=1e-14,
            epsrel=1e-13,
            limit=200,
        )[0]
        conv = np.exp(-beta * T_END) * theta0(T_END) + beta * integral
        direct = float(_series_w(r, T_END)) - T_INF
        worst = max(worst, abs(conv - direct))
        lines.append(f"r={r:.3f}: {direct:.10f} vs {conv:.10f}")
    _say(
        f"MAT13 (N4) convolution identity at t*: {'; '.join(lines)}; max diff = "
        f"{worst:.3e} K (bar 1e-10)"
    )
    assert worst <= 1e-10


# ------------------------------------------------------------------ FEM anchors


def test_fem_pw_e_controls():
    comm = MPI.COMM_WORLD
    msh, _, ft = _ball()
    n_glob = msh.topology.index_map(3).size_global
    t0 = time.perf_counter()
    s = _transient(msh, ft, 40, snapshot_steps=(39, 40))
    t_solve = time.perf_counter() - t0
    theta_ex = _exact_theta(msh, _series_w)
    e_pw = _interior_error(msh, s.temperature, theta_ex)
    pts = np.array([0.0, R])
    pred = _series_w(pts, T_END) - T_INF
    _say(
        f"MAT13 (P_w) CN+Rannacher N=40 res={RES} cells={n_glob}: interior_rel_L2="
        f"{e_pw:.4e} (predicted ≈ 1.7e-3 - 1.8e-3) max_balance={s.max_balance_residual:.3e} "
        f"perfused_power={s.perfused_power:.4e} W source={s.source_power:.4e} W "
        f"solve {t_solve:.1f}s; series theta(0,t*)={pred[0]:.4f}, theta(R,t*)={pred[1]:.4f} "
        f"(hand 4.1988 / 2.6009)"
    )

    # (E) recomputed from snapshots 39, 40 (CN step), and control (c).
    T39, T40 = s.snapshots[39], s.snapshots[40]
    ds = ufl.Measure("ds", domain=msh, subdomain_data=ft)
    dt = T_END / 40
    stor = float(np.real(_allsum(comm, RHO_C * (T40 - T39) / dt * ufl.dx)))
    conv = float(np.real(_allsum(comm, H * (0.5 * T40 + 0.5 * T39 - T_INF) * ds(ROBIN_TAG))))
    perf = float(np.real(_allsum(comm, W * (0.5 * T40 + 0.5 * T39 - T_INF) * ufl.dx)))
    src = s.source_power
    bal = abs(stor + conv + perf - src) / src
    bal_c = abs(stor + conv - src) / src
    _say(
        f"MAT13 (E) recomputed last-step balance = {bal:.3e}; split storage/convected/"
        f"perfused = {stor / src:.4f}/{conv / src:.4f}/{perf / src:.4f} (predicted "
        f"0.21/0.39/0.40); perf vs result.perfused_power {perf:.6e}/{s.perfused_power:.6e}"
    )
    _say(
        f"MAT13 control (c) perfused term dropped from (E): {bal_c:.4e} (asserted > 1e-3; "
        f"predicted ≈ 0.398)"
    )

    # (a) w dropped, and (w→0): the same run is MAT-9's.
    s0 = _transient(msh, ft, 40, w=0.0)
    e_a = _interior_error(msh, s0.temperature, theta_ex)
    e_mat9 = _interior_error(msh, s0.temperature, _exact_theta(msh, _series_0))
    _say(
        f"MAT13 control (a) w dropped: interior_rel_L2 vs perfused series = {e_a:.4e} "
        f"(asserted > 1e-2; predicted ≈ 0.39); (w->0) vs MAT-9 series = {e_mat9:.4e} "
        f"(record {MAT9_RECORD}), perfused_power = {s0.perfused_power!r}"
    )

    # (b) sign flipped (+wθ), through the private body only.
    sb = _transient(msh, ft, 40, w=-W, solver=_solve_transient)
    e_b = _interior_error(msh, sb.temperature, theta_ex)
    _say(
        f"MAT13 control (b) sign flipped: interior_rel_L2 = {e_b:.4e} (asserted > 1e-2; "
        f"predicted ≈ 1.03)"
    )

    # Quasi-static θ_ss^w for θ(t*): printed.
    st = solve_steady_conduction(
        msh, ft, k=K, q=Q, h=H, t_inf=T_INF, robin_tag=ROBIN_TAG, degree=1, perfusion_w_per_m3_k=W
    )
    e_qs = _interior_error(msh, st.temperature, theta_ex)
    _say(f"MAT13 [printed] quasi-static theta_ss^w for theta(t*): {e_qs:.4e} (predicted ≈ 0.28)")

    assert s.max_balance_residual <= BALANCE_BAND, "(E) on the (P_w) run"
    assert bal <= BALANCE_BAND, f"(E) recomputed {bal:.3e}"
    assert e_pw <= INTERIOR_BAND, f"(P_w) {e_pw:.4e} > 1 %"
    assert f"{e_mat9:.4e}" == MAT9_RECORD, f"(w->0) {e_mat9:.4e} != {MAT9_RECORD}"
    assert s0.perfused_power == 0.0
    assert e_a > INTERIOR_BAND, "(a) w dropped still inside (P_w)'s band"
    assert e_b > INTERIOR_BAND, "(b) sign flipped still inside (P_w)'s band"
    assert bal_c > 1e-3, "(c) balance does not see the perfused term"


def test_fem_t_rate_and_e_on_every_run():
    msh, _, ft = _ball()
    sols, worst = [], 0.0
    t0 = time.perf_counter()
    for n in LADDER_N:
        s = _transient(msh, ft, n)
        worst = max(worst, s.max_balance_residual)
        sols.append(s.temperature)
    d = [_rel_diff(msh, sols[i], sols[i + 1]) for i in range(len(sols) - 1)]
    dts = T_END / np.asarray(LADDER_N[: len(d)], dtype=float)
    lsq = float(np.polyfit(np.log(dts), np.log(d), 1)[0])
    pair = [float(np.log(d[i] / d[i + 1]) / np.log(2.0)) for i in range(len(d) - 1)]
    _say(
        f"MAT13 (T) CN+Rannacher d_N (N=10,20,40,80) = {[f'{x:.4e}' for x in d]} LSQ rate="
        f"{lsq:.3f} (predicted 2.00) pairwise={[f'{p:.3f}' for p in pair]}; (E) worst "
        f"per-step balance {worst:.3e} ({time.perf_counter() - t0:.1f}s)"
    )
    assert worst <= BALANCE_BAND, "(E) over the ladder"
    assert lsq >= 1.8, f"(T) CN LSQ rate {lsq:.3f}"


def test_fem_inf_long_be_equals_perfused_steady():
    msh, _, ft = _ball()
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
        perfusion_w_per_m3_k=W,
    )
    st = solve_steady_conduction(
        msh, ft, k=K, q=Q, h=H, t_inf=T_INF, robin_tag=ROBIN_TAG, degree=1, perfusion_w_per_m3_k=W
    )
    diff = _rel_diff(msh, s.temperature, st.temperature)
    _say(
        f"MAT13 (inf) BE 20 x 5 t_R vs perfused steady: rel diff = {diff:.3e} "
        f"max_balance = {s.max_balance_residual:.3e}; perfused power transient/steady = "
        f"{s.perfused_power:.8e}/{st.perfused_power:.8e}"
    )
    assert s.max_balance_residual <= BALANCE_BAND
    assert diff <= 1e-8


def test_fem_s_first_step_is_q_over_rho_c():
    comm = MPI.COMM_WORLD
    msh, _, ft = _ball()
    dt1 = 1e-5 * T_R
    s = solve_transient_conduction(
        msh,
        ft,
        k=K,
        rho_c=RHO_C,
        q=Q,
        h=H,
        t_inf=0.0,
        t_init=0.0,
        robin_tag=ROBIN_TAG,
        dt=dt1,
        n_steps=1,
        theta=1.0,
        perfusion_w_per_m3_k=W,
    )
    chi = _interior_indicator(msh)
    T1 = s.temperature
    dx = ufl.dx(domain=msh, metadata=QDEG)
    e = _allsum(comm, chi * (RHO_C * T1 / dt1 - Q) ** 2 * dx)
    n = _allsum(comm, chi * Q**2 * dx)
    s1 = float(np.sqrt(np.real(e) / np.real(n)))
    vol = float(np.real(_allsum(comm, 1.0 * dx)))
    mean_def = 1.0 - float(np.real(_allsum(comm, RHO_C * T1 / dt1 * dx))) / (Q * vol)
    _say(
        f"MAT13 (S) one BE step dt1 = 1e-5 t_R at T_inf = 0: interior ||rho_c theta1/dt1 - q||"
        f"/||q|| = {s1:.3e} (bar 1e-3; predicted ≈ 0.6e-4 - 0.9e-4); volume-mean deficit "
        f"{mean_def:.3e} (printed; predicted ≈ 9.0e-5)"
    )
    assert s.max_balance_residual <= BALANCE_BAND
    assert s1 <= S1_BAND
