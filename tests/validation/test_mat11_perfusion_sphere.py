"""`MAT-11` — the Pennes perfusion sink on `MAT-9`'s Robin ball (§9 item 97).

Anchor: −k∇²θ + wθ = q in r < R, −kθ′(R) = hθ(R), arterial = ambient
(Pennes 1948): θ(r) = (q/w)[1 − A·sinh(mr)/(mr)], m² = w/k, A = Bi·mR/D,
D = mR cosh mR + (Bi − 1) sinh mR
(`utils.analytical.robin_sphere_perfusion_temperature`; derivation in its
docstring, no textbook equation pinned).

Fixture: `MAT-9` step 1's by name (R 0.1, k 0.6, h 10, q 2000), with
**w = 240 W/m³/K** — mR = 2, 49 % of the power convected and 51 % perfused,
so k, h and w are all exercised.  A numerical choice, not a tissue claim.
w = 2 400 (mR = 6.32) is printed only.  **Errors are taken on θ, never on T.**

Numpy anchors (asserted): (N1) Robin condition, (N2) ODE residual, (N3) global
balance, (N4) w → 0 recovers `MAT-9`, (N5) mR = 40 limits.

FEM anchors (asserted, P1): (i) 0.01 m interior (|x|² < 0.81R²) relative L2
≤ 1 % (`MAT-9`'s (i) band); (ii) whole-ball LSQ rate on 0.02 / 0.01 / 0.005 m
≥ 1.7; (iii) |∮hθ_h dS + ∫wθ_h dV − ∫q dV|/∫q dV ≤ 1e-9 on every rung.

Negative controls (asserted only as "misses (i)'s 1 %"): (a) perfusion
dropped; (b) the sink's sign flipped — `MAT-5`'s growth form past runaway.
"""

from __future__ import annotations

import time

import numpy as np
import pytest
import ufl
from dolfinx import default_scalar_type, fem
from mpi4py import MPI

from fem_em_solver.thermal.steady import _solve, solve_steady_conduction
from fem_em_solver.utils.analytical import (
    robin_sphere_perfusion_coefficients,
    robin_sphere_perfusion_temperature,
    robin_sphere_steady_temperature,
)
from tests.validation.test_mat9_steady_conduction_sphere import (
    H,
    K,
    LADDER,
    Q,
    R,
    ROBIN_TAG,
    T_INF,
    _allsum,
    _ball,
    _interior_indicator,
    _say,
)

W = 240.0
W_STRONG = 2400.0
# `MAT-9` (i)'s band, a literal at test_mat9_steady_conduction_sphere.py:156
# (no module constant to import), restated by value.
INTERIOR_BAND = 1.0e-2
RATE_FLOOR = 1.7
BALANCE_BAND = 1.0e-9
QDEG = {"quadrature_degree": 6}
BI = H * R / K
V_BALL = 4.0 / 3.0 * np.pi * R**3
KW = dict(radius=R, k=K, q=Q, h=H)


def _theta(r, w):
    return robin_sphere_perfusion_temperature(r, w=w, **KW)


# ----------------------------------------------------------------- numpy ---


def test_n1_robin_condition():
    """(N1) −kθ′(R) = hθ(R), θ′ from the derivation, not from the function."""
    m, D, _A, _mean = robin_sphere_perfusion_coefficients(w=W, **KW)
    x = m * R
    c = BI * R * (Q / W) / D
    dtheta = -c * (x * np.cosh(x) - np.sinh(x)) / R**2
    lhs, rhs = -K * dtheta, H * float(_theta(R, W))
    print(f"\nMAT11 N1 -k theta'(R)={lhs:.12e} h theta(R)={rhs:.12e} rel={abs(lhs / rhs - 1):.3e}")
    assert abs(lhs / rhs - 1.0) <= 1e-12


def test_n2_ode_residual():
    """(N2) −k(θ″ + 2θ′/r) + wθ − q at r ∈ {0.2, 0.5, 0.8}R."""
    m, D, _A, _mean = robin_sphere_perfusion_coefficients(w=W, **KW)
    c = BI * R * (Q / W) / D
    worst = 0.0
    for frac in (0.2, 0.5, 0.8):
        r = frac * R
        sh, ch = np.sinh(m * r), np.cosh(m * r)
        du = (m * r * ch - sh) / r**2
        d2u = (m * m * r * r * sh - 2.0 * m * r * ch + 2.0 * sh) / r**3
        lap = -c * (d2u + 2.0 * du / r)
        res = -K * lap + W * float(_theta(r, W)) - Q
        worst = max(worst, abs(res) / Q)
    print(f"\nMAT11 N2 worst |residual|/q = {worst:.3e}")
    assert worst <= 1e-10


def test_n3_global_balance_and_predicted_values():
    """(N3) h·4πR²θ(R) + wVθ̄ = qV."""
    _m, _D, _A, mean = robin_sphere_perfusion_coefficients(w=W, **KW)
    conv = H * 4.0 * np.pi * R**2 * float(_theta(R, W))
    perf = W * V_BALL * mean
    rel = abs((conv + perf) / (Q * V_BALL) - 1.0)
    print(
        f"\nMAT11 N3 theta(0)={float(_theta(0.0, W)):.4f} theta(R)={float(_theta(R, W)):.4f} "
        f"mean={mean:.4f} K (predicted 5.5394 / 3.2668 / 4.2498); "
        f"convected={conv:.4f} W perfused={perf:.4f} W total={conv + perf:.4f} W "
        f"(predicted 4.1052 + 4.2724 = 8.3776); rel={rel:.3e}"
    )
    assert rel <= 1e-12


def test_n4_small_w_recovers_mat9():
    """(N4) w = 6e-3: |θ/θ_MAT9 − 1| ≤ 1e-4 at r ∈ {0, R/2, R}."""
    worst = 0.0
    for r in (0.0, 0.5 * R, R):
        ref = float(robin_sphere_steady_temperature(r, **KW))
        worst = max(worst, abs(float(_theta(r, 6e-3)) / ref - 1.0))
    print(f"\nMAT11 N4 worst |theta/theta_MAT9 - 1| = {worst:.3e} (hand 2.89e-5)")
    assert worst <= 1e-4


def test_n5_strong_perfusion_limits():
    """(N5) mR = 40: interior wθ/q → 1; wθ(R)/q = (mR − 1)/(mR + Bi − 1)."""
    w = K * (40.0 / R) ** 2
    interior = max(abs(w * float(_theta(r, w)) / Q - 1.0) for r in (0.0, 0.5 * R))
    surface = w * float(_theta(R, w)) / Q
    limit = (40.0 - 1.0) / (40.0 + BI - 1.0)
    print(f"\nMAT11 N5 interior |w theta/q - 1| = {interior:.3e}; "
          f"surface w theta(R)/q = {surface:.12f} vs {limit:.12f} (0.959016)")
    assert interior <= 1e-8
    assert abs(surface - limit) <= 1e-12
    # No overflow far past sinh's range.
    assert np.isfinite(_theta(np.array([0.0, 0.5 * R, R]), K * (2000.0 / R) ** 2)).all()
    with pytest.raises(ValueError):
        robin_sphere_perfusion_temperature(0.0, w=0.0, **KW)
    with pytest.raises(ValueError):
        robin_sphere_perfusion_coefficients(w=-1.0, **KW)


# ------------------------------------------------------------------- FEM ---


def _errors(msh, T_h, w):
    """(interior rel. L2 on θ, whole-ball rel. L2 on θ, V_h, N_global)."""
    comm = msh.comm
    m, _D, A, _mean = robin_sphere_perfusion_coefficients(w=w, **KW)
    x = ufl.SpatialCoordinate(msh)
    mr = m * ufl.sqrt(ufl.inner(x, x))
    theta = Q / w * (1.0 - A * ufl.sinh(mr) / mr)
    theta_h = T_h - T_INF
    dx = ufl.dx(domain=msh, metadata=QDEG)
    chi = _interior_indicator(msh)
    # Complex-safe (`MAT-10`'s `_interior_rel` idiom): squares as inner(d, d)
    # = |d|², and np.real before every float — float() raises on a Python
    # complex (the allreduced sum in the complex build).
    d = theta_h - theta
    e_in = np.real(_allsum(comm, chi * ufl.inner(d, d) * dx))
    n_in = np.real(_allsum(comm, chi * ufl.inner(theta, theta) * dx))
    e_all = np.real(_allsum(comm, ufl.inner(d, d) * dx))
    n_all = np.real(_allsum(comm, ufl.inner(theta, theta) * dx))
    vol = np.real(_allsum(comm, fem.Constant(msh, default_scalar_type(1.0)) * dx))
    n_glob = msh.topology.index_map(msh.topology.dim).size_global
    return float(np.sqrt(e_in / n_in)), float(np.sqrt(e_all / n_all)), float(vol), int(n_glob)


def _balance(sol):
    return abs(sol.convected_power + sol.perfused_power - sol.source_power) / abs(sol.source_power)


def _centre_surface(msh, T_h):
    """(θ_h at the vertex nearest the centre, mean θ_h over the Robin surface)."""
    comm = msh.comm
    V = T_h.function_space
    n_owned = V.dofmap.index_map.size_local
    xyz = V.tabulate_dof_coordinates()[:n_owned]
    r = np.linalg.norm(xyz, axis=1)
    local = (
        (float(r.min()), float(np.real(T_h.x.array[:n_owned][int(np.argmin(r))] - T_INF)))
        if n_owned else (np.inf, 0.0)
    )
    centre = min(comm.allgather(local))[1]
    ds = ufl.Measure("ds", domain=msh)
    area = np.real(_allsum(comm, fem.Constant(msh, default_scalar_type(1.0)) * ds))
    surf = np.real(_allsum(comm, (T_h - T_INF) * ds)) / area
    return float(centre), float(surf)


def _max_abs_imag(T_h):
    """Global max|Im θ_h| over owned dofs (real floats, so MAX is safe)."""
    V = T_h.function_space
    n_owned = V.dofmap.index_map.size_local
    loc = float(np.max(np.abs(np.imag(T_h.x.array[:n_owned])))) if n_owned else 0.0
    return V.mesh.comm.allreduce(loc, op=MPI.MAX)


def test_p1_ladder_interior_rate_and_balance():
    comm = MPI.COMM_WORLD
    rows = []
    for res in LADDER:
        t0 = time.perf_counter()
        msh, _, ft = _ball(res)
        sol = solve_steady_conduction(
            msh, ft, k=K, q=Q, h=H, t_inf=T_INF, robin_tag=ROBIN_TAG, degree=1,
            perfusion_w_per_m3_k=W,
        )
        e_in, e_all, vol, n = _errors(msh, sol.temperature, W)
        bal = _balance(sol)
        h_eff = (vol / n) ** (1.0 / 3.0)
        centre, surf = _centre_surface(msh, sol.temperature)
        rows.append((res, n, h_eff, e_in, e_all, bal))
        if res == 0.01:
            _say(comm, f"MAT11 max|Im theta_h| at 0.01 m = {_max_abs_imag(sol.temperature):.3e} "
                       f"(predicted 0.0; ScalarType={np.dtype(default_scalar_type).name})")
        _say(
            comm,
            f"MAT11 P1 res={res:.3f} cells={n} h_eff={h_eff:.5e} "
            f"interior_rel_L2={e_in:.4e} whole_rel_L2={e_all:.4e} "
            f"src={sol.source_power:.6f}W conv={sol.convected_power:.6f}W "
            f"perf={sol.perfused_power:.6f}W balance={bal:.3e} "
            f"theta_h centre={centre:.4f} surface-mean={surf:.4f} (cf 5.5394 / 3.2668) "
            f"Vh/V={vol / V_BALL:.6f} ({time.perf_counter() - t0:.1f}s)",
        )
        # (iii) on every rung.
        assert bal <= BALANCE_BAND, f"(iii) balance {bal:.3e} at res={res}"

    h = np.array([r[2] for r in rows])
    e = np.array([r[4] for r in rows])
    for i in range(len(rows) - 1):
        _say(comm, f"MAT11 pairwise rate {LADDER[i]}->{LADDER[i + 1]}: "
                   f"{np.log(e[i] / e[i + 1]) / np.log(h[i] / h[i + 1]):.3f}")
    rate = float(np.polyfit(np.log(h), np.log(e), 1)[0])
    interior_001 = rows[LADDER.index(0.01)][3]
    _say(comm, f"MAT11 least-squares P1 L2 rate (whole ball) = {rate:.3f}  (predicted ≈ 2.0; floor {RATE_FLOOR})")
    _say(comm, f"MAT11 (i) interior rel L2 at 0.01 m = {interior_001:.4e}  "
               f"(band {INTERIOR_BAND:.0e}; predicted ≈ 2.1e-3; MAT-9 read 3.21e-3)")
    assert rate >= RATE_FLOOR, f"(ii) P1 rate {rate:.3f}"
    assert interior_001 <= INTERIOR_BAND, f"(i) interior {interior_001:.4e}"


def test_strong_perfusion_printed_only():
    comm = MPI.COMM_WORLD
    msh, _, ft = _ball(0.01)
    sol = solve_steady_conduction(
        msh, ft, k=K, q=Q, h=H, t_inf=T_INF, robin_tag=ROBIN_TAG, degree=1,
        perfusion_w_per_m3_k=W_STRONG,
    )
    e_in, e_all, _vol, _n = _errors(msh, sol.temperature, W_STRONG)
    _say(comm, f"MAT11 printed-only w={W_STRONG:.0f}: interior_rel_L2={e_in:.4e} whole_rel_L2={e_all:.4e} "
               f"perfused fraction={sol.perfused_power / sol.source_power:.4f} (predicted 0.905) "
               f"balance={_balance(sol):.3e}")
    with pytest.raises(ValueError):
        solve_steady_conduction(
            msh, ft, k=K, q=Q, h=H, t_inf=T_INF, robin_tag=ROBIN_TAG, perfusion_w_per_m3_k=-W,
        )


def test_negative_controls():
    comm = MPI.COMM_WORLD
    msh, _, ft = _ball(0.01)

    # (a) perfusion dropped: MAT-9's field against the perfused closed form.
    sol_a = solve_steady_conduction(
        msh, ft, k=K, q=Q, h=H, t_inf=T_INF, robin_tag=ROBIN_TAG, degree=1
    )
    e_in_a, _, _, _ = _errors(msh, sol_a.temperature, W)
    _say(comm, f"MAT11 control (a) w = 0: interior_rel_L2={e_in_a:.4e} "
               f"(asserted > {INTERIOR_BAND:.0e}; predicted ≈ 1.11); perfused_power={sol_a.perfused_power!r}")
    assert e_in_a > INTERIOR_BAND
    assert sol_a.perfused_power == 0.0

    # (b) sign flipped: −wθ is MAT-5's growth form at mR = 2 > λ₁ (indefinite).
    sol_b = _solve(
        msh, ft, k=K, q=Q, h=H, t_inf=T_INF, robin_tag=ROBIN_TAG, degree=1,
        perfusion_w_per_m3_k=-W, petsc_prefix="mat11_ctrl_b_",
    )
    e_in_b, _, _, _ = _errors(msh, sol_b.temperature, W)
    centre_b, _surf_b = _centre_surface(msh, sol_b.temperature)
    _say(comm, f"MAT11 control (b) -w: interior_rel_L2={e_in_b:.4e} theta_h centre={centre_b:.2f} K "
               f"(asserted > {INTERIOR_BAND:.0e}; predicted ≈ 22, theta(0) ≈ -131 K, mesh-sensitive)")
    assert e_in_b > INTERIOR_BAND
