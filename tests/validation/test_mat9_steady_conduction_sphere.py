"""`MAT-9` step 1 — steady conduction in a Robin-cooled ball (§9 item 66).

Fixture: R 0.1 m, k 0.6 W/m/K, Q 2000 W/m³ (2 W/kg × 1000 kg/m³),
h 10 W/m²/K (Bi = hR/k = 1.67), T∞ 20 °C.  Closed form
(`utils.analytical.robin_sphere_steady_temperature`):
θ = T − T∞ = QR/(3h) + Q(R² − r²)/(6k): 12.222 K at the centre, 6.667 K at
the surface; ∫Q dV = 8.378 W.  **Errors are taken on θ, never on T** (T∞
would dilute every relative figure).

Anchors (asserted): (i) interior ‖θ_h − θ‖/‖θ‖ ≤ 1 % on cells whose
midpoint has |x|² < 0.81R², finest P1 rung; (ii) least-squares P1 L2 rate
(whole ball) vs h_eff = (V_h/N_global)^⅓ ≥ 1.7 on the 0.02/0.01/0.005 m
ladder, one-sided (`MAG-18` precedent); (iii) |∮hθ_h dS − ∫Q dV|/∫Q dV ≤ 1e-9
on every rung and at P2.  Negative controls: (a) Dirichlet θ = 0 on r = R
must read (i) > 1 % (*asserted*; *predicted* ≈ 69.5 %); (b) a flipped Robin
term must read (iii) > 1e-3 (*asserted*; *predicted* 2.000).  Printed only:
a flipped Q passes (iii) and reads ≈ 200 % on (i).
"""

from __future__ import annotations

import time

import numpy as np
import pytest
import ufl
from dolfinx import default_scalar_type, fem
from dolfinx import mesh as dmesh
from mpi4py import MPI

from fem_em_solver.io.mesh import MeshGenerator
from fem_em_solver.thermal.steady import _solve, solve_steady_conduction
from fem_em_solver.utils.analytical import robin_sphere_steady_temperature

R = 0.1
K = 0.6
Q = 2000.0
H = 10.0
T_INF = 20.0
ROBIN_TAG = 1
LADDER = (0.02, 0.01, 0.005)
QDEG = {"quadrature_degree": 4}

_MESHES: dict = {}


def _ball(res):
    if res not in _MESHES:
        _MESHES[res] = MeshGenerator.ball_domain(radius=R, resolution=res)
    return _MESHES[res]


def _allsum(comm, form):
    return comm.allreduce(fem.assemble_scalar(fem.form(form)), op=MPI.SUM)


def _interior_indicator(msh):
    """DG0 indicator of cells with midpoint |x|² < 0.81 R² (numpy mask)."""
    tdim = msh.topology.dim
    imap = msh.topology.index_map(tdim)
    n = imap.size_local + imap.num_ghosts
    mids = dmesh.compute_midpoints(msh, tdim, np.arange(n, dtype=np.int32))
    W = fem.functionspace(msh, ("DG", 0))
    chi = fem.Function(W)
    chi.x.array[:] = (np.sum(mids**2, axis=1) < 0.81 * R**2).astype(float)
    chi.x.scatter_forward()
    return chi


def _errors(msh, T_h):
    """(interior rel. L2 on θ, whole-ball rel. L2 on θ, V_h, N_global)."""
    comm = msh.comm
    x = ufl.SpatialCoordinate(msh)
    r2 = ufl.inner(x, x)
    theta = Q * R / (3.0 * H) + Q * (R**2 - r2) / (6.0 * K)
    theta_h = T_h - T_INF
    dx = ufl.dx(domain=msh, metadata=QDEG)
    chi = _interior_indicator(msh)
    e_in = _allsum(comm, chi * (theta_h - theta) ** 2 * dx)
    n_in = _allsum(comm, chi * theta**2 * dx)
    e_all = _allsum(comm, (theta_h - theta) ** 2 * dx)
    n_all = _allsum(comm, theta**2 * dx)
    vol = _allsum(comm, fem.Constant(msh, default_scalar_type(1.0)) * dx)
    n_glob = msh.topology.index_map(msh.topology.dim).size_global
    return (
        float(np.sqrt(e_in / n_in)),
        float(np.sqrt(e_all / n_all)),
        float(vol),
        int(n_glob),
    )


def _balance(res):
    return abs(res.convected_power - res.source_power) / abs(res.source_power)


def _say(comm, msg):
    if comm.rank == 0:
        print(msg, flush=True)


def test_closed_form_values():
    """The fixture's quoted numbers follow from the closed form."""
    tc = robin_sphere_steady_temperature(0.0, radius=R, k=K, q=Q, h=H, t_inf=T_INF)
    ts = robin_sphere_steady_temperature(R, radius=R, k=K, q=Q, h=H, t_inf=T_INF)
    assert tc - T_INF == pytest.approx(12.2222222, rel=1e-6)
    assert ts - T_INF == pytest.approx(6.6666667, rel=1e-6)
    assert 4.0 / 3.0 * np.pi * R**3 * Q == pytest.approx(8.37758, rel=1e-5)
    # Surface flux -k dT/dr = QR/3 matches h(T_s - T∞).
    r = np.array([R - 1e-7, R])
    dTdr = (
        np.diff(robin_sphere_steady_temperature(r, radius=R, k=K, q=Q, h=H, t_inf=T_INF))[0] / 1e-7
    )
    assert -K * dTdr == pytest.approx(H * (ts - T_INF), rel=1e-5)


def test_p1_ladder_interior_rate_and_balance():
    comm = MPI.COMM_WORLD
    V_exact = 4.0 / 3.0 * np.pi * R**3
    rows = []
    for res in LADDER:
        t0 = time.perf_counter()
        msh, _, ft = _ball(res)
        sol = solve_steady_conduction(
            msh, ft, k=K, q=Q, h=H, t_inf=T_INF, robin_tag=ROBIN_TAG, degree=1
        )
        e_in, e_all, vol, n = _errors(msh, sol.temperature)
        bal = _balance(sol)
        h_eff = (vol / n) ** (1.0 / 3.0)
        rows.append((res, n, h_eff, e_in, e_all, bal, vol))
        _say(
            comm,
            f"MAT9 P1 res={res:.3f} cells={n} h_eff={h_eff:.5e} "
            f"interior_rel_L2={e_in:.4e} whole_rel_L2={e_all:.4e} "
            f"src={sol.source_power:.6f}W conv={sol.convected_power:.6f}W "
            f"balance={bal:.3e} Vh/V={vol / V_exact:.6f} "
            f"({time.perf_counter() - t0:.1f}s)",
        )
        # (iii) on every rung.
        assert bal <= 1e-9, f"(iii) balance {bal:.3e} > 1e-9 at res={res}"
        assert sol.source_power == pytest.approx(Q * vol, rel=1e-12)

    h = np.array([r[2] for r in rows])
    e = np.array([r[4] for r in rows])
    for i in range(len(rows) - 1):
        _say(
            comm,
            f"MAT9 pairwise rate {LADDER[i]}->{LADDER[i + 1]}: "
            f"{np.log(e[i] / e[i + 1]) / np.log(h[i] / h[i + 1]):.3f}",
        )
    rate = float(np.polyfit(np.log(h), np.log(e), 1)[0])
    _say(comm, f"MAT9 least-squares P1 L2 rate (whole ball) = {rate:.3f}  (predicted ≈ 2)")
    _say(comm, f"MAT9 finest interior rel L2 = {rows[-1][3]:.4e}  (predicted ≈ 1e-3)")
    # (ii) one-sided rate.
    assert rate >= 1.7, f"(ii) P1 rate {rate:.3f} < 1.7"
    # (i) interior accuracy, finest P1 rung.
    assert rows[-1][3] <= 1e-2, f"(i) interior {rows[-1][3]:.4e} > 1 %"


def test_p2_balance_and_faceting():
    comm = MPI.COMM_WORLD
    res = 0.01
    msh, _, ft = _ball(res)
    sol = solve_steady_conduction(
        msh, ft, k=K, q=Q, h=H, t_inf=T_INF, robin_tag=ROBIN_TAG, degree=2
    )
    e_in, e_all, vol, n = _errors(msh, sol.temperature)
    bal = _balance(sol)
    _say(
        comm,
        f"MAT9 P2 res={res} cells={n} interior_rel_L2={e_in:.4e} "
        f"whole_rel_L2={e_all:.4e} balance={bal:.3e} "
        f"Vh/V={vol / (4.0 / 3.0 * np.pi * R**3):.6f} (predicted: faceting-limited)",
    )
    assert bal <= 1e-9, f"(iii) P2 balance {bal:.3e} > 1e-9"


def test_negative_controls():
    comm = MPI.COMM_WORLD
    msh, _, ft = _ball(0.01)

    # (a) Dirichlet θ = 0 on r = R: drops the QR/(3h) offset.
    sol_a = _solve(
        msh,
        ft,
        k=K,
        q=Q,
        h=H,
        t_inf=T_INF,
        robin_tag=ROBIN_TAG,
        degree=1,
        dirichlet_value=T_INF,
        petsc_prefix="mat9_ctrl_a_",
    )
    e_in_a, _, _, _ = _errors(msh, sol_a.temperature)
    _say(
        comm,
        f"MAT9 control (a) Dirichlet: interior_rel_L2={e_in_a:.4e} "
        "(asserted > 1e-2; predicted ≈ 0.695)",
    )
    assert e_in_a > 1e-2, f"control (a) reads {e_in_a:.4e} <= 1 %"

    # (b) Robin term with flipped sign.
    sol_b = _solve(
        msh,
        ft,
        k=K,
        q=Q,
        h=H,
        t_inf=T_INF,
        robin_tag=ROBIN_TAG,
        degree=1,
        robin_sign=-1.0,
        petsc_prefix="mat9_ctrl_b_",
    )
    bal_b = _balance(sol_b)
    _say(
        comm,
        f"MAT9 control (b) flipped Robin: balance={bal_b:.6f} "
        "(asserted > 1e-3; predicted 2.000)",
    )
    assert bal_b > 1e-3, f"control (b) balance {bal_b:.3e} <= 1e-3"

    # Printed only: flipped Q passes (iii) and reads ≈ 200 % on (i).
    sol_c = _solve(
        msh,
        ft,
        k=K,
        q=-Q,
        h=H,
        t_inf=T_INF,
        robin_tag=ROBIN_TAG,
        degree=1,
        petsc_prefix="mat9_ctrl_c_",
    )
    e_in_c, _, _, _ = _errors(msh, sol_c.temperature)
    _say(
        comm,
        f"MAT9 printed-only flipped Q: balance={_balance(sol_c):.3e} "
        f"interior_rel_L2={e_in_c:.4e} (predicted passes (iii), ≈ 2.0 on (i))",
    )
