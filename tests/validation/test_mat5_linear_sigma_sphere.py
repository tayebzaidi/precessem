"""`MAT-5` step 1 — fixed-|E| Joule-heated Robin ball, σ(T) linear (§9 item 85).

`MAT-9` step 1's fixture (imported by name): R 0.1, k 0.6, h 10 (Bi 5/3),
T∞ = T₀ = 20 °C, q₀ = 2000 W/m³; α = 0.02 /K (order of saline's; not a
material claim).  Holding |E| fixed, Q = q₀·σ(T)/σ₀ = q₀[1 + αθ] with the
factor taken from ``LinearTemperatureConductivity.factor_ufl``.  Closed form
``robin_sphere_linear_source_temperature``; m² = 66.667 m⁻², mR = 0.81650.
**Errors on θ, never on T.**

FEM anchors (asserted; P1, the 0.01 m rung, Picard from T₀ = T∞):
(i) interior (|x|² < 0.81R²) relative L2 on θ ≤ 1 % (predicted 0.3–0.4 %;
`MAT-9` read 0.321 % at this rung); (ii) increment ratio δ₆/δ₅,
δₙ = ‖T_{n+1} − T_n‖_L2, within 1 % of κ = (mR/λ₁)²; (iii) final-iterate
balance |∮hθ_N dS − ∫Q(T_{N−1}) dV|/∫Q ≤ 1e-9.  Printed: iteration count
(≈ 17), implied λ₁,h = mR/√ratio.
Negative controls (predicted; asserted only as "misses (i)'s 1 %"):
(a) α dropped — miss ≈ 18.7 %; (b) σ evaluated at T instead of T − T₀ —
θ_ctrl = 1.4θ exactly, miss 40.0 %.
"""

from __future__ import annotations

import time

import numpy as np
import ufl
from dolfinx import fem
from mpi4py import MPI

from fem_em_solver.materials import LinearTemperatureConductivity
from fem_em_solver.thermal.steady import solve_steady_conduction_picard
from fem_em_solver.utils.analytical import robin_sphere_linear_source_coefficients
from tests.validation.test_mat9_steady_conduction_sphere import (  # MAT-9 fixture
    H,
    K,
    QDEG,
    R,
    ROBIN_TAG,
    T_INF,
    Q,
    _allsum,
    _ball,
    _interior_indicator,
)

Q0 = Q
ALPHA = 0.02
T0 = T_INF
RES = 0.01
SIGMA0 = 0.57  # cancels in Q = q0 σ(T)/σ0; any positive value

M, D, A, LAM1 = robin_sphere_linear_source_coefficients(
    radius=R, k=K, q0=Q0, alpha=ALPHA, h=H
)
KAPPA = (M * R / LAM1) ** 2


def _say(comm, msg):
    if comm.rank == 0:
        print(msg, flush=True)


def _interior_rel(msh, T_h, chi):
    """Interior relative L2 on θ against the closed form (allreduced)."""
    comm = msh.comm
    x = ufl.SpatialCoordinate(msh)
    r = ufl.sqrt(ufl.inner(x, x))
    theta = (A * ufl.sin(M * r) / r) - 1.0 / ALPHA  # quadrature points never at r = 0
    dx = ufl.dx(domain=msh, metadata=QDEG)
    e = _allsum(comm, chi * (T_h - T_INF - theta) ** 2 * dx)
    n = _allsum(comm, chi * theta**2 * dx)
    return float(np.sqrt(np.real(e) / np.real(n)))


def _picard(msh, ft, source):
    return solve_steady_conduction_picard(
        msh, ft, k=K, source=source, h=H, t_inf=T_INF, robin_tag=ROBIN_TAG,
        degree=1, rtol=1e-12, max_iter=60,
    )


def test_picard_fixed_point_rate_and_balance():
    comm = MPI.COMM_WORLD
    t0 = time.perf_counter()
    msh, _, ft = _ball(RES)
    chi = _interior_indicator(msh)
    n_cells = msh.topology.index_map(msh.topology.dim).size_global
    cond = LinearTemperatureConductivity(sigma_ref=SIGMA0, alpha_per_k=ALPHA, t_ref=T0)

    sol = _picard(msh, ft, lambda T: Q0 * cond.sigma_ufl(T) / SIGMA0)
    inc = np.array(sol.increments)
    ratios = inc[1:] / inc[:-1]
    for i, d in enumerate(inc):
        rr = f" ratio δ{i}/δ{i - 1}={ratios[i - 1]:.6f}" if i else ""
        _say(comm, f"MAT5 Picard n={i} δ={d:.6e}{rr}")
    e_in = _interior_rel(msh, sol.temperature, chi)
    bal = abs(sol.convected_power - sol.source_power) / abs(sol.source_power)
    ratio6 = float(inc[6] / inc[5])
    lam1_h = M * R / np.sqrt(ratio6)
    _say(
        comm,
        f"MAT5 P1 res={RES} cells={n_cells} iterations={sol.iterations} "
        f"converged={sol.converged} interior_rel_L2(θ)={e_in:.4e} "
        f"(pred. 0.3–0.4 %) src={sol.source_power:.6f}W "
        f"conv={sol.convected_power:.6f}W (closed-form P 10.1998 W) "
        f"balance={bal:.3e} ratio6={ratio6:.6f} κ={KAPPA:.6f} "
        f"ratio6/κ−1={ratio6 / KAPPA - 1:+.4e} implied λ1,h={lam1_h:.6f} "
        f"(λ1={LAM1:.6f}) ({time.perf_counter() - t0:.1f}s)",
    )
    assert sol.converged, f"Picard did not converge in {sol.iterations}"
    # (i)
    assert e_in <= 1e-2, f"(i) interior {e_in:.4e} > 1 %"
    # (ii)
    assert abs(ratio6 / KAPPA - 1.0) <= 1e-2, (
        f"(ii) ratio6 {ratio6:.6f} vs κ {KAPPA:.6f}"
    )
    # (iii)
    assert bal <= 1e-9, f"(iii) balance {bal:.3e} > 1e-9"


def test_negative_controls():
    comm = MPI.COMM_WORLD
    msh, _, ft = _ball(RES)
    chi = _interior_indicator(msh)

    # (a) α dropped: Q = q0 (the MAT-9 source).
    q0_c = fem.Constant(msh, Q0)
    sol_a = _picard(msh, ft, lambda T: q0_c)
    e_a = _interior_rel(msh, sol_a.temperature, chi)
    _say(comm, f"MAT5 control (a) α dropped: interior_rel_L2={e_a:.4e} "
         f"iterations={sol_a.iterations} (asserted > 1e-2; predicted ≈ 0.187)")
    assert e_a > 1e-2, f"control (a) reads {e_a:.4e} <= 1 %"

    # (b) σ evaluated at T instead of T − T0: factor 1 + αT = 1.4 + αθ.
    sol_b = _picard(msh, ft, lambda T: Q0 * (1.0 + ALPHA * T))
    e_b = _interior_rel(msh, sol_b.temperature, chi)
    _say(comm, f"MAT5 control (b) σ at T: interior_rel_L2={e_b:.4e} "
         f"iterations={sol_b.iterations} (asserted > 1e-2; predicted 0.400)")
    assert e_b > 1e-2, f"control (b) reads {e_b:.4e} <= 1 %"
