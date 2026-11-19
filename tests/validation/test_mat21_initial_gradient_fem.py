"""`MAT-21` step 19 FEM (§9 item 669, §10 P47b) — the source-free Robin gel box
from the non-uniform initial temperature T(0) = T∞ + g·z, P1 CN on step 1's
finest rung, against item 659's closed form
(`utils.analytical.robin_box_initial_gradient_temperature`).

No EM solve; real build.  Mesh: step 1's ``create_box`` at m = ``RUNGS[-1]`` = 6,
``BASE`` = (13, 8, 2) (269 568 tetrahedra).  CN N = ``N_MAIN`` = 40 to
t* = 0.2 t_C with two BE half-steps (step 1's (T) schedule), q = 0, the initial
state interpolated through the additive keyword ``t_init_expr`` of
``solve_transient_conduction`` (a P1 interpolant of g·z is exact).  With no
source the per-step balance is normalised by the initial stored energy per
step, ∫ρc|θ₀| dV / Δtₙ (the new branch of that function).

Asserted: (P0) the series' constants = the FEM's; (I) interior rel L2 at t*
against the closed form (P2 interpolant, n_terms ``SERIES_TERMS``) ≤
``INTERIOR_BAND``; (E) the source-free balance ≤ ``E_BAND``; (Sym) at the five
F2182 probe points on z = 0, |θ(t*)| ≤ 1e-10·g C (a symmetry identity).
Control (Dz, asserted as a miss of ``INTERIOR_BAND``): the same FEM against
item 659's Dirichlet-z closed form (odd roots γ = mπ).  Printed: P5's drift
θ(0, 0, C, t*) beside the series and P37's **steady sourced** value
3.127848 K (``SERIES_DRAFT[5]``), with their ratio.

(Sym) re-posed as (Sym′) by the 2026-10-27 18:00 daily review ruling (§9 item
676, ``docs/planning/chunks/MAT-21.md``): item 669's z = 0 bar went red at
2.2e-3 g C because ``create_box``'s tetrahedra are not mirrored under z → −z;
they do map onto themselves under point reflection through the centre, and
θ₀ = g·z is odd under it, so the exact discrete identity is
θ_h(−p) = −θ_h(p).  Asserted: |θ(P1) + θ(P2)|, |θ(P3) + θ(P4)| and |θ(P0)|
each ≤ ``SYM_BAR``·g C (the bar unchanged); the four wall values print as the
O(h²) mid-plane record.  Not claimed: a z-mirror symmetry.
"""

from __future__ import annotations

import time

import numpy as np
import pytest
import ufl
from dolfinx import default_scalar_type, fem
from dolfinx import mesh as dmesh
from mpi4py import MPI

from fem_em_solver.post.evaluation import evaluate_vector_field_parallel
from fem_em_solver.thermal import solve_transient_conduction
from fem_em_solver.utils.analytical import robin_box_initial_gradient_temperature
from tests.validation.test_mat10_sar_to_heat_sphere import ROBIN_TAG, _allsum, _robin_tags, _say
from tests.validation.test_mat12_transient_sar_heating_sphere import E_BAND, RHO_C
from tests.validation.test_mat13_transient_perfusion_sphere import INTERIOR_BAND
from tests.validation.test_mat15_coil_driven_heating_phantom import H_GEL, K_GEL, N_MAIN, T_INF
from tests.validation.test_mat15_finite_cylinder_fem import _max_imag, _rel
from tests.validation.test_mat21_box_closed_form import (
    A_BOX,
    B_BOX,
    C_BOX,
    HALF,
)
from tests.validation.test_mat21_box_closed_form import RHO_C as RHO_C_0
from tests.validation.test_mat21_box_closed_form import SERIES_TERMS, T_C, H, K
from tests.validation.test_mat21_box_fem import BASE, RUNGS, T_STAR, _exact, _interior_chi
from tests.validation.test_mat21_box_probe_points_fem import PROBE_LABELS, PROBE_PTS, SERIES_DRAFT
from tests.validation.test_mat21_initial_gradient_closed_form import G_GRAD, GC

SYM_BAR = 1e-10  # × g C, the brief's symmetry-identity bar
DZ_ROOTS = np.arange(1, 41) * np.pi  # item 659's (Dz) control roots
P37_STEADY_P5 = SERIES_DRAFT[5]  # 3.127848 K — steady, sourced (P37), not a transient rise


class _Clock19:
    def __init__(self):
        self.t = time.perf_counter()

    def lap(self, phase):
        now = time.perf_counter()
        _say(f"MAT21-19F clock {phase}: {now - self.t:.1f} s")
        self.t = now


def _grad_ex(t, roots_z=None):
    def f(x, y, z):
        return robin_box_initial_gradient_temperature(
            x,
            y,
            z,
            t,
            half_widths=HALF,
            k=K,
            h=H,
            rho_c=RHO_C_0,
            g=G_GRAD,
            t_inf=0.0,
            n_terms=SERIES_TERMS,
            roots_z=roots_z,
        )

    return f


@pytest.fixture(scope="module")
def gradient_run():
    comm = MPI.COMM_WORLD
    t0 = time.perf_counter()
    clk = _Clock19()
    m = RUNGS[-1]
    ns = [b * m for b in BASE]
    msh = dmesh.create_box(
        comm,
        [np.array([-A_BOX, -B_BOX, -C_BOX]), np.array([A_BOX, B_BOX, C_BOX])],
        ns,
        cell_type=dmesh.CellType.tetrahedron,
    )
    ft = _robin_tags(msh)
    n = msh.topology.index_map(3).size_global
    chi = _interior_chi(msh)
    _say(
        f"\nMAT21-19F rung m = {m} {tuple(ns)}: cells {n}; g = {G_GRAD} K/m, g C = {GC} K; "
        f"t_C = {T_C:.3f} s, t* = {T_STAR:.3f} s, CN N = {N_MAIN} + 2 BE half-steps; ranks {comm.size}"
    )
    clk.lap("build")

    res = solve_transient_conduction(
        msh,
        ft,
        k=K_GEL,
        rho_c=RHO_C,
        q=0.0,
        h=H_GEL,
        t_inf=T_INF,
        t_init=T_INF,
        robin_tag=ROBIN_TAG,
        dt=T_STAR / N_MAIN,
        n_steps=N_MAIN,
        theta=0.5,
        startup_half_steps=2,
        degree=1,
        snapshot_steps=(0,),
        t_init_expr=lambda x: T_INF + G_GRAD * x[2],
    )
    clk.lap(f"CN N={N_MAIN}")
    out = {
        "E": res.max_balance_residual,
        "src": res.source_power,
        "imag": _max_imag(res.temperature),
    }

    # The initial state: exact P1 interpolant of g·z, and ∫θ₀ = 0 by symmetry (printed).
    th0 = res.snapshots[0] - T_INF
    xs = ufl.SpatialCoordinate(msh)
    e_ic = _rel(msh, th0, G_GRAD * xs[2])
    int0 = float(np.real(_allsum(comm, th0 * ufl.dx(domain=msh))))
    one = fem.Constant(msh, default_scalar_type(1.0))
    v_h = float(np.real(_allsum(comm, one * ufl.dx(domain=msh))))
    _say(
        f"MAT21-19F initial state: rel L2 |theta_0 - g z| = {e_ic:.3e} (P1 exact, printed); "
        f"int theta_0 dV/(g C V) = {int0 / (GC * v_h):.3e} (printed, symmetry)"
    )

    th = res.temperature - T_INF
    ex, npts = _exact(msh, _grad_ex(T_STAR))
    clk.lap(f"series interpolation at t* ({npts} P2 points)")
    out["I"] = _rel(msh, th, ex, chi)
    whole = _rel(msh, th, ex)
    _say(
        f"MAT21-19F (I) CN at t*: interior rel L2 = {out['I']:.6e} (bound {INTERIOR_BAND:.0e}); "
        f"whole {whole:.6e} (printed)"
    )
    nb = len(res.balance_residuals)
    _say(
        f"MAT21-19F (E) source power int Q = {res.source_power!r} W; balance (new normalisation "
        f"|stor + conv| dt_n / int rho c |theta_0|) max {res.max_balance_residual:.3e} over {nb} "
        f"(half-)steps (bound {E_BAND:.0e}); first {res.balance_residuals[0]:.3e}, "
        f"last {res.balance_residuals[-1]:.3e}"
    )

    ex_dz, _ = _exact(msh, _grad_ex(T_STAR, roots_z=DZ_ROOTS))
    out["Dz"] = _rel(msh, th, ex_dz, chi)
    _say(
        f"MAT21-19F control (Dz) Dirichlet-z closed form at t*: interior miss = {out['Dz']:.4e} "
        f"(asserted > {INTERIOR_BAND:.0e})"
    )
    clk.lap("control (Dz)")

    pts = np.array(PROBE_PTS, dtype=float)
    vals, ok = evaluate_vector_field_parallel(res.temperature, pts)
    fem_p = np.array(
        [
            float(np.real(np.ravel(vals[i])[0])) - T_INF if bool(ok[i]) else np.nan
            for i in range(len(PROBE_PTS))
        ]
    )
    ser_p = np.real(_grad_ex(T_STAR)(pts[:, 0], pts[:, 1], pts[:, 2]))
    out["probes"] = fem_p
    out["ok"] = [bool(o) for o in ok]
    for i, lab in enumerate(PROBE_LABELS):
        _say(
            f"MAT21-19F probe {lab} {PROBE_PTS[i]}: theta_h(t*) = {fem_p[i]:+.6e} K, series "
            f"{ser_p[i]:+.6e} K, |theta_h|/(g C) = {abs(fem_p[i]) / GC:.3e}"
        )
    p5, s5 = fem_p[5], ser_p[5]
    _say(
        f"MAT21-19F P5 drift theta(0, 0, C, t*) = {p5:.6f} K (series {s5:.6f} K, "
        f"rel {p5 / s5 - 1.0:+.3e}); beside it P37 = {P37_STEADY_P5:.6f} K (steady, sourced (P37), "
        f"not a transient rise); ratio drift / P37 = {p5 / P37_STEADY_P5:.6f}"
    )
    _say(f"MAT21-19F max|Im theta| = {out['imag']:.3e}; elapsed {time.perf_counter() - t0:.1f} s")
    return out


def test_precondition_constants():
    assert K == K_GEL
    assert H == H_GEL
    assert RHO_C_0 == RHO_C


def test_i_interior_against_closed_form(gradient_run):
    assert gradient_run["I"] <= INTERIOR_BAND, f"(I) {gradient_run['I']:.4e} > {INTERIOR_BAND}"


def test_e_source_free_balance(gradient_run):
    assert gradient_run["src"] == 0.0
    assert gradient_run["E"] <= E_BAND, f"(E) {gradient_run['E']:.3e} > {E_BAND:.0e}"


def test_sym_point_reflection_pairs(gradient_run):
    assert all(gradient_run["ok"]), gradient_run["ok"]
    p = gradient_run["probes"]
    s12, s34, s0 = abs(p[1] + p[2]), abs(p[3] + p[4]), abs(p[0])
    for i in range(1, 5):
        _say(
            f"MAT21-19F O(h^2) mid-plane record on a point-symmetric, not z-mirrored, mesh: "
            f"{PROBE_LABELS[i]} theta_h(t*) = {p[i]:+.6e} K, |theta_h|/(g C) = {abs(p[i]) / GC:.3e}"
        )
    _say(
        f"MAT21-19F (Sym') |P1 + P2|/(g C) = {s12 / GC:.3e}, |P3 + P4|/(g C) = {s34 / GC:.3e}, "
        f"|P0|/(g C) = {s0 / GC:.3e} (bar {SYM_BAR:.0e})"
    )
    assert s12 <= SYM_BAR * GC, ("P1 + P2", p[1], p[2])
    assert s34 <= SYM_BAR * GC, ("P3 + P4", p[3], p[4])
    assert s0 <= SYM_BAR * GC, ("P0", p[0])


def test_control_dz_misses(gradient_run):
    assert gradient_run["Dz"] > INTERIOR_BAND, f"(Dz) {gradient_run['Dz']:.4e} inside band"
