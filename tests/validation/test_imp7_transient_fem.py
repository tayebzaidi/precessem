"""`IMP-7` step 5 (§9 item 668, §10 I46): Crank–Nicolson transient on I43's
curved P2 adiabatic-hole fixture against I45's closed form.

Brief: docs/planning/chunks/IMP-7.md, "§9 item 668 — 2026-10-27".  Fixture:
the gel shell a < r < b (tag-2 sub-mesh of the geometry-order-2 two-layer
ball, `test_imp7_curved_p2_hole.py`), natural (zero-flux) condition on r = a,
Robin h on r = b, the β = 1 dipole SAR as a DG3 source, P2, T(0) = T∞.
Comparand: `pec_sphere_implant_transient_from_modes` (item 658; the same
series `pec_sphere_implant_transient_temperature` evaluates, modes computed
once) interpolated into P3, quadrature degree 6 (`DEG_ROUTE[2]`).

Gates, imported under `MAT-21` step 1's names (`test_mat21_box_fem.py:43–49`):
(T) CN interior rel L2 at t* ≤ INTERIOR_BAND; (S-slope) one BE step of DT_1,
|ρc∫θ¹/(DT_1∫Q) − 1| ≤ SLOPE_BAND (= `MAT-12`'s S2_BAND 1e-4, aliased: IMP-7's
own S2_BAND is 1e-2); (E) every run's balance ≤ E_BAND; (∞) BE 20 × 5 t_s
against the steady FEM solve ≤ INF_BAND (`MAT-21`'s (∞) form).  Negative
control (asserted): the CN transient misses the Dirichlet-at-a series by
> INTERIOR_BAND at t*.  Δt / n_steps / startup as `MAT-21` step 1 (CN N_MAIN
steps, two BE half-steps; DT_1 = 1e-5 of the time scale; BE 20 × 5 scale),
with the time scale t_s = (b − a)²/κ of this fixture (t* = 0.2 t_s).  Real
build.

Item 672 (review ruling 2026-10-27 10:30, IMP-7.md "Review ruling 2026-10-27
10:30"): item 668 found the DT_1 BE run's balance at T∞ = 20 K to be the
cancellation floor of a ≈ 8 µK step on a 20 K level (1.3e-9 vs 1.6e-14 in θ).
The `be_1` / slope run is therefore solved in θ (t_inf = t_init = 0) and its
(E) gated at the unchanged E_BAND; the T-variable `be_1` balance is a printed
record, beside a printed discriminator at 2·DT_1 (cancellation ⇒ ratio ≈ 0.5,
facet quadrature ⇒ ≈ 1).  CN and be_inf stay at T∞ = 20 K — they are the runs
that exercise the T∞ wiring.  `_gates` evaluates every gate before asserting.
"""

from __future__ import annotations

import time

import numpy as np
import ufl
from mpi4py import MPI

from fem_em_solver.thermal import solve_transient_conduction
from fem_em_solver.utils.analytical import (
    pec_sphere_implant_transient_from_modes,
    pec_sphere_implant_transient_modes,
)
from tests.validation import test_imp7_implant_heating_fem as m7
from tests.validation.test_imp7_curved_p2 import _ball_curved, _cmap_degree
from tests.validation.test_imp7_implant_heating_closed_form import BETA, FIX, H, K_GEL, T_INF
from tests.validation.test_imp7_transient_closed_form import (
    KAPPA,
    N_DEFAULT,
    SPAN,
    T_STAR,
)
from tests.validation.test_mat12_transient_sar_heating_sphere import (
    E_BAND,
    INF_BAND,
    N_MAIN,
    RHO_C,
)
from tests.validation.test_mat12_transient_sar_heating_sphere import S2_BAND as SLOPE_BAND
from tests.validation.test_mat13_transient_perfusion_sphere import INTERIOR_BAND

A, B, LADDER, ROBIN_TAG = m7.A, m7.B, m7.LADDER, m7.ROBIN_TAG
RT = m7.DEG_ROUTE[2]  # DG3 source, P3 comparand, quadrature 6
T_S = SPAN**2 / KAPPA  # this fixture's time scale; T_STAR = 0.2 T_S
DT_1 = 1.0e-5 * T_S  # MAT-21's DT_1 = 1e-5 t_C pattern
INTERP_CHUNK = 20000  # points per series call (bounds the (n_terms, N) mode matrix)


def _say(msg):
    if MPI.COMM_WORLD.rank == 0:
        print(msg, flush=True)


def _series_fn(t, inner):
    modes = pec_sphere_implant_transient_modes(beta=BETA, n_terms=N_DEFAULT, inner=inner,
                                               **FIX)

    def ev(x):
        n = x.shape[1]
        out = np.empty(n)
        for s in range(0, n, INTERP_CHUNK):
            out[s:s + INTERP_CHUNK] = pec_sphere_implant_transient_from_modes(
                x[:, s:s + INTERP_CHUNK].T, t, modes, rho_c=RHO_C, t_inf=T_INF) - T_INF
        return out
    return ev


def _rel(msh, num, den):
    dx = ufl.dx(domain=msh, metadata=RT["qmeta"])
    return float(np.sqrt(m7._allsum(msh.comm, num**2 * dx) / m7._allsum(msh.comm, den**2 * dx)))


def _fixture(res):
    msh, ct, _, t_build = _ball_curved(res)
    sub = m7._shell_submesh(msh, ct)
    deg = _cmap_degree(sub)
    assert deg == 2, f"curved hole sub-mesh cmap degree {deg}"
    sft = m7._robin_outer(sub, A, B)
    q_f = m7._source_dg2(sub, None, m7._q_shell, dg=RT["q_dg"])
    return sub, sft, q_f, t_build


def _transient(sub, sft, q_f, t_inf=T_INF, **kw):
    t0 = time.perf_counter()
    r = solve_transient_conduction(sub, sft, k=K_GEL, rho_c=RHO_C, q=q_f, h=H, t_inf=t_inf,
                                   t_init=t_inf, robin_tag=ROBIN_TAG, degree=2, **kw)
    return r, time.perf_counter() - t0


def _case(res):
    t0 = time.perf_counter()
    sub, sft, q_f, t_build = _fixture(res)
    nsub = sub.topology.index_map(3).size_global
    _say(f"\nIMP7-5 res={res} sub_cells={nsub} cmap={_cmap_degree(sub)} build={t_build:.1f}s; "
         f"rho c={RHO_C}, kappa={KAPPA:.6e}, t_s=(b-a)^2/kappa={T_S:.4f} s, t*={T_STAR:.4f} s, "
         f"CN dt=t*/{N_MAIN}={T_STAR / N_MAIN:.4f} s, DT_1={DT_1:.6f} s, n_terms={N_DEFAULT}")
    out = {"E": {}}

    # (T) CN N_MAIN steps with two BE half-steps, to t*.
    cn, tc = _transient(sub, sft, q_f, dt=T_STAR / N_MAIN, n_steps=N_MAIN, theta=0.5,
                        startup_half_steps=2)
    out["E"]["cn"] = cn.max_balance_residual
    tq = time.perf_counter()
    theta = m7._p2(sub, _series_fn(T_STAR, "neumann"), p=RT["comparand"])
    nan = sub.comm.allreduce(int(np.count_nonzero(~np.isfinite(theta.x.array))), op=MPI.SUM)
    assert nan == 0, f"transient comparand holds {nan} NaN dofs"
    e_all, e_in, _ = m7._errors(sub, None, cn.temperature, theta, B, qmeta=RT["qmeta"])
    out["T"], out["T_whole"] = e_in, e_all
    _say(f"IMP7-5 (T) res={res} CN N={N_MAIN} at t*: interior rel L2 = {e_in:.4e} "
         f"(INTERIOR_BAND {INTERIOR_BAND:.0e}); whole {e_all:.4e} (printed); "
         f"balance {cn.max_balance_residual:.3e}; solve {tc:.1f}s, "
         f"comparand {time.perf_counter() - tq:.1f}s")

    # Negative control: the Dirichlet-at-a series at t*.
    tq = time.perf_counter()
    th_d = m7._p2(sub, _series_fn(T_STAR, "dirichlet"), p=RT["comparand"])
    _, c_in, _ = m7._errors(sub, None, cn.temperature, th_d, B, qmeta=RT["qmeta"])
    out["cD"] = c_in
    _say(f"IMP7-5 control (c-D) res={res} Dirichlet-at-a series at t*: interior miss = "
         f"{c_in:.4e} (closed-form miss 0.431, item 658; asserted > INTERIOR_BAND "
         f"{INTERIOR_BAND:.0e}); {time.perf_counter() - tq:.1f}s")

    # (S-slope) one BE step of DT_1, solved in the θ variable (t_inf = t_init
    # = 0; item 672 ruling (i)): the storage integral carries no 20 K offset.
    be1, t1 = _transient(sub, sft, q_f, t_inf=0.0, dt=DT_1, n_steps=1, theta=1.0)
    out["E"]["be_1"] = be1.max_balance_residual
    stor = m7._allsum(sub.comm, be1.temperature * ufl.dx(domain=sub, metadata=RT["qmeta"]))
    ratio = RHO_C * stor / (DT_1 * be1.source_power)
    out["slope_dev"] = abs(ratio - 1.0)
    _say(f"IMP7-5 (S-slope) res={res} [theta, t_inf = 0] rho c int theta^1/(DT_1 int Q) = "
         f"{ratio:.12f}; |ratio - 1| = {out['slope_dev']:.4e} (SLOPE_BAND = MAT-12 S2_BAND "
         f"{SLOPE_BAND:.0e}; predicted ≈ h A_b DT_1 q_b/(rho c int Q) ≈ 5e-5); balance "
         f"{be1.max_balance_residual:.3e} (gated); solve {t1:.1f}s")
    # Record (printed, no gate): the same BE step in T at T∞ = 20 K.  Its
    # balance is the cancellation floor of a ≈ q·DT_1/ρc step on a 20 K level
    # (item 668: 1.300e-9 / 1.161e-9).
    beT, _ = _transient(sub, sft, q_f, dt=DT_1, n_steps=1, theta=1.0)
    storT = m7._allsum(sub.comm, (beT.temperature - T_INF)
                       * ufl.dx(domain=sub, metadata=RT["qmeta"]))
    ratio_T = RHO_C * storT / (DT_1 * beT.source_power)
    out["be_1_T"] = beT.max_balance_residual
    _say(f"IMP7-5 (E-record) res={res} BE DT_1 in T (T_inf = {T_INF}): balance "
         f"{beT.max_balance_residual:.3e} vs {be1.max_balance_residual:.3e} in theta; "
         f"|ratio - 1| = {abs(ratio_T - 1.0):.4e}; slope ratio T - theta = "
         f"{ratio_T - ratio:.3e} (printed)")
    # Discriminator (printed, ruling (iv)): BE in T at 2·DT_1.  Cancellation
    # predicts a balance ≈ 0.5 × the DT_1 one; the facet-quadrature channel ≈ 1.
    be2, _ = _transient(sub, sft, q_f, dt=2.0 * DT_1, n_steps=1, theta=1.0)
    out["be_2_T"] = be2.max_balance_residual
    out["disc_ratio"] = be2.max_balance_residual / beT.max_balance_residual
    _say(f"IMP7-5 (E-disc) res={res} BE in T: balance at 2 DT_1 = {be2.max_balance_residual:.3e}, "
         f"at DT_1 = {beT.max_balance_residual:.3e}; ratio = {out['disc_ratio']:.4f} "
         f"(cancellation predicts ≈ 0.5, facet quadrature ≈ 1; printed)")

    # (inf) BE 20 x 5 t_s against the steady FEM solve on the same mesh.
    bei, ti = _transient(sub, sft, q_f, dt=5.0 * T_S, n_steps=20, theta=1.0)
    out["E"]["be_inf"] = bei.max_balance_residual
    st, ts = m7._solve(sub, sft, K_GEL, q_f, degree=2)
    out["inf"] = _rel(sub, bei.temperature - st.temperature, st.temperature - T_INF)
    th_ss = m7._p2(sub, lambda x: m7._theta_shell(x, m7._coeffs(0.0)), p=RT["comparand"])
    inf_series = _rel(sub, (bei.temperature - T_INF) - th_ss, th_ss)
    _say(f"IMP7-5 (inf) res={res} BE 20 x 5 t_s vs steady FEM: rel L2 = {out['inf']:.3e} "
         f"(INF_BAND {INF_BAND:.0e}); vs steady series {inf_series:.4e} (printed; "
         f"the steady P2 discretisation error); balance {bei.max_balance_residual:.3e}; "
         f"solve {ti:.1f}s + steady {ts:.1f}s")

    worst = max(out["E"], key=out["E"].get)
    _say(f"IMP7-5 (E) res={res} per run " + ", ".join(f"{k}: {v:.3e}" for k, v in out["E"].items())
         + f"; worst {out['E'][worst]:.3e} on {worst} (E_BAND {E_BAND:.0e}); "
         f"elapsed {time.perf_counter() - t0:.1f}s")
    return out


def _gates(out, res):
    """Evaluate every gate, print the collected list, then assert (item 672 (iii))."""
    g = [(f"(E) {name}", r <= E_BAND, f"{r:.3e} <= {E_BAND:.0e}")
         for name, r in out["E"].items()]
    g += [
        ("(T)", out["T"] <= INTERIOR_BAND, f"{out['T']:.4e} <= {INTERIOR_BAND:.0e}"),
        ("(S-slope)", out["slope_dev"] <= SLOPE_BAND,
         f"{out['slope_dev']:.4e} <= {SLOPE_BAND:.0e}"),
        ("(inf)", out["inf"] <= INF_BAND, f"{out['inf']:.3e} <= {INF_BAND:.0e}"),
        ("(c-D)", out["cD"] > INTERIOR_BAND, f"{out['cD']:.4e} > {INTERIOR_BAND:.0e}"),
    ]
    _say(f"IMP7-5 gates res={res}: " + "; ".join(
        f"{n} {'PASS' if ok else 'FAIL'} [{m}]" for n, ok, m in g))
    bad = [f"{n} [{m}]" for n, ok, m in g if not ok]
    assert not bad, f"res={res} failing gates: " + ", ".join(bad)


# ---- tests --------------------------------------------------------------------


def test_imp7_5_w0_solver_accepts_curved_p2():
    """W0: `solve_transient_conduction` on a coordinate-degree-2 sub-mesh with
    P2 (the brief's check-first), two CN steps on the coarse 4 mm ball."""
    from fem_em_solver.io.mesh import MeshGenerator
    msh, ct, _ = MeshGenerator.two_layer_ball_domain(
        inner_radius=A, outer_radius=B, resolution=4e-3, geometry_order=2)
    sub = m7._shell_submesh(msh, ct)
    assert _cmap_degree(sub) == 2
    sft = m7._robin_outer(sub, A, B)
    q_f = m7._source_dg2(sub, None, m7._q_shell, dg=RT["q_dg"])
    r, t = _transient(sub, sft, q_f, dt=T_STAR / N_MAIN, n_steps=2, theta=0.5,
                      startup_half_steps=2)
    deg = r.temperature.function_space.ufl_element().basix_element.degree
    _say(f"IMP7-5 W0: sub cmap={_cmap_degree(sub)} cells={sub.topology.index_map(3).size_global} "
         f"T degree={deg}; 2 steps balance {r.max_balance_residual:.3e} (E_BAND {E_BAND:.0e}); "
         f"{t:.1f}s; t*={T_STAR:.4f} s")
    assert deg == 2
    assert r.max_balance_residual <= E_BAND


def test_imp7_5_w1_ladder0():
    """W1: LADDER[0]."""
    _gates(_case(LADDER[0]), LADDER[0])


def test_imp7_5_w2_ladder1():
    """W2: LADDER[1]."""
    _gates(_case(LADDER[1]), LADDER[1])
