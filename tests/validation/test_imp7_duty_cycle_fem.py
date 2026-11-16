"""`IMP-7` step 6 (§9 item 675, §10 I47): the duty-cycled Crank–Nicolson
transient on I43's curved P2 adiabatic hole against I45's duty series.

Brief: docs/planning/chunks/IMP-7.md, "§9 item 675 — 2026-10-27".  Fixture,
transient driver and the rel-L2 helper are imported from item 672's module
(`tests/validation/test_imp7_transient_fem.py`: `_fixture`, `_transient`,
`_rel`, with its constants), so the gel shell a < r < b, the natural condition
on r = a, Robin h on r = b, the β = 1 dipole SAR (DG3) and P2 on the
geometry-order-2 sub-mesh are item 672's exactly.  The source is switched by
χ(t) = 1 if (t mod period) < duty·period else 0 — `MAT-17`'s `_chi`
(`test_mat17_duty_cycle_fem.py:80–82`), copied, not imported — with I45's
(Duh) parameters period = 0.5 t*, duty 0.5, so t* spans two periods.  CN with
two BE start-up half-steps, dt = t*/N_MAIN (40 steps: a period is 20 steps, an
on-phase 10, so every switch lies on a step boundary), at T∞ = 20 K so the run
exercises the T∞ wiring (item 672 ruling (ii)).

Comparand: `pec_sphere_implant_duty_cycle_from_modes` (`analytical.py:3686`)
on `pec_sphere_implant_transient_modes(beta=BETA, n_terms=N_DEFAULT,
inner="neumann", **FIX)`, interpolated in chunks of INTERP_CHUNK points into
P3, quadrature degree 6.  W0 cross-checks it (numpy) at five axis radii
against I45's Duhamel sum of shifted step responses, ≤ B1_BAND.

Gates, imported and unmoved, collected then asserted: (T) interior rel L2 at
t* ≤ INTERIOR_BAND; (E) max balance ≤ E_BAND (the solver normalises by the
full-on ∫Q, so the off steps are defined); (S0) |Σχ·dt − on-time|/on-time
≤ 1e-12 (`MAT-17`'s census).  Negative control (asserted): the continuous-
source series at t* misses the duty FEM by > INTERIOR_BAND.  Printed only:
the no-startup CN run.  Real build.
"""

from __future__ import annotations

import time

import numpy as np
from mpi4py import MPI

from fem_em_solver.utils.analytical import (
    pec_sphere_implant_duty_cycle_from_modes,
    pec_sphere_implant_transient_modes,
)
from tests.validation import test_imp7_implant_heating_fem as m7
from tests.validation.test_imp4_pec_cylinder_closed_form import B1_BAND
from tests.validation.test_imp7_implant_heating_closed_form import BETA, FIX, T_INF
from tests.validation.test_imp7_transient_closed_form import (
    N_DEFAULT,
    T_STAR,
    _modes,
    _theta,
)
from tests.validation.test_imp7_transient_fem import (
    INTERP_CHUNK,
    LADDER,
    RT,
    _fixture,
    _series_fn,
    _transient,
)
from tests.validation.test_mat12_transient_sar_heating_sphere import E_BAND, N_MAIN, RHO_C
from tests.validation.test_mat13_transient_perfusion_sphere import INTERIOR_BAND

A, B = m7.A, m7.B
PERIOD = 0.5 * T_STAR  # I45 (Duh): period 0.5 t*, duty 0.5
DUTY = 0.5
ON_TIME = 2 * DUTY * PERIOD  # t* = two periods
DT = T_STAR / N_MAIN
S0_BAND = 1e-12  # MAT-17's census bar (test_mat17_duty_cycle_fem.py:116)


def _say(msg):
    if MPI.COMM_WORLD.rank == 0:
        print(msg, flush=True)


def _chi(period, duty, shift=0.0):
    def chi(t):
        return 1.0 if ((t - shift) % period) < duty * period else 0.0

    return chi


def _duty_fn(t):
    modes = pec_sphere_implant_transient_modes(beta=BETA, n_terms=N_DEFAULT, inner="neumann",
                                               **FIX)

    def ev(x):
        n = x.shape[1]
        out = np.empty(n)
        for s in range(0, n, INTERP_CHUNK):
            out[s:s + INTERP_CHUNK] = pec_sphere_implant_duty_cycle_from_modes(
                x[:, s:s + INTERP_CHUNK].T, t, modes, rho_c=RHO_C, t_inf=T_INF,
                period=PERIOD, duty=DUTY) - T_INF
        return out
    return ev


def _census(s, startup):
    dts = [0.5 * DT] * startup + [DT] * (len(s.source_multipliers) - startup)
    tot = float(np.dot(np.asarray(s.source_multipliers), np.asarray(dts)))
    n_on = int(sum(1 for x in s.source_multipliers if x == 1.0))
    n_off = int(sum(1 for x in s.source_multipliers if x == 0.0))
    return abs(tot - ON_TIME) / ON_TIME, tot, n_on, n_off


def _case(res):
    t0 = time.perf_counter()
    sub, sft, q_f, t_build = _fixture(res)
    nsub = sub.topology.index_map(3).size_global
    _say(f"\nIMP7-6 res={res} sub_cells={nsub} build={t_build:.1f}s; t*={T_STAR:.4f} s, "
         f"period={PERIOD:.4f} s, duty={DUTY}, on-time={ON_TIME:.4f} s, CN dt=t*/{N_MAIN}="
         f"{DT:.4f} s (period {PERIOD / DT:.6f} steps), n_terms={N_DEFAULT}")
    out = {}
    chi = _chi(PERIOD, DUTY)

    cn, tc = _transient(sub, sft, q_f, dt=DT, n_steps=N_MAIN, theta=0.5, startup_half_steps=2,
                        source_schedule=chi)
    out["E"] = cn.max_balance_residual
    out["S0"], tot, n_on, n_off = _census(cn, 2)
    _say(f"IMP7-6 (S0) res={res}: sum chi*dt = {tot:.12f} s vs on-time {ON_TIME:.12f} s, "
         f"rel {out['S0']:.3e} (bar {S0_BAND:.0e}); (half-)steps on={n_on} off={n_off} of "
         f"{len(cn.source_multipliers)}; chi = {cn.source_multipliers}")
    tq = time.perf_counter()
    th = m7._p2(sub, _duty_fn(T_STAR), p=RT["comparand"])
    nan = sub.comm.allreduce(int(np.count_nonzero(~np.isfinite(th.x.array))), op=MPI.SUM)
    assert nan == 0, f"duty comparand holds {nan} NaN dofs"
    e_all, e_in, _ = m7._errors(sub, None, cn.temperature, th, B, qmeta=RT["qmeta"])
    out["T"] = e_in
    _say(f"IMP7-6 (T) res={res} CN duty N={N_MAIN} at t*: interior rel L2 = {e_in:.4e} "
         f"(INTERIOR_BAND {INTERIOR_BAND:.0e}); whole {e_all:.4e} (printed); balance "
         f"{cn.max_balance_residual:.3e} (E_BAND {E_BAND:.0e}); solve {tc:.1f}s, comparand "
         f"{time.perf_counter() - tq:.1f}s")

    # Negative control: the continuous-source series (χ ≡ 1) at t*.
    tq = time.perf_counter()
    th_c = m7._p2(sub, _series_fn(T_STAR, "neumann"), p=RT["comparand"])
    _, c_in, _ = m7._errors(sub, None, cn.temperature, th_c, B, qmeta=RT["qmeta"])
    out["cC"] = c_in
    _say(f"IMP7-6 control (c-cont) res={res} continuous-source series at t*: interior miss = "
         f"{c_in:.4e} (asserted > INTERIOR_BAND {INTERIOR_BAND:.0e}; predicted ≫ 1e-2); "
         f"{time.perf_counter() - tq:.1f}s")

    # Printed only: CN without the start-up half-steps (MAT-17's s_ns).
    ns, tn = _transient(sub, sft, q_f, dt=DT, n_steps=N_MAIN, theta=0.5, startup_half_steps=0,
                        source_schedule=chi)
    _, e_ns, _ = m7._errors(sub, None, ns.temperature, th, B, qmeta=RT["qmeta"])
    s0_ns = _census(ns, 0)[0]
    _say(f"IMP7-6 printed res={res}: CN no-startup interior rel L2 = {e_ns:.4e}, balance "
         f"{ns.max_balance_residual:.3e}, census rel {s0_ns:.3e}; solve {tn:.1f}s; "
         f"elapsed {time.perf_counter() - t0:.1f}s")
    return out


def _gates(out, res):
    g = [
        ("(T)", out["T"] <= INTERIOR_BAND, f"{out['T']:.4e} <= {INTERIOR_BAND:.0e}"),
        ("(E)", out["E"] <= E_BAND, f"{out['E']:.3e} <= {E_BAND:.0e}"),
        ("(S0)", out["S0"] <= S0_BAND, f"{out['S0']:.3e} <= {S0_BAND:.0e}"),
        ("(c-cont)", out["cC"] > INTERIOR_BAND, f"{out['cC']:.4e} > {INTERIOR_BAND:.0e}"),
    ]
    _say(f"IMP7-6 gates res={res}: " + "; ".join(
        f"{n} {'PASS' if ok else 'FAIL'} [{m}]" for n, ok, m in g))
    bad = [f"{n} [{m}]" for n, ok, m in g if not ok]
    assert not bad, f"res={res} failing gates: " + ", ".join(bad)


# ---- tests --------------------------------------------------------------------


def test_imp7_6_w0_duty_helper_vs_duhamel():
    """W0 (numpy): the duty helper at t* against I45's Duhamel sum of shifted
    step responses at five axis radii, ≤ B1_BAND."""
    r = np.linspace(A * (1.0 + 1e-9), B, 5)
    pts = np.column_stack([np.zeros(5), np.zeros(5), r])
    modes = _modes()
    t, t_on = T_STAR, DUTY * PERIOD
    direct = pec_sphere_implant_duty_cycle_from_modes(
        pts, t, modes, rho_c=RHO_C, t_inf=T_INF, period=PERIOD, duty=DUTY) - T_INF
    duh = np.zeros(len(pts))
    k = 0
    while k * PERIOD < t:
        duh += _theta(pts, t - k * PERIOD, modes)
        if t - k * PERIOD - t_on > 0.0:
            duh -= _theta(pts, t - k * PERIOD - t_on, modes)
        k += 1
    rel = np.max(np.abs(duh - direct)) / np.max(np.abs(direct))
    cont = _theta(pts, t, modes)
    _say(f"IMP7-6 W0 axis r = {r}: duty θ = {direct}; Duhamel = {duh}; {k} periods; "
         f"max|Duhamel − modal|/max = {rel:.3e} (B1_BAND {B1_BAND:.0e}); continuous θ = "
         f"{cont} (printed)")
    assert rel <= B1_BAND


def test_imp7_6_w1_ladder0():
    """W1: LADDER[0]."""
    _gates(_case(LADDER[0]), LADDER[0])


def test_imp7_6_w2_ladder1():
    """W2: LADDER[1]."""
    _gates(_case(LADDER[1]), LADDER[1])
