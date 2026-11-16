"""`MAT-21` step 2, FEM — the F2182 gel box on P25's finest rung (m = 6) under
`MAT-17` step 3's duty cycle re-timed to t* = 0.2 t_C, read through `MAT-20`'s
operators (§9 item 424, §10 P26; brief `docs/planning/chunks/MAT-21.md`);
attempt 2, §9 item 439, the 2026-10-13 18:00 daily's ruling (provisional).

Schedule (W1 module's pins): N = 40 CN steps of DT_BOX = t*/40 with two BE
start-up half-steps, HALF_T = 4 DT_BOX, χ = 1 on [jP, jP + P/2); switches
t_m = m·HALF_T, m = 0 … 9, all grid points.  The schedule is written locally
(``_chi_box`` / ``_heav_box``): ``test_mat17_duty_cycle_coil_driven``'s
``HALF`` / ``DT`` / ``_chi_duty`` / ``_heaviside`` carry `MAT-15`'s t*.

W2 (fixture ``runs``: duty, ten Heaviside H_m, χ ≡ 1, unscheduled default —
13 runs).  Asserted, bands imported: (P0′) the pins equal the originals and
T_STAR = 0.2 t_C; (S0) Σχₙ Δtₙ = D t* ≤ 1e-12 (``_census``); (Δ) θ_duty =
Σ(−1)^m θ_Hm at t* ≤ ``L_BAND``; (D′) χ ≡ 1 = default ≤ ``D_BAND``; (H0) the
per-step mean series obeys Duhamel ≤ ``L_BAND``; (FEM′) / (FEM′-2d)
``running_average`` commutes with the Duhamel sum ≤ ``D_BAND`` once the
per-step residual r = duty − Σ(−1)^m H_m is carried through it
(|a − b − ra(r)| / max|a|), and (FEM-c) the window-mean bound
max|a − b| ≤ max|r|(1 + 1e-12) + ``D_BAND``·max|a| (K, 1-D and 2-D); (FEM)
as posed (= (FEM-ctl)) is printed beside its ceiling max|r|/max|a|; (E) every run's
balance ≤ ``E_BAND``.  Controls (rule (e), predicted, only the miss
asserted): (Δ-ctl) m = 9 dropped > ``L_BAND``; (Wc-FEM) P/2 shift > ``D_BAND``.

W3 (fixture ``duty_cf``: the duty run alone).  (T-duty) interior rel L2 at t*
against the ten-run triple-series sum in P2 ≤ ``INTERIOR_BAND``; control (c3)
the unscheduled step response as comparand misses ``INTERIOR_BAND``.
"""

from __future__ import annotations

import time

import numpy as np
import pytest
import ufl
from dolfinx import fem
from mpi4py import MPI

from fem_em_solver.post.evaluation import evaluate_vector_field_parallel
from fem_em_solver.thermal import solve_transient_conduction
from fem_em_solver.thermal.dose import cem43, running_average

from tests.complex_mode import complex_only
from tests.validation.test_mat10_sar_to_heat_sphere import ROBIN_TAG, _allsum, _say
from tests.validation.test_mat12_transient_sar_heating_sphere import E_BAND, RHO_C
from tests.validation.test_mat13_transient_perfusion_sphere import INTERIOR_BAND
from tests.validation.test_mat15_coil_driven_heating_phantom import H_GEL, K_GEL, N_MAIN, T_INF
from tests.validation.test_mat15_finite_cylinder_fem import _max_imag, _rel
from tests.validation.test_mat17_duty_cycle_coil_driven import (
    D_BAND,
    DUTY,
    L_BAND,
    N_SWITCH,
    STEPS_PER_HALF,
)
from tests.validation.test_mat17_duty_cycle_fem import _census
from tests.validation.test_mat21_box_closed_form import Q, T_C
from tests.validation.test_mat21_box_closed_form import T_INF as T_INF_0
from tests.validation.test_mat21_box_duty_closed_form import (
    DT_BOX,
    HALF_T,
    N_STEPS,
    QV,
    T_STAR,
    W_AVG,
    _box_mean_duty,
    _duty_cf,
)
from tests.validation.test_mat21_box_duty_closed_form import DUTY as DUTY_PIN
from tests.validation.test_mat21_box_duty_closed_form import N_SWITCH as N_SWITCH_PIN
from tests.validation.test_mat21_box_duty_closed_form import STEPS_PER_HALF as SPH_PIN
from tests.validation.test_mat21_box_fem import (
    _Clock,
    _exact,
    _interior_chi,
    _rung,
    _trans_ex,
)
from tests.validation.test_mat21_box_fem import T_STAR as T_STAR_1

STEPS = tuple(range(N_MAIN + 1))
M_RUNG = 6


def _chi_box(t):
    """1 on [jP, jP + P/2), 0 otherwise; t is always a (half-)step midpoint."""
    return 1.0 if int(np.floor(t / HALF_T)) % 2 == 0 else 0.0


def _heav_box(m):
    t_m = m * HALF_T
    return lambda t: 1.0 if t > t_m else 0.0


def _cn(msh, ft, schedule):
    kw = {} if schedule is None else {"source_schedule": schedule}
    return solve_transient_conduction(
        msh, ft, k=K_GEL, rho_c=RHO_C, q=Q, h=H_GEL, t_inf=T_INF, t_init=T_INF,
        robin_tag=ROBIN_TAG, dt=DT_BOX, n_steps=N_MAIN, theta=0.5,
        startup_half_steps=2, degree=1, snapshot_steps=STEPS, **kw)


def _series(res, work, f_mean, comm, v_h):
    """`MAT-20`'s ``_series`` (``test_mat20_fem_history.py:60–67``)."""
    out = []
    for n in STEPS:
        work.x.array[:] = res.snapshots[n].x.array - T_INF
        work.x.scatter_forward()
        out.append(float(np.real(comm.allreduce(fem.assemble_scalar(f_mean), op=MPI.SUM))) / v_h)
    res.snapshots.clear()
    return np.asarray(out)


def _identities(t, duty, hs):
    """`MAT-20`'s ``_identities`` (``test_mat20_fem_history.py:124–150``), on W_AVG."""
    signs = [(-1.0) ** m for m in range(N_SWITCH)]
    sup = sum(s * h for s, h in zip(signs, hs))
    r = duty - sup  # the one per-step Duhamel residual: (H0), (FEM'), (FEM-c), the ceiling
    r_max = float(np.max(np.abs(r)))
    h0 = r_max / float(np.max(np.abs(duty)))
    a, a_mx, a_tmx = running_average(t, duty, window=W_AVG)
    ra = [running_average(t, h, window=W_AVG) for h in hs]
    b = sum(s * x[0] for s, x in zip(signs, ra))  # from the ten H_m runs, never from a
    rr = running_average(t, r, window=W_AVG)[0]
    ok = ~np.isnan(a)
    scale = np.max(np.abs(a[ok]))
    res = np.abs(a[ok] - b[ok]) / scale
    fem_v = float(np.max(res))  # (FEM) as posed = (FEM-ctl): printed only
    res_p = np.abs(a[ok] - b[ok] - rr[ok]) / scale
    fem_p = float(np.max(res_p))  # (FEM')
    fem_c_lhs = float(np.max(np.abs(a[ok] - b[ok])))  # (FEM-c), K
    arr = np.column_stack([duty] + hs + [r])
    m2, _, _ = running_average(t, arr, window=W_AVG)
    b2 = m2[:, 1:1 + N_SWITCH] @ np.asarray(signs)
    rr2 = m2[:, -1]
    scale2 = np.max(np.abs(m2[ok, 0]))
    fem2d = float(np.max(np.abs(m2[ok, 0] - b2[ok])) / scale2)
    fem2d_p = float(np.max(np.abs(m2[ok, 0] - b2[ok] - rr2[ok])) / scale2)
    fem_c_lhs2 = float(np.max(np.abs(m2[ok, 0] - b2[ok])))
    idx = np.nonzero(ok)[0]
    idx = idx[idx - STEPS_PER_HALF >= idx[0]]
    wc = float(np.max(np.abs(a[idx] - b[idx - STEPS_PER_HALF])) / scale)
    return {"h0": h0, "a": a, "a_mx": a_mx, "a_tmx": a_tmx, "ok": ok, "res": res,
            "fem": fem_v, "fem2d": fem2d, "idx": idx, "wc": wc,
            "r_max": r_max, "scale": float(scale), "scale2": float(scale2),
            "res_p": res_p, "fem_p": fem_p, "fem2d_p": fem2d_p,
            "fem_c_lhs": fem_c_lhs, "fem_c_lhs2": fem_c_lhs2}


@pytest.fixture(scope="module")
def runs():
    comm = MPI.COMM_WORLD
    t0 = time.perf_counter()
    clk = _Clock("step2-runs")
    msh, ft, ns, n, v_h, _st = _rung(M_RUNG, clk)
    V1 = fem.functionspace(msh, ("Lagrange", 1))
    work = fem.Function(V1)
    f_mean = fem.form(work * ufl.dx(domain=msh))
    _say(f"MAT21-2 [runs] rung m = {M_RUNG} {tuple(ns)}: cells {n}; P1 dofs "
         f"{V1.dofmap.index_map.size_global}; t* = {T_STAR:.4f} s, DT_BOX = {DT_BOX:.6f} s, "
         f"HALF_T = {HALF_T:.4f} s, W/DT = {W_AVG / DT_BOX:.4f} samples")
    out = {"E": {}, "imag": 0.0, "H": [], "H_fields": []}

    def _go(label, schedule):
        res = _cn(msh, ft, schedule)
        clk.lap(f"CN N={N_MAIN} {label}")
        out["E"][label] = res.max_balance_residual
        out["imag"] = max(out["imag"], _max_imag(res.temperature))
        return res

    r = _go("duty", _chi_box)
    out["mult"] = list(r.source_multipliers)
    _census(r, DT_BOX, 2, DUTY * T_STAR, "box duty")
    out["census_ok"] = True
    out["duty_field"] = r.temperature
    out["duty"] = _series(r, work, f_mean, comm, v_h)
    for m in range(N_SWITCH):
        r = _go(f"H{m}", _heav_box(m))
        out["H_fields"].append(r.temperature)
        out["H"].append(_series(r, work, f_mean, comm, v_h))
    r = _go("chi1", lambda t: 1.0)
    out["one_field"] = r.temperature
    r.snapshots.clear()
    r = _go("default", None)
    out["default_field"] = r.temperature
    r.snapshots.clear()
    clk.lap("13 runs + reductions")

    # (Δ) and (Δ-ctl) on fields at t*.
    sig = fem.Function(V1)
    sig.x.array[:] = sum((-1.0) ** m * (f.x.array - T_INF) for m, f in enumerate(out["H_fields"]))
    sig9 = fem.Function(V1)
    sig9.x.array[:] = sig.x.array - (-1.0) ** 9 * (out["H_fields"][9].x.array - T_INF)
    th_duty = out["duty_field"] - T_INF
    out["delta"] = _rel(msh, th_duty, sig)
    out["delta_ctl"] = _rel(msh, th_duty, sig9)
    out["dprime"] = _rel(msh, out["one_field"] - T_INF, out["default_field"] - T_INF)
    out["t"] = np.arange(N_MAIN + 1) * DT_BOX
    out["ids"] = _identities(out["t"], out["duty"], out["H"])
    out["msh"] = msh
    _say(f"MAT21-2 [runs] max|Im theta| = {out['imag']:.3e} (predicted 0); elapsed "
         f"{time.perf_counter() - t0:.1f} s")
    return out


def test_precondition_schedule_shape():
    assert (N_STEPS, SPH_PIN, N_SWITCH_PIN, DUTY_PIN) == (N_MAIN, STEPS_PER_HALF, N_SWITCH, DUTY)
    assert T_STAR == 0.2 * T_C
    assert T_STAR == T_STAR_1
    assert T_INF == T_INF_0
    assert HALF_T == STEPS_PER_HALF * DT_BOX


@complex_only
def test_s0_schedule_census(runs):
    s = "".join(str(int(x)) for x in runs["mult"])
    _say(f"MAT21-2 (S0) schedule (2 half-steps, then 39 steps): {s}")
    assert runs["census_ok"]


@complex_only
def test_delta_duhamel_fields(runs):
    _say(f"MAT21-2 (Δ) whole-box rel L2 |theta_duty - sum (-1)^m theta_Hm| at t* = "
         f"{runs['delta']:.3e} (bound L_BAND {L_BAND:.0e}, predicted <= 1e-13)")
    assert runs["delta"] <= L_BAND


@complex_only
def test_d_prime_duty_one_is_default(runs):
    _say(f"MAT21-2 (D') chi = 1 vs default at t*: rel L2 = {runs['dprime']:.3e} "
         f"(bound D_BAND {D_BAND:.0e}, predicted <= 5e-15)")
    assert runs["dprime"] <= D_BAND


@complex_only
def test_h0_series_duhamel(runs):
    h0 = runs["ids"]["h0"]
    t, duty = runs["t"], runs["duty"]
    cf = np.array([_box_mean_duty(tn) for tn in t])
    nz = np.abs(cf) > 0
    dev = np.max(np.abs(duty[nz] - cf[nz]) / np.abs(cf[nz]))
    _say("MAT21-2 record <theta_duty>_n FEM vs closed-form mean: "
         + "; ".join(f"n{i} {duty[i]:.6f}/{cf[i]:.6f}" for i in range(len(t))))
    _say(f"MAT21-2 record max_n |FEM mean - closed-form mean|/closed-form = {dev:.3e} "
         f"(predicted 1e-4 to 3e-3; printed, not asserted)")
    _say(f"MAT21-2 (H0) max_n |duty - sum (-1)^m H_m| / max|duty| = {h0:.3e} "
         f"(bound L_BAND {L_BAND:.0e}, predicted <= 1e-13)")
    assert h0 <= L_BAND


@complex_only
def test_fem_six_minute_operator(runs):
    ids = runs["ids"]
    t, duty = runs["t"], runs["duty"]
    nw = int(ids["ok"].sum())
    ceil_ = ids["r_max"] / ids["scale"]
    rhs = ids["r_max"] * (1.0 + 1e-12) + D_BAND * ids["scale"]
    rhs2 = ids["r_max"] * (1.0 + 1e-12) + D_BAND * ids["scale2"]
    _say(f"MAT21-2 (FEM) as posed, printed only: W = {W_AVG} s, max |a - b| / max|a| over {nw} "
         f"windows = {ids['fem']:.3e}; 2-D twin {ids['fem2d']:.3e}; ceiling r_max/max|a| = "
         f"{ceil_:.3e} (predicted ≈ 1.46e-13 against ≈ 3.85e-13)")
    _say(f"MAT21-2 (FEM-ctl) = (FEM) as posed, r omitted, predicted ≈ 1.46e-13, not asserted: "
         f"{ids['fem']:.3e}")
    _say("MAT21-2 (FEM) per-window residuals: " + " ".join(f"{x:.2e}" for x in ids["res"]))
    _say("MAT21-2 (FEM') per-window residuals: " + " ".join(f"{x:.2e}" for x in ids["res_p"]))
    _say(f"MAT21-2 (FEM') max |a - b - rr| / max|a| = {ids['fem_p']:.3e}; (FEM'-2d) "
         f"{ids['fem2d_p']:.3e} (bound D_BAND {D_BAND:.0e}, predicted <= 1e-14, worst-case ≈ 3e-14)")
    _say(f"MAT21-2 (FEM-c) max |a - b| = {ids['fem_c_lhs']:.4e} K <= r_max(1 + 1e-12) + D_BAND max|a| "
         f"= {rhs:.4e} K (r_max {ids['r_max']:.4e} K); 2-D {ids['fem_c_lhs2']:.4e} K <= {rhs2:.4e} K "
         f"(predicted 9.03e-14 <= 2.99e-13 K)")
    _say(f"MAT21-2 record FEM six-minute maximum {ids['a_mx']:.5f} K at t_max = {ids['a_tmx']:.4f} s "
         f"(predicted ≈ 0.6168 K at t* = {T_STAR:.4f} s)")
    d_w = cem43(t, T_INF + np.real(duty) / QV)
    _say(f"MAT21-2 record FEM cem43 of T_inf + <theta>/(qV) = {d_w:.4e} s per W of qV "
         f"(W1 closed-form mean 4.0494e-11 s, 20261013T183409Z_MAT-21-step2-W1.log)")
    assert ids["fem_p"] <= D_BAND, f"(FEM') {ids['fem_p']:.3e}"
    assert ids["fem2d_p"] <= D_BAND, f"(FEM'-2d) {ids['fem2d_p']:.3e}"
    assert ids["fem_c_lhs"] <= rhs, f"(FEM-c) {ids['fem_c_lhs']:.4e} > {rhs:.4e} K"
    assert ids["fem_c_lhs2"] <= rhs2, f"(FEM-c 2-D) {ids['fem_c_lhs2']:.4e} > {rhs2:.4e} K"


@complex_only
def test_e_every_run_balances(runs):
    worst = max(runs["E"], key=runs["E"].get)
    _say(f"MAT21-2 (E) per run " + ", ".join(f"{k_}: {v:.3e}" for k_, v in runs["E"].items())
         + f"; worst {runs['E'][worst]:.3e} on {worst} (bound {E_BAND:.0e}, predicted <= 1e-11)")
    for name, r in runs["E"].items():
        assert r <= E_BAND, f"(E) {name}: {r:.3e}"


@complex_only
def test_control_delta_m9_dropped_misses(runs):
    _say(f"MAT21-2 control (Δ-ctl) m = 9 dropped at t*: miss = {runs['delta_ctl']:.4e} "
         f"(predicted 0.20-0.22, range 0.18-0.26; asserted > L_BAND)")
    assert runs["delta_ctl"] > L_BAND


@complex_only
def test_control_wc_shifted_window_misses(runs):
    ids = runs["ids"]
    _say(f"MAT21-2 control (Wc-FEM) shift P/2 = {STEPS_PER_HALF} samples: miss {ids['wc']:.4e} "
         f"over {ids['idx'].size} windows (predicted 0.141, range 0.10-0.18; asserted > D_BAND)")
    assert ids["wc"] > D_BAND


@pytest.fixture(scope="module")
def duty_cf():
    t0 = time.perf_counter()
    clk = _Clock("step2-cf")
    msh, ft, ns, n, v_h, _st = _rung(M_RUNG, clk)
    res = _cn(msh, ft, _chi_box)
    clk.lap(f"CN N={N_MAIN} duty")
    res.snapshots.clear()
    th = res.temperature - T_INF
    chi = _interior_chi(msh)
    ex, npts = _exact(msh, lambda x, y, z: _duty_cf(x, y, z, T_STAR))
    clk.lap(f"ten-run series interpolation at t* ({npts} P2 points)")
    out = {"T": _rel(msh, th, ex, chi), "whole": _rel(msh, th, ex)}
    ex_s, _ = _exact(msh, _trans_ex(T_STAR))
    clk.lap("step-response interpolation (c3)")
    out["c3"] = _rel(msh, th, ex_s, chi)
    try:
        vals, ok = evaluate_vector_field_parallel(res.temperature, np.zeros((1, 3)))
        th_c = float(np.real(np.ravel(vals[0])[0])) - T_INF if bool(ok[0]) else float("nan")
    except Exception as exc:  # printed quantity only
        _say(f"MAT21-2 centre evaluation failed: {exc!r}")
        th_c = float("nan")
    th_c_ex = float(np.real(_duty_cf(np.zeros(1), np.zeros(1), np.zeros(1), T_STAR)[0]))
    _say(f"MAT21-2 (T-duty) interior rel L2 at t* = {out['T']:.6e} (bound {INTERIOR_BAND:.0e}, "
         f"predicted 5e-4 to 3e-3; step 1's (T) 9.072e-4); whole {out['whole']:.6e} (printed); "
         f"theta_h(0) = {th_c:.6f} K vs ten-run series {th_c_ex:.6f} K (predicted ≈ 0.6665); "
         f"balance {res.max_balance_residual:.3e}; max|Im| {_max_imag(res.temperature):.3e}; "
         f"elapsed {time.perf_counter() - t0:.1f} s")
    _say(f"MAT21-2 control (c3) unscheduled step response as comparand: interior miss = "
         f"{out['c3']:.4e} (predicted 0.50, range 0.45-0.55; asserted > {INTERIOR_BAND:.0e})")
    return out


@complex_only
def test_t_duty_interior_against_closed_form(duty_cf):
    assert duty_cf["T"] <= INTERIOR_BAND, f"(T-duty) {duty_cf['T']:.4e} > 1 %"


@complex_only
def test_control_c3_step_response_misses(duty_cf):
    assert duty_cf["c3"] > INTERIOR_BAND, f"(c3) {duty_cf['c3']:.4e} inside 1 %"
