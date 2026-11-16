"""`MAT-21` step 4 — step 2's duty-cycled box at Lagrange degree 2 against a P3
comparand (§9 item 451, §10 chain P29, daily 2026-10-14 — provisional; brief
`docs/planning/chunks/MAT-21.md`).

A sibling of step 2's ``test_mat21_box_duty_fem.py`` (no degree knob; its
``_cn`` hard-codes degree 1 and is never imported).  This module imports step
2's degree-free helpers (``_chi_box``, ``_heav_box``, ``_series``,
``_identities``) and step 3's degree-2 helpers (``_rung_p``, ``_exact_p``,
``_rel_p``, ``_dofs_pred``).  It sets no band and moves none: every band is
imported by name.  Every P2 reading goes through ``_rel_p`` (quadrature degree
6) against a P3 comparand.

W2 (fixture ``runs_p``, m = 3, P2, 13 runs): (P0″) pins; (Pdeg) global dofs;
(S0) census ≤ 1e-12; (Δ) and (H0) ≤ ``L_BAND``; (D′), (FEM′), (FEM′-2d) ≤
``D_BAND``; (FEM-c); (E) ≤ ``E_BAND``.  Controls (separations only):
(Δ-ctl) m = 9 dropped > ``L_BAND``; (Wc-FEM) P/2 shift > ``D_BAND``.

W3 (fixture ``duty_cf_p``, m = 4): (Pdeg) P2 and P1 dofs; (T-duty) interior
≤ ``INTERIOR_BAND``; controls (c3) step-response comparand > ``INTERIOR_BAND``
and (c-deg-duty) the degree-1 duty run on the same mesh, same P3 comparand
and χ, reads worse (e₁ > e₂; label *predicted*).
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
from fem_em_solver.thermal.dose import cem43

from tests.complex_mode import complex_only
from tests.validation.test_mat10_sar_to_heat_sphere import ROBIN_TAG, _allsum, _say  # noqa: F401
from tests.validation.test_mat12_transient_sar_heating_sphere import E_BAND, RHO_C
from tests.validation.test_mat13_transient_perfusion_sphere import INTERIOR_BAND
from tests.validation.test_mat15_coil_driven_heating_phantom import H_GEL, K_GEL, N_MAIN, T_INF
from tests.validation.test_mat15_finite_cylinder_fem import _max_imag
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
from tests.validation.test_mat21_box_duty_fem import (
    STEPS,
    _chi_box,
    _heav_box,
    _identities,
    _series,
)
from tests.validation.test_mat21_box_fem import _interior_chi, _trans_ex
from tests.validation.test_mat21_box_fem import T_STAR as T_STAR_1
from tests.validation.test_mat21_box_fem_p2 import (
    CMP_DEGREE,
    DEGREE,
    M_T,
    QDEG_P,
    RUNGS_P,
    _dofs_pred,
    _exact_p,
    _rel_p,
    _rung_p,
)

M_ID = RUNGS_P[1]  # 3, the identities rung
DEG_CTL = 1


def _cn_p(msh, ft, schedule, degree):
    """Step 2's ``_cn`` (``test_mat21_box_duty_fem.py:96–101``) with ``degree=degree``."""
    kw = {} if schedule is None else {"source_schedule": schedule}
    return solve_transient_conduction(
        msh, ft, k=K_GEL, rho_c=RHO_C, q=Q, h=H_GEL, t_inf=T_INF, t_init=T_INF,
        robin_tag=ROBIN_TAG, dt=DT_BOX, n_steps=N_MAIN, theta=0.5,
        startup_half_steps=2, degree=degree, snapshot_steps=STEPS, **kw)


class _Clock4:
    """Per-phase wall clock, printed on rank 0 with flush, prefix MAT21-4."""

    def __init__(self, tag):
        self.tag = tag
        self.t = time.perf_counter()

    def lap(self, phase):
        now = time.perf_counter()
        _say(f"MAT21-4 [{self.tag}] clock {phase}: {now - self.t:.1f} s")
        self.t = now


@pytest.fixture(scope="module")
def runs_p():
    comm = MPI.COMM_WORLD
    t0 = time.perf_counter()
    clk = _Clock4("runs")
    msh, ft, ns, n, v_h, st = _rung_p(M_ID, clk, DEGREE)
    V2 = st.temperature.function_space
    work = fem.Function(V2)
    f_mean = fem.form(work * ufl.dx(domain=msh))
    ndof = int(V2.dofmap.index_map.size_global)
    _say(f"\nMAT21-4 [runs] DEGREE = {DEGREE}, comparand P{CMP_DEGREE}, quadrature {QDEG_P}; "
         f"rung m = {M_ID} {tuple(ns)}: cells {n}; P2 dofs {ndof} (predicted {_dofs_pred(M_ID)}); "
         f"t* = {T_STAR:.4f} s, DT_BOX = {DT_BOX:.6f} s, HALF_T = {HALF_T:.4f} s, "
         f"W/DT = {W_AVG / DT_BOX:.4f} samples")
    out = {"E": {}, "imag": _max_imag(st.temperature), "H": [], "H_fields": [], "ndof": ndof}

    def _go(label, schedule):
        res = _cn_p(msh, ft, schedule, DEGREE)
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

    # (Δ) and (Δ-ctl) on fields at t*, on the P2 space.
    sig = fem.Function(V2)
    sig.x.array[:] = sum((-1.0) ** m * (f.x.array - T_INF) for m, f in enumerate(out["H_fields"]))
    sig9 = fem.Function(V2)
    sig9.x.array[:] = sig.x.array - (-1.0) ** 9 * (out["H_fields"][9].x.array - T_INF)
    th_duty = out["duty_field"] - T_INF
    out["delta"] = _rel_p(msh, th_duty, sig)
    out["delta_ctl"] = _rel_p(msh, th_duty, sig9)
    out["dprime"] = _rel_p(msh, out["one_field"] - T_INF, out["default_field"] - T_INF)
    out["t"] = np.arange(N_MAIN + 1) * DT_BOX
    out["ids"] = _identities(out["t"], out["duty"], out["H"])
    clk.lap("(Δ), (D′) readings and identities")
    _say(f"MAT21-4 [runs] max|Im theta| = {out['imag']:.3e} (predicted 0); elapsed "
         f"{time.perf_counter() - t0:.1f} s")
    return out


def test_precondition_schedule_and_degree():
    assert (N_STEPS, SPH_PIN, N_SWITCH_PIN, DUTY_PIN) == (N_MAIN, STEPS_PER_HALF, N_SWITCH, DUTY)
    assert T_STAR == 0.2 * T_C
    assert T_STAR == T_STAR_1
    assert T_INF == T_INF_0
    assert HALF_T == STEPS_PER_HALF * DT_BOX
    assert DEGREE == 2 and CMP_DEGREE == 3


@complex_only
def test_precondition_degree_dofs_identities(runs_p):
    assert runs_p["ndof"] == _dofs_pred(M_ID), f"(Pdeg) {runs_p['ndof']} != {_dofs_pred(M_ID)}"


@complex_only
def test_s0_p2_schedule_census(runs_p):
    s = "".join(str(int(x)) for x in runs_p["mult"])
    _say(f"MAT21-4 (S0) schedule (2 half-steps, then 39 steps): {s} "
         f"(must equal step 2's 11111000011110000111100001111000011110000)")
    assert runs_p["census_ok"]


@complex_only
def test_delta_p2_duhamel_fields(runs_p):
    _say(f"MAT21-4 (Δ) whole-box rel L2 |theta_duty - sum (-1)^m theta_Hm| at t* = "
         f"{runs_p['delta']:.3e} (bound L_BAND {L_BAND:.0e}, predicted <= 3e-13; P1 7.748e-14)")
    assert runs_p["delta"] <= L_BAND


@complex_only
def test_d_prime_p2_duty_one_is_default(runs_p):
    _say(f"MAT21-4 (D') chi = 1 vs default at t*: rel L2 = {runs_p['dprime']:.3e} "
         f"(bound D_BAND {D_BAND:.0e}, predicted <= 1e-14; P1 3.105e-15)")
    assert runs_p["dprime"] <= D_BAND


@complex_only
def test_h0_p2_series_duhamel(runs_p):
    ids = runs_p["ids"]
    h0 = ids["h0"]
    t, duty = runs_p["t"], runs_p["duty"]
    cf = np.array([_box_mean_duty(tn) for tn in t])
    nz = np.abs(cf) > 0
    dev = np.max(np.abs(duty[nz] - cf[nz]) / np.abs(cf[nz]))
    _say("MAT21-4 record <theta_duty>_n FEM vs closed-form mean: "
         + "; ".join(f"n{i} {duty[i]:.6f}/{cf[i]:.6f}" for i in range(len(t))))
    _say(f"MAT21-4 record max_n |FEM mean - closed-form mean|/closed-form = {dev:.3e} "
         f"(predicted 3e-5 to 1.5e-3; P1 1.222e-3; printed, not asserted)")
    r = runs_p["duty"] - sum((-1.0) ** m * h for m, h in enumerate(runs_p["H"]))
    _say("MAT21-4 (H0) per-step r series (K): " + " ".join(f"{x:.2e}" for x in np.abs(r))
         + f"; r_max {ids['r_max']:.4e} K")
    _say(f"MAT21-4 (H0) max_n |duty - sum (-1)^m H_m| / max|duty| = {h0:.3e} "
         f"(bound L_BAND {L_BAND:.0e}, predicted 1e-13 to 1e-12, central ≈ 4e-13; P1 3.801e-13)")
    assert h0 <= L_BAND


@complex_only
def test_fem_p2_six_minute_operator(runs_p):
    ids = runs_p["ids"]
    t, duty = runs_p["t"], runs_p["duty"]
    nw = int(ids["ok"].sum())
    ceil_ = ids["r_max"] / ids["scale"]
    rhs = ids["r_max"] * (1.0 + 1e-12) + D_BAND * ids["scale"]
    rhs2 = ids["r_max"] * (1.0 + 1e-12) + D_BAND * ids["scale2"]
    _say(f"MAT21-4 (FEM) as posed = (FEM-ctl), printed only: W = {W_AVG} s, max |a - b| / max|a| "
         f"over {nw} windows = {ids['fem']:.3e}; 2-D twin {ids['fem2d']:.3e}; ceiling r_max/max|a| = "
         f"{ceil_:.3e} (P1 1.460e-13 / 1.456e-13 / 3.849e-13)")
    _say("MAT21-4 (FEM') per-window residuals: " + " ".join(f"{x:.2e}" for x in ids["res_p"]))
    _say(f"MAT21-4 (FEM') max |a - b - rr| / max|a| = {ids['fem_p']:.3e}; (FEM'-2d) "
         f"{ids['fem2d_p']:.3e} (bound D_BAND {D_BAND:.0e}, predicted <= 1e-14, worst-case ≈ 3e-14; "
         f"P1 1.564e-15 / 1.562e-15)")
    _say(f"MAT21-4 (FEM-c) max |a - b| = {ids['fem_c_lhs']:.4e} K <= r_max(1 + 1e-12) + D_BAND max|a| "
         f"= {rhs:.4e} K (r_max {ids['r_max']:.4e} K); 2-D {ids['fem_c_lhs2']:.4e} K <= {rhs2:.4e} K "
         f"(P1 9.0039e-14 <= 2.9903e-13 K)")
    _say(f"MAT21-4 record FEM six-minute maximum {ids['a_mx']:.5f} K at t_max = {ids['a_tmx']:.4f} s "
         f"(predicted 0.6167 ± 0.0005 K at t* = {T_STAR:.4f} s; P1 0.61668 K)")
    d_w = cem43(t, T_INF + np.real(duty) / QV)
    _say(f"MAT21-4 record FEM cem43 of T_inf + <theta>/(qV) = {d_w:.4e} s per W of qV "
         f"(predicted 4.049e-11 ± 0.1 %; P1 4.0493e-11)")
    assert ids["fem_p"] <= D_BAND, f"(FEM') {ids['fem_p']:.3e}"
    assert ids["fem2d_p"] <= D_BAND, f"(FEM'-2d) {ids['fem2d_p']:.3e}"
    assert ids["fem_c_lhs"] <= rhs, f"(FEM-c) {ids['fem_c_lhs']:.4e} > {rhs:.4e} K"
    assert ids["fem_c_lhs2"] <= rhs2, f"(FEM-c 2-D) {ids['fem_c_lhs2']:.4e} > {rhs2:.4e} K"


@complex_only
def test_e_p2_every_run_balances(runs_p):
    worst = max(runs_p["E"], key=runs_p["E"].get)
    _say("MAT21-4 (E) per run " + ", ".join(f"{k_}: {v:.3e}" for k_, v in runs_p["E"].items())
         + f"; worst {runs_p['E'][worst]:.3e} on {worst} (bound {E_BAND:.0e}, predicted <= 1e-11)")
    for name, r in runs_p["E"].items():
        assert r <= E_BAND, f"(E) {name}: {r:.3e}"


@complex_only
def test_control_delta_p2_m9_dropped_misses(runs_p):
    _say(f"MAT21-4 control (Δ-ctl) m = 9 dropped at t*: miss = {runs_p['delta_ctl']:.4e} "
         f"(predicted 0.177-0.180; P1 1.7843e-1; asserted > L_BAND)")
    assert runs_p["delta_ctl"] > L_BAND


@complex_only
def test_control_wc_p2_shifted_window_misses(runs_p):
    ids = runs_p["ids"]
    _say(f"MAT21-4 control (Wc-FEM) shift P/2 = {STEPS_PER_HALF} samples: miss {ids['wc']:.4e} "
         f"over {ids['idx'].size} windows (predicted 0.139-0.143; P1 1.4084e-1; asserted > D_BAND)")
    assert ids["wc"] > D_BAND


@pytest.fixture(scope="module")
def duty_cf_p():
    t0 = time.perf_counter()
    clk = _Clock4("cf")
    msh, ft, ns, n, v_h, _st = _rung_p(M_T, clk, DEGREE)
    res = _cn_p(msh, ft, _chi_box, DEGREE)
    clk.lap(f"CN N={N_MAIN} duty degree={DEGREE}")
    res.snapshots.clear()
    out = {"ndof": int(res.temperature.function_space.dofmap.index_map.size_global)}
    th = res.temperature - T_INF
    chi = _interior_chi(msh)
    ex, npts = _exact_p(msh, lambda x, y, z: _duty_cf(x, y, z, T_STAR))
    clk.lap(f"ten-run series interpolation at t* ({npts} P3 points)")
    out["T"] = _rel_p(msh, th, ex, chi)
    out["whole"] = _rel_p(msh, th, ex)
    ex_s, _ = _exact_p(msh, _trans_ex(T_STAR))
    clk.lap("step-response interpolation (c3)")
    out["c3"] = _rel_p(msh, th, ex_s, chi)
    try:
        vals, ok = evaluate_vector_field_parallel(res.temperature, np.zeros((1, 3)))
        th_c = float(np.real(np.ravel(vals[0])[0])) - T_INF if bool(ok[0]) else float("nan")
    except Exception as exc:  # printed quantity only
        _say(f"MAT21-4 centre evaluation failed: {exc!r}")
        th_c = float("nan")
    th_c_ex = float(np.real(_duty_cf(np.zeros(1), np.zeros(1), np.zeros(1), T_STAR)[0]))
    _say(f"\nMAT21-4 [cf] rung m = {M_T} {tuple(ns)}: cells {n}; P2 dofs {out['ndof']} "
         f"(predicted {_dofs_pred(M_T)}); P3 points {npts} (predicted {_dofs_pred(M_T, CMP_DEGREE)})")
    _say(f"MAT21-4 (T-duty) interior rel L2 at t* = {out['T']:.6e} (bound {INTERIOR_BAND:.0e}, "
         f"predicted 1.4e-4, range 7e-5 to 6e-4; step 3 P2 (T) 1.3229e-4; step 2 P1 m = 6 "
         f"9.401399e-4); whole {out['whole']:.6e} (printed, predicted 1.5e-4 to 8e-4; P1 1.414345e-3); "
         f"theta_h(0) = {th_c:.6f} K vs ten-run series {th_c_ex:.6f} K "
         f"(rel {th_c / th_c_ex - 1.0:+.3e}, predicted |rel| <= 3e-4); "
         f"balance {res.max_balance_residual:.3e}; max|Im| {_max_imag(res.temperature):.3e}")
    _say(f"MAT21-4 control (c3) unscheduled step response as comparand: interior miss = "
         f"{out['c3']:.4e} (predicted 0.49-0.52; P1 m = 6 5.0282e-1; asserted > {INTERIOR_BAND:.0e})")
    clk.lap("P2 readings")

    # (c-deg-duty): the same duty schedule at degree 1 on the same mesh, the same P3
    # comparand and χ; no second interpolation.
    res1 = _cn_p(msh, ft, _chi_box, DEG_CTL)
    clk.lap(f"CN N={N_MAIN} duty degree={DEG_CTL}")
    res1.snapshots.clear()
    out["ndof1"] = int(res1.temperature.function_space.dofmap.index_map.size_global)
    th1 = res1.temperature - T_INF
    out["e1"] = _rel_p(msh, th1, ex, chi)
    e1_whole = _rel_p(msh, th1, ex)
    c3_p1 = _rel_p(msh, th1, ex_s, chi)
    _say(f"MAT21-4 control (c-deg-duty) m = {M_T}: P1 dofs {out['ndof1']} (predicted "
         f"{_dofs_pred(M_T, DEG_CTL)}); interior P1 e1 = {out['e1']:.6e} vs P2 e2 = {out['T']:.6e}, "
         f"ratio e1/e2 = {out['e1'] / out['T']:.4f} (predicted e1 ≈ 2.1e-3, range 1.2e-3 to 4e-3; "
         f"ratio ≈ 15, range 3-50; asserted e1 > e2, label predicted); whole P1 {e1_whole:.6e}; "
         f"balance {res1.max_balance_residual:.3e}")
    _say(f"MAT21-4 record P1 control (c3) on m = {M_T} (printed only, rung-vs-degree separation): "
         f"interior miss = {c3_p1:.4e} (P2 m = {M_T} {out['c3']:.4e}; P1 m = 6 5.0282e-1)")
    clk.lap("P1 control readings")
    _say(f"MAT21-4 [cf] elapsed {time.perf_counter() - t0:.1f} s")
    return out


@complex_only
def test_precondition_degree_dofs_t_duty(duty_cf_p):
    assert duty_cf_p["ndof"] == _dofs_pred(M_T), f"(Pdeg) P2 {duty_cf_p['ndof']}"
    assert duty_cf_p["ndof1"] == _dofs_pred(M_T, DEG_CTL), f"(Pdeg) P1 {duty_cf_p['ndof1']}"


@complex_only
def test_t_duty_p2_interior_against_closed_form(duty_cf_p):
    assert duty_cf_p["T"] <= INTERIOR_BAND, f"(T-duty) {duty_cf_p['T']:.4e} > 1 %"


@complex_only
def test_control_c3_p2_step_response_misses(duty_cf_p):
    assert duty_cf_p["c3"] > INTERIOR_BAND, f"(c3) {duty_cf_p['c3']:.4e} inside 1 %"


@complex_only
def test_control_cdeg_p1_duty_worse_on_same_rung(duty_cf_p):
    assert duty_cf_p["e1"] > duty_cf_p["T"], (
        f"(c-deg-duty) P1 {duty_cf_p['e1']:.4e} <= P2 {duty_cf_p['T']:.4e}")
