"""`MAT-21` step 13, W2 — `MAT-20`'s dose and six-minute window on the gel box's
duty-cycled FEM point histories (P38's schedule, P36's transverse quadratic
source; P2, m = 4, CN N = 40 to t*) against the continuous run (§9 item 567,
§10 chain P41; brief `docs/planning/chunks/MAT-21.md` "§9 item 567", review
ruling 2026-10-20 03:00).

Sets no band and moves none; every band is imported by name.  Three CN runs
(default, χ ≡ 1, duty), each read on item 565's one batch ``HIST_PTS`` on all
41 snapshots.  Asserted: preconditions ((Pdeg) 116 025 on four solves, (S0)
census + schedule string, validity, (E) balance ≤ ``E_BAND``, max T < 43);
(R39′) the default run's dose strings equal ``P39P_FEM_DOSE``; (C100-h)
|D(χ ≡ 1)/D(default) − 1| ≤ ``D_BAND`` (the FEM side of "bitwise", review
ruling 2026-10-20 (2)); (Q-hd) / (Win-hd) on the duty history at peak and
equal mean power at ``DOSE_BAR`` / ``AVG_BAR``; (Rc′-hd) the R swap misses.
"""

from __future__ import annotations

import math
import time

import numpy as np
import pytest
import ufl
from mpi4py import MPI

from fem_em_solver.thermal import solve_transient_conduction
from fem_em_solver.thermal.dose import cem43, running_average
from tests.complex_mode import complex_only
from tests.validation.test_mat10_sar_to_heat_sphere import ROBIN_TAG, _say
from tests.validation.test_mat12_transient_sar_heating_sphere import E_BAND, RHO_C
from tests.validation.test_mat13_transient_perfusion_sphere import INTERIOR_BAND
from tests.validation.test_mat15_coil_driven_heating_phantom import H_GEL, K_GEL, N_MAIN, T_INF
from tests.validation.test_mat15_finite_cylinder_fem import _max_imag
from tests.validation.test_mat17_duty_cycle_coil_driven import D_BAND, DUTY
from tests.validation.test_mat17_duty_cycle_fem import _census
from tests.validation.test_mat21_box_closed_form import B_BOX, C_BOX
from tests.validation.test_mat21_box_closed_form import T_INF as T_INF_0
from tests.validation.test_mat21_box_closed_form import Q
from tests.validation.test_mat21_box_dose_closed_form import (
    AVG_BAR,
    DOSE_BAR,
    DOSE_PTS,
    HIST_LABELS,
    HIST_PTS,
    T_GRID,
    T_STAR,
    W_B2,
    _q_check,
    _series_hist,
    _win_check,
)
from tests.validation.test_mat21_box_duty_closed_form import DT_BOX
from tests.validation.test_mat21_box_duty_dose_closed_form import P39P_SERIES_DOSE, _duty_hist
from tests.validation.test_mat21_box_duty_fem import _chi_box
from tests.validation.test_mat21_box_fem_p2 import DEGREE, _dofs_pred, _rung_p
from tests.validation.test_mat21_box_probe_points_fem import _probe_vals
from tests.validation.test_mat21_box_quadratic_duty_fem import SCHEDULE
from tests.validation.test_mat21_box_separable_fem import M4

LN4 = math.log(4.0)
# Item 565's W2 FEM doses at :.10e, DOSE_PTS order, copied from
# docs/testing/logs/20261020T140522Z_MAT-21-step12-W2-dose-fem.log:85-87.
P39P_FEM_DOSE = ("4.7060617108e-10", "4.7060617108e-10", "4.1364655592e-11")
PRED_C100H = 5e-14


class _Clock13:
    """Per-phase wall clock, printed on rank 0 with flush, prefix MAT21-13."""

    def __init__(self, tag):
        self.tag = tag
        self.t = time.perf_counter()

    def lap(self, phase):
        now = time.perf_counter()
        _say(f"MAT21-13 [{self.tag}] clock {phase}: {now - self.t:.1f} s")
        self.t = now


def _cn_snap(msh, ft, q_quad, schedule):
    kw = {} if schedule is None else {"source_schedule": schedule}
    return solve_transient_conduction(
        msh,
        ft,
        k=K_GEL,
        rho_c=RHO_C,
        q=q_quad,
        h=H_GEL,
        t_inf=T_INF,
        t_init=T_INF,
        robin_tag=ROBIN_TAG,
        dt=T_STAR / N_MAIN,
        n_steps=N_MAIN,
        theta=0.5,
        startup_half_steps=2,
        degree=DEGREE,
        snapshot_steps=tuple(range(N_MAIN + 1)),
        **kw,
    )


@pytest.fixture(scope="module")
def ddose_m4():
    comm = MPI.COMM_WORLD
    t0 = time.perf_counter()
    clk = _Clock13("ddose")
    out = {"imag": 0.0, "E": {}, "ndof": {}, "th": {}, "ok": {}}
    _say(
        f"\nMAT21-13 DEGREE = {DEGREE}; m = {M4}; N_MAIN = {N_MAIN}; ranks {comm.size}; "
        f"HIST_PTS {HIST_PTS}"
    )

    msh, ft, ns, n, v_h, st = _rung_p(M4, clk, DEGREE)
    out["ndof"]["uniform"] = int(st.temperature.function_space.dofmap.index_map.size_global)

    xs = ufl.SpatialCoordinate(msh)
    q_quad = Q * 3.0 * (xs[1] * xs[1] + xs[2] * xs[2]) / (B_BOX**2 + C_BOX**2)
    pts = np.asarray(HIST_PTS, dtype=float)

    for label, schedule in (("default", None), ("chi1", lambda t: 1.0), ("duty", _chi_box)):
        r = _cn_snap(msh, ft, q_quad, schedule)
        clk.lap(f"CN N={N_MAIN} {label} ({N_MAIN + 1} snapshots)")
        out["ndof"][label] = int(r.temperature.function_space.dofmap.index_map.size_global)
        out["imag"] = max(out["imag"], _max_imag(r.temperature))
        out["E"][label] = float(r.max_balance_residual)
        if label == "duty":
            out["sched"] = "".join(str(int(x)) for x in r.source_multipliers)
            _census(r, DT_BOX, 2, DUTY * T_STAR, "quad duty")
            out["census_ok"] = True
        th_h = np.empty((N_MAIN + 1, len(HIST_PTS)))
        ok_all = np.ones(len(HIST_PTS), dtype=bool)
        for k in range(N_MAIN + 1):
            th, ok = _probe_vals(r.snapshots[k], pts)
            th_h[k] = th
            ok_all &= ok
        r.snapshots.clear()
        del r
        out["th"][label] = th_h
        out["ok"][label] = ok_all
        clk.lap(f"{label}: {N_MAIN + 1} owner-route batches")

    out["th_s"] = _series_hist(HIST_PTS)
    out["th_sd"] = _duty_hist(HIST_PTS)
    clk.lap("series histories (continuous + duty)")
    _say(f"MAT21-13 fixture elapsed {time.perf_counter() - t0:.1f} s; max|Im T| {out['imag']:.3e}")
    return out


def _T(th):
    return T_INF + np.asarray(th, dtype=float)


@complex_only
def test_preconditions(ddose_m4):
    pred = _dofs_pred(M4)
    _say(f"MAT21-13 (Pdeg) P2 dofs {ddose_m4['ndof']} (predicted {pred} = 116025)")
    assert pred == 116025
    assert len(ddose_m4["ndof"]) == 4
    for key, v in ddose_m4["ndof"].items():
        assert v == pred, f"(Pdeg) {key}: {v}"
    _say(f"MAT21-13 (S0) schedule {ddose_m4['sched']} (must equal {SCHEDULE})")
    assert ddose_m4["census_ok"] and ddose_m4["sched"] == SCHEDULE
    for key, ok in ddose_m4["ok"].items():
        assert bool(np.all(ok)), (key, ok)
    _say(
        "MAT21-13 (E) max balance residual: "
        + ", ".join(f"{k_} {v:.3e}" for k_, v in ddose_m4["E"].items())
        + f" (bound {E_BAND:.0e})"
    )
    assert len(ddose_m4["E"]) == 3
    for key, v in ddose_m4["E"].items():
        assert v <= E_BAND, (key, v)
    assert T_INF == T_INF_0
    for key, arr in ddose_m4["th"].items():
        for j in range(len(HIST_PTS)):
            assert float(np.max(_T(arr[:, j]))) < 43.0, (key, j)
            assert float(np.max(_T(arr[:, j] / DUTY))) < 43.0, (key, j)


@complex_only
def test_r39p_default_reproduces_p39p_fem_dose(ddose_m4):
    th = ddose_m4["th"]["default"]
    got = tuple(f"{cem43(T_GRID, _T(th[:, j])):.10e}" for j in range(len(DOSE_PTS)))
    _say(f"MAT21-13 (R39') default FEM doses {got} (P39' {P39P_FEM_DOSE})")
    assert got == P39P_FEM_DOSE, got


@complex_only
def test_c100_h_chi_one_dose_at_d_band(ddose_m4):
    th_d0, th_1 = ddose_m4["th"]["default"], ddose_m4["th"]["chi1"]
    rows = []
    for j in range(len(HIST_PTS)):
        d0 = cem43(T_GRID, _T(th_d0[:, j]))
        d1 = cem43(T_GRID, _T(th_1[:, j]))
        rel = abs(d1 / d0 - 1.0)
        dth = float(np.max(np.abs(th_1[:, j] - th_d0[:, j])))
        tag = "ASSERTED" if j < 3 else "printed"
        _say(
            f"MAT21-13 (C100-h) {HIST_LABELS[j]} [{tag}]: D(chi=1) {d1:.15e} s, D(default) "
            f"{d0:.15e} s, |ratio - 1| {rel:.3e} (D_BAND {D_BAND:.0e}; predicted <= "
            f"{PRED_C100H:.0e}); max_n|dtheta| {dth:.3e} K, ln4*that {LN4 * dth:.3e}"
        )
        rows.append(rel)
    for j in range(len(DOSE_PTS)):
        assert rows[j] <= D_BAND, ("(C100-h)", HIST_LABELS[j], rows[j])


@complex_only
def test_q_hd_win_hd_and_rc_hd(ddose_m4):
    th_d = ddose_m4["th"]["duty"]
    for scale, lab in ((1.0, "peak"), (1.0 / DUTY, "equal-mean")):
        for j, pt in enumerate(DOSE_PTS):
            th = th_d[:, j] * scale
            dose, ref, rel, swap, miss = _q_check(_T(th))
            mx, t_mx, lo, hi, wref, wrel = _win_check(th)
            _say(
                f"MAT21-13 (Q-hd) {lab} {HIST_LABELS[j]}: FEM duty CEM43 {dose:.15e} s vs quad "
                f"{ref:.15e} s rel {rel:.3e} (bar {DOSE_BAR:g}, ASSERTED); (Rc'-hd) swap miss "
                f"{miss:.4e} (asserted > {DOSE_BAR:g})"
            )
            _say(
                f"MAT21-13 (Win-hd) {lab} {HIST_LABELS[j]}: max six-minute mean {mx:.15e} K at "
                f"{t_mx:.4f} s; (ii) [{lo:.12e}, {hi:.12e}]; (iii) quad {wref:.15e} rel {wrel:.3e} "
                f"(bar {AVG_BAR:g}, ASSERTED)"
            )
            assert rel <= DOSE_BAR, (lab, pt, rel)
            assert miss > DOSE_BAR, (lab, pt, miss)
            assert T_GRID[0] + W_B2 <= t_mx <= T_GRID[-1], (lab, t_mx)
            assert lo <= mx <= hi, (lab, lo, mx, hi)
            assert wrel <= AVG_BAR, (lab, pt, wrel)


@complex_only
def test_printed_records(ddose_m4):
    """Printed records (no band); asserts only that they were formed."""
    th_c, th_d = ddose_m4["th"]["default"], ddose_m4["th"]["duty"]
    th_s, th_sd = ddose_m4["th_s"], ddose_m4["th_sd"]
    W = float(W_B2)
    for j in range(len(HIST_PTS)):
        dc_h = cem43(T_GRID, _T(th_c[:, j]))
        dd_h = cem43(T_GRID, _T(th_d[:, j]))
        de_h = cem43(T_GRID, _T(th_d[:, j] / DUTY))
        dc_s = cem43(T_GRID, _T(th_s[:, j]))
        dd_s = cem43(T_GRID, _T(th_sd[:, j]))
        de_s = cem43(T_GRID, _T(th_sd[:, j] / DUTY))
        M = float(np.max(np.abs(th_d[:, j] - th_sd[:, j])))
        cap = INTERIOR_BAND * float(np.max(th_sd[:, j]))
        _a, mx_h, t_h = running_average(T_GRID, th_d[:, j], window=W)
        _a, mx_s, t_s = running_average(T_GRID, th_sd[:, j], window=W)
        _say(
            f"MAT21-13 DUTY DOSE {HIST_LABELS[j]}: FEM ratio equal-mean {de_h / dc_h:.6f} / peak "
            f"{dd_h / dc_h:.6f} (W1 closed form {de_s / dc_s:.6f} / {dd_s / dc_s:.6f}); theta ratio "
            f"at t* FEM {th_d[-1, j] / th_c[-1, j]:.5f}, series {th_sd[-1, j] / th_s[-1, j]:.5f}; "
            f"FEM duty dose {dd_h:.10e} s vs series {dd_s:.15e} s: |D_h/D_s - 1| "
            f"{abs(dd_h / dd_s - 1.0):.4e} beside ln4*max|dtheta| {LN4 * M:.4e} (M {M:.4e} K) and "
            f"INTERIOR_BAND*max theta_sd {cap:.4e} K ({cap / M:.1f}x) [printed, not P41's target]; "
            f"six-minute peak FEM {mx_h:.9f} K at {t_h:.2f} s vs series {mx_s:.9f} K at {t_s:.2f} s"
        )
        assert np.isfinite(de_h / dc_h)
    _say(f"MAT21-13 P39' series dose strings (W1 import) {P39P_SERIES_DOSE}")
