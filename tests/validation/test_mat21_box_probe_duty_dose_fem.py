"""`MAT-21` step 16, W2 — `MAT-20`'s dose and six-minute window and P39′'s (Bd)
bound on the gel box's FEM point histories at all six F2182 probe points under
P41's D = 0.5 and P42's D = 0.25 schedules (P36's source; P2, m = 4, CN N = 40
to t*) (§9 item 577, §10 chain P44; brief `docs/planning/chunks/MAT-21.md`
"§9 item 577", review ruling 2026-10-20 18:00).

Sets no band and moves none; every band is imported by name.  Three CN runs
(default, d50, d25), each read on two fixed batches per snapshot — ``HIST_PTS``
then ``NEW_PTS`` — never one merged batch.  Asserted: preconditions; (R39′),
(R43), (R41-h), (R42-h) string reproductions; (Q-hP) / (Win-hP) at P3–P5;
(Bd50) / (Bd25) (i) and (ii) at all six points (first-order form, ruling
2026-10-20 03:00 (1)), every row printed before any is asserted; (Rc′-hP) at
P0–P2.  Controls (rule (e)): (Rc′-hP) at P3–P5 printed only; (c-50-dose) /
(c-25-dose) asserted only > ``INTERIOR_BAND``.  Ruling 18:00 (4): the
``PREDICTED RED: (Bd25)(ii) at P5, M_p / (INTERIOR_BAND·max θ_s) range
[1.0, 1.4]``.
"""

from __future__ import annotations

import math
import time

import numpy as np
import pytest
import ufl
from mpi4py import MPI

from fem_em_solver.thermal.dose import cem43, running_average
from tests.complex_mode import complex_only
from tests.validation.test_mat10_sar_to_heat_sphere import _say
from tests.validation.test_mat12_transient_sar_heating_sphere import E_BAND
from tests.validation.test_mat13_transient_perfusion_sphere import INTERIOR_BAND
from tests.validation.test_mat15_coil_driven_heating_phantom import N_MAIN, T_INF
from tests.validation.test_mat15_finite_cylinder_fem import _max_imag
from tests.validation.test_mat17_duty_cycle_coil_driven import DUTY
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
from tests.validation.test_mat21_box_dose_fem import R36, _bd
from tests.validation.test_mat21_box_duty25_dose_closed_form import (
    DUTY25,
    SCHEDULE_25,
    T_SW25,
    _chi25,
)
from tests.validation.test_mat21_box_duty25_dose_fem import P41_FEM_EQ
from tests.validation.test_mat21_box_duty_closed_form import DT_BOX, T_SW
from tests.validation.test_mat21_box_duty_dose_closed_form import _duty_hist
from tests.validation.test_mat21_box_duty_dose_fem import LN4, P39P_FEM_DOSE, _cn_snap
from tests.validation.test_mat21_box_duty_fem import _chi_box
from tests.validation.test_mat21_box_fem_p2 import DEGREE, _dofs_pred, _rung_p
from tests.validation.test_mat21_box_probe_dose_closed_form import NEW_LABELS, NEW_PTS
from tests.validation.test_mat21_box_probe_points_fem import PROBE_PTS, _probe_vals
from tests.validation.test_mat21_box_quadratic_duty_fem import SCHEDULE
from tests.validation.test_mat21_box_separable_fem import M4

# P42's W2 FEM D = 0.25 equal-mean ratios at :.6f, HIST_PTS order, copied from
# docs/testing/logs/20261020T185616Z_MAT-21-step14-W2-duty25-dose-fem.log:104-107.
P42_FEM_EQ = ("1.275667", "1.275667", "1.006985", "1.316409")
# P43's W2 FEM doses at :.10e, NEW_PTS order, copied from
# docs/testing/logs/20261020T190411Z_MAT-21-step15-W2-probe-dose-fem.log:109-111.
P43_FEM_DOSE = ("4.1094392184e-11", "4.1094392184e-11", "4.3189685419e-11")
RUNS = (("default", None, 1.0), ("d50", _chi_box, DUTY), ("d25", _chi25, DUTY25))
SIX = tuple(("h", j, HIST_LABELS[j]) for j in range(len(DOSE_PTS))) + tuple(
    ("n", j, NEW_LABELS[j]) for j in range(len(NEW_PTS))
)
# Predictions written first (brief), printed only.
PRED_M = {
    "d50": {"x": (3.5e-5, 9e-5), "z": (1.0e-4, 2.5e-4)},
    "d25": {"x": (1.8e-5, 9e-5), "z": (5e-5, 2.5e-4)},
}
PRED_HEAD = {"d50": {"x": (3.0, 7.5), "z": (2.0, 5.5)}, "d25": {"x": (1.5, 7.5), "z": (1.0, 5.5)}}
PRED_RATIO = (0.01, 1.0)
PRED_RC = (1e6, 1e7)
PRED_C50 = (0.5, 1.3)
PRED_C25 = (0.25, 0.8)
PRED_ASYM = 1e-12
# Ruling 2026-10-20 18:00 (4):
# PREDICTED RED: (Bd25)(ii) at P5, M_p / (INTERIOR_BAND·max θ_s) range [1.0, 1.4]
PRED_RED_BD25_P5 = (1.0, 1.4)


class _Clock16:
    """Per-phase wall clock, printed on rank 0 with flush, prefix MAT21-16."""

    def __init__(self, tag):
        self.tag = tag
        self.t = time.perf_counter()

    def lap(self, phase):
        now = time.perf_counter()
        _say(f"MAT21-16 [{self.tag}] clock {phase}: {now - self.t:.1f} s")
        self.t = now


@pytest.fixture(scope="module")
def pdd_m4():
    comm = MPI.COMM_WORLD
    t0 = time.perf_counter()
    clk = _Clock16("pdd")
    out = {"imag": 0.0, "E": {}, "ndof": {}, "th": {}, "ok": {}, "sched": {}, "census": {}}
    _say(
        f"\nMAT21-16 DEGREE = {DEGREE}; m = {M4}; N_MAIN = {N_MAIN}; ranks {comm.size}; "
        f"HIST_PTS {HIST_PTS}; NEW_PTS {NEW_PTS} (two fixed batches)"
    )

    msh, ft, ns, n, v_h, st = _rung_p(M4, clk, DEGREE)
    out["ndof"]["uniform"] = int(st.temperature.function_space.dofmap.index_map.size_global)

    xs = ufl.SpatialCoordinate(msh)
    q_quad = Q * 3.0 * (xs[1] * xs[1] + xs[2] * xs[2]) / (B_BOX**2 + C_BOX**2)
    p_h = np.asarray(HIST_PTS, dtype=float)
    p_n = np.asarray(NEW_PTS, dtype=float)

    for label, schedule, D in RUNS:
        r = _cn_snap(msh, ft, q_quad, schedule)
        clk.lap(f"CN N={N_MAIN} {label} ({N_MAIN + 1} snapshots)")
        out["ndof"][label] = int(r.temperature.function_space.dofmap.index_map.size_global)
        out["imag"] = max(out["imag"], _max_imag(r.temperature))
        out["E"][label] = float(r.max_balance_residual)
        if schedule is not None:
            _census(r, DT_BOX, 2, D * T_STAR, f"quad {label}")
            out["census"][label] = True
            out["sched"][label] = "".join(str(int(x)) for x in r.source_multipliers)
        th_h = np.empty((N_MAIN + 1, len(HIST_PTS)))
        th_n = np.empty((N_MAIN + 1, len(NEW_PTS)))
        ok_h = np.ones(len(HIST_PTS), dtype=bool)
        ok_n = np.ones(len(NEW_PTS), dtype=bool)
        for k in range(N_MAIN + 1):
            a, oa = _probe_vals(r.snapshots[k], p_h)
            b, ob = _probe_vals(r.snapshots[k], p_n)
            th_h[k], th_n[k] = a, b
            ok_h &= oa
            ok_n &= ob
        if label == "default":
            fin_h, okf_h = _probe_vals(r.temperature, p_h)
            fin_n, okf_n = _probe_vals(r.temperature, p_n)
            ok_h &= okf_h
            ok_n &= okf_n
            out["fin_h"], out["fin_n"] = fin_h, fin_n
        r.snapshots.clear()
        del r
        out["th"][label] = {"h": th_h, "n": th_n}
        out["ok"][label] = (ok_h, ok_n)
        clk.lap(f"{label}: 2 x {N_MAIN + 1} owner-route batches")

    out["s"] = {
        "default": {"h": _series_hist(HIST_PTS), "n": _series_hist(NEW_PTS)},
        "d50": {
            "h": _duty_hist(HIST_PTS, switch_times=T_SW),
            "n": _duty_hist(NEW_PTS, switch_times=T_SW),
        },
        "d25": {
            "h": _duty_hist(HIST_PTS, switch_times=T_SW25),
            "n": _duty_hist(NEW_PTS, switch_times=T_SW25),
        },
    }
    clk.lap("series histories (continuous + D = 0.5 + D = 0.25, two batches)")
    _say(f"MAT21-16 fixture elapsed {time.perf_counter() - t0:.1f} s; max|Im T| {out['imag']:.3e}")
    return out


def _T(th):
    return T_INF + np.asarray(th, dtype=float)


@complex_only
def test_preconditions(pdd_m4):
    pred = _dofs_pred(M4)
    _say(f"MAT21-16 (Pdeg) P2 dofs {pdd_m4['ndof']} (predicted {pred} = 116025)")
    assert pred == 116025
    assert len(pdd_m4["ndof"]) == 4
    for key, v in pdd_m4["ndof"].items():
        assert v == pred, f"(Pdeg) {key}: {v}"
    assert tuple(NEW_PTS) == tuple(PROBE_PTS[3:6])
    _say(
        f"MAT21-16 (S0) schedules {pdd_m4['sched']} (must equal d50 {SCHEDULE}, d25 {SCHEDULE_25})"
    )
    assert pdd_m4["census"] == {"d50": True, "d25": True}
    assert pdd_m4["sched"]["d50"] == SCHEDULE and pdd_m4["sched"]["d25"] == SCHEDULE_25
    for key, (ok_h, ok_n) in pdd_m4["ok"].items():
        _say(f"MAT21-16 validity {key}: HIST_PTS {ok_h.tolist()} NEW_PTS {ok_n.tolist()}")
        assert bool(np.all(ok_h)) and bool(np.all(ok_n)), key
    _say(
        "MAT21-16 (E) max balance residual: "
        + ", ".join(f"{k_} {v:.3e}" for k_, v in pdd_m4["E"].items())
        + f" (bound {E_BAND:.0e}; predicted 2.3e-13)"
    )
    assert len(pdd_m4["E"]) == 3
    for key, v in pdd_m4["E"].items():
        assert v <= E_BAND, (key, v)
    assert T_INF == T_INF_0
    th_d = pdd_m4["th"]["default"]
    r36 = (f"{th_d['h'][40, 2]:.6f}", f"{th_d['h'][40, 3]:.6f}")
    same_h = bool(np.array_equal(th_d["h"][40], pdd_m4["fin_h"]))
    same_n = bool(np.array_equal(th_d["n"][40], pdd_m4["fin_n"]))
    _say(
        f"MAT21-16 (R36) theta_h(t*) centre / (0,0.9B,0) {r36} (P36 {R36}); snapshot n=40 "
        f"array_equal final: HIST_PTS {same_h}, NEW_PTS {same_n}"
    )
    assert r36 == R36, r36
    assert same_h and same_n
    for label, _sch, D in RUNS:
        for src in (pdd_m4["th"][label], pdd_m4["s"][label]):
            for arr in (src["h"], src["n"]):
                for j in range(arr.shape[1]):
                    assert float(np.max(_T(arr[:, j]))) < 43.0, (label, j)
                    assert float(np.max(_T(arr[:, j] / D))) < 43.0, (label, j)


@complex_only
def test_string_reproductions(pdd_m4):
    th = pdd_m4["th"]
    got39 = tuple(
        f"{cem43(T_GRID, _T(th['default']['h'][:, j])):.10e}" for j in range(len(DOSE_PTS))
    )
    got43 = tuple(
        f"{cem43(T_GRID, _T(th['default']['n'][:, j])):.10e}" for j in range(len(NEW_PTS))
    )
    eq = {}
    for label, _sch, D in RUNS[1:]:
        g = []
        for j in range(len(HIST_PTS)):
            dc_h = cem43(T_GRID, _T(th["default"]["h"][:, j]))
            de_h = cem43(T_GRID, _T(th[label]["h"][:, j] / D))
            g.append(f"{de_h / dc_h:.6f}")
        eq[label] = tuple(g)
    _say(f"MAT21-16 (R39') default HIST_PTS FEM doses {got39} (P39' {P39P_FEM_DOSE})")
    _say(f"MAT21-16 (R43) default NEW_PTS FEM doses {got43} (P43 {P43_FEM_DOSE})")
    _say(f"MAT21-16 (R41-h) D = 0.5 FEM equal-mean {eq['d50']} (P41 {P41_FEM_EQ})")
    _say(f"MAT21-16 (R42-h) D = 0.25 FEM equal-mean {eq['d25']} (P42 {P42_FEM_EQ})")
    assert got39 == tuple(P39P_FEM_DOSE), got39
    assert got43 == P43_FEM_DOSE, got43
    assert eq["d50"] == tuple(P41_FEM_EQ), eq["d50"]
    assert eq["d25"] == P42_FEM_EQ, eq["d25"]


@complex_only
def test_q_hp_win_hp_new_points(pdd_m4):
    for label, _sch, D in RUNS[1:]:
        th_d = pdd_m4["th"][label]["n"]
        for scale, lab in ((1.0, "peak"), (1.0 / D, "equal-mean")):
            for j, pt in enumerate(NEW_PTS):
                th = th_d[:, j] * scale
                dose, ref, rel, swap, miss = _q_check(_T(th))
                mx, t_mx, lo, hi, wref, wrel = _win_check(th)
                _say(
                    f"MAT21-16 (Q-hP) {label} {lab} {NEW_LABELS[j]}: FEM CEM43 {dose:.15e} s vs "
                    f"quad {ref:.15e} s rel {rel:.3e} (bar {DOSE_BAR:g}, ASSERTED)"
                )
                _say(
                    f"MAT21-16 control (Rc'-hP) PREDICTED {label} {lab} {NEW_LABELS[j]}: R swap "
                    f"miss {miss:.4e}; factor {swap / ref:.4e} (predicted range {PRED_RC}; "
                    f"printed only, ruling 18:00 (3))"
                )
                _say(
                    f"MAT21-16 (Win-hP) {label} {lab} {NEW_LABELS[j]}: max six-minute mean "
                    f"{mx:.15e} K at {t_mx:.4f} s; (ii) [{lo:.12e}, {hi:.12e}]; (iii) quad "
                    f"{wref:.15e} rel {wrel:.3e} (bar {AVG_BAR:g}, ASSERTED)"
                )
                assert rel <= DOSE_BAR, (label, lab, pt, rel)
                assert T_GRID[0] + W_B2 <= t_mx <= T_GRID[-1], (label, lab, t_mx)
                assert lo <= mx <= hi, (label, lab, lo, mx, hi)
                assert wrel <= AVG_BAR, (label, lab, pt, wrel)


@complex_only
def test_rc_hp_hist_points(pdd_m4):
    for label, _sch, D in RUNS[1:]:
        th_d = pdd_m4["th"][label]["h"]
        for scale, lab in ((1.0, "peak"), (1.0 / D, "equal-mean")):
            for j, pt in enumerate(DOSE_PTS):
                _dose, ref, _rel, swap, miss = _q_check(_T(th_d[:, j] * scale))
                _say(
                    f"MAT21-16 (Rc'-hP) {label} {lab} {HIST_LABELS[j]}: R swap miss {miss:.4e} "
                    f"(asserted > {DOSE_BAR:g}); factor {swap / ref:.4e}"
                )
                assert miss > DOSE_BAR, (label, lab, pt, miss)


def _bd_rows(pdd_m4, label, D):
    rows = []
    for scale, lab in ((1.0, "peak"), (1.0 / D, "equal-mean")):
        for b, j, name in SIX:
            h = pdd_m4["th"][label][b] * scale
            s = pdd_m4["s"][label][b] * scale
            d_h, d_s, m, M, cap = _bd(h, s, j)
            rows.append((lab, name, d_h, d_s, m, M, cap))
    return rows


def _bd_test(pdd_m4, label, D, tag):
    rows = _bd_rows(pdd_m4, label, D)
    for lab, name, d_h, d_s, m, M, cap in rows:
        kz = "z" if name.startswith("P5") else "x"
        new = name.startswith(("P3", "P4", "P5"))
        pred = (
            f"predicted M {PRED_M[label][kz]} K at peak, headroom {PRED_HEAD[label][kz]}"
            if new
            else "P0-P2: reprints P41 / P42"
        )
        extra = ""
        if label == "d25" and name.startswith("P5"):
            extra = (
                f"; M_p/(INTERIOR_BAND max theta_s) {M / cap:.4f} (PREDICTED RED range "
                f"{PRED_RED_BD25_P5})"
            )
        _say(
            f"MAT21-16 ({tag}) {lab} {name}: D_h {d_h:.10e} s, D_s {d_s:.15e} s; (i) m_p {m:.4e} "
            f"<= ln4*M_p {LN4 * M:.4e} (exact e^(ln4 M)-1 {math.expm1(LN4 * M):.4e}, printed "
            f"only); ratio {m / (LN4 * M):.4f} (predicted {PRED_RATIO}); (ii) M_p {M:.4e} K <= "
            f"INTERIOR_BAND*max theta_s {cap:.4e} K, headroom {cap / M:.2f}x ({pred}){extra}"
        )
    bad = [
        (lab, name, m, LN4 * M, M, cap)
        for lab, name, d_h, d_s, m, M, cap in rows
        if not (m <= LN4 * M and M <= cap)
    ]
    assert not bad, (f"({tag})", bad)


@complex_only
def test_bd50_dose_sensitivity_bound_six_points(pdd_m4):
    _bd_test(pdd_m4, "d50", DUTY, "Bd50")


@complex_only
def test_bd25_dose_sensitivity_bound_six_points(pdd_m4):
    _bd_test(pdd_m4, "d25", DUTY25, "Bd25")


@complex_only
def test_wrong_schedule_controls(pdd_m4):
    th, s = pdd_m4["th"], pdd_m4["s"]
    rows = []
    for j, name in enumerate(NEW_LABELS):
        Mw = float(np.max(np.abs(th["d25"]["n"][:, j] - s["d50"]["n"][:, j])))
        r = Mw / float(np.max(s["d25"]["n"][:, j]))
        _say(
            f"MAT21-16 control (c-50-dose) PREDICTED {name}: M_w / max theta_s25 {r:.4f} "
            f"(asserted > INTERIOR_BAND {INTERIOR_BAND:g}; predicted {PRED_C50})"
        )
        rows.append(("c-50", name, r))
    for b, j, name in SIX:
        Mw = float(np.max(np.abs(th["d50"][b][:, j] - s["d25"][b][:, j])))
        r = Mw / float(np.max(s["d50"][b][:, j]))
        _say(
            f"MAT21-16 control (c-25-dose) PREDICTED {name}: M_w / max theta_s50 {r:.4f} "
            f"(asserted > INTERIOR_BAND {INTERIOR_BAND:g}; predicted {PRED_C25}; ratio to band "
            f"{r / INTERIOR_BAND:.1f})"
        )
        rows.append(("c-25", name, r))
    bad = [x for x in rows if not x[2] > INTERIOR_BAND]
    assert not bad, bad


@complex_only
def test_printed_records(pdd_m4):
    """Printed records (no band); asserts only that they were formed."""
    th, s = pdd_m4["th"], pdd_m4["s"]
    W = float(W_B2)
    for b, j, name in SIX:
        cells = []
        for label, _sch, _D in RUNS:
            M = float(np.max(np.abs(th[label][b][:, j] - s[label][b][:, j])))
            cap = INTERIOR_BAND * float(np.max(s[label][b][:, j]))
            cells.append(f"{label} {cap / M:.2f}x")
        _say(f"MAT21-16 (Bd)(ii) HEADROOM MATRIX (peak) {name}: " + "; ".join(cells))
    for label, _sch, D in RUNS[1:]:
        d = [cem43(T_GRID, _T(th[label]["n"][:, j])) for j in range(len(NEW_PTS))]
        _say(
            f"MAT21-16 {label} +x / -x FEM dose asymmetry {d[0] / d[1] - 1.0:+.4e} (predicted "
            f"|.| <= {PRED_ASYM:g})"
        )
        for b, j, name in SIX + (("h", 3, "(0,0.9B,0) [printed only]"),):
            dc_h = cem43(T_GRID, _T(th["default"][b][:, j]))
            dd_h = cem43(T_GRID, _T(th[label][b][:, j]))
            de_h = cem43(T_GRID, _T(th[label][b][:, j] / D))
            dc_s = cem43(T_GRID, _T(s["default"][b][:, j]))
            dd_s = cem43(T_GRID, _T(s[label][b][:, j]))
            de_s = cem43(T_GRID, _T(s[label][b][:, j] / D))
            _a, mx_h, t_h = running_average(T_GRID, th[label][b][:, j], window=W)
            a_s, mx_s, t_s = running_average(T_GRID, s[label][b][:, j], window=W)
            a_h = _a
            good = ~np.isnan(a_h)
            da = float(np.max(np.abs(a_h[good] - a_s[good])))
            _say(
                f"MAT21-16 DUTY DOSE {label} {name}: FEM ratio equal-mean {de_h / dc_h:.6f} / peak "
                f"{dd_h / dc_h:.6f} (W1 closed form {de_s / dc_s:.6f} / {dd_s / dc_s:.6f}); FEM "
                f"dose peak {dd_h:.10e} s / equal-mean {de_h:.10e} s; six-minute peak FEM "
                f"{mx_h:.9f} K at {t_h:.2f} s vs series {mx_s:.9f} K at {t_s:.2f} s; "
                f"(FEM-c-win) max|a_h - a_s| {da:.4e} K (printed only)"
            )
            assert np.isfinite(de_h / dc_h)
    for norm in ("peak", "equal-mean"):
        for b, j, name in SIX:
            row = []
            for label, _sch, D in RUNS:
                sc = 1.0 if norm == "peak" else 1.0 / D
                row.append(f"{label} {cem43(T_GRID, _T(th[label][b][:, j] * sc)):.10e}")
            _say(
                f"MAT21-16 SIX-POINT x THREE-SCHEDULE FEM DOSE TABLE ({norm}) {name}: "
                + "; ".join(row)
                + " s"
            )
