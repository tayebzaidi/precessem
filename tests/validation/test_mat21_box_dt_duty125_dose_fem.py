"""`MAT-21` step 17, W2 — `MAT-20`'s dose and six-minute window and P39′'s (Bd)
bound on the gel box's FEM point histories at all six F2182 probe points under
a D = 0.125 schedule at P38's period (P36's source; P2, m = 4), at Δt = t*/40
(N = 40) and Δt = t*/80 (N = 80) (§9 item 592, §10 chain P45; brief
`docs/planning/chunks/MAT-21.md` "§9 item 592", review ruling 2026-10-22 03:00).

Sets no band and moves none; every band is imported by name.  Three CN runs:
``default40`` and ``d125_40`` through the imported ``_cn_snap`` (P44's bytes),
``d125_80`` through the local ``_cn_snap_n`` (N = 80 only).  Each run is read on
two fixed batches per snapshot — ``HIST_PTS`` then ``NEW_PTS`` — never merged.
Asserted: preconditions; (R39′) / (R43) on ``default40``; (Q-h) / (Win-h) at
both Δt; (Rc′-h) at t*/40; (Bd125-80) (i) and (ii) at all six points
(first-order form), every row printed before any is asserted; (c-25-on-125)
> ``INTERIOR_BAND`` at both Δt.  Printed only: (Bd125-40) with the weekly's
``PREDICTED RED:`` line as a print tag (ruling 03:00 (1), (2)); (Rc′-h) at
t*/80; the Δt headroom table; (FEM-X).
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
from tests.validation.test_mat17_duty_cycle_fem import _census
from tests.validation.test_mat21_box_closed_form import B_BOX, C_BOX, Q
from tests.validation.test_mat21_box_closed_form import T_INF as T_INF_0
from tests.validation.test_mat21_box_dose_closed_form import (
    AVG_BAR,
    DOSE_BAR,
    DOSE_PTS,
    HIST_LABELS,
    HIST_PTS,
    T_GRID,
    T_STAR,
    W_B2,
    _series_hist,
)
from tests.validation.test_mat21_box_dose_fem import R36, _bd
from tests.validation.test_mat21_box_duty25_dose_closed_form import T_SW25
from tests.validation.test_mat21_box_duty_closed_form import DT_BOX
from tests.validation.test_mat21_box_duty_dose_closed_form import _duty_hist
from tests.validation.test_mat21_box_duty_dose_fem import LN4, P39P_FEM_DOSE, _cn_snap
from tests.validation.test_mat21_box_dt_duty125_dose_closed_form import (
    DT80,
    DUTY125,
    N80,
    SCHEDULE_125,
    SCHEDULE_125_80,
    T_GRID80,
    T_SW125,
    _bd_g,
    _chi125,
    _duty_hist_g,
    _q_check_g,
    _win_check_g,
)
from tests.validation.test_mat21_box_fem_p2 import DEGREE, _dofs_pred, _rung_p
from tests.validation.test_mat21_box_probe_dose_closed_form import NEW_LABELS, NEW_PTS
from tests.validation.test_mat21_box_probe_duty_dose_fem import P43_FEM_DOSE
from tests.validation.test_mat21_box_probe_points_fem import PROBE_PTS, _probe_vals
from tests.validation.test_mat21_box_separable_fem import M4

SIX = tuple(("h", j, HIST_LABELS[j]) for j in range(len(DOSE_PTS))) + tuple(
    ("n", j, NEW_LABELS[j]) for j in range(len(NEW_PTS)))
SIX_KIND = (0, 0, 1, 2, 2, 3)  # +-y walls, centre, +-x, P5
KIND_NAME = ("+-y", "centre", "+-x", "P5")
# P44's (Bd)(ii) headroom matrix (peak), copied as text from
# docs/testing/logs/20261021T005339Z_MAT-21-step16-W2-probe-duty-dose-fem.log:183-188.
P44_HEAD = ("default 24.71x; d50 11.22x; d25 5.40x", "default 24.71x; d50 11.22x; d25 5.40x",
            "default 26.52x; d50 15.15x; d25 7.88x", "default 6.94x; d50 6.65x; d25 5.86x",
            "default 6.94x; d50 6.65x; d25 5.86x", "default 5.04x; d50 3.05x; d25 1.72x")
# Predictions written first (brief), printed only.
PRED_HEAD40 = ((1.5, 4.0), (2.5, 6.0), (3.0, 7.5), (0.7, 1.4))
PRED_HEAD80_FACTOR = (1.0, 4.0)
PRED_HEAD80_P5 = (0.8, 4.0)
PRED_RATIO_PX = (0.85, 1.0)
PRED_RC = (1e6, 1e7)
PRED_C25 = (0.6, 1.3)
PRED_ASYM = 1e-12
# The weekly's pre-registration (2026-10-21, P45), carried verbatim as a print tag on the
# printed (Bd)(ii) cell at P5 at t*/40 (ruling 2026-10-22 03:00 (1), (2): the quantity is
# M_p / (INTERIOR_BAND·max θ_s); headroom printed beside it).
PRED_RED_LINE = "PREDICTED RED: (Bd)(ii) at P5 at t*/40 range [0.8, 1.4]"
PRED_RED_RANGE = (0.8, 1.4)
RUNS = (("default40", 1.0), ("d125_40", DUTY125), ("d125_80", DUTY125))


class _Clock17:
    """Per-phase wall clock, printed on rank 0 with flush, prefix MAT21-17."""

    def __init__(self, tag):
        self.tag = tag
        self.t = time.perf_counter()

    def lap(self, phase):
        now = time.perf_counter()
        _say(f"MAT21-17 [{self.tag}] clock {phase}: {now - self.t:.1f} s")
        self.t = now


def _cn_snap_n(msh, ft, q_quad, schedule, n):
    """``_cn_snap``'s body with ``dt = T_STAR / n``; used only for N = 80."""
    kw = {} if schedule is None else {"source_schedule": schedule}
    return solve_transient_conduction(
        msh, ft, k=K_GEL, rho_c=RHO_C, q=q_quad, h=H_GEL, t_inf=T_INF, t_init=T_INF,
        robin_tag=ROBIN_TAG, dt=T_STAR / n, n_steps=n, theta=0.5,
        startup_half_steps=2, degree=DEGREE, snapshot_steps=tuple(range(n + 1)), **kw)


@pytest.fixture(scope="module")
def dt125_m4():
    comm = MPI.COMM_WORLD
    t0 = time.perf_counter()
    clk = _Clock17("dt125")
    out = {"imag": 0.0, "E": {}, "ndof": {}, "th": {}, "ok": {}, "sched": {}, "census": {},
           "fin": {}, "nsteps": {}}
    _say(f"\nMAT21-17 DEGREE = {DEGREE}; m = {M4}; N = {N_MAIN} / {N80}; ranks {comm.size}; "
         f"HIST_PTS {HIST_PTS}; NEW_PTS {NEW_PTS} (two fixed batches)")

    msh, ft, ns, n, v_h, st = _rung_p(M4, clk, DEGREE)
    out["ndof"]["uniform"] = int(st.temperature.function_space.dofmap.index_map.size_global)

    xs = ufl.SpatialCoordinate(msh)
    q_quad = Q * 3.0 * (xs[1] * xs[1] + xs[2] * xs[2]) / (B_BOX**2 + C_BOX**2)
    p_h = np.asarray(HIST_PTS, dtype=float)
    p_n = np.asarray(NEW_PTS, dtype=float)

    for label, _D in RUNS:
        if label == "default40":
            r, nst, dt = _cn_snap(msh, ft, q_quad, None), N_MAIN, DT_BOX
        elif label == "d125_40":
            r, nst, dt = _cn_snap(msh, ft, q_quad, _chi125), N_MAIN, DT_BOX
        else:
            r, nst, dt = _cn_snap_n(msh, ft, q_quad, _chi125, N80), N80, DT80
        clk.lap(f"CN N={nst} {label} ({nst + 1} snapshots)")
        out["nsteps"][label] = nst
        out["ndof"][label] = int(r.temperature.function_space.dofmap.index_map.size_global)
        out["imag"] = max(out["imag"], _max_imag(r.temperature))
        out["E"][label] = float(r.max_balance_residual)
        if label != "default40":
            _census(r, dt, 2, DUTY125 * T_STAR, f"quad {label}")
            out["census"][label] = True
            out["sched"][label] = "".join(str(int(x)) for x in r.source_multipliers)
        th_h = np.empty((nst + 1, len(HIST_PTS)))
        th_n = np.empty((nst + 1, len(NEW_PTS)))
        ok_h = np.ones(len(HIST_PTS), dtype=bool)
        ok_n = np.ones(len(NEW_PTS), dtype=bool)
        for k in range(nst + 1):
            a, oa = _probe_vals(r.snapshots[k], p_h)
            b, ob = _probe_vals(r.snapshots[k], p_n)
            th_h[k], th_n[k] = a, b
            ok_h &= oa
            ok_n &= ob
        if label in ("default40", "d125_80"):
            fin_h, okf_h = _probe_vals(r.temperature, p_h)
            fin_n, okf_n = _probe_vals(r.temperature, p_n)
            ok_h &= okf_h
            ok_n &= okf_n
            out["fin"][label] = (fin_h, fin_n)
        r.snapshots.clear()
        del r
        out["th"][label] = {"h": th_h, "n": th_n}
        out["ok"][label] = (ok_h, ok_n)
        clk.lap(f"{label}: 2 x {nst + 1} owner-route batches")

    out["s"] = {
        "default40": {"h": _series_hist(HIST_PTS), "n": _series_hist(NEW_PTS)},
        "d125_40": {"h": _duty_hist(HIST_PTS, switch_times=T_SW125),
                    "n": _duty_hist(NEW_PTS, switch_times=T_SW125)},
        "d25_40": {"h": _duty_hist(HIST_PTS, switch_times=T_SW25),
                   "n": _duty_hist(NEW_PTS, switch_times=T_SW25)},
        "d125_80": {"h": _duty_hist_g(HIST_PTS, T_GRID80, T_SW125),
                    "n": _duty_hist_g(NEW_PTS, T_GRID80, T_SW125)},
        "d25_80": {"h": _duty_hist_g(HIST_PTS, T_GRID80, T_SW25),
                   "n": _duty_hist_g(NEW_PTS, T_GRID80, T_SW25)},
    }
    clk.lap("series histories (continuous; D = 0.125, 0.25 at t*/40 and t*/80; two batches)")
    _say(f"MAT21-17 fixture elapsed {time.perf_counter() - t0:.1f} s; max|Im T| {out['imag']:.3e}")
    return out


def _T(th):
    return T_INF + np.asarray(th, dtype=float)


def _grid(label):
    return T_GRID80 if label.endswith("_80") else T_GRID


@complex_only
def test_preconditions(dt125_m4):
    d = dt125_m4
    pred = _dofs_pred(M4)
    _say(f"MAT21-17 (Pdeg) P2 dofs {d['ndof']} (predicted {pred} = 116025)")
    assert pred == 116025
    assert len(d["ndof"]) == 4
    for key, v in d["ndof"].items():
        assert v == pred, f"(Pdeg) {key}: {v}"
    assert tuple(NEW_PTS) == tuple(PROBE_PTS[3:6])
    _say(f"MAT21-17 (S0) schedules {d['sched']} (must equal d125_40 {SCHEDULE_125}, d125_80 "
         f"{SCHEDULE_125_80})")
    assert d["census"] == {"d125_40": True, "d125_80": True}
    assert d["sched"]["d125_40"] == SCHEDULE_125 and d["sched"]["d125_80"] == SCHEDULE_125_80
    for key, (ok_h, ok_n) in d["ok"].items():
        _say(f"MAT21-17 validity {key}: HIST_PTS {ok_h.tolist()} NEW_PTS {ok_n.tolist()}")
        assert bool(np.all(ok_h)) and bool(np.all(ok_n)), key
    _say("MAT21-17 (E) max balance residual: " + ", ".join(
        f"{k_} {v:.3e}" for k_, v in d["E"].items()) + f" (bound {E_BAND:.0e}; predicted 2.3e-13)")
    assert len(d["E"]) == 3
    for key, v in d["E"].items():
        assert v <= E_BAND, (key, v)
    assert T_INF == T_INF_0
    th_d = d["th"]["default40"]
    r36 = (f"{th_d['h'][40, 2]:.6f}", f"{th_d['h'][40, 3]:.6f}")
    _say(f"MAT21-17 (R36) default40 theta_h(t*) centre / (0,0.9B,0) {r36} (P36 {R36})")
    assert r36 == R36, r36
    for label in ("default40", "d125_80"):
        nst = d["nsteps"][label]
        fin_h, fin_n = d["fin"][label]
        same_h = bool(np.array_equal(d["th"][label]["h"][nst], fin_h))
        same_n = bool(np.array_equal(d["th"][label]["n"][nst], fin_n))
        _say(f"MAT21-17 {label} snapshot n={nst} array_equal final: HIST_PTS {same_h}, "
             f"NEW_PTS {same_n}")
        assert same_h and same_n, label
    for label, D in RUNS:
        for src in (d["th"][label], d["s"][label]):
            for arr in (src["h"], src["n"]):
                for j in range(arr.shape[1]):
                    assert float(np.max(_T(arr[:, j]))) < 43.0, (label, j)
                    assert float(np.max(_T(arr[:, j] / D))) < 43.0, (label, j)


@complex_only
def test_string_reproductions(dt125_m4):
    th = dt125_m4["th"]["default40"]
    got39 = tuple(f"{cem43(T_GRID, _T(th['h'][:, j])):.10e}" for j in range(len(DOSE_PTS)))
    got43 = tuple(f"{cem43(T_GRID, _T(th['n'][:, j])):.10e}" for j in range(len(NEW_PTS)))
    _say(f"MAT21-17 (R39') default40 HIST_PTS FEM doses {got39} (P39' {P39P_FEM_DOSE})")
    _say(f"MAT21-17 (R43) default40 NEW_PTS FEM doses {got43} (P43 {P43_FEM_DOSE})")
    assert got39 == tuple(P39P_FEM_DOSE), got39
    assert got43 == P43_FEM_DOSE, got43


@complex_only
def test_q_h_win_h_both_dt(dt125_m4):
    bad = []
    for label in ("d125_40", "d125_80"):
        t = _grid(label)
        for scale, lab in ((1.0, "peak"), (1.0 / DUTY125, "equal-mean")):
            for b, j, name in SIX:
                th = dt125_m4["th"][label][b][:, j] * scale
                dose, ref, rel, swap, miss = _q_check_g(_T(th), t)
                mx, t_mx, lo, hi, wref, wrel = _win_check_g(th, t)
                _say(f"MAT21-17 (Q-h) {label} {lab} {name}: FEM CEM43 {dose:.15e} s vs quad "
                     f"{ref:.15e} s rel {rel:.3e} (bar {DOSE_BAR:g}, ASSERTED); max T "
                     f"{np.max(_T(th)):.9f} C")
                _say(f"MAT21-17 (Win-h) {label} {lab} {name}: max six-minute mean {mx:.15e} K at "
                     f"{t_mx:.4f} s; (ii) [{lo:.12e}, {hi:.12e}]; (iii) quad {wref:.15e} rel "
                     f"{wrel:.3e} (bar {AVG_BAR:g}, ASSERTED)")
                if not (rel <= DOSE_BAR and t[0] + W_B2 <= t_mx <= t[-1] and lo <= mx <= hi
                        and wrel <= AVG_BAR):
                    bad.append((label, lab, name, rel, t_mx, lo, mx, hi, wrel))
    assert not bad, bad


@complex_only
def test_rc_h_controls(dt125_m4):
    bad = []
    for label in ("d125_40", "d125_80"):
        t = _grid(label)
        for scale, lab in ((1.0, "peak"), (1.0 / DUTY125, "equal-mean")):
            for b, j, name in SIX:
                th = dt125_m4["th"][label][b][:, j] * scale
                _d, ref, _r, swap, miss = _q_check_g(_T(th), t)
                tag = (f"asserted > {DOSE_BAR:g}" if label == "d125_40" else
                       f"PREDICTED factor {PRED_RC}; printed only, ruling 03:00 (6)")
                _say(f"MAT21-17 (Rc'-h) {label} {lab} {name}: R swap miss {miss:.4e} ({tag}); "
                     f"factor {swap / ref:.4e}")
                if label == "d125_40" and not miss > DOSE_BAR:
                    bad.append((lab, name, miss))
    assert not bad, bad


def _bd_rows(d, label):
    rows = []
    t = _grid(label)
    for scale, lab in ((1.0, "peak"), (1.0 / DUTY125, "equal-mean")):
        for b, j, name in SIX:
            h = d["th"][label][b] * scale
            s = d["s"][label][b] * scale
            if label == "d125_40":
                d_h, d_s, m, M, cap = _bd(h, s, j)
            else:
                d_h, d_s, m, M, cap = _bd_g(h, s, j, t)
            rows.append((lab, b, j, name, d_h, d_s, m, M, cap))
    return rows


def _bd_print(rows, label, tag):
    for lab, _b, j, name, d_h, d_s, m, M, cap in rows:
        extra = ""
        if label == "d125_40" and name.startswith("P5"):
            r = M / cap
            inside = "inside" if PRED_RED_RANGE[0] <= r <= PRED_RED_RANGE[1] else "outside"
            extra = (f"; M_p/(INTERIOR_BAND max theta_s) {r:.4f}, headroom {cap / M:.4f} "
                     f"[{PRED_RED_LINE}: {inside}]")
        _say(f"MAT21-17 ({tag}) {lab} {name}: D_h {d_h:.10e} s, D_s {d_s:.15e} s; (i) m_p "
             f"{m:.4e} <= ln4*M_p {LN4 * M:.4e} (exact e^(ln4 M)-1 {math.expm1(LN4 * M):.4e}, "
             f"printed only); ratio {m / (LN4 * M):.4f}; (ii) M_p {M:.4e} K <= "
             f"INTERIOR_BAND*max theta_s {cap:.4e} K, headroom {cap / M:.2f}x{extra}")


@complex_only
def test_bd125_40_printed_only(dt125_m4):
    """(Bd125-40): printed only (the weekly's text); carries the PREDICTED RED tag."""
    rows = _bd_rows(dt125_m4, "d125_40")
    _bd_print(rows, "d125_40", "Bd125-40, printed only")
    bad = [(lab, name) for lab, _b, _j, name, _dh, _ds, m, M, cap in rows
           if not (m <= LN4 * M and M <= cap)]
    _say(f"MAT21-17 (Bd125-40) cells outside the bound (printed only, not asserted): {bad}")
    assert len(rows) == 12


@complex_only
def test_bd125_80_dose_sensitivity_bound_six_points(dt125_m4):
    rows = _bd_rows(dt125_m4, "d125_80")
    _bd_print(rows, "d125_80", "Bd125-80")
    bad = [(lab, name, m, LN4 * M, M, cap) for lab, _b, _j, name, _dh, _ds, m, M, cap in rows
           if not (m <= LN4 * M and M <= cap)]
    assert not bad, ("(Bd125-80)", bad)


@complex_only
def test_c25_on_125_controls(dt125_m4):
    th, s = dt125_m4["th"], dt125_m4["s"]
    rows = []
    for label, s25 in (("d125_40", "d25_40"), ("d125_80", "d25_80")):
        for b, j, name in SIX:
            Mw = float(np.max(np.abs(th[label][b][:, j] - s[s25][b][:, j])))
            r = Mw / float(np.max(s[label][b][:, j]))
            _say(f"MAT21-17 control (c-25-on-125) PREDICTED {label} {name}: M_w / max theta_s125 "
                 f"{r:.4f} (asserted > INTERIOR_BAND {INTERIOR_BAND:g}; predicted {PRED_C25})")
            rows.append((label, name, r))
    bad = [x for x in rows if not x[2] > INTERIOR_BAND]
    assert not bad, bad


@complex_only
def test_printed_records(dt125_m4):
    """Printed records (no band); asserts only that they were formed."""
    th, s = dt125_m4["th"], dt125_m4["s"]
    W = float(W_B2)
    for i, (b, j, name) in enumerate(SIX):
        cells = {}
        for label in ("default40", "d125_40", "d125_80"):
            M = float(np.max(np.abs(th[label][b][:, j] - s[label][b][:, j])))
            mx = float(np.max(s[label][b][:, j]))
            cells[label] = (M, mx, INTERIOR_BAND * mx / M)
        k = SIX_KIND[i]
        p80 = PRED_HEAD80_P5 if k == 3 else PRED_HEAD80_FACTOR
        _say(f"MAT21-17 DT HEADROOM TABLE (peak) {name}: P44 t*/40 [{P44_HEAD[i]}]; default40 "
             f"{cells['default40'][2]:.2f}x; d125 t*/40 M_p {cells['d125_40'][0]:.4e} K, max "
             f"theta_s {cells['d125_40'][1]:.4e} K, headroom {cells['d125_40'][2]:.4f}x "
             f"(predicted {PRED_HEAD40[k]}); d125 t*/80 M_p {cells['d125_80'][0]:.4e} K, max "
             f"theta_s {cells['d125_80'][1]:.4e} K, headroom {cells['d125_80'][2]:.4f}x "
             f"(predicted {'' if k == 3 else 'factor '}{p80}); headroom80/headroom40 "
             f"{cells['d125_80'][2] / cells['d125_40'][2]:.4f}; M_p40/M_p80 "
             f"{cells['d125_40'][0] / cells['d125_80'][0]:.4f}")
    for b, j, name in SIX:
        a40 = th["d125_40"][b][:, j]
        a80 = th["d125_80"][b][:, j]
        x = float(np.max(np.abs(a80[::2] - a40)))
        _say(f"MAT21-17 (FEM-X) {name}: max_n |theta80[2n] - theta40[n]| {x:.4e} K "
             f"({x / float(np.max(a40)):.4e} of max theta40)")
    for label in ("d125_40", "d125_80"):
        t = _grid(label)
        d = [cem43(t, _T(th[label]["n"][:, j])) for j in range(len(NEW_PTS))]
        _say(f"MAT21-17 {label} +x / -x FEM dose asymmetry {d[0] / d[1] - 1.0:+.4e} (predicted "
             f"|.| <= {PRED_ASYM:g})")
        for b, j, name in SIX + (("h", 3, "(0,0.9B,0) [printed only]"),):
            dc_h = cem43(T_GRID, _T(th["default40"][b][:, j]))
            dc_s = cem43(T_GRID, _T(s["default40"][b][:, j]))
            dd_h = cem43(t, _T(th[label][b][:, j]))
            de_h = cem43(t, _T(th[label][b][:, j] / DUTY125))
            dd_s = cem43(t, _T(s[label][b][:, j]))
            de_s = cem43(t, _T(s[label][b][:, j] / DUTY125))
            _a, mx_h, t_h = running_average(t, th[label][b][:, j], window=W)
            _a, mx_s, t_s = running_average(t, s[label][b][:, j], window=W)
            _ae, mxe_h, _te = running_average(t, th[label][b][:, j] / DUTY125, window=W)
            _say(f"MAT21-17 DUTY DOSE {label} {name}: FEM ratio equal-mean {de_h / dc_h:.6f} / "
                 f"peak {dd_h / dc_h:.6f} (series {de_s / dc_s:.6f} / {dd_s / dc_s:.6f}; cont at "
                 f"t*/40); FEM dose peak {dd_h:.10e} s / equal-mean {de_h:.10e} s; six-minute "
                 f"peak FEM {mx_h:.9f} K at {t_h:.2f} s vs series {mx_s:.9f} K at {t_s:.2f} s; "
                 f"equal-mean max T {float(np.max(_T(th[label][b][:, j] / DUTY125))):.4f} C")
            assert np.isfinite(de_h / dc_h)
