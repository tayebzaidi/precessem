"""`MAT-20` step 4 (§9 item 557, §10 chain P40) — dose and the six-minute window
on the tuned copper F-human's coil-driven gel heating at the degree-2 root.

Input: ``logs/mat15_s5/theta-series-ccw-d2.json`` (gitignored), written by
``xl`` entry 42 (`MAT-15` step 6, ``20261019T070008Z_MAT-15-step6-P22-d2root.log``
:6887–6929): the ccw CN run (N = 40 to t* = 1 254 s) with every step kept, at
the degree-2 root (115.5209 pF) — volume-mean, 10 g ball-mean at c*_θ and 5 mm
probe θ(t) series, per accepted W (ruling 361-1; the JSON already carries G).
Absent ⇒ every test skips (fresh clone, CI).  The degree-1 twin
(``theta-series-ccw.json``, item 380) is read beside it for the d1 → d2 prints.

Asserted (bars imported, none restated), on the 10 g ball-mean
(``theta_ball10g_k``) and the volume-mean (``theta_mean_k``) series,
T = t_inf_c + θ:

- (P0-d2) the file is entry 42's series (degree, node knob, C_in, 41 samples,
  six printed log strings); its time grid and T_INF equal the degree-1 file's.
- (Q-d2) ``cem43`` equals step 3's independent per-segment ``quad`` at
  ``DOSE_BAR``; max T < 43 °C.
- (Win-d2) the W = 360 s maximum: (i) inside [t₀ + W, t_N]; (ii) the
  mean-value sandwich (exact); (iii) equal to an independent ``quad`` of the
  interpolant at ``AVG_BAR``.

Negative controls (rule (e)): (Rc′-d2) the R swap misses (Q-d2)'s reference by
> ``DOSE_BAR`` (factor predicted ≈ 8.2e6, printed); (P0-ctl) the provenance
check rejects the degree-1 file.

Printed only (a record): CEM43, the dose excess over the T_INF baseline and the
six-minute peak per series and normalisation, d1 → d2 moves beside θ(t*)'s,
(Rect) and (Sh).  The 5 mm probe is printed only (it reads the mesh).

Not claimed: a dose or six-minute limit, compliance or a safety statement;
tissue, perfusion, σ(T); the ASTM phantom; an absolute temperature; a
convergence statement (two element orders on one mesh); any external check.
"""

from __future__ import annotations

import json
import math

import numpy as np
import pytest
from mpi4py import MPI
from scipy.integrate import quad

from fem_em_solver.thermal.dose import cem43, running_average

from tests.validation.test_mat20_cem43_closed_form import BAR as DOSE_BAR
from tests.validation.test_mat20_running_average import BAR as AVG_BAR
from tests.validation.test_mat20_running_average import W_B2
from tests.validation.test_mat20_tuned_heating_dose import (
    R_BELOW,
    SERIES_JSON,
    T_REF,
    _quad_dose,
    _rect,
)

SERIES_D2 = SERIES_JSON.with_name("theta-series-ccw-d2.json")

# Entry 42's printed series (log 20261019T070008Z_MAT-15-step6-P22-d2root.log,
# :6889 for n = 1 and :6928 for n = N), each the log's ``:.12e``.
E42_LOG = {
    1: {"theta_mean_k": "9.408837938992e-04", "theta_ball10g_k": "1.293205450941e-03",
        "theta_probe_k": "1.287310587223e-03"},
    -1: {"theta_mean_k": "3.507742275686e-02", "theta_ball10g_k": "4.916453675956e-02",
         "theta_probe_k": "4.979676417684e-02"},
}
E42_C_IN_PF = 115.5209
E42_SAMPLES = 41

ASSERTED = (("theta_ball10g_k", "10 g ball-mean"), ("theta_mean_k", "volume mean"))
ALL_SERIES = ASSERTED + (("theta_probe_k", "5 mm probe"),)

# Predictions written first (brief, §9 item 557), printed only.
PRED_Q = 1e-14
PRED_WIN_III = 1e-15
PRED_RC_FACTOR = 8.2e6
PRED_TABLE = {  # (dose move, excess move, peak move, theta(t*) move)
    "theta_mean_k": (-3.28e-4, -1.313e-2, -1.280e-2, -1.2745e-2),
    "theta_ball10g_k": (-9.2e-6, -2.6e-4, -1.49e-3, -1.96e-3),
    "theta_probe_k": (+6.0e-6, +1.7e-4, -1.11e-3, -1.61e-3),
}


def _p(*a):
    if MPI.COMM_WORLD.rank == 0:
        print(*a, flush=True)


def _is_d2_entry42(d) -> bool:
    if d.get("degree") != 2 or d.get("node_knob") != "d2root":
        return False
    if d.get("c_in_pf") != E42_C_IN_PF or len(d.get("times_s", [])) != E42_SAMPLES:
        return False
    for idx, strs in E42_LOG.items():
        for key, s in strs.items():
            if f"{d[key][idx]:.12e}" != s:
                return False
    return True


def _diff_report(d):
    out = [f"degree={d.get('degree')!r}", f"node_knob={d.get('node_knob')!r}",
           f"c_in_pf={d.get('c_in_pf')!r}", f"samples={len(d.get('times_s', []))}"]
    for idx, strs in E42_LOG.items():
        for key, s in strs.items():
            got = f"{d[key][idx]:.12e}" if key in d else "absent"
            if got != s:
                out.append(f"{key}[{idx}] {got} != {s}")
    return "; ".join(out)


@pytest.fixture(scope="module")
def pair():
    if not SERIES_D2.exists():
        pytest.skip(f"{SERIES_D2} absent (written by xl entry 42, MAT-15 step 6; never "
                    "re-run in-slot)")
    with open(SERIES_D2) as fh:
        d2 = json.load(fh)
    d1 = None
    if SERIES_JSON.exists():
        with open(SERIES_JSON) as fh:
            d1 = json.load(fh)
    t = np.asarray(d2["times_s"], dtype=float)
    _p(f"[MAT-20 s4] {SERIES_D2.name}: {t.size} samples, dt {d2['dt_s']:.6e} s, t_N "
       f"{t[-1]:.3f} s; degree {d2.get('degree')}, node_knob {d2.get('node_knob')}, C_in "
       f"{d2.get('c_in_pf')} pF; ranks at write {d2['ranks']}; cells {d2['cells']}, n_sub "
       f"{d2['n_sub']}; c* {d2['c_star_label']}; G {d2['G_per_w']:.6f} /W")
    _p(f"[MAT-20 s4] degree-1 twin {SERIES_JSON.name}: "
       f"{'present' if d1 is not None else 'ABSENT (d1 -> d2 prints skipped)'}")
    return d1, d2


def _series(d, key):
    t = np.asarray(d["times_s"], dtype=float)
    th = np.asarray(d[key], dtype=float)
    return t, th, float(d["t_inf_c"]) + th


def test_p0_provenance(pair):
    d1, d2 = pair
    ok = _is_d2_entry42(d2)
    _p(f"[MAT-20 s4] (P0-d2) _is_d2_entry42(d2) = {ok} ASSERTED True; {_diff_report(d2)}")
    assert ok, _diff_report(d2)
    if d1 is not None:
        same_t = d2["times_s"] == d1["times_s"]
        same_inf = d2["t_inf_c"] == d1["t_inf_c"]
        _p(f"[MAT-20 s4] (P0-d2) times_s equal d1: {same_t}; t_inf_c equal: {same_inf} "
           f"({d2['t_inf_c']} C) ASSERTED")
        assert same_t and same_inf


def test_p0_control_rejects_d1(pair):
    d1, _d2 = pair
    if d1 is None:
        pytest.skip(f"{SERIES_JSON} absent: (P0-ctl) needs the degree-1 file")
    ok = _is_d2_entry42(d1)
    _p(f"[MAT-20 s4] (P0-ctl) _is_d2_entry42(d1) = {ok} ASSERTED False; d1: "
       f"{_diff_report(d1)}")
    assert not ok


@pytest.mark.parametrize("key,label", ASSERTED)
def test_q_dose_equals_quadrature_d2(pair, key, label):
    _d1, d2 = pair
    t, _th, T = _series(d2, key)
    assert float(np.max(T)) < T_REF, float(np.max(T))
    dose = cem43(t, T)
    ref = _quad_dose(t, T, math.log(R_BELOW))
    rel = abs(dose - ref) / abs(ref)
    swap = cem43(t, T, r_above=R_BELOW, r_below=0.5)
    miss = abs(swap - ref) / abs(ref)
    rect = _rect(t, T)
    _p(f"[MAT-20 s4] (Q-d2) {label} [{key}] CEM43 {dose:.15e} s vs per-segment quad "
       f"{ref:.15e} s: rel {rel:.3e} (bar {DOSE_BAR:g}, ASSERTED; predicted <= {PRED_Q:g}); "
       f"max T {np.max(T):.9f} C < 43 ASSERTED")
    _p(f"[MAT-20 s4] (Rc'-d2) {label} [{key}] R swap {swap:.6e} s: miss {miss:.4e} relative "
       f"(asserted > {DOSE_BAR:g}; predicted factor ~{PRED_RC_FACTOR:.1e}, printed)")
    _p(f"[MAT-20 s4] (Rect) {label} [{key}] rectangle sum {rect:.15e} s, rel to (Q) "
       f"{rect / ref - 1.0:+.4e} (printed)")
    assert rel <= DOSE_BAR, f"(Q-d2) {key} {rel:.3e}"
    assert miss > DOSE_BAR, f"(Rc'-d2) {key} R swap did not miss: {miss:.3e}"


@pytest.mark.parametrize("key,label", ASSERTED)
def test_win_six_minute_max_inside_window_d2(pair, key, label):
    _d1, d2 = pair
    t, th, _T = _series(d2, key)
    W = float(W_B2)
    _means, mx, t_mx = running_average(t, th, window=W)
    s = t_mx - W
    inside = t[(t > s) & (t <= t_mx)]
    pts = np.concatenate([[np.interp(s, t, th)], np.interp(inside, t, th)])
    lo, hi = float(np.min(pts)), float(np.max(pts))
    brk = [x for x in t if s < x < t_mx]
    ref, _err = quad(lambda x: np.interp(x, t, th), s, t_mx, points=brk or None,
                     epsabs=0.0, epsrel=1e-13, limit=200)
    ref /= W
    rel = abs(mx - ref) / abs(ref)
    dt = float(d2["dt_s"])
    _m2, mx_sh, t_sh = running_average(t, th, window=W - dt)
    _p(f"[MAT-20 s4] (Win-d2) {label} [{key}] W {W:g} s ({W / dt:.4f} samples): max "
       f"six-minute mean {mx:.15e} K at t_max {t_mx:.6f} s (predicted t* {t[-1]:.3f}); (i) "
       f"[{t[0] + W:.3f}, {t[-1]:.3f}] ASSERTED; (ii) interpolant on window [{lo:.12e}, "
       f"{hi:.12e}] ASSERTED; (iii) quad {ref:.15e}, rel {rel:.3e} (bar {AVG_BAR:g}, "
       f"ASSERTED; predicted <= {PRED_WIN_III:g})")
    _p(f"[MAT-20 s4] (Sh) {label} [{key}] window W - dt = {W - dt:.6f} s: max {mx_sh:.12e} K "
       f"at t {t_sh:.6f} s (printed)")
    assert t[0] + W <= t_mx <= t[-1], t_mx
    assert lo <= mx <= hi, (lo, mx, hi)
    assert rel <= AVG_BAR, f"(Win-d2)(iii) {key} {rel:.3e}"


def _figures(d, key, scale=1.0):
    t = np.asarray(d["times_s"], dtype=float)
    th = scale * np.asarray(d[key], dtype=float)
    t_inf = float(d["t_inf_c"])
    dose = cem43(t, t_inf + th)
    base = R_BELOW ** (T_REF - t_inf) * (t[-1] - t[0])
    _m, peak, _tp = running_average(t, th, window=float(W_B2))
    return dose, dose - base, peak, th[-1], base


def test_moves_printed(pair):
    """Printed record (no limit, no compliance); asserts only that it was formed."""
    d1, d2 = pair
    p_acc = d2["p_acc_w"]["ccw"]
    _p(f"[MAT-20 s4] d2 P_acc,ccw {p_acc:.12e} W, P_VI,ccw {d2['p_vi_w']['ccw']:.12e} W, "
       f"P_src {d2['p_src_w']:.12e} W; G_per_w d2 {d2['G_per_w']:.6f}"
       + (f", d1 {d1['G_per_w']:.6f}" if d1 is not None else "")
       + " (series already per accepted W; G never re-applied)")
    formed = 0
    _p("[MAT-20 s4] MOVE TABLE per accepted W (ruling 361-1), d1 -> d2 (printed record; "
       "predicted moves from the brief beside; not a safety statement)")
    for key, label in ALL_SERIES:
        dose2, ex2, pk2, th2, base = _figures(d2, key)
        if d1 is None:
            _p(f"[MAT-20 s4] MOVE {label:<14s} [{key}]: d2 CEM43 {dose2:.6e} s = "
               f"{dose2 / 60:.6e} min, excess {ex2:.6e} s, peak {pk2:.6e} K (d1 absent: "
               "moves skipped)")
            formed += int(np.isfinite(dose2) and np.isfinite(pk2))
            continue
        dose1, ex1, pk1, th1, _b = _figures(d1, key)
        pr = PRED_TABLE[key]
        _p(f"[MAT-20 s4] MOVE {label:<14s} [{key}]: CEM43 {dose1:.6e} -> {dose2:.6e} s "
           f"({dose1 / 60:.6e} -> {dose2 / 60:.6e} min), dose move {dose2 / dose1 - 1:+.3e} "
           f"(pred {pr[0]:+.2e}); excess over baseline R^23 t* = {base:.6e} s: {ex1:.6e} -> "
           f"{ex2:.6e} s, move {ex2 / ex1 - 1:+.4e} (pred {pr[1]:+.2e}); six-minute peak "
           f"{pk1:.6e} -> {pk2:.6e} K, move {pk2 / pk1 - 1:+.4e} (pred {pr[2]:+.3e}); theta(t*) "
           f"{th1:.6e} -> {th2:.6e} K, move {th2 / th1 - 1:+.4e} (pred {pr[3]:+.4e})")
        formed += int(all(np.isfinite(x) for x in (dose1, dose2, ex1, ex2, pk1, pk2)))
    norms = {
        "per W of 1/2 Re sum V I*": lambda d: d["p_acc_w"]["ccw"] / d["p_vi_w"]["ccw"],
        "per W of P_src (not the pre-registered normaliser)":
            lambda d: d["p_acc_w"]["ccw"] / d["p_src_w"],
    }
    for norm, fsc in norms.items():
        for key, label in ALL_SERIES:
            sc2 = fsc(d2)
            dose2, ex2, pk2, th2, _b = _figures(d2, key, sc2)
            line = (f"[MAT-20 s4] READING {label:<14s} [{key}] {norm}: d2 (x {sc2:.6f}) CEM43 "
                    f"{dose2:.6e} s, peak {pk2:.6e} K, theta(t*) {th2:.6e} K")
            if d1 is not None:
                sc1 = fsc(d1)
                dose1, ex1, pk1, th1, _b = _figures(d1, key, sc1)
                line += (f"; d1 (x {sc1:.6f}) CEM43 {dose1:.6e} s, peak {pk1:.6e} K; moves: "
                         f"dose {dose2 / dose1 - 1:+.3e}, excess {ex2 / ex1 - 1:+.4e}, peak "
                         f"{pk2 / pk1 - 1:+.4e}, theta(t*) {th2 / th1 - 1:+.4e}")
            _p(line)
            formed += int(np.isfinite(dose2) and np.isfinite(pk2))
    assert formed == 9
