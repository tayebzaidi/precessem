"""`MAT-20` step 3 (§9 item 380, §10 chain P20) — dose and the six-minute window
on the tuned copper F-human's coil-driven gel heating.

Input: ``logs/mat15_s5/theta-series-ccw.json`` (gitignored), written by
`MAT-15` step 5's module under ``FEM_EM_MAT15_STEP5_SERIES=1``: the ccw CN run
(N = 40 to t* = 1 254 s) with every step kept — volume-mean, 10 g ball-mean at
c*_θ and 5 mm probe θ(t) series, per accepted W (ruling 361-1,
G = 1 / P_acc,ccw).  Absent ⇒ every test skips (fresh clone, CI).

Asserted on the 10 g ball-mean series at P★ = 1 W accepted, T = T_INF + θ̄
(bars imported, none restated):

- (Q) ``cem43`` equals an independent per-segment ``scipy.integrate.quad`` of
  R^(43 − T_lin) (R = 0.25; T < 43 asserted) at the step-1 bar ``DOSE_BAR``.
- (Win) for the W = 360 s running mean: (i) t₀ + W ≤ t_max ≤ t_N; (ii) the
  mean-value sandwich, min ≤ max-mean ≤ max of the interpolant over the
  window (exact); (iii) the max-mean equals an independent ``quad`` of the
  interpolant over [t_max − W, t_max] at the step-2 bar ``AVG_BAR``.

Negative control (rule (e)): (Rc′) the R swap misses (Q)'s reference by
> ``DOSE_BAR`` (factor predicted ≈ 8e6, printed).  Printed only: (Sh) the
window one sample short, (Rect) the rectangle sum, the six-minute peak and
CEM43 per W accepted / per W of ½Re ΣV·I* / per W of P_src, and the same two
figures on the volume-mean and probe series.

Not claimed: a dose or six-minute limit, compliance, tissue, perfusion, σ(T),
the ASTM phantom, degree 2, an absolute temperature, ANS-9.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pytest
from mpi4py import MPI
from scipy.integrate import quad

from fem_em_solver.thermal.dose import cem43, running_average

from tests.validation.test_mat20_cem43_closed_form import BAR as DOSE_BAR
from tests.validation.test_mat20_running_average import BAR as AVG_BAR
from tests.validation.test_mat20_running_average import W_B2

SERIES_JSON = Path(__file__).resolve().parents[2] / "logs" / "mat15_s5" / "theta-series-ccw.json"
R_BELOW = 0.25
T_REF = 43.0

# Predictions written first (brief, §9 item 380), printed only.
PRED_Q = 1e-14
PRED_WIN_III = 1e-15
PRED_RC_FACTOR = 8e6
PRED_DOSE_S = 1.8e-11
PRED_PVI_SCALE = 0.991


def _p(*a):
    if MPI.COMM_WORLD.rank == 0:
        print(*a, flush=True)


@pytest.fixture(scope="module")
def series():
    if not SERIES_JSON.exists():
        pytest.skip(f"{SERIES_JSON} absent (regenerate: MAT-15 step 5 with "
                    "FEM_EM_MAT15_STEP5_SERIES=1)")
    with open(SERIES_JSON) as fh:
        d = json.load(fh)
    t = np.asarray(d["times_s"], dtype=float)
    _p(f"[MAT-20 s3] {SERIES_JSON.name}: {t.size} samples, dt {d['dt_s']:.6e} s, t_N "
       f"{t[-1]:.3f} s; ranks at write {d['ranks']}; cells {d['cells']}, n_sub {d['n_sub']}; "
       f"c* {d['c_star_label']}; a10 {d['a10_m'] * 1e3:.4f} mm, probe {d['a_probe_m'] * 1e3:.1f} mm")
    return d


def _quad_dose(t, T, ln_r):
    """Independent per-segment quadrature of R^(T_REF − T_lin(t))."""
    tot = 0.0
    for i in range(t.size - 1):
        t0, t1, a, b = t[i], t[i + 1], T[i], T[i + 1]

        def f(x, t0=t0, t1=t1, a=a, b=b):
            return math.exp(ln_r * (T_REF - (a + (b - a) * (x - t0) / (t1 - t0))))

        val, _err = quad(f, t0, t1, epsabs=0.0, epsrel=1e-13, limit=200)
        tot += val
    return tot


def _rect(t, T):
    return float(np.sum(R_BELOW ** (T_REF - T[:-1]) * np.diff(t)))


def test_q_dose_equals_quadrature(series):
    t = np.asarray(series["times_s"], dtype=float)
    th = np.asarray(series["theta_ball10g_k"], dtype=float)
    T = series["t_inf_c"] + th
    assert float(np.max(T)) < T_REF, float(np.max(T))
    dose = cem43(t, T)
    ref = _quad_dose(t, T, math.log(R_BELOW))
    rel = abs(dose - ref) / abs(ref)
    swap = cem43(t, T, r_above=R_BELOW, r_below=0.5)
    miss = abs(swap - ref) / abs(ref)
    rect = _rect(t, T)
    _p(f"[MAT-20 s3] (Q) ball-mean CEM43 {dose:.15e} s vs per-segment quad {ref:.15e} s: rel "
       f"{rel:.3e} (bar {DOSE_BAR:g}, ASSERTED; predicted <= {PRED_Q:g}); max T {np.max(T):.9f} C "
       f"< 43 ASSERTED")
    _p(f"[MAT-20 s3] (Rc') R swap {swap:.6e} s: miss {miss:.4e} relative (asserted > "
       f"{DOSE_BAR:g}; predicted factor ~{PRED_RC_FACTOR:.0e}, printed)")
    _p(f"[MAT-20 s3] (Rect) rectangle sum {rect:.15e} s, rel to (Q) {rect / ref - 1.0:+.4e} "
       f"(printed)")
    assert rel <= DOSE_BAR, f"(Q) {rel:.3e}"
    assert miss > DOSE_BAR, f"(Rc') R swap did not miss: {miss:.3e}"


def test_win_six_minute_max_inside_window(series):
    t = np.asarray(series["times_s"], dtype=float)
    th = np.asarray(series["theta_ball10g_k"], dtype=float)
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
    dt = float(series["dt_s"])
    _m2, mx_sh, t_sh = running_average(t, th, window=W - dt)
    _p(f"[MAT-20 s3] (Win) W {W:g} s ({W / dt:.4f} samples): max six-minute mean {mx:.15e} K "
       f"at t_max {t_mx:.6f} s (predicted t* {t[-1]:.3f}); (i) [{t[0] + W:.3f}, {t[-1]:.3f}] "
       f"ASSERTED; (ii) interpolant on window [{lo:.12e}, {hi:.12e}] ASSERTED; (iii) quad "
       f"{ref:.15e}, rel {rel:.3e} (bar {AVG_BAR:g}, ASSERTED; predicted <= {PRED_WIN_III:g})")
    _p(f"[MAT-20 s3] (Sh) window W - dt = {W - dt:.6f} s: max {mx_sh:.12e} K at t {t_sh:.6f} s "
       f"(printed)")
    assert t[0] + W <= t_mx <= t[-1], t_mx
    assert lo <= mx <= hi, (lo, mx, hi)
    assert rel <= AVG_BAR, f"(Win)(iii) {rel:.3e}"


def test_readings_printed(series):
    """Printed readings (no limit, no compliance); asserts only that they were formed."""
    t = np.asarray(series["times_s"], dtype=float)
    p_acc = series["p_acc_w"]["ccw"]
    scales = {
        "per accepted W (361-1)": 1.0,
        "per W of 1/2 Re sum V I*": p_acc / series["p_vi_w"]["ccw"],
        "per W of P_src": p_acc / series["p_src_w"],
    }
    _p(f"[MAT-20 s3] P_acc,ccw {p_acc:.12e} W, P_VI,ccw {series['p_vi_w']['ccw']:.12e} W, "
       f"P_src {series['p_src_w']:.12e} W; P_acc/P_VI {scales['per W of 1/2 Re sum V I*']:.6f} "
       f"(predicted ~{PRED_PVI_SCALE}), P_acc/P_src {scales['per W of P_src']:.6f}")
    formed = 0
    for key, label in (("theta_ball10g_k", "10 g ball-mean"), ("theta_mean_k", "volume mean"),
                       ("theta_probe_k", "5 mm probe")):
        th = np.asarray(series[key], dtype=float)
        for norm, sc in scales.items():
            _m, mx, t_mx = running_average(t, sc * th, window=float(W_B2))
            dose = cem43(t, series["t_inf_c"] + sc * th)
            _p(f"[MAT-20 s3] READING {label:<14s} {norm:<26s}: six-minute peak theta {mx:.6e} K "
               f"at {t_mx:.1f} s; theta(t*) {sc * th[-1]:.6e} K; CEM43 at P* = 1 W {dose:.6e} s = "
               f"{dose / 60.0:.6e} min (predicted ball ~{PRED_DOSE_S:g} s; not a safety "
               f"statement)")
            formed += int(np.isfinite(mx) and np.isfinite(dose))
    assert formed == 9
