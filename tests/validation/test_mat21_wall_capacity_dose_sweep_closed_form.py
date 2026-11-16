"""`MAT-21` step 41 — the dose under the capacitive wall as a function of the wait, closed form
(§9 item 875, §10 P73; brief `docs/planning/chunks/MAT-21.md` "§9 item 875").

Numpy.  Sets no band and moves none.  P70's D = 0.25 dose — P64's walled uniform sum plus item
823's composite drift ``_drift_cap_z(Z, t + t_w)`` — at the four waits
``W41 = WAITS[:4]`` = {0, 0.5, 1, 2} / ``LAM1_W`` (P67's own wait set; ``LAM1_W`` the composite's
slowest odd rate), at P5 (z = C), ``ADJ_PT`` (0.9 C) and P0 (z = 0, the odd drift's node).

Asserted (all collected, then asserted — item 668): (Lim) at t_w = 0 against P70's own
construction at ``B1_BAND`` and the reproduction of P70's printed digits (half the last digit);
(Bq) N_TERMS vs 2·N_TERMS on the dose at t_w > 0 at ``B1_BAND``; (Mono) strict monotonicity in
the wait at P5 / ``ADJ_PT`` with each |ΔD|/D > ``DOSE_BAR``, constancy at P0 within
``DOSE_BAR``; MAT-20's (Q-d) with the R-swap control and (Win-d) per point and wait; the
side table's d = ``D_WALL`` row against the main t_w = 0 row at ``B1_BAND`` (precondition).
Control (asserted, factor 1 at ``DOSE_BAR``): at P5, every wait, the massless-wall
(``_drift_1d``) dose misses the composite dose by > ``DOSE_BAR``.
Printed (🟡 record, no band): per point × wait the composite drift at t* + t_w, D, D/D_ref − 1,
D/D(0), the massless D and composite/massless − 1; (Bq) at t_w = 0 (row 0 is the unconverged IC
series); the thickness side table at t_w = 0 for d ∈ {0.5, 1, 2}·``D_WALL``.

Not claimed: a dose band, any FEM (P74), a recommended wait or thickness, a cited PMMA ρc.
"""

from __future__ import annotations

import numpy as np

from fem_em_solver.utils.analytical import (
    two_layer_slab_eigenvalues,
    two_layer_slab_transient_temperature,
)

from tests.validation.test_imp4_pec_cylinder_closed_form import B1_BAND
from tests.validation.test_mat15_finite_cylinder_closed_form import H, Q
from tests.validation.test_mat21_box_dose_closed_form import AVG_BAR, DOSE_BAR, T_GRID, W_B2
from tests.validation.test_mat21_box_dt_duty125_dose_closed_form import _q_check_g, _win_check_g
from tests.validation.test_mat21_box_duty25_dose_closed_form import T_SW25
from tests.validation.test_mat21_initial_gradient_closed_form import N_TERMS
from tests.validation.test_mat21_pre_equilibration_duty_dose_closed_form import LABS, PTS, T_INF
from tests.validation.test_mat21_wall_capacity_drift_closed_form import (
    LAM1_W,
    WAIT_LABELS,
    WAITS,
    _drift_cap_z,
    _odd_eigs,
)
from tests.validation.test_mat21_wall_capacity_duty25_dose_closed_form import _cap_hists
from tests.validation.test_mat21_wall_closed_form import D_WALL
from tests.validation.test_mat21_wall_pre_equilibration_closed_form import (
    T_STAR,
    WAIT_C,
    _drift_1d,
)
from tests.validation.test_mat21_wall_pre_equilibration_duty25_dose_closed_form import (
    Z,
    _dose,
    _hists,
    _step,
    _uniform_w,
)
from tests.validation.test_mat21_wall_transient_closed_form import N_SERIES, WALL

W41 = tuple(WAITS[:4])
W41_LABELS = tuple(WAIT_LABELS[:4])
D_FACTORS = (0.5, 1.0, 2.0)

# P70's printed record (`20261105T224306Z_MAT-21-P70-W1.log`), each bound half the last digit.
P70_D0 = {0: (7.219294348e-11, 5e-21),   # :52, P5 composite-drift dose
          1: (7.360929e-11, 5e-18),      # :57, ADJ_PT D(t_w=0)
          2: (5.280312e-11, 5e-18)}      # :58, P0 D(t_w=0)
P70_RATIO_P5 = (3.7292e-02, 5e-7)        # :52, composite/massless - 1 at P5, t_w = 0
P70_RATIO_TC = 3.7934e-03                # :53, printed beside only
P70_BQ0 = 1.594e-05                      # :62, printed beside only
P70_CEIL_P5 = 4.546295e-01               # :56, D(0)/D_ref - 1 at P5, the wait's ceiling
P67_DCAP_P5 = (2.082472e-02, 1.586026e-02, 1.163907e-02, 5.697480e-03)  # P67-W1.log:68–71

_CAPW = {}
_LAMD = {}
_UNID = {}
_DRD = {}


def _say(msg):
    print(msg, flush=True)


def _cap_hist_w(w, n=N_TERMS):
    """(41, 3) composite drift on T_GRID at Z, at wait w; cached per (w, n), never per point."""
    key = (float(w), int(n))
    if key not in _CAPW:
        _CAPW[key] = np.array([_drift_cap_z(Z, t + w, n=n) for t in T_GRID])
    return _CAPW[key]


def _d_comp(j, w, n=N_TERMS):
    return _dose(_hists()["25"][:, j] + _cap_hist_w(w, n)[:, j])


def _d_mass(j, w):
    """Massless-wall comparand: P70's `_z1d_hists` route at wait w."""
    p = PTS[j]
    return _dose(_hists()["25"][:, j] + np.array([_drift_1d(p, t + w) for t in T_GRID]))


# ------------------------------------------------------------------ thickness side table
def _lam_d(d):
    if d not in _LAMD:
        _LAMD[d] = two_layer_slab_eigenvalues(h=H, n=N_SERIES, **dict(WALL, wall=d))
    return _LAMD[d]


def _step_d(tau, d):
    return np.asarray(two_layer_slab_transient_temperature(
        Z, tau, q=Q, h=H, t_inf=0.0, n_terms=N_SERIES, eigenvalues=_lam_d(d),
        **dict(WALL, wall=d)), dtype=float)


def _uniform_wd(sw, d):
    """Copy of P64's `_uniform_w` with the wall thickness d; cached per d."""
    if d not in _UNID:
        out = np.zeros((T_GRID.size, Z.size))
        for n, t in enumerate(T_GRID):
            for m, tm in enumerate(sw):
                if t > tm:
                    out[n] += (-1.0) ** m * _step_d(t - tm, d)
        _UNID[d] = out
    return _UNID[d]


def _drift_d(d):
    if d not in _DRD:
        _DRD[d] = np.array([_drift_cap_z(Z, t, d=d) for t in T_GRID])
    return _DRD[d]


def _side_row(d):
    """(D at the three points, lambda_1 odd) at t_w = 0 for wall thickness d."""
    u = _uniform_wd(T_SW25, d)
    dr = _drift_d(d)
    return [_dose(u[:, j] + dr[:, j]) for j in range(Z.size)], float(_odd_eigs(1, d)[0])


# ------------------------------------------------------------------ tests
def test_preconditions_waits_and_side_table_anchor():
    gates = []
    gates.append((("W41[0]",), W41[0] == 0.0, W41[0]))
    for w, c in zip(W41, WAIT_C):
        gates.append((("w*LAM1_W", c), abs(w * LAM1_W - c) <= 1e-12, w * LAM1_W - c))
    # the local copy is P64's construction: step and uniform sum at d = D_WALL
    worst_step = 0.0
    for tau in (0.5 * (T_GRID[1] - T_GRID[0]), T_GRID[7], T_GRID[-1]):
        a, b = _step_d(tau, D_WALL), _step(tau)
        worst_step = max(worst_step, float(np.max(np.abs(a - b) / np.abs(b))))
    u_rel = float(np.max(np.abs(_uniform_wd(T_SW25, D_WALL) - _uniform_w(T_SW25))))
    row, _l1 = _side_row(D_WALL)
    side = max(abs(row[j] / _d_comp(j, 0.0) - 1.0) for j in range(Z.size))
    _say(f"MAT21-41 W41 {[round(w, 2) for w in W41]} s, labels {list(W41_LABELS)}; LAM1_W "
         f"{LAM1_W:.6e} /s; local step(d=D_WALL) vs P64 _step max rel {worst_step:.3e}; uniform "
         f"sum max abs diff {u_rel:.3e} K; side d = D_WALL row vs main t_w = 0 dose max rel "
         f"{side:.3e} (band {B1_BAND:g}, ASSERTED)")
    gates.append((("side-anchor",), side <= B1_BAND, side))
    bad = [g for g in gates if not g[1]]
    assert not bad, bad


def test_lim_reproduction_bq_mono():
    gates = []
    cap70 = _cap_hists()["0"]
    d25 = _hists()["25"]
    for j, lab in enumerate(LABS):
        d70 = _dose(d25[:, j] + cap70[:, j])
        d0 = _d_comp(j, 0.0)
        rel = abs(d0 - d70) / d70
        val, tol = P70_D0[j]
        _say(f"MAT21-41 (Lim) {lab}: D(t_w=0) {d0:.12e} s vs P70 construction {d70:.12e} s, rel "
             f"{rel:.3e} (band {B1_BAND:g}); vs P70 printed {val:.9e} |diff| {abs(d0 - val):.3e} s "
             f"(bound {tol:g})")
        gates.append(((lab, "Lim"), rel <= B1_BAND, rel))
        gates.append(((lab, "repro"), abs(d0 - val) <= tol, d0 - val))
    r_p5 = _d_comp(0, 0.0) / _d_mass(0, 0.0) - 1.0
    _say(f"MAT21-41 (Lim) P5 composite/massless-1 at t_w = 0 {r_p5:+.7e} vs P70 "
         f"{P70_RATIO_P5[0]:+.4e} (bound {P70_RATIO_P5[1]:g})")
    gates.append((("P5", "repro-ratio"), abs(r_p5 - P70_RATIO_P5[0]) <= P70_RATIO_P5[1], r_p5))
    # (Bq)
    for k, (w, wl) in enumerate(zip(W41, W41_LABELS)):
        bq = max(abs(_d_comp(j, w) / _d_comp(j, w, 2 * N_TERMS) - 1.0) for j in range(Z.size))
        if k == 0:
            _say(f"MAT21-41 (Bq) t_w {wl}: dose N vs 2N worst rel {bq:.3e} (printed; row 0 is the "
                 f"unconverged IC series; P70 read {P70_BQ0:.3e}, P70-W1.log:62)")
        else:
            _say(f"MAT21-41 (Bq) t_w {wl}: dose N vs 2N worst rel {bq:.3e} (band {B1_BAND:g}, "
                 f"ASSERTED)")
            gates.append(((wl, "Bq"), bq <= B1_BAND, bq))
    # (Mono)
    for j, lab in enumerate(LABS):
        ds = np.array([_d_comp(j, w) for w in W41])
        if j < 2:
            diff = np.diff(ds)
            sgn = np.sign(diff)
            sep = np.abs(diff) / ds[:-1]
            ok = bool(np.all(sgn == sgn[0]) and sgn[0] != 0 and np.all(sep > DOSE_BAR))
            _say(f"MAT21-41 (Mono) {lab}: doses {', '.join(f'{d:.9e}' for d in ds)} s; dD signs "
                 f"{sgn.tolist()}; |dD|/D {', '.join(f'{s:.4e}' for s in sep)} (asserted > "
                 f"{DOSE_BAR:g}, one sign); monotone sign {'decreasing' if sgn[0] < 0 else 'increasing'}")
            gates.append(((lab, "Mono"), ok, sep.min()))
        else:
            dev = float(np.max(np.abs(ds / ds[0] - 1.0)))
            _say(f"MAT21-41 (Mono) {lab}: doses {', '.join(f'{d:.9e}' for d in ds)} s; max "
                 f"|D/D(0)-1| {dev:.3e} (asserted <= {DOSE_BAR:g})")
            gates.append(((lab, "Const"), dev <= DOSE_BAR, dev))
    bad = [g for g in gates if not g[1]]
    assert not bad, bad


def test_q_d_win_d_and_rc_on_the_new_waits():
    d25 = _hists()["25"]
    worst = {"q": 0.0, "w": 0.0, "rc": np.inf}
    gates = []
    for j, lab in enumerate(LABS):
        for w, wl in zip(W41, W41_LABELS):
            th = d25[:, j] + _cap_hist_w(w)[:, j]
            _dose_v, _quad, rel, _swap, miss = _q_check_g(T_INF + th, T_GRID)
            mx, t_mx, lo, hi, _ref, wrel = _win_check_g(th, T_GRID)
            gates.append(((lab, wl, "Q"), rel <= DOSE_BAR, rel))
            gates.append(((lab, wl, "swap"), miss > DOSE_BAR, miss))
            gates.append(((lab, wl, "win-loc"),
                          lo <= mx <= hi and T_GRID[0] + W_B2 <= t_mx <= T_GRID[-1], t_mx))
            gates.append(((lab, wl, "Win"), wrel <= AVG_BAR, wrel))
            worst["q"] = max(worst["q"], rel)
            worst["w"] = max(worst["w"], wrel)
            worst["rc"] = min(worst["rc"], miss)
    _say(f"MAT21-41 (Q-d) worst rel {worst['q']:.3e} (bar {DOSE_BAR:g}, ASSERTED); (Win-d) worst "
         f"rel {worst['w']:.3e} (bar {AVG_BAR:g}, ASSERTED); R swap least miss "
         f"{worst['rc']:.4e} (asserted > {DOSE_BAR:g}); 3 points x {len(W41)} waits")
    bad = [g for g in gates if not g[1]]
    assert not bad, bad


def test_control_massless_wall_misses_composite_at_p5_every_wait():
    gates = []
    for k, (w, wl) in enumerate(zip(W41, W41_LABELS)):
        d_c, d_m = _d_comp(0, w), _d_mass(0, w)
        miss = abs(d_m / d_c - 1.0)
        _say(f"MAT21-41 (Ctrl) P5 t_w {wl}: composite {d_c:.9e} s, massless {d_m:.9e} s, rel miss "
             f"{miss:.4e} (ASSERTED > {DOSE_BAR:g}); composite/massless-1 {d_c / d_m - 1.0:+.4e} "
             f"(predicted to fall from {P70_RATIO_P5[0]:+.4e} toward {P70_RATIO_TC:+.4e} at T_C; "
             f"P67 Delta_cap at t*+t_w {P67_DCAP_P5[k]:.6e} K, printed)")
        gates.append(((wl, "Ctrl"), miss > DOSE_BAR, miss))
    bad = [g for g in gates if not g[1]]
    assert not bad, bad


def test_record_printed():
    """The record (no band); asserts only that it was formed."""
    d25 = _hists()["25"]
    _say("MAT21-41 TABLE point | t_w | t_w s | drift at t*+t_w K | D(t_w) s | D/D_ref-1 | "
         "D/D(0) | massless D s | composite/massless-1")
    for j, lab in enumerate(LABS):
        ref = _dose(d25[:, j])
        d0 = _d_comp(j, 0.0)
        signs = []
        prev = None
        for w, wl in zip(W41, W41_LABELS):
            dr = float(_drift_cap_z(Z, T_STAR + w)[j])
            d = _d_comp(j, w)
            dm = _d_mass(j, w)
            _say(f"MAT21-41 TABLE {lab} | {wl} | {w:.2f} | {dr:+.6e} | {d:.6e} | "
                 f"{d / ref - 1.0:+.6e} | {d / d0:.6f} | {dm:.6e} | {d / dm - 1.0:+.4e}")
            if prev is not None:
                signs.append(int(np.sign(d - prev)))
            prev = d
        _say(f"MAT21-41 TABLE {lab} D_ref {ref:.6e} s; signs of successive dD {signs}"
             + (f"; P70 composite/massless-1 at t_w=0 {P70_RATIO_P5[0]:+.4e}, ceiling "
                f"D(0)/D_ref-1 {P70_CEIL_P5:+.6e} (P70-W1.log:56)" if j == 0 else ""))
        assert ref > 0.0
    _say("MAT21-41 SIDE (thickness, t_w = 0) d/D_WALL | d m | lambda1_odd /s | D P5 | D ADJ_PT | "
         "D P0 | D/D(D_WALL)-1 P5 | ADJ_PT | P0")
    base, _l = _side_row(D_WALL)
    for f in D_FACTORS:
        d = f * D_WALL
        try:
            row, l1 = _side_row(d)
        except Exception as exc:  # the table never gates
            _say(f"MAT21-41 SIDE {f:g} | {d:.4f} | REFUSED: {exc!r} (row skipped)")
            continue
        _say(f"MAT21-41 SIDE {f:g} | {d:.4f} | {l1:.6e} | {row[0]:.6e} | {row[1]:.6e} | "
             f"{row[2]:.6e} | {row[0] / base[0] - 1.0:+.6e} | {row[1] / base[1] - 1.0:+.6e} | "
             f"{row[2] / base[2] - 1.0:+.6e}")
