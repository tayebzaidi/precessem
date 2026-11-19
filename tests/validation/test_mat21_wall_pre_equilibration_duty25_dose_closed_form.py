"""`MAT-21` step 35 — the wait in a dose inside the walled box under D = 0.25, closed form
(§9 item 789, §10 P64; brief `docs/planning/chunks/MAT-21.md` "§9 item 789").

Numpy.  Sets no band and moves none.  P59′'s uniform-sum dose
(`test_mat21_pre_equilibration_duty_dose_closed_form.py`) with the wall in both parts:
θ(t; t_w) = Σ_{t > t_m} (−1)^m θ_wall(t − t_m) + θ_drift(t + t_w), where θ_wall is the
two-layer even series' step rise (item 788's ``_sourced``; here vectorised over the three
heights with its eigenvalues cached, asserted equal to ``_sourced``) and θ_drift is item 788's
h_eff drift (the walled series has even modes only; the wall's heat capacity is not in the
drift).  Schedule ``T_SW25`` (D = 0.25) on the box's ``T_GRID``; waits {0, ``T_C``}; points
P5 (z = C), ``ADJ_PT`` (0.9 C), P0 (z = 0).

Asserted: (Q-d) ``cem43`` against its quadrature at ``DOSE_BAR`` with the R-swap control and
(Win-d) the running average at ``AVG_BAR``, per point and wait; (Sym-d) at P0 the dose equals
the drift-free dose at ``DOSE_BAR``; (D50-cf) the walled uniform sum at ``T_SW`` (D = 0.5)
against ``two_layer_slab_duty_cycle_temperature`` (P = 2 ``HALF_T``, D = 0.5) in dose at
``DOSE_BAR``, with the drop-run-0 control > ``DOSE_BAR``.  The duty series is read with
``side="left"``: the sum's strict ``t > t_m`` excludes the run that starts at a grid switch,
which is the left limit (``side="right"`` adds the N = 400 truncation residual there).
Control (asserted): the massless-wall sourced rise (``rho_c_w = MASSLESS * RHO_C_W``) misses
the walled dose at z = C by > ``DOSE_BAR`` at each wait.  Printed (🟡 record, no band): the
dose with and without the wait per height, their ratio, beside P59′'s wall-free values.

Not claimed: the wall's heat capacity in the drift, any FEM (P65), a recommended wait.
"""

from __future__ import annotations

import numpy as np

from fem_em_solver.utils.analytical import (
    two_layer_slab_duty_cycle_temperature,
    two_layer_slab_eigenvalues,
    two_layer_slab_transient_temperature,
)
from tests.validation.test_mat15_finite_cylinder_closed_form import H, Q
from tests.validation.test_mat21_box_closed_form import T_C
from tests.validation.test_mat21_box_dose_closed_form import AVG_BAR, DOSE_BAR, T_GRID, W_B2
from tests.validation.test_mat21_box_dt_duty125_dose_closed_form import _q_check_g, _win_check_g
from tests.validation.test_mat21_box_duty25_dose_closed_form import DUTY25, T_SW25
from tests.validation.test_mat21_box_duty_closed_form import DT_BOX, DUTY, HALF_T, T_SW
from tests.validation.test_mat21_pre_equilibration_duty_dose_closed_form import (
    LABS,
    PTS,
    T_INF,
)
from tests.validation.test_mat21_pre_equilibration_duty_dose_closed_form import (
    _hists as _hists_free,
)
from tests.validation.test_mat21_wall_closed_form import RHO_C_W
from tests.validation.test_mat21_wall_pre_equilibration_closed_form import H_EFF, _drift, _sourced
from tests.validation.test_mat21_wall_transient_closed_form import N_SERIES, WALL
from tests.validation.test_mat21_wall_transient_fem import MASSLESS

WAITS = (("0", 0.0), ("T_C", T_C))
Z = np.array([p[2] for p in PTS], dtype=float)
PERIOD_BOX = 2.0 * HALF_T
_LAM = {}
_CACHE = {}


def _say(msg):
    print(msg, flush=True)


def _wall(rho_c_w):
    return dict(WALL, rho_c_w=rho_c_w)


def _lam(rho_c_w=RHO_C_W):
    if rho_c_w not in _LAM:
        _LAM[rho_c_w] = two_layer_slab_eigenvalues(h=H, n=N_SERIES, **_wall(rho_c_w))
    return _LAM[rho_c_w]


def _step(tau, rho_c_w=RHO_C_W):
    """The two-layer step rise θ (K) at the three heights ``Z``."""
    return np.asarray(
        two_layer_slab_transient_temperature(
            Z,
            tau,
            q=Q,
            h=H,
            t_inf=0.0,
            n_terms=N_SERIES,
            eigenvalues=_lam(rho_c_w),
            **_wall(rho_c_w),
        ),
        dtype=float,
    )


def _uniform_w(sw, rho_c_w=RHO_C_W, drop=()):
    """(41, 3) walled uniform sum Σ_{t > t_m} (−1)^m θ_wall(t − t_m); row 0 = 0.0."""
    out = np.zeros((T_GRID.size, Z.size))
    for n, t in enumerate(T_GRID):
        for m, tm in enumerate(sw):
            if t > tm and m not in drop:
                out[n] += (-1.0) ** m * _step(t - tm, rho_c_w)
    return out


def _duty_series(duty):
    out = np.zeros((T_GRID.size, Z.size))
    for n, t in enumerate(T_GRID):
        if n:
            out[n] = two_layer_slab_duty_cycle_temperature(
                Z,
                t,
                q=Q,
                h=H,
                period=PERIOD_BOX,
                duty=duty,
                t_inf=0.0,
                n_terms=N_SERIES,
                eigenvalues=_lam(),
                side="left",
                **WALL,
            )
    return out


def _dose(th):
    return _q_check_g(T_INF + th, T_GRID)[0]


def _hists():
    if not _CACHE:
        _CACHE["25"] = _uniform_w(T_SW25)
        _CACHE["50"] = _uniform_w(T_SW)
        _CACHE["50drop"] = _uniform_w(T_SW, drop=(0,))
        _CACHE["25m"] = _uniform_w(T_SW25, rho_c_w=MASSLESS * RHO_C_W)
        _CACHE["drift"] = {
            (j, wl): np.array([_drift(p, t + w) for t in T_GRID])
            for j, p in enumerate(PTS)
            for wl, w in WAITS
        }
    return _CACHE


def test_preconditions():
    assert all(min(abs(T_GRID - t)) <= 1e-9 * DT_BOX for t in T_SW25)
    assert all(min(abs(T_GRID - t)) <= 1e-9 * DT_BOX for t in T_SW)
    assert abs(PERIOD_BOX / DT_BOX - 8.0) < 1e-12 and DUTY == 0.5 and DUTY25 == 0.25
    assert abs(T_SW25[1] - DUTY25 * PERIOD_BOX) < 1e-9 * DT_BOX
    worst = 0.0
    for tau in (0.5 * DT_BOX, 7.0 * DT_BOX, T_GRID[-1]):
        v = _step(tau)
        for j, p in enumerate(PTS):
            worst = max(worst, abs(v[j] - _sourced(p, tau)) / abs(_sourced(p, tau)))
    _say(
        f"MAT21-35 h_eff {H_EFF:.6f}; {T_GRID.size} rows, DT {DT_BOX:.4f} s, P {PERIOD_BOX:.4f} s; "
        f"{len(T_SW25)} + {len(T_SW)} switches on the grid; vectorised step vs item 788 _sourced "
        f"max rel {worst:.3e} (asserted <= 1e-10, item 788 adds and removes T_INF); bars DOSE "
        f"{DOSE_BAR:g}, AVG {AVG_BAR:g}"
    )
    assert worst <= 1e-10, worst
    assert np.all(_hists()["25"][0] == 0.0)


def test_q_d_win_d_sym_d_and_rc():
    h = _hists()
    d25 = h["25"]
    worst = {"q": 0.0, "w": 0.0, "rc": np.inf, "sym": 0.0}
    for j, lab in enumerate(LABS):
        d_free = _dose(d25[:, j])
        for wl, _w in WAITS:
            th = d25[:, j] + h["drift"][(j, wl)]
            dose, _quad, rel, _swap, miss = _q_check_g(T_INF + th, T_GRID)
            mx, t_mx, lo, hi, _ref, wrel = _win_check_g(th, T_GRID)
            assert rel <= DOSE_BAR, (lab, wl, rel)
            assert miss > DOSE_BAR, (lab, wl, miss)
            assert lo <= mx <= hi and T_GRID[0] + W_B2 <= t_mx <= T_GRID[-1], (lab, wl, t_mx)
            assert wrel <= AVG_BAR, (lab, wl, wrel)
            worst["q"] = max(worst["q"], rel)
            worst["w"] = max(worst["w"], wrel)
            worst["rc"] = min(worst["rc"], miss)
            if j == 2:
                sym = abs(dose / d_free - 1.0)
                worst["sym"] = max(worst["sym"], sym)
                assert sym <= DOSE_BAR, (wl, sym)
    _say(
        f"MAT21-35 (Q-d) worst rel {worst['q']:.3e} (bar {DOSE_BAR:g}, ASSERTED); (Win-d) worst "
        f"rel {worst['w']:.3e} (bar {AVG_BAR:g}, ASSERTED); (Sym-d) P0 worst rel {worst['sym']:.3e} "
        f"(ASSERTED); R swap least miss {worst['rc']:.4e} (asserted > {DOSE_BAR:g})"
    )


def test_d50_cf_against_the_two_layer_duty_series():
    h = _hists()
    cf50 = _duty_series(DUTY)
    cf25 = _duty_series(DUTY25)
    worst, least, w25 = 0.0, np.inf, 0.0
    for j, lab in enumerate(LABS):
        d_sum = _dose(h["50"][:, j])
        d_cf = _dose(cf50[:, j])
        d_drop = _dose(h["50drop"][:, j])
        rel = abs(d_sum / d_cf - 1.0)
        miss = abs(d_drop / d_cf - 1.0)
        r25 = abs(_dose(h["25"][:, j]) / _dose(cf25[:, j]) - 1.0)
        worst, least, w25 = max(worst, rel), min(least, miss), max(w25, r25)
        _say(
            f"MAT21-35 (D50-cf) {lab}: walled sum dose {d_sum:.6e} s, two-layer duty series dose "
            f"{d_cf:.6e} s, rel {rel:.3e} (bar {DOSE_BAR:g}, ASSERTED); drop=(0,) miss {miss:.4e} "
            f"(asserted > {DOSE_BAR:g}); D = 0.25 sum vs series rel {r25:.3e} (printed)"
        )
        assert rel <= DOSE_BAR, (lab, rel)
        assert miss > DOSE_BAR, (lab, miss)
    _say(
        f"MAT21-35 (D50-cf) worst {worst:.3e}; drop=(0,) least miss {least:.4e}; D = 0.25 worst "
        f"{w25:.3e} (printed, no band)"
    )


def test_control_massless_wall_misses_walled_dose_at_z_c():
    h = _hists()
    j = 0  # P5, z = C
    least = np.inf
    for wl, _w in WAITS:
        dr = h["drift"][(j, wl)]
        d_w = _dose(h["25"][:, j] + dr)
        d_m = _dose(h["25m"][:, j] + dr)
        miss = abs(d_m / d_w - 1.0)
        least = min(least, miss)
        _say(
            f"MAT21-35 (Mw) {LABS[j]} t_w {wl}: walled dose {d_w:.9e} s, massless-wall dose "
            f"{d_m:.9e} s, rel miss {miss:.4e} (asserted > {DOSE_BAR:g})"
        )
        assert miss > DOSE_BAR, (wl, miss)
    _say(f"MAT21-35 (Mw) least miss {least:.4e}")


def test_table_printed():
    """The record (no band); asserts only that it was formed."""
    h = _hists()
    _q, uni_free, drift_free = _hists_free()
    _say(
        "MAT21-35 TABLE point | D_ref (no drift) s | D(t_w=0) s | D(t_w=T_C) s | D(T_C)/D(0) | "
        "D(0)/D_ref-1 | D(T_C)/D_ref-1 || wall-free (P59') D_ref | D(0)/D_ref-1 | D(T_C)/D_ref-1"
    )
    for j, lab in enumerate(LABS):
        ref = _dose(h["25"][:, j])
        d0 = _dose(h["25"][:, j] + h["drift"][(j, "0")])
        dc = _dose(h["25"][:, j] + h["drift"][(j, "T_C")])
        fref = _dose(uni_free["0.25"][:, j])
        f0 = _dose(uni_free["0.25"][:, j] + drift_free[(j, "0")])
        fc = _dose(uni_free["0.25"][:, j] + drift_free[(j, "T_C")])
        _say(
            f"MAT21-35 TABLE {lab} | {ref:.6e} | {d0:.6e} | {dc:.6e} | {dc / d0:.6f} | "
            f"{d0 / ref - 1.0:+.6e} | {dc / ref - 1.0:+.6e} || {fref:.6e} | {f0 / fref - 1.0:+.6e} | "
            f"{fc / fref - 1.0:+.6e}"
        )
        assert ref > 0.0
