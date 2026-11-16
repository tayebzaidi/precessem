"""`MAT-21` step 28 — the pre-equilibration bias in dose, in closed form (§9 item 727, §10
chain P, P55, provisional under rule (6); brief `docs/planning/chunks/MAT-21.md` "§9 item
727").

Numpy.  Sets no band and moves none.  On P51's θ(t; t_w) = ``_combined(p, t, t_w)``, the
history T∞ + θ on 161 rows over [0, t\\*], at P5, ``ADJ_PT`` and P0, on P51's four waits.  g
= ``G_GRAD`` and q are the fixture's choices, so this is a sensitivity at a 20 °C gel's
scale: an identity on a history.

Asserted: (Q-w) ``cem43`` against its quadrature at ``DOSE_BAR`` and (Win-w) the running
average against its identity at ``AVG_BAR``, per point and wait; (Sym-D) at P0 every wait's
dose equals the sourced-only dose at ``DOSE_BAR``; (Lim-D) at t_w = 40 / λ₁ the dose at P5
does; (Jen) at P5, t_w = 0, D(+g) + D(−g) > 2 D(g = 0).  Negative controls: the Dirichlet-z
drift misses D at P5 by > ``DOSE_BAR``; the R swap misses by > ``DOSE_BAR``.  Printed (**no
band**): b(t_w) beside D(t_w) / D_ref − 1, and the grid-halved dose.

Not claimed: a recommended wait, F2182's wording, a dose accuracy, a thermal-damage
statement.
"""

from __future__ import annotations

import numpy as np

from tests.validation.test_mat21_box_closed_form import KW
from tests.validation.test_mat21_box_dose_closed_form import AVG_BAR, DOSE_BAR, W_B2
from tests.validation.test_mat21_box_dt_duty125_dose_closed_form import (
    _q_check_g,
    _trap_g,
    _win_check_g,
)
from tests.validation.test_mat21_box_probe_points_fem import PROBE_LABELS, PROBE_PTS
from tests.validation.test_mat21_initial_gradient_closed_form import G_GRAD, N_TERMS, T_STAR
from tests.validation.test_mat21_pre_equilibration_closed_form import (
    ADJ_PT,
    LAM_ODD,
    P5,
    T_INF,
    WAITS,
    _combined,
    _drift,
    _sourced,
)

N_ROWS = 161
GRID = np.linspace(0.0, T_STAR, N_ROWS)
P0 = PROBE_PTS[0]
LIM_WAIT = 40.0 / LAM_ODD
READ_PTS = ((PROBE_LABELS[5], P5), ("ADJ_PT", ADJ_PT), (PROBE_LABELS[0], P0))
WAIT_COEF = (0.0, 0.5, 1.0, 2.0)
DZ_ROOTS = np.arange(1, N_TERMS + 1) * np.pi  # P51's control roots


def _say(msg):
    print(msg, flush=True)


def _dose(th, t=GRID):
    return _q_check_g(T_INF + th, t)[0]


def _hist(p, t_w, g=G_GRAD):
    return np.array([_combined(p, t, t_w, kw=KW, g=g) for t in GRID])


def _ref(p):
    """The sourced-only history; no drift term at all."""
    return np.array([_sourced(p, t) for t in GRID])


def test_preconditions():
    assert GRID.size == N_ROWS and abs(GRID[-1] / T_STAR - 1.0) <= 1e-14
    assert W_B2 <= T_STAR, (W_B2, T_STAR)
    assert P0[2] == 0.0 and P5[2] > 0.0 and ADJ_PT[2] > 0.0
    _say(f"MAT21-28 grid: {N_ROWS} rows to t* {T_STAR:.2f} s, window {W_B2:g} s; waits "
         f"{[round(w, 2) for w in WAITS]} s, limit wait {LIM_WAIT:.1f} s; g {G_GRAD} K/m (fixture)")


def test_q_w_win_w_rc_table_and_sym_d():
    worst = {"q": 0.0, "w": 0.0, "rc": np.inf, "sym": 0.0}
    _say("MAT21-28 TABLE point | t_w lambda_1 | b(t_w) | D s | D_ref s | D/D_ref-1 | "
         "4^dtheta-1 | D halved rel")
    for lab, p in READ_PTS:
        ref = _ref(p)
        d_ref = _dose(ref)
        w4 = np.exp(np.log(4.0) * ref)
        for c, w in zip(WAIT_COEF, WAITS):
            th = _hist(p, w)
            dose, quad, rel, _swap, miss = _q_check_g(T_INF + th, GRID)
            mx, t_mx, lo, hi, _wref, wrel = _win_check_g(th, GRID)
            assert rel <= DOSE_BAR, (lab, c, rel)
            assert miss > DOSE_BAR, (lab, c, miss)
            assert GRID[0] + W_B2 <= t_mx <= GRID[-1], (lab, c, t_mx)
            assert lo <= mx <= hi, (lab, c, lo, mx, hi)
            assert wrel <= AVG_BAR, (lab, c, wrel)
            worst["q"] = max(worst["q"], rel)
            worst["w"] = max(worst["w"], wrel)
            worst["rc"] = min(worst["rc"], miss)
            s_t = _sourced(p, T_STAR)
            b = _drift(p, T_STAR + w) / s_t
            dth = _trap_g(w4 * (th - ref), GRID) / _trap_g(w4, GRID)
            halved = _dose(th[::2], GRID[::2])
            _say(f"MAT21-28 TABLE {lab} | {c} | {b:+.6e} | {dose:.6e} | {d_ref:.6e} | "
                 f"{dose / d_ref - 1.0:+.6e} | {4.0 ** dth - 1.0:+.6e} | "
                 f"{abs(halved / dose - 1.0):.3e}")
            if lab == PROBE_LABELS[0]:
                sym = abs(dose / d_ref - 1.0)
                worst["sym"] = max(worst["sym"], sym)
                assert sym <= DOSE_BAR, (c, sym)
    _say(f"MAT21-28 (Q-w) worst rel {worst['q']:.3e} (bar {DOSE_BAR:g}, ASSERTED); (Win-w) worst "
         f"rel {worst['w']:.3e} (bar {AVG_BAR:g}, ASSERTED); R swap least miss {worst['rc']:.4e} "
         f"(asserted > {DOSE_BAR:g}); (Sym-D) P0 worst rel {worst['sym']:.3e} (bar {DOSE_BAR:g}, "
         "ASSERTED)")


def test_lim_d_jensen_and_dirichlet_z():
    ref = _ref(P5)
    d_ref = _dose(ref)
    lim = abs(_dose(_hist(P5, LIM_WAIT)) / d_ref - 1.0)
    _say(f"MAT21-28 (Lim-D) P5 at t_w = 40/lambda_1: rel {lim:.3e} (bar {DOSE_BAR:g}, ASSERTED)")
    assert lim <= DOSE_BAR, lim

    d_plus = _dose(_hist(P5, 0.0))
    d_minus = _dose(_hist(P5, 0.0, g=-G_GRAD))
    d_zero = _dose(_hist(P5, 0.0, g=0.0))
    margin = (d_plus + d_minus) / (2.0 * d_zero) - 1.0
    _say(f"MAT21-28 (Jen) P5, t_w = 0: D(+g) {d_plus:.6e} s, D(-g) {d_minus:.6e} s, D(0) "
         f"{d_zero:.6e} s; (D+ + D-)/(2 D0) - 1 = {margin:+.4e} (asserted > 0)")
    assert d_plus + d_minus > 2.0 * d_zero, margin

    th_dz = ref + np.array([_drift(P5, t, roots_z=DZ_ROOTS) for t in GRID])
    d_dz = _dose(th_dz)
    miss = abs(d_dz / d_plus - 1.0)
    _say(f"MAT21-28 control (Dz) P5, t_w = 0: D Dirichlet {d_dz:.6e} s vs Robin {d_plus:.6e} s, "
         f"miss {miss:.4e} (asserted > {DOSE_BAR:g})")
    assert miss > DOSE_BAR, miss
