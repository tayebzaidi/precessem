"""`MAT-21` step 31 — the Dirichlet control at an interior point (§9 item 743, §10 chain P, P58,
provisional under rule (6); brief `docs/planning/chunks/MAT-21.md` "§9 item 743").

Numpy.  Sets no band and moves none.  Item 727's Dirichlet-z control reads P5 on z = C, where
that drift is zero by construction, so it shows a wrong boundary removing the bias at the face
and nothing about the interior.  This reads the same control at ``ADJ_PT`` on z = 0.9 C.

Asserted: (Dz-int) the Dirichlet history's dose at ``ADJ_PT`` misses D(t_w) by > ``DOSE_BAR`` at
t_w = 0; (Rep) item 727's face control reproduces 2.5379e-01 within 1e-5, the print's last
digit; (Q-w) ``cem43`` against its quadrature on the Dirichlet histories at ``DOSE_BAR``.
Negative control: the Dirichlet drift at P5 reads ≤ ``SYM0_BAR`` of g C at every grid row of
every wait, the stated reason the face control has no interior content, now measured.  Printed
(**no band**): the miss at ``ADJ_PT`` and the Robin bias against the wait, and the sign of the
Dirichlet-vs-Robin drift there.

Not claimed: anything about the walled box.
"""

from __future__ import annotations

import numpy as np

from tests.validation.test_mat21_box_closed_form import C_BOX
from tests.validation.test_mat21_box_dose_closed_form import DOSE_BAR
from tests.validation.test_mat21_box_dt_duty125_dose_closed_form import _q_check_g
from tests.validation.test_mat21_initial_gradient_closed_form import G_GRAD, T_STAR
from tests.validation.test_mat21_pre_equilibration_closed_form import (
    ADJ_PT,
    P5,
    SYM0_BAR,
    T_INF,
    WAITS,
    _drift,
)
from tests.validation.test_mat21_pre_equilibration_dose_closed_form import (
    DZ_ROOTS,
    GRID,
    WAIT_COEF,
    _dose,
    _hist,
    _ref,
)

GC = G_GRAD * C_BOX
FACE_MISS = 2.5379e-01  # item 727's control, `20261101T195309Z_MAT-21-step28-W1.log`
FACE_TOL = 1e-5  # the print's last digit, the brief's


def _say(msg):
    print(msg, flush=True)


def _dz_drift(p, t_w):
    return np.array([_drift(p, t + t_w, roots_z=DZ_ROOTS) for t in GRID])


def test_preconditions():
    assert ADJ_PT[2] / C_BOX == 0.9 and abs(P5[2] / C_BOX - 1.0) <= 1e-14
    assert abs(GRID[-1] / T_STAR - 1.0) <= 1e-14
    _say(f"MAT21-31 ADJ_PT z/C {ADJ_PT[2] / C_BOX:g}, P5 z/C {P5[2] / C_BOX:g}, {GRID.size} rows, "
         f"{DZ_ROOTS.size} z roots; bars DOSE {DOSE_BAR:g}, SYM0 {SYM0_BAR:g}")


def test_rep_face_control_and_drift_at_p5():
    ref = _ref(P5)
    miss = abs(_dose(ref + _dz_drift(P5, 0.0)) / _dose(_hist(P5, 0.0)) - 1.0)
    _say(f"MAT21-31 (Rep) P5, t_w = 0: miss {miss:.4e} vs item 727's {FACE_MISS:.4e}, |diff| "
         f"{abs(miss - FACE_MISS):.3e} (bound {FACE_TOL:g}, ASSERTED)")
    assert abs(miss - FACE_MISS) <= FACE_TOL, miss

    worst = 0.0
    for c, w in zip(WAIT_COEF, WAITS):
        row = np.max(np.abs(_dz_drift(P5, w))) / abs(GC)
        worst = max(worst, row)
        _say(f"MAT21-31 control P5 t_w lambda_1 = {c}: max |Dirichlet drift|/(g C) {row:.3e}")
    _say(f"MAT21-31 control P5 worst {worst:.3e} (bar {SYM0_BAR:g}, ASSERTED)")
    assert worst <= SYM0_BAR, worst


def test_dz_int_and_q_w_at_adj_pt():
    ref = _ref(ADJ_PT)
    d_ref = _dose(ref)
    worst_q = 0.0
    first_miss = None
    _say("MAT21-31 TABLE t_w lambda_1 | D Dirichlet s | D Robin s | miss | Robin D/D_ref-1 | "
         "Dirichlet D/D_ref-1 | Dirichlet-Robin drift at t* / (g C)")
    for c, w in zip(WAIT_COEF, WAITS):
        dz = _dz_drift(ADJ_PT, w)
        d_dz, _quad, rel, _swap, _rc = _q_check_g(T_INF + ref + dz, GRID)
        assert rel <= DOSE_BAR, (c, rel)
        worst_q = max(worst_q, rel)
        th = _hist(ADJ_PT, w)
        d_rb = _dose(th)
        miss = abs(d_dz / d_rb - 1.0)
        first_miss = miss if first_miss is None else first_miss
        diff = (dz[-1] - (th[-1] - ref[-1])) / GC
        _say(f"MAT21-31 TABLE {c} | {d_dz:.6e} | {d_rb:.6e} | {miss:.4e} | "
             f"{d_rb / d_ref - 1.0:+.6e} | {d_dz / d_ref - 1.0:+.6e} | {diff:+.6e}")
    _say(f"MAT21-31 (Q-w) Dirichlet histories worst rel {worst_q:.3e} (bar {DOSE_BAR:g}, ASSERTED)")
    _say(f"MAT21-31 (Dz-int) ADJ_PT, t_w = 0: miss {first_miss:.4e} (asserted > {DOSE_BAR:g})")
    assert first_miss > DOSE_BAR, first_miss
