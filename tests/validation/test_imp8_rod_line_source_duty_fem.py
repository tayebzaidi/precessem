"""`IMP-8` step 5 (§9 item 799, §10 I60): the duty-cycled line source in FEM, on
item 803's graded gel box (``_graded_box`` via its ``box`` fixture) and rod
source (``_source_and_s0``), both imported unchanged; real build, degree 2,
CN (θ = 0.5, ``startup_half_steps=2``) to t* = 900 s.

Brief: docs/planning/chunks/IMP-8.md, "§9 item 799" (with the 2026-11-03 10:30
ruling 6 and the 18:00 ruling 5: it depends on item 803 and runs at ``-n 8``).
The schedule goes through ``source_schedule`` of ``solve_transient_conduction``:
period 450 s, on for [jP, jP + D·P), off after.  Δt = 22.5 s, so a phase at
D = 0.5 is 10 steps; every switch time is asserted to lie on the step grid.

Windows are selected by ``FEM_EM_IMP8_I60_DUTY`` = ``0.5`` (W1, the gates) or
``1`` (W2, the control); unset, every solve test skips.

Asserted in both windows: (E) max balance residual ≤ ``E_BAND``; (S0) ∫q = P
within ``B_BAND`` with the selected-cell count > 0; (I) θ_FEM against item
787's Duhamel sum ``_duhamel(r, z, t, D)`` (per watt) at item 778's four
points within ``INTERIOR_BAND``, read at t* = 900 s (the end of an off phase at
D = 0.5) and t₂ = 675 s (the end of an on phase), both from one solve (t₂ is
the step-30 snapshot).  No (St).

Control (the weekly's, a reproduction): the D = 1 run reproduces item 803's θ
at (10, 0) mm to ``.6e`` — the string ``4.651377e+00`` printed at
`20261104T030134Z_IMP-8-I55-W2-item803.log:141`.  Both strings are printed.
"""

from __future__ import annotations

import os
import time

import numpy as np
import pytest

from fem_em_solver.thermal import solve_transient_conduction

from tests.validation.test_imp8_rod_line_source_closed_form import MM, T_PRINT
from tests.validation.test_imp8_rod_line_source_duty_closed_form import PERIOD, _duhamel
from tests.validation.test_imp8_rod_line_source_fem import (  # noqa: F401  (box: fixture)
    DEGREE,
    GATE_PTS,
    T_INF,
    W0C_CELLS,
    _source_and_s0,
    _theta_pts,
    box,
)
from tests.validation.test_mat10_sar_to_heat_sphere import B_BAND, ROBIN_TAG, _say
from tests.validation.test_mat12_transient_sar_heating_sphere import E_BAND, RHO_C
from tests.validation.test_mat13_transient_perfusion_sphere import INTERIOR_BAND
from tests.validation.test_mat15_coil_driven_heating_phantom import H_GEL, K_GEL

DUTY_ENV = "FEM_EM_IMP8_I60_DUTY"
DT = 22.5  # s
N_STEPS = 40  # to t* = 900 s
T2 = 675.0  # s, the end of an on phase at D = 0.5
STEP_T2 = 30
#: Item 803's θ_FEM at (10, 0) mm, t* = 900 s, D = 1, as printed
#: (`20261104T030134Z_IMP-8-I55-W2-item803.log:141`).
ITEM803_10_0 = "4.651377e+00"
CTRL_PT = (10.0, 0.0)


def _duty():
    raw = os.environ.get(DUTY_ENV, "").strip()
    if not raw:
        return None
    if raw not in ("0.5", "1"):
        raise ValueError(f"{DUTY_ENV} must be unset, '0.5' or '1', got {raw!r}")
    return float(raw)


def _schedule(duty):
    def chi(t):
        return 1.0 if (t % PERIOD) < duty * PERIOD else 0.0
    return chi


def _switches(duty, t_end=T_PRINT):
    out, j = [], 0
    while j * PERIOD <= t_end:
        out.append(j * PERIOD)
        if duty < 1.0 and j * PERIOD + duty * PERIOD <= t_end:
            out.append(j * PERIOD + duty * PERIOD)
        j += 1
    return out


def test_schedule_on_step_grid():
    """Every switch lies on the Δt grid; 10 steps per phase at D = 0.5."""
    assert DT * N_STEPS == T_PRINT and DT * STEP_T2 == T2
    assert 0.5 * PERIOD / DT == 10.0
    for d in (0.5, 1.0):
        for ts in _switches(d):
            k = ts / DT
            assert k == round(k), (d, ts)
    chi = _schedule(0.5)
    mids = [(i + 0.5) * DT for i in range(N_STEPS)]
    on = [chi(t) for t in mids]
    assert on == ([1.0] * 10 + [0.0] * 10) * 2
    assert all(_schedule(1.0)(t) == 1.0 for t in mids)


@pytest.fixture(scope="module")
def duty_run(box):
    duty = _duty()
    if duty is None:
        pytest.skip(f"{DUTY_ENV} unset: no solve in the coding window")
    msh, ct, ft, n = box
    assert n == W0C_CELLS, n  # item 803's grading, unchanged
    q, out = _source_and_s0(msh)
    out["duty"] = duty
    t1 = time.perf_counter()
    res = solve_transient_conduction(
        msh, ft, k=K_GEL, rho_c=RHO_C, q=q, h=H_GEL, t_inf=T_INF, t_init=T_INF,
        robin_tag=ROBIN_TAG, dt=DT, n_steps=N_STEPS, theta=0.5, startup_half_steps=2,
        degree=DEGREE, snapshot_steps=(STEP_T2,), source_schedule=_schedule(duty))
    out["t_call"] = time.perf_counter() - t1
    out["E"] = float(res.max_balance_residual)
    mult = list(res.source_multipliers)
    _say(f"\n[imp8-duty] D={duty:g} CN N={N_STEPS} dt={DT:g} s: t_call {out['t_call']:.1f} s; "
         f"n_factorizations {res.n_factorizations}; ranks {msh.comm.size}; multipliers "
         f"{len(mult)} (sum {sum(mult):g}); (S0) {out['S0']:.3e}; selected cells {out['n_sel']}")
    out["mult"] = mult
    out["I"], out["th"] = [], {}
    for label, t, T in (("t*", T_PRINT, res.temperature), ("t2", T2, res.snapshots[STEP_T2])):
        th = _theta_pts(T, GATE_PTS)
        for (r, z), v in zip(GATE_PTS, th):
            ex = _duhamel(r * MM, z * MM, t, duty)
            rel = abs(v / ex - 1.0)
            out["I"].append(rel)
            out["th"][(label, r, z)] = v
            _say(f"[imp8-duty] (I) D={duty:g} ({r:g}, {z:g}) mm {label}={t:g} s: theta_FEM "
                 f"{v:.6e} K/W, Duhamel {ex:.6e}, rel {rel:.3e} (gated, band "
                 f"{INTERIOR_BAND:.0e})")
    _say(f"[imp8-duty] (E) {out['E']:.3e} (band {E_BAND:.0e}); (I) max {max(out['I']):.3e}")
    if duty == 1.0:
        s = f"{out['th'][('t*',) + CTRL_PT]:.6e}"
        out["ctrl"] = s
        _say(f"[imp8-duty] control D=1 theta(10, 0) mm t*: this run {s}; item 803 "
             f"{ITEM803_10_0}; reproduced to .6e: {s == ITEM803_10_0}")
    return out


def test_e_balance(duty_run):
    assert duty_run["E"] <= E_BAND, duty_run["E"]


def test_s0_source_power(duty_run):
    assert duty_run["n_sel"] > 0
    assert duty_run["S0"] <= B_BAND, duty_run["S0"]


def test_schedule_applied(duty_run):
    """The (half-)step multipliers the solver used: 2 + 39 of them."""
    mult = duty_run["mult"]
    assert len(mult) == N_STEPS + 1
    if duty_run["duty"] == 1.0:
        assert all(m == 1.0 for m in mult)
    else:
        assert mult == [1.0, 1.0] + [1.0] * 9 + ([0.0] * 10 + [1.0] * 10) + [0.0] * 10


def test_i_interior_against_duhamel(duty_run):
    assert max(duty_run["I"]) <= INTERIOR_BAND, duty_run["I"]


def test_control_d1_reproduces_item803(duty_run):
    if duty_run["duty"] != 1.0:
        pytest.skip("the control is the D = 1 window")
    assert duty_run["ctrl"] == ITEM803_10_0, duty_run["ctrl"]
