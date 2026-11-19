"""`IMP-8` step 8 (§9 item 814, §10 I63): the two-tip source under a duty cycle in
FEM, on item 803's graded gel box (``_graded_box`` via its ``box`` fixture) with
item 796's two-ball source (``_tips_source``, P/2 on each ball), both imported
unchanged; real build, degree 2, CN (θ = 0.5, ``startup_half_steps=2``) to
t* = 900 s at D = 0.5, period 450 s, Δt = 22.5 s (item 799's ``_schedule``,
``_switches``, ``DT``, ``N_STEPS``, ``STEP_T2``), one solve, t₂ = 675 s the
step-30 snapshot.

Brief: docs/planning/chunks/IMP-8.md, "§9 item 814".

Solve tests run only when ``FEM_EM_IMP8_I63_DUTY`` = ``0.5``; unset, they skip
(the W0 coding window).

Asserted: (E) max balance residual ≤ ``E_BAND``; (S0) ∫q = P and each ball's
∫q = P/2 within ``B_BAND`` with selected-cell counts > 0; the multipliers as item
799 asserted them; (I) θ_FEM against item 811's Φ₂ Duhamel sum
``_duhamel(_phi_tips, r, z, t, 0.5)`` (per watt) at item 803's four gate points,
at t* and t₂, within ``INTERIOR_BAND``.

Control (asserted, factor 1): the continuous Φ₂ comparand (D = 1 step response
``SCALE·_phi_tips``) misses θ_FEM at (10, 0) mm at t* by > ``INTERIOR_BAND``;
item 811's closed-form control measured 6.9531e-01 there.  The same miss at the
other three points is printed only.

Printed (records, no band): θ at (2, 15) and (2, 0) mm at t* and t₂ beside item
811's closed form, with (a/R₊)²; the tips / line ratio at (10, 0) mm, t*, over
item 799's θ_FEM 2.240893 K/W (`20261104T044429Z_IMP-8-I60-W1-item799.log`);
t_build, t_call and ``n_factorizations``.
"""

from __future__ import annotations

import os
import time

import numpy as np
import pytest
import ufl

from fem_em_solver.thermal import solve_transient_conduction
from tests.mesh.test_birdcage_implant_rod import ROD_RADIUS
from tests.validation.test_imp8_rod_line_source_closed_form import MM, T_PRINT, L, _phi_tips
from tests.validation.test_imp8_rod_line_source_duty_fem import (
    DT,
    N_STEPS,
    STEP_T2,
    T2,
    _schedule,
    _switches,
)
from tests.validation.test_imp8_rod_line_source_fem import (  # noqa: F401  (box: fixture)
    DEGREE,
    GATE_PTS,
    T_INF,
    W0C_CELLS,
    P,
    _theta_pts,
    box,
)
from tests.validation.test_imp8_rod_tips_duty_closed_form import SCALE, _duhamel
from tests.validation.test_imp8_rod_tips_fem import SIGNS, _tips_source
from tests.validation.test_mat10_sar_to_heat_sphere import B_BAND, ROBIN_TAG, _allsum, _say
from tests.validation.test_mat12_transient_sar_heating_sphere import E_BAND, RHO_C
from tests.validation.test_mat13_transient_perfusion_sphere import INTERIOR_BAND
from tests.validation.test_mat15_coil_driven_heating_phantom import H_GEL, K_GEL

DUTY_ENV = "FEM_EM_IMP8_I63_DUTY"
DUTY = 0.5
CTRL_PT = (10.0, 0.0)
RECORD_PTS = ((2.0, 15.0), (2.0, 0.0))  # (r, z) mm, printed only
#: Item 799's line-source θ_FEM at (10, 0) mm, t*, D = 0.5, K/W
#: (`20261104T044429Z_IMP-8-I60-W1-item799.log`).
LINE_FEM_10_0 = 2.240893
#: Item 811's closed-form control at (10, 0) mm (`20261104T103134Z_IMP-8-item811-W1.log`),
#: printed beside the measured miss.
CTRL_ITEM811 = 6.9531e-01


def _enabled():
    raw = os.environ.get(DUTY_ENV, "").strip()
    if not raw:
        return False
    if raw != "0.5":
        raise ValueError(f"{DUTY_ENV} must be unset or '0.5', got {raw!r}")
    return True


def test_schedule_on_step_grid_and_comparand():
    """Item 799's grid at D = 0.5; the comparand is per watt and zero at rest."""
    assert DT * N_STEPS == T_PRINT and DT * STEP_T2 == T2
    for ts in _switches(DUTY):
        k = ts / DT
        assert k == round(k), ts
    mids = [(i + 0.5) * DT for i in range(N_STEPS)]
    assert [_schedule(DUTY)(t) for t in mids] == ([1.0] * 10 + [0.0] * 10) * 2
    assert _duhamel(_phi_tips, 10 * MM, 0.0, 0.0, DUTY) == 0.0
    # at D = 1 the sum is the step response (item 811's (D1), restated as a precondition)
    v1 = _duhamel(_phi_tips, 10 * MM, 0.0, T_PRINT, 1.0)
    assert abs(v1 / (SCALE * _phi_tips(10 * MM, 0.0, T_PRINT)) - 1.0) < 1e-13


@pytest.fixture(scope="module")
def run(box):
    if not _enabled():
        pytest.skip(f"{DUTY_ENV} unset: no solve in the coding window")
    msh, ct, ft, n = box
    assert n == W0C_CELLS, n  # item 803's grading, unchanged
    t0 = time.perf_counter()
    q, balls = _tips_source(msh)
    out = {"S0": abs(float(np.real(_allsum(msh.comm, q * ufl.dx))) / P - 1.0), "balls": []}
    for s, (n_sel, v_sel, qb) in zip(SIGNS, balls):
        pb = float(np.real(_allsum(msh.comm, qb * ufl.dx)))
        out["balls"].append((n_sel, abs(pb / (0.5 * P) - 1.0)))
        _say(
            f"[imp8-tdf] ball z={s * 0.5 * L / MM:+g} mm: cells {n_sel}; "
            f"int q / (P/2) - 1 = {pb / (0.5 * P) - 1.0:.3e}"
        )
    out["t_build"] = time.perf_counter() - t0
    t1 = time.perf_counter()
    res = solve_transient_conduction(
        msh,
        ft,
        k=K_GEL,
        rho_c=RHO_C,
        q=q,
        h=H_GEL,
        t_inf=T_INF,
        t_init=T_INF,
        robin_tag=ROBIN_TAG,
        dt=DT,
        n_steps=N_STEPS,
        theta=0.5,
        startup_half_steps=2,
        degree=DEGREE,
        snapshot_steps=(STEP_T2,),
        source_schedule=_schedule(DUTY),
    )
    out["t_call"] = time.perf_counter() - t1
    out["E"] = float(res.max_balance_residual)
    out["mult"] = list(res.source_multipliers)
    _say(
        f"\n[imp8-tdf] D={DUTY:g} CN N={N_STEPS} dt={DT:g} s: t_build {out['t_build']:.1f} s; "
        f"t_call {out['t_call']:.1f} s; n_factorizations {res.n_factorizations}; ranks "
        f"{msh.comm.size}; multipliers {len(out['mult'])}; (S0) {out['S0']:.3e}"
    )
    out["I"], out["ctrl"] = [], None
    pts = GATE_PTS + RECORD_PTS
    for label, t, T in (("t*", T_PRINT, res.temperature), ("t2", T2, res.snapshots[STEP_T2])):
        th = _theta_pts(T, pts)
        for (r, z), v in zip(pts, th):
            ex = _duhamel(_phi_tips, r * MM, z * MM, t, DUTY)
            rel = abs(v / ex - 1.0)
            gated = (r, z) in GATE_PTS
            a2 = (ROD_RADIUS / np.hypot(r * MM, z * MM - 0.5 * L)) ** 2
            if gated:
                out["I"].append(rel)
            _say(
                f"[imp8-tdf] (I) ({r:g}, {z:g}) mm {label}={t:g} s: theta_FEM {v:.6e} K/W, "
                f"Phi_2 Duhamel {ex:.6e}, rel {rel:.3e} "
                f"({'gated' if gated else 'printed only'}, band {INTERIOR_BAND:.0e}); "
                f"(a/R+)^2 {a2:.3e}"
            )
            if label == "t*" and gated:
                cont = SCALE * _phi_tips(r * MM, z * MM, T_PRINT)
                miss = abs(v / cont - 1.0)
                is_ctrl = (r, z) == CTRL_PT
                if is_ctrl:
                    out["ctrl"] = miss
                _say(
                    f"[imp8-tdf] (c-cont) ({r:g}, {z:g}) mm t*: continuous Phi_2 {cont:.6e}; "
                    f"miss |theta_FEM/cmp - 1| {miss:.4e} "
                    f"({'asserted > INTERIOR_BAND' if is_ctrl else 'printed only'}; item 811 "
                    f"closed-form control at (10, 0) mm {CTRL_ITEM811:.4e})"
                )
            if label == "t*" and (r, z) == CTRL_PT:
                _say(
                    f"[imp8-tdf] record (no band): tips / line at (10, 0) mm t*, FEM/FEM = "
                    f"{v / LINE_FEM_10_0:.6f} (over item 799's {LINE_FEM_10_0:.6f} K/W)"
                )
    _say(
        f"[imp8-tdf] (E) {out['E']:.3e} (band {E_BAND:.0e}); (I) max {max(out['I']):.3e}; "
        f"(c-cont) {out['ctrl']:.4e}"
    )
    return out


def test_e_balance(run):
    assert run["E"] <= E_BAND, run["E"]


def test_s0_source_power(run):
    assert run["S0"] <= B_BAND, run["S0"]
    for n_sel, rel in run["balls"]:
        assert n_sel > 0
        assert rel <= B_BAND, rel


def test_schedule_applied(run):
    mult = run["mult"]
    assert len(mult) == N_STEPS + 1
    assert mult == [1.0, 1.0] + [1.0] * 9 + ([0.0] * 10 + [1.0] * 10) + [0.0] * 10


def test_i_interior_against_duhamel(run):
    assert len(run["I"]) == 2 * len(GATE_PTS)
    assert max(run["I"]) <= INTERIOR_BAND, run["I"]


def test_control_continuous_misses(run):
    assert run["ctrl"] is not None
    assert run["ctrl"] > INTERIOR_BAND, run["ctrl"]
