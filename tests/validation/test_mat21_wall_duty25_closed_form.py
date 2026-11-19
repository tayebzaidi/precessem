"""`MAT-21` step 25 — the duty-cycled walled slab at the probe heights, D = 0.5 and 0.25
(§9 item 713, §10 lane 3 → chain P, P53; brief `docs/planning/chunks/MAT-21.md` "§9 item
713"; the comparand of item 721).

Numpy.  Sets no band and moves none; every bar is `MAT-17`'s, imported.  The walled
comparand is one-dimensional in z, so P37's six points take two heights, z = 0 (P0–P4) and
z = C (P5), and P50a's near-wall height 0.9 C is added so (D1) / (Duh) are not read at two
heights only.

Asserted, at every height: (D1) ``duty = 1.0`` equals the step response; at both duties,
(Duh) the duty form from rest equals the Duhamel sum of ``N_CYC`` on / off steps; (Edge′a)
modal amplitudes continuous at every switch; (Edge′b) the jump minus Δs·r_N; (Col) the five
z = 0 series values bitwise equal.  Negative controls: the continuous series misses at t\\*
> ``NEG_FLOOR`` at both duties; the massless-wall duty series is asserted at D = 0.5 and at
D = 0.25 only if it clears the floor (a record otherwise, per the brief).

Not claimed: a FEM run (P54), a dose, lateral structure.
"""

from __future__ import annotations

import numpy as np

from fem_em_solver.utils.analytical import (
    _two_layer_slab_duty_state,
    two_layer_slab_mode,
    two_layer_slab_steady_temperature,
    two_layer_slab_transient_coefficients,
)
from tests.validation.test_mat12_transient_sar_heating_sphere import N_MAIN
from tests.validation.test_mat15_finite_cylinder_closed_form import H, Q
from tests.validation.test_mat21_box_closed_form import C_BOX
from tests.validation.test_mat21_box_duty25_dose_closed_form import DUTY25
from tests.validation.test_mat21_box_probe_points_fem import PROBE_LABELS, PROBE_PTS
from tests.validation.test_mat21_wall_duty_closed_form import (
    D1_BAND,
    DUH_BAND,
    EDGE_BAND,
    N_CYC,
    PERIOD,
    T2,
    _duty,
    _lam,
    _ss0,
    _step,
)
from tests.validation.test_mat21_wall_transient_closed_form import (
    NEG_FLOOR,
    T_STAR,
    WALL,
    Z_GEL,
)
from tests.validation.test_mat21_wall_transient_fem import MASSLESS

DUTIES = (0.5, DUTY25)
HEIGHTS = np.array([0.0, 0.9 * C_BOX, C_BOX])
DT = T_STAR / N_MAIN


def _say(msg):
    print(msg, flush=True)


def _edges(duty):
    """Every switch instant in (0, t₂]."""
    return sorted(
        [j * PERIOD for j in range(1, N_CYC + 1)]
        + [j * PERIOD + duty * PERIOD for j in range(N_CYC)]
    )


def test_preconditions_switches_on_the_grid():
    assert WALL["half_thickness"] == C_BOX, (WALL["half_thickness"], C_BOX)
    assert abs(N_CYC * PERIOD / T2 - 1.0) <= 1e-14
    for d in DUTIES:
        off = max(abs(te / DT - round(te / DT)) for te in _edges(d))
        _say(
            f"MAT21-25 D={d}: on-phase {d * PERIOD / DT:.9f} DT; {len(_edges(d))} switches, max "
            f"off-grid {off:.3e} DT"
        )
        assert off <= 1e-9, (d, off)


def test_d1_duty_one_equals_step_response():
    ss0 = _ss0()
    worst = 0.0
    for f in (0.3, 1.0, 2.0, 8.0):
        t = f * T_STAR
        worst = max(
            worst, float(np.abs(_duty(HEIGHTS, t, duty=1.0) - _step(HEIGHTS, t)).max() / ss0)
        )
    _say(f"MAT21-25 (D1) max {worst:.3e} of theta_ss(0) (band {D1_BAND:g}, ASSERTED)")
    assert worst <= D1_BAND, worst


def test_duh_from_rest_equals_duhamel_superposition():
    ss0 = _ss0()
    for d in DUTIES:
        on_edges = [j * PERIOD for j in range(N_CYC)]
        off_edges = [j * PERIOD + d * PERIOD for j in range(N_CYC)]
        worst = 0.0
        for j in range(N_CYC):
            for i in range(8):
                t = (j + i / 8.0) * PERIOD
                duh = np.zeros(HEIGHTS.shape)
                for te in on_edges:
                    if te <= t:
                        duh += _step(HEIGHTS, t - te)
                for te in off_edges:
                    if te <= t:
                        duh -= _step(HEIGHTS, t - te)
                worst = max(worst, float(np.abs(_duty(HEIGHTS, t, duty=d) - duh).max() / ss0))
        _say(
            f"MAT21-25 (Duh) D={d}: cycles 0-{N_CYC - 1}, 8 fractions/cycle, {HEIGHTS.size} "
            f"heights: max {worst:.3e} of theta_ss(0) (band {DUH_BAND:g}, ASSERTED)"
        )
        assert worst <= DUH_BAND, (d, worst)


def test_edge_prime_at_every_switch():
    ss0 = _ss0()
    lam = _lam()
    a = two_layer_slab_transient_coefficients(lam, q=Q, **WALL)
    resid = (
        two_layer_slab_steady_temperature(
            Z_GEL,
            q=Q,
            h=H,
            half_thickness=WALL["half_thickness"],
            wall=WALL["wall"],
            k=WALL["k"],
            k_w=WALL["k_w"],
        )
        - two_layer_slab_mode(Z_GEL, lam, **WALL) @ a
    )
    for d in DUTIES:
        jumps, a_jumps, id_err = [], [], []
        for te in _edges(d):
            sL, AL = _two_layer_slab_duty_state(te, lam, a, period=PERIOD, duty=d, side="left")
            sR, AR = _two_layer_slab_duty_state(te, lam, a, period=PERIOD, duty=d, side="right")
            thL = _duty(Z_GEL, te, duty=d, side="left")
            thR = _duty(Z_GEL, te, duty=d, side="right")
            jumps.append(np.abs(thR - thL).max() / ss0)
            a_jumps.append(np.abs(AR - AL).max() / np.abs(a).max())
            id_err.append(np.abs((thR - thL) - (sR - sL) * resid).max() / ss0)
        _say(
            f"MAT21-25 (Edge') D={d}: {len(jumps)} switches; (a) max |A_R - A_L|/max a "
            f"{max(a_jumps):.3e}; (b) jump - ds*r_N {max(id_err):.3e} (band {EDGE_BAND:g}, "
            f"ASSERTED); raw jump {max(jumps):.3e} (record)"
        )
        assert max(a_jumps) <= EDGE_BAND, (d, max(a_jumps))
        assert max(id_err) <= EDGE_BAND, (d, max(id_err))


def test_col_and_printed_six_points():
    """(Col) is the 1D reduction read at P0–P4; the rows are records."""
    for d in DUTIES:
        for name, t in (("t*", T_STAR), ("t2", T2)):
            vals = [float(_duty(p[2], t, duty=d, side="left")[0]) for p in PROBE_PTS]
            _say(
                f"MAT21-25 (Col) D={d} {name}: "
                + "; ".join(f"{lab} {v:.6f} K" for lab, v in zip(PROBE_LABELS, vals))
            )
            assert len(set(vals[:5])) == 1, (d, name, vals)
            assert vals[5] != vals[0], (d, name, vals)


def test_negative_controls():
    ss0 = _ss0()
    for d in DUTIES:
        nom = _duty(Z_GEL, T_STAR, duty=d, side="left")
        cont = _step(Z_GEL, T_STAR)
        miss = float(np.abs(cont - nom).max() / ss0)
        ml = _duty(Z_GEL, T_STAR, duty=d, side="left", rho_c_w=MASSLESS * WALL["rho_c_w"])
        miss_ml = float(np.abs(ml - nom).max() / ss0)
        above = miss_ml > NEG_FLOOR
        _say(
            f"MAT21-25 (neg) D={d}: continuous miss {miss:.4e} (floor {NEG_FLOOR}, ASSERTED); "
            f"massless-wall miss {miss_ml:.4e} (floor {NEG_FLOOR}; item 689 read 1.9941e-3 at "
            f"D = 0.5): {'ABOVE floor, asserted' if above else 'BELOW floor, a record'}"
        )
        assert miss > NEG_FLOOR, (d, miss)
        if d == 0.5:
            assert above, (d, miss_ml)
