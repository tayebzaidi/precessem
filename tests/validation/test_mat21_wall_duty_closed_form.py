"""`MAT-21` step 21a — the duty-cycled two-layer (gel + wall) slab in closed
form (§9 item 687, §10 P49a; the FEM comparand of item 688).

`utils.analytical.two_layer_slab_duty_cycle_temperature`: source on while
0 ≤ (t mod P) < D·P, from rest; θ = s(t) θ_ss − Σ (aₙ s(t) − Aₙ) Xₙ with
Ȧₙ = λₙ(aₙ s − Aₙ) advanced exactly phase by phase (derivation in the
docstring).  Errors are on θ = T − T∞, relative to θ_ss(0).

Fixture constants imported by name from item 677's module (`WALL`, `T_STAR`,
`N_SERIES`, `Z_GEL`, `NEG_FLOOR`) and `MAT-15` (`H`, `Q`); the massless-wall
scale `MASSLESS` from item 682's module.  Schedule: item 675's **P = 0.5 t\\*,
D = 0.5** (`test_imp7_duty_cycle_fem.py:63–64`).  `DUTY` is imported by name;
item 675's `PERIOD` is 0.5 × *IMP-7's* t\\*, so the ratio `PERIOD / T_STAR`
is imported (asserted == 0.5) and applied to this module's t\\* — item 688's
"a period is 20 steps" at Δt = t\\*/40 requires P = 0.5 t\\* of `MAT-21`.

Bars copied from `MAT-17` step 1 (`test_mat17_duty_cycle_closed_form.py:206`,
`:219`): (D1) 1e-13, (Duh) 1e-10; (Edge′) 1e-13 on the jump identity and modal continuity (Review ruling 2026-10-28 10:30).
"""

from __future__ import annotations

import time

import numpy as np
import pytest

from fem_em_solver.utils.analytical import (
    _two_layer_slab_duty_state,
    two_layer_slab_duty_cycle_temperature,
    two_layer_slab_eigenvalues,
    two_layer_slab_mode,
    two_layer_slab_steady_temperature,
    two_layer_slab_transient_coefficients,
    two_layer_slab_transient_temperature,
)
from tests.validation.test_imp7_duty_cycle_fem import DUTY
from tests.validation.test_imp7_duty_cycle_fem import PERIOD as PERIOD_IMP7
from tests.validation.test_imp7_transient_closed_form import T_STAR as T_STAR_IMP7
from tests.validation.test_mat15_finite_cylinder_closed_form import H, Q
from tests.validation.test_mat21_wall_transient_closed_form import (
    N_SERIES,
    NEG_FLOOR,
    T_STAR,
    WALL,
    Z_GEL,
)
from tests.validation.test_mat21_wall_transient_fem import MASSLESS

PERIOD_RATIO = PERIOD_IMP7 / T_STAR_IMP7
PERIOD = PERIOD_RATIO * T_STAR
T2 = 8 * T_STAR
N_CYC = 16  # cycles 0–15, to t₂ = 8 t* = 16 P
D1_BAND = 1.0e-13  # MAT-17 P4 (test_mat17_duty_cycle_closed_form.py:219)
DUH_BAND = 1.0e-10  # MAT-17 P3 (test_mat17_duty_cycle_closed_form.py:206)
EDGE_BAND = 1.0e-13  # the brief's (Edge)


def _say(msg):
    print(msg, flush=True)


_CACHE = {}


def _lam(**over):
    key = tuple(sorted(over.items()))
    if key not in _CACHE:
        _CACHE[key] = two_layer_slab_eigenvalues(h=H, n=N_SERIES, **dict(WALL, **over))
    return _CACHE[key]


def _duty(z, t, duty=DUTY, side="right", **over):
    return two_layer_slab_duty_cycle_temperature(
        z,
        t,
        q=Q,
        h=H,
        period=PERIOD,
        duty=duty,
        t_inf=0.0,
        n_terms=N_SERIES,
        eigenvalues=_lam(**over),
        side=side,
        **dict(WALL, **over),
    )


def _step(z, tau):
    if tau < 0.0:
        return np.zeros(np.atleast_1d(z).shape)
    return two_layer_slab_transient_temperature(
        z, tau, q=Q, h=H, t_inf=0.0, n_terms=N_SERIES, eigenvalues=_lam(), **WALL
    )


def _ss0():
    return float(
        two_layer_slab_steady_temperature(
            0.0,
            q=Q,
            h=H,
            half_thickness=WALL["half_thickness"],
            wall=WALL["wall"],
            k=WALL["k"],
            k_w=WALL["k_w"],
        )
    )


def test_schedule_is_item_675s():
    _say(
        f"MAT21-21a schedule: imported PERIOD/T_STAR (IMP-7) = {PERIOD_RATIO!r}, DUTY = {DUTY!r}; "
        f"MAT-21 t* = {T_STAR:.4f} s, P = {PERIOD:.4f} s, t2 = 8 t* = {T2:.4f} s"
    )
    assert PERIOD_RATIO == 0.5
    assert DUTY == 0.5


# ---------------------------------------------------------------- (D1)
def test_d1_duty_one_equals_step_response_and_validation():
    ss0 = _ss0()
    worst = 0.0
    for f in (0.3, 1.0, 2.0, 8.0):
        t = f * T_STAR
        e = np.abs(_duty(Z_GEL, t, duty=1.0) - _step(Z_GEL, t)).max() / ss0
        _say(
            f"MAT21-21a (D1) duty 1 vs two_layer_slab_transient_temperature at {f} t*: "
            f"{e:.3e} of theta_ss(0) (band {D1_BAND:g})"
        )
        worst = max(worst, e)
    for bad in (0.0, -0.1, 1.5):
        with pytest.raises(ValueError):
            _duty(0.0, T_STAR, duty=bad)
    with pytest.raises(ValueError):
        two_layer_slab_duty_cycle_temperature(
            0.0, 1.0, q=Q, h=H, period=0.0, duty=0.5, eigenvalues=_lam(), **WALL
        )
    _say(f"MAT21-21a (D1) max {worst:.3e}; ValueError raised for duty 0, -0.1, 1.5 and P = 0")
    assert worst <= D1_BAND


# ---------------------------------------------------------------- (Duh)
def test_duh_from_rest_equals_duhamel_superposition():
    ss0 = _ss0()
    t0 = time.perf_counter()
    on_edges = [j * PERIOD for j in range(N_CYC)]
    off_edges = [j * PERIOD + DUTY * PERIOD for j in range(N_CYC)]
    worst, first_fail = 0.0, None
    per_cycle = []
    for j in range(N_CYC):
        cyc = 0.0
        for i in range(8):
            t = (j + i / 8.0) * PERIOD
            modal = _duty(Z_GEL, t)
            duh = np.zeros_like(modal)
            for te in on_edges:
                if te <= t:
                    duh += _step(Z_GEL, t - te)
            for te in off_edges:
                if te <= t:
                    duh -= _step(Z_GEL, t - te)
            e = np.abs(modal - duh).max() / ss0
            cyc = max(cyc, e)
            if e > DUH_BAND and first_fail is None:
                first_fail = (j, i, e, float(modal[20]), float(duh[20]))
        per_cycle.append(cyc)
        worst = max(worst, cyc)
    _say(
        f"MAT21-21a (Duh) from rest vs Duhamel step superposition, cycles 0-{N_CYC - 1}, "
        f"8 fractions/cycle, {Z_GEL.size} points: max {worst:.3e} of theta_ss(0) "
        f"(band {DUH_BAND:g}); per cycle "
        + " ".join(f"{c:.1e}" for c in per_cycle)
        + f"; {time.perf_counter() - t0:.1f} s"
    )
    if first_fail is not None:
        _say(
            f"MAT21-21a (Duh) first failing (cycle, fraction/8, err, modal centre, duh centre) "
            f"= {first_fail}"
        )
    assert worst <= DUH_BAND


# ---------------------------------------------------------------- (Edge)
def test_edge_continuity_and_off_phase_decay():
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
    jumps, a_jumps, id_err = [], [], []
    edges = sorted(
        [j * PERIOD for j in range(1, N_CYC + 1)]
        + [j * PERIOD + DUTY * PERIOD for j in range(N_CYC)]
    )
    for te in edges:
        sL, AL = _two_layer_slab_duty_state(te, lam, a, period=PERIOD, duty=DUTY, side="left")
        sR, AR = _two_layer_slab_duty_state(te, lam, a, period=PERIOD, duty=DUTY, side="right")
        thL = _duty(Z_GEL, te, side="left")
        thR = _duty(Z_GEL, te, side="right")
        jumps.append(np.abs(thR - thL).max() / ss0)
        a_jumps.append(np.abs(AR - AL).max() / np.abs(a).max())
        id_err.append(np.abs((thR - thL) - (sR - sL) * resid).max() / ss0)
    decay_ok, n_off = True, 0
    for j in range(N_CYC):
        ts = j * PERIOD + DUTY * PERIOD + np.linspace(0.0, (1 - DUTY) * PERIOD, 9)[1:-1]
        c = np.array([_duty(0.0, t)[0] for t in ts])
        n_off += 1
        decay_ok &= bool(np.all(np.diff(c) < 0.0))
    _say(
        f"MAT21-21a (Edge) {len(edges)} switch instants in (0, t2]: max |theta_R − theta_L| "
        f"(the N = 400 truncation residual, record) {max(jumps):.3e} of theta_ss(0); modal amplitudes max |A_R - A_L|/max a "
        f"{max(a_jumps):.3e}; jump - ds*(theta_ss - sum a X) {max(id_err):.3e} (band {EDGE_BAND:g}); "
        f"truncation residual max|theta_ss - sum a X| over Z_GEL {np.abs(resid).max() / ss0:.3e}"
    )
    _say(
        f"MAT21-21a (Edge) off phases with strictly decaying centre (7 interior samples each): "
        f"{'all' if decay_ok else 'NOT all'} of {n_off}"
    )
    assert max(a_jumps) == 0.0
    assert decay_ok
    assert max(id_err) <= EDGE_BAND


# ---------------------------------------------------------------- printed comparands
def test_printed_comparands_and_period_mean():
    ss0 = _ss0()
    vals = {}
    for name, t in (("t*", T_STAR), ("t2", T2)):
        vals[name] = (float(_duty(0.0, t, side="left")[0]), float(_duty(0.0, t)[0]))
    ts = np.linspace(T2 - PERIOD, T2, 801)
    c = np.array([_duty(0.0, t)[0] for t in ts])
    mean = float(np.sum((c[1:] + c[:-1]) * 0.5 * np.diff(ts)) / PERIOD)
    _say(
        f"MAT21-21a (print) theta(0, t*) left (off-phase end) {vals['t*'][0]:.9f} K, right "
        f"{vals['t*'][1]:.9f} K; theta(0, t2) left {vals['t2'][0]:.9f} K, right "
        f"{vals['t2'][1]:.9f} K; theta_ss(0) {ss0:.9f} K"
    )
    _say(
        f"MAT21-21a (print) period-mean of theta(0, t) over [t2 - P, t2] (trapezoid, 801 pts) "
        f"{mean:.6f} K vs D*theta_ss(0) {DUTY * ss0:.6f} K, ratio {mean / (DUTY * ss0):.6f} "
        f"(record, not a gate)"
    )
    assert np.isfinite(mean)  # print-only


# ---------------------------------------------------------------- negative controls
def test_negative_control_continuous_source_misses_duty_at_t_star():
    ss0 = _ss0()
    duty = _duty(Z_GEL, T_STAR, side="left")
    cont = _step(Z_GEL, T_STAR)
    miss = np.abs(cont - duty).max() / ss0
    _say(
        f"MAT21-21a (neg, asserted) continuous-source series vs duty series at t* "
        f"(off-phase end): max miss {miss:.4e} of theta_ss(0) (floor {NEG_FLOOR}); centre "
        f"{float(cont[20]):.6f} K vs {float(duty[20]):.6f} K"
    )
    assert miss > NEG_FLOOR


def test_second_control_massless_wall_duty_series():
    """Conditional record per the brief: asserted only if the miss clears the floor."""
    ss0 = _ss0()
    nom = _duty(Z_GEL, T_STAR, side="left")
    ml = _duty(Z_GEL, T_STAR, side="left", rho_c_w=MASSLESS * WALL["rho_c_w"])
    miss = np.abs(ml - nom).max() / ss0
    above = miss > NEG_FLOOR
    _say(
        f"MAT21-21a (ctl2) massless-wall (rho_w c_w x {MASSLESS:g}) duty series vs walled at t*: "
        f"max miss {miss:.4e} of theta_ss(0) (floor {NEG_FLOOR}); centre {float(ml[20]):.6f} K "
        f"vs {float(nom[20]):.6f} K; "
        + ("ABOVE floor -> asserted" if above else "BELOW floor -> recorded, not asserted")
    )
    if above:
        assert miss > NEG_FLOOR
