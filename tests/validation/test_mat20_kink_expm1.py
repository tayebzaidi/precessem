"""`MAT-20` R3 — `MAT-16`'s kink and ``expm1`` hazards read against ``cem43`` (§9 item 638).

`MAT-16`'s draft named two numerical hazards of a CEM43 integrator: the 43 °C
kink (R switches 0.25 → 0.5) and cancellation in ``(exp(x) − 1)/x`` at small
slope.  This module tests the public ``thermal.dose.cem43`` there, every
assertion at `MAT-20` step 1's ``DOSE_BAR`` (1e-10 relative, imported):

(K) the kink at ε = 1e-3: a single segment 43 − ε → 43 + ε (and its descent)
    over τ = 600 s against step 1's split two-piece closed form ``ramp_dose``;
    the mean dose/τ, its miss from 1 − (ln 2/4)ε and the predicted remainder
    5(ln 2)²ε²/12 are printed, not asserted;
(E) the small-slope limit: one non-crossing segment T₀ → T₀ + δ over 600 s,
    δ ∈ {1e-6, 1e-9, 6e-10}, T₀ ∈ {41, 45}, against the second-order Taylor
    form Δt·R^(43−T₀)·(1 + x/2 + x²/6), x = −δ ln R (truncation ≈ x³/24);
(Z) the constant segment (the ``x == 0`` branch): Δt·R^(43−T) at T = 41, 45.

Controls under rule (e), *predicted*, printed with the measured value and never
asserted: (Kc) R swapped on the ε-ramp moves the mean by ≈ +(ln 2/2)ε; (Ec) the
naive difference quotient ``ramp_dose`` misses the Taylor form at δ = 6e-10 by
≈ 2e-7 relative.  The asserted control is step 1's ``test_rc_r_swap_misses``,
re-run unedited in the same window.

Scope: a numerical-hazard read of one post-processor; no damage model,
threshold, biological or regulatory statement.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from fem_em_solver.thermal.dose import cem43
from tests.validation.test_mat20_cem43_closed_form import BAR as DOSE_BAR
from tests.validation.test_mat20_cem43_closed_form import ramp_dose

P = "[MAT-20 R3]"
EPS = 1e-3
TAU = 600.0
DT = 600.0
DELTAS = (1e-6, 1e-9, 6e-10)
T0S = {41.0: 0.25, 45.0: 0.5}  # T0 -> R on its side
LN2 = math.log(2.0)


def _rel(x, ref):
    return abs(x - ref) / abs(ref)


def _taylor(T0, T1, dt, r):
    """Δt·R^(43−T₀)·(1 + x/2 + x²/6), x = −(T₁ − T₀) ln R (float-exact ΔT)."""
    x = -(T1 - T0) * math.log(r)
    return dt * r ** (43.0 - T0) * (1.0 + x / 2.0 + x * x / 6.0), x


@pytest.mark.parametrize("direction", ["rise", "fall"])
def test_k_kink(direction):
    a, b = (43.0 - EPS, 43.0 + EPS) if direction == "rise" else (43.0 + EPS, 43.0 - EPS)
    got = cem43(np.array([0.0, TAU]), np.array([a, b]))
    ref = ramp_dose(a, b, TAU)
    err = _rel(got, ref)
    mean = got / TAU
    lin = 1.0 - (LN2 / 4.0) * EPS
    pred_rem = 5.0 * LN2**2 * EPS**2 / 12.0
    print(
        f"{P} K {direction} eps={EPS:g} op={got:.15e} ref={ref:.15e} rel={err:.3e} "
        f"(bar {DOSE_BAR:g})"
    )
    print(
        f"{P} K {direction} mean={mean:.15f} miss_from_1-(ln2/4)eps={mean - lin:.4e} "
        f"predicted_remainder=5(ln2)^2eps^2/12={pred_rem:.4e} (printed, not asserted)"
    )
    assert err <= DOSE_BAR, (direction, err)


@pytest.mark.parametrize("T0", list(T0S))
@pytest.mark.parametrize("delta", DELTAS)
def test_e_small_slope(T0, delta):
    r = T0S[T0]
    T1 = T0 + delta
    got = cem43(np.array([0.0, DT]), np.array([T0, T1]))
    ref, x = _taylor(T0, T1, DT, r)
    err = _rel(got, ref)
    print(
        f"{P} E T0={T0:g} delta={delta:g} (float dT={T1 - T0:.6e}) x={x:.4e} "
        f"trunc~x^3/24={abs(x) ** 3 / 24:.1e} op={got:.15e} taylor={ref:.15e} rel={err:.3e}"
    )
    assert err <= DOSE_BAR, (T0, delta, err)


@pytest.mark.parametrize("T", list(T0S))
def test_z_constant_segment(T):
    r = T0S[T]
    got = cem43(np.array([0.0, DT]), np.array([T, T]))
    ref = DT * r ** (43.0 - T)
    err = _rel(got, ref)
    print(f"{P} Z T={T:g} op={got:.15e} ref={ref:.15e} rel={err:.3e}")
    assert err <= DOSE_BAR, (T, err)


def test_kc_r_swap_printed():
    """(Kc) predicted, printed: swap moves the mean by ≈ +(ln 2/2)ε."""
    pred = (LN2 / 2.0) * EPS
    for a, b in ((43.0 - EPS, 43.0 + EPS), (43.0 + EPS, 43.0 - EPS)):
        t, T = np.array([0.0, TAU]), np.array([a, b])
        base = cem43(t, T)
        swap = cem43(t, T, r_above=0.25, r_below=0.5)
        move = (swap - base) / base
        print(
            f"{P} Kc {a:g}->{b:g} swap_rel_move={move:+.6e} predicted=+{pred:.6e} "
            f"ratio={move / pred:.6f} miss_vs_ramp_dose={_rel(swap, ramp_dose(a, b, TAU)):.3e} "
            f"(predicted, printed, not asserted)"
        )


def test_ec_naive_quotient_printed():
    """(Ec) predicted, printed: naive quotient misses Taylor by ≈ 2e-7 at δ = 6e-10."""
    delta = 6e-10
    for T0, r in T0S.items():
        T1 = T0 + delta
        ref, _ = _taylor(T0, T1, DT, r)
        naive = ramp_dose(T0, T1, DT)
        op = cem43(np.array([0.0, DT]), np.array([T0, T1]))
        miss = _rel(naive, ref)
        print(
            f"{P} Ec T0={T0:g} delta={delta:g} naive_miss={miss:.3e} "
            f"(predicted ~2e-7 at T0=41) ratio_to_bar={miss / DOSE_BAR:.1f} "
            f"op_miss={_rel(op, ref):.3e} (predicted, printed, not asserted)"
        )
