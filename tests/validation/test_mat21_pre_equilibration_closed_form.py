"""`MAT-21` step 24 — the pre-equilibration bias against the wait, in closed form (§9 item
712, §10 lane 3 → chain P, P51; brief `docs/planning/chunks/MAT-21.md` "§9 item 712").

Numpy.  Sets no band and moves none.  The problem is linear, so a phantom that waited t_w
from T∞ + g·z and then met the source is θ(t; t_w) = θ_drift(t + t_w) + θ_sourced(t).  The
bias at a point is b(t_w) = θ_drift(t* + t_w) / θ_sourced(t*), on t_w ∈ {0, 0.5, 1, 2} / λ₁,
λ₁ the slowest **odd-in-z** mode's rate: the drift is odd in z, so step 0's slowest mode is
not its slowest.  g = ``G_GRAD`` is the fixture's choice, not a measurement.

Five of P37's six points lie in z = 0, where the odd drift is exactly zero, so the ratio
holds at P5 only and the decay gates read P5 and ``ADJ_PT``.

Asserted: (Sym0) the drift at P0–P4 ≤ 1e-13 g C; (Sup) the two limits and the mode-by-mode
time shift ≤ ``B1_BAND``; (Dec∞) on the first odd mode alone, the ratio over {1, 2} / λ₁
equals e⁻¹ to ``B1_BAND``.  Negative controls: the Dirichlet-z form misses the drift at P5
by > 1e-3 g C; the lowest even z root as λ₁ misses (Dec∞) by > ``B1_BAND``.

(Dec), the weekly's gate on the full series, is a strict xfail per the brief's
``PREDICTED RED:`` line, with its range asserted; the measured R is the record.

Not claimed: a recommended wait, a FEM run (P52), F2182's own wording.
"""

from __future__ import annotations

import numpy as np
import pytest

from fem_em_solver.utils.analytical import (
    robin_box_initial_gradient_temperature,
    robin_box_transient_temperature,
    robin_slab_eigenvalues,
    robin_slab_odd_eigenvalues,
    robin_slab_odd_linear_coefficients,
)

from tests.validation.test_imp4_pec_cylinder_closed_form import B1_BAND
from tests.validation.test_mat15_finite_cylinder_closed_form import H, K, RHO_C, T_INF
from tests.validation.test_mat21_box_closed_form import (
    A_BOX,
    B_BOX,
    C_BOX,
    HALF,
    KW,
    SERIES_TERMS,
)
from tests.validation.test_mat21_box_probe_points_fem import PROBE_LABELS, PROBE_PTS
from tests.validation.test_mat21_initial_gradient_closed_form import (
    ALPHA,
    G_GRAD,
    GC,
    N_TERMS,
    T_STAR,
    _even,
    _theta,
)

P5 = PROBE_PTS[5]
ADJ_PT = (0.0, 0.0, 0.9 * C_BOX)
GATE_PTS = ((PROBE_LABELS[5], P5), ("ADJ_PT", ADJ_PT))
SYM0_BAR = 1e-13  # the brief's (Sym0), of g C
DZ_FLOOR = 1e-3  # the brief's control floor, of g C (as `MAT-21` step 19's)
DEC_RANGE = (1e-4, 3e-1)  # PREDICTED RED: (Dec) at P5 and ADJ_PT range [1e-4, 3e-1]
_NO_SOURCE = {**KW, "q": 0.0}


def _say(msg):
    print(msg, flush=True)


def _lam(gz, gx=None, gy=None):
    gx = robin_slab_eigenvalues(H * A_BOX / K, 1)[0] if gx is None else gx
    gy = robin_slab_eigenvalues(H * B_BOX / K, 1)[0] if gy is None else gy
    return ALPHA * ((gx / A_BOX) ** 2 + (gy / B_BOX) ** 2 + (gz / C_BOX) ** 2)


LAM_ODD = _lam(robin_slab_odd_eigenvalues(H * C_BOX / K, 1)[0])
LAM_EVEN = _lam(robin_slab_eigenvalues(H * C_BOX / K, 1)[0])
WAITS = tuple(c / LAM_ODD for c in (0.0, 0.5, 1.0, 2.0))


def _drift(p, t, **kw):
    return float(_theta(*p, t, **kw)) - T_INF


def _sourced(p, t, kw=KW):
    return float(robin_box_transient_temperature(*p, t, half_widths=HALF, rho_c=RHO_C,
                                                 t_inf=T_INF, n_terms=SERIES_TERMS, **kw)) - T_INF


def _combined(p, t, t_w, kw=KW, g=G_GRAD):
    drift = float(robin_box_initial_gradient_temperature(
        *p, t + t_w, half_widths=HALF, k=K, h=H, rho_c=RHO_C, g=g, t_inf=T_INF,
        n_terms=N_TERMS)) - T_INF
    return drift + _sourced(p, t, kw)


def _advanced(p, t, t_w):
    """The drift at t from amplitudes advanced by e^{−λ t_w}, mode by mode."""
    gx, cx = _even(A_BOX, N_TERMS)
    gy, cy = _even(B_BOX, N_TERMS)
    gz = robin_slab_odd_eigenvalues(H * C_BOX / K, N_TERMS)
    bz = GC * robin_slab_odd_linear_coefficients(gz)
    rate = ALPHA * (((gx / A_BOX) ** 2)[:, None, None] + ((gy / B_BOX) ** 2)[None, :, None]
                    + ((gz / C_BOX) ** 2)[None, None, :])
    amp = cx[:, None, None] * cy[None, :, None] * bz[None, None, :] * np.exp(-rate * t_w)
    x, y, z = p
    return float(np.einsum("lmn,lmn,l,m,n->", amp, np.exp(-rate * t), np.cos(gx * x / A_BOX),
                           np.cos(gy * y / B_BOX), np.sin(gz * z / C_BOX)))


def _ratio_r(p, w1, w2, **kw):
    """b(w2) / b(w1): the source term cancels, so it is the drift's own ratio."""
    return _drift(p, T_STAR + w2, **kw) / _drift(p, T_STAR + w1, **kw)


def test_rates_and_sym0():
    _say(f"MAT21-24 lambda_1 (odd) {LAM_ODD:.6e} /s, tau {1 / LAM_ODD:.1f} s; lambda (even, "
         f"step 0's) {LAM_EVEN:.6e} /s; lambda_1/alpha {LAM_ODD / ALPHA:.2f} /m^2; t* "
         f"{T_STAR:.2f} s; waits {[round(w, 2) for w in WAITS]} s; g {G_GRAD} K/m (fixture)")
    assert LAM_ODD > LAM_EVEN
    worst = 0.0
    for p in PROBE_PTS[:5]:
        for w in WAITS:
            worst = max(worst, abs(_drift(p, T_STAR + w)) / GC)
    _say(f"MAT21-24 (Sym0) P0-P4 max |theta_drift|/(g C) {worst:.3e} (bar {SYM0_BAR:.0e}, ASSERTED)")
    assert worst <= SYM0_BAR, worst


def test_sup_limits_and_time_shift():
    w_lim = w_shift = 0.0
    for _name, p in GATE_PTS:
        for w in WAITS:
            d, s = _drift(p, T_STAR + w), _sourced(p, T_STAR)
            w_lim = max(w_lim, abs(_combined(p, T_STAR, w, kw=_NO_SOURCE) - d) / abs(d),
                        abs(_combined(p, T_STAR, w, g=0.0) - s) / abs(s))
            w_shift = max(w_shift, abs(_advanced(p, T_STAR, w) - d) / abs(d))
    _say(f"MAT21-24 (Sup) limits max rel {w_lim:.3e}; time shift max rel {w_shift:.3e} (band "
         f"{B1_BAND:g}, ASSERTED)")
    assert w_lim <= B1_BAND, w_lim
    assert w_shift <= B1_BAND, w_shift


def test_dec_inf_first_odd_mode_and_wrong_lambda():
    w1, w2 = 1.0 / LAM_ODD, 2.0 / LAM_ODD
    e1, e2 = 1.0 / LAM_EVEN, 2.0 / LAM_EVEN
    pred = abs(np.expm1(-(LAM_ODD - LAM_EVEN) / LAM_EVEN))  # waits on the even rate
    for name, p in GATE_PTS:
        r = abs(_ratio_r(p, w1, w2, n_terms=1) * np.e - 1.0)
        c = abs(_ratio_r(p, e1, e2, n_terms=1) * np.e - 1.0)
        _say(f"MAT21-24 (Dec-inf) {name}: |R| {r:.3e} (band {B1_BAND:g}, ASSERTED); wrong lambda "
             f"(even root) miss {c:.4e} (asserted > {B1_BAND:g}; predicted scale {pred:.4e})")
        assert r <= B1_BAND, (name, r)
        assert c > B1_BAND, (name, c)


def _dec_full():
    w1, w2 = 1.0 / LAM_ODD, 2.0 / LAM_ODD
    return {name: abs(_ratio_r(p, w1, w2) * np.e - 1.0) for name, p in GATE_PTS}


def test_dec_inside_predicted_range():
    """The brief's stop condition: (Dec) outside its range stops the item."""
    for name, r in _dec_full().items():
        _say(f"MAT21-24 (Dec) {name}: |R| {r:.4e}; PREDICTED RED range [{DEC_RANGE[0]:.0e}, "
             f"{DEC_RANGE[1]:.0e}] (ASSERTED); band {B1_BAND:g}")
        assert DEC_RANGE[0] <= r <= DEC_RANGE[1], (name, r)


@pytest.mark.xfail(strict=True, reason="PREDICTED RED: (Dec) at P5 and ADJ_PT range [1e-4, 3e-1] "
                                       "(docs/planning/chunks/MAT-21.md, §9 item 712)")
def test_dec_full_series_as_the_weekly_wrote_it():
    r = _dec_full()
    assert max(r.values()) <= B1_BAND, r


def test_negative_control_dirichlet_z():
    rz = np.arange(1, N_TERMS + 1) * np.pi
    rob, dz = _drift(P5, T_STAR), _drift(P5, T_STAR, roots_z=rz)
    miss = abs(rob - dz) / GC
    _say(f"MAT21-24 control (Dz) at P5, t*: Robin {rob:.6f} K, Dirichlet {dz:.3e} K, miss/(g C) "
         f"{miss:.4e} (asserted > {DZ_FLOOR:.0e})")
    assert miss > DZ_FLOOR, miss


def test_printed_bias_table():
    """The record (no band); asserts only that it was formed."""
    _say("MAT21-24 TABLE point | t_w lambda_1 | t_w (s) | theta_drift K | theta_sourced K | bias")
    for lab, p in zip(PROBE_LABELS, PROBE_PTS):
        s = _sourced(p, T_STAR)
        for c, w in zip((0.0, 0.5, 1.0, 2.0), WAITS):
            d = _drift(p, T_STAR + w)
            _say(f"MAT21-24 TABLE {lab} | {c} | {w:.1f} | {d:+.6e} | {s:.6f} | {d / s:+.6e}")
        assert s > 0.0
