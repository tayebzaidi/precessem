"""`MAT-21` step 34 — the pre-equilibration drift inside the walled box, in closed form
(§9 item 788, §10 P62; brief `docs/planning/chunks/MAT-21.md` "§9 item 788").

Numpy.  Sets no band and moves none.  θ(t; t_w) = θ_drift(t + t_w) + θ_sourced(t), as
step 24 (P51, `test_mat21_pre_equilibration_closed_form.py`), with the wall in both parts:

- **h_eff drift.**  The walled (two-layer) series has even modes only, and the drift from
  θ₀ = g·z is odd in z, so the drift is the Robin box at h = h_eff
  (`wall_effective_heat_transfer`): the wall as a quasi-steady (massless) resistance.
  The wall's heat capacity is **not** in the drift.  g = ``G_GRAD`` is the fixture's choice.
- **Sourced rise.**  The full two-layer even series (`two_layer_slab_transient_temperature`,
  item 677's ``WALL`` with ``D_WALL``, ``K_W``, ``RHO_C_W`` — the PMMA values are nominal,
  **citation owed**, W3), 400 terms, in z.

λ₁(h_eff) is P51's ``_lam`` with the first odd z root at Bi = h_eff C / K **and** the x, y
roots at h_eff (``_lam``'s gx / gy arguments): the h_eff drift carries h_eff on every face,
so its slowest odd mode decays at that rate.  Waits t_w ∈ {0, 0.5, 1, 2} / λ₁(h_eff),
t* = 0.2 T_C.

Geometry, disclosed.  The h_eff box is walled on all six faces; item 682's walled FEM box
(and the two-layer series) is 1D-reduced, walled on ±z with the side faces insulated.  At
the gate points (x = y = 0) the two drifts differ only by the side faces' reach; the
difference is printed (no band) beside the gates as a record for P63.

Asserted: (Sup) the q = 0 and g = 0 limits and the mode-by-mode time shift ≤ ``B1_BAND``;
(Dec∞) on the h_eff slowest odd mode alone (``n_terms=1``) |R e − 1| ≤ ``B1_BAND`` at P5 and
``ADJ_PT``, the even-root λ control > ``B1_BAND``; (Sym0) the drift at z = 0 ≤ ``SYM0_BAR``
g C.  Control (asserted): the wall-free series (P51's ``_combined`` at H) misses the walled
bias at z = C by > 1e-3 relative at every wait.  Printed (🟡 record): the bias at t*
relative to the sourced rise, per wait and height.

Not claimed: the wall's heat capacity in the drift, any FEM (P63), a recommended wait.
"""

from __future__ import annotations

import numpy as np

from fem_em_solver.utils.analytical import (
    robin_box_initial_gradient_temperature,
    robin_slab_eigenvalues,
    robin_slab_odd_eigenvalues,
    robin_slab_odd_linear_coefficients,
    two_layer_slab_transient_temperature,
    wall_effective_heat_transfer,
)

from tests.validation.test_imp4_pec_cylinder_closed_form import B1_BAND
from tests.validation.test_mat15_finite_cylinder_closed_form import H, K, Q, RHO_C, T_INF
from tests.validation.test_mat21_box_closed_form import A_BOX, B_BOX, C_BOX, HALF, T_C
from tests.validation.test_mat21_initial_gradient_closed_form import ALPHA, G_GRAD, GC, N_TERMS
from tests.validation.test_mat21_pre_equilibration_closed_form import (
    ADJ_PT,
    P5,
    SYM0_BAR,
    _combined,
    _lam,
)
from tests.validation.test_mat21_pre_equilibration_closed_form import (
    _drift as _drift_free,
)
from tests.validation.test_mat21_pre_equilibration_closed_form import (
    _sourced as _sourced_free,
)
from tests.validation.test_mat21_wall_closed_form import D_WALL, K_W, RHO_C_W
from tests.validation.test_mat21_wall_transient_closed_form import N_SERIES, WALL

H_EFF = wall_effective_heat_transfer(H, D_WALL, K_W)
T_STAR = 0.2 * T_C
P0 = (0.0, 0.0, 0.0)
GATE_PTS = (("P5 top face", P5), ("ADJ_PT", ADJ_PT))
TABLE_PTS = (("P0 centre", P0), ("P5 top face", P5), ("ADJ_PT", ADJ_PT))
CTRL_FLOOR = 1e-3  # the brief's control floor, relative


def _say(msg):
    print(msg, flush=True)


def _side(L):
    return robin_slab_eigenvalues(H_EFF * L / K, 1)[0]


GX, GY = _side(A_BOX), _side(B_BOX)
GZ_ODD = robin_slab_odd_eigenvalues(H_EFF * C_BOX / K, 1)[0]
GZ_EVEN = robin_slab_eigenvalues(H_EFF * C_BOX / K, 1)[0]
LAM_W = _lam(GZ_ODD, gx=GX, gy=GY)
LAM_W_EVEN = _lam(GZ_EVEN, gx=GX, gy=GY)
WAIT_C = (0.0, 0.5, 1.0, 2.0)
WAITS = tuple(c / LAM_W for c in WAIT_C)


def _drift(p, t, g=G_GRAD, n_terms=N_TERMS):
    """The h_eff drift θ (K) at p, t: the Robin box at H = h_eff from θ₀ = g·z."""
    return float(robin_box_initial_gradient_temperature(
        *p, t, half_widths=HALF, k=K, h=H_EFF, rho_c=RHO_C, g=g, t_inf=T_INF,
        n_terms=n_terms)) - T_INF


def _sourced(p, t, q=Q):
    """The two-layer even series' rise (K) at p's z."""
    return float(two_layer_slab_transient_temperature(
        p[2], t, q=q, h=H, t_inf=T_INF, n_terms=N_SERIES, **WALL)[0]) - T_INF


def _walled(p, t, t_w, q=Q, g=G_GRAD):
    return _drift(p, t + t_w, g=g) + _sourced(p, t, q=q)


def _even_h_eff(L, n):
    gam = robin_slab_eigenvalues(H_EFF * L / K, n)
    c = 2.0 * np.sin(gam) / (gam + np.sin(gam) * np.cos(gam))
    return gam, c


def _advanced(p, t, t_w):
    """The h_eff drift at t from amplitudes advanced by e^{−λ t_w}, mode by mode."""
    gx, cx = _even_h_eff(A_BOX, N_TERMS)
    gy, cy = _even_h_eff(B_BOX, N_TERMS)
    gz = robin_slab_odd_eigenvalues(H_EFF * C_BOX / K, N_TERMS)
    bz = GC * robin_slab_odd_linear_coefficients(gz)
    rate = ALPHA * (((gx / A_BOX) ** 2)[:, None, None] + ((gy / B_BOX) ** 2)[None, :, None]
                    + ((gz / C_BOX) ** 2)[None, None, :])
    amp = cx[:, None, None] * cy[None, :, None] * bz[None, None, :] * np.exp(-rate * t_w)
    x, y, z = p
    return float(np.einsum("lmn,lmn,l,m,n->", amp, np.exp(-rate * t), np.cos(gx * x / A_BOX),
                           np.cos(gy * y / B_BOX), np.sin(gz * z / C_BOX)))


def _drift_1d(p, t):
    """Print-only: the z-only h_eff drift (side faces insulated, item 682's geometry)."""
    gz = robin_slab_odd_eigenvalues(H_EFF * C_BOX / K, N_TERMS)
    bz = GC * robin_slab_odd_linear_coefficients(gz)
    return float(np.sum(bz * np.exp(-ALPHA * (gz / C_BOX) ** 2 * t) * np.sin(gz * p[2] / C_BOX)))


def test_rates_and_sym0():
    _say(f"MAT21-34 h_eff {H_EFF:.6f} W/m^2/K (h {H}, d {D_WALL}, k_w {K_W} citation owed, "
         f"rho_c_w {RHO_C_W:.4e} citation owed); lambda_1(h_eff, odd) {LAM_W:.6e} /s, tau "
         f"{1 / LAM_W:.1f} s; even-root lambda {LAM_W_EVEN:.6e} /s; t* {T_STAR:.2f} s; waits "
         f"{[round(w, 2) for w in WAITS]} s; g {G_GRAD} K/m (fixture)")
    assert LAM_W > LAM_W_EVEN
    worst = 0.0
    for w in WAITS:
        worst = max(worst, abs(_drift(P0, T_STAR + w)) / GC)
    _say(f"MAT21-34 (Sym0) z = 0 max |theta_drift|/(g C) {worst:.3e} (bar {SYM0_BAR:.0e}, "
         f"ASSERTED)")
    assert worst <= SYM0_BAR, worst


def test_sup_limits_and_time_shift():
    w_lim = w_shift = 0.0
    for _name, p in GATE_PTS:
        s = _sourced(p, T_STAR)
        for w in WAITS:
            d = _drift(p, T_STAR + w)
            w_lim = max(w_lim, abs(_walled(p, T_STAR, w, q=0.0) - d) / abs(d),
                        abs(_walled(p, T_STAR, w, g=0.0) - s) / abs(s))
            w_shift = max(w_shift, abs(_advanced(p, T_STAR, w) - d) / abs(d))
    _say(f"MAT21-34 (Sup) limits max rel {w_lim:.3e}; time shift max rel {w_shift:.3e} (band "
         f"{B1_BAND:g}, ASSERTED)")
    assert w_lim <= B1_BAND, w_lim
    assert w_shift <= B1_BAND, w_shift


def test_dec_inf_slowest_odd_mode_and_wrong_lambda():
    for name, p in GATE_PTS:
        def ratio(lam):
            return (_drift(p, T_STAR + 2.0 / lam, n_terms=1)
                    / _drift(p, T_STAR + 1.0 / lam, n_terms=1))
        r = abs(ratio(LAM_W) * np.e - 1.0)
        c = abs(ratio(LAM_W_EVEN) * np.e - 1.0)
        _say(f"MAT21-34 (Dec-inf) {name}: |R e - 1| {r:.3e} (band {B1_BAND:g}, ASSERTED); "
             f"wrong lambda (even root) miss {c:.4e} (asserted > {B1_BAND:g})")
        assert r <= B1_BAND, (name, r)
        assert c > B1_BAND, (name, c)


def test_control_wall_free_series_misses_walled_bias():
    least = np.inf
    for c, w in zip(WAIT_C, WAITS):
        b_wall = _drift(P5, T_STAR + w) / _sourced(P5, T_STAR)
        b_free = (_combined(P5, T_STAR, w) - _sourced_free(P5, T_STAR)) / _sourced_free(P5, T_STAR)
        d_free = _drift_free(P5, T_STAR + w)
        miss = abs(b_free - b_wall) / abs(b_wall)
        dmiss = abs(d_free - _drift(P5, T_STAR + w)) / abs(_drift(P5, T_STAR + w))
        least = min(least, miss)
        _say(f"MAT21-34 control P5 t_w {c}/lambda_1: walled bias {b_wall:+.6e}, wall-free bias "
             f"{b_free:+.6e}, rel miss {miss:.4e}; drift alone rel miss {dmiss:.4e}")
    _say(f"MAT21-34 control least rel miss {least:.4e} (asserted > {CTRL_FLOOR:.0e})")
    assert least > CTRL_FLOOR, least


def test_printed_bias_table_and_side_face_record():
    """The record (no band); asserts only that it was formed."""
    _say("MAT21-34 TABLE point | t_w lambda_1 | t_w (s) | theta_drift K | theta_sourced K | bias")
    for lab, p in TABLE_PTS:
        s = _sourced(p, T_STAR)
        for c, w in zip(WAIT_C, WAITS):
            d = _drift(p, T_STAR + w)
            _say(f"MAT21-34 TABLE {lab} | {c} | {w:.1f} | {d:+.6e} | {s:.6f} | {d / s:+.6e}")
        assert s > 0.0
    for lab, p in GATE_PTS:
        worst = max(abs(_drift(p, T_STAR + w) - _drift_1d(p, T_STAR + w)) / GC for w in WAITS)
        _say(f"MAT21-34 SIDE {lab}: max over waits |six-face h_eff drift - z-only drift|/(g C) "
             f"{worst:.3e} (printed, no band; item 682's box insulates the sides)")
