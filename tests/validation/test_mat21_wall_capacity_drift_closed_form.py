"""`MAT-21` step 38 — the pre-equilibration drift through a wall **with heat capacity**, in
closed form (§9 item 823, §10 P67; brief `docs/planning/chunks/MAT-21.md` "§9 item 823",
bound by the "Review ruling 2026-11-05 03:00" above it).

Numpy only.  Sets no band and moves none.  Geometry: item 682's 1D two-layer slab — gel
|z| < C (K, RHO_C), acrylic C < |z| < C + d (K_W, RHO_C_W; PMMA nominal, citation owed, W3),
Robin H at |z| = C + d.  The drift from θ₀ = g·z on [0, C + d] (gel **and** wall, as P63's R1
sets it on every cell) is odd in z, so it is expanded in the **odd** composite modes, built
here from the same μ₁ / μ₂ construction as `two_layer_slab_mode` (whose family is even):

    gel  X = sin(μ₁z);   wall (s = |z| − C)  X = sign(z)·[A cos(μ₂s) + B̃ sin(μ₂s)/μ₂],
    A = sin(μ₁C),  B̃ = (kμ₁/k_w) cos(μ₁C),  roots of F = k_w X′(C+d) + h X(C+d).

Item 788's `_drift_1d` is the massless-wall (h_eff) limit; this family keeps ρc_w.

Asserted: anchors (Lim-w) d = 0 → P47a's z family at Bi = H·C/K, (Lim-c) ρc_w = 0 →
`_drift_1d`, (Bq) N vs 2N, (Sup) q = 0 / g = 0 limits and the mode-by-mode time shift, all at
``B1_BAND``; preconditions (Root-odd) (Int-odd) (Orth-odd) (IC-odd) at P48a′'s bands; the
control (the nominal composite misses `_drift_1d` and the d = 0 form by > ``B1_BAND``).
Printed only (never asserted): the prediction Δ_cap against P63's wall-capacity term
2.173640e-02 K, reading (a) against P63's FEM drift, the tables.

Not claimed: the six-wall box, any FEM (P68 / P69), a recommended wait, a cited PMMA ρc.
"""

from __future__ import annotations

from functools import lru_cache

import numpy as np

from fem_em_solver.utils.analytical import (
    _sin2_over,
    _sinc_x,
    _two_layer_slab_mu,
    robin_slab_odd_eigenvalues,
    robin_slab_odd_linear_coefficients,
    two_layer_slab_eigenvalues,
)

from tests.validation.test_imp4_pec_cylinder_closed_form import B1_BAND
from tests.validation.test_mat13_transient_perfusion_sphere import INTERIOR_BAND
from tests.validation.test_mat15_finite_cylinder_closed_form import H, K, RHO_C
from tests.validation.test_mat21_box_closed_form import C_BOX, T_C
from tests.validation.test_mat21_initial_gradient_closed_form import ALPHA, G_GRAD, GC, N_TERMS
from tests.validation.test_mat21_pre_equilibration_closed_form import ADJ_PT, P5
from tests.validation.test_mat21_wall_closed_form import D_WALL, K_W, RHO_C_W
from tests.validation.test_mat21_wall_pre_equilibration_closed_form import (
    CTRL_FLOOR,
    H_EFF,
    T_STAR,
    WAIT_C,
    _drift_1d,
    _sourced,
)
from tests.validation.test_mat21_wall_transient_closed_form import (
    INT_BAND,
    ORTH_BAND,
    ROOT_BAND,
)

P0 = (0.0, 0.0, 0.0)
GATE_PTS = (("P5", P5), ("ADJ_PT", ADJ_PT))
TABLE_PTS = (("P0", P0), ("P5", P5), ("ADJ_PT", ADJ_PT))
MASSLESS = 1e-9  # test_mat21_wall_transient_fem.py:75, the (K0) scale; printed beside (Lim-c)
N_ROOT = 60      # item 673: (Root) at 200 roots fell below float64 resolution

# P63's FEM record, hard-coded (`20261104T060527Z_MAT-21-P63-W1b.log:76–80`, columns 2 and 4):
# (label, point, t_w is T_C, FEM drift R_k − R0 K, wall-capacity term FEM − _drift_1d K)
P63_ROWS = (
    ("P5 t_w=0", P5, False, 1.828983e-01, 2.173640e-02),
    ("ADJ_PT t_w=0", ADJ_PT, False, 1.828894e-01, 1.544980e-02),
    ("P5 t_w=T_C", P5, True, 6.721520e-03, 1.803301e-03),
    ("ADJ_PT t_w=T_C", ADJ_PT, True, 6.729312e-03, 1.617858e-03),
)


def _say(msg):
    print(msg, flush=True)


# ------------------------------------------------------------------ the odd family
def _F_odd(x, *, d=D_WALL, rho_c_w=RHO_C_W, h=H, return_scale=False):
    """Robin residual of the odd mode at |z| = C + d, x = √λ; F(0) = 0 (trivial zero)."""
    C = float(C_BOX)
    m1, m2 = _two_layer_slab_mu(x, K, K_W, RHO_C, rho_c_w)
    A = np.sin(m1 * C)
    Bt = (K / K_W) * m1 * np.cos(m1 * C)
    c2 = np.cos(m2 * d)
    sd = d * _sinc_x(m2 * d)
    t1 = h * A * c2
    t2 = h * Bt * sd
    t3 = -K_W * A * m2 * m2 * sd
    t4 = K_W * Bt * c2
    F = t1 + t2 + t3 + t4
    if return_scale:
        return F, np.abs(t1) + np.abs(t2) + np.abs(t3) + np.abs(t4)
    return F


def _odd_eigs_raw(n, *, d=D_WALL, rho_c_w=RHO_C_W, refine=64):
    """Copy of `two_layer_slab_eigenvalues` on `_F_odd`; the scan starts at x = step."""
    from scipy.optimize import brentq

    kw = dict(d=d, rho_c_w=rho_c_w)
    lam_len = C_BOX * np.sqrt(RHO_C / K) + d * np.sqrt(rho_c_w / K_W)
    step = np.pi / (float(refine) * lam_len)
    n = int(n)
    roots = []
    x_lo = step
    while len(roots) < n:
        need = n - len(roots)
        x_hi = x_lo + (need + 8) * 1.25 * np.pi / lam_len
        xs = np.arange(x_lo, x_hi + step, step)
        fs = _F_odd(xs, **kw)
        idx = np.nonzero(np.sign(fs[:-1]) * np.sign(fs[1:]) < 0)[0]
        for i in idx:
            x0 = brentq(lambda x: float(_F_odd(x, **kw)), xs[i], xs[i + 1],
                        xtol=1e-300, rtol=4 * np.finfo(float).eps, maxiter=500)
            cand = x0 + np.spacing(x0) * np.arange(-16, 17)
            roots.append(float(cand[np.argmin(np.abs(_F_odd(cand, **kw)))]))
            if len(roots) == n:
                break
        x_lo = xs[-1]
    x = np.asarray(roots[:n])
    return x * x


@lru_cache(maxsize=None)
def _odd_eigs(n, d=D_WALL, rho_c_w=RHO_C_W, refine=64):
    return _odd_eigs_raw(n, d=d, rho_c_w=rho_c_w, refine=refine)


def _mode_odd(z, lam, *, d=D_WALL, rho_c_w=RHO_C_W, deriv=0):
    C = float(C_BOX)
    z = np.atleast_1d(np.asarray(z, dtype=float))[:, None]
    lam = np.atleast_1d(np.asarray(lam, dtype=float))[None, :]
    m1, m2 = _two_layer_slab_mu(np.sqrt(lam), K, K_W, RHO_C, rho_c_w)
    a = np.abs(z)
    gel = a <= C
    s = np.where(gel, 0.0, a - C)
    A = np.sin(m1 * C)
    Bt = (K / K_W) * m1 * np.cos(m1 * C)
    cw, sw = np.cos(m2 * s), s * _sinc_x(m2 * s)
    sg = np.sign(z)
    if deriv == 0:
        return np.where(gel, np.sin(m1 * z), sg * (A * cw + Bt * sw))
    if deriv == 1:
        return np.where(gel, m1 * np.cos(m1 * z), -A * m2 * m2 * sw + Bt * cw)
    raise ValueError("deriv must be 0 or 1")


def _norm_closed(lam, *, d=D_WALL, rho_c_w=RHO_C_W):
    C = float(C_BOX)
    m1, m2 = _two_layer_slab_mu(np.sqrt(lam), K, K_W, RHO_C, rho_c_w)
    A = np.sin(m1 * C)
    Bt = (K / K_W) * m1 * np.cos(m1 * C)
    x2 = m2 * d
    cc = 0.5 * d + 0.5 * d * _sinc_x(2 * x2)
    ss = _sin2_over(x2, d)
    cs = 0.5 * d * d * _sinc_x(x2) ** 2
    n_gel = RHO_C * 0.5 * C * (1.0 - _sinc_x(2 * m1 * C))
    return n_gel + rho_c_w * (A * A * cc + Bt * Bt * ss + 2.0 * A * Bt * cs)


def _gl_wall(d, n_panels=40, order=12):
    """Composite 12-point GL on its own d-parameterised panels of [C, C + d]."""
    xg, wg = np.polynomial.legendre.leggauss(order)
    e = np.linspace(C_BOX, C_BOX + d, n_panels + 1)
    lo, hi = e[:-1, None], e[1:, None]
    return (0.5 * (hi - lo) * xg + 0.5 * (hi + lo)).ravel(), (0.5 * (hi - lo) * wg).ravel()


def _gl_layers(d, n_gel, n_wall, order=12):
    xg, wg = np.polynomial.legendre.leggauss(order)
    zs, ws = [], []
    for a, b, m in ((0.0, C_BOX, n_gel), (C_BOX, C_BOX + d, n_wall)):
        e = np.linspace(a, b, m + 1)
        lo, hi = e[:-1, None], e[1:, None]
        zs.append((0.5 * (hi - lo) * xg + 0.5 * (hi + lo)).ravel())
        ws.append((0.5 * (hi - lo) * wg).ravel())
    return np.concatenate(zs), np.concatenate(ws)


def _proj_gel(lam, g=G_GRAD):
    m1 = np.sqrt(lam * RHO_C / K)
    return RHO_C * g * (np.sin(m1 * C_BOX) - m1 * C_BOX * np.cos(m1 * C_BOX)) / m1**2


def _proj_wall(lam, d, rho_c_w, g=G_GRAD):
    if d == 0.0 or rho_c_w == 0.0:
        return np.zeros_like(lam)
    zw, ww = _gl_wall(d)
    X = _mode_odd(zw, lam, d=d, rho_c_w=rho_c_w)
    return rho_c_w * g * ((zw * ww) @ X)


@lru_cache(maxsize=None)
def _coeffs(n, d=D_WALL, rho_c_w=RHO_C_W, gel_only=False, g=G_GRAD):
    lam = _odd_eigs(n, d, rho_c_w)
    p = _proj_gel(lam, g)
    if not gel_only:
        p = p + _proj_wall(lam, d, rho_c_w, g)
    return lam, p / _norm_closed(lam, d=d, rho_c_w=rho_c_w)


def _drift_cap(p, t, n=N_TERMS, d=D_WALL, rho_c_w=RHO_C_W, gel_only=False, g=G_GRAD):
    """The capacitive-wall drift θ − T∞ (K) at p's z, t, from θ₀ = g·z on [0, C + d]."""
    lam, c = _coeffs(int(n), float(d), float(rho_c_w), bool(gel_only), float(g))
    X = _mode_odd(p[2], lam, d=d, rho_c_w=rho_c_w)[0]
    return float(np.sum(c * np.exp(-lam * t) * X))


def _drift_cap_z(z, t, n=N_TERMS, d=D_WALL, rho_c_w=RHO_C_W):
    """Vectorised over z: mode matrix @ amplitudes (items 835 / 836)."""
    lam, c = _coeffs(int(n), float(d), float(rho_c_w), False, float(G_GRAD))
    return _mode_odd(z, lam, d=d, rho_c_w=rho_c_w) @ (c * np.exp(-lam * float(t)))


LAM1_W = float(_odd_eigs(1)[0])
WAITS = tuple(c / LAM1_W for c in WAIT_C) + (T_C,)
WAIT_LABELS = tuple(f"{c}/LAM1_W" for c in WAIT_C) + ("T_C",)


def _p47a_z(p, t):
    """P47a's z family at H (Bi = H·C/K): GC Σ bₙ sin(γₙz/C) e^{−α(γₙ/C)²t}."""
    gz = robin_slab_odd_eigenvalues(H * C_BOX / K, N_TERMS)
    bz = GC * robin_slab_odd_linear_coefficients(gz)
    return float(np.sum(bz * np.exp(-ALPHA * (gz / C_BOX) ** 2 * t) * np.sin(gz * p[2] / C_BOX)))


def _rel(a, b):
    return abs(a - b) / abs(b)


# ------------------------------------------------------------------ preconditions
def test_preconditions_root_int_orth_ic():
    gates = {}
    # (Root-odd)
    lam = _odd_eigs(N_ROOT)
    F, S = _F_odd(np.sqrt(lam), return_scale=True)
    froot = float((np.abs(F) / S).max())
    lam_f = _odd_eigs_raw(N_ROOT, refine=256)
    agree = float(np.abs(lam_f - lam).max() / lam.max())
    incr = bool(np.all(np.diff(lam) > 0)) and lam[0] > 0
    _say(f"MAT21-38 (Root-odd) n <= {N_ROOT}: max |F|/scale {froot:.3e}; refine-256 rescan "
         f"{agree:.3e} (band {ROOT_BAND:g}); strictly increasing {incr}")
    gates["Root"] = froot <= ROOT_BAND and agree <= ROOT_BAND and incr
    # (Int-odd)
    m1 = np.sqrt(lam * RHO_C / K)
    zo = np.nextafter(C_BOX, np.inf)
    Xg, Xw = _mode_odd(C_BOX, lam)[0], _mode_odd(zo, lam)[0]
    dg, dw = _mode_odd(C_BOX, lam, deriv=1)[0], _mode_odd(zo, lam, deriv=1)[0]
    t_jump = float(np.abs(Xg - Xw).max())
    f_scale = K * m1
    f_jump = float((np.abs(K * dg - K_W * dw) / f_scale).max())
    L = C_BOX + D_WALL
    Xe, de = _mode_odd(L, lam)[0], _mode_odd(L, lam, deriv=1)[0]
    robin = float((np.abs(K_W * de + H * Xe) / (np.abs(K_W * de) + np.abs(H * Xe) + f_scale)).max())
    _say(f"MAT21-38 (Int-odd) {N_ROOT} modes: |X(C-) - X(C+)| {t_jump:.3e}; flux jump/(k mu1) "
         f"{f_jump:.3e}; outer Robin/scale {robin:.3e} (band {INT_BAND:g})")
    gates["Int"] = max(t_jump, f_jump, robin) <= INT_BAND
    # (Orth-odd)
    lam20 = _odd_eigs(20)
    z, w = _gl_layers(D_WALL, 200, 40)
    rc = np.where(z <= C_BOX, RHO_C, RHO_C_W)
    X = _mode_odd(z, lam20)
    Gm = (X * (rc * w)[:, None]).T @ X
    dgn = np.diag(Gm)
    off = float((np.abs(Gm - np.diag(dgn)) / np.sqrt(np.outer(dgn, dgn))).max())
    n_err = float(np.abs(_norm_closed(lam20) - dgn).max() / dgn.max())
    pg = (X * (rc * w * G_GRAD * z)[:, None]).sum(0)
    p_err = float(np.abs(pg - (_proj_gel(lam20) + _proj_wall(lam20, D_WALL, RHO_C_W))).max()
                  / np.abs(pg).max())
    _say(f"MAT21-38 (Orth-odd) 20 modes, 12-pt GL 200 gel + 40 wall panels: max offdiag "
         f"{off:.3e}; N closed vs GL {n_err:.3e} (band {ORTH_BAND:g}); projection vs GL "
         f"{p_err:.3e} (printed)")
    gates["Orth"] = off <= ORTH_BAND and n_err <= ORTH_BAND
    # (IC-odd)
    z, w = _gl_layers(D_WALL, 1000, 100, order=8)
    rc = np.where(z <= C_BOX, RHO_C, RHO_C_W)
    th0 = G_GRAD * z
    nrm = np.sqrt((rc * w * th0**2).sum())
    errs = {}
    for n in (25, 100, 400):
        rec = _drift_cap_z(z, 0.0, n=n)
        errs[n] = float(np.sqrt((rc * w * (rec - th0) ** 2).sum()) / nrm)
    ratios = [errs[4 * n] / errs[n] for n in (25, 100)]
    _say("MAT21-38 (IC-odd) rho c-weighted L2 of g z on [0, C+d]: "
         + ", ".join(f"N={n}: {e:.3e}" for n, e in errs.items())
         + f"; ratios 4N/N {ratios[0]:.4f}, {ratios[1]:.4f} (bar 0.5)")
    gates["IC"] = ratios[0] <= 0.5 and ratios[1] <= 0.5
    _say(f"MAT21-38 preconditions gates: {gates}")
    assert all(gates.values()), gates


# ------------------------------------------------------------------ anchors
def test_anchors_lim_w_lim_c_bq_sup():
    worst = {"Lim-w": 0.0, "Lim-c": 0.0, "Bq": 0.0, "Sup-lim": 0.0, "Sup-shift": 0.0}
    massless = 0.0
    lam, c = _coeffs(N_TERMS)
    for name, p in GATE_PTS:
        s = _sourced(p, T_STAR)
        Xp = _mode_odd(p[2], lam)[0]
        for w in WAITS:
            t = T_STAR + w
            worst["Lim-w"] = max(worst["Lim-w"], _rel(_drift_cap(p, t, d=0.0), _p47a_z(p, t)))
            ref = _drift_1d(p, t)
            worst["Lim-c"] = max(worst["Lim-c"], _rel(_drift_cap(p, t, rho_c_w=0.0), ref))
            massless = max(massless, _rel(_drift_cap(p, t, rho_c_w=MASSLESS * RHO_C_W), ref))
            dc = _drift_cap(p, t)
            worst["Bq"] = max(worst["Bq"], _rel(dc, _drift_cap(p, t, n=2 * N_TERMS)))
            th_q0 = _drift_cap(p, t) + _sourced(p, T_STAR, q=0.0)
            th_g0 = _drift_cap(p, t, g=0.0) + s
            worst["Sup-lim"] = max(worst["Sup-lim"], _rel(th_q0, dc), _rel(th_g0, s))
            adv = float(np.sum(c * np.exp(-lam * w) * np.exp(-lam * T_STAR) * Xp))
            worst["Sup-shift"] = max(worst["Sup-shift"], _rel(adv, dc))
    for k_, v in worst.items():
        _say(f"MAT21-38 anchor ({k_}) max rel over P5, ADJ_PT x {len(WAITS)} waits {v:.3e} "
             f"(band {B1_BAND:g}, ASSERTED)")
    _say(f"MAT21-38 (Lim-c) at MASSLESS*RHO_C_W instead of 0.0: max rel {massless:.3e} "
         f"(printed; predicted approx 5e-11, not asserted)")
    assert all(v <= B1_BAND for v in worst.values()), worst


# ------------------------------------------------------------------ control
def test_control_composite_misses_massless_and_d0():
    dc = _drift_cap(P5, T_STAR)
    m_1d = _rel(dc, _drift_1d(P5, T_STAR))
    m_d0 = _rel(dc, _drift_cap(P5, T_STAR, d=0.0))
    _say(f"MAT21-38 control P5 t*: composite {dc:+.6e} K; miss vs _drift_1d {m_1d:.4e}; miss "
         f"vs d = 0 form {m_d0:.4e} (both asserted > {B1_BAND:g}); predicted miss vs _drift_1d "
         f"approx 0.135 (printed beside CTRL_FLOOR {CTRL_FLOOR:g}, not asserted)")
    assert m_1d > B1_BAND and m_d0 > B1_BAND, (m_1d, m_d0)


# ------------------------------------------------------------------ printed records
def test_printed_prediction_and_tables():
    """The record (no band); the prediction is printed, never asserted."""
    for lab, p, is_tc, fem, term in P63_ROWS:
        t = T_STAR + (T_C if is_tc else 0.0)
        dc = _drift_cap(p, t)
        d1 = _drift_1d(p, t)
        delta = dc - d1
        rel = abs(delta - term) / term
        ra = abs(dc - fem) / fem
        note = " (T_C row carries item 805's even CN-switch residue 1.778e-05 K at P0)" if is_tc else ""
        _say(f"MAT21-38 PREDICTION {lab}: _drift_cap {dc:+.6e}, _drift_1d {d1:+.6e}, Delta_cap "
             f"{delta:+.6e} vs P63 term {term:.6e}: literal rel {rel:.4e} -> "
             f"{'GREEN' if rel <= INTERIOR_BAND else 'RED'} (INTERIOR_BAND {INTERIOR_BAND:g}); "
             f"reading (a) |_drift_cap - FEM {fem:.6e}|/FEM {ra:.4e} -> "
             f"{'GREEN' if ra <= INTERIOR_BAND else 'RED'}{note} (printed, never asserted)")
    _say("MAT21-38 TABLE point | t_w | t_w s | _drift_cap K | _drift_1d K | Delta_cap K | "
         "bias cap | gel-only-IC drift K | wall-IC share K")
    for lab, p in TABLE_PTS:
        s = _sourced(p, T_STAR)
        for wl, w in zip(WAIT_LABELS, WAITS):
            t = T_STAR + w
            dc, d1, dg = _drift_cap(p, t), _drift_1d(p, t), _drift_cap(p, t, gel_only=True)
            _say(f"MAT21-38 TABLE {lab} | {wl} | {w:.1f} | {dc:+.6e} | {d1:+.6e} | "
                 f"{dc - d1:+.6e} | {dc / s:+.6e} | {dg:+.6e} | {dc - dg:+.6e}")
        assert s > 0.0
    z0 = _drift_cap(P0, T_STAR)
    lam_heff = ALPHA * (robin_slab_odd_eigenvalues(H_EFF * C_BOX / K, 1)[0] / C_BOX) ** 2
    lam_even = float(two_layer_slab_eigenvalues(half_thickness=C_BOX, wall=D_WALL, k=K, k_w=K_W,
                                                rho_c=RHO_C, rho_c_w=RHO_C_W, h=H, n=1)[0])
    _say(f"MAT21-38 drift at z = 0: {z0!r} (exactly 0 by sin(0))")
    _say(f"MAT21-38 LAM1_W {LAM1_W:.6e} /s (tau {1 / LAM1_W:.1f} s); h_eff odd 1D rate "
         f"{lam_heff:.6e} /s; P48a'' even lambda_1 {lam_even:.6e} /s; t* {T_STAR:.2f} s; waits "
         f"{[round(w, 2) for w in WAITS]} s")
    assert z0 == 0.0
