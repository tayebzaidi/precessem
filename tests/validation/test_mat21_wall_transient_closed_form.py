"""`MAT-21` step 20b — the two-layer (gel + acrylic wall) slab **transient** in
closed form (§9 item 673, §10 P48a′).

The symmetric 1D slab: gel |z| < C (k, ρc, uniform q switched on at t = 0),
source-free wall C < |z| < C + d (k_w, ρ_w c_w), Robin h at |z| = C + d,
θ(0) = 0: `utils.analytical.two_layer_slab_transient_temperature`
(θ = θ_ss − Σ aₙ Xₙ e^{−λₙ t}; derivation in its docstring, composite-slab
textbook form: citation owed).  aₙ by Green's identity in closed form;
Nₙ = ∫ρc Xₙ² in closed form; the Gram matrix of (Orth) and the (IC) L2 norms
by composite Gauss–Legendre per layer (stated in each test).

Constants imported: d, k_w, ρ_w c_w from item 660's module (`K_W`, `RHO_C_W`
**nominal handbook values, citation owed**), C and T_C from step 0, gel k, ρc,
q, h, T∞ from `MAT-15`.

Times.  t\\* = 0.2 T_C (step 0's convention, T_C = C²ρc/k).  The brief's
"t_s" is not defined in it; it is read here as the slowest time constant
t_s = 1/λ₁ of the walled slab (≈ 2.163 T_C, logged 2.162790 T_C in the
attempt-2 probe), and (Lim) at 50 T_C is printed beside it.  t_s = 1/λ₁ is
ratified provisionally by `MAT-21.md` "Review ruling 2026-10-27 18:00".

(Root) is `MAT-18`'s R0a / R0b, copied (same ruling, §9 item 677): (Root-a)
|F|/scale ≤ 1e-12 on n ≤ 60; (Root-b) on all 200 roots F changes sign within
±2 ulp and the stored double is the ±2-ulp |F|-minimiser; the +8-ulp control
must fail (Root-b).

Asserted: (Root) (Int) (PDE) (Orth) (Lim) (IC: truncation ratio) (Bal),
the reductions (K0) massless wall = bare Robin slab at h_eff and (K∞) the
k_w → ∞, ρ_w c_w → 0 corner = bare slab at h, and the negative control
(bare slab at h misses the walled centre at t\\*).  Printed only: the
weekly's (ii) max |two-layer − bare(h_eff)|/rise beside ε_w, the
k_w → ∞ corner with the wall's capacity kept.
"""

from __future__ import annotations

import numpy as np

from fem_em_solver.utils.analytical import (
    _robin_slab_coefficients,
    _two_layer_slab_F,
    robin_slab_eigenvalues,
    two_layer_slab_eigenvalues,
    two_layer_slab_mode,
    two_layer_slab_mode_integrals,
    two_layer_slab_steady_temperature,
    two_layer_slab_transient_coefficients,
    two_layer_slab_transient_temperature,
    wall_effective_heat_transfer,
)
from tests.validation.test_mat15_finite_cylinder_closed_form import H, K, Q, RHO_C, T_INF
from tests.validation.test_mat21_box_closed_form import C_BOX, T_C
from tests.validation.test_mat21_wall_closed_form import D_WALL, K_W, RHO_C_W

L_OUT = C_BOX + D_WALL
T_STAR = 0.2 * T_C
N_MODES = 200             # (Root) / (Int) / (PDE) mode count, the brief's
N_SERIES = 400            # transient series default
N_BAL = 4000              # (Bal): the identity is the IC truncation of ∫ρcθ(0)
ROOT_BAND = 1.0e-12
INT_BAND = 1.0e-12
PDE_BAND = 1.0e-12
ORTH_BAND = 1.0e-10
LIM_BAND = 1.0e-12
BAL_BAND = 1.0e-10
RED_BAND = 1.0e-9
NEG_FLOOR = 1.0e-3
Z_GEL = np.linspace(-C_BOX, C_BOX, 41)
WALL = dict(half_thickness=C_BOX, wall=D_WALL, k=K, k_w=K_W, rho_c=RHO_C, rho_c_w=RHO_C_W)


def _say(msg):
    print(msg, flush=True)


def _lam(n=N_SERIES, **over):
    p = dict(WALL, **over)
    return two_layer_slab_eigenvalues(h=H, n=n, **p)


def _theta(z, t, n=N_SERIES, lam=None, **over):
    p = dict(WALL, **over)
    return two_layer_slab_transient_temperature(z, t, q=Q, h=H, t_inf=T_INF, n_terms=n,
                                                eigenvalues=lam, **p)


def _rise_ss():
    return float(two_layer_slab_steady_temperature(0.0, half_thickness=C_BOX, wall=D_WALL,
                                                   k=K, k_w=K_W, q=Q, h=H) )


def _bare_transient(z, t, h, n=N_SERIES):
    """Bare Robin slab |z| < C at film h: T∞ + (qC²/k) Σ cₙ (1 − e^{−κγₙ²t/C²}) cos(γₙ z/C)/γₙ²
    (`robin_slab_eigenvalues` + `_robin_slab_coefficients`, item 660's comparand)."""
    gam = robin_slab_eigenvalues(h * C_BOX / K, n)
    c = _robin_slab_coefficients(gam)
    z = np.atleast_1d(np.asarray(z, dtype=float))
    rate = (K / RHO_C) * gam**2 / C_BOX**2
    steady_exact = T_INF + Q * C_BOX / h + Q * (C_BOX**2 - z**2) / (2 * K)
    decay = (Q * C_BOX**2 / K) * (np.cos(np.outer(z, gam / C_BOX)) * (c / gam**2 * np.exp(-rate * t))).sum(1)
    return steady_exact - decay


def _gl_layers(n_sub_gel, n_sub_wall, order=12):
    """Composite Gauss–Legendre on [0, C] and [C, C + d]."""
    xg, wg = np.polynomial.legendre.leggauss(order)
    zs, ws = [], []
    for a, b, m in ((0.0, C_BOX, n_sub_gel), (C_BOX, L_OUT, n_sub_wall)):
        e = np.linspace(a, b, m + 1)
        lo, hi = e[:-1, None], e[1:, None]
        zs.append((0.5 * (hi - lo) * xg + 0.5 * (hi + lo)).ravel())
        ws.append((0.5 * (hi - lo) * wg).ravel())
    z = np.concatenate(zs)
    # a node exactly at C would take the gel branch; GL nodes never hit endpoints
    return z, np.concatenate(ws), np.where(z <= C_BOX, RHO_C, RHO_C_W)


# ---------------------------------------------------------------- (Root)
def test_root_residual_and_one_root_per_bracket():
    lam = _lam(N_MODES)
    x = np.sqrt(lam)
    F, S = _two_layer_slab_F(x, h=H, **WALL, return_scale=True)
    rel = np.abs(F) / S
    # one root per bracket: a 4x finer scan (refine 256) finds the same roots
    lam_f = two_layer_slab_eigenvalues(h=H, n=N_MODES, refine=256, **WALL)
    agree = np.abs(lam_f - lam).max() / lam.max()
    # and the scan step at refine 64 is below a tenth of the smallest root gap
    lam_len = C_BOX * np.sqrt(RHO_C / K) + D_WALL * np.sqrt(RHO_C_W / K_W)
    gap = np.diff(np.concatenate([[0.0], x])).min()
    step = np.pi / (64 * lam_len)
    _say(f"MAT21-20b (Root) {N_MODES} modes: lambda_1 {lam[0]:.9e} 1/s (t_s = 1/lambda_1 "
         f"{1 / lam[0]:.1f} s = {1 / (lam[0] * T_C):.6f} T_C), lambda_200 {lam[-1]:.6e}; "
         f"max |F|/scale {rel.max():.3e}; refine-256 scan max |dlam|/lam_max {agree:.3e}; "
         f"min root gap in sqrt(lam) {gap:.4e} vs scan step {step:.4e} (ratio {gap / step:.1f})")
    # probe (print only, added after W1 read 5.49e-12): F/scale one ulp either
    # side of each root, i.e. the float64 resolution of this residual
    up = np.abs(_two_layer_slab_F(np.nextafter(x, np.inf), h=H, **WALL)) / S
    dn = np.abs(_two_layer_slab_F(np.nextafter(x, -np.inf), h=H, **WALL)) / S
    sign_ok = np.sign(_two_layer_slab_F(np.nextafter(x, np.inf), h=H, **WALL)) != np.sign(
        _two_layer_slab_F(np.nextafter(x, -np.inf), h=H, **WALL))
    best = rel <= np.minimum(up, dn)
    worst_n = np.argsort(rel)[-5:][::-1]
    for i in worst_n:
        _say(f"MAT21-20b (Root probe) n={i + 1}: |F|/scale {rel[i]:.3e}; +1 ulp {up[i]:.3e}, "
             f"-1 ulp {dn[i]:.3e}; sign change across +-1 ulp {bool(sign_ok[i])}; "
             f"root is |F|-min of the 3 doubles {bool(best[i])}")
    _say(f"MAT21-20b (Root probe) all {N_MODES}: |F|-min among 3 neighbouring doubles "
         f"{int(best.sum())}/{N_MODES}; sign change across +-1 ulp {int(sign_ok.sum())}/{N_MODES}; "
         f"max over n<=60 {rel[:60].max():.3e}; min one-ulp step max(up,dn) over the 5 worst "
         f"{np.maximum(up, dn)[worst_n].min():.3e}")
    # (Root) re-posed as MAT-18's ratified R0a / R0b (MAT-21.md "Review ruling
    # 2026-10-27 18:00", item 677).  Attempt 2 (W1
    # 20261027T201000Z_MAT-21-step20b-W1.log) read max |F|/scale 5.49e-12 at
    # n = 159 against 1e-12 over all 200 roots, and its probe
    # (20261027T201037Z_MAT-21-step20b-W1-probe.log) measured one ulp of x
    # moving F/scale by >= 1.33e-12 on the five worst modes and 2.9e-11 at
    # n = 159: the band sat below float64 resolution there.  (Root-a) keeps
    # ROOT_BAND = 1e-12 on n <= 60; (Root-b) certifies all 200 roots to
    # float64 resolution at +-2 ulp (the MAT-18 precedent; never +-3).
    fres = float(rel[:60].max())
    _say(f"MAT21-20b (Root-a) max |F|/scale over n <= 60 = {fres:.3e} at n = "
         f"{int(np.argmax(rel[:60])) + 1} (band {ROOT_BAND:g}); (record) over all "
         f"{N_MODES} = {rel.max():.3e} at n = {int(np.argmax(rel)) + 1} (printed only)")

    # _nbhd / _r0b copied from tests/validation/test_mat18_two_layer_sphere_closed_form.py
    # :141-156 (MAT-18 R0b, §9 item 199), with _two_layer_sphere_F -> _two_layer_slab_F.
    def _nbhd(x0, m):
        """The 2m+1 consecutive float64 values centred on x0 (np.nextafter)."""
        x0 = np.asarray(x0, dtype=np.float64)
        lo, hi = [x0], [x0]
        for _ in range(m):
            lo.append(np.nextafter(lo[-1], np.float64(-np.inf)))
            hi.append(np.nextafter(hi[-1], np.float64(np.inf)))
        return np.stack(lo[:0:-1] + hi, axis=-1)  # (..., 2m+1), ascending

    def _r0b(xr):
        cand = _nbhd(xr, 2)
        Fc = _two_layer_slab_F(cand, h=H, **WALL)
        sign_ok = (Fc.min(axis=-1) <= 0.0) & (Fc.max(axis=-1) >= 0.0)
        a_ = np.abs(Fc)
        min_ok = a_[..., 2] <= a_.min(axis=-1)
        return sign_ok, min_ok

    sign_ok2, min_ok2 = _r0b(x)
    n_nosign = int(np.sum(~sign_ok2))
    n_nonmin = int(np.sum(~min_ok2))
    _say(f"MAT21-20b (Root-b) {N_MODES} roots: no sign change within +-2 ulp: {n_nosign} "
         f"(must be 0); stored double not the +-2-ulp |F|-minimiser: {n_nonmin} (must be 0)")
    bad_b = np.nonzero(~(sign_ok2 & min_ok2))[0]
    iw = int(bad_b[0]) if bad_b.size else int(np.argmax(rel))
    nb = _nbhd(x[iw], 4)
    _say(f"MAT21-20b (Root-b) +-4-ulp neighbourhood of n = {iw + 1} "
         f"({'first failing root' if bad_b.size else 'worst |F|/scale, record'}): "
         f"F/scale = {np.array2string(_two_layer_slab_F(nb, h=H, **WALL) / S[iw], precision=3)}")
    if bad_b.size:
        _say(f"MAT21-20b (Root-b) failing n: {(bad_b + 1)[:20].tolist()}")
    xp = x.copy()
    for _ in range(8):
        xp = np.nextafter(xp, np.float64(np.inf))
    _, min_ok_p = _r0b(xp)
    n_ctl = int(np.sum(~min_ok_p))
    _say(f"MAT21-20b (Root-b control) roots +8 ulp failing the +-2-ulp minimiser check: "
         f"{n_ctl} of {N_MODES} (asserted > 0)")
    assert fres <= ROOT_BAND
    assert n_nosign == 0 and n_nonmin == 0
    assert n_ctl > 0
    assert np.all(np.diff(lam) > 0)
    assert agree <= ROOT_BAND
    assert gap > 4 * step


# ---------------------------------------------------------------- (Int)
def test_interface_continuity_and_outer_robin_per_mode():
    lam = _lam(N_MODES)
    m1 = np.sqrt(lam * RHO_C / K)
    zi = C_BOX
    zo = np.nextafter(zi, np.inf)
    Xg = two_layer_slab_mode(zi, lam, **WALL)[0]
    Xw = two_layer_slab_mode(zo, lam, **WALL)[0]
    dg = two_layer_slab_mode(zi, lam, deriv=1, **WALL)[0]
    dw = two_layer_slab_mode(zo, lam, deriv=1, **WALL)[0]
    t_jump = np.abs(Xg - Xw).max()                          # |X| ≤ 1 in the gel
    f_scale = K * m1
    f_jump = (np.abs(K * dg - K_W * dw) / f_scale).max()
    Xe = two_layer_slab_mode(L_OUT, lam, **WALL)[0]
    de = two_layer_slab_mode(L_OUT, lam, deriv=1, **WALL)[0]
    r_scale = np.abs(K_W * de) + np.abs(H * Xe) + f_scale
    robin = (np.abs(K_W * de + H * Xe) / r_scale).max()
    _say(f"MAT21-20b (Int) {N_MODES} modes: max |X(C-) - X(C+)| {t_jump:.3e}; "
         f"max |k X'(C-) - k_w X'(C+)|/(k mu1) {f_jump:.3e}; max Robin residual at C+d "
         f"/scale {robin:.3e}")
    assert t_jump <= INT_BAND
    assert f_jump <= INT_BAND
    assert robin <= INT_BAND


# ---------------------------------------------------------------- (PDE)
def test_mode_ode_in_each_layer():
    lam = _lam(N_MODES)
    zg = np.linspace(0.0, C_BOX, 23)[1:-1]
    zw = np.linspace(C_BOX, L_OUT, 23)[1:-1]
    worst = 0.0
    for zz, kk, rc in ((zg, K, RHO_C), (zw, K_W, RHO_C_W)):
        X = two_layer_slab_mode(zz, lam, **WALL)
        X2 = two_layer_slab_mode(zz, lam, deriv=2, **WALL)
        res = np.abs(kk * X2 + rc * lam * X)
        scale = np.abs(kk * X2) + np.abs(rc * lam * X)
        scale = np.maximum(scale, rc * lam * np.abs(X).max(0))   # avoid 0/0 at nodes
        worst = max(worst, float((res / scale).max()))
    _say(f"MAT21-20b (PDE) {N_MODES} modes x 21 + 21 points: max |k X'' + rho c lambda X|/scale {worst:.3e}")
    assert worst <= PDE_BAND


# ---------------------------------------------------------------- (Orth)
def test_rho_c_weighted_orthogonality_and_closed_form_integrals():
    lam = _lam(20)
    z, w, rc = _gl_layers(200, 40)            # 12-point GL, 200 gel + 40 wall panels
    X = two_layer_slab_mode(z, lam, **WALL)
    Gm = (X * (rc * w)[:, None]).T @ X
    dg = np.diag(Gm)
    off = np.abs(Gm - np.diag(dg)) / np.sqrt(np.outer(dg, dg))
    N, I, G = two_layer_slab_mode_integrals(lam, **WALL)
    n_err = np.abs(N - dg).max() / dg.max()
    # Green's-identity coefficient vs the quadrature projection of θ_ss
    ss = two_layer_slab_steady_temperature(z, half_thickness=C_BOX, wall=D_WALL, k=K,
                                           k_w=K_W, q=Q, h=H)
    a_q = (X * (rc * w * ss)[:, None]).sum(0) / dg
    a_g = two_layer_slab_transient_coefficients(lam, q=Q, **WALL)
    a_err = np.abs(a_q - a_g).max() / np.abs(a_g).max()
    i_q = (X * (rc * w)[:, None]).sum(0)
    i_err = np.abs(i_q - I).max() / np.abs(I).max()
    _say(f"MAT21-20b (Orth) 20 modes, composite 12-pt GL (200 gel + 40 wall panels): "
         f"max offdiag/sqrt(diag diag) {off.max():.3e}; closed-form N vs GL diag {n_err:.3e}; "
         f"a_n Green vs GL projection {a_err:.3e}; I_n closed form vs GL {i_err:.3e}; "
         f"a_1..3 {a_g[:3]}")
    assert off.max() <= ORTH_BAND
    assert n_err <= ORTH_BAND
    assert a_err <= ORTH_BAND
    assert i_err <= ORTH_BAND


# ---------------------------------------------------------------- (Lim)
def test_late_time_equals_steady():
    lam = _lam()
    t_s = 1.0 / lam[0]
    rise = _rise_ss()
    zz = np.concatenate([Z_GEL, np.linspace(C_BOX, L_OUT, 9)[1:]])
    ss = two_layer_slab_steady_temperature(zz, half_thickness=C_BOX, wall=D_WALL, k=K,
                                           k_w=K_W, q=Q, h=H, t_inf=T_INF)
    err = np.abs(_theta(zz, 50 * t_s, lam=lam) - ss).max() / rise
    err_tc = np.abs(_theta(zz, 50 * T_C, lam=lam) - ss).max() / rise
    _say(f"MAT21-20b (Lim) theta(50 t_s = {50 * t_s:.1f} s) vs steady: {err:.3e} of the rise "
         f"{rise:.9f} K; (print) at 50 T_C = {50 * T_C:.1f} s: {err_tc:.3e}")
    assert err <= LIM_BAND


# ---------------------------------------------------------------- (IC)
def test_initial_condition_truncation_converges():
    z, w, _ = _gl_layers(1000, 100, order=8)
    ss = two_layer_slab_steady_temperature(z, half_thickness=C_BOX, wall=D_WALL, k=K,
                                           k_w=K_W, q=Q, h=H)
    nrm = np.sqrt((w * ss**2).sum())
    lam_all = _lam(400)
    errs = {}
    for n in (25, 100, 400):
        th0 = _theta(z, 0.0, lam=lam_all[:n]) - T_INF
        errs[n] = float(np.sqrt((w * th0**2).sum()) / nrm)
    ratios = [errs[4 * n] / errs[n] for n in (25, 100)]
    _say(f"MAT21-20b (IC) ||theta(0)||/||theta_ss|| (unweighted L2 on [0, C+d], composite 8-pt GL "
         f"1000 + 100 panels): " + ", ".join(f"N={n}: {e:.3e}" for n, e in errs.items())
         + f"; ratios 4N/N {', '.join(f'{r:.4f}' for r in ratios)}")
    assert ratios[0] <= 0.5 and ratios[1] <= 0.5


# ---------------------------------------------------------------- (Bal)
def test_energy_balance_at_t_star():
    lam = _lam(N_BAL)
    a = two_layer_slab_transient_coefficients(lam, q=Q, **WALL)
    _, I, _ = two_layer_slab_mode_integrals(lam, **WALL)
    rise_c = Q * C_BOX / wall_effective_heat_transfer(H, D_WALL, K_W)
    int_ss = (RHO_C * (rise_c * C_BOX + Q * C_BOX**3 / (3 * K))
              + RHO_C_W * (rise_c * D_WALL - (Q * C_BOX / K_W) * D_WALL**2 / 2))
    e = np.exp(-lam * T_STAR)
    lhs = int_ss - (a * I * e).sum()
    theta_L_ss = Q * C_BOX / H
    XL = two_layer_slab_mode(L_OUT, lam, **WALL)[0]
    time_int = theta_L_ss * T_STAR - (a * XL * (-np.expm1(-lam * T_STAR)) / lam).sum()
    rhs = Q * C_BOX * T_STAR - H * time_int
    err = abs(lhs - rhs) / (Q * C_BOX * T_STAR)
    _say(f"MAT21-20b (Bal) t* = 0.2 T_C = {T_STAR:.1f} s, {N_BAL} modes: int rho c theta dz "
         f"{lhs:.9e} J/m2 vs qCt* - h int theta(C+d) {rhs:.9e}; |diff|/(qCt*) {err:.3e}")
    assert err <= BAL_BAND


# ---------------------------------------------------------------- (K0) / (K∞)
def test_reductions_massless_and_conducting_wall():
    h_eff = wall_effective_heat_transfer(H, D_WALL, K_W)
    rise = _rise_ss()
    th_k0 = _theta(Z_GEL, T_STAR, rho_c_w=1e-9 * RHO_C_W)
    bare_eff = _bare_transient(Z_GEL, T_STAR, h_eff)
    k0 = np.abs(th_k0 - bare_eff).max() / rise
    rise_h = Q * C_BOX / H + Q * C_BOX**2 / (2 * K)
    th_kinf = _theta(Z_GEL, T_STAR, k_w=1e9, rho_c_w=1e-9 * RHO_C_W)
    bare_h = _bare_transient(Z_GEL, T_STAR, H)
    kinf = np.abs(th_kinf - bare_h).max() / rise_h
    th_cap = _theta(Z_GEL, T_STAR, k_w=1e9)
    cap = np.abs(th_cap - bare_h).max() / rise_h
    _say(f"MAT21-20b (K0) rho_w c_w = 1e-9 x nominal vs bare slab at h_eff {h_eff:.6f}, t*: "
         f"{k0:.3e} of the rise; (Kinf) k_w 1e9, rho_w c_w 1e-9 x vs bare at h: {kinf:.3e}; "
         f"(print) k_w 1e9 with rho_w c_w kept vs bare at h: {cap:.4e} "
         f"(centre {th_cap[20] - T_INF:.6f} K vs {bare_h[20] - T_INF:.6f} K)")
    assert k0 <= RED_BAND
    assert kinf <= RED_BAND


# ---------------------------------------------------------------- (ii) printed
def test_weekly_ii_transient_difference_printed():
    lam = _lam()
    t_s = 1.0 / lam[0]
    h_eff = wall_effective_heat_transfer(H, D_WALL, K_W)
    rise = _rise_ss()
    eps_w = RHO_C_W * D_WALL / (RHO_C * C_BOX)
    worst, t_w = 0.0, 0.0
    for t in np.linspace(0.0, 5 * t_s, 41):
        d = np.abs(_theta(Z_GEL, t, lam=lam) - _bare_transient(Z_GEL, t, h_eff)).max() / rise
        if d > worst:
            worst, t_w = d, t
    _say(f"MAT21-20b (ii, print) max_t max_gel |two-layer - bare(h_eff)|/rise_ss over 41 times "
         f"in [0, 5 t_s = {5 * t_s:.0f} s]: {worst:.6e} at t = {t_w:.1f} s ({t_w / T_C:.4f} T_C); "
         f"eps_w = {eps_w:.6f}; <= eps_w: {worst <= eps_w} (printed, not asserted)")
    assert np.isfinite(worst)   # print-only


# ---------------------------------------------------------------- negative control
def test_negative_control_bare_slab_at_h_misses_walled_centre():
    th = float(_theta(0.0, T_STAR)[0] - T_INF)
    bare = float(_bare_transient(0.0, T_STAR, H)[0] - T_INF)
    miss = abs(th - bare) / abs(th)
    _say(f"MAT21-20b (neg) centre rise at t* walled {th:.9f} K vs bare at h {bare:.9f} K; "
         f"relative miss {miss:.4e} (floor {NEG_FLOOR})")
    assert miss > NEG_FLOOR
