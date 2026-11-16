"""`MAT-21` step 0 — the rectangular gel box (F2182-sized), uniform source,
Robin walls on all six faces: the triple eigen-series in numpy, verified by
in-repo identities before any solve (§9 item 406).

|x| ≤ A, |y| ≤ B, |z| ≤ C with A, B, C = 0.325, 0.21, 0.045 m (650 × 420 ×
90 mm halved), `MAT-15`'s gel k 0.6, ρc 4.18e6, h 10, q 2000 (imported, not
restated).  Closed forms: `utils.analytical.robin_box_steady_temperature` /
`robin_box_transient_temperature`.

Anchors (asserted): (R0) roots in brackets, residual/γ ≤ 1e-12, two C forms
agree ≤ 1e-12 in max|ΔC|/max|C|; (R1) six-face Robin residuals/(qC) ≤ 1e-10,
steady and Fo ∈ FOS; (R2) z-reduced double series vs the triple series at
(40, 40, 200) ≤ 1e-7·θ_ss(0); (R3) truncated Gauss |h∮θ_N/(qV) − ΠS_d| ≤
1e-10; (R4) Duhamel ≤ 1e-10·θ_ss(0); (R4s) |ΠU_d − 1| ≤ 1e-9 at the centre;
(R5) t = 100 t_C = steady ≤ 1e-12, A = B = 20C → slab ≤ 1e-6; (R6) cube
permutation spread ≤ 1e-12.  Negative control (rule (e), predicted, only the
separation > 1e-3 asserted): Dirichlet z-roots.
"""

from __future__ import annotations

import itertools

import numpy as np

from fem_em_solver.utils.analytical import (
    robin_box_steady_temperature,
    robin_box_transient_temperature,
    robin_slab_eigenvalues,
)
from tests.validation.test_mat15_finite_cylinder_closed_form import H, K, Q, RHO_C, T_INF

# PROJECT_PLAN.md:4814 — the F2182 phantom, 650 × 420 × 90 mm, halved.
A_BOX, B_BOX, C_BOX = 0.325, 0.21, 0.045
HALF = (A_BOX, B_BOX, C_BOX)
T_C = C_BOX**2 * RHO_C / K
FOS = (0.05, 0.2, 1.0)
SERIES_TERMS = (40, 40, 40)
ALPHA = K / RHO_C
KW = dict(k=K, q=Q, h=H)


def _say(msg):
    print(msg, flush=True)


def _steady(x, y, z, **kw):
    kw = {"half_widths": HALF, "n_terms": SERIES_TERMS, **kw}
    return robin_box_steady_temperature(x, y, z, t_inf=T_INF, **KW, **kw)


def _trans(x, y, z, t, **kw):
    kw = {"half_widths": HALF, "n_terms": SERIES_TERMS, **kw}
    return robin_box_transient_temperature(x, y, z, t, t_inf=T_INF, rho_c=RHO_C, **KW, **kw)


def _theta_ss0():
    return float(_steady(0.0, 0.0, 0.0) - T_INF)


def _axis_data(L, n):
    gam = robin_slab_eigenvalues(H * L / K, n)
    c = 2.0 * np.sin(gam) / (gam + np.sin(gam) * np.cos(gam))
    return gam, c


# ---------------------------------------------------------------- (R0)
def test_r0_roots_and_coefficients():
    m = 400
    n = np.arange(1, m + 1)
    for L in HALF:
        bi = H * L / K
        gam = robin_slab_eigenvalues(bi, m)
        assert np.all(gam > (n - 1) * np.pi) and np.all(gam < (n - 0.5) * np.pi)
        res = np.abs(gam * np.sin(gam) - bi * np.cos(gam)) / gam
        c1 = 2.0 * np.sin(gam) / (gam + np.sin(gam) * np.cos(gam))
        c2 = 2.0 * bi * np.cos(gam) * (gam**2 + bi**2) / (gam**2 * (gam**2 + bi**2 + bi))
        dc = np.max(np.abs(c1 - c2)) / np.max(np.abs(c1))
        dc_rel = np.max(np.abs(c1 - c2) / np.abs(c1))
        _say(f"MAT21-0 (R0) Bi {bi:.4f}: gamma_1 {gam[0]:.6f}; max residual/gamma {res.max():.3e}; "
             f"C two forms max|dC|/max|C| {dc:.3e} (per-root rel {dc_rel:.3e}, printed only)")
        assert res.max() <= 1e-12
        assert dc <= 1e-12


# ---------------------------------------------------------------- (R1)
def _face_residuals(fn, half=HALF):
    """Per-face max |−k∂ₙθ − hθ|/(qC), faces ordered (±x, ±y, ±z)."""
    s = np.linspace(-0.95, 0.95, 7)
    u, v = np.meshgrid(s, s, indexing="ij")
    u, v = u.ravel(), v.ravel()
    out = []
    for d in range(3):
        o = [i for i in range(3) if i != d]
        for sgn in (1.0, -1.0):
            pts = [None, None, None]
            pts[d] = np.full_like(u, sgn * half[d])
            pts[o[0]] = u * half[o[0]]
            pts[o[1]] = v * half[o[1]]
            T, g = fn(*pts)
            res = np.abs(-K * sgn * g[d] - H * (T - T_INF)) / (Q * C_BOX)
            out.append(res.max())
    return out


def test_r1_robin_residuals():
    res = _face_residuals(lambda x, y, z: _steady(x, y, z, return_gradient=True))
    _say(f"MAT21-0 (R1) steady per-face residual/(qC) (+x,-x,+y,-y,+z,-z): "
         + ", ".join(f"{r:.3e}" for r in res))
    assert max(res) <= 1e-10
    for fo in FOS:
        t = fo * T_C
        res = _face_residuals(lambda x, y, z: _trans(x, y, z, t, return_gradient=True))
        _say(f"MAT21-0 (R1) transient Fo {fo}: max residual/(qC) {max(res):.3e}")
        assert max(res) <= 1e-10


# ---------------------------------------------------------------- (R2)
def _z_reduced(x, y, z, nl, nm):
    """θ = (q/k) Σ_lm C_l C_m X_l Y_m g_lm(z)/μ², exact in z; scaled cosh ratio."""
    gx, cx = _axis_data(A_BOX, nl)
    gy, cy = _axis_data(B_BOX, nm)
    mu = np.sqrt((gx[:, None] / A_BOX) ** 2 + (gy[None, :] / B_BOX) ** 2)
    az = abs(z)
    e2c = np.exp(-2.0 * mu * C_BOX)
    ratio = np.exp(mu * (az - C_BOX)) * (1.0 + np.exp(-2.0 * mu * az)) / (
        (K * mu / H) * (1.0 - e2c) + 1.0 + e2c)
    g = 1.0 - ratio
    w = cx[:, None] * cy[None, :] * g / mu**2
    X = np.cos(gx * x / A_BOX)
    Y = np.cos(gy * y / B_BOX)
    return (Q / K) * float(X @ w @ Y)


R2_PTS = ((0.0, 0.0, 0.0), (0.5 * A_BOX, 0.5 * B_BOX, 0.5 * C_BOX),
          (0.9 * A_BOX, 0.9 * B_BOX, 0.9 * C_BOX), (0.9 * A_BOX, 0.0, 0.0),
          (0.0, 0.0, 0.95 * C_BOX))


def test_r2_z_reduced_route():
    th0 = _theta_ss0()
    worst = 0.0
    for p in R2_PTS:
        a = float(_steady(*p, n_terms=(40, 40, 200)) - T_INF)
        b = _z_reduced(*p, 40, 40)
        worst = max(worst, abs(a - b) / th0)
        _say(f"MAT21-0 (R2) {p}: triple(40,40,200) {a:.12f}  reduced(40,40) {b:.12f}")
    _say(f"MAT21-0 (R2) max |triple - reduced|/theta_ss(0) = {worst:.3e} (predicted <~ 1e-9)")
    red0 = _z_reduced(0.0, 0.0, 0.0, 160, 160)
    p9 = R2_PTS[2]
    red9 = _z_reduced(*p9, 160, 160)
    d0 = abs(th0 - red0)
    d9 = abs(float(_steady(*p9) - T_INF) - red9) / th0
    _say(f"MAT21-0 printed theta_ss(0) default (40,40,40) {th0:.6f} vs reduced (160,160) "
         f"{red0:.6f} K: |diff| centre {d0:.3e} K (predicted 8.5e-7), "
         f"at 0.9(A,B,C) {d9:.3e} theta_ss(0) (predicted 2.4e-6)")
    assert worst <= 1e-7


# ---------------------------------------------------------------- (R3)
def test_r3_gauss_truncated_identity():
    xg, wg = np.polynomial.legendre.leggauss(128)
    total = 0.0
    area = 0.0
    for d in range(3):
        o = [i for i in range(3) if i != d]
        u, v = np.meshgrid(xg, xg, indexing="ij")
        wuv = np.outer(wg, wg).ravel() * HALF[o[0]] * HALF[o[1]]
        for sgn in (1.0, -1.0):
            pts = [None, None, None]
            pts[d] = np.full(u.size, sgn * HALF[d])
            pts[o[0]] = u.ravel() * HALF[o[0]]
            pts[o[1]] = v.ravel() * HALF[o[1]]
            total += np.sum(wuv * (_steady(*pts) - T_INF))
            area += np.sum(wuv)
    vol = 8.0 * A_BOX * B_BOX * C_BOX
    prod = 1.0
    defs = []
    for L, n in zip(HALF, SERIES_TERMS):
        gam, c = _axis_data(L, n)
        s = np.sum(c * np.sin(gam) / gam)
        prod *= s
        defs.append(1.0 - s)
    prod60 = 1.0
    defs60 = []
    for L in HALF:
        gam, c = _axis_data(L, 60)
        s = np.sum(c * np.sin(gam) / gam)
        prod60 *= s
        defs60.append(1.0 - s)
    ratio = H * total / (Q * vol)
    err = abs(ratio - prod)
    _say(f"MAT21-0 (R3) h*surface integral/(qV) {ratio:.15f}, prod S_d(40) {prod:.15f}, "
         f"|diff| {err:.3e} (predicted <= 1e-13)")
    _say(f"MAT21-0 printed truncation deficit 1 - prod S_d: 40 {1 - prod:.4e} (predicted 4.671e-6), "
         f"60 {1 - prod60:.4e} (predicted 1.368e-6); S-deficits at 60 "
         + " / ".join(f"{x:.3e}" for x in defs60) + " (predicted 9.52e-7 / 3.98e-7 / 1.83e-8)")
    _say(f"MAT21-0 printed surface area {area:.4f} m^2 (predicted 0.7386); face-mean theta "
         f"{total / area:.4f} K beside qV/(hA_s) {Q * vol / (H * area):.4f} K (predicted 6.6531)")
    assert err <= 1e-10


# ---------------------------------------------------------------- (R4)
def _duhamel(x, y, z, t, n=SERIES_TERMS, roots_z=None):
    edges = np.concatenate([[0.0], np.geomspace(1e-10 * t, t, 80)])
    xg, wg = np.polynomial.legendre.leggauss(16)
    s = (0.5 * (edges[1:, None] - edges[:-1, None]) * (xg + 1.0) + edges[:-1, None]).ravel()
    ws = (0.5 * (edges[1:, None] - edges[:-1, None]) * wg).ravel()
    return Q / RHO_C * np.sum(ws * _uprod(x, y, z, s, n))


def _uprod(x, y, z, s, n):
    s = np.atleast_1d(s)
    out = np.ones_like(s, dtype=float)
    for coord, L, nd in zip((x, y, z), HALF, n):
        gam, c = _axis_data(L, nd)
        out = out * (np.exp(-ALPHA * np.outer(s, (gam / L) ** 2)) * (c * np.cos(gam * coord / L))).sum(1)
    return out


R4_PTS = ((0.0, 0.0, 0.0), (0.5 * A_BOX, 0.5 * B_BOX, 0.5 * C_BOX),
          (0.9 * A_BOX, 0.8 * B_BOX, 0.8 * C_BOX))


def _r4_error(**kw):
    th0 = _theta_ss0()
    worst = 0.0
    tag = "-ctl" if kw else ""
    for fo in FOS:
        t = fo * T_C
        for p in R4_PTS:
            a = float(_trans(*p, t, **kw) - T_INF)
            b = _duhamel(*p, t)
            worst = max(worst, abs(a - b) / th0)
            _say(f"MAT21-0 (R4{tag}) Fo {fo} {p}: series {a:.12f}  Duhamel {b:.12f}  "
                 f"|d|/theta_ss(0) {abs(a - b) / th0:.3e}")
    return worst


def test_r4_duhamel():
    worst = _r4_error()
    _say(f"MAT21-0 (R4) max |series - Duhamel|/theta_ss(0) = {worst:.3e} (predicted <= 1e-13)")
    th_star = float(_trans(0.0, 0.0, 0.0, 0.2 * T_C) - T_INF)
    th_c = float(_trans(0.0, 0.0, 0.0, T_C) - T_INF)
    _say(f"MAT21-0 printed theta(0, 0.2 t_C) = {th_star:.6f} K (predicted 1.335705; q t/rho c "
         f"{Q * 0.2 * T_C / RHO_C:.4f}); theta(0, t_C) = {th_c:.6f} K (predicted 5.511718)")
    assert worst <= 1e-10


def test_r4s_initial_slope():
    s = 0.01 * T_C
    n = (200, 200, 60)
    pc = float(_uprod(0.0, 0.0, 0.0, s, n)[0])
    ph = float(_uprod(0.5 * A_BOX, 0.5 * B_BOX, 0.5 * C_BOX, s, n)[0])
    t = s
    dt = 1e-2 * t
    fd = (float(_trans(0.0, 0.0, 0.0, t + dt, n_terms=n)) - float(_trans(0.0, 0.0, 0.0, t - dt, n_terms=n))) \
        / (2 * dt) / (Q / RHO_C)
    _say(f"MAT21-0 (R4s) centre prod U_d at s = 0.01 t_C: |prod - 1| = {abs(pc - 1):.3e} (predicted 4.4e-14)")
    _say(f"MAT21-0 printed prod U_d at 0.5(A,B,C): 1 - {1 - ph:.3e} (predicted 1 - 1.05e-5); "
         f"FD slope / (q/rho c) at centre {fd:.12f}, |.-1| {abs(fd - 1):.3e} (predicted within 1e-8)")
    assert abs(pc - 1.0) <= 1e-9


# ---------------------------------------------------------------- (R5)
def test_r5_limits():
    pts = (np.array([0.0, 0.5 * A_BOX, 0.9 * A_BOX]), np.array([0.0, 0.3 * B_BOX, 0.7 * B_BOX]),
           np.array([0.0, 0.5 * C_BOX, 0.9 * C_BOX]))
    ss = _steady(*pts)
    late = _trans(*pts, 100.0 * T_C)
    e_late = np.max(np.abs(late - ss) / np.abs(ss - T_INF))
    zs = np.array([0.0, 0.5 * C_BOX, 0.9 * C_BOX])
    slab = Q * (C_BOX**2 - zs**2) / (2 * K) + Q * C_BOX / H
    got = robin_box_steady_temperature(
        np.zeros_like(zs), np.zeros_like(zs), zs, half_widths=(20 * C_BOX, 20 * C_BOX, C_BOX),
        n_terms=(400, 400, 60), **KW)
    e = np.abs(got / slab - 1.0)
    _say(f"MAT21-0 (R5a) t = 100 t_C vs steady: {e_late:.3e} (predicted 0)")
    _say(f"MAT21-0 (R5b) A = B = 20C, n = (400,400,60): |theta/theta_slab - 1| at z = 0/0.5C/0.9C "
         + " / ".join(f"{v:.3e}" for v in e) + " (predicted 3.98e-7 / - / 3.94e-7)")
    assert e_late <= 1e-12
    assert e.max() <= 1e-6


# ---------------------------------------------------------------- (R6)
def test_r6_cube_permutation():
    cube = (C_BOX, C_BOX, C_BOX)
    p = np.array([0.2, 0.5, 0.8]) * C_BOX
    perms = np.array(list(itertools.permutations(p)))
    for label, fn in (("steady", lambda x, y, z: robin_box_steady_temperature(
            x, y, z, half_widths=cube, **KW)),
                      ("Fo 0.2", lambda x, y, z: robin_box_transient_temperature(
            x, y, z, 0.2 * T_C, half_widths=cube, rho_c=RHO_C, **KW))):
        v = fn(perms[:, 0], perms[:, 1], perms[:, 2])
        c0 = float(fn(0.0, 0.0, 0.0))
        spread = (v.max() - v.min()) / c0
        _say(f"MAT21-0 (R6) cube {label}: permutation spread/theta(0) {spread:.3e} (predicted <= 1e-14)")
        assert spread <= 1e-12


# ---------------------------------------------------------------- printed records
def test_printed_records():
    (g1x, _), (g1y, _), (g1z, _) = (_axis_data(L, 1) for L in HALF)
    _say(f"MAT21-0 printed gamma_1 x/y/z {g1x[0]:.6f} / {g1y[0]:.6f} / {g1z[0]:.6f} "
         f"(predicted 1.330018 / 1.232272 / 0.771359)")
    vals = [float(_steady(0.0, 0.0, 0.0, n_terms=(n, n, n)) - T_INF) for n in (40, 60, 120)]
    slab = Q * C_BOX**2 / (2 * K) + Q * C_BOX / H
    _say(f"MAT21-0 printed theta_ss(0) n = 40/60/120: {vals[0]:.6f} / {vals[1]:.6f} / {vals[2]:.6f} K "
         f"(predicted 12.003947 / 12.003956 / -; reduced 12.003958; slab {slab:.4f})")
    corner = float(_steady(A_BOX, B_BOX, C_BOX) - T_INF)
    fx = float(_steady(A_BOX, 0.0, 0.0) - T_INF)
    fy = float(_steady(0.0, B_BOX, 0.0) - T_INF)
    fz = float(_steady(0.0, 0.0, C_BOX) - T_INF)
    _say(f"MAT21-0 printed corner {corner:.4f} K (predicted 2.5526); face centres x/y/z "
         f"{fx:.4f} / {fy:.4f} / {fz:.4f} K (predicted 6.1505 / 6.2367 / 8.7340)")
    lam = (g1x[0] / A_BOX) ** 2 + (g1y[0] / B_BOX) ** 2 + (g1z[0] / C_BOX) ** 2
    _say(f"MAT21-0 printed Lambda_111 {lam:.3f} m^-2 (predicted 345.005); tau_1 "
         f"{1 / (ALPHA * lam):.0f} s (predicted 20193); t_C {T_C:.1f} s")
    d = 1e-2 * C_BOX
    th = lambda x, y, z: float(_steady(x, y, z) - T_INF)  # noqa: E731
    c = th(0.0, 0.0, 0.0)
    lap = (th(d, 0, 0) + th(-d, 0, 0) + th(0, d, 0) + th(0, -d, 0) + th(0, 0, d) + th(0, 0, -d)
           - 6 * c) / d**2
    _say(f"MAT21-0 printed centre FD Laplacian residual |k lap + q|/q (delta 1e-2 C) "
         f"{abs(K * lap + Q) / Q:.3e} (predicted 1e-4 to 1e-3; never asserted)")


# ---------------------------------------------------------------- negative control
def test_negative_control_dirichlet_z():
    rz = (np.arange(1, 41) - 0.5) * np.pi
    roots = (None, None, rz)
    t = 0.2 * T_C
    res = _face_residuals(lambda x, y, z: _trans(x, y, z, t, roots=roots, return_gradient=True))
    zres = max(res[4:])
    worst = _r4_error(roots=roots)
    ss_c = float(_steady(0.0, 0.0, 0.0, roots=roots) - T_INF)
    _say(f"MAT21-0 control (Dirichlet z): (R1-ctl) z-face residual/(qC) {zres:.3e} "
         f"(predicted 0.95-1.0); lateral faces {max(res[:4]):.3e} (printed, round-off)")
    _say(f"MAT21-0 control (R4-ctl) max miss/theta_ss(0) {worst:.3e} (predicted 0.235)")
    _say(f"MAT21-0 control printed steady Dirichlet-z centre {ss_c:.4f} K vs qC^2/(2k) "
         f"{Q * C_BOX**2 / (2 * K):.4f} K (predicted 3.3735)")
    assert zres > 1e-3
    assert worst > 1e-3
