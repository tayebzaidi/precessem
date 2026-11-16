"""`MAT-21` step 9, W1 — the Robin gel box under the transverse quadratic
source q = 3Q(y² + z²)/(B² + C²) = Q[w_y (y/B)² + w_z (z/C)²] (§9 item 497,
§10 chain P36; brief `docs/planning/chunks/MAT-21.md`, review ruling
2026-10-16 10:30).  Numpy only (real build, no dolfinx).

Bars are step 0's literals (`test_mat21_box_closed_form.py:10–17`) applied to
the same identities as P32's W1: (P0) w_y + w_z = 3; (Q1) projection E_n vs
composite Gauss–Legendre ≤ 1e-12 (R0); (Q1b) Robin-root identity ≤ 1e-12;
(Q2) unit weights bitwise = step 0's series; (Q3) six-face Robin residuals ≤
1e-10·QC (R1); (Q4) truncated Gauss ≤ 1e-10 (R3); (Q5) y- and z-reduced
routes ≤ 1e-7·θ_ref (R2); (Q6) Duhamel ≤ 1e-10·θ_ref (R4), t = 100 t_C =
steady ≤ 1e-12 (R5a).  Negative control (predicted; only the separation >
1e-3 asserted): the z² term dropped.
"""

from __future__ import annotations

import numpy as np

from fem_em_solver.utils.analytical import (
    robin_box_quadratic_steady_temperature,
    robin_box_quadratic_transient_temperature,
    robin_box_steady_temperature,
    robin_box_transient_temperature,
    robin_slab_eigenvalues,
    robin_slab_quadratic_coefficients,
)
from tests.validation.test_mat15_finite_cylinder_closed_form import T_INF
from tests.validation.test_mat21_box_closed_form import (
    A_BOX,
    B_BOX,
    C_BOX,
    FOS,
    H,
    HALF,
    K,
    Q,
    R2_PTS,
    R4_PTS,
    RHO_C,
    SERIES_TERMS,
    T_C,
    _face_residuals,
)

WEIGHTS = (0.0, 0.0, 3.0 * B_BOX**2 / (B_BOX**2 + C_BOX**2), 3.0 * C_BOX**2 / (B_BOX**2 + C_BOX**2))
WEIGHTS_NOZ = (0.0, 0.0, WEIGHTS[2], 0.0)
UNIT = (1.0, 0.0, 0.0, 0.0)
QUAD_TERMS = (80, 80, 40)
ALPHA = K / RHO_C
KW = dict(k=K, q=Q, h=H)


def _say(msg):
    print(msg, flush=True)


def _quad(x, y, z, weights=WEIGHTS, n_terms=SERIES_TERMS, **kw):
    return robin_box_quadratic_steady_temperature(
        x, y, z, half_widths=HALF, weights=weights, n_terms=n_terms, t_inf=T_INF, **KW, **kw)


def _quad_t(x, y, z, t, weights=WEIGHTS, n_terms=SERIES_TERMS, **kw):
    return robin_box_quadratic_transient_temperature(
        x, y, z, t, half_widths=HALF, weights=weights, n_terms=n_terms, t_inf=T_INF,
        rho_c=RHO_C, **KW, **kw)


def _axis(L, n):
    gam = robin_slab_eigenvalues(H * L / K, n)
    c = 2.0 * np.sin(gam) / (gam + np.sin(gam) * np.cos(gam))
    return gam, c, robin_slab_quadratic_coefficients(gam)


def _gauss_target(n_terms, weights=WEIGHTS):
    S, S2 = [], []
    for L, n in zip(HALF, n_terms):
        gam, c, e = _axis(L, n)
        S.append(np.sum(c * np.sin(gam) / gam))
        S2.append(np.sum(e * np.sin(gam) / gam))
    w0, wx, wy, wz = weights
    return (w0 * S[0] * S[1] * S[2] + wx * S2[0] * S[1] * S[2] + wy * S[0] * S2[1] * S[2]
            + wz * S[0] * S[1] * S2[2]), S2


# ---------------------------------------------------------------- (P0)
def test_p0_weights_sum_to_three():
    d = abs(WEIGHTS[2] + WEIGHTS[3] - 3.0)
    _say(f"MAT21-9 (P0) w_y {WEIGHTS[2]:.6f} (drafting 2.868293), w_z {WEIGHTS[3]:.6f} "
         f"(drafting 0.131707), |w_y + w_z - 3| {d:.3e}")
    assert d <= 1e-15


# ---------------------------------------------------------------- (Q1), (Q1b)
def test_q1_projection_vs_quadrature_and_root_identity():
    xg, wg = np.polynomial.legendre.leggauss(16)
    worst = 0.0
    worst_b = 0.0
    for L in HALF:
        bi = H * L / K
        gam = robin_slab_eigenvalues(bi, 80)
        edges = np.linspace(-L, L, 65)
        s = (0.5 * (edges[1:, None] - edges[:-1, None]) * (xg + 1.0) + edges[:-1, None]).ravel()
        ws = (0.5 * (edges[1:, None] - edges[:-1, None]) * wg).ravel()
        Xn = np.cos(np.outer(s, gam) / L)
        nrm = ws @ (Xn**2)
        eq = (ws * (s / L) ** 2) @ Xn / nrm
        e = robin_slab_quadratic_coefficients(gam)
        err = np.max(np.abs(e - eq)) / np.max(np.abs(e))
        c = 2.0 * np.sin(gam) / (gam + np.sin(gam) * np.cos(gam))
        eb = c * (1.0 + 2.0 / bi - 2.0 / gam**2)
        err_b = np.max(np.abs(e - eb)) / np.max(np.abs(e))
        worst = max(worst, err)
        worst_b = max(worst_b, err_b)
        _say(f"MAT21-9 (Q1) Bi {bi:.4f}: max|dE|/max|E| {err:.3e} (drafting <= 4.2e-15); "
             f"(Q1b) root identity {err_b:.3e} (predicted <= 1e-14)")
    _say(f"MAT21-9 (Q1) worst {worst:.3e}; (Q1b) worst {worst_b:.3e}")
    assert worst <= 1e-12
    assert worst_b <= 1e-12


# ---------------------------------------------------------------- (Q2)
def test_q2_unit_weights_bitwise():
    pts = np.array(R2_PTS)
    g = np.meshgrid(np.linspace(-0.9, 0.9, 7) * A_BOX, np.linspace(-0.9, 0.9, 7) * B_BOX,
                    np.linspace(-0.9, 0.9, 3) * C_BOX, indexing="ij")
    x = np.concatenate([pts[:, 0], g[0].ravel()])
    y = np.concatenate([pts[:, 1], g[1].ravel()])
    z = np.concatenate([pts[:, 2], g[2].ravel()])
    a = _quad(x, y, z, weights=UNIT)
    b = robin_box_steady_temperature(x, y, z, half_widths=HALF, n_terms=SERIES_TERMS, t_inf=T_INF, **KW)
    t = 0.2 * T_C
    at = _quad_t(x, y, z, t, weights=UNIT)
    bt = robin_box_transient_temperature(x, y, z, t, half_widths=HALF, n_terms=SERIES_TERMS,
                                         t_inf=T_INF, rho_c=RHO_C, **KW)
    _say(f"MAT21-9 (Q2) unit weights vs step 0: steady equal {np.array_equal(a, b)}, "
         f"Fo 0.2 equal {np.array_equal(at, bt)} over {x.size} points")
    assert np.array_equal(a, b)
    assert np.array_equal(at, bt)


# ---------------------------------------------------------------- (Q3)
def test_q3_robin_residuals():
    res = _face_residuals(lambda x, y, z: _quad(x, y, z, return_gradient=True))
    _say("MAT21-9 (Q3) steady per-face residual/(QC): " + ", ".join(f"{r:.3e}" for r in res))
    assert max(res) <= 1e-10
    for fo in FOS:
        t = fo * T_C
        res = _face_residuals(lambda x, y, z: _quad_t(x, y, z, t, return_gradient=True))
        _say(f"MAT21-9 (Q3) transient Fo {fo}: max residual/(QC) {max(res):.3e}")
        assert max(res) <= 1e-10


# ---------------------------------------------------------------- (Q4)
def test_q4_gauss_truncated_identity():
    xg, wg = np.polynomial.legendre.leggauss(128)
    total = 0.0
    for d in range(3):
        o = [i for i in range(3) if i != d]
        u, v = np.meshgrid(xg, xg, indexing="ij")
        wuv = np.outer(wg, wg).ravel() * HALF[o[0]] * HALF[o[1]]
        for sgn in (1.0, -1.0):
            pts = [None, None, None]
            pts[d] = np.full(u.size, sgn * HALF[d])
            pts[o[0]] = u.ravel() * HALF[o[0]]
            pts[o[1]] = v.ravel() * HALF[o[1]]
            total += np.sum(wuv * (_quad(*pts) - T_INF))
    vol = 8.0 * A_BOX * B_BOX * C_BOX
    ratio = H * total / (Q * vol)
    tgt, s2 = _gauss_target(SERIES_TERMS)
    err = abs(ratio - tgt)
    _say(f"MAT21-9 (Q4) h*surface/(QV) {ratio:.15f}, target (40,40,40) {tgt:.15f}, "
         f"|diff| {err:.3e} (predicted <= 1e-14)")
    _say("MAT21-9 S'' per axis (40 terms): " + " / ".join(f"{v:.10f}" for v in s2)
         + " (drafting 0.3333288835 / 0.3333311989 / 0.3333331044)")
    _say(f"MAT21-9 printed floor 1 - target: (40,40,40) {1 - tgt:.4e} (drafting 9.5217e-6); "
         f"QUAD_TERMS {1 - _gauss_target(QUAD_TERMS)[0]:.4e} (drafting 1.2484e-6); step 0 uniform "
         f"{1 - _gauss_target(SERIES_TERMS, UNIT)[0]:.4e} (4.6708e-6)")
    assert err <= 1e-10


# ---------------------------------------------------------------- (Q5)
def _hom(gp_b, gp_db, mu, s_abs, L):
    e2 = np.exp(-2.0 * mu * L)
    return -(K * gp_db + H * gp_b) * np.exp(mu * (s_abs - L)) * (1.0 + np.exp(-2.0 * mu * s_abs)) / (
        K * mu * (1.0 - e2) + H * (1.0 + e2))


def _y_reduced(x, y, z, nl, nn, weights=WEIGHTS):
    """θ = Σ_ln X_l Z_n g_ln(y), exact in y; E used only on z."""
    _, _, wy, wz = weights
    gx, cx, _ = _axis(A_BOX, nl)
    gz, cz, ez = _axis(C_BOX, nn)
    mu2 = (gx[:, None] / A_BOX) ** 2 + (gz[None, :] / C_BOX) ** 2
    mu = np.sqrt(mu2)
    a = Q * cx[:, None] * wy * cz[None, :] / K
    b = Q * cx[:, None] * wz * ez[None, :] / K
    al = a / (B_BOX**2 * mu2)
    be = (2.0 * al + b) / mu2
    g = al * y * y + be + _hom(al * B_BOX**2 + be, 2.0 * al * B_BOX, mu, abs(y), B_BOX)
    X = np.cos(gx * x / A_BOX)
    Z = np.cos(gz * z / C_BOX)
    return float(X @ g @ Z)


def _z_reduced(x, y, z, nl, nm, weights=WEIGHTS):
    """θ = Σ_lm X_l Y_m g_lm(z), exact in z; E used only on y."""
    _, _, wy, wz = weights
    gx, cx, _ = _axis(A_BOX, nl)
    gy, cy, ey = _axis(B_BOX, nm)
    mu2 = (gx[:, None] / A_BOX) ** 2 + (gy[None, :] / B_BOX) ** 2
    mu = np.sqrt(mu2)
    a = Q * cx[:, None] * wz * cy[None, :] / K
    b = Q * cx[:, None] * wy * ey[None, :] / K
    al = a / (C_BOX**2 * mu2)
    be = (2.0 * al + b) / mu2
    g = al * z * z + be + _hom(al * C_BOX**2 + be, 2.0 * al * C_BOX, mu, abs(z), C_BOX)
    X = np.cos(gx * x / A_BOX)
    Y = np.cos(gy * y / B_BOX)
    return float(X @ g @ Y)


def test_q5_reduced_routes():
    red_y = [_y_reduced(*p, 40, 40) for p in R2_PTS]
    red_z = [_z_reduced(*p, 40, 40) for p in R2_PTS]
    th_ref = max(abs(r) for r in red_y)
    wy_ = 0.0
    wz_ = 0.0
    for p, by, bz in zip(R2_PTS, red_y, red_z):
        ay = float(_quad(*p, n_terms=(40, 200, 40)) - T_INF)
        az = float(_quad(*p, n_terms=(40, 40, 200)) - T_INF)
        wy_ = max(wy_, abs(ay - by) / th_ref)
        wz_ = max(wz_, abs(az - bz) / th_ref)
        _say(f"MAT21-9 (Q5) {p}: triple(40,200,40) {ay:.12f} y-red {by:.12f}; "
             f"triple(40,40,200) {az:.12f} z-red {bz:.12f}")
    _say(f"MAT21-9 (Q5) theta_ref {th_ref:.6f} K (drafting 8.805072); y-route {wy_:.3e} "
         f"(drafting 2.48e-9), z-route {wz_:.3e} (drafting 1.46e-11)")
    assert wy_ <= 1e-7
    assert wz_ <= 1e-7


# ---------------------------------------------------------------- (Q6)
def _uprod_quad(x, y, z, s, n=SERIES_TERMS):
    s = np.atleast_1d(s)
    U, U2 = [], []
    for coord, L, nd in zip((x, y, z), HALF, n):
        gam, c, e = _axis(L, nd)
        ex = np.exp(-ALPHA * np.outer(s, (gam / L) ** 2))
        cs = np.cos(gam * coord / L)
        U.append((ex * (c * cs)).sum(1))
        U2.append((ex * (e * cs)).sum(1))
    _, _, wy, wz = WEIGHTS
    return wy * U[0] * U2[1] * U[2] + wz * U[0] * U[1] * U2[2]


def _duhamel_quad(x, y, z, t):
    edges = np.concatenate([[0.0], np.geomspace(1e-10 * t, t, 80)])
    xg, wg = np.polynomial.legendre.leggauss(16)
    s = (0.5 * (edges[1:, None] - edges[:-1, None]) * (xg + 1.0) + edges[:-1, None]).ravel()
    ws = (0.5 * (edges[1:, None] - edges[:-1, None]) * wg).ravel()
    return Q / RHO_C * np.sum(ws * _uprod_quad(x, y, z, s))


def test_q6_duhamel_and_late_time():
    th_ref = max(abs(_y_reduced(*p, 40, 40)) for p in R2_PTS)
    worst = 0.0
    for fo in FOS:
        t = fo * T_C
        for p in R4_PTS:
            a = float(_quad_t(*p, t) - T_INF)
            b = _duhamel_quad(*p, t)
            worst = max(worst, abs(a - b) / th_ref)
            _say(f"MAT21-9 (Q6) Fo {fo} {p}: series {a:.12f}  Duhamel {b:.12f}")
    _say(f"MAT21-9 (Q6) max |series - Duhamel|/theta_ref {worst:.3e} (predicted <= 1e-13)")
    pts = (np.array([0.0, 0.5 * A_BOX, 0.9 * A_BOX]), np.array([0.0, 0.3 * B_BOX, 0.7 * B_BOX]),
           np.array([0.0, 0.5 * C_BOX, 0.9 * C_BOX]))
    ss = _quad(*pts)
    late = _quad_t(*pts, 100.0 * T_C)
    e_late = np.max(np.abs(late - ss) / np.abs(ss - T_INF))
    _say(f"MAT21-9 (Q6) t = 100 t_C vs steady: {e_late:.3e} (predicted 0)")
    assert worst <= 1e-10
    assert e_late <= 1e-12


# ---------------------------------------------------------------- printed records
def test_printed_records():
    for L, lab in zip(HALF, "xyz"):
        gam, c, e = _axis(L, 3)
        _say(f"MAT21-9 printed {lab}: gamma {np.array2string(gam, precision=6)}; C "
             f"{np.array2string(c, precision=6)}; E {np.array2string(e, precision=6)}")
    _say("MAT21-9 printed drafting E: x 0.296788/-0.438925/0.219782; y 0.310445/-0.440810/0.195929; "
         "z 0.334851/-0.425321/0.129260")
    pts = ((0.0, 0.0, 0.0), (0.0, 0.9 * B_BOX, 0.0), (0.9 * A_BOX, 0.0, 0.0),
           (0.9 * A_BOX, 0.9 * B_BOX, 0.9 * C_BOX))
    vq = [float(_quad(*p, n_terms=QUAD_TERMS) - T_INF) for p in pts]
    vu = [float(_quad(*p, weights=UNIT) - T_INF) for p in pts[:3]]
    _say(f"MAT21-9 printed steady quad (80,80,40) {vq[0]:.6f} / {vq[1]:.6f} / {vq[2]:.6f} / "
         f"{vq[3]:.6f} K (drafting 4.276700 / 13.337014 / 2.682231 / 8.585064); uniform "
         f"{vu[0]:.6f} / {vu[1]:.6f} / {vu[2]:.6f} (drafting 12.003947 / 8.074905 / 8.706255)")
    ts = 0.2 * T_C
    vt = [float(_quad_t(*p, ts, n_terms=QUAD_TERMS) - T_INF) for p in pts[:3]]
    _say(f"MAT21-9 printed transient t* quad (80,80,40) {vt[0]:.6f} / {vt[1]:.6f} / {vt[2]:.6f} K "
         f"(drafting 0.063466 / 2.892661 / 0.061875)")
    ref = float(_y_reduced(0.9 * A_BOX, 0.9 * B_BOX, 0.9 * C_BOX, 200, 200))
    _say(f"MAT21-9 printed QUAD_TERMS truncation at 0.9(A,B,C): |triple - y-red(200,200)| "
         f"{abs(vq[3] - ref):.3e} K (drafting <= 3.8e-6)")


# ---------------------------------------------------------------- negative control
def test_q_ctl_z_term_dropped_misses():
    red = [_y_reduced(*p, 40, 40) for p in R2_PTS]
    th_ref = max(abs(r) for r in red)
    worst = 0.0
    for p, b in zip(R2_PTS, red):
        a = float(_quad(*p, weights=WEIGHTS_NOZ, n_terms=(40, 200, 40)) - T_INF)
        worst = max(worst, abs(a - b) / th_ref)
    _say(f"MAT21-9 control (Q-ctl) PREDICTED: WEIGHTS_NOZ vs full y-reduced: max miss/theta_ref "
         f"{worst:.4e} (drafting 0.0516; asserted only > 1e-3)")
    assert worst > 1e-3
