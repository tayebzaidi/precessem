"""`MAT-21` step 6, W1 — the Robin gel box under the separable source
q(x) = Q·[1 − cos(πx/A)] (§9 item 467, §10 chain P32): the triple Robin
series generalised by the cosine projection D⁽ʲ⁾, checked in numpy before any
solve.  Numpy only (real build, no dolfinx).

Bars are step 0's literals (`test_mat21_box_closed_form.py:10–17`) applied
to the same identities: (Q1) projection vs composite Gauss–Legendre ≤ 1e-12
(R0); (Q2) uniform profile bitwise = step 0's series; (Q3) six-face Robin
residuals/(QC) ≤ 1e-10 (R1); (Q4) truncated Gauss ≤ 1e-10 (R3); (Q5)
x-reduced route ≤ 1e-7·θ_ref (R2); (Q6) Duhamel ≤ 1e-10·θ_ref (R4), t = 100
t_C = steady ≤ 1e-12 (R5a).  Negative control (predicted, only the separation
> 1e-3 asserted): the flipped profile 1 + cos.
"""

from __future__ import annotations

import numpy as np

from fem_em_solver.utils.analytical import (
    robin_box_separable_steady_temperature,
    robin_box_separable_transient_temperature,
    robin_box_steady_temperature,
    robin_box_transient_temperature,
    robin_slab_cosine_coefficients,
    robin_slab_eigenvalues,
)
from tests.validation.test_mat15_finite_cylinder_closed_form import T_INF
from tests.validation.test_mat21_box_closed_form import (
    A_BOX,
    B_BOX,
    C_BOX,
    FOS,
    HALF,
    R2_PTS,
    R4_PTS,
    RHO_C,
    SERIES_TERMS,
    T_C,
    H,
    K,
    Q,
    _face_residuals,
)

PROFILES = ((1.0, -1.0, 1), (1.0, 0.0, 0), (1.0, 0.0, 0))
FLIPPED = ((1.0, 1.0, 1), (1.0, 0.0, 0), (1.0, 0.0, 0))
UNIFORM = ((1.0, 0.0, 0), (1.0, 0.0, 0), (1.0, 0.0, 0))
SEP_TERMS = (80, 40, 40)
ALPHA = K / RHO_C
KW = dict(k=K, q=Q, h=H)


def _say(msg):
    print(msg, flush=True)


def _sep(x, y, z, profiles=PROFILES, n_terms=SERIES_TERMS, **kw):
    return robin_box_separable_steady_temperature(
        x, y, z, half_widths=HALF, profiles=profiles, n_terms=n_terms, t_inf=T_INF, **KW, **kw
    )


def _sep_t(x, y, z, t, profiles=PROFILES, n_terms=SERIES_TERMS, **kw):
    return robin_box_separable_transient_temperature(
        x,
        y,
        z,
        t,
        half_widths=HALF,
        profiles=profiles,
        n_terms=n_terms,
        t_inf=T_INF,
        rho_c=RHO_C,
        **KW,
        **kw,
    )


def _axis(L, n, prof=(1.0, 0.0, 0)):
    gam = robin_slab_eigenvalues(H * L / K, n)
    c = 2.0 * np.sin(gam) / (gam + np.sin(gam) * np.cos(gam))
    a, b, j = prof
    p = c if (a == 1.0 and b == 0.0) else a * c + b * robin_slab_cosine_coefficients(gam, j)
    return gam, c, p


def _sprime(n_terms, profiles=PROFILES):
    out = 1.0
    for L, n, prof in zip(HALF, n_terms, profiles):
        gam, _, p = _axis(L, n, prof)
        out *= np.sum(p * np.sin(gam) / gam)
    return out


# ---------------------------------------------------------------- (Q1)
def test_q1_projection_vs_quadrature():
    xg, wg = np.polynomial.legendre.leggauss(16)
    worst = 0.0
    for L in HALF:
        bi = H * L / K
        gam = robin_slab_eigenvalues(bi, 80)
        edges = np.linspace(-L, L, 65)
        s = (0.5 * (edges[1:, None] - edges[:-1, None]) * (xg + 1.0) + edges[:-1, None]).ravel()
        ws = (0.5 * (edges[1:, None] - edges[:-1, None]) * wg).ravel()
        Xn = np.cos(np.outer(s, gam) / L)
        nrm = ws @ (Xn**2)
        for j in (1, 2):
            dq = (ws * np.cos(j * np.pi * s / L)) @ Xn / nrm
            d = robin_slab_cosine_coefficients(gam, j)
            err = np.max(np.abs(d - dq)) / np.max(np.abs(d))
            worst = max(worst, err)
            _say(f"MAT21-6 (Q1) Bi {bi:.4f} j {j}: max|dD|/max|D| {err:.3e} (predicted <= 1e-14)")
    _say(f"MAT21-6 (Q1) worst {worst:.3e}")
    assert worst <= 1e-12


# ---------------------------------------------------------------- (Q2)
def test_q2_uniform_bitwise():
    pts = np.array(R2_PTS)
    g = np.meshgrid(
        np.linspace(-0.9, 0.9, 7) * A_BOX,
        np.linspace(-0.9, 0.9, 7) * B_BOX,
        np.linspace(-0.9, 0.9, 3) * C_BOX,
        indexing="ij",
    )
    x = np.concatenate([pts[:, 0], g[0].ravel()])
    y = np.concatenate([pts[:, 1], g[1].ravel()])
    z = np.concatenate([pts[:, 2], g[2].ravel()])
    a = _sep(x, y, z, profiles=UNIFORM)
    b = robin_box_steady_temperature(
        x, y, z, half_widths=HALF, n_terms=SERIES_TERMS, t_inf=T_INF, **KW
    )
    t = 0.2 * T_C
    at = _sep_t(x, y, z, t, profiles=UNIFORM)
    bt = robin_box_transient_temperature(
        x, y, z, t, half_widths=HALF, n_terms=SERIES_TERMS, t_inf=T_INF, rho_c=RHO_C, **KW
    )
    _say(
        f"MAT21-6 (Q2) uniform profile vs step 0: steady equal {np.array_equal(a, b)}, "
        f"Fo 0.2 equal {np.array_equal(at, bt)} over {x.size} points"
    )
    assert np.array_equal(a, b)
    assert np.array_equal(at, bt)


# ---------------------------------------------------------------- (Q3)
def test_q3_robin_residuals():
    res = _face_residuals(lambda x, y, z: _sep(x, y, z, return_gradient=True))
    _say("MAT21-6 (Q3) steady per-face residual/(QC): " + ", ".join(f"{r:.3e}" for r in res))
    assert max(res) <= 1e-10
    for fo in FOS:
        t = fo * T_C
        res = _face_residuals(lambda x, y, z: _sep_t(x, y, z, t, return_gradient=True))
        _say(f"MAT21-6 (Q3) transient Fo {fo}: max residual/(QC) {max(res):.3e}")
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
            total += np.sum(wuv * (_sep(*pts) - T_INF))
    vol = 8.0 * A_BOX * B_BOX * C_BOX
    ratio = H * total / (Q * vol)
    prod = _sprime(SERIES_TERMS)
    err = abs(ratio - prod)
    _say(
        f"MAT21-6 (Q4) h*surface/(QV) {ratio:.15f}, prod S' (40,40,40) {prod:.15f}, "
        f"|diff| {err:.3e} (predicted <= 1e-14)"
    )
    _say(
        f"MAT21-6 printed floor 1 - prod S': (40,40,40) {1 - prod:.4e} (predicted 7.9221e-6); "
        f"SEP_TERMS {1 - _sprime(SEP_TERMS):.4e} (predicted 2.2196e-6); step 0 uniform "
        f"{1 - _sprime(SERIES_TERMS, UNIFORM):.4e} (4.6708e-6)"
    )
    assert err <= 1e-10


# ---------------------------------------------------------------- (Q5)
def _x_reduced(x, y, z, nm, nn, profile_sign=-1.0):
    """θ = Σ_mn C_m C_n Y_m Z_n g_mn(x), exact in x, for q = Q[1 + s cos(πx/A)]."""
    gy, cy, _ = _axis(B_BOX, nm)
    gz, cz, _ = _axis(C_BOX, nn)
    mu2 = (gy[:, None] / B_BOX) ** 2 + (gz[None, :] / C_BOX) ** 2
    mu = np.sqrt(mu2)
    p = np.pi / A_BOX
    s = profile_sign
    gpa = Q / (K * mu2) - s * Q / (K * (mu2 + p**2))  # g_p(A), cos(pA) = -1
    ax = abs(x)
    e2a = np.exp(-2.0 * mu * A_BOX)
    hom = (
        -H
        * gpa
        * np.exp(mu * (ax - A_BOX))
        * (1.0 + np.exp(-2.0 * mu * ax))
        / (K * mu * (1.0 - e2a) + H * (1.0 + e2a))
    )
    g = Q / (K * mu2) + s * Q * np.cos(p * x) / (K * (mu2 + p**2)) + hom
    Y = np.cos(gy * y / B_BOX)
    Z = np.cos(gz * z / C_BOX)
    return float(Y @ (cy[:, None] * cz[None, :] * g) @ Z)


def test_q5_x_reduced_route():
    red = [_x_reduced(*p, 40, 40) for p in R2_PTS]
    th_ref = max(abs(r) for r in red)
    worst = 0.0
    for p, b in zip(R2_PTS, red):
        a = float(_sep(*p, n_terms=(200, 40, 40)) - T_INF)
        worst = max(worst, abs(a - b) / th_ref)
        _say(f"MAT21-6 (Q5) {p}: triple(200,40,40) {a:.12f}  reduced(40,40) {b:.12f}")
    _say(
        f"MAT21-6 (Q5) theta_ref {th_ref:.6f} K (predicted 14.92); max |triple - reduced|/theta_ref "
        f"{worst:.3e} (predicted <= 3e-9)"
    )
    assert worst <= 1e-7


# ---------------------------------------------------------------- (Q6)
def _uprod_sep(x, y, z, s, n=SERIES_TERMS):
    s = np.atleast_1d(s)
    out = np.ones_like(s, dtype=float)
    for coord, L, nd, prof in zip((x, y, z), HALF, n, PROFILES):
        gam, _, p = _axis(L, nd, prof)
        out = out * (
            np.exp(-ALPHA * np.outer(s, (gam / L) ** 2)) * (p * np.cos(gam * coord / L))
        ).sum(1)
    return out


def _duhamel_sep(x, y, z, t):
    edges = np.concatenate([[0.0], np.geomspace(1e-10 * t, t, 80)])
    xg, wg = np.polynomial.legendre.leggauss(16)
    s = (0.5 * (edges[1:, None] - edges[:-1, None]) * (xg + 1.0) + edges[:-1, None]).ravel()
    ws = (0.5 * (edges[1:, None] - edges[:-1, None]) * wg).ravel()
    return Q / RHO_C * np.sum(ws * _uprod_sep(x, y, z, s))


def test_q6_duhamel_and_late_time():
    th_ref = max(abs(_x_reduced(*p, 40, 40)) for p in R2_PTS)
    worst = 0.0
    for fo in FOS:
        t = fo * T_C
        for p in R4_PTS:
            a = float(_sep_t(*p, t) - T_INF)
            b = _duhamel_sep(*p, t)
            worst = max(worst, abs(a - b) / th_ref)
            _say(f"MAT21-6 (Q6) Fo {fo} {p}: series {a:.12f}  Duhamel {b:.12f}")
    _say(f"MAT21-6 (Q6) max |series - Duhamel|/theta_ref {worst:.3e} (predicted <= 1e-13)")
    pts = (
        np.array([0.0, 0.5 * A_BOX, 0.9 * A_BOX]),
        np.array([0.0, 0.3 * B_BOX, 0.7 * B_BOX]),
        np.array([0.0, 0.5 * C_BOX, 0.9 * C_BOX]),
    )
    ss = _sep(*pts)
    late = _sep_t(*pts, 100.0 * T_C)
    e_late = np.max(np.abs(late - ss) / np.abs(ss - T_INF))
    _say(f"MAT21-6 (Q6) t = 100 t_C vs steady: {e_late:.3e} (predicted 0)")
    assert worst <= 1e-10
    assert e_late <= 1e-12


# ---------------------------------------------------------------- printed records
def test_printed_records():
    for L, lab in zip(HALF, "xyz"):
        gam, c, _ = _axis(L, 3)
        d = robin_slab_cosine_coefficients(gam, 1)
        _say(
            f"MAT21-6 printed {lab}: gamma {np.array2string(gam, precision=6)}; C "
            f"{np.array2string(c, precision=6)}; D1 {np.array2string(d, precision=6)}; "
            f"C - D1 {np.array2string(c - d, precision=6)}"
        )
    _say(
        "MAT21-6 printed drafting (x): gamma 1.330018/4.068213/6.945535; C 1.243794/-0.351594/0.165527; "
        "D 0.271608/0.871010/-0.208103; P 0.972185/-1.222604/0.373630"
    )
    pts = ((0.0, 0.0, 0.0), (0.9 * A_BOX, 0.0, 0.0), (0.9 * A_BOX, 0.9 * B_BOX, 0.9 * C_BOX))
    vs = [float(_sep(*p, n_terms=SEP_TERMS) - T_INF) for p in pts]
    vu = [float(_sep(*p, profiles=UNIFORM) - T_INF) for p in pts]
    _say(
        f"MAT21-6 printed steady sep (80,40,40) {vs[0]:.6f} / {vs[1]:.6f} / {vs[2]:.6f} K "
        f"(predicted 2.741147 / 14.919189 / 8.541895); uniform {vu[0]:.6f} / {vu[1]:.6f} / {vu[2]:.6f} "
        f"(predicted 12.003947 / 8.706255 / 4.852187)"
    )
    ts = 0.2 * T_C
    t0 = float(_sep_t(0.0, 0.0, 0.0, ts, n_terms=SEP_TERMS) - T_INF)
    t9 = float(_sep_t(0.9 * A_BOX, 0.0, 0.0, ts, n_terms=SEP_TERMS) - T_INF)
    _say(f"MAT21-6 printed transient t* sep {t0:.6f} / {t9:.6f} K (predicted 0.024797 / 2.537123)")


# ---------------------------------------------------------------- negative control
def test_q_ctl_flipped_profile_misses():
    red = [_x_reduced(*p, 40, 40) for p in R2_PTS]
    th_ref = max(abs(r) for r in red)
    worst = 0.0
    for p, b in zip(R2_PTS, red):
        a = float(_sep(*p, profiles=FLIPPED, n_terms=(200, 40, 40)) - T_INF)
        worst = max(worst, abs(a - b) / th_ref)
    c = float(_sep(0.0, 0.0, 0.0, profiles=FLIPPED, n_terms=(200, 40, 40)) - T_INF)
    _say(
        f"MAT21-6 control (1 + cos vs 1 - cos reduced): max miss/theta_ref {worst:.3e} "
        f"(predicted 1.24); centre {c:.4f} vs {red[0]:.4f} K (predicted 21.27 vs 2.74)"
    )
    assert worst > 1e-3
