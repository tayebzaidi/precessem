"""`MAT-21` step 19 — the source-free Robin box from a non-uniform initial
temperature θ₀ = T∞ + g·z (a pre-equilibration top-to-bottom gradient), in
closed form, verified by in-repo identities (§9 item 659, §10 P47a).

Box and gel imported from step 0 (`test_mat21_box_closed_form.py`) and
`MAT-15` (k, ρc, h, T∞).  g = 10 K/m (0.9 K over the 90 mm depth) is the
fixture's choice, not a measured phantom value.  Closed form:
`utils.analytical.robin_box_initial_gradient_temperature` (odd z family
`robin_slab_odd_eigenvalues`, projection `robin_slab_odd_linear_coefficients`).

Asserted: (Odd) odd roots in ((n − ½)π, nπ), |γ cos γ + Bi sin γ|/γ ≤ 1e-12;
(Rob) six-face Robin residuals/(k g) ≤ 1e-10 at Fo ∈ FOS; (IC) the L2 error of
the z-series reconstruction of z at 160 terms ≤ ½ × that at 40; (Dec)
|θ(100 t_C)|/(g C) ≤ 1e-12; (Bal) whole-box energy balance ≤ 1e-10 (both sides
vanish by symmetry); (Half) the upper-half-box balance with the flux through
z = 0, ≤ 1e-10 of ρc g C²/2·4AB; (U) the uniform ḡ path equals
`robin_box_transient_temperature(q=0, t_init=T∞+ḡ)` ≤ 1e-12 K.  Control (Dz):
Dirichlet odd z roots (γ = mπ) miss the Robin solution at t* by > 1e-3 g C.
"""

from __future__ import annotations

import numpy as np

from fem_em_solver.utils.analytical import (
    robin_box_initial_gradient_temperature,
    robin_box_transient_temperature,
    robin_slab_eigenvalues,
    robin_slab_odd_eigenvalues,
    robin_slab_odd_linear_coefficients,
)
from tests.validation.test_mat15_finite_cylinder_closed_form import H, K, RHO_C, T_INF
from tests.validation.test_mat21_box_closed_form import (
    A_BOX,
    B_BOX,
    C_BOX,
    FOS,
    HALF,
    T_C,
)

G_GRAD = 10.0          # K/m — fixture's choice (0.9 K over 90 mm), not a measured value
ALPHA = K / RHO_C
N_TERMS = 40
T_STAR = 0.2 * T_C
GC = G_GRAD * C_BOX


def _say(msg):
    print(msg, flush=True)


def _theta(x, y, z, t, **kw):
    kw = {"n_terms": N_TERMS, **kw}
    return robin_box_initial_gradient_temperature(
        x, y, z, t, half_widths=HALF, k=K, h=H, rho_c=RHO_C, g=G_GRAD, t_inf=T_INF, **kw)


def _even(L, n):
    gam = robin_slab_eigenvalues(H * L / K, n)
    c = 2.0 * np.sin(gam) / (gam + np.sin(gam) * np.cos(gam))
    return gam, c


# ---------------------------------------------------------------- (Odd)
def test_odd_roots():
    m = 400
    n = np.arange(1, m + 1)
    for L in HALF:
        bi = H * L / K
        gam = robin_slab_odd_eigenvalues(bi, m)
        inside = np.all(gam > (n - 0.5) * np.pi) and np.all(gam < n * np.pi)
        res = np.abs(gam * np.cos(gam) + bi * np.sin(gam)) / gam
        _say(f"MAT21-19 (Odd) Bi {bi:.4f}: gamma'_1 {gam[0]:.12f}, gamma'_400 {gam[-1]:.6f}; "
             f"one root per ((n-1/2)pi, n pi): {inside}; max residual/gamma {res.max():.3e}")
        assert inside
        assert res.max() <= 1e-12


# ---------------------------------------------------------------- (Rob)
def _face_residuals(fn):
    """Per-face max |−k∂ₙθ − h(T − T∞)|/(k g), faces ordered (±x, ±y, ±z)."""
    s = np.linspace(-0.95, 0.95, 7)
    u, v = np.meshgrid(s, s, indexing="ij")
    u, v = u.ravel(), v.ravel()
    out = []
    for d in range(3):
        o = [i for i in range(3) if i != d]
        for sgn in (1.0, -1.0):
            pts = [None, None, None]
            pts[d] = np.full_like(u, sgn * HALF[d])
            pts[o[0]] = u * HALF[o[0]]
            pts[o[1]] = v * HALF[o[1]]
            T, gr = fn(*pts)
            res = np.abs(-K * sgn * gr[d] - H * (T - T_INF)) / (K * G_GRAD)
            out.append(res.max())
    return out


def test_robin_residuals():
    worst = 0.0
    for fo in FOS:
        res = _face_residuals(lambda x, y, z: _theta(x, y, z, fo * T_C, return_gradient=True))
        _say(f"MAT21-19 (Rob) Fo {fo} per-face residual/(k g) (+x,-x,+y,-y,+z,-z): "
             + ", ".join(f"{r:.3e}" for r in res))
        worst = max(worst, max(res))
    _say(f"MAT21-19 (Rob) max {worst:.3e} (bar 1e-10)")
    assert worst <= 1e-10


# ---------------------------------------------------------------- (IC)
def _ic_l2_error(n):
    gam = robin_slab_odd_eigenvalues(H * C_BOX / K, n)
    b = robin_slab_odd_linear_coefficients(gam)
    xg, wg = np.polynomial.legendre.leggauss(32)
    edges = np.linspace(-C_BOX, C_BOX, 401)
    zq = (0.5 * (edges[1:, None] - edges[:-1, None]) * (xg + 1.0) + edges[:-1, None]).ravel()
    wq = (0.5 * (edges[1:, None] - edges[:-1, None]) * wg).ravel()
    rec = C_BOX * (np.sin(np.outer(zq, gam / C_BOX)) @ b)
    err = np.sqrt(np.sum(wq * (zq - rec) ** 2))
    ref = np.sqrt(np.sum(wq * zq**2))
    return err / ref


def test_ic_reconstruction():
    ns = (10, 20, 40, 80, 160, 320)
    errs = {n: _ic_l2_error(n) for n in ns}
    _say("MAT21-19 (IC) rel L2 |z - S_N z| / |z| on (-C, C): "
         + ", ".join(f"N={n} {errs[n]:.4e}" for n in ns))
    rates = [np.log(errs[a] / errs[b]) / np.log(b / a) for a, b in zip(ns[:-1], ns[1:])]
    _say("MAT21-19 (IC) printed successive decay exponents: " + ", ".join(f"{r:.3f}" for r in rates))
    gam = robin_slab_odd_eigenvalues(H * C_BOX / K, 4)
    zq, wq = np.polynomial.legendre.leggauss(64)
    num = [np.sum(wq * zq * np.sin(gg * zq)) / np.sum(wq * np.sin(gg * zq) ** 2) for gg in gam]
    _say("MAT21-19 (IC) printed b_1..4 closed form "
         + " / ".join(f"{v:.12f}" for v in robin_slab_odd_linear_coefficients(gam))
         + "; quadrature " + " / ".join(f"{v:.12f}" for v in num))
    _say(f"MAT21-19 (IC) ratio e(160)/e(40) = {errs[160] / errs[40]:.4e} (bar 0.5)")
    assert errs[160] <= 0.5 * errs[40]


# ---------------------------------------------------------------- (Dec)
def test_decay_to_ambient():
    pts = (np.array([0.0, 0.5 * A_BOX, 0.9 * A_BOX, 0.0]), np.array([0.0, 0.3 * B_BOX, 0.7 * B_BOX, 0.0]),
           np.array([0.5 * C_BOX, -0.5 * C_BOX, 0.9 * C_BOX, C_BOX]))
    late = _theta(*pts, 100.0 * T_C)
    e = np.max(np.abs(late - T_INF)) / GC
    _say(f"MAT21-19 (Dec) max |T(100 t_C) - T_inf|/(g C) = {e:.3e} (bar 1e-12)")
    assert e <= 1e-12


# ---------------------------------------------------------------- (Bal) / (Half)
def _gauss(lo, hi, n):
    xg, wg = np.polynomial.legendre.leggauss(n)
    return 0.5 * (hi - lo) * (xg + 1.0) + lo, 0.5 * (hi - lo) * wg


def _vol_integral(t, zlo):
    xq, wx = _gauss(-A_BOX, A_BOX, 128)
    yq, wy = _gauss(-B_BOX, B_BOX, 128)
    zq, wz = _gauss(zlo, C_BOX, 96 if zlo == 0.0 else 192)
    th = _theta(xq[:, None, None], yq[None, :, None], zq[None, None, :], t) - T_INF
    return RHO_C * float(np.einsum("ijk,i,j,k->", th, wx, wy, wz))


def _half_flux_modes(t):
    """Cumulative heat leaving the half-box z > 0 on [0, t], in mode space:
    Robin on the top and the four half side faces, k∂_zθ through z = 0."""
    gx, cx = _even(A_BOX, N_TERMS)
    gy, cy = _even(B_BOX, N_TERMS)
    gz = robin_slab_odd_eigenvalues(H * C_BOX / K, N_TERMS)
    bz = G_GRAD * C_BOX * robin_slab_odd_linear_coefficients(gz)
    lam = ((gx / A_BOX) ** 2)[:, None, None] + ((gy / B_BOX) ** 2)[None, :, None] \
        + ((gz / C_BOX) ** 2)[None, None, :]
    tau = -np.expm1(-ALPHA * lam * t) / (ALPHA * lam)
    ix, iy = 2 * A_BOX * np.sin(gx) / gx, 2 * B_BOX * np.sin(gy) / gy
    jz = C_BOX * (1.0 - np.cos(gz)) / gz
    amp = cx[:, None, None] * cy[None, :, None] * bz[None, None, :] * tau
    top = H * np.einsum("lmn,l,m,n->", amp, ix, iy, np.sin(gz))
    sx = H * np.einsum("lmn,l,m,n->", amp, 2 * np.cos(gx), iy, jz)
    sy = H * np.einsum("lmn,l,m,n->", amp, ix, 2 * np.cos(gy), jz)
    bot = K * np.einsum("lmn,l,m,n->", amp, ix, iy, gz / C_BOX)
    full = H * (np.einsum("lmn,l,m,n->", amp, ix, iy, np.sin(gz) - np.sin(gz))
                + np.einsum("lmn,l,m,n->", amp, 2 * np.cos(gx), iy, 0.0 * jz)
                + np.einsum("lmn,l,m,n->", amp, ix, 2 * np.cos(gy), 0.0 * jz))
    return top + sx + sy + bot, (top, sx, sy, bot), full


BAL_FOS = (0.05, 0.2, 1.0)


def test_balance_whole_and_half():
    e_half_scale = RHO_C * G_GRAD * C_BOX**2 / 2.0 * 4.0 * A_BOX * B_BOX
    e_scale = RHO_C * GC * 8.0 * A_BOX * B_BOX * C_BOX
    e0 = _vol_integral(0.0, -C_BOX)
    h0 = _vol_integral(0.0, 0.0)
    _say(f"MAT21-19 (Bal) series E(0) whole box {e0:.3e} J (exact 0 by symmetry); "
         f"half box {h0:.6f} J vs exact rho c g C^2/2 4AB {e_half_scale:.6f} J "
         f"(truncation rel {abs(h0 / e_half_scale - 1):.3e}, printed)")
    worst_b = worst_h = 0.0
    for fo in BAL_FOS:
        t = fo * T_C
        et = _vol_integral(t, -C_BOX)
        ht = _vol_integral(t, 0.0)
        fh, parts, ffull = _half_flux_modes(t)
        rb = abs(et + ffull - e0) / e_scale
        rh = abs(ht + fh - h0) / e_half_scale
        worst_b, worst_h = max(worst_b, rb), max(worst_h, rh)
        _say(f"MAT21-19 (Bal) Fo {fo}: whole E(t) {et:.3e} J, wall flux (modes) {ffull:.3e} J; "
             f"residual/(rho c g C V) {rb:.3e}")
        _say(f"MAT21-19 (Half) Fo {fo}: E_half(t) {ht:.6f} J, out top {parts[0]:.6f} / sides x "
             f"{parts[1]:.6f} / sides y {parts[2]:.6f} / through z=0 {parts[3]:.6f} J; "
             f"residual/E_half {rh:.3e}")
    _say(f"MAT21-19 (Bal) max {worst_b:.3e}; (Half) max {worst_h:.3e} (bars 1e-10)")
    assert worst_b <= 1e-10
    assert abs(e0) / e_scale <= 1e-10
    assert worst_h <= 1e-10


# ---------------------------------------------------------------- (U)
def test_uniform_reduction():
    gbar = 1.0
    pts = (np.array([0.0, 0.5 * A_BOX, 0.9 * A_BOX, A_BOX]), np.array([0.0, 0.3 * B_BOX, 0.7 * B_BOX, B_BOX]),
           np.array([0.0, 0.5 * C_BOX, -0.9 * C_BOX, C_BOX]))
    worst = 0.0
    for fo in (0.0,) + FOS:
        t = fo * T_C
        a = robin_box_initial_gradient_temperature(
            *pts, t, half_widths=HALF, k=K, h=H, rho_c=RHO_C, g=0.0, t_inf=T_INF,
            theta_uniform=gbar, n_terms=N_TERMS)
        b = robin_box_transient_temperature(
            *pts, t, half_widths=HALF, k=K, q=0.0, h=H, rho_c=RHO_C, t_inf=T_INF,
            t_init=T_INF + gbar, n_terms=(N_TERMS,) * 3)
        d = np.max(np.abs(a - b))
        worst = max(worst, d)
        _say(f"MAT21-19 (U) Fo {fo}: max |new even path - robin_box_transient(q=0)| {d:.3e} K; "
             f"centre {float(a[0]) - T_INF:.12f} K")
    _say(f"MAT21-19 (U) max {worst:.3e} K (bar 1e-12, gbar 1 K)")
    assert worst <= 1e-12


# ---------------------------------------------------------------- printed records
def test_printed_records():
    probes = ((0.0, 0.0, 0.0), (0.0, 0.0, C_BOX), (0.0, 0.0, -C_BOX), (0.0, 0.0, 0.5 * C_BOX))
    for fo in (0.05, 0.2, 1.0):
        vals = [float(_theta(*p, fo * T_C) - T_INF) for p in probes]
        _say(f"MAT21-19 printed theta at Fo {fo}: centre {vals[0]:.3e}, (0,0,C) {vals[1]:.6f}, "
             f"(0,0,-C) {vals[2]:.6f}, (0,0,C/2) {vals[3]:.6f} K (g C = {GC:.3f} K)")
    g1 = robin_slab_odd_eigenvalues(H * C_BOX / K, 1)[0]
    gx1 = robin_slab_eigenvalues(H * A_BOX / K, 1)[0]
    gy1 = robin_slab_eigenvalues(H * B_BOX / K, 1)[0]
    lam = (gx1 / A_BOX) ** 2 + (gy1 / B_BOX) ** 2 + (g1 / C_BOX) ** 2
    _say(f"MAT21-19 printed slowest odd mode Lambda {lam:.3f} m^-2, tau {1 / (ALPHA * lam):.0f} s; "
         f"g = {G_GRAD} K/m (fixture choice, not measured)")


# ---------------------------------------------------------------- negative control
def test_negative_control_dirichlet_z():
    rz = np.arange(1, N_TERMS + 1) * np.pi
    s = np.linspace(-0.9, 0.9, 7)
    X, Y, Z = np.meshgrid(s * A_BOX, s * B_BOX, s * C_BOX, indexing="ij")
    rob = _theta(X, Y, Z, T_STAR)
    dz = _theta(X, Y, Z, T_STAR, roots_z=rz)
    miss = np.max(np.abs(rob - dz)) / GC
    _say(f"MAT21-19 control (Dz) Dirichlet odd z roots: max miss at t* /(g C) {miss:.3e} "
         f"(bar > 1e-3); max |theta_Robin| {np.max(np.abs(rob - T_INF)):.4f} K")
    assert miss > 1e-3
