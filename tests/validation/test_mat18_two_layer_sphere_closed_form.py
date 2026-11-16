"""`MAT-18` step 1 — the two-layer Robin sphere in closed form (§9 item 190).

Numpy only.  Fixture: gel core r < a = 0.1 (k₁ 0.6, ρc₁ 4.18e6, q₁ 2000),
acrylic-like shell a < r < b = 0.11 (k₂ 0.19, ρc₂ 1.7e6, q₂ 0 — representative
parameters, not a material claim), h 10 on r = b, T∞ = T(0) = 20,
t_R = a²ρc₁/k₁.

Anchor: `utils.analytical.two_layer_sphere_*` (derivations in the docstrings).
Asserted: (S1) steady interface / Robin / balance / two reductions; (R0)
roots — |F| residual, Sturm count, analytic Nₙ, weighted Gram; (R1) Green-form
cₙ vs direct quadrature; (R2) transient interface and Robin residuals from
analytic derivatives; (R3) global balance per mode; (R4) Fo 50 limit and the
equal-property reduction to `MAT-9`.  Negative controls (rule (e), all
predicted, asserted only as exceeding their band): (c1) flux jump omitted;
(c2) homogeneous-ball roots; (c3) Gram with a uniform weight.
"""

from __future__ import annotations

import numpy as np
from mpi4py import MPI

from fem_em_solver.utils.analytical import (
    _two_layer_sphere_F,
    _two_layer_sphere_mode_integrals,
    _two_layer_sphere_modes,
    _two_layer_sphere_series,
    _two_layer_sphere_steady,
    robin_sphere_eigenvalues,
    robin_sphere_steady_temperature,
    robin_sphere_transient_temperature,
    two_layer_sphere_eigenvalues,
    two_layer_sphere_steady_temperature,
    two_layer_sphere_transient_temperature,
)

A_, B_ = 0.1, 0.11
K1, RC1, Q1 = 0.6, 4.18e6, 2000.0
K2, RC2, Q2 = 0.19, 1.7e6, 0.0
H, T_INF = 10.0, 20.0
T_R = A_**2 * RC1 / K1
FOS = (0.05, 0.2, 1.0)
GEO = dict(a=A_, b=B_, k1=K1, k2=K2, rho_c1=RC1, rho_c2=RC2)
STEADY = dict(a=A_, b=B_, k1=K1, k2=K2, q1=Q1, q2=Q2, h=H)
TSCALE = Q1 * A_**2 / K1  # q₁a²/k₁ [K]
NT = 400


def _p(*a):
    if MPI.COMM_WORLD.rank == 0:
        print(*a, flush=True)


def _gl(lo, hi, panels, order=40):
    """Composite Gauss–Legendre nodes / weights on [lo, hi]."""
    x, w = np.polynomial.legendre.leggauss(order)
    edges = np.linspace(lo, hi, panels + 1)
    mid = 0.5 * (edges[1:] + edges[:-1])[:, None]
    half = 0.5 * (edges[1:] - edges[:-1])[:, None]
    return (mid + half * x).ravel(), (half * w).ravel()


# Split at r = a (the kink).
RC_, WC_ = _gl(0.0, A_, 400)
RS_, WS_ = _gl(A_, B_, 40)


def _u_core_shell(s, **kw):
    pc, _ = _two_layer_sphere_modes(RC_, s, **GEO, layer="core", **kw)
    ps, _ = _two_layer_sphere_modes(RS_, s, **GEO, layer="shell", **kw)
    return pc * RC_[:, None], ps * RS_[:, None]  # u = rψ


def _quad_N(s):
    uc, us = _u_core_shell(s)
    return RC1 * (WC_ @ uc**2) + RC2 * (WS_ @ us**2)


S = two_layer_sphere_eigenvalues(**GEO, h=H, n=NT)


# ---------------------------------------------------------------------------


def test_s1_steady():
    th = lambda r, **kw: _two_layer_sphere_steady(r, **STEADY, **kw)  # noqa: E731
    jump_t = abs(th(A_, layer="core") - th(A_, layer="shell")) / th(A_, layer="core")
    fc = K1 * th(A_, layer="core", want="dr")
    fs = K2 * th(A_, layer="shell", want="dr")
    jump_f = abs(fc - fs) / abs(fc)
    rob = abs(-K2 * th(B_, want="dr") - H * th(B_)) / (H * th(B_))
    bal = abs(4 * np.pi * B_**2 * H * th(B_) - Q1 * 4 * np.pi * A_**3 / 3) / (Q1 * 4 * np.pi * A_**3 / 3)
    _p(f"[MAT-18 S1] θ jump {jump_t:.3e}, flux jump {jump_f:.3e}, Robin {rob:.3e}, "
       f"balance {bal:.3e} (bands 1e-13)")
    assert jump_t <= 1e-13 and jump_f <= 1e-13 and rob <= 1e-13 and bal <= 1e-13

    r = np.linspace(0.0, B_, 23)
    eq = two_layer_sphere_steady_temperature(r, a=A_, b=B_, k1=K1, k2=K1, q1=Q1, q2=Q1, h=H, t_inf=T_INF)
    ref = robin_sphere_steady_temperature(r, radius=B_, k=K1, q=Q1, h=H, t_inf=T_INF)
    e_eq = float(np.max(np.abs(eq - ref)) / np.max(np.abs(ref - T_INF)))
    bt = A_ * (1 + 1e-9)
    r2 = np.linspace(0.0, A_, 21)
    thin = two_layer_sphere_steady_temperature(r2, a=A_, b=bt, k1=K1, k2=K2, q1=Q1, h=H, t_inf=T_INF)
    ref2 = robin_sphere_steady_temperature(r2, radius=A_, k=K1, q=Q1, h=H, t_inf=T_INF)
    e_thin = float(np.max(np.abs(thin - ref2)) / np.max(np.abs(ref2 - T_INF)))
    _p(f"[MAT-18 S1] equal properties vs robin_sphere(radius=b): {e_eq:.3e} (band 1e-13)")
    _p(f"[MAT-18 S1] b = a(1+1e-9) vs robin_sphere(radius=a): {e_thin:.3e} (band 1e-7, pred ≈ 4e-9)")
    tb, ta, t0 = th(B_), th(A_), th(0.0)
    mat9 = float(robin_sphere_steady_temperature(0.0, radius=A_, k=K1, q=Q1, h=H))
    _p(f"[MAT-18 records] θ(b) {float(tb):.4f} K (pred 5.5096), θ(a) {float(ta):.4f} K "
       f"(pred 8.6994), θ(0) {float(t0):.4f} K (pred 14.2550); MAT-9 no wall {mat9:.4f} K "
       f"(pred 12.2222); wall adds {100 * (t0 - mat9) / t0:.2f} % (pred 14.3 %)")
    assert e_eq <= 1e-13
    assert e_thin <= 1e-7


def test_r0_roots_and_controls():
    x = np.sqrt(S)
    F = _two_layer_sphere_F(x, **GEO, h=H)
    sc = _two_layer_sphere_F(x, **GEO, h=H, want="scale")
    # Scale = constituent magnitudes before cancellation.  W1 (2026-10-03,
    # log 20261003T021446Z_MAT-18-step1-W1.log) used the post-cancellation
    # k₂|u₂′(b)| + |k₂/b − h||u₂(b)| and read 7.770e-11 against this same
    # 1e-12 band: u₂′(b) is ≈ 4e-3 at the roots while its constituents are O(1),
    # so that denominator understated F's rounding floor.  Band unchanged.
    # R0 re-posed (§9 item 199, ruling of the 2026-10-03 03:00 daily review):
    # W1b read 1.823e-11 at n = 364 against 1e-12 over all 400 roots, and its
    # probe (20261003T021825Z_MAT-18-step1-probe.log) measured one ulp of x
    # moving F/scale by 2.74e-11 there — the 1e-12 band sat below one ulp's
    # worth of F above n ≈ 300.  (R0a) keeps 1e-12 on n ≤ 60 (W1b max
    # 2.09e-13); (R0b) certifies all 400 roots to float64 resolution.
    rel = np.abs(F) / sc
    fres = float(np.max(rel[:60]))
    post = np.abs(F) / _two_layer_sphere_F(x, **GEO, h=H, want="scale_post")
    _p(f"[MAT-18 R0a] max |F(sₙ)|/scale over n ≤ 60 = {fres:.3e} at n = "
       f"{int(np.argmax(rel[:60])) + 1} (band 1e-12)")
    _p(f"[MAT-18 R0 record] max |F(sₙ)|/scale over {NT} roots = {np.max(rel):.3e} at n = "
       f"{int(np.argmax(rel)) + 1} (printed only); W1 post-cancellation "
       f"normalisation (printed only) {np.max(post):.3e} at n = {int(np.argmax(post)) + 1}")

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
        Fc = _two_layer_sphere_F(cand, **GEO, h=H)
        sign_ok = (Fc.min(axis=-1) <= 0.0) & (Fc.max(axis=-1) >= 0.0)
        a_ = np.abs(Fc)
        min_ok = a_[..., 2] <= a_.min(axis=-1)
        return sign_ok, min_ok

    sign_ok, min_ok = _r0b(x)
    n_nosign = int(np.sum(~sign_ok))
    n_nonmin = int(np.sum(~min_ok))
    _p(f"[MAT-18 R0b] {NT} roots: no sign change within ±2 ulp: {n_nosign} (must be 0); "
       f"stored double not the ±2-ulp |F|-minimiser: {n_nonmin} (must be 0)")
    bad_b = np.nonzero(~(sign_ok & min_ok))[0]
    iw = int(bad_b[0]) if bad_b.size else int(np.argmax(rel))
    nb = _nbhd(x[iw], 4)
    _p(f"[MAT-18 R0b] ±4-ulp neighbourhood of n = {iw + 1} "
       f"({'first failing root' if bad_b.size else 'worst |F|/scale, record'}): "
       f"x = {np.array2string(nb, precision=17)}, F/scale = "
       f"{np.array2string(_two_layer_sphere_F(nb, **GEO, h=H) / sc[iw], precision=3)}")
    if bad_b.size:
        _p(f"[MAT-18 R0b] failing n: {(bad_b + 1)[:20].tolist()}")
    xp = x.copy()
    for _ in range(8):
        xp = np.nextafter(xp, np.float64(np.inf))
    _, min_ok_p = _r0b(xp)
    n_ctl = int(np.sum(~min_ok_p))
    _p(f"[MAT-18 R0b control] roots +8 ulp failing the ±2-ulp minimiser check: {n_ctl} "
       f"of {NT} (pred ≥ 1; asserted > 0)")

    # Sturm count, n ≤ 60.
    rg = np.linspace(0.0, B_, 200001)[1:]
    psi, _ = _two_layer_sphere_modes(rg, S[:60], **GEO)
    changes = np.sum(np.sign(psi[:-1]) * np.sign(psi[1:]) < 0, axis=0)
    bad = [n + 1 for n in range(60) if changes[n] != n]
    _p(f"[MAT-18 R0] Sturm sign changes n=1..60 match n−1: {not bad} (mismatches {bad[:5]})")

    N, _, _ = _two_layer_sphere_mode_integrals(S, **GEO)
    Nq = _quad_N(S)
    eN = float(np.max(np.abs(N - Nq)) / np.max(N))
    _p(f"[MAT-18 R0] Nₙ analytic vs split GL ({NT} modes): {eN:.3e} (band 1e-12)")

    uc, us = _u_core_shell(S[:60])
    G = RC1 * (uc.T * WC_) @ uc + RC2 * (us.T * WS_) @ us
    d = np.sqrt(np.diag(G))
    Gn = G / np.outer(d, d)
    off = float(np.max(np.abs(Gn - np.diag(np.diag(Gn)))))
    G3 = RC1 * ((uc.T * WC_) @ uc + (us.T * WS_) @ us)
    d3 = np.sqrt(np.diag(G3))
    off3 = float(np.max(np.abs(G3 / np.outer(d3, d3) - np.eye(60))))
    _p(f"[MAT-18 R0] weighted Gram off-diag max (m≠n≤60) = {off:.3e} (band 1e-10)")
    _p(f"[MAT-18 c3] Gram with uniform weight ρc₁: {off3:.3e} (pred 1e-2…1e-1; asserted > 1e-6)")

    lam_h = robin_sphere_eigenvalues(H * B_ / K1, 60)
    xh = lam_h * np.sqrt(K1 / RC1) / B_
    c2 = np.abs(_two_layer_sphere_F(xh, **GEO, h=H)) / _two_layer_sphere_F(xh, **GEO, h=H, want="scale")
    _p(f"[MAT-18 c2] homogeneous-ball roots in composite F: max {np.max(c2):.3e}, "
       f"median {np.median(c2):.3e} (pred O(1); asserted max > 1e-3)")

    lam_len = A_ * np.sqrt(RC1 / K1) + (B_ - A_) * np.sqrt(RC2 / K2)
    dx = np.diff(x)
    _p(f"[MAT-18 records] Λ = {lam_len:.4f}, π/Λ = {np.pi / lam_len:.6e}; √s spacing "
       f"min {dx.min():.6e} max {dx.max():.6e} mean(last 100) {dx[-100:].mean():.6e}")
    lam_eff_bi = None
    r_tot = (1 / A_ - 1 / B_) / (4 * np.pi * K2) + 1 / (4 * np.pi * B_**2 * H)
    h_eff = 1.0 / (4 * np.pi * A_**2 * r_tot)
    lam_eff_bi = h_eff * A_ / K1
    lam1 = float(robin_sphere_eigenvalues(lam_eff_bi, 1)[0])
    _p(f"[MAT-18 records] effective Bi at a = {lam_eff_bi:.4f} (pred ≈ 1.28), λ₁(Bi_eff) = "
       f"{lam1:.4f} (pred ≈ 1.73), τ₁ = t_R/λ₁² = {T_R / lam1**2:.4e} s (pred ≈ 2.3e4); "
       f"composite 1/s₁ = {1 / S[0]:.4e} s, √(s₁ t_R) = {np.sqrt(S[0] * T_R):.4f}")
    assert fres <= 1e-12
    assert n_nosign == 0 and n_nonmin == 0
    assert n_ctl > 0
    assert not bad
    assert eN <= 1e-12
    assert off <= 1e-10
    assert off3 > 1e-6
    assert np.max(c2) > 1e-3


def test_r1_green_coefficients():
    N, I1, I2 = _two_layer_sphere_mode_integrals(S, **GEO)
    c_green = (Q1 * I1 + Q2 * I2) / (S * N)
    uc, us = _u_core_shell(S)
    thc = _two_layer_sphere_steady(RC_, **STEADY)
    ths = _two_layer_sphere_steady(RS_, **STEADY)
    # ∫ρc θ_ss ψ r² dr = ∫ρc θ_ss u r dr
    c_dir = (RC1 * ((WC_ * thc * RC_) @ uc) + RC2 * ((WS_ * ths * RS_) @ us)) / N
    e = float(np.max(np.abs(c_green - c_dir)) / np.max(np.abs(c_green)))
    _p(f"[MAT-18 R1] cₙ Green vs direct GL ({NT} modes): {e:.3e} (band 1e-10)")
    assert e <= 1e-10


def _series(r, t, **kw):
    return _two_layer_sphere_series(r, t, **GEO, q1=Q1, q2=Q2, h=H, theta_i=0.0, s=S, **kw)


def test_r2_transient_residuals_and_c1():
    worst = {"θ": 0.0, "flux": 0.0, "Robin": 0.0}
    for fo in FOS:
        t = fo * T_R
        thc = _two_layer_sphere_steady(A_, **STEADY, layer="core") + _series(A_, t, layer="core")
        ths = _two_layer_sphere_steady(A_, **STEADY, layer="shell") + _series(A_, t, layer="shell")
        dc = _two_layer_sphere_steady(A_, **STEADY, layer="core", want="dr") + _series(A_, t, layer="core", want="dr")
        ds = _two_layer_sphere_steady(A_, **STEADY, layer="shell", want="dr") + _series(A_, t, layer="shell", want="dr")
        tb = _two_layer_sphere_steady(B_, **STEADY) + _series(B_, t)
        db = _two_layer_sphere_steady(B_, **STEADY, want="dr") + _series(B_, t, want="dr")
        e_t = float(abs(thc - ths)) / TSCALE
        e_f = float(abs(K1 * dc - K2 * ds)) * A_ / K1 / TSCALE
        e_r = float(abs(-K2 * db - H * tb)) * A_ / K1 / TSCALE
        _p(f"[MAT-18 R2] Fo {fo}: θ jump {e_t:.3e}, flux jump {e_f:.3e}, Robin {e_r:.3e} "
           f"(·q₁a²/k₁; band 1e-10)")
        worst["θ"] = max(worst["θ"], e_t)
        worst["flux"] = max(worst["flux"], e_f)
        worst["Robin"] = max(worst["Robin"], e_r)
    # (c1) flux jump omitted: roots, modes and coefficients rebuilt with u₂′(a) = u₁′(a)
    kw = dict(**GEO, q1=Q1, q2=Q2, h=H, theta_i=0.0, n_terms=NT, flux_jump=False)
    c1 = 0.0
    for fo in FOS:
        t = fo * T_R
        dc = _two_layer_sphere_steady(A_, **STEADY, layer="core", want="dr") + _two_layer_sphere_series(A_, t, layer="core", want="dr", **kw)
        ds = _two_layer_sphere_steady(A_, **STEADY, layer="shell", want="dr") + _two_layer_sphere_series(A_, t, layer="shell", want="dr", **kw)
        c1 = max(c1, float(abs(K1 * dc - K2 * ds)) * A_ / K1 / TSCALE)
    _p(f"[MAT-18 c1] flux jump omitted: max flux residual {c1:.3e}·q₁a²/k₁ "
       f"(pred O(|1 − k₂/k₁|) ≈ 0.7; asserted > 1e-3)")
    assert worst["θ"] <= 1e-10
    assert worst["flux"] <= 1e-10
    assert worst["Robin"] <= 1e-10
    assert c1 > 1e-3


def test_r3_global_balance():
    N, I1, I2 = _two_layer_sphere_mode_integrals(S, **GEO)
    cn = (Q1 * I1 + Q2 * I2) / (S * N)
    dn = -cn  # θ_i = 0
    J = RC1 * I1 + RC2 * I2
    psi_b, _ = _two_layer_sphere_modes(B_, S, **GEO)
    src = Q1 * 4 * np.pi * A_**3 / 3 + Q2 * 4 * np.pi * (B_**3 - A_**3) / 3
    worst = 0.0
    for fo in FOS:
        e = np.exp(-S * fo * T_R)
        dEdt = float(np.sum(-S * dn * e * 4 * np.pi * J))
        theta_b = float(_two_layer_sphere_steady(B_, **STEADY)) + float(np.sum(dn * e * psi_b))
        out = 4 * np.pi * B_**2 * H * theta_b
        rel = abs(dEdt + out - src) / src
        _p(f"[MAT-18 R3] Fo {fo}: dE/dt {dEdt:.6e} W + film {out:.6e} W vs source "
           f"{src:.6e} W: rel {rel:.3e} (band 1e-10)")
        worst = max(worst, rel)
    assert worst <= 1e-10


def test_r4_limits():
    r = np.array([0.0, A_ / 2, A_, (A_ + B_) / 2, B_])
    th_ss = _two_layer_sphere_steady(r, **STEADY)
    late = two_layer_sphere_transient_temperature(
        r, 50 * T_R, **GEO, q1=Q1, q2=Q2, h=H, t_inf=T_INF, t_init=T_INF) - T_INF
    e50 = float(np.max(np.abs(late - th_ss)) / th_ss[0])
    _p(f"[MAT-18 R4] Fo 50 vs θ_ss: {e50:.3e}·θ_ss(0) (band 1e-12)")

    geq = dict(a=A_, b=B_, k1=K1, k2=K1, rho_c1=RC1, rho_c2=RC1)
    th0 = float(robin_sphere_steady_temperature(0.0, radius=B_, k=K1, q=Q1, h=H))
    rr = np.linspace(0.0, B_, 12)
    worst = 0.0
    for t_init in (T_INF, T_INF + 5.0):
        for fo in FOS:
            t = fo * T_R
            ours = two_layer_sphere_transient_temperature(
                rr, t, **geq, q1=Q1, q2=Q1, h=H, t_inf=T_INF, t_init=t_init)
            ref = robin_sphere_transient_temperature(
                rr, t, radius=B_, k=K1, rho_c=RC1, q=Q1, h=H, t_inf=T_INF, t_init=t_init)
            e = float(np.max(np.abs(ours - ref))) / th0
            _p(f"[MAT-18 R4] equal props, T_init {t_init}, Fo {fo}: {e:.3e}·θ_ss(0) (band 1e-10)")
            worst = max(worst, e)
    s_eq = two_layer_sphere_eigenvalues(**geq, h=H, n=NT)
    lam_eq = np.sqrt(s_eq * RC1 / K1) * B_
    lam_h = robin_sphere_eigenvalues(H * B_ / K1, NT)
    _p(f"[MAT-18 R4] equal-props roots vs robin_sphere_eigenvalues(hb/k): max |Δλ| "
       f"{np.max(np.abs(lam_eq - lam_h)):.3e} (printed)")

    # t = 0 tail (printed only, never asserted)
    rt = np.array([0.0, A_ / 2, A_, B_])
    for nt in (100, 200, 400):
        v = two_layer_sphere_transient_temperature(
            rt, 0.0, **GEO, q1=Q1, q2=Q2, h=H, t_inf=T_INF, t_init=T_INF, n_terms=nt) - T_INF
        _p(f"[MAT-18 records] t = 0 tail N={nt}: |θ(r,0)| at r = 0, a/2, a, b = "
           f"{np.array2string(np.abs(v), precision=3)} K")
    assert e50 <= 1e-12
    assert worst <= 1e-10
