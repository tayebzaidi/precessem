"""IMP-5 step 0 (§9 item 210): dielectric-coated PEC circular cylinder in
IMP-4's lossy host at finite k, closed form, numpy only.

Brief: docs/planning/chunks/IMP-5.md, "§9 item 210". PEC radius a = 1 mm,
lossless εr 3 coating to b = 1.5 mm, host IMP-1's (σ 0.57 S/m, εr 78,
64 MHz, k_host = K_FIXTURE), c = 10 mm. k_coat = K0·√3. Two single-mode
drives: TM_z n = 0 and TE_z n = 1 (E0 ŷ on the axis). e^{+jωt}, H^(2)
outgoing in the host.

Anchors (asserted): (Z) scipy Y at complex z against the ascending series,
(B1) Bessel's equation per region (Jin E.1), (B2) Wronskians (Jin E.19),
(P) PEC traces at ρ = a, (IF) interface continuity at ρ = b (Jin 1.46–1.48),
(I) wire current, (Lε) coating → host and (Lb) b → a limits onto IMP-4's
series, (Q) independent quasi-static TE form, (L) (3, N) adapter layout.
Controls (c1)–(c3) are predicted (rule (e)); each is asserted only as
missing its anchor's band, measured printed beside predicted.
"""

from __future__ import annotations

import numpy as np
from mpi4py import MPI
from scipy import special

from fem_em_solver.utils.analytical import (
    _coated_cylinder_coefficients,
    _coated_cylinder_radial,
    coated_pec_cylinder_te1_field,
    coated_pec_cylinder_tm0_current,
    coated_pec_cylinder_tm0_field,
    pec_cylinder_te1_field,
    pec_cylinder_tm0_current,
    pec_cylinder_tm0_field,
)
from tests.validation.test_imp1_pec_sphere_lossy_host import (
    E0,
    FREQUENCY_HZ,
    TH15_BETA_BAND,
)
from tests.validation.test_imp1_step3_exact_comparand import K0, K_FIXTURE

A = 1.0e-3  # PEC radius (IMP-4)
B = 1.5e-3  # coating outer radius
C = 10.0e-3  # outer radius (IMP-4)
K_HOST = K_FIXTURE
K_COAT = K0 * np.sqrt(3.0)  # lossless εr 3, μ0
R_EPS = K_COAT**2 / K_HOST**2  # ε_c,coat / ε_c,host
OMEGA = 2.0 * np.pi * FREQUENCY_HZ
MU_0 = 4 * np.pi * 1e-7
N_RAD = 64
N_RING = 64
S_STATIC = 1e-6

Z_BAND = 1e-13
B1_BAND = 1e-12
B2_BAND = 1e-13
P_BAND = 1e-14
IF_BAND = 1e-12
I_BAND = 1e-12
L_BAND = 1e-12
Q_BAND = 1e-9
CONTROL_MISS = 1e-3

KW = dict(k_host=K_HOST, k_coat=K_COAT, radius=A, coat_radius=B)


def _say(msg):
    if MPI.COMM_WORLD.rank == 0:
        print(msg, flush=True)


def _ring(rho, n, z=0.0):
    phi = 2.0 * np.pi * (np.arange(n) + 0.25) / n
    return np.column_stack([rho * np.cos(phi), rho * np.sin(phi), np.full(n, z)]), phi


def _polar_components(e, phi):
    e_rho = np.cos(phi) * e[:, 0] + np.sin(phi) * e[:, 1]
    e_phi = -np.sin(phi) * e[:, 0] + np.cos(phi) * e[:, 1]
    return e_rho, e_phi, e[:, 2]


def _rel_pair(lhs, rhs):
    """max |lhs − rhs| / (|lhs| + |rhs|): the scale before cancellation.
    Points where both sides are exactly zero contribute 0."""
    lhs, rhs = np.asarray(lhs), np.asarray(rhs)
    scale = np.abs(lhs) + np.abs(rhs)
    diff = np.abs(lhs - rhs)
    if np.any((scale == 0) & (diff != 0)):
        return np.inf
    with np.errstate(invalid="ignore", divide="ignore"):
        ratio = np.where(scale > 0, diff / np.where(scale > 0, scale, 1.0), 0.0)
    return float(np.max(ratio))


def _bessel_residual(f, fp, fpp, rho, k, n):
    res = fpp + fp / rho + (k**2 - n**2 / rho**2) * f
    scale = np.abs(fpp) + np.abs(fp / rho) + np.abs((k**2 - n**2 / rho**2) * f)
    return np.max(np.abs(res) / scale)


def _te(points, **over):
    kw = dict(KW, omega=OMEGA, e0=E0)
    kw.update(over)
    return coated_pec_cylinder_te1_field(points, **kw)


def _tm(points, **over):
    kw = dict(KW, e0=E0)
    kw.update(over)
    return coated_pec_cylinder_tm0_field(points, **kw)


def _incident(points, drive, k=K_HOST):
    """Incident field in the host (IMP-4's normalisation)."""
    pts = np.asarray(points, float)
    rho = np.hypot(pts[:, 0], pts[:, 1])
    cos_p, sin_p = pts[:, 0] / rho, pts[:, 1] / rho
    out = np.zeros((pts.shape[0], 3), complex)
    if drive == "tm0":
        out[:, 2] = E0 * special.jv(0, k * rho)
        return out
    pref = 2.0 * E0 / k
    e_rho = pref * special.jv(1, k * rho) * sin_p / rho
    e_phi = pref * k * special.jvp(1, k * rho, 1) * cos_p
    out[:, 0] = e_rho * cos_p - e_phi * sin_p
    out[:, 1] = e_rho * sin_p + e_phi * cos_p
    return out


def _probe_points():
    pts = [_ring(m * B, 24)[0] for m in (1.2, 1.5, 2.0)]
    return np.vstack(pts)


def _lsq_beta(basis, target):
    """Complex LSQ amplitude β minimising |β·basis − target| over all
    components (IMP-4's fit)."""
    b, t = basis.ravel(), target.ravel()
    return complex(np.vdot(b, t) / np.vdot(b, b))


def _annulus_mean(field_fn, r1, r2, n_rad=64, n_phi=64):
    """⟨|E|²⟩/E0² over r1 < ρ < r2: Gauss–Legendre in ρ, uniform in φ
    (exact for the cos²/sin² content)."""
    x, w = np.polynomial.legendre.leggauss(n_rad)
    rho = 0.5 * (r2 - r1) * x + 0.5 * (r2 + r1)
    wr = 0.5 * (r2 - r1) * w
    phi = 2.0 * np.pi * (np.arange(n_phi) + 0.5) / n_phi
    rr, pp = np.meshgrid(rho, phi, indexing="ij")
    pts = np.column_stack([(rr * np.cos(pp)).ravel(), (rr * np.sin(pp)).ravel(), np.zeros(rr.size)])
    e2 = np.sum(np.abs(field_fn(pts)) ** 2, axis=1).reshape(rr.shape)
    integral = np.sum(e2 * rr * wr[:, None]) * (2.0 * np.pi / n_phi)
    return integral / (np.pi * (r2**2 - r1**2)) / E0**2


def _y_series(n, z, terms=30):
    """Ascending series of Y_n (n = 0, 1) for complex z (A&S 9.1.11):
    Y_n = −(1/π)(z/2)^−n Σ_{k<n} (n−k−1)!/k! (z²/4)^k + (2/π) ln(z/2) J_n
          − (1/π)(z/2)^n Σ_k [ψ(k+1) + ψ(n+k+1)] (−z²/4)^k / (k!(n+k)!)."""
    j_n = sum(
        (-1) ** m / (special.factorial(m) * special.factorial(m + n)) * (z / 2) ** (2 * m + n)
        for m in range(terms)
    )
    head = 0.0
    for k in range(n):
        head += special.factorial(n - k - 1) / special.factorial(k) * (z**2 / 4) ** k
    head = -(1.0 / np.pi) * (z / 2) ** (-n) * head if n else 0.0
    tail = sum(
        (special.digamma(k + 1) + special.digamma(n + k + 1))
        * (-(z**2) / 4) ** k
        / (special.factorial(k) * special.factorial(n + k))
        for k in range(terms)
    )
    return head + (2.0 / np.pi) * np.log(z / 2) * j_n - (1.0 / np.pi) * (z / 2) ** n * tail


# ----------------------------------------------------------------- anchors


def test_z_scipy_y_complex_arguments():
    worst = 0.0
    for z in (K_HOST * A, K_HOST * C):
        y0, y1 = _y_series(0, z), _y_series(1, z)
        refs = {
            ("yv", 0): y0,
            ("yv", 1): y1,
            ("yvp", 0): -y1,  # Y0' = −Y1
            ("yvp", 1): y0 - y1 / z,  # Y1' = Y0 − Y1/z
        }
        for (name, n), ref in refs.items():
            val = special.yv(n, z) if name == "yv" else special.yvp(n, z, 1)
            assert np.iscomplexobj(val) and abs(val.imag) > 0
            worst = max(worst, abs(val - ref) / abs(ref))
    _say(
        f"[IMP-5 W1] (Z) scipy Y at complex z vs ascending series: max rel "
        f"{worst:.3e} (band {Z_BAND:g})"
    )
    assert worst <= Z_BAND


def test_printed_fixture_values():
    _say(
        f"[IMP-5 W1] k_host = {K_HOST:.6g}, k_coat = {K_COAT:.6g}, "
        f"|k_coat| b = {abs(K_COAT)*B:.4g} (pred 3.48e-3), r = k1^2/k2^2 = "
        f"{R_EPS:.4g} (pred 7.39e-3+1.515e-2j), |r| = {abs(R_EPS):.4g}"
    )
    # TM
    _, s = _coated_cylinder_coefficients("tm0", **KW)
    s0 = -special.jv(0, K_HOST * A) / special.hankel2(0, K_HOST * A)
    cur = coated_pec_cylinder_tm0_current(**KW, omega=OMEGA, e0=E0)
    cur0 = pec_cylinder_tm0_current(k=K_HOST, radius=A, omega=OMEGA, e0=E0)
    _say(
        f"[IMP-5 W1] TM: s = {s:.10g}, s0 = {s0:.10g}, |s/s0 - 1| = "
        f"{abs(s/s0 - 1):.4g} (pred 3.1e-6); I = {cur:.6g} A per V/m "
        f"(pred 6.9305e-4-2.8348e-3j), |I/I_IMP4 - 1| = {abs(cur/cur0 - 1):.4g} "
        f"(pred 4.3e-5)"
    )
    # TE
    _, t = _coated_cylinder_coefficients("te1", **KW)
    t0 = -special.jvp(1, K_HOST * A, 1) / special.h2vp(1, K_HOST * A, 1)
    q = A**2 / B**2
    delta = (1 - q) + R_EPS * (1 + q)
    b_qs = B**2 * (R_EPS * (1 + q) - (1 - q)) / delta
    c_qs = -2.0 / delta
    _say(
        f"[IMP-5 W1] TE: t/t0 = {t/t0:.5g} (pred -2.158+0.166j); quasi-static "
        f"B/(E0 a^2) = {b_qs/A**2:.5g} (pred -2.159+0.170j); C/E0 = {c_qs:.4g} "
        f"(pred -3.53+0.14j)"
    )
    up = np.array([[0.0, B, 0.0]])
    e_in = _te(up, region="coat")[0, 1]
    e_out = _te(up, region="host")[0, 1]
    e_phi = _te(np.array([[B, 0.0, 0.0]]), region="host")[0, 1]
    _say(
        f"[IMP-5 W1] TE at rho = b: |E_rho(b-)|/|E_rho(b+)| = "
        f"{abs(e_in)/abs(e_out):.4g} (pred 59.4), |E_rho(b+)|/E0 = "
        f"{abs(e_out):.4g} (pred 0.086), |E_phi(b)|/E0 = {abs(e_phi):.4g} "
        f"(pred 1.960)"
    )
    # Annulus means (step 1's comparands)
    te_coat = _annulus_mean(lambda p: _te(p, region="coat"), A, B)
    te_host = _annulus_mean(lambda p: _te(p, region="host"), B, 1.5 * B)
    te_bare = _annulus_mean(
        lambda p: pec_cylinder_te1_field(
            p, k=K_HOST, radius=A, omega=OMEGA, e0=E0, interior="continue"
        ),
        B,
        1.5 * B,
    )
    tm_coat = _annulus_mean(lambda p: _tm(p, region="coat"), A, B)
    tm_host = _annulus_mean(lambda p: _tm(p, region="host"), B, 1.5 * B)
    _say(
        f"[IMP-5 W1] annulus <|E|^2>/E0^2: TE coat {te_coat:.5g} (pred 17.99), "
        f"TE host (b,1.5b) {te_host:.5g} (pred 1.411), bare wire same annulus "
        f"{te_bare:.5g} (pred 1.088), ratio {te_host/te_bare:.4g}; TM coat "
        f"{tm_coat:.5g} (pred 3.640e-3), TM host {tm_host:.5g} (pred 2.296e-2)"
    )
    # Non-vacuity of the TE coating effect (independent of the band choice).
    assert abs(t / t0 - 1) > 0.5


def test_b1_bessel_equation_per_region():
    rho_c = np.linspace(A, B, N_RAD)
    rho_h = np.linspace(B, C, N_RAD)
    worst = {}
    for drive, n in (("tm0", 0), ("te1", 1)):
        worst[f"{drive} coat"] = _bessel_residual(
            *_coated_cylinder_radial(rho_c, drive, "coat", **KW), rho_c, K_COAT, n
        )
        worst[f"{drive} host"] = _bessel_residual(
            *_coated_cylinder_radial(rho_h, drive, "host", **KW), rho_h, K_HOST, n
        )
    _say(
        "[IMP-5 W1] (B1) max residual / term-sum: "
        + ", ".join(f"{key} {v:.3e}" for key, v in worst.items())
        + f" (band {B1_BAND:g})"
    )
    for key, v in worst.items():
        assert v <= B1_BAND, f"(B1) {key}: {v:.3e}"


def test_b2_wronskians():
    worst_y = worst_h = 0.0
    for n in (0, 1, 2):
        for z in (K_COAT * A, K_COAT * B):
            z = complex(z)
            w = special.jv(n, z) * special.yvp(n, z, 1) - special.jvp(n, z, 1) * special.yv(n, z)
            ref = 2.0 / (np.pi * z)
            worst_y = max(worst_y, abs(w - ref) / abs(ref))
        for z in (K_HOST * B, K_HOST * C):
            w = special.jv(n, z) * special.h2vp(n, z, 1) - special.jvp(n, z, 1) * special.hankel2(
                n, z
            )
            ref = -2j / (np.pi * z)
            worst_h = max(worst_h, abs(w - ref) / abs(ref))
    _say(
        f"[IMP-5 W1] (B2) Wronskian rel. error: J/Y (k1a, k1b) {worst_y:.3e}, "
        f"J/H2 (k2b, k2c) {worst_h:.3e} (band {B2_BAND:g})"
    )
    assert worst_y <= B2_BAND and worst_h <= B2_BAND


def _pec_trace(condition="exact"):
    pts, phi = _ring(A, N_RING, z=0.003)
    tm = _tm(pts, region="coat", _condition=condition)
    te = _te(pts, region="coat", _condition=condition)
    _, tm_phi, tm_z = _polar_components(tm, phi)
    te_rho, te_phi, te_z = _polar_components(te, phi)
    vals = {
        "TM phi": np.max(np.abs(tm_phi)),
        "TM z": np.max(np.abs(tm_z)),
        "TE phi": np.max(np.abs(te_phi)),
        "TE z": np.max(np.abs(te_z)),
    }
    return vals, np.max(np.abs(te_rho))


def test_p_pec_traces():
    vals, te_rho = _pec_trace()
    _say(
        "[IMP-5 W1] (P) max tangential |E|/E0 on rho = a (region='coat'): "
        + ", ".join(f"{k} {v:.3e}" for k, v in vals.items())
        + f"; |E_rho| max {te_rho:.4g} (band {P_BAND:g})"
    )
    assert te_rho > 1.0 * E0  # non-vacuity: the normal TE field is O(E0)
    for key, v in vals.items():
        assert v / E0 <= P_BAND, f"(P) {key}: {v:.3e}"


def _interface(condition="exact"):
    pts, phi = _ring(B, N_RING, z=-0.002)
    out = {}
    tm_c = _tm(pts, region="coat", _condition=condition)
    tm_h = _tm(pts, region="host", _condition=condition)
    out["TM E_z"] = _rel_pair(tm_c[:, 2], tm_h[:, 2])
    _, fp_c, _ = _coated_cylinder_radial(
        np.full(N_RING, B), "tm0", "coat", **KW, condition=condition
    )
    _, fp_h, _ = _coated_cylinder_radial(
        np.full(N_RING, B), "tm0", "host", **KW, condition=condition
    )
    out["TM H_phi"] = _rel_pair(E0 * fp_c / (1j * OMEGA * MU_0), E0 * fp_h / (1j * OMEGA * MU_0))
    te_c = _te(pts, region="coat", _condition=condition)
    te_h = _te(pts, region="host", _condition=condition)
    rc, pc, zc = _polar_components(te_c, phi)
    rh, ph, zh = _polar_components(te_h, phi)
    out["TE E_phi"] = _rel_pair(pc, ph)
    out["TE E_z"] = _rel_pair(zc, zh)
    h0 = -2j * OMEGA * (K_HOST**2 / (OMEGA**2 * MU_0)) * E0 / K_HOST
    f_c, _, _ = _coated_cylinder_radial(
        np.full(N_RING, B), "te1", "coat", **KW, condition=condition
    )
    f_h, _, _ = _coated_cylinder_radial(
        np.full(N_RING, B), "te1", "host", **KW, condition=condition
    )
    out["TE H_z"] = _rel_pair(h0 * f_c * np.cos(phi), h0 * f_h * np.cos(phi))
    out["TE D_n"] = _rel_pair(K_COAT**2 * rc, K_HOST**2 * rh)
    return out, (pc, ph, rc, rh)


def test_if_interface_continuity():
    res, _ = _interface()
    _say(
        "[IMP-5 W1] (IF) rho = b residual / (|coat| + |host|): "
        + ", ".join(f"{k} {v:.3e}" for k, v in res.items())
        + f" (band {IF_BAND:g})"
    )
    for key, v in res.items():
        assert v <= IF_BAND, f"(IF) {key}: {v:.3e}"


def test_i_wire_current():
    _, fp, _ = _coated_cylinder_radial(np.array([A]), "tm0", "coat", **KW)
    from_field = 2.0 * np.pi * A * E0 * fp[0] / (1j * OMEGA * MU_0)
    closed = coated_pec_cylinder_tm0_current(**KW, omega=OMEGA, e0=E0)
    rel = abs(from_field - closed) / abs(closed)
    _say(
        f"[IMP-5 W1] (I) 2*pi*a*H_phi(a) = {from_field:.12g}, closed = "
        f"{closed:.12g}, rel {rel:.3e} (band {I_BAND:g})"
    )
    assert rel <= I_BAND


def _limit_points(n=200, seed=210):
    rng = np.random.default_rng(seed)
    rho = np.sqrt(rng.uniform(A**2, C**2, n))
    rho[0], rho[1], rho[2] = A, B, C
    phi = rng.uniform(0.0, 2.0 * np.pi, n)
    return np.column_stack([rho * np.cos(phi), rho * np.sin(phi), rng.uniform(-C, C, n)])


def _limit_misses(**over):
    pts = _limit_points()
    kw = dict(KW)
    kw.update(over)
    tm = coated_pec_cylinder_tm0_field(pts, **kw, e0=E0)
    tm0 = pec_cylinder_tm0_field(pts, k=K_HOST, radius=A, e0=E0)
    te = coated_pec_cylinder_te1_field(pts, **kw, omega=OMEGA, e0=E0)
    te0 = pec_cylinder_te1_field(pts, k=K_HOST, radius=A, omega=OMEGA, e0=E0)
    cur = coated_pec_cylinder_tm0_current(**kw, omega=OMEGA, e0=E0)
    cur0 = pec_cylinder_tm0_current(k=K_HOST, radius=A, omega=OMEGA, e0=E0)
    c_tm, s = _coated_cylinder_coefficients("tm0", **kw)
    c_te, t = _coated_cylinder_coefficients("te1", **kw)
    s0 = -special.jv(0, K_HOST * A) / special.hankel2(0, K_HOST * A)
    t0 = -special.jvp(1, K_HOST * A, 1) / special.h2vp(1, K_HOST * A, 1)
    return {
        "TM field": np.max(np.abs(tm - tm0)) / np.max(np.abs(tm0)),
        "TE field": np.max(np.abs(te - te0)) / np.max(np.abs(te0)),
        "current": abs(cur - cur0) / abs(cur0),
    }, {
        "TM (c, s)": (c_tm, s),
        "s0": s0,
        "TE (c, t)": (c_te, t),
        "t0": t0,
        "s/s0-1": abs(s / s0 - 1),
        "t/t0-1": abs(t / t0 - 1),
    }


def test_l_eps_coating_to_host():
    miss, coef = _limit_misses(k_coat=K_HOST)
    _say(
        "[IMP-5 W1] (Leps) k_coat = k_host vs IMP-4: "
        + ", ".join(f"{k} {v:.3e}" for k, v in miss.items())
        + f"; |s/s0-1| {coef['s/s0-1']:.3e}, |t/t0-1| {coef['t/t0-1']:.3e} "
        f"(pred ~3e-16, band {L_BAND:g})"
    )
    for key, v in miss.items():
        assert v <= L_BAND, f"(Leps) {key}: {v:.3e}; coefficients {coef}"


def test_l_b_thickness_to_zero():
    miss, coef = _limit_misses(coat_radius=A)
    _say(
        "[IMP-5 W1] (Lb) coat_radius = radius vs IMP-4: "
        + ", ".join(f"{k} {v:.3e}" for k, v in miss.items())
        + f"; |s/s0-1| {coef['s/s0-1']:.3e}, |t/t0-1| {coef['t/t0-1']:.3e} "
        f"(band {L_BAND:g})"
    )
    for key, v in miss.items():
        assert v <= L_BAND, f"(Lb) {key}: {v:.3e}; coefficients {coef}"


def _quasi_static_te(pts):
    """Independent 2D potential-theory field of a coated PEC cylinder in a
    uniform E0 ŷ (derivation in the brief / analytical.py block)."""
    rho = np.hypot(pts[:, 0], pts[:, 1])
    phi = np.arctan2(pts[:, 1], pts[:, 0])
    q = A**2 / B**2
    delta = (1 - q) + R_EPS * (1 + q)
    bq = E0 * B**2 * (R_EPS * (1 + q) - (1 - q)) / delta
    cq = -2.0 * E0 / delta
    host = rho >= B
    e_rho = np.where(host, (E0 + bq / rho**2), -cq * (1 + A**2 / rho**2)) * np.sin(phi)
    e_phi = np.where(host, (E0 - bq / rho**2), -cq * (1 - A**2 / rho**2)) * np.cos(phi)
    return np.column_stack(
        [
            e_rho * np.cos(phi) - e_phi * np.sin(phi),
            e_rho * np.sin(phi) + e_phi * np.cos(phi),
            np.zeros(len(rho)),
        ]
    )


def test_q_te_quasi_static():
    # Stratified: 100 points in the coating and 100 in the host (an
    # area-uniform draw over [a, c] put 1 of 200 in the sleeve, W1).
    rng = np.random.default_rng(2101)
    rho = np.concatenate(
        [np.sqrt(rng.uniform(A**2, B**2, 100)), np.sqrt(rng.uniform(B**2, C**2, 100))]
    )
    rho[0], rho[100], rho[199] = A, B, C
    phi = rng.uniform(0.0, 2.0 * np.pi, 200)
    pts = np.column_stack([rho * np.cos(phi), rho * np.sin(phi), np.zeros(200)])
    e = _te(pts, k_host=S_STATIC * K_HOST, k_coat=S_STATIC * K_COAT, interior="continue")
    ref = _quasi_static_te(pts)
    miss = np.max(np.abs(e - ref)) / E0
    rho = np.hypot(pts[:, 0], pts[:, 1])
    _say(
        f"[IMP-5 W1] (Q) s = {S_STATIC:g}: max |E - E_qs|/E0 = {miss:.3e} "
        f"(pred <~1e-11, band {Q_BAND:g}); points coat {int(np.sum(rho < B))}, "
        f"host {int(np.sum(rho >= B))}"
    )
    assert np.sum(rho < B) >= 90 and np.sum(rho >= B) >= 90
    assert miss <= Q_BAND


def test_l_dolfinx_adapter_layout():
    rng = np.random.default_rng(7)
    rho = rng.uniform(A, C, 50)
    phi = rng.uniform(0, 2 * np.pi, 50)
    pts = np.column_stack([rho * np.cos(phi), rho * np.sin(phi), rng.uniform(-C, C, 50)])
    x = np.ascontiguousarray(pts.T)  # dolfinx interpolation layout (3, N)
    for fn in (_tm, _te):
        adapter = lambda xx, fn=fn: fn(xx.T).T  # noqa: E731
        got = adapter(x)
        assert got.shape == (3, 50)
        assert np.array_equal(got, fn(pts).T)


# ---------------------------------------------------------------- controls


def test_c1_te_uncoated_comparand_misses_beta():
    pts = _probe_points()
    inc = _incident(pts, "te1")
    coated_s = _te(pts, region="host") - inc
    bare_s = pec_cylinder_te1_field(pts, k=K_HOST, radius=A, omega=OMEGA, e0=E0) - inc
    beta = _lsq_beta(coated_s, bare_s)
    _say(
        f"[IMP-5 W1] (c1) TE uncoated comparand: beta = {beta:.5g} "
        f"(pred -0.461-0.035j), |beta-1| = {abs(beta-1):.4g} (pred 1.461, "
        f"asserted > {TH15_BETA_BAND:.4g})"
    )
    inc_tm = _incident(pts, "tm0")
    beta_tm = _lsq_beta(
        _tm(pts, region="host") - inc_tm,
        pec_cylinder_tm0_field(pts, k=K_HOST, radius=A, e0=E0) - inc_tm,
    )
    _say(
        f"[IMP-5 W1] (c1, TM record) uncoated comparand |beta-1| = "
        f"{abs(beta_tm-1):.4g} (pred 3.1e-6; not a control)"
    )
    assert abs(beta - 1) > TH15_BETA_BAND


def test_c2_tm_swapped_condition_misses_interface():
    res, _ = _interface(condition="swapped")
    _, s = _coated_cylinder_coefficients("tm0", **KW)
    _, s_w = _coated_cylinder_coefficients("tm0", **KW, condition="swapped")
    pred = abs(R_EPS - 1) / (abs(R_EPS) + 1)
    rho_c = np.linspace(A, B, N_RAD)
    b1 = _bessel_residual(
        *_coated_cylinder_radial(rho_c, "tm0", "coat", **KW, condition="swapped"), rho_c, K_COAT, 0
    )
    p, _ = _pec_trace(condition="swapped")
    _say(
        f"[IMP-5 W1] (c2) TM swapped: (IF) H_phi residual {res['TM H_phi']:.4g} "
        f"(pred {pred:.4g} ~0.98, asserted > {CONTROL_MISS:g}); |s_w/s - 1| = "
        f"{abs(s_w/s - 1):.4g} (pred 0.104); still (B1) {b1:.3e}, (P) TM z "
        f"{p['TM z']:.3e}"
    )
    assert res["TM H_phi"] > CONTROL_MISS


def test_c3_te_swapped_condition_misses_interface():
    res, _ = _interface(condition="swapped")
    _, t = _coated_cylinder_coefficients("te1", **KW)
    _, t_w = _coated_cylinder_coefficients("te1", **KW, condition="swapped")
    beta = t_w / t
    rho_c = np.linspace(A, B, N_RAD)
    b1 = _bessel_residual(
        *_coated_cylinder_radial(rho_c, "te1", "coat", **KW, condition="swapped"), rho_c, K_COAT, 1
    )
    p, _ = _pec_trace(condition="swapped")
    _say(
        f"[IMP-5 W1] (c3) TE swapped: (IF) E_phi residual {res['TE E_phi']:.4g} "
        f"(pred ~0.98, asserted > {CONTROL_MISS:g}); beta = t_w/t = {beta:.5g}, "
        f"|beta-1| = {abs(beta-1):.4g} (pred ~1.461); still (B1) {b1:.3e}, "
        f"(P) TE phi {p['TE phi']:.3e}"
    )
    assert res["TE E_phi"] > CONTROL_MISS
