"""IMP-4 step 0 (§9 item 202): PEC circular cylinder in IMP-1's lossy host
at finite k, closed form, numpy only.

Brief: docs/planning/chunks/IMP-4.md, "§9 item 202". Infinite PEC cylinder of
radius a = 1 mm along z in the IMP-1 host (σ 0.57 S/m, εr 78, 64 MHz,
k = K_FIXTURE), outer radius c = 10 mm. Two single-mode drives:
TM_z n = 0 (axial E) and TE_z n = 1 (transverse E, E0 ŷ on the axis).
Convention e^{+jωt}, H^(2) outgoing.

Anchors (asserted): (B1) Bessel's equation (Jin E.1), (B2) Wronskian
(Jin E.19 for H^(2)), (P) PEC tangential traces on the Cartesian output,
(I) wire current 2πa H_φ(a) = 4E0/(ωμ0 H0^(2)(ka)), (Q) TE static limit,
(L) dolfinx-style (3, N) adapter layout. Negative controls (c1)-(c3) are
predicted (rule (e)); each is asserted only as missing its anchor's band,
measured printed beside predicted.
"""

from __future__ import annotations

import numpy as np
from mpi4py import MPI
from scipy import special

from fem_em_solver.utils.analytical import (
    _pec_cylinder_te1_radial,
    _pec_cylinder_tm0_radial,
    pec_cylinder_te1_field,
    pec_cylinder_tm0_current,
    pec_cylinder_tm0_field,
)
from tests.validation.test_imp1_step3_exact_comparand import (
    EPSILON_R_HOST,
    FREQUENCY_HZ,
    K0,
    K_FIXTURE,
)

A = 1.0e-3  # GEO-36 rod radius
C = 10.0e-3  # step-1 outer radius, |k|c = 0.179
OMEGA = 2.0 * np.pi * FREQUENCY_HZ
MU_0 = 4 * np.pi * 1e-7
E0 = 1.0
N_RAD = 64
S_STATIC = 1e-6

B1_BAND = 1e-12
B2_BAND = 1e-13
P_BAND = 1e-14
I_BAND = 1e-12
Q_BAND = 1e-9
CONTROL_MISS = 1e-3


def _say(msg):
    if MPI.COMM_WORLD.rank == 0:
        print(msg, flush=True)


def _radii():
    return np.linspace(A, C, N_RAD)


def _ring(rho, n, z=0.0):
    phi = 2.0 * np.pi * (np.arange(n) + 0.25) / n
    return np.column_stack([rho * np.cos(phi), rho * np.sin(phi), np.full(n, z)]), phi


def _bessel_residual(f, fp, fpp, rho, k, n):
    res = fpp + fp / rho + (k**2 - n**2 / rho**2) * f
    scale = np.abs(fpp) + np.abs(fp / rho) + np.abs((k**2 - n**2 / rho**2) * f)
    return np.max(np.abs(res) / scale)


def _wire_current_from_field(k, hankel_kind=2):
    _, fp, _ = _pec_cylinder_tm0_radial(np.array([A]), k=k, radius=A, hankel_kind=hankel_kind)
    h_phi = E0 * fp[0] / (1j * OMEGA * MU_0)
    return 2.0 * np.pi * A * h_phi


def test_scipy_accepts_complex_arguments():
    """Check, not assume: scipy's Bessel routines take complex z and agree
    with the ascending series for J0, J1 at ka and kc."""
    for z in (K_FIXTURE * A, K_FIXTURE * C):
        for n in (0, 1):
            series = sum(
                (-1) ** m
                / (special.factorial(m) * special.factorial(m + n))
                * (z / 2) ** (2 * m + n)
                for m in range(30)
            )
            val = special.jv(n, z)
            assert np.iscomplexobj(val) and abs(val.imag) > 0
            assert abs(val - series) <= 1e-14 * abs(series)
        for fn in (special.hankel2, special.jvp, special.h2vp):
            assert np.iscomplexobj(fn(1, z)) and abs(fn(1, z).imag) > 0


def test_printed_fixture_values():
    k = K_FIXTURE
    _say(
        f"[IMP-4 W1] k = {k:.6g} 1/m, |k| = {abs(k):.6g}, |k|a = {abs(k)*A:.6g}, "
        f"|k|c = {abs(k)*C:.6g}"
    )
    coef = special.jv(0, k * A) / special.hankel2(0, k * A)
    _say(
        f"[IMP-4 W1] TM n=0 coefficient J0(ka)/H0(ka) = {coef:.6g}, "
        f"|.| = {abs(coef):.6g} (predicted ~0.37); H0(ka) = "
        f"{special.hankel2(0, k*A):.6g} (predicted ~0.64+2.64j)"
    )
    for m in (1.2, 1.5, 2.0):
        pts = np.array([[m * A, 0.0, 0.0]])
        ez = pec_cylinder_tm0_field(pts, k=k, radius=A, e0=E0)[0, 2]
        _say(f"[IMP-4 W1] |E_z|/E0 at rho = {m}a: {abs(ez)/E0:.6g}")
    pts = np.array([[0.0, A, 0.0]])
    e = pec_cylinder_te1_field(pts, k=k, radius=A, omega=OMEGA, e0=E0)[0]
    _say(
        f"[IMP-4 W1] TE surface |E_rho(a, pi/2)|/E0 = {abs(e[1])/E0:.6g} "
        f"(predicted ~2.00, static doubling)"
    )
    current = pec_cylinder_tm0_current(k=k, radius=A, omega=OMEGA, e0=E0)
    _say(
        f"[IMP-4 W1] |I| = {abs(current):.6g} A per V/m (predicted ~2.9e-3), "
        f"arg I = {np.angle(current):.6g} rad"
    )


def test_b1_bessel_equation():
    rho = _radii()
    k = K_FIXTURE
    worst = {}
    for n in (0, 1):
        z = k * rho
        jn = (special.jv(n, z), k * special.jvp(n, z, 1), k**2 * special.jvp(n, z, 2))
        hn = (special.hankel2(n, z), k * special.h2vp(n, z, 1), k**2 * special.h2vp(n, z, 2))
        worst[f"J{n}"] = _bessel_residual(*jn, rho, k, n)
        worst[f"H{n}"] = _bessel_residual(*hn, rho, k, n)
    worst["f_TM0"] = _bessel_residual(*_pec_cylinder_tm0_radial(rho, k=k, radius=A), rho, k, 0)
    worst["f_TE1"] = _bessel_residual(*_pec_cylinder_te1_radial(rho, k=k, radius=A), rho, k, 1)
    _say(
        "[IMP-4 W1] (B1) max residual / term-sum: "
        + ", ".join(f"{key} {v:.3e}" for key, v in worst.items())
        + f" (band {B1_BAND:g})"
    )
    for key, v in worst.items():
        assert v <= B1_BAND, f"(B1) {key}: {v:.3e}"


def test_b2_wronskian():
    worst = 0.0
    for z in (K_FIXTURE * A, K_FIXTURE * C):
        for n in (0, 1, 2):
            w = special.jv(n, z) * special.h2vp(n, z, 1) - special.jvp(n, z, 1) * special.hankel2(
                n, z
            )
            ref = -2j / (np.pi * z)
            worst = max(worst, abs(w - ref) / abs(ref))
    _say(f"[IMP-4 W1] (B2) max Wronskian rel. error: {worst:.3e} (band {B2_BAND:g})")
    assert worst <= B2_BAND


def _tangential(e, phi):
    e_phi = -np.sin(phi) * e[:, 0] + np.cos(phi) * e[:, 1]
    return np.abs(e_phi), np.abs(e[:, 2])


def test_p_pec_traces_cartesian():
    pts, phi = _ring(A, 64, z=0.003)
    tm = pec_cylinder_tm0_field(pts, k=K_FIXTURE, radius=A, e0=E0, interior="continue")
    te = pec_cylinder_te1_field(pts, k=K_FIXTURE, radius=A, omega=OMEGA, e0=E0, interior="continue")
    tm_phi, tm_z = _tangential(tm, phi)
    te_phi, te_z = _tangential(te, phi)
    worst = max(tm_phi.max(), tm_z.max(), te_phi.max(), te_z.max()) / E0
    _say(
        f"[IMP-4 W1] (P) max tangential |E|/E0 on rho = a: TM phi {tm_phi.max():.3e} "
        f"z {tm_z.max():.3e}; TE phi {te_phi.max():.3e} z {te_z.max():.3e} "
        f"(band {P_BAND:g})"
    )
    # Non-vacuity: the normal TE component is O(E0) there.
    e_rho = np.cos(phi) * te[:, 0] + np.sin(phi) * te[:, 1]
    assert np.max(np.abs(e_rho)) > 1.0 * E0
    assert worst <= P_BAND


def test_i_wire_current():
    from_field = _wire_current_from_field(K_FIXTURE)
    closed = pec_cylinder_tm0_current(k=K_FIXTURE, radius=A, omega=OMEGA, e0=E0)
    rel = abs(from_field - closed) / abs(closed)
    _say(
        f"[IMP-4 W1] (I) 2*pi*a*H_phi(a) = {from_field:.12g}, closed = {closed:.12g}, "
        f"ratio = {from_field/closed:.15g}, rel {rel:.3e} (band {I_BAND:g})"
    )
    assert rel <= I_BAND


def test_q_te_static_limit():
    k = S_STATIC * K_FIXTURE
    rng = np.random.default_rng(202)
    rho = np.sqrt(rng.uniform(A**2, C**2, 200))
    rho[0], rho[1] = A, C
    phi = rng.uniform(0.0, 2.0 * np.pi, 200)
    pts = np.column_stack([rho * np.cos(phi), rho * np.sin(phi), np.zeros(200)])
    e = pec_cylinder_te1_field(pts, k=k, radius=A, omega=OMEGA, e0=E0, interior="continue")
    e_rho_s = E0 * (1 + A**2 / rho**2) * np.sin(phi)
    e_phi_s = E0 * (1 - A**2 / rho**2) * np.cos(phi)
    ref = np.column_stack(
        [
            e_rho_s * np.cos(phi) - e_phi_s * np.sin(phi),
            e_rho_s * np.sin(phi) + e_phi_s * np.cos(phi),
            np.zeros(200),
        ]
    )
    miss = np.max(np.abs(e - ref)) / E0
    _say(
        f"[IMP-4 W1] (Q) s = {S_STATIC:g}: max |E - E_static|/E0 = {miss:.3e} "
        f"(predicted <~1e-11, band {Q_BAND:g})"
    )
    assert miss <= Q_BAND


def test_l_dolfinx_adapter_layout():
    rng = np.random.default_rng(7)
    rho = rng.uniform(A, C, 50)
    phi = rng.uniform(0, 2 * np.pi, 50)
    pts = np.column_stack([rho * np.cos(phi), rho * np.sin(phi), rng.uniform(-C, C, 50)])
    x = np.ascontiguousarray(pts.T)  # dolfinx interpolation layout (3, N)
    for fn in (
        lambda p: pec_cylinder_tm0_field(p, k=K_FIXTURE, radius=A, e0=E0),
        lambda p: pec_cylinder_te1_field(p, k=K_FIXTURE, radius=A, omega=OMEGA, e0=E0),
    ):
        adapter = lambda xx, fn=fn: fn(xx.T).T  # noqa: E731
        got = adapter(x)
        assert got.shape == (3, 50)
        assert np.array_equal(got, fn(pts).T)


def test_c1_hankel1_misses_current():
    wrong = _wire_current_from_field(K_FIXTURE, hankel_kind=1)
    closed = pec_cylinder_tm0_current(k=K_FIXTURE, radius=A, omega=OMEGA, e0=E0)
    miss = abs(wrong - closed) / abs(closed)
    _say(
        f"[IMP-4 W1] (c1) H^(1) in TM vs (I): miss {miss:.4g} (predicted ~0.7, "
        f"asserted > {CONTROL_MISS:g})"
    )
    # (B1) and (P) still pass for the wrong-convention field.
    rho = _radii()
    b1 = _bessel_residual(
        *_pec_cylinder_tm0_radial(rho, k=K_FIXTURE, radius=A, hankel_kind=1), rho, K_FIXTURE, 0
    )
    pts, _ = _ring(A, 64)
    p = np.max(
        np.abs(
            pec_cylinder_tm0_field(
                pts, k=K_FIXTURE, radius=A, e0=E0, interior="continue", _hankel_kind=1
            )[:, 2]
        )
    )
    _say(f"[IMP-4 W1] (c1) wrong-convention field still: (B1) {b1:.3e}, (P) {p:.3e}")
    assert miss > CONTROL_MISS


def test_c2_dirichlet_te_misses_pec_trace():
    pts, phi = _ring(A, 64)
    e = pec_cylinder_te1_field(
        pts, k=K_FIXTURE, radius=A, omega=OMEGA, e0=E0, interior="continue", _condition="dirichlet"
    )
    e_phi = np.max(np.abs(-np.sin(phi) * e[:, 0] + np.cos(phi) * e[:, 1])) / E0
    _say(
        f"[IMP-4 W1] (c2) Dirichlet TE vs (P): max |E_phi(a)|/E0 = {e_phi:.4g} "
        f"(predicted O(1), asserted > {CONTROL_MISS:g})"
    )
    assert e_phi > CONTROL_MISS


def test_c3_sigma_dropped_misses_current():
    k_lossless = K0 * np.sqrt(EPSILON_R_HOST)
    wrong = pec_cylinder_tm0_current(k=k_lossless, radius=A, omega=OMEGA, e0=E0)
    ref = _wire_current_from_field(K_FIXTURE)
    miss = abs(wrong - ref) / abs(ref)
    _say(
        f"[IMP-4 W1] (c3) sigma dropped (k = {k_lossless:.6g}) vs (I): miss "
        f"{miss:.4g} (predicted ~0.15, asserted > {CONTROL_MISS:g})"
    )
    assert miss > CONTROL_MISS
