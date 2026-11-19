"""`IMP-7` step 0 (§9 item 579, §10 I37): steady heating of a Robin-cooled gel
shell by the quasi-static dipole SAR of `IMP-1`'s PEC sphere, closed form,
numpy only.

Brief: docs/planning/chunks/IMP-7.md, "§9 item 579 — 2026-10-22", the 10-20
draft above it, the 2026-10-21 weekly ruling and the 2026-10-22 03:00 review
ruling.  Fixture (rule (f) parameters): a = `IMP-1`'s 5 mm, b = 4a, gel
k₂ = 0.6 and q_far = 2000 (`MAT-9`), h = 100 (`MAT-10`), T∞ = 20, core
k₁ = 20 (metal-like, representative — not a material claim) plus the exact
adiabatic (k₁ = 0) and isothermal (k₁ = ∞) branches.  θ = R₀ + R₂P₂, exact.

Asserted: (Src), (PDE), (I), (Rb), (FD) ratio, (Bal), (Φ), (R0).  Printed
only: (FD) level (review ruling 2026-10-22 03:00), (Rk) (weekly amendment 1),
and every prediction.  Negative controls (c1)–(c3) are predicted (rule (e))
and asserted only as missing their anchor by > `CONTROL_MISS`.
"""

from __future__ import annotations

import numpy as np
from mpi4py import MPI

from fem_em_solver.utils.analytical import (
    _pec_sphere_implant_terms,
    dielectric_sphere_quasistatic_field,
    pec_sphere_implant_coefficients,
    pec_sphere_implant_radial,
    pec_sphere_implant_source,
    pec_sphere_implant_steady_temperature,
    two_layer_sphere_steady_temperature,
)
from tests.validation.test_imp1_pec_sphere_lossy_host import (
    PHI_CLOSED_FORM_PEC,
    SHELL_OUTER_FACTOR,
    SPHERE_RADIUS,
    _pec_exterior_numpy,
)
from tests.validation.test_imp4_pec_cylinder_closed_form import B1_BAND, CONTROL_MISS
from tests.validation.test_mat18_two_layer_sphere_closed_form import _gl

# Bands that exist only as literals in `assert` lines: restated by value, not widened.
R0_BAND = 1e-13  # MAT-18 (S1): test_mat18_two_layer_sphere_closed_form.py:113 (also :95)
FD_RATIO_LO, FD_RATIO_HI = 3.8, 4.2  # IMP-1 (H): test_imp1_step3_exact_comparand.py:193
FD_LEVEL_PRINTED = 1e-4  # IMP-1 (H): test_imp1_step3_exact_comparand.py:194 — printed only
PHI_BAND = 1e-10  # IMP-1 (Φ0): test_imp1_step3_exact_comparand.py:217

# Fixture (rule (f) parameters, not bands).
A = SPHERE_RADIUS  # 0.005, test_imp1_pec_sphere_lossy_host.py:61
B = 4.0 * A
K_GEL = 0.6  # MAT-9 K, test_mat9_steady_conduction_sphere.py:35
Q_FAR = 2000.0  # MAT-9 Q, test_mat9_steady_conduction_sphere.py:36
H = 100.0  # MAT-10 H, test_mat10_sar_to_heat_sphere.py:53
T_INF = 20.0
K_CORE = 20.0
BETA = 1.0
BRANCHES = {"adiabatic": 0.0, "k1=20": K_CORE, "isothermal": np.inf}
FIX = dict(radius=A, outer_radius=B, k_gel=K_GEL, q_far=Q_FAR, h=H)
N_SRC, N_RAD, N_FD = 512, 64, 32
DELTA = 1e-2 * A  # rule (f), IMP-2 (F) ruling


def _p(*args):
    if MPI.COMM_WORLD.rank == 0:
        print(*args, flush=True)


def _coeffs(k_core, beta=BETA):
    return pec_sphere_implant_coefficients(k_core=k_core, beta=beta, **FIX)


def _shell_points(n, rlo, rhi, seed):
    rng = np.random.default_rng(seed)
    r = rng.uniform(rlo, rhi, n)
    mu = rng.uniform(-1.0, 1.0, n)
    phi = rng.uniform(0.0, 2.0 * np.pi, n)
    s = np.sqrt(1.0 - mu * mu)
    return np.column_stack([r * s * np.cos(phi), r * s * np.sin(phi), r * mu])


def _q_from_field(field, pts):
    e = np.real(field(pts.T))  # (3, N) convention
    return Q_FAR * np.sum(e * e, axis=0)


def _q_l(r, l, beta=BETA):
    u = (A / r) ** 3
    if l == 0:
        return Q_FAR * (1.0 + 2.0 * beta**2 * u**2)
    return Q_FAR * (4.0 * beta * u + 2.0 * beta**2 * u**2)


def _term_scale(c, l, r, derivative):
    """Σ|c n^(d) r^(n−d)| — the magnitude before cancellation."""
    s = 0.0
    for cc, n in _pec_sphere_implant_terms(c, l):
        fac = 1.0 if derivative == 0 else (n if derivative == 1 else n * (n - 1))
        s = s + np.abs(cc * fac * r ** (n - derivative))
    return s


def _gl_power(source_on_points, rlo, rhi):
    r, wr = _gl(rlo, rhi, 40, 40)
    mu, wm = np.polynomial.legendre.leggauss(40)
    R, M = np.meshgrid(r, mu, indexing="ij")
    W = np.outer(wr, wm)
    pts = np.column_stack([(R * np.sqrt(1.0 - M * M)).ravel(), np.zeros(R.size), (R * M).ravel()])
    q = source_on_points(pts).reshape(R.shape)
    return 2.0 * np.pi * float(np.sum(W * q * R * R))


# ---------------------------------------------------------------------------


def test_imp7_0_src_source_vs_imp1_field():
    pts = _shell_points(N_SRC, 1.001 * A, B, seed=579)
    q = pec_sphere_implant_source(pts, radius=A, q_far=Q_FAR, beta=1.0)
    ref = _q_from_field(_pec_exterior_numpy(), pts)
    e1 = float(np.max(np.abs(q - ref)) / np.max(q))
    errs = {}
    for beta, eps in ((-0.5, 0.0), (0.5, 4.0)):
        qb = pec_sphere_implant_source(pts, radius=A, q_far=Q_FAR, beta=beta)
        rb = _q_from_field(dielectric_sphere_quasistatic_field(radius=A, epsilon_c=eps), pts)
        errs[beta] = float(np.max(np.abs(qb - rb)) / np.max(qb))
    u = (A / np.linalg.norm(pts, axis=1)) ** 3
    q_iso = Q_FAR * (1.0 + 2.0 * u**2)
    c1 = float(np.max(np.abs(q_iso - q)) / np.max(q))
    _p(
        f"\n[IMP7-0 Src] beta=1 vs IMP-1 _pec_exterior_numpy: {e1:.3e}; beta=-1/2 (eps 0): "
        f"{errs[-0.5]:.3e}; beta=1/2 (eps 4): {errs[0.5]:.3e} (band B1_BAND {B1_BAND:.0e} x max q)"
    )
    _p(
        f"[IMP7-0 c1] isotropic source vs (Src): {c1:.4f} (predicted 0.6667 ceiling; "
        f"asserted > CONTROL_MISS {CONTROL_MISS:.0e})"
    )
    assert e1 <= B1_BAND
    assert errs[-0.5] <= B1_BAND and errs[0.5] <= B1_BAND
    assert c1 > CONTROL_MISS


def test_imp7_0_pde_radial_residual():
    r = np.linspace(A, B, N_RAD + 1)[1:]
    worst = 0.0
    for name, k1 in BRANCHES.items():
        c = _coeffs(k1)
        for l in (0, 2):
            R = pec_sphere_implant_radial(r, c, l)
            R1 = pec_sphere_implant_radial(r, c, l, derivative=1)
            R2 = pec_sphere_implant_radial(r, c, l, derivative=2)
            ql = _q_l(r, l)
            res = -K_GEL * (R2 + 2.0 * R1 / r - l * (l + 1) * R / r**2) - ql
            scale = K_GEL * (
                _term_scale(c, l, r, 2)
                + 2.0 * _term_scale(c, l, r, 1) / r
                + l * (l + 1) * _term_scale(c, l, r, 0) / r**2
            ) + np.abs(ql)
            e = float(np.max(np.abs(res) / scale))
            worst = max(worst, e)
            _p(f"[IMP7-0 PDE] {name:10s} l={l}: max residual / term sum {e:.3e}")
    _p(f"[IMP7-0 PDE] worst {worst:.3e} (band B1_BAND {B1_BAND:.0e})")
    assert worst <= B1_BAND


def test_imp7_0_interface_robin_and_c3():
    worst_i, worst_b = 0.0, 0.0
    for name, k1 in BRANCHES.items():
        c = _coeffs(k1)
        # l = 0: the core constant carries no flux, so k₂R₀′(a⁺) = 0 on every branch.
        f0 = K_GEL * pec_sphere_implant_radial(A, c, 0, derivative=1)
        e_f0 = abs(float(f0)) / float(K_GEL * _term_scale(c, 0, A, 1))
        jumps = [("l=0 flux", e_f0)]
        if k1 != 0.0:
            t0 = c["theta0_core"] - pec_sphere_implant_radial(A, c, 0)
            jumps.append(
                (
                    "l=0 theta",
                    abs(float(t0)) / float(abs(c["theta0_core"]) + _term_scale(c, 0, A, 0)),
                )
            )
            t2 = c["F"] * A**2 - pec_sphere_implant_radial(A, c, 2)
            jumps.append(
                ("l=2 theta", abs(float(t2)) / float(abs(c["F"] * A**2) + _term_scale(c, 2, A, 0)))
            )
        if np.isfinite(k1):
            core_flux = 0.0 if k1 == 0.0 else k1 * 2.0 * c["F"] * A
            f2 = core_flux - K_GEL * pec_sphere_implant_radial(A, c, 2, derivative=1)
            jumps.append(
                (
                    "l=2 flux",
                    abs(float(f2)) / float(abs(core_flux) + K_GEL * _term_scale(c, 2, A, 1)),
                )
            )
        for lab, e in jumps:
            worst_i = max(worst_i, e)
            _p(f"[IMP7-0 I] {name:10s} {lab}: {e:.3e}")
        for l in (0, 2):
            res = K_GEL * pec_sphere_implant_radial(
                B, c, l, derivative=1
            ) + H * pec_sphere_implant_radial(B, c, l)
            sc = K_GEL * _term_scale(c, l, B, 1) + H * _term_scale(c, l, B, 0)
            e = abs(float(res)) / float(sc)
            worst_b = max(worst_b, e)
            _p(f"[IMP7-0 Rb] {name:10s} l={l}: {e:.3e}")
    # (c3): the adiabatic R₂ used for k₁ = K_CORE, against (I)'s l = 2 flux jump.
    ca = _coeffs(0.0)
    r2a = float(pec_sphere_implant_radial(A, ca, 2))
    core_flux = 2.0 * K_CORE * r2a / A
    shell_flux = float(K_GEL * pec_sphere_implant_radial(A, ca, 2, derivative=1))
    c3 = abs(core_flux - shell_flux) / (abs(core_flux) + K_GEL * float(_term_scale(ca, 2, A, 1)))
    _p(
        f"[IMP7-0 I] worst {worst_i:.3e}; [IMP7-0 Rb] worst {worst_b:.3e} (band B1_BAND {B1_BAND:.0e})"
    )
    _p(
        f"[IMP7-0 c3] adiabatic R2 at k1 = {K_CORE}: flux jump {c3:.4f} (predicted 1.0; "
        f"note: the term-sum scale includes the shell's pre-cancellation terms, so the "
        f"reading may sit below 1; asserted > CONTROL_MISS {CONTROL_MISS:.0e})"
    )
    assert worst_i <= B1_BAND
    assert worst_b <= B1_BAND
    assert c3 > CONTROL_MISS


def _fd_residual(pts, delta):
    kw = dict(k_core=K_CORE, t_inf=T_INF, beta=BETA, **FIX)
    T0 = pec_sphere_implant_steady_temperature(pts, **kw)
    lap = -6.0 * T0
    for ax in range(3):
        for s in (1.0, -1.0):
            p = pts.copy()
            p[:, ax] += s * delta
            lap = lap + pec_sphere_implant_steady_temperature(p, **kw)
    lap = lap / delta**2
    q = pec_sphere_implant_source(pts, radius=A, q_far=Q_FAR, beta=BETA)
    res = np.abs(-K_GEL * lap - q) / (9.0 * Q_FAR)
    return float(np.max(res))


def test_imp7_0_fd_cartesian_laplacian():
    pts = _shell_points(N_FD, 1.1 * A, 0.95 * B, seed=5790)
    e1 = _fd_residual(pts, DELTA)
    e2 = _fd_residual(pts, 2.0 * DELTA)
    ratio = e2 / e1
    _p(
        f"\n[IMP7-0 FD] max|-k2 Lap_h T - q|/(9 q_far), 32 shell points: delta = 1e-2 a "
        f"{e1:.3e} (predicted ~1.1e-5; IMP-1 (H) level {FD_LEVEL_PRINTED:.0e} printed only); "
        f"2 delta {e2:.3e}; ratio {ratio:.3f} (predicted 4.00, asserted in "
        f"[{FD_RATIO_LO}, {FD_RATIO_HI}])"
    )
    assert FD_RATIO_LO <= ratio <= FD_RATIO_HI


def test_imp7_0_balance_phi_and_c2():
    p_gl = _gl_power(lambda x: _q_from_field(_pec_exterior_numpy(), x), A, B)
    c = _coeffs(K_CORE)
    p_cf = c["power"]
    p_rb = 4.0 * np.pi * B**2 * H * float(pec_sphere_implant_radial(B, c, 0))
    pairs = {
        "GL vs closed form": abs(p_gl - p_cf) / p_cf,
        "GL vs 4 pi b^2 h theta0(b)": abs(p_gl - p_rb) / p_cf,
        "closed form vs Robin": abs(p_cf - p_rb) / p_cf,
    }
    for k, v in pairs.items():
        _p(f"[IMP7-0 Bal] {k}: {v:.3e}")
    _p(
        f"[IMP7-0 Bal] P(beta=1) GL {p_gl:.9e} W (predicted 6.803512e-2); band B1_BAND {B1_BAND:.0e}"
    )
    c_out = SHELL_OUTER_FACTOR * A
    vol = 4.0 * np.pi * (c_out**3 - A**3) / 3.0
    phi = (
        _gl_power(lambda x: pec_sphere_implant_source(x, radius=A, q_far=1.0, beta=1.0), A, c_out)
        / vol
    )
    e_phi = abs(phi - PHI_CLOSED_FORM_PEC)
    _p(
        f"[IMP7-0 Phi] shell mean q/q_far over a<r<1.5a {phi:.15f} vs PHI_CLOSED_FORM_PEC "
        f"{PHI_CLOSED_FORM_PEC:.15f}: {e_phi:.3e} (predicted ~1e-15, band {PHI_BAND:.0e})"
    )
    p0 = _gl_power(lambda x: pec_sphere_implant_source(x, radius=A, q_far=Q_FAR, beta=0.0), A, B)
    c2 = 1.0 - p0 / p_gl
    _p(
        f"[IMP7-0 c2] beta=0 source vs (Bal) power: 1 - P(0)/P(1) = {c2:.6f} (P(0) {p0:.6e}, "
        f"predicted 6.597345e-2; predicted 0.030303; asserted > CONTROL_MISS {CONTROL_MISS:.0e})"
    )
    assert max(pairs.values()) <= B1_BAND
    assert e_phi <= PHI_BAND
    assert c2 > CONTROL_MISS


def test_imp7_0_r0_beta0_reduces_to_mat18():
    r = np.linspace(0.0, B, 81)
    rng = np.random.default_rng(57900)
    mu = rng.uniform(-1.0, 1.0, r.size)
    s = np.sqrt(1.0 - mu * mu)
    pts = np.column_stack([r * s, np.zeros_like(r), r * mu])
    ref = (
        two_layer_sphere_steady_temperature(
            r, a=A, b=B, k1=K_CORE, k2=K_GEL, q1=0.0, q2=Q_FAR, h=H, t_inf=T_INF
        )
        - T_INF
    )
    worst = 0.0
    for name, k1 in BRANCHES.items():
        th = (
            pec_sphere_implant_steady_temperature(pts, k_core=k1, t_inf=T_INF, beta=0.0, **FIX)
            - T_INF
        )
        m = r >= A if k1 == 0.0 else np.ones_like(r, dtype=bool)
        e = float(np.max(np.abs(th[m] - ref[m])) / np.max(np.abs(ref[m])))
        worst = max(worst, e)
        _p(f"[IMP7-0 R0] beta=0 {name:10s} vs MAT-18 two_layer_sphere_steady_temperature: {e:.3e}")
    _p(
        f"[IMP7-0 R0] worst {worst:.3e} (band {R0_BAND:.0e}); theta(a)|beta=0 = "
        f"{float(ref[np.argmin(np.abs(r - A))]):.6f} K (predicted 0.318750)"
    )
    assert worst <= R0_BAND


def test_imp7_0_rk_and_predictions_printed():
    """(Rk) and the step's readouts: printed beside their predictions, never asserted."""
    ad, iso = _coeffs(0.0), _coeffs(np.inf)

    def e_to(ref, k1):
        c = _coeffs(k1)
        return abs(float(pec_sphere_implant_radial(A, c, 2) - pec_sphere_implant_radial(A, ref, 2)))

    for eps in (1e-2, 1e-3, 1e-4):
        ea, ea2 = e_to(ad, eps * K_GEL), e_to(ad, 0.5 * eps * K_GEL)
        ei, ei2 = e_to(iso, K_GEL / eps), e_to(iso, 2.0 * K_GEL / eps)
        _p(
            f"[IMP7-0 Rk] eps {eps:.0e}: |dR2(a)| toward adiabatic {ea:.3e}/{ea2:.3e} ratio "
            f"{ea / ea2:.4f}; toward isothermal {ei:.3e}/{ei2:.3e} ratio {ei / ei2:.4f} "
            f"(predicted 2; [1.9, 2.1] a reading aid, not asserted)"
        )
    c20 = _coeffs(K_CORE)
    _p(
        f"[IMP7-0 rec] theta0(a) {c20['theta0_core']:.6f} K (pred 0.350684); theta0(b) "
        f"{c20['theta0_b']:.6f} K (pred 0.135352); P {c20['power']:.6e} W (pred 6.803512e-2)"
    )
    preds = {
        "adiabatic": (0.045660, 0.396344, 0.327854),
        "k1=20": (0.001966, 0.352650, 0.349700),
        "isothermal": (0.0, 0.350684, 0.350684),
    }
    for name, k1 in BRANCHES.items():
        c = _coeffs(k1)
        r2a = float(pec_sphere_implant_radial(A, c, 2))
        pole = (
            float(
                pec_sphere_implant_steady_temperature(
                    np.array([[0, 0, A]]), k_core=k1, t_inf=T_INF, **FIX
                )[0]
            )
            - T_INF
        )
        eq = (
            float(
                pec_sphere_implant_steady_temperature(
                    np.array([[A, 0, 0]]), k_core=k1, t_inf=T_INF, **FIX
                )[0]
            )
            - T_INF
        )
        pr = preds[name]
        _p(
            f"[IMP7-0 rec] {name:10s} R2(a) {r2a:.6f} (pred {pr[0]:.6f}); pole/equator at r=a "
            f"{pole:.6f}/{eq:.6f} (pred {pr[1]:.6f}/{pr[2]:.6f}), ratio {pole / eq:.4f}"
        )
    t_b0 = (
        float(
            two_layer_sphere_steady_temperature(
                A, a=A, b=B, k1=K_CORE, k2=K_GEL, q1=0.0, q2=Q_FAR, h=H, t_inf=T_INF
            )
        )
        - T_INF
    )
    _p(
        f"[IMP7-0 rec] heating factor theta0(a)|b=1 / theta(a)|b=0 = {c20['theta0_core']:.6f}/"
        f"{t_b0:.6f} = {c20['theta0_core'] / t_b0:.4f} (pred 0.350684/0.318750 = 1.1002)"
    )
    cv = _coeffs(0.0, beta=-0.5)
    pv = (
        float(
            pec_sphere_implant_steady_temperature(
                np.array([[0, 0, A]]), k_core=0.0, t_inf=T_INF, beta=-0.5, **FIX
            )[0]
        )
        - T_INF
    )
    ev = (
        float(
            pec_sphere_implant_steady_temperature(
                np.array([[A, 0, 0]]), k_core=0.0, t_inf=T_INF, beta=-0.5, **FIX
            )[0]
        )
        - T_INF
    )
    _p(
        f"[IMP7-0 rec] void beta=-1/2 adiabatic pole/equator {pv:.6f}/{ev:.6f} (pred 0.310844/0.334678); "
        f"R2(a) {float(pec_sphere_implant_radial(A, cv, 2)):.6f}"
    )
