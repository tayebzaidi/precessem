r"""`MAT-21` step 51 — the box volume mean ⟨θ⟩(t) of the six-face Robin box under two sources,
closed form (§9 item 991, re-posed as item 999, §10 P91; briefs `docs/planning/chunks/MAT-21.md`
"§9 item 991" and "§9 item 999").

Numpy / scipy; complex build only because P89′ imports C33 (``mpi4py``) and C6.  Sets no band
and moves none (``B1_BAND`` imported).  No landed module is edited; C6 is import-only.

The box |a| ≤ A, |b| ≤ B, |c| ≤ C (step 0's ``HALF``), ρc θ_t = k∇²θ + q, −k∂_nθ = hθ on all
six faces, θ(0) = 0.  Per axis, even modes cos(γ s/L), γ tan γ = Bi = hL/k
(``robin_slab_eigenvalues``), 𝒩 = ∫_{−L}^{L} cos² = L(γ + sin γ cos γ)/γ; odd modes sin(γ s/L).
Λ = Σ_d (γ_d/L_d)², V = 8ABC.

The box mean of one axis factor (its own antiderivative):
    (1/2L) ∫_{−L}^{L} cos(γ s/L) ds = sin γ/γ,        (1/2L) ∫ sin(γ s/L) ds = 0,
so **every mode with any odd factor has box mean exactly 0 and is dropped** (printed).  With
c_d = 2 sin γ/(γ + sin γ cos γ) (step 0's ``_axis_data``) one has (2L sin γ/γ)/𝒩 = c_d.

(a) Uniform q = P_dep/V.  θ(x, t) = (q/k) Σ_lmn Π_d c_d cos(γ_d x_d/L_d)·(−expm1(−αΛt))/Λ (step 1,
``robin_box_transient_temperature``), so
    ⟨θ⟩_u(t) = (q/k) Σ_lmn Π_d [c_d sin γ_d/γ_d] · (−expm1(−αΛt))/Λ.

(b) P89′'s lattice (1.5 T, entry 52's weights ``p89._data()["w"]["1.5T"]``), ball i of radius
r = ``c6._a10()`` about p_i with q_i = w_i/V_B.  P89 gives ∫_{B_i} Φ = V_B F(√Λ r) Φ(p_i), so
    θ_i(x, t) = (w_i/k) Σ_modes F(√Λ r) Φ(p_i) Φ(x)/𝒩 · (−expm1(−αΛt))/Λ,
and the box mean of Φ(x)/𝒩 = Π_d (sin γ_d/γ_d)/𝒩_d = Π_d c_d/(2L_d) = Π_d c_d / V:
    ⟨θ_i⟩(t) = (w_i/(kV)) Σ_even F(√Λ r) cos(γ_a p_a/L_a) cos(γ_b p_b/L_b) cos(γ_c p_c/L_c)
               · Π_d c_d · (−expm1(−αΛt))/Λ,
with p_c = 0.  (Executor's note: the brief's displayed line writes the per-axis factor as
c_d·sin γ_d/γ_d without the 1/V; the derivation above — the brief's own identity
(2L sin γ/γ)/𝒩 = c_d — gives Π c_d / V for a point-like ball, the sin γ/γ arising only after a
further mean over the source position, which recovers (a).  This module uses the derived form;
(Bq-split) and (N-h) test it.)  The contraction is built here: never ``p89._gmat``, which
weights the c-axis by 1/𝒩_c (the point value at c = 0), not c_c.

(Bq-split): the same mean through P89's split at s₀ = 100 s: the free-space piece integrates over
the box to s₀ Σw/(ρcV) (exact to the image factor e^{−69.7}); the modal piece carries
e^{−αΛs₀}·(−expm1(−αΛ(t − s₀)))/Λ in place of (−expm1(−αΛt))/Λ.

Attempt 2 (§9 item 999, the re-pose of item 991).  Diagnosis: at an even root γ sin γ = Bi cos γ
gives sin γ_n ≈ ±Bi/γ_n, so c_n ≈ 2Bi/γ² and sin γ/γ ≈ Bi/γ² (not ≈ 1/n as item 991 read), and
a_n = c_n sin γ_n/γ_n ∝ γ⁻⁴ > 0.  Every mode with αΛt* ≫ 1 carries the saturated factor 1/Λ, so
the direct [0, t*] series' tail is the steady series' algebraic tail, ∝ N⁻⁵ (rigorous bound
T(N), ``tail_bound``, printed).  The cure is a split in time at s₀, not a steady split (the 3D
Robin-box steady mean *is* the algebraic series).

(a) by the s₀ split rests on an exact product identity: the box survival under unit initial
temperature is Π_d S_d(x_d, τ), its box mean Π_d M_d(τ) with M_d = Σ_n a_{d,n} e^{−α k_n² τ};
Duhamel gives ⟨θ⟩_u(t) = (q/ρc)∫_0^t Π_d M_d dτ, so with q = 1/V
    ⟨θ⟩_u(t) = I_s/(ρcV) + (1/(kV)) Σ_{N₀} Π_d a_d e^{−αΛs₀}(−expm1(−αΛ(t − s₀)))/Λ,
    I_s = ∫_0^{s₀} Π_d M_d^sh dτ,  M_d^sh(τ) = 1 − [erfcx(x) − 1 + 2x/√π]/Bi_d,  x = (h/k)√(ατ),
each face a semi-infinite Robin solid (surface temperature erfcx(x), heat lost ∫hθ_s dτ), image
correction ≈ e^{−L_d²/(αs₀)} ≤ e^{−141}.  I_s is GL64 in u = √τ on [0, √s₀] (x ≤ 0.0632, analytic
integrand).  The split needs P89's N₀ = (180, 120, 30) (αk_N² s₀ = 43 / 46 / 59), not step 0's
(40, 40, 40) (long-axis exponent 2.04).  (S-a) is the seam: M_d^sh(s₀) against the eigen-series
at s₀, per axis.  (b) by P89's split as above; the lattice line is Π c_d / V per ball (item 999,
R1), as derived above.  (Bq-split) compares P89's split at N₀ with the direct mean at 2N₀, N sized
from W1's measured tail (4.07e-13; R3).  (N-h) stays on the direct route (m_short cancels to
noise at Bi ≈ 1e-12) on N_h = (40, 40, 40) for both sources, below the first brentq bracket flip
of ``robin_slab_eigenvalues`` at h = 1e-11 (i = 167 / 84 / 53; R2 (A), no ``src/`` edit).

Asserted (collected, then asserted): preconditions; (Lim-a) the termwise antiderivative vs
64-node Gauss–Legendre (composite panels, printed) per axis; (Lim-∞) on the direct form; (S-a);
(Bq-N) split(N₀) → split(2N₀) for both sources; (Bq-split); (Point) custody vs step 0's
``robin_box_transient_temperature`` at the centre; Π > 0 for the uniform source; negative
control (N-h) at h = 1e-11.  Printed (🟡 record, from the split routes): ⟨θ⟩(t*) per deposited
watt for both sources, the ratio, the stored and surface fractions, I_s.

Not claimed: entry 69's FEM ⟨θ⟩ or nearness to it; a 3 T figure; a degree-2 lattice; a wall
model; an external comparison.
"""

from __future__ import annotations

import math

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import erfcx

from fem_em_solver.utils.analytical import (
    robin_box_transient_temperature,
    robin_slab_eigenvalues,
)
from tests.validation import test_geo38_f_clinical_box_loaded as c6
from tests.validation import test_mat21_lattice_ball_superposition_closed_form as p89
from tests.validation.test_imp4_pec_cylinder_closed_form import B1_BAND
from tests.validation.test_mat21_box_closed_form import (
    A_BOX,
    B_BOX,
    C_BOX,
    HALF,
    RHO_C,
    SERIES_TERMS,
    T_C,
    H,
    K,
)
from tests.validation.test_mat21_box_dose_closed_form import T_GRID

V_BOX = 8.0 * A_BOX * B_BOX * C_BOX
EPS = float(np.finfo(float).eps)
N_U = tuple(SERIES_TERMS)  # (40, 40, 40) even, source (a)
N_L = tuple(p89.N0)  # (180, 120, 30) even, source (b)
H_SMALL = 1e-11
N_HCTRL = 1e-12  # (N-h) ceiling, the weekly's number (= B1_BAND)
# (N-h) truncation (item 999, ruling R2 (A)): below the first brentq bracket flip of
# robin_slab_eigenvalues at h = 1e-11 (i = 167 / 84 / 53 on L = 0.325 / 0.21 / 0.045 m).
N_H = N_U
FLIP_H_SMALL = (167, 84, 53)
N_T_NEEDED = (402, 284, 84)  # drafter's N for T(N) <= B1_BAND (item 999), not run
GL_NODES = 64
SIDE = "1.5T"
# Predictions (brief arithmetic), printed only.
PRED_U = (0.0248, 0.0240, 0.0255)
PRED_L = (0.0273, 0.0268, 0.027473)
PRED_R = (1.10, 1.06, 1.14)


def _say(msg):
    print(msg, flush=True)


def _fail(gates):
    return [g for g in gates if not g[1]]


def _in(v, lo, hi):
    return lo <= v <= hi


def _even(L, n, h=H):
    """(γ, k = γ/L, c, S = sin γ/γ) for the even family of one axis."""
    g = robin_slab_eigenvalues(h * L / K, n)
    c = 2.0 * np.sin(g) / (g + np.sin(g) * np.cos(g))
    return g, g / L, c, np.sin(g) / g


def _lam(kk):
    return (kk[0] ** 2)[:, None, None] + (kk[1] ** 2)[None, :, None] + (kk[2] ** 2)[None, None, :]


def _tfac(lam, t, alpha, s0=0.0, steady=False):
    if steady:
        return np.exp(-alpha * lam * s0) / lam
    out = -np.expm1(-alpha * lam * (t - s0)) / lam
    if s0 > 0.0:
        out *= np.exp(-alpha * lam * s0)
    return out


def mean_uniform(t, n=N_U, h=H, steady=False):
    """⟨θ⟩_u per deposited watt [K/W] (q = 1/V)."""
    ax = [_even(L, m, h) for L, m in zip(HALF, n)]
    alpha = K / RHO_C
    lam = _lam([a[1] for a in ax])
    pf = [a[2] * a[3] for a in ax]
    P = pf[0][:, None, None] * pf[1][None, :, None] * pf[2][None, None, :]
    val = float(np.sum(P * _tfac(lam, t, alpha, steady=steady))) / (K * V_BOX)
    return val, float(P.min())


def point_uniform(t, n=N_U):
    """Series (a) at the box centre through this module's sum (q = 1/V)."""
    ax = [_even(L, m) for L, m in zip(HALF, n)]
    lam = _lam([a[1] for a in ax])
    P = ax[0][2][:, None, None] * ax[1][2][None, :, None] * ax[2][2][None, None, :]
    return float(np.sum(P * _tfac(lam, t, K / RHO_C))) / (K * V_BOX)


def m_short(L, tau, h=H):
    """Short-time slab mean M_d^sh(τ) (unit initial temperature, Robin faces, each face a
    semi-infinite Robin solid): 1 − [erfcx(x) − 1 + 2x/√π]/Bi, x = (h/k)√(ατ).  Image
    correction ≈ e^{−L²/(ατ)}.  Never at h ≈ 1e-11 (the bracket cancels to noise)."""
    alpha = K / RHO_C
    bi = h * L / K
    x = (h / K) * np.sqrt(alpha * np.asarray(tau, dtype=float))
    return 1.0 - (erfcx(x) - 1.0 + 2.0 * x / math.sqrt(math.pi)) / bi


def i_short(s0, h=H, nodes=GL_NODES):
    """I_s = ∫_0^{s0} Π_d M_d^sh(τ) dτ by GL64 in u = √τ (dτ = 2u du)."""
    x, wq = leggauss(nodes)
    rs = math.sqrt(s0)
    u = 0.5 * rs * (x + 1.0)
    prod = np.ones_like(u)
    for L in HALF:
        prod = prod * m_short(L, u * u, h)
    return float(0.5 * rs * np.sum(wq * prod * 2.0 * u))


def mean_uniform_split(t, n=N_L, s0=None):
    """⟨θ⟩_u per deposited watt [K/W] through the s₀ split (item 999):
    I_s/(ρcV) + (1/(kV)) Σ Π a_d e^{−αΛs₀}(−expm1(−αΛ(t − s₀)))/Λ."""
    s0 = p89.S0 if s0 is None else s0
    ax = [_even(L, m) for L, m in zip(HALF, n)]
    alpha = K / RHO_C
    lam = _lam([a[1] for a in ax])
    pf = [a[2] * a[3] for a in ax]
    P = pf[0][:, None, None] * pf[1][None, :, None] * pf[2][None, None, :]
    mod = float(np.sum(P * _tfac(lam, t, alpha, s0=s0))) / (K * V_BOX)
    isv = i_short(s0)
    return isv / (RHO_C * V_BOX) + mod, isv, mod


def lattice_split(w, pts, t, n=N_L, s0=None):
    """P89's route: free-space s₀Σw/(ρcV) + the modal [s₀, t] box mean."""
    s0 = p89.S0 if s0 is None else s0
    free = s0 * float(np.sum(w)) / (RHO_C * V_BOX)
    mod, _ = mean_lattice(w, pts, t, n=n, s0=s0)
    return free + mod, free, mod


def tail_bound(n):
    """T(N): rigorous upper bound on direct(∞) − direct(N) for (a), per deposited watt
    (all a_n > 0, Σ a = 1, 1 − e^{−x} ≤ 1, 1/Λ ≤ L²/γ², γ_n ≥ (n − 1)π)."""
    tot = 0.0
    for L, m in zip(HALF, n):
        bi = H * L / K
        g = (np.arange(m + 1, m + 2_000_001, dtype=float) - 1.0) * math.pi
        tot += float(np.sum(2.0 * bi**2 * L**2 / (g**5 * (g - 0.5))))
    return tot / (K * V_BOX)


def mean_lattice(w, pts, t, n=N_L, h=H, s0=0.0, steady=False):
    """⟨θ⟩_lat per deposited watt [K/W] — the box-mean contraction (even modes only)."""
    r = float(c6._a10())
    ax = [_even(L, m, h) for L, m in zip(HALF, n)]
    alpha = K / RHO_C
    lam = _lam([a[1] for a in ax])
    T = p89._F(np.sqrt(lam) * r) * _tfac(lam, t, alpha, s0=s0, steady=steady)
    del lam
    T *= ax[2][2][None, None, :]
    Tc = T @ np.cos(np.outer(ax[2][1], pts[:, 2]))  # (l, m, i)
    del T
    pa = np.cos(np.outer(pts[:, 0], ax[0][1])) * ax[0][2][None, :]  # (i, l)
    pb = np.cos(np.outer(pts[:, 1], ax[1][1])) * ax[1][2][None, :]  # (i, m)
    per = np.einsum("il,lmi,im->i", pa, Tc, pb) / (K * V_BOX)
    return float(np.dot(w, per)), per


def _lat():
    d = p89._data()
    return d, np.asarray(d["w"][SIDE], dtype=float), np.asarray(d["pts"], dtype=float)


def _gl_mean(kk, L, nodes=GL_NODES):
    """(1/2L)∫_{−L}^{L} cos(k s) ds by composite 64-node Gauss–Legendre."""
    x, wq = leggauss(nodes)
    panels = int(math.ceil(float(np.max(kk)) * 2.0 * L / 40.0)) + 1
    edges = np.linspace(-L, L, panels + 1)
    tot = np.zeros(len(kk))
    for a, b in zip(edges[:-1], edges[1:]):
        s = 0.5 * (b - a) * x + 0.5 * (a + b)
        tot += 0.5 * (b - a) * (np.cos(np.outer(kk, s)) @ wq)
    return tot / (2.0 * L), panels


# ------------------------------------------------------------------ tests
def test_preconditions():
    gates = []
    d, w, pts = _lat()
    r = d["r"]
    gates.append((("n", len(w)), len(w) == 18, len(w)))
    gates.append((("sum w", w.sum()), abs(w.sum() - 1.0) <= 4 * EPS, w.sum() - 1.0))
    worst = max(max(abs(p[k]) + r - HALF[k] for k in range(3)) for p in pts)
    gates.append((("balls inside", worst), worst <= 0.0, worst))
    ts = float(T_GRID[-1])
    gates.append((("t* >= S0", ts), ts >= p89.S0, ts))
    # custody: this module's even axis against p89._axis's even family at step 0's H
    cust = 0.0
    for L, m in zip(HALF, N_L):
        kk, nn, od = p89._axis(L, m, odd_too=False)
        g, k2, c, S = _even(L, m)
        cust = max(
            cust,
            float(np.max(np.abs(kk - k2) / k2)),
            float(np.max(np.abs(2 * L * S / nn - c) / np.abs(c))),
        )
        gates.append((("even only", L), not od.any(), 0))
    gates.append((("axis custody", cust), cust <= B1_BAND, cust))
    _say(
        f"MAT21-51 lattice {len(w)} balls ({SIDE}, entry 52 weights), |sum w - 1| "
        f"{abs(w.sum() - 1.0):.3e} (<= 4 eps ASSERTED); r {r:.6f} m, worst max_d(|p_d| + r - L_d) "
        f"{worst:.6f} m (<= 0 ASSERTED); t* {ts:.1f} s >= S0 {p89.S0:g} s (ASSERTED)"
    )
    _say(
        f"MAT21-51 odd families dropped (box mean of sin(gamma s/L) is exactly 0); even axis vs "
        f"p89._axis and (2L sin g/g)/N = c_d: max rel {cust:.3e} (band {B1_BAND:g}, ASSERTED)"
    )
    _say(f"MAT21-51 V_sub = 8ABC {V_BOX:.6f} m^3, 1/(rho c V) {1.0 / (RHO_C * V_BOX):.4e} K/s/W")
    bad = _fail(gates)
    assert not bad, bad


def test_lim_a_antiderivative():
    gates = []
    worst = 0.0
    for tag, n in (("(a)", N_U), ("(b)", N_L)):
        for L, m in zip(HALF, n):
            g, kk, c, S = _even(L, m)
            gl, panels = _gl_mean(kk, L)
            res = float(np.max(np.abs(gl - S)))
            worst = max(worst, res)
            gates.append((("Lim-a", tag, L, m), res <= B1_BAND, res))
            _say(
                f"MAT21-51 (Lim-a) {tag} axis L {L:g} m, {m} even roots (gamma_max {g[-1]:.2f}): "
                f"max |GL{GL_NODES} x {panels} panels - sin g/g| {res:.3e} (band {B1_BAND:g}, "
                f"ASSERTED)"
            )
    bad = _fail(gates)
    assert not bad, bad


def test_lim_inf_bq_point():
    gates = []
    d, w, pts = _lat()
    ts = float(T_GRID[-1])
    tl = 1000.0 * T_C
    # (Lim-inf)
    u_tr, _ = mean_uniform(tl)
    u_ss, _ = mean_uniform(None, steady=True)
    l_tr, _ = mean_lattice(w, pts, tl)
    l_ss, _ = mean_lattice(w, pts, None, steady=True)
    for tag, a, b in (("(a)", u_tr, u_ss), ("(b)", l_tr, l_ss)):
        rel = abs(a - b) / abs(b)
        gates.append((("Lim-inf", tag, a, b), rel <= B1_BAND, rel))
        _say(
            f"MAT21-51 (Lim-inf) {tag} <theta>(1000 T_C) {a:.9e} vs steady {b:.9e}: rel "
            f"{rel:.3e} (band {B1_BAND:g}, ASSERTED)"
        )
    # (S-a) the seam: short-time slab mean vs the eigen-series at s0, per axis (item 999)
    alpha = K / RHO_C
    for L, m in zip(HALF, N_L):
        g, kk, c, S = _even(L, m)
        ser = float(np.sum(c * S * np.exp(-alpha * kk**2 * p89.S0)))
        msh = float(m_short(L, p89.S0))
        res = abs(msh - ser)
        gates.append((("S-a", L, msh, ser), res <= B1_BAND, res))
        _say(
            f"MAT21-51 (S-a) axis L {L:g} m: M_sh(s0) {msh:.15f} vs sum_{{n<={m}}} a_n "
            f"exp(-alpha k_n^2 s0) {ser:.15f}: abs {res:.3e} (band {B1_BAND:g}, ASSERTED; "
            f"alpha k_N^2 s0 {alpha * kk[-1] ** 2 * p89.S0:.1f})"
        )
    # (Bq-N) on the s0-split routes (item 999)
    n2 = tuple(2 * v for v in N_L)
    u0, isv, umod = mean_uniform_split(ts)
    u1, _, _ = mean_uniform_split(ts, n=n2)
    l0, free, lmod = lattice_split(w, pts, ts)
    l1, _, _ = lattice_split(w, pts, ts, n=n2)
    _say(
        f"MAT21-51 (a) split s0 {p89.S0:g} s: I_s {isv:.10f} s (I_s/s0 {isv / p89.S0:.9f}), "
        f"I_s/(rho c V) {isv / (RHO_C * V_BOX):.12e} + modal [s0, t*] {umod:.12e} K/W"
    )
    for tag, a, b in (("(a)", u0, u1), ("(b)", l0, l1)):
        rel = abs(b - a) / abs(a)
        gates.append((("Bq-N", tag, a, b), rel <= B1_BAND, rel))
        _say(
            f"MAT21-51 (Bq-N) {tag} split N {N_L} -> 2N {n2}: {a:.15e} -> {b:.15e} K/W, rel "
            f"{rel:.3e} (band {B1_BAND:g}, ASSERTED)"
        )
    _, pmin = mean_uniform(ts)
    gates.append((("Pi > 0", pmin), pmin > 0.0, pmin))
    _say(f"MAT21-51 (a) min over modes of prod c_d sin g_d/g_d {pmin:.3e} (> 0 ASSERTED)")
    # (a) direct-route tail, printed only: gap to split(N0) beside the rigorous bound T(N)
    for n in (N_U, tuple(2 * v for v in N_U)):
        dv, _ = mean_uniform(ts, n=n)
        gap = (u0 - dv) / u0
        tb = tail_bound(n) / u0
        _say(
            f"MAT21-51 (a) direct N {n}: {dv:.15e}, gap to split(N0) rel {gap:.3e} beside T(N) "
            f"rel {tb:.3e}, gap <= T {gap <= tb} (printed, not asserted)"
        )
    _say(
        f"MAT21-51 (a) N for T(N) <= B1_BAND on the direct route {N_T_NEEDED} (drafter, item "
        f"999; not run)"
    )
    # (Bq-split) P89's split at N0 vs the direct mean at 2N0 (item 999, R3)
    ld, _ = mean_lattice(w, pts, ts, n=n2)
    rel = abs(l0 - ld) / abs(ld)
    gates.append((("Bq-split", l0, ld), rel <= B1_BAND, rel))
    _say(
        f"MAT21-51 (Bq-split) s0 {p89.S0:g} s: free {free:.12e} + modal [s0, t*] {lmod:.12e} = "
        f"{l0:.16e} (N0) vs direct {ld:.16e} (2N0 {n2}): rel {rel:.3e} (band {B1_BAND:g}, "
        f"ASSERTED)"
    )
    # (Point)
    me = point_uniform(ts)
    ref = float(
        robin_box_transient_temperature(
            0.0,
            0.0,
            0.0,
            ts,
            half_widths=HALF,
            k=K,
            q=1.0 / V_BOX,
            h=H,
            rho_c=RHO_C,
            t_inf=0.0,
            n_terms=N_U,
        )
    )
    rel = abs(me - ref) / abs(ref)
    gates.append((("Point", me, ref), rel <= B1_BAND, rel))
    _say(
        f"MAT21-51 (Point) centre theta(t*) at q = 1/V: module {me:.9e} vs "
        f"robin_box_transient_temperature {ref:.9e}: rel {rel:.3e} (band {B1_BAND:g}, ASSERTED)"
    )
    bad = _fail(gates)
    assert not bad, bad


def test_negative_control_small_h():
    gates = []
    d, w, pts = _lat()
    ts = float(T_GRID[-1])
    stored = ts / (RHO_C * V_BOX)
    for L in HALF:
        bi = H_SMALL * L / K
        g0 = float(robin_slab_eigenvalues(bi, 1)[0])
        res = abs(g0 * math.sin(g0) - bi * math.cos(g0)) / bi
        _say(
            f"MAT21-51 (N-h) Bi {bi:.3e}: gamma_0 {g0:.6e} (sqrt Bi {math.sqrt(bi):.6e}), "
            f"gamma_0^2/Bi - 1 {g0 * g0 / bi - 1.0:.3e} (-Bi/3 {-bi / 3.0:.3e}), "
            f"|g sin g - Bi cos g|/Bi {res:.3e} (printed)"
        )
    _say(
        f"MAT21-51 (N-h) N_h {N_H} below first bracket flip i = "
        f"{' / '.join(str(v) for v in FLIP_H_SMALL)} (h = {H_SMALL:g})"
    )
    u, _ = mean_uniform(ts, n=N_H, h=H_SMALL)
    lt, _ = mean_lattice(w, pts, ts, n=N_H, h=H_SMALL)
    for tag, v in (("(a)", u), ("(b)", lt)):
        rel = abs(v / stored - 1.0)
        gates.append((("N-h", tag, v, stored), rel <= N_HCTRL, rel))
        _say(
            f"MAT21-51 (N-h) {tag} h {H_SMALL:g}: <theta>(t*) {v:.9e} vs P_dep t*/(rho c V) "
            f"{stored:.9e}: |ratio - 1| {rel:.3e} (ceiling {N_HCTRL:g}, ASSERTED; predicted "
            f"approx 1e-13)"
        )
    bad = _fail(gates)
    assert not bad, bad


def test_record_printed():
    d, w, pts = _lat()
    ts = float(T_GRID[-1])
    stored = ts / (RHO_C * V_BOX)
    u, isv, _ = mean_uniform_split(ts)
    lt, _, _ = lattice_split(w, pts, ts)
    ratio = lt / u
    _say(
        f"MAT21-51 records from the s0-split routes (N0 {N_L}, s0 {p89.S0:g} s); I_s {isv:.7g} s, "
        f"I_s/s0 {isv / p89.S0:.7g}"
    )
    _say(
        "MAT21-51 normaliser: <theta> per deposited watt (no P_acc factor); entry 69's row 5 is "
        "per deposited watt as well"
    )
    for tag, v, (pv, lo, hi) in (("uniform", u, PRED_U), ("lattice 1.5 T deg 1", lt, PRED_L)):
        _say(
            f"MAT21-51 RECORD <theta>_{tag}(t*) {v:.7g} K/W (predicted approx {pv}, range "
            f"[{lo}, {hi}] inside {_in(v, lo, hi)}); stored fraction {v / stored:.6f}, surface "
            f"fraction {1.0 - v / stored:.6f}"
        )
    pv, lo, hi = PRED_R
    _say(
        f"MAT21-51 RECORD ratio lattice / uniform {ratio:.6f} (predicted approx {pv}, range "
        f"[{lo}, {hi}] inside {_in(ratio, lo, hi)}); beside P_dep t*/(rho c V) {stored:.6f} K/W"
    )
    assert u > 0.0 and lt > 0.0
