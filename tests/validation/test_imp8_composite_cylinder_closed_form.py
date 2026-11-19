"""`IMP-8` step 9 (§9 item 824, re-posed as item 843; §10 I64): the composite infinite cylinder in
closed form, numpy / scipy only.

Brief: docs/planning/chunks/IMP-8.md, "§9 item 824 — 2026-11-05", under the
"Review ruling 2026-11-05 03:00" paragraph above it (fixture values).  A rod
core 0 < r < a (k_r, ρc_r) in perfect contact with a gel shell a < r < b
(K_GEL, RHO_C), θ = 0 at r = b, a uniform source q = Q_SHELL in the shell only
(1 W per metre of rod), switched on at t = 0.  Results are in K per W/m.
Carslaw & Jaeger ch. XIII (composite regions) is cited from memory; the text
is not in the repo.  Template: `MAT-18`'s two-layer sphere
(`src/fem_em_solver/utils/analytical.py`), cylinder functions for sin/r.

Modes, μᵢ = √(s/κᵢ): core J0(μ₁r); shell Z = A J0(μ₂r) + B Y0(μ₂r), with θ and
k∂ᵣθ continuous at a.  A and B from the Wronskian W[J0, Y0](x) = 2/(πx)
(Jin 3e App. E, (E.18)–(E.19)).  Eigencondition F(√s) = Z(b) = 0, scanned.
Nₙ by Lommel's integral, cₙ = q ∫ₐᵇ Z r dr / (sₙNₙ) (the Green-form
projection of θ_ss), θ = θ_ss − Σ cₙψₙ e^{−sₙt}.

Asserted (bands imported, unmoved): (St), (En), (Lim) against the homogeneous
Dirichlet J0 series with the source confined to the shell (written here; its
roots are jn_zeros, independent of the scan), (Bq), the precondition (|F| over
the envelope sqrt(Z0^2 + Z1^2)(b) <= 1e-12, review ruling 2026-11-05 10:30;
the attempt-1 scale |A J0| + |B Y0| is printed as its control), and the
control (ρc_r → RHO_C alone misses at t_c by > INTERIOR_BAND at r = 0 and a,
predicted ≈ 0.151 / 0.140, printed beside).  Printed (🟡 record): θ for three
cores — Ti (nominal pair), gel core, adiabatic hole.  Scope: infinite
cylinder, Dirichlet wall, uniform shell source, constant properties, perfect
contact; no finite rod, no SAR shape, no Ti property claim.
"""

from __future__ import annotations

from functools import lru_cache

import numpy as np
from mpi4py import MPI
from scipy import special

from tests.mesh.test_birdcage_implant_rod import ROD_RADIUS
from tests.validation.test_imp4_pec_cylinder_closed_form import B1_BAND
from tests.validation.test_mat12_transient_sar_heating_sphere import RHO_C
from tests.validation.test_mat13_transient_perfusion_sphere import INTERIOR_BAND
from tests.validation.test_mat15_coil_driven_heating_phantom import K_GEL

# Fixture (review ruling 2026-11-05 03:00; parameters, not bands).
A_ROD = ROD_RADIUS  # 1 mm, test_birdcage_implant_rod.py:93
B_WALL = 5.0 * A_ROD
KAPPA_GEL = K_GEL / RHO_C
T_C = A_ROD**2 / KAPPA_GEL  # 6.966667 s
T_STAR = 5.0 * T_C
T_ST = 50.0 * B_WALL**2 / KAPPA_GEL
ROD_K = 6.7  # W/(m·K), Ti-6Al-4V room temperature  # nominal (from memory, no source in repo)
ROD_RHO_C = 2.33018e6  # J/(m³·K), Ti-6Al-4V  # nominal (from memory, no source in repo)
Q_SHELL = 1.0 / (np.pi * (B_WALL**2 - A_ROD**2))  # W/m³ → 1 W/m
N_MODES = 400
EPS = np.finfo(float).eps
R_PTS = (0.0, A_ROD, 0.5 * (A_ROD + B_WALL))
PRED_CTRL = (0.151, 0.140)  # unlogged review arithmetic, predicted (rule (e))


def _p(*args):
    if MPI.COMM_WORLD.rank == 0:
        print(*args, flush=True)


def _theta_ss(r, q=Q_SHELL, a=A_ROD, b=B_WALL, k_g=K_GEL):
    r = np.maximum(np.asarray(r, dtype=float), a)
    return q * (b**2 - r**2) / (4 * k_g) - q * a**2 * np.log(b / r) / (2 * k_g)


def _theta_ss_d1(r, q=Q_SHELL, a=A_ROD, k_g=K_GEL):
    return -q * r / (2 * k_g) + q * a**2 / (2 * k_g * r)


def _theta_ss_d2(r, q=Q_SHELL, a=A_ROD, k_g=K_GEL):
    return -q / (2 * k_g) - q * a**2 / (2 * k_g * r**2)


def _coeffs(sig, k_r, rhoc_r, hole=False, a=A_ROD):
    """A, B of the shell function for √s = sig (arrays)."""
    kap_r, kap_g = k_r / rhoc_r, KAPPA_GEL
    mu1, mu2 = sig / np.sqrt(kap_r), sig / np.sqrt(kap_g)
    x = mu2 * a
    if hole:
        v, dn = np.ones_like(sig), np.zeros_like(sig)
    else:
        v = special.j0(mu1 * a)
        dn = -k_r * mu1 * special.j1(mu1 * a) / (K_GEL * mu2)  # Z'(a)/μ₂
    # W[J0, Y0](x) = 2/(πx); J0' = −J1, Y0' = −Y1.
    A = 0.5 * np.pi * x * (-v * special.y1(x) - dn * special.y0(x))
    B = 0.5 * np.pi * x * (dn * special.j0(x) + v * special.j1(x))
    return A, B, mu1, mu2


def _F(sig, k_r, rhoc_r, hole=False):
    A, B, _, mu2 = _coeffs(sig, k_r, rhoc_r, hole)
    xb = mu2 * B_WALL
    t1, t2 = A * special.j0(xb), B * special.y0(xb)
    return t1 + t2, np.abs(t1) + np.abs(t2)


def _lam(k_r, rhoc_r, hole):
    shell = (B_WALL - A_ROD) / np.sqrt(KAPPA_GEL)
    return shell if hole else A_ROD / np.sqrt(k_r / rhoc_r) + shell


@lru_cache(maxsize=None)
def _modes(k_r, rhoc_r, n, hole=False):
    """The first n roots √sₙ (scan + vectorised bisection) and mode data."""
    lam = _lam(k_r, rhoc_r, hole)
    step = np.pi / (32 * lam)
    lo_list, hi_list, x0 = [], [], step
    while len(lo_list) < n:
        grid = x0 + step * np.arange(0, 32 * (n - len(lo_list)) + 64)
        f, _ = _F(grid, k_r, rhoc_r, hole)
        idx = np.nonzero(np.sign(f[:-1]) * np.sign(f[1:]) < 0)[0]
        lo_list.extend(grid[idx])
        hi_list.extend(grid[idx + 1])
        x0 = grid[-1]
    lo = np.array(lo_list[:n])
    hi = np.array(hi_list[:n])
    flo, _ = _F(lo, k_r, rhoc_r, hole)
    for _ in range(200):
        if np.all(hi - lo <= 4 * EPS * hi):
            break
        mid = 0.5 * (lo + hi)
        fm, _ = _F(mid, k_r, rhoc_r, hole)
        left = np.sign(fm) == np.sign(flo)
        lo = np.where(left, mid, lo)
        flo = np.where(left, fm, flo)
        hi = np.where(left, hi, mid)
    fl, _ = _F(lo, k_r, rhoc_r, hole)
    fh, _ = _F(hi, k_r, rhoc_r, hole)
    sig = np.where(np.abs(fl) <= np.abs(fh), lo, hi)
    fval, fscale = _F(sig, k_r, rhoc_r, hole)
    A, B, mu1, mu2 = _coeffs(sig, k_r, rhoc_r, hole)
    # Precondition scale (review ruling 2026-11-05 10:30): the envelope
    # env = sqrt(Z0(mu2 b)^2 + Z1(mu2 b)^2), same A, B as F; B*Y1 only where B != 0.
    xb_e = mu2 * B_WALL
    z1b_e = A * special.j1(xb_e) + np.where(B != 0, B * special.y1(xb_e), 0.0)
    env = np.sqrt(fval**2 + z1b_e**2)
    s = sig**2
    a, b = A_ROD, B_WALL

    def z0(x):
        return A * special.j0(x) + B * special.y0(x)

    def z1(x):
        return A * special.j1(x) + B * special.y1(x)

    shell_n = RHO_C * (
        0.5 * b**2 * (z0(mu2 * b) ** 2 + z1(mu2 * b) ** 2)
        - 0.5 * a**2 * (z0(mu2 * a) ** 2 + z1(mu2 * a) ** 2)
    )
    core_n = (
        0.0 if hole else rhoc_r * 0.5 * a**2 * (special.j0(mu1 * a) ** 2 + special.j1(mu1 * a) ** 2)
    )
    norm = core_n + shell_n
    int_shell = (b * z1(mu2 * b) - a * z1(mu2 * a)) / mu2  # ∫ₐᵇ Z r dr
    int_core = 0.0 if hole else a * special.j1(mu1 * a) / mu1
    c = Q_SHELL * int_shell / (s * norm)
    return dict(
        sig=sig,
        s=s,
        A=A,
        B=B,
        mu1=mu1,
        mu2=mu2,
        norm=norm,
        c=c,
        fval=fval,
        fscale=fscale,
        env=env,
        lo=lo,
        hi=hi,
        int_shell=int_shell,
        int_core=int_core,
        hole=hole,
        rhoc_r=rhoc_r,
        k_r=k_r,
    )


def _psi(m, r):
    if r < A_ROD:
        if m["hole"]:
            return np.full_like(m["sig"], np.nan)
        return special.j0(m["mu1"] * r)
    x = m["mu2"] * r
    return m["A"] * special.j0(x) + m["B"] * special.y0(x)


def _theta(r, t, *, k_r=ROD_K, rhoc_r=ROD_RHO_C, n_modes=N_MODES, hole=False):
    """θ(r, t) in K per W/m; t = np.inf is the steady state.  Never t = 0."""
    if np.isinf(t):
        return float(_theta_ss(r))
    m = _modes(float(k_r), float(rhoc_r), int(n_modes), bool(hole))
    return float(_theta_ss(r) - np.sum(m["c"] * _psi(m, r) * np.exp(-m["s"] * t)))


def _theta_j0(r, t, n_modes=N_MODES):
    """(Lim)'s comparand: the homogeneous Dirichlet J0 series, gel throughout,
    source q in a < r < b only."""
    lam = special.jn_zeros(0, n_modes) / B_WALL
    s = KAPPA_GEL * lam**2
    norm = RHO_C * 0.5 * B_WALL**2 * special.j1(lam * B_WALL) ** 2
    integ = (B_WALL * special.j1(lam * B_WALL) - A_ROD * special.j1(lam * A_ROD)) / lam
    c = Q_SHELL * integ / (s * norm)
    return float(_theta_ss(r) - np.sum(c * special.j0(lam * r) * np.exp(-s * t)))


def _gram(m, n):
    """Gram matrix ∫ρcψₘψₙ r dr / √(NₘNₙ) by Lommel's cross integral."""
    a, b = A_ROD, B_WALL
    G = np.zeros((n, n))

    def cross(al, be, f0a, f1a, g0b, g1b, r):
        return r * (al * f1a * g0b - be * f0a * g1b) / (al**2 - be**2)

    for i in range(n):
        for j in range(n):
            if i == j:
                G[i, j] = 1.0
                continue
            val = 0.0
            if not m["hole"]:
                al, be = m["mu1"][i], m["mu1"][j]
                val += m["rhoc_r"] * cross(
                    al,
                    be,
                    special.j0(al * a),
                    special.j1(al * a),
                    special.j0(be * a),
                    special.j1(be * a),
                    a,
                )
            al, be = m["mu2"][i], m["mu2"][j]
            Ai, Bi, Aj, Bj = m["A"][i], m["B"][i], m["A"][j], m["B"][j]
            for r, sgn in ((b, 1.0), (a, -1.0)):
                f0 = Ai * special.j0(al * r) + Bi * special.y0(al * r)
                f1 = Ai * special.j1(al * r) + Bi * special.y1(al * r)
                g0 = Aj * special.j0(be * r) + Bj * special.y0(be * r)
                g1 = Aj * special.j1(be * r) + Bj * special.y1(be * r)
                val += sgn * RHO_C * cross(al, be, f0, f1, g0, g1, r)
            G[i, j] = val / np.sqrt(m["norm"][i] * m["norm"][j])
    return G


def test_imp8_composite_cylinder_closed_form():
    gates = {}
    scale = _theta_ss(A_ROD)
    _p(
        f"[I64] a = {A_ROD:.4g} m, b = {B_WALL:.4g} m, t_c = {T_C:.7g} s, t* = {T_STAR:.7g} s, "
        f"q = {Q_SHELL:.7g} W/m3, kappa_gel = {KAPPA_GEL:.7g}, kappa_Ti = {ROD_K / ROD_RHO_C:.7g} (nominal)"
    )
    _p(
        f"[I64] theta_ss(a) = {float(_theta_ss(A_ROD)):.7f} (review 0.1148409), "
        f"theta_ss(3 mm) = {float(_theta_ss(3e-3)):.7f} (review 0.0827736) K per W/m"
    )

    # Precondition.
    pre = {}
    for name, kr, rc, hole in (
        ("Ti", ROD_K, ROD_RHO_C, False),
        ("gel", K_GEL, RHO_C, False),
        ("hole", ROD_K, ROD_RHO_C, True),
        ("ctrl", ROD_K, RHO_C, False),
    ):
        m = _modes(float(kr), float(rc), N_MODES, hole)
        rel = np.abs(m["fval"]) / m["env"]
        rel_old = (
            np.abs(m["fval"]) / m["fscale"]
        )  # attempt-1 scale |A J0| + |B Y0|: printed control only
        pre[name] = (
            float(rel.max()),
            len(m["sig"]),
            float(np.min(np.diff(m["sig"])) / (np.pi / _lam(kr, rc, hole))),
        )
        _p(
            f"[pre] {name}: max |F|/env = {rel.max():.3e} (n = {int(np.argmax(rel)) + 1}), roots = {len(m['sig'])}, "
            f"min spacing / (pi/Lambda) = {pre[name][2]:.4g}, sqrt(s)_1 = {m['sig'][0]:.7g}; "
            f"control (printed): max |F|/(|A J0|+|B Y0|) = {rel_old.max():.3e}"
        )
    gates["pre_F"] = max(v[0] for v in pre.values())
    gates["pre_count_ok"] = all(v[1] == N_MODES for v in pre.values())
    G = _gram(_modes(float(ROD_K), float(ROD_RHO_C), 20, False), 20)
    _p(f"[pre] Ti Gram N=20 max off-diagonal = {np.max(np.abs(G - np.eye(20))):.3e}")
    Gh = _gram(_modes(float(ROD_K), float(ROD_RHO_C), 20, True), 20)
    _p(f"[pre] hole Gram N=20 max off-diagonal = {np.max(np.abs(Gh - np.eye(20))):.3e}")

    # (St).
    st = [abs(_theta(r, T_ST) - float(_theta_ss(r))) / scale for r in R_PTS]
    sc = Q_SHELL * B_WALL**2 / K_GEL
    r3 = np.array([1.5e-3, 3e-3, 4.5e-3])
    res = [
        abs(float(_theta_ss(B_WALL))) / sc,
        abs(K_GEL * _theta_ss_d1(A_ROD)) * B_WALL / (K_GEL * sc),
    ]
    res += list(
        np.abs(K_GEL * (_theta_ss_d2(r3) + _theta_ss_d1(r3) / r3) + Q_SHELL)
        * B_WALL**2
        / (K_GEL * sc)
    )
    gates["St_series"] = max(st)
    gates["St_resid"] = max(res)
    _p(
        f"[St] series at 50 b2/kappa vs theta_ss: {', '.join(f'{v:.3e}' for v in st)}; "
        f"residuals: {', '.join(f'{v:.3e}' for v in res)}"
    )

    # (En).
    m = _modes(float(ROD_K), float(ROD_RHO_C), N_MODES, False)
    lhs = m["s"] * (ROD_RHO_C * m["int_core"] + RHO_C * m["int_shell"])
    zb1 = m["A"] * special.j1(m["mu2"] * B_WALL) + m["B"] * special.y1(m["mu2"] * B_WALL)
    rhs = -K_GEL * B_WALL * (-m["mu2"] * zb1)
    mag = (
        np.abs(m["s"] * ROD_RHO_C * m["int_core"])
        + np.abs(m["s"] * RHO_C * m["int_shell"])
        + np.abs(rhs)
    )
    gates["En_modes"] = float(np.max(np.abs(lhs - rhs) / mag))
    bal = -K_GEL * B_WALL * _theta_ss_d1(B_WALL)
    gates["En_steady"] = abs(bal - Q_SHELL * (B_WALL**2 - A_ROD**2) / 2) / (
        Q_SHELL * (B_WALL**2 - A_ROD**2) / 2
    )
    _p(
        f"[En] modes max rel = {gates['En_modes']:.3e}, steady balance rel = {gates['En_steady']:.3e} "
        f"(q(b2-a2)/2 = {Q_SHELL * (B_WALL**2 - A_ROD**2) / 2:.7g} W/m)"
    )

    # (Lim).
    lim = []
    for t in (T_C, T_STAR):
        for r in R_PTS:
            comp = _theta(r, t, k_r=K_GEL, rhoc_r=RHO_C)
            ref = _theta_j0(r, t)
            lim.append(abs(comp - ref) / scale)
            _p(
                f"[Lim] r = {r * 1e3:.1f} mm, t = {t:.7g} s: composite {comp:.7g}, J0 series {ref:.7g}, "
                f"diff {lim[-1]:.3e}"
            )
    mg = _modes(float(K_GEL), float(RHO_C), N_MODES, False)
    lam_c = mg["mu2"][:50]
    lam_j = special.jn_zeros(0, 50) / B_WALL
    gates["Lim_theta"] = max(lim)
    gates["Lim_roots"] = float(np.max(np.abs(lam_c - lam_j) / lam_j))
    _p(f"[Lim] first 50 roots max rel = {gates['Lim_roots']:.3e}")

    # (Bq).
    bq = [
        abs(_theta(r, t, n_modes=200) - _theta(r, t)) / scale for t in (T_C, T_STAR) for r in R_PTS
    ]
    gates["Bq"] = max(bq)
    _p(f"[Bq] N=200 vs 400: {', '.join(f'{v:.3e}' for v in bq)}")

    # Control: ρc_r → RHO_C alone.
    ctrl = []
    for r in R_PTS:
        ti = _theta(r, T_C)
        cc = _theta(r, T_C, rhoc_r=RHO_C)
        ctrl.append(abs(cc - ti) / abs(ti))
        _p(
            f"[ctrl] r = {r * 1e3:.1f} mm, t_c: Ti {ti:.7g}, rhoc->RHO_C {cc:.7g}, rel miss {ctrl[-1]:.4g}, "
            f"abs/theta_ss(a) {abs(cc - ti) / scale:.4g}"
        )
    gates["ctrl_min"] = min(ctrl[0], ctrl[1])
    _p(
        f"[ctrl] r = 0 / a: measured {ctrl[0]:.4g} / {ctrl[1]:.4g}, predicted {PRED_CTRL[0]} / {PRED_CTRL[1]}; "
        f"(a+b)/2 printed only {ctrl[2]:.4g}"
    )

    # Printed record: three cores.
    cores = (("Ti", dict()), ("gel", dict(k_r=K_GEL, rhoc_r=RHO_C)), ("hole", dict(hole=True)))
    at_a = {}
    for name, kw in cores:
        for t, tl in ((T_C, "t_c"), (T_STAR, "t*"), (np.inf, "steady")):
            vals = []
            for r in R_PTS:
                if kw.get("hole") and r < A_ROD:
                    vals.append("   n/a   ")
                    continue
                v = _theta(r, t, **kw)
                vals.append(f"{v:.7g}")
                if r == A_ROD:
                    at_a[(name, tl)] = v
            _p(f"[rec] {name:4s} {tl:6s} theta(0, a, (a+b)/2) = {', '.join(vals)} K per W/m")
    for tl in ("t_c", "t*"):
        order = sorted(("Ti", "gel", "hole"), key=lambda n: -at_a[(n, tl)])
        _p(f"[rec] ordering at r = a, {tl}: {' > '.join(order)}")
    for n in (400, 800):
        tail = [abs(_theta(r, 0.0, n_modes=n)) / scale for r in R_PTS]
        _p(f"[rec] t = 0 tail, Ti, N = {n}: {', '.join(f'{v:.3e}' for v in tail)} (of theta_ss(a))")

    _p(
        "[gates] "
        + ", ".join(
            f"{k} = {v:.3e}" if isinstance(v, float) else f"{k} = {v}" for k, v in gates.items()
        )
    )
    assert gates["pre_count_ok"], pre
    assert gates["pre_F"] <= 1e-12, pre
    assert gates["St_series"] <= B1_BAND, gates
    assert gates["St_resid"] <= B1_BAND, gates
    assert gates["En_modes"] <= B1_BAND, gates
    assert gates["En_steady"] <= B1_BAND, gates
    assert gates["Lim_theta"] <= B1_BAND, gates
    assert gates["Lim_roots"] <= B1_BAND, gates
    assert gates["Bq"] <= B1_BAND, gates
    assert gates["ctrl_min"] > INTERIOR_BAND, gates
