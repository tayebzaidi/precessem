"""`IMP-8` step 12 (§9 item 847; §10 I70): the composite infinite cylinder with
the source in the core, closed form, numpy / scipy only.

Brief: docs/planning/chunks/IMP-8.md, "§9 item 847 — 2026-11-05", under the
"Review ruling 2026-11-05 18:00" paragraph.  Item 843's fixture and modes
(`test_imp8_composite_cylinder_closed_form`, not edited): a rod core
0 < r < a (k_r, ρc_r) in perfect contact with a gel shell a < r < b
(K_GEL, RHO_C), θ = 0 at r = b, from θ = 0.  The only change is the source:
q = Q_CORE = 1/(πa²) in r < a (1 W per metre of rod), zero in the shell.
Results in K per W/m.

Steady state: core θ_ss = q(a² − r²)/(4k_r) + q a² ln(b/a)/(2k_g); shell
θ_ss = q a² ln(b/r)/(2k_g).  Series: cₙ = q·∫₀ᵃ J0(μ₁r) r dr/(sₙNₙ) (Green-form
projection of θ_ss onto item 843's modes; `_modes` returns `int_core` and the
ρc-weighted `norm`), θ = θ_ss − Σ cₙψₙ e^{−sₙt}.

Asserted (bands imported, unmoved: B1_BAND 1e-12, INTERIOR_BAND 1e-2):
(St), (En), (Lim) against the homogeneous Dirichlet J0 series with the source
confined to r < a (roots jn_zeros, independent of the scan), (Bq), and the
control (k_r → K_GEL alone misses the composite's steady θ(0) by
> INTERIOR_BAND; predicted 0.2752, unlogged review arithmetic, printed beside).
Scope: infinite cylinder, Dirichlet wall, uniform core source, constant
properties, perfect contact; no realistic implant source, no finite rod, no
Ti property claim.
"""

from __future__ import annotations

import numpy as np
from mpi4py import MPI
from scipy import special

from tests.validation.test_imp4_pec_cylinder_closed_form import B1_BAND
from tests.validation.test_imp8_composite_cylinder_closed_form import (
    A_ROD,
    B_WALL,
    KAPPA_GEL,
    N_MODES,
    R_PTS,
    ROD_K,
    ROD_RHO_C,
    T_C,
    T_ST,
    T_STAR,
    _gram,
    _modes,
    _psi,
)
from tests.validation.test_imp8_composite_cylinder_closed_form import (
    _theta_ss as _theta_ss_shell,  # item 843's shell-source steady state, printed comparison only
)
from tests.validation.test_mat12_transient_sar_heating_sphere import RHO_C
from tests.validation.test_mat13_transient_perfusion_sphere import INTERIOR_BAND
from tests.validation.test_mat15_coil_driven_heating_phantom import K_GEL

Q_CORE = 1.0 / (np.pi * A_ROD**2)  # W/m³ → 1 W/m
PRED_SS0 = (0.43879, 0.55955)  # Ti / gel core θ_ss(0), unlogged review arithmetic, predicted
PRED_SSA = 0.42692  # θ_ss(a), unlogged review arithmetic, predicted
PRED_CTRL = 0.2752  # 0.55955 / 0.43879 − 1, unlogged review arithmetic, predicted (rule (e))


def _p(*args):
    if MPI.COMM_WORLD.rank == 0:
        print(*args, flush=True)


def _theta_ss(r, k_r=ROD_K, q=Q_CORE, a=A_ROD, b=B_WALL, k_g=K_GEL):
    r = float(r)
    if r < a:
        return q * (a**2 - r**2) / (4 * k_r) + q * a**2 * np.log(b / a) / (2 * k_g)
    return q * a**2 * np.log(b / r) / (2 * k_g)


def _d1(r, k_r=ROD_K, q=Q_CORE, a=A_ROD, k_g=K_GEL, side=None):
    inside = (r < a) if side is None else (side == "core")
    return -q * r / (2 * k_r) if inside else -q * a**2 / (2 * k_g * r)


def _d2(r, k_r=ROD_K, q=Q_CORE, a=A_ROD, k_g=K_GEL):
    return -q / (2 * k_r) if r < a else q * a**2 / (2 * k_g * r**2)


def _theta(r, t, *, k_r=ROD_K, rhoc_r=ROD_RHO_C, n_modes=N_MODES):
    """θ(r, t) in K per W/m; t = np.inf is the steady state.  Never t = 0."""
    if np.isinf(t):
        return float(_theta_ss(r, k_r=k_r))
    m = _modes(float(k_r), float(rhoc_r), int(n_modes), False)
    c = Q_CORE * m["int_core"] / (m["s"] * m["norm"])
    return float(_theta_ss(r, k_r=k_r) - np.sum(c * _psi(m, r) * np.exp(-m["s"] * t)))


def _theta_j0(r, t, n_modes=N_MODES):
    """(Lim)'s comparand: homogeneous Dirichlet J0 series, gel throughout,
    source Q_CORE in r < a only."""
    lam = special.jn_zeros(0, n_modes) / B_WALL
    s = KAPPA_GEL * lam**2
    norm = RHO_C * 0.5 * B_WALL**2 * special.j1(lam * B_WALL) ** 2
    c = Q_CORE * (A_ROD * special.j1(lam * A_ROD) / lam) / (s * norm)
    return float(_theta_ss(r, k_r=K_GEL) - np.sum(c * special.j0(lam * r) * np.exp(-s * t)))


def test_imp8_composite_cylinder_core_source_closed_form():
    gates = {}
    a, b = A_ROD, B_WALL
    scale = Q_CORE * a**2 / K_GEL
    _p(f"[I70] a = {a:.4g} m, b = {b:.4g} m, q_c = {Q_CORE:.7g} W/m3, t_c = {T_C:.7g} s, t* = {T_STAR:.7g} s, "
       f"T_ST = {T_ST:.7g} s, scale q_c a^2/k_g = {scale:.7g} K per W/m")
    _p(f"[I70] theta_ss(0) Ti = {_theta_ss(0.0):.7g} (predicted {PRED_SS0[0]}), gel core = "
       f"{_theta_ss(0.0, k_r=K_GEL):.7g} (predicted {PRED_SS0[1]}); theta_ss(a) = {_theta_ss(a):.7g} "
       f"(predicted {PRED_SSA})")
    m20 = _modes(float(ROD_K), float(ROD_RHO_C), 20, False)
    G = _gram(m20, 20)
    _p(f"[pre] Ti Gram N=20 max off-diagonal = {np.max(np.abs(G - np.eye(20))):.3e} (item 843's modes, printed)")

    # (St): the series at T_ST against θ_ss, and θ_ss's own residuals.
    st = [abs(_theta(r, T_ST) - _theta_ss(r)) / scale for r in R_PTS]
    eps_r = 1e-12 * a
    res = {
        "theta(b)": abs(_theta_ss(b)) / scale,
        "theta cont at a": abs(_theta_ss(a - eps_r) - _theta_ss(a)) / scale,
        "flux cont at a": abs(ROD_K * _d1(a, side="core") - K_GEL * _d1(a, side="shell")) / (Q_CORE * a / 2),
    }
    for r in (0.3e-3, 0.7e-3):
        res[f"pde core {r * 1e3:.1f} mm"] = abs(ROD_K * (_d2(r) + _d1(r) / r) + Q_CORE) / Q_CORE
    for r in (2e-3, 4e-3):
        res[f"pde shell {r * 1e3:.1f} mm"] = abs(K_GEL * (_d2(r) + _d1(r) / r) + 0.0) / Q_CORE
    gates["St_series"] = max(st)
    gates["St_resid"] = max(res.values())
    _p(f"[St] series at T_ST vs theta_ss: {', '.join(f'{v:.3e}' for v in st)}")
    _p("[St] residuals: " + ", ".join(f"{k} {v:.3e}" for k, v in res.items()))

    # (En): steady wall flux equals the core's source per metre.
    bal = -K_GEL * b * _d1(b)
    gates["En_steady"] = abs(bal - Q_CORE * a**2 / 2) / (Q_CORE * a**2 / 2)
    _p(f"[En] -k_g b theta'(b) = {bal:.7g}, q_c a^2/2 = {Q_CORE * a**2 / 2:.7g} W/m, rel = {gates['En_steady']:.3e}")

    # (Lim): gel core equals the homogeneous J0 series, source in r < a.
    lim = []
    for t in (T_C, T_STAR):
        for r in R_PTS:
            comp = _theta(r, t, k_r=K_GEL, rhoc_r=RHO_C)
            ref = _theta_j0(r, t)
            lim.append(abs(comp - ref) / scale)
            _p(f"[Lim] r = {r * 1e3:.1f} mm, t = {t:.7g} s: composite {comp:.7g}, J0 series {ref:.7g}, "
               f"diff {lim[-1]:.3e}")
    gates["Lim_theta"] = max(lim)

    # (Bq): N = 200 vs 400.
    bq = [abs(_theta(r, t, n_modes=200) - _theta(r, t)) / scale for t in (T_C, T_STAR) for r in R_PTS]
    gates["Bq"] = max(bq)
    _p(f"[Bq] N=200 vs 400: {', '.join(f'{v:.3e}' for v in bq)}")

    # Control: k_r → K_GEL alone (ρc_r = ROD_RHO_C), steady θ on the axis.
    ti0_series = _theta(0.0, T_ST)
    kc0_series = _theta(0.0, T_ST, k_r=K_GEL, rhoc_r=ROD_RHO_C)
    ctrl = abs(kc0_series - ti0_series) / abs(ti0_series)
    ctrl_cf = abs(_theta_ss(0.0, k_r=K_GEL) - _theta_ss(0.0)) / _theta_ss(0.0)
    miss_a = abs(_theta(a, T_ST, k_r=K_GEL, rhoc_r=ROD_RHO_C) - _theta(a, T_ST)) / scale
    gates["ctrl"] = ctrl
    _p(f"[ctrl] k_r->K_GEL alone, steady (series at T_ST) theta(0): Ti {ti0_series:.7g}, swapped {kc0_series:.7g}, "
       f"rel miss {ctrl:.4g} (closed form {ctrl_cf:.4g}), predicted {PRED_CTRL}")
    _p(f"[ctrl] at r = a the steady miss is {miss_a:.3e} of scale (0 by construction: theta_ss(a) has no k_r)")

    # Printed record.
    for name, kw in (("Ti", dict()), ("gel", dict(k_r=K_GEL, rhoc_r=RHO_C))):
        for t, tl in ((T_C, "t_c"), (T_STAR, "t*"), (np.inf, "steady")):
            vals = [f"{_theta(r, t, **kw):.7g}" for r in R_PTS]
            _p(f"[rec] {name:4s} {tl:6s} theta(0, a, (a+b)/2) = {', '.join(vals)} K per W/m")
    ti_c = _theta(0.0, T_C)
    rc_c = _theta(0.0, T_C, rhoc_r=RHO_C)
    _p(f"[rec] rhoc-only swap (rhoc_r -> RHO_C) at t_c on the axis: Ti {ti_c:.7g}, swapped {rc_c:.7g}, "
       f"rel miss {abs(rc_c - ti_c) / abs(ti_c):.4g}")
    _p(f"[rec] theta(a, steady) per W/m: source in core {_theta_ss(a):.7g}, source in shell (item 843) "
       f"{float(_theta_ss_shell(a)):.7g}, ratio {_theta_ss(a) / float(_theta_ss_shell(a)):.7g}")

    _p("[gates] " + ", ".join(f"{k} = {v:.3e}" for k, v in gates.items()))
    assert gates["St_series"] <= B1_BAND, gates
    assert gates["St_resid"] <= B1_BAND, gates
    assert gates["En_steady"] <= B1_BAND, gates
    assert gates["Lim_theta"] <= B1_BAND, gates
    assert gates["Bq"] <= B1_BAND, gates
    assert gates["ctrl"] > INTERIOR_BAND, gates
