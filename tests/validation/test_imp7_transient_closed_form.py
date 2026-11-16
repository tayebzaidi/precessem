"""`IMP-7` step 4 (§9 item 658, §10 I45): the transient of the step 0 implant
heating problem on the adiabatic hole, in closed form (hollow-sphere
eigen-expansion, Carslaw & Jaeger 2e ch. IX; Jin 3e App. E).  Numpy / scipy
only.

Brief: docs/planning/chunks/IMP-7.md, "§9 item 658 — 2026-10-27".  Fixture
constants are step 0's (`test_imp7_implant_heating_closed_form.py`), ρc is
`MAT-9`'s (via `MAT-12`), θ_init = T∞.

Asserted at `B1_BAND`: per mode (PDE), (N), (Rb), (Orth); series (Lim),
(Duh).  Asserted as truncation convergence (2026-10-27 03:00 review: they
cannot reach 1e-12 pointwise): (IC), (Slope) — L2(r²) error at 400 terms
≤ ½ × that at 100; all three truncations printed.  Controls (Dirichlet at a,
isotropic source) asserted only as missing by > `CONTROL_MISS`.

Root count: the brief's asymptote ``(n + l/2)π/(b − a)`` is replaced by the
measured one.  For large λ the Neumann-at-a mode is ``cos λ(r − a)/r`` (the
``lπ/2`` phases of ``j_l``, ``y_l`` cancel between the two ends) and Robin at b
tends to Neumann, so the 0-indexed n-th root tends to ``nπ/(b − a)`` for every
l; for l = 2 the brief's ladder is exactly one spacing off.  Asserted: within
½ spacing of ``nπ/(b − a)`` for n ≥ 10, and every gap > ``π/(b − a)/8`` (so
no grid bracket of the finder can hold two roots).
"""

from __future__ import annotations

import numpy as np
from mpi4py import MPI

from fem_em_solver.utils.analytical import (
    _imp7_gl_nodes,
    _imp7_mode_matrix,
    pec_sphere_implant_duty_cycle_from_modes,
    pec_sphere_implant_eigenvalues,
    pec_sphere_implant_mode,
    pec_sphere_implant_radial,
    pec_sphere_implant_source_radial,
    pec_sphere_implant_steady_temperature,
    pec_sphere_implant_transient_from_modes,
    pec_sphere_implant_transient_modes,
    pec_sphere_implant_transient_temperature,
)
from tests.validation.test_imp4_pec_cylinder_closed_form import B1_BAND, CONTROL_MISS
from tests.validation.test_imp7_implant_heating_closed_form import (
    A,
    B,
    BETA,
    FIX,
    H,
    K_GEL,
    Q_FAR,
    T_INF,
    _shell_points,
)
from tests.validation.test_mat12_transient_sar_heating_sphere import RHO_C

SPAN = B - A
KAPPA = K_GEL / RHO_C
T_STAR = 0.2 * SPAN**2 / KAPPA  # brief: t* = 0.2 (b − a)² / κ
N_DEFAULT = 200
TRUNCATIONS = (100, 200, 400)
CONV_RATIO = 0.5  # brief: error at 400 ≤ ½ × error at 100 (a record, not a band)
N_PDE, N_ORTH = 50, 20


def _p(*args):
    if MPI.COMM_WORLD.rank == 0:
        print(*args, flush=True)


def _modes(n=N_DEFAULT, inner="neumann"):
    return pec_sphere_implant_transient_modes(beta=BETA, n_terms=n, inner=inner, **FIX)


def _pts(n=256, seed=4):
    return _shell_points(n, A * (1.0 + 1e-9), B, seed)


def _theta(pts, t, modes, ls=(0, 2)):
    return pec_sphere_implant_transient_from_modes(pts, t, modes, rho_c=RHO_C,
                                                   t_inf=T_INF, ls=ls) - T_INF


def test_imp7_4_eigenvalue_count():
    for inner in ("neumann", "dirichlet"):
        for l in (0, 2):
            lam = pec_sphere_implant_eigenvalues(l, 400, radius=A, outer_radius=B,
                                                 k=K_GEL, h=H, inner=inner)
            x = lam * SPAN / np.pi
            n = np.arange(lam.size)
            off = 0.0 if inner == "neumann" else 0.5
            dev = np.abs(x[10:] - (n[10:] + off))
            gap = np.min(np.diff(x))
            _p(f"(Count) {inner} l={l}: λ₀(b−a)/π = {x[0]:.6f}, λ₁ {x[1]:.6f}, "
               f"λ₃₉₉ {x[-1]:.6f}; max |x_n − (n+{off})| (n≥10) = {dev.max():.3e}; "
               f"min gap = {gap:.6f} spacings; brief ladder (n+l/2) dev at n=399: "
               f"{x[-1] - (399 + l / 2):+.4f}")
            assert dev.max() < 0.5
            assert gap > 1.0 / 8.0
            assert lam[0] > 0.0


def test_imp7_4_per_mode_identities():
    r = np.linspace(A, B, 401)
    for l in (0, 2):
        lam = pec_sphere_implant_eigenvalues(l, N_PDE, radius=A, outer_radius=B,
                                             k=K_GEL, h=H)
        worst = {"PDE": 0.0, "N": 0.0, "Rb": 0.0}
        for lm in lam:
            R0, R1, R2 = (pec_sphere_implant_mode(r, lm, l, radius=A, derivative=d)
                          for d in (0, 1, 2))
            s = np.max(np.abs(R0))
            res = R2 + 2.0 * R1 / r + (lm**2 - l * (l + 1) / r**2) * R0
            worst["PDE"] = max(worst["PDE"], np.max(np.abs(res)) / (lm**2 * s))
            worst["N"] = max(worst["N"], abs(R1[0]) / (lm * s))
            worst["Rb"] = max(worst["Rb"], abs(K_GEL * R1[-1] + H * R0[-1])
                              / ((K_GEL * lm + H) * s))
        _p(f"(PDE)/(N)/(Rb) l={l}, first {N_PDE} modes: "
           + ", ".join(f"{k} {v:.3e}" for k, v in worst.items())
           + f"  [B1_BAND {B1_BAND:.0e}]")
        for v in worst.values():
            assert v <= B1_BAND

    modes = _modes(N_ORTH)
    rq, wq = _imp7_gl_nodes(A, B, 400, 16)
    wr = wq * rq * rq
    for l in (0, 2):
        R = _imp7_mode_matrix(rq, modes, l)
        G = (R * wr) @ R.T
        d = np.sqrt(np.diag(G))
        G = G / np.outer(d, d)
        off = np.max(np.abs(G - np.diag(np.diag(G))))
        _p(f"(Orth) l={l}, first {N_ORTH} modes: max |⟨R_n,R_m⟩|/norms = {off:.3e}")
        assert off <= B1_BAND


def test_imp7_4_limit_and_duhamel():
    pts = _pts()
    ss = pec_sphere_implant_steady_temperature(pts, k_core=0.0, t_inf=T_INF,
                                               beta=BETA, **FIX) - T_INF
    th_inf = pec_sphere_implant_transient_temperature(
        pts, 1e3 * T_STAR, rho_c=RHO_C, t_inf=T_INF, beta=BETA, n_terms=N_DEFAULT,
        **FIX) - T_INF
    lim = np.max(np.abs(th_inf - ss)) / np.max(np.abs(ss))
    _p(f"(Lim) t = 1e3 t* = {1e3 * T_STAR:.4g} s: max|θ − θ_ss|/max|θ_ss| = {lim:.3e}"
       f"  (θ_ss max {np.max(ss):.6g} K)")
    assert lim <= B1_BAND

    modes = _modes()
    period = 0.5 * T_STAR
    t_on = 0.5 * period
    for t in (1.6 * T_STAR, 1.4 * T_STAR):  # on-phase, off-phase
        direct = pec_sphere_implant_duty_cycle_from_modes(
            pts, t, modes, rho_c=RHO_C, t_inf=T_INF, period=period) - T_INF
        duh = np.zeros(len(pts))
        k = 0
        while k * period < t:
            duh += _theta(pts, t - k * period, modes)
            if t - k * period - t_on > 0.0:
                duh -= _theta(pts, t - k * period - t_on, modes)
            k += 1
        rel = np.max(np.abs(duh - direct)) / np.max(np.abs(direct))
        _p(f"(Duh) t = {t / T_STAR:.2f} t*, {k} periods: max|Duhamel − modal|/max = "
           f"{rel:.3e}  (θ max {np.max(direct):.6g} K)")
        assert rel <= B1_BAND


def test_imp7_4_ic_and_slope_converge():
    rq, wq = _imp7_gl_nodes(A, B, 1000, 16)
    wr = wq * rq * rq
    err = {(name, l): [] for name in ("IC", "Slope") for l in (0, 2)}
    for n in TRUNCATIONS:
        modes = _modes(n)
        for l in (0, 2):
            m = modes[l]
            R = _imp7_mode_matrix(rq, modes, l)
            ss = pec_sphere_implant_radial(rq, modes["coeffs"], l)
            ic = ss - m["c"] @ R  # θ_l(t = 0) − θ_init,l
            ql = pec_sphere_implant_source_radial(rq, l, radius=A, q_far=Q_FAR,
                                                  beta=BETA)
            slope = K_GEL * (m["c"] * m["lam"] ** 2) @ R  # ρc ∂ₜθ_l at 0⁺
            err[("IC", l)].append(np.sqrt(np.sum(ic**2 * wr) / np.sum(ss**2 * wr)))
            err[("Slope", l)].append(np.sqrt(np.sum((slope - ql) ** 2 * wr)
                                             / np.sum(ql**2 * wr)))
    for (name, l), e in err.items():
        _p(f"({name}) l={l}: rel L2(r²) at n = {TRUNCATIONS}: "
           + " / ".join(f"{v:.3e}" for v in e)
           + f"; ratio 400/100 = {e[-1] / e[0]:.4f} (≤ {CONV_RATIO})")
    for e in err.values():
        assert e[-1] <= CONV_RATIO * e[0]


def test_imp7_4_controls_miss():
    pts = _pts(seed=5)
    full = _theta(pts, T_STAR, _modes())
    dirich = _theta(pts, T_STAR, _modes(inner="dirichlet"))
    iso = _theta(pts, T_STAR, _modes(), ls=(0,))

    def rel(x):
        return np.sqrt(np.sum((x - full) ** 2) / np.sum(full**2))

    _p(f"t* = {T_STAR:.4f} s; adiabatic θ(t*) max {np.max(full):.6g} K")
    _p(f"(c-D) Dirichlet-at-a miss rel L2 = {rel(dirich):.4e} (> {CONTROL_MISS})")
    _p(f"(c-iso) isotropic-source miss rel L2 = {rel(iso):.4e} (> {CONTROL_MISS})")
    assert rel(dirich) > CONTROL_MISS
    assert rel(iso) > CONTROL_MISS
