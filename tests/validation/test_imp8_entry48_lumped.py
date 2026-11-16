"""`IMP-8` step 6 (§9 item 786, §10 I58): the lumped reading of entry 48's rod,
infinite-medium closed form, numpy / scipy only.

Brief: docs/planning/chunks/IMP-8.md, "§9 item 786 — 2026-11-03".  Step 0's
line kernel Φ and step 1's tip pair Φ₂ (`test_imp8_rod_line_source_closed_form`,
imported, never restated) read at entry 48's points and power.  The tip 1 g
ball's centre lies on the line's axis at its end, where Φ and Φ₂ are singular,
so the lumped θ at a tip is read two ways, both printed: as the 1 g ball
average (radius `averaging_ball_radius` at the gel ρ) and at the rod-surface
point (r = ROD_RADIUS, z = ±L/2).

The ball average is a product Gauss–Legendre rule in (ρ, μ = cos ϑ, φ) about
the ball centre, with polynomial clustering ρ = a u^Q_RHO at the centre and
δ = s^P_MU at each pole (μ = ±1), where the axis line's log singularity sits.

Asserted: (Bq) the ball average at n and 2n agrees to ≤ BQ_BAND, beside the
point kernel's ball average about its own source = 3 / (2a) to ≤ BQ_BAND (the
singular-integrand check that gives (Bq) its meaning); (Mono) θ(t*) < θ∞ at
every point read; (En) and (Far) re-run as imported (the kernel's own gates;
they do not take P).  Control: the 1 g ball moved one ball radius off the axis
changes the ball average by > OFF_AXIS_MIN against the on-axis one.

Printed, as a record (no band): θ_line and θ_tips at t* = 900 s per watt and at
P_rod, at both tip balls (average and surface point), at P37's six points in
the clinical box (≈ 0: √(κt*) ≈ 11 mm against ≥ 0.1 m), and the ratio
θ_tips / θ_line at the tip ball.

Scope: an infinite-medium lumped shape at entry 48's power.  Not claimed: the
box walls, the rod's own conduction, the FEM tip θ (that is entry 55).
"""

from __future__ import annotations

import numpy as np
from mpi4py import MPI

from fem_em_solver.post.sar import averaging_ball_radius

from tests.mesh.test_birdcage_implant_rod import ROD_LENGTH, ROD_RADIUS
from tests.validation import test_geo36_rod_f_clinical as i35
from tests.validation.test_birdcage_sar_mass_averaged import PHANTOM_RHO_KG_PER_M3
from tests.validation.test_geo36_rod_f_clinical_drive import _probe_points
from tests.validation.test_geo36_rod_f_clinical_heating import _p37_fixture
from tests.validation.test_imp8_rod_line_source_closed_form import (
    K_GEL,
    L,
    T_PRINT,
    _phi,
    _phi_steady,
    _phi_tips,
)
from tests.validation.test_imp8_rod_line_source_closed_form import (  # noqa: F401  (re-run as imported)
    test_en_energy_identities,
    test_far_steady_limit,
)
from tests.validation.test_mass_averaged_sar_standard_masses import ONE_GRAM_KG
from tests.validation.test_mat21_box_probe_points_fem import PROBE_LABELS

#: Entry 48's rod power on (1, +1j), arm M, the log's own `.9e` string
#: (`20261102T080008Z_GEO-36-step9-I36.log:12494`).  Leak guard: printed as is.
P_ROD_LOG = "9.495544579e-08"
P_ROD = float(P_ROD_LOG)
#: The rod centre and its tips, `20261031T003859Z_GEO-36-I36c-W0b.log:90`.
CENTRE_LOG = (0.105, 0.0, 0.0)
TIP_HALF_LOG = 0.015
APEX_HALF_LOG = 0.016

BQ_BAND = 1.0e-6  # brief, (Bq)
OFF_AXIS_MIN = 1.0e-2  # brief, the control
N_BQ = 16  # the rule's n; (Bq) compares n with 2n
Q_RHO = 2  # ρ = a u^Q_RHO: ρ² log ρ → u^(3Q-1) log u at the centre
P_MU = 4  # δ = 1 ∓ μ = s^P_MU at each pole: log δ → s^(P-1) log s
N_PHI_OFF = 16  # φ nodes on [0, π] for the off-axis (non-axisymmetric) ball
N_OFF = 8  # (ρ, μ-half) nodes of the off-axis ball: the control needs ≈ 1e-3, not 1e-6
PER_WATT = 1.0 / (4.0 * np.pi * K_GEL * L)  # θ = P / (4πkL) · Φ, both shapes


def _p(*args):
    if MPI.COMM_WORLD.rank == 0:
        print(*args, flush=True)


def _gl01(n):
    x, w = np.polynomial.legendre.leggauss(n)
    return 0.5 * (x + 1.0), 0.5 * w


def _radial_nodes(a, n):
    """ρ nodes / weights of ∫₀^a ρ² g dρ, clustered at the centre."""
    u, w = _gl01(n)
    rho = a * u**Q_RHO
    return rho, w * rho**2 * a * Q_RHO * u ** (Q_RHO - 1)


def _mu_nodes(n):
    """μ nodes / weights on [−1, 1], two halves, each clustered at its pole."""
    s, w = _gl01(n)
    delta = s**P_MU
    jac = w * P_MU * s ** (P_MU - 1)
    return np.concatenate([-1.0 + delta, 1.0 - delta]), np.concatenate([jac, jac])


def _ball_average_axis(f, z_c, a, n):
    """Average of f(r, z) over the ball of radius a centred on the axis at z_c."""
    rho, wr = _radial_nodes(a, n)
    mu, wm = _mu_nodes(n)
    total = 0.0
    for rr, wrr in zip(rho, wr, strict=True):
        sin_t = np.sqrt(np.clip(1.0 - mu * mu, 0.0, None))
        vals = [f(rr * st, z_c + rr * m) for st, m in zip(sin_t, mu, strict=True)]
        total += wrr * float(np.dot(wm, vals))
    return total * 2.0 * np.pi / (4.0 * np.pi * a**3 / 3.0)


def _ball_average_off(f, r_c, z_c, a, n, n_phi=N_PHI_OFF):
    """Average of f(r, z) over the ball of radius a centred at (r_c, 0, z_c)."""
    rho, wr = _radial_nodes(a, n)
    mu, wm = _mu_nodes(n)
    ph, wp = _gl01(n_phi)
    ph, wp = np.pi * ph, np.pi * wp  # φ ∈ [0, π], the other half by y → −y
    total = 0.0
    for rr, wrr in zip(rho, wr, strict=True):
        for m, wmm in zip(mu, wm, strict=True):
            st = np.sqrt(max(1.0 - m * m, 0.0))
            x = r_c + rr * st * np.cos(ph)
            y = rr * st * np.sin(ph)
            vals = [f(np.hypot(xx, yy), z_c + rr * m) for xx, yy in zip(x, y, strict=True)]
            total += wrr * wmm * float(np.dot(wp, vals))
    return total * 2.0 / (4.0 * np.pi * a**3 / 3.0)


def _a_1g():
    return averaging_ball_radius(mass_kg=ONE_GRAM_KG, rho=PHANTOM_RHO_KG_PER_M3)


def _line(r, z, t=T_PRINT):
    return _phi(r, z, t) if np.isfinite(t) else _phi_steady(r, z)


def _tips(r, z, t=T_PRINT):
    return _phi_tips(r, z, t)


def _tips_steady(r, z):
    return 0.5 * L * sum(1.0 / np.hypot(r, z - s * 0.5 * L) for s in (1.0, -1.0))


def test_inputs_from_the_logs():
    centre = i35._rod_centre()
    tips, _t1 = _probe_points(centre)
    _p(f"[imp8 I58] P_rod = {P_ROD_LOG} W (log string); centre {tuple(map(float, centre))}; "
       f"tips {tips['top']} / {tips['bottom']}; apexes {tips['top_apex']} / {tips['bottom_apex']}")
    assert np.allclose(centre, CENTRE_LOG, rtol=0.0, atol=1e-15)
    assert np.isclose(tips["top"][2] - centre[2], TIP_HALF_LOG, rtol=0.0, atol=1e-15)
    assert np.isclose(centre[2] - tips["bottom"][2], TIP_HALF_LOG, rtol=0.0, atol=1e-15)
    assert np.isclose(tips["top_apex"][2] - centre[2], APEX_HALF_LOG, rtol=0.0, atol=1e-15)
    assert L == ROD_LENGTH and 0.5 * L == TIP_HALF_LOG


def test_bq_kernel_identity():
    """1/R averaged over a ball centred on its source = 3 / (2a), by the same rule."""
    a = _a_1g()
    for n in (N_BQ, 2 * N_BQ):
        avg = _ball_average_axis(lambda r, z: 1.0 / np.hypot(r, z), 0.0, a, n)
        rel = abs(avg / (1.5 / a) - 1.0)
        _p(f"[imp8 I58] (Bq-k) n={n} a={a:.7g} m avg={avg:.9e} 3/(2a)={1.5 / a:.9e} rel={rel:.3e} "
           f"band={BQ_BAND:.0e}")
        assert rel <= BQ_BAND


def test_bq_mono_and_record():
    a = _a_1g()
    half = 0.5 * L
    rows = {}
    for label, z_c in (("top", half), ("bottom", -half)):
        for name, f, f_inf in (("line", _line, lambda r, z: _phi_steady(r, z)),
                               ("tips", _tips, _tips_steady)):
            lo = _ball_average_axis(f, z_c, a, N_BQ)
            hi = _ball_average_axis(f, z_c, a, 2 * N_BQ)
            rel = abs(lo / hi - 1.0)
            steady = _ball_average_axis(f_inf, z_c, a, 2 * N_BQ)
            _p(f"[imp8 I58] (Bq) {label} ball {name}: n={N_BQ} {lo:.9e} 2n={hi:.9e} rel={rel:.3e} "
               f"band={BQ_BAND:.0e}; steady {steady:.7g}")
            assert rel <= BQ_BAND
            assert hi < steady, (label, name, hi, steady)  # (Mono)
            rows[(label, name, "ball")] = hi
        r_s, z_s = ROD_RADIUS, z_c
        for name, f, f_inf in (("line", _line, _phi_steady), ("tips", _tips, _tips_steady)):
            v, s = f(r_s, z_s), f_inf(r_s, z_s)
            assert v < s, (label, name, v, s)  # (Mono)
            rows[(label, name, "surface")] = v
    _p(f"[imp8 I58] 1 g ball radius a = {a:.7g} m (ρ = {PHANTOM_RHO_KG_PER_M3:g} kg/m³)")
    for label in ("top", "bottom"):
        for where in ("ball", "surface"):
            ln, tp = rows[(label, "line", where)], rows[(label, "tips", where)]
            _p(f"[imp8 I58] (print) {label} tip {where}: Φ {ln:.7g} Φ₂ {tp:.7g}; "
               f"θ_line {ln * PER_WATT:.7g} K/W = {ln * PER_WATT * P_ROD:.7g} K at P_rod; "
               f"θ_tips {tp * PER_WATT:.7g} K/W = {tp * PER_WATT * P_ROD:.7g} K at P_rod; "
               f"θ_tips/θ_line {tp / ln:.7g}")


def test_six_box_points_and_mono():
    centre = i35._rod_centre()
    probes = _p37_fixture(centre[1])
    for lab in PROBE_LABELS:
        x = np.asarray(probes[lab], dtype=float)
        r = float(np.hypot(x[0] - centre[0], x[1] - centre[1]))
        z = float(x[2] - centre[2])
        ln, tp = _line(r, z), _tips(r, z)
        ln_s, tp_s = _phi_steady(r, z), _tips_steady(r, z)
        _p(f"[imp8 I58] (print) {lab} at {tuple(float(v) for v in x)} (r {r:.4g}, z {z:.4g}) m: "
           f"θ_line {ln * PER_WATT:.3e} K/W = {ln * PER_WATT * P_ROD:.3e} K; "
           f"θ_tips {tp * PER_WATT:.3e} K/W = {tp * PER_WATT * P_ROD:.3e} K (≈ 0); "
           f"steady {ln_s * PER_WATT:.4g} / {tp_s * PER_WATT:.4g} K/W")
        assert ln < ln_s and tp < tp_s  # (Mono)


def test_control_off_axis_ball():
    a = _a_1g()
    on = _ball_average_axis(_line, 0.5 * L, a, 2 * N_BQ)
    off = _ball_average_off(_line, a, 0.5 * L, a, N_OFF)
    change = abs(off / on - 1.0)
    _p(f"[imp8 I58] (c-off) on-axis {on:.7g} off-axis {off:.7g} change={change:.3e} "
       f"> {OFF_AXIS_MIN:.0e}")
    assert change > OFF_AXIS_MIN
