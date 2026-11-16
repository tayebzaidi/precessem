"""`IMP-8` step 16 (§9 item 874, §10 I78): the lumped F2182-shaped prediction,
closed form, numpy / scipy only.

Brief: docs/planning/chunks/IMP-8.md, "§9 item 874 — 2026-11-08".  The two
idealised source shapes of steps 0 / 1 (the zero-radius line Φ and the tip pair
Φ₂, same total power, 1 W per rod through ``SCALE`` = 1/(4πkL)) in an infinite
gel medium, read on the axis 1 mm beyond each tip, (r, z) = (0, ±Z_P),
Z_P = L/2 + 1 mm, at t = 900 s, continuous and at D = 0.5 on item 787's
``PERIOD`` (450 s; the two-tip duty runs I62 / I63 import the same period, and
the weekly's "period as I56's" is read as that).  At t = 900 s the edge at 900 s
is excluded by both the step sum (t > tⱼ) and the convolution (j·P < t).

Asserted (bands imported, unmoved):
- (Lim) Φ(0, Z_P, ∞) against the on-axis steady form ln((Z_P + L/2)/(Z_P − L/2)),
  relative ≤ ``B1_BAND``.  The axis form is the r → 0 limit of ``_phi_steady``,
  written here because ``_phi_steady`` / I58's ``_line(·, ·, inf)`` divide by r.
- (Bq) a Gauss–Legendre rule in z′ (I58's ``_gl01``) on ``_point`` at n = 64 vs
  2n = 128, t = ``T_PRINT`` and ∞: |GL_n/GL_2n − 1| ≤ ``BQ_BAND``.
- (Dut) the D = 0.5 step sum at t = 900 s against the from-rest convolution, line
  and twin, |difference| / θ_ss(probe) ≤ ``B1_BAND``.
- Control (the weekly's, factor 1 on ``INTERIOR_BAND``): the probe moved to
  Z_C = L/2 + 2 mm reads lower under the continuous source, both shapes:
  INTERIOR_BAND < 1 − θ(Z_C)/θ(Z_P) < 1.

Printed (🟡 record, no band, K per watt per rod): the F2182-shaped probe
(placement and time quoted from memory; W3 confirms or corrects).

Scope: infinite gel, a zero-radius line and two point tips, one placement.  Not
claimed: box walls, the rod's own conduction (I73 / I75), a ball-averaged read,
any F2182 figure or published value.  The twin's point source sits 1 mm from the
probe (a point's 1/R there); the FEM's 1 mm balls (I56) put their surface at the
probe, so the twin is not the FEM's source shape at that distance — printed,
never compared.
"""

from __future__ import annotations

import time

import numpy as np
from mpi4py import MPI

from tests.validation.test_imp4_pec_cylinder_closed_form import B1_BAND
from tests.validation.test_imp8_entry48_lumped import BQ_BAND, _gl01, _tips_steady
from tests.validation.test_imp8_rod_line_source_closed_form import (
    KAPPA,
    MM,
    QUAD,
    T_PRINT,
    L,
    _phi,
    _phi_steady,
    _phi_tips,
    _point,
)
from tests.validation.test_imp8_rod_line_source_duty_closed_form import PERIOD, SCALE, T_STAR, _phi_dot
from tests.validation.test_imp8_rod_tips_duty_closed_form import _convolution, _duhamel, _phi_tips_dot
from tests.validation.test_mat13_transient_perfusion_sphere import INTERIOR_BAND

Z_P = 0.5 * L + 1.0 * MM  # 16 mm, the probe 1 mm beyond each tip (quoted from memory)
Z_C = 0.5 * L + 2.0 * MM  # 17 mm, the control placement
DUTY = 0.5
N_GL = 64  # (Bq): n vs 2n


def _p(*args):
    if MPI.COMM_WORLD.rank == 0:
        print(*args, flush=True)


def _axis_steady(z, length=L):
    """r → 0 limit of ``_phi_steady`` for |z| > L/2."""
    a = abs(z)
    return np.log((a + 0.5 * length) / (a - 0.5 * length))


def _gl_line(z, t, n, length=L):
    """Gauss–Legendre rule in z′ on [−L/2, L/2] of the line kernel at r = 0."""
    u, w = _gl01(n)
    zp = -0.5 * length + length * u
    return length * float(np.dot(w, _point(np.abs(z - zp), t)))


SHAPES = {
    "line": (_phi, _phi_dot, lambda z: _axis_steady(z)),
    "twin": (_phi_tips, _phi_tips_dot, lambda z: _tips_steady(0.0, z)),
}
#: The control's steady-form prediction (review arithmetic: line ≈ 0.19, twin ≈ 0.49), printed only.
CTRL_PRED = {n: 1.0 - s(Z_C) / s(Z_P) for n, (_k, _kd, s) in SHAPES.items()}


def test_lim_steady_axis_form():
    phi_inf = _phi(0.0, Z_P, np.inf)
    ax = _axis_steady(Z_P)
    rel = abs(phi_inf / ax - 1.0)
    tips_rel = abs(_phi_tips(0.0, Z_P, np.inf) / _tips_steady(0.0, Z_P) - 1.0)
    off_rel = abs(_phi_steady(1e-9, Z_P) / ax - 1.0)
    _p(f"\n[imp8-16] Z_P {Z_P / MM:g} mm, Z_C {Z_C / MM:g} mm, PERIOD {PERIOD:g} s, D {DUTY:g}, "
       f"t {T_STAR:g} s, KAPPA {KAPPA:.7g} m^2/s, SCALE {SCALE:.7g} K/W per unit Φ")
    _p(f"[imp8-16] (Lim) Φ(0,Z_P,inf) {phi_inf:.7g} vs ln(31) {ax:.7g}: rel {rel:.3e} (band {B1_BAND:g}); "
       f"printed: Φ2(0,Z_P,inf)/tips_steady - 1 = {tips_rel:.3e}; _phi_steady(1e-9,Z_P)/axis - 1 = "
       f"{off_rel:.3e}")
    assert rel <= B1_BAND


def test_bq_gauss_legendre_n_vs_2n():
    worst = 0.0
    for t in (T_PRINT, np.inf):
        g_n, g_2n = _gl_line(Z_P, t, N_GL), _gl_line(Z_P, t, 2 * N_GL)
        e = abs(g_n / g_2n - 1.0)
        _p(f"[imp8-16] (Bq) t={t:g}: GL{N_GL} {g_n:.7g} GL{2 * N_GL} {g_2n:.7g} |n/2n-1| {e:.3e} "
           f"(band {BQ_BAND:g}); printed |GL2n/_phi-1| {abs(g_2n / _phi(0.0, Z_P, t) - 1.0):.3e}")
        worst = max(worst, e)
    assert worst <= BQ_BAND


def test_dut_step_sum_equals_convolution():
    t0 = time.perf_counter()
    worst = 0.0
    for name, (kern, kern_dot, steady) in SHAPES.items():
        step = _duhamel(kern, 0.0, Z_P, T_STAR, DUTY)
        conv = _convolution(kern_dot, 0.0, Z_P, T_STAR, DUTY)
        e = abs(step - conv) / (SCALE * steady(Z_P))
        _p(f"[imp8-16] (Dut) {name}: step sum {step:.7g} conv {conv:.7g} |diff|/θ_ss {e:.3e} "
           f"(band {B1_BAND:g})")
        worst = max(worst, e)
    _p(f"[imp8-16] (Dut) worst {worst:.3e}; {time.perf_counter() - t0:.1f} s")
    assert worst <= B1_BAND


def _reads(kern, steady, z):
    cont = SCALE * kern(0.0, z, T_STAR)
    duty = _duhamel(kern, 0.0, z, T_STAR, DUTY)
    return cont, duty, SCALE * steady(z)


def test_control_probe_moved_reads_lower():
    for name, (kern, _kd, _st) in SHAPES.items():
        hi = SCALE * kern(0.0, Z_P, T_STAR)
        lo = SCALE * kern(0.0, Z_C, T_STAR)
        drop = 1.0 - lo / hi
        _p(f"[imp8-16] (ctrl) {name}: θ(Z_C) {lo:.7g} θ(Z_P) {hi:.7g} 1-ratio {drop:.7g} "
           f"(band > {INTERIOR_BAND:g}, < 1; predicted steady {CTRL_PRED[name]:.7g})")
        assert lo < hi
        assert INTERIOR_BAND < drop < 1.0


def test_printed_record():
    _p("[imp8-16] the F2182-shaped probe (placement and time quoted from memory; W3 confirms or "
       "corrects) — K per watt per rod, infinite gel")
    rec = {}
    for name, (kern, _kd, steady) in SHAPES.items():
        for lab, z in (("Z_P", Z_P), ("Z_C", Z_C)):
            cont, duty, ss = _reads(kern, steady, z)
            rec[(name, lab)] = (cont, duty, ss)
            _p(f"[imp8-16] {name} at (0, +{lab}={z / MM:g} mm): θ(900 s) cont {cont:.7g}  D=0.5 {duty:.7g}  "
               f"θ_ss {ss:.7g}  duty/cont {duty / cont:.7g}")
        sym = abs(SCALE * kern(0.0, Z_P, T_STAR) - SCALE * kern(0.0, -Z_P, T_STAR))
        sym_d = abs(_duhamel(kern, 0.0, Z_P, T_STAR, DUTY) - _duhamel(kern, 0.0, -Z_P, T_STAR, DUTY))
        _p(f"[imp8-16] {name}: |θ(+Z_P)-θ(-Z_P)| cont {sym:.3e} D=0.5 {sym_d:.3e} (round-off by symmetry)")
    for lab in ("Z_P", "Z_C"):
        tw, ln = rec[("twin", lab)], rec[("line", lab)]
        _p(f"[imp8-16] twin/line at {lab}: cont {tw[0] / ln[0]:.7g}  D=0.5 {tw[1] / ln[1]:.7g}  "
           f"ss {tw[2] / ln[2]:.7g}")
    _p(f"[imp8-16] QUAD {QUAD}")
    for v in rec.values():
        assert all(np.isfinite(x) and x > 0.0 for x in v)
