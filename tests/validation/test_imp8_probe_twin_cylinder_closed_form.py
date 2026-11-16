"""IMP-8 step 21 (§9 item 934, §10 I86): the finite-radius two-tip twin at the probe, closed form.

Brief: docs/planning/chunks/IMP-8.md, "§9 item 934 — 2026-11-10".  I80's point twin Φ₂ puts P/2
at each tip as a point; here each tip is a uniform cylinder of radius a = ``ROD_RADIUS`` and length
h = 2a centred on z = ±L/2 (the fixture is the 03:00 review's execution detail, flagged for the
weekly), carrying P/2.  Per unit ``SCALE`` (θ = SCALE·Φ₂,a), on the axis,

    Φ₂,a(0, z, t) = (L/2)·[M(|z − L/2|, t) + M(|z + L/2|, t)],
    M(ζ, t) = ½∫₋₁¹ B(ζ − c u, t) du,   c = h/2 (the segment mean in the unit variable u),
    B(d, t) = (2s/a²)[ierfc(d/s) − ierfc(D/s)],  D = hypot(a, d),  s = 2√(κt),

I83's disc integral (``test_imp8_probe_cylinder_source_closed_form.py:5–13``), normalised by πa² so
that B → erfc(d/s)/d as a → 0 and Φ₂,a → ``_phi_tips``.  At t = ∞, B = 2/(D + d) and
M(ζ, ∞) = [g(ζ + c) − g(ζ − c)]/(2c) with I83's ``_g``.  The rate per tip is
(L/(2h))·``_phi_cyl_dot``(0, ζ, u, a=a, length=h), used at a = ``ROD_RADIUS`` only.

Small-a cancellation removed: the bracket is evaluated as ∫ₓʸ erfc(v) dv by 16-point GL with the
width y − x written as a²/((D + d)s), never D/s − d/s; the segment mean is integrated in u (not z′)
so ``QUAD``'s epsabs never meets an O(h) integral.

Gates (bands imported, none set or moved): (Lim) B1_BAND, (Bq) BQ_BAND, (Dut) B1_BAND, (Lim-a)
B1_BAND.  Control (the weekly's, factor 1 on INTERIOR_BAND): at 16 mm, continuous θ(900 s),
1 − Φ₂,a/Φ₂ > INTERIOR_BAND.  (Lim-a)'s own negative: at a = ROD_RADIUS/10 the deficit > B1_BAND.

Scope: infinite gel, two uniform tip cylinders, 1 W per rod, four axis placements.  Not claimed:
I56's FEM balls, the box walls or the rod's own conduction, a ball-averaged read, any F2182 figure.
"""

from __future__ import annotations

import time
from functools import partial

import numpy as np
from scipy import special
from scipy.integrate import quad

from tests.mesh.test_birdcage_implant_rod import ROD_RADIUS
from tests.validation.test_imp4_pec_cylinder_closed_form import B1_BAND
from tests.validation.test_imp8_entry48_lumped import BQ_BAND, _gl01, _tips_steady
from tests.validation.test_imp8_probe_cylinder_source_closed_form import _GL_B, _g, _phi_cyl, _phi_cyl_dot
from tests.validation.test_imp8_probe_f2182_closed_form import DUTY, N_GL, _axis_steady, _p
from tests.validation.test_imp8_probe_ladder_closed_form import LADDER_MM, RUNGS, _lab
from tests.validation.test_imp8_rod_line_source_closed_form import KAPPA, MM, QUAD, T_PRINT, L, _phi, _phi_tips
from tests.validation.test_imp8_rod_line_source_duty_closed_form import PERIOD, SCALE, T_STAR
from tests.validation.test_imp8_rod_tips_duty_closed_form import _convolution, _duhamel, _phi_tips_dot
from tests.validation.test_mat13_transient_perfusion_sphere import INTERIOR_BAND

A_LIM = ROD_RADIUS * 1e-7  # 1e-10 m, (Lim-a)
#: I83's logged Φ_a/Φ on the line (20261109T110006Z_IMP-8-item911-W1.log:119/:125/:130/:135), printed beside.
I83_RADIUS_FACTOR = ("0.9532823", "0.9811647", "0.9889726", "0.9941423")
#: I80's logged point twin / line (20261109T043829Z_IMP-8-item894-W1.log:100/:103/:106/:109), printed beside.
I80_TWIN_LINE = ("6.543912", "4.400661", "3.60769", "2.940895")
#: The brief's drafter predictions (unlogged arithmetic), printed only.
PRED_CTRL, PRED_CTRL_SS, PRED_H_A = 0.044, 0.0408, (0.60, 0.63)


def _bracket_cf(d, t, a):
    """B(d, t): (2s/a²)[ierfc(d/s) − ierfc(D/s)] with y − x = a²/((D + d)s); 2/(D + d) at t = inf."""
    d = np.asarray(d, dtype=float)
    big_d = np.hypot(a, d)
    if np.isinf(t):
        return 2.0 / (big_d + d)
    s = 2.0 * np.sqrt(KAPPA * t)
    x = d / s
    width = a * a / ((big_d + d) * s)
    u, w = _GL_B
    v = x[..., None] + width[..., None] * u
    return (2.0 * s / (a * a)) * width * (special.erfc(v) @ w)


def _seg_mean(zeta, t, a, n=None, half=None):
    """M(ζ, t) = ½∫₋₁¹ B(ζ − c u, t) du, c = ``half`` (default a, i.e. h = 2a)."""
    c = a if half is None else half
    assert zeta >= c * (1.0 - 1e-12), (zeta, c)  # the probe is on or beyond the face
    if n is None:
        return 0.5 * quad(lambda u: float(_bracket_cf(np.array([zeta - c * u]), t, a)[0]),
                          -1.0, 1.0, **QUAD)[0]
    x, w = _gl01(n)
    return float(np.dot(w, _bracket_cf(zeta - c * (2.0 * x - 1.0), t, a)))


def _centres(length):
    """Tip-cylinder centres: ±L/2 for h = 2a; ±(L/2 − h/2) for the h = a inset (inside the rod's end)."""
    off = 0.0 if length is None else 0.5 * length
    return (0.5 * L - off, -(0.5 * L - off))


def _twin_a(r, z, t, a=ROD_RADIUS, n=None, length=None):
    """Φ₂,a(0, z, t); ``length`` None is the fixture h = 2a centred on ±L/2."""
    assert r == 0.0
    half = None if length is None else 0.5 * length
    return 0.5 * L * sum(_seg_mean(abs(z - zc), t, a, n=n, half=half) for zc in _centres(length))


def _twin_a_dot(r, z, u, a=ROD_RADIUS):
    """∂Φ₂,a/∂t at age u, closed form, h = 2a (use at a = ROD_RADIUS only)."""
    assert r == 0.0
    h = 2.0 * a
    return sum((0.5 * L / h) * _phi_cyl_dot(0.0, abs(z - zc), u, a=a, length=h) for zc in _centres(None))


def _twin_a_steady(z, a=ROD_RADIUS, length=None):
    c = a if length is None else 0.5 * length
    return 0.5 * L * sum((_g(abs(z - zc) + c, a) - _g(abs(z - zc) - c, a)) / (2.0 * c)
                         for zc in _centres(length))


def _signed():
    for k, z in enumerate(RUNGS):
        for sg in (1.0, -1.0):
            yield k, sg, sg * z


def test_lim_every_rung():
    _p(f"\n[imp8-i86] a {ROD_RADIUS / MM:g} mm, h 2a, L {L / MM:g} mm, t {T_PRINT:g} s, PERIOD {PERIOD:g} s, "
       f"D {DUTY:g}, SCALE {SCALE:.7g} K/W per unit Φ, κ {KAPPA:.7g} m^2/s")
    worst = 0.0
    for k, sg, z in _signed():
        ph, st = _twin_a(0.0, z, np.inf), _twin_a_steady(z)
        rel = abs(ph / st - 1.0)
        _p(f"[imp8-i86] (Lim) {_lab(k)} sign {sg:+g}: Φ2a(0,z,inf) quad {ph:.7g} g-form {st:.7g} "
           f"rel {rel:.3e} (band {B1_BAND:g})")
        worst = max(worst, rel)
    _p(f"[imp8-i86] (Lim) worst {worst:.3e}")
    assert worst <= B1_BAND


def test_bq_every_rung():
    worst = 0.0
    for k, sg, z in _signed():
        for t in (T_PRINT, np.inf):
            g_n, g_2n = _twin_a(0.0, z, t, n=N_GL), _twin_a(0.0, z, t, n=2 * N_GL)
            e = abs(g_n / g_2n - 1.0)
            _p(f"[imp8-i86] (Bq) {_lab(k)} sign {sg:+g} t={t:g}: GL{N_GL} {g_n:.7g} GL{2 * N_GL} {g_2n:.7g} "
               f"|n/2n-1| {e:.3e} (band {BQ_BAND:g}); printed |GL2n/quad-1| "
               f"{abs(g_2n / _twin_a(0.0, z, t) - 1.0):.3e}")
            worst = max(worst, e)
    _p(f"[imp8-i86] (Bq) worst {worst:.3e}")
    assert worst <= BQ_BAND


def test_dut_every_rung():
    t0 = time.perf_counter()
    worst = 0.0
    for k, sg, z in _signed():
        step = _duhamel(_twin_a, 0.0, z, T_STAR, DUTY)
        conv = _convolution(_twin_a_dot, 0.0, z, T_STAR, DUTY)
        e = abs(step - conv) / (SCALE * _twin_a_steady(z))
        _p(f"[imp8-i86] (Dut) {_lab(k)} sign {sg:+g}: step sum {step:.7g} conv {conv:.7g} "
           f"|diff|/θ_ss {e:.3e} (band {B1_BAND:g})")
        worst = max(worst, e)
    _p(f"[imp8-i86] (Dut) worst {worst:.3e}; {time.perf_counter() - t0:.1f} s")
    assert worst <= B1_BAND


def test_lim_a_small_radius_to_point_twin():
    small = partial(_twin_a, a=A_LIM)
    worst = 0.0
    for k, sg, z in _signed():
        reads = []
        for t in (T_PRINT, np.inf):
            reads.append((f"t={t:g}", small(0.0, z, t), _phi_tips(0.0, z, t)))
        reads.append(("D=0.5", _duhamel(small, 0.0, z, T_STAR, DUTY), _duhamel(_phi_tips, 0.0, z, T_STAR, DUTY)))
        for lab, got, ref in reads:
            rel = abs(got / ref - 1.0)
            _p(f"[imp8-i86] (Lim-a) {_lab(k)} sign {sg:+g} {lab}: a=1e-10 m {got:.7g} point twin {ref:.7g} "
               f"rel {rel:.3e} (band {B1_BAND:g})")
            worst = max(worst, rel)
    _p(f"[imp8-i86] (Lim-a) worst {worst:.3e}")
    z0 = RUNGS[0]
    ref = _phi_tips(0.0, z0, T_PRINT)
    defs = [_twin_a(0.0, z0, T_PRINT, a=ROD_RADIUS / f) / ref - 1.0 for f in (10, 20, 40)]
    orders = [np.log2(defs[i] / defs[i + 1]) for i in range(2)]
    _p(f"[imp8-i86] (Lim-a) 16 mm, t={T_PRINT:g}: Φ2a/Φ2 - 1 at a/10, /20, /40 {defs[0]:+.4e} / {defs[1]:+.4e} / "
       f"{defs[2]:+.4e}; observed order {orders[0]:.4f} / {orders[1]:.4f} (expected ≈ 2); predicted at a/10 "
       f"≈ 8e-4 (drafter)")
    _p(f"[imp8-i86] (Lim-a) own negative: |deficit| at a/10 {abs(defs[0]):.3e} asserted > B1_BAND {B1_BAND:g}")
    assert worst <= B1_BAND
    assert abs(defs[0]) > B1_BAND


def test_control_twin_cylinder_below_point_twin():
    z = RUNGS[0]
    cyl, pt = SCALE * _twin_a(0.0, z, T_PRINT), SCALE * _phi_tips(0.0, z, T_PRINT)
    ctrl = 1.0 - cyl / pt
    ss_cyl, ss_pt = SCALE * _twin_a_steady(z), SCALE * _tips_steady(0.0, z)
    d_cyl, d_pt = _duhamel(_twin_a, 0.0, z, T_STAR, DUTY), _duhamel(_phi_tips, 0.0, z, T_STAR, DUTY)
    _p(f"[imp8-i86] (ctrl) {_lab(0)} continuous θ(900 s): cylinder twin {cyl:.7g} point twin {pt:.7g} K/W; "
       f"1 - ratio {ctrl:.4e} (asserted > INTERIOR_BAND {INTERIOR_BAND:g}, < 1; predicted ≈ {PRED_CTRL:g})")
    _p(f"[imp8-i86] (ctrl) steady {ss_cyl:.7g} vs {ss_pt:.7g}: 1 - ratio {1.0 - ss_cyl / ss_pt:.4e} "
       f"(printed; predicted {PRED_CTRL_SS:g})")
    _p(f"[imp8-i86] (ctrl) D=0.5 {d_cyl:.7g} vs {d_pt:.7g}: 1 - ratio {1.0 - d_cyl / d_pt:.4e} "
       f"(below the band; not the weekly's control)")
    assert INTERIOR_BAND < ctrl < 1.0


def test_printed_record():
    _p("[imp8-i86] record: placement and time quoted from memory; W3 confirms or corrects. "
       "K per watt per rod, infinite gel, tips h = 2a centred on ±L/2.")
    ans = None
    for k, z in enumerate(RUNGS):
        cont = SCALE * _twin_a(0.0, z, T_PRINT)
        duty = _duhamel(_twin_a, 0.0, z, T_STAR, DUTY)
        ss = SCALE * _twin_a_steady(z)
        rf = _twin_a(0.0, z, T_PRINT) / _phi_tips(0.0, z, T_PRINT)
        rf_ss = _twin_a_steady(z) / _tips_steady(0.0, z)
        rf_d = duty / _duhamel(_phi_tips, 0.0, z, T_STAR, DUTY)
        tw = _twin_a(0.0, z, T_PRINT) / _phi_cyl(0.0, z, T_PRINT)
        tw_ss = _twin_a(0.0, z, np.inf) / _phi_cyl(0.0, z, np.inf)
        pt_line = _phi_tips(0.0, z, T_PRINT) / _phi(0.0, z, T_PRINT)
        sym = abs(_twin_a(0.0, z, T_PRINT) - _twin_a(0.0, -z, T_PRINT)) * SCALE
        sym_d = abs(duty - _duhamel(_twin_a, 0.0, -z, T_STAR, DUTY))
        _p(f"[imp8-i86] {_lab(k)} (+{LADDER_MM[k]:g} mm): SCALE*Φ2a 900 s {cont:.7g}, D=0.5 {duty:.7g}, "
           f"steady {ss:.7g}; |θ(+z)-θ(-z)| cont {sym:.3e} D=0.5 {sym_d:.3e}")
        _p(f"[imp8-i86]   radius factor Φ2a/Φ2: 900 s {rf:.7g}, D=0.5 {rf_d:.7g}, steady {rf_ss:.7g} "
           f"(I83 line Φ_a/Φ logged {I83_RADIUS_FACTOR[k]})")
        _p(f"[imp8-i86]   tip weighting twin_a/line_a: 900 s {tw:.7g}, steady {tw_ss:.7g} "
           f"(point twin/line here {pt_line:.7g}; I80 logged {I80_TWIN_LINE[k]}); "
           f"line steady axis form {_axis_steady(z):.7g}")
        if k == 0:
            ans = (rf, tw)
    rf, tw = ans
    dom = "the tip weighting" if abs(np.log(tw)) > abs(np.log(rf)) else "the radius"
    _p(f"[imp8-i86] answer (16 mm, continuous 900 s): radius factor Φ2a/Φ2 {rf:.7g} (moves the read "
       f"{100.0 * (1.0 - rf):.3g} %) against tip weighting twin_a/line_a × {tw:.7g}: {dom} dominates "
       f"(predicted ≈ 4 % against × 6.5)")
    _p("[imp8-i86] sensitivity (the cylinder's length is an execution detail; the weekly's to fix): "
       "h = a, cylinders [L/2 - a, L/2] and mirror, inside the rod's end")
    for k, z in enumerate(RUNGS):
        h = ROD_RADIUS
        r9 = _twin_a(0.0, z, T_PRINT, length=h) / _phi_tips(0.0, z, T_PRINT)
        rs = _twin_a(0.0, z, np.inf, length=h) / _tips_steady(0.0, z)
        rs_g = _twin_a_steady(z, length=h) / _tips_steady(0.0, z)
        _p(f"[imp8-i86]   h=a {_lab(k)}: Φ2a/Φ2 900 s {r9:.7g}, steady {rs:.7g} (g-form {rs_g:.7g})"
           + (f" (predicted ≈ {PRED_H_A[0]:g} / {PRED_H_A[1]:g})" if k == 0 else ""))
    _p(f"[imp8-i86] (aux) point-twin Φ̇ closed form at 16 mm, u=900 s: {_phi_tips_dot(0.0, RUNGS[0], 900.0):.7g}; "
       f"cylinder {_twin_a_dot(0.0, RUNGS[0], 900.0):.7g}")
    assert np.isfinite(rf) and np.isfinite(tw)
