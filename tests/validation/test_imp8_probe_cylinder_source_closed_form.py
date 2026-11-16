"""IMP-8 step 20 (§9 item 911, §10 I83): the rod's finite radius at the probe, closed form.

I55's source is the uniform DG0 cylinder r ≤ a = ``ROD_RADIUS``, |z| ≤ L/2, not the zero-radius
line I75 / I81 printed beside.  On the axis beyond the end face the disc integral of the point
kernel is exact,

    2π∫₀ᵃ erfc(√(r′² + d²)/s)/√(r′² + d²) r′ dr′ = 2πs[ierfc(d/s) − ierfc(√(a² + d²)/s)],  s = 2√(κt),

so per unit ``SCALE`` (θ = SCALE·Φ_a; Φ_a → ``_phi`` as a → 0)

    Φ_a(0, z, t) = (2s/a²) ∫ [ierfc(d/s) − ierfc(D/s)] dz′,   d = |z − z′|, D = hypot(a, d),

with the steady kernel 2/(D + d) and the steady antiderivative g(u) = u/(√(a² + u²) + u) + asinh(u/a).

Time derivative, closed: d/dt[s·ierfc(x/s)] = (s/2t)·e^{−x²/s²}/√π and
e^{−d²/s²} − e^{−D²/s²} = e^{−d²/s²}(1 − e^{−a²/s²}), so

    Φ̇_a(0, z, u) = (s²/(2a²u))·(−expm1(−a²/s²))·[erf((|z| + L/2)/s) − erf((|z| − L/2)/s)],  s² = 4κu,

whose a → 0 limit is the line's (1/2u)[erf − erf].

Disclosure (implementation): the bracket ierfc(x) − ierfc(y), y = D/s, is evaluated as the
identical ∫ₓʸ erfc(v) dv (ierfc′ = −erfc) by a 16-point Gauss–Legendre rule on the short interval
[x, y] (width ≤ a²/(2ds) ≪ 1), which removes the cancellation of the two-ierfc difference;
``_ierfc`` itself is kept and the two forms are printed side by side.

Gates (bands imported): (Lim) B1_BAND, (Bq) BQ_BAND, (Dut) B1_BAND, (K) string custody, and (G)
— flagged for the weekly to ratify — I81's gel-arm z-even read against SCALE·Φ_a at
INTERIOR_BAND.  Controls (C1) / (C2) asserted at factor 1 on INTERIOR_BAND at 16 / 17 mm.
"""

from __future__ import annotations

import gzip
import re
import time
from pathlib import Path

import numpy as np
from scipy import special
from scipy.integrate import quad

from tests.mesh.test_birdcage_implant_rod import ROD_RADIUS
from tests.validation.test_imp4_pec_cylinder_closed_form import B1_BAND
from tests.validation.test_imp8_entry48_lumped import BQ_BAND, _gl01
from tests.validation.test_imp8_probe_f2182_closed_form import DUTY, N_GL, _axis_steady, _p
from tests.validation.test_imp8_probe_ladder_closed_form import LADDER_MM, RUNGS, _lab
from tests.validation.test_imp8_rod_line_source_closed_form import KAPPA, MM, QUAD, T_PRINT, L, _phi
from tests.validation.test_imp8_rod_line_source_duty_closed_form import PERIOD, SCALE, T_STAR
from tests.validation.test_imp8_rod_tips_duty_closed_form import _convolution, _duhamel
from tests.validation.test_mat13_transient_perfusion_sphere import INTERIOR_BAND

A_EFF = 1.0022e-3  # m, I55's effective radius (20261103T182153Z_IMP-8-I55-W1s.log:56), printed only
LOG_DIR = Path(__file__).resolve().parents[2] / "docs" / "testing" / "logs"
_GL_B = _gl01(16)

#: I81's FEM reads, custody by string (asserted against the logs in ``test_k_reads_from_logs``).
#: mean / dmean: ``ladder |z|=NN mm … +/- mean theta(900 s) V K/W … D=0.5 mean V``;
#: plus / minus: ``control sign ±1 arm X: theta 16/17/18/20 …``.
FEM = {
    "gel": dict(
        log="20261109T061107Z_IMP-8-item895-W2.log",  # :124/:127/:130/:133, :134–135
        mean=("9.161894e+00", "6.662160e+00", "5.171959e+00", "3.400665e+00"),
        dmean=("1.975098e+00", "1.818740e+00", "1.665104e+00", "1.371623e+00"),
        plus=("9.220827e+00", "6.699986e+00", "5.198702e+00", "3.416442e+00"),
        minus=("9.102961e+00", "6.624334e+00", "5.145215e+00", "3.384888e+00"),
    ),
    "ti": dict(
        log="20261109T060343Z_IMP-8-item895-W1.log",  # :124/:127/:133/:145, :149–150
        mean=("1.057039e+01", "7.453618e+00", "5.672034e+00", "3.653661e+00"),
        dmean=("2.104578e+00", "1.903297e+00", "1.730378e+00", "1.421099e+00"),
        plus=("1.064117e+01", "7.500761e+00", "5.704421e+00", "3.671967e+00"),
        minus=("1.049960e+01", "7.406474e+00", "5.639646e+00", "3.635356e+00"),
    ),
    "hole": dict(
        log="20261109T061832Z_IMP-8-item895-W3.log",  # :124/:127/:130/:133, :134–135
        mean=("8.783168e+00", "6.473925e+00", "5.068436e+00", "3.361168e+00"),
        dmean=("1.891937e+00", "1.765061e+00", "1.625632e+00", "1.348305e+00"),
        plus=("8.815610e+00", "6.502218e+00", "5.090192e+00", "3.374817e+00"),
        minus=("8.750725e+00", "6.445632e+00", "5.046681e+00", "3.347518e+00"),
    ),
}


def _f(arm, key):
    return np.array([float(v) for v in FEM[arm][key]])


def _ierfc(u):
    return np.exp(-u * u) / np.sqrt(np.pi) - u * special.erfc(u)


def _bracket(d, t, a):
    """(2s/a²)[ierfc(d/s) − ierfc(D/s)], or 2/(D + d) at t = inf; vectorised in d."""
    d = np.asarray(d, dtype=float)
    big_d = np.hypot(a, d)
    if np.isinf(t):
        return 2.0 / (big_d + d)
    s = 2.0 * np.sqrt(KAPPA * t)
    x, y = d / s, big_d / s
    u, w = _GL_B
    v = x[..., None] + (y - x)[..., None] * u
    return (2.0 * s / (a * a)) * (y - x) * (special.erfc(v) @ w)


def _phi_cyl(r, z, t, a=ROD_RADIUS, length=L):
    assert r == 0.0 and abs(z) > 0.5 * length  # the axis beyond the end only, one piece
    return quad(lambda zp: float(_bracket(np.array([abs(z - zp)]), t, a)[0]),
                -0.5 * length, 0.5 * length, **QUAD)[0]


def _phi_cyl_dot(r, z, u, a=ROD_RADIUS, length=L):
    if u <= 0.0:
        return 0.0
    s2 = 4.0 * KAPPA * u
    s = np.sqrt(s2)
    za = abs(z)
    return (s2 / (2.0 * a * a * u)) * (-np.expm1(-a * a / s2)) * (
        special.erf((za + 0.5 * length) / s) - special.erf((za - 0.5 * length) / s))


def _g(u, a=ROD_RADIUS):
    return u / (np.sqrt(a * a + u * u) + u) + np.arcsinh(u / a)


def _cyl_steady(z, a=ROD_RADIUS, length=L):
    za = abs(z)
    return _g(za + 0.5 * length, a) - _g(za - 0.5 * length, a)


def _gl_cyl(z, t, n, a=ROD_RADIUS, length=L):
    u, w = _gl01(n)
    zp = -0.5 * length + length * u
    return length * float(np.dot(w, _bracket(np.abs(z - zp), t, a)))


def test_lim_every_rung():
    worst = 0.0
    _p(f"\n[imp8-i83] a {ROD_RADIUS / MM:g} mm, L {L / MM:g} mm, t {T_PRINT:g} s, PERIOD {PERIOD:g} s, "
       f"D {DUTY:g}, SCALE {SCALE:.7g} K/W per unit Φ")
    for k, z in enumerate(RUNGS):
        for sg in (1.0, -1.0):
            ph, st = _phi_cyl(0.0, sg * z, np.inf), _cyl_steady(sg * z)
            rel = abs(ph / st - 1.0)
            _p(f"[imp8-i83] (Lim) {_lab(k)} sign {sg:+g}: Φ_a(0,z,inf) {ph:.7g} g-form {st:.7g} "
               f"rel {rel:.3e} (band {B1_BAND:g})")
            worst = max(worst, rel)
    _p(f"[imp8-i83] (Lim) worst {worst:.3e}")
    assert worst <= B1_BAND


def test_bq_every_rung():
    worst = 0.0
    for k, z in enumerate(RUNGS):
        for sg in (1.0, -1.0):
            for t in (T_PRINT, np.inf):
                g_n, g_2n = _gl_cyl(sg * z, t, N_GL), _gl_cyl(sg * z, t, 2 * N_GL)
                e = abs(g_n / g_2n - 1.0)
                _p(f"[imp8-i83] (Bq) {_lab(k)} sign {sg:+g} t={t:g}: GL{N_GL} {g_n:.7g} "
                   f"GL{2 * N_GL} {g_2n:.7g} |n/2n-1| {e:.3e} (band {BQ_BAND:g}); printed "
                   f"|GL2n/quad-1| {abs(g_2n / _phi_cyl(0.0, sg * z, t) - 1.0):.3e}")
                worst = max(worst, e)
    _p(f"[imp8-i83] (Bq) worst {worst:.3e}")
    assert worst <= BQ_BAND


def test_dut_every_rung():
    t0 = time.perf_counter()
    worst = 0.0
    for k, z in enumerate(RUNGS):
        for sg in (1.0, -1.0):
            step = _duhamel(_phi_cyl, 0.0, sg * z, T_STAR, DUTY)
            conv = _convolution(_phi_cyl_dot, 0.0, sg * z, T_STAR, DUTY)
            e = abs(step - conv) / (SCALE * _cyl_steady(sg * z))
            _p(f"[imp8-i83] (Dut) {_lab(k)} sign {sg:+g}: step sum {step:.7g} conv {conv:.7g} "
               f"|diff|/θ_ss {e:.3e} (band {B1_BAND:g})")
            worst = max(worst, e)
    _p(f"[imp8-i83] (Dut) worst {worst:.3e}; {time.perf_counter() - t0:.1f} s")
    assert worst <= B1_BAND


def _read_log(name):
    p = LOG_DIR / name
    if p.exists():
        return p.read_text().splitlines()
    with gzip.open(str(p) + ".gz", "rt") as fh:
        return fh.read().splitlines()


def test_k_reads_from_logs():
    for arm, rec in FEM.items():
        lines = _read_log(rec["log"])
        means = [m for m in (re.search(
            r"\[imp8-i81\] ladder \|z\|=(\d+) mm.*? arm " + arm
            + r": \+/- mean theta\(900 s\) (\S+) K/W.*?D=0\.5 mean (\S+),", ln) for ln in lines) if m]
        ctrls = [m for m in (re.search(
            r"\[imp8-i81\] control sign ([+-]1) arm " + arm
            + r": theta 16/17/18/20 (\S+) / (\S+) / (\S+) / (\S+);", ln) for ln in lines) if m]
        _p(f"[imp8-i83] (K) arm {arm}: {rec['log']} mean lines {len(means)}, control lines {len(ctrls)}")
        assert len(means) == 4 and len(ctrls) == 2
        assert tuple(int(m.group(1)) for m in means) == tuple(int(round(z / MM)) for z in RUNGS)
        assert tuple(m.group(2) for m in means) == rec["mean"]
        assert tuple(m.group(3) for m in means) == rec["dmean"]
        by_sign = {m.group(1): tuple(m.group(i) for i in range(2, 6)) for m in ctrls}
        assert by_sign["+1"] == rec["plus"] and by_sign["-1"] == rec["minus"]


def _closed(z):
    return SCALE * _phi_cyl(0.0, z, T_PRINT), _duhamel(_phi_cyl, 0.0, z, T_STAR, DUTY)


def test_g_gel_arm_z_even_against_cylinder():
    """(G), flagged for the weekly's ratification: the gel arm's ± mean vs SCALE·Φ_a, both duties."""
    mean, dmean = _f("gel", "mean"), _f("gel", "dmean")
    plus, minus = _f("gel", "plus"), _f("gel", "minus")
    worst = 0.0
    for k, z in enumerate(RUNGS):
        for sg in (1.0, -1.0):
            c, d = _closed(sg * z)
            rc, rd = mean[k] / c - 1.0, dmean[k] / d - 1.0
            _p(f"[imp8-i83] (G) {_lab(k)} sign {sg:+g}: gel mean {mean[k]:.7g} / SCALE*Φ_a {c:.7g} - 1 = "
               f"{rc:+.3e}; D=0.5 mean {dmean[k]:.7g} / _duhamel {d:.7g} - 1 = {rd:+.3e} "
               f"(band {INTERIOR_BAND:g})")
            worst = max(worst, abs(rc), abs(rd))
        c, _ = _closed(z)
        _p(f"[imp8-i83] (G) printed only, per sign (not gated; I81's mesh asymmetry): "
           f"+z {plus[k] / c - 1.0:+.3e}, -z {minus[k] / c - 1.0:+.3e}")
    _p(f"[imp8-i83] (G) worst {worst:.3e} (band {INTERIOR_BAND:g})")
    assert worst <= INTERIOR_BAND


def test_controls_c1_c2():
    gel = _f("gel", "mean")
    for k, z in enumerate(RUNGS):
        c, _ = _closed(z)
        c1 = 1.0 - gel[k] / (SCALE * _phi(0.0, z, T_PRINT))
        c1_exact = 1.0 - _phi_cyl(0.0, z, T_PRINT) / _phi(0.0, z, T_PRINT)
        c2 = {arm: abs(_f(arm, "mean")[k] / c - 1.0) for arm in ("ti", "hole")}
        asserted = k < 2
        _p(f"[imp8-i83] (C1) {_lab(k)}: 1 - gel/(SCALE*Φ_line) {c1:.6f} (exact-FEM ceiling "
           f"1 - Φ_a/Φ {c1_exact:.5f}); (C2) |ti/(SCALE*Φ_a)-1| {c2['ti']:.5f}, |hole/..-1| "
           f"{c2['hole']:.5f}; {'asserted > band' if asserted else 'printed only'} {INTERIOR_BAND:g}")
        if asserted:
            assert c1 > INTERIOR_BAND
            assert c2["ti"] > INTERIOR_BAND and c2["hole"] > INTERIOR_BAND


def _drop(v):
    return 1.0 - v[1] / v[0]


def test_printed_record():
    _p("[imp8-i83] record: placement and time quoted from memory; W3 confirms or corrects. "
       "K per watt, 1 W in the rod, infinite gel.")
    for k, z in enumerate(RUNGS):
        c, d = _closed(z)
        ss = SCALE * _cyl_steady(z)
        rf = _phi_cyl(0.0, z, T_PRINT) / _phi(0.0, z, T_PRINT)
        c_eff = SCALE * _phi_cyl(0.0, z, T_PRINT, a=A_EFF)
        rf_eff = c_eff / (SCALE * _phi(0.0, z, T_PRINT))
        _p(f"[imp8-i83] {_lab(k)} (+{LADDER_MM[k]:g} mm): SCALE*Φ_a 900 s {c:.7g}, D=0.5 {d:.7g}, "
           f"steady {ss:.7g}; Φ_a/Φ {rf:.7g}; at a=1.0022 mm {c_eff:.7g} (Φ_a/Φ {rf_eff:.7g})")
        for arm in ("gel", "ti", "hole"):
            m, pl, mi = _f(arm, "mean")[k], _f(arm, "plus")[k], _f(arm, "minus")[k]
            _p(f"[imp8-i83]   arm {arm}: mean/Φ_a {m / c:.6f}, +z/Φ_a {pl / c:.6f}, -z/Φ_a {mi / c:.6f}")
        g = _f("gel", "mean")[k]
        _p(f"[imp8-i83]   rod effect on the gated baseline: ti/gel - 1 {_f('ti', 'mean')[k] / g - 1.0:+.4f}, "
           f"hole/gel - 1 {_f('hole', 'mean')[k] / g - 1.0:+.4f}")
    line = [SCALE * _phi(0.0, z, T_PRINT) for z in RUNGS[:2]]
    cyl = [_closed(z)[0] for z in RUNGS[:2]]
    dl, dc = _drop(line), _drop(cyl)
    dr = {arm: _drop(_f(arm, "mean")) for arm in ("gel", "ti", "hole")}
    _p(f"[imp8-i83] first drop (16 -> 17 mm), means: line {dl:.6f}, cylinder {dc:.6f}, gel {dr['gel']:.6f}, "
       f"ti {dr['ti']:.6f}, hole {dr['hole']:.6f}")
    _p(f"[imp8-i83] decomposition: radius (cyl - line) {dc - dl:+.6f}; mesh (gel - cyl) {dr['gel'] - dc:+.6f}; "
       f"rod (ti - gel) {dr['ti'] - dr['gel']:+.6f}; ti - line {dr['ti'] - dl:+.6f}")
    for arm in ("gel", "ti", "hole"):
        dp, dm = _drop(_f(arm, "plus")), _drop(_f(arm, "minus"))
        _p(f"[imp8-i83] per sign arm {arm}: +z drop {dp:.6f}, -z drop {dm:.6f}, spread {abs(dp - dm):.6f}; "
           f"+z - cyl {dp - dc:+.6f}, -z - cyl {dm - dc:+.6f}")
    z0 = RUNGS[0]
    defs = [_phi_cyl(0.0, z0, T_PRINT, a=ROD_RADIUS / f) / _phi(0.0, z0, T_PRINT) - 1.0 for f in (10, 20, 40)]
    orders = [np.log2(defs[i] / defs[i + 1]) for i in range(2)]
    _p(f"[imp8-i83] a -> 0 (16 mm, printed only): Φ_a/Φ - 1 at a/10, /20, /40 "
       f"{defs[0]:.4e} / {defs[1]:.4e} / {defs[2]:.4e}; observed order {orders[0]:.4f} / {orders[1]:.4f} "
       f"(expected ≈ 2)")
    x, y = 0.3, 0.3 + 1e-3
    _p(f"[imp8-i83] bracket forms at x={x:g}, y={y:g}: ierfc(x)-ierfc(y) {_ierfc(x) - _ierfc(y):.15e} vs "
       f"GL16 ∫erfc {(y - x) * float(special.erfc(x + (y - x) * _GL_B[0]) @ _GL_B[1]):.15e}")
