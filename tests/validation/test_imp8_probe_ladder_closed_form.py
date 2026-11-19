"""`IMP-8` step 18 (§9 item 894, §10 I80): the lumped probe-placement ladder,
closed form, numpy / scipy only.

Brief: docs/planning/chunks/IMP-8.md, "§9 item 894 — 2026-11-08".  I78's two
idealised source shapes (the zero-radius line Φ and the tip pair Φ₂, 1 W per rod
through ``SCALE``) in an infinite gel medium, read on the axis at four placements
z_k = L/2 + d_k, d_k ∈ ``LADDER_MM`` = (1, 2, 3, 5) mm, both signs, at
t = 900 s, continuous and at D = 0.5 on ``PERIOD``.  Rung 1 (16 mm) is I75's
bracket probe (1 mm beyond a flat end); rung 2 (17 mm) is entry 63's clinical
probe (1 mm beyond the capsule apex, 2 mm beyond the line's end).

Asserted at every rung (bands imported, unmoved):
- (Lim) |Φ(0, z_k, ∞)/``_axis_steady``(z_k) − 1| ≤ ``B1_BAND``.
- (Bq) |GL_N/GL_2N − 1| ≤ ``BQ_BAND`` (I78's ``_gl_line``, N = ``N_GL``) at
  t = ``T_PRINT`` and ∞.
- (Dut) line and twin, |step sum − convolution| / θ_ss(z_k) ≤ ``B1_BAND``.
- (Rep) at rungs 1 and 2, I78's ``_reads`` reproduces its printed record
  (``20261108T133157Z_IMP-8-item874-W1.log:58–62``) by ``.7g`` string.
- Control (I78's, extended; factor 1 on ``INTERIOR_BAND``): continuous θ(900 s)
  and θ_ss strictly decrease rung to rung, each drop in (``INTERIOR_BAND``, 1).
  The D = 0.5 reads are printed only (the heat spreads outward after switch-off).

Printed (🟡 record, no band, K per watt per rod): placement and time quoted from
memory; W3 confirms or corrects.  The rungs are placements, not a convergence
ladder; no rate is read.  Not claimed: box walls, the rod's own conduction
(I75 / I81), a ball-averaged read, the capsule's caps as a source, any F2182
figure or published value.  The twin is printed beside the line, never compared.
"""

from __future__ import annotations

import numpy as np

from tests.validation.test_imp4_pec_cylinder_closed_form import B1_BAND
from tests.validation.test_imp8_entry48_lumped import BQ_BAND
from tests.validation.test_imp8_probe_f2182_closed_form import (
    DUTY,
    N_GL,
    SHAPES,
    Z_C,
    Z_P,
    _axis_steady,
    _gl_line,
    _p,
    _reads,
)
from tests.validation.test_imp8_rod_line_source_closed_form import MM, T_PRINT, L, _phi
from tests.validation.test_imp8_rod_line_source_duty_closed_form import SCALE, T_STAR
from tests.validation.test_imp8_rod_tips_duty_closed_form import _convolution, _duhamel
from tests.validation.test_mat13_transient_perfusion_sphere import INTERIOR_BAND

LADDER_MM = (1.0, 2.0, 3.0, 5.0)  # mm beyond L/2
RUNGS = tuple(0.5 * L + d * MM for d in LADDER_MM)  # 16, 17, 18, 20 mm
LABELS = {
    0: "the bracket's probe (I75: 1 mm beyond a flat end)",
    1: "the clinical probe (entry 63: 1 mm beyond the capsule apex, 2 mm beyond the line's end)",
}
#: I78's printed record (``20261108T133157Z_IMP-8-item874-W1.log:58–62``): cont / D=0.5 / θ_ss.
I78_RECORD = {
    ("line", 0): ("9.648662", "1.978384", "15.18156"),
    ("line", 1): ("6.810428", "1.822237", "12.25753"),
    ("twin", 0): ("63.14", "2.498341", "68.45374"),
    ("twin", 1): ("29.97039", "2.451249", "35.22961"),
}
#: The brief's steady-form drop predictions (unlogged review arithmetic), printed beside, never asserted.
CTRL_PRED = {"line": (0.1926, 0.1352, 0.1885), "twin": (0.4854, 0.3155, 0.3714)}


def _lab(k):
    return f"rung {k + 1} (z {RUNGS[k] / MM:g} mm" + (f", {LABELS[k]})" if k in LABELS else ")")


def test_rungs_match_i78_placements():
    # Rungs 1 and 2 are I78's Z_P / Z_C by the same expression, bit for bit.
    assert RUNGS[0] == Z_P
    assert RUNGS[1] == Z_C


def test_lim_every_rung():
    worst = 0.0
    for k, z in enumerate(RUNGS):
        phi_inf = _phi(0.0, z, np.inf)
        ax = _axis_steady(z)
        rel = abs(phi_inf / ax - 1.0)
        _p(
            f"[imp8-18] (Lim) {_lab(k)}: Φ(0,z,inf) {phi_inf:.7g} axis steady {ax:.7g} rel {rel:.3e} "
            f"(band {B1_BAND:g})"
        )
        worst = max(worst, rel)
    _p(f"[imp8-18] (Lim) worst {worst:.3e}")
    assert worst <= B1_BAND


def test_bq_every_rung():
    worst = 0.0
    for k, z in enumerate(RUNGS):
        for t in (T_PRINT, np.inf):
            g_n, g_2n = _gl_line(z, t, N_GL), _gl_line(z, t, 2 * N_GL)
            e = abs(g_n / g_2n - 1.0)
            _p(
                f"[imp8-18] (Bq) {_lab(k)} t={t:g}: GL{N_GL} {g_n:.7g} GL{2 * N_GL} {g_2n:.7g} "
                f"|n/2n-1| {e:.3e} (band {BQ_BAND:g}); printed |GL2n/_phi-1| "
                f"{abs(g_2n / _phi(0.0, z, t) - 1.0):.3e}"
            )
            worst = max(worst, e)
    _p(f"[imp8-18] (Bq) worst {worst:.3e}")
    assert worst <= BQ_BAND


def test_dut_every_rung():
    worst = 0.0
    for k, z in enumerate(RUNGS):
        for name, (kern, kern_dot, steady) in SHAPES.items():
            step = _duhamel(kern, 0.0, z, T_STAR, DUTY)
            conv = _convolution(kern_dot, 0.0, z, T_STAR, DUTY)
            e = abs(step - conv) / (SCALE * steady(z))
            _p(
                f"[imp8-18] (Dut) {_lab(k)} {name}: step sum {step:.7g} conv {conv:.7g} "
                f"|diff|/θ_ss {e:.3e} (band {B1_BAND:g})"
            )
            worst = max(worst, e)
    _p(f"[imp8-18] (Dut) worst {worst:.3e}")
    assert worst <= B1_BAND


def test_rep_i78_record_rungs_1_2():
    misses = 0
    for (name, k), want in I78_RECORD.items():
        kern, _kd, steady = SHAPES[name]
        got = tuple(f"{x:.7g}" for x in _reads(kern, steady, RUNGS[k]))
        ok = got == want
        misses += 0 if ok else 1
        _p(
            f"[imp8-18] (Rep) {name} {_lab(k)}: got {got} want {want} {'match' if ok else 'MISMATCH'}"
        )
    _p(
        f"[imp8-18] (Rep) {len(I78_RECORD) * 3 - 3 * misses}/{len(I78_RECORD) * 3} strings in matching rows; "
        f"{misses} row mismatches"
    )
    assert misses == 0


def test_control_ladder_monotone():
    for name, (kern, _kd, steady) in SHAPES.items():
        cont = [SCALE * kern(0.0, z, T_STAR) for z in RUNGS]
        ss = [SCALE * steady(z) for z in RUNGS]
        for k in range(len(RUNGS) - 1):
            d_c = 1.0 - cont[k + 1] / cont[k]
            d_s = 1.0 - ss[k + 1] / ss[k]
            _p(
                f"[imp8-18] (ctrl) {name} rung {k + 1}->{k + 2} ({RUNGS[k] / MM:g}->{RUNGS[k + 1] / MM:g} mm): "
                f"drop θ(900 s) {d_c:.7g}  drop θ_ss {d_s:.7g} (review's steady prediction "
                f"{CTRL_PRED[name][k]:.4f}; "
                f"band > {INTERIOR_BAND:g}, < 1)"
            )
            assert cont[k + 1] < cont[k]
            assert ss[k + 1] < ss[k]
            assert INTERIOR_BAND < d_c < 1.0
            assert INTERIOR_BAND < d_s < 1.0


def test_printed_record():
    _p(
        "[imp8-18] the lumped probe-placement ladder — placement and time quoted from memory; W3 confirms "
        "or corrects — K per watt per rod, infinite gel"
    )
    rec = {}
    for name, (kern, _kd, steady) in SHAPES.items():
        for k, z in enumerate(RUNGS):
            rec[(name, k)] = _reads(kern, steady, z)
    for k, z in enumerate(RUNGS):
        for name, (kern, _kd, _st) in SHAPES.items():
            cont, duty, ss = rec[(name, k)]
            c0, d0, s0 = rec[(name, 0)]
            sym = abs(SCALE * kern(0.0, z, T_STAR) - SCALE * kern(0.0, -z, T_STAR))
            sym_d = abs(
                _duhamel(kern, 0.0, z, T_STAR, DUTY) - _duhamel(kern, 0.0, -z, T_STAR, DUTY)
            )
            _p(
                f"[imp8-18] {name} {_lab(k)}: θ(900 s) cont {cont:.7g}  D=0.5 {duty:.7g}  θ_ss {ss:.7g}  "
                f"duty/cont {duty / cont:.7g}; /θ(16 mm) cont {cont / c0:.7g} D=0.5 {duty / d0:.7g} "
                f"ss {ss / s0:.7g}; |θ(+z)-θ(-z)| cont {sym:.3e} D=0.5 {sym_d:.3e}"
            )
        tw, ln = rec[("twin", k)], rec[("line", k)]
        _p(
            f"[imp8-18] twin/line {_lab(k)}: cont {tw[0] / ln[0]:.7g}  D=0.5 {tw[1] / ln[1]:.7g}  "
            f"ss {tw[2] / ln[2]:.7g}"
        )
    for v in rec.values():
        assert all(np.isfinite(x) and x > 0.0 for x in v)
