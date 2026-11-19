r"""`MAT-21` step 48 — the dose is not per-watt, closed form (§9 item 962, §10 P88; brief
`docs/planning/chunks/MAT-21.md` "§9 item 962").

Numpy (complex build only because P83 imports P79 → P61 / C6 and C33's reader).  Sets no band
and moves none; every band is imported.  P81's end-weighted closed-form history at c\* per
deposited gel watt (``p83._data()["sd"][side]["dep"]["n"]``, one ``SL`` batch, both c\* equal
by P83's (Sym)) is scaled to P_dep ∈ {0.1, 1, 10, 100} W — θ = P·h by the series' homogeneity
— and each P gets its own ``MAT-20`` ``cem43`` call on T_INF + P·h.  ``cem43`` is nonlinear in
θ, so the excess over the unheated 20 °C baseline D0 is not proportional to P.

Asserted (collected, then asserted — item 668): preconditions (grid, t\*, row 0, labels, the
t\* token, (Sym), (R) every sample below 43 °C so only the R = 0.25 branch is in play); (D0);
(P83-tie) by string at P83's own normaliser; (Lim) ``_q_check_g``'s independent adaptive
quadrature at ``B1_BAND`` with its R-swap miss > ``DOSE_BAR``; (Cem) ``p80._cem_bracket``;
controls (N-0) the zero-power dose is D0 and (N-lin) excess(10 W)/excess(1 W) differs from 10
by more than ``INTERIOR_BAND``.  Printed (🟡 record, no band): per P the max θ / max T, cem43,
excess, excess/D0, excess per watt, the linear term on the same interpolant and the
nonlinearity; successive ratios; the per-watt excess at 100 W over that at 0.1 W.

Not claimed: a dose band, a damage threshold, a compliance figure, a tissue model; which
picture the clinical heating follows, or any entry-59 / 60 figure; the R = 0.5 branch (never
reached); that 100 W in the gel box is a realistic drive.
"""

from __future__ import annotations

import math

import numpy as np

from fem_em_solver.thermal.dose import cem43
from tests.validation import test_mat21_eddy_source_probe_points_closed_form as p80
from tests.validation import test_mat21_field_strength_heating_closed_form as p83
from tests.validation import test_mat21_wall_uniform_source_probe_points_closed_form as p75
from tests.validation.test_imp4_pec_cylinder_closed_form import B1_BAND
from tests.validation.test_mat13_transient_perfusion_sphere import INTERIOR_BAND
from tests.validation.test_mat20_tuned_heating_dose import R_BELOW, T_REF
from tests.validation.test_mat21_box_closed_form import T_INF
from tests.validation.test_mat21_box_dose_closed_form import DOSE_BAR, LN_R, T_GRID
from tests.validation.test_mat21_box_dt_duty125_dose_closed_form import _q_check_g

PFX = "MAT21-48"
P_SET = (0.1, 1.0, 10.0, 100.0)  # W deposited in the gel (the step's)
# `20261020T190033Z_MAT-21-step15-W1-probe-dose-closed-form.log.gz:83`
D0_CITED = 4.009592657894245e-11
H_TSTAR = "0.00356031"  # `…P83-W1.log:81`, `…P81-W1.log:95` (end picture, K per dep. W at t*)
# P83's `.7g` cem43 and `.6e` excess of the end picture per accepted watt (`…P83-W1.log:88`, `:92`)
P83_DOSE = {"1.5T": ("4.010097e-11", "5.045885e-15"), "3T": ("4.011457e-11", "1.864127e-14")}
LABELS = {"1.5T": "x-105mm dz-50mm", "3T": "x+105mm dz+50mm"}
N0_BAND = 1e-12  # the step's number (= B1_BAND)

# Predictions (brief, printed only; a miss outside a range is disclosed, never red)
PRED_TABLE = {
    0.1: (2.3813e-04, 9.548e-14, 1.62e-04),
    1.0: (2.3848e-03, 9.562e-14, 1.62e-03),
    10.0: (2.4199e-02, 9.703e-14, 1.64e-02),
    100.0: (0.2818, 1.130e-13, 0.184),
}
PRED_NLIN = (0.147, 0.135, 0.160)
PRED_SUCC = (10.015, 10.147, 11.65)
PRED_PW_RATIO = 1.183
PRED_CEM100 = 5.14e-11
PRED_TMAX = 20.356031

_CACHE = {}


def _say(msg):
    print(msg, flush=True)


def _trap(y, t):
    return float(np.sum(0.5 * (y[1:] + y[:-1]) * np.diff(t)))


def _data():
    if _CACHE:
        return _CACHE
    d = p83._data()
    sd = d["sd"]
    h = sd["1.5T"]["dep"]["n"]
    h3 = sd["3T"]["dep"]["n"]
    d0 = cem43(T_GRID, np.full(T_GRID.size, T_INF))
    lin_unit = math.log(1.0 / R_BELOW) * R_BELOW ** (T_REF - T_INF) * _trap(h, T_GRID)
    rows = {}
    for P in P_SET:
        th = h * P
        cem = cem43(T_GRID, T_INF + th)
        exc = cem - d0
        rows[P] = {
            "th": th,
            "cem": cem,
            "exc": exc,
            "lin": P * lin_unit,
            "tmax": float(np.max(T_INF + th)),
        }
    _CACHE.update(d=d, sd=sd, h=h, h3=h3, d0=d0, rows=rows, lin_unit=lin_unit)
    return _CACHE


def _fail(gates):
    return [g for g in gates if not g[1]]


# ------------------------------------------------------------------ tests
def test_preconditions_d0_tie():
    gates = []
    c = _data()
    sd, h, h3 = c["sd"], c["h"], c["h3"]
    gates.append((("T_GRID[0]", float(T_GRID[0])), T_GRID[0] == 0.0, 0))
    dts = abs(p75.T_STAR - float(T_GRID[-1]))
    gates.append((("T_STAR", dts), dts <= 1e-9, dts))
    gates.append((("row0", float(h[0]), float(h3[0])), h[0] == 0.0 and h3[0] == 0.0, 0))
    gates.append((("rows", h.size), h.size == T_GRID.size == 41, h.size))
    for s in p83.SIDES:
        gates.append((("label", s, sd[s]["label"]), sd[s]["label"] == LABELS[s], 0))
    tok = f"{h[-1]:.7g}"
    gates.append((("H_TSTAR", tok), tok == H_TSTAR, 0))
    sym = float(np.max(np.abs(h3 - h)) / np.max(h))
    gates.append((("Sym", sym), sym <= B1_BAND, sym))
    tmax = max(r["tmax"] for r in c["rows"].values())
    gates.append((("R", tmax), tmax < T_REF, tmax))
    lnr_ok = LN_R == math.log(R_BELOW)
    gates.append((("LN_R", LN_R), lnr_ok, 0))
    _say(
        f"{PFX} T_GRID[0] {float(T_GRID[0])} (== 0 ASSERTED); |p75.T_STAR - T_GRID[-1]| {dts:.3e} s "
        f"(<= 1e-9 ASSERTED); {h.size} rows, row 0 of both histories 0 exactly (ASSERTED)"
    )
    _say(
        f"{PFX} c* 1.5 T {sd['1.5T']['label']!r}, 3 T {sd['3T']['label']!r} (ASSERTED); h(t*) "
        f"{tok} K per deposited W == {H_TSTAR} (string ASSERTED)"
    )
    _say(f"{PFX} (Sym) max|h3 - h| / max h {sym:.3e} (band {B1_BAND:g}, ASSERTED; P83 read 0)")
    _say(
        f"{PFX} (R) hottest sample over P_SET T_INF + P h {tmax:.6f} C < T_REF {T_REF:g} C "
        f"(ASSERTED; predicted {PRED_TMAX}); only the R = {R_BELOW} branch is reached, "
        f"_quad_dose LN_R == ln {R_BELOW} {lnr_ok} (ASSERTED); the R = 0.5 branch is stated, not "
        f"exercised"
    )
    d0 = c["d0"]
    rel = abs(d0 - D0_CITED) / D0_CITED
    gates.append((("D0", rel), rel <= B1_BAND, rel))
    _say(
        f"{PFX} (D0) cem43(T_GRID, T_INF) {d0:.7g} s vs cited (step 15 W1 :83): rel {rel:.3e} "
        f"(band {B1_BAND:g}, ASSERTED; P83 read 1.11e-16)"
    )
    for s in p83.SIDES:
        P = sd[s]["s"]
        hs = sd[s]["dep"]["n"]
        cem = cem43(T_GRID, T_INF + hs * P)
        got = (f"{cem:.7g}", f"{cem - d0:.6e}")
        gates.append((("P83-tie", s, got, P83_DOSE[s]), got == P83_DOSE[s], 0))
        _say(
            f"{PFX} (P83-tie) {s} at P_dep = s {P:.7g} W: cem43 {got[0]} s, excess {got[1]} s vs "
            f"P83 {P83_DOSE[s]} equal {got == P83_DOSE[s]} (strings ASSERTED)"
        )
    bad = _fail(gates)
    assert not bad, bad


def test_lim_cem():
    gates = []
    c = _data()
    worst_q, least_swap, least_m = 0.0, np.inf, np.inf
    for P in P_SET:
        th = c["rows"][P]["th"]
        _dose, _ref, rel, _swap, miss = _q_check_g(T_INF + th, T_GRID)
        gates.append((("Lim", P, rel), rel <= B1_BAND, rel))
        gates.append((("swap", P, miss), miss > DOSE_BAR, miss))
        cem, lo, hi = p80._cem_bracket(th)
        gates.append((("Cem", P, lo, cem, hi), lo <= cem <= hi, cem))
        m = min(cem - lo, hi - cem) / cem
        worst_q, least_swap, least_m = max(worst_q, rel), min(least_swap, miss), min(least_m, m)
        _say(
            f"{PFX} (Lim) P {P:g} W: cem43 vs _quad_dose rel {rel:.3e} (band {B1_BAND:g}, "
            f"ASSERTED); R-swap miss {miss:.4e} (> {DOSE_BAR:g} ASSERTED); (Cem) bracket "
            f"[{lo:.6e}, {hi:.6e}] holds {lo <= cem <= hi} (ASSERTED), margin {m:.3e}"
        )
    _say(
        f"{PFX} (Lim) worst rel {worst_q:.3e}; R-swap least miss {least_swap:.4e}; (Cem) least "
        f"relative margin {least_m:.3e}"
    )
    bad = _fail(gates)
    assert not bad, bad


def test_controls():
    gates = []
    c = _data()
    h = c["h"]
    z = cem43(T_GRID, T_INF + 0.0 * h)
    n0 = abs(z - D0_CITED) / D0_CITED
    gates.append((("N-0", n0), n0 <= N0_BAND, n0))
    _say(
        f"{PFX} (N-0) P_dep = 0: |cem43 - D0_cited| / D0_cited {n0:.3e} (band {N0_BAND:g}, "
        f"ASSERTED; predicted 1.1e-16)"
    )
    r = c["rows"][10.0]["exc"] / c["rows"][1.0]["exc"]
    dev = abs(r - 10.0)
    gates.append((("N-lin", r, dev), dev > INTERIOR_BAND, dev))
    pv, lo, hi = PRED_NLIN
    _say(
        f"{PFX} (N-lin) r = excess(10 W) / excess(1 W) {r:.6f}; |r - 10| {dev:.4f} (ASSERTED > "
        f"{INTERIOR_BAND:g}; predicted {pv}, range [{lo}, {hi}] inside {lo <= dev <= hi}); "
        f"separation {dev / INTERIOR_BAND:.1f}x; relative form |r/10 - 1| {abs(r / 10 - 1):.3e} "
        f"(printed, not asserted)"
    )
    bad = _fail(gates)
    assert not bad, bad


def test_record_printed():
    """The record (no band); asserts only that it was formed."""
    c = _data()
    d0, rows = c["d0"], c["rows"]
    _say(
        f"{PFX} P81's end-weighted picture at c* (both c* one history), per deposited gel watt; "
        f"D0 {d0:.7g} s; linear term per W {c['lin_unit']:.6e} s/W (ln 4 R^(T_REF-T_INF) "
        f"trapezoid(h), exact on the interpolant)"
    )
    _say(
        f"{PFX} TABLE P_dep W | max theta K | max T C | cem43 s | D0 s | excess s | excess/D0 | "
        f"excess per W s/W | excess_lin s | excess/excess_lin - 1 | predicted (exc/D0, per W, "
        f"nonlin) inside"
    )
    misses = []
    for P in P_SET:
        rw = rows[P]
        e = rw["exc"]
        vals = (e / d0, e / P, e / rw["lin"] - 1.0)
        pr = PRED_TABLE[P]
        tol = (0.03 if P == 100.0 else 0.01, 0.03 if P == 100.0 else 0.01, 0.10)
        ins = tuple(abs(v / p - 1.0) <= t for v, p, t in zip(vals, pr, tol))
        if not all(ins):
            misses.append((P, vals, pr))
        _say(
            f"{PFX} TABLE {P:g} | {float(np.max(rw['th'])):.6e} | {rw['tmax']:.6f} | "
            f"{rw['cem']:.6e} | {d0:.6e} | {e:.6e} | {vals[0]:.5e} | {vals[1]:.5e} | "
            f"{rw['lin']:.6e} | {vals[2]:.4e} | {pr} {ins}"
        )
    succ = [rows[b]["exc"] / rows[a]["exc"] for a, b in zip(P_SET[:-1], P_SET[1:])]
    _say(
        f"{PFX} successive ratios excess(P_k+1)/excess(P_k): "
        f"{' / '.join(f'{v:.5g}' for v in succ)} (predicted approx "
        f"{' / '.join(str(v) for v in PRED_SUCC)})"
    )
    pw = (rows[100.0]["exc"] / 100.0) / (rows[0.1]["exc"] / 0.1)
    _say(
        f"{PFX} per-watt excess at 100 W / at 0.1 W {pw:.5g} (predicted approx {PRED_PW_RATIO}): "
        f"one number per watt cannot carry a quantity whose per-watt value moves "
        f"{100 * (pw - 1):.1f} % over the range"
    )
    _say(
        f"{PFX} cem43 at 100 W {rows[100.0]['cem']:.4g} s (predicted approx {PRED_CEM100:g}); "
        f"max T {rows[100.0]['tmax']:.6f} C, the R = 0.5 branch (T >= {T_REF:g} C) not reached"
    )
    _say(
        f"{PFX} NOTE a cem43 per accepted watt (P76, P83, P86, C39) is the 1 W history's dose: "
        f"approx D0 plus a near-linear excess, never scaled by s or by P"
    )
    _say(
        f"{PFX} predictions outside their ranges (disclosed, never red): "
        f"{misses if misses else 'none'}"
    )
    assert len(rows) == len(P_SET) and all(np.isfinite(r["exc"]) for r in rows.values())
