"""`IMP-7` step 10 — the dose sensitivity bound on the hollow-sphere duty series, a 3 × 2
matrix (§9 item 711, §10 lane 3 → chain I, I49; brief `docs/planning/chunks/IMP-7.md` "§9
item 711").

Numpy (the imported fixture modules bring dolfinx in).  Sets no band and moves none; every
band is imported.  Histories: `pec_sphere_implant_duty_cycle_from_modes` on ``T_GRID`` at
D = 1, 0.5, 0.25 and period ``PERIOD``, at ``P_POLE`` and ``P_MID``.  Equal mean power is
θ / D, a scale, so the ``SCALES`` labels take the factor 1 / D per duty here: the imported
tuple's factor is bound to D = 0.5.

P39′'s (Bd) relates two histories and this step has one, so the pair is the ``N_DEFAULT``
series and the same series at ``N_DEFAULT // 2`` modes, a real, known perturbation of the
comparand (the 2026-11-01 03:00 review's execution detail).  (Bd)(ii) has no FEM history
here: ``INTERIOR_BAND``·max θ_s is printed per cell, the cap I50 reads against.

Asserted, per cell and scale: (Q-m) `cem43` against `quad` ≤ ``DOSE_BAR``, max T < 43 °C;
(Win-m) the six-minute window at ``AVG_BAR``; (Bd-i) m_p ≤ ln 4·M_p, and (Bd-i′) the same on
``N_COARSE``, where the pair differs by more than rounding.  (C100-m) D = 1 equals the
continuous history to ``B1_BAND``.  Negative controls: (Rc-m) the R swap on every cell,
(β0-m) β = 0 at the pole, (Mis-m) D = 0.25 against D = 0.5, each missing > ``DOSE_BAR``.

Not claimed: a FEM history (I50), (Bd)(ii), a rate, a thermal-damage figure.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from fem_em_solver.thermal.dose import cem43
from fem_em_solver.utils.analytical import (
    pec_sphere_implant_duty_cycle_from_modes,
    pec_sphere_implant_transient_modes,
)
from tests.validation.test_imp4_pec_cylinder_closed_form import B1_BAND
from tests.validation.test_imp7_duty_cycle_fem import PERIOD, A, B
from tests.validation.test_imp7_duty_dose_closed_form import (
    AVG_BAR,
    DOSE_BAR,
    DT,
    POLE,
    PTS,
    SCALES,
    T_END,
    T_GRID,
)
from tests.validation.test_imp7_implant_heating_closed_form import FIX, T_INF
from tests.validation.test_imp7_transient_closed_form import N_DEFAULT, T_STAR, _modes, _theta
from tests.validation.test_mat12_transient_sar_heating_sphere import RHO_C
from tests.validation.test_mat13_transient_perfusion_sphere import INTERIOR_BAND
from tests.validation.test_mat20_cem43_closed_form import _rel
from tests.validation.test_mat20_running_average import W_B2
from tests.validation.test_mat20_tuned_heating_dose import T_REF
from tests.validation.test_mat21_box_duty_dose_fem import LN4
from tests.validation.test_mat21_wall_duty_dose_closed_form import _q_check, _win_check

P_POLE = PTS[POLE]
P_MID = np.array([0.0, 0.0, 0.5 * (A + B)])
M_PTS = np.vstack([P_POLE, P_MID])
M_LABELS = ("pole", "mid")
DUTIES = (1.0, 0.5, 0.25)
#: Not the brief's pair.  ``N_DEFAULT // 2`` reads M_p at one ulp of 20 °C (the modes it
#: drops have decayed by the first grid row), so (Bd-i) is met trivially there; this
#: truncation is the pair that does perturb the history, (Bd-i′).
N_COARSE = 2


def _p(msg):
    print(msg, flush=True)


def _duty_hist_at(modes, pts, duty):
    """`_duty_hist` (bound to ``PTS``) at any points; row 0 = 0.0 exactly."""
    out = np.zeros((T_GRID.size, pts.shape[0]))
    for n in range(1, T_GRID.size):
        out[n] = (
            pec_sphere_implant_duty_cycle_from_modes(
                pts, T_GRID[n], modes, rho_c=RHO_C, t_inf=T_INF, period=PERIOD, duty=duty
            )
            - T_INF
        )
    return out


def _scales(duty):
    return ((1.0, SCALES[0][1]), (1.0 / duty, SCALES[1][1]))


@pytest.fixture(scope="module")
def hists():
    m = _modes()
    m2 = _modes(N_DEFAULT // 2)
    # The control builds its own modes: `_modes()` fixes beta (brief).
    iso = pec_sphere_implant_transient_modes(beta=0.0, n_terms=N_DEFAULT, inner="neumann", **FIX)
    mc = _modes(N_COARSE)
    out = {"full": {}, "half": {}, "coarse": {}, "iso": {}}
    for d in DUTIES:
        out["full"][d] = _duty_hist_at(m, M_PTS, d)
        out["half"][d] = _duty_hist_at(m2, M_PTS, d)
        out["coarse"][d] = _duty_hist_at(mc, M_PTS, d)
        out["iso"][d] = _duty_hist_at(iso, M_PTS, d)
    cont = np.zeros((T_GRID.size, M_PTS.shape[0]))
    for n in range(1, T_GRID.size):
        cont[n] = _theta(M_PTS, T_GRID[n], m)
    out["cont"] = cont
    return out


def test_precondition_grid():
    assert T_GRID.size == 161 and abs(T_GRID[-1] / T_END - 1.0) <= 1e-14
    assert T_STAR < W_B2 <= T_END, (T_STAR, W_B2, T_END)
    assert abs(PERIOD / DT - 20.0) < 1e-12
    for d, on in ((0.5, 10.0), (0.25, 5.0)):
        assert abs(d * PERIOD / DT - on) < 1e-12, d
    assert A < P_POLE[2] < P_MID[2] < B, (A, P_POLE, P_MID, B)
    _p(
        f"IMP7-10 grid: t* {T_STAR:.4f} s, T_END {T_END:.4f} s, dt {DT:.4f} s, period "
        f"{PERIOD:.4f} s; {dict(zip(M_LABELS, M_PTS.tolist()))}; N {N_DEFAULT} / {N_DEFAULT // 2}"
    )


def test_c100_m_full_duty_equals_continuous(hists):
    d1, cont = hists["full"][1.0], hists["cont"]
    worst = float(np.max(np.abs(d1[1:] - cont[1:]) / np.abs(cont[1:])))
    _p(
        f"IMP7-10 (C100-m) D = 1 vs continuous over {cont.shape}: worst rel {worst:.3e} (band "
        f"{B1_BAND:g}, ASSERTED)"
    )
    assert worst <= B1_BAND, worst


def test_q_m_win_m_bd_i_and_rc_m(hists):
    for d in DUTIES:
        for scale, lab in _scales(d):
            for j, name in enumerate(M_LABELS):
                th = hists["full"][d][:, j] * scale
                th2 = hists["half"][d][:, j] * scale
                T = T_INF + th
                dose, ref, rel, swap, miss = _q_check(T_GRID, T)
                mx, t_mx, lo, hi, wref, wrel = _win_check(T_GRID, th)
                M = float(np.max(np.abs(th2 - th)))
                m = abs(cem43(T_GRID, T_INF + th2) / dose - 1.0)
                cell = f"D={d} {lab} {name}"
                _p(
                    f"IMP7-10 (Q-m) {cell}: CEM43 {dose:.6e} s rel quad {rel:.3e} (bar "
                    f"{DOSE_BAR:g}, ASSERTED); max T {np.max(T):.6f} C; (Rc-m) miss {miss:.4e}"
                )
                _p(
                    f"IMP7-10 (Win-m) {cell}: max six-minute mean {mx:.6e} K at {t_mx:.2f} s; "
                    f"[{lo:.6e}, {hi:.6e}]; rel quad {wrel:.3e} (bar {AVG_BAR:g}, ASSERTED)"
                )
                _p(
                    f"IMP7-10 (Bd-i) {cell}: m_p {m:.4e} <= ln4*M_p {LN4 * M:.4e} (M_p {M:.4e} K; "
                    f"exact {math.expm1(LN4 * M):.4e}, printed only); ratio {m / (LN4 * M):.4f}"
                )
                assert float(np.max(T)) < T_REF, (cell, float(np.max(T)))
                assert rel <= DOSE_BAR, (cell, rel)
                assert miss > DOSE_BAR, (cell, miss)
                assert T_GRID[0] + W_B2 <= t_mx <= T_GRID[-1], (cell, t_mx)
                assert lo <= mx <= hi, (cell, lo, mx, hi)
                assert wrel <= AVG_BAR, (cell, wrel)
                assert m <= LN4 * M, (cell, m, LN4 * M)


def test_bd_i_prime_on_the_coarse_pair(hists):
    """(Bd-i′): the same inequality on the pair that moves the history; see `N_COARSE`."""
    for d in DUTIES:
        for scale, lab in _scales(d):
            for j, name in enumerate(M_LABELS):
                th = hists["full"][d][:, j] * scale
                thc = hists["coarse"][d][:, j] * scale
                M = float(np.max(np.abs(thc - th)))
                m = abs(cem43(T_GRID, T_INF + thc) / cem43(T_GRID, T_INF + th) - 1.0)
                _p(
                    f"IMP7-10 (Bd-i') D={d} {lab} {name} N={N_COARSE}: m_p {m:.4e} <= ln4*M_p "
                    f"{LN4 * M:.4e} (M_p {M:.4e} K); ratio {m / (LN4 * M):.4f}"
                )
                assert m <= LN4 * M, (d, lab, name, m, LN4 * M)


def test_beta0_m_and_mis_m_controls(hists):
    for d in DUTIES:
        d_b = cem43(T_GRID, T_INF + hists["full"][d][:, POLE])
        d_0 = cem43(T_GRID, T_INF + hists["iso"][d][:, POLE])
        miss = _rel(d_0, d_b)
        _p(
            f"IMP7-10 (b0-m) D={d} pole: ratio {d_b / d_0:.6f}, miss {miss:.4e} (asserted > "
            f"{DOSE_BAR:g})"
        )
        assert miss > DOSE_BAR, (d, miss)
    for j, name in enumerate(M_LABELS):
        d25 = cem43(T_GRID, T_INF + hists["full"][0.25][:, j])
        d50 = cem43(T_GRID, T_INF + hists["full"][0.5][:, j])
        miss = _rel(d25, d50)
        _p(
            f"IMP7-10 (Mis-m) {name}: D = 0.25 against D = 0.5 ratio {d25 / d50:.6f}, miss "
            f"{miss:.4e} (asserted > {DOSE_BAR:g})"
        )
        assert miss > DOSE_BAR, (name, miss)


def test_printed_matrix(hists):
    """The records (no band); asserts only that they were formed."""
    _p(
        "IMP7-10 MATRIX D | point | D(peak) s | D(equal mean) s | six-minute peak K | "
        "cap INTERIOR_BAND*max theta_s K | M_p K | headroom cap/M_p | m_p/(ln4*M_p)"
    )
    for d in DUTIES:
        for j, name in enumerate(M_LABELS):
            th, th2 = hists["full"][d][:, j], hists["half"][d][:, j]
            d_pk = cem43(T_GRID, T_INF + th)
            d_eq = cem43(T_GRID, T_INF + th / d)
            mx = _win_check(T_GRID, th)[0]
            cap = INTERIOR_BAND * float(np.max(th))
            M = float(np.max(np.abs(th2 - th)))
            ratio = abs(cem43(T_GRID, T_INF + th2) / d_pk - 1.0) / (LN4 * M)
            _p(
                f"IMP7-10 MATRIX {d} | {name} | {d_pk:.6e} | {d_eq:.6e} | {mx:.6f} | {cap:.4e} | "
                f"{M:.4e} | {cap / M:.3e} | {ratio:.4f}"
            )
            assert np.isfinite(d_eq / d_pk)
