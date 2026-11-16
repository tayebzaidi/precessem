"""`IMP-7` step 14 — the comparand's own floor: the hollow-sphere duty series at 200, 400
and 800 modes (§9 item 741, §10 chain I, I53; brief `docs/planning/chunks/IMP-7.md` "§9
item 741").

Numpy (the imported fixture modules bring dolfinx in).  Sets no band and moves none; every
band is imported.  Histories: `_duty_hist(_modes(n))` on ``T_GRID`` at ``PTS``, D = 0.5.

Item 725 left about four fifths of the pole's FEM-vs-series dose gap in something neither
h nor Δt moves at one refinement; item 745 then read it as not Δt.  Nothing had tested the
comparand itself, which ``N_DEFAULT`` = 200 modes truncates.

Asserted: (Root) the 800-mode root ladder passes `test_imp7_4_eigenvalue_count`'s own
checks; (Q-hd) / (Win-hd) `cem43` against `quad` at ``DOSE_BAR`` and the six-minute window
at ``AVG_BAR``, on each n, point and scale, max T < 43 °C; (Self) the 800-mode dose at the
pole against the 200-mode dose ≤ ``B1_BAND``.  Negative controls: (Tr8) n = 8 misses the
200-mode dose at the pole by > ``DOSE_BAR``; (Rc-hd) the R swap misses by > ``DOSE_BAR``.
n = 8 is from the probe's log, the brief's n = 20 read 1.493e-16.

Printed (no band): rel dose change 200 → 400 → 800 at each point on the peak scale, beside
item 725's gaps, and the row of ``T_GRID`` with the largest |Δθ| between 200 and 800.

Not claimed: a dose accuracy, a rate, or that the residual is the comparand unless the
printed change reaches the gap's size.
"""

from __future__ import annotations

import numpy as np
import pytest

from fem_em_solver.thermal.dose import cem43

from tests.validation.test_imp4_pec_cylinder_closed_form import B1_BAND
from tests.validation.test_imp7_duty_dose_closed_form import (
    AVG_BAR,
    DOSE_BAR,
    LABELS,
    POLE,
    SCALES,
    T_GRID,
    _duty_hist,
)
from tests.validation.test_imp7_implant_heating_closed_form import T_INF
from tests.validation.test_imp7_transient_closed_form import N_DEFAULT, SPAN, _modes
from tests.validation.test_mat20_cem43_closed_form import _rel
from tests.validation.test_mat20_running_average import W_B2
from tests.validation.test_mat20_tuned_heating_dose import T_REF
from tests.validation.test_mat21_wall_duty_dose_closed_form import _q_check, _win_check

N_LADDER = (N_DEFAULT, 400, 800)
N_CONTROL = 8
#: Item 725's peak-scale "both" row, `20261101T182633Z_IMP-7-step12-W1.log`, at the log's
#: ``.3e`` (the brief's list).
ITEM_725_GAP = (1.011e-04, 1.251e-05, 1.724e-07)


def _p(msg):
    print(msg, flush=True)


@pytest.fixture(scope="module")
def modes():
    return {n: _modes(n) for n in (N_CONTROL, *N_LADDER)}


@pytest.fixture(scope="module")
def hists(modes):
    return {n: _duty_hist(m) for n, m in modes.items()}


def test_root_ladder_at_800(modes):
    """First in W1: the checks of `test_imp7_4_eigenvalue_count`, on the modes in use."""
    n_max = N_LADDER[-1]
    for l in (0, 2):
        lam = modes[n_max][l]["lam"]
        x = lam * SPAN / np.pi
        idx = np.arange(lam.size)
        dev = float(np.max(np.abs(x[10:] - idx[10:])))
        gap = float(np.min(np.diff(x)))
        _p(f"IMP7-14 (Root) l={l}: {lam.size} roots, x0 {x[0]:.6f}, x_last {x[-1]:.6f}; max "
           f"|x_n - n| (n>=10) {dev:.3e} (asserted < 0.5); min gap {gap:.6f} spacings "
           f"(asserted > 0.125)")
        assert lam.size == n_max
        assert lam[0] > 0.0
        assert dev < 0.5, (l, dev)
        assert gap > 1.0 / 8.0, (l, gap)


def test_q_hd_win_hd_and_rc_hd(hists):
    for n in N_LADDER:
        for scale, lab in SCALES:
            for j, name in enumerate(LABELS):
                th = hists[n][:, j] * scale
                T = T_INF + th
                dose, _ref, rel, _swap, miss = _q_check(T_GRID, T)
                mx, t_mx, lo, hi, _wref, wrel = _win_check(T_GRID, th)
                cell = f"n={n} {lab} {name}"
                _p(f"IMP7-14 (Q-hd) {cell}: CEM43 {dose:.6e} s rel quad {rel:.3e} (bar "
                   f"{DOSE_BAR:g}, ASSERTED); max T {np.max(T):.6f} C; (Rc-hd) miss {miss:.4e} "
                   f"(asserted > {DOSE_BAR:g})")
                _p(f"IMP7-14 (Win-hd) {cell}: max six-minute mean {mx:.6e} K at {t_mx:.2f} s; "
                   f"rel quad {wrel:.3e} (bar {AVG_BAR:g}, ASSERTED)")
                assert float(np.max(T)) < T_REF, (cell, float(np.max(T)))
                assert rel <= DOSE_BAR, (cell, rel)
                assert miss > DOSE_BAR, (cell, miss)
                assert T_GRID[0] + W_B2 <= t_mx <= T_GRID[-1], (cell, t_mx)
                assert lo <= mx <= hi, (cell, lo, mx, hi)
                assert wrel <= AVG_BAR, (cell, wrel)


def test_self_800_against_200_at_the_pole(hists):
    d_200 = cem43(T_GRID, T_INF + hists[N_DEFAULT][:, POLE])
    d_800 = cem43(T_GRID, T_INF + hists[N_LADDER[-1]][:, POLE])
    rel = _rel(d_800, d_200)
    _p(f"IMP7-14 (Self) pole: rel dose 800 vs 200 modes {rel:.3e} (band {B1_BAND:g}, ASSERTED); "
       f"item 725's gap {ITEM_725_GAP[POLE]:.3e}")
    assert rel <= B1_BAND, rel


def test_tr8_truncation_control(hists):
    d_8 = cem43(T_GRID, T_INF + hists[N_CONTROL][:, POLE])
    d_200 = cem43(T_GRID, T_INF + hists[N_DEFAULT][:, POLE])
    miss = _rel(d_8, d_200)
    _p(f"IMP7-14 (Tr8) pole: n = {N_CONTROL} against n = {N_DEFAULT} miss {miss:.4e} (asserted "
       f"> {DOSE_BAR:g})")
    assert miss > DOSE_BAR, miss


def test_printed_changes(hists):
    """The records (no band); asserts only that they were formed."""
    th = {n: hists[n] for n in N_LADDER}
    for j, name in enumerate(LABELS):
        d = [cem43(T_GRID, T_INF + th[n][:, j]) for n in N_LADDER]
        a, b = _rel(d[1], d[0]), _rel(d[2], d[1])
        _p(f"IMP7-14 CHANGE {name}: peak-scale rel dose 200->400 {a:.3e}, 400->800 {b:.3e}, "
           f"200->800 {_rel(d[2], d[0]):.3e}; item 725's gap {ITEM_725_GAP[j]:.3e}")
        assert np.isfinite(a) and np.isfinite(b)
    diff = np.abs(th[N_LADDER[-1]] - th[N_DEFAULT])
    row, col = np.unravel_index(int(np.argmax(diff)), diff.shape)
    _p(f"IMP7-14 ROW largest |dtheta| 200 vs 800: {diff[row, col]:.3e} K at row {row}, t "
       f"{T_GRID[row]:.2f} s, {LABELS[col]}; per point max "
       + ", ".join(f"{name} {np.max(diff[:, j]):.3e} K (row {int(np.argmax(diff[:, j]))})"
                   for j, name in enumerate(LABELS)))
    assert np.isfinite(diff[row, col])
