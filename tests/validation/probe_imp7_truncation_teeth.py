"""Measure-only probe for §9 item 741's red (Tr20): at which n does truncation move the
pole's dose against 200 modes?  Asserts nothing, closes nothing; the review rules."""

import numpy as np

from fem_em_solver.thermal.dose import cem43
from tests.validation.test_imp7_dose_comparand_floor import T_GRID, T_INF, _duty_hist, _modes, _rel
from tests.validation.test_imp7_duty_dose_closed_form import LABELS

ref = _duty_hist(_modes(200))
for n in (1, 2, 3, 4, 6, 8, 12, 20):
    th = _duty_hist(_modes(n))
    rels = [
        _rel(cem43(T_GRID, T_INF + th[:, j]), cem43(T_GRID, T_INF + ref[:, j])) for j in range(3)
    ]
    diff = np.abs(th - ref)
    row = np.argmax(diff, axis=0)
    print(
        f"PROBE n={n}: rel dose vs 200 "
        + ", ".join(f"{name} {r:.3e}" for name, r in zip(LABELS, rels))
        + "; max |dtheta| "
        + ", ".join(f"{diff[row[j], j]:.3e} K (row {row[j]})" for j in range(3)),
        flush=True,
    )
