"""`MAT-21` step 23 — the dose and six-minute window on P49b's walled-box CN
history (§9 item 708, §10 P50b; brief `docs/planning/chunks/MAT-21.md` "§9 item 708").

Item 688's fixture by name (`test_mat21_wall_duty_fem`: the 1D-reduced walled box,
rung m = 4, CN with two BE start-up half-steps, `N2` steps to t₂ = 8 t*), with
snapshots at every step, and P50a's grid and checks imported from
`test_mat21_wall_duty_dose_closed_form`.  Each row of ``T_GRID`` reads the three
``Z_PTS`` in one fixed batch at (x, y) = (0, 0), `_probe_vals`'s owner route.

Asserted (bands imported, none set): (T) gel-interior rel L2 ≤ ``INTERIOR_BAND``
at 4 t* against `_exact_duty_t`; (E) max balance ≤ ``E_BAND``; (S0) ≤ ``S0_BAND``;
every point valid; max T < ``T_REF``; (Q-wf) / (Win-wf) on the FEM history at both
``SCALES``, at ``DOSE_BAR`` / ``AVG_BAR``, with the ``t_max`` range and the bracket.
Negative controls (asserted): the massless-wall series misses the FEM's interface
dose by > ``DOSE_BAR``; the R swap by > ``DOSE_BAR``; the continuous-source series
misses at 4 t* by > ``INTERIOR_BAND``.

Printed (records, no band): at each point the FEM dose beside P50a's series dose,
peak and equal-mean; each six-minute peak and t_max; max_n|Δθ|.

Not claimed: a dose band, a rate, the six-wall box, a coil source.  k_w and ρ_w c_w
remain nominal.
"""

from __future__ import annotations

import time

import numpy as np
import pytest

from fem_em_solver.thermal import solve_transient_conduction
from fem_em_solver.thermal.dose import cem43, running_average
from tests.validation.test_imp7_duty_cycle_fem import DUTY, S0_BAND, _chi
from tests.validation.test_mat10_sar_to_heat_sphere import _say
from tests.validation.test_mat12_transient_sar_heating_sphere import E_BAND
from tests.validation.test_mat13_transient_perfusion_sphere import INTERIOR_BAND
from tests.validation.test_mat15_coil_driven_heating_phantom import H_GEL, K_GEL, N_MAIN
from tests.validation.test_mat15_finite_cylinder_closed_form import T_INF, H, K, Q
from tests.validation.test_mat15_finite_cylinder_fem import _rel as _rel_l2
from tests.validation.test_mat20_cem43_closed_form import _rel
from tests.validation.test_mat20_running_average import W_B2
from tests.validation.test_mat20_tuned_heating_dose import T_REF
from tests.validation.test_mat21_box_closed_form import A_BOX, B_BOX, C_BOX
from tests.validation.test_mat21_box_probe_points_fem import _probe_vals
from tests.validation.test_mat21_wall_closed_form import D_WALL, RHO_C_W
from tests.validation.test_mat21_wall_duty_closed_form import PERIOD
from tests.validation.test_mat21_wall_duty_dose_closed_form import (
    AVG_BAR,
    DOSE_BAR,
    INTERFACE,
    LABELS,
    SCALES,
    T_GRID,
    Z_PTS,
    _duty_hist,
    _q_check,
    _win_check,
)
from tests.validation.test_mat21_wall_duty_fem import DT, _census, _exact_duty_t
from tests.validation.test_mat21_wall_fem import (
    ROBIN,
    _materials,
    _robin_facets,
    _walled_box,
)
from tests.validation.test_mat21_wall_transient_closed_form import T_STAR
from tests.validation.test_mat21_wall_transient_fem import (
    MASSLESS,
    N2,
    RUNG,
    _exact_t,
    _lam,
    _rho_c_field,
)

STEP_4T = 4 * N_MAIN
T_4 = 4.0 * T_STAR
BATCH = np.array([(0.0, 0.0, z) for z in Z_PTS], dtype=float)


def _six_minute(th):
    """``(mx, t_mx)`` of the six-minute window on one history."""
    _a, mx, t_mx = running_average(T_GRID, th, window=float(W_B2))
    return float(mx), float(t_mx)


def test_precondition_grid_and_batch():
    """Imports and grid (W0): no solve."""
    assert T_GRID.size == N2 + 1 == 321, T_GRID.size
    assert abs(T_GRID[1] - DT) <= 1e-12 * DT
    assert abs(T_GRID[STEP_4T] / T_4 - 1.0) <= 1e-14
    assert abs(PERIOD / DT - 20.0) < 1e-12 and abs(DUTY * PERIOD / DT - 10.0) < 1e-12, (PERIOD, DT)
    assert K == K_GEL and H == H_GEL, (K, K_GEL, H, H_GEL)
    assert BATCH.shape == (3, 3) and LABELS[INTERFACE] == "interface"
    assert BATCH[INTERFACE, 2] == C_BOX
    _say(
        f"\nMAT21-23 (W0) T_GRID {T_GRID.size} rows, dt {DT:.4f} s, 4 t* {T_4:.3f} s at row {STEP_4T}; "
        f"batch z {tuple(float(z) for z in Z_PTS)} m"
    )


@pytest.fixture(scope="module")
def wall_dose():
    t0 = time.perf_counter()
    msh, ns, nw = _walled_box(RUNG, (False, False, True))
    n = msh.topology.index_map(3).size_global
    chi, kf, inter, _ct = _materials(msh)
    rcf = _rho_c_field(msh)
    ft = _robin_facets(msh, (A_BOX, B_BOX, C_BOX + D_WALL), (2,))
    tl = time.perf_counter()
    res = solve_transient_conduction(
        msh,
        ft,
        k=kf,
        rho_c=rcf,
        q=Q * chi,
        h=H_GEL,
        t_inf=T_INF,
        t_init=T_INF,
        robin_tag=ROBIN,
        theta=0.5,
        startup_half_steps=2,
        dt=DT,
        n_steps=N2,
        snapshot_steps=tuple(range(N2 + 1)),
        source_schedule=_chi(PERIOD, DUTY),
    )
    t_cn = time.perf_counter() - tl
    s0, tot, n_on, n_off = _census(res)
    _say(
        f"\nMAT21-23 [W1] walled box m = {RUNG} cells {tuple(ns)} N_w = {nw}: n = {n}; CN {t_cn:.1f} s; "
        f"(E) max balance {res.max_balance_residual:.3e} (band {E_BAND:g}); (S0) sum chi dt {tot:.6f} s: "
        f"rel {s0:.3e} (band {S0_BAND:g}); {n_on} on / {n_off} off (half-)steps"
    )

    tp = time.perf_counter()
    hist = np.zeros((T_GRID.size, BATCH.shape[0]))
    valid = np.ones(BATCH.shape[0], dtype=bool)
    for k in range(T_GRID.size):
        hist[k], ok = _probe_vals(res.snapshots[k], BATCH)
        valid &= ok
    t_probe = time.perf_counter() - tp

    th = res.snapshots[STEP_4T] - T_INF
    lam_w = _lam(RHO_C_W)
    e_int = _rel_l2(msh, th, _exact_duty_t(msh, T_4, RHO_C_W, lam_w), inter)
    miss_c = _rel_l2(msh, th, _exact_t(msh, T_4, RHO_C_W, lam_w), inter)
    _say(
        f"MAT21-23 (T) 4 t* = {T_4:.3f} s: interior rel L2 vs duty series {e_int:.6e} (band "
        f"{INTERIOR_BAND:.0e}); continuous-source series miss {miss_c:.6e}; point reads {t_probe:.1f} s; "
        f"all valid {bool(np.all(valid))}"
    )
    series = _duty_hist(RHO_C_W)
    massless = _duty_hist(MASSLESS * RHO_C_W)
    _say(f"MAT21-23 [W1] elapsed {time.perf_counter() - t0:.1f} s")
    return {
        "hist": hist,
        "valid": valid,
        "bal": res.max_balance_residual,
        "s0": s0,
        "e_int": e_int,
        "miss_c": miss_c,
        "series": series,
        "massless": massless,
    }


def test_all_gates(wall_dose):
    """Collect every gate before asserting (item 668), so one red hides no other."""
    w = wall_dose
    hist = w["hist"]
    gates = {
        "T": w["e_int"] <= INTERIOR_BAND,
        "E": w["bal"] <= E_BAND,
        "S0": w["s0"] <= S0_BAND,
        "valid": bool(np.all(w["valid"])),
        "continuous_control": w["miss_c"] > INTERIOR_BAND,
    }
    for scale, lab in SCALES:
        for j, name in enumerate(LABELS):
            th = hist[:, j] * scale
            T = T_INF + th
            dose, ref, rel, _swap, miss = _q_check(T_GRID, T)
            mx, t_mx, lo, hi, wref, wrel = _win_check(T_GRID, th)
            _say(
                f"MAT21-23 (Q-wf) {lab} {name}: FEM CEM43 {dose:.6e} s vs quad {ref:.6e} s rel "
                f"{rel:.3e} (bar {DOSE_BAR:g}); max T {np.max(T):.6f} C (< {T_REF}); R swap miss "
                f"{miss:.4e} (> {DOSE_BAR:g})"
            )
            _say(
                f"MAT21-23 (Win-wf) {lab} {name}: max six-minute mean {mx:.6e} K at t_max {t_mx:.4f} s; "
                f"bracket [{lo:.6e}, {hi:.6e}]; quad {wref:.6e} rel {wrel:.3e} (bar {AVG_BAR:g})"
            )
            key = f"{lab} {name}"
            gates[f"Q {key}"] = rel <= DOSE_BAR
            gates[f"maxT {key}"] = float(np.max(T)) < T_REF
            gates[f"Win {key}"] = wrel <= AVG_BAR
            gates[f"tmax {key}"] = T_GRID[0] + W_B2 <= t_mx <= T_GRID[-1]
            gates[f"bracket {key}"] = lo <= mx <= hi
            gates[f"Rswap {key}"] = miss > DOSE_BAR
    d_fem = cem43(T_GRID, T_INF + hist[:, INTERFACE])
    d_m = cem43(T_GRID, T_INF + w["massless"][:, INTERFACE])
    miss_m = _rel(d_m, d_fem)
    gates["massless_control"] = miss_m > DOSE_BAR
    _say(
        f"MAT21-23 (control) interface dose: massless-wall series misses the FEM by {miss_m:.4e} "
        f"(separation {miss_m / DOSE_BAR:.3e} x DOSE_BAR; item 698's series-side miss 1.1165e-01); "
        f"continuous-source series at 4 t* {w['miss_c']:.6e} vs INTERIOR_BAND {INTERIOR_BAND:.0e} "
        f"(separation {w['miss_c'] / INTERIOR_BAND:.2f}x)"
    )
    for j, name in enumerate(LABELS):
        fem_h, ser_h = hist[:, j], w["series"][:, j]
        rows = [f"max_n |dtheta| {np.max(np.abs(fem_h - ser_h)):.6e} K"]
        for scale, lab in SCALES:
            df, ds = cem43(T_GRID, T_INF + fem_h * scale), cem43(T_GRID, T_INF + ser_h * scale)
            mf, tf = _six_minute(fem_h * scale)
            ms, ts = _six_minute(ser_h * scale)
            rows.append(
                f"{lab}: dose FEM {df:.6e} s, series {ds:.6e} s (rel {_rel(df, ds):.3e}); "
                f"six-minute FEM {mf:.6f} K at {tf:.2f} s, series {ms:.6f} K at {ts:.2f} s"
            )
        _say(f"MAT21-23 RECORD {name}: " + "; ".join(rows))
    _say(f"MAT21-23 gates red: {[k for k, v in gates.items() if not v] or 'none'} of {len(gates)}")
    assert all(gates.values()), [k for k, v in gates.items() if not v]
