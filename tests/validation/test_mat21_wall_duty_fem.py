"""`MAT-21` step 21b — the duty-cycled walled box, Crank–Nicolson against the
duty-cycled two-layer slab closed form (§9 item 688, §10 P49b; brief
`docs/planning/chunks/MAT-21.md` "§9 item 688").

Item 682's fixture (`test_mat21_wall_transient_fem`: the 1D-reduced walled box,
rung m = 4, CN θ = 0.5 with two backward-Euler start-up half-steps, Δt =
t*/N_MAIN) with one change: ``source_schedule=_chi(PERIOD, DUTY)`` (`_chi`,
`DUTY` from item 675's `test_imp7_duty_cycle_fem`).  `PERIOD` is item 689's
(`test_mat21_wall_duty_closed_form.PERIOD` = 0.5 × MAT-21's t* = 1410.75 s —
item 675's `PERIOD` is 0.5 × *IMP-7's* t*, so its ratio is applied here, as
item 689 recorded), so a period is 20 CN steps and an on-phase 10.

Comparand: item 689's `two_layer_slab_duty_cycle_temperature` (400 modes) in z,
interpolated into P2 in chunks; errors on θ = T − T∞.

Asserted (imported gates, unmoved): gel-interior rel L2 ≤ ``INTERIOR_BAND`` at
t* and t₂ = 8 t*; max per-step balance ≤ ``E_BAND``; tag volumes ≤
``_vol_bar(n)``; (S0) the schedule census Σχ·Δt equals the analytic on-time
D·t₂ to ≤ ``S0_BAND``.  Negative control (asserted): the continuous-source
series misses the FEM at t* by more than ``CTL_FACTOR`` × the FEM's interior
error there.  The massless-wall duty series' miss is printed, not asserted.
"""

from __future__ import annotations

import time

import numpy as np
import pytest
from dolfinx import fem

from fem_em_solver.thermal import solve_transient_conduction
from fem_em_solver.utils.analytical import two_layer_slab_duty_cycle_temperature

from tests.validation.test_imp7_duty_cycle_fem import DUTY, S0_BAND, _chi
from tests.validation.test_mat10_sar_to_heat_sphere import _say
from tests.validation.test_mat12_transient_sar_heating_sphere import E_BAND
from tests.validation.test_mat13_transient_perfusion_sphere import INTERIOR_BAND
from tests.validation.test_mat15_coil_driven_heating_phantom import H_GEL, K_GEL, N_MAIN
from tests.validation.test_mat15_finite_cylinder_closed_form import H, K, Q, RHO_C, T_INF
from tests.validation.test_mat15_finite_cylinder_fem import _rel
from tests.validation.test_mat21_box_closed_form import A_BOX, B_BOX, C_BOX, T_C
from tests.validation.test_mat21_wall_closed_form import D_WALL, K_W, RHO_C_W
from tests.validation.test_mat21_wall_duty_closed_form import PERIOD
from tests.validation.test_mat21_wall_fem import (
    _centre,
    _materials,
    _robin_facets,
    _vol_bar,
    _volumes,
    _walled_box,
)
from tests.validation.test_mat21_wall_transient_closed_form import N_SERIES, T_STAR
from tests.validation.test_mat21_wall_transient_fem import (
    CTL_FACTOR,
    INTERP_CHUNK,
    MASSLESS,
    N2,
    RUNG,
    SLAB,
    T2,
    _exact_t,
    _lam,
    _rho_c_field,
)
from tests.validation.test_mat21_wall_fem import ROBIN

DT = T_STAR / N_MAIN
ON_TIME = DUTY * T2  # 16 whole periods to t₂


def _duty_series(z, t, rho_c_w, lam):
    return two_layer_slab_duty_cycle_temperature(
        z, t, q=Q, h=H, t_inf=0.0, rho_c_w=rho_c_w, period=PERIOD, duty=DUTY,
        n_terms=N_SERIES, eigenvalues=lam, **SLAB)


def _exact_duty_t(msh, t, rho_c_w, lam):
    V2 = fem.functionspace(msh, ("Lagrange", 2))
    g = fem.Function(V2)

    def ev(x):
        n = x.shape[1]
        out = np.empty(n)
        for s in range(0, n, INTERP_CHUNK):
            out[s:s + INTERP_CHUNK] = _duty_series(x[2, s:s + INTERP_CHUNK], t, rho_c_w, lam)
        return out

    g.interpolate(ev)
    return g


def _census(res, startup=2):
    m = np.asarray(res.source_multipliers, dtype=float)
    dts = np.array([0.5 * DT] * startup + [DT] * (len(m) - startup))
    tot = float(np.dot(m, dts))
    return abs(tot - ON_TIME) / ON_TIME, tot, int(np.sum(m == 1.0)), int(np.sum(m == 0.0))


@pytest.fixture(scope="module")
def walled():
    t0 = time.perf_counter()
    assert K == K_GEL and H == H_GEL, (K, K_GEL, H, H_GEL)
    # Schedule shape the brief requires: a period is 20 steps, an on-phase 10.
    assert abs(PERIOD / DT - 20.0) < 1e-12 and abs(DUTY * PERIOD / DT - 10.0) < 1e-12, (PERIOD, DT)
    msh, ns, nw = _walled_box(RUNG, (False, False, True))
    n = msh.topology.index_map(3).size_global
    chi, kf, inter, ct = _materials(msh)
    rcf = _rho_c_field(msh)
    ft = _robin_facets(msh, (A_BOX, B_BOX, C_BOX + D_WALL), (2,))
    vg, vw = _volumes(msh, ct)
    v_gel = 8.0 * A_BOX * B_BOX * C_BOX
    v_wall = 4.0 * A_BOX * B_BOX * 2.0 * D_WALL
    out = {"n": n, "dvg": abs(vg / v_gel - 1.0), "dvw": abs(vw / v_wall - 1.0)}
    _say(f"\nMAT21-21b [W1] 1D-reduced walled box m = {RUNG} cells {tuple(ns)} N_w = {nw}: "
         f"n = {n}; V_gel rel dev {out['dvg']:.3e}, V_wall rel dev {out['dvw']:.3e} "
         f"(bar {_vol_bar(n):.3e}); T_C {T_C:.3f} s, t* {T_STAR:.3f} s, t2 {T2:.3f} s; "
         f"P {PERIOD:.4f} s = {PERIOD / DT:.1f} steps, D {DUTY}; dt {DT:.4f} s, {N2} CN steps")
    tl = time.perf_counter()
    res = solve_transient_conduction(
        msh, ft, k=kf, rho_c=rcf, q=Q * chi, h=H_GEL, t_inf=T_INF, t_init=T_INF,
        robin_tag=ROBIN, theta=0.5, startup_half_steps=2, dt=DT, n_steps=N2,
        snapshot_steps=(N_MAIN, N2), source_schedule=_chi(PERIOD, DUTY))
    _say(f"MAT21-21b [W1] CN solve {time.perf_counter() - tl:.1f} s; max balance "
         f"{res.max_balance_residual:.3e} over {len(res.balance_residuals)} (half-)steps")
    out["bal"] = res.max_balance_residual
    s0, tot, n_on, n_off = _census(res)
    out["s0"] = s0
    _say(f"MAT21-21b (S0) census sum chi dt = {tot:.6f} s vs D t2 = {ON_TIME:.6f} s: rel "
         f"{s0:.3e} (band {S0_BAND:g}); {n_on} on / {n_off} off (half-)steps")
    lam_w = _lam(RHO_C_W)
    lam_0 = _lam(MASSLESS * RHO_C_W)
    out["snap"] = {}
    for step, t in ((N_MAIN, T_STAR), (N2, T2)):
        th = res.snapshots[step] - T_INF
        ex = _exact_duty_t(msh, t, RHO_C_W, lam_w)
        e_int = _rel(msh, th, ex, inter)
        e_whole = _rel(msh, th, ex)
        exc = _exact_t(msh, t, RHO_C_W, lam_w)  # continuous source (negative control)
        miss_c = _rel(msh, th, exc, inter)
        ex0 = _exact_duty_t(msh, t, MASSLESS * RHO_C_W, lam_0)
        miss0 = _rel(msh, th, ex0, inter)
        c_h = _centre(res.snapshots[step]) - T_INF
        c_ex = float(np.atleast_1d(_duty_series(0.0, t, RHO_C_W, lam_w))[0])
        out["snap"][step] = {"e_int": e_int, "e_whole": e_whole, "miss_c": miss_c,
                             "miss0": miss0, "c_h": c_h, "c_ex": c_ex}
        _say(f"MAT21-21b [W1] step {step} t = {t:.3f} s: interior rel L2 vs duty series "
             f"{e_int:.6e} (bound {INTERIOR_BAND:.0e}); whole {e_whole:.6e} (printed); "
             f"continuous-source series miss {miss_c:.6e}; massless-wall duty series miss "
             f"{miss0:.6e} (printed)")
        _say(f"MAT21-21b [W1] step {step} centre theta: FEM {c_h:.9f} K, duty series "
             f"{c_ex:.9f} K (rel {c_h / c_ex - 1.0:+.3e})")
    _say(f"MAT21-21b [W1] elapsed {time.perf_counter() - t0:.1f} s")
    return out


def test_all_gates(walled):
    """Collect every gate before asserting (item 668), so one red hides no other."""
    w = walled
    s1 = w["snap"][N_MAIN]
    gates = {
        "vol_gel": w["dvg"] <= _vol_bar(w["n"]),
        "vol_wall": w["dvw"] <= _vol_bar(w["n"]),
        "balance": w["bal"] <= E_BAND,
        "S0": w["s0"] <= S0_BAND,
        "interior_tstar": s1["e_int"] <= INTERIOR_BAND,
        "interior_t2": w["snap"][N2]["e_int"] <= INTERIOR_BAND,
        "control_continuous_tstar": s1["miss_c"] > CTL_FACTOR * s1["e_int"],
    }
    _say(f"MAT21-21b (control) at t*: continuous-source series misses by {s1['miss_c']:.6e} vs "
         f"{CTL_FACTOR:.0f} x FEM interior error {CTL_FACTOR * s1['e_int']:.6e} "
         f"(separation {s1['miss_c'] / s1['e_int']:.2f}x)")
    _say(f"MAT21-21b gates: {gates}")
    assert all(gates.values()), (gates, w)
