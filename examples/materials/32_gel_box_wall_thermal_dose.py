"""Example (`EX-166`): the thermal dose in the walled gel box -- CEM43 and the six-minute
mean on FEM point histories, against the two-layer slab series.

`mat:31` steps the acrylic-walled gel box under a duty-cycled source and compares the
*field* with the series. This example turns the same Crank-Nicolson history into the
quantities a safety reader asks for: the CEM43 dose and the maximum six-minute running
mean of the temperature rise, read at three points on the box axis (centre, near the wall,
and the gel-wall interface), where the wall's heat capacity matters most. The box is the
1D-reduced walled box of `MAT-21` step 20 (walls on the z faces only, so the two-layer
series is exact for the 3D FEM), m = 4 (119 808 cells), dt = t*/40, 320 CN steps to
t2 = 8 t*, period 0.5 t*, duty 0.5. Snapshots are held at every step and read at the three
points with the owner-route point evaluation.

**It asserts, and it does not re-implement** (the `ANS-1` rule). Every solver, series,
band and helper is imported by name from the gate chain of
``tests/validation/test_mat21_wall_duty_dose_fem.py`` (`MAT-21` step 23, item 708). That
module's ``wall_dose`` fixture is a pytest fixture and is not called: its body is written
out here in ``main()`` (disclosed, as `EX-165` did), with ``np.real`` on each integral
because the `mat:` group sources the complex build. No gate edit, no ``src/`` edit.

Asserted (bands imported, all collected before asserting): (E) <= ``E_BAND``; (S0) <=
``S0_BAND``; every point valid; max T < ``T_REF``; (T) interior rel L2 at 4 t* <=
``INTERIOR_BAND``; (Q-wf) the FEM history's ``cem43`` against the quadrature of the same
piecewise history at ``DOSE_BAR`` and the six-minute window check at ``AVG_BAR`` at each
point and both scales (an identity of ``cem43`` on a history, not a dose accuracy);
(c-mass) the massless-wall series misses the FEM's interface dose by more than
``DOSE_BAR``; (c-cont) the continuous-source series at 4 t* misses by more than
``INTERIOR_BAND``.
Printed only: the FEM-vs-series dose table. No dose band exists on this fixture.

Run it::

    ./run_examples.sh -e mat:32 -n 2 -t 480

Output: ``examples/materials/paraview_output/
materials_32_gel_box_wall_thermal_dose_combined.xdmf`` (theta at 4 t* and t2, the duty
series at 4 t*, and the error at 4 t*).
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
from dolfinx import default_scalar_type, fem
from mpi4py import MPI

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from fem_em_solver.io.paraview_utils import (  # noqa: E402
    adopt_host_ownership,
    write_xdmf_with_tags,
)
from fem_em_solver.post.setup_figure import write_setup_figure  # noqa: E402
from fem_em_solver.thermal import solve_transient_conduction  # noqa: E402
from fem_em_solver.thermal.dose import cem43, running_average  # noqa: E402

from tests.validation.test_imp7_duty_cycle_fem import DUTY, S0_BAND, _chi  # noqa: E402
from tests.validation.test_mat12_transient_sar_heating_sphere import E_BAND  # noqa: E402
from tests.validation.test_mat13_transient_perfusion_sphere import INTERIOR_BAND  # noqa: E402
from tests.validation.test_mat15_coil_driven_heating_phantom import (  # noqa: E402
    H_GEL,
    K_GEL,
    N_MAIN,
)
from tests.validation.test_mat15_finite_cylinder_closed_form import H, K, Q, T_INF  # noqa: E402
from tests.validation.test_mat15_finite_cylinder_fem import _max_imag  # noqa: E402
from tests.validation.test_mat15_finite_cylinder_fem import _rel as _rel_l2  # noqa: E402
from tests.validation.test_mat20_cem43_closed_form import _rel  # noqa: E402
from tests.validation.test_mat20_running_average import W_B2  # noqa: E402
from tests.validation.test_mat20_tuned_heating_dose import T_REF  # noqa: E402
from tests.validation.test_mat21_box_closed_form import A_BOX, B_BOX, C_BOX  # noqa: E402
from tests.validation.test_mat21_box_probe_points_fem import _probe_vals  # noqa: E402
from tests.validation.test_mat21_wall_closed_form import D_WALL, RHO_C_W  # noqa: E402
from tests.validation.test_mat21_wall_duty_closed_form import PERIOD  # noqa: E402
from tests.validation.test_mat21_wall_duty_dose_closed_form import (  # noqa: E402
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
from tests.validation.test_mat21_wall_duty_fem import DT, _census, _exact_duty_t  # noqa: E402
from tests.validation.test_mat21_wall_fem import (  # noqa: E402
    GEL,
    ROBIN,
    WALL,
    _materials,
    _robin_facets,
    _vol_bar,
    _volumes,
    _walled_box,
)
from tests.validation.test_mat21_wall_transient_closed_form import T_STAR  # noqa: E402
from tests.validation.test_mat21_wall_transient_fem import (  # noqa: E402
    MASSLESS,
    N2,
    RUNG,
    _exact_t,
    _lam,
    _rho_c_field,
)

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"  # committed, EX-57
BASENAME = "materials_32_gel_box_wall_thermal_dose"

STEP_4T = 4 * N_MAIN
T_4 = 4.0 * T_STAR
BATCH = np.array([(0.0, 0.0, z) for z in Z_PTS], dtype=float)

#: Item 708's records, real build, `-n 2`:
#: `docs/testing/logs/20261101T104306Z_MAT-21-step23-W1.log`.  Printed beside this run;
#: the two control figures are also the "ceilings" the bars must stay below.
RECORD = {"e_int": 1.291e-03, "bal": 2.25e-13, "miss_mass": 1.0479e-01, "miss_c": 5.063e-01}
CEIL_MASS = 1.05e-01
CEIL_CONT = 5.06e-01


def _say(comm, msg):
    if comm.rank == 0:
        print(msg, flush=True)


def _six_minute(th):
    """``(mx, t_mx)`` of the six-minute window on one history."""
    _a, mx, t_mx = running_average(T_GRID, th, window=float(W_B2))
    return float(mx), float(t_mx)


def main() -> None:
    comm = MPI.COMM_WORLD
    started = time.perf_counter()

    if not np.issubdtype(np.dtype(default_scalar_type), np.complexfloating):
        raise RuntimeError(
            "this example expects the complex DolfinX build the runner's `mat:` group "
            "sources (source /usr/local/bin/dolfinx-complex-mode)"
        )

    _say(comm, "=" * 78)
    _say(comm, "EX-166 -- thermal dose in the walled gel box: CEM43 and the six-minute mean "
               "(MAT-21 step 23 / P50b)")
    _say(comm, "=" * 78)
    _say(comm, f"[gate] imported: INTERIOR_BAND {INTERIOR_BAND:.0e}, E_BAND {E_BAND:.0e}, "
               f"S0_BAND {S0_BAND:g}, DOSE_BAR {DOSE_BAR:g}, AVG_BAR {AVG_BAR:g}, T_REF {T_REF}; "
               f"P = {PERIOD:.4f} s, D = {DUTY}")
    assert T_GRID.size == N2 + 1 == 321, T_GRID.size
    assert abs(T_GRID[1] - DT) <= 1e-12 * DT
    assert abs(T_GRID[STEP_4T] / T_4 - 1.0) <= 1e-14
    assert abs(PERIOD / DT - 20.0) < 1e-12 and abs(DUTY * PERIOD / DT - 10.0) < 1e-12, (PERIOD, DT)
    assert K == K_GEL and H == H_GEL, (K, K_GEL, H, H_GEL)
    assert BATCH.shape == (3, 3) and LABELS[INTERFACE] == "interface"
    assert BATCH[INTERFACE, 2] == C_BOX

    # ---- the body of item 708's `wall_dose` fixture, written out -----------------------
    msh, ns, nw = _walled_box(RUNG, (False, False, True))
    n = msh.topology.index_map(3).size_global
    chi, kf, inter, ct = _materials(msh)
    rcf = _rho_c_field(msh)
    ft = _robin_facets(msh, (A_BOX, B_BOX, C_BOX + D_WALL), (2,))
    vg, vw = _volumes(msh, ct)
    v_gel = 8.0 * A_BOX * B_BOX * C_BOX
    v_wall = 4.0 * A_BOX * B_BOX * 2.0 * D_WALL
    dvg, dvw = abs(vg / v_gel - 1.0), abs(vw / v_wall - 1.0)
    write_setup_figure(
        msh,
        ct,
        FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names={GEL: "gel (k 0.6, rho c 4.18e6), q = 2000 W/m^3 x duty chi(t)",
                      WALL: "acrylic wall, 6 mm (k 0.19, rho c 1.75e6), no source"},
        title="mat:32 -- 1D-reduced walled gel box, duty-cycled source: y = 0 slice "
              "(dose read at x = y = 0, z = 0, 0.9 C, C)",
        slice_normal=(0.0, 1.0, 0.0),
        slice_origin=(0.0, 0.0, 0.0),
        comm=comm,
    )
    _say(comm, f"[mesh] m = {RUNG} cells {tuple(ns)} N_w = {nw}: n = {n}; V_gel rel dev "
               f"{dvg:.3e}, V_wall rel dev {dvw:.3e}, bar {_vol_bar(n):.3e}; dt {DT:.4f} s, "
               f"{N2} CN steps; ranks {comm.size}")
    ts = time.perf_counter()
    res = solve_transient_conduction(
        msh, ft, k=kf, rho_c=rcf, q=Q * chi, h=H_GEL, t_inf=T_INF, t_init=T_INF,
        robin_tag=ROBIN, theta=0.5, startup_half_steps=2, dt=DT, n_steps=N2,
        snapshot_steps=tuple(range(N2 + 1)), source_schedule=_chi(PERIOD, DUTY))
    t_cn = time.perf_counter() - ts
    s0, tot, n_on, n_off = _census(res)
    s0 = float(np.real(s0))
    bal = float(np.real(res.max_balance_residual))
    imag = _max_imag(res.temperature)
    _say(comm, f"[E] CN solve {t_cn:.1f} s; (E) max balance {bal:.3e} (record {RECORD['bal']:.2e}, "
               f"band {E_BAND:g}); (S0) sum chi dt {tot:.6f} s: rel {s0:.3e} (band {S0_BAND:g}); "
               f"{n_on} on / {n_off} off (half-)steps")

    tp = time.perf_counter()
    hist = np.zeros((T_GRID.size, BATCH.shape[0]))
    valid = np.ones(BATCH.shape[0], dtype=bool)
    for k in range(T_GRID.size):
        hist[k], ok = _probe_vals(res.snapshots[k], BATCH)
        valid &= ok
    t_probe = time.perf_counter() - tp

    th = res.snapshots[STEP_4T] - T_INF
    lam_w = _lam(RHO_C_W)
    ex_4 = _exact_duty_t(msh, T_4, RHO_C_W, lam_w)
    e_int = float(np.real(_rel_l2(msh, th, ex_4, inter)))
    miss_c = float(np.real(_rel_l2(msh, th, _exact_t(msh, T_4, RHO_C_W, lam_w), inter)))
    _say(comm, f"[T] 4 t* = {T_4:.3f} s: interior rel L2 vs duty series {e_int:.6e} (record "
               f"{RECORD['e_int']:.3e}, band {INTERIOR_BAND:.0e}); continuous-source series "
               f"miss {miss_c:.6e} (record {RECORD['miss_c']:.3e}); point reads {t_probe:.1f} s; "
               f"all valid {bool(np.all(valid))}")
    series = _duty_hist(RHO_C_W)
    massless = _duty_hist(MASSLESS * RHO_C_W)

    # ---- gates, collected first (item 668) ---------------------------------------------
    hist = np.real(hist)
    gates = {
        "T": e_int <= INTERIOR_BAND,
        "E": bal <= E_BAND,
        "S0": s0 <= S0_BAND,
        "valid": bool(np.all(valid)),
        "continuous_control": miss_c > INTERIOR_BAND,
        "no_imag": imag <= 1e-12,
        "bars_below_ceilings": DOSE_BAR < CEIL_MASS and INTERIOR_BAND < CEIL_CONT,
    }
    for scale, lab in SCALES:
        for j, name in enumerate(LABELS):
            thj = hist[:, j] * scale
            T = T_INF + thj
            dose, ref, rel, _swap, miss = _q_check(T_GRID, T)
            mx, t_mx, lo, hi, wref, wrel = _win_check(T_GRID, thj)
            dose, ref, rel, miss = (float(np.real(x)) for x in (dose, ref, rel, miss))
            wrel = float(np.real(wrel))
            _say(comm, f"(Q-wf) {lab} {name}: FEM CEM43 {dose:.6e} s vs quad {ref:.6e} s rel "
                       f"{rel:.3e} (bar {DOSE_BAR:g}); max T {np.max(T):.6f} C (< {T_REF}); "
                       f"R swap miss {miss:.4e} (> {DOSE_BAR:g})")
            _say(comm, f"(Win-wf) {lab} {name}: max six-minute mean {mx:.6e} K at t_max "
                       f"{t_mx:.4f} s; bracket [{lo:.6e}, {hi:.6e}]; quad {wref:.6e} rel "
                       f"{wrel:.3e} (bar {AVG_BAR:g})")
            key = f"{lab} {name}"
            gates[f"Q {key}"] = rel <= DOSE_BAR
            gates[f"maxT {key}"] = float(np.max(T)) < T_REF
            gates[f"Win {key}"] = wrel <= AVG_BAR
            gates[f"tmax {key}"] = T_GRID[0] + W_B2 <= t_mx <= T_GRID[-1]
            gates[f"bracket {key}"] = lo <= mx <= hi
            gates[f"Rswap {key}"] = miss > DOSE_BAR
    d_fem = float(np.real(cem43(T_GRID, T_INF + hist[:, INTERFACE])))
    d_m = float(np.real(cem43(T_GRID, T_INF + massless[:, INTERFACE])))
    miss_m = float(np.real(_rel(d_m, d_fem)))
    gates["massless_control"] = miss_m > DOSE_BAR
    _say(comm, f"(c-mass) interface dose: massless-wall series misses the FEM by {miss_m:.4e} "
               f"(record {RECORD['miss_mass']:.4e}; ceiling {CEIL_MASS:.2e}; bar {DOSE_BAR:g}, "
               f"separation {miss_m / DOSE_BAR:.3e} x)")
    _say(comm, f"(c-cont) continuous-source series at 4 t* misses by {miss_c:.6e} vs "
               f"INTERIOR_BAND {INTERIOR_BAND:.0e} (ceiling {CEIL_CONT:.2e}; separation "
               f"{miss_c / INTERIOR_BAND:.2f}x)")

    # ---- printed only: FEM vs series dose table ----------------------------------------
    for j, name in enumerate(LABELS):
        fem_h, ser_h = hist[:, j], series[:, j]
        rows = [f"max_n |dtheta| {np.max(np.abs(fem_h - ser_h)):.6e} K"]
        for scale, lab in SCALES:
            df = float(np.real(cem43(T_GRID, T_INF + fem_h * scale)))
            ds = float(np.real(cem43(T_GRID, T_INF + ser_h * scale)))
            mf, tf = _six_minute(fem_h * scale)
            ms, tsr = _six_minute(ser_h * scale)
            rows.append(f"{lab}: dose FEM {df:.6e} s, series {ds:.6e} s (rel "
                        f"{float(np.real(_rel(df, ds))):+.3e}); six-minute FEM {mf:.6f} K at "
                        f"{tf:.2f} s, series {ms:.6f} K at {tsr:.2f} s")
        _say(comm, f"RECORD {name}: " + "; ".join(rows))
    _say(comm, f"max|Im theta| = {imag:.3e}")

    # ---- output ------------------------------------------------------------------------
    V = res.temperature.function_space

    def _real_fn(name, arr):
        f = fem.Function(V, name=name)
        f.x.array[:] = np.real(arr)
        f.x.scatter_forward()
        return f

    th_4 = _real_fn("theta_4_t_star", res.snapshots[STEP_4T].x.array - T_INF)
    th_2 = _real_fn("theta_t2", res.snapshots[N2].x.array - T_INF)
    ser_fn = fem.Function(V)
    ser_fn.interpolate(ex_4)
    ser4 = _real_fn("theta_duty_series_4_t_star", ser_fn.x.array)
    err = _real_fn("theta_error_4_t_star", th_4.x.array - ser4.x.array)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_combined", msh, ct,
        {"theta_4_t_star": th_4, "theta_t2": th_2, "theta_duty_series_4_t_star": ser4,
         "theta_error_4_t_star": err},
        comm=comm,
    )
    adopt_host_ownership(OUTPUT_DIR, comm=comm)
    _say(comm, f"\n[output] {path} -- theta(4 t*), theta(t2), duty series and error at 4 t*")

    red = [k for k, v in gates.items() if not v]
    _say(comm, f"[gates] red: {red or 'none'} of {len(gates)}")
    assert all(gates.values()), red

    _say(comm, f"\n[timing] total {time.perf_counter() - started:.1f} s at -n {comm.size}")
    _say(comm, "\nAll imported gates hold.")


if __name__ == "__main__":
    main()
