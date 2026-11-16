"""Example (`EX-164`): the walled gel box under a duty-cycled source -- the acrylic wall's
heat capacity when the heating switches on and off, against the duty-cycled two-layer slab
series.

`mat:30` steps the 1D-reduced walled box of `MAT-21` step 20 with the source on from
t = 0. Here the *same* source is switched on and off: period P = 0.5 t* (20 CN steps), duty
D = 0.5 (10 steps on), through ``source_schedule`` (chi is evaluated at each (half-)step
midpoint, so a period of exactly 20 steps makes the on-time census exact). Crank-Nicolson
(theta = 0.5, two backward-Euler start-up half-steps, dt = t*/N_MAIN, 8 N_MAIN steps) at
m = 4 (119 808 cells); snapshots at t* and t2 = 8 t* are compared with item 689's
``two_layer_slab_duty_cycle_temperature`` (400 modes -- the step response superposed phase
by phase). The control is the *continuous-source* series, which is wrong by about 50 % at
t*.

**It asserts, and it does not re-implement** (the `ANS-1` rule). Every solver, series,
band and helper is imported by name from the gate chain of
``tests/validation/test_mat21_wall_duty_fem.py`` (`MAT-21` step 21b, item 688); ``PERIOD``
is item 689's (0.5 x MAT-21's t*), not item 675's. That module's ``walled`` fixture is a
pytest fixture and is not called: its body is written out here (disclosed, as `EX-162` did
for item 682's), with ``np.real`` on each integral because the `mat:` group sources the
complex build. No gate edit, no ``src/`` edit. One import beyond the brief's list:
``_series`` from item 682's module (the continuous-source series' centre value, printed).

Asserted (bands imported): (Sch) P = 20 dt and D P = 10 dt to 1e-12; (V) gel and wall tag
volumes within ``_vol_bar(n)``; (E) max per-step balance <= ``E_BAND``; (S0) the census
sum chi dt against D t2 <= ``S0_BAND``; (T) gel-interior rel L2 against the duty series
<= ``INTERIOR_BAND`` at both t* and t2; (c-c) the continuous-source series misses the FEM
at t* by more than ``CTL_FACTOR`` x the FEM's interior error there.
Printed only: the centre theta at both snapshots -- FEM, duty series, continuous-source
series, massless-wall duty series, and the relative difference between the duty series
with and without the wall's heat capacity.

Run it::

    ./run_examples.sh -e mat:31 -n 2 -t 480

Output: ``examples/materials/paraview_output/
materials_31_gel_box_wall_duty_cycled_transient_combined.xdmf`` (theta(t*), theta(t2), the
duty series at t* and t2 and the error at t2, on the m = 4 walled box).
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
from fem_em_solver.utils.analytical import two_layer_slab_duty_cycle_temperature  # noqa: E402,F401

from tests.validation.test_imp7_duty_cycle_fem import DUTY, S0_BAND, _chi  # noqa: E402
from tests.validation.test_mat12_transient_sar_heating_sphere import E_BAND  # noqa: E402
from tests.validation.test_mat13_transient_perfusion_sphere import INTERIOR_BAND  # noqa: E402
from tests.validation.test_mat15_coil_driven_heating_phantom import (  # noqa: E402
    H_GEL,
    K_GEL,
    N_MAIN,
)
from tests.validation.test_mat15_finite_cylinder_closed_form import (  # noqa: E402
    H,
    K,
    Q,
    RHO_C,
    T_INF,
)
from tests.validation.test_mat15_finite_cylinder_fem import _max_imag, _rel  # noqa: E402
from tests.validation.test_mat21_box_closed_form import A_BOX, B_BOX, C_BOX, T_C  # noqa: E402
from tests.validation.test_mat21_wall_closed_form import D_WALL, K_W, RHO_C_W  # noqa: E402,F401
from tests.validation.test_mat21_wall_duty_closed_form import PERIOD  # noqa: E402
from tests.validation.test_mat21_wall_duty_fem import (  # noqa: E402
    DT,
    ON_TIME,
    _census,
    _duty_series,
    _exact_duty_t,
)
from tests.validation.test_mat21_wall_fem import (  # noqa: E402
    GEL,
    ROBIN,
    WALL,
    _centre,
    _materials,
    _robin_facets,
    _vol_bar,
    _volumes,
    _walled_box,
)
from tests.validation.test_mat21_wall_transient_closed_form import N_SERIES, T_STAR  # noqa: E402,F401
from tests.validation.test_mat21_wall_transient_fem import (  # noqa: E402
    CTL_FACTOR,
    INTERP_CHUNK,  # noqa: F401
    MASSLESS,
    N2,
    RUNG,
    SLAB,  # noqa: F401
    T2,
    _exact_t,
    _lam,
    _rho_c_field,
    _series,
)

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"  # committed, EX-57
BASENAME = "materials_31_gel_box_wall_duty_cycled_transient"

#: Item 688's records, real build, `-n 2`: `20261028T172612Z_MAT-21-step21b-W1.log:52-60`.
#: Printed beside this run, a difference reported and never asserted.
RECORD = {
    "dvg": 2.507e-13, "dvw": 2.724e-13, "balance": 2.269e-13, "s0": 0.0,
    "int_tstar": 1.697155e-3, "int_t2": 1.447368e-3, "miss_c_tstar": 5.066571e-1,
    "miss0_tstar": 3.077e-2, "miss0_t2": 2.064e-2,
    "c_fem_tstar": 0.661673131, "c_ex_tstar": 0.659441298,
    "c_fem_t2": 3.937390665, "c_ex_t2": 3.930539409,
}


def _say(comm, msg):
    if comm.rank == 0:
        print(msg, flush=True)


def main() -> None:
    comm = MPI.COMM_WORLD
    started = time.perf_counter()

    if not np.issubdtype(np.dtype(default_scalar_type), np.complexfloating):
        raise RuntimeError(
            "this example expects the complex DolfinX build the runner's `mat:` group "
            "sources (source /usr/local/bin/dolfinx-complex-mode)"
        )

    _say(comm, "=" * 78)
    _say(comm, "EX-164 -- walled gel box under a duty-cycled source: the wall's rho c "
               "(MAT-21 step 21b / P49b)")
    _say(comm, "=" * 78)
    _say(comm, f"[gate] imported: INTERIOR_BAND {INTERIOR_BAND:.0e}, E_BAND {E_BAND:.0e}, "
               f"S0_BAND {S0_BAND:g}, CTL_FACTOR {CTL_FACTOR:.0f}; "
               f"(A, B, C) = {(A_BOX, B_BOX, C_BOX)}, d = {D_WALL}; P = {PERIOD:.4f} s, D = {DUTY}")
    assert K == K_GEL and H == H_GEL, (K, K_GEL, H, H_GEL)
    sch_ok = abs(PERIOD / DT - 20.0) < 1e-12 and abs(DUTY * PERIOD / DT - 10.0) < 1e-12

    # ---- mesh: item 670's 1D-reduced walled box, m = RUNG ------------------------------
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
        title="mat:31 -- 1D-reduced walled gel box, duty-cycled source, Robin h on the outer "
              "z faces: y = 0 slice",
        slice_normal=(0.0, 1.0, 0.0),
        slice_origin=(0.0, 0.0, 0.0),
        comm=comm,
    )
    _say(comm, f"[mesh] m = {RUNG} cells {tuple(ns)} N_w = {nw}: n = {n}; V_gel rel dev "
               f"{dvg:.3e} (record {RECORD['dvg']:.3e}), V_wall rel dev {dvw:.3e} (record "
               f"{RECORD['dvw']:.3e}), bar {_vol_bar(n):.3e}; T_C {T_C:.3f} s, t* {T_STAR:.3f} s, "
               f"t2 {T2:.3f} s; P {PERIOD:.4f} s = {PERIOD / DT:.1f} steps, D {DUTY}; "
               f"dt {DT:.4f} s, {N2} CN steps; ranks {comm.size}")
    _say(comm, f"[Sch] P = 20 dt and D P = 10 dt to 1e-12: {sch_ok}")

    # ---- CN solve with the duty schedule -----------------------------------------------
    ts = time.perf_counter()
    res = solve_transient_conduction(
        msh, ft, k=kf, rho_c=rcf, q=Q * chi, h=H_GEL, t_inf=T_INF, t_init=T_INF,
        robin_tag=ROBIN, theta=0.5, startup_half_steps=2, dt=DT, n_steps=N2,
        snapshot_steps=(N_MAIN, N2), source_schedule=_chi(PERIOD, DUTY))
    t_solve = time.perf_counter() - ts
    bal = float(np.real(res.max_balance_residual))
    nb = len(res.balance_residuals)
    imag = _max_imag(res.temperature)
    _say(comm, f"[E] CN solve {t_solve:.1f} s; max balance {bal:.3e} over {nb} (half-)steps "
               f"(record {RECORD['balance']:.3e}, bound {E_BAND:.0e})")
    s0, tot, n_on, n_off = _census(res)
    s0 = float(np.real(s0))
    _say(comm, f"[S0] census sum chi dt = {tot:.6f} s vs D t2 = {ON_TIME:.6f} s: rel {s0:.3e} "
               f"(record {RECORD['s0']:.0e}, band {S0_BAND:g}); {n_on} on / {n_off} off "
               f"(half-)steps (record 161 / 160)")

    # ---- snapshots against the series --------------------------------------------------
    lam_w = _lam(RHO_C_W)
    lam_0 = _lam(MASSLESS * RHO_C_W)
    snap = {}
    last = {}
    for step, t, tag in ((N_MAIN, T_STAR, "tstar"), (N2, T2, "t2")):
        th = res.snapshots[step] - T_INF
        ex = _exact_duty_t(msh, t, RHO_C_W, lam_w)
        e_int = float(np.real(_rel(msh, th, ex, inter)))
        e_whole = float(np.real(_rel(msh, th, ex)))
        exc = _exact_t(msh, t, RHO_C_W, lam_w)  # continuous source (negative control)
        miss_c = float(np.real(_rel(msh, th, exc, inter)))
        ex0 = _exact_duty_t(msh, t, MASSLESS * RHO_C_W, lam_0)
        miss0 = float(np.real(_rel(msh, th, ex0, inter)))
        c_h = float(np.real(_centre(res.snapshots[step]) - T_INF))
        c_ex = float(np.real(np.atleast_1d(_duty_series(0.0, t, RHO_C_W, lam_w))[0]))
        c_c = float(np.real(_series(0.0, t, RHO_C_W, lam_w)[0]))
        c_0 = float(np.real(np.atleast_1d(_duty_series(0.0, t, MASSLESS * RHO_C_W, lam_0))[0]))
        snap[step] = {"e_int": e_int, "miss_c": miss_c, "miss0": miss0}
        last[step] = ex
        rec_miss0 = RECORD["miss0_" + tag]
        _say(comm, f"[T] step {step} t = {t:.3f} s ({t / T_C:.4f} T_C): interior rel L2 vs duty "
                   f"series {e_int:.6e} (record {RECORD['int_' + tag]:.6e}, bound "
                   f"{INTERIOR_BAND:.0e}); whole {e_whole:.6e} (printed); continuous-source "
                   f"series miss {miss_c:.6e}; massless-wall duty series miss {miss0:.6e} "
                   f"(printed; record {rec_miss0:.3e})")
        _say(comm, f"[centre] step {step} theta: FEM {c_h:.9f} K (record {RECORD['c_fem_' + tag]}), "
                   f"duty series {c_ex:.9f} K (record {RECORD['c_ex_' + tag]}; rel "
                   f"{c_h / c_ex - 1.0:+.3e}), continuous-source series {c_c:.9f} K "
                   f"(rel {c_c / c_ex - 1.0:+.3e}), massless-wall duty series {c_0:.9f} K; "
                   f"wall-heat-capacity effect (massless / two-layer - 1) {c_0 / c_ex - 1.0:+.4e}")
    s1 = snap[N_MAIN]
    _say(comm, f"[c-c] at t*: continuous-source series misses the FEM by {s1['miss_c']:.6e} "
               f"(record {RECORD['miss_c_tstar']:.6e}) vs {CTL_FACTOR:.0f} x FEM interior error "
               f"{CTL_FACTOR * s1['e_int']:.6e} (separation {s1['miss_c'] / s1['e_int']:.2f}x; "
               f"record 298.5x)")
    _say(comm, f"max|Im theta| = {imag:.3e}")

    # ---- output ------------------------------------------------------------------------
    V = res.temperature.function_space

    def _real_fn(name, arr):
        f = fem.Function(V, name=name)
        f.x.array[:] = np.real(arr)
        f.x.scatter_forward()
        return f

    th_s = _real_fn("theta_t_star", res.snapshots[N_MAIN].x.array - T_INF)
    th_2 = _real_fn("theta_t2", res.snapshots[N2].x.array - T_INF)
    def _series_fn(name, g):  # the gates' P2 series, interpolated onto the output space
        f = fem.Function(V, name=name)
        f.interpolate(g)
        f.x.scatter_forward()
        return f

    ser_s = _series_fn("theta_duty_series_t_star", last[N_MAIN])
    ser_2 = _series_fn("theta_duty_series_t2", last[N2])
    err = _real_fn("theta_error_t2", th_2.x.array - ser_2.x.array)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_combined", msh, ct,
        {"theta_t_star": th_s, "theta_t2": th_2, "theta_duty_series_t_star": ser_s,
         "theta_duty_series_t2": ser_2, "theta_error_t2": err},
        comm=comm,
    )
    adopt_host_ownership(OUTPUT_DIR, comm=comm)
    _say(comm, f"\n[output] {path} -- theta(t*), theta(t2), duty series at both, error at t2")

    # ---- assertions (the gate's), all collected first (item 668) -----------------------
    gates = {
        "Sch": bool(sch_ok),
        "vol_gel": dvg <= _vol_bar(n),
        "vol_wall": dvw <= _vol_bar(n),
        "balance": bal <= E_BAND,
        "S0": s0 <= S0_BAND,
        "interior_tstar": snap[N_MAIN]["e_int"] <= INTERIOR_BAND,
        "interior_t2": snap[N2]["e_int"] <= INTERIOR_BAND,
        "control_continuous_tstar": s1["miss_c"] > CTL_FACTOR * s1["e_int"],
        "no_imag": imag <= 1e-12,
    }
    _say(comm, f"[gates] {gates}")
    assert all(gates.values()), gates

    _say(comm, f"\n[timing] total {time.perf_counter() - started:.1f} s at -n {comm.size}")
    _say(comm, "\nAll imported gates hold.")


if __name__ == "__main__":
    main()
