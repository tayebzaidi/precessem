"""Example (`EX-162`): the walled gel box heating up -- the acrylic wall's heat capacity
in the transient, against the two-layer slab series.

`mat:28` reads the acrylic-walled gel box in *steady* state, where the wall enters only
through its conductance k_w / d (the h_eff shortcut). In a transient the wall also
carries a heat capacity rho_wc_w that delays the rise. Here the 1D-reduced walled box of
`MAT-21` step 20 (item 670: walls on the z faces only, so the two-layer slab series is
exact for the 3D FEM) is stepped by Crank-Nicolson (theta = 0.5, two backward-Euler
start-up half-steps, dt = t*/N_MAIN, 8 N_MAIN steps) at m = 4 (119 808 cells) with the
source on from t = 0. Snapshots at t* = 0.2 T_C and t2 = 1.6 T_C are compared with
``two_layer_slab_transient_temperature`` (400 modes), and the same FEM field is then
compared with the *massless-wall* series (rho_wc_w x 1e-9), the control.

**It asserts, and it does not re-implement** (the `ANS-1` rule). Every solver, series,
band and helper is imported by name from the gate chain of
``tests/validation/test_mat21_wall_transient_fem.py`` (`MAT-21` step 20c, item 682). That
module's ``walled`` fixture is a pytest fixture and is not called: its body is written
out here (disclosed, as `EX-159` did for item 670's fixtures), with ``np.real`` on each
integral because the `mat:` group sources the complex build. No gate edit, no ``src/``
edit.

Asserted (bands imported): (V) gel and wall tag volumes within ``_vol_bar(n)``; (E) max
per-step balance <= ``E_BAND``; (T) gel-interior rel L2 <= ``INTERIOR_BAND`` at both t*
and t2; (c-m) the massless-wall series misses the FEM at t2 by more than ``CTL_FACTOR``
times the FEM's own interior error there.
Printed only: the centre theta at both snapshots -- FEM, two-layer series, massless-wall
series, and their relative difference (the size of the wall-heat-capacity effect a
phantom builder who used h_eff would ignore).

Run it::

    ./run_examples.sh -e mat:30 -n 2 -t 480

Output: ``examples/materials/paraview_output/
materials_30_gel_box_wall_heat_capacity_transient_combined.xdmf`` (theta(t*), theta(t2),
the two-layer series at t2 and the error, on the m = 4 walled box).
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
from tests.validation.test_mat21_wall_closed_form import D_WALL, K_W, RHO_C_W  # noqa: E402
from tests.validation.test_mat21_wall_fem import (  # noqa: E402
    BASE,
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
from tests.validation.test_mat21_wall_transient_closed_form import T_STAR  # noqa: E402
from tests.validation.test_mat21_wall_transient_fem import (  # noqa: E402
    CTL_FACTOR,
    MASSLESS,
    N2,
    RUNG,
    T2,
    _exact_t,
    _lam,
    _rho_c_field,
    _series,
)

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"  # committed, EX-57
BASENAME = "materials_30_gel_box_wall_heat_capacity_transient"

#: Item 682's records, real build, `-n 2`:
#: `20261028T021107Z_MAT-21-step20c-W1-fem.log:52-57`.  Printed beside this run, a
#: difference reported and never asserted.
RECORD = {
    "dvg": 2.507e-13, "dvw": 2.724e-13, "balance": 2.350e-13,
    "int_tstar": 2.244928e-3, "int_t2": 1.483557e-3, "miss0_t2": 2.059306e-2,
    "c_fem_tstar": 1.330815144, "c_ser_tstar": 1.327311508, "c_0_tstar": 1.338774894,
    "c_fem_t2": 7.958737121, "c_ser_t2": 7.945611514, "c_0_t2": 8.090989721,
}
#: `mat:28`'s bare-gel-box-at-h_eff steady centre record (K) -- steady, not a t* rise.
MAT28_BARE_HEFF_STEADY = 14.6001


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
    _say(comm, "EX-162 -- walled gel box heating up: the wall's rho c in the transient "
               "(MAT-21 step 20c / P48c)")
    _say(comm, "=" * 78)
    _say(comm, f"[gate] imported: INTERIOR_BAND {INTERIOR_BAND:.0e}, E_BAND {E_BAND:.0e}, "
               f"CTL_FACTOR {CTL_FACTOR:.0f}; (A, B, C) = {(A_BOX, B_BOX, C_BOX)}, d = {D_WALL}")
    assert K == K_GEL and H == H_GEL, (K, K_GEL, H, H_GEL)

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
        region_names={GEL: "gel (k 0.6, rho c 4.18e6), q = 2000 W/m^3",
                      WALL: "acrylic wall, 6 mm (k 0.19, rho c 1.75e6), no source"},
        title="mat:30 -- 1D-reduced walled gel box, Robin h on the outer z faces: y = 0 slice",
        slice_normal=(0.0, 1.0, 0.0),
        slice_origin=(0.0, 0.0, 0.0),
        comm=comm,
    )
    _say(comm, f"[mesh] m = {RUNG} cells {tuple(ns)} N_w = {nw} (BASE {BASE}, tags gel {GEL} / "
               f"wall {WALL}): n = {n}; V_gel rel dev {dvg:.3e} (record {RECORD['dvg']:.3e}), "
               f"V_wall rel dev {dvw:.3e} (record {RECORD['dvw']:.3e}), bar {_vol_bar(n):.3e}; "
               f"T_C {T_C:.3f} s, t* {T_STAR:.3f} s, t2 {T2:.3f} s = {T2 / T_C:.4f} T_C; "
               f"dt {T_STAR / N_MAIN:.4f} s, {N2} CN steps; ranks {comm.size}")

    # ---- CN solve ----------------------------------------------------------------------
    ts = time.perf_counter()
    res = solve_transient_conduction(
        msh, ft, k=kf, rho_c=rcf, q=Q * chi, h=H_GEL, t_inf=T_INF, t_init=T_INF,
        robin_tag=ROBIN, theta=0.5, startup_half_steps=2, dt=T_STAR / N_MAIN, n_steps=N2,
        snapshot_steps=(N_MAIN, N2))
    t_solve = time.perf_counter() - ts
    src = float(np.real(res.source_power))
    bal = float(np.real(res.max_balance_residual))
    nb = len(res.balance_residuals)
    imag = _max_imag(res.temperature)
    _say(comm, f"[E] CN solve {t_solve:.1f} s; source power {src:.9e} W (Q V_gel {Q * v_gel:.9e}); "
               f"max balance {bal:.3e} over {nb} (half-)steps (record {RECORD['balance']:.3e}, "
               f"bound {E_BAND:.0e})")

    # ---- snapshots against the series --------------------------------------------------
    lam_w = _lam(RHO_C_W)
    lam_0 = _lam(MASSLESS * RHO_C_W)
    snap = {}
    last = {}
    for step, t, tag in ((N_MAIN, T_STAR, "tstar"), (N2, T2, "t2")):
        th = res.snapshots[step] - T_INF
        ex = _exact_t(msh, t, RHO_C_W, lam_w)
        ex0 = _exact_t(msh, t, MASSLESS * RHO_C_W, lam_0)
        e_int = _rel(msh, th, ex, inter)
        e_whole = _rel(msh, th, ex)
        miss0 = _rel(msh, th, ex0, inter)
        c_h = _centre(res.snapshots[step]) - T_INF
        c_ex = float(np.real(_series(0.0, t, RHO_C_W, lam_w)[0]))
        c_0 = float(np.real(_series(0.0, t, MASSLESS * RHO_C_W, lam_0)[0]))
        snap[step] = {"e_int": e_int, "miss0": miss0}
        last[step] = ex
        _say(comm, f"[T] step {step} t = {t:.3f} s ({t / T_C:.4f} T_C): interior rel L2 vs "
                   f"two-layer series {e_int:.6e} (record {RECORD['int_' + tag]:.6e}, bound "
                   f"{INTERIOR_BAND:.0e}); whole {e_whole:.6e} (printed); massless-wall series "
                   f"vs FEM, interior {miss0:.6e}")
        _say(comm, f"[centre] step {step} theta: FEM {c_h:.9f} K (record {RECORD['c_fem_' + tag]}), "
                   f"two-layer series {c_ex:.9f} K (record {RECORD['c_ser_' + tag]}; rel "
                   f"{c_h / c_ex - 1.0:+.3e}), massless-wall series {c_0:.9f} K (record "
                   f"{RECORD['c_0_' + tag]}); wall-heat-capacity effect (massless / two-layer - 1) "
                   f"{c_0 / c_ex - 1.0:+.4e}")
    s2 = snap[N2]
    _say(comm, f"[c-m] at t2: massless-wall series misses the FEM by {s2['miss0']:.6e} "
               f"(record {RECORD['miss0_t2']:.6e}) vs {CTL_FACTOR:.0f} x FEM interior error "
               f"{CTL_FACTOR * s2['e_int']:.6e} (separation {s2['miss0'] / s2['e_int']:.2f}x; "
               f"record 13.9x)")
    _say(comm, f"[orientation] mat:28's bare gel box at h_eff, STEADY centre record "
               f"{MAT28_BARE_HEFF_STEADY} K (steady, not a t* rise; no bare transient series "
               f"exists in the gate chain)")
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
    ser_fn = fem.Function(V)
    ser_fn.interpolate(last[N2])
    ser2 = _real_fn("theta_series_t2", ser_fn.x.array)
    err = _real_fn("theta_error_t2", th_2.x.array - ser2.x.array)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_combined", msh, ct,
        {"theta_t_star": th_s, "theta_t2": th_2, "theta_series_t2": ser2,
         "theta_error_t2": err},
        comm=comm,
    )
    adopt_host_ownership(OUTPUT_DIR, comm=comm)
    _say(comm, f"\n[output] {path} -- theta(t*), theta(t2), two-layer series at t2 and error")

    # ---- assertions (the gate's) -------------------------------------------------------
    assert dvg <= _vol_bar(n), f"(V) gel {dvg:.3e}"
    assert dvw <= _vol_bar(n), f"(V) wall {dvw:.3e}"
    assert bal <= E_BAND, f"(E) {bal:.3e} > {E_BAND:.0e}"
    errs = {s: v["e_int"] for s, v in snap.items()}
    assert all(e <= INTERIOR_BAND for e in errs.values()), f"(T) {errs}"
    assert s2["miss0"] > CTL_FACTOR * s2["e_int"], f"(c-m) {s2}"
    assert imag <= 1e-12, f"complex build left Im theta = {imag:.3e}"

    _say(comm, f"\n[timing] total {time.perf_counter() - started:.1f} s at -n {comm.size}")
    _say(comm, "\nAll imported gates hold.")


if __name__ == "__main__":
    main()
