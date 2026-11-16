"""Example (`EX-183`): the drift through a wall with heat capacity -- the composite
comparand against the FEM, on both pre-equilibration schedules.

`mat:37` shows the drift inside the walled box against P62's massless-wall comparand
(``_drift_1d``). Here item 836's fixture (item 682's 1D-reduced walled box, m = 4, 119 808
cells, Crank-Nicolson with two backward-Euler half-steps) is run three times:

- **R0**, from theta_0 = 0, source on from t = 0, to t*;
- **R1**, t_w = 0: theta_0 = g z on every cell, to t*;
- **R2**, t_w = T_C: theta_0 = g z, source off on [0, T_C), to T_C + t*.

The comparand is item 682's series ``s_w`` plus P67's two-layer composite drift
``_drift_cap_z`` (the wall's heat capacity carried in the drift) at t* / T_C + t*. The
massless (h_eff) comparand is printed beside it.

**It asserts, and it does not re-implement** (the `ANS-1` rule). Bands, builders and
helpers are imported by name from ``tests/validation/test_mat21_wall_capacity_drift_fem.py``
(item 836) and the modules it imports from. That gate's ``runs`` fixture is not called: its
body is written out here, plus ``np.real`` on each integral (the group sources the complex
build). No ``test_*`` name is imported. No gate edit, no ``src/`` edit.

Asserted (collected first, then asserted): (I) R1 and R2 interior rel L2 against the
composite <= ``INTERIOR_BAND``; (E) <= ``E_BAND`` on all three runs, ``source_power`` > 0 and
identical; (S0) on R2 <= ``S0_BAND``; the unique-z precondition <= ``PRE_BAR``; R0's control
reproduces item 682's record within ``REPRO_BOUND``; (c-ser) the series alone (no drift)
misses R1 by more than ``INTERIOR_BAND`` (size *predicted* about 0.17 from the P69 pointwise
ratio, not logged before this run); max|Im theta| = 0 on all three runs.
Printed only: the massless comparand's (I) beside the composite's, item 836's table at P5,
``ADJ_PT``, P0, and the composite / FEM-floor ratios.

The wall's k and rho_c are nominal PMMA values, citation owed (W3). One rung, one time step,
this fixture's g. Not claimed: the six-wall box, a rate, a recommended wait, a cited PMMA
rho_c, the D = 0.25 dose.

Run it::

    ./run_examples.sh -e mat:38 -n 2 -t 480

Output: ``examples/materials/paraview_output/
materials_38_gel_box_wall_capacity_drift_combined.xdmf`` (theta at t* on R1 and R1's
composite comparand).
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

from tests.validation.test_imp7_duty_cycle_fem import S0_BAND  # noqa: E402
from tests.validation.test_mat12_transient_sar_heating_sphere import E_BAND  # noqa: E402
from tests.validation.test_mat13_transient_perfusion_sphere import INTERIOR_BAND  # noqa: E402
from tests.validation.test_mat15_coil_driven_heating_phantom import (  # noqa: E402
    H_GEL,
    K_GEL,
    N_MAIN,
)
from tests.validation.test_mat15_finite_cylinder_closed_form import H, K, Q, T_INF  # noqa: E402
from tests.validation.test_mat15_finite_cylinder_fem import _max_imag, _rel  # noqa: E402
from tests.validation.test_mat21_box_closed_form import A_BOX, B_BOX, C_BOX, T_C  # noqa: E402
from tests.validation.test_mat21_initial_gradient_closed_form import G_GRAD, GC  # noqa: E402
from tests.validation.test_mat21_wall_capacity_drift_closed_form import _drift_cap  # noqa: E402
from tests.validation.test_mat21_wall_capacity_drift_fem import _cap_z  # noqa: E402
from tests.validation.test_mat21_wall_closed_form import D_WALL, RHO_C_W  # noqa: E402
from tests.validation.test_mat21_wall_fem import (  # noqa: E402
    GEL,
    ROBIN,
    WALL,
    _materials,
    _robin_facets,
    _walled_box,
)
from tests.validation.test_mat21_wall_pre_equilibration_fem import (  # noqa: E402
    DT,
    I_ADJ,
    I_P0,
    I_P5,
    ITEM_682_TSTAR,
    N_OFF,
    N_R2,
    PRE_BAR,
    READ_PTS,
    REPRO_BOUND,
    _census,
    _chi_off_on,
    _drift_z,
    _read,
    _series_z,
    _ZField,
)
from tests.validation.test_mat21_wall_transient_closed_form import T_STAR  # noqa: E402
from tests.validation.test_mat21_wall_transient_fem import RUNG, _rho_c_field  # noqa: E402

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"  # committed, EX-57
BASENAME = "materials_38_gel_box_wall_capacity_drift"

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
    _say(comm, "EX-183 -- the drift through a wall with heat capacity (MAT-21 step 38 / P67-P69)")
    _say(comm, "=" * 78)
    _say(comm, f"[gate] imported: INTERIOR_BAND {INTERIOR_BAND:.0e}, E_BAND {E_BAND:.0e}, "
               f"S0_BAND {S0_BAND:g}, PRE_BAR {PRE_BAR:.0e}, REPRO_BOUND {REPRO_BOUND:g}; "
               f"g {G_GRAD} K/m, g C {GC} K")
    assert K == K_GEL and H == H_GEL, (K, K_GEL, H, H_GEL)

    msh, ns, nw = _walled_box(RUNG, (False, False, True))
    chi, kf, inter, ct = _materials(msh)
    rcf = _rho_c_field(msh)
    ft = _robin_facets(msh, (A_BOX, B_BOX, C_BOX + D_WALL), (2,))
    n = msh.topology.index_map(3).size_global
    write_setup_figure(
        msh,
        ct,
        FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names={GEL: "gel, q on, theta_0 = g z", WALL: "acrylic wall (nominal)"},
        title="mat:38 -- gel box with acrylic walls on z = +/-C, Robin on the outer wall "
              "faces, T(0) = T_inf + g z: y = 0 slice",
        slice_normal=(0.0, 1.0, 0.0),
        slice_origin=(0.0, 0.0, 0.0),
        comm=comm,
    )
    _say(comm, f"[mesh] walled box m = {RUNG} cells {tuple(ns)} N_w = {nw}: n = {n}; "
               f"T_C {T_C:.3f} s, t* {T_STAR:.3f} s, dt {DT:.4f} s; R2 {N_R2} steps "
               f"(off {N_OFF}); ranks {comm.size}")

    grad = lambda x: T_INF + G_GRAD * x[2]  # noqa: E731  (every cell, gel and wall)
    kw = dict(k=kf, q=Q * chi, h=H_GEL, t_inf=T_INF, t_init=T_INF, robin_tag=ROBIN, theta=0.5,
              startup_half_steps=2, dt=DT)
    specs = {"R0": dict(rho_c=rcf, n_steps=N_MAIN),
             "R1": dict(rho_c=rcf, n_steps=N_MAIN, t_init_expr=grad),
             "R2": dict(rho_c=rcf, n_steps=N_R2, t_init_expr=grad,
                        source_schedule=_chi_off_on)}
    E, src, th_pt, ok_pt, res, wall = {}, {}, {}, {}, {}, {}
    for key, sp in specs.items():
        tq = time.perf_counter()
        res[key] = solve_transient_conduction(msh, ft, **kw, **sp)
        wall[key] = time.perf_counter() - tq
        E[key] = float(np.real(res[key].max_balance_residual))
        src[key] = float(np.real(res[key].source_power))
        th_pt[key], ok_pt[key] = _read(res[key].temperature, READ_PTS)
        _say(comm, f"[run {key}] {wall[key]:.1f} s; (E) {E[key]:.3e} (E_BAND {E_BAND:.0e}); "
                   f"source_power {src[key]:.9e} W")
    imag = max(_max_imag(r.temperature) for r in res.values())
    s0, tot = _census(res["R2"].source_multipliers)
    _say(comm, f"[S0] R2 sum chi dt {tot:.9f} s vs t* {T_STAR:.9f} s rel {s0:.3e} "
               f"(band {S0_BAND:g})")

    tq = time.perf_counter()
    zf = _ZField()
    s_w = zf(msh, _series_z(RHO_C_W), T_STAR)
    c_0 = zf(msh, _cap_z, T_STAR)
    c_c = zf(msh, _cap_z, T_C + T_STAR)
    d_0 = zf(msh, _drift_z, T_STAR)
    d_c = zf(msh, _drift_z, T_C + T_STAR)
    pre = comm.allreduce(zf.worst, op=MPI.MAX)
    _say(comm, f"[series] five z-only comparands on the unique dof z ({zf.n_unique} values per "
               f"rank max): {time.perf_counter() - tq:.1f} s; precondition worst rel "
               f"{pre:.3e} (bar {PRE_BAR:.0e})")

    th = {key: r.temperature - T_INF for key, r in res.items()}
    ref = _rel(msh, th["R0"], s_w, inter)
    I1 = _rel(msh, th["R1"], s_w + c_0, inter)
    I2 = _rel(msh, th["R2"], s_w + c_c, inter)
    h1 = _rel(msh, th["R1"], s_w + d_0, inter)
    h2 = _rel(msh, th["R2"], s_w + d_c, inter)
    cser = _rel(msh, th["R1"], s_w, inter)
    _say(comm, f"[I] vs the composite: R1 interior rel L2 {I1:.6e}, R2 {I2:.6e} (bound "
               f"{INTERIOR_BAND:.0e}; item 836 records 2.639057e-03 / 2.258920e-03); whole "
               f"{_rel(msh, th['R1'], s_w + c_0):.6e} / {_rel(msh, th['R2'], s_w + c_c):.6e} "
               f"(printed)")
    _say(comm, f"[massless, printed] h_eff comparand (I): R1 {h1:.6e}, R2 {h2:.6e} "
               f"(item 836 records 8.105518e-03 / 2.415575e-03)")
    _say(comm, f"[control] R0 interior {ref:.6e} vs item 682 {ITEM_682_TSTAR:.6e}, |delta| "
               f"{abs(ref - ITEM_682_TSTAR):.3e} (bound {REPRO_BOUND:g})")
    _say(comm, f"[c-ser] series alone (no drift) misses R1 by {cser:.6e} (asserted > "
               f"{INTERIOR_BAND:.0e}; size PREDICTED about 0.17 from the P69 pointwise ratio, "
               f"no prior log)")
    _say(comm, f"[ceiling, printed] FEM floor R0 {ref:.6e}; composite/floor R1 "
               f"{I1 / ref:.3f}x, R2 {I2 / ref:.3f}x")

    _say(comm, "[TABLE] run | point | FEM drift (Rk - R0) K | composite K | FEM - composite K")
    for key, t_w in (("R1", 0.0), ("R2", T_C)):
        for i, lab in ((I_P5, "P5"), (I_ADJ, "ADJ_PT"), (I_P0, "P0")):
            p = tuple(READ_PTS[i])
            dh = th_pt[key][i] - th_pt["R0"][i]
            cap = _drift_cap(p, T_STAR + t_w)
            _say(comm, f"[TABLE] {key} | {lab} | {dh:+.6e} | {cap:+.6e} | {dh - cap:+.3e}")
    _say(comm, f"max|Im theta| = {imag:.3e}")

    V = res["R0"].temperature.function_space
    outs = {}
    f = fem.Function(V, name="theta_R1_t_star")
    f.x.array[:] = np.real(res["R1"].temperature.x.array) - T_INF
    f.x.scatter_forward()
    outs[f.name] = f
    comp = fem.Function(V, name="theta_comparand_R1")
    d_p1 = fem.Function(V)
    comp.interpolate(s_w)
    d_p1.interpolate(c_0)
    comp.x.array[:] = comp.x.array + d_p1.x.array
    comp.x.scatter_forward()
    outs[comp.name] = comp
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_combined", msh, ct, outs, comm=comm)
    adopt_host_ownership(OUTPUT_DIR, comm=comm)
    _say(comm, f"\n[output] {path} -- theta(t*) on R1 and R1's composite comparand")

    gates = {
        "valid_points": all(all(o) for o in ok_pt.values()),
        "precondition_unique_z": pre <= PRE_BAR,
        "I_R1_composite": I1 <= INTERIOR_BAND,
        "I_R2_composite": I2 <= INTERIOR_BAND,
        "E_source_power_identical": src["R0"] > 0.0 and len(set(src.values())) == 1,
        "S0_R2": s0 <= S0_BAND,
        "control_repro_682": abs(ref - ITEM_682_TSTAR) <= REPRO_BOUND,
        "control_series_alone_misses": cser > INTERIOR_BAND,
        "max_imag_zero": imag == 0.0,
    }
    for key, e in E.items():
        gates[f"E_{key}"] = e <= E_BAND
    _say(comm, f"[gates] {gates}")
    _say(comm, f"\n[timing] total {time.perf_counter() - started:.1f} s at -n {comm.size}")
    assert all(gates.values()), gates
    _say(comm, "\nAll imported gates hold.")


if __name__ == "__main__":
    main()


