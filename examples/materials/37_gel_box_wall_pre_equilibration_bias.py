"""Example (`EX-180`): the skipped wait inside the walled box -- what the wall's heat
capacity adds to the drift.

`mat:33` shows the skipped wait in the wall-free box. `mat:28`, `mat:30` and `mat:31` show
the walled box with the source alone. None shows the leftover gradient inside the wall.
Here item 682's 1D-reduced walled box (acrylic wall d on z = +/-C, Robin h on the outer
wall faces, m = 4, 119 808 cells) is run three times, each Crank-Nicolson with two
backward-Euler half-steps, q on from t = 0, to t*:

- **R0**, the reference, from theta_0 = 0;
- **R1**, t_w = 0, from theta_0 = g z on every cell (gel and wall);
- **R3**, R1 with the wall's rho_c scaled by ``MASSLESS`` (a massless wall).

The comparand of R1 (`MAT-21` step 33 / P62) is the h_eff drift ``_drift_1d`` plus item
682's two-layer series: it treats the wall as a massless resistance *in the drift*. R3 is
the case where that is exact up to discretisation; R1 - R3 and (R1 - R0) - ``_drift_1d``
are the printed measure of what the comparand leaves out.

**It asserts, and it does not re-implement** (the `ANS-1` rule). The bands, the box
builders, the z-only comparand helpers and the readers are imported by name from
``tests/validation/test_mat21_wall_pre_equilibration_fem.py`` (item 805) and the modules
that gate imports from. That gate's ``runs`` fixture is a pytest fixture and is not called:
its body, without R2, is written out here with every call it makes, plus ``np.real`` on
each integral (the group sources the complex build). No ``test_*`` name is imported. No
gate edit, no ``src/`` edit.

Asserted (bands imported; every gate value collected first, then asserted): (I) interior
rel L2 of R1 against ``_drift_1d`` + the two-layer series <= ``INTERIOR_BAND``; (I-mw) the
same for R3 against its exact comparand; (E) the balance <= ``E_BAND`` on all three runs
with ``source_power`` > 0 and identical on all three; the unique-z precondition <=
``PRE_BAR``; R0's control reproduces item 682's t* record within ``REPRO_BOUND``; (c-wf)
the wall-free series (P51's comparand at H) misses R1 by more than ``INTERIOR_BAND``;
max|Im theta| = 0 on all three runs.
Printed only (records, no band): at P5, ``ADJ_PT`` and P0 the FEM drift part (R1 - R0),
``_drift_1d``, their difference (the wall-capacity term) and R1 - R3, beside item 805's
real-build values.

The wall's k and rho_c are nominal PMMA values, citation owed (W3). One rung, one time step,
this fixture's g. Not claimed: the D = 0.25 case, a wall-capacity comparand, a rate, a
recommended wait.

Run it::

    ./run_examples.sh -e mat:37 -n 2 -t 480

Output: ``examples/materials/paraview_output/
materials_37_gel_box_wall_pre_equilibration_bias_combined.xdmf`` (theta at t* on R0, R1 and
R3, and R1's comparand, on the walled box).
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
from tests.validation.test_mat15_finite_cylinder_closed_form import H, K, Q, T_INF  # noqa: E402
from tests.validation.test_mat15_finite_cylinder_fem import _max_imag, _rel  # noqa: E402
from tests.validation.test_mat21_box_closed_form import A_BOX, B_BOX, C_BOX, T_C  # noqa: E402
from tests.validation.test_mat21_box_fem import _exact, _trans_ex  # noqa: E402
from tests.validation.test_mat21_initial_gradient_closed_form import G_GRAD, GC  # noqa: E402
from tests.validation.test_mat21_initial_gradient_fem import _grad_ex  # noqa: E402
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
    PRE_BAR,
    READ_LABELS,
    READ_PTS,
    REPRO_BOUND,
    _drift_z,
    _read,
    _series_z,
    _ZField,
)
from tests.validation.test_mat21_wall_transient_closed_form import T_STAR  # noqa: E402
from tests.validation.test_mat21_wall_transient_fem import MASSLESS, RUNG, _rho_c_field  # noqa: E402
from tests.validation.test_mat21_wall_pre_equilibration_closed_form import (  # noqa: E402
    _drift_1d,
)

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"  # committed, EX-57
BASENAME = "materials_37_gel_box_wall_pre_equilibration_bias"

#: Item 805's real-build records, `20261104T060527Z_MAT-21-P63-W1b.log:66-68,76-78` (printed
#: beside this run's values; records, not bands).
ITEM_805 = {
    "I1": 8.105518e-03, "Imw": 2.284679e-03, "ref": 2.244928e-03, "wf": 2.1130e-02,
    "drift": {"P5": 1.828983e-01, "ADJ_PT": 1.828894e-01},
    "d1d": {"P5": 1.611619e-01, "ADJ_PT": 1.674396e-01},
    "term": {"P5": 2.174e-02, "ADJ_PT": 1.545e-02},
    "r13": {"P5": -6.931e-02, "ADJ_PT": -5.664e-02, "P0": -1.042e-02},
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
    _say(comm, "EX-180 -- the skipped wait inside the walled box: the wall's heat capacity "
               "in the drift (MAT-21 step 36 / P63)")
    _say(comm, "=" * 78)
    _say(comm, f"[gate] imported: INTERIOR_BAND {INTERIOR_BAND:.0e}, E_BAND {E_BAND:.0e}, "
               f"PRE_BAR {PRE_BAR:.0e}, REPRO_BOUND {REPRO_BOUND:g}, MASSLESS {MASSLESS:g}; "
               f"g {G_GRAD} K/m, g C {GC} K")
    assert K == K_GEL and H == H_GEL, (K, K_GEL, H, H_GEL)

    # ---- mesh: item 682's 1D-reduced walled box ----------------------------------------
    msh, ns, nw = _walled_box(RUNG, (False, False, True))
    chi, kf, inter, ct = _materials(msh)
    rcf = _rho_c_field(msh)
    rc0 = fem.Function(rcf.function_space)
    rc0.x.array[:] = np.where(np.real(rcf.x.array) == RHO_C_W, MASSLESS * RHO_C_W, rcf.x.array)
    ft = _robin_facets(msh, (A_BOX, B_BOX, C_BOX + D_WALL), (2,))
    n = msh.topology.index_map(3).size_global
    write_setup_figure(
        msh,
        ct,
        FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names={GEL: "gel, q on from t = 0, theta_0 = g z", WALL: "acrylic wall (nominal)"},
        title="mat:37 -- gel box with acrylic walls on z = +/-C, Robin on the outer wall "
              "faces, T(0) = T_inf + g z: y = 0 slice",
        slice_normal=(0.0, 1.0, 0.0),
        slice_origin=(0.0, 0.0, 0.0),
        comm=comm,
    )
    _say(comm, f"[mesh] walled box m = {RUNG} cells {tuple(ns)} N_w = {nw}: n = {n}; "
               f"T_C {T_C:.3f} s, t* {T_STAR:.3f} s, dt {DT:.4f} s; ranks {comm.size}")

    # ---- three CN runs (the gate's `runs` body without R2, written out) ----------------
    grad = lambda x: T_INF + G_GRAD * x[2]  # noqa: E731  (every cell, gel and wall)
    kw = dict(k=kf, q=Q * chi, h=H_GEL, t_inf=T_INF, t_init=T_INF, robin_tag=ROBIN, theta=0.5,
              startup_half_steps=2, dt=DT)
    specs = {"R0": dict(rho_c=rcf, n_steps=N_MAIN),
             "R1": dict(rho_c=rcf, n_steps=N_MAIN, t_init_expr=grad),
             "R3": dict(rho_c=rc0, n_steps=N_MAIN, t_init_expr=grad)}
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

    # ---- comparands: z-only, on the unique dof z ---------------------------------------
    tq = time.perf_counter()
    zf = _ZField()
    s_w = zf(msh, _series_z(RHO_C_W), T_STAR)
    s_0 = zf(msh, _series_z(MASSLESS * RHO_C_W), T_STAR)
    d_0 = zf(msh, _drift_z, T_STAR)
    pre = comm.allreduce(zf.worst, op=MPI.MAX)
    _say(comm, f"[series] three z-only comparands on the unique dof z ({zf.n_unique} values "
               f"per rank max): {time.perf_counter() - tq:.1f} s; precondition worst rel "
               f"{pre:.3e} (bar {PRE_BAR:.0e})")
    tq = time.perf_counter()
    wf_d, npts = _exact(msh, _grad_ex(T_STAR))
    wf_s, _ = _exact(msh, _trans_ex(T_STAR))
    _say(comm, f"[series] wall-free control, two 3D interpolations ({npts} P2 points): "
               f"{time.perf_counter() - tq:.1f} s")

    th = {key: r.temperature - T_INF for key, r in res.items()}
    ref = _rel(msh, th["R0"], s_w, inter)
    I1 = _rel(msh, th["R1"], s_w + d_0, inter)
    Imw = _rel(msh, th["R3"], s_0 + d_0, inter)
    wf = _rel(msh, th["R1"], wf_d + wf_s, inter)
    whole1 = _rel(msh, th["R1"], s_w + d_0)
    wholemw = _rel(msh, th["R3"], s_0 + d_0)
    _say(comm, f"[I] R1 interior rel L2 {I1:.6e} (bound {INTERIOR_BAND:.0e}; item 805 "
               f"{ITEM_805['I1']:.6e}); whole {whole1:.6e} (printed)")
    _say(comm, f"[I-mw] R3 interior rel L2 vs its exact comparand {Imw:.6e} (bound "
               f"{INTERIOR_BAND:.0e}; item 805 {ITEM_805['Imw']:.6e}); whole {wholemw:.6e} "
               f"(printed)")
    _say(comm, f"[control] R0 interior {ref:.6e} vs item 682 {ITEM_682_TSTAR:.6e}, |delta| "
               f"{abs(ref - ITEM_682_TSTAR):.3e} (bound {REPRO_BOUND:g}); [c-wf] wall-free "
               f"series miss of R1 {wf:.4e} (asserted > {INTERIOR_BAND:.0e}; item 805 "
               f"{ITEM_805['wf']:.4e})")

    # ---- printed: the wall-capacity term -----------------------------------------------
    _say(comm, "[TABLE] point | FEM drift (R1 - R0) K | _drift_1d K | FEM - _drift_1d K "
               "(wall-capacity term) | R1 - R3 K | item 805: drift, _drift_1d, term, R1 - R3")
    for i, lab in ((I_P5, "P5"), (I_ADJ, "ADJ_PT"), (I_P0, "P0")):
        p = tuple(READ_PTS[i])
        dh = th_pt["R1"][i] - th_pt["R0"][i]
        d1 = _drift_1d(p, T_STAR)
        r13 = th_pt["R1"][i] - th_pt["R3"][i]
        if lab == "P0":
            rec = f"item 805: R1 - R3 {ITEM_805['r13'][lab]:+.3e}"
        else:
            rec = (f"item 805: {ITEM_805['drift'][lab]:+.6e}, {ITEM_805['d1d'][lab]:+.6e}, "
                   f"{ITEM_805['term'][lab]:+.3e}, {ITEM_805['r13'][lab]:+.3e}")
        _say(comm, f"[TABLE] {lab} | {dh:+.6e} | {d1:+.6e} | {dh - d1:+.3e} | {r13:+.3e} | {rec}")
    _say(comm, f"max|Im theta| = {imag:.3e}")

    # ---- output ------------------------------------------------------------------------
    V = res["R0"].temperature.function_space
    outs = {}
    for key in ("R0", "R1", "R3"):
        f = fem.Function(V, name=f"theta_{key}_t_star")
        f.x.array[:] = np.real(res[key].temperature.x.array) - T_INF
        f.x.scatter_forward()
        outs[f.name] = f
    comp = fem.Function(V, name="theta_comparand_R1")
    # P2 comparands -> P1 output by Function interpolation (no extra series evaluation).
    d_p1 = fem.Function(V)
    comp.interpolate(s_w)
    d_p1.interpolate(d_0)
    comp.x.array[:] = comp.x.array + d_p1.x.array
    comp.x.scatter_forward()
    outs[comp.name] = comp
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_combined", msh, ct, outs, comm=comm)
    adopt_host_ownership(OUTPUT_DIR, comm=comm)
    _say(comm, f"\n[output] {path} -- theta(t*) on R0, R1, R3 and R1's comparand on the "
               f"m = {RUNG} walled box")

    # ---- assertions (the gate's), collected first --------------------------------------
    gates = {
        "valid_points": all(all(o) for o in ok_pt.values()),
        "precondition_unique_z": pre <= PRE_BAR,
        "I_R1": I1 <= INTERIOR_BAND,
        "I_mw_R3": Imw <= INTERIOR_BAND,
        "E_source_power_identical": src["R0"] > 0.0 and len(set(src.values())) == 1,
        "control_repro_682": abs(ref - ITEM_682_TSTAR) <= REPRO_BOUND,
        "control_wall_free_misses": wf > INTERIOR_BAND,
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
