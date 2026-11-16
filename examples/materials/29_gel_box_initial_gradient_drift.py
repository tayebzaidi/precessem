"""Example (`EX-161`): the F2182 gel box relaxing from a pre-equilibration gradient --
the drift an F2182 measurement subtracts, against the triple Robin series.

`mat:24` reads the gel box in steady state under a uniform source; `mat:25`/`mat:27`
drive it in time. Every one of them starts from a *uniform* temperature. No earlier
example starts a transient from a **spatial field**. Here the source is off (q = 0),
the box (A, B, C) = (0.325, 0.21, 0.045) m is cooled by Robin on all six faces, and
the initial state is theta_0 = g z with g = 10 K/m (a 1 K gradient over 90 mm, the
kind a phantom carries if it was not left to equilibrate). The CN solution at
t* = 0.2 t_C is compared to the closed form
``utils.analytical.robin_box_initial_gradient_temperature`` (triple Robin series),
and to the same series with Dirichlet-z roots (the control).

**It asserts, and it does not re-implement** (the `ANS-1` rule). The solve wrapper
call, the series wrapper ``_grad_ex``, the control roots, the point-reflection bar and
the probe points are imported by name from
``tests/validation/test_mat21_initial_gradient_fem.py`` (`g`); the rung, interpolation
and interior-weight helpers from ``test_mat21_box_fem``; the bands from the modules
that gate imports them from. The gate's ``gradient_run`` fixture is a pytest fixture
and is not called: its body is written out here with every call it makes, plus
``np.real`` on each integral (the group sources the complex build). No gate edit, no
``src/`` edit.

Asserted (bands imported): (I) interior rel L2 at t* <= ``INTERIOR_BAND``; (E) the
source-free balance <= ``E_BAND``; (Sym') |theta(P1) + theta(P2)|, |theta(P3) +
theta(P4)| and |theta(P0)| each <= ``SYM_BAR`` * g C (point reflection through the
centre, not a z-mirror); (c-Dz) the Dirichlet-z closed form misses the same FEM field
by more than ``INTERIOR_BAND``.
Printed only: the P5 (top face) drift beside the series and beside P37's *steady,
sourced* value (not a transient rise); the four mid-plane values as the O(h^2) record.

Run it::

    ./run_examples.sh -e mat:29 -n 2 -t 300

Output: ``examples/materials/paraview_output/
materials_29_gel_box_initial_gradient_drift_combined.xdmf`` (theta(0), theta(t*), the
series at t* and the error on step 1's finest box).
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import ufl
from dolfinx import default_scalar_type, fem
from dolfinx import mesh as dmesh
from mpi4py import MPI

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from fem_em_solver.io.paraview_utils import (  # noqa: E402
    adopt_host_ownership,
    write_xdmf_with_tags,
)
from fem_em_solver.post.evaluation import evaluate_vector_field_parallel  # noqa: E402
from fem_em_solver.post.setup_figure import write_setup_figure  # noqa: E402
from fem_em_solver.thermal import solve_transient_conduction  # noqa: E402

from tests.validation.test_mat10_sar_to_heat_sphere import (  # noqa: E402
    ROBIN_TAG,
    _allsum,
    _robin_tags,
)
from tests.validation.test_mat12_transient_sar_heating_sphere import E_BAND, RHO_C  # noqa: E402
from tests.validation.test_mat13_transient_perfusion_sphere import INTERIOR_BAND  # noqa: E402
from tests.validation.test_mat15_coil_driven_heating_phantom import (  # noqa: E402
    H_GEL,
    K_GEL,
    N_MAIN,
    T_INF,
)
from tests.validation.test_mat15_finite_cylinder_fem import _max_imag, _rel  # noqa: E402
from tests.validation.test_mat21_box_closed_form import A_BOX, B_BOX, C_BOX, T_C  # noqa: E402
from tests.validation.test_mat21_box_fem import (  # noqa: E402
    BASE,
    RUNGS,
    T_STAR,
    _exact,
    _interior_chi,
)
from tests.validation.test_mat21_box_probe_points_fem import (  # noqa: E402
    PROBE_LABELS,
    PROBE_PTS,
)
from tests.validation.test_mat21_initial_gradient_closed_form import G_GRAD, GC  # noqa: E402
from tests.validation.test_mat21_initial_gradient_fem import (  # noqa: E402
    DZ_ROOTS,
    P37_STEADY_P5,
    SYM_BAR,
    _grad_ex,
)

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"  # committed, EX-57
BASENAME = "materials_29_gel_box_initial_gradient_drift"

#: Item 676's records, real build, `-n 2`:
#: `20261028T003056Z_MAT-21-step19-W1-fem-r2.log:59-69`.  Printed beside this run, a
#: difference reported and never asserted.
RECORD = {
    "interior": 7.543166e-3,
    "balance": 4.106e-15,
    "dz": 2.9493,
    "p5": 0.146994,
    "p5_series": 0.146779,
    "sym12": 7.895e-14,
    "sym34": 6.316e-14,
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
    _say(comm, "EX-161 -- gel box relaxing from theta_0 = g z (MAT-21 step 19 / P47b')")
    _say(comm, "=" * 78)
    _say(comm, f"[gate] imported: INTERIOR_BAND {INTERIOR_BAND:.0e}, E_BAND {E_BAND:.0e}, "
               f"SYM_BAR {SYM_BAR:.0e}; (A, B, C) = {(A_BOX, B_BOX, C_BOX)}, k = {K_GEL}, "
               f"h = {H_GEL}, g = {G_GRAD} K/m, g C = {GC} K")

    # ---- mesh: step 1's finest rung ----------------------------------------------------
    m = RUNGS[-1]
    ns = [b * m for b in BASE]
    msh = dmesh.create_box(comm, [np.array([-A_BOX, -B_BOX, -C_BOX]),
                                  np.array([A_BOX, B_BOX, C_BOX])],
                           ns, cell_type=dmesh.CellType.tetrahedron)
    ft = _robin_tags(msh)
    n = msh.topology.index_map(3).size_global
    chi = _interior_chi(msh)
    tdim = msh.topology.dim
    imap_c = msh.topology.index_map(tdim)
    n_all = imap_c.size_local + imap_c.num_ghosts
    cell_tags = dmesh.meshtags(
        msh, tdim, np.arange(n_all, dtype=np.int32), np.full(n_all, 1, dtype=np.int32))
    write_setup_figure(
        msh,
        cell_tags,
        FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names={1: "gel (650 x 420 x 90 mm), q = 0, theta_0 = g z"},
        title="mat:29 -- gel box, Robin on all six faces, no source, T(0) = T_inf + g z: y = 0 slice",
        slice_normal=(0.0, 1.0, 0.0),
        slice_origin=(0.0, 0.0, 0.0),
        comm=comm,
    )
    _say(comm, f"[mesh] rung m = {m} {tuple(ns)}: cells {n}; t_C = {T_C:.3f} s, "
               f"t* = {T_STAR:.3f} s, CN N = {N_MAIN} + 2 BE half-steps; ranks {comm.size}")

    # ---- CN from g z, q = 0 ------------------------------------------------------------
    ts = time.perf_counter()
    res = solve_transient_conduction(
        msh, ft, k=K_GEL, rho_c=RHO_C, q=0.0, h=H_GEL, t_inf=T_INF, t_init=T_INF,
        robin_tag=ROBIN_TAG, dt=T_STAR / N_MAIN, n_steps=N_MAIN, theta=0.5,
        startup_half_steps=2, degree=1, snapshot_steps=(0,),
        t_init_expr=lambda x: T_INF + G_GRAD * x[2])
    t_solve = time.perf_counter() - ts
    imag = _max_imag(res.temperature)
    src = float(np.real(res.source_power))
    e_bal = float(np.real(res.max_balance_residual))
    nb = len(res.balance_residuals)

    th0 = res.snapshots[0] - T_INF
    xs = ufl.SpatialCoordinate(msh)
    e_ic = _rel(msh, th0, G_GRAD * xs[2])
    int0 = float(np.real(_allsum(comm, th0 * ufl.dx(domain=msh))))
    one = fem.Constant(msh, default_scalar_type(1.0))
    v_h = float(np.real(_allsum(comm, one * ufl.dx(domain=msh))))
    _say(comm, f"[IC] rel L2 |theta_0 - g z| = {e_ic:.3e} (P1 exact, printed); "
               f"int theta_0 dV/(g C V) = {int0 / (GC * v_h):.3e} (printed, symmetry); "
               f"CN solve {t_solve:.1f} s")

    th = res.temperature - T_INF
    ex, npts = _exact(msh, _grad_ex(T_STAR))
    interior = _rel(msh, th, ex, chi)
    whole = _rel(msh, th, ex)
    _say(comm, f"[I] CN at t*: interior rel L2 = {interior:.6e} (record {RECORD['interior']:.6e}, "
               f"bound {INTERIOR_BAND:.0e}); whole {whole:.6e} ({npts} P2 points)")
    _say(comm, f"[E] source power int Q = {src!r} W; balance max {e_bal:.3e} over {nb} (half-)steps "
               f"(record {RECORD['balance']:.3e}, bound {E_BAND:.0e}); first "
               f"{res.balance_residuals[0]:.3e}, last {res.balance_residuals[-1]:.3e}")

    ex_dz, _ = _exact(msh, _grad_ex(T_STAR, roots_z=DZ_ROOTS))
    dz = _rel(msh, th, ex_dz, chi)
    _say(comm, f"[c-Dz] Dirichlet-z closed form at t*: interior miss = {dz:.4e} "
               f"(record {RECORD['dz']:.4e}; asserted > {INTERIOR_BAND:.0e})")

    # ---- probes: point values through the parallel evaluator ---------------------------
    pts = np.array(PROBE_PTS, dtype=float)
    vals, ok = evaluate_vector_field_parallel(res.temperature, pts)
    ok = [bool(o) for o in ok]
    fem_p = np.array([float(np.real(np.ravel(vals[i])[0])) - T_INF if ok[i] else np.nan
                      for i in range(len(PROBE_PTS))])
    ser_p = np.real(_grad_ex(T_STAR)(pts[:, 0], pts[:, 1], pts[:, 2]))
    for i, lab in enumerate(PROBE_LABELS):
        _say(comm, f"[probe] {lab} {PROBE_PTS[i]}: theta_h(t*) = {fem_p[i]:+.6e} K, series "
                   f"{ser_p[i]:+.6e} K, |theta_h|/(g C) = {abs(fem_p[i]) / GC:.3e}")
    p5, s5 = fem_p[5], ser_p[5]
    _say(comm, f"[P5] drift theta(0, 0, C, t*) = {p5:.6f} K (record {RECORD['p5']}); series "
               f"{s5:.6f} K (record {RECORD['p5_series']}); rel {p5 / s5 - 1.0:+.3e}; beside it "
               f"P37 = {P37_STEADY_P5:.6f} K (steady, sourced, not a transient rise); "
               f"ratio drift / P37 = {p5 / P37_STEADY_P5:.6f}")
    s12, s34, s0 = abs(fem_p[1] + fem_p[2]), abs(fem_p[3] + fem_p[4]), abs(fem_p[0])
    for i in range(1, 5):
        _say(comm, f"[O(h^2) record] {PROBE_LABELS[i]} theta_h(t*) = {fem_p[i]:+.6e} K, "
                   f"|theta_h|/(g C) = {abs(fem_p[i]) / GC:.3e} (point-symmetric, not z-mirrored mesh)")
    _say(comm, f"[Sym'] |P1 + P2|/(g C) = {s12 / GC:.3e} (record {RECORD['sym12']:.3e}), "
               f"|P3 + P4|/(g C) = {s34 / GC:.3e} (record {RECORD['sym34']:.3e}), "
               f"|P0|/(g C) = {s0 / GC:.3e} (bar {SYM_BAR:.0e})")
    _say(comm, f"max|Im theta| = {imag:.3e}")

    # ---- output ------------------------------------------------------------------------
    V = res.temperature.function_space
    th_t = fem.Function(V, name="theta_t_star")
    th_t.x.array[:] = np.real(res.temperature.x.array) - T_INF
    th_t.x.scatter_forward()
    th_i = fem.Function(V, name="theta_0")
    th_i.x.array[:] = np.real(res.snapshots[0].x.array) - T_INF
    th_i.x.scatter_forward()
    ser = fem.Function(V, name="theta_series_t_star")
    ser.interpolate(lambda x: np.real(_grad_ex(T_STAR)(x[0], x[1], x[2])))
    err = fem.Function(V, name="theta_error")
    err.x.array[:] = th_t.x.array - ser.x.array
    err.x.scatter_forward()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_combined", msh, cell_tags,
        {"theta_0": th_i, "theta_t_star": th_t, "theta_series_t_star": ser, "theta_error": err},
        comm=comm,
    )
    adopt_host_ownership(OUTPUT_DIR, comm=comm)
    _say(comm, f"\n[output] {path} -- theta(0), theta(t*), series at t*, error on the m = {m} box")

    # ---- assertions (the gate's) -------------------------------------------------------
    assert src == 0.0, f"source power {src}"
    assert imag <= 1e-12, f"complex build left Im theta = {imag:.3e}"
    assert interior <= INTERIOR_BAND, f"(I) {interior:.4e} > {INTERIOR_BAND}"
    assert e_bal <= E_BAND, f"(E) {e_bal:.3e} > {E_BAND:.0e}"
    assert all(ok), ok
    assert s12 <= SYM_BAR * GC, ("P1 + P2", fem_p[1], fem_p[2])
    assert s34 <= SYM_BAR * GC, ("P3 + P4", fem_p[3], fem_p[4])
    assert s0 <= SYM_BAR * GC, ("P0", fem_p[0])
    assert dz > INTERIOR_BAND, f"(c-Dz) {dz:.4e} inside band"

    _say(comm, f"\n[timing] total {time.perf_counter() - started:.1f} s at -n {comm.size}")
    _say(comm, "\nAll imported gates hold.")


if __name__ == "__main__":
    main()
