"""Example (`EX-158`): the implant-heated gel shell switched on, continuous and duty-cycled.

Every earlier `examples/implants/*.py` is steady. `IMP-7` step 5 (I46') and
step 6 (I47) gated the *time-dependent* version on the curved adiabatic-hole
fixture: the Crank-Nicolson transient of the Robin-cooled gel shell heated by
the beta = 1 dipole SAR, against the eigen-series of the hollow sphere. This
example runs both drives on curved `LADDER[0]` at P2:

1. the CN transient from theta = 0 under the *continuous* dipole SAR to t*;
2. the same under the 50 % duty (period 0.5 t*), each switch on a step boundary;
3. each against its series (`_series_fn` / `_duty_fn`, adiabatic hole).

It runs no BE, slope or 2*DT_1 runs and no `LADDER[1]`.

**It asserts, and it does not re-implement** (the `ANS-1` rule). Every name below
is imported from the gate modules: from ``test_imp7_transient_fem`` `_fixture`,
`_series_fn`, `_transient`, `LADDER`, `RT`, `INTERP_CHUNK`; from
``test_imp7_duty_cycle_fem`` `PERIOD`, `DUTY`, `S0_BAND`, `_chi`, `_duty_fn`,
`_census`; `INTERIOR_BAND` and `E_BAND` from the modules those import them from.
Disclosed local piece: ``m7._errors`` squares the error as ``d**2`` (not a norm in
the complex build), so ``_errors_c`` uses ``ufl.inner(d, d)``.

Asserted: (T) interior rel L2 at t* <= INTERIOR_BAND for both runs; (E) CN max
balance <= E_BAND for both; (S0) duty census sum(chi dt) = on-time within
S0_BAND; (c-cont) the continuous-source series misses the duty FEM by more than
INTERIOR_BAND.

Printed only: theta(t) at the pole and the equator of r = 1.1a for both runs beside
the series; the deposited-energy ratio duty / continuous; cells, build and solve
seconds. No rate in h or dt is shown.

Run it::

    ./scripts/run_examples.sh -e imp:17 -n 4 -t 480

Output: ``examples/implants/paraview_output/
implants_17_implant_heating_transient_duty_combined.xdmf``.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import ufl
from dolfinx import default_scalar_type, fem, geometry
from mpi4py import MPI

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from fem_em_solver.io.paraview_utils import (  # noqa: E402
    adopt_host_ownership,
    write_xdmf_with_tags,
)
from fem_em_solver.post.setup_figure import write_setup_figure  # noqa: E402

import tests.validation.test_imp7_implant_heating_fem as m7  # noqa: E402
from tests.validation.test_imp7_implant_heating_closed_form import T_INF  # noqa: E402
from tests.validation.test_imp7_transient_closed_form import T_STAR  # noqa: E402
from tests.validation.test_imp7_transient_fem import (  # noqa: E402
    INTERP_CHUNK,  # noqa: F401 -- used inside the imported series helpers
    LADDER,
    RT,
    _fixture,
    _series_fn,
    _transient,
)
from tests.validation.test_imp7_duty_cycle_fem import (  # noqa: E402
    DUTY,
    PERIOD,
    S0_BAND,
    _census,
    _chi,
    _duty_fn,
)
from tests.validation.test_mat12_transient_sar_heating_sphere import E_BAND, N_MAIN  # noqa: E402
from tests.validation.test_mat13_transient_perfusion_sphere import INTERIOR_BAND  # noqa: E402

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"  # committed, EX-57
BASENAME = "implants_17_implant_heating_transient_duty"
RES = LADDER[0]
DT = T_STAR / N_MAIN
SNAP = (N_MAIN // 4, N_MAIN // 2, 3 * N_MAIN // 4, N_MAIN)  # t*/4 ... t*


def _say(comm, msg):
    if comm.rank == 0:
        print(msg, flush=True)


def _errors_c(msh, T_h, theta, qmeta):
    """(whole, interior |x|^2 < 0.81 b^2) rel L2 on theta; `inner`, complex-safe."""
    comm = msh.comm
    dx = ufl.dx(domain=msh, metadata=qmeta)
    d = (T_h - T_INF) - theta
    chi_in = m7._indicator(msh, lambda m: np.sum(m**2, axis=1) < 0.81 * m7.B**2)
    whole = np.sqrt(max(float(np.real(m7._allsum(comm, ufl.inner(d, d) * dx))), 0.0)
                    / float(np.real(m7._allsum(comm, ufl.inner(theta, theta) * dx))))
    inner_ = np.sqrt(max(float(np.real(m7._allsum(comm, chi_in * ufl.inner(d, d) * dx))), 0.0)
                     / float(np.real(m7._allsum(comm, chi_in * ufl.inner(theta, theta) * dx))))
    return float(whole), float(inner_)


def _point_value(msh, f, pts):
    """Scalar P2 values at points, rank-safe: owners report, others -inf, allreduce MAX."""
    tree = geometry.bb_tree(msh, msh.topology.dim)
    cand = geometry.compute_collisions_points(tree, pts)
    coll = geometry.compute_colliding_cells(msh, cand, pts)
    out = np.full(len(pts), -np.inf)
    for i in range(len(pts)):
        cells = coll.links(i)
        if len(cells) > 0:
            out[i] = float(np.real(f.eval(pts[i:i + 1], cells[:1])[0]))
    red = np.empty_like(out)
    msh.comm.Allreduce(out, red, op=MPI.MAX)
    return red


def main() -> None:
    comm = MPI.COMM_WORLD
    started = time.perf_counter()

    if not np.issubdtype(np.dtype(default_scalar_type), np.complexfloating):
        raise RuntimeError(
            "this example needs the complex DolfinX build "
            "(source /usr/local/bin/dolfinx-complex-mode); the runner's `imp:` "
            "group sources it automatically"
        )

    _say(comm, "=" * 78)
    _say(comm, "EX-158 -- the implant-heated gel shell, switched on and duty-cycled (IMP-7 steps 5-6)")
    _say(comm, "=" * 78)
    _say(comm, f"[gate] bands imported: INTERIOR {INTERIOR_BAND:.0e}, E {E_BAND:.0e}, "
               f"S0 {S0_BAND:.0e}; rung {RES} m; route {RT}; t* = {T_STAR:.4f} s, "
               f"CN dt = t*/{N_MAIN} = {DT:.4f} s, duty {DUTY} of period {PERIOD:.4f} s")

    # Setup figure from the affine twin of the same rung (the curved mesh only adds
    # mid-edge geometry nodes; the solve meshes the shell alone, hole at r = a).
    ma, cta, _fta, _tb = m7._ball(RES)
    write_setup_figure(
        ma,
        cta,
        FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names={
            1: f"adiabatic hole a = 5 mm (cut out of the solve; k {m7.K_CORE} W/m/K in the steady twin)",
            2: f"gel shell to b = 20 mm (k {m7.K_GEL} W/m/K, dipole SAR, Robin-cooled)",
        },
        translucent_tags=(2,),
        slice_normal=(0.0, 1.0, 0.0),
        slice_origin=(0.0, 0.0, 0.0),
        title="imp:17 -- gel shell around an adiabatic hole, CN transient continuous and duty-cycled",
        window_size=(1400, 630),
        comm=comm,
    )

    sub, sft, q_f, t_build = _fixture(RES)
    ncell = sub.topology.index_map(3).size_global
    _say(comm, f"[mesh] curved shell sub-mesh {ncell} cells, build {t_build:.1f} s")

    # ---- the two CN runs ------------------------------------------------------------------
    cont, t_cont = _transient(sub, sft, q_f, dt=DT, n_steps=N_MAIN, theta=0.5,
                              startup_half_steps=2, snapshot_steps=SNAP)
    duty, t_duty = _transient(sub, sft, q_f, dt=DT, n_steps=N_MAIN, theta=0.5,
                              startup_half_steps=2, snapshot_steps=SNAP,
                              source_schedule=_chi(PERIOD, DUTY))
    s0, tot, n_on, n_off = _census(duty, 2)
    _say(comm, f"[solve] continuous {t_cont:.1f} s, duty {t_duty:.1f} s")
    _say(comm, f"[(S0)] census sum chi*dt = {tot:.12f} s vs on-time {2 * DUTY * PERIOD:.12f} s, "
               f"rel {s0:.3e}  ASSERTED <= {S0_BAND:.0e}; (half-)steps on={n_on} off={n_off} "
               f"(record 0)")

    # ---- series and errors at t* ---------------------------------------------------------
    th_cont = m7._p2(sub, _series_fn(T_STAR, "neumann"), p=RT["comparand"])
    th_duty = m7._p2(sub, _duty_fn(T_STAR), p=RT["comparand"])
    for nm, th in (("continuous", th_cont), ("duty", th_duty)):
        nan = comm.allreduce(int(np.count_nonzero(~np.isfinite(th.x.array))), op=MPI.SUM)
        assert nan == 0, f"{nm} comparand holds {nan} NaN dofs"
    e_all_c, e_in_c = _errors_c(sub, cont.temperature, th_cont, RT["qmeta"])
    e_all_d, e_in_d = _errors_c(sub, duty.temperature, th_duty, RT["qmeta"])
    _, e_ctrl = _errors_c(sub, duty.temperature, th_cont, RT["qmeta"])
    _say(comm, f"[(T)] continuous: interior rel L2 {e_in_c:.4e} (whole {e_all_c:.4e}); "
               f"ASSERTED <= {INTERIOR_BAND:.0e}  (record 1.8004e-4, W1r)")
    _say(comm, f"[(T)] duty:       interior rel L2 {e_in_d:.4e} (whole {e_all_d:.4e}); "
               f"ASSERTED <= {INTERIOR_BAND:.0e}  (record 1.0444e-4, I47 W1)")
    _say(comm, f"[(E)] CN max balance: continuous {cont.max_balance_residual:.3e}, "
               f"duty {duty.max_balance_residual:.3e}; ASSERTED <= {E_BAND:.0e} "
               f"(records 5.807e-12 / 5.809e-12)")
    _say(comm, f"[(c-cont)] continuous-source series vs the duty FEM: interior miss "
               f"{e_ctrl:.4e}; ASSERTED > {INTERIOR_BAND:.0e} "
               f"({e_ctrl / INTERIOR_BAND:.0f}x the band; record 0.53873)")

    # ---- printed: theta(t) at the pole and equator of r = 1.1a, FEM beside series ---------
    r0 = m7.R0_READ
    pts = np.array([[0.0, 0.0, r0], [r0, 0.0, 0.0]])
    _say(comm, f"[theta(t)] K rise at pole (0,0,1.1a) / equator (1.1a,0,0); step -> t = step*dt")
    for step in SNAP:
        t = step * DT
        cs = _series_fn(t, "neumann")(pts.T)
        ds = _duty_fn(t)(pts.T)
        fc = _point_value(sub, cont.snapshots[step], pts) - T_INF
        fd = _point_value(sub, duty.snapshots[step], pts) - T_INF
        _say(comm, f"  t = {t:8.3f} s  cont FEM {fc[0]:.6f}/{fc[1]:.6f} series {cs[0]:.6f}/{cs[1]:.6f}"
                   f"   duty FEM {fd[0]:.6f}/{fd[1]:.6f} series {ds[0]:.6f}/{ds[1]:.6f}")
    ratio = tot / T_STAR
    _say(comm, f"[energy] deposited duty / continuous = sum chi*dt / t* = {ratio:.6f} (printed; ~0.5)")

    # ---- assertions ----------------------------------------------------------------------
    assert e_in_c <= INTERIOR_BAND, f"(T) continuous {e_in_c:.4e}"
    assert e_in_d <= INTERIOR_BAND, f"(T) duty {e_in_d:.4e}"
    assert cont.max_balance_residual <= E_BAND, f"(E) continuous {cont.max_balance_residual:.3e}"
    assert duty.max_balance_residual <= E_BAND, f"(E) duty {duty.max_balance_residual:.3e}"
    assert s0 <= S0_BAND, f"(S0) {s0:.3e}"
    assert e_ctrl > INTERIOR_BAND, f"(c-cont) did not fire: {e_ctrl:.4e}"

    # ---- outputs: P2 Lagrange on the P2-geometry mesh ------------------------------------
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    V2 = fem.functionspace(sub, ("Lagrange", 2))

    def _field(name, src, shift=0.0):
        f = fem.Function(V2, name=name)
        f.interpolate(src)
        f.x.array[:] = np.real(f.x.array) - shift
        f.x.scatter_forward()
        return f

    fc = _field("theta_cont_fem", cont.temperature, T_INF)
    fd = _field("theta_duty_fem", duty.temperature, T_INF)
    sc = _field("theta_cont_series", th_cont)
    sd = _field("theta_duty_series", th_duty)
    ec = fem.Function(V2, name="error_cont")
    ec.x.array[:] = np.abs(fc.x.array - sc.x.array)
    ec.x.scatter_forward()
    ed = fem.Function(V2, name="error_duty")
    ed.x.array[:] = np.abs(fd.x.array - sd.x.array)
    ed.x.scatter_forward()
    xdmf = OUTPUT_DIR / f"{BASENAME}_combined"
    write_xdmf_with_tags(xdmf, sub, None,
                         {"theta_cont_fem": fc, "theta_duty_fem": fd,
                          "theta_cont_series": sc, "theta_duty_series": sd,
                          "error_cont": ec, "error_duty": ed},
                         comm=comm)
    adopt_host_ownership(OUTPUT_DIR, comm=comm)
    _say(comm, f"\n[output] {xdmf}.xdmf -- P2 fields on the {ncell}-cell curved shell mesh")
    _say(comm, f"\n[timing] total {time.perf_counter() - started:.1f} s at -n {comm.size}")
    _say(comm, "\nAll imported gates hold.")


if __name__ == "__main__":
    main()
