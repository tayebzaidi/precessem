"""Example (`EX-156`): the implant-heated gel ball on curved (P2-geometry) cells.

`imp:14` solves `IMP-7` step 1/2's metal-core case at P2 on *affine* cells,
where the error floor is the chord of the faceted r = a / r = b spheres
(interior rel L2 8.4803e-04 at `LADDER[1]`, rate ~ 2). `IMP-7` step 3 (I42)
meshed the same fixture with a degree-2 coordinate element
(``geometry_order=2``). This example solves the metal core on that curved mesh
at P2 on the same rung and narrates the cut, then runs the beta = 0 control on
the same curved mesh. It also builds the affine mesh at the same rung for the
volume comparison only -- no solve on it.

**It asserts, and it does not re-implement** (the `ANS-1` rule). Bands, degree
route and helpers are imported by name from the gate modules
(``test_imp7_implant_heating_fem`` as ``m7`` and ``test_imp7_curved_p2``).
Disclosed local pieces: the gate's ``_curved_metal_rung`` / ``m7._errors`` square
the error as ``d**2`` (not a norm in the complex build), so ``_rung`` writes the
rung body out with ``rt = DEG_ROUTE[2]`` and a local ``_errors_c`` using
``ufl.inner(d, d)``; ``_volumes_c`` writes the gate's two volume integrals out
with ``np.real``.

Asserted: (cmap) curved coordinate degree 2, affine 1; (S3c) balance <=
S3_BAND on the curved rung and the control; (S2c) interior rel L2 <= S2_BAND
on curved LADDER[1]; (c-beta0c) the beta = 0 source on curved LADDER[1] misses
S2_BAND.

Printed only: (Geo) volume errors on both meshes, whole / interior error and the
curved / affine ratio beside imp:14's record, the l = 2 readout, cells, build and
solve seconds. No rate is shown (one rung).

Run it::

    ./scripts/run_examples.sh -e imp:16 -n 4 -t 540

Output: ``examples/implants/paraview_output/
implants_16_implant_heating_curved_p2_combined.xdmf``.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import ufl
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

import tests.validation.test_imp7_implant_heating_fem as m7  # noqa: E402
from tests.validation.test_imp7_curved_p2 import (  # noqa: E402
    V_BALL,
    V_CORE,
    _ball_curved,
    _cmap_degree,
    _volumes,  # noqa: F401 -- imported by name; not called (complex-unsafe format), see _volumes_c
)
from tests.validation.test_imp7_implant_heating_closed_form import T_INF  # noqa: E402

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"  # committed, EX-57
BASENAME = "implants_16_implant_heating_curved_p2"
RES = m7.LADDER[1]
assert m7.CONTROL_RES == RES

# imp:14 affine P2 record at this rung, complex build, -n 2
# (20261025T203007Z_EX-154.log:51); widths differ, compare ratios only.
AFFINE_P2_WHOLE = 8.1490e-04
AFFINE_P2_INTERIOR = 8.4803e-04


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


def _volumes_c(msh, ct):
    """`_volumes`' two integrals written out with np.real (complex build)."""
    comm = msh.comm
    dx = ufl.dx(domain=msh, metadata={"quadrature_degree": 4})
    vball = float(np.real(m7._allsum(comm, fem.Constant(msh, default_scalar_type(1.0)) * dx)))
    vcore = float(np.real(m7._allsum(comm, m7._dg0(msh, ct, 1.0, 0.0) * dx)))
    return vcore / V_CORE - 1.0, vball / V_BALL - 1.0


def _rung(q_fn=m7._q_shell):
    """`_curved_metal_rung`'s body (test_imp7_curved_p2 :109-131), complex-safe."""
    rt = m7.DEG_ROUTE[2]
    msh, ct, ft, t_build = _ball_curved(RES)
    k_f = m7._dg0(msh, ct, m7.K_CORE, m7.K_GEL)
    q_f = m7._source_dg2(msh, ct.find(2), q_fn, dg=rt["q_dg"])
    sol, t_solve = m7._solve(msh, ft, k_f, q_f, degree=2)
    theta = m7._p2(msh, m7._theta_metal, p=rt["comparand"])
    e_all, e_in = _errors_c(msh, sol.temperature, theta, rt["qmeta"])
    src = float(np.real(sol.source_power))
    bal = abs(float(np.real(sol.convected_power)) - src) / src
    return msh, ct, sol, theta, t_build, t_solve, e_all, e_in, bal


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
    _say(comm, "EX-156 -- the implant-heated gel ball on curved P2-geometry cells (IMP-7 step 3)")
    _say(comm, "=" * 78)
    _say(comm, f"[gate] bands imported: S2 interior {m7.S2_BAND:.0e}, S3 balance {m7.S3_BAND:.0e}; "
               f"rung {RES} m; P2 route {m7.DEG_ROUTE[2]}")

    # ---- meshes: curved (solved) and affine (volumes only) -------------------------------------
    mc, ctc, ftc, t_build_c = _ball_curved(RES)
    ma, cta, fta, t_build_a = m7._ball(RES)
    dcm, dam = _cmap_degree(mc), _cmap_degree(ma)
    nc = mc.topology.index_map(3).size_global
    na = ma.topology.index_map(3).size_global
    _say(comm, f"[mesh] res {RES}: curved {nc} cells (build {t_build_c:.1f} s), "
               f"affine {na} cells (build {t_build_a:.1f} s); coordinate degree "
               f"curved {dcm}, affine {dam}")
    assert dcm == 2, f"(cmap) curved coordinate degree {dcm}"
    assert dam == 1, f"(cmap) affine coordinate degree {dam}"

    # The setup figure shows the affine twin (same fixture, same cells); the curved
    # mesh differs only by the mid-edge geometry nodes.
    write_setup_figure(
        ma,
        cta,
        FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names={
            1: f"metal core a = 5 mm (k {m7.K_CORE} W/m/K, no source)",
            2: f"gel shell to b = 20 mm (k {m7.K_GEL} W/m/K, dipole SAR)",
        },
        translucent_tags=(2,),
        slice_normal=(0.0, 1.0, 0.0),
        slice_origin=(0.0, 0.0, 0.0),
        title="imp:16 -- metal-core sphere in a gel ball, solved on curved P2-geometry cells",
        window_size=(1400, 630),
        comm=comm,
    )

    # ---- (Geo) printed -------------------------------------------------------------------------
    gc = _volumes_c(mc, ctc)
    ga = _volumes_c(ma, cta)
    _say(comm, f"[geo] Vcore/(4/3 pi a^3)-1: curved {gc[0]:.4e}, affine {ga[0]:.4e}; "
               f"Vball/(4/3 pi b^3)-1: curved {gc[1]:.4e}, affine {ga[1]:.4e} "
               f"(gate record, curved core / ball -2.56e-05 / -1.44e-02; printed, no band)")

    # ---- curved metal core at P2 ---------------------------------------------------------------
    msh, ct, sol, theta, tb, ts, e_all, e_in, bal = _rung()
    _say(comm, f"\n[curved P2] cells {nc}, solve {ts:.1f} s: whole rel L2 {e_all:.4e}, "
               f"interior rel L2 {e_in:.4e}  ASSERTED <= {m7.S2_BAND:.0e}; "
               f"balance {bal:.3e}  ASSERTED <= {m7.S3_BAND:.0e}")
    _say(comm, f"[ratio] curved / affine P2 (imp:14 record {AFFINE_P2_WHOLE:.4e} / "
               f"{AFFINE_P2_INTERIOR:.4e}): whole {e_all / AFFINE_P2_WHOLE:.4f}, "
               f"interior {e_in / AFFINE_P2_INTERIOR:.4f} (printed; gate read 0.072)")
    sr = m7.SERIES_READ["metal"]
    m7._l2_readout(sol.temperature, "metal curved P2", sr)
    th_pole = float(m7._theta_metal(np.array([[0.0], [0.0], [m7.R0_READ]]))[0])
    th_eq = float(m7._theta_metal(np.array([[m7.R0_READ], [0.0], [0.0]]))[0])
    _say(comm, f"[readout] series (R2, pole, equator) = {sr[0]:.4e}, {sr[1]:.6f}, {sr[2]:.6f} "
               f"(closed form pole / equator {th_pole:.6f} / {th_eq:.6f}); gate curved record "
               f"-1.6 % on R2. R2_h is a cancellation readout, no band.")

    # ---- (c-beta0c) control on the same curved mesh --------------------------------------------
    _, _, _, _, _, ts0, e_all0, e_in0, bal0 = _rung(q_fn=lambda x: m7._q_shell(x, beta=0.0))
    _say(comm, f"[control] (c-beta0c) curved P2, beta = 0 source on {RES} m: interior miss "
               f"{e_in0:.4e}; ASSERTED > {m7.S2_BAND:.0e} ({e_in0 / m7.S2_BAND:.1f}x the band; "
               f"gate record 5.9317e-02); whole {e_all0:.4e}; balance {bal0:.3e}; solve {ts0:.1f} s")

    # ---- assertions ----------------------------------------------------------------------------
    assert bal <= m7.S3_BAND, f"(S3c) balance {bal:.3e}"
    assert bal0 <= m7.S3_BAND, f"(S3c) control balance {bal0:.3e}"
    assert e_in <= m7.S2_BAND, f"(S2c) interior {e_in:.4e}"
    assert e_in0 > m7.S2_BAND, f"(c-beta0c) did not fire: {e_in0:.4e}"

    # ---- outputs: P2 Lagrange on the P2-geometry mesh ------------------------------------------
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    V2 = fem.functionspace(msh, ("Lagrange", 2))
    th_h = fem.Function(V2, name="theta_curved_p2")
    th_h.interpolate(sol.temperature)
    th_h.x.array[:] = np.real(th_h.x.array) - T_INF
    th_h.x.scatter_forward()
    th_ex = fem.Function(V2, name="theta_exact")
    th_ex.interpolate(theta)
    th_ex.x.array[:] = np.real(th_ex.x.array)
    th_ex.x.scatter_forward()
    err = fem.Function(V2, name="theta_error")
    err.x.array[:] = np.abs(np.real(th_h.x.array) - np.real(th_ex.x.array))
    err.x.scatter_forward()
    xdmf = OUTPUT_DIR / f"{BASENAME}_combined"
    write_xdmf_with_tags(xdmf, msh, ct,
                         {"theta_curved_p2": th_h, "theta_exact": th_ex, "theta_error": err},
                         comm=comm)
    adopt_host_ownership(OUTPUT_DIR, comm=comm)
    _say(comm, f"\n[output] {xdmf}.xdmf -- P2 fields on the {nc}-cell curved mesh")
    _say(comm, f"\n[timing] total {time.perf_counter() - started:.1f} s at -n {comm.size}")
    _say(comm, "\nAll imported gates hold.")


if __name__ == "__main__":
    main()
