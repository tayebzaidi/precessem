"""Example (`EX-141`): the same implant SAR source around an insulating implant.

`imp:11` (`EX-140`) heats the gel around a *metal* core (k = 20 W/m/K) that carries
heat across itself. This one runs `IMP-7` step 1's other boundary case, the
**adiabatic hole**, on one rung (the middle one, 1.0 mm, `LADDER[1]`): the
implant is removed, the gel shell a < r < b is cut out of the two-layer ball as its
own submesh, and the inner surface r = a is left a natural (zero-flux) boundary,
which is how an implant is treated when its own conduction is ignored. The
outer surface is a Robin film. The heat source is the beta = 1 quasi-static
dipole SAR of the PEC sphere (q_far = 2000 W/m^3 far away), solved at P1 on the
submesh and compared with the closed form's adiabatic branch (k_core = 0).

**It asserts, and it does not re-implement** (the `ANS-1` rule). The bands,
the constants and every helper are imported by name from `IMP-7` step 1
(``test_imp7_implant_heating_fem.py``), its closed-form module
(``test_imp7_implant_heating_closed_form.py``) and `MAT-18`'s FEM module
(``test_mat18_two_layer_sphere_fem.py``). The gate's ``_hole_rung`` is *not*
called whole: it ends in ``_errors``, which squares the error as ``d**2`` (not a
norm in the complex build the ``imp:`` group runs in; `EX-140` met the same).
This example therefore rebuilds the rung's pieces in ``_hole_rung``'s own order
(``_ball`` -> ``_shell_submesh`` -> ``_robin_outer`` -> ``_source_dg2`` ->
``_solve`` -> ``_p2(_theta_shell(., _coeffs(0.0)))``), keeps its (G) asserts, and
computes the identical error quantities with ``ufl.inner(d, d)`` (``_errors_c``).

Asserted: (G) the submesh census (sub-cells = parent tag-2 cells, shell r bounds,
Robin area within G_AREA_BAND); (S2) interior rel L2 <= S2_BAND; (S3) balance <=
S3_BAND; (Im) max |Im T_h| = 0; (c-iso) the isotropic source (l = 2 term dropped)
misses the interior band.

Printed only, never gated: pole / equator rise at r = 1.1a beside the closed form
and beside `imp:11`'s metal-core digits (quoted from its log, not re-solved);
the source power against the closed-form shell power; the inner-surface flux
against the source power.

Run it::

    ./run_examples.sh -e imp:12 -n 2 -t 180

Output: ``examples/implants/paraview_output/
implants_12_insulating_implant_heating_in_gel_shell_combined.xdmf``.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import ufl
from dolfinx import default_scalar_type, fem, mesh as dmesh
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
from fem_em_solver.utils.analytical import pec_sphere_implant_radial  # noqa: E402

from tests.validation.test_imp7_implant_heating_closed_form import (  # noqa: E402
    H,
    K_GEL,
    Q_FAR,
    T_INF,
)
from tests.validation.test_imp7_implant_heating_fem import (  # noqa: E402
    A,
    AREA_A,
    AREA_B,
    B,
    BETA,
    CONTROL_RES,
    G_AREA_BAND,
    LADDER,
    R0_READ,
    ROBIN_TAG,
    S2_BAND,
    S3_BAND,
    _ball,
    _coeffs,
    _l2_readout,
    _p2,
    _q_shell,
    _rmin_rmax,
    _robin_outer,
    _shell_submesh,
    _solve,
    _source_dg2,
    _theta_shell,
    _vertex_radii,
    _volume,
)
from tests.validation.test_mat18_two_layer_sphere_fem import (  # noqa: E402
    QDEG,
    _allsum,
    _indicator,
)

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"  # committed, EX-57
BASENAME = "implants_12_insulating_implant_heating_in_gel_shell"
RES = LADDER[1]

#: Gate records (IMP-7 step 1 hole case, real build), printed beside this run.
#: Source: 20261022T202503Z_IMP-7-step1-W2.log (:63-64 rung 1.0 mm, :68 finest rung, :75 control).
RECORD = {
    "cells": 150334,       # :63
    "area": -5.0690e-04,   # :63 robin_area/4 pi b^2 - 1
    "e_in": 8.6093e-04,    # :64
    "bal": 2.571e-14,      # :64
    "src": 6.802147e-02,   # :64 (W)
    "c_iso": 3.6062e-02,   # :75
    "pole": 0.393152,      # :68 finest rung
    "equator": 0.327785,   # :68 finest rung
}

#: imp:11 (EX-140) metal-core digits at r = 1.1a, quoted from 20261023T005150Z_EX-140.log:55.
IMP11 = {"fem": (0.358025, 0.345290), "cf": (0.361610, 0.343661)}


def _say(comm, msg):
    if comm.rank == 0:
        print(msg, flush=True)


def _errors_c(msh, T_h, theta):
    """(whole, interior |x|^2 < 0.81 b^2) rel L2 on theta; `inner`, complex-safe."""
    comm = msh.comm
    dx = ufl.dx(domain=msh, metadata=QDEG)
    d = (T_h - T_INF) - theta
    chi_in = _indicator(msh, lambda m: np.sum(m**2, axis=1) < 0.81 * B**2)
    whole = np.sqrt(max(float(np.real(_allsum(comm, ufl.inner(d, d) * dx))), 0.0)
                    / float(np.real(_allsum(comm, ufl.inner(theta, theta) * dx))))
    inner = np.sqrt(float(np.real(_allsum(comm, chi_in * ufl.inner(d, d) * dx)))
                    / float(np.real(_allsum(comm, chi_in * ufl.inner(theta, theta) * dx))))
    return float(whole), float(inner)


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
    _say(comm, "EX-141 -- the same SAR source around an insulating implant "
               "(IMP-7 step 1 hole case, one rung)")
    _say(comm, "=" * 78)
    _say(comm, f"[gate] bands imported: S2 interior {S2_BAND:.0e}, S3 balance {S3_BAND:.0e}, "
               f"G_AREA_BAND {G_AREA_BAND:.0e}; hole a = {A} m (zero flux), shell to b = {B} m "
               f"(k {K_GEL}), q_far {Q_FAR}, h {H}, T_inf {T_INF}, beta {BETA}; rung {RES} m")

    # ---- the rung, in _hole_rung's order -------------------------------------------------------
    msh, ct, ft, t_build = _ball(RES)
    c2 = ct.find(2)
    n2 = comm.allreduce(int(np.sum(c2 < msh.topology.index_map(3).size_local)), op=MPI.SUM)
    sub = _shell_submesh(msh, ct)
    sft = _robin_outer(sub, A, B)
    nsub = sub.topology.index_map(3).size_global
    n_loc = sub.topology.index_map(3).size_local
    _say(comm, f"[mesh] parent res {RES}: build {t_build:.1f} s; shell submesh {nsub} cells "
               f"(parent tag-2 cells {n2})")

    # Every submesh cell is the gel shell: tag them all 2 (the writer and the figure need tags).
    sub_ct = dmesh.meshtags(sub, 3, np.arange(n_loc, dtype=np.int32),
                            np.full(n_loc, 2, dtype=np.int32))
    write_setup_figure(
        sub,
        sub_ct,
        FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names={2: f"gel shell a = 5 mm to b = 20 mm (k {K_GEL} W/m/K, dipole SAR); "
                         "empty hole at r < a, zero flux"},
        slice_normal=(0.0, 1.0, 0.0),
        slice_origin=(0.0, 0.0, 0.0),
        title="imp:12 -- gel shell around an insulating hole, Robin outside (IMP-7 hole case)",
        window_size=(1400, 630),
        comm=comm,
    )

    # ---- (G) the submesh census ----------------------------------------------------------------
    ds = ufl.Measure("ds", domain=sub, subdomain_data=sft)
    one = fem.Constant(sub, default_scalar_type(1.0))
    area_b = float(np.real(_allsum(comm, one * ds(ROBIN_TAG))))
    area_a = float(np.real(_allsum(comm, one * ds(2))))
    srmin, srmax = _rmin_rmax(comm, _vertex_radii(sub, np.arange(n_loc, dtype=np.int32)))
    g_area = area_b / AREA_B - 1
    _say(comm, f"\n[census] (G) sub_cells {nsub} (record {RECORD['cells']}) = parent tag-2 {n2}; "
               f"shell rmin/a-1 {srmin / A - 1:.3e}, rmax/b-1 {srmax / B - 1:.3e}; "
               f"Robin area/4 pi b^2-1 {g_area:.4e} (record {RECORD['area']:.4e}) "
               f"ASSERTED <= {G_AREA_BAND:.0e}; r=a area/4 pi a^2-1 {area_a / AREA_A - 1:.4e} (printed)")

    # ---- the hole solve ------------------------------------------------------------------------
    q_f = _source_dg2(sub, None, _q_shell)
    sol, t_solve = _solve(sub, sft, K_GEL, q_f)
    c0 = _coeffs(0.0)
    theta = _p2(sub, lambda x: _theta_shell(x, c0))
    nan = comm.allreduce(int(np.count_nonzero(~np.isfinite(theta.x.array))), op=MPI.SUM)
    e_all, e_in = _errors_c(sub, sol.temperature, theta)
    src = float(np.real(sol.source_power))
    conv = float(np.real(sol.convected_power))
    bal = abs(conv - src) / src
    p_cf = float(_coeffs(0.0)["power"])
    _say(comm, f"\n[steady] (S2) interior rel L2 {e_in:.4e} (record {RECORD['e_in']:.4e}); "
               f"whole {e_all:.4e}  ASSERTED interior <= {S2_BAND:.0e}; solve {t_solve:.1f} s")
    _say(comm, f"[steady] (S3) balance {bal:.3e} (record {RECORD['bal']:.3e})  "
               f"ASSERTED <= {S3_BAND:.0e}")
    _say(comm, f"[source] assembled source power {src:.6e} W (gate record {RECORD['src']:.6e} W); "
               f"closed-form shell power {p_cf:.6e} W; ratio {src / p_cf:.6f} "
               f"(shell volume {float(np.real(_volume(sub))):.6e} m^3; printed)")

    # ---- (c-iso) negative control (same rung; CONTROL_RES == RES) -------------------------------
    assert CONTROL_RES == RES, "the control is meant to run on the same rung"
    q_iso = _source_dg2(sub, None, lambda x: _q_shell(x, drop_l2=True))
    sol_iso, _ = _solve(sub, sft, K_GEL, q_iso)
    _, e_iso = _errors_c(sub, sol_iso.temperature, theta)
    _say(comm, f"[control] (c-iso) isotropic source (l = 2 term dropped): interior miss {e_iso:.4e} "
               f"(record {RECORD['c_iso']:.4e}); ASSERTED > {S2_BAND:.0e}; "
               f"{e_iso / S2_BAND:.1f}x the band")

    # ---- printed only --------------------------------------------------------------------------
    r2_cf = float(pec_sphere_implant_radial(R0_READ, c0, 2))
    th_pole = float(_theta_shell(np.array([[0.0], [0.0], [R0_READ]]), c0)[0])
    th_eq = float(_theta_shell(np.array([[R0_READ], [0.0], [0.0]]), c0)[0])
    _say(comm, "")
    _l2_readout(sol.temperature, "hole", (r2_cf, th_pole, th_eq))
    pts = np.array([[0.0, 0.0, R0_READ], [R0_READ, 0.0, 0.0]])
    vals, ok = evaluate_vector_field_parallel(sol.temperature, pts)
    pole_h = eq_h = float("nan")
    if comm.rank == 0:
        v = np.real(np.ravel(vals[:, 0]) if np.ndim(vals) == 2 else np.ravel(vals)) - T_INF
        pole_h, eq_h = float(v[0]), float(v[1])
    pole_h = comm.bcast(pole_h, root=0)
    eq_h = comm.bcast(eq_h, root=0)
    split_hole_cf = th_pole - th_eq
    split_metal_cf = IMP11["cf"][0] - IMP11["cf"][1]
    _say(comm, f"[readout] pole / equator rise at r = 1.1a: solved {pole_h:.6f} / {eq_h:.6f} K; "
               f"closed form {th_pole:.6f} / {th_eq:.6f} K; gate finest-rung records "
               f"{RECORD['pole']:.6f} / {RECORD['equator']:.6f} K. A cancellation readout; no band.")
    _say(comm, f"[compare] imp:11 metal core (20261023T005150Z_EX-140.log:55): solved "
               f"{IMP11['fem'][0]:.6f} / {IMP11['fem'][1]:.6f} K, closed form "
               f"{IMP11['cf'][0]:.6f} / {IMP11['cf'][1]:.6f} K")
    _say(comm, f"[compare] closed-form pole - equator split: hole {split_hole_cf:.6f} K, "
               f"metal core {split_metal_cf:.6f} K, ratio {split_hole_cf / split_metal_cf:.2f}x "
               f"(the metal core carries heat across itself and flattens the l = 2 pattern; "
               f"an insulating implant does not)")

    n = ufl.FacetNormal(sub)
    flux = float(np.real(_allsum(comm, K_GEL * ufl.inner(ufl.grad(sol.temperature), n) * ds(2))))
    _say(comm, f"[flux] inner surface r = a: int k grad T_h . n dS = {flux:.6e} W; "
               f"flux / source power = {flux / src:.3e} (natural boundary; printed, no band)")

    im = max(comm.allreduce(float(np.max(np.abs(np.imag(f.x.array)))) if f.x.array.size else 0.0,
                            op=MPI.MAX) for f in (sol.temperature, sol_iso.temperature))
    _say(comm, f"[complex] max |Im T_h| over both solves = {im:.3e}  ASSERTED == 0")

    # ---- assertions ----------------------------------------------------------------------------
    assert nsub == n2, f"(G) hole cells {nsub} != parent tag-2 {n2}"
    assert srmin >= A * (1 - 1e-9), f"(G) hole vertex at r={srmin!r} < a"
    assert srmax <= B * (1 + 1e-9), f"(G) hole vertex at r={srmax!r} > b"
    assert abs(g_area) <= G_AREA_BAND, f"(G-hole) Robin area {g_area:.3e}"
    assert nan == 0, f"hole comparand holds {nan} NaN dofs"
    assert e_in <= S2_BAND, f"(S2) interior {e_in:.4e}"
    assert bal <= S3_BAND, f"(S3) balance {bal:.3e}"
    assert e_iso > S2_BAND, f"(c-iso) did not fire: {e_iso:.4e}"
    assert im == 0.0, f"max |Im T_h| = {im:.3e}"

    # ---- outputs -------------------------------------------------------------------------------
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    V1 = sol.temperature.function_space
    th = fem.Function(V1, name="theta_h")
    th.x.array[:] = np.real(sol.temperature.x.array) - T_INF
    th.x.scatter_forward()
    th_ex = fem.Function(V1, name="theta_exact")
    th_ex.interpolate(theta)
    th_ex.x.scatter_forward()
    err = fem.Function(V1, name="theta_abs_error")
    err.x.array[:] = np.abs(np.real(th.x.array - th_ex.x.array))
    err.x.scatter_forward()
    q_p1 = fem.Function(V1, name="q_source")
    q_p1.interpolate(q_f)
    q_p1.x.scatter_forward()
    xdmf = OUTPUT_DIR / f"{BASENAME}_combined"
    write_xdmf_with_tags(xdmf, sub, sub_ct, {"theta_h": th, "theta_exact": th_ex,
                                              "theta_abs_error": err, "q_source": q_p1}, comm=comm)
    adopt_host_ownership(OUTPUT_DIR, comm=comm)
    _say(comm, f"\n[output] {xdmf}.xdmf -- P1 fields on the {nsub}-cell shell submesh "
               f"(every cell tagged 2)")
    _say(comm, f"\n[timing] total {time.perf_counter() - started:.1f} s at -n {comm.size}")
    _say(comm, "\nAll imported gates hold.")


if __name__ == "__main__":
    main()
