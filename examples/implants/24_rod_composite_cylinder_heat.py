"""Example (`EX-184`): the rod's own heat capacity and conductivity -- the composite cylinder,
closed form and FEM.

`imp:15` and `imp:17`-`imp:23` treat the implant as gel or as a hole. Here the rod is a
material: a core 0 < r < a (``ROD_K``, ``ROD_RHO_C``, a **nominal** titanium pair taken from
memory, no source in the repo) in perfect contact with a gel shell a < r < b (``K_GEL``,
``RHO_C``), a uniform source ``Q_SHELL`` (1 W per metre of rod) in the shell only, theta = 0 at
r = b, switched on at t = 0. a = ``A_ROD`` = 1 mm, b = 5 mm.

Closed form: ``_theta`` (`IMP-8` step 9 / item 843, the composite-region series of Carslaw &
Jaeger ch. XIII, cited from memory) with the core's k and rho*c as keywords -- Ti, the gel core
(``k_r=K_GEL, rhoc_r=RHO_C``) and an adiabatic hole (``hole=True``). FEM: item 837's
construction (`IMP-8` step 10), written out in ``main`` because the gate module's fixtures are
pytest fixtures: ``_composite_cylinder`` (imported, never copied), DG0 k / rho*c / q through
``_dg0``, one Crank-Nicolson run (N = 80 steps of t*/80, two backward-Euler half-steps, a
snapshot at step 16 = t_c) and one steady solve, degree 2 on a curved P2 mesh, the wall a
near-Dirichlet Robin at ``H_DIR``.

**It asserts, and it does not re-implement** (the `ANS-1` rule). Every band, helper and fixture
constant is imported by name from ``tests/validation/test_imp8_composite_cylinder_closed_form.py``
and ``test_imp8_composite_cylinder_fem.py`` and the modules they import. No gate edit, no
``src/`` edit.

Asserted (collected first): (I) |theta_FEM / theta_cf - 1| <= ``INTERIOR_BAND`` at (a, 0) and
((a+b)/2, 0) at t* and steady; (I0) the same at (0, 0), t_c; (E) <= ``E_BAND``; (S0) the integral
of q equals ``Q_SHELL`` * V_shell,h within ``B_BAND`` and the integral over the core is 0; the
closed form re-checked on the page (St) and (Lim) at ``B1_BAND`` through the module's helpers;
the control (the k_r = ``K_GEL`` comparand misses theta_FEM(0, t_c) by more than
``INTERIOR_BAND``, factor 1, item 837 read 0.1521); max|Im theta| <= 1e-12.

Printed only: theta(r, t) for the three cores at r in {0, a, (a+b)/2}, t in {t_c, t*, steady}
beside theta_FEM for Ti; the ordering at r = a; z-invariance at (a, +/- L/4); the Robin offset
(predicted 2.8e-5).

Run it::

    ./scripts/run_examples.sh -e imp:24 -n 2 -t 480

Output: ``examples/implants/paraview_output/implants_24_rod_composite_cylinder_heat_combined.xdmf``
(theta_FEM at t* and steady on P2, with the k and rho*c DG0 fields).
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
from fem_em_solver.thermal import solve_steady_conduction, solve_transient_conduction  # noqa: E402

from tests.validation.test_imp4_pec_cylinder_closed_form import B1_BAND  # noqa: E402
from tests.validation.test_imp8_composite_cylinder_closed_form import (  # noqa: E402
    A_ROD,
    B_WALL,
    KAPPA_GEL,
    Q_SHELL,
    R_PTS,
    ROD_K,
    ROD_RHO_C,
    T_C,
    T_ST,
    T_STAR,
    _theta,
    _theta_j0,
    _theta_ss,
)
from tests.validation.test_imp8_composite_cylinder_fem import (  # noqa: E402
    CELL_CAP,
    DEGREE,
    H_DIR,
    LEN,
    MM,
    N_STEPS,
    READ_PTS,
    SNAP,
    T_INF,
    WIT_PTS,
    _composite_cylinder,
)
from tests.validation.test_imp8_rod_line_source_fem import _theta_pts  # noqa: E402
from tests.validation.test_mat10_sar_to_heat_sphere import B_BAND, ROBIN_TAG, _allsum  # noqa: E402
from tests.validation.test_mat12_transient_sar_heating_sphere import E_BAND, RHO_C  # noqa: E402
from tests.validation.test_mat13_transient_perfusion_sphere import INTERIOR_BAND  # noqa: E402
from tests.validation.test_mat15_coil_driven_heating_phantom import K_GEL  # noqa: E402
from tests.validation.test_mat15_finite_cylinder_fem import _max_imag  # noqa: E402
from tests.validation.test_mat18_two_layer_sphere_fem import _dg0  # noqa: E402

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"  # committed, EX-57
BASENAME = "implants_24_rod_composite_cylinder_heat"

IMAG_BAND = 1e-12
#: Item 837's W1 records (`20261105T225620Z_IMP-8-item837-W1.log`), printed beside the readings.
REC_837 = {"cells": 53_504, "I": 6.017e-04, "I0": 1.657e-03, "E": 3.382e-13, "S0": 2.776e-15,
           "ctrl": 0.1521, "robin": 2.7717e-05}
PRED_ROBIN = 2.8e-5  # predicted, printed (rule (e))

CORES = (("Ti", {}), ("gel", dict(k_r=K_GEL, rhoc_r=RHO_C)), ("hole", dict(hole=True)))


def _say(comm, msg):
    if comm.rank == 0:
        print(msg, flush=True)


def _one(msh):
    """1.0 in the build's scalar type: a bare ``fem.Constant(msh, 1.0)`` is float64 and the
    complex build's form constructor rejects it."""
    return fem.Constant(msh, np.asarray(1.0, dtype=default_scalar_type))


def _p1_real(msh, T, name):
    """Real copy of ``T`` in the mesh's own coordinate degree (2): XDMF refuses a
    Lagrange function of another degree on a curved mesh."""
    V1 = fem.functionspace(msh, ("Lagrange", msh.geometry.cmaps[0].degree))
    tmp = fem.Function(V1)
    tmp.interpolate(T)
    out = fem.Function(V1, name=name)
    out.x.array[:] = np.real(tmp.x.array)
    out.x.scatter_forward()
    return out


def _dg0_real(f, name):
    out = fem.Function(f.function_space, name=name)
    out.x.array[:] = np.real(f.x.array)
    out.x.scatter_forward()
    return out


def main() -> None:
    comm = MPI.COMM_WORLD
    started = time.perf_counter()

    if not np.issubdtype(np.dtype(default_scalar_type), np.complexfloating):
        raise RuntimeError(
            "this example expects the complex DolfinX build the runner's `imp:` group "
            "sources (source /usr/local/bin/dolfinx-complex-mode)"
        )

    _say(comm, "=" * 78)
    _say(comm, "EX-184 -- the rod's own heat capacity and conductivity: the composite cylinder "
               "(IMP-8 steps 9-10 / I64, I65)")
    _say(comm, "=" * 78)
    _say(comm, f"[gate] imported: INTERIOR_BAND {INTERIOR_BAND:.0e}, E_BAND {E_BAND:.0e}, "
               f"B_BAND {B_BAND:.0e}, B1_BAND {B1_BAND:.0e}; a {A_ROD / MM:g} mm, b {B_WALL / MM:g} "
               f"mm, L {LEN / MM:g} mm, H_DIR {H_DIR:.0e}, N {N_STEPS}, degree {DEGREE}")
    _say(comm, f"[fixture] t_c {T_C:.7g} s, t* {T_STAR:.7g} s, Q_SHELL {Q_SHELL:.7g} W/m3 (1 W/m), "
               f"kappa_gel {KAPPA_GEL:.7g}; Ti pair ROD_K {ROD_K:g}, ROD_RHO_C {ROD_RHO_C:g} "
               f"(nominal, from memory, no source in the repo)")

    # ---- the closed form on the page (numpy / scipy; helpers imported) -------------------
    scale = float(_theta_ss(A_ROD))
    st = [abs(_theta(r, T_ST) - float(_theta_ss(r))) / scale for r in R_PTS]
    lim = []
    for t in (T_C, T_STAR):
        for r in R_PTS:
            lim.append(abs(_theta(r, t, k_r=K_GEL, rhoc_r=RHO_C) - _theta_j0(r, t)) / scale)
    _say(comm, f"[cf] theta_ss(a) {scale:.7f} K per W/m; (St) max {max(st):.3e}; (Lim) max "
               f"{max(lim):.3e} (B1_BAND {B1_BAND:.0e})")

    # ---- mesh ----------------------------------------------------------------------------
    t0 = time.perf_counter()
    msh, ct, ft = _composite_cylinder(comm)
    n = msh.topology.index_map(msh.topology.dim).size_global
    t_build = time.perf_counter() - t0
    write_setup_figure(
        msh,
        ct,
        FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names={1: "rod core (r < a, own k and rho c)", 2: "gel shell (a < r < b, source)"},
        title="imp:24 -- composite cylinder, r = b wall (Robin tag): y = 0 slice "
              "through the axis",
        slice_normal=(0.0, 1.0, 0.0),
        slice_origin=(0.0, 0.0, 0.0),
        comm=comm,
    )
    _say(comm, f"[mesh] cells {n} (cap {CELL_CAP}; item 837 read {REC_837['cells']}); "
               f"coordinate degree {msh.geometry.cmaps[0].degree}; ranks {comm.size}; "
               f"t_build {t_build:.1f} s")

    # ---- properties and source -----------------------------------------------------------
    k = _dg0(msh, ct, ROD_K, K_GEL)
    rc = _dg0(msh, ct, ROD_RHO_C, RHO_C)
    q = _dg0(msh, ct, 0.0, Q_SHELL)
    nan_loc = int(sum(np.count_nonzero(np.isnan(np.real(g.x.array))) for g in (k, rc, q)))
    nan_left = comm.allreduce(nan_loc, op=MPI.SUM)
    dx = ufl.Measure("dx", domain=msh, subdomain_data=ct)
    v_shell = float(np.real(_allsum(comm, _one(msh) * dx(2))))
    s0 = abs(float(np.real(_allsum(comm, q * dx))) / (Q_SHELL * v_shell) - 1.0)
    s0_core = abs(float(np.real(_allsum(comm, q * dx(1)))))

    # ---- FEM: transient and steady -------------------------------------------------------
    dt = T_STAR / N_STEPS
    t1 = time.perf_counter()
    res = solve_transient_conduction(
        msh, ft, k=k, rho_c=rc, q=q, h=H_DIR, t_inf=T_INF, t_init=T_INF,
        robin_tag=ROBIN_TAG, dt=dt, n_steps=N_STEPS, theta=0.5, startup_half_steps=2,
        degree=DEGREE, snapshot_steps=(SNAP,))
    t_tr = time.perf_counter() - t1
    e_bal = float(np.real(res.max_balance_residual))
    t2 = time.perf_counter()
    st_res = solve_steady_conduction(msh, ft, k=k, q=q, h=H_DIR, t_inf=T_INF,
                                     robin_tag=ROBIN_TAG, degree=DEGREE)
    t_st = time.perf_counter() - t2
    _say(comm, f"[run] t_call transient {t_tr:.1f} s (CN N={N_STEPS}, dt {dt:.6g} s, snapshot "
               f"{SNAP} -> t = {SNAP * dt:.6g} s = t_c {T_C:.6g}), n_factorizations "
               f"{res.n_factorizations}; t_call steady {t_st:.1f} s; NaN left {nan_left}")
    _say(comm, f"[run] (E) {e_bal:.3e} (E_BAND {E_BAND:.0e}; item 837 {REC_837['E']:.3e}); "
               f"(S0) {s0:.3e} (B_BAND {B_BAND:.0e}; item 837 {REC_837['S0']:.3e}), "
               f"core integral {s0_core:.3e}")

    th_c = _theta_pts(res.snapshots[SNAP], READ_PTS)
    th_s = _theta_pts(res.temperature, READ_PTS + WIT_PTS)
    th_ss = _theta_pts(st_res.temperature, READ_PTS)
    fem_ = {"t_c": th_c, "t*": th_s[:3], "ss": th_ss}
    times = {"t_c": T_C, "t*": T_STAR, "ss": np.inf}

    # ---- (I), (I0) -----------------------------------------------------------------------
    i_gated, i0 = [], None
    for key in ("t_c", "t*", "ss"):
        for (r, z), v in zip(READ_PTS, fem_[key]):
            e = float(_theta(r * MM, times[key]))
            rel = abs(v / e - 1.0)
            gated = (key in ("t*", "ss") and r > 0.0) or (key == "t_c" and r == 0.0)
            tag = "(I0)" if key == "t_c" and r == 0.0 else "(I)"
            if gated and tag == "(I)":
                i_gated.append(rel)
            if gated and tag == "(I0)":
                i0 = rel
            _say(comm, f"[I] {tag} ({r:g}, {z:g}) mm {key}: theta_FEM {v:.7e}, composite {e:.7e} "
                       f"K per W/m, rel {rel:.3e} ({'gated' if gated else 'printed only'}, "
                       f"INTERIOR_BAND {INTERIOR_BAND:.0e})")
    _say(comm, f"[I] max over the four (I) reads {max(i_gated):.3e} (item 837 "
               f"{REC_837['I']:.3e}); (I0) {i0:.3e} (item 837 {REC_837['I0']:.3e})")

    # ---- control -------------------------------------------------------------------------
    c0 = float(_theta(0.0, T_C, k_r=K_GEL))
    ctrl = abs(th_c[2] / c0 - 1.0)
    _say(comm, f"[ctrl] k_r = K_GEL comparand at (0, t_c): {c0:.7e}, miss {ctrl:.4e} (asserted > "
               f"INTERIOR_BAND, factor 1; item 837 read {REC_837['ctrl']})")

    # ---- printed record: three cores beside FEM -------------------------------------------
    at_a = {}
    for name, kw in CORES:
        for key in ("t_c", "t*", "ss"):
            vals = []
            for r in R_PTS:
                if kw.get("hole") and r < A_ROD:
                    vals.append("n/a")
                    continue
                v = float(_theta(r, times[key], **kw))
                vals.append(f"{v:.7g}")
                if r == A_ROD:
                    at_a[(name, key)] = v
            _say(comm, f"[rec] {name:4s} {key:3s} theta(0, a, (a+b)/2) = {', '.join(vals)} "
                       f"K per W/m")
    for key in ("t_c", "t*"):
        order = sorted(("Ti", "gel", "hole"), key=lambda nm: -at_a[(nm, key)])
        _say(comm, f"[rec] ordering at r = a, {key}: {' > '.join(order)}")
    for key in ("t_c", "t*", "ss"):
        v3 = fem_[key]  # READ_PTS order: a, (a+b)/2, 0
        _say(comm, f"[rec] Ti theta_FEM {key:3s} at (0, a, (a+b)/2) = "
                   f"{v3[2]:.7g}, {v3[0]:.7g}, {v3[1]:.7g} K per W/m")
    _say(comm, f"[rec] z-invariance at t*: theta_FEM(a, 0) {th_s[0]:.7e}, (a, +L/4) "
               f"{th_s[3]:.7e}, (a, -L/4) {th_s[4]:.7e}")

    ds = ufl.Measure("ds", domain=msh, subdomain_data=ft)
    area = float(np.real(_allsum(comm, _one(msh) * ds(ROBIN_TAG))))
    t_wall = float(np.real(_allsum(comm, st_res.temperature * ds(ROBIN_TAG)))) / area
    _say(comm, f"[rec] Robin offset: mean wall theta / theta_ss(a) = {t_wall / scale:.4e} "
               f"(predicted {PRED_ROBIN:.1e}; item 837 read {REC_837['robin']:.4e}); steady "
               f"convected / source - 1 = "
               f"{st_res.convected_power / st_res.source_power - 1:.3e}")

    # ---- output --------------------------------------------------------------------------
    f_t = _p1_real(msh, res.temperature, "theta_t_star")
    f_s = _p1_real(msh, st_res.temperature, "theta_steady")
    f_k = _dg0_real(k, "k_W_per_mK")
    f_rc = _dg0_real(rc, "rho_c_J_per_m3K")
    imag = max(_max_imag(res.temperature), _max_imag(res.snapshots[SNAP]),
               _max_imag(st_res.temperature))
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_combined", msh, ct,
        {f_t.name: f_t, f_s.name: f_s, f_k.name: f_k, f_rc.name: f_rc}, comm=comm)
    adopt_host_ownership(OUTPUT_DIR, comm=comm)
    _say(comm, f"\n[output] {path} -- theta_FEM at t* = {T_STAR:.7g} s and steady (K per W/m, "
               f"P2) with the k and rho*c DG0 fields")
    _say(comm, f"max|Im theta| = {imag:.3e}")

    # ---- assertions, collected first -----------------------------------------------------
    gates = {
        "cells_le_cap": n <= CELL_CAP,
        "coord_degree_2": msh.geometry.cmaps[0].degree == 2,
        "no_nan": nan_left == 0,
        "E": e_bal <= E_BAND,
        "S0": s0 <= B_BAND,
        "S0_core_zero": s0_core == 0.0,
        "I_four_reads": len(i_gated) == 4 and max(i_gated) <= INTERIOR_BAND,
        "I0_axis_t_c": i0 is not None and i0 <= INTERIOR_BAND,
        "cf_St": max(st) <= B1_BAND,
        "cf_Lim": max(lim) <= B1_BAND,
        "control_k_gel": ctrl > INTERIOR_BAND,
        "max_imag": imag <= IMAG_BAND,
    }
    _say(comm, f"[gates] {gates}")
    _say(comm, f"\n[timing] total {time.perf_counter() - started:.1f} s at -n {comm.size}")
    assert all(gates.values()), gates
    _say(comm, "\nAll imported gates hold.")


if __name__ == "__main__":
    try:
        main()
    except BaseException:  # one rank raising must not leave the other waiting in a collective
        import traceback

        traceback.print_exc()
        sys.stdout.flush()
        sys.stderr.flush()
        MPI.COMM_WORLD.Abort(1)
