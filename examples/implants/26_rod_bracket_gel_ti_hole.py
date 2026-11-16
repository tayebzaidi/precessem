"""Example (`EX-188`): the rod as gel, as titanium and as a hole -- the bracket on one page.

`imp:23` runs the finite rod's two heat sources in gel; `imp:24` / `imp:25` run the infinite
composite cylinder. Here the same finite rod (r = ``ROD_RADIUS``, length ``L``) is run with three
materials in FEM on one coarse box, per watt (P = 1 W): a 0.20 m gel cube, Robin on all six
faces (``H_GEL``), degree 2, Crank-Nicolson with two backward-Euler half-steps, ``N_MAIN`` steps
to ``T_PRINT`` = 900 s, from theta_0 = 0. Five runs:

1. tips / gel   -- two DG0 balls at z = +/- L/2, P/2 each (`IMP-8` I56's source);
2. tips / Ti    -- the same source, the rod cells (cylinder minus the two balls) at the nominal
                   Ti pair (`IMP-8` I66);
3. tips / hole  -- the same source, the rod cells at ``EPS_HOLE`` x gel k and rho c (I69);
4. rod / gel    -- the DG0 line source on the rod's cells (I55's source);
5. rod / Ti     -- the same source, the body (the cells of q's positive dofs) at the Ti pair (I73).

The box is **not** the gates' box (485 291 cells): it is `imp:23`'s construction, the same
graded gmsh body with the three target sizes multiplied by ``SIZE_SCALE`` -- a **disclosed local
copy** of `imp:23`'s ``_coarse_box`` (examples do not import examples). Every other constant is
imported by name from the gate modules (the `ANS-1` rule); ``box`` is a pytest fixture and is
not imported. The gates' ``_rod_cells`` assembles ``one * w * ufl.dx``, which the complex build
rejects, so the rod and overlap cell selection here is a **local copy of its midpoint
selection without the volume form**. No gate edit, no ``src/`` edit.

Asserted (collected first): the cell count equals ``W0_CELLS`` and is <= ``CELL_LIMIT``; (E) <=
``E_BAND`` on all five runs; (S0): the line's integral of q within ``B_BAND`` of P, each tip ball's
within ``B_BAND`` of P/2, and the selected / rod / overlap / body cell counts all > 0; max|Im
theta| <= 1e-12; no non-finite coefficient dof; the control: each non-gel arm's theta(2, 0) moves
off its same-source gel run's by more than ``CONTROL_FLOOR``.

Printed only: theta(2, 0), theta(2, 15) and their ratio for each run, beside the records on the
gates' mesh; (I) on both gel runs at ``GATE_PTS`` beside ``INTERIOR_BAND``; the Ti / gel / hole
ordering; t_build, t_call and n_factorizations.

Run it::

    ./scripts/run_examples.sh -e imp:26 -n 2 -t 480

Output: ``examples/implants/paraview_output/implants_26_rod_bracket_gel_ti_hole_combined.xdmf``
(the five theta fields at t*, P1, and the DG0 k of the Ti and hole arms).
"""

from __future__ import annotations

import os
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

from fem_em_solver.io.mesh import _model_to_mesh  # noqa: E402
from fem_em_solver.io.paraview_utils import (  # noqa: E402
    adopt_host_ownership,
    write_xdmf_with_tags,
)
from fem_em_solver.post.setup_figure import write_setup_figure  # noqa: E402
from fem_em_solver.thermal import solve_transient_conduction  # noqa: E402

from tests.mesh.test_birdcage_implant_rod import ROD_RADIUS  # noqa: E402
from tests.validation.test_imp8_rod_conducting_body_fem import (  # noqa: E402
    CONTROL_FLOOR,
    EPS_HOLE,
    I56_THETA,
    TI_THETA,
    _coefficients,
)
from tests.validation.test_imp8_composite_cylinder_closed_form import (  # noqa: E402
    ROD_K,
    ROD_RHO_C,
)
from tests.validation.test_imp8_rod_line_source_closed_form import (  # noqa: E402
    L,
    MM,
    T_PRINT,
    _phi,
    _phi_tips,
)
from tests.validation.test_imp8_rod_line_source_fem import (  # noqa: E402
    AXIS_HALF,
    D_IN0,
    D_IN1,
    D_OUT1,
    DEGREE,
    GATE_PTS,
    H_IN,
    H_MID,
    H_OUT,
    HALF_WIDTH,
    N_AXIS_PTS,
    P,
    QP,
    T_INF,
    _cells_global,
    _h_near,
    _rod_source,
    _theta_pts,
)
from tests.validation.test_imp8_rod_self_heated_fem import I55_THETA  # noqa: E402
from tests.validation.test_imp8_rod_tips_fem import RATIO_PTS, SIGNS, _tips_source  # noqa: E402
from tests.validation.test_mat10_sar_to_heat_sphere import B_BAND, ROBIN_TAG, _allsum  # noqa: E402
from tests.validation.test_mat12_transient_sar_heating_sphere import E_BAND, RHO_C  # noqa: E402
from tests.validation.test_mat13_transient_perfusion_sphere import INTERIOR_BAND  # noqa: E402
from tests.validation.test_mat15_coil_driven_heating_phantom import (  # noqa: E402
    H_GEL,
    K_GEL,
    N_MAIN,
)
from tests.validation.test_mat15_finite_cylinder_fem import _max_imag  # noqa: E402

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"  # committed, EX-57
BASENAME = "implants_26_rod_bracket_gel_ti_hole"

#: The three sizes of the gates' box are scaled by this (set from W0 `EX-188-W0`).
SIZE_SCALE = 2.0
#: The W0 build's cell count at ``SIZE_SCALE``, at -n 2: 64 140 (`20261108T210314Z_EX-188-W0.log`;
#: 2.2 gave 48 858, 2.4 gave 38 020).
W0_CELLS = 64_140
CELL_LIMIT = 120_000
W0_SCALES = (2.0, 2.2, 2.4)
IMAG_BAND = 1e-12

#: Records on the gates' 485 291-cell box, printed beside the rows (never asserted). theta in
#: K/W at t* = 900 s, (2, 0) and (2, 15) mm.
REC_TIPS_GEL = ("20261104T043245Z_IMP-8-I56-W1-item796.log", I56_THETA[(2.0, 0.0)],
                I56_THETA[(2.0, 15.0)])
REC_TIPS_TI = ("20261106T013519Z_IMP-8-item846-W1.log:119-121", TI_THETA[(2.0, 0.0)],
               TI_THETA[(2.0, 15.0)])
REC_ROD_GEL = ("20261106T103928Z_IMP-8-item856-W2.log:117-119", 1.772315e01, 9.482226e00)
REC_ROD_TI = ("20261106T103149Z_IMP-8-item856-W1.log:117-119", 1.735615e01, 1.041661e01)


def _say(comm, msg):
    if comm.rank == 0:
        print(msg, flush=True)


def _coarse_box(comm, scale):
    """Disclosed local copy of `imp:23`'s ``_coarse_box`` (itself a copy of the gates'
    ``_graded_box`` with ``H_IN``, ``H_MID``, ``H_OUT`` multiplied by ``scale``)."""
    import gmsh

    rank = 0
    if comm.rank == rank:
        gmsh.initialize()
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("imp26_coarse_box")
        occ = gmsh.model.occ
        hw = HALF_WIDTH
        box = occ.addBox(-hw, -hw, -hw, 2 * hw, 2 * hw, 2 * hw)
        occ.synchronize()
        surfs = [s for _, s in gmsh.model.getBoundary([(3, box)], oriented=False)]
        gmsh.model.addPhysicalGroup(3, [box], tag=1)
        gmsh.model.addPhysicalGroup(2, surfs, tag=ROBIN_TAG)
        pts = [occ.addPoint(0.0, 0.0, float(z)) for z in
               np.linspace(-AXIS_HALF, AXIS_HALF, N_AXIS_PTS)]
        occ.synchronize()
        fm = gmsh.model.mesh.field
        fd = fm.add("Distance")
        fm.setNumbers(fd, "PointsList", pts)
        t_in = fm.add("Threshold")
        fm.setNumber(t_in, "InField", fd)
        fm.setNumber(t_in, "SizeMin", scale * H_IN)
        fm.setNumber(t_in, "SizeMax", scale * H_MID)
        fm.setNumber(t_in, "DistMin", D_IN0)
        fm.setNumber(t_in, "DistMax", D_IN1)
        fm.setNumber(t_in, "StopAtDistMax", 1)
        t_out = fm.add("Threshold")
        fm.setNumber(t_out, "InField", fd)
        fm.setNumber(t_out, "SizeMin", scale * H_MID)
        fm.setNumber(t_out, "SizeMax", scale * H_OUT)
        fm.setNumber(t_out, "DistMin", D_IN1)
        fm.setNumber(t_out, "DistMax", D_OUT1)
        fmin = fm.add("Min")
        fm.setNumbers(fmin, "FieldsList", [t_in, t_out])
        fm.setAsBackgroundMesh(fmin)
        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
        gmsh.model.mesh.generate(3)
    msh, ct, ft = _model_to_mesh(gmsh.model, comm, rank, gdim=3)
    if comm.rank == rank:
        gmsh.finalize()
    return msh, ct, ft


def _rod_cell_sets(msh, radius=ROD_RADIUS, length=L):
    """Local copy of the gates' ``_rod_cells`` selection, without its volume form (the complex
    build rejects ``one * w * ufl.dx``). Owned + ghost cell indices of the rod (cylinder minus
    the two balls), of the overlap (cylinder and ball), and of the balls; allreduced owned
    counts."""
    comm = msh.comm
    tdim = msh.topology.dim
    n_loc = msh.topology.index_map(tdim).size_local
    n_all = n_loc + msh.topology.index_map(tdim).num_ghosts
    mid = dmesh.compute_midpoints(msh, tdim, np.arange(n_all, dtype=np.int32))
    cyl = (np.hypot(mid[:, 0], mid[:, 1]) <= radius) & (np.abs(mid[:, 2]) <= 0.5 * length)
    ball = np.zeros(n_all, dtype=bool)
    for s in SIGNS:
        ball |= np.linalg.norm(mid - np.array([0.0, 0.0, s * 0.5 * length]), axis=1) <= radius
    rod = cyl & ~ball
    over = cyl & ball
    info = {}
    for name, sel in (("rod", rod), ("overlap", over), ("cyl", cyl), ("ball", ball)):
        info[f"n_{name}"] = comm.allreduce(int(np.count_nonzero(sel[:n_loc])), op=MPI.SUM)
    return (np.flatnonzero(rod).astype(np.int32), np.flatnonzero(over).astype(np.int32),
            np.flatnonzero(ball).astype(np.int32), info)


def _body_cells(msh, q):
    """Cell indices (owned + ghost) of q's positive DG0 dofs: I73's body."""
    tdim = msh.topology.dim
    n_loc = msh.topology.index_map(tdim).size_local
    n_all = n_loc + msh.topology.index_map(tdim).num_ghosts
    dofs = np.asarray(q.function_space.dofmap.list).reshape(-1)[:n_all]
    body = np.flatnonzero(np.real(q.x.array[dofs]) > 0.0).astype(np.int32)
    n_body = msh.comm.allreduce(int(np.count_nonzero(body < n_loc)), op=MPI.SUM)
    return body, n_body


def _w0(comm):
    """Build-only window: the cell count at each candidate scale; no solve."""
    _say(comm, "[W0] build-only: no source, no solve")
    for s in W0_SCALES:
        t0 = time.perf_counter()
        msh, _, _ = _coarse_box(comm, s)
        n = _cells_global(msh)
        _say(comm, f"[W0] SIZE_SCALE {s:g}: cells {n}; t_build {time.perf_counter() - t0:.1f} s; "
                   f"ranks {comm.size}; limit {CELL_LIMIT}; within limit {n <= CELL_LIMIT}")


def _p1_theta(msh, T):
    V1 = fem.functionspace(msh, ("Lagrange", 1))
    tmp = fem.Function(V1)
    tmp.interpolate(T)
    out = fem.Function(V1)
    out.x.array[:] = np.real(tmp.x.array) - T_INF
    out.x.scatter_forward()
    return out


def _real_dg0(k):
    """A real-valued copy of a DG0 coefficient for the XDMF."""
    out = fem.Function(k.function_space)
    out.x.array[:] = np.real(k.x.array)
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
    _say(comm, "EX-188 -- the rod as gel, as titanium and as a hole "
               "(IMP-8 I55, I56, I66, I69, I73)")
    _say(comm, "=" * 78)

    if os.environ.get("FEM_EM_IMP26_W0") == "1":
        _w0(comm)
        return
    if SIZE_SCALE is None or W0_CELLS is None:
        raise RuntimeError("SIZE_SCALE / W0_CELLS are set from the W0 window's log")

    _say(comm, f"[gate] imported: E_BAND {E_BAND:.0e}, B_BAND {B_BAND:.0e}, INTERIOR_BAND "
               f"{INTERIOR_BAND:.0e}, CONTROL_FLOOR {CONTROL_FLOOR:.0e}, N_MAIN {N_MAIN}, "
               f"T_PRINT {T_PRINT:g} s, DEGREE {DEGREE}, ROD_RADIUS {ROD_RADIUS / MM:g} mm, "
               f"L {L / MM:g} mm")
    _say(comm, f"[gate] ROD_K {ROD_K:g} W/(m K), ROD_RHO_C {ROD_RHO_C:.6g} J/(m3 K) (nominal); "
               f"K_GEL {K_GEL:g}, RHO_C {RHO_C:.6g}; EPS_HOLE {EPS_HOLE:.0e}")
    _say(comm, f"[local] SIZE_SCALE {SIZE_SCALE:g}; W0_CELLS {W0_CELLS}")

    # ---- mesh and the material cell sets ------------------------------------------------
    t0 = time.perf_counter()
    msh, ct, ft = _coarse_box(comm, SIZE_SCALE)
    n = _cells_global(msh)
    t_build = time.perf_counter() - t0
    rod_idx, over_idx, ball_idx, info = _rod_cell_sets(msh)
    _say(comm, f"[mesh] cells {n} (W0 {W0_CELLS}; limit {CELL_LIMIT}); ranks {comm.size}; "
               f"t_build {t_build:.1f} s")
    _say(comm, f"[cells] cylinder {info['n_cyl']}; rod (cylinder minus balls) {info['n_rod']}; "
               f"overlap (gel properties, carry the source) {info['n_overlap']}; "
               f"ball cells {info['n_ball']}")

    # the figure's tags: 1 gel, 2 rod cells, 3 tip-ball cells (a copy; ``ct`` is not touched)
    n_loc = msh.topology.index_map(msh.topology.dim).size_local
    fig_vals = np.ones(n_loc, dtype=np.int32)
    fig_vals[rod_idx[rod_idx < n_loc]] = 2
    fig_vals[ball_idx[ball_idx < n_loc]] = 3
    ct_fig = dmesh.meshtags(msh, msh.topology.dim, np.arange(n_loc, dtype=np.int32), fig_vals)
    write_setup_figure(
        msh,
        ct_fig,
        FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names={1: "phantom (gel box)", 2: "rod cells (Ti / hole arms)",
                      3: "tip-ball source cells"},
        title="imp:26 -- coarse gel box, rod and tip-ball cells: y = 0 slice through the axis",
        translucent_tags=(1,),
        slice_normal=(0.0, 1.0, 0.0),
        slice_origin=(0.0, 0.0, 0.0),
        comm=comm,
    )

    pts = GATE_PTS + RATIO_PTS  # I56's / I66's order; RATIO_PTS = ((2, 15), (2, 0))
    k15, k0 = pts.index((2.0, 15.0)), pts.index((2.0, 0.0))
    hs = _h_near(msh, pts)
    for (r, z), h in zip(pts, hs):
        _say(comm, f"[mesh] h within 3 mm of ({r:g}, {z:g}) mm: {h / MM:.3f} mm")

    dt = T_PRINT / N_MAIN
    kw = dict(h=H_GEL, t_inf=T_INF, t_init=T_INF, robin_tag=ROBIN_TAG, dt=dt, n_steps=N_MAIN,
              theta=0.5, startup_half_steps=2, degree=DEGREE)

    # ---- sources and their (S0) ---------------------------------------------------------
    q_t, balls = _tips_source(msh)
    s0_t = abs(float(np.real(_allsum(comm, q_t * ufl.dx))) / P - 1.0)
    ball_miss = []
    for s, (n_sel, v_sel, qb) in zip(SIGNS, balls):
        pb = float(np.real(_allsum(comm, qb * ufl.dx)))
        ball_miss.append((n_sel, abs(pb / (0.5 * P) - 1.0)))
        _say(comm, f"[S0] tip ball z={s * 0.5 * L / MM:+g} mm: cells {n_sel}; "
                   f"int q/(P/2) - 1 = {pb / (0.5 * P) - 1.0:.3e}")
    q_l, n_sel_l, v_sel_l = _rod_source(msh)
    s0_l = abs(float(np.real(_allsum(comm, q_l * ufl.dx))) - P) / P
    body_idx, n_body = _body_cells(msh, q_l)
    _say(comm, f"[S0] tips total |int q/P - 1| {s0_t:.3e}; line |int q - P|/P {s0_l:.3e} "
               f"(B_BAND {B_BAND:.0e}); line selected cells {n_sel_l}; body cells {n_body}")

    # ---- the five runs ------------------------------------------------------------------
    k_ti, rc_ti, bad_ti = _coefficients(msh, rod_idx, "ti")
    k_ho, rc_ho, bad_ho = _coefficients(msh, rod_idx, "hole")
    k_bd, rc_bd, bad_bd = _coefficients(msh, body_idx, "ti")
    n_bad = bad_ti + bad_ho + bad_bd
    _say(comm, f"[coef] non-finite coefficient dofs {n_bad}")

    plan = (
        ("tips_gel", "tips / gel", q_t, dict(k=K_GEL, rho_c=RHO_C)),
        ("tips_ti", "tips / Ti", q_t, dict(k=k_ti, rho_c=rc_ti)),
        ("tips_hole", "tips / hole", q_t, dict(k=k_ho, rho_c=rc_ho)),
        ("rod_gel", "rod / gel", q_l, dict(k=K_GEL, rho_c=RHO_C)),
        ("rod_ti", "rod / Ti", q_l, dict(k=k_bd, rho_c=rc_bd)),
    )
    runs = {}
    for key, label, q, coef in plan:
        t1 = time.perf_counter()
        res = solve_transient_conduction(msh, ft, q=q, **coef, **kw)
        call = time.perf_counter() - t1
        e = float(np.real(res.max_balance_residual))
        th = _theta_pts(res.temperature, pts)  # one batch per run
        runs[key] = dict(label=label, res=res, E=e, call=call, th=th)
        _say(comm, f"[{key}] t_call {call:.1f} s; n_factorizations {res.n_factorizations}; "
                   f"(E) {e:.3e} (E_BAND {E_BAND:.0e})")

    # ---- the table ----------------------------------------------------------------------
    recs = {"tips_gel": REC_TIPS_GEL, "tips_ti": REC_TIPS_TI, "rod_gel": REC_ROD_GEL,
            "rod_ti": REC_ROD_TI}
    _say(comm, "\n[table] theta(2, 0), theta(2, 15) mm (K/W at t* = %g s) and the tip-point "
               "ratio theta(2, 15)/theta(2, 0)" % T_PRINT)
    for key, run in runs.items():
        a0, a15 = run["th"][k0], run["th"][k15]
        _say(comm, f"[table] {run['label']:12s} theta(2,0) {a0:.6g}; theta(2,15) {a15:.6g}; "
                   f"ratio {a15 / a0:.6g}")
        if key in recs:
            name, r0, r15 = recs[key]
            _say(comm, f"[table]   record on the 485 291-cell box: theta(2,0) {r0:.6g}; "
                       f"theta(2,15) {r15:.6g}; ratio {r15 / r0:.6g} ({name})")
        else:
            _say(comm, "[table]   record on the 485 291-cell box: item 851, log leak-held")

    # ---- control ------------------------------------------------------------------------
    pairs = (("tips_ti", "tips_gel", 0.555), ("tips_hole", "tips_gel", None),
             ("rod_ti", "rod_gel", -0.0207))
    moves = {}
    for arm, gel, rec in pairs:
        mv = runs[arm]["th"][k0] / runs[gel]["th"][k0] - 1.0
        moves[arm] = mv
        note = (f"; gates' mesh {rec:+.4f}, predicted to repeat in sign" if rec is not None
                else "; gates' mesh: item 851, not reprinted")
        _say(comm, f"[control] {runs[arm]['label']} over {runs[gel]['label']} at (2, 0) mm: "
                   f"theta ratio - 1 = {mv:+.4e} (CONTROL_FLOOR {CONTROL_FLOOR:.0e}){note}")

    # ---- ordering -----------------------------------------------------------------------
    for kk, name in ((k0, "z = 0"), (k15, "z = 15 mm")):
        vals = [("Ti", runs["tips_ti"]["th"][kk]), ("gel", runs["tips_gel"]["th"][kk]),
                ("hole", runs["tips_hole"]["th"][kk])]
        order = " > ".join(nm for nm, _ in sorted(vals, key=lambda x: -x[1]))
        _say(comm, f"[order] tip source, r = 2 mm, {name}: {order} (gates' mesh: Ti > gel > hole "
                   f"at z = 0, hole > gel > Ti at z = 15 mm)")

    # ---- (I), printed -------------------------------------------------------------------
    for key, fn, tag in (("tips_gel", _phi_tips, "run 1 tips / gel vs QP*_phi_tips"),
                         ("rod_gel", _phi, "run 4 rod / gel vs QP*_phi")):
        rels = []
        for (r, z) in GATE_PTS:
            kk = pts.index((r, z))
            ex = QP * fn(r * MM, z * MM, T_PRINT)
            rel = abs(runs[key]["th"][kk] / ex - 1.0)
            rels.append(rel)
            _say(comm, f"[I] {tag} at ({r:g}, {z:g}) mm: theta {runs[key]['th'][kk]:.6e} vs "
                       f"{ex:.6e} rel {rel:.3e}")
        _say(comm, f"[I] {tag}: max {max(rels):.3e} (INTERIOR_BAND {INTERIOR_BAND:.0e}); inside "
                   f"at all four: {all(v <= INTERIOR_BAND for v in rels)} (printed, not asserted)")

    # ---- output -------------------------------------------------------------------------
    fields = {}
    for key, run in runs.items():
        f = _p1_theta(msh, run["res"].temperature)
        f.name = f"theta_{key}_t_star"
        fields[f.name] = f
    for nm, kf in (("k_tips_ti", k_ti), ("k_tips_hole", k_ho), ("k_rod_ti", k_bd)):
        g = _real_dg0(kf)
        g.name = nm
        fields[nm] = g
    imag = max(_max_imag(run["res"].temperature) for run in runs.values())
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_combined", msh, ct, fields, comm=comm)
    adopt_host_ownership(OUTPUT_DIR, comm=comm)
    _say(comm, f"\n[output] {path} -- five theta fields (K per W, t* = {T_PRINT:g} s, P1) and "
               f"the DG0 k of the Ti and hole arms")
    _say(comm, f"max|Im theta| = {imag:.3e}")

    # ---- assertions, collected first ----------------------------------------------------
    gates = {
        "cells_eq_W0": n == W0_CELLS,
        "cells_le_limit": n <= CELL_LIMIT,
        "cell_sets_nonempty": (n_sel_l > 0 and info["n_rod"] > 0 and info["n_overlap"] > 0
                               and n_body > 0),
        "coef_finite": n_bad == 0,
        "S0_line": s0_l <= B_BAND,
        "S0_tips": s0_t <= B_BAND and all(m <= B_BAND and ns > 0 for ns, m in ball_miss),
        "max_imag": imag <= IMAG_BAND,
    }
    for key, run in runs.items():
        gates[f"E_{key}"] = run["E"] <= E_BAND
    for arm, mv in moves.items():
        gates[f"control_{arm}"] = abs(mv) > CONTROL_FLOOR
    _say(comm, f"[gates] {gates}")
    for key, run in runs.items():
        _say(comm, f"[timing] {key}: t_call {run['call']:.1f} s")
    _say(comm, f"\n[timing] total {time.perf_counter() - started:.1f} s at -n {comm.size}; "
               f"t_build {t_build:.1f} s")
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
