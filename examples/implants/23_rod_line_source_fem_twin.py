"""Example (`EX-182`): the rod's heat in FEM -- the line source and its two-tip twin on a coarse box.

`imp:21` and `imp:22` show the finite line source and its two-tip twin in an infinite gel
medium, in closed form. Here both are run in FEM: a 0.20 m gel cube, Robin on all six faces
(``H_GEL``), degree 2, Crank-Nicolson with two backward-Euler half-steps, ``N_MAIN`` steps to
``T_PRINT`` = 900 s, from theta_0 = 0, at P = 1 W.

- **line**: a DG0 source on the cells whose midpoint has r <= ``ROD_RADIUS`` and |z| <= L/2
  (``_rod_source``, `IMP-8` step 2 / item 778);
- **tips**: two DG0 balls of radius ``ROD_RADIUS`` at z = +/- L/2, P/2 each (``_tips_source``,
  step 3 / item 796).

The box is **not** the gates' box. The gates' graded box has 485 291 cells; this one is the
same construction with the three target sizes ``H_IN``, ``H_MID`` and ``H_OUT`` multiplied by
``SIZE_SCALE``, so it is a coarser mesh by design (item 829's ruling). ``_graded_box``'s sizes
are module constants, so ``_coarse_box`` is a **disclosed local copy** of its gmsh body with
that one change; every other constant (``D_IN0``, ``D_IN1``, ``D_OUT1``, ``HALF_WIDTH``,
``AXIS_HALF``, ``N_AXIS_PTS``, the ``StopAtDistMax`` on the inner threshold) is imported or
copied line for line.

**It asserts, and it does not re-implement** (the `ANS-1` rule). Bands, sources, comparands
and readers are imported by name from ``tests/validation/test_imp8_rod_line_source_fem.py``,
``test_imp8_rod_tips_fem.py`` and the modules they import. ``box`` is a pytest fixture and is
not imported; no ``test_*`` name is imported. No gate edit, no ``src/`` edit.

Asserted (collected first): (E) the balance residual <= ``E_BAND`` on both runs; (S0) the line's
integral of q within ``B_BAND`` of P and each tip ball's within ``B_BAND`` of P/2, with every
selected-cell count > 0; the cell count equals ``W0_CELLS`` (the build-only window's figure)
and is <= ``CELL_LIMIT``; max|Im theta| <= 1e-12. (I) at ``GATE_PTS`` for both runs, and the
controls (c-2L) and (c-one), are asserted only where the flagged probe run read them (see the
``ASSERT_*`` flags); otherwise printed beside ``INTERIOR_BAND``.

Printed only: theta_tips / theta_line at (2, 15) mm, FEM over FEM, beside the closed forms'
3.179482; h near each point; t_build, t_call and n_factorizations for each run.

Run it::

    ./scripts/run_examples.sh -e imp:23 -n 2 -t 480

Output: ``examples/implants/paraview_output/implants_23_rod_line_source_fem_twin_combined.xdmf``
(theta_line, theta_tips and their ratio at t*, P1).
"""

from __future__ import annotations

import os
import sys
import time
from pathlib import Path

import numpy as np
import ufl
from dolfinx import default_scalar_type, fem, la
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
from tests.validation.test_imp8_rod_line_source_closed_form import (  # noqa: E402
    L,
    MM,
    T_PRINT,
    _phi,
    _phi_tips,
)
from tests.validation.test_imp8_rod_line_source_fem import (  # noqa: E402
    AXIS_HALF,
    C2L_PT,
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
    _theta_pts,
)
from tests.validation.test_imp8_rod_tips_fem import SIGNS  # noqa: E402
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
BASENAME = "implants_23_rod_line_source_fem_twin"

#: The three sizes of the gates' box are scaled by this (W0 `EX-182-W0` chose it).
SIZE_SCALE = 1.8
#: The W0 build's cell count at ``SIZE_SCALE``: 86 726 at -n 2
#: (`20261105T193505Z_EX-182-W0.log`; 1.6 gave 122 448 > limit, 2.0 gave 64 140).
W0_CELLS = 86_726
CELL_LIMIT = 120_000
W0_SCALES = (1.6, 1.8, 2.0)
#: Set from the flagged window `20261105T194455Z_EX-182-figure.log`: (I) at all four
#: ``GATE_PTS`` read inside ``INTERIOR_BAND`` on both runs (line max 8.730e-03, tips max
#: 6.810e-03), (c-2L) 8.1455e-02 and (c-one) min 1.0092e+00 read above it (rule (e)).
ASSERT_I_LINE = True
ASSERT_I_TIPS = True
ASSERT_C2L = True
ASSERT_CONE = True

RATIO_PT = (2.0, 15.0)  # (r, z) mm
CLOSED_RATIO = 3.179482  # imp:22 / item 811, the closed forms' tips / line at (2, 15) mm
I56_RATIO = 3.186418  # item 796's FEM tips over the line's closed form, the brief's figure
IMAG_BAND = 1e-12


def _say(comm, msg):
    if comm.rank == 0:
        print(msg, flush=True)


def _coarse_box(comm, scale):
    """Disclosed local copy of ``test_imp8_rod_line_source_fem._graded_box`` with
    ``H_IN``, ``H_MID``, ``H_OUT`` multiplied by ``scale``. Everything else is imported."""
    import gmsh

    rank = 0
    if comm.rank == rank:
        gmsh.initialize()
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("imp23_coarse_box")
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


def _cell_volumes(msh):
    """Per-cell DG0 volumes (the gates' ``bvol`` block), with the test function conjugated:
    the gates' ``one * w * ufl.dx`` is rejected by the complex build (``ArityMismatch``,
    `20261105T193522Z_EX-182-figure.log`), and ``ufl.conj`` is the identity in real mode."""
    tdim = msh.topology.dim
    n_loc = msh.topology.index_map(tdim).size_local
    n_all = n_loc + msh.topology.index_map(tdim).num_ghosts
    mid = dmesh.compute_midpoints(msh, tdim, np.arange(n_all, dtype=np.int32))
    W = fem.functionspace(msh, ("DG", 0))
    w = ufl.TestFunction(W)
    one = fem.Constant(msh, np.asarray(1.0, dtype=fem.Function(W).x.array.dtype))
    bvol = fem.assemble_vector(fem.form(one * ufl.conj(w) * ufl.dx))
    bvol.scatter_reverse(la.InsertMode.add)
    bvol.scatter_forward()
    dofs = np.asarray(W.dofmap.list).reshape(-1)[:n_all]
    return W, n_loc, mid, dofs, np.real(bvol.array[dofs])


def _rod_source(msh, power=P, length=L, radius=ROD_RADIUS):
    """Disclosed local copy of the gate's ``_rod_source`` (same selection and
    normalisation) with the conjugated volume form of ``_cell_volumes``."""
    comm = msh.comm
    W, n_loc, mid, dofs, vol = _cell_volumes(msh)
    sel = (np.hypot(mid[:, 0], mid[:, 1]) <= radius) & (np.abs(mid[:, 2]) <= 0.5 * length)
    own = sel[:n_loc]
    v_sel = comm.allreduce(float(np.sum(vol[:n_loc][own])), op=MPI.SUM)
    n_sel = comm.allreduce(int(np.count_nonzero(own)), op=MPI.SUM)
    q = fem.Function(W)
    q.x.array[:] = 0.0
    if v_sel > 0:
        q.x.array[dofs[sel]] = power / v_sel
    q.x.scatter_forward()
    return q, n_sel, v_sel


def _tips_source(msh, power=P, radius=ROD_RADIUS, length=L):
    """Disclosed local copy of the gate's ``_tips_source`` (two balls at z = +/- length/2,
    power/2 each) with the conjugated volume form of ``_cell_volumes``."""
    comm = msh.comm
    W, n_loc, mid, dofs, vol = _cell_volumes(msh)
    q = fem.Function(W)
    q.x.array[:] = 0.0
    balls = []
    for s in SIGNS:
        c = np.array([0.0, 0.0, s * 0.5 * length])
        sel = np.linalg.norm(mid - c, axis=1) <= radius
        own = sel[:n_loc]
        v_sel = comm.allreduce(float(np.sum(vol[:n_loc][own])), op=MPI.SUM)
        n_sel = comm.allreduce(int(np.count_nonzero(own)), op=MPI.SUM)
        qb = fem.Function(W)
        qb.x.array[:] = 0.0
        if v_sel > 0:
            qb.x.array[dofs[sel]] = 0.5 * power / v_sel
        qb.x.scatter_forward()
        q.x.array[:] += qb.x.array
        balls.append((n_sel, v_sel, qb))
    q.x.scatter_forward()
    return q, balls


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


def main() -> None:
    comm = MPI.COMM_WORLD
    started = time.perf_counter()

    if not np.issubdtype(np.dtype(default_scalar_type), np.complexfloating):
        raise RuntimeError(
            "this example expects the complex DolfinX build the runner's `imp:` group "
            "sources (source /usr/local/bin/dolfinx-complex-mode)"
        )

    _say(comm, "=" * 78)
    _say(comm, "EX-182 -- the rod's heat in FEM: the line source and its two-tip twin "
               "(IMP-8 steps 2-3 / I55, I56)")
    _say(comm, "=" * 78)

    if os.environ.get("FEM_EM_IMP23_W0") == "1":
        _w0(comm)
        return
    if SIZE_SCALE is None or W0_CELLS is None:
        raise RuntimeError("SIZE_SCALE / W0_CELLS are set from the W0 window's log")

    _say(comm, f"[gate] imported: E_BAND {E_BAND:.0e}, B_BAND {B_BAND:.0e}, INTERIOR_BAND "
               f"{INTERIOR_BAND:.0e}, N_MAIN {N_MAIN}, T_PRINT {T_PRINT:g} s, DEGREE {DEGREE}, "
               f"ROD_RADIUS {ROD_RADIUS / MM:g} mm, L {L / MM:g} mm")
    _say(comm, f"[local] SIZE_SCALE {SIZE_SCALE:g} on H_IN {H_IN / MM:g}, H_MID {H_MID / MM:g}, "
               f"H_OUT {H_OUT / MM:g} mm; W0_CELLS {W0_CELLS}")

    # ---- mesh ---------------------------------------------------------------------------
    t0 = time.perf_counter()
    msh, ct, ft = _coarse_box(comm, SIZE_SCALE)
    n = _cells_global(msh)
    t_build = time.perf_counter() - t0
    write_setup_figure(
        msh,
        ct,
        FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names={1: "gel box (the rod's line and tip balls are sources, not meshed)"},
        title="imp:23 -- coarse graded gel box, Robin on all faces: y = 0 slice "
              "through the rod axis",
        slice_normal=(0.0, 1.0, 0.0),
        slice_origin=(0.0, 0.0, 0.0),
        comm=comm,
    )
    pts = GATE_PTS + ((2.0, 0.0), RATIO_PT)
    hs = _h_near(msh, pts)
    _say(comm, f"[mesh] cells {n} (W0 {W0_CELLS}; limit {CELL_LIMIT}); ranks {comm.size}; "
               f"t_build {t_build:.1f} s")
    for (r, z), h in zip(pts, hs):
        _say(comm, f"[mesh] h within 3 mm of ({r:g}, {z:g}) mm: {h / MM:.3f} mm")

    dt = T_PRINT / N_MAIN
    kw = dict(k=K_GEL, rho_c=RHO_C, h=H_GEL, t_inf=T_INF, t_init=T_INF, robin_tag=ROBIN_TAG,
              dt=dt, n_steps=N_MAIN, theta=0.5, startup_half_steps=2, degree=DEGREE)

    # ---- run 1: the line -----------------------------------------------------------------
    q_l, n_sel_l, v_sel_l = _rod_source(msh)
    s0_l = abs(float(np.real(_allsum(comm, q_l * ufl.dx))) - P) / P
    t1 = time.perf_counter()
    res_l = solve_transient_conduction(msh, ft, q=q_l, **kw)
    call_l = time.perf_counter() - t1
    e_l = float(np.real(res_l.max_balance_residual))
    _say(comm, f"[line] t_call {call_l:.1f} s; n_factorizations {res_l.n_factorizations}; "
               f"selected cells {n_sel_l}; (E) {e_l:.3e} (E_BAND {E_BAND:.0e}); (S0) "
               f"{s0_l:.3e} (B_BAND {B_BAND:.0e})")
    th_l = _theta_pts(res_l.temperature, pts)

    # ---- run 2: the two tips -------------------------------------------------------------
    q_t, balls = _tips_source(msh)
    s0_t = abs(float(np.real(_allsum(comm, q_t * ufl.dx))) / P - 1.0)
    ball_miss = []
    for s, (n_sel, v_sel, qb) in zip(SIGNS, balls):
        pb = float(np.real(_allsum(comm, qb * ufl.dx)))
        ball_miss.append((n_sel, abs(pb / (0.5 * P) - 1.0)))
        _say(comm, f"[tips] ball z={s * 0.5 * L / MM:+g} mm: cells {n_sel}; int q/(P/2) - 1 = "
                   f"{pb / (0.5 * P) - 1.0:.3e}")
    t1 = time.perf_counter()
    res_t = solve_transient_conduction(msh, ft, q=q_t, **kw)
    call_t = time.perf_counter() - t1
    e_t = float(np.real(res_t.max_balance_residual))
    _say(comm, f"[tips] t_call {call_t:.1f} s; n_factorizations {res_t.n_factorizations}; "
               f"(E) {e_t:.3e} (E_BAND {E_BAND:.0e}); (S0) {s0_t:.3e} (B_BAND {B_BAND:.0e})")
    th_t = _theta_pts(res_t.temperature, pts)

    # ---- (I), controls -------------------------------------------------------------------
    i_l, i_t, c2l, cone = [], [], None, []
    for k, (r, z) in enumerate(pts):
        ex_l = QP * _phi(r * MM, z * MM, T_PRINT)
        ex_t = QP * _phi_tips(r * MM, z * MM, T_PRINT)
        rl, rt = abs(th_l[k] / ex_l - 1.0), abs(th_t[k] / ex_t - 1.0)
        gated = (r, z) in GATE_PTS
        if gated:
            i_l.append(rl)
            i_t.append(rt)
            ex1 = QP * _phi_tips(r * MM, z * MM, T_PRINT, signs=(-1.0,))
            cone.append(abs(th_t[k] / ex1 - 1.0))
        if (r, z) == C2L_PT:
            ex2 = P / (2 * L) / (4.0 * np.pi * K_GEL) * _phi(r * MM, z * MM, T_PRINT, length=2 * L)
            c2l = abs(th_l[k] / ex2 - 1.0)
        _say(comm, f"[I] ({r:g}, {z:g}) mm: line theta {th_l[k]:.6e} vs {ex_l:.6e} rel {rl:.3e}; "
                   f"tips theta {th_t[k]:.6e} vs {ex_t:.6e} rel {rt:.3e} "
                   f"({'gated' if gated else 'printed only'}; INTERIOR_BAND {INTERIOR_BAND:.0e})")
    _say(comm, f"[I] max over GATE_PTS: line {max(i_l):.3e}, tips {max(i_t):.3e} "
               f"(INTERIOR_BAND {INTERIOR_BAND:.0e}); line inside band at all four: "
               f"{all(v <= INTERIOR_BAND for v in i_l)}; tips inside band at all four: "
               f"{all(v <= INTERIOR_BAND for v in i_t)}")
    _say(comm, f"[c-2L] L -> 2L comparand at {C2L_PT} mm misses the line by {c2l:.4e} (asserted "
               f"> INTERIOR_BAND iff ASSERT_C2L = {ASSERT_C2L}; measured on the gates' mesh "
               f"7.517e-02)")
    _say(comm, f"[c-one] one-ball comparand miss, min over GATE_PTS {min(cone):.4e} (asserted "
               f"> INTERIOR_BAND iff ASSERT_CONE = {ASSERT_CONE}; measured on the gates' mesh "
               f"min 1.0045)")
    k_r = pts.index(RATIO_PT)
    _say(comm, f"[record] Phi_2/Phi at {RATIO_PT} mm, FEM over FEM: {th_t[k_r] / th_l[k_r]:.6f} "
               f"(closed forms {CLOSED_RATIO:.6f}; item 796's FEM tips over the line's closed "
               f"form {I56_RATIO:.6f})")

    # ---- output --------------------------------------------------------------------------
    f_l = _p1_theta(msh, res_l.temperature)
    f_t = _p1_theta(msh, res_t.temperature)
    f_l.name, f_t.name = "theta_line_t_star", "theta_tips_t_star"
    f_r = fem.Function(f_l.function_space, name="ratio_tips_over_line")
    a_l, a_t = f_l.x.array, f_t.x.array
    safe = np.abs(a_l) > 1e-12
    f_r.x.array[:] = np.where(safe, a_t / np.where(safe, a_l, 1.0), 0.0)
    f_r.x.scatter_forward()
    imag = max(_max_imag(res_l.temperature), _max_imag(res_t.temperature))
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_combined", msh, ct,
        {f_l.name: f_l, f_t.name: f_t, f_r.name: f_r}, comm=comm)
    adopt_host_ownership(OUTPUT_DIR, comm=comm)
    _say(comm, f"\n[output] {path} -- theta_line, theta_tips (K per W, t* = {T_PRINT:g} s) and "
               f"their ratio on P1")
    _say(comm, f"max|Im theta| = {imag:.3e}")

    # ---- assertions, collected first -----------------------------------------------------
    gates = {
        "cells_eq_W0": n == W0_CELLS,
        "cells_le_limit": n <= CELL_LIMIT,
        "E_line": e_l <= E_BAND,
        "E_tips": e_t <= E_BAND,
        "S0_line": s0_l <= B_BAND and n_sel_l > 0,
        "S0_tips": s0_t <= B_BAND and all(m <= B_BAND and ns > 0 for ns, m in ball_miss),
        "max_imag": imag <= IMAG_BAND,
    }
    if ASSERT_I_LINE:
        gates["I_line"] = max(i_l) <= INTERIOR_BAND
    if ASSERT_I_TIPS:
        gates["I_tips"] = max(i_t) <= INTERIOR_BAND
    if ASSERT_C2L:
        gates["control_c2L"] = c2l > INTERIOR_BAND
    if ASSERT_CONE:
        gates["control_c_one"] = min(cone) > INTERIOR_BAND
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
