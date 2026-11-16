"""Example (`EX-178`): the same watts on a rod as a line and as two tips, on a viewing grid.

A closed form, no solve.  `IMP-8` steps 0-1 (items 759, 760,
`tests/validation/test_imp8_rod_line_source_closed_form.py`) give the temperature rise in an
infinite gel medium from a uniform line source of length L = `ROD_LENGTH` and from two tip
sources of the same total power.  This example evaluates both, per watt of total power, at
t = `T_PRINT` (900 s) on the vertices of a 120 mm cube grid of 41 cells per axis (an odd
count on a symmetric box, so no vertex lies on the line), and compares them.  It is a
viewing grid, not a model: no FEM, no wall, no rod heat capacity, and neither shape is
claimed to be a real rod's.

**It asserts, and it does not re-implement** (the `ANS-1` rule).  Every band, constant and
helper (`_phi`, `_phi_tips`, `_phi_midpoint`, `_point`, `_far2_miss`, `RATE_LO` ...) is
imported by name from the gate module.  Disclosed local pieces: `theta_line` is the gate's
midpoint sum written with the imported `_point` so it broadcasts over the vertices
(`_phi_midpoint` is scalar); the gates' `test_*` bodies for (Rate), (Far2) and (c-one) are
written out because test functions are not imported.

Asserted (collected first): (Rate) the midpoint-sum slope on `RATE_LADDER` inside
[`RATE_LO`, `RATE_HI`] (item 760's, provisional); (Far2) `_far2_miss()` within `FAR2_A`**4
and the gap positive; (Node) the vertex-array values at the vertex nearest (20, 0, 0) mm
equal the gate's scalar forms to `B1_BAND`; (Geom) the nearest vertex to the line is
outside `ROD_RADIUS`; (c-one) `_far2_miss(signs=(1.0,))` misses by more than `FAR2_A`**4;
(c-field) dropping the near tip moves the value at the vertex nearest (0, 0, 15) mm by more
than `B1_BAND`.  max|Im theta| <= 1e-12.

Printed only: item 760's two ratio records re-formed; min and max of the ratio on the grid;
rho c * integral theta dV per watt against 900 J per watt, both shapes.

Run it::

    ./scripts/run_examples.sh -e imp:21 -n 2 -t 180

Output: ``examples/implants/paraview_output/implants_21_rod_line_and_tip_sources_combined.xdmf``.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import ufl
from dolfinx import default_scalar_type, fem
from dolfinx import mesh as dmesh
from dolfinx.fem import assemble_scalar, form
from mpi4py import MPI

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from fem_em_solver.io.paraview_utils import (  # noqa: E402
    adopt_host_ownership,
    write_xdmf_with_tags,
)
from fem_em_solver.post.setup_figure import write_setup_figure  # noqa: E402

from tests.mesh.test_birdcage_implant_rod import ROD_RADIUS  # noqa: E402
from tests.validation.test_imp8_rod_line_source_closed_form import (  # noqa: E402
    B1_BAND,
    FAR2_A,
    K_GEL,
    L,
    MM,
    RATE_HI,
    RATE_LADDER,
    RATE_LO,
    RATE_POINT,
    RHO_C,
    T_PRINT,
    _far2_miss,
    _phi,
    _phi_midpoint,
    _phi_tips,
    _point,
)
from tests.validation.test_mat15_finite_cylinder_fem import _max_imag  # noqa: E402

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"  # committed, EX-57
BASENAME = "implants_21_rod_line_and_tip_sources"
HALF = 0.060  # m, half-width of the viewing grid on each axis
N_CELLS = 41  # per axis; odd on a symmetric box, so no vertex lies on the line
N_FIELD = 64  # midpoint-sum points for theta_line
ITEM_760 = {15.0: 3.179482, 0.0: 0.1714577}  # (print2) records of item 760's window, by z in mm
IM_BAND = 1e-12


def _say(comm, msg):
    if comm.rank == 0:
        print(msg, flush=True)


def _nearest(comm, xyz, target):
    """(owner-reduced) index of the owned vertex nearest `target`, its distance (m)."""
    d = np.linalg.norm(xyz - np.asarray(target, dtype=float) * MM, axis=1)
    i = int(np.argmin(d)) if d.size else -1
    loc = float(d[i]) if i >= 0 else np.inf
    dmin = comm.allreduce(loc, op=MPI.MIN)
    # lowest rank holding the minimum owns the probe
    me = comm.rank if (i >= 0 and loc == dmin) else comm.size
    owner = comm.allreduce(me, op=MPI.MIN)
    return (i if comm.rank == owner else -1), dmin


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
    _say(comm, "EX-178 -- the same watts on a rod as a line and as two tips "
               "(IMP-8 steps 0-1, closed form on a viewing grid)")
    _say(comm, "=" * 78)
    _say(comm, f"[gate] bands imported: B1 {B1_BAND:.0e}, Far2 bound {FAR2_A**4:.0e}, Rate "
               f"[{RATE_LO}, {RATE_HI}] (item 760, provisional); L {L * 1e3:g} mm, rod radius "
               f"{ROD_RADIUS * 1e3:g} mm, t {T_PRINT:g} s, N_FIELD {N_FIELD}")

    # ---- the viewing grid ----------------------------------------------------------------
    msh = dmesh.create_box(comm, [np.array([-HALF] * 3), np.array([HALF] * 3)],
                           [N_CELLS] * 3, cell_type=dmesh.CellType.tetrahedron)
    tdim = msh.topology.dim
    imap_c = msh.topology.index_map(tdim)
    n_all = imap_c.size_local + imap_c.num_ghosts
    cell_tags = dmesh.meshtags(
        msh, tdim, np.arange(n_all, dtype=np.int32), np.full(n_all, 1, dtype=np.int32))
    write_setup_figure(
        msh,
        cell_tags,
        FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names={1: "gel viewing grid, 120 mm cube (a grid, not a model: closed form, no solve)"},
        title="imp:21 -- the same watts as a 30 mm line and as two tips, read on a grid: y = 0 slice",
        slice_normal=(0.0, 1.0, 0.0),
        slice_origin=(0.0, 0.0, 0.0),
        comm=comm,
    )
    ncell = msh.topology.index_map(3).size_global
    V = fem.functionspace(msh, ("Lagrange", 1))
    nd = V.dofmap.index_map.size_local
    xyz = V.tabulate_dof_coordinates()[:nd]
    nglob = V.dofmap.index_map.size_global
    _say(comm, f"[mesh] {N_CELLS}^3 box, half-width {HALF * 1e3:g} mm: {ncell} cells, "
               f"{nglob} P1 dofs; ranks {comm.size}")

    r = np.hypot(xyz[:, 0], xyz[:, 1])
    z = xyz[:, 2]
    r_min = comm.allreduce(float(r.min()) if nd else np.inf, op=MPI.MIN)

    # ---- the two fields, per watt of total power ------------------------------------------
    per_watt = 1.0 / (4.0 * np.pi * K_GEL * L)
    zp = L * ((np.arange(N_FIELD) + 0.5) / N_FIELD - 0.5)
    big_r = np.hypot(r[:, None], z[:, None] - zp[None, :])
    th_line = L / N_FIELD * np.sum(_point(big_r, T_PRINT), axis=1) * per_watt
    th_tips = _phi_tips(r, z, T_PRINT) * per_watt
    ratio = th_tips / th_line

    f_line = fem.Function(V, name="theta_line_per_watt")
    f_tips = fem.Function(V, name="theta_tips_per_watt")
    f_ratio = fem.Function(V, name="ratio_tips_over_line")
    for f, vals in ((f_line, th_line), (f_tips, th_tips), (f_ratio, ratio)):
        f.x.array[:nd] = vals
        f.x.scatter_forward()
    imag = max(_max_imag(f) for f in (f_line, f_tips, f_ratio))

    r_lo = comm.allreduce(float(ratio.min()) if nd else np.inf, op=MPI.MIN)
    r_hi = comm.allreduce(float(ratio.max()) if nd else -np.inf, op=MPI.MAX)

    # ---- energy in the box, rho c * integral theta dV per watt ----------------------------
    dx = ufl.dx(domain=msh)
    energy = {}
    for name, f in (("line", f_line), ("tips", f_tips)):
        loc = np.real(assemble_scalar(form(f * dx)))
        energy[name] = RHO_C * comm.allreduce(float(loc), op=MPI.SUM)

    # ---- (Node): vertex arrays against the gated scalar forms -----------------------------
    i_node, d_node = _nearest(comm, xyz, (20.0, 0.0, 0.0))
    gap_line = gap_tips = 0.0
    node_xyz = None
    if i_node >= 0:
        rn, zn = float(r[i_node]), float(z[i_node])
        ref_line = _phi_midpoint(rn, zn, T_PRINT, N_FIELD) * per_watt
        ref_tips = _phi_tips(rn, zn, T_PRINT) * per_watt
        gap_line = abs(th_line[i_node] / ref_line - 1.0)
        gap_tips = abs(th_tips[i_node] / ref_tips - 1.0)
        node_xyz = xyz[i_node]
    gap_line = comm.allreduce(gap_line, op=MPI.MAX)
    gap_tips = comm.allreduce(gap_tips, op=MPI.MAX)
    node_xyz = comm.bcast(node_xyz, root=comm.allreduce(
        comm.rank if node_xyz is not None else comm.size, op=MPI.MIN))

    # ---- (c-field): dropping a tip at the vertex nearest (0, 0, 15) mm --------------------
    i_tip, d_tip = _nearest(comm, xyz, (0.0, 0.0, 15.0))
    move_near = move_far = 0.0
    tip_xyz = None
    if i_tip >= 0:
        rt, zt = float(r[i_tip]), float(z[i_tip])
        full = th_tips[i_tip]
        near = _phi_tips(rt, zt, T_PRINT, signs=(-1.0,)) * per_watt  # +L/2 tip dropped
        far = _phi_tips(rt, zt, T_PRINT, signs=(1.0,)) * per_watt  # -L/2 tip dropped
        move_near, move_far = abs(near / full - 1.0), abs(far / full - 1.0)
        tip_xyz = xyz[i_tip]
    move_near = comm.allreduce(move_near, op=MPI.MAX)
    move_far = comm.allreduce(move_far, op=MPI.MAX)
    tip_xyz = comm.bcast(tip_xyz, root=comm.allreduce(
        comm.rank if tip_xyz is not None else comm.size, op=MPI.MIN))

    # ---- (Rate), (Far2), (c-one): the gates' bodies, written out --------------------------
    rr, zr = RATE_POINT[0] * MM, RATE_POINT[1] * MM
    phi_ref = _phi(rr, zr, T_PRINT)
    errors = [abs(_phi_midpoint(rr, zr, T_PRINT, n) / phi_ref - 1.0) for n in RATE_LADDER]
    slope = float(-np.polyfit(np.log(RATE_LADDER), np.log(errors), 1)[0])
    far2_miss, far2_gap = _far2_miss()
    c_one = _far2_miss(signs=(1.0,))[0]
    bound = FAR2_A**4

    # ---- printed records ------------------------------------------------------------------
    probe_rec = {}
    for zmm in (15.0, 0.0):
        tips, line = _phi_tips(2.0 * MM, zmm * MM, T_PRINT), _phi(2.0 * MM, zmm * MM, T_PRINT)
        probe_rec[zmm] = tips / line
        _say(comm, f"[print] (2, {zmm:g}) mm t {T_PRINT:g} s: Phi2 {tips:.7g}  Phi {line:.7g}  "
                   f"Phi2/Phi {tips / line:.7g}   (item 760's record {ITEM_760[zmm]})")
    _say(comm, f"[grid] min r on the grid {r_min * 1e3:.4f} mm (rod radius {ROD_RADIUS * 1e3:g} mm)")
    _say(comm, f"[print] ratio theta_tips / theta_line on the grid: min {r_lo:.7g}, max {r_hi:.7g}")
    _say(comm, f"[print] energy rho c int theta dV per watt (900 J/W is the whole run): line "
               f"{energy['line']:.7g} J/W, tips {energy['tips']:.7g} J/W; vs {T_PRINT:g}: "
               f"{energy['line'] / T_PRINT:.5f}, {energy['tips'] / T_PRINT:.5f}")
    _say(comm, f"[node] (20, 0, 0) mm probe: vertex {np.round(node_xyz * 1e3, 4).tolist()} mm "
               f"(distance {d_node * 1e3:.4f} mm); line gap {gap_line:.3e}, tips gap {gap_tips:.3e}")
    _say(comm, f"[c-field] vertex {np.round(tip_xyz * 1e3, 4).tolist()} mm (distance {d_tip * 1e3:.4f} "
               f"mm from (0, 0, 15)): dropping the near tip moves theta_tips by {move_near:.4e}, "
               f"the far tip by {move_far:.4e}")
    _say(comm, f"[Rate] rel {' '.join(f'{e:.3e}' for e in errors)} slope {slope:.4f}")
    _say(comm, f"[Far2] a {FAR2_A} miss {far2_miss:.3e} bound {bound:.3e} gap {far2_gap:.4e}; "
               f"[c-one] miss {c_one:.3e}")
    _say(comm, f"max|Im theta| = {imag:.3e}")

    # ---- gates: collected above, asserted here --------------------------------------------
    _say(comm, f"[gates] (Geom) {r_min * 1e3:.4f} mm > {ROD_RADIUS * 1e3:g}; (Rate) {slope:.4f} in "
               f"[{RATE_LO}, {RATE_HI}]; (Far2) {far2_miss:.3e} <= {bound:.0e}; (Node) "
               f"{gap_line:.3e} / {gap_tips:.3e} <= {B1_BAND:.0e}; (c-one) {c_one:.3e} > "
               f"{bound:.0e}; (c-field) {move_near:.4e} > {B1_BAND:.0e}")
    assert r_min > ROD_RADIUS, f"(Geom) a vertex at r = {r_min:.4e} m is inside the rod"
    assert RATE_LO <= slope <= RATE_HI, f"(Rate) slope {slope:.4f}"
    assert far2_miss <= bound, f"(Far2) {far2_miss:.3e}"
    assert far2_gap > 0.0, f"(Far2) gap {far2_gap:.4e}"
    assert gap_line <= B1_BAND, f"(Node) line {gap_line:.3e}"
    assert gap_tips <= B1_BAND, f"(Node) tips {gap_tips:.3e}"
    assert c_one > bound, f"(c-one) did not fire: {c_one:.3e}"
    assert move_near > B1_BAND, f"(c-field) did not fire: {move_near:.4e}"
    assert imag <= IM_BAND, f"complex build left Im theta = {imag:.3e}"
    assert np.isfinite(r_lo) and np.isfinite(r_hi) and r_lo > 0.0

    # ---- outputs --------------------------------------------------------------------------
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    xdmf = OUTPUT_DIR / f"{BASENAME}_combined"
    write_xdmf_with_tags(
        xdmf, msh, cell_tags,
        {"theta_line_per_watt": f_line, "theta_tips_per_watt": f_tips,
         "ratio_tips_over_line": f_ratio}, comm=comm)
    adopt_host_ownership(OUTPUT_DIR, comm=comm)
    _say(comm, f"\n[output] {xdmf}.xdmf -- theta per watt for both shapes and their ratio, "
               f"P1 on the {ncell}-cell viewing grid")
    _say(comm, f"\n[timing] total {time.perf_counter() - started:.1f} s at -n {comm.size}")
    _say(comm, "\nAll imported gates hold.")


if __name__ == "__main__":
    main()
