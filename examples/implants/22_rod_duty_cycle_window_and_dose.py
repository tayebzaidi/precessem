"""Example (`EX-181`): a rod under a duty cycle, as a line and as two tips -- the six-minute window and the dose.

A closed form, no solve.  `IMP-8` step 4 (item 787,
`tests/validation/test_imp8_rod_line_source_duty_closed_form.py`) and step 7 (item 811,
`tests/validation/test_imp8_rod_tips_duty_closed_form.py`) give the temperature rise in an
infinite gel medium from a pulsed source of PERIOD 450 s, for DUTIES (1, 0.5, 0.25), read to
T_STAR = 900 s, spread over the rod as a uniform line or as two tip sources of the same total
power.  This example reads both histories at the two gated points, (r, z) = (2, 15) and
(10, 0) mm, runs the gates' checks on them, prints the six-minute maximum and CEM43 for the
two shapes at (2, 15) mm, and evaluates both shapes at T_STAR under D = 0.5 on the vertices
of `imp:21`'s viewing grid.  Infinite gel medium in closed form: no wall, no rod heat
capacity, no FEM; neither shape is claimed to be a real rod's.

**It asserts, and it does not re-implement** (the `ANS-1` rule).  Every band, constant and
helper is imported by name from the gate modules.  Disclosed local pieces: the gates'
`test_*` bodies (D1), (Duh), (Q-d), (Win-d), (c-cont) are written out because test functions
are not importable; `_line_kernel` is the line's midpoint sum (`N_FIELD` points, written with
the imported `_point`) so it broadcasts over the grid, as `imp:21`'s `theta_line` does, and the
grid's Duhamel sum is the tips module's imported `_duhamel` with a kernel argument, which is
kernel-agnostic and broadcasts.  The (Node) reference at one vertex is the scalar
`_phi_midpoint` through the same `_duhamel`.

Asserted (collected first): (D1) the D = 1 sum equals the step response within `D1_BAND` of
theta_ss; (Duh) step sum against convolution within `DUH_BAND` of theta_ss; (Q-d)
`_q_check_g` within `DOSE_BAR` and (Win-d) `_win_check_g` within `AVG_BAR`, the line at 1 W
(item 787) and the tips at `P_DOSE_W` (item 811's ruling), max T < 43 C on every case; the
reproduction control (theta(t*) against the logs' printed digits, relative <= 1e-9);
(c-cont) the continuous form misses the D = 0.5 history by more than `INTERIOR_BAND`;
(Node), (Geom) as `imp:21`; max|Im theta| <= 1e-12.

Printed only: theta(t*), six-minute maximum and CEM43 for both shapes at (2, 15) mm, per
duty, beside the two items' printed records; theta tips / line; grid min / max of the ratio.
The two doses are at different powers (1 W and 0.5 W), so no dose ratio is printed.

Run it::

    ./scripts/run_examples.sh -e imp:22 -n 2 -t 180

Output: ``examples/implants/paraview_output/implants_22_rod_duty_cycle_window_and_dose_combined.xdmf``.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
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
from fem_em_solver.post.setup_figure import write_setup_figure  # noqa: E402

from tests.mesh.test_birdcage_implant_rod import ROD_RADIUS  # noqa: E402
from tests.validation.test_imp8_rod_line_source_closed_form import (  # noqa: E402
    B1_BAND,
    L,
    MM,
    _phi_midpoint,
    _phi_tips,
    _point,
)
from tests.validation.test_imp8_rod_line_source_duty_closed_form import (  # noqa: E402
    DUTIES,
    POINTS,
    T_GRID,
    T_STAR,
    _theta_step,
)
from tests.validation.test_imp8_rod_line_source_duty_closed_form import (  # noqa: E402
    _convolution as _convolution_line,
)
from tests.validation.test_imp8_rod_line_source_duty_closed_form import (  # noqa: E402
    _hist as _hist_line,
)
from tests.validation.test_imp8_rod_line_source_duty_closed_form import (  # noqa: E402
    _theta_ss as _theta_ss_line,
)
from tests.validation.test_imp8_rod_tips_duty_closed_form import (  # noqa: E402
    P_DOSE_W,
    _convolution as _convolution_tips,
    _duhamel,
    _hist as _hist_tips,
    _phi_tips_dot,
    _step,
    _theta_ss as _theta_ss_tips,
)
from tests.validation.test_mat13_transient_perfusion_sphere import INTERIOR_BAND  # noqa: E402
from tests.validation.test_mat15_finite_cylinder_fem import _max_imag  # noqa: E402
from tests.validation.test_mat21_box_closed_form import T_INF as T_INF_0  # noqa: E402
from tests.validation.test_mat21_box_dose_closed_form import AVG_BAR, DOSE_BAR, W_B2  # noqa: E402
from tests.validation.test_mat21_box_dt_duty125_dose_closed_form import (  # noqa: E402
    _q_check_g,
    _win_check_g,
)
from tests.validation.test_mat21_wall_duty_closed_form import D1_BAND, DUH_BAND  # noqa: E402

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"  # committed, EX-57
BASENAME = "implants_22_rod_duty_cycle_window_and_dose"
HALF = 0.060  # m, half-width of the viewing grid on each axis
N_CELLS = 41  # per axis; odd on a symmetric box, so no vertex lies on the line
N_FIELD = 64  # midpoint-sum points for the line kernel on the grid
IM_BAND = 1e-12
REPRO_BAND = 1e-9  # the brief's: twice the .9e print's worst half-ulp, with margin
T_SAFE = 43.0  # C, the precondition _q_check_g's R_below-only quadrature needs

# theta(t*) per watt, the logs' .9e prints, keyed (duty, point).
# Line: 20261103T211919Z_IMP-8-item787-W1.log :54-59.  Tips: 20261104T103134Z_IMP-8-item811-W1.log :55-60.
LOG_LINE = {
    (1.0, POINTS[0]): 9.438476919e+00, (1.0, POINTS[1]): 4.645356263e+00,
    (0.5, POINTS[0]): 2.096144980e+00, (0.5, POINTS[1]): 2.241722225e+00,
    (0.25, POINTS[0]): 8.986021128e-01, (0.25, POINTS[1]): 1.004369893e+00,
}
LOG_TIPS = {
    (1.0, POINTS[0]): 3.000946732e+01, (1.0, POINTS[1]): 1.927919823e+00,
    (0.5, POINTS[0]): 2.474673564e+00, (0.5, POINTS[1]): 1.137209732e+00,
    (0.25, POINTS[0]): 9.807666287e-01, (0.25, POINTS[1]): 5.680057583e-01,
}
# Printed records at (2, 15) mm beside the new lines: (six-minute max K/W, CEM43 s), by duty.
REC_LINE = {1.0: (8.948814, 1.571787e-06), 0.5: (5.135368, 3.548044e-08), 0.25: (2.656731, 1.816366e-09)}
REC_TIPS = {1.0: (29.54154, 0.005531217), 0.5: (17.81965, 0.0002971548), 0.25: (9.05116, 1.651446e-05)}


def _say(comm, msg):
    if comm.rank == 0:
        print(msg, flush=True)


def _nearest(comm, xyz, target):
    """(owner-reduced) index of the owned vertex nearest `target`, its distance (m)."""
    d = np.linalg.norm(xyz - np.asarray(target, dtype=float) * MM, axis=1)
    i = int(np.argmin(d)) if d.size else -1
    loc = float(d[i]) if i >= 0 else np.inf
    dmin = comm.allreduce(loc, op=MPI.MIN)
    me = comm.rank if (i >= 0 and loc == dmin) else comm.size
    owner = comm.allreduce(me, op=MPI.MIN)
    return (i if comm.rank == owner else -1), dmin


def _line_kernel(r, z, tau):
    """The line's step-response kernel, the gate's midpoint sum with the imported `_point`,
    written so it broadcasts over arrays (`_phi_midpoint` is scalar)."""
    zp = L * ((np.arange(N_FIELD) + 0.5) / N_FIELD - 0.5)
    big_r = np.hypot(np.asarray(r, float)[..., None], np.asarray(z, float)[..., None] - zp)
    return L / N_FIELD * np.sum(_point(big_r, tau), axis=-1)


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
    _say(comm, "EX-181 -- a rod under a duty cycle, as a line and as two tips: the six-minute "
               "window and the dose (IMP-8 steps 4 and 7, closed form)")
    _say(comm, "=" * 78)
    _say(comm, f"[gate] bands imported: D1 {D1_BAND:g}, Duh {DUH_BAND:g}, dose {DOSE_BAR:g}, "
               f"window {AVG_BAR:g}, c-cont {INTERIOR_BAND:g}, B1 {B1_BAND:.0e}; t* {T_STAR:g} s, "
               f"duties {DUTIES}, tips dose at {P_DOSE_W:g} W, rod radius {ROD_RADIUS * 1e3:g} mm")

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
        title="imp:22 -- a pulsed rod as a line and as two tips, read on a grid: y = 0 slice",
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

    # ---- the grid fields: theta at t*, D = 0.5, per watt ---------------------------------
    th_line = _duhamel(_line_kernel, r, z, T_STAR, 0.5)
    th_tips = _duhamel(_phi_tips, r, z, T_STAR, 0.5)
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

    # ---- (Node): vertex arrays against the scalar forms ----------------------------------
    i_node, d_node = _nearest(comm, xyz, (20.0, 0.0, 0.0))
    gap_line = gap_tips = 0.0
    node_xyz = None
    if i_node >= 0:
        rn, zn = float(r[i_node]), float(z[i_node])
        ref_line = _duhamel(lambda a, b, t: _phi_midpoint(a, b, t, N_FIELD), rn, zn, T_STAR, 0.5)
        ref_tips = _duhamel(_phi_tips, rn, zn, T_STAR, 0.5)
        gap_line = abs(th_line[i_node] / ref_line - 1.0)
        gap_tips = abs(th_tips[i_node] / ref_tips - 1.0)
        node_xyz = xyz[i_node]
    gap_line = comm.allreduce(gap_line, op=MPI.MAX)
    gap_tips = comm.allreduce(gap_tips, op=MPI.MAX)
    node_xyz = comm.bcast(node_xyz, root=comm.allreduce(
        comm.rank if node_xyz is not None else comm.size, op=MPI.MIN))

    # ---- the histories and their gates, collected ----------------------------------------
    kernels = {
        "line": dict(hist=_hist_line, ss=_theta_ss_line, p_dose=1.0, log=LOG_LINE,
                     step=lambda rr, zz, t: _theta_step(rr, zz, t),
                     conv=lambda rr, zz, t, d: _convolution_line(rr, zz, t, d)),
        "tips": dict(hist=_hist_tips, ss=_theta_ss_tips, p_dose=P_DOSE_W, log=LOG_TIPS,
                     step=lambda rr, zz, t: _step(_phi_tips, rr, zz, t),
                     conv=lambda rr, zz, t, d: _convolution_tips(_phi_tips_dot, rr, zz, t, d)),
    }
    res = {}  # (kernel, duty, pt) -> dict
    d1 = {}
    duh_worst = {k: 0.0 for k in kernels}
    qd_worst = {k: 0.0 for k in kernels}
    wd_worst = {k: 0.0 for k in kernels}
    repro_worst = 0.0
    window_ok = True
    t_ok = True
    t_max_all = {k: -np.inf for k in kernels}
    for name, K in kernels.items():
        for pt in POINTS:
            rr, zz = pt[0] * MM, pt[1] * MM
            step = np.array([K["step"](rr, zz, t) for t in T_GRID])
            d1[(name, pt)] = float(np.max(np.abs(K["hist"](1.0, pt) - step))) / K["ss"](pt)
        for duty in DUTIES:
            for pt in POINTS:
                rr, zz = pt[0] * MM, pt[1] * MM
                th = K["hist"](duty, pt)
                conv = np.array([K["conv"](rr, zz, t, duty) for t in T_GRID])
                e_duh = float(np.max(np.abs(th - conv))) / K["ss"](pt)
                duh_worst[name] = max(duh_worst[name], e_duh)
                T = T_INF_0 + K["p_dose"] * th
                t_max_all[name] = max(t_max_all[name], float(np.max(T)))
                t_ok = t_ok and float(np.max(T)) < T_SAFE
                dose, _ref, rel, _swap, _miss = _q_check_g(T, T_GRID)
                mx, t_mx, lo, hi, _wref, wrel = _win_check_g(th, T_GRID)
                qd_worst[name] = max(qd_worst[name], rel)
                wd_worst[name] = max(wd_worst[name], wrel)
                window_ok = window_ok and (lo <= mx <= hi and T_GRID[0] + W_B2 <= t_mx <= T_GRID[-1])
                repro = abs(th[-1] / K["log"][(duty, pt)] - 1.0)
                repro_worst = max(repro_worst, repro)
                res[(name, duty, pt)] = dict(th=float(th[-1]), mx=float(mx), dose=float(dose),
                                             t_mx=float(t_mx), repro=repro)
    # (c-cont): the continuous form (D = 1) against the D = 0.5 history at t*, factor 1
    ccont = {}
    for name, K in kernels.items():
        ccont[name] = min(abs(K["hist"](1.0, pt)[-1] / K["hist"](0.5, pt)[-1] - 1.0) for pt in POINTS)

    # ---- printed records ------------------------------------------------------------------
    p0 = POINTS[0]
    for duty in DUTIES:
        a, b = res[("line", duty, p0)], res[("tips", duty, p0)]
        _say(comm, f"[print] D={duty:g} (2, 15) mm: line theta(t*) {a['th']:.7g} K/W, six-min max "
                   f"{a['mx']:.7g} K/W, CEM43 {a['dose']:.4g} s at 1 W (item 787: {REC_LINE[duty][0]:.7g}, "
                   f"{REC_LINE[duty][1]:.4g})")
        _say(comm, f"[print] D={duty:g} (2, 15) mm: tips theta(t*) {b['th']:.7g} K/W, six-min max "
                   f"{b['mx']:.7g} K/W, CEM43 {b['dose']:.4g} s at {P_DOSE_W:g} W (item 811: "
                   f"{REC_TIPS[duty][0]:.7g}, {REC_TIPS[duty][1]:.4g}); tips / line theta(t*) "
                   f"{b['th'] / a['th']:.7g}")
    _say(comm, "[print] the two doses are at different powers (1 W line, "
               f"{P_DOSE_W:g} W tips); no dose ratio is printed")
    for duty in DUTIES:
        _say(comm, f"[print] D={duty:g} (10, 0) mm: theta(t*) line {res[('line', duty, POINTS[1])]['th']:.7g}, "
                   f"tips {res[('tips', duty, POINTS[1])]['th']:.7g} K/W")
    _say(comm, f"[grid] min r on the grid {r_min * 1e3:.4f} mm (rod radius {ROD_RADIUS * 1e3:g} mm); "
               f"D = 0.5 theta tips / line at t*: min {r_lo:.7g}, max {r_hi:.7g}")
    _say(comm, f"[node] (20, 0, 0) mm probe: vertex {np.round(node_xyz * 1e3, 4).tolist()} mm "
               f"(distance {d_node * 1e3:.4f} mm); line gap {gap_line:.3e}, tips gap {gap_tips:.3e}")
    _say(comm, f"[D1] line {d1[('line', POINTS[0])]:.3e} / {d1[('line', POINTS[1])]:.3e}; tips "
               f"{d1[('tips', POINTS[0])]:.3e} / {d1[('tips', POINTS[1])]:.3e} (band {D1_BAND:g})")
    _say(comm, f"[Duh] worst line {duh_worst['line']:.3e}, tips {duh_worst['tips']:.3e} (band {DUH_BAND:g})")
    _say(comm, f"[Q-d] worst rel line {qd_worst['line']:.3e}, tips {qd_worst['tips']:.3e} "
               f"(bar {DOSE_BAR:g}); [Win-d] worst rel line {wd_worst['line']:.3e}, tips "
               f"{wd_worst['tips']:.3e} (bar {AVG_BAR:g}); max T line {t_max_all['line']:.5f} C, "
               f"tips {t_max_all['tips']:.5f} C (< {T_SAFE:g})")
    _say(comm, f"[repro] worst |theta(t*) / log - 1| over 12 cases: {repro_worst:.3e} (<= {REPRO_BAND:g})")
    _say(comm, f"[c-cont] D=1 vs D=0.5 miss at t*, least over the points: line {ccont['line']:.4e}, "
               f"tips {ccont['tips']:.4e} (> {INTERIOR_BAND:g})")
    _say(comm, f"max|Im theta| = {imag:.3e}")

    # ---- gates: collected above, asserted here --------------------------------------------
    _say(comm, f"[gates] (Geom) {r_min * 1e3:.4f} mm > {ROD_RADIUS * 1e3:g}; (Node) {gap_line:.3e} / "
               f"{gap_tips:.3e} <= {B1_BAND:.0e}; (c-cont) {min(ccont.values()):.4e} > {INTERIOR_BAND:g}")
    assert r_min > ROD_RADIUS, f"(Geom) a vertex at r = {r_min:.4e} m is inside the rod"
    for k in kernels:
        for pt in POINTS:
            assert d1[(k, pt)] <= D1_BAND, f"(D1) {k} {pt}: {d1[(k, pt)]:.3e}"
        assert duh_worst[k] <= DUH_BAND, f"(Duh) {k}: {duh_worst[k]:.3e}"
        assert t_max_all[k] < T_SAFE, f"max T {k}: {t_max_all[k]:.5f}"
        assert qd_worst[k] <= DOSE_BAR, f"(Q-d) {k}: {qd_worst[k]:.3e}"
        assert wd_worst[k] <= AVG_BAR, f"(Win-d) {k}: {wd_worst[k]:.3e}"
        assert ccont[k] > INTERIOR_BAND, f"(c-cont) {k} did not fire: {ccont[k]:.4e}"
    assert window_ok, "(Win-d) six-minute max outside its bracket or window"
    assert t_ok
    assert repro_worst <= REPRO_BAND, f"(repro) {repro_worst:.3e}"
    assert gap_line <= B1_BAND, f"(Node) line {gap_line:.3e}"
    assert gap_tips <= B1_BAND, f"(Node) tips {gap_tips:.3e}"
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
    _say(comm, f"\n[output] {xdmf}.xdmf -- theta(t* = {T_STAR:g} s) per watt under D = 0.5 for both "
               f"shapes and their ratio, P1 on the {ncell}-cell viewing grid")
    _say(comm, f"\n[timing] total {time.perf_counter() - started:.1f} s at -n {comm.size}")
    _say(comm, "\nAll imported gates hold.")


if __name__ == "__main__":
    main()
