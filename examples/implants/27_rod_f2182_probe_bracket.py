"""Example (`EX-190`): the F2182-shaped probe on the rod bracket -- gel, titanium and hole rods
read 1 and 2 mm beyond the tip, continuous and at D = 0.5.

`imp:26` puts the finite rod's material bracket on one coarse box but reads it at (2, 0) and
(2, 15) mm only. Here the same box (64 140 cells), the same source (`IMP-8` I55's DG0 line source
on the rod's cells, P = 1 W), the same stepping (degree 2, Crank-Nicolson with two backward-Euler
half-steps, ``N_MAIN`` steps to ``T_PRINT`` = 900 s, Robin ``H_GEL`` on all six faces, from
theta_0 = 0) is run three times, with ``snapshot_steps`` = (10, 20, 30):

1. rod / gel  -- the body cells at the gel pair (`imp:26`'s run 4);
2. rod / Ti   -- the body cells at the nominal Ti pair (`IMP-8` I73; `imp:26`'s run 5);
3. rod / hole -- the body cells at ``EPS_HOLE`` x gel k and rho c (`IMP-8` I75's hole arm; new on
                 this box).

The probe (`IMP-8` I75) is on the axis at (0, +-16) mm, 1 mm beyond this flat-ended rod, and at
(0, +-17) mm, where "1 mm beyond the apex" lands on a clinical capsule rod 2 mm beyond the same
line. Placement and the 900 s are quoted from memory. D = 0.5 is I75's superposition
theta(900) - theta(675) + theta(450) - theta(225), not a switched solve. I78's infinite-gel closed
forms (line and twin) are printed beside.

The box is a **disclosed local copy** of `imp:26`'s ``_coarse_box`` (itself a copy of `imp:23`'s;
examples do not import examples), with ``_body_cells`` and the figure's tag copy local too.
Everything else is imported by name from the gate modules (the `ANS-1` rule); ``box`` is a pytest
fixture and is not imported. No gate edit, no ``src/`` edit.

Asserted (collected first): the cell count equals ``W0_CELLS`` and is <= ``CELL_LIMIT``; the body
count equals ``BODY_CELLS_IMP26`` and is > 0; (E) <= ``E_BAND`` on all three runs; (S0)
|int q - P|/P <= ``B_BAND``; max|Im theta| <= 1e-12 (field and probe reads); no non-finite
coefficient dof; the reproduction: rod / gel and rod / Ti theta(2, 0) and theta(2, 15) equal
`imp:26`'s closing-log prints within one unit of the sixth significant digit; the coefficient
control: the Ti and hole arms' theta(2, 0) are off the gel run's by more than ``CONTROL_FLOOR``;
I75's control: theta_probe(900 s) > theta_probe(675 s) at (0, +-16) mm on every arm; I78's
placement control on the FEM: 1 - theta(0, +-17)/theta(0, +-16) > ``INTERIOR_BAND`` on every arm.

Printed only: the probe table per arm, the closed forms and the records beside, the ordering.

Run it::

    ./scripts/run_examples.sh -e imp:27 -n 2 -t 480

Output: ``examples/implants/paraview_output/implants_27_rod_f2182_probe_bracket_combined.xdmf``
(the three theta fields at t*, P1, and the DG0 k of the Ti and hole arms).
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

from fem_em_solver.io.mesh import _model_to_mesh  # noqa: E402
from fem_em_solver.io.paraview_utils import (  # noqa: E402
    adopt_host_ownership,
    write_xdmf_with_tags,
)
from fem_em_solver.post.setup_figure import write_setup_figure  # noqa: E402
from fem_em_solver.thermal import solve_transient_conduction  # noqa: E402

from tests.mesh.test_birdcage_implant_rod import ROD_RADIUS  # noqa: E402
from tests.validation.test_imp8_composite_cylinder_closed_form import (  # noqa: E402
    ROD_K,
    ROD_RHO_C,
)
from tests.validation.test_imp8_probe_f2182_closed_form import (  # noqa: E402
    SHAPES,
    Z_C,
    _reads,
)
from tests.validation.test_imp8_probe_f2182_closed_form import Z_P as Z_P_M  # noqa: E402
from tests.validation.test_imp8_rod_conducting_body_fem import (  # noqa: E402
    CONTROL_FLOOR,
    EPS_HOLE,
    _coefficients,
)
from tests.validation.test_imp8_rod_line_source_closed_form import L, MM, T_PRINT  # noqa: E402
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
    T_INF,
    _cells_global,
    _h_near,
    _rod_source,
    _theta_pts,
)
from tests.validation.test_imp8_rod_probe_f2182_fem import (  # noqa: E402
    SNAP_STEPS,
    _sixth,
    _superpose,
)
from tests.validation.test_imp8_rod_probe_f2182_fem import Z_P as Z_P_MM  # noqa: E402
from tests.validation.test_imp8_rod_tips_fem import RATIO_PTS  # noqa: E402
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
BASENAME = "implants_27_rod_f2182_probe_bracket"

#: `imp:26`'s box: SIZE_SCALE 2.0 gives 64 140 cells at -n 2 (`20261108T210314Z_EX-188-W0.log`).
SIZE_SCALE = 2.0
W0_CELLS = 64_140
CELL_LIMIT = 120_000
#: `imp:26`'s body-cell count (`20261108T211230Z_EX-188.log:55`).
BODY_CELLS_IMP26 = 833
IMAG_BAND = 1e-12

#: `imp:26`'s closing-log prints (`20261108T211230Z_EX-188.log:70`, `:72`), K/W at 900 s:
#: (theta(2, 0), theta(2, 15)).
REC_IMP26 = {"gel": (17.7246, 9.77473), "ti": (17.1764, 10.6763)}
#: I75's 485 291-cell probe means at (0, +-16) mm, K/W: (continuous, D = 0.5) and the log.
REC_I75 = {
    "ti": (1.057039e01, 2.104578e00, "20261108T194323Z_IMP-8-item881-W1.log:122"),
    "gel": (9.161894e00, 1.975098e00, "20261108T195853Z_IMP-8-item881-W2b.log:122"),
    "hole": (8.783168e00, 1.891937e00, "20261108T200636Z_IMP-8-item881-W3.log:122"),
}
REC_I78 = "20261108T133157Z_IMP-8-item874-W1.log:58-62"


def _say(comm, msg):
    if comm.rank == 0:
        print(msg, flush=True)


def _coarse_box(comm, scale):
    """Disclosed local copy of `imp:26`'s ``_coarse_box`` (itself a copy of `imp:23`'s, of the
    gates' ``_graded_box`` with ``H_IN``, ``H_MID``, ``H_OUT`` multiplied by ``scale``)."""
    import gmsh

    rank = 0
    if comm.rank == rank:
        gmsh.initialize()
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("imp27_coarse_box")
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


def _body_cells(msh, q):
    """Local copy of `imp:26`'s: cell indices (owned + ghost) of q's positive DG0 dofs."""
    tdim = msh.topology.dim
    n_loc = msh.topology.index_map(tdim).size_local
    n_all = n_loc + msh.topology.index_map(tdim).num_ghosts
    dofs = np.asarray(q.function_space.dofmap.list).reshape(-1)[:n_all]
    body = np.flatnonzero(np.real(q.x.array[dofs]) > 0.0).astype(np.int32)
    n_body = msh.comm.allreduce(int(np.count_nonzero(body < n_loc)), op=MPI.SUM)
    return body, n_body


def _p1_theta(msh, T):
    V1 = fem.functionspace(msh, ("Lagrange", 1))
    tmp = fem.Function(V1)
    tmp.interpolate(T)
    out = fem.Function(V1)
    out.x.array[:] = np.real(tmp.x.array) - T_INF
    out.x.scatter_forward()
    return out


def _real_dg0(k):
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
    _say(comm, "EX-190 -- the F2182-shaped probe on the rod bracket (IMP-8 I55, I73, I75, I78)")
    _say(comm, "=" * 78)
    _say(comm, f"[gate] imported: E_BAND {E_BAND:.0e}, B_BAND {B_BAND:.0e}, INTERIOR_BAND "
               f"{INTERIOR_BAND:.0e}, CONTROL_FLOOR {CONTROL_FLOOR:.0e}, N_MAIN {N_MAIN}, "
               f"T_PRINT {T_PRINT:g} s, DEGREE {DEGREE}, ROD_RADIUS {ROD_RADIUS / MM:g} mm, "
               f"L {L / MM:g} mm, SNAP_STEPS {SNAP_STEPS}")
    _say(comm, f"[gate] ROD_K {ROD_K:g} W/(m K), ROD_RHO_C {ROD_RHO_C:.6g} J/(m3 K) (nominal); "
               f"K_GEL {K_GEL:g}, RHO_C {RHO_C:.6g}; EPS_HOLE {EPS_HOLE:.0e}")
    z_p, z_c = float(Z_P_MM), float(Z_C / MM)
    _say(comm, f"[probe] I75 Z_P {z_p:g} mm (I78 Z_P {Z_P_M / MM:g} mm); I78 Z_C {z_c:g} mm; "
               f"placement and 900 s quoted from memory")
    probe_pts = ((0.0, z_p), (0.0, -z_p), (0.0, z_c), (0.0, -z_c))

    # ---- mesh, source, body -------------------------------------------------------------
    t0 = time.perf_counter()
    msh, ct, ft = _coarse_box(comm, SIZE_SCALE)
    n = _cells_global(msh)
    t_build = time.perf_counter() - t0
    _say(comm, f"[mesh] cells {n} (W0 {W0_CELLS}; limit {CELL_LIMIT}); ranks {comm.size}; "
               f"t_build {t_build:.1f} s")

    q_l, n_sel_l, v_sel_l = _rod_source(msh)
    s0_l = abs(float(np.real(_allsum(comm, q_l * ufl.dx))) - P) / P
    body_idx, n_body = _body_cells(msh, q_l)
    _say(comm, f"[S0] line |int q - P|/P {s0_l:.3e} (B_BAND {B_BAND:.0e}); selected cells "
               f"{n_sel_l}; body cells {n_body} (imp:26 {BODY_CELLS_IMP26})")

    n_loc = msh.topology.index_map(msh.topology.dim).size_local
    fig_vals = np.ones(n_loc, dtype=np.int32)
    fig_vals[body_idx[body_idx < n_loc]] = 2
    ct_fig = dmesh.meshtags(msh, msh.topology.dim, np.arange(n_loc, dtype=np.int32), fig_vals)
    write_setup_figure(
        msh,
        ct_fig,
        FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names={1: "phantom (gel box)", 2: "rod body cells (Ti / hole arms)"},
        title="imp:27 -- coarse gel box and rod body: y = 0 slice through the axis; probes on "
              "the axis at z = +-16 and +-17 mm",
        translucent_tags=(1,),
        slice_normal=(0.0, 1.0, 0.0),
        slice_origin=(0.0, 0.0, 0.0),
        comm=comm,
    )

    pts = GATE_PTS + RATIO_PTS  # imp:26's batch (a), unchanged
    k15, k0 = pts.index((2.0, 15.0)), pts.index((2.0, 0.0))
    for (r, z), h in zip(probe_pts, _h_near(msh, probe_pts)):
        _say(comm, f"[mesh] h within 3 mm of ({r:g}, {z:+g}) mm: {h / MM:.3f} mm")

    dt = T_PRINT / N_MAIN
    kw = dict(h=H_GEL, t_inf=T_INF, t_init=T_INF, robin_tag=ROBIN_TAG, dt=dt, n_steps=N_MAIN,
              theta=0.5, startup_half_steps=2, degree=DEGREE, snapshot_steps=SNAP_STEPS)

    k_ti, rc_ti, bad_ti = _coefficients(msh, body_idx, "ti")
    k_ho, rc_ho, bad_ho = _coefficients(msh, body_idx, "hole")
    n_bad = bad_ti + bad_ho
    _say(comm, f"[coef] non-finite coefficient dofs {n_bad}")

    plan = (
        ("gel", "rod / gel", dict(k=K_GEL, rho_c=RHO_C)),
        ("ti", "rod / Ti", dict(k=k_ti, rho_c=rc_ti)),
        ("hole", "rod / hole", dict(k=k_ho, rho_c=rc_ho)),
    )
    runs = {}
    for key, label, coef in plan:
        t1 = time.perf_counter()
        res = solve_transient_conduction(msh, ft, q=q_l, **coef, **kw)
        call = time.perf_counter() - t1
        e = float(np.real(res.max_balance_residual))
        th = _theta_pts(res.temperature, pts)  # (a): one batch
        fields_t = {N_MAIN: res.temperature, **{j: res.snapshots[j] for j in SNAP_STEPS}}
        probe = {j: np.asarray(_theta_pts(f, probe_pts)) for j, f in sorted(fields_t.items())}
        runs[key] = dict(label=label, res=res, E=e, call=call, th=th, probe=probe,
                         nfac=res.n_factorizations)
        _say(comm, f"[{key}] t_call {call:.1f} s; n_factorizations {res.n_factorizations}; "
                   f"(E) {e:.3e} (E_BAND {E_BAND:.0e})")

    # ---- the reproduction (anchor) ------------------------------------------------------
    repro_ok = {}
    for key in ("gel", "ti"):
        a0, a15 = np.real(runs[key]["th"][k0]), np.real(runs[key]["th"][k15])
        r0, r15 = REC_IMP26[key]
        d0, d15 = abs(a0 - r0), abs(a15 - r15)
        ok = d0 <= _sixth(r0) and d15 <= _sixth(r15)
        repro_ok[key] = ok
        _say(comm, f"[repro] {runs[key]['label']}: theta(2,0) {a0:.7g} vs imp:26 {r0:.6g} "
                   f"(|diff| {d0:.2e}, 6th-digit unit {_sixth(r0):.0e}); theta(2,15) {a15:.7g} vs "
                   f"{r15:.6g} (|diff| {d15:.2e}, unit {_sixth(r15):.0e}); within: {ok}")

    # ---- the coefficient control --------------------------------------------------------
    moves = {}
    g0 = np.real(runs["gel"]["th"][k0])
    for arm in ("ti", "hole"):
        mv = np.real(runs[arm]["th"][k0]) / g0 - 1.0
        moves[arm] = float(mv)
        _say(comm, f"[control] {runs[arm]['label']} over rod / gel at (2, 0) mm: theta ratio - 1 "
                   f"= {mv:+.4e} (CONTROL_FLOOR {CONTROL_FLOOR:.0e})")

    # ---- the probe table ----------------------------------------------------------------
    probe_imag = 0.0
    tab = {}
    _say(comm, "\n[table] probe theta (K/W per rod, t* = %g s); points (0, +16), (0, -16), "
               "(0, +17), (0, -17) mm" % T_PRINT)
    for key, run in runs.items():
        pr = run["probe"]
        for arr in pr.values():
            probe_imag = max(probe_imag, float(np.max(np.abs(np.imag(arr)))))
        p9 = np.real(pr[N_MAIN])
        p675 = np.real(pr[30])
        duty = [float(np.real(_superpose({j: pr[j][i] for j in pr}))) for i in range(4)]
        m16, m17 = 0.5 * (p9[0] + p9[1]), 0.5 * (p9[2] + p9[3])
        d16, d17 = 0.5 * (duty[0] + duty[1]), 0.5 * (duty[2] + duty[3])
        drop = [1.0 - p9[2] / p9[0], 1.0 - p9[3] / p9[1]]
        grow = [p9[0] / p675[0], p9[1] / p675[1]]
        tab[key] = dict(p9=p9, m16=m16, m17=m17, d16=d16, d17=d17, drop=drop, grow=grow)
        _say(comm, f"[table] {run['label']:10s} theta(900): {p9[0]:.7g} {p9[1]:.7g} {p9[2]:.7g} "
                   f"{p9[3]:.7g}; mean 16 mm {m16:.7g}, mean 17 mm {m17:.7g}")
        _say(comm, f"[table] {run['label']:10s} D=0.5 superposed means: 16 mm {d16:.7g}, 17 mm "
                   f"{d17:.7g}; theta(900)/theta(675) at +16 {grow[0]:.7g}, -16 {grow[1]:.7g}; "
                   f"theta(17)/theta(16) {m17 / m16:.7g}; 1 - ratio +/- {drop[0]:.7g} / "
                   f"{drop[1]:.7g} (INTERIOR_BAND {INTERIOR_BAND:.0e})")
        c16, d16r, _ = REC_I75[key][0], REC_I75[key][1], None
        _say(comm, f"[table]   I75 485 291-cell mean at 16 mm: continuous {c16:.7g}, D=0.5 "
                   f"{d16r:.7g} ({REC_I75[key][2]})")
    for name, (kern, _kd, steady) in SHAPES.items():
        for tag, z in (("16", Z_P_M), ("17", Z_C)):
            cont, duty, ss = _reads(kern, steady, z)
            _say(comm, f"[I78] {name} infinite gel, 1 W per rod, {tag} mm: theta(900) {cont:.7g}; "
                       f"D=0.5 {duty:.7g}; steady {ss:.7g} ({REC_I78})")
    order = " > ".join(nm for nm, _ in sorted(
        (("Ti", tab["ti"]["m16"]), ("gel", tab["gel"]["m16"]), ("hole", tab["hole"]["m16"])),
        key=lambda x: -x[1]))
    _say(comm, f"[order] probe mean at 16 mm, continuous: {order} (gates' mesh, item 881: "
               f"Ti > gel > hole)")

    # ---- output -------------------------------------------------------------------------
    fields = {}
    for key, run in runs.items():
        f = _p1_theta(msh, run["res"].temperature)
        f.name = f"theta_rod_{key}_t_star"
        fields[f.name] = f
    for nm, kf in (("k_rod_ti", k_ti), ("k_rod_hole", k_ho)):
        g = _real_dg0(kf)
        g.name = nm
        fields[nm] = g
    imag = max(max(_max_imag(run["res"].temperature) for run in runs.values()), probe_imag)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_combined", msh, ct, fields, comm=comm)
    adopt_host_ownership(OUTPUT_DIR, comm=comm)
    _say(comm, f"\n[output] {path} -- three theta fields (K per W, t* = {T_PRINT:g} s, P1) and "
               f"the DG0 k of the Ti and hole arms")
    _say(comm, f"max|Im theta| = {imag:.3e}")

    # ---- assertions, collected first ----------------------------------------------------
    gates = {
        "cells_eq_W0": n == W0_CELLS,
        "cells_le_limit": n <= CELL_LIMIT,
        "body_eq_imp26": n_body == BODY_CELLS_IMP26,
        "body_positive": n_body > 0 and n_sel_l > 0,
        "coef_finite": n_bad == 0,
        "S0_line": s0_l <= B_BAND,
        "max_imag": imag <= IMAG_BAND,
        "repro_gel": repro_ok["gel"],
        "repro_ti": repro_ok["ti"],
    }
    for key, run in runs.items():
        gates[f"E_{key}"] = run["E"] <= E_BAND
        gates[f"probe_growth_{key}"] = all(g > 1.0 for g in tab[key]["grow"])
        gates[f"placement_{key}"] = all(d > INTERIOR_BAND for d in tab[key]["drop"])
    for arm, mv in moves.items():
        gates[f"control_{arm}"] = abs(mv) > CONTROL_FLOOR
    _say(comm, f"[gates] {gates}")
    for key, run in runs.items():
        _say(comm, f"[timing] {key}: t_call {run['call']:.1f} s; n_factorizations {run['nfac']}")
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
