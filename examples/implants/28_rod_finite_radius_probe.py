"""Example (`EX-191`): the rod's finite radius at the probe -- gel and Ti beside the uniform-cylinder
closed form (`IMP-8` I83), read on the four ladder rungs.

`imp:27` sets I78's zero-radius line beside its probe reads and calls the line "the like-for-like
closed form". `IMP-8` I83 (item 911) showed that 1 mm from the end of a 1 mm-radius source the
line is the wrong comparand, and the uniform finite cylinder Phi_a is the right one. This page
runs `imp:27`'s box and source again (64 140 cells, I55's DG0 line source on the rod's cells,
P = 1 W; degree 2, Crank-Nicolson with two backward-Euler half-steps, ``N_MAIN`` steps to
``T_PRINT`` = 900 s, Robin ``H_GEL`` on all six faces, from theta_0 = 0), with
``snapshot_steps`` = (10, 20, 30), twice:

1. rod / gel -- the body cells at the gel pair;
2. rod / Ti  -- the body cells at the nominal Ti pair (`IMP-8` I66 / I73).

Two read batches per field (the batch rule): (a) `imp:27`'s probe batch, (0, +-16) and (0, +-17)
mm, unchanged; (b) the ladder's other two rungs, (0, +-18) and (0, +-20) mm, built from I80's
``RUNGS``. D = 0.5 is I75's superposition, not a switched solve. Beside each rung: I83's
SCALE*Phi_a (900 s, D = 0.5, steady), I78's line, and the radius factor Phi_a/Phi.

The box is a **disclosed local copy** of `imp:27`'s ``_coarse_box`` (examples do not import
examples), with ``_body_cells`` local too. Everything else is imported by name from the gate
modules (the `ANS-1` rule); ``box`` is a pytest fixture and is not imported. No gate edit, no
``src/`` edit.

Asserted (collected first): the cell count equals ``W0_CELLS`` and is <= ``CELL_LIMIT``; body count
833 (allreduced); (E) <= ``E_BAND`` on both runs; (S0) |int q - P|/P <= ``B_BAND``; max|Im theta| <=
1e-12; no non-finite coefficient dof; the reproduction: gel and Ti theta(900 s) at +-16, +-17 mm and
the D = 0.5 means at 16 / 17 mm equal `imp:27`'s closing-log prints within one unit of the sixth
significant digit; (Lim): I83's identity at every rung and sign, <= ``B1_BAND``; (C1) 1 - gel mean /
(SCALE*Phi_line) > ``INTERIOR_BAND`` at 16 / 17 mm; (C2) |Ti mean / (SCALE*Phi_a) - 1| >
``INTERIOR_BAND`` at 16 / 17 mm; I78's placement control on the FEM; I81's strict ordering over
the rungs at each sign.

Printed only, never gated: the probe table, the closed forms, I83's (G) read on this box (flagged
for the weekly), the first-drop decomposition.

Run it::

    ./scripts/run_examples.sh -e imp:28 -n 2 -t 480

Output: ``examples/implants/paraview_output/implants_28_rod_finite_radius_probe_combined.xdmf``
(the two theta fields at t*, P1, and the DG0 k of the Ti arm).
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
from tests.validation.test_imp4_pec_cylinder_closed_form import B1_BAND  # noqa: E402
from tests.validation.test_imp8_composite_cylinder_closed_form import (  # noqa: E402
    ROD_K,
    ROD_RHO_C,
)
from tests.validation.test_imp8_probe_cylinder_source_closed_form import (  # noqa: E402
    FEM as I83_FEM,
)
from tests.validation.test_imp8_probe_cylinder_source_closed_form import (  # noqa: E402
    _closed,
    _cyl_steady,
    _phi_cyl,
)
from tests.validation.test_imp8_probe_f2182_closed_form import (  # noqa: E402
    SHAPES,
    Z_C,
    _reads,
)
from tests.validation.test_imp8_probe_f2182_closed_form import Z_P as Z_P_M  # noqa: E402
from tests.validation.test_imp8_probe_ladder_closed_form import RUNGS, _lab  # noqa: E402
from tests.validation.test_imp8_rod_conducting_body_fem import (  # noqa: E402
    CONTROL_FLOOR,
    _coefficients,
)
from tests.validation.test_imp8_rod_line_source_closed_form import L, MM, T_PRINT  # noqa: E402
from tests.validation.test_imp8_rod_line_source_duty_closed_form import SCALE  # noqa: E402
from tests.validation.test_imp8_rod_line_source_fem import (  # noqa: E402
    AXIS_HALF,
    D_IN0,
    D_IN1,
    D_OUT1,
    DEGREE,
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
from tests.validation.test_imp8_rod_tips_duty_closed_form import _duhamel  # noqa: E402
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
BASENAME = "implants_28_rod_finite_radius_probe"

#: `imp:27`'s box: SIZE_SCALE 2.0 gives 64 140 cells at -n 2 (`20261109T105242Z_EX-190.log`).
SIZE_SCALE = 2.0
W0_CELLS = 64_140
CELL_LIMIT = 120_000
BODY_CELLS_IMP27 = 833  # `20261109T105242Z_EX-190.log:55`
IMAG_BAND = 1e-12

#: `imp:27`'s closing-log prints (`20261109T105242Z_EX-190.log:61`, `:62`, `:64`, `:65`), K/W at
#: 900 s: theta at +16, -16, +17, -17 mm, then the D = 0.5 superposed means at 16 / 17 mm.
REC_IMP27 = {
    "gel": ((9.638764, 8.866075, 6.919993, 6.477902), (1.970653, 1.814308)),
    "ti": ((11.02901, 10.10473, 7.739184, 7.178166), (2.100275, 1.898788)),
}
REC_IMP27_LOG = "20261109T105242Z_EX-190.log"


def _say(comm, msg):
    if comm.rank == 0:
        print(msg, flush=True)


def _coarse_box(comm, scale):
    """Disclosed local copy of `imp:27`'s ``_coarse_box`` (itself a copy of `imp:26`'s and `imp:23`'s,
    of the gates' ``_graded_box`` with ``H_IN``, ``H_MID``, ``H_OUT`` multiplied by ``scale``)."""
    import gmsh

    rank = 0
    if comm.rank == rank:
        gmsh.initialize()
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("imp28_coarse_box")
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
    """Local copy of `imp:27`'s: cell indices (owned + ghost) of q's positive DG0 dofs."""
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


def _drop(v):
    return 1.0 - v[1] / v[0]


def main() -> None:
    comm = MPI.COMM_WORLD
    started = time.perf_counter()

    if not np.issubdtype(np.dtype(default_scalar_type), np.complexfloating):
        raise RuntimeError(
            "this example expects the complex DolfinX build the runner's `imp:` group "
            "sources (source /usr/local/bin/dolfinx-complex-mode)"
        )

    _say(comm, "=" * 78)
    _say(comm, "EX-191 -- the rod's finite radius at the probe (IMP-8 I55, I66, I75, I78, I80, I83)")
    _say(comm, "=" * 78)
    _say(comm, f"[gate] imported: E_BAND {E_BAND:.0e}, B_BAND {B_BAND:.0e}, B1_BAND {B1_BAND:.0e}, "
               f"INTERIOR_BAND {INTERIOR_BAND:.0e}, CONTROL_FLOOR {CONTROL_FLOOR:.0e}, "
               f"N_MAIN {N_MAIN}, T_PRINT {T_PRINT:g} s, DEGREE {DEGREE}, "
               f"ROD_RADIUS {ROD_RADIUS / MM:g} mm, L {L / MM:g} mm, SNAP_STEPS {SNAP_STEPS}, "
               f"SCALE {SCALE:.7g} K/W per unit Phi")
    _say(comm, f"[gate] ROD_K {ROD_K:g} W/(m K), ROD_RHO_C {ROD_RHO_C:.6g} J/(m3 K) (nominal); "
               f"K_GEL {K_GEL:g}, RHO_C {RHO_C:.6g}")
    z_p, z_c = float(Z_P_MM), float(Z_C / MM)
    rung_mm = [float(z / MM) for z in RUNGS]
    _say(comm, f"[probe] I75 Z_P {z_p:g} mm (I78 Z_P {Z_P_M / MM:g} mm); I78 Z_C {z_c:g} mm; "
               f"I80 RUNGS {rung_mm} mm; placement and 900 s quoted from memory")
    pts_a = ((0.0, z_p), (0.0, -z_p), (0.0, z_c), (0.0, -z_c))  # imp:27's probe batch
    pts_b = ((0.0, rung_mm[2]), (0.0, -rung_mm[2]), (0.0, rung_mm[3]), (0.0, -rung_mm[3]))
    all_pts = pts_a + pts_b

    # ---- mesh, source, body -------------------------------------------------------------
    t0 = time.perf_counter()
    msh, ct, ft = _coarse_box(comm, SIZE_SCALE)
    n = _cells_global(msh)
    t_build = time.perf_counter() - t0
    _say(comm, f"[mesh] cells {n} (W0 {W0_CELLS}; limit {CELL_LIMIT}); ranks {comm.size}; "
               f"t_build {t_build:.1f} s")

    q_l, n_sel_l, _v_sel = _rod_source(msh)
    s0_l = abs(float(np.real(_allsum(comm, q_l * ufl.dx))) - P) / P
    body_idx, n_body = _body_cells(msh, q_l)
    _say(comm, f"[S0] line |int q - P|/P {s0_l:.3e} (B_BAND {B_BAND:.0e}); selected cells "
               f"{n_sel_l}; body cells {n_body} (imp:27 {BODY_CELLS_IMP27})")

    n_loc = msh.topology.index_map(msh.topology.dim).size_local
    fig_vals = np.ones(n_loc, dtype=np.int32)
    fig_vals[body_idx[body_idx < n_loc]] = 2
    ct_fig = dmesh.meshtags(msh, msh.topology.dim, np.arange(n_loc, dtype=np.int32), fig_vals)
    write_setup_figure(
        msh,
        ct_fig,
        FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names={1: "phantom (gel box)", 2: "rod body cells (Ti arm)"},
        title="imp:28 -- coarse gel box and rod body: y = 0 slice through the axis; probes on "
              "the axis at z = +-16, +-17, +-18 and +-20 mm",
        translucent_tags=(1,),
        slice_normal=(0.0, 1.0, 0.0),
        slice_origin=(0.0, 0.0, 0.0),
        comm=comm,
    )

    for (r, z), h in zip(all_pts, _h_near(msh, all_pts)):
        _say(comm, f"[mesh] h within 3 mm of ({r:g}, {z:+g}) mm: {h / MM:.3f} mm")

    dt = T_PRINT / N_MAIN
    kw = dict(h=H_GEL, t_inf=T_INF, t_init=T_INF, robin_tag=ROBIN_TAG, dt=dt, n_steps=N_MAIN,
              theta=0.5, startup_half_steps=2, degree=DEGREE, snapshot_steps=SNAP_STEPS)

    k_ti, rc_ti, n_bad = _coefficients(msh, body_idx, "ti")
    _say(comm, f"[coef] non-finite coefficient dofs {n_bad}")

    plan = (
        ("gel", "rod / gel", dict(k=K_GEL, rho_c=RHO_C)),
        ("ti", "rod / Ti", dict(k=k_ti, rho_c=rc_ti)),
    )
    runs = {}
    for key, label, coef in plan:
        t1 = time.perf_counter()
        res = solve_transient_conduction(msh, ft, q=q_l, **coef, **kw)
        call = time.perf_counter() - t1
        e = float(np.real(res.max_balance_residual))
        fields_t = {N_MAIN: res.temperature, **{j: res.snapshots[j] for j in SNAP_STEPS}}
        # (a) imp:27's batch unchanged, then (b) the ladder batch: one batch per call
        pa = {j: np.asarray(_theta_pts(f, pts_a)) for j, f in sorted(fields_t.items())}
        pb = {j: np.asarray(_theta_pts(f, pts_b)) for j, f in sorted(fields_t.items())}
        runs[key] = dict(label=label, res=res, E=e, call=call, pa=pa, pb=pb,
                         nfac=res.n_factorizations)
        _say(comm, f"[{key}] t_call {call:.1f} s; n_factorizations {res.n_factorizations}; "
                   f"(E) {e:.3e} (E_BAND {E_BAND:.0e})")

    # ---- per-arm tables: theta at (+z, -z) for the four rungs ----------------------------
    probe_imag = 0.0
    tab = {}
    for key, run in runs.items():
        for arr in list(run["pa"].values()) + list(run["pb"].values()):
            probe_imag = max(probe_imag, float(np.max(np.abs(np.imag(arr)))))
        fin = {j: np.concatenate([np.real(run["pa"][j]), np.real(run["pb"][j])])
               for j in run["pa"]}  # order: +16 -16 +17 -17 +18 -18 +20 -20
        sup = np.real(np.array([_superpose({j: fin[j][i] for j in fin}) for i in range(8)]))
        plus = np.array([fin[N_MAIN][0], fin[N_MAIN][2], fin[N_MAIN][4], fin[N_MAIN][6]])
        minus = np.array([fin[N_MAIN][1], fin[N_MAIN][3], fin[N_MAIN][5], fin[N_MAIN][7]])
        dplus, dminus = sup[[0, 2, 4, 6]], sup[[1, 3, 5, 7]]
        tab[key] = dict(fin=fin[N_MAIN], plus=plus, minus=minus, mean=0.5 * (plus + minus),
                        dmean=0.5 * (dplus + dminus))

    # ---- the reproduction (anchor 1) ----------------------------------------------------
    repro_ok = {}
    for key in ("gel", "ti"):
        ref_th, ref_d = REC_IMP27[key]
        got = tab[key]["fin"][:4]
        ok = all(abs(g - r) <= _sixth(r) for g, r in zip(got, ref_th))
        got_d = tab[key]["dmean"][:2]
        ok_d = all(abs(g - r) <= _sixth(r) for g, r in zip(got_d, ref_d))
        repro_ok[key] = ok and ok_d
        _say(comm, f"[repro] {runs[key]['label']}: theta(900) at +16, -16, +17, -17 "
                   + " ".join(f"{g:.7g} vs {r:.7g} (|diff| {abs(g - r):.2e})"
                              for g, r in zip(got, ref_th))
                   + f"; D=0.5 means 16 / 17 mm {got_d[0]:.7g} / {got_d[1]:.7g} vs "
                     f"{ref_d[0]:.7g} / {ref_d[1]:.7g}; within the sixth-digit unit: "
                     f"{ok} / {ok_d} ({REC_IMP27_LOG})")

    # ---- the closed forms beside (I83, I78) ---------------------------------------------
    kern, _kd, steady = SHAPES["line"]
    cf = []  # per rung: dict
    lim_worst = 0.0
    for k, z in enumerate(RUNGS):
        cyl, cyl_d = _closed(z)
        line, line_d, line_ss = _reads(kern, steady, z)
        cyl_ss = SCALE * _cyl_steady(z)
        for sg in (1.0, -1.0):
            ph, st = _phi_cyl(0.0, sg * z, np.inf), _cyl_steady(sg * z)
            rel = abs(ph / st - 1.0)
            lim_worst = max(lim_worst, rel)
        cf.append(dict(cyl=cyl, cyl_d=cyl_d, cyl_ss=cyl_ss, line=line, line_d=line_d,
                       line_ss=line_ss, rf=cyl / line))
        _say(comm, f"[I83] {_lab(k)}: SCALE*Phi_a 900 s {cyl:.7g}, D=0.5 {cyl_d:.7g}, steady "
                   f"{cyl_ss:.7g}; line 900 s {line:.7g}, D=0.5 {line_d:.7g}, steady "
                   f"{line_ss:.7g}; radius factor Phi_a/Phi {cyl / line:.7g}")
    _say(comm, f"[Lim] I83's identity |Phi_a(0,z,inf)/g-form - 1|, every rung and sign: worst "
               f"{lim_worst:.3e} (B1_BAND {B1_BAND:.0e})")

    # ---- the table ----------------------------------------------------------------------
    for key, run in runs.items():
        t = tab[key]
        for k in range(4):
            _say(comm, f"[table] {run['label']:10s} {_lab(k)}: theta(900) +z {t['plus'][k]:.7g}, "
                       f"-z {t['minus'][k]:.7g}; mean {t['mean'][k]:.7g}, |diff| "
                       f"{abs(t['plus'][k] - t['minus'][k]):.7g}; D=0.5 mean {t['dmean'][k]:.7g}")
    gel, ti = tab["gel"], tab["ti"]
    for k in range(4):
        _say(comm, f"[table] Ti/gel - 1 at rung {k + 1} ({rung_mm[k]:g} mm): "
                   f"{ti['mean'][k] / gel['mean'][k] - 1.0:+.6f}")

    # ---- I83's (G), printed only ---------------------------------------------------------
    for k in range(4):
        rc = gel["mean"][k] / cf[k]["cyl"] - 1.0
        rd = gel["dmean"][k] / cf[k]["cyl_d"] - 1.0
        i83c = float(I83_FEM["gel"]["mean"][k]) / cf[k]["cyl"] - 1.0
        i83d = float(I83_FEM["gel"]["dmean"][k]) / cf[k]["cyl_d"] - 1.0
        _say(comm, f"[G] I83's (G), flagged for the weekly; not gated here. Rung {k + 1} "
                   f"({rung_mm[k]:g} mm): gel mean / (SCALE*Phi_a) - 1 = {rc:+.3e}; D=0.5 mean / "
                   f"_duhamel - 1 = {rd:+.3e} (INTERIOR_BAND {INTERIOR_BAND:.0e}); I83 on the "
                   f"gates' mesh {i83c:+.3e} / {i83d:+.3e} ({I83_FEM['gel']['log']})")

    # ---- the first-drop decomposition on the means --------------------------------------
    dl = _drop([cf[0]["line"], cf[1]["line"]])
    dc = _drop([cf[0]["cyl"], cf[1]["cyl"]])
    dg, dti = _drop(gel["mean"][:2]), _drop(ti["mean"][:2])
    _say(comm, f"[drop] first drop (16 -> 17 mm), means: line {dl:.6f}, cylinder {dc:.6f}, "
               f"gel {dg:.6f}, Ti {dti:.6f}")
    _say(comm, f"[drop] decomposition: radius (cyl - line) {dc - dl:+.6f}; mesh (gel - cyl) "
               f"{dg - dc:+.6f}; rod (Ti - gel) {dti - dg:+.6f}")

    # ---- the controls -------------------------------------------------------------------
    c1 = [1.0 - gel["mean"][k] / cf[k]["line"] for k in range(2)]
    c2 = [abs(ti["mean"][k] / cf[k]["cyl"] - 1.0) for k in range(2)]
    for k in range(4):
        _say(comm, f"[C1] rung {k + 1} ({rung_mm[k]:g} mm): 1 - gel mean/(SCALE*Phi_line) "
                   f"{1.0 - gel['mean'][k] / cf[k]['line']:.6f}; [C2] |Ti mean/(SCALE*Phi_a) - 1| "
                   f"{abs(ti['mean'][k] / cf[k]['cyl'] - 1.0):.6f}; "
                   f"{'asserted > INTERIOR_BAND' if k < 2 else 'printed only'}")
    place_ok, order_ok = {}, {}
    for key in runs:
        t = tab[key]
        drops = [1.0 - t["plus"][1] / t["plus"][0], 1.0 - t["minus"][1] / t["minus"][0]]
        place_ok[key] = all(d > INTERIOR_BAND for d in drops)
        order_ok[key] = all(np.all(np.diff(v) < 0.0) for v in (t["plus"], t["minus"]))
        _say(comm, f"[control] {runs[key]['label']}: 1 - theta(17)/theta(16) +z {drops[0]:.7g}, "
                   f"-z {drops[1]:.7g} (INTERIOR_BAND {INTERIOR_BAND:.0e}); strict decrease over "
                   f"16 > 17 > 18 > 20 mm at each sign: {order_ok[key]}")

    # ---- output -------------------------------------------------------------------------
    fields = {}
    for key, run in runs.items():
        f = _p1_theta(msh, run["res"].temperature)
        f.name = f"theta_rod_{key}_t_star"
        fields[f.name] = f
    g = _real_dg0(k_ti)
    g.name = "k_rod_ti"
    fields[g.name] = g
    imag = max(max(_max_imag(run["res"].temperature) for run in runs.values()), probe_imag)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_combined", msh, ct, fields, comm=comm)
    adopt_host_ownership(OUTPUT_DIR, comm=comm)
    _say(comm, f"\n[output] {path} -- two theta fields (K per W, t* = {T_PRINT:g} s, P1) and "
               f"the DG0 k of the Ti arm")
    _say(comm, f"max|Im theta| = {imag:.3e}")

    # ---- assertions, collected first ----------------------------------------------------
    gates = {
        "cells_eq_W0": n == W0_CELLS,
        "cells_le_limit": n <= CELL_LIMIT,
        "body_eq_imp27": n_body == BODY_CELLS_IMP27,
        "body_positive": n_body > 0 and n_sel_l > 0,
        "coef_finite": n_bad == 0,
        "S0_line": s0_l <= B_BAND,
        "max_imag": imag <= IMAG_BAND,
        "repro_gel": repro_ok["gel"],
        "repro_ti": repro_ok["ti"],
        "Lim": lim_worst <= B1_BAND,
        "C1_k16": c1[0] > INTERIOR_BAND,
        "C1_k17": c1[1] > INTERIOR_BAND,
        "C2_k16": c2[0] > INTERIOR_BAND,
        "C2_k17": c2[1] > INTERIOR_BAND,
    }
    for key, run in runs.items():
        gates[f"E_{key}"] = run["E"] <= E_BAND
        gates[f"placement_{key}"] = place_ok[key]
        gates[f"ladder_order_{key}"] = order_ok[key]
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
