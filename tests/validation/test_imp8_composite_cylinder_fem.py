"""IMP-8 step 10 / §10 I65 (§9 item 837): the composite cylinder in FEM.

A Ti core 0 < r < a (``ROD_K``, ``ROD_RHO_C``, nominal) in perfect contact with
a gel shell a < r < b (``K_GEL``, ``RHO_C``), a uniform shell source
``Q_SHELL`` (1 W per metre), θ ≈ 0 at r = b through a Robin wall at
``H_DIR`` = 1e7 W/(m²·K) (``solve_transient_conduction`` takes Robin only;
predicted offset 2.8e-5 of θ_ss(a)), insulated end faces, length 10 b.  Degree
2 on a curved P2 mesh (affine cells at r = b would move θ_ss(a) by ≈ −1.7 %;
review ruling 2026-11-05 03:00).  Comparand: item 843's composite series
(``test_imp8_composite_cylinder_closed_form._theta``), imported unmoved.

Asserted: (I) at (a, 0) and ((a+b)/2, 0), at t* and at the steady limit, and
(I0) at (0, 0), t_c, all ≤ ``INTERIOR_BAND``; (E) ≤ ``E_BAND``; (S0)
∫q dx = Q_SHELL·V_shell,h ≤ ``B_BAND`` and ∫q dx(1) = 0.  Control (asserted,
factor 1 on ``INTERIOR_BAND``): the k_r = ``K_GEL`` comparand misses
θ_FEM(0, t_c) by more than the band (predicted ≈ 0.133, printed).  Windows by
``FEM_EM_IMP8_I65_PHASE`` = ``build`` (W0) / ``solve`` (W1).
"""

from __future__ import annotations

import os
import time

import numpy as np
import pytest
import ufl
from dolfinx import fem
from dolfinx import mesh as dmesh
from mpi4py import MPI

from fem_em_solver.io.mesh import _model_to_mesh
from fem_em_solver.post.evaluation import evaluate_vector_field_parallel
from fem_em_solver.thermal import solve_steady_conduction, solve_transient_conduction

from tests.validation.test_imp8_composite_cylinder_closed_form import (
    A_ROD,
    B_WALL,
    Q_SHELL,
    ROD_K,
    ROD_RHO_C,
    T_C,
    T_STAR,
    _theta,
    _theta_ss,
)
from tests.validation.test_imp8_rod_line_source_fem import _theta_pts
from tests.validation.test_mat10_sar_to_heat_sphere import B_BAND, ROBIN_TAG, _allsum, _say
from tests.validation.test_mat12_transient_sar_heating_sphere import E_BAND, RHO_C
from tests.validation.test_mat13_transient_perfusion_sphere import INTERIOR_BAND
from tests.validation.test_mat15_coil_driven_heating_phantom import K_GEL
from tests.validation.test_mat18_two_layer_sphere_fem import _dg0

MM = 1e-3
LEN = 10.0 * B_WALL  # 50 mm, centred on z = 0
H_DIR = 1e7  # W/(m²·K): the near-Dirichlet wall (ruling 2026-11-05 03:00)
CELL_CAP = 250_000
DEGREE = 2
N_STEPS = 80
SNAP = 16  # 16 · t*/80 = t_c
T_INF = 0.0
PRED_CTRL = 0.133  # unlogged review arithmetic, predicted (rule (e)), not asserted
PRED_CTRL_OTHERS = "<= 7.4e-3 (a, t_c), approx 4e-4 (a, t*), 0 (a, steady)"
PRED_ROBIN = 2.8e-5  # predicted Robin offset / theta_ss(a), printed
RA, RM = A_ROD / MM, 0.5 * (A_ROD + B_WALL) / MM  # 1 mm, 3 mm
READ_PTS = ((RA, 0.0), (RM, 0.0), (0.0, 0.0))
ZW = 0.25 * LEN / MM  # ±L/4 z-invariance witness
WIT_PTS = ((RA, ZW), (RA, -ZW))

PHASE_ENV = "FEM_EM_IMP8_I65_PHASE"
PHASES = ("build", "solve")


def _phase():
    raw = os.environ.get(PHASE_ENV, "").strip()
    if raw == "":
        return None
    if raw in PHASES:
        return raw
    raise ValueError(f"{PHASE_ENV} must be unset or one of {PHASES}, got {raw!r}")


def test_k0_phase_knob(monkeypatch):
    monkeypatch.delenv(PHASE_ENV, raising=False)
    assert _phase() is None
    for ph in PHASES:
        monkeypatch.setenv(PHASE_ENV, ph)
        assert _phase() == ph
    for bad in ("1", "all", "steady"):
        monkeypatch.setenv(PHASE_ENV, bad)
        with pytest.raises(ValueError):
            _phase()


def _composite_cylinder(comm, *, h_a=0.30e-3, h_b=2.0e-3):
    """gmsh OCC fragment of the radius-a and radius-b cylinders, length 10 b;
    core tag 1, shell tag 2 (bbox r_max); r = b lateral face ``ROBIN_TAG``,
    end faces untagged (Neumann); curved P2 (``setOrder(2)``)."""
    import gmsh

    rank = 0
    if comm.rank == rank:
        gmsh.initialize()
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("imp8_composite_cylinder")
        occ = gmsh.model.occ
        z0 = -0.5 * LEN
        c_in = occ.addCylinder(0.0, 0.0, z0, 0.0, 0.0, LEN, A_ROD)
        c_out = occ.addCylinder(0.0, 0.0, z0, 0.0, 0.0, LEN, B_WALL)
        occ.fragment([(3, c_out)], [(3, c_in)])
        occ.synchronize()
        core, shell = [], []
        for _, v in gmsh.model.getEntities(3):
            bb = gmsh.model.getBoundingBox(3, v)
            r_max = max(abs(bb[0]), abs(bb[1]), abs(bb[3]), abs(bb[4]))
            (core if r_max < 0.5 * (A_ROD + B_WALL) else shell).append(v)
        assert len(core) == 1 and len(shell) == 1, (core, shell)
        surfs = sorted({s for _, s in gmsh.model.getBoundary(
            [(3, core[0]), (3, shell[0])], combined=False, oriented=False)})
        outer, lateral_a = [], []
        for s in surfs:
            bb = gmsh.model.getBoundingBox(2, s)
            # z first: a planar end face (Neumann).  OCC pads bounding boxes by
            # its tolerance, so a 1e-9 test let the end annuli through and they
            # were tagged as the wall (W0 20261105T225535Z: a tagged node at
            # |r - b| = 4 mm, a tagged midpoint at |z| = L/2).
            if abs(bb[5] - bb[2]) < 1e-3 * LEN:
                continue
            r_max = max(abs(bb[0]), abs(bb[1]), abs(bb[3]), abs(bb[4]))
            if abs(r_max - B_WALL) < 1e-6:
                outer.append(s)
            elif abs(r_max - A_ROD) < 1e-6:
                lateral_a.append(s)
        assert outer and lateral_a, (outer, lateral_a)
        gmsh.model.addPhysicalGroup(3, core, tag=1)
        gmsh.model.addPhysicalGroup(3, shell, tag=2)
        gmsh.model.addPhysicalGroup(2, outer, tag=ROBIN_TAG)
        fm = gmsh.model.mesh.field
        fd = fm.add("Distance")
        fm.setNumbers(fd, "SurfacesList", lateral_a)
        fm.setNumber(fd, "Sampling", 200)
        th = fm.add("Threshold")
        fm.setNumber(th, "InField", fd)
        fm.setNumber(th, "SizeMin", h_a)
        fm.setNumber(th, "SizeMax", h_b)
        fm.setNumber(th, "DistMin", 0.0)
        fm.setNumber(th, "DistMax", B_WALL - A_ROD)
        fm.setAsBackgroundMesh(th)
        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
        tg = time.perf_counter()
        gmsh.model.mesh.generate(3)
        gmsh.model.mesh.setOrder(2)
        print(f"[imp8-i65] gmsh generate + setOrder(2) {time.perf_counter() - tg:.1f} s "
              f"(h_a {h_a / MM:g} mm, h_b {h_b / MM:g} mm)", flush=True)
    msh, ct, ft = _model_to_mesh(gmsh.model, comm, rank, gdim=3)
    if comm.rank == rank:
        gmsh.finalize()
    return msh, ct, ft


@pytest.fixture(scope="module")
def cyl():
    if _phase() is None:
        pytest.skip(f"{PHASE_ENV} unset: no gmsh in the coding window")
    comm = MPI.COMM_WORLD
    t0 = time.perf_counter()
    msh, ct, ft = _composite_cylinder(comm)
    t_build = time.perf_counter() - t0
    n = msh.topology.index_map(msh.topology.dim).size_global
    _say(f"\n[imp8-i65] composite cylinder: cells {n} (cap {CELL_CAP}); ranks {comm.size}; "
         f"t_build {t_build:.1f} s; phase {_phase()}")
    return msh, ct, ft, n, t_build


def test_w0_preconditions(cyl):
    msh, ct, ft, n, _ = cyl
    comm = msh.comm
    deg = msh.geometry.cmap.degree
    # Tag-1 facets: every geometry node (vertices and P2 edge nodes) at r = b.
    fdim = msh.topology.dim - 1
    msh.topology.create_connectivity(fdim, msh.topology.dim)
    facets = ft.find(ROBIN_TAG)
    gidx = dmesh.entities_to_geometry(msh, fdim, facets)
    x = msh.geometry.x
    rr = np.hypot(x[gidx.ravel(), 0], x[gidx.ravel(), 1]) if facets.size else np.zeros(0)
    dev_loc = float(np.max(np.abs(rr - B_WALL))) if rr.size else 0.0
    dev = comm.allreduce(dev_loc, op=MPI.MAX)
    n_fac = comm.allreduce(int(np.count_nonzero(
        facets < msh.topology.index_map(fdim).size_local)), op=MPI.SUM)
    mid = dmesh.compute_midpoints(msh, fdim, facets) if facets.size else np.zeros((0, 3))
    zmax_loc = float(np.max(np.abs(mid[:, 2]))) if facets.size else 0.0
    zmax = comm.allreduce(zmax_loc, op=MPI.MAX)
    # Wiring: f = x + 2y + 3z in P2 at the three read points.
    V = fem.functionspace(msh, ("Lagrange", DEGREE))
    f = fem.Function(V)
    f.interpolate(lambda X: X[0] + 2.0 * X[1] + 3.0 * X[2])
    pts = np.array([[r * MM, 0.0, z * MM] for r, z in READ_PTS], dtype=float)
    vals, ok = evaluate_vector_field_parallel(f, pts)
    werr = max(abs(float(np.real(np.ravel(vals[i])[0])) - (p[0] + 2 * p[1] + 3 * p[2]))
               if bool(ok[i]) else np.inf for i, p in enumerate(pts))
    dx = ufl.Measure("dx", domain=msh, subdomain_data=ct)
    v_core = _allsum(comm, fem.Constant(msh, 1.0) * dx(1))
    v_shell = _allsum(comm, fem.Constant(msh, 1.0) * dx(2))
    loc_tags = set(int(t) for t in np.unique(ct.values))
    all_tags = set().union(*comm.allgather(loc_tags))
    _say(f"[imp8-i65] W0 cells {n}; coordinate degree {deg}; cell tags {sorted(all_tags)}; "
         f"tag-{ROBIN_TAG} owned facets {n_fac}, max |r - b| over their geometry nodes "
         f"{dev:.3e} m, max |z| of their midpoints {zmax / MM:.4f} mm (L/2 = {LEN / 2 / MM:g}); "
         f"P2 wiring max err {werr:.3e}")
    _say(f"[imp8-i65] W0 V_core,h/(pi a^2 L) - 1 = {v_core / (np.pi * A_ROD**2 * LEN) - 1:.4e}; "
         f"V_shell,h/(pi (b^2-a^2) L) - 1 = "
         f"{v_shell / (np.pi * (B_WALL**2 - A_ROD**2) * LEN) - 1:.4e}")
    assert deg == 2, deg
    assert n <= CELL_CAP, n
    assert all_tags == {1, 2}, all_tags
    assert n_fac > 0
    assert dev <= 1e-6 * B_WALL, dev  # geometry nodes on r = b (gmsh round-off)
    assert zmax < 0.5 * LEN - 1e-6, zmax  # no end face carries the Robin tag
    assert werr <= 1e-12, werr


@pytest.fixture(scope="module")
def run(cyl):
    if _phase() != "solve":
        pytest.skip(f"{PHASE_ENV} != 'solve'")
    msh, ct, ft, n, t_build = cyl
    comm = msh.comm
    k = _dg0(msh, ct, ROD_K, K_GEL)
    rc = _dg0(msh, ct, ROD_RHO_C, RHO_C)
    q = _dg0(msh, ct, 0.0, Q_SHELL)
    nan_loc = int(sum(np.count_nonzero(np.isnan(np.real(g.x.array))) for g in (k, rc, q)))
    out = {"nan": comm.allreduce(nan_loc, op=MPI.SUM)}
    dx = ufl.Measure("dx", domain=msh, subdomain_data=ct)
    v_shell = _allsum(comm, fem.Constant(msh, 1.0) * dx(2))
    out["S0"] = abs(_allsum(comm, q * dx) / (Q_SHELL * v_shell) - 1.0)
    out["S0core"] = abs(_allsum(comm, q * dx(1)))
    dt = T_STAR / N_STEPS
    t1 = time.perf_counter()
    res = solve_transient_conduction(
        msh, ft, k=k, rho_c=rc, q=q, h=H_DIR, t_inf=T_INF, t_init=T_INF,
        robin_tag=ROBIN_TAG, dt=dt, n_steps=N_STEPS, theta=0.5, startup_half_steps=2,
        degree=DEGREE, snapshot_steps=(SNAP,))
    t_tr = time.perf_counter() - t1
    out["E"] = float(res.max_balance_residual)
    t2 = time.perf_counter()
    st = solve_steady_conduction(msh, ft, k=k, q=q, h=H_DIR, t_inf=T_INF,
                                 robin_tag=ROBIN_TAG, degree=DEGREE)
    t_st = time.perf_counter() - t2
    _say(f"[imp8-i65] W1 cells {n}; t_build {t_build:.1f} s; t_call transient {t_tr:.1f} s "
         f"(CN N={N_STEPS}, dt={dt:.6g} s, snapshot {SNAP} -> t = {SNAP * dt:.6g} s = t_c "
         f"{T_C:.6g}), n_factorizations {res.n_factorizations}; t_call steady {t_st:.1f} s; "
         f"NaN left {out['nan']}")
    th_c = _theta_pts(res.snapshots[SNAP], READ_PTS)
    th_s = _theta_pts(res.temperature, READ_PTS + WIT_PTS)
    th_ss = _theta_pts(st.temperature, READ_PTS)
    ex = {
        "t_c": [float(_theta(r * MM, T_C)) for r, _ in READ_PTS],
        "t*": [float(_theta(r * MM, T_STAR)) for r, _ in READ_PTS],
        "ss": [float(_theta_ss(r * MM)) for r, _ in READ_PTS],
    }
    fem_ = {"t_c": th_c, "t*": th_s[:3], "ss": th_ss}
    out["I"] = []
    for key in ("t_c", "t*", "ss"):
        for (r, z), v, e in zip(READ_PTS, fem_[key], ex[key]):
            rel = abs(v / e - 1.0)
            gated = (key in ("t*", "ss") and r in (RA, RM)) or (key == "t_c" and r == 0.0)
            tag = "(I0)" if key == "t_c" and r == 0.0 else "(I)"
            if gated and tag == "(I)":
                out["I"].append(rel)
            if gated and tag == "(I0)":
                out["I0"] = rel
            _say(f"[imp8-i65] {tag} ({r:g}, {z:g}) mm {key}: theta_FEM {v:.6e}, composite "
                 f"{e:.6e} K per W/m, rel {rel:.3e} "
                 f"({'gated' if gated else 'printed only'}, band {INTERIOR_BAND:.0e})")
    _say(f"[imp8-i65] z-invariance at t*: theta_FEM(a, 0) {th_s[0]:.6e}, (a, +L/4) "
         f"{th_s[3]:.6e}, (a, -L/4) {th_s[4]:.6e}")
    # Control: the k_r = K_GEL comparand at (0, t_c); same swap at the weekly's points.
    c0 = float(_theta(0.0, T_C, k_r=K_GEL))
    out["ctrl"] = abs(th_c[2] / c0 - 1.0)
    _say(f"[imp8-i65] control (k_r = K_GEL) at (0, t_c): comparand {c0:.6e}, miss "
         f"{out['ctrl']:.4e} (asserted > {INTERIOR_BAND:.0e}; predicted approx {PRED_CTRL}, "
         f"not asserted)")
    oth = [abs(th_c[0] / float(_theta(A_ROD, T_C, k_r=K_GEL)) - 1.0),
           abs(th_s[0] / float(_theta(A_ROD, T_STAR, k_r=K_GEL)) - 1.0),
           abs(th_ss[0] / float(_theta_ss(A_ROD)) - 1.0)]
    _say(f"[imp8-i65] control swap at the weekly's points (printed): (a, t_c) {oth[0]:.3e}, "
         f"(a, t*) {oth[1]:.3e}, (a, steady) {oth[2]:.3e} (predicted {PRED_CTRL_OTHERS})")
    # Robin offset: the steady wall temperature over theta_ss(a).
    ds = ufl.Measure("ds", domain=msh, subdomain_data=ft)
    area = _allsum(comm, fem.Constant(msh, 1.0) * ds(ROBIN_TAG))
    t_wall = _allsum(comm, st.temperature * ds(ROBIN_TAG)) / area
    _say(f"[imp8-i65] Robin offset: mean wall theta / theta_ss(a) = "
         f"{t_wall / float(_theta_ss(A_ROD)):.4e} (predicted {PRED_ROBIN:.1e}); steady "
         f"convected / source {st.convected_power / st.source_power - 1:.3e}; (E) {out['E']:.3e} "
         f"(band {E_BAND:.0e}); (S0) {out['S0']:.3e} (band {B_BAND:.0e}), core {out['S0core']:.3e}")
    return out


def test_w1_s0_source(run):
    assert run["nan"] == 0, run["nan"]
    assert run["S0"] <= B_BAND, run["S0"]
    assert run["S0core"] == 0.0, run["S0core"]


def test_w1_e_balance(run):
    assert run["E"] <= E_BAND, run["E"]


def test_w1_i_interior(run):
    assert len(run["I"]) == 4
    assert max(run["I"]) <= INTERIOR_BAND, run["I"]


def test_w1_i0_axis(run):
    assert run["I0"] <= INTERIOR_BAND, run["I0"]


def test_w1_control_k_gel(run):
    """Asserted at factor 1 (> ``INTERIOR_BAND``); the ≈ 0.133 is predicted, printed."""
    assert run["ctrl"] > INTERIOR_BAND, run["ctrl"]
