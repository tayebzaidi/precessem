"""`IMP-8` step 2 (§9 item 778, §10 I55): the finite line source in FEM on a
graded gmsh gel box, real build, degree 2, CN to t* = 900 s.

Brief: docs/planning/chunks/IMP-8.md, "§9 item 778".  A DG0 rod source (cells
whose midpoint has r ≤ ``ROD_RADIUS`` and |z| ≤ L/2, normalised to ∫q dV = P
= 1 W) in a Robin-cooled box [−0.10, 0.10]³ m, against step 0's infinite-medium
θ = q′/(4πk)·Φ.  Asserted: (E) max balance residual ≤ ``E_BAND``; (S0) ∫q = P
within ``B_BAND`` with the selected-cell count > 0; (I) θ at t* at four points
≤ ``INTERIOR_BAND``; (St) the steady box against Φ∞ ≤ ``INTERIOR_BAND`` —
PRE-REGISTERED RED, range [0.03, 0.30] (the walls are felt at steady).
Controls (predicted): (c-2L) L → 2L in the comparand misses (I) at (10, 15) mm
by > 0.1; (c-E) a steady solve with q(1 + 1e-3) balanced against the nominal P
misses by ≈ 1e-3 ≫ ``E_BAND``.  No band is set or moved here.

Item 800 (attempt 2, Review ruling 2026-11-03 10:30, docs/planning/chunks/
IMP-8.md): the windows are split by ``FEM_EM_IMP8_I55_PHASE`` = ``steady``
(W1s: (S0), (St), the price) / ``transient`` (W2: (E), (S0), (I), (c-2L)) /
``ce`` (W3: (c-E) on one full step of two BE half-steps against the nominal P,
plus the nominal one-step balance ≤ ``E_BAND``).  Each builds the mesh (W0c's
grading; its 485 291 cells asserted) and prints t_build, the solve call's wall
time and, for transients, ``n_factorizations``.

Item 803 (attempt 3, Review ruling 2026-11-03 18:00): (c-2L) is *asserted* at
the gate's own band — the L → 2L comparand misses θ_FEM at (10, 15) mm by
> ``INTERIOR_BAND`` — backed by the prior measurement of this comparison on
this fixture, 7.517e-02 (`20261103T182723Z_IMP-8-I55-W2-n8.log`, the (c-2L)
line); item 778's 0.1 is printed as predicted, not asserted.  (c-E) runs the
perturbed one-step solve alone and is *asserted* > ``E_BAND`` (factor 1); the
≈ 1e-3 miss is *predicted* and printed.  The nominal one-step balance test is
dropped as redundant: its two half-steps are the first two of W2's run at the
same Δt and ``startup_half_steps``, and W2's (E) line (1.898e-14 ≤ ``E_BAND``,
`20261103T182723Z_IMP-8-I55-W2-n8.log`) is the maximum over every half-step.
"""

from __future__ import annotations

import os
import time

import numpy as np
import pytest
import ufl
from dolfinx import fem, la
from dolfinx import mesh as dmesh
from mpi4py import MPI

from fem_em_solver.io.mesh import _model_to_mesh
from fem_em_solver.post.evaluation import evaluate_vector_field_parallel
from fem_em_solver.thermal import solve_steady_conduction, solve_transient_conduction
from tests.mesh.test_birdcage_implant_rod import ROD_RADIUS
from tests.validation.test_geo36_rod_f_clinical_heating import PERTURBED_SOURCE
from tests.validation.test_imp8_rod_line_source_closed_form import (
    MM,
    T_PRINT,
    T_STEADY,
    L,
    _phi,
    _phi_steady,
)
from tests.validation.test_mat10_sar_to_heat_sphere import B_BAND, ROBIN_TAG, _allsum, _say
from tests.validation.test_mat12_transient_sar_heating_sphere import E_BAND, RHO_C
from tests.validation.test_mat13_transient_perfusion_sphere import INTERIOR_BAND
from tests.validation.test_mat15_coil_driven_heating_phantom import H_GEL, K_GEL, N_MAIN

# Fixture parameters (rule (f)), not bands.
HALF_WIDTH = 0.10
AXIS_HALF = 0.5 * L + 2e-3
N_AXIS_PTS = 80
H_IN, D_IN0, H_MID, D_IN1 = 0.4e-3, 1.5e-3, 2.0e-3, 25e-3
H_OUT, D_OUT1 = 12e-3, 100e-3
CELL_CAP = 800_000
P = 1.0  # W, results per watt
DEGREE = 2
T_INF = 0.0
GATE_PTS = ((10.0, 0.0), (20.0, 0.0), (10.0, 15.0), (20.0, 15.0))  # (r, z) mm
PRINT_PTS = ((2.0, 0.0),)
C2L_PT = (10.0, 15.0)
C2L_PREDICTED = 0.1  # item 778's predicted factor, printed only (rule (e))
ST_LO, ST_HI = 0.03, 0.30  # PREDICTED RED range, not a band
#: Item 800 (attempt 2): W0c's measured count (`20261103T105748Z_IMP-8-I55-W0c.log`).
W0C_CELLS = 485_291
PHASE_ENV = "FEM_EM_IMP8_I55_PHASE"
PHASES = ("steady", "transient", "ce")


def _graded_box(comm, *, half_width=HALF_WIDTH, axis_half=AXIS_HALF):
    """gmsh OCC box [−hw, hw]³, one volume (tag 1), outer facets tag 1, graded on
    the z-axis segment |z| ≤ axis_half.  No embedded rod."""
    import gmsh

    rank = 0
    if comm.rank == rank:
        gmsh.initialize()
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("imp8_graded_box")
        occ = gmsh.model.occ
        hw = half_width
        box = occ.addBox(-hw, -hw, -hw, 2 * hw, 2 * hw, 2 * hw)
        occ.synchronize()
        surfs = [s for _, s in gmsh.model.getBoundary([(3, box)], oriented=False)]
        gmsh.model.addPhysicalGroup(3, [box], tag=1)
        gmsh.model.addPhysicalGroup(2, surfs, tag=ROBIN_TAG)
        pts = [
            occ.addPoint(0.0, 0.0, float(z)) for z in np.linspace(-axis_half, axis_half, N_AXIS_PTS)
        ]
        occ.synchronize()
        fm = gmsh.model.mesh.field
        fd = fm.add("Distance")
        fm.setNumbers(fd, "PointsList", pts)
        t_in = fm.add("Threshold")
        fm.setNumber(t_in, "InField", fd)
        fm.setNumber(t_in, "SizeMin", H_IN)
        fm.setNumber(t_in, "SizeMax", H_MID)
        fm.setNumber(t_in, "DistMin", D_IN0)
        fm.setNumber(t_in, "DistMax", D_IN1)
        # Beyond DistMax a Threshold returns SizeMax (2 mm) unless stopped, and the
        # Min would then cap the whole box at 2 mm: 4 101 108 cells measured in
        # 20261103T105626Z_IMP-8-I55-W0b-probe.log.  Stop it so the outer one grades.
        fm.setNumber(t_in, "StopAtDistMax", 1)
        t_out = fm.add("Threshold")
        fm.setNumber(t_out, "InField", fd)
        fm.setNumber(t_out, "SizeMin", H_MID)
        fm.setNumber(t_out, "SizeMax", H_OUT)
        fm.setNumber(t_out, "DistMin", D_IN1)
        fm.setNumber(t_out, "DistMax", D_OUT1)
        fmin = fm.add("Min")
        fm.setNumbers(fmin, "FieldsList", [t_in, t_out])
        fm.setAsBackgroundMesh(fmin)
        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
        algo = os.environ.get("IMP8_ALGO3D")  # cost probe only; unset = gmsh default
        if algo:
            gmsh.option.setNumber("Mesh.Algorithm3D", int(algo))
            gmsh.option.setNumber("General.NumThreads", 2)
        tg = time.perf_counter()
        gmsh.model.mesh.generate(3)
        _, etags, _ = gmsh.model.mesh.getElements(3)
        ntet = sum(len(t) for t in etags)
        print(
            f"[imp8-fem] gmsh: algo3d {algo or 'default'}; tets {ntet}; generate "
            f"{time.perf_counter() - tg:.1f} s",
            flush=True,
        )
    msh, ct, ft = _model_to_mesh(gmsh.model, comm, rank, gdim=3)
    if comm.rank == rank:
        gmsh.finalize()
    return msh, ct, ft


def _cells_global(msh):
    return msh.topology.index_map(msh.topology.dim).size_global


def _rod_source(msh, power=P, length=L, radius=ROD_RADIUS):
    """DG0 indicator of cells with midpoint r ≤ radius, |z| ≤ length/2,
    normalised so the owned-volume sum gives ∫q dV = power exactly."""
    comm = msh.comm
    tdim = msh.topology.dim
    n_loc = msh.topology.index_map(tdim).size_local
    n_all = n_loc + msh.topology.index_map(tdim).num_ghosts
    mid = dmesh.compute_midpoints(msh, tdim, np.arange(n_all, dtype=np.int32))
    sel = (np.hypot(mid[:, 0], mid[:, 1]) <= radius) & (np.abs(mid[:, 2]) <= 0.5 * length)
    W = fem.functionspace(msh, ("DG", 0))
    w = ufl.TestFunction(W)
    one = fem.Constant(msh, np.asarray(1.0, dtype=fem.Function(W).x.array.dtype))
    bvol = fem.assemble_vector(fem.form(one * ufl.conj(w) * ufl.dx))
    bvol.scatter_reverse(la.InsertMode.add)
    bvol.scatter_forward()
    dofs = np.asarray(W.dofmap.list).reshape(-1)[:n_all]  # cell -> its DG0 dof
    vol = np.real(bvol.array[dofs])
    own = sel[:n_loc]  # owned cells own their DG0 dof
    v_sel = comm.allreduce(float(np.sum(vol[:n_loc][own])), op=MPI.SUM)
    n_sel = comm.allreduce(int(np.count_nonzero(own)), op=MPI.SUM)
    q = fem.Function(W)
    q.x.array[:] = 0.0
    if v_sel > 0:
        q.x.array[dofs[sel]] = power / v_sel
    q.x.scatter_forward()
    return q, n_sel, v_sel


def _theta_pts(T, pts_mm):
    pts = np.array([[r * MM, 0.0, z * MM] for r, z in pts_mm], dtype=float)
    vals, ok = evaluate_vector_field_parallel(T, pts)
    return [
        float(np.real(np.ravel(vals[i])[0])) - T_INF if bool(ok[i]) else float("nan")
        for i in range(len(pts))
    ]


def _h_near(msh, pts_mm):
    """Max edge length of the cells whose midpoint lies within 3 mm of each point
    (allreduced max)."""
    tdim = msh.topology.dim
    n_loc = msh.topology.index_map(tdim).size_local
    cells = np.arange(n_loc, dtype=np.int32)
    mid = dmesh.compute_midpoints(msh, tdim, cells)
    hs = msh.h(tdim, cells)
    out = []
    for r, z in pts_mm:
        d = np.linalg.norm(mid - np.array([r * MM, 0.0, z * MM]), axis=1)
        loc = float(np.max(hs[d <= 3e-3])) if np.any(d <= 3e-3) else 0.0
        out.append(msh.comm.allreduce(loc, op=MPI.MAX))
    return out


def _phase():
    raw = os.environ.get(PHASE_ENV, "").strip()
    if raw == "":
        return None
    if raw in PHASES:
        return raw
    raise ValueError(f"{PHASE_ENV} must be unset or one of {PHASES}, got {raw!r}")


def _need(phase):
    if _phase() != phase:
        pytest.skip(f"{PHASE_ENV} != {phase!r}: item 800's windows are selected by environment")


def test_k0_phase_knob(monkeypatch):
    monkeypatch.delenv(PHASE_ENV, raising=False)
    assert _phase() is None
    for ph in PHASES:
        monkeypatch.setenv(PHASE_ENV, ph)
        assert _phase() == ph
    for bad in ("1", "all", "St"):
        monkeypatch.setenv(PHASE_ENV, bad)
        with pytest.raises(ValueError):
            _phase()


@pytest.fixture(scope="module")
def box():
    comm = MPI.COMM_WORLD
    t0 = time.perf_counter()
    msh, ct, ft = _graded_box(comm)
    n = _cells_global(msh)
    t_build = time.perf_counter() - t0
    _say(
        f"\n[imp8-fem] graded box: cells {n} (W0c {W0C_CELLS}; cap {CELL_CAP}); "
        f"ranks {comm.size}; t_build {t_build:.1f} s; phase {_phase()}"
    )
    return msh, ct, ft, n


def test_w0_build(box):
    if _phase() is None:
        pytest.skip(f"{PHASE_ENV} unset: no gmsh in the coding window")
    msh, ct, ft, n = box
    q, n_sel, v_sel = _rod_source(msh)
    a_eff = np.sqrt(v_sel / (np.pi * L))
    hs = _h_near(msh, GATE_PTS + PRINT_PTS)
    _say(
        f"[imp8-fem] W0 selected cells {n_sel}; selected volume {v_sel:.6e} m^3; effective "
        f"radius {a_eff / MM:.4f} mm (ROD_RADIUS {ROD_RADIUS / MM:g} mm)"
    )
    for (r, z), h in zip(GATE_PTS + PRINT_PTS, hs):
        _say(f"[imp8-fem] W0 h within 3 mm of ({r:g}, {z:g}) mm: {h / MM:.3f} mm")
    assert n <= CELL_CAP, n
    # Item 800: a different count than W0c's means the grading moved; stop.
    assert n == W0C_CELLS, n
    assert n_sel > 0


def _source_and_s0(msh):
    q, n_sel, v_sel = _rod_source(msh)
    s0 = float(np.real(_allsum(msh.comm, q * ufl.dx)))
    return q, {"n_sel": n_sel, "S0": abs(s0 - P) / P, "a_eff": np.sqrt(v_sel / (np.pi * L))}


QP = P / L / (4.0 * np.pi * K_GEL)  # q′/(4πk), K per unit Φ


# --- W1s: the build and the (St) steady solve (the price). -------------------


@pytest.fixture(scope="module")
def steady_run(box):
    _need("steady")
    msh, ct, ft, n = box
    q, out = _source_and_s0(msh)
    t1 = time.perf_counter()
    st = solve_steady_conduction(
        msh, ft, k=K_GEL, q=q, h=H_GEL, t_inf=T_INF, robin_tag=ROBIN_TAG, degree=DEGREE
    )
    out["t_call"] = time.perf_counter() - t1
    pts = GATE_PTS + PRINT_PTS
    th_s = _theta_pts(st.temperature, pts)
    out["St"] = []
    for (r, z), v in zip(pts, th_s):
        ex = QP * _phi_steady(r * MM, z * MM)
        ex_t = QP * _phi(r * MM, z * MM, T_STEADY)
        rel = abs(v / ex - 1.0)
        if (r, z) in GATE_PTS:
            out["St"].append(rel)
        _say(
            f"[imp8-fem] (St) ({r:g}, {z:g}) mm steady box: theta {v:.6e} K/W, Phi_inf form "
            f"{ex:.6e} (t=1e9 Phi {ex_t:.6e}), offset {v - ex:+.4e} K/W, rel {rel:.3e} "
            f"(band {INTERIOR_BAND:.0e}; PREDICTED RED range [{ST_LO}, {ST_HI}])"
        )
    gauss = abs(float(st.convected_power) - P) / P
    _say(
        f"[imp8-fem] W1s t_call (steady) {out['t_call']:.1f} s; convected vs P {gauss:.3e} "
        f"(printed); (S0) {out['S0']:.3e} (band {B_BAND:.0e}), selected cells {out['n_sel']}, "
        f"effective radius {out['a_eff'] / MM:.4f} mm; W2 projection = the box line's "
        f"t_build + 3 x {out['t_call']:.1f} s (-n 2 if <= 480 s)"
    )
    return out


def test_w1s_s0_source_power(steady_run):
    assert steady_run["n_sel"] > 0
    assert steady_run["S0"] <= B_BAND, steady_run["S0"]


def test_w1s_st_steady_against_infinite_medium(steady_run):
    """PREDICTED RED: (St) range [0.03, 0.30] (box Robin walls at 0.10 m)."""
    st = steady_run["St"]
    _say(
        f"[imp8-fem] (St) max rel {max(st):.4e}; inside PREDICTED RED range "
        f"[{ST_LO}, {ST_HI}]: {ST_LO <= max(st) <= ST_HI and min(st) >= 0.0}"
    )
    assert max(st) <= INTERIOR_BAND, st


# --- W2: the transient ((E), (S0), (I), (c-2L)). ------------------------------


@pytest.fixture(scope="module")
def transient_run(box):
    _need("transient")
    msh, ct, ft, n = box
    q, out = _source_and_s0(msh)
    dt = T_PRINT / N_MAIN
    t1 = time.perf_counter()
    res = solve_transient_conduction(
        msh,
        ft,
        k=K_GEL,
        rho_c=RHO_C,
        q=q,
        h=H_GEL,
        t_inf=T_INF,
        t_init=T_INF,
        robin_tag=ROBIN_TAG,
        dt=dt,
        n_steps=N_MAIN,
        theta=0.5,
        startup_half_steps=2,
        degree=DEGREE,
    )
    out["t_call"] = time.perf_counter() - t1
    out["E"] = float(res.max_balance_residual)
    _say(
        f"[imp8-fem] W2 CN N={N_MAIN} dt={dt:g} s: t_call {out['t_call']:.1f} s; "
        f"n_factorizations {res.n_factorizations}"
    )
    pts = GATE_PTS + PRINT_PTS
    th = _theta_pts(res.temperature, pts)
    out["I"], out["c2L"] = [], None
    for (r, z), v in zip(pts, th):
        ex = QP * _phi(r * MM, z * MM, T_PRINT)
        rel = abs(v / ex - 1.0)
        a2 = (ROD_RADIUS / (r * MM)) ** 2
        gated = (r, z) in GATE_PTS
        if gated:
            out["I"].append(rel)
        _say(
            f"[imp8-fem] (I) ({r:g}, {z:g}) mm t*={T_PRINT:g} s: theta_FEM {v:.6e} K/W, "
            f"closed form {ex:.6e}, rel {rel:.3e} ({'gated' if gated else 'printed only'}, "
            f"band {INTERIOR_BAND:.0e}); (a/r)^2 {a2:.3e}"
        )
        if (r, z) == C2L_PT:
            ex2 = P / (2 * L) / (4.0 * np.pi * K_GEL) * _phi(r * MM, z * MM, T_PRINT, length=2 * L)
            out["c2L"] = abs(v / ex2 - 1.0)
    _say(
        f"[imp8-fem] (c-2L) L->2L comparand at {C2L_PT} mm misses by {out['c2L']:.4e} "
        f"(asserted > INTERIOR_BAND {INTERIOR_BAND:.0e}; predicted {C2L_PREDICTED} "
        f"(item 778), not asserted); (E) {out['E']:.3e} (band {E_BAND:.0e}); (S0) "
        f"{out['S0']:.3e}"
    )
    return out


def test_w2_e_balance(transient_run):
    assert transient_run["E"] <= E_BAND, transient_run["E"]


def test_w2_s0_source_power(transient_run):
    assert transient_run["n_sel"] > 0
    assert transient_run["S0"] <= B_BAND, transient_run["S0"]


def test_w2_i_interior_against_line_source(transient_run):
    assert max(transient_run["I"]) <= INTERIOR_BAND, transient_run["I"]


def test_w2_control_c2l(transient_run):
    """Asserted at (I)'s own band; backed by 7.517e-02 in
    `20261103T182723Z_IMP-8-I55-W2-n8.log` (item 803 ruling)."""
    assert transient_run["c2L"] > INTERIOR_BAND, transient_run["c2L"]


# --- W3: (c-E) on one full step (two backward-Euler half-steps). --------------


@pytest.fixture(scope="module")
def ce_run(box):
    _need("ce")
    msh, ct, ft, n = box
    q, out = _source_and_s0(msh)
    dt = T_PRINT / N_MAIN
    qe = fem.Function(q.function_space)
    qe.x.array[:] = q.x.array * (1.0 + PERTURBED_SOURCE)
    kw = dict(
        k=K_GEL,
        rho_c=RHO_C,
        h=H_GEL,
        t_inf=T_INF,
        t_init=T_INF,
        robin_tag=ROBIN_TAG,
        dt=dt,
        n_steps=1,
        theta=0.5,
        startup_half_steps=2,
        degree=DEGREE,
    )
    t1 = time.perf_counter()
    pert = solve_transient_conduction(msh, ft, q=qe, **kw)
    out["t_call"] = time.perf_counter() - t1
    # The perturbed run's residual r_n is against its own ∫Qe; its storage + loss
    # is ∫Qe (1 ± r_n), so against the nominal P it misses by
    # |(∫Qe/P)(1 ± r_n) − 1| — bounded below by |∫Qe/P − 1| − (∫Qe/P)|r_1| on the
    # first half-step.
    ratio = float(pert.source_power) / P
    r_1 = abs(float(pert.balance_residuals[0]))
    out["cE_lo"] = abs(ratio - 1.0) - ratio * r_1
    out["cE_hi"] = abs(ratio - 1.0) + ratio * r_1
    out["n_half"] = len(pert.balance_residuals)
    _say(
        f"[imp8-fem] W3 one step (2 BE half-steps), perturbed only: t_call {out['t_call']:.1f} s; "
        f"n_factorizations {pert.n_factorizations}; residual count {out['n_half']}\n"
        f"[imp8-fem] (c-E) perturbed: own first-half-step residual {r_1:.3e} (own max "
        f"{float(pert.max_balance_residual):.3e}); int Qe / P - 1 = {ratio - 1.0:.6e}; first "
        f"half-step miss against nominal P in [{out['cE_lo']:.6e}, {out['cE_hi']:.6e}] "
        f"(asserted > E_BAND {E_BAND:.0e}; predicted approx {PERTURBED_SOURCE:.0e}, not "
        f"asserted); nominal reference: W2's (E) line, 20261103T182723Z_IMP-8-I55-W2-n8.log; "
        f"(S0) {out['S0']:.3e}"
    )
    return out


def test_w3_control_ce(ce_run):
    """Asserted at factor 1 (> ``E_BAND``); the ≈ 1e-3 miss is predicted, printed."""
    assert ce_run["S0"] <= B_BAND, ce_run["S0"]
    assert ce_run["cE_lo"] > E_BAND, ce_run["cE_lo"]
