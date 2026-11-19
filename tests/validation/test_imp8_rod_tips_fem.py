"""`IMP-8` step 3 (§9 item 796, §10 I56): the two-tip twin in FEM, on item 803's
graded gel box (``_graded_box`` via its ``box`` fixture, imported unchanged),
real build, degree 2, CN to t* = 900 s.

Brief: docs/planning/chunks/IMP-8.md, "§9 item 796" (with the 2026-11-03 10:30
ruling 6 and the 18:00 ruling 5: it depends on item 803 and runs at ``-n 8``).
The source is two DG0 indicators — the cells whose midpoint lies within
``ROD_RADIUS`` of (0, 0, ±L/2) — each normalised to exactly P/2.  Asserted:
(E) max balance residual ≤ ``E_BAND``; (S0) ∫q = P, and each ball's ∫q = P/2,
within ``B_BAND`` with selected-cell counts > 0; (I) θ_FEM at t* against
q′/(4πk)·Φ₂ (step 1's ``_phi_tips``) at item 778's four points within
``INTERIOR_BAND``.  No (St).

Printed as a record, no band: the FEM Φ₂/Φ at (2, 15) and (2, 0) mm beside
step 1's 3.179482 / 0.1714577.  At (2, 0) mm the line's θ_FEM is item 803's
1.772315e+01 K/W (`20261104T030134Z_IMP-8-I55-W2-item803.log:145`).  Item 803
did not evaluate (2, 15) mm, so the ratio there is printed against the line's
closed form q′/(4πk)·Φ (disclosed in the line itself), not a second solve.

Control (c-one), asserted at the gate's own band, factor 1 (rule (e)): the
comparand with one ball dropped (P/2 in total; the z = −L/2 ball kept) misses
θ_FEM by > ``INTERIOR_BAND`` at every gate point — arithmetic on the same
solve.  The brief's "≈ 0.5" has no log behind it, so it is a *predicted*
factor, printed beside the measured miss and not asserted.
"""

from __future__ import annotations

import time

import numpy as np
import pytest
import ufl
from dolfinx import fem, la
from dolfinx import mesh as dmesh
from mpi4py import MPI

from fem_em_solver.thermal import solve_transient_conduction
from tests.mesh.test_birdcage_implant_rod import ROD_RADIUS
from tests.validation.test_imp8_rod_line_source_closed_form import (
    MM,
    T_PRINT,
    L,
    _phi,
    _phi_tips,
)
from tests.validation.test_imp8_rod_line_source_fem import (  # noqa: F401  (box: fixture)
    DEGREE,
    GATE_PTS,
    QP,
    T_INF,
    W0C_CELLS,
    P,
    _need,
    _theta_pts,
    box,
)
from tests.validation.test_mat10_sar_to_heat_sphere import B_BAND, ROBIN_TAG, _allsum, _say
from tests.validation.test_mat12_transient_sar_heating_sphere import E_BAND, RHO_C
from tests.validation.test_mat13_transient_perfusion_sphere import INTERIOR_BAND
from tests.validation.test_mat15_coil_driven_heating_phantom import H_GEL, K_GEL, N_MAIN

RATIO_PTS = ((2.0, 15.0), (2.0, 0.0))  # (r, z) mm, printed record
STEP1_RATIO = {(2.0, 15.0): 3.179482, (2.0, 0.0): 0.1714577}  # step 1's Φ₂/Φ
#: Item 803's line-source θ_FEM at (2, 0) mm, K/W
#: (`20261104T030134Z_IMP-8-I55-W2-item803.log:145`).
LINE_FEM_2_0 = 1.772315e01
C_ONE_PREDICTED = 0.5  # the brief's factor, printed only (rule (e))
SIGNS = (1.0, -1.0)


def _tips_source(msh, power=P, radius=ROD_RADIUS, length=L):
    """Two DG0 indicators: cells with midpoint within ``radius`` of
    (0, 0, s·length/2), s = ±1, each normalised so its owned-volume sum gives
    ∫q dV = power/2 exactly.  Returns (q, per-ball [(n_sel, v_sel, q_ball)])."""
    comm = msh.comm
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
    vol = np.real(bvol.array[dofs])
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


@pytest.fixture(scope="module")
def tips_run(box):
    _need("transient")
    msh, ct, ft, n = box
    assert n == W0C_CELLS, n  # item 803's grading, unchanged
    q, balls = _tips_source(msh)
    out = {"S0": abs(float(np.real(_allsum(msh.comm, q * ufl.dx)))) / P - 1.0}
    out["S0"] = abs(out["S0"])
    out["balls"] = []
    for s, (n_sel, v_sel, qb) in zip(SIGNS, balls):
        pb = float(np.real(_allsum(msh.comm, qb * ufl.dx)))
        a_eq = (3.0 * v_sel / (4.0 * np.pi)) ** (1.0 / 3.0)
        out["balls"].append((n_sel, abs(pb / (0.5 * P) - 1.0)))
        _say(
            f"[imp8-tips] ball z={s * 0.5 * L / MM:+g} mm: cells {n_sel}; volume "
            f"{v_sel:.6e} m^3; equivalent radius {a_eq / MM:.4f} mm (ROD_RADIUS "
            f"{ROD_RADIUS / MM:g} mm); int q / (P/2) - 1 = {pb / (0.5 * P) - 1.0:.3e}"
        )
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
        f"[imp8-tips] W1 CN N={N_MAIN} dt={dt:g} s: t_call {out['t_call']:.1f} s; "
        f"n_factorizations {res.n_factorizations}; ranks {msh.comm.size}"
    )
    pts = GATE_PTS + RATIO_PTS
    th = _theta_pts(res.temperature, pts)
    out["I"], out["c_one"] = [], []
    for (r, z), v in zip(pts, th):
        ex = QP * _phi_tips(r * MM, z * MM, T_PRINT)
        rel = abs(v / ex - 1.0)
        R_plus = np.hypot(r * MM, z * MM - 0.5 * L)
        a2 = (ROD_RADIUS / R_plus) ** 2
        if (r, z) in GATE_PTS:
            out["I"].append(rel)
            ex1 = QP * _phi_tips(r * MM, z * MM, T_PRINT, signs=(-1.0,))
            m_fem = abs(v / ex1 - 1.0)  # (I)'s convention, |θ_FEM / comparand − 1|
            m_cmp = abs(ex1 / v - 1.0)  # the comparand's shortfall relative to θ_FEM
            out["c_one"].append(m_fem)
            _say(
                f"[imp8-tips] (c-one) ({r:g}, {z:g}) mm: one-ball comparand {ex1:.6e} K/W; "
                f"miss |theta_FEM/cmp - 1| {m_fem:.4e} (asserted > INTERIOR_BAND "
                f"{INTERIOR_BAND:.0e}); |cmp/theta_FEM - 1| {m_cmp:.4e} (predicted approx "
                f"{C_ONE_PREDICTED}, not asserted)"
            )
        _say(
            f"[imp8-tips] (I) ({r:g}, {z:g}) mm t*={T_PRINT:g} s: theta_FEM {v:.6e} K/W, "
            f"closed form q'/(4 pi k) Phi_2 {ex:.6e}, rel {rel:.3e} "
            f"({'gated' if (r, z) in GATE_PTS else 'printed only'}, band "
            f"{INTERIOR_BAND:.0e}); (a/R+)^2 {a2:.3e}"
        )
        if (r, z) in RATIO_PTS:
            line_cf = QP * _phi(r * MM, z * MM, T_PRINT)
            if (r, z) == (2.0, 0.0):
                ratio = v / LINE_FEM_2_0
                src = f"item 803's line theta_FEM {LINE_FEM_2_0:.6e} (W2 log :145)"
            else:
                ratio = v / line_cf
                src = (
                    f"the line's closed form {line_cf:.6e} (item 803 did not evaluate "
                    f"this point; not a FEM/FEM ratio)"
                )
            _say(
                f"[imp8-tips] record (no band): Phi_2/Phi at ({r:g}, {z:g}) mm = {ratio:.6f} "
                f"(this run's theta over {src}); FEM tips / line closed form "
                f"{v / line_cf:.6f}; step 1 closed form {STEP1_RATIO[(r, z)]:.7g}"
            )
    _say(
        f"[imp8-tips] (E) {out['E']:.3e} (band {E_BAND:.0e}); (S0) {out['S0']:.3e} "
        f"(band {B_BAND:.0e}); (I) max {max(out['I']):.3e}; (c-one) min "
        f"{min(out['c_one']):.4e}"
    )
    return out


def test_w1_e_balance(tips_run):
    assert tips_run["E"] <= E_BAND, tips_run["E"]


def test_w1_s0_source_power(tips_run):
    assert tips_run["S0"] <= B_BAND, tips_run["S0"]
    for n_sel, miss in tips_run["balls"]:
        assert n_sel > 0
        assert miss <= B_BAND, miss


def test_w1_i_interior_against_tips(tips_run):
    assert max(tips_run["I"]) <= INTERIOR_BAND, tips_run["I"]


def test_w1_control_c_one(tips_run):
    """Asserted at (I)'s own band, factor 1; the ≈ 0.5 is predicted, printed."""
    assert min(tips_run["c_one"]) > INTERIOR_BAND, tips_run["c_one"]
