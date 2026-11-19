"""`IMP-8` step 11 (§9 item 846, §10 I66): I56's rod as a conducting body.

Brief: docs/planning/chunks/IMP-8.md, "§9 item 846", under the "Review ruling
2026-11-05 18:00" paragraph (which binds over it).  Item 803's graded gel box
(``box`` fixture, imported unchanged), real build, degree 2, CN to t* = 900 s,
I56's two-tip source (``_tips_source``, the same call, byte-identical).  The rod
cells (midpoint r ≤ ``ROD_RADIUS`` ∧ |z| ≤ L/2) minus the two ball sets (a cell
in both keeps the gel's properties and carries the source — the 18:00 ruling)
are given (``ROD_K``, ``ROD_RHO_C``, both nominal) through DG0 ``k`` / ``rho_c``.

Arm knob ``FEM_EM_IMP8_I66_ARM`` = ``ti`` (default) / ``gel`` / ``hole``.  Arm
``hole`` (§9 item 851, §10 I69, ``IMP-8`` step 13) gives the rod cells
k = ε·``K_GEL`` and ρc = ε·``RHO_C``, ε = ``EPS_HOLE`` = 1e-6 — I64's adiabatic-hole
core in the FEM's terms; the source cells keep gel properties.  Its negative
control: θ(2, 0) mm moves from I56's gel θ by more than 1e-9 relative; it prints
θ beside item 846's Ti arm and I56's gel θ (the bracket, a record).
Arm ``gel`` uses the same DG0
objects with gel values everywhere: the reproduction control (I56's six θ to
one unit of the sixth significant digit).  Arm ``ti``'s negative control: θ at
(2, 0) mm moves from I56's by more than 1e-9 relative (a DG0 path that ignored
its coefficients would return I56's θ).

Asserted (bands imported): (E) ``max_balance_residual`` ≤ ``E_BAND``; (S0)
∫q = P and each ball's ∫q = P/2 within ``B_BAND``, counts > 0; preconditions
cells = ``W0C_CELLS``, rod and overlap counts > 0 (allreduced).  Everything else
is a 🟡 record, printed, no band: no closed form has the finite rod.
"""

from __future__ import annotations

import math
import os
import time

import numpy as np
import pytest
import ufl
from dolfinx import fem, la
from dolfinx import mesh as dmesh
from mpi4py import MPI

from fem_em_solver.thermal import solve_transient_conduction
from tests.mesh.test_birdcage_implant_rod import ROD_RADIUS
from tests.validation.test_imp8_composite_cylinder_closed_form import ROD_K  # nominal
from tests.validation.test_imp8_composite_cylinder_closed_form import ROD_RHO_C  # nominal
from tests.validation.test_imp8_rod_line_source_closed_form import MM, T_PRINT, L, _phi_tips
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
from tests.validation.test_imp8_rod_tips_fem import RATIO_PTS, SIGNS, _tips_source
from tests.validation.test_mat10_sar_to_heat_sphere import B_BAND, ROBIN_TAG, _allsum, _say
from tests.validation.test_mat12_transient_sar_heating_sphere import E_BAND, RHO_C
from tests.validation.test_mat13_transient_perfusion_sphere import INTERIOR_BAND  # noqa: F401
from tests.validation.test_mat15_coil_driven_heating_phantom import H_GEL, K_GEL, N_MAIN

ARM_ENV = "FEM_EM_IMP8_I66_ARM"
ARMS = ("ti", "gel", "hole")
#: I56's printed θ, K/W (`20261104T043245Z_IMP-8-I56-W1-item796.log:106–115`),
#: in ``GATE_PTS + RATIO_PTS`` order.
I56_THETA = {
    (10.0, 0.0): 1.932331e00,
    (20.0, 0.0): 6.373868e-01,
    (10.0, 15.0): 3.645521e00,
    (20.0, 15.0): 7.543411e-01,
    (2.0, 15.0): 3.007493e01,
    (2.0, 0.0): 3.044228e00,
}
CONTROL_FLOOR = 1e-9  # the round-off floor (I56's (E) is 1.754e-14)
#: The 18:00 review's axial-conductance prediction (unlogged arithmetic), printed only.
MOVE_PREDICTED = 0.09
#: Item 851: the rod cells' k and ρc scale on arm ``hole`` (I64's adiabatic-hole core).
EPS_HOLE = 1e-6
#: Item 846's Ti-arm θ, K/W (`20261106T013519Z_IMP-8-item846-W1.log:115–120`), printed
#: beside arm ``hole`` only.
TI_THETA = {
    (10.0, 0.0): 2.365053e00,
    (20.0, 0.0): 7.194612e-01,
    (10.0, 15.0): 3.506584e00,
    (20.0, 15.0): 7.502096e-01,
    (2.0, 15.0): 2.525271e01,
    (2.0, 0.0): 4.735218e00,
}


def _arm():
    raw = os.environ.get(ARM_ENV, "").strip()
    if raw == "":
        return "ti"
    if raw in ARMS:
        return raw
    raise ValueError(f"{ARM_ENV} must be unset or one of {ARMS}, got {raw!r}")


def test_k0_arm_knob(monkeypatch):
    monkeypatch.delenv(ARM_ENV, raising=False)
    assert _arm() == "ti"
    for a in ARMS:
        monkeypatch.setenv(ARM_ENV, a)
        assert _arm() == a
    for bad in ("Ti", "pec", "1"):
        monkeypatch.setenv(ARM_ENV, bad)
        with pytest.raises(ValueError):
            _arm()


def _rod_cells(msh, radius=ROD_RADIUS, length=L):
    """Owned + ghost cell indices: rod (r ≤ radius ∧ |z| ≤ length/2) minus the
    two balls (midpoint within radius of (0, 0, ±length/2)), and the overlap.
    Returns (rod_idx, overlap_idx, info) with allreduced owned counts/volumes."""
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
    W = fem.functionspace(msh, ("DG", 0))
    w = ufl.TestFunction(W)
    one = fem.Constant(msh, np.asarray(1.0, dtype=fem.Function(W).x.array.dtype))
    bvol = fem.assemble_vector(fem.form(one * w * ufl.dx))
    bvol.scatter_reverse(la.InsertMode.add)
    bvol.scatter_forward()
    dofs = np.asarray(W.dofmap.list).reshape(-1)[:n_all]
    vol = np.real(bvol.array[dofs])
    info = {}
    for name, sel in (("rod", rod), ("overlap", over), ("cyl", cyl)):
        own = sel[:n_loc]
        info[f"n_{name}"] = comm.allreduce(int(np.count_nonzero(own)), op=MPI.SUM)
        info[f"v_{name}"] = comm.allreduce(float(np.sum(vol[:n_loc][own])), op=MPI.SUM)
    return np.flatnonzero(rod).astype(np.int32), np.flatnonzero(over).astype(np.int32), info


def _coefficients(msh, rod_idx, arm):
    W = fem.functionspace(msh, ("DG", 0))
    n_all = (
        msh.topology.index_map(msh.topology.dim).size_local
        + msh.topology.index_map(msh.topology.dim).num_ghosts
    )
    dofs = np.asarray(W.dofmap.list).reshape(-1)[:n_all]
    k = fem.Function(W)
    rho_c = fem.Function(W)
    k.x.array[:] = K_GEL
    rho_c.x.array[:] = RHO_C
    if arm == "ti":
        k.x.array[dofs[rod_idx]] = ROD_K
        rho_c.x.array[dofs[rod_idx]] = ROD_RHO_C
    elif arm == "hole":
        k.x.array[dofs[rod_idx]] = EPS_HOLE * K_GEL
        rho_c.x.array[dofs[rod_idx]] = EPS_HOLE * RHO_C
    k.x.scatter_forward()
    rho_c.x.scatter_forward()
    bad = int(
        np.count_nonzero(~np.isfinite(k.x.array)) + np.count_nonzero(~np.isfinite(rho_c.x.array))
    )
    return k, rho_c, msh.comm.allreduce(bad, op=MPI.SUM)


@pytest.fixture(scope="module")
def body_run(box):
    _need("transient")
    arm = _arm()
    msh, ct, ft, n = box
    rod_idx, over_idx, info = _rod_cells(msh)
    v_ideal = np.pi * ROD_RADIUS**2 * L
    _say(
        f"[imp8-i66] arm {arm}; ROD_K {ROD_K:g} W/(m K), ROD_RHO_C {ROD_RHO_C:.6g} J/(m3 K) "
        f"(nominal); K_GEL {K_GEL:g}, RHO_C {RHO_C:.6g}"
    )
    _say(
        f"[imp8-i66] cylinder cells {info['n_cyl']} volume {info['v_cyl']:.6e} m^3; rod cells "
        f"(cylinder minus balls) {info['n_rod']} volume {info['v_rod']:.6e} m^3 "
        f"(/(pi a^2 L) {info['v_rod'] / v_ideal:.6f}); overlap cells {info['n_overlap']} volume "
        f"{info['v_overlap']:.6e} m^3 (gel properties, carry the source)"
    )
    k, rho_c, n_bad = _coefficients(msh, rod_idx, arm)
    q, balls = _tips_source(msh)
    out = {"arm": arm, "n": n, "n_bad": n_bad, **info}
    out["S0"] = abs(float(np.real(_allsum(msh.comm, q * ufl.dx))) / P - 1.0)
    out["balls"] = []
    for s, (n_sel, v_sel, qb) in zip(SIGNS, balls):
        pb = float(np.real(_allsum(msh.comm, qb * ufl.dx)))
        out["balls"].append((n_sel, abs(pb / (0.5 * P) - 1.0)))
        _say(
            f"[imp8-i66] ball z={s * 0.5 * L / MM:+g} mm: cells {n_sel}; volume {v_sel:.6e} m^3; "
            f"int q / (P/2) - 1 = {pb / (0.5 * P) - 1.0:.3e}"
        )
    dt = T_PRINT / N_MAIN
    t1 = time.perf_counter()
    res = solve_transient_conduction(
        msh,
        ft,
        k=k,
        rho_c=rho_c,
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
        f"[imp8-i66] CN N={N_MAIN} dt={dt:g} s: t_call {out['t_call']:.1f} s; "
        f"n_factorizations {res.n_factorizations}; ranks {msh.comm.size}"
    )
    pts = GATE_PTS + RATIO_PTS
    th = _theta_pts(res.temperature, pts)  # one batch, I56's order
    out["theta"] = dict(zip(pts, th))
    for (r, z), v in zip(pts, th):
        ref = I56_THETA[(r, z)]
        cf = QP * _phi_tips(r * MM, z * MM, T_PRINT)
        tol = 10.0 ** (math.floor(math.log10(abs(ref))) - 5)
        _say(
            f"[imp8-i66] ({r:g}, {z:g}) mm t*={T_PRINT:g} s arm {arm}: theta {v:.6e} K/W; I56 "
            f"{ref:.6e}; theta/theta_56 {v / ref:.7f} (|diff| {abs(v - ref):.3e}, 6th-digit unit "
            f"{tol:.0e}); vs q'/(4 pi k) Phi_2 {cf:.6e}: ratio {v / cf:.6f} (not a gate: the "
            f"comparand has no rod)"
        )
    tip = out["theta"][(2.0, 15.0)] / out["theta"][(2.0, 0.0)]
    tip56 = I56_THETA[(2.0, 15.0)] / I56_THETA[(2.0, 0.0)]
    move = out["theta"][(2.0, 0.0)] / I56_THETA[(2.0, 0.0)] - 1.0
    out["move"] = move
    _say(
        f"[imp8-i66] tip-point ratio theta(2,15)/theta(2,0): arm {arm} {tip:.6f}; I56 (gel) "
        f"{tip56:.6f}"
    )
    _say(
        f"[imp8-i66] (2, 0) mm move theta/theta_56 - 1 = {move:.4e} (control floor "
        f"{CONTROL_FLOOR:.0e}, asserted on arm ti; axial-conductance prediction approx "
        f"{MOVE_PREDICTED} printed only)"
    )
    if arm == "hole":
        _say(
            f"[imp8-i69] arm hole: rod cells k = {EPS_HOLE:.0e} K_GEL, rho_c = {EPS_HOLE:.0e} "
            f"RHO_C"
        )
        for pt in pts:
            vh, vt, vg = out["theta"][pt], TI_THETA[pt], I56_THETA[pt]
            ranked = sorted((("hole", vh), ("Ti", vt), ("gel", vg)), key=lambda x: -x[1])
            order = " > ".join(name for name, _ in ranked)
            _say(
                f"[imp8-i69] ({pt[0]:g}, {pt[1]:g}) mm t*={T_PRINT:g} s: theta hole {vh:.6e}; Ti "
                f"(item 846 W1) {vt:.6e}; gel (I56) {vg:.6e} K/W; hole/gel {vh / vg:.6f}; "
                f"ordering {order}"
            )
        _say(
            "[imp8-i69] I64 closed form at r = a, t_c and t*: hole > Ti > gel "
            "(`20261105T181438Z_IMP-8-item843-W1.log:80-81`; infinite cylinder, shell source "
            "- not the same geometry)"
        )
    _say(
        f"[imp8-i66] (E) {out['E']:.3e} (band {E_BAND:.0e}); (S0) {out['S0']:.3e} (band "
        f"{B_BAND:.0e}); cells {n} (W0c {W0C_CELLS}); non-finite coefficient dofs {n_bad}"
    )
    return out


def test_w_preconditions(body_run):
    assert body_run["n"] == W0C_CELLS, body_run["n"]
    assert body_run["n_rod"] > 0
    assert body_run["n_overlap"] > 0
    assert body_run["n_bad"] == 0


def test_w_e_balance(body_run):
    assert body_run["E"] <= E_BAND, body_run["E"]


def test_w_s0_source_power(body_run):
    assert body_run["S0"] <= B_BAND, body_run["S0"]
    for n_sel, miss in body_run["balls"]:
        assert n_sel > 0
        assert miss <= B_BAND, miss


def test_w2_gel_reproduces_i56(body_run):
    if body_run["arm"] != "gel":
        pytest.skip("W2 only (arm gel): the reproduction control")
    for pt, ref in I56_THETA.items():
        v = body_run["theta"][pt]
        tol = 10.0 ** (math.floor(math.log10(abs(ref))) - 5)
        assert abs(v - ref) <= tol, (pt, v, ref, tol)


def test_w1_ti_control_moves(body_run):
    if body_run["arm"] != "ti":
        pytest.skip("W1 only (arm ti): the negative control")
    assert abs(body_run["move"]) > CONTROL_FLOOR, body_run["move"]


def test_w3_hole_control_moves(body_run):
    """Item 851's negative control: θ_hole(2, 0) mm differs from I56's gel θ."""
    if body_run["arm"] != "hole":
        pytest.skip("W3 only (arm hole): the negative control")
    assert abs(body_run["move"]) > CONTROL_FLOOR, body_run["move"]
