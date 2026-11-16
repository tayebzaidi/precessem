"""IMP-8 step 14 / §10 I71 (§9 item 852): the core-source composite cylinder in FEM.

Brief: docs/planning/chunks/IMP-8.md, "§9 item 852", under the "Review ruling
2026-11-05 18:00" paragraph.  Item 837's fixture, unchanged and imported
(``_composite_cylinder``: Ti core tag 1, gel shell tag 2, curved P2, Robin wall
at ``H_DIR`` on r = b, insulated ends; its CN settings and ``snapshot_steps``),
real build, degree 2.  The only physics change from item 837: the source is
q = ``Q_CORE`` = 1/(πa²) in the core (1 W per metre), zero in the shell.
Comparand: item 847's composite series (``_theta``, ``_theta_ss``), imported
unmoved.

Asserted (bands imported, unmoved): (I) at (a, 0) and ((a+b)/2, 0), at t* and
steady, and (I0) at (0, 0), t* and steady, all ≤ ``INTERIOR_BAND``; (E) ≤
``E_BAND``; (S0) ∫q dx = Q_CORE·V_core,h ≤ ``B_BAND`` and ∫q dx(2) = 0;
preconditions as item 837 (coordinate degree 2, 53 504 cells, tag-1 geometry
nodes at r = b, no tagged end face, P2 wiring).  Control (asserted, factor 1 on
``INTERIOR_BAND``): the k_r = ``K_GEL`` comparand misses θ_FEM(0, steady) by more
than the band (predicted ≈ 0.275 from item 847's arithmetic, printed).  Window
by ``FEM_EM_IMP8_I71_PHASE`` = ``solve`` (W1).
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

from fem_em_solver.post.evaluation import evaluate_vector_field_parallel
from fem_em_solver.thermal import solve_steady_conduction, solve_transient_conduction

from tests.validation.test_imp8_composite_cylinder_closed_form import (
    A_ROD,
    B_WALL,
    ROD_K,
    ROD_RHO_C,
    T_C,
    T_STAR,
)
from tests.validation.test_imp8_composite_cylinder_core_source_closed_form import (
    PRED_CTRL,
    Q_CORE,
    _theta,
    _theta_ss,
)
from tests.validation.test_imp8_composite_cylinder_fem import (
    CELL_CAP,
    DEGREE,
    H_DIR,
    LEN,
    MM,
    N_STEPS,
    READ_PTS,
    RA,
    RM,
    SNAP,
    T_INF,
    WIT_PTS,
    _composite_cylinder,
)
from tests.validation.test_imp8_rod_line_source_fem import _theta_pts
from tests.validation.test_mat10_sar_to_heat_sphere import B_BAND, ROBIN_TAG, _allsum, _say
from tests.validation.test_mat12_transient_sar_heating_sphere import E_BAND, RHO_C
from tests.validation.test_mat13_transient_perfusion_sphere import INTERIOR_BAND
from tests.validation.test_mat15_coil_driven_heating_phantom import K_GEL
from tests.validation.test_mat18_two_layer_sphere_fem import _dg0

#: Item 837's measured count (`20261105T225620Z_IMP-8-item837-W1.log`).
I65_CELLS = 53_504
PHASE_ENV = "FEM_EM_IMP8_I71_PHASE"
PHASES = ("solve",)


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
    monkeypatch.setenv(PHASE_ENV, "solve")
    assert _phase() == "solve"
    for bad in ("1", "build", "steady"):
        monkeypatch.setenv(PHASE_ENV, bad)
        with pytest.raises(ValueError):
            _phase()


@pytest.fixture(scope="module")
def run():
    if _phase() != "solve":
        pytest.skip(f"{PHASE_ENV} != 'solve'")
    comm = MPI.COMM_WORLD
    t0 = time.perf_counter()
    msh, ct, ft = _composite_cylinder(comm)
    t_build = time.perf_counter() - t0
    n = msh.topology.index_map(msh.topology.dim).size_global
    out = {"n": n, "deg": msh.geometry.cmap.degree}
    # Preconditions (item 837's): tag-1 geometry nodes at r = b, no tagged end face,
    # P2 wiring at the three read points.
    fdim = msh.topology.dim - 1
    msh.topology.create_connectivity(fdim, msh.topology.dim)
    facets = ft.find(ROBIN_TAG)
    gidx = dmesh.entities_to_geometry(msh, fdim, facets)
    x = msh.geometry.x
    rr = np.hypot(x[gidx.ravel(), 0], x[gidx.ravel(), 1]) if facets.size else np.zeros(0)
    out["dev"] = comm.allreduce(float(np.max(np.abs(rr - B_WALL))) if rr.size else 0.0,
                                op=MPI.MAX)
    out["n_fac"] = comm.allreduce(int(np.count_nonzero(
        facets < msh.topology.index_map(fdim).size_local)), op=MPI.SUM)
    mid = dmesh.compute_midpoints(msh, fdim, facets) if facets.size else np.zeros((0, 3))
    out["zmax"] = comm.allreduce(float(np.max(np.abs(mid[:, 2]))) if facets.size else 0.0,
                                 op=MPI.MAX)
    V = fem.functionspace(msh, ("Lagrange", DEGREE))
    f = fem.Function(V)
    f.interpolate(lambda X: X[0] + 2.0 * X[1] + 3.0 * X[2])
    pts = np.array([[r * MM, 0.0, z * MM] for r, z in READ_PTS], dtype=float)
    vals, ok = evaluate_vector_field_parallel(f, pts)
    out["werr"] = max(abs(float(np.real(np.ravel(vals[i])[0])) - (p[0] + 2 * p[1] + 3 * p[2]))
                      if bool(ok[i]) else np.inf for i, p in enumerate(pts))
    # Coefficients: item 837's, the source moved to the core.
    k = _dg0(msh, ct, ROD_K, K_GEL)
    rc = _dg0(msh, ct, ROD_RHO_C, RHO_C)
    q = _dg0(msh, ct, Q_CORE, 0.0)
    nan_loc = int(sum(np.count_nonzero(np.isnan(np.real(g.x.array))) for g in (k, rc, q)))
    out["nan"] = comm.allreduce(nan_loc, op=MPI.SUM)
    dx = ufl.Measure("dx", domain=msh, subdomain_data=ct)
    v_core = _allsum(comm, fem.Constant(msh, 1.0) * dx(1))
    out["S0"] = abs(_allsum(comm, q * dx) / (Q_CORE * v_core) - 1.0)
    out["S0shell"] = abs(_allsum(comm, q * dx(2)))
    _say(f"\n[imp8-i71] composite cylinder, core source: cells {n} (I65 {I65_CELLS}, cap "
         f"{CELL_CAP}); coordinate degree {out['deg']}; ranks {comm.size}; t_build {t_build:.1f} s; "
         f"tag-{ROBIN_TAG} owned facets {out['n_fac']}, max |r - b| {out['dev']:.3e} m, max |z| "
         f"{out['zmax'] / MM:.4f} mm; P2 wiring max err {out['werr']:.3e}")
    _say(f"[imp8-i71] V_core,h/(pi a^2 L) - 1 = {v_core / (np.pi * A_ROD**2 * LEN) - 1:.4e}; "
         f"Q_CORE {Q_CORE:.6e} W/m^3")
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
    _say(f"[imp8-i71] t_call transient {t_tr:.1f} s (CN N={N_STEPS}, dt={dt:.6g} s, snapshot "
         f"{SNAP} -> t = {SNAP * dt:.6g} s = t_c {T_C:.6g}), n_factorizations "
         f"{res.n_factorizations}; t_call steady {t_st:.1f} s; NaN left {out['nan']}")
    th_c = _theta_pts(res.snapshots[SNAP], READ_PTS)
    th_s = _theta_pts(res.temperature, READ_PTS + WIT_PTS)
    th_ss = _theta_pts(st.temperature, READ_PTS)
    ex = {
        "t_c": [_theta(r * MM, T_C) for r, _ in READ_PTS],
        "t*": [_theta(r * MM, T_STAR) for r, _ in READ_PTS],
        "ss": [float(_theta_ss(r * MM)) for r, _ in READ_PTS],
    }
    fem_ = {"t_c": th_c, "t*": th_s[:3], "ss": th_ss}
    out["I"], out["I0"] = [], []
    for key in ("t_c", "t*", "ss"):
        for (r, z), v, e in zip(READ_PTS, fem_[key], ex[key]):
            rel = abs(v / e - 1.0)
            gated = key in ("t*", "ss")
            tag = "(I0)" if r == 0.0 else "(I)"
            if gated:
                out["I0" if r == 0.0 else "I"].append(rel)
            _say(f"[imp8-i71] {tag} ({r:g}, {z:g}) mm {key}: theta_FEM {v:.6e}, composite "
                 f"{e:.6e} K per W/m, rel {rel:.3e} "
                 f"({'gated' if gated else 'printed only'}, band {INTERIOR_BAND:.0e})")
    _say(f"[imp8-i71] z-invariance at t*: theta_FEM(a, 0) {th_s[0]:.6e}, (a, +L/4) "
         f"{th_s[3]:.6e}, (a, -L/4) {th_s[4]:.6e}")
    c0 = float(_theta_ss(0.0, k_r=K_GEL))
    out["ctrl"] = abs(th_ss[2] / c0 - 1.0)
    _say(f"[imp8-i71] control (k_r = K_GEL) at (0, steady): comparand {c0:.6e}, miss "
         f"{out['ctrl']:.4e} (asserted > {INTERIOR_BAND:.0e}; predicted approx {PRED_CTRL}, "
         f"not asserted)")
    ds = ufl.Measure("ds", domain=msh, subdomain_data=ft)
    area = _allsum(comm, fem.Constant(msh, 1.0) * ds(ROBIN_TAG))
    t_wall = _allsum(comm, st.temperature * ds(ROBIN_TAG)) / area
    _say(f"[imp8-i71] Robin offset: mean wall theta / theta_ss(a) = "
         f"{t_wall / float(_theta_ss(A_ROD)):.4e}; steady convected / source "
         f"{st.convected_power / st.source_power - 1:.3e}; (E) {out['E']:.3e} (band "
         f"{E_BAND:.0e}); (S0) {out['S0']:.3e} (band {B_BAND:.0e}), shell {out['S0shell']:.3e}")
    return out


def test_w1_preconditions(run):
    assert run["deg"] == 2, run["deg"]
    assert run["n"] == I65_CELLS, run["n"]
    assert run["n_fac"] > 0
    assert run["dev"] <= 1e-6 * B_WALL, run["dev"]
    assert run["zmax"] < 0.5 * LEN - 1e-6, run["zmax"]
    assert run["werr"] <= 1e-12, run["werr"]
    assert run["nan"] == 0, run["nan"]


def test_w1_s0_source(run):
    assert run["S0"] <= B_BAND, run["S0"]
    assert run["S0shell"] == 0.0, run["S0shell"]


def test_w1_e_balance(run):
    assert run["E"] <= E_BAND, run["E"]


def test_w1_i_interior(run):
    assert len(run["I"]) == 4
    assert max(run["I"]) <= INTERIOR_BAND, run["I"]


def test_w1_i0_axis(run):
    assert len(run["I0"]) == 2
    assert max(run["I0"]) <= INTERIOR_BAND, run["I0"]


def test_w1_control_k_gel(run):
    """Asserted at factor 1 (> ``INTERIOR_BAND``); the ≈ 0.275 is predicted, printed."""
    assert run["ctrl"] > INTERIOR_BAND, run["ctrl"]
