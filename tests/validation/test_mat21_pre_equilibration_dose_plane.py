"""`MAT-21` step 30 — the FEM-vs-series dose gap on a mesh plane and midway (§9 item 742, §10
chain P, P57, provisional under rule (6); brief `docs/planning/chunks/MAT-21.md` "§9 item 742").

Sets no band and moves none; every band is imported.  Two of item 736's runs on ``RUNGS[-1]``,
reference and t_w = 0, with ``snapshot_steps=tuple(range(N_MAIN + 1))``.  Each run's snapshots
are read on one fixed batch, ``(PLANE_PT, ADJ_PT, MID_PT)`` (item 537), then cleared, so one
run's snapshots are held at a time.  Real build.

Asserted: W0, ``PLANE_PT``'s z on a mesh plane and ``MID_PT``'s midway, each to 1e-12 relative,
on Δz = C / 6; (I) ≤ ``INTERIOR_BAND`` on t_w = 0; (E) ≤ ``E_BAND`` on both runs; (Q-wf)
``cem43`` against ``_q_check_g``'s quadrature at ``DOSE_BAR`` on all six histories.  Negative
control: ``ADJ_PT``'s FEM doses reproduce item 736's records within rel 1e-6, and the t_w = 0
dose at ``PLANE_PT`` misses the reference's by > ``DOSE_BAR``.

Printed (**no band**): the signed gap FEM / series − 1 at the three heights on both runs, beside
½ s (Δz − s) T″ at the two off-plane heights, with s the distance to the lower plane and T″ the
series' z-curvature at t*.  Also printed, beyond the brief: the temperature gap at t*, and the
dose gap of the series' own P1 interpolant between the two bracketing planes.

Not claimed: a dose band, a rate.
"""

from __future__ import annotations

import time

import numpy as np
import pytest
from dolfinx import mesh as dmesh
from mpi4py import MPI

from fem_em_solver.thermal import solve_transient_conduction

from tests.validation.test_mat10_sar_to_heat_sphere import ROBIN_TAG, _robin_tags, _say
from tests.validation.test_mat12_transient_sar_heating_sphere import E_BAND, RHO_C
from tests.validation.test_mat13_transient_perfusion_sphere import INTERIOR_BAND
from tests.validation.test_mat15_coil_driven_heating_phantom import H_GEL, K_GEL, N_MAIN, T_INF
from tests.validation.test_mat15_finite_cylinder_fem import _rel
from tests.validation.test_mat21_box_closed_form import A_BOX, B_BOX, C_BOX, Q
from tests.validation.test_mat21_box_dose_closed_form import DOSE_BAR
from tests.validation.test_mat21_box_dt_duty125_dose_closed_form import _q_check_g
from tests.validation.test_mat21_box_fem import (
    BASE,
    RUNGS,
    T_STAR,
    _exact,
    _interior_chi,
    _trans_ex,
)
from tests.validation.test_mat21_initial_gradient_closed_form import G_GRAD
from tests.validation.test_mat21_initial_gradient_fem import _grad_ex
from tests.validation.test_mat21_pre_equilibration_closed_form import ADJ_PT, _combined, _sourced
from tests.validation.test_mat21_pre_equilibration_dose_closed_form import GRID, _dose
from tests.validation.test_mat21_pre_equilibration_dose_fem import T_FEM, _series
from tests.validation.test_mat21_pre_equilibration_fem import _read

DZ = C_BOX / 6.0
PLANE_PT = (0.0, 0.0, 5.0 * C_BOX / 6.0)
MID_PT = (0.0, 0.0, 11.0 * C_BOX / 12.0)
TOP_PT = (0.0, 0.0, C_BOX)
HIST_LABELS = ("PLANE_PT", "ADJ_PT", "MID_PT")
HIST_PTS = np.array((PLANE_PT, ADJ_PT, MID_PT), dtype=float)
KEYS = ("ref", "tw0")
#: Item 736's FEM doses at ``ADJ_PT``, `20261101T224626Z_MAT-21-step29-W1.log:90–98`.
ITEM736_ADJ_D = {"ref": 1.008673e-10, "tw0": 1.363467e-10}
REPRO_BOUND = 1e-6  # the brief's, one unit in the print's last digit
PLANE_BOUND = 1e-12  # the brief's W0 bound, relative to C


def _box(comm):
    ns = [b * RUNGS[-1] for b in BASE]
    return dmesh.create_box(comm, [np.array([-A_BOX, -B_BOX, -C_BOX]),
                                   np.array([A_BOX, B_BOX, C_BOX])],
                            ns, cell_type=dmesh.CellType.tetrahedron), ns


def _at_t(p, key, t):
    return _sourced(p, t) if key == "ref" else _combined(p, t, 0.0)


@pytest.fixture(scope="module")
def runs():
    comm = MPI.COMM_WORLD
    t0 = time.perf_counter()
    msh, ns = _box(comm)
    ft = _robin_tags(msh)
    chi = _interior_chi(msh)
    _say(f"\nMAT21-30 rung m = {RUNGS[-1]} {tuple(ns)}: cells "
         f"{msh.topology.index_map(3).size_global}; CN N = {N_MAIN} + 2 BE half-steps, "
         f"{N_MAIN + 1} snapshots per run; batch {HIST_LABELS}; ranks {comm.size}")

    kw = dict(k=K_GEL, rho_c=RHO_C, q=Q, h=H_GEL, t_inf=T_INF, t_init=T_INF, robin_tag=ROBIN_TAG,
              dt=T_STAR / N_MAIN, n_steps=N_MAIN, theta=0.5, startup_half_steps=2, degree=1,
              snapshot_steps=tuple(range(N_MAIN + 1)))
    starts = {"ref": None, "tw0": lambda x: T_INF + G_GRAD * x[2]}
    out = {"E": {}, "src": {}, "hist": {}, "hist_ok": {}}
    res = {}
    for key, expr in starts.items():
        tq = time.perf_counter()
        res[key] = solve_transient_conduction(msh, ft, t_init_expr=expr, **kw)
        ts = time.perf_counter()
        out["E"][key] = float(res[key].max_balance_residual)
        out["src"][key] = float(res[key].source_power)

        # One batch, one order, on every snapshot; read and clear.
        hist = np.empty((N_MAIN + 1, HIST_PTS.shape[0]))
        ok_all = [True] * HIST_PTS.shape[0]
        for n in range(N_MAIN + 1):
            hist[n], ok = _read(res[key].snapshots[n], HIST_PTS)
            ok_all = [a and b for a, b in zip(ok_all, ok)]
        res[key].snapshots.clear()
        out["hist"][key], out["hist_ok"][key] = hist, ok_all
        _say(f"MAT21-30 run {key}: solve {ts - tq:.1f} s, reads {time.perf_counter() - ts:.1f} s; "
             f"(E) {out['E'][key]:.3e} (E_BAND {E_BAND:.0e}); source_power "
             f"{out['src'][key]:.6e} W")

    tq = time.perf_counter()
    ex_s, npts = _exact(msh, _trans_ex(T_STAR))
    ex_0, _ = _exact(msh, _grad_ex(T_STAR))
    _say(f"MAT21-30 two series interpolations ({npts} P2 points): "
         f"{time.perf_counter() - tq:.1f} s")
    out["I0"] = _rel(msh, res["tw0"].temperature - T_INF, ex_s + ex_0, chi)
    _say(f"MAT21-30 (I) t_w = 0: interior rel L2 {out['I0']:.6e}; bound {INTERIOR_BAND:.0e}")
    _say(f"MAT21-30 elapsed {time.perf_counter() - t0:.1f} s")
    return out


def test_w0_heights_on_and_midway_between_mesh_planes():
    """W0: the mesh's own z coordinates, no solve, on every rank's vertices."""
    comm = MPI.COMM_WORLD
    msh, ns = _box(comm)
    z = np.unique(np.concatenate(comm.allgather(np.unique(msh.geometry.x[:, 2]))))
    d_plane = float(np.min(np.abs(z - PLANE_PT[2]))) / C_BOX
    lo, hi = float(z[z <= MID_PT[2]].max()), float(z[z >= MID_PT[2]].min())
    d_mid = abs(MID_PT[2] - 0.5 * (lo + hi)) / C_BOX
    d_dz = abs((hi - lo) / DZ - 1.0)
    _say(f"MAT21-30 W0 nz {ns[2]}: PLANE_PT off a plane {d_plane:.3e} C; MID_PT off midway "
         f"{d_mid:.3e} C; dz / (C / 6) - 1 {d_dz:.3e} (bound {PLANE_BOUND:g})")
    assert d_dz <= PLANE_BOUND, (lo, hi)
    assert d_plane <= PLANE_BOUND, d_plane
    assert d_mid <= PLANE_BOUND, d_mid
    assert lo < ADJ_PT[2] < hi and abs(lo - PLANE_PT[2]) / C_BOX <= PLANE_BOUND
    assert abs(hi / C_BOX - 1.0) <= PLANE_BOUND  # z = C is the upper plane of both off-plane reads


def test_preconditions_every_history_read(runs):
    for key, ok in runs["hist_ok"].items():
        assert all(ok), (key, ok)


def test_i_on_t_w_0(runs):
    assert runs["I0"] <= INTERIOR_BAND, runs["I0"]


def test_e_both_runs_on_the_source_branch(runs):
    assert runs["src"]["ref"] > 0.0 and runs["src"]["ref"] == runs["src"]["tw0"], runs["src"]
    for key, e in runs["E"].items():
        assert e <= E_BAND, (key, e)


def test_q_wf_on_all_six_histories(runs):
    rows = []
    for key, hist in runs["hist"].items():
        for j, lab in enumerate(HIST_LABELS):
            dose, quad, rel, _swap, miss = _q_check_g(T_INF + hist[:, j], T_FEM)
            _say(f"MAT21-30 (Q-wf) {key} {lab}: D {dose:.6e} s, quad {quad:.6e} s, rel {rel:.3e} "
                 f"(bar {DOSE_BAR:g}, ASSERTED); R swap miss {miss:.4e} (printed)")
            rows.append((key, lab, rel))
    for key, lab, rel in rows:
        assert rel <= DOSE_BAR, (key, lab, rel)


def test_control_adj_pt_reproduces_item_736(runs):
    j = HIST_LABELS.index("ADJ_PT")
    for key, rec in ITEM736_ADJ_D.items():
        d = _dose(runs["hist"][key][:, j], T_FEM)
        rel = abs(d / rec - 1.0)
        _say(f"MAT21-30 control ADJ_PT {key}: D {d:.6e} s vs item 736 {rec:.6e} s, rel {rel:.3e} "
             f"(bound {REPRO_BOUND:g})")
        assert rel <= REPRO_BOUND, (key, d)


def test_control_gradient_dose_misses_the_reference_at_plane(runs):
    j = HIST_LABELS.index("PLANE_PT")
    d_0 = _dose(runs["hist"]["tw0"][:, j], T_FEM)
    d_r = _dose(runs["hist"]["ref"][:, j], T_FEM)
    miss = abs(d_0 / d_r - 1.0)
    _say(f"MAT21-30 control PLANE_PT: t_w = 0 FEM D {d_0:.6e} s vs reference {d_r:.6e} s, miss "
         f"{miss:.4e} (asserted > {DOSE_BAR:g})")
    assert miss > DOSE_BAR, miss


def test_printed_gap_table(runs):
    """The record (no band); asserts only that it was formed."""
    h = DZ / 4.0  # the curvature's stencil; stays strictly inside the upper cell's z < C
    _say("MAT21-30 TABLE height | run | D FEM s | D series161 s | FEM/series161-1 | D series41 s | "
         "FEM/series41-1 | P1 interpolant/series41-1 | gap41 minus plane gap41 | T gap at t* K | "
         "1/2 s (dz-s) T'' K")
    for key in KEYS:
        lo = _series(PLANE_PT, key, T_FEM)
        hi = _series(TOP_PT, key, T_FEM)
        plane_gap = None
        for j, (lab, p) in enumerate(zip(HIST_LABELS, HIST_PTS)):
            p = tuple(float(c) for c in p)
            fem = runs["hist"][key][:, j]
            s41, s161 = _series(p, key, T_FEM), _series(p, key, GRID)
            d_h, d_41, d_161 = _dose(fem, T_FEM), _dose(s41, T_FEM), _dose(s161)
            gap41 = d_h / d_41 - 1.0
            plane_gap = gap41 if lab == "PLANE_PT" else plane_gap
            t_gap = float(fem[-1] - s41[-1])
            if lab == "PLANE_PT":
                est, interp = 0.0, 0.0
            else:
                s = p[2] - PLANE_PT[2]
                up, dn = (p[0], p[1], p[2] + h), (p[0], p[1], p[2] - h)
                curv = (_at_t(up, key, T_STAR) - 2.0 * s41[-1] + _at_t(dn, key, T_STAR)) / h ** 2
                est = 0.5 * s * (DZ - s) * curv
                a = s / DZ
                interp = _dose((1.0 - a) * lo + a * hi, T_FEM) / d_41 - 1.0
            _say(f"MAT21-30 TABLE {lab} | {key} | {d_h:.6e} | {d_161:.6e} | "
                 f"{d_h / d_161 - 1.0:+.3e} | {d_41:.6e} | {gap41:+.3e} | {interp:+.3e} | "
                 f"{gap41 - plane_gap:+.3e} | {t_gap:+.3e} | {est:+.3e}")
            assert np.isfinite(gap41)
