"""`MAT-21` step 29 — the pre-equilibration dose on the FEM history, the FEM half of P55
(§9 item 736, §10 chain P, P56, provisional under rule (6); brief
`docs/planning/chunks/MAT-21.md` "§9 item 736").

Sets no band and moves none; every band is imported.  P52's three runs on ``RUNGS[-1]`` (the
reference from uniform T∞, t_w = 0 from T∞ + g·z, t_w = ``T_C`` from T∞ + the drift series),
called as P52 calls them, with ``snapshot_steps=tuple(range(N_MAIN + 1))``.  Each run's
snapshots are read on one fixed batch, ``HIST_PTS`` = P0, P5, ``ADJ_PT`` (item 537), then
cleared, so one run's snapshots are held at a time; ``PROBE_PTS`` stays P52's own batch on the
final field for (Sym′).  Real build.

Asserted: (I) ≤ ``INTERIOR_BAND`` on both gradient runs; (E) ≤ ``E_BAND`` on all three, on the
source branch; (Sym′) as P52; (Q-wf) ``cem43`` against ``_q_check_g``'s quadrature at
``DOSE_BAR`` on each FEM history.  Negative control: the t_w = 0 FEM dose at P5 misses the
reference run's by > ``DOSE_BAR``.  Control: the reference reproduces ``STEP1_T`` at
``STEP1_T_BOUND``.  Printed (**no band**): D at P0, P5, ``ADJ_PT``, FEM vs item 727's series,
per run; the FEM-vs-series gap on one rung and one Δt, no rate.

Not claimed: a dose band, a rate, a recommended wait.
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
from tests.validation.test_mat21_box_closed_form import A_BOX, B_BOX, C_BOX, T_C, Q
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
from tests.validation.test_mat21_box_probe_points_fem import PROBE_PTS
from tests.validation.test_mat21_initial_gradient_closed_form import G_GRAD, GC
from tests.validation.test_mat21_initial_gradient_fem import SYM_BAR, _grad_ex
from tests.validation.test_mat21_pre_equilibration_closed_form import (
    ADJ_PT,
    P5,
)
from tests.validation.test_mat21_pre_equilibration_closed_form import T_INF as T_INF_51
from tests.validation.test_mat21_pre_equilibration_closed_form import _combined, _sourced
from tests.validation.test_mat21_pre_equilibration_dose_closed_form import GRID, _dose
from tests.validation.test_mat21_pre_equilibration_fem import (
    STEP1_T,
    STEP1_T_BOUND,
    _read,
    _wait_init,
)

HIST_LABELS = ("P0", "P5", "ADJ_PT")
HIST_PTS = np.array((PROBE_PTS[0], P5, ADJ_PT), dtype=float)
T_FEM = (T_STAR / N_MAIN) * np.arange(N_MAIN + 1)
WAIT = {"tw0": 0.0, "twC": T_C}
#: P55's D at P5, t_w = 0, `20261101T195309Z_MAT-21-step28-W1.log:61`.
P55_D_P5_T0 = "1.275743e-10"
#: P52's (I), `docs/planning/chunks/MAT-21.md` step 26 (printed beside, no assertion).
P52_I = {"I0": "1.125097e-03", "IC": "9.076702e-04"}


def _series(p, key, t):
    """Item 727's closed-form history at ``p`` on grid ``t``; the sourced-only one for "ref"."""
    if key == "ref":
        return np.array([_sourced(p, s) for s in t])
    return np.array([_combined(p, s, WAIT[key]) for s in t])


@pytest.fixture(scope="module")
def runs():
    comm = MPI.COMM_WORLD
    t0 = time.perf_counter()
    m = RUNGS[-1]
    ns = [b * m for b in BASE]
    msh = dmesh.create_box(
        comm,
        [np.array([-A_BOX, -B_BOX, -C_BOX]), np.array([A_BOX, B_BOX, C_BOX])],
        ns,
        cell_type=dmesh.CellType.tetrahedron,
    )
    ft = _robin_tags(msh)
    chi = _interior_chi(msh)
    pts = np.array(PROBE_PTS, dtype=float)
    _say(
        f"\nMAT21-29 rung m = {m} {tuple(ns)}: cells {msh.topology.index_map(3).size_global}; "
        f"g = {G_GRAD} K/m; CN N = {N_MAIN} + 2 BE half-steps, {N_MAIN + 1} snapshots per run; "
        f"history batch {HIST_LABELS}; ranks {comm.size}"
    )

    kw = dict(
        k=K_GEL,
        rho_c=RHO_C,
        q=Q,
        h=H_GEL,
        t_inf=T_INF,
        t_init=T_INF,
        robin_tag=ROBIN_TAG,
        dt=T_STAR / N_MAIN,
        n_steps=N_MAIN,
        theta=0.5,
        startup_half_steps=2,
        degree=1,
        snapshot_steps=tuple(range(N_MAIN + 1)),
    )
    starts = {"ref": None, "tw0": lambda x: T_INF + G_GRAD * x[2], "twC": _wait_init(T_C)}
    out = {"E": {}, "src": {}, "th": {}, "ok": {}, "hist": {}, "hist_ok": {}, "last": {}}
    res = {}
    for key, expr in starts.items():
        tq = time.perf_counter()
        res[key] = solve_transient_conduction(msh, ft, t_init_expr=expr, **kw)
        ts = time.perf_counter()
        out["E"][key] = float(res[key].max_balance_residual)
        out["src"][key] = float(res[key].source_power)
        out["th"][key], out["ok"][key] = _read(res[key].temperature, pts)

        # One batch, one order, on every snapshot; read and clear.
        hist = np.empty((N_MAIN + 1, HIST_PTS.shape[0]))
        ok_all = [True] * HIST_PTS.shape[0]
        for n in range(N_MAIN + 1):
            hist[n], ok = _read(res[key].snapshots[n], HIST_PTS)
            ok_all = [a and b for a, b in zip(ok_all, ok)]
        res[key].snapshots.clear()
        out["hist"][key], out["hist_ok"][key] = hist, ok_all
        out["last"][key] = _read(res[key].temperature, HIST_PTS)[0]
        _say(
            f"MAT21-29 run {key}: solve {ts - tq:.1f} s, reads {time.perf_counter() - ts:.1f} s; "
            f"(E) {out['E'][key]:.3e} (E_BAND {E_BAND:.0e}); source_power "
            f"{out['src'][key]:.6e} W"
        )

    tq = time.perf_counter()
    ex_s, npts = _exact(msh, _trans_ex(T_STAR))
    ex_0, _ = _exact(msh, _grad_ex(T_STAR))
    ex_c, _ = _exact(msh, _grad_ex(T_STAR + T_C))
    _say(
        f"MAT21-29 three series interpolations ({npts} P2 points): "
        f"{time.perf_counter() - tq:.1f} s"
    )

    th = {key: r.temperature - T_INF for key, r in res.items()}
    out["ref_T"] = _rel(msh, th["ref"], ex_s, chi)
    out["I0"] = _rel(msh, th["tw0"], ex_s + ex_0, chi)
    out["IC"] = _rel(msh, th["twC"], ex_s + ex_c, chi)
    _say(
        f"MAT21-29 (I) t_w = 0: interior rel L2 {out['I0']:.6e} (P52 {P52_I['I0']}); t_w = T_C: "
        f"{out['IC']:.6e} (P52 {P52_I['IC']}); bound {INTERIOR_BAND:.0e}"
    )
    _say(
        f"MAT21-29 control (ref) reference interior {out['ref_T']:.6e} vs step 1 {STEP1_T:.6e}, "
        f"|delta| {abs(out['ref_T'] - STEP1_T):.3e} (bound {STEP1_T_BOUND:.0e})"
    )
    _say(f"MAT21-29 elapsed {time.perf_counter() - t0:.1f} s")
    return out


def test_w0_imports_and_p55_dose_at_p5():
    """W0 (numpy): the imports and item 727's D at P5, t_w = 0, with no solve."""
    d = _dose(_series(P5, "tw0", GRID))
    _say(
        f"MAT21-29 W0 P5, t_w = 0: D {d:.6e} s (P55 {P55_D_P5_T0}); {N_MAIN} steps, "
        f"series rows {GRID.size}, T_C {T_C:.3f} s, DOSE_BAR {DOSE_BAR:g}"
    )
    assert f"{d:.6e}" == P55_D_P5_T0, d
    assert N_MAIN == 40 and RUNGS[-1] == 6
    assert T_INF == T_INF_51
    assert abs(T_FEM[-1] / GRID[-1] - 1.0) <= 1e-14
    assert HIST_PTS.shape == (3, 3) and HIST_PTS[0, 2] == 0.0 and HIST_PTS[1, 2] > 0.0


def test_preconditions_every_history_read_and_n40_is_the_final_field(runs):
    for key, ok in runs["hist_ok"].items():
        assert all(ok), (key, ok)
        assert np.array_equal(runs["hist"][key][N_MAIN], runs["last"][key]), key


def test_i_both_gradient_runs(runs):
    assert runs["I0"] <= INTERIOR_BAND, runs["I0"]
    assert runs["IC"] <= INTERIOR_BAND, runs["IC"]


def test_e_all_three_on_the_source_branch(runs):
    src = runs["src"]
    assert src["ref"] > 0.0 and src["ref"] == src["tw0"] == src["twC"], src
    for key, e in runs["E"].items():
        assert e <= E_BAND, (key, e)


def test_sym_prime_as_p52(runs):
    assert all(all(o) for o in runs["ok"].values()), runs["ok"]
    worst = {}
    for key in ("tw0", "twC"):
        d = runs["th"][key] - runs["th"]["ref"]
        s12, s34, s0 = abs(d[1] + d[2]), abs(d[3] + d[4]), abs(d[0])
        worst[key] = max(s12, s34, s0)
        _say(
            f"MAT21-29 (Sym') {key}: |P1 + P2|/(g C) {s12 / GC:.3e}, |P3 + P4|/(g C) "
            f"{s34 / GC:.3e}, |P0|/(g C) {s0 / GC:.3e} (bar {SYM_BAR:.0e})"
        )
    for key, w in worst.items():
        assert w <= SYM_BAR * GC, (key, w)


def test_q_wf_on_each_fem_history(runs):
    rows = []
    for key, hist in runs["hist"].items():
        for j, lab in enumerate(HIST_LABELS):
            dose, quad, rel, _swap, miss = _q_check_g(T_INF + hist[:, j], T_FEM)
            _say(
                f"MAT21-29 (Q-wf) {key} {lab}: D {dose:.6e} s, quad {quad:.6e} s, rel {rel:.3e} "
                f"(bar {DOSE_BAR:g}, ASSERTED); R swap miss {miss:.4e} (printed)"
            )
            rows.append((key, lab, rel))
    for key, lab, rel in rows:
        assert rel <= DOSE_BAR, (key, lab, rel)


def test_control_reference_reproduces_step_1(runs):
    assert abs(runs["ref_T"] - STEP1_T) <= STEP1_T_BOUND, runs["ref_T"]


def test_control_gradient_dose_misses_the_reference_at_p5(runs):
    j = HIST_LABELS.index("P5")
    d_0 = _dose(runs["hist"]["tw0"][:, j], T_FEM)
    d_r = _dose(runs["hist"]["ref"][:, j], T_FEM)
    miss = abs(d_0 / d_r - 1.0)
    _say(
        f"MAT21-29 control P5: t_w = 0 FEM D {d_0:.6e} s vs reference {d_r:.6e} s, miss "
        f"{miss:.4e} (asserted > {DOSE_BAR:g})"
    )
    assert miss > DOSE_BAR, miss


def test_printed_dose_table(runs):
    """The record (no band); asserts only that it was formed."""
    _say(
        "MAT21-29 TABLE point | run | D FEM s | D series 41 rows s | FEM/series41-1 | "
        "D series 161 rows s | FEM/series161-1 | FEM D/D_ref-1 | series41 D/D_ref-1 | "
        "series161 D/D_ref-1"
    )
    for j, (lab, p) in enumerate(zip(HIST_LABELS, HIST_PTS)):
        d_h = {key: _dose(runs["hist"][key][:, j], T_FEM) for key in runs["hist"]}
        d_41 = {key: _dose(_series(tuple(p), key, T_FEM), T_FEM) for key in runs["hist"]}
        d_161 = {key: _dose(_series(tuple(p), key, GRID)) for key in runs["hist"]}
        for key in runs["hist"]:
            _say(
                f"MAT21-29 TABLE {lab} | {key} | {d_h[key]:.6e} | {d_41[key]:.6e} | "
                f"{d_h[key] / d_41[key] - 1.0:+.3e} | {d_161[key]:.6e} | "
                f"{d_h[key] / d_161[key] - 1.0:+.3e} | {d_h[key] / d_h['ref'] - 1.0:+.6e} | "
                f"{d_41[key] / d_41['ref'] - 1.0:+.6e} | {d_161[key] / d_161['ref'] - 1.0:+.6e}"
            )
            assert np.isfinite(d_h[key] / d_161[key])
