"""Example (`EX-176`): the same skipped wait under two schedules -- the FEM dose of the
sourced gel box at D = 1 and at D = 0.25, reference against t_w = 0.

`mat:34` prices the skipped wait in dose at D = 1 (the source on for the whole run). Here the
four runs of `MAT-21` steps 27-33 are made on the same box (the finest rung of `MAT-21` step 1,
m = 6): **reference** from uniform T_inf and **t_w = 0** from T_inf + g z, each under the
continuous source (D = 1) and under the D = 0.25 schedule ``_chi25`` (item 753's, 40 steps of
``DT_BOX``). Each run's per-step history is read at P5 and ``ADJ_PT`` on one fixed batch, and the
CEM43 dose of each history is printed beside the closed-form dose, as a table of one wait under
two schedules. A dose at D = 0.25 is far smaller in absolute terms (it is a fourth of the
energy), so only the ratio D(t_w = 0) / D_ref - 1 is compared across the schedules.

**It asserts, and it does not re-implement** (the `ANS-1` rule). Bands, helpers, the series and
the schedule are imported by name from the modules the two FEM gates
``tests/validation/test_mat21_pre_equilibration_dose_fem.py`` (D = 1) and
``tests/validation/test_mat21_pre_equilibration_duty_dose_fem.py`` (D = 0.25, `MAT-21` step 33)
draw on. Those gates' ``runs`` fixtures are pytest fixtures and are not called: their bodies are
written out here, plus ``np.real`` on each integral (the group sources the complex build) and a
``write_setup_figure`` call after the mesh. The two modules' ``_dose`` differ (a grid argument
against none) and are imported under their own names. No gate edit, no ``src/`` edit.

Asserted (bands imported; every gate value collected first, then asserted): (I) interior rel L2
at t* against the superposition <= ``INTERIOR_BAND`` on both t_w = 0 runs and on the D = 0.25
reference (its comparand the duty sum without the drift); (E) the balance <= ``E_BAND`` on all
four, ``source_power`` equal within each schedule; (S0) ``_census`` at ``S0_BAND`` on the recorded
multipliers and the schedule string on the two D = 0.25 runs; (Q-wf) ``cem43`` against the
quadrature ``_q_check_g`` at ``DOSE_BAR`` on all eight histories; (c-dose) at each schedule the
t_w = 0 FEM dose at P5 misses the reference's by more than ``DOSE_BAR``; (c-cont) the continuous
superposition misses (I) on the D = 0.25 t_w = 0 run by more than ``INTERIOR_BAND``. Printed only
(records, no band): D_ref, D(t_w = 0), D / D_ref - 1 on the FEM and on the series. No dose band,
no rate, no recommended wait.

Run it::

    ./run_examples.sh -e mat:36 -n 2 -t 480

Output: ``examples/materials/paraview_output/
materials_36_gel_box_wait_across_schedules_combined.xdmf`` (theta at t* on the four runs, on
step 1's finest box).
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
from fem_em_solver.thermal import solve_transient_conduction  # noqa: E402

from tests.validation.test_imp7_duty_cycle_fem import S0_BAND  # noqa: E402
from tests.validation.test_mat10_sar_to_heat_sphere import ROBIN_TAG, _robin_tags  # noqa: E402
from tests.validation.test_mat12_transient_sar_heating_sphere import E_BAND, RHO_C  # noqa: E402
from tests.validation.test_mat13_transient_perfusion_sphere import INTERIOR_BAND  # noqa: E402
from tests.validation.test_mat15_coil_driven_heating_phantom import (  # noqa: E402
    H_GEL,
    K_GEL,
    N_MAIN,
    T_INF,
)
from tests.validation.test_mat15_finite_cylinder_fem import _max_imag, _rel  # noqa: E402
from tests.validation.test_mat17_duty_cycle_fem import _census  # noqa: E402
from tests.validation.test_mat21_box_closed_form import (  # noqa: E402
    A_BOX,
    B_BOX,
    C_BOX,
    Q,
)
from tests.validation.test_mat21_box_dose_closed_form import DOSE_BAR, T_GRID  # noqa: E402
from tests.validation.test_mat21_box_dt_duty125_dose_closed_form import _q_check_g  # noqa: E402
from tests.validation.test_mat21_box_duty25_dose_closed_form import (  # noqa: E402
    DUTY25,
    SCHEDULE_25,
    T_SW25,
    _chi25,
    _sched_string,
)
from tests.validation.test_mat21_box_duty_closed_form import DT_BOX  # noqa: E402
from tests.validation.test_mat21_box_fem import (  # noqa: E402
    BASE,
    RUNGS,
    T_STAR,
    _exact,
    _interior_chi,
    _trans_ex,
)
from tests.validation.test_mat21_initial_gradient_closed_form import G_GRAD  # noqa: E402
from tests.validation.test_mat21_initial_gradient_fem import _grad_ex  # noqa: E402
from tests.validation.test_mat21_pre_equilibration_closed_form import (  # noqa: E402
    ADJ_PT,
    P5,
    _combined,
    _sourced,
)
from tests.validation.test_mat21_pre_equilibration_dose_closed_form import (  # noqa: E402
    _dose as _dose_d1,
)
from tests.validation.test_mat21_pre_equilibration_duty_dose_closed_form import (  # noqa: E402
    _dose as _dose_d25,
)
from tests.validation.test_mat21_pre_equilibration_duty_dose_closed_form import (  # noqa: E402
    _drifts,
    _uniform,
)
from tests.validation.test_mat21_pre_equilibration_fem import _read  # noqa: E402

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"  # committed, EX-57
BASENAME = "materials_36_gel_box_wait_across_schedules"

# Re-formed from imported pieces, as the gates' module-level lines.
HIST_LABELS = ("P5", "ADJ_PT")
HIST_PTS = np.array((P5, ADJ_PT), dtype=float)
T_FEM = (T_STAR / N_MAIN) * np.arange(N_MAIN + 1)
RUN_KEYS = ("d1_ref", "d1_tw0", "d25_ref", "d25_tw0")
SCHED_OF = {"d1": "1", "d25": f"{DUTY25:g}"}


def _series_d1(p, key):
    """Item 727's closed-form history at ``p`` on ``T_FEM``; sourced-only for "ref"."""
    if key == "ref":
        return np.array([_sourced(p, s) for s in T_FEM])
    return np.array([_combined(p, s, 0.0) for s in T_FEM])


def _series_d25(p, key):
    """Item 753's history at ``p``: the uniform duty sum, plus the drift on "tw0"."""
    sum_h = _uniform(p, T_SW25)
    return sum_h if key == "ref" else sum_h + _drifts(p, 0.0)


def _say(comm, msg):
    if comm.rank == 0:
        print(msg, flush=True)


def main() -> None:
    comm = MPI.COMM_WORLD
    started = time.perf_counter()

    if not np.issubdtype(np.dtype(default_scalar_type), np.complexfloating):
        raise RuntimeError(
            "this example expects the complex DolfinX build the runner's `mat:` group "
            "sources (source /usr/local/bin/dolfinx-complex-mode)"
        )

    _say(comm, "=" * 78)
    _say(comm, "EX-176 -- the same skipped wait under two schedules "
               "(MAT-21 steps 32-33 / P59, P60)")
    _say(comm, "=" * 78)
    _say(comm, f"[gate] imported: INTERIOR_BAND {INTERIOR_BAND:.0e}, E_BAND {E_BAND:.0e}, "
               f"S0_BAND {S0_BAND:.0e}, DOSE_BAR {DOSE_BAR:g}; (A, B, C) = "
               f"{(A_BOX, B_BOX, C_BOX)}, q = {Q} W/m^3, g = {G_GRAD} K/m, D = {DUTY25}, "
               f"schedule {SCHEDULE_25}")
    assert S0_BAND == 1e-12  # `_census` carries the same bar as a literal
    assert N_MAIN == 40 and RUNGS[-1] == 6 and T_GRID.size == N_MAIN + 1
    assert np.max(np.abs(T_GRID - T_FEM)) <= 1e-9 * DT_BOX
    assert _sched_string(_chi25) == SCHEDULE_25

    # ---- mesh: step 1's finest rung ----------------------------------------------------
    m = RUNGS[-1]
    ns = [b * m for b in BASE]
    msh = dmesh.create_box(comm, [np.array([-A_BOX, -B_BOX, -C_BOX]),
                                  np.array([A_BOX, B_BOX, C_BOX])],
                           ns, cell_type=dmesh.CellType.tetrahedron)
    ft = _robin_tags(msh)
    chi = _interior_chi(msh)
    n = msh.topology.index_map(3).size_global
    tdim = msh.topology.dim
    imap_c = msh.topology.index_map(tdim)
    n_all = imap_c.size_local + imap_c.num_ghosts
    cell_tags = dmesh.meshtags(
        msh, tdim, np.arange(n_all, dtype=np.int32), np.full(n_all, 1, dtype=np.int32))
    write_setup_figure(
        msh,
        cell_tags,
        FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names={1: "gel (650 x 420 x 90 mm), q on from t = 0 (D = 1) or duty 0.25, "
                         "theta_0 = 0 or g z"},
        title="mat:36 -- gel box, Robin on all six faces, one wait under two source "
              "schedules: y = 0 slice",
        slice_normal=(0.0, 1.0, 0.0),
        slice_origin=(0.0, 0.0, 0.0),
        comm=comm,
    )
    _say(comm, f"[mesh] rung m = {m} {tuple(ns)}: cells {n}; t* = {T_STAR:.3f} s, CN N = "
               f"{N_MAIN} + 2 BE half-steps, {N_MAIN + 1} snapshots per run; batch "
               f"{HIST_LABELS}; four runs (D = 1, D = {DUTY25} x reference, t_w = 0); ranks "
               f"{comm.size}")

    # ---- four CN runs with histories (the gates' `runs` bodies, written out) -----------
    kw = dict(k=K_GEL, rho_c=RHO_C, q=Q, h=H_GEL, t_inf=T_INF, t_init=T_INF,
              robin_tag=ROBIN_TAG, dt=T_STAR / N_MAIN, n_steps=N_MAIN, theta=0.5,
              startup_half_steps=2, degree=1, snapshot_steps=tuple(range(N_MAIN + 1)))
    starts = {"ref": None, "tw0": lambda x: T_INF + G_GRAD * x[2]}
    E, src, hist, hist_ok, last, res, sched = {}, {}, {}, {}, {}, {}, {}
    for sk, extra in (("d1", {}), ("d25", {"source_schedule": _chi25})):
        for key, expr in starts.items():
            rk = f"{sk}_{key}"
            tq = time.perf_counter()
            res[rk] = solve_transient_conduction(msh, ft, t_init_expr=expr, **kw, **extra)
            ts = time.perf_counter()
            E[rk] = float(np.real(res[rk].max_balance_residual))
            src[rk] = float(np.real(res[rk].source_power))
            if sk == "d25":
                sched[rk] = "".join(str(int(x)) for x in res[rk].source_multipliers)
                _census(res[rk], DT_BOX, 2, DUTY25 * T_STAR, f"mat36 d25 {key}")

            # One batch, one order, on every snapshot; read and clear.
            h = np.empty((N_MAIN + 1, HIST_PTS.shape[0]))
            ok_all = [True] * HIST_PTS.shape[0]
            for i in range(N_MAIN + 1):
                h[i], ok = _read(res[rk].snapshots[i], HIST_PTS)
                ok_all = [a and b for a, b in zip(ok_all, ok)]
            res[rk].snapshots.clear()
            hist[rk], hist_ok[rk] = np.real(h), ok_all
            last[rk] = np.real(_read(res[rk].temperature, HIST_PTS)[0])
            _say(comm, f"[run {rk}] solve {ts - tq:.1f} s, reads "
                       f"{time.perf_counter() - ts:.1f} s; (E) {E[rk]:.3e} (E_BAND "
                       f"{E_BAND:.0e}); source_power {src[rk]:.6e} W")

    # ---- the series on the mesh: one interpolation per switch, as the D = 0.25 gate ----
    tq = time.perf_counter()
    ex_c, npts = _exact(msh, _trans_ex(T_STAR))  # m = 0: also the continuous (D = 1) series
    ex_d = fem.Function(ex_c.function_space)
    ex_d.x.array[:] = ex_c.x.array
    for i, tm in enumerate(T_SW25):
        if i:
            g_m, _ = _exact(msh, _trans_ex(T_STAR - tm))
            ex_d.x.array[:] += (-1.0) ** i * g_m.x.array
    ex_0, _ = _exact(msh, _grad_ex(T_STAR))
    _say(comm, f"[series] {len(T_SW25) + 1} interpolations ({npts} P2 points): "
               f"{time.perf_counter() - tq:.1f} s")

    th = {rk: r.temperature - T_INF for rk, r in res.items()}
    imag = max(_max_imag(r.temperature) for r in res.values())
    I = {
        "d1_tw0": float(np.real(_rel(msh, th["d1_tw0"], ex_c + ex_0, chi))),
        "d25_ref": float(np.real(_rel(msh, th["d25_ref"], ex_d, chi))),
        "d25_tw0": float(np.real(_rel(msh, th["d25_tw0"], ex_d + ex_0, chi))),
    }
    I_cont = float(np.real(_rel(msh, th["d25_tw0"], ex_c + ex_0, chi)))
    for k, v in I.items():
        _say(comm, f"[I] {k}: interior rel L2 {v:.6e} (bound {INTERIOR_BAND:.0e}, ASSERTED)")
    _say(comm, f"[c-cont] continuous superposition on d25_tw0: {I_cont:.4e} (asserted > "
               f"{INTERIOR_BAND:.0e}; gate ceiling 7.4919e-01)")

    # ---- (S0): the recorded schedule on both D = 0.25 runs ------------------------------
    for rk, s in sched.items():
        _say(comm, f"[S0] {rk}: schedule {s} (must equal {SCHEDULE_25}); `_census` at "
                   f"{S0_BAND:.0e} passed in the run")

    # ---- (Q-wf): cem43 against the quadrature, eight histories --------------------------
    qwf = []
    for rk, h in hist.items():
        for j, lab in enumerate(HIST_LABELS):
            dose, quad, rel, _swap, miss = _q_check_g(T_INF + h[:, j], T_FEM)
            rel = float(np.real(rel))
            qwf.append((rk, lab, rel))
            _say(comm, f"[Q-wf] {rk} {lab}: D {np.real(dose):.6e} s, quad {np.real(quad):.6e} s, "
                       f"rel {rel:.3e} (bar {DOSE_BAR:g}, ASSERTED); R swap miss "
                       f"{float(np.real(miss)):.4e} (printed)")

    # ---- dose per run, point and schedule; (c-dose) at P5 -------------------------------
    dose_f = {"d1": lambda x: float(np.real(_dose_d1(x, T_FEM))),
              "d25": lambda x: float(np.real(_dose_d25(x)))}
    series_f = {"d1": _series_d1, "d25": _series_d25}
    d_fem, d_ser = {}, {}
    for sk in ("d1", "d25"):
        for j, (lab, p) in enumerate(zip(HIST_LABELS, HIST_PTS)):
            p = tuple(float(c) for c in p)
            for key in ("ref", "tw0"):
                d_fem[(sk, lab, key)] = dose_f[sk](hist[f"{sk}_{key}"][:, j])
                d_ser[(sk, lab, key)] = dose_f[sk](series_f[sk](p, key))
    miss_c = {sk: abs(d_fem[(sk, "P5", "tw0")] / d_fem[(sk, "P5", "ref")] - 1.0)
              for sk in ("d1", "d25")}
    for sk, v in miss_c.items():
        _say(comm, f"[c-dose] {sk} P5: t_w = 0 FEM D {d_fem[(sk, 'P5', 'tw0')]:.6e} s vs "
                   f"reference {d_fem[(sk, 'P5', 'ref')]:.6e} s, miss {v:.4e} (asserted > "
                   f"{DOSE_BAR:g})")

    # ---- printed: one wait, two schedules -----------------------------------------------
    _say(comm, "[TABLE] point | D | D_ref FEM s | D(t_w=0) FEM s | FEM D/D_ref - 1 | "
               "series D/D_ref - 1")
    for lab in HIST_LABELS:
        for sk in ("d1", "d25"):
            fr, f0 = d_fem[(sk, lab, "ref")], d_fem[(sk, lab, "tw0")]
            sr, s0 = d_ser[(sk, lab, "ref")], d_ser[(sk, lab, "tw0")]
            _say(comm, f"[TABLE] {lab} | {SCHED_OF[sk]} | {fr:.6e} | {f0:.6e} | "
                       f"{f0 / fr - 1.0:+.6e} | {s0 / sr - 1.0:+.6e}")
    _say(comm, "[note] a dose at D = 0.25 is smaller in absolute terms; the ratio compares "
               "one wait under two schedules on this fixture's g, one rung, one dt, and is "
               "not a rate")
    _say(comm, f"max|Im theta| = {imag:.3e}")

    # ---- output ------------------------------------------------------------------------
    V = res["d1_ref"].temperature.function_space
    outs = {}
    for rk in RUN_KEYS:
        f_th = fem.Function(V, name=f"theta_{rk}_t_star")
        f_th.x.array[:] = np.real(res[rk].temperature.x.array) - T_INF
        f_th.x.scatter_forward()
        outs[f_th.name] = f_th
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_combined", msh, cell_tags, outs, comm=comm)
    adopt_host_ownership(OUTPUT_DIR, comm=comm)
    _say(comm, f"\n[output] {path} -- theta(t*) on the four runs (D = 1 and D = {DUTY25}, "
               f"reference and t_w = 0), on the m = {m} box")

    # ---- assertions (the gates'), collected first --------------------------------------
    assert imag <= 1e-12, f"complex build left Im theta = {imag:.3e}"
    assert all(all(o) for o in hist_ok.values()), hist_ok
    for rk in hist:
        assert np.array_equal(hist[rk][N_MAIN], last[rk]), rk
    for k, v in I.items():
        assert v <= INTERIOR_BAND, f"(I) {k}: {v:.4e} > {INTERIOR_BAND}"
    assert src["d1_ref"] > 0.0 and src["d1_ref"] == src["d1_tw0"], src
    assert src["d25_ref"] > 0.0 and src["d25_ref"] == src["d25_tw0"], src
    for rk, e in E.items():
        assert e <= E_BAND, f"(E) {rk}: {e:.3e} > {E_BAND:.0e}"
    for rk, s in sched.items():
        assert s == SCHEDULE_25, (rk, s)
    for rk, lab, rel in qwf:
        assert rel <= DOSE_BAR, f"(Q-wf) {rk} {lab}: {rel:.3e} > {DOSE_BAR:g}"
    for sk, v in miss_c.items():
        assert v > DOSE_BAR, f"(c-dose) {sk}: miss {v:.4e} inside {DOSE_BAR:g}"
    assert I_cont > INTERIOR_BAND, f"(c-cont) {I_cont:.4e} inside {INTERIOR_BAND}"

    _say(comm, f"\n[timing] total {time.perf_counter() - started:.1f} s at -n {comm.size}")
    _say(comm, "\nAll imported gates hold.")


if __name__ == "__main__":
    main()
