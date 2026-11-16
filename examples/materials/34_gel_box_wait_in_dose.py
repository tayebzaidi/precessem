"""Example (`EX-171`): what skipping the F2182 wait costs in thermal dose -- the FEM
histories of the sourced gel box against the series' dose.

`mat:33` stops at temperature: the t_w = 0 start biases the rise at the top-face probe P5
by about 14 %. Dose is not linear in temperature (CEM43 grows as 4^(T - 43) above 43 C
and 2^(T - 43) below, integrated over the history), so the same leftover gradient reads
as a larger relative bias in dose. Here the three runs of `mat:33` are made again on the
same box (the finest rung of `MAT-21` step 1; **reference** from uniform T_inf; **t_w = 0**
from T_inf + g z; **t_w = T_C** from T_inf plus the drift series at T_C), each with the
per-step history at P0, P5 and ADJ_PT, and the CEM43 dose of each history is printed
beside the closed-form dose of `MAT-21` step 27 (P55).

**It asserts, and it does not re-implement** (the `ANS-1` rule). Bands, helpers, the
series and the step-1 record are imported by name from
``tests/validation/test_mat21_pre_equilibration_dose_fem.py``'s own imports (the modules
it draws on). That gate's ``runs`` fixture is a pytest fixture and is not called: its body
is written out here with every call it makes, plus ``np.real`` on each integral (the group
sources the complex build) and a ``write_setup_figure`` call after the mesh. The gate's
``HIST_PTS`` / ``T_FEM`` / ``_series`` are re-formed from imported pieces (they live in the
test module, and no ``test_*`` module is imported). No gate edit, no ``src/`` edit.

Asserted (bands imported; every gate value collected first, then asserted): (I) interior
rel L2 at t* against the superposition <= ``INTERIOR_BAND`` on both gradient runs; (E) the
balance <= ``E_BAND`` on all three, ``source_power`` > 0 and identical; (Sym') the
point-reflection pair identity on the drift part at P0-P4, <= ``SYM_BAR`` * g C; (Q-wf)
``cem43`` against the quadrature ``_q_check_g`` at ``DOSE_BAR`` on each of the nine
histories; (c-dose) the t_w = 0 FEM dose at P5 misses the reference run's by more than
``DOSE_BAR``; (c-ref) the reference reproduces ``STEP1_T`` within ``STEP1_T_BOUND``.
Printed only (records, no band): D / D_ref - 1 at P5 and ADJ_PT, FEM beside the series,
beside the temperature bias b at t*; the signed FEM-vs-series dose gap. No dose band, no
rate, no recommended wait.

Run it::

    ./run_examples.sh -e mat:34 -n 2 -t 300

Output: ``examples/materials/paraview_output/
materials_34_gel_box_wait_in_dose_combined.xdmf`` (theta at t* on the reference,
t_w = 0 and t_w = T_C runs, on step 1's finest box).
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
from tests.validation.test_mat21_box_closed_form import (  # noqa: E402
    A_BOX,
    B_BOX,
    C_BOX,
    Q,
    T_C,
)
from tests.validation.test_mat21_box_dose_closed_form import DOSE_BAR  # noqa: E402
from tests.validation.test_mat21_box_dt_duty125_dose_closed_form import _q_check_g  # noqa: E402
from tests.validation.test_mat21_box_fem import (  # noqa: E402
    BASE,
    RUNGS,
    T_STAR,
    _exact,
    _interior_chi,
    _trans_ex,
)
from tests.validation.test_mat21_box_probe_points_fem import PROBE_PTS  # noqa: E402
from tests.validation.test_mat21_initial_gradient_closed_form import G_GRAD, GC  # noqa: E402
from tests.validation.test_mat21_initial_gradient_fem import SYM_BAR, _grad_ex  # noqa: E402
from tests.validation.test_mat21_pre_equilibration_closed_form import (  # noqa: E402
    ADJ_PT,
    P5,
    _combined,
    _sourced,
)
from tests.validation.test_mat21_pre_equilibration_dose_closed_form import GRID, _dose  # noqa: E402
from tests.validation.test_mat21_pre_equilibration_fem import (  # noqa: E402
    STEP1_T,
    STEP1_T_BOUND,
    _read,
    _wait_init,
)

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"  # committed, EX-57
BASENAME = "materials_34_gel_box_wait_in_dose"

# Re-formed from imported pieces, as the gate's module-level lines (that module is a test
# module and is not imported).
HIST_LABELS = ("P0", "P5", "ADJ_PT")
HIST_PTS = np.array((PROBE_PTS[0], P5, ADJ_PT), dtype=float)
T_FEM = (T_STAR / N_MAIN) * np.arange(N_MAIN + 1)
WAIT = {"tw0": 0.0, "twC": T_C}


def _series(p, key, t):
    """Item 727's closed-form history at ``p`` on grid ``t``; sourced-only for "ref"."""
    if key == "ref":
        return np.array([_sourced(p, s) for s in t])
    return np.array([_combined(p, s, WAIT[key]) for s in t])


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
    _say(comm, "EX-171 -- what skipping the wait costs in dose: FEM histories against the "
               "series (MAT-21 steps 27-29 / P55, P56)")
    _say(comm, "=" * 78)
    _say(comm, f"[gate] imported: INTERIOR_BAND {INTERIOR_BAND:.0e}, E_BAND {E_BAND:.0e}, "
               f"SYM_BAR {SYM_BAR:.0e}, DOSE_BAR {DOSE_BAR:g}, STEP1_T_BOUND "
               f"{STEP1_T_BOUND:.0e}; (A, B, C) = {(A_BOX, B_BOX, C_BOX)}, q = {Q} W/m^3, "
               f"g = {G_GRAD} K/m, g C = {GC} K")

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
        region_names={1: "gel (650 x 420 x 90 mm), q on from t = 0, theta_0 = g z or the drift"},
        title="mat:34 -- gel box, Robin on all six faces, sourced from t = 0, "
              "dose read at P0, P5, ADJ_PT: y = 0 slice",
        slice_normal=(0.0, 1.0, 0.0),
        slice_origin=(0.0, 0.0, 0.0),
        comm=comm,
    )
    pts = np.array(PROBE_PTS, dtype=float)
    _say(comm, f"[mesh] rung m = {m} {tuple(ns)}: cells {n}; t_C = {T_C:.3f} s, "
               f"t* = {T_STAR:.3f} s, CN N = {N_MAIN} + 2 BE half-steps, {N_MAIN + 1} "
               f"snapshots per run; history batch {HIST_LABELS}; ranks {comm.size}")

    # ---- three CN runs with histories (the gate's `runs` body, written out) ------------
    kw = dict(k=K_GEL, rho_c=RHO_C, q=Q, h=H_GEL, t_inf=T_INF, t_init=T_INF,
              robin_tag=ROBIN_TAG, dt=T_STAR / N_MAIN, n_steps=N_MAIN, theta=0.5,
              startup_half_steps=2, degree=1, snapshot_steps=tuple(range(N_MAIN + 1)))
    starts = {"ref": None, "tw0": lambda x: T_INF + G_GRAD * x[2], "twC": _wait_init(T_C)}
    E, src, th_pt, ok_pt, hist, hist_ok, last, res = {}, {}, {}, {}, {}, {}, {}, {}
    for key, expr in starts.items():
        tq = time.perf_counter()
        res[key] = solve_transient_conduction(msh, ft, t_init_expr=expr, **kw)
        ts = time.perf_counter()
        E[key] = float(np.real(res[key].max_balance_residual))
        src[key] = float(np.real(res[key].source_power))
        th_pt[key], ok_pt[key] = _read(res[key].temperature, pts)

        # One batch, one order, on every snapshot; read and clear.
        h = np.empty((N_MAIN + 1, HIST_PTS.shape[0]))
        ok_all = [True] * HIST_PTS.shape[0]
        for i in range(N_MAIN + 1):
            h[i], ok = _read(res[key].snapshots[i], HIST_PTS)
            ok_all = [a and b for a, b in zip(ok_all, ok)]
        res[key].snapshots.clear()
        hist[key], hist_ok[key] = np.real(h), ok_all
        last[key] = np.real(_read(res[key].temperature, HIST_PTS)[0])
        _say(comm, f"[run {key}] solve {ts - tq:.1f} s, reads {time.perf_counter() - ts:.1f} s; "
                   f"(E) {E[key]:.3e} (E_BAND {E_BAND:.0e}); source_power "
                   f"{src[key]:.6e} W")

    # ---- the series on the mesh: three interpolations, as the gate ---------------------
    tq = time.perf_counter()
    ex_s, npts = _exact(msh, _trans_ex(T_STAR))
    ex_0, _ = _exact(msh, _grad_ex(T_STAR))
    ex_c, _ = _exact(msh, _grad_ex(T_STAR + T_C))
    _say(comm, f"[series] three P2 interpolations ({npts} points): "
               f"{time.perf_counter() - tq:.1f} s")

    th = {key: r.temperature - T_INF for key, r in res.items()}
    imag = max(_max_imag(r.temperature) for r in res.values())
    ref_T = float(np.real(_rel(msh, th["ref"], ex_s, chi)))
    I0 = float(np.real(_rel(msh, th["tw0"], ex_s + ex_0, chi)))
    IC = float(np.real(_rel(msh, th["twC"], ex_s + ex_c, chi)))
    _say(comm, f"[I] interior rel L2, t_w = 0: {I0:.6e}; t_w = T_C: {IC:.6e} "
               f"(bound {INTERIOR_BAND:.0e}; gate P52 1.125097e-03 / 9.076702e-04)")
    _say(comm, f"[c-ref] reference interior {ref_T:.6e} vs step 1 {STEP1_T:.6e}, |delta| "
               f"{abs(ref_T - STEP1_T):.3e} (bound {STEP1_T_BOUND:.0e})")

    # ---- (Sym') on the drift part ------------------------------------------------------
    worst = {}
    for key in ("tw0", "twC"):
        d = np.real(th_pt[key]) - np.real(th_pt["ref"])
        s12, s34, s0 = abs(d[1] + d[2]), abs(d[3] + d[4]), abs(d[0])
        worst[key] = max(s12, s34, s0)
        _say(comm, f"[Sym'] {key}: |P1 + P2|/(g C) {s12 / GC:.3e}, |P3 + P4|/(g C) "
                   f"{s34 / GC:.3e}, |P0|/(g C) {s0 / GC:.3e} (bar {SYM_BAR:.0e})")

    # ---- (Q-wf): cem43 against the quadrature, nine histories --------------------------
    qwf = []
    for key, h in hist.items():
        for j, lab in enumerate(HIST_LABELS):
            dose, quad, rel, _swap, miss = _q_check_g(T_INF + h[:, j], T_FEM)
            rel = float(np.real(rel))
            qwf.append((key, lab, rel))
            _say(comm, f"[Q-wf] {key} {lab}: D {np.real(dose):.6e} s, quad {np.real(quad):.6e} s, "
                       f"rel {rel:.3e} (bar {DOSE_BAR:g}, ASSERTED); R swap miss "
                       f"{float(np.real(miss)):.4e} (printed)")

    # ---- (c-dose): the t_w = 0 dose at P5 misses the reference's -----------------------
    j5 = HIST_LABELS.index("P5")
    d_0 = float(np.real(_dose(hist["tw0"][:, j5], T_FEM)))
    d_r = float(np.real(_dose(hist["ref"][:, j5], T_FEM)))
    miss_c = abs(d_0 / d_r - 1.0)
    _say(comm, f"[c-dose] P5: t_w = 0 FEM D {d_0:.6e} s vs reference {d_r:.6e} s, miss "
               f"{miss_c:.4e} (asserted > {DOSE_BAR:g}; gate record 3.4217e-01)")

    # ---- printed: dose bias, FEM beside the series, beside the temperature bias --------
    _say(comm, "[TABLE] point | run | D FEM s | D series 41 rows s | FEM/series41 - 1 | "
               "D series 161 rows s | FEM/series161 - 1 | FEM D/D_ref - 1 | series41 D/D_ref - 1 "
               "| series161 D/D_ref - 1 | b FEM (theta bias at t*) | b series")
    for j, (lab, p) in enumerate(zip(HIST_LABELS, HIST_PTS)):
        d_h = {key: float(np.real(_dose(hist[key][:, j], T_FEM))) for key in hist}
        d_41 = {key: float(np.real(_dose(_series(tuple(p), key, T_FEM), T_FEM))) for key in hist}
        d_161 = {key: float(np.real(_dose(_series(tuple(p), key, GRID)))) for key in hist}
        for key in hist:
            b_f = hist[key][N_MAIN, j] / hist["ref"][N_MAIN, j] - 1.0
            s_t = _series(tuple(p), key, T_FEM)[-1] / _series(tuple(p), "ref", T_FEM)[-1] - 1.0
            _say(comm, f"[TABLE] {lab} | {key} | {d_h[key]:.6e} | {d_41[key]:.6e} | "
                       f"{d_h[key] / d_41[key] - 1.0:+.3e} | {d_161[key]:.6e} | "
                       f"{d_h[key] / d_161[key] - 1.0:+.3e} | {d_h[key] / d_h['ref'] - 1.0:+.6e} | "
                       f"{d_41[key] / d_41['ref'] - 1.0:+.6e} | "
                       f"{d_161[key] / d_161['ref'] - 1.0:+.6e} | {b_f:+.6e} | {s_t:+.6e}")
    _say(comm, "[note] the FEM-vs-series dose gap at ADJ_PT changes sign relative to P5; "
               "MAT-21 step 30 (item 742) attributes that to P1 interpolation between mesh "
               "planes at this rung; this example does not re-diagnose it")
    _say(comm, f"max|Im theta| = {imag:.3e}")

    # ---- output ------------------------------------------------------------------------
    V = res["ref"].temperature.function_space
    outs = {}
    for key in ("ref", "tw0", "twC"):
        f_th = fem.Function(V, name=f"theta_{key}_t_star")
        f_th.x.array[:] = np.real(res[key].temperature.x.array) - T_INF
        f_th.x.scatter_forward()
        outs[f_th.name] = f_th
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_combined", msh, cell_tags, outs, comm=comm)
    adopt_host_ownership(OUTPUT_DIR, comm=comm)
    _say(comm, f"\n[output] {path} -- theta(t*) on the reference, t_w = 0 and t_w = T_C runs, "
               f"on the m = {m} box")

    # ---- assertions (the gate's), collected first --------------------------------------
    assert imag <= 1e-12, f"complex build left Im theta = {imag:.3e}"
    assert all(all(o) for o in hist_ok.values()), hist_ok
    for key in hist:
        assert np.array_equal(hist[key][N_MAIN], last[key]), key
    assert I0 <= INTERIOR_BAND, f"(I) t_w = 0: {I0:.4e} > {INTERIOR_BAND}"
    assert IC <= INTERIOR_BAND, f"(I) t_w = T_C: {IC:.4e} > {INTERIOR_BAND}"
    assert src["ref"] > 0.0 and src["ref"] == src["tw0"] == src["twC"], src
    for key, e in E.items():
        assert e <= E_BAND, f"(E) {key}: {e:.3e} > {E_BAND:.0e}"
    assert all(all(o) for o in ok_pt.values()), ok_pt
    for key, w in worst.items():
        assert w <= SYM_BAR * GC, f"(Sym') {key}: {w:.3e}"
    for key, lab, rel in qwf:
        assert rel <= DOSE_BAR, f"(Q-wf) {key} {lab}: {rel:.3e} > {DOSE_BAR:g}"
    assert miss_c > DOSE_BAR, f"(c-dose) miss {miss_c:.4e} inside {DOSE_BAR:g}"
    assert abs(ref_T - STEP1_T) <= STEP1_T_BOUND, f"(c-ref) {ref_T:.6e} vs {STEP1_T:.6e}"

    _say(comm, f"\n[timing] total {time.perf_counter() - started:.1f} s at -n {comm.size}")
    _say(comm, "\nAll imported gates hold.")


if __name__ == "__main__":
    main()
