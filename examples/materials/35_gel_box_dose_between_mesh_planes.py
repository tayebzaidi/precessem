"""Example (`EX-173`): a dose read between two mesh planes is the interpolant's -- the FEM
dose of the sourced gel box at a probe on a mesh plane, one plane away and midway.

`mat:32` and `mat:34` print FEM-vs-series dose gaps at probes without saying where a probe
sits against the mesh. Here the box is `MAT-21` step 1's finest rung (m = 6, 78 x 48 x 12
cells, dz = C / 6 along z), and two runs are made -- **reference** from uniform T_inf and
**t_w = 0** from T_inf + g z -- each with the per-step history read on one fixed batch
``(PLANE_PT, ADJ_PT, MID_PT)``: ``PLANE_PT`` on a mesh plane (z = 5 C / 6), ``ADJ_PT`` one
off-plane height, ``MID_PT`` midway between that plane and z = C. The CEM43 dose of each
history is printed beside the closed-form dose (41-row and 161-row grids) and beside the dose
of the series' own P1 interpolant between the two bracketing planes.

**It asserts, and it does not re-implement** (the `ANS-1` rule). Bands, constants, helpers
and the series are imported by name from the step-30 gate
``tests/validation/test_mat21_pre_equilibration_dose_plane.py`` (``PLANE_PT``, ``MID_PT``,
``TOP_PT``, ``ITEM736_ADJ_D``, ``REPRO_BOUND``, ``PLANE_BOUND``, ``DZ``, ``_at_t``,
``_box``) and the modules it draws on. That gate's ``runs`` fixture is a pytest fixture and is
not called: its body is written out here, plus ``np.real`` on each integral (the group sources
the complex build) and a ``write_setup_figure`` call after the mesh. ``_series`` and ``T_FEM``
are re-formed from imported pieces, as in `mat:34`. No gate edit, no ``src/`` edit.

Asserted (bands imported; every gate value collected first, then asserted): (I) interior rel
L2 at t* against the superposition <= ``INTERIOR_BAND`` on t_w = 0; (E) the balance <=
``E_BAND`` on both runs, ``source_power`` > 0 and identical; (Q-wf) ``cem43`` against the
quadrature ``_q_check_g`` at ``DOSE_BAR`` on all six histories; (Grid) ``PLANE_PT``'s z on a
mesh plane and ``MID_PT``'s midway, each to ``PLANE_BOUND`` (relative to C); (c-rep) ``ADJ_PT``'s
FEM doses reproduce item 736's at ``REPRO_BOUND``; (c-plane) the t_w = 0 dose at ``PLANE_PT``
misses the reference's by more than ``DOSE_BAR``. Printed only (records, no band): the signed
gap FEM / series - 1 at the three heights on both runs, beside the interpolant's own gap and
the gap less the plane's. No dose band, no rate, no recommended probe height.

Run it::

    ./run_examples.sh -e mat:35 -n 2 -t 300

Output: ``examples/materials/paraview_output/
materials_35_gel_box_dose_between_mesh_planes_combined.xdmf`` (theta at t* on the reference
and t_w = 0 runs, on step 1's finest box).
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
from tests.validation.test_mat21_initial_gradient_closed_form import G_GRAD  # noqa: E402
from tests.validation.test_mat21_initial_gradient_fem import _grad_ex  # noqa: E402
from tests.validation.test_mat21_pre_equilibration_closed_form import (  # noqa: E402
    ADJ_PT,
    _combined,
    _sourced,
)
from tests.validation.test_mat21_pre_equilibration_dose_closed_form import GRID, _dose  # noqa: E402
from tests.validation.test_mat21_pre_equilibration_dose_plane import (  # noqa: E402
    DZ,
    ITEM736_ADJ_D,
    MID_PT,
    PLANE_BOUND,
    PLANE_PT,
    REPRO_BOUND,
    TOP_PT,
)
from tests.validation.test_mat21_pre_equilibration_fem import _read  # noqa: E402

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"  # committed, EX-57
BASENAME = "materials_35_gel_box_dose_between_mesh_planes"

# Re-formed from imported pieces, as the gate's module-level lines.
HIST_LABELS = ("PLANE_PT", "ADJ_PT", "MID_PT")
HIST_PTS = np.array((PLANE_PT, ADJ_PT, MID_PT), dtype=float)
KEYS = ("ref", "tw0")
T_FEM = (T_STAR / N_MAIN) * np.arange(N_MAIN + 1)


def _series(p, key, t):
    """Item 727's closed-form history at ``p`` on grid ``t``; sourced-only for "ref"."""
    if key == "ref":
        return np.array([_sourced(p, s) for s in t])
    return np.array([_combined(p, s, 0.0) for s in t])


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
    _say(comm, "EX-173 -- a dose read between two mesh planes is the interpolant's "
               "(MAT-21 step 30 / P57)")
    _say(comm, "=" * 78)
    _say(comm, f"[gate] imported: INTERIOR_BAND {INTERIOR_BAND:.0e}, E_BAND {E_BAND:.0e}, "
               f"DOSE_BAR {DOSE_BAR:g}, REPRO_BOUND {REPRO_BOUND:g}, PLANE_BOUND "
               f"{PLANE_BOUND:g}; (A, B, C) = {(A_BOX, B_BOX, C_BOX)}, q = {Q} W/m^3, "
               f"g = {G_GRAD} K/m, T_C = {T_C:.3f} s")

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
        region_names={1: "gel (650 x 420 x 90 mm), q on from t = 0, theta_0 = 0 or g z"},
        title="mat:35 -- gel box, Robin on all six faces, dose read on a mesh plane, "
              "one plane away and midway: y = 0 slice",
        slice_normal=(0.0, 1.0, 0.0),
        slice_origin=(0.0, 0.0, 0.0),
        comm=comm,
    )
    _say(comm, f"[mesh] rung m = {m} {tuple(ns)}: cells {n}; t* = {T_STAR:.3f} s, CN N = "
               f"{N_MAIN} + 2 BE half-steps, {N_MAIN + 1} snapshots per run; batch "
               f"{HIST_LABELS}; ranks {comm.size}")

    # ---- (Grid): the mesh's own z coordinates, on every rank's vertices ----------------
    z = np.unique(np.concatenate(comm.allgather(np.unique(msh.geometry.x[:, 2]))))
    d_plane = float(np.min(np.abs(z - PLANE_PT[2]))) / C_BOX
    lo, hi = float(z[z <= MID_PT[2]].max()), float(z[z >= MID_PT[2]].min())
    d_mid = abs(MID_PT[2] - 0.5 * (lo + hi)) / C_BOX
    d_dz = abs((hi - lo) / DZ - 1.0)
    grid_adj = lo < ADJ_PT[2] < hi and abs(lo - PLANE_PT[2]) / C_BOX <= PLANE_BOUND
    grid_top = abs(hi / C_BOX - 1.0) <= PLANE_BOUND
    _say(comm, f"[Grid] nz {ns[2]}: PLANE_PT off a plane {d_plane:.3e} C; MID_PT off midway "
               f"{d_mid:.3e} C; dz / (C / 6) - 1 {d_dz:.3e} (bound {PLANE_BOUND:g}); "
               f"bracketing planes z = {lo:.6e}, {hi:.6e} m")

    # ---- two CN runs with histories (the gate's `runs` body, written out) --------------
    kw = dict(k=K_GEL, rho_c=RHO_C, q=Q, h=H_GEL, t_inf=T_INF, t_init=T_INF,
              robin_tag=ROBIN_TAG, dt=T_STAR / N_MAIN, n_steps=N_MAIN, theta=0.5,
              startup_half_steps=2, degree=1, snapshot_steps=tuple(range(N_MAIN + 1)))
    starts = {"ref": None, "tw0": lambda x: T_INF + G_GRAD * x[2]}
    E, src, hist, hist_ok, last, res = {}, {}, {}, {}, {}, {}
    for key, expr in starts.items():
        tq = time.perf_counter()
        res[key] = solve_transient_conduction(msh, ft, t_init_expr=expr, **kw)
        ts = time.perf_counter()
        E[key] = float(np.real(res[key].max_balance_residual))
        src[key] = float(np.real(res[key].source_power))

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

    # ---- the series on the mesh and (I) on t_w = 0 --------------------------------------
    tq = time.perf_counter()
    ex_s, npts = _exact(msh, _trans_ex(T_STAR))
    ex_0, _ = _exact(msh, _grad_ex(T_STAR))
    _say(comm, f"[series] two P2 interpolations ({npts} points): "
               f"{time.perf_counter() - tq:.1f} s")
    imag = max(_max_imag(r.temperature) for r in res.values())
    I0 = float(np.real(_rel(msh, res["tw0"].temperature - T_INF, ex_s + ex_0, chi)))
    _say(comm, f"[I] interior rel L2, t_w = 0: {I0:.6e} (bound {INTERIOR_BAND:.0e}; gate "
               f"P52 1.125097e-03)")

    # ---- (Q-wf): cem43 against the quadrature, six histories ----------------------------
    qwf = []
    for key, h in hist.items():
        for j, lab in enumerate(HIST_LABELS):
            dose, quad, rel, _swap, miss = _q_check_g(T_INF + h[:, j], T_FEM)
            rel = float(np.real(rel))
            qwf.append((key, lab, rel))
            _say(comm, f"[Q-wf] {key} {lab}: D {np.real(dose):.6e} s, quad {np.real(quad):.6e} s, "
                       f"rel {rel:.3e} (bar {DOSE_BAR:g}, ASSERTED); R swap miss "
                       f"{float(np.real(miss)):.4e} (printed)")

    # ---- (c-rep): ADJ_PT's FEM doses reproduce item 736's -------------------------------
    jadj = HIST_LABELS.index("ADJ_PT")
    c_rep = {}
    for key, rec in ITEM736_ADJ_D.items():
        d = float(np.real(_dose(hist[key][:, jadj], T_FEM)))
        c_rep[key] = abs(d / rec - 1.0)
        _say(comm, f"[c-rep] ADJ_PT {key}: D {d:.6e} s vs item 736 {rec:.6e} s, rel "
                   f"{c_rep[key]:.3e} (bound {REPRO_BOUND:g}, ASSERTED)")

    # ---- (c-plane): the t_w = 0 dose at PLANE_PT misses the reference's -----------------
    jp = HIST_LABELS.index("PLANE_PT")
    d_0 = float(np.real(_dose(hist["tw0"][:, jp], T_FEM)))
    d_r = float(np.real(_dose(hist["ref"][:, jp], T_FEM)))
    miss_c = abs(d_0 / d_r - 1.0)
    _say(comm, f"[c-plane] PLANE_PT: t_w = 0 FEM D {d_0:.6e} s vs reference {d_r:.6e} s, miss "
               f"{miss_c:.4e} (asserted > {DOSE_BAR:g}; gate record 3.5830e-01)")

    # ---- printed: the signed gap at three heights, beside the interpolant's -------------
    _say(comm, "[TABLE] height | run | D FEM s | D series161 s | FEM/series161-1 | D series41 s "
               "| FEM/series41-1 | P1 interpolant/series41-1 | FEM/interpolant-1 | "
               "gap41 minus plane gap41")
    for key in KEYS:
        s_lo = _series(PLANE_PT, key, T_FEM)
        s_hi = _series(TOP_PT, key, T_FEM)
        plane_gap = None
        for j, (lab, p) in enumerate(zip(HIST_LABELS, HIST_PTS)):
            p = tuple(float(c) for c in p)
            fem_h = hist[key][:, j]
            s41, s161 = _series(p, key, T_FEM), _series(p, key, GRID)
            d_h = float(np.real(_dose(fem_h, T_FEM)))
            d_41 = float(np.real(_dose(s41, T_FEM)))
            d_161 = float(np.real(_dose(s161)))
            gap41 = d_h / d_41 - 1.0
            if lab == "PLANE_PT":
                plane_gap = gap41
                interp, d_int = 0.0, d_41
            else:
                a = (p[2] - PLANE_PT[2]) / DZ
                d_int = float(np.real(_dose((1.0 - a) * s_lo + a * s_hi, T_FEM)))
                interp = d_int / d_41 - 1.0
            _say(comm, f"[TABLE] {lab} | {key} | {d_h:.6e} | {d_161:.6e} | "
                       f"{d_h / d_161 - 1.0:+.3e} | {d_41:.6e} | {gap41:+.3e} | "
                       f"{interp:+.3e} | {d_h / d_int - 1.0:+.3e} | {gap41 - plane_gap:+.3e}")
    _say(comm, "[note] a dose at a point is cem43 on a history read through P1; a sign "
               "there is not by itself a solver error; what sets the on-plane gap is not "
               "addressed here")
    _say(comm, f"max|Im theta| = {imag:.3e}")

    # ---- output ------------------------------------------------------------------------
    V = res["ref"].temperature.function_space
    outs = {}
    for key in KEYS:
        f_th = fem.Function(V, name=f"theta_{key}_t_star")
        f_th.x.array[:] = np.real(res[key].temperature.x.array) - T_INF
        f_th.x.scatter_forward()
        outs[f_th.name] = f_th
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_combined", msh, cell_tags, outs, comm=comm)
    adopt_host_ownership(OUTPUT_DIR, comm=comm)
    _say(comm, f"\n[output] {path} -- theta(t*) on the reference and t_w = 0 runs, on the "
               f"m = {m} box")

    # ---- assertions (the gate's), collected first --------------------------------------
    assert imag <= 1e-12, f"complex build left Im theta = {imag:.3e}"
    assert all(all(o) for o in hist_ok.values()), hist_ok
    for key in hist:
        assert np.array_equal(hist[key][N_MAIN], last[key]), key
    assert d_dz <= PLANE_BOUND, (lo, hi)
    assert d_plane <= PLANE_BOUND, d_plane
    assert d_mid <= PLANE_BOUND, d_mid
    assert grid_adj and grid_top, (grid_adj, grid_top)
    assert I0 <= INTERIOR_BAND, f"(I) t_w = 0: {I0:.4e} > {INTERIOR_BAND}"
    assert src["ref"] > 0.0 and src["ref"] == src["tw0"], src
    for key, e in E.items():
        assert e <= E_BAND, f"(E) {key}: {e:.3e} > {E_BAND:.0e}"
    for key, lab, rel in qwf:
        assert rel <= DOSE_BAR, f"(Q-wf) {key} {lab}: {rel:.3e} > {DOSE_BAR:g}"
    for key, rel in c_rep.items():
        assert rel <= REPRO_BOUND, f"(c-rep) {key}: {rel:.3e} > {REPRO_BOUND:g}"
    assert miss_c > DOSE_BAR, f"(c-plane) miss {miss_c:.4e} inside {DOSE_BAR:g}"

    _say(comm, f"\n[timing] total {time.perf_counter() - started:.1f} s at -n {comm.size}")
    _say(comm, "\nAll imported gates hold.")


if __name__ == "__main__":
    main()
