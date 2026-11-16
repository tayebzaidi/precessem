"""Example (`EX-169`): what skipping the F2182 wait costs -- the gel box's initial
gradient biasing the sourced rise, against the triple Robin series.

`mat:29` shows the leftover gradient relaxing with the source off. No example switches
a source on after a wait. Here the source is on (q from t = 0) and three runs are made
on the same box, the finest rung of `MAT-21` step 1: the **reference** from uniform
T_inf; **t_w = 0**, from T_inf + g z (the wait skipped); **t_w = T_C**, from T_inf plus
the drift series evaluated at t = T_C on the mesh (the wait kept). By linearity the
t_w run is the sourced rise from uniform plus the drift advanced by t_w,

    theta(t; t_w) = theta_drift(t + t_w) + theta_sourced(t),

and the closed form of `MAT-21` step 24 (P51) is the sum of the two series.

**It asserts, and it does not re-implement** (the `ANS-1` rule). The bands, the
series wrappers, the probe points, the wait-initial helper and the step-1 record are
imported by name from ``tests/validation/test_mat21_pre_equilibration_fem.py`` and the
modules that gate imports from. That gate's ``runs`` fixture is a pytest fixture and is
not called: its body is written out here with every call it makes, plus ``np.real`` on
each integral (the group sources the complex build). No ``test_*`` name is imported. No
gate edit, no ``src/`` edit.

Asserted (bands imported; every gate value collected first, then asserted): (I)
interior rel L2 at t* against the superposition <= ``INTERIOR_BAND`` on both gradient
runs; (E) the balance <= ``E_BAND`` on all three runs, with ``source_power`` > 0 and
identical on all three; (Sym') the point-reflection pair identity on the drift part
(run minus reference) at P0-P4, <= ``SYM_BAR`` * g C; (c-ref) the reference reproduces
step 1's interior record within ``STEP1_T_BOUND``; (c-Dz) the Dirichlet-z series misses
the t_w = 0 field by more than ``INTERIOR_BAND``.
Printed only (records, no band): the bias at P37's six points, FEM beside P51, at
t_w = 0 and T_C; P51's closed-form column at t_w = 0, 0.5, 1, 2 / lambda_1. (Dec),
P51's parked gate, is neither asserted nor printed as a pass.

Run it::

    ./run_examples.sh -e mat:33 -n 2 -t 300

Output: ``examples/materials/paraview_output/
materials_33_gel_box_pre_equilibration_bias_combined.xdmf`` (theta at t* at t_w = 0 and
t_w = T_C, the superpositions, the errors, on step 1's finest box).
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
from fem_em_solver.post.evaluation import evaluate_vector_field_parallel  # noqa: E402
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
from tests.validation.test_mat21_box_fem import (  # noqa: E402
    BASE,
    RUNGS,
    T_STAR,
    _exact,
    _interior_chi,
    _trans_ex,
)
from tests.validation.test_mat21_box_probe_points_fem import (  # noqa: E402
    PROBE_LABELS,
    PROBE_PTS,
)
from tests.validation.test_mat21_initial_gradient_closed_form import G_GRAD, GC  # noqa: E402
from tests.validation.test_mat21_initial_gradient_fem import (  # noqa: E402
    DZ_ROOTS,
    SYM_BAR,
    _grad_ex,
)
from tests.validation.test_mat21_pre_equilibration_closed_form import (  # noqa: E402
    LAM_ODD,
    WAITS,
    _drift,
    _sourced,
)
from tests.validation.test_mat21_pre_equilibration_fem import (  # noqa: E402
    STEP1_T,
    STEP1_T_BOUND,
    _read,
    _wait_init,
)

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"  # committed, EX-57
BASENAME = "materials_33_gel_box_pre_equilibration_bias"


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
    _say(comm, "EX-169 -- what skipping the wait costs: the sourced gel box from g z "
               "(MAT-21 steps 24, 26 / P51, P52)")
    _say(comm, "=" * 78)
    _say(comm, f"[gate] imported: INTERIOR_BAND {INTERIOR_BAND:.0e}, E_BAND {E_BAND:.0e}, "
               f"SYM_BAR {SYM_BAR:.0e}, STEP1_T_BOUND {STEP1_T_BOUND:.0e}; (A, B, C) = "
               f"{(A_BOX, B_BOX, C_BOX)}, q = {Q} W/m^3, g = {G_GRAD} K/m, g C = {GC} K")

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
        title="mat:33 -- gel box, Robin on all six faces, sourced from t = 0, "
              "T(0) = T_inf + g z (or waited): y = 0 slice",
        slice_normal=(0.0, 1.0, 0.0),
        slice_origin=(0.0, 0.0, 0.0),
        comm=comm,
    )
    pts = np.array(PROBE_PTS, dtype=float)
    _say(comm, f"[mesh] rung m = {m} {tuple(ns)}: cells {n}; t_C = {T_C:.3f} s, "
               f"t* = {T_STAR:.3f} s, CN N = {N_MAIN} + 2 BE half-steps; ranks {comm.size}")

    # ---- three CN runs (the gate's `runs` body, written out) ---------------------------
    kw = dict(k=K_GEL, rho_c=RHO_C, q=Q, h=H_GEL, t_inf=T_INF, t_init=T_INF,
              robin_tag=ROBIN_TAG, dt=T_STAR / N_MAIN, n_steps=N_MAIN, theta=0.5,
              startup_half_steps=2, degree=1)
    starts = {"ref": None, "tw0": lambda x: T_INF + G_GRAD * x[2], "twC": _wait_init(T_C)}
    E, src, th_pt, ok_pt, res = {}, {}, {}, {}, {}
    for key, expr in starts.items():
        tq = time.perf_counter()
        res[key] = solve_transient_conduction(msh, ft, t_init_expr=expr, **kw)
        E[key] = float(np.real(res[key].max_balance_residual))
        src[key] = float(np.real(res[key].source_power))
        th_pt[key], ok_pt[key] = _read(res[key].temperature, pts)
        _say(comm, f"[run {key}] {time.perf_counter() - tq:.1f} s; (E) {E[key]:.3e} "
                   f"(E_BAND {E_BAND:.0e}); source_power {src[key]:.9e} W, / q V_box "
                   f"{src[key] / (Q * 8.0 * A_BOX * B_BOX * C_BOX) - 1.0:+.3e}")

    # ---- the series on the mesh: four interpolations, as the gate (no fifth) -----------
    tq = time.perf_counter()
    ex_s, npts = _exact(msh, _trans_ex(T_STAR))
    ex_0, _ = _exact(msh, _grad_ex(T_STAR))
    ex_c, _ = _exact(msh, _grad_ex(T_STAR + T_C))
    ex_dz, _ = _exact(msh, _grad_ex(T_STAR, roots_z=DZ_ROOTS))
    _say(comm, f"[series] four P2 interpolations ({npts} points): {time.perf_counter() - tq:.1f} s")

    th = {key: r.temperature - T_INF for key, r in res.items()}
    imag = max(_max_imag(r.temperature) for r in res.values())
    ref_T = _rel(msh, th["ref"], ex_s, chi)
    I0 = _rel(msh, th["tw0"], ex_s + ex_0, chi)
    IC = _rel(msh, th["twC"], ex_s + ex_c, chi)
    Dz = _rel(msh, th["tw0"], ex_s + ex_dz, chi)
    whole0 = _rel(msh, th["tw0"], ex_s + ex_0)
    wholeC = _rel(msh, th["twC"], ex_s + ex_c)
    _say(comm, f"[I] interior rel L2, t_w = 0: {I0:.6e}; t_w = T_C: {IC:.6e} "
               f"(bound {INTERIOR_BAND:.0e}); whole {whole0:.6e} / {wholeC:.6e} (printed)")
    _say(comm, f"[c-ref] reference interior {ref_T:.6e} vs step 1 {STEP1_T:.6e}, |delta| "
               f"{abs(ref_T - STEP1_T):.3e} (bound {STEP1_T_BOUND:.0e}); [c-Dz] Dirichlet-z "
               f"interior miss {Dz:.4e} (asserted > {INTERIOR_BAND:.0e}; gate record 6.76e-02)")

    # ---- (Sym') on the drift part ------------------------------------------------------
    worst = {}
    for key in ("tw0", "twC"):
        d = th_pt[key] - th_pt["ref"]
        s12, s34, s0 = abs(d[1] + d[2]), abs(d[3] + d[4]), abs(d[0])
        worst[key] = max(s12, s34, s0)
        _say(comm, f"[Sym'] {key}: |P1 + P2|/(g C) {s12 / GC:.3e}, |P3 + P4|/(g C) "
                   f"{s34 / GC:.3e}, |P0|/(g C) {s0 / GC:.3e} (bar {SYM_BAR:.0e})")

    # ---- printed: the bias table, FEM beside P51 ---------------------------------------
    _say(comm, "[TABLE] point | t_w | FEM drift K | P51 drift K | FEM sourced K | P51 sourced K "
               "| FEM bias | P51 bias | |FEM drift|/(g C)")
    ref = th_pt["ref"]
    for key, t_w, tag in (("tw0", 0.0, "0"), ("twC", T_C, "T_C")):
        for i, (lab, p) in enumerate(zip(PROBE_LABELS, PROBE_PTS)):
            d_h = th_pt[key][i] - ref[i]
            d_s, s_s = _drift(p, T_STAR + t_w), _sourced(p, T_STAR)
            _say(comm, f"[TABLE] {lab} | {tag} | {d_h:+.6e} | {d_s:+.6e} | {ref[i]:.6f} | "
                       f"{s_s:.6f} | {d_h / ref[i]:+.6e} | {d_s / s_s:+.6e} | "
                       f"{abs(d_h) / GC:.3e}")
    p5 = PROBE_PTS[5]
    s5 = _sourced(p5, T_STAR)
    _say(comm, f"[P51] lambda_1 (odd) {LAM_ODD:.6e} /s, 1/lambda_1 {1.0 / LAM_ODD:.1f} s; "
               f"P5 closed-form bias against the wait:")
    for c, w in zip((0.0, 0.5, 1.0, 2.0), WAITS):
        d = _drift(p5, T_STAR + w)
        _say(comm, f"[P51] t_w = {c} / lambda_1 = {w:.1f} s: drift {d:+.6e} K, sourced "
                   f"{s5:.6f} K, bias {d / s5:+.6e}")
    _say(comm, f"max|Im theta| = {imag:.3e}")

    # ---- output ------------------------------------------------------------------------
    V = res["ref"].temperature.function_space
    # P2 series -> P1 output by Function interpolation (no extra series evaluation).
    def to_p1(name, f2):
        f = fem.Function(V, name=name)
        f.interpolate(f2)
        return f

    s1, z1, c1 = to_p1("s", ex_s), to_p1("z", ex_0), to_p1("c", ex_c)
    outs = {}
    for key, tag, drift1 in (("tw0", "tw0", z1), ("twC", "twC", c1)):
        f_th = fem.Function(V, name=f"theta_{tag}_t_star")
        f_th.x.array[:] = np.real(res[key].temperature.x.array) - T_INF
        f_sup = fem.Function(V, name=f"theta_superposition_{tag}")
        f_sup.x.array[:] = s1.x.array + drift1.x.array
        f_err = fem.Function(V, name=f"theta_error_{tag}")
        f_err.x.array[:] = f_th.x.array - f_sup.x.array
        for f in (f_th, f_sup, f_err):
            f.x.scatter_forward()
            outs[f.name] = f
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_combined", msh, cell_tags, outs, comm=comm)
    adopt_host_ownership(OUTPUT_DIR, comm=comm)
    _say(comm, f"\n[output] {path} -- theta(t*) at t_w = 0 and T_C, superpositions, errors "
               f"on the m = {m} box")

    # ---- assertions (the gate's), collected first --------------------------------------
    assert imag <= 1e-12, f"complex build left Im theta = {imag:.3e}"
    assert I0 <= INTERIOR_BAND, f"(I) t_w = 0: {I0:.4e} > {INTERIOR_BAND}"
    assert IC <= INTERIOR_BAND, f"(I) t_w = T_C: {IC:.4e} > {INTERIOR_BAND}"
    assert src["ref"] > 0.0 and src["ref"] == src["tw0"] == src["twC"], src
    for key, e in E.items():
        assert e <= E_BAND, f"(E) {key}: {e:.3e} > {E_BAND:.0e}"
    assert all(all(o) for o in ok_pt.values()), ok_pt
    for key, w in worst.items():
        assert w <= SYM_BAR * GC, f"(Sym') {key}: {w:.3e}"
    assert abs(ref_T - STEP1_T) <= STEP1_T_BOUND, f"(c-ref) {ref_T:.6e} vs {STEP1_T:.6e}"
    assert Dz > INTERIOR_BAND, f"(c-Dz) {Dz:.4e} inside band"

    _say(comm, f"\n[timing] total {time.perf_counter() - started:.1f} s at -n {comm.size}")
    _say(comm, "\nAll imported gates hold.")


if __name__ == "__main__":
    main()
