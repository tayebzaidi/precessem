"""Example (`EX-168`): the dose around the implant hole against the duty cycle.

`imp:18` reads the CEM43 dose and the six-minute window at one duty, D = 0.5.  This
example asks what the *drive* does: the same hole, the same three-fold gel shell, and
the same mean power delivered at different duties.  It narrates two gate modules:

1. `IMP-7` step 11 (I50, `tests/validation/test_imp7_duty25_dose_fem.py`, item 719): the
   Crank-Nicolson FEM history of the Robin-cooled gel shell on curved `LADDER[0]` (P2)
   under the D = 0.25 source (on-phase of 5 dt), 161 snapshots read at the pole
   (1.05 a) and the mid-shell point on the z axis in one fixed batch, against the
   hollow-sphere D = 0.25 series;
2. `IMP-7` step 10 (I49, `tests/validation/test_imp7_duty_matrix_closed_form.py`, item
   711): the closed-form matrix, D = 1, 0.5, 0.25 at both points, at the peak-power
   scale and the equal-mean-power scale (theta / D).

**It asserts, and it does not re-implement** (the `ANS-1` rule).  Every band, grid and
helper is imported by name; the `fem_dose` fixture and every `test_*` name are never
imported, so its body is written out below (disclosed in the guide).  Disclosed local
pieces: `_err_c` (``ufl.inner(d, d)`` with the real part taken, as `imp:17` / `imp:18`
do); every integral, norm and series value is passed through ``np.real`` in this
complex build.

Asserted: (E) balance <= E_BAND; (S0) census to the on-time <= S0_BAND; every point
valid; max T < T_REF; (T) interior rel L2 at 4 t* <= INTERIOR_BAND; (Q-hd) and (Win-hd) at
DOSE_BAR / AVG_BAR at both points and both scales; (Q-m) and (Win-m) on the matrix's 12
cells; negative controls (c-R) the R swap, (c-cont) the continuous series, (c-50) the
D = 0.5 series, each missing.

Printed only: the matrix (D, point, D(peak), D(equal mean), six-minute peak, cap) and the
FEM-vs-series gap at D = 0.25 (no band exists).  The doses are tiny (about 1e-11 s at the
fixture's source scale): an identity on a history, not a thermal-damage statement.  The
fixture has no metal core, no coil, no tissue.

Run it::

    ./scripts/run_examples.sh -e imp:19 -n 4 -t 480

Output: ``examples/implants/paraview_output/implants_19_implant_dose_against_duty_combined.xdmf``.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import ufl
from dolfinx import default_scalar_type, fem
from mpi4py import MPI

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from fem_em_solver.io.paraview_utils import (  # noqa: E402
    adopt_host_ownership,
    write_xdmf_with_tags,
)
from fem_em_solver.post.setup_figure import write_setup_figure  # noqa: E402
from fem_em_solver.thermal.dose import cem43, running_average  # noqa: E402

from tests.validation import test_imp7_implant_heating_fem as m7  # noqa: E402
from tests.validation.test_imp7_duty_cycle_fem import DUTY, PERIOD, S0_BAND, _chi  # noqa: E402,F401
from tests.validation.test_imp7_duty_dose_closed_form import (  # noqa: E402
    AVG_BAR,
    DOSE_BAR,
    DT,
    SCALES,  # noqa: F401  (labels only: its equal-mean factor is bound to D = 0.5)
    T_END,
    T_GRID,
)
from tests.validation.test_imp7_duty_dose_fem import (  # noqa: E402
    N_RUN,
    STARTUP,
    _probe,
    _shell_fixture,
)
from tests.validation.test_imp7_duty_matrix_closed_form import (  # noqa: E402
    DUTIES,
    M_LABELS,
    M_PTS,
    _duty_hist_at,
    _scales,
)
from tests.validation.test_imp7_duty25_dose_fem import (  # noqa: E402
    CAP_50_POLE,
    D25,
    ITEM_700_POLE,
    ON_END_25,
    _census_25,
    _duty_fn_at,
)
from tests.validation.test_imp7_implant_heating_closed_form import T_INF  # noqa: E402
from tests.validation.test_imp7_transient_closed_form import _modes  # noqa: E402
from tests.validation.test_imp7_transient_fem import (  # noqa: E402
    LADDER,
    RT,
    _series_fn,
    _transient,
)
from tests.validation.test_mat12_transient_sar_heating_sphere import E_BAND  # noqa: E402
from tests.validation.test_mat13_transient_perfusion_sphere import INTERIOR_BAND  # noqa: E402
from tests.validation.test_mat20_running_average import W_B2  # noqa: E402
from tests.validation.test_mat20_tuned_heating_dose import T_REF  # noqa: E402
from tests.validation.test_mat21_wall_duty_dose_closed_form import _q_check, _win_check  # noqa: E402

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"  # committed, EX-57
BASENAME = "implants_19_implant_dose_against_duty"
RES = LADDER[0]


def _say(comm, msg):
    if comm.rank == 0:
        print(msg, flush=True)


def _err_c(sub, T_h, theta):
    """(whole, interior |x|^2 < 0.81 b^2) rel L2 on theta; `inner`, real part (complex build)."""
    comm = sub.comm
    dx = ufl.dx(domain=sub, metadata=RT["qmeta"])
    d = (T_h - T_INF) - theta
    chi_in = m7._indicator(sub, lambda m: np.sum(m**2, axis=1) < 0.81 * m7.B**2)

    def s(form):
        return float(np.real(m7._allsum(comm, form)))

    whole = np.sqrt(max(s(ufl.inner(d, d) * dx), 0.0) / s(ufl.inner(theta, theta) * dx))
    inner_ = np.sqrt(max(s(chi_in * ufl.inner(d, d) * dx), 0.0)
                     / s(chi_in * ufl.inner(theta, theta) * dx))
    return float(whole), float(inner_)


def _comparand(sub, fn):
    th = m7._p2(sub, fn, p=RT["comparand"])
    nan = sub.comm.allreduce(int(np.count_nonzero(~np.isfinite(th.x.array))), op=MPI.SUM)
    assert nan == 0, f"comparand holds {nan} NaN dofs"
    return th


def main() -> None:
    comm = MPI.COMM_WORLD
    started = time.perf_counter()

    if not np.issubdtype(np.dtype(default_scalar_type), np.complexfloating):
        raise RuntimeError(
            "this example needs the complex DolfinX build "
            "(source /usr/local/bin/dolfinx-complex-mode); the runner's `imp:` "
            "group sources it automatically"
        )

    _say(comm, "=" * 78)
    _say(comm, "EX-168 -- the dose around the implant hole against the duty cycle "
               "(IMP-7 steps 10-11, I49 / I50)")
    _say(comm, "=" * 78)
    _say(comm, f"[gate] bands imported: INTERIOR {INTERIOR_BAND:.0e}, E {E_BAND:.0e}, "
               f"S0 {S0_BAND:.0e}, DOSE {DOSE_BAR:g}, AVG {AVG_BAR:g}, T_REF {T_REF} C; "
               f"rung {RES} m; dt {DT:.4f} s, {N_RUN} steps to {T_END:.4f} s (= 4 t*), "
               f"period {PERIOD:.4f} s, duties {DUTIES}, window {float(W_B2):g} s")

    # Setup figure from the affine twin of the same rung (the curved mesh only adds
    # mid-edge geometry nodes; the solve meshes the shell alone, hole at r = a).
    ma, cta, _fta, _tb = m7._ball(RES)
    write_setup_figure(
        ma,
        cta,
        FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names={
            1: f"adiabatic hole a = 5 mm (cut out of the solve; k {m7.K_CORE} W/m/K in the steady twin)",
            2: f"gel shell to b = 20 mm (k {m7.K_GEL} W/m/K, dipole SAR, Robin-cooled)",
        },
        translucent_tags=(2,),
        slice_normal=(0.0, 1.0, 0.0),
        slice_origin=(0.0, 0.0, 0.0),
        title="imp:19 -- gel shell around an adiabatic hole, the dose against the duty",
        window_size=(1400, 630),
        comm=comm,
    )

    # ---- the fixture's body (`fem_dose`, test_imp7_duty25_dose_fem.py), written out ----------
    sub, sft, q_f, t_build, parent = _shell_fixture(RES)
    ncell = sub.topology.index_map(3).size_global
    _say(comm, f"[mesh] curved shell sub-mesh {ncell} cells, build {t_build:.1f} s; "
               f"ranks {comm.size}")
    cn, t_cn = _transient(sub, sft, q_f, dt=DT, n_steps=N_RUN, theta=0.5,
                          startup_half_steps=STARTUP, source_schedule=_chi(PERIOD, D25),
                          snapshot_steps=tuple(range(N_RUN + 1)))
    e_bal = float(np.real(cn.max_balance_residual))
    s0, tot = _census_25(cn.source_multipliers)
    _say(comm, f"[solve] CN {t_cn:.1f} s; (E) {e_bal:.3e}; (S0) sum chi*dt {tot:.9f} s vs "
               f"{ON_END_25:.9f} s rel {s0:.3e}; {len(cn.source_multipliers)} (half-)steps")

    tq = time.perf_counter()
    th_h = np.empty((N_RUN + 1, M_PTS.shape[0]))
    ok_all = np.ones(M_PTS.shape[0], dtype=bool)
    for k in range(N_RUN + 1):
        th, ok = _probe(cn.snapshots[k], M_PTS, parent)
        th_h[k] = np.real(th)
        ok_all &= ok
    _say(comm, f"[reads] {N_RUN + 1} batches of {M_PTS.shape[0]} in "
               f"{time.perf_counter() - tq:.1f} s; valid {ok_all.tolist()}")
    cn.snapshots.clear()

    # ---- (T) at 4 t*, (c-cont), (c-50) ------------------------------------------------------
    th_duty = _comparand(sub, _duty_fn_at(T_END, D25))
    e_all, e_in = _err_c(sub, cn.temperature, th_duty)
    th_cont = _comparand(sub, _series_fn(T_END, "neumann"))
    _, c_cont = _err_c(sub, cn.temperature, th_cont)
    th_50 = _comparand(sub, _duty_fn_at(T_END, 0.5))
    _, c_50 = _err_c(sub, cn.temperature, th_50)
    _say(comm, f"[(T)] interior rel L2 at 4 t* {e_in:.4e} (whole {e_all:.4e}, printed); "
               f"ASSERTED <= {INTERIOR_BAND:.0e}")
    _say(comm, f"[(c-cont)] continuous series at T_END: interior miss {c_cont:.4e}; "
               f"ASSERTED > {INTERIOR_BAND:.0e} (gate record 0.779)")
    _say(comm, f"[(c-50)] D = 0.5 series at T_END: interior miss {c_50:.4e}; "
               f"ASSERTED > {INTERIOR_BAND:.0e} (gate record 0.520)")

    # ---- series histories: the matrix's closed-form side -----------------------------------
    tq = time.perf_counter()
    m_full = _modes()
    hists = {d: _duty_hist_at(m_full, M_PTS, d) for d in DUTIES}
    th_s = hists[D25]
    _say(comm, f"[series] matrix histories {th_s.shape} x {len(DUTIES)} duties in "
               f"{time.perf_counter() - tq:.1f} s")

    # ---- collect every reading first, then assert ------------------------------------------
    rows = []  # FEM at D = 0.25
    for scale, lab in _scales(D25):
        mx_T = T_INF + np.max(th_h, axis=0) * scale
        _say(comm, f"[max T] FEM D=0.25 {lab}: {np.array2string(mx_T, precision=6)} C "
                   f"(< {T_REF} ASSERTED)")
        for j, name in enumerate(M_LABELS):
            th = th_h[:, j] * scale
            dose, ref, rel, swap, miss = _q_check(T_GRID, T_INF + th)
            mx, t_mx, lo, hi, wref, wrel = _win_check(T_GRID, th)
            rows.append((lab, name, dose, ref, rel, swap, miss, mx, t_mx, lo, hi, wref, wrel,
                         float(mx_T[j])))
            _say(comm, f"[(Q-hd)] {lab} {name}: FEM CEM43 {dose:.6e} s vs quad {ref:.6e} s rel "
                       f"{rel:.3e} (bar {DOSE_BAR:g}, ASSERTED); (c-R) swap miss {miss:.4e} "
                       f"(asserted > {DOSE_BAR:g}), factor {swap / ref:.3e}")
            _say(comm, f"[(Win-hd)] {lab} {name}: max six-minute mean {mx:.6e} K at {t_mx:.2f} s; "
                       f"[{lo:.6e}, {hi:.6e}]; quad {wref:.6e} rel {wrel:.3e} "
                       f"(bar {AVG_BAR:g}, ASSERTED)")

    mrows = []  # the closed-form matrix, 12 cells
    for d in DUTIES:
        for scale, lab in _scales(d):
            for j, name in enumerate(M_LABELS):
                th = hists[d][:, j] * scale
                T = T_INF + th
                dose, ref, rel, swap, miss = _q_check(T_GRID, T)
                mx, t_mx, lo, hi, wref, wrel = _win_check(T_GRID, th)
                mrows.append((d, lab, name, dose, rel, miss, float(np.max(T)), mx, t_mx, lo, hi,
                              wrel))
                _say(comm, f"[(Q-m)] D={d} {lab} {name}: CEM43 {dose:.6e} s rel quad {rel:.3e} "
                           f"(bar {DOSE_BAR:g}, ASSERTED); max T {np.max(T):.6f} C; (c-R) miss "
                           f"{miss:.4e}")
                _say(comm, f"[(Win-m)] D={d} {lab} {name}: max six-minute mean {mx:.6e} K at "
                           f"{t_mx:.2f} s; [{lo:.6e}, {hi:.6e}]; rel quad {wrel:.3e} "
                           f"(bar {AVG_BAR:g}, ASSERTED)")

    # ---- printed only: the matrix and the FEM-vs-series gap ---------------------------------
    W = float(W_B2)
    _say(comm, "[MATRIX] D | point | D(peak) s | D(equal mean) s | six-minute peak K (peak scale) | "
               "six-minute peak K (equal mean) | cap INTERIOR_BAND*max theta_s K")
    for d in DUTIES:
        for j, name in enumerate(M_LABELS):
            th = hists[d][:, j]
            d_pk = cem43(T_GRID, T_INF + th)
            d_eq = cem43(T_GRID, T_INF + th / d)
            mx_pk = _win_check(T_GRID, th)[0]
            mx_eq = _win_check(T_GRID, th / d)[0]
            cap = INTERIOR_BAND * float(np.max(th))
            _say(comm, f"[MATRIX] {d} | {name} | {d_pk:.6e} | {d_eq:.6e} | {mx_pk:.6f} | "
                       f"{mx_eq:.6f} | {cap:.4e}")
            assert np.isfinite(d_eq / d_pk)

    p = ITEM_700_POLE
    _say(comm, f"[GAP] item 700 D=0.5 pole (record, not re-run): peak-scale dose rel {p[0]:.3e}; "
               f"equal-mean rel {p[1]:.3e}; six-minute peak {p[2]:.6f} K vs {p[3]:.6f} K; "
               f"max_n|dtheta| {p[4]:.3e} K; cap {CAP_50_POLE:.4e} K")
    _say(comm, "[GAP] FEM against series at D = 0.25, printed, no band")
    for j, name in enumerate(M_LABELS):
        d_h = cem43(T_GRID, T_INF + th_h[:, j])
        d_s = cem43(T_GRID, T_INF + th_s[:, j])
        de_h = cem43(T_GRID, T_INF + th_h[:, j] / D25)
        de_s = cem43(T_GRID, T_INF + th_s[:, j] / D25)
        _a, mx_h, t_h = running_average(T_GRID, th_h[:, j], window=W)
        _a, mx_s, t_s = running_average(T_GRID, th_s[:, j], window=W)
        M = float(np.max(np.abs(th_h[:, j] - th_s[:, j])))
        cap = INTERIOR_BAND * float(np.max(th_s[:, j]))
        _say(comm, f"  {name}: peak-scale dose {d_h:.6e} s vs {d_s:.6e} s, rel "
                   f"{abs(d_h / d_s - 1.0):.3e}; equal-mean {de_h:.6e} s vs {de_s:.6e} s, rel "
                   f"{abs(de_h / de_s - 1.0):.3e}; six-minute peak {mx_h:.6f} K at {t_h:.2f} s vs "
                   f"{mx_s:.6f} K at {t_s:.2f} s; max_n|dtheta| {M:.3e} K")
        _say(comm, f"  {name}: M_p {M:.4e} K vs cap {cap:.4e} K, cap/M_p {cap / M:.3f} "
                   f"(a reading, no gate)")
        assert np.isfinite(d_h / d_s)

    # ---- assertions ----------------------------------------------------------------------
    assert e_bal <= E_BAND, f"(E) {e_bal:.3e}"
    assert s0 <= S0_BAND, f"(S0) {s0:.3e}"
    assert bool(np.all(ok_all)), f"invalid points {ok_all}"
    assert np.all(np.isfinite(th_h))
    assert e_in <= INTERIOR_BAND, f"(T) {e_in:.4e}"
    assert c_cont > INTERIOR_BAND, f"(c-cont) did not fire: {c_cont:.4e}"
    assert c_50 > INTERIOR_BAND, f"(c-50) did not fire: {c_50:.4e}"
    for (lab, name, dose, ref, rel, swap, miss, mx, t_mx, lo, hi, wref, wrel, tmax) in rows:
        assert tmax < T_REF, (lab, name, tmax)
        assert rel <= DOSE_BAR, ("(Q-hd)", lab, name, rel)
        assert miss > DOSE_BAR, ("(c-R)", lab, name, miss)
        assert T_GRID[0] + W_B2 <= t_mx <= T_GRID[-1], ("(Win-hd) t_max", lab, name, t_mx)
        assert lo <= mx <= hi, ("(Win-hd) bracket", lab, name, lo, mx, hi)
        assert wrel <= AVG_BAR, ("(Win-hd)", lab, name, wrel)
    assert len(mrows) == 12, len(mrows)
    for (d, lab, name, dose, rel, miss, tmax, mx, t_mx, lo, hi, wrel) in mrows:
        cell = (d, lab, name)
        assert tmax < T_REF, (cell, tmax)
        assert rel <= DOSE_BAR, ("(Q-m)", cell, rel)
        assert miss > DOSE_BAR, ("(c-R) matrix", cell, miss)
        assert T_GRID[0] + W_B2 <= t_mx <= T_GRID[-1], ("(Win-m) t_max", cell, t_mx)
        assert lo <= mx <= hi, ("(Win-m) bracket", cell, lo, mx, hi)
        assert wrel <= AVG_BAR, ("(Win-m)", cell, wrel)

    # ---- outputs: P2 Lagrange on the P2-geometry mesh ------------------------------------
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    V2 = fem.functionspace(sub, ("Lagrange", 2))
    ft = fem.Function(V2, name="theta_fem_d025_4tstar")
    ft.interpolate(cn.temperature)
    ft.x.array[:] = np.real(ft.x.array) - T_INF
    ft.x.scatter_forward()
    fs = fem.Function(V2, name="theta_series_d025_4tstar")
    fs.interpolate(th_duty)
    fs.x.array[:] = np.real(fs.x.array)
    fs.x.scatter_forward()
    fe = fem.Function(V2, name="error_d025_4tstar")
    fe.x.array[:] = np.abs(ft.x.array - fs.x.array)
    fe.x.scatter_forward()
    xdmf = OUTPUT_DIR / f"{BASENAME}_combined"
    write_xdmf_with_tags(xdmf, sub, None,
                         {"theta_fem_d025_4tstar": ft, "theta_series_d025_4tstar": fs,
                          "error_d025_4tstar": fe}, comm=comm)
    adopt_host_ownership(OUTPUT_DIR, comm=comm)
    _say(comm, f"\n[output] {xdmf}.xdmf -- P2 fields on the {ncell}-cell curved shell mesh")
    _say(comm, f"\n[timing] total {time.perf_counter() - started:.1f} s at -n {comm.size}")
    _say(comm, "\nAll imported gates hold.")


if __name__ == "__main__":
    main()
