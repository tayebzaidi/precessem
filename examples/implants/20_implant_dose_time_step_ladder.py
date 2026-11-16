"""Example (`EX-174`): the dose gap around the implant hole on a time-step ladder.

`imp:19` reads the CEM43 dose at one time step.  This example refines the step on one
mesh until a reading stops moving: the curved `LADDER[0]` (P2) gel shell around the
adiabatic hole, Crank-Nicolson at dt = DT, DT / 2 and DT / 4 under the D = 0.5 duty
source, the same 161 rows of `T_GRID` read at pole, equator and outer pole, against the
hollow-sphere duty series at `N_DEFAULT` modes.  It narrates `IMP-7` step 9 (item 710,
`tests/validation/test_imp7_duty_dose_ladder.py`) and step 15 (item 745,
`tests/validation/test_imp7_duty_dose_quarter_dt.py`).

**It asserts, and it does not re-implement** (the `ANS-1` rule).  Every band, grid and
helper is imported by name (`_read_times`, `_census_dt`, `_gaps`, `GRID_BAR`, the item 700
and 710 records); the `run` fixtures and every `test_*` name are never imported, so their
body is written out below (disclosed in the guide).  Disclosed local pieces: `_err_c`
(``ufl.inner(d, d)`` with the real part taken, as `imp:17` - `imp:19` do); every integral
and norm is passed through ``np.real`` in this complex build.

Asserted, on each of the three runs: (T) interior rel L2 at 4 t* <= INTERIOR_BAND; (E) <=
E_BAND; (S0) <= S0_BAND; (Grid) <= GRID_BAR; every point valid; max T < T_REF.  Asserted
control (c-cont): the continuous series at T_END misses by more than INTERIOR_BAND.

Printed only: the three-row table at pole, equator and outer pole (peak-scale dose rel,
equal-mean dose rel, six-minute peak delta, max_n |dtheta|), beside the gate modules'
records.  No band and no rate.  The doses are about 1e-11 s: an identity on a history, not
a thermal-damage statement.  The fixture has no metal core, no coil, no tissue.

Run it::

    ./scripts/run_examples.sh -e imp:20 -n 4 -t 560

Output: ``examples/implants/paraview_output/implants_20_implant_dose_time_step_ladder_combined.xdmf``.
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

from tests.validation import test_imp7_implant_heating_fem as m7  # noqa: E402
from tests.validation.test_imp7_duty_cycle_fem import (  # noqa: E402
    DUTY,
    PERIOD,
    S0_BAND,
    _chi,
    _duty_fn,
)
from tests.validation.test_imp7_duty_dose_both import ITEM_710_DT  # noqa: E402
from tests.validation.test_imp7_duty_dose_closed_form import (  # noqa: E402
    DT,
    LABELS,
    PTS,
    SCALES,
    T_END,
    T_GRID,
    _duty_hist,
)
from tests.validation.test_imp7_duty_dose_fem import (  # noqa: E402
    N_RUN,
    STARTUP,
    _probe,
    _shell_fixture,
)
from tests.validation.test_imp7_duty_dose_ladder import (  # noqa: E402
    GRID_BAR,
    ITEM_700,
    _census_dt,
    _gaps,
    _read_times,
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
from tests.validation.test_mat20_tuned_heating_dose import T_REF  # noqa: E402

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"  # committed, EX-57
BASENAME = "implants_20_implant_dose_time_step_ladder"
RES = LADDER[0]
#: (tag, dt as a multiple of DT, step stride that lands on T_GRID)
RUNS = (("DT", 1.0, 1), ("DT/2", 0.5, 2), ("DT/4", 0.25, 4))


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


def _one_run(comm, tag, fac, stride, th_s):
    """The gates' `run` body at dt = fac * DT, written out; returns (record, sub, cn)."""
    dt = fac * DT
    n_steps = N_RUN * stride
    t0 = time.perf_counter()
    sub, sft, q_f, t_build, parent = _shell_fixture(RES)
    ncell = sub.topology.index_map(3).size_global
    _say(comm, f"[{tag}] mesh {ncell} cells, build {t_build:.1f} s; ranks {comm.size}; dt "
               f"{dt:.6f} s, {n_steps} steps to {T_END:.4f} s, reads every {stride} step(s)")
    steps = tuple(range(0, n_steps + 1, stride))
    cn, t_cn = _transient(sub, sft, q_f, dt=dt, n_steps=n_steps, theta=0.5,
                          startup_half_steps=STARTUP, source_schedule=_chi(PERIOD, DUTY),
                          snapshot_steps=steps)
    rec = {"tag": tag, "ncell": ncell}
    rec["E"] = float(np.real(cn.max_balance_residual))
    rec["S0"], tot = _census_dt(cn.source_multipliers, dt)
    _say(comm, f"[{tag}] CN {t_cn:.1f} s; (E) {rec['E']:.3e} (E_BAND {E_BAND:.0e}); (S0) sum "
               f"chi*dt {tot:.9f} s rel {rec['S0']:.3e} (bar {S0_BAND:.0e}); "
               f"{len(cn.source_multipliers)} (half-)steps")

    tq = time.perf_counter()
    th_h = np.empty((len(steps), PTS.shape[0]))
    ok_all = np.ones(PTS.shape[0], dtype=bool)
    for row, k in enumerate(steps):
        th, ok = _probe(cn.snapshots[k], PTS, parent)
        th_h[row] = np.real(th)
        ok_all &= ok
    rec["th"], rec["ok"] = th_h, ok_all
    rec["grid"] = float(np.max(np.abs(_read_times(dt, n_steps, stride)[1:] / T_GRID[1:] - 1.0)))
    _say(comm, f"[{tag}] reads: {len(steps)} batches of {PTS.shape[0]} in "
               f"{time.perf_counter() - tq:.1f} s; valid {ok_all.tolist()}; (Grid) "
               f"{rec['grid']:.3e} (bar {GRID_BAR:.0e})")

    tq = time.perf_counter()
    th_duty = _comparand(sub, _duty_fn(T_END))
    e_all, rec["T"] = _err_c(sub, cn.snapshots[n_steps], th_duty)
    th_cont = _comparand(sub, _series_fn(T_END, "neumann"))
    _, rec["cC"] = _err_c(sub, cn.snapshots[n_steps], th_cont)
    _say(comm, f"[{tag}] (T) interior rel L2 at 4 t* {rec['T']:.4e} (whole {e_all:.4e}, printed); "
               f"ASSERTED <= {INTERIOR_BAND:.0e}; (c-cont) continuous series interior miss "
               f"{rec['cC']:.4e}, ASSERTED > {INTERIOR_BAND:.0e}; "
               f"{time.perf_counter() - tq:.1f} s")
    cn.snapshots.clear()
    rec["gaps"] = _gaps(th_h, th_s)
    rec["elapsed"] = time.perf_counter() - t0
    return rec, sub, cn


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
    _say(comm, "EX-174 -- the dose gap around the implant hole on a time-step ladder "
               "(IMP-7 steps 9 and 15)")
    _say(comm, "=" * 78)
    _say(comm, f"[gate] bands imported: INTERIOR {INTERIOR_BAND:.0e}, E {E_BAND:.0e}, "
               f"S0 {S0_BAND:.0e}, Grid {GRID_BAR:.0e}, T_REF {T_REF} C; rung {RES} m; "
               f"DT {DT:.4f} s, T_END {T_END:.4f} s, duty {DUTY}, period {PERIOD:.4f} s")

    # Setup figure from the affine twin of the same rung (the solve meshes the shell alone).
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
        title="imp:20 -- gel shell around an adiabatic hole, the dose on a time-step ladder",
        window_size=(1400, 630),
        comm=comm,
    )

    tq = time.perf_counter()
    th_s = _duty_hist(_modes())
    _say(comm, f"[series] duty history {th_s.shape} in {time.perf_counter() - tq:.1f} s")

    recs = []
    sub = cn = None
    for tag, fac, stride in RUNS:
        del cn
        rec, sub, cn = _one_run(comm, tag, fac, stride, th_s)
        recs.append(rec)
        _say(comm, f"[{tag}] run elapsed {rec['elapsed']:.1f} s")

    # ---- collect every reading first, then assert ------------------------------------------
    for rec in recs:
        for scale, lab in SCALES:
            mx = T_INF + np.max(rec["th"], axis=0) * scale
            rec[f"maxT_{lab}"] = mx
            _say(comm, f"[{rec['tag']}] max T {lab}: {np.array2string(mx, precision=6)} C "
                       f"(< {T_REF} ASSERTED)")
        _say(comm, f"[gates] {rec['tag']}: (T) {rec['T']:.4e} <= {INTERIOR_BAND:.0e}; (E) "
                   f"{rec['E']:.3e} <= {E_BAND:.0e}; (S0) {rec['S0']:.3e} <= {S0_BAND:.0e}; "
                   f"(Grid) {rec['grid']:.3e} <= {GRID_BAR:.0e}; (c-cont) {rec['cC']:.4e} > "
                   f"{INTERIOR_BAND:.0e}")

    # ---- printed only: the three-row table, beside the gate modules' records ----------------
    gate_rec = {"DT": ITEM_700, "DT/2": ITEM_710_DT}
    _say(comm, "[TABLE] point | run | peak-scale dose rel | equal-mean dose rel | six-minute "
               "peak delta (K) | max_n|dtheta| (K)")
    for name in LABELS:
        for rec in recs:
            r = rec["gaps"][name]
            line = (f"[TABLE] {name} | {rec['tag']} | {r[0]:.3e} | {r[1]:.3e} | {r[2]:+.3e} | "
                    f"{r[3]:.3e}")
            if rec["tag"] in gate_rec:
                g = gate_rec[rec["tag"]][name]
                same = all(f"{a:.3e}" == f"{b:.3e}" for a, b in zip(r, g))
                line += f"   (gate record {g[0]:.3e} | {g[1]:.3e} | {g[2]:+.3e} | {g[3]:.3e}; same at printed digits: {same})"
            _say(comm, line)
            assert np.all(np.isfinite(r))

    # ---- assertions ----------------------------------------------------------------------
    for rec in recs:
        tag = rec["tag"]
        assert rec["E"] <= E_BAND, f"(E) {tag} {rec['E']:.3e}"
        assert rec["S0"] <= S0_BAND, f"(S0) {tag} {rec['S0']:.3e}"
        assert rec["grid"] <= GRID_BAR, f"(Grid) {tag} {rec['grid']:.3e}"
        assert bool(np.all(rec["ok"])), f"invalid points {tag} {rec['ok']}"
        assert np.all(np.isfinite(rec["th"])), tag
        assert rec["th"].shape == T_GRID.shape + (PTS.shape[0],), (tag, rec["th"].shape)
        for scale, lab in SCALES:
            assert np.all(rec[f"maxT_{lab}"] < T_REF), (tag, lab, rec[f"maxT_{lab}"])
        assert rec["T"] <= INTERIOR_BAND, f"(T) {tag} {rec['T']:.4e}"
        assert rec["cC"] > INTERIOR_BAND, f"(c-cont) {tag} did not fire: {rec['cC']:.4e}"

    # ---- outputs: theta at T_END on the finest dt, P2 Lagrange on the curved shell ----------
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    V2 = fem.functionspace(sub, ("Lagrange", 2))
    ft = fem.Function(V2, name="theta_fem_dt4_4tstar")
    ft.interpolate(cn.temperature)
    ft.x.array[:] = np.real(ft.x.array) - T_INF
    ft.x.scatter_forward()
    xdmf = OUTPUT_DIR / f"{BASENAME}_combined"
    write_xdmf_with_tags(xdmf, sub, None, {"theta_fem_dt4_4tstar": ft}, comm=comm)
    adopt_host_ownership(OUTPUT_DIR, comm=comm)
    _say(comm, f"\n[output] {xdmf}.xdmf -- theta at 4 t* on the finest dt, P2, "
               f"{recs[-1]['ncell']}-cell curved shell mesh")
    _say(comm, f"\n[timing] total {time.perf_counter() - started:.1f} s at -n {comm.size}")
    _say(comm, "\nAll imported gates hold.")


if __name__ == "__main__":
    main()
