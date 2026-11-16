"""Example (`EX-151`): a time-step convergence rate under a switched source --
backward Euler on the duty-cycled gel ball converges at order about 1.00, and
Crank-Nicolson's order 2 falls outside that band.

`mat:8` (`EX-74`) prints `MAT-9`'s BE / CN rates on a *constant* source;
`mat:11` (`EX-85`) shows the duty cycle at one time step. This example runs
`MAT-17` step 2 gate (ii), the Delta-t ladder 300 / 150 / 75 / 37.5 / 18.75 s
on `MAT-9`'s 0.01 m ball under case A's schedule (P = 1200 s, D = 0.5, to
t_end = 14 100 s), and fits the whole-ball self-difference rate.

**It asserts, and it does not re-implement** (the `ANS-1` rule). Bands, ladder,
fixture and helpers are imported by name from
``tests/validation/test_mat17_duty_cycle_fem.py`` (`m17`) and the complex-safe
norms of ``test_mat13_transient_perfusion_sphere.py`` (`m13`). The gate's
``_ladder`` is *not* called: it uses `MAT-9`'s ``_rel_diff``, which takes
``float()`` of a complex value and raises in the complex build. The loop is
written out here with every call the gate makes; ``m13._rel_diff`` is
`MAT-9`'s formula with ``np.real`` before the square root. No gate edit, no
``src/`` edit.

Asserted: (S0) sum chi*dt = on-time on every run (inside ``m17._census``);
(iii) every BE and CN run's balance <= ``BALANCE_BAND``; (ii) BE least-squares
rate within ``BE_RATE_BAND``; (c-CN) the CN + Rannacher rate on the same
ladder lies above ``BE_RATE_BAND[1]``. Printed only: per-run N, dt, wall,
balance; per-pair d and pairwise rates; the BE N = 752 interior error against
the closed form beside ``INTERIOR_BAND``.

Run it::

    ./run_examples.sh -e mat:27 -n 2 -t 420

Output: ``examples/materials/paraview_output/
materials_27_duty_cycle_time_step_rate_finest_combined.xdmf`` (finest BE
theta_h beside the closed form at t_end) and
``materials_27_duty_cycle_time_step_rate_ladder.csv``.
"""

from __future__ import annotations

import csv
import sys
import time
from pathlib import Path

import numpy as np
from mpi4py import MPI

from dolfinx import default_scalar_type, fem

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from fem_em_solver.io.paraview_utils import (  # noqa: E402
    adopt_host_ownership,
    write_xdmf_with_tags,
)
from fem_em_solver.post.setup_figure import write_setup_figure  # noqa: E402

from tests.validation import test_mat13_transient_perfusion_sphere as m13  # noqa: E402
from tests.validation import test_mat17_duty_cycle_fem as m17  # noqa: E402

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"  # committed, EX-57
BASENAME = "materials_27_duty_cycle_time_step_rate"
RES = 0.01

#: Gate records, real build -- `20261003T171113Z_MAT-17-step2-W1b.log:110-144`.
#: Printed beside this run, a difference reported and never asserted.
RECORD = {
    "be_d": (2.6051e-03, 1.3085e-03, 6.5583e-04, 3.2832e-04),
    "be_rate": 0.996,
    "be_pair": (0.993, 0.997, 0.998),
    "cn_rate": 2.104,
    "cn_pair": (2.298, 2.034, 2.005),
    "be_balance": 9.960e-14,
    "cn_balance": 2.358e-13,
    "e752": 1.8819e-03,
}


def _say(comm, msg):
    if comm.rank == 0:
        print(msg, flush=True)


def _ladder(comm, msh, ft, theta, startup, label):
    """The gate's `_ladder` (:295-317), with the complex-safe `m13._rel_diff`."""
    chi = m17._chi(m17.P_A, m17.D_A)
    sols, rows, worst, last = [], [], 0.0, None
    t0 = time.perf_counter()
    for dt in m17.BE_DTS:
        n = m17._n(m17.T_END_A, dt)
        s = m17._run(msh, ft, dt, n, theta, startup, schedule=chi)
        m17._census(s, dt, startup, m17.ON_TIME_A, f"(ii) {label} dt={dt}")
        worst = max(worst, s.max_balance_residual)
        _say(comm, f"[run] {label} dt={dt} N={n} max_balance={s.max_balance_residual:.3e} "
                   f"wall {s.wall:.1f}s")
        rows.append((label, dt, n, s.max_balance_residual, s.wall))
        sols.append(s.temperature)
        last = (s, n)
    d = [m13._rel_diff(msh, sols[i], sols[i + 1]) for i in range(len(sols) - 1)]
    lsq = float(np.polyfit(np.log(np.asarray(m17.BE_DTS[: len(d)])), np.log(d), 1)[0])
    pair = [float(np.log(d[i] / d[i + 1]) / np.log(2.0)) for i in range(len(d) - 1)]
    _say(comm, f"[ladder] {label} d (dt=300,150,75,37.5) = {[f'{x:.4e}' for x in d]} "
               f"LSQ rate={lsq:.3f} pairwise={[f'{p:.3f}' for p in pair]} "
               f"({time.perf_counter() - t0:.1f}s)")
    return lsq, d, pair, worst, last, rows


def main() -> None:
    comm = MPI.COMM_WORLD
    started = time.perf_counter()

    if not np.issubdtype(np.dtype(default_scalar_type), np.complexfloating):
        raise RuntimeError(
            "this example needs the complex DolfinX build "
            "(source /usr/local/bin/dolfinx-complex-mode); the runner's `mat:` "
            "group sources it automatically"
        )

    _say(comm, "=" * 78)
    _say(comm, "EX-151 -- time-step rate under a switched source (MAT-17 step 2 gate (ii))")
    _say(comm, "=" * 78)
    _say(comm, f"\n[gate] imported from MAT-17: BE_DTS {m17.BE_DTS}, BE_RATE_BAND "
               f"{m17.BE_RATE_BAND}, BALANCE_BAND {m17.BALANCE_BAND:.0e}, INTERIOR_BAND "
               f"{m17.INTERIOR_BAND:.0e}; case A P={m17.P_A:g} s D={m17.D_A:g} "
               f"t_end={m17.T_END_A:g} s on-time={m17.ON_TIME_A:g} s")

    msh, cell_tags, ft = m17._ball(RES)
    n_cells = msh.topology.index_map(3).size_global
    _say(comm, f"[mesh] MAT-9's _ball({RES}): {n_cells} cells")
    write_setup_figure(
        msh,
        cell_tags,
        FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names={1: "gel ball (Robin h = 10, Bi = 5/3; duty-cycled uniform source)"},
        translucent_tags=(1,),
        slice_normal=(0.0, 0.0, 1.0),
        slice_origin=(0.0, 0.0, 0.0),
        title="mat:27 -- 10 cm Robin gel ball, 0.01 m mesh, duty-cycled source (case A)",
        comm=comm,
    )

    lsq_be, d_be, pair_be, w_be, (s_fin, n_fin), rows_be = _ladder(
        comm, msh, ft, 1.0, 0, "BE")
    ex = m17._exact_duty(msh, m17.T_END_A, m17.P_A, m17.D_A)
    e752 = m13._interior_error(msh, s_fin.temperature, ex)
    _say(comm, f"[closed form] BE N={n_fin} interior_rel_L2 vs closed form = {e752:.4e} "
               f"(record {RECORD['e752']:.4e}; band {m17.INTERIOR_BAND:.0e}; printed only)")
    lsq_cn, d_cn, pair_cn, w_cn, _last, rows_cn = _ladder(
        comm, msh, ft, 0.5, 2, "CN+Rannacher")

    _say(comm, f"\n[records] BE d {RECORD['be_d']}, rate {RECORD['be_rate']}, pairwise "
               f"{RECORD['be_pair']}; CN rate {RECORD['cn_rate']}, pairwise {RECORD['cn_pair']}")
    _say(comm, f"[rate] (ii) BE LSQ rate = {lsq_be:.3f} (record {RECORD['be_rate']:.3f}, diff "
               f"{lsq_be - RECORD['be_rate']:+.1e}) ASSERTED in {m17.BE_RATE_BAND}")
    _say(comm, f"[rate] (c-CN) CN+Rannacher LSQ rate = {lsq_cn:.3f} (record "
               f"{RECORD['cn_rate']:.3f}, diff {lsq_cn - RECORD['cn_rate']:+.1e}) ASSERTED > "
               f"{m17.BE_RATE_BAND[1]} (margin {lsq_cn - m17.BE_RATE_BAND[1]:+.3f})")
    _say(comm, f"[balance] (iii) worst BE {w_be:.3e} (record {RECORD['be_balance']:.3e}), "
               f"CN {w_cn:.3e} (record {RECORD['cn_balance']:.3e}) ASSERTED <= "
               f"{m17.BALANCE_BAND:.0e}")

    # ---- output ---------------------------------------------------------------------
    V = s_fin.temperature.function_space
    th_h = fem.Function(V, name="theta_h")
    th_h.x.array[:] = s_fin.temperature.x.array - m17.T_INF
    th_h.x.scatter_forward()
    th_cf = fem.Function(V, name="theta_closed_form")
    th_cf.interpolate(ex)
    th_cf.x.scatter_forward()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_finest_combined", msh, None,
        {"theta_h": th_h, "theta_closed_form": th_cf}, comm=comm,
    )
    if comm.rank == 0:
        csv_path = OUTPUT_DIR / f"{BASENAME}_ladder.csv"
        with csv_path.open("w", newline="") as fh:
            wr = csv.writer(fh)
            wr.writerow(["scheme", "dt_s", "N", "max_balance", "wall_s"])
            for r in rows_be + rows_cn:
                wr.writerow([r[0], f"{r[1]:g}", r[2], f"{r[3]:.6e}", f"{r[4]:.3f}"])
            wr.writerow([])
            wr.writerow(["scheme", "pair", "d", "pairwise_rate_to_next"])
            for lab, d, pr in (("BE", d_be, pair_be), ("CN+Rannacher", d_cn, pair_cn)):
                for i, x in enumerate(d):
                    wr.writerow([lab, f"{m17.BE_DTS[i]:g}->{m17.BE_DTS[i + 1]:g}",
                                 f"{x:.6e}", f"{pr[i]:.4f}" if i < len(pr) else ""])
            wr.writerow([])
            wr.writerow(["scheme", "lsq_rate"])
            wr.writerow(["BE", f"{lsq_be:.4f}"])
            wr.writerow(["CN+Rannacher", f"{lsq_cn:.4f}"])
    adopt_host_ownership(OUTPUT_DIR, comm=comm)
    _say(comm, f"\n[output] {path} -- finest BE run: P1 theta_h beside P1 closed form at t_end")
    _say(comm, f"[output] {OUTPUT_DIR / (BASENAME + '_ladder.csv')}")

    # ---- assertions (the gate's) ----------------------------------------------------
    assert w_be <= m17.BALANCE_BAND and w_cn <= m17.BALANCE_BAND, (
        f"(iii) BE {w_be:.3e}, CN {w_cn:.3e}")
    assert m17.BE_RATE_BAND[0] <= lsq_be <= m17.BE_RATE_BAND[1], (
        f"(ii) BE LSQ rate {lsq_be:.3f}, pairwise {pair_be}")
    assert lsq_cn > m17.BE_RATE_BAND[1], f"(c-CN) CN rate {lsq_cn:.3f} inside the BE band"

    _say(comm, f"\n[timing] total {time.perf_counter() - started:.1f} s at -n {comm.size}")
    _say(comm, "\nAll imported gates hold.")


if __name__ == "__main__":
    main()
