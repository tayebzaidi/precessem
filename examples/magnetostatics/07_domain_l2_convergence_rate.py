#!/usr/bin/env python3
"""Sampler-free convergence rate of the straight-wire solve (`EX-148`).

`mag:6` (`EX-9`) measures the h-convergence of the straight wire with a *sampled*
statistic (ten points on a line); that statistic swung 34% with its sample count
and `MAG-19` retired its rate band. The rate duty moved to `MAG-18`, whose
statistic is the annulus-restricted relative domain L2 error of ``|B|``

    E_Omega = || |B_h| - |B_ana| ||_{L2(Omega)} / || |B_ana| ||_{L2(Omega)},
    Omega = {2a <= r <= 0.8 R, |z| <= 0.25 L},

an integral over cells with no sample count in it. This example runs that gate's
own measurement as a narrated example:

* the ladder h = 0.004 / 0.0025 / 0.0018 with the analytic Dirichlet wall;
* the fitted rate, asserted against the imported one-sided floor;
* the natural (n x H = 0) wall at h = 0.0025, asserted strictly worse.

Everything -- ladder, floor, record, helpers, the rate fit -- is *imported by
name* from the gate modules (`ANS-1`), never restated. No ``test_*`` name is
imported. The h = 0.0025 record is **printed beside its reading and not
asserted** (it is version-tagged; `test_domain_l2_record` owns it).

Real build only (magnetostatics); do not source the complex mode. Heavy tier,
``-n 2``. Degree 1 only: see the `MAG-10` note in ``_solve_straight_wire``.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
from dolfinx import fem
from dolfinx import mesh as dmesh
from mpi4py import MPI

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from fem_em_solver.io.paraview_utils import (  # noqa: E402
    write_combined_paraview_output,
)
from fem_em_solver.post.setup_figure import write_setup_figure  # noqa: E402
from fem_em_solver.utils.analytical import AnalyticalSolutions  # noqa: E402
from tests.validation.test_convergence import fit_convergence_rate  # noqa: E402
from tests.validation.test_straight_wire import (  # noqa: E402
    CURRENT,
    DOMAIN_RADIUS,
    E_OMEGA_H0025_RECORD,
    E_OMEGA_LADDER,
    E_OMEGA_RATE_MIN,
    WIRE_LENGTH,
    WIRE_RADIUS,
    _annulus_indicator,
    _domain_l2_error,
    _solve_straight_wire,
)

# Measured when `MAG-18` ran green (20260824T003650Z_MAG-18-regate-run2.log.gz).
# Printed beside this run's readings; never asserted.
R2_ERRORS = [2.5286824423e-01, 1.0617169964e-01, 6.6458072979e-02]
R2_CELLS = [38740, 147235, 383146]
R2_RATE = 1.6854
R2_NATURAL_H0025 = 0.32315493

CONTROL_H = 0.0025

OUTPUT_DIR = "paraview_output"
BASENAME = "magnetostatics_07_domain_l2_convergence_rate"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"


def _wire_air_tags(mesh):
    """Cell tags 1 = wire (r <= a), 2 = air, from owned-cell midpoints.

    `_solve_straight_wire` returns no tags; the picture needs them.
    """
    tdim = mesh.topology.dim
    n_owned = mesh.topology.index_map(tdim).size_local
    owned = np.arange(n_owned, dtype=np.int32)
    mid = dmesh.compute_midpoints(mesh, tdim, owned)
    r = np.hypot(mid[:, 0], mid[:, 1])
    values = np.where(r <= WIRE_RADIUS, 1, 2).astype(np.int32)
    return dmesh.meshtags(mesh, tdim, owned, values)


def main() -> None:
    comm = MPI.COMM_WORLD
    rank0 = comm.rank == 0
    t_start = time.time()

    if rank0:
        print("=" * 74)
        print("Annulus domain-L2 convergence rate of the straight-wire solve")
        print("=" * 74)
        print("  fixture        : tests/validation/test_straight_wire.py (MAG-18)")
        print(f"  current        : {CURRENT:.3f} A, a = {WIRE_RADIUS} m, R = {DOMAIN_RADIUS} m")
        print(f"  wire length    : {WIRE_LENGTH} m")
        print(f"  ladder h       : {E_OMEGA_LADDER} m")
        print(f"  rate floor     : {E_OMEGA_RATE_MIN} (imported, one-sided)")
        print(f"  MPI ranks      : {comm.size}")
        print(flush=True)

    errors = []
    n_cells = []
    times = []
    mesh = None
    b_field = None
    for i, h in enumerate(E_OMEGA_LADDER):
        t0 = time.time()
        mesh, b_field = _solve_straight_wire(h, comm)
        e = _domain_l2_error(mesh, b_field)
        dt = time.time() - t0
        nc = int(mesh.topology.index_map(mesh.topology.dim).size_global)
        errors.append(e)
        n_cells.append(nc)
        times.append(dt)
        if rank0:
            print(f"  h = {h:<7.4f} {nc:>8d} cells   E_Omega = {e:>9.4%}   ({dt:.1f} s)", flush=True)
        if i < len(E_OMEGA_LADDER) - 1:
            # Only the finest rung is kept for export.
            mesh, b_field = None, None

    rate = fit_convergence_rate(np.array(E_OMEGA_LADDER), np.array(errors))

    # ---- control: natural wall at the middle rung ---------------------------
    t0 = time.time()
    mesh_n, b_nat = _solve_straight_wire(CONTROL_H, comm, analytic_bc=False)
    e_natural = _domain_l2_error(mesh_n, b_nat)
    dt_nat = time.time() - t0
    e_analytic = errors[E_OMEGA_LADDER.index(CONTROL_H)]
    del mesh_n, b_nat

    # ---- setup figure + export of the finest rung ---------------------------
    cell_tags = _wire_air_tags(mesh)
    write_setup_figure(
        mesh,
        cell_tags,
        FIGURE_DIR / "magnetostatics_07_domain_l2_convergence_rate_setup.png",
        region_names={1: "wire (conductor)", 2: "air"},
        hide_tags=(2,),
        slice_normal=(0.0, 0.0, 1.0),
        slice_origin=(0.0, 0.0, 0.0),
        title="mag:7 -- Omega-restricted L2 rate, finest rung (h = %.4f m)" % E_OMEGA_LADDER[-1],
        comm=comm,
    )

    v_cg = fem.functionspace(mesh, ("Lagrange", 1, (3,)))
    b_h_cg = fem.Function(v_cg, name="B_numeric")
    b_h_cg.interpolate(b_field)

    # The analytic field interpolated as `_domain_l2_error` does it.
    b_ana_cg = fem.Function(v_cg, name="B_analytic")

    def _b_ana_interp(x):
        points = np.ascontiguousarray(x[:3].T)
        return AnalyticalSolutions.straight_wire_magnetic_field(points, CURRENT).T

    b_ana_cg.interpolate(_b_ana_interp)
    chi = _annulus_indicator(mesh)

    written_files = write_combined_paraview_output(
        OUTPUT_DIR,
        BASENAME,
        mesh,
        cell_tags,
        {
            "B_numeric": (b_h_cg, b_h_cg),
            "B_analytic": (b_ana_cg, b_ana_cg),
            "chi_Omega": (chi, chi),
        },
        comm=comm,
    )

    record = {
        "E_Omega": errors,
        "n_cells": n_cells,
        "rate": rate,
        "E_Omega_natural_h0025": e_natural,
        "natural_over_analytic": e_natural / e_analytic,
        "E_Omega_h0025": e_analytic,
    }

    if rank0:
        print()
        print("-" * 74)
        print("Reproduction records (printed beside the readings; not asserted)")
        print("-" * 74)
        print(f"{'h (m)':>9} {'cells':>9} {'E_Omega':>13} {'MAG-18 cells':>13} {'MAG-18 E_Omega':>16}")
        for i, h in enumerate(E_OMEGA_LADDER):
            print(
                f"{h:>9.4f} {n_cells[i]:>9d} {errors[i]:>13.10e} "
                f"{R2_CELLS[i]:>13d} {R2_ERRORS[i]:>16.10e}"
            )
        print()
        print(f"  fitted rate        : {rate:.4f}   (MAG-18 record {R2_RATE}; floor {E_OMEGA_RATE_MIN})")
        print(f"  rate / floor       : {rate / E_OMEGA_RATE_MIN:.2f}x")
        print(
            f"  h = 0.0025 record  : reading {errors[1]:.6e} beside E_OMEGA_H0025_RECORD "
            f"{E_OMEGA_H0025_RECORD:.6e} (diff {errors[1] - E_OMEGA_H0025_RECORD:+.3e}; "
            "not asserted)"
        )
        print(
            f"  natural wall h={CONTROL_H}: {e_natural:.6%} vs analytic {e_analytic:.6%}, "
            f"ratio analytic/natural {e_analytic / e_natural:.4f} "
            f"(MAG-18 record {R2_NATURAL_H0025:.6f}; solve {dt_nat:.1f} s)"
        )
        print(f"  RECORD             : {record}")
        print(f"  wall time per rung : {[round(t, 1) for t in times]} s")
        print(flush=True)

    # ---- gates ---------------------------------------------------------------
    for i in range(1, len(errors)):
        assert errors[i] < errors[i - 1], (
            f"(M) E_Omega rose from {errors[i - 1]:.4%} at h = {E_OMEGA_LADDER[i - 1]} to "
            f"{errors[i]:.4%} at h = {E_OMEGA_LADDER[i]}"
        )
    assert rate >= E_OMEGA_RATE_MIN, (
        f"(Rate) fitted rate {rate:.4f} below the floor {E_OMEGA_RATE_MIN}"
    )
    assert e_analytic < e_natural, (
        f"(c-BC) analytic-wall E_Omega {e_analytic:.4%} is not below the natural wall's "
        f"{e_natural:.4%} at h = {CONTROL_H}: the analytic wall did not reach the solve"
    )

    if rank0:
        print("=" * 74)
        print("ParaView output (finest mesh, h = %.4f m)" % E_OMEGA_LADDER[-1])
        print("=" * 74)
        for name, path in sorted(written_files.items()):
            print(f"  {name:>14}: {path}")
        print()
        print(
            f"All assertions hold. Total elapsed {time.time() - t_start:.1f} s "
            f"({comm.size} ranks, {n_cells[-1]} cells at the finest h)."
        )


if __name__ == "__main__":
    main()
