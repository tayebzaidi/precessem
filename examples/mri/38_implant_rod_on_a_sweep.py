"""Example (`EX-170`): the PEC implant rod ``Study`` on a two-frequency sweep, against its one-frequency solves.

`mri:37` runs the rod ``Study`` at one frequency (64 MHz).  `WF-15` step 11 (item 728) put the same rod
``Study`` through ``Study.solve_sweep``.  No example sweeps a rod.  This one runs

* ``Study.solve_sweep((64 MHz, 128 MHz))`` on ``FourPortBirdcageGeometry(implant_rod=ROD)`` (degree 1,
  four-port ``LumpedPortDrive``, exact-rectangle ports, kappa off): one mesh, one assembly route;
* the two one-frequency ``Study.solve`` calls on the same geometry, one per frequency (each builds its own mesh);
* the sweep keeps the 64 MHz drive fields (``solve_sweep(..., fields_at=64 MHz)``, the public route) for the
  ParaView file only; there is no separate view solve.

Gates, collected first and then asserted (every band imported by name):

* (Cells) 168 507 on every build (sweep and point, both frequencies), equal to ``ROD_EXACT_CELL_RECORD``;
* (S) sweep S equals the point-solve S within ``IDENTITY_BAND`` at each frequency (matrix norm and worst entry);
* (i) reciprocity within ``RECIPROCITY_BAND``; (ii) sigma_max <= 1 + ``PASSIVITY_SIGMA_TOLERANCE`` and
  column power <= 1, at each frequency;
* (c-f), the negative control: the sweep's 64 MHz matrix misses the 128 MHz point solve by more than
  ``IDENTITY_BAND`` (the gate's record is 5.52e-01, so the bar is reachable).

Printed only: sigma_max and the worst column power per frequency; the sweep time against the two point solves
(a record, not a promise); item 709's 64 MHz sigma_max beside this run's.

Run it through the example runner::

    ./scripts/run_examples.sh -e mri:38 -n 2 -t 560

Not claimed: a tuned coil, a resonance, a tip Q ordering across frequency, the two-port drive, a tissue model.
"""

from __future__ import annotations

import dataclasses
import gc
import sys
import time
from pathlib import Path

import numpy as np
import ufl
from mpi4py import MPI

from dolfinx import fem

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from fem_em_solver import FourPortBirdcageGeometry, LumpedPortDrive, Study  # noqa: E402
from fem_em_solver.io.paraview_utils import (  # noqa: E402
    adopt_host_ownership,
    write_xdmf_with_tags,
)
from fem_em_solver.post.setup_figure import write_setup_figure  # noqa: E402

from tests.mesh.test_birdcage_implant_rod import ROD  # noqa: E402
from tests.mesh.test_birdcage_leg_offset import IDENTITY_BAND  # noqa: E402
from tests.package.test_wf15_study_implant_rod import (  # noqa: E402
    RECORD_WIDTH,
    ROD_EXACT_CELL_RECORD,
    _phantom,
    _study,
)
from tests.package.test_wf15_study_rod_sweep import (  # noqa: E402
    FREQS,
    ITEM_709_SIGMA,
    _rel,
)
from tests.validation.test_imp3_rod_birdcage import _matrix_quantities  # noqa: E402
from tests.validation.test_port_birdcage_four_port import (  # noqa: E402
    PASSIVITY_SIGMA_TOLERANCE,
    POOLED_SEPARATION_FLOOR,
)
from tests.validation.test_port_lumped_sheet_sweep import RECIPROCITY_BAND  # noqa: E402

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"
BASENAME = "mri_38_implant_rod_on_a_sweep"

# The gate's records, printed beside the readings (never asserted by equality).
REC_GATE_LOG = "20261101T195622Z_WF-15-s11-W1.log"
REC_CONTROL_MISS = 5.5237e-01


def _p0(comm, text):
    if comm.rank == 0:
        print(text, flush=True)


def main() -> None:
    comm = MPI.COMM_WORLD
    started = time.perf_counter()
    assert comm.size == RECORD_WIDTH, f"the record holds at -n {RECORD_WIDTH}, got -n {comm.size}"
    OUTPUT_DIR.mkdir(exist_ok=True)
    FIGURE_DIR.mkdir(exist_ok=True)
    comm.Barrier()

    # ---- the sweep, then the two one-frequency solves (the gate's route, written out) -------------
    geom = FourPortBirdcageGeometry(implant_rod=ROD)
    t0 = time.perf_counter()
    sweep = _study(geom).solve_sweep(FREQS, comm, fields_at=FREQS[0])
    t_sweep = time.perf_counter() - t0

    t0 = time.perf_counter()
    points = {}
    for f in FREQS:
        material = dataclasses.replace(_phantom(), frequency_hz=f)
        points[f] = Study(geom, material, LumpedPortDrive(), frequency_hz=f, degree=1).solve(comm)
    t_points = time.perf_counter() - t0
    _p0(comm, f"[time] -n {comm.size}: sweep {t_sweep:.1f} s (mesh {sweep.mesh_time:.1f}, solves "
              f"{[round(t, 1) for t in sweep.sweep_times]}); two point Studies {t_points:.1f} s "
              f"(a record, not a promise)")

    # ---- setup figure (no-op unless FEM_EM_SETUP_FIGURES=1) ---------------------------------------
    first = sweep.results[0]
    msh, cell_tags = first.mesh, first.cell_tags
    local_tags = np.unique(cell_tags.values)
    all_tags = sorted({int(x) for part in comm.allgather(local_tags.tolist()) for x in part})
    phantom_tag = geom.phantom_tag
    other = tuple(x for x in all_tags if x != phantom_tag)
    names = {phantom_tag: "gel phantom (rod hole)"}
    for x in other:
        names[x] = "air / copper regions"
    write_setup_figure(
        msh, cell_tags, FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names=names, hide_tags=other, translucent_tags=(phantom_tag,),
        slice_normal=(0.0, 1.0, 0.0), slice_origin=(0.0, 0.0, 0.0),
        title="mri:38 -- F-small phantom with the PEC rod hole, swept (degree 1)",
        comm=comm)
    gc.collect()
    comm.Barrier()

    # ---- the gates: collect everything, then assert ----------------------------------------------
    gates = {}
    for f, res in zip(FREQS, sweep.results):
        mhz = f / 1e6
        point = points[f]
        mq = _matrix_quantities(res)
        sigma = float(mq["sigma"][0])
        colpow_max = float(np.max(mq["colpow"]))
        rel = _rel(res.s_matrix, point.s_matrix)
        entry = float(np.max(np.abs(np.asarray(res.s_matrix) - np.asarray(point.s_matrix))
                             / np.abs(np.asarray(point.s_matrix))))
        gates[f"{mhz:g} S"] = rel <= IDENTITY_BAND and entry <= IDENTITY_BAND
        gates[f"{mhz:g} i"] = mq["recip"] <= RECIPROCITY_BAND
        gates[f"{mhz:g} ii_sigma"] = sigma <= 1.0 + PASSIVITY_SIGMA_TOLERANCE
        gates[f"{mhz:g} ii_colpow"] = bool(np.all(mq["colpow"] <= 1.0))
        gates[f"{mhz:g} cells sweep"] = res.cells == ROD_EXACT_CELL_RECORD
        gates[f"{mhz:g} cells point"] = point.cells == ROD_EXACT_CELL_RECORD
        _p0(comm, f"[(S)] {mhz:g} MHz: sweep vs point rel {rel:.3e}, worst entry {entry:.3e} "
                  f"(ASSERTED <= IDENTITY_BAND {IDENTITY_BAND:.0e}); record {REC_GATE_LOG}")
        _p0(comm, f"[(i)] {mhz:g} MHz: {mq['recip']:.3e} (ASSERTED <= RECIPROCITY_BAND "
                  f"{RECIPROCITY_BAND:.0e}); [(ii)] sigma_max {sigma:.12f} (ASSERTED <= 1 + "
                  f"{PASSIVITY_SIGMA_TOLERANCE:g}), column power max {colpow_max:.7f} (ASSERTED <= 1)")
        _p0(comm, f"[(Cells)] {mhz:g} MHz: sweep {res.cells}, point {point.cells}, record "
                  f"{ROD_EXACT_CELL_RECORD} (ASSERTED equal)")
    sig64 = float(_matrix_quantities(sweep.results[0])["sigma"][0])
    _p0(comm, f"[709] 64 MHz sigma_max {sig64:.12f} vs item 709's {ITEM_709_SIGMA:.12f}: |diff| "
              f"{abs(sig64 - ITEM_709_SIGMA):.3e} (printed, not asserted here)")

    miss = _rel(sweep.results[0].s_matrix, points[FREQS[1]].s_matrix)
    gates["c_f"] = miss > IDENTITY_BAND
    _p0(comm, f"[(c-f)] the sweep's 64 MHz matrix misses the 128 MHz point solve by {miss:.4e} = "
              f"{miss / IDENTITY_BAND:.3e} x IDENTITY_BAND {IDENTITY_BAND:.0e} (ASSERTED >); gate record "
              f"{REC_CONTROL_MISS:.4e}; POOLED_SEPARATION_FLOOR {POOLED_SEPARATION_FLOOR:g} (printed)")

    red = sorted(comm.allreduce([k for k, v in gates.items() if not v], op=MPI.SUM))
    _p0(comm, f"[gates] red: {red or 'none'} of {len(gates)}")

    # ---- ParaView: the sweep's own 64 MHz result, kept by fields_at (not a gate comparand) ---------
    f_view = FREQS[0]
    view = sweep.results[0]
    assert view.fields is not None and float(sweep.frequencies_hz[0]) == f_view
    _p0(comm, f"[paraview fields] kept by solve_sweep(fields_at={f_view:g}) on the sweep's own solve "
              f"(cells {view.cells}; the file only, no gate reads it)")
    vmsh = view.mesh
    e = view.fields["P1"].e_complex
    p1_space = fem.functionspace(vmsh, ("Lagrange", 1))
    e_abs = fem.Function(p1_space, name="E_abs")
    e_abs.interpolate(
        fem.Expression(ufl.sqrt(ufl.real(ufl.inner(e, e))), p1_space.element.interpolation_points))
    e_abs.x.array[:] = np.real(e_abs.x.array)
    e_abs.x.scatter_forward()
    dg0 = fem.functionspace(vmsh, ("DG", 0))
    sigma_out = fem.Function(dg0, name="sigma")
    sigma_out.interpolate(view.fields["P1"].sigma_field)
    sigma_out.x.array[:] = np.real(sigma_out.x.array)
    sigma_out.x.scatter_forward()
    xdmf_path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_combined", vmsh, view.cell_tags,
        {"E_abs": e_abs, "sigma": sigma_out}, comm=comm)
    adopt_host_ownership(OUTPUT_DIR, comm=comm)
    _p0(comm, f"[paraview] wrote {xdmf_path} (|E| of the P1 drive at 64 MHz with the rod, sigma, cell tags)")

    assert not red, red
    _p0(comm, f"\nAll assertions hold. Total elapsed {time.perf_counter() - started:.1f} s.")


if __name__ == "__main__":
    main()
