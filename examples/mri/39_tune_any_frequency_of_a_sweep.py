"""Example (`EX-172`): sweep once, tune at any of its frequencies, with no second solve -- one
``Study.solve_sweep`` written to one ``.s4p``, then ``Study.load_result`` + ``tune(1, at_hz=f)`` and
``femem run``'s ``tune.at_hz`` on that file, on a clean clone.

`WF-16` step 17 (U108) let a ``Study`` load a multi-frequency four-port Touchstone and select a
column with ``tune(k, at_hz=f)``; `WF-1` step 27 (U111) put ``at_hz`` behind the ``tune:`` block of a
``femem run`` file.  Both gates read a gitignored file left by another gate's window, so nothing on
a clean clone shows the route.  This example makes its own file:

* one F-small four-port sweep (degree 1, PEC proxy, exact-rectangle ports, 10 / 64 / 128 MHz,
  ``-n 2``) through ``Study.solve_sweep`` and ``write_touchstone``;
* a fresh ``Study`` + ``load_result`` + ``tune(1, at_hz=f)`` at each frequency, and
  ``run_study_file`` on a YAML written beside the file, once per ``at_hz``;
* (C=) each ``C_mode`` against ``ITEM722_C_MODE_F`` at ``REDUCTION_BAND``;
* (I=) API = YAML by ``float.hex`` at each frequency;
* (f=) the selected frequency equals the requested one exactly;
* (S0) 0 mesh builds and 0 solves after the sweep (a raising wrapper, restored in ``finally``);
* (c-next) the column one frequency along misses by more than ``REDUCTION_BAND`` (asserted);
* (Ref) an ``at_hz`` outside the span raises ``ValueError`` through the API (asserted).

Every band and record is imported by name from the gate modules (the ANS-1 rule); the control
record printed beside the reading sits in the example-local ``REC_*`` with its log and is never
asserted by equality.  Run it through the example runner::

    ./scripts/run_examples.sh -e mri:39 -n 2 -t 300

Not claimed: interpolation between frequencies, a resonance, a match or Q, degree 2, copper.
"""

from __future__ import annotations

import copy
import gc
import shutil
import sys
import time
from pathlib import Path

import gmsh
import numpy as np
import yaml
from mpi4py import MPI

from dolfinx import fem

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from fem_em_solver import LumpedPortDrive, Study, run_study_file  # noqa: E402
from fem_em_solver.io.paraview_utils import (  # noqa: E402
    adopt_host_ownership,
    write_xdmf_with_tags,
)
from fem_em_solver.post.setup_figure import write_setup_figure  # noqa: E402
from fem_em_solver.study import FourPortBirdcageGeometry  # noqa: E402

from tests.package.test_wf16_load_result_sweep import FREQS  # noqa: E402
from tests.package.test_wf16_study_sweep import _port_study  # noqa: E402
from tests.package.test_wf1_yaml_load_result import _hex, _tune_block_hex_equal  # noqa: E402
from tests.package.test_wf1_yaml_load_result_sweep import ITEM722_C_MODE_F  # noqa: E402
from tests.validation.test_lossy_sphere_fullwave import SALINE_SIGMA  # noqa: E402
from tests.validation.test_port_lumped_rlc_termination import REDUCTION_BAND  # noqa: E402

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"
BASENAME = "mri_39_tune_any_frequency_of_a_sweep"
RUN_DIR = OUTPUT_DIR / "mri_39_run"
SWEEP_YAML = _REPO_ROOT / "tests" / "package" / "wf1" / "f_small_sweep.yaml"

# Gate-log records (printed beside, never asserted by equality).
REC_CTL_MISS = 0.754562            # 20261101T210812Z_WF-1-s27-W1.log (128 MHz against 64 MHz)
REC_EX139_SWEEP_S = 109            # 20261021T020252Z_EX-139.log.gz (whole-example window, s)


def _p0(comm, text):
    if comm.rank == 0:
        print(text, flush=True)


def _boom(*a, **k):
    raise AssertionError("(S0) a loaded-sweep tune reached a mesh build or a solve")


def main() -> None:
    comm = MPI.COMM_WORLD
    started = time.perf_counter()
    assert comm.size == 2, f"asserted at -n 2; got -n {comm.size}"
    OUTPUT_DIR.mkdir(exist_ok=True)
    FIGURE_DIR.mkdir(exist_ok=True)
    if comm.rank == 0:
        if RUN_DIR.exists():
            shutil.rmtree(RUN_DIR)
        RUN_DIR.mkdir(parents=True)
    comm.Barrier()

    # ---- the one sweep, and the file -------------------------------------------------------------
    study0, _phantom = _port_study(1)
    t0 = time.perf_counter()
    sweep = study0.solve_sweep(FREQS, comm)
    t_sweep = time.perf_counter() - t0
    assert tuple(float(f) for f in sweep.frequencies_hz) == FREQS
    ts_path = Path(sweep.write_touchstone(RUN_DIR / "f_small_sweep_saved.s4p"))
    res0 = sweep.results[0]
    _p0(comm, f"[sweep] {res0.cells} cells, solve_sweep {t_sweep:.1f} s (mesh {sweep.mesh_time:.1f}, "
              f"solves {[round(t, 1) for t in sweep.sweep_times]}); wrote {ts_path.name} under "
              f"{RUN_DIR.name}/")

    base = yaml.safe_load(SWEEP_YAML.read_text())
    del base["study"]["sweep_hz"]            # load_result and sweep_hz are mutually exclusive
    ymls = []
    for i, f in enumerate(FREQS):
        m = {**copy.deepcopy(base), "load_result": ts_path.name, "tune": {"k": 1, "at_hz": f}}
        path = RUN_DIR / f"f{i}.yaml"
        if comm.rank == 0:
            path.write_text(yaml.safe_dump(m, sort_keys=False))
        ymls.append(path)
    comm.Barrier()

    # ---- the loaded routes, with a raising wrapper on every mesh and solve entry (S0) -----------
    calls = {"n": 0}
    saved = (Study.solve, Study._solve_ports, FourPortBirdcageGeometry.build, gmsh.initialize)

    def _counting_boom(*a, **k):
        calls["n"] += 1
        return _boom(*a, **k)

    Study.solve = _counting_boom
    Study._solve_ports = _counting_boom
    FourPortBirdcageGeometry.build = _counting_boom
    gmsh.initialize = _counting_boom
    apis, runs = {}, {}
    ref_raised = None
    try:
        t0 = time.perf_counter()
        fresh, _ = _port_study(1)
        fresh.load_result(ts_path)
        for f in FREQS:
            apis[f] = fresh.tune(k=1, comm=comm, at_hz=f)
        t_api = time.perf_counter() - t0
        t0 = time.perf_counter()
        for i, f in enumerate(FREQS):
            runs[f] = run_study_file(ymls[i], output_dir=RUN_DIR / f"out{i}", comm=comm)
        t_yaml = time.perf_counter() - t0
        # the column one frequency along, through the API (c-next)
        nxt = {FREQS[0]: apis[FREQS[1]], FREQS[1]: apis[FREQS[2]],
               FREQS[2]: None}
        # (Ref)
        try:
            fresh.tune(k=1, comm=comm, at_hz=FREQS[2] * 1.5)
        except ValueError as exc:
            ref_raised = str(exc)
    finally:
        Study.solve, Study._solve_ports, FourPortBirdcageGeometry.build, gmsh.initialize = saved
    n_calls = comm.allreduce(calls["n"], op=MPI.SUM)

    # ---- collect every reading first, then assert -------------------------------------------------
    rows = []
    for i, f in enumerate(FREQS):
        api, block = apis[f], runs[f].summary["tune"]
        rec = ITEM722_C_MODE_F[i]
        rel = abs(float(api.c_mode_farad) / rec - 1.0)
        differing = _tune_block_hex_equal(block, api)
        rows.append((f, api, block, rec, rel, differing))
    miss_next = [abs(float(nxt[f].c_mode_farad) / ITEM722_C_MODE_F[i] - 1.0)
                 for i, f in enumerate(FREQS) if nxt[f] is not None]
    miss_6428 = abs(float(apis[FREQS[2]].c_mode_farad) / float(apis[FREQS[1]].c_mode_farad) - 1.0)

    _p0(comm, f"[(S0)] mesh-build / solve / gmsh.initialize calls after the sweep: {n_calls} "
              f"(ASSERTED == 0); API load + 3 tunes {t_api:.3f} s, run_study_file x3 {t_yaml:.3f} s "
              f"against the sweep's {t_sweep:.1f} s")
    for f, api, block, rec, rel, differing in rows:
        _p0(comm, f"[(C=)] {f / 1e6:g} MHz: C_mode {float(api.c_mode_farad):.12e} F, record "
                  f"{rec:.12e} F, rel {rel:.3e} (ASSERTED <= REDUCTION_BAND {REDUCTION_BAND:.0e})")
        _p0(comm, f"[(I=)] {f / 1e6:g} MHz: API {_hex(api.c_mode_farad)} vs YAML "
                  f"{block['c_mode_farad'].hex()}; differing {differing} (ASSERTED []); "
                  f"selected {block['frequency_hz']!r} Hz vs requested {f!r} Hz (ASSERTED equal)")
    _p0(comm, f"[(c-next)] column one frequency along against the record: miss "
              f"{' / '.join(f'{x:.6f}' for x in miss_next)} (ASSERTED each > REDUCTION_BAND "
              f"{REDUCTION_BAND:.0e}); 128 MHz against 64 MHz {miss_6428:.6f} (gate record "
              f"{REC_CTL_MISS})")
    _p0(comm, f"[(Ref)] at_hz {FREQS[2] * 1.5:.3e} Hz through the API: "
              f"{'ValueError' if ref_raised else 'no error'}: {ref_raised} (ASSERTED ValueError)")
    _p0(comm, f"[time] sweep {t_sweep:.1f} s vs load + three tunes {t_api:.3f} s (API) and "
              f"{t_yaml:.3f} s (three run_study_file); EX-139's whole sweep example took "
              f"{REC_EX139_SWEEP_S} s (printed)")

    assert n_calls == 0
    for f, api, block, rec, rel, differing in rows:
        assert differing == [], f
        assert rel <= REDUCTION_BAND, f
        assert api.frequency_hz == f and block["frequency_hz"] == f, f
    assert all(x > REDUCTION_BAND for x in miss_next)
    assert ref_raised is not None and "at_hz" in ref_raised

    # ---- setup figure (no-op unless FEM_EM_SETUP_FIGURES=1) -------------------------------------
    geometry = study0.geometry
    msh, cell_tags = res0.mesh, res0.cell_tags
    local_tags = np.unique(cell_tags.values)
    all_tags = sorted({int(x) for part in comm.allgather(local_tags.tolist()) for x in part})
    phantom_tag = geometry.phantom_tag
    other = tuple(x for x in all_tags if x != phantom_tag)
    names = {phantom_tag: "gel phantom"}
    for x in other:
        names[x] = "air / coil regions"
    write_setup_figure(
        msh, cell_tags, FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names=names, hide_tags=other, translucent_tags=(phantom_tag,),
        slice_normal=(0.0, 0.0, 1.0), slice_origin=(0.0, 0.0, 0.0),
        title="mri:39 -- the one sweep behind the saved file (F-small four-port, degree 1)",
        comm=comm)
    gc.collect()
    comm.Barrier()

    # ---- ParaView: cell tags + sigma (DG0); a loaded result has no field --------------------------
    dg0 = fem.functionspace(msh, ("DG", 0))
    sigma_out = fem.Function(dg0, name="sigma")
    sigma_out.x.array[:] = 0.0
    cells = cell_tags.indices
    in_phantom = cells[cell_tags.values == phantom_tag]
    for c in in_phantom:
        sigma_out.x.array[dg0.dofmap.cell_dofs(int(c))[0]] = SALINE_SIGMA
    sigma_out.x.scatter_forward()
    xdmf_path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_combined", msh, cell_tags, {"sigma": sigma_out}, comm=comm)
    adopt_host_ownership(OUTPUT_DIR, comm=comm)
    _p0(comm, f"[paraview] wrote {xdmf_path} (cell tags + sigma; no field by design)")
    _p0(comm, f"\nAll assertions hold. Total elapsed {time.perf_counter() - started:.1f} s.")


if __name__ == "__main__":
    main()
