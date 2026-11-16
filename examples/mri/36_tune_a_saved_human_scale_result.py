"""Example (`EX-163`): tune a saved human-scale birdcage with no second solve -- the 32-port copper
F-human through ``Study.load_result``.

`WF-16` step 15 (U100) let a Study read a single-frequency 32-port ``.s32p`` back into
``Study.result`` and ring-tune it.  ``mri:35`` does this for a four-port; this is the same
capability at human scale (32 ring-gap ports, Leontovich copper, ``n_slots=16``), where one solve
costs minutes at 8 ranks and the load costs seconds.

* one copper F-human solve (degree 1, hole mesh, 64 MHz, kappa off, ``-n 8``), the in-memory
  ``tune(1, n_slots=16)`` and ``write_touchstone``;
* a fresh Study: ``load_result(path)`` + ``tune(1, n_slots=16)``;
* (CM) the loaded ``C_mode`` equals ``WF12_STEP3_C_MODE_CU_F`` within ``REDUCTION_BAND``;
* (T=) the loaded ``C_mode`` equals the in-memory tune's within ``REDUCTION_BAND``;
* (S0) class-level spies on ``Study.solve`` / ``Study._solve_ports`` read 0 calls;
* (c-sw) the file with ports 1 and 2 swapped in the data (names kept), loaded into a third fresh
  Study, misses ``WF12_STEP3_C_MODE_CU_F`` by more than ``LEG_D0_REPRODUCTION_BAND`` (asserted).

Every record and band is imported by name from the gate modules (the ANS-1 rule).  Run it through
the example runner::

    ./scripts/run_examples.sh -e mri:36 -n 8 -t 540

Not claimed: a 64-port (F-clinical) load, a field from a loaded file, degree 2, any resonance /
Q / match statement.
"""

from __future__ import annotations

import gc
import sys
import time
from pathlib import Path

import numpy as np
from mpi4py import MPI

from dolfinx import fem

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from fem_em_solver import Study  # noqa: E402
from fem_em_solver import study as study_mod  # noqa: E402
from fem_em_solver.io.paraview_utils import (  # noqa: E402
    adopt_host_ownership,
    write_xdmf_with_tags,
)
from fem_em_solver.post.setup_figure import write_setup_figure  # noqa: E402

from tests.package.test_wf16_load_result_f_human import (  # noqa: E402
    N_SLOTS,
    _swap12_n,
    _write_named,
)
from tests.package.test_wf16_study_tune_f_human_copper import (  # noqa: E402
    RECORD_RANK_WIDTH,
    WF12_STEP3_C_MODE_CU_F,
    _study,
)
from tests.validation.test_lossy_sphere_fullwave import SALINE_SIGMA  # noqa: E402
from tests.validation.test_port_birdcage_four_port import LEG_D0_REPRODUCTION_BAND  # noqa: E402
from tests.validation.test_port_lumped_rlc_termination import REDUCTION_BAND  # noqa: E402

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"
BASENAME = "mri_36_tune_a_saved_human_scale_result"
RUN_DIR = OUTPUT_DIR / "mri_36_run"

# Gate-log records (printed beside, never asserted by equality).
REC_SWAP_MISS = 0.4254907          # 20261028T004356Z_WF-16-s15-W1.log:7356 (U100)
REC_U100_SOLVE_S = 227.9           # same log, :7348


def _p0(comm, text):
    if comm.rank == 0:
        print(text, flush=True)


def main() -> None:
    comm = MPI.COMM_WORLD
    started = time.perf_counter()
    assert comm.size == RECORD_RANK_WIDTH, (comm.size, RECORD_RANK_WIDTH)
    OUTPUT_DIR.mkdir(exist_ok=True)
    FIGURE_DIR.mkdir(exist_ok=True)
    if comm.rank == 0:
        RUN_DIR.mkdir(exist_ok=True)
    comm.Barrier()

    # ---- the one solve, the in-memory tune, and the file ------------------------------------------
    study0 = _study()
    fhz = study0.frequency_hz
    t0 = time.perf_counter()
    res0 = study0.solve(comm)
    t_solve = time.perf_counter() - t0
    tune_mem = study0.tune(1, n_slots=N_SLOTS, comm=comm)
    ts_path = Path(res0.write_touchstone(RUN_DIR / f"f_human_cu_d1_n{comm.size}.s32p"))
    _p0(comm, f"[solve] {res0.cells} cells, solve {t_solve:.1f} s (U100 solve {REC_U100_SOLVE_S} s), "
              f"-n {comm.size}; wrote {ts_path.name} under {RUN_DIR.name}/")

    # ---- loaded Study, with class-level spies (S0) ------------------------------------------------
    calls = {"solve": 0, "_solve_ports": 0}

    def _spy(name):
        def _boom(self, *a, **k):
            calls[name] += 1
            raise AssertionError(f"(S0) the loaded Study called {name}")
        return _boom

    real_solve, real_ports = Study.solve, Study._solve_ports
    Study.solve, Study._solve_ports = _spy("solve"), _spy("_solve_ports")
    try:
        t0 = time.perf_counter()
        fresh = _study()
        loaded = fresh.load_result(ts_path)
        assert isinstance(loaded, study_mod.LoadedPortResult)
        assert len(loaded.port_ids) == 32
        tune_ld = fresh.tune(1, n_slots=N_SLOTS, comm=comm)
        t_load = time.perf_counter() - t0

        sw_path = _write_named(RUN_DIR / f"f_human_swap12_n{comm.size}.s32p", fhz,
                               _swap12_n(loaded.s_matrix), list(loaded.port_ids), comm)
        st_sw = _study()
        st_sw.load_result(sw_path)
        try:
            c_sw = st_sw.tune(1, n_slots=N_SLOTS, comm=comm).c_mode_farad
            miss = abs(c_sw / WF12_STEP3_C_MODE_CU_F - 1.0)
            sw_note = f"C {c_sw:.12e} F"
        except ValueError as exc:  # a non-capacitive read is a miss too; printed
            miss = float("inf")
            sw_note = f"tune raised {exc}"
    finally:
        Study.solve, Study._solve_ports = real_solve, real_ports

    _p0(comm, f"[(S0)] class-level spies: solve {calls['solve']}, _solve_ports "
              f"{calls['_solve_ports']} (ASSERTED 0 / 0); load + tune {t_load:.3f} s against the "
              f"solve's {t_solve:.1f} s")
    assert calls == {"solve": 0, "_solve_ports": 0}

    cl, cm = float(tune_ld.c_mode_farad), float(tune_mem.c_mode_farad)
    g_cm = comm.allreduce(abs(cl / WF12_STEP3_C_MODE_CU_F - 1.0), op=MPI.MAX)
    g_t = comm.allreduce(abs(cl / cm - 1.0), op=MPI.MAX)
    _p0(comm, f"[(CM)] loaded C_mode {cl:.12e} F vs WF12_STEP3_C_MODE_CU_F "
              f"{WF12_STEP3_C_MODE_CU_F:.12e} F: rel {g_cm:.3e} (ASSERTED <= REDUCTION_BAND "
              f"{REDUCTION_BAND:.1e})")
    _p0(comm, f"[(T=)] loaded {cl:.12e} F vs in-memory {cm:.12e} F: rel {g_t:.3e} (ASSERTED <= "
              f"REDUCTION_BAND {REDUCTION_BAND:.1e}); off-block residual loaded "
              f"{tune_ld.off_diagonal_residual:.6e} vs in-memory "
              f"{tune_mem.off_diagonal_residual:.6e} (printed)")
    _p0(comm, f"[ports] file order ({len(loaded.port_ids)}): {list(loaded.port_ids)}")
    assert g_cm <= REDUCTION_BAND
    assert g_t <= REDUCTION_BAND

    miss_g = comm.allreduce(float(miss), op=MPI.MAX)
    _p0(comm, f"[(c-sw)] P1<->P2 swapped file: {sw_note}, miss {miss_g:.6e} (record "
              f"{REC_SWAP_MISS}); ASSERTED > LEG_D0_REPRODUCTION_BAND "
              f"{LEG_D0_REPRODUCTION_BAND:.1e}")
    assert miss_g > LEG_D0_REPRODUCTION_BAND

    # ---- setup figure (no-op unless FEM_EM_SETUP_FIGURES=1) ---------------------------------------
    geometry = study0.geometry
    msh, cell_tags = res0.mesh, res0.cell_tags
    local_tags = np.unique(cell_tags.values)
    all_tags = sorted({int(x) for part in comm.allgather(local_tags.tolist()) for x in part})
    phantom_tag = geometry.phantom_tag
    other = tuple(x for x in all_tags if x != phantom_tag)
    names = {phantom_tag: "gel phantom"}
    for x in other:
        names[x] = "air / copper-hole regions"
    write_setup_figure(
        msh, cell_tags, FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names=names, hide_tags=other, translucent_tags=(phantom_tag,),
        slice_normal=(0.0, 0.0, 1.0), slice_origin=(0.0, 0.0, 0.0),
        title="mri:36 -- the one solve behind the saved 32-port (copper F-human, degree 1)",
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
