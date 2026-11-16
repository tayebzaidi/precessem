"""Example (`EX-160`): tune a saved four-port with no field solve -- ``Study.load_result`` and
``femem run`` on a Touchstone file, and a swapped file that misses.

`WF-16` step 14 (U94) let a Study read a single-frequency four-port ``.s4p`` back into
``Study.result`` (a ``LoadedPortResult``: S only, no mesh, no fields) and `WF-1` step 24 (U99)
put the same read behind ``load_result:`` in a ``femem run`` file.  Every earlier example that
reads a Touchstone does so to compare matrices; this is the first whose tuned capacitance comes
from a file rather than a field solve.

* one F-small four-port solve (degree 1, PEC hole, exact-rectangle ports, 64 MHz, ``-n 2``) and
  ``write_touchstone``;
* (API) a fresh ``Study.load_result(path)`` + ``tune(k=1)``;
* (YAML) ``tests/package/wf1/f_small_hole_pec_tune.yaml`` plus ``load_result:`` through
  ``run_study_file``;
* (CM) both routes' ``C_mode`` equal ``ITEM429_C_MODE_F`` within ``LEG_D0_REPRODUCTION_BAND``;
* (S0) a counter on ``Study.solve`` / ``Study._solve_ports`` reads 0 for the loaded routes;
* (=) the YAML summary's ``c_mode_farad`` is ``float.hex``-equal to the API's;
* (c-sw) the P1 / P2 swapped copy of the file, through both routes, misses by more than
  ``LEG_D0_REPRODUCTION_BAND`` (asserted).

Every record and band is imported by name from the gate modules (the ANS-1 rule); the records
printed beside the readings sit in the example-local ``REC_*`` with their log, never asserted by
equality.  Run it through the example runner::

    ./scripts/run_examples.sh -e mri:35 -n 2 -t 300

Not claimed: F-human or F-clinical loads; fields, SAR or B1+ from a file; an absolute tuning
claim; any external check.
"""

from __future__ import annotations

import copy
import gc
import shutil
import sys
import time
from pathlib import Path

import numpy as np
import ufl
import yaml
from mpi4py import MPI

from dolfinx import fem

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from fem_em_solver import LumpedPortDrive, Study, run_study_file  # noqa: E402
from fem_em_solver import study as study_mod  # noqa: E402
from fem_em_solver.io.paraview_utils import (  # noqa: E402
    adopt_host_ownership,
    write_xdmf_with_tags,
)
from fem_em_solver.ports.sparameters import run_n_port_sparameter_sweep  # noqa: E402
from fem_em_solver.post.setup_figure import write_setup_figure  # noqa: E402

from tests.package.test_wf16_study_load_result import _swap12, _write  # noqa: E402
from tests.package.test_wf16_study_tune_capacitors import ITEM429_C_MODE_F, _study  # noqa: E402
from tests.validation.test_port_birdcage_four_port import LEG_D0_REPRODUCTION_BAND  # noqa: E402

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"
BASENAME = "mri_35_tune_a_saved_result"
RUN_DIR = OUTPUT_DIR / "mri_35_run"
TUNE_YAML = _REPO_ROOT / "tests" / "package" / "wf1" / "f_small_hole_pec_tune.yaml"

# Gate-log records (printed beside, never asserted by equality).
REC_SWAP_MISS = 0.19679            # 20261027T214345Z_WF-1-s24-W1.log (U99)
REC_U99_SOLVE_S = 29.8             # U99 W1, same log


def _p0(comm, text):
    if comm.rank == 0:
        print(text, flush=True)


def _fresh(path, comm):
    if comm.rank == 0 and path.exists():
        shutil.rmtree(path)
    comm.Barrier()


def _write_yaml(path, mapping, comm):
    if comm.rank == 0:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(yaml.safe_dump(mapping, sort_keys=False))
    comm.Barrier()
    return path


def _hex(x):
    return float(x).hex()


def main() -> None:
    comm = MPI.COMM_WORLD
    started = time.perf_counter()
    assert comm.size == 2, f"asserted at -n 2; got -n {comm.size}"
    OUTPUT_DIR.mkdir(exist_ok=True)
    FIGURE_DIR.mkdir(exist_ok=True)
    _fresh(RUN_DIR, comm)

    # ---- the one solve, and the files ----------------------------------------------------------
    study0 = _study(LumpedPortDrive())
    fhz = study0.frequency_hz
    t0 = time.perf_counter()
    res0 = study0.solve(comm)
    t_solve = time.perf_counter() - t0
    s_solved = np.asarray(res0.s_matrix, dtype=np.complex128).copy()
    RUN_DIR.mkdir(exist_ok=True)
    ts_path = Path(res0.write_touchstone(RUN_DIR / "f_small_saved.s4p"))
    swap_path = _write(RUN_DIR / "f_small_swap12.s4p", [(fhz, _swap12(s_solved))], comm)
    base = yaml.safe_load(TUNE_YAML.read_text())
    ok_map = {**copy.deepcopy(base), "load_result": ts_path.name}
    sw_map = {**copy.deepcopy(base), "load_result": swap_path.name}
    yml_ok = _write_yaml(RUN_DIR / "f_small_saved_load.yaml", ok_map, comm)
    yml_sw = _write_yaml(RUN_DIR / "f_small_swap12_load.yaml", sw_map, comm)
    _p0(comm, f"[solve] {res0.cells} cells, {t_solve:.1f} s (U99 W1 solve {REC_U99_SOLVE_S} s); "
              f"wrote {ts_path.name} and {swap_path.name} under {RUN_DIR.name}/")

    # ---- loaded routes, with the solve counter (S0) ---------------------------------------------
    calls = {"n": 0}
    real_solve, real_ports = Study.solve, Study._solve_ports

    def _count(real):
        def _wrapped(*a, **k):
            calls["n"] += 1
            return real(*a, **k)
        return _wrapped

    Study.solve = _count(real_solve)
    Study._solve_ports = _count(real_ports)
    try:
        t0 = time.perf_counter()
        st_api = _study(LumpedPortDrive())
        loaded = st_api.load_result(ts_path)
        api = st_api.tune(k=1, comm=comm)
        t_api = time.perf_counter() - t0
        t0 = time.perf_counter()
        run = run_study_file(yml_ok, output_dir=RUN_DIR / "out_ok", comm=comm)
        t_yaml = time.perf_counter() - t0
        st_sw = _study(LumpedPortDrive())
        st_sw.load_result(swap_path)
        api_sw = st_sw.tune(k=1, comm=comm)
        run_sw = run_study_file(yml_sw, output_dir=RUN_DIR / "out_sw", comm=comm)
    finally:
        Study.solve, Study._solve_ports = real_solve, real_ports
    n_calls = comm.allreduce(calls["n"], op=MPI.SUM)

    # (S0)
    _p0(comm, f"[(S0)] Study.solve / _solve_ports calls during the four loaded runs: {n_calls} "
              f"(ASSERTED == 0); load + tune {t_api:.3f} s (API), run_study_file {t_yaml:.3f} s "
              f"against the solve's {t_solve:.1f} s")
    assert n_calls == 0
    assert isinstance(loaded, study_mod.LoadedPortResult)
    assert isinstance(run.study.result, study_mod.LoadedPortResult)
    rt = comm.allreduce(float(np.max(np.abs(loaded.s_matrix - s_solved))), op=MPI.MAX)
    _p0(comm, f"[RT] max|S_loaded - S_solved| = {rt:.3e} (printed)")

    # (CM)
    c_api = float(api.c_mode_farad)
    c_yaml = float(run.summary["tune"]["c_mode_farad"])
    rel_api = abs(c_api / ITEM429_C_MODE_F - 1.0)
    rel_yaml = abs(c_yaml / ITEM429_C_MODE_F - 1.0)
    d9_api = f"{c_api:.9e}" == f"{ITEM429_C_MODE_F:.9e}"
    d9_yaml = f"{c_yaml:.9e}" == f"{ITEM429_C_MODE_F:.9e}"
    _p0(comm, f"[(CM)] C_mode API {c_api:.9e} F, YAML {c_yaml:.9e} F, item 429 "
              f"{ITEM429_C_MODE_F:.9e} F; rel {rel_api:.3e} / {rel_yaml:.3e} "
              f"(ASSERTED <= {LEG_D0_REPRODUCTION_BAND:.0e}); .9e digits equal "
              f"{d9_api} / {d9_yaml} (printed)")
    assert rel_api <= LEG_D0_REPRODUCTION_BAND
    assert rel_yaml <= LEG_D0_REPRODUCTION_BAND

    # (=)
    block = run.summary["tune"]
    differing = []
    if block["c_mode_farad"].hex() != _hex(api.c_mode_farad):
        differing.append("c_mode_farad")
    if block["off_diagonal_residual"].hex() != _hex(api.off_diagonal_residual):
        differing.append("off_diagonal_residual")
    _p0(comm, f"[(=)] YAML summary c_mode_farad {block['c_mode_farad'].hex()} vs API "
              f"{_hex(api.c_mode_farad)}; residual {block['off_diagonal_residual']:.3e}; "
              f"differing {differing} (ASSERTED [])")
    assert differing == []

    _p0(comm, "[z_k] " + ", ".join(f"{complex(z).real:+.6e}{complex(z).imag:+.6e}j"
                                   for z in api.z_k) + "  (printed)")

    # (c-sw)
    c_sw_api = float(api_sw.c_mode_farad)
    c_sw_yaml = float(run_sw.summary["tune"]["c_mode_farad"])
    miss_api = abs(c_sw_api / ITEM429_C_MODE_F - 1.0)
    miss_yaml = abs(c_sw_yaml / ITEM429_C_MODE_F - 1.0)
    _p0(comm, f"[(c-sw)] P1/P2 swapped file: C_mode API {c_sw_api:.9e} F (miss {miss_api:.5f}), "
              f"YAML {c_sw_yaml:.9e} F (miss {miss_yaml:.5f}); record {REC_SWAP_MISS}; ASSERTED > "
              f"{LEG_D0_REPRODUCTION_BAND:.0e} ({miss_api / LEG_D0_REPRODUCTION_BAND:.2e} x the bar)")
    assert miss_api > LEG_D0_REPRODUCTION_BAND
    assert miss_yaml > LEG_D0_REPRODUCTION_BAND

    # ---- setup figure (no-op unless FEM_EM_SETUP_FIGURES=1) -------------------------------------
    geometry = study0.geometry
    msh, cell_tags = res0.mesh, res0.cell_tags
    local_tags = np.unique(cell_tags.values)
    all_tags = sorted({int(x) for part in comm.allgather(local_tags.tolist()) for x in part})
    phantom_tag = geometry.phantom_tag
    other = tuple(x for x in all_tags if x != phantom_tag)
    names = {phantom_tag: "gel phantom"}
    for x in other:
        names[x] = "air / PEC-proxy hole region"
    write_setup_figure(
        msh, cell_tags, FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names=names, hide_tags=other, translucent_tags=(phantom_tag,),
        slice_normal=(0.0, 0.0, 1.0), slice_origin=(0.0, 0.0, 0.0),
        title="mri:35 -- the one solve behind the saved four-port (F-small, PEC hole, degree 1)",
        comm=comm)
    gc.collect()
    comm.Barrier()

    # ---- ParaView: the solved Study's P1-drive |E| and sigma (a loaded result has no field) ----
    t0 = time.perf_counter()
    sweep = run_n_port_sparameter_sweep(
        res0.problem, res0.port_definitions, lumped_sheet_ports=res0.specs,
        lumped_sheet_facet_tags=res0.facet_tags, degree=res0.degree, keep_fields=True)
    fields = sweep.fields["P1"]
    e = fields.e_complex
    p1_space = fem.functionspace(msh, ("Lagrange", 1))
    e_abs = fem.Function(p1_space, name="E_abs")
    e_abs.interpolate(fem.Expression(ufl.sqrt(ufl.real(ufl.inner(e, e))),
                                     p1_space.element.interpolation_points))
    e_abs.x.array[:] = np.real(e_abs.x.array)
    e_abs.x.scatter_forward()
    dg0 = fem.functionspace(msh, ("DG", 0))
    sigma_out = fem.Function(dg0, name="sigma")
    sigma_out.interpolate(fields.sigma_field)
    sigma_out.x.array[:] = np.real(sigma_out.x.array)
    sigma_out.x.scatter_forward()
    xdmf_path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_combined", msh, cell_tags,
        {"E_abs": e_abs, "sigma": sigma_out}, comm=comm)
    adopt_host_ownership(OUTPUT_DIR, comm=comm)
    t_x = time.perf_counter() - t0
    rel_x = float(np.linalg.norm(np.asarray(sweep.s_matrix) - s_solved) / np.linalg.norm(s_solved))
    _p0(comm, f"[paraview] wrote {xdmf_path}; field sweep {t_x:.1f} s; keep_fields S vs Study S "
              f"{rel_x:.3e} (printed)")
    _p0(comm, f"\nAll assertions hold. Total elapsed {time.perf_counter() - started:.1f} s.")


if __name__ == "__main__":
    main()
