"""Example (`EX-107`): figures that equal their records -- write the records, draw
them, and prove the drawing is the record.

`WF-8` (``post.figures`` and ``femem figure``) draws three kinds of record a
``Study`` user has: an S-parameter Touchstone file, a cell-centred field in an
XDMF/HDF5 bundle, and a probe-series CSV.  This example regenerates all three
sources in-process, under its own file stem (no example reads another's output),
through the same gated fixtures `mri:11`, `mri:10` and `mat:13` call:

* (a) S: ``_port_study(1)`` + ``solve_sweep`` at ``FREQUENCIES`` +
  ``write_touchstone`` (`WF-16`);
* (b) SAR: ``_study(degree=1)`` + ``temperature(RobinCooling(...))`` and the
  ``theta`` / ``SAR`` / ``q`` fields written by ``write_xdmf_with_tags``
  (`WF-15` / `MAT-10`);
* (c) a probe series CSV from ``MAT-19``'s ``_build_run()`` (the lifted fixture
  `mat:13` calls; ``_check_probe_series`` is never called).

It then draws each through ``plot_sparameters`` / ``plot_xdmf_slice`` /
``plot_probe_series``, reads every plotted array back from the figure object and
asserts it bitwise equal (float64, ``np.array_equal``) to the record parsed
independently from disk by `WF-8`'s own source parsers (imported, now taking a
path keyword).

Every helper comes from ``post.figures`` or from the gate modules by import (the
ANS-1 rule).  Cited inputs, not records or bands (the `EX-101` ruling (3)
precedent): the plant indices ``[1, 0, 0]`` (S) and ``5`` (probe), the slice
plant rule "first non-zero selected value", the factor ``1 + 1e-6`` and the wrong
field ``"real_q"`` -- each as in ``tests/post/test_wf8_figures.py``.

Asserted (all exact, no tolerance exists):

* (S)   16 / 16 (i, j) pairs: ``|S|`` and ``arg S`` lines equal the source, with
        4 ports and 3 frequencies;
* (X)   slice colours = ``val[idx]`` and slice points = ``cent[idx][:, :2]``,
        ``len(idx) > 0``;
* (P)   probe line x / y equal the CSV parse;
* (F)   the three PNGs exist and are non-empty on rank 0; every other rank's
        ``plot_*`` returned ``None``;
* (CLI) ``femem figure sparams|slice|probe`` returns ``[0, 0, 0]`` on every rank
        and rank 0's three PNGs are non-empty.

Negative controls (asserted):

* (c1)  a planted (1 + 1e-6) change is caught at exactly one index per figure:
        ``[1]`` (S), ``[k]`` (slice), ``[5]`` (probe);
* (c2)  a wrong column (S12 under the S11 label, ``real_q`` for ``real_SAR``,
        ``probe_c1_K`` for ``probe_c0_K``) gives ``mismatch_indices(...).size > 0``.

Printed, not asserted: each source's path, bytes and sha256; frequencies; cells
and slab count; rows; the (c1) / (c2) indices and counts; PNG sizes; wall time
per source build; the display slab's owned-cell count.

Run it through the example runner (the ``mri:`` group sources the complex build)::

    ./scripts/run_examples.sh -e mri:15 -n 2 -t 300

Not claimed: that the sources are correct (they are records of other gates), any
S, SAR or temperature accuracy, a ``|B1+|`` or interpolated slice, style, or
ParaView opening the XDMF.
"""

from __future__ import annotations

import csv
import gc
import hashlib
import os
import sys
import time
from pathlib import Path

import numpy as np
from mpi4py import MPI

from dolfinx import default_scalar_type, fem
from dolfinx.mesh import meshtags

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from fem_em_solver import RobinCooling  # noqa: E402
from fem_em_solver.cli import main as femem_main  # noqa: E402
from fem_em_solver.io.paraview_utils import (  # noqa: E402
    adopt_host_ownership,
    write_xdmf_with_tags,
)
from fem_em_solver.post.figures import (  # noqa: E402
    draw_series,
    draw_slice,
    draw_sparameters,
    mismatch_indices,
    plot_probe_series,
    plot_sparameters,
    plot_xdmf_slice,
    select_slab,
)
from fem_em_solver.post.setup_figure import write_setup_figure  # noqa: E402

from tests.package.test_wf15_study_lossy_sphere import _study  # noqa: E402
from tests.package.test_wf16_study_sweep import _port_study  # noqa: E402
from tests.post.test_wf8_figures import (  # noqa: E402
    COLUMN,
    FIELD,
    WRONG_COLUMN,
    _coll_values,
    _coll_xy,
    _exact,
    _line_xy,
    _source_csv,
    _source_s,
    _source_slice,
)
from tests.validation.test_ans4_exact_strip_records import (  # noqa: E402
    ANS4_EXACT_STRIP_RECORD_RANK_WIDTH,
    FREQUENCIES,
)
from tests.validation.test_mat10_sar_to_heat_sphere import (  # noqa: E402
    G,
    H,
    K,
    T_INF,
)
from tests.validation.test_mat19_virtual_temperature_probe import (  # noqa: E402
    DT,
    N_MAIN,
    _build_run,
)
from tests.validation.test_lossy_sphere_sar import RHO_KG_M3  # noqa: E402

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"
BASENAME = "mri_15_figures_from_records"
DEGREE = 1  # the records' degree; the Study default is 2.

# Cited inputs (not records or bands), tests/post/test_wf8_figures.py:
PLANT_S = (1, 0, 0)  # :143  mag[1, 0, 0]
PLANT_PROBE = 5  # :223  y[5]
PLANT_FACTOR = 1.0 + 1e-6  # :143, :184, :223
WRONG_FIELD = "real_q"  # :193


def _all(comm, flag) -> bool:
    """Rank 0's verdict, broadcast so every rank asserts the same thing."""
    return bool(comm.bcast(bool(flag), root=0))


def _sha(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _say(comm, msg) -> None:
    if comm.rank == 0:
        print(msg, flush=True)


def main() -> None:
    comm = MPI.COMM_WORLD
    started = time.perf_counter()

    if not np.issubdtype(np.dtype(default_scalar_type), np.complexfloating):
        raise RuntimeError(
            "this example needs the complex DolfinX build "
            "(source /usr/local/bin/dolfinx-complex-mode); the runner's `mri:` "
            "group sources it automatically"
        )
    assert comm.size == ANS4_EXACT_STRIP_RECORD_RANK_WIDTH, (
        f"the S source is asserted at -n {ANS4_EXACT_STRIP_RECORD_RANK_WIDTH}; got -n {comm.size}"
    )
    OUTPUT_DIR.mkdir(exist_ok=True)

    _say(comm, "=" * 78)
    _say(comm, "EX-107 -- figures that equal their records (WF-8 on this example's own sources)")
    _say(comm, "=" * 78)

    # ---- (a) S: a three-frequency .s4p -----------------------------------------
    s4p = OUTPUT_DIR / f"{BASENAME}.s4p"
    t0 = time.perf_counter()
    study_s, _phantom = _port_study(DEGREE)
    freqs = tuple(FREQUENCIES.values())
    sw = study_s.solve_sweep(freqs, comm)
    sw.write_touchstone(s4p)
    t_a = time.perf_counter() - t0
    del sw, study_s, _phantom
    gc.collect()

    # ---- (b) SAR: a combined XDMF on the sphere submesh ------------------------
    t0 = time.perf_counter()
    study = _study(degree=DEGREE)
    res = study.solve(comm)
    th = study.temperature(RobinCooling(K, H, T_INF, source_scale=G))
    sub = th.submesh
    theta_out = fem.Function(th.temperature.function_space, name="theta")
    theta_out.x.array[:] = np.real(th.temperature.x.array) - T_INF
    theta_out.x.scatter_forward()
    q_out = fem.Function(th.heat_source.function_space, name="q")
    q_out.x.array[:] = np.real(th.heat_source.x.array)
    q_out.x.scatter_forward()
    sar_out = fem.Function(th.heat_source.function_space, name="SAR")
    sar_out.x.array[:] = np.real(th.heat_source.x.array) / (G * RHO_KG_M3)
    sar_out.x.scatter_forward()
    write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_sar_combined",
        sub,
        None,
        {"theta": theta_out, "SAR": sar_out, "q": q_out},
        comm=comm,
    )  # returns (None, None) off rank 0, so the paths are named here
    xdmf = OUTPUT_DIR / f"{BASENAME}_sar_combined.xdmf"
    h5 = OUTPUT_DIR / f"{BASENAME}_sar_combined.h5"
    t_b = time.perf_counter() - t0

    # ---- setup figure (display-only tags; never asserted on) --------------------
    # Tag 2: cells whose centroid lies in the slab the slice figure draws;
    # tag 1: the rest.  select_slab's own defaults (mid-plane, 5 % of the range),
    # made global.  Every allreduce is outside any rank branch.
    dim = sub.topology.dim
    n_local = sub.topology.index_map(dim).size_local
    n_all = n_local + sub.topology.index_map(dim).num_ghosts
    cen = sub.geometry.x[sub.geometry.dofmap].mean(axis=1)
    z_own = cen[:n_local, 2]
    zmin = comm.allreduce(float(z_own.min()), op=MPI.MIN)
    zmax = comm.allreduce(float(z_own.max()), op=MPI.MAX)
    centre = 0.5 * (zmin + zmax)
    half = 0.05 * (zmax - zmin)
    slab_local = select_slab(cen[:n_all], axis="z", center=centre, half_width=half)
    vals = np.ones(n_all, dtype=np.int32)
    vals[slab_local] = 2
    display_tags = meshtags(sub, dim, np.arange(n_all, dtype=np.int32), vals)
    n_display = comm.allreduce(int(np.count_nonzero(vals[:n_local] == 2)), op=MPI.SUM)
    write_setup_figure(
        sub,
        display_tags,
        FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names={1: "gel sphere (SAR source)", 2: "slab drawn by the slice figure"},
        translucent_tags=(1,),
        slice_normal=(0.0, 1.0, 0.0),
        slice_origin=(0.0, 0.0, 0.0),
        title="mri:15 -- the slice figure is a selection: sphere cells on the z mid-plane slab (WF-8)",
        comm=comm,
    )
    del th, theta_out, q_out, sar_out, res, study
    gc.collect()

    # ---- (c) a probe series CSV -------------------------------------------------
    probe_csv = OUTPUT_DIR / f"{BASENAME}_probe_series.csv"
    t0 = time.perf_counter()
    run = _build_run()
    if comm.rank == 0:
        with probe_csv.open("w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["n", "t_s", "probe_c0_K", "probe_c1_K", "probe_c2_K"])
            for n in range(N_MAIN + 1):
                w.writerow(
                    [n, f"{n * DT:.6f}", *(f"{np.real(run['s_cf'].values[n, j]):.9e}" for j in range(3))]
                )
    comm.Barrier()
    t_c = time.perf_counter() - t0
    del run
    gc.collect()
    adopt_host_ownership(OUTPUT_DIR, comm=comm)

    # ---- the source audit trail (rule (i)) --------------------------------------
    _say(comm, f"[source] build wall: S {t_a:.1f} s, SAR {t_b:.1f} s, probe {t_c:.1f} s")
    if comm.rank == 0:
        for p in (s4p, xdmf, h5, probe_csv):
            print(
                f"[source] {p.relative_to(_REPO_ROOT)} bytes={p.stat().st_size} sha256={_sha(p)}",
                flush=True,
            )
    _say(comm, "")

    out_s = OUTPUT_DIR / f"{BASENAME}_sparams.png"
    out_x = OUTPUT_DIR / f"{BASENAME}_slice.png"
    out_p = OUTPUT_DIR / f"{BASENAME}_probe.png"

    # ---- the three plot_* calls ---------------------------------------------------
    t0 = time.perf_counter()
    res_s = plot_sparameters(s4p, out_s, comm=comm)
    res_x = plot_xdmf_slice(xdmf, FIELD, out_x, axis="z", comm=comm)
    res_p = plot_probe_series(probe_csv, COLUMN, out_p, comm=comm)
    t_plot = time.perf_counter() - t0
    none_ok = comm.rank == 0 or (res_s is None and res_x is None and res_p is None)
    none_ok = bool(comm.allreduce(int(none_ok), op=MPI.MIN))

    # ---- (S), (X), (P): read back and compare (rank 0), verdicts broadcast ---------
    ok_s = ok_x = ok_p = True
    pairs = 0
    n_ports = n_freq = n_cells = n_sel = n_rows = 0
    first = []
    if comm.rank == 0:
        f, mag, deg = _source_s(s4p)
        n_ports, n_freq = mag.shape[1], f.shape[0]
        for i in range(n_ports):
            for j in range(n_ports):
                xm, ym = _line_xy(res_s.lines[(i, j, "mag")])
                xp, yp = _line_xy(res_s.lines[(i, j, "deg")])
                good = _exact(xm, f) and _exact(ym, mag[:, i, j]) and _exact(xp, f) and _exact(yp, deg[:, i, j])
                pairs += int(good)
        ok_s = pairs == 16 and n_ports == 4 and n_freq == 3
        print(f"[(S)] frequencies={f.tolist()} ports={n_ports}; bitwise-equal (i,j) pairs {pairs}/16", flush=True)

        cent, val = _source_slice(h5, FIELD)
        idx = res_x.indices
        n_cells, n_sel, first = int(val.shape[0]), int(len(idx)), idx[:5].tolist()
        ok_x = (
            len(idx) > 0
            and _exact(_coll_values(res_x.collection), val[idx])
            and _exact(_coll_xy(res_x.collection), cent[idx][:, :2])
        )
        print(
            f"[(X)] cells={n_cells} selected len(idx)={n_sel} first={first} "
            f"(WF-8: 19496 / 2772, first [137, 228, 258, 263, 423]); display slab owned cells "
            f"{n_display} (predicted = len(idx))",
            flush=True,
        )

        t_src, y_src = _source_csv(COLUMN, probe_csv)
        n_rows = int(t_src.shape[0])
        x_pl, y_pl = _line_xy(res_p.lines[COLUMN])
        ok_p = _exact(x_pl, t_src) and _exact(y_pl, y_src) and n_rows == N_MAIN + 1
        print(f"[(P)] rows={n_rows} (N_MAIN + 1 = {N_MAIN + 1})", flush=True)

    # ---- (c1) / (c2) controls -------------------------------------------------------
    ok_c1 = ok_c2 = True
    if comm.rank == 0:
        # S
        f, mag, deg = _source_s(s4p)
        mag_p = mag.copy()
        mag_p[PLANT_S] *= PLANT_FACTOR
        r = draw_sparameters(f, mag_p, deg)
        _, y = _line_xy(r.lines[(0, 0, "mag")])
        bad_s = mismatch_indices(y, mag[:, 0, 0])
        r = draw_sparameters(f, mag, deg)
        _, y = _line_xy(r.lines[(0, 1, "mag")])
        bad2_s = mismatch_indices(y, mag[:, 0, 0])
        # slice
        cent, val = _source_slice(h5, FIELD)
        idx = select_slab(cent, axis="z")
        k = int(np.flatnonzero(val[idx] != 0.0)[0])  # first non-zero selected value (test :182)
        val_p = val.copy()
        val_p[idx[k]] *= PLANT_FACTOR
        r = draw_slice(cent, val_p, idx)
        bad_x = mismatch_indices(_coll_values(r.collection), val[idx])
        other = _source_slice(h5, WRONG_FIELD)[1]
        r = draw_slice(cent, other, idx)
        bad2_x = mismatch_indices(_coll_values(r.collection), val[idx])
        # probe
        t_src, y_src = _source_csv(COLUMN, probe_csv)
        y_pp = y_src.copy()
        y_pp[PLANT_PROBE] *= PLANT_FACTOR
        r = draw_series(t_src, y_pp, label=COLUMN)
        _, y = _line_xy(r.lines[COLUMN])
        bad_p = mismatch_indices(y, y_src)
        t_w, y_w = _source_csv(WRONG_COLUMN, probe_csv)
        r = draw_series(t_w, y_w, label=COLUMN)
        _, y = _line_xy(r.lines[COLUMN])
        bad2_p = mismatch_indices(y, y_src)
        print(
            f"[(c1)] planted (1+1e-6) mismatch indices: S {bad_s.tolist()} (WF-8 [1]); slice "
            f"{bad_x.tolist()} (planted at {k}; WF-8 [0]); probe {bad_p.tolist()} (WF-8 [5])",
            flush=True,
        )
        print(
            f"[(c2)] wrong column mismatching: S12 vs S11 {bad2_s.size} (WF-8 3); {WRONG_FIELD} vs "
            f"{FIELD} {bad2_x.size} (WF-8 2772); {WRONG_COLUMN} vs {COLUMN} {bad2_p.size} (WF-8 40)",
            flush=True,
        )
        ok_c1 = bad_s.tolist() == [1] and bad_x.tolist() == [k] and bad_p.tolist() == [PLANT_PROBE]
        ok_c2 = bad2_s.size > 0 and bad2_x.size > 0 and bad2_p.size > 0

    # ---- (F) and (CLI) ------------------------------------------------------------------
    ok_f = none_ok
    if comm.rank == 0:
        sizes = {p.name: p.stat().st_size for p in (out_s, out_x, out_p)}
        ok_f = ok_f and all(v > 0 for v in sizes.values())
        print(f"[(F)] files={sizes}; other ranks returned None: {none_ok}", flush=True)
    cli_out = {k_: OUTPUT_DIR / f"{BASENAME}_cli_{k_}.png" for k_ in ("sparams", "slice", "probe")}
    cases = [
        ("sparams", s4p, []),
        ("slice", xdmf, ["--name", FIELD]),
        ("probe", probe_csv, ["--column", COLUMN]),
    ]
    rcs = [femem_main(["figure", kind, str(src), str(cli_out[kind]), *extra]) for kind, src, extra in cases]
    ok_cli = rcs == [0, 0, 0]
    if comm.rank == 0:
        cli_sizes = {k_: v.stat().st_size for k_, v in cli_out.items()}
        ok_cli = ok_cli and all(v > 0 for v in cli_sizes.values())
        print(f"[(CLI)] rc={rcs} sizes={cli_sizes}", flush=True)
    adopt_host_ownership(OUTPUT_DIR, comm=comm)

    verdicts = {
        "S": _all(comm, ok_s),
        "X": _all(comm, ok_x),
        "P": _all(comm, ok_p),
        "F": _all(comm, ok_f),
        "CLI": _all(comm, ok_cli),
        "c1": _all(comm, ok_c1),
        "c2": _all(comm, ok_c2),
    }
    _say(comm, f"[verdicts] {verdicts}; plot_* wall {t_plot:.2f} s")
    for name, v in verdicts.items():
        assert v, f"({name}) failed"

    _say(
        comm,
        f"\nAll anchors and controls hold. Total elapsed {time.perf_counter() - started:.1f} s.",
    )


if __name__ == "__main__":
    main()
