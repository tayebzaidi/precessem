"""Example (`EX-138`): ``femem run`` on the four-port F-small EM file with point probes --
the per-port point SAR and ``|E|`` table written by the CLI, read back from disk, tied to an
independent Python Study, and a 1 mm displaced probe that misses.

`WF-1` step 20 (item 564) let ``outputs.probes: {points: [[x, y, z], ...]}`` ride on a
``four_port_birdcage`` file: ``run_study_file`` solves the four single-port 1 V drives once,
evaluates `POST-8`'s ``point_sar`` (owner-rule default) at the listed points on the Study's own
kept fields, and writes ``probes.{points, owner, sigma_s_per_m, rho_kg_m3, values}`` into the
run's JSON summary.  Every earlier probe example (``mri:29``-``mri:32``) is a thermal ball
average; none ran point probes on an EM file.  This one runs the tracked
``tests/package/wf1/f_small_em_probes.yaml`` (degree 1, PEC, 64 MHz, exact-rectangle ports,
``-n 2``) through ``cli.main(["run", ...])`` with a call-through wrapper on
``fem_em_solver.study_yaml.run_study_file`` (restored in ``finally``):

* (B) complex build, ``-n 2``; (K) the file's Study equals ``f_small_four_port.yaml``'s, the
  parsed points are hex-equal to the gate's, ``femem check`` rc (0, 0); (X-out) two probes
  outside the phantom are refused naming their index, before any mesh;
* (R) cells and S tie to item 192's records; (U') the summary is ``point_sar`` on the run's
  own fields, bitwise across ranks, JSON == memory, ``|E|`` identity; (J-disk) the probe table
  read from disk alone carries the gated values; (I-Py) file = Python Study at ``L0_BAND``;
* (c-disp) a point displaced 1 mm radially reads differently (asserted).

Every record, band and helper is imported from the gate module (the ANS-1 rule); record
strings sit only in the example-local ``REC`` with their ``log:line`` and are PRINTED beside
the readings, never asserted by equality.  Values print at <= 6 significant digits except the
gate's own ``.6e`` per-port rows.

Run it through the example runner (the ``mri:`` group sources the complex build)::

    ./scripts/run_examples.sh -e mri:33 -n 2 -t 300

Not claimed: absolute, quadrature-drive or compliance SAR (per-1 V single-port values only);
a combined-drive SAR rebuilt from the table; C4 symmetry of the table; probes on F-human,
copper or capacitor files; degree 2; other widths; a peak search or 10 g averaging; the
display fields' values at the probes; any external check.
"""

from __future__ import annotations

import contextlib
import csv
import gc
import io
import json
import sys
import time
from pathlib import Path

import numpy as np
import ufl
from mpi4py import MPI

from dolfinx import default_scalar_type, fem
from dolfinx import mesh as dmesh

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

import fem_em_solver.study_yaml as study_yaml  # noqa: E402
from fem_em_solver.io.paraview_utils import (  # noqa: E402
    adopt_host_ownership,
    write_xdmf_with_tags,
)
from fem_em_solver.post.setup_figure import write_setup_figure  # noqa: E402
from fem_em_solver.utils.instrumentation import report_peak_rss  # noqa: E402

import tests.package.test_wf1_yaml_em_probes as t  # noqa: E402
from tests.validation.test_port_birdcage_lumped_column import PHANTOM_CELL_TAG  # noqa: E402

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"
BASENAME = "mri_33_femem_run_em_file_point_probes"
RUN_DIR = OUTPUT_DIR / "mri_33_run_em_probes"
CSV_PATH = OUTPUT_DIR / f"{BASENAME}_probe_table.csv"
SITE_RADIUS_M = 2.0e-3

# Gate-log records, docs/testing/logs/20261020T123448Z_WF-1-s20-W1-em-probes-n2.log
# (complex, -n 2): per-port rows at .6e, :3688-3695; printed beside, never asserted.
REC = {
    "sar": {
        "P1": ['2.950345e-05', '5.331334e-07', '3.204894e-06', '3.925963e-07', '8.587866e-06'],
        "P2": ['1.617152e-07', '2.678039e-05', '1.100258e-07', '5.252647e-06', '9.167913e-06'],
        "P3": ['3.566888e-06', '4.383337e-07', '2.633419e-05', '4.130965e-07', '1.561175e-06'],
        "P4": ['1.130206e-07', '3.208626e-06', '4.783794e-07', '4.852049e-05', '1.700377e-06'],
    },
    "eabs": {
        "P1": ['3.435314e-01', '4.617936e-02', '1.132236e-01', '3.962809e-02', '1.853415e-01'],
        "P2": ['2.543346e-02', '3.272943e-01', '2.097864e-02', '1.449503e-01', '1.914984e-01'],
        "P3": ['1.194469e-01', '4.187284e-02', '3.245562e-01', '4.064955e-02', '7.902342e-02'],
        "P4": ['2.126223e-02', '1.132895e-01', '4.374377e-02', '4.405474e-01', '8.247125e-02'],
    },
}
REC_RUN_S = 68.2        # :3686
REC_TWIN_S = 34.3       # :5509
REC_CDISP = [6.3865e-02, 2.1232e-01, 8.3482e-02, 5.6184e-02]   # :5510
REC_S_TIE = 9.558e-11   # :3696
REC_IDENT = 2.213e-16   # :3697


def _p0(comm, text):
    if comm.rank == 0:
        print(text, flush=True)


def _run_file(path, outdir, comm):
    """``cli.main(['run', ...])`` with a call-through wrapper; exactly one call asserted."""
    t._fresh_dir(outdir, comm)
    captured = []
    original = study_yaml.run_study_file

    def _call_through(*args, **kwargs):
        out = original(*args, **kwargs)
        captured.append(out)
        return out

    study_yaml.run_study_file = _call_through
    t0 = time.perf_counter()
    try:
        rc = t.cli.main(["run", str(path), "--output-dir", str(outdir)])
    finally:
        study_yaml.run_study_file = original
    wall = time.perf_counter() - t0
    assert study_yaml.run_study_file is original
    rc_all = t._rc_all(rc, comm)
    assert rc_all == (0, 0), rc_all
    assert len(captured) == 1, f"wrapper calls {len(captured)}"
    comm.Barrier()
    return captured[0], wall


def main() -> None:
    comm = MPI.COMM_WORLD
    rank0 = comm.rank == 0
    started = time.perf_counter()

    # ---- (B) ------------------------------------------------------------------------
    assert t.IS_COMPLEX, "complex build required (source dolfinx-complex-mode)"
    assert np.issubdtype(default_scalar_type, np.complexfloating)
    assert comm.size == 2, f"this example is for -n 2; got -n {comm.size}"

    # ---- (K) parse, no mesh -----------------------------------------------------------
    diff = t._spec_diff(t.load_study(t.EM_PROBES), t.load_study(t.FOUR_PORT))
    _, kind, outputs, cooling = t._load(t.EM_PROBES)
    got, want = t._hex(outputs["probes"]), t._hex(t.POINTS)
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        rc_check = t.cli.main(["check", str(t.EM_PROBES)])
    check_all = t._rc_all(rc_check, comm)
    family_lines = [ln for ln in buf.getvalue().splitlines() if ln.startswith("family: ")]
    assert diff == [], diff
    assert kind == "four_port_birdcage" and cooling is None
    assert set(outputs) == {"probes"}
    assert got == want, (got, want)
    assert check_all == (0, 0), check_all
    if rank0:
        assert family_lines == ["family: four_port_birdcage"], family_lines
    _p0(comm, f"[(K)] spec diff {diff}; kind {kind!r}; five points hex-equal {got == want}; "
              f"check rc {check_all}")

    # ---- (X-out) probes outside the phantom are refused before any mesh -------------------
    s = t._mapping(t.EM_PROBES)
    base = {k: v for k, v in s.items() if k != "outputs"}
    rows = [
        ([[0.02, 0.0, 0.0], [0.0, 0.0, 0.05]], "outputs.probes.points[1]"),
        ([[0.031, 0.0, 0.0]], "outputs.probes.points[0]"),
    ]
    fields = [t._field(lambda pts=pts: t._study_from_mapping(
        {**base, "outputs": {"probes": {"points": pts}}})) for pts, _ in rows]
    _p0(comm, "[(X-out)] " + "; ".join(
        f"{pts} -> {f!r}" for (pts, _), f in zip(rows, fields)))
    assert fields == [f for _, f in rows], fields

    # ---- femem run on the EM probe file ----------------------------------------------------
    OUTPUT_DIR.mkdir(exist_ok=True)
    run, wall_run = _run_file(t.EM_PROBES, RUN_DIR, comm)
    summ = run.summary
    pr = summ["probes"]
    order, sar, eab = t._arrays(summ)
    res = run.study.fields_result
    sigma, rho = float(run.study.material.sigma), float(run.study.density_kg_m3)
    pts = np.asarray(t.POINTS)
    _p0(comm, f"[run] femem run {wall_run:.1f} s (gate {REC_RUN_S} s, W1:3686); "
              f"{summ['cells']} cells; ports {order}; sigma {sigma!r}; rho {rho!r}; "
              f"owner {pr['owner']!r}")

    # printed: per-port rows beside REC with digits-equal flags
    n_equal = 0
    n_total = 0
    for i, pid in enumerate(order):
        s_str = [f"{v:.6e}" for v in sar[i]]
        e_str = [f"{v:.6e}" for v in eab[i]]
        eq_s = [a == b for a, b in zip(s_str, REC["sar"][pid])]
        eq_e = [a == b for a, b in zip(e_str, REC["eabs"][pid])]
        n_equal += sum(eq_s) + sum(eq_e)
        n_total += len(eq_s) + len(eq_e)
        _p0(comm, f"  {pid}: SAR [W/kg per 1 V] {s_str}  digits == REC {eq_s}")
        _p0(comm, f"  {pid}: |E| [V/m per 1 V]  {e_str}  digits == REC {eq_e}")
    _p0(comm, f"[REC digits] {n_equal} / {n_total} per-port values equal the W1 strings "
              f"(predicted 40 / 40; printed, not asserted)")
    near = [int(np.argmax(sar[i])) for i in range(len(order))]
    _p0(comm, f"[largest-SAR point per port] {dict(zip(order, near))} "
              f"(predicted P1 0, P2 1, P3 2, P4 3; printed)")

    # ---- (R) ---------------------------------------------------------------------------
    s_mat = np.array([[complex(re, im) for re, im in row] for row in summ["s_matrix"]])
    rel_s = t._rel(s_mat, t._record(1))
    cells = int(summ["cells"])
    _p0(comm, f"[(R)] cells {cells} (record {t.ANS4_EXACT_STRIP_CELL_RECORD}); "
              f"||S - P1||/||P1|| {rel_s:.3e} (gate {REC_S_TIE:.3e}, W1:3696; ASSERTED <= "
              f"{t.LEG_D0_REPRODUCTION_BAND:.0e})")
    assert cells == t.ANS4_EXACT_STRIP_CELL_RECORD
    assert rel_s <= t.LEG_D0_REPRODUCTION_BAND
    assert [t._hex([p])[0] for p in pr["points"]] == t._hex(t.POINTS)
    assert pr["owner"] == "lowest_original_cell"

    # ---- (U') the summary is point_sar on the run's own fields, same five-point batch -------
    own = np.array([
        t.point_sar(res.fields[pid].e_real, res.fields[pid].e_imag, pts, sigma=sigma,
                    rho=rho, comm=res.comm)
        for pid in order
    ])
    exact = bool(np.array_equal(own, sar))
    bitwise = t._bitwise_equal_across_ranks(comm, sar)
    js = json.loads(run.paths["summary"].read_text())
    rest = {k: v for k, v in summ.items() if k != "elapsed_s"}
    js_rest = {k: v for k, v in js.items() if k != "elapsed_s"}
    same_json = t._all(comm, js_rest == rest and (comm.rank != 0 or js == summ))
    ident = float(np.max(np.abs(eab ** 2 * sigma / (2.0 * rho) - sar) / sar))
    _p0(comm, f"[(U')] summary == point_sar on run fields {exact}; bitwise across ranks "
              f"{bitwise}; JSON == memory {same_json}; max |E|^2 sigma/(2 rho) vs SAR rel "
              f"{ident:.3e} (gate {REC_IDENT:.3e}, W1:3697; ASSERTED <= "
              f"{t.EXACT_IDENTITY_RTOL:.0e})")
    assert exact
    assert bitwise
    assert same_json
    assert ident <= t.EXACT_IDENTITY_RTOL
    assert run.probes is None
    assert np.all(sar > 0.0)

    # ---- (J-disk) the probe table read from disk alone ---------------------------------------
    comm.Barrier()
    with open(run.paths["summary"]) as fh:
        disk = json.load(fh)
    dp = disk["probes"]
    assert t._hex(dp["points"]) == t._hex(t.POINTS)
    assert float(dp["sigma_s_per_m"]).hex() == float(run.study.material.sigma).hex()
    assert float(dp["rho_kg_m3"]).hex() == float(run.study.density_kg_m3).hex()
    assert list(dp["values"]) == order
    d_sar = np.array([dp["values"][p]["sar_w_per_kg"] for p in order])
    d_eab = np.array([dp["values"][p]["e_abs_v_per_m"] for p in order])
    assert np.array_equal(d_sar, sar) and np.array_equal(d_eab, eab)
    d_ident = float(np.max(np.abs(d_eab ** 2 * float(dp["sigma_s_per_m"])
                                  / (2.0 * float(dp["rho_kg_m3"])) - d_sar) / d_sar))
    assert d_ident <= t.EXACT_IDENTITY_RTOL
    ok_disk = t._all(comm, True)
    _p0(comm, f"[(J-disk)] points hex-equal, sigma / rho hex-equal, SAR and |E| rows "
              f"array_equal the in-memory rows, on every rank {ok_disk}; disk |E|^2 sigma/"
              f"(2 rho) vs disk SAR rel {d_ident:.3e} (ASSERTED <= {t.EXACT_IDENTITY_RTOL:.0e})")
    comm.Barrier()

    # ---- (c-disp) asserted: point 0 displaced 1 mm radially misses ----------------------------
    disp = np.array([
        t.point_sar(res.fields[pid].e_real, res.fields[pid].e_imag,
                    np.asarray([t.DISPLACED_0]), sigma=sigma, rho=rho, comm=res.comm)[0]
        for pid in order
    ])
    c_disp = np.abs(disp - sar[:, 0]) / sar[:, 0]
    _p0(comm, f"[(c-disp)] 1 mm displaced point 0 per port rel diff "
              f"{[f'{v:.4e}' for v in c_disp]} (gate {REC_CDISP}, W1:5510); min "
              f"{float(c_disp.min()):.4e} = {float(c_disp.min()) / t.L0_BAND:.3e} x L0_BAND "
              f"(ASSERTED > {t.WIDTH_RTOL:.0e})")
    assert float(c_disp.min()) > t.WIDTH_RTOL

    # ---- owner=None equals the owner rule (printed only) ---------------------------------------
    eq_none = []
    for i, pid in enumerate(order):
        leg = t.point_sar(res.fields[pid].e_real, res.fields[pid].e_imag, pts, sigma=sigma,
                          rho=rho, comm=res.comm, owner=None)
        eq_none.append(bool(np.array_equal(leg, sar[i])))
    _p0(comm, f"[owner=None] equals the owner-rule row per port {dict(zip(order, eq_none))} "
              f"(printed; W1:5511-5518)")

    # ---- (I-Py) file = Python: an independent Study, its own field solve -----------------------
    t0 = time.perf_counter()
    twin = t._port_study(degree=1)
    f_twin = float(twin.frequency_hz)  # the twin's own frequency, never typed (WF-17 step 29)
    tsweep = twin.solve_sweep((f_twin,), comm, fields_at=f_twin)  # the public field route (WF-16 step 18)
    tres = tsweep.results[0]
    assert tres.fields is not None and float(tsweep.frequencies_hz[0]) == f_twin
    py =np.array([
        t.point_sar(tres.fields[pid].e_real, tres.fields[pid].e_imag, pts,
                    sigma=float(twin.material.sigma), rho=float(twin.density_kg_m3),
                    comm=tres.comm)
        for pid in order
    ])
    wall_twin = time.perf_counter() - t0
    i_py = float(np.max(np.abs(sar - py) / np.abs(py)))
    _p0(comm, f"[(I-Py)] Python twin {wall_twin:.1f} s (gate {REC_TWIN_S} s, W1:5509), "
              f"{tres.cells} cells; max |SAR_file - SAR_py|/SAR_py {i_py:.3e} (gate 0.000e+00; "
              f"ASSERTED <= {t.L0_BAND:.0e})")
    assert i_py <= t.L0_BAND
    del twin, tres, tsweep
    gc.collect()
    comm.Barrier()

    # ---- setup figure (no-op unless FEM_EM_SETUP_FIGURES=1) -------------------------------------
    msh = res.mesh
    local_tags = np.unique(res.cell_tags.values)
    all_tags = sorted({int(x) for part in comm.allgather(local_tags.tolist()) for x in part})
    air_tags = tuple(x for x in all_tags if x != PHANTOM_CELL_TAG)
    names = {PHANTOM_CELL_TAG: "gel phantom"}
    for x in air_tags:
        names[x] = "air (outer box PEC)"
    FIGURE_DIR.mkdir(exist_ok=True)
    write_setup_figure(
        msh, res.cell_tags, FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names=names, hide_tags=air_tags, translucent_tags=(PHANTOM_CELL_TAG,),
        slice_normal=(0.0, 0.0, 1.0), slice_origin=(0.0, 0.0, 0.0),
        title=("mri:33 -- femem run on the F-small EM file: five point probes in the gel "
               "phantom, per-port point SAR / |E| (WF-1 step 20)"),
        comm=comm)
    gc.collect()
    comm.Barrier()

    # ---- ParaView fields (display only; never compared with the table) --------------------------
    p1_space = fem.functionspace(msh, ("Lagrange", 1))
    dg0 = fem.functionspace(msh, ("DG", 0))
    out_fields = {}
    for pid in order:
        f = res.fields[pid]
        e = f.e_complex
        e_abs = fem.Function(p1_space, name=f"E_abs_{pid}")
        e_abs.interpolate(fem.Expression(ufl.sqrt(ufl.real(ufl.inner(e, e))),
                                         p1_space.element.interpolation_points))
        e_abs.x.array[:] = np.real(e_abs.x.array)
        e_abs.x.scatter_forward()
        out_fields[f"E_abs_{pid}"] = e_abs
        lsar = fem.Function(dg0, name=f"local_sar_{pid}")
        lsar.interpolate(fem.Expression(
            0.5 * ufl.real(f.sigma_field * ufl.inner(e, e)) / rho,
            dg0.element.interpolation_points))
        lsar.x.array[:] = np.real(lsar.x.array)
        lsar.x.scatter_forward()
        out_fields[f"local_sar_{pid}"] = lsar
    tdim = msh.topology.dim
    n_cells = int(msh.topology.index_map(tdim).size_local
                  + msh.topology.index_map(tdim).num_ghosts)
    mids = dmesh.compute_midpoints(msh, tdim, np.arange(n_cells, dtype=np.int32))
    site_vals = np.zeros(n_cells)
    all_sites = [np.asarray(p) for p in t.POINTS] + [np.asarray(t.DISPLACED_0)]
    for k, p in enumerate(all_sites):
        near_k = np.linalg.norm(mids - p[None, :], axis=1) <= SITE_RADIUS_M
        if k == 5:
            near_k &= site_vals == 0.0   # the displaced site never overwrites a real probe
        site_vals[near_k] = float(k + 1)
    sites = fem.Function(dg0, name="probe_sites")
    dofs = dg0.dofmap.list.reshape(-1)
    sites.x.array[:] = 0.0
    sites.x.array[dofs] = site_vals
    sites.x.scatter_forward()
    out_fields["probe_sites"] = sites
    for key, fn in out_fields.items():
        assert fn.name == key, (key, fn.name)
    t0 = time.perf_counter()
    xdmf_path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_combined", msh, res.cell_tags, out_fields, comm=comm)
    adopt_host_ownership(OUTPUT_DIR, comm=comm)
    wall_xdmf = time.perf_counter() - t0
    peak = report_peak_rss(comm, label="EX-138 after outputs")

    # ---- CSV and prints (rank 0 only, no collective inside) ------------------------------------
    if rank0:
        with open(CSV_PATH, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["port", "point", "x", "y", "z", "sar_disk_w_per_kg",
                        "e_abs_disk_v_per_m", "sar_python_twin_w_per_kg",
                        "rel_diff_file_vs_twin", "rec_digits_equal"])
            for i, pid in enumerate(order):
                for j, p in enumerate(t.POINTS):
                    rel = abs(d_sar[i, j] - py[i, j]) / abs(py[i, j])
                    eq = (f"{d_sar[i, j]:.6e}" == REC["sar"][pid][j]
                          and f"{d_eab[i, j]:.6e}" == REC["eabs"][pid][j])
                    w.writerow([pid, j, p[0], p[1], p[2], f"{d_sar[i, j]:.6g}",
                                f"{d_eab[i, j]:.6g}", f"{py[i, j]:.6g}", f"{rel:.3e}", eq])
            for i, pid in enumerate(order):
                w.writerow([pid, "displaced_0", *t.DISPLACED_0, f"{disp[i]:.6g}", "",
                            f"{d_sar[i, 0]:.6g}", f"{c_disp[i]:.3e}", "miss"])
        print("=" * 72)
        print("EX-138 -- femem run on the F-small EM file with point probes (mri:33)")
        print("=" * 72)
        print(f"[wall] femem run {wall_run:.1f}, twin {wall_twin:.1f}, xdmf {wall_xdmf:.1f} s; "
              f"summed peak RSS {peak / 2**30:.2f} GiB")
        print(f"[paraview] wrote {xdmf_path}")
        print(f"[csv] wrote {CSV_PATH}")

    _p0(comm, f"\nAll assertions hold. Total elapsed {time.perf_counter() - started:.1f} s.")


if __name__ == "__main__":
    main()
