"""Example (`EX-139`): ``femem run`` on the four-port F-small sweep file -- one multi-frequency
Touchstone ``.s4p`` and a JSON summary written by the CLI, read back from disk, each frequency's
S tied to the `ANS-4` ``P1`` record, the file tied to an independent Python ``solve_sweep``, and
the 10 <-> 128 MHz exchange that misses.

`WF-1` step 21 (item 572) let ``study.sweep_hz`` ride on a ``four_port_birdcage`` file:
``run_study_file`` meshes once, runs `WF-16`'s ``Study.solve_sweep`` at each listed frequency
and writes ``f_small_sweep.s4p`` plus the JSON summary.  Every earlier ``femem run`` example
writes one single-frequency S (or a thermal result).  This one runs the tracked
``tests/package/wf1/f_small_sweep.yaml`` (degree 1, PEC proxy, 10 / 64 / 128 MHz, exact-rectangle
ports, ``-n 2``) through ``cli.main(["run", ...])`` with a call-through wrapper on
``fem_em_solver.study_yaml.run_study_file`` (restored in ``finally``):

* (B) complex build, ``-n 2``; (K) the file's Study equals ``f_small_four_port.yaml``'s, the
  parsed sweep is hex-equal to the gate's ``FREQUENCIES``, ``femem check`` rc (0, 0);
  (X-sweep) four conflicting files and the Palace backend are refused by field, before any mesh;
* (R') cells and per-frequency S tie to item 192's ``P1`` records on one mesh; (T') the one
  ``.s4p`` ties to memory; (G') reciprocity and passivity per frequency; (U') JSON == memory;
  (I-Py) file = Python ``solve_sweep`` at ``L0_BAND``;
* (D-disk) the two files on disk alone carry the gated per-frequency records;
* (c-swap) the 10 <-> 128 MHz exchange misses (asserted).

Every record, band and helper is imported from the gate module through ``t`` (the ANS-1
rule); record strings sit only in the example-local ``REC`` with their ``log:line`` and are
PRINTED beside the readings, never asserted by equality.  New values print at <= 6 significant
digits; no ``P1`` record and no Touchstone text is printed.

Run it through the example runner (the ``mri:`` group sources the complex build)::

    ./scripts/run_examples.sh -e mri:34 -n 2 -t 300

Not claimed: resonance, tuning or an absolute S; frequency dispersion; degree 2; copper, hole,
capacitor, F-human or sphere sweeps (refused); SAR / B1+ / probes per frequency (refused);
fields per frequency; ``--backend palace`` with a sweep (refused); other widths; any external
check.
"""

from __future__ import annotations

import contextlib
import csv
import gc
import hashlib
import io
import json
import sys
import time
from pathlib import Path

import numpy as np
from mpi4py import MPI

from dolfinx import default_scalar_type, fem

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

import tests.package.test_wf1_yaml_sweep as t  # noqa: E402
from tests.validation.test_port_birdcage_lumped_column import PHANTOM_CELL_TAG  # noqa: E402

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"
BASENAME = "mri_34_femem_run_sweep_file_touchstone"
RUN_DIR = OUTPUT_DIR / "mri_34_run_sweep"
PALACE_DIR = OUTPUT_DIR / "mri_34_palace_refused"
CSV_PATH = OUTPUT_DIR / f"{BASENAME}_s_table.csv"

# Gate-log records (complex, -n 2): W1a = docs/testing/logs/20261020T201454Z_WF-1-s21-W1a-sweep-file.log,
# W1b = ...20261020T201555Z_WF-1-s21-W1b-sweep-python.log; printed beside, never asserted.
REC_RUN_S = 50.3            # W1a:1933
REC_MESH_S = 25.28          # W1a:1933
REC_SWEEP_S = (8.11, 8.02, 8.34)   # W1a:1933
REC_CELLS = 117206          # W1a:1933
REC_R = ("6.601e-11", "9.557e-11", "9.449e-11")    # W1a:1934 (R')
REC_SWAP = (1.2367, 1.2763)                        # W1a:1935 (c-swap)
REC_DS = ("5.541e-14", "6.873e-14", "4.969e-14")   # W1a:1936 (T')
REC_RECIP = ("4.593e-15", "4.509e-16", "6.624e-16")  # W1a:1937-1939 (G')
REC_SIGMA = ("0.999994", "0.999760", "0.999136")     # W1a:1937-1939 (G')
REC_PY_S = 49.7             # W1b:1933
REC_PY_MESH_S = 25.72       # W1b:1933
REC_IPY = ("9.010e-15", "9.110e-16", "0.000e+00")    # W1b:1934 (I-Py)


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


def _db(z):
    return 20.0 * np.log10(abs(z))


def main() -> None:
    comm = MPI.COMM_WORLD
    rank0 = comm.rank == 0
    started = time.perf_counter()

    # ---- (B) ------------------------------------------------------------------------------
    assert t.IS_COMPLEX, "complex build required (source dolfinx-complex-mode)"
    assert np.issubdtype(default_scalar_type, np.complexfloating)
    assert comm.size == t.ANS4_EXACT_STRIP_RECORD_RANK_WIDTH, (
        comm.size, t.ANS4_EXACT_STRIP_RECORD_RANK_WIDTH)
    keys = t.KEYS

    # ---- (K) parse, no mesh (gate :140-163) ---------------------------------------------------
    diff = t._spec_diff(t.load_study(t.SWEEP), t.load_study(t.FOUR_PORT))
    study0, kind, outputs, cooling = t._load(t.SWEEP)
    got = [float(f).hex() for f in outputs["sweep_hz"]]
    want = [float(f).hex() for f in t.FREQUENCIES.values()]
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        rc_check = t.cli.main(["check", str(t.SWEEP)])
    check_all = t._rc_all(rc_check, comm)
    family_lines = [ln for ln in buf.getvalue().splitlines() if ln.startswith("family: ")]
    assert diff == [], diff
    assert kind == "four_port_birdcage"
    assert cooling is None
    assert set(outputs) == {"sweep_hz"}
    assert isinstance(outputs["sweep_hz"], tuple)
    assert got == want, (got, want)
    assert check_all == (0, 0), check_all
    if rank0:
        assert family_lines == ["family: four_port_birdcage"], family_lines
    _p0(comm, f"[(K)] spec diff {diff}; kind {kind!r}; outputs keys {sorted(outputs)}; "
              f"sweep_hz hex-equal {got == want}; frequency_hz {study0.frequency_hz!r}; "
              f"check rc {check_all}")

    # ---- (X-sweep) refused by field, before any mesh ---------------------------------------------
    base = t._mapping(t.SWEEP)
    probes = t._mapping(t.WF1 / "f_small_em_probes.yaml")["outputs"]["probes"]
    rows = [
        (t._with_sweep(base, t.BAD_LISTS[2][0]), "study.sweep_hz"),
        (t._with_sweep(base, t.BAD_LISTS[6][0]), "study.sweep_hz[0]"),
        (t._with_sweep(t._mapping(t.WF1 / "f_small_copper.yaml"), t.GOOD_SWEEP),
         "study.sweep_hz"),
        ({**base, "outputs": {"probes": probes}}, "study.sweep_hz"),
    ]
    fields = [t._field(lambda m=m: t._study_from_mapping(m)) for m, _ in rows]
    labels = ["decreasing list", '"10 MHz" entry', "copper file", "beside outputs.probes"]
    _p0(comm, "[(X-sweep)] " + "; ".join(
        f"{lab}: want {w!r} got {g!r}" for lab, (_, w), g in zip(labels, rows, fields)))
    assert fields == [w for _, w in rows], fields
    OUTPUT_DIR.mkdir(exist_ok=True)
    t._fresh_dir(PALACE_DIR, comm)
    g_pal = t._field(lambda: t._run_study_file_palace(t.SWEEP, output_dir=PALACE_DIR, comm=comm))
    absent = t._all(comm, not PALACE_DIR.exists())
    _p0(comm, f"[(X-sweep)] --backend palace on {t.SWEEP.name}: field {g_pal!r}; nothing written "
              f"on every rank {absent}")
    assert g_pal == "--backend"
    assert absent

    # ---- femem run on the sweep file ---------------------------------------------------------------
    run, wall_run = _run_file(t.SWEEP, RUN_DIR, comm)
    s = run.summary
    sw = run.result
    _p0(comm, f"[run] femem run {wall_run:.1f} s (gate {REC_RUN_S} s, W1a:1933; summary elapsed "
              f"{s['elapsed_s']:.1f} s); mesh {s['mesh_time_s']:.2f} s (gate {REC_MESH_S} s); "
              f"sweeps " + " / ".join(f"{k} MHz {x:.2f} s (gate {g} s)" for k, x, g in
                                      zip(keys, s["sweep_times_s"], REC_SWEEP_S))
              + f"; {s['cells']} cells (gate {REC_CELLS}); frequency_hz {s['frequency_hz']!r} "
                f"beside sweep_hz {s['sweep_hz']}")
    src_sha = hashlib.sha256(t.SWEEP.read_bytes()).hexdigest()
    _p0(comm, f"[sha256] summary source digest == sha256 of the file bytes "
              f"{s['sha256'] == src_sha} (printed)")

    # ---- (R') per frequency, one mesh (gate :244-257) ------------------------------------------------
    rec = t._records(1)
    s_file = np.array([t._unpairs(m) for m in s["s_matrices"]])
    rels = {k: t._rel(s_file[i], rec[k]) for i, k in enumerate(keys)}
    one_mesh = all(r.mesh is sw.results[0].mesh for r in sw.results)
    r_flags = [f"{rels[k]:.3e}" == g for k, g in zip(keys, REC_R)]
    _p0(comm, "[(R')] ||S_f - P1_f||/||P1_f|| " + " / ".join(
        f"{k} MHz {v:.3e}" for k, v in rels.items())
        + f" (gate {' / '.join(REC_R)}, W1a:1934; digits == {r_flags}; ASSERTED <= "
          f"{t.LEG_D0_REPRODUCTION_BAND:.0e}); cells {s['cells']} (record "
          f"{t.ANS4_EXACT_STRIP_CELL_RECORD}); one mesh {one_mesh}")
    assert s["cells"] == t.ANS4_EXACT_STRIP_CELL_RECORD
    assert one_mesh
    assert len(sw.results) == len(keys)
    assert [float(f).hex() for f in s["sweep_hz"]] == want
    for k, v in rels.items():
        assert v <= t.LEG_D0_REPRODUCTION_BAND, k

    # ---- (T') the one Touchstone (gate :269-280) ------------------------------------------------------
    ts = run.paths["touchstone"]
    freqs, mats, ports = t.load_touchstone(ts)
    loaded = np.asarray(mats, dtype=np.complex128)
    d_ts = [float(np.max(np.abs(loaded[i] - sw.s_matrices[i]))) for i in range(len(keys))]
    d_f = [abs(float(freqs[i]) - f) for i, f in enumerate(s["sweep_hz"])]
    _p0(comm, f"[(T')] {ts.name}: {len(freqs)} points; |df| {d_f}; max|dS| "
              + " / ".join(f"{v:.3e}" for v in d_ts)
              + f" (gate {' / '.join(REC_DS)}, W1a:1936; ASSERTED <= {t.ANCHOR_TOL:.0e})")
    assert ts.name == f"{t.SWEEP.stem}.s4p" and s["touchstone"] == ts.name
    assert len(freqs) == len(keys) and loaded.shape == (len(keys), 4, 4)
    assert all(d <= t.FREQ_ATOL_HZ for d in d_f)
    assert all(d <= t.ANCHOR_TOL for d in d_ts)

    # ---- (G') PORT-11's identities per frequency (gate :282-290) --------------------------------------
    for i, k in enumerate(keys):
        S = sw.s_matrices[i]
        recip = float(np.linalg.norm(S - S.T) / np.linalg.norm(S))
        smax = float(np.linalg.svd(S, compute_uv=False).max())
        _p0(comm, f"[(G')] {k} MHz reciprocity {recip:.3e} (gate {REC_RECIP[i]}; ASSERTED <= "
                  f"{t.RECIPROCITY_BAND:.0e}); sigma_max {smax:.6f} (gate {REC_SIGMA[i]}, "
                  f"W1a:1937-1939; ASSERTED <= 1 + {t.PASSIVITY_SIGMA_TOLERANCE:.0e})")
        assert recip <= t.RECIPROCITY_BAND, k
        assert smax <= 1.0 + t.PASSIVITY_SIGMA_TOLERANCE, k

    # ---- (U') JSON = memory; summary = the sweep exactly; no scalar solve (gate :292-310) -------------
    # Timings are rank-local (each rank's own clock); rank 0's are the ones on disk.
    timing = ("elapsed_s", "mesh_time_s", "sweep_times_s")
    js = json.loads(run.paths["summary"].read_text())
    rest = {k: v for k, v in s.items() if k not in timing}
    js_rest = {k: v for k, v in js.items() if k not in timing}
    diff_keys = sorted(k for k in set(js) | set(s) if js.get(k) != s.get(k))
    _p0(comm, f"[(U')] rank 0 keys differing JSON vs memory: {diff_keys} (gate [], W1a:1940)")
    same_json = t._all(comm, js_rest == rest and (comm.rank != 0 or js == s))
    exact = bool(np.array_equal(s_file, sw.s_matrices))
    _p0(comm, f"[(U')] JSON == memory {same_json}; summary S == sweep S exactly {exact}; "
              f"study.result is None {run.study.result is None}")
    assert same_json
    assert exact
    assert run.study.result is None
    agree = all(np.array_equal(a, s_file) for a in comm.allgather(s_file))
    assert agree

    # ---- (D-disk) the two files on disk alone ------------------------------------------------------------
    comm.Barrier()
    with open(run.paths["summary"]) as fh:
        disk = json.load(fh)
    assert [float(f).hex() for f in disk["sweep_hz"]] == want
    assert disk["cells"] == t.ANS4_EXACT_STRIP_CELL_RECORD
    assert disk["touchstone"] == "f_small_sweep.s4p"
    for i in range(len(keys)):
        assert np.array_equal(t._unpairs(disk["s_matrices"][i]), s_file[i])
    d_freqs, d_mats, _ = t.load_touchstone(RUN_DIR / disk["touchstone"])
    s_ts = np.asarray(d_mats, dtype=np.complex128)
    assert s_ts.shape == (len(keys), 4, 4)
    rel_ts = {k: t._rel(s_ts[i], rec[k]) for i, k in enumerate(keys)}
    _p0(comm, "[(D-disk)] JSON from disk: sweep_hz hex-equal, cells == record, touchstone name, "
              "S blocks array_equal memory; .s4p alone vs P1 "
              + " / ".join(f"{k} MHz {v:.3e}" for k, v in rel_ts.items())
              + f" (EX-95 read <= 9.556e-11; ASSERTED <= {t.LEG_D0_REPRODUCTION_BAND:.0e})")
    for k, v in rel_ts.items():
        assert v <= t.LEG_D0_REPRODUCTION_BAND, k
    ok_disk = t._all(comm, True)
    _p0(comm, f"[(D-disk)] on every rank {ok_disk}")
    comm.Barrier()

    # ---- (c-swap) asserted: the 10 <-> 128 MHz exchange misses on the .s4p blocks (gate :259-267) ------
    i10, i128 = keys.index("10"), keys.index("128")
    sw_10 = t._rel(s_ts[i10], rec["128"])
    sw_128 = t._rel(s_ts[i128], rec["10"])
    _p0(comm, f"[(c-swap)] file 10 MHz block vs P1_128 {sw_10:.4f} (hand "
              f"{t.C4_HAND_TABLE[1][('10', '128')]}, gate {REC_SWAP[0]}); file 128 MHz block vs "
              f"P1_10 {sw_128:.4f} (hand {t.C4_HAND_TABLE[1][('128', '10')]}, gate "
              f"{REC_SWAP[1]}, W1a:1935); ASSERTED > {t.LEG_D0_REPRODUCTION_BAND:.0e} "
              f"({sw_10 / t.LEG_D0_REPRODUCTION_BAND:.2e} x / "
              f"{sw_128 / t.LEG_D0_REPRODUCTION_BAND:.2e} x the bar)")
    assert sw_10 > t.LEG_D0_REPRODUCTION_BAND
    assert sw_128 > t.LEG_D0_REPRODUCTION_BAND
    assert sw_10 > 1.0e-6 and sw_128 > 1.0e-6, "blind instrument"

    # ---- (I-Py) file = Python solve_sweep (gate :319-337, substitution 4) ------------------------------
    t0 = time.perf_counter()
    twin, _phantom = t._port_study(1)
    tsw = twin.solve_sweep(tuple(t.FREQUENCIES.values()), comm)
    wall_twin = time.perf_counter() - t0
    _p0(comm, f"[(I-Py)] Python solve_sweep {wall_twin:.1f} s (gate {REC_PY_S} s, W1b:1933); "
              f"mesh {tsw.mesh_time:.2f} s (gate {REC_PY_MESH_S} s); sweeps "
              + " / ".join(f"{k} MHz {x:.2f} s" for k, x in zip(keys, tsw.sweep_times))
              + f"; {tsw.results[0].cells} cells (file {s['cells']})")
    assert np.array_equal(np.asarray(tsw.frequencies_hz), np.asarray(s["sweep_hz"]))
    assert int(tsw.results[0].cells) == int(s["cells"])
    py_s = np.asarray(tsw.s_matrices)
    rels_py = [t._rel(py_s[i], s_file[i]) for i in range(len(keys))]
    _p0(comm, "[(I-Py)] ||S_py - S_file||/||S_file|| " + " / ".join(
        f"{k} MHz {v:.3e}" for k, v in zip(keys, rels_py))
        + f" (gate {' / '.join(REC_IPY)}, W1b:1934; predicted <= 1e-12; ASSERTED <= "
          f"{t.L0_BAND:.0e})")
    for k, v in zip(keys, rels_py):
        assert v <= t.L0_BAND, k
    rel_py_ts = [t._rel(py_s[i], s_ts[i]) for i in range(len(keys))]
    del twin, tsw
    gc.collect()
    comm.Barrier()

    # ---- setup figure (no-op unless FEM_EM_SETUP_FIGURES=1) -------------------------------------------------
    res0 = sw.results[0]
    msh, cell_tags, facet_tags = res0.mesh, res0.cell_tags, res0.facet_tags
    local_tags = np.unique(cell_tags.values)
    all_tags = sorted({int(x) for part in comm.allgather(local_tags.tolist()) for x in part})
    air_tags = tuple(x for x in all_tags if x != PHANTOM_CELL_TAG)
    names = {PHANTOM_CELL_TAG: "gel phantom"}
    for x in air_tags:
        names[x] = "air (outer box PEC)"
    FIGURE_DIR.mkdir(exist_ok=True)
    write_setup_figure(
        msh, cell_tags, FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names=names, hide_tags=air_tags, translucent_tags=(PHANTOM_CELL_TAG,),
        slice_normal=(0.0, 0.0, 1.0), slice_origin=(0.0, 0.0, 0.0),
        title=("mri:34 -- femem run on the F-small sweep file: 10 / 64 / 128 MHz on one mesh, "
               "one .s4p (WF-1 step 21)"),
        comm=comm)
    gc.collect()
    comm.Barrier()

    # ---- ParaView: cell tags, facet tags (port sheets), DG0 sigma (display only) ----------------------------
    geometry = run.study.geometry
    dg0 = fem.functionspace(msh, ("DG", 0))
    sigma_out = fem.Function(dg0, name="sigma_s_per_m")
    sigma_out.x.array[:] = 0.0
    sigma_by_tag = {
        geometry.conductor_tag: float(geometry.conductor_sigma_s_per_m),
        geometry.phantom_tag: float(np.real(run.study.material.sigma)),
    }
    idx = np.asarray(cell_tags.indices)
    val = np.asarray(cell_tags.values)
    for c, tag in zip(idx, val):
        if int(tag) in sigma_by_tag:
            sigma_out.x.array[dg0.dofmap.cell_dofs(int(c))[0]] = sigma_by_tag[int(tag)]
    sigma_out.x.array[:] = np.real(sigma_out.x.array)
    sigma_out.x.scatter_forward()
    out_fields = {"sigma_s_per_m": sigma_out}
    for key, fn in out_fields.items():
        assert fn.name == key, (key, fn.name)
    t0 = time.perf_counter()
    xdmf_path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_combined", msh, cell_tags, out_fields, comm=comm,
        facet_tags=facet_tags)
    adopt_host_ownership(OUTPUT_DIR, comm=comm)
    wall_xdmf = time.perf_counter() - t0
    peak = report_peak_rss(comm, label="EX-139 after outputs")

    # ---- CSV and prints (rank 0 only, no collective inside) ------------------------------------------------
    if rank0:
        with open(CSV_PATH, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["f_hz", "i", "j", "re", "im", "abs_db", "block_rel_vs_P1",
                        "block_rel_vs_python_twin"])
            for k_i, k in enumerate(keys):
                for i in range(4):
                    for j in range(4):
                        z = s_ts[k_i][i, j]
                        w.writerow([f"{float(d_freqs[k_i]):.6g}", i + 1, j + 1,
                                    f"{z.real:.6g}", f"{z.imag:.6g}", f"{_db(z):.4f}",
                                    f"{rel_ts[k]:.3e}", f"{rel_py_ts[k_i]:.3e}"])
            w.writerow(["swap_10_block_vs_P1_128", "", "", "", "", "", f"{sw_10:.4f}", ""])
            w.writerow(["swap_128_block_vs_P1_10", "", "", "", "", "", f"{sw_128:.4f}", ""])
        print("=" * 72)
        print("EX-139 -- femem run on the F-small sweep file (mri:34)")
        print("=" * 72)
        for k_i, k in enumerate(keys):
            print(f"[S] {k} MHz from the .s4p: |S11| {_db(s_ts[k_i][0, 0]):.4f} dB, "
                  f"|S21| {_db(s_ts[k_i][1, 0]):.4f} dB (no resonance or tuning reading)")
        print(f"[wall] femem run {wall_run:.1f}, python twin {wall_twin:.1f}, xdmf "
              f"{wall_xdmf:.1f} s; summed peak RSS {peak / 2**30:.2f} GiB")
        print(f"[paraview] wrote {xdmf_path}")
        print(f"[csv] wrote {CSV_PATH}")

    _p0(comm, f"\nAll assertions hold. Total elapsed {time.perf_counter() - started:.1f} s.")


if __name__ == "__main__":
    main()
