"""`WF-8` step 1 (§9 item 328, chain U22): figures that equal their records.

Every plotted array, read back from the figure object, is bitwise equal
(``np.array_equal``, float64) to the source record it draws; a planted
(1 + 1e-6) perturbation is detected at exactly one index on every figure; a
deliberately wrong source is detected; files exist and are non-empty and only
rank 0 writes; ``--check-links`` over ``docs/guide`` reads 0 dead.  All exact.

Sources are gitignored build products on this box (review ruling (1)); each
sha256 is printed.  Source arrays are parsed here independently of the
plotter (skrf for S, h5py + vertex mean for the slice, the ``csv`` module +
``float`` for the probe series).  Run at ``-n 2``.
"""

from __future__ import annotations

import csv
import hashlib
import os
from pathlib import Path

import numpy as np
import pytest
from mpi4py import MPI

from fem_em_solver.post import figures
from tests.unit.test_doc_reference_exit_codes import run_checker

REPO = Path(__file__).resolve().parents[2]
S4P = REPO / "examples/mri/paraview_output/mri_11_study_frequency_sweep.s4p"
XDMF = REPO / "examples/mri/paraview_output/mri_10_study_temperature_combined.xdmf"
H5 = XDMF.with_suffix(".h5")
CSV = REPO / "examples/materials/paraview_output/materials_13_virtual_probe_in_gel_ball_series.csv"
FIELD = "real_SAR"
COLUMN = "probe_c0_K"
WRONG_COLUMN = "probe_c1_K"

comm = MPI.COMM_WORLD


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _all(flag: bool) -> bool:
    """Rank 0's verdict, broadcast so every rank asserts the same thing."""
    return bool(comm.bcast(bool(flag), root=0))


def _exact(plotted, source) -> bool:
    p, s = np.asarray(plotted), np.asarray(source)
    return p.dtype == np.float64 and s.dtype == np.float64 and np.array_equal(p, s)


@pytest.fixture(scope="module")
def outdir(tmp_path_factory):
    local = str(tmp_path_factory.mktemp("wf8")) if comm.rank == 0 else None
    path = Path(comm.bcast(local, root=0))
    return path


@pytest.fixture(scope="module", autouse=True)
def sources():
    present = all(p.is_file() for p in (S4P, XDMF, H5, CSV))
    assert _all(present), "a WF-8 source record is absent: regenerate it with its example (ruling (1))"
    if comm.rank == 0:
        for p in (S4P, XDMF, H5, CSV):
            print(f"[wf8] source {p.relative_to(REPO)} sha256={_sha(p)} bytes={p.stat().st_size}", flush=True)


def _source_s(path=S4P):
    import skrf

    nw = skrf.Network(str(path))
    return nw.f, nw.s_mag, nw.s_deg


def _source_slice(h5=H5, field=FIELD):
    import h5py

    with h5py.File(h5, "r") as f:
        geo = f["/Mesh/mesh/geometry"][...].astype(np.float64)
        top = f["/Mesh/mesh/topology"][...].astype(np.int64)
        val = f[f"/Function/{field}/0"][...].ravel()
    return geo[top].mean(axis=1), val


def _source_csv(column, path=CSV):
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))
    return (
        np.array([float(r["t_s"]) for r in rows], dtype=np.float64),
        np.array([float(r[column]) for r in rows], dtype=np.float64),
    )


def _line_xy(line):
    return np.asarray(line.get_xdata(orig=True)), np.asarray(line.get_ydata(orig=True))


def _coll_values(coll):
    arr = coll.get_array()
    assert not np.ma.is_masked(arr)
    return np.asarray(np.ma.getdata(arr))


def _coll_xy(coll):
    return np.asarray(np.ma.getdata(coll.get_offsets()))


# ------------------------------------------------------------------------- (S)


def test_sparameters_bitwise(outdir):
    out = outdir / "wf8_sparams.png"
    res = figures.plot_sparameters(S4P, out, comm=comm)
    ok = True
    if comm.rank == 0:
        f, mag, deg = _source_s()
        n = mag.shape[1]
        print(f"[wf8] (S) frequencies={f.tolist()} ports={n}", flush=True)
        assert n == 4 and f.shape[0] == 3
        count = 0
        for i in range(n):
            for j in range(n):
                xm, ym = _line_xy(res.lines[(i, j, "mag")])
                xp, yp = _line_xy(res.lines[(i, j, "deg")])
                good = _exact(xm, f) and _exact(ym, mag[:, i, j]) and _exact(xp, f) and _exact(yp, deg[:, i, j])
                ok = ok and good
                count += int(good)
        print(f"[wf8] (S) bitwise-equal (i,j) pairs: {count}/16", flush=True)
    else:
        assert res is None
    assert _all(ok)


def test_sparameters_controls():
    ok = True
    if comm.rank == 0:
        f, mag, deg = _source_s()
        # (c1) one planted perturbation -> exactly one mismatching index
        mag_p = mag.copy()
        mag_p[1, 0, 0] *= 1.0 + 1e-6
        res = figures.draw_sparameters(f, mag_p, deg)
        _, y = _line_xy(res.lines[(0, 0, "mag")])
        bad = figures.mismatch_indices(y, mag[:, 0, 0])
        print(f"[wf8] (S c1) mismatch indices: {bad.tolist()}", flush=True)
        ok = ok and bad.tolist() == [1]
        # (c2) wrong source under the S11 label
        res = figures.draw_sparameters(f, mag, deg)
        _, y = _line_xy(res.lines[(0, 1, "mag")])
        bad2 = figures.mismatch_indices(y, mag[:, 0, 0])
        print(f"[wf8] (S c2) S12 drawn vs S11 source: {bad2.size} mismatching", flush=True)
        ok = ok and bad2.size > 0 and not np.array_equal(y, mag[:, 0, 0])
    assert _all(ok)


# ------------------------------------------------------------------------- (X)


def test_slice_bitwise(outdir):
    out = outdir / "wf8_slice.png"
    res = figures.plot_xdmf_slice(XDMF, FIELD, out, axis="z", comm=comm)
    ok = True
    if comm.rank == 0:
        cent, val = _source_slice()
        idx = res.indices
        print(f"[wf8] (X) cells={val.shape[0]} selected len(idx)={len(idx)} first={idx[:5].tolist()}", flush=True)
        ok = len(idx) > 0
        ok = ok and _exact(_coll_values(res.collection), val[idx])
        ok = ok and _exact(_coll_xy(res.collection), cent[idx][:, :2])
    else:
        assert res is None
    assert _all(ok)


def test_slice_controls():
    ok = True
    if comm.rank == 0:
        cent, val = _source_slice()
        idx = figures.select_slab(cent, axis="z")
        k = int(np.flatnonzero(val[idx] != 0.0)[0])  # a non-zero value, so the plant is not a no-op
        val_p = val.copy()
        val_p[idx[k]] *= 1.0 + 1e-6
        res = figures.draw_slice(cent, val_p, idx)
        bad = figures.mismatch_indices(_coll_values(res.collection), val[idx])
        print(f"[wf8] (X c1) planted at {k}; mismatch indices: {bad.tolist()}", flush=True)
        ok = bad.tolist() == [k]
        # (c2) a wrong field (real_q) drawn against the real_SAR source
        import h5py

        with h5py.File(H5, "r") as f:
            other = f["/Function/real_q/0"][...].ravel()
        res = figures.draw_slice(cent, other, idx)
        bad2 = figures.mismatch_indices(_coll_values(res.collection), val[idx])
        print(f"[wf8] (X c2) real_q drawn vs real_SAR source: {bad2.size} mismatching", flush=True)
        ok = ok and bad2.size > 0
    assert _all(ok)


# ------------------------------------------------------------------------- (P)


def test_probe_bitwise(outdir):
    out = outdir / "wf8_probe.png"
    res = figures.plot_probe_series(CSV, COLUMN, out, comm=comm)
    ok = True
    if comm.rank == 0:
        t, y = _source_csv(COLUMN)
        x_pl, y_pl = _line_xy(res.lines[COLUMN])
        print(f"[wf8] (P) rows={t.shape[0]}", flush=True)
        ok = _exact(x_pl, t) and _exact(y_pl, y)
    else:
        assert res is None
    assert _all(ok)


def test_probe_controls():
    ok = True
    if comm.rank == 0:
        t, y = _source_csv(COLUMN)
        y_p = y.copy()
        y_p[5] *= 1.0 + 1e-6
        res = figures.draw_series(t, y_p, label=COLUMN)
        _, y_pl = _line_xy(res.lines[COLUMN])
        bad = figures.mismatch_indices(y_pl, y)
        print(f"[wf8] (P c1) mismatch indices: {bad.tolist()}", flush=True)
        ok = bad.tolist() == [5]
        t_w, y_w = _source_csv(WRONG_COLUMN)
        res = figures.draw_series(t_w, y_w, label=COLUMN)
        _, y_pl = _line_xy(res.lines[COLUMN])
        bad2 = figures.mismatch_indices(y_pl, y)
        print(f"[wf8] (P c2) {WRONG_COLUMN} drawn vs {COLUMN}: {bad2.size} mismatching", flush=True)
        ok = ok and bad2.size > 0
    assert _all(ok)


# ------------------------------------------------------------------------- (F)


def test_files_written_by_rank0_only(outdir):
    comm.Barrier()
    ok = True
    if comm.rank == 0:
        names = sorted(os.listdir(outdir))
        sizes = {n: os.path.getsize(outdir / n) for n in names}
        print(f"[wf8] (F) files={sizes}", flush=True)
        ok = names == ["wf8_probe.png", "wf8_slice.png", "wf8_sparams.png"] and all(s > 0 for s in sizes.values())
    assert _all(ok)


def test_femem_figure_cli(tmp_path_factory):
    """`femem figure` is thin over the three: each kind returns 0 on every rank
    and rank 0 writes exactly the three non-empty files."""
    from fem_em_solver.cli import main

    local = str(tmp_path_factory.mktemp("wf8cli")) if comm.rank == 0 else None
    d = Path(comm.bcast(local, root=0))
    cases = [("sparams", S4P, []), ("slice", XDMF, ["--name", FIELD]), ("probe", CSV, ["--column", COLUMN])]
    rcs = [main(["figure", kind, str(src), str(d / f"cli_{kind}.png"), *extra]) for kind, src, extra in cases]
    assert rcs == [0, 0, 0]
    ok = True
    if comm.rank == 0:
        sizes = {kind: os.path.getsize(d / f"cli_{kind}.png") for kind, _, _ in cases}
        print(f"[wf8] (CLI) rc={rcs} sizes={sizes}", flush=True)
        ok = sorted(os.listdir(d)) == ["cli_probe.png", "cli_slice.png", "cli_sparams.png"] and all(
            s > 0 for s in sizes.values()
        )
    assert _all(ok)


# ------------------------------------------------------------------------- (D)


def test_guide_links():
    dead = None
    if comm.rank == 0:
        _rc, counts = run_checker("--docs-root", "docs/guide", "--check-links")
        dead = counts["dead"]
        print(f"[wf8] (D) docs/guide --check-links dead={dead} counts={counts}", flush=True)
    dead = comm.bcast(dead, root=0)
    assert dead == 0
