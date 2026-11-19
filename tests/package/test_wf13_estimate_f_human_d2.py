"""`WF-13` step 11 (§9 item 585, chain U77): the F-human degree-2 price row FH2.

Row FH2 is `WF-10`'s degree-2 copper F-human hole solve (`xl` entry 41): the
node pair (the EM solve, build excluded) and the summed ``ru_maxrss``.  A
degree-2 F-human Study (`WF-16`'s refusal lifted for ``estimate`` only) is
priced from the ten-row degree-2 fit that includes FH2.

Brief: docs/planning/chunks/WF-13.md, "§9 item 585" and the review ruling
2026-10-22 03:00 beneath it.  Run under ``mpiexec -n 2``.

Anchors (asserted; bands imported by name, none set here):
  (FH)      FH2's seven fields equal the transcribed literals; unknowns ==
            N_DOFS_D2_RECORD, cells == HOLE_CELL_RECORD; exactly one FH2 row.
  (FH-log)  the row is entry 41's log: one ``price:`` line whose node pair,
            RSS summed and ``-n`` equal the row; every ``global unknowns`` line
            equals the row's N; ``- Status: 0``.
  (LOO)     ``_predict`` at FH2's N and ``-n 16`` on the table without FH2: the
            seconds ratio inside ``[1/BRACKET, BRACKET]``; rows used are the
            nine other degree-2 rows.
  (U7-d2)   (knob ``FEM_EM_WF13_S11=1``, one real build, memoized) a
            ``tmp_path`` copy of the tune file at ``degree: 2``: unknowns ==
            N_DOFS_D2_RECORD, cells == HOLE_CELL_RECORD, degree 2, "FH2" used,
            not_priced == ("thermal transient",), no solve; (CLI = API) rc
            [0, 0], printed fields equal the API's, no output directory.

Negative controls (asserted, by construction): (c-N) an N-blind estimator
misses FH2's seconds by more than ×2; (c-deg) the degree-1 fit at FH2's N and
width misses by more than ×2; (c-row) FH2's unknowns ×1.1 moves the in-sample
price by more than ``EXACT_IDENTITY_RTOL``; (c-plant) FH2 carries the largest
hat weight in P6's LOO fit on the live table, and ``PLANT_K`` on FH2 pushes
P6's LOO memory and time outside ×2.

Printed only: the LOO memory ratio (and against the cgroup memory.peak), the
in-sample ratios, the time-definition sensitivity, the LOO without P6, the LOO
against every other degree-2 F-human hole window on record (duplicates, never
fitted), the ``-n 2`` price of the degree-2 file and the width-mismatch reading.
"""

from __future__ import annotations

import dataclasses
import os
import re
import time
from pathlib import Path

import numpy as np
import pytest
from mpi4py import MPI

from fem_em_solver.study import PRICE_TABLE, FHumanBirdcageGeometry, _predict
from tests.log_records import read_log_text
from tests.package.test_wf1_study_yaml import WF1
from tests.package.test_wf13_cli_estimate import _diff, _parse, _patch_solves

KNOB = "FEM_EM_WF13_S11"
REPO = Path(__file__).resolve().parents[2]
ENTRY41_LOG = REPO / "docs/testing/logs/20261022T070008Z_WF-10-step4-H28.log"
TUNE_FILE = WF1 / "f_human_hole_copper_tune.yaml"
# FH2 transcribed: 20261022T070008Z_WF-10-step4-H28.log:6855 (node pair, RSS
# summed), N :6613 (2E + 2F); cells = the hole record.
FH2_LITERALS = dict(
    row="FH2", degree=2, ranks=16, cells=369387, unknowns=2570024, seconds=298.9, peak_gib=88.58
)
# Entry 41's cgroup memory.peak, :14681 (printed only; the P5 precedent).
ENTRY41_MEMORY_PEAK_GIB = 87.12
# Entry 41's other time definitions (printed sensitivity only): the factorising
# solve pair :6827 (293.22 s), [solve] done :6788 (281.1 s), the second
# module's node pair :13598 (290.0 s).
SENSITIVITY_SECONDS = (
    ("solve pair :6827", 293.22),
    ("[solve] done :6788", 281.1),
    ("second-module node pair :13598", 290.0),
)
# Other degree-2 F-human hole windows on record, -n 16 — duplicates of FH2's
# fixture, printed, never fitted (item 212's precedent).
DUPLICATE_WINDOWS = (
    ("entry 35 20261016T070008Z_WF-10-step4-H20.log:6854", 367.1, 88.72),
    ("entry 35 second module :13465", 333.4, None),
    ("entry 42 20261019T070008Z_MAT-15-step6-P22-d2root.log:6886", 300.5, 88.40),
    ("H26 xxl 20261017T070008Z_WF-10-step1e-H26.log:8641 (low)", 375.2, 104.90),
    ("H26 xxl (high)", 436.5, None),
)
PRICE_RE = re.compile(
    r"\[WF-10 step2\] price: build ([\d.]+) s, .*node pair ([\d.]+) s, .*RSS summed "
    r"([\d.]+) GiB; window elapsed ([\d.]+) s at -n (\d+)"
)
UNKNOWNS_RE = re.compile(r"global unknowns (\d+) vs 2E\+2F")
NINE_ROWS = ("P2", "F16", "Cu7", "Cu7h", "A6", "S4c", "P4", "P5", "P6")


def _say(msg):
    if MPI.COMM_WORLD.rank == 0:
        print(msg, flush=True)


def _row(name, table=PRICE_TABLE):
    rows = [r for r in table if r["row"] == name]
    assert len(rows) == 1, name
    return rows[0]


def _no_fh2():
    return tuple(r for r in PRICE_TABLE if r["row"] != "FH2")


def _in(r, bracket):
    return 1.0 / bracket <= r <= bracket


def test_wf13_s11_fh2_row_is_its_log():
    """(FH), (FH-log).  No mesh."""
    from tests.validation.test_wf10_f_human_tuned_b1_sar import N_DOFS_D2_RECORD
    from tests.validation.test_wf12_f_human_copper_matrix import HOLE_CELL_RECORD

    fh2 = _row("FH2")
    # (FH)
    assert {k: fh2[k] for k in FH2_LITERALS} == FH2_LITERALS
    assert set(fh2) == set(FH2_LITERALS)
    assert fh2["unknowns"] == N_DOFS_D2_RECORD
    assert fh2["cells"] == HOLE_CELL_RECORD

    # (FH-log)
    text = read_log_text(ENTRY41_LOG)
    prices = PRICE_RE.findall(text)
    unknowns = [int(u) for u in UNKNOWNS_RE.findall(text)]
    status_ok = any(line.strip() == "- Status: 0" for line in text.splitlines())
    _say(
        f"\n[WF-13 s11 FH-log] {ENTRY41_LOG.name}: price lines {prices}; "
        f"global unknowns {unknowns}; Status 0 footer {status_ok}"
    )
    assert len(prices) == 1
    build_s, node_s, rss_gib, window_s, ranks = prices[0]
    assert float(node_s) == fh2["seconds"]
    assert float(rss_gib) == fh2["peak_gib"]
    assert int(ranks) == fh2["ranks"]
    assert len(unknowns) >= 1
    assert all(u == fh2["unknowns"] for u in unknowns)
    assert status_ok
    _say(
        f"[WF-13 s11 FH-log] not the row (printed): build {build_s} s, module window "
        f"{window_s} s; cgroup memory.peak {ENTRY41_MEMORY_PEAK_GIB} GiB"
    )


def test_wf13_s11_fh2_loo_and_controls():
    """(LOO), (c-N), (c-deg), (c-row), (c-plant); printed readings.  No mesh."""
    from tests.package.test_wf13_study_estimate import BRACKET, PLANT_K, _hat_weights, _ratios
    from tests.validation.test_birdcage_sar_mass_averaged import EXACT_IDENTITY_RTOL

    fh2 = _row("FH2")
    n, w = fh2["unknowns"], fh2["ranks"]
    nine = _no_fh2()

    # (LOO) the weekly's target.
    m_loo, s_loo, used = _predict(n, 2, w, table=nine)
    rt = s_loo / fh2["seconds"]
    rm = m_loo / fh2["peak_gib"]
    _say(
        f"\n[WF-13 s11 LOO] FH2 at N={n}, -n {w}: predicted {s_loo:.3f} s vs "
        f"{fh2['seconds']} s -> ratio {rt:.4f} (asserted in [{1/BRACKET}, {BRACKET}]); "
        f"rows {used}"
    )
    _say(
        f"[WF-13 s11 LOO] memory (printed only): {m_loo:.3f} GiB vs {fh2['peak_gib']} GiB "
        f"-> ratio {rm:.4f} ({'inside' if _in(rm, BRACKET) else 'OUTSIDE'} x2); vs cgroup "
        f"memory.peak {ENTRY41_MEMORY_PEAK_GIB} GiB -> {m_loo / ENTRY41_MEMORY_PEAK_GIB:.4f}"
    )
    hw = _hat_weights(PRICE_TABLE, "FH2")
    _say("[WF-13 s11 LOO] hat weights at FH2: " + ", ".join(f"{k} {v:.3f}" for k, v in hw.items()))
    assert used == NINE_ROWS
    assert _in(rt, BRACKET)

    # Printed: in-sample, sensitivity, LOO without P6, duplicates, -n 2 price.
    m_in, s_in, used_in = _predict(n, 2, w)
    _say(
        f"[WF-13 s11 in-sample] rows {used_in}: FH2 time {s_in / fh2['seconds']:.4f} "
        f"memory {m_in / fh2['peak_gib']:.4f}"
    )
    for label, secs in SENSITIVITY_SECONDS:
        _say(f"[WF-13 s11 sensitivity] LOO vs {label} {secs} s -> {s_loo / secs:.4f}")
    _, s_np6, used_np6 = _predict(n, 2, w, table=tuple(r for r in nine if r["row"] != "P6"))
    _say(
        f"[WF-13 s11 LOO without P6] {s_np6:.2f} s -> ratio {s_np6 / fh2['seconds']:.4f} "
        f"(rows {used_np6})"
    )
    for label, secs, gib in DUPLICATE_WINDOWS:
        mem = f"; memory {m_loo / gib:.4f} vs {gib} GiB" if gib else ""
        _say(
            f"[WF-13 s11 duplicate, not fitted] {label}: {secs} s -> LOO ratio "
            f"{s_loo / secs:.4f} ({'inside' if _in(s_loo / secs, BRACKET) else 'OUTSIDE'} x2){mem}"
        )
    m_n2, s_n2, _ = _predict(n, 2, 2)
    _say(
        f"[WF-13 s11 -n 2 price] degree-2 F-human file, ten-row fit: {s_n2:.1f} s / "
        f"{m_n2:.2f} GiB (no -n 2 degree-2 F-human solve exists)"
    )
    _, s_w2, _ = _predict(n, 2, 2, table=nine)
    _say(
        f"[WF-13 s11 width mismatch, printed] LOO at -n 2: {s_w2:.1f} s -> ratio "
        f"{s_w2 / fh2['seconds']:.2f}"
    )

    # (c-N) N-blind estimator over the nine rows.
    gm = float(
        np.exp(np.mean(np.log([r["seconds"] * r["ranks"] for r in nine if r["degree"] == 2])))
    )
    gmm = float(np.exp(np.mean(np.log([r["peak_gib"] for r in nine if r["degree"] == 2]))))
    rn = (gm / w) / fh2["seconds"]
    _say(
        f"[WF-13 s11 c-N] N-blind {gm / w:.2f} s vs {fh2['seconds']} s -> ratio {rn:.4f}; "
        f"memory twin {gmm / fh2['peak_gib']:.4f} (printed)"
    )
    assert not _in(rn, BRACKET)

    # (c-deg) degree confusion: the degree-1 fit at FH2's N and width.
    m_d1, s_d1, used_d1 = _predict(n, 1, w)
    rd = s_d1 / fh2["seconds"]
    _say(
        f"[WF-13 s11 c-deg] degree-1 fit {used_d1} at FH2 N/-n {w}: {s_d1:.2f} s -> ratio "
        f"{rd:.4f}; memory {m_d1 / fh2['peak_gib']:.4f} (printed)"
    )
    assert not _in(rd, BRACKET)

    # (c-row) FH2 is load-bearing in the fit.
    bumped = tuple(
        dict(r, unknowns=int(r["unknowns"] * 1.1)) if r["row"] == "FH2" else r for r in PRICE_TABLE
    )
    _, s_b, _ = _predict(n, 2, w, table=bumped)
    rel = abs(s_b / s_in - 1.0)
    _say(f"[WF-13 s11 c-row] FH2 unknowns x1.1: {s_b!r} vs {s_in!r} (rel {rel:.3e})")
    assert rel > EXACT_IDENTITY_RTOL

    # (c-plant) on the live ten-row table.
    weights = _hat_weights(PRICE_TABLE, "P6")
    top = max(weights, key=weights.get)
    wt = weights["FH2"]
    rm0, rt0 = _ratios(PRICE_TABLE, "P6")
    kmin_m = (BRACKET / rm0) ** (1.0 / wt) if rm0 < BRACKET else 1.0
    kmin_t = (BRACKET / rt0) ** (1.0 / wt) if rt0 < BRACKET else 1.0
    planted = tuple(
        (
            dict(r, peak_gib=r["peak_gib"] * PLANT_K, seconds=r["seconds"] * PLANT_K)
            if r["row"] == "FH2"
            else r
        )
        for r in PRICE_TABLE
    )
    rmp, rtp = _ratios(planted, "P6")
    _say(
        "[WF-13 s11 c-plant] P6 LOO hat weights: "
        + ", ".join(f"{k} {v:.3f}" for k, v in weights.items())
        + f"; largest {top}; unplanted mem {rm0:.3f} time {rt0:.3f}; min failing k mem "
        f"{kmin_m:.2f} time {kmin_t:.2f}; FH2 x{PLANT_K:g} -> mem {rmp:.3f} time {rtp:.3f}"
    )
    assert top == "FH2"
    assert PLANT_K >= max(kmin_m, kmin_t)
    assert not _in(rmp, BRACKET)
    assert not _in(rtp, BRACKET)


def test_wf13_s11_f_human_degree2_priced(monkeypatch, capsys, tmp_path):
    """(U7-d2), (CLI = API) on the degree-2 tune copy; one real F-human build."""
    if os.environ.get(KNOB) != "1":
        pytest.skip(f"{KNOB}=1 runs the one-build window")
    from fem_em_solver.cli import main
    from fem_em_solver.study import ResourceEstimate
    from fem_em_solver.study_yaml import load_study
    from tests.validation.test_wf10_f_human_tuned_b1_sar import N_DOFS_D2_RECORD
    from tests.validation.test_wf12_f_human_copper_matrix import HOLE_CELL_RECORD

    comm = MPI.COMM_WORLD
    assert comm.size == 2
    _patch_solves(monkeypatch)

    text = TUNE_FILE.read_text()
    lines = text.splitlines(keepends=True)
    hits_d1 = [i for i, ln in enumerate(lines) if ln.rstrip("\n") == "  degree: 1"]
    assert len(hits_d1) == 1
    lines[hits_d1[0]] = "  degree: 2\n"
    d2_file = tmp_path / "f_human_hole_copper_tune_d2.yaml"
    d2_file.write_text("".join(lines))

    orig = FHumanBirdcageGeometry.build
    memo, count = {}, {"builds": 0, "hits": 0}

    def _memo_build(self, c=MPI.COMM_WORLD):
        key = (repr(self), c.size)
        if key in memo:
            count["hits"] += 1
            return memo[key]
        count["builds"] += 1
        t0 = time.perf_counter()
        memo[key] = orig(self, c)
        _say(f"[WF-13 s11 U7-d2] real F-human build {time.perf_counter() - t0:.1f} s")
        return memo[key]

    monkeypatch.setattr(FHumanBirdcageGeometry, "build", _memo_build)
    names = {f.name for f in dataclasses.fields(ResourceEstimate)}
    t_start = time.perf_counter()
    study = load_study(d2_file)
    est = study.estimate(comm)
    out_dir = tmp_path / "out-d2"
    capsys.readouterr()
    rc = main(["run", "--estimate", str(d2_file), "--output-dir", str(out_dir)])
    out = capsys.readouterr().out
    rcs = comm.allgather(rc)
    diff = _diff(_parse(out, names), est) if comm.rank == 0 else set()
    diff = comm.bcast(diff, root=0)
    out_exists = comm.allgather(out_dir.exists())
    builds = comm.allgather(count["builds"])
    hits = comm.allgather(count["hits"])
    with capsys.disabled():
        _say(
            f"\n[WF-13 s11 U7-d2] {d2_file.name}: cells {est.cells} unknowns {est.unknowns} "
            f"degree {est.degree} ranks {est.ranks} rows {est.rows_used}"
        )
        _say(
            f"[WF-13 s11 U7-d2] predicted_seconds {est.predicted_seconds:.3f} "
            f"predicted_peak_gib {est.predicted_peak_gib:.4f} at -n {est.ranks} "
            f"(printed only; no -n 2 degree-2 F-human solve to compare)"
        )
        _say(
            f"[WF-13 s11 CLI] rc {rcs} diff {diff} not_priced {est.not_priced} "
            f"output dir exists {out_exists}"
        )
        _say(
            f"[WF-13 s11 U7-d2] real builds per rank {builds}, memo hits {hits}; "
            f"node {time.perf_counter() - t_start:.1f} s"
        )
    assert est.unknowns == N_DOFS_D2_RECORD
    assert est.cells == HOLE_CELL_RECORD
    assert est.degree == 2
    assert "FH2" in est.rows_used
    assert est.not_priced == ("thermal transient",)
    assert study.result is None
    assert rcs == [0] * comm.size
    assert diff == set()
    assert out_exists == [False] * comm.size
    assert builds == [1] * comm.size
    assert hits == [1] * comm.size
