"""`WF-13` step 10 (§9 item 562, chain U68): the F-human estimate branch.

``Study.estimate`` / ``femem run --estimate`` on the three tracked F-human
files (`WF-16`'s refusal lifted for ``estimate`` only) price the EM
factorisation from row FH1 at degree 1 and state the thermal transient
unpriced.  FH1's N is `WF-10`'s logged edge count.

Brief: docs/planning/chunks/WF-13.md, "§9 item 562".  Run under ``mpiexec -n 2``.

Anchors (asserted; bands imported by name, none set here):
  (FH)  FH1's six numbers equal the transcribed literals; unknowns == E_RECORD,
        cells == HOLE_CELL_RECORD.
  (LOO) ``_predict`` at FH1's N and width on the table without FH1: the
        seconds ratio inside ``[1/BRACKET, BRACKET]``; rows used P1, P3, T49.
  (B)   with the F-human ``build`` patched to raise, ``estimate`` reaches it
        (not ``TypeError``) on all three files; ``_refuse_f_human`` still
        refuses ``mean_sar``.
  (U7)  (knob ``FEM_EM_WF13_S10=1``, one real build, memoized) per file:
        unknowns == E_RECORD, cells == HOLE_CELL_RECORD, degree 1, "FH1" used,
        no solve; (T) ``not_priced == ("thermal transient",)``; (CLI = API)
        rc [0, 0], printed fields equal the API's, no output directory.

Negative controls (asserted, by construction): (c-N) an N-blind estimator
misses FH1's seconds by more than ×2; (c-row) FH1's unknowns ×1.1 moves the
in-sample price at FH1's N by more than ``EXACT_IDENTITY_RTOL``.

Printed only: the LOO memory ratio (predicted outside ×2, a disclosed finding),
the in-sample ratios with FH1 fitted, the degree-2 fit at FH1's N and width,
and the ``-n 2`` predictions per file.
"""

from __future__ import annotations

import dataclasses
import os
import time

import numpy as np
import pytest
from mpi4py import MPI

from fem_em_solver.study import PRICE_TABLE, FHumanBirdcageGeometry, _predict, _refuse_f_human
from tests.package.test_wf1_study_yaml import WF1
from tests.package.test_wf13_cli_estimate import _diff, _parse, _patch_solves

KNOB = "FEM_EM_WF13_S10"
FILES = (
    WF1 / "f_human_hole_copper_tune.yaml",
    WF1 / "f_human_tuned_copper.yaml",
    WF1 / "f_human_tuned_copper_thermal.yaml",
)
THERMAL = WF1 / "f_human_tuned_copper_thermal.yaml"
# FH1 transcribed: 20261015T213534Z_WF-16-s12-W1.log:7343 (seconds), :7344
# (memory); N 20261011T140313Z_WF-10-step4a-W1.log:6215; cells = the hole record.
FH1_LITERALS = dict(
    row="FH1", degree=1, ranks=8, cells=369387, unknowns=500568, seconds=93.4, peak_gib=9.7652
)


class _BuildReached(Exception):
    pass


def _say(msg):
    if MPI.COMM_WORLD.rank == 0:
        print(msg, flush=True)


def _row(name):
    rows = [r for r in PRICE_TABLE if r["row"] == name]
    assert len(rows) == 1, name
    return rows[0]


def test_wf13_s10_fh1_loo_and_controls():
    """(FH), (LOO), (c-N), (c-row); printed readings.  No mesh."""
    from tests.package.test_wf13_study_estimate import BRACKET
    from tests.validation.test_birdcage_sar_mass_averaged import EXACT_IDENTITY_RTOL
    from tests.validation.test_wf10_f_human_tuned_b1_sar import E_RECORD
    from tests.validation.test_wf12_f_human_copper_matrix import HOLE_CELL_RECORD

    fh1 = _row("FH1")
    # (FH)
    assert {k: fh1[k] for k in FH1_LITERALS} == FH1_LITERALS
    assert set(fh1) == set(FH1_LITERALS)
    assert fh1["unknowns"] == E_RECORD
    assert fh1["cells"] == HOLE_CELL_RECORD
    n, w = fh1["unknowns"], fh1["ranks"]

    # (LOO)
    loo = tuple(r for r in PRICE_TABLE if r["row"] != "FH1")
    m_loo, s_loo, used = _predict(n, 1, w, table=loo)
    rt = s_loo / fh1["seconds"]
    rm = m_loo / fh1["peak_gib"]
    _say(
        f"\n[WF-13 s10 LOO] FH1 at N={n}, -n {w}: predicted {s_loo:.3f} s vs "
        f"{fh1['seconds']} s -> ratio {rt:.4f} (asserted in [{1/BRACKET}, {BRACKET}]); "
        f"rows {used}"
    )
    _say(
        f"[WF-13 s10 LOO] memory (printed only): {m_loo:.4f} GiB vs {fh1['peak_gib']} GiB "
        f"-> ratio {rm:.4f} ({'inside' if 1/BRACKET <= rm <= BRACKET else 'OUTSIDE'} x2)"
    )
    assert used == ("P1", "P3", "T49")
    assert 1.0 / BRACKET <= rt <= BRACKET

    # Printed: in-sample with FH1 fitted, and the degree-2 fit at FH1's N / width.
    m_in, s_in, used_in = _predict(n, 1, w)
    _say(
        f"[WF-13 s10 in-sample] rows {used_in}: FH1 time {s_in / fh1['seconds']:.4f} "
        f"memory {m_in / fh1['peak_gib']:.4f}"
    )
    t49 = _row("T49")
    _, s49, _ = _predict(t49["unknowns"], 1, t49["ranks"])
    _say(f"[WF-13 s10 in-sample] T49 time {s49 / t49['seconds']:.4f}")
    m2, s2, used2 = _predict(n, 2, w)
    _say(
        f"[WF-13 s10 degree-2 fit at FH1 N/width] {s2:.3f} s -> ratio "
        f"{s2 / fh1['seconds']:.4f} ({'inside' if 1/BRACKET <= s2 / fh1['seconds'] <= BRACKET else 'outside'} x2; "
        f"degree confusion is not a control here); memory {m2:.3f} GiB; rows {used2}"
    )

    # (c-N) N-blind estimator, the three small rows.
    small = [r for r in loo if r["degree"] == 1]
    gm = float(np.exp(np.mean(np.log([r["seconds"] * r["ranks"] for r in small]))))
    rn = (gm / w) / fh1["seconds"]
    _say(f"[WF-13 s10 c-N] N-blind {gm / w:.3f} s vs {fh1['seconds']} s -> ratio {rn:.4f}")
    assert not (1.0 / BRACKET <= rn <= BRACKET)

    # (c-row) FH1 is load-bearing.
    bumped = tuple(
        dict(r, unknowns=int(r["unknowns"] * 1.1)) if r["row"] == "FH1" else r for r in PRICE_TABLE
    )
    _, s_b, _ = _predict(n, 1, w, table=bumped)
    rel = abs(s_b / s_in - 1.0)
    _say(f"[WF-13 s10 c-row] FH1 unknowns x1.1: {s_b!r} vs {s_in!r} (rel {rel:.3e})")
    assert rel > EXACT_IDENTITY_RTOL


def test_wf13_s10_branch_reaches_build(monkeypatch):
    """(B): the refusal is lifted before any mesh; other refusals stand.  No mesh."""
    from fem_em_solver.study_yaml import load_study

    comm = MPI.COMM_WORLD

    def _raise(*a, **k):
        raise _BuildReached()

    monkeypatch.setattr(FHumanBirdcageGeometry, "build", _raise)
    for f in FILES:
        study = load_study(f)
        assert isinstance(study.geometry, FHumanBirdcageGeometry)
        with pytest.raises(_BuildReached):
            study.estimate(comm)
        with pytest.raises(TypeError):
            _refuse_f_human(study, "mean_sar")
        _say(f"[WF-13 s10 B] {f.name}: estimate reached build; mean_sar still refused")


def test_wf13_s10_f_human_files_priced(monkeypatch, capsys, tmp_path):
    """(U7), (T), (CLI = API) on the three files; one real F-human build."""
    if os.environ.get(KNOB) != "1":
        pytest.skip(f"{KNOB}=1 runs the one-build window")
    from fem_em_solver.cli import main
    from fem_em_solver.study import ResourceEstimate
    from fem_em_solver.study_yaml import load_study
    from tests.validation.test_wf10_f_human_tuned_b1_sar import E_RECORD
    from tests.validation.test_wf12_f_human_copper_matrix import HOLE_CELL_RECORD

    comm = MPI.COMM_WORLD
    assert comm.size == 2
    _patch_solves(monkeypatch)

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
        _say(f"[WF-13 s10 U7] real F-human build {time.perf_counter() - t0:.1f} s")
        return memo[key]

    monkeypatch.setattr(FHumanBirdcageGeometry, "build", _memo_build)
    names = {f.name for f in dataclasses.fields(ResourceEstimate)}
    t_start = time.perf_counter()
    for f in FILES:
        study = load_study(f)
        est = study.estimate(comm)
        out_dir = tmp_path / f"out-{f.stem}"
        capsys.readouterr()
        rc = main(["run", "--estimate", str(f), "--output-dir", str(out_dir)])
        out = capsys.readouterr().out
        rcs = comm.allgather(rc)
        diff = _diff(_parse(out, names), est) if comm.rank == 0 else set()
        diff = comm.bcast(diff, root=0)
        has_line = comm.bcast(
            "not_priced: ('thermal transient',)" in out.splitlines() if comm.rank == 0 else None,
            root=0,
        )
        _say(
            f"\n[WF-13 s10 U7] {f.name}: cells {est.cells} unknowns {est.unknowns} "
            f"degree {est.degree} ranks {est.ranks} rows {est.rows_used}"
        )
        _say(
            f"[WF-13 s10 U7] {f.name}: predicted_seconds {est.predicted_seconds:.3f} "
            f"predicted_peak_gib {est.predicted_peak_gib:.4f} at -n {est.ranks} "
            f"(printed only; no -n 2 F-human solve to compare)"
        )
        _say(
            f"[WF-13 s10 CLI] {f.name}: rc {rcs} diff {diff} not_priced {est.not_priced} "
            f"stdout line present {has_line}"
        )
        assert est.unknowns == E_RECORD
        assert est.cells == HOLE_CELL_RECORD
        assert est.degree == 1
        assert "FH1" in est.rows_used
        assert study.result is None
        assert est.not_priced == ("thermal transient",)
        assert rcs == [0] * comm.size
        assert diff == set()
        assert not out_dir.exists()
        if f == THERMAL:
            assert has_line
    builds = comm.allgather(count["builds"])
    hits = comm.allgather(count["hits"])
    _say(
        f"[WF-13 s10 U7] real builds per rank {builds}, memo hits {hits}; "
        f"node {time.perf_counter() - t_start:.1f} s"
    )
    assert builds == [1] * comm.size
    assert hits == [5] * comm.size
