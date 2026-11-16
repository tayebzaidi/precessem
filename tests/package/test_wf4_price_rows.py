"""`WF-13` step 6 (§9 item 471a, chain U42, `WF-4` folded in): the degree-1 price row.

Brief: docs/planning/chunks/WF-4.md, "§9 item 471" part 471a.  ``PRICE_TABLE``
P3 (item 192's exact-rectangle F-small, degree 1, ``-n 2`` — the `ANS-4`
fixture `WF-4`'s ``dev`` preset uses) was time-only; this module measures its
N1curl unknown count and summed peak RSS on the Study route (Wm), and after the
row is completed from that log checks (W2): (U3) the estimator's N equals the
row's; (Pd) ``femem run --estimate`` prices both degree-1 files and equals the
API field by field; (Int) the two-row power law passes through both rows
(``EXACT_IDENTITY_RTOL`` — by construction, so the degree-1 ``×2`` checks of
item 471b are reproductions, not tests of a fit); (D2) the degree-2 estimate
is unmoved against ``20261008T200148Z_WF-13-s4-W1.log:212``.

Controls (asserted, by construction): (c1) a copy of the table with P3's
unknowns ×1.1 moves the F-small price off P3's seconds; (c2) a P1-only table
still refuses degree 1 with the API message.
"""

from __future__ import annotations

import dataclasses

import pytest
from mpi4py import MPI

from tests.complex_mode import complex_only
from tests.package.test_wf13_cli_estimate import DEG1_MESSAGE, DEG2_FILE, _diff, _parse, _patch_solves
from tests.package.test_wf16_study_four_port import _record, _rel
from tests.package.test_wf1_study_yaml import BIRDCAGE_FILE, SPHERE_FILE
from tests.validation.test_birdcage_sar_mass_averaged import EXACT_IDENTITY_RTOL
from tests.validation.test_port_birdcage_four_port import LEG_D0_REPRODUCTION_BAND

# 20261008T200148Z_WF-13-s4-W1.log:212 — the degree-2 sphere estimate at -n 2
# before this item (a record of the same comparison on the same fixture/width).
# FH2: §9 item 585 — re-recorded under (1*) to the ten-row table (review ruling
# 2026-10-22 03:00 (2)): FH2 joins the degree-2 fit and moves this estimate by
# rel 2.226e-03 (peak) / 8.595e-03 (seconds); the nine-row values are kept
# below and the move is printed by both D2 nodes.  No band moves.
D2_ESTIMATE_RECORD_NINE_ROW = dict(
    predicted_peak_gib=15.361991136892966,
    predicted_seconds=246.50021705192188,
)
D2_ESTIMATE_RECORD = dict(
    cells=74020,
    unknowns=482158,
    degree=2,
    ranks=2,
    predicted_peak_gib=15.32780177373478,
    predicted_seconds=244.38150702034363,
)


def _p0(msg):
    if MPI.COMM_WORLD.rank == 0:
        print(msg, flush=True)


def _row(name, table=None):
    from fem_em_solver.study import PRICE_TABLE

    return [r for r in (table or PRICE_TABLE) if r["row"] == name][0]


def _width2():
    assert MPI.COMM_WORLD.size == 2, "P3 is a -n 2 row"


def _two_row(monkeypatch):
    """Item 479 re-point: the P1 / P3 construction, T49 removed (``_predict`` patched)."""
    import fem_em_solver.study as study_mod

    orig = study_mod._predict
    two = tuple(r for r in study_mod.PRICE_TABLE if r["row"] not in ("T49", "FH1"))  # FH1: §9 item 562
    monkeypatch.setattr(study_mod, "_predict", lambda u, d, r, n_rhs=1: orig(u, d, r, n_rhs, table=two))


@complex_only
def test_wf4_p3_measure():
    """Wm: the F-small degree-1 Study solve; N and summed peak printed for P3."""
    import numpy as np
    from dolfinx import fem

    from fem_em_solver.study_yaml import load_study
    from fem_em_solver.utils.instrumentation import report_peak_rss

    _width2()
    comm = MPI.COMM_WORLD
    study = load_study(BIRDCAGE_FILE)
    assert study.degree == 1
    res = study.solve(comm)
    V = fem.functionspace(res.mesh, ("N1curl", 1))
    n = int(V.dofmap.index_map.size_global * V.dofmap.index_map_bs)
    s = np.asarray(res.s_matrix, dtype=np.complex128).reshape(4, 4)
    r1 = _rel(s, _record(1))
    total = report_peak_rss(comm, label="[WF-4 P3] summed peak ru_maxrss")
    _p0(
        f"\n[WF-4 P3] cells {res.cells} (P3 row {_row('P3')['cells']}); "
        f"unknowns {n}; mesh {res.mesh_time:.2f} s, sweep {res.sweep_time:.2f} s; "
        f"peak {total / 2**30:.2f} GiB summed over {comm.size} ranks; "
        f"S rel _record(1) {r1:.3e} (band {LEG_D0_REPRODUCTION_BAND:.0e})"
    )
    assert res.cells == _row("P3")["cells"]
    assert r1 <= LEG_D0_REPRODUCTION_BAND


@complex_only
def test_wf4_u3_unknowns_equal_row():
    """(U3) the estimator's N and cells equal P3's, exactly."""
    from fem_em_solver.study_yaml import load_study

    _width2()
    est = load_study(BIRDCAGE_FILE).estimate(MPI.COMM_WORLD)
    p3 = _row("P3")
    _p0(f"\n[WF-4 U3] estimate {est!r}; P3 cells {p3['cells']} unknowns {p3['unknowns']}")
    assert est.unknowns == p3["unknowns"]
    assert est.cells == p3["cells"]


@complex_only
def test_wf4_pd_cli_prices_degree1(monkeypatch, capsys, tmp_path):
    """(Pd) both degree-1 files priced through the CLI, equal to the API."""
    from fem_em_solver.cli import main
    from fem_em_solver.study import ResourceEstimate
    from fem_em_solver.study_yaml import load_study

    _width2()
    _two_row(monkeypatch)
    comm = MPI.COMM_WORLD
    names = [f.name for f in dataclasses.fields(ResourceEstimate)]
    apis = {f: load_study(f).estimate(comm) for f in (SPHERE_FILE, BIRDCAGE_FILE)}
    _patch_solves(monkeypatch)
    for f, api in apis.items():
        out_dir = tmp_path / f.stem
        capsys.readouterr()
        rc = main(["run", "--estimate", str(f), "--output-dir", str(out_dir)])
        out = capsys.readouterr().out
        verdict = None
        if comm.rank == 0:
            parsed = _parse(out, names)
            complete = set(parsed) == set(names) and all(len(v) == 1 for v in parsed.values())
            verdict = dict(complete=complete, differ=_diff(parsed, api) if complete else None)
            with capsys.disabled():  # the next readouterr() would swallow it
                print(f"\n[WF-4 Pd] {f.name}: rc {rc}; API {api!r}; {verdict}", flush=True)
        verdict = comm.bcast(verdict, root=0)
        assert comm.allgather(rc) == [0] * comm.size
        assert verdict["complete"]
        assert verdict["differ"] == set()
        assert api.rows_used == ("P1", "P3")
        assert not out_dir.exists()


@complex_only
def test_wf4_int_two_row_interpolation(monkeypatch):
    """(Int) the two-row power law reproduces each degree-1 row (by construction)."""
    from fem_em_solver.study_yaml import load_study

    _width2()
    _two_row(monkeypatch)
    comm = MPI.COMM_WORLD
    for f, name in ((SPHERE_FILE, "P1"), (BIRDCAGE_FILE, "P3")):
        est = load_study(f).estimate(comm)
        row = _row(name)
        dt = abs(est.predicted_seconds / row["seconds"] - 1.0)
        dm = abs(est.predicted_peak_gib / row["peak_gib"] - 1.0)
        _p0(
            f"\n[WF-4 Int] {f.name} vs {name}: seconds {est.predicted_seconds!r} vs "
            f"{row['seconds']} (rel {dt:.3e}); peak {est.predicted_peak_gib!r} vs "
            f"{row['peak_gib']} (rel {dm:.3e}); band {EXACT_IDENTITY_RTOL:.0e}"
        )
        assert dt <= EXACT_IDENTITY_RTOL
        assert dm <= EXACT_IDENTITY_RTOL


@complex_only
def test_wf4_d2_degree2_fit_unmoved():
    """(D2) the degree-2 sphere estimate equals the step-4 log's, floats to rtol."""
    from fem_em_solver.study_yaml import load_study

    _width2()
    est = load_study(DEG2_FILE).estimate(MPI.COMM_WORLD)
    rec = D2_ESTIMATE_RECORD
    dp = abs(est.predicted_peak_gib / rec["predicted_peak_gib"] - 1.0)
    ds = abs(est.predicted_seconds / rec["predicted_seconds"] - 1.0)
    _p0(f"\n[WF-4 D2] {est!r}; rel peak {dp:.3e}, rel seconds {ds:.3e}")
    old = D2_ESTIMATE_RECORD_NINE_ROW  # FH2: §9 item 585 — the move, printed
    _p0(f"[WF-4 D2] FH2 move vs nine-row record: rel peak "
        f"{abs(est.predicted_peak_gib / old['predicted_peak_gib'] - 1.0):.3e}, rel seconds "
        f"{abs(est.predicted_seconds / old['predicted_seconds'] - 1.0):.3e}")
    for k in ("cells", "unknowns", "degree", "ranks"):
        assert getattr(est, k) == rec[k]
    assert est.rows_used[0] == "P2" and "P1" not in est.rows_used and "P3" not in est.rows_used
    assert dp <= EXACT_IDENTITY_RTOL
    assert ds <= EXACT_IDENTITY_RTOL


@complex_only
def test_wf4_controls():
    """(c1) the interpolation reads P3; (c2) a P1-only table still refuses degree 1."""
    from fem_em_solver.study import PRICE_TABLE, _predict

    _width2()
    p3 = _row("P3")
    bumped = tuple(dict(r, unknowns=int(r["unknowns"] * 1.1)) if r["row"] == "P3" else r for r in PRICE_TABLE)
    _, sec, _ = _predict(p3["unknowns"], 1, p3["ranks"], table=bumped)
    c1 = abs(sec / p3["seconds"] - 1.0)
    _p0(f"\n[WF-4 c1] P3 unknowns x1.1: F-small seconds {sec!r} vs {p3['seconds']} (rel {c1:.3e})")
    assert c1 > EXACT_IDENTITY_RTOL
    p1_only = tuple(r for r in PRICE_TABLE if r["row"] not in ("P3", "T49", "FH1"))  # T49: item 479; FH1: §9 item 562
    with pytest.raises(ValueError) as exc:
        _predict(p3["unknowns"], 1, p3["ranks"], table=p1_only)
    _p0(f"[WF-4 c2] P1-only table: {exc.value}")
    assert str(exc.value) == DEG1_MESSAGE
