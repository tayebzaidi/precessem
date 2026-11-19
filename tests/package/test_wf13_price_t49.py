"""`WF-13` step 7 (§9 item 479, chain U55a): the tune YAML's degree-1 price row T49.

Brief: docs/planning/chunks/WF-13.md, "Draft brief — §9 item 479 (U55a)" and
"Review ruling 2026-10-15 18:00" (the ruling overrides the draft).  Wm (this
node, run before any ``PRICE_TABLE`` edit) solves ``TUNE_FILE`` (the F-small
exact-rectangle PEC hole, degree 1, ``-n 2``) through ``run_study_file`` and
prints the N1curl unknown count and the summed peak RSS that a T49 row would
carry.  Backed fixture checks: ``cells == ANS6_EXACT_STRIP_HOLE_CELL_RECORD``
exactly and ``|C_mode / TH17_STEP3_C_MODE_RECORD_F - 1| <= REDUCTION_BAND``.
The leave-one-out reading (the tune file priced from P1 / P3 alone, predicted
ratio ≈ 0.53 against item 442's 29.668 s) is printed, never asserted.
"""

from __future__ import annotations

import time

from mpi4py import MPI

from tests.complex_mode import complex_only
from tests.package.test_wf1_study_yaml import LOGS, _fresh_dir
from tests.package.test_wf1_yaml_tune import TUNE_FILE
from tests.package.test_wf16_study_tune import TH17_STEP3_C_MODE_RECORD_F
from tests.validation.test_ans6_exact_strip import ANS6_EXACT_STRIP_HOLE_CELL_RECORD
from tests.validation.test_port_lumped_rlc_termination import REDUCTION_BAND

# 20261014T110326Z_WF-1-s11-W1.log:1830 — item 442's run_study_file on TUNE_FILE, -n 2.
ITEM442_RUN_SECONDS = 29.668


def _p0(msg):
    if MPI.COMM_WORLD.rank == 0:
        print(msg, flush=True)


@complex_only
def test_wf13_t49_measure():
    """Wm: the tune YAML solved; N, summed peak and the leave-one-out price printed."""
    from dolfinx import fem

    from fem_em_solver import load_study, run_study_file
    from fem_em_solver.utils.instrumentation import report_peak_rss

    comm = MPI.COMM_WORLD
    assert comm.size == 2, "T49 is a -n 2 row"
    est = load_study(TUNE_FILE).estimate(comm)
    _p0(
        f"\n[WF-13 T49 loo] P1/P3-only estimate {est!r}; ratio predicted / item 442 "
        f"{est.predicted_seconds / ITEM442_RUN_SECONDS:.4f} (printed, never asserted)"
    )
    d = LOGS / "wf13_t49_measure"
    _fresh_dir(d, comm)
    t0 = time.perf_counter()
    run = run_study_file(TUNE_FILE, output_dir=d, comm=comm)
    t_run = time.perf_counter() - t0
    res = run.result
    V = fem.functionspace(res.mesh, ("N1curl", 1))
    n = int(V.dofmap.index_map.size_global * V.dofmap.index_map_bs)
    c = run.summary["tune"]["c_mode_farad"]
    gap = c / TH17_STEP3_C_MODE_RECORD_F - 1.0
    total = report_peak_rss(comm, label="[WF-13 T49] summed peak ru_maxrss")
    _p0(
        f"\n[WF-13 T49] cells {res.cells}; unknowns {n}; mesh {res.mesh_time:.2f} s, "
        f"sweep {res.sweep_time:.2f} s; run_study_file {t_run:.3f} s; "
        f"peak {total / 2**30:.4f} GiB summed over {comm.size} ranks; "
        f"C_mode gap {gap:+.6e} (band {REDUCTION_BAND:.1e})"
    )
    assert res.cells == ANS6_EXACT_STRIP_HOLE_CELL_RECORD
    assert abs(gap) <= REDUCTION_BAND


def _t49(table=None):
    from fem_em_solver.study import PRICE_TABLE

    return [r for r in (table or PRICE_TABLE) if r["row"] == "T49"][0]


@complex_only
def test_wf13_t49_u7_e_d2():
    """W2: (U7) N exact vs Wm; (E) ×2 price (by construction); (D2) degree 2 unmoved."""
    from fem_em_solver import load_study
    from tests.package.test_wf4_price_rows import D2_ESTIMATE_RECORD
    from tests.package.test_wf13_cli_estimate import DEG2_FILE
    from tests.package.test_wf13_study_estimate import BRACKET
    from tests.validation.test_birdcage_sar_mass_averaged import EXACT_IDENTITY_RTOL

    comm = MPI.COMM_WORLD
    assert comm.size == 2, "T49 is a -n 2 row"
    t49 = _t49()
    est = load_study(TUNE_FILE).estimate(comm)
    ratio = est.predicted_seconds / ITEM442_RUN_SECONDS
    _p0(
        f"\n[WF-13 T49 U7/E] {est!r}; T49 {t49}; (E, by construction) ratio {ratio:.4f} "
        f"(BRACKET {BRACKET})"
    )
    # (U7) 20261016T035745Z_WF-13-s7-Wm.log:1788 measured N 111 946, cells 80 798.
    assert est.unknowns == t49["unknowns"] == 111946
    assert est.cells == t49["cells"] == ANS6_EXACT_STRIP_HOLE_CELL_RECORD
    assert "T49" in est.rows_used
    # (E) by construction: T49 sits inside the three-row degree-1 fit (rule (h)).
    assert 1.0 / BRACKET <= ratio <= BRACKET
    d2 = load_study(DEG2_FILE).estimate(comm)
    rec = D2_ESTIMATE_RECORD
    dp = abs(d2.predicted_peak_gib / rec["predicted_peak_gib"] - 1.0)
    ds = abs(d2.predicted_seconds / rec["predicted_seconds"] - 1.0)
    _p0(
        f"[WF-13 T49 D2] {d2!r}; rel peak {dp:.3e}, rel seconds {ds:.3e} "
        f"(EXACT_IDENTITY_RTOL {EXACT_IDENTITY_RTOL:.0e})"
    )
    # FH2: §9 item 585 — D2_ESTIMATE_RECORD re-recorded under (1*) to the
    # ten-row table (review ruling 2026-10-22 03:00 (2)); the move is printed.
    from tests.package.test_wf4_price_rows import D2_ESTIMATE_RECORD_NINE_ROW as old

    _p0(
        f"[WF-13 T49 D2] FH2 move vs nine-row record: rel peak "
        f"{abs(d2.predicted_peak_gib / old['predicted_peak_gib'] - 1.0):.3e}, rel seconds "
        f"{abs(d2.predicted_seconds / old['predicted_seconds'] - 1.0):.3e}"
    )
    for k in ("cells", "unknowns", "degree", "ranks"):
        assert getattr(d2, k) == rec[k]
    assert "T49" not in d2.rows_used
    assert dp <= EXACT_IDENTITY_RTOL
    assert ds <= EXACT_IDENTITY_RTOL


@complex_only
def test_wf13_t49_controls():
    """Controls (by construction): T49 N ×1.1 moves the price; a P1-only table refuses."""
    import pytest

    from fem_em_solver.study import PRICE_TABLE, _predict
    from tests.package.test_wf13_cli_estimate import DEG1_MESSAGE
    from tests.validation.test_birdcage_sar_mass_averaged import EXACT_IDENTITY_RTOL

    t49 = _t49()
    _, sec0, _ = _predict(t49["unknowns"], 1, 2)
    bumped = tuple(
        dict(r, unknowns=int(r["unknowns"] * 1.1)) if r["row"] == "T49" else r for r in PRICE_TABLE
    )
    _, sec1, _ = _predict(t49["unknowns"], 1, 2, table=bumped)
    c1 = abs(sec1 / sec0 - 1.0)
    _p0(f"\n[WF-13 T49 c1] T49 unknowns x1.1: seconds {sec1!r} vs {sec0!r} (rel {c1:.3e})")
    assert c1 > EXACT_IDENTITY_RTOL
    p1_only = tuple(
        r for r in PRICE_TABLE if r["row"] not in ("P3", "T49", "FH1")
    )  # FH1: §9 item 562
    with pytest.raises(ValueError) as exc:
        _predict(t49["unknowns"], 1, 2, table=p1_only)
    _p0(f"[WF-13 T49 c2] P1-only table: {exc.value}")
    assert str(exc.value) == DEG1_MESSAGE
