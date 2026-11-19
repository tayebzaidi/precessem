"""`WF-13` step 2 (§9 item 268, re-posed as item 287): ``Study.estimate()``.

(U) unknowns exact against the record / the solve's own space; (U4) the Study
route builds P4's mesh; (K) no solve; (B) leave-one-out ×2 bracket on P2 / P5 /
P6 over the nine-row degree-2 table; (c1) degree confusion; (c2) a planted
price on the row that sets P6's level (×8); (c3) an N-blind estimator.  No band
is set here: ×2 is the brief's (items 268 / 287), never widened.
"""

from __future__ import annotations

import numpy as np
from mpi4py import MPI

from fem_em_solver import FourPortBirdcageGeometry, LumpedPortDrive, Study
from fem_em_solver.materials import GelledSalinePhantomMaterial
from fem_em_solver.study import PRICE_TABLE, _predict
from tests.complex_mode import complex_only
from tests.mesh.test_birdcage_port_sheet_prerequisite import CONDUCTOR_RESOLUTION
from tests.package.test_wf15_study_lossy_sphere import _study
from tests.package.test_wf16_study_four_port import (
    FREQUENCY_64_HZ,
    SALINE_EPSILON_R,
    SALINE_SIGMA,
)
from tests.validation.test_ops61_exact_strip_ports import ANS4_EXACT_STRIP_CELL_RECORD

# `20261005T170537Z_WF-15-r2-W2a-deg2-study.log:158`: the degree-2 sphere Study.
WF15_DEG2_UNKNOWNS_RECORD = 482158
WF15_CELLS_RECORD = 74020
BRACKET = 2.0  # item 268 (B): 0.5 <= pred/rec <= 2, never widened
LOO_ROWS = ("P2", "P5", "P6")
# `20261006T185702Z_WF-13-s2-W1.log:1989`: F-small Study, degree 2.
F_SMALL_DEG2_UNKNOWNS_RECORD = 753062
# Item 225's printed count (`20261004T183338Z_ANS-4-item225-W3-p2-cf0p5-128-koff.log:2083`).
P4_CELLS_RECORD = 259524
# (c2) item 287: P5 carries the largest hat weight in P6's LOO fit; minimum
# failing k ≈ 4.0 (mem) / 2.9 (time) on nine rows, so x8 is twice it.
PLANT_K = 8.0
# Item 212: P4's mesh at 10 / 64 MHz — duplicates of P4's fixture, never fitted.
ITEM212_RANK_ROWS = (
    dict(row="212-W1", ranks=8, peak_gib=36.12, seconds=186.93),
    dict(row="212-W2", ranks=12, peak_gib=38.40, seconds=183.12),
)


def _say(msg):
    if MPI.COMM_WORLD.rank == 0:
        print(msg, flush=True)


@complex_only
def test_wf13_estimate_sphere_unknowns_exact_and_no_solve():
    comm = MPI.COMM_WORLD
    s2 = _study()
    assert s2.degree == 2
    e2 = s2.estimate(comm)
    _say(
        f"WF13-S2 (U) sphere deg2 cells={e2.cells} unknowns={e2.unknowns} "
        f"pred_peak={e2.predicted_peak_gib:.2f} GiB pred_s={e2.predicted_seconds:.1f} "
        f"ranks={e2.ranks} rows={e2.rows_used}"
    )
    assert s2.result is None  # (K)
    assert e2.cells == WF15_CELLS_RECORD
    assert e2.unknowns == WF15_DEG2_UNKNOWNS_RECORD
    # Degree 1 carries < 2 priced rows, so its count is read from the space the
    # Study's own solve() would build (TimeHarmonicSolver.function_space).
    from fem_em_solver.core import HomogeneousMaterial, TimeHarmonicProblem, TimeHarmonicSolver

    msh, cell_tags, _ = s2.geometry.build(comm)
    problem = TimeHarmonicProblem(
        mesh=msh,
        frequency_hz=s2.frequency_hz,
        material=HomogeneousMaterial(sigma=0.0, epsilon_r=1.0),
        cell_tags=cell_tags,
        phantom_material=s2.material,
        phantom_tag=s2.geometry.sphere_tag,
        boundary_condition="pec_zero_tangential_a",
        dirichlet_e_field=s2.dirichlet_field(),
    )
    V2 = TimeHarmonicSolver(problem, degree=2).function_space()
    assert V2.dofmap.index_map.size_global * V2.dofmap.index_map_bs == e2.unknowns
    V1 = TimeHarmonicSolver(problem, degree=1).function_space()
    n1 = V1.dofmap.index_map.size_global * V1.dofmap.index_map_bs
    _say(f"WF13-S2 (c1) sphere deg1 unknowns={n1} ratio deg2/deg1={e2.unknowns / n1:.4f}")
    assert n1 != WF15_DEG2_UNKNOWNS_RECORD  # (c1)


@complex_only
def test_wf13_estimate_f_small_cells():
    comm = MPI.COMM_WORLD
    phantom = GelledSalinePhantomMaterial(SALINE_SIGMA, SALINE_EPSILON_R, FREQUENCY_64_HZ)
    st = Study(
        FourPortBirdcageGeometry.f_small(), phantom, LumpedPortDrive(), frequency_hz=FREQUENCY_64_HZ
    )
    est = st.estimate(comm)
    _say(
        f"WF13-S2 (U) f_small deg{est.degree} cells={est.cells} unknowns={est.unknowns} "
        f"pred_peak={est.predicted_peak_gib:.2f} GiB pred_s={est.predicted_seconds:.1f}"
    )
    assert st.result is None
    assert est.cells == ANS4_EXACT_STRIP_CELL_RECORD
    assert est.unknowns == F_SMALL_DEG2_UNKNOWNS_RECORD


@complex_only
def test_wf13_p4_cells_and_unknowns_read():
    """(U4): the Study route builds P4's mesh (item 225's cells); prints P4's N."""
    comm = MPI.COMM_WORLD
    phantom = GelledSalinePhantomMaterial(SALINE_SIGMA, SALINE_EPSILON_R, FREQUENCY_64_HZ)
    geom = FourPortBirdcageGeometry(conductor_resolution=0.5 * CONDUCTOR_RESOLUTION)
    st = Study(geom, phantom, LumpedPortDrive(), frequency_hz=FREQUENCY_64_HZ)
    est = st.estimate(comm)
    _say(
        f"WF13-S2 (U4) P4 route deg{est.degree} cells={est.cells} "
        f"(record {P4_CELLS_RECORD}) unknowns={est.unknowns} "
        f"({est.unknowns / est.cells:.4f} per cell)"
    )
    assert st.result is None
    assert est.degree == 2
    assert est.cells == P4_CELLS_RECORD


def _row(table, name):
    return next(r for r in table if r["row"] == name)


def _deg2(table):
    return [r for r in table if r["degree"] == 2 and r["unknowns"]]


def _ratios(table, name, drop=()):
    rec = _row(table, name)
    rest = tuple(r for r in table if r["row"] != name and r["row"] not in drop)
    peak, secs, _ = _predict(rec["unknowns"], rec["degree"], rec["ranks"], table=rest)
    return peak / rec["peak_gib"], secs / rec["seconds"]


def _hat_weights(table, name):
    """Log-space hat weights of the LOO fit's prediction at ``name``'s N."""
    x0 = np.log(_row(table, name)["unknowns"])
    rest = [r for r in _deg2(table) if r["row"] != name]
    xs = np.log([r["unknowns"] for r in rest])
    xb = xs.mean()
    sxx = float(((xs - xb) ** 2).sum())
    w = 1.0 / len(xs) + (xs - xb) * (x0 - xb) / sxx
    return {r["row"]: float(wi) for r, wi in zip(rest, w)}


def _inside(r):
    return 1.0 / BRACKET <= r <= BRACKET


def test_wf13_price_leave_one_out_bracket():
    rows = _deg2(PRICE_TABLE)
    _say(f"WF13-S2 table: {len(rows)} degree-2 rows fitted: " f"{[r['row'] for r in rows]}")
    for r in PRICE_TABLE:
        _say(f"WF13-S2 row {r}")
    for r in ITEM212_RANK_ROWS:
        _say(f"WF13-S2 item-212 rank-sensitivity (not fitted) {r}")
    for r in rows:
        peak, secs, _ = _predict(r["unknowns"], 2, r["ranks"])
        _say(
            f"WF13-S2 in-sample {r['row']} mem {peak / r['peak_gib']:.3f} "
            f"time {secs / r['seconds']:.3f}"
        )
    for r in rows:
        if r["row"] in LOO_ROWS:
            continue
        rm, rt = _ratios(PRICE_TABLE, r["row"])
        _say(f"WF13-S2 LOO (printed) {r['row']} mem {rm:.3f} time {rt:.3f}")
    rm, rt = _ratios(PRICE_TABLE, "P5", drop=("A6",))
    _say(f"WF13-S2 LOO (printed) P5 with A6 removed mem {rm:.3f} time {rt:.3f}")
    bad = []
    for name in LOO_ROWS:
        w = _hat_weights(PRICE_TABLE, name)
        _say(f"WF13-S2 hat weights LOO {name}: " + " ".join(f"{k}={v:.3f}" for k, v in w.items()))
        rm, rt = _ratios(PRICE_TABLE, name)
        _say(f"WF13-S2 (B) LOO {name} mem pred/rec={rm:.3f} time pred/rec={rt:.3f}")
        for q, r in (("mem", rm), ("time", rt)):
            if not _inside(r):
                bad.append(f"{name} {q} {r:.3f}")
    assert not bad, f"(B) outside x2: {bad}"


def test_wf13_planted_price_fails_bracket():
    """(c2): P5 (largest hat weight in P6's LOO fit) memory and seconds x8."""
    # FH2: §9 item 585 — the nine-row construction this node and PLANT_K were
    # written on; the live ten-row table's planted control is
    # test_wf13_estimate_f_human_d2.py (c-plant).
    NO_FH2 = tuple(r for r in PRICE_TABLE if r["row"] != "FH2")
    w = _hat_weights(NO_FH2, "P6")["P5"]
    rm0, rt0 = _ratios(NO_FH2, "P6")
    # held-out ratio moves by k**w: the minimum failing scale per quantity
    kmin_m = (BRACKET / rm0) ** (1.0 / w) if rm0 < BRACKET else 1.0
    kmin_t = (BRACKET / rt0) ** (1.0 / w) if rt0 < BRACKET else 1.0
    planted = tuple(
        (
            dict(r, peak_gib=r["peak_gib"] * PLANT_K, seconds=r["seconds"] * PLANT_K)
            if r["row"] == "P5"
            else r
        )
        for r in NO_FH2
    )
    rm, rt = _ratios(planted, "P6")
    _say(
        f"WF13-S2 (c2) P5 hat weight in LOO P6 w={w:.3f}; min failing k mem "
        f"{kmin_m:.2f} time {kmin_t:.2f}; P5 x{PLANT_K:g} -> LOO P6 mem "
        f"{rm:.3f} time {rt:.3f}"
    )
    assert PLANT_K >= max(kmin_m, kmin_t)
    assert not _inside(rm)
    assert not _inside(rt)


def test_wf13_n_blind_estimator_fails_bracket():
    """(c3): a geometric-mean estimator with no N term fails (B) on every target."""
    rows = _deg2(PRICE_TABLE)
    bad = []
    for name in LOO_ROWS:
        rec = _row(rows, name)
        rest = [r for r in rows if r["row"] != name]
        mem = float(np.exp(np.mean(np.log([r["peak_gib"] for r in rest]))))
        rs = float(np.exp(np.mean(np.log([r["seconds"] * r["ranks"] for r in rest]))))
        rm = mem / rec["peak_gib"]
        rt = rs / (rec["seconds"] * rec["ranks"])
        _say(f"WF13-S2 (c3) N-blind LOO {name} mem {rm:.3f} time {rt:.3f}")
        for q, r in (("mem", rm), ("time", rt)):
            if _inside(r):
                bad.append(f"{name} {q} {r:.3f}")
    assert not bad, f"(c3) N-blind estimator inside x2: {bad}"
