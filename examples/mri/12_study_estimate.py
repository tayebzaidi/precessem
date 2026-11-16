"""Example (`EX-96`): price a ``Study`` before you solve it -- ``Study.estimate()``
on the lossy sphere, then the solve it priced (`WF-13` step 2 on `WF-15`'s sphere).

``study.estimate(comm)`` builds the mesh and the N1curl space and returns a
``ResourceEstimate``: cells, the **exact** global unknown count, degree, ranks, and
a priced peak memory and wall time from ``PRICE_TABLE``. This example builds
`WF-15`'s sphere ``Study`` at the **default degree** (2), prints the estimate,
solves the same ``Study`` and asserts `WF-13` step 2's identity: the solve's global
N1curl unknowns and cells equal the estimate's, exactly.

On this fixture the printed price is **in-sample**: the sphere at degree 2 and
``-n 2`` is ``PRICE_TABLE`` row P2. The gated, out-of-sample claim is `WF-13`'s
leave-one-out factor-of-2 bracket; (B-P2) re-asserts it on P2's row.

Every constant is imported from the gate modules (the ANS-1 rule); the example
states no sigma, epsilon_r, frequency, field, radius, unknown count or price of
its own.

Asserted:

* (K)   ``estimate()`` does not solve (``study.result is None``); after ``solve()``
        ``study.result is res``;
* (U)   estimate unknowns / cells == the solve's global N1curl unknowns / cells
        (integer ``==``), degree 2, ranks == comm.size, and both equal the
        imported records;
* (In)  ``"P2" in est.rows_used`` and P2's registered unknowns == ``est.unknowns``;
* (R)   |SAR - WF15_DEG2_RECORD| / record <= ``L0_BAND``;
* (G)   |SAR - closed form| / closed form <= ``MAT4_BAND``;
* (B-P2) the leave-one-out ratios on P2 (memory, time) sit inside ``BRACKET``.

Negative controls (all asserted):

* (c1) a degree-1 N1curl space on the solved mesh has a different global count
       than ``est.unknowns`` and equals P1's registered unknowns;
* (c2) on a P1-only table (re-pointed 2026-10-16 after 471a added degree-1 rows P3
       and T49) the degree-1 price raises ``ValueError`` ("fewer than two priced rows
       at degree 1"); asserted by construction (deterministic raise). The full table
       prices degree 1 since 471a; that price is printed, not asserted;
* (c3) the (R) comparator separates degrees: ``MAT4_DEG1_RECORD`` against
       ``WF15_DEG2_RECORD`` differs by more than ``L0_BAND`` (no solve).

Printed, not asserted (this machine): the estimate's predicted peak (GiB summed
over ranks) and seconds beside the measured summed ``ru_maxrss`` and solve
seconds, the estimate call's own seconds, the rows used.

Run it through the example runner (the ``mri:`` group sources the complex build)::

    ./scripts/run_examples.sh -e mri:12 -n 2 -t 590

Not claimed: a live-run memory or time accuracy for ``estimate()``; an
out-of-sample price on this fixture beyond `WF-13`'s leave-one-out factor of 2;
degree-1 pricing; a port ``Study``'s estimate; anything `WF-13` step 2 or `WF-15`
do not claim.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import ufl
from mpi4py import MPI

from dolfinx import fem

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from fem_em_solver import Study  # noqa: E402,F401  (the API under test, via _study)
from fem_em_solver.io.paraview_utils import (  # noqa: E402
    adopt_host_ownership,
    write_xdmf_with_tags,
)
from fem_em_solver.post.setup_figure import write_setup_figure  # noqa: E402
from fem_em_solver.study import PRICE_TABLE, _predict  # noqa: E402

from tests.materials.test_phantom_material_model import MAT4_BAND  # noqa: E402
from tests.package.test_wf13_study_estimate import (  # noqa: E402
    BRACKET,
    WF15_CELLS_RECORD,
    WF15_DEG2_UNKNOWNS_RECORD,
    _inside,
    _ratios,
)
from tests.package.test_wf15_study_lossy_sphere import (  # noqa: E402
    MAT4_DEG1_RECORD,
    _closed,
    _study,
    _summed_rss_gib,
)
from tests.package.test_wf17_first_study_tutorial import WF15_DEG2_RECORD  # noqa: E402
from tests.validation.test_mat5_em_coupled_sphere import L0_BAND  # noqa: E402

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"
BASENAME = "mri_12_study_estimate"
AIR_TAG = 2
# Printed beside the live reading (WF-13 step 2, 20261007T033535Z_WF-13-s2-r-W1.log:228).
WF13_DEG_RATIO_RECORD = 5.3449
# Printed beside (c3) (EX-91 control, 20261006T202424Z_EX-91.log:191).
EX91_C3_RECORD = 3.4937e-02


def _global_dofs(V) -> int:
    return int(V.dofmap.index_map.size_global * V.dofmap.index_map_bs)


def _row(name):
    return next(r for r in PRICE_TABLE if r["row"] == name)


def main() -> None:
    comm = MPI.COMM_WORLD
    started = time.perf_counter()

    # ---- price it first (no solve) ------------------------------------------
    study = _study()  # no degree= keyword
    t0 = time.perf_counter()
    est = study.estimate(comm)
    estimate_seconds = time.perf_counter() - t0
    no_solve_after_estimate = study.result is None

    # ---- then solve the same Study ------------------------------------------
    t0 = time.perf_counter()
    res = study.solve(comm)
    solve_seconds = time.perf_counter() - t0
    rss_gib = _summed_rss_gib(comm)  # collective: every rank
    solved_is_result = study.result is res
    sar = float(np.real(study.mean_sar()["mean_sar_w_per_kg"]))

    msh = res.mesh
    sigma_value = study.material.sigma
    eps_value = study.material.epsilon_r
    freq = study.frequency_hz

    # EX-57 setup figure; no-op unless FEM_EM_SETUP_FIGURES=1.
    write_setup_figure(
        msh,
        res.cell_tags,
        FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names={
            study.geometry.sphere_tag: f"gel sphere (sigma {sigma_value:g} S/m)",
            AIR_TAG: "air box (imposed-field walls)",
        },
        hide_tags=(AIR_TAG,),
        translucent_tags=(study.geometry.sphere_tag,),
        slice_normal=(0.0, 0.0, 1.0),
        slice_origin=(0.0, 0.0, 0.0),
        title="mri:12 -- Study.estimate() then the solve it priced (WF-13 step 2)",
        comm=comm,
    )

    # ---- (U): the identity ---------------------------------------------------
    V = res.fields.e_complex.function_space
    solve_unknowns = _global_dofs(V)
    solve_cells = int(msh.topology.index_map(msh.topology.dim).size_global)

    # ---- (In): the fixture is a priced row ----------------------------------
    p2 = _row("P2")

    # ---- (R), (G) -------------------------------------------------------------
    closed = _closed(sigma_value, eps_value)["sar_w_per_kg"]
    rel_record = abs(sar - WF15_DEG2_RECORD) / WF15_DEG2_RECORD
    rel_closed = abs(sar - closed) / closed

    # ---- (B-P2): WF-13's leave-one-out bracket on this fixture's row --------
    rm, rt = _ratios(PRICE_TABLE, "P2")

    # ---- (c1): a degree-1 space on the solved mesh (no form, no solve) ------
    V1 = fem.functionspace(msh, ("N1curl", 1))
    n1 = _global_dofs(V1)
    p1 = _row("P1")

    # ---- (c2): a P1-only table refuses degree 1 (every rank raises) ----------
    # Since 471a the full table prices degree 1 (P1, P3, T49 since 479), so the
    # refusal is exercised on a P1-only table through _predict(table=...).
    p1_only = tuple(r for r in PRICE_TABLE if r["row"] == "P1")
    c2_message = None
    try:
        _predict(n1, degree=1, ranks=comm.size, table=p1_only)
    except ValueError as exc:
        c2_message = str(exc)
    full_d1 = _predict(n1, degree=1, ranks=comm.size)

    # ---- (c3): the (R) comparator separates degrees (no solve) --------------
    c3 = abs(MAT4_DEG1_RECORD / WF15_DEG2_RECORD - 1.0)

    # ---- output fields (every collective on every rank) ---------------------
    e = res.fields.e_complex
    p1_space = fem.functionspace(msh, ("Lagrange", 1))
    e_abs = fem.Function(p1_space, name="E_abs")
    e_abs.interpolate(
        fem.Expression(ufl.sqrt(ufl.real(ufl.inner(e, e))), p1_space.element.interpolation_points)
    )
    e_abs.x.array[:] = np.real(e_abs.x.array)
    e_abs.x.scatter_forward()

    sigma_field = res.fields.sigma_field
    dg0 = fem.functionspace(msh, ("DG", 0))
    sar_field = fem.Function(dg0, name="SAR")
    sar_field.interpolate(
        fem.Expression(
            sigma_field * ufl.real(ufl.inner(e, e)) / (2.0 * study.density_kg_m3),
            dg0.element.interpolation_points,
        )
    )
    sar_field.x.array[:] = np.real(sar_field.x.array)
    sar_field.x.scatter_forward()

    rank_field = fem.Function(dg0, name="rank")
    n_owned = dg0.dofmap.index_map.size_local
    rank_field.x.array[:] = -1.0
    rank_field.x.array[:n_owned] = float(comm.rank)
    rank_field.x.scatter_forward()

    OUTPUT_DIR.mkdir(exist_ok=True)
    xdmf_path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_combined",
        msh,
        res.cell_tags,
        {"E_abs": e_abs, "SAR": sar_field, "rank": rank_field},
        comm=comm,
    )
    adopt_host_ownership(OUTPUT_DIR, comm=comm)

    # ---- prints (rank 0, no collective inside) -----------------------------
    if comm.rank == 0:
        print("=" * 72)
        print("EX-96 -- Study.estimate() on the lossy sphere, then the solve it priced")
        print("=" * 72)
        print(
            f"[study] sigma = {sigma_value} S/m, eps_r = {eps_value}, f = {freq:.6e} Hz; "
            f"Study.degree = {study.degree}, result.degree = {res.degree}"
        )
        print(
            f"[estimate] cells={est.cells} unknowns={est.unknowns} degree={est.degree} "
            f"ranks={est.ranks} pred_peak={est.predicted_peak_gib:.2f} GiB (summed over ranks) "
            f"pred_s={est.predicted_seconds:.1f} rows={est.rows_used}"
        )
        print(f"[estimate] the call took {estimate_seconds:.1f} s (mesh + space, no solve)")
        print(
            f"[solve] cells={solve_cells} unknowns={solve_unknowns}; Study solve+mesh "
            f"{solve_seconds:.1f} s; summed ru_maxrss {rss_gib:.2f} GiB"
        )
        print(
            f"[printed] measured/estimate: memory {rss_gib / est.predicted_peak_gib:.3f}, "
            f"time {solve_seconds / est.predicted_seconds:.3f} "
            f"(this machine; in-sample on row P2; not asserted)"
        )
        print(f"[(K)] result is None after estimate: {no_solve_after_estimate}; "
              f"study.result is res after solve: {solved_is_result}")
        print(
            f"[(U)] est.unknowns {est.unknowns} == solve {solve_unknowns} "
            f"(record {WF15_DEG2_UNKNOWNS_RECORD}); est.cells {est.cells} == solve "
            f"{solve_cells} (record {WF15_CELLS_RECORD})"
        )
        print(f"[(In)] 'P2' in rows_used: {'P2' in est.rows_used}; P2 unknowns {p2['unknowns']}")
        print(
            f"[(R)] SAR = {sar!r} ({sar:.10e}); record {WF15_DEG2_RECORD:.10e}; "
            f"rel {rel_record:.3e} (band L0_BAND {L0_BAND:g})"
        )
        print(
            f"[(G)] closed form {closed:.6e}; rel {rel_closed:.4%} (band MAT4_BAND {MAT4_BAND:g})"
        )
        print(
            f"[(B-P2)] leave-one-out P2: mem {rm:.3f}, time {rt:.3f} "
            f"(inside 1/{BRACKET:g}..{BRACKET:g}: {_inside(rm)}, {_inside(rt)}; "
            f"WF-13 record 0.635 / 0.795 (FH2, §9 item 585))"
        )
        print(
            f"[(c1)] degree-1 N1curl on the solved mesh: {n1} unknowns (P1 row {p1['unknowns']}); "
            f"deg2/deg1 ratio {est.unknowns / n1:.4f} (WF-13 record {WF13_DEG_RATIO_RECORD}); "
            f"an estimator that ignored degree would miss (U) by {est.unknowns - n1}"
        )
        print(f"[(c2)] P1-only table, degree-1 price raised ValueError: {c2_message!r}")
        print(
            f"[(c2-full)] full table, degree-1 priced since 471a (P1, P3, T49 since 479), "
            f"in-sample on its own rows: {full_d1[0]:.3f} GiB / {full_d1[1]:.1f} s "
            f"for {n1} unknowns (rows {full_d1[2]}); printed, not asserted"
        )
        print(
            f"[(c3)] |MAT4_DEG1_RECORD/WF15_DEG2_RECORD - 1| = {c3:.4e} "
            f"(must exceed L0_BAND {L0_BAND:g}; EX-91 read {EX91_C3_RECORD:.4e})"
        )
        print(f"[paraview] wrote {xdmf_path}")

    assert no_solve_after_estimate, "(K) estimate() solved"
    assert solved_is_result, "(K) study.result is not the solve's result"
    assert est.unknowns == solve_unknowns, f"(U) {est.unknowns} != {solve_unknowns}"
    assert est.cells == solve_cells, f"(U) {est.cells} != {solve_cells}"
    assert est.degree == study.degree == res.degree == 2, "(U) degree"
    assert est.ranks == comm.size == 2, f"(U) ranks {est.ranks} / {comm.size}"
    assert est.unknowns == WF15_DEG2_UNKNOWNS_RECORD, f"(U) record {est.unknowns}"
    assert est.cells == WF15_CELLS_RECORD, f"(U) record {est.cells}"
    assert "P2" in est.rows_used, f"(In) rows {est.rows_used}"
    assert p2["unknowns"] == est.unknowns, "(In) P2 unknowns"
    assert rel_record <= L0_BAND, f"(R) rel {rel_record:.3e} > {L0_BAND:g}"
    assert rel_closed <= MAT4_BAND, f"(G) rel {rel_closed:.4%} > {MAT4_BAND:g}"
    assert _inside(rm) and _inside(rt), f"(B-P2) {rm:.3f} / {rt:.3f}"
    assert n1 != est.unknowns, "(c1) degree-1 count equals the degree-2 count"
    assert n1 == p1["unknowns"], f"(c1) n1 {n1} != P1 {p1['unknowns']}"
    assert c2_message is not None, "(c2) P1-only degree-1 price did not raise"
    assert "fewer than two priced rows at degree 1" in c2_message, f"(c2) {c2_message!r}"
    assert c3 > L0_BAND, f"(c3) the (R) comparator does not separate degrees: {c3:.3e}"

    if comm.rank == 0:
        print(
            f"\nAll anchors and controls hold. Total elapsed "
            f"{time.perf_counter() - started:.1f} s.",
            flush=True,
        )


if __name__ == "__main__":
    main()
