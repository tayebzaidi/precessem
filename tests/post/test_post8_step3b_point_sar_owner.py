"""`POST-8` step 3b-i — ``point_sar`` defaults to the owner-cell rule.

PROJECT_PLAN §9 item 526; brief ``docs/planning/chunks/POST-8.md``
(``### §9 item 526`` and its 2026-10-18 03:00 review ruling).

Fixture: step 2's unit cube (``create_unit_cube(comm, 4, 4, 4)``) and its
three points (interior / facet / vertex).  Two non-polynomial real fields are
interpolated into N1curl degree 1 and then into vector DG1, exactly as
``core/time_harmonic.py`` splits the phasor, so the samples at the facet and
the vertex depend on the cell chosen.

Asserted: (D0) the default is the rule and ``owner=None`` is the legacy path;
(K) the known answer σ|E_c|²/(2ρ) in the oracle's minimum cell, by a route
independent of ``point_owner_cells``; (P) ``phantom_fields.py:163``
(``_evaluate_on_cells`` on owned-cell centroids) is owner-compliant; (W)
width invariance against the ``-n 2`` npz (complex build only); (N1) the
legacy cell differs from the owner at the widths step 2 measured.
"""

from pathlib import Path

import numpy as np
import pytest
from mpi4py import MPI

from dolfinx import fem

from fem_em_solver.post.evaluation import (
    evaluate_vector_field_parallel,
    point_owner_cells,
)
from fem_em_solver.post.phantom_fields import _cell_centroids, _evaluate_on_cells
from fem_em_solver.post.sar import point_sar
from tests.complex_mode import complex_only
from tests.post.test_point_owner_rule import (
    CLOSED_FORM_TOL,
    LABELS,
    POINTS,
    _bitwise_equal_across_ranks,
    _mesh,
    _oracle,
    _original,
)

SIGMA = 0.5
RHO = 1000.0
NPZ_DIR = Path("/workspace/logs/post8")
NPZ = NPZ_DIR / "step3b-n2.npz"  # never step3a-n2.npz
#: (N1) — step 2's measured legacy-vs-owner separations (W2 :165–166, W3 :79).
N1_CASES = {4: (1, 2), 1: (2,)}


def _f_r(x):
    return np.vstack((np.sin(3.0 * x[1]) + x[2], x[0] * np.cos(2.0 * x[2]), x[0] * x[1]))


def _f_i(x):
    return np.vstack((x[1] ** 2, np.sin(2.0 * x[0]), np.cos(x[1]) + x[2]))


def _fixture(comm):
    msh = _mesh(comm)
    v_n1 = fem.functionspace(msh, ("N1curl", 1))
    dg = fem.functionspace(msh, ("DG", 1, (3,)))
    out = []
    for f in (_f_r, _f_i):
        n1 = fem.Function(v_n1)
        n1.interpolate(f)
        d = fem.Function(dg)
        d.interpolate(n1)
        d.x.scatter_forward()
        out.append(d)
    return msh, out[0], out[1]


def _dg0_index(msh):
    v_dg = fem.functionspace(msh, ("DG", 0))
    f_dg = fem.Function(v_dg)
    n_owned = msh.topology.index_map(3).size_local
    f_dg.x.array[:n_owned] = _original(msh)[:n_owned]
    f_dg.x.scatter_forward()
    return f_dg


def _sar_in_cell(msh, e_r, e_i, point, original_index, comm):
    """σ|E|²/(2ρ) evaluated in the cell with the given original index.

    The owning rank calls ``.eval`` in that cell; the 6-component result is
    shared by buffer ``Allreduce`` SUM (independent of ``point_owner_cells``).
    """
    n_owned = msh.topology.index_map(3).size_local
    orig = _original(msh)[:n_owned]
    hit = np.flatnonzero(orig == int(original_index))
    buf = np.zeros(6, dtype=np.float64)
    if hit.size:
        cell = np.array([int(hit[0])], dtype=np.int32)
        p = np.asarray(point, dtype=np.float64).reshape(1, 3)
        buf[:3] = np.real(np.asarray(e_r.eval(p, cell)).reshape(-1))
        buf[3:] = np.real(np.asarray(e_i.eval(p, cell)).reshape(-1))
    tot = np.zeros_like(buf)
    comm.Allreduce(buf, tot, op=MPI.SUM)
    e = tot[:3] + 1j * tot[3:]
    return SIGMA * float(np.sum(np.abs(e) ** 2)) / (2.0 * RHO)


def _hand_sar(e_r, e_i, comm, owner):
    rv, vr = evaluate_vector_field_parallel(e_r, POINTS, comm, owner=owner)
    iv, vi = evaluate_vector_field_parallel(e_i, POINTS, comm, owner=owner)
    assert bool(np.all(vr & vi))
    e = np.real(rv) + 1j * np.real(iv)
    return float(SIGMA) * np.sum(np.abs(e) ** 2, axis=1) / (2.0 * float(RHO))


def _readings(comm):
    msh, e_r, e_i = _fixture(comm)
    sets = _oracle(msh, comm)
    owner_original, _ = point_owner_cells(msh, POINTS, comm)
    sar = point_sar(e_r, e_i, POINTS, sigma=SIGMA, rho=RHO)
    return msh, e_r, e_i, sets, owner_original, sar


def test_d0_k_default_is_rule_and_known_answer():
    comm = MPI.COMM_WORLD
    msh, e_r, e_i, sets, owner_original, sar_default = _readings(comm)
    sar_rule = point_sar(
        e_r, e_i, POINTS, sigma=SIGMA, rho=RHO, owner="lowest_original_cell"
    )
    sar_legacy = point_sar(e_r, e_i, POINTS, sigma=SIGMA, rho=RHO, owner=None)
    hand_legacy = _hand_sar(e_r, e_i, comm, None)

    # (K) known answer in the oracle's minimum cell
    sar_c = np.array(
        [_sar_in_cell(msh, e_r, e_i, POINTS[i], sets[i][0], comm) for i in range(3)]
    )
    rel_k = np.abs(sar_default - sar_c) / np.abs(sar_c)

    # incident-set spread at the vertex (printed, predicted > 1.01)
    vertex_vals = np.array(
        [_sar_in_cell(msh, e_r, e_i, POINTS[2], c, comm) for c in sets[2]]
    )

    # (N1) legacy cell via the DG0 index field
    f_dg = _dg0_index(msh)
    leg_dg, leg_ok = evaluate_vector_field_parallel(f_dg, POINTS, comm, owner=None)
    leg_dg = np.asarray(leg_dg).reshape(-1)
    assert bool(np.all(leg_ok))
    assert np.all(np.imag(leg_dg) == 0.0)
    legacy_cell = np.real(leg_dg).astype(np.int64)

    if comm.rank == 0:
        print(f"\n[POST-8 3b-i D0/K] -n {comm.size}")
        for i, lab in enumerate(LABELS):
            print(
                f"  {lab:8s} SAR(rule) {sar_default[i]:.17e}  SAR(legacy) "
                f"{sar_legacy[i]:.17e}  closed form {sar_c[i]:.17e}  rel {rel_k[i]:.3e}; "
                f"owner {int(owner_original[i])} (oracle min {sets[i][0]}), "
                f"legacy cell {int(legacy_cell[i])}"
            )
        print(
            f"  vertex incident set size {len(sets[2])}: SAR max/min "
            f"{vertex_vals.max() / vertex_vals.min():.6f} (predicted > 1.01, not asserted)"
        )
        rel_n2 = abs(sar_legacy[2] - sar_default[2]) / sar_default[2]
        print(
            f"  (N2) vertex |legacy - rule| / rule = {rel_n2:.3e} "
            f"(predicted >= 1e-2 where (N1) fires, printed only)"
        )

    # (D0)
    assert sar_default.tobytes() == sar_rule.tobytes()
    assert sar_legacy.tobytes() == hand_legacy.tobytes()
    for arr in (sar_default, sar_legacy, owner_original):
        assert _bitwise_equal_across_ranks(comm, np.asarray(arr))

    # (K)
    for i in range(3):
        assert int(owner_original[i]) == sets[i][0]
    assert float(np.max(rel_k)) <= CLOSED_FORM_TOL, rel_k

    # (N1) asserted only at the widths step 2 measured a separation
    for i in N1_CASES.get(comm.size, ()):
        assert int(legacy_cell[i]) != int(owner_original[i]), (
            LABELS[i], int(legacy_cell[i]), int(owner_original[i])
        )


def test_p_phantom_fields_centroids_owner_compliant():
    comm = MPI.COMM_WORLD
    msh, e_r, _ = _fixture(comm)
    n_owned = msh.topology.index_map(3).size_local
    cells = np.arange(n_owned, dtype=np.int32)
    local_centroids = _cell_centroids(msh, cells)
    local_orig = _original(msh)[:n_owned]

    counts = comm.allgather(int(n_owned))
    centroids = np.concatenate(comm.allgather(local_centroids))
    orig_all = np.concatenate(comm.allgather(local_orig))
    offset = int(sum(counts[: comm.rank]))

    owner_original, _ = point_owner_cells(msh, centroids, comm)
    rule_vals, rule_ok = evaluate_vector_field_parallel(
        e_r, centroids, comm, owner="lowest_original_cell"
    )
    pf_vals, _, _, n_fallback = _evaluate_on_cells(e_r, local_centroids, cells)

    rows = np.asarray(rule_vals)[offset : offset + n_owned]
    pf_vals = np.asarray(pf_vals).reshape(rows.shape)
    local_eq = bool(np.ascontiguousarray(rows).tobytes() == np.ascontiguousarray(pf_vals).astype(rows.dtype).tobytes())
    all_eq = comm.allreduce(int(local_eq), op=MPI.MIN)
    n_owner_ok = int(np.sum(owner_original == orig_all))

    if comm.rank == 0:
        print(
            f"\n[POST-8 3b-i P] -n {comm.size}: {centroids.shape[0]} centroids; "
            f"owner == own original index on {n_owner_ok}; "
            f"_evaluate_on_cells == rule bitwise on every rank: {bool(all_eq)}; "
            f"rank-0 fallback count {n_fallback}"
        )

    assert centroids.shape[0] == 384
    assert bool(np.all(rule_ok))
    assert n_owner_ok == 384
    assert all_eq == 1


@complex_only
def test_w_width_invariance():
    from tests.validation.test_post8_step3a_averaged_row1 import WIDTH_RTOL

    comm = MPI.COMM_WORLD
    _, _, _, _, owner_original, sar = _readings(comm)
    if comm.size == 2:
        if comm.rank == 0:
            NPZ_DIR.mkdir(parents=True, exist_ok=True)
            np.savez(NPZ, sar=sar, owner=np.asarray(owner_original, dtype=np.int64))
            print(f"\n[POST-8 3b-i W] -n 2 wrote {NPZ}: SAR {[f'{v:.17e}' for v in sar]}")
        comm.Barrier()
        return
    comm.Barrier()
    if not NPZ.exists():
        pytest.fail(f"{NPZ} missing — run the -n 2 window (W1) first; never recompute")
    ref = np.load(NPZ)
    rel = np.abs(sar - ref["sar"]) / np.abs(ref["sar"])
    if comm.rank == 0:
        print(
            f"\n[POST-8 3b-i W] -n {comm.size} vs -n 2: owner {owner_original.tolist()} "
            f"vs {ref['owner'].tolist()}; max rel SAR {float(rel.max()):.3e} "
            f"(asserted <= {WIDTH_RTOL:.0e}; predicted <= 1e-15)"
        )
    assert np.array_equal(np.asarray(owner_original, dtype=np.int64), ref["owner"])
    assert float(rel.max()) <= WIDTH_RTOL, rel
