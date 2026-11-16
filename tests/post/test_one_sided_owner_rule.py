"""`POST-1` step 7a (K) — the one-sided owner rule on a known answer.

PROJECT_PLAN §9 item 187; brief ``docs/planning/chunks/POST-1.md``.
``candidate_cells`` restricts ``point_owner_cells`` /
``evaluate_vector_field_parallel(owner="lowest_original_cell")`` to owned
cells ∩ the given local set; the tie-break is unchanged.

Fixture: ``create_unit_cube(comm, 4, 4, 4)`` and the face point
(0.5, 0.3, 0.6) on the hex-face plane x = 0.5.  Restricted to cells with
centroid x < 0.5 and then x > 0.5, each call must return the minimum of a DG0
``original_cell_index`` field over the incident set *on that side* (oracle:
allgathered vertices, barycentrics ≥ −1e-12, independent of DolfinX's
collision search), imaginary part 0.  Also: ``candidate_cells`` with
``owner=None`` raises ``ValueError``; an all-owned candidate set returns the
unrestricted values bitwise.  Runs in both builds.
"""

import numpy as np
import pytest
from mpi4py import MPI

from dolfinx import fem, mesh as dmesh

from fem_em_solver.post.evaluation import (
    evaluate_vector_field_parallel,
    point_owner_cells,
)

POINT = np.array([[0.5, 0.3, 0.6]], dtype=np.float64)
BARY_TOL = 1e-12


def _original(msh):
    oci = msh.topology.original_cell_index
    if isinstance(oci, (list, tuple)):
        oci = oci[0]
    return np.asarray(oci, dtype=np.int64)


def _cell_vertices(msh, cells):
    dofmap = msh.geometry.dofmap
    return np.stack([msh.geometry.x[dofmap[c]][:, :3] for c in cells])


def _oracle_side_sets(msh, comm):
    """Incident original indices at POINT, split by centroid side, rank 0 → all."""
    n_owned = msh.topology.index_map(3).size_local
    orig = _original(msh)[:n_owned]
    verts = _cell_vertices(msh, range(n_owned))
    all_orig = comm.gather(orig, root=0)
    all_verts = comm.gather(verts, root=0)
    out = None
    if comm.rank == 0:
        o = np.concatenate(all_orig)
        v = np.concatenate(all_verts)
        p = POINT[0]
        t = np.stack([v[:, 1] - v[:, 0], v[:, 2] - v[:, 0], v[:, 3] - v[:, 0]], axis=2)
        lam = np.linalg.solve(t, (p - v[:, 0])[:, :, None])[:, :, 0]
        bary = np.column_stack([1.0 - lam.sum(axis=1), lam])
        inside = np.all(bary >= -BARY_TOL, axis=1)
        cx = v[:, :, 0].mean(axis=1)
        out = {
            "left": sorted(int(x) for x in o[inside & (cx < 0.5)]),
            "right": sorted(int(x) for x in o[inside & (cx > 0.5)]),
        }
    return comm.bcast(out, root=0)


def _side_cells(msh, side):
    imap = msh.topology.index_map(3)
    n_local = imap.size_local + imap.num_ghosts
    cx = _cell_vertices(msh, range(n_local))[:, :, 0].mean(axis=1)
    sel = cx < 0.5 if side == "left" else cx > 0.5
    return np.flatnonzero(sel).astype(np.int32)


def _dg0_original(msh):
    v_dg = fem.functionspace(msh, ("DG", 0))
    f_dg = fem.Function(v_dg)
    n_owned = msh.topology.index_map(3).size_local
    f_dg.x.array[:n_owned] = _original(msh)[:n_owned]
    f_dg.x.scatter_forward()
    return f_dg


def test_k_one_sided_owner_returns_side_minimum():
    comm = MPI.COMM_WORLD
    msh = dmesh.create_unit_cube(comm, 4, 4, 4)
    sets = _oracle_side_sets(msh, comm)
    f_dg = _dg0_original(msh)

    results = {}
    for side in ("left", "right"):
        cells = _side_cells(msh, side)
        owner_original, owner_local = point_owner_cells(
            msh, POINT, comm, candidate_cells=cells
        )
        vals, valid = evaluate_vector_field_parallel(
            f_dg, POINT, comm, owner="lowest_original_cell", candidate_cells=cells
        )
        vals = np.asarray(vals).reshape(-1)
        claims = comm.allreduce((owner_local >= 0).astype(np.int64), op=MPI.SUM)
        results[side] = (owner_original, vals, valid, claims)

    unrestricted, _ = evaluate_vector_field_parallel(
        f_dg, POINT, comm, owner="lowest_original_cell"
    )
    unrestricted = np.asarray(unrestricted).reshape(-1)

    if comm.rank == 0:
        print(f"\n[POST-1 7a K] -n {comm.size}, point {POINT[0].tolist()}")
        for side in ("left", "right"):
            oo, vals, valid, claims = results[side]
            print(
                f"  {side:5s}: oracle set {sets[side]} (min "
                f"{sets[side][0] if sets[side] else None}); rule owner_original "
                f"{int(oo[0])}; DG0 value {vals[0]!r}; valid {bool(valid[0])}; "
                f"claims {int(claims[0])}"
            )
        print(f"  unrestricted DG0 {unrestricted[0]!r}")

    for side in ("left", "right"):
        oo, vals, valid, claims = results[side]
        assert len(sets[side]) >= 1, (side, sets)
        assert bool(valid[0])
        assert int(claims[0]) == 1
        assert int(oo[0]) == sets[side][0], (side, oo, sets[side])
        assert vals[0].real == float(sets[side][0]), (side, vals, sets[side])
        assert np.imag(vals[0]) == 0.0
    # The two sides are disjoint sets, and the unrestricted rule is their min.
    assert set(sets["left"]).isdisjoint(sets["right"])
    assert unrestricted[0].real == float(min(sets["left"][0], sets["right"][0]))


def test_k_candidate_cells_requires_owner_rule():
    comm = MPI.COMM_WORLD
    msh = dmesh.create_unit_cube(comm, 4, 4, 4)
    f_dg = _dg0_original(msh)
    with pytest.raises(ValueError, match="candidate_cells is valid only"):
        evaluate_vector_field_parallel(
            f_dg, POINT, comm, candidate_cells=np.arange(1, dtype=np.int32)
        )


def test_k_all_owned_candidates_equal_unrestricted_and_empty_side_invalid():
    comm = MPI.COMM_WORLD
    msh = dmesh.create_unit_cube(comm, 4, 4, 4)
    f_dg = _dg0_original(msh)
    n_owned = msh.topology.index_map(3).size_local
    pts = np.array([[0.13, 0.27, 0.61], [0.5, 0.3, 0.6], [0.5, 0.5, 0.5]])
    ref, ref_valid = evaluate_vector_field_parallel(
        f_dg, pts, comm, owner="lowest_original_cell"
    )
    got, got_valid = evaluate_vector_field_parallel(
        f_dg, pts, comm, owner="lowest_original_cell",
        candidate_cells=np.arange(n_owned, dtype=np.int32),
    )
    assert np.asarray(ref).tobytes() == np.asarray(got).tobytes()
    assert np.array_equal(ref_valid, got_valid)
    # A point with no candidate on its side: x > 0.5 cells only, point at x = 0.13.
    right = _side_cells(msh, "right")
    vals, valid = evaluate_vector_field_parallel(
        f_dg, pts[:1], comm, owner="lowest_original_cell", candidate_cells=right
    )
    assert not bool(valid[0])
