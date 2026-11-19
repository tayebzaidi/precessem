"""`POST-10` step 1 — the point-owner rule on a ``create_submesh`` mesh.

PROJECT_PLAN §9 item 254; brief ``docs/planning/chunks/POST-10.md``.  A
``create_submesh`` mesh carries no original cell index (DolfinX 0.11 raises
``Mixed topology unsupported.``), so ``owner="lowest_original_cell"`` takes
``parent=(parent_mesh, parent_cells)`` and keys each submesh cell by its
parent's original cell index.

Fixture: ``create_unit_cube(comm, 4, 4, 4)`` (384 cells, the `POST-8` K
fixture) cut at midpoint x < 0.75 (288 cells), and six points P0–P5
(interior, facet, vertex, boundary facet, boundary vertex, outside).  The
oracle is independent of DolfinX's collision search (``test_point_owner_rule``
idiom): owned parent cells' original index, vertices and keep flag are
gathered and rank 0 takes the kept incident set at barycentrics ≥ −1e-12.

Asserted: (S0) fixture identity; (K) the submesh owner equals the oracle's
minimum exactly, −1 outside, one claiming rank; (I) owner integers equal and
DG0 values / validity **bitwise** equal to the parent's ``candidate_cells=keep``
rule; (N) N1curl ``a + b×x`` to ≤ 1e-12; (W) bitwise equal on every rank (the
three windows run -n 1 / 2 / 4); (C1) the kept incident sets at P1 / P2 / P4
separate a wrong owner by ≥ 1; (C2) no ``parent`` raises ``RuntimeError`` on
every rank; (C3) ``parent`` without the rule raises ``ValueError``; (C4) a
short ``parent_cells`` raises ``ValueError`` on every rank.  Printed only:
the same-rank ownership fraction, (C5) the naive ``local_to_global`` key, and
the unrestricted parent owner at P3 / P4.  Runs in both builds.
"""

import numpy as np
import pytest
from dolfinx import fem
from dolfinx import mesh as dmesh
from mpi4py import MPI

from fem_em_solver.post.evaluation import (
    evaluate_vector_field_parallel,
    point_owner_cells,
)

POINTS = np.array(
    [
        [0.13, 0.27, 0.61],  # P0 interior
        [0.5, 0.3, 0.6],  # P1 facet plane x = 0.5
        [0.5, 0.5, 0.5],  # P2 interior vertex
        [0.75, 0.3, 0.6],  # P3 submesh-boundary facet
        [0.75, 0.5, 0.5],  # P4 submesh-boundary vertex
        [0.9, 0.3, 0.6],  # P5 outside the submesh
    ],
    dtype=np.float64,
)
LABELS = ("P0 interior", "P1 facet", "P2 vertex", "P3 bdry facet", "P4 bdry vertex", "P5 outside")
INSIDE = (0, 1, 2, 3, 4)
OUTSIDE = 5
X_CUT = 0.75
SUB_CELLS_GLOBAL = 288  # 6 * 4**3 * 3/4: three of four x-slabs
BARY_TOL = 1e-12
#: The guard module's bound (``test_evaluation_collective_guard.py:83``).
CLOSED_FORM_TOL = 1e-12
A_VEC = np.array([1.0, -2.0, 0.5])
B_VEC = np.array([0.3, 0.1, -0.2])
RULE = "lowest_original_cell"


def _original(msh):
    oci = msh.topology.original_cell_index
    if isinstance(oci, (list, tuple)):
        oci = oci[0]
    return np.asarray(oci, dtype=np.int64)


def _fixture(comm):
    msh = dmesh.create_unit_cube(comm, 4, 4, 4)
    tdim = msh.topology.dim
    imap = msh.topology.index_map(tdim)
    n_local = imap.size_local + imap.num_ghosts
    mid = dmesh.compute_midpoints(msh, tdim, np.arange(n_local, dtype=np.int32))
    keep = np.flatnonzero(mid[:, 0] < X_CUT).astype(np.int32)
    # The idiom, verbatim (test_mat10_sar_to_heat_sphere.py:82-88).
    sub, emap, _vmap, _gmap = dmesh.create_submesh(msh, tdim, keep)
    simap = sub.topology.index_map(tdim)
    n_sub = simap.size_local + simap.num_ghosts
    parent_ids = np.asarray(
        emap.sub_topology_to_topology(np.arange(n_sub, dtype=np.int32), False),
        dtype=np.int32,
    )
    return msh, keep, sub, parent_ids


def _oracle(msh, keep, comm):
    """Kept incident sets of parent original indices, from gathered vertices."""
    n_owned = msh.topology.index_map(3).size_local
    imap = msh.topology.index_map(3)
    keep_mask = np.zeros(imap.size_local + imap.num_ghosts, dtype=bool)
    keep_mask[keep] = True
    orig = _original(msh)[:n_owned]
    dofmap = msh.geometry.dofmap
    verts = np.stack([msh.geometry.x[dofmap[c]][:, :3] for c in range(n_owned)])
    all_orig = comm.gather(orig, root=0)
    all_verts = comm.gather(verts, root=0)
    all_keep = comm.gather(keep_mask[:n_owned], root=0)
    sets = None
    if comm.rank == 0:
        o = np.concatenate(all_orig)
        v = np.concatenate(all_verts)
        k = np.concatenate(all_keep)
        sets = []
        for p in POINTS:
            t = np.stack([v[:, 1] - v[:, 0], v[:, 2] - v[:, 0], v[:, 3] - v[:, 0]], axis=2)
            lam = np.linalg.solve(t, (p - v[:, 0])[:, :, None])[:, :, 0]
            bary = np.column_stack([1.0 - lam.sum(axis=1), lam])
            inside = np.all(bary >= -BARY_TOL, axis=1)
            sets.append(sorted(int(x) for x in o[inside & k]))
    return comm.bcast(sets, root=0)


def _bitwise_equal_across_ranks(comm, arr):
    b = np.ascontiguousarray(arr).tobytes()
    return all(x == b for x in comm.allgather(b))


def _dg0(msh, cell_values):
    """Scalar DG0 holding ``cell_values[c]`` on owned cells, then scatter."""
    v = fem.functionspace(msh, ("DG", 0))
    f = fem.Function(v)
    n_owned = msh.topology.index_map(3).size_local
    for c in range(n_owned):
        f.x.array[v.dofmap.cell_dofs(c)[0]] = cell_values[c]
    f.x.scatter_forward()
    return f


def test_s0_fixture_identity():
    comm = MPI.COMM_WORLD
    msh, keep, sub, parent_cells = _fixture(comm)
    simap = sub.topology.index_map(3)
    pimap = msh.topology.index_map(3)
    n_sub_owned = simap.size_local
    kept_orig = _original(msh)[parent_cells[:n_sub_owned]]
    gathered = np.concatenate(comm.allgather(kept_orig))
    same_rank = int(np.count_nonzero(parent_cells[:n_sub_owned] < pimap.size_local))
    tot_same = comm.allreduce(same_rank, op=MPI.SUM)
    if comm.rank == 0:
        print(
            f"\n[POST-10 S0] -n {comm.size}: sub size_global {simap.size_global}; "
            f"kept parent originals gathered {gathered.size}, unique "
            f"{np.unique(gathered).size}; owned sub cells whose parent is owned "
            f"on the same rank: {tot_same}/{simap.size_global} = "
            f"{tot_same / max(simap.size_global, 1):.4f} (predicted 1.0, not asserted)"
        )
    per_rank = comm.allgather((simap.size_local, simap.num_ghosts, len(parent_cells)))
    if comm.rank == 0:
        print(f"  per rank (sub size_local, num_ghosts, len(parent_cells)): {per_rank}")
    assert simap.size_global == SUB_CELLS_GLOBAL
    assert len(parent_cells) == simap.size_local + simap.num_ghosts
    assert gathered.size == SUB_CELLS_GLOBAL
    assert np.unique(gathered).size == SUB_CELLS_GLOBAL


def _naive_key_values(sub, f_sub_dg, comm):
    """(C5): key by the sub index map's partition-dependent global number."""
    from dolfinx import geometry

    simap = sub.topology.index_map(3)
    n_owned = simap.size_local
    l2g = simap.local_to_global(np.arange(n_owned, dtype=np.int32))
    tree = geometry.bb_tree(sub, 3)
    hits = geometry.compute_collisions_points(tree, POINTS)
    coll = geometry.compute_colliding_cells(sub, hits, POINTS)
    dofmap = f_sub_dg.function_space.dofmap
    local = []
    for i in range(POINTS.shape[0]):
        links = np.asarray(coll.links(i))
        owned = links[links < n_owned]
        if owned.size == 0:
            local.append(None)
            continue
        j = int(np.argmin(l2g[owned]))
        c = int(owned[j])
        local.append((int(l2g[c]), complex(f_sub_dg.x.array[dofmap.cell_dofs(c)[0]])))
    gathered = comm.allgather(local)
    out = []
    for i in range(POINTS.shape[0]):
        cands = [g[i] for g in gathered if g[i] is not None]
        out.append(min(cands)[1].real if cands else None)
    return out


def test_k_i_n_w_submesh_rule():
    comm = MPI.COMM_WORLD
    msh, keep, sub, parent_cells = _fixture(comm)
    sets = _oracle(msh, keep, comm)
    parent = (msh, parent_cells)
    orig = _original(msh)

    # (K) owner integers
    sub_oo, sub_local = point_owner_cells(sub, POINTS, comm, parent=parent)
    par_oo, par_local = point_owner_cells(msh, POINTS, comm, candidate_cells=keep)
    unr_oo, _ = point_owner_cells(msh, POINTS, comm)

    # planted DG0 fields
    f_par = _dg0(msh, orig)
    f_sub = _dg0(sub, orig[parent_cells])
    sub_vals, sub_valid = evaluate_vector_field_parallel(
        f_sub, POINTS, comm, owner=RULE, parent=parent
    )
    par_vals, par_valid = evaluate_vector_field_parallel(
        f_par, POINTS, comm, owner=RULE, candidate_cells=keep
    )
    sub_vals = np.asarray(sub_vals).reshape(-1)
    par_vals = np.asarray(par_vals).reshape(-1)

    # (N) N1curl a + b x x on the submesh
    v_n1 = fem.functionspace(sub, ("N1curl", 1))
    f_n1 = fem.Function(v_n1)
    f_n1.interpolate(lambda x: A_VEC[:, None] + np.cross(B_VEC, x.T).T)
    n1_vals, n1_valid = evaluate_vector_field_parallel(
        f_n1, POINTS, comm, owner=RULE, parent=parent
    )
    n1_vals = np.asarray(n1_vals)
    n1_exact = A_VEC[None, :] + np.cross(B_VEC[None, :], POINTS)
    err_n1 = float(np.max(np.abs(n1_vals[list(INSIDE)] - n1_exact[list(INSIDE)])))

    # (C5) printed only
    naive = _naive_key_values(sub, f_sub, comm)

    claims = comm.allreduce((sub_local >= 0).astype(np.int64), op=MPI.SUM)

    if comm.rank == 0:
        print(f"\n[POST-10 K/I/N] -n {comm.size}")
        for i, lab in enumerate(LABELS):
            s = sets[i]
            print(
                f"  {lab:15s} {POINTS[i].tolist()}: kept set size {len(s)}, "
                f"min {s[0] if s else None}, max {s[-1] if s else None}, "
                f"sep {(s[-1] - s[0]) if s else None}; submesh owner {int(sub_oo[i])}; "
                f"parent keep owner {int(par_oo[i])}; unrestricted parent owner "
                f"{int(unr_oo[i])}; sub DG0 {sub_vals[i].real:.1f} "
                f"(valid {bool(sub_valid[i])}); naive l2g-key DG0 {naive[i]}"
            )
        print(f"  owner integers -n {comm.size}: {[int(x) for x in sub_oo]}")
        print(f"  claims {claims.tolist()}; N1curl max |err| P0-P4 {err_n1:.3e}")
        diffs = [i for i in (1, 2, 4) if naive[i] is not None and naive[i] != sub_vals[i].real]
        print(
            f"  (C5) naive key differs from the rule at points {diffs} (predicted >= 1 at some width)"
        )

    # (C1) discrimination ceiling
    for i in (1, 2, 4):
        assert len(sets[i]) >= 2, (i, sets[i])
        assert sets[i][-1] - sets[i][0] >= 1, (i, sets[i])
    for i in INSIDE:
        assert len(sets[i]) >= 1 and len(set(sets[i])) == len(sets[i])
    assert sets[OUTSIDE] == []

    # (K)
    for i in INSIDE:
        assert int(sub_oo[i]) == sets[i][0], (i, sub_oo, sets[i])
        assert bool(sub_valid[i])
        assert np.imag(sub_vals[i]) == 0.0
        assert sub_vals[i].real == float(sets[i][0]), (i, sub_vals[i], sets[i])
        assert int(claims[i]) == 1, (i, claims)
    assert int(sub_oo[OUTSIDE]) == -1
    assert not bool(sub_valid[OUTSIDE])
    assert int(claims[OUTSIDE]) == 0

    # (I) the identity with the parent's one-sided rule, band 0
    assert np.array_equal(sub_oo, par_oo), (sub_oo, par_oo)
    assert np.array_equal(sub_vals, par_vals), (sub_vals, par_vals)
    assert np.array_equal(sub_valid, par_valid), (sub_valid, par_valid)

    # (N)
    assert all(bool(n1_valid[i]) for i in INSIDE)
    assert not bool(n1_valid[OUTSIDE])
    assert err_n1 <= CLOSED_FORM_TOL, err_n1

    # (W) bitwise equal on every rank
    for arr in (sub_oo, sub_vals, n1_vals):
        assert _bitwise_equal_across_ranks(comm, np.asarray(arr))


def test_c2_submesh_without_parent_raises_on_every_rank():
    comm = MPI.COMM_WORLD
    msh, keep, sub, parent_cells = _fixture(comm)
    f_sub = _dg0(sub, _original(msh)[parent_cells])
    raised = False
    message = ""
    try:
        evaluate_vector_field_parallel(f_sub, POINTS, comm, owner=RULE)
    except RuntimeError as exc:
        raised = True
        message = str(exc)
    flags = comm.allgather(raised)
    if comm.rank == 0:
        print(f"\n[POST-10 C2] -n {comm.size}: flags {flags}; message {message!r}")
    assert flags == [True] * comm.size, flags
    assert "parent=" in message and "create_submesh" in message


def test_c3_parent_requires_owner_rule():
    comm = MPI.COMM_WORLD
    msh, keep, sub, parent_cells = _fixture(comm)
    f_sub = _dg0(sub, _original(msh)[parent_cells])
    with pytest.raises(ValueError, match="parent is valid only"):
        evaluate_vector_field_parallel(f_sub, POINTS, comm, parent=(msh, parent_cells))


def test_c4_short_parent_cells_raises_on_every_rank():
    comm = MPI.COMM_WORLD
    msh, keep, sub, parent_cells = _fixture(comm)
    f_sub = _dg0(sub, _original(msh)[parent_cells])
    raised = False
    message = ""
    try:
        evaluate_vector_field_parallel(
            f_sub, POINTS, comm, owner=RULE, parent=(msh, parent_cells[:-1])
        )
    except ValueError as exc:
        raised = True
        message = str(exc)
    flags = comm.allgather(raised)
    if comm.rank == 0:
        print(f"\n[POST-10 C4] -n {comm.size}: flags {flags}; message {message!r}")
    assert flags == [True] * comm.size, flags
    assert "parent" in message
