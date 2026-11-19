"""`POST-8` step 2 — the deterministic point-owner rule on a known answer.

PROJECT_PLAN §9 item 177; brief ``docs/planning/chunks/POST-8.md``.  The rule
(09-30 weekly): raw point evaluation takes the cell with the **lowest original
cell index** among the cells containing the point
(``post.evaluation.point_owner_cells``; opt-in as
``evaluate_vector_field_parallel(..., owner="lowest_original_cell")``).

Fixture: ``create_unit_cube(comm, 4, 4, 4)`` and three points — interior
(one incident cell), on a hex-face plane (a shared facet, two cells) and an
interior vertex (many cells).  The oracle is independent of DolfinX's
collision search: every owned cell's original index and four vertices are
allgathered and rank 0 takes the incident set at barycentrics ≥ −1e-12.

Asserted: K0 (original indices unique, count = ``size_global``), K1 / K2
(closed forms of a CG1 and an N1curl field under the rule to ≤ 1e-12), K3 (the
tie-break returns the oracle's minimum, imaginary part 0), K4 (bitwise-equal
on every rank), K5 (a rank-mismatched ``owner`` raises on every rank; the
``-n 1`` control does not raise), K6 (an unknown rule raises ``ValueError``),
CK1 (the incident set's max ≠ min at the facet and vertex), CK2 (the legacy
vertex value is a member of the incident set).  Runs in both builds.
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
        [0.13, 0.27, 0.61],  # interior
        [0.5, 0.3, 0.6],  # hex-face plane x = 0.5 (a shared facet)
        [0.5, 0.5, 0.5],  # interior vertex
    ],
    dtype=np.float64,
)
LABELS = ("interior", "facet", "vertex")
BARY_TOL = 1e-12
#: The guard module's bound (``test_evaluation_collective_guard.py:83``).
CLOSED_FORM_TOL = 1e-12
A_VEC = np.array([1.0, -2.0, 0.5])
B_VEC = np.array([0.3, 0.1, -0.2])


def _mesh(comm):
    return dmesh.create_unit_cube(comm, 4, 4, 4)


def _original(msh):
    oci = msh.topology.original_cell_index
    if isinstance(oci, (list, tuple)):
        oci = oci[0]
    return np.asarray(oci, dtype=np.int64)


def _oracle(msh, comm):
    """Incident sets of original indices per point, from allgathered vertices."""
    n_owned = msh.topology.index_map(3).size_local
    orig = _original(msh)[:n_owned]
    dofmap = msh.geometry.dofmap
    verts = np.stack([msh.geometry.x[dofmap[c]][:, :3] for c in range(n_owned)])
    all_orig = comm.gather(orig, root=0)
    all_verts = comm.gather(verts, root=0)
    sets = None
    if comm.rank == 0:
        o = np.concatenate(all_orig)
        v = np.concatenate(all_verts)
        sets = []
        for p in POINTS:
            t = np.stack([v[:, 1] - v[:, 0], v[:, 2] - v[:, 0], v[:, 3] - v[:, 0]], axis=2)
            lam = np.linalg.solve(t, (p - v[:, 0])[:, :, None])[:, :, 0]
            bary = np.column_stack([1.0 - lam.sum(axis=1), lam])
            inside = np.all(bary >= -BARY_TOL, axis=1)
            sets.append(sorted(int(x) for x in o[inside]))
    return comm.bcast(sets, root=0)


def _bitwise_equal_across_ranks(comm, arr):
    b = np.ascontiguousarray(arr).tobytes()
    return all(x == b for x in comm.allgather(b))


def test_k0_original_cell_index_unique_and_complete():
    comm = MPI.COMM_WORLD
    msh = _mesh(comm)
    raw = msh.topology.original_cell_index
    imap = msh.topology.index_map(3)
    if comm.rank == 0:
        print(
            f"\n[POST-8 K0] -n {comm.size}: original_cell_index type {type(raw).__name__}, "
            f"dtype {getattr(raw, 'dtype', None)}, shape {np.shape(raw)}; "
            f"size_local {imap.size_local}, num_ghosts {imap.num_ghosts}, "
            f"size_global {imap.size_global}"
        )
    owned = _original(msh)[: imap.size_local]
    gathered = np.concatenate(comm.allgather(owned))
    assert gathered.size == imap.size_global
    assert np.unique(gathered).size == imap.size_global


def test_k1_k2_k3_k4_closed_forms_and_tie_break():
    comm = MPI.COMM_WORLD
    msh = _mesh(comm)
    sets = _oracle(msh, comm)
    owner_original, owner_local = point_owner_cells(msh, POINTS, comm)

    # K1: CG1 (x, 2y, 3z)
    v_cg = fem.functionspace(msh, ("Lagrange", 1, (3,)))
    f_cg = fem.Function(v_cg)
    f_cg.interpolate(lambda x: np.vstack((x[0], 2.0 * x[1], 3.0 * x[2])))
    cg_vals, cg_valid = evaluate_vector_field_parallel(
        f_cg, POINTS, comm, owner="lowest_original_cell"
    )
    cg_exact = POINTS * np.array([1.0, 2.0, 3.0])
    err_cg = float(np.max(np.abs(np.asarray(cg_vals) - cg_exact)))

    # K2: N1curl a + b x x
    v_n1 = fem.functionspace(msh, ("N1curl", 1))
    f_n1 = fem.Function(v_n1)
    f_n1.interpolate(lambda x: A_VEC[:, None] + np.cross(B_VEC, x.T).T)
    n1_vals, n1_valid = evaluate_vector_field_parallel(
        f_n1, POINTS, comm, owner="lowest_original_cell"
    )
    n1_exact = A_VEC[None, :] + np.cross(B_VEC[None, :], POINTS)
    err_n1 = float(np.max(np.abs(np.asarray(n1_vals) - n1_exact)))

    # K3: DG0 holding original_cell_index
    v_dg = fem.functionspace(msh, ("DG", 0))
    f_dg = fem.Function(v_dg)
    imap = msh.topology.index_map(3)
    n_owned = imap.size_local
    orig = _original(msh)
    f_dg.x.array[:n_owned] = orig[:n_owned]
    f_dg.x.scatter_forward()
    dg_vals, dg_valid = evaluate_vector_field_parallel(
        f_dg, POINTS, comm, owner="lowest_original_cell"
    )
    dg_vals = np.asarray(dg_vals).reshape(-1)

    # CK2: the legacy path at the vertex
    legacy_vals, _ = evaluate_vector_field_parallel(f_dg, POINTS, comm)
    legacy_vals = np.asarray(legacy_vals).reshape(-1)

    if comm.rank == 0:
        print(f"\n[POST-8 K] -n {comm.size}")
        for i, lab in enumerate(LABELS):
            s = sets[i]
            print(
                f"  {lab:8s} {POINTS[i].tolist()}: oracle incident set size {len(s)}, "
                f"min {s[0] if s else None}, max {s[-1] if s else None}; "
                f"rule owner_original {int(owner_original[i])}; "
                f"rule DG0 {dg_vals[i].real:.1f}; legacy DG0 {legacy_vals[i].real:.1f}"
                f" (legacy == min: {bool(s) and legacy_vals[i].real == s[0]})"
            )
        print(f"  K1 CG1 max |err| {err_cg:.3e}; K2 N1curl max |err| {err_n1:.3e}")
        print(f"  oracle vertex set (predicted 24, not asserted): {sets[2]}")

    assert bool(np.all(cg_valid)) and bool(np.all(n1_valid)) and bool(np.all(dg_valid))
    assert err_cg <= CLOSED_FORM_TOL, err_cg
    assert err_n1 <= CLOSED_FORM_TOL, err_n1

    # K3
    assert len(sets[0]) == 1
    assert len(sets[1]) >= 2 and len(sets[2]) >= 2
    for i in range(3):
        assert len(set(sets[i])) == len(sets[i])
        assert int(owner_original[i]) == sets[i][0], (i, owner_original[i], sets[i])
        assert dg_vals[i].real == float(sets[i][0]), (i, dg_vals[i], sets[i])
        assert np.imag(dg_vals[i]) == 0.0
    # owner_local: exactly one rank claims each point
    claims = comm.allreduce((owner_local >= 0).astype(np.int64), op=MPI.SUM)
    assert np.all(claims == 1), claims
    # CK1: by construction once K3 holds with distinct indices
    for i in (1, 2):
        assert max(sets[i]) != min(sets[i])
    # CK2: the legacy value is a member of the incident set
    assert int(round(legacy_vals[2].real)) in sets[2]

    # K4: bitwise-equal on every rank
    for arr in (cg_vals, n1_vals, dg_vals, owner_original):
        assert _bitwise_equal_across_ranks(comm, np.asarray(arr))


def test_k5_mismatched_owner_raises_on_every_rank():
    comm = MPI.COMM_WORLD
    msh = _mesh(comm)
    v_cg = fem.functionspace(msh, ("Lagrange", 1, (3,)))
    f = fem.Function(v_cg)
    owner = None if comm.rank == 0 else "lowest_original_cell"
    raised = False
    message = ""
    try:
        evaluate_vector_field_parallel(f, POINTS, comm, owner=owner)
    except ValueError as exc:
        raised = True
        message = str(exc)
    flags = comm.allgather(raised)
    if comm.rank == 0:
        print(f"\n[POST-8 K5] -n {comm.size}: flags {flags}; message {message!r}")
    if comm.size == 1:
        assert flags == [False]
        return
    assert flags == [True] * comm.size, flags
    assert "owner" in message and "collective" in message


def test_k6_unknown_owner_raises():
    comm = MPI.COMM_WORLD
    msh = _mesh(comm)
    v_cg = fem.functionspace(msh, ("Lagrange", 1, (3,)))
    f = fem.Function(v_cg)
    with pytest.raises(ValueError, match="unknown owner rule"):
        evaluate_vector_field_parallel(f, POINTS, comm, owner="first_cell")
