"""Utilities for robust field evaluation on distributed meshes."""

from __future__ import annotations

import hashlib

import numpy as np
from mpi4py import MPI
from dolfinx import geometry

#: The accepted values of ``evaluate_vector_field_parallel``'s ``owner``
#: keyword.  ``None`` is the legacy path (first colliding cell on every rank
#: that finds one, last gathered rank wins), byte for byte.
OWNER_RULES = (None, "lowest_original_cell")

_NO_OWNER = np.iinfo(np.int64).max


def _original_cell_index(mesh) -> np.ndarray:
    """``mesh.topology.original_cell_index`` as one int64 array.

    DolfinX 0.11 returns an ndarray; a per-cell-type list is accepted too
    (single cell type ⇒ entry 0).
    """
    oci = mesh.topology.original_cell_index
    if isinstance(oci, (list, tuple)):
        oci = oci[0]
    return np.asarray(oci, dtype=np.int64)


_KEYS_OK = 0
_KEYS_NO_INDEX = 1
_KEYS_BAD_PARENT = 2


def _owner_keys(mesh, parent, comm) -> np.ndarray:
    """The owner rule's per-cell keys, collectively checked (`POST-10` step 1).

    ``parent=None``: the mesh's own original cell index (unchanged).
    ``parent=(parent_mesh, parent_cells)``: ``parent_cells`` is the
    parent-local cell index of every local (owned + ghost) cell of ``mesh``
    (``EntityMap.sub_topology_to_topology(arange(n_sub), False)`` for a
    ``create_submesh`` mesh); the key is the parent's original cell index of
    that cell.

    Every rank computes its keys or a local error code, then **one**
    ``allreduce(MAX)`` of the code runs on both paths, and a failure on any
    rank raises on every rank — a rank-local raise ahead of the owner rule's
    collectives would hang the others.
    """
    status = _KEYS_OK
    keys = None
    detail = ""
    cause = None
    if parent is None:
        try:
            keys = _original_cell_index(mesh)
        except RuntimeError as exc:  # DolfinX 0.11 on a create_submesh mesh
            status = _KEYS_NO_INDEX
            cause = exc
            detail = str(exc)
    else:
        try:
            parent_mesh, parent_cells = parent
        except (TypeError, ValueError):
            status = _KEYS_BAD_PARENT
            detail = "parent must be a (parent_mesh, parent_cells) pair"
        if status == _KEYS_OK:
            tdim = mesh.topology.dim
            imap = mesh.topology.index_map(tdim)
            n_local = imap.size_local + imap.num_ghosts
            cells = np.asarray(parent_cells, dtype=np.int64).reshape(-1)
            pimap = parent_mesh.topology.index_map(parent_mesh.topology.dim)
            n_parent = pimap.size_local + pimap.num_ghosts
            if cells.size != n_local:
                status = _KEYS_BAD_PARENT
                detail = (
                    f"parent_cells has length {cells.size}; the mesh has "
                    f"{n_local} local (owned + ghost) cells"
                )
            elif cells.size and (cells.min() < 0 or cells.max() >= n_parent):
                status = _KEYS_BAD_PARENT
                detail = (
                    f"parent_cells must be parent-local cell indices in "
                    f"[0, {n_parent}); got range [{cells.min()}, {cells.max()}]"
                )
            else:
                try:
                    keys = _original_cell_index(parent_mesh)[cells]
                except RuntimeError as exc:  # the parent is itself a submesh
                    status = _KEYS_NO_INDEX
                    cause = exc
                    detail = str(exc)

    worst = comm.allreduce(status, op=MPI.MAX)
    if worst == _KEYS_OK:
        return keys
    ranks = [r for r, s in enumerate(comm.allgather(status)) if s != _KEYS_OK]
    if worst == _KEYS_NO_INDEX:
        raise RuntimeError(
            "point_owner_cells: the mesh carries no original cell index "
            f"(failing ranks {ranks}; DolfinX: {detail or 'see another rank'}). "
            "A create_submesh mesh has none — pass "
            "parent=(parent_mesh, parent_cells) with parent_cells the "
            "parent-local index of every local sub cell (POST-10)."
        ) from cause
    raise ValueError(
        f"point_owner_cells: malformed parent= on ranks {ranks}: "
        f"{detail or 'see another rank'}"
    )


def _candidate_mask(mesh, candidate_cells) -> np.ndarray | None:
    """Boolean mask over local cells (owned + ghost) for ``candidate_cells``.

    ``None`` -> ``None`` (no restriction).  `POST-1` step 7a.
    """
    if candidate_cells is None:
        return None
    imap = mesh.topology.index_map(mesh.topology.dim)
    n_local = imap.size_local + imap.num_ghosts
    cells = np.asarray(candidate_cells, dtype=np.int64).reshape(-1)
    if cells.size and (cells.min() < 0 or cells.max() >= n_local):
        raise ValueError(
            "candidate_cells must be local cell indices in "
            f"[0, {n_local}); got range [{cells.min()}, {cells.max()}]"
        )
    mask = np.zeros(n_local, dtype=bool)
    mask[cells] = True
    return mask


def point_owner_cells(
    mesh,
    points: np.ndarray,
    comm: MPI.Intracomm | None = None,
    *,
    candidate_cells=None,
    parent=None,
):
    """Deterministic, partition-free owner cell of each point (`POST-8` step 2).

    The rule (09-30 weekly): among all cells containing a point, the owner is
    the one with the **lowest original cell index** — the cell's position in
    the array given to ``create_mesh``, which is set by the mesher and not by
    the partition.  The index-map global number (``local_to_global``) is
    partition-dependent and is deliberately not used (step 1 measured the same
    cell as global 85735 at ``-n 2`` and 68816 at ``-n 4``).

    Collective.  Only *owned* cells are candidates (ghosts are owned elsewhere
    and would make two ranks claim one point).

    ``candidate_cells`` (keyword only, `POST-1` step 7a): local cell indices
    restricting the candidates to *owned cells ∩ candidate_cells* — the
    one-sided rule at a material interface (pass ``cell_tags.find(tag)``).
    The tie-break is unchanged.  A point with no candidate on any rank gets
    ``-1``.  ``None`` (default) is the unrestricted rule, byte for byte.

    ``parent`` (keyword only, `POST-10` step 1): ``None`` or
    ``(parent_mesh, parent_cells)``.  A ``create_submesh`` mesh carries no
    original cell index (DolfinX 0.11 raises ``Mixed topology unsupported``),
    so the caller passes the parent map: ``parent_cells`` holds the
    parent-local cell index of every local (owned + ghost) cell of ``mesh``,
    and each cell is keyed by its parent's original cell index.  The rule is
    otherwise unchanged.  A mesh without an index and no ``parent`` raises
    ``RuntimeError`` on every rank; malformed ``parent_cells`` raises
    ``ValueError`` on every rank (one collectively checked status code).

    Returns
    -------
    owner_original : np.ndarray[int64], shape (N,)
        The owner's original cell index, identical on every rank; ``-1`` for a
        point no cell contains.
    owner_local_cell : np.ndarray[int32], shape (N,)
        The owner's local cell index on the owning rank; ``-1`` on every other
        rank and for unowned points.
    """
    comm = comm or mesh.comm
    points = np.ascontiguousarray(np.asarray(points, dtype=np.float64))
    if points.ndim != 2 or points.shape[1] != 3:
        raise ValueError(f"points must have shape (N, 3), got {points.shape}")
    n_points = points.shape[0]
    tdim = mesh.topology.dim
    n_owned = mesh.topology.index_map(tdim).size_local
    original = _owner_keys(mesh, parent, comm)
    mask = _candidate_mask(mesh, candidate_cells)

    bb_tree = geometry.bb_tree(mesh, tdim)
    bb_hits = geometry.compute_collisions_points(bb_tree, points)
    colliding_cells = geometry.compute_colliding_cells(mesh, bb_hits, points)

    local_min = np.full(n_points, _NO_OWNER, dtype=np.int64)
    local_cell = np.full(n_points, -1, dtype=np.int32)
    for i in range(n_points):
        links = np.asarray(colliding_cells.links(i))
        owned = links[links < n_owned]
        if mask is not None:
            owned = owned[mask[owned]]
        if owned.size == 0:
            continue
        orig = original[owned]
        j = int(np.argmin(orig))
        local_min[i] = orig[j]
        local_cell[i] = int(owned[j])

    global_min = np.empty_like(local_min)
    comm.Allreduce(local_min, global_min, op=MPI.MIN)

    claims = ((local_min == global_min) & (global_min != _NO_OWNER)).astype(np.int64)
    claim_count = np.empty_like(claims)
    comm.Allreduce(claims, claim_count, op=MPI.SUM)
    if np.any(claim_count > 1):
        bad = [int(i) for i in np.flatnonzero(claim_count > 1)]
        raise RuntimeError(
            "point_owner_cells: more than one rank claims the owner of point(s) "
            f"{bad[:10]} — an original cell index is owned on two ranks"
        )

    owner_local_cell = np.where(claims == 1, local_cell, -1).astype(np.int32)
    owner_original = np.where(global_min == _NO_OWNER, -1, global_min).astype(np.int64)
    return owner_original, owner_local_cell


def evaluate_vector_field_parallel(
    function,
    points: np.ndarray,
    comm: MPI.Intracomm | None = None,
    *,
    owner: str | None = None,
    candidate_cells=None,
    parent=None,
):
    """Evaluate a vector-valued dolfinx function at points across MPI ranks.

    Parameters
    ----------
    function:
        dolfinx function to evaluate.
    points:
        Array with shape (N, 3).
    comm:
        MPI communicator. Defaults to ``function.function_space.mesh.comm``.
    owner:
        Keyword only.  ``None`` (default) is the legacy path, unchanged.
        ``"lowest_original_cell"`` evaluates each point in the cell chosen by
        :func:`point_owner_cells` — the same cell at every MPI width, so a
        field discontinuous across cell faces (N1curl) returns a
        width-invariant value (`POST-8` step 2).  Must be identical on every
        rank (validated with the OPS-40/55 key).
    candidate_cells:
        Keyword only (`POST-1` step 7a).  Local cell indices; the owner is
        chosen among *owned cells ∩ candidate_cells* (one-sided evaluation at
        a material interface).  Valid only with
        ``owner="lowest_original_cell"`` — anything else raises
        ``ValueError`` before any collective.  Whether it is ``None`` must
        agree on every rank (it rides the OPS-40/55 key).  A point with no
        candidate returns ``valid=False``.  ``None`` is byte-identical to the
        call without it.
    parent:
        Keyword only (`POST-10` step 1).  ``None`` or
        ``(parent_mesh, parent_cells)`` for a ``create_submesh`` mesh:
        ``parent_cells`` is the parent-local cell index of every local
        (owned + ghost) cell of the function's mesh, and the owner rule keys
        each cell by its parent's original cell index (see
        :func:`point_owner_cells`).  Valid only with
        ``owner="lowest_original_cell"`` — anything else raises
        ``ValueError`` before any collective.  Whether it is ``None`` must
        agree on every rank (it rides the OPS-40/55 key).

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        values:
            Array of shape (N, value_size) with point evaluations.
            Invalid/outside points are filled with zeros.
        valid_mask:
            Boolean mask with shape (N,) marking which points were evaluated.
    """
    mesh = function.function_space.mesh
    comm = comm or mesh.comm

    points = np.asarray(points, dtype=np.float64)
    if points.ndim != 2 or points.shape[1] != 3:
        raise ValueError(f"points must have shape (N, 3), got {points.shape}")

    n_points = points.shape[0]

    if owner not in OWNER_RULES:
        raise ValueError(
            f"unknown owner rule {owner!r}; accepted: {OWNER_RULES}"
        )
    if candidate_cells is not None and owner != "lowest_original_cell":
        raise ValueError(
            "candidate_cells is valid only with owner='lowest_original_cell' "
            f"(got owner={owner!r})"
        )
    if parent is not None and owner != "lowest_original_cell":
        raise ValueError(
            "parent is valid only with owner='lowest_original_cell' "
            f"(got owner={owner!r})"
        )

    # OPS-40: this routine is collective — rank 0's output buffer is sized by
    # ``n_points`` and every rank's hits are scattered into it by index, so a
    # caller that hands each rank its own *local* points corrupts the gather
    # (it died with ``IndexError: index 36378 is out of bounds for axis 0 with
    # size 36378`` at -n 2, EX-52, 20260907T033551Z_EX-52.log:400) and passes
    # silently at -n 1. Refuse it on EVERY rank — a rank-0-only raise ahead of
    # the gather below would hang the other ranks.
    #
    # OPS-55: equal counts are not enough — with equal counts and different
    # coordinates each rank locates ITS OWN rows and rank 0 scatters the hits
    # by positional index, so ``values[i]`` silently becomes the field at some
    # other rank's row ``i``. A digest rides the same allgather (no extra round
    # on the valid path). Validate, never broadcast: a broadcast would hide a
    # caller that still indexes the result against its own local ``points``.
    # The contract is bitwise identity.
    points = np.ascontiguousarray(points)
    digest = hashlib.blake2b(points.tobytes(), digest_size=16).digest()
    # POST-8 step 2: ``owner`` rides the same allgather — the rule is
    # collective, so a rank-mismatched ``owner`` must raise on every rank
    # rather than hang.
    # POST-1 step 7a: whether ``candidate_cells`` is given joins the key (the
    # set itself is rank-local by construction and cannot be compared).
    # POST-10 step 1: so does whether ``parent`` is given.
    gathered_keys = comm.allgather(
        (n_points, digest, (owner, candidate_cells is None, parent is None))
    )
    counts = [c for c, _, _ in gathered_keys]
    if len(set(counts)) > 1:
        per_rank = ", ".join(f"rank {r}: {c}" for r, c in enumerate(counts))
        raise ValueError(
            "evaluate_vector_field_parallel received a different number of "
            f"points on different ranks ({per_rank}). This routine is "
            "collective: the point list must be identical on every rank — "
            "allgather your points first."
        )
    if len({d for _, d, _ in gathered_keys}) > 1:
        # Every rank knows the digests differ, so further collectives are safe
        # here: compare against rank 0's list and raise on EVERY rank.
        ref = comm.bcast(points if comm.rank == 0 else None, root=0)
        diff = np.abs(points - ref)
        rows = np.flatnonzero(np.any(points != ref, axis=1))
        first_row = int(rows[0]) if rows.size else -1
        max_dev = float(np.max(diff)) if diff.size else 0.0
        stats = comm.allgather((first_row, max_dev))
        per_rank = "; ".join(
            f"rank {r}: first differing row {row}, max abs deviation {dev:.3e}"
            for r, (row, dev) in enumerate(stats)
        )
        raise ValueError(
            "evaluate_vector_field_parallel received point lists of equal "
            "length but different coordinates on different ranks (vs rank 0 — "
            f"{per_rank}). This routine is collective: the point list must be "
            "bitwise identical on every rank — broadcast the points from "
            "rank 0 (a gross deviation means rank-local points; ~1e-16 means "
            "rank-dependent floating-point ordering)."
        )
    owners = [o for _, _, o in gathered_keys]
    if len(set(owners)) > 1:
        per_rank = ", ".join(f"rank {r}: {o!r}" for r, o in enumerate(owners))
        raise ValueError(
            "evaluate_vector_field_parallel received different owner rules on "
            f"different ranks ({per_rank}). This routine is collective: "
            "``owner`` (and whether ``candidate_cells`` and ``parent`` are "
            "given) must be identical on every rank."
        )

    value_shape = function.function_space.element.value_shape
    value_size = int(np.prod(value_shape)) if len(value_shape) > 0 else 1

    # Follow the function's own scalar dtype: in the complex DolfinX build
    # ``function.eval`` returns complex128, and gathering it into a float64
    # buffer raises a casting error (or, worse, would drop the phase).
    scalar_dtype = np.asarray(function.x.array).dtype

    local_point_idx: list[int] = []
    local_cells: list[int] = []
    if owner is None:
        bb_tree = geometry.bb_tree(mesh, mesh.topology.dim)
        candidate_cells = geometry.compute_collisions_points(bb_tree, points)
        colliding_cells = geometry.compute_colliding_cells(mesh, candidate_cells, points)

        for i in range(n_points):
            links = colliding_cells.links(i)
            if len(links) > 0:
                local_point_idx.append(i)
                local_cells.append(links[0])
    else:  # "lowest_original_cell"
        _, owner_local_cell = point_owner_cells(
            mesh, points, comm, candidate_cells=candidate_cells, parent=parent
        )
        for i in range(n_points):
            if owner_local_cell[i] >= 0:
                local_point_idx.append(i)
                local_cells.append(int(owner_local_cell[i]))

    local_point_idx_arr = np.asarray(local_point_idx, dtype=np.int32)

    if local_point_idx_arr.size > 0:
        local_points = points[local_point_idx_arr]
        local_cells_arr = np.asarray(local_cells, dtype=np.int32)
        local_values = function.eval(local_points, local_cells_arr)
    else:
        local_values = np.zeros((0, value_size), dtype=scalar_dtype)

    gathered_values = comm.gather(local_values, root=0)
    gathered_indices = comm.gather(local_point_idx_arr, root=0)

    if comm.rank == 0:
        values = np.zeros((n_points, value_size), dtype=scalar_dtype)
        valid_mask = np.zeros(n_points, dtype=bool)
        for rank_values, rank_indices in zip(gathered_values, gathered_indices):
            if rank_indices.size > 0:
                values[rank_indices] = rank_values
                valid_mask[rank_indices] = True
    else:
        values = None
        valid_mask = None

    values = comm.bcast(values, root=0)
    valid_mask = comm.bcast(valid_mask, root=0)
    return values, valid_mask
