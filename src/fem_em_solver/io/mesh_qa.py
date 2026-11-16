"""Mesh QA helpers for quick tag-integrity diagnostics."""

from __future__ import annotations

from typing import Dict, Mapping, Optional

import numpy as np
from mpi4py import MPI


def cell_tag_counts(cell_tags, comm: Optional[MPI.Intracomm] = None) -> Dict[int, int]:
    """Return global cell counts per cell-tag value.

    Parameters
    ----------
    cell_tags : dolfinx.mesh.MeshTags
        Cell tag data returned from mesh generation.
    comm : MPI.Intracomm, optional
        MPI communicator. Defaults to ``cell_tags.topology.comm`` when available.
    """
    if comm is None:
        comm = cell_tags.topology.comm

    local_values = np.asarray(cell_tags.values, dtype=np.int64)

    local_max = np.array([-1], dtype=np.int64)
    if local_values.size:
        local_max[0] = int(np.max(local_values))

    global_max = np.array([-1], dtype=np.int64)
    comm.Allreduce(local_max, global_max, op=MPI.MAX)

    if global_max[0] < 0:
        return {}

    local_hist = np.zeros(int(global_max[0]) + 1, dtype=np.int64)
    if local_values.size:
        local_hist = np.bincount(local_values, minlength=local_hist.size).astype(np.int64)

    global_hist = np.zeros_like(local_hist)
    comm.Allreduce(local_hist, global_hist, op=MPI.SUM)

    return {int(tag): int(count) for tag, count in enumerate(global_hist) if count > 0}


def format_cell_tag_summary(
    counts: Mapping[int, int],
    tag_names: Optional[Mapping[int, str]] = None,
) -> str:
    """Build a compact deterministic summary string for tag counts."""
    if not counts:
        return "(no tagged cells)"

    parts = []
    for tag in sorted(counts):
        name = tag_names.get(tag, f"tag_{tag}") if tag_names else f"tag_{tag}"
        parts.append(f"{name}={counts[tag]}")

    return ", ".join(parts)


def print_cell_tag_summary(
    cell_tags,
    tag_names: Optional[Mapping[int, str]] = None,
    comm: Optional[MPI.Intracomm] = None,
    prefix: str = "[mesh] ",
) -> Dict[int, int]:
    """Print global cell counts by tag on rank 0 and return the counts."""
    if comm is None:
        comm = cell_tags.topology.comm

    counts = cell_tag_counts(cell_tags, comm=comm)

    if comm.rank == 0:
        print(f"{prefix}cell-tag counts: {format_cell_tag_summary(counts, tag_names=tag_names)}")

    return counts


def format_expected_tag_counts(
    counts: Mapping[int, int],
    required_tags: Mapping[int, str],
    minimum_count: int = 1,
) -> str:
    """Format expected-vs-actual counts for required tags.

    Parameters
    ----------
    counts : Mapping[int, int]
        Observed global cell counts by tag.
    required_tags : Mapping[int, str]
        Required tags and display names.
    minimum_count : int, default=1
        Minimum required count for each required tag.
    """
    if not required_tags:
        return "(no required tags)"

    rendered = []
    for tag in sorted(required_tags):
        name = required_tags[tag]
        actual = int(counts.get(tag, 0))
        status = "OK" if actual >= minimum_count else "MISSING"
        rendered.append(
            f"{name}(tag={tag}): expected>={minimum_count}, actual={actual} [{status}]"
        )

    return "; ".join(rendered)


def print_required_tag_failure_summary(
    counts: Mapping[int, int],
    required_tags: Mapping[int, str],
    comm: Optional[MPI.Intracomm] = None,
    prefix: str = "[mesh-qa] ",
) -> None:
    """Print compact diagnostics for required-tag failures on rank 0.

    Includes expected-vs-actual required tag counts and a compact summary of
    all observed tags to speed up debugging for missing/empty tag failures.
    """
    if comm is not None and comm.rank != 0:
        return

    print(
        f"{prefix}required-tag expected vs actual: "
        f"{format_expected_tag_counts(counts, required_tags)}"
    )
    print(
        f"{prefix}observed-tag summary: "
        f"{format_cell_tag_summary(counts, tag_names=required_tags)}"
    )


# `GEO-2` / `GEO-6` / `GEO-7` (§9 item 184, 2026-10-03): one CAD-identity QA
# utility.  The default band is the exact-arithmetic tagged-volume partition
# band of tests/mesh/helpers.py (``VOLUME_PARTITION_BAND``); the test asserts
# the two are equal (anchor A0) rather than this module importing from tests.
CAD_IDENTITY_RTOL = 1.0e-9


from dataclasses import dataclass, field  # noqa: E402
from typing import FrozenSet, Iterable  # noqa: E402


@dataclass(frozen=True)
class CadIdentityReport:
    """Volume identities of a tagged mesh against its CAD (see ``cad_identity_report``)."""

    tag_volumes: Dict[int, float]
    total_volume: float
    missing_tags: FrozenSet[int]
    unexpected_tags: FrozenSet[int]
    partition_rel: float
    box_rel: Optional[float]
    exact_rel: Dict[int, float]
    record_rel: Dict[int, float]
    flags: FrozenSet[str]
    rtol: float = CAD_IDENTITY_RTOL
    passed: bool = field(default=False)


def _owned_tag_set(mesh, cell_tags, comm) -> set:
    """Global tag set from *owned* cells only (ghosted meshes repeat ghosts)."""
    tdim = mesh.topology.dim
    size_local = mesh.topology.index_map(tdim).size_local
    owned = np.asarray(cell_tags.indices) < size_local
    local = set(int(v) for v in np.unique(np.asarray(cell_tags.values)[owned]).tolist())
    return set().union(*comm.allgather(local))


def cad_identity_report(
    mesh,
    cell_tags,
    *,
    expected_tags: Iterable[int],
    box_volume: Optional[float] = None,
    exact_volumes: Optional[Mapping[int, float]] = None,
    record_volumes: Optional[Mapping[int, float]] = None,
    rtol: float = CAD_IDENTITY_RTOL,
    comm: Optional[MPI.Intracomm] = None,
) -> CadIdentityReport:
    """Check a tagged mesh's cell-volume identities; return a frozen report.

    Flags (a subset of ``tag_set`` / ``partition`` / ``box`` /
    ``exact_volume`` / ``record``):

    * ``tag_set`` — the global owned-cell tag set differs from ``expected_tags``;
    * ``partition`` — ``|sum_expected V_t / V_mesh - 1| > rtol``;
    * ``box`` — ``|V_mesh / box_volume - 1| > rtol``;
    * ``exact_volume`` — some ``|V_t / exact_volumes[t] - 1| > rtol``;
    * ``record`` — some ``|V_t / record_volumes[t] - 1| > rtol``.

    Volumes are assembled with the expressions of tests/mesh/helpers.py
    (``tag_volume`` / ``total_volume``) and reduced over ``comm``.
    """
    import ufl
    from dolfinx import default_scalar_type, fem

    comm = mesh.comm if comm is None else comm
    expected = sorted(int(t) for t in expected_tags)

    present = _owned_tag_set(mesh, cell_tags, comm)
    missing = frozenset(set(expected) - present)
    unexpected = frozenset(present - set(expected))

    tag_volumes: Dict[int, float] = {}
    for tag in expected:
        dx_tag = ufl.Measure(
            "dx", domain=mesh, subdomain_data=cell_tags, subdomain_id=(int(tag),)
        )
        one = fem.Constant(mesh, default_scalar_type(1.0))
        local = fem.assemble_scalar(fem.form(one * dx_tag))
        tag_volumes[tag] = float(np.real(comm.allreduce(local, op=MPI.SUM)))

    one = fem.Constant(mesh, default_scalar_type(1.0))
    local = fem.assemble_scalar(fem.form(one * ufl.dx))
    total = float(np.real(comm.allreduce(local, op=MPI.SUM)))

    partition_rel = abs(sum(tag_volumes.values()) / total - 1.0)
    box_rel = None if box_volume is None else abs(total / float(box_volume) - 1.0)

    def _rel(reference):
        return {
            int(t): abs(tag_volumes.get(int(t), 0.0) / float(v) - 1.0)
            for t, v in (reference or {}).items()
        }

    exact_rel = _rel(exact_volumes)
    record_rel = _rel(record_volumes)

    flags = set()
    if missing or unexpected:
        flags.add("tag_set")
    if partition_rel > rtol:
        flags.add("partition")
    if box_rel is not None and box_rel > rtol:
        flags.add("box")
    if any(r > rtol for r in exact_rel.values()):
        flags.add("exact_volume")
    if any(r > rtol for r in record_rel.values()):
        flags.add("record")

    return CadIdentityReport(
        tag_volumes=tag_volumes,
        total_volume=total,
        missing_tags=missing,
        unexpected_tags=unexpected,
        partition_rel=partition_rel,
        box_rel=box_rel,
        exact_rel=exact_rel,
        record_rel=record_rel,
        flags=frozenset(flags),
        rtol=rtol,
        passed=not flags,
    )
