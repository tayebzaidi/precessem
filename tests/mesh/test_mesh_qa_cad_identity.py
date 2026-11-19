"""`GEO-2` / `GEO-6` / `GEO-7` (§9 item 184): one CAD-identity QA utility.

``fem_em_solver.io.mesh_qa.cad_identity_report`` gathers the volume identities
the mesh gates already assert one by one — tag set, tagged-volume partition,
analytic box, exact (planar) region volumes, recorded volumes — into one
report. This file checks that the report reproduces those gates' own readings
on the three fixtures they gate (anchors A1/A2), and that a planted defect on
each built mesh raises exactly its own flag (negative controls N1-N3).

Every constant, generator call argument and analytic reference is imported from
the module that owns it; none is restated here.
"""

from __future__ import annotations

import dolfinx
import numpy as np
import pytest
from mpi4py import MPI

from fem_em_solver.io import CAD_IDENTITY_RTOL, cad_identity_report
from fem_em_solver.io.mesh import MeshGenerator
from tests.mesh import helpers
from tests.mesh import test_birdcage_port_tags as birdcage
from tests.mesh import test_coil_phantom_conforming as coil_phantom_conforming
from tests.mesh import test_two_torus_conforming as two_torus
from tests.mesh.helpers import REQUIRED_COIL_PHANTOM_TAGS, VOLUME_PARTITION_BAND
from tests.mesh.test_birdcage_conductor_sizing import BASELINE_CONTROL_RESOLUTION
from tests.mesh.test_mesh_tag_integrity import GEOMETRY, UNIFORM_VOLUMES_RECORD

# Recorded cell count of the birdcage box-port mesh at the baseline control
# resolution (`20260826T093202Z_GEO-21-step2-gate.log:11075`): printed, not
# asserted.
BIRDCAGE_BASELINE_CELLS_RECORD = 33185


def _global_cells(mesh) -> int:
    return mesh.topology.index_map(mesh.topology.dim).size_global


def _print(comm, msg):
    if comm.rank == 0:
        print(msg, flush=True)


def _print_report(comm, label, report):
    _print(
        comm,
        f"\n[cad-identity] {label}: flags={sorted(report.flags)} passed={report.passed}"
        f"\n  tag_volumes={ {t: f'{v:.9e}' for t, v in report.tag_volumes.items()} }"
        f"\n  total_volume={report.total_volume:.9e}"
        f" missing={sorted(report.missing_tags)} unexpected={sorted(report.unexpected_tags)}"
        f"\n  partition_rel={report.partition_rel:.3e} box_rel={report.box_rel}"
        f"\n  exact_rel={ {t: f'{r:.3e}' for t, r in report.exact_rel.items()} }"
        f"\n  record_rel={ {t: f'{r:.3e}' for t, r in report.record_rel.items()} }",
    )


def _replant(mesh, cell_tags, indices, values):
    order = np.argsort(indices)
    idx = np.asarray(indices, dtype=np.int32)[order]
    val = np.asarray(values, dtype=np.int32)[order]
    assert np.all(np.diff(idx) > 0), "planted indices must be sorted and unique"
    return dolfinx.mesh.meshtags(mesh, mesh.topology.dim, idx, val)


def _owned_mask(mesh, cell_tags):
    size_local = mesh.topology.index_map(mesh.topology.dim).size_local
    return np.asarray(cell_tags.indices) < size_local


# --- fixtures: each mesh built once -------------------------------------------------


@pytest.fixture(scope="module")
def coil_phantom():
    mesh, cell_tags, _ = MeshGenerator.coil_phantom_domain(comm=MPI.COMM_WORLD, **GEOMETRY)
    return mesh, cell_tags


@pytest.fixture(scope="module")
def two_torus_mesh():
    mesh, cell_tags, _ = MeshGenerator.two_torus_domain(
        separation=two_torus.SEPARATION,
        major_radius=two_torus.MAJOR_RADIUS,
        minor_radius=two_torus.MINOR_RADIUS,
        resolution=two_torus.RESOLUTION,
        air_padding=two_torus.AIR_PADDING,
        wire_resolution=two_torus.WIRE_RESOLUTION,
        far_resolution=two_torus.RESOLUTION,
        comm=MPI.COMM_WORLD,
    )
    return mesh, cell_tags


@pytest.fixture(scope="module")
def birdcage_mesh():
    mesh, cell_tags, _, _ = MeshGenerator.birdcage_port_domain(
        leg_count=birdcage.LEG_COUNT,
        ring_radius=birdcage.RING_RADIUS,
        leg_width=birdcage.LEG_WIDTH,
        leg_spacing=birdcage.LEG_SPACING,
        coil_length=birdcage.COIL_LENGTH,
        ring_minor_radius=birdcage.RING_MINOR_RADIUS,
        phantom_radius=birdcage.PHANTOM_RADIUS,
        phantom_height=birdcage.PHANTOM_HEIGHT,
        port_box_size=birdcage.PORT_BOX_SIZE,
        air_padding=birdcage.AIR_PADDING,
        resolution=birdcage.RESOLUTION,
        conductor_resolution=BASELINE_CONTROL_RESOLUTION,
        comm=MPI.COMM_WORLD,
        return_diagnostics=True,
    )
    return mesh, cell_tags


def _birdcage_port_tags():
    return [100 + i for i in range(1, birdcage.LEG_COUNT + 1)]


def _birdcage_expected_tags():
    return [*birdcage.BIRDCAGE_CORE_TAGS, *_birdcage_port_tags()]


def _birdcage_exact_volumes():
    v_port = birdcage.PORT_BOX_SIZE[0] * birdcage.PORT_BOX_SIZE[1] * birdcage.PORT_BOX_SIZE[2]
    return {tag: v_port for tag in _birdcage_port_tags()}


def _r1_report(mesh, cell_tags, **kwargs):
    return cad_identity_report(
        mesh,
        cell_tags,
        expected_tags=REQUIRED_COIL_PHANTOM_TAGS,
        box_volume=coil_phantom_conforming._analytic_box_volume(),
        **kwargs,
    )


# --- A0 -----------------------------------------------------------------------------


def test_default_rtol_is_the_partition_band():
    """(A0) The utility's default band is the mesh gates' partition band."""
    import inspect

    default = inspect.signature(cad_identity_report).parameters["rtol"].default
    assert default == VOLUME_PARTITION_BAND
    assert CAD_IDENTITY_RTOL == VOLUME_PARTITION_BAND


# --- A1 + A2 ------------------------------------------------------------------------


def _assert_bitwise_with_helpers(mesh, cell_tags, report):
    comm = mesh.comm
    for tag, v in report.tag_volumes.items():
        assert v == helpers.tag_volume(mesh, cell_tags, tag, comm), f"tag {tag}"
    assert report.total_volume == helpers.total_volume(mesh, comm)


def test_r1_coil_phantom_identities(coil_phantom):
    """(A1, A2-R1) coil+phantom: tags, record at 1e-9, box, partition."""
    mesh, cell_tags = coil_phantom
    report = _r1_report(mesh, cell_tags, record_volumes=UNIFORM_VOLUMES_RECORD)
    _print_report(mesh.comm, "R1 coil_phantom_domain", report)
    _print(mesh.comm, f"  cells={_global_cells(mesh)}")
    _assert_bitwise_with_helpers(mesh, cell_tags, report)
    assert report.flags == frozenset(), report
    assert report.passed


def test_r2_two_torus_identities(two_torus_mesh):
    """(A1, A2-R2) two-torus: tag set {1,2,3}, box, partition."""
    mesh, cell_tags = two_torus_mesh
    report = cad_identity_report(
        mesh,
        cell_tags,
        expected_tags=[*two_torus.WIRE_TAGS, two_torus.AIR_TAG],
        box_volume=two_torus._analytic_box_volume(
            two_torus.SEPARATION,
            two_torus.MAJOR_RADIUS,
            two_torus.MINOR_RADIUS,
            two_torus.AIR_PADDING,
        ),
    )
    _print_report(mesh.comm, "R2 two_torus_domain", report)
    _print(mesh.comm, f"  cells={_global_cells(mesh)}")
    _assert_bitwise_with_helpers(mesh, cell_tags, report)
    assert report.flags == frozenset(), report
    assert report.passed


def test_r3_birdcage_box_port_identities(birdcage_mesh):
    """(A1, A2-R3) birdcage box-port at the baseline control resolution."""
    mesh, cell_tags = birdcage_mesh
    report = cad_identity_report(
        mesh,
        cell_tags,
        expected_tags=_birdcage_expected_tags(),
        box_volume=birdcage._analytic_box_volume(),
        exact_volumes=_birdcage_exact_volumes(),
    )
    _print_report(mesh.comm, "R3 birdcage_port_domain", report)
    _print(
        mesh.comm,
        f"  cells={_global_cells(mesh)} (record {BIRDCAGE_BASELINE_CELLS_RECORD}, printed not asserted)",
    )
    _assert_bitwise_with_helpers(mesh, cell_tags, report)
    assert report.flags == frozenset(), report
    assert report.passed


# --- negative controls --------------------------------------------------------------


def test_n1_merged_tag_raises_only_tag_set(coil_phantom):
    """(N1) Every tag-3 value relabelled 4: flags == {tag_set}, missing {3}."""
    mesh, cell_tags = coil_phantom
    values = np.array(cell_tags.values, copy=True)
    values[values == 3] = 4
    planted = _replant(mesh, cell_tags, np.array(cell_tags.indices, copy=True), values)
    report = _r1_report(mesh, planted)
    _print_report(mesh.comm, "N1 tag merged (3 -> 4)", report)
    assert report.flags == frozenset({"tag_set"})
    assert report.missing_tags == frozenset({3})
    assert report.unexpected_tags == frozenset()
    assert report.partition_rel <= VOLUME_PARTITION_BAND
    assert report.box_rel <= VOLUME_PARTITION_BAND
    assert not report.passed


def _drop_or_retag_one_owned_cell(mesh, cell_tags, tag, new_value):
    """Plant on the lowest rank owning a ``tag`` cell; return (meshtags, rank)."""
    comm = mesh.comm
    owned = _owned_mask(mesh, cell_tags)
    hits = np.flatnonzero(owned & (np.asarray(cell_tags.values) == tag))
    owners = comm.allgather(hits.size > 0)
    rank = min(r for r, has in enumerate(owners) if has)
    indices = np.array(cell_tags.indices, copy=True)
    values = np.array(cell_tags.values, copy=True)
    if comm.rank == rank:
        k = hits[0]
        if new_value is None:
            indices = np.delete(indices, k)
            values = np.delete(values, k)
        else:
            values[k] = new_value
    return _replant(mesh, cell_tags, indices, values), rank


def test_n2_ungrouped_cell_raises_only_partition(coil_phantom):
    """(N2) One owned tag-3 cell dropped: flags == {partition}."""
    mesh, cell_tags = coil_phantom
    planted, rank = _drop_or_retag_one_owned_cell(mesh, cell_tags, 3, None)
    report = _r1_report(mesh, planted)
    _print_report(mesh.comm, f"N2 ungrouped cell (rank {rank})", report)
    predicted_floor = 100.0 * VOLUME_PARTITION_BAND
    _print(
        mesh.comm,
        f"  N2 miss partition_rel={report.partition_rel:.3e} "
        f"PREDICTED >= {predicted_floor:.1e}: "
        f"{'met' if report.partition_rel >= predicted_floor else 'MISSED'}",
    )
    assert report.flags == frozenset({"partition"})
    assert not report.passed


def test_n3_port_box_retagged_raises_only_exact_volume(birdcage_mesh):
    """(N3, `GEO-2`) One owned tag-101 cell re-tagged 2: flags == {exact_volume}, 101 only."""
    mesh, cell_tags = birdcage_mesh
    planted, rank = _drop_or_retag_one_owned_cell(mesh, cell_tags, 101, 2)
    report = cad_identity_report(
        mesh,
        planted,
        expected_tags=_birdcage_expected_tags(),
        box_volume=birdcage._analytic_box_volume(),
        exact_volumes=_birdcage_exact_volumes(),
    )
    _print_report(mesh.comm, f"N3 port box re-tagged (rank {rank})", report)
    predicted_floor = 1.0e-3
    miss = report.exact_rel[101]
    _print(
        mesh.comm,
        f"  N3 miss exact_rel[101]={miss:.3e} PREDICTED >= {predicted_floor:.1e}: "
        f"{'met' if miss >= predicted_floor else 'MISSED'}",
    )
    assert report.flags == frozenset({"exact_volume"})
    over = {t for t, r in report.exact_rel.items() if r > report.rtol}
    assert over == {101}
    assert not report.passed
