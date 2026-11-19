"""Mesh tests for coarse coil + phantom + air geometry presets."""

import numpy as np
import pytest
from mpi4py import MPI

from fem_em_solver.io.mesh import MeshGenerator
from tests.mesh.helpers import (
    REQUIRED_COIL_PHANTOM_TAGS,
    assert_required_tags_nonempty,
    compute_tag_cell_centroid,
)


def _generate_default_mesh(**kwargs):
    return MeshGenerator.coil_phantom_domain(
        coil_major_radius=0.08,
        coil_minor_radius=0.01,
        coil_separation=0.08,
        phantom_radius=0.04,
        phantom_height=0.10,
        air_padding=0.04,
        resolution=0.015,
        comm=MPI.COMM_WORLD,
        **kwargs,
    )


def test_coil_phantom_mesh_generates_required_tags_centered_preset():
    """Centered preset should generate required tags, including non-empty phantom."""
    mesh, cell_tags, _ = _generate_default_mesh(phantom_placement_preset="centered")

    n_cells = mesh.topology.index_map(3).size_global
    assert n_cells > 0, "Mesh must contain cells"
    assert n_cells < 50000, f"Mesh must stay under 50k cells, got {n_cells}"

    assert_required_tags_nonempty(cell_tags, REQUIRED_COIL_PHANTOM_TAGS, comm=MPI.COMM_WORLD)

    phantom_centroid = compute_tag_cell_centroid(mesh, cell_tags, tag=3, comm=MPI.COMM_WORLD)
    assert abs(phantom_centroid[0]) < 1.0e-2
    assert abs(phantom_centroid[1]) < 1.0e-2


def test_coil_phantom_mesh_off_center_preset_moves_phantom_without_overlap():
    """Off-center preset should keep phantom non-empty and displaced in +x direction."""
    mesh, cell_tags, _ = _generate_default_mesh(phantom_placement_preset="off_center")

    assert_required_tags_nonempty(cell_tags, REQUIRED_COIL_PHANTOM_TAGS, comm=MPI.COMM_WORLD)

    phantom_centroid = compute_tag_cell_centroid(mesh, cell_tags, tag=3, comm=MPI.COMM_WORLD)
    assert (
        phantom_centroid[0] > 5.0e-3
    ), "Off-center preset should produce positive x displacement for phantom centroid"


_PHANTOM_RADIUS = 0.04
_PHANTOM_HEIGHT = 0.10

# GEO-3 (§9 item 73): placements and the cylinder centre each one places.
_GEO3_PLACEMENTS = {
    "centered": (dict(phantom_placement_preset="centered"), (0.0, 0.0)),
    "off_center": (dict(phantom_placement_preset="off_center"), (0.35 * _PHANTOM_RADIUS, 0.0)),
    "override_neg_y": (
        dict(phantom_placement_preset="centered", phantom_offset_xy=(0.0, -0.35 * _PHANTOM_RADIUS)),
        (0.0, -0.35 * _PHANTOM_RADIUS),
    ),
}


@pytest.fixture(scope="module")
def geo3_placed_meshes():
    """Mesh the three GEO-3 placements once for the module."""
    out = {}
    for name, (kwargs, centre_xy) in _GEO3_PLACEMENTS.items():
        mesh, cell_tags, _ = _generate_default_mesh(**kwargs)
        out[name] = (mesh, cell_tags, centre_xy)
    return out


def test_coil_phantom_placed_phantom_first_moment_is_cylinder_centre(geo3_placed_meshes):
    """GEO-3: tag 3's volume-weighted first moment is the placed cylinder's centre.

    Anchor: C (the placed cylinder) is convex and the mesh nodes lie on its
    surface, so M (tag 3) ⊂ C.  With D = C∖M and δ = V_C/V_M − 1,
    m_M − c_C = δ(c_C − c_D), hence |m_xy − c_xy| ≤ δ·r and |m_z| ≤ δ·h/2.
    """
    import ufl
    from dolfinx import fem

    comm = MPI.COMM_WORLD
    r, h = _PHANTOM_RADIUS, _PHANTOM_HEIGHT
    v_c = np.pi * r**2 * h
    results = {}
    for name, (mesh, cell_tags, (cx, cy)) in geo3_placed_meshes.items():
        dx3 = ufl.Measure(
            "dx", domain=mesh, subdomain_data=cell_tags, metadata={"quadrature_degree": 2}
        )(3)
        x = ufl.SpatialCoordinate(mesh)
        one = fem.Constant(mesh, 1.0)
        v_m = comm.allreduce(fem.assemble_scalar(fem.form(one * dx3)), op=MPI.SUM).real
        moment = np.array(
            [
                comm.allreduce(fem.assemble_scalar(fem.form(x[i] * dx3)), op=MPI.SUM).real
                for i in range(3)
            ]
        )
        m = moment / v_m
        delta = v_c / v_m - 1.0
        bound_xy = delta * r
        bound_z = delta * h / 2.0
        miss_x, miss_y, miss_z = abs(m[0] - cx), abs(m[1] - cy), abs(m[2])
        cell_centroid = compute_tag_cell_centroid(mesh, cell_tags, tag=3, comm=comm)
        if comm.rank == 0:
            print(
                f"\nGEO-3 [{name}] c=({cx:.6e},{cy:.6e},0) V_C={v_c:.9e} V_M={v_m:.9e} delta={delta:.6f} "
                f"m=({m[0]:.6e},{m[1]:.6e},{m[2]:.6e}) bound_xy={bound_xy:.3e} bound_z={bound_z:.3e} "
                f"ratios x={miss_x / bound_xy:.4f} y={miss_y / bound_xy:.4f} z={miss_z / bound_z:.4f} "
                f"cell_count_centroid=({cell_centroid[0]:.6e},{cell_centroid[1]:.6e},{cell_centroid[2]:.6e})"
            )
        results[name] = (delta, m, bound_xy, cell_centroid)

        # Premise: M ⊂ C with a small discretisation deficit.
        assert (
            0.0 < delta < 0.1
        ), f"[{name}] premise failed: delta={delta:.6f} (V_M={v_m:.9e}, V_C={v_c:.9e})"
        # Band: every moment inside its own in-run bound.
        assert miss_x <= bound_xy, f"[{name}] |m_x - c_x|={miss_x:.3e} > {bound_xy:.3e}"
        assert miss_y <= bound_xy, f"[{name}] |m_y - c_y|={miss_y:.3e} > {bound_xy:.3e}"
        assert miss_z <= bound_z, f"[{name}] |m_z|={miss_z:.3e} > {bound_z:.3e}"

    # Negative control (comparand swap): centred m_x against the off-centre 0.014 m.
    delta_c, m_c, bound_c, cc_c = results["centered"]
    swap_miss = abs(m_c[0] - 0.35 * r)
    if comm.rank == 0:
        print(
            f"\nGEO-3 control: miss={swap_miss:.6e} bound={bound_c:.3e} separation={swap_miss / bound_c:.2f}x "
            f"(predicted delta 0.016744, bound 6.70e-4 m, separation ~20.9x, ceiling 21.9x); "
            f"measured delta={delta_c:.6f}; cell-count centroid x={cc_c[0]:.6e}"
        )
    assert (
        swap_miss > bound_c
    ), f"control: swapped miss {swap_miss:.3e} not above bound {bound_c:.3e}"


def test_coil_phantom_mesh_rejects_overlapping_off_center_placement():
    """Explicit offsets that overlap the coil envelope should be rejected."""
    with pytest.raises(ValueError, match="Phantom overlaps coil conductor envelope"):
        _generate_default_mesh(
            phantom_placement_preset="off_center",
            phantom_offset_xy=(0.04, 0.0),
        )
