"""Private sheet-geometry helpers for ``study.FourPortBirdcageGeometry`` (`WF-16`).

Verbatim copies, in operation order, of the test helpers the F-small four-port
rung composes (``src/`` never imports ``tests/``, `WF-13`'s rule):

* ``_sheet_nodes``, ``_projected_extents``, ``_sheet_azimuth_deg``,
  ``_sheet_areas`` — ``tests/mesh/test_birdcage_leg_offset.py``;
* ``_port_frame``, ``_narrowed_radial`` —
  ``tests/validation/test_port_birdcage_leg_offset_sweep.py``;
* ``_sheet_facet_count`` — ``tests/mesh/test_two_torus_port_sheet.py``;
* ``_facet_group_area`` — ``tests/mesh/test_two_torus_port_facets.py``;
* ``_global_facet_count``, ``_interface_area_or_zero`` —
  ``tests/mesh/test_birdcage_port_terminals.py``.

Extents are reduced elementwise through ``allgather`` (never ``allreduce`` on
arrays); facet counts and areas are reduced. Nothing here is public.
"""

from __future__ import annotations

import dolfinx
import numpy as np
import ufl
from dolfinx import default_scalar_type, fem
from mpi4py import MPI

from ..io.mesh import _interface_facet_tags


def _sheet_nodes(msh, facet_tags, tag):
    fdim = msh.topology.dim - 1
    facets = np.asarray(facet_tags.find(tag), dtype=np.int32)
    if not facets.size:
        return np.zeros((0, 3))
    nodes = dolfinx.cpp.mesh.entities_to_geometry(msh._cpp_object, fdim, facets, False)
    return msh.geometry.x[nodes.reshape(-1)]


def _projected_extents(msh, facet_tags, tag, comm, directions):
    pts = _sheet_nodes(msh, facet_tags, tag)
    directions = np.asarray(directions, dtype=float)
    if pts.shape[0]:
        proj = pts @ directions.T
        lo, hi = proj.min(axis=0), proj.max(axis=0)
    else:
        lo = np.full(directions.shape[0], np.inf)
        hi = np.full(directions.shape[0], -np.inf)
    lo = np.min(np.array(comm.allgather(lo)), axis=0)
    hi = np.max(np.array(comm.allgather(hi)), axis=0)
    return hi - lo, 0.5 * (hi + lo)


def _sheet_azimuth_deg(msh, facet_tags, tag, comm):
    pts = _sheet_nodes(msh, facet_tags, tag)
    if pts.shape[0]:
        lo, hi = pts.min(axis=0), pts.max(axis=0)
    else:
        lo, hi = np.full(3, np.inf), np.full(3, -np.inf)
    lo = np.min(np.array(comm.allgather(lo)), axis=0)
    hi = np.max(np.array(comm.allgather(hi)), axis=0)
    centre = 0.5 * (lo + hi)
    return float(np.degrees(np.arctan2(centre[1], centre[0])) % 360.0)


def _global_facet_count(mesh, facet_tags, tag, comm) -> int:
    fdim = mesh.topology.dim - 1
    n_owned = mesh.topology.index_map(fdim).size_local
    indices = np.asarray(facet_tags.indices)
    values = np.asarray(facet_tags.values)
    local = int(np.count_nonzero((values == tag) & (indices < n_owned)))
    return int(comm.allreduce(local, op=MPI.SUM))


def _facet_group_area(msh, facet_tags, tag, comm):
    msh.topology.create_entity_permutations()

    dS_tag = ufl.Measure("dS", domain=msh, subdomain_data=facet_tags, subdomain_id=(tag,))
    one = fem.Constant(msh, default_scalar_type(1.0))
    local = fem.assemble_scalar(fem.form(ufl.avg(one) * dS_tag))
    return float(np.real(comm.allreduce(local, op=MPI.SUM)))


def _interface_area_or_zero(mesh, facet_tags, tag, comm) -> float:
    if _global_facet_count(mesh, facet_tags, tag, comm) == 0:
        return 0.0
    return _facet_group_area(mesh, facet_tags, tag, comm)


def _sheet_areas(msh, cell_tags, ports, comm, sheet_iface, port_lower, port_upper):
    tags_f = _interface_facet_tags(
        msh,
        cell_tags,
        {sheet_iface + i: (port_lower + i, port_upper + i) for i in ports},
    )
    areas = {i: _interface_area_or_zero(msh, tags_f, sheet_iface + i, comm) for i in ports}
    return tags_f, areas


def _sheet_facet_count(msh, facet_tags, tag, comm):
    fdim = msh.topology.dim - 1
    n_local = msh.topology.index_map(fdim).size_local
    facets = facet_tags.find(tag)
    return int(comm.allreduce(int(np.count_nonzero(facets < n_local)), op=MPI.SUM))


def _port_frame(azimuth_deg):
    beta = np.radians(azimuth_deg)
    radial = np.array([np.cos(beta), np.sin(beta), 0.0])
    azimuthal = np.array([-np.sin(beta), np.cos(beta), 0.0])
    axial = np.array([0.0, 0.0, 1.0])
    return radial, azimuthal, axial


def _narrowed_radial(msh, facet_tags, tag, fraction, centre, radial, half_width):
    fdim = msh.topology.dim - 1
    idx = np.asarray(facet_tags.indices, dtype=np.int32)
    val = np.asarray(facet_tags.values, dtype=np.int32)
    on_sheet = val == int(tag)
    sheet_facets = idx[on_sheet]
    keep = np.ones(sheet_facets.size, dtype=bool)
    if sheet_facets.size:
        mid = dolfinx.mesh.compute_midpoints(msh, fdim, sheet_facets)
        keep = np.abs(mid @ np.asarray(radial, dtype=float) - centre) <= (fraction * half_width)
    new_idx = np.concatenate([idx[~on_sheet], sheet_facets[keep]]).astype(np.int32)
    new_val = np.concatenate(
        [
            val[~on_sheet],
            np.full(int(np.count_nonzero(keep)), int(tag), dtype=np.int32),
        ]
    ).astype(np.int32)
    order = np.argsort(new_idx, kind="stable")
    return dolfinx.mesh.meshtags(msh, fdim, new_idx[order], new_val[order])
