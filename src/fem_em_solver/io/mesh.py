"""Mesh generation utilities for EM simulations."""

import math
import time
from typing import Dict, List, Optional, Sequence, Tuple

import dolfinx
import gmsh
import numpy as np
from dolfinx.io import gmsh as dolfinx_gmsh
from mpi4py import MPI


def _model_to_mesh(model, comm, rank, **kwargs):
    """`OPS-18` step 2: the 0.7 → 0.11 gmsh-interop shim.

    ``dolfinx.io.gmshio.model_to_mesh`` moved to ``dolfinx.io.gmsh`` and now
    returns a six-field :class:`~dolfinx.io.gmsh.MeshData` named tuple
    (``mesh, cell_tags, facet_tags, ridge_tags, peak_tags, physical_groups``)
    instead of the 0.7 three-tuple, so every call site's unpacking broke.  This
    module is the project's only gmsh-interop boundary; keeping the conversion
    here means one definition to revisit at the next upgrade, and the call
    sites keep the ``(mesh, cell_tags, facet_tags)`` contract the rest of
    ``io/mesh.py`` is written against.  ``ridge_tags``/``peak_tags`` are 0.11
    additions this project has no consumer for; ``physical_groups`` duplicates
    tag bookkeeping the generators already track by name.
    """
    data = dolfinx_gmsh.model_to_mesh(model, comm, rank, **kwargs)
    return data.mesh, data.cell_tags, data.facet_tags


def _raise_geometry_failure_on_every_rank(
    comm: MPI.Intracomm,
    rank: int,
    generator: str,
    build_error: Optional[BaseException],
    resolution: Optional[float] = None,
) -> None:
    """`GEO-23` step 2a: turn a rank-local gmsh throw into a collective raise.

    Only the building rank runs the gmsh model construction, so a throw there
    (``Invalid boundary mesh (overlapping facets)`` is this project's recurring
    one) leaves every other rank blocked forever inside the collective
    ``_model_to_mesh``: the command dies on the harness timeout instead of
    reporting the geometry failure.  Step 1 measured the price — the three
    unwrapped sites cost Status 124 at 120–121 s each where the already-wrapped
    ``birdcage_port_domain`` footers at Status 1 in 5 s.

    The pattern (originally ``birdcage_port_domain``, `GEO-9` step 2b): catch
    ``BaseException`` on the building rank, ``gmsh.finalize()`` so the failure
    does not poison later ``occ`` calls in the same process (known-issues 7),
    then ``bcast`` the flag and re-raise on *every* rank before the collective.
    The building rank re-raises the original exception; the others raise a
    ``RuntimeError`` naming the generator, the building rank and the
    ``resolution`` that failed.
    """
    if comm.bcast(build_error is not None, root=rank):
        if build_error is not None:
            raise build_error
        raise RuntimeError(
            f"{generator} geometry generation failed on rank {rank} "
            f"(resolution={resolution!r}); this is rank {comm.rank}"
        )


#: `GEO-13`: ``cylindrical_domain``'s surface-classification tolerance, as a
#: fraction of the radial gap ``outer_radius - inner_radius``.  Set from
#: `20260807T033127Z_GEO-13-probe.log`; see the sizing note at its use site.
_WALL_TOL_FRACTION = 0.01


def _cylindrical_hole_domain(
    a, c, length, resolution_wall, resolution_far, comm, rank, coat_radius=None
):
    """`IMP-4` step 1 (§9 item 206): ``cylindrical_domain(as_hole=True)``.

    The radius-``a`` cylinder is *cut* from the radius-``c`` one (both of
    ``length`` centred on z = 0): volume tag ``2`` is the only cell group.
    Facet groups come from ``getBoundary`` of the **meshed** volume (`TH-15`'s
    rule), classified z first — a face with z-extent < tol at |z| = L/2 is an
    end annulus (tag ``3``; the bbox r_max of an annulus is c, so an r-first
    rule would put it in tag 1); then bbox r_max within tol of c is the outer
    lateral wall (tag ``1``); the rest is the cavity wall (tag ``2``).  Each
    group must be non-empty.  tol = ``_WALL_TOL_FRACTION``·(c − a).

    Sizing: Distance from the cavity faces → Threshold (``resolution_wall`` at
    distance 0, ``resolution_far`` at distance c − a) as the background mesh,
    i.e. h linear in ρ.  No Netgen optimise on this branch.

    ``coat_radius`` (`IMP-5` step 1, §9 item 211; a < b < c): the sleeve
    ``cut(cyl_b, cyl_a)`` and the host ``cut(cyl_c, cyl_b)`` are
    ``occ.fragment``-ed so the ρ = b face is shared (conforming); the host
    keeps volume tag ``2`` and the sleeve is cell tag ``4`` "coating"
    (identified by bbox r_max, b against c).  Facet groups come from
    ``getBoundary`` of both volumes combined, so the shared face drops out;
    no interface facet group is emitted.  Sizing is unchanged.  ``None``
    takes the single-volume branch untouched.
    """
    if not 0.0 < a < c:
        raise ValueError(f"need 0 < inner_radius < outer_radius, got {a!r}, {c!r}")
    if length <= 0.0 or resolution_wall <= 0.0 or resolution_far <= 0.0:
        raise ValueError("length and resolutions must be positive")
    if coat_radius is not None and not a < coat_radius < c:
        raise ValueError(
            f"need inner_radius < coat_radius < outer_radius, got {a!r}, " f"{coat_radius!r}, {c!r}"
        )
    build_error: Optional[BaseException] = None
    if comm.rank == rank:
        try:
            gmsh.initialize()
            gmsh.model.add("cylindrical_hole_domain")
            if coat_radius is None:
                inner = gmsh.model.occ.addCylinder(0, 0, -length / 2, 0, 0, length, a)
                outer = gmsh.model.occ.addCylinder(0, 0, -length / 2, 0, 0, length, c)
                cut_out, _ = gmsh.model.occ.cut([(3, outer)], [(3, inner)])
                gmsh.model.occ.synchronize()
                volumes = gmsh.model.getEntities(dim=3)
                if len(volumes) != 1:
                    raise RuntimeError(
                        "cylindrical_domain(as_hole=True): cut did not leave a single "
                        f"volume (got {volumes}, cut returned {cut_out})"
                    )
                vol = volumes[0][1]
                gmsh.model.addPhysicalGroup(3, [vol], tag=2)
                gmsh.model.setPhysicalName(3, 2, "host")
                boundary_of = [(3, vol)]
            else:
                b = float(coat_radius)
                occ = gmsh.model.occ
                cyl_a = occ.addCylinder(0, 0, -length / 2, 0, 0, length, a)
                cyl_b = occ.addCylinder(0, 0, -length / 2, 0, 0, length, b)
                sleeve, _ = occ.cut([(3, cyl_b)], [(3, cyl_a)], removeTool=True)
                cyl_b2 = occ.addCylinder(0, 0, -length / 2, 0, 0, length, b)
                cyl_c = occ.addCylinder(0, 0, -length / 2, 0, 0, length, c)
                host, _ = occ.cut([(3, cyl_c)], [(3, cyl_b2)], removeTool=True)
                occ.fragment(sleeve, host)
                occ.synchronize()
                volumes = gmsh.model.getEntities(dim=3)
                if len(volumes) != 2:
                    raise RuntimeError(
                        "cylindrical_domain(as_hole=True, coat_radius): fragment did "
                        f"not leave exactly two volumes (got {volumes})"
                    )
                tol_v = _WALL_TOL_FRACTION * (c - a)
                host_vols, coat_vols = [], []
                for _, v in volumes:
                    x0, y0, _z0, x1, y1, _z1 = gmsh.model.getBoundingBox(3, v)
                    r_max = math.sqrt(max(x0**2, x1**2, y0**2, y1**2))
                    if abs(r_max - c) < tol_v:
                        host_vols.append(v)
                    elif abs(r_max - b) < tol_v:
                        coat_vols.append(v)
                if len(host_vols) != 1 or len(coat_vols) != 1:
                    raise RuntimeError(
                        "cylindrical_domain(as_hole=True, coat_radius): could not "
                        f"identify host/coating volumes (host {host_vols}, coat "
                        f"{coat_vols})"
                    )
                gmsh.model.addPhysicalGroup(3, host_vols, tag=2)
                gmsh.model.setPhysicalName(3, 2, "host")
                gmsh.model.addPhysicalGroup(3, coat_vols, tag=4)
                gmsh.model.setPhysicalName(3, 4, "coating")
                boundary_of = [(3, host_vols[0]), (3, coat_vols[0])]

            tol = _WALL_TOL_FRACTION * (c - a)
            walls, cavity, ends = [], [], []
            for dim, surf in gmsh.model.getBoundary(
                boundary_of, combined=True, oriented=False, recursive=False
            ):
                x0, y0, z0, x1, y1, z1 = gmsh.model.getBoundingBox(dim, surf)
                if (z1 - z0) < tol and abs(abs(0.5 * (z0 + z1)) - length / 2) < tol:
                    ends.append(surf)
                    continue
                r_max = math.sqrt(max(x0**2, x1**2, y0**2, y1**2))
                if abs(r_max - c) < tol:
                    walls.append(surf)
                else:
                    cavity.append(surf)
            print(
                f"[cylindrical_domain(as_hole)] faces: outer wall {len(walls)}, "
                f"cavity {len(cavity)}, ends {len(ends)}",
                flush=True,
            )
            if not walls or not cavity or not ends:
                raise RuntimeError(
                    "cylindrical_domain(as_hole=True): empty facet group (walls "
                    f"{walls}, cavity {cavity}, ends {ends})"
                )
            gmsh.model.addPhysicalGroup(2, walls, tag=1)
            gmsh.model.setPhysicalName(2, 1, "outer_boundary")
            gmsh.model.addPhysicalGroup(2, cavity, tag=2)
            gmsh.model.setPhysicalName(2, 2, "cavity_surface")
            gmsh.model.addPhysicalGroup(2, ends, tag=3)
            gmsh.model.setPhysicalName(2, 3, "end_planes")

            dist = gmsh.model.mesh.field.add("Distance")
            gmsh.model.mesh.field.setNumbers(dist, "SurfacesList", cavity)
            gmsh.model.mesh.field.setNumber(dist, "Sampling", 200)
            thr = gmsh.model.mesh.field.add("Threshold")
            gmsh.model.mesh.field.setNumber(thr, "InField", dist)
            gmsh.model.mesh.field.setNumber(thr, "SizeMin", resolution_wall)
            gmsh.model.mesh.field.setNumber(thr, "SizeMax", resolution_far)
            gmsh.model.mesh.field.setNumber(thr, "DistMin", 0.0)
            gmsh.model.mesh.field.setNumber(thr, "DistMax", c - a)
            gmsh.model.mesh.field.setAsBackgroundMesh(thr)
            gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
            gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
            gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
            gmsh.model.mesh.generate(3)

            etypes, etags, _ = gmsh.model.mesh.getElements(3)
            qmin = float("inf")
            for tags in etags:
                q = gmsh.model.mesh.getElementQualities(tags, "minSICN")
                if len(q):
                    qmin = min(qmin, float(np.min(q)))
            print(
                f"[cylindrical_domain(as_hole)] minimum SICN quality {qmin:.4f}",
                flush=True,
            )
            if coat_radius is not None:
                per_tag = {}
                for name, vols_ in (("host(2)", host_vols), ("coating(4)", coat_vols)):
                    n = 0
                    for v in vols_:
                        _, et, _ = gmsh.model.mesh.getElements(3, v)
                        n += sum(len(t) for t in et)
                    per_tag[name] = n
                print(
                    f"[cylindrical_domain(as_hole, coat)] cells per volume tag {per_tag}",
                    flush=True,
                )
        except BaseException as exc:  # noqa: BLE001 - re-raised on every rank
            build_error = exc
            if gmsh.isInitialized():
                gmsh.finalize()
    _raise_geometry_failure_on_every_rank(
        comm, rank, "cylindrical_domain(as_hole=True)", build_error, resolution_wall
    )
    mesh, cell_tags, facet_tags = _model_to_mesh(gmsh.model, comm, rank, gdim=3)
    if comm.rank == rank:
        gmsh.finalize()
    return mesh, cell_tags, facet_tags


#: `PORT-1` step 3b-xvi: ceiling [m] on the ``Thickness`` of
#: ``two_torus_domain``'s gap-box refinement field, i.e. how far outside each
#: gap box the refinement is allowed to ramp back to ``far_resolution``.  See
#: the sizing note at its use site — the slope-0.3 rule the other size fields
#: use puts a 10 cm shell of refinement in the air here, and the locality
#: control measured it at +35.4560% of the outside cell count.
GAP_BOX_THICKNESS_CAP_M = 5.0e-3

#: `TH-15` step 2a: the dim-2 physical-group tag ``two_torus_domain(
#: as_hole=True)`` puts on the conductor cavity wall — the surface a PEC
#: boundary condition is imposed on when the conductor is a hole rather than a
#: meshed solid.  It is an *exterior* facet group (the cavity is outside the
#: mesh), unlike the interior port/sheet groups ``201``/``202``/``211``/``212``,
#: so it is emitted as a gmsh physical group directly and not rebuilt from cell
#: tags.  Value carried over from `TH-15` step 0's probe.
TWO_TORUS_CONDUCTOR_SURFACE_TAG = 301

#: `TH-15` step 3a: the same thing for ``birdcage_port_domain(as_hole=True)`` —
#: the dim-2 physical group on the coil's cavity wall once the ring and leg
#: solids have been cut out of the air.  Distinct value from the two-torus one
#: so a mesh that ever carried both could not confuse them; 401 is free of the
#: birdcage's cell groups (1/2/3, 100+i, 200+i) and of its only other facet
#: group (``1``, the outer boundary).
BIRDCAGE_CONDUCTOR_SURFACE_TAG = 401

#: `GEO-36` step 2: the dim-2 physical group on the wall of a capped rod cut
#: out of the phantom as a hole (``birdcage_port_domain(implant_rod=...)``).
#: Exterior facet group, like 401; 402 is free of the birdcage's cell groups
#: (1/2/3, 100+i, 200+i) and of its facet groups (1, 401, 501, 502).
BIRDCAGE_IMPLANT_SURFACE_TAG = 402

#: `GEO-36` step 3 (§9 item 219): the cell group of the rod's own cells when
#: ``implant_rod={..., "keep_cells": True}`` keeps them (the rod is then a
#: material region, and 402 is the *interior* phantom/implant interface).
#: 403 is free of the birdcage's cell groups (1/2/3, 100+i, 200+i; i <= 32 on
#: F-human) and of its facet groups (1, 401, 402, 501, 502) — grep of
#: ``src/`` for ``403`` was empty on 2026-10-04.
BIRDCAGE_IMPLANT_CELL_TAG = 403

#: `GEO-36` step 5 (§9 item 316): the cell group of the shaft-only dielectric
#: sleeve when ``implant_rod={..., "keep_cells": True, "coat_radius": b}``
#: (the cylinder of radius b between the two cap-centre planes, minus the
#: rod; the hemispherical caps stay bare).  No new facet tag: 402 stays the
#: (3, 403) interface, which then covers the bare caps only.  A grep of
#: ``src/`` for ``404`` was empty on 2026-10-08.
BIRDCAGE_IMPLANT_COAT_CELL_TAG = 404

#: `TH-16` step 3a: the two symmetry cut planes of
#: ``birdcage_port_domain(symmetry_quadrant=True)`` — dim-2 physical groups on
#: the x = 0 and y = 0 faces of the quarter air box.  Kept out of the outer
#: boundary group (``1``) so a later step can assign PEC / PMC per plane.
BIRDCAGE_CUT_X0_TAG = 501
BIRDCAGE_CUT_Y0_TAG = 502


def _interface_facet_tags(
    mesh: "dolfinx.mesh.Mesh",
    cell_tags: "dolfinx.mesh.MeshTags",
    interfaces: Dict[int, Tuple[int, int]],
    existing: Optional["dolfinx.mesh.MeshTags"] = None,
) -> "dolfinx.mesh.MeshTags":
    """Tag the facets shared by two differently-tagged cell regions.

    ``interfaces`` maps a facet tag to the pair of *cell* tags whose common
    facets carry it (order irrelevant), or to a *sequence* of such pairs when
    one facet group is the union of several cell-tag interfaces (`GEO-16`: the
    gap box split at its port-sheet mid-plane meets the conductor with both
    halves, so port group ``201`` is ``(101, 1)`` and ``(111, 1)`` together).
    The result merges those facets into
    ``existing``'s tags when one is given.

    Why this is done here and not in gmsh: dim-2 physical groups on facets
    that are **interior** to the partitioned mesh hang ``model_to_mesh`` at
    ``-n 2`` inside ``distribute_entity_data`` (known-issues 9, measured
    2026-08-05; ``-n 1`` completes the identical case in 22.5 s). Cell tags
    distribute fine, and an interface is exactly derivable from them, so the
    interface is reconstructed on the dolfinx side where every rank already
    holds the data it needs.

    Rank-safety: the classification needs the tag of *both* cells behind a
    facet, and one of them is a ghost on a partition-boundary facet. Ghost
    cells are not necessarily carried by ``cell_tags``, so the tag is pushed
    through a DG0 function and ``scatter_forward``ed rather than read from
    ``cell_tags.values`` directly.
    """
    tdim = mesh.topology.dim
    fdim = tdim - 1
    mesh.topology.create_entities(fdim)
    mesh.topology.create_connectivity(fdim, tdim)

    marker_space = dolfinx.fem.functionspace(mesh, ("DG", 0))
    marker = dolfinx.fem.Function(marker_space)
    marker.x.array[:] = 0.0
    cell_to_dof = marker_space.dofmap.list.reshape(-1)
    marker.x.array[cell_to_dof[cell_tags.indices]] = cell_tags.values
    marker.x.scatter_forward()
    cell_value = np.rint(np.real(marker.x.array[cell_to_dof])).astype(np.int32)

    facet_to_cell = mesh.topology.connectivity(fdim, tdim)
    offsets = facet_to_cell.offsets
    links = facet_to_cell.array
    counts = offsets[1:] - offsets[:-1]
    interior = np.flatnonzero(counts == 2)
    side_a = cell_value[links[offsets[interior]]]
    side_b = cell_value[links[offsets[interior] + 1]]

    indices: List[np.ndarray] = []
    values: List[np.ndarray] = []
    if existing is not None and existing.indices.size:
        indices.append(np.asarray(existing.indices, dtype=np.int32))
        values.append(np.asarray(existing.values, dtype=np.int32))

    for facet_tag, spec in sorted(interfaces.items()):
        pairs = [spec] if len(spec) == 2 and np.isscalar(spec[0]) else list(spec)
        on_interface = np.zeros(interior.size, dtype=bool)
        for tag_a, tag_b in pairs:
            on_interface |= ((side_a == tag_a) & (side_b == tag_b)) | (
                (side_a == tag_b) & (side_b == tag_a)
            )
        found = interior[on_interface]
        indices.append(found.astype(np.int32))
        values.append(np.full(found.size, facet_tag, dtype=np.int32))

    if indices:
        all_indices = np.concatenate(indices)
        all_values = np.concatenate(values)
    else:
        all_indices = np.empty(0, dtype=np.int32)
        all_values = np.empty(0, dtype=np.int32)

    order = np.argsort(all_indices, kind="stable")
    all_indices = all_indices[order]
    all_values = all_values[order]
    if np.unique(all_indices).size != all_indices.size:
        raise RuntimeError(
            "_interface_facet_tags: a facet was claimed by more than one tag; "
            f"{all_indices.size} entries, {np.unique(all_indices).size} unique"
        )
    return dolfinx.mesh.meshtags(mesh, fdim, all_indices, all_values)


def _z_rotation_affine(angle: float) -> Tuple[float, ...]:
    """Row-major 4x4 rotation about `ẑ`, snapped exact on the coordinate axes.

    `GEO-19` step B places the birdcage port boxes and their sheets by rotating
    an azimuth-0 construction, so at a multiple of 90 degrees the transform must
    be the *exact* identity or the exact coordinate swap — `occ.rotate` would
    instead apply `cos(pi/2) = 6.1e-17` and move every vertex by a few ulps.
    That is not cosmetic: the ulp is below any geometric tolerance, but gmsh's
    tie-breaking is not tolerant, and last-bit coordinates change the cell count
    at the 1e-3 level (measured, `GEO-19` step B, 2026-08-24).
    """
    c, s = float(np.cos(angle)), float(np.sin(angle))
    for name, value in (("c", c), ("s", s)):
        for exact in (-1.0, 0.0, 1.0):
            if abs(value - exact) < 1.0e-12:
                if name == "c":
                    c = exact
                else:
                    s = exact
    return (c, -s, 0.0, 0.0, s, c, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0)


def _rotate_about_z(dim: int, tag: int, angle: float) -> None:
    """Rotate one entity about `ẑ` through the origin; a no-op at azimuth 0."""
    matrix = _z_rotation_affine(angle)
    if matrix[0] == 1.0 and matrix[4] == 0.0:
        return
    gmsh.model.occ.affineTransform([(dim, tag)], list(matrix))


def _place_sheet_in_leg_frame(tag: int, angle: float) -> None:
    """Take an xy-plane rectangle to the leg's mid-plane at `angle`.

    One transform, `R_z(angle) . R_x(pi/2)`, so no intermediate rounding: the
    in-plane extent that must become the gap goes to `ẑ` and the `y`-normal
    plane of a leg on `+x` turns to the leg's own azimuth. Snapped by
    `_z_rotation_affine`, so at 0 / 90 / 180 / 270 degrees every entry is
    exactly 0 or +-1 and the rectangle's corners stay exactly on their
    coordinate values.
    """
    c, _ms, _z0, _z1, s, *_ = _z_rotation_affine(angle)
    gmsh.model.occ.affineTransform(
        [(2, tag)],
        [
            c,
            0.0,
            s,
            0.0,
            s,
            0.0,
            -c,
            0.0,
            0.0,
            1.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
        ],
    )


# `GEO-4` step 2 (§9 item 185): the box-term padding heuristic for the
# two-torus fixture.  Provenance: `PORT-1` step 3b-xi measured the signed
# ``Im Z₁₂/ωM₁₂ − 1`` deficit (percentage points, projected drive, h_far 0.03,
# d = 0.04) at three air paddings (`20260807T124038Z_PORT-1-step3bxi-gate.log`);
# `PORT-1` decision (4) fitted deficit(W) = D∞ + C·W^(−p) through them
# (`scripts/probes/port1_padding_fit.py`, `20260812T213337Z_PORT-1-dec4-fit.log`).
# The fit is recomputed from these rungs on every call — never hard-code the
# printed C / p (rounded, they reproduce the 0.08 m rung only to 5e-4 pp).
PORT1_PADDING_RUNGS_M = (0.08, 0.10, 0.12)
PORT1_PADDING_DEFICITS_PP = (-8.0324, -5.0256, -3.2733)


def _free_exponent_fit(widths, values):
    """Fit ``value(W) = v_inf − C·W**(−p)`` through exactly three points.

    Copied verbatim from ``scripts/probes/port1_padding_fit.py``
    ``power_law_extrapolation`` (`PORT-1` decision (4)), itself carried over
    from ``tests/validation/test_dodd_deeds_reactance_box_truncation.py``
    (`MAT-6` step 9) unchanged in method: writing ``e_i = v_inf − value_i`` the
    model is linear in logs — ``ln e_i = ln C − p ln W_i`` — so the correct
    ``v_inf`` is the one making the three ``(ln W_i, ln e_i)`` collinear, i.e.
    the one where the two adjacent-pair slopes agree.  Bisect on that residual.

    Note the sign convention: ``v_inf − C·W**(−p)`` here is the commissioning
    text's ``D_inf + C·W**(−p)`` with ``C`` negated.

    Returns ``(v_inf, p)``, or ``(None, None)`` when the trend does not admit
    such a fit (non-monotone rungs, or no sign change in the bracket).
    """
    w = [float(x) for x in widths]
    r = [float(x) for x in values]
    if not (r[0] < r[1] < r[2]):
        return None, None

    def slope_mismatch(v_inf):
        e = [v_inf - x for x in r]
        if min(e) <= 0.0:
            return float("nan")
        p_lo = (math.log(e[0]) - math.log(e[1])) / (math.log(w[1]) - math.log(w[0]))
        p_hi = (math.log(e[1]) - math.log(e[2])) / (math.log(w[2]) - math.log(w[1]))
        return p_lo - p_hi

    lo = r[2] + 1.0e-12
    hi = r[2] + 10.0 * (abs(r[2]) + 1.0)
    f_lo, f_hi = slope_mismatch(lo), slope_mismatch(hi)
    if not (math.isfinite(f_lo) and math.isfinite(f_hi)) or f_lo * f_hi > 0.0:
        return None, None
    for _ in range(400):
        mid = 0.5 * (lo + hi)
        f_mid = slope_mismatch(mid)
        if not math.isfinite(f_mid):
            return None, None
        if f_lo * f_mid <= 0.0:
            hi, f_hi = mid, f_mid
        else:
            lo, f_lo = mid, f_mid
    v_inf = 0.5 * (lo + hi)
    e = [v_inf - x for x in r]
    p = (math.log(e[0]) - math.log(e[2])) / (math.log(w[2]) - math.log(w[0]))
    return v_inf, p


def _padding_for_box_term(c_pp: float, p: float, tolerance_pp: float, quantum_m: float):
    """Smallest padding, on a ``quantum_m`` grid, with ``|C|·W^(−p) ≤ τ``.

    Returns ``(padding_m, unrounded_padding_m)``.  The 1e-12 guard keeps an
    unrounded value that lands on a grid point (to round-off) from being
    pushed up one quantum.
    """
    if tolerance_pp <= 0.0:
        raise ValueError("tolerance_pp must be > 0")
    if quantum_m <= 0.0:
        raise ValueError("quantum_m must be > 0")
    if not (p > 0.0):
        raise ValueError("box-term exponent must be > 0")
    unrounded = (abs(c_pp) / tolerance_pp) ** (1.0 / p)
    steps = math.ceil(unrounded / quantum_m - 1.0e-12)
    return round(steps * quantum_m, 12), unrounded


class MeshGenerator:
    """Generate meshes for common geometries using Gmsh."""

    @staticmethod
    def straight_wire_domain(
        wire_length: float = 1.0,
        wire_radius: float = 0.001,
        domain_radius: float = 0.1,
        resolution: float = 0.005,
        comm: MPI.Intracomm = MPI.COMM_WORLD,
        rank: int = 0,
    ) -> Tuple[dolfinx.mesh.Mesh, dolfinx.mesh.MeshTags, dolfinx.mesh.MeshTags]:
        """Generate mesh for straight wire in cylindrical domain.

        Creates a cylindrical domain with a thin wire along the z-axis.
        Useful for validating against the analytical solution for an
        infinite straight wire.

        Parameters
        ----------
        wire_length : float
            Length of wire [m]
        wire_radius : float
            Radius of wire [m] (should be small for thin wire approximation)
        domain_radius : float
            Radius of surrounding cylindrical domain [m]
        resolution : float
            Characteristic mesh size [m]
        comm : MPI.Intracomm
            MPI communicator
        rank : int
            Rank for Gmsh model (usually 0)

        Returns
        -------
        mesh : dolfinx.mesh.Mesh
            The generated mesh
        cell_tags : dolfinx.mesh.MeshTags
            Cell tags for subdomains (wire, surrounding)
        facet_tags : dolfinx.mesh.MeshTags
            Facet tags for boundaries
        """
        build_error: Optional[BaseException] = None
        if comm.rank == rank:
            try:
                # Initialize Gmsh
                gmsh.initialize()
                gmsh.model.add("straight_wire")

                # Create wire (cylinder along z-axis)
                wire_tag = gmsh.model.occ.addCylinder(
                    0,
                    0,
                    -wire_length / 2,  # center of bottom face
                    0,
                    0,
                    wire_length,  # axis direction and height
                    wire_radius,
                )

                # Create surrounding domain (hollow cylinder)
                domain_tag = gmsh.model.occ.addCylinder(
                    0, 0, -wire_length / 2, 0, 0, wire_length, domain_radius
                )

                # Cut wire out of domain to create separate volumes
                # Actually, we want both as separate volumes, so we fragment
                ov, ovv = gmsh.model.occ.fragment([(3, domain_tag)], [(3, wire_tag)])
                gmsh.model.occ.synchronize()

                # Get volumes
                volumes = gmsh.model.getEntities(dim=3)
                wire_volume = None
                domain_volume = None

                # Tag volumes based on their bounding box
                for vol in volumes:
                    bbox = gmsh.model.getBoundingBox(vol[0], vol[1])
                    # Check if this is the wire (small radius)
                    x_min, y_min, z_min, x_max, y_max, z_max = bbox
                    r_max = np.sqrt(max(x_max**2, y_max**2))
                    if r_max < 2 * wire_radius:
                        wire_volume = vol[1]
                    else:
                        domain_volume = vol[1]

                # Add physical groups
                if wire_volume:
                    gmsh.model.addPhysicalGroup(3, [wire_volume], tag=1)
                    gmsh.model.setPhysicalName(3, 1, "wire")

                if domain_volume:
                    gmsh.model.addPhysicalGroup(3, [domain_volume], tag=2)
                    gmsh.model.setPhysicalName(3, 2, "domain")

                # Tag boundaries
                surfaces = gmsh.model.getEntities(dim=2)
                boundary_surfaces = []
                wire_surfaces = []

                for surf in surfaces:
                    bbox = gmsh.model.getBoundingBox(surf[0], surf[1])
                    x_min, y_min, z_min, x_max, y_max, z_max = bbox
                    r_max = np.sqrt(max(x_max**2, y_max**2))

                    # Cylindrical boundary of domain
                    if abs(r_max - domain_radius) < resolution:
                        boundary_surfaces.append(surf[1])
                    # Wire surface
                    elif r_max < 2 * wire_radius:
                        wire_surfaces.append(surf[1])

                if boundary_surfaces:
                    gmsh.model.addPhysicalGroup(2, boundary_surfaces, tag=1)
                    gmsh.model.setPhysicalName(2, 1, "boundary")

                if wire_surfaces:
                    gmsh.model.addPhysicalGroup(2, wire_surfaces, tag=2)
                    gmsh.model.setPhysicalName(2, 2, "wire_surface")

                # Set mesh size
                gmsh.model.mesh.setSize(gmsh.model.getEntities(0), resolution)

                # Generate mesh
                gmsh.model.mesh.generate(3)
                gmsh.model.mesh.optimize("Netgen")
            except BaseException as exc:  # noqa: BLE001 - re-raised on every rank
                build_error = exc
                if gmsh.isInitialized():
                    gmsh.finalize()

        # Collective: every rank learns whether the geometry exists before
        # any of them enters the collective `_model_to_mesh` below (GEO-23
        # step 2a; the unwrapped throw deadlocked at 120 s instead of
        # failing in seconds).
        _raise_geometry_failure_on_every_rank(
            comm, rank, "straight_wire_domain", build_error, resolution
        )

        # Convert to dolfinx mesh
        mesh, cell_tags, facet_tags = _model_to_mesh(gmsh.model, comm, rank, gdim=3)

        if comm.rank == rank:
            gmsh.finalize()

        return mesh, cell_tags, facet_tags

    @staticmethod
    def circular_loop_domain(
        loop_radius: float = 0.05,
        wire_radius: float = 0.001,
        domain_radius: float = 0.15,
        resolution: float = 0.005,
        comm: MPI.Intracomm = MPI.COMM_WORLD,
        rank: int = 0,
    ) -> Tuple[dolfinx.mesh.Mesh, dolfinx.mesh.MeshTags, dolfinx.mesh.MeshTags]:
        """Generate mesh for circular current loop in spherical domain.

        Creates a torus (ring) for the wire surrounded by a spherical air domain.
        The loop lies in the xy-plane centered at origin.

        Parameters
        ----------
        loop_radius : float
            Major radius of loop (distance from center to wire center) [m]
        wire_radius : float
            Minor radius of wire cross-section [m]
        domain_radius : float
            Radius of surrounding spherical domain [m]
        resolution : float
            Characteristic mesh size [m]
        comm : MPI.Intracomm
            MPI communicator
        rank : int
            Rank for Gmsh model (usually 0)

        Returns
        -------
        mesh : dolfinx.mesh.Mesh
            The generated mesh
        cell_tags : dolfinx.mesh.MeshTags
            Cell tags for subdomains (wire=1, air=2)
        facet_tags : dolfinx.mesh.MeshTags
            Facet tags for boundaries
        """
        if comm.rank == rank:
            # Initialize Gmsh
            gmsh.initialize()
            gmsh.model.add("circular_loop")

            # Create wire as torus in xy-plane
            # addTorus(x, y, z, major_radius, minor_radius)
            wire_tag = gmsh.model.occ.addTorus(0, 0, 0, loop_radius, wire_radius)

            # Create surrounding spherical domain
            domain_tag = gmsh.model.occ.addSphere(0, 0, 0, domain_radius)

            # Fragment to get separate volumes
            ov, ovv = gmsh.model.occ.fragment([(3, domain_tag)], [(3, wire_tag)])
            gmsh.model.occ.synchronize()

            # Get volumes and tag them
            volumes = gmsh.model.getEntities(dim=3)
            wire_volume = None
            air_volume = None

            for vol in volumes:
                bbox = gmsh.model.getBoundingBox(vol[0], vol[1])
                x_min, y_min, z_min, x_max, y_max, z_max = bbox

                # Compute bounding box dimensions
                dx = x_max - x_min
                dy = y_max - y_min
                dz = z_max - z_min

                # Wire has smaller bounding box
                max_dim = max(dx, dy, dz)
                if max_dim < 4 * loop_radius:  # Wire is smaller
                    wire_volume = vol[1]
                else:
                    air_volume = vol[1]

            # Add physical groups
            if wire_volume:
                gmsh.model.addPhysicalGroup(3, [wire_volume], tag=1)
                gmsh.model.setPhysicalName(3, 1, "wire")

            if air_volume:
                gmsh.model.addPhysicalGroup(3, [air_volume], tag=2)
                gmsh.model.setPhysicalName(3, 2, "air")

            # Tag boundaries
            surfaces = gmsh.model.getEntities(dim=2)
            boundary_surfaces = []

            for surf in surfaces:
                bbox = gmsh.model.getBoundingBox(surf[0], surf[1])
                x_min, y_min, z_min, x_max, y_max, z_max = bbox

                # Check if on outer spherical boundary
                r_max = np.sqrt(x_max**2 + y_max**2 + z_max**2)
                if abs(r_max - domain_radius) < resolution:
                    boundary_surfaces.append(surf[1])

            if boundary_surfaces:
                gmsh.model.addPhysicalGroup(2, boundary_surfaces, tag=1)
                gmsh.model.setPhysicalName(2, 1, "outer_boundary")

            # Set mesh size
            gmsh.model.mesh.setSize(gmsh.model.getEntities(0), resolution)

            # Generate mesh
            gmsh.model.mesh.generate(3)
            gmsh.model.mesh.optimize("Netgen")

        # Convert to dolfinx mesh
        mesh, cell_tags, facet_tags = _model_to_mesh(gmsh.model, comm, rank, gdim=3)

        if comm.rank == rank:
            gmsh.finalize()

        return mesh, cell_tags, facet_tags

    @staticmethod
    def helmholtz_coil_domain(
        loop_radius: float = 0.05,
        wire_radius: float = 0.002,  # Increased from 0.001 for simpler mesh
        domain_radius: float = 0.12,  # Reduced from 0.15
        resolution: float = 0.008,  # Coarser mesh
        comm: MPI.Intracomm = MPI.COMM_WORLD,
        rank: int = 0,
    ) -> Tuple[dolfinx.mesh.Mesh, dolfinx.mesh.MeshTags, dolfinx.mesh.MeshTags]:
        """Generate mesh for Helmholtz coil in spherical domain.

        A Helmholtz coil consists of two identical circular loops separated
        by a distance equal to their radius. This configuration creates a
        highly uniform magnetic field in the central region.

        The loops are positioned at z = -R/2 and z = +R/2, where R is the
        loop radius.

        Parameters
        ----------
        loop_radius : float
            Major radius of each loop (distance from center to wire center) [m]
        wire_radius : float
            Minor radius of wire cross-section [m]
        domain_radius : float
            Radius of surrounding spherical domain [m]
        resolution : float
            Characteristic mesh size [m]
        comm : MPI.Intracomm
            MPI communicator
        rank : int
            Rank for Gmsh model (usually 0)

        Returns
        -------
        mesh : dolfinx.mesh.Mesh
            The generated mesh
        cell_tags : dolfinx.mesh.MeshTags
            Cell tags for subdomains (wire=1, air=2)
        facet_tags : dolfinx.mesh.MeshTags
            Facet tags for boundaries
        """
        if comm.rank == rank:
            # Initialize Gmsh
            gmsh.initialize()
            gmsh.model.add("helmholtz_coil")

            # Helmholtz condition: separation = loop radius
            separation = loop_radius
            z1 = -separation / 2
            z2 = separation / 2

            # Create two wire tori at z = -R/2 and z = +R/2
            # addTorus(x, y, z, major_radius, minor_radius)
            wire1_tag = gmsh.model.occ.addTorus(0, 0, z1, loop_radius, wire_radius)
            wire2_tag = gmsh.model.occ.addTorus(0, 0, z2, loop_radius, wire_radius)

            # Create surrounding spherical domain
            domain_tag = gmsh.model.occ.addSphere(0, 0, 0, domain_radius)

            # Fragment to get separate volumes
            # First fragment domain with wire1
            ov1, ovv1 = gmsh.model.occ.fragment([(3, domain_tag)], [(3, wire1_tag)])
            gmsh.model.occ.synchronize()

            # Then fragment result with wire2
            # Get the air volume from previous fragmentation
            volumes_after_first = gmsh.model.getEntities(dim=3)
            air_tag = None
            for vol in volumes_after_first:
                bbox = gmsh.model.getBoundingBox(vol[0], vol[1])
                x_min, y_min, z_min, x_max, y_max, z_max = bbox
                max_dim = max(x_max - x_min, y_max - y_min, z_max - z_min)
                # Wire is smaller than air domain
                if not max_dim < 4 * loop_radius:
                    air_tag = vol[1]

            # Fragment air with second wire
            if air_tag:
                ov2, ovv2 = gmsh.model.occ.fragment([(3, air_tag)], [(3, wire2_tag)])
            gmsh.model.occ.synchronize()

            # Get final volumes and tag them
            volumes = gmsh.model.getEntities(dim=3)
            wire_volumes = []
            air_volume = None

            for vol in volumes:
                bbox = gmsh.model.getBoundingBox(vol[0], vol[1])
                x_min, y_min, z_min, x_max, y_max, z_max = bbox

                # Compute bounding box dimensions
                dx = x_max - x_min
                dy = y_max - y_min
                dz = z_max - z_min
                max_dim = max(dx, dy, dz)

                # Wire volumes are smaller and have z-extent ~ 2*wire_radius
                if max_dim < 4 * loop_radius and dz < 4 * wire_radius:
                    wire_volumes.append(vol[1])
                else:
                    air_volume = vol[1]

            # Add physical groups
            if wire_volumes:
                # Tag both wires as "wire" (tag=1)
                gmsh.model.addPhysicalGroup(3, wire_volumes, tag=1)
                gmsh.model.setPhysicalName(3, 1, "wire")

            if air_volume:
                gmsh.model.addPhysicalGroup(3, [air_volume], tag=2)
                gmsh.model.setPhysicalName(3, 2, "air")

            # Tag outer boundary
            surfaces = gmsh.model.getEntities(dim=2)
            boundary_surfaces = []

            for surf in surfaces:
                bbox = gmsh.model.getBoundingBox(surf[0], surf[1])
                x_min, y_min, z_min, x_max, y_max, z_max = bbox

                # Check if on outer spherical boundary
                r_max = np.sqrt(x_max**2 + y_max**2 + z_max**2)
                if abs(r_max - domain_radius) < resolution:
                    boundary_surfaces.append(surf[1])

            if boundary_surfaces:
                gmsh.model.addPhysicalGroup(2, boundary_surfaces, tag=1)
                gmsh.model.setPhysicalName(2, 1, "outer_boundary")

            # Set mesh size
            gmsh.model.mesh.setSize(gmsh.model.getEntities(0), resolution)

            # Generate mesh
            gmsh.model.mesh.generate(3)
            # Skip optimization for Helmholtz coil to avoid Gmsh hangs
            # gmsh.model.mesh.optimize("Netgen")

        # Convert to dolfinx mesh
        mesh, cell_tags, facet_tags = _model_to_mesh(gmsh.model, comm, rank, gdim=3)

        if comm.rank == rank:
            gmsh.finalize()

        return mesh, cell_tags, facet_tags

    @staticmethod
    def rectangular_domain(
        bounds: Tuple[float, float, float, float, float, float],
        resolution: float,
        comm: MPI.Intracomm = MPI.COMM_WORLD,
        rank: int = 0,
    ) -> dolfinx.mesh.Mesh:
        """Generate simple rectangular domain.

        Parameters
        ----------
        bounds : tuple
            (xmin, xmax, ymin, ymax, zmin, zmax)
        resolution : float
            Mesh size
        comm : MPI.Intracomm
            MPI communicator
        rank : int
            Rank for Gmsh

        Returns
        -------
        mesh : dolfinx.mesh.Mesh
            The generated mesh
        """
        if comm.rank == rank:
            gmsh.initialize()
            gmsh.model.add("rectangular")

            xmin, xmax, ymin, ymax, zmin, zmax = bounds

            # Create box
            box = gmsh.model.occ.addBox(xmin, ymin, zmin, xmax - xmin, ymax - ymin, zmax - zmin)
            gmsh.model.occ.synchronize()

            gmsh.model.addPhysicalGroup(3, [box], tag=1)

            # Set mesh size
            gmsh.model.mesh.setSize(gmsh.model.getEntities(0), resolution)
            gmsh.model.mesh.generate(3)

        mesh, _, _ = _model_to_mesh(gmsh.model, comm, rank, gdim=3)

        if comm.rank == rank:
            gmsh.finalize()

        return mesh

    @staticmethod
    def create_simple_box(
        L: float = 1.0, n: int = 10, comm: MPI.Intracomm = MPI.COMM_WORLD
    ) -> dolfinx.mesh.Mesh:
        """Create simple box mesh using dolfinx built-in generator.

        Parameters
        ----------
        L : float
            Box half-length (domain is [-L, L]³)
        n : int
            Number of cells in each direction
        comm : MPI.Intracomm
            MPI communicator

        Returns
        -------
        mesh : dolfinx.mesh.Mesh
            Box mesh
        """
        from dolfinx.mesh import CellType, create_box

        domain = [(-L, -L, -L), (L, L, L)]
        mesh = create_box(comm, domain, [n, n, n], CellType.tetrahedron)

        return mesh

    @staticmethod
    def cylindrical_domain(
        inner_radius: float = 0.01,
        outer_radius: float = 0.1,
        length: float = 0.2,
        resolution: float = 0.02,
        comm: MPI.Intracomm = MPI.COMM_WORLD,
        rank: int = 0,
        *,
        as_hole: bool = False,
        resolution_wall: Optional[float] = None,
        resolution_far: Optional[float] = None,
        coat_radius: Optional[float] = None,
    ) -> Tuple[dolfinx.mesh.Mesh, dolfinx.mesh.MeshTags, dolfinx.mesh.MeshTags]:
        """Generate mesh with cylindrical inner volume inside cylindrical domain.

        ``as_hole=True`` (`IMP-4` step 1, §9 item 206) takes the separate
        `TH-15` hole route in ``_cylindrical_hole_domain``: the inner cylinder
        is cut, volume tag ``2`` only, facet tags ``1`` outer lateral wall,
        ``2`` cavity wall, ``3`` end annuli, sized h(ρ) = ``resolution_wall``
        at the cavity rising linearly to ``resolution_far`` at ρ = c.
        ``as_hole=False`` (the default) leaves the fragment branch below
        untouched; passing either resolution without ``as_hole`` raises.
        ``coat_radius`` (`IMP-5` step 1, keyword-only, ``as_hole=True`` only)
        adds a conforming sleeve inner_radius < ρ < coat_radius as cell tag
        ``4``; ``None`` leaves the hole branch untouched.
        """
        if as_hole:
            if resolution_wall is None or resolution_far is None:
                raise ValueError(
                    "cylindrical_domain(as_hole=True) needs resolution_wall and " "resolution_far"
                )
            return _cylindrical_hole_domain(
                float(inner_radius),
                float(outer_radius),
                float(length),
                float(resolution_wall),
                float(resolution_far),
                comm,
                rank,
                coat_radius=None if coat_radius is None else float(coat_radius),
            )
        if resolution_wall is not None or resolution_far is not None:
            raise ValueError("resolution_wall / resolution_far apply only with as_hole=True")
        if coat_radius is not None:
            raise ValueError("coat_radius applies only with as_hole=True")
        return MeshGenerator._cylindrical_fragment_domain(
            inner_radius, outer_radius, length, resolution, comm, rank
        )

    @staticmethod
    def _cylindrical_fragment_domain(
        inner_radius, outer_radius, length, resolution, comm, rank
    ) -> Tuple[dolfinx.mesh.Mesh, dolfinx.mesh.MeshTags, dolfinx.mesh.MeshTags]:
        """The pre-`IMP-4` ``cylindrical_domain`` body, unchanged.

        Creates two concentric cylinders along the z-axis. Useful for practicing
        multi-volume meshing and for problems with cylindrical symmetry.

        Parameters
        ----------
        inner_radius : float
            Radius of inner cylinder [m]
        outer_radius : float
            Radius of outer cylinder [m]
        length : float
            Length of cylinders along z-axis [m]
        resolution : float
            Characteristic mesh size [m]
        comm : MPI.Intracomm
            MPI communicator
        rank : int
            Rank for Gmsh model (usually 0)

        Returns
        -------
        mesh : dolfinx.mesh.Mesh
            The generated mesh
        cell_tags : dolfinx.mesh.MeshTags
            Cell tags for subdomains (inner=1, outer=2)
        facet_tags : dolfinx.mesh.MeshTags
            Facet tags for boundaries
        """
        build_error: Optional[BaseException] = None
        if comm.rank == rank:
            try:
                # Initialize Gmsh
                gmsh.initialize()
                gmsh.model.add("cylindrical_domain")

                # Create inner cylinder along z-axis
                inner_tag = gmsh.model.occ.addCylinder(
                    0,
                    0,
                    -length / 2,  # center of bottom face
                    0,
                    0,
                    length,  # axis direction and height
                    inner_radius,
                )

                # Create outer cylinder
                outer_tag = gmsh.model.occ.addCylinder(
                    0, 0, -length / 2, 0, 0, length, outer_radius
                )

                # Fragment to create separate volumes (inner and outer region)
                ov, ovv = gmsh.model.occ.fragment([(3, outer_tag)], [(3, inner_tag)])
                gmsh.model.occ.synchronize()

                # Get volumes and tag them
                volumes = gmsh.model.getEntities(dim=3)
                inner_volume = None
                outer_volume = None

                for vol in volumes:
                    bbox = gmsh.model.getBoundingBox(vol[0], vol[1])
                    x_min, y_min, z_min, x_max, y_max, z_max = bbox

                    # Check if this is the inner cylinder by radius
                    r_max = np.sqrt(max(x_max**2, y_max**2))
                    if r_max < (inner_radius + outer_radius) / 2:
                        inner_volume = vol[1]
                    else:
                        outer_volume = vol[1]

                # Add physical groups
                if inner_volume:
                    gmsh.model.addPhysicalGroup(3, [inner_volume], tag=1)
                    gmsh.model.setPhysicalName(3, 1, "inner")

                if outer_volume:
                    gmsh.model.addPhysicalGroup(3, [outer_volume], tag=2)
                    gmsh.model.setPhysicalName(3, 2, "outer")

                # Tag boundaries
                surfaces = gmsh.model.getEntities(dim=2)
                outer_boundary_surfaces = []
                inner_boundary_surfaces = []

                # `GEO-13`: the classification tolerance is a fraction of the radial
                # gap, never `resolution` — a tolerance keyed to mesh size makes the
                # margin geometry-over-mesh-size, and at `resolution >= 0.09` (the
                # gap itself) the inner cylinder was swept into `outer_boundary`
                # (known-issues 13, 6 of 6 surfaces accepted;
                # `20260807T033127Z_GEO-13-probe.log`).  0.01 sits in the middle of
                # the measured window [1e-4, 0.05] where both sides of the `GEO-11`
                # two-sided margin hold on every geometry this generator is called
                # with: worst accepted 1.1e-4 x tol (ceiling 0.1, the gmsh OCC
                # bounding-box padding is 1.000e-07) and nearest rejected 1.0e+02 x
                # tol (floor 10).  Precondition: the gap must exceed ~1e-4 m, or the
                # tolerance stops clearing that padding by 10x.
                tol = _WALL_TOL_FRACTION * (outer_radius - inner_radius)

                for surf in surfaces:
                    bbox = gmsh.model.getBoundingBox(surf[0], surf[1])
                    x_min, y_min, z_min, x_max, y_max, z_max = bbox
                    r_max = np.sqrt(max(x_max**2, y_max**2))

                    # Outer cylindrical boundary
                    if abs(r_max - outer_radius) < tol:
                        outer_boundary_surfaces.append(surf[1])
                    # Inner cylinder surface
                    elif abs(r_max - inner_radius) < tol:
                        inner_boundary_surfaces.append(surf[1])

                if outer_boundary_surfaces:
                    gmsh.model.addPhysicalGroup(2, outer_boundary_surfaces, tag=1)
                    gmsh.model.setPhysicalName(2, 1, "outer_boundary")

                if inner_boundary_surfaces:
                    gmsh.model.addPhysicalGroup(2, inner_boundary_surfaces, tag=2)
                    gmsh.model.setPhysicalName(2, 2, "inner_boundary")

                # Set mesh size
                gmsh.model.mesh.setSize(gmsh.model.getEntities(0), resolution)

                # Generate mesh
                gmsh.model.mesh.generate(3)
                gmsh.model.mesh.optimize("Netgen")
            except BaseException as exc:  # noqa: BLE001 - re-raised on every rank
                build_error = exc
                if gmsh.isInitialized():
                    gmsh.finalize()

        # Collective: every rank learns whether the geometry exists before
        # any of them enters the collective `_model_to_mesh` below (GEO-23
        # step 2a; the unwrapped throw deadlocked at 120 s instead of
        # failing in seconds).
        _raise_geometry_failure_on_every_rank(
            comm, rank, "cylindrical_domain", build_error, resolution
        )

        # Convert to dolfinx mesh
        mesh, cell_tags, facet_tags = _model_to_mesh(gmsh.model, comm, rank, gdim=3)

        if comm.rank == rank:
            gmsh.finalize()

        return mesh, cell_tags, facet_tags

    @staticmethod
    def two_cylinder_domain(
        separation: float = 0.05,
        radius: float = 0.01,
        length: float = 0.1,
        resolution: float = 0.02,
        comm: MPI.Intracomm = MPI.COMM_WORLD,
        rank: int = 0,
    ) -> Tuple[dolfinx.mesh.Mesh, dolfinx.mesh.MeshTags, dolfinx.mesh.MeshTags]:
        """Generate mesh with two side-by-side cylinders inside a box domain.

        This intentionally avoids boolean operations/fragmentation and keeps
        all three volumes explicitly tagged.

        Parameters
        ----------
        separation : float
            Center-to-center distance between the two cylinders [m]
        radius : float
            Radius of each cylinder [m]
        length : float
            Cylinder length along z-axis [m]
        resolution : float
            Characteristic mesh size [m]
        comm : MPI.Intracomm
            MPI communicator
        rank : int
            Rank for Gmsh model (usually 0)

        Returns
        -------
        mesh : dolfinx.mesh.Mesh
            The generated mesh
        cell_tags : dolfinx.mesh.MeshTags
            Cell tags for subdomains (cylinder_1=1, cylinder_2=2, domain=3)
        facet_tags : dolfinx.mesh.MeshTags
            Facet tags for boundaries
        """
        if comm.rank == rank:
            gmsh.initialize()
            gmsh.model.add("two_cylinder_domain")

            x_offset = separation / 2

            cylinder_1 = gmsh.model.occ.addCylinder(-x_offset, 0, -length / 2, 0, 0, length, radius)
            cylinder_2 = gmsh.model.occ.addCylinder(x_offset, 0, -length / 2, 0, 0, length, radius)

            box_half_x = x_offset + 2.0 * radius
            box_half_y = 2.0 * radius
            domain = gmsh.model.occ.addBox(
                -box_half_x,
                -box_half_y,
                -length / 2,
                2.0 * box_half_x,
                2.0 * box_half_y,
                length,
            )

            gmsh.model.occ.synchronize()

            gmsh.model.addPhysicalGroup(3, [cylinder_1], tag=1)
            gmsh.model.setPhysicalName(3, 1, "cylinder_1")

            gmsh.model.addPhysicalGroup(3, [cylinder_2], tag=2)
            gmsh.model.setPhysicalName(3, 2, "cylinder_2")

            gmsh.model.addPhysicalGroup(3, [domain], tag=3)
            gmsh.model.setPhysicalName(3, 3, "domain")

            # Tag all box boundary surfaces as outer boundary
            boundary_surfaces = []
            x_tol = box_half_x + resolution
            y_tol = box_half_y + resolution
            z_tol = (length / 2) + resolution

            for dim, surf in gmsh.model.getEntities(dim=2):
                x_min, y_min, z_min, x_max, y_max, z_max = gmsh.model.getBoundingBox(dim, surf)
                if (
                    abs(abs(x_min) - box_half_x) < x_tol
                    or abs(abs(x_max) - box_half_x) < x_tol
                    or abs(abs(y_min) - box_half_y) < y_tol
                    or abs(abs(y_max) - box_half_y) < y_tol
                    or abs(abs(z_min) - (length / 2)) < z_tol
                    or abs(abs(z_max) - (length / 2)) < z_tol
                ):
                    boundary_surfaces.append(surf)

            if boundary_surfaces:
                gmsh.model.addPhysicalGroup(2, boundary_surfaces, tag=1)
                gmsh.model.setPhysicalName(2, 1, "outer_boundary")

            gmsh.model.mesh.setSize(gmsh.model.getEntities(0), resolution)
            gmsh.model.mesh.generate(3)
            gmsh.model.mesh.optimize("Netgen")

        mesh, cell_tags, facet_tags = _model_to_mesh(gmsh.model, comm, rank, gdim=3)

        if comm.rank == rank:
            gmsh.finalize()

        return mesh, cell_tags, facet_tags

    @staticmethod
    def two_torus_domain(
        separation: float = 0.05,
        major_radius: float = 0.02,
        minor_radius: float = 0.005,
        resolution: float = 0.02,
        comm: MPI.Intracomm = MPI.COMM_WORLD,
        rank: int = 0,
        *,
        air_padding: Optional[float] = None,
        wire_resolution: Optional[float] = None,
        far_resolution: Optional[float] = None,
        port_gap: bool = False,
        gap_angle: float = 0.30,
        gap_clearance: float = 1.0e-3,
        gap_burial: Optional[float] = None,
        gap_overhang: Optional[float] = None,
        gap_arc_resolution: Optional[float] = None,
        gap_arc_tube_radius: Optional[float] = None,
        gap_box_resolution: Optional[float] = None,
        emit_port_sheet: bool = False,
        as_hole: bool = False,
    ) -> Tuple[dolfinx.mesh.Mesh, dolfinx.mesh.MeshTags, dolfinx.mesh.MeshTags]:
        """Generate mesh with two tori inside a box domain.

        The box is ``occ.fragment``-ed against both tori, so the three volumes
        share conforming interfaces and the mesh is a single connected
        component. Cell tags: ``1`` wire 1 (``z < 0``), ``2`` wire 2
        (``z > 0``), ``3`` the air box minus the two tori. Facet tag ``1`` is
        the outer boundary.

        With ``port_gap=True`` each torus additionally carries a lumped-port
        gap (``PORT-1`` step 3b-i): the torus becomes a partial torus of
        opening ``2*pi - gap_angle`` centred on the ``+x`` axis, and a
        rectangular gap box bridges the two arc ends, tagged ``101`` (wire 1)
        and ``102`` (wire 2). Default is **off**, so every existing caller
        meshes exactly as before.

        This fixture is partitioned with ``GhostMode.shared_facet`` (unlike
        the rest of `io/mesh.py`, which takes gmshio's ghost-free default):
        the port facets are interior, so both the tag reconstruction and any
        ``dS`` integral over them need the cell on the far side of a
        partition-boundary facet. See :func:`_interface_facet_tags` and
        known-issues 9.

        Until 2026-08-01 this fixture added the volumes without fragmenting
        (`GEO-8`): gmsh then meshed the box *solid* through the torus regions
        and the tori as two disconnected islands, so a source restricted to a
        torus tag produced no field in the air (``PORT-1`` measured
        ``Z12 == 0``). Do not remove the fragment call.

        Parameters
        ----------
        air_padding:
            Clearance between the coil envelope and the outer box [m]. When
            ``None`` (default) this is ``2 * minor_radius``, which preserves the
            historical box size but couples the air gap to the wire radius --
            making the wire thinner shrinks the air box.

            That coupling matters physically. The outer boundary carries the
            natural condition ``n x (mu^-1 curl A) = 0``, i.e. ``n x H = 0``,
            which behaves as a perfect *magnetic* conductor and mirrors flux
            back into the domain, inflating the on-axis field. Measured centre-
            field error against the analytic Helmholtz solution:

                box half-width 1.45R -> 43.7%
                box half-width 1.75R -> 20.5%   (historical default)
                box half-width 2.50R ->  4.4%

            For free-space comparisons pass ``air_padding >= 2 * major_radius``.
        wire_resolution:
            Target mesh size on the torus surfaces [m]. When ``None``, a single
            uniform size (``resolution``) is used everywhere, as before. Setting
            this enables graded refinement so a large air box stays affordable:
            cell count is otherwise driven by the box volume, not the wire.
        far_resolution:
            Target mesh size at the outer boundary [m]. Defaults to
            ``resolution`` when grading is enabled.
        port_gap:
            Opt in to the gapped (lumped-port) variant. Default ``False``
            reproduces the ungapped fixture byte for byte.
        gap_angle:
            Angular opening of the removed wedge [rad], centred on ``+x``.
        gap_clearance:
            Legacy single margin [m]: the default for **both** ``gap_burial``
            and ``gap_overhang`` below, so a call that passes only this one
            meshes exactly as it did before the split (``PORT-1`` step 3b-i's
            gate is written against that geometry).
        gap_burial:
            How far the gap box buries past each arc end-face centre along
            ``ŷ`` [m], i.e. the box half-length is ``major_radius *
            sin(gap_angle/2) + gap_burial``. Must stay **strictly positive**:
            the two end planes meet at ``gap_angle``, so a box face flush with
            a tilted arc end is not constructible (step 3b-i's recorded
            deviation). Defaults to ``gap_clearance``.
        gap_overhang:
            Transverse margin by which the box overhangs the tube in the
            ``xz`` plane [m], i.e. the box half-size there is ``minor_radius +
            gap_overhang``. Unlike the burial this one is free to vanish, and
            it is the one that costs cross-section: at overhang ``o`` the box
            face is ``4(r+o)²`` against the tube's ``πr²``, so the fraction of
            the face that is *not* conductor shadow is ``1 - πr²/(4(r+o)²)``
            — 45.5% at the 1 mm legacy value, and never below ``1 - π/4 =
            21.5%`` however small ``o`` gets (a rectangle circumscribing a
            circle keeps its corners). ``PORT-1`` step 3b-ii traced a +72.12%
            mutual-impedance error to that fringe; step 3b-iii varied this
            parameter and found the box-volume voltage *sign-unstable* in it
            (+1.7210 / −0.2391 / +0.3317 × ωM₁₂ at o = 1e-3 / 5e-4 / 2e-4),
            which excludes the box-average family outright. Steps 3b-v and
            3b-vi both run at ``o = 2e-4``; note that below ``o ≈ 6e-4`` the
            tube protrudes through the box's ``−x`` face, so the port facet
            tags pick up lateral strips (known-issues 11). Defaults to
            ``gap_clearance``.
        gap_arc_resolution:
            Target mesh size in a tube around each gap's *centreline arc* [m].
            ``None`` (default) leaves the sizing exactly as before. Requires
            ``wire_resolution`` (it composes with that graded field through a
            ``Min``), and is ignored unless ``port_gap``.

            `PORT-1` step 3b-vii: the lumped-port voltage estimator integrates
            ``E·t̂`` along that arc, and N1curl guarantees continuity of the
            *facet*-tangential component only — the arc's own tangent is not
            facet-tangential, so the integrand jumps at every cell crossing.
            Step 3b-vi measured the consequence: with ``h_wire = 2.5e-3``
            against an arc length ``a·g = 1.2e-2`` only ~5 cells span the path,
            and the quadrature converges at ``O(1/n)`` to a ~1e-3 plateau
            instead of the 1e-3 precondition. The fix is to resolve the
            integrand, not the quadrature; ``3e-4`` puts ~40 cells across the
            arc.
        gap_arc_tube_radius:
            Radius of that refined tube [m]. Defaults to
            ``4 * gap_arc_resolution``, which is wide enough that the cells the
            arc actually passes through are all at ``gap_arc_resolution``. The
            size then ramps linearly to ``far_resolution`` at slope 0.3 (about
            1.3x growth per cell, gmsh's comfortable gradation), which is what
            keeps the refinement affordable: the tube itself is a few 1e-8 m^3.
        gap_box_resolution:
            Target mesh size *inside each gap box* [m], overhang included.
            ``None`` (default) leaves the sizing exactly as before — every
            digit-string the `PORT-1` gates pin is measured at that default.
            Requires ``wire_resolution``; ignored unless ``port_gap``.

            `PORT-1` step 3b-xvi. ``gap_arc_resolution`` refines a tube of
            radius ``4*gap_arc_resolution`` around the gap's centreline arc,
            which is where the path-integral estimator samples — but that tube
            never reaches the conductor wall at ``minor_radius``, so the box's
            transverse structure, and the ``gap_overhang`` shell in
            particular, is still sized by the conductor field at
            ``wire_resolution``. At the landed fixture that leaves a 2e-4
            overhang inside 2.5e-3 cells, i.e. sub-cell, while step 3b-ix
            located 45% of the loop EMF there. This knob is the feed-region
            h-refinement that makes the estimator's mesh convergence
            answerable; it is deliberately **local** (a ``Box`` field with a
            slope-0.3 ramp, composed through the same ``Min``), so the PEC-box
            deficit stays common-mode with the unrefined record.
        emit_port_sheet:
            Opt in to the longitudinal **port-sheet mid-plane** inside each gap
            box (`GEO-16`, the `PORT-9` step-1 mesh prerequisite). Default
            ``False`` reproduces the gapped fixture byte for byte — every
            `PORT-1` / `PORT-10` record was measured without it. Requires
            ``port_gap``.

            A lumped-element port sheet spans terminal to terminal with the
            port current flowing **in** its plane. The fixture's only tagged
            surfaces (facet tags ``201`` / ``202``) are the gap<->conductor
            cross-sections, which are *normal* to that current — the wrong
            constitutive law for ``R = Z_p * w / h``. This kwarg fragments each
            gap box with its own mid-plane ``z = +/- separation/2``, the plane
            that contains both the gap's centreline arc and the current
            direction ``y``, so the tet mesh conforms to it.

            Consequences for the tags, all opt-in:

            * each gap box becomes **two** cell groups — ``101`` / ``111`` for
              gap 1 (below / above its mid-plane) and ``102`` / ``112`` for
              gap 2. A caller that selects the gap volume by tag must take
              both halves.
            * the sheet itself is facet tag ``211`` (gap 1) / ``212`` (gap 2),
              rebuilt from the distributed cell tags exactly like the port
              facets (known-issues 9 forbids a dim-2 gmsh group on an interior
              surface).
            * the port facet groups ``201`` / ``202`` are unchanged as *sets*:
              the mid-plane also cuts the arc-end discs, so each is now the
              union of two cell-tag interfaces.

            The sheet's measured extents (area, length ``h`` along the current
            direction, transverse width ``w``) are printed, never gated: the
            gap box crosses a round arc, so the "number of squares" that
            ``R = Z_p * w / h`` needs is a measured quantity on this fixture,
            not the box's nominal dimensions.
        as_hole:
            Opt in to the **PEC-hole** variant (`TH-15` step 2a): the two
            conductors are *cut* out of the air box instead of fragmented into
            it, their volumes are dropped from the model, and the cavity wall
            arrives as facet tag :data:`TWO_TORUS_CONDUCTOR_SURFACE_TAG`
            (``301``).  Default ``False`` leaves every existing caller — every
            `PORT-1` / `PORT-10` / `PORT-12` record — meshing exactly as before.

            The meshed cell set is then air (tag ``3``) plus the four gap-box
            halves (``101`` / ``111`` / ``102`` / ``112``) only: cell tags ``1``
            and ``2`` do not exist, and neither do the gap<->conductor port
            facet groups ``201`` / ``202``, which are cell-tag interfaces of a
            conductor that is no longer meshed.  The port-sheet mid-planes
            (``211`` / ``212``) are unchanged in construction and in area.

            Requires ``port_gap`` and ``emit_port_sheet``: the hole route exists
            to carry the lumped-sheet ports of `TH-15` step 2, and a hole with
            no sheet has no port at all.

            Two rules, both measured by `TH-15` step 0
            (`tests/mesh/probe_two_torus_conductor_hole.py`):

            * ``removeTool=False`` on both cuts.  With ``True`` the cut deletes
              the tool's surfaces, and any group standing on them goes with it.
            * the conductor-surface group is built from ``getBoundary`` of the
              **meshed** volumes, never from the retained tool's own faces.  The
              retained solid's faces are *not* the cavity's: tagging them made
              ``_model_to_mesh`` abort with 716 ``nodes not attached to any
              tet``.
        """
        if as_hole and not (port_gap and emit_port_sheet):
            raise ValueError(
                "as_hole needs port_gap and emit_port_sheet: the hole route is "
                "the lumped-sheet port fixture with the conductors removed "
                f"(got port_gap={port_gap!r}, emit_port_sheet={emit_port_sheet!r})"
            )
        if emit_port_sheet and not port_gap:
            raise ValueError(
                "emit_port_sheet needs port_gap: the sheet is the gap box's "
                "mid-plane, and the ungapped fixture has no gap box"
            )

        if comm.rank == rank:
            gmsh.initialize()
            gmsh.model.add("two_torus_domain")

            z_offset = separation / 2

            if port_gap:
                if not 0.0 < gap_angle < np.pi:
                    raise ValueError(f"gap_angle must be in (0, pi), got {gap_angle!r}")
                burial = gap_clearance if gap_burial is None else float(gap_burial)
                overhang = gap_clearance if gap_overhang is None else float(gap_overhang)
                if burial <= 0.0:
                    raise ValueError(
                        f"gap_burial must be positive, got {burial!r} "
                        f"(gap_clearance={gap_clearance!r})"
                    )
                if overhang < 0.0:
                    raise ValueError(
                        f"gap_overhang must be non-negative, got {overhang!r} "
                        f"(gap_clearance={gap_clearance!r})"
                    )
                # Partial torus: OCC sweeps from phi = 0 to phi = angle, so
                # rotating by +gap_angle/2 leaves the wedge centred on +x.
                wire_1 = gmsh.model.occ.addTorus(
                    0,
                    0,
                    -z_offset,
                    major_radius,
                    minor_radius,
                    angle=2.0 * np.pi - gap_angle,
                )
                wire_2 = gmsh.model.occ.addTorus(
                    0,
                    0,
                    z_offset,
                    major_radius,
                    minor_radius,
                    angle=2.0 * np.pi - gap_angle,
                )
                for wire in (wire_1, wire_2):
                    gmsh.model.occ.rotate([(3, wire)], 0, 0, 0, 0, 0, 1, 0.5 * gap_angle)
                # The box must *cross* both arc ends, not stop short of them:
                # a box face flush with a tilted arc end is not constructible
                # (the two end planes meet at gap_angle), so the box buries
                # `gap_burial` past the end-face centres and the gap group
                # takes precedence over the conductor below. That makes the
                # gap the box exactly -- planar faces, meshed to roundoff --
                # and the conductor the arc minus what the box swallowed.
                #
                # The two margins are independent (step 3b-iii): `burial` is
                # the along-arc one that must stay positive for the box to be
                # constructible at all; `overhang` is the transverse one, which
                # buys nothing and costs cross-section (see the docstring).
                gap_half_xz = minor_radius + overhang
                gap_half_y = major_radius * np.sin(0.5 * gap_angle) + burial
                gap_size = (2.0 * gap_half_xz, 2.0 * gap_half_y, 2.0 * gap_half_xz)
                gap_1 = gmsh.model.occ.addBox(
                    major_radius - gap_half_xz,
                    -gap_half_y,
                    -z_offset - gap_half_xz,
                    *gap_size,
                )
                gap_2 = gmsh.model.occ.addBox(
                    major_radius - gap_half_xz,
                    -gap_half_y,
                    z_offset - gap_half_xz,
                    *gap_size,
                )
                # `GEO-16`: the longitudinal mid-plane of each gap box, added
                # as a *tool* to the same fragment below so the tets conform
                # to it. Exactly the box's own cross-section at ``z = z_c`` —
                # no larger, or the fragment would also cut the conductor and
                # the air box along that plane.
                sheet_1 = sheet_2 = None
                if emit_port_sheet:
                    sheet_1 = gmsh.model.occ.addRectangle(
                        major_radius - gap_half_xz,
                        -gap_half_y,
                        -z_offset,
                        2.0 * gap_half_xz,
                        2.0 * gap_half_y,
                    )
                    sheet_2 = gmsh.model.occ.addRectangle(
                        major_radius - gap_half_xz,
                        -gap_half_y,
                        z_offset,
                        2.0 * gap_half_xz,
                        2.0 * gap_half_y,
                    )
            else:
                wire_1 = gmsh.model.occ.addTorus(0, 0, -z_offset, major_radius, minor_radius)
                wire_2 = gmsh.model.occ.addTorus(0, 0, z_offset, major_radius, minor_radius)

            padding = 2.0 * minor_radius if air_padding is None else float(air_padding)
            if padding <= 0.0:
                raise ValueError(f"air_padding must be positive, got {padding!r}")

            radial_extent = major_radius + minor_radius
            box_half_x = radial_extent + padding
            box_half_y = radial_extent + padding
            box_half_z = z_offset + minor_radius + padding

            domain = gmsh.model.occ.addBox(
                -box_half_x,
                -box_half_y,
                -box_half_z,
                2.0 * box_half_x,
                2.0 * box_half_y,
                2.0 * box_half_z,
            )

            cond_of: Dict[int, int] = {}
            if as_hole:
                # `TH-15` step 2a, step 0's measured sequence verbatim:
                #   conductors = cut(arcs, gap boxes)   # gap wins over metal
                #   air        = cut(box, conductors)   # the cavity
                #   fragment(air, [gap boxes, sheets])  # gap groups + mid-plane
                # `removeTool=False` throughout — the first cut keeps the gap
                # boxes for the fragment below, the second keeps the conductor
                # solids alive so their surfaces stay addressable until the
                # cavity wall has been identified from the *meshed* volumes.
                conductors, _ = gmsh.model.occ.cut(
                    [(3, wire_1), (3, wire_2)],
                    [(3, gap_1), (3, gap_2)],
                    removeObject=True,
                    removeTool=False,
                )
                gmsh.model.occ.synchronize()
                for dim, tag in conductors:
                    _, _, zc = gmsh.model.occ.getCenterOfMass(dim, tag)
                    cond_of[tag] = 1 if zc < 0.0 else 2
                if sorted(cond_of.values()) != [1, 2]:
                    raise RuntimeError(
                        "two_torus_domain(as_hole=True): the arc cut did not "
                        f"leave one conductor solid per side (got {conductors}, "
                        f"z-classification {cond_of})"
                    )
                air_pieces, _ = gmsh.model.occ.cut(
                    [(3, domain)],
                    conductors,
                    removeObject=True,
                    removeTool=False,
                )
                tool_dimtags = [(3, gap_1), (3, gap_2), (2, sheet_1), (2, sheet_2)]
                _, fragment_map = gmsh.model.occ.fragment(air_pieces, tool_dimtags)
                gmsh.model.occ.synchronize()
                input_dimtags = list(air_pieces) + tool_dimtags
            else:
                tool_tags = [wire_1, wire_2]
                if port_gap:
                    tool_tags += [gap_1, gap_2]
                tool_dimtags = [(3, tag) for tag in tool_tags]
                if port_gap and emit_port_sheet:
                    tool_dimtags += [(2, sheet_1), (2, sheet_2)]
                _, fragment_map = gmsh.model.occ.fragment([(3, domain)], tool_dimtags)
                gmsh.model.occ.synchronize()
                input_dimtags = [(3, domain)] + tool_dimtags

            if port_gap:
                # Groups re-derived from the fragment out-map (`GEO-9` step 2b
                # machinery): fragment renumbers, so absolute tags from before
                # the call mean nothing after it. The out-map is positional,
                # objects first then tools.
                if len(fragment_map) != len(input_dimtags):
                    raise RuntimeError(
                        "two_torus_domain: occ.fragment returned an out-map of "
                        f"{len(fragment_map)} entries for {len(input_dimtags)} "
                        "inputs"
                    )
                ancestors: Dict[int, set] = {}
                for input_dimtag, pieces in zip(input_dimtags, fragment_map):
                    for dim, piece in pieces:
                        if dim == 3:
                            ancestors.setdefault(piece, set()).add(input_dimtag)

                # With `as_hole` the conductors were cut away before the
                # fragment, so no piece descends from a wire and every non-gap
                # piece is air.
                wire_of = {} if as_hole else {(3, wire_1): 1, (3, wire_2): 2}
                # `GEO-16`: with the sheet in, each gap box is two pieces —
                # below and above its own mid-plane. The sheet is a dim-2
                # input, so it is never an *ancestor* of a 3-D piece; the
                # halves are told apart by centroid z against the box centre,
                # the same mass/centroid discipline the ungapped branch uses.
                gap_of = {(3, gap_1): (101, 111, -z_offset), (3, gap_2): (102, 112, z_offset)}
                group_of_piece: Dict[int, int] = {}
                for piece, sources in ancestors.items():
                    # Gap wins over metal: the box IS the dielectric gap, and
                    # its exact volume is what the step-3b-i anchor measures.
                    hit_gap = sources & gap_of.keys()
                    hit_wire = sources & wire_of.keys()
                    if hit_gap:
                        lower, upper, z_c = gap_of[min(hit_gap)]
                        if emit_port_sheet:
                            _, _, zc_piece = gmsh.model.occ.getCenterOfMass(3, piece)
                            group_of_piece[piece] = upper if zc_piece > z_c else lower
                        else:
                            group_of_piece[piece] = lower
                    elif hit_wire:
                        group_of_piece[piece] = min(wire_of[t] for t in hit_wire)
                    else:
                        group_of_piece[piece] = 3

                pieces_by_group: Dict[int, List[int]] = {}
                for piece, group in sorted(group_of_piece.items()):
                    pieces_by_group.setdefault(group, []).append(piece)

                group_names = {
                    1: "wire_1",
                    2: "wire_2",
                    3: "domain",
                    101: "gap_1",
                    102: "gap_2",
                    111: "gap_1_upper",
                    112: "gap_2_upper",
                }
                required_groups = (3, 101, 102) if as_hole else (1, 2, 3, 101, 102)
                if emit_port_sheet:
                    required_groups = required_groups + (111, 112)
                # The retained conductor solids are still in the model at this
                # point (`removeTool=False`); they are not meshed volumes and
                # must not enter the mass census or the ungrouped-entity check.
                volumes = [vt for vt in gmsh.model.getEntities(dim=3) if vt[1] not in cond_of]
                masses = {tag: gmsh.model.occ.getMass(3, tag) for _, tag in volumes}
                missing = [group_names[g] for g in required_groups if g not in pieces_by_group]
                if missing:
                    raise RuntimeError(
                        "two_torus_domain: occ.fragment left no piece for "
                        f"{', '.join(missing)}; {len(volumes)} volumes, "
                        "per-volume masses [m^3]: "
                        + ", ".join(f"tag {t}: {m:.6e}" for t, m in sorted(masses.items()))
                    )

                for group, pieces in sorted(pieces_by_group.items()):
                    gmsh.model.addPhysicalGroup(3, pieces, tag=group)
                    gmsh.model.setPhysicalName(3, group, group_names[group])

                # Every 3-D entity must carry a marker; a cell with none is
                # what `gmshio.py:118` asserts on. Checked against the model,
                # not the out-map, so an unmentioned piece cannot hide.
                grouped_volumes = set()
                for _, group_tag in gmsh.model.getPhysicalGroups(dim=3):
                    grouped_volumes.update(gmsh.model.getEntitiesForPhysicalGroup(3, group_tag))
                ungrouped = sorted({tag for _, tag in volumes} - grouped_volumes)
                if ungrouped:
                    raise RuntimeError(
                        "two_torus_domain: 3-D entities carry no physical group: "
                        f"{ungrouped}; {len(volumes)} volumes, masses [m^3]: "
                        + ", ".join(f"tag {t}: {masses[t]:.6e}" for t in ungrouped)
                    )

                print(
                    f"[two-torus-mesh] gapped fragment volumes={len(volumes)} "
                    + " ".join(
                        f"{group_names[g]}={sum(masses[p] for p in pieces):.6e}" f"({len(pieces)}p)"
                        for g, pieces in sorted(pieces_by_group.items())
                    )
                    + f" gap_box_analytic={gap_size[0] * gap_size[1] * gap_size[2]:.6e}",
                    flush=True,
                )

                wire_volumes = (
                    {1: [], 2: []} if as_hole else {1: pieces_by_group[1], 2: pieces_by_group[2]}
                )

                # `PORT-1` step 3b-iv: the port facet groups (`201` / `202`)
                # are the surfaces each gap piece shares with its conductor
                # piece -- exactly the two planar cuts where the gap box's
                # y-faces cross the arc, since the box overhangs the tube in x
                # and z (`gap_half_xz > minor_radius`) and the tube can only
                # leave it through those faces. That is the gap-conductor
                # interface a lumped-port voltage integrates over (step 3b-v).
                #
                # They are NOT emitted as gmsh physical groups: dim-2 groups on
                # facets interior to the partition hang `model_to_mesh` at
                # `-n 2` (known-issues 9). The identical facet set is rebuilt
                # from the distributed *cell* tags after `model_to_mesh`, via
                # `_interface_facet_tags`. What is computed here is the CAD
                # cross-check -- OCC's own area for the surfaces, printed so
                # the meshed area has an independent number to be scored
                # against. Derived from the fragment's own boundaries;
                # absolute tags from before the fragment call mean nothing
                # after it.
                def _boundary_surfaces(groups):
                    return {
                        surf
                        for group in groups
                        for vol in pieces_by_group[group]
                        for _, surf in gmsh.model.getBoundary(
                            [(3, vol)], oriented=False, recursive=False
                        )
                    }

                # `GEO-16`: with the sheet in, "the gap" is both halves.
                gap_groups = {101: (101, 111), 102: (102, 112)}

                cad_areas = {}
                for gap_group, wire_group, facet_group in (
                    () if as_hole else ((101, 1, 201), (102, 2, 202))
                ):
                    gap_boundary = _boundary_surfaces(
                        gap_groups[gap_group] if emit_port_sheet else (gap_group,)
                    )
                    conductor_boundary = _boundary_surfaces((wire_group,))
                    shared = sorted(gap_boundary & conductor_boundary)
                    if not shared:
                        raise RuntimeError(
                            "two_torus_domain: gap group "
                            f"{group_names[gap_group]} shares no surface with "
                            f"{group_names[wire_group]}; the gap box is not "
                            "crossing the arc ends"
                        )
                    cad_areas[facet_group] = (
                        len(shared),
                        sum(gmsh.model.occ.getMass(2, s) for s in shared),
                    )

                print(
                    "[two-torus-mesh] port interfaces (CAD) "
                    + " ".join(
                        f"{tag}: {n} surface(s) area={area:.6e}"
                        for tag, (n, area) in sorted(cad_areas.items())
                    ),
                    flush=True,
                )

                if emit_port_sheet:
                    # `GEO-16`: the sheet is what the two halves of one gap box
                    # share -- the mid-plane, and nothing else, because the
                    # halves touch only there. Same CAD-denominator pattern as
                    # the port interfaces above: OCC's own area, plus the
                    # bounding-box extents `R = Z_p * w / h` needs (h along the
                    # current direction y, w transverse in x).
                    for lower, upper, facet_group in ((101, 111, 211), (102, 112, 212)):
                        shared = sorted(_boundary_surfaces((lower,)) & _boundary_surfaces((upper,)))
                        if not shared:
                            raise RuntimeError(
                                "two_torus_domain: gap halves "
                                f"{group_names[lower]}/{group_names[upper]} "
                                "share no surface; the fragment did not "
                                "conform to the port-sheet mid-plane"
                            )
                        area = sum(gmsh.model.occ.getMass(2, s) for s in shared)
                        boxes = np.array([gmsh.model.getBoundingBox(2, s) for s in shared])
                        w_cad = boxes[:, 3].max() - boxes[:, 0].min()
                        h_cad = boxes[:, 4].max() - boxes[:, 1].min()
                        print(
                            f"[two-torus-mesh] port sheet (CAD) {facet_group}: "
                            f"{len(shared)} surface(s) area={area:.9e} "
                            f"h_along_current={h_cad:.9e} w_transverse={w_cad:.9e} "
                            f"squares_w_over_h={w_cad / h_cad:.9f} "
                            f"nominal_area={4.0 * gap_half_xz * gap_half_y:.9e}",
                            flush=True,
                        )

                cavity_surfaces: List[int] = []
                if as_hole:
                    # `TH-15` step 2a.  The cavity wall must be identified from
                    # the *meshed* volumes alone: the retained conductor solids
                    # are still in the model, and testing a face against a set
                    # that includes them is self-satisfying (step 0's abort:
                    # 716 `nodes not attached to any tet`).  A cavity wall is a
                    # face that bounds exactly ONE meshed volume and is not flat
                    # against an outer wall — the air/gap interfaces and the
                    # sheet mid-planes bound two, the outer walls bound one but
                    # lie on a wall.
                    def _face_on_wall(surf: int) -> bool:
                        bb = gmsh.model.getBoundingBox(2, surf)
                        return any(
                            abs(bb[i] - lo) < 1e-6 and abs(bb[i + 3] - lo) < 1e-6
                            for i, lo in (
                                (0, -box_half_x),
                                (0, box_half_x),
                                (1, -box_half_y),
                                (1, box_half_y),
                                (2, -box_half_z),
                                (2, box_half_z),
                            )
                        )

                    face_use: Dict[int, int] = {}
                    for _, vol in volumes:
                        for _, surf in gmsh.model.getBoundary(
                            [(3, vol)], oriented=False, recursive=False
                        ):
                            face_use[surf] = face_use.get(surf, 0) + 1
                    cavity_surfaces = sorted(
                        s for s, n in face_use.items() if n == 1 and not _face_on_wall(s)
                    )
                    if not cavity_surfaces:
                        raise RuntimeError(
                            "two_torus_domain(as_hole=True): the meshed volumes "
                            "have no single-use non-wall boundary face, so the "
                            "conductor cavity surface is missing"
                        )
                    cavity_area = sum(gmsh.model.occ.getMass(2, s) for s in cavity_surfaces)

                    # Drop the conductor solids: the cavity IS the hole.
                    # `recursive=False` so the faces survive for the group.
                    gmsh.model.occ.remove([(3, t) for t in cond_of], recursive=False)
                    gmsh.model.occ.synchronize()
                    still_present = [t for t in cond_of if (3, t) in gmsh.model.getEntities(dim=3)]
                    if still_present:
                        # occ.remove + synchronize did not propagate; drop them
                        # model-side (step 0 needed this fallback).
                        gmsh.model.removeEntities([(3, t) for t in still_present], recursive=False)
                    cavity_surfaces = [
                        s for s in cavity_surfaces if (2, s) in gmsh.model.getEntities(dim=2)
                    ]
                    if not cavity_surfaces:
                        raise RuntimeError(
                            "two_torus_domain(as_hole=True): every cavity "
                            "surface was destroyed with the conductor solids"
                        )
                    gmsh.model.addPhysicalGroup(
                        2, cavity_surfaces, tag=TWO_TORUS_CONDUCTOR_SURFACE_TAG
                    )
                    gmsh.model.setPhysicalName(
                        2, TWO_TORUS_CONDUCTOR_SURFACE_TAG, "conductor_surface"
                    )
                    print(
                        "[two-torus-mesh] conductor cavity (CAD) tag "
                        f"{TWO_TORUS_CONDUCTOR_SURFACE_TAG}: "
                        f"{len(cavity_surfaces)} surface(s) "
                        f"area={cavity_area:.9e}; volumes after the drop="
                        f"{len(gmsh.model.getEntities(dim=3))} "
                        f"(conductor solids removed: {len(cond_of)})",
                        flush=True,
                    )

            else:
                # Fragment renumbers volumes; identify them by mass and centroid
                # rather than by the tags it hands back (same discipline as
                # loop_over_half_space_domain). Each torus is orders of magnitude
                # smaller than the air region, and the two tori are told apart by
                # the sign of their centroid z.
                torus_mass = 2.0 * np.pi**2 * major_radius * minor_radius**2
                wire_1_volume = None
                wire_2_volume = None
                air_volume = None
                for dim, tag in gmsh.model.getEntities(dim=3):
                    mass = gmsh.model.occ.getMass(dim, tag)
                    _, _, zc = gmsh.model.occ.getCenterOfMass(dim, tag)
                    if mass > 10.0 * torus_mass:
                        air_volume = tag
                    elif zc < 0.0:
                        wire_1_volume = tag
                    else:
                        wire_2_volume = tag

                if wire_1_volume is None or wire_2_volume is None or air_volume is None:
                    raise RuntimeError(
                        "two_torus_domain: fragment did not produce the expected two "
                        f"tori plus air volume (got {gmsh.model.getEntities(dim=3)})"
                    )

                gmsh.model.addPhysicalGroup(3, [wire_1_volume], tag=1)
                gmsh.model.setPhysicalName(3, 1, "wire_1")

                gmsh.model.addPhysicalGroup(3, [wire_2_volume], tag=2)
                gmsh.model.setPhysicalName(3, 2, "wire_2")

                gmsh.model.addPhysicalGroup(3, [air_volume], tag=3)
                gmsh.model.setPhysicalName(3, 3, "domain")

                wire_volumes = {1: [wire_1_volume], 2: [wire_2_volume]}

            # A face is on the outer boundary only if it is *flat against* a
            # wall (both bounding-box extremes on it), not merely within one
            # mesh size of it: fragment introduces interior faces, and the old
            # `< resolution` test would have swept some of them into the BC.
            #
            # `GEO-10`: the tolerance may not be roundoff-tight. gmsh inflates
            # an OCC entity's bounding box by its geometric tolerance, measured
            # 2026-08-06 at exactly **1.000e-07** on all six walls of this box
            # (`20260806T050143Z_GEO-10-probe.log`) — so `tol = 1e-9` rejected
            # every wall, `boundary_surfaces` came out empty, and the group was
            # never declared at all (known-issues 10). 1e-6 clears that padding
            # by 10x while staying four orders below the nearest interior face,
            # whose residual in the same probe is 2.000e-02.
            tol = 1e-6
            boundary_surfaces = []
            for dim, surf in gmsh.model.getEntities(dim=2):
                x_min, y_min, z_min, x_max, y_max, z_max = gmsh.model.getBoundingBox(dim, surf)
                on_wall = (
                    abs(x_min + box_half_x) < tol
                    and abs(x_max + box_half_x) < tol
                    or abs(x_min - box_half_x) < tol
                    and abs(x_max - box_half_x) < tol
                    or abs(y_min + box_half_y) < tol
                    and abs(y_max + box_half_y) < tol
                    or abs(y_min - box_half_y) < tol
                    and abs(y_max - box_half_y) < tol
                    or abs(z_min + box_half_z) < tol
                    and abs(z_max + box_half_z) < tol
                    or abs(z_min - box_half_z) < tol
                    and abs(z_max - box_half_z) < tol
                )
                if on_wall:
                    boundary_surfaces.append(surf)

            if boundary_surfaces:
                gmsh.model.addPhysicalGroup(2, boundary_surfaces, tag=1)
                gmsh.model.setPhysicalName(2, 1, "outer_boundary")

            if wire_resolution is None:
                if gap_arc_resolution is not None:
                    raise ValueError(
                        "gap_arc_resolution needs wire_resolution: it composes "
                        "with the graded field, and uniform sizing has none"
                    )
                if gap_box_resolution is not None:
                    raise ValueError(
                        "gap_box_resolution needs wire_resolution: it composes "
                        "with the graded field, and uniform sizing has none"
                    )
                gmsh.model.mesh.setSize(gmsh.model.getEntities(0), resolution)
            else:
                # Graded sizing: fine on the tori, coarsening with distance so a
                # large air box does not blow up the cell count.
                h_wire = float(wire_resolution)
                h_far = float(resolution if far_resolution is None else far_resolution)
                if h_wire <= 0.0 or h_far <= 0.0:
                    raise ValueError("wire_resolution and far_resolution must be positive")

                # Refine on the conductors and, when gapped, on the gap boxes
                # too -- the gap is the smallest feature in the model. Pieces
                # share surfaces after fragment, so de-duplicate.
                refine_volumes = list(wire_volumes[1]) + list(wire_volumes[2])
                if port_gap:
                    for group in required_groups:
                        if group >= 101:
                            refine_volumes += pieces_by_group[group]
                # `TH-15` step 2a: with the conductors gone there is no wire
                # volume to refine towards, so the cavity wall — the same
                # surface, now a boundary — carries the conductor field.
                wire_surfaces = sorted(
                    (set(cavity_surfaces) if as_hole else set())
                    | {
                        surf
                        for vol in refine_volumes
                        for _, surf in gmsh.model.getBoundary(
                            [(3, vol)], oriented=False, recursive=False
                        )
                    }
                )

                dist = gmsh.model.mesh.field.add("Distance")
                gmsh.model.mesh.field.setNumbers(dist, "SurfacesList", wire_surfaces)
                gmsh.model.mesh.field.setNumber(dist, "Sampling", 200)

                thr = gmsh.model.mesh.field.add("Threshold")
                gmsh.model.mesh.field.setNumber(thr, "InField", dist)
                gmsh.model.mesh.field.setNumber(thr, "SizeMin", h_wire)
                gmsh.model.mesh.field.setNumber(thr, "SizeMax", h_far)
                gmsh.model.mesh.field.setNumber(thr, "DistMin", minor_radius)
                gmsh.model.mesh.field.setNumber(thr, "DistMax", major_radius + padding)

                # Extra *local* refinement fields, composed with the conductor
                # field through a `Min` below (finest wins, nothing coarsens).
                local_fields = []
                if port_gap and gap_arc_resolution is not None:
                    # `PORT-1` step 3b-vii: refine a tube around each gap's
                    # centreline arc, where the path-integral estimator samples.
                    #
                    # The field is defined by *coordinates*, evaluated on the
                    # fragmented model: fragment renumbers, and the arc is not
                    # a model entity at all (it runs through the gap box's
                    # interior), so there is nothing to point a `Distance` field
                    # at. `MathEval` is the distance to the arc written out --
                    # distance to the centreline circle, plus an out-of-band
                    # penalty in `y` so only the wedge is refined and not the
                    # whole 2*pi*a of conductor, plus one in `-x` so the
                    # circle's far branch is excluded. `max(0, u)` is spelled
                    # `(u + sqrt(u^2))/2` and `|y|` as `sqrt(y^2)`: gmsh's
                    # expression parser is not guaranteed to carry `fabs`/`max`,
                    # and sqrt of a square is exact for both.
                    h_gap = float(gap_arc_resolution)
                    if h_gap <= 0.0:
                        raise ValueError(f"gap_arc_resolution must be positive, got {h_gap!r}")
                    tube = (
                        4.0 * h_gap if gap_arc_tube_radius is None else float(gap_arc_tube_radius)
                    )
                    if tube <= 0.0:
                        raise ValueError(f"gap_arc_tube_radius must be positive, got {tube!r}")
                    # Slope 0.3 in size per unit distance ~ 1.3x growth per
                    # cell. `SizeMax = h_far` (not h_wire) is what keeps the
                    # `Min` below from clamping the whole air box to the wire
                    # size; the conductor field still wins where it is finer.
                    dist_max = tube + (h_far - h_gap) / 0.3
                    # OPS-18 step 3: every number interpolated into a gmsh
                    # `MathEval` string must be a plain Python float. numpy 2
                    # renders `repr(np.float64(x))` as `np.float64(x)`, a token
                    # gmsh's expression parser rejects with SIGABRT at meshing
                    # ("invalid token on expression"); numpy 1.x's bare-number
                    # repr masked this. Coerce at the boundary, never `!r` a
                    # value of unknown type into a foreign grammar.
                    arc_half_y = float(major_radius * np.sin(0.5 * gap_angle))
                    r_major = float(major_radius)
                    gap_fields = []
                    for z_c in (-float(z_offset), float(z_offset)):
                        band = (
                            f"(0.5*((sqrt(y^2)-{arc_half_y!r})"
                            f"+sqrt((sqrt(y^2)-{arc_half_y!r})^2)))"
                        )
                        expr = (
                            f"sqrt((sqrt(x^2+y^2)-{r_major!r})^2"
                            f"+(z-({z_c!r}))^2"
                            f"+{band}^2"
                            f"+(0.5*(sqrt(x^2)-x))^2)"
                        )
                        arc_dist = gmsh.model.mesh.field.add("MathEval")
                        gmsh.model.mesh.field.setString(arc_dist, "F", expr)
                        arc_thr = gmsh.model.mesh.field.add("Threshold")
                        gmsh.model.mesh.field.setNumber(arc_thr, "InField", arc_dist)
                        gmsh.model.mesh.field.setNumber(arc_thr, "SizeMin", h_gap)
                        gmsh.model.mesh.field.setNumber(arc_thr, "SizeMax", h_far)
                        gmsh.model.mesh.field.setNumber(arc_thr, "DistMin", tube)
                        gmsh.model.mesh.field.setNumber(arc_thr, "DistMax", dist_max)
                        gap_fields.append(arc_thr)

                    local_fields.extend(gap_fields)
                    print(
                        "[two-torus-mesh] gap-arc refinement h_gap="
                        f"{h_gap:.3e} tube={tube:.3e} dist_max={dist_max:.3e} "
                        f"arc_half_y={arc_half_y:.6e} "
                        f"cells_across_arc={major_radius * gap_angle / h_gap:.1f}",
                        flush=True,
                    )

                if port_gap and gap_box_resolution is not None:
                    # `PORT-1` step 3b-xvi: refine the *whole gap box* -- the
                    # feed region -- not just the centreline tube.  The arc
                    # field above resolves where the estimator's path integral
                    # samples; it does nothing for the transverse structure the
                    # terminal fields have, because its tube radius
                    # (4*h_gap = 1.2e-3) never reaches the tube wall at
                    # r_wire = 5e-3, let alone the `gap_overhang` shell outside
                    # it.  Jin (3rd ed., ch. 12) refines the feed region
                    # locally as routine practice for exactly this reason.
                    #
                    # A `Box` field, not `MathEval`: the gap solids *are*
                    # axis-aligned boxes, so the refined region is the feed
                    # region exactly, and `Thickness` ramps back out to
                    # `far_resolution` so the refinement stays affordable.
                    #
                    # That ramp is *capped*, and the cap is measured, not
                    # guessed.  A slope-0.3 thickness — `(h_far - h_box)/0.3`,
                    # the rule the arc field above uses — is 0.098 m at
                    # `h_box = 6e-4`: a 10 cm shell of refinement around each
                    # gap box.  The arc field survives the same slope because
                    # it ramps away from a thin tube; a `Box` field ramps away
                    # from a *surface*, so the same slope buys vastly more
                    # volume, and the step 3b-xvi mesh arm measured the
                    # consequence — cells outside the gap boxes +35.4560%
                    # against the < 5% locality band
                    # (`20260812T094005Z_PORT-1-step3bxvi-mesh6e4.log`).  Past
                    # a few mm the conductor field (`h_wire` near the wire) is
                    # the finer of the two anyway and the `Min` takes it back,
                    # so the cap removes no resolution the refinement was
                    # buying: it only stops the shell leaking into the air.
                    h_box = float(gap_box_resolution)
                    if h_box <= 0.0:
                        raise ValueError(f"gap_box_resolution must be positive, got {h_box!r}")
                    box_thickness = min(max(h_far - h_box, 0.0) / 0.3, GAP_BOX_THICKNESS_CAP_M)
                    for z_c in (-z_offset, z_offset):
                        box_field = gmsh.model.mesh.field.add("Box")
                        gmsh.model.mesh.field.setNumber(
                            box_field, "XMin", major_radius - gap_half_xz
                        )
                        gmsh.model.mesh.field.setNumber(
                            box_field, "XMax", major_radius + gap_half_xz
                        )
                        gmsh.model.mesh.field.setNumber(box_field, "YMin", -gap_half_y)
                        gmsh.model.mesh.field.setNumber(box_field, "YMax", gap_half_y)
                        gmsh.model.mesh.field.setNumber(box_field, "ZMin", z_c - gap_half_xz)
                        gmsh.model.mesh.field.setNumber(box_field, "ZMax", z_c + gap_half_xz)
                        gmsh.model.mesh.field.setNumber(box_field, "VIn", h_box)
                        gmsh.model.mesh.field.setNumber(box_field, "VOut", h_far)
                        gmsh.model.mesh.field.setNumber(box_field, "Thickness", box_thickness)
                        local_fields.append(box_field)
                    print(
                        "[two-torus-mesh] gap-box refinement h_box="
                        f"{h_box:.3e} thickness={box_thickness:.3e} "
                        f"(uncapped {max(h_far - h_box, 0.0) / 0.3:.3e}, cap "
                        f"{GAP_BOX_THICKNESS_CAP_M:.3e}) "
                        f"half_xz={gap_half_xz:.6e} half_y={gap_half_y:.6e} "
                        f"cells_across_overhang={overhang / h_box:.3f}",
                        flush=True,
                    )

                if local_fields:
                    background = gmsh.model.mesh.field.add("Min")
                    gmsh.model.mesh.field.setNumbers(background, "FieldsList", [thr] + local_fields)
                else:
                    background = thr

                gmsh.model.mesh.field.setAsBackgroundMesh(background)
                # Background field must win over point/curve-derived sizing.
                gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
                gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
                gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)

            gmsh.model.mesh.generate(3)
            gmsh.model.mesh.optimize("Netgen")

        # `PORT-1` step 3b-iv: gmshio's default partitioner is
        # `GhostMode.none`, which leaves `cells_ghost = 0` on every rank. The
        # port facets rebuilt below are *interior* to the mesh, so classifying
        # one needs the tag of the cell on both sides — and on a partition
        # boundary one of those cells lives on the other rank. `shared_facet`
        # is what makes that cell present as a ghost; it is also what a `dS`
        # integral over those facets needs. Plumbed here only, so no other
        # fixture changes partition.
        # `OPS-18` step 3: 0.11 made `max_facet_to_cell_links` a *required*
        # argument of `create_cell_partitioner` (0.7.2 took the mode alone).
        # 2 is dolfinx's own default elsewhere (`create_mesh`) and its
        # documented value for a non-branching manifold mesh, which every
        # fixture here is — tetrahedra with each interior facet shared by
        # exactly two cells. Passing `None` would mean "no upper bound" and
        # is not what this mesh is.
        partitioner = dolfinx.mesh.create_cell_partitioner(dolfinx.mesh.GhostMode.shared_facet, 2)
        mesh, cell_tags, facet_tags = _model_to_mesh(
            gmsh.model, comm, rank, gdim=3, partitioner=partitioner
        )

        if comm.rank == rank:
            gmsh.finalize()

        if port_gap:
            # `PORT-1` step 3b-iv: rebuild the port facet groups from the
            # distributed cell tags (see the CAD-side comment above and
            # known-issues 9). gap_1 <-> wire_1 is port 201, gap_2 <-> wire_2
            # is port 202.
            if as_hole:
                # `TH-15` step 2a: the conductor is not meshed, so `201`/`202`
                # — gap<->conductor cell-tag interfaces — do not exist. The
                # cavity wall arrives as the exterior gmsh group `301`; only the
                # sheets still have to be rebuilt from the distributed tags.
                interfaces = {211: ((101, 111),), 212: ((102, 112),)}
            elif emit_port_sheet:
                # `GEO-16`: the mid-plane cuts the arc-end discs too, so each
                # port group is the union of its two half-box interfaces; the
                # sheet groups (`211` / `212`) are the halves against each
                # other. Same rebuild-from-cell-tags rule as above.
                interfaces = {
                    201: ((101, 1), (111, 1)),
                    202: ((102, 2), (112, 2)),
                    211: ((101, 111),),
                    212: ((102, 112),),
                }
            else:
                interfaces = {201: (101, 1), 202: (102, 2)}
            facet_tags = _interface_facet_tags(mesh, cell_tags, interfaces, facet_tags)

        return mesh, cell_tags, facet_tags

    @staticmethod
    def loop_over_half_space_domain(
        loop_radius: float = 0.04,
        wire_radius: float = 0.005,
        liftoff: float = 0.015,
        box_half_width: float = 0.10,
        resolution_wire: float = 0.004,
        resolution_near: float = 0.005,
        resolution_far: float = 0.02,
        near_half_width: float = 0.06,
        near_depth: float = 0.05,
        near_height: float = 0.025,
        comm: MPI.Intracomm = MPI.COMM_WORLD,
        rank: int = 0,
    ) -> Tuple[dolfinx.mesh.Mesh, dolfinx.mesh.MeshTags, dolfinx.mesh.MeshTags]:
        """Filamentary loop at height ``liftoff`` over a conductive half-space.

        `MAT-6` step 2: the FEM counterpart of
        :func:`~fem_em_solver.utils.dodd_deeds.coil_impedance_change`.  A torus
        of major radius ``loop_radius`` and minor radius ``wire_radius`` sits in
        the plane ``z = liftoff``; the cube ``[-W, W]³`` around it is split at
        ``z = 0`` into an upper air region and a lower slab that stands in for
        the half-space.  The slab fills the whole lower half of the box, so the
        PEC truncation plane below it sits many skin depths into the conductor
        and the half-space is only truncated where the fields are already dead.

        Cell tags: ``1`` wire, ``2`` air (``z > 0``), ``3`` slab (``z < 0``).
        Facet tag ``1`` is the outer boundary of the cube.

        Sizing is **graded**, deliberately: a single global ``setSize`` fine
        enough to put 3–4 cells across the skin depth would be unaffordable over
        the whole air box, and the tight-padding/uniform-size pattern is exactly
        what cost 20% on Helmholtz (docs/testing/known-issues.md).  Three scales:
        ``resolution_wire`` on the torus surface, ``resolution_near`` inside the
        near-field box ``|x|,|y| ≤ near_half_width``,
        ``−near_depth ≤ z ≤ near_height`` (which must contain the coil and the
        skin layer), and ``resolution_far`` everywhere else.

        Parameters
        ----------
        loop_radius, wire_radius, liftoff : float
            Loop major radius ``a``, wire minor radius, and the height ``h`` of
            the loop plane above the conductor surface [m].  The wire must clear
            the interface: ``liftoff > wire_radius``.
        box_half_width : float
            Half side ``W`` of the cubic truncation box [m].  This is the knob
            the `MAT-6` step-2a probe sweeps: ``ΔZ`` is extracted from a
            loaded/free *difference*, so the coil self-impedance cancels, but
            the PEC wall still images the induced currents.
        resolution_wire, resolution_near, resolution_far : float
            Target mesh sizes [m] for the three zones described above.
        near_half_width, near_depth, near_height : float
            Extent of the fine near-field box [m]: lateral half-width, depth
            below the interface, height above it.
        """
        if liftoff <= wire_radius:
            raise ValueError(
                f"liftoff={liftoff!r} must exceed wire_radius={wire_radius!r} so the "
                "wire clears the half-space interface"
            )
        if box_half_width <= loop_radius + wire_radius:
            raise ValueError(
                f"box_half_width={box_half_width!r} must exceed the loop outer radius "
                f"{loop_radius + wire_radius!r}"
            )
        for name, value in (
            ("resolution_wire", resolution_wire),
            ("resolution_near", resolution_near),
            ("resolution_far", resolution_far),
        ):
            if value <= 0.0:
                raise ValueError(f"{name} must be positive, got {value!r}")

        W = float(box_half_width)

        if comm.rank == rank:
            gmsh.initialize()
            gmsh.model.add("loop_over_half_space")

            wire = gmsh.model.occ.addTorus(0, 0, liftoff, loop_radius, wire_radius)
            air = gmsh.model.occ.addBox(-W, -W, 0.0, 2 * W, 2 * W, W)
            slab = gmsh.model.occ.addBox(-W, -W, -W, 2 * W, 2 * W, W)

            gmsh.model.occ.fragment([(3, air), (3, slab)], [(3, wire)])
            gmsh.model.occ.synchronize()

            # Identify the three volumes by mass and centroid rather than by the
            # tags fragment happens to hand back: the wire is orders of magnitude
            # the smallest, and air/slab are told apart by the sign of z̄.
            wire_volume = None
            air_volume = None
            slab_volume = None
            torus_mass = 2.0 * np.pi**2 * loop_radius * wire_radius**2
            for dim, tag in gmsh.model.getEntities(dim=3):
                mass = gmsh.model.occ.getMass(dim, tag)
                _, _, zc = gmsh.model.occ.getCenterOfMass(dim, tag)
                if mass < 10.0 * torus_mass:
                    wire_volume = tag
                elif zc > 0.0:
                    air_volume = tag
                else:
                    slab_volume = tag

            if wire_volume is None or air_volume is None or slab_volume is None:
                raise RuntimeError(
                    "loop_over_half_space_domain: fragment did not produce the expected "
                    f"wire/air/slab volumes (got {gmsh.model.getEntities(dim=3)})"
                )

            gmsh.model.addPhysicalGroup(3, [wire_volume], tag=1)
            gmsh.model.setPhysicalName(3, 1, "wire")
            gmsh.model.addPhysicalGroup(3, [air_volume], tag=2)
            gmsh.model.setPhysicalName(3, 2, "air")
            gmsh.model.addPhysicalGroup(3, [slab_volume], tag=3)
            gmsh.model.setPhysicalName(3, 3, "slab")

            # `GEO-12`: 1e-9 -> 1e-6, exactly `GEO-10`'s fix and for the same
            # measured reason. gmsh inflates an OCC entity's bounding box by its
            # geometric tolerance — 1.000e-07 here as on `two_torus_domain`
            # (`20260806T140325Z_GEO-11-probe.log`) — so `tol = 1e-9` sat 100x
            # *below* the padding, accepted 0 of 12 surfaces, and the group was
            # never declared (known-issues 12). The nearest interior face sits at
            # 9.000e-02, five orders above 1e-6, so the interior-face protection
            # the tight test existed for is intact.
            tol = 1e-6
            boundary_surfaces = []
            for dim, surf in gmsh.model.getEntities(dim=2):
                x0, y0, z0, x1, y1, z1 = gmsh.model.getBoundingBox(dim, surf)
                on_wall = (
                    abs(x0 + W) < tol
                    and abs(x1 + W) < tol
                    or abs(x0 - W) < tol
                    and abs(x1 - W) < tol
                    or abs(y0 + W) < tol
                    and abs(y1 + W) < tol
                    or abs(y0 - W) < tol
                    and abs(y1 - W) < tol
                    or abs(z0 + W) < tol
                    and abs(z1 + W) < tol
                    or abs(z0 - W) < tol
                    and abs(z1 - W) < tol
                )
                if on_wall:
                    boundary_surfaces.append(surf)
            if boundary_surfaces:
                gmsh.model.addPhysicalGroup(2, boundary_surfaces, tag=1)
                gmsh.model.setPhysicalName(2, 1, "outer_boundary")

            wire_surfaces = [
                surf
                for _, surf in gmsh.model.getBoundary(
                    [(3, wire_volume)], oriented=False, recursive=False
                )
            ]

            dist = gmsh.model.mesh.field.add("Distance")
            gmsh.model.mesh.field.setNumbers(dist, "SurfacesList", wire_surfaces)
            gmsh.model.mesh.field.setNumber(dist, "Sampling", 200)

            wire_thr = gmsh.model.mesh.field.add("Threshold")
            gmsh.model.mesh.field.setNumber(wire_thr, "InField", dist)
            gmsh.model.mesh.field.setNumber(wire_thr, "SizeMin", resolution_wire)
            gmsh.model.mesh.field.setNumber(wire_thr, "SizeMax", resolution_far)
            gmsh.model.mesh.field.setNumber(wire_thr, "DistMin", wire_radius)
            gmsh.model.mesh.field.setNumber(wire_thr, "DistMax", W)

            near = gmsh.model.mesh.field.add("Box")
            gmsh.model.mesh.field.setNumber(near, "VIn", resolution_near)
            gmsh.model.mesh.field.setNumber(near, "VOut", resolution_far)
            gmsh.model.mesh.field.setNumber(near, "XMin", -near_half_width)
            gmsh.model.mesh.field.setNumber(near, "XMax", near_half_width)
            gmsh.model.mesh.field.setNumber(near, "YMin", -near_half_width)
            gmsh.model.mesh.field.setNumber(near, "YMax", near_half_width)
            gmsh.model.mesh.field.setNumber(near, "ZMin", -near_depth)
            gmsh.model.mesh.field.setNumber(near, "ZMax", near_height)
            # Blend over one far-cell so the near/far jump is not a mesh shock.
            gmsh.model.mesh.field.setNumber(near, "Thickness", resolution_far)

            combined = gmsh.model.mesh.field.add("Min")
            gmsh.model.mesh.field.setNumbers(combined, "FieldsList", [wire_thr, near])
            gmsh.model.mesh.field.setAsBackgroundMesh(combined)

            gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
            gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
            gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)

            gmsh.model.mesh.generate(3)

        mesh, cell_tags, facet_tags = _model_to_mesh(gmsh.model, comm, rank, gdim=3)

        if comm.rank == rank:
            gmsh.finalize()

        return mesh, cell_tags, facet_tags

    @staticmethod
    def sphere_in_box_domain(
        sphere_radius: float = 0.05,
        box_half_width: float = 0.10,
        resolution_sphere: float = 0.0125,
        resolution_far: float = 0.025,
        comm: MPI.Intracomm = MPI.COMM_WORLD,
        rank: int = 0,
        *,
        as_hole: bool = False,
    ) -> Tuple[dolfinx.mesh.Mesh, dolfinx.mesh.MeshTags, dolfinx.mesh.MeshTags]:
        """Dielectric sphere centred in a cubic air box (`TH-8`).

        Cell tags: ``1`` sphere, ``2`` air.  Facet tag ``1`` is the outer
        boundary of the cube.

        With ``as_hole=True`` (`TH-15` step 1) the sphere is *cut* from the box
        instead of fragmented into it: the meshed cell set is air only (tag
        ``2``, no tag ``1``), and the cavity wall arrives as facet tag ``2``.
        That surface group is built from the **meshed** volume's boundary, never
        from the cut tool — `TH-15` step 0 measured that a retained tool's faces
        are not the cavity's, and tagging them makes ``_model_to_mesh`` abort
        with ``nodes not attached to any tet`` / ``Invalid rank``.  ``False``
        (the default) leaves the `TH-8` generator bit-identical.

        Sizing is graded through a ``Ball`` field rather than a single global
        ``setSize``: the interior of the sphere is where the gate measures, and
        an unsigned ``Distance`` field from the sphere *surface* would coarsen
        towards the centre, which is exactly the region that must stay
        resolved.  The far zone only has to carry a smooth exterior field, so
        ``resolution_far`` may be several times ``resolution_sphere`` (the
        tight-padding/uniform-size pattern that cost 20% on Helmholtz is
        avoided by keeping the box independent of the sphere resolution).

        Parameters
        ----------
        sphere_radius : float
            Sphere radius ``R`` [m].
        box_half_width : float
            Half side ``W`` of the cubic truncation box [m].  `TH-8` imposes the
            exact exterior (uniform + dipole) field on this wall, so ``W`` sets
            how much air is discretised, not how wrong the truncation is.
        resolution_sphere, resolution_far : float
            Target mesh sizes [m] inside ``1.2 R`` and in the remaining air.
        """
        R = float(sphere_radius)
        W = float(box_half_width)
        if R <= 0.0:
            raise ValueError(f"sphere_radius must be positive, got {sphere_radius!r}")
        if W <= R:
            raise ValueError(
                f"box_half_width={box_half_width!r} must exceed sphere_radius={sphere_radius!r}"
            )
        for name, value in (
            ("resolution_sphere", resolution_sphere),
            ("resolution_far", resolution_far),
        ):
            if value <= 0.0:
                raise ValueError(f"{name} must be positive, got {value!r}")

        if comm.rank == rank:
            gmsh.initialize()
            gmsh.model.add("sphere_in_box")

            sphere = gmsh.model.occ.addSphere(0.0, 0.0, 0.0, R)
            box = gmsh.model.occ.addBox(-W, -W, -W, 2 * W, 2 * W, 2 * W)

            # `GEO-12`: 1e-9 -> 1e-6, same measured defect as above and as
            # `GEO-10`. The OCC bounding-box padding is 1.000e-07, so the old
            # tolerance accepted 0 of 7 surfaces and never declared the group
            # (known-issues 12); the nearest interior face is at 1.500e-01.
            tol = 1e-6

            def _on_wall(dim: int, surf: int) -> bool:
                x0, y0, z0, x1, y1, z1 = gmsh.model.getBoundingBox(dim, surf)
                return (
                    abs(x0 + W) < tol
                    and abs(x1 + W) < tol
                    or abs(x0 - W) < tol
                    and abs(x1 - W) < tol
                    or abs(y0 + W) < tol
                    and abs(y1 + W) < tol
                    or abs(y0 - W) < tol
                    and abs(y1 - W) < tol
                    or abs(z0 + W) < tol
                    and abs(z1 + W) < tol
                    or abs(z0 - W) < tol
                    and abs(z1 - W) < tol
                )

            if as_hole:
                cut_out, _ = gmsh.model.occ.cut([(3, box)], [(3, sphere)])
                gmsh.model.occ.synchronize()

                volumes = gmsh.model.getEntities(dim=3)
                if len(volumes) != 1:
                    raise RuntimeError(
                        "sphere_in_box_domain(as_hole=True): cut did not leave a "
                        f"single air volume (got {volumes}, cut returned {cut_out})"
                    )
                air_volume = volumes[0][1]

                gmsh.model.addPhysicalGroup(3, [air_volume], tag=2)
                gmsh.model.setPhysicalName(3, 2, "air")

                # The surface groups come from the *meshed* volume's boundary
                # (`TH-15` step 0's rule): every face of the air volume is either
                # an outer wall or the cavity wall, and both are faces of tets.
                boundary_surfaces = []
                cavity_surfaces = []
                for dim, surf in gmsh.model.getBoundary(
                    [(3, air_volume)], oriented=False, recursive=False
                ):
                    if _on_wall(dim, surf):
                        boundary_surfaces.append(surf)
                    else:
                        cavity_surfaces.append(surf)
                if not cavity_surfaces:
                    raise RuntimeError(
                        "sphere_in_box_domain(as_hole=True): the air volume has no "
                        "non-wall boundary face, so the cavity surface is missing"
                    )
                gmsh.model.addPhysicalGroup(2, cavity_surfaces, tag=2)
                gmsh.model.setPhysicalName(2, 2, "cavity_surface")
            else:
                gmsh.model.occ.fragment([(3, box)], [(3, sphere)])
                gmsh.model.occ.synchronize()

                # Identify by mass, not by the tags fragment hands back.
                sphere_mass = 4.0 / 3.0 * np.pi * R**3
                sphere_volume = None
                air_volume = None
                for dim, tag in gmsh.model.getEntities(dim=3):
                    mass = gmsh.model.occ.getMass(dim, tag)
                    if abs(mass - sphere_mass) < 0.05 * sphere_mass:
                        sphere_volume = tag
                    else:
                        air_volume = tag

                if sphere_volume is None or air_volume is None:
                    raise RuntimeError(
                        "sphere_in_box_domain: fragment did not produce the expected "
                        f"sphere/air volumes (got {gmsh.model.getEntities(dim=3)})"
                    )

                gmsh.model.addPhysicalGroup(3, [sphere_volume], tag=1)
                gmsh.model.setPhysicalName(3, 1, "sphere")
                gmsh.model.addPhysicalGroup(3, [air_volume], tag=2)
                gmsh.model.setPhysicalName(3, 2, "air")

                boundary_surfaces = [
                    surf for dim, surf in gmsh.model.getEntities(dim=2) if _on_wall(dim, surf)
                ]

            if boundary_surfaces:
                gmsh.model.addPhysicalGroup(2, boundary_surfaces, tag=1)
                gmsh.model.setPhysicalName(2, 1, "outer_boundary")

            ball = gmsh.model.mesh.field.add("Ball")
            gmsh.model.mesh.field.setNumber(ball, "XCenter", 0.0)
            gmsh.model.mesh.field.setNumber(ball, "YCenter", 0.0)
            gmsh.model.mesh.field.setNumber(ball, "ZCenter", 0.0)
            gmsh.model.mesh.field.setNumber(ball, "Radius", 1.2 * R)
            gmsh.model.mesh.field.setNumber(ball, "VIn", resolution_sphere)
            gmsh.model.mesh.field.setNumber(ball, "VOut", resolution_far)
            # Blend over one far-cell so the near/far jump is not a mesh shock.
            gmsh.model.mesh.field.setNumber(ball, "Thickness", resolution_far)
            gmsh.model.mesh.field.setAsBackgroundMesh(ball)

            gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
            gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
            gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)

            gmsh.model.mesh.generate(3)

        mesh, cell_tags, facet_tags = _model_to_mesh(gmsh.model, comm, rank, gdim=3)

        if comm.rank == rank:
            gmsh.finalize()

        return mesh, cell_tags, facet_tags

    @staticmethod
    def spheroid_in_box_domain(
        semi_major: float,
        semi_minor: float,
        box_half_extents: Tuple[float, float, float],
        resolution_tip: float,
        resolution_surface: float,
        resolution_far: float,
        comm: MPI.Intracomm = MPI.COMM_WORLD,
        rank: int = 0,
        *,
        as_hole: bool = True,
        oblate: bool = False,
        oblate_rim: str = "balls",
        oblate_netgen: bool = True,
    ) -> Tuple[dolfinx.mesh.Mesh, dolfinx.mesh.MeshTags, dolfinx.mesh.MeshTags]:
        """PEC prolate spheroid cut from a rectangular air box (`IMP-2` step 1,
        §9 item 144) — `sphere_in_box_domain`'s ``as_hole`` route with the unit
        sphere dilated to semi-axes (a, a, c), c along z.

        Cell tag ``2`` air (the only volume); facet tag ``1`` the outer walls,
        facet tag ``2`` the cavity wall, built from the **meshed** volume's
        boundary (`TH-15` step 0's rule).  Sizing: ``Min`` of two ``Ball``
        fields at the tips (0, 0, ±c) (radius 3ρ_t, ρ_t = a²/c, value
        ``resolution_tip``) and a Distance→Threshold from the cavity
        (``resolution_surface`` at the wall, ``resolution_far`` beyond 2a).
        Only ``as_hole=True`` is implemented.  Prints (on ``rank``) the air
        volume's OCC mass against 8·bx·by·bz − 4πa²c/3 and the minimum SICN
        element quality.

        ``oblate=True`` (`IMP-6` step 1, §9 item 404; additive, the default
        path is unchanged): the **oblate** spheroid with ``semi_major`` the
        equatorial radius (x, y) and ``semi_minor`` along z.  The unit sphere
        is rotated 45° about z before dilation so the OCC seam stays off the
        rim points (±a, 0, 0); two ``Ball`` fields at (±a, 0, 0) (radius 3ρ_r,
        ρ_r = c²/a, value ``resolution_tip``) plus a rim-tube ``MathEval`` →
        ``Threshold`` (SizeMin ρ_r within ρ_r of the rim circle, SizeMax
        ``resolution_far`` beyond 4ρ_r) join the ``Min``.

        ``oblate_rim`` / ``oblate_netgen`` (`IMP-6` step 2 attempt 2, §9 item
        426; rule (c)) are forwarded to the oblate builder as ``rim`` /
        ``netgen``; their defaults reproduce step 1's mesh, and a non-default
        value without ``oblate=True`` raises ``ValueError``.
        """
        if not as_hole:
            raise NotImplementedError("spheroid_in_box_domain: only as_hole=True is implemented")
        if not oblate and (oblate_rim != "balls" or oblate_netgen is not True):
            raise ValueError("oblate_rim / oblate_netgen need oblate=True")
        if oblate:
            return MeshGenerator._oblate_spheroid_in_box_domain(
                semi_major,
                semi_minor,
                box_half_extents,
                resolution_tip,
                resolution_surface,
                resolution_far,
                comm,
                rank,
                rim=oblate_rim,
                netgen=oblate_netgen,
            )
        c = float(semi_major)
        a = float(semi_minor)
        bx, by, bz = (float(v) for v in box_half_extents)
        if not c > a > 0.0:
            raise ValueError("spheroid_in_box_domain needs semi_major > semi_minor > 0")
        if not (bx > a and by > a and bz > c):
            raise ValueError("box_half_extents must enclose the spheroid")
        for name, value in (
            ("resolution_tip", resolution_tip),
            ("resolution_surface", resolution_surface),
            ("resolution_far", resolution_far),
        ):
            if value <= 0.0:
                raise ValueError(f"{name} must be positive, got {value!r}")
        rho_t = a * a / c

        if comm.rank == rank:
            gmsh.initialize()
            gmsh.model.add("spheroid_in_box")
            body = gmsh.model.occ.addSphere(0.0, 0.0, 0.0, 1.0)
            gmsh.model.occ.dilate([(3, body)], 0.0, 0.0, 0.0, a, a, c)
            box = gmsh.model.occ.addBox(-bx, -by, -bz, 2 * bx, 2 * by, 2 * bz)
            # Absolute, as `sphere_in_box_domain` (`GEO-12`): the OCC bounding
            # box pads by 1e-7, so a scale-relative 1e-6·6 mm tolerance missed
            # every wall (first W0, 2026-10-01).
            tol = 1e-6

            def _on_wall(dim: int, surf: int) -> bool:
                x0, y0, z0, x1, y1, z1 = gmsh.model.getBoundingBox(dim, surf)
                return (
                    abs(x0 + bx) < tol
                    and abs(x1 + bx) < tol
                    or abs(x0 - bx) < tol
                    and abs(x1 - bx) < tol
                    or abs(y0 + by) < tol
                    and abs(y1 + by) < tol
                    or abs(y0 - by) < tol
                    and abs(y1 - by) < tol
                    or abs(z0 + bz) < tol
                    and abs(z1 + bz) < tol
                    or abs(z0 - bz) < tol
                    and abs(z1 - bz) < tol
                )

            cut_out, _ = gmsh.model.occ.cut([(3, box)], [(3, body)])
            gmsh.model.occ.synchronize()
            volumes = gmsh.model.getEntities(dim=3)
            if len(volumes) != 1:
                raise RuntimeError(
                    "spheroid_in_box_domain: cut did not leave a single air volume "
                    f"(got {volumes}, cut returned {cut_out})"
                )
            air_volume = volumes[0][1]
            mass = gmsh.model.occ.getMass(3, air_volume)
            mass_ref = 8.0 * bx * by * bz - 4.0 * np.pi * a * a * c / 3.0
            print(
                f"[spheroid_in_box_domain] a = {a:g}, c = {c:g}: air OCC mass "
                f"{mass:.9e} vs 8 bx by bz - 4 pi a^2 c/3 = {mass_ref:.9e} "
                f"(rel {abs(mass / mass_ref - 1.0):.3e})",
                flush=True,
            )
            gmsh.model.addPhysicalGroup(3, [air_volume], tag=2)
            gmsh.model.setPhysicalName(3, 2, "air")
            boundary_surfaces = []
            cavity_surfaces = []
            for dim, surf in gmsh.model.getBoundary(
                [(3, air_volume)], oriented=False, recursive=False
            ):
                if _on_wall(dim, surf):
                    boundary_surfaces.append(surf)
                else:
                    cavity_surfaces.append(surf)
            if not cavity_surfaces or len(boundary_surfaces) != 6:
                raise RuntimeError(
                    "spheroid_in_box_domain: expected 6 wall faces and a cavity "
                    f"(walls {boundary_surfaces}, cavity {cavity_surfaces})"
                )
            gmsh.model.addPhysicalGroup(2, boundary_surfaces, tag=1)
            gmsh.model.setPhysicalName(2, 1, "outer_boundary")
            gmsh.model.addPhysicalGroup(2, cavity_surfaces, tag=2)
            gmsh.model.setPhysicalName(2, 2, "cavity_surface")

            fields = []
            for zc in (c, -c):
                ball = gmsh.model.mesh.field.add("Ball")
                gmsh.model.mesh.field.setNumber(ball, "XCenter", 0.0)
                gmsh.model.mesh.field.setNumber(ball, "YCenter", 0.0)
                gmsh.model.mesh.field.setNumber(ball, "ZCenter", zc)
                gmsh.model.mesh.field.setNumber(ball, "Radius", 3.0 * rho_t)
                gmsh.model.mesh.field.setNumber(ball, "VIn", resolution_tip)
                gmsh.model.mesh.field.setNumber(ball, "VOut", resolution_far)
                fields.append(ball)
            dist = gmsh.model.mesh.field.add("Distance")
            gmsh.model.mesh.field.setNumbers(dist, "SurfacesList", cavity_surfaces)
            gmsh.model.mesh.field.setNumber(dist, "Sampling", 200)
            thr = gmsh.model.mesh.field.add("Threshold")
            gmsh.model.mesh.field.setNumber(thr, "InField", dist)
            gmsh.model.mesh.field.setNumber(thr, "SizeMin", resolution_surface)
            gmsh.model.mesh.field.setNumber(thr, "SizeMax", resolution_far)
            gmsh.model.mesh.field.setNumber(thr, "DistMin", 0.0)
            gmsh.model.mesh.field.setNumber(thr, "DistMax", 2.0 * a)
            fields.append(thr)
            fmin = gmsh.model.mesh.field.add("Min")
            gmsh.model.mesh.field.setNumbers(fmin, "FieldsList", fields)
            gmsh.model.mesh.field.setAsBackgroundMesh(fmin)

            gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
            gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
            gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
            gmsh.model.mesh.generate(3)

            etypes, etags, _ = gmsh.model.mesh.getElements(3)
            qmin = float("nan")
            for et, tags in zip(etypes, etags):
                q = gmsh.model.mesh.getElementQualities(tags, "minSICN")
                qmin = float(np.min(q)) if np.isnan(qmin) else min(qmin, float(np.min(q)))
            print(
                f"[spheroid_in_box_domain] minimum SICN quality {qmin:.4f}",
                flush=True,
            )

        mesh, cell_tags, facet_tags = _model_to_mesh(gmsh.model, comm, rank, gdim=3)
        if comm.rank == rank:
            gmsh.finalize()
        return mesh, cell_tags, facet_tags

    @staticmethod
    def _oblate_spheroid_in_box_domain(
        semi_major: float,
        semi_minor: float,
        box_half_extents: Tuple[float, float, float],
        resolution_tip: float,
        resolution_surface: float,
        resolution_far: float,
        comm: MPI.Intracomm = MPI.COMM_WORLD,
        rank: int = 0,
        *,
        rim: str = "balls",
        netgen: bool = True,
    ) -> Tuple[dolfinx.mesh.Mesh, dolfinx.mesh.MeshTags, dolfinx.mesh.MeshTags]:
        """``spheroid_in_box_domain(..., oblate=True)`` (`IMP-6` step 1).

        ``rim`` / ``netgen`` (`IMP-6` step 2 attempt 2, §9 item 426; rule (c),
        defaults reproduce the step-1 path): ``rim="edge"`` builds the body by
        revolving a half-ellipse so the equator is a model edge; ``rim="tube"``
        replaces each rim ``Ball`` by ``Max(Threshold(rim distance), Ball)`` so
        the size grows with distance from the rim circle; ``netgen=False``
        skips the Netgen optimisation pass.
        """
        a = float(semi_major)
        c = float(semi_minor)
        bx, by, bz = (float(v) for v in box_half_extents)
        if not a > c > 0.0:
            raise ValueError("spheroid_in_box_domain needs semi_major > semi_minor > 0")
        if not (bx > a and by > a and bz > c):
            raise ValueError("box_half_extents must enclose the spheroid")
        for name, value in (
            ("resolution_tip", resolution_tip),
            ("resolution_surface", resolution_surface),
            ("resolution_far", resolution_far),
        ):
            if value <= 0.0:
                raise ValueError(f"{name} must be positive, got {value!r}")
        if rim not in ("balls", "edge", "tube"):
            raise ValueError(f"rim must be 'balls', 'edge' or 'tube', got {rim!r}")
        rho_r = c * c / a

        if comm.rank == rank:
            gmsh.initialize()
            gmsh.model.add("oblate_spheroid_in_box")
            if rim == "edge":
                # IMP-6 step 2 attempt 2 (§9 item 426, variant (1a)): the body
                # by revolution, so the equator is a model edge shared by the
                # upper and lower faces and no surface triangle bridges it.
                occ = gmsh.model.occ
                p_ctr = occ.addPoint(0.0, 0.0, 0.0)
                p_top = occ.addPoint(0.0, 0.0, c)
                p_rim = occ.addPoint(a, 0.0, 0.0)
                p_bot = occ.addPoint(0.0, 0.0, -c)
                arc_up = occ.addEllipseArc(p_top, p_ctr, p_rim, p_rim)
                arc_lo = occ.addEllipseArc(p_rim, p_ctr, p_rim, p_bot)
                axis = occ.addLine(p_bot, p_top)
                prof = occ.addPlaneSurface([occ.addCurveLoop([arc_up, arc_lo, axis])])
                rev = occ.revolve([(2, prof)], 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0 * np.pi)
                rev_vols = [t for d, t in rev if d == 3]
                if len(rev_vols) != 1:
                    raise RuntimeError(
                        f"spheroid_in_box_domain(oblate, rim='edge'): revolve gave {rev}"
                    )
                body = rev_vols[0]
            else:
                body = gmsh.model.occ.addSphere(0.0, 0.0, 0.0, 1.0)
            # Seam off the drive's rim points (±a, 0, 0).
            gmsh.model.occ.rotate([(3, body)], 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, np.pi / 4.0)
            if rim != "edge":
                gmsh.model.occ.dilate([(3, body)], 0.0, 0.0, 0.0, a, a, c)
            box = gmsh.model.occ.addBox(-bx, -by, -bz, 2 * bx, 2 * by, 2 * bz)
            tol = 1e-6  # absolute, as the prolate path

            def _on_wall(dim: int, surf: int) -> bool:
                x0, y0, z0, x1, y1, z1 = gmsh.model.getBoundingBox(dim, surf)
                return (
                    abs(x0 + bx) < tol
                    and abs(x1 + bx) < tol
                    or abs(x0 - bx) < tol
                    and abs(x1 - bx) < tol
                    or abs(y0 + by) < tol
                    and abs(y1 + by) < tol
                    or abs(y0 - by) < tol
                    and abs(y1 - by) < tol
                    or abs(z0 + bz) < tol
                    and abs(z1 + bz) < tol
                    or abs(z0 - bz) < tol
                    and abs(z1 - bz) < tol
                )

            cut_out, _ = gmsh.model.occ.cut([(3, box)], [(3, body)])
            gmsh.model.occ.synchronize()
            volumes = gmsh.model.getEntities(dim=3)
            if len(volumes) != 1:
                raise RuntimeError(
                    "spheroid_in_box_domain(oblate): cut did not leave a single air "
                    f"volume (got {volumes}, cut returned {cut_out})"
                )
            air_volume = volumes[0][1]
            if rim == "edge":
                # The revolve's profile face can survive the cut as a free
                # surface inside the removed body; drop anything not bounding air.
                keep = {
                    t
                    for _, t in gmsh.model.getBoundary(
                        [(3, air_volume)], oriented=False, recursive=False
                    )
                }
                stray = [(2, t) for _, t in gmsh.model.getEntities(2) if t not in keep]
                if stray:
                    gmsh.model.occ.remove(stray, recursive=True)
                    gmsh.model.occ.synchronize()
                print(
                    f"[spheroid_in_box_domain oblate] rim='edge': stray surfaces "
                    f"removed {stray}",
                    flush=True,
                )
            mass = gmsh.model.occ.getMass(3, air_volume)
            mass_ref = 8.0 * bx * by * bz - 4.0 * np.pi * a * a * c / 3.0
            print(
                f"[spheroid_in_box_domain oblate] a = {a:g}, c = {c:g}: air OCC mass "
                f"{mass:.9e} vs 8 bx by bz - 4 pi a^2 c/3 = {mass_ref:.9e} "
                f"(rel {abs(mass / mass_ref - 1.0):.3e})",
                flush=True,
            )
            gmsh.model.addPhysicalGroup(3, [air_volume], tag=2)
            gmsh.model.setPhysicalName(3, 2, "air")
            boundary_surfaces = []
            cavity_surfaces = []
            for dim, surf in gmsh.model.getBoundary(
                [(3, air_volume)], oriented=False, recursive=False
            ):
                if _on_wall(dim, surf):
                    boundary_surfaces.append(surf)
                else:
                    cavity_surfaces.append(surf)
            if not cavity_surfaces or len(boundary_surfaces) != 6:
                raise RuntimeError(
                    "spheroid_in_box_domain(oblate): expected 6 wall faces and a "
                    f"cavity (walls {boundary_surfaces}, cavity {cavity_surfaces})"
                )
            gmsh.model.addPhysicalGroup(2, boundary_surfaces, tag=1)
            gmsh.model.setPhysicalName(2, 1, "outer_boundary")
            gmsh.model.addPhysicalGroup(2, cavity_surfaces, tag=2)
            gmsh.model.setPhysicalName(2, 2, "cavity_surface")
            if rim == "edge":
                print(
                    f"[spheroid_in_box_domain oblate] rim='edge': cavity faces "
                    f"{len(cavity_surfaces)} (2 predicted)",
                    flush=True,
                )

            fields = []
            balls = []
            for xc in (a, -a):
                ball = gmsh.model.mesh.field.add("Ball")
                gmsh.model.mesh.field.setNumber(ball, "XCenter", xc)
                gmsh.model.mesh.field.setNumber(ball, "YCenter", 0.0)
                gmsh.model.mesh.field.setNumber(ball, "ZCenter", 0.0)
                gmsh.model.mesh.field.setNumber(ball, "Radius", 3.0 * rho_r)
                gmsh.model.mesh.field.setNumber(ball, "VIn", resolution_tip)
                gmsh.model.mesh.field.setNumber(ball, "VOut", resolution_far)
                balls.append(ball)
            if rim != "tube":
                fields.extend(balls)
            rim_expr = f"Sqrt((Sqrt(x*x+y*y)-{a!r})^2+z*z)"
            print(f"[spheroid_in_box_domain oblate] rim MathEval {rim_expr}", flush=True)
            rim_dist = gmsh.model.mesh.field.add("MathEval")
            gmsh.model.mesh.field.setString(rim_dist, "F", rim_expr)
            if rim == "tube":
                # Variant (1b): inside the old Ball the size grows with
                # distance from the rim circle; outside it Max gives far.
                for ball in balls:
                    graded = gmsh.model.mesh.field.add("Threshold")
                    gmsh.model.mesh.field.setNumber(graded, "InField", rim_dist)
                    gmsh.model.mesh.field.setNumber(graded, "SizeMin", resolution_tip)
                    gmsh.model.mesh.field.setNumber(graded, "SizeMax", resolution_far)
                    gmsh.model.mesh.field.setNumber(graded, "DistMin", 0.0)
                    gmsh.model.mesh.field.setNumber(graded, "DistMax", 3.0 * rho_r)
                    fmax = gmsh.model.mesh.field.add("Max")
                    gmsh.model.mesh.field.setNumbers(fmax, "FieldsList", [graded, ball])
                    fields.append(fmax)
            rim_thr = gmsh.model.mesh.field.add("Threshold")
            gmsh.model.mesh.field.setNumber(rim_thr, "InField", rim_dist)
            gmsh.model.mesh.field.setNumber(rim_thr, "SizeMin", rho_r)
            gmsh.model.mesh.field.setNumber(rim_thr, "SizeMax", resolution_far)
            gmsh.model.mesh.field.setNumber(rim_thr, "DistMin", rho_r)
            gmsh.model.mesh.field.setNumber(rim_thr, "DistMax", 4.0 * rho_r)
            fields.append(rim_thr)
            dist = gmsh.model.mesh.field.add("Distance")
            gmsh.model.mesh.field.setNumbers(dist, "SurfacesList", cavity_surfaces)
            gmsh.model.mesh.field.setNumber(dist, "Sampling", 200)
            thr = gmsh.model.mesh.field.add("Threshold")
            gmsh.model.mesh.field.setNumber(thr, "InField", dist)
            gmsh.model.mesh.field.setNumber(thr, "SizeMin", resolution_surface)
            gmsh.model.mesh.field.setNumber(thr, "SizeMax", resolution_far)
            gmsh.model.mesh.field.setNumber(thr, "DistMin", 0.0)
            gmsh.model.mesh.field.setNumber(thr, "DistMax", 2.0 * c)
            fields.append(thr)
            fmin = gmsh.model.mesh.field.add("Min")
            gmsh.model.mesh.field.setNumbers(fmin, "FieldsList", fields)
            gmsh.model.mesh.field.setAsBackgroundMesh(fmin)

            gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
            gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
            gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
            gmsh.model.mesh.generate(3)
            # IMP-6 W0 (20261012T183906Z): aspect 2 read min SICN 0.0247 with
            # the default optimiser; a Netgen pass is added on this path only.
            if netgen:
                gmsh.model.mesh.optimize("Netgen")

            etypes, etags, _ = gmsh.model.mesh.getElements(3)
            qmin = float("nan")
            worst = None
            for et, tags in zip(etypes, etags):
                q = gmsh.model.mesh.getElementQualities(tags, "minSICN")
                i = int(np.argmin(q))
                if np.isnan(qmin) or float(q[i]) < qmin:
                    qmin = float(q[i])
                    worst = int(tags[i])
            loc = ""
            if worst is not None:
                _, nodes, _, _ = gmsh.model.mesh.getElement(worst)
                xyz = np.array([gmsh.model.mesh.getNode(int(n))[0] for n in nodes])
                loc = f" at centroid {np.array2string(xyz.mean(axis=0), precision=5)}"
            print(
                f"[spheroid_in_box_domain oblate] minimum SICN quality {qmin:.4f}{loc}",
                flush=True,
            )

        mesh, cell_tags, facet_tags = _model_to_mesh(gmsh.model, comm, rank, gdim=3)
        if comm.rank == rank:
            gmsh.finalize()
        return mesh, cell_tags, facet_tags

    @staticmethod
    def ball_domain(
        radius: float = 0.1,
        resolution: float = 0.01,
        comm: MPI.Intracomm = MPI.COMM_WORLD,
        rank: int = 0,
    ) -> Tuple[dolfinx.mesh.Mesh, dolfinx.mesh.MeshTags, dolfinx.mesh.MeshTags]:
        """Homogeneous solid ball centred at the origin (`MAT-9` step 1).

        One OCC ``addSphere`` meshed at a uniform target size ``resolution``
        (``Mesh.MeshSizeMin`` = ``Mesh.MeshSizeMax``, no curvature or point
        sizing), so a ratio-2 ``resolution`` ladder is a ratio-2 ``h`` ladder.
        Cell tag ``1`` is the ball; facet tag ``1`` is its exterior surface —
        the only surface, so nothing interior is tagged (known-issues 9).
        """
        R = float(radius)
        if R <= 0.0:
            raise ValueError(f"radius must be positive, got {radius!r}")
        if resolution <= 0.0:
            raise ValueError(f"resolution must be positive, got {resolution!r}")

        if comm.rank == rank:
            gmsh.initialize()
            gmsh.option.setNumber("General.Terminal", 0)
            gmsh.model.add("ball")
            ball = gmsh.model.occ.addSphere(0.0, 0.0, 0.0, R)
            gmsh.model.occ.synchronize()

            gmsh.model.addPhysicalGroup(3, [ball], tag=1)
            gmsh.model.setPhysicalName(3, 1, "ball")
            surfaces = [
                surf
                for _, surf in gmsh.model.getBoundary([(3, ball)], oriented=False, recursive=False)
            ]
            gmsh.model.addPhysicalGroup(2, surfaces, tag=1)
            gmsh.model.setPhysicalName(2, 1, "ball_surface")

            gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
            gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
            gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
            gmsh.option.setNumber("Mesh.MeshSizeMin", float(resolution))
            gmsh.option.setNumber("Mesh.MeshSizeMax", float(resolution))
            gmsh.model.mesh.generate(3)

        mesh, cell_tags, facet_tags = _model_to_mesh(gmsh.model, comm, rank, gdim=3)

        if comm.rank == rank:
            gmsh.finalize()

        return mesh, cell_tags, facet_tags

    @staticmethod
    def two_layer_ball_domain(
        inner_radius: float = 0.1,
        outer_radius: float = 0.11,
        resolution: float = 0.01,
        comm: MPI.Intracomm = MPI.COMM_WORLD,
        rank: int = 0,
        geometry_order: int = 1,
    ) -> Tuple[dolfinx.mesh.Mesh, dolfinx.mesh.MeshTags, dolfinx.mesh.MeshTags]:
        """Conforming core + shell ball centred at the origin (`MAT-18` step 2(a)).

        ``geometry_order`` (`IMP-7` step 3, additive): 1 (default) is the
        affine path, unchanged; 2 calls ``gmsh.model.mesh.setOrder(2)`` after
        ``generate(3)`` so the coordinate element is degree 2 (10-node tets
        whose edge midpoints lie on the r = a / r = b spheres).  Anything else
        raises.

        Two OCC spheres (radii ``inner_radius`` = a < ``outer_radius`` = b)
        fragmented so the r = a interface is shared, meshed at a uniform
        target size ``resolution`` exactly as :meth:`ball_domain`.  The two
        volumes are identified by ``occ.getMass`` (the core is 4πa³/3), never
        by tag order.  Cell tag ``1`` is the core, ``2`` the shell; facet tag
        ``1`` is the r = b surface only.  No physical group is placed on the
        r = a interface (interior facet groups hang ``model_to_mesh`` at
        ``-n 2``, known-issues 9); use ``_interface_facet_tags`` if needed.
        """
        a = float(inner_radius)
        b = float(outer_radius)
        if not (0.0 < a < b):
            raise ValueError(
                f"need 0 < inner_radius < outer_radius, got {inner_radius!r}, " f"{outer_radius!r}"
            )
        if resolution <= 0.0:
            raise ValueError(f"resolution must be positive, got {resolution!r}")
        if geometry_order not in (1, 2):
            raise ValueError(f"geometry_order must be 1 or 2, got {geometry_order!r}")

        if comm.rank == rank:
            gmsh.initialize()
            gmsh.option.setNumber("General.Terminal", 0)
            gmsh.model.add("two_layer_ball")
            occ = gmsh.model.occ
            outer = occ.addSphere(0.0, 0.0, 0.0, b)
            inner = occ.addSphere(0.0, 0.0, 0.0, a)
            occ.fragment([(3, outer)], [(3, inner)])
            occ.synchronize()

            vols = [tag for dim, tag in gmsh.model.getEntities(3)]
            if len(vols) != 2:
                raise RuntimeError(f"expected 2 volumes after fragment, got {vols}")
            v_core = 4.0 / 3.0 * np.pi * a**3
            masses = {v: occ.getMass(3, v) for v in vols}
            core = min(vols, key=lambda v: abs(masses[v] - v_core))
            shell = [v for v in vols if v != core][0]

            # Outer surface: the boundary surface whose area is ≈ 4πb².
            area_b = 4.0 * np.pi * b**2
            shell_surfs = [
                s for _, s in gmsh.model.getBoundary([(3, shell)], oriented=False, recursive=False)
            ]
            outer_surfs = [
                s for s in shell_surfs if abs(occ.getMass(2, s) - area_b) <= 1e-6 * area_b
            ]
            if len(outer_surfs) != 1:
                raise RuntimeError(f"could not identify the r = b surface among {shell_surfs}")

            gmsh.model.addPhysicalGroup(3, [core], tag=1)
            gmsh.model.setPhysicalName(3, 1, "core")
            gmsh.model.addPhysicalGroup(3, [shell], tag=2)
            gmsh.model.setPhysicalName(3, 2, "shell")
            gmsh.model.addPhysicalGroup(2, outer_surfs, tag=1)
            gmsh.model.setPhysicalName(2, 1, "outer_surface")

            gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
            gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
            gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
            gmsh.option.setNumber("Mesh.MeshSizeMin", float(resolution))
            gmsh.option.setNumber("Mesh.MeshSizeMax", float(resolution))
            gmsh.model.mesh.generate(3)
            if geometry_order == 2:
                gmsh.model.mesh.setOrder(2)

        mesh, cell_tags, facet_tags = _model_to_mesh(gmsh.model, comm, rank, gdim=3)

        if comm.rank == rank:
            gmsh.finalize()

        return mesh, cell_tags, facet_tags

    @staticmethod
    def coil_phantom_domain_sizing_diagnostics(
        *,
        coil_major_radius: float,
        coil_minor_radius: float,
        coil_separation: float,
        phantom_radius: float,
        phantom_height: float,
        air_padding: float,
        phantom_offset_xy: Tuple[float, float] = (0.0, 0.0),
        min_air_padding_ratio: float = 0.35,
    ) -> Dict[str, float | bool]:
        """Compute air-box sizing diagnostics for the coil+phantom fixture.

        The returned values are lightweight geometry heuristics (no meshing).
        They are used to flag undersized domains that can amplify boundary artifacts.

        Sizing rule (the intent the numbers encode, `GEO-4` step 1): the air box is
        **centred on the origin** (see `coil_phantom_domain`, which builds it from
        these extents), so its half-width must contain every object's largest |x|:

            half_width = max(coil_major + coil_minor, |offset| + phantom_radius)
                         + padding

        An off-centre phantom therefore enters through the second term of the max
        and only *governs* the box once it reaches outside the coil envelope. Below
        that it is masked by the coil, and the correct domain does not grow — the
        clearance to the wall shrinks instead, which is what
        `phantom_boundary_clearance_m` reports.
        """
        if min_air_padding_ratio <= 0.0:
            raise ValueError("min_air_padding_ratio must be > 0")
        if air_padding < 0.0:
            raise ValueError("air_padding must be >= 0")

        phantom_cx, phantom_cy = float(phantom_offset_xy[0]), float(phantom_offset_xy[1])
        phantom_offset_radius = np.sqrt(phantom_cx**2 + phantom_cy**2)

        radial_extent_without_padding = max(
            coil_major_radius + coil_minor_radius,
            phantom_offset_radius + phantom_radius,
        )
        z_offset = coil_separation / 2.0
        z_extent_without_padding = max(z_offset + coil_minor_radius, phantom_height / 2.0)

        reference_extent = max(radial_extent_without_padding, z_extent_without_padding)
        recommended_min_padding = min_air_padding_ratio * reference_extent
        effective_air_padding = max(air_padding, recommended_min_padding)

        # Clearance from the phantom's outermost radial point to the recommended
        # wall. By construction of the max above this is >= recommended_min_padding,
        # with equality exactly when the phantom governs the box.
        phantom_outer_radial_extent = phantom_offset_radius + phantom_radius
        recommended_domain_half_width = radial_extent_without_padding + recommended_min_padding
        phantom_boundary_clearance = recommended_domain_half_width - phantom_outer_radial_extent

        return {
            "phantom_offset_radius_m": float(phantom_offset_radius),
            "phantom_outer_radial_extent_m": float(phantom_outer_radial_extent),
            "phantom_boundary_clearance_m": float(phantom_boundary_clearance),
            "phantom_governs_radial_extent": bool(
                phantom_outer_radial_extent >= coil_major_radius + coil_minor_radius
            ),
            "radial_extent_without_padding_m": float(radial_extent_without_padding),
            "z_extent_without_padding_m": float(z_extent_without_padding),
            "reference_extent_m": float(reference_extent),
            "provided_air_padding_m": float(air_padding),
            "recommended_min_air_padding_m": float(recommended_min_padding),
            "effective_air_padding_m": float(effective_air_padding),
            "is_domain_undersized": bool(air_padding < recommended_min_padding),
            "recommended_domain_half_width_m": float(recommended_domain_half_width),
            "recommended_domain_half_height_m": float(
                z_extent_without_padding + recommended_min_padding
            ),
        }

    @staticmethod
    def two_torus_box_padding(tolerance_pp: float, *, quantum_m: float = 0.01) -> Dict[str, float]:
        """Air padding whose PEC-box term on the two-torus fixture is ≤ τ.

        `GEO-4` step 2.  The box term is ``B(W) = C·W^(−p)`` from `PORT-1`
        decision (4)'s fit ``deficit(W) = D∞ + C·W^(−p)`` through the step
        3b-xi rungs (``PORT1_PADDING_RUNGS_M`` / ``PORT1_PADDING_DEFICITS_PP``),
        recomputed here.  Returns the smallest padding on a ``quantum_m`` grid
        with ``|B| ≤ tolerance_pp``.

        Scope: two-torus only, an effective-range fit (fitted on 0.08–0.12 m),
        not a converged box-free deficit.  The coil+phantom 0.35·extent rule
        in `coil_phantom_domain_sizing_diagnostics` is unaffected.
        """
        if tolerance_pp <= 0.0:
            raise ValueError("tolerance_pp must be > 0")
        d_inf, p = _free_exponent_fit(PORT1_PADDING_RUNGS_M, PORT1_PADDING_DEFICITS_PP)
        if d_inf is None:
            raise ValueError("the PORT-1 padding rungs admit no power-law box term")
        c_pp = (PORT1_PADDING_DEFICITS_PP[0] - d_inf) * PORT1_PADDING_RUNGS_M[0] ** p
        padding, unrounded = _padding_for_box_term(c_pp, p, tolerance_pp, quantum_m)
        return {
            "padding_m": float(padding),
            "unrounded_padding_m": float(unrounded),
            "predicted_box_term_pp": float(abs(c_pp) * padding ** (-p)),
            "d_inf_pp": float(d_inf),
            "c_pp": float(c_pp),
            "p": float(p),
        }

    @staticmethod
    def coil_phantom_geometry_sanity_report(
        *,
        cell_tags: dolfinx.mesh.MeshTags,
        coil_major_radius: float,
        coil_minor_radius: float,
        coil_separation: float,
        phantom_radius: float,
        phantom_height: float,
        air_padding: float,
        phantom_offset_xy: Tuple[float, float] = (0.0, 0.0),
        comm: MPI.Intracomm = MPI.COMM_WORLD,
    ) -> Dict[str, object]:
        """Build a compact geometry/tag sanity report for the coil+phantom fixture.

        The report combines:
        - required tag counts (global)
        - expected analytic region-volume ratios (coil/phantom/air)
        - observed tag cell-count ratios
        - lightweight warnings for obviously suspicious setups
        """
        required_tags = {
            1: "coil_1",
            2: "coil_2",
            3: "phantom",
            4: "air",
        }

        counts = {}
        for tag in required_tags:
            local = int(np.count_nonzero(cell_tags.values == tag))
            counts[tag] = int(comm.allreduce(local, op=MPI.SUM))

        total_cells = int(sum(counts.values()))
        observed_cell_ratios = {
            tag: (counts[tag] / total_cells if total_cells > 0 else 0.0) for tag in required_tags
        }

        sizing = MeshGenerator.coil_phantom_domain_sizing_diagnostics(
            coil_major_radius=coil_major_radius,
            coil_minor_radius=coil_minor_radius,
            coil_separation=coil_separation,
            phantom_radius=phantom_radius,
            phantom_height=phantom_height,
            air_padding=air_padding,
            phantom_offset_xy=phantom_offset_xy,
        )

        effective_air_padding = sizing["effective_air_padding_m"]
        radial_extent = sizing["radial_extent_without_padding_m"]
        z_extent = sizing["z_extent_without_padding_m"]

        box_volume = (
            (2.0 * (radial_extent + effective_air_padding))
            * (2.0 * (radial_extent + effective_air_padding))
            * (2.0 * (z_extent + effective_air_padding))
        )
        coil_single_volume = 2.0 * np.pi**2 * coil_major_radius * (coil_minor_radius**2)
        phantom_volume = np.pi * (phantom_radius**2) * phantom_height
        occupied_volume = 2.0 * coil_single_volume + phantom_volume
        air_volume = max(0.0, box_volume - occupied_volume)

        expected_volumes = {
            1: float(coil_single_volume),
            2: float(coil_single_volume),
            3: float(phantom_volume),
            4: float(air_volume),
        }
        expected_total_volume = float(sum(expected_volumes.values()))
        expected_volume_ratios = {
            tag: (
                expected_volumes[tag] / expected_total_volume
                if expected_total_volume > 0.0
                else 0.0
            )
            for tag in required_tags
        }

        warnings: List[str] = []
        missing = [required_tags[tag] for tag, count in counts.items() if count <= 0]
        if missing:
            warnings.append(f"missing required tags: {', '.join(missing)}")

        if sizing["is_domain_undersized"]:
            warnings.append(
                "air padding below recommended minimum "
                f"({sizing['provided_air_padding_m']:.6e} < {sizing['recommended_min_air_padding_m']:.6e} m)"
            )

        for tag, name in required_tags.items():
            expected_ratio = expected_volume_ratios[tag]
            observed_ratio = observed_cell_ratios[tag]
            if expected_ratio > 0.0 and observed_ratio > 0.0:
                mismatch = max(observed_ratio / expected_ratio, expected_ratio / observed_ratio)
                if mismatch > 5.0:
                    warnings.append(
                        f"{name} observed/expected ratio mismatch is large "
                        f"(observed={observed_ratio:.3f}, expected={expected_ratio:.3f})"
                    )

        return {
            "required_tag_counts": {required_tags[tag]: counts[tag] for tag in required_tags},
            "expected_volume_ratios": {
                required_tags[tag]: expected_volume_ratios[tag] for tag in required_tags
            },
            "observed_cell_ratios": {
                required_tags[tag]: observed_cell_ratios[tag] for tag in required_tags
            },
            "effective_air_padding_m": float(effective_air_padding),
            "warnings": warnings,
            "ok": len(warnings) == 0,
        }

    @staticmethod
    def print_coil_phantom_geometry_sanity_report(
        *,
        report: Dict[str, object],
        prefix: str = "[coil-phantom-sanity] ",
        comm: MPI.Intracomm = MPI.COMM_WORLD,
    ) -> None:
        """Print a deterministic, compact sanity report on rank 0."""
        if comm.rank != 0:
            return

        print(f"{prefix}geometry sanity report:")

        counts = report["required_tag_counts"]
        print(f"{prefix}required tag counts:")
        for name in ("coil_1", "coil_2", "phantom", "air"):
            print(f"{prefix}  {name}: {counts[name]}")

        expected = report["expected_volume_ratios"]
        observed = report["observed_cell_ratios"]
        print(f"{prefix}volume ratio check (expected_volume_ratio vs observed_cell_ratio):")
        for name in ("coil_1", "coil_2", "phantom", "air"):
            print(f"{prefix}  {name}: expected={expected[name]:.6f}, observed={observed[name]:.6f}")

        print(f"{prefix}effective air padding: {report['effective_air_padding_m']:.6e} m")

        warnings = report["warnings"]
        if warnings:
            print(f"{prefix}warnings:")
            for warning in warnings:
                print(f"{prefix}  - {warning}")
        else:
            print(f"{prefix}warnings: none")

    @staticmethod
    def coil_phantom_region_resolution_policy(
        *,
        resolution: float,
        coil_resolution: Optional[float] = None,
        phantom_resolution: Optional[float] = None,
        air_resolution: Optional[float] = None,
        region_resolutions: Optional[Dict[str, float]] = None,
    ) -> Dict[str, float]:
        """Resolve per-region mesh resolution policy for coil+phantom fixtures.

        Resolution precedence (highest to lowest):
        1) Explicit keyword args (``coil_resolution``, ``phantom_resolution``, ``air_resolution``)
        2) ``region_resolutions`` mapping keys: ``coil``, ``phantom``, ``air``
        3) Global ``resolution`` fallback
        """
        if resolution <= 0.0:
            raise ValueError("resolution must be > 0")

        mapping = region_resolutions or {}
        allowed_mapping_keys = {"coil", "phantom", "air"}
        unknown_keys = sorted(set(mapping.keys()) - allowed_mapping_keys)
        if unknown_keys:
            raise ValueError(
                "region_resolutions contains unsupported keys: "
                f"{unknown_keys}; allowed keys are {sorted(allowed_mapping_keys)}"
            )

        coil_h = float(
            coil_resolution if coil_resolution is not None else mapping.get("coil", resolution)
        )
        phantom_h = float(
            phantom_resolution
            if phantom_resolution is not None
            else mapping.get("phantom", resolution)
        )
        air_h = float(
            air_resolution if air_resolution is not None else mapping.get("air", resolution)
        )

        for name, h in (("coil", coil_h), ("phantom", phantom_h), ("air", air_h)):
            if h <= 0.0:
                raise ValueError(f"{name}_resolution must be > 0 (got {h:.6e})")

        return {
            "coil_resolution_m": coil_h,
            "phantom_resolution_m": phantom_h,
            "air_resolution_m": air_h,
            "min_resolution_m": min(coil_h, phantom_h, air_h),
            "max_resolution_m": max(coil_h, phantom_h, air_h),
        }

    @staticmethod
    def coil_phantom_domain(
        coil_major_radius: float = 0.08,
        coil_minor_radius: float = 0.01,
        coil_separation: float = 0.08,
        phantom_radius: float = 0.04,
        phantom_height: float = 0.10,
        air_padding: float = 0.04,
        resolution: float = 0.015,
        coil_resolution: Optional[float] = None,
        phantom_resolution: Optional[float] = None,
        air_resolution: Optional[float] = None,
        region_resolutions: Optional[Dict[str, float]] = None,
        phantom_placement_preset: str = "centered",
        phantom_offset_xy: Optional[Tuple[float, float]] = None,
        comm: MPI.Intracomm = MPI.COMM_WORLD,
        rank: int = 0,
    ) -> Tuple[dolfinx.mesh.Mesh, dolfinx.mesh.MeshTags, dolfinx.mesh.MeshTags]:
        """Generate a coarse two-coil + cylindrical phantom + air mesh.

        Parameters
        ----------
        coil_resolution, phantom_resolution, air_resolution : float | None
            Optional per-region characteristic mesh sizes [m].
            When omitted, ``resolution`` is used. Smaller values refine locally.
        region_resolutions : dict[str, float] | None
            Optional mapping with keys ``coil``, ``phantom``, ``air`` as an
            alternative to individual keyword arguments.
        phantom_placement_preset : str
            Placement preset for phantom center in the xy-plane.
            Supported values: ``centered`` and ``off_center``.
        phantom_offset_xy : tuple[float, float] | None
            Optional explicit ``(x, y)`` phantom center offset [m]. When provided,
            this overrides the preset offset.

        Cell tags:
        - 1: coil_1
        - 2: coil_2
        - 3: phantom
        - 4: air
        """
        allowed_presets = {"centered", "off_center"}
        if phantom_placement_preset not in allowed_presets:
            raise ValueError(
                "phantom_placement_preset must be one of "
                f"{sorted(allowed_presets)}, got '{phantom_placement_preset}'"
            )

        if phantom_offset_xy is None:
            if phantom_placement_preset == "centered":
                phantom_cx, phantom_cy = 0.0, 0.0
            else:
                phantom_cx, phantom_cy = 0.35 * phantom_radius, 0.0
        else:
            phantom_cx, phantom_cy = float(phantom_offset_xy[0]), float(phantom_offset_xy[1])

        phantom_offset_radius = np.sqrt(phantom_cx**2 + phantom_cy**2)
        coil_inner_radius = coil_major_radius - coil_minor_radius
        radial_clearance = coil_inner_radius - (phantom_offset_radius + phantom_radius)
        if radial_clearance <= 0.0:
            raise ValueError(
                "Phantom overlaps coil conductor envelope for selected placement: "
                f"offset_radius={phantom_offset_radius:.6e} m, "
                f"phantom_radius={phantom_radius:.6e} m, "
                f"available_inner_radius={coil_inner_radius:.6e} m"
            )

        sizing = MeshGenerator.coil_phantom_domain_sizing_diagnostics(
            coil_major_radius=coil_major_radius,
            coil_minor_radius=coil_minor_radius,
            coil_separation=coil_separation,
            phantom_radius=phantom_radius,
            phantom_height=phantom_height,
            air_padding=air_padding,
            phantom_offset_xy=(phantom_cx, phantom_cy),
        )
        effective_air_padding = sizing["effective_air_padding_m"]

        if comm.rank == rank and sizing["is_domain_undersized"]:
            print(
                "[coil-phantom-domain] WARNING: requested air_padding is below recommended minimum; "
                f"provided={sizing['provided_air_padding_m']:.6e} m, "
                f"recommended_min={sizing['recommended_min_air_padding_m']:.6e} m. "
                f"Using effective air_padding={effective_air_padding:.6e} m to reduce boundary artifacts."
            )

        region_h = MeshGenerator.coil_phantom_region_resolution_policy(
            resolution=resolution,
            coil_resolution=coil_resolution,
            phantom_resolution=phantom_resolution,
            air_resolution=air_resolution,
            region_resolutions=region_resolutions,
        )

        build_error: Optional[BaseException] = None
        if comm.rank == rank:
            try:
                gmsh.initialize()
                gmsh.model.add("coil_phantom_domain")

                z_offset = coil_separation / 2
                coil_1 = gmsh.model.occ.addTorus(
                    0, 0, -z_offset, coil_major_radius, coil_minor_radius
                )
                coil_2 = gmsh.model.occ.addTorus(
                    0, 0, z_offset, coil_major_radius, coil_minor_radius
                )
                phantom = gmsh.model.occ.addCylinder(
                    phantom_cx,
                    phantom_cy,
                    -phantom_height / 2,
                    0,
                    0,
                    phantom_height,
                    phantom_radius,
                )

                radial_extent = sizing["radial_extent_without_padding_m"]
                z_extent = sizing["z_extent_without_padding_m"]
                air = gmsh.model.occ.addBox(
                    -(radial_extent + effective_air_padding),
                    -(radial_extent + effective_air_padding),
                    -(z_extent + effective_air_padding),
                    2 * (radial_extent + effective_air_padding),
                    2 * (radial_extent + effective_air_padding),
                    2 * (z_extent + effective_air_padding),
                )

                gmsh.model.occ.fragment(
                    [(3, air)],
                    [(3, coil_1), (3, coil_2), (3, phantom)],
                )
                gmsh.model.occ.synchronize()

                volumes = gmsh.model.getEntities(dim=3)
                masses = {tag: gmsh.model.occ.getMass(3, tag) for _, tag in volumes}

                # GEO-9 step 1: the group re-derivation below assumes fragment
                # returned exactly four volumes (air + two coils + phantom). An
                # extra piece (an overlap split off) would silently receive no
                # physical group and surface far downstream as a dolfinx gmshio
                # `assert len(entity_types) == 1`; a merged pair would raise
                # IndexError on the phantom_tag line. Fail here instead, with the
                # count and the masses that identify which happened.
                if len(volumes) != 4:
                    raise RuntimeError(
                        "coil_phantom_domain: occ.fragment returned "
                        f"{len(volumes)} volumes, expected exactly 4 "
                        "(air + coil_1 + coil_2 + phantom); per-volume masses [m^3]: "
                        + ", ".join(
                            f"tag {tag}: {mass:.6e}" for tag, mass in sorted(masses.items())
                        )
                    )

                air_tag = max(masses, key=masses.get)

                remaining = [tag for _, tag in volumes if tag != air_tag]
                z_centers = {tag: gmsh.model.occ.getCenterOfMass(3, tag)[2] for tag in remaining}

                coil_1_tag = min(remaining, key=lambda tag: z_centers[tag])
                coil_2_tag = max(remaining, key=lambda tag: z_centers[tag])
                phantom_tag = [tag for tag in remaining if tag not in (coil_1_tag, coil_2_tag)][0]

                gmsh.model.addPhysicalGroup(3, [coil_1_tag], tag=1)
                gmsh.model.setPhysicalName(3, 1, "coil_1")
                gmsh.model.addPhysicalGroup(3, [coil_2_tag], tag=2)
                gmsh.model.setPhysicalName(3, 2, "coil_2")
                gmsh.model.addPhysicalGroup(3, [phantom_tag], tag=3)
                gmsh.model.setPhysicalName(3, 3, "phantom")
                gmsh.model.addPhysicalGroup(3, [air_tag], tag=4)
                gmsh.model.setPhysicalName(3, 4, "air")

                # Every 3-D entity must carry a marker; a cell with none is exactly
                # what gmshio asserts on. Checked against the model, not against the
                # four tags we just wrote, so a renumbering by fragment cannot hide.
                grouped_volumes = set()
                for _, group_tag in gmsh.model.getPhysicalGroups(dim=3):
                    grouped_volumes.update(gmsh.model.getEntitiesForPhysicalGroup(3, group_tag))
                ungrouped = sorted({tag for _, tag in volumes} - grouped_volumes)
                if ungrouped:
                    raise RuntimeError(
                        "coil_phantom_domain: 3-D entities carry no physical group: "
                        f"{ungrouped}; masses [m^3]: "
                        + ", ".join(f"tag {tag}: {masses[tag]:.6e}" for tag in ungrouped)
                    )

                print(
                    f"[coil-phantom-mesh] fragment volumes={len(volumes)} masses[m^3]: "
                    + ", ".join(f"{tag}:{mass:.6e}" for tag, mass in sorted(masses.items()))
                    + f" | air={air_tag} coil_1={coil_1_tag} coil_2={coil_2_tag} phantom={phantom_tag}",
                    flush=True,
                )

                outer_boundary_surfaces = []
                box_half_x = radial_extent + effective_air_padding
                box_half_z = z_extent + effective_air_padding
                for dim, surf in gmsh.model.getEntities(dim=2):
                    x_min, y_min, z_min, x_max, y_max, z_max = gmsh.model.getBoundingBox(dim, surf)
                    if (
                        abs(abs(x_min) - box_half_x) < resolution
                        or abs(abs(x_max) - box_half_x) < resolution
                        or abs(abs(y_min) - box_half_x) < resolution
                        or abs(abs(y_max) - box_half_x) < resolution
                        or abs(abs(z_min) - box_half_z) < resolution
                        or abs(abs(z_max) - box_half_z) < resolution
                    ):
                        outer_boundary_surfaces.append(surf)

                if outer_boundary_surfaces:
                    gmsh.model.addPhysicalGroup(2, outer_boundary_surfaces, tag=1)
                    gmsh.model.setPhysicalName(2, 1, "outer_boundary")

                # GEO-17 step 1 (2026-08-20): the region sizes are carried by gmsh
                # size *fields*, composed with `Min`, not by `mesh.setSize` on CAD
                # points. The point path this replaces was inert: it walked
                # volume -> surfaces -> curves -> points with
                # `gmsh.model.getBoundary`'s default `combined=True`, and the
                # boundary of the *combined* closed shell of a volume is empty, so
                # every region collected zero points and `setSize` was never called
                # at all (measured `20260820T110127Z_GEO-17-step1-diag.log`: "air:
                # 0 pts -> NO SIZE SET" for all four regions, both sizings). Mesh
                # size then came only from the CharacteristicLength clamps, so
                # asking for a finer coil (0.012) while any region asked for a
                # coarser one (air 0.020) meshed the coil at the 0.020 ceiling and
                # LOST 22% of its volume — known-issues "Four defects" §1. With a
                # `Min` over per-volume Constant fields, a region's request bounds
                # the size on its own boundary, so a shared curved interface takes
                # the finer of the two neighbours rather than whatever the clamps
                # allow.
                region_size_requests = (
                    ("air", air_tag, region_h["air_resolution_m"]),
                    ("coil_1", coil_1_tag, region_h["coil_resolution_m"]),
                    ("coil_2", coil_2_tag, region_h["coil_resolution_m"]),
                    ("phantom", phantom_tag, region_h["phantom_resolution_m"]),
                )

                region_field_ids = []
                for _, volume_tag, size_value in region_size_requests:
                    field_id = gmsh.model.mesh.field.add("Constant")
                    gmsh.model.mesh.field.setNumbers(field_id, "VolumesList", [volume_tag])
                    gmsh.model.mesh.field.setNumber(field_id, "IncludeBoundary", 1)
                    gmsh.model.mesh.field.setNumber(field_id, "VIn", size_value)
                    # Outside its own volume a region must not constrain anything;
                    # `Min` then reduces to the other fields there.
                    gmsh.model.mesh.field.setNumber(field_id, "VOut", 1.0e22)
                    region_field_ids.append(field_id)

                size_field = gmsh.model.mesh.field.add("Min")
                gmsh.model.mesh.field.setNumbers(size_field, "FieldsList", region_field_ids)
                gmsh.model.mesh.field.setAsBackgroundMesh(size_field)

                # gmsh's own sizing heuristics stay at their defaults on purpose.
                # Measured both ways (`20260820T110302Z_GEO-17-step1-fieldfix.log`
                # with MeshSizeExtendFromBoundary/FromPoints/FromCurvature forced to
                # 0, `20260820T110407Z_GEO-17-step1-probe-defaults.log` with them
                # left alone): forcing them off moves the *uniform* mesh's tagged
                # volumes (coil_1 1.191750413e-04 -> 1.154535949e-04 m^3, -3.12%),
                # while at the defaults the uniform table reproduces exactly and the
                # policy mesh gains more coil volume (+10.72% vs +13.21% off a
                # 3% lower base: 1.319468693e-04 vs 1.307098011e-04 m^3). With no
                # point sizes set anywhere, the boundary-extension heuristic has
                # nothing of its own to extend and does not compete with the field.

                gmsh.option.setNumber("Mesh.CharacteristicLengthMin", region_h["min_resolution_m"])
                gmsh.option.setNumber("Mesh.CharacteristicLengthMax", region_h["max_resolution_m"])

                print(
                    "[coil-phantom-mesh] sizing ownership: "
                    + ", ".join(
                        f"{name}: field {field_id} on volume {volume_tag} -> {size_value:.6e} m"
                        for (name, volume_tag, size_value), field_id in zip(
                            region_size_requests, region_field_ids
                        )
                    )
                    + f" | Min field {size_field}"
                    + f" | clamps min={region_h['min_resolution_m']:.6e} "
                    f"max={region_h['max_resolution_m']:.6e} m",
                    flush=True,
                )

                print(
                    "[coil-phantom-mesh] region resolution policy: "
                    f"coil={region_h['coil_resolution_m']:.6e} m, "
                    f"phantom={region_h['phantom_resolution_m']:.6e} m, "
                    f"air={region_h['air_resolution_m']:.6e} m"
                )

                gmsh.model.mesh.generate(3)
                gmsh.model.mesh.optimize("Netgen")
            except BaseException as exc:  # noqa: BLE001 - re-raised on every rank
                build_error = exc
                if gmsh.isInitialized():
                    gmsh.finalize()

        # Collective: every rank learns whether the geometry exists before
        # any of them enters the collective `_model_to_mesh` below (GEO-23
        # step 2a; the unwrapped throw deadlocked at 120 s instead of
        # failing in seconds).
        _raise_geometry_failure_on_every_rank(
            comm, rank, "coil_phantom_domain", build_error, resolution
        )

        mesh, cell_tags, facet_tags = _model_to_mesh(gmsh.model, comm, rank, gdim=3)

        sanity_report = MeshGenerator.coil_phantom_geometry_sanity_report(
            cell_tags=cell_tags,
            coil_major_radius=coil_major_radius,
            coil_minor_radius=coil_minor_radius,
            coil_separation=coil_separation,
            phantom_radius=phantom_radius,
            phantom_height=phantom_height,
            air_padding=air_padding,
            phantom_offset_xy=(phantom_cx, phantom_cy),
            comm=comm,
        )
        MeshGenerator.print_coil_phantom_geometry_sanity_report(report=sanity_report, comm=comm)

        if comm.rank == rank:
            gmsh.finalize()

        return mesh, cell_tags, facet_tags

    @staticmethod
    def birdcage_port_layout_diagnostics(
        *,
        leg_count: int,
        ring_radius: float,
        leg_width: float,
        ring_minor_radius: float,
        phantom_radius: float,
        port_box_size: Tuple[float, float, float],
        port_clearance: float = 1.0e-3,
        min_port_face_area: float = 2.5e-5,
        min_port_center_separation: Optional[float] = None,
    ) -> Dict[str, float]:
        """Compute and validate birdcage port geometry diagnostics.

        This helper is intentionally lightweight and does not require meshing.
        It enforces a minimum port face area, checks center-to-center
        separation, and guards against radial overlap with conductor/phantom
        bulk regions.
        """
        if leg_count < 3:
            raise ValueError("leg_count must be >= 3 for port layout checks")
        if leg_width <= 0.0:
            raise ValueError("leg_width must be > 0")
        if ring_radius <= 0.0:
            raise ValueError("ring_radius must be > 0")

        port_dx, port_dy, port_dz = port_box_size
        if port_dx <= 0.0 or port_dy <= 0.0 or port_dz <= 0.0:
            raise ValueError("port_box_size components must be > 0")
        if ring_minor_radius <= 0.0:
            raise ValueError("ring_minor_radius must be > 0")
        if phantom_radius <= 0.0:
            raise ValueError("phantom_radius must be > 0")
        if port_clearance < 0.0:
            raise ValueError("port_clearance must be >= 0")

        port_face_area = port_dx * port_dz
        if port_face_area < min_port_face_area:
            raise ValueError(
                "Port face area too small for robust tagging: "
                f"area={port_face_area:.6e} m^2, required>={min_port_face_area:.6e} m^2"
            )

        leg_radius_eff = 0.5 * leg_width
        conductor_outer_radius = ring_radius + max(leg_radius_eff, ring_minor_radius)
        port_radius = conductor_outer_radius + 0.5 * port_dy + port_clearance

        theta = np.linspace(0.0, 2.0 * np.pi, leg_count, endpoint=False)
        port_centers = []
        for idx, angle in enumerate(theta):
            next_angle = theta[(idx + 1) % leg_count]
            midpoint_angle = np.arctan2(
                np.sin(angle) + np.sin(next_angle),
                np.cos(angle) + np.cos(next_angle),
            )
            port_centers.append(
                np.array(
                    [
                        port_radius * np.cos(midpoint_angle),
                        port_radius * np.sin(midpoint_angle),
                    ]
                )
            )

        min_center_separation = float("inf")
        for idx in range(len(port_centers)):
            for jdx in range(idx + 1, len(port_centers)):
                distance = float(np.linalg.norm(port_centers[idx] - port_centers[jdx]))
                min_center_separation = min(min_center_separation, distance)

        if min_port_center_separation is None:
            min_port_center_separation = max(5.0e-4, 1.25 * max(port_dx, port_dy))

        if min_center_separation < min_port_center_separation:
            raise ValueError(
                "Port center separation too small: "
                f"min={min_center_separation:.6e} m, "
                f"required>={min_port_center_separation:.6e} m"
            )

        conductor_radial_clearance = port_radius - 0.5 * port_dy - conductor_outer_radius
        phantom_radial_clearance = port_radius - 0.5 * port_dy - phantom_radius

        if conductor_radial_clearance <= 0.0:
            raise ValueError(
                "Port/conductor radial overlap detected: "
                f"clearance={conductor_radial_clearance:.6e} m"
            )
        if phantom_radial_clearance <= 0.0:
            raise ValueError(
                "Port/phantom radial overlap detected: "
                f"clearance={phantom_radial_clearance:.6e} m"
            )

        return {
            "port_face_area_m2": float(port_face_area),
            "min_port_face_area_m2": float(min_port_face_area),
            "min_port_center_separation_m": float(min_center_separation),
            "required_port_center_separation_m": float(min_port_center_separation),
            "conductor_radial_clearance_m": float(conductor_radial_clearance),
            "phantom_radial_clearance_m": float(phantom_radial_clearance),
            "port_radius_m": float(port_radius),
            "conductor_outer_radius_m": float(conductor_outer_radius),
        }

    @staticmethod
    def _birdcage_leg_gap_layout(
        *,
        leg_count: int,
        ring_radius: float,
        leg_radius_eff: float,
        leg_spacing: float,
        coil_length: float,
        ring_minor_radius: float,
        phantom_radius: float,
        leg_gap_length: float,
        port_clearance: float,
        min_port_face_area: float,
        min_port_center_separation: Optional[float],
        leg_azimuth_offsets_rad: Optional[Sequence[float]] = None,
    ) -> Tuple[Tuple[float, float, float], Dict[str, float]]:
        """Size and validate the `GEO-18` gapped-leg port layout.

        Returns the derived ``port_box_size`` and the diagnostics dict
        `birdcage_port_domain` consumes. The box is square in transverse
        section and exactly as tall as the gap, so its z-faces are the leg
        stubs' cut faces.
        """
        if leg_gap_length <= 0.0:
            raise ValueError("leg_gap_length must be > 0")
        if port_clearance <= 0.0:
            raise ValueError("leg_gap_length requires port_clearance > 0")
        stub_length = 0.5 * (coil_length - leg_gap_length)
        if stub_length <= 0.0:
            raise ValueError(
                "leg_gap_length leaves no leg: "
                f"gap={leg_gap_length:.6e} m, coil_length={coil_length:.6e} m"
            )
        # The cut must stay clear of the rings, or the removed segment is not a
        # cylinder and the closed-form terminal/mass forms below are wrong.
        ring_inner_z = 0.5 * leg_spacing - ring_minor_radius
        ring_clearance = ring_inner_z - 0.5 * leg_gap_length
        if ring_clearance <= 0.0:
            raise ValueError("Leg gap reaches the rings: " f"clearance={ring_clearance:.6e} m")

        box_width = 2.0 * leg_radius_eff + 2.0 * port_clearance
        port_box_size = (box_width, box_width, leg_gap_length)

        port_face_area = box_width * leg_gap_length
        if port_face_area < min_port_face_area:
            raise ValueError(
                "Port face area too small for robust tagging: "
                f"area={port_face_area:.6e} m^2, required>={min_port_face_area:.6e} m^2"
            )

        # Ports now sit on the legs, so their separation is the leg separation.
        # `PORT-9` step 3 leg (d1): with a per-leg azimuth offset the spacing is
        # no longer uniform, so the pairwise minimum below — which already
        # scanned every pair rather than assuming one gap — is the check that has
        # to see the displaced layout, and that is why the offsets reach here at
        # all rather than only the builder.
        theta = np.linspace(0.0, 2.0 * np.pi, leg_count, endpoint=False)
        if leg_azimuth_offsets_rad is not None:
            theta = theta + np.asarray(leg_azimuth_offsets_rad, dtype=float)
        centers = [np.array([ring_radius * np.cos(a), ring_radius * np.sin(a)]) for a in theta]
        min_center_separation = float("inf")
        for idx in range(len(centers)):
            for jdx in range(idx + 1, len(centers)):
                min_center_separation = min(
                    min_center_separation,
                    float(np.linalg.norm(centers[idx] - centers[jdx])),
                )
        if min_port_center_separation is None:
            min_port_center_separation = max(5.0e-4, 1.25 * box_width)
        if min_center_separation < min_port_center_separation:
            raise ValueError(
                "Port center separation too small: "
                f"min={min_center_separation:.6e} m, "
                f"required>={min_port_center_separation:.6e} m"
            )

        # The box is axis-aligned, so its nearest corner to the axis is at
        # radius >= ring_radius - (sqrt(2)/2)*box_width; keep the phantom out of
        # it, or the port region is not the clean air+two-disks volume assumed.
        phantom_radial_clearance = ring_radius - 0.5 * np.sqrt(2.0) * box_width - phantom_radius
        if phantom_radial_clearance <= 0.0:
            raise ValueError(
                "Port/phantom radial overlap detected: "
                f"clearance={phantom_radial_clearance:.6e} m"
            )

        return port_box_size, {
            "port_face_area_m2": float(port_face_area),
            "min_port_face_area_m2": float(min_port_face_area),
            "min_port_center_separation_m": float(min_center_separation),
            "required_port_center_separation_m": float(min_port_center_separation),
            "conductor_radial_clearance_m": 0.0,
            "phantom_radial_clearance_m": float(phantom_radial_clearance),
            "port_radius_m": float(ring_radius),
            "conductor_outer_radius_m": float(ring_radius + max(leg_radius_eff, ring_minor_radius)),
            "leg_gap_length_m": float(leg_gap_length),
            "leg_stub_length_m": float(stub_length),
            "leg_gap_ring_clearance_m": float(ring_clearance),
        }

    @staticmethod
    def _birdcage_ring_gap_layout(
        *,
        leg_count: int,
        ring_radius: float,
        leg_radius_eff: float,
        leg_spacing: float,
        ring_minor_radius: float,
        phantom_radius: float,
        ring_gap_length: float,
        port_clearance: float,
        min_port_face_area: float,
    ) -> Dict[str, float]:
        """Size and validate the `GEO-20` ring-gap port layout (high-pass).

        The gap is an *arc* of length ``g`` centred on the mid-azimuth between
        each adjacent leg pair, on both end rings — ``2·leg_count`` of them. It
        is cut by the two half-planes ``phi = phi_c +- alpha`` with
        ``alpha = g / (2·ring_radius)``, so both cut faces are exact planar
        disks of area ``pi·ring_minor_radius²``: the closed form only exists
        because the cut planes are *radial*, which is why the box that spans
        the gap must live in the gap's own ``(r_hat, phi_hat, z_hat)`` frame
        rather than being axis-aligned (`GEO-18`'s leg construction rotated).

        That port solid is the box ``|z − z_ring| <= w/2``,
        ``|u − ring_radius| <= w/2`` (``u`` = the radial coordinate *of the
        gap's own frame*, i.e. ``rho·cos(phi − phi_c)``) intersected with the
        wedge ``|phi − phi_c| <= alpha``. Every one of its six faces is planar,
        so its volume, surface area and mid-plane section all have exact closed
        forms and a linear mesh reproduces them to machine precision — a
        constant-``rho`` face would not, being curved.

        Returns the closed forms and clearances; the geometry itself is built
        in `_build_birdcage_port_model`.
        """
        if ring_gap_length <= 0.0:
            raise ValueError("ring_gap_length must be > 0")
        if port_clearance <= 0.0:
            raise ValueError("ring_gap_length requires port_clearance > 0")

        half_angle = 0.5 * ring_gap_length / ring_radius
        leg_half_angle = np.pi / leg_count
        if half_angle >= leg_half_angle:
            raise ValueError(
                "ring_gap_length spans past the neighbouring leg: "
                f"half-angle={half_angle:.6e} rad, leg half-spacing="
                f"{leg_half_angle:.6e} rad"
            )

        box_width = 2.0 * ring_minor_radius + 2.0 * port_clearance
        # Arc clearance between the gap's end plane and the nearest leg's
        # surface, measured along the ring centreline. The port solid spans the
        # gap exactly, so this is also the metal the terminals sit on.
        leg_arc_clearance = ring_radius * (leg_half_angle - half_angle) - leg_radius_eff
        if leg_arc_clearance <= 0.0:
            raise ValueError("Ring gap reaches the legs: " f"clearance={leg_arc_clearance:.6e} m")

        port_face_area = np.pi * ring_minor_radius**2
        if port_face_area < min_port_face_area:
            raise ValueError(
                "Ring-gap terminal area too small for robust tagging: "
                f"area={port_face_area:.6e} m^2, "
                f"required>={min_port_face_area:.6e} m^2"
            )

        # The port solid's innermost point sits at u = ring_radius - w/2 on the
        # gap's own radial ray; the phantom must stay outside it or the port
        # region is not the clean air+two-disks volume assumed.
        phantom_radial_clearance = ring_radius - 0.5 * box_width - phantom_radius
        if phantom_radial_clearance <= 0.0:
            raise ValueError(
                "Ring port/phantom radial overlap detected: "
                f"clearance={phantom_radial_clearance:.6e} m"
            )
        # The gap solid straddles the ring plane; it must not reach the legs'
        # own gap region or the ring tube's own z extent is not what is cut.
        if 0.5 * box_width <= ring_minor_radius:
            raise ValueError(
                "Ring port box is thinner than the ring tube: "
                f"w/2={0.5 * box_width:.6e} m, r={ring_minor_radius:.6e} m"
            )

        tan_half = float(np.tan(half_angle))
        return {
            "ring_gap_length_m": float(ring_gap_length),
            "ring_gap_half_angle_rad": float(half_angle),
            "ring_port_box_width_m": float(box_width),
            "ring_port_count": int(2 * leg_count),
            # V = 2·R·w²·tan(alpha) — the wedge∩slab∩slab solid, exact.
            "ring_port_volume_m3": float(2.0 * ring_radius * box_width**2 * tan_half),
            # A = 2·w²/cos(alpha) + 8·R·w·tan(alpha) — two radial caps, two
            # z-faces, two u-faces, all planar.
            "ring_port_surface_m2": float(
                2.0 * box_width**2 / float(np.cos(half_angle))
                + 8.0 * ring_radius * box_width * tan_half
            ),
            # The mid-plane section (phi = phi_c) is the w x w rectangle.
            "ring_port_sheet_area_m2": float(box_width**2),
            # `GEO-26`: the gap's *chord* at u = ring_radius, 2·R·tan(alpha).
            # `ring_gap_length_m` above is the *arc* g = 2·R·alpha; the port
            # box's two radial cap faces are planar, so what a longitudinal
            # sheet in the plane u = R actually spans is the chord, not the
            # arc — +0.10% at the fixture's g = 0.008 m / R = 0.07 m. Both are
            # emitted so a caller never has to pick by re-deriving.
            "ring_port_gap_chord_m": float(2.0 * ring_radius * tan_half),
            # The longitudinal sheet (`ring_sheet_orientation="longitudinal"`):
            # the chord x w rectangle in the plane u = R, normal u_hat(phi_c).
            # This is the surface the lumped-sheet port model needs — its
            # extent along the drive direction phi_hat is the chord, i.e. the
            # `h` that the transverse mid-section above offers as exactly 0
            # (`PORT-13` step 1's blocker).
            "ring_port_sheet_longitudinal_area_m2": float(2.0 * ring_radius * tan_half * box_width),
            "ring_terminal_area_m2": float(2.0 * port_face_area),
            "ring_removed_mass_m3": float(
                2 * leg_count * np.pi * ring_minor_radius**2 * ring_gap_length
            ),
            "ring_leg_arc_clearance_m": float(leg_arc_clearance),
            "ring_phantom_radial_clearance_m": float(phantom_radial_clearance),
            "ring_gap_z_m": float(0.5 * leg_spacing),
        }

    @staticmethod
    def birdcage_port_domain(
        leg_count: int = 4,
        ring_radius: float = 0.07,
        leg_width: float = 0.012,
        leg_spacing: float = 0.11,
        coil_length: float = 0.14,
        ring_minor_radius: float = 0.004,
        phantom_radius: float = 0.03,
        phantom_height: float = 0.08,
        port_box_size: Tuple[float, float, float] = (0.010, 0.008, 0.010),
        port_clearance: float = 1.0e-3,
        min_port_face_area: float = 2.5e-5,
        min_port_center_separation: Optional[float] = None,
        leg_gap_length: Optional[float] = None,
        ring_gap_length: Optional[float] = None,
        leg_azimuth_offsets_rad: Optional[Sequence[float]] = None,
        emit_port_sheets: bool = False,
        air_padding: float = 0.03,
        resolution: float = 0.015,
        conductor_resolution: Optional[float] = None,
        conductor_refine_distance: Optional[float] = None,
        phantom_resolution: Optional[float] = None,
        comm: MPI.Intracomm = MPI.COMM_WORLD,
        rank: int = 0,
        n_legs: Optional[int] = None,
        leg_radius: Optional[float] = None,
        leg_height: Optional[float] = None,
        return_diagnostics: bool = False,
        *,
        ring_sheet_orientation: str = "transverse",
        as_hole: bool = False,
        c4_congruent_sheets: bool = False,
        step_export_path: Optional[str] = None,
        symmetry_quadrant: bool = False,
        leg_sheet_strip_fraction: Optional[float] = None,
        implant_rod: Optional[Dict[str, object]] = None,
        emit_leg_ports: bool = True,
        implant_spheroid: Optional[Dict[str, object]] = None,
        shield_radius: Optional[float] = None,
        phantom_box: Optional[Dict[str, object]] = None,
    ):
        """Generate a coarse, parametric birdcage-like geometry fixture with port tags.

        ``phantom_box`` (`GEO-38` step 2) is additive and keyword-only; ``None``
        (default) executes no new line. ``{"size_m": (Lx, Ly, Lz),
        "table_offset_m": off}`` replaces the phantom cylinder by the
        axis-aligned box ``[-Lx/2, Lx/2] x [-Ly/2 + off, Ly/2 + off] x
        [-Lz/2, Lz/2]`` (cell tag 3, kept as cells on the hole route);
        ``z_extent`` uses ``Lz/2`` and the layout clearances use the box's xy
        circumradius ``hypot(Lx/2, Ly/2 + |off|)`` as ``phantom_radius``.
        Refused with ``symmetry_quadrant`` and
        ``phantom_resolution``. With it the build
        also reports ``ring_sheet_gmsh_triangles`` ({ordinal: triangles on the
        post-fragment ring sheet}) and, on ``as_hole`` with ``shield_radius``,
        ``conductor_cad_mass_m3`` (the can's CAD mass minus the cut air, minus
        the implant body's CAD mass when one was cut from that air).
        ``implant_spheroid`` composes with ``phantom_box`` (`IMP-6` step 4),
        its containment then tested against the box's half-extents; so does
        ``implant_rod`` (`GEO-36` step 8), for rod, coat and tip ball.

        ``shield_radius`` (`GEO-38` step 1a) is additive and keyword-only;
        ``None`` (default) executes no new line. A float replaces the PEC air
        box by a closed PEC can: the cylinder of that radius spanning
        ``±z_extent`` (lateral face plus the two end caps, all outer
        boundary). Refused with ``symmetry_quadrant`` and ``leg_gap_length``;
        must clear the coil, the ring-gap box corners and the floating leg
        box corners by ``port_clearance`` (diagnostic
        ``leg_box_corner_radius_m``).

        ``emit_leg_ports`` (`GEO-35` step 1) is additive and keyword-only;
        ``True`` (default) executes no new line. ``False`` on an unbroken-leg
        layout (``leg_gap_length=None``) does not build the ``leg_count``
        floating port boxes at the mid-leg azimuths; the ring-port ordinals
        stay ``leg_count+1 .. 3·leg_count`` (ordinals ``1 .. leg_count`` are
        reserved gaps), so no ring tag is renumbered. Refused with
        ``leg_gap_length`` (the leg boxes are the drive there), ``as_hole`` or
        ``symmetry_quadrant`` (later steps).

        ``implant_rod`` (`GEO-36` step 2) is additive and keyword-only; ``None``
        (default) executes no new line.  A mapping with ``radius`` and
        ``length`` [m] (``length`` is the cylindrical section between the two
        hemispherical cap centres, so tip to tip is ``length + 2*radius``),
        optional ``centre`` (default the origin) and ``axis`` (default ``ẑ``)
        cuts that capped rod out of the phantom **and** the air box before the
        fragment, so it is a hole; its wall (read off the meshed volumes'
        boundary) carries the exterior facet group
        :data:`BIRDCAGE_IMPLANT_SURFACE_TAG`.  The wall is graded (Distance ->
        Threshold ``r/3`` -> ``resolution`` over ``r`` -> ``10 r``, plus a
        ``Ball`` of radius ``2 r`` and size ``r/4`` at each cap centre).
        Refused with ``symmetry_quadrant``, with ``as_hole`` unless kept and
        uncoated (`GEO-36` step 6), for ``radius`` or
        ``length`` <= 0, and for a rod not strictly inside the phantom.
        `GEO-36` step 3: ``"keep_cells": True`` (a bool; anything else is
        refused before gmsh) keeps the rod's cells instead — cell group
        :data:`BIRDCAGE_IMPLANT_CELL_TAG`, and 402 becomes the *interior*
        phantom/implant interface, rebuilt from the cell tags after
        ``model_to_mesh``.  Absent or ``False`` is the hole above, unchanged.

        ``leg_sheet_strip_fraction`` (`OPS-61` step 1b) is additive and
        keyword-only: when set to ``f`` in (0, 1), each leg-port sheet is
        fragmented at azimuth 0 against its centre rectangle of width
        ``f * port_dx`` (full height) before placement, so the strip edges are
        CAD edges and a midpoint filter selects exactly that rectangle. Needs
        ``emit_port_sheets=True`` and ``leg_gap_length``; refuses
        ``c4_congruent_sheets`` and ``symmetry_quadrant`` (composes with
        ``as_hole`` since `ANS-6` step 4).
        ``None`` (default) executes no new line.

        ``symmetry_quadrant`` (`TH-16` step 3a) is additive, keyword-only and
        defaults to ``False``, which executes no new gmsh call.  ``True`` meshes
        only the x >= 0, y >= 0 quadrant: the air box becomes
        [0, r_e] x [0, r_e] x [-z_e, z_e], the leg-port sheet tools are dropped
        (the kept port-box halves' faces on y = 0 / x = 0 *are* the sheets), and
        every fragment piece whose centroid lies at x < 0 or y < 0 is removed.
        The cut faces carry :data:`BIRDCAGE_CUT_X0_TAG` / :data:`BIRDCAGE_CUT_Y0_TAG`
        and are excluded from the outer-boundary group.  Needs
        ``emit_port_sheets`` and one of two layouts: L4 (``leg_gap_length``,
        ``leg_count == 4``, no ring gaps, not ``as_hole``) or LH (`TH-20` step
        0a: ``ring_gap_length``, no leg gaps, ``leg_count % 4 == 0``,
        ``as_hole=True``, ``ring_sheet_orientation="longitudinal"``; the ring
        sheets stay fragment tools and every kept ring port lies wholly inside
        the quadrant). Refuses ``c4_congruent_sheets``, azimuth offsets, both
        gap kinds together, the solid ring-gap quadrant and transverse sheets.

        ``step_export_path`` (`GEO-34` step 1a) is additive, keyword-only and
        defaults to ``None``, which executes no new gmsh call: every existing
        caller gets the identical mesh.  A path makes the building rank write
        the finished OCC model (after the last ``occ.synchronize()``, sizing
        fields set, before ``generate(3)``) with gmsh's own STEP writer — the
        emitter `io/step_import.import_step` re-imports.

        ``c4_congruent_sheets`` (`GEO-32`) is additive, keyword-only and defaults
        to ``False``, which executes no new gmsh call: every existing caller gets
        the identical mesh.  With ``True`` the four leg-port sheets are meshed as
        exact C4 copies of sheet 1 — after the fragment, sheets 2..4 are bound to
        sheet 1 with ``gmsh.model.mesh.setPeriodic`` under the snapped
        ``_z_rotation_affine(theta_k - theta_1)``, so gmsh copies sheet 1's
        triangulation instead of cutting each sheet on its own (`GEO-31`
        measured the four sheets as the same patch cut four ways, facet-centroid
        Hausdorff 0.44-0.99 of a facet edge).  Valid only with
        ``emit_port_sheets=True``, `leg_gap_length` set and ``leg_count == 4``;
        anything else raises ``ValueError``.

        ``as_hole`` (`TH-15` step 3a) is additive and defaults to ``False``, so
        every existing caller builds the identical mesh it always has.  With
        ``True`` the ring and leg solids are **cut out** of the air box
        (``removeTool=False``, then dropped from the model) instead of being
        meshed: the coil becomes a cavity whose wall carries the dim-2 physical
        group :data:`BIRDCAGE_CONDUCTOR_SURFACE_TAG`, which is where a PEC
        (``n x E = 0``) condition goes.  The phantom, the port boxes and their
        sheets stay exactly what they were.  It needs ``emit_port_sheets=True``:
        without terminals the hole has no drive at all.  Same trap as
        `two_torus_domain(as_hole=True)` — the retained tools' faces are *not*
        the cavity's, so the wall is identified from the **meshed** volumes'
        boundary (step 0's abort: 716 ``nodes not attached to any tet``).

        Parameters
        ----------
        leg_count : int
            Number of birdcage legs distributed uniformly around the ring.
        ring_radius : float
            Radius of the birdcage rings [m].
        leg_width : float
            Leg diameter [m].
        leg_spacing : float
            Center-to-center spacing between bottom and top ring planes [m].
        coil_length : float
            Axial conductor span used for vertical legs [m].
        leg_gap_length : float, optional
            `GEO-18` step 1. ``None`` (default) keeps today's geometry
            bit-for-bit: unbroken legs, and port boxes floating in the air at
            the *midpoint* azimuth between adjacent legs, touching no metal
            (measured: conductor-facing area exactly 0 m² per port, `PORT-9`
            step 3 leg (b)). When given, the segment ``|z| <= g/2`` is removed
            from every leg **before** the fragment, and each port box is
            re-placed centred on its own leg axis spanning exactly the gap
            (``dz = g``, square transverse ``dx = dy = 2·r_leg +
            2·port_clearance``), so the two leg stubs' planar cut faces lie in
            the box's z-faces: the conductor↔port interface is two disks of
            closed-form area ``π·r_leg²`` each, the drive direction is ``ẑ``
            for every port, and the four-port layout is exactly C4-invariant.
            `port_box_size` is ignored in this mode — the gap sets it.
        ring_gap_length : float, optional
            `GEO-20` step 1 — the **high-pass** layout, independent of
            `leg_gap_length` (a high-pass birdcage has uncut legs, so ``None``
            for that one is the normal pairing). ``None`` (default) keeps the
            rings whole and the geometry bit-for-bit. When given, the arc
            ``|phi − phi_mid| <= g/(2·ring_radius)`` is removed from **both**
            end rings at the mid-azimuth between every adjacent leg pair —
            ``2·leg_count`` gaps — and one port solid is placed in each,
            spanning the gap exactly in the gap's own ``(r_hat, phi_hat,
            z_hat)`` frame. The cut planes are radial, so the two terminals per
            port are exact planar disks of area ``pi·ring_minor_radius²`` and
            the drive direction is ``phi_hat``. The ring ports take ordinals
            ``leg_count+1 .. leg_count+2·leg_count``, after the leg ports, so
            every existing tag keeps its meaning.
        ring_sheet_orientation : {'transverse', 'longitudinal'}, optional
            `GEO-26` step 1, keyword-only, meaningful only with
            `ring_gap_length` + `emit_port_sheets`. ``'transverse'`` (default)
            is `GEO-20`'s emission byte-for-byte: the gap's mid-*section* at
            ``phi = phi_c``, the ``w x w`` rectangle spanned by ``û(phi_c)``
            and ``ẑ``, normal ``phi_hat``. Its extent along the drive
            direction is zero (measured ``<= 1.43e-17`` m, `PORT-13` step 1),
            so the lumped-sheet port model — ``R_s = Z_p·w/h``,
            ``I = (1/R_s)∫E·ĥ dS / h`` — has no ``h`` on it and ``E·ĥ`` is the
            normal trace of an H(curl) field on an interior facet.
            ``'longitudinal'`` emits instead the planar rectangle in the plane
            ``u = ring_radius`` (normal ``û(phi_c)``), spanning the gap
            **chord** ``2·R·tan(alpha)`` along ``phi_hat`` and
            `ring_port_box_width_m` along ``ẑ``: a section *through* the gap,
            so ``h`` is the chord and the port halves become the inner
            (``u < R``, tag ``100+i``) and outer (``u > R``, ``200+i``) pieces
            with closed-form volumes ``w·tan(alpha)·(R·w ∓ w²/4)``. Its four
            edges lie on the box's two radial caps and two ``z`` faces, so it
            spans the box exactly. The horizontal trapezoid at ``z = z_ring``
            was considered and rejected: also planar and box-spanning, but its
            ``h(u) = 2u·tan(alpha)`` varies by ``±w/(2R)`` across the sheet and
            the port model needs one ``h``.
        leg_azimuth_offsets_rad : sequence of float, optional
            `PORT-9` step 3 leg (d1), the geometric negative control of the C4
            symmetry gate. ``None`` (the default) and an all-zero sequence both
            keep the uniform layout **bit-for-bit** — the offsets are *added* to
            the nominal azimuths, and adding an exact zero is exact, so the
            all-zero build is the same construction and not a zero-angle
            rotation of it. Otherwise entry ``i`` is added to leg ``i``'s
            azimuth, and leg ``i`` **with everything that belongs to it** — its
            two stubs, its gap, its terminals, its port box and (with
            `emit_port_sheets`) its sheet — moves rigidly about the z axis by
            that angle, while the rings, the phantom and every other leg stay
            put. The coil loses C4; nothing else about it changes. Since
            `GEO-19` step B the box and sheet are already built in the leg's own
            local frame and rotated onto its azimuth, so a displaced port stays
            square about its leg axis and its sheet stays the plane containing
            that axis and its own radial direction. Requires `leg_gap_length`
            (an ungapped layout floats its boxes *between* legs, so a leg has no
            port to carry) and is refused with `ring_gap_length` (the ring gaps
            sit at mid-azimuths that a displaced leg would no longer bisect).
            Length must be `leg_count`.
        emit_port_sheets : bool, optional
            `GEO-18` step 2, requires `leg_gap_length` or `ring_gap_length`.
            A port only gets a sheet if it has terminals: with
            `ring_gap_length` alone, the four floating leg boxes stay single
            pieces at ``100+i`` and only the ring ports split.
            ``False`` (default)
            leaves the gapped mesh bit-for-bit as step 1 built it. When ``True``,
            each port box is split by the axis-aligned coordinate plane
            containing its own leg axis (y-normal for a leg on the x-axis,
            x-normal for one on the y-axis), the rectangle entering the fragment
            as a dim-2 tool so the tets conform to it. The halves become separate
            cell tags — ``100+i`` below the plane, ``200+i`` above — and the sheet
            is what they share, so a caller rebuilds it dolfinx-side as an
            interface facet group over that tag pair. Mesh-side only: no port
            model and no drive live here.
            Mesh-side only: a gapped birdcage carries no lumped elements and
            therefore still cannot resonate.
        conductor_resolution : float, optional
            Target element size **on and near the conductor surfaces** [m]
            (`GEO-15` step 1). ``None`` keeps the historical behaviour: one
            global ``setSize`` at `resolution` everywhere. When given, a gmsh
            Distance→Threshold background field grades from
            `conductor_resolution` at the conductor boundary up to `resolution`
            beyond `conductor_refine_distance`; air/box sizing is unchanged
            because the field's ``SizeMax`` *is* `resolution`.
        conductor_refine_distance : float, optional
            Width of the graded shell around the conductor [m]. Defaults to
            ``3 x ring_minor_radius``, which keeps the fine zone a thin skin
            rather than flooding the air box with small cells.
        phantom_resolution : float, optional
            `WF-6` step 3f₀. Target element size **inside the phantom
            cylinder** [m]. ``None`` (default) leaves the sizing exactly as it
            is — no field is created, so the mesh is bit-for-bit what it was.
            When given, a gmsh ``Box`` field over the phantom's own bounding
            box (``|z| <= phantom_height/2``, ``|x|,|y| <= phantom_radius``,
            plus a hair of margin so the cylinder's own boundary is inside)
            is set to `phantom_resolution` in and `resolution` out, and is
            combined with the `conductor_resolution` Threshold through a
            ``Min`` field, so neither knob can coarsen what the other
            refines. Setting it equal to `resolution` is a no-op by
            construction (the ``Min`` sees the same numbers) and is the
            negative control the step 3f₀ test asserts.
        return_diagnostics : bool
            When True, return a fourth element: a dict with the **CAD (occ)
            mass** per physical group — the junction-double-count-free
            denominator `GEO-15` measures against — and the gmsh meshing wall
            time. Broadcast from `rank`, so it is identical on every rank.

        Notes
        -----
        `n_legs`, `leg_radius`, and `leg_height` are accepted as backward-compatible
        aliases and map to `leg_count`, `leg_width/2`, and `coil_length` respectively.

        Cell tags:
        - 1: conductor (rings + legs)
        - 2: air
        - 3: phantom
        - 101..(100+leg_count): per-port regions between adjacent legs
        """
        if n_legs is not None:
            leg_count = n_legs
        if leg_radius is not None:
            leg_width = 2.0 * leg_radius
        if leg_height is not None:
            coil_length = leg_height

        if as_hole and not emit_port_sheets:
            raise ValueError(
                "as_hole needs emit_port_sheets: with the coil cut out of the "
                "mesh the port sheets are the only thing left to drive it"
            )
        if phantom_box is not None:
            # `GEO-38` step 2: checked before any gmsh state exists.
            # `GEO-36` step 8 (§9 item 530): the rod composes with the box;
            # its containment is tested against the box's half-extents.
            if symmetry_quadrant:
                raise ValueError("phantom_box cannot be combined with symmetry_quadrant")
            if phantom_resolution is not None:
                raise ValueError(
                    "phantom_box cannot be combined with phantom_resolution "
                    "(its Box field is cylinder-sized)"
                )
            try:
                box_size = tuple(float(s) for s in phantom_box["size_m"])
            except (KeyError, TypeError, ValueError) as exc:
                raise ValueError("phantom_box needs 'size_m': (Lx, Ly, Lz) [m]") from exc
            box_off = float(phantom_box.get("table_offset_m", 0.0))
            if len(box_size) != 3 or not all(np.isfinite(s) and s > 0.0 for s in box_size):
                raise ValueError(f"phantom_box size_m must be three finite > 0, got {box_size}")
            if not np.isfinite(box_off):
                raise ValueError(f"phantom_box table_offset_m must be finite, got {box_off}")
            phantom_box = {"size_m": box_size, "table_offset_m": box_off}
            phantom_radius = float(np.hypot(0.5 * box_size[0], 0.5 * box_size[1] + abs(box_off)))
            phantom_height = box_size[2]
        if implant_rod is not None:
            # `GEO-36` step 2. Contradictions in the request, checked before
            # any gmsh state exists.
            implant_rod = MeshGenerator._normalise_implant_rod(
                implant_rod,
                phantom_radius=phantom_radius,
                phantom_height=phantom_height,
                as_hole=as_hole,
                symmetry_quadrant=symmetry_quadrant,
                phantom_box=phantom_box,
            )
        if implant_spheroid is not None:
            # `GEO-36` step 7: the oblate body; refused before any gmsh state.
            implant_spheroid = MeshGenerator._normalise_implant_spheroid(
                implant_spheroid,
                phantom_radius=phantom_radius,
                phantom_height=phantom_height,
                as_hole=as_hole,
                symmetry_quadrant=symmetry_quadrant,
                implant_rod=implant_rod,
                phantom_box=phantom_box,
            )
        if not emit_leg_ports:
            # `GEO-35` step 1. Contradictions in the request, checked before
            # any gmsh state exists (a failed build deadlocks at rank > 1,
            # `GEO-23`).
            if leg_gap_length is not None:
                raise ValueError(
                    "emit_leg_ports=False requires leg_gap_length=None: in the "
                    "gapped layout the leg boxes are the drive"
                )
            if as_hole or symmetry_quadrant:
                raise ValueError(
                    "emit_leg_ports=False is not yet gated with as_hole or "
                    "symmetry_quadrant (`GEO-35` later steps)"
                )
        if leg_sheet_strip_fraction is not None:
            # `OPS-61` step 1b. Contradictions in the request.
            if not emit_port_sheets or leg_gap_length is None:
                raise ValueError(
                    "leg_sheet_strip_fraction requires emit_port_sheets=True " "and leg_gap_length"
                )
            if not (0.0 < float(leg_sheet_strip_fraction) < 1.0):
                raise ValueError(
                    "leg_sheet_strip_fraction must lie in (0, 1), got "
                    f"{leg_sheet_strip_fraction!r}"
                )
            if c4_congruent_sheets or symmetry_quadrant:
                raise ValueError(
                    "leg_sheet_strip_fraction refuses c4_congruent_sheets "
                    "and symmetry_quadrant (each needs one surface per "
                    "sheet or its own sheet construction)"
                )
        if c4_congruent_sheets:
            # `GEO-32`. Checked before any gmsh state exists: each is a
            # contradiction in the request, not a build failure.
            if not emit_port_sheets:
                raise ValueError(
                    "c4_congruent_sheets requires emit_port_sheets=True: it "
                    "constrains the port sheets' triangulation, and without "
                    "emit_port_sheets there are no sheets"
                )
            if leg_gap_length is None:
                raise ValueError(
                    "c4_congruent_sheets requires leg_gap_length: it binds the "
                    "leg-port sheets, which exist only in the gapped layout"
                )
            if leg_count != 4:
                raise ValueError(
                    "c4_congruent_sheets requires leg_count == 4: the "
                    f"constraint is the C4 copy of sheet 1, got leg_count={leg_count}"
                )
        if symmetry_quadrant:
            # `TH-16` step 3a. Contradictions in the request, checked before
            # any gmsh state exists. Two layouts are admitted (`TH-20` step
            # 0a): L4 (4 legs, leg gaps, solid) and LH (ring gaps, no leg
            # gaps, leg_count % 4 == 0, as_hole, longitudinal sheets).
            if not emit_port_sheets:
                raise ValueError("symmetry_quadrant requires emit_port_sheets=True")
            if c4_congruent_sheets or leg_azimuth_offsets_rad is not None:
                raise ValueError(
                    "symmetry_quadrant cannot be combined with "
                    "c4_congruent_sheets or leg_azimuth_offsets_rad"
                )
            if leg_gap_length is not None and ring_gap_length is not None:
                raise ValueError(
                    "symmetry_quadrant: ring and leg gaps together are not "
                    "gated (`TH-20` step 0a)"
                )
            if leg_gap_length is not None:
                if leg_count != 4 or as_hole:
                    raise ValueError(
                        "symmetry_quadrant with leg_gap_length requires "
                        "leg_count == 4 and as_hole=False"
                    )
            elif ring_gap_length is not None:
                if not as_hole:
                    raise ValueError(
                        "symmetry_quadrant: the solid ring-gap quadrant is not "
                        "gated (`TH-20` step 0a); use as_hole=True"
                    )
                if leg_count % 4 != 0:
                    raise ValueError(
                        "symmetry_quadrant with ring_gap_length requires "
                        f"leg_count % 4 == 0, got leg_count={leg_count}"
                    )
                if ring_sheet_orientation != "longitudinal":
                    raise ValueError(
                        "symmetry_quadrant: transverse ring sheets are not "
                        "gated (`TH-20` step 0a)"
                    )
            else:
                raise ValueError(
                    "symmetry_quadrant requires leg_gap_length (L4) or " "ring_gap_length (LH)"
                )
        if leg_count < 3:
            raise ValueError("leg_count must be >= 3 for a birdcage-like fixture")
        if leg_width <= 0.0:
            raise ValueError("leg_width must be > 0")
        if leg_spacing <= 0.0:
            raise ValueError("leg_spacing must be > 0")
        if coil_length <= 0.0:
            raise ValueError("coil_length must be > 0")
        if ring_radius <= 0.0:
            raise ValueError("ring_radius must be > 0")

        if emit_port_sheets and leg_gap_length is None and ring_gap_length is None:
            # Checked here, before any gmsh state exists: the mode is a
            # contradiction, not a build failure.
            raise ValueError(
                "emit_port_sheets requires leg_gap_length or ring_gap_length: "
                "the sheet is a gap box's mid-plane, and without a gap the port "
                "boxes float in the air with no terminals to span"
            )

        if ring_sheet_orientation not in ("transverse", "longitudinal"):
            raise ValueError(
                "ring_sheet_orientation must be 'transverse' or 'longitudinal', "
                f"got {ring_sheet_orientation!r}"
            )
        if ring_sheet_orientation == "longitudinal" and ring_gap_length is None:
            # The mode names a surface inside a ring gap; without ring gaps it
            # has no referent, and silently ignoring it would let a caller
            # believe it got longitudinal sheets on the leg ports.
            raise ValueError(
                "ring_sheet_orientation='longitudinal' requires ring_gap_length: "
                "it selects the ring gap's own sheet, and the leg-port sheets "
                "(`GEO-16`/`GEO-18`) are longitudinal already"
            )

        if leg_azimuth_offsets_rad is not None:
            # `PORT-9` step 3 leg (d1). Checked before gmsh exists, like the
            # `emit_port_sheets` contradiction above.
            offsets = np.asarray(leg_azimuth_offsets_rad, dtype=float)
            if offsets.shape != (leg_count,):
                raise ValueError(
                    "leg_azimuth_offsets_rad must carry one angle per leg: got "
                    f"shape {tuple(offsets.shape)} for leg_count={leg_count}"
                )
            if leg_gap_length is None:
                # Without a gap the port boxes sit at the *midpoint* azimuths
                # between legs, so "rotate the leg with its port" has no
                # referent; refusing is honest, silently moving legs out from
                # under unmoved boxes is not.
                raise ValueError(
                    "leg_azimuth_offsets_rad requires leg_gap_length: the "
                    "ungapped layout floats its port boxes between legs, so a "
                    "leg cannot carry its port around with it"
                )
            if ring_gap_length is not None:
                # `GEO-20` cuts each ring at the mid-azimuth between adjacent
                # legs. A displaced leg no longer bisects those arcs, so the
                # ring gaps would silently stop being centred on anything.
                raise ValueError(
                    "leg_azimuth_offsets_rad cannot be combined with "
                    "ring_gap_length: the ring gaps are centred on the "
                    "mid-azimuths of the *uniform* layout, which a displaced "
                    "leg no longer bisects"
                )

        leg_radius_eff = 0.5 * leg_width
        ring_port_layout: Optional[Dict[str, float]] = None
        if ring_gap_length is not None:
            ring_port_layout = MeshGenerator._birdcage_ring_gap_layout(
                leg_count=leg_count,
                ring_radius=ring_radius,
                leg_radius_eff=leg_radius_eff,
                leg_spacing=leg_spacing,
                ring_minor_radius=ring_minor_radius,
                phantom_radius=phantom_radius,
                ring_gap_length=ring_gap_length,
                port_clearance=port_clearance,
                min_port_face_area=min_port_face_area,
            )
        if leg_gap_length is None:
            port_diagnostics = MeshGenerator.birdcage_port_layout_diagnostics(
                leg_count=leg_count,
                ring_radius=ring_radius,
                leg_width=leg_width,
                ring_minor_radius=ring_minor_radius,
                phantom_radius=phantom_radius,
                port_box_size=port_box_size,
                port_clearance=port_clearance,
                min_port_face_area=min_port_face_area,
                min_port_center_separation=min_port_center_separation,
            )
        else:
            # `GEO-18` step 1. The gapped layout is a different problem from the
            # floating-box one `birdcage_port_layout_diagnostics` validates: the
            # box is *supposed* to meet metal (that is the point — terminals),
            # so its radial-clearance guards do not apply. What still has to
            # hold is checked here.
            port_box_size, port_diagnostics = MeshGenerator._birdcage_leg_gap_layout(
                leg_count=leg_count,
                ring_radius=ring_radius,
                leg_radius_eff=leg_radius_eff,
                leg_spacing=leg_spacing,
                coil_length=coil_length,
                ring_minor_radius=ring_minor_radius,
                phantom_radius=phantom_radius,
                leg_gap_length=leg_gap_length,
                port_clearance=port_clearance,
                min_port_face_area=min_port_face_area,
                min_port_center_separation=min_port_center_separation,
                leg_azimuth_offsets_rad=leg_azimuth_offsets_rad,
            )

        if shield_radius is not None:
            # `GEO-38` step 1a: checked before gmsh exists, on every rank.
            if symmetry_quadrant:
                raise ValueError("shield_radius cannot be combined with symmetry_quadrant")
            if leg_gap_length is not None:
                # Chain C is high-pass only; the gap-mode leg boxes are rotated
                # and would need their own corner formula (review 2026-10-15 10:30).
                raise ValueError("shield_radius cannot be combined with leg_gap_length")
            # Corner radii, not face centres (item 468 W2: an axis-aligned
            # leg box's outer corner sat at 0.37486 m > R_s on a table the
            # face-centre check passed). Leg boxes are axis-aligned at the
            # leg-pair mid-azimuths, exactly as `_build_birdcage_port_model`
            # places them; ring-gap wedges reach (R + w/2)/cos(alpha).
            r_port = float(port_diagnostics["port_radius_m"])
            half_bx = 0.5 * float(port_box_size[0])
            half_by = 0.5 * float(port_box_size[1])
            leg_box_corner = 0.0
            if emit_leg_ports:
                theta_s = np.linspace(0.0, 2.0 * np.pi, leg_count, endpoint=False)
                for k_s in range(leg_count):
                    a_s, b_s = theta_s[k_s], theta_s[(k_s + 1) % leg_count]
                    phi_s = np.arctan2(np.sin(a_s) + np.sin(b_s), np.cos(a_s) + np.cos(b_s))
                    leg_box_corner = max(
                        leg_box_corner,
                        float(
                            np.hypot(
                                abs(r_port * np.cos(phi_s)) + half_bx,
                                abs(r_port * np.sin(phi_s)) + half_by,
                            )
                        ),
                    )
            ring_box_corner = 0.0
            if ring_gap_length is not None:
                ring_box_w = 2.0 * ring_minor_radius + 2.0 * port_clearance
                ring_box_corner = (ring_radius + 0.5 * ring_box_w) / float(
                    np.cos(0.5 * ring_gap_length / ring_radius)
                )
            shield_need = (
                max(
                    ring_radius + max(leg_radius_eff, ring_minor_radius),
                    ring_box_corner,
                    leg_box_corner,
                )
                + port_clearance
            )
            if not float(shield_radius) > shield_need:
                raise ValueError(
                    f"shield_radius {float(shield_radius)!r} must exceed {shield_need!r} "
                    f"(coil, ring-gap box corner {ring_box_corner!r} and leg-box corner "
                    f"{leg_box_corner!r} + port_clearance)"
                )

        # The generator below can raise (overlapping facets, GEO-9 step 2b). Two
        # things must happen when it does, or the failure poisons the rest of the
        # process: gmsh must be finalized (otherwise it stays initialised and
        # mid-command, and every later `occ` call in this process is refused —
        # that is known-issues 7), and the *other* ranks must be told, or they
        # block forever in the collective `model_to_mesh` below that the raising
        # rank never reaches. Both were measured: 20260803T123116Z log, 5 failed
        # 2 passed in 3.16 s of pytest but harness exit 124 at the 180 s ceiling.
        build_error: Optional[BaseException] = None
        build_diagnostics: Optional[Dict[str, object]] = None
        if comm.rank == rank:
            try:
                build_diagnostics = MeshGenerator._build_birdcage_port_model(
                    leg_count=leg_count,
                    ring_radius=ring_radius,
                    leg_radius_eff=leg_radius_eff,
                    leg_spacing=leg_spacing,
                    coil_length=coil_length,
                    ring_minor_radius=ring_minor_radius,
                    phantom_radius=phantom_radius,
                    phantom_height=phantom_height,
                    port_box_size=port_box_size,
                    port_radius=port_diagnostics["port_radius_m"],
                    leg_gap_length=leg_gap_length,
                    ring_gap_length=ring_gap_length,
                    ring_sheet_orientation=ring_sheet_orientation,
                    leg_azimuth_offsets_rad=leg_azimuth_offsets_rad,
                    port_clearance=port_clearance,
                    emit_port_sheets=emit_port_sheets,
                    as_hole=as_hole,
                    c4_congruent_sheets=c4_congruent_sheets,
                    air_padding=air_padding,
                    resolution=resolution,
                    conductor_resolution=conductor_resolution,
                    conductor_refine_distance=conductor_refine_distance,
                    phantom_resolution=phantom_resolution,
                    step_export_path=step_export_path,
                    symmetry_quadrant=symmetry_quadrant,
                    leg_sheet_strip_fraction=leg_sheet_strip_fraction,
                    implant_rod=implant_rod,
                    emit_leg_ports=emit_leg_ports,
                    implant_spheroid=implant_spheroid,
                    shield_radius=shield_radius,
                    phantom_box=phantom_box,
                )
            except BaseException as exc:  # noqa: BLE001 — re-raised below, on every rank
                build_error = exc
                if gmsh.isInitialized():
                    gmsh.finalize()

        # Collective: every rank learns whether the geometry exists before any of
        # them enters `model_to_mesh`.
        if comm.bcast(build_error is not None, root=rank):
            if build_error is not None:
                raise build_error
            raise RuntimeError(f"birdcage_port_domain geometry generation failed on rank {rank}")

        # `GEO-24` step 2a, the same reason `two_torus_domain` was plumbed in
        # `PORT-1` step 3b-iv: gmshio's default partitioner is
        # `GhostMode.none`, which leaves `cells_ghost = 0` on every rank. The
        # port facets rebuilt below are *interior* to the mesh, so classifying
        # one needs the tag of the cell on both sides — and on a partition
        # boundary one of those cells lives on the other rank. `shared_facet`
        # is what makes that cell present as a ghost; it is also what a `dS`
        # integral over those facets needs. `max_facet_to_cell_links = 2` is
        # dolfinx 0.11's required second argument and the documented value for
        # a non-branching manifold mesh (`OPS-18` step 3).
        partitioner = dolfinx.mesh.create_cell_partitioner(dolfinx.mesh.GhostMode.shared_facet, 2)
        mesh, cell_tags, facet_tags = _model_to_mesh(
            gmsh.model, comm, rank, gdim=3, partitioner=partitioner
        )

        if comm.rank == rank:
            gmsh.finalize()

        if (implant_rod is not None and implant_rod.get("keep_cells", False)) or (
            implant_spheroid is not None and implant_spheroid.get("keep_cells", False)
        ):
            # `GEO-36` step 3 (and step 7's kept spheroid): the kept rod's wall (402) is interior — the
            # phantom (3) / implant (`BIRDCAGE_IMPLANT_CELL_TAG`) interface —
            # so it is rebuilt here from the cell tags, never as a gmsh group
            # (known-issues 9; see `_interface_facet_tags`).
            facet_tags = _interface_facet_tags(
                mesh,
                cell_tags,
                {BIRDCAGE_IMPLANT_SURFACE_TAG: (3, BIRDCAGE_IMPLANT_CELL_TAG)},
                existing=facet_tags,
            )

        if not return_diagnostics:
            return mesh, cell_tags, facet_tags

        # Collective: the CAD masses only exist on the building rank, and every
        # rank must see the same denominator or the ratio assertions disagree
        # across ranks (the rank-local trap `GEO-9` already paid once).
        diagnostics = dict(comm.bcast(build_diagnostics, root=rank))
        # The realised port geometry, so a caller never has to restate it: with
        # `leg_gap_length` the box size is *derived*, not the one passed in.
        diagnostics["port_box_size_m"] = tuple(float(s) for s in port_box_size)
        diagnostics["leg_gap_length_m"] = None if leg_gap_length is None else float(leg_gap_length)
        diagnostics["port_layout"] = dict(port_diagnostics)
        diagnostics["leg_azimuth_offsets_rad"] = (
            None
            if leg_azimuth_offsets_rad is None
            else tuple(float(a) for a in leg_azimuth_offsets_rad)
        )
        diagnostics["ring_gap_length_m"] = (
            None if ring_gap_length is None else float(ring_gap_length)
        )
        diagnostics["ring_port_layout"] = (
            None if ring_port_layout is None else dict(ring_port_layout)
        )
        diagnostics["ring_sheet_orientation"] = str(ring_sheet_orientation)
        diagnostics["leg_ports_emitted"] = bool(emit_leg_ports)
        if shield_radius is not None:
            # `GEO-38` step 1a re-pose: the refusal's own corner radii.
            diagnostics["leg_box_corner_radius_m"] = float(leg_box_corner)
            diagnostics["ring_box_corner_radius_m"] = float(ring_box_corner)
        if phantom_box is not None:
            diagnostics["phantom_box"] = dict(phantom_box)
            diagnostics["phantom_box_circumradius_m"] = float(phantom_radius)
        return mesh, cell_tags, facet_tags, diagnostics

    @staticmethod
    def _normalise_implant_spheroid(
        implant_spheroid: Dict[str, object],
        *,
        phantom_radius: float,
        phantom_height: float,
        as_hole: bool = False,
        symmetry_quadrant: bool = False,
        implant_rod: Optional[Dict[str, object]] = None,
        phantom_box: Optional[Dict[str, object]] = None,
    ) -> Dict[str, object]:
        """`GEO-36` step 7: validate an ``implant_spheroid`` request; no gmsh call.

        With ``phantom_box`` (the normalised box dict) the containment test is
        the box's half-extents (`IMP-6` step 4); without it, the cylinder test.

        An oblate PEC body of revolution: equatorial semi-axis ``semi_major``
        a [m], polar semi-axis ``semi_minor`` c [m] along the unit ``axis``,
        centred at ``centre``. One implant per build (``implant_rod`` refused
        beside it); ``symmetry_quadrant`` refused; under ``as_hole`` only the
        kept body (``keep_cells=True``), the rod's step-6 rule. The body must
        lie strictly inside the phantom (bounding-sphere test). Returns
        ``{"semi_major", "semi_minor", "centre", "axis", "keep_cells"}``.
        """
        if implant_rod is not None:
            raise ValueError(
                "implant_spheroid cannot be combined with implant_rod " "(one implant per build)"
            )
        if symmetry_quadrant:
            raise ValueError("implant_spheroid cannot be combined with symmetry_quadrant")
        keep_cells = implant_spheroid.get("keep_cells", False)
        if not isinstance(keep_cells, bool):
            raise ValueError(f"implant_spheroid 'keep_cells' must be a bool, got {keep_cells!r}")
        if as_hole and keep_cells is not True:
            raise ValueError(
                "implant_spheroid under as_hole needs 'keep_cells': True "
                "(a hole body under as_hole is out of scope)"
            )
        try:
            a = float(implant_spheroid["semi_major"])
            c = float(implant_spheroid["semi_minor"])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError("implant_spheroid needs 'semi_major' and 'semi_minor' [m]") from exc
        if not (np.isfinite(a) and np.isfinite(c)) or c <= 0.0 or a <= c:
            raise ValueError(
                "implant_spheroid needs finite semi_major > semi_minor > 0, got "
                f"semi_major={a!r} semi_minor={c!r}"
            )
        if "centre" not in implant_spheroid or "axis" not in implant_spheroid:
            raise ValueError("implant_spheroid needs 'centre' and 'axis'")
        centre = np.asarray(implant_spheroid["centre"], dtype=float)
        axis = np.asarray(implant_spheroid["axis"], dtype=float)
        if centre.shape != (3,) or axis.shape != (3,):
            raise ValueError("implant_spheroid centre and axis must be 3-vectors")
        if not (np.all(np.isfinite(centre)) and np.all(np.isfinite(axis))):
            raise ValueError("implant_spheroid centre and axis must be finite")
        norm = float(np.linalg.norm(axis))
        if norm == 0.0:
            raise ValueError("implant_spheroid axis must be non-zero")
        axis = axis / norm
        if phantom_box is not None:
            # `IMP-6` step 4 (item 482): the box's own half-extents, not its
            # xy circumradius (bounding sphere of radius a, any axis).
            bx, by, bz = (float(s) for s in phantom_box["size_m"])
            boff = float(phantom_box.get("table_offset_m", 0.0))
            if not (
                abs(float(centre[0])) + a < 0.5 * bx
                and abs(float(centre[1]) - boff) + a < 0.5 * by
                and abs(float(centre[2])) + a < 0.5 * bz
            ):
                raise ValueError(
                    "implant_spheroid is not strictly inside the phantom_box: "
                    f"centre {tuple(float(v) for v in centre)} with semi_major {a} "
                    f"against size_m=({bx}, {by}, {bz}), table_offset_m={boff}"
                )
        elif not (
            float(np.hypot(centre[0], centre[1])) + a < phantom_radius
            and abs(float(centre[2])) + a < 0.5 * phantom_height
        ):
            raise ValueError(
                "implant_spheroid is not strictly inside the phantom: centre "
                f"{tuple(float(v) for v in centre)} with semi_major {a} against "
                f"R={phantom_radius}, H={phantom_height}"
            )
        return {
            "semi_major": a,
            "semi_minor": c,
            "centre": tuple(float(v) for v in centre),
            "axis": tuple(float(v) for v in axis),
            "keep_cells": keep_cells,
        }

    @staticmethod
    def _normalise_implant_rod(
        implant_rod: Dict[str, object],
        *,
        phantom_radius: float,
        phantom_height: float,
        as_hole: bool = False,
        symmetry_quadrant: bool = False,
        phantom_box: Optional[Dict[str, object]] = None,
    ) -> Dict[str, object]:
        """`GEO-36` step 2: validate an ``implant_rod`` request; no gmsh call.

        Returns ``{"radius", "length", "centre", "axis"}`` with ``axis`` a unit
        vector.  The rod is the set of points within ``radius`` of the segment
        ``centre ± axis*length/2``; it is strictly inside the phantom cylinder
        iff both segment end points are at radial distance ``< R - r`` and at
        ``|z| < H/2 - r`` (the radial distance is convex along the segment, so
        the end points bound it).  With ``phantom_box`` set (`GEO-36` step 8)
        each body (rod, coat, tip ball) is tested against the box's own
        half-extents instead: exact, since the capsule is the convex hull of
        its end balls and the box is convex.
        """

        def _box_refuse(end: np.ndarray, rho: float, body: str) -> None:
            # `GEO-36` step 8 (§9 item 530): strict box containment of a ball
            # of radius ``rho`` at cap centre ``end``.
            bx, by, bz = (float(v) for v in phantom_box["size_m"])
            boff = float(phantom_box.get("table_offset_m", 0.0))
            if not (
                abs(float(end[0])) + rho < 0.5 * bx
                and abs(float(end[1]) - boff) + rho < 0.5 * by
                and abs(float(end[2])) + rho < 0.5 * bz
            ):
                raise ValueError(
                    "implant_rod is not strictly inside the phantom_box: "
                    f"{body} cap centre {tuple(float(v) for v in end)} with "
                    f"radius {rho} against size_m=({bx}, {by}, {bz}), "
                    f"table_offset_m={boff}"
                )

        # `GEO-36` step 6 (§9 item 372): a kept, uncoated rod composes with
        # ``as_hole`` (its wall bounds two meshed volumes, so the single-use
        # cavity rule cannot take it); the hole rod and the coat stay refused.
        if symmetry_quadrant or (
            as_hole
            and (
                implant_rod.get("keep_cells", False) is not True
                or implant_rod.get("coat_radius", None) is not None
            )
        ):
            raise ValueError(
                "implant_rod cannot be combined with symmetry_quadrant, nor "
                "with as_hole unless keep_cells is True and no coat_radius "
                "(the coil-cavity rule would swallow a hole rod's wall; the "
                "quadrant cut would split the rod)"
            )
        try:
            radius = float(implant_rod["radius"])
            length = float(implant_rod["length"])
        except (KeyError, TypeError) as exc:
            raise ValueError("implant_rod needs 'radius' and 'length' [m]") from exc
        if not (radius > 0.0) or not (length > 0.0):
            raise ValueError(
                f"implant_rod radius and length must be > 0, got radius={radius!r} "
                f"length={length!r}"
            )
        centre = np.asarray(implant_rod.get("centre", (0.0, 0.0, 0.0)), dtype=float)
        axis = np.asarray(implant_rod.get("axis", (0.0, 0.0, 1.0)), dtype=float)
        if centre.shape != (3,) or axis.shape != (3,):
            raise ValueError("implant_rod centre and axis must be 3-vectors")
        norm = float(np.linalg.norm(axis))
        if not (norm > 0.0):
            raise ValueError("implant_rod axis must be non-zero")
        axis = axis / norm
        # `GEO-36` step 3: additive; absent or False is the step-2 hole.
        keep_cells = implant_rod.get("keep_cells", False)
        if not isinstance(keep_cells, bool):
            raise ValueError(f"implant_rod 'keep_cells' must be a bool, got {keep_cells!r}")
        # `GEO-36` step 4: additive cap-Ball ``VIn`` [m]; absent / None keeps
        # today's ``radius / 4``.  Refused before gmsh: bool, non-number,
        # non-finite, <= 0, or > radius / 4 (the knob only refines).
        tip_size = implant_rod.get("tip_size", None)
        if tip_size is not None:
            if isinstance(tip_size, bool) or not isinstance(
                tip_size, (int, float, np.integer, np.floating)
            ):
                raise ValueError(f"implant_rod 'tip_size' must be a number [m], got {tip_size!r}")
            tip_size = float(tip_size)
            if not np.isfinite(tip_size) or not (tip_size > 0.0) or (tip_size > radius / 4.0):
                raise ValueError(
                    "implant_rod 'tip_size' must be finite and in "
                    f"(0, radius/4 = {radius / 4.0!r}] m, got {tip_size!r}"
                )
        for end in (centre - 0.5 * length * axis, centre + 0.5 * length * axis):
            if phantom_box is not None:
                _box_refuse(end, radius, "rod")
            elif not (
                float(np.hypot(end[0], end[1])) + radius < phantom_radius
                and abs(end[2]) + radius < 0.5 * phantom_height
            ):
                raise ValueError(
                    "implant_rod is not strictly inside the phantom: cap centre "
                    f"{tuple(float(v) for v in end)} with radius {radius} against "
                    f"R={phantom_radius}, H={phantom_height}"
                )
        # `GEO-36` step 5: additive shaft-only sleeve radius b [m]; absent /
        # None is no coat.  Refused before gmsh: bool, non-number, non-finite,
        # b <= radius, without ``keep_cells`` (a coat around a hole is out of
        # scope), or a sleeve not strictly inside the phantom.
        coat_radius = implant_rod.get("coat_radius", None)
        if coat_radius is not None:
            if isinstance(coat_radius, bool) or not isinstance(
                coat_radius, (int, float, np.integer, np.floating)
            ):
                raise ValueError(
                    "implant_rod 'coat_radius' must be a number [m], got " f"{coat_radius!r}"
                )
            coat_radius = float(coat_radius)
            if not np.isfinite(coat_radius) or not (coat_radius > radius):
                raise ValueError(
                    "implant_rod 'coat_radius' must be finite and > radius "
                    f"= {radius!r} m, got {coat_radius!r}"
                )
            if keep_cells is not True:
                raise ValueError(
                    "implant_rod 'coat_radius' needs 'keep_cells': True (a coat "
                    "around a hole is out of scope)"
                )
            for end in (centre - 0.5 * length * axis, centre + 0.5 * length * axis):
                if phantom_box is not None:
                    _box_refuse(end, coat_radius, "coat")
                elif not (
                    float(np.hypot(end[0], end[1])) + coat_radius < phantom_radius
                    and abs(end[2]) + coat_radius < 0.5 * phantom_height
                ):
                    raise ValueError(
                        "implant_rod coat is not strictly inside the phantom: "
                        f"cap centre {tuple(float(v) for v in end)} with "
                        f"coat_radius {coat_radius} against R={phantom_radius}, "
                        f"H={phantom_height}"
                    )
        # `IMP-3` I17 (§9 item 464): additive cap-Ball radius [m]; absent /
        # None keeps today's ``2 * radius``.  Refused before gmsh: bool,
        # non-number, non-finite, <= radius, > 4 * radius (runaway guard), or
        # a ball not strictly inside the phantom.
        tip_ball_radius = implant_rod.get("tip_ball_radius", None)
        if tip_ball_radius is not None:
            if isinstance(tip_ball_radius, bool) or not isinstance(
                tip_ball_radius, (int, float, np.integer, np.floating)
            ):
                raise ValueError(
                    "implant_rod 'tip_ball_radius' must be a number [m], got "
                    f"{tip_ball_radius!r}"
                )
            tip_ball_radius = float(tip_ball_radius)
            if (
                not np.isfinite(tip_ball_radius)
                or not (tip_ball_radius > radius)
                or tip_ball_radius > 4.0 * radius
            ):
                raise ValueError(
                    "implant_rod 'tip_ball_radius' must be finite and in "
                    f"(radius = {radius!r}, 4*radius = {4.0 * radius!r}] m, "
                    f"got {tip_ball_radius!r}"
                )
            for end in (centre - 0.5 * length * axis, centre + 0.5 * length * axis):
                if phantom_box is not None:
                    _box_refuse(end, tip_ball_radius, "tip ball")
                elif not (
                    float(np.hypot(end[0], end[1])) + tip_ball_radius < phantom_radius
                    and abs(end[2]) + tip_ball_radius < 0.5 * phantom_height
                ):
                    raise ValueError(
                        "implant_rod tip ball is not strictly inside the "
                        f"phantom: cap centre {tuple(float(v) for v in end)} "
                        f"with tip_ball_radius {tip_ball_radius} against "
                        f"R={phantom_radius}, H={phantom_height}"
                    )
        return {
            "radius": radius,
            "length": length,
            "centre": tuple(float(v) for v in centre),
            "axis": tuple(float(v) for v in axis),
            "keep_cells": keep_cells,
            "tip_size": tip_size,
            "coat_radius": coat_radius,
            "tip_ball_radius": tip_ball_radius,
        }

    @staticmethod
    def _build_birdcage_port_model(
        *,
        leg_count: int,
        ring_radius: float,
        leg_radius_eff: float,
        leg_spacing: float,
        coil_length: float,
        ring_minor_radius: float,
        phantom_radius: float,
        phantom_height: float,
        port_box_size: Tuple[float, float, float],
        port_radius: float,
        air_padding: float,
        resolution: float,
        leg_gap_length: Optional[float] = None,
        ring_gap_length: Optional[float] = None,
        ring_sheet_orientation: str = "transverse",
        leg_azimuth_offsets_rad: Optional[Sequence[float]] = None,
        port_clearance: float = 1.0e-3,
        emit_port_sheets: bool = False,
        as_hole: bool = False,
        c4_congruent_sheets: bool = False,
        conductor_resolution: Optional[float] = None,
        conductor_refine_distance: Optional[float] = None,
        phantom_resolution: Optional[float] = None,
        step_export_path: Optional[str] = None,
        symmetry_quadrant: bool = False,
        leg_sheet_strip_fraction: Optional[float] = None,
        implant_rod: Optional[Dict[str, object]] = None,
        emit_leg_ports: bool = True,
        implant_spheroid: Optional[Dict[str, object]] = None,
        shield_radius: Optional[float] = None,
        phantom_box: Optional[Dict[str, object]] = None,
    ) -> Dict[str, object]:
        """Build the birdcage gmsh model on the calling rank (see `birdcage_port_domain`)."""
        if implant_rod is not None:
            implant_rod = MeshGenerator._normalise_implant_rod(
                implant_rod,
                phantom_radius=phantom_radius,
                phantom_height=phantom_height,
                as_hole=as_hole,
                symmetry_quadrant=symmetry_quadrant,
                phantom_box=phantom_box,
            )
        if implant_spheroid is not None:
            implant_spheroid = MeshGenerator._normalise_implant_spheroid(
                implant_spheroid,
                phantom_radius=phantom_radius,
                phantom_height=phantom_height,
                as_hole=as_hole,
                symmetry_quadrant=symmetry_quadrant,
                implant_rod=implant_rod,
                phantom_box=phantom_box,
            )
        gmsh.initialize()
        gmsh.model.add("birdcage_port_domain")

        # Simplified conductor scaffold: two rings plus vertical legs.
        z_ring_offset = 0.5 * leg_spacing
        theta = np.linspace(0.0, 2.0 * np.pi, leg_count, endpoint=False)
        # Mid-azimuth between each adjacent leg pair — the gap centres, and the
        # same angles the ungapped port boxes have always floated at. Deliberately
        # off `theta`, the *nominal* layout: `birdcage_port_domain` refuses to
        # combine leg offsets with ring gaps, so the two never disagree here.
        phi_mid = theta + np.pi / leg_count
        # `PORT-9` step 3 leg (d1): the legs and the leg ports (and only those)
        # ride the displaced azimuths. Adding an exact zero is exact, so `None`
        # and an all-zero vector both give back `theta` bit-for-bit — the
        # identity case is the same construction, not a zero-angle rotation.
        theta_leg = theta + (
            0.0
            if leg_azimuth_offsets_rad is None
            else np.asarray(leg_azimuth_offsets_rad, dtype=float)
        )

        ring_tags: List[int] = []
        # Each ring gap centre, in build order: bottom ring first, then top, so
        # port ordinal `leg_count + 1 + j` and `leg_count + 1 + leg_count + j`
        # are the mirror pair the top/bottom gate compares.
        ring_gap_centres: List[Tuple[float, float]] = []
        # Pappus on the ring primitives, before the fragment touches them: the
        # swept volume is `pi·r²·R·angle` whether the sweep is 2*pi or an arc.
        # Measured here so a mass discrepancy after the union can be attributed —
        # OCC's own quadrature on this primitive vs the arcs being the wrong
        # arcs. `GEO-18`'s leg cut had no equivalent (its removed segment is a
        # plain cylinder, integrated exactly).
        ring_analytic_mass = 0.0
        if ring_gap_length is None:
            for z_sign in (+1.0, -1.0):
                ring_tags.append(
                    gmsh.model.occ.addTorus(
                        0, 0, z_sign * z_ring_offset, ring_radius, ring_minor_radius
                    )
                )
                ring_analytic_mass += np.pi * ring_minor_radius**2 * ring_radius * 2.0 * np.pi
        else:
            # `GEO-20` step 1: `leg_count` arcs per ring instead of one torus.
            # The removed arcs are bounded by the radial half-planes
            # `phi = phi_mid +- alpha`, so each cut face is an exact disk of
            # area pi*r^2 and the removed mass is exactly pi*r^2*g per gap.
            # `addTorus`'s partial sweep starts on the +x half-plane, so build
            # the arc at phi = 0 and rotate its start into place.
            half_angle = 0.5 * ring_gap_length / ring_radius
            arc_sweep = 2.0 * np.pi / leg_count - 2.0 * half_angle
            for z_sign in (-1.0, +1.0):
                z_ring = z_sign * z_ring_offset
                for j in range(leg_count):
                    ring_gap_centres.append((float(phi_mid[j]), float(z_ring)))
                    arc = gmsh.model.occ.addTorus(
                        0,
                        0,
                        z_ring,
                        ring_radius,
                        ring_minor_radius,
                        angle=arc_sweep,
                    )
                    gmsh.model.occ.rotate(
                        [(3, arc)],
                        0.0,
                        0.0,
                        0.0,
                        0.0,
                        0.0,
                        1.0,
                        float(phi_mid[j]) + half_angle,
                    )
                    ring_tags.append(arc)
                    ring_analytic_mass += np.pi * ring_minor_radius**2 * ring_radius * arc_sweep

        ring_cad_mass = float(sum(gmsh.model.occ.getMass(3, tag) for tag in ring_tags))

        leg_tags: List[int] = []
        for angle in theta_leg:
            x = ring_radius * np.cos(angle)
            y = ring_radius * np.sin(angle)
            if leg_gap_length is None:
                leg_tags.append(
                    gmsh.model.occ.addCylinder(
                        x, y, -0.5 * coil_length, 0.0, 0.0, coil_length, leg_radius_eff
                    )
                )
            else:
                # `GEO-18` step 1: two stubs instead of one leg. The removed
                # segment |z| <= g/2 is a plain cylinder (the layout helper has
                # already checked it stays clear of the rings), so the analytic
                # conductor mass is the uncut mass minus `leg_count` times
                # pi*r^2*g exactly, and each stub end face is a planar disk.
                stub_length = 0.5 * (coil_length - leg_gap_length)
                leg_tags.append(
                    gmsh.model.occ.addCylinder(
                        x, y, -0.5 * coil_length, 0.0, 0.0, stub_length, leg_radius_eff
                    )
                )
                leg_tags.append(
                    gmsh.model.occ.addCylinder(
                        x,
                        y,
                        0.5 * leg_gap_length,
                        0.0,
                        0.0,
                        stub_length,
                        leg_radius_eff,
                    )
                )

        if phantom_box is None:
            phantom_tag = gmsh.model.occ.addCylinder(
                0.0,
                0.0,
                -0.5 * phantom_height,
                0.0,
                0.0,
                phantom_height,
                phantom_radius,
            )
        else:
            # `GEO-38` step 2: the ASTM F2182 box in place of the cylinder.
            pb_x, pb_y, pb_z = phantom_box["size_m"]
            pb_off = float(phantom_box["table_offset_m"])
            phantom_tag = gmsh.model.occ.addBox(
                -0.5 * pb_x, -0.5 * pb_y + pb_off, -0.5 * pb_z, pb_x, pb_y, pb_z
            )
        # `GEO-38` step 2: ring-sheet surface per ordinal (pre-fragment tag).
        ring_sheet_tag_of_ordinal: Dict[int, int] = {}

        port_dx, port_dy, port_dz = port_box_size
        port_tags: List[int] = []
        # `GEO-18` step 2: one rectangle per port, the gap box's own mid-section
        # in the plane containing that leg's axis. Added as dim-2 *tools* to the
        # same fragment below so the tets conform to it — exactly the box's
        # cross-section, no larger, or the fragment would cut the conductor and
        # the air box along the same plane. Keyed by port ordinal: the plane as
        # a (unit normal, in-plane point) pair in the xy plane. `GEO-19` step B
        # replaced the old (pinned axis, coordinate) encoding, which could only
        # name a plane whose normal was x̂ or ŷ — i.e. only legs on a coordinate
        # axis, `leg_count <= 4`.
        sheet_tags: List[int] = []
        sheet_of_ordinal: Dict[int, Tuple[Tuple[float, float], Tuple[float, float]]] = {}
        # `GEO-20`: the ring gaps' mid-planes are radial, so their half-space
        # test is a `(normal, point)` pair in the gap's own frame — the
        # `GEO-18` lower/upper coordinate convention is not C_N-covariant and
        # must not be reused here.
        ring_sheet_of_ordinal: Dict[int, Tuple[np.ndarray, np.ndarray]] = {}
        n_ports_total = leg_count + (0 if ring_gap_length is None else 2 * leg_count)
        if emit_port_sheets:
            if leg_gap_length is None and ring_gap_length is None:
                raise ValueError(
                    "emit_port_sheets requires leg_gap_length or ring_gap_length: "
                    "the sheet is a gap box's mid-plane, and without a gap there "
                    "is no box on the metal to split"
                )
            if n_ports_total > 99:
                raise ValueError(
                    "emit_port_sheets encodes the port halves as cell tags "
                    f"100+i and 200+i, so the port count must be <= 99, got "
                    f"{n_ports_total}"
                )
        for idx, angle in enumerate(theta_leg):
            if not emit_leg_ports:
                # `GEO-35` step 1: no floating leg box (unbroken legs only).
                continue
            if leg_gap_length is None:
                next_angle = theta_leg[(idx + 1) % leg_count]
                box_angle = np.arctan2(
                    np.sin(angle) + np.sin(next_angle),
                    np.cos(angle) + np.cos(next_angle),
                )
            else:
                # On the leg, not between legs: the drive element of a low-pass
                # birdcage lives in the leg (`GEO-18`).
                box_angle = angle
            cx = port_radius * np.cos(box_angle)
            cy = port_radius * np.sin(box_angle)
            if leg_gap_length is None:
                port = gmsh.model.occ.addBox(
                    cx - 0.5 * port_dx,
                    cy - 0.5 * port_dy,
                    -0.5 * port_dz,
                    port_dx,
                    port_dy,
                    port_dz,
                )
            else:
                # `GEO-19` step B: on the leg, the box is the leg's *own* local
                # box — built at azimuth 0 and rotated about `ẑ` by the leg
                # azimuth, so its transverse faces stay normal to `φ̂` at any
                # count. It has to move with the sheet: a rotated mid-plane
                # rectangle of width `dx` only spans an axis-aligned square
                # section when the azimuth is a multiple of 90 degrees, and a
                # sheet that does not span its box leaves the box one piece,
                # i.e. no port halves at all. In gap mode the section *is* a
                # square (`_birdcage_leg_gap_layout` sets dx = dy = box_width),
                # so at the four axis azimuths this rotation maps the box onto
                # exactly the box the previous axis-aligned construction built.
                port = gmsh.model.occ.addBox(
                    port_radius - 0.5 * port_dx,
                    -0.5 * port_dy,
                    -0.5 * port_dz,
                    port_dx,
                    port_dy,
                    port_dz,
                )
                _rotate_about_z(3, port, float(box_angle))
            port_tags.append(port)

            if emit_port_sheets and leg_gap_length is not None:
                # The plane contains the leg axis *and* the radial direction, so
                # its normal is azimuthal, `φ̂ = (-sin θ, cos θ)`. `addRectangle`
                # only builds in the xy plane, so build the azimuth-0 rectangle
                # there and take it to the leg's mid-plane with the same
                # snapped rotation the box above takes — the sheet is that
                # box's mid-section by construction at every azimuth. The
                # rectangle is symmetric about both rotation axes, so neither
                # sign convention matters.
                ordinal = idx + 1
                sheet = gmsh.model.occ.addRectangle(
                    port_radius - 0.5 * port_dx, -0.5 * port_dz, 0.0, port_dx, port_dz
                )
                if leg_sheet_strip_fraction is None:
                    _place_sheet_in_leg_frame(sheet, float(box_angle))
                    sheet_tags.append(sheet)
                else:
                    # `OPS-61` step 1b: make the centre strip's edges CAD
                    # edges. Fragment at azimuth 0, *before* placement, so the
                    # cut is exact; every piece then takes the same placement.
                    f_strip = float(leg_sheet_strip_fraction)
                    centre = gmsh.model.occ.addRectangle(
                        port_radius - 0.5 * f_strip * port_dx,
                        -0.5 * port_dz,
                        0.0,
                        f_strip * port_dx,
                        port_dz,
                    )
                    pieces, _ = gmsh.model.occ.fragment([(2, sheet)], [(2, centre)])
                    for _, piece in pieces:
                        _place_sheet_in_leg_frame(piece, float(box_angle))
                        sheet_tags.append(piece)
                cos_t, _, _, _, sin_t, *_ = _z_rotation_affine(float(box_angle))
                sheet_of_ordinal[ordinal] = (
                    (-sin_t, cos_t),
                    (port_radius * cos_t, port_radius * sin_t),
                )

        # `GEO-20` step 1: one port solid per ring gap, in the gap's own frame.
        # Cross-section in that frame is the trapezoid bounded by the two radial
        # half-planes and by u = ring_radius +- w/2 (u = rho·cos(phi - phi_c)),
        # extruded through w in z. Six planar faces — so volume, surface and the
        # mid-plane section are exact under a linear mesh, which a constant-rho
        # (curved) face would not be. Corners are evaluated directly in global
        # coordinates rather than built at phi = 0 and rotated: the rotation is
        # what puts the C4 images ulps apart (`GEO-19` ruling (4*)).
        if ring_gap_length is not None:
            half_angle = 0.5 * ring_gap_length / ring_radius
            tan_half = float(np.tan(half_angle))
            ring_w = 2.0 * ring_minor_radius + 2.0 * port_clearance
            for phi_c, z_ring in ring_gap_centres:
                cos_c, sin_c = float(np.cos(phi_c)), float(np.sin(phi_c))
                corners = []
                for u, v in (
                    (ring_radius - 0.5 * ring_w, -(ring_radius - 0.5 * ring_w) * tan_half),
                    (ring_radius + 0.5 * ring_w, -(ring_radius + 0.5 * ring_w) * tan_half),
                    (ring_radius + 0.5 * ring_w, +(ring_radius + 0.5 * ring_w) * tan_half),
                    (ring_radius - 0.5 * ring_w, +(ring_radius - 0.5 * ring_w) * tan_half),
                ):
                    corners.append(
                        gmsh.model.occ.addPoint(
                            u * cos_c - v * sin_c,
                            u * sin_c + v * cos_c,
                            z_ring - 0.5 * ring_w,
                        )
                    )
                lines = [gmsh.model.occ.addLine(corners[k], corners[(k + 1) % 4]) for k in range(4)]
                loop = gmsh.model.occ.addCurveLoop(lines)
                face = gmsh.model.occ.addPlaneSurface([loop])
                extruded = gmsh.model.occ.extrude([(2, face)], 0.0, 0.0, ring_w)
                solids = [tag for dim, tag in extruded if dim == 3]
                if len(solids) != 1:
                    raise RuntimeError(
                        "birdcage_port_domain: ring-gap port extrusion produced "
                        f"{len(solids)} solids, expected 1"
                    )
                port_tags.append(solids[0])

                if emit_port_sheets:
                    ordinal = len(port_tags)
                    if not emit_leg_ports:
                        # `GEO-35` step 1: ordinals 1..leg_count stay reserved.
                        ordinal += leg_count
                    if ring_sheet_orientation == "longitudinal":
                        # `GEO-26` step 1: the section *through* the gap, in the
                        # plane u = ring_radius, normal û(phi_c). It spans the
                        # gap chord `2·R·tan(alpha)` along phi_hat (the drive
                        # direction, where the transverse section below spans
                        # exactly 0 — `PORT-13` step 1) and `w` along ẑ. Its
                        # four edges lie on the box's two radial cap planes
                        # (v = ±u·tan(alpha), evaluated here at u = R) and its
                        # two z faces, so it spans the box and splits it into
                        # an inner (u < R) and an outer (u > R) half. Corners
                        # evaluated in global coordinates like the box above —
                        # a build-at-0-and-rotate puts the C_N images ulps
                        # apart (`GEO-19` ruling (4*)).
                        v_half = ring_radius * tan_half
                        sheet_corners = [
                            gmsh.model.occ.addPoint(
                                ring_radius * cos_c - v * sin_c,
                                ring_radius * sin_c + v * cos_c,
                                z,
                            )
                            for v, z in (
                                (-v_half, z_ring - 0.5 * ring_w),
                                (+v_half, z_ring - 0.5 * ring_w),
                                (+v_half, z_ring + 0.5 * ring_w),
                                (-v_half, z_ring + 0.5 * ring_w),
                            )
                        ]
                        sheet_normal = np.array([cos_c, sin_c, 0.0])
                    else:
                        sheet_corners = [
                            gmsh.model.occ.addPoint(u * cos_c, u * sin_c, z)
                            for u, z in (
                                (ring_radius - 0.5 * ring_w, z_ring - 0.5 * ring_w),
                                (ring_radius + 0.5 * ring_w, z_ring - 0.5 * ring_w),
                                (ring_radius + 0.5 * ring_w, z_ring + 0.5 * ring_w),
                                (ring_radius - 0.5 * ring_w, z_ring + 0.5 * ring_w),
                            )
                        ]
                        sheet_normal = np.array([-sin_c, cos_c, 0.0])
                    sheet_lines = [
                        gmsh.model.occ.addLine(sheet_corners[k], sheet_corners[(k + 1) % 4])
                        for k in range(4)
                    ]
                    sheet = gmsh.model.occ.addPlaneSurface(
                        [gmsh.model.occ.addCurveLoop(sheet_lines)]
                    )
                    sheet_tags.append(sheet)
                    ring_sheet_tag_of_ordinal[ordinal] = sheet
                    ring_sheet_of_ordinal[ordinal] = (
                        sheet_normal,
                        np.array([ring_radius * cos_c, ring_radius * sin_c, z_ring]),
                    )

        radial_extent = ring_radius + max(leg_radius_eff, ring_minor_radius) + port_dy + air_padding
        z_extent = (
            max(0.5 * coil_length, 0.5 * leg_spacing + ring_minor_radius, 0.5 * phantom_height)
            + air_padding
        )
        if symmetry_quadrant:
            # `TH-16` step 3a: the quarter box. The leg-port sheets are not
            # fragment tools here — the kept port-box halves' faces on y = 0
            # (P1) and x = 0 (P2) are the sheets — so the rectangles built
            # above are dropped before the fragment.
            air_tag = gmsh.model.occ.addBox(
                0.0,
                0.0,
                -z_extent,
                radial_extent,
                radial_extent,
                2.0 * z_extent,
            )
            if leg_gap_length is not None:
                # L4 only; in LH (`TH-20` step 0a) the ring sheets stay
                # fragment tools and `sheet_tags` is kept whole.
                gmsh.model.occ.remove([(2, s) for s in sheet_tags], recursive=True)
                sheet_tags = []
        elif shield_radius is not None:
            # `GEO-38` step 1a: the closed PEC can replaces the box.
            air_tag = gmsh.model.occ.addCylinder(
                0.0, 0.0, -z_extent, 0.0, 0.0, 2.0 * z_extent, float(shield_radius)
            )
            shield_air_cad_mass = float(gmsh.model.occ.getMass(3, air_tag))

            def _on_shield(bb) -> bool:
                # The can's lateral face spans x, y = ±R_s and z = ±z_extent
                # (the caps are flat at ±z_extent). 1e-5 absorbs the OCC bbox
                # enlargement measured on curved faces (`TH-20` step 0a).
                r_s = float(shield_radius)
                return (
                    all(abs(bb[i] + r_s) < 1.0e-5 and abs(bb[i + 3] - r_s) < 1.0e-5 for i in (0, 1))
                    and abs(bb[2] + z_extent) < 1.0e-5
                    and abs(bb[5] - z_extent) < 1.0e-5
                )

        else:
            air_tag = gmsh.model.occ.addBox(
                -radial_extent,
                -radial_extent,
                -z_extent,
                2.0 * radial_extent,
                2.0 * radial_extent,
                2.0 * z_extent,
            )

        conductor_tags = ring_tags + leg_tags

        rod_occ_mass = float("nan")
        rod_bbox: Optional[Tuple[float, ...]] = None
        # `GEO-36` step 7: one implant per build, so the rod's keep flag and
        # body tag (`rod_tag`) carry the spheroid's too; every rod line keyed
        # on them behaves as before when ``implant_spheroid`` is None.
        keep_rod_cells = bool(
            (implant_rod is not None and implant_rod.get("keep_cells", False))
            or (implant_spheroid is not None and implant_spheroid.get("keep_cells", False))
        )
        spheroid_seam_bbox: Optional[Tuple[float, ...]] = None
        rod_tag: Optional[int] = None
        # `GEO-36` step 5: the optional shaft-only sleeve (tag 404).
        coat_b: Optional[float] = (
            None
            if implant_rod is None or implant_rod.get("coat_radius") is None
            else float(implant_rod["coat_radius"])
        )
        coat_tag: Optional[int] = None
        coat_occ_mass = float("nan")
        if implant_rod is not None:
            # `GEO-36` step 2: the capped rod, cut first (the `as_hole`
            # sequence). Shaft + two spheres fused (the probe's construction),
            # then cut from the phantom (tool kept) and from the air box (tool
            # dropped). Cut from the phantom alone, the air box would fill the
            # hole in the fragment below.
            rod_r = float(implant_rod["radius"])
            rod_l = float(implant_rod["length"])
            rod_c = np.asarray(implant_rod["centre"], dtype=float)
            rod_a = np.asarray(implant_rod["axis"], dtype=float)
            end0 = rod_c - 0.5 * rod_l * rod_a
            end1 = rod_c + 0.5 * rod_l * rod_a
            shaft = gmsh.model.occ.addCylinder(*end0.tolist(), *(rod_l * rod_a).tolist(), rod_r)
            cap0 = gmsh.model.occ.addSphere(*end0.tolist(), rod_r)
            cap1 = gmsh.model.occ.addSphere(*end1.tolist(), rod_r)
            rod_dimtags, _ = gmsh.model.occ.fuse([(3, shaft)], [(3, cap0), (3, cap1)])
            gmsh.model.occ.synchronize()
            if len(rod_dimtags) != 1 or rod_dimtags[0][0] != 3:
                raise RuntimeError(
                    f"implant_rod: the shaft/cap fuse gave {rod_dimtags}, not one solid"
                )
            rod_occ_mass = float(gmsh.model.occ.getMass(*rod_dimtags[0]))
            # `GEO-36` step 5: the shaft-only sleeve's outer cylinder (radius
            # b between the two cap-centre planes) is cut from the phantom and
            # the air box beside the rod; the sleeve itself is cyl_b - rod.
            cut_tools = list(rod_dimtags)
            if coat_b is not None:
                cyl_b = gmsh.model.occ.addCylinder(
                    *end0.tolist(), *(rod_l * rod_a).tolist(), coat_b
                )
                cut_tools = cut_tools + [(3, cyl_b)]
            phantom_out, _ = gmsh.model.occ.cut(
                [(3, phantom_tag)], cut_tools, removeObject=True, removeTool=False
            )
            # `GEO-36` step 3: with ``keep_cells`` the rod solid survives this
            # cut too and enters the fragment below as its own tool (group
            # `BIRDCAGE_IMPLANT_CELL_TAG`); without it, today's hole.
            air_out, _ = gmsh.model.occ.cut(
                [(3, air_tag)],
                cut_tools,
                removeObject=True,
                removeTool=not keep_rod_cells,
            )
            gmsh.model.occ.synchronize()
            if keep_rod_cells:
                rod_tag = int(rod_dimtags[0][1])
            if coat_b is not None:
                coat_out, _ = gmsh.model.occ.cut(
                    [(3, cyl_b)], rod_dimtags, removeObject=True, removeTool=False
                )
                gmsh.model.occ.synchronize()
                if len(coat_out) != 1 or coat_out[0][0] != 3:
                    raise RuntimeError(
                        "implant_rod: cutting the rod from the coat cylinder gave "
                        f"{coat_out}, not one solid"
                    )
                coat_tag = int(coat_out[0][1])
                coat_occ_mass = float(gmsh.model.occ.getMass(3, coat_tag))
            if len(phantom_out) != 1 or len(air_out) != 1:
                raise RuntimeError(
                    "implant_rod: cutting the rod left "
                    f"{len(phantom_out)} phantom / {len(air_out)} air volumes, not 1 / 1"
                )
            phantom_tag = phantom_out[0][1]
            air_tag = air_out[0][1]
            pad = 1.0e-6
            rod_bbox = tuple(
                [float(min(end0[i], end1[i]) - rod_r - pad) for i in range(3)]
                + [float(max(end0[i], end1[i]) + rod_r + pad) for i in range(3)]
            )
        elif implant_spheroid is not None:
            # `GEO-36` step 7: the exact oblate body of revolution. A half
            # ellipse (two quarter arcs, pole -> rim -> pole, major axis local
            # x) closed by the axis segment, revolved 2*pi about local z, so
            # the rim circle is a model edge (no `dilate`: its B-spline body
            # misses 4*pi*a^2*c/3 by 4.37e-4, `IMP-6` step 1 W0). Then rotate
            # pi/4 about local z (seam between rim points), local z -> axis,
            # translate to centre; cut first exactly as the rod.
            sph_a = float(implant_spheroid["semi_major"])
            sph_c = float(implant_spheroid["semi_minor"])
            sph_ctr = np.asarray(implant_spheroid["centre"], dtype=float)
            sph_ax = np.asarray(implant_spheroid["axis"], dtype=float)
            occ = gmsh.model.occ
            p_o = occ.addPoint(0.0, 0.0, 0.0)
            p_top = occ.addPoint(0.0, 0.0, sph_c)
            p_rim = occ.addPoint(sph_a, 0.0, 0.0)
            p_bot = occ.addPoint(0.0, 0.0, -sph_c)
            arc_up = occ.addEllipseArc(p_top, p_o, p_rim, p_rim)
            arc_dn = occ.addEllipseArc(p_rim, p_o, p_rim, p_bot)
            axis_line = occ.addLine(p_bot, p_top)
            loop = occ.addCurveLoop([arc_up, arc_dn, axis_line])
            profile = occ.addPlaneSurface([loop])
            revolved = occ.revolve([(2, profile)], 0, 0, 0, 0, 0, 1, 2.0 * np.pi)
            body = [dt for dt in revolved if dt[0] == 3]
            # The revolve leaves the profile face behind as a free surface at
            # the origin (inside the phantom); drop it, as `IMP-6`'s
            # `oblate_rim="edge"` does — left in, the 3D boundary recovery
            # fails ("PLC Error: A segment and a facet intersect",
            # 20261013T201959Z_GEO-36-step7-W1a.log:2171).
            occ.remove([(2, profile)], recursive=True)
            occ.synchronize()
            if len(body) != 1:
                raise RuntimeError(f"implant_spheroid: the revolve gave {revolved}, not one solid")
            occ.rotate(body, 0, 0, 0, 0, 0, 1, 0.25 * np.pi)
            z_hat = np.array([0.0, 0.0, 1.0])
            rot_k = np.cross(z_hat, sph_ax)
            if float(np.linalg.norm(rot_k)) > 1.0e-12:
                rot_ang = float(np.arccos(np.clip(float(z_hat @ sph_ax), -1.0, 1.0)))
                occ.rotate(body, 0, 0, 0, *rot_k.tolist(), rot_ang)
            elif float(z_hat @ sph_ax) < 0.0:
                occ.rotate(body, 0, 0, 0, 1, 0, 0, np.pi)
            occ.translate(body, *sph_ctr.tolist())
            occ.synchronize()
            rod_dimtags = body
            rod_occ_mass = float(occ.getMass(*rod_dimtags[0]))
            # The seam: the body's curves other than the rim (any curve whose
            # bbox spans the polar extent).
            body_faces = gmsh.model.getBoundary(
                rod_dimtags, combined=False, oriented=False, recursive=False
            )
            body_curves = sorted(
                {
                    c
                    for d, c in gmsh.model.getBoundary(
                        body_faces, combined=False, oriented=False, recursive=False
                    )
                    if d == 1
                }
            )
            for cur in body_curves:
                cbb = gmsh.model.getBoundingBox(1, cur)
                print(
                    f"[birdcage-mesh] implant spheroid curve {cur} bbox "
                    f"{tuple(float(v) for v in cbb)}",
                    flush=True,
                )
                if (
                    spheroid_seam_bbox is None
                    and abs(
                        float(np.dot(np.asarray(cbb[3:]) - np.asarray(cbb[:3]), np.abs(sph_ax)))
                        - 2.0 * sph_c
                    )
                    < 1.0e-6
                ):
                    spheroid_seam_bbox = tuple(float(v) for v in cbb)
            phantom_out, _ = gmsh.model.occ.cut(
                [(3, phantom_tag)], rod_dimtags, removeObject=True, removeTool=False
            )
            air_out, _ = gmsh.model.occ.cut(
                [(3, air_tag)],
                rod_dimtags,
                removeObject=True,
                removeTool=not keep_rod_cells,
            )
            gmsh.model.occ.synchronize()
            if keep_rod_cells:
                rod_tag = int(rod_dimtags[0][1])
            if len(phantom_out) != 1 or len(air_out) != 1:
                raise RuntimeError(
                    "implant_spheroid: cutting the body left "
                    f"{len(phantom_out)} phantom / {len(air_out)} air volumes, not 1 / 1"
                )
            phantom_tag = phantom_out[0][1]
            air_tag = air_out[0][1]
            pad = 1.0e-6
            rod_bbox = tuple(
                [float(sph_ctr[i] - sph_a - pad) for i in range(3)]
                + [float(sph_ctr[i] + sph_a + pad) for i in range(3)]
            )

        # GEO-9 step 2b. The previous `occ.cut(..., removeTool=False)` carved the
        # air box around the tools but never booleaned the tools against *each
        # other*: the legs pierce both rings by construction, so those solids
        # overlapped and gmsh failed with "Invalid boundary mesh (overlapping
        # facets) on surface 3 surface 49" (known-issues 7). Fragment the box
        # against ALL tools at once, so every pairwise overlap becomes its own
        # conforming piece, and re-derive the physical groups from the fragment
        # out-map — fragment renumbers and reorders, so absolute tags from before
        # the call mean nothing afterwards.
        if as_hole:
            # `TH-15` step 3a, `two_torus_domain(as_hole=True)`'s sequence:
            #   air = cut(box, conductors)     # the cavity
            #   fragment(air, [phantom, port boxes, sheets])
            # `removeTool=False` keeps the ring/leg solids alive (and so their
            # own faces addressable) until the cavity wall has been read off the
            # *meshed* volumes; they are dropped from the model below, before
            # anything is meshed.  The tools overlap each other by construction
            # (every leg pierces both rings) — a cut subtracts their union, so
            # that is exactly the coil.
            air_pieces, _ = gmsh.model.occ.cut(
                [(3, air_tag)],
                [(3, tag) for tag in conductor_tags],
                removeObject=True,
                removeTool=False,
            )
            gmsh.model.occ.synchronize()
            if not air_pieces:
                raise RuntimeError(
                    "birdcage_port_domain(as_hole=True): cutting the coil out " "left no air volume"
                )
            if phantom_box is not None and shield_radius is not None:
                # `GEO-38` step 2: the conductor union's CAD mass, read as the
                # can minus the cut air (the cut subtracts the union, so the
                # leg∩ring overlaps count once).
                conductor_cad_mass = shield_air_cad_mass - float(
                    sum(gmsh.model.occ.getMass(d, t) for d, t in air_pieces)
                )
                if implant_rod is not None or implant_spheroid is not None:
                    # `IMP-6` step 4 (item 482): the implant cut already took
                    # the body out of this air, so the difference above read
                    # conductor + body; remove the body (definition repair of
                    # item 470's diagnostic, not a band).
                    conductor_cad_mass = conductor_cad_mass - rod_occ_mass
            tool_dimtags = (
                [(3, phantom_tag)]
                + [(3, tag) for tag in port_tags]
                + [(2, tag) for tag in sheet_tags]
            )
            if keep_rod_cells:
                # `GEO-36` step 6: the kept rod, once (`GEO-9`).
                tool_dimtags = tool_dimtags + [(3, rod_tag)]
            _, fragment_map = gmsh.model.occ.fragment(air_pieces, tool_dimtags)
            gmsh.model.occ.synchronize()
            # out-map is positional: one entry per input, objects first then tools.
            input_dimtags = list(air_pieces) + tool_dimtags
        else:
            tool_tags = conductor_tags + [phantom_tag] + port_tags
            if keep_rod_cells:
                # `GEO-36` step 3: the kept rod, once (entering it twice is the
                # `GEO-9` overlapping-facet failure).
                tool_tags = tool_tags + [rod_tag]
            if coat_tag is not None:
                # `GEO-36` step 5: the sleeve, once, beside the rod.
                tool_tags = tool_tags + [coat_tag]
            tool_dimtags = [(3, tag) for tag in tool_tags] + [(2, tag) for tag in sheet_tags]
            _, fragment_map = gmsh.model.occ.fragment(
                [(3, air_tag)],
                tool_dimtags,
            )
            gmsh.model.occ.synchronize()

            # out-map is positional: one entry per input, objects first then tools.
            input_dimtags = [(3, air_tag)] + tool_dimtags
        if len(fragment_map) != len(input_dimtags):
            raise RuntimeError(
                "birdcage_port_domain: occ.fragment returned an out-map of "
                f"{len(fragment_map)} entries for {len(input_dimtags)} inputs"
            )

        removed_pieces: set = set()
        if symmetry_quadrant:
            # `TH-16` step 3a: the fragment is the full one (every tool, the
            # quarter box); every dim-3 piece outside the kept quadrant is
            # dropped here, before `ancestors`. The side test is the piece's
            # centroid, not its bounding box: OCC bounding boxes are enlarged
            # by the shape tolerance (≈ 1e-7 m), so a kept piece touching a cut
            # plane reports x_min < 0 and a bbox test at -1e-9 would drop it.
            # After the fragment every piece lies wholly on one side of each
            # cut face, so the centroid sign is the side.
            # `TH-20` step 0a: under `as_hole` the retained conductor solids
            # are dim-3 entities too; they are dropped wholesale after the
            # cavity read, and the ring arc on θ = 0 has a centroid y of order
            # 1e-17, so they are skipped here.
            hole_tools = set(conductor_tags) if as_hole else set()
            for _, piece in gmsh.model.getEntities(dim=3):
                if piece in hole_tools:
                    continue
                cx, cy, _ = gmsh.model.occ.getCenterOfMass(3, piece)
                if cx < 0.0 or cy < 0.0:
                    removed_pieces.add(piece)
            gmsh.model.occ.remove([(3, p) for p in sorted(removed_pieces)], recursive=True)
            gmsh.model.occ.synchronize()
            kept_surfaces = {
                surf
                for _, vol in gmsh.model.getEntities(dim=3)
                for d, surf in gmsh.model.getBoundary([(3, vol)], oriented=False, recursive=False)
                if d == 2
            }
            orphan_surfaces = [
                s for _, s in gmsh.model.getEntities(dim=2) if s not in kept_surfaces
            ]
            if orphan_surfaces:
                gmsh.model.occ.remove([(2, s) for s in orphan_surfaces], recursive=True)
                gmsh.model.occ.synchronize()

        # Only dim-3 inputs are ancestors of a dim-3 piece; the sheets are dim-2
        # tools, and gmsh numbers surfaces and volumes in separate spaces, so a
        # sheet tag entering `ancestors` could collide with a volume tag.
        ancestors: Dict[int, set] = {}
        for (input_dim, input_tag), pieces in zip(input_dimtags, fragment_map):
            if input_dim != 3:
                continue
            for dim, piece in pieces:
                if dim == 3 and piece not in removed_pieces:
                    ancestors.setdefault(piece, set()).add(input_tag)

        # With `as_hole` no piece descends from a conductor — they were cut away
        # before the fragment — so the metal-wins rule below is simply empty.
        conductor_set = set() if as_hole else set(conductor_tags)
        port_ordinal = {tag: idx for idx, tag in enumerate(port_tags, start=1)}
        if not emit_leg_ports:
            # `GEO-35` step 1: explicit map, ring ordinals unchanged at
            # leg_count+1 .. 3·leg_count; 1..leg_count are reserved gaps, so
            # `group_names` / `expected_groups` skip them below.
            port_ordinal = {tag: idx for idx, tag in enumerate(port_tags, start=leg_count + 1)}

        # Piece policy. Metal wins over everything (a leg∩ring piece is conductor
        # either way, and a port box grazing a conductor is metal, not an
        # integration volume); the phantom outranks a port box for the same
        # reason. Only pieces descended from a port box alone carry `100+i`, and
        # what is left — descended from the air box only — is air.
        group_of_piece: Dict[int, int] = {}
        for piece, sources in ancestors.items():
            if sources & conductor_set:
                group_of_piece[piece] = 1
            elif keep_rod_cells and rod_tag in sources:
                # `GEO-36` step 3: the rod was cut out of the phantom and the
                # air box, so its pieces descend from the rod alone.
                group_of_piece[piece] = BIRDCAGE_IMPLANT_CELL_TAG
            elif coat_tag is not None and coat_tag in sources:
                # `GEO-36` step 5: the sleeve was cut out of the phantom and
                # the air box too, so its pieces descend from the sleeve alone.
                group_of_piece[piece] = BIRDCAGE_IMPLANT_COAT_CELL_TAG
            elif phantom_tag in sources:
                group_of_piece[piece] = 3
            elif sources & port_ordinal.keys():
                ordinal = min(port_ordinal[tag] for tag in sources if tag in port_ordinal)
                if ordinal in sheet_of_ordinal:
                    # `GEO-18` step 2 / `GEO-16`: with the sheet in, each port box
                    # is two pieces, one either side of its own mid-plane. The
                    # sheet is a dim-2 input and so never an ancestor; the halves
                    # are told apart by the signed projection of the centroid on
                    # the plane's own normal, the same centroid discipline the
                    # two-torus branch uses. `GEO-19` step B: the projection is
                    # C4-covariant where the old x-or-y coordinate test was not,
                    # so "upper" is now the `+φ̂` side at *every* azimuth. On the
                    # four axis azimuths that relabels two of the four ports
                    # (the old rule always took the positive coordinate axis);
                    # nothing downstream reads the labels' sense — the sheet is
                    # the interface *between* the pair and the drive polarity is
                    # `ẑ` — so the relabel is inert.
                    (nx, ny), (px, py) = sheet_of_ordinal[ordinal]
                    centre = gmsh.model.occ.getCenterOfMass(3, piece)
                    value = (centre[0] - px) * nx + (centre[1] - py) * ny
                    group_of_piece[piece] = 200 + ordinal if value > 0.0 else 100 + ordinal
                elif ordinal in ring_sheet_of_ordinal:
                    # `GEO-20`: same discipline, but the sheet's plane is radial,
                    # so the side is the sign of (centroid - gap centre)·n. In
                    # `transverse` mode n = phi_hat and the halves are the two
                    # azimuthal sides; in `GEO-26`'s `longitudinal` mode n = û
                    # and they are the inner (u < R, 100+i) and outer (u > R,
                    # 200+i) halves. Nothing downstream reads the labels' sense.
                    normal, point = ring_sheet_of_ordinal[ordinal]
                    centre = np.asarray(gmsh.model.occ.getCenterOfMass(3, piece))
                    group_of_piece[piece] = (
                        200 + ordinal
                        if float(np.dot(centre - point, normal)) > 0.0
                        else 100 + ordinal
                    )
                else:
                    group_of_piece[piece] = 100 + ordinal
            else:
                group_of_piece[piece] = 2

        pieces_by_group: Dict[int, List[int]] = {}
        for piece, group in sorted(group_of_piece.items()):
            pieces_by_group.setdefault(group, []).append(piece)

        # `TH-15` step 3a — the cavity wall, and the drop of the retained tools.
        # Identified from the *meshed* volumes alone (`pieces_by_group`): the
        # ring and leg solids are still in the model at this point, and testing a
        # face against a set that includes them is self-satisfying (step 0's
        # abort: 716 `nodes not attached to any tet`).  A cavity wall is a face
        # that bounds exactly ONE meshed volume and does not lie flat on an outer
        # wall — every air/phantom/port interface bounds two, the outer walls
        # bound one but lie on a wall.
        cavity_surfaces: List[int] = []
        cavity_cad_area = float("nan")
        if as_hole:

            def _face_on_wall(surf: int) -> bool:
                bb = gmsh.model.getBoundingBox(2, surf)
                if symmetry_quadrant and any(
                    # `TH-20` step 0a: the cut faces are single-use faces too;
                    # flat on x = 0 / y = 0 at the measured OCC bbox
                    # enlargement (`…T171233Z_TH-16-step3a-probe.log:36`).
                    abs(bb[i]) < 1.0e-5 and abs(bb[i + 3]) < 1.0e-5
                    for i in (0, 1)
                ):
                    return True
                if shield_radius is not None and _on_shield(bb):
                    return True
                return any(
                    abs(bb[i] - lo) < 1e-6 and abs(bb[i + 3] - lo) < 1e-6
                    for i, lo in (
                        (0, -radial_extent),
                        (0, radial_extent),
                        (1, -radial_extent),
                        (1, radial_extent),
                        (2, -z_extent),
                        (2, z_extent),
                    )
                )

            face_use: Dict[int, int] = {}
            for pieces in pieces_by_group.values():
                for piece in pieces:
                    for dim, surf in gmsh.model.getBoundary(
                        [(3, piece)], oriented=False, recursive=False
                    ):
                        if dim == 2:
                            face_use[surf] = face_use.get(surf, 0) + 1
            cavity_surfaces = sorted(
                s for s, n in face_use.items() if n == 1 and not _face_on_wall(s)
            )
            if not cavity_surfaces:
                raise RuntimeError(
                    "birdcage_port_domain(as_hole=True): the meshed volumes have "
                    "no single-use non-wall boundary face, so the conductor "
                    "cavity surface is missing"
                )
            cavity_cad_area = float(sum(gmsh.model.occ.getMass(2, s) for s in cavity_surfaces))

            # Drop the conductor solids: the cavity IS the hole.
            # `recursive=False` so the cavity faces (which belong to the meshed
            # volumes, not to the tools) are untouched either way.
            gmsh.model.occ.remove([(3, t) for t in conductor_tags], recursive=False)
            gmsh.model.occ.synchronize()
            still_present = [t for t in conductor_tags if (3, t) in gmsh.model.getEntities(dim=3)]
            if still_present:
                # occ.remove + synchronize did not propagate; drop them
                # model-side (the two-torus route needed this fallback).
                gmsh.model.removeEntities([(3, t) for t in still_present], recursive=False)
            cavity_surfaces = [
                s for s in cavity_surfaces if (2, s) in gmsh.model.getEntities(dim=2)
            ]
            if not cavity_surfaces:
                raise RuntimeError(
                    "birdcage_port_domain(as_hole=True): every cavity surface "
                    "was destroyed with the conductor solids"
                )
            gmsh.model.addPhysicalGroup(2, cavity_surfaces, tag=BIRDCAGE_CONDUCTOR_SURFACE_TAG)
            gmsh.model.setPhysicalName(2, BIRDCAGE_CONDUCTOR_SURFACE_TAG, "conductor_surface")
            print(
                "[birdcage-mesh] conductor cavity (CAD) tag "
                f"{BIRDCAGE_CONDUCTOR_SURFACE_TAG}: {len(cavity_surfaces)} "
                f"surface(s) area={cavity_cad_area:.9e}; volumes after the drop="
                f"{len(gmsh.model.getEntities(dim=3))} "
                f"(conductor solids removed: {len(conductor_tags)})",
                flush=True,
            )

        implant_rod_cad: Optional[Dict[str, object]] = None
        rod_surfaces: List[int] = []
        implant_kept_mass = float("nan")
        implant_spheroid_cad: Optional[Dict[str, object]] = None
        if implant_rod is not None or implant_spheroid is not None:
            # `GEO-36` step 2 (step 7: the spheroid shares it): the rod wall from the *meshed* volumes' boundary
            # (the `TH-15` step 0 trap): a face bounding exactly one meshed
            # volume whose bounding box lies inside the rod's padded one.
            rod_face_use: Dict[int, int] = {}
            for pieces in pieces_by_group.values():
                for piece in pieces:
                    for dim, surf in gmsh.model.getBoundary(
                        [(3, piece)], oriented=False, recursive=False
                    ):
                        if dim == 2:
                            rod_face_use[surf] = rod_face_use.get(surf, 0) + 1
            if keep_rod_cells:
                # `GEO-36` step 3: the wall is the boundary of the kept
                # `"implant"` pieces, every face of which is shared with the
                # phantom — a two-use face, so the single-use filter below
                # would find nothing. A face of any other use count means the
                # rod touches something other than the phantom: a finding.
                implant_pieces = pieces_by_group.get(BIRDCAGE_IMPLANT_CELL_TAG, [])
                if not implant_pieces:
                    raise RuntimeError(
                        "implant_rod(keep_cells=True): occ.fragment left no " "implant piece"
                    )
                rod_surfaces = sorted(
                    {
                        surf
                        for piece in implant_pieces
                        for dim, surf in gmsh.model.getBoundary(
                            [(3, piece)], oriented=False, recursive=False
                        )
                        if dim == 2
                    }
                )
                bad_use = {
                    s: rod_face_use.get(s, 0) for s in rod_surfaces if rod_face_use.get(s, 0) != 2
                }
                if bad_use:
                    raise RuntimeError(
                        "implant_rod(keep_cells=True): rod wall faces not "
                        f"bounded by exactly two meshed volumes: {bad_use}"
                    )
                implant_kept_mass = float(sum(gmsh.model.occ.getMass(3, p) for p in implant_pieces))
            for surf, uses in sorted(rod_face_use.items()):
                if keep_rod_cells or uses != 1:
                    continue
                bb = gmsh.model.getBoundingBox(2, surf)
                if all(bb[i] >= rod_bbox[i] for i in range(3)) and all(
                    bb[i + 3] <= rod_bbox[i + 3] for i in range(3)
                ):
                    rod_surfaces.append(surf)
            if not rod_surfaces:
                raise RuntimeError(
                    "implant_rod: the meshed volumes have no single-use face "
                    "inside the rod's bounding box, so the rod wall is missing"
                )
            rod_areas = [float(gmsh.model.occ.getMass(2, s)) for s in rod_surfaces]
            rod_wall_area = float(sum(rod_areas))
            rod_centroid = np.zeros(3)
            for s, a in zip(rod_surfaces, rod_areas):
                rod_centroid += a * np.asarray(gmsh.model.occ.getCenterOfMass(2, s))
            rod_centroid /= rod_wall_area
            if not keep_rod_cells:
                gmsh.model.addPhysicalGroup(2, rod_surfaces, tag=BIRDCAGE_IMPLANT_SURFACE_TAG)
                gmsh.model.setPhysicalName(2, BIRDCAGE_IMPLANT_SURFACE_TAG, "implant_surface")
            # With ``keep_cells`` the wall is an *interior* facet group, and a
            # dim-2 gmsh group on interior facets hangs `model_to_mesh` at
            # `-n 2` (`_interface_facet_tags`' docstring, known-issues 9), so
            # tag 402 is rebuilt dolfinx-side from the 3 / 403 cell tags in
            # `birdcage_port_domain` instead — the module's interior-facet rule.
            if implant_spheroid is not None:
                sph_a = float(implant_spheroid["semi_major"])
                sph_c = float(implant_spheroid["semi_minor"])
                implant_spheroid_cad = {
                    "surfaces": list(rod_surfaces),
                    "wall_area_m2": rod_wall_area,
                    "body_mass_m3": rod_occ_mass,
                    "wall_centroid_m": tuple(float(v) for v in rod_centroid),
                    "facet_tag": BIRDCAGE_IMPLANT_SURFACE_TAG,
                    "spheroid": dict(implant_spheroid),
                    "keep_cells": keep_rod_cells,
                    "wall_size_m": sph_c / 3.0,
                    "rim_size_m": sph_c * sph_c / sph_a,
                    "seam_bbox_m": spheroid_seam_bbox,
                }
                if keep_rod_cells:
                    implant_spheroid_cad["cell_tag"] = BIRDCAGE_IMPLANT_CELL_TAG
                    implant_spheroid_cad["kept_mass_m3"] = implant_kept_mass
                print(
                    "[birdcage-mesh] implant spheroid (CAD) tag "
                    f"{BIRDCAGE_IMPLANT_SURFACE_TAG}: surfaces {rod_surfaces} "
                    f"wall area={rod_wall_area:.12e} body mass={rod_occ_mass:.12e} "
                    f"kept mass={implant_kept_mass:.12e} wall centroid="
                    f"({rod_centroid[0]:.12e}, {rod_centroid[1]:.12e}, "
                    f"{rod_centroid[2]:.12e}) seam bbox={spheroid_seam_bbox}",
                    flush=True,
                )
        if implant_rod is not None:
            implant_rod_cad = {
                "surfaces": list(rod_surfaces),
                "wall_area_m2": rod_wall_area,
                "rod_mass_m3": rod_occ_mass,
                "wall_centroid_m": tuple(float(v) for v in rod_centroid),
                "facet_tag": BIRDCAGE_IMPLANT_SURFACE_TAG,
                "rod": dict(implant_rod),
                "keep_cells": keep_rod_cells,
                # `GEO-36` step 4: the cap-Ball sizing actually used.
                "tip_size_m": (
                    float(implant_rod["radius"]) / 4.0
                    if implant_rod.get("tip_size") is None
                    else float(implant_rod["tip_size"])
                ),
                # `IMP-3` I17 (§9 item 464): the cap-Ball radius actually used.
                "tip_ball_radius_m": (
                    2.0 * float(implant_rod["radius"])
                    if implant_rod.get("tip_ball_radius") is None
                    else float(implant_rod["tip_ball_radius"])
                ),
            }
            if keep_rod_cells:
                implant_rod_cad["cell_tag"] = BIRDCAGE_IMPLANT_CELL_TAG
                implant_rod_cad["kept_mass_m3"] = implant_kept_mass
            # `GEO-36` step 5: the sleeve (``None`` when unset).
            coat_pieces = pieces_by_group.get(BIRDCAGE_IMPLANT_COAT_CELL_TAG, [])
            if coat_b is not None and not coat_pieces:
                raise RuntimeError("implant_rod(coat_radius): occ.fragment left no coat piece")
            implant_rod_cad["coat_radius_m"] = coat_b
            implant_rod_cad["coat_cell_tag"] = (
                None if coat_b is None else BIRDCAGE_IMPLANT_COAT_CELL_TAG
            )
            implant_rod_cad["coat_kept_mass_m3"] = (
                None
                if coat_b is None
                else float(sum(gmsh.model.occ.getMass(3, p) for p in coat_pieces))
            )
            implant_rod_cad["coat_occ_mass_m3"] = None if coat_b is None else coat_occ_mass
            coat_wall_area = None
            coat_centroid = None
            if coat_b is not None:
                coat_faces = sorted(
                    {
                        surf
                        for piece in coat_pieces
                        for dim, surf in gmsh.model.getBoundary(
                            [(3, piece)], oriented=False, recursive=False
                        )
                        if dim == 2
                    }
                )
                coat_wall_area = float(sum(gmsh.model.occ.getMass(2, s) for s in coat_faces))
                c_mass = [float(gmsh.model.occ.getMass(3, p)) for p in coat_pieces]
                c_sum = np.zeros(3)
                for p, m in zip(coat_pieces, c_mass):
                    c_sum += m * np.asarray(gmsh.model.occ.getCenterOfMass(3, p))
                coat_centroid = tuple(float(v) for v in c_sum / sum(c_mass))
            implant_rod_cad["coat_wall_area_m2"] = coat_wall_area
            implant_rod_cad["coat_centroid_m"] = coat_centroid
            print(
                "[birdcage-mesh] implant rod (CAD) tag "
                f"{BIRDCAGE_IMPLANT_SURFACE_TAG}: surfaces {rod_surfaces} "
                f"wall area={rod_wall_area:.9e} rod mass={rod_occ_mass:.9e} "
                f"wall centroid=({rod_centroid[0]:.9e}, {rod_centroid[1]:.9e}, "
                f"{rod_centroid[2]:.9e})",
                flush=True,
            )

        # A port is split only if it actually has a sheet: with `ring_gap_length`
        # alone the floating leg boxes have no terminals and stay single pieces.
        sheeted_ordinals = sorted(set(sheet_of_ordinal) | set(ring_sheet_of_ordinal))
        group_names = {1: "conductor", 2: "air", 3: "phantom"}
        if keep_rod_cells:
            group_names[BIRDCAGE_IMPLANT_CELL_TAG] = "implant"
        if coat_b is not None:
            group_names[BIRDCAGE_IMPLANT_COAT_CELL_TAG] = "implant_coat"
        for idx in port_ordinal.values():
            has_sheet = idx in sheet_of_ordinal or idx in ring_sheet_of_ordinal
            group_names[100 + idx] = f"port_P{idx}_lower" if has_sheet else f"port_P{idx}"
            if has_sheet:
                group_names[200 + idx] = f"port_P{idx}_upper"

        if symmetry_quadrant:
            # `TH-16` step 3a: the quadrant keeps one half of P1 (y = 0 sheet)
            # and one of P2 (x = 0 sheet), and nothing of P3 / P4.
            sheeted_ordinals = [
                o
                for o in sheeted_ordinals
                if (100 + o) in pieces_by_group or (200 + o) in pieces_by_group
            ]
            # `TH-20` step 0a: classify each kept ordinal. Straddling (one
            # half kept) is legal only for L4 leg-gap ports; interior (both
            # halves kept) only for LH ring ports.
            quadrant_interior_ordinals: List[int] = []
            for o in sheeted_ordinals:
                both = (100 + o) in pieces_by_group and (200 + o) in pieces_by_group
                if both:
                    if leg_gap_length is not None or o not in ring_sheet_of_ordinal:
                        raise RuntimeError(
                            f"symmetry_quadrant: both halves of P{o} survived the cut"
                        )
                    quadrant_interior_ordinals.append(o)
                elif leg_gap_length is None or o not in sheet_of_ordinal:
                    raise RuntimeError(
                        f"symmetry_quadrant: only one half of P{o} survived the "
                        "cut, legal only for an L4 leg-gap port"
                    )
            if as_hole:
                quadrant_kept_ordinals = sorted(
                    o
                    for o in port_ordinal.values()
                    if (100 + o) in pieces_by_group or (200 + o) in pieces_by_group
                )
                expected_groups = [2, 3] + [
                    g
                    for o in quadrant_kept_ordinals
                    for g in ((100 + o, 200 + o) if o in quadrant_interior_ordinals else (100 + o,))
                ]
            else:
                quadrant_kept_ordinals = list(sheeted_ordinals)
                expected_groups = [1, 2, 3]
        else:
            expected_groups = ([2, 3] if as_hole else [1, 2, 3]) + [
                100 + i for i in port_ordinal.values()
            ]
            expected_groups += [200 + i for i in sheeted_ordinals]
        if keep_rod_cells:
            expected_groups.append(BIRDCAGE_IMPLANT_CELL_TAG)
        if coat_b is not None:
            expected_groups.append(BIRDCAGE_IMPLANT_COAT_CELL_TAG)
        missing_groups = [group_names[tag] for tag in expected_groups if tag not in pieces_by_group]
        volumes = gmsh.model.getEntities(dim=3)
        masses = {tag: gmsh.model.occ.getMass(3, tag) for _, tag in volumes}
        if missing_groups:
            raise RuntimeError(
                "birdcage_port_domain: occ.fragment left no piece for "
                f"{', '.join(missing_groups)}; {len(volumes)} volumes, "
                "per-volume masses [m^3]: "
                + ", ".join(f"tag {tag}: {mass:.6e}" for tag, mass in sorted(masses.items()))
            )

        for group, pieces in sorted(pieces_by_group.items()):
            gmsh.model.addPhysicalGroup(3, pieces, tag=group)
            gmsh.model.setPhysicalName(3, group, group_names[group])

        # Every 3-D entity must carry a marker; a cell with none is exactly what
        # `gmshio.py:118` asserts on. Checked against the model rather than
        # against the out-map, so a piece the out-map never mentioned cannot
        # hide (the `GEO-9` step-1 guard, same reasoning).
        grouped_volumes = set()
        for _, group_tag in gmsh.model.getPhysicalGroups(dim=3):
            grouped_volumes.update(gmsh.model.getEntitiesForPhysicalGroup(3, group_tag))
        ungrouped = sorted({tag for _, tag in volumes} - grouped_volumes)
        if ungrouped:
            raise RuntimeError(
                "birdcage_port_domain: 3-D entities carry no physical group: "
                f"{ungrouped}; {len(volumes)} volumes, masses [m^3]: "
                + ", ".join(f"tag {tag}: {masses[tag]:.6e}" for tag in ungrouped)
            )

        print(
            f"[birdcage-mesh] fragment volumes={len(volumes)} "
            + " ".join(
                f"{group_names[group]}={sum(masses[p] for p in pieces):.6e}" f"({len(pieces)}p)"
                for group, pieces in sorted(pieces_by_group.items())
            ),
            flush=True,
        )

        # `GEO-18` step 2 / `GEO-16`: the sheet is what the two halves of one port
        # box share — the mid-plane, and nothing else, because the halves touch
        # only there. OCC's own area and bounding box give the CAD denominator the
        # meshed reading is scored against, and the `PORT-9` step-2b effective
        # width `w = A/h` its extents.
        port_sheet_cad: Dict[str, object] = {}
        quadrant_layout_diag: Optional[Dict[str, object]] = None
        cut_surfaces: Dict[int, List[int]] = {BIRDCAGE_CUT_X0_TAG: [], BIRDCAGE_CUT_Y0_TAG: []}
        if symmetry_quadrant:
            # `TH-16` step 3a: a cut face is one whose bounding box is flat on
            # x = 0 (y = 0), read off the post-fragment, post-removal model.
            # The tolerance is the OCC bbox enlargement, measured: after the
            # fragment the kept P2 half's x = 0 face reports x in
            # [-5.000e-06, 5.000e-06] (`…T171233Z_TH-16-step3a-probe.log:36`),
            # so 1e-6 missed it; 1e-5 is still far below any other face's
            # thickness. The meshed facets are checked flat to 1e-12 m.
            plane_tol = 1.0e-5
            for _, surf in gmsh.model.getEntities(dim=2):
                bb = gmsh.model.getBoundingBox(2, surf)
                if abs(bb[0]) < plane_tol and abs(bb[3]) < plane_tol:
                    cut_surfaces[BIRDCAGE_CUT_X0_TAG].append(surf)
                elif abs(bb[1]) < plane_tol and abs(bb[4]) < plane_tol:
                    cut_surfaces[BIRDCAGE_CUT_Y0_TAG].append(surf)
            for ordinal in sheeted_ordinals:
                group = 100 + ordinal if (100 + ordinal) in pieces_by_group else 200 + ordinal
                if ordinal in quadrant_interior_ordinals:
                    # `TH-20` step 0a: an interior LH ring port's sheet is the
                    # face its two halves share (the full-mesh path).
                    on_cut = {
                        surf
                        for vol in pieces_by_group[200 + ordinal]
                        for _, surf in gmsh.model.getBoundary(
                            [(3, vol)], oriented=False, recursive=False
                        )
                    }
                else:
                    on_cut = set(cut_surfaces[BIRDCAGE_CUT_X0_TAG]) | set(
                        cut_surfaces[BIRDCAGE_CUT_Y0_TAG]
                    )
                shared = sorted(
                    {
                        surf
                        for vol in pieces_by_group[group]
                        for _, surf in gmsh.model.getBoundary(
                            [(3, vol)], oriented=False, recursive=False
                        )
                    }
                    & on_cut
                )
                if not shared:
                    raise RuntimeError(
                        f"symmetry_quadrant: port P{ordinal}'s kept half has no "
                        "face on a cut plane; its faces' bboxes: "
                        + "; ".join(
                            f"{surf}: "
                            + ",".join(f"{v:.3e}" for v in gmsh.model.getBoundingBox(2, surf))
                            for vol in pieces_by_group[group]
                            for _, surf in gmsh.model.getBoundary(
                                [(3, vol)], oriented=False, recursive=False
                            )
                        )
                    )
                boxes = np.array([gmsh.model.getBoundingBox(2, surf) for surf in shared])
                port_sheet_cad[f"P{ordinal}"] = {
                    "surfaces": len(shared),
                    "area_m2": float(sum(gmsh.model.occ.getMass(2, surf) for surf in shared)),
                    "extent_x_m": float(boxes[:, 3].max() - boxes[:, 0].min()),
                    "extent_y_m": float(boxes[:, 4].max() - boxes[:, 1].min()),
                    "extent_z_m": float(boxes[:, 5].max() - boxes[:, 2].min()),
                }
            for cut_tag, name in (
                (BIRDCAGE_CUT_X0_TAG, "symmetry_cut_x0"),
                (BIRDCAGE_CUT_Y0_TAG, "symmetry_cut_y0"),
            ):
                if not cut_surfaces[cut_tag]:
                    raise RuntimeError(f"symmetry_quadrant: no surface on {name}")
                gmsh.model.addPhysicalGroup(2, cut_surfaces[cut_tag], tag=cut_tag)
                gmsh.model.setPhysicalName(2, cut_tag, name)
            quadrant_layout_diag = {
                "layout": "L4" if leg_gap_length is not None else "LH",
                "kept_ordinals": list(quadrant_kept_ordinals),
                "interior_ordinals": list(quadrant_interior_ordinals),
                "cut_cad_area_m2": {
                    int(t): float(sum(gmsh.model.occ.getMass(2, s) for s in cut_surfaces[t]))
                    for t in (BIRDCAGE_CUT_X0_TAG, BIRDCAGE_CUT_Y0_TAG)
                },
            }
            print(
                "[birdcage-mesh] symmetry quadrant: removed "
                f"{len(removed_pieces)} pieces; cut x0 {len(cut_surfaces[BIRDCAGE_CUT_X0_TAG])} "
                f"surface(s) area={sum(gmsh.model.occ.getMass(2, s) for s in cut_surfaces[BIRDCAGE_CUT_X0_TAG]):.9e}; "
                f"cut y0 {len(cut_surfaces[BIRDCAGE_CUT_Y0_TAG])} surface(s) "
                f"area={sum(gmsh.model.occ.getMass(2, s) for s in cut_surfaces[BIRDCAGE_CUT_Y0_TAG]):.9e}; "
                "port sheets (CAD) "
                + " ".join(
                    f"{name}: area={d['area_m2']:.9e} extents=({d['extent_x_m']:.6e}, "
                    f"{d['extent_y_m']:.6e}, {d['extent_z_m']:.6e})"
                    for name, d in sorted(port_sheet_cad.items())
                ),
                flush=True,
            )
        elif sheeted_ordinals:
            for ordinal in sheeted_ordinals:
                shared = sorted(
                    {
                        surf
                        for vol in pieces_by_group[100 + ordinal]
                        for _, surf in gmsh.model.getBoundary(
                            [(3, vol)], oriented=False, recursive=False
                        )
                    }
                    & {
                        surf
                        for vol in pieces_by_group[200 + ordinal]
                        for _, surf in gmsh.model.getBoundary(
                            [(3, vol)], oriented=False, recursive=False
                        )
                    }
                )
                if not shared:
                    raise RuntimeError(
                        f"birdcage_port_domain: port P{ordinal}'s halves share no "
                        "surface; the fragment did not conform to the port-sheet "
                        "mid-plane"
                    )
                boxes = np.array([gmsh.model.getBoundingBox(2, surf) for surf in shared])
                port_sheet_cad[f"P{ordinal}"] = {
                    "surfaces": len(shared),
                    "area_m2": float(sum(gmsh.model.occ.getMass(2, surf) for surf in shared)),
                    "extent_x_m": float(boxes[:, 3].max() - boxes[:, 0].min()),
                    "extent_y_m": float(boxes[:, 4].max() - boxes[:, 1].min()),
                    "extent_z_m": float(boxes[:, 5].max() - boxes[:, 2].min()),
                }
            print(
                "[birdcage-mesh] port sheets (CAD) "
                + " ".join(
                    f"{name}: {d['surfaces']} surface(s) area={d['area_m2']:.9e} "
                    f"extents=({d['extent_x_m']:.6e}, {d['extent_y_m']:.6e}, "
                    f"{d['extent_z_m']:.6e})"
                    for name, d in sorted(port_sheet_cad.items())
                ),
                flush=True,
            )

        outer_boundary_surfaces = []
        # `TH-16` step 3a: the cut faces reach r_e, so the wall rule below
        # would tag them; they were classified first and are excluded here.
        cut_surface_set = set(cut_surfaces[BIRDCAGE_CUT_X0_TAG]) | set(
            cut_surfaces[BIRDCAGE_CUT_Y0_TAG]
        )
        for dim, surf in gmsh.model.getEntities(dim=2):
            if surf in cut_surface_set:
                continue
            x_min, y_min, z_min, x_max, y_max, z_max = gmsh.model.getBoundingBox(dim, surf)
            if shield_radius is not None:
                # `GEO-38` step 1a: the box rule's `resolution` window would
                # catch leg boxes ≈ 5.7 mm inside the can; the can's own faces only.
                bb = (x_min, y_min, z_min, x_max, y_max, z_max)
                if _on_shield(bb) or (
                    abs(z_max - z_min) < 1.0e-5
                    and min(abs(z_min - z_extent), abs(z_min + z_extent)) < 1.0e-5
                ):
                    outer_boundary_surfaces.append(surf)
                continue
            if (
                abs(abs(x_min) - radial_extent) < resolution
                or abs(abs(x_max) - radial_extent) < resolution
                or abs(abs(y_min) - radial_extent) < resolution
                or abs(abs(y_max) - radial_extent) < resolution
                or abs(abs(z_min) - z_extent) < resolution
                or abs(abs(z_max) - z_extent) < resolution
            ):
                outer_boundary_surfaces.append(surf)

        if outer_boundary_surfaces:
            gmsh.model.addPhysicalGroup(2, outer_boundary_surfaces, tag=1)
            gmsh.model.setPhysicalName(2, 1, "outer_boundary")

        gmsh.model.mesh.setSize(gmsh.model.getEntities(0), resolution)

        # `GEO-15` step 1 — conductor grading. `mesh.setSize` only binds
        # dimension-0 entities, and an OCC torus carries a single seam point, so
        # a per-point constraint cannot resolve a 0.004 m minor radius no matter
        # what value it is given. The size *field* is the mechanism that can: a
        # Distance field sampled over the conductor's own boundary surfaces fed
        # into a Threshold that is `conductor_resolution` at the metal and
        # relaxes to the untouched global `resolution` past `refine_distance`.
        # The three MeshSizeFrom* switches must be off or gmsh takes the minimum
        # of the field and the point/curvature constraints, which would silently
        # re-impose the coarse global size inside the shell.
        # `WF-6` step 3f₀ adds a second, independent size field (the phantom
        # Box), so the two are collected here and combined through a `Min`
        # below. With `phantom_resolution=None` the list holds at most the
        # conductor Threshold and the tail reduces to the single
        # `setAsBackgroundMesh(threshold_field)` this code has always made —
        # the `None` path is a true no-op, asserted at 0.000e+00.
        size_fields: List[int] = []
        if conductor_resolution is not None:
            refine_distance = (
                3.0 * ring_minor_radius
                if conductor_refine_distance is None
                else conductor_refine_distance
            )
            # With `as_hole` the metal is gone and its boundary *is* the cavity
            # wall, so the Distance field is sampled over exactly the same
            # surface set the solid route grades on — the grading is the one
            # thing that must not change between the two routes.
            conductor_surfaces = (
                list(cavity_surfaces)
                if as_hole
                else sorted(
                    {
                        surf
                        for dim, surf in gmsh.model.getBoundary(
                            [(3, piece) for piece in pieces_by_group[1]],
                            combined=True,
                            oriented=False,
                            recursive=False,
                        )
                        if dim == 2 and surf not in cut_surface_set
                    }
                )
            )
            if not conductor_surfaces:
                raise RuntimeError(
                    "birdcage_port_domain: conductor grading requested but the "
                    "conductor group has no boundary surfaces "
                    f"({len(pieces_by_group.get(1, ()))} pieces, as_hole={as_hole})"
                )
            distance_field = gmsh.model.mesh.field.add("Distance")
            gmsh.model.mesh.field.setNumbers(distance_field, "SurfacesList", conductor_surfaces)
            gmsh.model.mesh.field.setNumber(distance_field, "Sampling", 20)
            threshold_field = gmsh.model.mesh.field.add("Threshold")
            gmsh.model.mesh.field.setNumber(threshold_field, "InField", distance_field)
            gmsh.model.mesh.field.setNumber(threshold_field, "SizeMin", conductor_resolution)
            gmsh.model.mesh.field.setNumber(threshold_field, "SizeMax", resolution)
            gmsh.model.mesh.field.setNumber(threshold_field, "DistMin", 0.0)
            gmsh.model.mesh.field.setNumber(threshold_field, "DistMax", refine_distance)
            size_fields.append(threshold_field)
            print(
                f"[birdcage-mesh] conductor grading: h_c={conductor_resolution:.6e} m "
                f"over {len(conductor_surfaces)} surfaces, shell={refine_distance:.6e} m, "
                f"h_global={resolution:.6e} m",
                flush=True,
            )

        # `WF-6` step 3f₀ — phantom sizing. A Box field is the cheapest thing
        # that can hold: the phantom is an axis-aligned cylinder, so its own
        # bounding box contains it exactly and nothing else of the geometry
        # (the nearest metal is the legs at r = 0.07 m, well outside a 0.03 m
        # phantom). The margin is one part in a thousand of the extent so the
        # cylinder's own boundary points fall strictly inside VIn rather than
        # on the Box face. VOut is the untouched global `resolution`, so with
        # `phantom_resolution == resolution` the field is numerically constant
        # at `resolution` and the `Min` below cannot move a single size.
        if phantom_resolution is not None:
            box_margin = 1.0e-3 * max(phantom_radius, 0.5 * phantom_height)
            phantom_field = gmsh.model.mesh.field.add("Box")
            gmsh.model.mesh.field.setNumber(phantom_field, "VIn", phantom_resolution)
            gmsh.model.mesh.field.setNumber(phantom_field, "VOut", resolution)
            gmsh.model.mesh.field.setNumber(phantom_field, "XMin", -phantom_radius - box_margin)
            gmsh.model.mesh.field.setNumber(phantom_field, "XMax", phantom_radius + box_margin)
            gmsh.model.mesh.field.setNumber(phantom_field, "YMin", -phantom_radius - box_margin)
            gmsh.model.mesh.field.setNumber(phantom_field, "YMax", phantom_radius + box_margin)
            gmsh.model.mesh.field.setNumber(
                phantom_field, "ZMin", -0.5 * phantom_height - box_margin
            )
            gmsh.model.mesh.field.setNumber(
                phantom_field, "ZMax", 0.5 * phantom_height + box_margin
            )
            gmsh.model.mesh.field.setNumber(phantom_field, "Thickness", 0.0)
            size_fields.append(phantom_field)
            print(
                f"[birdcage-mesh] phantom sizing: h_p={phantom_resolution:.6e} m over "
                f"r<={phantom_radius:.6e} m, |z|<={0.5 * phantom_height:.6e} m "
                f"(+{box_margin:.3e} m margin), h_global={resolution:.6e} m",
                flush=True,
            )

        if implant_rod is not None:
            # `GEO-36` step 2: the probe's grading (`tests/mesh/
            # probe_geo36_rod_in_phantom.py:114-133`), `RESOLUTION` -> the
            # generator's own `resolution`.
            rod_r = float(implant_rod["radius"])
            rod_l = float(implant_rod["length"])
            rod_c = np.asarray(implant_rod["centre"], dtype=float)
            rod_a = np.asarray(implant_rod["axis"], dtype=float)
            rod_distance = gmsh.model.mesh.field.add("Distance")
            gmsh.model.mesh.field.setNumbers(rod_distance, "SurfacesList", rod_surfaces)
            gmsh.model.mesh.field.setNumber(rod_distance, "Sampling", 50)
            rod_threshold = gmsh.model.mesh.field.add("Threshold")
            gmsh.model.mesh.field.setNumber(rod_threshold, "InField", rod_distance)
            gmsh.model.mesh.field.setNumber(rod_threshold, "SizeMin", rod_r / 3.0)
            gmsh.model.mesh.field.setNumber(rod_threshold, "SizeMax", resolution)
            gmsh.model.mesh.field.setNumber(rod_threshold, "DistMin", rod_r)
            gmsh.model.mesh.field.setNumber(rod_threshold, "DistMax", 10.0 * rod_r)
            size_fields.append(rod_threshold)
            # `GEO-36` step 4: the cap-Ball ``VIn`` knob (None -> r/4, today's).
            rod_tip = (
                rod_r / 4.0
                if implant_rod.get("tip_size") is None
                else float(implant_rod["tip_size"])
            )
            # `IMP-3` I17 (§9 item 464): the cap-Ball radius knob (None -> 2r).
            ball_r = (
                2.0 * rod_r
                if implant_rod.get("tip_ball_radius") is None
                else float(implant_rod["tip_ball_radius"])
            )
            for sign in (-1.0, 1.0):
                cap_centre = rod_c + sign * 0.5 * rod_l * rod_a
                ball = gmsh.model.mesh.field.add("Ball")
                gmsh.model.mesh.field.setNumber(ball, "XCenter", float(cap_centre[0]))
                gmsh.model.mesh.field.setNumber(ball, "YCenter", float(cap_centre[1]))
                gmsh.model.mesh.field.setNumber(ball, "ZCenter", float(cap_centre[2]))
                gmsh.model.mesh.field.setNumber(ball, "Radius", ball_r)
                gmsh.model.mesh.field.setNumber(ball, "VIn", rod_tip)
                gmsh.model.mesh.field.setNumber(ball, "VOut", resolution)
                size_fields.append(ball)
            print(
                f"[birdcage-mesh] implant rod sizing: r/3={rod_r / 3.0:.6e} m over "
                f"{len(rod_surfaces)} surfaces (r -> 10r), cap balls r/4="
                f"{rod_r / 4.0:.6e} m radius 2r, h_global={resolution:.6e} m, "
                f"tip_size(VIn)={rod_tip:.6e} m, ball radius={ball_r:.6e} m",
                flush=True,
            )
        elif implant_spheroid is not None:
            # `GEO-36` step 7: wall grading c/3 over c -> 10c (the rod's r/3
            # with r <-> c, the thin dimension) plus `IMP-6` step 1's rim
            # tube (MathEval distance to the rim circle); no Balls (item 408).
            sph_a = float(implant_spheroid["semi_major"])
            sph_c = float(implant_spheroid["semi_minor"])
            sph_ctr = np.asarray(implant_spheroid["centre"], dtype=float)
            sph_ax = np.asarray(implant_spheroid["axis"], dtype=float)
            rho_r = sph_c * sph_c / sph_a
            sph_distance = gmsh.model.mesh.field.add("Distance")
            gmsh.model.mesh.field.setNumbers(sph_distance, "SurfacesList", rod_surfaces)
            gmsh.model.mesh.field.setNumber(sph_distance, "Sampling", 50)
            sph_threshold = gmsh.model.mesh.field.add("Threshold")
            gmsh.model.mesh.field.setNumber(sph_threshold, "InField", sph_distance)
            gmsh.model.mesh.field.setNumber(sph_threshold, "SizeMin", sph_c / 3.0)
            gmsh.model.mesh.field.setNumber(sph_threshold, "SizeMax", resolution)
            gmsh.model.mesh.field.setNumber(sph_threshold, "DistMin", sph_c)
            gmsh.model.mesh.field.setNumber(sph_threshold, "DistMax", 10.0 * sph_c)
            size_fields.append(sph_threshold)
            # Plain-float `repr` (a numpy scalar's repr is `np.float64(...)`,
            # which mathex rejects — W1a's first run, 20261013T201740Z).
            cx0, cy0, cz0 = (repr(float(v)) for v in sph_ctr)
            ax0, ay0, az0 = (repr(float(v)) for v in sph_ax)
            dx = f"(x-({cx0}))"
            dy = f"(y-({cy0}))"
            dz = f"(z-({cz0}))"
            zl = f"({dx}*({ax0})+{dy}*({ay0})+{dz}*({az0}))"
            # rho_l = |d x axis| (a sum of squares, never negative).
            ax_, ay_, az_ = (repr(float(v)) for v in sph_ax)
            cx_ = f"({dy}*({az_})-{dz}*({ay_}))"
            cy_ = f"({dz}*({ax_})-{dx}*({az_}))"
            cz_ = f"({dx}*({ay_})-{dy}*({ax_}))"
            rl = f"Sqrt({cx_}^2+{cy_}^2+{cz_}^2)"
            rim_expr = f"Sqrt(({rl}-({float(sph_a)!r}))^2+{zl}^2)"
            rim_eval = gmsh.model.mesh.field.add("MathEval")
            gmsh.model.mesh.field.setString(rim_eval, "F", rim_expr)
            rim_threshold = gmsh.model.mesh.field.add("Threshold")
            gmsh.model.mesh.field.setNumber(rim_threshold, "InField", rim_eval)
            gmsh.model.mesh.field.setNumber(rim_threshold, "SizeMin", rho_r)
            gmsh.model.mesh.field.setNumber(rim_threshold, "SizeMax", resolution)
            gmsh.model.mesh.field.setNumber(rim_threshold, "DistMin", rho_r)
            gmsh.model.mesh.field.setNumber(rim_threshold, "DistMax", 4.0 * rho_r)
            size_fields.append(rim_threshold)
            print(
                f"[birdcage-mesh] implant spheroid sizing: c/3={sph_c / 3.0:.6e} m over "
                f"{len(rod_surfaces)} surfaces (c -> 10c), rim tube rho_r="
                f"{rho_r:.6e} m (rho_r -> 4 rho_r), h_global={resolution:.6e} m, "
                f"rim MathEval F={rim_expr}",
                flush=True,
            )

        if size_fields:
            if len(size_fields) == 1:
                background_field = size_fields[0]
            else:
                background_field = gmsh.model.mesh.field.add("Min")
                gmsh.model.mesh.field.setNumbers(background_field, "FieldsList", size_fields)
            gmsh.model.mesh.field.setAsBackgroundMesh(background_field)
            gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
            gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
            gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)

        if c4_congruent_sheets:
            # `GEO-32`. After every `occ.synchronize()` that follows the
            # fragment (the `as_hole` drop included) and before `generate(3)`.
            # The tags are the sheets' *post-fragment* surfaces, read off the
            # positional out-map: a pre-fragment tag names an entity the
            # fragment deleted. The leg sheets are the first `leg_count` dim-2
            # tools, in ordinal order. The affine is the snapped one, so at 90
            # degrees the copy is exact (last-bit coordinates move the cell
            # count at 1e-3, `GEO-19` step B). gmsh matches the boundary
            # curves under the affine itself; if they do not coincide it
            # fails, and that failure is a finding — no tolerance is set here.
            if not (emit_port_sheets and leg_gap_length is not None and leg_count == 4):
                raise ValueError(
                    "c4_congruent_sheets requires emit_port_sheets=True, "
                    "leg_gap_length and leg_count == 4"
                )
            post_fragment_sheets: List[int] = []
            for ordinal, sheet in enumerate(sheet_tags[:leg_count], start=1):
                pieces = fragment_map[input_dimtags.index((2, sheet))]
                surfaces = [tag for dim, tag in pieces if dim == 2]
                if len(surfaces) != 1:
                    raise RuntimeError(
                        f"c4_congruent_sheets: port P{ordinal}'s sheet became "
                        f"{len(surfaces)} surfaces in the fragment "
                        f"({surfaces}), so a one-to-one periodic copy is undefined"
                    )
                post_fragment_sheets.append(surfaces[0])
            master = post_fragment_sheets[0]
            for k in range(1, leg_count):
                gmsh.model.mesh.setPeriodic(
                    2,
                    [post_fragment_sheets[k]],
                    [master],
                    list(_z_rotation_affine(float(theta_leg[k] - theta_leg[0]))),
                )
            print(
                "[birdcage-mesh] c4_congruent_sheets: post-fragment sheet "
                f"surfaces {post_fragment_sheets} bound to master {master} "
                "under the snapped z-rotation",
                flush=True,
            )

        if step_export_path is not None:
            # `GEO-34` step 1a: the OCC model as built, no mesh in it yet.
            gmsh.write(str(step_export_path))

        mesh_start = time.perf_counter()
        gmsh.model.mesh.generate(3)
        gmsh.model.mesh.optimize("Netgen")
        mesh_wall_time = time.perf_counter() - mesh_start

        # CAD (occ) masses per physical group. `masses` above is per fragment
        # piece and is what `pieces_by_group` indexes, so summing it group-wise
        # gives the exact CAD volume of the group — the leg∩ring junctions are
        # separate pieces after the fragment, so each is counted once. That is
        # precisely the denominator the analytic ring+leg sum cannot supply.
        cad_mass_by_group = {
            group_names[group]: float(sum(masses[piece] for piece in pieces))
            for group, pieces in sorted(pieces_by_group.items())
        }
        print(
            f"[birdcage-mesh] mesh wall time {mesh_wall_time:.2f} s; CAD masses [m^3]: "
            + " ".join(f"{name}={mass:.9e}" for name, mass in cad_mass_by_group.items()),
            flush=True,
        )
        build_out = {
            "cad_mass_by_group": cad_mass_by_group,
            "mesh_wall_time_s": float(mesh_wall_time),
            "conductor_resolution_m": conductor_resolution,
            "phantom_resolution_m": phantom_resolution,
            "ring_cad_mass_m3": float(ring_cad_mass),
            "ring_analytic_mass_m3": float(ring_analytic_mass),
            "port_sheet_cad": port_sheet_cad,
            "as_hole": bool(as_hole),
            "symmetry_quadrant": bool(symmetry_quadrant),
            "symmetry_quadrant_layout": quadrant_layout_diag,
            "conductor_cavity_cad": (
                None
                if not as_hole
                else {
                    "surfaces": len(cavity_surfaces),
                    "area_m2": cavity_cad_area,
                    "facet_tag": BIRDCAGE_CONDUCTOR_SURFACE_TAG,
                }
            ),
            "implant_rod_cad": implant_rod_cad,
            "implant_spheroid_cad": implant_spheroid_cad,
        }
        if shield_radius is not None:
            # `GEO-38` step 1a diagnostics (absent on the default path).
            build_out["shield_radius_m"] = float(shield_radius)
            build_out["shield_length_m"] = float(2.0 * z_extent)
            build_out["air_domain_cad_mass_m3"] = shield_air_cad_mass
        if phantom_box is not None:
            # `GEO-38` step 2 diagnostics (absent on the default path):
            # `GEO-37`'s N_gmsh — triangles gmsh put on each post-fragment
            # ring sheet (positional out-map), read after `generate(3)`.
            n_gmsh: Dict[int, int] = {}
            for ordinal, sheet in sorted(ring_sheet_tag_of_ordinal.items()):
                pieces = fragment_map[input_dimtags.index((2, sheet))]
                n_tri = 0
                for dim, surf in pieces:
                    if dim != 2:
                        continue
                    _, elem_tags, _ = gmsh.model.mesh.getElements(2, surf)
                    n_tri += int(sum(len(t) for t in elem_tags))
                n_gmsh[int(ordinal)] = n_tri
            build_out["ring_sheet_gmsh_triangles"] = n_gmsh
            if as_hole and shield_radius is not None:
                build_out["conductor_cad_mass_m3"] = float(conductor_cad_mass)
        return build_out
