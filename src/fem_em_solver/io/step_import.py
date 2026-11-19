"""`GEO-34` step 1a — STEP import through Gmsh/OpenCASCADE with a group map.

``import_step(path, config, comm)`` reads a STEP file with
``gmsh.model.occ.importShapes``, restores shared faces with
``occ.removeAllDuplicates()`` (a fragment of everything, so two solids that
touch share one face again and the mesh is conforming), then claims **every**
imported volume for exactly one configured group by geometric selectors, tags
it with the group's physical tag, sizes it, meshes it and distributes it
through the project's ``_model_to_mesh``.  It returns the same
``(mesh, cell_tags, facet_tags, diagnostics)`` tuple the primitive generators
in ``io/mesh.py`` return.

Config (a dict, or a path to a JSON file)::

    {
      "occ_target_unit": null,          # optional Geometry.OCCTargetUnit
      "resolution": 0.015,              # global size [m] on every CAD point
      "outer_boundary_tag": 1,          # optional dim-2 group on the exterior boundary
      "volume_groups": [
        {"name": "conductor", "tag": 1,
         "selectors": [{"centroid_in_box": [x0, y0, z0, x1, y1, z1]}, ...],
         "grading": {"size": 0.0016, "distance": 0.012}},   # optional
        ...
      ],
      "congruent_sheets": {               # optional (`GEO-34` step 1b′)
        "pairs": [["port_P1_lower", "port_P1_upper"], ...],
        "rotation_deg": 90.0}
    }

Two additive keys (`GEO-34` step 1d, §9 item 985); absent, the import is
byte-identical:

* ``"highest_dim_only": false`` calls ``importShapes(path,
  highestDimOnly=False)`` so a free STEP face (its own product, an
  ``OPEN_SHELL``) is read at all — the default drops it (known-issues
  2026-10-08).
* ``"surface_groups": [{"name": "port_1", "step_name": "PortSheet",
  "role": "port", "tag": 201}, ...]`` resolves each ``step_name`` by
  ``getEntityName`` (exact, or the last ``/`` component; exactly one dim-2
  match, read after ``removeAllDuplicates``), embeds the named faces with
  ``occ.fragment(<all volumes>, <named faces>)``, traces each through the
  ``out_map`` (exactly one child, area equal to the pre-fragment ``getMass``
  to 1e-12, listed by ``mesh.getEmbedded`` of a volume) and tags the child
  natively.  Identity is name plus lineage, never a geometric selector.
  Diagnostics gain ``surface_groups[name] = {tag, role, step_name, surface,
  cad_area, n_triangles, route}``.  With named surfaces the Netgen
  optimisation pass is skipped (it SEGV'd on the embedded face, W1 of item
  985); ``generate(3)``'s default optimisation still runs.

``congruent_sheets`` binds pair k's shared face to pair 1's with
``gmsh.model.mesh.setPeriodic`` under the snapped z-rotation by
``rotation_deg·(k − 1)`` — `GEO-32`'s congruent cut on an imported model.
Absent, the import is unchanged.

A *selector* is a dict of predicates that must all hold; a group matches a
volume when any of its selectors does.  Predicates: ``centroid_in_box``
(OCC centre of mass inside the box), ``bbox_within`` (the volume's bounding
box inside the box), ``bbox_contains`` (the volume's bounding box contains
the box).  A volume matched by no group, a volume matched by two, and a group
that matches no volume are each a ``ValueError`` naming the offender — a
silent untagged port is the failure a config invites.

Grading mirrors ``birdcage_port_domain``'s conductor grading: a Distance field
over the group's boundary surfaces into a Threshold that is ``size`` at the
surface and the global ``resolution`` past ``distance``; the MeshSizeFrom*
switches are turned off exactly as the generator does, then
``generate(3)`` + ``optimize("Netgen")``.

gmsh runs on ``rank`` only (one session per call, finalised on the way out,
also on failure); a failure on the building rank is re-raised on every rank
with the same type and message before the collective ``_model_to_mesh``.
"""

from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import dolfinx
import gmsh
import numpy as np
from mpi4py import MPI

from fem_em_solver.io.mesh import _model_to_mesh, _z_rotation_affine

_PREDICATES = ("centroid_in_box", "bbox_within", "bbox_contains")


def load_step_config(config) -> Dict[str, object]:
    """A config dict, from a dict or a JSON file path."""
    if isinstance(config, (str, Path)):
        with open(config, "r", encoding="utf-8") as handle:
            return json.load(handle)
    return dict(config)


def _in_box(point, box) -> bool:
    x0, y0, z0, x1, y1, z1 = box
    x, y, z = point
    return x0 <= x <= x1 and y0 <= y <= y1 and z0 <= z <= z1


def _selector_matches(selector, centroid, bbox) -> bool:
    unknown = set(selector) - set(_PREDICATES)
    if unknown or not selector:
        raise ValueError(
            f"import_step: selector {selector!r} has unknown or no predicates "
            f"(known: {', '.join(_PREDICATES)})"
        )
    lo, hi = bbox[:3], bbox[3:]
    if "centroid_in_box" in selector and not _in_box(centroid, selector["centroid_in_box"]):
        return False
    if "bbox_within" in selector:
        box = selector["bbox_within"]
        if not (_in_box(lo, box) and _in_box(hi, box)):
            return False
    if "bbox_contains" in selector:
        box = selector["bbox_contains"]
        if not (_in_box(box[:3], bbox) and _in_box(box[3:], bbox)):
            return False
    return True


def claim_volumes(volumes, groups) -> Dict[str, List[int]]:
    """Map every volume to exactly one group, or raise ``ValueError``.

    ``volumes`` is ``{tag: (centroid, bbox)}``; ``groups`` the config's
    ``volume_groups`` list.  Pure function: no gmsh call, so the claim rules
    are checkable without a CAD model.
    """
    claimed: Dict[str, List[int]] = {g["name"]: [] for g in groups}
    unclaimed, doubly = [], []
    for tag, (centroid, bbox) in sorted(volumes.items()):
        owners = [
            g["name"]
            for g in groups
            if any(_selector_matches(s, centroid, bbox) for s in g["selectors"])
        ]
        if not owners:
            unclaimed.append((tag, centroid))
        elif len(owners) > 1:
            doubly.append((tag, centroid, owners))
        else:
            claimed[owners[0]].append(tag)
    if unclaimed:
        raise ValueError(
            "import_step: volume(s) claimed by no config group: "
            + "; ".join(
                f"volume {tag} centroid=({c[0]:.6e}, {c[1]:.6e}, {c[2]:.6e})"
                for tag, c in unclaimed
            )
        )
    if doubly:
        raise ValueError(
            "import_step: volume(s) claimed by more than one config group: "
            + "; ".join(
                f"volume {tag} centroid=({c[0]:.6e}, {c[1]:.6e}, {c[2]:.6e}) " f"groups {owners}"
                for tag, c, owners in doubly
            )
        )
    empty = [name for name, tags in claimed.items() if not tags]
    if empty:
        raise ValueError("import_step: config group(s) matched no volume: " + ", ".join(empty))
    return claimed


def _build(path, config) -> Dict[str, object]:
    """Import, claim, tag, size and mesh on the calling rank."""
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("step_import")
    unit = config.get("occ_target_unit")
    if unit:
        gmsh.option.setString("Geometry.OCCTargetUnit", str(unit))

    started = time.perf_counter()
    if "highest_dim_only" in config:
        gmsh.model.occ.importShapes(str(path), highestDimOnly=bool(config["highest_dim_only"]))
    else:
        gmsh.model.occ.importShapes(str(path))
    gmsh.model.occ.removeAllDuplicates()
    gmsh.model.occ.synchronize()
    surface_spec = config.get("surface_groups")
    named_surfaces = _embed_named_surfaces(surface_spec, config) if surface_spec else None
    import_wall_time = time.perf_counter() - started

    volumes = {}
    masses = {}
    for _, tag in gmsh.model.getEntities(dim=3):
        centroid = tuple(float(v) for v in gmsh.model.occ.getCenterOfMass(3, tag))
        bbox = tuple(float(v) for v in gmsh.model.getBoundingBox(3, tag))
        volumes[tag] = (centroid, bbox)
        masses[tag] = float(gmsh.model.occ.getMass(3, tag))
    model_bbox = tuple(float(v) for v in gmsh.model.getBoundingBox(-1, -1))
    print(
        f"[step-import] {path}: {len(volumes)} volumes in {import_wall_time:.2f} s, "
        f"model bbox={tuple(f'{v:.6e}' for v in model_bbox)}",
        flush=True,
    )
    for tag, (c, b) in sorted(volumes.items()):
        print(
            f"[step-import]   volume {tag}: mass={masses[tag]:.9e} "
            f"centroid=({c[0]:.6e}, {c[1]:.6e}, {c[2]:.6e}) "
            f"bbox=({b[0]:.4e}, {b[1]:.4e}, {b[2]:.4e}, {b[3]:.4e}, {b[4]:.4e}, {b[5]:.4e})",
            flush=True,
        )

    groups = config["volume_groups"]
    claimed = claim_volumes(volumes, groups)
    for group in groups:
        gmsh.model.addPhysicalGroup(3, claimed[group["name"]], tag=int(group["tag"]))
        gmsh.model.setPhysicalName(3, int(group["tag"]), group["name"])

    outer_tag = config.get("outer_boundary_tag")
    outer_surfaces: List[int] = []
    outer_bboxes: List[Tuple[float, ...]] = []
    outer_cad_area = 0.0
    if outer_tag is not None:
        # The exterior boundary is the combined boundary of every volume:
        # after removeAllDuplicates + synchronize a face shared by two volumes
        # appears twice and cancels, so what remains is the model's skin.  No
        # tolerance: the bbox-plane predicate this replaces compared
        # gmsh.model.getBoundingBox against 1e-9 * extent, but OCC pads
        # bounding boxes by the shape tolerance (±1.200001e-01 against the CAD
        # box's ±0.12, `20260930T201605Z_GEO-34-step2-W1.log:1848`), so it
        # matched no surface (known-issues 2026-09-30).  All dim-3 entities in
        # one call with recursive=False — a per-volume call cancels nothing.
        walls = {
            abs(s)
            for d, s in gmsh.model.getBoundary(
                gmsh.model.getEntities(3),
                combined=True,
                oriented=False,
                recursive=False,
            )
            if d == 2
        }
        if not walls:
            raise ValueError(
                f"import_step: outer_boundary_tag={outer_tag} is set but the "
                "model's exterior boundary is empty"
            )
        if named_surfaces:
            # W0b (`20261114T181729Z_GEO-34-U12p-W0b-read.log:182`) read the
            # embedded child outside the wall set; removing by tag keeps the
            # rule if another CAD route ever puts it there.
            walls -= {g["surface"] for g in named_surfaces.values()}
        outer_surfaces = sorted(walls)
        outer_bboxes = [
            tuple(float(v) for v in gmsh.model.getBoundingBox(2, s)) for s in outer_surfaces
        ]
        outer_cad_area = float(sum(gmsh.model.occ.getMass(2, s) for s in outer_surfaces))
        gmsh.model.addPhysicalGroup(2, outer_surfaces, tag=int(outer_tag))
        gmsh.model.setPhysicalName(2, int(outer_tag), "outer_boundary")
        print(
            f"[step-import] outer_boundary_tag={outer_tag}: "
            f"{len(outer_surfaces)} exterior surface(s) {outer_surfaces}, "
            f"CAD area {outer_cad_area:.12e} m^2",
            flush=True,
        )

    resolution = float(config["resolution"])
    gmsh.model.mesh.setSize(gmsh.model.getEntities(0), resolution)
    size_fields: List[int] = []
    for group in groups:
        grading = group.get("grading")
        if not grading:
            continue
        surfaces = sorted(
            {
                s
                for d, s in gmsh.model.getBoundary(
                    [(3, v) for v in claimed[group["name"]]],
                    combined=True,
                    oriented=False,
                    recursive=False,
                )
                if d == 2
            }
        )
        distance = gmsh.model.mesh.field.add("Distance")
        gmsh.model.mesh.field.setNumbers(distance, "SurfacesList", surfaces)
        gmsh.model.mesh.field.setNumber(distance, "Sampling", 20)
        threshold = gmsh.model.mesh.field.add("Threshold")
        gmsh.model.mesh.field.setNumber(threshold, "InField", distance)
        gmsh.model.mesh.field.setNumber(threshold, "SizeMin", float(grading["size"]))
        gmsh.model.mesh.field.setNumber(threshold, "SizeMax", resolution)
        gmsh.model.mesh.field.setNumber(threshold, "DistMin", 0.0)
        gmsh.model.mesh.field.setNumber(threshold, "DistMax", float(grading["distance"]))
        size_fields.append(threshold)
    if size_fields:
        if len(size_fields) == 1:
            background = size_fields[0]
        else:
            background = gmsh.model.mesh.field.add("Min")
            gmsh.model.mesh.field.setNumbers(background, "FieldsList", size_fields)
        gmsh.model.mesh.field.setAsBackgroundMesh(background)
        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)

    bound_sheets = _bind_congruent_sheets(config.get("congruent_sheets"), claimed)

    mesh_start = time.perf_counter()
    gmsh.model.mesh.generate(3)
    if named_surfaces:
        # Netgen's optimiser does not respect an embedded face: on this
        # fixture it reported "456 illegal tets" and the rank died with SEGV
        # (`20261114T181903Z_GEO-34-U12p-W1.log:110–116`).  With named
        # surfaces only generate(3)'s default tet optimisation runs; without
        # the key this branch is never taken.
        pass
    else:
        gmsh.model.mesh.optimize("Netgen")
    mesh_wall_time = time.perf_counter() - mesh_start

    if named_surfaces:
        for name, record in named_surfaces.items():
            _, elem_tags, _ = gmsh.model.mesh.getElements(2, record["surface"])
            record["n_triangles"] = int(sum(len(t) for t in elem_tags))
            print(
                f"[step-import] surface_groups {name}: surface {record['surface']} "
                f"({record['step_name']!r}) tag {record['tag']} role "
                f"{record['role']} cad_area={record['cad_area']!r} "
                f"n_triangles={record['n_triangles']} route={record['route']}",
                flush=True,
            )

    diagnostics = {
        "cad_mass_by_group": {
            g["name"]: float(sum(masses[v] for v in claimed[g["name"]])) for g in groups
        },
        "volumes_by_group": {g["name"]: len(claimed[g["name"]]) for g in groups},
        "tag_by_group": {g["name"]: int(g["tag"]) for g in groups},
        "model_bbox": model_bbox,
        "import_wall_time_s": float(import_wall_time),
        "mesh_wall_time_s": float(mesh_wall_time),
        "n_volumes": len(volumes),
    }
    if outer_tag is not None:
        diagnostics["outer_boundary_surfaces"] = outer_surfaces
        diagnostics["outer_boundary_bboxes"] = outer_bboxes
        diagnostics["outer_boundary_cad_area"] = outer_cad_area
    if bound_sheets is not None:
        diagnostics["congruent_sheet_surfaces"] = bound_sheets
    if named_surfaces is not None:
        diagnostics["surface_groups"] = named_surfaces
    return diagnostics


def _step_name_matches(entity_name: str, step_name: str) -> bool:
    return bool(entity_name) and (
        entity_name == step_name or entity_name.split("/")[-1] == step_name
    )


def _embed_named_surfaces(spec, config) -> Dict[str, Dict[str, object]]:
    """`GEO-34` step 1d: embed STEP-named faces and tag them natively.

    Identity is the STEP product name (``getEntityName``, read after
    ``removeAllDuplicates``, where W0b read it surviving) plus the topological
    lineage through ``occ.fragment``'s ``out_map`` — never a geometric
    selector.  Raises ``ValueError`` for an absent name (listing the names
    present), a non-unique match, a role other than ``port``, a tag collision,
    or a child that is not one face of unchanged area embedded in a volume.
    """
    outer = config.get("outer_boundary_tag")
    taken = set() if outer is None else {int(outer)}
    present = sorted({gmsh.model.getEntityName(d, t) for d, t in gmsh.model.getEntities()} - {""})
    resolved = []
    for group in spec:
        name, step_name = group["name"], group["step_name"]
        role, tag = group.get("role"), int(group["tag"])
        if role != "port":
            raise ValueError(
                f"import_step: surface group {name!r} has role {role!r}; only "
                "'port' is supported"
            )
        if tag in taken:
            raise ValueError(
                f"import_step: surface group {name!r} tag {tag} collides with "
                "another dim-2 tag in the config"
            )
        taken.add(tag)
        matches = [
            s
            for _, s in gmsh.model.getEntities(2)
            if _step_name_matches(gmsh.model.getEntityName(2, s), step_name)
        ]
        if not matches:
            raise ValueError(
                f"import_step: surface group {name!r}: no dim-2 entity named "
                f"{step_name!r} (names present: {present}); a free STEP face "
                'needs "highest_dim_only": false'
            )
        if len(matches) != 1:
            raise ValueError(
                f"import_step: surface group {name!r}: {len(matches)} dim-2 "
                f"entities named {step_name!r} ({matches}); exactly one required"
            )
        resolved.append((group, matches[0], float(gmsh.model.occ.getMass(2, matches[0]))))

    volumes = gmsh.model.getEntities(3)
    tools = [(2, s) for _, s, _ in resolved]
    _, out_map = gmsh.model.occ.fragment(volumes, tools)
    gmsh.model.occ.synchronize()
    embedded = {
        s
        for _, v in gmsh.model.getEntities(3)
        for d, s in gmsh.model.mesh.getEmbedded(3, v)
        if d == 2
    }
    records: Dict[str, Dict[str, object]] = {}
    for (group, source, cad_area), children in zip(resolved, out_map[len(volumes) :]):
        name = group["name"]
        faces = [c for d, c in children if d == 2]
        if len(faces) != 1:
            raise ValueError(
                f"import_step: surface group {name!r}: fragment split "
                f"{group['step_name']!r} into {len(faces)} faces ({children}); "
                "exactly one required"
            )
        child = faces[0]
        child_area = float(gmsh.model.occ.getMass(2, child))
        if abs(child_area - cad_area) > 1e-12:
            raise ValueError(
                f"import_step: surface group {name!r}: fragment child area "
                f"{child_area!r} != pre-fragment {cad_area!r}"
            )
        if child not in embedded:
            raise ValueError(
                f"import_step: surface group {name!r}: fragment child {child} of "
                f"{group['step_name']!r} is not embedded in any volume"
            )
        gmsh.model.addPhysicalGroup(2, [child], tag=int(group["tag"]))
        gmsh.model.setPhysicalName(2, int(group["tag"]), name)
        records[name] = {
            "tag": int(group["tag"]),
            "role": group["role"],
            "step_name": group["step_name"],
            "surface": int(child),
            "cad_area": cad_area,
            "n_triangles": None,
            "route": f"name -> surface {source} -> fragment child {child} (embedded)",
        }
    return records


def _shared_face(group_a: List[int], group_b: List[int]) -> List[int]:
    """Surfaces on the boundary of both volume sets (post-heal tags)."""

    def faces(vols):
        return {
            abs(s)
            for v in vols
            for d, s in gmsh.model.getBoundary(
                [(3, v)], combined=False, oriented=False, recursive=False
            )
            if d == 2
        }

    return sorted(faces(group_a) & faces(group_b))


def _bind_congruent_sheets(spec, claimed) -> Optional[List[int]]:
    """`GEO-34` step 1b′: bind each pair's shared face to the first pair's.

    ``spec`` is ``{"pairs": [[lower, upper], ...], "rotation_deg": θ}`` or
    ``None`` (the default: nothing happens, the import is byte-identical).
    Pair k's shared face — the one surface in the boundary of both groups'
    volumes, read *after* ``removeAllDuplicates`` so the tag names a live
    entity — is made a periodic copy of pair 1's under the snapped
    ``_z_rotation_affine(θ·(k − 1))`` (`GEO-32`'s congruent cut, the same
    call as ``io/mesh.py``'s ``c4_congruent_sheets``).  More or fewer than one
    shared face is a ``RuntimeError`` naming the pair.  gmsh matches the
    boundary curves under the affine itself; if they do not coincide it
    fails, and that failure is a finding — no tolerance is set here.
    """
    if spec is None:
        return None
    pairs = spec["pairs"]
    rotation = np.radians(float(spec["rotation_deg"]))
    sheets: List[int] = []
    for lower, upper in pairs:
        for name in (lower, upper):
            if name not in claimed:
                raise ValueError(
                    f"import_step: congruent_sheets names group {name!r}, which "
                    "is not a config volume group"
                )
        shared = _shared_face(claimed[lower], claimed[upper])
        if len(shared) != 1:
            raise RuntimeError(
                f"import_step: congruent_sheets pair ({lower}, {upper}) shares "
                f"{len(shared)} surfaces ({shared}), so a one-to-one periodic "
                "copy is undefined"
            )
        sheets.append(shared[0])
    master = sheets[0]
    for k in range(1, len(sheets)):
        gmsh.model.mesh.setPeriodic(
            2, [sheets[k]], [master], list(_z_rotation_affine(float(rotation * k)))
        )
    print(
        f"[step-import] congruent_sheets: post-heal sheet surfaces {sheets} bound "
        f"to master {master} under the snapped z-rotation "
        f"({float(spec['rotation_deg'])} deg steps)",
        flush=True,
    )
    return sheets


def import_step(
    path,
    config,
    comm: MPI.Intracomm = MPI.COMM_WORLD,
    rank: int = 0,
) -> Tuple[object, object, object, Dict[str, object]]:
    """``(mesh, cell_tags, facet_tags, diagnostics)`` from a STEP file + config.

    Collective over ``comm``.  ``path`` need only be readable on ``rank``.
    """
    config = load_step_config(config)
    build_error: Optional[BaseException] = None
    build_diagnostics = None
    if comm.rank == rank:
        try:
            build_diagnostics = _build(path, config)
        except BaseException as exc:  # noqa: BLE001 — re-raised below, on every rank
            build_error = exc
            if gmsh.isInitialized():
                gmsh.finalize()

    failure = comm.bcast(
        None if build_error is None else (isinstance(build_error, ValueError), str(build_error)),
        root=rank,
    )
    if failure is not None:
        if build_error is not None:
            raise build_error
        is_value_error, message = failure
        raise (ValueError if is_value_error else RuntimeError)(
            f"{message} (raised on rank {rank}; this is rank {comm.rank})"
        )

    partitioner = dolfinx.mesh.create_cell_partitioner(dolfinx.mesh.GhostMode.shared_facet, 2)
    try:
        mesh, cell_tags, facet_tags = _model_to_mesh(
            gmsh.model, comm, rank, gdim=3, partitioner=partitioner
        )
    finally:
        if comm.rank == rank and gmsh.isInitialized():
            gmsh.finalize()
    diagnostics = dict(comm.bcast(build_diagnostics, root=rank))
    return mesh, cell_tags, facet_tags, diagnostics
