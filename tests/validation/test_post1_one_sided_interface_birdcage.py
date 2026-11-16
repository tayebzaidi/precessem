"""`POST-1` step 7b′ — one-sided interface residuals, D2-gated, on the mirror-built birdcage.

PROJECT_PLAN §9 item 245; brief ``docs/planning/chunks/POST-1.md`` (§9 item
245).  The 10-03 weekly's re-pose of 7b: the 7b version (item 203, C4 pairs on
the 1 g rung, red) is at ``cbd32287`` on ``attempt/POST-1-step7b-20261003T141500Z``.

Fixture: `TH-16` step 3c's mirror-built F-small mesh — the quadrant
``_build(True, symmetry_quadrant=True)`` mirrored across y = 0 then x = 0 by
vertex index, so the D2 images exist cell for cell (144 136 cells, phantom
tag 3 with 1 256 cells).  The port sheets are the **un-narrowed** cut-plane
box faces (the whole CAD face; the 2026-10-02 ruling forbids a sub-sheet
triangulated out of a box), width A/h, κ off.  64 MHz, N1curl degree 1, one
four-drive ``run_n_port_sparameter_sweep(..., keep_fields=True)``.

Step 7a's one-sided extraction (``candidate_cells``, item 187) unchanged.
Residual r = ε_c E⁻·n − E⁺·n; N_j^{(k)} = Σ_f (A_f/3) Σ_q w_j |r|², D the same
with |E⁺·n|².  Drive k = port number − 1 sits at θ_k = kπ/2.

D2 maps (mod 4): σ_x k ↦ 2 − k, j ↦ 2 − j; σ_y k ↦ −k, j ↦ −j; σ_π k ↦ k + 2,
j ↦ j + 2.  Identity: N_{σj}^{(σk)} = N_j^{(k)}, ρ_{σf}^{(σk)} = ρ_f^{(k)}.

Asserted: (M0) the mesh is `TH-16`'s; (A1) uniform σ; (I3) ``S_f``'s D2 pairs;
(X0a) F = 2V − 4; (X0b) one tag-3, one tag-2 owned cell per point; (X1)
tangential continuity; (P) partition identity; (B) legacy route gives
|ε_c − 1|; (M1) the face set is D2-closed; (G1′) / (G2′) / (G3′) the D2 pairs
at the imported ``C4_COVARIANCE_BAND``; (W0) rank-invariant face set (``-n 4``
against the ``-n 2`` npz).  Everything else, including item 203's dropped C4
pairing on this mesh, is printed.  Scope: degree 1, 64 MHz, one fixture.
"""

from __future__ import annotations

import os
import time

import dolfinx
import numpy as np
import pytest
from dolfinx import geometry
from mpi4py import MPI

from fem_em_solver.core import HomogeneousMaterial, TimeHarmonicProblem
from fem_em_solver.ports.definitions import PortDefinition
from fem_em_solver.ports.lumped import LumpedSheetPortSpec
from fem_em_solver.ports.sparameters import run_n_port_sparameter_sweep
from fem_em_solver.post import mean_sar
from fem_em_solver.post.evaluation import evaluate_vector_field_parallel
from tests.complex_mode import complex_only
from tests.mesh.test_birdcage_port_sheets import _build
from tests.mesh.test_birdcage_port_tags import PHANTOM_HEIGHT, PHANTOM_RADIUS
from tests.mesh.test_birdcage_symmetry_quadrant import EXACT_IDENTITY_RTOL
from tests.mesh.test_two_torus_port_sheet import _sheet_extents, _sheet_facet_count
from tests.validation.test_birdcage_b1_plus_map import C4_COVARIANCE_BAND
from tests.validation.test_birdcage_sar_integral import (
    PARTITION_EPS_M,
    PARTITION_RTOL,
)
from tests.validation.test_birdcage_sar_map import (
    PHANTOM_SIGMA_RTOL,
    _phantom_sigma_range,
)
from tests.validation.test_lossy_sphere_fullwave import (
    FREQUENCY_64_HZ,
    SALINE_EPSILON_R,
    SALINE_SIGMA,
    _series,
)
from tests.validation.test_mat15_coil_driven_heating_phantom import CAD_AREA_M2
from tests.validation.test_port_birdcage_four_port import (
    LEG_D0_REPRODUCTION_BAND,
    PASSIVITY_SIGMA_TOLERANCE,
)
from tests.validation.test_port_birdcage_lumped_column import (
    CONDUCTOR_CELL_TAG,
    PHANTOM_CELL_TAG,
)
from tests.validation.test_port_gap_voltage_impedance import SIGMA_WIRE_S_PER_M
from tests.validation.test_th16_quarter_birdcage import (
    FULL_CELL_RECORD,
    FULL_SHEET_TAG,
    FULL_SHEETS,
    FULL_Z0_OHM,
    KEPT_HALF,
    MIDPOINT_MATCH_TOL_M,
    PLANE_TOL,
    QUARTER_CELL_RECORD,
    SHEET_AREA_RTOL,
    SHEET_FACET_RECORD,
    SHEET_PLANE_AXIS,
    SHEET_TAG,
    SHEET_TRANSVERSE_AXIS,
    _facet_midpoints,
    _gather_owned,
    _mirror_by_index,
    _owned_cell_tags,
    _owned_count,
    _owned_exterior_on_plane,
    _sheet_integral_area,
    _tag_counts,
    _tri_area_owned,
    _volume,
)
from tests.validation.test_th16_symmetry_plane_driven import _distribute

# Air is tag 2 (test_port_birdcage_lumped_column.py:96-99 comment).
AIR_CELL_TAG = 2

ON_SURFACE_TOL = 1e-9  # relative; never tuned (brief trap)
TANGENTIAL_TOL = 1e-10
LEGACY_TOL = 1e-12
# Reporting rule (brief, not a band): a D2 reading above this is journaled.
UNEXPLAINED_LEVEL = 1e-8
NPZ_PATH = "/workspace/logs/post1/step7b2-n2.npz"
_FACES = ((1, 2, 3), (0, 2, 3), (0, 1, 3), (0, 1, 2))
# Barycentric (2/3, 1/6, 1/6) and its permutations: rows are points.
_BARY = np.array(
    [[2 / 3, 1 / 6, 1 / 6], [1 / 6, 2 / 3, 1 / 6], [1 / 6, 1 / 6, 2 / 3]],
    dtype=np.float64,
)

# D2 maps on drive / quadrant indices and on points.
SIGMA_MAPS = {
    "sx": (lambda i: (2 - i) % 4, np.array([-1.0, 1.0, 1.0])),
    "sy": (lambda i: (-i) % 4, np.array([1.0, -1.0, 1.0])),
    "spi": (lambda i: (i + 2) % 4, np.array([-1.0, -1.0, 1.0])),
}
# Sheet k's mean midpoint side (P1 x > 0, P2 y > 0, P3 x < 0, P4 y < 0).
SHEET_SIDE = {1: (0, 1.0), 2: (1, 1.0), 3: (0, -1.0), 4: (1, -1.0)}


def _w(points, j, eps=PARTITION_EPS_M):
    """Numpy twin of ``test_birdcage_sar_integral._quadrant_weight``."""
    theta = j * np.pi / 2.0
    x, y = points[:, 0], points[:, 1]
    radius = np.sqrt(x * x + y * y + eps * eps)
    c = (x * float(np.cos(theta)) + y * float(np.sin(theta))) / radius
    return ((c + np.sqrt(c * c + eps * eps)) / 2.0) ** 2


def _interface_faces(msh, cell_tags, comm):
    """Allgathered, sorted, deduplicated interface faces (identical on ranks)."""
    R, h2 = PHANTOM_RADIUS, PHANTOM_HEIGHT / 2.0
    n_owned = msh.topology.index_map(3).size_local
    cells = cell_tags.find(PHANTOM_CELL_TAG)
    cells = cells[cells < n_owned]
    dofmap = msh.geometry.dofmap
    x = msh.geometry.x[:, :3]
    input_ids = np.asarray(msh.geometry.input_global_indices, dtype=np.int64)

    verts_out, normals, vids = [], [], []
    for c in cells:
        nodes = np.asarray(dofmap[c])[:4]
        xc = x[nodes]
        centroid = xc.mean(axis=0)
        for face in _FACES:
            fn = nodes[list(face)]
            v = x[fn]
            rho = np.linalg.norm(v[:, :2], axis=1)
            lateral = bool(
                np.all(np.abs(rho - R) <= ON_SURFACE_TOL * R)
                and np.all(np.abs(v[:, 2]) <= h2 * (1.0 + ON_SURFACE_TOL))
            )
            cap = bool(
                np.all(np.abs(np.abs(v[:, 2]) - h2) <= ON_SURFACE_TOL * h2)
                and np.all(rho <= R * (1.0 + ON_SURFACE_TOL))
            )
            if not (lateral or cap):
                continue
            order = np.lexsort((v[:, 2], v[:, 1], v[:, 0]))
            v = v[order]
            p = v.mean(axis=0)
            n = np.cross(v[1] - v[0], v[2] - v[0])
            n = n / np.linalg.norm(n)
            if float(np.dot(n, p - centroid)) < 0.0:
                n = -n
            verts_out.append(v.reshape(-1))
            normals.append(n)
            vids.append([int(input_ids[k]) for k in fn[order]])

    def _cat(lst, shape, dtype=np.float64):
        return np.asarray(lst, dtype=dtype).reshape(shape)

    g_verts = np.concatenate([_cat(v, (-1, 9)) for v in comm.allgather(verts_out)])
    g_normals = np.concatenate([_cat(v, (-1, 3)) for v in comm.allgather(normals)])
    g_vids = np.concatenate(
        [_cat(v, (-1, 3), np.int64) for v in comm.allgather(vids)]
    )
    order = np.lexsort(tuple(g_verts[:, j] for j in range(8, -1, -1)))
    g_verts, g_normals, g_vids = g_verts[order], g_normals[order], g_vids[order]

    same = np.all(g_verts[1:] == g_verts[:-1], axis=1)
    dup = np.zeros(g_verts.shape[0], dtype=bool)
    dup[1:] |= same
    dup[:-1] |= same
    keep = ~dup
    return (
        np.ascontiguousarray(g_verts[keep]),
        np.ascontiguousarray(g_normals[keep]),
        g_vids[keep],
        int(np.sum(dup)),
    )


def _side_counts(msh, cell_tags, points, comm):
    """Global count of owned containing cells of tag 3, tag 2 and in total."""
    tdim = msh.topology.dim
    imap = msh.topology.index_map(tdim)
    n_owned = imap.size_local
    tag_of = np.zeros(n_owned + imap.num_ghosts, dtype=np.int32)
    tag_of[cell_tags.indices] = cell_tags.values
    bb = geometry.bb_tree(msh, tdim)
    hits = geometry.compute_collisions_points(bb, points)
    coll = geometry.compute_colliding_cells(msh, hits, points)
    local = np.zeros((points.shape[0], 3), dtype=np.int64)
    for i in range(points.shape[0]):
        links = np.asarray(coll.links(i))
        owned = links[links < n_owned]
        local[i, 0] = int(np.sum(tag_of[owned] == PHANTOM_CELL_TAG))
        local[i, 1] = int(np.sum(tag_of[owned] == AIR_CELL_TAG))
        local[i, 2] = int(owned.size)
    total = np.empty_like(local)
    comm.Allreduce(local, total, op=MPI.SUM)
    return total


def _face_image_map(centroids, normals, areas, mirror):
    """(M1): for each face the index of the face at its mirror image, or -1."""
    img = centroids * mirror[None, :]
    d = np.linalg.norm(img[:, None, :] - centroids[None, :, :], axis=2)
    hits = d <= MIDPOINT_MATCH_TOL_M
    n_hit = hits.sum(axis=1)
    perm = np.where(n_hit == 1, np.argmax(hits, axis=1), -1)
    ok = perm >= 0
    area_rel = np.full(perm.size, np.inf)
    area_rel[ok] = np.abs(areas[perm[ok]] / areas[ok] - 1.0)
    normal_abs = np.full(perm.size, np.inf)
    normal_abs[ok] = np.max(np.abs(normals[perm[ok]] - normals[ok] * mirror[None, :]), axis=1)
    return {
        "perm": perm,
        "n_hit_min": int(n_hit.min()),
        "n_hit_max": int(n_hit.max()),
        "worst_dist": float(np.max(np.where(hits, d, 0.0))),
        "area_rel": float(area_rel.max()),
        "normal_abs": float(normal_abs.max()),
    }


def _quarter(comm):
    """`TH-16`'s quarter with un-narrowed sheets: records and the rank-0 gather."""
    out = {}
    t0 = time.perf_counter()
    msh, cells, _facets, diag, t_mesh = _build(True, symmetry_quadrant=True)
    tdim = msh.topology.dim
    fdim = tdim - 1
    # Hoisted on every rank before any facet form (known-issues 9).
    msh.topology.create_connectivity(fdim, tdim)
    msh.topology.create_entity_permutations()
    n_owned_f = msh.topology.index_map(fdim).size_local
    out["cells"] = int(msh.topology.index_map(tdim).size_global)

    # `test_th16_quarter_birdcage.py:304-317`, then `:328-341`; the narrowing
    # (`:318-327`) is omitted: each sheet is the whole cut-plane box face.
    sheet_idx, sheet_val, unnarrowed = [], [], {}
    for p in (1, 2):
        f = _owned_exterior_on_plane(msh, cells, SHEET_PLANE_AXIS[p], KEPT_HALF[p])
        unnarrowed[p] = comm.allreduce(int(f.size), op=MPI.SUM)
        sheet_idx.append(f)
        sheet_val.append(np.full(f.size, SHEET_TAG[p], dtype=np.int32))
    idx = np.concatenate(sheet_idx).astype(np.int32)
    val = np.concatenate(sheet_val).astype(np.int32)
    order = np.argsort(idx, kind="stable")
    idx, val = idx[order], val[order]
    assert np.unique(idx).size == idx.size
    tags_s = dolfinx.mesh.meshtags(msh, fdim, idx, val)
    out["sheets_unnarrowed"] = unnarrowed

    sheets = {}
    for p in (1, 2):
        tag = SHEET_TAG[p]
        n = _sheet_facet_count(msh, tags_s, tag, comm)
        f_owned = np.asarray(tags_s.find(tag))
        f_owned = f_owned[f_owned < n_owned_f]
        area = comm.allreduce(_tri_area_owned(msh, f_owned), op=MPI.SUM)
        ext = _sheet_extents(msh, tags_s, tag, comm)
        h = float(ext[2])
        mids = np.vstack(comm.allgather(_facet_midpoints(msh, f_owned)))
        sheets[p] = {"facets": int(n), "area": area, "h": h, "w": area / h, "mid": mids,
                     "ds_area": _sheet_integral_area(msh, tags_s, tag, False),
                     "w_bbox": float(ext[SHEET_TRANSVERSE_AXIS[p]]),
                     "flat": float(ext[SHEET_PLANE_AXIS[p]])}
    out["sheets"] = sheets
    out["tag_counts"] = _tag_counts(comm, _owned_cell_tags(msh, cells))
    out["volume"] = _volume(msh)
    out["gathered"] = _gather_owned(msh, cells)
    out["t_mesh"] = t_mesh
    out["mesh_wall"] = float(diag["mesh_wall_time_s"])
    out["t_quarter"] = time.perf_counter() - t0
    return out


def _mirror(comm, q):
    """`test_th16_quarter_birdcage.py:513-598` with the un-narrowed sheets."""
    out = {}
    comm.Barrier()
    t0 = time.perf_counter()
    x = c = t = None
    if comm.rank == 0:
        x, c, t = q["gathered"]
        for axis in (1, 0):
            x, c, t = _mirror_by_index(x, c, t, axis)
    out["builder_vertices"] = comm.bcast(None if x is None else int(x.shape[0]), root=0)
    msh = _distribute(x, c, comm)
    tdim = msh.topology.dim
    fdim = tdim - 1
    if comm.rank == 0:
        ents = np.ascontiguousarray(c, dtype=np.int64)
        vals = np.ascontiguousarray(t, dtype=np.int32)
    else:
        ents = np.empty((0, 4), dtype=np.int64)
        vals = np.empty(0, dtype=np.int32)
    local_entities, local_values = dolfinx.io.distribute_entity_data(msh, tdim, ents, vals)
    msh.topology.create_connectivity(tdim, tdim)
    cells = dolfinx.mesh.meshtags_from_entities(
        msh, tdim, dolfinx.graph.adjacencylist(local_entities),
        local_values.astype(np.int32, copy=False),
    )
    del x, c, t, ents, vals
    # Hoisted on every rank before any facet form (known-issues 9).
    msh.topology.create_connectivity(fdim, tdim)
    msh.topology.create_entity_permutations()
    comm.Barrier()
    out["t_build"] = time.perf_counter() - t0

    n_owned_f = msh.topology.index_map(fdim).size_local
    out["cells"] = int(msh.topology.index_map(tdim).size_global)
    out["vertices"] = int(msh.geometry.index_map().size_global)
    out["tag_counts"] = _tag_counts(comm, _owned_cell_tags(msh, cells))
    out["volume"] = _volume(msh)

    ext = np.asarray(dolfinx.mesh.exterior_facet_indices(msh.topology), dtype=np.int32)
    ext = ext[ext < n_owned_f]
    out["exterior_on_plane"] = {}
    for axis in (0, 1):
        n_on = 0
        if ext.size:
            nodes = dolfinx.cpp.mesh.entities_to_geometry(msh._cpp_object, fdim, ext, False)
            n_on = int(np.count_nonzero(
                np.all(np.abs(msh.geometry.x[nodes][:, :, axis]) < PLANE_TOL, axis=1)
            ))
        out["exterior_on_plane"][axis] = comm.allreduce(n_on, op=MPI.SUM)

    plane = {}
    for axis in (0, 1):
        f = dolfinx.mesh.locate_entities(
            msh, fdim, lambda p, a=axis: np.abs(p[a]) < PLANE_TOL
        ).astype(np.int32)
        f = f[f < n_owned_f]
        plane[axis] = (f, _facet_midpoints(msh, f))
    idx, val, worst = [], [], 0.0
    side = {}
    for k, (p, flip) in FULL_SHEETS.items():
        target = q["sheets"][p]["mid"].copy()
        if flip is not None:
            target[:, flip] = -target[:, flip]
        f, mid = plane[SHEET_PLANE_AXIS[p]]
        loc_sum, loc_n = np.zeros(3), 0
        if f.size:
            d = np.linalg.norm(mid[:, None, :] - target[None, :, :], axis=2).min(axis=1)
            hit = d <= MIDPOINT_MATCH_TOL_M
            if hit.any():
                worst = max(worst, float(d[hit].max()))
                loc_sum = mid[hit].sum(axis=0)
                loc_n = int(np.count_nonzero(hit))
            idx.append(f[hit])
            val.append(np.full(int(np.count_nonzero(hit)), FULL_SHEET_TAG[k], dtype=np.int32))
        tot = comm.allreduce(loc_sum, op=MPI.SUM)
        n_tot = comm.allreduce(loc_n, op=MPI.SUM)
        side[k] = tot / max(n_tot, 1)
    idx = np.concatenate(idx).astype(np.int32) if idx else np.empty(0, dtype=np.int32)
    val = np.concatenate(val).astype(np.int32) if val else np.empty(0, dtype=np.int32)
    assert np.unique(idx).size == idx.size
    order = np.argsort(idx, kind="stable")
    tags_s = dolfinx.mesh.meshtags(msh, fdim, idx[order], val[order])
    out["midpoint_worst"] = comm.allreduce(worst, op=MPI.MAX)
    out["sheet_mean_mid"] = side
    out["sheet_facets"] = {
        k: _owned_count(comm, tags_s, FULL_SHEET_TAG[k], n_owned_f) for k in FULL_SHEETS
    }
    out["sheet_area"] = {
        k: _sheet_integral_area(msh, tags_s, FULL_SHEET_TAG[k], True) for k in FULL_SHEETS
    }
    out["mesh"], out["cell_tags"], out["tags_s"] = msh, cells, tags_s
    return out


@pytest.fixture(scope="module")
def step7b2():
    comm = MPI.COMM_WORLD
    t_window = time.perf_counter()
    eps_c = _series(FREQUENCY_64_HZ).epsilon_c

    q = _quarter(comm)
    m = _mirror(comm, q)
    msh, cell_tags, tags_s = m["mesh"], m["cell_tags"], m["tags_s"]

    # Problem as `test_th16_quarter_birdcage.py:603-641` at 64 MHz.
    problem = TimeHarmonicProblem(
        mesh=msh,
        frequency_hz=FREQUENCY_64_HZ,
        material=HomogeneousMaterial(sigma=0.0, epsilon_r=1.0, mu_r=1.0),
        cell_tags=cell_tags,
        material_map={
            CONDUCTOR_CELL_TAG: HomogeneousMaterial(
                sigma=SIGMA_WIRE_S_PER_M, epsilon_r=1.0, mu_r=1.0
            ),
            PHANTOM_CELL_TAG: HomogeneousMaterial(
                sigma=SALINE_SIGMA, epsilon_r=SALINE_EPSILON_R, mu_r=1.0
            ),
        },
        boundary_condition="pec_zero_tangential_a",
    )
    defs, specs = [], []
    for k, (p, _flip) in FULL_SHEETS.items():
        pid = f"P{k}"
        defs.append(
            PortDefinition(
                port_id=pid,
                positive_tag=FULL_SHEET_TAG[k],
                negative_tag=CONDUCTOR_CELL_TAG,
                orientation="leg_gap_axial_plus_z",
                z0_ohm=FULL_Z0_OHM,
            )
        )
        specs.append(
            LumpedSheetPortSpec(
                port_id=pid,
                facet_tag=FULL_SHEET_TAG[k],
                port_impedance_ohm=FULL_Z0_OHM,
                gap_height_m=q["sheets"][p]["h"],
                sheet_width_m=q["sheets"][p]["w"],
                drive_direction=(0.0, 0.0, 1.0),
                drive_voltage_v=1.0 + 0.0j,
                interior=True,
            )
        )
    comm.Barrier()
    t0 = time.perf_counter()
    result = run_n_port_sparameter_sweep(
        problem, defs, lumped_sheet_ports=specs, lumped_sheet_facet_tags=tags_s,
        degree=1, keep_fields=True,
    )
    comm.Barrier()
    t_sweep = time.perf_counter() - t0
    assert tuple(result.port_ids) == tuple(d.port_id for d in defs)
    sf = np.asarray(result.s_matrix, dtype=np.complex128).copy()
    fields = {k: result.fields[f"P{k + 1}"] for k in range(4)}
    off = ~np.eye(4, dtype=bool)
    recip = float(np.max(np.abs(sf - sf.T)[off] / np.abs(sf.T)[off]))
    adj = np.array([sf[i, j] for i in range(4) for j in range(4) if (i - j) % 2 == 1])
    d2s = {
        "S11/S33": float(abs(sf[0, 0] - sf[2, 2]) / abs(0.5 * (sf[0, 0] + sf[2, 2]))),
        "S22/S44": float(abs(sf[1, 1] - sf[3, 3]) / abs(0.5 * (sf[1, 1] + sf[3, 3]))),
        "adjacent": float(np.max(np.abs(adj - adj.mean())) / abs(adj.mean())),
    }

    # (A1).
    n_owned = msh.topology.index_map(3).size_local
    tag3 = cell_tags.find(PHANTOM_CELL_TAG)
    tag3_global = comm.allreduce(int(np.sum(tag3 < n_owned)), op=MPI.SUM)
    sigma_lo, sigma_hi, sigma_count = _phantom_sigma_range(
        {"mesh": msh, "cell_tags": cell_tags}, fields[0].sigma_field
    )

    # Extraction (item 203's loop, unchanged).
    t0 = time.perf_counter()
    verts, normals, vids, n_dropped = _interface_faces(msh, cell_tags, comm)
    F = int(verts.shape[0])
    V = int(np.unique(vids).size)
    v3 = verts.reshape(-1, 3, 3)
    areas = 0.5 * np.linalg.norm(np.cross(v3[:, 1] - v3[:, 0], v3[:, 2] - v3[:, 0]), axis=1)
    centroids = v3.mean(axis=1)
    pts = np.ascontiguousarray(np.einsum("qi,fij->fqj", _BARY, v3).reshape(-1, 3))
    nq = np.repeat(normals, 3, axis=0)
    wq = np.repeat(areas / 3.0, 3)
    w_j = np.stack([_w(pts, j) for j in range(4)])  # (4, 3F)
    counts = _side_counts(msh, cell_tags, pts, comm)

    cand_in = cell_tags.find(PHANTOM_CELL_TAG)
    cand_out = cell_tags.find(AIR_CELL_TAG)
    drives = {}
    for k in range(4):
        e = fields[k].e_complex
        e_in, v_in = evaluate_vector_field_parallel(
            e, pts, comm, owner="lowest_original_cell", candidate_cells=cand_in
        )
        e_out, v_out = evaluate_vector_field_parallel(
            e, pts, comm, owner="lowest_original_cell", candidate_cells=cand_out
        )
        e_leg, v_leg = evaluate_vector_field_parallel(e, pts, comm)
        e_in, e_out, e_leg = (np.asarray(z, dtype=complex) for z in (e_in, e_out, e_leg))
        en_in = np.sum(e_in * nq, axis=1)
        en_out = np.sum(e_out * nq, axis=1)
        en_leg = np.sum(e_leg * nq, axis=1)
        r2 = np.abs(eps_c * en_in - en_out) ** 2
        d2 = np.abs(en_out) ** 2
        N = float(np.sum(wq * r2))
        D = float(np.sum(wq * d2))
        Nj = [float(np.sum(wq * w_j[j] * r2)) for j in range(4)]
        Dj = [float(np.sum(wq * w_j[j] * d2)) for j in range(4)]

        def _rn(eps):
            return float(np.sqrt(np.sum(wq * np.abs(eps * en_in - en_out) ** 2) / D))

        r_leg = float(
            np.sqrt(
                np.sum(wq * np.abs(eps_c * en_leg - en_leg) ** 2)
                / np.sum(wq * np.abs(en_leg) ** 2)
            )
        )
        legacy_tag3 = float(
            np.mean(
                np.linalg.norm(e_leg - e_in, axis=1)
                < np.linalg.norm(e_leg - e_out, axis=1)
            )
        )
        x1 = float(np.max(np.linalg.norm(np.cross(nq, e_in - e_out), axis=1)))
        x1_scale = float(np.max(np.linalg.norm(e_out, axis=1)))
        power = float(
            mean_sar(
                e,
                sigma=fields[k].sigma_field,
                rho=1000.0,
                subdomain_ids=PHANTOM_CELL_TAG,
                cell_tags=cell_tags,
                comm=comm,
            )["dissipated_power_w"]
        )
        drives[k] = {
            "N": N, "D": D, "Nj": Nj, "Dj": Dj,
            "r_n": float(np.sqrt(N / D)),
            "r_n_blind": _rn(eps_c.real),
            "r_n_conj": _rn(np.conj(eps_c)),
            "r_n_legacy": r_leg,
            "legacy_tag3": legacy_tag3,
            "valid": (np.asarray(v_in), np.asarray(v_out), np.asarray(v_leg)),
            "x1": x1, "x1_scale": x1_scale,
            "rho_face": r2.reshape(-1, 3).mean(axis=1),
            "power": power,
        }
    t_extract = time.perf_counter() - t0
    del result
    t_window = time.perf_counter() - t_window

    # (M1) face image maps, identical on every rank (allgathered arrays).
    m1 = {s: _face_image_map(centroids, normals, areas, SIGMA_MAPS[s][1]) for s in ("sx", "sy")}

    out = {
        "comm": comm, "eps_c": eps_c, "q": q, "m": m,
        "S": sf, "recip": recip, "d2s": d2s,
        "sigma_max": float(np.linalg.svd(sf, compute_uv=False).max()),
        "tag3_global": tag3_global,
        "sigma": (sigma_lo, sigma_hi, sigma_count),
        "verts": verts, "F": F, "V": V, "n_dropped": n_dropped,
        "n_points": int(pts.shape[0]), "area": float(np.sum(areas)),
        "counts": counts, "drives": drives, "m1": m1,
        "t_sweep": t_sweep, "t_extract": t_extract, "t_window": t_window,
    }
    _report(out)

    if comm.size == 2:
        if comm.rank == 0:
            os.makedirs(os.path.dirname(NPZ_PATH), exist_ok=True)
            np.savez(
                NPZ_PATH,
                verts=verts, F=F, V=V, area=out["area"],
                N=np.array([drives[k]["N"] for k in range(4)]),
                D=np.array([drives[k]["D"] for k in range(4)]),
                r_n=np.array([drives[k]["r_n"] for k in range(4)]),
            )
        comm.Barrier()
    return out


# ---- D2 readings --------------------------------------------------------------


def _g1p(d, s, k):
    kk = SIGMA_MAPS[s][0](k)
    return abs(d[kk]["N"] - d[k]["N"]) / d[k]["N"]


def _g2p_pairs():
    pairs = []
    for s in ("sx", "sy", "spi"):
        f = SIGMA_MAPS[s][0]
        for k in range(4):
            for j in range(4):
                if (f(k), f(j)) != (k, j):
                    pairs.append((s, k, j))
    return pairs


G2P_PAIRS = _g2p_pairs()
G1P_PAIRS = [("sx", 0), ("sx", 2), ("sy", 1), ("sy", 3)]
G3P_PAIRS = [(s, k) for s in ("sx", "sy") for k in range(4)]


def _g2p(d, s, k, j, key="Nj"):
    f = SIGMA_MAPS[s][0]
    ref = d[k][key][j]
    return abs(d[f(k)][key][f(j)] - ref) / ref


def _g3p(o, s, k):
    perm = o["m1"][s]["perm"]
    kk = SIGMA_MAPS[s][0](k)
    rho_k = o["drives"][k]["rho_face"]
    rho_img = o["drives"][kk]["rho_face"][perm]
    return float(np.max(np.abs(rho_img - rho_k) / rho_k))


def _report(o):
    if o["comm"].rank != 0:
        return
    eps_c = o["eps_c"]
    em1 = abs(eps_c - 1.0)
    d = o["drives"]
    q, m = o["q"], o["m"]
    tag = "[POST-1 7b']"
    print(
        f"\n{tag} -n {o['comm'].size}, TH-16 mirror-built F-small birdcage, un-narrowed "
        f"sheets, f {FREQUENCY_64_HZ:.6e} Hz, degree 1; eps_c = {eps_c!r}, |eps_c - 1| = {em1!r}"
    )
    print(
        f"  (M0) quarter cells {q['cells']} (record {QUARTER_CELL_RECORD}); mesh "
        f"{q['mesh_wall']:.2f} s, rung {q['t_mesh']:.2f} s, quarter total {q['t_quarter']:.2f} s; "
        f"un-narrowed sheets {q['sheets_unnarrowed']} (record {SHEET_FACET_RECORD})"
    )
    for p in (1, 2):
        s = q["sheets"][p]
        print(
            f"  P{p}q whole face: {s['facets']} facets area={s['area']:.9e} m^2 "
            f"h={s['h']:.9e} m w=A/h={s['w']:.9e} m w_bbox={s['w_bbox']:.9e} m "
            f"flat={s['flat']:.3e} m ds/geom-1={s['ds_area'] / s['area'] - 1.0:.2e}"
        )
    qc, fc = q["tag_counts"], m["tag_counts"]
    print(
        f"  (M0) full cells {m['cells']} (record {FULL_CELL_RECORD}); vertices {m['vertices']} "
        f"(builder {m['builder_vertices']}); mirror build {m['t_build']:.2f} s; per-tag "
        + " ".join(f"{k}:{fc.get(k)}/{qc[k]}" for k in sorted(qc))
        + f"; |V_f/(4V_q)-1| {abs(m['volume'] / (4.0 * q['volume']) - 1.0):.3e}; exterior "
        f"on planes {m['exterior_on_plane']}"
    )
    print(
        "  (M0) full sheets owned facets "
        + " / ".join(str(m["sheet_facets"][k]) for k in FULL_SHEETS)
        + f" (50 / 56 / 50 / 56); worst midpoint {m['midpoint_worst']:.3e} m; dS area rel "
        + " ".join(
            f"{abs(m['sheet_area'][k] / q['sheets'][p]['area'] - 1.0):.2e}"
            for k, (p, _f) in FULL_SHEETS.items()
        )
        + "; mean midpoints "
        + " ".join(
            f"P{k}({v[0]:+.4e},{v[1]:+.4e},{v[2]:+.4e})" for k, v in m["sheet_mean_mid"].items()
        )
    )
    print(
        f"  sweep {o['t_sweep']:.2f} s; S_f=\n" + np.array2string(o["S"], precision=9)
        + f"\n  (I3) D2 pairs " + " ".join(f"{k} {v:.3e}" for k, v in o["d2s"].items())
        + f" (band {LEG_D0_REPRODUCTION_BAND:.0e}, pred <= 1e-11); (I1) reciprocity "
        f"{o['recip']:.3e}; (I2) sigma_max {o['sigma_max']:.9f}"
    )
    print(f"  (A1) sigma tag-3 [{o['sigma'][0]!r}, {o['sigma'][1]!r}] over {o['sigma'][2]} "
          f"cells; tag-3 global {o['tag3_global']}")
    c = o["counts"]
    one_sided = int(np.sum((c[:, 0] == 0) ^ (c[:, 1] == 0)))
    print(
        f"  F {o['F']} (pred [400,900], mod 4 = {o['F'] % 4}), V {o['V']}, 2V-4 {2 * o['V'] - 4}, "
        f"dropped (chords) {o['n_dropped']}, points {o['n_points']}; side counts tag3 "
        f"[{int(c[:, 0].min())},{int(c[:, 0].max())}] tag2 [{int(c[:, 1].min())},"
        f"{int(c[:, 1].max())}] total [{int(c[:, 2].min())},{int(c[:, 2].max())}], "
        f"one-sided points {one_sided}"
    )
    print(
        f"  sum A_f {o['area']!r} vs CAD {CAD_AREA_M2} (rel "
        f"{(o['area'] - CAD_AREA_M2) / CAD_AREA_M2:+.4e})"
    )
    for s, r in o["m1"].items():
        perm = r["perm"]
        inv = bool(np.all(perm >= 0) and np.array_equal(perm[perm], np.arange(perm.size)))
        fixed = int(np.sum(perm == np.arange(perm.size)))
        print(
            f"  (M1) {s}: hits/face [{r['n_hit_min']},{r['n_hit_max']}], worst dist "
            f"{r['worst_dist']:.3e} m, involution {inv}, fixed points {fixed}, area rel "
            f"{r['area_rel']:.3e}, normal abs {r['normal_abs']:.3e}"
        )
    for k in range(4):
        r = d[k]
        print(
            f"  drive k={k} (P{k + 1}): R_n {r['r_n']!r} (pred [0.8,1.5]); N {r['N']!r} "
            f"D {r['D']!r}; partition N {abs(sum(r['Nj']) - r['N']) / r['N']:.3e} "
            f"D {abs(sum(r['Dj']) - r['D']) / r['D']:.3e}"
        )
        print(
            f"    X1 {r['x1']:.3e} vs max|E+| {r['x1_scale']:.4e} (ratio "
            f"{r['x1'] / r['x1_scale']:.3e}); legacy R_n {r['r_n_legacy']!r} (rel "
            f"{abs(r['r_n_legacy'] - em1) / em1:.3e}); (B) separation {em1 / r['r_n']:.1f}x "
            f"(pred 100-200, ceiling 160.155/R_n)"
        )
        print(
            f"    sigma-blind {r['r_n_blind']:.4f} (pred [0.95,1.10]), eps_c* "
            f"{r['r_n_conj']:.4f} (pred [1.6,2.6]), legacy tag-3 fraction "
            f"{r['legacy_tag3']:.3f} (~0.5); phantom power {r['power']!r} W"
        )
    for s in ("sx", "sy", "spi"):
        f = SIGMA_MAPS[s][0]
        print(
            f"  phantom power D2 pairs {s}: "
            + " ".join(
                f"k{k}->{f(k)} {abs(d[f(k)]['power'] - d[k]['power']) / d[k]['power']:.3e}"
                for k in range(4) if f(k) != k
            )
        )
    for key, name in (("Nj", "N_j^(k)"), ("Dj", "D_j^(k)")):
        print(f"  {name} table (rows k, cols j):")
        for k in range(4):
            print("    " + "  ".join(f"{v:.9e}" for v in d[k][key]))
    print("  quadrant R_n table:")
    for k in range(4):
        print("    " + "  ".join(
            f"{np.sqrt(n / dd):.5f}" for n, dd in zip(d[k]["Nj"], d[k]["Dj"])
        ))
    g1 = {pk: _g1p(d, *pk) for pk in G1P_PAIRS}
    print(
        "  (G1') whole-interface D2 pairs: "
        + ", ".join(f"{s} k={k} {v:.3e}" for (s, k), v in g1.items())
        + f" (band {C4_COVARIANCE_BAND:.0e}, pred <= 1e-9)"
    )
    g2 = {pk: _g2p(d, *pk) for pk in G2P_PAIRS}
    g2d = {pk: _g2p(d, *pk, key="Dj") for pk in G2P_PAIRS}
    print(f"  (G2') quadrant D2 pairs, N then D (40 readings; band {C4_COVARIANCE_BAND:.0e}, "
          "pred <= 1e-9):")
    for pk in G2P_PAIRS:
        print(f"    {pk[0]} k={pk[1]} j={pk[2]}: N {g2[pk]:.3e}  D {g2d[pk]:.3e}")
    g3 = {pk: _g3p(o, *pk) for pk in G3P_PAIRS}
    print(
        "  (G3') face-level D2 max rel: "
        + ", ".join(f"{s} k={k} {v:.3e}" for (s, k), v in g3.items())
        + f" (band {C4_COVARIANCE_BAND:.0e}, pred <= 1e-7)"
    )
    for k in range(4):
        rho = d[k]["rho_face"]
        print(f"    k={k}: min rho_f / mean rho {rho.min() / rho.mean():.3e}")
    worst = max(max(g1.values()), max(g2.values()), max(g3.values()))
    print(
        f"  D2 worst reading {worst:.3e}; reporting rule (> {UNEXPLAINED_LEVEL:.0e} is "
        f"journaled as unexplained): {'ABOVE' if worst > UNEXPLAINED_LEVEL else 'below'}"
    )
    # (N1) half-mapped pairing on drives sigma moves, printed only.
    g2max = max(g2.values())
    for s, ks, js in (("sx", (0, 2), (0, 2)), ("sy", (1, 3), (1, 3))):
        f = SIGMA_MAPS[s][0]
        for k in ks:
            nk = d[k]["Nj"]
            ceiling = max(nk) / min(nk) - 1.0
            vals = [abs(nk[f(j)] - nk[j]) / nk[j] for j in js]
            print(
                f"  (N1) {s} k={k}: half-mapped " + " ".join(f"j={j} {v * 100:.3f}%"
                                                             for j, v in zip(js, vals))
                + f" (pred [50,400]%), ceiling {ceiling * 100:.3f}%, separation over G2' max "
                f"{min(vals) / g2max if g2max > 0 else float('inf'):.3e}x (pred >= 1e6)"
            )
    # (N2) face-level half-mapped pairing.
    for s, ks in (("sx", (0, 2)), ("sy", (1, 3))):
        perm = o["m1"][s]["perm"]
        for k in ks:
            rho = d[k]["rho_face"]
            rel = np.abs(rho[perm] - rho) / rho
            print(
                f"  (N2) {s} k={k}: face half-mapped median {np.median(rel) * 100:.2f}%, "
                f"p90 {np.percentile(rel, 90) * 100:.2f}% (pred median >= 10%)"
            )
    # Item 203's dropped C4 statement, as a record only (this mesh's C4 contrast).
    c4g1 = [abs(d[(k + 1) % 4]["N"] - d[k]["N"]) / d[k]["N"] for k in range(4)]
    print(
        "  dropped C4 record (never gated): G1 "
        + ", ".join(f"{v * 100:.4f}%" for v in c4g1)
        + " (pred 0-20%); G2 (rows k, cols j), N then D:"
    )
    for k in range(4):
        n_row = [abs(d[(k + 1) % 4]["Nj"][(j + 1) % 4] - d[k]["Nj"][j]) / d[k]["Nj"][j]
                 for j in range(4)]
        d_row = [abs(d[(k + 1) % 4]["Dj"][(j + 1) % 4] - d[k]["Dj"][j]) / d[k]["Dj"][j]
                 for j in range(4)]
        print("    N: " + "  ".join(f"{v * 100:7.4f}%" for v in n_row)
              + "   D: " + "  ".join(f"{v * 100:7.4f}%" for v in d_row))
    print(
        f"  elapsed: quarter build {q['t_quarter']:.1f}s (mesh {q['mesh_wall']:.1f}s), mirror "
        f"build {m['t_build']:.1f}s, sweep {o['t_sweep']:.1f}s, extraction "
        f"{o['t_extract']:.1f}s, fixture {o['t_window']:.1f}s",
        flush=True,
    )


# ---- Tests ----------------------------------------------------------------------


@complex_only
def test_m0_mesh_is_th16s(step7b2):
    q, m = step7b2["q"], step7b2["m"]
    assert q["cells"] == QUARTER_CELL_RECORD
    assert m["cells"] == FULL_CELL_RECORD
    assert m["vertices"] == m["builder_vertices"]
    assert m["tag_counts"] == {k: 4 * v for k, v in q["tag_counts"].items()}
    assert abs(m["volume"] / (4.0 * q["volume"]) - 1.0) <= EXACT_IDENTITY_RTOL
    assert m["exterior_on_plane"] == {0: 0, 1: 0}
    assert q["sheets_unnarrowed"] == SHEET_FACET_RECORD
    assert {p: q["sheets"][p]["facets"] for p in (1, 2)} == SHEET_FACET_RECORD
    for k, (p, _flip) in FULL_SHEETS.items():
        assert m["sheet_facets"][k] == SHEET_FACET_RECORD[p], k
        assert abs(m["sheet_area"][k] / q["sheets"][p]["area"] - 1.0) <= SHEET_AREA_RTOL, k
        axis, sign = SHEET_SIDE[k]
        assert sign * m["sheet_mean_mid"][k][axis] > 0.0, (k, m["sheet_mean_mid"][k])


@complex_only
def test_a1_phantom_sigma_uniform(step7b2):
    lo, hi, count = step7b2["sigma"]
    assert count == step7b2["tag3_global"]
    assert count == step7b2["m"]["tag_counts"][PHANTOM_CELL_TAG]
    assert abs(lo - SALINE_SIGMA) <= PHANTOM_SIGMA_RTOL * SALINE_SIGMA, lo
    assert abs(hi - SALINE_SIGMA) <= PHANTOM_SIGMA_RTOL * SALINE_SIGMA, hi


@complex_only
def test_i3_drives_are_d2_images(step7b2):
    for k, v in step7b2["d2s"].items():
        assert v <= LEG_D0_REPRODUCTION_BAND, (k, v)


@complex_only
def test_x0a_closed_genus0_surface(step7b2):
    assert step7b2["F"] == 2 * step7b2["V"] - 4, (step7b2["F"], step7b2["V"])
    assert step7b2["n_points"] == 3 * step7b2["F"]


@complex_only
def test_x0b_one_cell_each_side(step7b2):
    c = step7b2["counts"]
    assert np.all(c[:, 0] == 1), (int(c[:, 0].min()), int(c[:, 0].max()))
    assert np.all(c[:, 1] == 1), (int(c[:, 1].min()), int(c[:, 1].max()))
    assert np.all(c[:, 2] == 2), (int(c[:, 2].min()), int(c[:, 2].max()))
    for k in range(4):
        v_in, v_out, _ = step7b2["drives"][k]["valid"]
        assert bool(np.all(v_in)) and bool(np.all(v_out)), k


@complex_only
@pytest.mark.parametrize("k", range(4))
def test_x1_tangential_continuity(step7b2, k):
    r = step7b2["drives"][k]
    assert r["x1"] <= TANGENTIAL_TOL * r["x1_scale"], (r["x1"], r["x1_scale"])


@complex_only
@pytest.mark.parametrize("k", range(4))
def test_p_partition_identity(step7b2, k):
    r = step7b2["drives"][k]
    assert abs(sum(r["Nj"]) - r["N"]) <= PARTITION_RTOL * r["N"], (r["Nj"], r["N"])
    assert abs(sum(r["Dj"]) - r["D"]) <= PARTITION_RTOL * r["D"], (r["Dj"], r["D"])


@complex_only
@pytest.mark.parametrize("k", range(4))
def test_b_legacy_route_reads_eps_c_minus_one(step7b2, k):
    em1 = abs(step7b2["eps_c"] - 1.0)
    r = step7b2["drives"][k]["r_n_legacy"]
    assert abs(r - em1) <= LEGACY_TOL * em1, (r, em1)


@complex_only
@pytest.mark.parametrize("s", ["sx", "sy"])
def test_m1_face_set_is_d2_closed(step7b2, s):
    r = step7b2["m1"][s]
    perm = r["perm"]
    assert r["n_hit_min"] == 1 and r["n_hit_max"] == 1, (r["n_hit_min"], r["n_hit_max"])
    assert np.array_equal(perm[perm], np.arange(perm.size))
    assert not np.any(perm == np.arange(perm.size))
    assert step7b2["F"] % 4 == 0, step7b2["F"]
    assert r["area_rel"] <= EXACT_IDENTITY_RTOL, r["area_rel"]
    assert r["normal_abs"] <= EXACT_IDENTITY_RTOL, r["normal_abs"]


@complex_only
@pytest.mark.parametrize("s,k", G1P_PAIRS)
def test_g1p_whole_interface_d2(step7b2, s, k):
    reading = _g1p(step7b2["drives"], s, k)
    assert reading <= C4_COVARIANCE_BAND, (s, k, reading)


@complex_only
@pytest.mark.parametrize("s,k,j", G2P_PAIRS)
def test_g2p_quadrant_d2(step7b2, s, k, j):
    reading = _g2p(step7b2["drives"], s, k, j)
    assert reading <= C4_COVARIANCE_BAND, (s, k, j, reading)


@complex_only
@pytest.mark.parametrize("s,k", G3P_PAIRS)
def test_g3p_face_level_d2(step7b2, s, k):
    reading = _g3p(step7b2, s, k)
    assert reading <= C4_COVARIANCE_BAND, (s, k, reading)


@complex_only
def test_w0_rank_invariant_face_set(step7b2):
    comm = step7b2["comm"]
    if comm.size == 2:
        pytest.skip("-n 2 is the writer width; W0 is read at -n 4")
    if not os.path.exists(NPZ_PATH):
        pytest.fail(f"{NPZ_PATH} missing: run the -n 2 window first, never recompute")
    ref = np.load(NPZ_PATH)
    d = step7b2["drives"]
    if comm.rank == 0:
        dn = max(abs(d[k]["N"] - ref["N"][k]) / ref["N"][k] for k in range(4))
        dd = max(abs(d[k]["D"] - ref["D"][k]) / ref["D"][k] for k in range(4))
        dr = max(abs(d[k]["r_n"] - ref["r_n"][k]) / ref["r_n"][k] for k in range(4))
        print(
            f"\n[POST-1 7b' W0] -n {comm.size} vs -n 2: max rel N {dn:.3e}, D {dd:.3e}, "
            f"R_n {dr:.3e} (pred <= 1e-8); area diff "
            f"{abs(step7b2['area'] - float(ref['area'])):.3e}",
            flush=True,
        )
    assert int(ref["F"]) == step7b2["F"]
    assert int(ref["V"]) == step7b2["V"]
    assert ref["verts"].shape == step7b2["verts"].shape
    assert np.array_equal(ref["verts"], step7b2["verts"])
