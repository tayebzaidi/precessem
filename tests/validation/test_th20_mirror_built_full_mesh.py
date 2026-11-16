"""TH-20 step 0b (§9 item 198, chain H13-0b): the mirror-built full F-human mesh.

The step-0a quarter (``symmetry_quadrant=True``, LH layout, 1.6 mm fixture) is
gathered to rank 0, mirrored across y = 0 with the σ_y port-tag map and then
across x = 0 with the σ_x map applied to the whole half (so the far quadrant
carries σ_xσ_y), distributed (cells only — interior facet tags hang at -n 2,
known-issues 9), and its facets are re-derived geometrically.

Anchors (asserted; bands imported):
  (M0) cells = 4 × quarter; global vertices = the builder's count; owned census
       sums to size_global; volume = 4·V_q at EXACT_IDENTITY_RTOL.
  (M1) tag set == FULL_TAGS; tags 2 / 3 carry 4 × the quarter's counts; every
       port tag's count = its unique preimage's; preimage map a bijection onto
       the 80 port tags.
  (M2) all 80 port volumes = CAD (imported) and = the direct full mesh at EXACT.
  (M3) D2 port map on cell-volume-weighted centroids: mirror vs direct full mesh
       within EXACT·R; the realised (quarter tag, g) -> t is the closed form.
  (M4) 0 exterior facets on x = 0 / y = 0; exterior = outer ∪ 401;
       N401 = 4 × quarter; A401 = 4·A401_q and outer = 8r_e² + 16 r_e z_e at
       EXACT_IDENTITY_RTOL; every 401 neighbour in {2} ∪ ring-port tags.
  (M5) all 32 ring sheets: facet count = preimage's, area = chord·w at
       SHEET_AREA_BAND, centroid = direct full mesh's within EXACT·R.
Negative controls (rule (e); printed in W1, asserted under
FEM_EM_TH20_STEP0B_CONTROL=1):
  (C4) identity relabel (TH-16's unpermuted copy): 22 tags, and each kept port
       tag's pooled centroid misses the direct one by exactly its radius ρ_f.
  (C5) no vertex reuse on the planes: 4 × 1825 / 4 × 1811 once-seen faces on
       x = 0 / y = 0; the real build reads 0 / 0.

`_gather_owned` and `_distribute` are lifted verbatim from
tests/validation/test_th16_quarter_birdcage.py:234–252 and
tests/validation/test_th16_symmetry_plane_driven.py:145–154 (that module's
imports are complex-only); the mirror is TH-16's `_mirror_by_index` with the
tag copy permuted per copy (the `np.concatenate([tags, tags])` is C4).
"""

from __future__ import annotations

import math
import os
import resource
import time

import numpy as np
import pytest
from mpi4py import MPI

from tests.validation.test_th20_quarter_f_human_mesh import (
    CAD_LEG_BOX,
    CAD_RING_LOWER,
    CAD_RING_UPPER,
    CAVITY_TAG,
    FULL_HOLE_CELL_RECORD,
    FULL_TAGS,
    PLANE_TOL,
    QUARTER_TAGS,
    _extents,
    _meshes,
    _rank0_print,
)

STEP0B_ENV = "FEM_EM_TH20_STEP0B"
CONTROL = os.environ.get("FEM_EM_TH20_STEP0B_CONTROL") == "1"
FACE_TOL = 1.0e-9  # box-face classification tolerance (m)
PORT_TAGS = FULL_TAGS - {2, 3}
RING_PORT_TAGS = set(range(117, 149)) | set(range(217, 249))
# 0a cut records (20261002T185456Z_TH-20-step0a-W1.log:8160).
CUT_FACES_X0 = 1825
CUT_FACES_Y0 = 1811
CLOSED = {
    "sy": lambda j: 15 - j,
    "sx": lambda j: (7 - j) % 16,
    "sxsy": lambda j: (j + 8) % 16,
}
OPS = {
    "e": np.array([1.0, 1.0, 1.0]),
    "sy": np.array([1.0, -1.0, 1.0]),
    "sx": np.array([-1.0, 1.0, 1.0]),
    "sxsy": np.array([-1.0, -1.0, 1.0]),
}


def _skip_unless_on():
    if os.environ.get(STEP0B_ENV, "").strip() != "1":
        pytest.skip(f"{STEP0B_ENV} unset: TH-20 step 0b is selected by environment only")


def _d2_tag(t, g):
    """Closed-form D2 port-tag map; inner/outer (1xx / 2xx) never swap."""
    t = int(t)
    if t in (2, 3) or g == "e":
        return t
    b = 100 if t < 200 else 200
    o = t - b
    s = 1 if o <= 16 else (17 if o <= 32 else 33)
    return b + s + CLOSED[g](o - s)


def _map_tags(tags, g):
    lut = np.arange(256, dtype=np.int64)
    for t in range(256):
        if t in (2, 3) or (101 <= t <= 148) or (217 <= t <= 248):
            lut[t] = _d2_tag(t, g)
    return lut[tags].astype(np.int32)


def _gather_owned(msh, cells):
    """Rank 0: (x, cells, tags) of the whole mesh in global geometry indices."""
    comm = msh.comm
    gmap = msh.geometry.index_map()
    tdim = msh.topology.dim
    imap = msh.topology.index_map(tdim)
    n_cells = imap.size_local
    dm = np.asarray(msh.geometry.dofmap, dtype=np.int32)[:n_cells]
    gcells = np.asarray(gmap.local_to_global(dm.reshape(-1)), dtype=np.int64).reshape(-1, 4)
    tag_of_cell = np.full(imap.size_local + imap.num_ghosts, -1, dtype=np.int32)
    tag_of_cell[np.asarray(cells.indices)] = np.asarray(cells.values)
    own_tags = tag_of_cell[:n_cells]
    xs = comm.gather(np.array(msh.geometry.x[: gmap.size_local, :3], dtype=np.float64), root=0)
    cs = comm.gather(gcells, root=0)
    ts = comm.gather(own_tags, root=0)
    if comm.rank != 0:
        return None
    return np.vstack(xs), np.vstack(cs), np.concatenate(ts)


def _mirror(x, cells, tags, axis, g, reuse=True):
    """Mirror across {x_axis = 0}; on-plane vertices reused by index unless
    ``reuse`` is False (C5); the copy's tags mapped by ``g`` (None = identity, C4)."""
    on = np.abs(x[:, axis]) < PLANE_TOL if reuse else np.zeros(x.shape[0], dtype=bool)
    n = x.shape[0]
    new_index = np.arange(n, dtype=np.int64)
    new_index[~on] = n + np.arange(int(np.count_nonzero(~on)), dtype=np.int64)
    xm = x[~on].copy()
    xm[:, axis] = -xm[:, axis]
    copy_tags = tags if g is None else _map_tags(tags, g)
    return (
        np.vstack([x, xm]),
        np.vstack([cells, new_index[cells][:, [1, 0, 2, 3]]]),
        np.concatenate([tags, copy_tags]),
    )


def _distribute(x, cells, comm):
    import basix.ufl
    import ufl
    from dolfinx import mesh as dmesh

    if comm.rank != 0:
        x = np.empty((0, 3), dtype=np.float64)
        cells = np.empty((0, 4), dtype=np.int64)
    domain = ufl.Mesh(basix.ufl.element("Lagrange", "tetrahedron", 1, shape=(3,)))
    partitioner = dmesh.create_cell_partitioner(dmesh.GhostMode.shared_facet, 2)
    return dmesh.create_mesh(
        comm, cells=np.ascontiguousarray(cells), x=np.ascontiguousarray(x), e=domain,
        partitioner=partitioner,
    )


def _array_moments(p, tags):
    """(256, 4) per-tag [volume, V·cx, V·cy, V·cz] from cell node coords p (n, 4, 3)."""
    v = np.abs(np.einsum("ij,ij->i", np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]),
                         p[:, 3] - p[:, 0])) / 6.0
    c = p.mean(axis=1)
    out = np.zeros((256, 4))
    if len(tags):
        out[:, 0] = np.bincount(tags, weights=v, minlength=256)[:256]
        for k in range(3):
            out[:, 1 + k] = np.bincount(tags, weights=v * c[:, k], minlength=256)[:256]
    return out


def _moments(msh, cells, comm):
    tdim = msh.topology.dim
    imap = msh.topology.index_map(tdim)
    n = imap.size_local
    tag_of_cell = np.full(imap.size_local + imap.num_ghosts, -1, dtype=np.int32)
    tag_of_cell[np.asarray(cells.indices)] = np.asarray(cells.values)
    own = tag_of_cell[:n]
    dm = np.asarray(msh.geometry.dofmap, dtype=np.int32)[:n]
    p = msh.geometry.x[dm][:, :, :3]
    ok = own >= 0
    m = comm.allreduce(_array_moments(p[ok], own[ok]), op=MPI.SUM)
    counts = comm.allreduce(np.bincount(own[ok], minlength=256)[:256], op=MPI.SUM)
    return m, counts


def _centroid(m, t):
    return m[t, 1:] / m[t, 0]


def _sheet_stats(msh, cells, ordinals, comm):
    """Per ring ordinal: (area, area-weighted centroid, facet count), owned, reduced."""
    import dolfinx

    from fem_em_solver.io.mesh import _interface_facet_tags
    from tests.mesh.test_birdcage_port_sheets import PORT_LOWER, PORT_UPPER, SHEET_IFACE

    tags = _interface_facet_tags(
        msh, cells, {SHEET_IFACE + o: (PORT_LOWER + o, PORT_UPPER + o) for o in ordinals}
    )
    fdim = msh.topology.dim - 1
    n_owned = msh.topology.index_map(fdim).size_local
    idx = np.asarray(tags.indices)
    val = np.asarray(tags.values)
    out = {}
    for o in ordinals:
        f = idx[(val == SHEET_IFACE + o) & (idx < n_owned)]
        s = np.zeros(5)
        if len(f):
            g = dolfinx.cpp.mesh.entities_to_geometry(
                msh._cpp_object, fdim, f.astype(np.int32), False)
            p = msh.geometry.x[g][:, :3, :]
            a = 0.5 * np.linalg.norm(np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]), axis=1)
            s = np.array([a.sum(), *(a[:, None] * p.mean(axis=1)).sum(axis=0), len(f)])
        s = comm.allreduce(s, op=MPI.SUM)
        out[o] = (float(s[0]), s[1:4] / s[0] if s[0] > 0 else np.full(3, np.nan), int(s[4]))
    return out


def _plane_face_census(x, cells):
    """Once-seen (exterior) faces with all three vertices on x = 0 / y = 0."""
    res = []
    for axis in (0, 1):
        on = np.abs(x[:, axis]) < PLANE_TOL
        faces = []
        for drop in range(4):
            f = np.delete(cells, drop, axis=1)
            faces.append(f[np.all(on[f], axis=1)])
        f = np.sort(np.vstack(faces), axis=1)
        if len(f) == 0:
            res.append(0)
            continue
        _, cnt = np.unique(f, axis=0, return_counts=True)
        res.append(int(np.count_nonzero(cnt == 1)))
    return tuple(res)


def _derived_facets(msh, cells, r_e, z_e, comm):
    """M4 on the distributed mirror mesh: owned exterior facets classified."""
    import dolfinx

    from tests.validation.test_wf12_f_human_hole import _cell_value

    tdim = msh.topology.dim
    fdim = tdim - 1
    n_owned = msh.topology.index_map(fdim).size_local
    ext = np.asarray(dolfinx.mesh.exterior_facet_indices(msh.topology), dtype=np.int32)
    ext = ext[ext < n_owned]
    n_ext = len(ext)
    on_plane = outer_n = cav_n = 0
    outer_a = cav_a = 0.0
    nbr = set()
    if n_ext:
        g = dolfinx.cpp.mesh.entities_to_geometry(msh._cpp_object, fdim, ext, False)
        p = msh.geometry.x[g][:, :3, :]
        a = 0.5 * np.linalg.norm(np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]), axis=1)
        on_plane = int(np.count_nonzero(np.all(np.abs(p[:, :, 0]) < PLANE_TOL, axis=1))
                       + np.count_nonzero(np.all(np.abs(p[:, :, 1]) < PLANE_TOL, axis=1)))
        outer = (np.all(np.abs(np.abs(p[:, :, 0]) - r_e) < FACE_TOL, axis=1)
                 | np.all(np.abs(np.abs(p[:, :, 1]) - r_e) < FACE_TOL, axis=1)
                 | np.all(np.abs(np.abs(p[:, :, 2]) - z_e) < FACE_TOL, axis=1))
        cav = ~outer
        outer_n, cav_n = int(outer.sum()), int(cav.sum())
        outer_a, cav_a = float(a[outer].sum()), float(a[cav].sum())
        f2c = msh.topology.connectivity(fdim, tdim)
        cv = _cell_value(msh, cells)
        nbr = set(int(t) for t in np.unique(cv[f2c.array[f2c.offsets[ext[cav]]]]).tolist())
    nbr = set().union(*comm.allgather(nbr))
    s = lambda v: comm.allreduce(v, op=MPI.SUM)  # noqa: E731
    return {"n_ext": s(n_ext), "on_plane": s(on_plane), "outer_n": s(outer_n),
            "cav_n": s(cav_n), "outer_a": s(outer_a), "cav_a": s(cav_a), "nbr": nbr}


def test_th20_m0_to_m5_mirror_built_full_mesh():
    _skip_unless_on()
    import dolfinx

    from tests.mesh.probe_two_torus_conductor_hole import _exterior_facet_area
    from tests.mesh.test_birdcage_conductor_hole import EXACT
    from tests.mesh.test_birdcage_leg_offset import SHEET_AREA_BAND
    from tests.mesh.test_birdcage_symmetry_quadrant import EXACT_IDENTITY_RTOL
    from tests.validation.test_birdcage_f_human_rung import _ring_ports
    from tests.validation.test_wf12_f_human_hole import _owned_census

    comm = MPI.COMM_WORLD
    m = _meshes()
    full, q = m["full"], m["quarter"]
    kw = q["kwargs"]
    radius = float(kw["ring_radius"])
    tol_c = EXACT * radius
    fails = []

    # ---- build
    t0 = time.perf_counter()
    gathered = _gather_owned(q["msh"], q["cells"])
    x = c = t = None
    c4 = c5 = real_census = None
    if comm.rank == 0:
        xq, cq, tq = gathered
        x, c, t = xq, cq, tq
        for axis, g in ((1, "sy"), (0, "sx")):
            x, c, t = _mirror(x, c, t, axis, g)
        real_census = _plane_face_census(x, c)
        # (C4) identity relabel
        xi, ci, ti = xq, cq, tq
        for axis in (1, 0):
            xi, ci, ti = _mirror(xi, ci, ti, axis, None)
        c4 = (set(int(v) for v in np.unique(ti).tolist()),
              _array_moments(xi[ci], ti))
        del xi, ci, ti
        # (C5) no vertex reuse
        xn, cn, tn = xq, cq, tq
        for axis in (1, 0):
            xn, cn, tn = _mirror(xn, cn, tn, axis, None, reuse=False)
        c5 = _plane_face_census(xn, cn)
        del xn, cn, tn
    builder_vertices = comm.bcast(None if x is None else int(x.shape[0]), root=0)
    real_census, c5 = comm.bcast((real_census, c5), root=0)
    c4 = comm.bcast(c4, root=0)
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
    msh.topology.create_connectivity(fdim, tdim)
    msh.topology.create_entity_permutations()
    comm.Barrier()
    t_build = time.perf_counter() - t0
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0**2
    n_cells = int(msh.topology.index_map(tdim).size_global)
    _rank0_print(comm, f"\n[TH-20 0b] mirror build {t_build:.2f} s, peak RSS rank0 {rss:.2f} "
                 f"GiB, -n {comm.size}; mirror cells {n_cells} vs direct full "
                 f"{full['n_cells']} (record {FULL_HOLE_CELL_RECORD})")

    # ---- (M0)
    census = _owned_census(msh, cells, comm)
    n_vert = int(msh.geometry.index_map().size_global)
    mm, cnt_m = _moments(msh, cells, comm)
    mq, cnt_q = _moments(q["msh"], q["cells"], comm)
    mf, cnt_f = _moments(full["msh"], full["cells"], comm)
    vol_m, vol_q = float(mm[:, 0].sum()), float(mq[:, 0].sum())
    dv = abs(vol_m / (4.0 * vol_q) - 1.0)
    _rank0_print(comm, f"[TH-20 0b] (M0) cells {n_cells} vs 4×{q['n_cells']} = "
                 f"{4 * q['n_cells']}; vertices {n_vert} vs builder {builder_vertices}; census "
                 f"sum {sum(census.values())}; V {vol_m:.12e} vs 4 V_q dev {dv:.3e} "
                 f"(band {EXACT_IDENTITY_RTOL:g})")
    if n_cells != 4 * q["n_cells"]:
        fails.append(f"M0 cells {n_cells}")
    if n_vert != builder_vertices:
        fails.append(f"M0 vertices {n_vert} != {builder_vertices}")
    if sum(census.values()) != n_cells:
        fails.append(f"M0 census {sum(census.values())}")
    if not dv <= EXACT_IDENTITY_RTOL:
        fails.append(f"M0 volume dev {dv:.3e}")

    # ---- (M1) tag census + preimage map
    tagset = set(int(v) for v in np.nonzero(cnt_m)[0].tolist())
    _rank0_print(comm, f"[TH-20 0b] (M1) {len(tagset)} tags; == FULL_TAGS "
                 f"{tagset == FULL_TAGS}; tag 2 {cnt_m[2]} vs 4×{cnt_q[2]}, tag 3 {cnt_m[3]} "
                 f"vs 4×{cnt_q[3]}")
    if tagset != FULL_TAGS:
        fails.append(f"M1 tag set diff {sorted(tagset ^ FULL_TAGS)}")
    for tt in (2, 3):
        if cnt_m[tt] != 4 * cnt_q[tt]:
            fails.append(f"M1 tag {tt} count {cnt_m[tt]} != 4×{cnt_q[tt]}")
    pre = {}
    for tq_ in sorted(QUARTER_TAGS - {2, 3}):
        for g in OPS:
            pre.setdefault(_d2_tag(tq_, g), []).append((tq_, g))
    if set(pre) != PORT_TAGS or any(len(v) != 1 for v in pre.values()):
        fails.append("M1 preimage map not a bijection onto the 80 port tags")
    bad_cnt = [tt for tt, ((tq_, _g),) in pre.items() if cnt_m[tt] != cnt_q[tq_]]
    _rank0_print(comm, f"[TH-20 0b] (M1) preimage map: {len(pre)} images, bijection "
                 f"{set(pre) == PORT_TAGS and all(len(v) == 1 for v in pre.values())}; "
                 f"count mismatches {bad_cnt}")
    if bad_cnt:
        fails.append(f"M1 counts off preimage at {bad_cnt}")

    # ---- (M2) volumes, (M3) centroids
    cad_f = full["diag"]["cad_mass_by_group"]
    worst2 = worst3 = worst3_map = 0.0
    for tt in sorted(PORT_TAGS):
        o = tt - (100 if tt < 200 else 200)
        if tt <= 116:
            name, imp = f"port_P{o}", CAD_LEG_BOX
        elif tt < 200:
            name, imp = f"port_P{o}_lower", CAD_RING_LOWER
        else:
            name, imp = f"port_P{o}_upper", CAD_RING_UPPER
        vm, vf = mm[tt, 0], mf[tt, 0]
        cadv = float(cad_f.get(name, float("nan")))
        d = max(abs(vm / cadv - 1.0), abs(vm / vf - 1.0), abs(cadv / imp - 1.0))
        worst2 = max(worst2, d)
        if not d <= EXACT:
            fails.append(f"M2 tag {tt} dev {d:.3e}")
        miss = float(np.linalg.norm(_centroid(mm, tt) - _centroid(mf, tt)))
        worst3 = max(worst3, miss)
        if not miss <= tol_c:
            fails.append(f"M3 tag {tt} miss {miss:.3e}")
        ((tq_, g),) = pre.get(tt, [(None, None)]) if tt in pre else ((None, None),)
        if tq_ is not None:
            dmap = float(np.linalg.norm(_centroid(mq, tq_) * OPS[g] - _centroid(mm, tt)))
            worst3_map = max(worst3_map, dmap)
            if not dmap <= tol_c:
                fails.append(f"M3 map ({tq_},{g})->{tt} miss {dmap:.3e}")
    _rank0_print(comm, f"[TH-20 0b] (M2) 80 port volumes worst dev {worst2:.3e} (band "
                 f"{EXACT:g}); (M3) worst centroid miss vs direct {worst3:.3e} m, vs "
                 f"g·c_q(preimage) {worst3_map:.3e} m (tol {tol_c:.3e} m)")
    dv3 = mm[3, 0] / mf[3, 0] - 1.0
    _rank0_print(comm, f"[TH-20 0b] (printed) tag-3 volume mirror/direct - 1 = {dv3:+.3e}")

    # ---- (M4) cavity, re-derived
    r_e, z_e = _extents(kw, float(q["diag"]["port_box_size_m"][1]))
    fd = _derived_facets(msh, cells, r_e, z_e, comm)
    a_q = _exterior_facet_area(q["msh"], q["facets"], CAVITY_TAG, comm)
    a_f = _exterior_facet_area(full["msh"], full["facets"], CAVITY_TAG, comm)
    n_q401 = comm.allreduce(int(np.count_nonzero(
        (np.asarray(q["facets"].values) == CAVITY_TAG)
        & (np.asarray(q["facets"].indices) < q["msh"].topology.index_map(fdim).size_local))),
        op=MPI.SUM)
    a_outer = 8.0 * r_e**2 + 16.0 * r_e * z_e
    d_cav = abs(fd["cav_a"] / (4.0 * a_q) - 1.0)
    d_out = abs(fd["outer_a"] / a_outer - 1.0)
    bad_nbr = sorted(fd["nbr"] - ({2} | RING_PORT_TAGS))
    _rank0_print(comm, f"[TH-20 0b] (M4) exterior {fd['n_ext']} = outer {fd['outer_n']} + "
                 f"401 {fd['cav_n']}; on planes {fd['on_plane']}; N401 {fd['cav_n']} vs "
                 f"4×{n_q401} = {4 * n_q401}; A401 {fd['cav_a']:.9e} vs 4 A_q dev {d_cav:.3e}; "
                 f"outer {fd['outer_a']:.12e} vs {a_outer:.12e} dev {d_out:.3e}; 401 "
                 f"neighbours outside {{2}}∪ring-port {bad_nbr}")
    _rank0_print(comm, f"[TH-20 0b] (printed) |A401 mirror/direct - 1| = "
                 f"{abs(fd['cav_a'] / a_f - 1.0):.3e} (Q3 record 1.498e-4, two triangulations)")
    if fd["on_plane"] != 0:
        fails.append(f"M4 {fd['on_plane']} exterior facets on a plane")
    if fd["outer_n"] + fd["cav_n"] != fd["n_ext"]:
        fails.append("M4 partition")
    if fd["cav_n"] != 4 * n_q401:
        fails.append(f"M4 N401 {fd['cav_n']} != {4 * n_q401}")
    if not d_cav <= EXACT_IDENTITY_RTOL:
        fails.append(f"M4 A401 dev {d_cav:.3e}")
    if not d_out <= EXACT_IDENTITY_RTOL:
        fails.append(f"M4 outer dev {d_out:.3e}")
    if bad_nbr:
        fails.append(f"M4 401 neighbours {bad_nbr}")

    # ---- (M5) sheets
    ring_ports = _ring_ports()
    sm = _sheet_stats(msh, cells, ring_ports, comm)
    sf = _sheet_stats(full["msh"], full["cells"], ring_ports, comm)
    kept_ring = sorted(o for o in (tq_ - 100 for tq_ in QUARTER_TAGS if 117 <= tq_ <= 148))
    sq = _sheet_stats(q["msh"], q["cells"], kept_ring, comm)
    sheet_analytic = float(full["diag"]["ring_port_layout"]["ring_port_sheet_longitudinal_area_m2"])
    worst5a = worst5c = 0.0
    for o in ring_ports:
        ((tq_, _g),) = pre[100 + o]
        if sm[o][2] != sq[tq_ - 100][2]:
            fails.append(f"M5 sheet {o} facets {sm[o][2]} != preimage {sq[tq_ - 100][2]}")
        da = abs(sm[o][0] / sheet_analytic - 1.0)
        dc = float(np.linalg.norm(sm[o][1] - sf[o][1]))
        worst5a, worst5c = max(worst5a, da), max(worst5c, dc)
        if not da <= SHEET_AREA_BAND:
            fails.append(f"M5 sheet {o} area dev {da:.3e}")
        if not dc <= tol_c:
            fails.append(f"M5 sheet {o} centroid miss {dc:.3e}")
    _rank0_print(comm, f"[TH-20 0b] (M5) 32 sheets: worst area dev {worst5a:.3e} (band "
                 f"{SHEET_AREA_BAND:g}); worst centroid miss {worst5c:.3e} m (tol {tol_c:.3e})")

    # ---- (C4) identity relabel
    c4_tags, c4_m = c4
    c4_worst = 0.0
    for tq_ in sorted(QUARTER_TAGS - {2, 3}):
        cf = _centroid(mf, tq_)
        rho = math.hypot(cf[0], cf[1])
        miss = float(np.linalg.norm(_centroid(c4_m, tq_) - cf))
        c4_worst = max(c4_worst, abs(miss / rho - 1.0))
    _rank0_print(comm, f"[TH-20 0b] (C4) identity relabel: {len(c4_tags)} tags, "
                 f"{len(FULL_TAGS - c4_tags)} of FULL_TAGS missing; worst |miss/ρ_f - 1| "
                 f"{c4_worst:.3e} (band {EXACT:g}; M3 tol {tol_c:.3e} m)")
    # ---- (C5) no vertex reuse
    _rank0_print(comm, f"[TH-20 0b] (C5) once-seen plane faces x=0 / y=0: no reuse {c5} "
                 f"(predicted {(4 * CUT_FACES_X0, 4 * CUT_FACES_Y0)}); real build "
                 f"{real_census}")
    if real_census != (0, 0):
        fails.append(f"M4 builder plane faces {real_census}")
    if CONTROL:
        if not (len(c4_tags) == 22 and len(FULL_TAGS - c4_tags) == 60 and c4_worst <= EXACT):
            fails.append(f"C4 tags {len(c4_tags)} worst {c4_worst:.3e}")
        if c5 != (4 * CUT_FACES_X0, 4 * CUT_FACES_Y0):
            fails.append(f"C5 {c5}")

    _rank0_print(comm, f"[TH-20 0b] failures: {fails}")
    assert not fails, fails
