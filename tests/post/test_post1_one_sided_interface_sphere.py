"""`POST-1` step 7a — one-sided interface residuals on the lossy sphere.

PROJECT_PLAN §9 item 187; brief ``docs/planning/chunks/POST-1.md``.  The
`TH-10` fixture (saline sphere a = 0.05 m, ε_r = 78, σ = 0.5 S/m, imposed
series field on the box wall), 64 MHz, N1curl degree 1, ``fields.e_complex``.

Interface points are the centroids of every face of an owned tag-1 cell whose
three vertices lie on |x| = a (to 1e-9·a); they are allgathered and sorted by
(inside cell's ``original_cell_index``, sorted vertex coordinates), so the
array is bitwise identical on every rank.  E⁻ / E⁺ come from the one-sided
owner rule (``candidate_cells = cell_tags.find(1)`` / ``find(2)``).

Series values at the centroid p: S⁻ = ``internal_field(p)``, S⁺ =
``incident_field(p) + scattered_field(p)`` — never ``total_field`` (it
switches on r < a and every centroid sits inside by the chord sagitta).

Residuals: R_n = ‖ε_c E⁻·n − E⁺·n‖₂ / ‖E⁺·n‖₂;  e∓ = ‖E∓ − S∓‖₂ / ‖S∓‖₂.

Asserted: (S0) the series premise at p̂ = a·p/|p| to 1e-9; (X0a) F = 2V − 4;
(X0b) one owned containing cell per side per point (global SUM); (X1) the
tangential trace is continuous to 1e-10·max|E⁺|; (C) R_n, e⁻, e⁺ strictly
decrease over the three rungs; (B) the legacy ``owner=None`` route on both
sides gives R_n = |ε_c − 1| to 1e-12 relative.  Printed, not asserted: levels,
pairwise rates, separation |ε_c − 1|/R_n, σ-blind and ε_c* controls, legacy
side fraction.  Scope: no peak, SAR, birdcage or ✅ claim; no band on levels.
"""

from __future__ import annotations

import numpy as np
from mpi4py import MPI

from fem_em_solver.post.evaluation import evaluate_vector_field_parallel
from tests.complex_mode import complex_only
from tests.validation.test_lossy_sphere_fullwave import (
    RESOLUTIONS_64,
    RESOLUTIONS_128,
    SPHERE_TAG,
    _mesh_and_solve,
    _series,
)
from dolfinx import geometry

AIR_TAG = 2  # io/mesh.py sphere_in_box_domain: 1 sphere, 2 air
FREQUENCY_HZ = 64.0e6
RUNGS = list(RESOLUTIONS_64) + [RESOLUTIONS_128[-1]]
ON_SPHERE_TOL = 1e-9  # × a; never tuned (brief trap)
SERIES_PREMISE_TOL = 1e-9
TANGENTIAL_TOL = 1e-10
LEGACY_TOL = 1e-12
_FACES = ((1, 2, 3), (0, 2, 3), (0, 1, 3), (0, 1, 2))


def _original(msh):
    oci = msh.topology.original_cell_index
    if isinstance(oci, (list, tuple)):
        oci = oci[0]
    return np.asarray(oci, dtype=np.int64)


def _interface_points(msh, cell_tags, radius, comm):
    """Allgathered, sorted face centroids, normals and surface-vertex count.

    Returns ``(points, normals, n_faces, n_vertices)``; identical on all ranks.
    """
    n_owned = msh.topology.index_map(3).size_local
    sphere_cells = cell_tags.find(SPHERE_TAG)
    sphere_cells = sphere_cells[sphere_cells < n_owned]
    original = _original(msh)
    dofmap = msh.geometry.dofmap
    x = msh.geometry.x[:, :3]
    input_ids = np.asarray(msh.geometry.input_global_indices, dtype=np.int64)

    keys, verts_out, normals, vertex_ids = [], [], [], []
    for c in sphere_cells:
        nodes = np.asarray(dofmap[c])[:4]
        on = np.abs(np.linalg.norm(x[nodes], axis=1) - radius) <= ON_SPHERE_TOL * radius
        for face in _FACES:
            if not all(on[j] for j in face):
                continue
            fn = nodes[list(face)]
            v = x[fn]
            order = np.lexsort((v[:, 2], v[:, 1], v[:, 0]))
            v = v[order]
            p = v.mean(axis=0)
            n = np.cross(v[1] - v[0], v[2] - v[0])
            n = n / np.linalg.norm(n)
            if float(np.dot(n, p)) < 0.0:
                n = -n
            keys.append(int(original[c]))
            verts_out.append(v.reshape(-1))
            normals.append(n)
            vertex_ids.extend(int(input_ids[k]) for k in fn)

    def _cat(lst, shape):
        return np.asarray(lst, dtype=np.float64).reshape(shape)

    g_keys = np.concatenate(
        [np.asarray(k, dtype=np.int64) for k in comm.allgather(keys)]
    )
    g_verts = np.concatenate([_cat(v, (-1, 9)) for v in comm.allgather(verts_out)])
    g_normals = np.concatenate([_cat(v, (-1, 3)) for v in comm.allgather(normals)])
    g_vids = np.concatenate(
        [np.asarray(v, dtype=np.int64) for v in comm.allgather(vertex_ids)]
    )

    sort_keys = tuple(g_verts[:, j] for j in range(8, -1, -1)) + (g_keys,)
    order = np.lexsort(sort_keys)
    g_verts = g_verts[order]
    g_normals = np.ascontiguousarray(g_normals[order])
    points = np.ascontiguousarray(g_verts.reshape(-1, 3, 3).mean(axis=1))
    return points, g_normals, int(g_keys.size), int(np.unique(g_vids).size)


def _side_counts(msh, cell_tags, points, comm):
    """Global count of owned containing cells of tag 1 and tag 2 per point."""
    tdim = msh.topology.dim
    imap = msh.topology.index_map(tdim)
    n_owned = imap.size_local
    tag_of = np.zeros(n_owned + imap.num_ghosts, dtype=np.int32)
    tag_of[cell_tags.indices] = cell_tags.values
    bb = geometry.bb_tree(msh, tdim)
    hits = geometry.compute_collisions_points(bb, points)
    coll = geometry.compute_colliding_cells(msh, hits, points)
    local = np.zeros((points.shape[0], 2), dtype=np.int64)
    for i in range(points.shape[0]):
        links = np.asarray(coll.links(i))
        owned = links[links < n_owned]
        local[i, 0] = int(np.sum(tag_of[owned] == SPHERE_TAG))
        local[i, 1] = int(np.sum(tag_of[owned] == AIR_TAG))
    total = np.empty_like(local)
    comm.Allreduce(local, total, op=MPI.SUM)
    return total


def _rel(a, b):
    return float(np.linalg.norm(a - b) / np.linalg.norm(b))


def _rung(series, hs, hf, comm):
    msh, cell_tags, fields, ncells = _mesh_and_solve(series, hs, hf)
    a = series.radius
    eps_c = series.epsilon_c
    points, normals, n_faces, n_vertices = _interface_points(msh, cell_tags, a, comm)
    counts = _side_counts(msh, cell_tags, points, comm)

    e = fields.e_complex
    e_in, v_in = evaluate_vector_field_parallel(
        e, points, comm, owner="lowest_original_cell",
        candidate_cells=cell_tags.find(SPHERE_TAG),
    )
    e_out, v_out = evaluate_vector_field_parallel(
        e, points, comm, owner="lowest_original_cell",
        candidate_cells=cell_tags.find(AIR_TAG),
    )
    e_leg, v_leg = evaluate_vector_field_parallel(e, points, comm)
    e_in, e_out, e_leg = (np.asarray(z, dtype=complex) for z in (e_in, e_out, e_leg))

    s_in = series.internal_field(points)
    s_out = series.incident_field(points) + series.scattered_field(points)

    # (S0) premise at the projected point p̂ = a·p/|p|.
    p_hat = a * points / np.linalg.norm(points, axis=1)[:, None]
    r_hat = p_hat / a
    s_in_hat = series.internal_field(p_hat)
    s_out_hat = series.incident_field(p_hat) + series.scattered_field(p_hat)
    s0_normal = float(
        np.linalg.norm(np.sum((eps_c * s_in_hat - s_out_hat) * r_hat, axis=1))
        / np.linalg.norm(np.sum(s_out_hat * r_hat, axis=1))
    )
    s0_tangent = float(
        np.linalg.norm(np.cross(r_hat, s_in_hat - s_out_hat))
        / np.linalg.norm(s_out_hat)
    )

    en_in = np.sum(e_in * normals, axis=1)
    en_out = np.sum(e_out * normals, axis=1)
    en_leg = np.sum(e_leg * normals, axis=1)
    r_n = float(np.linalg.norm(eps_c * en_in - en_out) / np.linalg.norm(en_out))
    r_n_blind = float(
        np.linalg.norm(eps_c.real * en_in - en_out) / np.linalg.norm(en_out)
    )
    r_n_conj = float(
        np.linalg.norm(np.conj(eps_c) * en_in - en_out) / np.linalg.norm(en_out)
    )
    r_n_legacy = float(np.linalg.norm(eps_c * en_leg - en_leg) / np.linalg.norm(en_leg))
    legacy_tag1 = float(
        np.mean(
            np.linalg.norm(e_leg - e_in, axis=1) < np.linalg.norm(e_leg - e_out, axis=1)
        )
    )
    x1 = float(np.max(np.linalg.norm(np.cross(normals, e_in - e_out), axis=1)))
    x1_scale = float(np.max(np.linalg.norm(e_out, axis=1)))

    return {
        "hs": hs,
        "ncells": ncells,
        "n_points": points.shape[0],
        "F": n_faces,
        "V": n_vertices,
        "counts": counts,
        "valid": (v_in, v_out, v_leg),
        "s0": (s0_normal, s0_tangent),
        "r_n": r_n,
        "r_n_blind": r_n_blind,
        "r_n_conj": r_n_conj,
        "r_n_legacy": r_n_legacy,
        "legacy_tag1": legacy_tag1,
        "e_in": _rel(e_in, s_in),
        "e_out": _rel(e_out, s_out),
        "x1": x1,
        "x1_scale": x1_scale,
        "eps_c": eps_c,
    }


@complex_only
def test_post1_one_sided_interface_residuals_sphere_64mhz():
    comm = MPI.COMM_WORLD
    series = _series(FREQUENCY_HZ)
    eps_c = series.epsilon_c
    eps_minus_one = abs(eps_c - 1.0)
    runs = [_rung(series, hs, hf, comm) for hs, hf in RUNGS]

    if comm.rank == 0:
        print(
            f"\n[POST-1 7a] -n {comm.size}, lossy sphere 64 MHz, degree 1; "
            f"eps_c = {eps_c!r}, |eps_c - 1| = {eps_minus_one!r}"
        )
        for r in runs:
            c = r["counts"]
            one_sided = int(np.sum((c[:, 0] == 0) ^ (c[:, 1] == 0)))
            print(
                f"  h_sphere {r['hs']:.5f} ({r['ncells']:6d} cells): points {r['n_points']}, "
                f"F {r['F']}, V {r['V']}, 2V-4 {2 * r['V'] - 4}; "
                f"side counts tag1 [{int(c[:, 0].min())},{int(c[:, 0].max())}] "
                f"tag2 [{int(c[:, 1].min())},{int(c[:, 1].max())}], "
                f"one-sided points {one_sided}"
            )
            print(
                f"    S0 normal {r['s0'][0]:.3e}, tangential {r['s0'][1]:.3e}; "
                f"X1 max|n x (E- - E+)| {r['x1']:.3e} vs max|E+| {r['x1_scale']:.4e} "
                f"(ratio {r['x1'] / r['x1_scale']:.3e})"
            )
            print(
                f"    R_n {r['r_n']!r}  e- {r['e_in']!r}  e+ {r['e_out']!r}"
            )
            print(
                f"    controls: sigma-blind {r['r_n_blind']:.4f} (pred 0.874), "
                f"eps_c* {r['r_n_conj']:.4f} (pred 1.748), legacy R_n "
                f"{r['r_n_legacy']!r} (rel to |eps_c-1| "
                f"{abs(r['r_n_legacy'] - eps_minus_one) / eps_minus_one:.3e}), "
                f"legacy tag-1 fraction {r['legacy_tag1']:.3f} (pred ~0.5), "
                f"separation |eps_c-1|/R_n {eps_minus_one / r['r_n']:.1f}x (pred >= 500)"
            )
        for lo, hi in zip(runs[:-1], runs[1:]):
            h_ratio = lo["hs"] / hi["hs"]
            rates = {
                k: float(np.log(lo[k] / hi[k]) / np.log(h_ratio))
                for k in ("r_n", "e_in", "e_out")
            }
            print(
                f"  rate h {lo['hs']:.5f} -> {hi['hs']:.5f}: "
                + ", ".join(f"{k} {v:.3f}" for k, v in rates.items())
                + " (pred ~1, printed)"
            )

    for r in runs:
        # (S0)
        assert r["s0"][0] <= SERIES_PREMISE_TOL, r["s0"]
        assert r["s0"][1] <= SERIES_PREMISE_TOL, r["s0"]
        # (X0a)
        assert r["F"] == 2 * r["V"] - 4, (r["F"], r["V"])
        assert r["n_points"] == r["F"]
        # (X0b)
        assert np.all(r["counts"][:, 0] == 1), r["counts"][:, 0]
        assert np.all(r["counts"][:, 1] == 1), r["counts"][:, 1]
        for v in r["valid"]:
            assert bool(np.all(v))
        # (X1)
        assert r["x1"] <= TANGENTIAL_TOL * r["x1_scale"], (r["x1"], r["x1_scale"])
        # (B)
        assert abs(r["r_n_legacy"] - eps_minus_one) <= LEGACY_TOL * eps_minus_one, (
            r["r_n_legacy"], eps_minus_one,
        )
    # (C) strictly decreasing over the three rungs
    for key in ("r_n", "e_in", "e_out"):
        vals = [r[key] for r in runs]
        assert all(b < a for a, b in zip(vals[:-1], vals[1:])), (key, vals)
