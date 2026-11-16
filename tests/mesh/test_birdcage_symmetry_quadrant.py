"""`TH-16` step 3a — the quadrant cutter on the `GEO-19` gapped birdcage.

``birdcage_port_domain(symmetry_quadrant=True)`` meshes only x >= 0, y >= 0:
the air box is [0, r_e] x [0, r_e] x [-z_e, z_e], and the kept port-box
halves' faces on y = 0 (P1) and x = 0 (P2) are the leg-port sheets.  Every face
involved here is planar, so the anchors are exact identities, not bands:

(0) cell tags {1, 2, 3, 201, 102} — one half of P1, one of P2, nothing of
    P3 / P4 ("upper" is the +phi side, `io/mesh.py` port-half rule); every
    vertex at x, y >= -1e-12 (allreduced);
(1) meshed volume r_e² · 2 z_e = 2.880e-3 m³, each kept half
    0.014 · 0.014 · 0.008 / 2 = 7.840e-7 m³, to ``EXACT_IDENTITY_RTOL``;
(2) each cut group has area r_e · 2 z_e = 2.400e-2 m², lies on its plane to
    1e-12 m, and is exactly the owned exterior facets on that plane; the outer
    group (tag 1) is 2 r_e · 2 z_e + 2 r_e² = 7.680e-2 m² and shares no facet
    with a cut group;
(3) the sheets — the kept P1 (P2) half's exterior facets on y = 0 (x = 0) —
    have area dx · g = 1.120e-4 m², and A / h equals their bbox extent;
(4) flag off: the `GEO-19` record, 116 085 cells exactly.

Negative control: the same plane checker on the flag-off full mesh.  Predicted
|V/V_q - 1| = 3 (four quadrants), > 0 vertices outside the quadrant, 0 exterior
facets on either plane.  Printed always; asserted only under
``FEM_EM_TH16_STEP3A_CONTROL=1``.

Mesh only: no solve, no port spec, no PEC / PMC assignment, no memory claim.
"""

from __future__ import annotations

import os

import dolfinx
import numpy as np
from mpi4py import MPI

from fem_em_solver.io.mesh import BIRDCAGE_CUT_X0_TAG, BIRDCAGE_CUT_Y0_TAG
from tests.mesh.helpers import global_cell_tag_set
from tests.mesh.test_birdcage_port_sheets import _build
from tests.mesh.test_coil_phantom_conforming import _tag_volume, _total_volume
from tests.mesh.test_birdcage_leg_gaps import LEG_GAP_LENGTH

CONTROL = os.environ.get("FEM_EM_TH16_STEP3A_CONTROL") == "1"

# Imported value (`tests/validation/test_birdcage_sar_mass_averaged.py:163`),
# restated to avoid importing a validation module into a mesh test.
EXACT_IDENTITY_RTOL = 1.0e-10
# `GEO-19` record of the flag-off mesh at `-n 2`
# (`20260929T110818Z_PORT-8-W1.log:1866`).
FULL_CELL_RECORD = 116085
PLANE_TOL = 1.0e-12

# Fixture closed forms (`io/mesh.py` extents: r_e = 0.07 + 0.006 + 0.014 +
# 0.03, z_e = 0.07 + 0.03; gap box 0.014 x 0.014 x 0.008).
R_E = 0.120
Z_E = 0.100
PORT_DX = 0.014
V_QUARTER = R_E * R_E * 2.0 * Z_E            # 2.880e-3 m^3
V_HALF_PORT = 0.5 * PORT_DX * PORT_DX * LEG_GAP_LENGTH  # 7.840e-7 m^3
A_CUT = R_E * 2.0 * Z_E                      # 2.400e-2 m^2
A_OUTER = 2.0 * R_E * 2.0 * Z_E + 2.0 * R_E * R_E  # 7.680e-2 m^2
A_SHEET = PORT_DX * LEG_GAP_LENGTH           # 1.120e-4 m^2

KEPT_HALVES = {1: 201, 2: 102}  # predicted: P1 upper (+y), P2 lower (+x)

_CACHE: dict = {}


def _mesh(quadrant: bool):
    if quadrant not in _CACHE:
        _CACHE[quadrant] = _build(True, symmetry_quadrant=quadrant)
    return _CACHE[quadrant]


def _facet_nodes(mesh, facets):
    fdim = mesh.topology.dim - 1
    if len(facets) == 0:
        return np.zeros((0, 3, 3))
    nodes = dolfinx.cpp.mesh.entities_to_geometry(
        mesh._cpp_object, fdim, np.asarray(facets, dtype=np.int32), False
    )
    return mesh.geometry.x[nodes][:, :3, :]


def _tri_area(pts):
    if len(pts) == 0:
        return 0.0
    return float(
        0.5 * np.linalg.norm(np.cross(pts[:, 1] - pts[:, 0], pts[:, 2] - pts[:, 0]), axis=1).sum()
    )


def _plane_readings(mesh, cells, facet_tags):
    """Collective. Owned exterior facets on x = 0 / y = 0, tagged sets, areas."""
    comm = mesh.comm
    tdim = mesh.topology.dim
    fdim = tdim - 1
    mesh.topology.create_connectivity(fdim, tdim)
    n_owned = mesh.topology.index_map(fdim).size_local
    ext = np.asarray(dolfinx.mesh.exterior_facet_indices(mesh.topology), dtype=np.int32)
    ext = ext[ext < n_owned]
    pts = _facet_nodes(mesh, ext)
    f2c = mesh.topology.connectivity(fdim, tdim)
    n_cells_local = mesh.topology.index_map(tdim).size_local + mesh.topology.index_map(tdim).num_ghosts
    tag_of_cell = np.full(n_cells_local, -1, dtype=np.int64)
    tag_of_cell[np.asarray(cells.indices)] = np.asarray(cells.values)
    ext_cell_tag = np.array(
        [tag_of_cell[f2c.links(int(f))[0]] for f in ext], dtype=np.int64
    )

    tag_idx = np.asarray(facet_tags.indices)
    tag_val = np.asarray(facet_tags.values)
    out = {}
    for axis, cut_tag in ((0, BIRDCAGE_CUT_X0_TAG), (1, BIRDCAGE_CUT_Y0_TAG)):
        on_plane = (
            np.all(np.abs(pts[:, :, axis]) < PLANE_TOL, axis=1)
            if len(ext) else np.zeros(0, dtype=bool)
        )
        geo_set = set(ext[on_plane].tolist())
        tagged = tag_idx[(tag_val == cut_tag) & (tag_idx < n_owned)]
        tagged_set = set(tagged.tolist())
        tagged_pts = _facet_nodes(mesh, tagged)
        flat = float(np.abs(tagged_pts[:, :, axis]).max()) if len(tagged) else 0.0
        out[cut_tag] = {
            "n_geo": comm.allreduce(len(geo_set), op=MPI.SUM),
            "n_tagged": comm.allreduce(len(tagged_set), op=MPI.SUM),
            "mismatch": comm.allreduce(len(geo_set ^ tagged_set), op=MPI.SUM),
            "area": comm.allreduce(_tri_area(tagged_pts), op=MPI.SUM),
            "flat": comm.allreduce(flat, op=MPI.MAX),
            "on_plane_mask": on_plane,
        }
    outer = tag_idx[(tag_val == 1) & (tag_idx < n_owned)]
    cut_all = tag_idx[
        ((tag_val == BIRDCAGE_CUT_X0_TAG) | (tag_val == BIRDCAGE_CUT_Y0_TAG))
        & (tag_idx < n_owned)
    ]
    out["outer_area"] = comm.allreduce(_tri_area(_facet_nodes(mesh, outer)), op=MPI.SUM)
    out["outer_cut_overlap"] = comm.allreduce(
        len(set(outer.tolist()) & set(cut_all.tolist())), op=MPI.SUM
    )
    x = mesh.geometry.x
    n_geom_owned = mesh.geometry.index_map().size_local
    xo = x[:n_geom_owned]
    out["n_outside"] = comm.allreduce(
        int(np.count_nonzero((xo[:, 0] < -PLANE_TOL) | (xo[:, 1] < -PLANE_TOL))),
        op=MPI.SUM,
    )
    out["min_xy"] = comm.allreduce(float(xo[:, :2].min()) if len(xo) else np.inf, op=MPI.MIN)
    out["_ext"] = ext
    out["_pts"] = pts
    out["_ext_cell_tag"] = ext_cell_tag
    return out


def _sheet(mesh, readings, axis, cut_tag, cell_tag):
    """Collective. Area and bbox extents of the kept half's facets on the plane."""
    comm = mesh.comm
    mask = readings[cut_tag]["on_plane_mask"] & (readings["_ext_cell_tag"] == cell_tag)
    pts = readings["_pts"][mask]
    area = comm.allreduce(_tri_area(pts), op=MPI.SUM)
    flat = pts.reshape(-1, 3)
    lo = flat.min(axis=0) if len(flat) else np.full(3, np.inf)
    hi = flat.max(axis=0) if len(flat) else np.full(3, -np.inf)
    lo = np.min(np.array(comm.allgather(lo)), axis=0)
    hi = np.max(np.array(comm.allgather(hi)), axis=0)
    n = comm.allreduce(int(mask.sum()), op=MPI.SUM)
    return area, hi - lo, n


def test_quadrant_cutter_planar_identities():
    comm = MPI.COMM_WORLD
    q_mesh, q_cells, q_facets, q_diag, q_elapsed = _mesh(True)
    f_mesh, f_cells, _, f_diag, f_elapsed = _mesh(False)

    q_tags = global_cell_tag_set(q_mesh, q_cells)
    n_q = q_mesh.topology.index_map(3).size_global
    n_f = f_mesh.topology.index_map(3).size_global
    v_q = _total_volume(q_mesh, comm)
    v_half = {o: _tag_volume(q_mesh, q_cells, t, comm) for o, t in KEPT_HALVES.items()}
    r = _plane_readings(q_mesh, q_cells, q_facets)
    sheets = {
        1: _sheet(q_mesh, r, 1, BIRDCAGE_CUT_Y0_TAG, KEPT_HALVES[1]),
        2: _sheet(q_mesh, r, 0, BIRDCAGE_CUT_X0_TAG, KEPT_HALVES[2]),
    }

    if comm.rank == 0:
        print(
            f"\n[TH-16 step 3a] quarter cells={n_q}  full cells={n_f} "
            f"(record {FULL_CELL_RECORD})  ratio={n_q / FULL_CELL_RECORD:.6f} "
            f"(predicted 0.25-0.30)  quarter mesh={q_diag['mesh_wall_time_s']:.2f} s "
            f"rung={q_elapsed:.2f} s  full mesh={f_diag['mesh_wall_time_s']:.2f} s "
            f"rung={f_elapsed:.2f} s"
            f"\n[TH-16 step 3a] (0) cell tags {sorted(q_tags)} (predicted "
            f"[1, 2, 3, 102, 201]); vertices outside quadrant={r['n_outside']} "
            f"min(x,y)={r['min_xy']:.3e}"
            f"\n[TH-16 step 3a] (1) V={v_q:.12e} m^3 rel={v_q / V_QUARTER - 1:.3e}; "
            + " ".join(
                f"P{o} half tag {KEPT_HALVES[o]} V={v:.12e} rel={v / V_HALF_PORT - 1:.3e}"
                for o, v in v_half.items()
            ),
            flush=True,
        )
        for cut_tag in (BIRDCAGE_CUT_X0_TAG, BIRDCAGE_CUT_Y0_TAG):
            d = r[cut_tag]
            print(
                f"[TH-16 step 3a] (2) cut {cut_tag}: area={d['area']:.12e} "
                f"rel={d['area'] / A_CUT - 1:.3e} flat={d['flat']:.3e} m "
                f"owned exterior on plane={d['n_geo']} tagged={d['n_tagged']} "
                f"symmetric difference={d['mismatch']}",
                flush=True,
            )
        print(
            f"[TH-16 step 3a] (2) outer tag 1 area={r['outer_area']:.12e} "
            f"rel={r['outer_area'] / A_OUTER - 1:.3e}; facets in both outer and cut="
            f"{r['outer_cut_overlap']}",
            flush=True,
        )
        for o, (area, ext, n) in sheets.items():
            h = ext[2]
            w_bbox = ext[0] if o == 1 else ext[1]
            print(
                f"[TH-16 step 3a] (3) P{o} sheet: {n} facets area={area:.12e} "
                f"rel={area / A_SHEET - 1:.3e} h={h:.9e} w_bbox={w_bbox:.9e} "
                f"A/h/w_bbox-1={area / h / w_bbox - 1:.3e}",
                flush=True,
            )
        cad_q = q_diag["cad_mass_by_group"]
        cad_f = f_diag["cad_mass_by_group"]
        print(
            "[TH-16 step 3a] CAD masses x4 / flag-off: "
            + " ".join(
                f"{k}={4.0 * cad_q[k] / cad_f[k]:.9f}"
                for k in ("conductor", "air", "phantom")
            )
            + f"; generator port sheets (CAD) {q_diag['port_sheet_cad']}",
            flush=True,
        )

    # (4) flag off is the record mesh.
    assert n_f == FULL_CELL_RECORD
    # (0)
    assert q_tags == {1, 2, 3, 201, 102}
    assert r["n_outside"] == 0
    # (1)
    assert abs(v_q / V_QUARTER - 1.0) < EXACT_IDENTITY_RTOL
    for o, v in v_half.items():
        assert abs(v / V_HALF_PORT - 1.0) < EXACT_IDENTITY_RTOL, (o, v)
    # (2)
    for cut_tag in (BIRDCAGE_CUT_X0_TAG, BIRDCAGE_CUT_Y0_TAG):
        d = r[cut_tag]
        assert abs(d["area"] / A_CUT - 1.0) < EXACT_IDENTITY_RTOL, (cut_tag, d["area"])
        assert d["flat"] <= PLANE_TOL
        assert d["n_tagged"] > 0 and d["n_geo"] == d["n_tagged"] and d["mismatch"] == 0
    assert abs(r["outer_area"] / A_OUTER - 1.0) < EXACT_IDENTITY_RTOL
    assert r["outer_cut_overlap"] == 0
    # (3)
    for o, (area, ext, n) in sheets.items():
        assert n > 0
        assert abs(area / A_SHEET - 1.0) < EXACT_IDENTITY_RTOL, (o, area)
        h = ext[2]
        w_bbox = ext[0] if o == 1 else ext[1]
        assert abs(h / LEG_GAP_LENGTH - 1.0) < 1e-9
        assert abs(area / h / w_bbox - 1.0) < 1e-9


def test_negative_control_full_mesh_fails_the_checker():
    """The plane checker on the flag-off full mesh: predicted to fail all three."""
    comm = MPI.COMM_WORLD
    f_mesh, f_cells, f_facets, _, _ = _mesh(False)
    v_f = _total_volume(f_mesh, comm)
    r = _plane_readings(f_mesh, f_cells, f_facets)
    ratio_dev = abs(v_f / V_QUARTER - 1.0)
    n_on_planes = r[BIRDCAGE_CUT_X0_TAG]["n_geo"] + r[BIRDCAGE_CUT_Y0_TAG]["n_geo"]
    if comm.rank == 0:
        print(
            f"\n[TH-16 step 3a] control (flag-off full mesh): |V/V_q - 1|="
            f"{ratio_dev:.12f} (predicted 3); vertices outside quadrant="
            f"{r['n_outside']} (predicted > 0); owned exterior facets on x=0 "
            f"{r[BIRDCAGE_CUT_X0_TAG]['n_geo']}, on y=0 "
            f"{r[BIRDCAGE_CUT_Y0_TAG]['n_geo']} (predicted 0); "
            f"asserted={'yes' if CONTROL else 'no (W1: predicted, printed)'}",
            flush=True,
        )
    if CONTROL:
        assert abs(ratio_dev - 3.0) < EXACT_IDENTITY_RTOL * 3.0
        assert r["n_outside"] > 0
        assert n_on_planes == 0
