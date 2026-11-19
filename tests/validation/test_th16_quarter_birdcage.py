"""`TH-16` step 3c — the quarter birdcage's three solves equal its mirrored full 4x4.

Step 3a's quadrant cutter (``birdcage_port_domain(symmetry_quadrant=True)``)
keeps x >= 0, y >= 0 of the `GEO-19` gapped birdcage: 36 034 cells at ``-n 2``,
cut groups 501 (x = 0) and 502 (y = 0), the outer box on tag 1.  Every port
sheet lies *in* a cut plane (P1, P3 in y = 0; P2, P4 in x = 0; the drive ẑ is
tangential), so the quarter keeps P1q (the y = 0 face of cell tag 201) and P2q
(the x = 0 face of cell tag 102) as exterior sheets.

**Symmetry map** (class notation (x = 0 cut, y = 0 cut), + = PMC, − = PEC; a
sheet in a PEC plane is shorted, i.e. dead):

* (+,+): P1q, P2q live, 2x2;
* (−,+): P1q only; (+,−): P2q only;
* (−,−): no live port, not solved.

Quarter port: ``V_q = V``, ``I_q = I/2``, so ``Z_p`` doubles (100 Ω) and the
S reference is ``z0_q = 2 z0 = 100 Ω``; the quarter S entries are then exactly
the full mode ratios.  The sheet keeps its full width ``w`` (it lies in the
cut plane; ``halve_sheet_spec(…, "pmc")`` would halve it and is **not** used).

**Anchor (re-anchored 2026-09-30 03:00 review; §9 item 93).**  Step 3b read the
quarter against the independently meshed full record and missed its two-mesh
band (s 3.682e-2 / a 1.081e-2 / o 1.024e-2 vs 1e-3): a statement about two
meshes, not about the symmetry planes.  Here the quarter is mirrored across
y = 0, then x = 0, by vertex index (plane vertices reused), so the full mesh
is *exactly* D2-symmetric — not C4: P1q / P2q carry 24 / 23 facets — and its
4x4 ``S_f`` has five independent entries the three quarter solves must
reproduce to solver round-off.  The reconstruction ``R``:

``S11 = S33 = ½(S⁺⁺11 + S⁻⁺11)``, ``S13 = ½(S⁺⁺11 − S⁻⁺11)``,
``S22 = S44 = ½(S⁺⁺22 + S⁺⁻22)``, ``S24 = ½(S⁺⁺22 − S⁺⁻22)``,
all eight adjacent ``½ S⁺⁺12``.

Reconstructed on S, not Z (``z_matrix`` is normalised per column and does not
superpose).  The circulant built from the quarter's class means sits 1.05e-3
off ``R`` and is never compared entrywise.

Gate: premise (P1) 36 034 cells, (P2) owned cut facets 1 416 / 1 402 and
un-narrowed sheets 50 / 56, (P3) PEC classes constrain more DOFs than (+,+);
(M0) the mirrored mesh conforms — 144 136 cells, global vertices equal to the
builder's array, per-tag cell counts and volume x4, no exterior facet on
either plane, owned sheet facets 24 / 23 / 24 / 23 with the quarter's areas;
(I1) S⁺⁺ reciprocity ≤ 1e-12, ``S_f`` reciprocal ≤ 1e-12; (I2) passivity of
every S; (I3) ``S_f``'s own D2 pairs; **(A) max |S_f − R|_ij/|R_ij|** and
(B) R's class means against ``S_f``'s, all ≤ `LEG_D0_REPRODUCTION_BAND`
(1e-9); summed quarter solve-phase memory ΔM_q ≤ 0.425 GiB, read before the
mirrored mesh exists.  Negative control (wrong BC on x = 0: σ₋ := σ₊) against
``S_f``'s class means: predicted self miss 1.1447, printed in W1, asserted
> 0.5 under ``FEM_EM_TH16_STEP3C_CONTROL=1``.

Printed only: the two-mesh record (R's and ``S_f``'s class misses against
``S_10MHZ_EPS0_RECORD``), (A) per entry, ΔM_f.

Scope: one fixture, 10 MHz, degree 1, uncorrected width, ``-n 2``.  Not: any
two-mesh band, width matching, 64 / 128 MHz, degree 2, absolute S.
"""

from __future__ import annotations

import dataclasses
import gc
import os
import time

import dolfinx
import numpy as np
import pytest
import ufl
from dolfinx import default_scalar_type, fem
from mpi4py import MPI

from fem_em_solver.core import HomogeneousMaterial, TimeHarmonicProblem
from fem_em_solver.io.mesh import BIRDCAGE_CUT_X0_TAG, BIRDCAGE_CUT_Y0_TAG
from fem_em_solver.ports.definitions import PortDefinition
from fem_em_solver.ports.lumped import LumpedSheetPortSpec
from fem_em_solver.ports.sparameters import run_n_port_sparameter_sweep
from tests.complex_mode import complex_only
from tests.mesh.test_birdcage_port_sheets import _build
from tests.mesh.test_birdcage_symmetry_quadrant import EXACT_IDENTITY_RTOL
from tests.mesh.test_two_torus_port_sheet import _sheet_extents, _sheet_facet_count
from tests.validation.test_lossy_sphere_fullwave import SALINE_EPSILON_R, SALINE_SIGMA
from tests.validation.test_port21_feed_sensitivity import _class_means, _rss_report
from tests.validation.test_port_birdcage_four_port import (
    LEG_D0_REPRODUCTION_BAND,
    PASSIVITY_SIGMA_TOLERANCE,
    _circulant_classes,
)
from tests.validation.test_port_birdcage_lumped_column import (
    CONDUCTOR_CELL_TAG,
    PHANTOM_CELL_TAG,
    _narrowed_transverse,
    _sheet_bbox_centre,
)
from tests.validation.test_port_circuit_layer_field import S_10MHZ_EPS0_RECORD
from tests.validation.test_port_gap_voltage_impedance import (
    FREQUENCY_HZ,
    SIGMA_WIRE_S_PER_M,
)
from tests.validation.test_port_lumped_narrowed_sheet import GATED_WIDTH_FRACTION
from tests.validation.test_th16_symmetry_plane_driven import _distribute

CONTROL = os.environ.get("FEM_EM_TH16_STEP3C_CONTROL") == "1"

# Premise records, `TH-16` step 3a (`20260929T171313Z_TH-16-step3a-W1.log:2646–2650`).
QUARTER_CELL_RECORD = 36034
CUT_FACET_RECORD = {BIRDCAGE_CUT_X0_TAG: 1416, BIRDCAGE_CUT_Y0_TAG: 1402}
SHEET_FACET_RECORD = {1: 50, 2: 56}

# Full-mesh record: 50 Ω port at 50 Ω reference; the quarter doubles both.
FULL_Z0_OHM = 50.0
QUARTER_Z_OHM = 2.0 * FULL_Z0_OHM

# (I1): predicted ≈ 1e-14 (the full 4x4's own 1.04e-14).
RECIPROCITY_Q_BAND = 1.0e-12
# Memory pin: 0.310410 (3a cell ratio) x 1.369 GiB (full-mesh summed solve-phase
# ΔRSS, `20260927T110318Z_PORT-21-step3-W1.log:86,1886`: 1.870 − 0.501).
FULL_DELTA_M_GIB = 1.369
QUARTER_CELL_RATIO = 0.310410
MEMORY_BAR_GIB = 0.425
# Negative control (wrong BC on x = 0): predicted |o|/|s| of the quarter's own
# classes, 0.438294/0.382886 (`20260930T003844Z_TH-16-step3b-W1.log:942–943`).
CONTROL_PREDICTED = 1.1447
CONTROL_WINDOW = 0.001
CONTROL_FLOOR = 0.5

# (M0): 4 x the quarter; narrowed sheets P1q / P2q.
FULL_CELL_RECORD = 4 * QUARTER_CELL_RECORD
NARROWED_FACET_RECORD = {1: 24, 2: 23}
SHEET_AREA_RTOL = 1.0e-12
MIDPOINT_MATCH_TOL_M = 1.0e-10
# Full sheet k: (quarter sheet it copies, axis negated to reach it or None).
FULL_SHEETS = {1: (1, None), 2: (2, None), 3: (1, 0), 4: (2, 1)}
FULL_SHEET_TAG = {k: 210 + k for k in FULL_SHEETS}

PLANE_TOL = 1.0e-12
SHEET_TAG = {1: 211, 2: 212}
KEPT_HALF = {1: 201, 2: 102}
SHEET_PLANE_AXIS = {1: 1, 2: 0}  # P1 in y = 0, P2 in x = 0
SHEET_TRANSVERSE_AXIS = {1: 0, 2: 1}

CLASSES = {
    "(+,+)": ((1,), (1, 2)),
    "(-,+)": ((1, BIRDCAGE_CUT_X0_TAG), (1,)),
    "(+,-)": ((1, BIRDCAGE_CUT_Y0_TAG), (2,)),
}


def _say(comm, msg):
    if comm.rank == 0:
        print(f"[TH-16 step 3c] {msg}", flush=True)


def _owned_exterior_on_plane(msh, cells, axis, cell_tag):
    """Owned exterior facets on {x_axis = 0} whose cell carries ``cell_tag``."""
    tdim = msh.topology.dim
    fdim = tdim - 1
    n_owned = msh.topology.index_map(fdim).size_local
    ext = np.asarray(dolfinx.mesh.exterior_facet_indices(msh.topology), dtype=np.int32)
    ext = ext[ext < n_owned]
    if ext.size == 0:
        return ext
    nodes = dolfinx.cpp.mesh.entities_to_geometry(msh._cpp_object, fdim, ext, False)
    pts = msh.geometry.x[nodes]
    on = np.all(np.abs(pts[:, :, axis]) < PLANE_TOL, axis=1)
    imap = msh.topology.index_map(tdim)
    tag_of_cell = np.full(imap.size_local + imap.num_ghosts, -1, dtype=np.int64)
    tag_of_cell[np.asarray(cells.indices)] = np.asarray(cells.values)
    f2c = msh.topology.connectivity(fdim, tdim)
    ctag = np.array([tag_of_cell[f2c.links(int(f))[0]] for f in ext], dtype=np.int64)
    return np.sort(ext[on & (ctag == cell_tag)]).astype(np.int32)


def _tri_area_owned(msh, facets):
    if len(facets) == 0:
        return 0.0
    nodes = dolfinx.cpp.mesh.entities_to_geometry(
        msh._cpp_object, msh.topology.dim - 1, np.asarray(facets, dtype=np.int32), False
    )
    p = msh.geometry.x[nodes][:, :3, :]
    return float(0.5 * np.linalg.norm(np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]), axis=1).sum())


def _owned_count(comm, facet_tags, tag, n_owned):
    f = np.asarray(facet_tags.find(tag))
    return comm.allreduce(int(np.count_nonzero(f < n_owned)), op=MPI.SUM)


def _constrained_dofs(msh, facet_tags, pec_tags):
    """Rank-summed owned N1curl-1 DOFs on the union of ``pec_tags`` facets."""
    V = fem.functionspace(msh, ("N1curl", 1))
    found = [np.asarray(facet_tags.find(int(t)), dtype=np.int32) for t in pec_tags]
    facets = np.unique(np.concatenate(found)).astype(np.int32)
    dofs = fem.locate_dofs_topological(V, msh.topology.dim - 1, facets)
    n_owned = V.dofmap.index_map.size_local * V.dofmap.index_map_bs
    return msh.comm.allreduce(int(np.count_nonzero(dofs < n_owned)), op=MPI.SUM)


def _facet_midpoints(msh, facets):
    if len(facets) == 0:
        return np.empty((0, 3), dtype=np.float64)
    nodes = dolfinx.cpp.mesh.entities_to_geometry(
        msh._cpp_object, msh.topology.dim - 1, np.asarray(facets, dtype=np.int32), False
    )
    return np.asarray(msh.geometry.x[nodes][:, :3, :].mean(axis=1), dtype=np.float64)


def _owned_cell_tags(msh, cells):
    imap = msh.topology.index_map(msh.topology.dim)
    tag_of_cell = np.full(imap.size_local + imap.num_ghosts, -1, dtype=np.int32)
    tag_of_cell[np.asarray(cells.indices)] = np.asarray(cells.values)
    return tag_of_cell[: imap.size_local]


def _tag_counts(comm, owned_tags):
    """Rank-summed owned cell count per tag."""
    vals, counts = np.unique(np.concatenate(comm.allgather(owned_tags)), return_counts=True)
    return {int(v): int(c) for v, c in zip(vals, counts)}


def _volume(msh):
    one = fem.Constant(msh, default_scalar_type(1.0))
    local = fem.assemble_scalar(fem.form(one * ufl.dx(domain=msh)))
    return float(np.real(msh.comm.allreduce(local, op=MPI.SUM)))


def _sheet_integral_area(msh, tags, tag, interior):
    one = fem.Constant(msh, default_scalar_type(1.0))
    meas = ufl.Measure("dS" if interior else "ds", domain=msh, subdomain_data=tags)
    local = fem.assemble_scalar(fem.form(one * meas(int(tag))))
    return float(np.real(msh.comm.allreduce(local, op=MPI.SUM)))


def _gather_owned(msh, cells):
    """Rank 0: (x, cells, tags) of the whole mesh, cells in global geometry indices.

    Owned geometry nodes are globally numbered contiguously in rank order, so
    the rank-ordered concatenation of the owned coordinates *is* the global
    numbering.
    """
    comm = msh.comm
    gmap = msh.geometry.index_map()
    n_cells = msh.topology.index_map(msh.topology.dim).size_local
    dm = np.asarray(msh.geometry.dofmap, dtype=np.int32)[:n_cells]
    gcells = np.asarray(gmap.local_to_global(dm.reshape(-1)), dtype=np.int64).reshape(-1, 4)
    own_tags = _owned_cell_tags(msh, cells)
    xs = comm.gather(np.array(msh.geometry.x[: gmap.size_local, :3], dtype=np.float64), root=0)
    cs = comm.gather(gcells, root=0)
    ts = comm.gather(own_tags, root=0)
    if comm.rank != 0:
        return None
    return np.vstack(xs), np.vstack(cs), np.concatenate(ts)


def _mirror_by_index(x, cells, tags, axis):
    """Mirror across {x_axis = 0}; vertices on the plane are reused by index."""
    on = np.abs(x[:, axis]) < PLANE_TOL
    n = x.shape[0]
    new_index = np.arange(n, dtype=np.int64)
    new_index[~on] = n + np.arange(int(np.count_nonzero(~on)), dtype=np.int64)
    xm = x[~on].copy()
    xm[:, axis] = -xm[:, axis]
    return (
        np.vstack([x, xm]),
        np.vstack([cells, new_index[cells][:, [1, 0, 2, 3]]]),
        np.concatenate([tags, tags]),
    )


def _full_matrix_from_quarter(spp, s_mp, s_pm):
    """R: the full 4x4 the three quarter solves imply under exact D2 symmetry."""
    r = np.empty((4, 4), dtype=np.complex128)
    r[:, :] = 0.5 * spp[0, 1]
    r[0, 0] = r[2, 2] = 0.5 * (spp[0, 0] + s_mp)
    r[0, 2] = r[2, 0] = 0.5 * (spp[0, 0] - s_mp)
    r[1, 1] = r[3, 3] = 0.5 * (spp[1, 1] + s_pm)
    r[1, 3] = r[3, 1] = 0.5 * (spp[1, 1] - s_pm)
    return r


def _record_classes():
    means = _class_means(S_10MHZ_EPS0_RECORD)
    return {"s": means["self"], "a": means["adjacent"], "o": means["opposite"]}


@pytest.fixture(scope="module")
def run():
    comm = MPI.COMM_WORLD
    out = {"rss": {"start": _rss_report(comm, "start (baseline B)")}}

    msh, cells, facets, diag, t_mesh = _build(True, symmetry_quadrant=True)
    tdim = msh.topology.dim
    fdim = tdim - 1
    # Hoisted on every rank before any facet form (known-issues 9).
    msh.topology.create_connectivity(fdim, tdim)
    msh.topology.create_entity_permutations()
    n_owned_f = msh.topology.index_map(fdim).size_local
    out["cells"] = int(msh.topology.index_map(tdim).size_global)
    out["cut_facets"] = {t: _owned_count(comm, facets, t, n_owned_f) for t in CUT_FACET_RECORD}

    # Sheet meshtags, separate from the problem's facet tags.
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
        ext = _sheet_extents(msh, tags_s, tag, comm)
        ax = SHEET_TRANSVERSE_AXIS[p]
        w_bbox = float(ext[ax])
        centre = _sheet_bbox_centre(msh, tags_s, tag, comm)
        tags_s = _narrowed_transverse(
            msh, tags_s, tag, GATED_WIDTH_FRACTION, float(centre[ax]), ax, 0.5 * w_bbox
        )
    for p in (1, 2):
        tag = SHEET_TAG[p]
        n = _sheet_facet_count(msh, tags_s, tag, comm)
        f_owned = np.asarray(tags_s.find(tag))
        f_owned = f_owned[f_owned < n_owned_f]
        area = comm.allreduce(_tri_area_owned(msh, f_owned), op=MPI.SUM)
        ext = _sheet_extents(msh, tags_s, tag, comm)
        h = float(ext[2])
        mids = np.vstack(comm.allgather(_facet_midpoints(msh, f_owned)))
        sheets[p] = {
            "facets": int(n),
            "area": area,
            "h": h,
            "w": area / h,
            "mid": mids,
            "ds_area": _sheet_integral_area(msh, tags_s, tag, False),
            "w_bbox": float(ext[SHEET_TRANSVERSE_AXIS[p]]),
            "flat": float(ext[SHEET_PLANE_AXIS[p]]),
        }
    out["sheets"] = sheets

    base = TimeHarmonicProblem(
        mesh=msh,
        frequency_hz=FREQUENCY_HZ,
        material=HomogeneousMaterial(sigma=0.0, epsilon_r=1.0, mu_r=1.0),
        cell_tags=cells,
        material_map={
            CONDUCTOR_CELL_TAG: HomogeneousMaterial(
                sigma=SIGMA_WIRE_S_PER_M, epsilon_r=1.0, mu_r=1.0
            ),
            PHANTOM_CELL_TAG: HomogeneousMaterial(
                sigma=SALINE_SIGMA, epsilon_r=SALINE_EPSILON_R, mu_r=1.0
            ),
        },
        boundary_condition="pec_zero_tangential_a",
        facet_tags=facets,
    )

    def _port(p):
        pid = f"P{p}q"
        d = PortDefinition(
            port_id=pid,
            positive_tag=SHEET_TAG[p],
            negative_tag=CONDUCTOR_CELL_TAG,
            orientation="leg_gap_axial_plus_z",
            z0_ohm=QUARTER_Z_OHM,
        )
        s = LumpedSheetPortSpec(
            port_id=pid,
            facet_tag=SHEET_TAG[p],
            port_impedance_ohm=QUARTER_Z_OHM,
            gap_height_m=sheets[p]["h"],
            sheet_width_m=sheets[p]["w"],
            drive_direction=(0.0, 0.0, 1.0),
            drive_voltage_v=1.0 + 0.0j,
            interior=False,
        )
        return d, s

    out["dofs"] = {}
    out["S"] = {}
    out["t"] = {}
    for name, (pec_tags, live) in CLASSES.items():
        out["dofs"][name] = _constrained_dofs(msh, facets, pec_tags)
        problem = dataclasses.replace(base, pec_facet_tags=tuple(pec_tags))
        defs, specs = zip(*[_port(p) for p in live])
        comm.Barrier()
        t0 = time.perf_counter()
        result = run_n_port_sparameter_sweep(
            problem,
            list(defs),
            lumped_sheet_ports=list(specs),
            lumped_sheet_facet_tags=tags_s,
            degree=1,
        )
        comm.Barrier()
        out["t"][name] = time.perf_counter() - t0
        assert tuple(result.port_ids) == tuple(d.port_id for d in defs)
        out["S"][name] = np.asarray(result.s_matrix, dtype=np.complex128).copy()
        del result, problem
        gc.collect()
        out["rss"][name] = _rss_report(comm, f"{name} sweep")

    out["t_mesh"] = t_mesh
    out["diag"] = diag

    spp = out["S"]["(+,+)"]
    s_mp = complex(out["S"]["(-,+)"][0, 0])
    s_pm = complex(out["S"]["(+,-)"][0, 0])
    sig_p = 0.5 * (spp[0, 0] + spp[1, 1])
    sig_m = 0.5 * (s_mp + s_pm)
    alpha = 0.25 * (spp[0, 1] + spp[1, 0])
    q = {"s": 0.5 * (sig_p + sig_m), "a": alpha, "o": 0.5 * (sig_p - sig_m)}
    rec = _record_classes()
    miss = {k: abs(q[k] - rec[k]) / abs(rec[k]) for k in q}

    circ = np.empty((4, 4), dtype=np.complex128)
    for r in range(4):
        for c in range(4):
            step = (r - c) % 4
            circ[r, c] = q["s"] if step == 0 else (q["o"] if step == 2 else q["a"])
    rec_classes = _circulant_classes(S_10MHZ_EPS0_RECORD)
    rec_circ_floor = float(
        np.linalg.norm(
            S_10MHZ_EPS0_RECORD
            - np.array(
                [
                    [
                        (
                            rec["s"]
                            if (r - c) % 4 == 0
                            else (rec["o"] if (r - c) % 4 == 2 else rec["a"])
                        )
                        for c in range(4)
                    ]
                    for r in range(4)
                ]
            )
        )
        / np.linalg.norm(S_10MHZ_EPS0_RECORD)
    )
    circ_rel = float(
        np.linalg.norm(circ - S_10MHZ_EPS0_RECORD) / np.linalg.norm(S_10MHZ_EPS0_RECORD)
    )
    peak = max(v["peak_sum_gib"] for v in out["rss"].values())
    # ΔM_q is fixed here, before the mirrored mesh exists (ru_maxrss is a
    # high-water mark); the full mesh's RSS goes in the `full` fixture's dict.
    out.update(
        sig_p=sig_p,
        sig_m=sig_m,
        alpha=alpha,
        q=q,
        rec=rec,
        miss=miss,
        circ_rel=circ_rel,
        s_mp=s_mp,
        s_pm=s_pm,
        R=_full_matrix_from_quarter(spp, s_mp, s_pm),
        rec_circ_floor=rec_circ_floor,
        recip=float(abs(spp[0, 1] - spp[1, 0]) / abs(spp[1, 0])),
        sigma_pp=float(np.linalg.svd(spp, compute_uv=False).max()),
        abs_mp=abs(s_mp),
        abs_pm=abs(s_pm),
        mirror_pp=abs(spp[0, 0] - spp[1, 1]) / abs(0.5 * (spp[0, 0] + spp[1, 1])),
        mirror_m=abs(s_mp - s_pm) / abs(0.5 * (s_mp + s_pm)),
        delta_m=peak - out["rss"]["start"]["peak_sum_gib"],
    )

    if comm.rank == 0:
        print(
            f"\n[TH-16 step 3c] quarter cells={out['cells']} (record {QUARTER_CELL_RECORD}); "
            f"mesh {diag['mesh_wall_time_s']:.2f} s, rung {t_mesh:.2f} s; -n {comm.size}"
            f"\n[TH-16 step 3c] owned cut facets 501={out['cut_facets'][BIRDCAGE_CUT_X0_TAG]} "
            f"502={out['cut_facets'][BIRDCAGE_CUT_Y0_TAG]} (records 1416 / 1402); "
            f"un-narrowed sheets P1q={unnarrowed[1]} P2q={unnarrowed[2]} (records 50 / 56)",
            flush=True,
        )
        for p in (1, 2):
            s = sheets[p]
            print(
                f"[TH-16 step 3c] P{p}q narrowed (f={GATED_WIDTH_FRACTION}): {s['facets']} facets "
                f"area={s['area']:.9e} m^2 h={s['h']:.9e} m w=A/h={s['w']:.9e} m "
                f"w_bbox={s['w_bbox']:.9e} m flat={s['flat']:.3e} m "
                f"(full mesh: 26 facets, w=7.294123600e-03 m)",
                flush=True,
            )
        for name in CLASSES:
            S = out["S"][name]
            print(
                f"[TH-16 step 3c] {name}: constrained DOFs={out['dofs'][name]} "
                f"sweep {out['t'][name]:.2f} s; S=" + np.array2string(S, precision=9),
                flush=True,
            )
        print(
            f"[TH-16 step 3c] sigma+={sig_p:.9f} sigma-={sig_m:.9f} alpha={alpha:.9f}"
            f"\n[TH-16 step 3c] classes: s_q={q['s']:.9f} a_q={q['a']:.9f} o_q={q['o']:.9f}"
            f"\n[TH-16 step 3c] record:  s={rec['s']:.9f} a={rec['a']:.9f} o={rec['o']:.9f}"
            f"\n[TH-16 step 3c] two-mesh record, printed only: quarter class misses s {miss['s']:.3e}  "
            f"a {miss['a']:.3e}  o {miss['o']:.3e}  (predicted 3.682e-2 / 1.081e-2 / 1.024e-2)"
            f"\n[TH-16 step 3c] (I1) |S12-S21|/|S21|={out['recip']:.3e} (band {RECIPROCITY_Q_BAND:.0e})"
            f"\n[TH-16 step 3c] (I2) sigma_max(S++)={out['sigma_pp']:.9f} |S-+11|={out['abs_mp']:.9f} "
            f"|S+-22|={out['abs_pm']:.9f} (ceiling 1 + {PASSIVITY_SIGMA_TOLERANCE:.0e})"
            f"\n[TH-16 step 3c] mirror spreads: |S++11-S++22|/mean={out['mirror_pp']:.3e} "
            f"|S-+11-S+-22|/mean={out['mirror_m']:.3e} (predicted <~ 6e-4)"
            f"\n[TH-16 step 3c] mode eigenvalues: s+o+2a={q['s'] + q['o'] + 2 * q['a']:.6f} "
            f"s+o-2a={q['s'] + q['o'] - 2 * q['a']:.6f} s-o={q['s'] - q['o']:.6f} "
            "(record +0.999904-0.013299j, -0.868303+0.153771j, -0.808126+0.213282j)"
            f"\n[TH-16 step 3c] ||circ - record||_F/||record||_F={circ_rel:.3e} "
            f"(record's own circulant floor {rec_circ_floor:.3e}); record class spreads "
            + " ".join(
                f"{k}={max(abs(v - np.mean(v))) / abs(np.mean(v)):.3e}"
                for k, v in rec_classes.items()
            )
            + f"\n[TH-16 step 3c] MEMORY dM_q={out['delta_m']:.3f} GiB (bar {MEMORY_BAR_GIB} = "
            f"{QUARTER_CELL_RATIO} x {FULL_DELTA_M_GIB}; literal 1/4 = {0.25 * FULL_DELTA_M_GIB:.3f}; "
            f"predicted 0.29-0.43)",
            flush=True,
        )

    # For the mirror build: the whole quarter on rank 0, and its own counts.
    out["tag_counts"] = _tag_counts(comm, _owned_cell_tags(msh, cells))
    out["volume"] = _volume(msh)
    out["gathered"] = _gather_owned(msh, cells)
    return out


@pytest.fixture(scope="module")
def full(run):
    """The quarter mirrored across y = 0, then x = 0: four interior sheets, one sweep."""
    comm = MPI.COMM_WORLD
    out = {"rss": {"start": _rss_report(comm, "before the mirror build")}}
    comm.Barrier()
    t0 = time.perf_counter()
    x = c = t = None
    if comm.rank == 0:
        x, c, t = run["gathered"]
        for axis in (1, 0):
            x, c, t = _mirror_by_index(x, c, t, axis)
    out["builder_vertices"] = comm.bcast(None if x is None else int(x.shape[0]), root=0)
    msh = _distribute(x, c, comm)
    tdim = msh.topology.dim
    fdim = tdim - 1

    # Cell tags the way the installed `dolfinx.io.gmsh.model_to_mesh` does it
    # (cells only: interior facets hang at -n 2, known-issues 9).  Entities are
    # in *input* vertex indices; non-root ranks pass empty arrays.
    if comm.rank == 0:
        ents = np.ascontiguousarray(c, dtype=np.int64)
        vals = np.ascontiguousarray(t, dtype=np.int32)
    else:
        ents = np.empty((0, 4), dtype=np.int64)
        vals = np.empty(0, dtype=np.int32)
    local_entities, local_values = dolfinx.io.distribute_entity_data(msh, tdim, ents, vals)
    msh.topology.create_connectivity(tdim, tdim)
    cells = dolfinx.mesh.meshtags_from_entities(
        msh,
        tdim,
        dolfinx.graph.adjacencylist(local_entities),
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
            n_on = int(
                np.count_nonzero(
                    np.all(np.abs(msh.geometry.x[nodes][:, :, axis]) < PLANE_TOL, axis=1)
                )
            )
        out["exterior_on_plane"][axis] = comm.allreduce(n_on, op=MPI.SUM)

    # Sheets 211-214: the quarter's narrowed facets (now interior) or their
    # mirror images, matched by midpoint.  No second narrowing.
    plane = {}
    for axis in (0, 1):
        f = dolfinx.mesh.locate_entities(
            msh, fdim, lambda p, a=axis: np.abs(p[a]) < PLANE_TOL
        ).astype(np.int32)
        f = f[f < n_owned_f]
        plane[axis] = (f, _facet_midpoints(msh, f))
    idx, val, worst = [], [], 0.0
    for k, (p, flip) in FULL_SHEETS.items():
        target = run["sheets"][p]["mid"].copy()
        if flip is not None:
            target[:, flip] = -target[:, flip]
        f, mid = plane[SHEET_PLANE_AXIS[p]]
        if f.size:
            d = np.linalg.norm(mid[:, None, :] - target[None, :, :], axis=2).min(axis=1)
            hit = d <= MIDPOINT_MATCH_TOL_M
            if hit.any():
                worst = max(worst, float(d[hit].max()))
            idx.append(f[hit])
            val.append(np.full(int(np.count_nonzero(hit)), FULL_SHEET_TAG[k], dtype=np.int32))
    idx = np.concatenate(idx).astype(np.int32) if idx else np.empty(0, dtype=np.int32)
    val = np.concatenate(val).astype(np.int32) if val else np.empty(0, dtype=np.int32)
    assert np.unique(idx).size == idx.size
    order = np.argsort(idx, kind="stable")
    tags_s = dolfinx.mesh.meshtags(msh, fdim, idx[order], val[order])
    out["midpoint_worst"] = comm.allreduce(worst, op=MPI.MAX)
    out["sheet_facets"] = {
        k: _owned_count(comm, tags_s, FULL_SHEET_TAG[k], n_owned_f) for k in FULL_SHEETS
    }
    out["sheet_area"] = {
        k: _sheet_integral_area(msh, tags_s, FULL_SHEET_TAG[k], True) for k in FULL_SHEETS
    }
    out["rss"]["built"] = _rss_report(comm, "mirror build + tags")

    # Problem as `test_port_birdcage_four_port.py:406-424`: no facet_tags, so
    # every exterior facet is PEC.
    problem = TimeHarmonicProblem(
        mesh=msh,
        frequency_hz=FREQUENCY_HZ,
        material=HomogeneousMaterial(sigma=0.0, epsilon_r=1.0, mu_r=1.0),
        cell_tags=cells,
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
                gap_height_m=run["sheets"][p]["h"],
                sheet_width_m=run["sheets"][p]["w"],
                drive_direction=(0.0, 0.0, 1.0),
                drive_voltage_v=1.0 + 0.0j,
                interior=True,
            )
        )
    comm.Barrier()
    t0 = time.perf_counter()
    result = run_n_port_sparameter_sweep(
        problem, defs, lumped_sheet_ports=specs, lumped_sheet_facet_tags=tags_s, degree=1
    )
    comm.Barrier()
    out["t_sweep"] = time.perf_counter() - t0
    assert tuple(result.port_ids) == tuple(d.port_id for d in defs)
    sf = np.asarray(result.s_matrix, dtype=np.complex128).copy()
    del result, problem
    gc.collect()
    out["rss"]["sweep"] = _rss_report(comm, "full four-drive sweep")

    R = run["R"]
    entry = np.abs(sf - R) / np.abs(R)
    off = ~np.eye(4, dtype=bool)
    recip = float(np.max(np.abs(sf - sf.T)[off] / np.abs(sf.T)[off]))
    adj = np.array([sf[i, j] for i in range(4) for j in range(4) if (i - j) % 2 == 1])
    d2 = {
        "S11/S33": float(abs(sf[0, 0] - sf[2, 2]) / abs(0.5 * (sf[0, 0] + sf[2, 2]))),
        "S22/S44": float(abs(sf[1, 1] - sf[3, 3]) / abs(0.5 * (sf[1, 1] + sf[3, 3]))),
        "adjacent": float(np.max(np.abs(adj - adj.mean())) / abs(adj.mean())),
    }
    cm_f, cm_r = _class_means(sf), _class_means(R)
    class_miss = {k: float(abs(cm_r[k] - cm_f[k]) / abs(cm_f[k])) for k in cm_f}
    rec = _class_means(S_10MHZ_EPS0_RECORD)
    # Negative control: PMC for PEC on x = 0 makes the (−,+) solve the (+,+)
    # one, so s_w = σ₊ and o_w = 0.
    ctrl_self = float(abs(run["sig_p"] - cm_f["self"]) / abs(cm_f["self"]))
    ctrl_opp = float(abs(0.0 - cm_f["opposite"]) / abs(cm_f["opposite"]))
    out.update(
        S=sf,
        entry=entry,
        A=float(entry.max()),
        recip=recip,
        d2=d2,
        sigma=float(np.linalg.svd(sf, compute_uv=False).max()),
        class_miss=class_miss,
        ctrl_self=ctrl_self,
        ctrl_opp=ctrl_opp,
        delta_m=out["rss"]["sweep"]["peak_sum_gib"] - out["rss"]["start"]["peak_sum_gib"],
    )

    if comm.rank == 0:
        qc, fc = run["tag_counts"], out["tag_counts"]
        print(
            f"\n[TH-16 step 3c] (M0) mirrored full: cells={out['cells']} (4 x {run['cells']} = "
            f"{FULL_CELL_RECORD}); vertices={out['vertices']} (builder {out['builder_vertices']}); "
            f"build {out['t_build']:.2f} s"
            f"\n[TH-16 step 3c] (M0) per-tag cells full/quarter: "
            + " ".join(f"{k}:{fc.get(k)}/{qc[k]}" for k in sorted(qc))
            + f"\n[TH-16 step 3c] (M0) volume full={out['volume']:.12e} quarter={run['volume']:.12e} "
            f"|V_f/(4 V_q) - 1|={abs(out['volume'] / (4.0 * run['volume']) - 1.0):.3e} "
            f"(rtol {EXACT_IDENTITY_RTOL:.0e})"
            f"\n[TH-16 step 3c] (M0) owned exterior facets on x=0: {out['exterior_on_plane'][0]}, "
            f"on y=0: {out['exterior_on_plane'][1]} (both 0)"
            f"\n[TH-16 step 3c] (M0) sheets 211-214 owned facets: "
            + " / ".join(str(out["sheet_facets"][k]) for k in FULL_SHEETS)
            + f" (24 / 23 / 24 / 23); worst midpoint match {out['midpoint_worst']:.3e} m; dS areas rel. "
            + " ".join(
                f"{abs(out['sheet_area'][k] / run['sheets'][p]['area'] - 1.0):.2e}"
                for k, (p, _f) in FULL_SHEETS.items()
            )
            + f" (rtol {SHEET_AREA_RTOL:.0e}); quarter ds/geometric "
            + " ".join(
                f"{abs(run['sheets'][p]['ds_area'] / run['sheets'][p]['area'] - 1.0):.2e}"
                for p in (1, 2)
            )
            + f"\n[TH-16 step 3c] full sweep {out['t_sweep']:.2f} s; S_f=\n"
            + np.array2string(sf, precision=9)
            + "\n[TH-16 step 3c] R (from the three quarter solves)=\n"
            + np.array2string(R, precision=9)
            + "\n[TH-16 step 3c] predicted: S11 -0.356438+0.138737j, S22 -0.357219+0.138930j, "
            "S13 +0.432256-0.070029j, S24 +0.433036-0.070238j, adjacent +0.462043-0.040980j"
            + "\n[TH-16 step 3c] (A) |S_f - R|_ij/|R_ij| per entry=\n"
            + np.array2string(entry, precision=3)
            + f"\n[TH-16 step 3c] (A) max={out['A']:.3e} (band {LEG_D0_REPRODUCTION_BAND:.0e}; "
            "predicted 1e-13 - 1e-10)"
            "\n[TH-16 step 3c] (B) class means R vs S_f: "
            + " ".join(f"{k} {v:.3e}" for k, v in class_miss.items())
            + f" (band {LEG_D0_REPRODUCTION_BAND:.0e})"
            f"\n[TH-16 step 3c] (I1) S_f reciprocity max={recip:.3e} (band {RECIPROCITY_Q_BAND:.0e}); "
            f"(I2) sigma_max(S_f)={out['sigma']:.9f}"
            "\n[TH-16 step 3c] (I3) S_f D2 pairs: "
            + " ".join(f"{k} {v:.3e}" for k, v in d2.items())
            + f" (band {LEG_D0_REPRODUCTION_BAND:.0e})"
            "\n[TH-16 step 3c] two-mesh record, printed only: class misses against "
            "S_10MHZ_EPS0_RECORD: R "
            + " ".join(f"{k} {abs(cm_r[k] - rec[k]) / abs(rec[k]):.3e}" for k in rec)
            + "; S_f "
            + " ".join(f"{k} {abs(cm_f[k] - rec[k]) / abs(rec[k]):.3e}" for k in rec)
            + " (predicted unchanged: 3.682e-2 / 1.081e-2 / 1.024e-2)"
            f"\n[TH-16 step 3c] MEMORY dM_f={out['delta_m']:.3f} GiB (printed only; predicted ~1.8), "
            f"peak sum {out['rss']['sweep']['peak_sum_gib']:.3f} GiB (predicted ~2.4)"
            "\n[TH-16 step 3c] CONTROL (PMC for PEC on x = 0) against S_f's class means: "
            f"self miss={ctrl_self:.4f} (predicted {CONTROL_PREDICTED} +/- {CONTROL_WINDOW}), "
            f"opposite miss={ctrl_opp:.4f} (predicted 1); "
            f"asserted={'yes' if CONTROL else 'no (W1: predicted, printed)'}",
            flush=True,
        )
    return out


@complex_only
def test_premise(run):
    assert run["cells"] == QUARTER_CELL_RECORD
    assert run["cut_facets"] == CUT_FACET_RECORD
    assert run["sheets_unnarrowed"] == SHEET_FACET_RECORD
    assert run["dofs"]["(-,+)"] > run["dofs"]["(+,+)"]
    assert run["dofs"]["(+,-)"] > run["dofs"]["(+,+)"]


@complex_only
def test_one_mesh_identities(run):
    assert run["recip"] <= RECIPROCITY_Q_BAND
    ceiling = 1.0 + PASSIVITY_SIGMA_TOLERANCE
    assert run["sigma_pp"] <= ceiling
    assert run["abs_mp"] <= ceiling
    assert run["abs_pm"] <= ceiling


@complex_only
def test_mirror_build_conforms(run, full):
    """(M0)."""
    assert full["cells"] == FULL_CELL_RECORD
    assert full["vertices"] == full["builder_vertices"]
    assert full["tag_counts"] == {k: 4 * v for k, v in run["tag_counts"].items()}
    assert abs(full["volume"] / (4.0 * run["volume"]) - 1.0) <= EXACT_IDENTITY_RTOL
    assert full["exterior_on_plane"] == {0: 0, 1: 0}
    assert {p: run["sheets"][p]["facets"] for p in (1, 2)} == NARROWED_FACET_RECORD
    for k, (p, _flip) in FULL_SHEETS.items():
        assert full["sheet_facets"][k] == NARROWED_FACET_RECORD[p], k
        assert abs(full["sheet_area"][k] / run["sheets"][p]["area"] - 1.0) <= SHEET_AREA_RTOL, k


@complex_only
def test_full_identities(full):
    """(I1)-(I3) on the mirrored full mesh's own 4x4."""
    assert full["recip"] <= RECIPROCITY_Q_BAND
    assert full["sigma"] <= 1.0 + PASSIVITY_SIGMA_TOLERANCE
    for k, v in full["d2"].items():
        assert v <= LEG_D0_REPRODUCTION_BAND, (k, v)


@complex_only
def test_quarter_reproduces_mirrored_full(full):
    """(A) per entry and (B) per class."""
    assert full["A"] <= LEG_D0_REPRODUCTION_BAND
    for k, m in full["class_miss"].items():
        assert m <= LEG_D0_REPRODUCTION_BAND, (k, m)


@complex_only
def test_memory_bar(run):
    assert run["delta_m"] <= MEMORY_BAR_GIB


@complex_only
def test_negative_control_wrong_bc(full):
    if CONTROL:
        assert full["ctrl_self"] > CONTROL_FLOOR
