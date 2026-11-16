"""`POST-8` step 3b-ii-t1 — the time-harmonic test sites under the owner rule.

PROJECT_PLAN §9 item 536; brief ``docs/planning/chunks/POST-8.md``
(``### §9 item 536`` and its 2026-10-18 10:30 review ruling).

Fixture: ``_run_at_degree(d)`` for d = 1, 2 — `TH-12` step 1's coarse-rung
lossy saline sphere (5 866 cells), whose ``e_real`` / ``e_imag`` are the DG
splits of the N1curl phasor sampled at `TH-10`'s 25 interior probes
(``test_lossy_sphere_degree2.py`` / ``test_lossy_sphere_fullwave.py`` now pass
``owner="lowest_original_cell"`` there).

Asserted: (S) the rule is a bitwise no-op at every single-incidence probe and
matches ``.eval`` in the owner cell at any multi-incidence probe; (O) owner
original indices equal the ``-n 2`` npz at ``-n 4`` / ``-n 1``; (W) every
sampled E within ``WIDTH_RTOL`` of the ``-n 2`` npz; (N0) at the vertex V
(lexicographically smallest vertex of probe 0's degree-1 owner cell)
``n_incident >= 4`` and ``max/min |E|`` over the incident cells > 1 + 1e-6.
Printed: the single-incidence count, (N1) the legacy cell at V.
"""

from pathlib import Path

import numpy as np
import pytest
from mpi4py import MPI

from dolfinx import geometry

from fem_em_solver.post.evaluation import (
    evaluate_vector_field_parallel,
    point_owner_cells,
)
from tests.complex_mode import complex_only
from tests.post.test_point_owner_rule import (
    CLOSED_FORM_TOL,
    _bitwise_equal_across_ranks,
)
from tests.post.test_post8_step3b_point_sar_owner import _dg0_index
from tests.validation.test_lossy_sphere_degree2 import _run_at_degree
from tests.validation.test_lossy_sphere_fullwave import _interior_probe_points

NPZ_DIR = Path("/workspace/logs/post8")
NPZ = NPZ_DIR / "step3b-ii-t1-n2.npz"  # never step3a-n2.npz / step3b-n2.npz
RULE = "lowest_original_cell"
#: (N0) asserted floors (rule (e)); predicted max/min > 1.01.
N0_MIN_INCIDENT = 4
N0_SPREAD_FLOOR = 1.0e-6

pytestmark = complex_only


@pytest.fixture(scope="module")
def rows():
    return {d: _run_at_degree(d) for d in (1, 2)}


def _mesh_of(row):
    return row["fields"].e_real.function_space.mesh


def _n_incident(msh, points, comm):
    """Number of *owned* cells containing each point, summed over ranks."""
    tdim = msh.topology.dim
    n_owned = msh.topology.index_map(tdim).size_local
    tree = geometry.bb_tree(msh, tdim)
    hits = geometry.compute_collisions_points(tree, points)
    coll = geometry.compute_colliding_cells(msh, hits, points)
    local = np.zeros(points.shape[0], dtype=np.int64)
    for i in range(points.shape[0]):
        links = np.asarray(coll.links(i))
        local[i] = int(np.sum(links < n_owned))
    tot = np.empty_like(local)
    comm.Allreduce(local, tot, op=MPI.SUM)
    return tot


def _owned_incident_cells(msh, point):
    tdim = msh.topology.dim
    n_owned = msh.topology.index_map(tdim).size_local
    p = np.asarray(point, dtype=np.float64).reshape(1, 3)
    tree = geometry.bb_tree(msh, tdim)
    hits = geometry.compute_collisions_points(tree, p)
    coll = geometry.compute_colliding_cells(msh, hits, p)
    links = np.asarray(coll.links(0))
    return links[links < n_owned].astype(np.int32)


def _e_in_cell(row, point, local_cell, comm):
    """Complex E at ``point`` in ``local_cell`` (on the rank where it is >= 0).

    The owning rank evaluates; the result is shared by buffer ``Allreduce``
    SUM over real/imag parts — independent of the rule's gather path.
    """
    f = row["fields"]
    buf = np.zeros(6, dtype=np.float64)
    if local_cell >= 0:
        p = np.asarray(point, dtype=np.float64).reshape(1, 3)
        c = np.array([int(local_cell)], dtype=np.int32)
        buf[:3] = np.real(np.asarray(f.e_real.eval(p, c)).reshape(-1))
        buf[3:] = np.real(np.asarray(f.e_imag.eval(p, c)).reshape(-1))
    tot = np.zeros_like(buf)
    comm.Allreduce(buf, tot, op=MPI.SUM)
    return tot[:3] + 1j * tot[3:]


def _sample(row, points, comm, owner):
    f = row["fields"]
    rv, vr = evaluate_vector_field_parallel(f.e_real, points, comm, owner=owner)
    iv, vi = evaluate_vector_field_parallel(f.e_imag, points, comm, owner=owner)
    assert bool(np.all(vr & vi))
    return np.asarray(rv), np.asarray(iv)


def _vertex_v(row, comm):
    """(V): lexicographically smallest vertex of probe 0's owner cell."""
    msh = _mesh_of(row)
    p0 = _interior_probe_points()[:1]
    _, owner_local = point_owner_cells(msh, p0, comm)
    buf = np.zeros(3, dtype=np.float64)
    if owner_local[0] >= 0:
        verts = msh.geometry.x[msh.geometry.dofmap[int(owner_local[0])]][:, :3]
        order = np.lexsort((verts[:, 2], verts[:, 1], verts[:, 0]))
        buf[:] = verts[order[0]]
    v = np.zeros_like(buf)
    comm.Allreduce(buf, v, op=MPI.SUM)
    return v.reshape(1, 3)


def _probe_set(row, comm):
    """The 25 probes plus V (V defined on the degree-1 mesh; same mesh at d=2)."""
    return np.ascontiguousarray(
        np.vstack([_interior_probe_points(), _vertex_v(row, comm)])
    )


def test_s_single_incidence_noop(rows):
    comm = MPI.COMM_WORLD
    points = np.ascontiguousarray(_interior_probe_points())
    for d in (1, 2):
        row = rows[d]
        msh = _mesh_of(row)
        n_inc = _n_incident(msh, points, comm)
        rr, ri = _sample(row, points, comm, RULE)
        lr, li = _sample(row, points, comm, None)
        _, owner_local = point_owner_cells(msh, points, comm)
        single = np.flatnonzero(n_inc == 1)
        multi = np.flatnonzero(n_inc > 1)
        multi_rel = []
        for i in multi:
            e_c = _e_in_cell(row, points[i], int(owner_local[i]), comm)
            e_r = np.real(rr[i]) + 1j * np.real(ri[i])
            multi_rel.append(float(np.linalg.norm(e_r - e_c) / np.linalg.norm(e_c)))
        if comm.rank == 0:
            print(
                f"\n[POST-8 3b-ii-t1 S] -n {comm.size} degree {d}: "
                f"{single.size}/{points.shape[0]} probes single-incidence "
                f"(predicted 25/25); n_incident {n_inc.tolist()}; "
                f"multi-incidence rel vs owner-cell eval {multi_rel}"
            )
        assert np.all(n_inc >= 1)
        for arr in (rr, ri):
            assert _bitwise_equal_across_ranks(comm, arr)
        for i in single:
            assert rr[i].tobytes() == lr[i].tobytes(), (d, int(i))
            assert ri[i].tobytes() == li[i].tobytes(), (d, int(i))
        for rel in multi_rel:
            assert rel <= CLOSED_FORM_TOL, (d, multi_rel)


def test_n0_vertex_incident_spread(rows):
    comm = MPI.COMM_WORLD
    row = rows[1]
    msh = _mesh_of(row)
    v = _vertex_v(row, comm)
    n_inc = int(_n_incident(msh, v, comm)[0])
    cells = _owned_incident_cells(msh, v[0])
    local_mags = []
    for c in cells:
        p = v.reshape(1, 3)
        cc = np.array([int(c)], dtype=np.int32)
        e = np.real(np.asarray(row["fields"].e_real.eval(p, cc)).reshape(-1)) + 1j * np.real(
            np.asarray(row["fields"].e_imag.eval(p, cc)).reshape(-1)
        )
        local_mags.append(float(np.linalg.norm(e)))
    mags = np.concatenate([np.asarray(m, dtype=np.float64) for m in comm.allgather(local_mags)])
    spread = float(mags.max() / mags.min())

    owner_original, _ = point_owner_cells(msh, v, comm)
    f_dg = _dg0_index(msh)
    leg, ok = evaluate_vector_field_parallel(f_dg, v, comm, owner=None)
    leg = np.asarray(leg).reshape(-1)
    assert bool(np.all(ok))
    assert np.all(np.imag(leg) == 0.0)
    legacy_cell = int(np.real(leg[0]))
    if comm.rank == 0:
        print(
            f"\n[POST-8 3b-ii-t1 N0] -n {comm.size} degree 1: V = {v[0].tolist()}; "
            f"n_incident {n_inc} (asserted >= {N0_MIN_INCIDENT}); "
            f"max/min |E| over {mags.size} incident cells {spread:.6f} "
            f"(asserted > 1 + {N0_SPREAD_FLOOR:.0e}; predicted > 1.01); "
            f"(N1, printed) owner {int(owner_original[0])} vs legacy cell {legacy_cell}"
        )
    assert mags.size == n_inc
    assert n_inc >= N0_MIN_INCIDENT
    assert spread > 1.0 + N0_SPREAD_FLOOR


def test_o_w_width_invariance(rows):
    from tests.validation.test_post8_step3a_averaged_row1 import WIDTH_RTOL

    comm = MPI.COMM_WORLD
    points = _probe_set(rows[1], comm)
    owners = {}
    samples = {}
    for d in (1, 2):
        row = rows[d]
        msh = _mesh_of(row)
        owners[d], _ = point_owner_cells(msh, points, comm)
        rr, ri = _sample(row, points, comm, RULE)
        samples[d] = np.real(rr) + 1j * np.real(ri)
    if comm.size == 2:
        if comm.rank == 0:
            NPZ_DIR.mkdir(parents=True, exist_ok=True)
            np.savez(
                NPZ,
                points=points,
                owner1=np.asarray(owners[1], dtype=np.int64),
                owner2=np.asarray(owners[2], dtype=np.int64),
                e1=samples[1],
                e2=samples[2],
            )
            print(f"\n[POST-8 3b-ii-t1 W] -n 2 wrote {NPZ} ({points.shape[0]} points)")
        comm.Barrier()
        return
    comm.Barrier()
    if not NPZ.exists():
        pytest.fail(f"{NPZ} missing — run the -n 2 window (W1) first; never recompute")
    ref = np.load(NPZ)
    spreads = {}
    for d in (1, 2):
        e_ref = ref[f"e{d}"]
        spreads[d] = np.linalg.norm(samples[d] - e_ref, axis=1) / np.linalg.norm(e_ref, axis=1)
    if comm.rank == 0:
        print(
            f"\n[POST-8 3b-ii-t1 O/W] -n {comm.size} vs -n 2: V {points[-1].tolist()} vs "
            f"{ref['points'][-1].tolist()}; owners equal d1 "
            f"{bool(np.array_equal(owners[1], ref['owner1']))}, d2 "
            f"{bool(np.array_equal(owners[2], ref['owner2']))}; max spread d1 "
            f"{float(spreads[1].max()):.3e}, d2 {float(spreads[2].max()):.3e} "
            f"(asserted <= {WIDTH_RTOL:.0e}; predicted <= 5e-10)"
        )
    assert np.array_equal(points, ref["points"])
    assert np.array_equal(np.asarray(owners[1], dtype=np.int64), ref["owner1"])
    assert np.array_equal(np.asarray(owners[2], dtype=np.int64), ref["owner2"])
    for d in (1, 2):
        assert float(spreads[d].max()) <= WIDTH_RTOL, (d, spreads[d])
