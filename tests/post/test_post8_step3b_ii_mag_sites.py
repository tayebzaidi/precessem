"""`POST-8` step 3b-ii-t2 — the magnetostatic DG test sites under the owner rule.

PROJECT_PLAN §9 item 537; brief ``docs/planning/chunks/POST-8.md``
(``### §9 item 537`` and its 2026-10-18 18:00 review ruling, which binds).

Sites (each now passes ``owner="lowest_original_cell"``):
``test_coil_phantom_bfield_metrics.py`` centerline / ±x symmetry samples of
``b_dg0``, and ``test_straight_wire.py::_sample_radial``.  The CG1 ``b_lagrange``
print-only calls of the coil-phantom module stay on the legacy path (ruling).

Fixtures are spies: ``evaluate_vector_field_parallel`` is wrapped *on the
importing module* (``cpm`` / ``sw``), recording the field, a copy of the points
and the ``owner`` kwarg, and forwarding every argument unchanged.  The
coil-phantom spy runs the module's own test body, so its assertions are (X).

Asserted: (E) the edit is the one ruled (four ``cpm`` calls ``B_dg0`` / ``B`` /
``B_dg0`` / ``B`` with owners rule / None / rule / None; every ``sw`` call on the
rule); (S) owner sample bitwise equal across ranks, bytes = legacy at every
single-incidence probe, = ``.eval`` in the owner cell to ``CLOSED_FORM_TOL`` at
any multi-incidence probe; (O) coil-phantom owner original indices equal the
``-n 2`` npz exactly at ``-n 4`` / ``-n 1``; (N0) at V ``n_incident >= 4`` and
``max/min |B_dg0|`` over V's incident cells > 1 + 1e-6.
Printed only (ruling): (R) owner vs legacy per probe, the field-value width
spread beside the `MAG-18` LU cross-width floor, (N1) legacy cell at V.
"""

from pathlib import Path

import numpy as np
import pytest
from mpi4py import MPI

import tests.validation.test_coil_phantom_bfield_metrics as cpm
import tests.validation.test_straight_wire as sw
from fem_em_solver.post.evaluation import (
    evaluate_vector_field_parallel,
    point_owner_cells,
)
from tests.post.test_point_owner_rule import (
    CLOSED_FORM_TOL,
    _bitwise_equal_across_ranks,
)
from tests.post.test_post8_step3b_ii_th_sites import (
    N0_MIN_INCIDENT,
    N0_SPREAD_FLOOR,
    _n_incident,
    _owned_incident_cells,
)
from tests.post.test_post8_step3b_point_sar_owner import _dg0_index
from tests.tolerances import FIELD_SCALE_FLOOR

NPZ_DIR = Path("/workspace/logs/post8")
#: New name; never step3a-n2.npz / step3b-n2.npz / step3b-ii-t1-*.npz.
NPZ = NPZ_DIR / "step3b-ii-t2-n2.npz"
RULE = "lowest_original_cell"


def _spy(calls):
    real = evaluate_vector_field_parallel

    def spy(function, points, *args, **kwargs):
        values, valid = real(function, points, *args, **kwargs)
        calls.append(
            {
                "f": function,
                "points": np.array(points, dtype=np.float64, copy=True),
                "owner": kwargs.get("owner"),
                "values": np.array(values, copy=True),
            }
        )
        return values, valid

    return spy


@pytest.fixture(scope="module")
def cp():
    """Run the unedited coil-phantom test body under the spy (its asserts = (X))."""
    calls = []
    mp = pytest.MonkeyPatch()
    try:
        mp.setattr(cpm, "evaluate_vector_field_parallel", _spy(calls))
        cpm.test_coil_phantom_bfield_metrics_are_finite_smooth_and_symmetric()
    finally:
        mp.undo()
    return calls


def _b_in_cell(f, point, local_cell, comm):
    """Real field ``f`` at ``point`` in ``local_cell`` (>= 0 on the owning rank).

    Real-field analogue of t1's ``_e_in_cell``: the owning rank evaluates, the
    result is shared by buffer ``Allreduce`` SUM, independent of the rule's
    gather path.
    """
    buf = np.zeros(3, dtype=np.float64)
    if local_cell >= 0:
        p = np.asarray(point, dtype=np.float64).reshape(1, 3)
        c = np.array([int(local_cell)], dtype=np.int32)
        buf[:] = np.asarray(f.eval(p, c), dtype=np.float64).reshape(-1)
    tot = np.zeros_like(buf)
    comm.Allreduce(buf, tot, op=MPI.SUM)
    return tot


def _legacy_cells(msh, points, comm):
    leg, ok = evaluate_vector_field_parallel(_dg0_index(msh), points, comm, owner=None)
    assert bool(np.all(ok))
    return np.rint(np.asarray(leg).reshape(-1)).astype(np.int64)


def _s_and_r(tag, f, points, comm):
    """(S) asserted, (R) printed, on one captured field and point set."""
    msh = f.function_space.mesh
    n_inc = _n_incident(msh, points, comm)  # measured first
    own, ok_o = evaluate_vector_field_parallel(f, points, comm, owner=RULE)
    leg, ok_l = evaluate_vector_field_parallel(f, points, comm, owner=None)
    assert bool(np.all(ok_o)) and bool(np.all(ok_l))
    own = np.asarray(own, dtype=np.float64)
    leg = np.asarray(leg, dtype=np.float64)
    owner_orig, owner_local = point_owner_cells(msh, points, comm)
    legacy = _legacy_cells(msh, points, comm)
    single = np.flatnonzero(n_inc == 1)
    multi = np.flatnonzero(n_inc > 1)
    multi_rel = []
    r_rows = []
    for i in multi:
        b_c = _b_in_cell(f, points[i], int(owner_local[i]), comm)
        multi_rel.append(float(np.linalg.norm(own[i] - b_c) / np.linalg.norm(b_c)))
        delta = float(
            np.linalg.norm(own[i] - leg[i]) / max(np.linalg.norm(leg[i]), FIELD_SCALE_FLOOR)
        )
        r_rows.append(
            f"p{int(i)} n_inc {int(n_inc[i])} owner {int(owner_orig[i])} "
            f"legacy {int(legacy[i])} delta {delta:.3e}"
        )
    r_equal = sum(own[i].tobytes() == leg[i].tobytes() for i in range(points.shape[0]))
    if comm.rank == 0:
        print(
            f"\n[POST-8 3b-ii-t2 S] {tag} -n {comm.size}: "
            f"{single.size}/{points.shape[0]} probes single-incidence; "
            f"n_incident {n_inc.tolist()}; multi-incidence rel vs owner-cell eval "
            f"{[f'{x:.3e}' for x in multi_rel]}"
        )
        print(
            f"[POST-8 3b-ii-t2 R, printed] {tag} -n {comm.size}: owner bytes == legacy "
            f"bytes at {r_equal}/{points.shape[0]} probes; multi-incidence: "
            f"{r_rows if r_rows else 'none'}"
        )
    assert np.all(n_inc >= 1)
    assert _bitwise_equal_across_ranks(comm, own)
    for i in single:
        assert own[i].tobytes() == leg[i].tobytes(), (tag, int(i))
    for rel in multi_rel:
        assert rel <= CLOSED_FORM_TOL, (tag, multi_rel)
    return own, leg


def _mirror_metric(b):
    """The ±x pair metric as the coil-phantom module forms it (``:211–218``)."""
    mag = np.linalg.norm(b, axis=1).reshape(-1, 2)
    ref = np.maximum(np.maximum(mag[:, 0], mag[:, 1]), FIELD_SCALE_FLOOR)
    return float(np.max(np.abs(mag[:, 0] - mag[:, 1]) / ref))


def _check_e_cp(cp):
    names = [c["f"].name for c in cp]
    owners = [c["owner"] for c in cp]
    assert names == ["B_dg0", "B", "B_dg0", "B"], names
    assert owners == [RULE, None, RULE, None], owners
    assert cp[0]["f"] is cp[2]["f"]
    return cp[0]["f"], cp[0]["points"], cp[2]["points"]


def _vertex_v(b_dg0, centerline, comm):
    """(V): lexicographically smallest vertex of centerline probe 0's owner cell."""
    msh = b_dg0.function_space.mesh
    _, owner_local = point_owner_cells(msh, centerline[:1], comm)
    buf = np.zeros(3, dtype=np.float64)
    if owner_local[0] >= 0:
        verts = msh.geometry.x[msh.geometry.dofmap[int(owner_local[0])]][:, :3]
        order = np.lexsort((verts[:, 2], verts[:, 1], verts[:, 0]))
        buf[:] = verts[order[0]]
    v = np.zeros_like(buf)
    comm.Allreduce(buf, v, op=MPI.SUM)
    return v.reshape(1, 3)


def test_cp_s_r_n0(cp):
    comm = MPI.COMM_WORLD
    b_dg0, centerline, symmetry = _check_e_cp(cp)
    msh = b_dg0.function_space.mesh

    own_c, _ = _s_and_r("coil-phantom centerline", b_dg0, centerline, comm)
    own_s, leg_s = _s_and_r("coil-phantom symmetry", b_dg0, symmetry, comm)
    # The module's own gated samples are the owner-path samples re-taken here.
    assert own_c.tobytes() == np.asarray(cp[0]["values"], dtype=np.float64).tobytes()
    assert own_s.tobytes() == np.asarray(cp[2]["values"], dtype=np.float64).tobytes()
    if comm.rank == 0:
        print(
            f"[POST-8 3b-ii-t2 R, printed] -n {comm.size}: mirror metric owner path "
            f"{_mirror_metric(own_s):.6f}, legacy path {_mirror_metric(leg_s):.6f} "
            "(record: cpm docstring :30-31 / message :234, printed by the module above)"
        )

    # (N0) asserted; (N1) printed.
    v = _vertex_v(b_dg0, centerline, comm)
    n_inc = int(_n_incident(msh, v, comm)[0])
    local_mags = [
        float(np.linalg.norm(_b_local(b_dg0, v, c))) for c in _owned_incident_cells(msh, v[0])
    ]
    mags = np.concatenate([np.asarray(m, dtype=np.float64) for m in comm.allgather(local_mags)])
    spread = float(mags.max() / mags.min())
    owner_original, _ = point_owner_cells(msh, v, comm)
    legacy_cell = int(_legacy_cells(msh, v, comm)[0])
    if comm.rank == 0:
        print(
            f"\n[POST-8 3b-ii-t2 N0] -n {comm.size}: V = "
            f"{[float(f'{x:.6g}') for x in v[0]]}; n_incident {n_inc} "
            f"(asserted >= {N0_MIN_INCIDENT}); max/min |B_dg0| over {mags.size} "
            f"incident cells {spread:.6f} (asserted > 1 + {N0_SPREAD_FLOOR:.0e}; "
            f"predicted > 1.01); (N1, printed) owner {int(owner_original[0])} vs "
            f"legacy cell {legacy_cell}"
        )
    assert mags.size == n_inc
    assert n_inc >= N0_MIN_INCIDENT
    assert spread > 1.0 + N0_SPREAD_FLOOR


def _b_local(f, v, local_cell):
    p = np.asarray(v, dtype=np.float64).reshape(1, 3)
    c = np.array([int(local_cell)], dtype=np.int32)
    return np.asarray(f.eval(p, c), dtype=np.float64).reshape(-1)


def test_cp_o_width(cp):
    comm = MPI.COMM_WORLD
    b_dg0, centerline, symmetry = _check_e_cp(cp)
    msh = b_dg0.function_space.mesh
    probes = np.ascontiguousarray(np.vstack([centerline, symmetry]))
    points = np.ascontiguousarray(np.vstack([probes, _vertex_v(b_dg0, centerline, comm)]))
    owners, _ = point_owner_cells(msh, points, comm)
    owners = np.asarray(owners, dtype=np.int64)
    n_inc = np.asarray(_n_incident(msh, points, comm), dtype=np.int64)
    b, ok = evaluate_vector_field_parallel(b_dg0, probes, comm, owner=RULE)
    assert bool(np.all(ok))
    b = np.asarray(b, dtype=np.float64)
    if comm.size == 2:
        if comm.rank == 0:
            NPZ_DIR.mkdir(parents=True, exist_ok=True)
            np.savez(NPZ, points=points, owner=owners, n_incident=n_inc, b=b)
            print(f"\n[POST-8 3b-ii-t2 O] -n 2 wrote {NPZ} ({points.shape[0]} points)")
        comm.Barrier()
        return
    comm.Barrier()
    if not NPZ.exists():
        pytest.fail(f"{NPZ} missing — run the -n 2 window (W1) first; never recompute")
    ref = np.load(NPZ)
    spread = np.linalg.norm(b - ref["b"], axis=1) / np.linalg.norm(ref["b"], axis=1)
    if comm.rank == 0:
        print(
            f"\n[POST-8 3b-ii-t2 O] -n {comm.size} vs -n 2: points equal "
            f"{bool(np.array_equal(points, ref['points']))}; owners equal "
            f"{bool(np.array_equal(owners, ref['owner']))}; n_incident equal "
            f"{bool(np.array_equal(n_inc, ref['n_incident']))}; (printed only) max "
            f"field spread {float(spread.max()):.3e} beside the MAG-18 LU cross-width "
            "floor 7.28e-08 (test_straight_wire.py:113-122, a different fixture); "
            "predicted <= 1e-6, no band"
        )
    assert np.array_equal(points, ref["points"])
    assert np.array_equal(owners, ref["owner"])


def test_sw_s_r():
    comm = MPI.COMM_WORLD
    calls = []
    mp = pytest.MonkeyPatch()
    fields = []
    try:
        mp.setattr(sw, "evaluate_vector_field_parallel", _spy(calls))
        _, b_fine = sw._solve_straight_wire(sw.RESOLUTION, comm)
        n0 = len(calls)
        for n in sorted(sw._npoints_control()[1]):
            sw._sample_radial(b_fine, n, comm, sw.R_MIN, sw.R_MAX_BC)
        sw._sample_radial(b_fine, 8, comm, sw.R_MIN, sw.R_MAX)
        fields.append((f"straight wire h={sw.RESOLUTION}", b_fine, calls[n0:]))
        _, b_coarse = sw._solve_straight_wire(sw.E_OMEGA_LADDER[0], comm)
        n1 = len(calls)
        sw._sample_radial(b_coarse, 10, comm, sw.R_MIN, sw.R_MAX_BC)
        sw._sample_radial(b_coarse, 8, comm, sw.R_MIN, sw.R_MAX)
        fields.append((f"straight wire h={sw.E_OMEGA_LADDER[0]}", b_coarse, calls[n1:]))
    finally:
        mp.undo()

    # (E): every sw call on the rule.
    assert len(calls) == 6, len(calls)
    assert all(c["owner"] == RULE for c in calls), [c["owner"] for c in calls]
    for tag, f, cs in fields:
        assert all(c["f"] is f for c in cs)
        pts = np.ascontiguousarray(np.vstack([c["points"] for c in cs]))
        _s_and_r(tag, f, pts, comm)
        # The module's gated samples are the owner-path samples, compared call
        # by call on the *same* batch.  Measured 2026-10-19
        # (20261019T020649Z_POST-8-step3b-ii-t2-W4-wire.log): re-sampling the
        # four h = 0.0025 batches as one 46-point batch differs from the
        # captured values in the last ULP at some probes — ``.eval`` bytes
        # depend on the batch a point is evaluated in — so the first form of
        # this check (one concatenated batch) was a wrong premise, not a rule
        # defect.  (S) compares owner vs legacy on one batch and is unaffected.
        for c in cs:
            own_c, ok = evaluate_vector_field_parallel(f, c["points"], comm, owner=RULE)
            assert bool(np.all(ok))
            assert (
                np.asarray(own_c, dtype=np.float64).tobytes()
                == np.asarray(c["values"], dtype=np.float64).tobytes()
            )
