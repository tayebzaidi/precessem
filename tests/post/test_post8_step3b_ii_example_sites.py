"""`POST-8` step 3b-ii-e1 — seven example call sites under the owner rule.

PROJECT_PLAN §9 item 547; brief ``docs/planning/chunks/POST-8.md``
(``### §9 item 547`` and its 2026-10-19 03:00 review ruling, which binds:
R1 declined, ``mag:1`` takes the keyword on ``B_lag`` only).

Each test runs one edited example as ``__main__`` (``runpy.run_path``) under a
spy on ``evaluate_vector_field_parallel`` (patched on both
``fem_em_solver.post.evaluation`` and ``fem_em_solver.post``).  The example's
own assertions run unmodified inside the run: that is (X).  The spy records
only calls whose caller frame is the example's own file, forwards every
argument unchanged, and on every rank adds side evaluations on the *identical*
point array of that call.

Asserted: (X) the example exits 0; (E) every recorded call carries
``owner == "lowest_original_cell"``, at least one call, forwarded bytes = the
rule's bytes; (S) rule sample bitwise equal across ranks, bytes = legacy bytes at
every single-incidence probe of each rank's owner batch (S-sb, same batch; e2 sites:
cross-batch (S) printed per the weekly 2026-10-19 ruling), = in-cell ``.eval`` to
``CLOSED_FORM_TOL`` at a multi-incidence probe; (O) at ``-n 4`` / ``-n 1`` the
call count, every ``n_incident`` array and every owner original-index array
equal the ``-n 2`` npz exactly; (N0) at V (lexicographically smallest vertex of
probe 0's owner cell) ``n_incident >= 4`` and ``max/min |f|`` over V's incident
cells > 1 + 1e-6 — ``th:3``'s first ``e_real`` call and ``mag:5``'s penalty
call.  Printed only: (R) owner vs legacy per multi-incidence probe, the
field-value width spread, (N1) the legacy cell at V.

Item 568 Part B: ``test_ports1_owner_measure`` measures ``ports:1`` under the
owner rule by injecting the keyword through the spy (the example on ``main``
stays on the legacy path); nothing is re-recorded.

Item 586 (e4): ``ports:1`` carries the keyword on ``main``; ``test_ports1``
checks it at ``-n 1 / 2 / 4``; ``test_ports1_owner_measure`` is item 568's
record (its injection is now a no-op).

Gated by ``FEM_EM_POST8_EXAMPLES=1`` so no CI job runs whole examples.
"""

import os
import runpy
import sys
from pathlib import Path

import numpy as np
import pytest
from mpi4py import MPI

import fem_em_solver.post as post_pkg
import fem_em_solver.post.evaluation as post_eval
from fem_em_solver.post.evaluation import (
    evaluate_vector_field_parallel,
    point_owner_cells,
)
from tests.complex_mode import IS_COMPLEX, complex_only
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

pytestmark = [
    pytest.mark.slow,
    pytest.mark.skipif(
        os.environ.get("FEM_EM_POST8_EXAMPLES") != "1",
        reason="whole-example windows only: set FEM_EM_POST8_EXAMPLES=1",
    ),
]

REPO = Path(__file__).resolve().parents[2]
NPZ_DIR = Path("/workspace/logs/post8")
RULE = "lowest_original_cell"
#: tag -> (path, record width = digits of the guide's records)
EXAMPLES = {
    "mag1": (REPO / "examples/magnetostatics/01_straight_wire.py", 2),
    "mag5": (REPO / "examples/magnetostatics/05_gauge_cross_check.py", 2),
    "th3": (REPO / "examples/time_harmonic/03_dielectric_sphere_in_uniform_field.py", 2),
    "th4": (REPO / "examples/time_harmonic/04_evanescent_waveguide_decay.py", 2),
    "th6": (REPO / "examples/time_harmonic/06_larmor_lossy_sphere.py", 2),
    "th7": (REPO / "examples/time_harmonic/07_element_order_lossy_sphere.py", 2),
    "th9": (REPO / "examples/time_harmonic/09_pec_sphere_as_hole.py", 2),
}
PREDICTED_CALLS = {"mag1": 1, "mag5": 2, "th3": 4, "th4": 9, "th6": 4, "th7": 4, "th9": 4}
#: calls (by index) carrying the asserted (N0): th3 first e_real, mag5 penalty (first).
N0_CALL = {"th3": 0, "mag5": 0}


def _comm_of(args, kwargs):
    if kwargs.get("comm") is not None:
        return kwargs["comm"]
    if args:
        return args[0]
    return MPI.COMM_WORLD


def _family(f):
    try:
        return str(f.function_space.ufl_element().family_name)
    except Exception:  # pragma: no cover - printing aid only
        return "?"


def _ncomp(f):
    shape = f.function_space.element.value_shape
    return int(np.prod(shape)) if len(shape) > 0 else 1


def _in_cell(f, point, local_cell, comm):
    """``f`` at ``point`` in ``local_cell`` (>= 0 on the owning rank only).

    Dtype-preserving analogue of t2's ``_b_in_cell`` (test_post8_step3b_ii_mag_sites.py
    :89–103, which casts to float64) and t1's ``_e_in_cell``; the buffer dtype is
    ``f.x.array.dtype``; shared by buffer ``Allreduce`` SUM.
    """
    buf = np.zeros(_ncomp(f), dtype=f.x.array.dtype)
    if local_cell >= 0:
        p = np.asarray(point, dtype=np.float64).reshape(1, 3)
        c = np.array([int(local_cell)], dtype=np.int32)
        buf[:] = np.asarray(f.eval(p, c)).reshape(-1)
    tot = np.zeros_like(buf)
    comm.Allreduce(buf, tot, op=MPI.SUM)
    return tot


def _legacy_cells(msh, points, comm):
    leg, ok = evaluate_vector_field_parallel(_dg0_index(msh), points, comm, owner=None)
    assert bool(np.all(ok))
    leg = np.asarray(leg).reshape(-1)
    assert np.all(np.imag(leg) == 0.0)
    return np.rint(np.real(leg)).astype(np.int64)


def _make_spy(path, calls, inject_owner=None):
    real = evaluate_vector_field_parallel
    target = path.resolve()

    def spy(function, points, *args, **kwargs):
        caller = sys._getframe(1)
        is_target = Path(caller.f_code.co_filename).resolve() == target
        if inject_owner is not None and is_target and "owner" not in kwargs:
            kwargs = {**kwargs, "owner": inject_owner}
        fwd = real(function, points, *args, **kwargs)
        if not is_target:
            return fwd
        comm = _comm_of(args, kwargs)
        pts = np.array(points, dtype=np.float64, copy=True)
        msh = function.function_space.mesh
        n_inc = _n_incident(msh, pts, comm)  # measured first
        own, ok_o = real(function, points, comm, owner=RULE)
        leg, ok_l = real(function, points, comm, owner=None)
        assert bool(np.all(ok_o)) and bool(np.all(ok_l))
        own = np.array(own, copy=True)
        leg = np.array(leg, copy=True)
        owner_orig, owner_local = point_owner_cells(msh, pts, comm)
        legacy = _legacy_cells(msh, pts, comm)
        multi = np.flatnonzero(n_inc > 1)
        in_cell = {int(i): _in_cell(function, pts[i], int(owner_local[i]), comm) for i in multi}
        calls.append(
            {
                "line": int(caller.f_lineno),
                "name": function.name,
                "family": _family(function),
                "owner": kwargs.get("owner"),
                "f": function,
                "pts": pts,
                "n_inc": np.asarray(n_inc, dtype=np.int64),
                "owner_orig": np.asarray(owner_orig, dtype=np.int64),
                "legacy": legacy,
                "own": own,
                "leg": leg,
                "fwd": np.array(fwd[0], copy=True),
                "in_cell": in_cell,
            }
        )
        return fwd

    return spy


def _run_spied(tag):
    path, _ = EXAMPLES[tag]
    calls = []
    mp = pytest.MonkeyPatch()
    try:
        spy = _make_spy(path, calls)
        mp.setattr(post_eval, "evaluate_vector_field_parallel", spy)
        mp.setattr(post_pkg, "evaluate_vector_field_parallel", spy)
        mp.setattr(sys, "argv", [str(path)])
        mp.syspath_prepend(str(path.parent))
        mp.syspath_prepend(str(REPO))
        try:
            runpy.run_path(str(path), run_name="__main__")
        except SystemExit as exc:
            if exc.code not in (0, None):
                raise
    finally:
        mp.undo()
    return calls


def _run_spied_measure(tag):
    """`_run_spied` for item 568 Part B: the spy injects ``owner=RULE`` on the
    example's own calls (the example on ``main`` stays on the legacy path), two
    pass-through recorders capture the ladder and S calls, and the example's
    0.11-image record AssertionError (only) is caught by its message."""
    import fem_em_solver.ports.sparameters as sp_mod
    import fem_em_solver.ports.systematics as sys_mod

    path, _ = EXAMPLES[tag]
    calls, ladder_calls, s_calls = [], [], []
    record_msg = None
    real_ladder = sys_mod.mutual_systematics_ladder
    real_s = sp_mod.sparameters_from_impedance

    def rec_ladder(*args, **kwargs):
        out = real_ladder(*args, **kwargs)
        ladder_calls.append((args, kwargs, out))
        return out

    def rec_s(*args, **kwargs):
        out = real_s(*args, **kwargs)
        s_calls.append((args, kwargs, out))
        return out

    mp = pytest.MonkeyPatch()
    try:
        spy = _make_spy(path, calls, inject_owner=RULE)
        mp.setattr(post_eval, "evaluate_vector_field_parallel", spy)
        mp.setattr(post_pkg, "evaluate_vector_field_parallel", spy)
        mp.setattr(sys_mod, "mutual_systematics_ladder", rec_ladder)
        mp.setattr(sp_mod, "sparameters_from_impedance", rec_s)
        mp.setattr(sys, "argv", [str(path)])
        mp.syspath_prepend(str(path.parent))
        mp.syspath_prepend(str(REPO))
        try:
            runpy.run_path(str(path), run_name="__main__")
        except SystemExit as exc:
            if exc.code not in (0, None):
                raise
        except AssertionError as exc:
            if "does not reproduce the 0.11-image record" not in str(exc):
                raise
            record_msg = str(exc)
    finally:
        mp.undo()
    return calls, ladder_calls, s_calls, record_msg


def _check_e_s(tag, calls, comm):
    assert len(calls) >= 1, f"{tag}: no call recorded"
    n_single = n_total = 0
    for k, c in enumerate(calls):
        assert c["owner"] == RULE, (tag, k, c["owner"])
        assert c["fwd"].tobytes() == c["own"].tobytes(), (tag, k)
        assert np.all(c["n_inc"] >= 1), (tag, k)
        assert _bitwise_equal_across_ranks(comm, c["own"]), (tag, k)
        single = np.flatnonzero(c["n_inc"] == 1)
        for i in single:
            assert c["own"][i].tobytes() == c["leg"][i].tobytes(), (tag, k, int(i))
        multi_rel, r_rows = [], []
        for i, b_c in c["in_cell"].items():
            rel = float(
                np.linalg.norm(c["own"][i] - b_c) / max(np.linalg.norm(b_c), FIELD_SCALE_FLOOR)
            )
            multi_rel.append(rel)
            delta = float(
                np.linalg.norm(c["own"][i] - c["leg"][i])
                / max(np.linalg.norm(c["leg"][i]), FIELD_SCALE_FLOOR)
            )
            r_rows.append(
                f"p{i} n_inc {int(c['n_inc'][i])} owner {int(c['owner_orig'][i])} "
                f"legacy {int(c['legacy'][i])} delta {delta:.3e}"
            )
        n_single += single.size
        n_total += c["n_inc"].size
        if comm.rank == 0:
            print(
                f"[POST-8 3b-ii-e1 S] {tag} -n {comm.size} call {k} (line {c['line']}, "
                f"{c['name']}, {c['family']}, owner {c['owner']}): "
                f"{single.size}/{c['n_inc'].size} probes single-incidence; "
                f"multi rel vs owner-cell eval {[f'{x:.3e}' for x in multi_rel]}; "
                f"(R, printed) {r_rows if r_rows else 'none'}",
                flush=True,
            )
        for rel in multi_rel:
            assert rel <= CLOSED_FORM_TOL, (tag, k, multi_rel)
    if comm.rank == 0:
        print(
            f"[POST-8 3b-ii-e1 E] {tag} -n {comm.size}: {len(calls)} calls recorded "
            f"(predicted {PREDICTED_CALLS[tag]}); {n_single}/{n_total} probes "
            "single-incidence across calls",
            flush=True,
        )


def _vertex_v(f, pts0, comm):
    msh = f.function_space.mesh
    _, owner_local = point_owner_cells(msh, pts0, comm)
    buf = np.zeros(3, dtype=np.float64)
    if owner_local[0] >= 0:
        verts = msh.geometry.x[msh.geometry.dofmap[int(owner_local[0])]][:, :3]
        order = np.lexsort((verts[:, 2], verts[:, 1], verts[:, 0]))
        buf[:] = verts[order[0]]
    v = np.zeros_like(buf)
    comm.Allreduce(buf, v, op=MPI.SUM)
    return v.reshape(1, 3)


def _check_n0(tag, calls, comm, asserted):
    f = calls[N0_CALL.get(tag, 0)]["f"]
    pts0 = calls[N0_CALL.get(tag, 0)]["pts"][:1]
    msh = f.function_space.mesh
    v = _vertex_v(f, pts0, comm)
    n_inc = int(_n_incident(msh, v, comm)[0])
    local = []
    for c in _owned_incident_cells(msh, v[0]):
        p = v.reshape(1, 3)
        cc = np.array([int(c)], dtype=np.int32)
        local.append(float(np.linalg.norm(np.asarray(f.eval(p, cc)).reshape(-1))))
    mags = np.concatenate([np.asarray(m, dtype=np.float64) for m in comm.allgather(local)])
    spread = float(mags.max() / mags.min())
    owner_original, _ = point_owner_cells(msh, v, comm)
    legacy_cell = int(_legacy_cells(msh, v, comm)[0])
    if comm.rank == 0:
        print(
            f"[POST-8 3b-ii-e1 N0] {tag} -n {comm.size}: V = "
            f"{[float(f'{x:.6g}') for x in v[0]]}; n_incident {n_inc} "
            f"(asserted >= {N0_MIN_INCIDENT}); max/min |f| over {mags.size} incident cells "
            f"{spread:.6f} (asserted > 1 + {N0_SPREAD_FLOOR:.0e}; predicted > 1.007); "
            f"(N1, printed) owner {int(owner_original[0])} vs legacy cell {legacy_cell}",
            flush=True,
        )
    if asserted:
        assert mags.size == n_inc
        assert n_inc >= N0_MIN_INCIDENT
        assert spread > 1.0 + N0_SPREAD_FLOOR


def _check_o(tag, calls, comm):
    npz = NPZ_DIR / f"step3b-ii-e1-{tag}-n2.npz"  # new names only
    if comm.size == 2:
        if comm.rank == 0:
            NPZ_DIR.mkdir(parents=True, exist_ok=True)
            data = {"ncalls": np.array(len(calls))}
            for k, c in enumerate(calls):
                data[f"line_{k}"] = np.array(c["line"])
                data[f"owner_orig_{k}"] = c["owner_orig"]
                data[f"n_inc_{k}"] = c["n_inc"]
                data[f"own_{k}"] = c["own"]
            np.savez(npz, **data)
            print(f"[POST-8 3b-ii-e1 O] {tag} -n 2 wrote {npz} ({len(calls)} calls)", flush=True)
        comm.Barrier()
        return
    comm.Barrier()
    if not npz.exists():
        pytest.fail(f"{npz} missing — run the -n 2 window first; never recompute")
    ref = np.load(npz)
    assert len(calls) == int(ref["ncalls"]), (tag, len(calls), int(ref["ncalls"]))
    spread = 0.0
    for k, c in enumerate(calls):
        assert c["line"] == int(ref[f"line_{k}"]), (tag, k)
        assert np.array_equal(c["n_inc"], ref[f"n_inc_{k}"]), (tag, k)
        assert np.array_equal(c["owner_orig"], ref[f"owner_orig_{k}"]), (tag, k)
        a, b = c["own"], ref[f"own_{k}"]
        assert a.shape == b.shape, (tag, k)
        num = np.linalg.norm(a - b, axis=1)
        den = np.maximum(np.linalg.norm(b, axis=1), FIELD_SCALE_FLOOR)
        spread = max(spread, float(np.max(num / den)))
    if comm.rank == 0:
        print(
            f"[POST-8 3b-ii-e1 O] {tag} -n {comm.size} vs -n 2: {len(calls)} calls, "
            f"n_incident and owner original indices exactly equal on every call; "
            f"field spread (printed) {spread:.3e}",
            flush=True,
        )


def _drive(tag):
    comm = MPI.COMM_WORLD
    calls = _run_spied(tag)
    _check_e_s(tag, calls, comm)
    if tag in N0_CALL:  # th:3 first e_real call and mag:5 penalty call (both DG)
        _check_n0(tag, calls, comm, asserted=True)
    _check_o(tag, calls, comm)


@pytest.mark.skipif(IS_COMPLEX, reason="magnetostatic example: real build")
def test_mag1():
    _drive("mag1")


@pytest.mark.skipif(IS_COMPLEX, reason="magnetostatic example: real build")
def test_mag5():
    _drive("mag5")


@complex_only
def test_th3():
    _drive("th3")


@complex_only
def test_th4():
    _drive_e3("th4")


@complex_only
def test_th6():
    _drive("th6")


@complex_only
def test_th7():
    _drive("th7")


@complex_only
def test_th9():
    _drive("th9")


# ---------------------------------------------------------------------------
# `POST-8` step 3b-ii-e2 (§9 item 548): five expensive example sites, record
# width only, additive to the e1 block above (no e1 name is changed).  There is
# no (O) here: each run writes its own ``step3b-ii-e2-<tag>-n<w>.npz`` for a
# future width check.  (N0) asserted on mag:2 (DG B) and ports:12 (raw N1curl
# e_complex, first call); printed for mag:4 / ports:1; skipped for mag:6 (CG1).
# ---------------------------------------------------------------------------
EXAMPLES.update(
    {
        "mag2": (REPO / "examples/magnetostatics/02_circular_loop.py", 2),
        "mag4": (REPO / "examples/magnetostatics/04_helmholtz_analytic_comparison.py", 2),
        "mag6": (REPO / "examples/magnetostatics/06_h_convergence_rate.py", 2),
        "ports1": (REPO / "examples/ports/01_two_torus_port_pair.py", 2),
        "ports12": (REPO / "examples/ports/12_birdcage_ring_quarter_turn.py", 4),
    }
)
PREDICTED_CALLS.update({"mag2": 1, "mag4": 3, "mag6": 1, "ports1": 8, "ports12": 2})
# POST-8 3b-ii-e3 (item 556): mri:1 keyword-only (2 calls)
EXAMPLES["mri1"] = (REPO / "examples/mri/01_coil_phantom_fields.py", 2)
PREDICTED_CALLS["mri1"] = 2
N0_CALL.update({"mag2": 0})
N0_PRINTED_E2 = ("mag4", "ports1", "ports12")  # printed only; mag6 (CG1) skipped


def _write_e2_npz(tag, calls, comm):
    npz = NPZ_DIR / f"step3b-ii-e2-{tag}-n{comm.size}.npz"
    if comm.rank == 0:
        NPZ_DIR.mkdir(parents=True, exist_ok=True)
        data = {"ncalls": np.array(len(calls))}
        for k, c in enumerate(calls):
            data[f"line_{k}"] = np.array(c["line"])
            data[f"owner_orig_{k}"] = c["owner_orig"]
            data[f"n_inc_{k}"] = c["n_inc"]
            data[f"own_{k}"] = c["own"]
        np.savez(npz, **data)
        print(
            f"[POST-8 3b-ii-e2 npz] {tag} -n {comm.size} wrote {npz} ({len(calls)} calls)",
            flush=True,
        )
    comm.Barrier()


def _legacy_hits(msh, pts):
    """Legacy-path locator on this rank: {point: (first colliding cell, owned?)}.

    Mirrors ``evaluate_vector_field_parallel(owner=None)`` (evaluation.py:380-388):
    ``links[0]`` of ``compute_colliding_cells``; ``owned`` is ``cell < size_local``.
    """
    from dolfinx import geometry

    tree = geometry.bb_tree(msh, msh.topology.dim)
    cand = geometry.compute_collisions_points(tree, pts)
    coll = geometry.compute_colliding_cells(msh, cand, pts)
    n_owned = msh.topology.index_map(msh.topology.dim).size_local
    hits = {}
    for i in range(pts.shape[0]):
        links = coll.links(i)
        if len(links) > 0:
            hits[i] = (int(links[0]), bool(links[0] < n_owned))
    return hits


def _ulp_distance(a, b):
    """Per-component distance in float64 ULP steps (complex: real and imag parts)."""
    a = np.ascontiguousarray(np.asarray(a)).view(np.float64).reshape(-1)
    b = np.ascontiguousarray(np.asarray(b)).view(np.float64).reshape(-1)
    return np.frombuffer(a.tobytes(), np.int64) - np.frombuffer(b.tobytes(), np.int64)


def _check_e2_s(tag, calls, comm):
    """`_check_e_s`'s logic (item 551) in the ruled (S) form (weekly 2026-10-19,
    PROJECT_PLAN.md:4694-4705; item 555).

    (S-sb), asserted: on each rank's owner batch (``owner_local >= 0``) the legacy
    cell evaluated on exactly that batch is byte-equal to the rule's bytes at the
    single-incidence rows; multi-incidence rows print their same-batch count only
    (gated by the in-cell ``CLOSED_FORM_TOL`` check).  (S), cross-batch, printed:
    every byte mismatch against the full-array legacy bytes with its ULP distance,
    per-call and per-site counts, max |ULP| and the A1 / A2 hit-pattern split.
    The cross-batch byte assertion is retired (epitaph PROJECT_PLAN.md:4782-4785)."""
    assert len(calls) >= 1, f"{tag}: no call recorded"
    msh0 = calls[0]["f"].function_space.mesh
    ghosts = comm.allgather(int(msh0.topology.index_map(msh0.topology.dim).num_ghosts))
    if comm.rank == 0:
        print(
            f"[POST-8 3b-ii-e2a G] {tag} -n {comm.size}: cells_ghost per rank {ghosts} "
            "(mesh-class pre-condition: all 0; printed, not asserted)",
            flush=True,
        )
    n_single = n_total = 0
    sb_failures = []
    site_mism = site_a1 = site_a2 = site_ulp = 0
    for k, c in enumerate(calls):
        assert c["owner"] == RULE, (tag, k, c["owner"])
        assert c["fwd"].tobytes() == c["own"].tobytes(), (tag, k)
        assert np.all(c["n_inc"] >= 1), (tag, k)
        assert _bitwise_equal_across_ranks(comm, c["own"]), (tag, k)
        msh = c["f"].function_space.mesh
        pts = c["pts"]
        single = np.flatnonzero(c["n_inc"] == 1)
        mism = [int(i) for i in single if c["own"][i].tobytes() != c["leg"][i].tobytes()]
        hits = _legacy_hits(msh, pts)
        _, owner_local = point_owner_cells(msh, pts, comm)
        all_hits = comm.allgather({i: hits[i] for i in mism if i in hits})
        all_owner_local = comm.allgather({i: int(owner_local[i]) for i in mism})
        call_a1 = call_a2 = call_ulp = 0
        for i in mism:
            owning = [r for r in range(comm.size) if all_owner_local[r][i] >= 0]
            per_rank = {r: all_hits[r][i] for r in range(comm.size) if i in all_hits[r]}
            ulp = _ulp_distance(c["own"][i], c["leg"][i]).tolist()
            a2 = any(r not in owning for r in per_rank)
            call_a2 += int(a2)
            call_a1 += int(not a2)
            call_ulp = max([call_ulp] + [abs(int(u)) for u in ulp])
            if comm.rank == 0:
                print(
                    f"[POST-8 3b-ii-e2c S (S, cross-batch, printed)] {tag} call {k} probe {i}: n_inc "
                    f"{int(c['n_inc'][i])}, owner_orig {int(c['owner_orig'][i])}, "
                    f"legacy_orig {int(c['legacy'][i])}, owning rank {owning}, legacy "
                    f"hits (rank: (cell, owned)) {per_rank}, ULP distance per component {ulp}, "
                    f"pattern {'A2' if a2 else 'A1'}",
                    flush=True,
                )
        # (S-sb): legacy cells on exactly the owner batch's points and order.
        idx = np.flatnonzero(owner_local >= 0)
        n_sb = n_sb_multi = 0
        if idx.size:
            cells = np.array([hits[int(i)][0] for i in idx], dtype=np.int32)
            vals = np.asarray(c["f"].eval(pts[idx], cells)).reshape(idx.size, -1)
            own_rows = np.asarray(c["own"]).reshape(c["own"].shape[0], -1)[idx]
            inc_b = np.asarray(c["n_inc"])[idx]
            for j in range(idx.size):
                if vals[j].tobytes() != own_rows[j].tobytes():
                    if inc_b[j] == 1:
                        n_sb += 1
                    else:
                        n_sb_multi += 1
        n_sb_tot = comm.allreduce(n_sb, op=MPI.SUM)
        n_sb_multi_tot = comm.allreduce(n_sb_multi, op=MPI.SUM)
        if n_sb_tot:
            sb_failures.append((k, n_sb_tot))
        site_mism += len(mism)
        site_a1 += call_a1
        site_a2 += call_a2
        site_ulp = max(site_ulp, call_ulp)
        multi_rel, r_rows = [], []
        for i, b_c in c["in_cell"].items():
            rel = float(
                np.linalg.norm(c["own"][i] - b_c) / max(np.linalg.norm(b_c), FIELD_SCALE_FLOOR)
            )
            multi_rel.append(rel)
            delta = float(
                np.linalg.norm(c["own"][i] - c["leg"][i])
                / max(np.linalg.norm(c["leg"][i]), FIELD_SCALE_FLOOR)
            )
            r_rows.append(
                f"p{i} n_inc {int(c['n_inc'][i])} owner {int(c['owner_orig'][i])} "
                f"legacy {int(c['legacy'][i])} delta {delta:.3e}"
            )
        n_single += single.size
        n_total += c["n_inc"].size
        if comm.rank == 0:
            print(
                f"[POST-8 3b-ii-e2a S] {tag} -n {comm.size} call {k} (line {c['line']}, "
                f"{c['name']}, {c['family']}, owner {c['owner']}): "
                f"{single.size}/{c['n_inc'].size} probes single-incidence; "
                f"(S, cross-batch, printed) {len(mism)} byte mismatches, max |ULP| {call_ulp}, "
                f"A1 {call_a1} / A2 {call_a2}; (S-sb, asserted) single-incidence same-batch "
                f"byte mismatches {n_sb_tot}; multi-incidence same-batch (printed) "
                f"{n_sb_multi_tot}; multi rel vs owner-cell eval {[f'{x:.3e}' for x in multi_rel]}; "
                f"(R, printed) {r_rows if r_rows else 'none'}",
                flush=True,
            )
        for rel in multi_rel:
            assert rel <= CLOSED_FORM_TOL, (tag, k, multi_rel)
    if comm.rank == 0:
        print(
            f"[POST-8 3b-ii-e2a E] {tag} -n {comm.size}: {len(calls)} calls recorded "
            f"(predicted {PREDICTED_CALLS[tag]}); {n_single}/{n_total} probes "
            f"single-incidence across calls; (S, cross-batch, printed) site total "
            f"{site_mism} mismatches, max |ULP| {site_ulp}, A1 {site_a1} / A2 {site_a2}; "
            f"(S-sb) failing calls {sb_failures}",
            flush=True,
        )
    assert not sb_failures, (tag, sb_failures)


def _drive_e2(tag):
    comm = MPI.COMM_WORLD
    assert comm.size == EXAMPLES[tag][1], f"{tag}: record width is -n {EXAMPLES[tag][1]}"
    calls = _run_spied(tag)
    _check_e2_s(tag, calls, comm)
    if tag in N0_CALL:
        _check_n0(tag, calls, comm, asserted=True)
    elif tag in N0_PRINTED_E2:
        _check_n0(tag, calls, comm, asserted=False)
    _write_e2_npz(tag, calls, comm)


def _drive_e3(tag):
    comm = MPI.COMM_WORLD
    calls = _run_spied(tag)
    _check_e2_s(tag, calls, comm)
    if tag in N0_CALL:
        _check_n0(tag, calls, comm, asserted=True)
    _check_o(tag, calls, comm)


@pytest.mark.skipif(IS_COMPLEX, reason="magnetostatic example: real build")
def test_mag2():
    _drive_e2("mag2")


@pytest.mark.skipif(IS_COMPLEX, reason="magnetostatic example: real build")
def test_mag4():
    _drive_e2("mag4")


@pytest.mark.skipif(IS_COMPLEX, reason="magnetostatic example: real build")
def test_mag6():
    _drive_e2("mag6")


@complex_only
def test_ports1():
    _drive_e3("ports1")


@complex_only
def test_ports12():
    _drive_e2("ports12")


@complex_only
def test_mri1():
    _drive_e2("mri1")


@complex_only
def test_ports1_owner_measure():
    """Item 568 Part B: ``ports:1`` under the owner rule, by injection; no width
    assertion.  The four 0.11-image record values and misses are printed; the
    example's pre-empted gates are asserted on the captured values."""
    comm = MPI.COMM_WORLD
    K = runpy.run_path(str(EXAMPLES["ports1"][0]), run_name="ports1_constants")
    calls, ladder_calls, s_calls, record_msg = _run_spied_measure("ports1")
    _check_e2_s("ports1", calls, comm)
    _check_n0("ports1", calls, comm, asserted=False)
    _check_o("ports1m", calls, comm)

    assert len(ladder_calls) >= 2 and len(s_calls) >= 1
    ladder = ladder_calls[0][2]
    blind = ladder_calls[1][2]
    s_args, s_kwargs, s_mat = s_calls[0]
    z_mat = np.asarray(s_args[0] if s_args else s_kwargs["z_matrix"])
    s_mat = np.asarray(s_mat)
    s_sym = float(np.linalg.norm(s_mat - s_mat.T) / np.linalg.norm(s_mat))
    s_spec = float(np.linalg.norm(s_mat, 2))
    vals = {
        "raw": float(ladder["raw"]),
        "corrected": float(ladder["corrected"]),
        "symmetry": s_sym,
        "spectral": s_spec,
    }
    rel, ab = K["_relative_miss"], K["_absolute_miss"]
    rec = {
        "raw": K["RECORDED_RAW_RATIO"],
        "corrected": K["RECORDED_CORRECTED_RATIO"],
        "symmetry": K["RECORDED_S_SYMMETRY_RESIDUAL"],
        "spectral": K["RECORDED_S_SPECTRAL_NORM"],
    }
    misses = {
        "raw": rel(vals["raw"], rec["raw"]),
        "corrected": rel(vals["corrected"], rec["corrected"]),
        "symmetry": ab(vals["symmetry"], rec["symmetry"]),
        "spectral": rel(vals["spectral"], rec["spectral"]),
    }
    sup = K["SUPERSEDED_V072_RECORD"]
    sup_miss = {n: rel(vals[n], sup[n]) for n in ("raw", "corrected")}
    band_r = K["REPRODUCTION_BAND_RELATIVE"]
    if comm.rank == 0:
        print(
            f"[POST-8 3b-ii-e2d B] ports1 -n {comm.size} owner-rule values: "
            + ", ".join(f"{n} {vals[n]:.10g}" for n in vals),
            flush=True,
        )
        print(
            f"[POST-8 3b-ii-e2d B] misses vs 0.11-image records (bands {band_r:.0e} relative, "
            f"{K['REPRODUCTION_BAND_SYMMETRY_ABSOLUTE']:.0e} absolute on symmetry): "
            f"raw {misses['raw']:.3e}, corrected {misses['corrected']:.3e}, "
            f"symmetry {misses['symmetry']:.3e} (abs), spectral {misses['spectral']:.3e}",
            flush=True,
        )
        print(f"[POST-8 3b-ii-e2d B] record_msg: {record_msg or 'records reproduced'}", flush=True)
        print(
            f"[POST-8 3b-ii-e2d B] pre-empted gates: superseded misses raw "
            f"{sup_miss['raw']:.3e} / corrected {sup_miss['corrected']:.3e} (> {band_r:.0e}); "
            f"corrected deviation {float(ladder['deviation']):+.4%} "
            f"(< {K['MUTUAL_TOLERANCE']}); blind deviation {float(blind['deviation']):+.4%} "
            f"(>= {K['MUTUAL_TOLERANCE']}); symmetry {s_sym:.4e} (< {K['S_SYMMETRY_BAND']}); "
            f"spectral {s_spec:.6f} (<= {K['S_SPECTRAL_NORM_CEILING']})",
            flush=True,
        )
    for n in ("raw", "corrected"):
        assert sup_miss[n] > band_r, (n, sup_miss[n])
    assert abs(float(ladder["deviation"])) < K["MUTUAL_TOLERANCE"]
    assert abs(float(blind["deviation"])) >= K["MUTUAL_TOLERANCE"]
    assert s_sym < K["S_SYMMETRY_BAND"]
    assert s_spec <= K["S_SPECTRAL_NORM_CEILING"]

    rec_npz = NPZ_DIR / "step3b-ii-e2d-ports1-records-n2.npz"
    if comm.size == 2:
        if comm.rank == 0:
            NPZ_DIR.mkdir(parents=True, exist_ok=True)
            np.savez(
                rec_npz,
                raw=vals["raw"],
                corrected=vals["corrected"],
                symmetry=vals["symmetry"],
                spectral=vals["spectral"],
                z=z_mat,
            )
            print(f"[POST-8 3b-ii-e2d B] -n 2 wrote {rec_npz}", flush=True)
        comm.Barrier()
        return
    comm.Barrier()
    if not rec_npz.exists():
        pytest.fail(f"{rec_npz} missing - run the -n 2 window first; never recompute")
    ref = np.load(rec_npz)
    if comm.rank == 0:
        spreads = {}
        for n in vals:
            r = float(ref[n])
            spreads[n] = abs(vals[n] - r) / max(abs(r), FIELD_SCALE_FLOOR)
        sym_abs = abs(vals["symmetry"] - float(ref["symmetry"]))
        z_spread = float(np.max(np.abs(z_mat - ref["z"])) / np.max(np.abs(ref["z"])))
        print(
            f"[POST-8 3b-ii-e2d B] -n {comm.size} vs -n 2 cross-width spread (printed beside "
            f"the Done-when 1e-9): relative raw {spreads['raw']:.3e}, corrected "
            f"{spreads['corrected']:.3e}, symmetry {spreads['symmetry']:.3e} "
            f"(absolute {sym_abs:.3e}), spectral {spreads['spectral']:.3e}; Z max-element "
            f"relative {z_spread:.3e}",
            flush=True,
        )
