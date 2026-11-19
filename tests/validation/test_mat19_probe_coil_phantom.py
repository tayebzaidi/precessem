"""`MAT-19` step 2 — the probe on `MAT-15` P1's coil-driven gel phantom (§9 item 253).

Brief: ``docs/planning/chunks/MAT-19.md`` (§9 item 253, chain P8).  Step 1's
ball-averaged probe (``thermal.probe``, unchanged) read on `MAT-15` step 1's
fixture (``_run_heated(return_fields=True)``: `ANS-2`'s 10 MHz loaded F-small
birdcage, 199 920 / 58 866 cells, drive k = 0 = P1, G = 1 W / P₁, P1 thermal
on the tag-3 submesh, Robin on every face) at `ANS-9`'s SPEC probe centres
A–E (loaded from the ``ans:9`` script, never restated), plus C′ (the σ_y image
of C), the whole-phantom ball Ω and the displaced control B⁺.  a = 5 mm.

Asserted: (A0) fixture identity, (M) every ball inside the phantom, snapshot 0
reads 0.0; (R0) the re-run is the gated run ≤ ``E_BAND``; (A) P_Ω = ∫θ/V_sub
≤ ``B_BAND`` at n = 0 … 40 and for the BE θ¹, V_Ω/V_sub ≤ ``B_BAND``; (U4)
bitwise rank agreement at ``-n 2``; (W) width invariance ≤ ``WIDTH_BAND`` at
``-n 4`` / ``-n 1`` against the ``-n 2`` npz (hash-guarded); (Σ) σ_y pair
C / C′ at t* ≤ ``C4_COVARIANCE_BAND``.  Negative controls (rule (e)): (U3c)
ghost-included Ω misses ``B_BAND``; (C1″) B⁺ vs B > ``WIDTH_BAND``; (f″) C vs
B > ``F_CONTROL_FLOOR``.  Printed: (P4) probe volume vs row-4 point, (P5)
placement sensitivity at B, (S) local calorimetric slope, series tables.

Scope: a record on the SPEC's (ragged-sheet, degree-1) fixture — no AED value,
no adjudication of row 4, not a real probe's sensing volume.

Re-pose, §9 item 258 (attempt 2 of 4): attempt 1 measured the tag-3 submesh
unpartitioned at ``-n 2`` (owned 58 866 / 0, ghosts 0 / 0,
``20261005T184012Z_MAT-19-step2-ghost-diag.log:1871–1873``), so (U3c) had
nothing to double-count.  Pre-registered applicability rule: (U3c) is asserted
iff G_w := Σ_ranks submesh ``num_ghosts`` > 0; the "submesh partition" reading
of (A) / (U4) holds iff O_w := #{ranks with submesh ``size_local`` > 0} ≥ 2.
(U4) is asserted at every width > 1; (W) at ``-n 4`` / ``-n 1`` only (printed
at any other width).  The (S) helper pre-filters cells by bounding radius as
``probe._geometry`` does (E6; a dropped cell has v_K = 0 exactly).
"""

from __future__ import annotations

import hashlib
import importlib.util
import os
import sys
import time
from pathlib import Path

import numpy as np
import pytest
import ufl
from dolfinx import fem
from mpi4py import MPI

from fem_em_solver.thermal import ProbeBall, ball_average  # noqa: F401
from fem_em_solver.thermal import probe as _probe_mod
from fem_em_solver.thermal import probe_series
from fem_em_solver.thermal.probe import _ball_average_on_cells, _lattice, _owned_cells
from tests.complex_mode import complex_only
from tests.mesh.test_birdcage_port_tags import PHANTOM_HEIGHT
from tests.validation.test_mat15_coil_driven_heating_phantom import (
    B_BAND,
    C4_COVARIANCE_BAND,
    DT_1,
    E_BAND,
    F_CONTROL_FLOOR,
    N_MAIN,
    RHO_C,
    T_INF,
    T_STAR,
    V_BAND,
    V_SUB_RECORD_M3,
    A,
    _allsum,
    _be_one,
    _run_heated,
    _say,
    _transient,
)
from tests.validation.test_mat19_virtual_temperature_probe import WIDTH_BAND

_REPO = Path(__file__).resolve().parents[2]
_ANS9 = (
    _REPO
    / "examples"
    / "ansys_benchmarks"
    / "ans9_coil_driven_gel_heating_10MHz"
    / "09_coil_driven_gel_heating_10MHz.py"
)

A_PROBE = 5.0e-3  # 2 h_p, a fixture parameter, not a fibre-optic sensing volume
OMEGA_RADIUS = 0.1  # whole-phantom ball (farthest vertex 0.05 m)
B_SHIFT = 1.0e-4  # B⁺: B moved 0.1 mm radially
N_SUB = 58866  # 20261001T123309Z_MAT-15-step1-W1.log:1926
DT = T_STAR / N_MAIN
NPZ = "/workspace/logs/mat19/step2-n2.npz"
NPZ4 = "/workspace/logs/mat19/step2-n4.npz"  # printed pairwise -n 1 vs -n 4 only

# Records printed beside (never asserted): ANS-9 W1 log
# 20261004T142653Z_ANS-9-step1a-W1.log:1943 (row 4 points) and :1936 (row 2).
ROW4_POINT_K = {"A": 0.905798, "B": 3.564854, "C": 0.685305, "D": 0.697444, "E": 0.958639}
ROW2_MEAN_K = {10: 0.318096, 20: 0.614315, 30: 0.892104, 40: 1.153314}
S2_OMEGA_RECORD = 5.754e-6  # MAT-15 (S2) on P1, row 3


def _load_spec_probes():
    name = "ans9_coil_driven_gel_heating_10mhz_spec"
    if name in sys.modules:
        mod = sys.modules[name]
    else:
        spec = importlib.util.spec_from_file_location(name, _ANS9)
        mod = importlib.util.module_from_spec(spec)
        sys.modules[name] = mod
        spec.loader.exec_module(mod)
    return dict(mod.SPEC_THERMAL["probes"])


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _probe_set():
    spec = _load_spec_probes()
    names = list(spec)  # A, B, C, D, E
    c = spec["C"]
    b = np.asarray(spec["B"], dtype=float)
    bp = b * (np.linalg.norm(b) + B_SHIFT) / np.linalg.norm(b)
    centres = {k: tuple(float(v) for v in spec[k]) for k in names}
    centres["C'"] = (float(c[0]), -float(c[1]), float(c[2]))
    centres["Omega"] = (0.0, 0.0, 0.0)
    centres["B+"] = tuple(float(v) for v in bp)
    order = ["A", "B", "C", "C'", "D", "E", "Omega", "B+"]
    radii = {k: (OMEGA_RADIUS if k == "Omega" else A_PROBE) for k in order}
    return names, order, [ProbeBall(centres[k], radii[k]) for k in order], centres


def _lattice_q_mean(V, q, probe, lam, chunk=4096):
    """(S) helper: Σ_K q_K v_K(B) / Σ_K v_K(B) over owned cells, one allreduce on every rank."""
    xdof = np.real(V.tabulate_dof_coordinates()).astype(np.float64)
    owned = _owned_cells(V)
    c = np.asarray(probe.centre, dtype=np.float64)
    a = probe.radius
    qarr = np.real(q.x.array)
    qdm = q.function_space.dofmap.list
    if (
        owned.size
    ):  # E6: bounding-radius pre-filter, as probe._geometry (v_K = 0 for a dropped cell)
        Xa = xdof[V.dofmap.list[owned]]
        cen = Xa.mean(axis=1)
        rho = np.sqrt(np.max(np.sum((Xa - cen[:, None, :]) ** 2, axis=2), axis=1))
        owned = owned[np.linalg.norm(cen - c, axis=1) - rho <= a]
    num, den = 0.0, 0.0
    for s in range(0, owned.size, chunk):
        cells = owned[s : s + chunk]
        X = xdof[V.dofmap.list[cells]]
        vol = np.abs(np.linalg.det(X[:, 1:, :] - X[:, :1, :])) / 6.0
        P = np.einsum("pk,mkd->mpd", lam, X)
        inside = np.sum((P - c) ** 2, axis=2) <= a * a
        vK = vol * np.count_nonzero(inside, axis=1) / lam.shape[0]
        qK = qarr[qdm[cells, 0]]
        num += float(np.sum(qK * vK))
        den += float(np.sum(vK))
    tot = V.mesh.comm.allreduce(np.array([num, den]), op=MPI.SUM)
    return float(tot[0] / tot[1]), float(tot[1])


@pytest.fixture(scope="module")
def run():
    comm = MPI.COMM_WORLD
    clock = {}
    t0 = time.perf_counter()
    out = _run_heated(return_fields=True)
    clock["_run_heated"] = time.perf_counter() - t0
    sub, ft = out["sub"], out["ft"]
    q_h = out["q_h"][0]
    v_sub = out["v_sub"]
    n_sub = sub.topology.index_map(3).size_global

    t1 = time.perf_counter()
    res = _transient(
        sub,
        ft,
        q_h,
        dt=DT,
        n_steps=N_MAIN,
        theta=0.5,
        startup_half_steps=2,
        snapshot_steps=range(N_MAIN + 1),
    )
    res_be = _be_one(sub, ft, q_h)
    clock["snapshot re-run + BE"] = time.perf_counter() - t1

    # (R0) in ans:9's form.
    d = res.snapshots[N_MAIN] - out["T_star"][0]
    th_star = out["T_star"][0] - T_INF
    r0 = float(
        np.sqrt(np.real(_allsum(comm, ufl.inner(d, d) * ufl.dx(domain=sub))))
        / np.sqrt(np.real(_allsum(comm, ufl.inner(th_star, th_star) * ufl.dx(domain=sub))))
    )

    V = res.snapshots[N_MAIN].function_space
    th = {}
    for n in range(N_MAIN + 1):
        f = fem.Function(V)
        f.x.array[:] = res.snapshots[n].x.array - T_INF
        f.x.scatter_forward()
        th[n] = f
    th_be = fem.Function(V)
    th_be.x.array[:] = res_be.temperature.x.array - T_INF
    th_be.x.scatter_forward()

    names, order, probes, centres = _probe_set()
    t2 = time.perf_counter()
    lam = _lattice(3)  # the same lattice probe._geometry uses at levels=3, for (S)
    s = probe_series(th, DT, probes)
    s_be = probe_series({1: th_be}, DT_1, probes)
    clock["probe geometry + series sums"] = time.perf_counter() - t2

    t3 = time.perf_counter()
    qbar = {}
    for k, p in zip(order, probes):
        qbar[k], _vol = _lattice_q_mean(V, q_h, p, lam)
    clock["(S) helper"] = time.perf_counter() - t3

    means = np.array(
        [
            float(np.real(_allsum(comm, th[n] * ufl.dx(domain=sub)))) / v_sub
            for n in range(N_MAIN + 1)
        ]
    )
    mean_be = float(np.real(_allsum(comm, th_be * ufl.dx(domain=sub)))) / v_sub

    # E2: the submesh partition table, every rank, every width.
    imap = sub.topology.index_map(3)
    part = comm.allgather((int(imap.size_local), int(imap.num_ghosts)))
    o_w = sum(1 for so, _g in part if so > 0)
    g_w = sum(g for _so, g in part)
    _say(
        f"MAT19-2 (part) -n {comm.size}: sub owned {[p[0] for p in part]} ghosts "
        f"{[p[1] for p in part]} O_w {o_w} G_w {g_w}; partition exercised: {o_w >= 2}"
    )

    # (U3c) planted ghost control on Ω (computed always at size > 1; asserted iff G_w > 0).
    ghost = None
    if comm.size > 1:
        allc = np.arange(imap.size_local + imap.num_ghosts, dtype=np.int32)
        gval, gvol = _ball_average_on_cells(th[N_MAIN], probes[order.index("Omega")], allc)
        ghost = (abs(gvol / v_sub - 1.0), abs(gval / means[N_MAIN] - 1.0))

    # (U4) in-run bitwise agreement across ranks at n = 40.
    allrows = comm.allgather(s.values[N_MAIN].copy())
    u4 = all(np.array_equal(allrows[0], r) for r in allrows)

    me = Path(__file__).resolve()
    hashes = (_sha(me), _sha(_probe_mod.__file__))
    clock["fixture total"] = time.perf_counter() - t0
    return dict(
        out=out,
        n_sub=n_sub,
        v_sub=v_sub,
        r0=r0,
        s=s,
        s_be=s_be,
        qbar=qbar,
        means=means,
        mean_be=mean_be,
        ghost=ghost,
        u4=u4,
        names=names,
        order=order,
        centres=centres,
        hashes=hashes,
        clock=clock,
        comm=comm,
        o_w=o_w,
        g_w=g_w,
    )


def _idx(run, k):
    return run["order"].index(k)


@complex_only
def test_mat19_step2_probe_coil_phantom(run):  # noqa: C901
    comm = run["comm"]
    out, s, s_be = run["out"], run["s"], run["s_be"]
    order, names, centres = run["order"], run["names"], run["centres"]
    v = s.values
    j = {k: _idx(run, k) for k in order}
    a_name = "Omega"

    _say(
        f"\nMAT19-2 drive by_k[0] = {out['by_k'][0]}; P_1 = {out['P'][0]:.12e} W  "
        f"G = {out['G']:.12e} 1/W; width -n {comm.size}"
    )
    err_V = abs(run["v_sub"] / V_SUB_RECORD_M3 - 1.0)
    _say(
        f"MAT19-2 (A0) sub cells = {run['n_sub']} (record {N_SUB}); V_sub = {run['v_sub']:.12e} "
        f"rel {err_V:.3e} (bound {V_BAND:.0e}); dt = {DT:.4f} s; steps {list(s.steps)[:2]}..{s.steps[-1]}"
    )
    _say(
        "MAT19-2 centres [m]: "
        + "  ".join(
            f"{k} ({centres[k][0]:+.4f},{centres[k][1]:+.4f},{centres[k][2]:+.4f})" for k in order
        )
    )
    _say(
        "MAT19-2 ball V_h/(4 pi a^3/3) - 1: "
        + "  ".join(
            f"{k} {s.volumes[j[k]] / s.nominal_volumes[j[k]] - 1.0:+.3e}"
            for k in order
            if k != a_name
        )
        + " (predicted |.| <= 5e-3)"
    )
    m_ok = all(
        np.hypot(centres[k][0], centres[k][1]) + A_PROBE <= A
        and abs(centres[k][2]) + A_PROBE <= PHANTOM_HEIGHT / 2.0
        for k in order
        if k != a_name
    )
    _say(f"MAT19-2 (M) every probe ball inside the phantom: {m_ok}")
    _say(
        f"MAT19-2 (R0) |theta40_rerun - theta*_gated| / |theta*_gated| = {run['r0']:.3e} "
        f"(bound {E_BAND:.0e}, predicted <= 1e-14)"
    )

    # (A)
    PO = v[:, j[a_name]]
    relA = np.abs(PO[1:] - run["means"][1:]) / np.abs(run["means"][1:])
    relA_be = abs(s_be.values[0, j[a_name]] - run["mean_be"]) / abs(run["mean_be"])
    volA = abs(s.volumes[j[a_name]] / run["v_sub"] - 1.0)
    _say(
        f"MAT19-2 (A) P_Omega vs int theta/V_sub: max over n=1..40 {np.max(relA):.3e}, BE theta^1 "
        f"{relA_be:.3e}; V_Omega/V_sub - 1 = {volA:.3e} (bound {B_BAND:.0e}, predicted <= 1e-14); "
        f"n=0: P_Omega {PO[0]!r} mean {run['means'][0]!r}"
    )
    _say(
        "MAT19-2 (A) SPEC times P_Omega vs row 2: "
        + "  ".join(f"n={n} {PO[n]:.6f} (row2 {ROW2_MEAN_K[n]:.6f})" for n in (10, 20, 30, 40))
    )

    # Series table
    _say("MAT19-2 series theta_bar [K] (n: " + " ".join(order) + "):")
    for n in (1, 10, 20, 30, 40):
        _say(f"MAT19-2   n={n:2d} t={n * DT:7.2f}  " + " ".join(f"{v[n, j[k]]:.6e}" for k in order))
    _say(
        f"MAT19-2 max_imag_ratio {s.max_imag_ratio:.3e} (predicted 0); BE {s_be.max_imag_ratio:.3e}"
    )

    # (P4)
    _say(
        "MAT19-2 (P4) ball(t*) / row-4 point - 1: "
        + "  ".join(f"{k} {v[N_MAIN, j[k]] / ROW4_POINT_K[k] - 1.0:+.4e}" for k in names)
        + "  (predicted B in [-8 %, 0], A/C/D/E within +-3 %)"
    )
    # (P5) / (C1″)
    c1 = abs(v[N_MAIN, j["B+"]] - v[N_MAIN, j["B"]]) / v[N_MAIN, j["B"]]
    _say(
        f"MAT19-2 (P5)/(C1'') |B+ - B|/B at t* = {c1:.4e} (predicted ≈ 3.3e-3, range [1e-3, 1e-2]; "
        f"asserted > {WIDTH_BAND:.0e})"
    )
    # (Σ)
    sig = abs(v[N_MAIN, j["C"]] - v[N_MAIN, j["C'"]]) / v[N_MAIN, j["C"]]
    qc, qcp = run["qbar"]["C"], run["qbar"]["C'"]
    _say(
        f"MAT19-2 (Sigma) |C - C'|/C at t* = {sig:.4e} (bound {C4_COVARIANCE_BAND:.0e}, predicted "
        f"<= 5e-3); q_h ball pair C {qc:.6e} C' {qcp:.6e} rel {abs(qc - qcp) / abs(qc):.4e} (printed)"
    )
    jc, jcp = j["C"], j["C'"]
    _say(
        "MAT19-2 (Sigma) C/C' series rel: "
        + " ".join(
            f"n={n}:{abs(v[n, jc] - v[n, jcp]) / v[n, jc]:.3e}" for n in (1, 5, 10, 20, 30, 40)
        )
    )
    # (f″)
    ff = abs(v[N_MAIN, j["C"]] - v[N_MAIN, j["B"]]) / v[N_MAIN, j["B"]]
    _say(
        f"MAT19-2 (f'') |C - B|/B at t* = {ff:.4f} (predicted ≈ 0.81, asserted > {F_CONTROL_FLOOR:.0e})"
    )
    # (S)
    _say(
        "MAT19-2 (S) rho c theta_bar(theta^1)/(dt_1 q_bar) - 1: "
        + "  ".join(
            f"{k} {RHO_C * s_be.values[0, j[k]] / (DT_1 * run['qbar'][k]) - 1.0:+.4e}"
            for k in order
        )
        + f"  (printed; Omega record (S2) {S2_OMEGA_RECORD:.3e})"
    )
    # (U3c)
    if run["ghost"] is not None:
        _say(
            f"MAT19-2 (U3c) ghost-included Omega: volume miss {run['ghost'][0]:.4e}, value miss "
            f"{run['ghost'][1]:.4e} (asserted volume miss > {B_BAND:.0e} iff G_w > 0; predicted "
            f"[2e-2, 2e-1] at -n 4; sharp G_w/{N_SUB} = {run['g_w'] / N_SUB:.4e} within x[0.5, 2])"
        )
        if run["g_w"] == 0:
            _say(f"MAT19-2 (U3c) inapplicable at -n {comm.size}: G_w = 0")
    _say(f"MAT19-2 (U4) bitwise rank agreement at n=40: {run['u4']} (-n {comm.size})")

    # (W)
    w = None
    if comm.size == 2:
        if comm.rank == 0:
            os.makedirs(os.path.dirname(NPZ), exist_ok=True)
            np.savez(
                NPZ,
                series=v,
                volumes=s.volumes,
                cells=np.array([run["n_sub"]]),
                sha_module=run["hashes"][0],
                sha_probe=run["hashes"][1],
            )
            print(
                f"MAT19-2 (W) wrote {NPZ} (sha module {run['hashes'][0][:12]}, probe "
                f"{run['hashes'][1][:12]})",
                flush=True,
            )
        comm.Barrier()
    else:
        if not os.path.exists(NPZ):
            pytest.fail(f"(W) {NPZ} missing: run W1 (-n 2) first; never recompute")
        z = np.load(NPZ)
        if str(z["sha_module"]) != run["hashes"][0] or str(z["sha_probe"]) != run["hashes"][1]:
            pytest.fail("(W) npz hash differs from this module / thermal/probe.py: re-run W1")
        cells_eq = int(z["cells"][0]) == run["n_sub"]
        rel = np.abs(v[1:] / z["series"][1:] - 1.0)
        w = float(np.max(rel))
        wvol = float(np.max(np.abs(s.volumes / z["volumes"] - 1.0)))
        _say(
            f"MAT19-2 (W) -n {comm.size} vs -n 2: cells equal {cells_eq}; series max rel {w:.3e} "
            f"(band {WIDTH_BAND:.0e}, predicted <= 1e-10); per probe max: "
            + "  ".join(f"{k} {np.max(rel[:, j[k]]):.2e}" for k in order)
            + f"; volumes {wvol:.3e}"
        )
        if comm.size not in (4, 1):  # E5: (W) printed only at other widths
            w = None
        if comm.size == 4:
            if comm.rank == 0:
                np.savez(NPZ4, series=v, sha_module=run["hashes"][0], sha_probe=run["hashes"][1])
            comm.Barrier()
        elif comm.size == 1 and os.path.exists(NPZ4):
            z4 = np.load(NPZ4)
            if (
                str(z4["sha_module"]) == run["hashes"][0]
                and str(z4["sha_probe"]) == run["hashes"][1]
            ):
                _say(
                    f"MAT19-2 (W) pairwise -n 1 vs -n 4: series max rel "
                    f"{float(np.max(np.abs(v[1:] / z4['series'][1:] - 1.0))):.3e} (printed)"
                )
    _say("MAT19-2 phase clocks: " + "  ".join(f"{k} {t:.1f} s" for k, t in run["clock"].items()))

    assert out["by_k"][0] == "P1", f"(A0) by_k[0] = {out['by_k'][0]}"
    assert run["n_sub"] == N_SUB, f"(A0) sub cells {run['n_sub']} != {N_SUB}"
    assert err_V <= V_BAND, f"(A0) V_sub rel {err_V:.3e}"
    assert m_ok, "(M) a probe ball leaves the phantom"
    assert np.all(v[0] == 0.0), "(A0) snapshot 0 does not read 0.0"
    assert run["r0"] <= E_BAND, f"(R0) {run['r0']:.3e}"
    assert PO[0] == 0.0 and run["means"][0] == 0.0
    assert np.max(relA) <= B_BAND, f"(A) {np.max(relA):.3e}"
    assert relA_be <= B_BAND, f"(A) BE {relA_be:.3e}"
    assert volA <= B_BAND, f"(A) volume {volA:.3e}"
    if comm.size > 1:
        assert run["u4"], "(U4) ranks disagree bitwise"
    if run["ghost"] is not None and run["g_w"] > 0:
        assert (
            run["ghost"][0] > B_BAND
        ), f"(U3c) ghost-included volume did not miss: {run['ghost'][0]:.3e}"
    assert sig <= C4_COVARIANCE_BAND, f"(Sigma) {sig:.4e}"
    assert c1 > WIDTH_BAND, f"(C1'') {c1:.3e}"
    assert ff > F_CONTROL_FLOOR, f"(f'') {ff:.4f}"
    if w is not None:
        assert cells_eq, "(W) cells differ"
        assert w <= WIDTH_BAND, f"(W) {w:.3e} > {WIDTH_BAND:.0e}"
