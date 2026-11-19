"""`TH-20` step 2a (§10 chain H13b, §9 item 350) — the four degree-2 quarter
sectors of the F-human copper birdcage, priced, with the D2 reassembly R saved.

Selected by environment only, never by ``-k``:
``FEM_EM_TH20_STEP2=1`` runs test 1 (the sector solves) on the sectors named in
``FEM_EM_TH20_STEP2_SECTORS`` (comma list from EE, EM, ME, MM; unset = all);
``FEM_EM_TH20_STEP2_ASSEMBLE=1`` runs test 2 (no mesh, no solve: reads the
four saved sector npz files, rebuilds R and saves it for step 2b).

Fixture: degree 2, 64 MHz, `TH-14`'s Leontovich copper wall (σ = 5.8e7) on tag
401, 50 Ω exact ring sheets (κ off), the 1.6 mm hole quarter of step 0a (built
after the direct full mesh, `GEO-23`). Sector (b_x, b_y): BC on x = 0 (501)
and y = 0 (502); E = PEC (pinned), M = PMC (natural).

Anchors (asserted, bands imported): (P0) quarter cell / cut-facet records;
(Deg) solver.degree == 2 and unknowns == 2E + 2F; (P2) P_surf > 0 on every
sector operator; (I1) per-sector reciprocity ≤ ``RECIPROCITY_BAND``; (I2)
per-sector σ_max ≤ ``MATRIX_PASSIVITY_CEILING``; in test 2, (P1) constrained-DOF
ordering from the saved counts and (P3) the image bijection / orientation.

(C2) the collapsed-sectors control is PRINTED, NOT ASSERTED (2026-10-09 18:00
review ruling (2) on item 350: its only backing measurement is degree 1, and a
change of element order is a different comparison under rule (e)). The brief's
``FEM_EM_TH20_STEP2_CONTROL`` assert path is therefore not implemented.

Scope: no reassembly identity (A) at degree 2 (step 2b, `xl`), no comparison
with H15b's σ = 800 solid matrix, no C_mode / tuning / B1+ / SAR, no absolute S.
"""

from __future__ import annotations

import gc
import os
import time

import numpy as np
import pytest
from mpi4py import MPI

from tests.complex_mode import complex_only
from tests.validation.test_th20_quarter_solves import (
    C2_FLOOR,
    COPPER_SIGMA,
    GROUP,
    QUARTER_CELL_RECORD,
    SECTORS,
    TAU,
    _make_drive,
    _owned_tag_count,
    _pi,
    _recip,
    _reconstruct,
    _sigma_max,
    _specs,
)

STEP2_ENV = "FEM_EM_TH20_STEP2"
SECTORS_ENV = "FEM_EM_TH20_STEP2_SECTORS"
ASSEMBLE_ENV = "FEM_EM_TH20_STEP2_ASSEMBLE"
DEGREE_2 = 2
OUT_DIR = "/workspace/logs/th20"
NAME = {s: f"{s[0]}{s[1]}" for s in SECTORS}
BY_NAME = {v: k for k, v in NAME.items()}
PEC_OF = {("E", "E"): (1, 501, 502), ("E", "M"): (1, 501), ("M", "E"): (1, 502), ("M", "M"): (1,)}
# Records, printed beside the degree-2 readings, never compared (brief):
# step 1 W2 degree-1 S^s(P17, P17) (`20261005T020907Z_TH-20-step1-W2.log:8307/8316/8325/8334`)
D1_S1717 = {
    ("E", "E"): 0.087643 + 0.432975j,
    ("E", "M"): 0.584504 + 0.435851j,
    ("M", "E"): 0.275464 + 0.258848j,
    ("M", "M"): 0.586144 + 0.438849j,
}
# step 1 W2 degree-1 R(P17, P17) (`:8381`) and degree-1 constrained DOFs (`:8342`)
D1_R1717 = 0.383439 + 0.391631j
D1_DOFS = {("E", "E"): 10_662, ("E", "M"): 7_874, ("M", "E"): 7_853, ("M", "M"): 5_035}
# step 1's degree-1 C2 (`:8385`); degree-2 prediction band (brief), printed only
D1_C2 = 1.1346e-1
C2_PREDICTED = (0.08, 0.15)
# H15b's full σ = 800 peak summed RSS at -n 16 (`docs/testing/xl-ledger.md:55`), printed
H15B_FULL_RSS_GIB = 103.95
# Pre-stated per-sector price bracket at -n 8 (brief), printed only
PRED_UNKNOWNS = (0.66e6, 0.72e6)
PRED_FACTOR_S = (60.0, 130.0)
PRED_SECTOR_S = (75.0, 160.0)
PRED_RSS_GIB = (15.0, 30.0)


def _path(name):
    return os.path.join(OUT_DIR, f"th20-step2-d2-{name}.npz")


def _sectors_from_env():
    raw = os.environ.get(SECTORS_ENV, "").strip()
    if not raw:
        return list(SECTORS)
    names = [t.strip() for t in raw.split(",") if t.strip()]
    bad = [t for t in names if t not in BY_NAME]
    if bad or not names or len(set(names)) != len(names):
        raise ValueError(f"{SECTORS_ENV} must be a comma list from EE, EM, ME, MM; got {raw!r}")
    return [BY_NAME[t] for t in names]


def _constrained_dofs_degree(msh, facet_tags, pec_tags, degree):
    """`test_th16_quarter_birdcage._constrained_dofs` on an N1curl-``degree`` space."""
    from dolfinx import fem

    V = fem.functionspace(msh, ("N1curl", degree))
    found = [np.asarray(facet_tags.find(int(t)), dtype=np.int32) for t in pec_tags]
    facets = np.unique(np.concatenate(found)).astype(np.int32)
    dofs = fem.locate_dofs_topological(V, msh.topology.dim - 1, facets)
    n_owned = V.dofmap.index_map.size_local * V.dofmap.index_map_bs
    return msh.comm.allreduce(int(np.count_nonzero(dofs < n_owned)), op=MPI.SUM)


def _preimage_and_p3(kept, ordinals):
    """Step 1's (P3) block (`test_th20_quarter_solves.py:303–321`), pure arithmetic."""
    from tests.mesh.test_birdcage_ring_gaps_scaleup import _ring_gap_frame
    from tests.validation.test_birdcage_f_human_rung import LEG_COUNT
    from tests.validation.test_th20_mirror_built_full_mesh import OPS

    preimage = {}
    for o in kept:
        for g in GROUP:
            preimage.setdefault(_pi(g, o), []).append((o, g))
    bij = set(preimage) == set(ordinals) and all(len(v) == 1 for v in preimage.values())
    preimage = {a: v[0] for a, v in preimage.items()}
    worst_tau = 0.0
    for o in kept:
        d_o, _c = _ring_gap_frame(o, LEG_COUNT)
        for g in ("sx", "sy", "sxsy"):
            d_img, _c2 = _ring_gap_frame(_pi(g, o), LEG_COUNT)
            worst_tau = max(worst_tau, abs(float(np.dot(OPS[g] * d_o, d_img)) - TAU[g]))
    return preimage, bij, worst_tau


@complex_only
def test_th20_step2_degree2_sectors():
    if os.environ.get(STEP2_ENV, "").strip() != "1":
        pytest.skip(f"{STEP2_ENV} unset: TH-20 step 2a sectors are selected by environment only")
    from fem_em_solver.core import TimeHarmonicSolver
    from fem_em_solver.core.cavity import surface_resistance_ohm
    from fem_em_solver.io.mesh import _interface_facet_tags
    from tests.mesh.test_birdcage_port_sheets import PORT_LOWER, PORT_UPPER, SHEET_IFACE
    from tests.mesh.test_birdcage_port_terminals import _interface_area_or_zero
    from tests.validation.test_port_birdcage_ring_matrix import MATRIX_PASSIVITY_CEILING
    from tests.validation.test_port_lumped_sheet_sweep import RECIPROCITY_BAND
    from tests.validation.test_th14_birdcage_copper import _problem
    from tests.validation.test_th16_quarter_birdcage import RECIPROCITY_Q_BAND, _constrained_dofs
    from tests.validation.test_th20_mirror_built_full_mesh import CUT_FACES_X0, CUT_FACES_Y0
    from tests.validation.test_th20_quarter_f_human_mesh import KEPT_RING, _meshes
    from tests.validation.test_wf7_f_human_b1_quadrature import FREQUENCY_HZ
    from tests.validation.test_wf12_f_human_hole import _rss_gib

    comm = MPI.COMM_WORLD
    t_window = time.perf_counter()
    run = _sectors_from_env()  # raises before the build on a bad list
    fails = []

    def say(msg):
        if comm.rank == 0:
            print(f"[TH-20 step2a] {msg}", flush=True)

    # Stale-file refusal before the build (rank 0 decides, broadcast).
    stale = (
        [_path(NAME[s]) for s in run if os.path.exists(_path(NAME[s]))] if comm.rank == 0 else None
    )
    stale = comm.bcast(stale, root=0)
    if stale:
        raise FileExistsError(f"stale step-2a npz present (refusing to mix attempts): {stale}")
    if comm.rank == 0:
        os.makedirs(OUT_DIR, exist_ok=True)

    say(
        f"start: -n {comm.size}, degree {DEGREE_2}, f {FREQUENCY_HZ:.3e} Hz, sectors "
        f"{[NAME[s] for s in run]}; C2 printed in test 2, never asserted (review ruling 2)"
    )
    omega = 2.0 * np.pi * FREQUENCY_HZ
    z_s = (1.0 + 1.0j) * surface_resistance_ohm(omega, COPPER_SIGMA)
    kept = list(KEPT_RING)

    # ---- meshes (full first, then the quarter: GEO-23) ---------------------
    rss0 = _rss_gib(comm)[2]
    m = _meshes()
    full_d, q = m["full"], m["quarter"]
    msh_q, cells_q, facets_q = q["msh"], q["cells"], q["facets"]
    n501 = _owned_tag_count(msh_q, facets_q, 501, comm)
    n502 = _owned_tag_count(msh_q, facets_q, 502, comm)
    # Edges / faces exist only once created (W1 attempt 1, 20261010T003223Z: RuntimeError
    # "Missing IndexMap in Topology" here; W1 re-run under rule (i)).
    msh_q.topology.create_entities(1)
    msh_q.topology.create_entities(2)
    n_edges = int(msh_q.topology.index_map(1).size_global)
    n_faces = int(msh_q.topology.index_map(2).size_global)
    say(
        f"(P0) quarter cells {q['n_cells']} (record {QUARTER_CELL_RECORD}, ASSERTED); owned "
        f"501/502 {n501}/{n502} (records {CUT_FACES_X0}/{CUT_FACES_Y0}, ASSERTED); direct full "
        f"build {float(full_d['t_build']):.1f} s; t={time.perf_counter() - t_window:.1f} s"
    )
    if q["n_cells"] != QUARTER_CELL_RECORD:
        fails.append(f"P0 quarter cells {q['n_cells']}")
    if (n501, n502) != (CUT_FACES_X0, CUT_FACES_Y0):
        fails.append(f"P0 cut facets {(n501, n502)}")
    assert not fails, fails  # P0 off its record stops the window before any solve

    tdim = msh_q.topology.dim
    msh_q.topology.create_connectivity(tdim - 1, tdim)
    msh_q.topology.create_entity_permutations()
    tags_q = _interface_facet_tags(
        msh_q, cells_q, {SHEET_IFACE + o: (PORT_LOWER + o, PORT_UPPER + o) for o in kept}
    )
    chord = float(q["diag"]["ring_port_layout"]["ring_port_gap_chord_m"])
    width_of = {
        o: _interface_area_or_zero(msh_q, tags_q, SHEET_IFACE + o, comm) / chord for o in kept
    }
    specs_q = _specs(kept, width_of, chord)
    drive_q = _make_drive(msh_q, tags_q, facets_q, specs_q, omega, z_s, comm)

    t_before_first = time.perf_counter() - t_window
    _mine, per, tot = _rss_gib(comm)
    say(
        f"phase mesh+tags: t {t_before_first:.1f} s, RSS summed {tot:.2f} GiB (start "
        f"{rss0:.2f}), per rank {', '.join(f'{r:.2f}' for r in (per or []))}"
    )

    t_sector_max = 0.0
    deg_checked = False
    for sector in run:
        name = NAME[sector]
        r_before = _rss_gib(comm)[2]
        t0 = time.perf_counter()
        pec = PEC_OF[sector]
        dofs_d1 = _constrained_dofs(msh_q, facets_q, pec)
        dofs_d2 = _constrained_dofs_degree(msh_q, facets_q, pec, DEGREE_2)
        problem = _problem(msh_q, cells_q, FREQUENCY_HZ, facets_q, pec)
        solver = TimeHarmonicSolver(problem, degree=DEGREE_2)
        v_space = solver.function_space()
        n_dofs = int(v_space.dofmap.index_map.size_global * v_space.dofmap.index_map_bs)
        if not deg_checked:
            say(
                f"(Deg) solver.degree {solver.degree}; global unknowns {n_dofs} vs 2E+2F = "
                f"2*{n_edges}+2*{n_faces} = {2 * n_edges + 2 * n_faces} (ASSERTED exact; "
                f"predicted {PRED_UNKNOWNS[0]:.2e}-{PRED_UNKNOWNS[1]:.2e})"
            )
            say(
                f"(Deg) control (PREDICTED, printed; rule (e)): a degree-1 space would read E = "
                f"{n_edges} unknowns, not 2E+2F (the `WF-15` step 1 degree-1 default trap)"
            )
            assert solver.degree == DEGREE_2, solver.degree
            assert n_dofs == 2 * n_edges + 2 * n_faces, (n_dofs, n_edges, n_faces)
            deg_checked = True
        say(
            f"sector {name}: pec {pec}; constrained DOFs degree-2 {dofs_d2} (degree-1 "
            f"{dofs_d1}; step 1 read {D1_DOFS[sector]})"
        )
        cols = {}
        for k, o in enumerate(kept):
            d = drive_q(solver, f"P{o}", held=(k > 0))
            cols[o] = d
            _mine, per, tot = _rss_gib(comm)
            say(
                f"  sector {name} drive {k + 1}/8 P{o}: {'held' if k else 'factor'} "
                f"{d['t_solve']:.2f} + post {d['t_post']:.2f} s; RSS summed {tot:.2f} GiB; "
                f"t={time.perf_counter() - t_window:.1f} s"
            )
        del solver, problem, v_space
        gc.collect()
        s8 = np.array([[cols[oj]["s"][f"P{oi}"] for oj in kept] for oi in kept], dtype=complex)
        t_sec = time.perf_counter() - t0
        t_sector_max = max(t_sector_max, t_sec)
        f0 = cols[kept[0]]
        p_surf = float(f0["p_surf"])
        t_factor = f0["t_solve"] + f0["t_post"]
        held = [cols[o]["t_solve"] + cols[o]["t_post"] for o in kept[1:]]
        rec, rec_max = _recip(s8)
        smax = _sigma_max(s8)
        _mine, per, tot = _rss_gib(comm)
        drss = tot - r_before
        say(
            f"sector {name}: factorising drive {f0['t_solve']:.2f} + post {f0['t_post']:.2f} = "
            f"{t_factor:.2f} s (predicted {PRED_FACTOR_S[0]:.0f}-{PRED_FACTOR_S[1]:.0f}); held "
            f"max {max(held):.2f} s mean {np.mean(held):.2f} s; sector total {t_sec:.1f} s "
            f"(predicted {PRED_SECTOR_S[0]:.0f}-{PRED_SECTOR_S[1]:.0f}); P_surf {p_surf:.9e} W "
            f"(> 0, ASSERTED)"
        )
        say(
            f"sector {name}: (I1) {rec:.3e} (band {RECIPROCITY_BAND:g}, ASSERTED; entrywise max "
            f"{rec_max:.3e} beside {RECIPROCITY_Q_BAND:g}, printed); (I2) sigma_max "
            f"{smax:.12f} (ceiling {MATRIX_PASSIVITY_CEILING}, ASSERTED)"
        )
        say(
            f"sector {name}: RSS summed (peak) {tot:.2f} GiB (predicted {PRED_RSS_GIB[0]:.0f}-"
            f"{PRED_RSS_GIB[1]:.0f}; beside H15b full sigma=800 {H15B_FULL_RSS_GIB} GiB at -n 16 "
            f"and the degree-2 copper full column '>445 s unreturned at -n 8'), per rank "
            f"{', '.join(f'{r:.2f}' for r in (per or []))}, dRSS {drss:+.2f} GiB; "
            f"S^{name}(P17,P17) d2 {s8[0, 0]:.6f} vs d1 {D1_S1717[sector]:.6f} (record, printed); "
            f"t={time.perf_counter() - t_window:.1f} s"
        )
        if comm.rank == 0:
            for ii, oi in enumerate(kept):
                print(
                    f"    S^{name} row P{oi}: "
                    + " ".join(f"{v.real:+.6f}{v.imag:+.6f}j" for v in s8[ii]),
                    flush=True,
                )
            with open(_path(name), "xb") as fh:
                np.savez(
                    fh,
                    s=s8,
                    kept=np.asarray(kept),
                    dofs_d2=dofs_d2,
                    dofs_d1=dofs_d1,
                    p_surf=p_surf,
                    t_factor=t_factor,
                    t_held_max=max(held),
                    t_held_mean=float(np.mean(held)),
                    t_sector=t_sec,
                    rss_sum_gib=tot,
                    n_dofs=n_dofs,
                    n_ranks=comm.size,
                )
            print(f"    saved {_path(name)}", flush=True)
        comm.Barrier()
        if not rec <= RECIPROCITY_BAND:
            fails.append(f"I1 sector {name} {rec:.3e}")
        if not smax <= MATRIX_PASSIVITY_CEILING:
            fails.append(f"I2 sector {name} {smax:.12f}")
        if not p_surf > 0.0:
            fails.append(f"P2 sector {name} P_surf {p_surf}")

    for k in (1, 2, 3):
        say(
            f"PROJECTION({k}) = {t_before_first:.1f} + {k} x {t_sector_max:.1f} + 30 = "
            f"{t_before_first + k * t_sector_max + 30.0:.0f} s (gate <= 540 s)"
        )
    say(f"window {time.perf_counter() - t_window:.1f} s; failures: {fails or 'none'}")
    assert not fails, fails


@complex_only
def test_th20_step2_assemble():
    if os.environ.get(ASSEMBLE_ENV, "").strip() != "1":
        pytest.skip(f"{ASSEMBLE_ENV} unset: TH-20 step 2a assembly is selected by environment only")
    from tests.validation.test_birdcage_f_human_rung import _ring_ports
    from tests.validation.test_port_birdcage_ring_matrix import MATRIX_PASSIVITY_CEILING
    from tests.validation.test_port_lumped_sheet_sweep import RECIPROCITY_BAND
    from tests.validation.test_th20_quarter_f_human_mesh import KEPT_RING

    comm = MPI.COMM_WORLD
    assert comm.size == 1, "test 2 is the -n 1 assembly window"
    fails = []

    def say(msg):
        print(f"[TH-20 step2a W3] {msg}", flush=True)

    ordinals = list(_ring_ports())
    kept = list(KEPT_RING)
    assert len(ordinals) == 32
    missing = [_path(NAME[s]) for s in SECTORS if not os.path.exists(_path(NAME[s]))]
    assert not missing, f"missing sector npz: {missing}"
    if os.path.exists(_path("R")):
        raise FileExistsError(f"stale {_path('R')} present")
    data = {}
    for s in SECTORS:
        with np.load(_path(NAME[s])) as z:
            data[s] = {k: z[k] for k in z.files}
        assert list(data[s]["kept"]) == kept, (NAME[s], data[s]["kept"])
    ss = {s: data[s]["s"] for s in SECTORS}

    # (P3)
    preimage, bij, worst_tau = _preimage_and_p3(kept, ordinals)
    say(
        f"(P3) image map bijection onto 32 ordinals: {bij}; worst |(R_g d_o)·d_pi(o) - tau_g| "
        f"{worst_tau:.3e} (<= 1e-12, ASSERTED)"
    )
    if not bij:
        fails.append("P3 bijection")
    if not worst_tau <= 1e-12:
        fails.append(f"P3 orientation {worst_tau:.3e}")

    # (P1) from the saved counts
    d = {s: int(data[s]["dofs_d2"]) for s in SECTORS}
    d_ee, d_em, d_me, d_mm = (d[s] for s in SECTORS)
    p1 = d_ee > d_em and d_ee > d_me and d_em > d_mm and d_me > d_mm
    say(
        f"(P1) degree-2 constrained DOFs EE {d_ee} > EM {d_em}, ME {d_me} > MM {d_mm}: {p1} "
        f"(ASSERTED); ratio to step 1's degree-1 counts "
        f"{', '.join(f'{NAME[s]} {d[s] / D1_DOFS[s]:.3f}' for s in SECTORS)} (predicted ≈2, printed)"
    )
    if not p1:
        fails.append(f"P1 {d_ee, d_em, d_me, d_mm}")

    # (I1)/(I2) re-read from the saved matrices, (P2) from the saved P_surf
    for s in SECTORS:
        rec, rec_max = _recip(ss[s])
        smax = _sigma_max(ss[s])
        ps = float(data[s]["p_surf"])
        say(
            f"sector {NAME[s]}: (I1) {rec:.3e} (<= {RECIPROCITY_BAND:g}, ASSERTED; entrywise "
            f"{rec_max:.3e}); (I2) sigma_max {smax:.12f} (<= {MATRIX_PASSIVITY_CEILING}, "
            f"ASSERTED); P_surf {ps:.9e} W (> 0, ASSERTED); price: unknowns {int(data[s]['n_dofs'])}, "
            f"factor {float(data[s]['t_factor']):.2f} s, held max {float(data[s]['t_held_max']):.2f} "
            f"/ mean {float(data[s]['t_held_mean']):.2f} s, sector {float(data[s]['t_sector']):.1f} s, "
            f"RSS summed {float(data[s]['rss_sum_gib']):.2f} GiB at -n {int(data[s]['n_ranks'])}"
        )
        say(
            f"sector {NAME[s]}: S^s(P17,P17) d1 {D1_S1717[s]:.6f} -> d2 {ss[s][0, 0]:.6f} "
            f"(|d2-d1|/|d1| {abs(ss[s][0, 0] - D1_S1717[s]) / abs(D1_S1717[s]):.3e}; record, printed)"
        )
        if not rec <= RECIPROCITY_BAND:
            fails.append(f"I1 {NAME[s]} {rec:.3e}")
        if not smax <= MATRIX_PASSIVITY_CEILING:
            fails.append(f"I2 {NAME[s]} {smax:.12f}")
        if not ps > 0.0:
            fails.append(f"P2 {NAME[s]} {ps}")

    # R (32×32), saved for step 2b
    r = _reconstruct(ss, kept, ordinals, preimage)
    with open(_path("R"), "xb") as fh:
        np.savez(
            fh,
            R=r,
            ordinals=np.asarray(ordinals),
            kept=np.asarray(kept),
            **{f"S_{NAME[s]}": ss[s] for s in SECTORS},
        )
    say(
        f"saved {_path('R')} (step 2b's comparand); R(P17,P17) d2 {r[0, 0]:.6f} beside step 1 "
        f"d1 {D1_R1717:.6f} (record, printed); R reciprocity {_recip(r)[0]:.3e}, sigma_max "
        f"{_sigma_max(r):.12f} (printed)"
    )

    # (C2) printed, never asserted (review ruling 2)
    c2 = min(
        float(
            abs(ss[("M", "E")][i, i] - ss[("E", "E")][i, i])
            / abs(sum(ss[s][i, i] for s in SECTORS))
        )
        for i in range(len(kept))
    )
    say(
        f"(C2) min_i |S^ME_ii - S^EE_ii| / |sum_s S^s_ii| = {c2:.4e} (floor {C2_FLOOR:g}, "
        f"PRINTED NOT ASSERTED per 2026-10-09 18:00 review ruling 2; predicted "
        f"[{C2_PREDICTED[0]}, {C2_PREDICTED[1]}]; step 1 degree-1 {D1_C2:.4e}); "
        f"{'above' if c2 >= C2_FLOOR else 'BELOW'} floor"
    )

    # step 2b projection (printed)
    fmax = max(float(data[s]["t_factor"]) for s in SECTORS)
    say(
        f"PROJECTION step 2b full degree-2 factor at -n 8: {fmax:.1f} s x 8 ... 16 = "
        f"{8 * fmax:.0f} ... {16 * fmax:.0f} s (projection, printed)"
    )
    say(f"failures: {fails or 'none'}")
    assert not fails, fails
