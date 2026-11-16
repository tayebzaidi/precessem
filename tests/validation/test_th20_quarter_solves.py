"""`TH-20` step 1 (§10 chain H13a, §9 item 233) — four quarter solves of the
F-human copper birdcage reconstruct the mirror-built full 32×32.

Selected by environment only (``FEM_EM_TH20_STEP1=1``), never by ``-k``.
Knobs: ``FEM_EM_TH20_STEP1_FULL_DRIVES=N`` drives the first N ring ports on the
full mesh (unset = 32; below 32 the matrix-wide gates print SKIPPED);
``FEM_EM_TH20_STEP1_CONTROL=1`` asserts the negative controls C1 / C1x / C2.
Step 2b-a (§9 item 412): ``FEM_EM_TH20_DEGREE=2`` runs the sectors and the
mirror-built full solve at degree 2 (unset = the degree-1 path below; any other
value raises); ``FEM_EM_TH20_STEP1_BUILD_ONLY=1`` builds the five spaces and
asserts (Deg) — unknowns = E at degree 1, 2E + 2F at degree 2 — with no form,
assembly or factorisation.

Fixture: degree 1, 64 MHz, `TH-14`'s Leontovich copper wall (σ = 5.8e7) on tag
401, 50 Ω exact ring sheets (κ off), the 1.6 mm hole quarter of step 0a and
its mirror-built full mesh of step 0b.  Sector (b_x, b_y): BC on x = 0 (501)
and y = 0 (502); E = PEC (pinned), M = PMC (natural).

Sector character (brief, §9 item 233): azimuthal ports flip under a
reflection, so χ_s(g) = ε_x^[σ_x∈g] ε_y^[σ_y∈g] with ε = +1 for E, −1 for M.
Reconstruction: S_{π_h(i), π_k(j)} = ¼ Σ_s χ_s(hk) S^s_ij, i, j ∈ K.

Anchors (asserted, bands imported): (P0) cell / facet records; (P1)
constrained-DOF ordering; (P2) P_surf > 0 on every operator; (P3) the
orientation premise and image bijection; (A) max |S_f − R|/|R| ≤
``LEG_D0_REPRODUCTION_BAND`` on the driven columns; (I1) reciprocity ≤
``RECIPROCITY_BAND``; (I2) σ_max ≤ ``MATRIX_PASSIVITY_CEILING`` (S_f in W2 and
every S^s); (I3) S_f's D2 invariance ≤ 1e-9 (W2).
Negative controls (rule (e); predicted and printed in W1, asserted under
``FEM_EM_TH20_STEP1_CONTROL=1``): C1 orientation trap, C1x PMC-for-PEC on
x = 0, C2 one mis-set plane in one sector.

Scope: one same-mesh, same-code identity; no two-mesh band, no C16, no
C_mode / tuning / B1+ / SAR, no degree 2 or memory claim (H13b), 64 MHz only.
"""

from __future__ import annotations

import gc
import os
import time

import numpy as np
import pytest
from mpi4py import MPI

from tests.complex_mode import complex_only

STEP1_ENV = "FEM_EM_TH20_STEP1"
DRIVES_ENV = "FEM_EM_TH20_STEP1_FULL_DRIVES"
CONTROL = os.environ.get("FEM_EM_TH20_STEP1_CONTROL") == "1"
COPPER_SIGMA = 5.8e7
OUTER_BOX_TAG = 499

# (P0) records: 0a quarter (`20261003T130342Z_TH-20-step0b-W2.log:8179`), 0b
# mirror-built cells and owned 401 facets (0b M4).
QUARTER_CELL_RECORD = 97_414
MIRROR_CELL_RECORD = 389_656
MIRROR_401_FACET_RECORD = 79_344
# (I3) D2 invariance bound (brief, §9 item 233).
D2_INVARIANCE_BAND = 1.0e-9
# C1 / C1x exact-2 signature tolerance and C2 floor (brief).
CONTROL_EXACT_TOL = 1.0e-9
C2_FLOOR = 1.0e-3
# Printed beside S_f(P17, P17), never compared (two meshes):
# WF-12 step 2's direct-mesh copper value.
WF12_S1717 = 0.383271 + 0.391623j
TH16_RATIO = 0.29
W2_PROJECTION_CEILING_S = 560.0
NORMAL_BUILD_S = 108.0

SECTORS = (("E", "E"), ("E", "M"), ("M", "E"), ("M", "M"))
GROUP = {"e": (0, 0), "sx": (1, 0), "sy": (0, 1), "sxsy": (1, 1)}
BITS_TO_G = {v: k for k, v in GROUP.items()}
TAU = {"sx": -1.0, "sy": -1.0, "sxsy": 1.0}


def _skip_unless_on():
    if os.environ.get(STEP1_ENV, "").strip() != "1":
        pytest.skip(f"{STEP1_ENV} unset: TH-20 step 1 is selected by environment only")


def _n_drives():
    raw = os.environ.get(DRIVES_ENV, "").strip()
    n = int(raw) if raw else 32
    if not 1 <= n <= 32:
        raise ValueError(f"{DRIVES_ENV} must be in 1..32, got {raw!r}")
    return n


DEGREE_ENV = "FEM_EM_TH20_DEGREE"
BUILD_ONLY_ENV = "FEM_EM_TH20_STEP1_BUILD_ONLY"
# Step 2a's degree-2 quarter count (`20261010T003627Z_TH-20-step2-W1-EE.log:8305`,
# `-n 8`): 684 904 = 2·134 230 + 2·208 222; printed beside (Deg-q).
Q_D2_UNKNOWNS_RECORD = 684_904
Q_EDGES_RECORD = 134_230
Q_FACES_RECORD = 208_222
SAVED_R_PATH = "/workspace/logs/th20/th20-step2-d2-R.npz"


def _degree_knob():
    """``FEM_EM_TH20_DEGREE`` (H20a's parser shape, `TH-20` step 2b-a): unset /
    blank ⇒ None (the degree-1 path), "2" ⇒ 2, anything else — "1" included,
    one degree-1 path only — ⇒ ValueError."""
    raw = os.environ.get(DEGREE_ENV, "").strip()
    if not raw:
        return None
    if raw == "2":
        return 2
    raise ValueError(f"{DEGREE_ENV} must be unset or '2', got {raw!r}")


def _effective_degree(knob, base):
    return base if knob is None else 2


def _build_only():
    raw = os.environ.get(BUILD_ONLY_ENV, "").strip()
    if not raw:
        return False
    if raw == "1":
        return True
    raise ValueError(f"{BUILD_ONLY_ENV} must be unset or '1', got {raw!r}")


def _deg_count(solver, msh, deg):
    """(Deg): global unknowns against E (degree 1) or 2E + 2F (degree 2)."""
    msh.topology.create_entities(1)
    msh.topology.create_entities(2)
    v_space = solver.function_space()
    n_dofs = int(v_space.dofmap.index_map.size_global * v_space.dofmap.index_map_bs)
    n_e = int(msh.topology.index_map(1).size_global)
    n_f = int(msh.topology.index_map(2).size_global)
    expect = n_e if deg == 1 else 2 * n_e + 2 * n_f
    ok = solver.degree == deg and n_dofs == expect
    return ok, n_e, n_f, n_dofs, expect


def _mirror_built_full(msh_q, cells_q, q, comm):
    """Step 0b's mirror build + full facet tags (lifted from the step-1 test,
    `TH-20` step 2b-a, rule (a)); connectivity hoisted before any facet query."""
    import dolfinx

    from tests.validation.test_th20_mirror_built_full_mesh import (
        _distribute,
        _gather_owned,
        _mirror,
    )
    from tests.validation.test_th20_quarter_f_human_mesh import _extents

    tdim = msh_q.topology.dim
    t0 = time.perf_counter()
    gathered = _gather_owned(msh_q, cells_q)
    x = c = t = None
    if comm.rank == 0:
        x, c, t = gathered
        for axis, g in ((1, "sy"), (0, "sx")):
            x, c, t = _mirror(x, c, t, axis, g)
    del gathered
    msh = _distribute(x, c, comm)
    fdim = msh.topology.dim - 1
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
        local_values.astype(np.int32, copy=False))
    del x, c, t, ents, vals
    msh.topology.create_connectivity(fdim, tdim)
    msh.topology.create_entity_permutations()
    comm.Barrier()
    t_mirror = time.perf_counter() - t0
    n_cells = int(msh.topology.index_map(tdim).size_global)
    r_e, z_e = _extents(q["kwargs"], float(q["diag"]["port_box_size_m"][1]))
    ftags, fcensus = _full_facet_tags(msh, r_e, z_e, comm)
    return {"msh": msh, "cells": cells, "ftags": ftags, "fcensus": fcensus,
            "n_cells": n_cells, "t_mirror": t_mirror}


def _pi(g, o):
    """Image of ring ordinal o under g (bottom 17 + j, top 33 + j)."""
    from tests.validation.test_th20_mirror_built_full_mesh import CLOSED

    if g == "e":
        return o
    base = 17 if o <= 32 else 33
    return base + CLOSED[g](o - base)


def _mul(h, k):
    a, b = GROUP[h], GROUP[k]
    return BITS_TO_G[(a[0] ^ b[0], a[1] ^ b[1])]


def _chi(sector, g, eps_e=1.0, flip_x=False):
    """χ_s(g).  ``eps_e`` = −1 is C1's axial convention; ``flip_x`` is C1x."""
    bx, by = sector
    if flip_x:
        bx = "M" if bx == "E" else "E"
    ex = eps_e if bx == "E" else -eps_e
    ey = eps_e if by == "E" else -eps_e
    gx, gy = GROUP[g]
    return (ex if gx else 1.0) * (ey if gy else 1.0)


def _reconstruct(ss, kept, ordinals, preimage, **chi_kw):
    """R (32×32, ordinals order) from the four 8×8 sector matrices."""
    n = len(ordinals)
    r = np.zeros((n, n), dtype=complex)
    pos = {o: k for k, o in enumerate(kept)}
    for ia, a in enumerate(ordinals):
        i, h = preimage[a]
        for ib, b in enumerate(ordinals):
            j, k = preimage[b]
            g = _mul(h, k)
            r[ia, ib] = 0.25 * sum(_chi(s, g, **chi_kw) * ss[s][pos[i], pos[j]]
                                   for s in SECTORS)
    return r


def _full_facet_tags(msh, r_e, z_e, comm):
    """0b's geometric rule on every local exterior facet (owned and ghost):
    box faces → 499, the rest → 401; owned counts reduced."""
    import dolfinx

    from tests.validation.test_th20_mirror_built_full_mesh import FACE_TOL
    from tests.validation.test_th20_quarter_f_human_mesh import PLANE_TOL

    tdim = msh.topology.dim
    fdim = tdim - 1
    n_owned = msh.topology.index_map(fdim).size_local
    ext = np.asarray(dolfinx.mesh.exterior_facet_indices(msh.topology), dtype=np.int32)
    vals = np.zeros(len(ext), dtype=np.int32)
    on_plane = 0
    if len(ext):
        g = dolfinx.cpp.mesh.entities_to_geometry(msh._cpp_object, fdim, ext, False)
        p = msh.geometry.x[g][:, :3, :]
        outer = (np.all(np.abs(np.abs(p[:, :, 0]) - r_e) < FACE_TOL, axis=1)
                 | np.all(np.abs(np.abs(p[:, :, 1]) - r_e) < FACE_TOL, axis=1)
                 | np.all(np.abs(np.abs(p[:, :, 2]) - z_e) < FACE_TOL, axis=1))
        vals = np.where(outer, OUTER_BOX_TAG, 401).astype(np.int32)
        own = ext < n_owned
        on_plane = int(np.count_nonzero(own & np.all(np.abs(p[:, :, 0]) < PLANE_TOL, axis=1))
                       + np.count_nonzero(own & np.all(np.abs(p[:, :, 1]) < PLANE_TOL, axis=1)))
    order = np.argsort(ext)
    tags = dolfinx.mesh.meshtags(msh, fdim, ext[order], vals[order])
    own = ext < n_owned
    census = {
        "outer": comm.allreduce(int(np.count_nonzero(own & (vals == OUTER_BOX_TAG))), op=MPI.SUM),
        "cav": comm.allreduce(int(np.count_nonzero(own & (vals == 401))), op=MPI.SUM),
        "on_plane": comm.allreduce(on_plane, op=MPI.SUM),
    }
    return tags, census


def _owned_tag_count(msh, facets, tag, comm):
    n_owned = msh.topology.index_map(msh.topology.dim - 1).size_local
    f = np.asarray(facets.find(int(tag)))
    return comm.allreduce(int(np.count_nonzero(f < n_owned)), op=MPI.SUM)


def _make_drive(msh, tags_f, wall_facets, specs, omega, z_s, comm):
    """Local copy of WF-12's ``_drive`` (`test_wf12_f_human_copper_matrix.py:
    170–220`): wall term first in the factorising call; S column by its
    formula (`:197–199`); P_surf on the factorising drive only."""
    from fem_em_solver.ports.lumped import (
        lumped_port_bilinear_term,
        lumped_port_linear_term,
        sheet_terminal_current,
    )
    from tests.validation.test_port_birdcage_four_port import TERMINATED_PORT_IMPEDANCE_OHM
    from tests.validation.test_th14_birdcage_copper import _leontovich_term, _surface_loss_w

    z_p = complex(TERMINATED_PORT_IMPEDANCE_OHM)
    wall = _leontovich_term(msh, wall_facets, z_s, omega)

    def drive(solver, driven_id, held):
        sheets = [sp.sheet(driven=(sp.port_id == driven_id)) for sp in specs]
        dsheet = next(s for s in sheets if s.port_id == driven_id)
        linear = [lambda test, _s=dsheet: lumped_port_linear_term(
            msh, tags_f, _s, test, omega_rad_per_s=omega)]
        comm.Barrier()
        t0 = time.perf_counter()
        if held:
            fields = solver.solve_with_held_factorization(extra_linear_terms=linear)
        else:
            bil = [wall] + [
                (lambda trial, test, _s=s: lumped_port_bilinear_term(
                    msh, tags_f, _s, trial, test, omega_rad_per_s=omega))
                for s in sheets
            ]
            fields = solver.solve(current_density=None, project_source=False,
                                  extra_bilinear_terms=bil, extra_linear_terms=linear)
        comm.Barrier()
        t_solve = time.perf_counter() - t0
        t1 = time.perf_counter()
        e = fields.e_complex
        currents = {s.port_id: sheet_terminal_current(msh, tags_f, s, e, comm) for s in sheets}
        v_src = complex(dsheet.source_voltage_v)
        s_col = {}
        for s in sheets:
            v = complex(s.source_voltage_v) - currents[s.port_id] * z_p
            s_col[s.port_id] = (v - z_p * currents[s.port_id]) / v_src
        p_surf = None if held else _surface_loss_w(msh, wall_facets, z_s, e, comm)
        del fields, e
        t_post = time.perf_counter() - t1
        return {"s": s_col, "held": held, "t_solve": t_solve, "t_post": t_post,
                "p_surf": p_surf}

    return drive


def _specs(ordinals, width_of, chord):
    from fem_em_solver.ports.lumped import LumpedSheetPortSpec
    from tests.mesh.test_birdcage_port_sheets import SHEET_IFACE
    from tests.mesh.test_birdcage_ring_gaps_scaleup import _ring_gap_frame
    from tests.validation.test_birdcage_f_human_rung import LEG_COUNT
    from tests.validation.test_port_birdcage_four_port import TERMINATED_PORT_IMPEDANCE_OHM

    out = []
    for o in ordinals:
        phi_hat, _c = _ring_gap_frame(o, LEG_COUNT)
        out.append(LumpedSheetPortSpec(
            port_id=f"P{o}", facet_tag=int(SHEET_IFACE + o),
            port_impedance_ohm=TERMINATED_PORT_IMPEDANCE_OHM, gap_height_m=chord,
            sheet_width_m=width_of[o], drive_direction=tuple(float(c) for c in phi_hat),
            drive_voltage_v=1.0 + 0.0j, interior=True,
        ))
    return out


def _sigma_max(s):
    return float(np.linalg.norm(s, 2))


def _recip(s):
    return (float(np.linalg.norm(s - s.T, "fro") / np.linalg.norm(s, "fro")),
            float(np.max(np.abs(s - s.T))))


@complex_only
def test_th20_step1_quarter_solves_reconstruct_full():
    _skip_unless_on()
    from fem_em_solver.core import TimeHarmonicSolver
    from fem_em_solver.core.cavity import surface_resistance_ohm
    from fem_em_solver.io.mesh import _interface_facet_tags
    from tests.mesh.test_birdcage_port_sheets import PORT_LOWER, PORT_UPPER, SHEET_IFACE
    from tests.mesh.test_birdcage_port_terminals import _interface_area_or_zero
    from tests.mesh.test_birdcage_ring_gaps_scaleup import _ring_gap_frame
    from tests.validation.test_birdcage_f_human_rung import LEG_COUNT, _ring_ports
    from tests.validation.test_port_birdcage_four_port import LEG_D0_REPRODUCTION_BAND
    from tests.validation.test_port_birdcage_ring_matrix import MATRIX_PASSIVITY_CEILING
    from tests.validation.test_port_lumped_sheet_sweep import RECIPROCITY_BAND
    from tests.validation.test_th14_birdcage_copper import _problem
    from tests.validation.test_th16_quarter_birdcage import (
        RECIPROCITY_Q_BAND,
        _constrained_dofs,
    )
    from tests.validation.test_th20_mirror_built_full_mesh import (
        CUT_FACES_X0,
        CUT_FACES_Y0,
        OPS,
    )
    from tests.validation.test_th20_quarter_degree2 import _constrained_dofs_degree
    from tests.validation.test_th20_quarter_f_human_mesh import (
        KEPT_RING,
        _meshes,
    )
    from tests.validation.test_wf7_f_human_b1_quadrature import DEGREE, FREQUENCY_HZ
    from tests.validation.test_wf12_f_human_hole import _rss_gib

    comm = MPI.COMM_WORLD
    t_window = time.perf_counter()
    n_drives = _n_drives()
    knob = _degree_knob()
    deg = _effective_degree(knob, DEGREE)
    build_only = _build_only()
    fails = []
    deg_q_counts = []

    def say(msg):
        if comm.rank == 0:
            print(f"[TH-20 step1] {msg}", flush=True)

    def rss_sum():
        return _rss_gib(comm)[2]

    say(f"start: -n {comm.size}, degree {deg} ({DEGREE_ENV} "
        f"{'unset' if knob is None else knob}), build-only {build_only}, f "
        f"{FREQUENCY_HZ:.3e} Hz, full drives {n_drives}, CONTROL "
        f"{'asserted' if CONTROL else 'printed'}")

    def deg_check(solver, msh, where):
        ok, n_e, n_f, n_dofs, expect = _deg_count(solver, msh, deg)
        say(f"(Deg) {where}: solver.degree {solver.degree}; global unknowns {n_dofs} vs "
            f"{'E' if deg == 1 else '2E+2F'} = {expect} (E {n_e}, F {n_f}; ASSERTED exact); "
            f"n/E {n_dofs / n_e:.4f}")
        if not ok:
            fails.append(f"Deg {where} degree {solver.degree} n {n_dofs} expect {expect}")
        return n_e, n_f, n_dofs
    omega = 2.0 * np.pi * FREQUENCY_HZ
    z_s = (1.0 + 1.0j) * surface_resistance_ohm(omega, COPPER_SIGMA)
    ordinals = list(_ring_ports())
    kept = list(KEPT_RING)
    assert len(ordinals) == 32

    # ---- (P3) orientation premise + image bijection (pure arithmetic) -----
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
    say(f"(P3) image map bijection onto 32 ordinals: {bij}; worst |(R_g d_o)·d_pi(o) - tau_g| "
        f"{worst_tau:.3e} (<= 1e-12, ASSERTED)")
    if not bij:
        fails.append("P3 image map not a bijection")
    if not worst_tau <= 1e-12:
        fails.append(f"P3 orientation {worst_tau:.3e}")

    # ---- step 1: meshes ---------------------------------------------------
    rss0 = rss_sum()
    m = _meshes()
    full_d, q = m["full"], m["quarter"]
    msh_q, cells_q, facets_q = q["msh"], q["cells"], q["facets"]
    t_full_build = float(full_d["t_build"])
    speed = t_full_build / NORMAL_BUILD_S
    n401_q = _owned_tag_count(msh_q, facets_q, 401, comm)
    n501 = _owned_tag_count(msh_q, facets_q, 501, comm)
    n502 = _owned_tag_count(msh_q, facets_q, 502, comm)
    say(f"(P0) quarter cells {q['n_cells']} (record {QUARTER_CELL_RECORD}); owned 501/502 "
        f"{n501}/{n502} (records {CUT_FACES_X0}/{CUT_FACES_Y0}); quarter owned 401 {n401_q}; "
        f"direct full build {t_full_build:.1f} s -> box speed s = {speed:.3f} "
        f"(normal 109.9 / 107.2 s); t={time.perf_counter() - t_window:.1f} s")
    if q["n_cells"] != QUARTER_CELL_RECORD:
        fails.append(f"P0 quarter cells {q['n_cells']}")
    if (n501, n502) != (CUT_FACES_X0, CUT_FACES_Y0):
        fails.append(f"P0 cut facets {(n501, n502)}")

    pec_of = {("E", "E"): (1, 501, 502), ("E", "M"): (1, 501), ("M", "E"): (1, 502),
              ("M", "M"): (1,)}

    def check_quarter_deg(counts):
        """(Deg-q): the four sector counts equal; at degree 2 the backed negative
        control n(2) != E and n(2)·E == E·(2E + 2F) in exact integers."""
        if len(set(counts)) != 1:
            fails.append(f"Deg-q sector counts differ {counts}")
        n_e, n_f, n_d = counts[0]
        say(f"(Deg-q) four sector counts equal: {len(set(counts)) == 1} (ASSERTED); "
            f"unknowns {n_d}, E {n_e}, F {n_f}"
            + (f" beside step 2a's record {Q_D2_UNKNOWNS_RECORD} = 2*{Q_EDGES_RECORD} + "
               f"2*{Q_FACES_RECORD} (-n 8; printed)" if deg == 2 else
               f" beside quarter E record {Q_EDGES_RECORD} (printed)"))
        if deg == 2:
            ctl = n_d != n_e and n_d * n_e == n_e * (2 * n_e + 2 * n_f)
            say(f"(Deg-q) negative control: degree-2 count {n_d} != degree-1 count E {n_e} "
                f"and n*E == E*(2E+2F): {ctl} (ASSERTED; ratio {n_d / n_e:.3f}, record 5.102)")
            if not ctl:
                fails.append(f"Deg-q control {n_d} vs E {n_e}")

    def check_mirror_p0(full):
        fc = full["fcensus"]
        say(f"(P0) mirror-built cells {full['n_cells']} (record {MIRROR_CELL_RECORD}); build "
            f"{full['t_mirror']:.2f} s; owned 401 {fc['cav']} (record {MIRROR_401_FACET_RECORD}; "
            f"4x quarter {4 * n401_q}); owned 499 {fc['outer']}; exterior on planes "
            f"{fc['on_plane']} (record 0)")
        if full["n_cells"] != MIRROR_CELL_RECORD:
            fails.append(f"P0 mirror cells {full['n_cells']}")
        if fc["cav"] != MIRROR_401_FACET_RECORD:
            fails.append(f"P0 mirror 401 {fc['cav']}")
        if fc["on_plane"] != 0:
            fails.append(f"P0 on-plane {fc['on_plane']}")

    tdim = msh_q.topology.dim
    msh_q.topology.create_connectivity(tdim - 1, tdim)
    msh_q.topology.create_entity_permutations()

    # ---- build-only branch (`TH-20` step 2b-a, W1): no drive, no sheet tags, no
    # form, no assembly, no factorisation -----------------------------------
    if build_only:
        for sector in SECTORS:
            problem = _problem(msh_q, cells_q, FREQUENCY_HZ, facets_q, pec_of[sector])
            solver = TimeHarmonicSolver(problem, degree=deg)
            deg_q_counts.append(deg_check(solver, msh_q, f"sector {sector}"))
            del solver, problem
            gc.collect()
        check_quarter_deg(deg_q_counts)
        full = _mirror_built_full(msh_q, cells_q, q, comm)
        check_mirror_p0(full)
        problem_f = _problem(full["msh"], full["cells"], FREQUENCY_HZ, full["ftags"],
                             (OUTER_BOX_TAG,))
        solver_f = TimeHarmonicSolver(problem_f, degree=deg)
        n_e, _n_f, n_d = deg_check(solver_f, full["msh"], "mirror-built full")
        if deg == 2:
            say(f"(Deg-f) mirror-built degree-2 unknowns {n_d} (predicted 2 711 047, range "
                f"[2.65e6, 2.76e6], printed); n/E {n_d / n_e:.4f} (predicted 5.10-5.14, printed)")
        del solver_f, problem_f
        gc.collect()
        _mine, per, tot = _rss_gib(comm)
        say(f"build-only RSS per rank {', '.join(f'{r:.2f}' for r in (per or []))} GiB, summed "
            f"{tot:.2f} GiB (start {rss0:.2f}); window {time.perf_counter() - t_window:.1f} s; "
            f"failures: {fails or 'none'}")
        assert not fails, fails
        return

    # ---- step 2: quarter sectors -----------------------------------------
    tags_q = _interface_facet_tags(
        msh_q, cells_q, {SHEET_IFACE + o: (PORT_LOWER + o, PORT_UPPER + o) for o in kept})
    chord = float(q["diag"]["ring_port_layout"]["ring_port_gap_chord_m"])
    width_of = {}
    for o in kept:
        width_of[o] = _interface_area_or_zero(msh_q, tags_q, SHEET_IFACE + o, comm) / chord
    for a in ordinals:  # each image sheet takes its preimage's spec numbers
        width_of[a] = width_of[preimage[a][0]]
    say(f"quarter sheet widths (m) {', '.join(f'P{o} {width_of[o]:.9e}' for o in kept)}; "
        f"chord {chord:.9e} m")
    specs_q = _specs(kept, width_of, chord)
    drive_q = _make_drive(msh_q, tags_q, facets_q, specs_q, omega, z_s, comm)
    ss, dofs, p_surf_s, t_sector, drss_sector = {}, {}, {}, {}, {}
    t_quarter_solve = 0.0
    for sector in SECTORS:
        r_before = rss_sum()
        t0 = time.perf_counter()
        pec = pec_of[sector]
        if deg == 2:
            dofs[sector] = _constrained_dofs_degree(msh_q, facets_q, pec, 2)
        else:
            dofs[sector] = _constrained_dofs(msh_q, facets_q, pec)
        problem = _problem(msh_q, cells_q, FREQUENCY_HZ, facets_q, pec)
        solver = TimeHarmonicSolver(problem, degree=deg)
        deg_q_counts.append(deg_check(solver, msh_q, f"sector {sector}"))
        cols = {}
        tsum = 0.0
        for k, o in enumerate(kept):
            d = drive_q(solver, f"P{o}", held=(k > 0))
            cols[o] = d
            tsum += d["t_solve"] + d["t_post"]
            if k == 0:
                p_surf_s[sector] = d["p_surf"]
        del solver, problem
        gc.collect()
        s8 = np.array([[cols[oj]["s"][f"P{oi}"] for oj in kept] for oi in kept], dtype=complex)
        ss[sector] = s8
        t_sector[sector] = time.perf_counter() - t0
        t_quarter_solve += tsum
        drss_sector[sector] = rss_sum() - r_before
        f0 = cols[kept[0]]
        held = [cols[o]["t_solve"] + cols[o]["t_post"] for o in kept[1:]]
        rec, rec_max = _recip(s8)
        smax = _sigma_max(s8)
        say(f"sector {sector}: pec {pec}, constrained DOFs {dofs[sector]}; factor "
            f"{f0['t_solve']:.2f} + post {f0['t_post']:.2f} s, held max {max(held):.2f} s mean "
            f"{np.mean(held):.2f} s; sector total {t_sector[sector]:.1f} s; P_surf "
            f"{p_surf_s[sector]:.9e} W; (I1) {rec:.3e} (band {RECIPROCITY_BAND:g}, ASSERTED; "
            f"entrywise max {rec_max:.3e} beside {RECIPROCITY_Q_BAND:g}, printed); (I2) "
            f"sigma_max {smax:.12f} (ceiling {MATRIX_PASSIVITY_CEILING}, ASSERTED); dRSS "
            f"{drss_sector[sector]:+.2f} GiB; t={time.perf_counter() - t_window:.1f} s")
        if comm.rank == 0:
            for ii, oi in enumerate(kept):
                print(f"    S^{sector[0]}{sector[1]} row P{oi}: "
                      + " ".join(f"{v.real:+.6f}{v.imag:+.6f}j" for v in s8[ii]), flush=True)
        if not rec <= RECIPROCITY_BAND:
            fails.append(f"I1 sector {sector} {rec:.3e}")
        if not smax <= MATRIX_PASSIVITY_CEILING:
            fails.append(f"I2 sector {sector} {smax:.12f}")
        if not p_surf_s[sector] > 0.0:
            fails.append(f"P2 sector {sector} P_surf {p_surf_s[sector]}")
    d_ee, d_em, d_me, d_mm = (dofs[s] for s in SECTORS)
    p1 = d_ee > d_em and d_ee > d_me and d_em > d_mm and d_me > d_mm
    say(f"(P1) constrained DOFs EE {d_ee} > EM {d_em}, ME {d_me} > MM {d_mm}: {p1} (ASSERTED)")
    if not p1:
        fails.append(f"P1 dofs {d_ee, d_em, d_me, d_mm}")
    check_quarter_deg(deg_q_counts)

    # ---- steps 3-4: mirror-built full mesh + full facet tags (lifted) ------
    full = _mirror_built_full(msh_q, cells_q, q, comm)
    msh, cells, ftags = full["msh"], full["cells"], full["ftags"]
    check_mirror_p0(full)

    # ---- step 5: full sheets -----------------------------------------------
    tags_full = _interface_facet_tags(
        msh, cells, {SHEET_IFACE + o: (PORT_LOWER + o, PORT_UPPER + o) for o in ordinals})
    specs_f = _specs(ordinals, width_of, chord)

    # ---- step 6: full sweep ------------------------------------------------
    drive_f = _make_drive(msh, tags_full, ftags, specs_f, omega, z_s, comm)
    problem_f = _problem(msh, cells, FREQUENCY_HZ, ftags, (OUTER_BOX_TAG,))
    solver_f = TimeHarmonicSolver(problem_f, degree=deg)
    deg_check(solver_f, msh, "mirror-built full")
    t_before_sweep = time.perf_counter() - t_window
    r_before = rss_sum()
    driven = ordinals[:n_drives]
    cols_f = {}
    for k, o in enumerate(driven):
        d = drive_f(solver_f, f"P{o}", held=(k > 0))
        cols_f[o] = d
        say(f"full drive {k + 1:2d}/{n_drives} P{o}: {'held' if k else 'factor'} "
            f"{d['t_solve']:.2f} s + post {d['t_post']:.2f} s; S_driven "
            f"{d['s'][f'P{o}']:.6f}"
            + (f"; P_surf {d['p_surf']:.9e} W" if k == 0 else "")
            + f"; t={time.perf_counter() - t_window:.1f} s")
    drss_full = rss_sum() - r_before
    f0 = cols_f[driven[0]]
    t_factor = f0["t_solve"] + f0["t_post"]
    held = [cols_f[o]["t_solve"] + cols_f[o]["t_post"] for o in driven[1:]]
    max_held = max(held) if held else float("nan")
    mean_held = float(np.mean(held)) if held else float("nan")
    projection = t_before_sweep + t_factor + 31 * max_held + 30.0
    full_solve_32 = t_factor + 31 * mean_held
    full_solve_meas = t_factor + sum(held)
    p_surf_f = f0["p_surf"]
    rss_mine, rss_per, rss_total = _rss_gib(comm)
    say(f"(P2) P_surf sectors {', '.join(f'{p_surf_s[s]:.6e}' for s in SECTORS)} W, full "
        f"{p_surf_f:.9e} W (all > 0, ASSERTED)")
    if not p_surf_f > 0.0:
        fails.append(f"P2 full P_surf {p_surf_f}")
    say(f"PRICE: elapsed before full sweep {t_before_sweep:.1f} s; full factor {t_factor:.2f} s; "
        f"held max {max_held:.2f} s mean {mean_held:.2f} s; box speed s = {speed:.3f}; "
        f"PROJECTED W2 = {t_before_sweep:.1f} + {t_factor:.1f} + 31 x {max_held:.2f} + 30 = "
        f"{projection:.0f} s (gate <= {W2_PROJECTION_CEILING_S:.0f} s; predicted W2 330-480 s)")
    say(f"RSS per rank {', '.join(f'{r:.2f}' for r in (rss_per or []))} GiB, summed "
        f"{rss_total:.2f} GiB (start {rss0:.2f}); dRSS sectors "
        f"{', '.join(f'{drss_sector[s]:+.2f}' for s in SECTORS)} GiB, full sweep "
        f"{drss_full:+.2f} GiB")
    say(f"quarter-to-full solve-phase ratio: four sectors {t_quarter_solve:.1f} s / full 32-drive "
        f"{full_solve_32:.1f} s (factor + 31 x mean held) = "
        f"{t_quarter_solve / full_solve_32:.3f} (TH-16 read ≈{TH16_RATIO}; printed); measured "
        f"full ({n_drives} drives) {full_solve_meas:.1f} s")

    # ---- (A) the anchor ----------------------------------------------------
    r = _reconstruct(ss, kept, ordinals, preimage)
    sf = np.full((32, 32), np.nan + 0j, dtype=complex)
    for o in driven:
        jb = ordinals.index(o)
        for ia, a in enumerate(ordinals):
            sf[ia, jb] = cols_f[o]["s"][f"P{a}"]
    dcols = [ordinals.index(o) for o in driven]
    sfd, rd = sf[:, dcols], r[:, dcols]
    rel = np.abs(sfd - rd) / np.abs(rd)
    a_max = float(rel.max())
    frob = float(np.linalg.norm(sfd - rd) / np.linalg.norm(sfd))
    min_r = float(np.abs(r).min())
    # by-g pattern
    gpat = {}
    for ia, a in enumerate(ordinals):
        for kb, jb in enumerate(dcols):
            g = _mul(preimage[a][1], preimage[ordinals[jb]][1])
            gpat[g] = max(gpat.get(g, 0.0), float(rel[ia, kb]))
    say(f"(A) max |S_f - R|/|R| over {rel.size} entries {a_max:.3e} (band "
        f"{LEG_D0_REPRODUCTION_BAND:g}, ASSERTED; predicted 1e-13 - 1e-10); by hk "
        f"{ {g: f'{v:.3e}' for g, v in gpat.items()} }; Frobenius ||S_f - R||/||S_f|| "
        f"{frob:.3e}; min |R_ab| (all 1024) {min_r:.4e} (predicted >= 1e-3)")
    if not a_max <= LEG_D0_REPRODUCTION_BAND:
        fails.append(f"A {a_max:.3e}")
    say(f"S_f(P17,P17) {sf[0, 0]:.6f} beside WF-12 direct-mesh copper {WF12_S1717:.6f} "
        f"(two meshes; printed); R(P17,P17) {r[0, 0]:.6f}")
    if deg == 2 and comm.rank == 0:  # rank-0 numpy only, after every collective
        if os.path.exists(SAVED_R_PATH):
            with np.load(SAVED_R_PATH, allow_pickle=False) as z:
                ords_saved = [int(v) for v in np.asarray(z["ordinals"]).ravel()]
                r_saved = np.asarray(z["R"]) if "R" in z.files else None
            if ords_saved != [int(o) for o in ordinals] or r_saved is None:
                say(f"saved R {SAVED_R_PATH}: ordinals/keys do not match the window (printed)")
            else:
                rs = r_saved[:, dcols]
                dev = float(np.max(np.abs(rd - rs) / np.abs(rs)))
                say(f"saved R vs window R over driven columns: max |dR|/|R_saved| {dev:.3e} "
                    f"(predicted <= 1e-10; printed, never asserted)")
        else:
            say(f"saved R {SAVED_R_PATH}: absent")

    # ---- matrix-wide gates (W2) ---------------------------------------------
    if n_drives == 32:
        rec, rec_max = _recip(sf)
        smax = _sigma_max(sf)
        i3 = 0.0
        idx = {o: k for k, o in enumerate(ordinals)}
        for g in ("sx", "sy", "sxsy"):
            perm = [idx[_pi(g, a)] for a in ordinals]
            i3 = max(i3, float(np.max(np.abs(sf[np.ix_(perm, perm)] - sf) / np.abs(sf))))
        say(f"(I1) S_f reciprocity {rec:.3e} (band {RECIPROCITY_BAND:g}, ASSERTED; entrywise max "
            f"{rec_max:.3e} beside {RECIPROCITY_Q_BAND:g}, printed); (I2) S_f sigma_max "
            f"{smax:.12f} (ceiling {MATRIX_PASSIVITY_CEILING}, ASSERTED); max sigma_max(S^s) "
            f"{max(_sigma_max(ss[s]) for s in SECTORS):.12f}; (I3) D2 invariance {i3:.3e} "
            f"(<= {D2_INVARIANCE_BAND:g}, ASSERTED)")
        if not rec <= RECIPROCITY_BAND:
            fails.append(f"I1 S_f {rec:.3e}")
        if not smax <= MATRIX_PASSIVITY_CEILING:
            fails.append(f"I2 S_f {smax:.12f}")
        if not i3 <= D2_INVARIANCE_BAND:
            fails.append(f"I3 {i3:.3e}")
    else:
        say(f"(I1)/(I2)/(I3) on S_f: SKIPPED ({n_drives} < 32 full drives)")

    # ---- negative controls --------------------------------------------------
    def _exact2(r_ctl, label):
        miss = np.abs(r_ctl[:, dcols] - sfd) / np.abs(rd)
        two = np.abs(miss / 2.0 - 1.0) <= CONTROL_EXACT_TOL
        same = miss <= CONTROL_EXACT_TOL
        n_two, n_same = int(two.sum()), int(same.sum())
        ok = n_two == miss.size // 2 and n_same == miss.size - n_two
        worst2 = float(np.max(np.abs(miss[two] / 2.0 - 1.0))) if n_two else float("nan")
        say(f"({label}) entries at miss = 2: {n_two}, unchanged: {n_same}, of {miss.size} "
            f"(predicted {miss.size // 2} / {miss.size // 2}); worst |miss/2 - 1| {worst2:.3e}, "
            f"worst unchanged {float(miss[same].max()) if n_same else float('nan'):.3e} "
            f"({'ASSERTED' if CONTROL else 'printed'}): {'PASS' if ok else 'FAIL'}")
        return ok

    ok_c1 = _exact2(_reconstruct(ss, kept, ordinals, preimage, eps_e=-1.0), "C1")
    ok_c1x = _exact2(_reconstruct(ss, kept, ordinals, preimage, flip_x=True), "C1x")
    ss_c2 = dict(ss)
    ss_c2[("E", "E")] = ss[("M", "E")]
    r_c2 = _reconstruct(ss_c2, kept, ordinals, preimage)
    diag_idx = [ordinals.index(o) for o in driven]
    c2 = min(float(abs(r_c2[a, a] - r[a, a]) / abs(sf[a, a])) for a in diag_idx)
    say(f"(C2) min over {len(diag_idx)} self entries |dR_aa|/|S_f,aa| {c2:.4e} (floor "
        f"{C2_FLOOR:g}, {'ASSERTED' if CONTROL else 'printed'}; predicted [1e-2, 0.91])")
    if CONTROL:
        if not ok_c1:
            fails.append("C1 not exact-2 on half")
        if not ok_c1x:
            fails.append("C1x not exact-2 on half")
        if not c2 >= C2_FLOOR:
            fails.append(f"C2 {c2:.3e}")

    del solver_f, problem_f
    gc.collect()
    say(f"window {time.perf_counter() - t_window:.1f} s; failures: {fails or 'none'}")
    assert not fails, fails


def test_th20_step2ba_knob_coding(monkeypatch):
    """(K0) `TH-20` step 2b-a, ungated, numpy only: the degree and build-only
    parsers (rejected values are the parser's negative control — "1" above all,
    the second degree-1 path H20a forbade) and the records' own arithmetic."""
    from tests.validation.test_wf7_f_human_b1_quadrature import DEGREE

    for raw in (None, "", "  "):
        if raw is None:
            monkeypatch.delenv(DEGREE_ENV, raising=False)
        else:
            monkeypatch.setenv(DEGREE_ENV, raw)
        assert _degree_knob() is None
    for raw in ("2", " 2 "):
        monkeypatch.setenv(DEGREE_ENV, raw)
        assert _degree_knob() == 2
    for raw in ("1", "3", "x", "2.0"):
        monkeypatch.setenv(DEGREE_ENV, raw)
        with pytest.raises(ValueError):
            _degree_knob()
    assert DEGREE == 1
    assert _effective_degree(None, DEGREE) == 1
    assert _effective_degree(2, DEGREE) == 2

    monkeypatch.delenv(BUILD_ONLY_ENV, raising=False)
    assert _build_only() is False
    monkeypatch.setenv(BUILD_ONLY_ENV, "")
    assert _build_only() is False
    monkeypatch.setenv(BUILD_ONLY_ENV, "1")
    assert _build_only() is True
    for raw in ("0", "yes", "2"):
        monkeypatch.setenv(BUILD_ONLY_ENV, raw)
        with pytest.raises(ValueError):
            _build_only()

    assert Q_D2_UNKNOWNS_RECORD == 2 * Q_EDGES_RECORD + 2 * Q_FACES_RECORD == 684_904
    n_e, n_d = Q_EDGES_RECORD, Q_D2_UNKNOWNS_RECORD
    assert n_d != n_e and n_d * n_e == n_e * (2 * n_e + 2 * Q_FACES_RECORD)
