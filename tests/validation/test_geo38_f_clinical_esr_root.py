"""GEO-38 C4a (§9 item 561, §10 chain C4a): the F-clinical-1.5T in-model
Larmor root module — H26a's pattern (`test_wf10_f_human_copper_esr_root.py`,
§9 items 334 / 354 / 459) transposed from F-human's 16 legs / 32 sheets to
F-clinical's 32 legs / 64 sheets.

Brief: ``docs/planning/chunks/GEO-38.md`` "§9 item 561" and its review ruling
2026-10-20 03:00 (which binds over the brief): this module carries H26a's
root-finding route — (Z50), (P50), NEG, (P_RC) — not the gates entry 46's
readout names (circuit <-> field closure, (P_w), the feed-bisector mirror);
H17's (R) gate is dropped (C3's npz carries no per-drive power).

Selected by ``FEM_EM_GEO38_C4`` only (never by ``-k``):

* unset — the gated test skips; the numpy coding tests (K0)–(K3), (Y0-64),
  (N1), (N2) and the planted (Z) control run;
* ``build`` — build-only proof: C3a's fixture asserts (A0) / (iv) / outer
  census / (S) / (C) / (W) / (Az), (U1), the two-port feed (F), the ESR
  series-RC specs (Z), then degree 2 in the same build (Deg-2); returns before
  the Leontovich term, any assembly or factorisation; the npz is never opened;
* ``1`` — the node run (``xl-pending.md`` entry 46 only), nodes from
  ``FEM_EM_GEO38_C4_NODES_PF`` (required);
* anything else — ``ValueError``.

``FEM_EM_GEO38_C4_DEGREE``: unset ⇒ degree 1; ``2`` ⇒ degree 2; else
``ValueError`` (H20a's ``_degree_knob`` semantics under its own name, so a
GEO-38 window never reads a WF-10 variable).

F-human's 16-leg helpers (``_ring_quadrature_slots``, ``QUADRATURE_SLOT_OFFSET``,
``_k1_shares``, ``_apply_degree``) are never used on a production path here;
(N1) / (N2) import two of them inside their bodies only, as the asserted
negative controls that show why.

Scope: the module, its coding tests and a build-only proof.  Not claimed: any
Y_in, root, C_in, C_mode or tuning statement; any S entry or (Z50) / (P_RC)
reading on the fixture; degree 2 beyond the unknown count.
"""

from __future__ import annotations

import json
import os
import time

import numpy as np
import pytest
from mpi4py import MPI

import tests.validation.test_geo38_f_clinical_column as col
import tests.validation.test_geo38_f_clinical_matrix as c3
import tests.validation.test_wf10_f_human_copper_esr_root as h17
from tests.complex_mode import complex_only

C4_ENV = "FEM_EM_GEO38_C4"
C4_NODES_ENV = "FEM_EM_GEO38_C4_NODES_PF"
C4_DEGREE_ENV = "FEM_EM_GEO38_C4_DEGREE"
DRIVEN = c3.DRIVEN  # "P33"
#: Capacitor ESR (§9 binding ruling "capacitor ESR"), imported, never retyped.
ESR_OHM = h17.ESR_OHM
LEG_COUNT = c3.LEG_COUNT_32
AZIMUTH_STEP_DEG = c3.AZIMUTH_STEP_32_DEG
#: 90° on the 11.25° grid — a derived count (H17's QUADRATURE_SLOT_OFFSET is 4 on 22.5°).
QUADRATURE_SLOT_OFFSET_32 = LEG_COUNT // 4
JSON_DIR = "/workspace/logs/geo38"
#: (Z)'s wiring value only — not a node, not a prediction.
SPEC_WIRING_C_PF = 100.0
#: F-clinical-1.5T Larmor frequency (brief: 64 MHz); asserted against ctx["omega"].
FREQUENCY_HZ_C4 = 64.0e6

# Predictions written first (brief, printed only, never asserted).
PREDICTED_T_BUILD_S = (320.0, 560.0)
PREDICTED_T_D2_S = (10.0, 60.0)
PREDICTED_RSS_BUILD_GIB = 2.0
PREDICTED_RSS_D2_GIB = (4.0, 6.0)
PREDICTED_FEED2 = "P41"
PREDICTED_AZ_P33_DEG = 5.625
PREDICTED_BISECTOR_DEG = 50.625
PREDICTED_Z50 = 1.0e-10
PREDICTED_P_RC = 1.0e-10


def _say(comm, msg):
    if comm.rank == 0:
        print(f"[GEO-38 C4] {msg}", flush=True)


# ---------------------------------------------------------------------------
# Knobs.
# ---------------------------------------------------------------------------
def _mode():
    raw = os.environ.get(C4_ENV, "").strip()
    if raw == "":
        return None
    if raw in ("build", "1"):
        return raw
    raise ValueError(f"{C4_ENV} must be unset, 'build' or '1', got {raw!r}")


def _degree_knob_c4():
    """Unset / empty ⇒ ``None`` (degree 1); ``"2"`` ⇒ 2; anything else ⇒ ValueError."""
    raw = os.environ.get(C4_DEGREE_ENV, "").strip()
    if raw == "":
        return None
    if raw == "2":
        return 2
    raise ValueError(f"{C4_DEGREE_ENV}={raw!r}: only unset/empty (degree 1) or '2' are accepted")


def _nodes_pf():
    """``C4_NODES_ENV`` (comma-separated pF), required: every node > 0, sorted
    ascending, distinct and round-tripping at 4 decimals (the node JSON name)."""
    raw = os.environ.get(C4_NODES_ENV, "").strip()
    if raw == "":
        raise ValueError(f"{C4_NODES_ENV} is required under {C4_ENV}=1 (C10 sets the grid)")
    nodes = [float(x) for x in raw.split(",") if x.strip()]
    if not nodes:
        raise ValueError(f"{C4_NODES_ENV}={raw!r}: no nodes")
    if any(not np.isfinite(c) or c <= 0.0 for c in nodes):
        raise ValueError(f"{C4_NODES_ENV}={raw!r}: every node must be > 0")
    nodes = sorted(nodes)
    names = [f"{c:.4f}" for c in nodes]
    if len(set(names)) != len(names):
        raise ValueError(f"{C4_NODES_ENV}={raw!r}: node names collide at 4 decimals")
    if any(abs(float(n) - c) > 1e-12 * max(1.0, abs(c)) for n, c in zip(names, nodes)):
        raise ValueError(f"{C4_NODES_ENV}={raw!r}: nodes must round-trip at 4 decimals")
    return nodes


def _c4_node_json_path(c_pf, degree):
    """``c4-node-<pF:.4f>.json`` at degree 1 (``degree is None``), a ``-d2``
    sibling at degree 2; anything else raises."""
    if degree is None:
        return os.path.join(JSON_DIR, f"c4-node-{c_pf:.4f}.json")
    if degree == 2:
        return os.path.join(JSON_DIR, f"c4-node-{c_pf:.4f}-d2.json")
    raise ValueError(degree)


# ---------------------------------------------------------------------------
# 32-leg helpers.
# ---------------------------------------------------------------------------
def _slot_of(az_deg):
    return int(round(az_deg / AZIMUTH_STEP_DEG - 0.5)) % LEG_COUNT


def _slots_32(ids, geo):
    """``{pid: (slot k, ring sign)}`` from measured azimuth / z (geo[pid] = (n, az, z))."""
    return {pid: (_slot_of(geo[pid][1]), float(np.sign(geo[pid][2]))) for pid in ids}


def _ang(a):
    """Wrap to (-180, 180]."""
    return (a + 180.0) % 360.0 - 180.0


def _feeds_32(ids, geo, offset_slots=QUADRATURE_SLOT_OFFSET_32):
    """feed2 = the one port on DRIVEN's ring at az(DRIVEN) + offset·step within
    ``AZIMUTH_MATCH_DEG``.  Returns ``(feeds, az_feeds, az_bisector, d_az)``,
    ``d_az`` the measured feed2 − DRIVEN azimuth in [0, 360); raises unless
    exactly one candidate exists."""
    from tests.validation.test_port_birdcage_ring_column import AZIMUTH_MATCH_DEG

    _n, az0, z0 = geo[DRIVEN]
    target = az0 + offset_slots * AZIMUTH_STEP_DEG
    cands = [
        pid
        for pid in ids
        if pid != DRIVEN
        and np.sign(geo[pid][2]) == np.sign(z0)
        and abs(_ang(geo[pid][1] - target)) < AZIMUTH_MATCH_DEG
    ]
    if len(cands) != 1:
        raise AssertionError(
            f"(F) {len(cands)} candidates for feed2 at az {target % 360.0:.9f} deg on {DRIVEN}'s "
            f"ring: {cands}; measured azimuths "
            + ", ".join(f"{p} {geo[p][1]:.6f}" for p in ids if np.sign(geo[p][2]) == np.sign(z0))
        )
    feed2 = cands[0]
    d_az = (geo[feed2][1] - az0) % 360.0
    return [DRIVEN, feed2], (az0, geo[feed2][1]), (az0 + 45.0) % 360.0, d_az


def _k1_shares_32(currents, order, slots):
    """``{ring sign: k=±1 share}`` of the 32 slot currents per ring (Parseval)."""
    out = {}
    for sign in (1.0, -1.0):
        vec = np.zeros(LEG_COUNT, dtype=np.complex128)
        for idx, pid in enumerate(order):
            k, sg = slots[pid]
            if sg == sign:
                vec[k] = currents[idx]
        spec = np.fft.fft(vec)
        total = float(np.sum(np.abs(spec) ** 2))
        out[sign] = (
            float((abs(spec[1]) ** 2 + abs(spec[LEG_COUNT - 1]) ** 2) / total)
            if total
            else float("nan")
        )
    return out


def _esr_specs(base_specs, feeds, c_f, frequency_hz, z0, r_ohm=ESR_OHM):
    """Z_RC on the non-feed sheets, Z_th = z0·Z_RC/(z0 + Z_RC) on the feeds, κ off."""
    from tests.validation.test_wf10_f_human_tuned_circuit import _respec

    z_rc = h17._z_rc(c_f, frequency_hz, r_ohm=r_ohm)
    z_th = h17._thevenin_feed_impedance(z_rc, z0)
    imp = {sp.port_id: (z_th if sp.port_id in feeds else z_rc) for sp in base_specs}
    return _respec(base_specs, None, imp), z_rc, z_th


def _z_failures(specs, base_ids, feeds, c_f, frequency_hz, z0, r_ohm=ESR_OHM):
    """(Z): returns (failures, worst Im rel, worst feed rel)."""
    omega = 2.0 * np.pi * frequency_hz
    x_c = -1.0 / (omega * c_f)
    fails = []
    im_worst = 0.0
    th_worst = 0.0
    if [sp.port_id for sp in specs] != list(base_ids):
        fails.append("ids changed")
    for sp in specs:
        z = complex(sp.port_impedance_ohm)
        if sp.width_correction_kappa is not None:
            fails.append(f"{sp.port_id} kappa {sp.width_correction_kappa}")
        if sp.port_id in feeds:
            z_rc = complex(r_ohm, x_c)
            want = z0 * z_rc / (z0 + z_rc)
            th_worst = max(th_worst, abs(z - want) / abs(want))
        else:
            if z.real != r_ohm:
                fails.append(f"{sp.port_id} Re Z {z.real!r} != {r_ohm}")
            im_worst = max(im_worst, abs(z.imag / x_c - 1.0))
    if im_worst > 1e-12:
        fails.append(f"Im Z_RC rel {im_worst:.3e} > 1e-12")
    if th_worst > 1e-15:
        fails.append(f"Z_th rel {th_worst:.3e} > 1e-15")
    n_rc = sum(1 for sp in specs if sp.port_id not in feeds)
    if n_rc != len(specs) - 2:
        fails.append(f"{n_rc} RC specs")
    return fails, im_worst, th_worst


# ---------------------------------------------------------------------------
# Fixture.
# ---------------------------------------------------------------------------
def _apply_degree_c4(ctx, comm):
    """`_apply_degree`'s body under ``C4_DEGREE_ENV`` (H20a), printing the
    F-clinical records; (Deg) asserted inside.  Returns ``(solver, degree, counts)``;
    on the unset path ``ctx`` is not touched."""
    from tests.validation.test_geo38_f_clinical_build import F_CLINICAL_R1_UNKNOWNS_RECORD

    knob = _degree_knob_c4()
    if knob is None:
        solver, degree = ctx["solver"], 1
    else:
        from fem_em_solver.core import TimeHarmonicSolver

        problem = ctx["solver"].problem
        del ctx["solver"]
        ctx["solver"] = TimeHarmonicSolver(problem, degree=2)
        solver, degree = ctx["solver"], 2
    msh = solver.mesh
    v_space = solver.function_space()
    n_dofs = int(v_space.dofmap.index_map.size_global * v_space.dofmap.index_map_bs)
    for dim in (1, 2):
        if msh.topology.index_map(dim) is None:
            msh.topology.create_entities(dim)
    n_edges = int(msh.topology.index_map(1).size_global)
    n_faces = int(msh.topology.index_map(2).size_global)
    expected = n_edges if degree == 1 else 2 * n_edges + 2 * n_faces
    _say(
        comm,
        f"(Deg) knob {knob}; solver.degree {solver.degree}; global unknowns {n_dofs} vs "
        f"{'E' if degree == 1 else '2E+2F'} = {expected} (E {n_edges}, F {n_faces}; records "
        f"U1 {F_CLINICAL_R1_UNKNOWNS_RECORD[0]}, U2 {F_CLINICAL_R1_UNKNOWNS_RECORD[1]}; "
        f"ASSERTED exact)",
    )
    assert solver.degree == degree, (solver.degree, degree)
    assert n_dofs == expected, (n_dofs, expected, n_edges, n_faces)
    return solver, degree, {"n_dofs": n_dofs, "n_edges": n_edges, "n_faces": n_faces}


def _build_c4_fixture(comm):
    """`col._build_f_clinical`, then C3a's fixture asserts (copied from
    `_run_c3`'s fixture block, C3a's helpers by import), (U1), (F), then the
    degree knob.  Compiles nothing beyond `_context_from_mesh`'s width form."""
    from fem_em_solver.io.mesh import BIRDCAGE_CONDUCTOR_SURFACE_TAG
    from tests.mesh.test_birdcage_port_sheet_prerequisite import CELL_COUNT_BAND
    from tests.validation.test_geo38_f_clinical_build import (
        F_CLINICAL_R1_CELL_RECORD,
        F_CLINICAL_R1_UNKNOWNS_RECORD,
        _f_clinical_kwargs,
        _facet_coords,
        _owned,
        _ring_ports,
    )
    from tests.validation.test_ops61_exact_strip_ports import STRIP_BAND
    from tests.validation.test_port_birdcage_ring_column import AZIMUTH_MATCH_DEG
    from tests.validation.test_wf12_f_human_hole import _census_checks, _owned_census

    t0 = time.perf_counter()
    built = col._build_f_clinical(comm)
    comm.Barrier()
    t_build = time.perf_counter() - t0
    ctx = built["ctx"]
    msh, omega, cell_tags = ctx["msh"], ctx["omega"], ctx["cell_tags"]
    facets = built["facets"]
    n_cells = built["cells"]
    f_rel = omega / (2.0 * np.pi * FREQUENCY_HZ_C4) - 1.0
    census = _owned_census(msh, cell_tags, comm)
    census_fails = _census_checks(census, n_cells, _ring_ports())
    oc = built["outer_census"]
    cell_rel = n_cells / F_CLINICAL_R1_CELL_RECORD - 1.0
    R_s = float(built["shield_radius_m"])
    fv = np.asarray(facets.values)
    fi = np.asarray(facets.indices)
    cav_p = _facet_coords(msh, _owned(msh, fi[fv == BIRDCAGE_CONDUCTOR_SURFACE_TAG])).reshape(-1, 3)
    cav_hmax = comm.allreduce(
        float(np.hypot(cav_p[:, 0], cav_p[:, 1]).max()) if len(cav_p) else 0.0, op=MPI.MAX
    )
    _say(
        comm,
        f"fixture: t_build {t_build:.1f} s at -n {comm.size} (PREDICTED "
        f"{PREDICTED_T_BUILD_S[0]:.0f}-{PREDICTED_T_BUILD_S[1]:.0f} s, printed); (A0) cells "
        f"{n_cells} vs record {F_CLINICAL_R1_CELL_RECORD}, rel {cell_rel:.3e} (band "
        f"{CELL_COUNT_BAND}, ASSERTED; PREDICTED rel 0); (iv) failures "
        f"{census_fails or 'none'}; outer census {oc}; (S) cavity max hypot {cav_hmax:.9f} m "
        f"< R_s {R_s} (ASSERTED); f = omega/2pi rel to 64 MHz {f_rel:.1e} (ASSERTED <= 1e-12)",
    )
    assert abs(f_rel) <= 1e-12, f_rel
    assert abs(cell_rel) <= CELL_COUNT_BAND, cell_rel
    assert not census_fails, census_fails
    assert oc["tag_collision"] == 0 and oc["cavity"] > 0, oc
    assert (
        oc["cavity_exterior"] == oc["cavity"] and oc["outer"] + oc["cavity"] == oc["exterior"]
    ), oc
    assert cav_hmax < R_s, (cav_hmax, R_s)

    # ---- (C) census, (W) widths, (Az) layout, (U1) — C3a's helpers ------
    ids = [sp.port_id for sp in ctx["specs"]]
    geo = c3._sheet_geometry(ctx, comm)
    owned_counts = {pid: geo[pid][0] for pid in ids}
    c_fails = c3._census_failures(ids, owned_counts)
    kw = _f_clinical_kwargs("r1")
    w_expected = 2.0 * kw["ring_minor_radius"] + 2.0 * 1.0e-3
    assert abs(w_expected - c3.EXPECTED_SHEET_WIDTH_M) < 1e-15, w_expected
    w_rel = {sp.port_id: abs(sp.sheet_width_m / w_expected - 1.0) for sp in ctx["specs"]}
    w_worst = max(w_rel.values()) if w_rel else float("inf")
    az_fails, az_worst = c3._layout_failures(ids, geo)
    v_space = ctx["solver"].function_space()
    u1 = int(v_space.dofmap.index_map.size_global * v_space.dofmap.index_map_bs)
    _say(
        comm,
        f"(C) {len(ids)} specs (= 2*LEG_COUNT {2 * LEG_COUNT}), owned facets per sheet min "
        f"{min(owned_counts.values())} max {max(owned_counts.values())}; failures "
        f"{c_fails or 'none'} (ASSERTED)",
    )
    _say(
        comm,
        f"(W) sheet widths vs 2r+2c = {w_expected:.6f} m: worst rel {w_worst:.3e} "
        f"(band STRIP_BAND {STRIP_BAND:g}, ASSERTED)",
    )
    _say(
        comm,
        f"(Az) worst slot residual {az_worst:.3e} deg (< {AZIMUTH_MATCH_DEG:.0e}, ASSERTED); "
        f"layout failures {az_fails or 'none'}",
    )
    _say(
        comm,
        f"(U1) N1curl global size x bs = {u1} vs record {F_CLINICAL_R1_UNKNOWNS_RECORD[0]} "
        f"(EXACT, ASSERTED)",
    )
    assert not c_fails, c_fails
    assert w_worst <= STRIP_BAND, sorted(w_rel.items(), key=lambda kv: -kv[1])[:4]
    assert not az_fails, az_fails
    assert az_worst < AZIMUTH_MATCH_DEG, az_worst
    assert u1 == F_CLINICAL_R1_UNKNOWNS_RECORD[0], u1

    # ---- (F) the two-port feed -------------------------------------------
    feeds, az_feeds, az_bis, d_az = _feeds_32(ids, geo)
    f_res = abs(_ang(d_az - 90.0))
    _say(
        comm,
        f"(F) feeds {feeds}: az {az_feeds[0]:.9f} / {az_feeds[1]:.9f} deg, z sign "
        f"{np.sign(geo[feeds[0]][2]):+.0f} / {np.sign(geo[feeds[1]][2]):+.0f}; delta az "
        f"{d_az:.9f} deg (90 within {AZIMUTH_MATCH_DEG:.0e}, ASSERTED; one candidate "
        f"ASSERTED); bisector az(P33) + 45 = {az_bis:.9f} deg (PREDICTED feed2 "
        f"{PREDICTED_FEED2}, az(P33) {PREDICTED_AZ_P33_DEG}, bisector "
        f"{PREDICTED_BISECTOR_DEG}; printed)",
    )
    assert f_res < AZIMUTH_MATCH_DEG, (d_az, f_res)
    assert np.sign(geo[feeds[0]][2]) == np.sign(geo[feeds[1]][2])

    solver, degree, counts = _apply_degree_c4(ctx, comm)
    del solver
    return {
        "built": built,
        "ctx": ctx,
        "cells": n_cells,
        "t_build": t_build,
        "ids": ids,
        "geo": geo,
        "feeds": feeds,
        "az_bisector": az_bis,
        "degree": degree,
        "counts": counts,
        "u1": u1,
    }


# ---------------------------------------------------------------------------
# The build-only branch (H26a's K1 shape).
# ---------------------------------------------------------------------------
def _run_build(comm, monkeypatch):
    from tests.validation.test_geo38_f_clinical_build import F_CLINICAL_R1_UNKNOWNS_RECORD
    from tests.validation.test_port_package_sparameters import REFERENCE_IMPEDANCE_OHM
    from tests.validation.test_wf12_f_human_hole import _rss_gib

    t_window = time.perf_counter()
    z0 = float(REFERENCE_IMPEDANCE_OHM)
    monkeypatch.delenv(C4_DEGREE_ENV, raising=False)
    npz = c3.SAVE_PATH
    _say(comm, f"start: -n {comm.size}, BUILD-ONLY; npz {npz} never opened on this path")
    out = _build_c4_fixture(comm)
    ctx = out["ctx"]
    n1, e, f = (out["counts"][k] for k in ("n_dofs", "n_edges", "n_faces"))
    assert out["degree"] == 1 and ctx["solver"].degree == 1
    assert n1 == e == out["u1"], (n1, e, out["u1"])
    _m, rss_per_b, rss_build = _rss_gib(comm)

    # ---- (Z) the ESR series-RC sheets at the wiring value ----------------
    frequency_hz = ctx["omega"] / (2.0 * np.pi)
    base_specs = list(ctx["specs"])
    base_ids = [sp.port_id for sp in base_specs]
    c_f = SPEC_WIRING_C_PF * 1e-12
    specs, z_rc, z_th = _esr_specs(base_specs, out["feeds"], c_f, frequency_hz, z0)
    z_fails, im_worst, th_worst = _z_failures(specs, base_ids, out["feeds"], c_f, frequency_hz, z0)
    _say(
        comm,
        f"(Z) at C = {SPEC_WIRING_C_PF} pF (wiring value, not a node, not a prediction): "
        f"Z_RC {z_rc:.9e} Ohm on {len(specs) - 2} sheets (Re = {ESR_OHM} exact, Im rel "
        f"{im_worst:.1e} <= 1e-12), Z_th {z_th:.9e} Ohm on {out['feeds']} (rel {th_worst:.1e} "
        f"<= 1e-15); kappa None on all {len(specs)}; ids unchanged; failures "
        f"{z_fails or 'none'} (ASSERTED)",
    )
    assert not z_fails, z_fails
    assert len(specs) == 64

    # ---- (Deg-2) in the same build ---------------------------------------
    monkeypatch.setenv(C4_DEGREE_ENV, "2")
    t0 = time.perf_counter()
    solver2, deg2, counts2 = _apply_degree_c4(ctx, comm)
    t_d2 = time.perf_counter() - t0
    n2 = counts2["n_dofs"]
    assert counts2["n_edges"] == e and counts2["n_faces"] == f
    _m, rss_per, rss_sum = _rss_gib(comm)
    _say(
        comm,
        f"(Deg-2) knob 2: degree {deg2} (solver.degree {solver2.degree}); unknowns {n2} vs "
        f"2E+2F = {2 * e + 2 * f} vs record {F_CLINICAL_R1_UNKNOWNS_RECORD[1]} (ASSERTED "
        f"exact); measured F {f} beside predicted (U2 - 2 U1)/2 = "
        f"{(F_CLINICAL_R1_UNKNOWNS_RECORD[1] - 2 * F_CLINICAL_R1_UNKNOWNS_RECORD[0]) // 2}; "
        f"degree-2 space {t_d2:.1f} s (PREDICTED {PREDICTED_T_D2_S[0]:.0f}-"
        f"{PREDICTED_T_D2_S[1]:.0f} s, printed)",
    )
    _say(
        comm,
        f"NEG (asserted, rule (e)): n(2) {n2} != n(1) {n1}; n(2)/n(1) = {n2 / n1:.6f} vs "
        f"(2E+2F)/E = {(2 * e + 2 * f) / e:.6f} (exact integer identity ASSERTED)",
    )
    _say(
        comm,
        f"BUILD-ONLY return: t_build {out['t_build']:.1f} s, t_d2 {t_d2:.1f} s; RSS summed "
        f"after build {rss_build:.2f} GiB (PREDICTED ≈{PREDICTED_RSS_BUILD_GIB:g}), after the "
        f"degree-2 space per rank {', '.join(f'{r:.2f}' for r in (rss_per or []))} GiB, summed "
        f"{rss_sum:.2f} GiB (PREDICTED {PREDICTED_RSS_D2_GIB[0]:g}-{PREDICTED_RSS_D2_GIB[1]:g}; "
        f"printed); window {time.perf_counter() - t_window:.1f} s (no wall term, no assembly, "
        f"no factorisation)",
    )
    assert deg2 == 2 and solver2.degree == 2 and ctx["solver"] is solver2
    assert n2 == 2 * e + 2 * f, (n2, e, f)
    assert n2 == F_CLINICAL_R1_UNKNOWNS_RECORD[1], n2
    assert n2 != n1
    assert n2 * e == n1 * (2 * e + 2 * f), (n1, n2, e, f)


# ---------------------------------------------------------------------------
# The node branch (entry 46 only): H17's node body transposed to 64 ports.
# ---------------------------------------------------------------------------
def _run_nodes(comm):
    import ufl
    from dolfinx import fem

    from fem_em_solver.core.cavity import surface_resistance_ohm
    from fem_em_solver.ports.circuit import s_to_z
    from fem_em_solver.ports.lumped import (
        lumped_port_bilinear_term,
        lumped_port_linear_term,
        sheet_terminal_current,
    )
    from tests.validation.test_birdcage_power_identity import (
        DISCRETE_IDENTITY_RTOL,
        _sheet_field_dissipation_w,
        _source_power_w,
    )
    from tests.validation.test_port_birdcage_four_port import TERMINATED_PORT_IMPEDANCE_OHM
    from tests.validation.test_port_birdcage_lumped_column import PHANTOM_CELL_TAG
    from tests.validation.test_port_birdcage_ring_matrix import REUSE_REPRODUCTION_RTOL
    from tests.validation.test_port_circuit_layer_field import _residual
    from tests.validation.test_port_package_sparameters import REFERENCE_IMPEDANCE_OHM
    from tests.validation.test_th14_birdcage_copper import _leontovich_term, _surface_loss_w
    from tests.validation.test_wf10_f_human_tuned_circuit import _respec
    from tests.validation.test_wf12_f_human_hole import _rss_gib

    z0 = float(REFERENCE_IMPEDANCE_OHM)
    t_window = time.perf_counter()
    nodes_pf = _nodes_pf()
    knob = _degree_knob_c4()
    d2 = knob is not None

    # ---- the stored comparand (absent ⇒ stop, never re-solve) -------------
    npz_path = str(c3.SAVE_PATH)
    assert os.path.exists(
        npz_path
    ), f"{npz_path} absent: C4 stops (item 334's ruling (2)); never re-solve"
    with np.load(npz_path, allow_pickle=False) as npz:
        s64_npz = np.array(npz["s64"], dtype=np.complex128)
        order_npz = [str(p) for p in npz["order"]]
        cells_npz = int(npz["cells"])
    _say(
        comm,
        f"start: -n {comm.size}, nodes {', '.join(f'{c:.4f}' for c in nodes_pf)} pF, degree "
        f"{knob or 1} ({C4_DEGREE_ENV} {knob if d2 else 'unset'}), kappa off, r = {ESR_OHM} "
        f"Ohm on all 64 gaps; npz {npz_path}",
    )

    fixture = _build_c4_fixture(comm)
    assert fixture["degree"] == (knob or 1), (fixture["degree"], knob)
    built, t_build = fixture["built"], fixture["t_build"]
    feeds, geo, ids = fixture["feeds"], fixture["geo"], fixture["ids"]
    ctx = built["ctx"]
    msh, tags_f, omega, cell_tags = ctx["msh"], ctx["tags_f"], ctx["omega"], ctx["cell_tags"]
    frequency_hz = omega / (2.0 * np.pi)
    _say(comm, f"(A0) npz cells {cells_npz} (ASSERTED equal to {built['cells']})")
    assert built["cells"] == cells_npz, (built["cells"], cells_npz)

    z_s = (1.0 + 1.0j) * surface_resistance_ohm(omega, c3.COPPER_SIGMA)
    wall = _leontovich_term(msh, built["facets"], z_s, omega)
    z_p = complex(TERMINATED_PORT_IMPEDANCE_OHM)
    dx = ufl.Measure("dx", domain=msh, subdomain_data=cell_tags)
    acct = {"mesh": msh, "facet_tags": tags_f}
    solver = ctx["solver"]
    base_specs = list(ctx["specs"])

    # H17's `_drive` closure (`test_wf10_f_human_copper_esr_root.py:402–449`),
    # wall term first in the factorising call.
    def _drive(specs, driven_id, *, held, partition):
        sheets = [sp.sheet(driven=(sp.port_id == driven_id)) for sp in specs]
        dsheet = next(s for s in sheets if s.port_id == driven_id)
        linear = [
            lambda test, _s=dsheet: lumped_port_linear_term(
                msh, tags_f, _s, test, omega_rad_per_s=omega
            )
        ]
        comm.Barrier()
        t0 = time.perf_counter()
        if held:
            fields = solver.solve_with_held_factorization(extra_linear_terms=linear)
        else:
            bil = [wall] + [
                (
                    lambda trial, test, _s=s: lumped_port_bilinear_term(
                        msh, tags_f, _s, trial, test, omega_rad_per_s=omega
                    )
                )
                for s in sheets
            ]
            fields = solver.solve(
                current_density=None,
                project_source=False,
                extra_bilinear_terms=bil,
                extra_linear_terms=linear,
            )
        comm.Barrier()
        t_solve = time.perf_counter() - t0
        t1 = time.perf_counter()
        e = fields.e_complex
        currents = {s.port_id: sheet_terminal_current(msh, tags_f, s, e, comm) for s in sheets}
        out = {"currents": currents, "kind": "held" if held else "solve", "t_solve": t_solve}
        if partition:
            sigma_f = fields.sigma_field

            # F821: pyflakes sees the `del` at the end of this function; the closure runs before it.
            def _vol(measure):
                form = fem.form(0.5 * sigma_f * ufl.inner(e, e) * measure)  # noqa: F821
                return float(np.real(comm.allreduce(fem.assemble_scalar(form), op=MPI.SUM)))

            by_id = {
                s.port_id: float(_sheet_field_dissipation_w(acct, s, e, omega)) for s in sheets
            }
            out["p_sheet_by_id"] = by_id
            out["p_sheets"] = float(sum(by_id.values()))
            out["p_src"] = float(_source_power_w(acct, dsheet, e, omega))
            p_total = _vol(dx)
            out["p_phantom"] = _vol(dx(PHANTOM_CELL_TAG))
            out["p_rest"] = p_total - out["p_phantom"]
            out["p_surf"] = _surface_loss_w(msh, built["facets"], z_s, e, comm)
        del fields, e
        out["t_post"] = time.perf_counter() - t1
        return out

    # ---- order: [P33, feed2] + others by number --------------------------
    others = sorted((pid for pid in ids if pid not in feeds), key=lambda p: int(p[1:]))
    order = list(feeds) + others
    feeds_idx = [0, 1]
    slots = _slots_32(ids, geo)
    assert len(order) == 64 and sorted(order) == sorted(order_npz), "npz order differs in ids"
    npz_ix = [order_npz.index(p) for p in order]
    s64_ord = s64_npz[np.ix_(npz_ix, npz_ix)]
    _say(
        comm,
        f"feeds {feeds} (slots {slots[feeds[0]][0]}, {slots[feeds[1]][0]}; ring "
        f"{slots[feeds[0]][1]:+.0f}); bisector {fixture['az_bisector']:.6f} deg",
    )

    def _pair(imp_x, z_feed, v_feed):
        specs_x = _respec(base_specs, None, imp_x)
        assert all(sp.width_correction_kappa is None for sp in specs_x)
        cur, knd, tt, parts = [], [], [], []
        for n_drive, pid in enumerate(feeds):
            d = _drive(specs_x, pid, held=(n_drive > 0), partition=(n_drive == 0))
            knd.append(d["kind"])
            tt.append(d["t_solve"] + d["t_post"])
            cur.append(np.array([d["currents"][p] for p in order], dtype=np.complex128))
            if n_drive == 0:
                parts.append({k: v for k, v in d.items() if k.startswith("p_")})
            del d
        return h17._current_route_s(cur, feeds_idx, z_feed, v_feed), cur, knd, tt, parts[0]

    def _show(name, m):
        _say(
            comm,
            f"    {name} " + "; ".join(", ".join(f"{complex(v):.10e}" for v in row) for row in m),
        )

    def _partition(p):
        p_vol = p["p_phantom"] + p["p_rest"]
        rel = abs(p["p_src"] - p["p_sheets"] - p_vol - p["p_surf"]) / abs(p["p_src"])
        blind = abs(p["p_src"] - p["p_sheets"] - p_vol) / abs(p["p_src"])
        return rel, blind

    # ---- the Z50 pair: (P50), NEG, (Z50) ---------------------------------
    s_50, _c50, kinds_50, t_50, p50 = _pair({pid: z0 for pid in order}, z0, 1.0)
    del _c50
    rel_p50, blind50 = _partition(p50)
    d1_only = "degree-1 npz comparand, printed only at degree 2 (H26a)" if d2 else "ASSERTED"
    _say(comm, f"Z50 pair {kinds_50}, {', '.join(f'{t:.2f}' for t in t_50)} s; z_p {z_p}")
    _say(
        comm,
        f"(P50) P_src {p50['p_src']:.9e} sumP_sheet {p50['p_sheets']:.9e} P_phantom "
        f"{p50['p_phantom']:.9e} P_rest {p50['p_rest']:.9e} P_surf {p50['p_surf']:.9e} W; rel "
        f"{rel_p50:.3e} (band {DISCRETE_IDENTITY_RTOL:g}, ASSERTED)",
    )
    _say(
        comm,
        f"NEG P_surf-blind miss {blind50:.4e} = {blind50 / DISCRETE_IDENTITY_RTOL:.1f}x the "
        f"band (bar {c3.BLIND_FACTOR:.0f}x, ASSERTED)",
    )
    s_ff = s64_ord[:2, :2]
    r_z50 = _residual(s_50, s_ff)
    _show("S64npz[feeds] ", s_ff)
    _show("MEAS50        ", s_50)
    _say(
        comm,
        f"(Z50) ||MEAS50 - S64npz[feeds]||/||S64npz[feeds]|| = {r_z50:.6e} (band "
        f"{REUSE_REPRODUCTION_RTOL:g}, {d1_only}; PREDICTED <= {PREDICTED_Z50:g}); "
        f"t={time.perf_counter() - t_window:.1f} s",
    )

    # ---- per node --------------------------------------------------------
    z64_npz = s_to_z(s64_ord, z0)
    node_results, node_pair_s, y_by_node = [], [], []
    for c_pf in nodes_pf:
        c_f = c_pf * 1e-12
        z_rc = h17._z_rc(c_f, frequency_hz)
        z_c = h17._z_rc(c_f, frequency_hz, r_ohm=0.0)
        z_th = h17._thevenin_feed_impedance(z_rc, z0)
        v_s = h17._source_voltage(z_rc, z0)
        _say(
            comm,
            f"node C = {c_f:.12e} F ({c_pf:.4f} pF): Z_RC {z_rc:.9e} Ohm, Z_th {z_th:.9e} "
            f"Ohm, V_s {v_s:.9e} V",
        )
        imp = {pid: z_rc for pid in others}
        imp.update({pid: z_th for pid in feeds})
        s_node, cur_node, kinds_n, t_n, pn = _pair(imp, z_th, v_s)
        _show("S_node", s_node)
        rel_prc, blind_n = _partition(pn)
        p_rc = sum(v for pid, v in pn["p_sheet_by_id"].items() if pid not in feeds)
        p_feeds = sum(v for pid, v in pn["p_sheet_by_id"].items() if pid in feeds)
        p_vol = pn["p_phantom"] + pn["p_rest"]
        dz = abs(pn["p_src"] - (pn["p_sheets"] - p_rc) - p_vol - pn["p_surf"]) / abs(pn["p_src"])
        _say(
            comm,
            f"(P_RC) P_src {pn['p_src']:.9e} sumP_sheet {pn['p_sheets']:.9e} (62 RC "
            f"{p_rc:.9e}, 2 feeds {p_feeds:.9e}) P_phantom {pn['p_phantom']:.9e} P_rest "
            f"{pn['p_rest']:.9e} P_surf {pn['p_surf']:.9e} W; rel {rel_prc:.3e} (band "
            f"{DISCRETE_IDENTITY_RTOL:g}, ASSERTED; PREDICTED <= {PREDICTED_P_RC:g}); (Dz) "
            f"without the 62 RC sheets {dz:.4e}; P_surf-blind {blind_n:.4e} (printed)",
        )
        y_in = {}
        for label, a in h17.SENSES.items():
            y_in[label] = h17._active_y_in(s_node, a, z0)
            ya, yb = y_in[label]
            spread = abs(ya - yb) / abs(ya)
            _say(
                comm,
                f"Y_in a = {label}: {feeds[0]} {ya.real:+.9e} {ya.imag:+.9e}j S; {feeds[1]} "
                f"{yb.real:+.9e} {yb.imag:+.9e}j S; C4/mirror spread {spread:.4e} "
                f"(READING, never asserted)",
            )
        for lbl, z_load in (("Z_RC (r=0.2)", z_rc), ("Z_C (r=0)", z_c)):
            yc = h17._current_route_y_in(
                z64_npz, feeds_idx, z_load, z_load, h17.SENSES["(1,+1j)"], z0
            )
            _say(
                comm,
                f"circuit-only (degree-1 npz s64, local loads {lbl}; not a prediction): "
                f"Y_in (1,+1j) {feeds[0]} {yc[0]:.9e} S, {feeds[1]} {yc[1]:.9e} S",
            )
        sh = _k1_shares_32(cur_node[0], order, slots)
        ring = slots[feeds[0]][1]
        _say(
            comm,
            f"k=+-1 share (32-point), {feeds[0]} column: feed ring {sh[ring]:.6f} (ring+ "
            f"{sh[1.0]:.6f} ring- {sh[-1.0]:.6f}); node pair {kinds_n}, "
            f"{', '.join(f'{t:.2f}' for t in t_n)} s; t={time.perf_counter() - t_window:.1f} s",
        )
        rec = {
            "c_f": c_f,
            "r_ohm": ESR_OHM,
            "freq_hz": float(frequency_hz),
            "order": order,
            "degree": knob or 1,
            "s_node": [[[complex(v).real, complex(v).imag] for v in row] for row in s_node],
            "y_in": {lbl: [[y.real, y.imag] for y in ys] for lbl, ys in y_in.items()},
            "partition": dict(pn),
            "rel_p_rc": rel_prc,
            "dz": dz,
            "k1_share_feed_ring": sh[ring],
            "times_s": {"pair": t_n, "z50_pair": t_50, "build": t_build},
        }
        if comm.rank == 0:
            os.makedirs(JSON_DIR, exist_ok=True)
            path = _c4_node_json_path(c_pf, knob)
            with open(path, "w") as fh:
                json.dump(rec, fh, indent=1)
            print(f"[GEO-38 C4] wrote {path}", flush=True)
        comm.Barrier()
        node_results.append((c_pf, kinds_n, rel_prc))
        node_pair_s.append(sum(t_n))
        y_by_node.append(y_in)

    # ---- the root reader (printed; C_in is a record, never asserted) -----
    c_arr = np.array(nodes_pf) * 1e-12
    for lbl in h17.SENSES:
        for k, name in enumerate(feeds):
            b = np.array([y[lbl][k].imag for y in y_by_node])
            found = h17._sign_change_roots(c_arr, b) if len(c_arr) > 1 else []
            _say(
                comm,
                f"{name} a = {lbl}: Im Y_in "
                + ", ".join(f"{x * 1e12:.4f}:{v:+.6e}" for x, v in zip(c_arr, b))
                + " S; sign-change roots "
                + (
                    ", ".join(
                        f"[{a * 1e12:.4f}, {bb * 1e12:.4f}] -> {r * 1e12:.6f} pF"
                        for a, bb, r in found
                    )
                    or "none"
                )
                + (
                    " (C_in,F-clinical: RECORD, never asserted)"
                    if (name == DRIVEN and lbl == "(1,+1j)")
                    else ""
                ),
            )

    _m, rss_per, rss_sum = _rss_gib(comm)
    _say(
        comm,
        f"price: build {t_build:.1f} s, Z50 pair {sum(t_50):.1f} s, node pairs "
        f"{', '.join(f'{t:.1f}' for t in node_pair_s)} s; RSS per rank "
        f"{', '.join(f'{r:.2f}' for r in (rss_per or []))} GiB, summed {rss_sum:.2f} GiB; "
        f"window {time.perf_counter() - t_window:.1f} s at -n {comm.size}",
    )

    assert kinds_50 == ["solve", "held"], kinds_50
    assert rel_p50 <= DISCRETE_IDENTITY_RTOL, rel_p50
    assert blind50 >= c3.BLIND_FACTOR * DISCRETE_IDENTITY_RTOL, blind50
    if not d2:
        assert r_z50 <= REUSE_REPRODUCTION_RTOL, f"(Z50) closure {r_z50:.6e}"
    for c_pf, kinds_n, rel_prc in node_results:
        assert kinds_n == ["solve", "held"], kinds_n
        assert (
            rel_prc <= DISCRETE_IDENTITY_RTOL
        ), f"(P_RC) at {c_pf:.4f} pF: {rel_prc:.3e} — never re-band"


# ---------------------------------------------------------------------------
# Coding tests (numpy only); skip when FEM_EM_GEO38_C4 is set.
# ---------------------------------------------------------------------------
_coding = pytest.mark.skipif(
    os.environ.get(C4_ENV, "").strip() != "",
    reason=f"{C4_ENV} set: coding tests skip in build/solve windows",
)


def _ideal_geo():
    az, zc = c3._ideal_layout()
    return {f"P{33 + k}": (1, float(az[k]), float(zc[k])) for k in range(64)}


@_coding
def test_c4_k0_knob_parsing(monkeypatch):
    """(K0): all three knobs."""
    monkeypatch.delenv(C4_ENV, raising=False)
    assert _mode() is None
    for good in ("build", "1"):
        monkeypatch.setenv(C4_ENV, good)
        assert _mode() == good
    for bad in ("2", "x"):
        monkeypatch.setenv(C4_ENV, bad)
        with pytest.raises(ValueError):
            _mode()
    monkeypatch.delenv(C4_NODES_ENV, raising=False)
    with pytest.raises(ValueError):
        _nodes_pf()
    monkeypatch.setenv(C4_NODES_ENV, "60.25,58.5, 59.1234")
    assert _nodes_pf() == [58.5, 59.1234, 60.25]
    for bad in (
        "",
        " , ",
        "-1.0,58.5",
        "0,58.5",
        "58.50001,58.50002",
        "58.5,58.5",
        "58.12345",
        "x",
    ):
        monkeypatch.setenv(C4_NODES_ENV, bad)
        with pytest.raises(ValueError):
            _nodes_pf()
    monkeypatch.delenv(C4_DEGREE_ENV, raising=False)
    assert _degree_knob_c4() is None
    monkeypatch.setenv(C4_DEGREE_ENV, "2")
    assert _degree_knob_c4() == 2
    for bad in ("1", "x"):
        monkeypatch.setenv(C4_DEGREE_ENV, bad)
        with pytest.raises(ValueError):
            _degree_knob_c4()
    print(
        "[GEO-38 C4 K0] knobs: C4 unset/build/1 accepted, 2/x refused; nodes required, "
        "sorted, >0, distinct at 4 decimals; degree unset/2 accepted, 1/x refused"
    )


@_coding
def test_c4_k1_feeds_on_ideal_layout():
    """(K1): P33 and P41, Δaz = 90° ≤ 1e-12, bisector 50.625°."""
    geo = _ideal_geo()
    feeds, az_f, bis, d_az = _feeds_32(list(geo), geo)
    print(f"[GEO-38 C4 K1] feeds {feeds}, az {az_f}, delta {d_az:.15f}, bisector {bis:.15f}")
    assert QUADRATURE_SLOT_OFFSET_32 == 8
    assert feeds == ["P33", "P41"], feeds
    assert abs(d_az - 90.0) <= 1e-12, d_az
    assert abs(bis - 50.625) <= 1e-12, bis
    assert abs(az_f[0] - 5.625) <= 1e-12


@_coding
def test_c4_k2_k1_shares_32():
    """(K2): a planted k = 1 pattern reads 1.0, a planted k = 2 pattern reads 0."""
    geo = _ideal_geo()
    ids = list(geo)
    slots = _slots_32(ids, geo)
    out = {}
    for m in (1, 2):
        cur = np.array([np.exp(2j * np.pi * m * slots[p][0] / LEG_COUNT) for p in ids])
        out[m] = _k1_shares_32(cur, ids, slots)
    print(f"[GEO-38 C4 K2] planted k=1 -> {out[1]}; k=2 -> {out[2]}")
    for sign in (1.0, -1.0):
        assert abs(out[1][sign] - 1.0) <= 1e-12, out[1]
        assert abs(out[2][sign]) <= 1e-12, out[2]


@_coding
def test_c4_k3_node_json_paths():
    """(K3): the -d2 sibling differs, both under JSON_DIR, never an h17-node name."""
    for c_pf in (58.5, 100.0, 120.1234):
        p1, p2 = _c4_node_json_path(c_pf, None), _c4_node_json_path(c_pf, 2)
        assert p1 != p2 and p2.endswith("-d2.json")
        assert os.path.dirname(p1) == JSON_DIR == os.path.dirname(p2)
        for p in (p1, p2):
            assert os.path.basename(p).startswith("c4-node-")
            assert p not in (h17._h17_node_json_path(c_pf, None), h17._h17_node_json_path(c_pf, 2))
    with pytest.raises(ValueError):
        _c4_node_json_path(58.5, 1)
    print(f"[GEO-38 C4 K3] {_c4_node_json_path(58.5, None)} / {_c4_node_json_path(58.5, 2)}")


@_coding
def test_c4_y0_64_routes_agree():
    """(Y0-64): H17's (Y0) at 64 ports — S-route Y_in = current-route Y_in
    within ``h17.Y0_CODING_BAND``."""
    from fem_em_solver.ports.circuit import z_to_s

    z0 = 50.0
    f = FREQUENCY_HZ_C4
    rng = np.random.default_rng(20261020)
    n = 64
    b = rng.standard_normal((n, n))
    r_part = 5.0 * (b @ b.T) / n + 0.5 * np.eye(n)
    x = rng.standard_normal((n, n))
    z = r_part + 1j * 40.0 * (x + x.T)
    assert np.min(np.linalg.eigvalsh(0.5 * (z + z.conj().T))) > 0.0
    s = z_to_s(z, z0)
    z_rc = h17._z_rc(SPEC_WIRING_C_PF * 1e-12, f)
    s_rt = h17._s_route_2x2(s, [0, 1], z_rc, z_rc, z0)
    s_cm = h17._circuit_current_route_s(z, [0, 1], z_rc, z_rc, z0)
    r_s = h17._rel(s_cm, s_rt)
    worst = 0.0
    for a in h17.SENSES.values():
        y_s = h17._active_y_in(s_rt, a, z0)
        y_c = h17._current_route_y_in(z, [0, 1], z_rc, z_rc, a, z0)
        worst = max(worst, h17._rel(y_s, y_c))
    print(
        f"[GEO-38 C4 Y0-64] S routes {r_s:.3e}, Y_in routes {worst:.3e} "
        f"(band {h17.Y0_CODING_BAND:g}, ASSERTED)"
    )
    assert r_s <= h17.Y0_CODING_BAND, r_s
    assert worst <= h17.Y0_CODING_BAND, worst


@_coding
def test_c4_z_specs_and_planted_r0():
    """(Z) on synthetic specs, and the planted control: r = 0 gives Re Z_RC = 0 ≠ 0.2."""
    from fem_em_solver.ports.lumped import LumpedSheetPortSpec

    base = [
        LumpedSheetPortSpec(
            port_id=f"P{i}",
            facet_tag=1000 + i,
            port_impedance_ohm=50.0,
            gap_height_m=0.01,
            sheet_width_m=0.01,
            drive_direction=(0.0, 0.0, 1.0),
        )
        for i in range(33, 97)
    ]
    ids = [sp.port_id for sp in base]
    c_f = SPEC_WIRING_C_PF * 1e-12
    specs, z_rc, z_th = _esr_specs(base, ["P33", "P41"], c_f, FREQUENCY_HZ_C4, 50.0)
    fails, im_w, th_w = _z_failures(specs, ids, ["P33", "P41"], c_f, FREQUENCY_HZ_C4, 50.0)
    specs0, z_rc0, _z = _esr_specs(base, ["P33", "P41"], c_f, FREQUENCY_HZ_C4, 50.0, r_ohm=0.0)
    fails0, _a, _b = _z_failures(specs0, ids, ["P33", "P41"], c_f, FREQUENCY_HZ_C4, 50.0)
    print(
        f"[GEO-38 C4 Z] Z_RC {z_rc:.9e}, Z_th {z_th:.9e}, fails {fails}; planted r=0: Re Z_RC "
        f"{z_rc0.real} (ASSERTED != {ESR_OHM}), {len(fails0)} failures"
    )
    assert not fails, fails
    assert z_rc.real == ESR_OHM
    assert z_rc0.real == 0.0 and z_rc0.real != ESR_OHM
    assert len(fails0) >= 62


@_coding
def test_c4_n1_f_human_offset_misses_quadrature():
    """(N1), asserted by construction: F-human's QUADRATURE_SLOT_OFFSET (4)
    on the 32-leg layout puts feed2 at Δaz = 45°, missing 90° by exactly 45°.
    The F-human constant is imported here only, as the negative control."""
    from tests.validation.test_port_birdcage_ring_column import AZIMUTH_MATCH_DEG
    from tests.validation.test_wf10_f_human_tuned_circuit import QUADRATURE_SLOT_OFFSET

    geo = _ideal_geo()
    feeds, _az, _bis, d_az = _feeds_32(list(geo), geo, offset_slots=QUADRATURE_SLOT_OFFSET)
    miss = abs(_ang(d_az - 90.0))
    print(
        f"[GEO-38 C4 N1] F-human offset {QUADRATURE_SLOT_OFFSET}: feed2 {feeds[1]}, delta az "
        f"{d_az:.12f} deg, miss {miss:.12f} deg (ASSERTED 45, > {AZIMUTH_MATCH_DEG:.0e})"
    )
    assert QUADRATURE_SLOT_OFFSET == 4
    assert abs(miss - 45.0) <= 1e-12, miss
    assert miss > AZIMUTH_MATCH_DEG


@_coding
def test_c4_n2_f_human_k1_reader_refuses_32_slots():
    """(N2), asserted: F-human's 16-point `_k1_shares` on 32-slot input raises
    IndexError.  Imported here only, as the negative control."""
    from tests.validation.test_wf10_f_human_tuned_circuit import _k1_shares

    geo = _ideal_geo()
    ids = list(geo)
    slots = _slots_32(ids, geo)
    cur = np.ones(len(ids), dtype=np.complex128)
    with pytest.raises(IndexError):
        _k1_shares(cur, ids, slots)
    print("[GEO-38 C4 N2] F-human _k1_shares on 32 slots raised IndexError (ASSERTED)")


# ---------------------------------------------------------------------------
# The gated test.
# ---------------------------------------------------------------------------
@complex_only
def test_geo38_c4_esr_root(monkeypatch):
    mode = _mode()
    if mode is None:
        pytest.skip(f"{C4_ENV} unset: GEO-38 C4 is selected by environment only")
    comm = MPI.COMM_WORLD
    if mode == "build":
        _run_build(comm, monkeypatch)
        return
    _run_nodes(comm)
