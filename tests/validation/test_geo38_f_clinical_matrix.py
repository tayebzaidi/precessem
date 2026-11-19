"""GEO-38 C3a (§9 item 519, §10 chain C3): the F-clinical-1.5T 64×64 copper
sweep module — `PORT-13`'s ring-port sweep with `PORT-19`'s held factor,
transposed from `WF-12` step 2's copper 32×32
(`test_wf12_f_human_copper_matrix.py:86–382`) to 32 legs / 64 ring sheets.

Selected by ``FEM_EM_GEO38_C3`` only (never by ``-k``):

* unset — the gated test skips; the numpy coding tests (K0)–(K2), N1–N3 run;
* ``build`` — build-only proof: the column module's ``_build_f_clinical``
  (entry 40's in-memory route), then (A0), (iv), the outer census, (S), (C),
  (W), (Az) and (U1); returns before the Leontovich wall term, any assembly
  or factorisation;
* ``1`` — the full 64-drive sweep (``xl-pending.md`` entry 44 only);
* anything else — ``ValueError``.

``FEM_EM_GEO38_C3_DRIVES`` (1..64, default 64): below 64 the matrix gates
print SKIPPED (`WF-12`'s ``_knobs`` pattern).

``FEM_EM_GEO38_C3_DEGREE`` (§9 item 705, C5c; H26a's and C4a's pattern): unset /
empty ⇒ degree 1, ``2`` ⇒ the degree-2 solver swapped in after (U1), (Deg)
asserted, and ``_save_path()`` ⇒ ``SAVE_PATH_D2``.  Claimed by item 705: the
knob, the save path and the count.  ``_sweep`` has never run at degree 2.

``FEM_EM_GEO38_C3_RUNG`` (§9 item 889, C35): unset / empty / ``r1`` ⇒ r1;
``r2`` ⇒ ``_r2_patches`` rebinds s1a's kwargs and R1 records to r2's and
``SAVE_PATH`` to ``SAVE_PATH_R2`` (degree 1 only; r2 with the degree knob set is
refused before gmsh).  Read only by the gated test, never by ``_run_c3``.
Claimed by item 889: the knob, the save path and the r2 build-only record.

Since §9 item 706 (C17) ``build`` asserts the save path's state unchanged, a
before / after identity (``_build_keeps_save_state``), not its absence: entry 44
wrote the degree-1 npz, which turned the old assert red by construction.

Counts derived from ``LEG_COUNT`` = 32 (not bands): step 11.25°, 17 azimuth
classes, 34 symmetry classes.  F-human's 16-leg constants (``_azimuth_class``,
``N_SYMMETRY_CLASSES`` = 18, ``MATRIX_CONTROL_MARGIN``) are never used for
classification here; N1 asserts why.

Bands, all imported: ``CELL_COUNT_BAND``, ``STRIP_BAND``, ``AZIMUTH_MATCH_DEG``,
``DISCRETE_IDENTITY_RTOL``, ``RECIPROCITY_BAND``, ``MATRIX_PASSIVITY_CEILING``,
``OPPOSITE_SPREAD_BAND``.  Entry 44's readout names ``PASSIVITY_SIGMA_TOLERANCE``
and ``ADJACENT_SPREAD_BAND``; the worst spread is printed against the latter,
never asserted (brief "Band reconciliation", ruling (e)).

Scope: the module, its coding tests and a build-only proof.  No S entry, no
reciprocity/passivity/class reading on the fixture is claimed by item 519.
"""

from __future__ import annotations

import os
import time
from pathlib import Path

import numpy as np
import pytest
from mpi4py import MPI

from tests.complex_mode import complex_only

C3_ENV = "FEM_EM_GEO38_C3"
C3_DRIVES_ENV = "FEM_EM_GEO38_C3_DRIVES"
C3_DEGREE_ENV = "FEM_EM_GEO38_C3_DEGREE"
C3_RUNG_ENV = "FEM_EM_GEO38_C3_RUNG"
COPPER_SIGMA = 5.8e7
DRIVEN = "P33"
N_PORTS = 64
SAVE_PATH = Path("logs/geo38/c3-cu-s64-d1.npz")
SAVE_PATH_D2 = Path("logs/geo38/c3-cu-s64-d2.npz")
SAVE_PATH_R2 = Path("logs/geo38/c3-cu-s64-d1-r2.npz")
# (S-r2): entry 44's and entry 54's records under their literal paths (never rebound).
ENTRY44_D1_PATH = Path("logs/geo38/c3-cu-s64-d1.npz")
ENTRY54_D2_PATH = Path("logs/geo38/c3-cu-s64-d2.npz")
# Printed prediction (rule (e)): C5c's r1 t_build 338.6 s × s1a's r2 / r1 build 245.4 / 323.3.
PREDICTED_T_BUILD_R2_S = (338.6 * 245.4 / 323.3, (210.0, 400.0))
EXPECTED_SHEET_WIDTH_M = 0.010  # 2·ring_minor_radius (4 mm) + 2·port_clearance (1 mm)
BLIND_FACTOR = 1.0  # binary NEG: the identity with P_surf dropped misses the band

# Counts derived from LEG_COUNT (not bands).
LEG_COUNT_32 = 32
AZIMUTH_STEP_32_DEG = 360.0 / LEG_COUNT_32
N_AZIMUTH_CLASSES_32 = LEG_COUNT_32 // 2 + 1
N_SYMMETRY_CLASSES_32 = 2 * N_AZIMUTH_CLASSES_32

# N2 ceiling, computed first: equal-norm columns, one column ×1.01 lifts the
# reciprocity ratio by ≈ 0.01·√2/√64 = 1.77e-3 = 1.77× RECIPROCITY_BAND.
N2_CEILING = 0.01 * np.sqrt(2.0) / np.sqrt(N_PORTS)

# Printed predictions for the full branch (rule (e); none measured on this fixture).
PREDICTED_RECIPROCITY = 1e-14
PREDICTED_SPREAD_F_HUMAN = 0.0327  # `WF-12.md:142`, F-human copper


def _say(comm, msg):
    if comm.rank == 0:
        print(f"[GEO-38 C3] {msg}", flush=True)


def _mode():
    raw = os.environ.get(C3_ENV, "").strip()
    if raw == "":
        return None
    if raw in ("build", "1"):
        return raw
    raise ValueError(f"{C3_ENV} must be unset, 'build' or '1', got {raw!r}")


def _n_drives():
    raw = os.environ.get(C3_DRIVES_ENV, "").strip()
    n = int(raw) if raw else N_PORTS
    if not 1 <= n <= N_PORTS:
        raise ValueError(f"{C3_DRIVES_ENV} must be in 1..{N_PORTS}, got {raw!r}")
    return n


def _degree_knob_c3():
    """Unset / empty ⇒ ``None`` (degree 1); ``"2"`` ⇒ 2; anything else ⇒ ValueError."""
    raw = os.environ.get(C3_DEGREE_ENV, "").strip()
    if raw == "":
        return None
    if raw == "2":
        return 2
    raise ValueError(f"{C3_DEGREE_ENV}={raw!r}: only unset/empty (degree 1) or '2' are accepted")


def _rung_knob_c3():
    """Unset / empty / ``"r1"`` ⇒ ``"r1"``; ``"r2"`` ⇒ ``"r2"``; anything else ⇒ ValueError."""
    raw = os.environ.get(C3_RUNG_ENV, "").strip()
    if raw in ("", "r1"):
        return "r1"
    if raw == "r2":
        return "r2"
    raise ValueError(f"{C3_RUNG_ENV}={raw!r}: only unset/empty, 'r1' or 'r2' are accepted")


def _rung_checked():
    """The rung, with r2 refused at degree 2 before gmsh (C35 builds no r2 degree-2 path)."""
    rung = _rung_knob_c3()
    if rung == "r2" and os.environ.get(C3_DEGREE_ENV, "").strip() != "":
        raise ValueError(f"{C3_RUNG_ENV}=r2 with {C3_DEGREE_ENV} set: r2 is degree 1 only (C35)")
    return rung


def _r2_patches(monkeypatch):
    """The r2 rung through `_run_c3`'s route (§9 item 889, C35; `_c11_patches`'
    pattern): rebind the s1a names `_run_c3` / `_build_f_clinical` /
    `_apply_degree_c3` import at call time, and this module's ``SAVE_PATH``."""
    import sys

    import tests.validation.test_geo38_f_clinical_build as s1a

    c3 = sys.modules[__name__]
    orig = s1a._f_clinical_kwargs
    monkeypatch.setattr(s1a, "_f_clinical_kwargs", lambda _rung: orig("r2"))
    monkeypatch.setattr(s1a, "F_CLINICAL_R1_CELL_RECORD", s1a.F_CLINICAL_R2_CELL_RECORD)
    monkeypatch.setattr(s1a, "F_CLINICAL_R1_UNKNOWNS_RECORD", s1a.F_CLINICAL_R2_UNKNOWNS_RECORD)
    monkeypatch.setattr(c3, "SAVE_PATH", SAVE_PATH_R2)


def _save_path():
    """``SAVE_PATH_D2`` at degree 2; otherwise the module global ``SAVE_PATH``,
    read at call time: C11's ``_c11_patches`` rebinds it."""
    return SAVE_PATH_D2 if _degree_knob_c3() == 2 else SAVE_PATH


def _check_unknowns(n_dofs, degree, n_edges, n_faces):
    """(Deg): E at degree 1, 2E + 2F at degree 2, exact.  Returns the expected count."""
    expected = n_edges if degree == 1 else 2 * n_edges + 2 * n_faces
    assert n_dofs == expected, (n_dofs, expected, n_edges, n_faces)
    return expected


def _apply_degree_c3(ctx, comm):
    """`_apply_degree_c4`'s body under ``C3_DEGREE_ENV``, (Deg) asserted inside.
    Returns ``(solver, degree, counts)``; on the unset path ``ctx`` is not touched."""
    from tests.validation.test_geo38_f_clinical_build import F_CLINICAL_R1_UNKNOWNS_RECORD

    knob = _degree_knob_c3()
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
    assert solver.degree == degree, (solver.degree, degree)
    expected = _check_unknowns(n_dofs, degree, n_edges, n_faces)
    _say(
        comm,
        f"(Deg) knob {knob}; solver.degree {solver.degree}; global unknowns {n_dofs} vs "
        f"{'E' if degree == 1 else '2E+2F'} = {expected} (E {n_edges}, F {n_faces}; records "
        f"U1 {F_CLINICAL_R1_UNKNOWNS_RECORD[0]}, U2 {F_CLINICAL_R1_UNKNOWNS_RECORD[1]}; "
        f"ASSERTED exact)",
    )
    assert n_dofs == F_CLINICAL_R1_UNKNOWNS_RECORD[degree - 1], (
        n_dofs,
        F_CLINICAL_R1_UNKNOWNS_RECORD,
    )
    return solver, degree, {"n_dofs": n_dofs, "n_edges": n_edges, "n_faces": n_faces}


def _azimuth_class_32(az_i, az_j):
    """``(steps, residual_deg)`` with the 32-leg 11.25° step (fold of
    `test_port_birdcage_ring_matrix.py:206–215`)."""
    delta = abs((az_i - az_j + 180.0) % 360.0 - 180.0)
    steps = int(round(delta / AZIMUTH_STEP_32_DEG))
    return steps, abs(delta - steps * AZIMUTH_STEP_32_DEG)


def _expected_class_sizes():
    sizes = {}
    for same in (True, False):
        for k in range(N_AZIMUTH_CLASSES_32):
            sizes[(same, k)] = 64 if k in (0, LEG_COUNT_32 // 2) else 128
    return sizes


def _class_table(az, zc, classify=_azimuth_class_32):
    """{(same_ring, steps): [(i, j), ...]} and the worst residual."""
    n = len(az)
    classes = {}
    worst = 0.0
    for i in range(n):
        for j in range(n):
            steps, res = classify(az[i], az[j])
            worst = max(worst, res)
            classes.setdefault((bool(zc[i] * zc[j] > 0.0), steps), []).append((i, j))
    return classes, worst


def _class_spreads(s, classes):
    out = {}
    for key, pairs in classes.items():
        m = np.array([abs(s[i, j]) for i, j in pairs])
        out[key] = (float((m.max() - m.min()) / m.mean()), len(m), float(m.mean()))
    return out


def _ideal_layout():
    """Ideal azimuths (j + ½)·11.25° and z signs: bottom ring (P33…P64) first."""
    az = np.array([(j + 0.5) * AZIMUTH_STEP_32_DEG for j in range(LEG_COUNT_32)] * 2)
    zc = np.array([-1.0] * LEG_COUNT_32 + [1.0] * LEG_COUNT_32)
    return az, zc


def _synthetic_s(az, zc):
    """Reciprocal C32 × mirror-symmetric S_ij = f(ring relation, steps)."""
    n = len(az)
    s = np.zeros((n, n), dtype=complex)
    for i in range(n):
        for j in range(n):
            steps, _r = _azimuth_class_32(az[i], az[j])
            same = zc[i] * zc[j] > 0.0
            s[i, j] = (
                (0.08 if same else 0.06)
                * np.exp(-0.15 * steps)
                * np.exp(0.3j * steps + (0.0 if same else 0.7j))
            )
    return s


def _census_failures(spec_ids, owned_counts):
    """(C): returns a list of failure strings (empty ⇒ green)."""
    from tests.validation.test_geo38_f_clinical_build import LEG_COUNT, _ring_ports

    expected = [f"P{i}" for i in _ring_ports()]
    fails = []
    if len(spec_ids) != N_PORTS or N_PORTS != 2 * LEG_COUNT:
        fails.append(f"{len(spec_ids)} specs != {N_PORTS} = 2*{LEG_COUNT}")
    if sorted(spec_ids) != sorted(expected) or len(set(spec_ids)) != len(spec_ids):
        fails.append(f"spec ids are not exactly P33..P96: {sorted(set(expected) ^ set(spec_ids))}")
    empty = [pid for pid in spec_ids if owned_counts.get(pid, 0) <= 0]
    if empty:
        fails.append(f"sheets with no owned facets: {empty}")
    return fails


def _sheet_geometry(ctx, comm):
    """Per spec: (global owned facet count, centroid azimuth deg, centroid z),
    from the owned facets of ``ctx["tags_f"]``, area-weighted, reduced once."""
    from tests.validation.test_geo38_f_clinical_build import _area, _facet_coords, _owned

    msh, tags_f = ctx["msh"], ctx["tags_f"]
    fv = np.asarray(tags_f.values)
    fi = np.asarray(tags_f.indices)
    local = []
    for sp in ctx["specs"]:
        xyz = _facet_coords(msh, _owned(msh, fi[fv == sp.facet_tag]))
        a = _area(xyz)
        c = xyz.mean(axis=1) if len(xyz) else np.zeros((0, 3))
        local.append(
            [float(len(a)), float(a.sum()), *(a[:, None] * c).sum(axis=0).tolist()]
            if len(a)
            else [0.0, 0.0, 0.0, 0.0, 0.0]
        )
    tot = comm.allreduce(np.array(local), op=MPI.SUM)
    geo = {}
    for sp, row in zip(ctx["specs"], tot):
        n, area, mx, my, mz = row
        cx, cy, cz = (mx / area, my / area, mz / area) if area > 0 else (np.nan,) * 3
        geo[sp.port_id] = (int(n), float(np.degrees(np.arctan2(cy, cx)) % 360.0), float(cz))
    return geo


def _layout_failures(ids, geo):
    """(Az): measured azimuths on (j + ½)·step, z sign per ring, 32 distinct
    slots per ring.  Returns (failures, worst slot residual deg)."""
    fails = []
    worst = 0.0
    slots = {-1: [], 1: []}
    for pid in ids:
        _n, az, z = geo[pid]
        j = int(round(az / AZIMUTH_STEP_32_DEG - 0.5)) % LEG_COUNT_32
        res = abs(((az - (j + 0.5) * AZIMUTH_STEP_32_DEG) + 180.0) % 360.0 - 180.0)
        worst = max(worst, res)
        ordinal = int(pid[1:])
        want = -1 if ordinal <= 2 * LEG_COUNT_32 else 1  # build order: bottom ring first
        if np.sign(z) != want:
            fails.append(f"{pid} z {z:+.6e} sign != {want:+d}")
        slots[want].append(j)
    for ring, js in slots.items():
        if sorted(js) != list(range(LEG_COUNT_32)):
            fails.append(f"ring z{ring:+d} slots not one-per-gap: {sorted(js)}")
    return fails, worst


def _run_c3(comm, *, build_only, n_drives=N_PORTS):
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
    from tests.validation.test_geo38_f_clinical_column import _build_f_clinical
    from tests.validation.test_ops61_exact_strip_ports import STRIP_BAND
    from tests.validation.test_port_birdcage_ring_column import AZIMUTH_MATCH_DEG
    from tests.validation.test_wf12_f_human_hole import _census_checks, _owned_census, _rss_gib

    t_window = time.perf_counter()
    _say(
        comm,
        f"start: -n {comm.size}, {'BUILD-ONLY' if build_only else f'full sweep, {n_drives} drives'}",
    )

    # ---- fixture (entry 40's in-memory route), (A0), (iv), outer census, (S)
    t0 = time.perf_counter()
    built = _build_f_clinical(comm)
    comm.Barrier()
    t_build = time.perf_counter() - t0
    ctx = built["ctx"]
    msh, cell_tags = ctx["msh"], ctx["cell_tags"]
    facets = built["facets"]
    n_cells = built["cells"]
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
        f"fixture: t_build {t_build:.1f} s at -n {comm.size}; (A0) cells {n_cells} vs record "
        f"{F_CLINICAL_R1_CELL_RECORD}, rel {cell_rel:.3e} (band {CELL_COUNT_BAND}, ASSERTED); "
        f"(iv) failures {census_fails or 'none'}; outer census {oc}; (S) cavity max hypot "
        f"{cav_hmax:.9f} m < R_s {R_s} (ASSERTED)",
    )
    assert abs(cell_rel) <= CELL_COUNT_BAND, cell_rel
    assert not census_fails, census_fails
    assert oc["tag_collision"] == 0 and oc["cavity"] > 0, oc
    assert (
        oc["cavity_exterior"] == oc["cavity"] and oc["outer"] + oc["cavity"] == oc["exterior"]
    ), oc
    assert cav_hmax < R_s, (cav_hmax, R_s)

    # ---- (C) census, (W) widths, (Az) layout, (U1) ----------------------
    t1 = time.perf_counter()
    ids = [sp.port_id for sp in ctx["specs"]]
    geo = _sheet_geometry(ctx, comm)
    owned_counts = {pid: geo[pid][0] for pid in ids}
    c_fails = _census_failures(ids, owned_counts)
    kw = _f_clinical_kwargs("r1")
    w_expected = 2.0 * kw["ring_minor_radius"] + 2.0 * 1.0e-3
    assert abs(w_expected - EXPECTED_SHEET_WIDTH_M) < 1e-15, w_expected
    w_rel = {sp.port_id: abs(sp.sheet_width_m / w_expected - 1.0) for sp in ctx["specs"]}
    w_worst = max(w_rel.values()) if w_rel else float("inf")
    az_fails, az_worst = _layout_failures(ids, geo)
    order = [f"P{i}" for i in _ring_ports()]
    az = np.array([geo[p][1] for p in order])
    zc = np.array([geo[p][2] for p in order])
    classes, cls_res = _class_table(az, zc)
    sizes = {k: len(v) for k, v in classes.items()}
    V = ctx["solver"].function_space()
    u1 = int(V.dofmap.index_map.size_global * V.dofmap.index_map_bs)
    t_checks = time.perf_counter() - t1
    _say(
        comm,
        f"(C) {len(ids)} specs (= 2*LEG_COUNT), owned facets per sheet min "
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
        f"layout failures {az_fails or 'none'}; pair classes {len(classes)} (expected "
        f"{N_SYMMETRY_CLASSES_32}, ASSERTED), worst pair residual {cls_res:.3e} deg; sizes "
        f"equal expected {sizes == _expected_class_sizes()} (ASSERTED), total {sum(sizes.values())}",
    )
    _say(
        comm,
        f"(U1) N1curl global size x bs = {u1} vs record {F_CLINICAL_R1_UNKNOWNS_RECORD[0]} "
        f"(EXACT, ASSERTED); checks {t_checks:.1f} s",
    )
    assert not c_fails, c_fails
    assert w_worst <= STRIP_BAND, sorted(w_rel.items(), key=lambda kv: -kv[1])[:4]
    assert not az_fails, az_fails
    assert az_worst < AZIMUTH_MATCH_DEG, az_worst
    assert cls_res < AZIMUTH_MATCH_DEG, cls_res
    assert len(classes) == N_SYMMETRY_CLASSES_32, len(classes)
    assert sizes == _expected_class_sizes(), sizes
    assert u1 == F_CLINICAL_R1_UNKNOWNS_RECORD[0], u1
    del V  # the local holds the degree-1 space
    _apply_degree_c3(ctx, comm)

    if build_only:
        _m, rss_per, rss_sum = _rss_gib(comm)
        _say(
            comm,
            f"BUILD-ONLY return: t_build {t_build:.1f} s; RSS per rank "
            f"{', '.join(f'{r:.2f}' for r in (rss_per or []))} GiB, summed {rss_sum:.2f} GiB; "
            f"window {time.perf_counter() - t_window:.1f} s (no wall term, no assembly, "
            f"no factorisation)",
        )
        return {
            "t_build": t_build,
            "cells": n_cells,
            "u1": u1,
            "rss_sum": rss_sum,
            "geo": geo,
            "classes": len(classes),
        }

    return _sweep(comm, built, geo, order, az, zc, classes, t_build, t_window, n_drives)


def _sweep(comm, built, geo, order, az, zc, classes, t_build, t_window, n_drives):
    """`WF-12` step 2's Cu drive loop and matrix block (`:164–345`),
    transposed to 64 ports.  Never run by item 519 (entry 44 only)."""
    import ufl
    from dolfinx import fem

    from fem_em_solver.core.cavity import surface_resistance_ohm
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
    from tests.validation.test_port_birdcage_four_port import (
        PASSIVITY_SIGMA_TOLERANCE,
        TERMINATED_PORT_IMPEDANCE_OHM,
    )
    from tests.validation.test_port_birdcage_lumped_column import (
        ADJACENT_SPREAD_BAND,
        PHANTOM_CELL_TAG,
    )
    from tests.validation.test_port_birdcage_ring_column import (
        CONTROL_COLUMN_SCALE,
        OPPOSITE_SPREAD_BAND,
        _reciprocity_ratio,
    )
    from tests.validation.test_port_birdcage_ring_matrix import MATRIX_PASSIVITY_CEILING
    from tests.validation.test_port_lumped_sheet_sweep import RECIPROCITY_BAND
    from tests.validation.test_th14_birdcage_copper import _leontovich_term, _surface_loss_w
    from tests.validation.test_wf12_f_human_hole import _rss_gib

    ctx = built["ctx"]
    msh, tags_f, omega, cell_tags = ctx["msh"], ctx["tags_f"], ctx["omega"], ctx["cell_tags"]
    facets = built["facets"]
    z_s = (1.0 + 1.0j) * surface_resistance_ohm(omega, COPPER_SIGMA)
    wall = _leontovich_term(msh, facets, z_s, omega)
    z_p = complex(TERMINATED_PORT_IMPEDANCE_OHM)
    dx = ufl.Measure("dx", domain=msh, subdomain_data=cell_tags)
    acct = {"mesh": msh, "facet_tags": tags_f}
    solver = ctx["solver"]

    def _drive(driven_id, *, held):
        sheets = [sp.sheet(driven=(sp.port_id == driven_id)) for sp in ctx["specs"]]
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
        v_src = complex(dsheet.source_voltage_v)
        s_col = {}
        for s in sheets:
            v = complex(s.source_voltage_v) - currents[s.port_id] * z_p
            s_col[s.port_id] = (v - z_p * currents[s.port_id]) / v_src
        sigma_f = fields.sigma_field

        def _vol(measure):
            form = fem.form(0.5 * sigma_f * ufl.inner(e, e) * measure)
            return float(np.real(comm.allreduce(fem.assemble_scalar(form), op=MPI.SUM)))

        p_sheets = float(sum(_sheet_field_dissipation_w(acct, s, e, omega) for s in sheets))
        p_src = float(_source_power_w(acct, dsheet, e, omega))
        p_total = _vol(dx)
        p_phantom = _vol(dx(PHANTOM_CELL_TAG))
        p_surf = _surface_loss_w(msh, facets, z_s, e, comm)
        del fields
        t_post = time.perf_counter() - t1
        d = {
            "s": s_col,
            "kind": "held" if held else "solve",
            "t_solve": t_solve,
            "t_post": t_post,
            "p_src": p_src,
            "p_sheets": p_sheets,
            "p_phantom": p_phantom,
            "p_rest": p_total - p_phantom,
            "p_surf": p_surf,
        }
        d["rel"] = abs(p_src - p_sheets - p_total - p_surf) / abs(p_src)
        d["blind"] = abs(p_src - p_sheets - p_total) / abs(p_src)
        return d

    assert order.index(DRIVEN) == 0, order[:3]
    driven_ids = order[:n_drives]
    cols = {}
    for k, pid in enumerate(driven_ids):
        d = _drive(pid, held=(k > 0))
        cols[pid] = d
        _say(
            comm,
            f"Cu drive {k + 1:2d}/{n_drives} {pid}: {d['kind']} {d['t_solve']:.2f} s + post "
            f"{d['t_post']:.2f} s; S_driven {d['s'][pid]:.6f}; sum|S|^2 "
            f"{sum(abs(v) ** 2 for v in d['s'].values()):.6f} | P_src {d['p_src']:.9e} "
            f"P_surf {d['p_surf']:.9e} W; (P) rel {d['rel']:.3e}; NEG blind "
            f"{d['blind']:.4e} t={time.perf_counter() - t_window:.1f} s",
        )
    kinds = [cols[p]["kind"] for p in driven_ids]
    _m, rss_per, rss_sum = _rss_gib(comm)
    first = cols[DRIVEN]
    held_cad = [cols[p]["t_solve"] + cols[p]["t_post"] for p in driven_ids[1:]]
    _say(
        comm,
        f"PRICE: t_b {t_build:.1f} s; t_col {first['t_solve']:.1f} s (+ post "
        f"{first['t_post']:.1f} s); held cadence mean "
        f"{(np.mean(held_cad) if held_cad else float('nan')):.2f} s max "
        f"{(max(held_cad) if held_cad else float('nan')):.2f} s; RSS per rank "
        f"{', '.join(f'{r:.2f}' for r in (rss_per or []))} GiB, summed {rss_sum:.2f} GiB "
        f"(beside entry 40's t_b / t_col / t_held)",
    )
    worst_rel = max(cols[p]["rel"] for p in driven_ids)
    min_blind = min(cols[p]["blind"] for p in driven_ids)
    _say(
        comm,
        f"kinds {kinds.count('solve')} solve / {kinds.count('held')} held (ASSERTED 1 / "
        f"{n_drives - 1}); (P) worst rel {worst_rel:.3e} (band {DISCRETE_IDENTITY_RTOL:g}, "
        f"ASSERTED); NEG min blind {min_blind:.4e} > band (ASSERTED)",
    )
    assert kinds.count("solve") == 1 and kinds.count("held") == n_drives - 1, kinds
    assert worst_rel <= DISCRETE_IDENTITY_RTOL, worst_rel
    assert all(cols[p]["p_surf"] > 0.0 for p in driven_ids)
    assert min_blind > BLIND_FACTOR * DISCRETE_IDENTITY_RTOL, min_blind

    if n_drives < N_PORTS:
        _say(
            comm,
            f"(ii), (iii), (iv') and the column-scale control: SKIPPED "
            f"({C3_DRIVES_ENV}={n_drives} < {N_PORTS})",
        )
        return {"cols": cols, "kinds": kinds, "t_build": t_build}

    s = np.array([[cols[pj]["s"][pi] for pj in order] for pi in order], dtype=complex)
    ratio = _reciprocity_ratio(s)
    control = s.copy()
    control[:, 0] *= CONTROL_COLUMN_SCALE
    control_ratio = _reciprocity_ratio(control)
    sigma_max = float(np.linalg.norm(s, 2))
    col_norms = np.sum(np.abs(s) ** 2, axis=0)
    spreads = _class_spreads(s, classes)
    worst = max(spreads, key=lambda key: spreads[key][0])
    _say(
        comm,
        f"(ii) reciprocity {ratio:.6e} (band {RECIPROCITY_BAND:.0e}, ASSERTED; PREDICTED "
        f"≈{PREDICTED_RECIPROCITY:.0e})",
    )
    _say(
        comm,
        f"(iii) sigma_max {sigma_max:.12f} (ceiling {MATRIX_PASSIVITY_CEILING}, ASSERTED; "
        f"entry 44 names 1 + {PASSIVITY_SIGMA_TOLERANCE:g}, printed); column sum|S|^2 min "
        f"{col_norms.min():.6f} max {col_norms.max():.6f} (ASSERTED <= ceiling)",
    )
    for same_ring in (True, False):
        for steps in range(N_AZIMUTH_CLASSES_32):
            sp, cnt, mean = spreads[(same_ring, steps)]
            _say(
                comm,
                f"     {'same ' if same_ring else 'other'} ring, {steps:2d} steps  n = {cnt:4d}"
                f"  mean |S| = {mean:.9e}  spread {sp * 100:8.4f}%"
                + ("   <-- worst" if (same_ring, steps) == worst else ""),
            )
    _say(
        comm,
        f"(iv') worst class spread {spreads[worst][0]:.4e} (band OPPOSITE_SPREAD_BAND "
        f"{OPPOSITE_SPREAD_BAND}, ASSERTED; vs ADJACENT_SPREAD_BAND {ADJACENT_SPREAD_BAND} "
        f"= {spreads[worst][0] / ADJACENT_SPREAD_BAND:.2f}x, PRINTED; F-human "
        f"{PREDICTED_SPREAD_F_HUMAN})",
    )
    _say(
        comm,
        f"control (predicted, printed): column {order[0]} x {CONTROL_COLUMN_SCALE}: ratio "
        f"{control_ratio:.6e} (ceiling ≈{N2_CEILING:.2e}; band {RECIPROCITY_BAND:.0e})",
    )
    save_path = _save_path()
    if comm.rank == 0:
        save_path.parent.mkdir(parents=True, exist_ok=True)
        np.savez(
            save_path,
            s64=s.astype(np.complex128),
            order=np.array(order),
            azimuth_deg=az,
            z_m=zc,
            cells=built["cells"],
            n_ranks=comm.size,
        )
    comm.Barrier()
    _say(comm, f"saved {save_path}; window {time.perf_counter() - t_window:.1f} s")
    assert ratio <= RECIPROCITY_BAND, ratio
    assert sigma_max <= MATRIX_PASSIVITY_CEILING, sigma_max
    assert float(col_norms.max()) <= MATRIX_PASSIVITY_CEILING, col_norms.max()
    assert spreads[worst][0] <= OPPOSITE_SPREAD_BAND, spreads[worst]
    return {
        "cols": cols,
        "kinds": kinds,
        "t_build": t_build,
        "s": s,
        "ratio": ratio,
        "sigma_max": sigma_max,
        "spreads": spreads,
    }


# ---------------------------------------------------------------------------
# Coding tests (numpy only); skip when FEM_EM_GEO38_C3 is set.
# ---------------------------------------------------------------------------
_coding = pytest.mark.skipif(
    os.environ.get(C3_ENV, "").strip() != "",
    reason=f"{C3_ENV} set: coding tests skip in solve windows",
)


@_coding
def test_k0_knob_parsing(monkeypatch):
    monkeypatch.delenv(C3_ENV, raising=False)
    assert _mode() is None
    monkeypatch.setenv(C3_ENV, "build")
    assert _mode() == "build"
    monkeypatch.setenv(C3_ENV, "1")
    assert _mode() == "1"
    for bad in ("2", "x"):
        monkeypatch.setenv(C3_ENV, bad)
        with pytest.raises(ValueError):
            _mode()
    monkeypatch.delenv(C3_DRIVES_ENV, raising=False)
    assert _n_drives() == 64
    for good in ("1", "64", "7"):
        monkeypatch.setenv(C3_DRIVES_ENV, good)
        assert _n_drives() == int(good)
    for bad in ("0", "65"):
        monkeypatch.setenv(C3_DRIVES_ENV, bad)
        with pytest.raises(ValueError):
            _n_drives()


@_coding
def test_k1_ideal_layout_has_34_classes():
    az, zc = _ideal_layout()
    classes, res = _class_table(az, zc)
    sizes = {k: len(v) for k, v in classes.items()}
    print(f"[GEO-38 C3 K1] classes {len(classes)}, residual {res:.3e}, total {sum(sizes.values())}")
    assert N_SYMMETRY_CLASSES_32 == 34 and N_AZIMUTH_CLASSES_32 == 17
    assert abs(AZIMUTH_STEP_32_DEG - 11.25) < 1e-15
    assert len(classes) == 34
    assert res <= 1e-12, res
    assert sizes == _expected_class_sizes()
    assert sum(sizes.values()) == 4096
    geo = {f"P{33 + k}": (1, float(az[k]), float(zc[k])) for k in range(64)}
    fails, worst = _layout_failures(list(geo), geo)
    assert not fails and worst <= 1e-12, (fails, worst)


@_coding
def test_k2_synthetic_matrix_gates():
    from tests.validation.test_port_birdcage_ring_column import _reciprocity_ratio

    az, zc = _ideal_layout()
    s = _synthetic_s(az, zc)
    classes, _res = _class_table(az, zc)
    spreads = _class_spreads(s, classes)
    ratio = _reciprocity_ratio(s)
    sigma_max = float(np.linalg.norm(s, 2))
    print(
        f"[GEO-38 C3 K2] reciprocity {ratio:.3e}, worst spread "
        f"{max(v[0] for v in spreads.values()):.3e}, sigma_max {sigma_max:.6f}"
    )
    assert ratio <= 1e-15, ratio
    assert all(v[0] == 0.0 for v in spreads.values()), spreads


@_coding
def test_n1_f_human_classifier_misses_odd_separations():
    """N1: F-human's 22.5° `_azimuth_class` on the 32-leg layout leaves an
    11.25° residual on odd separations (the import trap); ceiling exactly 11.25°."""
    from tests.validation.test_port_birdcage_ring_column import AZIMUTH_MATCH_DEG
    from tests.validation.test_port_birdcage_ring_matrix import _azimuth_class

    az, _zc = _ideal_layout()
    worst = max(_azimuth_class(az[i], az[j])[1] for i in range(32) for j in range(32))
    odd = max(_azimuth_class(az[0], az[j])[1] for j in range(1, 32, 2))
    print(f"[GEO-38 C3 N1] F-human classifier worst residual {worst:.6f} deg (odd {odd:.6f})")
    assert abs(worst - 11.25) <= 1e-9 and abs(odd - 11.25) <= 1e-9, (worst, odd)
    assert worst > AZIMUTH_MATCH_DEG


@_coding
def test_n2_column_scale_lifts_reciprocity_above_band():
    from tests.validation.test_port_birdcage_ring_column import (
        CONTROL_COLUMN_SCALE,
        _reciprocity_ratio,
    )
    from tests.validation.test_port_lumped_sheet_sweep import RECIPROCITY_BAND

    az, zc = _ideal_layout()
    s = _synthetic_s(az, zc)
    s[:, 0] *= CONTROL_COLUMN_SCALE
    ratio = _reciprocity_ratio(s)
    print(
        f"[GEO-38 C3 N2] controlled ratio {ratio:.4e} = {ratio / RECIPROCITY_BAND:.3f}x band "
        f"(ceiling ≈{N2_CEILING:.3e})"
    )
    assert ratio > RECIPROCITY_BAND, ratio
    assert ratio <= N2_CEILING * (1.0 + 1e-9), ratio


@_coding
def test_n3_census_with_p96_removed_fails():
    from tests.validation.test_geo38_f_clinical_build import _ring_ports

    ids = [f"P{i}" for i in _ring_ports()]
    counts = {p: 10 for p in ids}
    assert _census_failures(ids, counts) == []
    short = [p for p in ids if p != "P96"]
    fails = _census_failures(short, counts)
    print(f"[GEO-38 C3 N3] 63 specs -> {fails}")
    assert len(short) == 63 and fails
    counts0 = dict(counts, P50=0)
    assert _census_failures(ids, counts0)


@_coding
def test_k0b_degree_knob_parsing(monkeypatch):
    monkeypatch.delenv(C3_DEGREE_ENV, raising=False)
    assert _degree_knob_c3() is None
    monkeypatch.setenv(C3_DEGREE_ENV, "")
    assert _degree_knob_c3() is None
    monkeypatch.setenv(C3_DEGREE_ENV, "2")
    assert _degree_knob_c3() == 2
    for bad in ("1", "x"):
        monkeypatch.setenv(C3_DEGREE_ENV, bad)
        with pytest.raises(ValueError):
            _degree_knob_c3()


@_coding
def test_k4_save_path_route(monkeypatch):
    """(K4): the unset knob reads the module global at call time (the module
    the body reads, past pytest's double import); C11's patch is honoured."""
    import sys

    import tests.validation.test_geo38_f_clinical_3t_matrix as c11

    this = sys.modules[__name__]
    sentinel = Path("logs/geo38/k4-sentinel.npz")
    monkeypatch.delenv(C3_DEGREE_ENV, raising=False)
    original = this.SAVE_PATH
    assert this._save_path() == original
    monkeypatch.setattr(this, "SAVE_PATH", sentinel)
    assert this._save_path() == sentinel
    monkeypatch.setenv(C3_DEGREE_ENV, "2")
    assert this._save_path() == SAVE_PATH_D2
    assert SAVE_PATH_D2 not in (original, sentinel, c11.C11_SAVE_PATH)
    monkeypatch.undo()

    c3 = c11._modules()[0]
    with pytest.MonkeyPatch.context() as mp:
        mp.delenv(C3_DEGREE_ENV, raising=False)
        c11._c11_patches(mp)
        assert c3._save_path() == c11.C11_SAVE_PATH
    print(f"\n[GEO-38 C3 K4] d1 {original}, d2 {SAVE_PATH_D2}, C11 {c11.C11_SAVE_PATH}", flush=True)


@_coding
def test_n4_deg_refuses_degree_1_count():
    """N4 (asserted, integers only): (Deg) against 2E + 2F refuses the degree-1
    size; the separation is U2 / U1 = 5.17×."""
    from tests.validation.test_geo38_f_clinical_build import F_CLINICAL_R1_UNKNOWNS_RECORD

    u1, u2 = F_CLINICAL_R1_UNKNOWNS_RECORD
    n_edges, n_faces = 1_525_193, 2_417_644
    assert n_edges == u1
    assert _check_unknowns(u2, 2, n_edges, n_faces) == u2 == 7_885_674
    assert _check_unknowns(u1, 1, n_edges, n_faces) == u1
    with pytest.raises(AssertionError):
        _check_unknowns(u1, 2, n_edges, n_faces)
    print(f"\n[GEO-38 C3 N4] U2 / U1 = {u2 / u1:.4f} (> 5, ASSERTED)", flush=True)
    assert u2 / u1 > 5.0, (u1, u2)


@_coding
def test_k5_build_keeps_save_state(tmp_path):
    """(S) and its control (§9 item 706, C17), on a file present and on one
    absent: a build that leaves the path alone passes; one that touches or
    creates it fails on every rank.  All file work is rank 0's, on rank 0's
    ``tmp_path``, so the other rank's verdict can only come from the bcast."""
    comm = MPI.COMM_WORLD
    root = Path(comm.bcast(str(tmp_path), root=0))
    present, absent = root / "present.npz", root / "absent.npz"
    if comm.rank == 0:
        present.write_bytes(b"s64")
    comm.Barrier()

    def nothing():
        return "built"

    def touch():
        if comm.rank == 0:
            st = present.stat()
            os.utime(present, ns=(st.st_atime_ns, st.st_mtime_ns + 1_000_000_000))

    def create():
        if comm.rank == 0:
            absent.write_bytes(b"s64")

    out, state = _build_keeps_save_state(comm, present, nothing)
    assert out == "built"
    assert state == (3, comm.bcast(present.stat().st_mtime_ns if comm.rank == 0 else None, root=0))
    out, state = _build_keeps_save_state(comm, absent, nothing)
    assert out == "built"
    assert state is None
    with pytest.raises(AssertionError):
        _build_keeps_save_state(comm, present, touch)
    with pytest.raises(AssertionError):
        _build_keeps_save_state(comm, absent, create)
    print(
        "\n[GEO-38 C3 K5] (S) identical present and absent; touch and create both fail it",
        flush=True,
    )


@_coding
def test_k6_rung_knob_parsing(monkeypatch):
    """(K6): the rung knob, and r2 refused at degree 2 before gmsh."""
    monkeypatch.delenv(C3_RUNG_ENV, raising=False)
    monkeypatch.delenv(C3_DEGREE_ENV, raising=False)
    assert _rung_knob_c3() == "r1"
    for raw, want in (("", "r1"), ("r1", "r1"), ("r2", "r2")):
        monkeypatch.setenv(C3_RUNG_ENV, raw)
        assert _rung_knob_c3() == want and _rung_checked() == want, raw
    for bad in ("r3", "R2", "2"):
        monkeypatch.setenv(C3_RUNG_ENV, bad)
        with pytest.raises(ValueError):
            _rung_knob_c3()
        with pytest.raises(ValueError):
            _rung_checked()
    monkeypatch.setenv(C3_RUNG_ENV, "r2")
    monkeypatch.setenv(C3_DEGREE_ENV, "2")
    with pytest.raises(ValueError, match="degree 1 only"):
        _rung_checked()
    monkeypatch.setenv(C3_RUNG_ENV, "r1")
    assert _rung_checked() == "r1"
    print(
        "\n[GEO-38 C3 K6] rung knob r1 / r2 accepted, r3 / R2 / 2 refused; r2 at degree 2 refused",
        flush=True,
    )


@_coding
def test_k7_r2_patches_route_and_restore():
    """(K7): under `_r2_patches` only the rung moves at the builder; the records
    and the save path are r2's; all of it is restored on exit."""
    import sys

    import tests.validation.test_geo38_f_clinical_3t_matrix as c11
    import tests.validation.test_geo38_f_clinical_build as s1a

    this = sys.modules[__name__]
    r1_kw = s1a._f_clinical_kwargs("r1")
    with pytest.MonkeyPatch.context() as mp:
        mp.delenv(C3_DEGREE_ENV, raising=False)
        _r2_patches(mp)
        seen = c11._tripwire_kwargs(mp)
        assert seen["conductor_resolution"] == s1a.CONDUCTOR_RESOLUTION["r2"] == 2.4e-3, seen
        assert set(seen) - {"comm", "return_diagnostics"} == set(r1_kw), sorted(
            set(seen) ^ set(r1_kw)
        )
        moved = sorted(k for k in r1_kw if seen[k] != r1_kw[k])
        assert moved == ["conductor_resolution"], moved
        assert s1a.F_CLINICAL_R1_CELL_RECORD == 912_804
        assert s1a.F_CLINICAL_R1_UNKNOWNS_RECORD == (1_184_518, 6_179_214)
        assert this._save_path() == SAVE_PATH_R2
        # the bare name SAVE_PATH is rebound here too (this module's globals); compare the literal d1
        assert SAVE_PATH_R2 not in (ENTRY44_D1_PATH, SAVE_PATH_D2, c11.C11_SAVE_PATH)
    with pytest.MonkeyPatch.context() as mp:
        restored = c11._tripwire_kwargs(mp)
    assert restored["conductor_resolution"] == 1.6e-3, restored
    assert s1a.F_CLINICAL_R1_CELL_RECORD == 1_147_366
    assert s1a.F_CLINICAL_R1_UNKNOWNS_RECORD == (1_525_193, 7_885_674)
    assert this.SAVE_PATH == ENTRY44_D1_PATH
    print(
        f"\n[GEO-38 C3 K7] r2 patches: builder sees only conductor_resolution "
        f"{seen['conductor_resolution']} (r1 {r1_kw['conductor_resolution']}), save "
        f"{SAVE_PATH_R2}; restored on exit; c11's c3 is this module: {c11._modules()[0] is this}",
        flush=True,
    )


@_coding
def test_n5_r1_records_refuse_r2_mesh():
    """N5 (asserted, integers only): r1's records refuse the r2 mesh."""
    import tests.validation.test_geo38_f_clinical_build as s1a
    from tests.mesh.test_birdcage_port_sheet_prerequisite import CELL_COUNT_BAND

    cells_rel = abs(s1a.F_CLINICAL_R2_CELL_RECORD / s1a.F_CLINICAL_R1_CELL_RECORD - 1.0)
    u1_rel = abs(s1a.F_CLINICAL_R2_UNKNOWNS_RECORD[0] / s1a.F_CLINICAL_R1_UNKNOWNS_RECORD[0] - 1.0)
    print(
        f"\n[GEO-38 C3 N5] |912804 / 1147366 - 1| = {cells_rel:.4f} = "
        f"{cells_rel / CELL_COUNT_BAND:.1f}x CELL_COUNT_BAND; |1184518 / 1525193 - 1| = "
        f"{u1_rel:.4f}",
        flush=True,
    )
    assert cells_rel > CELL_COUNT_BAND, cells_rel
    assert s1a.F_CLINICAL_R2_UNKNOWNS_RECORD[0] != s1a.F_CLINICAL_R1_UNKNOWNS_RECORD[0]
    assert u1_rel > 0.2, u1_rel


@_coding
def test_n6_rung_knob_does_not_leak(monkeypatch):
    """N6: `_run_c3`'s route never reads the rung knob, and a C11 window under
    the knob sees exactly the kwargs it sees with it unset."""
    import tests.validation.test_geo38_f_clinical_3t_matrix as c11

    monkeypatch.delenv(C3_DEGREE_ENV, raising=False)
    monkeypatch.setenv(C3_RUNG_ENV, "r2")
    with pytest.MonkeyPatch.context() as mp:
        bare = c11._tripwire_kwargs(mp)
    assert bare["conductor_resolution"] == 1.6e-3, bare
    with pytest.MonkeyPatch.context() as mp:
        c11._c11_patches(mp)
        c11_knob = c11._tripwire_kwargs(mp)
    monkeypatch.delenv(C3_RUNG_ENV)
    with pytest.MonkeyPatch.context() as mp:
        c11._c11_patches(mp)
        c11_unset = c11._tripwire_kwargs(mp)
    assert c11_knob == c11_unset, sorted(k for k in c11_unset if c11_knob.get(k) != c11_unset[k])
    print(
        f"\n[GEO-38 C3 N6] knob r2, no patches: conductor_resolution {bare['conductor_resolution']}; "
        f"C11 kwargs with knob == unset: {c11_knob == c11_unset}",
        flush=True,
    )


# ---------------------------------------------------------------------------
# The gated test.
# ---------------------------------------------------------------------------
def _size_mtime(path):
    if not path.exists():
        return None
    st = path.stat()
    return (st.st_size, st.st_mtime_ns)


def _build_keeps_save_state(comm, path, build):
    """(S): ``path``'s (size, mtime_ns), or its absence, is identical before and
    after ``build()``.  Rank 0's stat, bcast, so two ranks' calls are never
    compared.  Returns ``build()``'s value and the state."""
    before = comm.bcast(_size_mtime(path) if comm.rank == 0 else None, root=0)
    out = build()
    after = comm.bcast(_size_mtime(path) if comm.rank == 0 else None, root=0)
    assert after == before, (path, after, before)
    return out, before


@complex_only
def test_geo38_c3_f_clinical_copper_matrix(monkeypatch):
    rung = _rung_checked()
    mode = _mode()
    if mode is None:
        pytest.skip(f"{C3_ENV} unset: GEO-38 C3 is selected by environment only")
    comm = MPI.COMM_WORLD
    literal = (ENTRY44_D1_PATH, ENTRY54_D2_PATH)
    if rung == "r2":
        lit_before = comm.bcast(
            [_size_mtime(p) for p in literal] if comm.rank == 0 else None, root=0
        )
        lo, hi = PREDICTED_T_BUILD_R2_S[1]
        _say(
            comm,
            f"(C35) rung r2: save path {SAVE_PATH_R2}; PREDICTED t_build ≈ "
            f"{PREDICTED_T_BUILD_R2_S[0]:.0f} s in [{lo:.0f}, {hi:.0f}] (printed only)",
        )
        _r2_patches(monkeypatch)
    if mode == "build":
        save_path = _save_path()
        d1_before = _size_mtime(SAVE_PATH)
        _, save_state = _build_keeps_save_state(
            comm, save_path, lambda: _run_c3(comm, build_only=True)
        )
        d1_after = _size_mtime(SAVE_PATH)
        _say(
            comm,
            f"(S) save path {save_path} (size, mtime_ns) or absence unchanged by build-only: "
            f"{save_state} (ASSERTED)",
        )
        _say(comm, f"(S) {SAVE_PATH} (size, mtime_ns) {d1_after} == before {d1_before} (ASSERTED)")
        assert d1_after == d1_before, (d1_after, d1_before)
        if rung == "r2":
            lit_after = comm.bcast(
                [_size_mtime(p) for p in literal] if comm.rank == 0 else None, root=0
            )
            _say(
                comm,
                f"(S-r2) {save_path} absent before and after: {save_state is None}; "
                f"{ENTRY44_D1_PATH} / {ENTRY54_D2_PATH} (size, mtime_ns) {lit_after} == before "
                f"{lit_before} (ASSERTED)",
            )
            assert save_path == SAVE_PATH_R2 and save_state is None, (save_path, save_state)
            assert lit_after == lit_before, (lit_after, lit_before)
        return
    _run_c3(comm, build_only=False, n_drives=_n_drives())
