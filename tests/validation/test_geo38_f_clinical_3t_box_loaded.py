"""GEO-38 C21 (§9 item 775, §10 chain C): the 3 T box-loaded body — C6a's
module (``test_geo38_f_clinical_box_loaded.py``) composed with C20's mesh route
and C12's frequency patches, stub-run (no solve, no gmsh).

Brief: ``docs/planning/chunks/GEO-38.md`` "§9 item 775".

Composition, never an edit: C6a's `_run_build` / `_run_box` run unchanged under
`_c21_patches`, which sets

1. ``s2._f_clinical_kwargs`` -> ``t3._f_clinical_3t_kwargs`` (C20's route,
   `test_geo38_f_clinical_3t_rod_box.py:106–110`).  ``_box_kwargs`` is a
   name-bound import in I35 (`test_geo36_rod_f_clinical.py:55`) and
   ``_f_clinical_kwargs`` one in the box module, so C11's ``s1a`` patch never
   reaches the box build — K7-3T's negative control shows exactly that.  The rod
   lands at I35's `_rod_centre()`, which on the 3 T table is C20's `_kwargs_3t`
   centre (K7-3T asserts the whole kwargs against `_kwargs_3t`).
2. ``i35.F_CLINICAL_ROD_CELL_RECORD`` / ``_UNKNOWNS_RECORD`` -> C20's
   1 128 430 / (1 508 256, 7 783 028), read by `i36._route_drive` at call time
   into ``s1a``'s records, so (A0) and (U1) assert C20's values.
3. ``wf7.FREQUENCY_HZ`` and ``c4.FREQUENCY_HZ_C4`` -> ``t3.FREQUENCY_3T_HZ``
   (C12's `_c12_patches`, `test_geo38_f_clinical_3t_esr_root.py:139`).
4. ``box._box_readouts`` -> a wrapper that calls it unchanged and prints the
   lattice mean |B₁⁺| beside the centroid mean (both labelled; neither asserted;
   which is the record is the 11-04 weekly's call).  The lattice values are
   already returned by `_box_readouts` (``lattice_b_plus``), so C6a's module is
   not edited.

``c4._feeds_32`` is unchanged (P33 / P41 at 5.625° / 95.625° at both field
strengths).  The JSON is C6a's ``c6-node-<pF>.json``; the node differs from
entry 52's, so the name cannot collide.

Selected by ``FEM_EM_GEO38_C21`` only (never by ``-k``): unset ⇒ the gated test
skips and the coding tests run; ``build`` ⇒ C6a's `_run_build` under the patch
set; ``1`` ⇒ C6a's `_run_box` under the patch set (node knob
``FEM_EM_GEO36_S8_DRIVE_C_PF`` required, refused under ``build``); anything
else ⇒ ``ValueError``.

Scope: the module and its stub-run coding tests.  Not claimed: that the 3 T
build reaches C20's cells (item 791), any field, any solve.

C52 (§9 item 929): ``FEM_EM_GEO38_C21_DEGREE`` — unset / empty ⇒ degree 1,
``2`` ⇒ 2, anything else ⇒ ``ValueError`` — read by the gated test only.  Under
``1`` it is routed to C6a's `_run_box(degree=…)`; under ``build`` at ``2`` the
route is I36a's (`i36._route_drive` after `_c21_patches`, then C4a's
`_run_build`, which builds both degrees and asserts (Deg-2) against C20's
record), the first degree-2 space on C20's mesh.  Unlike C6a's knob it is
accepted under ``build``.  Not claimed: any 3 T degree-2 solve.
"""

from __future__ import annotations

import os

import numpy as np
import pytest
from mpi4py import MPI

from tests.complex_mode import complex_only

C21_ENV = "FEM_EM_GEO38_C21"
#: C52 (§9 item 929): the 3 T degree knob (C34's pattern, accepted under ``build``).
C21_DEGREE_ENV = "FEM_EM_GEO38_C21_DEGREE"
#: The 3 T node (``xl`` entry 49 (C16): C_in(F-clinical-3T, d1) = 44.0615 pF).
K_NODE_3T_PF = 44.0615
#: Variables a C21 window must never read: C6a's list plus C6's own knob and the
#: C11 / C12 / C20 knobs.
FOREIGN_ENVS_C21_EXTRA = (
    "FEM_EM_GEO38_C6",
    "FEM_EM_GEO38_C11",
    "FEM_EM_GEO38_C11_DRIVES",
    "FEM_EM_GEO38_C12",
    "FEM_EM_GEO38_C12_NODES_PF",
    "FEM_EM_GEO38_C12_DEGREE",
    "FEM_EM_GEO38_C20",
)

#: Filled by the lattice-mean wrapper: ``{sense: (centroid_mean, lattice_mean, n)}``.
_LATTICE_MEANS: list = []


def _modules():
    import tests.validation.test_geo36_rod_f_clinical as i35
    import tests.validation.test_geo36_rod_f_clinical_drive as i36
    import tests.validation.test_geo38_f_clinical_3t as t3
    import tests.validation.test_geo38_f_clinical_3t_rod_box as c20
    import tests.validation.test_geo38_f_clinical_box as s2
    import tests.validation.test_geo38_f_clinical_box_loaded as box
    import tests.validation.test_geo38_f_clinical_build as s1a
    import tests.validation.test_geo38_f_clinical_esr_root as c4
    import tests.validation.test_wf7_f_human_b1_quadrature as wf7

    return {"i35": i35, "i36": i36, "t3": t3, "c20": c20, "s2": s2, "box": box, "s1a": s1a,
            "c4": c4, "wf7": wf7}


def _foreign_envs():
    return tuple(_modules()["box"].FOREIGN_ENVS) + FOREIGN_ENVS_C21_EXTRA


def _say(comm, msg):
    if comm.rank == 0:
        print(f"[GEO-38 C21] {msg}", flush=True)


# ---------------------------------------------------------------------------
# Knobs.
# ---------------------------------------------------------------------------
def _c21_mode():
    raw = os.environ.get(C21_ENV, "").strip()
    if raw == "":
        return None
    if raw in ("build", "1"):
        return raw
    raise ValueError(f"{C21_ENV} must be unset, 'build' or '1', got {raw!r}")


def _knobs():
    """Foreign guard, then ``(mode, c_pf)`` with I36b's node-knob semantics."""
    i36 = _modules()["i36"]
    bad = [k for k in _foreign_envs() if os.environ.get(k) is not None]
    if bad:
        raise ValueError(f"C21 never reads these variables; unset them: {bad}")
    mode = _c21_mode()
    return mode, i36._node_knob(mode)


def _c21_degree(mode):
    """C52: `_c6_degree`'s parsing under ``C21_DEGREE_ENV`` — unset / empty ⇒
    ``None`` (degree 1), ``"2"`` ⇒ 2, anything else ⇒ ValueError — **accepted under
    ``build``** (no build has made the degree-2 space on C20's mesh; item 929
    premise 3).  Separate from `_knobs`, so K1's equalities are unchanged."""
    raw = os.environ.get(C21_DEGREE_ENV, "").strip()
    if raw == "":
        return None
    if raw == "2":
        return 2
    raise ValueError(f"{C21_DEGREE_ENV}={raw!r}: only unset/empty (degree 1) or '2' are accepted")


# ---------------------------------------------------------------------------
# The 3 T patch set.
# ---------------------------------------------------------------------------
def _c21_patches(monkeypatch):
    m = _modules()
    s2, t3, c20, i35, wf7, c4, box = (m[k] for k in ("s2", "t3", "c20", "i35", "wf7", "c4", "box"))
    _LATTICE_MEANS.clear()
    # 1. mesh: C20's route (the box module's name-bound `_f_clinical_kwargs`).
    monkeypatch.setattr(s2, "_f_clinical_kwargs", t3._f_clinical_3t_kwargs)
    # 2. records: read by `i36._route_drive` at call time.
    monkeypatch.setattr(i35, "F_CLINICAL_ROD_CELL_RECORD", c20.F_CLINICAL_3T_ROD_CELL_RECORD)
    monkeypatch.setattr(i35, "F_CLINICAL_ROD_UNKNOWNS_RECORD", c20.F_CLINICAL_3T_ROD_UNKNOWNS_RECORD)
    # 3. frequency: C12's precedent.
    monkeypatch.setattr(wf7, "FREQUENCY_HZ", t3.FREQUENCY_3T_HZ)
    monkeypatch.setattr(c4, "FREQUENCY_HZ_C4", t3.FREQUENCY_3T_HZ)
    # 4. the lattice mean beside the centroid mean (print-only wrapper).
    orig = box._box_readouts

    def _readouts_with_lattice_mean(*a, **k):
        out = orig(*a, **k)
        comm = a[9] if len(a) > 9 else k["comm"]
        for sense in ("ccw", "cw"):
            s = out["senses"][sense]
            lat = np.asarray(s["lattice_b_plus"], dtype=float)
            _LATTICE_MEANS.append((sense, float(s["mean_b_plus"]), float(np.mean(lat)), int(lat.size)))
            box._say(comm, f"(C21) mean |B1+|_{sense}: centroid mean {s['mean_b_plus']:.9e} T over "
                           f"{s['n_valid']} tag-3 centroids; lattice mean {np.mean(lat):.9e} T over the "
                           f"{lat.size}-point lattice (both printed, neither asserted; the record is "
                           f"the 11-04 weekly's call)")
        return out

    monkeypatch.setattr(box, "_box_readouts", _readouts_with_lattice_mean)


@complex_only
def test_geo38_c21_box_loaded_3t(monkeypatch):
    mode, c_pf = _knobs()  # node knob missing under 1 / set under build => ValueError, before gmsh
    if mode is None:
        pytest.skip(f"{C21_ENV} unset: GEO-38 C21 is selected by environment only")
    deg = _c21_degree(mode)  # C52: before gmsh
    m = _modules()
    comm = MPI.COMM_WORLD
    _c21_patches(monkeypatch)
    _say(comm, f"3 T patch set: rung Ø {m['t3'].RUNG_DIAMETER_3T_M} m, f {m['wf7'].FREQUENCY_HZ:.6e} Hz "
               f"(C4 {m['c4'].FREQUENCY_HZ_C4:.6e}); records cells "
               f"{m['i35'].F_CLINICAL_ROD_CELL_RECORD}, U {m['i35'].F_CLINICAL_ROD_UNKNOWNS_RECORD}; "
               f"degree {deg or 1} ({C21_DEGREE_ENV} {deg if deg else 'unset'})")
    if mode == "1":
        m["box"]._run_box(comm, monkeypatch, c_pf, degree=deg)
        return
    # (Bore), §9 item 791: C20's `_bore_gap` on the kwargs the patched route builds
    # (K7-3T asserts they are C20's `_kwargs_3t`); the box sits strictly inside the rungs.
    gap, corner = m["c20"]._bore_gap(m["i35"]._kwargs())
    _say(comm, f"(Bore) corner {corner:.9f} m vs rung radius {m['t3'].RUNG_DIAMETER_3T_M / 2.0} m, "
               f"gap {gap:.9f} m (> 0, ASSERTED)")
    assert gap > 0.0, (gap, corner)
    if deg is None:
        m["box"]._run_build(comm, monkeypatch)
        return
    # C52, degree 2: I36a's build route (`test_geo36_rod_f_clinical_drive.py:656–659`).
    # `_route_drive` after `_c21_patches` (it copies C20's records into s1a at call time).
    records = m["i36"]._route_drive(monkeypatch)
    m["c4"]._run_build(comm, monkeypatch)  # (A0) / (U1) / (Z) / (Deg-2) + NEG asserted inside
    m2 = m["i36"]._m2_failures(records)
    _say(comm, f"(M2) degrees {[d for d, _m in records]}; failures {m2 or 'none'} (ASSERTED)")
    assert not m2, m2


# ---------------------------------------------------------------------------
# Coding tests (no gmsh); skip whenever C21_ENV is set.
# ---------------------------------------------------------------------------
_coding = pytest.mark.skipif(os.environ.get(C21_ENV, "").strip() != "",
                             reason=f"{C21_ENV} set: coding tests skip in windows")


@_coding
def test_c21_k1_knob_parsing(monkeypatch):
    from tests.validation.test_geo36_rod_f_clinical_drive import S8_DRIVE_C_PF_ENV

    foreign = _foreign_envs()
    for k in (C21_ENV, S8_DRIVE_C_PF_ENV) + foreign:
        monkeypatch.delenv(k, raising=False)
    assert _knobs() == (None, None)
    monkeypatch.setenv(C21_ENV, "build")
    assert _knobs() == ("build", None)
    monkeypatch.setenv(S8_DRIVE_C_PF_ENV, "44.0615")
    with pytest.raises(ValueError, match="refused"):
        _knobs()
    monkeypatch.setenv(C21_ENV, "1")
    assert _knobs() == ("1", K_NODE_3T_PF)
    for bad in ("2", "x"):
        monkeypatch.setenv(C21_ENV, bad)
        with pytest.raises(ValueError):
            _c21_mode()
    monkeypatch.setenv(C21_ENV, "build")
    monkeypatch.delenv(S8_DRIVE_C_PF_ENV)
    for k in foreign:
        with monkeypatch.context() as m:
            m.setenv(k, "1")
            with pytest.raises(ValueError):
                _knobs()
    assert set(FOREIGN_ENVS_C21_EXTRA) <= set(foreign) and "FEM_EM_GEO38_C6" in foreign
    print(f"\n[GEO-38 C21 K1] mode unset/build/1 accepted, 2/x refused; node knob refused under "
          f"build, read under 1 ({K_NODE_3T_PF}); {len(foreign)} foreign variables refused "
          f"(C6a's {len(foreign) - len(FOREIGN_ENVS_C21_EXTRA)} + {FOREIGN_ENVS_C21_EXTRA}) (ASSERTED)")


@_coding
def test_c21_k2_one_without_node_stops_before_gmsh(monkeypatch):
    from fem_em_solver.io.mesh import MeshGenerator
    from tests.validation.test_geo36_rod_f_clinical_drive import S8_DRIVE_C_PF_ENV

    def _sentinel(*a, **k):
        raise AssertionError("gmsh reached")

    with monkeypatch.context() as m:
        for k in _foreign_envs() + (S8_DRIVE_C_PF_ENV,):
            m.delenv(k, raising=False)
        m.setenv(C21_ENV, "1")
        m.setattr(MeshGenerator, "birdcage_port_domain", staticmethod(_sentinel))
        with pytest.raises(ValueError, match=S8_DRIVE_C_PF_ENV):
            test_geo38_c21_box_loaded_3t(m)
    print("\n[GEO-38 C21 K2] C21=1 without the node knob raised ValueError before gmsh (ASSERTED)")


@_coding
def test_c21_k5_sense_weights():
    """(K5): ccw = (1, −1j), cw = (1, +1j)."""
    from tests.validation.test_wf10_f_human_tuned_b1_sar import _sense_weights

    w = _sense_weights()
    d = [max(abs(np.asarray(w[s]) - np.array(t))) for s, t in (("ccw", (1.0, -1.0j)), ("cw", (1.0, 1.0j)))]
    print(f"\n[GEO-38 C21 K5] ccw {w['ccw']}, cw {w['cw']}; misses {d} (<= 1e-15, ASSERTED)")
    assert max(d) <= 1.0e-15, d


@complex_only
@_coding
def test_c21_k6_readouts_at_128mhz():
    """(K6-3T): C6a's K6 at ω = 2π·128e6 (C6a hard-codes 64 MHz at `:769`; not
    edited).  (K6-B) |B₁^±| = |c|/(2ω), cw twice; cw/ccw SAR at c* = 4; (K6-i)
    coverage power and mass; (K6-c) tag 403 zeroed on one side misses."""
    from fem_em_solver.post.sar import build_density_field
    from tests.validation.test_birdcage_power_identity import DISCRETE_IDENTITY_RTOL
    from tests.validation.test_birdcage_sar_mass_averaged import (
        EXACT_IDENTITY_RTOL,
        PHANTOM_RHO_KG_PER_M3,
    )
    from tests.validation.test_wf7_f_human_sar10g_hotspot import COVERAGE_BALL_PAD_M

    m = _modules()
    box, t3 = m["box"], m["t3"]
    comm = MPI.COMM_WORLD
    msh, cell_tags, e_ccw, e_cw, sigma_f, c, half = box._k6_fixture(comm)
    omega = 2.0 * np.pi * t3.FREQUENCY_3T_HZ
    centre = np.zeros(3)
    radius = float(np.linalg.norm(half) + COVERAGE_BALL_PAD_M)
    rho_field = box._gel_density(msh, cell_tags)
    points = box._box_sample_points(msh, cell_tags, comm)
    assert np.array_equal(points, comm.bcast(points, root=0))
    rows = [(f"k6-{k}", (0.0, 0.0, dz), None, None, True)
            for k, dz in enumerate((0.0, 0.05, 0.1, -0.05))]
    out = box._box_readouts(msh, cell_tags, e_ccw, e_cw, sigma_f, omega, rho_field, rows, points, comm,
                            coverage=(centre, radius))
    ccw, cw = out["senses"]["ccw"], out["senses"]["cw"]
    b0 = abs(c) / (2.0 * omega)
    miss_b = {
        "ccw+": float(np.max(np.abs(ccw["b_plus"][ccw["valid"]] / b0 - 1.0))),
        "ccw-": float(np.max(np.abs(ccw["b_minus"][ccw["valid"]] / b0 - 1.0))),
        "cw+": float(np.max(np.abs(cw["b_plus"][cw["valid"]] / (2.0 * b0) - 1.0))),
    }
    # Printed beside: the 64 MHz closed form, which the 128 MHz field must miss by 1/2.
    b0_64 = abs(c) / (2.0 * 2.0 * np.pi * 64.0e6)
    miss_64 = float(np.max(np.abs(ccw["b_plus"][ccw["valid"]] / b0_64 - 1.0)))
    r_star = abs(out["sar_star_cw"] / out["sar_star"] - 4.0) / 4.0

    rho_403_zero = build_density_field(msh, PHANTOM_RHO_KG_PER_M3, cell_tags=cell_tags,
                                       density_map={403: 0.0})
    _w, tag3, p_miss_c, m_miss_zero = box._coverage_identity(e_ccw, sigma_f, rho_403_zero, cell_tags, 3,
                                                             centre, radius, comm)
    _w, _t, _p, m_miss_c = box._coverage_identity(e_ccw, sigma_f, rho_field, cell_tags, 3, centre,
                                                  radius, comm)
    v_all = out["tagged"]["volume_m3"]
    v_403 = v_all - tag3["volume_m3"]
    if comm.rank == 0:
        print(f"\n[GEO-38 C21 K6] omega 2pi x {omega / (2.0 * np.pi):.6e} Hz; {len(points)} tag-3 "
              f"centroids, {ccw['n_valid']} valid; b0 = |c|/(2 omega) {b0:.9e} T; (K6-B) misses "
              f"{miss_b} (band {DISCRETE_IDENTITY_RTOL:g}, ASSERTED); against the 64 MHz b0 "
              f"{b0_64:.9e} T the miss is {miss_64:.6f} (printed)")
        print(f"[GEO-38 C21 K6] c* {out['star_label']} {out['c_star']}, SAR10g {out['sar10g'].tolist()} "
              f"W/kg; cw/ccw ball at c* rel to 4 {r_star:.3e} (band {EXACT_IDENTITY_RTOL:g}, ASSERTED)")
        print(f"[GEO-38 C21 K6] (K6-i) radius {radius:.9f} m: power miss {out['power_miss']:.3e}, mass "
              f"miss {out['mass_miss']:.3e} (band {EXACT_IDENTITY_RTOL:g}, ASSERTED)")
        print(f"[GEO-38 C21 K6] (K6-c) tag 3 alone: ball rho(403) = 0, power miss {p_miss_c:.6e} (mass "
              f"then matches, {m_miss_zero:.3e}); ball rho(403) = gel, mass miss {m_miss_c:.6e}; V403/V "
              f"{v_403 / v_all:.6e} (> band, ASSERTED)", flush=True)
    assert len(points) > 0 and ccw["n_valid"] == len(points) and cw["n_valid"] == len(points)
    assert max(miss_b.values()) <= DISCRETE_IDENTITY_RTOL, miss_b
    assert out["star_label"] == "k6-2", out["star_label"]
    assert r_star <= EXACT_IDENTITY_RTOL, r_star
    assert out["power_miss"] <= EXACT_IDENTITY_RTOL, out["power_miss"]
    assert out["mass_miss"] <= EXACT_IDENTITY_RTOL, out["mass_miss"]
    assert p_miss_c > EXACT_IDENTITY_RTOL, p_miss_c
    assert m_miss_c > EXACT_IDENTITY_RTOL, m_miss_c


class _BuilderReached(Exception):
    pass


def _builder_sentinel(m, seen):
    """Replace the gmsh builder with a recorder that stores the kwargs plus the
    frequencies and records C4a's fixture would read, then raises (no gmsh)."""
    from fem_em_solver.io.mesh import MeshGenerator

    mods = _modules()

    def _sentinel(*_a, **kw):
        seen.append({"kw": dict(kw), "f_wf7": float(mods["wf7"].FREQUENCY_HZ),
                     "f_c4": float(mods["c4"].FREQUENCY_HZ_C4),
                     "cells": mods["s1a"].F_CLINICAL_R1_CELL_RECORD,
                     "unknowns": tuple(mods["s1a"].F_CLINICAL_R1_UNKNOWNS_RECORD)})
        raise _BuilderReached

    m.setattr(MeshGenerator, "birdcage_port_domain", staticmethod(_sentinel))


@_coding
def test_c21_k7_one_with_node_reaches_builder_with_3t_kwargs(monkeypatch):
    """(K7-3T): ``1`` + node 44.0615 runs `_run_box` into the gmsh builder with the
    3 T kwargs (rung Ø 0.704, C20's rod), f = 128 MHz and C20's records.  Negative
    control (ASSERTED, the brief's "must"): C11's ``s1a`` route alone leaves the box
    build on the 1.5 T kwargs and 64 MHz — the name-binding trap."""
    from tests.validation.test_geo36_rod_f_clinical_drive import S8_DRIVE_C_PF_ENV

    mods = _modules()
    t3, c20, s1a = mods["t3"], mods["c20"], mods["s1a"]
    kw15 = mods["i35"]._kwargs()  # the 1.5 T rod kwargs, unpatched
    kw3 = c20._kwargs_3t()
    comm = MPI.COMM_WORLD

    seen = []
    with monkeypatch.context() as m:
        for k in _foreign_envs():
            m.delenv(k, raising=False)
        m.setenv(C21_ENV, "1")
        m.setenv(S8_DRIVE_C_PF_ENV, f"{K_NODE_3T_PF}")
        _builder_sentinel(m, seen)
        with pytest.raises(_BuilderReached):
            test_geo38_c21_box_loaded_3t(m)

    neg = []
    with monkeypatch.context() as m:
        for k in _foreign_envs():
            m.delenv(k, raising=False)
        m.setattr(s1a, "_f_clinical_kwargs", t3._f_clinical_3t_kwargs)  # C11's route only
        _builder_sentinel(m, neg)
        with pytest.raises(_BuilderReached):
            mods["box"]._run_box(comm, m, K_NODE_3T_PF)

    assert len(seen) == 1 and len(neg) == 1, (len(seen), len(neg))
    got, got_neg = seen[0], neg[0]
    kw, kwn = got["kw"], got_neg["kw"]
    rung_d, rung_d_neg = 2.0 * kw["ring_radius"], 2.0 * kwn["ring_radius"]
    over = {k: (kw.get(k), v) for k, v in t3.F_CLINICAL_3T_OVERRIDES.items() if kw.get(k) != v}
    call_only = ("comm", "return_diagnostics")  # `_build_rod_f_clinical`'s own call arguments
    keys_same = sorted(k for k in kw if k not in call_only) == sorted(kw3)
    diff3 = [k for k in kw3 if k != "comm" and kw.get(k) != kw3[k]]
    print(f"\n[GEO-38 C21 K7] C21=1 node {K_NODE_3T_PF}: builder reached with rung Ø {rung_d:.12f} m "
          f"(3 T {t3.RUNG_DIAMETER_3T_M}), coil_length {kw['coil_length']}, shield_radius "
          f"{kw['shield_radius']}; rod {kw['implant_rod']}; f wf7 {got['f_wf7']:.6e} / C4 "
          f"{got['f_c4']:.6e} Hz; records cells {got['cells']}, U {got['unknowns']}; overrides "
          f"missed {over or 'none'}; vs C20 _kwargs_3t: keys same {keys_same}, differing "
          f"{diff3 or 'none'} (ASSERTED)")
    print(f"[GEO-38 C21 K7] control (C11's s1a route only): rung Ø {rung_d_neg:.12f} m (1.5 T "
          f"{2.0 * kw15['ring_radius']:.12f}), f wf7 {got_neg['f_wf7']:.6e} / C4 {got_neg['f_c4']:.6e} "
          f"Hz, records cells {got_neg['cells']} (ASSERTED: the 1.5 T kwargs)", flush=True)
    # (K7-3T)
    assert abs(rung_d / t3.RUNG_DIAMETER_3T_M - 1.0) <= 1e-15, rung_d
    assert not over, over
    assert keys_same and not diff3, (sorted(kw), diff3)
    assert kw["implant_rod"] == kw3["implant_rod"] and kw["phantom_box"] == kw3["phantom_box"]
    assert kw["ring_radius"] != kw15["ring_radius"]
    assert got["f_wf7"] == t3.FREQUENCY_3T_HZ == 128e6 and got["f_c4"] == t3.FREQUENCY_3T_HZ
    assert got["cells"] == c20.F_CLINICAL_3T_ROD_CELL_RECORD == 1_128_430
    assert got["unknowns"] == c20.F_CLINICAL_3T_ROD_UNKNOWNS_RECORD == (1_508_256, 7_783_028)
    # Negative control: the name-binding trap.
    assert kwn["ring_radius"] == kw15["ring_radius"] != kw["ring_radius"], kwn["ring_radius"]
    assert kwn["coil_length"] == kw15["coil_length"]
    assert got_neg["f_wf7"] == 64.0e6
    # Teardown restored the 1.5 T namespace.
    assert mods["s2"]._box_kwargs()["ring_radius"] == kw15["ring_radius"]


def _k8_stubs_3t(m, fx, json_dir, *, p_src_p2=1.0):
    """C6a's `_k8_stubs` (`:850`), then the 3 T overrides: env ω = 2π·128e6, the
    ESR specs at node 44.0615 pF and 128 MHz, and `_drive` on the 3 T loads (the
    same reciprocal synthetic circuit ``z_sym``).  Returns ``(lines, seen, z_sym, z_rc)``."""
    from dolfinx import fem

    from tests.validation import test_wf10_f_human_copper_esr_root as h17

    mods = _modules()
    box, i36, c4, t3 = mods["box"], mods["i36"], mods["c4"], mods["t3"]
    comm = MPI.COMM_WORLD
    lines, seen, z_sym = box._k8_stubs(m, fx, json_dir, p_src_p2=p_src_p2)
    env = dict(i36._drive_env(None, comm))
    env["omega"] = 2.0 * np.pi * t3.FREQUENCY_3T_HZ
    specs, _zrc, _zth = c4._esr_specs()
    z0 = 50.0
    z_rc = h17._z_rc(K_NODE_3T_PF * 1e-12, t3.FREQUENCY_3T_HZ)
    z_th = h17._thevenin_feed_impedance(z_rc, z0)
    loads = h17._parallel_loads(4, [0, 1], z_rc, z_th)
    _msh, _ct, e, _e2, sigma_f, _c, _half = fx
    e_2 = fem.Function(e.function_space)
    e_2.x.array[:] = 0.5j * e.x.array
    e_2.x.scatter_forward()
    columns = (e, e_2)
    feeds, order = ["P1", "P2"], ["P1", "P2", "P3", "P4"]

    def _drive(_env, _solver, _specs, driven_id, *, held):
        j = feeds.index(driven_id)
        v = np.zeros(4, dtype=np.complex128)
        v[j] = 1.0
        cur = h17._parallel_currents(z_sym, loads, v)
        return {"currents": dict(zip(order, cur)), "kind": "held" if held else "solve",
                "t_solve": 0.0, "t_post": 0.0, "e": columns[j].copy(), "sigma": sigma_f,
                "p_src": p_src_p2 if j == 1 else 1.0, "p_sheets": 0.5, "p_phantom": 0.2,
                "p_rest": 0.2, "p_surf": 0.1, "p_rod": 0.0,
                "p_sheet_by_id": {"P1": 0.1, "P2": 0.1, "P3": 0.15, "P4": 0.15}}

    m.setattr(i36, "_drive_env", lambda _fixture, _comm: env)
    m.setattr(c4, "_esr_specs", lambda *_a: (specs, z_rc, z_th))
    m.setattr(i36, "_drive", _drive)
    return lines, seen, z_sym, z_rc


def _k8_json_3t(comm, directory):
    import json

    rec = None
    if comm.rank == 0:
        with open(os.path.join(directory, f"c6-node-{K_NODE_3T_PF:.4f}.json")) as fh:
            rec = json.load(fh)
    return comm.bcast(rec, root=0)


@complex_only
@_coding
def test_c21_k8_run_box_end_to_end_on_stubs_at_128mhz(monkeypatch, tmp_path):
    """(K8-3T): `_run_box` under the 3 T patch set on C6a's K8 stubs at ω =
    2π·128e6 and node 44.0615 pF.  Column 2 = 0.5j × column 1, so ccw is 1.5 E and
    cw 0.5 E: B ratio 3, SAR ratio 9; mean |B₁⁺|_ccw = 1.5|c|/(2ω) at 128 MHz (the
    64 MHz value misses by 1/2, printed); the JSON keys and S2 against the circuit
    route at 128 MHz; the lattice-mean print equals the JSON lattice mean.  Control:
    column 2's P_src × (1 + 1e-3) raises naming (P) P2."""
    from tests.validation import test_wf10_f_human_copper_esr_root as h17
    from tests.validation.test_birdcage_power_identity import DISCRETE_IDENTITY_RTOL
    from tests.validation.test_birdcage_sar_mass_averaged import EXACT_IDENTITY_RTOL
    from tests.validation.test_wf10_f_human_tuned_b1_sar import _sense_weights

    mods = _modules()
    box, t3 = mods["box"], mods["t3"]
    comm = MPI.COMM_WORLD
    fx = box._k6_fixture(comm)
    c = fx[5]
    z0 = 50.0
    root = comm.bcast(str(tmp_path), root=0)
    good_dir, neg_dir = os.path.join(root, "good"), os.path.join(root, "neg")

    with monkeypatch.context() as m:
        _c21_patches(m)
        lines, seen, z_sym, z_rc = _k8_stubs_3t(m, fx, good_dir)
        box._run_box(comm, m, K_NODE_3T_PF)
        means = list(_LATTICE_MEANS)
    rec = _k8_json_3t(comm, good_dir)

    with monkeypatch.context() as m:
        _c21_patches(m)
        _k8_stubs_3t(m, fx, neg_dir, p_src_p2=1.0 + 1.0e-3)
        with pytest.raises(AssertionError, match=r"\(P\) P2") as neg:
            box._run_box(comm, m, K_NODE_3T_PF)
    rec_neg = _k8_json_3t(comm, neg_dir)

    omega = 2.0 * np.pi * t3.FREQUENCY_3T_HZ
    r_b = abs(rec["mean_b_plus"]["ccw"] / rec["mean_b_plus"]["cw"] / 3.0 - 1.0)
    r_s = abs(rec["sar_star"] / rec["sar_star_cw"] / 9.0 - 1.0)
    r_abs = abs(rec["mean_b_plus"]["ccw"] / (1.5 * abs(c) / (2.0 * omega)) - 1.0)
    r_abs_64 = abs(rec["mean_b_plus"]["ccw"] / (1.5 * abs(c) / (2.0 * 2.0 * np.pi * 64.0e6)) - 1.0)
    w = _sense_weights()["ccw"]
    expected_seen = [("P1", complex(w[0])), ("P2", complex(w[1]))]
    s_route = h17._circuit_current_route_s(z_sym, [0, 1], z_rc, z_rc, z0)
    s_json = np.array([[complex(*v) for v in row] for row in rec["s2"]])
    r_s2 = float(np.linalg.norm(s_json - s_route) / np.linalg.norm(s_route))
    z_rc_64 = h17._z_rc(K_NODE_3T_PF * 1e-12, 64.0e6)
    r_s2_64 = float(np.linalg.norm(s_json - h17._circuit_current_route_s(z_sym, [0, 1], z_rc_64,
                                                                         z_rc_64, z0))
                    / np.linalg.norm(s_route))

    by_sense = {s: (cm, lm, n) for s, cm, lm, n in means}
    json_lat = float(np.mean([p["b_plus_ccw"] for p in rec["lattice"]]))
    r_lat = abs(by_sense["ccw"][1] / json_lat - 1.0)
    r_cen = abs(by_sense["ccw"][0] / rec["mean_b_plus"]["ccw"] - 1.0)
    keys = {"c_pf", "feeds", "az_feeds", "s2", "y_in", "mean_b_plus", "mean_b_minus", "lattice",
            "c_star", "star_label", "sar_star", "sar_star_cw", "p_box", "p_src", "p_acc", "parts"}
    heads = ("(i) ball radius", "normalisers:", "price:", "(C21) mean |B1+|_ccw", "(C21) mean |B1+|_cw")
    missing_heads = [h for h in heads if not any(ln.startswith(h) for ln in lines)]
    if comm.rank == 0:
        print(f"\n[GEO-38 C21 K8] (K8-B) ccw/cw miss from 3 {r_b:.3e} (band "
              f"{DISCRETE_IDENTITY_RTOL:g}, ASSERTED); (K8-S) miss from 9 {r_s:.3e} (band "
              f"{EXACT_IDENTITY_RTOL:g}, ASSERTED)")
        print(f"[GEO-38 C21 K8] (K8-w) mean |B1+|_ccw {rec['mean_b_plus']['ccw']:.9e} T vs 1.5|c|/(2 omega) "
              f"at 128 MHz: miss {r_abs:.3e} (band {DISCRETE_IDENTITY_RTOL:g}, ASSERTED); vs the 64 MHz "
              f"value {r_abs_64:.6f} (printed)")
        print(f"[GEO-38 C21 K8] (K8-Ps) seen {seen}, JSON p_src {rec['p_src']!r}, c_pf {rec['c_pf']!r} "
              f"(ASSERTED)")
        print(f"[GEO-38 C21 K8] (K8-J) keys {sorted(rec)}; s2 vs circuit route at 128 MHz {r_s2:.3e} (band "
              f"{EXACT_IDENTITY_RTOL:g}, ASSERTED); vs the 64 MHz route {r_s2_64:.3e} (printed); (K8-L) "
              f"missing {missing_heads or 'none'}")
        print(f"[GEO-38 C21 K8] (Lat) printed lattice mean ccw {by_sense['ccw'][1]:.9e} T over "
              f"{by_sense['ccw'][2]} points vs JSON lattice mean {json_lat:.9e} T: rel {r_lat:.3e}; "
              f"printed centroid mean vs JSON mean_b_plus rel {r_cen:.3e} (<= 1e-12, ASSERTED)")
        print(f"[GEO-38 C21 K8] control: raised {str(neg.value).splitlines()[0]}; JSON still "
              f"written, parts[1] p_src {rec_neg['parts'][1]['p_src']!r} (ASSERTED)", flush=True)
    assert r_b <= DISCRETE_IDENTITY_RTOL, r_b
    assert r_s <= EXACT_IDENTITY_RTOL, r_s
    assert r_abs <= DISCRETE_IDENTITY_RTOL, r_abs
    assert seen == expected_seen, seen
    assert rec["p_src"] == 2.0, rec["p_src"]
    assert rec["c_pf"] == K_NODE_3T_PF, rec["c_pf"]
    assert set(rec) == keys, sorted(rec)
    assert r_s2 <= EXACT_IDENTITY_RTOL, r_s2
    assert not missing_heads, missing_heads
    assert r_lat <= 1.0e-12 and r_cen <= 1.0e-12, (r_lat, r_cen)
    assert rec_neg["parts"][1]["p_src"] == 1.0 + 1.0e-3, rec_neg["parts"][1]
    assert "(P) P2" in str(neg.value)


@complex_only
@_coding
def test_c21_k8_grid_3t_on_stubs(monkeypatch, tmp_path):
    """(K8-grid-3T), §9 item 819: C6b's (G1) on the 3 T box under `_c21_patches` (the
    same 90 of 91), then C6a's (K8-grid) under `_c21_patches` + `_k8_stubs_3t` at
    128 MHz, with control (G-off)."""
    mods = _modules()
    box = mods["box"]
    comm = MPI.COMM_WORLD
    with monkeypatch.context() as m:
        _c21_patches(m)
        box._g1_asserts(box._g1_check(" 3T"))
    fx = box._k6_fixture(comm)
    root = comm.bcast(str(tmp_path), root=0)
    keys = {"c_pf", "feeds", "az_feeds", "s2", "y_in", "mean_b_plus", "mean_b_minus", "lattice",
            "c_star", "star_label", "sar_star", "sar_star_cw", "p_box", "p_src", "p_acc", "parts"}
    rec, _off = box._k8_grid_check(comm, monkeypatch, fx, os.path.join(root, "on"),
                                   os.path.join(root, "off"), keys, "-3T", _c21_patches,
                                   lambda m, d: _k8_stubs_3t(m, fx, d), K_NODE_3T_PF, _k8_json_3t)
    assert rec["c_pf"] == K_NODE_3T_PF, rec["c_pf"]


# ---------------------------------------------------------------------------
# C52 (§9 item 929): the 3 T degree knob, routed on stubs.
# ---------------------------------------------------------------------------
#: C20's 3 T rod-mesh counts (`20261101T181021Z_GEO-38-C20-W1-build.log:11929`): E = U₁, F.
C20_E_F = (1_508_256, 2_383_258)


@_coding
def test_c52_kdeg_3t_knob_parsing_and_build_route(monkeypatch):
    """(K-deg-3T): unset / '' / 2 accepted, 1 / x / 3 refused, accepted under
    ``build``.  Under ``build`` at 2 the gmsh builder is reached through C4a's
    `_run_build` with s1a's records C20's and f 128 MHz, and C6a's `_run_build` is
    never called; with the knob unset, C6a's `_run_build` is called and C4a's is
    not.  A stray C6 degree knob under C21 ``1`` (C21 knob unset) is refused by
    `_run_box`.  (N-order): `_route_drive` before `_c21_patches` leaves s1a on the
    1.5 T rod pair."""
    from tests.validation.test_geo36_rod_f_clinical_drive import S8_DRIVE_C_PF_ENV

    mods = _modules()
    box, c4, s1a, i36, t3, c20 = (mods[k] for k in ("box", "c4", "s1a", "i36", "t3", "c20"))
    for k in _foreign_envs() + (C21_ENV, C21_DEGREE_ENV, S8_DRIVE_C_PF_ENV, box.C6_DEGREE_ENV):
        monkeypatch.delenv(k, raising=False)
    assert _c21_degree(None) is None and _c21_degree("1") is None and _c21_degree("build") is None
    monkeypatch.setenv(C21_DEGREE_ENV, "")
    assert _c21_degree("build") is None
    monkeypatch.setenv(C21_DEGREE_ENV, "2")
    assert _c21_degree("1") == 2 and _c21_degree("build") == 2 and _c21_degree(None) == 2
    for bad in ("1", "x", "3"):
        monkeypatch.setenv(C21_DEGREE_ENV, bad)
        with pytest.raises(ValueError, match=C21_DEGREE_ENV):
            _c21_degree("build")
    monkeypatch.delenv(C21_DEGREE_ENV)

    # build at 2: C4a's `_run_build` reaches the builder; C6a's `_run_build` never called.
    seen, box_calls = [], []
    with monkeypatch.context() as m:
        m.setenv(C21_ENV, "build")
        m.setenv(C21_DEGREE_ENV, "2")
        _builder_sentinel(m, seen)
        m.setattr(box, "_run_build", lambda *a, **k: box_calls.append(1))
        with pytest.raises(_BuilderReached):
            test_geo38_c21_box_loaded_3t(m)
    c4_env_after = os.environ.get(c4.C4_DEGREE_ENV)
    # build unset: C6a's `_run_build` called, C4a's not.
    box_calls_u, c4_calls_u = [], []
    with monkeypatch.context() as m:
        m.setenv(C21_ENV, "build")
        m.setattr(box, "_run_build", lambda *a, **k: box_calls_u.append(1))
        m.setattr(c4, "_run_build", lambda *a, **k: c4_calls_u.append(1))
        test_geo38_c21_box_loaded_3t(m)
    # A stray C6 degree knob under C21 = 1, the C21 knob unset.
    with monkeypatch.context() as m:
        m.setenv(C21_ENV, "1")
        m.setenv(S8_DRIVE_C_PF_ENV, f"{K_NODE_3T_PF}")
        m.setenv(box.C6_DEGREE_ENV, "2")
        _builder_sentinel(m, [])
        with pytest.raises(ValueError, match=box.C6_DEGREE_ENV) as stray:
            test_geo38_c21_box_loaded_3t(m)
    # (N-order): the reverse patch order.
    with monkeypatch.context() as m:
        i36._route_drive(m)
        _c21_patches(m)
        rec_wrong = tuple(s1a.F_CLINICAL_R1_UNKNOWNS_RECORD)
    got = seen[0] if seen else None
    if MPI.COMM_WORLD.rank == 0:
        print(f"\n[GEO-38 C52 K-deg-3T] {C21_DEGREE_ENV}: unset / '' / 2 accepted (build included), "
              f"1 / x / 3 refused; build+2: builder reached {len(seen)}x, s1a records U "
              f"{got and got['unknowns']}, cells {got and got['cells']}, f {got and got['f_wf7']}, "
              f"rung Ø {got and 2.0 * got['kw']['ring_radius']}; C6a _run_build calls {len(box_calls)}; "
              f"build unset: C6a {len(box_calls_u)} / C4a {len(c4_calls_u)}; stray C6 knob: "
              f"{stray.value}; (N-order) reverse order leaves s1a U {rec_wrong} (ASSERTED)", flush=True)
    assert len(seen) == 1, len(seen)
    assert got["unknowns"] == c20.F_CLINICAL_3T_ROD_UNKNOWNS_RECORD == (1_508_256, 7_783_028)
    assert got["cells"] == 1_128_430
    assert got["f_wf7"] == got["f_c4"] == t3.FREQUENCY_3T_HZ == 128e6
    assert abs(2.0 * got["kw"]["ring_radius"] / 0.704 - 1.0) <= 1e-15
    assert box_calls == [] and c4_env_after is None, (box_calls, c4_env_after)
    assert box_calls_u == [1] and c4_calls_u == [], (box_calls_u, c4_calls_u)
    assert rec_wrong != (1_508_256, 7_783_028) and rec_wrong == (1_606_744, 8_327_394), rec_wrong


def _k8_d2_3t(comm, monkeypatch, fx, json_dir, *, patched=True, report_degree=2):
    """(K8-3T-d2)'s run: [`_c21_patches`], `_k8_stubs_3t`, `box._c34_recorders`, then
    `box._run_box(degree=2)`.  Returns ``(lines, env_seen, degrees, env_after, exc)``."""
    from tests.validation import test_geo38_f_clinical_esr_root as c4

    box = _modules()["box"]
    exc = None
    with monkeypatch.context() as m:
        m.delenv(box.C6_DEGREE_ENV, raising=False)
        m.delenv(c4.C4_DEGREE_ENV, raising=False)
        if patched:
            _c21_patches(m)
        lines, _seen, _z_sym, _z_rc = _k8_stubs_3t(m, fx, json_dir)
        env_seen, degrees = box._c34_recorders(m, degree=report_degree)
        try:
            box._run_box(comm, m, K_NODE_3T_PF, degree=2)
        except AssertionError as e:
            exc = e
        env_after = os.environ.get(c4.C4_DEGREE_ENV)
    return lines, env_seen, degrees, env_after, exc


@complex_only
@_coding
def test_c52_k8_3t_d2_run_box_on_stubs(monkeypatch, tmp_path):
    """(K8-3T-d2): `_run_box(degree=2)` under `_c21_patches` on the 3 T K8 stubs —
    C4a's variable '2' at the fixture call, unset after; solver degrees [2]; the
    ``-d2`` JSON with K8's exact keys, the degree-1 name not written; (Deg) prints
    7783028; (K8-B) 3, (K8-S) 9, (K8-w) at 128 MHz.  Controls: (N-patch) without
    `_c21_patches` (Deg) prints 8327394, missing 7783028 by > 6e-2; (d2-c-3T) a
    fixture reporting degree 1 raises at the degree assert before any solver or file."""
    import json

    from tests.validation.test_birdcage_power_identity import DISCRETE_IDENTITY_RTOL
    from tests.validation.test_birdcage_sar_mass_averaged import EXACT_IDENTITY_RTOL

    mods = _modules()
    box, t3 = mods["box"], mods["t3"]
    comm = MPI.COMM_WORLD
    fx = box._k6_fixture(comm)
    c = fx[5]
    root = comm.bcast(str(tmp_path), root=0)
    dirs = {k: os.path.join(root, k) for k in ("d2", "npatch", "d2c")}
    lines, env_seen, degrees, env_after, exc = _k8_d2_3t(comm, monkeypatch, fx, dirs["d2"])
    lines_n, _e, _d, _a, exc_n = _k8_d2_3t(comm, monkeypatch, fx, dirs["npatch"], patched=False)
    _l, _e3, degrees_c, _a3, exc_c = _k8_d2_3t(comm, monkeypatch, fx, dirs["d2c"], report_degree=1)

    rec, d1_exists, files_c = None, None, None
    if comm.rank == 0:
        with open(os.path.join(dirs["d2"], f"c6-node-{K_NODE_3T_PF:.4f}-d2.json")) as fh:
            rec = json.load(fh)
        d1_exists = os.path.exists(os.path.join(dirs["d2"], f"c6-node-{K_NODE_3T_PF:.4f}.json"))
        files_c = sorted(os.listdir(dirs["d2c"])) if os.path.isdir(dirs["d2c"]) else []
    rec, d1_exists, files_c = comm.bcast((rec, d1_exists, files_c), root=0)
    omega = 2.0 * np.pi * t3.FREQUENCY_3T_HZ
    r_b = abs(rec["mean_b_plus"]["ccw"] / rec["mean_b_plus"]["cw"] / 3.0 - 1.0)
    r_s = abs(rec["sar_star"] / rec["sar_star_cw"] / 9.0 - 1.0)
    r_w = abs(rec["mean_b_plus"]["ccw"] / (1.5 * abs(c) / (2.0 * omega)) - 1.0)

    def _deg(ls):
        d = [ln for ln in ls if ln.startswith("(Deg) degree 2")]
        return d, (int(d[0].split("unknowns ")[1].split()[0]) if d else None)

    deg_lines, n_deg = _deg(lines)
    deg_lines_n, n_deg_n = _deg(lines_n)
    miss_n = abs(n_deg_n / 7_783_028 - 1.0) if n_deg_n else float("nan")
    if comm.rank == 0:
        print(f"\n[GEO-38 C52 K8-3T-d2] C4 variable at the fixture call {env_seen}, after {env_after!r}; "
              f"solver degrees {degrees}; JSON -d2 keys exact {set(rec) == box.K8_JSON_KEYS}, degree-1 "
              f"name written {d1_exists}; {deg_lines[0] if deg_lines else '(Deg) MISSING'} (ASSERTED)")
        print(f"[GEO-38 C52 K8-3T-d2] (K8-B) miss from 3 {r_b:.3e}, (K8-S) miss from 9 {r_s:.3e}, (K8-w) "
              f"at 128 MHz {r_w:.3e} (ASSERTED)")
        print(f"[GEO-38 C52 N-patch] without _c21_patches: "
              f"{deg_lines_n[0] if deg_lines_n else '(Deg) MISSING'}; |{n_deg_n}/7783028 - 1| = "
              f"{miss_n:.3e} (> 6e-02, ASSERTED)")
        print(f"[GEO-38 C52 d2-c-3T] degree 1 under degree=2 raised {exc_c!s}; solvers {degrees_c}, "
              f"files {files_c} (ASSERTED)", flush=True)
    assert exc is None and exc_n is None, (exc, exc_n)
    assert env_seen == ["2"] and env_after is None, (env_seen, env_after)
    assert degrees == [2], degrees
    assert set(rec) == box.K8_JSON_KEYS, sorted(rec)
    assert not d1_exists
    assert len(deg_lines) == 1 and n_deg == 7_783_028, deg_lines
    assert r_b <= DISCRETE_IDENTITY_RTOL and r_s <= EXACT_IDENTITY_RTOL and r_w <= DISCRETE_IDENTITY_RTOL
    assert n_deg_n == 8_327_394 != 7_783_028 and miss_n > 6.0e-2, (n_deg_n, miss_n)
    assert exc_c is not None and "(1, 2)" in str(exc_c), exc_c
    assert degrees_c == [] and files_c == [], (degrees_c, files_c)


def test_c52_deg_c_3t_degree1_count_is_not_degree2():
    """(Deg-c-3T): `_check_unknowns` refuses U₁ at degree 2 on C20's E / F and accepts
    2E + 2F = 7 783 028; U₂/U₁ > 5 on integers."""
    from tests.validation.test_geo38_f_clinical_matrix import _check_unknowns

    c20 = _modules()["c20"]
    u1, u2 = c20.F_CLINICAL_3T_ROD_UNKNOWNS_RECORD
    e, f = C20_E_F
    with pytest.raises(AssertionError):
        _check_unknowns(1_508_256, 2, 1_508_256, 2_383_258)
    assert _check_unknowns(u2, 2, e, f) == u2 == 2 * e + 2 * f == 7_783_028
    if MPI.COMM_WORLD.rank == 0:
        print(f"\n[GEO-38 C52 Deg-c-3T] 2E + 2F = {2 * e + 2 * f} = record {u2}; U1 {u1} at degree 2 "
              f"refused; U2/U1 {u2 / u1:.4f} > 5 (ASSERTED on integers)", flush=True)
    assert u1 == e
    assert u2 > 5 * u1


@_coding
def test_c52_k66_c21_call_alone(monkeypatch):
    """(K-66), C21's half: C21 ``1`` + ``C21_DEGREE=2`` + the node is accepted alone."""
    from tests.validation.test_geo36_rod_f_clinical_drive import S8_DRIVE_C_PF_ENV

    for k in _foreign_envs() + (C21_ENV, C21_DEGREE_ENV, S8_DRIVE_C_PF_ENV):
        monkeypatch.delenv(k, raising=False)
    monkeypatch.setenv(C21_ENV, "1")
    monkeypatch.setenv(S8_DRIVE_C_PF_ENV, f"{K_NODE_3T_PF:.4f}")
    monkeypatch.setenv(C21_DEGREE_ENV, "2")
    got = (_knobs(), _c21_degree("1"))
    if MPI.COMM_WORLD.rank == 0:
        print(f"\n[GEO-38 C52 K-66] C21=1 + {C21_DEGREE_ENV}=2 + node: {got} (ASSERTED)", flush=True)
    assert got == (("1", K_NODE_3T_PF), 2), got
