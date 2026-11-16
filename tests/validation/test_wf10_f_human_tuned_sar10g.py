"""`WF-10` step 2 part 2 (§9 item 359, §10 H9′): the 10 g SAR hotspot of the
tuned co-rotating drive on the copper F-human at C_in,F-human.

Item 358's builder (``build_tuned_quadrature``) supplies the fixture, (A0),
the node pair at ``C_IN_F_HUMAN_PF`` with (P_RC) and its P_surf-blind NEG, and
the two sense superpositions (``ccw`` = (1, −1j) from the convention code).
This module adds, on the 0.015 m phantom (≈ 1.8 cells across the 10 g ball —
the brief's named regime, never refined in-slot):

* (S_rep) — ``s_node`` reproduces 358's node JSON (read *before* the builder,
  which rewrites it) within ``REUSE_REPRODUCTION_RTOL``;
* (P_w) on both senses and (Σ_P), as in 358;
* (i) coverage — the whole-phantom ball equals ``mean_sar`` on the tag;
* (iii-a) mirror at c* and (iii-b) relative L2 on the lattice, in the
  feed-bisector plane az₁₇ + 45°, at ``C4_COVARIANCE_BAND``;
* (iv) containment — every ball 1 mm clear of the measured phantom bbox.

Printed only (rule (e), no SAR prior on this fixture): the single-feed P17
control at c* vs Mc* and its exact M-partner P21(Mc*); the azimuthal spread of
the four z_c rows; SAR per watt; the price.  Scope: c* is a lattice argmax; no
absolute SAR claim; no compliance claim.
"""

from __future__ import annotations

import json
import os
import time

import numpy as np
import pytest
from mpi4py import MPI

from tests.complex_mode import complex_only
from tests.validation.test_birdcage_b1_quadrature import _mirror_xy
from tests.validation.test_wf7_f_human_sar10g_hotspot import (
    CONTAINMENT_MARGIN_M,
    LATTICE_DZ_M,
    LATTICE_R_M,
)
from tests.validation.test_wf10_f_human_tuned_b1_sar import (
    C_IN_F_HUMAN_D2_PF,
    C_IN_F_HUMAN_PF,
    DEGREE_ENV,
    NODE_ENV,
    NODE_JSON,
    _degree_knob,
    _node_json_path,
    _node_pf,
)

STEP2_SAR_ENV = "FEM_EM_WF10_STEP2_SAR"
#: az₁₇ of 358's node JSON (``az_feeds[0]``), used by the ungated coding test
#: only; the gated test takes az₁₇ from the build.
AZ17_RECORD_DEG = 11.25
#: Nominal phantom (``…s2bi-W2.log:10769``), for (L0) only.
NOMINAL_R_M = 0.0642857
NOMINAL_HALF_H_M = 0.0857143
#: The extra azimuths of the z = z_c row (relative to az₁₇).
ROW_AZ_OFFSETS_DEG = (90.0, 180.0, 270.0)
# Predictions written first (brief), printed only.
PREDICTED_MIRROR = (0.5e-2, 3.0e-2)
PREDICTED_AZ_SPREAD = 0.10
PREDICTED_CONTROL_MISS = 0.30
PREDICTED_SREP = 1.0e-10


def _lattice(az17_deg, z_c):
    """The 18 candidate points: ``[(label, r, dz, d_az, (x, y, z))]``.

    z = z_c rows first (az₁₇, then +90 / +180 / +270), then az₁₇ at
    dz = ±50 mm, so the STOP-rule reduction keeps a prefix-independent subset.
    """
    out = []
    for d_az in (0.0,) + ROW_AZ_OFFSETS_DEG:
        phi = np.radians(az17_deg + d_az)
        for r in LATTICE_R_M:
            out.append((f"r={r * 1e3:.0f} dz=+0 az17{d_az:+.0f}", r, 0.0, d_az,
                        (r * np.cos(phi), r * np.sin(phi), z_c)))
    phi = np.radians(az17_deg)
    for dz in LATTICE_DZ_M:
        if dz == 0.0:
            continue
        for r in LATTICE_R_M:
            out.append((f"r={r * 1e3:.0f} dz={dz * 1e3:+.0f} az17+0", r, dz, 0.0,
                        (r * np.cos(phi), r * np.sin(phi), z_c + dz)))
    return out


def _plane_distance(points, az_deg):
    n = np.array([-np.sin(np.radians(az_deg)), np.cos(np.radians(az_deg)), 0.0])
    return np.abs(points @ n)


def _containment(c, centre, radius, half_height, a):
    """<= 0 means contained with the 1 mm margin (WF-7's ``contained``)."""
    rr = float(np.hypot(c[0] - centre[0], c[1] - centre[1]))
    radial = rr + a - (radius - CONTAINMENT_MARGIN_M)
    axial = abs(c[2] - centre[2]) + a - (half_height - CONTAINMENT_MARGIN_M)
    return max(radial, axial)


def _read_node_record(comm):
    """358's node record of the run's own degree, read BEFORE the builder (which
    rewrites it).  ``WF-10`` step 4b (§9 item 414).

    Unset knob: ``NODE_JSON`` (the degree-1 record), every check as item 359.
    ``FEM_EM_WF10_DEGREE=2``: the ``-d2`` sibling the b1_sar module's gated test
    writes earlier in the same window; absent ⇒ a BLOCKED assertion, never a
    fallback to the degree-1 file.  Returns ``(path, rec, s_rec)``.

    At degree 2 (Σ_P) stays exact: |E|² is quartic per cell and both the
    builder's ``_vol`` and ``mean_sar`` integrate with UFL's estimated degree
    (≥ 4), so the identity holds to round-off (predicted ≤ 1e-12).
    """
    knob = _degree_knob()
    c_node = _node_pf()  # step 4c (§9 item 511): the effective node, never a fallback
    path = _node_json_path(knob, c_pf=c_node)
    if knob is None:
        assert os.path.exists(path), f"{path} absent: item 359 stops (BLOCKED)"
    else:
        assert os.path.exists(path), (
            f"BLOCKED: {path} absent — at degree 2 (S_rep) reads the -d2 node record of node "
            f"{c_node} pF the b1_sar module's gated test (FEM_EM_WF10_STEP2=1) writes earlier "
            f"in the same window; never the degree-1 {NODE_JSON} nor another node's -d2 record")
    rec = None
    age_s = None
    if comm.rank == 0:
        with open(path) as fh:
            rec = json.load(fh)
        age_s = time.time() - os.path.getmtime(path)
    rec, age_s = comm.bcast((rec, age_s), root=0)
    if knob is not None:
        assert rec.get("degree") == 2, f"{path}: degree {rec.get('degree')!r}, expected 2"
        if comm.rank == 0:
            print(f"[WF-10 step4b] degree-2 node record {path}: age {age_s:.1f} s", flush=True)
    assert abs(rec["c_in_pf"] - c_node) == 0.0, (rec["c_in_pf"], c_node)
    s_rec = np.array([[complex(a, b) for a, b in row] for row in rec["s_node"]],
                     dtype=np.complex128)
    return path, rec, s_rec


def test_wf10_step4b_node_record_coding(monkeypatch, tmp_path):
    """(K0′): the 10 g module reads the node record of its own degree; a missing
    or wrong-degree ``-d2`` record raises, with a degree-1 file present beside it
    (no fallback).  numpy + JSON only, never the real ``logs/wf10/``."""
    import sys

    b1 = sys.modules["tests.validation.test_wf10_f_human_tuned_b1_sar"]
    me = sys.modules[__name__]
    monkeypatch.delenv(NODE_ENV, raising=False)  # step 4c: a window's knob never leaks in
    assert _node_json_path(None) == NODE_JSON  # unpatched, as imported
    d1 = str(tmp_path / f"h9p-node-{C_IN_F_HUMAN_PF:.4f}.json")
    monkeypatch.setattr(b1, "JSON_DIR", str(tmp_path))
    monkeypatch.setattr(b1, "NODE_JSON", d1)
    monkeypatch.setattr(me, "NODE_JSON", d1)
    comm = MPI.COMM_WORLD

    def plant(p, s, degree):
        with open(p, "w") as fh:
            json.dump({"c_in_pf": C_IN_F_HUMAN_PF, "degree": degree,
                       "s_node": [[[complex(v).real, complex(v).imag] for v in row]
                                  for row in s]}, fh)

    a = np.array([[0.1 + 0.2j, -0.3 + 0.05j], [-0.3 + 0.05j, 0.15 - 0.4j]])
    b = np.array([[0.11 + 0.19j, -0.29 + 0.06j], [-0.29 + 0.06j, 0.16 - 0.41j]])
    plant(d1, a, 1)

    # (a) knob unset: the degree-1 path, bit for bit.
    monkeypatch.delenv("FEM_EM_WF10_DEGREE", raising=False)
    path, _rec, s = _read_node_record(comm)
    assert path == d1 and np.array_equal(s, a), (path, s)

    # (b) knob 2, no -d2 file, degree-1 file present: BLOCKED, never a fallback.
    monkeypatch.setenv("FEM_EM_WF10_DEGREE", "2")
    with pytest.raises(AssertionError) as exc:
        _read_node_record(comm)
    assert "BLOCKED" in str(exc.value) and "-d2" in str(exc.value), str(exc.value)

    # (c) knob 2, a planted -d2 record: it is the one read.
    d2 = str(tmp_path / f"h9p-node-{C_IN_F_HUMAN_PF:.4f}-d2.json")
    plant(d2, b, 2)
    path, rec, s = _read_node_record(comm)
    assert path == d2 and path.endswith("-d2.json") and os.path.dirname(path) == str(tmp_path)
    assert rec["degree"] == 2 and np.array_equal(s, b), s

    # (d) knob 2, a -d2 record that says degree 1: the wrong-degree guard fires.
    plant(d2, b, 1)
    with pytest.raises(AssertionError) as exc:
        _read_node_record(comm)
    assert "degree" in str(exc.value), str(exc.value)
    print(f"[WF-10 step4b K0'] (a) unset -> {os.path.basename(d1)} bitwise; (b) missing -d2 "
          f"BLOCKED with the degree-1 file present; (c) -d2 read; (d) wrong degree raises "
          f"(all ASSERTED)", flush=True)


def test_wf10_step4c_node_record_coding(monkeypatch, tmp_path):
    """(N0′) step 4c (§9 item 511, H28): with ``FEM_EM_WF10_NODE=d2root`` and
    degree 2 the 10 g module reads only the 115.5209 ``-d2`` record — a missing
    one is BLOCKED (never the 118.2437 ``-d2`` or degree-1 record), a wrong-node
    one raises; with the knob unset K0′'s cases are unchanged.  numpy + JSON
    only, never the real ``logs/wf10/``."""
    import sys

    b1 = sys.modules["tests.validation.test_wf10_f_human_tuned_b1_sar"]
    me = sys.modules[__name__]
    monkeypatch.delenv(NODE_ENV, raising=False)
    monkeypatch.delenv(DEGREE_ENV, raising=False)
    d1 = str(tmp_path / f"h9p-node-{C_IN_F_HUMAN_PF:.4f}.json")
    monkeypatch.setattr(b1, "JSON_DIR", str(tmp_path))
    monkeypatch.setattr(b1, "NODE_JSON", d1)
    monkeypatch.setattr(me, "NODE_JSON", d1)
    comm = MPI.COMM_WORLD

    def plant(p, s, degree, c_pf):
        with open(p, "w") as fh:
            json.dump({"c_in_pf": c_pf, "degree": degree,
                       "s_node": [[[complex(v).real, complex(v).imag] for v in row]
                                  for row in s]}, fh)

    a = np.array([[0.1 + 0.2j, -0.3 + 0.05j], [-0.3 + 0.05j, 0.15 - 0.4j]])
    b = np.array([[0.11 + 0.19j, -0.29 + 0.06j], [-0.29 + 0.06j, 0.16 - 0.41j]])
    c = np.array([[0.12 + 0.18j, -0.28 + 0.07j], [-0.28 + 0.07j, 0.17 - 0.42j]])
    d2_old = str(tmp_path / f"h9p-node-{C_IN_F_HUMAN_PF:.4f}-d2.json")
    d2_new = str(tmp_path / f"h9p-node-{C_IN_F_HUMAN_D2_PF:.4f}-d2.json")
    plant(d1, a, 1, C_IN_F_HUMAN_PF)
    plant(d2_old, b, 2, C_IN_F_HUMAN_PF)

    monkeypatch.setenv(NODE_ENV, "d2root")
    monkeypatch.setenv(DEGREE_ENV, "2")
    # (a) only the 118.2437 -d2 (and degree-1) records present: BLOCKED, no fallback.
    with pytest.raises(AssertionError) as exc:
        _read_node_record(comm)
    msg = str(exc.value)
    assert "BLOCKED" in msg and "115.5209" in msg, msg

    # (b) a 115.5209 -d2 record whose c_in_pf is 118.2437: raises.
    plant(d2_new, c, 2, C_IN_F_HUMAN_PF)
    with pytest.raises(AssertionError):
        _read_node_record(comm)

    # (c) the correct record is the one read.
    plant(d2_new, c, 2, C_IN_F_HUMAN_D2_PF)
    path, rec, s = _read_node_record(comm)
    assert path == d2_new and path.endswith("115.5209-d2.json"), path
    assert rec["c_in_pf"] == C_IN_F_HUMAN_D2_PF and np.array_equal(s, c), s

    # (d) knob unset: K0′'s cases unchanged (degree 1 -> d1; degree 2 -> 118.2437 -d2).
    monkeypatch.delenv(NODE_ENV, raising=False)
    monkeypatch.delenv(DEGREE_ENV, raising=False)
    path, _rec, s = _read_node_record(comm)
    assert path == d1 and np.array_equal(s, a), path
    monkeypatch.setenv(DEGREE_ENV, "2")
    path, _rec, s = _read_node_record(comm)
    assert path == d2_old and np.array_equal(s, b), path
    print(f"[WF-10 step4c N0'] (a) only the {C_IN_F_HUMAN_PF} -d2 present -> BLOCKED; (b) "
          f"wrong-node record raises; (c) {os.path.basename(d2_new)} read; (d) unset -> K0' "
          f"paths unchanged (all ASSERTED)", flush=True)


def test_wf10_step2_sar_lattice_coding():
    """(L0): lattice size, mirror pairing, off-plane, nominal containment."""
    from fem_em_solver.post.sar import averaging_ball_radius

    from tests.validation.test_birdcage_b1_plus_map import PHANTOM_RHO_KG_PER_M3
    from tests.validation.test_mass_averaged_sar_standard_masses import TEN_GRAM_KG

    a10 = averaging_ball_radius(mass_kg=TEN_GRAM_KG, rho=PHANTOM_RHO_KG_PER_M3)
    for az17 in (AZ17_RECORD_DEG, 0.0, 33.7):
        lat = _lattice(az17, 0.0)
        assert len(lat) == 18
        pts = np.array([p for *_x, p in lat])
        az_m = az17 + 45.0
        img = _mirror_xy(pts, az_m)
        # z_c rows: index d*3 + i, d over (0, 90, 180, 270): 0 <-> 90, 180 <-> 270.
        rows = pts[:12].reshape(4, 3, 3)
        irow = img[:12].reshape(4, 3, 3)
        d_pair = max(np.max(np.abs(irow[0] - rows[1])), np.max(np.abs(irow[1] - rows[0])),
                     np.max(np.abs(irow[2] - rows[3])), np.max(np.abs(irow[3] - rows[2])))
        dist = _plane_distance(pts, az_m)
        bar = min(LATTICE_R_M) * np.sin(np.radians(45.0))
        worst = max(_containment(c, (0.0, 0.0, 0.0), NOMINAL_R_M, NOMINAL_HALF_H_M, a10)
                    for c in np.vstack([pts, img]))
        print(f"[WF-10 step2s L0] az17 {az17}: 18 points; mirror pairing {d_pair:.3e} m "
              f"(<= 1e-12 ASSERTED); min plane distance {dist.min() * 1e3:.4f} mm (>= "
              f"{bar * 1e3:.4f} mm ASSERTED); worst containment {worst * 1e3:.4f} mm (<= 0 "
              f"ASSERTED, a10 {a10 * 1e3:.4f} mm)")
        assert d_pair <= 1e-12, d_pair
        assert dist.min() >= bar - 1e-12, dist.min()
        assert worst <= 0.0, worst
    # The reduced (STOP) set is exactly the 12 z_c points.
    lat = _lattice(AZ17_RECORD_DEG, 0.0)
    assert sum(1 for x in lat if x[2] == 0.0) == 12


@complex_only
def test_wf10_step2_tuned_sar10g():
    if os.environ.get(STEP2_SAR_ENV, "").strip() != "1":
        pytest.skip(f"{STEP2_SAR_ENV} unset: heavy -n 8 F-human window, environment-selected")

    from fem_em_solver.post.sar import (
        averaging_ball_radius,
        build_density_field,
        mass_averaged_sar,
        mean_sar,
    )

    from tests.mesh.helpers import global_cell_tag_set
    from tests.validation.test_birdcage_b1_plus_map import (
        C4_COVARIANCE_BAND,
        PHANTOM_RHO_KG_PER_M3,
        _relative_l2,
    )
    from tests.validation.test_birdcage_f_human_rung import F_HUMAN_RING_RADIUS
    from tests.validation.test_birdcage_power_identity import DISCRETE_IDENTITY_RTOL
    from tests.validation.test_birdcage_sar_mass_averaged import EXACT_IDENTITY_RTOL
    from tests.validation.test_mass_averaged_sar_standard_masses import (
        QUADRATURE_DEGREE,
        TEN_GRAM_KG,
    )
    from tests.validation.test_port_birdcage_lumped_column import PHANTOM_CELL_TAG
    from tests.validation.test_port_birdcage_ring_matrix import REUSE_REPRODUCTION_RTOL
    from tests.validation.test_port_circuit_layer_field import _residual
    from tests.validation.test_port_drive_superposition import UNIT_WEIGHT_CONTROL_RTOL
    from tests.validation.test_wf7_f_human_sar10g_hotspot import (
        CONDUCTOR_CLEARANCE_M,
        COVERAGE_BALL_PAD_M,
        FIRST_BALL_STOP_S,
        _phantom_geometry,
    )
    from tests.validation.test_wf10_f_human_copper_esr_root import Y0_CODING_BAND
    from tests.validation.test_wf10_f_human_tuned_b1_sar import build_tuned_quadrature
    from tests.validation.test_wf12_f_human_hole import _rss_gib

    comm = MPI.COMM_WORLD

    # ---- 358's record of this degree, read BEFORE the builder (which rewrites it) --
    node_json, rec, s_rec = _read_node_record(comm)

    out = build_tuned_quadrature(comm, with_z50=False)
    say, ctx, times = out["say"], out["ctx"], out["times"]
    msh, cell_tags, sigma = ctx["msh"], ctx["cell_tags"], out["sigma"]

    # ---- (S_rep) -----------------------------------------------------------------
    s_node = np.asarray(out["s_node"], dtype=np.complex128)
    r_srep = _residual(s_node, s_rec)
    say(f"(S_rep) s_node vs {os.path.basename(node_json)} (degree {out.get('degree')}) at "
        f"{_node_pf()} pF: residual {r_srep:.6e} (band "
        f"{REUSE_REPRODUCTION_RTOL:g}, ASSERTED; PREDICTED <= {PREDICTED_SREP:g}); JSON "
        f"az_feeds {rec['az_feeds']} vs run {list(out['az_feeds'])}")

    # ---- geometry, density, ball ---------------------------------------------------
    centre, radius, half_height, lo, hi = _phantom_geometry(msh, cell_tags, PHANTOM_CELL_TAG, comm)
    z_c = float(centre[2])
    a10 = averaging_ball_radius(mass_kg=TEN_GRAM_KG, rho=PHANTOM_RHO_KG_PER_M3)
    az17 = float(out["az_feeds"][0])
    az_m = float(out["az_mirror"])
    say(f"phantom bbox lo {lo} hi {hi}: centre ({centre[0]:.6e}, {centre[1]:.6e}, "
        f"{centre[2]:.6e}) m, R {radius * 1e3:.4f} mm, H/2 {half_height * 1e3:.4f} mm; 10 g "
        f"a = {a10 * 1e3:.4f} mm; az17 {az17:.9f} deg, mirror plane {az_m:.9f} deg")
    all_tags = global_cell_tag_set(msh, cell_tags)
    rho_field = build_density_field(
        msh, PHANTOM_RHO_KG_PER_M3, cell_tags=cell_tags,
        density_map={int(t): 0.0 for t in sorted(all_tags) if int(t) != PHANTOM_CELL_TAG})

    fields = {"ccw": out["e"]["ccw"], "cw": out["e"]["cw"], "P17": out["e17"],
              "P21": out["e21"]}
    containment = []
    walls = []

    def ball(key, c, label):
        c = tuple(float(v) for v in c)
        containment.append((label, c, _containment(c, centre, radius, half_height, a10)))
        t = time.perf_counter()
        r = mass_averaged_sar(fields[key], sigma=sigma, rho=rho_field, center=c, radius=a10,
                              comm=comm, quadrature_degree=QUADRATURE_DEGREE)
        r["wall_s"] = comm.allreduce(time.perf_counter() - t, op=MPI.MAX)
        walls.append(r["wall_s"])
        say(f"  ball {key:<3s} {label:<30s} at ({c[0] * 1e3:8.3f}, {c[1] * 1e3:8.3f}, "
            f"{c[2] * 1e3:8.3f}) mm: SAR10g {r['averaged_sar_w_per_kg']:.9e} W/kg, mass "
            f"{r['mass_kg'] * 1e3:.6f} g, {r['wall_s']:.2f} s")
        return r

    # ---- 1. ccw lattice (STOP rule) -------------------------------------------------
    comm.Barrier()
    t_sar = time.perf_counter()
    lattice = _lattice(az17, z_c)
    ccw = {}
    reduced = False
    for idx, (label, _r, dz, _d, c) in enumerate(lattice):
        if reduced and dz != 0.0:
            continue
        ccw[label] = ball("ccw", c, label)
        if idx == 0 and ccw[label]["wall_s"] > FIRST_BALL_STOP_S:
            reduced = True
            say(f"REDUCED LATTICE: first ball {ccw[label]['wall_s']:.2f} s > "
                f"{FIRST_BALL_STOP_S:.0f} s — z = z_c rows only (pre-registered, not a red)")
    coords = {label: np.array(c) for label, *_x, c in lattice}
    used = [label for label, *_x in lattice if label in ccw]
    # ---- 2. c* -------------------------------------------------------------------
    key_star = max(used, key=lambda k: ccw[k]["averaged_sar_w_per_kg"])
    sar_star = ccw[key_star]["averaged_sar_w_per_kg"]
    c_star = coords[key_star]
    # ---- 3. cw at mirrored points ---------------------------------------------------
    cw = {}
    for label in used:
        mc = _mirror_xy(coords[label][None, :], az_m)[0]
        cw[label] = ball("cw", mc, f"M({label})")
    # ---- 4. controls -------------------------------------------------------------
    mc_star = _mirror_xy(c_star[None, :], az_m)[0]
    ctl17_c = ball("P17", c_star, "control P17 c*")["averaged_sar_w_per_kg"]
    ctl17_m = ball("P17", mc_star, "control P17 Mc*")["averaged_sar_w_per_kg"]
    ctl21_m = ball("P21", mc_star, "control P21 Mc*")["averaged_sar_w_per_kg"]
    # ---- 5. coverage + mean_sar ------------------------------------------------------
    coverage_radius = float(np.hypot(radius, half_height) + COVERAGE_BALL_PAD_M)
    conductor_limit = F_HUMAN_RING_RADIUS - CONDUCTOR_CLEARANCE_M
    t = time.perf_counter()
    whole = mass_averaged_sar(fields["ccw"], sigma=sigma, rho=rho_field,
                              center=tuple(float(v) for v in centre), radius=coverage_radius,
                              comm=comm, quadrature_degree=QUADRATURE_DEGREE)
    tagged = mean_sar(fields["ccw"], sigma=sigma, rho=PHANTOM_RHO_KG_PER_M3,
                      cell_tags=cell_tags, subdomain_ids=PHANTOM_CELL_TAG, comm=comm)
    t_cov = comm.allreduce(time.perf_counter() - t, op=MPI.MAX)
    t_sar = comm.allreduce(time.perf_counter() - t_sar, op=MPI.MAX)
    power_miss = abs(whole["dissipated_power_w"] / tagged["dissipated_power_w"] - 1.0)
    mass_miss = abs(whole["mass_kg"] / (PHANTOM_RHO_KG_PER_M3 * tagged["volume_m3"]) - 1.0)
    p_ph = out["parts"]["ccw"]["p_phantom"]
    p_src = out["parts"]["ccw"]["p_src"]
    r_sigma = abs(tagged["dissipated_power_w"] - p_ph) / abs(p_ph)

    # ---- readings --------------------------------------------------------------------
    sar_c = np.array([ccw[k]["averaged_sar_w_per_kg"] for k in used])
    sar_m = np.array([cw[k]["averaged_sar_w_per_kg"] for k in used])
    mask = np.ones(len(used), dtype=bool)
    mirror_star = abs(cw[key_star]["averaged_sar_w_per_kg"] - sar_star) / sar_star
    mirror_l2 = _relative_l2(sar_m, sar_c, mask)
    ctl_miss = abs(ctl17_m - ctl17_c) / ctl17_c
    ctl_partner = abs(ctl21_m - ctl17_c) / ctl17_c
    spreads = {}
    for r in LATTICE_R_M:
        vals = [ccw[lab]["averaged_sar_w_per_kg"] for lab, rr, dz, _d, _c in lattice
                if rr == r and dz == 0.0]
        spreads[r] = (max(vals) - min(vals)) / np.mean(vals)
    worst_cont = max(x[2] for x in containment)
    _m, _rss_per, rss_sum = _rss_gib(comm)

    if comm.rank == 0:
        print(f"[WF-10 step2s] ccw 10 g table ({len(used)} points"
              f"{', REDUCED' if reduced else ''}); cw at the feed-bisector image:", flush=True)
        for k in used:
            d = abs(cw[k]["averaged_sar_w_per_kg"] / ccw[k]["averaged_sar_w_per_kg"] - 1.0)
            print(f"[WF-10 step2s]   {k:<26s} ccw {ccw[k]['averaged_sar_w_per_kg']:.9e} cw(M) "
                  f"{cw[k]['averaged_sar_w_per_kg']:.9e} W/kg, rel {d * 100:.4f}%"
                  f"{'   <- c*' if k == key_star else ''}", flush=True)
        print(f"[WF-10 step2s] c* = {key_star} ({c_star[0] * 1e3:.3f}, {c_star[1] * 1e3:.3f}, "
              f"{c_star[2] * 1e3:.3f}) mm, SAR10g(c*) {sar_star:.9e} W/kg (PREDICTED r = 45 mm, "
              f"z_c; reading)", flush=True)
        print(f"[WF-10 step2s] (iii-a) mirror at c*: {mirror_star * 100:.4f}% (band "
              f"{C4_COVARIANCE_BAND * 100:.1f}%, ASSERTED; PREDICTED {PREDICTED_MIRROR[0] * 100:g}-"
              f"{PREDICTED_MIRROR[1] * 100:g}%)", flush=True)
        print(f"[WF-10 step2s] (iii-b) mirror lattice rel L2 over {len(used)} points: "
              f"{mirror_l2 * 100:.4f}% (band {C4_COVARIANCE_BAND * 100:.1f}%, ASSERTED; PREDICTED "
              f"{PREDICTED_MIRROR[0] * 100:g}-{PREDICTED_MIRROR[1] * 100:g}%)", flush=True)
        print("[WF-10 step2s] azimuthal spread of the four z_c rows (reading, not C16; PREDICTED "
              f"<= {PREDICTED_AZ_SPREAD * 100:g}%): "
              + ", ".join(f"r={r * 1e3:.0f} mm {v * 100:.4f}%" for r, v in spreads.items()),
              flush=True)
        print(f"[WF-10 step2s] control (printed, PREDICTED >= {PREDICTED_CONTROL_MISS * 100:g}%): "
              f"P17 SAR10g(c*) {ctl17_c:.6e}, P17 SAR10g(Mc*) {ctl17_m:.6e} W/kg, miss "
              f"{ctl_miss * 100:.4f}%; M-partner P21(Mc*) {ctl21_m:.6e} vs P17(c*): "
              f"{ctl_partner * 100:.4f}% (PREDICTED inside {C4_COVARIANCE_BAND * 100:.1f}%)",
              flush=True)
        print(f"[WF-10 step2s] (i) coverage ball {coverage_radius * 1e3:.4f} mm (< conductor "
              f"limit {conductor_limit * 1e3:.1f} mm ASSERTED): power "
              f"{whole['dissipated_power_w']:.12e} vs tagged {tagged['dissipated_power_w']:.12e} W, "
              f"rel {power_miss:.3e}; mass {whole['mass_kg']:.12e} vs rho*V "
              f"{PHANTOM_RHO_KG_PER_M3 * tagged['volume_m3']:.12e} kg, rel {mass_miss:.3e} (band "
              f"{EXACT_IDENTITY_RTOL:g}, ASSERTED); {t_cov:.2f} s", flush=True)
        print(f"[WF-10 step2s] (Sigma_P) P_phantom {p_ph:.12e} vs mean_sar dissipated "
              f"{tagged['dissipated_power_w']:.12e} W, rel {r_sigma:.3e} (band "
              f"{EXACT_IDENTITY_RTOL:g}, ASSERTED)", flush=True)
        for sense in ("ccw", "cw"):
            rel_w, blind_w = out["anchors"]["p_w"][sense]
            r_lin, r_y = out["anchors"]["lin"][sense]
            print(f"[WF-10 step2s] [{sense}] (P_w) {rel_w:.3e} (band {DISCRETE_IDENTITY_RTOL:g}, "
                  f"ASSERTED), blind {blind_w:.4e}; (Lin) {r_lin:.3e} / Y {r_y:.3e} (ASSERTED)",
                  flush=True)
        print(f"[WF-10 step2s] (P_RC) {out['p_rc_node']:.3e}; node NEG blind "
              f"{out['blind_node']:.4e} (both ASSERTED in the builder)", flush=True)
        print(f"[WF-10 step2s] (iv) worst containment {worst_cont * 1e3:.4f} mm over "
              f"{len(containment)} balls (<= 0 ASSERTED, measured bbox)", flush=True)
        print(f"[WF-10 step2s] SAR per watt (ccw): SAR10g(c*)/P_src {sar_star / p_src:.6e}, "
              f"/P_phantom {sar_star / p_ph:.6e} W/kg/W; whole-phantom mean "
              f"{tagged['mean_sar_w_per_kg']:.6e} W/kg, /P_src "
              f"{tagged['mean_sar_w_per_kg'] / p_src:.6e}, /P_phantom "
              f"{tagged['mean_sar_w_per_kg'] / p_ph:.6e} W/kg/W; peak-to-mean "
              f"{sar_star / tagged['mean_sar_w_per_kg']:.4f}; P_phantom/P_src {p_ph / p_src:.6e}",
              flush=True)
        print(f"[WF-10 step2s] PRICE: {len(walls)} 10 g balls + coverage, per ball mean "
              f"{np.mean(walls):.2f} s max {np.max(walls):.2f} s; SAR phase {t_sar:.1f} s; build "
              f"{times['build']:.1f} s, node pair {times['node_pair']:.1f} s, superposition "
              f"{times['superposition']:.1f} s; RSS summed {rss_sum:.2f} GiB; window "
              f"{time.perf_counter() - out['t_window']:.1f} s at -n {comm.size}", flush=True)

    # ---- asserts ---------------------------------------------------------------------
    assert r_srep <= REUSE_REPRODUCTION_RTOL, f"(S_rep) {r_srep:.6e}"
    for sense in ("ccw", "cw"):
        rel_w, _b = out["anchors"]["p_w"][sense]
        r_lin, r_y = out["anchors"]["lin"][sense]
        assert rel_w <= DISCRETE_IDENTITY_RTOL, f"(P_w) {sense} {rel_w:.3e} — never re-band"
        assert r_lin <= UNIT_WEIGHT_CONTROL_RTOL, f"(Lin) {sense} {r_lin:.3e}"
        assert r_y <= Y0_CODING_BAND, f"(Lin) Y_in {sense} {r_y:.3e}"
    assert r_sigma <= EXACT_IDENTITY_RTOL, f"(Sigma_P) {r_sigma:.3e}"
    assert coverage_radius < conductor_limit, (coverage_radius, conductor_limit)
    assert power_miss <= EXACT_IDENTITY_RTOL, f"(i) power {power_miss:.3e}"
    assert mass_miss <= EXACT_IDENTITY_RTOL, f"(i) mass {mass_miss:.3e}"
    assert worst_cont <= 0.0, f"(iv) containment {worst_cont:.4e} m"
    assert mirror_star <= C4_COVARIANCE_BAND, f"(iii-a) {mirror_star:.4e} — never re-band"
    assert mirror_l2 <= C4_COVARIANCE_BAND, f"(iii-b) {mirror_l2:.4e} — never re-band"
