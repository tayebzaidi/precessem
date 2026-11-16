"""`WF-16` step 9 part A (§9 item 382, chain U36): the tuned copper F-human
through the ``Study``.

``Study(FHumanBirdcageGeometry(as_hole=True), gel,
TwoPortQuadratureDrive(C_in, esr_ohm=0.2), conductor=LeontovichWall(),
frequency_hz=64e6, degree=1)`` reproduces H9′'s node (`WF-10` step 2, §9
item 358): the feeds' node S, the active input admittance per sense and the
two superposed sense partitions, against the stored node JSON at
``REUSE_REPRODUCTION_RTOL``.  The Study never calls ``build_tuned_quadrature``
(that builder rewrites the JSON); the JSON is read on rank 0 and broadcast
before any solve.

Anchors (asserted, every band imported, none set here):
  (K0)  no mesh: the hole route's generator kwargs equal `WF-7`'s
        ``_f_human_kwargs(as_hole=True)``; the default and boxed dicts are
        unchanged; the drive's refusals; ``Z_RC`` / ``Z_th`` / ``V_s`` equal the
        H17 helpers' values on the same inputs exactly.
  (A0)  cells == ``HOLE_CELL_RECORD`` within ``CELL_COUNT_BAND`` and == the
        npz ``cells``; az21 − az17 = 90° within ``AZIMUTH_MATCH_DEG``;
        ``az_mirror`` against the JSON's.
  (U)   the Study's 32 specs against the builder's construction re-done
        in-test on the Study's mesh (``h``, ``w`` within ``IDENTITY_BAND``,
        every other field ``==``).
  (W)   the sense weights against the JSON ``w`` and ``SENSES["(1,-1j)"]``
        within ``Y0_CODING_BAND``.
  (S_rep) ``_residual(s_node, s_json)`` <= ``REUSE_REPRODUCTION_RTOL``.
  (Y)   per sense, ||y_in − y_json|| / ||y_json|| <= ``REUSE_REPRODUCTION_RTOL``.
  (P_w) per sense, the superposed-drive partition <= ``DISCRETE_IDENTITY_RTOL``;
        P_src and P_phantom against the JSON at ``REUSE_REPRODUCTION_RTOL``.

Negative controls (rule (e)):
  asserted, backed: the P_surf-blind miss of (P_w) ccw >= ``BLIND_FACTOR`` x
  ``DISCRETE_IDENTITY_RTOL`` (record 7.0397e-03).
  predicted, printed only: (S_rep) against the H17 node JSON at 118.2000 pF.

Opt-in: ``FEM_EM_WF16_S9=1`` (unset ⇒ the solve nodes skip; any other value
raises); width 8 asserted when set.  Scope: degree 1, 64 MHz, κ off, one
copper hole fixture, one phantom; not f₀ / Q / match / S₁₁ / tuning, not
|B₁⁺| / SAR (part B), not degree 2, not AED.
"""

from __future__ import annotations

import dataclasses
import json
import os
import time

import numpy as np
import pytest
from mpi4py import MPI

from fem_em_solver import (
    FHumanBirdcageGeometry,
    FourPortBirdcageGeometry,
    LeontovichWall,
    Study,
    TwoPortQuadratureDrive,
)
from fem_em_solver.materials import GelledSalinePhantomMaterial

from tests.complex_mode import complex_only
from tests.mesh.test_birdcage_leg_offset import IDENTITY_BAND
from tests.mesh.test_birdcage_port_sheet_prerequisite import CELL_COUNT_BAND
from tests.validation.test_birdcage_power_identity import (
    DISCRETE_IDENTITY_RTOL,
    _sheet_field_dissipation_w,
    _source_power_w,
)
from tests.validation.test_lossy_sphere_fullwave import SALINE_EPSILON_R, SALINE_SIGMA
from tests.validation.test_port_birdcage_lumped_column import PHANTOM_CELL_TAG
from tests.validation.test_port_birdcage_ring_column import AZIMUTH_MATCH_DEG
from tests.validation.test_port_birdcage_ring_matrix import REUSE_REPRODUCTION_RTOL
from tests.validation.test_port_circuit_layer_field import _residual
from tests.validation.test_th14_birdcage_copper import _surface_loss_w
from tests.validation.test_wf7_f_human_b1_quadrature import (
    DEGREE,
    FREQUENCY_HZ,
    _f_human_kwargs,
)
from tests.validation.test_wf10_f_human_copper_esr_root import (
    ESR_OHM,
    JSON_DIR,
    NPZ_PATH,
    SENSES,
    Y0_CODING_BAND,
    _active_y_in,
    _source_voltage,
    _thevenin_feed_impedance,
    _z_rc,
)
from tests.validation.test_wf10_f_human_tuned_b1_sar import C_IN_F_HUMAN_PF, NODE_JSON
from tests.validation.test_wf12_f_human_copper_matrix import (
    BLIND_FACTOR,
    COPPER_SIGMA,
    HOLE_CELL_RECORD,
)

S9_ENV = "FEM_EM_WF16_S9"
RECORD_RANK_WIDTH = 8  # the width H9′'s node JSON was recorded at
H17_NODE_JSON = os.path.join(JSON_DIR, "h17-node-118.2000.json")
#: 358's node line (`20261010T170356Z_WF-10-step2-W1.log`), printed beside.
NODE_LINE_RECORD = ("Z_RC 2.000000000e-01-2.103110767e+01j, Z_th 7.634949838e+00"
                    "-1.774868389e+01j, V_s 1.022606662e+00+2.377215668e+00j")
#: (S_rep)'s same-comparison prior (`20261010T171705Z_WF-10-step2s-W1.log:6312`).
S_REP_PRIOR = 6.027027e-13
#: (P_w) priors (`20261010T170356Z_WF-10-step2-W1.log`, [ccw] / [cw] lines).
P_W_PRIOR = {"ccw": 5.065e-13, "cw": 5.893e-13}


def _knob() -> bool:
    raw = os.environ.get(S9_ENV, "").strip()
    if raw == "":
        return False
    if raw == "1":
        return True
    raise ValueError(f"{S9_ENV}={raw!r}: only unset/empty or '1' are accepted")


def _c_f() -> float:
    # The builder's own expression (one ulp from a literal 118.2437e-12).
    return C_IN_F_HUMAN_PF * 1e-12


def _phantom():
    return GelledSalinePhantomMaterial(SALINE_SIGMA, SALINE_EPSILON_R, FREQUENCY_HZ)


def _study(geom=None, drive=None, conductor="cu"):
    return Study(
        FHumanBirdcageGeometry(as_hole=True) if geom is None else geom,
        _phantom(),
        TwoPortQuadratureDrive(_c_f(), esr_ohm=ESR_OHM) if drive is None else drive,
        frequency_hz=FREQUENCY_HZ,
        degree=DEGREE,
        conductor=LeontovichWall(COPPER_SIGMA) if conductor == "cu" else conductor,
    )


def _cx(pair):
    return complex(pair[0], pair[1])


def test_wf16_s9a_api_and_refusals():
    """(K0): no mesh anywhere in this node."""
    _knob()  # a bad knob value raises here too
    assert FHumanBirdcageGeometry(as_hole=True).generator_kwargs() == \
        _f_human_kwargs(as_hole=True)
    assert FHumanBirdcageGeometry(as_hole=True, emit_leg_ports=True).generator_kwargs() == \
        _f_human_kwargs(as_hole=True)
    assert FHumanBirdcageGeometry().generator_kwargs() == _f_human_kwargs()
    assert "as_hole" not in FHumanBirdcageGeometry().generator_kwargs()
    assert FHumanBirdcageGeometry(emit_leg_ports=True).generator_kwargs() == \
        _f_human_kwargs(emit_leg_ports=True)
    assert FHumanBirdcageGeometry().as_hole is False

    drv = TwoPortQuadratureDrive(_c_f(), esr_ohm=ESR_OHM)
    assert drv.feeds == ("P17", "P21") and drv.z0_ohm == 50.0
    z_rc = _z_rc(_c_f(), FREQUENCY_HZ)
    assert drv.gap_impedance_ohm(FREQUENCY_HZ) == z_rc
    assert drv.feed_impedance_ohm(FREQUENCY_HZ) == _thevenin_feed_impedance(z_rc, 50.0)
    assert drv.source_voltage_v(FREQUENCY_HZ) == _source_voltage(z_rc, 50.0)
    print(f"[WF-16 s9a K0] Z_RC {drv.gap_impedance_ohm(FREQUENCY_HZ):.9e} Z_th "
          f"{drv.feed_impedance_ohm(FREQUENCY_HZ):.9e} V_s "
          f"{drv.source_voltage_v(FREQUENCY_HZ):.9e}; record {NODE_LINE_RECORD}")
    for bad in (dict(c_in_f=0.0), dict(c_in_f=float("nan")), dict(c_in_f=-1e-12),
                dict(c_in_f=1e-10, esr_ohm=-0.1), dict(c_in_f=1e-10, esr_ohm=float("inf")),
                dict(c_in_f=1e-10, feeds="P17"), dict(c_in_f=1e-10, feeds=("P17",)),
                dict(c_in_f=1e-10, feeds=("P17", "P17")), dict(c_in_f=1e-10, z0_ohm=0.0)):
        with pytest.raises(ValueError):
            TwoPortQuadratureDrive(**bad)
    with pytest.raises(dataclasses.FrozenInstanceError):
        drv.esr_ohm = 1.0
    # Refusals before any mesh, through item 373's hook.
    with pytest.raises(TypeError):
        _study(geom=FourPortBirdcageGeometry(as_hole=True)).solve()
    with pytest.raises(ValueError):
        _study(geom=FHumanBirdcageGeometry()).solve()
    with pytest.raises(TypeError):
        _study(conductor=None).solve()
    with pytest.raises(ValueError):
        _study(drive=TwoPortQuadratureDrive(_c_f(), feeds=("P17", "P1"))).solve()


@pytest.fixture(scope="module")
def tuned():
    if not _knob():
        pytest.skip(f"{S9_ENV} unset: the tuned F-human node is a heavy -n 8 window, "
                    "environment-selected")
    comm = MPI.COMM_WORLD
    assert comm.size == RECORD_RANK_WIDTH, (comm.size, RECORD_RANK_WIDTH)
    rec, npz_cells, h17 = None, None, None
    if comm.rank == 0:
        assert os.path.exists(NODE_JSON), f"{NODE_JSON} absent: stop (BLOCKED, rule (d))"
        assert os.path.exists(NPZ_PATH), f"{NPZ_PATH} absent: stop (BLOCKED, rule (d))"
        with open(NODE_JSON) as fh:
            rec = json.load(fh)
        with np.load(NPZ_PATH, allow_pickle=False) as npz:
            npz_cells = int(npz["cells"])
        if os.path.exists(H17_NODE_JSON):
            with open(H17_NODE_JSON) as fh:
                h17 = json.load(fh)
        st = os.stat(NODE_JSON)
        print(f"[WF-16 s9a] node JSON {NODE_JSON}: {st.st_size} bytes, mtime "
              f"{time.strftime('%Y-%m-%d %H:%M:%S %Z', time.localtime(st.st_mtime))}; "
              f"c_f {rec['c_f']!r}", flush=True)
    rec, npz_cells, h17 = comm.bcast((rec, npz_cells, h17), root=0)
    assert rec["c_f"] == _c_f(), (rec["c_f"], _c_f())
    assert rec["r_ohm"] == ESR_OHM and rec["degree"] == 1
    t0 = time.perf_counter()
    st = _study()
    res = st.solve(comm)
    elapsed = time.perf_counter() - t0
    assert st.result is res
    if comm.rank == 0:
        print(f"[WF-16 s9a] solve: build {res.times['build']:.1f} s, node pair "
              f"{res.times['node_pair']:.2f} s, superposition {res.times['superposition']:.2f} s; "
              f"total {elapsed:.1f} s (H9′ JSON times {rec['times_s']})", flush=True)
    return {"st": st, "res": res, "rec": rec, "npz_cells": npz_cells, "h17": h17, "comm": comm}


@complex_only
def test_wf16_s9a_fixture_and_specs(tuned):
    """(A0) and (U)."""
    from fem_em_solver.io.mesh import _interface_facet_tags
    from tests.mesh.test_birdcage_port_sheets import PORT_LOWER, PORT_UPPER, SHEET_IFACE
    from tests.mesh.test_birdcage_port_terminals import _interface_area_or_zero
    from tests.mesh.test_birdcage_ring_gaps_scaleup import _ring_gap_frame
    from tests.validation.test_birdcage_f_human_rung import LEG_COUNT, _ring_ports
    from tests.validation.test_wf10_f_human_tuned_circuit import _respec

    res, rec, comm = tuned["res"], tuned["rec"], tuned["comm"]
    cell_rel = res.cells / HOLE_CELL_RECORD - 1.0
    d_az = (res.az_feeds[1] - res.az_feeds[0]) % 360.0
    d_mirror = abs(res.az_mirror - rec["az_mirror"])
    d_feeds = max(abs(a - b) for a, b in zip(res.az_feeds, rec["az_feeds"]))
    if comm.rank == 0:
        print(f"[WF-16 s9a A0] cells {res.cells} vs record {HOLE_CELL_RECORD} (rel {cell_rel:.3e}, "
              f"band {CELL_COUNT_BAND}, ASSERTED), npz cells {tuned['npz_cells']} (== ASSERTED); "
              f"feeds {res.feeds} at az {res.az_feeds[0]:.12f} / {res.az_feeds[1]:.12f} deg, "
              f"d_az {d_az:.12f} (90 within {AZIMUTH_MATCH_DEG:g}, ASSERTED); az_mirror "
              f"{res.az_mirror:.12f} vs JSON {rec['az_mirror']:.12f} (|d| {d_mirror:.3e}, band "
              f"AZIMUTH_MATCH_DEG ASSERTED); feed azimuths vs JSON |d| {d_feeds:.3e}", flush=True)
    assert abs(cell_rel) <= CELL_COUNT_BAND, cell_rel
    assert res.cells == tuned["npz_cells"], (res.cells, tuned["npz_cells"])
    assert abs(d_az - 90.0) < AZIMUTH_MATCH_DEG, d_az
    assert d_mirror < AZIMUTH_MATCH_DEG, d_mirror
    assert d_feeds < AZIMUTH_MATCH_DEG, d_feeds
    assert list(res.feeds) == rec["feeds"]
    assert res.currents["order"] == rec["order"]

    # (U): the builder's spec construction (`_build_f_human`, test_wf7 module),
    # re-done on the Study's mesh, then `_respec(…, None, imp)` as `_pair` does.
    from fem_em_solver.ports.lumped import LumpedSheetPortSpec

    msh, cell_tags = res.mesh, res.cell_tags
    ring = _ring_ports()
    tags_f = _interface_facet_tags(
        msh, cell_tags, {SHEET_IFACE + i: (PORT_LOWER + i, PORT_UPPER + i) for i in ring})
    chord = float(res.specs[0].gap_height_m)
    base = []
    for i in ring:
        area = _interface_area_or_zero(msh, tags_f, SHEET_IFACE + i, comm)
        phi_hat, _c = _ring_gap_frame(i, LEG_COUNT)
        base.append(LumpedSheetPortSpec(
            port_id=f"P{i}", facet_tag=int(SHEET_IFACE + i), port_impedance_ohm=50.0,
            gap_height_m=chord, sheet_width_m=area / chord,
            drive_direction=tuple(float(c) for c in phi_hat), drive_voltage_v=1.0 + 0.0j,
            interior=True))
    z_rc = _z_rc(_c_f(), FREQUENCY_HZ)
    z_th = _thevenin_feed_impedance(z_rc, 50.0)
    imp = {sp.port_id: (z_th if sp.port_id in res.feeds else z_rc) for sp in base}
    ref = {sp.port_id: sp for sp in _respec(base, None, imp)}
    got = {sp.port_id: sp for sp in res.specs}
    assert sorted(got) == sorted(ref) and len(got) == 32
    worst_w, worst_h, other = 0.0, 0.0, []
    for pid, r in ref.items():
        g = got[pid]
        worst_w = max(worst_w, abs(g.sheet_width_m / r.sheet_width_m - 1.0))
        worst_h = max(worst_h, abs(g.gap_height_m / (res.sheets[0]["h"]) - 1.0))
        for f in dataclasses.fields(r):
            if f.name in ("sheet_width_m", "gap_height_m"):
                continue
            if getattr(g, f.name) != getattr(r, f.name):
                other.append((pid, f.name, getattr(g, f.name), getattr(r, f.name)))
    if comm.rank == 0:
        print(f"[WF-16 s9a U] 32 specs: worst w rel {worst_w:.3e}, h rel {worst_h:.3e} (band "
              f"IDENTITY_BAND {IDENTITY_BAND:g}, ASSERTED); other-field mismatches "
              f"{other or 'none'} (== ASSERTED)", flush=True)
    assert worst_w <= IDENTITY_BAND and worst_h <= IDENTITY_BAND, (worst_w, worst_h)
    assert not other, other


@complex_only
def test_wf16_s9a_node_reproduces_h9p(tuned):
    """(W), (S_rep), (Y); the H17-node control printed."""
    res, rec, comm = tuned["res"], tuned["rec"], tuned["comm"]
    say = (lambda m: print(f"[WF-16 s9a] {m}", flush=True)) if comm.rank == 0 else (lambda m: None)
    say(f"node line: Z_RC {res.gap_impedance_ohm:.9e} Z_th {res.feed_impedance_ohm:.9e} V_s "
        f"{res.source_voltage_v:.9e}; record {NODE_LINE_RECORD}")
    d_w = {}
    for sense in ("ccw", "cw"):
        w_json = np.array([_cx(p) for p in rec["w"][sense]])
        d_w[sense] = float(np.max(np.abs(np.asarray(res.weights[sense]) - w_json)))
    d_ccw_senses = float(np.max(np.abs(np.asarray(res.weights["ccw"])
                                       - np.asarray(SENSES["(1,-1j)"]))))
    say(f"(W) weights {res.weights}; vs JSON ccw {d_w['ccw']:.3e} cw {d_w['cw']:.3e}; ccw vs "
        f"SENSES['(1,-1j)'] {d_ccw_senses:.3e} (band Y0_CODING_BAND {Y0_CODING_BAND:g}, ASSERTED)")
    s_json = np.array([[_cx(v) for v in row] for row in rec["s_node"]], dtype=np.complex128)
    r_s = _residual(res.s_node, s_json)
    say(f"(S_rep) s_node {np.array2string(res.s_node, precision=12)}; residual vs node JSON "
        f"{r_s:.6e} (band {REUSE_REPRODUCTION_RTOL:g}, ASSERTED; prior {S_REP_PRIOR:.6e})")
    for i in range(2):
        for j in range(2):
            say(f"(S_rep) S[{i},{j}] study {res.s_node[i, j]:.12e} json {s_json[i, j]:.12e} "
                f"|d| {abs(res.s_node[i, j] - s_json[i, j]):.3e}")
    r_y = {}
    for sense in ("ccw", "cw"):
        y_json = np.array([_cx(p) for p in rec["y_in"][sense]])
        y = np.asarray(res.y_in[sense])
        r_y[sense] = float(np.linalg.norm(y - y_json) / np.linalg.norm(y_json))
        y_ref = _active_y_in(res.s_node, res.weights[sense], 50.0)
        assert max(abs(a - b) for a, b in zip(y, y_ref)) <= Y0_CODING_BAND * max(abs(b) for b in y_ref)
        say(f"(Y) [{sense}] Y_in P17 {y[0]:.10e} P21 {y[1]:.10e} S; JSON {y_json[0]:.10e} / "
            f"{y_json[1]:.10e}; rel {r_y[sense]:.3e} (band {REUSE_REPRODUCTION_RTOL:g}, ASSERTED)")
    h17 = tuned["h17"]
    if h17 is not None and "s_node" in h17:
        try:
            s17 = np.array([[_cx(v) for v in row] for row in h17["s_node"]], dtype=np.complex128)
            say(f"control (PREDICTED miss > {REUSE_REPRODUCTION_RTOL:g}, printed only): s_node vs "
                f"H17 node JSON at {h17.get('c_f')!r} F residual {_residual(res.s_node, s17):.6e}")
        except Exception as exc:  # printed only, never asserted
            say(f"control: H17 node JSON s_node not comparable ({exc!r})")
    else:
        say("control: H17 node JSON absent or without s_node (printed only)")
    assert max(d_w.values()) <= Y0_CODING_BAND and d_ccw_senses <= Y0_CODING_BAND, (d_w,
                                                                                  d_ccw_senses)
    assert r_s <= REUSE_REPRODUCTION_RTOL, f"(S_rep) {r_s:.6e}"
    for sense, r in r_y.items():
        assert r <= REUSE_REPRODUCTION_RTOL, (sense, r)


@complex_only
def test_wf16_s9a_power_partition(tuned):
    """(P_w) per sense; the P_surf-blind control (ccw) asserted."""
    import ufl
    from dolfinx import fem

    res, rec, comm = tuned["res"], tuned["rec"], tuned["comm"]
    msh = res.mesh
    omega = 2.0 * np.pi * FREQUENCY_HZ
    z_s = LeontovichWall(COPPER_SIGMA).surface_impedance_ohm(omega)
    acct = {"mesh": msh, "facet_tags": res.facet_tags}
    dx = ufl.Measure("dx", domain=msh, subdomain_data=res.cell_tags)
    specs_by = {sp.port_id: sp for sp in res.specs}
    order = res.currents["order"]

    def _vol(sigma_f, e, measure):
        form = fem.form(0.5 * sigma_f * ufl.inner(e, e) * measure)
        return float(np.real(comm.allreduce(fem.assemble_scalar(form), op=MPI.SUM)))

    out = {}
    for sense in ("ccw", "cw"):
        f = res.fields[sense]
        ef, sigma = f.e_complex, f.sigma_field
        wf = dict(zip(res.feeds, res.weights[sense]))
        sheets = {pid: (dataclasses.replace(specs_by[pid].sheet(driven=True),
                                            source_voltage_v=complex(wf[pid]))
                        if pid in res.feeds else specs_by[pid].sheet(driven=False))
                  for pid in order}
        p_src = float(sum(_source_power_w(acct, sheets[pid], ef, omega) for pid in res.feeds))
        p_sheets = float(sum(_sheet_field_dissipation_w(acct, sheets[pid], ef, omega)
                             for pid in order))
        p_total = _vol(sigma, ef, dx)
        p_ph = _vol(sigma, ef, dx(PHANTOM_CELL_TAG))
        p_rest = p_total - p_ph
        p_surf = _surface_loss_w(msh, res.conductor_facet_tags, z_s, ef, comm)
        rel = abs(p_src - p_sheets - p_total - p_surf) / abs(p_src)
        blind = abs(p_src - p_sheets - p_total) / abs(p_src)
        jp = rec["partitions"][sense]
        r_src = abs(p_src / jp["p_src"] - 1.0)
        r_ph = abs(p_ph / jp["p_phantom"] - 1.0)
        out[sense] = (rel, blind, r_src, r_ph)
        if comm.rank == 0:
            print(f"[WF-16 s9a P_w] [{sense}] P_src {p_src:.9e} sumP_sheet {p_sheets:.9e} "
                  f"P_phantom {p_ph:.9e} P_rest {p_rest:.9e} P_surf {p_surf:.9e} W; rel "
                  f"{rel:.3e} (band {DISCRETE_IDENTITY_RTOL:g}, ASSERTED; prior "
                  f"{P_W_PRIOR[sense]:.3e}); P_surf-blind miss {blind:.4e} = "
                  f"{blind / DISCRETE_IDENTITY_RTOL:.1f}x (JSON {jp['blind']:.4e}; bar "
                  f"{BLIND_FACTOR:.0f}x, {'ASSERTED' if sense == 'ccw' else 'printed'}); vs JSON "
                  f"P_src rel {r_src:.3e}, P_phantom rel {r_ph:.3e} (band "
                  f"{REUSE_REPRODUCTION_RTOL:g}, ASSERTED)", flush=True)
    for sense, (rel, _blind, r_src, r_ph) in out.items():
        assert rel <= DISCRETE_IDENTITY_RTOL, (sense, rel)
        assert r_src <= REUSE_REPRODUCTION_RTOL, (sense, r_src)
        assert r_ph <= REUSE_REPRODUCTION_RTOL, (sense, r_ph)
    assert out["ccw"][1] >= BLIND_FACTOR * DISCRETE_IDENTITY_RTOL, out["ccw"][1]
