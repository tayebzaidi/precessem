"""GEO-38 C31 (§9 item 821, §10 chain C): P33-3T's body — P61's module
(``test_geo38_f_clinical_box_heating.py``) composed with C21's 3 T patch set,
stub-run before ``xl`` entry 60 (no solve, no gmsh).

Brief: ``docs/planning/chunks/GEO-38.md`` "§9 item 821".  Composition, never an
edit: P61's `_run_heat_box` (C6a's `_run_box`, then the gel twin's thermal arm on
each sense, the dose at C6's c*, per accepted watt) runs unchanged under C21's
`_c21_patches` (C20's mesh route, C20's records, 128 MHz on ``wf7`` and ``c4``,
the lattice-mean print).  P61, C21, C6a and ``src`` stay byte-identical.

Selected by ``FEM_EM_GEO38_C31`` only: unset ⇒ the gated test skips and the coding
tests run (no gmsh); ``1`` ⇒ the body, with the node knob
``FEM_EM_GEO36_S8_DRIVE_C_PF`` required (I36b's `_node_knob`; C16's 3 T node is
44.0615 pF, never C6a's 193.5996); anything else ⇒ ``ValueError``.  The senses are
P61b's knob ``FEM_EM_GEO38_C6_HEAT_SENSE`` (`_sense_knob`, unchanged).  Degree 1
by default (I32a's `_refuse_degree` still refuses the rod's S8 knob); P87 (§9 item
969): C52's ``FEM_EM_GEO38_C21_DEGREE=2`` routes degree 2 through `_c31_degree`
to P61's `_run_heat_box(degree=2)`, and the 1.5 T knob ``FEM_EM_GEO38_C6_DEGREE``
is refused here.  The foreign guard is P61's list plus P61's own
knob, C21's knob and C21's extras.

**The 3 T gel σ / ε_r and the thermal constants are 1.5 T's** — provisional until
W3; the body prints that and claims nothing about a 3 T gel.

Scope: the body and its stub proof.  Not claimed: any 128 MHz temperature (that is
entry 60), a 3 T gel value, or the 1.5 T / 3 T ratio.
"""

from __future__ import annotations

import os
import sys

import numpy as np
import pytest
from mpi4py import MPI

from tests.complex_mode import complex_only
from tests.validation import test_geo36_rod_f_clinical_heating as i32
from tests.validation import test_geo38_f_clinical_3t_box_loaded as c21
from tests.validation import test_geo38_f_clinical_box_heating as p61
from tests.validation.test_mat10_sar_to_heat_sphere import B_BAND

C31_ENV = "FEM_EM_GEO38_C31"
#: Variables a C31 window never reads: P61's list, P61's own knob, C21's knob and
#: C21's extras (order kept, duplicates dropped).
FOREIGN_ENVS = tuple(
    dict.fromkeys(
        tuple(p61.FOREIGN_ENVS) + (p61.HEAT_ENV, c21.C21_ENV) + tuple(c21.FOREIGN_ENVS_C21_EXTRA)
    )
)


def _say(comm, msg):
    if comm.rank == 0:
        print(f"[GEO-38 C31] {msg}", flush=True)


def _c31_mode():
    raw = os.environ.get(C31_ENV, "").strip()
    if raw == "":
        return None
    if raw == "1":
        return raw
    raise ValueError(f"{C31_ENV} must be unset or '1', got {raw!r}")


def _knobs():
    """Foreign guard, degree refusal, then ``(mode, c_pf)`` (I36b's node knob)."""
    from tests.validation import test_geo36_rod_f_clinical_drive as i36

    bad = [k for k in FOREIGN_ENVS if os.environ.get(k) is not None]
    if bad:
        raise ValueError(f"C31 never reads these variables; unset them: {bad}")
    mode = _c31_mode()
    i32._refuse_degree()
    return mode, i36._node_knob(mode)


def _c31_degree(mode):
    """P87 (§9 item 969): C52's ``C21_DEGREE_ENV`` through `c21._c21_degree` (unset /
    empty ⇒ ``None``, ``"2"`` ⇒ 2, anything else ⇒ ValueError); the 1.5 T knob
    ``FEM_EM_GEO38_C6_DEGREE`` refused."""
    if os.environ.get(p61.c6.C6_DEGREE_ENV) is not None:
        raise ValueError(f"C31 never routes a 1.5 T degree; unset {p61.c6.C6_DEGREE_ENV}")
    return c21._c21_degree(mode)


def _run_heat_box_3t(comm, monkeypatch, c_pf, senses, *, degree=None):
    """Entry 60's body: C21's patch set, then P61's `_run_heat_box` unchanged.
    (A0) / U₁ assert C20's 1 128 430 / 1 508 256 inside C4a's fixture (the records
    `_c21_patches` sets, read by `i36._route_drive` at call time)."""
    m = c21._modules()
    c21._c21_patches(monkeypatch)
    _say(
        comm,
        f"3 T patch set: rung Ø {m['t3'].RUNG_DIAMETER_3T_M} m, f {m['wf7'].FREQUENCY_HZ:.6e} Hz "
        f"(C4 {m['c4'].FREQUENCY_HZ_C4:.6e}); records cells "
        f"{m['i35'].F_CLINICAL_ROD_CELL_RECORD}, U {m['i35'].F_CLINICAL_ROD_UNKNOWNS_RECORD}; "
        f"node {c_pf:.4f} pF; senses {list(senses)}" + (f"; degree {degree}" if degree else ""),
    )
    _say(
        comm,
        "the gel sigma / eps_r and the thermal constants are the 1.5 T values (provisional "
        "until W3; no 3 T gel claim)",
    )
    return p61._run_heat_box(comm, monkeypatch, c_pf, senses=senses, degree=degree)


@complex_only
def test_geo38_c31_box_heating_3t(monkeypatch):
    mode, c_pf = _knobs()  # node knob missing under 1 => ValueError, before gmsh
    senses = p61._sense_knob()  # a bad sense raises before the mode skip
    deg = _c31_degree(mode)  # P87: a bad degree / the 1.5 T knob raises before the mode skip
    if mode is None:
        pytest.skip(f"{C31_ENV} unset: GEO-38 C31 is selected by environment only")
    _run_heat_box_3t(MPI.COMM_WORLD, monkeypatch, c_pf, senses, degree=deg)


# ---------------------------------------------------------------------------
# Coding tests (no gmsh); skip whenever C31_ENV is set.
# ---------------------------------------------------------------------------
_coding = pytest.mark.skipif(
    os.environ.get(C31_ENV, "").strip() != "", reason=f"{C31_ENV} set: coding tests skip in windows"
)


@_coding
def test_c31_k0_knobs_foreign_sense_and_box_time_base(monkeypatch):
    """(K0-3T): modes unset / 1 accepted, others refused; the node required under 1;
    degree 2 refused; every foreign variable refused; the sense knob P61b's; the
    box time base (`i32._box_preconditions`)."""
    from tests.validation.test_geo36_rod_f_clinical_drive import (
        S8_DRIVE_C_PF_ENV,
        S8_DRIVE_DEGREE_ENV,
    )

    for k in (C31_ENV, S8_DRIVE_C_PF_ENV, S8_DRIVE_DEGREE_ENV, p61.HEAT_SENSE_ENV) + FOREIGN_ENVS:
        monkeypatch.delenv(k, raising=False)
    assert _knobs() == (None, None)
    monkeypatch.setenv(C31_ENV, "1")
    with pytest.raises(ValueError, match=S8_DRIVE_C_PF_ENV):
        _knobs()
    monkeypatch.setenv(S8_DRIVE_C_PF_ENV, f"{c21.K_NODE_3T_PF}")
    assert _knobs() == ("1", c21.K_NODE_3T_PF)
    for bad in ("0", "build", "x"):
        monkeypatch.setenv(C31_ENV, bad)
        with pytest.raises(ValueError):
            _c31_mode()
    monkeypatch.setenv(C31_ENV, "1")
    monkeypatch.setenv(S8_DRIVE_DEGREE_ENV, "2")
    # The degree variable is in the foreign list (P61's FOREIGN_ENVS via C6a), so
    # `_knobs()` refuses it through the foreign guard first; `_refuse_degree` is
    # asserted on its own, as P61's K0 does (review ruling 2026-11-05 10:30, pt 1).
    with pytest.raises(ValueError, match="never reads"):
        _knobs()
    with pytest.raises(ValueError, match="refused"):
        i32._refuse_degree()
    monkeypatch.delenv(S8_DRIVE_DEGREE_ENV)
    for k in FOREIGN_ENVS:
        with monkeypatch.context() as m:
            m.setenv(k, "1")
            with pytest.raises(ValueError, match="C31 never reads"):
                _knobs()
    required = (
        {p61.HEAT_ENV, c21.C21_ENV, "FEM_EM_GEO38_C6"}
        | set(p61.FOREIGN_ENVS)
        | set(c21.FOREIGN_ENVS_C21_EXTRA)
    )
    assert set(FOREIGN_ENVS) == required and len(FOREIGN_ENVS) == len(set(FOREIGN_ENVS))
    assert p61.HEAT_SENSE_ENV not in FOREIGN_ENVS
    monkeypatch.setenv(p61.HEAT_SENSE_ENV, "both")
    assert p61._sense_knob() == ("ccw", "cw")
    monkeypatch.setenv(p61.HEAT_SENSE_ENV, "x")
    with pytest.raises(ValueError):
        test_geo38_c31_box_heating_3t(monkeypatch)  # a bad sense raises before the skip / body
    monkeypatch.delenv(p61.HEAT_SENSE_ENV)
    i32._box_preconditions()
    print(
        f"\n[GEO-38 C31 K0] modes unset/1 accepted, 0/build/x refused; node required under 1 "
        f"({c21.K_NODE_3T_PF}); degree 2 refused; {len(FOREIGN_ENVS)} foreign variables refused "
        f"(P61's {len(p61.FOREIGN_ENVS)} + {p61.HEAT_ENV} + {c21.C21_ENV} + C21's extras, deduplicated); "
        f"sense knob P61b's (both -> ccw, cw; x refused); box time base (ASSERTED)",
        flush=True,
    )


@_coding
def test_c31_k10p_point_map_under_the_3t_patches(monkeypatch):
    """(K10-p-3T): under `_c21_patches`, P37's six points equal P61's mapping (the
    same box and ``TABLE_OFFSET_M`` at both field strengths)."""
    pts15 = p61._probes()
    with monkeypatch.context() as m:
        c21._c21_patches(m)
        pts3 = p61._probes()
        box3 = c21._modules()["box"]._box()
    box15 = c21._modules()["box"]._box()
    same_box = all(np.array_equal(a, b) for a, b in zip(box15, box3))
    print(
        f"\n[GEO-38 C31 K10-p] P37's points 1.5 T {pts15}; under the 3 T patches {pts3}; box "
        f"centre / half equal {same_box} (ASSERTED)",
        flush=True,
    )
    assert pts3 == pts15, (pts15, pts3)
    assert same_box, (box15, box3)


@_coding
def test_c31_k7_one_with_node_reaches_builder_with_3t_kwargs(monkeypatch):
    """(K7-3T): ``1`` + node 44.0615 runs the gated test into C21's
    `_builder_sentinel` with rung Ø 0.704, f 128e6 on ``wf7`` and ``c4``, and C20's
    records 1 128 430 / (1 508 256, 7 783 028).  Negative control (asserted,
    categorical): P61's `_run_heat_box` without C21's patches reaches rung Ø 0.710,
    64 MHz and 1 217 025 cells — the name-binding trap
    (`20261103T103420Z_GEO-38-C21-W0.log:83`)."""
    from tests.validation.test_geo36_rod_f_clinical_drive import S8_DRIVE_C_PF_ENV

    mods = c21._modules()
    t3, c20, i35 = mods["t3"], mods["c20"], mods["i35"]
    kw15 = i35._kwargs()
    comm = MPI.COMM_WORLD

    seen = []
    with monkeypatch.context() as m:
        for k in FOREIGN_ENVS + (p61.HEAT_SENSE_ENV,):
            m.delenv(k, raising=False)
        m.setenv(C31_ENV, "1")
        m.setenv(S8_DRIVE_C_PF_ENV, f"{c21.K_NODE_3T_PF}")
        c21._builder_sentinel(m, seen)
        with pytest.raises(c21._BuilderReached):
            test_geo38_c31_box_heating_3t(m)

    neg = []
    with monkeypatch.context() as m:
        for k in FOREIGN_ENVS + (p61.HEAT_SENSE_ENV,):
            m.delenv(k, raising=False)
        c21._builder_sentinel(m, neg)
        with pytest.raises(c21._BuilderReached):
            p61._run_heat_box(comm, m, c21.K_NODE_3T_PF)

    assert len(seen) == 1 and len(neg) == 1, (len(seen), len(neg))
    got, got_neg = seen[0], neg[0]
    rung_d, rung_d_neg = 2.0 * got["kw"]["ring_radius"], 2.0 * got_neg["kw"]["ring_radius"]
    if comm.rank == 0:
        print(
            f"\n[GEO-38 C31 K7] C31=1 node {c21.K_NODE_3T_PF}: builder reached with rung Ø {rung_d:.12f} m "
            f"(3 T {t3.RUNG_DIAMETER_3T_M}); f wf7 {got['f_wf7']:.6e} / C4 {got['f_c4']:.6e} Hz; records "
            f"cells {got['cells']}, U {got['unknowns']} (ASSERTED)"
        )
        print(
            f"[GEO-38 C31 K7] control (P61 without C21's patches): rung Ø {rung_d_neg:.12f} m, f wf7 "
            f"{got_neg['f_wf7']:.6e} / C4 {got_neg['f_c4']:.6e} Hz, records cells {got_neg['cells']} "
            f"(ASSERTED: the 1.5 T fixture)",
            flush=True,
        )
    assert (
        abs(rung_d / t3.RUNG_DIAMETER_3T_M - 1.0) <= 1e-15 and abs(rung_d - 0.704) <= 1e-12
    ), rung_d
    assert got["f_wf7"] == t3.FREQUENCY_3T_HZ == 128e6 and got["f_c4"] == t3.FREQUENCY_3T_HZ
    assert got["cells"] == c20.F_CLINICAL_3T_ROD_CELL_RECORD == 1_128_430
    assert got["unknowns"] == c20.F_CLINICAL_3T_ROD_UNKNOWNS_RECORD == (1_508_256, 7_783_028)
    # Negative control: the name-binding trap, categorical.
    assert rung_d_neg == 2.0 * kw15["ring_radius"] and abs(rung_d_neg - 0.710) <= 1e-12, rung_d_neg
    assert got_neg["f_wf7"] == 64.0e6 and got_neg["f_c4"] == 64.0e6, got_neg
    assert got_neg["cells"] == 1_217_025, got_neg["cells"]
    # Teardown restored the 1.5 T namespace.
    assert i35.F_CLINICAL_ROD_CELL_RECORD == 1_217_025
    assert mods["s2"]._box_kwargs()["ring_radius"] == kw15["ring_radius"]


def _stub_3t(m, fx, json_dir):
    """C21's `_k8_stubs_3t` plus P61's K10-R stub settings (every face Robin, the stub
    probes).  Returns C6a's recorded lines and P61's module."""
    lines, _seen, _z, _zrc = c21._k8_stubs_3t(m, fx, json_dir)
    this = sys.modules[p61._run_heat_box.__module__]
    m.setattr(i32, "WALL_SELECT_M", -1.0)  # no rod wall on K6's box: every face is Robin
    m.setattr(this, "_probes", lambda: i32.STUB_PROBES)
    return lines, this


@complex_only
@_coding
def test_c31_k10r_run_heat_box_3t_on_stubs(monkeypatch, tmp_path):
    """(K10-R-3T): `_run_heat_box_3t` on C21's `_k8_stubs_3t` at 128 MHz and node
    44.0615 pF (column 2 = 0.5j × column 1, so ccw = 1.5 E, cw = 0.5 E): P_gel ccw/cw
    = 9; G·P_gel = 1; P_acc ccw against the 3 T JSON to 1e-12; the scale applied
    exactly; the thermal gates clean; the ``(C21) mean |B1+|`` heads present."""
    comm = MPI.COMM_WORLD
    box = c21._modules()["box"]
    fx = box._k6_fixture(comm)
    root = comm.bcast(str(tmp_path), root=0)
    with monkeypatch.context() as m:
        lines, _this = _stub_3t(m, fx, root)
        out = _run_heat_box_3t(comm, m, c21.K_NODE_3T_PF, ("cw",))
        f_seen = float(c21._modules()["wf7"].FREQUENCY_HZ)
    arm, dose = out["arm"], out["dose"]
    r9 = abs(out["p_gel_ccw"] / arm["P"] / 9.0 - 1.0)
    r_gp = abs(arm["G"] * arm["P"] - 1.0)
    r_scale = float(np.max(np.abs(dose["th_acc"] - out["scale"] * arm["hist"])))
    rec = c21._k8_json_3t(comm, root)
    r_pacc = abs(out["p_acc"]["ccw"] / rec["p_acc"] - 1.0)
    heads = (
        "(i) ball radius",
        "normalisers:",
        "price:",
        "(C21) mean |B1+|_ccw",
        "(C21) mean |B1+|_cw",
    )
    missing = [h for h in heads if not any(ln.startswith(h) for ln in lines)]
    fails = i32._thermal_failures(arm)
    if comm.rank == 0:
        print(
            f"\n[GEO-38 C31 K10-R] f {f_seen:.6e} Hz, node {rec['c_pf']:.4f} pF; P_gel ccw/cw miss from 9 "
            f"{r9:.3e}; G*P_gel miss {r_gp:.3e} (band {B_BAND:.0e}); P_acc ccw vs the 3 T JSON {r_pacc:.3e} "
            f"(<= 1e-12); scale {out['scale']:.9e}, applied miss {r_scale:.3e}; cem43 {dose['cem43']:.6e}, "
            f"running max {dose['ra_max']:.6e} K; thermal failures {fails or 'none'}; heads missing "
            f"{missing or 'none'} (ASSERTED)",
            flush=True,
        )
    assert f_seen == 128e6, f_seen
    assert rec["c_pf"] == c21.K_NODE_3T_PF, rec["c_pf"]
    assert r9 <= B_BAND, r9
    assert r_gp <= B_BAND, r_gp
    assert r_pacc <= 1e-12, r_pacc
    assert r_scale == 0.0, r_scale
    assert out["merged_rod"] == 0
    assert arm["counts"]["wall"] == 0 and all(arm["probe_ok"]) and all(arm["hist_ok"])
    assert not fails, fails
    assert np.isfinite(dose["cem43"]) and dose["cem43"] > 0.0
    assert not missing, missing


@complex_only
@_coding
def test_c31_k10r2_both_senses_3t_on_stubs(monkeypatch, tmp_path):
    """(K10-R2-3T): P61b's (R2-p) / (R2-h) / (R2-s) / (R2-o) / (R2-g) at 128 MHz with
    ``senses=("ccw", "cw")``.  Control (R2-c-3T): the ccw arm normalised by P_acc(cw)
    misses (R2-s) by |P_acc(ccw)/P_acc(cw) − 1|, asserted > ``B_BAND`` at factor 1
    (no log backs a 128 MHz figure; the value is printed)."""
    comm = MPI.COMM_WORLD
    box = c21._modules()["box"]
    fx = box._k6_fixture(comm)
    root = comm.bcast(str(tmp_path), root=0)
    said = []
    with monkeypatch.context() as m:
        _lines, this = _stub_3t(m, fx, root)
        real_say = this._say

        def _rec(c, msg):
            said.append(msg)
            real_say(c, msg)

        m.setattr(this, "_say", _rec)
        out = _run_heat_box_3t(comm, m, c21.K_NODE_3T_PF, ("ccw", "cw"))
    bs, p_acc = out["by_sense"], out["p_acc"]
    a_cc, a_cw = bs["ccw"]["arm"], bs["cw"]["arm"]
    d_cc, d_cw = bs["ccw"]["dose"], bs["cw"]["dose"]
    r_p = abs(a_cc["P"] / a_cw["P"] / 9.0 - 1.0)
    r_h = float(np.max(np.abs(a_cc["hist"] - a_cw["hist"])) / np.max(np.abs(a_cw["hist"])))
    applied = {
        s: float(np.max(np.abs(bs[s]["dose"]["th_acc"] - bs[s]["scale"] * bs[s]["arm"]["hist"])))
        for s in ("ccw", "cw")
    }
    scale_def = {s: abs(bs[s]["scale"] - bs[s]["arm"]["P"] / p_acc[s]) for s in ("ccw", "cw")}
    predicted = 9.0 * p_acc["cw"] / p_acc["ccw"]
    r_s = abs(d_cc["th_acc"][-1] / d_cw["th_acc"][-1] / predicted - 1.0)
    wrong = a_cc["hist"][-1] * a_cc["P"] / p_acc["cw"]  # control: ccw arm on cw's normaliser
    miss_c = abs(wrong / d_cw["th_acc"][-1] / predicted - 1.0)
    pred_c = abs(p_acc["ccw"] / p_acc["cw"] - 1.0)
    thermal = [i for i, ln in enumerate(said) if ln.startswith(("hotspot ", "probe "))]
    first = {t: next((i for i in thermal if t in said[i]), None) for t in ("(1,-1j)", "(1,+1j)")}
    order_ok = (
        first["(1,-1j)"] is not None
        and first["(1,+1j)"] is not None
        and first["(1,-1j)"] < first["(1,+1j)"]
    )
    g = {s: i32._thermal_failures(bs[s]["arm"]) for s in ("ccw", "cw")}
    if comm.rank == 0:
        print(
            f"\n[GEO-38 C31 K10-R2] (R2-p) P_gel ccw {a_cc['P']:.12e} / cw {a_cw['P']:.12e} W miss from 9 "
            f"{r_p:.3e}; (R2-h) max|h_ccw - h_cw|/max|h_cw| {r_h:.3e}; (R2-s) applied {applied}, scale "
            f"defs {scale_def}, theta_acc(t*) ccw {d_cc['th_acc'][-1]:.9e} / cw {d_cw['th_acc'][-1]:.9e} K "
            f"vs 9 P_acc(cw)/P_acc(ccw) {predicted:.12e}: miss {r_s:.3e} (band {B_BAND:.0e}); (R2-o) first "
            f"thermal line (1,-1j) {first['(1,-1j)']} < (1,+1j) {first['(1,+1j)']}: {order_ok}; (R2-g) "
            f"failures {g} (ASSERTED)"
        )
        print(
            f"[GEO-38 C31 K10-R2] control (R2-c-3T) ccw on P_acc(cw) at 128 MHz: miss {miss_c:.6e} vs "
            f"|P_acc(ccw)/P_acc(cw) - 1| {pred_c:.6e} (1.5 T stub: {p61.SENSE_NORMALISER_MISS:.4e}, "
            f"printed); {miss_c / B_BAND:.4e} x B_BAND (> B_BAND, ASSERTED)",
            flush=True,
        )
    assert r_p <= B_BAND, r_p
    assert r_h <= B_BAND, r_h
    assert all(v == 0.0 for v in applied.values()), applied
    assert all(v == 0.0 for v in scale_def.values()), scale_def
    assert r_s <= B_BAND, r_s
    assert order_ok, first
    assert not g["ccw"] and not g["cw"], g
    assert all(a_cc["hist_ok"]) and all(a_cw["hist_ok"])
    assert out["arm"] is a_cw
    assert miss_c > B_BAND, miss_c


# ---------------------------------------------------------------------------
# P87 (§9 item 969): C52's 3 T degree knob routed through C31.
# ---------------------------------------------------------------------------
#: (Deg-arith): C20's build, E and F (`20261101T181021Z_GEO-38-C20-W1-build.log:11929`); U₂
#: measured at C52's W1-build (`20261112T180748Z_GEO-38-C52-W1-build.log:11936`).
E_3T, F_3T = 1_508_256, 2_383_258
#: (N-patch-3T): |8327394 / 7783028 - 1| = 6.994e-02 (review arithmetic); bar.
N_PATCH_BAR = 6.0e-02


@_coding
def test_p87_kdeg_c31_knob_parsing(monkeypatch):
    """(K-deg-C31): C52's semantics (unset / '' / 2 accepted; 1 / x / 3 refused); the
    1.5 T knob ``FEM_EM_GEO38_C6_DEGREE=2`` refused, naming it."""
    c6_deg = p61.c6.C6_DEGREE_ENV
    monkeypatch.delenv(c6_deg, raising=False)
    monkeypatch.delenv(c21.C21_DEGREE_ENV, raising=False)
    got = {"(unset)": _c31_degree("1")}
    for raw in ("", "2"):
        monkeypatch.setenv(c21.C21_DEGREE_ENV, raw)
        got[repr(raw)] = _c31_degree("1")
    refused = []
    for bad in ("1", "x", "3"):
        monkeypatch.setenv(c21.C21_DEGREE_ENV, bad)
        with pytest.raises(ValueError, match=c21.C21_DEGREE_ENV):
            _c31_degree("1")
        refused.append(bad)
    monkeypatch.delenv(c21.C21_DEGREE_ENV)
    monkeypatch.setenv(c6_deg, "2")
    with pytest.raises(ValueError, match=c6_deg):
        _c31_degree("1")
    print(
        f"\n[GEO-38 P87 K-deg-C31] {c21.C21_DEGREE_ENV}: {got}; refused {refused}; {c6_deg}=2 refused "
        f"(ASSERTED)",
        flush=True,
    )
    assert got == {"(unset)": None, "''": None, "'2'": 2}, got


@_coding
def test_p87_k70_entry70_environment_shape(monkeypatch):
    """(K-70): entry 70's environment (node 44.0525, a stub value) accepted alone."""
    from tests.validation.test_geo36_rod_f_clinical_drive import (
        S8_DRIVE_C_PF_ENV,
        S8_DRIVE_DEGREE_ENV,
    )

    for k in (
        C31_ENV,
        S8_DRIVE_C_PF_ENV,
        S8_DRIVE_DEGREE_ENV,
        p61.HEAT_SENSE_ENV,
        c21.C21_DEGREE_ENV,
        p61.c6.C6_DEGREE_ENV,
    ) + FOREIGN_ENVS:
        monkeypatch.delenv(k, raising=False)
    monkeypatch.setenv(C31_ENV, "1")
    monkeypatch.setenv(c21.C21_DEGREE_ENV, "2")
    monkeypatch.setenv(p61.HEAT_SENSE_ENV, "both")
    monkeypatch.setenv(S8_DRIVE_C_PF_ENV, "44.0525")
    got = (_knobs(), p61._sense_knob(), _c31_degree("1"))
    print(
        f"\n[GEO-38 P87 K-70] {C31_ENV}=1 {c21.C21_DEGREE_ENV}=2 {p61.HEAT_SENSE_ENV}=both "
        f"{S8_DRIVE_C_PF_ENV}=44.0525: {got} (ASSERTED)",
        flush=True,
    )
    assert got == (("1", 44.0525), ("ccw", "cw"), 2), got


@_coding
def test_p87_n_c6deg_refused_under_c31(monkeypatch):
    """(N-C6deg-C31), categorical: ``FEM_EM_GEO38_C6_DEGREE=2`` under C31 ``1`` (node
    set) is refused before gmsh — the builder sentinel is never reached."""
    from tests.validation.test_geo36_rod_f_clinical_drive import S8_DRIVE_C_PF_ENV

    seen = []
    with monkeypatch.context() as m:
        for k in FOREIGN_ENVS + (p61.HEAT_SENSE_ENV, c21.C21_DEGREE_ENV):
            m.delenv(k, raising=False)
        m.setenv(C31_ENV, "1")
        m.setenv(S8_DRIVE_C_PF_ENV, f"{c21.K_NODE_3T_PF}")
        m.setenv(p61.c6.C6_DEGREE_ENV, "2")
        c21._builder_sentinel(m, seen)
        with pytest.raises(ValueError, match=p61.c6.C6_DEGREE_ENV) as err:
            test_geo38_c31_box_heating_3t(m)
    print(
        f"\n[GEO-38 P87 N-C6deg-C31] C31=1 + {p61.c6.C6_DEGREE_ENV}=2: refused ({err.value}); builder "
        f"reached {len(seen)} (ASSERTED 0)",
        flush=True,
    )
    assert seen == [], seen


@complex_only
@_coding
def test_p87_k10r_3t_d2_on_degree2_stubs(monkeypatch, tmp_path):
    """(K10-R-3T-d2): `_run_heat_box_3t(…, degree=2)` on `_k6_fixture_d2` with C34's
    recorders under C21's patches — f 128e6, (Deg) 7783028 exact, (Deg-arith), the
    ``-d2`` JSON, (B) / (R2-*) / (V-acc) / (V-s), the (C21) heads.  Control
    (N-patch-3T): the same stubs through P61 without `_c21_patches` print 8327394."""
    import json

    mods = c21._modules()
    box, c20 = mods["box"], mods["c20"]
    comm = MPI.COMM_WORLD
    fx2 = p61._k6_fixture_d2(comm)
    root = comm.bcast(str(tmp_path), root=0)
    d2_dir, neg_dir = os.path.join(root, "d2"), os.path.join(root, "neg")
    said = []
    with monkeypatch.context() as m:
        m.delenv(box.C6_DEGREE_ENV, raising=False)
        m.delenv(mods["c4"].C4_DEGREE_ENV, raising=False)
        lines, this = _stub_3t(m, fx2, d2_dir)
        real_say = this._say

        def _rec(c, msg):
            said.append(msg)
            real_say(c, msg)

        m.setattr(this, "_say", _rec)
        env_seen, degrees = box._c34_recorders(m, degree=2)
        out = _run_heat_box_3t(comm, m, c21.K_NODE_3T_PF, ("ccw", "cw"), degree=2)
        f_seen = float(mods["wf7"].FREQUENCY_HZ)
    with monkeypatch.context() as m:
        m.delenv(box.C6_DEGREE_ENV, raising=False)
        lines_n, this_n = _stub_3t(m, fx2, neg_dir)
        box._c34_recorders(m, degree=2)
        this_n._run_heat_box(comm, m, c21.K_NODE_3T_PF, senses=("cw",), degree=2)
    rec = None
    if comm.rank == 0:
        with open(os.path.join(d2_dir, f"c6-node-{c21.K_NODE_3T_PF:.4f}-d2.json")) as fh:
            rec = json.load(fh)
    rec = comm.bcast(rec, root=0)
    u2 = c20.F_CLINICAL_3T_ROD_UNKNOWNS_RECORD[1]
    deg_lines = [ln for ln in lines if ln.startswith("(Deg)")]
    deg_want = f"(Deg) degree 2 global unknowns {u2} vs record {u2} (EXACT, ASSERTED)"
    arith = 2 * (E_3T + F_3T)
    r2 = p61._r2_checks(out, said)
    heads = ("(C21) mean |B1+|_ccw", "(C21) mean |B1+|_cw")
    missing = [h for h in heads if not any(ln.startswith(h) for ln in lines)]
    neg_deg = [ln for ln in lines_n if ln.startswith("(Deg)")]
    n_neg = int(neg_deg[0].split("global unknowns ")[1].split(" ")[0]) if neg_deg else None
    r_neg = abs(n_neg / u2 - 1.0) if n_neg else None
    if comm.rank == 0:
        print(
            f"\n[GEO-38 P87 K10-R-3T-d2] f {f_seen:.6e} Hz; element degree {p61._element_degree(fx2[2])}; C4 "
            f"variable at the fixture call {env_seen}; solver degrees {degrees}; "
            f"{deg_lines[0] if deg_lines else '(Deg) line MISSING'}; (Deg-arith) 2*({E_3T} + {F_3T}) = "
            f"{arith} vs C20's record {u2} (ASSERTED exact)"
        )
        print(
            f"[GEO-38 P87 K10-R-3T-d2] JSON -d2 keys == K8_JSON_KEYS {set(rec) == box.K8_JSON_KEYS}; (B) ccw "
            f"{r2['B']['ccw']:.3e} cw {r2['B']['cw']:.3e} (band {B_BAND:.0e}); failures {r2['g']}; (R2-p) "
            f"{r2['r_p']:.3e}, (R2-h) {r2['r_h']:.3e}, (R2-s) applied {r2['applied']} miss {r2['r_s']:.3e}, "
            f"(R2-o) {r2['order_ok']}; (V-acc) {r2['v_acc']}; (V-s) vol_mean ccw {r2['vm']['ccw']:.12e} / "
            f"cw {r2['vm']['cw']:.12e} miss {r2['r_vs']:.3e}; heads missing {missing or 'none'} (ASSERTED)"
        )
        print(
            f"[GEO-38 P87 N-patch-3T] without _c21_patches: {neg_deg[0] if neg_deg else 'MISSING'}; "
            f"|{n_neg}/{u2} - 1| {r_neg if r_neg is None else f'{r_neg:.3e}'} (> {N_PATCH_BAR:.0e}, ASSERTED)",
            flush=True,
        )
    assert f_seen == 128e6, f_seen
    assert p61._element_degree(fx2[2]) == 2
    assert env_seen == ["2"], env_seen
    assert degrees == [2], degrees
    assert out["degree"] == 2
    assert deg_lines == [deg_want], deg_lines
    assert arith == 7_783_028 == u2, (arith, u2)
    assert set(rec) == box.K8_JSON_KEYS, sorted(rec)
    p61._r2_asserts(r2)
    assert len(r2["vlines"]) == 2, r2["vlines"]
    assert not missing, missing
    assert n_neg == 8_327_394 and n_neg != u2, n_neg
    assert r_neg > N_PATCH_BAR, r_neg
