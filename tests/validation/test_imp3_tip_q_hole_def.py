"""`GEO-36` I67 (§9 item 825): the rod model's lever on one definition.

Entry 48's tip-Q / η definition (Q = SAR_1g(rod) / SAR_1g(gel twin) on the four
balls of ``_probe_points``, η_rod(T1) = |E_rod(T1)| / |E_gel(T1)|) run with the
**PEC-hole** rod model on I41's matched pair (review ruling 2026-11-05 03:00):

- **hole arm** — ``_rb._arm(_build_rod(), rod=True)``: the rod is a hole with
  the PEC wall on tag 402, ``ROD_CELL_RECORD_W1`` cells;
- **gel twin** — I57's ``_gel_seen`` on ``_build_kept_rod()``,
  ``KEPT_ROD_CELL_RECORD`` cells, tag 403 → gel.

The two share their tag-3 cells exactly (I41's (A0m)).  Arms are labelled by
content.  No landed module is edited; every gate and band is imported unmoved.

Gates: (A0) the two cell records exact; (A0m) owned tag-3 counts equal and the
kept build's tag-403 count = 15 000; per arm (i) ``RECIPROCITY_BAND``, (ii)
σ_max ≤ 1 + ``PASSIVITY_SIGMA_TOLERANCE`` and every column power ≤ 1, (P) all
cells at ``DISCRETE_IDENTITY_RTOL`` (the two-tag residual printed beside it;
on the hole arm the two are expected equal — a hole build has no lossy tag
outside 3 and 1, `OPS-80`); (E) on the hole arm ``wall_max`` ≤
``PEC_WALL_MAX_E`` for every drive, wall dofs > 0.

Controls (asserted): (C1) ``_control_c1({"rod": hole, "off": gel})`` reads
max|E| on the 402 dofs > ``PEC_WALL_MAX_E`` with its (P) ≤
``DISCRETE_IDENTITY_RTOL`` (`IMP-3` W1 measured 1.29e-4 on this build and drive,
`20261002T170428Z_IMP-3-step1-W1.log`; C1's "natural/without" predicted ranges
were written for the no-rod mesh and are mislabelled here — printed only).
(Rep) the apex Q reproduces I41's 4.090342 / 4.113695
(`20261025T201613Z_IMP-3-I41-W1b.log:3776–3777`) within
``DISCRETE_IDENTITY_RTOL`` relative; the ``.6f`` print's half-ulp is ≈ 1.2e-7
relative, so the bound is not the print's.

Record (no band): Q on the four balls beside I57's kept-cell values, the
hole / kept ratio per ball, η_rod(T1) beside I57's 2.157532, `IMP-3`'s
interval and entry 48's figures; the ball masses on both arms.

Scope: one F-small pair, degree 1, 64 MHz, the ccw drive.  Not claimed:
convergence, which rod model is right, anything clinical.
"""

from __future__ import annotations

import gc
import os
import time

import numpy as np
import pytest
from mpi4py import MPI

import tests.validation.test_imp3_rod_birdcage as _rb
from tests.complex_mode import complex_only
from tests.mesh.test_birdcage_implant_rod import (
    KEPT_ROD_CELL_RECORD,
    ROD_CENTRE,
    ROD_LENGTH,
    ROD_RADIUS,
    _build_kept_rod,
    _build_rod,
)
from tests.validation.test_birdcage_power_identity import DISCRETE_IDENTITY_RTOL
from tests.validation.test_geo36_rod_f_clinical_drive import _probe_points
from tests.validation.test_imp3_rod_birdcage import (
    PEC_WALL_MAX_E,
    PORT_IDS,
    ROD_CELL_RECORD_W1,
    TIPS,
    _arm,
    _control_c1,
)
from tests.validation.test_imp3_rod_heating_matched import _census
from tests.validation.test_imp3_rod_heating_phantom import SOURCE
from tests.validation.test_imp3_tip_ladder_1g import ENTRY43_READOUT_INTERVAL
from tests.validation.test_imp3_tip_q_entry48_def import (
    BALLS,
    ENTRY48_ETA_T1,
    ENTRY48_Q,
    _gate_failures,
    _gel_seen,
    _make_hook,
    _print_gates,
    _two_tag_resid,
)

KNOB = "FEM_EM_IMP3_I67"
ROD_CELLS = KEPT_ROD_CELL_RECORD - ROD_CELL_RECORD_W1  # 15 000
# I41's apex Q, hole / gel, ccw (`20261025T201613Z_IMP-3-I41-W1b.log:3776–3777`).
I41_APEX_Q = {"top_apex": 4.090342, "bottom_apex": 4.113695}
# I57's kept-cell Q and η (`20261104T104820Z_GEO-36-I57-W1d.log:1999–2003`).
I57_KEPT_Q = {"top": 3.951922, "bottom": 3.973311, "top_apex": 3.899322, "bottom_apex": 3.921665}
I57_KEPT_ETA_T1 = 2.157532


def _say(msg):
    if MPI.COMM_WORLD.rank == 0:
        print(msg, flush=True)


def _knob():
    v = os.environ.get(KNOB)
    if v is None:
        return False
    if v != "1":
        raise ValueError(f"{KNOB} must be unset or '1', got {v!r}")
    return True


def _run():
    comm = MPI.COMM_WORLD
    t0 = time.perf_counter()
    tips, t1 = _probe_points(ROD_CENTRE)
    _say(f"\nI67 start: -n {comm.size}; source {SOURCE}; balls "
         + ", ".join(f"{n} {tuple(round(x, 6) for x in tips[n])}" for n in BALLS)
         + f"; T1 {tuple(round(x, 6) for x in t1)}")
    store = {}

    # 1. hole build + hole arm (PEC wall on 402).
    built_h = _build_rod()
    c_h = _census(built_h)
    _say(f"I67 (A0) hole build: {c_h['cells']} cells (record {ROD_CELL_RECORD_W1}, ASSERTED exact); "
         f"tag-3 {c_h['t3']}, tag-403 {c_h['t403']}; build {built_h[4]:.1f} s")
    assert c_h["cells"] == ROD_CELL_RECORD_W1, c_h
    hole = _arm(built_h, "hole (PEC rod, wall 402)", rod=True,
                field_hook=_make_hook(store, "hole", tips, t1))
    # 2.
    gc.collect()
    # 3. kept build + gel twin.
    built_k = _build_kept_rod()
    c_k = _census(built_k)
    _say(f"I67 (A0) kept build: {c_k['cells']} cells (record {KEPT_ROD_CELL_RECORD}, ASSERTED exact); "
         f"tag-3 {c_k['t3']}, tag-403 {c_k['t403']} (expected {ROD_CELLS}); build {built_k[4]:.1f} s")
    assert c_k["cells"] == KEPT_ROD_CELL_RECORD, c_k
    gel = _gel_seen(built_k, "gel twin (403 -> gel)", _make_hook(store, "gel", tips, t1))
    del built_k
    gc.collect()
    # 4. the all-cells losses from the store (I57's `_run`, :265–266).
    for arm, key in ((hole, "hole"), (gel, "gel")):
        arm["losses"] = store[key]["losses"]

    _say(f"I67 gel twin: build map had 403 {gel['spy']['build_map_has_403']}, mapped "
         f"{gel['material_403']['mapped']}")
    for arm in (hole, gel):
        _print_gates(arm)
    for pid in PORT_IDS:
        lo = hole["losses"][pid]
        _say(f"I67 hole {pid}: two-tag resid {_two_tag_resid(hole, pid):.3e} vs all-cells; "
             f"P_vol 403 {lo['403']:.3e} W, rest(all-3-403) {lo['rest']:.6e} W, "
             f"all - 3 - cond {lo['all'] - lo['3'] - lo['cond']:.3e} W")
    _say(f"I67 (E) hole arm: tag-402 dofs (owned, reduced) {hole['wall_dofs']}; max |E| per drive "
         + ", ".join(f"{p} {v:.3e}" for p, v in hole["wall_max"].items())
         + f" (bound {PEC_WALL_MAX_E:g}, ASSERTED)")

    # (C1) the rod wall left natural, P1, fresh solver.  Its "natural/without" prints
    # carry predicted ranges written for the no-rod mesh: mislabelled here, print only.
    _say("I67 (C1) follows; its 'natural/without' and 'PEC/natural' predicted ranges were "
         "written for the no-rod 'without' arm and are mislabelled on the gel twin (printed only)")
    c1 = _control_c1({"rod": hole, "off": gel})
    _say(f"I67 (C1) natural wall P1: max |E| on 402 dofs {c1['wall_max']:.3e} (> {PEC_WALL_MAX_E:g}, "
         f"ASSERTED; IMP-3 W1 1.29e-4); (P) {c1['resid']:.3e} (<= {DISCRETE_IDENTITY_RTOL:g}, ASSERTED)")

    hh, hg = store["hole"], store["gel"]
    _say(f"I67 (A0m) owned tag-3 hole {c_h['t3']} / kept {c_k['t3']} (ASSERTED equal); kept tag-403 "
         f"{c_k['t403']} (ASSERTED {ROD_CELLS}); hook census hole {hh['census']} / gel {hg['census']}")
    q = {n: hh["sar"][n] / hg["sar"][n] for n in BALLS}
    eta = float(np.linalg.norm(hh["e_t1"]) / np.linalg.norm(hg["e_t1"]))
    q_arm = {t: hole["sar"][(t, SOURCE)] / gel["sar"][(t, SOURCE)] for t in TIPS}
    _say(f"I67 1 g ball a = {hh['a_1g'] * 1e3:.4f} mm; masses hole / gel twin: "
         + ", ".join(f"{n} {hh['mass'][n] * 1e3:.6f} / {hg['mass'][n] * 1e3:.6f} g" for n in BALLS))
    for n in BALLS:
        _say(f"I67 record Q[{n}] = SAR_1g(hole)/SAR_1g(gel twin) = {q[n]:.6f} (hole {hh['sar'][n]:.6e}, "
             f"gel twin {hg['sar'][n]:.6e} W/kg); I57 kept-cell {I57_KEPT_Q[n]:.6f}; hole/kept "
             f"{q[n] / I57_KEPT_Q[n]:.6f}; entry 48 F-clinical {ENTRY48_Q[n]:.6f}")
    rep = {n: abs(q[n] / I41_APEX_Q[n] - 1.0) for n in I41_APEX_Q}
    _say("I67 (Rep) apex Q vs I41: "
         + ", ".join(f"{n} {q[n]:.6f} vs {I41_APEX_Q[n]:.6f} rel {rep[n]:.3e}" for n in I41_APEX_Q)
         + f" (band {DISCRETE_IDENTITY_RTOL:g}, ASSERTED); `_arm`'s own TIPS ratio (I41's computation) "
         + ", ".join(f"{t} {v:.6f}" for t, v in q_arm.items()))
    _say(f"I67 record eta_rod(T1) = |E_hole|/|E_gel| = {eta:.6f} (valid {hh['t1_valid'] and hg['t1_valid']}); "
         f"I57 kept-cell {I57_KEPT_ETA_T1:.6f}; IMP-3 interval {ENTRY43_READOUT_INTERVAL}; "
         f"entry 48 {ENTRY48_ETA_T1:.6f}")
    _say(f"I67 timing: hole {hole['t_arm']:.1f} / gel twin {gel['t_arm']:.1f} s; total "
         f"{time.perf_counter() - t0:.1f} s at -n {comm.size}")
    return {"hole": hole, "gel": gel, "store": store, "c_h": c_h, "c_k": c_k, "c1": c1,
            "q": q, "rep": rep, "eta": eta}


@pytest.fixture(scope="module")
def i67():
    if not _knob():
        pytest.skip(f"{KNOB}=1 not set (heavy window, ≈ 4-5 min at -n 2)")
    return _run()


# --- coding tests (W0) -------------------------------------------------------

def test_knob_refuses_other_values(monkeypatch):
    monkeypatch.setenv(KNOB, "yes")
    with pytest.raises(ValueError, match=KNOB):
        _knob()
    monkeypatch.delenv(KNOB)
    assert _knob() is False
    monkeypatch.setenv(KNOB, "1")
    assert _knob() is True


def test_records_are_the_pair():
    assert ROD_CELL_RECORD_W1 == 167_371 and KEPT_ROD_CELL_RECORD == 182_371
    assert ROD_CELLS == 15_000


def test_apex_balls_are_i41s_tips_and_t1_outside_rod():
    tips, t1 = _probe_points(ROD_CENTRE)
    assert np.allclose(tips["top_apex"], TIPS["top"], atol=1e-15)
    assert np.allclose(tips["bottom_apex"], TIPS["bottom"], atol=1e-15)
    # T1 lies beyond the rod's apex, so it is in the gel on the hole build too.
    assert t1[2] - ROD_CENTRE[2] > 0.5 * ROD_LENGTH + ROD_RADIUS


# --- W1 ------------------------------------------------------------------------

@complex_only
def test_a0_and_a0m(i67):
    c_h, c_k = i67["c_h"], i67["c_k"]
    assert c_h["cells"] == ROD_CELL_RECORD_W1
    assert c_k["cells"] == KEPT_ROD_CELL_RECORD
    assert c_h["t3"] == c_k["t3"], (c_h, c_k)
    assert c_k["t403"] == ROD_CELLS, c_k
    assert i67["gel"]["material_403"]["mapped"]


@complex_only
@pytest.mark.parametrize("arm", ["hole", "gel"])
def test_arm_gates(i67, arm):
    fails = _gate_failures(i67[arm])
    assert not fails, fails


@complex_only
def test_e_hole_wall(i67):
    h = i67["hole"]
    assert h["wall_dofs"] > 0
    for pid in PORT_IDS:
        assert h["wall_max"][pid] <= PEC_WALL_MAX_E, (pid, h["wall_max"][pid])


@complex_only
def test_control_c1_natural_wall(i67):
    c1 = i67["c1"]
    assert c1["wall_max"] > PEC_WALL_MAX_E, c1
    assert c1["resid"] <= DISCRETE_IDENTITY_RTOL, c1


@complex_only
def test_control_rep_apex_q_is_i41s(i67):
    for n, r in i67["rep"].items():
        assert r <= DISCRETE_IDENTITY_RTOL, (n, i67["q"][n], I41_APEX_Q[n], r)


@complex_only
def test_record_t1_valid(i67):
    # The record carries no band; this is wiring only (η must be a point value).
    st = i67["store"]
    assert st["hole"]["t1_valid"] and st["gel"]["t1_valid"]
