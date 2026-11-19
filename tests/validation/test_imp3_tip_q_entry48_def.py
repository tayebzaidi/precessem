"""`GEO-36` I57 (§9 item 779): the F-small rod's tip Q by entry 48's definition.

Entry 48 (`test_geo36_rod_f_clinical_drive.py`, I36) defines the tip ratio on one
mesh: **rod** arm = tag 403 at ``SIGMA_WIRE_S_PER_M`` with the gel's ε_r / μ_r
(`_m2_failures`, `:160–186`), **gel twin** = tag 403 → the gel's material
(`_gel_twin_map`), Q = SAR_1g(rod) / SAR_1g(gel twin) at the four balls of
``_probe_points`` and η_rod(T1) = |E_rod(T1)| / |E_gel(T1)|.  Here that
definition runs on F-small's only like-for-like fixture: I41's kept build
(``_build_kept_rod``, ``KEPT_ROD_CELL_RECORD`` cells, the rod's cells tag 403),
with I41's four-port drive at 64 MHz, degree 1, source ``SOURCE`` (ccw).

**Premise correction (measured, not assumed).** The brief's rod arm was
``_arm_with_403(map_gel=False)`` on the reading "403 keeps its built σ".  It does
not: ``build_four_port_sweep``'s map has no tag 403 (I41 W1b,
`20261025T201613Z_IMP-3-I41-W1b.log:3775`, "had False"), so the unmapped arm is
**vacuum** — I41's own control label.  Entry 48's rod is σ_wire, so the rod arm
here maps 403 to ``HomogeneousMaterial(SIGMA_WIRE_S_PER_M, gel ε_r, gel μ_r)``
exactly as I35 does (`test_geo36_rod_f_clinical.py:481–483`), through a wrapper
of the same shape as I41's ``_arm_with_403`` (no landed module edited).  Arms
are named by content ("rod", "gel twin"): I41's letters are the reverse.

Gates (imported, unmoved): (A0) global cells == ``KEPT_ROD_CELL_RECORD``;
(A0m) per-tag owned census equal between arms; per arm (i) reciprocity ≤
``RECIPROCITY_BAND``, (ii) σ_max(S) ≤ 1 + ``PASSIVITY_SIGMA_TOLERANCE`` and every
column power ≤ 1, (P) P_src = P_vol + P_sheet within ``DISCRETE_IDENTITY_RTOL``
on every column, with P_vol the all-cells loss integral (§9 item 804, review
ruling 2026-11-03 18:00: ``_exact_shares``' two-tag sum omits tag 403's loss —
item 801 measured that omission as the whole 2.882e-06 P1 residual; the
two-tag residual is still printed and asserted as a control on rod P1).  (E) ``PEC_WALL_MAX_E``: I41 asserts it on no arm, and the kept
build has no PEC wall (402 is interior), so it is printed n/a.  Control: the rod
arm run a second time on the same build reads max_n |Q_repeat/Q_rod − 1| ≤
``DISCRETE_IDENTITY_RTOL`` on the four balls — *asserted at the (P) band, factor 1*
(review ruling 2026-11-04 03:00, §9 item 810: item 779's "exactly" was a bound
no rule backs; item 804 measured 2.63e-13, a leak's signature is ≈ 2.9).

Record (no band): Q on the four balls, η_rod(T1), beside entry 48's F-clinical
readings, ``ENTRY43_READOUT_INTERVAL`` and I41's hole-vs-gel Q.

Scope: one F-small mesh, degree 1, one drive; not converged, no tip-Q band,
nothing clinical.
"""

from __future__ import annotations

import gc
import os
import time
from dataclasses import replace

import numpy as np
import pytest
from mpi4py import MPI

import tests.validation.test_imp3_rod_birdcage as _rb
from fem_em_solver.core import HomogeneousMaterial
from fem_em_solver.io.mesh import BIRDCAGE_IMPLANT_CELL_TAG
from fem_em_solver.post.evaluation import evaluate_vector_field_parallel
from fem_em_solver.post.sar import averaging_ball_radius, mass_averaged_sar
from tests.complex_mode import complex_only
from tests.mesh import test_birdcage_implant_rod as _rodmesh
from tests.mesh.test_birdcage_implant_rod import (
    KEPT_ROD_CELL_RECORD,
    ROD_CENTRE,
    _build_kept_rod,
)
from tests.validation.test_birdcage_power_identity import DISCRETE_IDENTITY_RTOL
from tests.validation.test_birdcage_sar_mass_averaged import PHANTOM_RHO_KG_PER_M3
from tests.validation.test_geo36_rod_f_clinical_drive import _probe_points
from tests.validation.test_imp3_rod_birdcage import PEC_WALL_MAX_E, PORT_IDS
from tests.validation.test_imp3_rod_heating_matched import _arm_with_403, _census
from tests.validation.test_imp3_rod_heating_phantom import SOURCE
from tests.validation.test_imp3_tip_ladder_1g import ENTRY43_READOUT_INTERVAL
from tests.validation.test_mass_averaged_sar_standard_masses import (
    ONE_GRAM_KG,
    QUADRATURE_DEGREE,
)
from tests.validation.test_port_birdcage_four_port import PASSIVITY_SIGMA_TOLERANCE
from tests.validation.test_port_birdcage_lumped_column import PHANTOM_CELL_TAG
from tests.validation.test_port_drive_superposition import TWO_TAG_ACK_KEY
from tests.validation.test_port_gap_voltage_impedance import SIGMA_WIRE_S_PER_M
from tests.validation.test_port_lumped_sheet_sweep import RECIPROCITY_BAND

KNOB = "FEM_EM_IMP3_I57"
# Entry 48's F-clinical readings (`20261102T080008Z_GEO-36-step9-I36.log:12501`).
ENTRY48_Q = {"top": 4.460568, "bottom": 4.428040, "top_apex": 4.413534, "bottom_apex": 4.341688}
ENTRY48_ETA_T1 = 3.642324
# I41's hole-vs-gel Q, ccw (`20261025T201613Z_IMP-3-I41-W1b.log:3776–3777`).
I41_HOLE_Q = {"top": 4.090, "bottom": 4.114}
BALLS = ("top", "bottom", "top_apex", "bottom_apex")


def _say(msg):
    if MPI.COMM_WORLD.rank == 0:
        print(msg, flush=True)


def _rod_material(gel):
    return HomogeneousMaterial(sigma=SIGMA_WIRE_S_PER_M, epsilon_r=gel.epsilon_r, mu_r=gel.mu_r)


def _arm_rod(built, label, field_hook):
    """``_rb._arm`` (rod=False) on the kept build with 403 → σ_wire (I35's map)."""
    orig = _rb.build_four_port_sweep
    seen = {}

    def _wrapped(*a, **k):
        case = orig(*a, **k)
        mmap = dict(case["problem"].material_map)
        seen["had_403"] = BIRDCAGE_IMPLANT_CELL_TAG in mmap
        gel = mmap[PHANTOM_CELL_TAG]
        case["problem"] = replace(
            case["problem"], material_map={**mmap, BIRDCAGE_IMPLANT_CELL_TAG: _rod_material(gel)}
        )
        case[TWO_TAG_ACK_KEY] = "(P) gated over all cells by _all_cells_resid"
        seen["mat"] = case["problem"].material_map[BIRDCAGE_IMPLANT_CELL_TAG]
        seen["gel"] = gel
        return case

    _rb.build_four_port_sweep = _wrapped
    try:
        out = _rb._arm(built, label, rod=False, field_hook=field_hook)
    finally:
        _rb.build_four_port_sweep = orig
    out["material_403"] = seen
    return out


def _gel_seen(built, label, field_hook):
    """I41's ``_arm_with_403(map_gel=True)``, with the mapped material recorded."""
    orig = _rb.build_four_port_sweep
    seen = {}

    def _spy(*a, **k):
        case = orig(*a, **k)
        seen["build_map_has_403"] = BIRDCAGE_IMPLANT_CELL_TAG in case["problem"].material_map
        return case

    _rb.build_four_port_sweep = _spy
    try:
        out = _arm_with_403(built, label, field_hook, map_gel=True)
    finally:
        _rb.build_four_port_sweep = orig
    out["spy"] = seen
    return out


def _owned_tag_counts(msh, cell_tags):
    n_owned = msh.topology.index_map(msh.topology.dim).size_local
    idx = np.asarray(cell_tags.indices)
    vals = np.asarray(cell_tags.values)[idx < n_owned]
    tags = sorted(msh.comm.allreduce(set(int(t) for t in np.unique(vals)), op=MPI.BOR))
    return {t: msh.comm.allreduce(int(np.count_nonzero(vals == t)), op=MPI.SUM) for t in tags}


def _make_hook(store, key, tips, t1):
    def hook(label, rod, msh, cell_tags, e_of, sigma_field, rho_field):
        comm = msh.comm
        e = e_of[SOURCE]
        a_1g = averaging_ball_radius(mass_kg=ONE_GRAM_KG, rho=PHANTOM_RHO_KG_PER_M3)
        sar, mass = {}, {}
        for name in BALLS:
            o = mass_averaged_sar(
                e,
                sigma=sigma_field,
                rho=rho_field,
                center=tips[name],
                radius=a_1g,
                quadrature_degree=QUADRATURE_DEGREE,
            )
            sar[name] = float(o["averaged_sar_w_per_kg"])
            mass[name] = float(o["mass_kg"])
        pts = np.ascontiguousarray(np.array([t1], dtype=np.float64))
        vals, mask = evaluate_vector_field_parallel(e, pts, owner="lowest_original_cell")
        vals = comm.bcast(vals if comm.rank == 0 else None, root=0)
        mask = comm.bcast(mask if comm.rank == 0 else None, root=0)
        # §9 item 804: (P)'s P_vol over all cells (entry 48's form), per drive, with the
        # per-tag split so the two-tag helper's omitted 403 term stays visible.
        losses = {
            pid: _tag_losses(msh, cell_tags, e_of[pid], sigma_field, None)[0] for pid in PORT_IDS
        }
        store[key] = {
            "losses": losses,
            "sar": sar,
            "mass": mass,
            "a_1g": a_1g,
            "e_t1": np.asarray(vals).reshape(1, 3)[0],
            "t1_valid": bool(np.all(mask)),
            "census": _owned_tag_counts(msh, cell_tags),
            "cells": int(msh.topology.index_map(msh.topology.dim).size_global),
        }

    return hook


def _gate_failures(arm):
    mq = arm["mq"]
    fails = []
    if not mq["recip"] <= RECIPROCITY_BAND:
        fails.append(f"{arm['label']} (i) {mq['recip']:.3e}")
    if not float(mq["sigma"][0]) <= 1.0 + PASSIVITY_SIGMA_TOLERANCE:
        fails.append(f"{arm['label']} (ii) sigma_max {float(mq['sigma'][0]):.12f}")
    if not np.all(mq["colpow"] <= 1.0):
        fails.append(f"{arm['label']} (ii) column powers {mq['colpow']}")
    for pid in PORT_IDS:
        resid, form = _p_resid(arm, pid)
        if not resid <= DISCRETE_IDENTITY_RTOL:
            fails.append(f"{arm['label']} (P) {pid} {resid:.3e} ({form})")
    return fails


def _two_tag_resid(arm, pid):
    """The helper's form: P_vol = `_exact_shares` p_vol (tags 3 + conductor only)."""
    ex = arm["shares"][pid]
    return abs(ex["p_src"] - ex["p_vol"] - ex["sheet_field_total"]) / abs(ex["p_src"])


def _all_cells_resid(arm, pid):
    """Entry 48's form: P_vol = ½∫σ|E|² over all cells (`test_geo36_rod_f_clinical_drive.py:381`)."""
    ex = arm["shares"][pid]
    p_all = arm["losses"][pid]["all"]
    return abs(ex["p_src"] - p_all - ex["sheet_field_total"]) / abs(ex["p_src"])


def _p_resid(arm, pid):
    """(P) residual asserted: all cells when the arm carries per-tag losses (item 804), else
    the two-tag form (the item 801 diag arm, whose pre-registered red is the two-tag one)."""
    if "losses" in arm:
        return _all_cells_resid(arm, pid), "all cells"
    return _two_tag_resid(arm, pid), "two-tag"


def _print_gates(arm):
    mq = arm["mq"]
    _say(
        f"[I57 {arm['label']}] (i) ||S-S^T||/||S|| = {mq['recip']:.3e} (band {RECIPROCITY_BAND:g}, "
        f"ASSERTED); (ii) sigma_max(S) = {float(mq['sigma'][0]):.9f} (<= 1 + "
        f"{PASSIVITY_SIGMA_TOLERANCE:g}, ASSERTED), column powers "
        + ", ".join(f"{x:.6f}" for x in mq["colpow"])
        + " (<= 1, ASSERTED)"
    )
    for pid in PORT_IDS:
        ex = arm["shares"][pid]
        resid, form = _p_resid(arm, pid)
        _say(
            f"    (P) {pid}: rel {resid:.3e} ({form}, band {DISCRETE_IDENTITY_RTOL:g}, ASSERTED); "
            f"phantom share {ex['phantom'] / ex['p_src']:.4e} of P_src"
        )
        if "losses" in arm:
            lo = arm["losses"][pid]
            _say(
                f"        P_src {ex['p_src']:.9e} W, sheets {ex['sheet_field_total']:.9e} W; P_vol by tag: "
                f"3 {lo['3']:.9e}, 403 {lo['403']:.9e}, conductor {lo['cond']:.9e}, all {lo['all']:.9e} W; "
                f"two-tag (_exact_shares p_vol {ex['p_vol']:.9e}) resid {_two_tag_resid(arm, pid):.3e} "
                f"(printed, not asserted) beside all-cells {_all_cells_resid(arm, pid):.3e}"
            )
    _say(f"    times: arm {arm['t_arm']:.1f} s, sweep {arm['t_sweep']:.2f} s")


def _run():
    comm = MPI.COMM_WORLD
    t0 = time.perf_counter()
    tips, t1 = _probe_points(ROD_CENTRE)
    _say(
        f"\nI57 start: -n {comm.size}; source {SOURCE}; balls "
        + ", ".join(f"{n} {tuple(round(x, 6) for x in tips[n])}" for n in BALLS)
        + f"; T1 {tuple(round(x, 6) for x in t1)}"
    )
    _say(
        "I57 sheets (the kept build's arguments, `_build_rod` + implant_rod keep_cells=True): "
        f"emit_port_sheets=True, port_box_size={_rodmesh.PORT_BOX_SIZE!r}, "
        f"leg_gap_length={_rodmesh.LEG_GAP_LENGTH!r}, leg_width={_rodmesh.LEG_WIDTH!r}, "
        f"leg_count={_rodmesh.LEG_COUNT!r}, resolution={_rodmesh.RESOLUTION!r}, "
        f"conductor_resolution={_rodmesh.CONDUCTOR_RESOLUTION!r}"
    )

    built = _build_kept_rod()
    c_b = _census(built)
    _say(
        f"I57 (A0) kept build: {c_b['cells']} cells (record {KEPT_ROD_CELL_RECORD}, ASSERTED "
        f"exact); tag-3 {c_b['t3']}, tag-403 {c_b['t403']}; build {built[4]:.1f} s"
    )
    assert c_b["cells"] == KEPT_ROD_CELL_RECORD, c_b

    store = {}
    rod = _arm_rod(built, "rod (403 sigma_wire)", _make_hook(store, "rod", tips, t1))
    gc.collect()
    gel = _gel_seen(built, "gel twin (403 -> gel)", _make_hook(store, "gel", tips, t1))
    gc.collect()
    rod2 = _arm_rod(built, "rod, repeat (control)", _make_hook(store, "rod2", tips, t1))
    gc.collect()
    for arm, key in ((rod, "rod"), (gel, "gel"), (rod2, "rod2")):
        arm["losses"] = store[key]["losses"]

    m = rod["material_403"]
    _say(
        f"I57 tag 403 material: build map had 403 {m['had_403']} (gel twin spy "
        f"{gel['spy']['build_map_has_403']}; unmapped = vacuum, so map_gel=False is not the rod); "
        f"rod sigma {m['mat'].sigma!r} S/m (SIGMA_WIRE_S_PER_M {SIGMA_WIRE_S_PER_M!r}), eps_r "
        f"{m['mat'].epsilon_r!r}, mu_r {m['mat'].mu_r!r}; gel sigma {m['gel'].sigma!r}, eps_r "
        f"{m['gel'].epsilon_r!r}; gel twin mapped {gel['material_403']['mapped']}"
    )
    for arm in (rod, gel, rod2):
        _print_gates(arm)
    _say(
        f"I57 (E) PEC_WALL_MAX_E {PEC_WALL_MAX_E:g}: n/a — I41 asserts it on no arm and the "
        "kept build has no PEC wall (tag 402 interior)"
    )

    hr, hg, hr2 = store["rod"], store["gel"], store["rod2"]
    _say(
        f"I57 (A0m) owned census rod {hr['census']} / gel twin {hg['census']} / repeat "
        f"{hr2['census']} (ASSERTED equal); cells {hr['cells']} / {hg['cells']} / {hr2['cells']}"
    )
    q = {n: hr["sar"][n] / hg["sar"][n] for n in BALLS}
    q_ctl = {n: hr2["sar"][n] / hr["sar"][n] for n in BALLS}
    eta = float(np.linalg.norm(hr["e_t1"]) / np.linalg.norm(hg["e_t1"]))
    _say(
        f"I57 1 g ball a = {hr['a_1g'] * 1e3:.4f} mm; masses "
        + ", ".join(f"{n} {hr['mass'][n] * 1e3:.6f} g" for n in BALLS)
    )
    for n in BALLS:
        hole = f"{I41_HOLE_Q[n]:.3f}" if n in I41_HOLE_Q else "-"
        _say(
            f"I57 record Q[{n}] = SAR_1g(rod)/SAR_1g(gel twin) = {q[n]:.6f} (rod {hr['sar'][n]:.6e}, "
            f"gel twin {hg['sar'][n]:.6e} W/kg); entry 48 F-clinical {ENTRY48_Q[n]:.6f}; I41 "
            f"hole-vs-gel {hole}"
        )
    _say(
        f"I57 record eta_rod(T1) = |E_rod|/|E_gel| = {eta:.6f} (valid {hr['t1_valid'] and hg['t1_valid']}); "
        f"entry 48 {ENTRY48_ETA_T1:.6f}; IMP-3 hole-vs-gel interval {ENTRY43_READOUT_INTERVAL}"
    )
    _say(
        "I57 control (ASSERTED exact): rod repeat / rod Q = "
        + ", ".join(f"{n} {q_ctl[n]!r}" for n in BALLS)
    )
    _say(
        f"I57 timing: rod {rod['t_arm']:.1f} / gel twin {gel['t_arm']:.1f} / repeat "
        f"{rod2['t_arm']:.1f} s; total {time.perf_counter() - t0:.1f} s at -n {comm.size}"
    )
    return {
        "rod": rod,
        "gel": gel,
        "rod2": rod2,
        "store": store,
        "q": q,
        "q_ctl": q_ctl,
        "eta": eta,
        "c_b": c_b,
    }


@pytest.fixture(scope="module")
def i57():
    if os.environ.get(KNOB) != "1":
        pytest.skip(f"{KNOB}=1 not set (heavy window, ≈ 4-5 min at -n 2)")
    return _run()


# --- coding tests (W0) -------------------------------------------------------


def test_knob_refuses_without_value():
    assert os.environ.get(KNOB) in (None, "1"), f"{KNOB} must be unset or 1"


def test_rod_material_is_entry48s():
    gel = HomogeneousMaterial(sigma=0.5, epsilon_r=77.0, mu_r=1.0)
    rod = _rod_material(gel)
    assert rod.sigma == SIGMA_WIRE_S_PER_M
    assert rod.epsilon_r == gel.epsilon_r and rod.mu_r == gel.mu_r


def test_probe_points_on_f_small():
    tips, t1 = _probe_points(ROD_CENTRE)
    assert tips["top"][2] == pytest.approx(0.5 * _rodmesh.ROD_LENGTH)
    assert tips["top_apex"][2] == pytest.approx(0.5 * _rodmesh.ROD_LENGTH + _rodmesh.ROD_RADIUS)
    assert tips["bottom"][0] == pytest.approx(ROD_CENTRE[0])
    assert len(t1) == 3


# --- W1 ------------------------------------------------------------------------


@complex_only
def test_a0_and_a0m_census(i57):
    st = i57["store"]
    assert i57["c_b"]["cells"] == KEPT_ROD_CELL_RECORD
    assert st["rod"]["cells"] == st["gel"]["cells"] == KEPT_ROD_CELL_RECORD
    assert st["rod"]["census"] == st["gel"]["census"] == st["rod2"]["census"]
    assert st["rod"]["census"].get(BIRDCAGE_IMPLANT_CELL_TAG, 0) > 0
    assert i57["rod"]["material_403"]["mat"].sigma == SIGMA_WIRE_S_PER_M
    assert i57["gel"]["material_403"]["mapped"]


@complex_only
@pytest.mark.parametrize("arm", ["rod", "gel", "rod2"])
def test_arm_gates(i57, arm):
    fails = _gate_failures(i57[arm])
    assert not fails, fails


@complex_only
def test_control_two_tag_p1_reproduces_omitted_term(i57):
    """§9 item 804 control: the rod arm's two-tag (P) residual on P1 reproduces items 779/801's
    2.882e-06 to three digits (`20261103T195044Z_GEO-36-I57-diag-W1.log:1942`), so the term the
    all-cells form adds is the omitted tag-403 loss."""
    r = _two_tag_resid(i57["rod"], "P1")
    _say(
        f"I57 control: rod P1 two-tag |resid| {r:.4e} vs {W1_ROD_P1_RESID:.3e} (ASSERTED 3 digits)"
    )
    assert f"{r:.3e}" == f"{W1_ROD_P1_RESID:.3e}", r


@complex_only
def test_control_rod_twice_q_is_one(i57):
    dev = {n: abs(float(i57["q_ctl"][n]) - 1.0) for n in BALLS}
    worst = max(dev.values())
    _say(
        f"    control repeat/rod: max |Q - 1| = {worst:.3e} (band {DISCRETE_IDENTITY_RTOL:g}, "
        "ASSERTED); item 779 0.0, item 804 2.63e-13 "
        "(`20261104T013204Z_GEO-36-I57-W1c.log:2004`)"
    )
    assert worst <= DISCRETE_IDENTITY_RTOL, dev


# --- §9 item 801: the rod arm's (P) red separated (quadrature or solve floor) --------

DIAG_KNOB = "FEM_EM_IMP3_I57_DIAG"
#: Item 779's W1 rod-arm P1 residual (`20261103T120306Z_GEO-36-I57-W1.log:1968`).
W1_ROD_P1_RESID = 2.882e-06
#: The identity's P_vol tags (`_loss_power_w`: PHANTOM_CELL_TAG + CONDUCTOR_CELL_TAG).
DIAG_DEGREES = (2, 4)


def _tag_losses(msh, cell_tags, e, sigma, degree):
    """½∫σ|E|² over tag 3, tag 403, the identity's conductor tag and all cells, all-reduced."""
    import ufl
    from dolfinx import fem

    from tests.validation.test_port_birdcage_lumped_column import CONDUCTOR_CELL_TAG

    md = {} if degree is None else {"quadrature_degree": int(degree)}
    dx = ufl.Measure("dx", domain=msh, subdomain_data=cell_tags, metadata=md)
    integrand = 0.5 * sigma * ufl.inner(e, e)
    comm = msh.comm

    def _r(form):
        return float(np.real(comm.allreduce(fem.assemble_scalar(fem.form(form)), op=MPI.SUM)))

    out = {
        "3": _r(integrand * dx(PHANTOM_CELL_TAG)),
        "403": _r(integrand * dx(BIRDCAGE_IMPLANT_CELL_TAG)),
        "cond": _r(integrand * dx(CONDUCTOR_CELL_TAG)),
        "all": _r(integrand * dx),
    }
    out["rest"] = out["all"] - out["3"] - out["403"]
    est = ufl.algorithms.estimate_total_polynomial_degree(integrand)
    return out, int(est)


def _residual_recorder(monkeypatch, rows):
    """Wrap both solve paths: after each call, ‖A e − b‖/‖b‖ on the held LinearProblem."""
    from fem_em_solver.core import TimeHarmonicSolver

    def _rec(self, fields, kind):
        lp = self._linear_problem
        a_op, b = lp.A, lp.b
        x = fields.e_complex.x.petsc_vec
        y = b.duplicate()
        a_op.mult(x, y)
        y.axpy(-1.0, b)
        rows.append({"kind": kind, "r": float(y.norm()) / float(b.norm()), "b": float(b.norm())})
        y.destroy()

    orig_solve = TimeHarmonicSolver.solve
    orig_held = TimeHarmonicSolver.solve_with_held_factorization

    def _solve(self, *a, **k):
        f = orig_solve(self, *a, **k)
        _rec(self, f, "factorise")
        return f

    def _held(self, *a, **k):
        f = orig_held(self, *a, **k)
        _rec(self, f, "held")
        return f

    monkeypatch.setattr(TimeHarmonicSolver, "solve", _solve)
    monkeypatch.setattr(TimeHarmonicSolver, "solve_with_held_factorization", _held)


def _run_diag(monkeypatch):
    comm = MPI.COMM_WORLD
    t0 = time.perf_counter()
    built = _build_kept_rod()
    c_b = _census(built)
    _say(
        f"\nI57-diag (A0) kept build: {c_b['cells']} cells (record {KEPT_ROD_CELL_RECORD}, ASSERTED "
        f"exact); tag-403 {c_b['t403']}; build {built[4]:.1f} s; -n {comm.size}; rod arm only"
    )
    keep = {}

    def hook(label, rod, msh, cell_tags, e_of, sigma_field, rho_field):
        keep.update(msh=msh, cell_tags=cell_tags, e_of=dict(e_of), sigma=sigma_field)

    rows = []
    with monkeypatch.context() as m:
        _residual_recorder(m, rows)
        rod = _arm_rod(built, "rod (403 sigma_wire), diag", hook)
    _print_gates(rod)
    mat = rod["material_403"]["mat"]
    _say(
        f"I57-diag tag 403 material: sigma {mat.sigma!r} S/m, eps_r {mat.epsilon_r!r}; solves recorded "
        f"{[r['kind'] for r in rows]}"
    )
    msh, ct, sig = keep["msh"], keep["cell_tags"], keep["sigma"]
    diag = {}
    for j, pid in enumerate(PORT_IDS):
        ex = rod["shares"][pid]
        p_src, p_sh = ex["p_src"], ex["sheet_field_total"]
        resid_id = (p_src - ex["p_vol"] - p_sh) / p_src
        per = {}
        for d in (None,) + DIAG_DEGREES:
            losses, est = _tag_losses(msh, ct, keep["e_of"][pid], sig, d)
            per[d] = {
                "losses": losses,
                "est": est,
                "r_id": (p_src - losses["3"] - losses["cond"] - p_sh) / p_src,
                "r_all": (p_src - losses["all"] - p_sh) / p_src,
            }
        r = rows[j] if j < len(rows) else None
        diag[pid] = {"resid_identity": resid_id, "per": per, "r": r}
        _say(
            f"I57-diag {pid}: P_src {p_src:.9e} W, sheets {p_sh:.9e} W; identity (_exact_shares) "
            f"P_vol {ex['p_vol']:.9e} W (phantom {ex['phantom']:.9e} + conductor {ex['conductor']:.9e}); "
            f"(P) rel {resid_id:+.3e}"
        )
        for d, v in per.items():
            lo = v["losses"]
            tag = "unpinned" if d is None else f"d={d}"
            _say(
                f"    (q) {tag} (estimated degree {v['est']}): P_vol tag 3 {lo['3']:.9e}, tag 403 "
                f"{lo['403']:.9e}, conductor tag {lo['cond']:.9e}, rest(all-3-403) {lo['rest']:.9e}, all "
                f"{lo['all']:.9e} W; resid on identity tags (3 + cond) {v['r_id']:+.3e}; resid on all "
                f"cells {v['r_all']:+.3e}"
            )
        if r is None:
            _say(f"    (r) n/a: {len(rows)} solves recorded for {len(PORT_IDS)} drives")
        else:
            _say(f"    (r) ||A e - b||/||b|| = {r['r']:.3e} ({r['kind']}, ||b|| {r['b']:.6e})")
    _say(f"I57-diag total {time.perf_counter() - t0:.1f} s; arm {rod['t_arm']:.1f} s")
    return {"rod": rod, "c_b": c_b, "diag": diag, "rows": rows}


@pytest.fixture(scope="module")
def i57_diag():
    if os.environ.get(DIAG_KNOB) != "1":
        pytest.skip(f"{DIAG_KNOB}=1 not set (heavy window, rod arm only)")
    mp = pytest.MonkeyPatch()
    try:
        yield _run_diag(mp)
    finally:
        mp.undo()


@complex_only
def test_diag_gates_except_p(i57_diag):
    """(A0), (i), (ii) on the rod arm, imported bands unmoved."""
    rod = i57_diag["rod"]
    assert i57_diag["c_b"]["cells"] == KEPT_ROD_CELL_RECORD
    fails = [f for f in _gate_failures(rod) if " (P) " not in f]
    assert not fails, fails


@complex_only
def test_diag_control_reproduces_w1(i57_diag):
    """The run is the parked run: P1's identity residual reproduces W1's to three digits."""
    r = abs(i57_diag["diag"]["P1"]["resid_identity"])
    _say(f"I57-diag control: P1 |resid| {r:.4e} vs W1 {W1_ROD_P1_RESID:.3e}")
    assert f"{r:.3e}" == f"{W1_ROD_P1_RESID:.3e}", r


@complex_only
def test_diag_p_rod_arm(i57_diag):
    """(P) at DISCRETE_IDENTITY_RTOL. PREDICTED RED: (P) rod arm P1 range [2.5e-06, 3.3e-06]."""
    fails = [f for f in _gate_failures(i57_diag["rod"]) if " (P) " in f]
    assert not fails, fails


@complex_only
def test_record_finite_and_t1_valid(i57):
    # The record carries no band (brief: "a record, no band"); this is wiring only.
    st = i57["store"]
    assert st["rod"]["t1_valid"] and st["gel"]["t1_valid"]
    assert all(np.isfinite(i57["q"][n]) and i57["q"][n] > 0 for n in BALLS)
