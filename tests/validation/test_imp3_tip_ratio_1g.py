"""`IMP-3` step 6 (§9 item 263, chain I14): the tip ratio on the 1 g rung at degree 2.

Item 235's exact-rectangle rod fixture (``leg_sheet_strip_fraction =
GATED_WIDTH_FRACTION``, κ off) built at ``PHANTOM_RESOLUTION_1G_RUNG``, 64 MHz,
off and rod-PEC arms, solved at degree 1 and degree 2.  One env knob
``FEM_EM_IMP3_STEP6`` in {build, d2off, d2rod, d1}; one test dispatches on it.

Asserted (bands imported, never restated):
    W1 build: (G1) rod cells < ROD_CELL_GUARD_W2; (G2) N_i areas within
        STRIP_BAND of 5.6e-5, both arms; (G3) tag 402 > 0 and all exterior on
        the rod arm, 0 on the off arm; (G4) the 1 g ball fits.
    W2/W3/W4 per arm and degree: (i) reciprocity, (ii) passivity, (PP) the
        PORT-16 exact power partition, (M) rod σ_y pairs, (S) off separation,
        (E) rod PEC wall, (V) six points valid, (D) dofs = 2E + 2F (degree 2) /
        E (degree 1).
    W4 only: (Kw) the rod wall left natural (P1, degree 1, fresh solver) has
        max |E| on the 402 dofs > PEC_WALL_MAX_E.

Printed, never gated: η_rod at both orders, its p-move beside −4.961 % (record
mesh), the record-mesh -> 1 g rung move, Q at both tips and orders, the closed
form, (C-p), (C2), cells, dofs, sweep seconds, summed peak RSS.  No accuracy
claim (`ANS-10`).
"""

from __future__ import annotations

import gc
import json
import os
import time
from pathlib import Path

import numpy as np
import pytest
from mpi4py import MPI

from fem_em_solver.ports.sparameters import run_n_port_sparameter_sweep
from fem_em_solver.ports.superposition import superpose_drives
from fem_em_solver.post.evaluation import evaluate_vector_field_parallel
from fem_em_solver.post.sar import (
    averaging_ball_radius,
    build_density_field,
    mass_averaged_sar,
)
from fem_em_solver.utils.analytical import pec_prolate_spheroid_field
from tests.complex_mode import complex_only
from tests.mesh.helpers import global_cell_tag_set
from tests.mesh.test_birdcage_implant_rod import ROD_CENTRE, _build_rod
from tests.mesh.test_birdcage_port_sheet_prerequisite import CELL_COUNT_BAND
from tests.mesh.test_birdcage_port_sheets import SHEET_IFACE, _build
from tests.mesh.test_birdcage_port_tags import PHANTOM_HEIGHT, PHANTOM_RADIUS
from tests.validation.test_birdcage_b1_quadrature import (
    QUADRATURE_STEP_DEG,
    _port_index,
    quadrature_phase_weights,
)
from tests.validation.test_birdcage_power_identity import (
    DISCRETE_IDENTITY_RTOL,
    STEP2D_FREQUENCY_HZ,
    _exact_shares,
)
from tests.validation.test_birdcage_sar_1g_rung import (
    ONE_GRAM_RUNG_CELL_RECORD,
    PHANTOM_RESOLUTION_1G_RUNG,
)
from tests.validation.test_birdcage_sar_mass_averaged import PHANTOM_RHO_KG_PER_M3
from tests.validation.test_imp3_palace_fields import PNAMES, PTS
from tests.validation.test_imp3_palace_twin import _sheet_areas
from tests.validation.test_imp3_rod_birdcage import (
    DRIVES,
    PEC_WALL_MAX_E,
    PORT_IDS,
    ROD_CELL_GUARD_W2,
    TIPS,
    _control_c1,
    _matrix_quantities,
    _outer_box_tags,
    _peak_rss_gib,
    _print_arm,
    _wall_dofs,
    _wall_max_e,
)
from tests.validation.test_mass_averaged_sar_standard_masses import (
    ONE_GRAM_KG,
    QUADRATURE_DEGREE,
)
from tests.validation.test_ops61_exact_strip_ports import STRIP_BAND
from tests.validation.test_port_birdcage_four_port import (
    PASSIVITY_SIGMA_TOLERANCE,
    POOLED_SEPARATION_FLOOR,
    build_four_port_sweep,
)
from tests.validation.test_port_birdcage_lumped_column import (
    ADJACENT_SPREAD_BAND,
    PHANTOM_CELL_TAG,
)
from tests.validation.test_port_lumped_narrowed_sheet import GATED_WIDTH_FRACTION
from tests.validation.test_port_lumped_sheet_sweep import RECIPROCITY_BAND

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "logs" / "imp3"
KNOB = os.environ.get("FEM_EM_IMP3_STEP6")
KNOBS = ("build", "d2off", "d2rod", "d1")
NOMINAL_STRIP_AREA = 5.6e-5
TAG = "[IMP-3 6]"
# Record-mesh readings (item 252 W5, `20261005T173453Z_IMP-3-item252-W5-ours-p2.log:3910-3911`).
RECORD_ETA_O2_T1 = 3.594534
RECORD_PMOVE_T1 = -4.961114e-02
# Step 1's ragged 1 g rung, degree 1, P1 top (`20261002T170843Z_IMP-3-step1-W2.log:3864`).
STEP1_RAGGED_Q_P1_TOP = 4.496087
PRED_CP_FLOOR = 1.0e-3
PRED_C2_FLOOR = 3.0


def _say(msg):
    if MPI.COMM_WORLD.rank == 0:
        print(msg, flush=True)


def _built(arm):
    kw = {
        "phantom_resolution": PHANTOM_RESOLUTION_1G_RUNG,
        "leg_sheet_strip_fraction": GATED_WIDTH_FRACTION,
    }
    if arm == "off":
        return _build(True, **kw)
    return _build_rod(**kw)


def _edge_face_counts(msh):
    msh.topology.create_entities(1)
    msh.topology.create_entities(2)
    return (int(msh.topology.index_map(1).size_global), int(msh.topology.index_map(2).size_global))


def _eval_points(e_of, comm):
    e = np.zeros((len(PTS), 3, len(PORT_IDS)), dtype=complex)
    valid = True
    for j, p in enumerate(PORT_IDS):
        vals, mask = evaluate_vector_field_parallel(e_of[p], PTS, owner="lowest_original_cell")
        vals = comm.bcast(vals if comm.rank == 0 else None, root=0)
        mask = comm.bcast(mask if comm.rank == 0 else None, root=0)
        valid = valid and bool(np.all(mask))
        e[:, :, j] = np.asarray(vals).reshape(len(PTS), 3)
    return e, valid


def _arm_d(built, label, rod, degree):
    """`_arm` (test_imp3_rod_birdcage.py:295) with ``degree`` passed to the build."""
    comm = MPI.COMM_WORLD
    t_arm = time.perf_counter()
    case = build_four_port_sweep(
        frequency_hz=STEP2D_FREQUENCY_HZ, built=built, build_only=True, degree=degree
    )
    assert int(case["degree"]) == degree, case["degree"]
    msh = case["mesh"]
    ncells = int(case["cells"])
    areas = _sheet_areas(msh, case["facet_tags"], comm)
    _say(
        f"\n{TAG} {label}: {ncells} cells, mesh wall {built[3]['mesh_wall_time_s']:.2f} s, "
        f"build {built[4]:.2f} s; N_i areas "
        + ", ".join(f"P{i} {a:.12e}" for i, a in areas.items())
    )
    if rod:
        assert ncells < ROD_CELL_GUARD_W2, (ncells, ROD_CELL_GUARD_W2)  # (G1) before the solve
    comm.Barrier()
    t0 = time.perf_counter()
    result = run_n_port_sparameter_sweep(
        case["problem"],
        case["port_defs"],
        lumped_sheet_ports=case["specs"],
        lumped_sheet_facet_tags=case["facet_tags"],
        keep_fields=True,
        degree=degree,
    )
    comm.Barrier()
    t_sweep = time.perf_counter() - t0
    assert tuple(result.port_ids) == PORT_IDS, result.port_ids
    mq = _matrix_quantities(result)
    omega = 2.0 * np.pi * float(STEP2D_FREQUENCY_HZ)
    shares = {
        pid: _exact_shares(case, {"driven": pid, "omega": omega, "fields": result.fields[pid]})
        for pid in PORT_IDS
    }
    azimuths = {f"P{s['tag'] - SHEET_IFACE}": float(s["azimuth_deg"]) for s in case["sheets"]}
    ks = np.array(
        [_port_index(azimuths[p], azimuths["P1"], QUADRATURE_STEP_DEG) for p in PORT_IDS],
        dtype=float,
    )
    weights = {sense: quadrature_phase_weights(ks, sense) for sense in ("ccw", "cw")}
    superposed = {sense: superpose_drives(result, weights[sense]) for sense in weights}
    e_of = {pid: result.fields[pid].e_complex for pid in PORT_IDS}
    e_of.update({sense: superposed[sense].e_complex for sense in superposed})
    sigma_field = result.fields["P1"].sigma_field

    cell_tags = case["cell_tags"]
    all_tags = global_cell_tag_set(msh, cell_tags)
    rho_field = build_density_field(
        msh,
        PHANTOM_RHO_KG_PER_M3,
        cell_tags=cell_tags,
        density_map={int(t): 0.0 for t in sorted(all_tags) if int(t) != PHANTOM_CELL_TAG},
    )
    a_1g = averaging_ball_radius(mass_kg=ONE_GRAM_KG, rho=PHANTOM_RHO_KG_PER_M3)
    t1 = time.perf_counter()
    sar, mass = {}, {}
    for tip, centre in TIPS.items():
        for d in DRIVES:
            o = mass_averaged_sar(
                e_of[d],
                sigma=sigma_field,
                rho=rho_field,
                center=centre,
                radius=a_1g,
                quadrature_degree=QUADRATURE_DEGREE,
            )
            sar[(tip, d)] = float(o["averaged_sar_w_per_kg"])
            mass[tip] = float(o["mass_kg"])
    t_sar = time.perf_counter() - t1

    v_space = e_of["P1"].function_space
    ndofs = int(v_space.dofmap.index_map.size_global * v_space.dofmap.index_map_bs)
    n_e, n_f = _edge_face_counts(msh)
    e_pts, valid = _eval_points(e_of, comm)
    out = {
        "label": label,
        "cells": ncells,
        "areas": areas,
        "mq": mq,
        "shares": shares,
        "sar": sar,
        "mass": mass,
        "a_1g": a_1g,
        "t_sweep": t_sweep,
        "t_sar": t_sar,
        "ndofs": ndofs,
        "edges": n_e,
        "faces": n_f,
        "E": e_pts,
        "valid": valid,
        "degree": degree,
    }
    if rod:
        owned, n_wall = _wall_dofs(msh, built[2], v_space)
        out["wall_dofs"] = n_wall
        out["wall_max"] = {pid: _wall_max_e(e_of[pid], owned, comm) for pid in PORT_IDS}
        out["case"] = case
        out["built"] = built
        out["owned_wall"] = owned
        out["rho_field"] = rho_field
    out["t_arm"] = time.perf_counter() - t_arm
    del superposed, e_of, result, v_space
    gc.collect()
    return out


def _identities(a, rod):
    """Print, then return the list of failed anchors for one arm at one degree."""
    mq, deg = a["mq"], a["degree"]
    _print_arm(a)
    pred_d = 2 * a["edges"] + 2 * a["faces"] if deg == 2 else a["edges"]
    _say(
        f"    (D) degree {deg}: global dofs {a['ndofs']} vs "
        f"{'2E + 2F' if deg == 2 else 'E'} = {pred_d} (E {a['edges']}, F {a['faces']})"
    )
    _say(f"    (V) valid_mask all True at {', '.join(PNAMES)}: {a['valid']}")
    fails = []
    if not mq["recip"] <= RECIPROCITY_BAND:
        fails.append(("i", mq["recip"]))
    if not float(mq["sigma"][0]) <= 1.0 + PASSIVITY_SIGMA_TOLERANCE:
        fails.append(("ii sigma", float(mq["sigma"][0])))
    if not np.all(mq["colpow"] <= 1.0):
        fails.append(("ii colpow", mq["colpow"].tolist()))
    for pid in PORT_IDS:
        ex = a["shares"][pid]
        r = abs(ex["p_src"] - ex["p_vol"] - ex["sheet_field_total"]) / abs(ex["p_src"])
        if not r <= DISCRETE_IDENTITY_RTOL:
            fails.append(("PP", pid, r))
    if rod:
        if not max(mq["mirror"]) <= ADJACENT_SPREAD_BAND:
            fails.append(("M", mq["mirror"]))
        _say(
            f"    (E) tag-402 dofs {a['wall_dofs']}; max |E| "
            + ", ".join(f"{p} {v:.3e}" for p, v in a["wall_max"].items())
            + f" (bound {PEC_WALL_MAX_E:g})"
        )
        if not a["wall_dofs"] > 0 or any(v > PEC_WALL_MAX_E for v in a["wall_max"].values()):
            fails.append(("E", a["wall_dofs"], a["wall_max"]))
    else:
        if not mq["separation"] >= POOLED_SEPARATION_FLOOR:
            fails.append(("S", mq["separation"]))
    if not a["valid"]:
        fails.append(("V",))
    if a["ndofs"] != pred_d:
        fails.append(("D", a["ndofs"], pred_d))
    for k, v in a["areas"].items():
        if not abs(v - NOMINAL_STRIP_AREA) <= STRIP_BAND:
            fails.append(("G2", k, v))
    _say(f"{TAG} {a['label']} anchors failed: {fails if fails else 'none'}")
    return fails


def _ser(a):
    return {
        "label": a["label"],
        "degree": a["degree"],
        "cells": a["cells"],
        "ndofs": a["ndofs"],
        "edges": a["edges"],
        "faces": a["faces"],
        "t_sweep": a["t_sweep"],
        "E": [
            [[[z.real, z.imag] for z in a["E"][k, c, :]] for c in range(3)] for k in range(len(PTS))
        ],
        "sar": {f"{t}|{d}": v for (t, d), v in a["sar"].items()},
        "mass": a["mass"],
        "s": [[[z.real, z.imag] for z in row] for row in a["mq"]["s"]],
    }


def _save(knob, payload):
    if MPI.COMM_WORLD.rank == 0:
        with open(OUT / f"ours_imp3g_{knob}.json", "x") as fh:
            json.dump(payload, fh)
    MPI.COMM_WORLD.Barrier()


def _cplx(arr):
    a = np.array(arr)
    return a[..., 0] + 1j * a[..., 1]


def _census(arm):
    """W1: build serially, print cells, (G1)-(G3)."""
    comm = MPI.COMM_WORLD
    built = _built(arm)
    case = build_four_port_sweep(frequency_hz=STEP2D_FREQUENCY_HZ, built=built, build_only=True)
    msh = case["mesh"]
    ncells = int(case["cells"])
    areas = _sheet_areas(msh, case["facet_tags"], comm)
    _, census = _outer_box_tags(msh, built[2], comm)
    _say(
        f"\n{TAG} W1 {arm}: {ncells} cells, mesh wall {built[3]['mesh_wall_time_s']:.2f} s; "
        f"N_i areas "
        + ", ".join(f"P{i} {a:.12e}" for i, a in areas.items())
        + f"; 402 facets {census['wall']} (exterior {census['wall_exterior']}); "
        f"exterior {census['exterior']}"
    )
    return {"cells": ncells, "areas": areas, "census": census}


def _w1():
    off = _census("off")
    gc.collect()
    rod = _census("rod")
    a_1g = averaging_ball_radius(mass_kg=ONE_GRAM_KG, rho=PHANTOM_RHO_KG_PER_M3)
    rel = off["cells"] / ONE_GRAM_RUNG_CELL_RECORD - 1
    _say(
        f"{TAG} (G0) off cells {off['cells']} vs ragged record {ONE_GRAM_RUNG_CELL_RECORD} "
        f"({rel:+.4%}; CELL_COUNT_BAND {CELL_COUNT_BAND:.2%} shown, printed not asserted; "
        f"predicted ≈ 201 000 / ≈ 252 000); rod cells {rod['cells']}"
    )
    _say(
        f"{TAG} (G4) a_1g {a_1g * 1e3:.4f} mm; x {ROD_CENTRE[0] + a_1g:.6f} < R "
        f"{PHANTOM_RADIUS}; z {0.016 + a_1g:.6f} < H/2 {0.5 * PHANTOM_HEIGHT}"
    )
    _save(
        "build",
        {
            "off": {"cells": off["cells"], "census": off["census"]},
            "rod": {"cells": rod["cells"], "census": rod["census"]},
        },
    )
    assert rod["cells"] < ROD_CELL_GUARD_W2, (rod["cells"], ROD_CELL_GUARD_W2)  # (G1)
    for arm in (off, rod):  # (G2)
        for i, a in arm["areas"].items():
            assert abs(a - NOMINAL_STRIP_AREA) <= STRIP_BAND, (i, a)
    assert rod["census"]["wall"] > 0  # (G3)
    assert rod["census"]["wall_exterior"] == rod["census"]["wall"]
    assert off["census"]["wall"] == 0
    assert ROD_CENTRE[0] + a_1g < PHANTOM_RADIUS  # (G4)
    assert 0.016 + a_1g < 0.5 * PHANTOM_HEIGHT


def _w2_w3(knob):
    comm = MPI.COMM_WORLD
    t0 = time.perf_counter()
    rod = knob == "d2rod"
    arm = "rod" if rod else "off"
    built = _built(arm)
    a = _arm_d(built, f"degree 2 {arm}", rod, 2)
    fails = _identities(a, rod)
    rss = _peak_rss_gib()
    _say(
        f"{TAG} {knob}: cells {a['cells']}, dofs {a['ndofs']}, sweep {a['t_sweep']:.2f} s, "
        f"SARs {a['t_sar']:.2f} s, arm {a['t_arm']:.1f} s, window body "
        f"{time.perf_counter() - t0:.1f} s, summed peak RSS {rss:.2f} GiB at -n {comm.size}"
    )
    payload = _ser(a)
    payload["rss_gib"] = rss
    payload["ranks"] = comm.size
    _save(knob, payload)
    assert not fails, fails


def _w4():
    comm = MPI.COMM_WORLD
    t0 = time.perf_counter()
    d2 = {}
    for arm, knob in (("off", "d2off"), ("rod", "d2rod")):
        d2[arm] = json.loads((OUT / f"ours_imp3g_{knob}.json").read_text())
    fails = []
    off = _arm_d(_built("off"), "degree 1 off", False, 1)
    fails += _identities(off, False)
    for k in ("case",):
        off.pop(k, None)
    gc.collect()
    rod = _arm_d(_built("rod"), "degree 1 rod", True, 1)
    fails += _identities(rod, True)
    kw = _control_c1({"off": off, "rod": rod})
    _say(
        f"{TAG} (Kw) natural wall P1 degree 1: max |E| on 402 dofs {kw['wall_max']:.3e} "
        f"> {PEC_WALL_MAX_E:g}: {kw['wall_max'] > PEC_WALL_MAX_E}; "
        f"ratio {kw['wall_max'] / PEC_WALL_MAX_E:.3e}x"
    )
    if not kw["wall_max"] > PEC_WALL_MAX_E:
        fails.append(("Kw", kw["wall_max"]))
    rss = _peak_rss_gib()
    _say(
        f"{TAG} d1: cells {off['cells']} / {rod['cells']}, dofs {off['ndofs']} / "
        f"{rod['ndofs']}, sweeps {off['t_sweep']:.2f} / {rod['t_sweep']:.2f} s, summed "
        f"peak RSS {rss:.2f} GiB at -n {comm.size}"
    )
    payload = {
        "off": _ser(off),
        "rod": _ser(rod),
        "kw_wall_max": kw["wall_max"],
        "rss_gib": rss,
        "ranks": comm.size,
    }
    _save("d1", payload)

    # ------------------------------------------------------------- readout
    j = 0
    iT1, iF1, iF2 = PNAMES.index("T1"), PNAMES.index("F1"), PNAMES.index("F2")
    E1 = {"off": off["E"], "rod": rod["E"]}
    E2 = {a: _cplx(d2[a]["E"]) for a in ("off", "rod")}
    m1 = {a: np.linalg.norm(E1[a], axis=1) for a in E1}
    m2 = {a: np.linalg.norm(E2[a], axis=1) for a in E2}
    eta1 = m1["rod"] / m1["off"]
    eta2 = m2["rod"] / m2["off"]
    for jj in range(len(PORT_IDS)):
        _say(
            f"{TAG} η_rod P{jj + 1} o1 / o2 / o2/o1-1: "
            + "  ".join(
                f"{n} {eta1[k, jj]:.6f}/{eta2[k, jj]:.6f}/{eta2[k, jj] / eta1[k, jj] - 1:+.6e}"
                for k, n in enumerate(PNAMES)
            )
        )
    pmove = eta2[iT1, j] / eta1[iT1, j] - 1
    _say(
        f"{TAG} p-move η_rod^o2(T1)/η_rod^o1(T1) - 1 (P1) on the 1 g rung {pmove:+.6e} "
        f"({pmove:+.4%}) beside record mesh {RECORD_PMOVE_T1:+.6e} ({RECORD_PMOVE_T1:+.4%})"
    )
    _say(
        f"{TAG} (C-p) |p-move| {abs(pmove):.6e} (predicted >= {PRED_CP_FLOOR:g}, ceiling ≈ 1): "
        f"{'IN' if abs(pmove) >= PRED_CP_FLOOR else 'OUT'}"
    )
    for lab, eta in (("o1", eta1), ("o2", eta2)):
        c2 = abs(eta[iT1, j] - 1) / max(abs(eta[iF1, j] - 1), abs(eta[iF2, j] - 1))
        _say(
            f"{TAG} (C2) {lab}: |η(T1)-1| / max_F |η(F)-1| = {abs(eta[iT1, j] - 1):.6f} / "
            f"{max(abs(eta[iF1, j] - 1), abs(eta[iF2, j] - 1)):.6f} = {c2:.3f} "
            f"(predicted >= {PRED_C2_FLOOR:g}: {'IN' if c2 >= PRED_C2_FLOOR else 'OUT'})"
        )
    rec1 = json.loads((OUT / "ours_imp3p_fields.json").read_text())
    rec2 = json.loads((OUT / "ours_imp3q_fields.json").read_text())
    for lab, rec, eta in (("o1", rec1, eta1), ("o2", rec2, eta2)):
        mr = {a: np.linalg.norm(_cplx(rec["E"][a]), axis=1) for a in ("off", "rod")}
        er = mr["rod"] / mr["off"]
        _say(
            f"{TAG} record mesh -> 1 g rung at {lab}, η_rod(T1) P1: {er[iT1, j]:.6f} -> "
            f"{eta[iT1, j]:.6f}, move {eta[iT1, j] / er[iT1, j] - 1:+.6e} (remesh reading)"
        )
    _say(f"{TAG} record-mesh η_rod^o2(T1) P1 of record {RECORD_ETA_O2_T1}")
    sar2 = {a: d2[a]["sar"] for a in ("off", "rod")}
    for d in DRIVES:
        for tip in TIPS:
            q1 = rod["sar"][(tip, d)] / off["sar"][(tip, d)]
            q2 = sar2["rod"][f"{tip}|{d}"] / sar2["off"][f"{tip}|{d}"]
            _say(
                f"{TAG} Q {d:>3s} {tip:>6s}: o1 {q1:.6f}  o2 {q2:.6f}  Q^o2/Q^o1 - 1 "
                f"{q2 / q1 - 1:+.6e}"
            )
    q1top = rod["sar"][("top", "P1")] / off["sar"][("top", "P1")]
    _say(f"{TAG} Q^o1 P1 top {q1top:.6f} beside step 1's ragged {STEP1_RAGGED_Q_P1_TOP}")
    _say(f"{TAG} degree-2 masses off {d2['off']['mass']} rod {d2['rod']['mass']}")
    cf = np.linalg.norm(
        pec_prolate_spheroid_field(
            PTS - np.array(ROD_CENTRE), semi_major=0.016, semi_minor=0.001, e0=1.0
        ),
        axis=1,
    )
    for lab, eta in (("o1", eta1), ("o2", eta2)):
        _say(
            f"{TAG} closed form |E_cf|/E0 vs η_rod(p)/η_rod(F1) ({lab}, P1), ratio: "
            + "  ".join(
                f"{PNAMES[k]} {cf[k]:.4f}/{eta[k, j] / eta[iF1, j]:.4f}="
                f"{cf[k] / (eta[k, j] / eta[iF1, j]):.3f}"
                for k in (PNAMES.index(n) for n in ("T1", "T2", "T3", "S1"))
            )
        )
    _say(
        f"{TAG} degree-2 prices: off sweep {d2['off']['t_sweep']:.2f} s "
        f"({d2['off']['ndofs']} dofs, {d2['off']['rss_gib']:.2f} GiB at -n "
        f"{d2['off']['ranks']}); rod sweep {d2['rod']['t_sweep']:.2f} s "
        f"({d2['rod']['ndofs']} dofs, {d2['rod']['rss_gib']:.2f} GiB)"
    )
    _say(f"{TAG} W4 body {time.perf_counter() - t0:.1f} s")
    assert not fails, fails


@complex_only
def test_imp3_tip_ratio_1g():
    if KNOB not in KNOBS:
        pytest.skip(f"FEM_EM_IMP3_STEP6 must be one of {KNOBS}")
    if KNOB == "build":
        assert MPI.COMM_WORLD.size == 1, "W1 builds serially (-n 1)"
        _w1()
    elif KNOB in ("d2off", "d2rod"):
        _w2_w3(KNOB)
    else:
        _w4()
