"""`GEO-35` step 1 — the inert leg-port boxes dropped behind a knob.

§9 item 275 (chain U17), brief ``docs/planning/chunks/GEO-35.md``.
``birdcage_port_domain(emit_leg_ports=False)`` on an unbroken-leg layout does
not build the ``leg_count`` floating port boxes at the mid-leg azimuths; the
ring-port ordinals stay ``LEG_COUNT + 1 .. 3·LEG_COUNT``.

Anchors (asserted):
  (K) contract, no mesh — default ``True``; three refusals before gmsh.
  (D) the default build is the `WF-7` record exactly (507 266,
      ``20261005T003402Z_WF-7-H15a-r3-W2.log:10583``), leg tags present.
  (N) knob-on: zero cells on ``PORT_LOWER + i`` / ``PORT_UPPER + i``,
      i = 1..16, allreduced; ``diag["leg_ports_emitted"] is False``.
  (R) knob-on: every ring ordinal 17..48 has both halves with cells, every
      sheet interface has area > 0, ring tag list equals the default's.
  (C) knob-on CAD identities at imported bands (partition over the
      ring-only tag list, terminal ratios, conductor mass).
Negative controls (asserted): (c1) the leg-tag census reads > 0 on all 16
tags of the default build; (c2) on the default build, the ring-only tag list
breaks the partition by more than ``EXACT`` (≈ 16 box volumes / V_total).

Every band imported, none restated. Tier: (K) smoke, mesh test heavy, ``-n 2``.
"""

from __future__ import annotations

import inspect
import time

import gmsh
import numpy as np
import pytest
from mpi4py import MPI

from fem_em_solver.io.mesh import MeshGenerator, _interface_facet_tags

from tests.mesh.helpers import global_cell_tag_set
from tests.mesh.test_birdcage_conductor_sizing import CAD_MASS_GATE
from tests.mesh.test_birdcage_port_sheet_prerequisite import CELL_COUNT_BAND
from tests.mesh.test_birdcage_port_sheets import PORT_LOWER, PORT_UPPER, SHEET_IFACE
from tests.mesh.test_birdcage_port_terminals import (
    CONDUCTOR_IFACE,
    _interface_area_or_zero,
)
from tests.mesh.test_birdcage_ring_gaps import EXACT, _port_boundary_partition
from tests.mesh.test_coil_phantom_conforming import _tag_volume, _total_volume
from tests.validation.test_birdcage_f_human_rung import (
    F_HUMAN_BRANCH_A_CELL_RECORD,
    F_HUMAN_BRANCH_B_CELL_RECORD,
    F_HUMAN_RING_RADIUS,
    LEG_COUNT,
    TERMINAL_RATIO_BAND,
    _build_rung,
    _params,
    _ring_ports,
)

# `WF-7`'s degree-1 F-human cell record on the 0.11 image, same kwargs
# (``20261005T003402Z_WF-7-H15a-r3-W2.log:10583``). (D) asserts it exactly:
# the default path must be digit-identical.
WF7_CELL_RECORD = 507266

LEG_IDS = list(range(1, LEG_COUNT + 1))


def _wf7_kwargs() -> dict:
    """`WF-7`'s build (``test_wf7_f_human_ring_matrix.py:451–452``)."""
    kwargs = _params(F_HUMAN_RING_RADIUS, scale_sizing=False)
    kwargs["ring_sheet_orientation"] = "longitudinal"
    return kwargs


def _owned_census(msh, cell_tags, tags, comm) -> dict[int, int]:
    """Global count of *owned* cells per tag (ghosts excluded), allreduced."""
    n_local = msh.topology.index_map(3).size_local
    owned = cell_tags.indices < n_local
    vals = cell_tags.values[owned]
    local = np.array([np.count_nonzero(vals == t) for t in tags], dtype=np.int64)
    total = np.zeros_like(local)
    comm.Allreduce(local, total, op=MPI.SUM)
    return {t: int(c) for t, c in zip(tags, total)}


def _build(comm, **extra) -> dict:
    kwargs = _wf7_kwargs()
    kwargs.update(extra)
    started = time.perf_counter()
    msh, cells, _facets, diag = MeshGenerator.birdcage_port_domain(
        comm=comm, return_diagnostics=True, **kwargs
    )
    build_s = time.perf_counter() - started
    ring_ports = _ring_ports()
    ring_cell_tags = {i: (PORT_LOWER + i, PORT_UPPER + i) for i in ring_ports}
    ring_tag_list = sorted(t for v in ring_cell_tags.values() for t in v)

    v_total = _total_volume(msh, comm)
    base_vol = sum(_tag_volume(msh, cells, t, comm) for t in (1, 2, 3))
    ring_vol = sum(_tag_volume(msh, cells, t, comm) for t in ring_tag_list)
    leg_vol = sum(
        _tag_volume(msh, cells, t, comm)
        for i in LEG_IDS
        for t in (PORT_LOWER + i, PORT_UPPER + i)
    )

    _counts, areas = _port_boundary_partition(msh, cells, comm, ring_cell_tags)
    terminal_analytic = float(diag["ring_port_layout"]["ring_terminal_area_m2"])
    ratios = np.array(
        [areas[CONDUCTOR_IFACE + i] / terminal_analytic for i in ring_ports]
    )

    sheet_tags = _interface_facet_tags(
        msh, cells,
        {SHEET_IFACE + i: (PORT_LOWER + i, PORT_UPPER + i) for i in ring_ports},
    )
    sheet_areas = np.array(
        [_interface_area_or_zero(msh, sheet_tags, SHEET_IFACE + i, comm)
         for i in ring_ports]
    )

    leg_tags = [PORT_LOWER + i for i in LEG_IDS] + [PORT_UPPER + i for i in LEG_IDS]
    return {
        "kwargs": kwargs,
        "diag": diag,
        "n_cells": int(msh.topology.index_map(3).size_global),
        "build_s": build_s,
        "v_total": v_total,
        "partition_ring_only": (base_vol + ring_vol) / v_total,
        "partition_with_legs": (base_vol + ring_vol + leg_vol) / v_total,
        "terminal_ratios": ratios,
        "sheet_areas": sheet_areas,
        "cad_ratio": _tag_volume(msh, cells, 1, comm)
        / diag["cad_mass_by_group"]["conductor"],
        "leg_census": _owned_census(msh, cells, leg_tags, comm),
        "ring_census": _owned_census(msh, cells, ring_tag_list, comm),
        "ring_tags_present": sorted(
            t for t in global_cell_tag_set(msh, cells)
            if PORT_LOWER + LEG_COUNT < t <= PORT_LOWER + 3 * LEG_COUNT
            or PORT_UPPER + LEG_COUNT < t <= PORT_UPPER + 3 * LEG_COUNT
        ),
        "ring_tag_list": ring_tag_list,
    }


def test_geo35_knob_contract():
    """(K): default ``True``; ``False`` refused before any gmsh state exists."""
    sig = inspect.signature(MeshGenerator.birdcage_port_domain)
    param = sig.parameters["emit_leg_ports"]
    assert param.default is True
    assert param.kind is inspect.Parameter.KEYWORD_ONLY

    base = _wf7_kwargs()
    refusals = {
        "leg_gap_length": dict(base, leg_gap_length=0.01, ring_gap_length=None,
                               ring_sheet_orientation="transverse"),
        "as_hole": dict(base, as_hole=True),
        "symmetry_quadrant": dict(base, symmetry_quadrant=True, as_hole=False),
    }
    for name, kwargs in refusals.items():
        assert not gmsh.isInitialized()
        with pytest.raises(ValueError, match="emit_leg_ports"):
            MeshGenerator.birdcage_port_domain(
                comm=MPI.COMM_WORLD, emit_leg_ports=False, **kwargs
            )
        assert not gmsh.isInitialized(), f"{name}: gmsh state existed at refusal"
        if MPI.COMM_WORLD.rank == 0:
            print(f"\n[GEO-35 (K)] emit_leg_ports=False with {name}: ValueError "
                  "before gmsh", flush=True)


def test_geo35_f_human_builder_contract():
    """(K, step 2a): ``_f_human_kwargs`` resolution, the pin and the hole exemption.

    Negative control (c1), asserted: "key present, False" (solid default) is
    separated from "no key" (pinned boxed, or ``as_hole``).
    """
    from tests.validation.test_wf7_f_human_b1_quadrature import _f_human_kwargs

    default = _f_human_kwargs()
    assert default == dict(_wf7_kwargs(), emit_leg_ports=False)
    assert default["emit_leg_ports"] is False

    pinned = _f_human_kwargs(emit_leg_ports=True)
    assert "emit_leg_ports" not in pinned
    assert pinned == _wf7_kwargs()
    pinned_res = _f_human_kwargs(phantom_resolution=0.0075, emit_leg_ports=True)
    assert "emit_leg_ports" not in pinned_res
    assert pinned_res == dict(_wf7_kwargs(), phantom_resolution=0.0075)

    hole = _f_human_kwargs(as_hole=True)
    assert "emit_leg_ports" not in hole
    assert hole == dict(_wf7_kwargs(), as_hole=True)

    assert not gmsh.isInitialized()
    with pytest.raises(ValueError, match="emit_leg_ports"):
        _f_human_kwargs(as_hole=True, emit_leg_ports=False)
    assert not gmsh.isInitialized()
    if MPI.COMM_WORLD.rank == 0:
        print("\n[GEO-35 s2a (K)] default drops leg boxes (key False); pinned and "
              "as_hole add no key; as_hole+False refused before gmsh", flush=True)


def test_geo35_f_human_drops_leg_ports_keeps_identities():
    """(D), (N), (R), (C), (c1), (c2) on the `WF-7` F-human fixture."""
    comm = MPI.COMM_WORLD
    default = _build(comm)
    knob = _build(comm, emit_leg_ports=False)
    lo, hi = TERMINAL_RATIO_BAND

    box = default["kwargs"]["port_box_size"]
    c2_predicted = LEG_COUNT * float(np.prod(box)) / default["v_total"]
    c2_measured = 1.0 - default["partition_ring_only"]
    if comm.rank == 0:
        for name, b in (("default (knob off)", default), ("emit_leg_ports=False", knob)):
            r = b["terminal_ratios"]
            print(
                f"\n[GEO-35 step 1] {name}"
                f"\n  cells                 {b['n_cells']}  vs WF-7 {WF7_CELL_RECORD} "
                f"(rel {b['n_cells'] / WF7_CELL_RECORD - 1.0:+.4e})  vs GEO-25 "
                f"{F_HUMAN_BRANCH_B_CELL_RECORD} "
                f"(rel {b['n_cells'] / F_HUMAN_BRANCH_B_CELL_RECORD - 1.0:+.4e}, "
                f"band {CELL_COUNT_BAND})"
                f"\n  build_s               {b['build_s']:.2f}"
                f"\n  leg_ports_emitted     {b['diag']['leg_ports_emitted']}"
                f"\n  partition ring-only   {b['partition_ring_only']:.12f}  "
                f"with legs {b['partition_with_legs']:.12f}  (EXACT {EXACT})"
                f"\n  terminal ratio        min {r.min():.9f} max {r.max():.9f} "
                f"n {len(r)}"
                f"\n  sheet area            min {b['sheet_areas'].min():.6e} m^2"
                f"\n  meshed/CAD conductor  {b['cad_ratio']:.6f} (gate {CAD_MASS_GATE})"
                f"\n  leg census            min {min(b['leg_census'].values())} "
                f"max {max(b['leg_census'].values())}"
                f"\n  ring census           min {min(b['ring_census'].values())}",
                flush=True,
            )
        print(
            f"\n[GEO-35 (c2)] default build, ring-only tag list: partition deficit "
            f"{c2_measured:.6e} (predicted 16 box volumes / V_total = "
            f"{c2_predicted:.6e}; ratio {c2_measured / c2_predicted:.6f}); "
            f"{c2_measured / EXACT:.3e}x EXACT",
            flush=True,
        )

    # (D) the default path is digit-identical to the `WF-7` record.
    assert default["n_cells"] == WF7_CELL_RECORD, default["n_cells"]
    assert abs(default["n_cells"] / F_HUMAN_BRANCH_B_CELL_RECORD - 1.0) < CELL_COUNT_BAND
    assert default["diag"]["leg_ports_emitted"] is True
    # (c1) the leg-tag detector sees all 16 boxes on the default build.
    assert all(default["leg_census"][PORT_LOWER + i] > 0 for i in LEG_IDS), (
        default["leg_census"]
    )
    # (c2) the partition comparator sees 16 missing boxes.
    assert c2_measured > EXACT, c2_measured
    assert abs(default["partition_with_legs"] - 1.0) < EXACT

    # (N) no leg-port tags on the knob-on build.
    assert knob["diag"]["leg_ports_emitted"] is False
    assert all(c == 0 for c in knob["leg_census"].values()), knob["leg_census"]

    # (R) ring-port census unchanged.
    assert all(c > 0 for c in knob["ring_census"].values()), knob["ring_census"]
    assert len(knob["sheet_areas"]) == 2 * LEG_COUNT
    assert np.all(knob["sheet_areas"] > 0.0), knob["sheet_areas"]
    assert knob["ring_tags_present"] == default["ring_tags_present"]
    assert knob["ring_tags_present"] == knob["ring_tag_list"]

    # (C) CAD identities on the knob-on build, imported bands.
    assert abs(knob["partition_ring_only"] - 1.0) < EXACT, knob["partition_ring_only"]
    r = knob["terminal_ratios"]
    assert len(r) == 2 * LEG_COUNT
    assert r.min() >= lo and r.max() <= hi, (r.min(), r.max())
    assert knob["cad_ratio"] >= CAD_MASS_GATE, knob["cad_ratio"]


# Consumers of the boxed `F_HUMAN_BRANCH_B_CELL_RECORD` (the longitudinal
# `WF-7` build, 507 266 cells, asserted against it at `CELL_COUNT_BAND`).
B_RECORD_CONSUMERS = (
    "test_wf7_f_human_ring_matrix.py:1092",
    "test_wf12_f_human_hole.py:380",
    "test_geo35_leg_port_boxes.py:243",
)
# Step 1's measured knob-on move on the longitudinal build (−0.899 %), the
# prediction for branch B here; printed, never asserted (rule (e)).
PREDICTED_B_MOVE = -0.009


def test_geo35_transverse_rung_measure():
    """`GEO-35` step 2b-i part 2 (§9 item 357): `GEO-25` rung, knob on.

    Branches B and A of the transverse `GEO-25` 0.15 m rung built with
    ``emit_leg_ports=False``. Asserted (imported bands): partition within
    ``EXACT`` over the ring-only tag list, 32 terminal ratios in
    ``TERMINAL_RATIO_BAND``, the ``CAD_MASS_GATE`` separation (B above, A
    below), (N′) zero leg-tag cells, and each knob-on count ≠ its boxed
    record (negative control: the switch took effect). Printed only (10-10
    03:00 ruling (3)): each count's distance from its landed record against
    ``CELL_COUNT_BAND``, B's predicted ≈ −0.9 % move, and the distance of
    the boxed 507 266 consumers from the knob-on B count. Pre-registration
    for 2b-ii; nothing re-recorded.
    """
    comm = MPI.COMM_WORLD
    leg_tags = [PORT_LOWER + i for i in LEG_IDS] + [PORT_UPPER + i for i in LEG_IDS]
    records = {False: F_HUMAN_BRANCH_B_CELL_RECORD, True: F_HUMAN_BRANCH_A_CELL_RECORD}
    names = {False: "branch B (fixed absolute sizing)", True: "branch A (scaled — CONTROL)"}
    lo, hi = TERMINAL_RATIO_BAND

    rungs = {}
    for scale_sizing in (False, True):
        rung = _build_rung(
            F_HUMAN_RING_RADIUS, scale_sizing=scale_sizing, keep_mesh=True,
            emit_leg_ports=False,
        )
        rung["leg_census"] = _owned_census(
            rung["mesh"], rung["cell_tags"], leg_tags, comm
        )
        rungs[scale_sizing] = rung
        rec = records[scale_sizing]
        rel = rung["n_cells"] / rec - 1.0
        r = rung["terminal_ratios"]
        if comm.rank == 0:
            print(
                f"\n[GEO-35 s2b-i p2] {names[scale_sizing]}, emit_leg_ports=False"
                f"\n  cells                 {rung['n_cells']} vs boxed record {rec} "
                f"rel {rel:+.4e}  in CELL_COUNT_BAND {CELL_COUNT_BAND}: "
                f"{abs(rel) < CELL_COUNT_BAND} (printed, ruling (3))"
                f"\n  build_s               {rung['build_s']:.2f}  mesh_wall "
                f"{rung['mesh_wall_time_s']:.2f}"
                f"\n  leg_ports_emitted     {rung['diag']['leg_ports_emitted']}"
                f"\n  partition ring-only   {rung['partition']:.12f} (EXACT {EXACT})"
                f"\n  terminal ratio        min {r.min():.9f} max {r.max():.9f} n {len(r)}"
                f"\n  meshed/CAD conductor  {rung['cad_ratio']:.6f} (gate {CAD_MASS_GATE})"
                f"\n  (N') leg census       total {sum(rung['leg_census'].values())} "
                f"over {len(leg_tags)} tags",
                flush=True,
            )

    b, a = rungs[False], rungs[True]
    if comm.rank == 0:
        print(
            f"\n[GEO-35 s2b-i p2] branch B move {b['n_cells'] / records[False] - 1.0:+.4%}"
            f" (predicted ≈ {PREDICTED_B_MOVE:+.1%}, printed); branch A move "
            f"{a['n_cells'] / records[True] - 1.0:+.4%} (unpredicted)",
            flush=True,
        )
        for site in B_RECORD_CONSUMERS:
            d = WF7_CELL_RECORD / b["n_cells"] - 1.0
            print(
                f"  consumer {site}: boxed {WF7_CELL_RECORD} vs knob-on B "
                f"{b['n_cells']} rel {d:+.4%} (band {CELL_COUNT_BAND}, inside: "
                f"{abs(d) < CELL_COUNT_BAND}; predicted ≈ +1.4 % outside)",
                flush=True,
            )
        print(
            f"  CAD separation: A {a['cad_ratio']:.6f} < gate {CAD_MASS_GATE} <= "
            f"B {b['cad_ratio']:.6f}",
            flush=True,
        )

    for scale_sizing, rung in rungs.items():
        assert rung["diag"]["leg_ports_emitted"] is False
        # (N′) no leg-port cells, allreduced over owned cells.
        assert all(c == 0 for c in rung["leg_census"].values()), rung["leg_census"]
        assert abs(rung["partition"] - 1.0) < EXACT, rung["partition"]
        r = rung["terminal_ratios"]
        assert len(r) == 2 * LEG_COUNT
        assert r.min() >= lo and r.max() <= hi, (r.min(), r.max())
        # Negative control: the switch took effect on the count.
        assert rung["n_cells"] != records[scale_sizing], rung["n_cells"]
    assert a["cad_ratio"] < CAD_MASS_GATE <= b["cad_ratio"], (
        a["cad_ratio"], b["cad_ratio"],
    )
