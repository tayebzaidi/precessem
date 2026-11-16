"""`PORT-2` step 1 (§9 item 247, lane 1) — the port tagging contract as a
round trip on the exact-strip x1 birdcage.

``src/fem_em_solver/ports/definitions.py`` states the contract (facet tags
primary; positive/negative terminal tag per port) and validates it on tag
*sets* only.  ``tests/ports/test_port_definition.py`` exercises it on literal
integers.  This module reads a real mesh's tags back through it: on the
`OPS-61` exact-strip fixture (``leg_sheet_strip_fraction =
GATED_WIDTH_FRACTION``), the ``PortDefinition`` objects built exactly as
``_four_port_rung`` builds them (``test_port_birdcage_leg_offset_sweep.py``
:380-405) must round-trip:

* (T0) premise: 0.5*dx*dz = ``SPEC_AREA_M2`` within ``STRIP_BAND``; cells vs
  ``ANS4_EXACT_STRIP_CELL_RECORD`` within ``CELL_COUNT_BAND``;
* (T1) ``validate_required_port_tags_exist`` passes on the globally reduced
  facet + cell tag set, positive tags among facet tags, negative among cell;
* (T2) every port's owned-facet count n_i > 0 and equals the size of the union
  over ranks of the global indices of every local facet carrying the tag
  (integer identity: a ghost double-count or a lost facet fails it);
* (T3) every area = ``SPEC_AREA_M2`` within ``STRIP_BAND``;
* (T4) port_id <-> tag bijection, spec facet_tag == positive_tag, quarter-turn
  index of P_i's set == i - 1 (integer, no band), ordering and calibration
  checks pass.

Negative controls (asserted, every rank): (C1) a dropped tag raises "missing
required port tags"; (C2) a swapped pair passes existence and area but fails
the round trip naming both ports.  (c3) — the rank-local read `OPS-14` paid
for — is printed only.

No solve; the only form is ``_interface_area_or_zero``'s ``dS``.
"""

from __future__ import annotations

import re

import dolfinx
import numpy as np
import pytest
from mpi4py import MPI

from fem_em_solver.ports.definitions import (
    PortDefinition,
    run_port_calibration_checks,
    validate_port_ordering,
    validate_required_port_tags_exist,
)
from fem_em_solver.ports.lumped import LumpedSheetPortSpec
from tests.complex_mode import complex_only
from tests.mesh.test_birdcage_leg_offset import (
    _build,
    _projected_extents,
    _sheet_areas,
    _sheet_azimuth_deg,
)
from tests.mesh.test_birdcage_port_sheet_prerequisite import CELL_COUNT_BAND
from tests.mesh.test_birdcage_port_sheets import SHEET_IFACE
from tests.mesh.test_birdcage_port_tags import LEG_COUNT
from tests.mesh.test_birdcage_port_terminals import _interface_area_or_zero
from tests.validation.test_ans4_exact_strip_records import SPEC_AREA_M2
from tests.validation.test_ops61_exact_strip_ports import (
    ANS4_EXACT_STRIP_CELL_RECORD,
    STRIP_BAND,
    _owned_facet_count,
)
from tests.validation.test_port_birdcage_four_port import TERMINATED_PORT_IMPEDANCE_OHM
from tests.validation.test_port_birdcage_leg_offset_sweep import (
    _narrowed_radial,
    _port_frame,
)
from tests.validation.test_port_birdcage_lumped_column import CONDUCTOR_CELL_TAG
from tests.validation.test_port_lumped_narrowed_sheet import GATED_WIDTH_FRACTION
from tests.validation.test_port_package_sparameters import REFERENCE_IMPEDANCE_OHM

PORTS = list(range(1, LEG_COUNT + 1))

# Printed beside n_i, never asserted: the mesher's record, 40 triangles per
# port at -n 1 and -n 2 (`20261001T215942Z_OPS-61-F10a-W1.log`:3576-3579;
# `20261001T220051Z_OPS-61-F10a-W2.log`:3591-3594).
RECORDED_STRIP_TRIANGLES = 40


def _reduced_tag_set(values, comm):
    """Sorted union over ranks of a meshtags' local values (excitation.py idiom)."""
    local = np.unique(np.asarray(values, dtype=np.int64)).tolist()
    return sorted({int(t) for part in comm.allgather(local) for t in part})


def _global_facet_union_size(msh, facet_tags, tag, comm):
    """|union over ranks of local_to_global(every local facet carrying tag)|."""
    fdim = msh.topology.dim - 1
    local = np.asarray(facet_tags.find(tag), dtype=np.int32)
    glob = msh.topology.index_map(fdim).local_to_global(local).tolist()
    return len({int(g) for part in comm.allgather(glob) for g in part})


def _relabel(msh, facet_tags, mapping=None, drop=None):
    """Rebuild facet meshtags with values relabelled / one tag dropped."""
    fdim = msh.topology.dim - 1
    idx = np.asarray(facet_tags.indices, dtype=np.int32)
    val = np.asarray(facet_tags.values, dtype=np.int32).copy()
    keep = np.ones(idx.size, dtype=bool)
    if drop is not None:
        keep = val != int(drop)
    new_val = val.copy()
    if mapping:
        for a, b in mapping.items():
            new_val[val == int(a)] = int(b)
    idx, new_val = idx[keep], new_val[keep]
    order = np.argsort(idx, kind="stable")
    return dolfinx.mesh.meshtags(msh, fdim, idx[order], new_val[order])


def _roundtrip(port_defs, specs, msh, tags, comm):
    """Per-port tag/geometry records; raise naming every mismatched port.

    Every collective runs on every rank before any raise.
    """
    spec_by_id = {s.port_id: s for s in specs}
    records = {}
    for idx, pd in enumerate(port_defs):
        azimuth = _sheet_azimuth_deg(msh, tags, int(pd.positive_tag), comm)
        records[pd.port_id] = {
            "leg": idx + 1,
            "azimuth_deg": azimuth,
            "quarter": int(round(azimuth / 90.0)) % 4,
            "spec_tag": int(spec_by_id[pd.port_id].facet_tag),
        }
    bad = []
    for pd in port_defs:
        r = records[pd.port_id]
        if r["spec_tag"] != int(pd.positive_tag) or r["quarter"] != r["leg"] - 1:
            bad.append(
                f"{pd.port_id} (tag {pd.positive_tag}, spec {r['spec_tag']}, "
                f"quarter {r['quarter']} vs leg index {r['leg'] - 1}, "
                f"az {r['azimuth_deg']:.4f} deg)"
            )
    if bad:
        raise ValueError("port tag/geometry mismatch: " + "; ".join(bad))
    return records


@pytest.fixture(scope="module")
def rung():
    comm = MPI.COMM_WORLD
    msh, cell_tags, _facets, diag, elapsed = _build(
        np.zeros(LEG_COUNT), leg_sheet_strip_fraction=GATED_WIDTH_FRACTION
    )
    tdim = msh.topology.dim
    msh.topology.create_connectivity(tdim - 1, tdim)
    msh.topology.create_entity_permutations()
    tags_f, _full = _sheet_areas(msh, cell_tags, PORTS, comm)
    dx, _dy, dz = (float(v) for v in diag["port_box_size_m"])

    geometry = {}
    for i in PORTS:
        tag = SHEET_IFACE + i
        azimuth = _sheet_azimuth_deg(msh, tags_f, tag, comm)
        radial, azimuthal, axial = _port_frame(azimuth)
        spans, centres = _projected_extents(
            msh, tags_f, tag, comm, [radial, azimuthal, axial]
        )
        geometry[i] = {
            "tag": tag,
            "azimuth_deg": float(azimuth),
            "radial": radial,
            "w_full": float(spans[0]),
            "centre_radial": float(centres[0]),
        }
    for i in PORTS:
        g = geometry[i]
        tags_f = _narrowed_radial(
            msh,
            tags_f,
            g["tag"],
            GATED_WIDTH_FRACTION,
            g["centre_radial"],
            g["radial"],
            0.5 * g["w_full"],
        )

    sheets = []
    for i in PORTS:
        g = geometry[i]
        area = float(_interface_area_or_zero(msh, tags_f, g["tag"], comm))
        radial, azimuthal, axial = _port_frame(g["azimuth_deg"])
        spans, _c = _projected_extents(
            msh, tags_f, g["tag"], comm, [radial, azimuthal, axial]
        )
        h = float(spans[2])
        sheets.append({**g, "area": area, "h": h, "w": area / h if h > 0 else 0.0})

    # Exactly as `_four_port_rung` (:380-405), kappa None; no solve.
    specs = [
        LumpedSheetPortSpec(
            port_id=f"P{s['tag'] - SHEET_IFACE}",
            facet_tag=int(s["tag"]),
            port_impedance_ohm=TERMINATED_PORT_IMPEDANCE_OHM,
            gap_height_m=s["h"],
            sheet_width_m=s["w"],
            drive_direction=(0.0, 0.0, 1.0),
            drive_voltage_v=1.0 + 0.0j,
            interior=True,
        )
        for s in sheets
    ]
    port_defs = [
        PortDefinition(
            port_id=f"P{s['tag'] - SHEET_IFACE}",
            positive_tag=int(s["tag"]),
            negative_tag=CONDUCTOR_CELL_TAG,
            orientation="leg_gap_axial_plus_z",
            z0_ohm=REFERENCE_IMPEDANCE_OHM,
        )
        for s in sheets
    ]

    owned = {i: _owned_facet_count(msh, tags_f, SHEET_IFACE + i, comm) for i in PORTS}
    union = {
        i: _global_facet_union_size(msh, tags_f, SHEET_IFACE + i, comm) for i in PORTS
    }
    facet_set = _reduced_tag_set(tags_f.values, comm)
    cell_set = _reduced_tag_set(cell_tags.values, comm)

    # (c3) printed only: does each rank's *local* tag set alone hold all four
    # sheet tags?  No validator call here; gathered, then printed.
    local_vals = set(np.unique(np.asarray(tags_f.values)).tolist())
    local_has_all = all((SHEET_IFACE + i) in local_vals for i in PORTS)
    c3 = comm.allgather(bool(local_has_all))

    rec = {
        "comm": comm,
        "msh": msh,
        "cell_tags": cell_tags,
        "tags": tags_f,
        "dx": dx,
        "dz": dz,
        "cells": int(msh.topology.index_map(tdim).size_global),
        "elapsed": elapsed,
        "sheets": sheets,
        "specs": specs,
        "port_defs": port_defs,
        "owned": owned,
        "union": union,
        "facet_set": facet_set,
        "cell_set": cell_set,
        "c3": c3,
    }
    if comm.rank == 0:
        print(
            f"\n[PORT-2 step1] ranks {comm.size}; rung build {elapsed:.2f} s; "
            f"cells {rec['cells']} (record {ANS4_EXACT_STRIP_CELL_RECORD}, "
            f"rel {rec['cells'] / ANS4_EXACT_STRIP_CELL_RECORD - 1:+.3e})"
        )
        print(f"[PORT-2 step1] facet tag set {facet_set}")
        print(f"[PORT-2 step1] cell tag set {cell_set}")
        for s in sheets:
            i = s["tag"] - SHEET_IFACE
            print(
                f"[PORT-2 step1] P{i}: tag {s['tag']} az {s['azimuth_deg']:.4f} deg "
                f"quarter {int(round(s['azimuth_deg'] / 90.0)) % 4}; "
                f"A {s['area']:.12e} (rel {s['area'] / SPEC_AREA_M2 - 1:+.3e}); "
                f"n_owned {owned[i]} union {union[i]} "
                f"(recorded {RECORDED_STRIP_TRIANGLES}, printed only)"
            )
        print(f"[PORT-2 step1] (c3) per-rank local set holds all four sheet tags: {c3}")
    return rec


@complex_only
def test_t0_premise(rung):
    half = 0.5 * rung["dx"] * rung["dz"]
    assert abs(half / SPEC_AREA_M2 - 1.0) <= STRIP_BAND, (half, SPEC_AREA_M2)
    assert abs(rung["cells"] / ANS4_EXACT_STRIP_CELL_RECORD - 1.0) < CELL_COUNT_BAND


@complex_only
def test_t1_existence_and_dimension(rung):
    facet_set, cell_set = rung["facet_set"], rung["cell_set"]
    validate_required_port_tags_exist(rung["port_defs"], set(facet_set) | set(cell_set))
    for pd in rung["port_defs"]:
        assert pd.positive_tag in facet_set, (pd.port_id, pd.positive_tag)
        assert pd.negative_tag in cell_set, (pd.port_id, pd.negative_tag)


@complex_only
def test_t2_owned_count_identity(rung):
    for i in PORTS:
        n_i, u_i = rung["owned"][i], rung["union"][i]
        assert n_i > 0, f"P{i}: no owned facets"
        assert n_i == u_i, f"P{i}: owned {n_i} != global union {u_i}"


@complex_only
def test_t3_area(rung):
    for s in rung["sheets"]:
        assert abs(s["area"] / SPEC_AREA_M2 - 1.0) <= STRIP_BAND, (s["tag"], s["area"])


@complex_only
def test_t4_bijection_and_placement(rung):
    comm, msh, tags = rung["comm"], rung["msh"], rung["tags"]
    port_defs, specs = rung["port_defs"], rung["specs"]
    positives = [pd.positive_tag for pd in port_defs]
    assert len(set(positives)) == len(positives)
    assert set(positives) == {SHEET_IFACE + i for i in PORTS}
    assert set(positives) <= set(rung["facet_set"])
    assert [(s.port_id, s.facet_tag) for s in specs] == [
        (pd.port_id, pd.positive_tag) for pd in port_defs
    ]
    records = _roundtrip(port_defs, specs, msh, tags, comm)
    for i, pd in enumerate(port_defs, start=1):
        assert records[pd.port_id]["quarter"] == i - 1
    validate_port_ordering(port_defs, expected_port_order=["P1", "P2", "P3", "P4"])
    areas = {f"P{s['tag'] - SHEET_IFACE}": s["area"] for s in rung["sheets"]}
    checks = run_port_calibration_checks(
        port_defs, expected_port_order=["P1", "P2", "P3", "P4"], port_face_areas=areas
    )
    if comm.rank == 0:
        print(
            f"\n[PORT-2 step1] (T4) quarters "
            f"{[records[pd.port_id]['quarter'] for pd in port_defs]}; "
            f"calibration max face-area ratio {checks.max_face_area_ratio:.15f}"
        )


@complex_only
def test_c1_missing_tag_raises(rung):
    comm, msh = rung["comm"], rung["msh"]
    dropped = SHEET_IFACE + 3
    tags_c1 = _relabel(msh, rung["tags"], drop=dropped)
    facet_set = _reduced_tag_set(tags_c1.values, comm)
    assert dropped not in facet_set
    with pytest.raises(
        ValueError, match=re.escape(f"missing required port tags: [{dropped}]")
    ):
        validate_required_port_tags_exist(
            rung["port_defs"], set(facet_set) | set(rung["cell_set"])
        )


@complex_only
def test_c2_swap_passes_existence_fails_roundtrip(rung):
    comm, msh = rung["comm"], rung["msh"]
    a, b = SHEET_IFACE + 2, SHEET_IFACE + 4
    tags_c2 = _relabel(msh, rung["tags"], mapping={a: b, b: a})
    facet_set = _reduced_tag_set(tags_c2.values, comm)
    validate_required_port_tags_exist(
        rung["port_defs"], set(facet_set) | set(rung["cell_set"])
    )
    areas = {i: float(_interface_area_or_zero(msh, tags_c2, SHEET_IFACE + i, comm)) for i in PORTS}
    for i in PORTS:
        assert abs(areas[i] / SPEC_AREA_M2 - 1.0) <= STRIP_BAND, (i, areas[i])
    with pytest.raises(ValueError, match="port tag/geometry mismatch") as excinfo:
        _roundtrip(rung["port_defs"], rung["specs"], msh, tags_c2, comm)
    msg = str(excinfo.value)
    assert "P2 " in msg and "P4 " in msg, msg
    assert "P1 " not in msg and "P3 " not in msg, msg
    if comm.rank == 0:
        print(f"\n[PORT-2 step1] (C2) raised: {msg}")
