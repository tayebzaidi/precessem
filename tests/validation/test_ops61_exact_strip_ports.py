"""`OPS-61` step 1b (§9 item 148, chain F10a) — the exact f = 0.5 port strip.

On the `ANS-4` x1 rung (``_build(np.zeros(4))``), step 2b's midpoint filter
(``_narrowed_radial``, ``GATED_WIDTH_FRACTION`` = 0.5) selects a *ragged* set
N_i: the mesh does not conform to the strip edges, so step 1a's W1 read
5.835298880e-05 m^2 (fraction 0.521009) on every port, +4.20e-2 relative off
the exact half rectangle (``20260929T110345Z_OPS-61-step1a-W1.log``), and
Palace refuses that set as not aligned with its bounding box.

``MeshGenerator.birdcage_port_domain(leg_sheet_strip_fraction=0.5)`` fragments
each leg sheet against its centre rectangle at azimuth 0 before placement, so
the strip edges are CAD edges and the same filter selects exactly the centre
rectangle.  Anchors (bands imported, never restated):

* (A1) flag on, every F_i = dx*dz within ``SHEET_AREA_BAND``;
* (A2) every N_i = 0.5*dx*dz within 1e-9 (relative), width A/h = 0.5*dx
  within 1e-9;
* (A3) N_i is a subset of F_i and F_i minus N_i = 0.5*dx*dz within 1e-9;
* (A4) flag off, cells = ``STEP2_CELL_COUNT`` within ``CELL_COUNT_BAND``;
* (A5) flag on, `GEO-18`'s box-half volume partition within ``EXACT``;
* (A6) flag-on cell count ``ANS4_EXACT_STRIP_CELL_RECORD`` within
  ``CELL_COUNT_BAND`` (registered after W1; ``None`` means print-only).

Negative control (asserted, same fixture and quantity): the flag-off N_i
misses 0.5*dx*dz by more than 1e-3 relative.

Areas, facet counts and extents are rank-local and reduced before asserting.
"""

from __future__ import annotations

import numpy as np
import pytest
from mpi4py import MPI

from tests.mesh.test_birdcage_conductor_hole import EXACT
from tests.mesh.test_birdcage_leg_offset import (
    SHEET_AREA_BAND,
    _build,
    _projected_extents,
    _sheet_areas,
    _sheet_azimuth_deg,
)
from tests.mesh.test_birdcage_port_sheet_prerequisite import CELL_COUNT_BAND
from tests.mesh.test_birdcage_port_sheets import PORT_LOWER, PORT_UPPER, SHEET_IFACE
from tests.mesh.test_birdcage_port_tags import LEG_COUNT
from tests.mesh.test_birdcage_port_terminals import _interface_area_or_zero
from tests.mesh.test_coil_phantom_conforming import _tag_volume
from tests.validation.test_port_birdcage_leg_offset_sweep import (
    _narrowed_radial,
    _port_frame,
)
from tests.validation.test_port_birdcage_lumped_column import STEP2_CELL_COUNT
from tests.validation.test_port_lumped_narrowed_sheet import GATED_WIDTH_FRACTION

# (A2)/(A3)'s band, from the brief: an exact CAD strip leaves only
# floating-point reassociation in the area reduction.
STRIP_BAND = 1.0e-9

# (A6) — the flag-on cell count, version-tagged: 117 206 measured at -n 1 on
# the complex 0.11 image `main` boots (`20261001T215942Z_OPS-61-F10a-W1.log`
# :3575; ratio 1.009657 to STEP2_CELL_COUNT).  Mesher-dependent — re-record
# under the (1*) licence on an image change, never re-band.
ANS4_EXACT_STRIP_CELL_RECORD = 117206

# The negative control's floor: step 1a W1 read +4.20e-2 on the ragged set.
NEGATIVE_CONTROL_FLOOR = 1.0e-3

PORTS = list(range(1, LEG_COUNT + 1))


def _owned_facet_count(msh, facet_tags, tag, comm):
    fdim = msh.topology.dim - 1
    n_owned = msh.topology.index_map(fdim).size_local
    facets = np.asarray(facet_tags.find(tag), dtype=np.int64)
    return int(comm.allreduce(int(np.count_nonzero(facets < n_owned)), op=MPI.SUM))


def _rung(strip_fraction):
    comm = MPI.COMM_WORLD
    msh, cell_tags, _facets, diag, elapsed = _build(
        np.zeros(LEG_COUNT), leg_sheet_strip_fraction=strip_fraction
    )
    tdim = msh.topology.dim
    msh.topology.create_connectivity(tdim - 1, tdim)
    tags_full, full_areas = _sheet_areas(msh, cell_tags, PORTS, comm)
    dx, _dy, dz = (float(v) for v in diag["port_box_size_m"])
    box_volume = float(np.prod(diag["port_box_size_m"]))

    tags_n = tags_full
    ports = {}
    for i in PORTS:
        tag = SHEET_IFACE + i
        azimuth = _sheet_azimuth_deg(msh, tags_full, tag, comm)
        radial, azimuthal, axial = _port_frame(azimuth)
        spans, centres = _projected_extents(msh, tags_full, tag, comm, [radial, azimuthal, axial])
        ports[i] = {
            "tag": tag,
            "frame": (radial, azimuthal, axial),
            "w_full": float(spans[0]),
            "centre_radial": float(centres[0]),
            "F": float(full_areas[i]),
            "F_facets": _owned_facet_count(msh, tags_full, tag, comm),
        }
    for i in PORTS:
        p = ports[i]
        tags_n = _narrowed_radial(
            msh,
            tags_n,
            p["tag"],
            GATED_WIDTH_FRACTION,
            p["centre_radial"],
            p["frame"][0],
            0.5 * p["w_full"],
        )
    for i in PORTS:
        p = ports[i]
        full = np.asarray(tags_full.find(p["tag"]), dtype=np.int64)
        narrow = np.asarray(tags_n.find(p["tag"]), dtype=np.int64)
        p["subset_local"] = bool(np.all(np.isin(narrow, full)))
        p["N"] = float(_interface_area_or_zero(msh, tags_n, p["tag"], comm))
        p["N_facets"] = _owned_facet_count(msh, tags_n, p["tag"], comm)
        spans_n, _c = _projected_extents(msh, tags_n, p["tag"], comm, list(p["frame"]))
        p["N_span_radial"] = float(spans_n[0])
        p["N_h"] = float(spans_n[2])
        p["halves"] = (
            _tag_volume(msh, cell_tags, PORT_LOWER + i, comm),
            _tag_volume(msh, cell_tags, PORT_UPPER + i, comm),
        )
    subset_all = bool(comm.allreduce(all(ports[i]["subset_local"] for i in PORTS), op=MPI.LAND))
    return {
        "cells": int(msh.topology.index_map(tdim).size_global),
        "dx": dx,
        "dz": dz,
        "box_volume": box_volume,
        "ports": ports,
        "subset_all": subset_all,
        "elapsed": elapsed,
    }


@pytest.fixture(scope="module")
def rungs():
    comm = MPI.COMM_WORLD
    built = {"on": _rung(GATED_WIDTH_FRACTION), "off": _rung(None)}
    if comm.rank == 0:
        for name, r in built.items():
            half = 0.5 * r["dx"] * r["dz"]
            print(
                f"\n[OPS-61 1b] flag {name}: cells={r['cells']} "
                f"(STEP2_CELL_COUNT {STEP2_CELL_COUNT}, ratio "
                f"{r['cells'] / STEP2_CELL_COUNT:.6f}) rung={r['elapsed']:.2f} s "
                f"dx*dz={r['dx'] * r['dz']:.9e} half={half:.9e}",
                flush=True,
            )
            for i in PORTS:
                p = r["ports"][i]
                comp = p["F"] - p["N"]
                print(
                    f"    P{i}: F={p['F']:.9e} ({p['F_facets']} tri) "
                    f"N={p['N']:.9e} ({p['N_facets']} tri, rel "
                    f"{p['N'] / half - 1.0:+.3e}) width A/h={p['N'] / p['N_h']:.9e} "
                    f"N span radial={p['N_span_radial']:.9e} h={p['N_h']:.9e} "
                    f"F\\N={comp:.9e} (rel {comp / half - 1.0:+.3e}) halves sum="
                    f"{sum(p['halves']) / r['box_volume']:.12f}",
                    flush=True,
                )
        print(
            "[OPS-61 1b] predicted (printed, not asserted): flag-on cells "
            "116085 x (1.001-1.01); F_i triangles +10-30 % over flag-off; "
            "N_i facet count equal on all four ports",
            flush=True,
        )
    return built


def test_a1_full_sheet_is_the_rectangle(rungs):
    r = rungs["on"]
    for i in PORTS:
        ratio = r["ports"][i]["F"] / (r["dx"] * r["dz"])
        assert abs(ratio - 1.0) < SHEET_AREA_BAND, f"P{i}: F/dxdz={ratio:.12f}"


def test_a2_strip_is_exactly_half(rungs):
    r = rungs["on"]
    half = 0.5 * r["dx"] * r["dz"]
    for i in PORTS:
        p = r["ports"][i]
        assert abs(p["N"] / half - 1.0) < STRIP_BAND, (
            f"P{i}: N={p['N']:.12e} vs {half:.12e} ({p['N_facets']} tri, "
            f"radial span {p['N_span_radial']:.9e})"
        )
        width = p["N"] / p["N_h"]
        target = GATED_WIDTH_FRACTION * r["dx"]
        assert (
            abs(width / target - 1.0) < STRIP_BAND
        ), f"P{i}: width A/h={width:.12e} vs {target:.12e}"


def test_a3_strip_is_a_subset_and_the_complement_is_half(rungs):
    r = rungs["on"]
    half = 0.5 * r["dx"] * r["dz"]
    assert r["subset_all"], "some N_i is not a subset of F_i"
    for i in PORTS:
        p = r["ports"][i]
        comp = p["F"] - p["N"]
        assert abs(comp / half - 1.0) < STRIP_BAND, f"P{i}: F\\N={comp:.12e} vs {half:.12e}"


def test_a4_default_path_cell_count_unchanged(rungs):
    cells = rungs["off"]["cells"]
    assert (
        abs(cells / STEP2_CELL_COUNT - 1.0) < CELL_COUNT_BAND
    ), f"flag off meshed {cells} cells against {STEP2_CELL_COUNT}"


def test_a5_box_halves_partition_the_box(rungs):
    r = rungs["on"]
    for i in PORTS:
        partition = sum(r["ports"][i]["halves"]) / r["box_volume"]
        assert abs(partition - 1.0) < EXACT, f"P{i}: halves sum {partition:.12f}"


def test_a6_flag_on_cell_record(rungs):
    cells = rungs["on"]["cells"]
    if MPI.COMM_WORLD.rank == 0:
        print(
            f"\n[OPS-61 1b] A6 flag-on cells = {cells} " f"(record {ANS4_EXACT_STRIP_CELL_RECORD})",
            flush=True,
        )
    if ANS4_EXACT_STRIP_CELL_RECORD is None:
        pytest.skip("A6 record not yet registered (W1 measures it)")
    assert abs(cells / ANS4_EXACT_STRIP_CELL_RECORD - 1.0) < CELL_COUNT_BAND, (
        f"flag-on meshed {cells} cells against the record " f"{ANS4_EXACT_STRIP_CELL_RECORD}"
    )


def test_negative_control_flag_off_strip_is_ragged(rungs):
    r = rungs["off"]
    half = 0.5 * r["dx"] * r["dz"]
    for i in PORTS:
        miss = r["ports"][i]["N"] / half - 1.0
        if MPI.COMM_WORLD.rank == 0:
            print(f"[OPS-61 1b] negative control P{i}: miss {miss:+.3e}", flush=True)
        assert (
            abs(miss) > NEGATIVE_CONTROL_FLOOR
        ), f"P{i}: flag-off N_i misses only {miss:+.3e} — the control is not one"


def test_strip_fraction_refusals():
    with pytest.raises(ValueError):
        _build(np.zeros(LEG_COUNT), leg_sheet_strip_fraction=1.0)
