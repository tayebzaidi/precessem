"""`GEO-37` step 2 (§9 item 267, chain U11): asserting reader for the native
interior-tag ``dS`` probe.

Reads ``sf_n{1,2,4}.json`` (ghost mode ``shared_facet``, asserted) and
``none_n{1,2,4}.json`` (ghost mode ``none``, the 2026-08 configuration,
printed only) written by ``tests/mesh/probe_geo37_interior_dS.py`` into
``$GEO37_DS_RECORD_DIR``.  Skipped when that variable is unset.

Anchors (every width, ``shared_facet``): (D1) no hang; (D2) I1 = 0.48;
(D3) I2 = 2.88, I3 = 1.92; (D4) instrument I4 = 6 -- all within the imported
``STRIP_BAND``; (D5) the port-tag validator raises on an absent tag and passes
on 201.  Controls: (c1) ``dS(999)`` assembles to exactly 0.0 (asserted -- why
(D5) needs the validator); (c2) one planted dropped facet leaves the band
(asserted); (c3) the unghosted leg is printed, never asserted.
"""

import json
import os

import pytest

from fem_em_solver.ports.definitions import (
    PortDefinition,
    validate_required_port_tags_exist,
)
from tests.validation.test_ops61_exact_strip_ports import STRIP_BAND

RECORD_DIR = os.environ.get("GEO37_DS_RECORD_DIR")
pytestmark = pytest.mark.skipif(
    RECORD_DIR is None, reason="GEO37_DS_RECORD_DIR unset (probe records not produced)"
)

WIDTHS = (1, 2, 4)
I1_CF = 0.8 * 0.6  # rectangle area
I2_CF = I1_CF * (2 + 3 * 0.5 + 5 * 0.5)  # g = 2 + 3x + 5y, centroid (0.5, 0.5)
I3_CF = I1_CF * (1 + 3)  # u- + u+ = 1 + 3
I4_CF = 6.0  # unit-box surface


def _load(mode, w):
    path = os.path.join(RECORD_DIR, f"{mode}_n{w}.json")
    if not os.path.exists(path):
        return None
    with open(path) as f:
        return json.load(f)


def _rel(val, cf):
    return None if val is None else val / cf - 1.0


def _fmt(x):
    return "None" if x is None else f"{x:+.3e}"


def _table(mode, recs):
    print(f"\n[GEO-37 s2] ghost={mode}: width | cells | cells_ghost | count_201 | "
          "I1/0.48-1 | I2/2.88-1 | I3/1.92-1 | u+ only | I4/6-1 | dS(999) | "
          "c2 drop/0.48-1 | per-rank owned / one-sided 201 | s_mesh/import/assemble")
    for w, r in recs.items():
        if r is None:
            print(f"[GEO-37 s2] ghost={mode} {w} | MISSING RECORD (hang or crash)")
            continue
        print(
            f"[GEO-37 s2] ghost={mode} {w} | {r['cells']} | {r['cells_ghost']} | "
            f"{r['count_201']} | {_fmt(_rel(r['I1'], I1_CF))} | "
            f"{_fmt(_rel(r['I2'], I2_CF))} | {_fmt(_rel(r['I3'], I3_CF))} | "
            f"{r['I3_plus_only']!r} | {_fmt(_rel(r['I4'], I4_CF))} | "
            f"{r['c1_absent']!r} | {_fmt(_rel(r['c2_drop'], I1_CF))} | "
            f"{r['per_rank_owned_201']} / {r['per_rank_one_sided_201']} | "
            f"{r['seconds_mesh']:.2f}/{r['seconds_import']:.2f}/"
            f"{r['seconds_assemble']:.2f}"
        )
        errs = {k: v for k, v in r.items() if k.endswith("_error") and v}
        if errs:
            print(f"[GEO-37 s2] ghost={mode} {w} | errors: {errs}")


def _green(recs):
    def ok(r, w):
        if r is None or not r.get("completed") or r["comm_size"] != w:
            return False
        for k, cf in (("I1", I1_CF), ("I2", I2_CF), ("I3", I3_CF), ("I4", I4_CF)):
            if r[k] is None or abs(r[k] / cf - 1) > STRIP_BAND:
                return False
        return True
    return all(ok(r, w) for w, r in recs.items())


@pytest.fixture(scope="module")
def records():
    sf = {w: _load("sf", w) for w in WIDTHS}
    _table("shared_facet", sf)
    print(f"[GEO-37 step2] VERDICT dS over native tag (ghost shared_facet): "
          f"{'GREEN' if _green(sf) else 'RED'}")
    return sf


@pytest.mark.parametrize("w", WIDTHS)
def test_d1_no_hang(records, w):
    r = records[w]
    assert r is not None, f"sf_n{w}.json missing (width {w} hung or crashed)"
    assert r.get("completed") is True
    assert r["comm_size"] == w
    assert r["ghost"] == "shared_facet"


@pytest.mark.parametrize("w", WIDTHS)
def test_d2_area(records, w):
    r = records[w]
    assert r["I1_error"] is None, r["I1_error"]
    assert abs(r["I1"] / I1_CF - 1) <= STRIP_BAND, r["I1"]


@pytest.mark.parametrize("w", WIDTHS)
def test_d3_linear_and_two_sided(records, w):
    r = records[w]
    assert r["I2_error"] is None and r["I3_error"] is None
    assert abs(r["I2"] / I2_CF - 1) <= STRIP_BAND, r["I2"]
    assert abs(r["I3"] / I3_CF - 1) <= STRIP_BAND, r["I3"]


@pytest.mark.parametrize("w", WIDTHS)
def test_d4_exterior_instrument(records, w):
    r = records[w]
    assert r["I4_error"] is None, r["I4_error"]
    assert abs(r["I4"] / I4_CF - 1) <= STRIP_BAND, r["I4"]


@pytest.mark.parametrize("w", WIDTHS)
def test_d5_missing_tag_raises_via_validator(records, w):
    tags = records[w]["global_facet_tag_set"]
    validate_required_port_tags_exist(
        [PortDefinition("P1", 201, 101, "native_interior_sheet")], tags
    )
    with pytest.raises(ValueError, match="missing required port tags: \\[999\\]"):
        validate_required_port_tags_exist(
            [PortDefinition("P1", 999, 101, "absent_tag_control")], tags
        )


@pytest.mark.parametrize("w", WIDTHS)
def test_c1_absent_tag_assembles_to_zero(records, w):
    r = records[w]
    print(f"[GEO-37 s2] (c1) n={w} dS(999) = {r['c1_absent']!r}")
    assert r["c1_absent_error"] is None, r["c1_absent_error"]
    assert r["c1_absent"] == 0.0


@pytest.mark.parametrize("w", WIDTHS)
def test_c2_planted_drop_detected(records, w):
    r = records[w]
    assert r["c2_drop_error"] is None, r["c2_drop_error"]
    print(f"[GEO-37 s2] (c2) n={w} drop_rank={r['drop_rank']} "
          f"drop/0.48-1 = {r['c2_drop'] / I1_CF - 1:+.3e} "
          f"(one-facet share ~1/N_gmsh = {1.0 / r['N_gmsh']:.3e})")
    assert abs(r["c2_drop"] / I1_CF - 1) > STRIP_BAND


def test_zz_c3_unghosted_leg_printed(records):
    """(c3) ghost mode ``none`` (2026-08 configuration): printed, never asserted."""
    none = {w: _load("none", w) for w in WIDTHS}
    _table("none", none)
    print("[GEO-37 s2] (c3) ghost=none is a record, never a verdict")
