"""`GEO-37` step 1 (§9 item 222): asserting reader for the interior-facet probe.

Reads the JSON records ``n1.json`` / ``n2.json`` / ``n4.json`` written by
``tests/mesh/probe_geo37_interior_facet.py`` into ``$GEO37_RECORD_DIR``.  The
module is skipped when that variable is unset, so ordinary suite runs never
fail on a missing record.

Anchors: (A1) no hang (record complete, right width); (A2) owned 201 count
equal across widths and equal to gmsh's triangle count; (A3) owned 201 area
= 0.48 at 1e-9; (A4) instrument check on the exterior group 101 (area 6 at
1e-9, count equal across widths); (c1) planted dropped facet is detected.
(c2), the unfiltered local-copy excess, is printed, never asserted.
"""

import json
import os

import pytest

RECORD_DIR = os.environ.get("GEO37_RECORD_DIR")
pytestmark = pytest.mark.skipif(
    RECORD_DIR is None, reason="GEO37_RECORD_DIR unset (probe records not produced)"
)

WIDTHS = (1, 2, 4)
RECT_AREA = 0.48
BOX_AREA = 6.0
BAND = 1e-9


def _load(w):
    path = os.path.join(RECORD_DIR, f"n{w}.json")
    if not os.path.exists(path):
        return None
    with open(path) as f:
        return json.load(f)


@pytest.fixture(scope="module")
def records():
    recs = {w: _load(w) for w in WIDTHS}
    print(
        "\n[GEO-37] width | cells | N_gmsh | count_201 | area_201-0.48 | "
        "local excess (c2) | per-rank owned 201 | count_101 | area_101-6 | "
        "s_mesh/s_import/s_measure"
    )
    for w, r in recs.items():
        if r is None:
            print(f"[GEO-37] {w} | MISSING RECORD")
            continue
        print(
            f"[GEO-37] {w} | {r['cells']} | {r['N_gmsh']} | {r['count_201']} | "
            f"{r['area_201'] - RECT_AREA:+.3e} | "
            f"{r['area_local_201'] - RECT_AREA:+.3e} "
            f"(count {r['count_local_201'] - r['count_201']:+d}) | "
            f"{r['per_rank_owned_201']} | {r['count_101']} | "
            f"{r['area_101'] - BOX_AREA:+.3e} | "
            f"{r['seconds_mesh']:.2f}/{r['seconds_import']:.2f}/"
            f"{r['seconds_measure']:.2f}"
        )
        print(
            f"[GEO-37] {w} | (c1) drop_rank={r['drop_rank']} count_drop="
            f"{r['count_drop']} area_drop/0.48-1={r['area_drop'] / RECT_AREA - 1:+.3e} "
            f"(one-facet share ~1/N_gmsh={1.0 / r['N_gmsh']:.3e}) | "
            f"cells_ghost={r['cells_ghost']}"
        )
    green = (
        all(
            r is not None
            and r.get("completed")
            and r["comm_size"] == w
            and r["count_201"] == r["N_gmsh"]
            and abs(r["area_201"] / RECT_AREA - 1) <= BAND
            for w, r in recs.items()
        )
        and len({recs[w]["count_201"] for w in WIDTHS if recs[w]}) == 1
    )
    print(f"[GEO-37] VERDICT native interior facet tags on 0.11: " f"{'GREEN' if green else 'RED'}")
    return recs


@pytest.mark.parametrize("w", WIDTHS)
def test_a1_no_hang(records, w):
    r = records[w]
    assert r is not None, f"n{w}.json missing (width {w} hung or crashed)"
    assert r.get("completed") is True
    assert r["comm_size"] == w


def test_a2_count_equal_across_widths_and_to_gmsh(records):
    counts = {w: records[w]["count_201"] for w in WIDTHS}
    n_gmsh = {w: records[w]["N_gmsh"] for w in WIDTHS}
    assert len(set(n_gmsh.values())) == 1, n_gmsh
    assert len(set(counts.values())) == 1, counts
    for w in WIDTHS:
        assert counts[w] == n_gmsh[w], (w, counts[w], n_gmsh[w])


@pytest.mark.parametrize("w", WIDTHS)
def test_a3_area_closed_form(records, w):
    r = records[w]
    assert abs(r["area_201"] / RECT_AREA - 1) <= BAND, r["area_201"]


@pytest.mark.parametrize("w", WIDTHS)
def test_a4_exterior_instrument(records, w):
    r = records[w]
    assert abs(r["area_101"] / BOX_AREA - 1) <= BAND, r["area_101"]
    assert r["count_101"] == records[1]["count_101"]


@pytest.mark.parametrize("w", WIDTHS)
def test_c1_planted_drop_detected(records, w):
    r = records[w]
    assert r["count_drop"] == r["N_gmsh"] - 1
    assert abs(r["area_drop"] / RECT_AREA - 1) > BAND
