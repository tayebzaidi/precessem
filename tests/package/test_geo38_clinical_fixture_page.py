"""`GEO-38` C9 (§9 item 527, §10 chain C) — the clinical-fixture guide page is
the code.

Brief: docs/planning/chunks/GEO-38.md, "§9 item 527" and "Review ruling
2026-10-18 03:00". Ungated; no mesh is built.

Anchors (asserted):
  (P1) every compared cell of the table under ``## Geometry`` in
       docs/guide/clinical_coil_fixtures.md equals the kwargs constant
       formatted in the page's convention (2·ring_radius, 2·shield_radius,
       leg_spacing to 3 decimals in m; coil_length to 4; the generator's can
       length 2·z_extent to 3; legs and the census 2·legs; cell / unknown
       records with space-grouped thousands);
  (P2) the three source URLs of GEO-38.md "Nominal geometry" appear verbatim;
  (P3) the live mesh:15 guide's can sentence (corrected by §9 item 543)
       carries the generator's 2·z_extent on the 1.5 T kwargs at
       IDENTITY_BAND and formats to the page's "0.954 m" (no mismatch).
Negative controls (asserted):
  (N1) the same regex and comparator on the mesh:15 guide's sentence as it
       read before §9 item 543 (pinned literal MESH15_CAN_SENTENCE_PRE_543,
       "length 0.745 m") report exactly that one mismatch;
  (N2) a planted wrong digit in the page's 3 T can-length cell (0.700 ->
       0.770) reports exactly that one mismatch.
"""

from __future__ import annotations

import re
from pathlib import Path

import tests.validation.test_geo38_f_clinical_3t as c8
import tests.validation.test_geo38_f_clinical_build as s1a

ROOT = Path(__file__).resolve().parents[2]
PAGE = ROOT / "docs" / "guide" / "clinical_coil_fixtures.md"
CHUNK = ROOT / "docs" / "planning" / "chunks" / "GEO-38.md"
MESH15 = ROOT / "examples" / "meshing" / "15_f_clinical_coil_astm_box.md"

# the `mesh:15` guide's sentence before §9 item 543
# (`a373abc0:examples/meshing/15_f_clinical_coil_astm_box.md:11`)
MESH15_CAN_SENTENCE_PRE_543 = "PEC can (`shield_radius` 0.3725 m, length 0.745 m)"
MESH15_CAN_RE = r"PEC can \(`shield_radius` [0-9.]+ m, length ([0-9.]+) m\)"

SOURCE_URLS = (
    "https://mri-q.com/uploads/3/4/5/7/34572113/body_coil_distriubtion_mrm.20177.pdf",
    "https://pmc.ncbi.nlm.nih.gov/articles/PMC12137786/",
    "https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/11442125",
)


def _grouped(n: int) -> str:
    return f"{n:,}".replace(",", " ")


def can_length_m(kw: dict) -> float:
    """The generator's rule (`io/mesh.py:5636`): the can spans ±z_extent."""
    z_e = (
        max(
            0.5 * kw["coil_length"],
            0.5 * kw["leg_spacing"] + kw["ring_minor_radius"],
            0.5 * kw["phantom_height"],
        )
        + kw["air_padding"]
    )
    return 2.0 * z_e


def expected_cells(kw: dict, cells: tuple[int, int], unknowns_r1: tuple[int, int]) -> dict:
    legs = kw["leg_count"]
    return {
        "legs": str(legs),
        "ring-gap sheets (port census)": str(2 * legs),
        "rung (conductor) diameter": f"{2 * kw['ring_radius']:.3f} m",
        "RF shield diameter": f"{2 * kw['shield_radius']:.3f} m",
        "shield (can) length": f"{can_length_m(kw):.3f} m",
        "rung length (ring-centreline separation)": f"{kw['leg_spacing']:.3f} m",
        "leg span": f"{kw['coil_length']:.4f} m",
        "cells r1 / r2 (records)": f"{_grouped(cells[0])} / {_grouped(cells[1])}",
        "U₁ / U₂ at r1 (records)": f"{_grouped(unknowns_r1[0])} / {_grouped(unknowns_r1[1])}",
    }


EXPECTED = {
    "1.5": expected_cells(
        s1a.F_CLINICAL_15T_KWARGS,
        (s1a.F_CLINICAL_R1_CELL_RECORD, s1a.F_CLINICAL_R2_CELL_RECORD),
        s1a.F_CLINICAL_R1_UNKNOWNS_RECORD,
    ),
    "3": expected_cells(
        c8.F_CLINICAL_3T_KWARGS,
        (c8.F_CLINICAL_3T_R1_CELL_RECORD, c8.F_CLINICAL_3T_R2_CELL_RECORD),
        c8.F_CLINICAL_3T_R1_UNKNOWNS_RECORD,
    ),
}


def parse_geometry_table(text: str) -> dict[str, tuple[str, str]]:
    """{row label: (1.5 T cell, 3 T cell)} from the table under ``## Geometry``."""
    m = re.search(r"^## Geometry\s*$(.*?)(?=^## )", text, flags=re.M | re.S)
    assert m, "no '## Geometry' section"
    rows = {}
    for line in m.group(1).splitlines():
        line = line.strip()
        if not line.startswith("|") or set(line) <= set("|-: "):
            continue
        cells = [c.strip() for c in line.strip("|").split("|")]
        if len(cells) != 3 or cells[0] == "Quantity":
            continue
        rows[cells[0]] = (cells[1], cells[2])
    return rows


def mismatches(rows: dict[str, tuple[str, str]]) -> list[tuple[str, str, str, str]]:
    out = []
    for label in EXPECTED["1.5"]:
        got = rows.get(label)
        for col, key in ((0, "1.5"), (1, "3")):
            want = EXPECTED[key][label]
            have = None if got is None else got[col]
            if have != want:
                out.append((label, key, str(have), want))
    return out


def test_geo38_c9_p1_page_table_is_the_code():
    rows = parse_geometry_table(PAGE.read_text(encoding="utf-8"))
    for label in EXPECTED["1.5"]:
        print(
            f"[GEO-38 C9 P1] {label}: page {rows.get(label)} | code "
            f"({EXPECTED['1.5'][label]!r}, {EXPECTED['3'][label]!r})"
        )
    bad = mismatches(rows)
    print(f"[GEO-38 C9 P1] compared {2 * len(EXPECTED['1.5'])} cells, mismatches {len(bad)}")
    assert not bad, bad
    # Spot values from the brief, independent of the formatter.
    assert rows["shield (can) length"] == ("0.954 m", "0.700 m")
    assert rows["RF shield diameter"] == ("0.745 m", "0.743 m")
    assert rows["rung (conductor) diameter"] == ("0.710 m", "0.704 m")


def test_geo38_c9_p2_source_urls_verbatim():
    page = PAGE.read_text(encoding="utf-8")
    chunk = CHUNK.read_text(encoding="utf-8")
    for url in SOURCE_URLS:
        assert url in chunk, f"URL not in GEO-38.md: {url}"
        assert url in page, f"URL missing from page: {url}"
    print(f"[GEO-38 C9 P2] {len(SOURCE_URLS)} / {len(SOURCE_URLS)} source URLs verbatim")


def test_geo38_c9_n1_mesh15_can_sentence_is_caught():
    m = re.search(MESH15_CAN_RE, MESH15_CAN_SENTENCE_PRE_543)
    assert m, "pre-543 mesh:15 can sentence not matched"
    have = f"{m.group(1)} m"
    want = EXPECTED["1.5"]["shield (can) length"]
    rows = {"shield (can) length": (have, EXPECTED["3"]["shield (can) length"])}
    bad = [b for b in mismatches(rows) if b[0] == "shield (can) length"]
    print(f"[GEO-38 C9 N1] pre-543 mesh:15 can length {have!r} vs code {want!r}: {bad}")
    assert bad == [("shield (can) length", "1.5", "0.745 m", "0.954 m")]


def test_geo38_c9_p3_mesh15_can_sentence_is_the_rule():
    m = re.search(MESH15_CAN_RE, MESH15.read_text(encoding="utf-8"))
    assert m, "mesh:15 can sentence not found"
    v = float(m.group(1))
    rule = can_length_m(s1a.F_CLINICAL_15T_KWARGS)
    rel = abs(v / rule - 1.0)
    have = f"{v:.3f} m"
    rows = {"shield (can) length": (have, EXPECTED["3"]["shield (can) length"])}
    bad = [b for b in mismatches(rows) if b[0] == "shield (can) length"]
    print(
        f"[GEO-38 C9 P3] mesh:15 can length {v!r} vs rule {rule!r}: rel {rel:.3e} "
        f"(band {s1a.IDENTITY_BAND:.0e}); formatted {have!r}; mismatches {bad}"
    )
    assert rel <= s1a.IDENTITY_BAND, rel
    assert have == EXPECTED["1.5"]["shield (can) length"]
    assert bad == []


def test_geo38_c9_n2_planted_digit_is_the_only_mismatch():
    text = PAGE.read_text(encoding="utf-8")
    row = "| shield (can) length | 0.954 m | 0.700 m |"
    assert text.count(row) == 1
    planted = text.replace(row, "| shield (can) length | 0.954 m | 0.770 m |")
    bad = mismatches(parse_geometry_table(planted))
    print(f"[GEO-38 C9 N2] planted 0.700 -> 0.770: {bad}")
    assert bad == [("shield (can) length", "3", "0.770 m", "0.700 m")]
