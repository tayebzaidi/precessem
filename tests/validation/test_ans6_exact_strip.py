"""`ANS-6` step 4 (§9 item 193, chain F16) — exact SPEC port rectangles on the
hole / copper route, and the copper four-port re-recorded at 10 / 64 / 128 MHz.

The `ANS-6` / `ANS-7` fixture is `TH-14`'s hole mesh (``_sheets_build(True,
as_hole=True)``, 80 181 cells).  Its ports were a ragged midpoint-filtered
sub-sheet of each whole F_i face.  ``leg_sheet_strip_fraction`` (`OPS-61` step
1b) now composes with ``as_hole`` (the refusal was conservative, not
structural), so the strip edges are CAD edges on the hole route as well and
``_build_ports``' unchanged midpoint filter selects exactly the 0.007 x 0.008
SPEC rectangle.  Operator ruling 2026-10-02: exact CAD rectangles, κ off;
``measure_kappa`` only reads the terminal-form deficit as a diagnostic.

Opt-in through ``FEM_EM_EXACT_ANS6`` ∈ {``mesh``, ``solve``}; unset ⇒ every
test skips, any other value raises.  (Not ``FEM_EM_ANS6_*``: ``ans:7``'s
``main()`` refuses every such variable.)

Knob notes.  ``FEM_EM_EXACT_ANS6_DEGREE`` ∈ {1, 2} (§9 item 216),
``FEM_EM_EXACT_ANS6_PART`` ∈ {sweeps, identity} (§9 item 226), and
``FEM_EM_EXACT_ANS6_RESOLUTION`` (§9 item 266, `ANS-6` step 5b′, xl entry 31's
contract): a finite float > 0 forwarded verbatim as ``resolution`` to both
``_hole_rung`` calls in ``meshes`` and to ``_build_ladder`` in ``ladder``;
unset or blank forwards nothing (the F19 records' call, digit-identical).
Knobbed: (M6)/(M7) become (F′) (cells > ``KNOB_CELL_FACTOR`` × the ×1 count),
(R) and (N3) skip, the h-move against the ×1 record is printed only, and the
JSON name gains ``_r<res>``.  Parsed in fixtures only, never at import.

Anchors (asserted, bands imported, never restated): (M1)–(M7) on the mesh,
(G) / (b) / (σ) / (c) / (R) on the solve, negative controls (N1) printed and
(N2) floor asserted.  Areas, facet counts and extents are rank-local and
reduced before asserting.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pytest
from mpi4py import MPI

from fem_em_solver.utils.instrumentation import report_peak_rss
from tests.complex_mode import complex_only
from tests.mesh.test_birdcage_conductor_hole import EXACT
from tests.mesh.test_birdcage_leg_offset import (
    SHEET_AREA_BAND,
    _projected_extents,
    _sheet_areas,
    _sheet_azimuth_deg,
)
from tests.mesh.test_birdcage_port_sheet_prerequisite import CELL_COUNT_BAND
from tests.mesh.test_birdcage_port_sheets import (
    PORT_LOWER,
    PORT_UPPER,
    SHEET_IFACE,
    _build as _sheets_build,
)
from tests.mesh.test_birdcage_port_tags import LEG_COUNT, RING_RADIUS
from tests.mesh.test_birdcage_port_terminals import _interface_area_or_zero
from tests.mesh.test_coil_phantom_conforming import _tag_volume
from tests.validation import test_th14_birdcage_copper as th14
from tests.validation.test_ans4_exact_strip_records import (
    ANS4_EXACT_STRIP_S_P1_RECORDS,
    KNOB_CELL_FACTOR,
    SPEC_AREA_M2,
    SPEC_HEIGHT_M,
    SPEC_WIDTH_M,
)
from tests.validation.test_ans4_resolution_ladder import STEP3A_RECORD_RTOL
from tests.validation.test_birdcage_power_identity import DISCRETE_IDENTITY_RTOL
from tests.validation.test_ops61_exact_strip_ports import (
    NEGATIVE_CONTROL_FLOOR,
    STRIP_BAND,
)
from tests.validation.test_port21_feed_sensitivity import (
    KNOB_REACHES_THE_SOLVE,
    _class_means,
    _relative_move,
)
from tests.validation.test_port_birdcage_four_port import PASSIVITY_SIGMA_TOLERANCE
from tests.validation.test_port_birdcage_larmor_gate import LEG_D_S_MATRIX_10MHZ
from tests.validation.test_port_birdcage_leg_offset_sweep import (
    _narrowed_radial,
    _port_frame,
)
from tests.validation.test_port_birdcage_lumped_column import ADJACENT_SPREAD_BAND
from tests.validation.test_port_lumped_narrowed_sheet import GATED_WIDTH_FRACTION
from tests.validation.test_port_lumped_sheet_sweep import RECIPROCITY_BAND

REPO = Path(__file__).resolve().parents[2]
OUT_DIR = REPO / "logs" / "ans6_exact"
MODE_ENV = "FEM_EM_EXACT_ANS6"
ANS6_METRICS = (
    REPO / "examples" / "ansys_benchmarks"
    / "ans6_copper_birdcage_four_port_10_64_128MHz" / "metrics.json"
)
ANS7_METRICS = (
    REPO / "examples" / "ansys_benchmarks"
    / "ans7_tuned_copper_birdcage_driven_sweep" / "metrics.json"
)
PORTS = list(range(1, LEG_COUNT + 1))
CLASSES = ("self", "adjacent", "opposite")
FREQ_KEYS = ("10", "64", "128")
# The ragged P_coil/P_in of `ans:6`'s runnable half (printed comparand only).
RAGGED_P_COIL_OVER_P_IN = {"10": 0.9288392, "64": 0.4483664, "128": 0.2186740}

# (M7) — flag-on hole cell count, version-tagged: 80 798 measured at -n 1 on
# the complex 0.11 image `main` boots (`20261003T095457Z_ANS-6-item193-W1-mesh.log`
# :1740; ratio 1.007695 to the 80 181 flag-off mesh).  Mesher-dependent —
# re-record under the (1*) licence on an image change, never re-band.
ANS6_EXACT_STRIP_HOLE_CELL_RECORD = 80798

# (R) — degree-1 κ-off class means on the exact-rectangle hole mesh (80 798
# cells), per frequency per column, 10 digits.  Version tag: complex 0.11
# image, `-n 2`, `20261003T095650Z_ANS-6-item193-W2-solve.log:3091–3101`;
# (1*) licence — re-record on an image change, never re-band.
ANS6_EXACT_STRIP_CLASS_MEANS_P1 = {
    "10": {
        "cu": {"self": -4.835946573e-01 + 1.496865370e-01j,
               "adjacent": +4.968908987e-01 - 4.219890799e-02j,
               "opposite": +4.897195031e-01 - 7.859240439e-02j},
        "pec": {"self": -4.837960870e-01 + 1.495553853e-01j,
                "adjacent": +4.969410045e-01 - 4.216223166e-02j,
                "opposite": +4.898207240e-01 - 7.853460271e-02j},
    },
    "64": {
        "cu": {"self": -3.894212853e-02 + 6.590775812e-01j,
               "adjacent": +3.977072490e-01 - 2.222710707e-01j,
               "opposite": +2.397020121e-01 - 3.002075581e-01j},
        "pec": {"self": -3.933329848e-02 + 6.592411299e-01j,
                "adjacent": +3.978340023e-01 - 2.222897411e-01j,
                "opposite": +2.398399832e-01 - 3.003335088e-01j},
    },
    "128": {
        "cu": {"self": +4.654952469e-01 + 6.579036879e-01j,
               "adjacent": +2.385171187e-01 - 2.935241150e-01j,
               "opposite": +4.167270495e-02 - 2.459835329e-01j},
        "pec": {"self": +4.653444453e-01 + 6.581856746e-01j,
                "adjacent": +2.385966529e-01 - 2.936075081e-01j,
                "opposite": +4.166643522e-02 - 2.460973517e-01j},
    },
}
ANS6_EXACT_STRIP_RECORD_RANK_WIDTH = 2

# (R) at degree 2 (§9 item 216, `ANS-6` step 5) — κ-off class means on the same
# 80 798-cell exact-rectangle hole mesh, per frequency per column, 10 digits.
# Version tag: complex 0.11 image, `-n 8`, registered by §9 item 226 and item 232
# (128 MHz) from the sweeps windows: 10 MHz
# `20261004T184300Z_ANS-6-item226-W1a-p2-10-sweeps.log:1256–1257`,
# 64 MHz `20261004T184757Z_ANS-6-item226-W2a-p2-64-sweeps.log:1256–1257`, 128 MHz
# `20261005T004117Z_ANS-6-item232-W3a-p2-128-sweeps.log:1256–1257`.  (1*) licence —
# re-record on an image change, never re-band.
ANS6_EXACT_STRIP_CLASS_MEANS_P2 = {
    "10": {
        "cu": {"self": -4.821901317e-01 + 1.563310044e-01j,
               "adjacent": +4.966236707e-01 - 4.361263927e-02j,
               "opposite": +4.888645793e-01 - 8.135292126e-02j},
        "pec": {"self": -4.823981130e-01 + 1.561975831e-01j,
                "adjacent": +4.966752944e-01 - 4.357522874e-02j,
                "opposite": +4.889693156e-01 - 8.129431878e-02j},
    },
    "64": {
        "cu": {"self": -1.276864349e-02 + 6.709270564e-01j,
               "adjacent": +3.908392789e-01 - 2.255913789e-01j,
               "opposite": +2.278736968e-01 - 2.986429559e-01j},
        "pec": {"self": -1.315445856e-02 + 6.711050243e-01j,
                "adjacent": +3.909665706e-01 - 2.256145287e-01j,
                "opposite": +2.280051972e-01 - 2.987743957e-01j},
    },
    "128": {
        "cu": {"self": +4.918927876e-01 + 6.533106200e-01j,
               "adjacent": +2.288719133e-01 - 2.895746771e-01j,
               "opposite": +3.702806662e-02 - 2.355310451e-01j},
        "pec": {"self": +4.917553754e-01 + 6.535906126e-01j,
                "adjacent": +2.289466054e-01 - 2.896597086e-01j,
                "opposite": +3.701783028e-02 - 2.356397332e-01j},
    },
}
ANS6_EXACT_STRIP_P2_RECORD_RANK_WIDTH = 8

DEGREE_ENV = "FEM_EM_EXACT_ANS6_DEGREE"
# §9 item 226 (`ANS-6` step 5): split each frequency into a ``sweeps`` window
# (PEC + the whole σ ladder, no field solve) and an ``identity`` window (PEC +
# copper + P1).  Not ``FEM_EM_ANS6_*``: ``ans:7``'s ``main()`` refuses that prefix.
PART_ENV = "FEM_EM_EXACT_ANS6_PART"
PARTS = ("sweeps", "identity")
# §9 item 266 (`ANS-6` step 5b′, xl entry 31's contract): mesh resolution knob.
RESOLUTION_ENV = "FEM_EM_EXACT_ANS6_RESOLUTION"
# Predicted (printed, never asserted) h = 0.0075 cell brackets: flag-off 226 642
# exactly (`20260923T093956Z_ANS-6-step2b-b-h0p0075-d1-64MHz.log:950`), flag-on
# in [226 642, 232 000] (×1 strip ratio 1.007695 ⇒ ≈ 228 400).
PREDICTED_R0P0075_OFF_CELLS = 226642
PREDICTED_R0P0075_ON_BRACKET = (226642, 232000)


def _resolution():
    """``FEM_EM_EXACT_ANS6_RESOLUTION``: unset or blank ⇒ ``None`` (no
    ``resolution`` forwarded — the F19 records' call); otherwise a finite
    float > 0, anything else raises.  Parsed in fixtures only, never at import
    (`test_ans7_exact_strip.py` imports this module)."""
    raw = os.environ.get(RESOLUTION_ENV, "").strip()
    if not raw:
        return None
    try:
        value = float(raw)
    except ValueError:
        raise ValueError(f"{RESOLUTION_ENV}={raw!r}: not a float") from None
    if not np.isfinite(value) or value <= 0.0:
        raise ValueError(f"{RESOLUTION_ENV}={raw!r}: must be finite and > 0")
    return value


def _resolution_set():
    """Non-raising read for the skip markers."""
    return bool(os.environ.get(RESOLUTION_ENV, "").strip())


def _part():
    """``FEM_EM_EXACT_ANS6_PART`` ∈ {unset, ``sweeps``, ``identity``}; unset ⇒
    ``None`` (the item-193 / item-216 call).  Parsed in fixtures only, never at
    import (`test_ans7_exact_strip.py` imports this module)."""
    raw = os.environ.get(PART_ENV, "").strip()
    if not raw:
        return None
    if raw not in PARTS:
        raise ValueError(f"{PART_ENV}={raw!r}: only 'sweeps' or 'identity'")
    return raw


def _part_set():
    """Non-raising read for the skip markers."""
    return bool(os.environ.get(PART_ENV, "").strip())


def _degree():
    """``FEM_EM_EXACT_ANS6_DEGREE`` ∈ {1, 2}; unset ⇒ ``None`` (no ``degree``
    forwarded, ``_build_ladder``'s default 1).  Parsed in fixtures only, never
    at import (`test_ans7_exact_strip.py` imports this module)."""
    raw = os.environ.get(DEGREE_ENV, "").strip()
    if not raw:
        return None
    if raw not in ("1", "2"):
        raise ValueError(f"{DEGREE_ENV}={raw!r}: only '1' or '2'")
    return int(raw)


def _mode():
    raw = os.environ.get(MODE_ENV, "").strip()
    if not raw:
        return None
    if raw not in ("mesh", "solve"):
        raise ValueError(f"{MODE_ENV}={raw!r}: only 'mesh' or 'solve'")
    return raw


needs_mode = pytest.mark.skipif(_mode() is None, reason=f"{MODE_ENV} unset (opt-in window)")
needs_solve = pytest.mark.skipif(_mode() != "solve", reason=f"{MODE_ENV} != solve")
skip_if_part = pytest.mark.skipif(_part_set(), reason="asserted in window P and W0")


def _owned_facet_count(msh, facet_tags, tag, comm):
    fdim = msh.topology.dim - 1
    n_owned = msh.topology.index_map(fdim).size_local
    facets = np.asarray(facet_tags.find(tag), dtype=np.int64)
    return int(comm.allreduce(int(np.count_nonzero(facets < n_owned)), op=MPI.SUM))


def _hole_rung(strip_fraction, resolution=None):
    """`test_ops61_exact_strip_ports._rung`'s measurement on the hole mesh.

    ``resolution`` (§9 item 214, `ANS-7` step 6, rule (a), additive): ``None``
    -- the default -- forwards nothing, so the call is byte-identical to
    before; a float is passed to ``_sheets_build``'s ``resolution`` keyword.
    """
    comm = MPI.COMM_WORLD
    extra = {} if strip_fraction is None else {"leg_sheet_strip_fraction": strip_fraction}
    if resolution is not None:
        extra["resolution"] = float(resolution)
    msh, cell_tags, _facets, diag, elapsed = _sheets_build(True, as_hole=True, **extra)
    tdim = msh.topology.dim
    msh.topology.create_connectivity(tdim - 1, tdim)
    tags_full, full_areas = _sheet_areas(msh, cell_tags, PORTS, comm)
    dx, _dy, dz = (float(v) for v in diag["port_box_size_m"])
    box_volume = float(np.prod(diag["port_box_size_m"]))
    ports = {}
    for i in PORTS:
        tag = SHEET_IFACE + i
        azimuth = _sheet_azimuth_deg(msh, tags_full, tag, comm)
        frame = _port_frame(azimuth)
        spans, centres = _projected_extents(msh, tags_full, tag, comm, list(frame))
        ports[i] = {"tag": tag, "frame": frame, "w_full": float(spans[0]),
                    "centre_radial": float(centres[0]), "F": float(full_areas[i]),
                    "F_facets": _owned_facet_count(msh, tags_full, tag, comm)}
    tags_n = tags_full
    for i in PORTS:
        p = ports[i]
        tags_n = _narrowed_radial(msh, tags_n, p["tag"], GATED_WIDTH_FRACTION,
                                  p["centre_radial"], p["frame"][0], 0.5 * p["w_full"])
    for i in PORTS:
        p = ports[i]
        full = np.asarray(tags_full.find(p["tag"]), dtype=np.int64)
        narrow = np.asarray(tags_n.find(p["tag"]), dtype=np.int64)
        p["subset_local"] = bool(np.all(np.isin(narrow, full)))
        p["N"] = float(_interface_area_or_zero(msh, tags_n, p["tag"], comm))
        p["N_facets"] = _owned_facet_count(msh, tags_n, p["tag"], comm)
        spans_n, centres_n = _projected_extents(msh, tags_n, p["tag"], comm, list(p["frame"]))
        p["N_span_radial"] = float(spans_n[0])
        p["N_h"] = float(spans_n[2])
        p["N_centre_radial"] = float(centres_n[0])
        p["halves"] = (_tag_volume(msh, cell_tags, PORT_LOWER + i, comm),
                       _tag_volume(msh, cell_tags, PORT_UPPER + i, comm))
    subset_all = bool(comm.allreduce(all(ports[i]["subset_local"] for i in PORTS), op=MPI.LAND))
    return {"cells": int(msh.topology.index_map(tdim).size_global), "dx": dx, "dz": dz,
            "box_volume": box_volume, "ports": ports, "subset_all": subset_all,
            "cavity_cad": diag["conductor_cavity_cad"], "elapsed": elapsed}


@pytest.fixture(scope="module")
def meshes():
    comm = MPI.COMM_WORLD
    res = _resolution()  # §9 item 266: a bad value raises on every rank before any mesh
    _degree()  # §9 item 216: a bad value raises on every rank before any mesh
    built = {"on": _hole_rung(GATED_WIDTH_FRACTION, resolution=res),
             "off": _hole_rung(None, resolution=res)}
    if comm.rank == 0:
        for name, r in built.items():
            print(f"\n[ANS-6 step4] hole flag {name}: cells={r['cells']} build "
                  f"{r['elapsed']:.2f} s dx*dz={r['dx'] * r['dz']:.9e} cavity CAD "
                  f"{r['cavity_cad']}", flush=True)
            for i in PORTS:
                p = r["ports"][i]
                comp = p["F"] - p["N"]
                print(f"    P{i}: F={p['F']:.9e} ({p['F_facets']} tri) N={p['N']:.9e} "
                      f"({p['N_facets']} tri, rel {p['N'] / SPEC_AREA_M2 - 1.0:+.3e}) "
                      f"w=A/h={p['N'] / p['N_h']:.9e} radial span={p['N_span_radial']:.9e} "
                      f"h={p['N_h']:.9e} centre r={p['N_centre_radial']:.9e} "
                      f"F\\N={comp:.9e} halves sum={sum(p['halves']) / r['box_volume']:.12f}",
                      flush=True)
        on = built["on"]["cells"]
        if res is None:
            print(f"[ANS-6 step4] strip hole cells {on} (PREDICTED [80600, 81800]: "
                  f"{'met' if 80600 <= on <= 81800 else 'NOT met'}; ratio to flag-off "
                  f"{on / built['off']['cells']:.6f})", flush=True)
        else:
            lo, hi = PREDICTED_R0P0075_ON_BRACKET
            print(f"[ANS-6 step5b′] resolution {res:g}: strip hole cells {on} (printed; "
                  f"PREDICTED at 0.0075 [{lo}, {hi}]: {'met' if lo <= on <= hi else 'NOT met'};"
                  f" ratio to flag-off {on / built['off']['cells']:.6f}); flag-off "
                  f"{built['off']['cells']} (printed; PREDICTED at 0.0075 "
                  f"{PREDICTED_R0P0075_OFF_CELLS} exactly)", flush=True)
    return built


@complex_only
@needs_mode
@skip_if_part
def test_m1_full_sheet_is_the_box_face(meshes):
    r = meshes["on"]
    for i in PORTS:
        ratio = r["ports"][i]["F"] / (r["dx"] * r["dz"])
        assert abs(ratio - 1.0) < SHEET_AREA_BAND, f"P{i}: F/dxdz={ratio:.12f}"


@complex_only
@needs_mode
@skip_if_part
def test_m2_strip_is_the_spec_rectangle(meshes):
    for i in PORTS:
        p = meshes["on"]["ports"][i]
        checks = {
            "area": (p["N"], SPEC_AREA_M2),
            "w=A/h": (p["N"] / p["N_h"], SPEC_WIDTH_M),
            "radial span": (p["N_span_radial"], SPEC_WIDTH_M),
            "h": (p["N_h"], SPEC_HEIGHT_M),
            "centre radius": (p["N_centre_radial"], RING_RADIUS),
        }
        for name, (got, want) in checks.items():
            assert abs(got / want - 1.0) < STRIP_BAND, (
                f"P{i} {name}: {got:.12e} vs {want:.12e} ({p['N_facets']} tri)")


@complex_only
@needs_mode
@skip_if_part
def test_m3_complement_and_subset(meshes):
    r = meshes["on"]
    assert r["subset_all"], "some N_i is not a subset of F_i"
    for i in PORTS:
        comp = r["ports"][i]["F"] - r["ports"][i]["N"]
        assert abs(comp / SPEC_AREA_M2 - 1.0) < STRIP_BAND, f"P{i}: F\\N={comp:.12e}"


@complex_only
@needs_mode
@skip_if_part
def test_m4_box_halves_partition_the_box(meshes):
    r = meshes["on"]
    for i in PORTS:
        partition = sum(r["ports"][i]["halves"]) / r["box_volume"]
        assert abs(partition - 1.0) < EXACT, f"P{i}: halves sum {partition:.12f}"


@complex_only
@needs_mode
@skip_if_part
def test_m5_strip_leaves_the_cavity_alone(meshes):
    on = meshes["on"]["cavity_cad"]
    off = meshes["off"]["cavity_cad"]
    assert on["surfaces"] == off["surfaces"], (on, off)
    assert abs(on["area_m2"] / off["area_m2"] - 1.0) < EXACT, (on, off)


@complex_only
@needs_mode
@skip_if_part
def test_m6_default_hole_cells_unchanged(meshes):
    tracked = int(json.loads(ANS7_METRICS.read_text())["n_cells"])
    cells = meshes["off"]["cells"]
    if _resolution() is not None:
        # (F′) §9 item 266: the knob reaches the flag-off hole mesh.
        if MPI.COMM_WORLD.rank == 0:
            print(f"\n[ANS-6 step5b′] (F′) M6 flag-off hole cells = {cells} > "
                  f"{KNOB_CELL_FACTOR} x {tracked} (asserted); PREDICTED at 0.0075 "
                  f"{PREDICTED_R0P0075_OFF_CELLS} (printed)", flush=True)
        assert cells > KNOB_CELL_FACTOR * tracked, f"flag off {cells} vs {tracked}"
        return
    assert abs(cells / tracked - 1.0) < CELL_COUNT_BAND, f"flag off {cells} vs {tracked}"


@complex_only
@needs_mode
@skip_if_part
def test_m7_strip_hole_cell_record(meshes):
    cells = meshes["on"]["cells"]
    if MPI.COMM_WORLD.rank == 0:
        print(f"\n[ANS-6 step4] M7 flag-on hole cells = {cells} "
              f"(record {ANS6_EXACT_STRIP_HOLE_CELL_RECORD})", flush=True)
    if ANS6_EXACT_STRIP_HOLE_CELL_RECORD is None:
        pytest.skip("M7 record not yet registered (W1 measures it)")
    if _resolution() is not None:
        # (F′) §9 item 266: the knob reaches the flag-on hole mesh.
        lo, hi = PREDICTED_R0P0075_ON_BRACKET
        if MPI.COMM_WORLD.rank == 0:
            print(f"[ANS-6 step5b′] (F′) M7 flag-on {cells} > {KNOB_CELL_FACTOR} x "
                  f"{ANS6_EXACT_STRIP_HOLE_CELL_RECORD} (asserted); PREDICTED at 0.0075 "
                  f"[{lo}, {hi}] (printed)", flush=True)
        assert cells > KNOB_CELL_FACTOR * ANS6_EXACT_STRIP_HOLE_CELL_RECORD
        return
    assert abs(cells / ANS6_EXACT_STRIP_HOLE_CELL_RECORD - 1.0) < CELL_COUNT_BAND


@complex_only
@needs_mode
@skip_if_part
def test_n1_flag_off_hole_strip_is_ragged(meshes):
    """(N1) predicted, printed: the flag-off hole N_i misses the SPEC area."""
    for i in PORTS:
        miss = meshes["off"]["ports"][i]["N"] / SPEC_AREA_M2 - 1.0
        if MPI.COMM_WORLD.rank == 0:
            print(f"[ANS-6 step4] N1 (PREDICTED ≈ +4e-2, printed) P{i}: miss {miss:+.3e} "
                  f"({'exceeds' if abs(miss) > NEGATIVE_CONTROL_FLOOR else 'DOES NOT exceed'}"
                  f" floor {NEGATIVE_CONTROL_FLOOR:.0e})", flush=True)


def _cmeans_json(means):
    return {c: [float(z.real), float(z.imag)] for c, z in means.items()}


def _ragged_cu(key):
    m = json.loads(ANS6_METRICS.read_text())
    s = m["rungs"][f"{key} MHz"]["cu"]["s_matrix"]
    return np.array([[complex(e["re"], e["im"]) for e in row] for row in s])


@pytest.fixture(scope="module")
def ladder(request):
    comm = MPI.COMM_WORLD
    res = _resolution()  # §9 item 266: a bad value raises on every rank before any mesh
    degree = _degree()
    part = _part()  # §9 item 226: a bad value raises on every rank before any mesh
    if part is None:
        # Unset: the branch's ``ladder(meshes)`` (the M-tests built it already).
        request.getfixturevalue("meshes")
    # §9 item 216: unset forwards no ``degree`` (identical call to item 193's).
    extra = {} if degree is None else {"degree": degree}
    if part == "sweeps":
        extra["field_solve"] = False
    elif part == "identity":
        extra["sigmas"] = (th14.COPPER_SIGMA,)
    if res is not None:
        extra["resolution"] = res  # §9 item 266: unset adds no key (item 232's call)
    if comm.rank == 0:
        res_txt = "unset" if res is None else f"{res:g}"
        print(f"\n[ANS-6 step5b′] knobs: resolution {res_txt}, "
              f"degree {'unset' if degree is None else degree}, part "
              f"{'unset' if part is None else part}, frequencies "
              f"{os.environ.get(th14.FREQS_ENV, 'all')}", flush=True)
    lad = th14._build_ladder(leg_sheet_strip_fraction=GATED_WIDTH_FRACTION, measure_kappa=True,
                             **extra)
    if res is not None:
        report_peak_rss(comm, label="ANS-6 step5b′ after ladder")
    lad["degree"] = 1 if degree is None else degree
    lad["part"] = part
    lad["resolution"] = res
    means = {k: {"cu": _class_means(lad["freqs"][k]["sigma"][th14.COPPER_SIGMA]["s"]),
                 "pec": _class_means(lad["freqs"][k]["pec"]["s"])} for k in lad["freqs"]}
    moves = {k: {"cu": _relative_move(means[k]["cu"], _class_means(_ragged_cu(k)))}
             for k in lad["freqs"]}
    # (N3), degree 2 only: the move of every class mean from the degree-1 record.
    order_moves = {}
    # §9 item 266: knobbed, (N3) would mix h and order — suppressed.
    if lad["degree"] == 2 and res is None:
        order_moves = {k: {col: _relative_move(means[k][col],
                                               ANS6_EXACT_STRIP_CLASS_MEANS_P1[k][col])
                           for col in ("cu", "pec")} for k in lad["freqs"]}
    lad["means"] = means
    lad["moves"] = moves
    lad["order_moves"] = order_moves
    if comm.rank == 0:
        print(f"\n[ANS-6 step4] solve: degree {lad['degree']}, part {part}, cells "
              f"{lad['cells']}, unknowns {lad['unknowns']}, -n {comm.size}, total "
              f"{lad['t_total']:.1f} s", flush=True)
        for k, rec in lad["freqs"].items():
            # §9 item 226 rank-0 guards: power / kappa exist only with the field
            # solve, the PEC-limit rung only with the whole σ ladder.
            txt = []
            if "kappa" in rec:
                txt.append(f"kappa pooled {rec['kappa']['pooled']:.6e} (PREDICTED [5e-6, 5e-5])")
            if "power" in rec:
                p = rec["power"]
                p_in = p["p_src"] - p["p_sheets"]
                ratio = p["p_surf"] / p_in
                rag = RAGGED_P_COIL_OVER_P_IN[k]
                txt.append(f"identity residual {rec['identity_residual']:.3e} (PREDICTED <~ "
                           f"1e-11); P_coil/P_in {ratio:.7f} vs ragged {rag:.7f} "
                           f"({ratio / rag - 1.0:+.3%}, PREDICTED within ±2 %)")
            last = rec["sigma"].get(th14.SIGMA_LADDER[-1])
            if last is not None:
                txt.append(f"PEC-limit max|dS| {last['max_abs_ds']:.3e} (PREDICTED <= "
                           f"{th14.PREDICTED_PEC_LIMIT_MAX_ABS_DS:.0e})")
            print(f"[ANS-6 step4] {k} MHz: " + "; ".join(txt), flush=True)
            for col in ("cu", "pec"):
                print(f"    {col} class means (10 digits): " + ", ".join(
                    f"{c} {means[k][col][c].real:+.9e}{means[k][col][c].imag:+.9e}j"
                    for c in CLASSES), flush=True)
            print(f"    cu exact vs ragged ans:6 metrics.json class move: " + "  ".join(
                f"{c} {moves[k]['cu'][c] * 100:.4f}%" for c in CLASSES)
                + (" (N2, floor asserted; PREDICTED self 0.5-2 %)" if k == "64" else
                   " (printed)")
                + (" [degree 2 vs degree-1 ragged: measures order plus sheet]"
                   if lad["degree"] == 2 else ""), flush=True)
            if order_moves:
                for col in ("cu", "pec"):
                    print(f"    (N3) {k} MHz {col} degree-2 vs degree-1 exact class move "
                          f"(floor {KNOB_REACHES_THE_SOLVE:.0e} asserted; PREDICTED 64 MHz "
                          f"self 5-9 % couplings 2-6 %, 10 MHz self 1-4 %, 128 MHz self "
                          f"4-9 %): " + "  ".join(
                              f"{c} {order_moves[k][col][c] * 100:.4f}%" for c in CLASSES),
                          flush=True)
            cu_m, pec_m = means[k]["cu"], means[k]["pec"]
            print(f"    cu vs pec class-mean relative gap |cu-pec|/|pec| (printed; degree 1 "
                  f"read <= 0.0642 %, item 216 W0): " + "  ".join(
                f"{c} {abs(cu_m[c] - pec_m[c]) / abs(pec_m[c]) * 100:.4f}%" for c in CLASSES),
                flush=True)
            if res is not None:
                # §9 item 266: entry 31's readout — h-move against the ×1 record.
                x1 = (ANS6_EXACT_STRIP_CLASS_MEANS_P2 if lad["degree"] == 2
                      else ANS6_EXACT_STRIP_CLASS_MEANS_P1)
                if k in x1:
                    for col in ("cu", "pec"):
                        hmv = _relative_move(means[k][col], x1[k][col])
                        print(f"    h-move {k} MHz {col} r={res:g} vs x1 record degree "
                              f"{lad['degree']} (printed, never asserted): " + "  ".join(
                                  f"{c} {hmv[c] * 100:.4f}%" for c in CLASSES), flush=True)
                    x_cu, x_pec = x1[k]["cu"], x1[k]["pec"]
                    print(f"    Cu-PEC gap |cu-pec|/|pec| rung r={res:g} vs x1 record "
                          f"(printed, never asserted): " + "  ".join(
                              f"{c} {abs(cu_m[c] - pec_m[c]) / abs(pec_m[c]) * 100:.4f}% vs "
                              f"{abs(x_cu[c] - x_pec[c]) / abs(x_pec[c]) * 100:.4f}%"
                              for c in CLASSES), flush=True)
        OUT_DIR.mkdir(parents=True, exist_ok=True)
        if degree is None and th14.FREQS_ENV not in os.environ:
            out = OUT_DIR / "ans6_exact_p1.json"
        else:
            out = OUT_DIR / f"ans6_exact_p{lad['degree']}_f{'-'.join(lad['freqs'])}.json"
        if part is not None:
            out = OUT_DIR / f"ans6_exact_p{lad['degree']}_f{'-'.join(lad['freqs'])}_{part}.json"
        if res is not None:
            tail = f"_{part}" if part is not None else ""
            out = OUT_DIR / (f"ans6_exact_p{lad['degree']}_r{res:g}_f"
                             f"{'-'.join(lad['freqs'])}{tail}.json")

        def _cj(s):
            return [[[float(z.real), float(z.imag)] for z in row] for row in np.asarray(s)]

        payload_extra = {} if degree is None else {
            "degree": lad["degree"], "order_moves_vs_p1": order_moves}
        out.write_text(json.dumps({
            **payload_extra,
            "ranks": comm.size, "cells": lad["cells"], "unknowns": lad["unknowns"],
            "census": lad["census"],
            "s": {k: {"cu": _cj(r["sigma"][th14.COPPER_SIGMA]["s"]), "pec": _cj(r["pec"]["s"])}
                  for k, r in lad["freqs"].items()},
            "class_means": {k: {col: _cmeans_json(v) for col, v in m.items()}
                            for k, m in means.items()},
            "kappa": {k: r.get("kappa") for k, r in lad["freqs"].items()},
            "power": {k: r.get("power") for k, r in lad["freqs"].items()},
            "identity_residual": {k: r.get("identity_residual")
                                  for k, r in lad["freqs"].items()},
            "moves_vs_ragged": moves,
        }, indent=2) + "\n")
        print(f"    wrote {out.relative_to(REPO)}", flush=True)
    return lad


def _require_key(ladder, key):
    """§9 item 216: a frequency absent from ``TH14_STEP2_FREQS_MHZ`` skips."""
    if key not in ladder["freqs"]:
        pytest.skip(f"{key} MHz not selected ({th14.FREQS_ENV})")


@complex_only
@needs_solve
@pytest.mark.parametrize("key", FREQ_KEYS)
def test_g_gates_both_columns(ladder, key):
    _require_key(ladder, key)
    rec = ladder["freqs"][key]
    for col, sw in (("cu", rec["sigma"][th14.COPPER_SIGMA]), ("pec", rec["pec"])):
        assert sw["reciprocity"] <= RECIPROCITY_BAND, (col, sw["reciprocity"])
        assert float(np.max(sw["sigma"])) <= 1.0 + PASSIVITY_SIGMA_TOLERANCE, (col, sw["sigma"])
        for name, value in sw["spreads"].items():
            assert value <= ADJACENT_SPREAD_BAND, (col, name, value)


@complex_only
@needs_solve
@pytest.mark.parametrize("key", FREQ_KEYS)
def test_b_copper_surface_loss_identity(ladder, key):
    _require_key(ladder, key)
    if ladder["part"] == "sweeps":
        pytest.skip("(b) is asserted in the identity windows")
    rec = ladder["freqs"][key]
    assert rec["power"]["p_surf"] > 0.0, rec["power"]
    assert rec["identity_residual"] <= DISCRETE_IDENTITY_RTOL, rec["identity_residual"]


@complex_only
@needs_solve
@pytest.mark.parametrize("key", FREQ_KEYS)
def test_sigma_ladder_monotone(ladder, key):
    _require_key(ladder, key)
    if ladder["part"] == "identity":
        pytest.skip("(sigma) is asserted in the sweeps windows")
    rungs = ladder["freqs"][key]["sigma"]
    if MPI.COMM_WORLD.rank == 0:
        print(f"[ANS-6 step4] (sigma) {key} MHz degree {ladder['degree']} per class "
              f"max|S_sigma - S_PEC| over {th14.SIGMA_LADDER}: " + "; ".join(
                  f"{c} " + " > ".join(f"{rungs[sig]['class_ds'][c]:.3e}"
                                       for sig in th14.SIGMA_LADDER) for c in CLASSES),
              flush=True)
    for c in CLASSES:
        seq = [rungs[sig]["class_ds"][c] for sig in th14.SIGMA_LADDER]
        assert all(b < a for a, b in zip(seq, seq[1:])), (c, seq)


@complex_only
@needs_solve
def test_c_bracket_10mhz_against_exact_solid_record(ladder):
    """(c) item 192 landed: per class max|S_cu − S_PEC| ≤ max|S_800 − S_PEC|,
    S_800 = ``ANS4_EXACT_STRIP_S_P1_RECORDS["10"]`` (same exact rectangles)."""
    _require_key(ladder, "10")
    rec = ladder["freqs"]["10"]
    solid = th14._class_max_abs(
        np.asarray(ANS4_EXACT_STRIP_S_P1_RECORDS["10"], dtype=np.complex128) - rec["pec"]["s"])
    ragged = th14._class_max_abs(
        np.asarray(LEG_D_S_MATRIX_10MHZ, dtype=np.complex128) - rec["pec"]["s"])
    cu = rec["sigma"][th14.COPPER_SIGMA]["class_ds"]
    if MPI.COMM_WORLD.rank == 0:
        mismatch = (" [ORDER MISMATCHED: S_800 record is degree 1, S_PEC is degree "
                    f"{ladder['degree']}]" if ladder["degree"] != 1 else "")
        if ladder["resolution"] is not None:
            mismatch += " [MESH MISMATCHED: comparand is ×1]"
        print("[ANS-6 step4] (c) 10 MHz per class max|dS|: " + "; ".join(
            f"{c}: copper {cu[c]:.9e} vs exact solid-800 {solid[c]:.9e} (ragged "
            f"LEG_D, printed only {ragged[c]:.9e})" for c in CLASSES) + mismatch, flush=True)
    for c in CLASSES:
        assert cu[c] <= solid[c], (c, cu[c], solid[c])


@complex_only
@needs_solve
def test_n2_exact_vs_ragged_64mhz_moves(ladder):
    """(N2) floor asserted: every 64 MHz Cu class mean moves > floor."""
    _require_key(ladder, "64")
    if ladder["resolution"] is not None and MPI.COMM_WORLD.rank == 0:
        print("[ANS-6 step5b′] (N2) [MESH MISMATCHED: comparand is ×1] 64 MHz cu moves: "
              + "  ".join(f"{c} {v * 100:.4f}%" for c, v in ladder["moves"]["64"]["cu"].items()),
              flush=True)
    for c, value in ladder["moves"]["64"]["cu"].items():
        assert value > KNOB_REACHES_THE_SOLVE, f"{c} move {value:.3e}"


@complex_only
@needs_solve
def test_n3_degree_knob_reaches_the_solve(ladder):
    """(N3), degree 2 only (§9 item 216): every class mean of each present
    frequency and column moves from the degree-1 record by > floor."""
    if ladder["degree"] != 2:
        pytest.skip("(N3) is a degree-2 control")
    if ladder["resolution"] is not None:
        pytest.skip("(N3) would mix h and order")
    for k, cols in ladder["order_moves"].items():
        for col, mv in cols.items():
            for c, value in mv.items():
                assert value > KNOB_REACHES_THE_SOLVE, f"{k} {col} {c} move {value:.3e}"


@complex_only
@needs_solve
@pytest.mark.parametrize("key", FREQ_KEYS)
def test_r_class_mean_records(ladder, key):
    _require_key(ladder, key)
    if ladder["resolution"] is not None:
        pytest.skip(f"(R) records are ×1 (h = 0.015); not asserted under {RESOLUTION_ENV}")
    if ladder["degree"] == 2:
        record, width = ANS6_EXACT_STRIP_CLASS_MEANS_P2, ANS6_EXACT_STRIP_P2_RECORD_RANK_WIDTH
    else:
        record, width = ANS6_EXACT_STRIP_CLASS_MEANS_P1, ANS6_EXACT_STRIP_RECORD_RANK_WIDTH
    if MPI.COMM_WORLD.size != width:
        pytest.skip(f"(R) at degree {ladder['degree']} is asserted at -n {width} only")
    if record is None or key not in record:
        pytest.skip(f"(R) degree-{ladder['degree']} record for {key} MHz not yet registered")
    for col in ("cu", "pec"):
        for c in CLASSES:
            got = ladder["means"][key][col][c]
            want = record[key][col][c]
            rel = abs(got - want) / abs(want)
            assert rel <= STEP3A_RECORD_RTOL, (key, col, c, got, want, rel)


@complex_only
@needs_solve
@pytest.mark.parametrize("key", FREQ_KEYS)
def test_p_prime_part_reaches_the_ladder(ladder, key):
    """(P′) asserted (§9 item 226): the part knob reaches ``_build_ladder``."""
    if ladder["part"] is None:
        pytest.skip(f"{PART_ENV} unset")
    _require_key(ladder, key)
    rec = ladder["freqs"][key]
    if MPI.COMM_WORLD.rank == 0:
        print(f"[ANS-6 step5] (P') {key} MHz part {ladder['part']}: sigmas "
              f"{sorted(rec['sigma'])}, power {'present' if 'power' in rec else 'absent'}",
              flush=True)
    if ladder["part"] == "sweeps":
        assert set(rec["sigma"]) == set(th14.SIGMA_LADDER), sorted(rec["sigma"])
        assert "power" not in rec
    else:
        assert set(rec["sigma"]) == {th14.COPPER_SIGMA}, sorted(rec["sigma"])
        assert "power" in rec


@complex_only
@needs_solve
def test_f_prime_resolution_reaches_the_ladder(ladder):
    """(F′) ladder level (§9 item 266): knobbed, the solved mesh is refined."""
    if ladder["resolution"] is None:
        pytest.skip(f"{RESOLUTION_ENV} unset")
    assert ladder["cells"] > KNOB_CELL_FACTOR * ANS6_EXACT_STRIP_HOLE_CELL_RECORD, ladder["cells"]
