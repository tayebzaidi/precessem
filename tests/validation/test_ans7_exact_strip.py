"""`ANS-7` step 5 (§9 item 195, chain T17) -- the tuned fixture on exact port and
RLC rectangles at degree 1, κ off: Route A = Route B, rows re-recorded.

Operator ruling 2026-10-02: port and lumped-element sheets are exact CAD
rectangles and the κ width correction is not applied.  ``_build_ports`` builds
the port sheet (P1, 50 Ω) and the three capacitor sheets (P2–P4, Z_C in Route
A) with one call, so the hole-route ``leg_sheet_strip_fraction`` of item 193
makes all four exact.  κ is read only through ``_kappa_in_run`` on uncorrected
specs, printed.

Knobs: ``FEM_EM_ANS7_EXACT_ROUTE`` ∈ {``B``, ``A``} (unset ⇒ every test skips;
anything else raises); ``FEM_EM_ANS7_EXACT_DEGREE`` (default 1; Route A is
degree 1 only).  Writes only under ``logs/ans7_exact/``.  Never calls ``ans:7``'s
``_stage_*`` / ``_assemble`` / ``_write_outputs``.

Anchors (bands imported, never restated): (F) fixture tie, (B) Route B per
node + Δ_LOO, (AB) Route A = Route B at all 102 points, (R) rows record at
``-n 4``; negative control (N) ragged κ-off control (ii), asserted.

§9 item 196 (`ANS-7` step 5b, chain T18): ``FEM_EM_ANS7_EXACT_MHZ`` (a comma
list ⊂ ``DEFAULT_NODES_MHZ``) with ``ROUTE=B`` switches to the per-node split
(item 132's): the listed nodes are solved at ``FEM_EM_ANS7_EXACT_DEGREE``, κ
off, and rank 0 writes ``logs/ans7_exact/node_<mhz>_p<deg>.json``; the
whole-grid ``route_b`` fixture skips.  A subset containing 64.5 MHz runs with
``field_solve=True`` and asserts the loss identity (L).
``FEM_EM_ANS7_EXACT_HELDKAPPA=1`` adds, at 64.5 MHz only, one sweep + P1 field
solve with the ragged degree-2 κ of item 157 -- the negative control, printed,
never asserted, never fed to the rows (``pec_column=True`` there too: item
157's trap forbids ``pec_column=False`` with ``field_solve=True``).  Assembly
(``sparse_rows``) runs once all seven node JSONs exist with one degree and
one cell count.

§9 item 214 (`ANS-7` step 6, chain T19): ``FEM_EM_ANS7_EXACT_RESOLUTION`` (a
float in metres; unset ⇒ nothing forwarded, item 196's path byte-identical)
reaches ``_build_ladder(resolution=…)`` and ``_hole_rung(resolution=…)``;
``FEM_EM_ANS7_EXACT_COLUMNS`` (unset ⇒ Cu and PEC; ``cu`` ⇒ Cu only, no PEC
sweep, no field solve) is the Cu-only five-node path of ``xl`` entry 28.  Node
JSONs at a resolution or Cu-only carry ``_h<res:g>`` / ``_cu`` suffixes, so
item 196's unsuffixed files are never overwritten by a rung run.  The Cu-only
assembly runs on ``ASSEMBLY5_NODES_MHZ``.  (R2) asserts the unset path's rows
against item 196's record; (S5) licenses the five-node Cu subset at ×1 on the
stored JSONs (pure); (N-wire) asserts a rung's cells moved off the ×1 record.
"""

from __future__ import annotations

import json
import math
import os
import resource
import time
from pathlib import Path

import numpy as np
import pytest
from mpi4py import MPI

from tests.complex_mode import complex_only
from tests.mesh.test_birdcage_port_sheet_prerequisite import CELL_COUNT_BAND
from tests.validation.test_ans4_exact_strip_records import SPEC_AREA_M2, SPEC_WIDTH_M
from tests.validation.test_ans6_exact_strip import (
    ANS6_EXACT_STRIP_HOLE_CELL_RECORD,
    _hole_rung,
)
from tests.validation.test_ans7_sparse_route_b import (
    _ANS7,
    ASSEMBLY5_NODES_MHZ,
    BANDS,
    CONTROL_NODES_MHZ,
    DEFAULT_NODES_MHZ,
    GRID_MHZ,
    K_SHARE_FLOOR,
    N_FACTOR,
    NODE_S_RTOL,
    bary_interp,
    sparse_rows,
    z_from_s,
)
from tests.validation.test_ops61_exact_strip_ports import STRIP_BAND
from tests.validation.test_port_birdcage_four_port import PASSIVITY_SIGMA_TOLERANCE
from tests.validation.test_port_circuit_layer_field import REDUCTION_BAND
from tests.validation.test_port_lumped_narrowed_sheet import GATED_WIDTH_FRACTION
from tests.validation.test_port_lumped_sheet_sweep import RECIPROCITY_BAND

route_b_s = _ANS7.route_b_s
route_b_z = _ANS7.route_b_z
_route_a_solve = _ANS7._route_a_solve
_kappa_in_run = _ANS7._kappa_in_run
_s11_from_z = _ANS7._s11_from_z
parallel_rows = _ANS7.parallel_rows
C_F = _ANS7.C_F
ANCHOR1_HZ = _ANS7.ANCHOR1_FREQUENCY_HZ

REPO = Path(__file__).resolve().parents[2]
CASE_DIR = REPO / "examples" / "ansys_benchmarks" / "ans7_tuned_copper_birdcage_driven_sweep"
OUT_DIR = REPO / "logs" / "ans7_exact"
ROUTE_ENV = "FEM_EM_ANS7_EXACT_ROUTE"
DEGREE_ENV = "FEM_EM_ANS7_EXACT_DEGREE"
ROWS = ("f0_hz", "r_in_ohm", "q_loaded")
PORTS = (1, 2, 3, 4)
RECORD_RANKS = 4

# Ragged κ-corrected degree-1 rows (metrics.json "rows", our numbers) -- the
# printed comparand of the predictions.
RAGGED_ROWS_P1 = {
    "cu": {"f0_hz": 63983031.49581605, "r_in_ohm": 16998.7066800671, "q_loaded": 918.440964550853},
    "pec": {
        "f0_hz": 63997103.30801084,
        "r_in_ohm": 28425.152089802326,
        "q_loaded": 1535.9515756955145,
    },
}

# (R) record -- Route B rows on the exact-rectangle hole mesh, κ off, degree 1.
# Registered from W1 (`20261003T112739Z_ANS-7-item195-W1-routeB.log:2258–2262`,
# 80 798 cells).  Version tag: complex 0.11 image, -n 4; (1*) licence --
# re-record on an image change, never re-band.
ANS7_EXACT_STRIP_ROUTE_B_ROWS_P1 = {
    "cu": {
        "f0_hz": 63971137.8564212,
        "r_in_ohm": 16907.699520475457,
        "q_loaded": 912.7173367759675,
    },
    "pec": {
        "f0_hz": 63985193.342905395,
        "r_in_ohm": 28088.54262947128,
        "q_loaded": 1516.4394852018208,
    },
}


# §9 item 196 knobs and printed comparands.
MHZ_ENV = "FEM_EM_ANS7_EXACT_MHZ"
HELDKAPPA_ENV = "FEM_EM_ANS7_EXACT_HELDKAPPA"
LOSS_NODE_MHZ = 64.5
# The ragged degree-2 κ of item 157 (its W0 log `:1124`, as cited and restated
# by item 196's brief) -- negative control only.
HELD_KAPPA_P2_RAGGED = 5.2094420914264195e-02
# Item 157's ragged in-run p2 rows (printed comparand only).  Cu at full repr
# (`20261002T005913Z_ANS-7-step3iii-W2.log:1325-1327`, imported constant's
# digits); PEC as quoted in item 196's brief.
ITEM157_P2_RAGGED_ROWS = {
    "cu": {
        "f0_hz": 62704215.74044876,
        "r_in_ohm": 18569.60744330647,
        "q_loaded": 984.2226364364154,
    },
    "pec": {"f0_hz": 62.717977e6, "r_in_ohm": 32545.0, "q_loaded": 1724.8},
}


# §9 item 214 knobs (with _ROUTE / _DEGREE / _MHZ: xl entry 28's contract).
RESOLUTION_ENV = "FEM_EM_ANS7_EXACT_RESOLUTION"
COLUMNS_ENV = "FEM_EM_ANS7_EXACT_COLUMNS"
BOTH_COLUMNS = ("cu", "pec")

# (R2) record -- the degree-2 seven-node Route B rows on the exact-rectangle
# hole mesh, κ off, ×1 (h = 0.015).  Registered from item 196 W3
# (`20261003T124814Z_ANS-7-item196-W3-p2.log:2189,2193`, full repr, 80 798
# cells).  Version tag: complex 0.11 image, -n 8; (1*) licence -- re-record on
# an image change, never re-band.
ANS7_EXACT_STRIP_ROUTE_B_ROWS_P2 = {
    "cu": {
        "f0_hz": 62701081.93936459,
        "r_in_ohm": 18441.754600739256,
        "q_loaded": 977.2613699232847,
    },
    "pec": {
        "f0_hz": 62714806.5966287,
        "r_in_ohm": 32123.932870177392,
        "q_loaded": 1701.733130785279,
    },
}


def _resolution():
    raw = os.environ.get(RESOLUTION_ENV, "").strip()
    if not raw:
        return None
    return float(raw)


def _cols():
    raw = os.environ.get(COLUMNS_ENV, "").strip()
    if not raw:
        return BOTH_COLUMNS
    if raw == "cu":
        return ("cu",)
    raise ValueError(f"{COLUMNS_ENV}={raw!r}: only unset (cu + pec) or 'cu'")


def _ceil_1sf(x):
    """Ceiling to one significant figure (B₅ = ceil_1sf(max(3·Δ_LOO, 1e-6)))."""
    if not math.isfinite(x) or x <= 0:
        return x
    e = math.floor(math.log10(x))
    m = round(x / 10.0**e, 9)
    return math.ceil(m) * 10.0**e


def _nodes_env():
    raw = os.environ.get(MHZ_ENV, "").strip()
    if not raw:
        return None
    nodes = tuple(round(float(v), 6) for v in raw.split(",") if v.strip())
    bad = [m for m in nodes if m not in DEFAULT_NODES_MHZ]
    if bad or not nodes:
        raise ValueError(f"{MHZ_ENV}={raw!r}: nodes must be a subset of {DEFAULT_NODES_MHZ}")
    return nodes


def _node_json(mhz, deg, res=None, cols=BOTH_COLUMNS):
    """Unset resolution with both columns keeps item 196's names; otherwise
    ``node_<mhz>_p<deg>[_h<res:g>][_cu].json``."""
    suffix = ""
    if res is not None:
        suffix += f"_h{res:g}"
    if tuple(cols) == ("cu",):
        suffix += "_cu"
    return OUT_DIR / f"node_{mhz:.1f}_p{deg}{suffix}.json"


def _route():
    raw = os.environ.get(ROUTE_ENV, "").strip()
    if not raw:
        return None
    if raw not in ("A", "B"):
        raise ValueError(f"{ROUTE_ENV}={raw!r}: only 'A' or 'B'")
    return raw


def _degree():
    return int(os.environ.get(DEGREE_ENV, "1") or 1)


def _say(msg):
    if MPI.COMM_WORLD.rank == 0:
        print(msg, flush=True)


def _need(route):
    r = _route()
    if r != route:
        pytest.skip(f"{ROUTE_ENV}={r!r}: this test runs on route {route}")


def _cj(m):
    return [
        [[float(complex(v).real), float(complex(v).imag)] for v in row] for row in np.asarray(m)
    ]


def _uj(m):
    return np.array([[complex(v[0], v[1]) for v in row] for row in m])


def _route_b_json(deg):
    return OUT_DIR / f"route_b_p{deg}.json"


def _fixture_tie(res=None):
    """(F): cells and the four N_i rectangles on the strip hole mesh (reduced)."""
    r = _hole_rung(GATED_WIDTH_FRACTION, resolution=res)
    for i in PORTS:
        p = r["ports"][i]
        _say(
            f"[ans7-exact] F P{i}: N={p['N']:.12e} (rel {p['N'] / SPEC_AREA_M2 - 1:+.3e}) "
            f"w=A/h={p['N'] / p['N_h']:.12e} (rel {p['N'] / p['N_h'] / SPEC_WIDTH_M - 1:+.3e})"
        )
    if res is None:
        _say(
            f"[ans7-exact] F cells {r['cells']} (PREDICTED = record "
            f"{ANS6_EXACT_STRIP_HOLE_CELL_RECORD})"
        )
    else:
        _say(
            f"[ans7-exact] F cells {r['cells']} at h = {res:g} m (record x1 "
            f"{ANS6_EXACT_STRIP_HOLE_CELL_RECORD}, ratio "
            f"{r['cells'] / ANS6_EXACT_STRIP_HOLE_CELL_RECORD:.4f}; (N-wire) ASSERTED to differ "
            f"by > CELL_COUNT_BAND {CELL_COUNT_BAND:g})"
        )
    return r


def _crossings(f_hz, z_in):
    b = (1.0 / np.asarray(z_in)).imag
    return sum(1 for k in range(len(b) - 1) if b[k] < 0 <= b[k + 1])


# -- Route B ------------------------------------------------------------------


def _loo(nodes, z_nodes, columns=BOTH_COLUMNS, c_f=C_F, grid_mhz=GRID_MHZ):
    # ``grid_mhz`` (§9 item 228, additive; default ``GRID_MHZ`` byte-identical):
    # forwarded to both ``sparse_rows`` calls so a 127.74 MHz node set is read on
    # its own grid, never extrapolated onto 62–67 MHz.
    full = sparse_rows(nodes, z_nodes, grid_mhz=grid_mhz, columns=columns, c_f=c_f)
    worst = {k: 0.0 for k in ROWS}
    for i in range(1, len(nodes) - 1):
        keep = [m for j, m in enumerate(nodes) if j != i]
        sub = {c: [z_nodes[c][j] for j in range(len(nodes)) if j != i] for c in z_nodes}
        loo = sparse_rows(keep, sub, grid_mhz=grid_mhz, columns=columns, c_f=c_f)
        for c in columns:
            for k in ROWS:
                if k not in loo[c]["rows"] or k not in full[c]["rows"]:
                    worst[k] = np.inf
                    continue
                d = abs(loo[c]["rows"][k] - full[c]["rows"][k]) / abs(full[c]["rows"][k])
                worst[k] = max(worst[k], d)
    return full, worst


@pytest.fixture(scope="module")
def route_b():
    from tests.validation import test_th14_birdcage_copper as th14

    _need("B")
    if _nodes_env() is not None:
        pytest.skip(f"{MHZ_ENV} set: per-node split (item 196) -- see test_p2_*")
    if _cols() != BOTH_COLUMNS:
        raise ValueError(f"{COLUMNS_ENV}=cu needs the per-node split ({MHZ_ENV})")
    comm = MPI.COMM_WORLD
    deg = _degree()
    res = _resolution()
    nodes = DEFAULT_NODES_MHZ
    tie = _fixture_tie(res)
    t0 = time.perf_counter()
    lad = th14._build_ladder(
        deg,
        sigmas=(th14.COPPER_SIGMA,),
        resolution=res,
        frequencies_hz=[m * 1e6 for m in nodes],
        width_correction_kappa=None,
        measure_kappa=False,
        pec_column=True,
        field_solve=False,
        leg_sheet_strip_fraction=GATED_WIDTH_FRACTION,
    )
    comm.Barrier()
    t_solve = time.perf_counter() - t0
    recs = sorted(lad["freqs"].values(), key=lambda r: r["f_hz"])
    assert len(recs) == len(nodes)
    data = {
        "nodes": nodes,
        "s": {"cu": [], "pec": []},
        "z": {"cu": [], "pec": []},
        "recip": {"cu": [], "pec": []},
        "sigma_max": {"cu": [], "pec": []},
        "cells": int(lad["cells"]),
        "tie": tie,
        "deg": deg,
        "res": res,
    }
    for mhz, rec in zip(nodes, recs):
        assert abs(rec["f_hz"] - mhz * 1e6) < 1.0
        for c, sw in (("cu", rec["sigma"][th14.COPPER_SIGMA]), ("pec", rec["pec"])):
            s4 = np.asarray(sw["s"], dtype=np.complex128)
            data["s"][c].append(s4)
            data["z"][c].append(z_from_s(s4))
            data["recip"][c].append(float(sw["reciprocity"]))
            data["sigma_max"][c].append(float(np.max(sw["sigma"])))
    full, worst = _loo(list(nodes), data["z"])
    data["full"], data["worst"] = full, worst
    for c in ("cu", "pec"):
        r = full[c]["rows"]
        _say(
            f"[ans7-exact] B p{deg} {c} rows: f0 {r.get('f0_hz')!r} Hz, R_in "
            f"{r.get('r_in_ohm')!r} ohm, Q {r.get('q_loaded')!r}, k-share "
            f"{r.get('k_share')!r}; up-crossings {r.get('crossings_up_hz')}"
        )
        for k in ROWS:
            if k in r:
                rag = RAGGED_ROWS_P1[c][k]
                _say(
                    f"    {c} {k}: exact {r[k]:.10e} vs ragged κ-corrected p1 {rag:.10e} "
                    f"(move {r[k] / rag - 1:+.4%})"
                )
    for k in ROWS:
        _say(
            f"[ans7-exact] B Delta_LOO {k} = {worst[k]:.3e} (band {BANDS[k]:.0e}; "
            f"PREDICTED <= 1e-9)"
        )
    _say(
        f"[ans7-exact] B {len(nodes)} nodes in {t_solve:.1f} s, cells {data['cells']}, "
        f"-n {comm.size}"
    )
    if comm.rank == 0:
        OUT_DIR.mkdir(parents=True, exist_ok=True)
        _route_b_json(deg).write_text(
            json.dumps(
                {
                    "degree": deg,
                    "ranks": comm.size,
                    "cells": data["cells"],
                    "nodes_mhz": nodes,
                    "s": {c: [_cj(m) for m in data["s"][c]] for c in data["s"]},
                    "z": {c: [_cj(m) for m in data["z"][c]] for c in data["z"]},
                    "rows": {c: {k: v for k, v in full[c]["rows"].items()} for c in full},
                    "delta_loo": worst,
                }
            )
        )
    comm.Barrier()
    return data


def _assert_tie(cells, tie, res=None):
    """(F).  Unset resolution: today's tie, unchanged.  At a resolution: the
    tie measured the mesh the ladder solved (cells equal exactly) and (N-wire)
    -- the cells moved off the ×1 record by more than ``CELL_COUNT_BAND``; the
    CAD rectangles are independent of h."""
    if res is None:
        assert abs(cells / ANS6_EXACT_STRIP_HOLE_CELL_RECORD - 1.0) < CELL_COUNT_BAND, cells
        assert abs(tie["cells"] / ANS6_EXACT_STRIP_HOLE_CELL_RECORD - 1.0) < CELL_COUNT_BAND
    else:
        assert int(tie["cells"]) == int(cells), (tie["cells"], cells)
        assert abs(cells / ANS6_EXACT_STRIP_HOLE_CELL_RECORD - 1.0) > CELL_COUNT_BAND, (
            "(N-wire)",
            res,
            cells,
        )
    for i in PORTS:
        p = tie["ports"][i]
        assert abs(p["N"] / SPEC_AREA_M2 - 1.0) < STRIP_BAND, (i, p["N"])
        assert abs(p["N"] / p["N_h"] / SPEC_WIDTH_M - 1.0) < STRIP_BAND, (i, p["N"] / p["N_h"])


@complex_only
def test_f_fixture_tie_route_b(route_b):
    _need("B")
    _assert_tie(route_b["cells"], route_b["tie"], route_b["res"])


@complex_only
def test_b_route_b_per_node(route_b):
    _need("B")
    for i, mhz in enumerate(route_b["nodes"]):
        for c in ("cu", "pec"):
            _say(
                f"[ans7-exact] B node {mhz:.1f} {c}: recip {route_b['recip'][c][i]:.3e} "
                f"sigma_max {route_b['sigma_max'][c][i]:.12f}"
            )
            assert route_b["recip"][c][i] <= RECIPROCITY_BAND, (mhz, c)
            assert route_b["sigma_max"][c][i] <= 1.0 + PASSIVITY_SIGMA_TOLERANCE, (mhz, c)
    for c in ("cu", "pec"):
        r = route_b["full"][c]["rows"]
        assert len(r["crossings_up_hz"]) == 1, (c, r["crossings_up_hz"])
        assert r["k_share"] > K_SHARE_FLOOR, (c, r["k_share"])
    for k in ROWS:
        assert route_b["worst"][k] <= BANDS[k], (k, route_b["worst"][k])


@complex_only
def test_r_route_b_rows_record(route_b):
    _need("B")
    if ANS7_EXACT_STRIP_ROUTE_B_ROWS_P1 is None:
        pytest.skip("(R) record not yet registered (W1 measures it)")
    if MPI.COMM_WORLD.size != RECORD_RANKS or route_b["deg"] != 1:
        pytest.skip(f"(R) asserted at -n {RECORD_RANKS}, degree 1 only")
    for c in ("cu", "pec"):
        for k in ROWS:
            got = route_b["full"][c]["rows"][k]
            want = ANS7_EXACT_STRIP_ROUTE_B_ROWS_P1[c][k]
            assert abs(got / want - 1.0) <= BANDS[k], (c, k, got, want)


# -- Route A ------------------------------------------------------------------


def _exact_geo(th14, comm, strip, res=None):
    """The hole mesh + uncorrected ports for Route A (lifted from ``route_a``'s
    inner ``_geo``, unchanged, so §9 item 217's ``tune_inmodel`` reuses it).
    ``res`` (§9 item 310, additive; ``None`` forwards nothing, byte-identical)
    reaches ``_sheets_build``'s ``resolution`` keyword, as ``_hole_rung`` does."""
    extra = {} if strip is None else {"leg_sheet_strip_fraction": strip}
    if res is not None:
        extra["resolution"] = float(res)
    msh, cell_tags, facet_tags, _d, _t = th14._sheets_build(True, as_hole=True, **extra)
    tdim = msh.topology.dim
    msh.topology.create_connectivity(tdim - 1, tdim)
    msh.topology.create_entity_permutations()
    cells = int(msh.topology.index_map(tdim).size_global)
    outer_tags, _census = th14._outer_box_tags(msh, facet_tags, comm)
    tags_f, port_defs, specs = th14._build_ports(msh, cell_tags, comm, width_correction_kappa=None)
    return (msh, cell_tags, facet_tags, outer_tags, tags_f, port_defs), specs, cells


@pytest.fixture(scope="module")
def route_a():
    from tests.validation import test_th14_birdcage_copper as th14

    _need("A")
    comm = MPI.COMM_WORLD
    deg = _degree()
    if deg != 1:
        pytest.skip("Route A is degree 1 only")
    jb = _route_b_json(1)
    if not jb.exists():
        pytest.fail(f"{jb} absent: run Route B (W1) first")
    rb = json.loads(jb.read_text())
    tie = _fixture_tie()

    def _geo(strip):
        return _exact_geo(th14, comm, strip)

    t0 = time.perf_counter()
    geo, specs, cells = _geo(GATED_WIDTH_FRACTION)
    kappa_exact = _kappa_in_run(th14, geo, specs, comm)
    _say(
        f"[ans7-exact] A kappa_exact (64 MHz Cu, uncorrected, printed) = {kappa_exact!r} "
        f"(PREDICTED [5e-6, 5e-5]: {'met' if 5e-6 <= kappa_exact <= 5e-5 else 'NOT met'})"
    )
    xn = np.asarray(rb["nodes_mhz"], dtype=float) - 64.5
    z_nodes = {c: [_uj(m) for m in rb["z"][c]] for c in ("cu", "pec")}
    res = {"cu": [], "pec": []}
    z_in_a = {"cu": [], "pec": []}
    for g in GRID_MHZ:
        f = g * 1e6
        for c in ("cu", "pec"):
            s_a, _ = _route_a_solve(th14, geo, specs, f, c)
            zb = route_b_z(bary_interp(xn, z_nodes[c], g - 64.5), f, C_F)[1]
            s_b = _s11_from_z(zb)
            res[c].append(float(abs(s_a - s_b) / abs(s_b)))
            z_bare_a = _ANS7._z_from_s11(s_a)
            z_in_a[c].append(_ANS7.add_c1(z_bare_a, 2.0 * np.pi * f, C_F))
    t_sweep = time.perf_counter() - t0
    for c in ("cu", "pec"):
        for g, r in zip(GRID_MHZ, res[c]):
            _say(f"[ans7-exact] AB {c} {g:.1f} MHz: |S11_A - S11_B|/|S11_B| = {r:.3e}")
        _say(
            f"[ans7-exact] AB {c} max {max(res[c]):.3e} min {min(res[c]):.3e} "
            f"(band {REDUCTION_BAND:.0e}; PREDICTED [5e-6, 3e-4] every point)"
        )
    rows_a = {
        c: parallel_rows(np.array(GRID_MHZ) * 1e6, np.array(z_in_a[c])) for c in ("cu", "pec")
    }
    for c in ("cu", "pec"):
        r = rows_a[c]
        _say(
            f"[ans7-exact] A {c} rows: f0 {r.get('f0_hz')!r} R_in {r.get('r_in_ohm')!r} "
            f"Q {r.get('q_loaded')!r}; up-crossings {r['crossings_up_hz']}"
        )
        for k in ROWS:
            if k in r:
                rag = RAGGED_ROWS_P1[c][k]
                _say(
                    f"    {c} {k}: exact A {r[k]:.10e} vs ragged p1 {rag:.10e} "
                    f"(move {r[k] / rag - 1:+.4%})"
                )
    # (N) ragged control (ii), verbatim: default hole mesh, uncorrected, 64 MHz Cu.
    geo_r, specs_r, _cells_r = _geo(None)
    s11_unc_a, _ = _route_a_solve(th14, geo_r, specs_r, ANCHOR1_HZ, "cu")
    tracked = json.loads((CASE_DIR / "metrics.json").read_text())
    s4_unc = _ANS7._unmat(tracked["check_64mhz_uncorrected"]["cu_s"])
    _zin, zb_unc = route_b_s(s4_unc, ANCHOR1_HZ)
    s_b_unc = _s11_from_z(zb_unc)
    ctrl = float(abs(s11_unc_a - s_b_unc) / abs(s_b_unc))
    j64 = list(np.round(GRID_MHZ, 6)).index(64.0)
    _say(
        f"[ans7-exact] N ragged kappa-off control (ii) = {ctrl:.6e} = "
        f"{ctrl / REDUCTION_BAND:.2f} x REDUCTION_BAND (ASSERTED > {N_FACTOR:g} x; "
        f"backing 1.657e-02); exact-sheet 64.0 MHz Cu (AB) = {res['cu'][j64]:.3e}"
    )
    _say(
        f"[ans7-exact] A total {time.perf_counter() - t0:.1f} s (sweep+kappa {t_sweep:.1f} s), "
        f"cells {cells}, -n {comm.size}"
    )
    if comm.rank == 0:
        OUT_DIR.mkdir(parents=True, exist_ok=True)
        (OUT_DIR / "route_a_p1.json").write_text(
            json.dumps(
                {
                    "ranks": comm.size,
                    "cells": cells,
                    "kappa_exact": kappa_exact,
                    "residual": res,
                    "control_ii": ctrl,
                    "rows_a": rows_a,
                }
            )
        )
    comm.Barrier()
    return {"cells": cells, "tie": tie, "res": res, "rows_a": rows_a, "ctrl": ctrl}


@complex_only
def test_f_fixture_tie_route_a(route_a):
    _need("A")
    _assert_tie(route_a["cells"], route_a["tie"])


@complex_only
def test_ab_route_a_equals_route_b(route_a):
    _need("A")
    n = 0
    for c in ("cu", "pec"):
        assert len(route_a["res"][c]) == len(GRID_MHZ)
        for g, r in zip(GRID_MHZ, route_a["res"][c]):
            assert r <= REDUCTION_BAND, (c, g, r)
            n += 1
        assert len(route_a["rows_a"][c]["crossings_up_hz"]) == 1, (c, route_a["rows_a"][c])
    assert n == 2 * len(GRID_MHZ)


@complex_only
def test_n_ragged_kappa_off_control(route_a):
    _need("A")
    assert route_a["ctrl"] > N_FACTOR * REDUCTION_BAND, route_a["ctrl"]


# -- §9 item 196: per-node Route B split (degree 2), (L) at 64.5 MHz -----------


def _tie_json(tie):
    return {
        "cells": int(tie["cells"]),
        "ports": {
            str(i): {"N": float(tie["ports"][i]["N"]), "N_h": float(tie["ports"][i]["N_h"])}
            for i in PORTS
        },
    }


def _tie_from_json(t):
    return {"cells": t["cells"], "ports": {int(k): v for k, v in t["ports"].items()}}


def _im_y_in(z4, mhz):
    z_in = route_b_z(z4, mhz * 1e6, C_F)[0]
    return float((1.0 / z_in).imag), z_in


@pytest.fixture(scope="module")
def p2_nodes():
    from tests.validation import test_th14_birdcage_copper as th14
    from tests.validation.test_ans7_sparse_route_b import _column_power, _loss_pair, _p_inc_w

    _need("B")
    nodes = _nodes_env()
    if nodes is None:
        pytest.skip(f"{MHZ_ENV} unset: whole-grid route_b fixture instead")
    comm = MPI.COMM_WORLD
    deg = _degree()
    res = _resolution()
    cols = _cols()
    pec_column = "pec" in cols
    # Item 157's trap: never pec_column=False with field_solve=True.
    field = (LOSS_NODE_MHZ in nodes) and pec_column
    held = os.environ.get(HELDKAPPA_ENV, "").strip() == "1"
    if held and not pec_column:
        raise ValueError(f"{HELDKAPPA_ENV}=1 with {COLUMNS_ENV}=cu")
    if held and not field:
        raise ValueError(f"{HELDKAPPA_ENV}=1 needs {LOSS_NODE_MHZ} MHz in {MHZ_ENV}")
    tie = _fixture_tie(res)
    t0 = time.perf_counter()
    lad = th14._build_ladder(
        deg,
        sigmas=(th14.COPPER_SIGMA,),
        resolution=res,
        frequencies_hz=[m * 1e6 for m in nodes],
        width_correction_kappa=None,
        measure_kappa=False,
        pec_column=pec_column,
        field_solve=field,
        leg_sheet_strip_fraction=GATED_WIDTH_FRACTION,
    )
    comm.Barrier()
    t_solve = time.perf_counter() - t0
    rss_gib = (
        comm.allreduce(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, op=MPI.SUM) / 1024**2
    )
    recs = sorted(lad["freqs"].values(), key=lambda r: r["f_hz"])
    assert len(recs) == len(nodes)
    cells = int(lad["cells"])
    unknowns = int(lad["unknowns"])
    out = {
        "nodes": nodes,
        "deg": deg,
        "cells": cells,
        "tie": tie,
        "field": field,
        "per": {},
        "held": None,
        "res": res,
        "cols": cols,
        "unknowns": unknowns,
        "rss_gib": rss_gib,
    }
    if comm.rank == 0:
        OUT_DIR.mkdir(parents=True, exist_ok=True)
    for mhz, rec in zip(nodes, recs):
        assert abs(rec["f_hz"] - mhz * 1e6) < 1.0
        per = {"s": {}, "z": {}, "recip": {}, "sigma_max": {}, "im_y": {}, "colpow": {}}
        pairs = [("cu", rec["sigma"][th14.COPPER_SIGMA])]
        if pec_column:
            pairs.append(("pec", rec["pec"]))
        per["t_cu_s"] = float(rec["sigma"][th14.COPPER_SIGMA]["t"])
        per["t_pec_s"] = float(rec["t_pec"]) if pec_column else None
        _say(
            f"[ans7-exact] p{deg} node {mhz:.1f} MHz: Cu sweep {per['t_cu_s']:.1f} s"
            + (f", PEC sweep {per['t_pec_s']:.1f} s" if pec_column else "")
            + f" (h {res if res is not None else 'x1'}, columns {cols}, -n {comm.size})"
        )
        for c, sw in pairs:
            s4 = np.asarray(sw["s"], dtype=np.complex128)
            per["s"][c], per["z"][c] = s4, z_from_s(s4)
            per["recip"][c] = float(sw["reciprocity"])
            per["sigma_max"][c] = float(np.max(sw["sigma"]))
            per["im_y"][c], z_in = _im_y_in(per["z"][c], mhz)
            per["colpow"][c] = _column_power(s4)
            _say(
                f"[ans7-exact] p{deg} node {mhz:.1f} MHz {c}: Im Y_in {per['im_y'][c]:+.9e} S "
                f"(Z_in {z_in.real:.6e}{z_in.imag:+.6e}j ohm); recip {per['recip'][c]:.3e}; "
                f"sigma_max {per['sigma_max'][c]:.12f}; column power {per['colpow'][c]:.9e}"
            )
        per["power"] = dict(rec["power"]) if field else None
        per["identity"] = float(rec["identity_residual"]) if field else None
        out["per"][mhz] = per
        if comm.rank == 0:
            _node_json(mhz, deg, res, cols).write_text(
                json.dumps(
                    {
                        "f_hz": rec["f_hz"],
                        "degree": deg,
                        "ranks": comm.size,
                        "cells": cells,
                        "tie": _tie_json(tie),
                        "kappa": None,
                        "field_solve": field,
                        "s": {c: _cj(per["s"][c]) for c in per["s"]},
                        "z": {c: _cj(per["z"][c]) for c in per["z"]},
                        "recip": per["recip"],
                        "sigma_max": per["sigma_max"],
                        "im_y_in": per["im_y"],
                        "power": per["power"],
                        "identity_residual": per["identity"],
                        "resolution": res,
                        "columns": list(cols),
                        "unknowns": unknowns,
                        "t_cu_s": per["t_cu_s"],
                        "t_pec_s": per["t_pec_s"],
                        "summed_ru_maxrss_gib": rss_gib,
                    }
                )
            )
    comm.Barrier()
    _say(
        f"[ans7-exact] p{deg} nodes {nodes} solved in {t_solve:.1f} s (field_solve {field}), "
        f"cells {cells}, unknowns {unknowns}, h {res if res is not None else 'x1'}, columns "
        f"{cols}, -n {comm.size}; summed ru_maxrss {rss_gib:.2f} GiB; JSON in {OUT_DIR}"
    )
    if field:
        p = out["per"][LOSS_NODE_MHZ]
        out["L"] = _loss_pair(p["s"]["cu"], p["power"])
        out["p_inc"] = _p_inc_w()
    del lad, recs
    if held:
        import gc

        gc.collect()
        t1 = time.perf_counter()
        hl = th14._build_ladder(
            deg,
            sigmas=(th14.COPPER_SIGMA,),
            frequencies_hz=[LOSS_NODE_MHZ * 1e6],
            width_correction_kappa=HELD_KAPPA_P2_RAGGED,
            measure_kappa=False,
            pec_column=True,
            field_solve=True,
            leg_sheet_strip_fraction=GATED_WIDTH_FRACTION,
        )
        comm.Barrier()
        hrec = list(hl["freqs"].values())[0]
        s_h = np.asarray(hrec["sigma"][th14.COPPER_SIGMA]["s"], dtype=np.complex128)
        out["held"] = {
            "L": _loss_pair(s_h, dict(hrec["power"])),
            "colpow": _column_power(s_h),
            "power": dict(hrec["power"]),
            "t_s": time.perf_counter() - t1,
        }
        del hl, hrec
    return out


@complex_only
def test_p2_fixture_tie(p2_nodes):
    _need("B")
    _assert_tie(p2_nodes["cells"], p2_nodes["tie"], p2_nodes["res"])


@complex_only
def test_p2_node_anchors(p2_nodes):
    """Per node: reciprocity, passivity; at 64.5 MHz the copper identity residual."""
    from tests.validation import test_th14_birdcage_copper as th14

    _need("B")
    for mhz, p in p2_nodes["per"].items():
        for c in p2_nodes["cols"]:
            assert p["recip"][c] <= RECIPROCITY_BAND, (mhz, c, p["recip"][c])
            assert p["sigma_max"][c] <= 1.0 + PASSIVITY_SIGMA_TOLERANCE, (mhz, c)
        if p["identity"] is not None:
            _say(
                f"[ans7-exact] p{p2_nodes['deg']} node {mhz:.1f} MHz copper identity residual "
                f"{p['identity']:.3e} (ASSERTED <= {th14.DISCRETE_IDENTITY_RTOL:g})"
            )
            assert p["identity"] <= th14.DISCRETE_IDENTITY_RTOL, (mhz, p["identity"])


@complex_only
def test_p2_window_edge_guard(p2_nodes):
    """Item 157's guard: Im Y_in at 62.0 MHz must be < 0 in every solved column."""
    _need("B")
    if 62.0 not in p2_nodes["per"]:
        pytest.skip("62.0 MHz not in this subset")
    for c in p2_nodes["cols"]:
        b = p2_nodes["per"][62.0]["im_y"][c]
        _say(
            f"[ans7-exact] window-edge guard {c}: Im Y_in(62.0 MHz) = {b:+.9e} S "
            f"(stop before W2 if > 0)"
        )
        assert b < 0.0, (c, b)


@complex_only
def test_p2_loss_identity(p2_nodes):
    """(L) at 64.5 MHz, κ off: |L_S - L_F| <= REDUCTION_BAND.  Held-κ control printed."""
    _need("B")
    if not p2_nodes["field"]:
        pytest.skip(f"(L) needs {LOSS_NODE_MHZ} MHz in the subset")
    deg = p2_nodes["deg"]
    p = p2_nodes["per"][LOSS_NODE_MHZ]
    l_s, l_f = p2_nodes["L"]
    gap = abs(l_s - l_f)
    cp_cu, cp_pec = p["colpow"]["cu"], p["colpow"]["pec"]
    _say(
        f"[ans7-exact] (L) {LOSS_NODE_MHZ} MHz p{deg}, kappa off: column power Cu {cp_cu:.9e}, "
        f"PEC {cp_pec:.9e}; L_S {l_s:.9e}, L_F {l_f:.9e} (P_inc {p2_nodes['p_inc']:.3e} W; "
        f"power {p['power']}); |L_S - L_F| = {gap:.6e} (ASSERTED <= {REDUCTION_BAND:.0e})"
    )
    _say(
        f"[ans7-exact] PREDICTED (printed): (L) gap in [1e-7, 1e-4] -> {gap:.3e} "
        f"({'met' if 1e-7 <= gap <= 1e-4 else 'NOT met'}); Cu column power in "
        f"[0.995, 0.9995] -> {cp_cu:.6f} ({'met' if 0.995 <= cp_cu <= 0.9995 else 'NOT met'})"
    )
    h = p2_nodes["held"]
    if h is not None:
        hs, hf = h["L"]
        hgap = abs(hs - hf)
        _say(
            f"[ans7-exact] control (printed, never asserted) held ragged kappa_deg2 "
            f"{HELD_KAPPA_P2_RAGGED!r} on exact sheets, {LOSS_NODE_MHZ} MHz p{deg}: column power "
            f"Cu {h['colpow']:.9e}; L_S {hs:.9e}, L_F {hf:.9e} (power {h['power']}); "
            f"|L_S - L_F| = {hgap:.6e} (PREDICTED [2e-2, 2e-1]: "
            f"{'met' if 2e-2 <= hgap <= 2e-1 else 'NOT met'}); separation from (L) "
            f"{hgap / gap if gap > 0 else float('inf'):.3e} x; {h['t_s']:.1f} s"
        )
    assert gap <= REDUCTION_BAND, (gap, l_s, l_f, cp_cu, cp_pec)


def _required_nodes(cols):
    return DEFAULT_NODES_MHZ if tuple(cols) == BOTH_COLUMNS else ASSEMBLY5_NODES_MHZ


def _load_nodes(nodes, deg, res, cols):
    """The node JSONs at (degree, resolution, columns), or the missing names."""
    paths = [_node_json(m, deg, res, cols) for m in nodes]
    missing = [p.name for p in paths if not p.exists()]
    if missing:
        return None, missing
    return [json.loads(p.read_text()) for p in paths], []


def _control_report(tag, set_nodes, z_nodes, full, cols):
    """(N3), predicted and printed: the 3-node control against the set's rows."""
    idx = [list(set_nodes).index(m) for m in CONTROL_NODES_MHZ]
    sub = {c: [z_nodes[c][i] for i in idx] for c in cols}
    ctl = sparse_rows(CONTROL_NODES_MHZ, sub, columns=cols)
    for c in cols:
        parts = []
        miss = False
        for k in ROWS:
            a, b = ctl[c]["rows"].get(k), full[c]["rows"].get(k)
            d = math.inf if a is None or b is None else abs(a - b) / abs(b)
            if k != "f0_hz" and d > BANDS[k]:
                miss = True
            parts.append(f"{k} {d:.3e}")
        _say(
            f"[ans7-exact] {tag} (N3) 3-node control {CONTROL_NODES_MHZ} {c} vs the "
            f"{len(set_nodes)}-node rows: {', '.join(parts)} (band {BANDS['r_in_ohm']:.0e}; "
            f"PREDICTED miss on R_in or Q -> {'met' if miss else 'NOT met'}; printed, never "
            f"asserted)"
        )


def _b5_report(tag, worst):
    for k in ROWS:
        _say(
            f"[ans7-exact] {tag} Delta_LOO {k} = {worst[k]:.3e}; B5 = ceil_1sf(max(3*Delta_LOO, "
            f"1e-6)) = {_ceil_1sf(max(3.0 * worst[k], 1e-6)):.0e} (printed)"
        )


@pytest.fixture(scope="module")
def p2_asm(p2_nodes):
    """The assembly on this run's (degree, resolution, columns) node JSONs."""
    deg, res, cols = p2_nodes["deg"], p2_nodes["res"], p2_nodes["cols"]
    set_nodes = _required_nodes(cols)
    recs, missing = _load_nodes(set_nodes, deg, res, cols)
    if recs is None:
        return {"missing": missing}
    z_nodes = {c: [_uj(r["z"][c]) for r in recs] for c in cols}
    full, worst = _loo(list(set_nodes), z_nodes, columns=cols)
    return {
        "missing": [],
        "recs": recs,
        "nodes": set_nodes,
        "z": z_nodes,
        "full": full,
        "worst": worst,
        "deg": deg,
        "res": res,
        "cols": cols,
    }


def _x1_cu5_rows():
    """×1 five-node Cu rows from item 196's stored seven two-column p2 JSONs."""
    recs, missing = _load_nodes(DEFAULT_NODES_MHZ, 2, None, BOTH_COLUMNS)
    if recs is None:
        return None, missing
    idx = [list(DEFAULT_NODES_MHZ).index(m) for m in ASSEMBLY5_NODES_MHZ]
    z5 = {"cu": [_uj(recs[i]["z"]["cu"]) for i in idx]}
    return sparse_rows(ASSEMBLY5_NODES_MHZ, z5, columns=("cu",)), recs


@complex_only
def test_p2_assembly(p2_asm):
    """Once the node set's JSONs exist at this run's (degree, resolution,
    columns) -- seven nodes for both columns, ``ASSEMBLY5_NODES_MHZ`` for Cu
    only: one degree, one cell count, κ None, one resolution, (F) on every
    JSON; one - -> + zero per column, k-share > K_SHARE_FLOOR, Δ_LOO <= BANDS."""
    _need("B")
    if p2_asm["missing"]:
        pytest.skip(f"assembly waits for {p2_asm['missing']}")
    deg, res, cols = p2_asm["deg"], p2_asm["res"], p2_asm["cols"]
    recs, full, worst = p2_asm["recs"], p2_asm["full"], p2_asm["worst"]
    assert all(r["degree"] == deg for r in recs), [r["degree"] for r in recs]
    assert len({r["cells"] for r in recs}) == 1, [r["cells"] for r in recs]
    assert all(r["kappa"] is None for r in recs)
    assert {r.get("resolution") for r in recs} == {res}, [r.get("resolution") for r in recs]
    assert all(tuple(r.get("columns") or BOTH_COLUMNS) == tuple(cols) for r in recs)
    for r in recs:
        _assert_tie(r["cells"], _tie_from_json(r["tie"]), res)
    tag = f"p{deg} h {res if res is not None else 'x1'} {'+'.join(cols)}"
    for c in cols:
        r = full[c]["rows"]
        _say(
            f"[ans7-exact] {tag} assembly ({len(p2_asm['nodes'])} nodes) {c} rows: f0 "
            f"{r.get('f0_hz')!r} Hz, R_in {r.get('r_in_ohm')!r} ohm, Q {r.get('q_loaded')!r}, "
            f"k-share {r.get('k_share')!r}; up-crossings {r.get('crossings_up_hz')}"
        )
        if res is not None:
            continue
        for k in ROWS:
            if k not in r:
                continue
            rag = ITEM157_P2_RAGGED_ROWS[c][k]
            p1 = ANS7_EXACT_STRIP_ROUTE_B_ROWS_P1[c][k]
            lim = 0.01 if k == "f0_hz" else 0.20
            mv = r[k] / rag - 1
            _say(
                f"    {c} {k}: exact p{deg} {r[k]:.10e} vs item 157 ragged in-run p2 {rag:.10e} "
                f"(move {mv:+.4%}; PREDICTED |move| <= {lim:.0%}: "
                f"{'met' if abs(mv) <= lim else 'NOT met'}); vs exact p1 {p1:.10e} "
                f"(order move {r[k] / p1 - 1:+.4%})"
            )
    for m, rec in zip(p2_asm["nodes"], recs):
        ims = ", ".join(f"{c} {rec['im_y_in'][c]:+.9e} S" for c in cols)
        t_cu, t_pec = rec.get("t_cu_s"), rec.get("t_pec_s")
        _say(
            f"[ans7-exact] {tag} node {m:.1f} MHz: Im Y_in {ims}; Cu sweep "
            f"{'n/a' if t_cu is None else f'{t_cu:.1f} s'}, PEC sweep "
            f"{'n/a' if t_pec is None else f'{t_pec:.1f} s'}, -n {rec['ranks']}; summed "
            f"ru_maxrss {rec.get('summed_ru_maxrss_gib', 'n/a')} GiB"
        )
    _b5_report(tag, worst)
    for k in ROWS:
        _say(
            f"[ans7-exact] {tag} Delta_LOO {k} = {worst[k]:.3e} (ASSERTED <= {BANDS[k]:.0e}; "
            f"PREDICTED <= {'1e-9' if res is None else '1e-7'})"
        )
    _control_report(tag, p2_asm["nodes"], p2_asm["z"], full, cols)
    if res is not None and "cu" in cols:
        rung = full["cu"]["rows"]
        x1_5, _x1_recs = _x1_cu5_rows()
        lim = {"f0_hz": 0.005, "r_in_ohm": 0.10, "q_loaded": 0.10}
        for k in ROWS:
            if k not in rung:
                continue
            rec7 = ANS7_EXACT_STRIP_ROUTE_B_ROWS_P2["cu"][k]
            mv7 = rung[k] / rec7 - 1
            msg = (
                f"[ans7-exact] h-move (h {res:g} - x1)/x1 Cu {k}: rung {rung[k]:.10e}; vs "
                f"ANS7_EXACT_STRIP_ROUTE_B_ROWS_P2 (7-node) {rec7:.10e} -> {mv7:+.4%}"
            )
            if x1_5 is not None and k in x1_5["cu"]["rows"]:
                v5 = x1_5["cu"]["rows"][k]
                msg += f"; vs x1 5-node Cu {v5:.10e} -> {rung[k] / v5 - 1:+.4%}"
            msg += (
                f" (PREDICTED |move| <= {lim[k]:.1%}: "
                f"{'met' if abs(mv7) <= lim[k] else 'NOT met'}; a record, never asserted)"
            )
            _say(msg)
        _say(
            f"[ans7-exact] h {res:g} cells {recs[0]['cells']} unknowns {recs[0].get('unknowns')} "
            f"(PREDICTED at h 0.01: cells [135000, 175000], unknowns ~ 6.6 x cells); k-share "
            f"{rung.get('k_share')!r} (PREDICTED >= 0.9999)"
        )
    for c in cols:
        r = full[c]["rows"]
        assert len(r["crossings_up_hz"]) == 1, (c, r["crossings_up_hz"])
        assert r["k_share"] > K_SHARE_FLOOR, (c, r["k_share"])
    for k in ROWS:
        assert worst[k] <= BANDS[k], (k, worst[k])


@complex_only
def test_r2_p2_rows_record(p2_asm):
    """(R2): degree 2, resolution unset, both columns -- the assembly's rows
    reproduce item 196's record within ``NODE_S_RTOL`` ("unset digit-identical")."""
    _need("B")
    if p2_asm["missing"]:
        pytest.skip(f"(R2) waits for the assembly ({p2_asm['missing']})")
    if p2_asm["deg"] != 2 or p2_asm["res"] is not None or tuple(p2_asm["cols"]) != BOTH_COLUMNS:
        pytest.skip("(R2) is asserted at degree 2, resolution unset, both columns")
    fails = []
    for c in BOTH_COLUMNS:
        for k in ROWS:
            got = p2_asm["full"][c]["rows"][k]
            want = ANS7_EXACT_STRIP_ROUTE_B_ROWS_P2[c][k]
            d = abs(got / want - 1.0)
            _say(
                f"[ans7-exact] (R2) {c} {k}: {got!r} vs record {want!r}: rel {d:.3e} "
                f"(ASSERTED <= NODE_S_RTOL {NODE_S_RTOL:.0e}; PREDICTED <= 1e-11: "
                f"{'met' if d <= 1e-11 else 'NOT met'})"
            )
            if d > NODE_S_RTOL:
                fails.append((c, k, got, want, d))
    assert not fails, fails


@complex_only
def test_p2_cu5_subset_anchor():
    """(S5), pure: the Cu-only five-node path at ×1 on item 196's seven stored
    two-column degree-2 JSONs -- five-node Cu rows vs seven-node Cu rows <=
    BANDS, five-node Δ_LOO <= BANDS; B₅ and the 3-node control printed."""
    _need("B")
    recs, missing = _load_nodes(DEFAULT_NODES_MHZ, 2, None, BOTH_COLUMNS)
    if recs is None:
        pytest.skip(f"(S5) needs the seven x1 p2 JSONs ({missing})")
    assert all(r["degree"] == 2 and r["kappa"] is None for r in recs)
    assert len({r["cells"] for r in recs}) == 1, [r["cells"] for r in recs]
    z7 = {"cu": [_uj(r["z"]["cu"]) for r in recs]}
    rows7 = sparse_rows(DEFAULT_NODES_MHZ, z7, columns=("cu",))["cu"]["rows"]
    idx = [list(DEFAULT_NODES_MHZ).index(m) for m in ASSEMBLY5_NODES_MHZ]
    z5 = {"cu": [z7["cu"][i] for i in idx]}
    full5, worst5 = _loo(list(ASSEMBLY5_NODES_MHZ), z5, columns=("cu",))
    rows5 = full5["cu"]["rows"]
    pred_rel = {"f0_hz": 1e-12, "r_in_ohm": 1e-10, "q_loaded": 3e-10}
    pred_loo = {"f0_hz": 1e-9, "r_in_ohm": 4e-8, "q_loaded": 3e-8}
    rel = {}
    for k in ROWS:
        rel[k] = abs(rows5[k] - rows7[k]) / abs(rows7[k])
        _say(
            f"[ans7-exact] (S5) x1 Cu {k}: 5-node {rows5[k]!r} vs 7-node {rows7[k]!r}: rel "
            f"{rel[k]:.3e} (ASSERTED <= {BANDS[k]:.0e}; PREDICTED ~ {pred_rel[k]:.0e})"
        )
    for k in ROWS:
        _say(
            f"[ans7-exact] (S5) x1 Cu 5-node Delta_LOO {k} = {worst5[k]:.3e} (ASSERTED <= "
            f"{BANDS[k]:.0e}; PREDICTED ~ {pred_loo[k]:.0e})"
        )
    _b5_report("(S5) x1 Cu 5-node", worst5)
    _control_report("(S5) x1", ASSEMBLY5_NODES_MHZ, z5, full5, ("cu",))
    _say(
        f"[ans7-exact] (S5) x1 Cu 5-node k-share {rows5.get('k_share')!r}, up-crossings "
        f"{rows5.get('crossings_up_hz')}"
    )
    for k in ROWS:
        assert rel[k] <= BANDS[k], (k, rel[k])
        assert worst5[k] <= BANDS[k], (k, worst5[k])


# -- §9 item 217 (`ANS-7` step 7, chain T21): tune to 64.000 MHz at degree 2 ----
#
# ``FEM_EM_ANS7_EXACT_TUNE`` ∈ {``cmode``, ``inmodel``} (unset ⇒ these tests
# skip; anything else raises); ``ROUTE`` unset in both windows.  W1 (cmode)
# solves the degree-2 ×1 4×4 at 64.0 MHz, applies `TH-17` step 3's C_mode
# construction (``test_th17_birdcage_eigenmodes.py:1093–1100``: ``_modal_z``
# imported, z₁ = ½(z_k[1] + z_k[3]), C_mode = 1/(ω₆₄ Im z₁)), prints C_in (the
# full-network root of Im Y_in(64 MHz)), and reads Route B at C_mode on item
# 196's seven stored nodes.  W2 (inmodel) drives Route A at degree 2 with
# C_mode on ``TUNE_MHZ``.  Anchors (C), (T1), (AB), (T2); controls (N1), (W)
# asserted, (N2) predicted/printed.  Reads item 196's node JSONs, never writes
# them; writes ``logs/ans7_exact/tune_p2.json`` / ``tune_inmodel_p2.json``.

TUNE_ENV = "FEM_EM_ANS7_EXACT_TUNE"
TUNE_TARGET_HZ = 64.0e6
TUNE_MHZ = (63.9, 64.0, 64.1)

# -- §9 item 228 (`ANS-7` step 8, chain T22): T21's construction at 127.74 MHz --
#
# ``FEM_EM_ANS7_EXACT_TUNE_MHZ`` ∈ {unset, ``127.74``}: unset ⇒ T21's objects
# verbatim; ``127.74`` ⇒ the literal constants below (registered before W1; the
# nodes are DEFAULT_NODES_MHZ + 63.74, the same offsets from the target as
# T21's from 64.0).  ``FEM_EM_ANS7_EXACT_TUNE=nodes`` solves a comma subset
# (``FEM_EM_ANS7_EXACT_TUNE_NODES``) of the 127.74 MHz nodes and writes
# ``logs/ans7_exact/t127.74/node_<mhz:.2f>_p2[_cu].json``; it never touches
# item 196's 64 MHz node files or T21's tune JSONs.
TUNE_MHZ_ENV = "FEM_EM_ANS7_EXACT_TUNE_MHZ"
TUNE_NODES_ENV = "FEM_EM_ANS7_EXACT_TUNE_NODES"
TUNE128_TARGET_HZ = 127.74e6
TUNE128_MHZ = (127.64, 127.74, 127.84)
TUNE128_NODES_MHZ = (125.74, 126.04, 126.94, 128.24, 129.54, 130.44, 130.74)
TUNE128_GRID_MHZ = tuple(round(125.74 + 0.1 * k, 6) for k in range(51))
TUNE128_OUT_DIR = OUT_DIR / "t127.74"


def _tune_mhz():
    raw = os.environ.get(TUNE_MHZ_ENV, "").strip()
    if not raw:
        return None
    if raw != "127.74":
        raise ValueError(f"{TUNE_MHZ_ENV}={raw!r}: only unset or '127.74' (not a free frequency)")
    return 127.74


def _active():
    """(target_hz, tune points, nodes, grid, out dir, label) for the knob."""
    if _tune_mhz() is None:
        if _h_tune():
            # §9 item 310: item 214's five h = 0.01 Cu nodes, outputs suffixed.
            return {
                "target_hz": TUNE_TARGET_HZ,
                "tune": TUNE_MHZ,
                "nodes": H_NODES_MHZ,
                "grid": GRID_MHZ,
                "out": OUT_DIR,
                "lab": "T28",
                "tmhz": "64.0",
                "res": H_RES,
                "sfx": H_SFX,
            }
        return {
            "target_hz": TUNE_TARGET_HZ,
            "tune": TUNE_MHZ,
            "nodes": DEFAULT_NODES_MHZ,
            "grid": GRID_MHZ,
            "out": OUT_DIR,
            "lab": "T21",
            "tmhz": "64.0",
        }
    if _h_tune():
        # §9 item 365 (T34): T22′'s seven-node set on the h = 0.01 rung.
        return {
            "target_hz": TUNE128_TARGET_HZ,
            "tune": TUNE128_MHZ,
            "nodes": TUNE128_NODES_MHZ,
            "grid": TUNE128_GRID_MHZ,
            "out": TUNE128_OUT_DIR,
            "lab": "T34",
            "tmhz": "127.74",
            "res": H_RES,
            "sfx": H_SFX,
        }
    return {
        "target_hz": TUNE128_TARGET_HZ,
        "tune": TUNE128_MHZ,
        "nodes": TUNE128_NODES_MHZ,
        "grid": TUNE128_GRID_MHZ,
        "out": TUNE128_OUT_DIR,
        "lab": "T22",
        "tmhz": "127.74",
    }


def _node128_json(mhz, cols, res=None):
    if res is not None:
        # §9 item 365 (T34): the h = 0.01 Cu node files.
        return TUNE128_OUT_DIR / f"node_{mhz:.2f}_p2{H_SFX}_cu.json"
    suffix = "_cu" if tuple(cols) == ("cu",) else ""
    return TUNE128_OUT_DIR / f"node_{mhz:.2f}_p2{suffix}.json"


def _tune_nodes_env():
    raw = os.environ.get(TUNE_NODES_ENV, "").strip()
    nodes = tuple(round(float(v), 6) for v in raw.split(",") if v.strip())
    bad = [m for m in nodes if m not in TUNE128_NODES_MHZ]
    if bad or not nodes:
        raise ValueError(f"{TUNE_NODES_ENV}={raw!r}: nodes must be a subset of {TUNE128_NODES_MHZ}")
    return nodes


def _tune():
    raw = os.environ.get(TUNE_ENV, "").strip()
    if not raw:
        return None
    if raw not in ("cmode", "inmodel", "nodes"):
        raise ValueError(f"{TUNE_ENV}={raw!r}: only 'cmode', 'inmodel' or 'nodes'")
    if raw == "nodes" and _tune_mhz() != 127.74:
        raise ValueError(
            f"{TUNE_ENV}=nodes needs {TUNE_MHZ_ENV}=127.74 (item 196's 64 MHz "
            f"node files are never re-solved here)"
        )
    return raw


# -- §9 item 262 (`ANS-7` step 7′, chain T21′): drive at C_in, read the in-model zero --
#
# ``FEM_EM_ANS7_EXACT_TUNE_C`` ∈ {unset, ``cin``}: unset ⇒ T21's C_mode paths
# byte-identical; ``cin`` (needs ``TUNE=inmodel``, ``TUNE_MHZ`` unset) swaps the
# Route A / Route B capacitor to C_in read from ``tune_p2.json`` (item 217 W1,
# never recomputed) and writes ``tune_inmodel_cin_p2.json``.
TUNE_C_ENV = "FEM_EM_ANS7_EXACT_TUNE_C"
# The grid's own frequency resolution, pre-registered before W1 (literal, never
# recomputed from this run): the max of T21's (AB) node residuals at C_mode at
# 63.9 / 64.0 / 64.1 MHz, 1.237e-5 / 1.258e-5 / 1.287e-5
# (``20261004T140738Z_ANS-7-item217-W2-inmodel.log:2046–2048``).  It exceeds
# |r| (−2.531e-6 at 64 MHz, ``…140738Z…:2050``) ×5.1.
T21P_ZERO_RESOLUTION = 1.287e-5
# §9 item 283 (`ANS-7` step 8′, chain T22′): the 127.74 MHz band, pre-registered
# before W1 (literal, never recomputed from this run): the max of T22's (AB) node
# residuals at C_mode at 127.64 / 127.74 / 127.84 MHz, 1.589e-5 / 1.601e-5 /
# 1.616e-5 (``20261004T204206Z_ANS-7-item228-W4-inmodel.log:2054–2056``).  It
# exceeds |r| (−2.263e-6, ``…204206Z…:2058``) ×7.1.  Chosen by the knob, never
# by a reading.
T22P_ZERO_RESOLUTION = 1.616e-5


def _zero_resolution():
    """The active (T2′)/(T2″) band: T21P at TUNE_MHZ unset, T22P at 127.74."""
    return T21P_ZERO_RESOLUTION if _tune_mhz() is None else T22P_ZERO_RESOLUTION


def _tune_c():
    raw = os.environ.get(TUNE_C_ENV, "").strip()
    if not raw:
        return None
    if raw != "cin":
        raise ValueError(f"{TUNE_C_ENV}={raw!r}: only unset or 'cin'")
    # §9 item 283: TUNE_MHZ ∈ {unset, 127.74} (``_tune_mhz`` raises on anything else).
    if _tune() != "inmodel" or _tune_mhz() not in (None, 127.74):
        raise ValueError(
            f"{TUNE_C_ENV}=cin needs {TUNE_ENV}=inmodel and {TUNE_MHZ_ENV} unset " f"or '127.74'"
        )
    return raw


# -- §9 item 310 (`ANS-7` step 12, chain T28): C_in on step 6's h = 0.01 rung --
#
# ``FEM_EM_ANS7_EXACT_RESOLUTION=0.01`` with ``TUNE`` ∈ {``cmode``, ``inmodel``}
# is legal only at 64 MHz (``TUNE_MHZ`` unset) and Cu only (``COLUMNS=cu``); any
# other tune combination with a resolution raises.  ``cmode`` solves the rung's
# Cu 4×4 at 64 MHz (``pec_column=False``, ``field_solve=False``: outside item
# 157's trap, which is ``pec_column=False`` *with* ``field_solve=True``) and
# writes ``tune_p2_h0.01.json`` (C_mode(h), C_in(h), z64).  ``inmodel`` with
# ``TUNE_C=cin`` needs ``FEM_EM_ANS7_EXACT_TUNE_POINT`` ∈ {63.9, 64.0, 64.1,
# ``assemble``}: one Route A drive per window (rank 0 writes
# ``tune_inmodel_cin_p2_h0.01_<mhz>.json``; the (W) ×1 degree-1 wiring drive runs
# in the 64.0 window only), or the pure assembly (``tune_inmodel_cin_p2_h0.01.json``).
# Unset resolution / unset point: every existing path unchanged.
H_RES = 0.01
H_SFX = "_h0.01"
H_NODES_MHZ = ASSEMBLY5_NODES_MHZ
TUNE_POINT_ENV = "FEM_EM_ANS7_EXACT_TUNE_POINT"
# Step 6's rung rows at C_F (item 214, PROJECT_PLAN §7 ANS-7 row): printed comparands.
STEP6_RUNG_F0_MOVE = 4.11e-4
STEP6_RUNG_ROWS = {"f0_hz": 62.726847e6, "r_in_ohm": 18313.602, "q_loaded": 970.439}
# §9 item 365 (T34): the 64 MHz C_in h-move, a printed record literal
# (``20261008T111045Z_ANS-7-item310-W1-cin-h001.log:2234``), never asserted.
T28_C_IN_H_REL_64 = 8.346526e-4


def _h_tune():
    """True on the T28 path (resolution set with a 64 MHz Cu tune); raises on any
    other tune combination with a resolution; False when either is unset."""
    res = _resolution()
    if res is None:
        return False
    t = _tune()
    if t is None:
        return False
    # §9 item 365 (T34): at 127.74 MHz the h = 0.01 rung also takes ``nodes``
    # (T22′'s seven-node set); the 64 MHz set stays {cmode, inmodel}.
    legal = ("cmode", "inmodel") if _tune_mhz() is None else ("nodes", "cmode", "inmodel")
    if t not in legal or res != H_RES or _cols() != ("cu",):
        raise ValueError(
            f"{RESOLUTION_ENV}={res!r} with {TUNE_ENV}={t!r}: only {H_RES:g} with "
            f"cmode/inmodel ({TUNE_MHZ_ENV} unset) or nodes/cmode/inmodel "
            f"({TUNE_MHZ_ENV}=127.74), and {COLUMNS_ENV}=cu"
        )
    return True


def _tune_point():
    raw = os.environ.get(TUNE_POINT_ENV, "").strip()
    h_in = _h_tune() and _tune() == "inmodel"
    if not raw:
        # §9 item 395 (T37): under the ESR knob the rung runs in ``tune_esr``, unset point.
        if h_in and _esr_ohm() is None:
            raise ValueError(
                f"{TUNE_POINT_ENV} unset on the h = {H_RES:g} rung: one point per "
                f"window (63.9 / 64.0 / 64.1) or 'assemble'"
            )
        return None
    if not h_in or _tune_c() != "cin":
        raise ValueError(
            f"{TUNE_POINT_ENV}={raw!r} needs {RESOLUTION_ENV}={H_RES:g}, "
            f"{TUNE_ENV}=inmodel, {TUNE_C_ENV}=cin"
        )
    if raw == "assemble":
        return raw
    try:
        v = round(float(raw), 6)
    except ValueError:
        v = None
    if v not in _active()["tune"]:
        raise ValueError(f"{TUNE_POINT_ENV}={raw!r}: only {_active()['tune']} or 'assemble'")
    return v


def _need_tune(mode):
    t = _tune()
    if t != mode:
        pytest.skip(f"{TUNE_ENV}={t!r}: this test runs on {mode}")
    if _route() is not None:
        raise ValueError(f"{TUNE_ENV} runs with {ROUTE_ENV} unset")


def _tune_json():
    act = _active()
    return act["out"] / f"tune_p{_degree()}{act.get('sfx', '')}.json"


def _tune_inmodel_json():
    return _active()["out"] / f"tune_inmodel_p{_degree()}.json"


def _met(ok):
    return "met" if ok else "NOT met"


def _stored_x1_p2():
    """Item 196's seven ×1 two-column degree-2 node JSONs, checked (fail if absent)."""
    recs, missing = _load_nodes(DEFAULT_NODES_MHZ, 2, None, BOTH_COLUMNS)
    if recs is None:
        pytest.fail(
            f"item 196's x1 p2 node JSONs absent: {missing} (stop; re-running item 196 "
            f"costs ~21 min)"
        )
    bad = [
        (r["f_hz"], r["degree"], r["cells"], r["kappa"], r.get("resolution"))
        for r in recs
        if r["degree"] != 2 or r["kappa"] is not None or r.get("resolution") is not None
    ]
    if bad or len({r["cells"] for r in recs}) != 1:
        pytest.fail(
            f"item 196's x1 p2 node JSONs inconsistent: {bad} "
            f"cells {[r['cells'] for r in recs]}"
        )
    for r in recs:
        _assert_tie(r["cells"], _tie_from_json(r["tie"]))
    return recs


def _stored_active():
    """(recs, columns) of the active node set: item 196's at unset (both
    columns), the ``t127.74/`` files at ``_cols()`` at 127.74 (checked); item
    214's five ``_h0.01_cu`` files on the T28 path (§9 item 310, checked); the
    seven ``t127.74/node_*_p2_h0.01_cu.json`` files on the T34 path (§9 item 365)."""
    if _h_tune() and _tune_mhz() == 127.74:
        paths = [_node128_json(m, ("cu",), H_RES) for m in TUNE128_NODES_MHZ]
        missing = [p.name for p in paths if not p.exists()]
        if missing:
            pytest.fail(
                f"T34's h = {H_RES:g} 127.74 MHz Cu node JSONs absent: {missing} "
                f"(run TUNE=nodes with {RESOLUTION_ENV}={H_RES:g} first)"
            )
        recs = [json.loads(p.read_text()) for p in paths]
        bad = [
            (
                r["f_hz"],
                r["degree"],
                r["cells"],
                r["kappa"],
                r.get("resolution"),
                r.get("field_solve"),
                r.get("columns"),
                r.get("target_mhz"),
            )
            for r in recs
            if r["degree"] != 2
            or r["kappa"] is not None
            or r.get("resolution") != H_RES
            or r.get("field_solve") is not False
            or r.get("columns") != ["cu"]
            or r.get("target_mhz") != 127.74
        ]
        if bad or len({r["cells"] for r in recs}) != 1:
            pytest.fail(
                f"T34's node JSONs inconsistent: {bad} " f"cells {[r['cells'] for r in recs]}"
            )
        for m, r in zip(TUNE128_NODES_MHZ, recs):
            assert abs(r["f_hz"] - m * 1e6) < 1.0, (m, r["f_hz"])
            _assert_tie(r["cells"], _tie_from_json(r["tie"]), res=H_RES)
        return recs, ("cu",)
    if _h_tune():
        recs, missing = _load_nodes(H_NODES_MHZ, 2, H_RES, ("cu",))
        if recs is None:
            pytest.fail(f"item 214's h = {H_RES:g} Cu node JSONs absent: {missing}")
        bad = [
            (
                r["f_hz"],
                r["degree"],
                r["cells"],
                r["kappa"],
                r.get("resolution"),
                r.get("field_solve"),
                r.get("columns"),
            )
            for r in recs
            if r["degree"] != 2
            or r["kappa"] is not None
            or r.get("resolution") != H_RES
            or r.get("field_solve") is not False
            or r.get("columns") != ["cu"]
        ]
        if bad or len({r["cells"] for r in recs}) != 1:
            pytest.fail(
                f"item 214's node JSONs inconsistent: {bad} " f"cells {[r['cells'] for r in recs]}"
            )
        for m, r in zip(H_NODES_MHZ, recs):
            assert abs(r["f_hz"] - m * 1e6) < 1.0, (m, r["f_hz"])
            _assert_tie(r["cells"], _tie_from_json(r["tie"]), res=H_RES)
        return recs, ("cu",)
    if _tune_mhz() is None:
        return _stored_x1_p2(), BOTH_COLUMNS
    cols = _cols()
    paths = [_node128_json(m, cols) for m in TUNE128_NODES_MHZ]
    missing = [p.name for p in paths if not p.exists()]
    if missing:
        pytest.fail(f"item 228's 127.74 MHz node JSONs absent: {missing} (run TUNE=nodes first)")
    recs = [json.loads(p.read_text()) for p in paths]
    bad = [
        (r["f_hz"], r["degree"], r["cells"], r["kappa"], r.get("resolution"), r.get("target_mhz"))
        for r in recs
        if r["degree"] != 2
        or r["kappa"] is not None
        or r.get("resolution") is not None
        or r.get("target_mhz") != 127.74
    ]
    if bad or len({r["cells"] for r in recs}) != 1:
        pytest.fail(f"item 228's node JSONs inconsistent: {bad} cells {[r['cells'] for r in recs]}")
    for m, r in zip(TUNE128_NODES_MHZ, recs):
        assert abs(r["f_hz"] - m * 1e6) < 1.0, (m, r["f_hz"])
        _assert_tie(r["cells"], _tie_from_json(r["tie"]))
    return recs, cols


@pytest.fixture(scope="module")
def tune_nodes():
    from tests.validation import test_th14_birdcage_copper as th14
    from tests.validation.test_ans7_sparse_route_b import _column_power

    _need_tune("nodes")
    nodes = _tune_nodes_env()
    comm = MPI.COMM_WORLD
    deg = _degree()
    if deg != 2:
        raise ValueError(f"{TUNE_ENV}=nodes runs at degree 2 only, got {deg}")
    cols = _cols()
    pec_column = "pec" in cols
    # §9 item 365 (T34): the h = 0.01 rung at 127.74 MHz; ``None`` is the default.
    res = H_RES if _h_tune() else None
    rlab = "x1" if res is None else f"h {res:g}"
    tie = _fixture_tie(res)
    t0 = time.perf_counter()
    lad = th14._build_ladder(
        deg,
        sigmas=(th14.COPPER_SIGMA,),
        frequencies_hz=[m * 1e6 for m in nodes],
        width_correction_kappa=None,
        measure_kappa=False,
        pec_column=pec_column,
        field_solve=False,
        leg_sheet_strip_fraction=GATED_WIDTH_FRACTION,
        resolution=res,
    )
    comm.Barrier()
    t_solve = time.perf_counter() - t0
    rss_gib = (
        comm.allreduce(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, op=MPI.SUM) / 1024**2
    )
    recs = sorted(lad["freqs"].values(), key=lambda r: r["f_hz"])
    assert len(recs) == len(nodes)
    cells = int(lad["cells"])
    unknowns = int(lad["unknowns"])
    out = {"nodes": nodes, "cells": cells, "tie": tie, "per": {}, "res": res}
    if comm.rank == 0:
        TUNE128_OUT_DIR.mkdir(parents=True, exist_ok=True)
    for mhz, rec in zip(nodes, recs):
        assert abs(rec["f_hz"] - mhz * 1e6) < 1.0
        per = {"s": {}, "z": {}, "recip": {}, "sigma_max": {}, "im_y": {}, "colpow": {}}
        pairs = [("cu", rec["sigma"][th14.COPPER_SIGMA])]
        if pec_column:
            pairs.append(("pec", rec["pec"]))
        per["t_cu_s"] = float(rec["sigma"][th14.COPPER_SIGMA]["t"])
        per["t_pec_s"] = float(rec["t_pec"]) if pec_column else None
        _say(
            f"[ans7-exact] T22 p{deg} node {mhz:.2f} MHz: Cu sweep {per['t_cu_s']:.1f} s"
            + (f", PEC sweep {per['t_pec_s']:.1f} s" if pec_column else "")
            + f" ({rlab}, columns {cols}, -n {comm.size})"
        )
        for c, sw in pairs:
            s4 = np.asarray(sw["s"], dtype=np.complex128)
            per["s"][c], per["z"][c] = s4, z_from_s(s4)
            per["recip"][c] = float(sw["reciprocity"])
            per["sigma_max"][c] = float(np.max(sw["sigma"]))
            per["im_y"][c], z_in = _im_y_in(per["z"][c], mhz)
            per["colpow"][c] = _column_power(s4)
            _say(
                f"[ans7-exact] T22 p{deg} node {mhz:.2f} MHz {c}: Im Y_in at C_F "
                f"{per['im_y'][c]:+.9e} S (PREDICTED > 0: {_met(per['im_y'][c] > 0)}) (Z_in "
                f"{z_in.real:.6e}{z_in.imag:+.6e}j ohm); recip {per['recip'][c]:.3e}; sigma_max "
                f"{per['sigma_max'][c]:.12f}; column power {per['colpow'][c]:.9e}"
            )
        out["per"][mhz] = per
        if comm.rank == 0:
            _node128_json(mhz, cols, res).write_text(
                json.dumps(
                    {
                        "f_hz": rec["f_hz"],
                        "degree": deg,
                        "ranks": comm.size,
                        "cells": cells,
                        "tie": _tie_json(tie),
                        "kappa": None,
                        "field_solve": False,
                        "s": {c: _cj(per["s"][c]) for c in per["s"]},
                        "z": {c: _cj(per["z"][c]) for c in per["z"]},
                        "recip": per["recip"],
                        "sigma_max": per["sigma_max"],
                        "im_y_in": per["im_y"],
                        "power": None,
                        "identity_residual": None,
                        "resolution": res,
                        "columns": list(cols),
                        "unknowns": unknowns,
                        "t_cu_s": per["t_cu_s"],
                        "t_pec_s": per["t_pec_s"],
                        "summed_ru_maxrss_gib": rss_gib,
                        "target_mhz": 127.74,
                    }
                )
            )
    comm.Barrier()
    _say(
        f"[ans7-exact] T22 p{deg} nodes {nodes} solved in {t_solve:.1f} s, cells {cells}, "
        f"unknowns {unknowns}, columns {cols}, -n {comm.size}; summed ru_maxrss "
        f"{rss_gib:.2f} GiB; JSON in {TUNE128_OUT_DIR}"
    )
    del lad, recs
    return out


@complex_only
def test_t22_node_anchors(tune_nodes):
    """Per 127.74 MHz node: reciprocity, passivity, and the (F) tie on the ×1 record."""
    _need_tune("nodes")
    _assert_tie(tune_nodes["cells"], tune_nodes["tie"], res=tune_nodes["res"])
    for mhz, per in tune_nodes["per"].items():
        for c in per["recip"]:
            assert per["recip"][c] <= RECIPROCITY_BAND, (mhz, c, per["recip"][c])
            assert per["sigma_max"][c] <= 1.0 + PASSIVITY_SIGMA_TOLERANCE, (
                mhz,
                c,
                per["sigma_max"][c],
            )


# Printed predictions, keyed by target (64 MHz: T21's text; 127.74: item 228's).
_PRED = {
    "64.0": {
        "cmode": "[63.4, 63.8] pF",
        "cmode_lo": 63.4e-12,
        "cmode_hi": 63.8e-12,
        "pec": "~ +4e-4",
        "off": 1e-6,
        "cin": "|.| in [1e-7, 5e-5]",
        "cin_lo": 1e-7,
        "cin_hi": 5e-5,
        "f0b": "[1e-7, 5e-5]",
        "f0b_lo": 1e-7,
        "f0b_hi": 5e-5,
        "ab_gap": "~ -3e-6",
    },
    "127.74": {
        "cmode": "[14.8, 15.6] pF",
        "cmode_lo": 14.8e-12,
        "cmode_hi": 15.6e-12,
        "pec": "in [+1e-4, +1e-3]",
        "off": None,
        "cin": "negative, |.| in [3e-6, 1e-4]",
        "cin_lo": 3e-6,
        "cin_hi": 1e-4,
        "f0b": "[3e-6, 5e-5]",
        "f0b_lo": 3e-6,
        "f0b_hi": 5e-5,
        "ab_gap": "in [-1.5e-5, -1e-6]",
    },
}


def _c_mode(z4, omega):
    """`TH-17` step 3's construction (``test_th17_birdcage_eigenmodes.py:1095–1100``)."""
    from tests.validation.test_th17_birdcage_eigenmodes import _modal_z

    zk, off = _modal_z(np.asarray(z4, dtype=np.complex128))
    z1 = complex(0.5 * (zk[1] + zk[3]))
    split = float(abs(zk[1] - zk[3]) / abs(z1))
    c = 1.0 / (omega * z1.imag) if z1.imag > 0.0 else float("nan")
    return {"zk": zk, "off": off, "z1": z1, "split": split, "c_mode": c}


def _c_in(z4, f_hz, c0):
    """Root of g(C) = Im[1/Z_in(C)] at ``f_hz`` on the full 4×4, secant from ``c0``
    on the real scalar g, to |g| <= 1e-15 S in <= 20 iterations."""

    def g(c):
        return float((1.0 / complex(route_b_z(z4, f_hz, c)[0])).imag)

    ca, cb = float(c0), float(c0) * (1.0 + 1e-6)
    ga, gb = g(ca), g(cb)
    it = 0
    while abs(gb) > 1e-15 and it < 20:
        if gb == ga:
            break
        ca, cb = cb, cb - gb * (cb - ca) / (gb - ga)
        ga, gb = gb, g(cb)
        it += 1
    return cb, gb, it


@pytest.fixture(scope="module")
def tune_cmode():
    from tests.validation import test_th14_birdcage_copper as th14

    _need_tune("cmode")
    comm = MPI.COMM_WORLD
    deg = _degree()
    act = _active()
    tgt, lab, tm, pr = act["target_hz"], act["lab"], act["tmhz"], _PRED[act["tmhz"]]
    recs, cols = _stored_active()
    z_nodes = {c: [_uj(r["z"][c]) for r in recs] for c in cols}
    # §9 item 310: on the h rung Cu only, no PEC sweep, field solve off (item 214's
    # flags; outside item 157's pec_column=False + field_solve=True trap).
    h = _h_tune()
    hres = H_RES if h else None
    scols = ("cu",) if h else BOTH_COLUMNS
    tie = _fixture_tie(hres)
    t0 = time.perf_counter()
    lad = th14._build_ladder(
        deg,
        sigmas=(th14.COPPER_SIGMA,),
        resolution=hres,
        frequencies_hz=[tgt],
        width_correction_kappa=None,
        measure_kappa=False,
        pec_column=not h,
        field_solve=False,
        leg_sheet_strip_fraction=GATED_WIDTH_FRACTION,
    )
    comm.Barrier()
    t_solve = time.perf_counter() - t0
    rec = list(lad["freqs"].values())[0]
    assert abs(rec["f_hz"] - tgt) < 1.0
    cells = int(lad["cells"])
    omega = 2.0 * np.pi * tgt
    out = {
        "cells": cells,
        "tie": tie,
        "deg": deg,
        "recip": {},
        "sigma_max": {},
        "z64": {},
        "cm": {},
        "stored_cells": recs[0]["cells"],
        "cols": cols,
        "res": hres,
        "scols": scols,
    }
    xn = np.asarray(act["nodes"], dtype=float) - 64.5
    pairs = [("cu", rec["sigma"][th14.COPPER_SIGMA])]
    if not h:
        pairs.append(("pec", rec["pec"]))
    for c, sw in pairs:
        s4 = np.asarray(sw["s"], dtype=np.complex128)
        z64 = z_from_s(s4)
        out["z64"][c] = z64
        out["recip"][c] = float(sw["reciprocity"])
        out["sigma_max"][c] = float(np.max(sw["sigma"]))
        cm = _c_mode(z64, omega)
        out["cm"][c] = cm
        if pr["off"] is not None:
            off_txt = f"(PREDICTED <= 1e-6: {_met(cm['off'] <= 1e-6)})"
        else:
            off_txt = f"(PREDICTED [5e-6, 3e-5]: {_met(5e-6 <= cm['off'] <= 3e-5)})"
        _say(
            f"[ans7-exact] {lab} {tm} MHz p{deg} {c}: recip {out['recip'][c]:.3e}, sigma_max "
            f"{out['sigma_max'][c]:.12f}; off-diag {cm['off']:.3e} {off_txt}; z_k = "
            + ", ".join(f"{complex(z):.9e}" for z in cm["zk"])
            + f"; z1 = {cm['z1']:.12e}, |z1-z3|/|z1| {cm['split']:.3e}; C_mode {cm['c_mode']!r} F"
        )
        if c not in z_nodes:
            continue
        zi = bary_interp(xn, z_nodes[c], tgt / 1e6 - 64.5)
        d = np.abs(zi - z64)
        ent = float(np.max(d / np.abs(z64)))
        glob_rel = float(np.max(d) / np.max(np.abs(z64)))
        _say(
            f"[ans7-exact] {lab} {c} interpolated-vs-solved Z64: entrywise max rel "
            f"{ent:.3e}, max|dZ|/max|Z| {glob_rel:.3e} (PREDICTED <= 1e-8: "
            f"{_met(max(ent, glob_rel) <= 1e-8)})"
        )
    c_cu = out["cm"]["cu"]["c_mode"]
    out["c_mode"] = c_cu
    if h:
        _say(
            f"[ans7-exact] {lab} C_mode,Cu(h = {H_RES:g}) = {c_cu!r} F = {c_cu * 1e12:.6f} pF "
            f"(printed; x1 PREDICTED band {pr['cmode']}: "
            f"{_met(pr['cmode_lo'] <= c_cu <= pr['cmode_hi'])})"
        )
    c_pec = out["cm"]["pec"]["c_mode"] if not h else float("nan")
    _say(
        f"[ans7-exact] {lab} C_mode,Cu = {c_cu!r} F = {c_cu * 1e12:.6f} pF (PREDICTED "
        f"{pr['cmode']}: {_met(pr['cmode_lo'] <= c_cu <= pr['cmode_hi'])}); C_mode,PEC = "
        f"{c_pec!r} F; (PEC - Cu)/Cu = {c_pec / c_cu - 1:+.3e} (PREDICTED {pr['pec']}); C_F "
        f"{C_F!r}, C_F*(62.701082/{'64' if pr['off'] is not None else tm})^2 ={C_F * (62.701082 / (tgt / 1e6)) ** 2:.6e}"
    )
    if pr["off"] is None:
        _say(
            f"[ans7-exact] {lab} Im z1,Cu = {out['cm']['cu']['z1'].imag:.6f} ohm (PREDICTED "
            f"[79.8, 84.2])"
        )
    c_in, g_in, it = _c_in(out["z64"]["cu"], tgt, c_cu)
    out["c_in"], out["g_in"], out["c_in_iter"] = c_in, g_in, it
    off_ci = c_cu / c_in - 1.0
    ci_ok = pr["cin_lo"] <= abs(off_ci) <= pr["cin_hi"] and (pr["off"] is not None or off_ci < 0)
    _say(
        f"[ans7-exact] {lab} C_in (root of Im Y_in({tm} MHz) on the solved Cu 4x4) = {c_in!r} F "
        f"(|g| {abs(g_in):.3e} S after {it} secant iterations; target <= 1e-15 S in <= 20: "
        f"{_met(abs(g_in) <= 1e-15 and it <= 20)}); (C_mode - C_in)/C_in = {off_ci:+.6e} "
        f"(PREDICTED {pr['cin']}: {_met(ci_ok)})"
    )
    if h:
        # §9 item 310: the first h-move of C_in, against the x1 tune_p2.json (read only).
        # §9 item 365: the x1 comparand is the active target's (OUT_DIR at 64 MHz).
        tx1 = json.loads((act["out"] / f"tune_p{deg}.json").read_text())
        c_x1 = float(tx1["c_in"])
        c_h_rel = c_in / c_x1 - 1.0
        # (N-h): the secant stops at |g| <= 1e-15 S; dg/dC = omega at the root, so its
        # C tolerance is 1e-15/omega, relative 1e-15/(omega C_in).
        sec_rel = 1e-15 / (omega * c_in)
        out["c_x1"], out["c_h_rel"], out["sec_rel"] = c_x1, c_h_rel, sec_rel
        if _tune_mhz() is None:
            _say(
                f"[ans7-exact] T28 C_in(h = {H_RES:g}) = {c_in!r} F; C_in(x1) = {c_x1!r} F (tune_p{deg}"
                f".json); C_in(h)/C_in(x1) - 1 = {c_h_rel:+.6e} (PREDICTED [-1.5e-3, -0.3e-3], ~ -8.2e-4 "
                f"from f ~ C^-1/2: {_met(-1.5e-3 <= c_h_rel <= -0.3e-3)}); beside step 6's f0 move "
                f"{STEP6_RUNG_F0_MOVE:+.4%} at fixed C (-2 x = {-2.0 * STEP6_RUNG_F0_MOVE:+.3e})"
            )
        else:
            # §9 item 365b's prediction (T34): positive, as at 64 MHz.
            _say(
                f"[ans7-exact] {lab} C_in(h = {H_RES:g}) = {c_in!r} F; C_in(x1) = {c_x1!r} F "
                f"(t127.74/tune_p{deg}.json); C_in(h)/C_in(x1) - 1 at {tm} MHz = {c_h_rel:+.6e} "
                f"(PREDICTED [+2e-4, +2e-3], sign positive: "
                f"{_met(2e-4 <= c_h_rel <= 2e-3)}); beside 64 MHz's record "
                f"{T28_C_IN_H_REL_64:+.6e}"
            )
        _say(
            f"[ans7-exact] {lab} (N-h) |C_in(h)/C_in(x1) - 1| = {abs(c_h_rel):.3e} vs 10 x secant "
            f"tolerance 10 x {sec_rel:.3e} = {10 * sec_rel:.3e} (PREDICTED larger: "
            f"{_met(abs(c_h_rel) >= 10 * sec_rel)}; printed, never asserted)"
        )
    full, worst = _loo(list(act["nodes"]), z_nodes, columns=cols, c_f=c_cu, grid_mhz=act["grid"])
    out["full"], out["worst"] = full, worst
    rows196 = ANS7_EXACT_STRIP_ROUTE_B_ROWS_P2["cu"]
    for c in cols:
        r = full[c]["rows"]
        _say(
            f"[ans7-exact] {lab} (T1) Route B at C_mode,Cu {c} rows: f0 {r.get('f0_hz')!r} Hz, "
            f"R_in {r.get('r_in_ohm')!r} ohm, Q {r.get('q_loaded')!r}, k-share "
            f"{r.get('k_share')!r}; up-crossings {r.get('crossings_up_hz')}"
        )
    rc = full["cu"]["rows"]
    if "f0_hz" in rc and pr["off"] is not None:
        df = abs(rc["f0_hz"] / tgt - 1.0)
        dr = rc["r_in_ohm"] / rows196["r_in_ohm"] - 1.0
        dq = rc["q_loaded"] / rows196["q_loaded"] - 1.0
        _say(
            f"[ans7-exact] T21 |f0,B/64 MHz - 1| (Cu) = {df:.6e} (printed; PREDICTED [1e-7, "
            f"5e-5]: {_met(1e-7 <= df <= 5e-5)}); R_in vs item 196 {dr:+.4%}, Q vs item 196 "
            f"{dq:+.4%} (PREDICTED within +-10%: {_met(abs(dr) <= 0.1 and abs(dq) <= 0.1)}); "
            f"k-share {rc.get('k_share')!r} (PREDICTED >= 0.9999: "
            f"{_met(rc.get('k_share', 0.0) >= 0.9999)})"
        )
    elif "f0_hz" in rc:
        sf0 = rc["f0_hz"] / tgt - 1.0
        q = rc["q_loaded"]
        _say(
            f"[ans7-exact] {lab} f0,B/{tm} MHz - 1 (Cu) = {sf0:+.6e} (printed; PREDICTED |.| "
            f"{pr['f0b']}: {_met(pr['f0b_lo'] <= abs(sf0) <= pr['f0b_hi'])}); R_in "
            f"{rc['r_in_ohm']!r} ohm, Q {q!r} (PREDICTED [400, 2500]: "
            f"{_met(400 <= q <= 2500)}); k-share {rc.get('k_share')!r} (PREDICTED >= 0.9999: "
            f"{_met(rc.get('k_share', 0.0) >= 0.9999)})"
        )
        # Sensitivity: the band is resolvable.
        sens = sparse_rows(
            act["nodes"], z_nodes, grid_mhz=act["grid"], columns=("cu",), c_f=c_cu * (1.0 + 1e-5)
        )["cu"]["rows"]
        if "f0_hz" in sens:
            _say(
                f"[ans7-exact] {lab} sensitivity |f0,B(C_mode*(1+1e-5))/f0,B(C_mode) - 1| = "
                f"{abs(sens['f0_hz'] / rc['f0_hz'] - 1.0):.6e} (PREDICTED ~ 5e-6; printed)"
            )
    for k in ROWS:
        _say(
            f"[ans7-exact] {lab} (T1) Delta_LOO {k} at C_mode = {worst[k]:.3e} (ASSERTED <= "
            f"{BANDS[k]:.0e}"
            + ("" if pr["off"] is not None else f"; PREDICTED <= 1e-9: {_met(worst[k] <= 1e-9)}")
            + ")"
        )
    # (N0), predicted / printed, 127.74 only: at C_F no − → + zero on the new grid.
    if pr["off"] is None:
        n0 = sparse_rows(act["nodes"], z_nodes, grid_mhz=act["grid"], columns=("cu",), c_f=C_F)
        n0x = n0["cu"]["rows"].get("crossings_up_hz", [])
        out["n0"] = n0x
        _say(
            f"[ans7-exact] {lab} (N0) at C_F on the {tm} MHz Cu grid: up-crossings {n0x} "
            f"(PREDICTED none: {_met(len(n0x) == 0)}; printed, never asserted)"
        )
    # (N1): c_f wiring at C_F on item 196's stored 64 MHz nodes (either target).
    recs196 = recs if (_tune_mhz() is None and not h) else _stored_x1_p2()
    z196 = {c: [_uj(r["z"][c]) for r in recs196] for c in BOTH_COLUMNS}
    n1 = sparse_rows(DEFAULT_NODES_MHZ, z196, columns=BOTH_COLUMNS, c_f=C_F)
    out["n1"] = n1
    for c in BOTH_COLUMNS:
        for k in ROWS:
            got, want = n1[c]["rows"][k], ANS7_EXACT_STRIP_ROUTE_B_ROWS_P2[c][k]
            _say(
                f"[ans7-exact] {lab} (N1) c_f=C_F {c} {k}: {got!r} vs record {want!r}: rel "
                f"{abs(got / want - 1.0):.3e} (ASSERTED <= NODE_S_RTOL {NODE_S_RTOL:.0e})"
            )
    nf0 = abs(n1["cu"]["rows"]["f0_hz"] / tgt - 1.0)
    out["n1_f0_off"] = nf0
    _say(
        f"[ans7-exact] {lab} (N1) |f0(C_F)/{tm} MHz - 1| = {nf0:.6e} = "
        f"{nf0 / BANDS['f0_hz']:.3e} x BANDS (ASSERTED > N_FACTOR {N_FACTOR:g} x "
        f"{BANDS['f0_hz']:.0e}; backing 2.03e-2)"
    )
    _say(
        f"[ans7-exact] {lab} W1 {tm} MHz node solved in {t_solve:.1f} s, cells {cells}, "
        f"-n {comm.size}"
    )
    if comm.rank == 0:
        act["out"].mkdir(parents=True, exist_ok=True)
        payload = {
            "degree": deg,
            "ranks": comm.size,
            "cells": cells,
            "f_hz": tgt,
            "c_mode": {c: out["cm"][c]["c_mode"] for c in scols},
            "c_in": c_in,
            "c_in_g": g_in,
            "c_in_iter": it,
            "z_k": {
                c: [[complex(z).real, complex(z).imag] for z in out["cm"][c]["zk"]] for c in scols
            },
            "off_diag": {c: out["cm"][c]["off"] for c in scols},
            "z64": {c: _cj(out["z64"][c]) for c in scols},
            "recip": out["recip"],
            "sigma_max": out["sigma_max"],
            "rows_at_c_mode": {c: full[c]["rows"] for c in cols},
            "delta_loo_at_c_mode": worst,
            "t_solve_s": t_solve,
        }
        if h:
            payload.update(
                {
                    "resolution": H_RES,
                    "columns": list(scols),
                    "c_in_x1": out["c_x1"],
                    "c_in_h_rel": out["c_h_rel"],
                    "tie": _tie_json(tie),
                }
            )
        _tune_json().write_text(json.dumps(payload))
    comm.Barrier()
    return out


@complex_only
def test_t21_c_construction(tune_cmode):
    """(C): Im z₁(64 MHz) > 0 (TH-17's assert); reciprocity / passivity per column;
    (F) on the solved mesh against item 196's cells."""
    _need_tune("cmode")
    _assert_tie(tune_cmode["cells"], tune_cmode["tie"], tune_cmode["res"])
    assert tune_cmode["cells"] == tune_cmode["stored_cells"], (
        tune_cmode["cells"],
        tune_cmode["stored_cells"],
    )
    for c in tune_cmode["scols"]:
        assert tune_cmode["cm"][c]["z1"].imag > 0.0, (c, tune_cmode["cm"][c]["z1"])
        assert tune_cmode["recip"][c] <= RECIPROCITY_BAND, (c, tune_cmode["recip"][c])
        assert tune_cmode["sigma_max"][c] <= 1.0 + PASSIVITY_SIGMA_TOLERANCE, (
            c,
            tune_cmode["sigma_max"][c],
        )


@complex_only
def test_t21_t1_route_b_at_c_mode(tune_cmode):
    """(T1): at C_mode,Cu both columns have one − → + zero, k-share >
    ``K_SHARE_FLOOR``, Δ_LOO <= ``BANDS``."""
    _need_tune("cmode")
    for c in tune_cmode["cols"]:
        r = tune_cmode["full"][c]["rows"]
        assert len(r["crossings_up_hz"]) == 1, (c, r["crossings_up_hz"])
        assert r["k_share"] > K_SHARE_FLOOR, (c, r["k_share"])
    for k in ROWS:
        assert tune_cmode["worst"][k] <= BANDS[k], (k, tune_cmode["worst"][k])


@complex_only
def test_t21_n1_cf_wiring(tune_cmode):
    """(N1), asserted: c_f=C_F reproduces item 196's record; untuned f0 is
    > N_FACTOR bands from 64 MHz (so (T2) is not vacuous)."""
    _need_tune("cmode")
    fails = []
    for c in BOTH_COLUMNS:
        for k in ROWS:
            got = tune_cmode["n1"][c]["rows"][k]
            want = ANS7_EXACT_STRIP_ROUTE_B_ROWS_P2[c][k]
            if abs(got / want - 1.0) > NODE_S_RTOL:
                fails.append((c, k, got, want))
    assert not fails, fails
    assert tune_cmode["n1_f0_off"] > N_FACTOR * BANDS["f0_hz"], tune_cmode["n1_f0_off"]


@pytest.fixture(scope="module")
def tune_inmodel():
    from tests.validation import test_th14_birdcage_copper as th14

    _need_tune("inmodel")
    # §9 item 310: the h = 0.01 rung runs per point (``tune_h``); a point knob off
    # the rung raises (``_tune_point``).
    if _tune_point() is not None or _h_tune():
        pytest.skip(f"{RESOLUTION_ENV}={H_RES:g}: T28 runs per point (test_t28_*)")
    # §9 item 311: under the ESR knob the drives run in ``tune_esr`` (test_t26_*) only.
    if _esr_ohm() is not None:
        pytest.skip(f"{ESR_ENV} set: T26 runs in tune_esr (test_t26_*)")
    comm = MPI.COMM_WORLD
    deg = _degree()
    act = _active()
    tgt, lab, tm, pr, tune = (
        act["target_hz"],
        act["lab"],
        act["tmhz"],
        _PRED[act["tmhz"]],
        act["tune"],
    )
    jt = _tune_json()
    if not jt.exists():
        pytest.fail(f"{jt} absent: run W1 (cmode) first")
    tj = json.loads(jt.read_text())
    if tj["degree"] != deg:
        pytest.fail(f"{jt}: degree {tj['degree']} != {deg}")
    if abs(tj["f_hz"] - tgt) > 1.0:
        pytest.fail(f"{jt}: f_hz {tj['f_hz']} != active target {tgt}")
    c_mode = float(tj["c_mode"]["cu"])
    cin = _tune_c() is not None
    c_drv = float(tj["c_in"]) if cin else c_mode
    cl = "C_in" if cin else "C_mode"
    z64_cu = _uj(tj["z64"]["cu"])
    recs, _ = _stored_active()
    z_nodes_cu = [_uj(r["z"]["cu"]) for r in recs]
    xn = np.asarray(act["nodes"], dtype=float) - 64.5
    _say(
        f"[ans7-exact] {lab} W2: C_mode,Cu read from {jt.name if _tune_mhz() is None else jt} ={c_mode!r} F; C_in "
        f"{tj['c_in']!r} F"
    )
    tie = _fixture_tie()
    t0 = time.perf_counter()
    geo, specs, cells = _exact_geo(th14, comm, GATED_WIDTH_FRACTION)
    t_mesh = time.perf_counter() - t0
    # (W) first: defaults (degree 1, C_F) against item 195's stored 64.0 MHz Cu point.
    t1 = time.perf_counter()
    s_w, _ = _route_a_solve(th14, geo, specs, TUNE_TARGET_HZ, "cu")
    rb = json.loads(_route_b_json(1).read_text())
    ra = json.loads((OUT_DIR / "route_a_p1.json").read_text())
    xn1 = np.asarray(rb["nodes_mhz"], dtype=float) - 64.5
    zb1 = route_b_z(
        bary_interp(xn1, [_uj(m) for m in rb["z"]["cu"]], 64.0 - 64.5), TUNE_TARGET_HZ, C_F
    )[1]
    s_b1 = _s11_from_z(zb1)
    res_w = float(abs(s_w - s_b1) / abs(s_b1))
    j64 = list(np.round(GRID_MHZ, 6)).index(64.0)
    res_195 = float(ra["residual"]["cu"][j64])
    dw = abs(res_w - res_195)
    t_w = time.perf_counter() - t1
    _say(
        f"[ans7-exact] T21 (W) degree-1 default drive 64.0 MHz Cu: residual vs stored Route B "
        f"{res_w!r}; item 195 stored {res_195!r}; |delta| = {dw:.3e} (ASSERTED <= NODE_S_RTOL "
        f"{NODE_S_RTOL:.0e}); {t_w:.1f} s"
    )
    res, s_as, s_bs, z_in_a, t_drive = [], [], [], [], []
    tfmt = ".1f" if _tune_mhz() is None else ".2f"
    for g in tune:
        f = g * 1e6
        omega = 2.0 * np.pi * f
        t2 = time.perf_counter()
        s_a, _ = _route_a_solve(th14, geo, specs, f, "cu", c_f=c_drv, degree=deg)
        t_drive.append(time.perf_counter() - t2)
        s_b = _s11_from_z(route_b_z(bary_interp(xn, z_nodes_cu, g - 64.5), f, c_drv)[1])
        r = float(abs(s_a - s_b) / abs(s_b))
        res.append(r)
        s_as.append(complex(s_a))
        s_bs.append(complex(s_b))
        z_in_a.append(_ANS7.add_c1(_ANS7._z_from_s11(s_a), omega, c_drv))
        if cin and _tune_mhz() is not None:
            plo, phi, ptxt = 1.4e-5, 1.8e-5, "[1.4e-5, 1.8e-5]"
        else:
            plo, phi = (1.0e-5, 1.5e-5) if cin else (5e-6, 3e-4)
            ptxt = "[1.0e-5, 1.5e-5]" if cin else "[5e-6, 3e-4]"
        msg = (
            f"[ans7-exact] {lab}{chr(8242) if cin else ''} (AB{chr(8242) if cin else ''}) "
            f"p{deg} Cu {g:{tfmt}} MHz at {cl}: S11_A "
            f"{complex(s_a):.12e}, "
            f"S11_B {complex(s_b):.12e}; |S11_A - S11_B|/|S11_B| = {r:.3e} (ASSERTED <= "
            f"{REDUCTION_BAND:.0e}; PREDICTED {ptxt}: {_met(plo <= r <= phi)}); "
            f"drive {t_drive[-1]:.1f} s"
        )
        if abs(f - tgt) < 1.0:
            s_b64 = _s11_from_z(route_b_z(z64_cu, f, c_drv)[1])
            msg += (
                f"; comparand from solved Z64 {complex(s_b64):.12e} (residual "
                f"{float(abs(s_a - s_b64) / abs(s_b64)):.3e}, printed)"
            )
        _say(msg)
    rows_a = parallel_rows(np.array(tune) * 1e6, np.array(z_in_a))
    b = (1.0 / np.asarray(z_in_a)).imag
    _say(
        f"[ans7-exact] {lab} (T2) Im Y_in,A at {tune} MHz = "
        + ", ".join(f"{v:+.9e}" for v in b)
        + f" S; rows A: f0 {rows_a.get('f0_hz')!r} Hz, R_in "
        f"{rows_a.get('r_in_ohm')!r} ohm, Q {rows_a.get('q_loaded')!r}; up-crossings "
        f"{rows_a['crossings_up_hz']}"
    )
    f0b = tj["rows_at_c_mode"]["cu"].get("f0_hz")
    t21p = None
    if cin:
        # §9 item 262: c, f0,B(C_in) on the stored nodes, r, (N1') and (N3').
        zc = {"cu": z_nodes_cu}
        c_rel = c_mode / c_drv - 1.0
        rb_ci = sparse_rows(act["nodes"], zc, grid_mhz=act["grid"], columns=("cu",), c_f=c_drv)[
            "cu"
        ]["rows"]
        f0b_ci = rb_ci.get("f0_hz")
        f0a = rows_a.get("f0_hz")
        sb = (f0b_ci / tgt - 1.0) if f0b_ci else float("nan")
        sa = (f0a / tgt - 1.0) if f0a else float("nan")
        r_ab = (f0a / f0b_ci - 1.0) if (f0a and f0b_ci) else float("nan")
        # §9 item 283: the 127.74 MHz set's tag, band and predictions (64 MHz bytes unchanged).
        is64 = _tune_mhz() is None
        pt = "T21'" if is64 else "T22'"
        bn = "T21P_ZERO_RESOLUTION" if is64 else "T22P_ZERO_RESOLUTION"
        band = _zero_resolution()
        c_pred = -1.494697e-5 if is64 else -7.754048e-6
        tp = "T2'" if is64 else "T2''"
        tmc = "64" if is64 else tm
        _say(
            f"[ans7-exact] {pt} C_in {c_drv!r} F (read from {jt.name}); c = C_mode/C_in - 1 = "
            f"{c_rel:+.6e} (PREDICTED {c_pred:+.6e}: {_met(abs(c_rel - c_pred) <= 1e-11)})"
        )
        _say(
            f"[ans7-exact] {pt} f0,B(C_in) {f0b_ci!r} Hz on the stored nodes; f0,B/{tmc} MHz - 1 = "
            f"{sb:+.6e} (PREDICTED within +-1e-9: {_met(abs(sb) <= 1e-9)}); up-crossings "
            f"{rb_ci.get('crossings_up_hz')}"
        )
        _say(
            f"[ans7-exact] {pt} ({tp}) f0,A(C_in) {f0a!r} Hz; f0,A/{tmc} MHz - 1 = {sa:+.6e} "
            f"(ASSERTED |.| <= {bn} {band:.3e}: "
            f"{_met(abs(sa) <= band)}; = {abs(sa) / band:.3f} "
            f"bands); r = f0,A/f0,B(C_in) - 1 = {r_ab:+.6e} (PREDICTED [-4e-6, -1.5e-6]: "
            f"{_met(-4e-6 <= r_ab <= -1.5e-6)})"
        )
        if is64:
            prev = (
                json.loads(_tune_inmodel_json().read_text())["rows_a"]
                if _tune_inmodel_json().exists()
                else {}
            )
            for k in ("r_in_ohm", "q_loaded"):
                if k in rows_a and k in prev:
                    d = rows_a[k] / prev[k] - 1.0
                    _say(
                        f"[ans7-exact] T21' {k} at C_in {rows_a[k]!r} vs T21 Route A at C_mode "
                        f"{prev[k]!r}: {d:+.4%} (PREDICTED within 0.1%: {_met(abs(d) <= 1e-3)})"
                    )
        else:
            # Item 283's comparand: T22's Route B rows at C_mode (``tune_p2.json``), 0.3 %.
            prev = tj["rows_at_c_mode"]["cu"]
            for k in ("r_in_ohm", "q_loaded"):
                if k in rows_a and k in prev:
                    d = rows_a[k] / prev[k] - 1.0
                    _say(
                        f"[ans7-exact] {pt} {k} at C_in {rows_a[k]!r} vs T22 Route B at C_mode "
                        f"{prev[k]!r}: {d:+.4%} (PREDICTED within 0.3%: {_met(abs(d) <= 3e-3)})"
                    )
        # (N1'), asserted at 64 MHz: the untuned zero on the stored nodes is >= 10 bands off.
        n1 = sparse_rows(act["nodes"], zc, grid_mhz=act["grid"], columns=("cu",), c_f=C_F)["cu"][
            "rows"
        ]
        n1_off = abs(n1["f0_hz"] / tgt - 1.0) if "f0_hz" in n1 else float("nan")
        # A missing key is None (fails (N0′)), never a vacuous [].
        n0x = list(n1["crossings_up_hz"]) if "crossings_up_hz" in n1 else None
        if is64:
            _say(
                f"[ans7-exact] T21' (N1') |f0,B(C_F)/64 MHz - 1| = {n1_off:.6e} = "
                f"{n1_off / T21P_ZERO_RESOLUTION:.1f} x T21P_ZERO_RESOLUTION (ASSERTED >= 10 x; "
                f"backing 2.029559e-2)"
            )
        else:
            # (N0'), asserted at 127.74: at C_F the stored Cu grid has no − → + zero.
            _say(
                f"[ans7-exact] {pt} (N0') at C_F on the {tm} MHz Cu grid: up-crossings {n0x} "
                f"(ASSERTED == []; backing item 228 W3 :2104)"
            )
        # (N3'), predicted, printed: a 1e-4 capacitor error is resolved by the band.
        n3 = sparse_rows(
            act["nodes"], zc, grid_mhz=act["grid"], columns=("cu",), c_f=c_drv * (1.0 + 1e-4)
        )["cu"]["rows"]
        n3_d = (n3["f0_hz"] / f0b_ci - 1.0) if ("f0_hz" in n3 and f0b_ci) else float("nan")
        if is64:
            _say(
                f"[ans7-exact] T21' (N3') f0,B(C_in*(1+1e-4))/f0,B(C_in) - 1 = {n3_d:+.6e} = "
                f"{abs(n3_d) / T21P_ZERO_RESOLUTION:.2f} x T21P_ZERO_RESOLUTION (PREDICTED ~ -5e-5, "
                f"~3.9 x: {_met(-7e-5 <= n3_d <= -3e-5)}; printed, never asserted)"
            )
        else:
            _say(
                f"[ans7-exact] {pt} (N3'') f0,B(C_in*(1+1e-4))/f0,B(C_in) - 1 = {n3_d:+.6e} = "
                f"{abs(n3_d) / band:.2f} x {bn} (PREDICTED ~ -4.7e-5, ~2.9 x: "
                f"{_met(-7e-5 <= n3_d <= -3e-5)}; printed, never asserted)"
            )
        t21p = {
            "c_in": c_drv,
            "c_rel": c_rel,
            "f0b_ci": f0b_ci,
            "sb": sb,
            "sa": sa,
            "r": r_ab,
            "n1_off": n1_off,
            "n3": n3_d,
        }
        if not is64:
            t21p["n0_crossings_up_hz"] = n0x
            t21p["band"] = band
    elif "f0_hz" in rows_a:
        da = abs(rows_a["f0_hz"] / tgt - 1.0)
        _say(
            f"[ans7-exact] {lab} (T2) |f0,A/{tm} MHz - 1| = {da:.6e} (ASSERTED <= "
            f"{BANDS['f0_hz']:.0e}; PREDICTED RED in [1e-6, 5e-5]: "
            f"{_met(1e-6 <= da <= 5e-5)}); f0,A - f0,B = "
            f"{(rows_a['f0_hz'] - f0b) if f0b else float('nan'):+.3f} Hz, relative "
            f"{(rows_a['f0_hz'] / f0b - 1.0) if f0b else float('nan'):+.3e} (PREDICTED "
            f"{pr['ab_gap']}); (C_mode - C_in)/C_in = {c_mode / tj['c_in'] - 1.0:+.3e}"
        )
        if _tune_mhz() is not None:
            sa = rows_a["f0_hz"] / tgt - 1.0
            sb = (f0b / tgt - 1.0) if f0b else float("nan")
            hyp = (
                "c flat"
                if da < 1e-5
                else "c tracks the off-diagonal" if da > 1.2e-5 else "undecided"
            )
            _say(
                f"[ans7-exact] {lab} (T2) split: f0,A/{tm} - 1 = {sa:+.6e} = c (f0,B/{tm} - 1 = "
                f"{sb:+.6e}) + r (f0,A/f0,B - 1 = "
                f"{(rows_a['f0_hz'] / f0b - 1.0) if f0b else float('nan'):+.6e}); hypothesis "
                f"readout: {hyp} (|f0,A/{tm} - 1| < 1e-5 c flat; > 1.2e-5 c tracks)"
            )
    # (N2) predicted, printed: circuit at 0.5*C_mode vs in-model at C_mode, target, Cu.
    j = [round(g, 6) for g in tune].index(round(tgt / 1e6, 6))
    s_n2 = _s11_from_z(route_b_z(z64_cu, tgt, 0.5 * c_drv)[1])
    n2 = float(abs(s_n2 - s_as[j]) / abs(s_as[j]))
    _say(
        f"[ans7-exact] {lab} (N2{chr(8242) if cin else ''}) circuit at 0.5*{cl} vs in-model at "
        f"{cl}, {tm} MHz Cu: "
        f"|dS|/|S| = {n2:.6e} = {n2 / REDUCTION_BAND:.2f} x REDUCTION_BAND (PREDICTED >= 10 x: "
        f"{_met(n2 >= 10.0 * REDUCTION_BAND)}; printed, never asserted)"
    )
    cin_rec = None
    # Item 283: under ``cin`` the three drives above already run at C_in, so item
    # 228's one-drive C_in record is not repeated.
    if _tune_mhz() is not None and not cin:
        # C_in record (printed, never asserted, never a target): one drive at C_in.
        c_in = float(tj["c_in"])
        t3 = time.perf_counter()
        s_ci, _ = _route_a_solve(th14, geo, specs, tgt, "cu", c_f=c_in, degree=deg)
        t_ci = time.perf_counter() - t3
        y_ci = complex(1.0 / _ANS7.add_c1(_ANS7._z_from_s11(s_ci), 2.0 * np.pi * tgt, c_in))
        dbdf = rows_a.get("dbdf_s_per_hz")
        dff = (-y_ci.imag / (tgt * dbdf)) if dbdf else float("nan")
        s_b_ci = _s11_from_z(route_b_z(z64_cu, tgt, c_in)[1])
        ab_ci = float(abs(s_ci - s_b_ci) / abs(s_b_ci))
        cin_rec = {
            "c_in": c_in,
            "im_y_a": y_ci.imag,
            "df_over_f_est": dff,
            "ab": ab_ci,
            "s11_a": [complex(s_ci).real, complex(s_ci).imag],
            "t_s": t_ci,
        }
        _say(
            f"[ans7-exact] {lab} C_in record: drive at C_in {c_in!r} F, {tm} MHz Cu: Im Y_in,A "
            f"{y_ci.imag:+.9e} S; estimate (one drive + slope) delta f/f = {dff:+.6e} (slope "
            f"dbdf {dbdf!r} S/Hz from the C_mode rows; PREDICTED |.| in [1e-7, 2e-5]: "
            f"{_met(1e-7 <= abs(dff) <= 2e-5)}; and outside BANDS {BANDS['f0_hz']:.0e}: "
            f"{_met(abs(dff) > BANDS['f0_hz'])}); |S11_A(C_in) - S11_B(C_in)|/|S11_B| = "
            f"{ab_ci:.3e}; drive {t_ci:.1f} s (printed, never asserted, never a target)"
        )
    _say(
        f"[ans7-exact] {lab} W2 total {time.perf_counter() - t0:.1f} s (mesh {t_mesh:.1f} s, (W) "
        f"{t_w:.1f} s, degree-{deg} drives {', '.join(f'{t:.1f}' for t in t_drive)} s), "
        f"cells {cells}, -n {comm.size}"
    )
    if comm.rank == 0:
        act["out"].mkdir(parents=True, exist_ok=True)
        payload = {
            "degree": deg,
            "ranks": comm.size,
            "cells": cells,
            "c_mode": c_mode,
            "tune_mhz": tune,
            "residual_ab": res,
            "s11_a": [[s.real, s.imag] for s in s_as],
            "s11_b": [[s.real, s.imag] for s in s_bs],
            "z_in_a": [[complex(z).real, complex(z).imag] for z in z_in_a],
            "rows_a": rows_a,
            "w_residual": res_w,
            "w_item195": res_195,
            "n2": n2,
            "t_drive_s": t_drive,
        }
        if cin_rec is not None:
            payload["c_in_record"] = cin_rec
            payload["target_hz"] = tgt
        if cin:
            payload["c_in"] = c_drv
            payload["t21p"] = t21p
            (act["out"] / f"tune_inmodel_cin_p{deg}.json").write_text(json.dumps(payload))
        else:
            _tune_inmodel_json().write_text(json.dumps(payload))
    comm.Barrier()
    return {
        "cells": cells,
        "tie": tie,
        "res": res,
        "rows_a": rows_a,
        "dw": dw,
        "n2": n2,
        "stored_cells": recs[0]["cells"],
        "tune": tune,
        "tgt": tgt,
        "t21p": t21p,
    }


@complex_only
def test_t21_w_degree_default_wiring(tune_inmodel):
    """(W), asserted: the degree-1 default drive reproduces item 195's stored
    64.0 MHz Cu residual within ``NODE_S_RTOL`` (absolute)."""
    _need_tune("inmodel")
    _assert_tie(tune_inmodel["cells"], tune_inmodel["tie"])
    assert tune_inmodel["dw"] <= NODE_S_RTOL, tune_inmodel["dw"]


@complex_only
def test_t21_ab_circuit_equals_field_p2(tune_inmodel):
    """(AB): degree 2, Cu, at C_mode -- |S11_A - S11_B|/|S11_B| <=
    ``REDUCTION_BAND`` at all three ``TUNE_MHZ``."""
    _need_tune("inmodel")
    if _tune_c() is not None:
        pytest.skip(f"{TUNE_C_ENV}=cin: (AB) asserts at C_mode")
    assert tune_inmodel["cells"] == tune_inmodel["stored_cells"]
    assert len(tune_inmodel["res"]) == len(tune_inmodel["tune"])
    for g, r in zip(tune_inmodel["tune"], tune_inmodel["res"]):
        assert r <= REDUCTION_BAND, (g, r)


@complex_only
def test_t21_t2_inmodel_zero(tune_inmodel):
    """(T2): exactly one − → + zero of Im Y_in,A in (63.9, 64.1) MHz and
    |f₀,A/64 MHz − 1| <= ``BANDS["f0_hz"]`` (the band half is the brief's
    PREDICTED RED, range [1e-6, 5e-5])."""
    _need_tune("inmodel")
    if _tune_c() is not None:
        pytest.skip(f"{TUNE_C_ENV}=cin: (T2) asserts at C_mode")
    r = tune_inmodel["rows_a"]
    assert len(r["crossings_up_hz"]) == 1, r["crossings_up_hz"]
    lo, hi = tune_inmodel["tune"][0] * 1e6, tune_inmodel["tune"][-1] * 1e6
    assert lo < r["f0_hz"] < hi, r["f0_hz"]
    assert abs(r["f0_hz"] / tune_inmodel["tgt"] - 1.0) <= BANDS["f0_hz"], r["f0_hz"]


def _need_cin():
    if _tune_c() is None:
        pytest.skip(f"{TUNE_C_ENV} unset: item 262's tests run on 'cin'")
    _need_tune("inmodel")


@complex_only
def test_t21p_ab_circuit_equals_field_cin(tune_inmodel):
    """(AB′): degree 2, Cu, at C_in -- |S11_A - S11_B|/|S11_B| <=
    ``REDUCTION_BAND`` at all three ``TUNE_MHZ``."""
    _need_cin()
    assert tune_inmodel["cells"] == tune_inmodel["stored_cells"]
    assert len(tune_inmodel["res"]) == len(tune_inmodel["tune"]) == 3
    for g, r in zip(tune_inmodel["tune"], tune_inmodel["res"]):
        assert r <= REDUCTION_BAND, (g, r)


@complex_only
def test_t21p_inmodel_zero_cin(tune_inmodel):
    """(T2′): exactly one − → + zero of Im Y_in,A in (63.9, 64.1) MHz at C_in and
    |f₀,A/64 MHz − 1| <= ``T21P_ZERO_RESOLUTION`` (pre-registered literal).
    §9 item 283 (T2″) at 127.74: bracket (127.64, 127.84) MHz, target 127.74 MHz,
    band ``T22P_ZERO_RESOLUTION`` -- all from the knob, never from a reading."""
    _need_cin()
    r = tune_inmodel["rows_a"]
    assert len(r["crossings_up_hz"]) == 1, r["crossings_up_hz"]
    lo, hi = tune_inmodel["tune"][0] * 1e6, tune_inmodel["tune"][-1] * 1e6
    assert lo < r["f0_hz"] < hi, r["f0_hz"]
    assert abs(r["f0_hz"] / _active()["target_hz"] - 1.0) <= _zero_resolution(), r["f0_hz"]


@complex_only
def test_t21p_n1_untuned(tune_inmodel):
    """(N1′), asserted: the untuned (C_F) zero on the stored nodes is >= 10 ×
    ``T21P_ZERO_RESOLUTION`` from 64 MHz, so (T2′) is not vacuous.  §9 item 283
    (N0′) at 127.74: at C_F the stored Cu grid has no − → + zero at all."""
    _need_cin()
    if _tune_mhz() is not None:
        n0 = tune_inmodel["t21p"]["n0_crossings_up_hz"]
        assert n0 == [], n0
        return
    off = tune_inmodel["t21p"]["n1_off"]
    assert off >= 10.0 * T21P_ZERO_RESOLUTION, off


# -- §9 item 310 (`ANS-7` step 12, chain T28): the per-point drives + assembly --


def _h_point_json(mhz):
    if _tune_mhz() == 127.74:
        # §9 item 365 (T34).
        return TUNE128_OUT_DIR / f"tune_inmodel_cin_p2{H_SFX}_{mhz:.2f}.json"
    return OUT_DIR / f"tune_inmodel_cin_p2{H_SFX}_{mhz:.1f}.json"


def _h_assemble(tj, c_drv, z64_cu, z_nodes_cu, recs):
    """Pure: read the three point files, rows, (T2‴) crossing, printed controls."""
    comm = MPI.COMM_WORLD
    # §9 item 365 (T34): every 64 MHz literal is read from ``_active()``; the
    # 64 MHz print bytes are unchanged (lab "T28", ".1f" points, "64 MHz").
    act = _active()
    at128 = _tune_mhz() == 127.74
    lab, tune = act["lab"], act["tune"]
    nodes, grid = act["nodes"], act["grid"]
    tml = act["tmhz"] if at128 else "64"
    zres = _zero_resolution()
    zres_name = "T22P_ZERO_RESOLUTION" if at128 else "T21P_ZERO_RESOLUTION"

    def _pf(g):
        return f"{g:.2f}" if at128 else f"{g:.1f}"

    paths = [_h_point_json(g) for g in tune]
    missing = [p.name for p in paths if not p.exists()]
    if missing:
        pytest.fail(f"{lab} point JSONs absent: {missing} (run the point windows first)")
    pts = [json.loads(p.read_text()) for p in paths]
    for g, p in zip(tune, pts):
        assert abs(p["f_hz"] - g * 1e6) < 1.0, (g, p["f_hz"])
        assert p["c_in"] == c_drv, (g, p["c_in"], c_drv)
        assert p["cells"] == recs[0]["cells"], (g, p["cells"], recs[0]["cells"])
        assert p["degree"] == 2 and p["resolution"] == H_RES, (g, p["degree"], p["resolution"])
    res = [float(p["residual_ab"]) for p in pts]
    s_as = [complex(*p["s11_a"]) for p in pts]
    z_in_a = np.array([complex(*p["z_in_a"]) for p in pts])
    tgt = act["target_hz"]
    rows_a = parallel_rows(np.array(tune) * 1e6, z_in_a)
    b = (1.0 / z_in_a).imag
    n_up = sum(1 for k in range(len(b) - 1) if b[k] < 0 <= b[k + 1])
    n_dn = sum(1 for k in range(len(b) - 1) if b[k] >= 0 > b[k + 1])
    for g, r in zip(tune, res):
        _say(
            f"[ans7-exact] {lab} (AB''') p2 Cu {_pf(g)} MHz at C_in(h): |S11_A - S11_B|/|S11_B| = "
            f"{r:.3e} (ASSERTED <= REDUCTION_BAND {REDUCTION_BAND:.0e}; PREDICTED [5e-6, 3e-5]: "
            f"{_met(5e-6 <= r <= 3e-5)})"
        )
    _say(
        f"[ans7-exact] {lab} (T2''') Im Y_in,A at {tune} MHz = "
        + ", ".join(f"{v:+.9e}" for v in b)
        + f" S; - -> + sign changes {n_up}, + -> - {n_dn} "
        f"(ASSERTED exactly one - -> +, none + -> -); rows A: f0 {rows_a.get('f0_hz')!r} Hz, "
        f"R_in {rows_a.get('r_in_ohm')!r} ohm, Q {rows_a.get('q_loaded')!r}; up-crossings "
        f"{rows_a['crossings_up_hz']}"
    )
    f0a = rows_a.get("f0_hz")
    sa = (f0a / tgt - 1.0) if f0a else float("nan")
    _say(
        f"[ans7-exact] {lab} |f0,A/{tml} MHz - 1| = {abs(sa):.6e} (printed, never asserted; "
        f"{zres_name} {zres:.3e} beside it for scale only = "
        f"{abs(sa) / zres:.3f} x)"
    )
    zc = {"cu": z_nodes_cu}
    rb_ci = sparse_rows(nodes, zc, grid_mhz=grid, columns=("cu",), c_f=c_drv)["cu"]["rows"]
    f0b_ci = rb_ci.get("f0_hz")
    sb = (f0b_ci / tgt - 1.0) if f0b_ci else float("nan")
    r_ab = (f0a / f0b_ci - 1.0) if (f0a and f0b_ci) else float("nan")
    _say(
        f"[ans7-exact] {lab} f0,B(C_in(h)) {f0b_ci!r} Hz on the h nodes; f0,B/{tml} MHz - 1 = "
        f"{sb:+.6e} (printed); r(h) = f0,A/f0,B(C_in(h)) - 1 = {r_ab:+.6e} (PREDICTED "
        f"[-5e-6, -1e-6]: {_met(-5e-6 <= r_ab <= -1e-6)}); up-crossings "
        f"{rb_ci.get('crossings_up_hz')}"
    )
    if not at128:
        for k in ("r_in_ohm", "q_loaded"):
            if k in rows_a:
                d = rows_a[k] / STEP6_RUNG_ROWS[k] - 1.0
                _say(
                    f"[ans7-exact] T28 {k} at C_in(h) {rows_a[k]!r} vs step 6 rung row "
                    f"{STEP6_RUNG_ROWS[k]!r}: {d:+.4%} (PREDICTED within 1%: {_met(abs(d) <= 1e-2)}; "
                    f"printed)"
                )
    else:
        # §9 item 365: the x1 127.74 rows at C_in (item 283), read, never recomputed.
        x1rows = json.loads((TUNE128_OUT_DIR / "tune_inmodel_cin_p2.json").read_text())["rows_a"]
        for k in ("r_in_ohm", "q_loaded"):
            if k in rows_a and k in x1rows:
                d = rows_a[k] / x1rows[k] - 1.0
                _say(
                    f"[ans7-exact] {lab} {k} at C_in(h) {rows_a[k]!r} vs x1 127.74 row at C_in "
                    f"{x1rows[k]!r}: {d:+.4%} (PREDICTED within 2%: {_met(abs(d) <= 2e-2)}; "
                    f"printed)"
                )
    if not at128:
        n1 = sparse_rows(nodes, zc, grid_mhz=grid, columns=("cu",), c_f=C_F)["cu"]["rows"]
        n1_off = abs(n1["f0_hz"] / tgt - 1.0) if "f0_hz" in n1 else float("nan")
        _say(
            f"[ans7-exact] T28 (N1''') |f0,B(C_F)/64 MHz - 1| on the h grid = {n1_off:.6e} = "
            f"{n1_off / T21P_ZERO_RESOLUTION:.1f} x T21P_ZERO_RESOLUTION (PREDICTED >= 10 x, ~1500 x: "
            f"{_met(n1_off >= 10.0 * T21P_ZERO_RESOLUTION)}; printed, never asserted); f0,B(C_F) "
            f"{n1.get('f0_hz')!r} Hz (step 6 {STEP6_RUNG_ROWS['f0_hz']!r})"
        )
    else:
        # (N0⁗): the up-crossings at C_F on the h 127.74 grid, printed, never asserted.
        n0 = sparse_rows(nodes, zc, grid_mhz=grid, columns=("cu",), c_f=C_F)["cu"]["rows"]
        n0x = n0.get("crossings_up_hz", [])
        n1_off = float("nan")
        _say(
            f"[ans7-exact] {lab} (N0'''') at C_F on the h {tml} MHz Cu grid: up-crossings {n0x} "
            f"(PREDICTED none: {_met(len(n0x) == 0)}; printed, never asserted)"
        )
    mid = tune[1]
    j = list(tune).index(mid)
    s_n2 = _s11_from_z(route_b_z(z64_cu, tgt, 0.5 * c_drv)[1])
    n2 = float(abs(s_n2 - s_as[j]) / abs(s_as[j]))
    _say(
        f"[ans7-exact] {lab} (N2''') circuit at 0.5*C_in(h) vs in-model at C_in(h), {tml} MHz Cu: "
        f"|dS|/|S| = {n2:.6e} = {n2 / REDUCTION_BAND:.2f} x REDUCTION_BAND (PREDICTED >= 10 x: "
        f"{_met(n2 >= 10.0 * REDUCTION_BAND)}; printed, never asserted)"
    )
    c_h_rel = float(tj["c_in_h_rel"])
    if not at128:
        _say(
            f"[ans7-exact] T28 C_in(h)/C_in(x1) - 1 = {c_h_rel:+.6e} (from {_tune_json().name}); "
            f"step 6 f0 move at fixed C {STEP6_RUNG_F0_MOVE:+.4%}"
        )
    else:
        _say(
            f"[ans7-exact] {lab} C_in(h)/C_in(x1) - 1 at {tml} MHz = {c_h_rel:+.6e} (from "
            f"{_tune_json().name}); beside 64 MHz's record {T28_C_IN_H_REL_64:+.6e}"
        )
    t28 = {
        "c_in": c_drv,
        "c_in_x1": float(tj["c_in_x1"]),
        "c_in_h_rel": c_h_rel,
        "residual_ab": res,
        "im_y_a": [float(v) for v in b],
        "n_up": n_up,
        "n_dn": n_dn,
        "rows_a": rows_a,
        "f0a_rel": sa,
        "f0b_ci": f0b_ci,
        "sb": sb,
        "r": r_ab,
        "n1_off": n1_off,
        "n2": n2,
        "w_residual": pts[j].get("w_residual"),
        "w_item195": pts[j].get("w_item195"),
        "cells": recs[0]["cells"],
        "t_point_s": [p["t_window_s"] for p in pts],
    }
    if comm.rank == 0:
        (act["out"] / f"tune_inmodel_cin_p2{H_SFX}.json").write_text(json.dumps(t28))
    comm.Barrier()
    return {
        "mode": "assemble",
        "res": res,
        "b": list(b),
        "n_up": n_up,
        "n_dn": n_dn,
        "rows_a": rows_a,
        "dw": abs(pts[j]["w_residual"] - pts[j]["w_item195"]),
    }


@pytest.fixture(scope="module")
def tune_h():
    from tests.validation import test_th14_birdcage_copper as th14

    _need_tune("inmodel")
    if not _h_tune():
        pytest.skip(f"{RESOLUTION_ENV} unset: T28 runs on the h = {H_RES:g} rung")
    # §9 item 395 (T37): under the ESR knob the rung's drives run in ``tune_esr``.
    if _esr_ohm() is not None:
        pytest.skip(f"{ESR_ENV} set: T37 runs in tune_esr (test_t26_*)")
    pt = _tune_point()
    if _tune_c() != "cin":
        raise ValueError(f"T28 needs {TUNE_C_ENV}=cin")
    comm = MPI.COMM_WORLD
    deg = _degree()
    if deg != 2:
        raise ValueError(f"T28 runs at degree 2 only, got {deg}")
    jt = _tune_json()
    if not jt.exists():
        pytest.fail(f"{jt} absent: run W1 (cmode on the rung) first")
    tj = json.loads(jt.read_text())
    # §9 item 365 (T34): the active target, tune points and nodes (64 MHz unchanged).
    act = _active()
    lab, mid = act["lab"], act["tune"][1]
    at128 = _tune_mhz() == 127.74

    def _pf(g):
        return f"{g:.2f}" if at128 else f"{g:.1f}"

    if (
        tj["degree"] != 2
        or abs(tj["f_hz"] - act["target_hz"]) > 1.0
        or tj.get("resolution") != H_RES
    ):
        pytest.fail(f"{jt}: degree {tj['degree']} f_hz {tj['f_hz']} res {tj.get('resolution')}")
    c_drv = float(tj["c_in"])
    z64_cu = _uj(tj["z64"]["cu"])
    recs, _ = _stored_active()
    z_nodes_cu = [_uj(r["z"]["cu"]) for r in recs]
    if tj["cells"] != recs[0]["cells"]:
        pytest.fail(f"{jt}: cells {tj['cells']} != node cells {recs[0]['cells']}")
    _say(
        f"[ans7-exact] {lab} C_in(h) read from {jt.name} = {c_drv!r} F; point {pt!r}; -n {comm.size}"
    )
    if pt == "assemble":
        return _h_assemble(tj, c_drv, z64_cu, z_nodes_cu, recs)
    t0 = time.perf_counter()
    f = pt * 1e6
    omega = 2.0 * np.pi * f
    w = None
    if abs(pt - mid) < 1e-9:
        # (W): item 262's wiring drive, unchanged -- degree 1, C_F, on the x1 mesh
        # (64 MHz, x1, in the middle-point window at either target; §9 item 365).
        t1 = time.perf_counter()
        geo1, specs1, _c1 = _exact_geo(th14, comm, GATED_WIDTH_FRACTION)
        s_w, _ = _route_a_solve(th14, geo1, specs1, TUNE_TARGET_HZ, "cu")
        del geo1, specs1
        rb = json.loads(_route_b_json(1).read_text())
        ra = json.loads((OUT_DIR / "route_a_p1.json").read_text())
        xn1 = np.asarray(rb["nodes_mhz"], dtype=float) - 64.5
        zb1 = route_b_z(
            bary_interp(xn1, [_uj(m) for m in rb["z"]["cu"]], 64.0 - 64.5), TUNE_TARGET_HZ, C_F
        )[1]
        s_b1 = _s11_from_z(zb1)
        res_w = float(abs(s_w - s_b1) / abs(s_b1))
        j64 = list(np.round(GRID_MHZ, 6)).index(64.0)
        res_195 = float(ra["residual"]["cu"][j64])
        w = {
            "res_w": res_w,
            "res_195": res_195,
            "dw": abs(res_w - res_195),
            "t": time.perf_counter() - t1,
        }
        _say(
            f"[ans7-exact] {lab} (W) degree-1 default drive 64.0 MHz Cu (x1 mesh): residual vs "
            f"stored Route B {res_w!r}; item 195 stored {res_195!r}; |delta| = {w['dw']:.3e} "
            f"(ASSERTED <= NODE_S_RTOL {NODE_S_RTOL:.0e}); {w['t']:.1f} s"
        )
    t2 = time.perf_counter()
    geo, specs, cells = _exact_geo(th14, comm, GATED_WIDTH_FRACTION, res=H_RES)
    t_mesh = time.perf_counter() - t2
    t3 = time.perf_counter()
    s_a, _ = _route_a_solve(th14, geo, specs, f, "cu", c_f=c_drv, degree=deg)
    t_drive = time.perf_counter() - t3
    xn = np.asarray(act["nodes"], dtype=float) - 64.5
    s_b = _s11_from_z(route_b_z(bary_interp(xn, z_nodes_cu, pt - 64.5), f, c_drv)[1])
    r = float(abs(s_a - s_b) / abs(s_b))
    z_in = _ANS7.add_c1(_ANS7._z_from_s11(s_a), omega, c_drv)
    im_y = float((1.0 / complex(z_in)).imag)
    msg = (
        f"[ans7-exact] {lab} (AB''') p{deg} Cu {_pf(pt)} MHz at C_in(h): S11_A {complex(s_a):.12e}, "
        f"S11_B {complex(s_b):.12e}; |S11_A - S11_B|/|S11_B| = {r:.3e} (ASSERTED <= "
        f"{REDUCTION_BAND:.0e}; PREDICTED [5e-6, 3e-5]: {_met(5e-6 <= r <= 3e-5)}); Im Y_in,A "
        f"{im_y:+.9e} S; mesh {t_mesh:.1f} s, drive {t_drive:.1f} s, cells {cells}"
    )
    if abs(pt - mid) < 1e-9:
        s_b64 = _s11_from_z(route_b_z(z64_cu, f, c_drv)[1])
        msg += (
            f"; comparand from solved Z64(h) {complex(s_b64):.12e} (residual "
            f"{float(abs(s_a - s_b64) / abs(s_b64)):.3e}, printed)"
        )
    _say(msg)
    t_win = time.perf_counter() - t0
    if comm.rank == 0:
        p = {
            "f_hz": f,
            "degree": deg,
            "ranks": comm.size,
            "cells": cells,
            "resolution": H_RES,
            "c_in": c_drv,
            "residual_ab": r,
            "s11_a": [complex(s_a).real, complex(s_a).imag],
            "s11_b": [complex(s_b).real, complex(s_b).imag],
            "z_in_a": [complex(z_in).real, complex(z_in).imag],
            "im_y_a": im_y,
            "t_mesh_s": t_mesh,
            "t_drive_s": t_drive,
            "t_window_s": t_win,
        }
        if w is not None:
            p.update({"w_residual": w["res_w"], "w_item195": w["res_195"], "t_w_s": w["t"]})
        _h_point_json(pt).write_text(json.dumps(p))
    comm.Barrier()
    _say(f"[ans7-exact] {lab} point {_pf(pt)} MHz window {t_win:.1f} s, -n {comm.size}")
    return {
        "mode": "point",
        "pt": pt,
        "cells": cells,
        "stored_cells": recs[0]["cells"],
        "res": [r],
        "dw": None if w is None else w["dw"],
    }


@complex_only
def test_t28_rung_and_wiring(tune_h):
    """(F) on the rung: the drive mesh is item 214's (cells equal the node files'
    and moved off the x1 record, (N-wire)); (W), asserted, in the 64.0 window."""
    if tune_h["mode"] == "point":
        assert tune_h["cells"] == tune_h["stored_cells"], (tune_h["cells"], tune_h["stored_cells"])
        assert abs(tune_h["cells"] / ANS6_EXACT_STRIP_HOLE_CELL_RECORD - 1.0) > CELL_COUNT_BAND
    if tune_h["dw"] is not None:
        assert tune_h["dw"] <= NODE_S_RTOL, tune_h["dw"]


@complex_only
def test_t28_ab_circuit_equals_field(tune_h):
    """(AB‴): |S11_A − S11_B|/|S11_B| <= ``REDUCTION_BAND`` at C_in(h) -- the
    window's point, or all three at assembly."""
    if tune_h["mode"] == "assemble":
        assert len(tune_h["res"]) == len(_active()["tune"])
    for r in tune_h["res"]:
        assert r <= REDUCTION_BAND, tune_h["res"]


@complex_only
def test_t28_t2_zero_in_bracket(tune_h):
    """(T2‴), ruling 310-Q1: Im Y_in,A changes sign − → + exactly once between
    63.9 and 64.1 MHz (and never + → −); the offset is printed, never asserted."""
    if tune_h["mode"] != "assemble":
        pytest.skip("(T2''') asserts at assembly")
    b = tune_h["b"]
    assert b[0] < 0.0 < b[-1], b
    assert tune_h["n_up"] == 1 and tune_h["n_dn"] == 0, b
    r = tune_h["rows_a"]
    assert len(r["crossings_up_hz"]) == 1, r["crossings_up_hz"]
    assert _active()["tune"][0] * 1e6 < r["f0_hz"] < _active()["tune"][-1] * 1e6, r["f0_hz"]


# -- §9 item 311 (`ANS-7` step 11, chain T26): tune on C_in with r = 0.2 Ω ESR --
#
# ``FEM_EM_ANS7_EXACT_ESR_OHM`` ∈ {unset, ``0.2``}: unset ⇒ every path above is
# unchanged (``tune_inmodel`` gains only a skip under the knob).  ``0.2`` (needs
# ``TUNE=inmodel``, ``TUNE_C=cin``, ``TUNE_MHZ`` / ``RESOLUTION`` / ``TUNE_POINT``
# unset) puts r in series with every ring-gap capacitor on both sides at once
# (operator ESR ruling 2026-10-06; ruling 311-Q5: the feed's branch carries r
# too).  Route A's terminated sheets are item 309's series-RC sheet
# ``series_rlc_impedance(f, r_ohm=r, l_h=0.0, c_f=C)``; Route B and the feed's
# (r + C) branch are formed here, test-side (2026-10-11 review amendment: item
# 363 lands ``esr_ohm`` only in ``ports/ring_modes.py``; ``ans:7`` is not edited).
# Reads ``tune_p2.json`` / ``tune_inmodel_cin_p2.json`` only; writes
# ``tune_p2_esr0.2.json`` and ``tune_inmodel_cin_esr0.2_p2.json``.
ESR_ENV = "FEM_EM_ANS7_EXACT_ESR_OHM"
ESR_VALUES = ("0.2",)
# §9 item 387 (T35): one record literal, printed at 127.74 MHz only and never
# asserted -- T26's C_in(r)/C_in(0) - 1 at 64 MHz
# (``20261011T094028Z_ANS-7-item311-W2-cin64-esr.log:359``).
T26_C_IN_R_REL_64 = 4.602610e-5
# §9 item 395 (T37): x1 64 MHz ESR record literals, printed on the h path only and
# never asserted (W0 rewrites item 311's x1 ESR JSONs in the same slot, so the
# comparands come from the registered log, not the files).
# C_in(r, x1): ``20261011T094028Z_ANS-7-item311-W2-cin64-esr.log:359``.
T26_C_IN_R_X1 = 6.359301781809668e-11
# Q(r, x1): ``20261011T094028Z_ANS-7-item311-W2-cin64-esr.log:2130``.
T26_Q_R_X1 = 164.5434331710199
# Q(r, x1)/closed form - 1: ``20261011T094028Z_ANS-7-item311-W2-cin64-esr.log:2133``.
T26_Q_CLOSED_EXCESS_X1 = 1.0996e-2
# Q_0 at the x1 C_in: ``20261011T094028Z_ANS-7-item311-W2-cin64-esr.log:2133``.
T21P_Q0_X1 = 971.0667834101893
# T30's Q h-move at C_F (h = 0.01 vs x1): ``20261011T094910Z_ANS-7-item364-T30.log:76``
# prints the magnitude 6.9811e-03; the sign (h below x1) is the item 395 brief's.
T30_Q_H_REL_CF = -6.9811e-3
# R_in(r, x1): ``20261011T094028Z_ANS-7-item311-W2-cin64-esr.log:2130``.
T26_R_IN_R_X1 = 3175.979466188352

# §9 item 419 (T38): x1 127.74 MHz ESR record literals, printed on the h38 path only
# and never asserted (W0 rewrites item 387's x1 ESR JSONs in the same slot, so the
# comparands come from the registered logs, not the files).
# C_in(r, x1) at 127.74 MHz: ``20261012T004412Z_ANS-7-item387-W2-esr128.log:291``.
T35_C_IN_R_X1 = 1.5275048999871143e-11
# Q(r, x1) at 127.74 MHz: ``20261012T004412Z_ANS-7-item387-W2-esr128.log:2130``.
T35_Q_R_X1 = 254.46925911792314
# R_in(r, x1) at 127.74 MHz: ``20261012T004412Z_ANS-7-item387-W2-esr128.log:2130``.
T35_R_IN_R_X1 = 9790.703605847839
# Q(r, x1)/closed form - 1 at 127.74 MHz: ``20261012T004412Z_ANS-7-item387-W2-esr128.log:2133``.
T35_Q_CLOSED_EXCESS_X1 = 3.5494e-2
# Q_0 at the x1 127.74 MHz C_in: ``20261012T004412Z_ANS-7-item387-W2-esr128.log:2133``.
T22P_Q0_X1 = 618.338591513781
# Lossless C_in(h)/C_in(x1) - 1 at 127.74 MHz (§10's "+8.968e-4"):
# ``20261011T112734Z_ANS-7-item365b-A-assemble.log:113``.
T34_C_IN_H_REL_128 = 8.967920e-4


def _esr_ohm():
    raw = os.environ.get(ESR_ENV, "").strip()
    if not raw:
        return None
    if raw not in ESR_VALUES:
        raise ValueError(f"{ESR_ENV}={raw!r}: only unset or {ESR_VALUES}")
    # §9 item 387 (T35): 127.74 MHz is the second legal set, and needs COLUMNS=cu
    # (the stored 127.74 nodes are ``_cu``-suffixed; item 283's trap).
    # §9 item 395 (T37): the third legal set is the h = 0.01 rung at 64 MHz, Cu only
    # (item 214's nodes are ``_cu``-suffixed); TUNE_POINT stays unset in every set.
    res_raw = os.environ.get(RESOLUTION_ENV, "").strip()
    h37 = bool(res_raw) and _tune_mhz() is None and _resolution() == H_RES and _cols() == ("cu",)
    # §9 item 419 (T38): the fourth legal set is the h = 0.01 rung at 127.74 MHz, Cu only
    # (T34's h nodes are ``_cu``-suffixed).
    h38 = bool(res_raw) and _tune_mhz() == 127.74 and _resolution() == H_RES and _cols() == ("cu",)
    if (
        _tune() != "inmodel"
        or _tune_c() != "cin"
        or _tune_mhz() not in (None, 127.74)
        or (_tune_mhz() == 127.74 and _cols() != ("cu",))
        or (res_raw and not (h37 or h38))
        or os.environ.get(TUNE_POINT_ENV, "").strip()
    ):
        raise ValueError(
            f"{ESR_ENV}={raw!r} needs {TUNE_ENV}=inmodel, {TUNE_C_ENV}=cin, "
            f"{TUNE_POINT_ENV} unset, and one of: {RESOLUTION_ENV} and "
            f"{TUNE_MHZ_ENV} unset (64 MHz, T26); {RESOLUTION_ENV} unset and "
            f"{TUNE_MHZ_ENV}=127.74 with {COLUMNS_ENV}=cu (T35); "
            f"{RESOLUTION_ENV}={H_RES:g}, {TUNE_MHZ_ENV} unset, {COLUMNS_ENV}=cu "
            f"(64 MHz on the h rung, T37); or {RESOLUTION_ENV}={H_RES:g}, "
            f"{TUNE_MHZ_ENV}=127.74, {COLUMNS_ENV}=cu (127.74 MHz on the h rung, T38)"
        )
    return float(raw)


def _z_branch(omega, c_f, r):
    """One ring-gap branch, r in series with C: r + 1/(jωC) (ESR ruling (iii): never r ‖ C)."""
    return r + 1.0 / (1j * omega * c_f)


def _route_b_z_esr(z4, f_hz, c_f, r):
    """``ans:7``'s ``route_b_z`` (``07_…_driven_sweep.py:181–188``, commit 2f0f83d1) with
    each Z_C replaced by r + Z_C, and the feed's branch added in parallel likewise.
    Returns Z_in, Z_bare, gap currents."""
    omega = 2.0 * np.pi * f_hz
    zb = _z_branch(omega, c_f, r)
    z4 = np.asarray(z4, dtype=np.complex128)
    z_tt = z4[1:, 1:] + zb * np.eye(3)
    i_t = -np.linalg.solve(z_tt, z4[1:, 0])
    z_bare = z4[0, 0] + z4[0, 1:] @ i_t
    return 1.0 / (1.0 / z_bare + 1.0 / zb), z_bare, np.concatenate(([1.0 + 0j], i_t))


def _c_in_esr(z4, f_hz, c0, r):
    """``_c_in`` (above) on the ESR Route B: the root of Im[1/Z_in(C; r)] at
    ``f_hz``, real-scalar secant from ``c0``, to |g| <= 1e-15 S in <= 20 iterations."""

    def g(c):
        return float((1.0 / complex(_route_b_z_esr(z4, f_hz, c, r)[0])).imag)

    ca, cb = float(c0), float(c0) * (1.0 + 1e-6)
    ga, gb = g(ca), g(cb)
    it = 0
    while abs(gb) > 1e-15 and it < 20:
        if gb == ga:
            break
        ca, cb = cb, cb - gb * (cb - ca) / (gb - ga)
        ga, gb = gb, g(cb)
        it += 1
    return cb, gb, it


def _route_a_solve_esr(th14, geo, specs, f_hz, c_f, r, degree):
    """Copied from ``ans:7``'s ``_route_a_solve`` (``07_tuned_copper_birdcage_driven_sweep.py
    :541–574``, commit 2f0f83d1), Cu column only, no fields; the one change is the
    terminated sheets' impedance: item 309's series-RC sheet
    ``series_rlc_impedance(f, r_ohm=r, l_h=0.0, c_f=C)`` in place of ``z_cap``."""
    import dataclasses

    from fem_em_solver.core.cavity import surface_resistance_ohm
    from fem_em_solver.ports.lumped import run_lumped_sheet_port_case, series_rlc_impedance
    from fem_em_solver.ports.sparameters import _power_waves

    msh, cell_tags, facet_tags, outer_tags, tags_f, port_defs = geo
    omega = 2.0 * np.pi * f_hz
    z_c = series_rlc_impedance(f_hz, r_ohm=r, l_h=0.0, c_f=c_f)
    sp = [
        (
            dataclasses.replace(s, port_impedance_ohm=z_c)
            if i in _ANS7.TUNING_TERMINATED_INDICES
            else s
        )
        for i, s in enumerate(specs)
    ]
    problem = th14._problem(msh, cell_tags, f_hz, outer_tags, (th14.OUTER_BOX_TAG,))
    z_s = (1.0 + 1.0j) * surface_resistance_ohm(omega, th14.COPPER_SIGMA)
    extra = [th14._leontovich_term(msh, facet_tags, z_s, omega)]
    result = run_lumped_sheet_port_case(
        problem,
        port_defs,
        sp,
        facet_tags=tags_f,
        driven_port_id=th14.DRIVEN,
        verbose=False,
        return_fields=False,
        extra_bilinear_terms=extra,
        degree=degree,
    )
    drive = result.responses[th14.DRIVEN]
    a1, b1 = _power_waves(drive.voltage_v, drive.current_a, float(_ANS7.REFERENCE_IMPEDANCE_OHM))
    assert abs(a1) > 0.0, "incident wave at P1 vanished"
    return complex(b1 / a1)


@pytest.fixture(scope="module")
def tune_esr():
    from tests.validation import test_th14_birdcage_copper as th14

    r_ohm = _esr_ohm()
    if r_ohm is None:
        pytest.skip(f"{ESR_ENV} unset: T26 runs under the ESR knob")
    comm = MPI.COMM_WORLD
    deg = _degree()
    act = _active()
    tgt, tune = act["target_hz"], act["tune"]
    jt = _tune_json()
    if not jt.exists():
        pytest.fail(f"{jt} absent")
    tj = json.loads(jt.read_text())
    if tj["degree"] != deg or abs(tj["f_hz"] - tgt) > 1.0:
        pytest.fail(f"{jt}: degree {tj['degree']} / f_hz {tj['f_hz']} != {deg} / {tgt}")
    # §9 item 395 (T37): the h = 0.01 rung at 64 MHz (``_active`` / ``_tune_json`` /
    # ``_stored_active`` already return the T28 set; the name below is unchanged on
    # the x1 and 127.74 paths, whose ``act`` carries no ``sfx``).
    h37 = _h_tune()
    if h37 and tj.get("resolution") != H_RES:
        pytest.fail(f"{jt}: resolution {tj.get('resolution')} != {H_RES}")
    jl = act["out"] / f"tune_inmodel_cin_p{deg}{act.get('sfx', '')}.json"
    lossless = json.loads(jl.read_text()) if jl.exists() else None
    if h37 and lossless is None:
        pytest.fail(f"{jl} absent: the rung's lossless Q_0 comparand (item 310 W5)")
    c_mode = float(tj["c_mode"]["cu"])
    c_in0 = float(tj["c_in"])
    z64_cu = _uj(tj["z64"]["cu"])
    omega_t = 2.0 * np.pi * tgt
    # C_in(r), pure: secant from the lossless C_in on the stored z64.
    c_in_r, g_r, it_r = _c_in_esr(z64_cu, tgt, c_in0, r_ohm)
    c_rel = c_in_r / c_in0 - 1.0
    # §9 item 387 (T35): the 127.74 MHz path's tag, frequency text and predictions;
    # on the 64 MHz path every print string is byte-identical to item 311's.
    # §9 item 419 (T38): the h = 0.01 rung at 127.74 MHz takes its own branch (h38) first
    # in every print chain; t35 stays the x1 127.74 path (h37 False there).
    h38 = h37 and _tune_mhz() == 127.74
    t35 = _tune_mhz() == 127.74 and not h37
    tag = "T38" if h38 else ("T35" if t35 else ("T37" if h37 else "T26"))
    if h38:
        hm = c_in_r / T35_C_IN_R_X1 - 1.0
        _say(
            f"[ans7-exact] T38 C_in(r = {r_ohm:g} ohm, h {H_RES:g}) {c_in_r!r} F (secant from "
            f"C_in(0, h) {c_in0!r} F read from {jt.parent.name}/{jt.name}; g {g_r:+.3e} S, "
            f"{it_r} iterations); C_in(r, h)/C_in(0, h) - 1 = {c_rel:+.6e} (PREDICTED in "
            f"[+0.8e-5, +2.0e-5], central +1.351e-5: {_met(0.8e-5 <= c_rel <= 2.0e-5)}; x1 "
            f"127.74 record +1.344645e-05 (item 387 W2); printed)"
        )
        _say(
            f"[ans7-exact] T38 C_in(r, h)/C_in(r, x1) - 1 = {hm:+.6e} vs T35_C_IN_R_X1 "
            f"{T35_C_IN_R_X1!r} F (PREDICTED in [+8.6e-4, +9.3e-4], central +8.9685e-4: "
            f"{_met(8.6e-4 <= hm <= 9.3e-4)}); lossless h-move T34_C_IN_H_REL_128 "
            f"{T34_C_IN_H_REL_128:+.6e}; difference {hm - T34_C_IN_H_REL_128:+.3e}; (N-h') "
            f"|.| >= 10 x secant tol (1e-11): {_met(abs(hm) >= 1e-11)} (PREDICTED met; "
            f"printed, never asserted)"
        )
    elif h37:
        hm = c_in_r / T26_C_IN_R_X1 - 1.0
        _say(
            f"[ans7-exact] T37 C_in(r = {r_ohm:g} ohm, h {H_RES:g}) {c_in_r!r} F (secant from "
            f"C_in(0, h) {c_in0!r} F read from {jt.name}; g {g_r:+.3e} S, {it_r} iterations); "
            f"C_in(r, h)/C_in(0, h) - 1 = {c_rel:+.6e} (PREDICTED in [+3e-5, +6e-5], central "
            f"+4.610e-5: {_met(3e-5 <= c_rel <= 6e-5)}; x1 record {T26_C_IN_R_REL_64:+.6e}; "
            f"printed)"
        )
        _say(
            f"[ans7-exact] T37 C_in(r, h)/C_in(r, x1) - 1 = {hm:+.6e} vs T26_C_IN_R_X1 "
            f"{T26_C_IN_R_X1!r} F (PREDICTED in [+8.0e-4, +8.7e-4], central +8.3473e-4: "
            f"{_met(8.0e-4 <= hm <= 8.7e-4)}); lossless h-move T28_C_IN_H_REL_64 "
            f"{T28_C_IN_H_REL_64:+.6e}; (N-h) |.| >= 10 x secant tol (1e-11): "
            f"{_met(abs(hm) >= 1e-11)} (PREDICTED met; printed, never asserted)"
        )
    elif t35:
        _say(
            f"[ans7-exact] T35 C_in(r = {r_ohm:g} ohm) {c_in_r!r} F (secant from C_in(0) "
            f"{c_in0!r} F read from {jt.parent.name}/{jt.name}; g {g_r:+.3e} S, {it_r} "
            f"iterations); C_in(r)/C_in(0) - 1 = {c_rel:+.6e} (PREDICTED in [+2e-6, +5e-5], "
            f"central ~+1.06e-5: {_met(2e-6 <= c_rel <= 5e-5)}; printed); 64 MHz record "
            f"T26_C_IN_R_REL_64 {T26_C_IN_R_REL_64:+.6e} (item 311 W2; ratio 127.74/64 = "
            f"{c_rel / T26_C_IN_R_REL_64:.4f}; printed, never asserted)"
        )
    else:
        _say(
            f"[ans7-exact] T26 C_in(r = {r_ohm:g} ohm) {c_in_r!r} F (secant from C_in(0) {c_in0!r} F "
            f"read from {jt.name}; g {g_r:+.3e} S, {it_r} iterations); C_in(r)/C_in(0) - 1 = "
            f"{c_rel:+.6e} (PREDICTED in [-1e-3, +1e-3]: {_met(abs(c_rel) <= 1e-3)}; printed)"
        )
    # (N-r0), predicted, printed: r = 0 through the ESR path from _c_in's start (C_mode).
    c_r0, _, _ = _c_in_esr(z64_cu, tgt, c_mode, 0.0)
    _say(
        f"[ans7-exact] {tag} (N-r0) ESR-path root at r = 0 from C_mode: {c_r0!r} F vs "
        f"tune_p2.json c_in {c_in0!r} F; bit-equal {c_r0 == c_in0}, rel "
        f"{c_r0 / c_in0 - 1.0:+.3e} (PREDICTED bit-equal: {_met(c_r0 == c_in0)}; printed)"
    )
    recs, _ = _stored_active()
    z_nodes_cu = [_uj(rc["z"]["cu"]) for rc in recs]
    xn = np.asarray(act["nodes"], dtype=float) - 64.5
    if h37:
        # §9 item 395: (F-Q0) -- the lossless comparand is the rung's, not x1's.
        q0_ok = lossless["c_in"] == tj["c_in"] and lossless["cells"] == recs[0]["cells"]
        _say(
            f"[ans7-exact] {tag} (F-Q0) {jl.name}: c_in {lossless['c_in']!r} vs {jt.name} "
            f"{tj['c_in']!r}; cells {lossless['cells']} vs node cells {recs[0]['cells']} "
            f"(ASSERTED equal: {_met(q0_ok)})"
        )
        # The rung's tie from tune_p2_h0.01.json (no x1 mesh, no x1 record assert).
        tie = _tie_from_json(tj["tie"])
    else:
        tie = _fixture_tie()
    t0 = time.perf_counter()
    if h37:
        # (W) first, on the x1 mesh, then free it before meshing the rung (tune_h).
        geo1, specs1, _c1 = _exact_geo(th14, comm, GATED_WIDTH_FRACTION)
        t1 = time.perf_counter()
        s_w, _ = _route_a_solve(th14, geo1, specs1, TUNE_TARGET_HZ, "cu")
        del geo1, specs1
    else:
        geo, specs, cells = _exact_geo(th14, comm, GATED_WIDTH_FRACTION)
        t_mesh = time.perf_counter() - t0
        # (W) unchanged: item 262's degree-1 default drive against item 195's stored point.
        t1 = time.perf_counter()
        s_w, _ = _route_a_solve(th14, geo, specs, TUNE_TARGET_HZ, "cu")
    rb = json.loads(_route_b_json(1).read_text())
    ra = json.loads((OUT_DIR / "route_a_p1.json").read_text())
    xn1 = np.asarray(rb["nodes_mhz"], dtype=float) - 64.5
    zb1 = route_b_z(
        bary_interp(xn1, [_uj(m) for m in rb["z"]["cu"]], 64.0 - 64.5), TUNE_TARGET_HZ, C_F
    )[1]
    s_b1 = _s11_from_z(zb1)
    res_w = float(abs(s_w - s_b1) / abs(s_b1))
    j64 = list(np.round(GRID_MHZ, 6)).index(64.0)
    res_195 = float(ra["residual"]["cu"][j64])
    dw = abs(res_w - res_195)
    t_w = time.perf_counter() - t1
    _say(
        f"[ans7-exact] {tag} (W) degree-1 default drive 64.0 MHz Cu: residual {res_w!r}; item 195 "
        f"stored {res_195!r}; |delta| = {dw:.3e} (ASSERTED <= NODE_S_RTOL {NODE_S_RTOL:.0e}); "
        f"{t_w:.1f} s"
    )
    if h37:
        t4 = time.perf_counter()
        geo, specs, cells = _exact_geo(th14, comm, GATED_WIDTH_FRACTION, res=H_RES)
        t_mesh = time.perf_counter() - t4
        _say(
            f"[ans7-exact] {tag} rung mesh h {H_RES:g}: cells {cells} (node cells "
            f"{recs[0]['cells']}, {jt.name} {tj['cells']}); {t_mesh:.1f} s"
        )
    res, s_as, s_bs, z_in_a, z_in_b, t_drive = [], [], [], [], [], []
    for g in tune:
        f = g * 1e6
        omega = 2.0 * np.pi * f
        t2 = time.perf_counter()
        s_a = _route_a_solve_esr(th14, geo, specs, f, c_in_r, r_ohm, deg)
        t_drive.append(time.perf_counter() - t2)
        zin_b, zbare_b, _ = _route_b_z_esr(bary_interp(xn, z_nodes_cu, g - 64.5), f, c_in_r, r_ohm)
        s_b = _s11_from_z(zbare_b)
        rr = float(abs(s_a - s_b) / abs(s_b))
        res.append(rr)
        s_as.append(complex(s_a))
        s_bs.append(complex(s_b))
        zbr = _z_branch(omega, c_in_r, r_ohm)
        z_in_a.append(1.0 / (1.0 / _ANS7._z_from_s11(s_a) + 1.0 / zbr))
        z_in_b.append(complex(zin_b))
        if h38:
            _say(
                f"[ans7-exact] T38 (AB-rh') p{deg} Cu {g:.2f} MHz at (r, C_in(r, h)), h "
                f"{H_RES:g}: S11_A {complex(s_a):.12e}, S11_B {complex(s_b):.12e}; |S11_A - "
                f"S11_B|/|S11_B| = {rr:.3e} (ASSERTED <= REDUCTION_BAND {REDUCTION_BAND:.0e}; "
                f"PREDICTED [8e-6, 4e-5]: {_met(8e-6 <= rr <= 4e-5)}; lossless h 1.593 / "
                f"1.605 / 1.619e-5, x1 ESR 1.587 / 1.600 / 1.614e-5); drive "
                f"{t_drive[-1]:.1f} s"
            )
        elif t35:
            _say(
                f"[ans7-exact] T35 (AB-r') p{deg} Cu {g:.2f} MHz at (r, C_in(r)): S11_A "
                f"{complex(s_a):.12e}, S11_B {complex(s_b):.12e}; |S11_A - S11_B|/|S11_B| = "
                f"{rr:.3e} (ASSERTED <= REDUCTION_BAND {REDUCTION_BAND:.0e}; PREDICTED "
                f"[1.0e-5, 3e-5]: {_met(1.0e-5 <= rr <= 3e-5)}; lossless 1.590 / 1.601 / "
                f"1.616e-5); drive {t_drive[-1]:.1f} s"
            )
        elif h37:
            _say(
                f"[ans7-exact] T37 (AB-rh) p{deg} Cu {g:.1f} MHz at (r, C_in(r, h)), h {H_RES:g}: "
                f"S11_A {complex(s_a):.12e}, S11_B {complex(s_b):.12e}; |S11_A - S11_B|/|S11_B| = "
                f"{rr:.3e} (ASSERTED <= REDUCTION_BAND {REDUCTION_BAND:.0e}; PREDICTED "
                f"[5e-6, 3e-5]: {_met(5e-6 <= rr <= 3e-5)}; lossless h 1.239 / 1.262 / "
                f"1.289e-5); drive {t_drive[-1]:.1f} s"
            )
        else:
            _say(
                f"[ans7-exact] T26 (AB-r) p{deg} Cu {g:.1f} MHz at (r, C_in(r)): S11_A "
                f"{complex(s_a):.12e}, S11_B {complex(s_b):.12e}; |S11_A - S11_B|/|S11_B| = "
                f"{rr:.3e} (ASSERTED <= {REDUCTION_BAND:.0e}; PREDICTED [5e-6, 3e-5]: "
                f"{_met(5e-6 <= rr <= 3e-5)}); drive {t_drive[-1]:.1f} s"
            )
    rows_a = parallel_rows(np.array(tune) * 1e6, np.array(z_in_a))
    rows_b = parallel_rows(np.array(tune) * 1e6, np.array(z_in_b))
    b = (1.0 / np.asarray(z_in_a)).imag
    n_up = sum(1 for k in range(len(b) - 1) if b[k] < 0 <= b[k + 1])
    n_dn = sum(1 for k in range(len(b) - 1) if b[k] >= 0 > b[k + 1])
    _say(
        f"[ans7-exact] {tag} (T2-r) Im Y_in,A at {tune} MHz = "
        + ", ".join(f"{v:+.9e}" for v in b)
        + f" S; rows A: f0 {rows_a.get('f0_hz')!r} Hz, R_in "
        f"{rows_a.get('r_in_ohm')!r} ohm, Q {rows_a.get('q_loaded')!r}; up-crossings "
        f"{rows_a['crossings_up_hz']}"
    )
    f0a, f0b = rows_a.get("f0_hz"), rows_b.get("f0_hz")
    sa = (f0a / tgt - 1.0) if f0a else float("nan")
    sb = (f0b / tgt - 1.0) if f0b else float("nan")
    r_ab = (f0a / f0b - 1.0) if (f0a and f0b) else float("nan")
    if h38:
        zr = T22P_ZERO_RESOLUTION
        _say(
            f"[ans7-exact] T38 (T2-rh') crossings up {rows_a['crossings_up_hz']}; Im Y_in,A sign "
            f"changes - -> + {n_up}, + -> - {n_dn} (ASSERTED exactly one up, none down, f0,A in "
            f"(127.64, 127.84) MHz); f0,A(C_in(r, h))/127.74 MHz - 1 = {sa:+.6e} (PRINTED for "
            f"scale only beside T22P_ZERO_RESOLUTION {zr:.3e}, = {abs(sa) / zr:.3f} x, never "
            f"asserted: rulings 310-Q1 / 365-Q1 / 395-Q3; PREDICTED ~2.29e-6, ~0.142 x, in "
            f"[1e-6, 5e-6]: {_met(1e-6 <= abs(sa) <= 5e-6)}; lossless h 2.271224e-6, x1 ESR "
            f"2.301805e-6)"
        )
        _say(
            f"[ans7-exact] T38 f0,B(C_in(r, h)) {f0b!r} Hz (Route B on the three interpolated "
            f"tune points); f0,B/127.74 MHz - 1 = {sb:+.6e} (PREDICTED within +-1e-9: "
            f"{_met(abs(sb) <= 1e-9)}); r_AB = f0,A/f0,B - 1 = {r_ab:+.6e} (PREDICTED "
            f"[-5e-6, -1e-6]: {_met(-5e-6 <= r_ab <= -1e-6)}; lossless h -2.271224e-6); R_in,B "
            f"{rows_b.get('r_in_ohm')!r} ohm, Q_B {rows_b.get('q_loaded')!r}"
        )
    elif t35:
        zr = _zero_resolution()
        _say(
            f"[ans7-exact] T35 (T2-r') f0,A(C_in(r))/127.74 MHz - 1 = {sa:+.6e} (ASSERTED |.| "
            f"<= T22P_ZERO_RESOLUTION {zr:.3e}: {_met(abs(sa) <= zr)}; = {abs(sa) / zr:.3f} "
            f"bands; PREDICTED ~2.3e-6, ~0.14 bands: {_met(1e-6 <= abs(sa) <= 5e-6)})"
        )
        _say(
            f"[ans7-exact] T35 f0,B(C_in(r)) {f0b!r} Hz (Route B on the three interpolated tune "
            f"points); f0,B/127.74 MHz - 1 = {sb:+.6e} (PREDICTED within +-1e-9: "
            f"{_met(abs(sb) <= 1e-9)}); r_AB = f0,A/f0,B - 1 = {r_ab:+.6e} (PREDICTED "
            f"[-5e-6, -1e-6]: {_met(-5e-6 <= r_ab <= -1e-6)}; lossless -2.265e-6); R_in,B "
            f"{rows_b.get('r_in_ohm')!r} ohm, Q_B {rows_b.get('q_loaded')!r}"
        )
    elif h37:
        zr = T21P_ZERO_RESOLUTION
        _say(
            f"[ans7-exact] T37 (T2-rh) crossings up {rows_a['crossings_up_hz']}; Im Y_in,A sign "
            f"changes - -> + {n_up}, + -> - {n_dn} (ASSERTED exactly one up, none down, f0,A in "
            f"(63.9, 64.1) MHz); f0,A(C_in(r, h))/64 MHz - 1 = {sa:+.6e} (PRINTED for scale "
            f"only beside T21P_ZERO_RESOLUTION {zr:.3e}, = {abs(sa) / zr:.3f} x, never "
            f"asserted: ruling 310-Q1; PREDICTED ~2.59e-6, ~0.20 x: "
            f"{_met(1e-6 <= abs(sa) <= 5e-6)})"
        )
        _say(
            f"[ans7-exact] T37 f0,B(C_in(r, h)) {f0b!r} Hz (Route B on the three interpolated "
            f"tune points); f0,B/64 MHz - 1 = {sb:+.6e} (PREDICTED within +-1e-9: "
            f"{_met(abs(sb) <= 1e-9)}); r_AB = f0,A/f0,B - 1 = {r_ab:+.6e} (PREDICTED "
            f"[-5e-6, -1e-6]: {_met(-5e-6 <= r_ab <= -1e-6)}; lossless h -2.545621e-6); R_in,B "
            f"{rows_b.get('r_in_ohm')!r} ohm, Q_B {rows_b.get('q_loaded')!r}"
        )
    else:
        _say(
            f"[ans7-exact] T26 (T2-r) f0,A(C_in(r))/64 MHz - 1 = {sa:+.6e} (ASSERTED |.| <= "
            f"T21P_ZERO_RESOLUTION {T21P_ZERO_RESOLUTION:.3e}: "
            f"{_met(abs(sa) <= T21P_ZERO_RESOLUTION)}; = {abs(sa) / T21P_ZERO_RESOLUTION:.3f} bands; "
            f"PREDICTED ~2.5e-6)"
        )
        _say(
            f"[ans7-exact] T26 f0,B(C_in(r)) {f0b!r} Hz (Route B on the three interpolated tune "
            f"points); f0,B/64 MHz - 1 = {sb:+.6e} (PREDICTED within +-1e-9: "
            f"{_met(abs(sb) <= 1e-9)}); r_AB = f0,A/f0,B - 1 = {r_ab:+.6e} (PREDICTED [-5e-6, -1e-6]: "
            f"{_met(-5e-6 <= r_ab <= -1e-6)}; lossless -2.536e-6); R_in,B {rows_b.get('r_in_ohm')!r} "
            f"ohm, Q_B {rows_b.get('q_loaded')!r}"
        )
    # Q(r)/Q(0) beside the closed form 1/Q = 1/Q_0 + 1/Q_C, Q_C = 1/(omega C_in r).
    q0 = lossless["rows_a"].get("q_loaded") if lossless else None
    r0 = lossless["rows_a"].get("r_in_ohm") if lossless else None
    q_c = 1.0 / (omega_t * c_in_r * r_ohm)
    qr = rows_a.get("q_loaded")
    if q0 and qr:
        q_cf = 1.0 / (1.0 / q0 + 1.0 / q_c)
        if h38:
            ex = qr / q_cf - 1.0
            q0m = q0 / T22P_Q0_X1 - 1.0
            qm = qr / T35_Q_R_X1 - 1.0
            _say(
                f"[ans7-exact] T38 Q(r, h) {qr!r} vs Q(0, h) = Q_0(h) {q0!r} ({jl.parent.name}/"
                f"{jl.name}, item 365b A): Q(r, h)/Q_0(h) = {qr / q0:.6f} (PREDICTED in "
                f"[0.405, 0.425]: {_met(0.405 <= qr / q0 <= 0.425)}); closed form Q_C(h) = "
                f"1/(omega C_in r) = {q_c:.4f} (PREDICTED ~407.47, x1 407.83), 1/Q = 1/Q_0 + "
                f"1/Q_C -> Q {q_cf:.4f}, ratio {q_cf / q0:.6f} (PREDICTED ~243.756, ~0.40178); "
                f"Q(r, h) PREDICTED in [246, 258], central 252.41: {_met(246.0 <= qr <= 258.0)}; "
                f"Q(r, h)/closed form - 1 = {ex:+.4%} (PREDICTED in [+2.5 %, +4.5 %]: "
                f"{_met(0.025 <= ex <= 0.045)}; x1 T35_Q_CLOSED_EXCESS_X1 "
                f"{T35_Q_CLOSED_EXCESS_X1:+.4%}; printed, never asserted)"
            )
            _say(
                f"[ans7-exact] T38 step reading: Q(r, h)/Q(r, x1) - 1 = {qm:+.6e} vs "
                f"T35_Q_R_X1 {T35_Q_R_X1!r} (PREDICTED in [-1.2e-2, -4e-3], central -8.10e-3: "
                f"{_met(-1.2e-2 <= qm <= -4e-3)}); Q_0's h-move at C_in Q_0(h)/Q_0(x1) - 1 = "
                f"{q0m:+.6e} (T22P_Q0_X1 {T22P_Q0_X1!r}; brief -1.8827e-2); ratio Q(r) move / "
                f"Q_0 move = {qm / q0m:.4f} (PREDICTED ~0.43, in [0.2, 0.7]: "
                f"{_met(0.2 <= qm / q0m <= 0.7)}; 64 MHz rung 0.2915, item 395); |Q(r) move| < "
                f"|Q_0 move|: {_met(abs(qm) < abs(q0m))} (printed, never asserted)"
            )
        elif t35:
            _say(
                f"[ans7-exact] T35 Q(r) {qr!r} vs Q(0) = Q_0 {q0!r} ({jl.parent.name}/{jl.name}, "
                f"item 283): Q(r)/Q(0) = {qr / q0:.6f}; closed form Q_C = 1/(omega C_in r) = "
                f"{q_c:.4f}, 1/Q = 1/Q_0 + 1/Q_C -> Q {q_cf:.4f}, ratio {q_cf / q0:.6f} "
                f"(PREDICTED Q(r) in [240, 255]: {_met(240.0 <= qr <= 255.0)}; ratio in "
                f"[0.38, 0.42]: {_met(0.38 <= qr / q0 <= 0.42)}); Q(r)/closed form - 1 = "
                f"{qr / q_cf - 1.0:+.4%} (PREDICTED in [0, +3 %]: "
                f"{_met(0.0 <= qr / q_cf - 1.0 <= 0.03)}; printed, never asserted)"
            )
        elif h37:
            ex = qr / q_cf - 1.0
            q0m = q0 / T21P_Q0_X1 - 1.0
            qm = qr / T26_Q_R_X1 - 1.0
            _say(
                f"[ans7-exact] T37 Q(r, h) {qr!r} vs Q(0, h) = Q_0(h) {q0!r} ({jl.name}, item "
                f"310 W5): Q(r, h)/Q_0(h) = {qr / q0:.6f} (PREDICTED in [0.166, 0.174]: "
                f"{_met(0.166 <= qr / q0 <= 0.174)}); closed form Q_C(h) = 1/(omega C_in r) = "
                f"{q_c:.4f} (PREDICTED ~195.36), 1/Q = 1/Q_0 + 1/Q_C -> Q {q_cf:.4f}, ratio "
                f"{q_cf / q0:.6f} (PREDICTED ~162.454, ~0.16844); Q(r, h) PREDICTED in [161, 167], "
                f"central 164.24: {_met(161.0 <= qr <= 167.0)}; Q(r, h)/closed form - 1 = "
                f"{ex:+.4%} (PREDICTED in [+0.7 %, +1.5 %]: {_met(0.007 <= ex <= 0.015)}; x1 "
                f"T26_Q_CLOSED_EXCESS_X1 {T26_Q_CLOSED_EXCESS_X1:+.4%}; printed, never asserted)"
            )
            _say(
                f"[ans7-exact] T37 step reading: Q(r, h)/Q(r, x1) - 1 = {qm:+.6e} vs "
                f"T26_Q_R_X1 {T26_Q_R_X1!r} (PREDICTED in [-4e-3, 0], central -1.84e-3: "
                f"{_met(-4e-3 <= qm <= 0.0)}); Q_0's h-move at C_in Q_0(h)/Q_0(x1) - 1 = "
                f"{q0m:+.6e} (T21P_Q0_X1 {T21P_Q0_X1!r}; brief -6.832e-3); T30's move at C_F "
                f"T30_Q_H_REL_CF {T30_Q_H_REL_CF:+.4e}; ratio Q(r) move / Q_0 move = "
                f"{qm / q0m:.4f} (PREDICTED ~0.27, in [0, 0.6]: {_met(0.0 <= qm / q0m <= 0.6)}); "
                f"|Q(r) move| < |Q_0 move|: {_met(abs(qm) < abs(q0m))} (printed, never asserted)"
            )
        else:
            _say(
                f"[ans7-exact] T26 Q(r) {qr!r} vs Q(0) = Q_0 {q0!r} ({jl.name}, item 262): "
                f"Q(r)/Q(0) = {qr / q0:.6f}; closed form Q_C = 1/(omega C_in r) = {q_c:.4f}, "
                f"1/Q = 1/Q_0 + 1/Q_C -> Q {q_cf:.4f}, ratio {q_cf / q0:.6f} (PREDICTED Q(r) ~163, "
                f"ratio ~0.168: {_met(0.1 <= qr / q0 <= 0.25)}); Q(r)/closed form - 1 = "
                f"{qr / q_cf - 1.0:+.4%} (printed, never asserted)"
            )
    if r0 and rows_a.get("r_in_ohm"):
        if h38:
            rh = rows_a["r_in_ohm"] / T35_R_IN_R_X1 - 1.0
            _say(
                f"[ans7-exact] T38 R_in(r, h) {rows_a['r_in_ohm']!r} ohm (PREDICTED ~9703) vs "
                f"R_in(0, h) {r0!r} ohm: ratio {rows_a['r_in_ohm'] / r0:.6f}; R_in(r, h)/R_in(r, "
                f"x1) - 1 = {rh:+.6e} vs T35_R_IN_R_X1 {T35_R_IN_R_X1!r} (PREDICTED in "
                f"[-1.4e-2, -4e-3], central -9.0e-3: {_met(-1.4e-2 <= rh <= -4e-3)}; printed)"
            )
        elif t35:
            rq = rows_a["r_in_ohm"] / r0
            _say(
                f"[ans7-exact] T35 R_in(r) {rows_a['r_in_ohm']!r} ohm vs R_in(0) {r0!r} ohm: "
                f"ratio {rq:.6f} (PREDICTED in [0.37, 0.43]: {_met(0.37 <= rq <= 0.43)}; printed)"
            )
        elif h37:
            rh = rows_a["r_in_ohm"] / T26_R_IN_R_X1 - 1.0
            _say(
                f"[ans7-exact] T37 R_in(r, h) {rows_a['r_in_ohm']!r} ohm (PREDICTED ~3166) vs "
                f"R_in(0, h) {r0!r} ohm: ratio {rows_a['r_in_ohm'] / r0:.6f}; R_in(r, h)/R_in(r, "
                f"x1) - 1 = {rh:+.6e} vs T26_R_IN_R_X1 {T26_R_IN_R_X1!r} (PREDICTED in [-6e-3, "
                f"-1e-3]: {_met(-6e-3 <= rh <= -1e-3)}; printed)"
            )
        else:
            _say(
                f"[ans7-exact] T26 R_in(r) {rows_a['r_in_ohm']!r} ohm vs R_in(0) {r0!r} ohm: ratio "
                f"{rows_a['r_in_ohm'] / r0:.6f} (printed)"
            )
    j = [round(g, 6) for g in tune].index(round(tgt / 1e6, 6))
    # (N-L), predicted, printed: lossless circuit (r = 0 at C_in(0)) vs in-model at (r, C_in(r)).
    s_nl = _s11_from_z(route_b_z(z64_cu, tgt, c_in0)[1])
    n_l = float(abs(s_nl - s_as[j]) / abs(s_as[j]))
    zin_nl = complex(route_b_z(z64_cu, tgt, c_in0)[0])
    n_l_zin = float(abs(zin_nl - z_in_a[j]) / abs(z_in_a[j]))
    if h38:
        _say(
            f"[ans7-exact] T38 (N-Lh') lossless circuit at C_in(0, h) vs in-model at (r, "
            f"C_in(r, h)), 127.74 MHz Cu: |dS11_bare|/|S11| = {n_l:.6e} = "
            f"{n_l / REDUCTION_BAND:.2f} x REDUCTION_BAND (PREDICTED in [1, 4] x, i.e. NOT "
            f"separating at 10 x: {_met(1.0 * REDUCTION_BAND <= n_l <= 4.0 * REDUCTION_BAND)}; "
            f">= 10 x: {_met(n_l >= 10.0 * REDUCTION_BAND)}; x1 2.17 x; printed, never "
            f"asserted); on Z_in |dZ|/|Z| = {n_l_zin:.6e} (PREDICTED ~Q_0(h)/Q_closed - 1 "
            f"~1.489, in [1.3, 1.7]: {_met(1.3 <= n_l_zin <= 1.7)}; x1 1.419; printed)"
        )
    elif t35:
        _say(
            f"[ans7-exact] T35 (N-L') lossless circuit at C_in(0) vs in-model at (r, C_in(r)), "
            f"127.74 MHz Cu: |dS11_bare|/|S11| = {n_l:.6e} = {n_l / REDUCTION_BAND:.2f} x "
            f"REDUCTION_BAND (PREDICTED in [0.5, 5] x, i.e. NOT separating at 10 x: "
            f"{_met(0.5 * REDUCTION_BAND <= n_l <= 5.0 * REDUCTION_BAND)}; >= 10 x: "
            f"{_met(n_l >= 10.0 * REDUCTION_BAND)}; printed, never asserted); on Z_in |dZ|/|Z| "
            f"= {n_l_zin:.6e} (PREDICTED ~Q_0/Q_closed - 1 ~1.52, in [1.3, 1.8]: "
            f"{_met(1.3 <= n_l_zin <= 1.8)}; printed)"
        )
    elif h37:
        _say(
            f"[ans7-exact] T37 (N-Lh) lossless circuit at C_in(0, h) vs in-model at (r, C_in(r, h)), "
            f"64.0 MHz Cu: |dS11_bare|/|S11| = {n_l:.6e} = {n_l / REDUCTION_BAND:.2f} x "
            f"REDUCTION_BAND (PREDICTED in [3, 7] x, i.e. NOT separating at 10 x: "
            f"{_met(3.0 * REDUCTION_BAND <= n_l <= 7.0 * REDUCTION_BAND)}; >= 10 x: "
            f"{_met(n_l >= 10.0 * REDUCTION_BAND)}; x1 4.97 x; printed, never asserted); on Z_in "
            f"|dZ|/|Z| = {n_l_zin:.6e} (PREDICTED ~4.94, in [4.6, 5.2]: "
            f"{_met(4.6 <= n_l_zin <= 5.2)}; x1 4.878; printed)"
        )
    else:
        _say(
            f"[ans7-exact] T26 (N-L) lossless circuit at C_in(0) vs in-model at (r, C_in(r)), "
            f"64.0 MHz Cu: |dS11_bare|/|S11| = {n_l:.6e} = {n_l / REDUCTION_BAND:.2f} x "
            f"REDUCTION_BAND (PREDICTED >= 10 x: {_met(n_l >= 10.0 * REDUCTION_BAND)}; printed, "
            f"never asserted); on Z_in |dZ|/|Z| = {n_l_zin:.6e}"
        )
    # (N2-r), predicted, printed: circuit at 0.5*C_in(r) (with r) vs in-model at C_in(r).
    s_n2 = _s11_from_z(_route_b_z_esr(z64_cu, tgt, 0.5 * c_in_r, r_ohm)[1])
    n2 = float(abs(s_n2 - s_as[j]) / abs(s_as[j]))
    if h38:
        _say(
            f"[ans7-exact] T38 (N2-rh') circuit at 0.5*C_in(r, h) vs in-model at C_in(r, h), "
            f"127.74 MHz Cu: |dS|/|S| = {n2:.6e} = {n2 / REDUCTION_BAND:.2f} x REDUCTION_BAND "
            f"(PREDICTED >= 10 x: {_met(n2 >= 10.0 * REDUCTION_BAND)}; lossless h 494.00 x, x1 "
            f"ESR 494.07 x; printed, never asserted)"
        )
    elif t35:
        _say(
            f"[ans7-exact] T35 (N2-r') circuit at 0.5*C_in(r) vs in-model at C_in(r), 127.74 MHz "
            f"Cu: |dS|/|S| = {n2:.6e} = {n2 / REDUCTION_BAND:.2f} x REDUCTION_BAND (PREDICTED "
            f">= 10 x: {_met(n2 >= 10.0 * REDUCTION_BAND)}; lossless 493.82 x; printed, never "
            f"asserted)"
        )
    elif h37:
        _say(
            f"[ans7-exact] T37 (N2-rh) circuit at 0.5*C_in(r, h) vs in-model at C_in(r, h), 64.0 "
            f"MHz Cu: |dS|/|S| = {n2:.6e} = {n2 / REDUCTION_BAND:.2f} x REDUCTION_BAND (PREDICTED "
            f">= 10 x: {_met(n2 >= 10.0 * REDUCTION_BAND)}; lossless h 493.62 x, x1 ESR 494.39 x; "
            f"printed, never asserted)"
        )
    else:
        _say(
            f"[ans7-exact] T26 (N2-r) circuit at 0.5*C_in(r) vs in-model at C_in(r), 64.0 MHz Cu: "
            f"|dS|/|S| = {n2:.6e} = {n2 / REDUCTION_BAND:.2f} x REDUCTION_BAND (PREDICTED >= 10 x: "
            f"{_met(n2 >= 10.0 * REDUCTION_BAND)}; lossless 493.78 x; printed, never asserted)"
        )
    _say(
        f"[ans7-exact] {tag} W2 total {time.perf_counter() - t0:.1f} s (mesh {t_mesh:.1f} s, (W) "
        f"{t_w:.1f} s, degree-{deg} drives {', '.join(f'{t:.1f}' for t in t_drive)} s), "
        f"cells {cells}, -n {comm.size}"
    )
    if comm.rank == 0:
        act["out"].mkdir(parents=True, exist_ok=True)
        p_tune = {
            "c_in_r": c_in_r,
            "c_in_0": c_in0,
            "esr_ohm": r_ohm,
            "f_hz": tgt,
            "g": g_r,
            "iterations": it_r,
            "c_r0": c_r0,
        }
        p_inm = {
            "degree": deg,
            "ranks": comm.size,
            "cells": cells,
            "esr_ohm": r_ohm,
            "c_in_r": c_in_r,
            "c_in_0": c_in0,
            "tune_mhz": list(tune),
            "residual_ab": res,
            "s11_a": [[s.real, s.imag] for s in s_as],
            "s11_b": [[s.real, s.imag] for s in s_bs],
            "z_in_a": [[z.real, z.imag] for z in z_in_a],
            "rows_a": rows_a,
            "rows_b": rows_b,
            "w_residual": res_w,
            "w_item195": res_195,
            "n_l": n_l,
            "n2": n2,
            "q_c": q_c,
            "t_drive_s": t_drive,
        }
        if h37:
            # §9 item 395: the rung's own names; item 311's x1 ESR files are never written here.
            p_tune["resolution"] = H_RES
            p_inm["resolution"] = H_RES
            (act["out"] / f"tune_p{deg}{H_SFX}_esr{r_ohm:g}.json").write_text(json.dumps(p_tune))
            (act["out"] / f"tune_inmodel_cin_esr{r_ohm:g}_p{deg}{H_SFX}.json").write_text(
                json.dumps(p_inm)
            )
        else:
            (act["out"] / f"tune_p{deg}_esr{r_ohm:g}.json").write_text(json.dumps(p_tune))
            (act["out"] / f"tune_inmodel_cin_esr{r_ohm:g}_p{deg}.json").write_text(
                json.dumps(p_inm)
            )
    comm.Barrier()
    return {
        "cells": cells,
        "tie": tie,
        "res": res,
        "rows_a": rows_a,
        "dw": dw,
        "n_l": n_l,
        "stored_cells": recs[0]["cells"],
        "tune": tune,
        "tgt": tgt,
        "h": h37,
        "hres": H_RES if h37 else None,
        "n_up": n_up,
        "n_dn": n_dn,
        "q0_ok": (
            (lossless["c_in"] == tj["c_in"] and lossless["cells"] == recs[0]["cells"])
            if h37
            else None
        ),
        "tj_cells": tj.get("cells"),
    }


@complex_only
def test_t26_ab_circuit_equals_field_esr(tune_esr):
    """(AB-r) + (W): at (r, C_in(r)), |S11_A − S11_B|/|S11_B| <= ``REDUCTION_BAND``
    at all three ``TUNE_MHZ``; (W) <= ``NODE_S_RTOL``.  On the h = 0.01 rung: T37 at
    64 MHz (§9 item 395) and T38 at 127.74 MHz (§9 item 419)."""
    _assert_tie(tune_esr["cells"], tune_esr["tie"], res=tune_esr["hres"])
    assert tune_esr["dw"] <= NODE_S_RTOL, tune_esr["dw"]
    assert tune_esr["cells"] == tune_esr["stored_cells"]
    if tune_esr["h"]:
        # §9 item 395 (T37): (F) the drive mesh is tune_p2_h0.01.json's; (F-Q0) the
        # lossless Q_0 comparand is the rung's (c_in and cells).
        assert tune_esr["cells"] == tune_esr["tj_cells"], (tune_esr["cells"], tune_esr["tj_cells"])
        assert tune_esr["q0_ok"], "(F-Q0)"
    assert len(tune_esr["res"]) == len(tune_esr["tune"]) == 3
    for g, r in zip(tune_esr["tune"], tune_esr["res"]):
        assert r <= REDUCTION_BAND, (g, r)


@complex_only
def test_t26_inmodel_zero_esr(tune_esr):
    """(T2-r) / (T2-r′): exactly one − → + zero of Im Y_in,A inside the tune
    points at (r, C_in(r)) and |f₀,A/target − 1| <= ``_zero_resolution()``
    (``T21P_ZERO_RESOLUTION`` at 64 MHz, ``T22P_ZERO_RESOLUTION`` at 127.74 MHz;
    imported, chosen by the knob).  (T2-rh), §9 item 395 on the h = 0.01 rung:
    the crossing only -- no band is registered for the rung (ruling 310-Q1), so
    the band line is skipped and the offset is printed by the fixture; and T38 at
    127.74 MHz (§9 item 419)."""
    r = tune_esr["rows_a"]
    assert len(r["crossings_up_hz"]) == 1, r["crossings_up_hz"]
    lo, hi = tune_esr["tune"][0] * 1e6, tune_esr["tune"][-1] * 1e6
    assert lo < r["f0_hz"] < hi, r["f0_hz"]
    if tune_esr["h"]:
        assert tune_esr["n_up"] == 1 and tune_esr["n_dn"] == 0, (tune_esr["n_up"], tune_esr["n_dn"])
        return
    assert abs(r["f0_hz"] / tune_esr["tgt"] - 1.0) <= _zero_resolution(), r["f0_hz"]


@complex_only
def test_t26_n_lossless_separates(tune_esr):
    """(N-L): predicted, printed, never asserted (rule (e)) -- the reading is the
    fixture's print; this test records only that the control was computed."""
    assert np.isfinite(tune_esr["n_l"])
