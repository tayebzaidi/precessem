"""`GEO-34` step 1b-probe — are the imported gap sheets congruent the way the
primitive's are?

Measurement only.  No solve, real build, no `src/` change, no band moved, no
record re-registered, no generator touched.  The two self-checks below are the
only things that can make the script exit non-zero; every symmetry reading is
printed and asserted nowhere.

**Import, never rebuild.**  The primitive rung and its STEP export come from
`GEO-34` step 1a's route —
``tests/mesh/test_birdcage_port_sheets._build(True, step_export_path=...)`` —
and the imported mesh from ``io/step_import.import_step(path, CONFIG_PATH,
comm)`` with ``CONFIG_PATH`` from ``tests/mesh/test_step_import.py`` (the
``round_trip`` fixture body).  The STEP goes to a ``tempfile.mkdtemp`` made on
rank 0 and broadcast, outside the tree.  The measurement is `GEO-31`'s own
``_measure_rung`` / ``_analyse``, loaded from
``scripts/probes/geo31_birdcage_port_sheet_congruence.py`` the way that script
loads ``geo30``, so both tables compare directly with `GEO-31`'s.

**Self-checks (exit non-zero on failure):**

(A) on the primitive, `GEO-31`'s x1 reading reproduced — counts 58/58/58/58,
    areas agreeing <= 1e-3 relative, d_i = 0 / 2.570e-14 / 1.114339e-03 /
    1.114339e-03 m (`20260910T003835Z_GEO-31.log.gz`).  d_1 must be exactly 0,
    d_2 round-off (<= 1e-12 m), d_3 and d_4 within 1e-6 relative of 1.114339e-03.
(B) `GEO-31`'s metric self-test on the *import's* sheet 1: self-distance
    exactly 0, one-l displacement along (1,1,1)/sqrt(3) >= l/2.

**Also printed, only because it is reachable without a solve:** the f = 0.5
narrowed facet counts and w = A/h on both meshes, through the four-port
module's own helpers (``_narrowed_transverse``, ``_sheet_bbox_centre``,
``_sheet_extents``, ``_sheet_axes``, ``_facet_group_area``,
``_sheet_facet_count``, ``GATED_WIDTH_FRACTION``), composed exactly as
``tests/validation/test_port_birdcage_four_port.py:335-381`` composes them.

Run (real build, no complex mode)::

    mpiexec -n 2 python3 scripts/probes/geo34_step1b_import_sheet_congruence.py
"""

from __future__ import annotations

import os
import sys
import tempfile
import time
from pathlib import Path

import numpy as np
from mpi4py import MPI

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))
if str(_REPO_ROOT / "scripts" / "probes") not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT / "scripts" / "probes"))

import geo31_birdcage_port_sheet_congruence as geo31  # noqa: E402

from fem_em_solver.io.step_import import import_step  # noqa: E402
from tests.mesh.test_birdcage_port_sheets import _build  # noqa: E402
from tests.mesh.test_step_import import CONFIG_PATH  # noqa: E402

geo30 = geo31.geo30
PORTS = geo31.PORTS

# `GEO-31`'s x1 reading (`20260910T003835Z_GEO-31.log.gz`; docs/planning/chunks/GEO-31.md).
GEO31_X1_COUNTS = (58, 58, 58, 58)
GEO31_X1_D = (0.0, 2.570e-14, 1.114339e-03, 1.114339e-03)
D_ROUNDOFF = 1.0e-12
D_REL = 1.0e-6

# 1b's import narrowed reading (`20260925T094541Z_GEO-34-step1b.log:3667-3670`).
STEP1B_IMPORT_NARROWED = ((26, 26, 29, 26), (7.264e-3, 7.294e-3, 7.508e-3, 7.294e-3))


def _print(text, comm):
    if comm.rank == 0:
        print(text, flush=True)


def _narrowed(mesh, cell_tags, comm):
    """f = 0.5 narrowed facet counts and w = A/h, or None if not importable."""
    try:
        from tests.mesh.test_birdcage_port_sheets import _sheet_axes
        from tests.mesh.test_two_torus_port_facets import _facet_group_area
        from tests.mesh.test_two_torus_port_sheet import (
            _sheet_extents,
            _sheet_facet_count,
        )
        from tests.validation.test_port_birdcage_lumped_column import (
            _narrowed_transverse,
            _sheet_bbox_centre,
        )
        from tests.validation.test_port_lumped_narrowed_sheet import (
            GATED_WIDTH_FRACTION,
        )
    except Exception as exc:  # noqa: BLE001 — reported, not swallowed
        return {"error": f"{type(exc).__name__}: {exc}"}

    tdim = mesh.topology.dim
    mesh.topology.create_connectivity(tdim - 1, tdim)
    mesh.topology.create_entity_permutations()
    tags_f = geo30._interface_facet_tags(
        mesh,
        cell_tags,
        {
            geo30.SHEET_IFACE + i: (geo30.PORT_LOWER + i, geo30.PORT_UPPER + i)
            for i in PORTS
        },
    )
    geometry = {}
    for i in PORTS:
        tag = geo30.SHEET_IFACE + i
        extents = _sheet_extents(mesh, tags_f, tag, comm)
        w_bbox, _h, _s = _sheet_axes(extents, None, i)
        centroid = _sheet_bbox_centre(mesh, tags_f, tag, comm)
        geometry[i] = {
            "tag": tag,
            "axis": 0 if extents[0] >= extents[1] else 1,
            "centroid": centroid,
            "w_full": float(w_bbox),
        }
    for i in PORTS:
        g = geometry[i]
        tags_f = _narrowed_transverse(
            mesh,
            tags_f,
            g["tag"],
            GATED_WIDTH_FRACTION,
            float(g["centroid"][g["axis"]]),
            g["axis"],
            0.5 * g["w_full"],
        )
    out = {"fraction": GATED_WIDTH_FRACTION, "facets": [], "area": [], "h": [],
           "w": [], "w_full": []}
    for i in PORTS:
        g = geometry[i]
        n = _sheet_facet_count(mesh, tags_f, g["tag"], comm)
        a = _facet_group_area(mesh, tags_f, g["tag"], comm)
        ext = _sheet_extents(mesh, tags_f, g["tag"], comm)
        _wb, h, _s = _sheet_axes(ext, None, i)
        out["facets"].append(int(n))
        out["area"].append(float(a))
        out["h"].append(float(h))
        out["w"].append(float(a / h) if h else float("nan"))
        out["w_full"].append(g["w_full"])
    return out


def _print_table(r, name, size_global, narrowed):
    print(f"\n  mesh {name}  (size_global={size_global}, l={r['ell']:.6e} m)")
    print("    (1) facet count            " + "  ".join(f"{c:13d}" for c in r["counts"])
          + f"   spread={geo31._spread(r['counts']):.4e}")
    print("    (1) gathered facets        "
          + "  ".join(f"{s['n_gathered']:13d}" for s in r["sheets"]))
    print("    (1) total area [m^2]       " + "  ".join(f"{a:.6e}" for a in r["areas"])
          + f"   spread={geo31._spread(r['areas']):.4e}")
    for key in ("min", "mean", "max"):
        vals = [getattr(s["facet_area"], key)() for s in r["sheets"]]
        print(f"    (2) facet area {key:<4s} [m^2]   " + "  ".join(f"{v:.6e}" for v in vals)
              + f"   spread={geo31._spread(vals):.4e}")
    print("    (2) |sorted(Pi)-sorted(P1)|_inf/mean(P1) "
          + "  ".join(f"{v:.6e}" for v in r["sorted_area_delta"])
          + "   (nan = counts differ)")
    print("    (3) d_i centroids [m]      " + "  ".join(f"{v:.6e}" for v in r["d_centroid"]))
    print("    (3) d_i / l                " + "  ".join(f"{v / r['ell']:13.6f}" for v in r["d_centroid"]))
    print("    (4) d_i boundary verts [m] " + "  ".join(f"{v:.6e}" for v in r["d_boundary"]))
    print("    (4) boundary vertex count  " + "  ".join(f"{v:13d}" for v in r["n_boundary"]))
    print("    (5) bbox delta, all verts  " + "  ".join(f"{v:.6e}" for v in r["bbox_delta_all"]))
    print("    (5) bbox delta, boundary   " + "  ".join(f"{v:.6e}" for v in r["bbox_delta_boundary"]))
    if "error" in narrowed:
        print(f"    (6) f=0.5 narrowing NOT REACHABLE: {narrowed['error']}")
    else:
        print(f"    (6) f={narrowed['fraction']} narrowed facets  "
              + "  ".join(f"{v:13d}" for v in narrowed["facets"]))
        print("    (6) narrowed area [m^2]    " + "  ".join(f"{v:.9e}" for v in narrowed["area"]))
        print("    (6) narrowed h [m]         " + "  ".join(f"{v:.9e}" for v in narrowed["h"]))
        print("    (6) narrowed w = A/h [m]   " + "  ".join(f"{v:.9e}" for v in narrowed["w"]))
        print("    (6) full-sheet w_bbox [m]  " + "  ".join(f"{v:.9e}" for v in narrowed["w_full"]))


def main():
    comm = MPI.COMM_WORLD
    started = time.perf_counter()
    _print(
        "\n[GEO-34 1b-probe] gap-sheet C4 congruence, primitive vs STEP import "
        f"(no solve; ranks={comm.size}; measurement is GEO-31's _measure_rung/_analyse)",
        comm,
    )

    step_dir = tempfile.mkdtemp(prefix="geo34p_") if comm.rank == 0 else None
    step_dir = comm.bcast(step_dir, root=0)
    step_path = os.path.join(step_dir, "fsmall_birdcage.step")

    t0 = time.perf_counter()
    g_mesh, g_cells, _gf, _gdiag, g_elapsed = _build(True, step_export_path=step_path)
    t_prim = time.perf_counter() - t0
    g_size = int(g_mesh.topology.index_map(3).size_global)
    _print(f"[GEO-34 1b-probe] primitive built + STEP exported: {g_size} cells, "
           f"{t_prim:.2f} s (STEP -> {step_path})", comm)

    t0 = time.perf_counter()
    i_mesh, i_cells, _if, _idiag = import_step(step_path, CONFIG_PATH, comm=comm)
    t_imp = time.perf_counter() - t0
    i_size = int(i_mesh.topology.index_map(3).size_global)
    _print(f"[GEO-34 1b-probe] STEP imported + meshed: {i_size} cells, {t_imp:.2f} s", comm)

    t0 = time.perf_counter()
    rows = {}
    for name, mesh, cells in (("primitive", g_mesh, g_cells), ("import", i_mesh, i_cells)):
        row = geo31._measure_rung(mesh, cells, comm)
        if comm.rank == 0:
            row = geo31._analyse(row, name)
        rows[name] = row
    narrowed = {
        "primitive": _narrowed(g_mesh, g_cells, comm),
        "import": _narrowed(i_mesh, i_cells, comm),
    }
    t_meas = time.perf_counter() - t0

    payload = None
    if comm.rank == 0:
        print("\n[GEO-34 1b-probe] === per-sheet table, columns P1..P4 (Pi vs P1 "
              "rotated (i-1)*90 deg about z) ===")
        sizes = {"primitive": g_size, "import": i_size}
        for name in ("primitive", "import"):
            _print_table(rows[name], name, sizes[name], narrowed[name])

        print("\n[GEO-34 1b-probe] === import vs primitive, sheet by sheet (no rotation) ===")
        p, q = rows["primitive"], rows["import"]
        for k in range(len(PORTS)):
            dc = geo31._hausdorff(q["sheets"][k]["centroids"], p["sheets"][k]["centroids"])
            db = geo31._hausdorff(q["sheets"][k]["boundary"], p["sheets"][k]["boundary"])
            bb = geo31._bbox_discrepancy(q["sheets"][k]["vertices"], p["sheets"][k]["vertices"])
            print(f"    P{k + 1}: d_centroid(imp,prim)={dc:.6e} m ({dc / p['ell']:.4f} l)  "
                  f"d_boundary(imp,prim)={db:.6e} m  bbox delta={bb:.6e} m")

        print("\n[GEO-34 1b-probe] === discriminator (GEO-31's scale-free thresholds, "
              "1e-3 l; printed, asserted nowhere) ===")
        for name in ("primitive", "import"):
            r = rows[name]
            ell = r["ell"]
            sup = max(r["bbox_delta_all"] + r["bbox_delta_boundary"])
            bdry = max(r["d_boundary"])
            cut = max(r["d_centroid"])
            if sup > 1.0e-3 * ell:
                verdict = "DIFFERENT SUPPORT"
            elif bdry > 1.0e-3 * ell:
                verdict = "SAME SUPPORT, DIFFERENT BOUNDARY"
            elif cut > 1.0e-3 * ell:
                verdict = "SAME SUPPORT AND BOUNDARY, DIFFERENT INTERIOR CUT"
            else:
                verdict = "C4-CONGRUENT TO ROUND-OFF"
            print(f"  {name:<10s} max bbox delta={sup:.6e} m  max d_boundary={bdry:.6e} m  "
                  f"max d_centroid={cut:.6e} m ({cut / ell:.3f} l)  => {verdict}")

        ctl = rows["import"]
        print("\n[GEO-34 1b-probe] === self-checks ===")
        print(f"  (A) primitive counts={tuple(p['counts'])} vs GEO-31 {GEO31_X1_COUNTS}; "
              f"area spread={geo31._spread(p['areas']):.6e} vs 1e-3; "
              f"d_i={tuple(f'{v:.6e}' for v in p['d_centroid'])} vs {GEO31_X1_D}")
        print(f"  (B) import sheet P1: l={ctl['ell']:.6e} m; d(P1,P1)={ctl['ctl_zero']:.6e} m "
              f"(must be 0.0); d(P1,P1+l*(1,1,1)/sqrt3)={ctl['ctl_shift']:.6e} m = "
              f"{ctl['ctl_shift'] / ctl['ell']:.6f} l (bar 0.5 l)")
        print(f"[GEO-34 1b-probe] wall: primitive+export {t_prim:.2f} s, import {t_imp:.2f} s, "
              f"measurement {t_meas:.2f} s, total {time.perf_counter() - started:.2f} s",
              flush=True)

        d = p["d_centroid"]
        fails = []
        if tuple(int(c) for c in p["counts"]) != GEO31_X1_COUNTS:
            fails.append(f"(A) counts {tuple(p['counts'])}")
        if geo31._spread(p["areas"]) > geo31.SHEET_AREA_BAND:
            fails.append("(A) area spread")
        if d[0] != 0.0:
            fails.append(f"(A) d_1={d[0]!r} not exactly 0")
        if not d[1] <= D_ROUNDOFF:
            fails.append(f"(A) d_2={d[1]:.3e} not round-off")
        for k in (2, 3):
            if not abs(d[k] - GEO31_X1_D[k]) <= D_REL * GEO31_X1_D[k]:
                fails.append(f"(A) d_{k + 1}={d[k]:.9e} vs {GEO31_X1_D[k]:.6e}")
        if ctl["ctl_zero"] != 0.0:
            fails.append("(B) self-distance not 0")
        if not ctl["ctl_shift"] >= 0.5 * ctl["ell"]:
            fails.append("(B) shifted distance under l/2")
        payload = fails

    fails = comm.bcast(payload, root=0)
    if fails:
        _print("[GEO-34 1b-probe] SELF-CHECK FAILED — table uncountable: " + "; ".join(fails), comm)
        sys.exit(1)
    _print("[GEO-34 1b-probe] self-checks (A) and (B) green; every symmetry reading "
           "above is printed and asserted nowhere.", comm)


if __name__ == "__main__":
    main()
