"""`IMP-3` step 3 (§9 item 235, chain I9) — the Palace twin of the rod fixture.

Two subcommands, both serial in the complex ``fem-em-solver`` container:

* ``export`` (W1) — builds the exact-rectangle sibling of `IMP-3`'s fixture
  (``_build(True, leg_sheet_strip_fraction=GATED_WIDTH_FRACTION)`` then
  ``_build_rod(leg_sheet_strip_fraction=...)``, in that order), takes the
  driven sheets N_i from ``build_four_port_sweep(..., build_only=True)`` and
  writes Gmsh 2.2 meshes ``imp3x_off.msh`` / ``imp3x_rod.msh`` to the
  gitignored ``logs/imp3/``: box exterior -> 1, N_i -> 20 + i, rod wall -> 402.
  Configs ``imp3x_off`` (PEC [1]), ``imp3x_rod`` (PEC [1, 402]) and
  ``imp3x_rodnat`` (PEC [1], 402 left unassigned = Palace's natural wall, the
  form item 234's ``h05x_pmc`` ran).  Anchors (X1)-(X5) asserted.
* ``read`` (W5) — (P) asserted on ``imp3x_off`` / ``imp3x_rod`` via the OPS-61
  probe's ``_p_check``; ΔS(rod), ΔS_wall, σ_y pairs and the f05x tie printed.

Imports helpers from ``ops61_step1a_palace_ans4.py`` by importlib; never edits
it and never writes to ``logs/ops61/``.
"""

from __future__ import annotations

import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

OUT_DIR = REPO / "logs" / "imp3"
ROD_WALL_TAG = 402
ARMS = ("off", "rod")
RUNS = ("imp3x_off", "imp3x_rod", "imp3x_rodnat")
NOMINAL_STRIP_AREA = 5.6e-5


def _ops61():
    path = REPO / "scripts" / "probes" / "ops61_step1a_palace_ans4.py"
    spec = importlib.util.spec_from_file_location("ops61_step1a_palace_ans4", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def build_arm(arm):
    """The build tuple for ``arm`` ('off' | 'rod'), exact-rectangle ports."""
    from tests.mesh.test_birdcage_implant_rod import _build_rod
    from tests.mesh.test_birdcage_port_sheets import _build
    from tests.validation.test_port_lumped_narrowed_sheet import GATED_WIDTH_FRACTION

    if arm == "off":
        return _build(True, leg_sheet_strip_fraction=GATED_WIDTH_FRACTION)
    return _build_rod(leg_sheet_strip_fraction=GATED_WIDTH_FRACTION)


def export():
    import dolfinx
    from mpi4py import MPI

    from tests.mesh.helpers import global_cell_tag_set
    from tests.mesh.test_birdcage_implant_rod import _rod_area
    from tests.mesh.test_birdcage_leg_gaps import TERMINAL_AREA_BAND
    from tests.mesh.test_birdcage_port_sheet_prerequisite import CELL_COUNT_BAND
    from tests.mesh.test_birdcage_port_sheets import SHEET_IFACE
    from tests.mesh.test_birdcage_port_tags import LEG_COUNT
    from tests.validation.test_imp3_rod_birdcage import ROD_CELL_GUARD_W2
    from tests.validation.test_lossy_sphere_fullwave import SALINE_EPSILON_R, SALINE_SIGMA
    from tests.validation.test_ops61_exact_strip_ports import (
        ANS4_EXACT_STRIP_CELL_RECORD,
        STRIP_BAND,
    )
    from tests.validation.test_birdcage_power_identity import STEP2D_FREQUENCY_HZ
    from tests.validation.test_port_birdcage_four_port import build_four_port_sweep
    from tests.validation.test_port_birdcage_lumped_column import (
        CONDUCTOR_CELL_TAG,
        PHANTOM_CELL_TAG,
    )
    from tests.validation.test_port_gap_voltage_impedance import SIGMA_WIRE_S_PER_M
    from tests.validation.test_port_package_sparameters import REFERENCE_IMPEDANCE_OHM

    o61 = _ops61()
    comm = MPI.COMM_WORLD
    assert comm.size == 1, "export is a serial writer: run it at mpiexec -n 1"
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    targets = [OUT_DIR / f"imp3x_{a}.msh" for a in ARMS] + [
        OUT_DIR / f"{r}.json" for r in RUNS
    ] + [OUT_DIR / "imp3x_meta.json"]
    for t in targets:
        assert not t.exists(), f"{t} exists: never overwrite"
    ports_idx = list(range(1, LEG_COUNT + 1))
    meta = {}
    fails = []
    tags_by_arm = {}
    t_all = time.perf_counter()
    for arm in ARMS:
        built = build_arm(arm)
        msh, cell_tags, build_facets, diag, t_build = built
        case = build_four_port_sweep(
            frequency_hz=STEP2D_FREQUENCY_HZ, built=built, build_only=True
        )
        tdim = msh.topology.dim
        fdim = tdim - 1
        ncells = int(msh.topology.index_map(tdim).size_global)
        print(f"\n[IMP-3 3 W1] arm {arm}: {ncells} cells, build {t_build:.2f} s; "
              f"sweep case cells {int(case['cells'])}", flush=True)
        assert int(case["cells"]) == ncells
        if arm == "off":
            rel = abs(ncells - ANS4_EXACT_STRIP_CELL_RECORD) / ANS4_EXACT_STRIP_CELL_RECORD
            ok = rel <= CELL_COUNT_BAND
            print(f"[IMP-3 3 W1] (X1) {ncells} vs ANS4_EXACT_STRIP_CELL_RECORD "
                  f"{ANS4_EXACT_STRIP_CELL_RECORD}: rel {rel:.3e} (band {CELL_COUNT_BAND}); "
                  f"exact equality {ncells == ANS4_EXACT_STRIP_CELL_RECORD}: "
                  f"{'PASS' if ok else 'FAIL'}", flush=True)
            if not ok:
                fails.append(f"X1 {ncells}")
        else:
            ok = ncells < ROD_CELL_GUARD_W2
            print(f"[IMP-3 3 W1] (X2) {ncells} < ROD_CELL_GUARD_W2 {ROD_CELL_GUARD_W2}: "
                  f"{'PASS' if ok else 'FAIL'}", flush=True)
            assert ok, "X2 red: stop before writing"
        msh.topology.create_connectivity(fdim, tdim)
        x = np.asarray(msh.geometry.x, dtype=float)

        def facet_nodes(facets):
            facets = np.asarray(facets, dtype=np.int32)
            if not facets.size:
                return np.zeros((0, 3), dtype=np.int64)
            return np.asarray(
                dolfinx.cpp.mesh.entities_to_geometry(msh._cpp_object, fdim, facets, False)
            ).reshape(-1, 3)

        ext = np.asarray(dolfinx.mesh.exterior_facet_indices(msh.topology), dtype=np.int32)
        wall = np.asarray(build_facets.find(ROD_WALL_TAG), dtype=np.int32)
        outer = np.setdiff1d(ext, wall)
        wall_ext = np.intersect1d(wall, ext)
        outer_nodes, wall_nodes = facet_nodes(outer), facet_nodes(wall)
        outer_area = float(np.sum(o61._tri_areas(x, outer_nodes)))
        wall_area = float(np.sum(o61._tri_areas(x, wall_nodes)))
        tri_blocks = [(outer_nodes, o61.EXTERIOR_TAG)]
        port_meta = {}
        sheet_facets = case["facet_tags"]
        for i in ports_idx:
            n_i = np.asarray(sheet_facets.find(SHEET_IFACE + i), dtype=np.int32)
            nodes = facet_nodes(n_i)
            a = float(np.sum(o61._tri_areas(x, nodes)))
            zs = x[nodes.ravel(), 2] if nodes.size else np.zeros(1)
            h = float(zs.max() - zs.min())
            ok = abs(a - NOMINAL_STRIP_AREA) <= STRIP_BAND
            print(f"    P{i}: N_{i} {n_i.size} facets area {a:.12e} m^2 h {h:.9e} "
                  f"w = A/h {a / h:.9e}; (X3) |A - 5.6e-5| {abs(a - NOMINAL_STRIP_AREA):.3e} "
                  f"<= {STRIP_BAND}: {'PASS' if ok else 'FAIL'}", flush=True)
            if not ok:
                fails.append(f"X3 {arm} P{i} {a}")
            port_meta[f"P{i}"] = {"facets": int(n_i.size), "area": a, "h": h}
            tri_blocks.append((nodes, o61.NARROW_BASE + i))
        ext_err = abs(outer_area - o61.EXTERIOR_AREA_M2)
        ok = ext_err <= o61.EXTERIOR_AREA_BAND
        print(f"[IMP-3 3 W1] (X4) outer exterior {outer.size} tris area {outer_area:.12e} "
              f"(|diff| {ext_err:.3e}, band {o61.EXTERIOR_AREA_BAND}): "
              f"{'PASS' if ok else 'FAIL'}", flush=True)
        if not ok:
            fails.append(f"X4 {arm} {outer_area}")
        if arm == "rod":
            ratio = wall_area / _rod_area()
            ok = (wall.size > 0 and wall_ext.size == wall.size
                  and TERMINAL_AREA_BAND[0] <= ratio <= TERMINAL_AREA_BAND[1])
            print(f"[IMP-3 3 W1] (X5) rod wall {wall.size} facets ({wall_ext.size} exterior), "
                  f"area {wall_area:.9e} / analytic {_rod_area():.9e} = {ratio:.6f} in "
                  f"{TERMINAL_AREA_BAND}: {'PASS' if ok else 'FAIL'}", flush=True)
            tri_blocks.append((wall_nodes, ROD_WALL_TAG))
        else:
            ok = wall.size == 0
            print(f"[IMP-3 3 W1] (X5) without-rod 402 count {wall.size} == 0: "
                  f"{'PASS' if ok else 'FAIL'}", flush=True)
        if not ok:
            fails.append(f"X5 {arm}")

        # --- Gmsh 2.2 ASCII (copied from the OPS-61 export, :257-292) --------
        cells = np.asarray(msh.geometry.dofmap, dtype=np.int64).reshape(-1, 4)
        ncl = msh.topology.index_map(tdim).size_local
        cells = cells[:ncl]
        ctag = np.zeros(ncl, dtype=np.int64)
        ctag[np.asarray(cell_tags.indices)] = np.asarray(cell_tags.values)
        assert np.all(ctag > 0), "an untagged tet"
        p = x[cells]
        vol6 = np.einsum(
            "ij,ij->i", np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]), p[:, 3] - p[:, 0]
        )
        neg = vol6 < 0
        cells[neg, 2], cells[neg, 3] = cells[neg, 3].copy(), cells[neg, 2].copy()
        mesh_path = OUT_DIR / f"imp3x_{arm}.msh"
        ntri = sum(b.shape[0] for b, _ in tri_blocks)
        with open(mesh_path, "x") as fh:
            fh.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n")
            fh.write(f"$Nodes\n{x.shape[0]}\n")
            for k, (a, b, c) in enumerate(x, start=1):
                fh.write(f"{k} {a:.17g} {b:.17g} {c:.17g}\n")
            fh.write("$EndNodes\n")
            fh.write(f"$Elements\n{ntri + ncl}\n")
            eid = 1
            for nodes, attr in tri_blocks:
                for t in nodes:
                    fh.write(f"{eid} 2 2 {attr} {attr} {t[0]+1} {t[1]+1} {t[2]+1}\n")
                    eid += 1
            for t, a in zip(cells, ctag):
                fh.write(f"{eid} 4 2 {a} {a} {t[0]+1} {t[1]+1} {t[2]+1} {t[3]+1}\n")
                eid += 1
            fh.write("$EndElements\n")
        tags = sorted(global_cell_tag_set(msh, cell_tags))
        tags_by_arm[arm] = tags
        print(f"[IMP-3 3 W1] wrote {mesh_path.relative_to(REPO)}: {x.shape[0]} nodes, "
              f"{ncl} tets ({int(neg.sum())} reoriented), {ntri} tris; cell tags "
              f"{ {t: int(np.count_nonzero(ctag == t)) for t in tags} }", flush=True)
        meta[arm] = {"cells": ncells, "ports": port_meta, "wall_facets": int(wall.size),
                     "wall_area": wall_area, "outer_area": outer_area}
        del case, built, msh, cell_tags, build_facets

    assert SIGMA_WIRE_S_PER_M == 800.0 and SALINE_SIGMA == 0.5
    assert SALINE_EPSILON_R == 78.0 and REFERENCE_IMPEDANCE_OHM == 50.0

    def config(name, arm, pec):
        tags = tags_by_arm[arm]
        air = [t for t in tags if t not in (CONDUCTOR_CELL_TAG, PHANTOM_CELL_TAG)]
        materials = [
            {"Attributes": [CONDUCTOR_CELL_TAG], "Permeability": 1.0,
             "Permittivity": 1.0, "Conductivity": SIGMA_WIRE_S_PER_M},
            {"Attributes": [PHANTOM_CELL_TAG], "Permeability": 1.0,
             "Permittivity": SALINE_EPSILON_R, "Conductivity": SALINE_SIGMA},
            {"Attributes": air, "Permeability": 1.0, "Permittivity": 1.0},
        ]
        return {
            "Problem": {"Type": "Driven", "Verbose": 2, "Output": f"postpro_{name}"},
            "Model": {"Mesh": f"imp3x_{arm}.msh", "L0": 1.0},
            "Domains": {"Materials": materials},
            "Boundaries": {
                "PEC": {"Attributes": pec},
                "LumpedPort": [
                    {"Index": i, "R": REFERENCE_IMPEDANCE_OHM, "Excitation": i,
                     "Attributes": [o61.NARROW_BASE + i], "Direction": "+Z"}
                    for i in ports_idx
                ],
            },
            "Solver": {
                "Order": 1,
                "Device": "CPU",
                "Driven": {"Samples": [{"Type": "Point", "Freq": [o61.FREQUENCY_GHZ]}]},
                "Linear": {"Type": "MUMPS", "KSPType": "GMRES", "Tol": 1.0e-10,
                           "MaxIts": 100},
            },
        }

    configs = {
        "imp3x_off": config("imp3x_off", "off", [o61.EXTERIOR_TAG]),
        "imp3x_rod": config("imp3x_rod", "rod", [o61.EXTERIOR_TAG, ROD_WALL_TAG]),
        # 402 left unassigned: Palace's natural (PMC) default, item 234's form.
        "imp3x_rodnat": config("imp3x_rodnat", "rod", [o61.EXTERIOR_TAG]),
    }
    for name, cfg in configs.items():
        with open(OUT_DIR / f"{name}.json", "x") as fh:
            fh.write(json.dumps(cfg, indent=2) + "\n")
        print(f"[IMP-3 3 W1] wrote logs/imp3/{name}.json: PEC "
              f"{cfg['Boundaries']['PEC']['Attributes']}", flush=True)
    print(json.dumps(configs["imp3x_rod"], indent=2), flush=True)
    with open(OUT_DIR / "imp3x_meta.json", "x") as fh:
        fh.write(json.dumps(meta, indent=2) + "\n")
    print(f"[IMP-3 3 W1] meta {json.dumps(meta)}", flush=True)
    print(f"[IMP-3 3 W1] export total {time.perf_counter() - t_all:.2f} s", flush=True)
    assert not fails, f"(X1)-(X5) red: {fails}"
    print("[IMP-3 3 W1] ANCHORS PASS: (X1)-(X5)", flush=True)


def csv_of(name):
    return OUT_DIR / f"postpro_{name}" / "port-S.csv"


def class_means(s):
    from tests.validation.test_port_birdcage_four_port import _circulant_classes

    return {k: complex(np.mean(v)) for k, v in _circulant_classes(s).items()}


def rel_class_gap(a, b, ref=None):
    """|mean_k(a) - mean_k(b)| / |mean_k(ref or b)| per class."""
    ma, mb = class_means(a), class_means(b)
    mr = class_means(ref) if ref is not None else mb
    return {k: abs(ma[k] - mb[k]) / abs(mr[k]) for k in ma}


def delta_vs(d_ours, d_pal):
    """Per-class |Δ_ours - Δ_pal| / |Δ_pal| on class-mean deltas."""
    return {k: abs(d_ours[k] - d_pal[k]) / abs(d_pal[k]) for k in d_pal}


def class_delta(a, b):
    ma, mb = class_means(a), class_means(b)
    return {k: ma[k] - mb[k] for k in ma}


def line(d):
    return "  ".join(f"{k} {v:.3e}" for k, v in d.items())


def read_palace():
    """The three Palace matrices (raises if a CSV is missing)."""
    o61 = _ops61()
    return {name: o61._read_port_s(csv_of(name))[0] for name in RUNS}


def read():
    from tests.validation.test_imp3_rod_birdcage import MIRROR_PAIRS
    from tests.validation.test_port_birdcage_four_port import _class_spread
    from tests.validation.test_port_birdcage_lumped_column import ADJACENT_SPREAD_BAND

    o61 = _ops61()
    tag = "[IMP-3 3 W5]"
    s, fails = {}, []
    for name in RUNS:
        s[name], _ident, f = o61._p_check(name, csv_of(name), tag)
        if name != "imp3x_rodnat":
            fails += [f"{name}: {x}" for x in f]
        else:
            print(f"{tag} imp3x_rodnat (P) printed only (control)", flush=True)
    off, rod, nat = s["imp3x_off"], s["imp3x_rod"], s["imp3x_rodnat"]
    d_rod = class_delta(rod, off)
    d_wall = class_delta(rod, nat)
    mo = class_means(off)
    print(f"{tag} Palace ΔS(rod)/|mean_k(off)|  {line({k: abs(v) / abs(mo[k]) for k, v in d_rod.items()})}"
          f"  (predicted [1e-5, 1e-3])", flush=True)
    print(f"{tag} Palace ΔS_wall/|mean_k(off)| {line({k: abs(v) / abs(mo[k]) for k, v in d_wall.items()})}"
          f"  (predicted [1e-6, 1e-3])", flush=True)
    for name in RUNS:
        mp = [_class_spread([s[name][a], s[name][b]]) for a, b in MIRROR_PAIRS]
        print(f"{tag} {name} σ_y pairs on S: " + ", ".join(f"{v:.3e}" for v in mp)
              + f" (ADJACENT_SPREAD_BAND {ADJACENT_SPREAD_BAND})", flush=True)
    f05x = REPO / "logs" / "ops61" / "postpro_f05x" / "port-S.csv"
    if f05x.exists():
        sf = o61._read_port_s(f05x)[0]
        print(f"{tag} Palace imp3x_off vs f05x per class {line(rel_class_gap(off, sf))}"
              f" (free tie, predicted <= 1e-9)", flush=True)
    else:
        print(f"{tag} f05x tie: not measured (no {f05x.relative_to(REPO)})", flush=True)
    assert not fails, f"(P) red: {fails}"
    print(f"{tag} (P) PASS on imp3x_off and imp3x_rod", flush=True)


if __name__ == "__main__":
    if len(sys.argv) != 2 or sys.argv[1] not in ("export", "read"):
        raise SystemExit("usage: imp3_step3_palace_rod.py export|read")
    {"export": export, "read": read}[sys.argv[1]]()
