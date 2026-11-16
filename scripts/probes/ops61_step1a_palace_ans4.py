"""`OPS-61` step 1a (§9 item 75, F10) — Palace on the `ANS-4` x1 rung at 64 MHz.

Two subcommands, both run at ``mpiexec -n 1`` in the complex ``fem-em-solver``
container (the Palace windows W2-W4 run in ``fem-em-palace`` between them):

* ``export`` (W1) — builds the gated four-leg birdcage exactly as
  ``_four_port_rung`` does (``_build(np.zeros(4))`` then its sheet steps,
  ``test_port_birdcage_leg_offset_sweep.py:271-307``), so it carries both the
  full sheet F_i and the narrowed N_i (``GATED_WIDTH_FRACTION`` = 0.5) we drive.
  The generator finalizes gmsh and the sheets are dolfinx tags, so the mesh is
  written as **Gmsh 2.2 ASCII from dolfinx** (``entities_to_geometry`` /
  the geometry dofmap, never topology vertex ids) to the gitignored
  ``logs/ops61/``.  Tets carry their cell tag; triangles carry 1 (exterior),
  20+i (N_i) or 30+i (F_i minus N_i).  Three Palace configs are written:
  ``f1`` (port i = {20+i, 30+i}, the exact rectangle), ``f05`` (port i = {20+i},
  the ragged set we drive) and ``flip`` (``f1`` with P1 along -Z).
  Asserted: exactly ``STEP2_CELL_COUNT`` tets; each F_i area equal to its
  nominal rectangle within ``SHEET_AREA_BAND``; exterior area 0.3072 m^2 to
  1e-9.

* ``export --exact-strip`` (`OPS-61` step 1b, W3) — the same export on
  ``_build(..., leg_sheet_strip_fraction=GATED_WIDTH_FRACTION)``, so N_i is the
  exact centre rectangle. Writes ``ops61_ans4_x1_strip.msh`` and one config,
  ``f05x`` (``Model.Mesh`` that file; port i = {20+i}, ``+Z``); never touches
  step 1a's mesh or ``f1`` / ``f05`` / ``flip``. Asserts the tet count equals
  ``ANS4_EXACT_STRIP_CELL_RECORD`` instead of 116 085, plus the F_i and
  exterior-area anchors above.

* ``read`` (W5) — turns each run's ``port-S.csv`` (dB / degrees, parsed by
  header) into the 4x4 and reads it with ``_circulant_classes`` /
  ``_class_spread``.  Asserted on ``f1``: reciprocity <= ``RECIPROCITY_BAND``,
  sigma_max <= 1 + ``PASSIVITY_SIGMA_TOLERANCE``, every S-class spread <=
  ``ADJACENT_SPREAD_BAND``.  The ``flip`` negative control (W4 = D S D,
  D = diag(-1, 1, 1, 1)) is *predicted*, printed beside the measurement, never
  asserted (§9 rule (e)).  Class means against ``ans4`` ``metrics.json``
  ``rungs["64 MHz"]["s_matrix"]`` are printed, never gated, and confounded by
  port width for ``f1``.  Never reads ``z_matrix_ohm`` (terminated
  transimpedance, 2026-09-29 ruling) nor the ``aed`` key.

* ``export --hole-pec`` / ``read --hole-pec`` (`OPS-61` step 2a, §9 item 234,
  F22) — ``export_hole_pec()`` writes `ANS-6`'s exact-rectangle hole mesh
  (``_sheets_build(True, as_hole=True, leg_sheet_strip_fraction=...)``) as
  ``ops61_ans6_x1_hole_strip.msh`` with box -> 1, cavity wall -> 401, N_i ->
  20 + i, plus ``h05x_pec`` / ``h05x_pmc`` / ``h05x_pec_o2`` configs from
  ``f05x.json``; never overwrites, never calls ``export()``.  Anchors (E1)-(E5).
  ``read_hole_pec()`` asserts (P) on ``h05x_pec`` (+o2 if present) and the
  (K') PMC-wall floor; ours-vs-Palace (T1)/(T2)/(C3)/(C2''') printed only.

* ``cu-configs`` / ``read --hole-cu`` (`OPS-61` step 2b, §9 item 241, F23) —
  ``cu_configs()`` writes ``h05x_cu`` / ``h05x_cu_o2`` / ``h05x_cu11`` from
  ``h05x_pec.json`` (PEC -> [1] only, ``Boundaries.Conductivity`` on the wall
  with ``th14.COPPER_SIGMA`` / ``SIGMA_LADDER[-1]``, mu_r 1, no Thickness);
  never overwrites; asserts (Z), the thick-conductor (1 + j)/(sigma delta) equal
  to our ``(1 + j) surface_resistance_ohm`` within ``STRIP_BAND``.
  ``read_hole_cu()`` asserts (P) on ``h05x_cu`` (+o2 if present), the (K1)
  floor against ``h05x_pec`` and (K3) sigma_max(cu) < sigma_max(pec); ours-vs-
  Palace (T3)/(T3m)/(T4)/(C2'''') and the (K2) sigma-scaling ratio printed only.

Run (W1)::

    scripts/testing/run_and_log.sh OPS-61-step1a-W1 "docker compose exec -T fem-em-solver \\
      bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && \\
      PYTHONPATH=/workspace/src timeout -k 30 300 mpiexec -n 1 \\
      python3 scripts/probes/ops61_step1a_palace_ans4.py export'"
"""

from __future__ import annotations

import csv
import json
import re
import sys
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

OUT_DIR = REPO / "logs" / "ops61"
MESH_NAME = "ops61_ans4_x1.msh"
STRIP_MESH_NAME = "ops61_ans4_x1_strip.msh"  # step 1b W3 (``export --exact-strip``)
FREQUENCY_GHZ = 0.064
EXTERIOR_TAG = 1
NARROW_BASE = 20  # N_i  -> 20 + i
REMAINDER_BASE = 30  # F_i \ N_i -> 30 + i
EXTERIOR_AREA_M2 = 0.3072
EXTERIOR_AREA_BAND = 1.0e-9
FLIP_CONTROL_BAND = 1.0e-6  # predicted, never asserted
ANS4_METRICS = (
    REPO
    / "examples"
    / "ansys_benchmarks"
    / "ans4_birdcage_four_port_10_64_128MHz"
    / "metrics.json"
)


def _tri_areas(x, tris):
    a, b, c = x[tris[:, 0]], x[tris[:, 1]], x[tris[:, 2]]
    return 0.5 * np.linalg.norm(np.cross(b - a, c - a), axis=1)


def export(exact_strip=False):
    """W1 export; ``exact_strip`` is `OPS-61` step 1b's W3 mode (see docstring)."""
    import dolfinx
    from mpi4py import MPI

    from tests.mesh.helpers import global_cell_tag_set
    from tests.mesh.test_birdcage_leg_offset import (
        SHEET_AREA_BAND,
        _build,
        _projected_extents,
        _sheet_areas,
        _sheet_azimuth_deg,
    )
    from tests.mesh.test_birdcage_port_sheets import SHEET_IFACE
    from tests.mesh.test_birdcage_port_tags import LEG_COUNT
    from tests.validation.test_lossy_sphere_fullwave import (
        SALINE_EPSILON_R,
        SALINE_SIGMA,
    )
    from tests.validation.test_port_birdcage_leg_offset_sweep import (
        _narrowed_radial,
        _port_frame,
    )
    from tests.validation.test_port_birdcage_lumped_column import (
        CONDUCTOR_CELL_TAG,
        PHANTOM_CELL_TAG,
        STEP2_CELL_COUNT,
    )
    from tests.validation.test_port_gap_voltage_impedance import SIGMA_WIRE_S_PER_M
    from tests.validation.test_port_lumped_narrowed_sheet import GATED_WIDTH_FRACTION
    from tests.validation.test_port_package_sparameters import REFERENCE_IMPEDANCE_OHM

    comm = MPI.COMM_WORLD
    assert comm.size == 1, "export is a serial writer: run it at mpiexec -n 1"
    t_start = time.perf_counter()
    ports_idx = list(range(1, LEG_COUNT + 1))

    if exact_strip:
        from tests.validation.test_ops61_exact_strip_ports import (
            ANS4_EXACT_STRIP_CELL_RECORD,
        )

        mesh_name = STRIP_MESH_NAME
        expected_cells = ANS4_EXACT_STRIP_CELL_RECORD
        msh, cell_tags, _facet_tags, diag, t_mesh = _build(
            np.zeros(LEG_COUNT), leg_sheet_strip_fraction=GATED_WIDTH_FRACTION
        )
    else:
        mesh_name = MESH_NAME
        expected_cells = STEP2_CELL_COUNT
        msh, cell_tags, _facet_tags, diag, t_mesh = _build(np.zeros(LEG_COUNT))
    tdim = msh.topology.dim
    fdim = tdim - 1
    ncells = int(msh.topology.index_map(tdim).size_global)
    print(
        f"[OPS-61 1a W1] mesh built in {t_mesh:.2f} s: {ncells} tets "
        f"(expected {expected_cells}; exact_strip={exact_strip}); diag "
        f"port_box_size_m = {diag['port_box_size_m']}",
        flush=True,
    )
    # Anchor 1 — negative-result rule: a count other than the record stops
    # here, before any config is written and before any Palace window.
    assert ncells == expected_cells, (
        f"W1 count {ncells} != {expected_cells}: stop before any Palace window"
    )
    msh.topology.create_connectivity(fdim, tdim)
    msh.topology.create_entity_permutations()

    # --- `_four_port_rung`'s sheet steps, verbatim (lines 271-307) ----------
    tags_f, full_areas = _sheet_areas(msh, cell_tags, ports_idx, comm)
    dx, _dy, dz = diag["port_box_size_m"]
    sheet_analytic = float(dx) * float(dz)

    geometry = {}
    for i in ports_idx:
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
            "out_of_plane_full": float(spans[1]),
            "h_full": float(spans[2]),
            "centre_radial": float(centres[0]),
            "area_full": float(full_areas[i]),
            "area_ratio_full": float(full_areas[i] / sheet_analytic),
        }

    tags_full = tags_f
    for i in ports_idx:
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
    # ------------------------------------------------------------------------

    x = np.asarray(msh.geometry.x, dtype=float)

    def facet_nodes(facets):
        facets = np.asarray(facets, dtype=np.int32)
        if not facets.size:
            return np.zeros((0, 3), dtype=np.int64)
        return np.asarray(
            dolfinx.cpp.mesh.entities_to_geometry(msh._cpp_object, fdim, facets, False)
        ).reshape(-1, 3)

    tri_blocks = []  # (node array, attribute)
    ext = dolfinx.mesh.exterior_facet_indices(msh.topology)
    ext_nodes = facet_nodes(ext)
    ext_area = float(np.sum(_tri_areas(x, ext_nodes)))
    tri_blocks.append((ext_nodes, EXTERIOR_TAG))

    print(
        f"[OPS-61 1a W1] nominal sheet rectangle dx*dz = {dx:.6g} x {dz:.6g} = "
        f"{sheet_analytic:.9e} m^2; GATED_WIDTH_FRACTION = {GATED_WIDTH_FRACTION}",
        flush=True,
    )
    area_ok = True
    for i in ports_idx:
        tag = SHEET_IFACE + i
        full = np.asarray(tags_full.find(tag), dtype=np.int32)
        narrow = np.asarray(tags_f.find(tag), dtype=np.int32)
        assert np.all(np.isin(narrow, full)), f"N_{i} not a subset of F_{i}"
        rest = np.setdiff1d(full, narrow)
        n_nodes, r_nodes = facet_nodes(narrow), facet_nodes(rest)
        a_n = float(np.sum(_tri_areas(x, n_nodes)))
        a_r = float(np.sum(_tri_areas(x, r_nodes)))
        a_f_numpy = a_n + a_r
        ratio = geometry[i]["area_ratio_full"]
        ok = abs(ratio - 1.0) < SHEET_AREA_BAND
        area_ok &= ok
        print(
            f"    P{i}: azimuth {geometry[i]['azimuth_deg']:8.4f} deg; F_{i} "
            f"{full.size} tris area {geometry[i]['area_full']:.12e} m^2 "
            f"(numpy {a_f_numpy:.12e}), meshed/nominal {ratio:.12f} "
            f"[{'PASS' if ok else 'FAIL'} band {SHEET_AREA_BAND:.0e}]; "
            f"w_full {geometry[i]['w_full']:.9e} h_full {geometry[i]['h_full']:.9e}; "
            f"N_{i} {narrow.size} tris area {a_n:.9e} m^2 "
            f"(fraction {a_n / a_f_numpy:.6f}); F_{i}\\N_{i} {rest.size} tris",
            flush=True,
        )
        tri_blocks.append((n_nodes, NARROW_BASE + i))
        tri_blocks.append((r_nodes, REMAINDER_BASE + i))
    ext_err = abs(ext_area - EXTERIOR_AREA_M2)
    print(
        f"[OPS-61 1a W1] exterior: {ext.size} tris, area {ext_area:.12e} m^2 vs "
        f"{EXTERIOR_AREA_M2} (|diff| {ext_err:.3e}, band {EXTERIOR_AREA_BAND:.0e})",
        flush=True,
    )

    # --- Gmsh 2.2 ASCII -----------------------------------------------------
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    cells = np.asarray(msh.geometry.dofmap, dtype=np.int64).reshape(-1, 4)
    ncl = msh.topology.index_map(tdim).size_local
    cells = cells[:ncl]
    ctag = np.zeros(ncl, dtype=np.int64)
    ctag[np.asarray(cell_tags.indices)] = np.asarray(cell_tags.values)
    assert np.all(ctag > 0), "an untagged tet"
    # Positive orientation (MFEM is happiest with it): swap two nodes where the
    # signed volume is negative.
    p = x[cells]
    vol6 = np.einsum(
        "ij,ij->i", np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]), p[:, 3] - p[:, 0]
    )
    neg = vol6 < 0
    cells[neg, 2], cells[neg, 3] = cells[neg, 3].copy(), cells[neg, 2].copy()
    total_vol = float(np.sum(np.abs(vol6)) / 6.0)

    mesh_path = OUT_DIR / mesh_name
    ntri = sum(b.shape[0] for b, _ in tri_blocks)
    with open(mesh_path, "w") as fh:
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

    lo, hi = x.min(axis=0), x.max(axis=0)
    tags = sorted(global_cell_tag_set(msh, cell_tags))
    counts = {t: int(np.count_nonzero(ctag == t)) for t in tags}
    print(
        f"[OPS-61 1a W1] wrote {mesh_path.relative_to(REPO)}: {x.shape[0]} nodes, "
        f"{ncl} tets ({int(neg.sum())} reoriented), {ntri} tris; volume "
        f"{total_vol:.9e} m^3; bbox lo {lo} hi {hi}",
        flush=True,
    )
    print(f"[OPS-61 1a W1] cell tags (count): {counts}", flush=True)

    # --- Palace configs (key spellings from W0: config-schema.json v0.18.1,
    #     examples/cpw/cpw_lumped_uniform.json, examples/coaxial) -----------
    assert SIGMA_WIRE_S_PER_M == 800.0 and SALINE_SIGMA == 0.5
    assert SALINE_EPSILON_R == 78.0 and REFERENCE_IMPEDANCE_OHM == 50.0
    air = [t for t in tags if t not in (CONDUCTOR_CELL_TAG, PHANTOM_CELL_TAG)]
    materials = [
        {
            "Attributes": [CONDUCTOR_CELL_TAG],
            "Permeability": 1.0,
            "Permittivity": 1.0,
            "Conductivity": SIGMA_WIRE_S_PER_M,
        },
        {
            "Attributes": [PHANTOM_CELL_TAG],
            "Permeability": 1.0,
            "Permittivity": SALINE_EPSILON_R,
            "Conductivity": SALINE_SIGMA,
        },
        {"Attributes": air, "Permeability": 1.0, "Permittivity": 1.0},
    ]

    def config(name, port_attrs, directions):
        return {
            "Problem": {"Type": "Driven", "Verbose": 2, "Output": f"postpro_{name}"},
            "Model": {"Mesh": mesh_name, "L0": 1.0},
            "Domains": {"Materials": materials},
            "Boundaries": {
                "PEC": {"Attributes": [EXTERIOR_TAG]},
                "LumpedPort": [
                    {
                        "Index": i,
                        "R": REFERENCE_IMPEDANCE_OHM,
                        "Excitation": i,
                        "Attributes": port_attrs(i),
                        "Direction": directions(i),
                    }
                    for i in ports_idx
                ],
            },
            "Solver": {
                "Order": 1,
                "Device": "CPU",
                "Driven": {"Samples": [{"Type": "Point", "Freq": [FREQUENCY_GHZ]}]},
                "Linear": {"Type": "MUMPS", "KSPType": "GMRES", "Tol": 1.0e-10,
                           "MaxIts": 100},
            },
        }

    if exact_strip:
        # Step 1b W3: only the exact-strip drive set; f1/f05/flip untouched.
        configs = {"f05x": config("f05x", lambda i: [NARROW_BASE + i], lambda i: "+Z")}
    else:
        configs = None
    configs = configs or {
        "f1": config(
            "f1", lambda i: [NARROW_BASE + i, REMAINDER_BASE + i], lambda i: "+Z"
        ),
        "f05": config("f05", lambda i: [NARROW_BASE + i], lambda i: "+Z"),
        "flip": config(
            "flip",
            lambda i: [NARROW_BASE + i, REMAINDER_BASE + i],
            lambda i: "-Z" if i == 1 else "+Z",
        ),
    }
    for name, cfg in configs.items():
        path = OUT_DIR / f"{name}.json"
        path.write_text(json.dumps(cfg, indent=2) + "\n")
        print(f"[OPS-61 1a W1] wrote {path.relative_to(REPO)}", flush=True)
    print(json.dumps(configs["f05x" if exact_strip else "f1"], indent=2), flush=True)

    print(
        f"[OPS-61 1a W1] export total {time.perf_counter() - t_start:.2f} s "
        f"(mesh {t_mesh:.2f} s)",
        flush=True,
    )
    # Anchors 2 and 3 (after writing, so a failing area still leaves the
    # artefacts to inspect — but the harness exit is non-zero either way).
    assert area_ok, "an F_i area is outside SHEET_AREA_BAND of its nominal rectangle"
    assert ext_err <= EXTERIOR_AREA_BAND, (
        f"exterior area {ext_area:.12e} != {EXTERIOR_AREA_M2} to {EXTERIOR_AREA_BAND}"
    )
    print("[OPS-61 1a W1] ANCHORS PASS: tet count, F_i areas, exterior area", flush=True)


# ---------------------------------------------------------------------------
def _read_port_s(csv_path):
    """4x4 complex S from Palace's port-S.csv (|S| in dB, arg in degrees)."""
    with open(csv_path) as fh:
        rows = list(csv.reader(fh))
    header = [h.strip() for h in rows[0]]
    data = [float(v) for v in rows[1]]
    if len(rows) > 2:
        print(f"    note: {len(rows) - 1} data rows; reading the first", flush=True)
    mag, arg = {}, {}
    for h, v in zip(header, data):
        m = re.search(r"S\[(\d+)\]\[(\d+)\]", h)
        if not m:
            print(f"    column {h!r} = {v}", flush=True)
            continue
        key = (int(m.group(1)), int(m.group(2)))
        if "dB" in h:
            mag[key] = 10.0 ** (v / 20.0)
        elif "arg" in h:
            arg[key] = np.radians(v)
        else:
            raise ValueError(f"unparsed column {h!r}")
    n = max(max(k) for k in mag)
    s = np.full((n, n), np.nan + 0j)
    for key, a in mag.items():
        s[key[0] - 1, key[1] - 1] = a * np.exp(1j * arg[key])
    return s, header


def _identities(s):
    from tests.validation.test_port_birdcage_four_port import (
        _circulant_classes,
        _class_spread,
    )

    classes = _circulant_classes(s)
    return {
        "reciprocity": float(np.linalg.norm(s - s.T) / np.linalg.norm(s)),
        "sigma_max": float(np.linalg.svd(s, compute_uv=False).max()),
        "spreads": {k: _class_spread(v) for k, v in classes.items()},
        "means": {k: complex(np.mean(v)) for k, v in classes.items()},
    }


def _fmt(z):
    return f"{z.real:+.9e}{z.imag:+.9e}j"


def read():
    from tests.validation.test_port_birdcage_four_port import (
        PASSIVITY_SIGMA_TOLERANCE,
        _circulant_classes,
    )
    from tests.validation.test_port_birdcage_lumped_column import ADJACENT_SPREAD_BAND
    from tests.validation.test_port_lumped_sheet_sweep import RECIPROCITY_BAND

    z0 = 50.0
    rung = json.loads(ANS4_METRICS.read_text())["rungs"]["64 MHz"]
    s_ours = np.array(
        [[complex(e["re"], e["im"]) for e in row] for row in rung["s_matrix"]]
    )
    ours_means = {k: complex(np.mean(v)) for k, v in _circulant_classes(s_ours).items()}

    results = {}
    for name in ("f1", "f05", "flip"):
        path = OUT_DIR / f"postpro_{name}" / "port-S.csv"
        if not path.exists():
            print(f"[OPS-61 1a W5] {name}: no {path.relative_to(REPO)}", flush=True)
            continue
        s, header = _read_port_s(path)
        ident = _identities(s)
        results[name] = (s, ident)
        print(f"\n[OPS-61 1a W5] {name}: header {header}", flush=True)
        for r in range(s.shape[0]):
            print("    S row " + "  ".join(_fmt(v) for v in s[r]), flush=True)
        print(
            f"    reciprocity {ident['reciprocity']:.3e}  sigma_max "
            f"{ident['sigma_max']:.12f}  spreads "
            + "  ".join(f"{k} {v * 100:.4f}%" for k, v in ident["spreads"].items()),
            flush=True,
        )
        eye = np.eye(s.shape[0])
        zoc = z0 * (eye + s) @ np.linalg.inv(eye - s)
        print("    open-circuit Z = 50 (I+S)(I-S)^-1 [ohm] (printed, never gated):",
              flush=True)
        for r in range(zoc.shape[0]):
            print("      " + "  ".join(_fmt(v) for v in zoc[r]), flush=True)
        for k, m in ident["means"].items():
            ref = ours_means[k]
            print(
                f"    class {k:8s} mean Palace {_fmt(m)} vs ans4 64 MHz s_matrix "
                f"{_fmt(ref)}: |dS|/|S| = {abs(m - ref) / abs(ref) * 100:.3f}% "
                f"(printed, never gated"
                + ("; confounded by port width f=1 vs f=0.5)" if name != "f05"
                   else ")"),
                flush=True,
            )

    # --- Negative control (predicted, rule (e)) — printed, never asserted --
    if "f1" in results and "flip" in results:
        s1, i1 = results["f1"]
        sf, iflip = results["flip"]
        d = np.diag([-1.0, 1.0, 1.0, 1.0])
        dev = float(np.linalg.norm(sf - d @ s1 @ d) / np.linalg.norm(s1))
        print(
            f"\n[OPS-61 1a W5] NEGATIVE CONTROL (predicted, not asserted): "
            f"||S_flip - D S_f1 D||/||S_f1|| predicted <= {FLIP_CONTROL_BAND:.0e}, "
            f"measured {dev:.3e} -> {'as predicted' if dev <= FLIP_CONTROL_BAND else 'NOT as predicted'}",
            flush=True,
        )
        print(
            "    flip spreads predicted adjacent/opposite ~200% (sign-flip ceiling), "
            "measured "
            + "  ".join(f"{k} {v * 100:.4f}%" for k, v in iflip["spreads"].items()),
            flush=True,
        )
        print(
            f"    flip reciprocity {iflip['reciprocity']:.3e} (f1 {i1['reciprocity']:.3e}); "
            f"sigma_max {iflip['sigma_max']:.12f} (f1 {i1['sigma_max']:.12f}) "
            f"— predicted unchanged",
            flush=True,
        )

    # --- Anchors (asserted) on W2's S ---------------------------------------
    assert "f1" in results, "no f1 port-S.csv: W2 did not produce S"
    _s1, i1 = results["f1"]
    fails = []
    if not i1["reciprocity"] <= RECIPROCITY_BAND:
        fails.append(f"reciprocity {i1['reciprocity']:.3e} > {RECIPROCITY_BAND}")
    if not i1["sigma_max"] <= 1.0 + PASSIVITY_SIGMA_TOLERANCE:
        fails.append(f"sigma_max {i1['sigma_max']:.12f} > 1 + {PASSIVITY_SIGMA_TOLERANCE}")
    for k, v in i1["spreads"].items():
        if not v <= ADJACENT_SPREAD_BAND:
            fails.append(f"{k} spread {v:.3e} > {ADJACENT_SPREAD_BAND}")
    print(
        f"\n[OPS-61 1a W5] f1 anchors: reciprocity <= {RECIPROCITY_BAND:.0e}, "
        f"sigma_max <= 1 + {PASSIVITY_SIGMA_TOLERANCE:.0e}, spreads <= "
        f"{ADJACENT_SPREAD_BAND}: {'PASS' if not fails else 'FAIL ' + '; '.join(fails)}",
        flush=True,
    )
    assert not fails, "; ".join(fails)


# ---------------------------------------------------------------------------
# `OPS-61` step 1c (F10b-i): the exact-strip flip control and its reader.
def flip_config():
    """Write ``flipx.json`` from ``f05x.json``: port 1 along ``-Z`` only."""
    src = json.loads((OUT_DIR / "f05x.json").read_text())
    assert src["Model"]["Mesh"] == STRIP_MESH_NAME, src["Model"]["Mesh"]
    ports = src["Boundaries"]["LumpedPort"]
    assert len(ports) == 4, f"{len(ports)} ports in f05x.json"
    for p in ports:
        i = int(p["Index"])
        assert p["Attributes"] == [NARROW_BASE + i], (i, p["Attributes"])
        assert p["Direction"] == "+Z", (i, p["Direction"])
    out = json.loads(json.dumps(src))
    for p in out["Boundaries"]["LumpedPort"]:
        if int(p["Index"]) == 1:
            p["Direction"] = "-Z"
    out["Problem"]["Output"] = "postpro_flipx"
    path = OUT_DIR / "flipx.json"
    path.write_text(json.dumps(out, indent=2) + "\n")
    print(f"[OPS-61 1c W1] wrote {path.relative_to(REPO)}", flush=True)
    print(json.dumps(out["Boundaries"]["LumpedPort"], indent=2), flush=True)


def read_exact_strip():
    """W4 of item 159: (P) asserted on ``f05x``; C1 (``flipx``) printed only."""
    from tests.validation.test_port_birdcage_four_port import (
        PASSIVITY_SIGMA_TOLERANCE,
        POOLED_SEPARATION_FLOOR,
        _circulant_classes,
        _class_spread,
    )
    from tests.validation.test_port_birdcage_lumped_column import ADJACENT_SPREAD_BAND
    from tests.validation.test_port_lumped_sheet_sweep import RECIPROCITY_BAND

    rung = json.loads(ANS4_METRICS.read_text())["rungs"]["64 MHz"]
    s_ours = np.array(
        [[complex(e["re"], e["im"]) for e in row] for row in rung["s_matrix"]]
    )
    ours_means = {k: complex(np.mean(v)) for k, v in _circulant_classes(s_ours).items()}

    results = {}
    for name in ("f05x", "flipx"):
        path = OUT_DIR / f"postpro_{name}" / "port-S.csv"
        if not path.exists():
            print(f"[OPS-61 1c W4] {name}: no {path.relative_to(REPO)}", flush=True)
            continue
        s, header = _read_port_s(path)
        ident = _identities(s)
        results[name] = (s, ident)
        print(f"\n[OPS-61 1c W4] {name}: header {header}", flush=True)
        print(
            f"    reciprocity {ident['reciprocity']:.3e}  sigma_max "
            f"{ident['sigma_max']:.12f}  spreads "
            + "  ".join(f"{k} {v * 100:.4f}%" for k, v in ident["spreads"].items()),
            flush=True,
        )
        if name == "f05x":
            for k, m in ident["means"].items():
                ref = ours_means[k]
                print(
                    f"    class {k:8s} Palace f05x vs ans4 64 MHz s_matrix (ours, "
                    f"ragged): |dS|/|S| = {abs(m - ref) / abs(ref) * 100:.3f}% "
                    "(printed, never gated)",
                    flush=True,
                )
            f1_csv = OUT_DIR / "postpro_f1" / "port-S.csv"
            if f1_csv.exists():
                i_f1 = _identities(_read_port_s(f1_csv)[0])
                for k, m in ident["means"].items():
                    ref = i_f1["means"][k]
                    print(
                        f"    class {k:8s} Palace f05x vs Palace f1: |dS|/|S| = "
                        f"{abs(m - ref) / abs(ref) * 100:.3f}% (printed, never "
                        "gated; width and mesh confounded)",
                        flush=True,
                    )

    # --- C1 negative control (predicted, rule (e)) — printed, never asserted
    if "f05x" in results and "flipx" in results:
        sx, ix = results["f05x"]
        sf, iflip = results["flipx"]
        d = np.diag([-1.0, 1.0, 1.0, 1.0])
        dev = float(np.linalg.norm(sf - d @ sx @ d) / np.linalg.norm(sx))
        print(
            f"\n[OPS-61 1c W4] C1 (predicted, not asserted): "
            f"||S_flipx - D S_f05x D||/||S_f05x|| predicted <= "
            f"{FLIP_CONTROL_BAND:.0e}, measured {dev:.3e} -> "
            f"{'as predicted' if dev <= FLIP_CONTROL_BAND else 'NOT as predicted'}",
            flush=True,
        )
        print(
            "    flipx spreads predicted adjacent/opposite ~200% (sign-flip "
            "ceiling), measured "
            + "  ".join(f"{k} {v * 100:.4f}%" for k, v in iflip["spreads"].items()),
            flush=True,
        )
        print(
            f"    flipx reciprocity {iflip['reciprocity']:.3e} (f05x "
            f"{ix['reciprocity']:.3e}); sigma_max {iflip['sigma_max']:.12f} "
            f"(f05x {ix['sigma_max']:.12f})",
            flush=True,
        )

    # --- (P) anchors (asserted) on f05x --------------------------------------
    assert "f05x" in results, "no f05x port-S.csv: W2 did not produce S"
    sx, ix = results["f05x"]
    classes = _circulant_classes(sx)
    column_power = np.sum(np.abs(sx) ** 2, axis=0)
    worst = max(ix["spreads"].values())
    pooled = _class_spread(np.concatenate([classes["adjacent"], classes["opposite"]]))
    fails = []
    if not ix["reciprocity"] <= RECIPROCITY_BAND:
        fails.append(f"reciprocity {ix['reciprocity']:.3e} > {RECIPROCITY_BAND}")
    if not ix["sigma_max"] <= 1.0 + PASSIVITY_SIGMA_TOLERANCE:
        fails.append(
            f"sigma_max {ix['sigma_max']:.12f} > 1 + {PASSIVITY_SIGMA_TOLERANCE}"
        )
    for k, v in enumerate(column_power, start=1):
        if not v <= 1.0 + 1.0e-9:
            fails.append(f"column {k} power {v:.12f} > 1 + 1e-9")
    for k, v in ix["spreads"].items():
        if not v <= ADJACENT_SPREAD_BAND:
            fails.append(f"{k} spread {v:.3e} > {ADJACENT_SPREAD_BAND}")
    if not pooled >= POOLED_SEPARATION_FLOOR * worst:
        fails.append(
            f"pooled {pooled:.3e} < {POOLED_SEPARATION_FLOOR} x worst {worst:.3e}"
        )
    print(
        f"\n[OPS-61 1c W4] f05x column power sums "
        + ", ".join(f"{v:.12f}" for v in column_power)
        + f"; pooled off-diagonal {pooled * 100:.4f}% vs floor "
        f"{POOLED_SEPARATION_FLOOR} x worst {worst * 100:.4f}%",
        flush=True,
    )
    print(
        f"[OPS-61 1c W4] (P) f05x anchors: reciprocity <= {RECIPROCITY_BAND:.0e}, "
        f"sigma_max <= 1 + {PASSIVITY_SIGMA_TOLERANCE:.0e}, column power <= 1 + "
        f"1e-9, spreads <= {ADJACENT_SPREAD_BAND}, pooled >= "
        f"{POOLED_SEPARATION_FLOOR} x worst: "
        f"{'PASS' if not fails else 'FAIL ' + '; '.join(fails)}",
        flush=True,
    )
    assert not fails, "; ".join(fails)


# ---------------------------------------------------------------------------
# `OPS-61` step 1d (§9 item 194, F17): Palace f05x at 10 / 128 MHz (order 1)
# and at 64 MHz order 2, from the same exported strip mesh.
FREQ_SET = {
    # name: (Freq [GHz], Order)
    "f05x_f10": (0.010, 1),
    "f05x_f128": (0.128, 1),
    "f05x_o2": (0.064, 2),
}

# `OPS-61` step 1e (§9 item 229, F21): Palace order 2 at 10 / 128 MHz.
FREQ_SET_O2 = {
    # name: (Freq [GHz], Order)
    "f05x_o2_f10": (0.010, 2),
    "f05x_o2_f128": (0.128, 2),
}

# `OPS-61` step 2a (§9 item 234, F22): the exact-rectangle hole mesh + PEC column.
HOLE_MESH_NAME = "ops61_ans6_x1_hole_strip.msh"

# `OPS-61` step 1f (§9 item 418, F34): Palace order 2 on `ANS-4`'s exact-strip
# h = 0.0075 rung at 10 MHz, one excitation per Palace command.
H0075_RES = 0.0075
H0075_MESH_NAME = "ops61_ans4_h0075_strip.msh"
# Literal, provenance `docs/testing/logs/20261005T070009Z_ANS-4-step6b.log:2167,
# :3700` (xl entry 27's degree-2 run); no tracked constant holds this count.
H0075_CELL_RECORD = 282198
H0075_SET = {f"f05x_h0075_o2_f10_e{k}": k for k in (1, 2, 3, 4)}
H0075_OURS_LOG = (
    REPO / "docs" / "testing" / "logs" / "20261005T070009Z_ANS-4-step6b.log"
)
H0075_OURS_HEADER = "[ANS-4 step5] 10 MHz kappa-off S (10 digits):"
H0075_OURS_SIGMA_MAX = 0.999995526884  # printed at the log's :3707 (Tie target)
H0075_TIE_BAND = 1.0e-9  # parse check, not a physics band (brief, item 418)
# Ours' printed h-move ×1 -> 0.0075, degree 2, 10 MHz (log :3726), printed beside (C-h).
H0075_OURS_H_MOVE = {"self": 2.527e-3, "adjacent": 4.47e-4, "opposite": 1.375e-3}

# `OPS-61` steps 1g / 1h (§9 items 434 / 435, F34): the same rung at 64 / 128 MHz.
# Shared code written once by whichever item ran first (item 434); additive.
H0075_FREQS = {"64": 0.064, "128": 0.128}  # GHz
H0075_FREQ_SET = {
    f: {f"f05x_h0075_o2_f{f}_e{k}": k for k in (1, 2, 3, 4)} for f in H0075_FREQS
}
H0075_OURS_BLOCKS = {
    # header at 20261005T070009Z_ANS-4-step6b.log:3708 (S lines :3709-3712, :3709
    # PASSED-prefixed); sigma_max printed :3713 (and (G) :3851); h-move :3727;
    # Palace x1 order 2 at 64 MHz = postpro_f05x_o2 (item 194).
    "64": {
        "header": "[ANS-4 step5] 64 MHz kappa-off S (10 digits):",
        "sigma_max": 0.999823496773,
        "h_move": {"self": 2.345e-3, "adjacent": 9.50e-4, "opposite": 2.098e-3},
        "h_move_line": ":3727",
        "pal_x1": "postpro_f05x_o2",
        "proj": (189.8, 193.7),  # item 418 W6's projected single excitation, s
    },
    # header :3714 (S lines :3715-3718); sigma_max printed on the pytest-prefixed
    # :3720 (and (G) :3859); h-move :3728; Palace x1 order 2 at 128 MHz =
    # postpro_f05x_o2_f128 (item 229).
    "128": {
        "header": "[ANS-4 step5] 128 MHz kappa-off S (10 digits):",
        "sigma_max": 0.999320852426,
        "h_move": {"self": 2.204e-3, "adjacent": 1.590e-3, "opposite": 2.863e-3},
        "h_move_line": ":3728",
        "pal_x1": "postpro_f05x_o2_f128",
        "proj": (265.9, 271.3),
    },
}
H0075_ITEM = {"64": 434, "128": 435}


# `OPS-61` step 1i (§9 item 437, F36, provisional): Palace order 2 on `ANS-4`'s
# exact-strip h = 0.005 rung at 10 MHz (xl entry 30), one excitation per command.
H0050_RES = 0.005
H0050_MESH_NAME = "ops61_ans4_h0050_strip.msh"
# Literal, provenance `docs/testing/logs/20261011T070008Z_ANS-4-step6c.log:2280,
# :3814` (xl entry 30's degree-2 run); no tracked constant holds this count.
H0050_CELL_RECORD = 593478
H0050_SET = {f"f05x_h0050_o2_f10_e{k}": k for k in (1, 2, 3, 4)}
H0050_OURS_LOG = (
    REPO / "docs" / "testing" / "logs" / "20261011T070008Z_ANS-4-step6c.log"
)
H0050_OURS_HEADER = "[ANS-4 step5] 10 MHz kappa-off S (10 digits):"  # :3815
H0050_OURS_SIGMA_MAX = 0.999995535227  # printed at the log's :3820 (Tie target)
# Ours' printed move x1 -> 0.005, degree 2, 10 MHz (log :3840), printed beside (C2-h).
H0050_OURS_X1_MOVE = {"self": 2.519e-3, "adjacent": 3.77e-4, "opposite": 1.528e-3}

# `OPS-61` step 1j (§9 item 458 / 458b, §10 F39): Palace order 2 on the same
# h = 0.005 exact-strip rung at 64 MHz.  Additive; nothing above is edited.
H0050_F64_SET = {f"f05x_h0050_o2_f64_e{k}": k for k in (1, 2, 3, 4)}
# Header at 20261011T070008Z_ANS-4-step6c.log:3821; rows :3822, :3823, :3825, :3827
# (blank :3824 / :3826, pytest node-ID prefixes on :3825 / :3827).
H0050_F64_OURS_HEADER = "[ANS-4 step5] 64 MHz kappa-off S (10 digits):"
H0050_F64_OURS_SIGMA_MAX = 0.999823581365  # printed at :3828 (and (G) :3907); (Tie) target
# Ours' printed move x1 -> 0.005, degree 2, 64 MHz (log :3841), printed beside (C2-h).
H0050_F64_OURS_X1_MOVE = {"self": 2.868e-3, "adjacent": 1.144e-3, "opposite": 2.699e-3}
H0050_F64_PROJ = (380.0, 435.0)  # item 458 brief's bracket, s (item 437 W6 printed 392.1-429.3)


def _hole_wall_attr():
    from fem_em_solver.io.mesh import BIRDCAGE_CONDUCTOR_SURFACE_TAG

    return int(BIRDCAGE_CONDUCTOR_SURFACE_TAG)


def _hole_set():
    """name -> (PEC attributes, Order); the wall attribute imported, never restated."""
    wall = _hole_wall_attr()
    return {
        "h05x_pec": ([EXTERIOR_TAG, wall], 1),
        "h05x_pmc": ([EXTERIOR_TAG], 1),
        "h05x_pec_o2": ([EXTERIOR_TAG, wall], 2),
    }


def freq_configs():
    """Write ``f05x_f10`` / ``f05x_f128`` / ``f05x_o2`` from ``f05x.json`` only."""
    src = json.loads((OUT_DIR / "f05x.json").read_text())
    assert src["Model"]["Mesh"] == STRIP_MESH_NAME, src["Model"]["Mesh"]
    ports = src["Boundaries"]["LumpedPort"]
    assert len(ports) == 4, f"{len(ports)} ports in f05x.json"
    for p in ports:
        i = int(p["Index"])
        assert p["Attributes"] == [NARROW_BASE + i], (i, p["Attributes"])
        assert p["Direction"] == "+Z", (i, p["Direction"])
    assert src["Solver"]["Order"] == 1, src["Solver"]["Order"]
    samples = src["Solver"]["Driven"]["Samples"]
    assert len(samples) == 1 and samples[0]["Freq"] == [0.064], samples
    for name, (freq_ghz, order) in FREQ_SET.items():
        out = json.loads(json.dumps(src))
        out["Solver"]["Order"] = order
        out["Solver"]["Driven"]["Samples"][0]["Freq"] = [freq_ghz]
        out["Problem"]["Output"] = f"postpro_{name}"
        path = OUT_DIR / f"{name}.json"
        assert not path.exists(), f"{path.relative_to(REPO)} exists: never overwrite"
        path.write_text(json.dumps(out, indent=2) + "\n")
        print(
            f"[OPS-61 1d W1] wrote {path.relative_to(REPO)}: Order {order}, "
            f"Freq [{freq_ghz}] GHz, Output postpro_{name}",
            flush=True,
        )


def _p_check(name, path, tag):
    """Item 159's (P) set on one Palace ``port-S.csv`` (lifted, §9 item 229).

    Prints with the caller's ``tag``; returns ``(s, ident, fails)``.
    """
    from tests.validation.test_port_birdcage_four_port import (
        PASSIVITY_SIGMA_TOLERANCE,
        POOLED_SEPARATION_FLOOR,
        _circulant_classes,
        _class_spread,
    )
    from tests.validation.test_port_birdcage_lumped_column import ADJACENT_SPREAD_BAND
    from tests.validation.test_port_lumped_sheet_sweep import RECIPROCITY_BAND

    def line(d):
        return "  ".join(f"{k} {v * 100:.4f}%" for k, v in d.items())

    s, header = _read_port_s(path)
    ident = _identities(s)
    classes = _circulant_classes(s)
    column_power = np.sum(np.abs(s) ** 2, axis=0)
    worst = max(ident["spreads"].values())
    pooled = _class_spread(
        np.concatenate([classes["adjacent"], classes["opposite"]])
    )
    f = []
    if not ident["reciprocity"] <= RECIPROCITY_BAND:
        f.append(f"reciprocity {ident['reciprocity']:.3e} > {RECIPROCITY_BAND}")
    if not ident["sigma_max"] <= 1.0 + PASSIVITY_SIGMA_TOLERANCE:
        f.append(
            f"sigma_max {ident['sigma_max']:.12f} > 1 + {PASSIVITY_SIGMA_TOLERANCE}"
        )
    for k, v in enumerate(column_power, start=1):
        if not v <= 1.0 + 1.0e-9:
            f.append(f"column {k} power {v:.12f} > 1 + 1e-9")
    for k, v in ident["spreads"].items():
        if not v <= ADJACENT_SPREAD_BAND:
            f.append(f"{k} spread {v:.3e} > {ADJACENT_SPREAD_BAND}")
    if not pooled >= POOLED_SEPARATION_FLOOR * worst:
        f.append(
            f"pooled {pooled:.3e} < {POOLED_SEPARATION_FLOOR} x worst {worst:.3e}"
        )
    print(
        f"\n{tag} {name}: header {header}\n"
        f"    reciprocity {ident['reciprocity']:.3e}  sigma_max "
        f"{ident['sigma_max']:.12f}  spreads {line(ident['spreads'])}\n"
        f"    column power sums "
        + ", ".join(f"{v:.12f}" for v in column_power)
        + f"; pooled off-diagonal {pooled * 100:.4f}% vs floor "
        f"{POOLED_SEPARATION_FLOOR} x worst {worst * 100:.4f}% "
        f"(separation {pooled / worst:.2f}x)",
        flush=True,
    )
    print(
        f"{tag} (P) {name} anchors: reciprocity <= "
        f"{RECIPROCITY_BAND:.0e}, sigma_max <= 1 + {PASSIVITY_SIGMA_TOLERANCE:.0e}, "
        f"column power <= 1 + 1e-9, spreads <= {ADJACENT_SPREAD_BAND}, pooled >= "
        f"{POOLED_SEPARATION_FLOOR} x worst: "
        f"{'PASS' if not f else 'FAIL ' + '; '.join(f)}",
        flush=True,
    )
    return s, ident, f


def read_exact_strip_set():
    """W5 of item 194: (P) asserted on each present run; o2 readings printed."""
    from tests.validation.test_port_birdcage_four_port import (
        PASSIVITY_SIGMA_TOLERANCE,
        POOLED_SEPARATION_FLOOR,
        _circulant_classes,
        _class_spread,
    )
    from tests.validation.test_port_birdcage_lumped_column import ADJACENT_SPREAD_BAND
    from tests.validation.test_port_lumped_sheet_sweep import RECIPROCITY_BAND

    def rel(a, b):
        return {k: abs(a[k] - b[k]) / abs(b[k]) for k in b}

    def line(d):
        return "  ".join(f"{k} {v * 100:.4f}%" for k, v in d.items())

    def ours_means(name):
        path = OUT_DIR / name
        if not path.exists():
            return None
        m = np.array(
            [[complex(a, b) for a, b in row]
             for row in json.loads(path.read_text())["s_kappa_off"]]
        )
        return {k: complex(np.mean(v)) for k, v in _circulant_classes(m).items()}

    results, fails = {}, []
    for name in FREQ_SET:
        path = OUT_DIR / f"postpro_{name}" / "port-S.csv"
        if not path.exists():
            print(f"[OPS-61 1d W5] {name}: no {path.relative_to(REPO)} "
                  "(not measured)", flush=True)
            continue
        s, ident, f = _p_check(name, path, "[OPS-61 1d W5]")
        results[name] = (s, ident)
        fails += [f"{name}: {x}" for x in f]

    # --- Order-2 readings (printed, never gated) -----------------------------
    if "f05x_o2" in results:
        o2 = results["f05x_o2"][1]["means"]
        o1_csv = OUT_DIR / "postpro_f05x" / "port-S.csv"
        if o1_csv.exists():
            o1 = _identities(_read_port_s(o1_csv)[0])["means"]
            print(
                f"\n[OPS-61 1d W5] Palace o2 vs Palace o1 (f05x), 64 MHz: |d mean|/|mean| "
                f"{line(rel(o2, o1))} (printed, never gated; predicted ~6.7 / 2.5 / 4.5 %)",
                flush=True,
            )
        p2 = ours_means("ours_f05x_p2.json")
        p1 = ours_means("ours_f05x_p1.json")
        if p2 is None:
            print("[OPS-61 1d W5] ours p2: absent", flush=True)
        else:
            g2 = rel(p2, o2)
            print(
                f"[OPS-61 1d W5] ours p2 (kappa off) vs Palace o2, matched order: "
                f"|d mean|/|mean| {line(g2)} (printed, never gated; predicted <= 1e-3)",
                flush=True,
            )
            if p1 is not None:
                g1 = rel(p1, o2)
                print(
                    f"[OPS-61 1d W5] C2' ours p1 (kappa off) vs Palace o2, order "
                    f"mismatch: {line(g1)} (predicted ~6.67 / 2.54 / 4.52 %)",
                    flush=True,
                )
                print(
                    "[OPS-61 1d W5] C2' separation |ours p1 - o2| / |ours p2 - o2| "
                    + "  ".join(
                        f"{k} {g1[k] / g2[k]:.4g}x" if g2[k] > 0 else f"{k} inf"
                        for k in g1
                    )
                    + " (predicted >= 25x)",
                    flush=True,
                )
            else:
                print("[OPS-61 1d W5] ours p1: absent", flush=True)
    else:
        print("\n[OPS-61 1d W5] Palace o2: not measured", flush=True)

    assert results, "no item-194 port-S.csv present: W2-W4 produced no S"
    assert not fails, "; ".join(fails)
    print(f"[OPS-61 1d W5] (P) PASS on {sorted(results)}", flush=True)


def o2_configs():
    """Item 229 W1: write ``f05x_o2_f10`` / ``f05x_o2_f128`` from ``f05x.json`` only."""
    src = json.loads((OUT_DIR / "f05x.json").read_text())
    assert src["Model"]["Mesh"] == STRIP_MESH_NAME, src["Model"]["Mesh"]
    ports = src["Boundaries"]["LumpedPort"]
    assert len(ports) == 4, f"{len(ports)} ports in f05x.json"
    for p in ports:
        i = int(p["Index"])
        assert p["Attributes"] == [NARROW_BASE + i], (i, p["Attributes"])
        assert p["Direction"] == "+Z", (i, p["Direction"])
    assert src["Solver"]["Order"] == 1, src["Solver"]["Order"]
    samples = src["Solver"]["Driven"]["Samples"]
    assert len(samples) == 1 and samples[0]["Freq"] == [0.064], samples
    for name in FREQ_SET_O2:
        path = OUT_DIR / f"{name}.json"
        assert not path.exists(), f"{path.relative_to(REPO)} exists: never overwrite"
    for name, (freq_ghz, order) in FREQ_SET_O2.items():
        out = json.loads(json.dumps(src))
        out["Solver"]["Order"] = order
        out["Solver"]["Driven"]["Samples"][0]["Freq"] = [freq_ghz]
        out["Problem"]["Output"] = f"postpro_{name}"
        path = OUT_DIR / f"{name}.json"
        path.write_text(json.dumps(out, indent=2) + "\n")
        print(
            f"[OPS-61 1e W1] wrote {path.relative_to(REPO)}: Order {order}, "
            f"Freq [{freq_ghz}] GHz, Output postpro_{name}",
            flush=True,
        )


def read_exact_strip_o2():
    """Item 229 W4: (P) on Palace o2 at 10 / 128 MHz; ours p2 record vs Palace o2."""
    from tests.validation.test_ans4_exact_strip_records import (
        ANS4_EXACT_STRIP_S_P1_RECORDS,
        ANS4_EXACT_STRIP_S_P2_X1,
    )
    from tests.validation.test_port_birdcage_four_port import _circulant_classes

    tag = "[OPS-61 1e W4]"
    TIE_BAND = 1.0e-6  # (T), predicted; a miss stops the 10 / 128 reading

    def rel(a, b):
        return {k: abs(a[k] - b[k]) / abs(b[k]) for k in b}

    def line(d):
        return "  ".join(f"{k} {v * 100:.4f}%" for k, v in d.items())

    def sci(d):
        return "  ".join(f"{k} {v:.3e}" for k, v in d.items())

    def means(m):
        return {k: complex(np.mean(v)) for k, v in _circulant_classes(m).items()}

    # --- (T) the 64 MHz tie: record vs Palace o2 (f05x_o2) --------------------
    tie_csv = OUT_DIR / "postpro_f05x_o2" / "port-S.csv"
    if not tie_csv.exists():
        print(f"{tag} (T) no {tie_csv.relative_to(REPO)}: not measured", flush=True)
        tie_ok = True
    else:
        s64, ident64, f64 = _p_check("f05x_o2", tie_csv, tag)
        assert not f64, "f05x_o2: " + "; ".join(f64)
        t = rel(means(ANS4_EXACT_STRIP_S_P2_X1["64"]), ident64["means"])
        tie_ok = max(t.values()) <= TIE_BAND
        print(
            f"{tag} (T) ours p2 record (ANS4_EXACT_STRIP_S_P2_X1['64']) vs Palace o2 "
            f"f05x_o2, |d mean|/|mean| {sci(t)} (predicted <= {TIE_BAND:.0e}): "
            f"{'GREEN' if tie_ok else 'MISS'}",
            flush=True,
        )
    assert tie_ok, "(T) missed: the record is not Palace's discrete system; stop"

    results, fails = {}, []
    for name, (freq_ghz, _order) in FREQ_SET_O2.items():
        path = OUT_DIR / f"postpro_{name}" / "port-S.csv"
        if not path.exists():
            print(f"{tag} {name}: no {path.relative_to(REPO)} (not measured)",
                  flush=True)
            continue
        s, ident, f = _p_check(name, path, tag)
        results[name] = (s, ident)
        fails += [f"{name}: {x}" for x in f]
        mhz = str(int(round(freq_ghz * 1000)))
        o2 = ident["means"]
        p2 = means(ANS4_EXACT_STRIP_S_P2_X1[mhz])
        p1 = means(ANS4_EXACT_STRIP_S_P1_RECORDS[mhz])
        g2 = rel(p2, o2)
        g1 = rel(p1, o2)
        print(
            f"{tag} {mhz} MHz ours p2 record (kappa off) vs Palace o2, matched "
            f"order: |d mean|/|mean| {line(g2)} = {sci(g2)} (printed, never gated; "
            f"predicted <= 1e-6)",
            flush=True,
        )
        print(
            f"{tag} {mhz} MHz C2'' ours p1 record vs Palace o2, order mismatch: "
            f"{line(g1)} (predicted = ours p2-vs-p1 move)",
            flush=True,
        )
        print(
            f"{tag} {mhz} MHz C2'' separation |ours p1 - o2| / |ours p2 - o2| "
            + "  ".join(
                f"{k} {g1[k] / g2[k]:.4g}x" if g2[k] > 0 else f"{k} inf" for k in g1
            )
            + " (predicted >= 25x)",
            flush=True,
        )
        print(
            f"{tag} {mhz} MHz ours p2 record vs ours p1 record: "
            f"{line(rel(p2, p1))}",
            flush=True,
        )
        o1_csv = OUT_DIR / f"postpro_f05x_f{mhz}" / "port-S.csv"
        if o1_csv.exists():
            o1 = _identities(_read_port_s(o1_csv)[0])["means"]
            print(
                f"{tag} {mhz} MHz Palace o2 vs Palace o1 (f05x_f{mhz}): "
                f"{line(rel(o2, o1))} (printed, never gated)",
                flush=True,
            )
        else:
            print(f"{tag} {mhz} MHz Palace o1: absent", flush=True)

    assert results, "no item-229 port-S.csv present: W2/W3 produced no S"
    assert not fails, "; ".join(fails)
    print(f"{tag} (P) PASS on {sorted(results)}", flush=True)


def export_hole_pec():
    """Item 234 W1: the exact-rectangle hole mesh (Gmsh 2.2) + three Palace configs.

    Box -> attribute 1, cavity wall (tag 401) -> its own attribute, N_i -> 20 + i.
    Never overwrites; ``export()`` is not called.  Anchors (E1)-(E5) after writing.
    """
    import dolfinx
    from mpi4py import MPI

    from tests.mesh.helpers import global_cell_tag_set
    from tests.mesh.test_birdcage_port_sheets import SHEET_IFACE
    from tests.mesh.test_birdcage_port_sheets import _build as _sheets_build
    from tests.mesh.test_birdcage_port_tags import LEG_COUNT
    from tests.validation import test_th14_birdcage_copper as th14
    from tests.validation.test_ans4_exact_strip_records import (
        SPEC_AREA_M2,
        SPEC_HEIGHT_M,
        SPEC_WIDTH_M,
    )
    from tests.validation.test_ans6_exact_strip import ANS6_EXACT_STRIP_HOLE_CELL_RECORD
    from tests.validation.test_lossy_sphere_fullwave import (
        SALINE_EPSILON_R,
        SALINE_SIGMA,
    )
    from tests.validation.test_ops61_exact_strip_ports import STRIP_BAND
    from tests.validation.test_port_birdcage_lumped_column import (
        CONDUCTOR_CELL_TAG,
        PHANTOM_CELL_TAG,
    )
    from tests.validation.test_port_lumped_narrowed_sheet import GATED_WIDTH_FRACTION

    tag = "[OPS-61 2a W1]"
    comm = MPI.COMM_WORLD
    assert comm.size == 1, "export_hole_pec is a serial writer: run it at mpiexec -n 1"
    t_start = time.perf_counter()
    wall_attr = _hole_wall_attr()
    hole_set = _hole_set()
    ports_idx = list(range(1, LEG_COUNT + 1))

    # 1. never overwrite; f05x.json checked as o2_configs does.
    mesh_path = OUT_DIR / HOLE_MESH_NAME
    assert not mesh_path.exists(), f"{mesh_path.relative_to(REPO)} exists: never overwrite"
    for name in hole_set:
        path = OUT_DIR / f"{name}.json"
        assert not path.exists(), f"{path.relative_to(REPO)} exists: never overwrite"
    src = json.loads((OUT_DIR / "f05x.json").read_text())
    assert src["Model"]["Mesh"] == STRIP_MESH_NAME, src["Model"]["Mesh"]
    ports = src["Boundaries"]["LumpedPort"]
    assert len(ports) == 4, f"{len(ports)} ports in f05x.json"
    for p in ports:
        i = int(p["Index"])
        assert p["Attributes"] == [NARROW_BASE + i], (i, p["Attributes"])
        assert p["Direction"] == "+Z", (i, p["Direction"])
    assert src["Solver"]["Order"] == 1, src["Solver"]["Order"]
    samples = src["Solver"]["Driven"]["Samples"]
    assert len(samples) == 1 and samples[0]["Freq"] == [0.064], samples

    # 2. the hole route with exact rectangles, as th14._build_ladder builds it.
    msh, cell_tags, facet_tags, diag, t_mesh = _sheets_build(
        True, as_hole=True, leg_sheet_strip_fraction=float(GATED_WIDTH_FRACTION)
    )
    tdim = msh.topology.dim
    fdim = tdim - 1
    msh.topology.create_connectivity(fdim, tdim)
    msh.topology.create_entity_permutations()
    ncells = int(msh.topology.index_map(tdim).size_global)
    print(
        f"{tag} hole mesh built in {t_mesh:.2f} s: {ncells} tets (record "
        f"{ANS6_EXACT_STRIP_HOLE_CELL_RECORD})",
        flush=True,
    )

    # 3. box / wall split and the ports our PEC column was solved with.
    _outer, census = th14._outer_box_tags(msh, facet_tags, comm)
    tags_f, _defs, specs = th14._build_ports(msh, cell_tags, comm)
    print(f"{tag} census {census}", flush=True)

    x = np.asarray(msh.geometry.x, dtype=float)

    def facet_nodes(facets):
        facets = np.asarray(facets, dtype=np.int32)
        if not facets.size:
            return np.zeros((0, 3), dtype=np.int64)
        return np.asarray(
            dolfinx.cpp.mesh.entities_to_geometry(msh._cpp_object, fdim, facets, False)
        ).reshape(-1, 3)

    # 4. triangle blocks: box -> 1, wall -> 401, N_i -> 20 + i.
    ext = np.asarray(dolfinx.mesh.exterior_facet_indices(msh.topology), dtype=np.int32)
    wall = np.asarray(facet_tags.find(wall_attr), dtype=np.int32)
    box = np.setdiff1d(ext, wall).astype(np.int32)
    box_nodes, wall_nodes = facet_nodes(box), facet_nodes(wall)
    box_area = float(np.sum(_tri_areas(x, box_nodes)))
    wall_area = float(np.sum(_tri_areas(x, wall_nodes)))
    tri_blocks = [(box_nodes, EXTERIOR_TAG), (wall_nodes, wall_attr)]
    port_rows = {}
    for i, spec in zip(ports_idx, specs):
        ptag = SHEET_IFACE + i
        assert int(spec.facet_tag) == ptag, (i, spec.facet_tag)
        narrow = np.asarray(tags_f.find(ptag), dtype=np.int32)
        n_nodes = facet_nodes(narrow)
        a_n = float(np.sum(_tri_areas(x, n_nodes)))
        port_rows[i] = {
            "facets": int(narrow.size),
            "area": a_n,
            "width": float(spec.sheet_width_m),
            "height": float(spec.gap_height_m),
        }
        tri_blocks.append((n_nodes, NARROW_BASE + i))

    OUT_DIR.mkdir(parents=True, exist_ok=True)
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
    total_vol = float(np.sum(np.abs(vol6)) / 6.0)

    ntri = sum(b.shape[0] for b, _ in tri_blocks)
    with open(mesh_path, "w") as fh:
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

    lo, hi = x.min(axis=0), x.max(axis=0)
    tags = sorted(global_cell_tag_set(msh, cell_tags))
    counts = {t: int(np.count_nonzero(ctag == t)) for t in tags}
    print(
        f"{tag} wrote {mesh_path.relative_to(REPO)}: {x.shape[0]} nodes, {ncl} tets "
        f"({int(neg.sum())} reoriented), {ntri} tris; volume {total_vol:.9e} m^3; "
        f"bbox lo {lo} hi {hi}",
        flush=True,
    )
    print(f"{tag} cell tags (count): {counts}", flush=True)
    cad = diag.get("conductor_cavity_cad", {}) if isinstance(diag, dict) else {}
    cad_area = cad.get("area_m2") if isinstance(cad, dict) else None
    ratio = f"{wall_area / float(cad_area):.9f}" if cad_area else "n/a"
    print(
        f"{tag} box (attr {EXTERIOR_TAG}): {box.size} tris area {box_area:.12e} m^2; "
        f"wall (attr {wall_attr}): {wall.size} tris area {wall_area:.12e} m^2, ratio "
        f"to conductor_cavity_cad area_m2 {cad_area}: {ratio} (printed only)",
        flush=True,
    )
    for i, r in port_rows.items():
        print(
            f"    P{i}: N_{i} {r['facets']} tris (attr {NARROW_BASE + i}) area "
            f"{r['area']:.12e} m^2 (rel {r['area'] / SPEC_AREA_M2 - 1.0:+.3e}); spec "
            f"sheet_width_m {r['width']:.12e} (rel {r['width'] / SPEC_WIDTH_M - 1.0:+.3e}) "
            f"gap_height_m {r['height']:.12e} (rel {r['height'] / SPEC_HEIGHT_M - 1.0:+.3e})",
            flush=True,
        )

    # 5. three configs from the f05x template.
    assert SALINE_SIGMA == 0.5 and SALINE_EPSILON_R == 78.0
    air = [t for t in tags if t not in (CONDUCTOR_CELL_TAG, PHANTOM_CELL_TAG)]
    materials = [
        {
            "Attributes": [PHANTOM_CELL_TAG],
            "Permeability": 1.0,
            "Permittivity": SALINE_EPSILON_R,
            "Conductivity": SALINE_SIGMA,
        },
        {"Attributes": air, "Permeability": 1.0, "Permittivity": 1.0},
    ]
    for name, (pec_attrs, order) in hole_set.items():
        out = json.loads(json.dumps(src))
        out["Model"]["Mesh"] = HOLE_MESH_NAME
        out["Domains"]["Materials"] = materials
        out["Boundaries"]["PEC"]["Attributes"] = list(pec_attrs)
        out["Solver"]["Order"] = order
        out["Problem"]["Output"] = f"postpro_{name}"
        path = OUT_DIR / f"{name}.json"
        path.write_text(json.dumps(out, indent=2) + "\n")
        print(
            f"{tag} wrote {path.relative_to(REPO)}: Order {order}\n"
            f"  Boundaries {json.dumps(out['Boundaries'])}\n"
            f"  Domains {json.dumps(out['Domains'])}",
            flush=True,
        )
    print(f"{tag} export total {time.perf_counter() - t_start:.2f} s (mesh {t_mesh:.2f} s)",
          flush=True)

    # 6. anchors (E1)-(E5), after writing.
    fails = []
    if ncl != ANS6_EXACT_STRIP_HOLE_CELL_RECORD or ncells != ncl:
        fails.append(f"(E1) tets {ncl} != {ANS6_EXACT_STRIP_HOLE_CELL_RECORD}")
    if CONDUCTOR_CELL_TAG in tags or PHANTOM_CELL_TAG not in tags:
        fails.append(f"(E2) cell tags {tags}")
    if not (census["cavity_exterior"] == census["cavity"] > 0
            and census["outer"] == census["exterior"] - census["cavity"]
            and census["tag_collision"] == 0):
        fails.append(f"(E3) census {census}")
    if not abs(box_area - EXTERIOR_AREA_M2) <= EXTERIOR_AREA_BAND:
        fails.append(f"(E4) box area {box_area:.12e} vs {EXTERIOR_AREA_M2}")
    for i, r in port_rows.items():
        for k, (got, want) in {"area": (r["area"], SPEC_AREA_M2),
                               "width": (r["width"], SPEC_WIDTH_M),
                               "height": (r["height"], SPEC_HEIGHT_M)}.items():
            if not abs(got / want - 1.0) < STRIP_BAND:
                fails.append(f"(E5) P{i} {k} {got:.12e} vs {want:.12e}")
    print(
        f"{tag} anchors (E1)-(E5) [box band {EXTERIOR_AREA_BAND:.0e}, strip band "
        f"{STRIP_BAND}]: {'PASS' if not fails else 'FAIL ' + '; '.join(fails)}",
        flush=True,
    )
    assert not fails, "; ".join(fails)


def read_hole_pec():
    """Item 234 W5: (P) on h05x_pec (+o2); (K') floor; ours-vs-Palace printed."""
    from tests.validation.test_ans4_exact_strip_records import (
        ANS4_EXACT_STRIP_S_P1_RECORDS,
    )
    from tests.validation.test_ans6_exact_strip import (
        ANS6_EXACT_STRIP_CLASS_MEANS_P1,
        ANS6_EXACT_STRIP_CLASS_MEANS_P2,
    )
    from tests.validation.test_port21_feed_sensitivity import KNOB_REACHES_THE_SOLVE
    from tests.validation.test_port_birdcage_four_port import _circulant_classes

    tag = "[OPS-61 2a W5]"

    def rel(a, b):
        return {k: abs(a[k] - b[k]) / abs(b[k]) for k in b}

    def line(d):
        return "  ".join(f"{k} {v * 100:.4f}%" for k, v in d.items())

    def sci(d):
        return "  ".join(f"{k} {v:.3e}" for k, v in d.items())

    def means(m):
        return {k: complex(np.mean(v)) for k, v in _circulant_classes(np.asarray(m)).items()}

    def csv_of(name):
        return OUT_DIR / f"postpro_{name}" / "port-S.csv"

    p1 = dict(ANS6_EXACT_STRIP_CLASS_MEANS_P1["64"]["pec"])
    p2 = dict(ANS6_EXACT_STRIP_CLASS_MEANS_P2["64"]["pec"])

    fails = []
    assert csv_of("h05x_pec").exists(), f"no {csv_of('h05x_pec').relative_to(REPO)}"
    s_pec, id_pec, f = _p_check("h05x_pec", csv_of("h05x_pec"), tag)
    fails += [f"h05x_pec: {x}" for x in f]
    o1 = id_pec["means"]
    o2 = None
    if csv_of("h05x_pec_o2").exists():
        _s2, id_o2, f2 = _p_check("h05x_pec_o2", csv_of("h05x_pec_o2"), tag)
        fails += [f"h05x_pec_o2: {x}" for x in f2]
        o2 = id_o2["means"]
    else:
        print(f"{tag} h05x_pec_o2: not measured", flush=True)

    # --- (T1) / (T2) ours vs Palace, matched order (printed) -----------------
    t1 = rel(p1, o1)
    print(
        f"\n{tag} (T1) ours ANS6_EXACT_STRIP_CLASS_MEANS_P1['64']['pec'] vs Palace "
        f"h05x_pec (o1): |d mean|/|mean| {line(t1)} = {sci(t1)} (printed; predicted <= 1e-6)",
        flush=True,
    )
    t2 = None
    if o2 is not None:
        t2 = rel(p2, o2)
        print(
            f"{tag} (T2) ours ANS6_EXACT_STRIP_CLASS_MEANS_P2['64']['pec'] vs Palace "
            f"h05x_pec_o2: |d mean|/|mean| {line(t2)} = {sci(t2)} (printed; predicted <= 1e-6)",
            flush=True,
        )
    else:
        print(f"{tag} (T2) not measured", flush=True)

    # --- (C3) conductor-model move as each code sees it ---------------------
    f05x_csv = OUT_DIR / "postpro_f05x" / "port-S.csv"
    ours_solid = means(ANS4_EXACT_STRIP_S_P1_RECORDS["64"])
    ours_move = rel(p1, ours_solid)
    if f05x_csv.exists():
        pal_solid = _identities(_read_port_s(f05x_csv)[0])["means"]
        pal_move = rel(o1, pal_solid)
        print(
            f"{tag} (C3) Palace h05x_pec vs Palace f05x (sigma 800 solid): {line(pal_move)}\n"
            f"{tag} (C3) ours P1 pec hole vs ours ANS4_EXACT_STRIP_S_P1_RECORDS['64']: "
            f"{line(ours_move)} (predicted equal to the printed digit)",
            flush=True,
        )
    else:
        print(f"{tag} (C3) postpro_f05x: absent; ours move {line(ours_move)}", flush=True)

    # --- (K') PMC-wall control: floor asserted, size predicted ---------------
    if csv_of("h05x_pmc").exists():
        s_pmc, id_pmc, fk = _p_check("h05x_pmc", csv_of("h05x_pmc"), tag + " [control]")
        print(f"{tag} h05x_pmc (P) printed only (control, never asserted): "
              f"{'clean' if not fk else '; '.join(fk)}", flush=True)
        k = rel(id_pmc["means"], o1)
        k_ok = min(k.values()) > KNOB_REACHES_THE_SOLVE
        print(
            f"{tag} (K') h05x_pmc vs h05x_pec class move {line(k)} = {sci(k)}; floor "
            f"> {KNOB_REACHES_THE_SOLVE:.0e} on every class: {'PASS' if k_ok else 'FAIL'}; "
            f"predicted >= 10 % on at least one class: "
            f"{'met' if max(k.values()) >= 0.10 else 'NOT met'}",
            flush=True,
        )
        if not k_ok:
            fails.append(f"(K') min class move {min(k.values()):.3e}")
    else:
        print(f"{tag} (K') h05x_pmc: not measured", flush=True)

    # --- (C2''') order-mismatch separation (printed) -------------------------
    g_num = rel(p2, o1)
    print(
        f"{tag} (C2''') ours p2 vs Palace o1: {line(g_num)} (predicted = ours p1->p2 move "
        f"{line(rel(p2, p1))})\n"
        f"{tag} (C2''') separation |ours p2 - o1| / |ours p1 - o1| "
        + "  ".join(f"{c} {g_num[c] / t1[c]:.4g}x" if t1[c] > 0 else f"{c} inf"
                    for c in g_num)
        + " (predicted >= 25x)",
        flush=True,
    )
    if t2 is not None:
        g_rev = rel(p1, o2)
        print(
            f"{tag} (C2''') reverse |ours p1 - o2| / |ours p2 - o2| "
            + "  ".join(f"{c} {g_rev[c] / t2[c]:.4g}x" if t2[c] > 0 else f"{c} inf"
                        for c in g_rev),
            flush=True,
        )
        print(f"{tag} Palace o2 vs Palace o1 (h05x): {line(rel(o2, o1))}", flush=True)

    # --- gitignored source, printed only ------------------------------------
    p1_json = REPO / "logs" / "ans6_exact" / "ans6_exact_p1.json"
    if p1_json.exists():
        ours_s = np.array(
            [[complex(a, b) for a, b in row]
             for row in json.loads(p1_json.read_text())["s"]["64"]["pec"]]
        )
        print(
            f"{tag} gitignored source, printed only: ours p1 64 MHz pec 4x4 vs Palace "
            f"h05x_pec entry max|dS| {np.max(np.abs(ours_s - s_pec)):.3e}; sigma_max ours "
            f"{np.linalg.svd(ours_s, compute_uv=False).max():.12f} Palace "
            f"{id_pec['sigma_max']:.12f}",
            flush=True,
        )
    else:
        print(f"{tag} logs/ans6_exact/ans6_exact_p1.json: absent", flush=True)

    assert not fails, "; ".join(fails)
    print(f"{tag} (P) PASS on h05x_pec" + (" + h05x_pec_o2" if o2 is not None else "")
          + ("; (K') floor PASS" if csv_of("h05x_pmc").exists()
             else "; (K') not measured"), flush=True)


# `OPS-61` step 2b (§9 item 241, F23): Palace's copper column on the hole mesh.
CU_CONFIG_NAMES = ("h05x_cu", "h05x_cu_o2", "h05x_cu11")


def _cu_set():
    """name -> (sigma [S/m], Order); both sigmas imported from th14, never restated."""
    from tests.validation.test_th14_birdcage_copper import COPPER_SIGMA, SIGMA_LADDER

    return {
        "h05x_cu": (COPPER_SIGMA, 1),
        "h05x_cu_o2": (COPPER_SIGMA, 2),
        "h05x_cu11": (SIGMA_LADDER[-1], 1),
    }


def cu_configs():
    """Item 241 W1: write the three copper configs from ``h05x_pec.json``; assert (Z)."""
    from fem_em_solver.core.cavity import surface_resistance_ohm
    from fem_em_solver.utils.constants import MU_0
    from tests.validation.test_ops61_exact_strip_ports import STRIP_BAND
    from tests.validation.test_th14_birdcage_copper import (
        COPPER_SIGMA,
        FREQUENCIES_HZ,
        SIGMA_LADDER,
    )

    tag = "[OPS-61 2b W1]"
    try:
        from mpi4py import MPI

        assert MPI.COMM_WORLD.size == 1, f"run at -n 1, not {MPI.COMM_WORLD.size}"
    except ImportError:
        pass
    cu_set = _cu_set()
    assert tuple(cu_set) == CU_CONFIG_NAMES
    for name in cu_set:
        path = OUT_DIR / f"{name}.json"
        assert not path.exists(), f"{path.relative_to(REPO)} exists: never overwrite"
    wall = _hole_wall_attr()
    src = json.loads((OUT_DIR / "h05x_pec.json").read_text())
    assert src["Model"]["Mesh"] == HOLE_MESH_NAME, src["Model"]["Mesh"]
    ports = src["Boundaries"]["LumpedPort"]
    assert len(ports) == 4, f"{len(ports)} ports in h05x_pec.json"
    for p in ports:
        i = int(p["Index"])
        assert p["Attributes"] == [NARROW_BASE + i], (i, p["Attributes"])
        assert p["Direction"] == "+Z", (i, p["Direction"])
    assert src["Solver"]["Order"] == 1, src["Solver"]["Order"]
    samples = src["Solver"]["Driven"]["Samples"]
    assert len(samples) == 1 and samples[0]["Freq"] == [0.064], samples
    assert src["Boundaries"]["PEC"]["Attributes"] == [EXTERIOR_TAG, wall], src[
        "Boundaries"
    ]["PEC"]
    assert "Conductivity" not in src["Boundaries"], src["Boundaries"].keys()

    for name, (sigma, order) in cu_set.items():
        out = json.loads(json.dumps(src))
        out["Boundaries"]["PEC"]["Attributes"] = [EXTERIOR_TAG]
        out["Boundaries"]["Conductivity"] = [
            {"Attributes": [wall], "Conductivity": sigma, "Permeability": 1.0}
        ]
        out["Solver"]["Order"] = order
        out["Problem"]["Output"] = f"postpro_{name}"
        path = OUT_DIR / f"{name}.json"
        path.write_text(json.dumps(out, indent=2) + "\n")
        print(
            f"{tag} wrote {path.relative_to(REPO)}: Order {order}, sigma {sigma:.3e} S/m, "
            f"Output postpro_{name}\n{tag} {name} Boundaries (PEC + Conductivity): "
            + json.dumps({k: out["Boundaries"][k] for k in ("PEC", "Conductivity")}),
            flush=True,
        )

    # --- (Z) the surface-impedance tie, read back from each written JSON -------
    f64 = FREQUENCIES_HZ["64"]
    omega_ours = 2.0 * np.pi * f64
    mu0_codata = 1.25663706212e-6
    fails = []
    for name, (sigma_set, order) in cu_set.items():
        cfg = json.loads((OUT_DIR / f"{name}.json").read_text())
        blocks = cfg["Boundaries"]["Conductivity"]
        assert len(blocks) == 1, blocks
        b = blocks[0]
        assert b["Attributes"] == [wall], b
        assert "Thickness" not in b and "External" not in b, b
        assert cfg["Boundaries"]["PEC"]["Attributes"] == [EXTERIOR_TAG], cfg["Boundaries"]
        assert cfg["Solver"]["Order"] == order, cfg["Solver"]["Order"]
        sigma_json = float(b["Conductivity"])
        mur_json = float(b["Permeability"])
        freq_json = cfg["Solver"]["Driven"]["Samples"][0]["Freq"]
        assert len(freq_json) == 1, freq_json
        omega_json = 2.0 * np.pi * float(freq_json[0]) * 1.0e9
        expect_sigma = COPPER_SIGMA if name != "h05x_cu11" else SIGMA_LADDER[-1]
        if sigma_json != expect_sigma:
            fails.append(f"{name}: sigma_json {sigma_json!r} != {expect_sigma!r}")
        if mur_json != 1.0:
            fails.append(f"{name}: mu_r_json {mur_json!r} != 1.0")
        d_omega = abs(omega_json - omega_ours) / omega_ours
        if not d_omega <= STRIP_BAND:
            fails.append(f"{name}: |d omega|/omega {d_omega:.3e} > {STRIP_BAND}")
        z_ours = (1.0 + 1.0j) * surface_resistance_ohm(omega_ours, sigma_json)
        delta = np.sqrt(2.0 / (omega_json * MU_0 * mur_json * sigma_json))
        z_doc = (1.0 + 1.0j) / (sigma_json * delta)
        z_rel = abs(z_ours - z_doc) / abs(z_ours)
        if not z_rel <= STRIP_BAND:
            fails.append(f"{name}: (Z) |Z_ours - Z_doc|/|Z_ours| {z_rel:.3e} > {STRIP_BAND}")
        delta_codata = np.sqrt(2.0 / (omega_json * mu0_codata * mur_json * sigma_json))
        z_codata = (1.0 + 1.0j) / (sigma_json * delta_codata)
        print(
            f"{tag} (Z) {name}: sigma_json {sigma_json:.6e} S/m, mu_r {mur_json}, "
            f"Freq {freq_json[0]} GHz (|d omega|/omega {d_omega:.3e}); delta "
            f"{delta:.6e} m, R_s {z_doc.real:.6e} ohm, Z_doc {z_doc.real:.6e}(1+j) ohm, "
            f"|Z_s| {abs(z_doc):.6e} ohm; Z_ours {z_ours.real:.6e}(1+j) ohm "
            f"(ours logged 2.087161e-03(1+j) at 5.8e7); |Z_ours - Z_doc|/|Z_ours| "
            f"{z_rel:.3e} vs STRIP_BAND {STRIP_BAND}; CODATA-mu0 variant "
            f"|Z_codata - Z_doc|/|Z_doc| {abs(z_codata - z_doc) / abs(z_doc):.3e} (printed)",
            flush=True,
        )
    print(
        f"{tag} (Z) anchors on {list(cu_set)}: "
        f"{'PASS' if not fails else 'FAIL ' + '; '.join(fails)}",
        flush=True,
    )
    assert not fails, "; ".join(fails)


def read_hole_cu():
    """Item 241 W6: (P) on h05x_cu (+o2); (K1)/(K3) asserted; (T3)/(T3m)/(T4)/(K2) printed."""
    from tests.validation.test_ans6_exact_strip import (
        ANS6_EXACT_STRIP_CLASS_MEANS_P1,
        ANS6_EXACT_STRIP_CLASS_MEANS_P2,
    )
    from tests.validation.test_port21_feed_sensitivity import KNOB_REACHES_THE_SOLVE

    tag = "[OPS-61 2b W6]"

    def rel(a, b):
        return {k: abs(a[k] - b[k]) / abs(b[k]) for k in b}

    def line(d):
        return "  ".join(f"{k} {v * 100:.4f}%" for k, v in d.items())

    def sci(d):
        return "  ".join(f"{k} {v:.3e}" for k, v in d.items())

    def csv_of(name):
        return OUT_DIR / f"postpro_{name}" / "port-S.csv"

    def resolved_block(name):
        path = OUT_DIR / f"postpro_{name}" / f"{name}_resolved.json"
        if not path.exists():
            print(f"{tag} {path.relative_to(REPO)}: absent", flush=True)
            return
        b = json.loads(path.read_text()).get("Boundaries", {})
        print(
            f"{tag} {name} resolved Boundaries.Conductivity (file evidence): "
            f"{json.dumps(b.get('Conductivity'))}; PEC {json.dumps(b.get('PEC'))}",
            flush=True,
        )

    fails = []
    assert csv_of("h05x_cu").exists(), f"no {csv_of('h05x_cu').relative_to(REPO)}"
    assert csv_of("h05x_pec").exists(), f"no {csv_of('h05x_pec').relative_to(REPO)}"
    for name in CU_CONFIG_NAMES:
        resolved_block(name)
    s_cu, id_cu, f = _p_check("h05x_cu", csv_of("h05x_cu"), tag)
    fails += [f"h05x_cu: {x}" for x in f]
    id_pec = _identities(_read_port_s(csv_of("h05x_pec"))[0])
    o2 = None
    if csv_of("h05x_cu_o2").exists():
        _s2, o2, f2 = _p_check("h05x_cu_o2", csv_of("h05x_cu_o2"), tag)
        fails += [f"h05x_cu_o2: {x}" for x in f2]
    else:
        print(f"{tag} h05x_cu_o2: not measured", flush=True)

    p1_cu = dict(ANS6_EXACT_STRIP_CLASS_MEANS_P1["64"]["cu"])
    p1_pec = dict(ANS6_EXACT_STRIP_CLASS_MEANS_P1["64"]["pec"])
    p2_cu = dict(ANS6_EXACT_STRIP_CLASS_MEANS_P2["64"]["cu"])
    p2_pec = dict(ANS6_EXACT_STRIP_CLASS_MEANS_P2["64"]["pec"])
    pal_cu, pal_pec = id_cu["means"], id_pec["means"]

    def t3m(pal_c, pal_p, our_c, our_p, label):
        d_pal = {k: pal_c[k] - pal_p[k] for k in our_c}
        d_ours = {k: our_c[k] - our_p[k] for k in our_c}
        mag = {k: abs(d_pal[k] - d_ours[k]) / abs(d_ours[k]) for k in d_ours}
        ang = {k: float(np.degrees(np.angle(d_pal[k] / d_ours[k]))) for k in d_ours}
        print(
            f"{tag} {label} copper move, Palace vs ours: |d_pal - d_ours|/|d_ours| "
            f"{sci(mag)}; arg(d_pal/d_ours) "
            + "  ".join(f"{k} {v:+.4f} deg" for k, v in ang.items())
            + f"\n{tag} {label} moves |d|/|pec mean|: Palace {sci(rel(pal_c, pal_p))}; "
            f"ours {sci(rel(our_c, our_p))} (printed; predicted <= 1e-3 and |arg| < 0.1 deg)",
            flush=True,
        )
        return mag, ang

    # --- (T3) / (T3m) ours vs Palace at order 1 (printed) --------------------
    t3 = rel(p1_cu, pal_cu)
    print(
        f"\n{tag} (T3) ours ANS6_EXACT_STRIP_CLASS_MEANS_P1['64']['cu'] vs Palace h05x_cu "
        f"(o1): |d mean|/|mean| {line(t3)} = {sci(t3)} (printed; predicted <= 1e-6)",
        flush=True,
    )
    t3m(pal_cu, pal_pec, p1_cu, p1_pec, "(T3m) o1")

    # --- (T4) order 2 (printed) -----------------------------------------------
    if o2 is not None and csv_of("h05x_pec_o2").exists():
        pal_pec_o2 = _identities(_read_port_s(csv_of("h05x_pec_o2"))[0])["means"]
        t4 = rel(p2_cu, o2["means"])
        print(
            f"{tag} (T4) ours ANS6_EXACT_STRIP_CLASS_MEANS_P2['64']['cu'] vs Palace "
            f"h05x_cu_o2: |d mean|/|mean| {line(t4)} = {sci(t4)} (printed; predicted <= 1e-6)",
            flush=True,
        )
        t3m(o2["means"], pal_pec_o2, p2_cu, p2_pec, "(T4m) o2")
    else:
        print(f"{tag} (T4) not measured", flush=True)

    # --- (C2'''') order-mismatch separation (printed) -------------------------
    g_num = rel(p2_cu, pal_cu)
    print(
        f"{tag} (C2'''') ours p2 cu vs Palace h05x_cu: {line(g_num)} (predicted = ours "
        f"p1->p2 cu move {line(rel(p2_cu, p1_cu))})\n"
        f"{tag} (C2'''') separation |ours p2 - o1| / |ours p1 - o1| "
        + "  ".join(f"{c} {g_num[c] / t3[c]:.4g}x" if t3[c] > 0 else f"{c} inf"
                    for c in g_num)
        + " (predicted >= 25x)",
        flush=True,
    )

    # --- (K1) the sigma boundary reaches Palace's operator: floor asserted ----
    k1 = rel(pal_cu, pal_pec)
    k1_ok = min(k1.values()) > KNOB_REACHES_THE_SOLVE
    print(
        f"{tag} (K1) h05x_cu vs h05x_pec class move {sci(k1)}; floor > "
        f"{KNOB_REACHES_THE_SOLVE:.0e} on every class: {'PASS' if k1_ok else 'FAIL'} "
        f"(ours 6.420e-4 / 2.811e-4 / 4.861e-4)",
        flush=True,
    )
    if not k1_ok:
        fails.append(f"(K1) min class move {min(k1.values()):.3e}")

    # --- (K3) loss direction: asserted -----------------------------------------
    k3_ok = id_cu["sigma_max"] < id_pec["sigma_max"]
    print(
        f"{tag} (K3) sigma_max h05x_cu {id_cu['sigma_max']:.12f} < h05x_pec "
        f"{id_pec['sigma_max']:.12f} (gap {id_pec['sigma_max'] - id_cu['sigma_max']:.3e}): "
        f"{'PASS' if k3_ok else 'FAIL'}; predicted Palace cu = ours 0.999851545 to ~1e-9",
        flush=True,
    )
    if not k3_ok:
        fails.append(f"(K3) sigma_max cu {id_cu['sigma_max']!r} >= pec {id_pec['sigma_max']!r}")

    # --- (K2) sigma scaling: predicted, printed -------------------------------
    if csv_of("h05x_cu11").exists():
        s11, id11, fk = _p_check("h05x_cu11", csv_of("h05x_cu11"), tag + " [control]")
        print(f"{tag} h05x_cu11 (P) printed only (control, never asserted): "
              f"{'clean' if not fk else '; '.join(fk)}", flush=True)
        k2 = rel(id11["means"], pal_pec)
        s_pec = _read_port_s(csv_of("h05x_pec"))[0]
        ent_cu = np.max(np.abs(s_cu - s_pec))
        ent_11 = np.max(np.abs(s11 - s_pec))
        print(
            f"{tag} (K2) h05x_cu11 vs h05x_pec class move {sci(k2)}; ratio (K1)/(K2) "
            + "  ".join(f"{c} {k1[c] / k2[c]:.4g}" if k2[c] > 0 else f"{c} inf" for c in k2)
            + f"; entry-level max|dS| cu {ent_cu:.3e} / cu11 {ent_11:.3e} = "
            f"{ent_cu / ent_11 if ent_11 > 0 else float('inf'):.4g} (predicted ~100; "
            f"ours 99.83; sigma_max cu11 {id11['sigma_max']:.12f}, ours 0.999900690)",
            flush=True,
        )
    else:
        print(f"{tag} (K2) h05x_cu11: not measured", flush=True)

    # --- gitignored source, printed only ------------------------------------
    p1_json = REPO / "logs" / "ans6_exact" / "ans6_exact_p1.json"
    if p1_json.exists():
        ours_s = np.array(
            [[complex(a, b) for a, b in row]
             for row in json.loads(p1_json.read_text())["s"]["64"]["cu"]]
        )
        print(
            f"{tag} gitignored source, printed only: ours p1 64 MHz cu 4x4 vs Palace "
            f"h05x_cu entry max|dS| {np.max(np.abs(ours_s - s_cu)):.3e}; sigma_max ours "
            f"{np.linalg.svd(ours_s, compute_uv=False).max():.12f} Palace "
            f"{id_cu['sigma_max']:.12f}",
            flush=True,
        )
    else:
        print(f"{tag} logs/ans6_exact/ans6_exact_p1.json: absent", flush=True)

    assert not fails, "; ".join(fails)
    print(f"{tag} (P) PASS on h05x_cu" + (" + h05x_cu_o2" if o2 is not None else "")
          + "; (K1) floor PASS; (K3) PASS", flush=True)


# ---------------------------------------------------------------------------
# `OPS-61` step 1f (§9 item 418, F34): the h = 0.0075 exact-strip rung.
def export_exact_strip_h0075():
    """W1 of item 418: the h = 0.0075 exact-strip mesh + four single-excitation configs.

    The sheet steps, the N_i strip selection, the triangle / tet blocks and the
    Gmsh 2.2 writer are COPIED from ``export(exact_strip=True)`` (commit
    4989ebdd, lines 143-383 of this file); ``export()`` is never called or
    edited. Never overwrites. Anchors (E1)-(E4) of item 418.
    """
    import dolfinx
    from mpi4py import MPI

    from tests.mesh.helpers import global_cell_tag_set
    from tests.mesh.test_birdcage_leg_offset import (
        _build,
        _projected_extents,
        _sheet_areas,
        _sheet_azimuth_deg,
    )
    from tests.mesh.test_birdcage_port_sheets import SHEET_IFACE
    from tests.mesh.test_birdcage_port_tags import LEG_COUNT
    from tests.validation.test_ans4_exact_strip_records import (
        SPEC_AREA_M2,
        SPEC_HEIGHT_M,
        SPEC_WIDTH_M,
    )
    from tests.validation.test_ops61_exact_strip_ports import STRIP_BAND
    from tests.validation.test_port_birdcage_leg_offset_sweep import (
        _narrowed_radial,
        _port_frame,
    )
    from tests.validation.test_port_lumped_narrowed_sheet import GATED_WIDTH_FRACTION

    tag = "[OPS-61 1f W1]"
    comm = MPI.COMM_WORLD
    assert comm.size == 1, "export is a serial writer: run it at mpiexec -n 1"
    t_start = time.perf_counter()
    ports_idx = list(range(1, LEG_COUNT + 1))

    # 1. never overwrite
    mesh_path = OUT_DIR / H0075_MESH_NAME
    cfg_paths = {name: OUT_DIR / f"{name}.json" for name in H0075_SET}
    for p in [mesh_path, *cfg_paths.values()]:
        assert not p.exists(), f"{p.relative_to(REPO)} exists: never overwrite"

    # 2. the template, checked the way o2_configs checks its own
    src = json.loads((OUT_DIR / "f05x_o2_f10.json").read_text())
    assert src["Model"]["Mesh"] == STRIP_MESH_NAME, src["Model"]["Mesh"]
    ports = src["Boundaries"]["LumpedPort"]
    assert len(ports) == 4, f"{len(ports)} ports in f05x_o2_f10.json"
    for p in ports:
        i = int(p["Index"])
        assert p["Attributes"] == [NARROW_BASE + i], (i, p["Attributes"])
        assert p["Direction"] == "+Z", (i, p["Direction"])
        assert float(p["R"]) == 50.0, (i, p["R"])
    assert src["Solver"]["Order"] == 2, src["Solver"]["Order"]
    samples = src["Solver"]["Driven"]["Samples"]
    assert len(samples) == 1 and samples[0]["Freq"] == [0.010], samples

    # 3. the build entry 27 measured (exact_rungs under RESOLUTION=0.0075)
    msh, cell_tags, _facet_tags, diag, t_mesh = _build(
        np.zeros(LEG_COUNT),
        resolution=H0075_RES,
        leg_sheet_strip_fraction=GATED_WIDTH_FRACTION,
    )
    tdim = msh.topology.dim
    fdim = tdim - 1
    ncells = int(msh.topology.index_map(tdim).size_global)
    print(
        f"{tag} mesh built in {t_mesh:.2f} s: {ncells} tets (record "
        f"{H0075_CELL_RECORD}); resolution {H0075_RES}; diag port_box_size_m = "
        f"{diag['port_box_size_m']}",
        flush=True,
    )
    # (E1) before anything is written
    assert ncells == H0075_CELL_RECORD, (
        f"(E1) W1 count {ncells} != {H0075_CELL_RECORD}: stop before any Palace window"
    )
    msh.topology.create_connectivity(fdim, tdim)
    msh.topology.create_entity_permutations()

    # 4. --- copied from export(exact_strip=True): sheet steps -----------------
    tags_f, full_areas = _sheet_areas(msh, cell_tags, ports_idx, comm)
    dx, _dy, dz = diag["port_box_size_m"]
    sheet_analytic = float(dx) * float(dz)
    geometry = {}
    for i in ports_idx:
        stag = SHEET_IFACE + i
        azimuth = _sheet_azimuth_deg(msh, tags_f, stag, comm)
        radial, azimuthal, axial = _port_frame(azimuth)
        spans, centres = _projected_extents(
            msh, tags_f, stag, comm, [radial, azimuthal, axial]
        )
        geometry[i] = {
            "tag": stag,
            "azimuth_deg": float(azimuth),
            "radial": radial,
            "frame": [radial, azimuthal, axial],
            "w_full": float(spans[0]),
            "h_full": float(spans[2]),
            "centre_radial": float(centres[0]),
            "area_full": float(full_areas[i]),
            "area_ratio_full": float(full_areas[i] / sheet_analytic),
        }
    tags_full = tags_f
    for i in ports_idx:
        g = geometry[i]
        tags_f = _narrowed_radial(
            msh, tags_f, g["tag"], GATED_WIDTH_FRACTION, g["centre_radial"],
            g["radial"], 0.5 * g["w_full"],
        )

    x = np.asarray(msh.geometry.x, dtype=float)

    def facet_nodes(facets):
        facets = np.asarray(facets, dtype=np.int32)
        if not facets.size:
            return np.zeros((0, 3), dtype=np.int64)
        return np.asarray(
            dolfinx.cpp.mesh.entities_to_geometry(msh._cpp_object, fdim, facets, False)
        ).reshape(-1, 3)

    tri_blocks = []
    ext = dolfinx.mesh.exterior_facet_indices(msh.topology)
    ext_nodes = facet_nodes(ext)
    ext_area = float(np.sum(_tri_areas(x, ext_nodes)))
    tri_blocks.append((ext_nodes, EXTERIOR_TAG))
    strips = {}
    for i in ports_idx:
        stag = SHEET_IFACE + i
        full = np.asarray(tags_full.find(stag), dtype=np.int32)
        narrow = np.asarray(tags_f.find(stag), dtype=np.int32)
        subset = bool(np.all(np.isin(narrow, full)))
        rest = np.setdiff1d(full, narrow)
        n_nodes, r_nodes = facet_nodes(narrow), facet_nodes(rest)
        a_n = float(np.sum(_tri_areas(x, n_nodes)))
        pts = x[np.unique(n_nodes.ravel())]
        proj = pts @ np.asarray(geometry[i]["frame"], dtype=float).T
        span = proj.max(axis=0) - proj.min(axis=0)
        h = float(span[2])
        strips[i] = {"area": a_n, "h": h, "w": a_n / h, "w_bbox": float(span[0]),
                     "subset": subset, "n": int(narrow.size), "rest": int(rest.size)}
        print(
            f"    P{i}: azimuth {geometry[i]['azimuth_deg']:8.4f} deg; F_{i} "
            f"{full.size} tris (F/dxdz {geometry[i]['area_ratio_full']:.12f}); N_{i} "
            f"{narrow.size} tris area {a_n:.12e} (rel {a_n / SPEC_AREA_M2 - 1:+.3e}); "
            f"h {h:.12e} (rel {h / SPEC_HEIGHT_M - 1:+.3e}); A/h {a_n / h:.12e} "
            f"(rel {a_n / h / SPEC_WIDTH_M - 1:+.3e}); radial span {span[0]:.12e}; "
            f"F\\N {rest.size} tris; N subset F: {subset}",
            flush=True,
        )
        tri_blocks.append((n_nodes, NARROW_BASE + i))
        tri_blocks.append((r_nodes, REMAINDER_BASE + i))
    ext_err = abs(ext_area - EXTERIOR_AREA_M2)
    print(
        f"{tag} exterior: {ext.size} tris, area {ext_area:.12e} m^2 vs "
        f"{EXTERIOR_AREA_M2} (|diff| {ext_err:.3e}, band {EXTERIOR_AREA_BAND:.0e})",
        flush=True,
    )

    # --- copied: Gmsh 2.2 ASCII writer -----------------------------------------
    OUT_DIR.mkdir(parents=True, exist_ok=True)
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
    total_vol = float(np.sum(np.abs(vol6)) / 6.0)
    ntri = sum(b.shape[0] for b, _ in tri_blocks)
    with open(mesh_path, "w") as fh:
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
    lo, hi = x.min(axis=0), x.max(axis=0)
    tags = sorted(global_cell_tag_set(msh, cell_tags))
    counts = {t: int(np.count_nonzero(ctag == t)) for t in tags}
    print(
        f"{tag} wrote {mesh_path.relative_to(REPO)}: {x.shape[0]} nodes, {ncl} tets "
        f"({int(neg.sum())} reoriented), {ntri} tris; volume {total_vol:.9e} m^3; "
        f"bbox lo {lo} hi {hi}",
        flush=True,
    )
    print(f"{tag} cell tags (count): {counts}", flush=True)
    template_air = sorted(src["Domains"]["Materials"][2]["Attributes"])
    print(f"{tag} template air attributes {template_air}; mesh tags {tags}", flush=True)

    # 5. four single-excitation configs from the template
    for name, k in H0075_SET.items():
        out = json.loads(json.dumps(src))
        out["Model"]["Mesh"] = H0075_MESH_NAME
        for p in out["Boundaries"]["LumpedPort"]:
            p["Excitation"] = int(p["Index"]) == k
        out["Problem"]["Output"] = f"postpro_{name}"
        cfg_paths[name].write_text(json.dumps(out, indent=2) + "\n")
        print(f"{tag} wrote {cfg_paths[name].relative_to(REPO)}: excitation port {k}, "
              f"Output postpro_{name}", flush=True)
    print(json.dumps(json.loads(cfg_paths["f05x_h0075_o2_f10_e1"].read_text())
                     ["Boundaries"], indent=2), flush=True)
    print(f"{tag} export total {time.perf_counter() - t_start:.2f} s (mesh {t_mesh:.2f} s)",
          flush=True)

    # 6. (E2)-(E4) after writing, so a failure leaves the artefacts to inspect
    fails = []
    if not ext_err <= EXTERIOR_AREA_BAND:
        fails.append(f"(E2) exterior area {ext_area:.12e}")
    for i, s in strips.items():
        if not abs(s["area"] / SPEC_AREA_M2 - 1.0) <= STRIP_BAND:
            fails.append(f"(E3) N_{i} area {s['area']:.12e}")
        if not abs(s["w"] / SPEC_WIDTH_M - 1.0) <= STRIP_BAND:
            fails.append(f"(E3) N_{i} A/h {s['w']:.12e}")
        if not s["subset"]:
            fails.append(f"(E4) N_{i} not a subset of F_{i}")
    print(
        f"{tag} ANCHORS (E1) tets == {H0075_CELL_RECORD}: PASS; (E2) exterior area "
        f"within {EXTERIOR_AREA_BAND:.0e}; (E3) N_i area / A/h within STRIP_BAND "
        f"{STRIP_BAND:.0e}; (E4) N_i subset F_i: "
        f"{'PASS' if not fails else 'FAIL ' + '; '.join(fails)}",
        flush=True,
    )
    assert not fails, "; ".join(fails)


def _p_check_s(name, s, tag):
    """``_p_check``'s assertions on an assembled matrix (``_p_check`` not edited)."""
    from tests.validation.test_port_birdcage_four_port import (
        PASSIVITY_SIGMA_TOLERANCE,
        POOLED_SEPARATION_FLOOR,
        _circulant_classes,
        _class_spread,
    )
    from tests.validation.test_port_birdcage_lumped_column import ADJACENT_SPREAD_BAND
    from tests.validation.test_port_lumped_sheet_sweep import RECIPROCITY_BAND

    ident = _identities(s)
    classes = _circulant_classes(s)
    column_power = np.sum(np.abs(s) ** 2, axis=0)
    worst = max(ident["spreads"].values())
    pooled = _class_spread(np.concatenate([classes["adjacent"], classes["opposite"]]))
    f = []
    if not ident["reciprocity"] <= RECIPROCITY_BAND:
        f.append(f"reciprocity {ident['reciprocity']:.3e} > {RECIPROCITY_BAND}")
    if not ident["sigma_max"] <= 1.0 + PASSIVITY_SIGMA_TOLERANCE:
        f.append(f"sigma_max {ident['sigma_max']:.12f} > 1 + {PASSIVITY_SIGMA_TOLERANCE}")
    for k, v in enumerate(column_power, start=1):
        if not v <= 1.0 + 1.0e-9:
            f.append(f"column {k} power {v:.12f} > 1 + 1e-9")
    for k, v in ident["spreads"].items():
        if not v <= ADJACENT_SPREAD_BAND:
            f.append(f"{k} spread {v:.3e} > {ADJACENT_SPREAD_BAND}")
    if not pooled >= POOLED_SEPARATION_FLOOR * worst:
        f.append(f"pooled {pooled:.3e} < {POOLED_SEPARATION_FLOOR} x worst {worst:.3e}")
    spreads = "  ".join(f"{k} {v * 100:.4f}%" for k, v in ident["spreads"].items())
    print(
        f"\n{tag} {name}: reciprocity {ident['reciprocity']:.3e}  sigma_max "
        f"{ident['sigma_max']:.12f}  spreads {spreads}\n    column power sums "
        + ", ".join(f"{v:.12f}" for v in column_power)
        + f"; pooled off-diagonal {pooled * 100:.4f}% vs floor {POOLED_SEPARATION_FLOOR}"
        f" x worst {worst * 100:.4f}% (separation {pooled / worst:.2f}x)",
        flush=True,
    )
    print(
        f"{tag} (P) {name} anchors: reciprocity <= {RECIPROCITY_BAND:.0e}, sigma_max <= "
        f"1 + {PASSIVITY_SIGMA_TOLERANCE:.0e}, column power <= 1 + 1e-9, spreads <= "
        f"{ADJACENT_SPREAD_BAND}, pooled >= {POOLED_SEPARATION_FLOOR} x worst: "
        f"{'PASS' if not f else 'FAIL ' + '; '.join(f)}",
        flush=True,
    )
    return ident, f


def _ours_h0075_s():
    """Ours at h = 0.0075, degree 2, kappa off, 10 MHz, parsed from the tracked log."""
    import gzip

    path = H0075_OURS_LOG
    if path.exists():
        lines = path.read_text().splitlines()
    else:
        with gzip.open(str(path) + ".gz", "rt") as fh:
            lines = fh.read().splitlines()
    idx = [n for n, ln in enumerate(lines) if ln.strip() == H0075_OURS_HEADER]
    assert len(idx) == 1, f"{len(idx)} headers {H0075_OURS_HEADER!r}"
    tok = re.compile(r"[+-]\d\.\d+e[+-]\d+[+-]\d\.\d+e[+-]\d+j")
    rows = []
    for k, ln in enumerate(lines[idx[0] + 1: idx[0] + 5], start=1):
        body = re.sub(r"^(PASSED)+", "", ln.strip()).strip()
        assert body.startswith(f"S_{k}k ="), (k, ln)
        vals = tok.findall(body)
        assert len(vals) == 4, (k, ln)
        rows.append([complex(v) for v in vals])
    return np.array(rows), idx[0] + 1


def read_exact_strip_h0075():
    """W6 of item 418: assemble Palace's 4x4, (P) asserted, (Tie) asserted;
    (T5) / (C-h) / (C2-h) and the price rule printed, never gated."""
    from tests.validation.test_ans4_exact_strip_records import ANS4_EXACT_STRIP_S_P2_X1
    from tests.validation.test_port_birdcage_four_port import _circulant_classes

    tag = "[OPS-61 1f W6]"

    def means(m):
        return {c: complex(np.mean(v)) for c, v in _circulant_classes(m).items()}

    def rel(a, b):
        return {c: abs(a[c] - b[c]) / abs(b[c]) for c in b}

    def line(d):
        return "  ".join(f"{c} {v * 100:.4f}% ({v:.3e})" for c, v in d.items())

    # --- ours, parsed; (Tie) ---------------------------------------------------
    s_ours, hdr_line = _ours_h0075_s()
    sig_ours = float(np.linalg.svd(s_ours, compute_uv=False).max())
    tie = abs(sig_ours - H0075_OURS_SIGMA_MAX)
    print(
        f"{tag} ours parsed from {H0075_OURS_LOG.relative_to(REPO)}:{hdr_line} "
        f"(+4 lines): sigma_max {sig_ours:.12f} vs printed {H0075_OURS_SIGMA_MAX} "
        f"|diff| {tie:.3e} (Tie band {H0075_TIE_BAND:.0e}): "
        f"{'PASS' if tie <= H0075_TIE_BAND else 'FAIL'}",
        flush=True,
    )
    for r in s_ours:
        print("    ours " + "  ".join(_fmt(z) for z in r), flush=True)
    assert tie <= H0075_TIE_BAND, f"(Tie) |diff| {tie:.3e}: wrong block parsed"

    # --- Palace columns ---------------------------------------------------------
    s = np.full((4, 4), np.nan + 0j)
    present = []
    for name, k in H0075_SET.items():
        path = OUT_DIR / f"postpro_{name}" / "port-S.csv"
        if not path.exists():
            print(f"{tag} {name}: no {path.relative_to(REPO)} (not measured)", flush=True)
            continue
        sk, header = _read_port_s(path)
        cols = [c for c in range(sk.shape[1]) if not np.all(np.isnan(sk[:, c]))]
        labels = [h for h in header if re.search(r"S\[\d+\]\[\d+\]", h)]
        print(f"{tag} {name}: S columns present {[c + 1 for c in cols]}; labels "
              f"{labels[:2]} ...", flush=True)
        if (k - 1) in cols and sk.shape[0] == 4:
            col = sk[:, k - 1]
            print(f"{tag} {name}: kept column labelled S[i][{k}] (asserted)", flush=True)
        else:
            assert len(cols) == 1 and sk.shape[0] == 4, (name, cols, sk.shape)
            col = sk[:, cols[0]]
            print(f"{tag} {name}: READING -- Palace labelled the single excitation "
                  f"S[i][{cols[0] + 1}], kept by position as column {k}", flush=True)
        assert not np.any(np.isnan(col)), (name, col)
        s[:, k - 1] = col
        present.append(k)
    print(f"{tag} Palace 4x4 (columns present {present}):", flush=True)
    for r in s:
        print("    Palace " + "  ".join(_fmt(z) for z in r), flush=True)

    # --- (T5) per entry of the present columns ----------------------------------
    for k in present:
        g = np.abs(s[:, k - 1] - s_ours[:, k - 1]) / np.abs(s_ours[:, k - 1])
        print(f"{tag} (T5) column {k} per entry |Palace - ours|/|ours|: "
              + "  ".join(f"{v * 100:.4f}% ({v:.3e})" for v in g), flush=True)

    fails = []
    if len(present) == 4:
        ident, fails = _p_check_s("f05x_h0075_o2_f10 (assembled)", s, tag)
        pal = means(s)
        t5 = rel(pal, means(s_ours))
        print(f"{tag} (T5) per class |Palace o2 h - ours p2 h|/|ours p2 h|: {line(t5)} "
              "(predicted <= 1e-6; printed, never gated)", flush=True)
        x1_path = OUT_DIR / "postpro_f05x_o2_f10" / "port-S.csv"
        pal_x1 = means(_read_port_s(x1_path)[0])
        ch = rel(pal, pal_x1)
        print(f"{tag} (C-h) Palace o2 h=0.0075 vs Palace o2 x1 per class: {line(ch)}; "
              "ours' h-move " + "  ".join(f"{c} {v * 100:.4f}%"
                                         for c, v in H0075_OURS_H_MOVE.items())
              + " (log :3726; printed)", flush=True)
        ours_x1 = means(ANS4_EXACT_STRIP_S_P2_X1["10"])
        num = {c: abs(ours_x1[c] - pal[c]) for c in pal}
        den = {c: abs(means(s_ours)[c] - pal[c]) for c in pal}
        print(
            f"{tag} (C2-h) |ours x1 p2 - Palace o2 h| / |ours h p2 - Palace o2 h| per "
            "class: " + "  ".join(
                f"{c} {num[c]:.3e}/{den[c]:.3e} = "
                f"{(num[c] / den[c]) if den[c] > 0 else float('inf'):.4g}x" for c in pal)
            + " (predicted >= 447x on adjacent, vs item 229's 25x bar; printed)",
            flush=True,
        )
    else:
        print(f"{tag} (P), per-class (T5), (C-h), (C2-h): not measured "
              f"(columns present {present})", flush=True)

    # --- price readout from the W2-W5 harness logs --------------------------------
    logs = sorted((REPO / "docs" / "testing" / "logs").glob(
        "*_OPS-61-item418-W[2-5]-palace-h0075-e*.log"))
    r_f = {"64": 37 / 26, "128": 55 / 26}
    proj128 = []
    for lg in logs:
        txt = lg.read_text()
        # Elapsed Time Report columns Min / Max / Avg: keep Avg (brief).
        tot = re.findall(r"^\s*Total\s+[\d.]+\s+[\d.]+\s+([\d.]+)\s*$", txt, re.M)
        setup = re.findall(r"^\s*Setup\s+[\d.]+\s+[\d.]+\s+([\d.]+)\s*$", txt, re.M)
        its = re.findall(r"GMRES solver converged in (\d+) iteration", txt)
        if not tot or not setup:
            print(f"{tag} price {lg.name}: no Total/Setup line", flush=True)
            continue
        t_e, s_e = float(tot[-1]), float(setup[-1])
        p = {f: s_e + (t_e - s_e) * r for f, r in r_f.items()}
        p4 = {f: s_e + 4 * (t_e - s_e) * r for f, r in r_f.items()}
        proj128.append(p["128"])
        print(
            f"{tag} price {lg.name}: Total {t_e:.1f} s, Setup {s_e:.1f} s, GMRES its "
            f"{its}; projected single-excitation 64 MHz {p['64']:.1f} s, 128 MHz "
            f"{p['128']:.1f} s; projected four-excitation one-command 64 MHz "
            f"{p4['64']:.1f} s, 128 MHz {p4['128']:.1f} s",
            flush=True,
        )
    if proj128:
        heavy = max(proj128) <= 450.0
        print(f"{tag} price rule: max_k projected T_e(128) {max(proj128):.1f} s vs 450 s "
              f"-> 64 / 128 MHz {'HEAVY (four single-excitation windows per frequency)' if heavy else 'xl cost probe (10-14 weekly)'}",
              flush=True)

    assert not fails, "; ".join(fails)
    print(f"{tag} done: (Tie) PASS" + ("; (P) PASS" if len(present) == 4 else ""),
          flush=True)


def h0075_freq_configs(f):
    """W1 of items 434 / 435: (E1') on the reused mesh; four single-excitation
    configs at frequency ``f`` from item 418's 10 MHz configs (Freq / Output only)."""
    import hashlib

    try:
        from mpi4py import MPI

        assert MPI.COMM_WORLD.size == 1, "h0075-configs runs at -n 1"
    except ImportError:
        pass
    tag = f"[OPS-61 h0075 {f} MHz W1]"
    targets = {name: OUT_DIR / f"{name}.json" for name in H0075_FREQ_SET[f]}
    for p in targets.values():
        assert not p.exists(), f"{p} exists: never overwrite"

    # (E1') tet count parsed from the Gmsh 2.2 $Elements block
    mesh_path = OUT_DIR / H0075_MESH_NAME
    raw = mesh_path.read_bytes()
    lines = raw.decode().splitlines()
    i0 = lines.index("$Elements")
    i1 = lines.index("$EndElements")
    nel = int(lines[i0 + 1])
    ntet = sum(1 for ln in lines[i0 + 2: i1] if ln.split()[1] == "4")
    assert i1 - i0 - 2 == nel, (nel, i1 - i0 - 2)
    print(f"{tag} {mesh_path.relative_to(REPO)}: sha256 "
          f"{hashlib.sha256(raw).hexdigest()}, size {len(raw)} B; {nel} elements, "
          f"{ntet} type-4 tets", flush=True)
    print(f"{tag} (E1') tets {ntet} == {H0075_CELL_RECORD}: "
          f"{'PASS' if ntet == H0075_CELL_RECORD else 'FAIL'}", flush=True)
    assert ntet == H0075_CELL_RECORD, f"(E1') {ntet} != {H0075_CELL_RECORD}: stop"

    for name, k in H0075_FREQ_SET[f].items():
        tpl_path = OUT_DIR / f"f05x_h0075_o2_f10_e{k}.json"
        tpl = json.loads(tpl_path.read_text())
        assert tpl["Model"]["Mesh"] == H0075_MESH_NAME, tpl["Model"]["Mesh"]
        ports = tpl["Boundaries"]["LumpedPort"]
        assert [int(p["Index"]) for p in ports] == [1, 2, 3, 4], ports
        for p in ports:
            assert p["Attributes"] == [20 + int(p["Index"])], p
            assert p["Direction"] == "+Z" and float(p["R"]) == 50.0, p
            assert bool(p["Excitation"]) == (int(p["Index"]) == k), p
        assert tpl["Solver"]["Order"] == 2, tpl["Solver"]["Order"]
        samples = tpl["Solver"]["Driven"]["Samples"]
        assert len(samples) == 1 and samples[0]["Freq"] == [0.010], samples
        out = json.loads(json.dumps(tpl))
        out["Solver"]["Driven"]["Samples"][0]["Freq"] = [H0075_FREQS[f]]
        out["Problem"]["Output"] = f"postpro_{name}"
        targets[name].write_text(json.dumps(out, indent=2) + "\n")

        # template diff: exactly the two keys
        back = json.loads(targets[name].read_text())
        t2 = json.loads(json.dumps(tpl))
        b2 = json.loads(json.dumps(back))
        assert b2["Solver"]["Driven"]["Samples"][0]["Freq"] == [H0075_FREQS[f]]
        assert b2["Problem"]["Output"] == f"postpro_{name}"
        for d in (t2, b2):
            d["Solver"]["Driven"]["Samples"][0]["Freq"] = None
            d["Problem"]["Output"] = None
        assert t2 == b2, f"{name}: differs from its template beyond Freq / Output"
        print(f"{tag} wrote {targets[name].relative_to(REPO)} from {tpl_path.name}: "
              f"excitation port {k}, Freq [{H0075_FREQS[f]}] GHz, Output "
              f"postpro_{name}; template diff = Freq, Output only: PASS", flush=True)
    print(f"{tag} done: (E1') PASS; four configs written", flush=True)


def _ours_rung_s(log, header):
    """``_ours_h0075_s`` with the log and header as arguments (that one not edited)."""
    import gzip

    path = log
    if path.exists():
        lines = path.read_text().splitlines()
    else:
        with gzip.open(str(path) + ".gz", "rt") as fh:
            lines = fh.read().splitlines()
    idx = [n for n, ln in enumerate(lines) if ln.strip() == header]
    assert len(idx) == 1, f"{len(idx)} headers {header!r}"
    tok = re.compile(r"[+-]\d\.\d+e[+-]\d+[+-]\d\.\d+e[+-]\d+j")
    rows = []
    for k, ln in enumerate(lines[idx[0] + 1: idx[0] + 5], start=1):
        body = re.sub(r"^(PASSED)+", "", ln.strip()).strip()
        assert body.startswith(f"S_{k}k ="), (k, ln)
        vals = tok.findall(body)
        assert len(vals) == 4, (k, ln)
        rows.append([complex(v) for v in vals])
    return np.array(rows), idx[0] + 1


def read_exact_strip_h0075_freq(f):
    """W6 of items 434 / 435: a copy of ``read_exact_strip_h0075`` (item 418's,
    not edited) with the frequency's set, ours' block, sigma_max, h-move and the
    x1 Palace output substituted; (P) and (Tie) asserted, the rest printed."""
    from tests.validation.test_ans4_exact_strip_records import ANS4_EXACT_STRIP_S_P2_X1
    from tests.validation.test_port_birdcage_four_port import _circulant_classes

    blk = H0075_OURS_BLOCKS[f]
    item = H0075_ITEM[f]
    tag = f"[OPS-61 h0075 {f} MHz W6]"

    def means(m):
        return {c: complex(np.mean(v)) for c, v in _circulant_classes(m).items()}

    def rel(a, b):
        return {c: abs(a[c] - b[c]) / abs(b[c]) for c in b}

    def line(d):
        return "  ".join(f"{c} {v * 100:.4f}% ({v:.3e})" for c, v in d.items())

    s_ours, hdr_line = _ours_rung_s(H0075_OURS_LOG, blk["header"])
    sig_ours = float(np.linalg.svd(s_ours, compute_uv=False).max())
    tie = abs(sig_ours - blk["sigma_max"])
    print(
        f"{tag} ours parsed from {H0075_OURS_LOG.relative_to(REPO)}:{hdr_line} "
        f"(+4 lines): sigma_max {sig_ours:.12f} vs printed {blk['sigma_max']} "
        f"|diff| {tie:.3e} (Tie band {H0075_TIE_BAND:.0e}): "
        f"{'PASS' if tie <= H0075_TIE_BAND else 'FAIL'}",
        flush=True,
    )
    for r in s_ours:
        print("    ours " + "  ".join(_fmt(z) for z in r), flush=True)
    assert tie <= H0075_TIE_BAND, f"(Tie) |diff| {tie:.3e}: wrong block parsed"

    s = np.full((4, 4), np.nan + 0j)
    present = []
    for name, k in H0075_FREQ_SET[f].items():
        path = OUT_DIR / f"postpro_{name}" / "port-S.csv"
        if not path.exists():
            print(f"{tag} {name}: no {path.relative_to(REPO)} (not measured)", flush=True)
            continue
        sk, header = _read_port_s(path)
        cols = [c for c in range(sk.shape[1]) if not np.all(np.isnan(sk[:, c]))]
        labels = [h for h in header if re.search(r"S\[\d+\]\[\d+\]", h)]
        print(f"{tag} {name}: S columns present {[c + 1 for c in cols]}; labels "
              f"{labels[:2]} ...", flush=True)
        if (k - 1) in cols and sk.shape[0] == 4:
            col = sk[:, k - 1]
            print(f"{tag} {name}: kept column labelled S[i][{k}] (asserted)", flush=True)
        else:
            assert len(cols) == 1 and sk.shape[0] == 4, (name, cols, sk.shape)
            col = sk[:, cols[0]]
            print(f"{tag} {name}: READING -- Palace labelled the single excitation "
                  f"S[i][{cols[0] + 1}], kept by position as column {k}", flush=True)
        assert not np.any(np.isnan(col)), (name, col)
        s[:, k - 1] = col
        present.append(k)
    print(f"{tag} Palace 4x4 (columns present {present}):", flush=True)
    for r in s:
        print("    Palace " + "  ".join(_fmt(z) for z in r), flush=True)

    for k in present:
        g = np.abs(s[:, k - 1] - s_ours[:, k - 1]) / np.abs(s_ours[:, k - 1])
        print(f"{tag} (T5) column {k} per entry |Palace - ours|/|ours|: "
              + "  ".join(f"{v * 100:.4f}% ({v:.3e})" for v in g), flush=True)

    fails = []
    if len(present) == 4:
        ident, fails = _p_check_s(f"f05x_h0075_o2_f{f} (assembled)", s, tag)
        pal = means(s)
        t5 = rel(pal, means(s_ours))
        print(f"{tag} (T5) per class |Palace o2 h - ours p2 h|/|ours p2 h|: {line(t5)} "
              "(predicted <= 1e-6; printed, never gated)", flush=True)
        x1_path = OUT_DIR / blk["pal_x1"] / "port-S.csv"
        pal_x1 = means(_read_port_s(x1_path)[0])
        ch = rel(pal, pal_x1)
        print(f"{tag} (C-h) Palace o2 h=0.0075 vs Palace o2 x1 ({blk['pal_x1']}) per "
              f"class: {line(ch)}; ours' h-move "
              + "  ".join(f"{c} {v * 100:.4f}%" for c, v in blk["h_move"].items())
              + f" (log {blk['h_move_line']}; printed)", flush=True)
        ours_x1 = means(ANS4_EXACT_STRIP_S_P2_X1[f])
        num = {c: abs(ours_x1[c] - pal[c]) for c in pal}
        den = {c: abs(means(s_ours)[c] - pal[c]) for c in pal}
        print(
            f"{tag} (C2-h) |ours x1 p2 - Palace o2 h| / |ours h p2 - Palace o2 h| per "
            "class: " + "  ".join(
                f"{c} {num[c]:.3e}/{den[c]:.3e} = "
                f"{(num[c] / den[c]) if den[c] > 0 else float('inf'):.4g}x" for c in pal)
            + " (vs item 229's 25x bar; printed)",
            flush=True,
        )
    else:
        print(f"{tag} (P), per-class (T5), (C-h), (C2-h): not measured "
              f"(columns present {present})", flush=True)

    logs = sorted((REPO / "docs" / "testing" / "logs").glob(
        f"*_OPS-61-item{item}-W[2-5]-palace-h0075-f{f}-e*.log"))
    for lg in logs:
        txt = lg.read_text()
        tot = re.findall(r"^\s*Total\s+[\d.]+\s+[\d.]+\s+([\d.]+)\s*$", txt, re.M)
        setup = re.findall(r"^\s*Setup\s+[\d.]+\s+[\d.]+\s+([\d.]+)\s*$", txt, re.M)
        its = re.findall(r"GMRES solver converged in (\d+) iteration", txt)
        if not tot or not setup:
            print(f"{tag} price {lg.name}: no Total/Setup line", flush=True)
            continue
        print(f"{tag} price {lg.name}: measured Total {float(tot[-1]):.1f} s, Setup "
              f"{float(setup[-1]):.1f} s, GMRES its {its}; item 418's projection "
              f"{blk['proj'][0]:.1f}-{blk['proj'][1]:.1f} s", flush=True)

    assert not fails, "; ".join(fails)
    print(f"{tag} done: (Tie) PASS" + ("; (P) PASS" if len(present) == 4 else ""),
          flush=True)


# ---------------------------------------------------------------------------
# `OPS-61` step 3a (§9 item 420, §10 T39): Palace twin of the ESR-tuned F-small
# at degree 2, 64.000 MHz.  Port 1 (feed) is the bare 50 ohm port; ports 2-4 are
# unexcited LumpedPorts carrying the gap capacitor.  Three configs:
#   t39_esr      -- ports 2-4 R = R_p, C = C_p (ESR ruling (iii)'s exact
#                   single-frequency series -> parallel map of r + C_in(r));
#   t39_lossless -- ports 2-4 C = C_in(0) only (convention-free C path);
#   t39_short    -- ports 2-4 R = 0.2, C = C_in(r): ruling (iii)'s misread,
#                   the labelled negative control / convention discriminator.
# Only S[1][1] is compared (S[i][1], i >= 2, are referenced to Palace's port
# impedances, not 50 ohm).  Never reads z_matrix_ohm, the aed key, docs/private/
# or aed_results/.
T39_NAMES = ("t39_esr", "t39_lossless", "t39_short")
T39_ESR_OHM = 0.2
# Ours' T26 Route A S11_A, 64.0 MHz, 12-digit prints (brief §9 item 420):
T39_S11_A_ESR = -2.392625126044e-01 + 9.637828488809e-01j  # ...item395-W0-esr64.log:2128
T39_S11_A_LOSSLESS = -2.404062978868e-01 + 9.685876774136e-01j  # ...item311-W1-cin64-lossless.log:2055
T39_C_IN_0_LITERAL = 6.359009101426113e-11  # ...item311-W1-cin64-lossless.log:2058
T39_TIE_PREDICTED = 1.0e-6  # predicted, printed, never asserted (no cross-solver band, Q1)


def _t39_map(omega, c_f, r):
    """ESR ruling (iii): series r + C -> parallel (R_p, C_p), exact at omega."""
    q = 1.0 / (omega * c_f * r)
    return q, r * (1.0 + q * q), c_f * q * q / (1.0 + q * q)


def ans7_configs():
    """Item 420 W1: write t39_esr / t39_lossless / t39_short from h05x_cu_o2.json; assert (Z)."""
    from fem_em_solver.ports.lumped import series_rlc_impedance
    from tests.validation.test_ans7_exact_strip import T26_C_IN_R_X1, TUNE_TARGET_HZ
    from tests.validation.test_ops61_exact_strip_ports import STRIP_BAND
    from tests.validation.test_th14_birdcage_copper import COPPER_SIGMA

    tag = "[OPS-61 3a W1]"
    try:
        from mpi4py import MPI

        assert MPI.COMM_WORLD.size == 1, f"run at -n 1, not {MPI.COMM_WORLD.size}"
    except ImportError:
        pass
    for name in T39_NAMES:
        path = OUT_DIR / f"{name}.json"
        assert not path.exists(), f"{path.relative_to(REPO)} exists: never overwrite"
    wall = _hole_wall_attr()
    src = json.loads((OUT_DIR / "h05x_cu_o2.json").read_text())
    assert src["Model"]["Mesh"] == HOLE_MESH_NAME, src["Model"]["Mesh"]
    ports = src["Boundaries"]["LumpedPort"]
    assert len(ports) == 4, f"{len(ports)} ports in h05x_cu_o2.json"
    for p in ports:
        i = int(p["Index"])
        assert p["Attributes"] == [NARROW_BASE + i], (i, p["Attributes"])
        assert p["Direction"] == "+Z", (i, p["Direction"])
        assert p["R"] == 50.0, (i, p["R"])
    assert src["Solver"]["Order"] == 2, src["Solver"]["Order"]
    samples = src["Solver"]["Driven"]["Samples"]
    assert len(samples) == 1 and samples[0]["Freq"] == [0.064], samples
    assert src["Boundaries"]["PEC"]["Attributes"] == [EXTERIOR_TAG], src["Boundaries"]["PEC"]
    cond = src["Boundaries"]["Conductivity"]
    assert len(cond) == 1 and cond[0]["Attributes"] == [wall], cond
    assert cond[0]["Conductivity"] == COPPER_SIGMA, cond

    ans7 = REPO / "logs" / "ans7_exact"
    c_in_r_json = json.loads((ans7 / "tune_p2_esr0.2.json").read_text())["c_in_r"]
    c_in_0 = json.loads((ans7 / "tune_p2.json").read_text())["c_in"]
    assert c_in_r_json == T26_C_IN_R_X1, (c_in_r_json, T26_C_IN_R_X1)
    assert c_in_0 == T39_C_IN_0_LITERAL, (c_in_0, T39_C_IN_0_LITERAL)
    c_in_r = float(T26_C_IN_R_X1)

    omega_set = 2.0 * np.pi * float(samples[0]["Freq"][0]) * 1.0e9
    q, r_p, c_p = _t39_map(omega_set, c_in_r, T39_ESR_OHM)
    gap = {
        "t39_esr": {"R": r_p, "C": c_p},
        "t39_lossless": {"C": c_in_0},
        "t39_short": {"R": T39_ESR_OHM, "C": c_in_r},  # the labelled misread control
    }
    for name in T39_NAMES:
        out = json.loads(json.dumps(src))
        for p in out["Boundaries"]["LumpedPort"]:
            if int(p["Index"]) == 1:
                continue  # feed: R 50, Excitation as the template (excited)
            p.pop("R")
            p["Excitation"] = False
            p.update(gap[name])
        out["Problem"]["Output"] = f"postpro_{name}"
        path = OUT_DIR / f"{name}.json"
        path.write_text(json.dumps(out, indent=2) + "\n")
        print(f"{tag} wrote {path.relative_to(REPO)}; Boundaries.LumpedPort: "
              + json.dumps(out["Boundaries"]["LumpedPort"]), flush=True)

    # --- (Z) the map's own identity, read back from each written JSON ----------
    fails = []
    f_t = float(TUNE_TARGET_HZ)
    z_ser = series_rlc_impedance(f_t, r_ohm=T39_ESR_OHM, c_f=c_in_r)
    for name in T39_NAMES:
        cfg = json.loads((OUT_DIR / f"{name}.json").read_text())
        freq_json = cfg["Solver"]["Driven"]["Samples"][0]["Freq"]
        assert len(freq_json) == 1, freq_json
        omega_json = 2.0 * np.pi * float(freq_json[0]) * 1.0e9
        d_omega = abs(omega_json - 2.0 * np.pi * f_t) / (2.0 * np.pi * f_t)
        if not d_omega <= STRIP_BAND:
            fails.append(f"{name}: |d omega|/omega {d_omega:.3e} > {STRIP_BAND}")
        lp = cfg["Boundaries"]["LumpedPort"]
        p1 = [p for p in lp if int(p["Index"]) == 1][0]
        if p1["R"] != 50.0 or p1["Excitation"] in (False, 0):
            fails.append(f"{name}: port 1 {p1}")
        for p in lp:
            if int(p["Index"]) == 1:
                continue
            if p["Excitation"] is not False:
                fails.append(f"{name}: port {p['Index']} Excitation {p['Excitation']!r}")
            if {k: p.get(k) for k in ("R", "C")} != {k: gap[name].get(k) for k in ("R", "C")}:
                fails.append(f"{name}: port {p['Index']} R/C {p.get('R')!r}/{p.get('C')!r}")
        if name == "t39_esr":
            p2 = [p for p in lp if int(p["Index"]) == 2][0]
            z_par = 1.0 / (1.0 / float(p2["R"]) + 1j * omega_json * float(p2["C"]))
            z_rel = abs(z_par - z_ser) / abs(z_ser)
            if not z_rel <= STRIP_BAND:
                fails.append(f"(Z) {name}: |Z_par - Z_ser|/|Z_ser| {z_rel:.3e} > {STRIP_BAND}")
            print(
                f"{tag} (Z) {name}: Freq {freq_json[0]} GHz (|d omega|/omega {d_omega:.3e}); "
                f"Q {q:.8f}, R_p {r_p:.8f} ohm, C_p {c_p:.10e} F (C_p/C - 1 {c_p / c_in_r - 1:.4e}); "
                f"Z_par {_fmt(z_par)} ohm, Z_ser {_fmt(z_ser)} ohm; |Z_par - Z_ser|/|Z_ser| "
                f"{z_rel:.3e} vs STRIP_BAND {STRIP_BAND} (predicted ~4e-16)",
                flush=True,
            )
            w, l = 7.0e-3, 8.0e-3
            print(f"{tag} predicted sheet values (w 7 mm, l 8 mm): R_p w/l {r_p * w / l:.3e} ohm/sq, "
                  f"C_p l/w {c_p * l / w:.3e} F/sq", flush=True)
        elif name == "t39_short":
            print(f"{tag} {name} (control, printed): |Z| as parallel 0.2 ohm || C = "
                  f"{abs(1.0 / (1.0 / T39_ESR_OHM + 1j * omega_json * c_in_r)):.4f} ohm "
                  f"vs series |Z_ser| {abs(z_ser):.3f} ohm", flush=True)
    print(f"{tag} (Z) anchors on {list(T39_NAMES)}: "
          f"{'PASS' if not fails else 'FAIL ' + '; '.join(fails)}", flush=True)
    assert not fails, "; ".join(fails)


def read_ans7_esr():
    """Item 420 W6: S[1][1] per t39 run against ours' Route A; (T6)/(T6-0)/(K-*) printed."""
    from tests.validation.test_ans7_exact_strip import _ANS7

    tag = "[OPS-61 3a W6]"

    def rel(a, b):
        return abs(a - b) / abs(b)

    s11 = {}
    for name in T39_NAMES:
        post = OUT_DIR / f"postpro_{name}"
        csv_path = post / "port-S.csv"
        res = post / f"{name}_resolved.json"
        if res.exists():
            b = json.loads(res.read_text()).get("Boundaries", {})
            print(f"{tag} {name} resolved Boundaries.LumpedPort (file evidence): "
                  f"{json.dumps(b.get('LumpedPort'))}", flush=True)
        else:
            print(f"{tag} {res.relative_to(REPO)}: absent", flush=True)
        if not csv_path.exists():
            print(f"{tag} {name}: not measured ({csv_path.relative_to(REPO)} absent)", flush=True)
            continue
        s, _hdr = _read_port_s(csv_path)
        s11[name] = complex(s[0, 0])
        total = None
        pj = post / "palace.json"
        if pj.exists():
            total = json.loads(pj.read_text())["ElapsedTime"]["Durations"]["Total"]
        print(f"{tag} {name}: S11_Pal {_fmt(s11[name])}, |S11| {abs(s11[name]):.12f}, "
              f"Z_bare {_fmt(complex(_ANS7._z_from_s11(s11[name])))} ohm; Palace Total "
              f"{total if total is None else f'{total:.3f}'} s", flush=True)

    print(f"{tag} ours S11_A (r, C_in(r)) {_fmt(T39_S11_A_ESR)}, |S11_A| "
          f"{abs(T39_S11_A_ESR):.12f}, Z_bare {_fmt(complex(_ANS7._z_from_s11(T39_S11_A_ESR)))} ohm; "
          f"lossless {_fmt(T39_S11_A_LOSSLESS)}, |S11| {abs(T39_S11_A_LOSSLESS):.12f}, "
          f"Z_bare {_fmt(complex(_ANS7._z_from_s11(T39_S11_A_LOSSLESS)))} ohm", flush=True)
    t6 = rel(s11["t39_esr"], T39_S11_A_ESR) if "t39_esr" in s11 else None
    t60 = rel(s11["t39_lossless"], T39_S11_A_LOSSLESS) if "t39_lossless" in s11 else None
    ks = rel(s11["t39_short"], T39_S11_A_ESR) if "t39_short" in s11 else None

    def show(x):
        return "not measured" if x is None else f"{x:.3e}"

    print(f"{tag} (T6) |S11_Pal(t39_esr) - S11_A|/|S11_A| {show(t6)} (printed; predicted "
          f"<= {T39_TIE_PREDICTED:.0e} if Palace's lumped R || C is parallel)", flush=True)
    print(f"{tag} (T6-0) |S11_Pal(t39_lossless) - S11_A,lossless|/|S11_A,lossless| {show(t60)} "
          f"(printed; predicted <= {T39_TIE_PREDICTED:.0e}, convention-free)", flush=True)
    ceiling = (1.0 + abs(T39_S11_A_ESR)) / abs(T39_S11_A_ESR)
    print(f"{tag} (K-short) |S11_Pal(t39_short) - S11_A|/|S11_A| {show(ks)} (control; "
          f"predicted in [5e-2, 2.0] if parallel; ceiling {ceiling:.3f})", flush=True)
    if "t39_esr" in s11 and "t39_lossless" in s11:
        k_pal = rel(s11["t39_esr"], s11["t39_lossless"])
        k_ours = rel(T39_S11_A_ESR, T39_S11_A_LOSSLESS)
        print(f"{tag} (K-esr) Palace esr-vs-lossless move {k_pal:.6e}; ours {k_ours:.6e}; "
              f"|Pal/ours - 1| {abs(k_pal / k_ours - 1.0):.3e} (predicted <= 1e-6; a dropped R "
              f"reads ~1e-10)", flush=True)
    if t6 is not None and ks is not None:
        par = t6 <= T39_TIE_PREDICTED
        ser = ks <= T39_TIE_PREDICTED
        reading = ("parallel" if par and not ser else "series" if ser and not par
                   else "neither" if not par and not ser else "both (degenerate)")
        print(f"{tag} convention reading: {reading} ((T6) {t6:.3e}, (K-short) {ks:.3e} vs "
              f"{T39_TIE_PREDICTED:.0e})", flush=True)
    print(f"{tag} done (all cross-solver figures printed, none asserted)", flush=True)


# --- OPS-61 step 3b (§9 item 421, §10 T40): Palace's driven sweep on T39 ------
# T26's Route A triple (63.9 / 64.0 / 64.1 MHz, TUNE_MHZ); 64.000 is item 420's
# t39_esr run (never re-run).  Item 420's convention reading is "parallel"
# (20261013T140547Z_OPS-61-item420-W6-read.log:49), so every frequency gets its
# own exact (R_p(f), C_p(f)) map at C_in(r) -- never t39_esr's values reused
# (ESR ruling (iii): R_p ~ 1/f^2 off the conversion frequency).  The feed's
# (r + C) branch is post-processing (ruling 311-Q5), as T26's Route A.
T40_SWEEP = (("t40_esr_f63.9", 0.0639, 63.9e6), ("t39_esr", 0.064, 64.0e6),
             ("t40_esr_f64.1", 0.0641, 64.1e6))
# Ours' S11_A at (r, C_in(r)), 12-digit prints, ...item395-W0-esr64.log:2127-2129:
T40_S11_A = (-2.437299889957e-01 + 9.627593952511e-01j,
             T39_S11_A_ESR,
             -2.347525671823e-01 + 9.647916492670e-01j)
# T26 (T2-r) rows A, ...item395-W0-esr64.log:2130:
T40_F0_A = 63999834.680890724
T40_R_IN_A = 3175.9794661876795
T40_Q_A = 164.5434331272106
T40_IM_Y_A = (-1.616350144e-04, +2.676560526e-07, +1.624710261e-04)
# T40's Q target, quoted from PROJECT_PLAN §10 T40 ("Q within the 5 % the copper
# contribution is quoted at"): the copper contribution's AGREE level as the
# 2026-09-30 weekly quoted it (§6 Phase 6 row, "AGREE on the mode and the copper
# contribution (5 %)").  Not a code band chosen here; review ruling 2026-10-13
# 03:00 (3) accepted it provisionally for the 10-14 weekly (Q2).
T40_Q_TARGET = 0.05


def ans7_sweep_configs():
    """Item 421 W1: write t40_esr_f63.9 / t40_esr_f64.1 from t39_esr.json; assert (Z')."""
    from fem_em_solver.ports.lumped import series_rlc_impedance
    from tests.validation.test_ans7_exact_strip import T26_C_IN_R_X1
    from tests.validation.test_ops61_exact_strip_ports import STRIP_BAND

    tag = "[OPS-61 3b W1]"
    try:
        from mpi4py import MPI

        assert MPI.COMM_WORLD.size == 1, f"run at -n 1, not {MPI.COMM_WORLD.size}"
    except ImportError:
        pass
    new = [(n, fg, fh) for n, fg, fh in T40_SWEEP if n.startswith("t40_")]
    for name, _fg, _fh in new:
        path = OUT_DIR / f"{name}.json"
        assert not path.exists(), f"{path.relative_to(REPO)} exists: never overwrite"
    src = json.loads((OUT_DIR / "t39_esr.json").read_text())
    samples = src["Solver"]["Driven"]["Samples"]
    assert len(samples) == 1 and samples[0]["Freq"] == [0.064], samples
    c_in_r = float(T26_C_IN_R_X1)
    for name, f_ghz, _fh in new:
        out = json.loads(json.dumps(src))
        out["Solver"]["Driven"]["Samples"][0]["Freq"] = [f_ghz]
        omega = 2.0 * np.pi * f_ghz * 1.0e9
        _q, r_p, c_p = _t39_map(omega, c_in_r, T39_ESR_OHM)
        for p in out["Boundaries"]["LumpedPort"]:
            if int(p["Index"]) == 1:
                continue
            p["R"] = r_p
            p["C"] = c_p
        out["Problem"]["Output"] = f"postpro_{name}"
        path = OUT_DIR / f"{name}.json"
        path.write_text(json.dumps(out, indent=2) + "\n")
        print(f"{tag} wrote {path.relative_to(REPO)}; Freq {[f_ghz]} GHz; "
              f"Boundaries.LumpedPort: " + json.dumps(out["Boundaries"]["LumpedPort"]), flush=True)

    # --- (Z') per config, read back; and the template diff ---------------------
    fails = []
    for name, _fg, f_hz in new:
        cfg = json.loads((OUT_DIR / f"{name}.json").read_text())
        freq_json = cfg["Solver"]["Driven"]["Samples"][0]["Freq"]
        assert len(freq_json) == 1, freq_json
        omega_json = 2.0 * np.pi * float(freq_json[0]) * 1.0e9
        d_omega = abs(omega_json - 2.0 * np.pi * f_hz) / (2.0 * np.pi * f_hz)
        if not d_omega <= STRIP_BAND:
            fails.append(f"{name}: |d omega|/omega {d_omega:.3e} > {STRIP_BAND}")
        z_ser = series_rlc_impedance(f_hz, r_ohm=T39_ESR_OHM, c_f=c_in_r)
        for p in cfg["Boundaries"]["LumpedPort"]:
            if int(p["Index"]) == 1:
                if p["R"] != 50.0 or p["Excitation"] in (False, 0) or "C" in p:
                    fails.append(f"{name}: port 1 {p}")
                continue
            if p["Excitation"] is not False:
                fails.append(f"{name}: port {p['Index']} Excitation {p['Excitation']!r}")
            z_par = 1.0 / (1.0 / float(p["R"]) + 1j * omega_json * float(p["C"]))
            z_rel = abs(z_par - z_ser) / abs(z_ser)
            if not z_rel <= STRIP_BAND:
                fails.append(f"(Z') {name} port {p['Index']}: {z_rel:.3e} > {STRIP_BAND}")
            if int(p["Index"]) == 2:
                q, r_p, c_p = _t39_map(omega_json, c_in_r, T39_ESR_OHM)
                print(f"{tag} (Z') {name}: Freq {freq_json[0]} GHz (|d omega|/omega "
                      f"{d_omega:.3e}); Q {q:.8f}, R_p {float(p['R']):.8f} ohm, C_p "
                      f"{float(p['C']):.10e} F; Z_par {_fmt(z_par)} ohm, Z_ser {_fmt(z_ser)} ohm; "
                      f"|Z_par - Z_ser|/|Z_ser| {z_rel:.3e} vs STRIP_BAND {STRIP_BAND} "
                      f"(predicted ~1e-16)", flush=True)
        # only Freq, ports 2-4 R/C and Problem.Output differ from t39_esr.json
        a = json.loads(json.dumps(src))
        b = json.loads(json.dumps(cfg))
        for d in (a, b):
            d["Solver"]["Driven"]["Samples"][0].pop("Freq")
            d["Problem"].pop("Output")
            for p in d["Boundaries"]["LumpedPort"]:
                if int(p["Index"]) != 1:
                    p.pop("R")
                    p.pop("C")
        if a != b:
            fails.append(f"{name}: differs from t39_esr.json beyond Freq / ports 2-4 R, C / Output")
    # (N-map), printed: t39_esr's 64.000 MHz (R_p, C_p) misused at 63.9 / 64.1 MHz
    _q64, r_p64, c_p64 = _t39_map(2.0 * np.pi * 64.0e6, c_in_r, T39_ESR_OHM)
    for _n, _fg, f_hz in new:
        z_ser = series_rlc_impedance(f_hz, r_ohm=T39_ESR_OHM, c_f=c_in_r)
        z_mis = 1.0 / (1.0 / r_p64 + 1j * 2.0 * np.pi * f_hz * c_p64)
        print(f"{tag} (N-map, printed) 64-MHz map at {f_hz / 1e6:.1f} MHz: |Z_par(64-map) - "
              f"Z_ser|/|Z_ser| {abs(z_mis - z_ser) / abs(z_ser):.3e} (predicted 1.60e-5; "
              f"STRIP_BAND {STRIP_BAND})", flush=True)
    print(f"{tag} (Z') anchors on {[n for n, _f, _h in new]}: "
          f"{'PASS' if not fails else 'FAIL ' + '; '.join(fails)}", flush=True)
    assert not fails, "; ".join(fails)


def read_ans7_sweep():
    """Item 421 W4: Palace's f0 / R_in / Q via parallel_rows against T26's Route A rows."""
    from tests.validation.test_ans7_exact_strip import (
        T21P_ZERO_RESOLUTION,
        T26_C_IN_R_X1,
        _ANS7,
        _z_branch,
        parallel_rows,
    )

    tag = "[OPS-61 3b W4]"
    c_in_r = float(T26_C_IN_R_X1)
    freqs = np.array([fh for _n, _fg, fh in T40_SWEEP])
    s11, totals = [], []
    for name, _fg, f_hz in T40_SWEEP:
        post = OUT_DIR / f"postpro_{name}"
        csv_path = post / "port-S.csv"
        if not csv_path.exists():
            print(f"{tag} {name}: not measured ({csv_path.relative_to(REPO)} absent)", flush=True)
            s11.append(None)
            totals.append(None)
            continue
        s, _hdr = _read_port_s(csv_path)
        s11.append(complex(s[0, 0]))
        pj = post / "palace.json"
        totals.append(json.loads(pj.read_text())["ElapsedTime"]["Durations"]["Total"]
                      if pj.exists() else None)
    fails = []
    if any(x is None for x in s11):
        print(f"{tag} a point is not measured: (C1) cannot be read", flush=True)
        assert False, "point not measured"

    def rel(a, b):
        return abs(a - b) / abs(b)

    z_bare, z_in, z_in_nofeed_r = [], [], []
    for k, (name, _fg, f_hz) in enumerate(T40_SWEEP):
        omega = 2.0 * np.pi * f_hz
        zb = complex(_ANS7._z_from_s11(s11[k]))
        zbr = _z_branch(omega, c_in_r, T39_ESR_OHM)
        z_bare.append(zb)
        z_in.append(1.0 / (1.0 / zb + 1.0 / zbr))
        z_in_nofeed_r.append(1.0 / (1.0 / zb + 1j * omega * c_in_r))
        tie = rel(s11[k], T40_S11_A[k])
        tot = totals[k]
        print(f"{tag} {f_hz / 1e6:.1f} MHz ({name}): S11_Pal {_fmt(s11[k])}, ours S11_A "
              f"{_fmt(T40_S11_A[k])}; |S11_Pal - S11_A|/|S11_A| {tie:.3e} (printed; predicted <= "
              f"{T39_TIE_PREDICTED:.0e}); Z_bare {_fmt(zb)} ohm; Palace Total "
              f"{tot if tot is None else f'{tot:.3f}'} s", flush=True)
    rows = parallel_rows(freqs, np.array(z_in))
    b = (1.0 / np.array(z_in)).imag
    print(f"{tag} Palace Im Y_in at (63.9, 64.0, 64.1) MHz = {b[0]:+.9e}, {b[1]:+.9e}, "
          f"{b[2]:+.9e} S; ours {T40_IM_Y_A[0]:+.9e}, {T40_IM_Y_A[1]:+.9e}, "
          f"{T40_IM_Y_A[2]:+.9e} S", flush=True)
    for k in range(3):
        print(f"{tag}   Im Y_in rel gap at {freqs[k] / 1e6:.1f} MHz: "
              f"{abs(b[k] - T40_IM_Y_A[k]) / abs(T40_IM_Y_A[k]):.3e} (printed)", flush=True)
    print(f"{tag} rows (Palace): {json.dumps(rows)}", flush=True)
    # ours' rows recomputed from the 12-digit S11_A prints (truncation check, printed)
    z_in_ours = [1.0 / (1.0 / complex(_ANS7._z_from_s11(T40_S11_A[k]))
                        + 1.0 / _z_branch(2.0 * np.pi * freqs[k], c_in_r, T39_ESR_OHM))
                 for k in range(3)]
    rows_o = parallel_rows(freqs, np.array(z_in_ours))
    print(f"{tag} ours' rows from the 12-digit S11_A prints: f0 {rows_o.get('f0_hz')!r} Hz, "
          f"R_in {rows_o.get('r_in_ohm')!r} ohm, Q {rows_o.get('q_loaded')!r} (vs logged "
          f"{T40_F0_A!r} / {T40_R_IN_A!r} / {T40_Q_A!r})", flush=True)

    # --- (C1) one crossing ------------------------------------------------------
    up, down = rows["crossings_up_hz"], rows["crossings_down_hz"]
    c1 = len(up) == 1 and not down and 63.9e6 < up[0] < 64.1e6
    print(f"{tag} (C1) up-crossings {up}, down-crossings {down}: "
          f"{'PASS' if c1 else 'FAIL'} (ASSERTED: exactly one - -> +, none + -> -, inside "
          f"(63.9, 64.1) MHz)", flush=True)
    if not c1:
        fails.append("(C1)")
    else:
        f0, r_in, q = rows["f0_hz"], rows["r_in_ohm"], rows["q_loaded"]
        df0 = f0 / T40_F0_A - 1.0
        dq = q / T40_Q_A - 1.0
        drin = r_in / T40_R_IN_A - 1.0
        f0_ok = abs(df0) <= T21P_ZERO_RESOLUTION
        q_ok = abs(dq) <= T40_Q_TARGET
        print(f"{tag} (f0) f0,Pal {f0!r} Hz vs f0,A {T40_F0_A!r} Hz: f0,Pal/f0,A - 1 = "
              f"{df0:+.3e} (ASSERTED |.| <= T21P_ZERO_RESOLUTION {T21P_ZERO_RESOLUTION:.3e}: "
              f"{'PASS' if f0_ok else 'FAIL'}; predicted <= 1e-6); f0,Pal/64 MHz - 1 = "
              f"{f0 / 64.0e6 - 1.0:+.6e}", flush=True)
        print(f"{tag} (Q) Q_Pal {q!r} vs Q_A {T40_Q_A!r}: Q_Pal/Q_A - 1 = {dq:+.3e} (ASSERTED "
              f"|.| <= T40_Q_TARGET {T40_Q_TARGET} (§10 T40's 5 %): {'PASS' if q_ok else 'FAIL'}; "
              f"predicted <= 1e-6)", flush=True)
        print(f"{tag} R_in,Pal {r_in!r} ohm vs R_in,A {T40_R_IN_A!r} ohm: R_in,Pal/R_in,A - 1 = "
              f"{drin:+.3e} (printed; predicted <= 1e-6)", flush=True)
        if not f0_ok:
            fails.append("(f0)")
        if not q_ok:
            fails.append("(Q)")

    # --- negative controls (rule (e); predicted, printed) -----------------------
    rows_nf = parallel_rows(freqs, np.array(z_in_nofeed_r))
    if "q_loaded" in rows_nf:
        qn = rows_nf["q_loaded"]
        in_band = 0.20 <= qn / T40_Q_A - 1.0 <= 0.80
        print(f"{tag} (N-feed) feed's r omitted: Q' {qn!r}, Q'/Q_A - 1 = {qn / T40_Q_A - 1.0:+.3e} "
              f"(predicted in [+20 %, +80 %]: {'met' if in_band else 'NOT met'}; ceiling Q_0 "
              f"971.07 -> +490 %); f0' {rows_nf['f0_hz']!r} Hz (printed)", flush=True)
    else:
        print(f"{tag} (N-feed) feed's r omitted: no - -> + crossing ({json.dumps(rows_nf)}) "
              f"(predicted a crossing with Q' in [+20 %, +80 %]: NOT met)", flush=True)
    rows_nb = parallel_rows(freqs, np.array(z_bare))
    bb = (1.0 / np.array(z_bare)).imag
    if "f0_hz" in rows_nb:
        print(f"{tag} (N-nofeed) feed branch dropped: crossing at {rows_nb['f0_hz']!r} Hz "
              f"(predicted no crossing: NOT met); Im Y_bare {bb[0]:+.4e}, {bb[1]:+.4e}, "
              f"{bb[2]:+.4e} S", flush=True)
    else:
        print(f"{tag} (N-nofeed) feed branch dropped: no crossing / f0 undefined (predicted: met); "
              f"Im Y_bare {bb[0]:+.4e}, {bb[1]:+.4e}, {bb[2]:+.4e} S; rows {json.dumps(rows_nb)}",
              flush=True)
    print(f"{tag} anchors (C1)/(f0)/(Q): {'PASS' if not fails else 'FAIL ' + ', '.join(fails)}; "
          f"ties and R_in printed, none asserted (ruling 2026-10-13 03:00 (1))", flush=True)
    assert not fails, ", ".join(fails)


# --- OPS-61 step 3c (§9 item 436, §10 T42, provisional): items 420 + 421 at ----
# the second Larmor frequency.  Palace's driven read of the ESR-tuned F-small at
# 127.74 MHz against T35's Route A (item 387; reproduced by item 419 W0, the
# window that wrote the files read here).  Each frequency gets its own exact
# parallel map of (r, C_in(r)) with C_in(r) = T35_C_IN_R_X1 (never a 64 MHz
# literal); the feed's (r + C) branch is post-processing (ruling 311-Q5).
T42_SWEEP = (("t42_esr_f127.64", 0.12764, 127.64e6), ("t42_esr_f127.74", 0.12774, 127.74e6),
             ("t42_esr_f127.84", 0.12784, 127.84e6))
T42_LOSSLESS = ("t42_lossless_f127.74", 0.12774, 127.74e6)
# Ours' S11_A at (r, C_in(r)), 12-digit prints, 20261013T124951Z_ANS-7-item419-W0-esr128.log:2127-2129:
T42_S11_A = (4.493753890384e-01 + 8.875155598136e-01j,
             4.513684481906e-01 + 8.864778722343e-01j,
             4.533642518309e-01 + 8.854328851717e-01j)
# T35 (T2-r) rows A, ...item419-W0-esr128.log:2130:
T42_F0_A = 127739705.96749091
T42_R_IN_A = 9790.70360584452
T42_Q_A = 254.46925909529418
T42_IM_Y_A = (-4.057389094e-05, +1.196522461e-07, +4.085113315e-05)
# Lossless C path at 127.74 MHz, C_in(0): 20261006T171745Z_ANS-7-item283-W2-cin128.log:2055
T42_S11_A_LOSSLESS = 4.523660701161e-01 + 8.884107356215e-01j
# C_in(0) at 127.74 MHz: ...item283-W2-cin128.log:2058; ...item419-W0-esr128.log:393
T42_C_IN_0_LITERAL = 1.5274843607471887e-11


def _t42_names():
    return [n for n, _g, _h in T42_SWEEP] + [T42_LOSSLESS[0]]


def ans7_128_configs():
    """Item 436 W1: write t42_esr_f127.{64,74,84} / t42_lossless_f127.74; assert (Z') / (Z0)."""
    from fem_em_solver.ports.lumped import series_rlc_impedance
    from tests.validation.test_ans7_exact_strip import T35_C_IN_R_X1
    from tests.validation.test_ops61_exact_strip_ports import STRIP_BAND

    tag = "[OPS-61 3c W1]"
    try:
        from mpi4py import MPI

        assert MPI.COMM_WORLD.size == 1, f"run at -n 1, not {MPI.COMM_WORLD.size}"
    except ImportError:
        pass
    for name in _t42_names():
        path = OUT_DIR / f"{name}.json"
        assert not path.exists(), f"{path.relative_to(REPO)} exists: never overwrite"
    t128 = REPO / "logs" / "ans7_exact" / "t127.74"
    c_in_r_json = json.loads((t128 / "tune_p2_esr0.2.json").read_text())["c_in_r"]
    c_in_0 = json.loads((t128 / "tune_p2.json").read_text())["c_in"]
    assert c_in_r_json == T35_C_IN_R_X1, (c_in_r_json, T35_C_IN_R_X1)
    assert c_in_0 == T42_C_IN_0_LITERAL, (c_in_0, T42_C_IN_0_LITERAL)
    c_in_r = float(T35_C_IN_R_X1)
    print(f"{tag} C_in(r) {c_in_r!r} F (T35_C_IN_R_X1, == t127.74/tune_p2_esr0.2.json c_in_r); "
          f"C_in(0) {c_in_0!r} F (== t127.74/tune_p2.json c_in)", flush=True)

    # --- templates, checked as item 421 / 420 check them ------------------------
    src = json.loads((OUT_DIR / "t39_esr.json").read_text())
    samples = src["Solver"]["Driven"]["Samples"]
    assert len(samples) == 1 and samples[0]["Freq"] == [0.064], samples
    for p in src["Boundaries"]["LumpedPort"]:
        if int(p["Index"]) == 1:
            assert p["R"] == 50.0 and p["Excitation"] not in (False, 0), p
        else:
            assert p["Excitation"] is False and set(p) >= {"R", "C"}, p
    src0 = json.loads((OUT_DIR / "t39_lossless.json").read_text())
    assert src0["Solver"]["Driven"]["Samples"][0]["Freq"] == [0.064], src0["Solver"]
    for p in src0["Boundaries"]["LumpedPort"]:
        if int(p["Index"]) != 1:
            assert "C" in p and "R" not in p and p["Excitation"] is False, p

    for name, f_ghz, _fh in T42_SWEEP:
        out = json.loads(json.dumps(src))
        out["Solver"]["Driven"]["Samples"][0]["Freq"] = [f_ghz]
        omega = 2.0 * np.pi * f_ghz * 1.0e9
        _q, r_p, c_p = _t39_map(omega, c_in_r, T39_ESR_OHM)
        for p in out["Boundaries"]["LumpedPort"]:
            if int(p["Index"]) == 1:
                continue
            p["R"] = r_p
            p["C"] = c_p
        out["Problem"]["Output"] = f"postpro_{name}"
        path = OUT_DIR / f"{name}.json"
        path.write_text(json.dumps(out, indent=2) + "\n")
        print(f"{tag} wrote {path.relative_to(REPO)}; Freq {[f_ghz]} GHz; "
              f"Boundaries.LumpedPort: " + json.dumps(out["Boundaries"]["LumpedPort"]), flush=True)
    name0, f0_ghz, _f0h = T42_LOSSLESS
    out = json.loads(json.dumps(src0))
    out["Solver"]["Driven"]["Samples"][0]["Freq"] = [f0_ghz]
    for p in out["Boundaries"]["LumpedPort"]:
        if int(p["Index"]) != 1:
            p["C"] = c_in_0
    out["Problem"]["Output"] = f"postpro_{name0}"
    path = OUT_DIR / f"{name0}.json"
    path.write_text(json.dumps(out, indent=2) + "\n")
    print(f"{tag} wrote {path.relative_to(REPO)}; Freq {[f0_ghz]} GHz; "
          f"Boundaries.LumpedPort: " + json.dumps(out["Boundaries"]["LumpedPort"]), flush=True)

    # --- (Z') / (Z0), read back; template diffs ---------------------------------
    fails = []

    def strip(d, drop_r):
        d = json.loads(json.dumps(d))
        d["Solver"]["Driven"]["Samples"][0].pop("Freq")
        d["Problem"].pop("Output")
        for p in d["Boundaries"]["LumpedPort"]:
            if int(p["Index"]) != 1:
                if drop_r:
                    p.pop("R")
                p.pop("C")
        return d

    for name, _fg, f_hz in T42_SWEEP:
        cfg = json.loads((OUT_DIR / f"{name}.json").read_text())
        freq_json = cfg["Solver"]["Driven"]["Samples"][0]["Freq"]
        assert len(freq_json) == 1, freq_json
        omega_json = 2.0 * np.pi * float(freq_json[0]) * 1.0e9
        d_omega = abs(omega_json - 2.0 * np.pi * f_hz) / (2.0 * np.pi * f_hz)
        if not d_omega <= STRIP_BAND:
            fails.append(f"{name}: |d omega|/omega {d_omega:.3e} > {STRIP_BAND}")
        z_ser = series_rlc_impedance(f_hz, r_ohm=T39_ESR_OHM, c_f=c_in_r)
        for p in cfg["Boundaries"]["LumpedPort"]:
            if int(p["Index"]) == 1:
                if p["R"] != 50.0 or p["Excitation"] in (False, 0) or "C" in p:
                    fails.append(f"{name}: port 1 {p}")
                continue
            if p["Excitation"] is not False:
                fails.append(f"{name}: port {p['Index']} Excitation {p['Excitation']!r}")
            z_par = 1.0 / (1.0 / float(p["R"]) + 1j * omega_json * float(p["C"]))
            z_rel = abs(z_par - z_ser) / abs(z_ser)
            if not z_rel <= STRIP_BAND:
                fails.append(f"(Z') {name} port {p['Index']}: {z_rel:.3e} > {STRIP_BAND}")
            if int(p["Index"]) == 2:
                q, r_p, c_p = _t39_map(omega_json, c_in_r, T39_ESR_OHM)
                w, l = 7.0e-3, 8.0e-3
                print(f"{tag} (Z') {name}: Freq {freq_json[0]} GHz (|d omega|/omega "
                      f"{d_omega:.3e}); Q_C {q:.8f}, R_p {float(p['R']):.8f} ohm, C_p "
                      f"{float(p['C']):.10e} F (C_p/C - 1 {float(p['C']) / c_in_r - 1:.4e}); "
                      f"Z_par {_fmt(z_par)} ohm, Z_ser {_fmt(z_ser)} ohm (|Z_ser| {abs(z_ser):.4f}); "
                      f"|Z_par - Z_ser|/|Z_ser| {z_rel:.3e} vs STRIP_BAND {STRIP_BAND} (predicted "
                      f"~1e-16); predicted sheet Rs {float(p['R']) * w / l:.3e} ohm/sq, Cs "
                      f"{float(p['C']) * l / w:.3e} F/sq", flush=True)
        if strip(src, True) != strip(cfg, True):
            fails.append(f"{name}: differs from t39_esr.json beyond Freq / ports 2-4 R, C / Output")
    cfg0 = json.loads((OUT_DIR / f"{name0}.json").read_text())
    fq0 = cfg0["Solver"]["Driven"]["Samples"][0]["Freq"]
    d0 = abs(fq0[0] * 1.0e9 - T42_LOSSLESS[2]) / T42_LOSSLESS[2]
    if not (len(fq0) == 1 and d0 <= STRIP_BAND):
        fails.append(f"{name0}: Freq {fq0}")
    for p in cfg0["Boundaries"]["LumpedPort"]:
        if int(p["Index"]) == 1:
            if p["R"] != 50.0 or p["Excitation"] in (False, 0) or "C" in p:
                fails.append(f"{name0}: port 1 {p}")
            continue
        if "R" in p or p["C"] != T42_C_IN_0_LITERAL or p["Excitation"] is not False:
            fails.append(f"(Z0) {name0} port {p['Index']}: {p}")
    if strip(src0, False) != strip(cfg0, False):
        fails.append(f"{name0}: differs from t39_lossless.json beyond Freq / ports 2-4 C / Output")
    print(f"{tag} (Z0) {name0}: Freq {fq0} GHz, ports 2-4 C == C_in(0) {T42_C_IN_0_LITERAL!r} F bit "
          f"for bit, no R key; lossless Cs {T42_C_IN_0_LITERAL * 8.0 / 7.0:.3e} F/sq (predicted)",
          flush=True)
    # (N-map), printed: the 127.74 MHz map misused at 127.64 / 127.84 MHz
    _q, r_pc, c_pc = _t39_map(2.0 * np.pi * 127.74e6, c_in_r, T39_ESR_OHM)
    for _n, _fg, f_hz in (T42_SWEEP[0], T42_SWEEP[2]):
        z_ser = series_rlc_impedance(f_hz, r_ohm=T39_ESR_OHM, c_f=c_in_r)
        z_mis = 1.0 / (1.0 / r_pc + 1j * 2.0 * np.pi * f_hz * c_pc)
        print(f"{tag} (N-map, printed) 127.74-MHz map at {f_hz / 1e6:.2f} MHz: |Z_par(127.74-map) - "
              f"Z_ser|/|Z_ser| {abs(z_mis - z_ser) / abs(z_ser):.3e} (predicted 3.84e-6; "
              f"STRIP_BAND {STRIP_BAND})", flush=True)
    print(f"{tag} (Z')/(Z0) anchors on {_t42_names()}: "
          f"{'PASS' if not fails else 'FAIL ' + '; '.join(fails)}", flush=True)
    assert not fails, "; ".join(fails)


def read_ans7_128():
    """Item 436 W6: Palace's f0 / R_in / Q at 127.74 MHz via parallel_rows against T35's rows A."""
    from tests.validation.test_ans7_exact_strip import (
        T21P_ZERO_RESOLUTION,
        T22P_ZERO_RESOLUTION,
        T35_C_IN_R_X1,
        T35_Q_R_X1,
        _ANS7,
        _z_branch,
        parallel_rows,
    )

    tag = "[OPS-61 3c W6]"
    c_in_r = float(T35_C_IN_R_X1)
    freqs = np.array([fh for _n, _fg, fh in T42_SWEEP])

    def load(name):
        post = OUT_DIR / f"postpro_{name}"
        csv_path = post / "port-S.csv"
        if not csv_path.exists():
            print(f"{tag} {name}: not measured ({csv_path.relative_to(REPO)} absent)", flush=True)
            return None, None
        s, _hdr = _read_port_s(csv_path)
        pj = post / "palace.json"
        tot = (json.loads(pj.read_text())["ElapsedTime"]["Durations"]["Total"]
               if pj.exists() else None)
        return complex(s[0, 0]), tot

    s11, totals = [], []
    for name, _fg, _fh in T42_SWEEP:
        a, t = load(name)
        s11.append(a)
        totals.append(t)
    s11_0, tot_0 = load(T42_LOSSLESS[0])
    fails = []
    if any(x is None for x in s11):
        print(f"{tag} a point is not measured: (C1) cannot be read", flush=True)
        assert False, "point not measured"

    def rel(a, b):
        return abs(a - b) / abs(b)

    z_bare, z_in, z_in_nofeed_r = [], [], []
    for k, (name, _fg, f_hz) in enumerate(T42_SWEEP):
        omega = 2.0 * np.pi * f_hz
        zb = complex(_ANS7._z_from_s11(s11[k]))
        zbr = _z_branch(omega, c_in_r, T39_ESR_OHM)
        z_bare.append(zb)
        z_in.append(1.0 / (1.0 / zb + 1.0 / zbr))
        z_in_nofeed_r.append(1.0 / (1.0 / zb + 1j * omega * c_in_r))
        tie = rel(s11[k], T42_S11_A[k])
        tot = totals[k]
        print(f"{tag} {f_hz / 1e6:.2f} MHz ({name}): S11_Pal {_fmt(s11[k])}, ours S11_A "
              f"{_fmt(T42_S11_A[k])}; |S11_Pal - S11_A|/|S11_A| {tie:.3e} (printed; predicted <= "
              f"{T39_TIE_PREDICTED:.0e}); Z_bare {_fmt(zb)} ohm; Palace Total "
              f"{tot if tot is None else f'{tot:.3f}'} s", flush=True)
    if s11_0 is not None:
        t60 = rel(s11_0, T42_S11_A_LOSSLESS)
        print(f"{tag} (T6-0') 127.74 MHz lossless C path: S11_Pal {_fmt(s11_0)}, ours "
              f"{_fmt(T42_S11_A_LOSSLESS)}; |S11_Pal - S11_A,lossless|/|S11_A,lossless| {t60:.3e} "
              f"(printed; predicted <= {T39_TIE_PREDICTED:.0e}, convention-free); Palace Total "
              f"{tot_0 if tot_0 is None else f'{tot_0:.3f}'} s", flush=True)
        k_pal = rel(s11[1], s11_0)
        k_ours = rel(T42_S11_A[1], T42_S11_A_LOSSLESS)
        print(f"{tag} (K-esr') Palace esr-vs-lossless move at 127.74 MHz {k_pal:.6e}; ours "
              f"{k_ours:.6e}; |Pal/ours - 1| {abs(k_pal / k_ours - 1.0):.3e} (predicted <= 1e-6; "
              f"a dropped R reads ~1e-10 move)", flush=True)
    rows = parallel_rows(freqs, np.array(z_in))
    b = (1.0 / np.array(z_in)).imag
    print(f"{tag} Palace Im Y_in at (127.64, 127.74, 127.84) MHz = {b[0]:+.9e}, {b[1]:+.9e}, "
          f"{b[2]:+.9e} S; ours {T42_IM_Y_A[0]:+.9e}, {T42_IM_Y_A[1]:+.9e}, "
          f"{T42_IM_Y_A[2]:+.9e} S", flush=True)
    for k in range(3):
        print(f"{tag}   Im Y_in rel gap at {freqs[k] / 1e6:.2f} MHz: "
              f"{abs(b[k] - T42_IM_Y_A[k]) / abs(T42_IM_Y_A[k]):.3e} (printed)", flush=True)
    print(f"{tag} rows (Palace): {json.dumps(rows)}", flush=True)
    z_in_ours = [1.0 / (1.0 / complex(_ANS7._z_from_s11(T42_S11_A[k]))
                        + 1.0 / _z_branch(2.0 * np.pi * freqs[k], c_in_r, T39_ESR_OHM))
                 for k in range(3)]
    rows_o = parallel_rows(freqs, np.array(z_in_ours))
    print(f"{tag} ours' rows from the 12-digit S11_A prints: f0 {rows_o.get('f0_hz')!r} Hz, "
          f"R_in {rows_o.get('r_in_ohm')!r} ohm, Q {rows_o.get('q_loaded')!r} (vs logged "
          f"{T42_F0_A!r} / {T42_R_IN_A!r} / {T42_Q_A!r}; T35_Q_R_X1 {T35_Q_R_X1!r}, item 387)",
          flush=True)

    # --- (C1) one crossing ------------------------------------------------------
    up, down = rows["crossings_up_hz"], rows["crossings_down_hz"]
    c1 = len(up) == 1 and not down and 127.64e6 < up[0] < 127.84e6
    print(f"{tag} (C1) up-crossings {up}, down-crossings {down}: "
          f"{'PASS' if c1 else 'FAIL'} (ASSERTED: exactly one - -> +, none + -> -, inside "
          f"(127.64, 127.84) MHz)", flush=True)
    if not c1:
        fails.append("(C1)")
    else:
        f0, r_in, q = rows["f0_hz"], rows["r_in_ohm"], rows["q_loaded"]
        df0 = f0 / T42_F0_A - 1.0
        dq = q / T42_Q_A - 1.0
        drin = r_in / T42_R_IN_A - 1.0
        f0_ok = abs(df0) <= T21P_ZERO_RESOLUTION
        q_ok = abs(dq) <= T40_Q_TARGET
        print(f"{tag} (f0) f0,Pal {f0!r} Hz vs f0,A {T42_F0_A!r} Hz: f0,Pal/f0,A - 1 = "
              f"{df0:+.3e} (ASSERTED |.| <= T21P_ZERO_RESOLUTION {T21P_ZERO_RESOLUTION:.3e}: "
              f"{'PASS' if f0_ok else 'FAIL'}; printed beside T22P_ZERO_RESOLUTION "
              f"{T22P_ZERO_RESOLUTION:.3e}: {'within' if abs(df0) <= T22P_ZERO_RESOLUTION else 'outside'}; "
              f"predicted <= 1e-6); f0,Pal/127.74 MHz - 1 = {f0 / 127.74e6 - 1.0:+.6e} "
              f"(ours -2.301805e-06)", flush=True)
        print(f"{tag} (Q) Q_Pal {q!r} vs Q_A {T42_Q_A!r}: Q_Pal/Q_A - 1 = {dq:+.3e} (ASSERTED "
              f"|.| <= T40_Q_TARGET {T40_Q_TARGET} (§10 T40's 5 %): {'PASS' if q_ok else 'FAIL'}; "
              f"predicted <= 1e-6); Q_Pal/T35_Q_R_X1 - 1 = {q / T35_Q_R_X1 - 1.0:+.3e} (printed)",
              flush=True)
        print(f"{tag} R_in,Pal {r_in!r} ohm vs R_in,A {T42_R_IN_A!r} ohm: R_in,Pal/R_in,A - 1 = "
              f"{drin:+.3e} (printed; predicted <= 1e-6)", flush=True)
        if not f0_ok:
            fails.append("(f0)")
        if not q_ok:
            fails.append("(Q)")

    # --- negative controls (rule (e); predicted, printed) -----------------------
    rows_nf = parallel_rows(freqs, np.array(z_in_nofeed_r))
    if "q_loaded" in rows_nf:
        qn = rows_nf["q_loaded"]
        in_band = 0.20 <= qn / T42_Q_A - 1.0 <= 0.70
        print(f"{tag} (N-feed) feed's r omitted: Q' {qn!r}, Q'/Q_A - 1 = {qn / T42_Q_A - 1.0:+.3e} "
              f"(predicted in [+20 %, +70 %]: {'met' if in_band else 'NOT met'}; ceiling Q_0 "
              f"618.34 -> +143 %); f0' {rows_nf['f0_hz']!r} Hz (printed)", flush=True)
    else:
        print(f"{tag} (N-feed) feed's r omitted: no - -> + crossing ({json.dumps(rows_nf)}) "
              f"(predicted a crossing with Q' in [+20 %, +70 %]: NOT met)", flush=True)
    rows_nb = parallel_rows(freqs, np.array(z_bare))
    bb = (1.0 / np.array(z_bare)).imag
    if "f0_hz" in rows_nb:
        print(f"{tag} (N-nofeed) feed branch dropped: crossing at {rows_nb['f0_hz']!r} Hz "
              f"(predicted no crossing: NOT met); Im Y_bare {bb[0]:+.4e}, {bb[1]:+.4e}, "
              f"{bb[2]:+.4e} S", flush=True)
    else:
        print(f"{tag} (N-nofeed) feed branch dropped: no crossing / f0 undefined (predicted: met); "
              f"Im Y_bare {bb[0]:+.4e}, {bb[1]:+.4e}, {bb[2]:+.4e} S; rows {json.dumps(rows_nb)}",
              flush=True)
    print(f"{tag} anchors (C1)/(f0)/(Q): {'PASS' if not fails else 'FAIL ' + ', '.join(fails)}; "
          f"ties and R_in printed, none asserted (ruling 2026-10-13 10:30 (4))", flush=True)
    assert not fails, ", ".join(fails)


def export_exact_strip_h0050():
    """W1 of item 437: the h = 0.005 exact-strip mesh + four single-excitation configs.

    A COPY of ``export_exact_strip_h0075`` (commit 81e3cc17, lines 1743-1980 of
    this file) with the H0050 constants substituted; that function is neither
    edited nor called. Its own provenance note follows.

    The sheet steps, the N_i strip selection, the triangle / tet blocks and the
    Gmsh 2.2 writer are COPIED from ``export(exact_strip=True)`` (commit
    4989ebdd, lines 143-383 of this file); ``export()`` is never called or
    edited. Never overwrites. Anchors (E1)-(E4) of item 418.
    """
    import dolfinx
    from mpi4py import MPI

    from tests.mesh.helpers import global_cell_tag_set
    from tests.mesh.test_birdcage_leg_offset import (
        _build,
        _projected_extents,
        _sheet_areas,
        _sheet_azimuth_deg,
    )
    from tests.mesh.test_birdcage_port_sheets import SHEET_IFACE
    from tests.mesh.test_birdcage_port_tags import LEG_COUNT
    from tests.validation.test_ans4_exact_strip_records import (
        SPEC_AREA_M2,
        SPEC_HEIGHT_M,
        SPEC_WIDTH_M,
    )
    from tests.validation.test_ops61_exact_strip_ports import STRIP_BAND
    from tests.validation.test_port_birdcage_leg_offset_sweep import (
        _narrowed_radial,
        _port_frame,
    )
    from tests.validation.test_port_lumped_narrowed_sheet import GATED_WIDTH_FRACTION

    tag = "[OPS-61 1i W1]"
    comm = MPI.COMM_WORLD
    assert comm.size == 1, "export is a serial writer: run it at mpiexec -n 1"
    t_start = time.perf_counter()
    ports_idx = list(range(1, LEG_COUNT + 1))

    # 1. never overwrite
    mesh_path = OUT_DIR / H0050_MESH_NAME
    cfg_paths = {name: OUT_DIR / f"{name}.json" for name in H0050_SET}
    for p in [mesh_path, *cfg_paths.values()]:
        assert not p.exists(), f"{p.relative_to(REPO)} exists: never overwrite"

    # 2. the template, checked the way o2_configs checks its own
    src = json.loads((OUT_DIR / "f05x_o2_f10.json").read_text())
    assert src["Model"]["Mesh"] == STRIP_MESH_NAME, src["Model"]["Mesh"]
    ports = src["Boundaries"]["LumpedPort"]
    assert len(ports) == 4, f"{len(ports)} ports in f05x_o2_f10.json"
    for p in ports:
        i = int(p["Index"])
        assert p["Attributes"] == [NARROW_BASE + i], (i, p["Attributes"])
        assert p["Direction"] == "+Z", (i, p["Direction"])
        assert float(p["R"]) == 50.0, (i, p["R"])
    assert src["Solver"]["Order"] == 2, src["Solver"]["Order"]
    samples = src["Solver"]["Driven"]["Samples"]
    assert len(samples) == 1 and samples[0]["Freq"] == [0.010], samples

    # 3. the build entry 30 measured (exact_rungs under RESOLUTION=0.005)
    msh, cell_tags, _facet_tags, diag, t_mesh = _build(
        np.zeros(LEG_COUNT),
        resolution=H0050_RES,
        leg_sheet_strip_fraction=GATED_WIDTH_FRACTION,
    )
    tdim = msh.topology.dim
    fdim = tdim - 1
    ncells = int(msh.topology.index_map(tdim).size_global)
    print(
        f"{tag} mesh built in {t_mesh:.2f} s: {ncells} tets (record "
        f"{H0050_CELL_RECORD}); resolution {H0050_RES}; diag port_box_size_m = "
        f"{diag['port_box_size_m']}",
        flush=True,
    )
    # (E1) before anything is written
    assert ncells == H0050_CELL_RECORD, (
        f"(E1) W1 count {ncells} != {H0050_CELL_RECORD}: stop before any Palace window"
    )
    msh.topology.create_connectivity(fdim, tdim)
    msh.topology.create_entity_permutations()

    # 4. --- copied from export(exact_strip=True): sheet steps -----------------
    tags_f, full_areas = _sheet_areas(msh, cell_tags, ports_idx, comm)
    dx, _dy, dz = diag["port_box_size_m"]
    sheet_analytic = float(dx) * float(dz)
    geometry = {}
    for i in ports_idx:
        stag = SHEET_IFACE + i
        azimuth = _sheet_azimuth_deg(msh, tags_f, stag, comm)
        radial, azimuthal, axial = _port_frame(azimuth)
        spans, centres = _projected_extents(
            msh, tags_f, stag, comm, [radial, azimuthal, axial]
        )
        geometry[i] = {
            "tag": stag,
            "azimuth_deg": float(azimuth),
            "radial": radial,
            "frame": [radial, azimuthal, axial],
            "w_full": float(spans[0]),
            "h_full": float(spans[2]),
            "centre_radial": float(centres[0]),
            "area_full": float(full_areas[i]),
            "area_ratio_full": float(full_areas[i] / sheet_analytic),
        }
    tags_full = tags_f
    for i in ports_idx:
        g = geometry[i]
        tags_f = _narrowed_radial(
            msh, tags_f, g["tag"], GATED_WIDTH_FRACTION, g["centre_radial"],
            g["radial"], 0.5 * g["w_full"],
        )

    x = np.asarray(msh.geometry.x, dtype=float)

    def facet_nodes(facets):
        facets = np.asarray(facets, dtype=np.int32)
        if not facets.size:
            return np.zeros((0, 3), dtype=np.int64)
        return np.asarray(
            dolfinx.cpp.mesh.entities_to_geometry(msh._cpp_object, fdim, facets, False)
        ).reshape(-1, 3)

    tri_blocks = []
    ext = dolfinx.mesh.exterior_facet_indices(msh.topology)
    ext_nodes = facet_nodes(ext)
    ext_area = float(np.sum(_tri_areas(x, ext_nodes)))
    tri_blocks.append((ext_nodes, EXTERIOR_TAG))
    strips = {}
    for i in ports_idx:
        stag = SHEET_IFACE + i
        full = np.asarray(tags_full.find(stag), dtype=np.int32)
        narrow = np.asarray(tags_f.find(stag), dtype=np.int32)
        subset = bool(np.all(np.isin(narrow, full)))
        rest = np.setdiff1d(full, narrow)
        n_nodes, r_nodes = facet_nodes(narrow), facet_nodes(rest)
        a_n = float(np.sum(_tri_areas(x, n_nodes)))
        pts = x[np.unique(n_nodes.ravel())]
        proj = pts @ np.asarray(geometry[i]["frame"], dtype=float).T
        span = proj.max(axis=0) - proj.min(axis=0)
        h = float(span[2])
        strips[i] = {"area": a_n, "h": h, "w": a_n / h, "w_bbox": float(span[0]),
                     "subset": subset, "n": int(narrow.size), "rest": int(rest.size)}
        print(
            f"    P{i}: azimuth {geometry[i]['azimuth_deg']:8.4f} deg; F_{i} "
            f"{full.size} tris (F/dxdz {geometry[i]['area_ratio_full']:.12f}); N_{i} "
            f"{narrow.size} tris area {a_n:.12e} (rel {a_n / SPEC_AREA_M2 - 1:+.3e}); "
            f"h {h:.12e} (rel {h / SPEC_HEIGHT_M - 1:+.3e}); A/h {a_n / h:.12e} "
            f"(rel {a_n / h / SPEC_WIDTH_M - 1:+.3e}); radial span {span[0]:.12e}; "
            f"F\\N {rest.size} tris; N subset F: {subset}",
            flush=True,
        )
        tri_blocks.append((n_nodes, NARROW_BASE + i))
        tri_blocks.append((r_nodes, REMAINDER_BASE + i))
    ext_err = abs(ext_area - EXTERIOR_AREA_M2)
    print(
        f"{tag} exterior: {ext.size} tris, area {ext_area:.12e} m^2 vs "
        f"{EXTERIOR_AREA_M2} (|diff| {ext_err:.3e}, band {EXTERIOR_AREA_BAND:.0e})",
        flush=True,
    )

    # --- copied: Gmsh 2.2 ASCII writer -----------------------------------------
    OUT_DIR.mkdir(parents=True, exist_ok=True)
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
    total_vol = float(np.sum(np.abs(vol6)) / 6.0)
    ntri = sum(b.shape[0] for b, _ in tri_blocks)
    with open(mesh_path, "w") as fh:
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
    lo, hi = x.min(axis=0), x.max(axis=0)
    tags = sorted(global_cell_tag_set(msh, cell_tags))
    counts = {t: int(np.count_nonzero(ctag == t)) for t in tags}
    print(
        f"{tag} wrote {mesh_path.relative_to(REPO)}: {x.shape[0]} nodes, {ncl} tets "
        f"({int(neg.sum())} reoriented), {ntri} tris; volume {total_vol:.9e} m^3; "
        f"bbox lo {lo} hi {hi}",
        flush=True,
    )
    print(f"{tag} cell tags (count): {counts}", flush=True)
    template_air = sorted(src["Domains"]["Materials"][2]["Attributes"])
    print(f"{tag} template air attributes {template_air}; mesh tags {tags}", flush=True)

    # 5. four single-excitation configs from the template
    for name, k in H0050_SET.items():
        out = json.loads(json.dumps(src))
        out["Model"]["Mesh"] = H0050_MESH_NAME
        for p in out["Boundaries"]["LumpedPort"]:
            p["Excitation"] = int(p["Index"]) == k
        out["Problem"]["Output"] = f"postpro_{name}"
        cfg_paths[name].write_text(json.dumps(out, indent=2) + "\n")
        print(f"{tag} wrote {cfg_paths[name].relative_to(REPO)}: excitation port {k}, "
              f"Output postpro_{name}", flush=True)
    print(json.dumps(json.loads(cfg_paths["f05x_h0050_o2_f10_e1"].read_text())
                     ["Boundaries"], indent=2), flush=True)
    print(f"{tag} export total {time.perf_counter() - t_start:.2f} s (mesh {t_mesh:.2f} s)",
          flush=True)

    # 6. (E2)-(E4) after writing, so a failure leaves the artefacts to inspect
    fails = []
    if not ext_err <= EXTERIOR_AREA_BAND:
        fails.append(f"(E2) exterior area {ext_area:.12e}")
    for i, s in strips.items():
        if not abs(s["area"] / SPEC_AREA_M2 - 1.0) <= STRIP_BAND:
            fails.append(f"(E3) N_{i} area {s['area']:.12e}")
        if not abs(s["w"] / SPEC_WIDTH_M - 1.0) <= STRIP_BAND:
            fails.append(f"(E3) N_{i} A/h {s['w']:.12e}")
        if not s["subset"]:
            fails.append(f"(E4) N_{i} not a subset of F_{i}")
    print(
        f"{tag} ANCHORS (E1) tets == {H0050_CELL_RECORD}: PASS; (E2) exterior area "
        f"within {EXTERIOR_AREA_BAND:.0e}; (E3) N_i area / A/h within STRIP_BAND "
        f"{STRIP_BAND:.0e}; (E4) N_i subset F_i: "
        f"{'PASS' if not fails else 'FAIL ' + '; '.join(fails)}",
        flush=True,
    )
    assert not fails, "; ".join(fails)


def read_exact_strip_h0050():
    """W6 of item 437: a copy of ``read_exact_strip_h0075`` (commit 81e3cc17, lines
    2054-2182; not edited) on the h = 0.005 rung. (P) and (Tie) asserted; (T5),
    (C-h), (C2-h) and the price printed, never gated."""
    from tests.validation.test_ans4_exact_strip_records import ANS4_EXACT_STRIP_S_P2_X1
    from tests.validation.test_port_birdcage_four_port import _circulant_classes

    tag = "[OPS-61 1i W6]"

    def means(m):
        return {c: complex(np.mean(v)) for c, v in _circulant_classes(m).items()}

    def rel(a, b):
        return {c: abs(a[c] - b[c]) / abs(b[c]) for c in b}

    def line(d):
        return "  ".join(f"{c} {v * 100:.4f}% ({v:.3e})" for c, v in d.items())

    # --- ours, parsed; (Tie) ---------------------------------------------------
    s_ours, hdr_line = _ours_rung_s(H0050_OURS_LOG, H0050_OURS_HEADER)
    sig_ours = float(np.linalg.svd(s_ours, compute_uv=False).max())
    tie = abs(sig_ours - H0050_OURS_SIGMA_MAX)
    print(
        f"{tag} ours parsed from {H0050_OURS_LOG.relative_to(REPO)}:{hdr_line} "
        f"(+4 lines): sigma_max {sig_ours:.12f} vs printed {H0050_OURS_SIGMA_MAX} "
        f"|diff| {tie:.3e} (Tie band {H0075_TIE_BAND:.0e}): "
        f"{'PASS' if tie <= H0075_TIE_BAND else 'FAIL'}",
        flush=True,
    )
    for r in s_ours:
        print("    ours " + "  ".join(_fmt(z) for z in r), flush=True)
    assert tie <= H0075_TIE_BAND, f"(Tie) |diff| {tie:.3e}: wrong block parsed"

    # --- Palace columns (h = 0.005, and h = 0.0075 for (C-h)) -------------------
    def assemble(names, label):
        s = np.full((4, 4), np.nan + 0j)
        present = []
        for name, k in names.items():
            path = OUT_DIR / f"postpro_{name}" / "port-S.csv"
            if not path.exists():
                print(f"{tag} {name}: no {path.relative_to(REPO)} (not measured)",
                      flush=True)
                continue
            sk, header = _read_port_s(path)
            cols = [c for c in range(sk.shape[1]) if not np.all(np.isnan(sk[:, c]))]
            labels = [h for h in header if re.search(r"S\[\d+\]\[\d+\]", h)]
            print(f"{tag} {name}: S columns present {[c + 1 for c in cols]}; labels "
                  f"{labels[:2]} ...", flush=True)
            if (k - 1) in cols and sk.shape[0] == 4:
                col = sk[:, k - 1]
                print(f"{tag} {name}: kept column labelled S[i][{k}] (asserted)",
                      flush=True)
            else:
                assert len(cols) == 1 and sk.shape[0] == 4, (name, cols, sk.shape)
                col = sk[:, cols[0]]
                print(f"{tag} {name}: READING -- Palace labelled the single excitation "
                      f"S[i][{cols[0] + 1}], kept by position as column {k}", flush=True)
            assert not np.any(np.isnan(col)), (name, col)
            s[:, k - 1] = col
            present.append(k)
        print(f"{tag} Palace {label} 4x4 (columns present {present}):", flush=True)
        for r in s:
            print("    Palace " + "  ".join(_fmt(z) for z in r), flush=True)
        return s, present

    s, present = assemble(H0050_SET, "h=0.005")

    # --- (T5) per entry of the present columns ----------------------------------
    for k in present:
        g = np.abs(s[:, k - 1] - s_ours[:, k - 1]) / np.abs(s_ours[:, k - 1])
        print(f"{tag} (T5) column {k} per entry |Palace - ours|/|ours|: "
              + "  ".join(f"{v * 100:.4f}% ({v:.3e})" for v in g), flush=True)

    fails = []
    if len(present) == 4:
        ident, fails = _p_check_s("f05x_h0050_o2_f10 (assembled)", s, tag)
        pal = means(s)
        ours_h = means(s_ours)
        t5 = rel(pal, ours_h)
        print(f"{tag} (T5) per class |Palace o2 h - ours p2 h|/|ours p2 h|: {line(t5)} "
              "(predicted <= 1e-6; printed, never gated)", flush=True)
        s75, p75 = assemble(H0075_SET, "h=0.0075")
        s_ours75, hdr75 = _ours_h0075_s()
        ours75 = means(s_ours75)
        ours_move = rel(ours_h, ours75)
        print(f"{tag} ours' 0.0075 -> 0.005 move, computed from "
              f"{H0075_OURS_LOG.name}:{hdr75} and {H0050_OURS_LOG.name}:{hdr_line}: "
              f"{line(ours_move)} (drafter's arithmetic 0.0144 / 0.0088 / 0.0161 %)",
              flush=True)
        if len(p75) == 4:
            pal75 = means(s75)
            ch = rel(pal, pal75)
            print(f"{tag} (C-h) Palace o2 h=0.005 vs Palace o2 h=0.0075 per class: "
                  f"{line(ch)}; beside ours' same move {line(ours_move)} (predicted "
                  "equal to the printed digit; printed)", flush=True)
        else:
            print(f"{tag} (C-h) not measured: h=0.0075 columns present {p75}", flush=True)
        num = {c: abs(ours75[c] - pal[c]) for c in pal}
        den = {c: abs(ours_h[c] - pal[c]) for c in pal}
        print(
            f"{tag} (C2-h) |ours h=0.0075 p2 - Palace o2 h=0.005| / |ours h=0.005 p2 - "
            "Palace o2 h=0.005| per class: " + "  ".join(
                f"{c} {num[c]:.3e}/{den[c]:.3e} = "
                f"{(num[c] / den[c]) if den[c] > 0 else float('inf'):.4g}x" for c in pal)
            + " (predicted >= 88x, vs item 229's 25x bar; printed)",
            flush=True,
        )
        ours_x1 = means(ANS4_EXACT_STRIP_S_P2_X1["10"])
        num1 = {c: abs(ours_x1[c] - pal[c]) for c in pal}
        print(
            f"{tag} (C2-h) against ours x1 p2: " + "  ".join(
                f"{c} {num1[c]:.3e}/{den[c]:.3e} = "
                f"{(num1[c] / den[c]) if den[c] > 0 else float('inf'):.4g}x" for c in pal)
            + " (predicted >= 380x; printed); ours' x1 -> 0.005 move (log :3840) "
            + "  ".join(f"{c} {v * 100:.4f}%" for c, v in H0050_OURS_X1_MOVE.items()),
            flush=True,
        )
    else:
        print(f"{tag} (P), per-class (T5), (C-h), (C2-h): not measured "
              f"(columns present {present})", flush=True)

    # --- price readout from the W2-W5 harness logs --------------------------------
    logs = sorted((REPO / "docs" / "testing" / "logs").glob(
        "*_OPS-61-item437-W[2-5]-palace-h0050-e*.log"))
    r_f = {"64": 37 / 26, "128": 55 / 26}
    for lg in logs:
        txt = lg.read_text()
        # Elapsed Time Report columns Min / Max / Avg: keep Avg (brief).
        tot = re.findall(r"^\s*Total\s+[\d.]+\s+[\d.]+\s+([\d.]+)\s*$", txt, re.M)
        setup = re.findall(r"^\s*Setup\s+[\d.]+\s+[\d.]+\s+([\d.]+)\s*$", txt, re.M)
        its = re.findall(r"GMRES solver converged in (\d+) iteration", txt)
        if not tot or not setup:
            print(f"{tag} price {lg.name}: no Total/Setup line", flush=True)
            continue
        t_e, s_e = float(tot[-1]), float(setup[-1])
        p = {f: s_e + (t_e - s_e) * r for f, r in r_f.items()}
        print(
            f"{tag} price {lg.name}: Total {t_e:.1f} s, Setup {s_e:.1f} s, GMRES its "
            f"{its}; projected single-excitation 64 MHz {p['64']:.1f} s, 128 MHz "
            f"{p['128']:.1f} s (item 418's iteration ratios; printed for the weekly)",
            flush=True,
        )

    assert not fails, "; ".join(fails)
    print(f"{tag} done: (Tie) PASS" + ("; (P) PASS" if len(present) == 4 else ""),
          flush=True)


# --- OPS-61 step 1j (§9 item 458 / 458b, §10 F39): 64 MHz on the h = 0.005 ------
# exact-strip rung.  Additive copies; h0075_freq_configs, _ours_rung_s and
# read_exact_strip_h0050 are not edited.
def h0050_f64_configs():
    """W1 of item 458: (E1') on the h = 0.005 mesh; four single-excitation 64 MHz
    configs from item 437's 10 MHz configs (Freq / Output only).  A copy of
    ``h0075_freq_configs`` (commit 970466b9, lines 2202-2265) with the h = 0.005
    names; the sha256 / byte-size print is dropped (private-leak hook trap)."""
    try:
        from mpi4py import MPI

        assert MPI.COMM_WORLD.size == 1, "h0050-configs runs at -n 1"
    except ImportError:
        pass
    tag = "[OPS-61 h0050 64 MHz W1]"
    fghz = 0.064
    targets = {name: OUT_DIR / f"{name}.json" for name in H0050_F64_SET}
    for p in targets.values():
        assert not p.exists(), f"{p} exists: never overwrite"

    # (E1') tet count parsed from the Gmsh 2.2 $Elements block
    mesh_path = OUT_DIR / H0050_MESH_NAME
    lines = mesh_path.read_bytes().decode().splitlines()
    i0 = lines.index("$Elements")
    i1 = lines.index("$EndElements")
    nel = int(lines[i0 + 1])
    ntet = sum(1 for ln in lines[i0 + 2: i1] if ln.split()[1] == "4")
    assert i1 - i0 - 2 == nel, (nel, i1 - i0 - 2)
    print(f"{tag} {mesh_path.relative_to(REPO)}: {nel} elements, {ntet} type-4 tets",
          flush=True)
    print(f"{tag} (E1') tets {ntet} == {H0050_CELL_RECORD}: "
          f"{'PASS' if ntet == H0050_CELL_RECORD else 'FAIL'}", flush=True)
    assert ntet == H0050_CELL_RECORD, f"(E1') {ntet} != {H0050_CELL_RECORD}: stop"

    for name, k in H0050_F64_SET.items():
        tpl_path = OUT_DIR / f"f05x_h0050_o2_f10_e{k}.json"
        tpl = json.loads(tpl_path.read_text())
        assert tpl["Model"]["Mesh"] == H0050_MESH_NAME, tpl["Model"]["Mesh"]
        ports = tpl["Boundaries"]["LumpedPort"]
        assert [int(p["Index"]) for p in ports] == [1, 2, 3, 4], ports
        for p in ports:
            assert p["Attributes"] == [20 + int(p["Index"])], p
            assert p["Direction"] == "+Z" and float(p["R"]) == 50.0, p
            assert bool(p["Excitation"]) == (int(p["Index"]) == k), p
        assert tpl["Solver"]["Order"] == 2, tpl["Solver"]["Order"]
        samples = tpl["Solver"]["Driven"]["Samples"]
        assert len(samples) == 1 and samples[0]["Freq"] == [0.010], samples
        out = json.loads(json.dumps(tpl))
        out["Solver"]["Driven"]["Samples"][0]["Freq"] = [fghz]
        out["Problem"]["Output"] = f"postpro_{name}"
        targets[name].write_text(json.dumps(out, indent=2) + "\n")

        # template diff: exactly the two keys
        back = json.loads(targets[name].read_text())
        t2 = json.loads(json.dumps(tpl))
        b2 = json.loads(json.dumps(back))
        assert b2["Solver"]["Driven"]["Samples"][0]["Freq"] == [fghz]
        assert b2["Problem"]["Output"] == f"postpro_{name}"
        for d in (t2, b2):
            d["Solver"]["Driven"]["Samples"][0]["Freq"] = None
            d["Problem"]["Output"] = None
        assert t2 == b2, f"{name}: differs from its template beyond Freq / Output"
        print(f"{tag} wrote {targets[name].relative_to(REPO)} from {tpl_path.name}: "
              f"excitation port {k}, Freq [{fghz}] GHz, Output "
              f"postpro_{name}; template diff = Freq, Output only: PASS", flush=True)
    print(f"{tag} done: (E1') PASS; four configs written", flush=True)


def _ours_rung_s_scan(log, header):
    """Gap-tolerant copy of ``_ours_rung_s`` (not edited): after the single header,
    walk at most 8 lines, skip blanks, find ``S_{k}k =`` by re.search (tolerating
    PASSED / pytest node-ID prefixes), k = 1..4 in order, no new ``[`` header first."""
    import gzip

    path = log
    if path.exists():
        lines = path.read_text().splitlines()
    else:
        with gzip.open(str(path) + ".gz", "rt") as fh:
            lines = fh.read().splitlines()
    idx = [n for n, ln in enumerate(lines) if ln.strip() == header]
    assert len(idx) == 1, f"{len(idx)} headers {header!r}"
    tok = re.compile(r"[+-]\d\.\d+e[+-]\d+[+-]\d\.\d+e[+-]\d+j")
    rows = []
    for ln in lines[idx[0] + 1: idx[0] + 9]:
        if len(rows) == 4:
            break
        if not ln.strip():
            continue
        assert not ln.lstrip().startswith("["), ("new header before row 4", ln)
        k = len(rows) + 1
        m = re.search(rf"S_{k}k =", ln)
        assert m is not None, (k, ln)
        vals = tok.findall(ln[m.start():])
        assert len(vals) == 4, (k, ln)
        rows.append([complex(v) for v in vals])
    assert len(rows) == 4, f"{len(rows)} rows within 8 lines of {header!r}"
    return np.array(rows), idx[0] + 1


def read_exact_strip_h0050_f64():
    """W4 of item 458 / W3 of item 458b: a copy of ``read_exact_strip_h0050``
    (commit 970466b9, lines 3417-3562; not edited) at 64 MHz.  (Parse-0), (Tie),
    (Tie75) and, with four columns, (P) asserted; (T5), (C-h), (C2-h) and the
    price printed, never gated."""
    from tests.validation.test_ans4_exact_strip_records import ANS4_EXACT_STRIP_S_P2_X1
    from tests.validation.test_port_birdcage_four_port import _circulant_classes

    tag = "[OPS-61 1j read]"

    def means(m):
        return {c: complex(np.mean(v)) for c, v in _circulant_classes(m).items()}

    def rel(a, b):
        return {c: abs(a[c] - b[c]) / abs(b[c]) for c in b}

    def line(d):
        return "  ".join(f"{c} {v * 100:.4f}% ({v:.3e})" for c, v in d.items())

    # --- (Parse-0): the new parser reproduces the old on the contiguous block ------
    s_old, l_old = _ours_rung_s(H0050_OURS_LOG, H0050_OURS_HEADER)
    s_new, l_new = _ours_rung_s_scan(H0050_OURS_LOG, H0050_OURS_HEADER)
    p0 = bool(l_old == l_new and s_old.shape == s_new.shape
              and np.array_equal(s_old, s_new))
    print(f"{tag} (Parse-0) _ours_rung_s_scan == _ours_rung_s bitwise on the 10 MHz "
          f"block ({H0050_OURS_LOG.name}:{l_old}): {'PASS' if p0 else 'FAIL'}",
          flush=True)
    assert p0, "(Parse-0) new parser differs from _ours_rung_s"

    # --- ours at h = 0.005, 64 MHz; (Tie) ---------------------------------------
    s_ours, hdr_line = _ours_rung_s_scan(H0050_OURS_LOG, H0050_F64_OURS_HEADER)
    sig_ours = float(np.linalg.svd(s_ours, compute_uv=False).max())
    tie = abs(sig_ours - H0050_F64_OURS_SIGMA_MAX)
    print(
        f"{tag} ours parsed from {H0050_OURS_LOG.relative_to(REPO)}:{hdr_line} "
        f"(gap-tolerant): sigma_max {sig_ours:.12f} vs printed "
        f"{H0050_F64_OURS_SIGMA_MAX} |diff| {tie:.3e} (Tie band {H0075_TIE_BAND:.0e}): "
        f"{'PASS' if tie <= H0075_TIE_BAND else 'FAIL'}",
        flush=True,
    )
    for r in s_ours:
        print("    ours " + "  ".join(_fmt(z) for z in r), flush=True)
    assert tie <= H0075_TIE_BAND, f"(Tie) |diff| {tie:.3e}: wrong block parsed"

    # --- ours at h = 0.0075, 64 MHz; (Tie75) ------------------------------------
    b64 = H0075_OURS_BLOCKS["64"]
    s_ours75, hdr75 = _ours_rung_s(H0075_OURS_LOG, b64["header"])
    sig75 = float(np.linalg.svd(s_ours75, compute_uv=False).max())
    tie75 = abs(sig75 - b64["sigma_max"])
    print(f"{tag} ours h=0.0075 parsed from {H0075_OURS_LOG.name}:{hdr75}: sigma_max "
          f"{sig75:.12f} vs printed {b64['sigma_max']} |diff| {tie75:.3e} (Tie75 band "
          f"{H0075_TIE_BAND:.0e}): {'PASS' if tie75 <= H0075_TIE_BAND else 'FAIL'}",
          flush=True)
    assert tie75 <= H0075_TIE_BAND, f"(Tie75) |diff| {tie75:.3e}: wrong block parsed"

    # --- Palace columns (h = 0.005, and h = 0.0075 for (C-h)) -------------------
    def assemble(names, label):
        s = np.full((4, 4), np.nan + 0j)
        present = []
        for name, k in names.items():
            path = OUT_DIR / f"postpro_{name}" / "port-S.csv"
            if not path.exists():
                print(f"{tag} {name}: no {path.relative_to(REPO)} (not measured)",
                      flush=True)
                continue
            sk, header = _read_port_s(path)
            cols = [c for c in range(sk.shape[1]) if not np.all(np.isnan(sk[:, c]))]
            labels = [h for h in header if re.search(r"S\[\d+\]\[\d+\]", h)]
            print(f"{tag} {name}: S columns present {[c + 1 for c in cols]}; labels "
                  f"{labels[:2]} ...", flush=True)
            if (k - 1) in cols and sk.shape[0] == 4:
                col = sk[:, k - 1]
                print(f"{tag} {name}: kept column labelled S[i][{k}] (asserted)",
                      flush=True)
            else:
                assert len(cols) == 1 and sk.shape[0] == 4, (name, cols, sk.shape)
                col = sk[:, cols[0]]
                print(f"{tag} {name}: READING -- Palace labelled the single excitation "
                      f"S[i][{cols[0] + 1}], kept by position as column {k}", flush=True)
            assert not np.any(np.isnan(col)), (name, col)
            s[:, k - 1] = col
            present.append(k)
        print(f"{tag} Palace {label} 4x4 (columns present {present}):", flush=True)
        for r in s:
            print("    Palace " + "  ".join(_fmt(z) for z in r), flush=True)
        return s, present

    s, present = assemble(H0050_F64_SET, "h=0.005 64 MHz")

    # --- (T5) per entry of the present columns ----------------------------------
    for k in present:
        g = np.abs(s[:, k - 1] - s_ours[:, k - 1]) / np.abs(s_ours[:, k - 1])
        print(f"{tag} (T5) column {k} per entry |Palace - ours|/|ours|: "
              + "  ".join(f"{v * 100:.4f}% ({v:.3e})" for v in g)
              + " (predicted <= 1e-6; printed, never gated)", flush=True)

    fails = []
    if len(present) == 4:
        ident, fails = _p_check_s("f05x_h0050_o2_f64 (assembled)", s, tag)
        pal = means(s)
        ours_h = means(s_ours)
        t5 = rel(pal, ours_h)
        print(f"{tag} (T5) per class |Palace o2 h - ours p2 h|/|ours p2 h|: {line(t5)} "
              "(predicted <= 1e-6; printed, never gated)", flush=True)
        sp_p = _identities(s)["spreads"]
        sp_o = _identities(s_ours)["spreads"]
        print(f"{tag} _class_spread per class, Palace / ours: " + "  ".join(
            f"{c} {sp_p[c] * 100:.4f}% / {sp_o[c] * 100:.4f}%" for c in sp_p), flush=True)
        s75, p75 = assemble(H0075_FREQ_SET["64"], "h=0.0075 64 MHz")
        ours75 = means(s_ours75)
        ours_move = rel(ours_h, ours75)
        print(f"{tag} ours' 0.0075 -> 0.005 move at 64 MHz, computed from "
              f"{H0075_OURS_LOG.name}:{hdr75} and {H0050_OURS_LOG.name}:{hdr_line}: "
              f"{line(ours_move)} (drafter's arithmetic 0.0569 / 0.0235 / 0.0611 %)",
              flush=True)
        if len(p75) == 4:
            pal75 = means(s75)
            ch = rel(pal, pal75)
            print(f"{tag} (C-h) Palace o2 h=0.005 vs Palace o2 h=0.0075 per class: "
                  f"{line(ch)}; beside ours' same move {line(ours_move)} (predicted "
                  "equal to the printed digit; printed)", flush=True)
        else:
            print(f"{tag} (C-h) not measured: h=0.0075 columns present {p75}", flush=True)
        num = {c: abs(ours75[c] - pal[c]) for c in pal}
        den = {c: abs(ours_h[c] - pal[c]) for c in pal}
        print(
            f"{tag} (C2-h) |ours h=0.0075 p2 - Palace o2 h=0.005| / |ours h=0.005 p2 - "
            "Palace o2 h=0.005| per class: " + "  ".join(
                f"{c} {num[c]:.3e}/{den[c]:.3e} = "
                f"{(num[c] / den[c]) if den[c] > 0 else float('inf'):.4g}x" for c in pal)
            + " (predicted >= 235x; printed)",
            flush=True,
        )
        ours_x1 = means(ANS4_EXACT_STRIP_S_P2_X1["64"])
        num1 = {c: abs(ours_x1[c] - pal[c]) for c in pal}
        print(
            f"{tag} (C2-h) against ours x1 p2: " + "  ".join(
                f"{c} {num1[c]:.3e}/{den[c]:.3e} = "
                f"{(num1[c] / den[c]) if den[c] > 0 else float('inf'):.4g}x" for c in pal)
            + " (predicted >= 1144x; printed); ours' x1 -> 0.005 move (log :3841) "
            + "  ".join(f"{c} {v * 100:.4f}%" for c, v in H0050_F64_OURS_X1_MOVE.items()),
            flush=True,
        )
    else:
        print(f"{tag} (P), per-class (T5), (C-h), (C2-h): not measured "
              f"(columns present {present})", flush=True)

    # --- price readout from both halves' Palace harness logs ---------------------
    logs = sorted((REPO / "docs" / "testing" / "logs").glob(
        "*_OPS-61-item458*-W*-palace-h0050-f64-e*.log"))
    for lg in logs:
        txt = lg.read_text()
        tot = re.findall(r"^\s*Total\s+[\d.]+\s+[\d.]+\s+([\d.]+)\s*$", txt, re.M)
        setup = re.findall(r"^\s*Setup\s+[\d.]+\s+[\d.]+\s+([\d.]+)\s*$", txt, re.M)
        its = re.findall(r"GMRES solver converged in (\d+) iteration", txt)
        peak = re.findall(r"peak per-node memory usage is:.*Total ([\d.]+[GM])", txt)
        if not tot or not setup:
            print(f"{tag} price {lg.name}: no Total/Setup line (> 590 s or aborted)",
                  flush=True)
            continue
        print(f"{tag} price {lg.name}: Total {float(tot[-1]):.1f} s, Setup "
              f"{float(setup[-1]):.1f} s, GMRES its {its}, peak "
              f"{peak[-1] if peak else 'n/a'}; predicted "
              f"Total {H0050_F64_PROJ[0]:.0f}-{H0050_F64_PROJ[1]:.0f} s", flush=True)

    assert not fails, "; ".join(fails)
    print(f"{tag} done: (Parse-0) PASS; (Tie) PASS; (Tie75) PASS"
          + ("; (P) PASS" if len(present) == 4 else ""), flush=True)


# --- OPS-61 step 3d (§9 item 444, §10 T43, provisional): items 420 + 421 on ----
# T37's h = 0.01 hole rung.  Palace's twin + driven read of the r = 0.2 ohm tune
# at 64.000 MHz on the exact-rectangle hole mesh at resolution 0.01, against
# T37's Route A rows (20261012T094711Z_ANS-7-item395-W1-esr64-h001.log, the
# window that wrote tune_p2_h0.01_esr0.2.json).  C_in(r, h) is the rung's own
# value, never a x1 literal (item 419's silent misread).  Additive: no earlier
# mode is edited or called.
H001_RES = 0.01  # == H_RES, test_ans7_exact_strip.py:1059 (asserted on use)
H001_HOLE_MESH_NAME = "ops61_ans6_h0010_hole_strip.msh"
# Rung mesh cells: ...item395-W1-esr64-h001.log:2203; ...item214-W3-rung.log:2236
# (no tracked constant holds it).
H001_HOLE_CELL_RECORD = 140270
T43_SWEEP = (("t43_esr_f63.9", 0.0639, 63.9e6), ("t43_esr_f64.0", 0.064, 64.0e6),
             ("t43_esr_f64.1", 0.0641, 64.1e6))
T43_LOSSLESS = ("t43_lossless_f64.0", 0.064, 64.0e6)
# ...item395-W1-esr64-h001.log:486 (== tune_p2_h0.01_esr0.2.json c_in_r / c_in_0)
T43_C_IN_R_H_LITERAL = 6.364610922203437e-11
T43_C_IN_0_H_LITERAL = 6.364316664701229e-11
# Ours' S11_A at (r, C_in(r, h)), 12-digit prints, ...item395-W1-esr64-h001.log:2204-2206:
T43_S11_A = (-2.445046675335e-01 + 9.625434769724e-01j,
             -2.400384573551e-01 + 9.635701508391e-01j,
             -2.355297117830e-01 + 9.645820886247e-01j)
# T37 rows A, ...item395-W1-esr64-h001.log:2207:
T43_F0_A = 63999834.331344604
T43_R_IN_A = 3166.7479872416206
T43_Q_A = 164.21575983485775
T43_IM_Y_A = (-1.617828160e-04, +2.684681835e-07, +1.626210783e-04)
T43_Q_H_MOVE_OURS = -1.991409e-03  # ...item395-W1-esr64-h001.log:2211
T43_Q_0_H = 964.432729560342  # ...item395-W1-esr64-h001.log:2210 (printed ceiling only)
# Lossless C path on the rung at C_in(0, h), 64.0 MHz: 20261008T111904Z_ANS-7-item310-W3-p64.0.log:2288
T43_S11_A_LOSSLESS = -2.411870945081e-01 + 9.683792879359e-01j
# Palace's x1 rows (f0, Q, R_in): 20261013T141248Z_OPS-61-item421-W4-read.log:48-50
T43_PAL_X1_ROWS = (63999834.68507022, 164.54343232572575, 3175.979463819766)
T43_GLOBAL_TAGS = (2, 3, 101, 102, 103, 104, 201, 202, 203, 204)  # == t39_esr.json Materials (asserted)


def _t43_names():
    return [n for n, _g, _h in T43_SWEEP] + [T43_LOSSLESS[0]]


def export_hole_h001():
    """Item 444 W1: the h = 0.01 exact-rectangle hole mesh (Gmsh 2.2), mesh only.

    A copy of ``export_hole_pec`` (commit b776becf, lines 1095-1319) with the
    rung's ``resolution=H001_RES`` -- the call ``_exact_geo(th14, comm,
    GATED_WIDTH_FRACTION, res=H_RES)`` makes inside ``tune_h``
    (test_ans7_exact_strip.py:409-425, :2091).  No configs are written.
    """
    import dolfinx
    from mpi4py import MPI

    from tests.mesh.helpers import global_cell_tag_set
    from tests.mesh.test_birdcage_port_sheets import SHEET_IFACE
    from tests.mesh.test_birdcage_port_sheets import _build as _sheets_build
    from tests.mesh.test_birdcage_port_tags import LEG_COUNT
    from tests.validation import test_th14_birdcage_copper as th14
    from tests.validation.test_ans4_exact_strip_records import (
        SPEC_AREA_M2,
        SPEC_HEIGHT_M,
        SPEC_WIDTH_M,
    )
    from tests.validation.test_ans7_exact_strip import H_RES
    from tests.validation.test_ops61_exact_strip_ports import STRIP_BAND
    from tests.validation.test_port_birdcage_lumped_column import (
        CONDUCTOR_CELL_TAG,
        PHANTOM_CELL_TAG,
    )
    from tests.validation.test_port_lumped_narrowed_sheet import GATED_WIDTH_FRACTION

    tag = "[OPS-61 3d W1]"
    comm = MPI.COMM_WORLD
    assert comm.size == 1, "export_hole_h001 is a serial writer: run it at mpiexec -n 1"
    assert H001_RES == H_RES, (H001_RES, H_RES)
    assert th14._sheets_build is _sheets_build
    t_start = time.perf_counter()
    wall_attr = _hole_wall_attr()
    ports_idx = list(range(1, LEG_COUNT + 1))

    # 1. never overwrite
    mesh_path = OUT_DIR / H001_HOLE_MESH_NAME
    assert not mesh_path.exists(), f"{mesh_path.relative_to(REPO)} exists: never overwrite"
    tmpl = json.loads((OUT_DIR / "t39_esr.json").read_text())
    tmpl_tags = sorted(t for m in tmpl["Domains"]["Materials"] for t in m["Attributes"])
    assert tuple(tmpl_tags) == T43_GLOBAL_TAGS, (tmpl_tags, T43_GLOBAL_TAGS)

    # 2. the rung's hole route with exact rectangles, as tune_h builds it.
    msh, cell_tags, facet_tags, diag, t_mesh = _sheets_build(
        True, as_hole=True, leg_sheet_strip_fraction=float(GATED_WIDTH_FRACTION),
        resolution=float(H001_RES),
    )
    tdim = msh.topology.dim
    fdim = tdim - 1
    msh.topology.create_connectivity(fdim, tdim)
    msh.topology.create_entity_permutations()
    ncells = int(msh.topology.index_map(tdim).size_global)
    print(f"{tag} rung hole mesh (resolution {H001_RES}) built in {t_mesh:.2f} s: {ncells} tets "
          f"(record {H001_HOLE_CELL_RECORD})", flush=True)
    # 3. (E1) before anything is written
    assert ncells == H001_HOLE_CELL_RECORD, f"(E1) tets {ncells} != {H001_HOLE_CELL_RECORD}"

    _outer, census = th14._outer_box_tags(msh, facet_tags, comm)
    tags_f, _defs, specs = th14._build_ports(msh, cell_tags, comm)
    print(f"{tag} census {census}", flush=True)

    x = np.asarray(msh.geometry.x, dtype=float)

    def facet_nodes(facets):
        facets = np.asarray(facets, dtype=np.int32)
        if not facets.size:
            return np.zeros((0, 3), dtype=np.int64)
        return np.asarray(
            dolfinx.cpp.mesh.entities_to_geometry(msh._cpp_object, fdim, facets, False)
        ).reshape(-1, 3)

    # 4. triangle blocks: box -> 1, wall -> 401, N_i -> 20 + i.
    ext = np.asarray(dolfinx.mesh.exterior_facet_indices(msh.topology), dtype=np.int32)
    wall = np.asarray(facet_tags.find(wall_attr), dtype=np.int32)
    box = np.setdiff1d(ext, wall).astype(np.int32)
    box_nodes, wall_nodes = facet_nodes(box), facet_nodes(wall)
    box_area = float(np.sum(_tri_areas(x, box_nodes)))
    wall_area = float(np.sum(_tri_areas(x, wall_nodes)))
    tri_blocks = [(box_nodes, EXTERIOR_TAG), (wall_nodes, wall_attr)]
    port_rows = {}
    for i, spec in zip(ports_idx, specs):
        ptag = SHEET_IFACE + i
        assert int(spec.facet_tag) == ptag, (i, spec.facet_tag)
        narrow = np.asarray(tags_f.find(ptag), dtype=np.int32)
        n_nodes = facet_nodes(narrow)
        a_n = float(np.sum(_tri_areas(x, n_nodes)))
        port_rows[i] = {
            "facets": int(narrow.size),
            "area": a_n,
            "width": float(spec.sheet_width_m),
            "height": float(spec.gap_height_m),
        }
        tri_blocks.append((n_nodes, NARROW_BASE + i))

    OUT_DIR.mkdir(parents=True, exist_ok=True)
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
    total_vol = float(np.sum(np.abs(vol6)) / 6.0)

    ntri = sum(b.shape[0] for b, _ in tri_blocks)
    with open(mesh_path, "w") as fh:
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

    # (E1) again, parsed from the written file's $Elements type-4 count
    n_tet_file = 0
    with open(mesh_path) as fh:
        inside = False
        for line in fh:
            if line.startswith("$Elements"):
                inside = True
                next(fh)
                continue
            if line.startswith("$EndElements"):
                break
            if inside and line.split(None, 2)[1] == "4":
                n_tet_file += 1

    lo, hi = x.min(axis=0), x.max(axis=0)
    tags = sorted(global_cell_tag_set(msh, cell_tags))
    counts = {t: int(np.count_nonzero(ctag == t)) for t in tags}
    print(
        f"{tag} wrote {mesh_path.relative_to(REPO)}: {ncl} tets ({int(neg.sum())} reoriented), "
        f"$Elements type-4 count {n_tet_file}; volume {total_vol:.9e} m^3; bbox lo {lo} hi {hi}",
        flush=True,
    )
    print(f"{tag} cell tags (count): {counts}", flush=True)
    print(
        f"{tag} box (attr {EXTERIOR_TAG}): {box.size} tris area {box_area:.12e} m^2; "
        f"wall (attr {wall_attr}): {wall.size} tris area {wall_area:.12e} m^2",
        flush=True,
    )
    for i, r in port_rows.items():
        print(
            f"    P{i}: N_{i} {r['facets']} tris (attr {NARROW_BASE + i}) area "
            f"{r['area']:.12e} m^2 (rel {r['area'] / SPEC_AREA_M2 - 1.0:+.3e}); spec "
            f"sheet_width_m {r['width']:.12e} (rel {r['width'] / SPEC_WIDTH_M - 1.0:+.3e}) "
            f"gap_height_m {r['height']:.12e} (rel {r['height'] / SPEC_HEIGHT_M - 1.0:+.3e})",
            flush=True,
        )
    print(f"{tag} export total {time.perf_counter() - t_start:.2f} s (mesh {t_mesh:.2f} s)",
          flush=True)

    # 5. anchors (E1)-(E6), after writing.
    fails = []
    if ncl != H001_HOLE_CELL_RECORD or ncells != ncl or n_tet_file != H001_HOLE_CELL_RECORD:
        fails.append(f"(E1) tets {ncl} / file {n_tet_file} != {H001_HOLE_CELL_RECORD}")
    if CONDUCTOR_CELL_TAG in tags or PHANTOM_CELL_TAG not in tags:
        fails.append(f"(E2) cell tags {tags}")
    if not (census["cavity_exterior"] == census["cavity"] > 0
            and census["outer"] == census["exterior"] - census["cavity"]
            and census["tag_collision"] == 0):
        fails.append(f"(E3) census {census}")
    if not abs(box_area - EXTERIOR_AREA_M2) <= EXTERIOR_AREA_BAND:
        fails.append(f"(E4) box area {box_area:.12e} vs {EXTERIOR_AREA_M2}")
    for i, r in port_rows.items():
        for k, (got, want) in {"area": (r["area"], SPEC_AREA_M2),
                               "width": (r["width"], SPEC_WIDTH_M),
                               "height": (r["height"], SPEC_HEIGHT_M)}.items():
            if not abs(got / want - 1.0) < STRIP_BAND:
                fails.append(f"(E5) P{i} {k} {got:.12e} vs {want:.12e}")
    if tuple(tags) != tuple(tmpl_tags):
        fails.append(f"(E6) cell-tag set {tags} != t39_esr.json Materials {tmpl_tags}")
    print(f"{tag} (E1) {ncl} == file {n_tet_file} == {H001_HOLE_CELL_RECORD}; (E4) box area "
          f"|d| {abs(box_area - EXTERIOR_AREA_M2):.3e}; (E6) tags {tags} vs template {tmpl_tags}",
          flush=True)
    print(
        f"{tag} anchors (E1)-(E6) [box band {EXTERIOR_AREA_BAND:.0e}, strip band "
        f"{STRIP_BAND}]: {'PASS' if not fails else 'FAIL ' + '; '.join(fails)}",
        flush=True,
    )
    assert not fails, "; ".join(fails)


def ans7_h001_configs():
    """Item 444 W2: write t43_esr_f{63.9,64.0,64.1} / t43_lossless_f64.0; assert (Z') / (Z0)."""
    from fem_em_solver.ports.lumped import series_rlc_impedance
    from tests.validation.test_ops61_exact_strip_ports import STRIP_BAND

    tag = "[OPS-61 3d W2]"
    try:
        from mpi4py import MPI

        assert MPI.COMM_WORLD.size == 1, f"run at -n 1, not {MPI.COMM_WORLD.size}"
    except ImportError:
        pass
    for name in _t43_names():
        path = OUT_DIR / f"{name}.json"
        assert not path.exists(), f"{path.relative_to(REPO)} exists: never overwrite"
    assert (OUT_DIR / H001_HOLE_MESH_NAME).exists(), "W1's rung mesh is absent"
    ans7 = REPO / "logs" / "ans7_exact"
    tr = json.loads((ans7 / "tune_p2_h0.01_esr0.2.json").read_text())
    t0 = json.loads((ans7 / "tune_p2_h0.01.json").read_text())
    assert tr["c_in_r"] == T43_C_IN_R_H_LITERAL, (tr["c_in_r"], T43_C_IN_R_H_LITERAL)
    assert tr["c_in_0"] == T43_C_IN_0_H_LITERAL, (tr["c_in_0"], T43_C_IN_0_H_LITERAL)
    assert tr["resolution"] == 0.01 and tr["esr_ohm"] == 0.2, (tr["resolution"], tr["esr_ohm"])
    assert t0["c_in"] == T43_C_IN_0_H_LITERAL, (t0["c_in"], T43_C_IN_0_H_LITERAL)
    assert t0["cells"] == H001_HOLE_CELL_RECORD, (t0["cells"], H001_HOLE_CELL_RECORD)
    c_in_r = float(T43_C_IN_R_H_LITERAL)
    c_in_0 = float(T43_C_IN_0_H_LITERAL)
    print(f"{tag} C_in(r, h) {c_in_r!r} F (== tune_p2_h0.01_esr0.2.json c_in_r); C_in(0, h) "
          f"{c_in_0!r} F (== tune_p2_h0.01.json c_in); rung cells {t0['cells']}", flush=True)

    # --- templates, checked as item 421 / 436 check them ------------------------
    src = json.loads((OUT_DIR / "t39_esr.json").read_text())
    assert src["Model"]["Mesh"] == HOLE_MESH_NAME, src["Model"]["Mesh"]
    samples = src["Solver"]["Driven"]["Samples"]
    assert len(samples) == 1 and samples[0]["Freq"] == [0.064], samples
    for p in src["Boundaries"]["LumpedPort"]:
        if int(p["Index"]) == 1:
            assert p["R"] == 50.0 and p["Excitation"] not in (False, 0), p
        else:
            assert p["Excitation"] is False and set(p) >= {"R", "C"}, p
    src0 = json.loads((OUT_DIR / "t39_lossless.json").read_text())
    assert src0["Model"]["Mesh"] == HOLE_MESH_NAME, src0["Model"]["Mesh"]
    assert src0["Solver"]["Driven"]["Samples"][0]["Freq"] == [0.064], src0["Solver"]
    for p in src0["Boundaries"]["LumpedPort"]:
        if int(p["Index"]) != 1:
            assert "C" in p and "R" not in p and p["Excitation"] is False, p

    for name, f_ghz, _fh in T43_SWEEP:
        out = json.loads(json.dumps(src))
        out["Model"]["Mesh"] = H001_HOLE_MESH_NAME
        out["Solver"]["Driven"]["Samples"][0]["Freq"] = [f_ghz]
        omega = 2.0 * np.pi * f_ghz * 1.0e9
        _q, r_p, c_p = _t39_map(omega, c_in_r, T39_ESR_OHM)
        for p in out["Boundaries"]["LumpedPort"]:
            if int(p["Index"]) == 1:
                continue
            p["R"] = r_p
            p["C"] = c_p
        out["Problem"]["Output"] = f"postpro_{name}"
        path = OUT_DIR / f"{name}.json"
        path.write_text(json.dumps(out, indent=2) + "\n")
        print(f"{tag} wrote {path.relative_to(REPO)}; Mesh {out['Model']['Mesh']}; Freq {[f_ghz]} "
              f"GHz; Boundaries.LumpedPort: " + json.dumps(out["Boundaries"]["LumpedPort"]),
              flush=True)
    name0, _f0g, _f0h = T43_LOSSLESS
    out = json.loads(json.dumps(src0))
    out["Model"]["Mesh"] = H001_HOLE_MESH_NAME
    for p in out["Boundaries"]["LumpedPort"]:
        if int(p["Index"]) != 1:
            p["C"] = c_in_0
    out["Problem"]["Output"] = f"postpro_{name0}"
    path = OUT_DIR / f"{name0}.json"
    path.write_text(json.dumps(out, indent=2) + "\n")
    print(f"{tag} wrote {path.relative_to(REPO)}; Mesh {out['Model']['Mesh']}; Freq "
          f"{out['Solver']['Driven']['Samples'][0]['Freq']} GHz; Boundaries.LumpedPort: "
          + json.dumps(out["Boundaries"]["LumpedPort"]), flush=True)

    # --- (Z') / (Z0), read back; template diffs ---------------------------------
    fails = []

    def strip(d, drop_r):
        d = json.loads(json.dumps(d))
        d["Model"].pop("Mesh")
        d["Solver"]["Driven"]["Samples"][0].pop("Freq")
        d["Problem"].pop("Output")
        for p in d["Boundaries"]["LumpedPort"]:
            if int(p["Index"]) != 1:
                if drop_r:
                    p.pop("R")
                p.pop("C")
        return d

    for name, _fg, f_hz in T43_SWEEP:
        cfg = json.loads((OUT_DIR / f"{name}.json").read_text())
        if cfg["Model"]["Mesh"] != H001_HOLE_MESH_NAME:
            fails.append(f"{name}: Model.Mesh {cfg['Model']['Mesh']}")
        freq_json = cfg["Solver"]["Driven"]["Samples"][0]["Freq"]
        assert len(freq_json) == 1, freq_json
        omega_json = 2.0 * np.pi * float(freq_json[0]) * 1.0e9
        d_omega = abs(omega_json - 2.0 * np.pi * f_hz) / (2.0 * np.pi * f_hz)
        if not d_omega <= STRIP_BAND:
            fails.append(f"{name}: |d omega|/omega {d_omega:.3e} > {STRIP_BAND}")
        z_ser = series_rlc_impedance(f_hz, r_ohm=T39_ESR_OHM, c_f=c_in_r)
        for p in cfg["Boundaries"]["LumpedPort"]:
            if int(p["Index"]) == 1:
                if p["R"] != 50.0 or p["Excitation"] in (False, 0) or "C" in p:
                    fails.append(f"{name}: port 1 {p}")
                continue
            if p["Excitation"] is not False:
                fails.append(f"{name}: port {p['Index']} Excitation {p['Excitation']!r}")
            z_par = 1.0 / (1.0 / float(p["R"]) + 1j * omega_json * float(p["C"]))
            z_rel = abs(z_par - z_ser) / abs(z_ser)
            if not z_rel <= STRIP_BAND:
                fails.append(f"(Z') {name} port {p['Index']}: {z_rel:.3e} > {STRIP_BAND}")
            if int(p["Index"]) == 2:
                q, _r_p, _c_p = _t39_map(omega_json, c_in_r, T39_ESR_OHM)
                w, l = 7.0e-3, 8.0e-3
                print(f"{tag} (Z') {name}: Freq {freq_json[0]} GHz (|d omega|/omega "
                      f"{d_omega:.3e}); Q_C {q:.8f}, R_p {float(p['R']):.8f} ohm, C_p "
                      f"{float(p['C']):.10e} F (C_p/C - 1 {float(p['C']) / c_in_r - 1:.4e}); "
                      f"Z_par {_fmt(z_par)} ohm, Z_ser {_fmt(z_ser)} ohm (|Z_ser| {abs(z_ser):.4f}); "
                      f"|Z_par - Z_ser|/|Z_ser| {z_rel:.3e} vs STRIP_BAND {STRIP_BAND} (predicted "
                      f"<= 2e-16); predicted sheet Rs {float(p['R']) * w / l:.3e} ohm/sq, Cs "
                      f"{float(p['C']) * l / w:.3e} F/sq", flush=True)
        if strip(src, True) != strip(cfg, True):
            fails.append(f"{name}: differs from t39_esr.json beyond Mesh / Freq / ports 2-4 R, C / Output")
    cfg0 = json.loads((OUT_DIR / f"{name0}.json").read_text())
    if cfg0["Model"]["Mesh"] != H001_HOLE_MESH_NAME:
        fails.append(f"{name0}: Model.Mesh {cfg0['Model']['Mesh']}")
    fq0 = cfg0["Solver"]["Driven"]["Samples"][0]["Freq"]
    d0 = abs(fq0[0] * 1.0e9 - T43_LOSSLESS[2]) / T43_LOSSLESS[2]
    if not (len(fq0) == 1 and d0 <= STRIP_BAND):
        fails.append(f"{name0}: Freq {fq0}")
    for p in cfg0["Boundaries"]["LumpedPort"]:
        if int(p["Index"]) == 1:
            if p["R"] != 50.0 or p["Excitation"] in (False, 0) or "C" in p:
                fails.append(f"{name0}: port 1 {p}")
            continue
        if "R" in p or p["C"] != T43_C_IN_0_H_LITERAL or p["Excitation"] is not False:
            fails.append(f"(Z0) {name0} port {p['Index']}: {p}")
    if strip(src0, False) != strip(cfg0, False):
        fails.append(f"{name0}: differs from t39_lossless.json beyond Mesh / ports 2-4 C / Output")
    print(f"{tag} (Z0) {name0}: Freq {fq0} GHz, ports 2-4 C == C_in(0, h) {T43_C_IN_0_H_LITERAL!r} F "
          f"bit for bit, no R key; lossless Cs {T43_C_IN_0_H_LITERAL * 8.0 / 7.0:.4e} F/sq (predicted)",
          flush=True)
    # (N-map), printed: the 64.0 MHz map misused at 63.9 / 64.1 MHz
    _q, r_pc, c_pc = _t39_map(2.0 * np.pi * 64.0e6, c_in_r, T39_ESR_OHM)
    for _n, _fg, f_hz in (T43_SWEEP[0], T43_SWEEP[2]):
        z_ser = series_rlc_impedance(f_hz, r_ohm=T39_ESR_OHM, c_f=c_in_r)
        z_mis = 1.0 / (1.0 / r_pc + 1j * 2.0 * np.pi * f_hz * c_pc)
        print(f"{tag} (N-map, printed) 64.0-MHz map at {f_hz / 1e6:.1f} MHz: |Z_par(64.0-map) - "
              f"Z_ser|/|Z_ser| {abs(z_mis - z_ser) / abs(z_ser):.3e} (predicted 1.6008e-5 / "
              f"1.5983e-5; STRIP_BAND {STRIP_BAND})", flush=True)
    print(f"{tag} (Z')/(Z0) anchors on {_t43_names()}: "
          f"{'PASS' if not fails else 'FAIL ' + '; '.join(fails)}", flush=True)
    assert not fails, "; ".join(fails)


def read_ans7_h001():
    """Item 444 W7: Palace's f0 / R_in / Q on the h = 0.01 rung via parallel_rows vs T37's rows A."""
    from tests.validation.test_ans7_exact_strip import (
        T21P_ZERO_RESOLUTION,
        _ANS7,
        _z_branch,
        parallel_rows,
    )

    tag = "[OPS-61 3d W7]"
    c_in_r = float(T43_C_IN_R_H_LITERAL)
    freqs = np.array([fh for _n, _fg, fh in T43_SWEEP])

    def load(name):
        post = OUT_DIR / f"postpro_{name}"
        csv_path = post / "port-S.csv"
        if not csv_path.exists():
            print(f"{tag} {name}: not measured ({csv_path.relative_to(REPO)} absent)", flush=True)
            return None, None
        s, _hdr = _read_port_s(csv_path)
        pj = post / "palace.json"
        tot = (json.loads(pj.read_text())["ElapsedTime"]["Durations"]["Total"]
               if pj.exists() else None)
        return complex(s[0, 0]), tot

    s11, totals = [], []
    for name, _fg, _fh in T43_SWEEP:
        a, t = load(name)
        s11.append(a)
        totals.append(t)
    s11_0, tot_0 = load(T43_LOSSLESS[0])
    fails = []
    if any(x is None for x in s11):
        print(f"{tag} a point is not measured: (C1) cannot be read", flush=True)
        assert False, "point not measured"

    def rel(a, b):
        return abs(a - b) / abs(b)

    z_bare, z_in, z_in_nofeed_r = [], [], []
    for k, (name, _fg, f_hz) in enumerate(T43_SWEEP):
        omega = 2.0 * np.pi * f_hz
        zb = complex(_ANS7._z_from_s11(s11[k]))
        zbr = _z_branch(omega, c_in_r, T39_ESR_OHM)
        z_bare.append(zb)
        z_in.append(1.0 / (1.0 / zb + 1.0 / zbr))
        z_in_nofeed_r.append(1.0 / (1.0 / zb + 1j * omega * c_in_r))
        tie = rel(s11[k], T43_S11_A[k])
        tot = totals[k]
        print(f"{tag} (T6') {f_hz / 1e6:.1f} MHz ({name}): S11_Pal {_fmt(s11[k])}, ours S11_A "
              f"{_fmt(T43_S11_A[k])}; |S11_Pal - S11_A|/|S11_A| {tie:.3e} (printed; predicted <= "
              f"{T39_TIE_PREDICTED:.0e}); Z_bare {_fmt(zb)} ohm; Palace Total "
              f"{tot if tot is None else f'{tot:.3f}'} s", flush=True)
    if s11_0 is not None:
        t60 = rel(s11_0, T43_S11_A_LOSSLESS)
        print(f"{tag} (T6-0'') 64.0 MHz lossless C path: S11_Pal {_fmt(s11_0)}, ours "
              f"{_fmt(T43_S11_A_LOSSLESS)}; |S11_Pal - S11_A,lossless|/|S11_A,lossless| {t60:.3e} "
              f"(printed; predicted <= {T39_TIE_PREDICTED:.0e}, convention-free); Palace Total "
              f"{tot_0 if tot_0 is None else f'{tot_0:.3f}'} s", flush=True)
        k_pal = rel(s11[1], s11_0)
        k_ours = rel(T43_S11_A[1], T43_S11_A_LOSSLESS)
        print(f"{tag} (K-esr'') Palace esr-vs-lossless move at 64.0 MHz {k_pal:.6e}; ours "
              f"{k_ours:.6e}; |Pal/ours - 1| {abs(k_pal / k_ours - 1.0):.3e} (predicted <= 1e-6; "
              f"a dropped R reads ~1e-10 move)", flush=True)
    rows = parallel_rows(freqs, np.array(z_in))
    b = (1.0 / np.array(z_in)).imag
    print(f"{tag} Palace Im Y_in at (63.9, 64.0, 64.1) MHz = {b[0]:+.9e}, {b[1]:+.9e}, "
          f"{b[2]:+.9e} S; ours {T43_IM_Y_A[0]:+.9e}, {T43_IM_Y_A[1]:+.9e}, "
          f"{T43_IM_Y_A[2]:+.9e} S", flush=True)
    for k in range(3):
        print(f"{tag}   Im Y_in rel gap at {freqs[k] / 1e6:.1f} MHz: "
              f"{abs(b[k] - T43_IM_Y_A[k]) / abs(T43_IM_Y_A[k]):.3e} (printed)", flush=True)
    print(f"{tag} rows (Palace): {json.dumps(rows)}", flush=True)
    z_in_ours = [1.0 / (1.0 / complex(_ANS7._z_from_s11(T43_S11_A[k]))
                        + 1.0 / _z_branch(2.0 * np.pi * freqs[k], c_in_r, T39_ESR_OHM))
                 for k in range(3)]
    rows_o = parallel_rows(freqs, np.array(z_in_ours))
    print(f"{tag} ours' rows from the 12-digit S11_A prints: f0 {rows_o.get('f0_hz')!r} Hz, "
          f"R_in {rows_o.get('r_in_ohm')!r} ohm, Q {rows_o.get('q_loaded')!r} (vs logged "
          f"{T43_F0_A!r} / {T43_R_IN_A!r} / {T43_Q_A!r})", flush=True)

    # --- (C1) one crossing ------------------------------------------------------
    up, down = rows["crossings_up_hz"], rows["crossings_down_hz"]
    c1 = len(up) == 1 and not down and 63.9e6 < up[0] < 64.1e6
    print(f"{tag} (C1) up-crossings {up}, down-crossings {down}: "
          f"{'PASS' if c1 else 'FAIL'} (ASSERTED: exactly one - -> +, none + -> -, inside "
          f"(63.9, 64.1) MHz)", flush=True)
    if not c1:
        fails.append("(C1)")
    else:
        f0, r_in, q = rows["f0_hz"], rows["r_in_ohm"], rows["q_loaded"]
        df0 = f0 / T43_F0_A - 1.0
        dq = q / T43_Q_A - 1.0
        drin = r_in / T43_R_IN_A - 1.0
        f0_ok = abs(df0) <= T21P_ZERO_RESOLUTION
        q_ok = abs(dq) <= T40_Q_TARGET
        print(f"{tag} (f0) f0,Pal {f0!r} Hz vs f0,A {T43_F0_A!r} Hz: f0,Pal/f0,A - 1 = "
              f"{df0:+.3e} (ASSERTED |.| <= T21P_ZERO_RESOLUTION {T21P_ZERO_RESOLUTION:.3e}: "
              f"{'PASS' if f0_ok else 'FAIL'}; predicted <= 1e-6); f0,Pal/64 MHz - 1 = "
              f"{f0 / 64.0e6 - 1.0:+.6e} (ours {T43_F0_A / 64.0e6 - 1.0:+.6e})", flush=True)
        print(f"{tag} (Q) Q_Pal {q!r} vs Q_A {T43_Q_A!r}: Q_Pal/Q_A - 1 = {dq:+.3e} (ASSERTED "
              f"|.| <= T40_Q_TARGET {T40_Q_TARGET} (§10 T40's 5 %): {'PASS' if q_ok else 'FAIL'}; "
              f"predicted <= 1e-6)", flush=True)
        print(f"{tag} R_in,Pal {r_in!r} ohm vs R_in,A {T43_R_IN_A!r} ohm: R_in,Pal/R_in,A - 1 = "
              f"{drin:+.3e} (printed; predicted <= 1e-6)", flush=True)
        if not f0_ok:
            fails.append("(f0)")
        if not q_ok:
            fails.append("(Q)")
        # --- (C-h) Palace's own x1 -> 0.01 move beside ours' ------------------
        pf0, pq, prin = T43_PAL_X1_ROWS
        mv_pal = (f0 / pf0 - 1.0, q / pq - 1.0, r_in / prin - 1.0)
        mv_ours = (T43_F0_A / T40_F0_A - 1.0, T43_Q_A / T40_Q_A - 1.0, T43_R_IN_A / T40_R_IN_A - 1.0)
        print(f"{tag} (C-h) x1 -> 0.01 move, Palace (vs item 421 W4 rows) / ours (T26 -> T37 rows): "
              f"f0 {mv_pal[0]:+.6e} / {mv_ours[0]:+.6e}; Q {mv_pal[1]:+.6e} / {mv_ours[1]:+.6e} "
              f"(ours logged {T43_Q_H_MOVE_OURS:+.6e}); R_in {mv_pal[2]:+.6e} / {mv_ours[2]:+.6e} "
              f"(printed; predicted equal to the printed digit)", flush=True)
        # --- (C2-h) mesh-mismatch separation on Q -----------------------------
        num_q = abs(T40_Q_A - q) / q
        den_q = abs(T43_Q_A - q) / q
        print(f"{tag} (C2-h) Q: |Q_A(x1) - Q_Pal(h)| / |Q_A(h) - Q_Pal(h)| = "
              f"{num_q:.3e} / {den_q:.3e} = "
              f"{(num_q / den_q) if den_q > 0 else float('inf'):.3e} (control; predicted >= 2000; "
              f"bar 25)", flush=True)

    # (C2-h) on S11 at 64.0 MHz
    num = abs(T39_S11_A_ESR - s11[1])
    den = abs(T43_S11_A[1] - s11[1])
    print(f"{tag} (C2-h) S11 64.0 MHz: |S11_A(x1) - S11_Pal(h)| / |S11_A(h) - S11_Pal(h)| = "
          f"{num / abs(s11[1]):.3e} / {den / abs(s11[1]):.3e} = "
          f"{(num / den) if den > 0 else float('inf'):.3e} (control; predicted >= 810, central ~3e6, "
          f"ceiling ~8e8; bar 25)", flush=True)

    # --- negative controls (rule (e); predicted, printed) -----------------------
    rows_nf = parallel_rows(freqs, np.array(z_in_nofeed_r))
    if "q_loaded" in rows_nf:
        qn = rows_nf["q_loaded"]
        in_band = 0.20 <= qn / T43_Q_A - 1.0 <= 0.80
        print(f"{tag} (N-feed) feed's r omitted: Q' {qn!r}, Q'/Q_A - 1 = {qn / T43_Q_A - 1.0:+.3e} "
              f"(predicted in [+20 %, +80 %]: {'met' if in_band else 'NOT met'}; ceiling Q_0(h) "
              f"{T43_Q_0_H:.2f} -> {T43_Q_0_H / T43_Q_A - 1.0:+.0%}); f0' {rows_nf['f0_hz']!r} Hz "
              f"(printed)", flush=True)
    else:
        print(f"{tag} (N-feed) feed's r omitted: no - -> + crossing ({json.dumps(rows_nf)}) "
              f"(predicted a crossing with Q' in [+20 %, +80 %]: NOT met)", flush=True)
    rows_nb = parallel_rows(freqs, np.array(z_bare))
    bb = (1.0 / np.array(z_bare)).imag
    if "f0_hz" in rows_nb:
        print(f"{tag} (N-nofeed) feed branch dropped: crossing at {rows_nb['f0_hz']!r} Hz "
              f"(predicted no crossing: NOT met); Im Y_bare {bb[0]:+.4e}, {bb[1]:+.4e}, "
              f"{bb[2]:+.4e} S", flush=True)
    else:
        print(f"{tag} (N-nofeed) feed branch dropped: no crossing / f0 undefined (predicted: met); "
              f"Im Y_bare {bb[0]:+.4e}, {bb[1]:+.4e}, {bb[2]:+.4e} S; rows {json.dumps(rows_nb)}",
              flush=True)
    print(f"{tag} anchors (C1)/(f0)/(Q): {'PASS' if not fails else 'FAIL ' + ', '.join(fails)}; "
          f"ties, R_in and controls printed, none asserted (rulings 2026-10-13 (1) / (4))", flush=True)
    assert not fails, ", ".join(fails)


# --- OPS-61 step 3e (§9 item 448, §10 T44, provisional): item 444 at 127.74 MHz -
# Palace's twin + driven read of T38's r = 0.2 ohm tune at 127.74 MHz on the SAME
# h = 0.01 exact-rectangle hole mesh item 444 exported (never re-exported here),
# against T38's Route A rows (20261013T125603Z_ANS-7-item419-W1-esr128-h001.log,
# the window that wrote t127.74/tune_p2_h0.01_esr0.2.json).  C_in(r, h) is the
# 127.74 MHz rung's own value -- never a x1 (T35_*) or a 64 MHz (T43_*) literal.
# Additive: ans7_h001_configs / read_ans7_h001 / read_ans7_128 are not edited.
T44_SWEEP = (("t44_esr_f127.64", 0.12764, 127.64e6), ("t44_esr_f127.74", 0.12774, 127.74e6),
             ("t44_esr_f127.84", 0.12784, 127.84e6))
T44_LOSSLESS = ("t44_lossless_f127.74", 0.12774, 127.74e6)
# ...item419-W1-esr128-h001.log:476 (== t127.74/tune_p2_h0.01_esr0.2.json c_in_r / c_in_0)
T44_C_IN_R_H_LITERAL = 1.5288749690061608e-11
T44_C_IN_0_H_LITERAL = 1.5288541964455146e-11
# Ours' S11_A at (r, C_in(r, h)), 12-digit prints, ...item419-W1-esr128-h001.log:2204-2206:
T44_S11_A = (4.486361641616e-01 + 8.878179861686e-01j,
             4.506307628802e-01 + 8.867812811584e-01j,
             4.526280544402e-01 + 8.857371566245e-01j)
# T38 rows A, ...item419-W1-esr128-h001.log:2207:
T44_F0_A = 127739705.26328915
T44_R_IN_A = 9698.113320373346
T44_Q_A = 252.30120548497052
T44_IM_Y_A = (-4.061198730e-05, +1.200522737e-07, +4.089005599e-05)
T44_Q_H_MOVE_OURS = -8.519904e-03  # ...item419-W1-esr128-h001.log:2211 (vs T35_Q_R_X1)
T44_R_IN_H_MOVE_OURS = -9.456959e-03  # ...item419-W1-esr128-h001.log:2212 (vs T35_R_IN_R_X1)
T44_Q_0_H = 606.6969050131794  # ...item419-W1-esr128-h001.log:2210 (printed ceiling only)
# Lossless C path on the rung at C_in(0, h), 127.74 MHz: 20261011T111940Z_ANS-7-item365b-D2-p127.74.log:2288
T44_S11_A_LOSSLESS = 4.516287177353e-01 + 8.887181496060e-01j
# Palace's x1 rows at 127.74 MHz (f0, Q, R_in): 20261013T215115Z_OPS-61-item436-W6-read.log:48
T44_PAL_X1_ROWS = (127739705.97524849, 254.46926048070927, 9790.703634300547)


def _t44_names():
    return [n for n, _g, _h in T44_SWEEP] + [T44_LOSSLESS[0]]


def ans7_h001_128_configs():
    """Item 448 W1: write t44_esr_f{127.64,127.74,127.84} / t44_lossless_f127.74; (E1') / (Z') / (Z0).

    A copy of ``ans7_h001_configs`` (commit 2877ca28, lines 3795-3949) with the
    T44 sets and literals, the t127.74/ tune JSONs, an (E1') re-read of the
    reused rung mesh before any write, and the lossless ``Freq`` set to 127.74 MHz.
    """
    from fem_em_solver.ports.lumped import series_rlc_impedance
    from tests.validation.test_ops61_exact_strip_ports import STRIP_BAND

    tag = "[OPS-61 3e W1]"
    try:
        from mpi4py import MPI

        assert MPI.COMM_WORLD.size == 1, f"run at -n 1, not {MPI.COMM_WORLD.size}"
    except ImportError:
        pass
    for name in _t44_names():
        path = OUT_DIR / f"{name}.json"
        assert not path.exists(), f"{path.relative_to(REPO)} exists: never overwrite"

    # --- (E1') the reused rung mesh's $Elements type-4 count, before any write ---
    mesh_path = OUT_DIR / H001_HOLE_MESH_NAME
    assert mesh_path.exists(), "item 444's rung mesh is absent (never re-export here)"
    n_tet_file = 0
    with open(mesh_path) as fh:
        inside = False
        for line in fh:
            if line.startswith("$Elements"):
                inside = True
                next(fh)
                continue
            if line.startswith("$EndElements"):
                break
            if inside and line.split(None, 2)[1] == "4":
                n_tet_file += 1
    e1 = n_tet_file == H001_HOLE_CELL_RECORD
    print(f"{tag} (E1') {H001_HOLE_MESH_NAME}: $Elements type-4 count {n_tet_file} vs "
          f"H001_HOLE_CELL_RECORD {H001_HOLE_CELL_RECORD}: {'PASS' if e1 else 'FAIL'} (ASSERTED)",
          flush=True)
    assert e1, (n_tet_file, H001_HOLE_CELL_RECORD)

    t128 = REPO / "logs" / "ans7_exact" / "t127.74"
    tr = json.loads((t128 / "tune_p2_h0.01_esr0.2.json").read_text())
    t0 = json.loads((t128 / "tune_p2_h0.01.json").read_text())
    assert tr["c_in_r"] == T44_C_IN_R_H_LITERAL, (tr["c_in_r"], T44_C_IN_R_H_LITERAL)
    assert tr["c_in_0"] == T44_C_IN_0_H_LITERAL, (tr["c_in_0"], T44_C_IN_0_H_LITERAL)
    assert tr["resolution"] == 0.01 and tr["esr_ohm"] == 0.2, (tr["resolution"], tr["esr_ohm"])
    assert tr["f_hz"] == 127740000.0, tr["f_hz"]
    assert t0["c_in"] == T44_C_IN_0_H_LITERAL, (t0["c_in"], T44_C_IN_0_H_LITERAL)
    assert t0["cells"] == H001_HOLE_CELL_RECORD, (t0["cells"], H001_HOLE_CELL_RECORD)
    assert t0["f_hz"] == 127740000.0, t0["f_hz"]
    c_in_r = float(T44_C_IN_R_H_LITERAL)
    c_in_0 = float(T44_C_IN_0_H_LITERAL)
    print(f"{tag} C_in(r, h) {c_in_r!r} F (== t127.74/tune_p2_h0.01_esr0.2.json c_in_r, f_hz "
          f"{tr['f_hz']!r}); C_in(0, h) {c_in_0!r} F (== t127.74/tune_p2_h0.01.json c_in, f_hz "
          f"{t0['f_hz']!r}); rung cells {t0['cells']}", flush=True)

    # --- templates, checked as item 444 checks them -----------------------------
    src = json.loads((OUT_DIR / "t39_esr.json").read_text())
    assert src["Model"]["Mesh"] == HOLE_MESH_NAME, src["Model"]["Mesh"]
    samples = src["Solver"]["Driven"]["Samples"]
    assert len(samples) == 1 and samples[0]["Freq"] == [0.064], samples
    for p in src["Boundaries"]["LumpedPort"]:
        if int(p["Index"]) == 1:
            assert p["R"] == 50.0 and p["Excitation"] not in (False, 0), p
        else:
            assert p["Excitation"] is False and set(p) >= {"R", "C"}, p
    src0 = json.loads((OUT_DIR / "t39_lossless.json").read_text())
    assert src0["Model"]["Mesh"] == HOLE_MESH_NAME, src0["Model"]["Mesh"]
    assert src0["Solver"]["Driven"]["Samples"][0]["Freq"] == [0.064], src0["Solver"]
    for p in src0["Boundaries"]["LumpedPort"]:
        if int(p["Index"]) != 1:
            assert "C" in p and "R" not in p and p["Excitation"] is False, p

    for name, f_ghz, _fh in T44_SWEEP:
        out = json.loads(json.dumps(src))
        out["Model"]["Mesh"] = H001_HOLE_MESH_NAME
        out["Solver"]["Driven"]["Samples"][0]["Freq"] = [f_ghz]
        omega = 2.0 * np.pi * f_ghz * 1.0e9
        _q, r_p, c_p = _t39_map(omega, c_in_r, T39_ESR_OHM)
        for p in out["Boundaries"]["LumpedPort"]:
            if int(p["Index"]) == 1:
                continue
            p["R"] = r_p
            p["C"] = c_p
        out["Problem"]["Output"] = f"postpro_{name}"
        path = OUT_DIR / f"{name}.json"
        path.write_text(json.dumps(out, indent=2) + "\n")
        print(f"{tag} wrote {path.relative_to(REPO)}; Mesh {out['Model']['Mesh']}; Freq {[f_ghz]} "
              f"GHz; Boundaries.LumpedPort: " + json.dumps(out["Boundaries"]["LumpedPort"]),
              flush=True)
    name0, f0_ghz, f0_hz = T44_LOSSLESS
    out = json.loads(json.dumps(src0))
    out["Model"]["Mesh"] = H001_HOLE_MESH_NAME
    out["Solver"]["Driven"]["Samples"][0]["Freq"] = [f0_ghz]
    for p in out["Boundaries"]["LumpedPort"]:
        if int(p["Index"]) != 1:
            p["C"] = c_in_0
    out["Problem"]["Output"] = f"postpro_{name0}"
    path = OUT_DIR / f"{name0}.json"
    path.write_text(json.dumps(out, indent=2) + "\n")
    print(f"{tag} wrote {path.relative_to(REPO)}; Mesh {out['Model']['Mesh']}; Freq "
          f"{out['Solver']['Driven']['Samples'][0]['Freq']} GHz; Boundaries.LumpedPort: "
          + json.dumps(out["Boundaries"]["LumpedPort"]), flush=True)

    # --- (Z') / (Z0), read back; template diffs ---------------------------------
    fails = []

    def strip(d, drop_r):
        d = json.loads(json.dumps(d))
        d["Model"].pop("Mesh")
        d["Solver"]["Driven"]["Samples"][0].pop("Freq")
        d["Problem"].pop("Output")
        for p in d["Boundaries"]["LumpedPort"]:
            if int(p["Index"]) != 1:
                if drop_r:
                    p.pop("R")
                p.pop("C")
        return d

    for name, _fg, f_hz in T44_SWEEP:
        cfg = json.loads((OUT_DIR / f"{name}.json").read_text())
        if cfg["Model"]["Mesh"] != H001_HOLE_MESH_NAME:
            fails.append(f"{name}: Model.Mesh {cfg['Model']['Mesh']}")
        freq_json = cfg["Solver"]["Driven"]["Samples"][0]["Freq"]
        assert len(freq_json) == 1, freq_json
        omega_json = 2.0 * np.pi * float(freq_json[0]) * 1.0e9
        d_omega = abs(omega_json - 2.0 * np.pi * f_hz) / (2.0 * np.pi * f_hz)
        if not d_omega <= STRIP_BAND:
            fails.append(f"{name}: |d omega|/omega {d_omega:.3e} > {STRIP_BAND}")
        z_ser = series_rlc_impedance(f_hz, r_ohm=T39_ESR_OHM, c_f=c_in_r)
        for p in cfg["Boundaries"]["LumpedPort"]:
            if int(p["Index"]) == 1:
                if p["R"] != 50.0 or p["Excitation"] in (False, 0) or "C" in p:
                    fails.append(f"{name}: port 1 {p}")
                continue
            if p["Excitation"] is not False:
                fails.append(f"{name}: port {p['Index']} Excitation {p['Excitation']!r}")
            z_par = 1.0 / (1.0 / float(p["R"]) + 1j * omega_json * float(p["C"]))
            z_rel = abs(z_par - z_ser) / abs(z_ser)
            if not z_rel <= STRIP_BAND:
                fails.append(f"(Z') {name} port {p['Index']}: {z_rel:.3e} > {STRIP_BAND}")
            if int(p["Index"]) == 2:
                q, _r_p, _c_p = _t39_map(omega_json, c_in_r, T39_ESR_OHM)
                w, l = 7.0e-3, 8.0e-3
                print(f"{tag} (Z') {name}: Freq {freq_json[0]} GHz (|d omega|/omega "
                      f"{d_omega:.3e}); Q_C {q:.8f}, R_p {float(p['R']):.8f} ohm, C_p "
                      f"{float(p['C']):.10e} F (C_p/C - 1 {float(p['C']) / c_in_r - 1:.4e}); "
                      f"Z_par {_fmt(z_par)} ohm, Z_ser {_fmt(z_ser)} ohm (|Z_ser| {abs(z_ser):.4f}); "
                      f"|Z_par - Z_ser|/|Z_ser| {z_rel:.3e} vs STRIP_BAND {STRIP_BAND} (predicted "
                      f"<= 2e-16); predicted sheet Rs {float(p['R']) * w / l:.4e} ohm/sq, Cs "
                      f"{float(p['C']) * l / w:.4e} F/sq", flush=True)
        if strip(src, True) != strip(cfg, True):
            fails.append(f"{name}: differs from t39_esr.json beyond Mesh / Freq / ports 2-4 R, C / Output")
    cfg0 = json.loads((OUT_DIR / f"{name0}.json").read_text())
    if cfg0["Model"]["Mesh"] != H001_HOLE_MESH_NAME:
        fails.append(f"{name0}: Model.Mesh {cfg0['Model']['Mesh']}")
    fq0 = cfg0["Solver"]["Driven"]["Samples"][0]["Freq"]
    if fq0 != [0.12774]:
        fails.append(f"(Z0) {name0}: Freq {fq0} != [0.12774]")
    d0 = abs(fq0[0] * 1.0e9 - f0_hz) / f0_hz
    if not (len(fq0) == 1 and d0 <= STRIP_BAND):
        fails.append(f"{name0}: Freq {fq0}")
    for p in cfg0["Boundaries"]["LumpedPort"]:
        if int(p["Index"]) == 1:
            if p["R"] != 50.0 or p["Excitation"] in (False, 0) or "C" in p:
                fails.append(f"{name0}: port 1 {p}")
            continue
        if "R" in p or p["C"] != T44_C_IN_0_H_LITERAL or p["Excitation"] is not False:
            fails.append(f"(Z0) {name0} port {p['Index']}: {p}")
    if strip(src0, False) != strip(cfg0, False):
        fails.append(f"{name0}: differs from t39_lossless.json beyond Mesh / Freq / ports 2-4 C / Output")
    print(f"{tag} (Z0) {name0}: Freq {fq0} GHz (== [0.12774]), ports 2-4 C == C_in(0, h) "
          f"{T44_C_IN_0_H_LITERAL!r} F bit for bit, no R key; lossless Cs "
          f"{T44_C_IN_0_H_LITERAL * 8.0 / 7.0:.4e} F/sq (predicted)", flush=True)
    # (N-map), printed: the 127.74 MHz map misused at 127.64 / 127.84 MHz
    _q, r_pc, c_pc = _t39_map(2.0 * np.pi * 127.74e6, c_in_r, T39_ESR_OHM)
    for _n, _fg, f_hz in (T44_SWEEP[0], T44_SWEEP[2]):
        z_ser = series_rlc_impedance(f_hz, r_ohm=T39_ESR_OHM, c_f=c_in_r)
        z_mis = 1.0 / (1.0 / r_pc + 1j * 2.0 * np.pi * f_hz * c_pc)
        print(f"{tag} (N-map, printed) 127.74-MHz map at {f_hz / 1e6:.2f} MHz: |Z_par(127.74-map) - "
              f"Z_ser|/|Z_ser| {abs(z_mis - z_ser) / abs(z_ser):.3e} (predicted 3.8440e-6 / "
              f"3.8410e-6; STRIP_BAND {STRIP_BAND})", flush=True)
    print(f"{tag} (E1')/(Z')/(Z0) anchors on {_t44_names()}: "
          f"{'PASS' if not fails else 'FAIL ' + '; '.join(fails)}", flush=True)
    assert not fails, "; ".join(fails)


def read_ans7_h001_128():
    """Item 448 W6: Palace's f0 / R_in / Q on the h = 0.01 rung at 127.74 MHz vs T38's rows A.

    A copy of ``read_ans7_h001`` (commit 2877ca28, lines 3952-4109) with every
    64 MHz comparand replaced by its 127.74 MHz counterpart: T44_* (T38's rows),
    T44_PAL_X1_ROWS (item 436 W6), T42_* (ours x1, 127.74 MHz).
    """
    from tests.validation.test_ans7_exact_strip import (
        T21P_ZERO_RESOLUTION,
        T22P_ZERO_RESOLUTION,
        _ANS7,
        _z_branch,
        parallel_rows,
    )

    tag = "[OPS-61 3e W6]"
    c_in_r = float(T44_C_IN_R_H_LITERAL)
    freqs = np.array([fh for _n, _fg, fh in T44_SWEEP])

    def load(name):
        post = OUT_DIR / f"postpro_{name}"
        csv_path = post / "port-S.csv"
        if not csv_path.exists():
            print(f"{tag} {name}: not measured ({csv_path.relative_to(REPO)} absent)", flush=True)
            return None, None
        s, _hdr = _read_port_s(csv_path)
        pj = post / "palace.json"
        tot = (json.loads(pj.read_text())["ElapsedTime"]["Durations"]["Total"]
               if pj.exists() else None)
        return complex(s[0, 0]), tot

    s11, totals = [], []
    for name, _fg, _fh in T44_SWEEP:
        a, t = load(name)
        s11.append(a)
        totals.append(t)
    s11_0, tot_0 = load(T44_LOSSLESS[0])
    fails = []
    if any(x is None for x in s11):
        print(f"{tag} a point is not measured: (C1) cannot be read", flush=True)
        assert False, "point not measured"

    def rel(a, b):
        return abs(a - b) / abs(b)

    z_bare, z_in, z_in_nofeed_r = [], [], []
    for k, (name, _fg, f_hz) in enumerate(T44_SWEEP):
        omega = 2.0 * np.pi * f_hz
        zb = complex(_ANS7._z_from_s11(s11[k]))
        zbr = _z_branch(omega, c_in_r, T39_ESR_OHM)
        z_bare.append(zb)
        z_in.append(1.0 / (1.0 / zb + 1.0 / zbr))
        z_in_nofeed_r.append(1.0 / (1.0 / zb + 1j * omega * c_in_r))
        tie = rel(s11[k], T44_S11_A[k])
        tot = totals[k]
        print(f"{tag} (T6') {f_hz / 1e6:.2f} MHz ({name}): S11_Pal {_fmt(s11[k])}, ours S11_A "
              f"{_fmt(T44_S11_A[k])}; |S11_Pal - S11_A|/|S11_A| {tie:.3e} (printed; predicted <= "
              f"{T39_TIE_PREDICTED:.0e}); Z_bare {_fmt(zb)} ohm; Palace Total "
              f"{tot if tot is None else f'{tot:.3f}'} s", flush=True)
    if s11_0 is not None:
        t60 = rel(s11_0, T44_S11_A_LOSSLESS)
        print(f"{tag} (T6-0''') 127.74 MHz lossless C path: S11_Pal {_fmt(s11_0)}, ours "
              f"{_fmt(T44_S11_A_LOSSLESS)}; |S11_Pal - S11_A,lossless|/|S11_A,lossless| {t60:.3e} "
              f"(printed; predicted <= {T39_TIE_PREDICTED:.0e}, convention-free); Palace Total "
              f"{tot_0 if tot_0 is None else f'{tot_0:.3f}'} s", flush=True)
        k_pal = rel(s11[1], s11_0)
        k_ours = rel(T44_S11_A[1], T44_S11_A_LOSSLESS)
        print(f"{tag} (K-esr''') Palace esr-vs-lossless move at 127.74 MHz {k_pal:.6e}; ours "
              f"{k_ours:.6e}; |Pal/ours - 1| {abs(k_pal / k_ours - 1.0):.3e} (predicted <= 1e-6; "
              f"a dropped R reads ~1e-10 move)", flush=True)
    rows = parallel_rows(freqs, np.array(z_in))
    b = (1.0 / np.array(z_in)).imag
    print(f"{tag} Palace Im Y_in at (127.64, 127.74, 127.84) MHz = {b[0]:+.9e}, {b[1]:+.9e}, "
          f"{b[2]:+.9e} S; ours {T44_IM_Y_A[0]:+.9e}, {T44_IM_Y_A[1]:+.9e}, "
          f"{T44_IM_Y_A[2]:+.9e} S", flush=True)
    for k in range(3):
        print(f"{tag}   Im Y_in rel gap at {freqs[k] / 1e6:.2f} MHz: "
              f"{abs(b[k] - T44_IM_Y_A[k]) / abs(T44_IM_Y_A[k]):.3e} (printed)", flush=True)
    print(f"{tag} rows (Palace): {json.dumps(rows)}", flush=True)
    z_in_ours = [1.0 / (1.0 / complex(_ANS7._z_from_s11(T44_S11_A[k]))
                        + 1.0 / _z_branch(2.0 * np.pi * freqs[k], c_in_r, T39_ESR_OHM))
                 for k in range(3)]
    rows_o = parallel_rows(freqs, np.array(z_in_ours))
    print(f"{tag} ours' rows from the 12-digit S11_A prints: f0 {rows_o.get('f0_hz')!r} Hz, "
          f"R_in {rows_o.get('r_in_ohm')!r} ohm, Q {rows_o.get('q_loaded')!r} (vs logged "
          f"{T44_F0_A!r} / {T44_R_IN_A!r} / {T44_Q_A!r})", flush=True)

    # --- (C1) one crossing ------------------------------------------------------
    up, down = rows["crossings_up_hz"], rows["crossings_down_hz"]
    c1 = len(up) == 1 and not down and 127.64e6 < up[0] < 127.84e6
    print(f"{tag} (C1) up-crossings {up}, down-crossings {down}: "
          f"{'PASS' if c1 else 'FAIL'} (ASSERTED: exactly one - -> +, none + -> -, inside "
          f"(127.64, 127.84) MHz)", flush=True)
    if not c1:
        fails.append("(C1)")
    else:
        f0, r_in, q = rows["f0_hz"], rows["r_in_ohm"], rows["q_loaded"]
        df0 = f0 / T44_F0_A - 1.0
        dq = q / T44_Q_A - 1.0
        drin = r_in / T44_R_IN_A - 1.0
        f0_ok = abs(df0) <= T21P_ZERO_RESOLUTION
        q_ok = abs(dq) <= T40_Q_TARGET
        print(f"{tag} (f0) f0,Pal {f0!r} Hz vs f0,A {T44_F0_A!r} Hz: f0,Pal/f0,A - 1 = "
              f"{df0:+.3e} (ASSERTED |.| <= T21P_ZERO_RESOLUTION {T21P_ZERO_RESOLUTION:.3e}: "
              f"{'PASS' if f0_ok else 'FAIL'}; printed beside T22P_ZERO_RESOLUTION "
              f"{T22P_ZERO_RESOLUTION:.3e}: {'within' if abs(df0) <= T22P_ZERO_RESOLUTION else 'outside'}; "
              f"predicted <= 1e-6); f0,Pal/127.74 MHz - 1 = {f0 / 127.74e6 - 1.0:+.6e} "
              f"(ours {T44_F0_A / 127.74e6 - 1.0:+.6e})", flush=True)
        print(f"{tag} (Q) Q_Pal {q!r} vs Q_A {T44_Q_A!r}: Q_Pal/Q_A - 1 = {dq:+.3e} (ASSERTED "
              f"|.| <= T40_Q_TARGET {T40_Q_TARGET} (§10 T40's 5 %): {'PASS' if q_ok else 'FAIL'}; "
              f"predicted <= 1e-6)", flush=True)
        print(f"{tag} R_in,Pal {r_in!r} ohm vs R_in,A {T44_R_IN_A!r} ohm: R_in,Pal/R_in,A - 1 = "
              f"{drin:+.3e} (printed; predicted <= 1e-6)", flush=True)
        if not f0_ok:
            fails.append("(f0)")
        if not q_ok:
            fails.append("(Q)")
        # --- (C-h) Palace's own x1 -> 0.01 move at 127.74 MHz beside ours' ----
        pf0, pq, prin = T44_PAL_X1_ROWS
        mv_pal = (f0 / pf0 - 1.0, q / pq - 1.0, r_in / prin - 1.0)
        mv_ours = (T44_F0_A / T42_F0_A - 1.0, T44_Q_A / T42_Q_A - 1.0, T44_R_IN_A / T42_R_IN_A - 1.0)
        print(f"{tag} (C-h) x1 -> 0.01 move at 127.74 MHz, Palace (vs item 436 W6 rows) / ours "
              f"(T42_* -> T38 rows): f0 {mv_pal[0]:+.6e} / {mv_ours[0]:+.6e}; Q {mv_pal[1]:+.6e} / "
              f"{mv_ours[1]:+.6e} (ours logged {T44_Q_H_MOVE_OURS:+.6e} vs T35_Q_R_X1); R_in "
              f"{mv_pal[2]:+.6e} / {mv_ours[2]:+.6e} (ours logged {T44_R_IN_H_MOVE_OURS:+.6e}); "
              f"|dQ-move| {abs(mv_pal[1] - mv_ours[1]):.3e} (printed; predicted equal to ~5 "
              f"significant digits, |dQ-move| <~ 1e-8)", flush=True)
        # --- (C2-h) mesh-mismatch separation on Q -----------------------------
        num_q = abs(T42_Q_A - q) / q
        den_q = abs(T44_Q_A - q) / q
        print(f"{tag} (C2-h) Q: |Q_A(x1) - Q_Pal(h)| / |Q_A(h) - Q_Pal(h)| = "
              f"{num_q:.3e} / {den_q:.3e} = "
              f"{(num_q / den_q) if den_q > 0 else float('inf'):.3e} (control; predicted >= 8000)",
              flush=True)

    # (C2-h) on S11 at 127.74 MHz, numerator from ours' x1 127.74 MHz S11 (T42_S11_A[1])
    num = abs(T42_S11_A[1] - s11[1])
    den = abs(T44_S11_A[1] - s11[1])
    print(f"{tag} (C2-h) S11 127.74 MHz: |S11_A(x1) - S11_Pal(h)| / |S11_A(h) - S11_Pal(h)| = "
          f"{num / abs(s11[1]):.3e} / {den / abs(s11[1]):.3e} = "
          f"{(num / den) if den > 0 else float('inf'):.3e} (control; predicted >= 800, central ~3e6)",
          flush=True)

    # --- negative controls (rule (e); predicted, printed) -----------------------
    rows_nf = parallel_rows(freqs, np.array(z_in_nofeed_r))
    if "q_loaded" in rows_nf:
        qn = rows_nf["q_loaded"]
        in_band = 0.20 <= qn / T44_Q_A - 1.0 <= 0.70
        print(f"{tag} (N-feed) feed's r omitted: Q' {qn!r}, Q'/Q_A - 1 = {qn / T44_Q_A - 1.0:+.3e} "
              f"(predicted in [+20 %, +70 %], central +41.3 %: {'met' if in_band else 'NOT met'}; "
              f"ceiling Q_0(h) {T44_Q_0_H:.2f} -> {T44_Q_0_H / T44_Q_A - 1.0:+.0%}); f0' "
              f"{rows_nf['f0_hz']!r} Hz (printed)", flush=True)
    else:
        print(f"{tag} (N-feed) feed's r omitted: no - -> + crossing ({json.dumps(rows_nf)}) "
              f"(predicted a crossing with Q' in [+20 %, +70 %]: NOT met)", flush=True)
    rows_nb = parallel_rows(freqs, np.array(z_bare))
    bb = (1.0 / np.array(z_bare)).imag
    if "f0_hz" in rows_nb:
        print(f"{tag} (N-nofeed) feed branch dropped: crossing at {rows_nb['f0_hz']!r} Hz "
              f"(predicted no crossing: NOT met); Im Y_bare {bb[0]:+.4e}, {bb[1]:+.4e}, "
              f"{bb[2]:+.4e} S", flush=True)
    else:
        print(f"{tag} (N-nofeed) feed branch dropped: no crossing / f0 undefined (predicted: met); "
              f"Im Y_bare {bb[0]:+.4e}, {bb[1]:+.4e}, {bb[2]:+.4e} S; rows {json.dumps(rows_nb)}",
              flush=True)
    print(f"{tag} anchors (C1)/(f0)/(Q): {'PASS' if not fails else 'FAIL ' + ', '.join(fails)}; "
          f"ties, R_in and controls printed, none asserted (rulings 2026-10-13 (1) / (4))", flush=True)
    assert not fails, ", ".join(fails)


# ---------------------------------------------------------------------------
# Item 454 (OPS-61 step 2c, section 10 F37): the h = 0.0075 exact-rectangle
# hole mesh (ANS-6 step 5b's rung, xl entry 31) and Palace order 2 on it.
# Comparand: docs/testing/logs/20261012T070009Z_ANS-6-step5b.log (not re-solved).
H0075_HOLE_MESH_NAME = "ops61_ans6_h0075_hole_strip.msh"
# entry 31 log :2093 "hole flag on: cells=226904", :3101, :4587 (no tracked constant holds it)
H0075_HOLE_CELL_RECORD = 226904
# entry 31 log :3102 / :4587 degree-2 unknowns (a predicted Palace reading, not an anchor)
H0075_HOLE_ND_P2 = 1538030
H0075_HOLE_SETS = {col: {f"h05x_h0075_{col}_o2_e{k}": k for k in (1, 2, 3, 4)}
                   for col in ("cu", "pec")}
H0075_HOLE_TEMPLATES = {"cu": "h05x_cu_o2.json", "pec": "h05x_pec_o2.json"}
H0075_HOLE_OURS_LOG = REPO / "docs" / "testing" / "logs" / "20261012T070009Z_ANS-6-step5b.log"
H0075_HOLE_OURS = {
    "cu": {
        "header": "sigma = 5.8e+07 S/m: Z_s =",  # :4560
        "label": "S_",  # rows :4561-4564
        # :4589 "cu class means (10 digits)"
        "means": {"self": -1.353440825e-02 + 6.706648129e-01j,
                  "adjacent": +3.910222169e-01 - 2.254546627e-01j,
                  "opposite": +2.282814942e-01 - 2.986006208e-01j},
        "sigma_max": 0.999904939998,  # :4565
        "h_move": {"self": 1.206e-3, "adjacent": 5.06e-4, "opposite": 1.091e-3},  # :4593
        "pal_x1": "postpro_h05x_cu_o2",
    },
    "pec": {
        "header": "[TH-14 step2] f = 6.400000e+07 Hz: PEC sweep",  # :4555
        "label": "S_PEC_",  # rows :4556-4559
        # :4590 "pec class means (10 digits)"
        "means": {"self": -1.390952768e-02 + 6.708375558e-01j,
                  "adjacent": +3.911462934e-01 - 2.254772062e-01j,
                  "opposite": +2.284087254e-01 - 2.987280517e-01j},
        "sigma_max": None,  # not printed at 10 digits for PEC
        "h_move": {"self": 1.193e-3, "adjacent": 5.01e-4, "opposite": 1.081e-3},  # :4594
        "pal_x1": "postpro_h05x_pec_o2",
    },
}
# :4592 cu vs pec class-mean relative gap (ours, this rung)
H0075_HOLE_CU_PEC_GAP_OURS = {"self": 6.15e-4, "adjacent": 2.79e-4, "opposite": 4.79e-4}


def export_hole_h0075():
    """Item 454 W1: the h = 0.0075 exact-rectangle hole mesh (Gmsh 2.2), mesh only.

    A copy of ``export_hole_h001`` (commit d17b69d2, lines 3595-3792; not edited,
    not called) with ``resolution=H0075_RES`` -- ANS-6 step 5b's rung, the call
    ``_sheets_build(True, as_hole=True, resolution=0.0075,
    leg_sheet_strip_fraction=GATED_WIDTH_FRACTION)`` (test_th14_birdcage_copper.py
    :310-314).  (E6) reads the template actually used (h05x_cu_o2.json).
    """
    import dolfinx
    from mpi4py import MPI

    from tests.mesh.helpers import global_cell_tag_set
    from tests.mesh.test_birdcage_port_sheets import SHEET_IFACE
    from tests.mesh.test_birdcage_port_sheets import _build as _sheets_build
    from tests.mesh.test_birdcage_port_tags import LEG_COUNT
    from tests.validation import test_th14_birdcage_copper as th14
    from tests.validation.test_ans4_exact_strip_records import (
        SPEC_AREA_M2,
        SPEC_HEIGHT_M,
        SPEC_WIDTH_M,
    )
    from tests.validation.test_ops61_exact_strip_ports import STRIP_BAND
    from tests.validation.test_port_birdcage_lumped_column import (
        CONDUCTOR_CELL_TAG,
        PHANTOM_CELL_TAG,
    )
    from tests.validation.test_port_lumped_narrowed_sheet import GATED_WIDTH_FRACTION

    tag = "[OPS-61 2c W1]"
    comm = MPI.COMM_WORLD
    assert comm.size == 1, "export_hole_h0075 is a serial writer: run it at mpiexec -n 1"
    assert th14._sheets_build is _sheets_build
    t_start = time.perf_counter()
    wall_attr = _hole_wall_attr()
    ports_idx = list(range(1, LEG_COUNT + 1))

    # 1. never overwrite
    mesh_path = OUT_DIR / H0075_HOLE_MESH_NAME
    assert not mesh_path.exists(), f"{mesh_path.relative_to(REPO)} exists: never overwrite"
    tmpl = json.loads((OUT_DIR / H0075_HOLE_TEMPLATES["cu"]).read_text())
    tmpl_tags = sorted(t for m in tmpl["Domains"]["Materials"] for t in m["Attributes"])
    assert tuple(tmpl_tags) == T43_GLOBAL_TAGS, (tmpl_tags, T43_GLOBAL_TAGS)

    # 2. the rung's hole route with exact rectangles, as ANS-6 step 5b builds it.
    msh, cell_tags, facet_tags, diag, t_mesh = _sheets_build(
        True, as_hole=True, leg_sheet_strip_fraction=float(GATED_WIDTH_FRACTION),
        resolution=float(H0075_RES),
    )
    tdim = msh.topology.dim
    fdim = tdim - 1
    msh.topology.create_connectivity(fdim, tdim)
    msh.topology.create_entity_permutations()
    ncells = int(msh.topology.index_map(tdim).size_global)
    print(f"{tag} rung hole mesh (resolution {H0075_RES}) built in {t_mesh:.2f} s: {ncells} tets "
          f"(record {H0075_HOLE_CELL_RECORD})", flush=True)
    # 3. (E1) before anything is written
    assert ncells == H0075_HOLE_CELL_RECORD, f"(E1) tets {ncells} != {H0075_HOLE_CELL_RECORD}"

    _outer, census = th14._outer_box_tags(msh, facet_tags, comm)
    tags_f, _defs, specs = th14._build_ports(msh, cell_tags, comm)
    print(f"{tag} census {census}", flush=True)

    x = np.asarray(msh.geometry.x, dtype=float)

    def facet_nodes(facets):
        facets = np.asarray(facets, dtype=np.int32)
        if not facets.size:
            return np.zeros((0, 3), dtype=np.int64)
        return np.asarray(
            dolfinx.cpp.mesh.entities_to_geometry(msh._cpp_object, fdim, facets, False)
        ).reshape(-1, 3)

    # 4. triangle blocks: box -> 1, wall -> 401, N_i -> 20 + i.
    ext = np.asarray(dolfinx.mesh.exterior_facet_indices(msh.topology), dtype=np.int32)
    wall = np.asarray(facet_tags.find(wall_attr), dtype=np.int32)
    box = np.setdiff1d(ext, wall).astype(np.int32)
    box_nodes, wall_nodes = facet_nodes(box), facet_nodes(wall)
    box_area = float(np.sum(_tri_areas(x, box_nodes)))
    wall_area = float(np.sum(_tri_areas(x, wall_nodes)))
    tri_blocks = [(box_nodes, EXTERIOR_TAG), (wall_nodes, wall_attr)]
    port_rows = {}
    for i, spec in zip(ports_idx, specs):
        ptag = SHEET_IFACE + i
        assert int(spec.facet_tag) == ptag, (i, spec.facet_tag)
        narrow = np.asarray(tags_f.find(ptag), dtype=np.int32)
        n_nodes = facet_nodes(narrow)
        a_n = float(np.sum(_tri_areas(x, n_nodes)))
        port_rows[i] = {
            "facets": int(narrow.size),
            "area": a_n,
            "width": float(spec.sheet_width_m),
            "height": float(spec.gap_height_m),
        }
        tri_blocks.append((n_nodes, NARROW_BASE + i))

    OUT_DIR.mkdir(parents=True, exist_ok=True)
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
    total_vol = float(np.sum(np.abs(vol6)) / 6.0)

    ntri = sum(b.shape[0] for b, _ in tri_blocks)
    with open(mesh_path, "w") as fh:
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

    # (E1) again, parsed from the written file's $Elements type-4 count
    n_tet_file = 0
    n_el_file = 0
    with open(mesh_path) as fh:
        inside = False
        for line in fh:
            if line.startswith("$Elements"):
                inside = True
                n_el_file = int(next(fh))
                continue
            if line.startswith("$EndElements"):
                break
            if inside and line.split(None, 2)[1] == "4":
                n_tet_file += 1

    lo, hi = x.min(axis=0), x.max(axis=0)
    tags = sorted(global_cell_tag_set(msh, cell_tags))
    counts = {t: int(np.count_nonzero(ctag == t)) for t in tags}
    print(
        f"{tag} wrote {mesh_path.relative_to(REPO)}: {ncl} tets ({int(neg.sum())} reoriented), "
        f"$Elements count {n_el_file}, type-4 count {n_tet_file}; volume {total_vol:.9e} m^3; "
        f"bbox lo {lo} hi {hi}",
        flush=True,
    )
    print(f"{tag} cell tags (count): {counts}", flush=True)
    print(
        f"{tag} box (attr {EXTERIOR_TAG}): {box.size} tris area {box_area:.12e} m^2; "
        f"wall (attr {wall_attr}): {wall.size} tris area {wall_area:.12e} m^2",
        flush=True,
    )
    for i, r in port_rows.items():
        print(
            f"    P{i}: N_{i} {r['facets']} tris (attr {NARROW_BASE + i}) area "
            f"{r['area']:.12e} m^2 (rel {r['area'] / SPEC_AREA_M2 - 1.0:+.3e}); spec "
            f"sheet_width_m {r['width']:.12e} (rel {r['width'] / SPEC_WIDTH_M - 1.0:+.3e}) "
            f"gap_height_m {r['height']:.12e} (rel {r['height'] / SPEC_HEIGHT_M - 1.0:+.3e})",
            flush=True,
        )
    print(f"{tag} export total {time.perf_counter() - t_start:.2f} s (mesh {t_mesh:.2f} s)",
          flush=True)

    # 5. anchors (E1)-(E6), after writing.
    fails = []
    if ncl != H0075_HOLE_CELL_RECORD or ncells != ncl or n_tet_file != H0075_HOLE_CELL_RECORD:
        fails.append(f"(E1) tets {ncl} / file {n_tet_file} != {H0075_HOLE_CELL_RECORD}")
    if CONDUCTOR_CELL_TAG in tags or PHANTOM_CELL_TAG not in tags:
        fails.append(f"(E2) cell tags {tags}")
    if not (census["cavity_exterior"] == census["cavity"] > 0
            and census["outer"] == census["exterior"] - census["cavity"]
            and census["tag_collision"] == 0):
        fails.append(f"(E3) census {census}")
    if not abs(box_area - EXTERIOR_AREA_M2) <= EXTERIOR_AREA_BAND:
        fails.append(f"(E4) box area {box_area:.12e} vs {EXTERIOR_AREA_M2}")
    for i, r in port_rows.items():
        for k, (got, want) in {"area": (r["area"], SPEC_AREA_M2),
                               "width": (r["width"], SPEC_WIDTH_M),
                               "height": (r["height"], SPEC_HEIGHT_M)}.items():
            if not abs(got / want - 1.0) < STRIP_BAND:
                fails.append(f"(E5) P{i} {k} {got:.12e} vs {want:.12e}")
    if tuple(tags) != tuple(tmpl_tags) or tuple(tags) != T43_GLOBAL_TAGS:
        fails.append(f"(E6) cell-tag set {tags} != {H0075_HOLE_TEMPLATES['cu']} Materials "
                     f"{tmpl_tags} / T43_GLOBAL_TAGS")
    print(f"{tag} (E1) {ncl} == file {n_tet_file} == {H0075_HOLE_CELL_RECORD}; (E4) box area "
          f"|d| {abs(box_area - EXTERIOR_AREA_M2):.3e}; (E6) tags {tags} vs template {tmpl_tags} "
          f"vs T43_GLOBAL_TAGS {list(T43_GLOBAL_TAGS)}", flush=True)
    print(
        f"{tag} anchors (E1)-(E6) [box band {EXTERIOR_AREA_BAND:.0e}, strip band "
        f"{STRIP_BAND}]: {'PASS' if not fails else 'FAIL ' + '; '.join(fails)}",
        flush=True,
    )
    assert not fails, "; ".join(fails)


def hole_h0075_configs(col):
    """Item 454 W2 (cu) / the PEC half's W1 (pec): (E1') on the exported rung mesh;
    four single-excitation configs from item 241 / 234's x1 order-2 template
    (Model.Mesh, each port's Excitation and Problem.Output change, nothing else)."""
    from tests.validation.test_th14_birdcage_copper import COPPER_SIGMA

    try:
        from mpi4py import MPI

        assert MPI.COMM_WORLD.size == 1, "hole-h0075-configs runs at -n 1"
    except ImportError:
        pass
    tag = f"[OPS-61 2c configs {col}]"
    wall_attr = _hole_wall_attr()
    targets = {name: OUT_DIR / f"{name}.json" for name in H0075_HOLE_SETS[col]}
    for p in targets.values():
        assert not p.exists(), f"{p} exists: never overwrite"

    # (E1') tet count parsed from the Gmsh 2.2 $Elements block
    mesh_path = OUT_DIR / H0075_HOLE_MESH_NAME
    lines = mesh_path.read_text().splitlines()
    i0 = lines.index("$Elements")
    i1 = lines.index("$EndElements")
    nel = int(lines[i0 + 1])
    ntet = sum(1 for ln in lines[i0 + 2: i1] if ln.split()[1] == "4")
    assert i1 - i0 - 2 == nel, (nel, i1 - i0 - 2)
    print(f"{tag} {mesh_path.relative_to(REPO)}: {nel} elements, {ntet} type-4 tets", flush=True)
    print(f"{tag} (E1') tets {ntet} == {H0075_HOLE_CELL_RECORD}: "
          f"{'PASS' if ntet == H0075_HOLE_CELL_RECORD else 'FAIL'}", flush=True)
    assert ntet == H0075_HOLE_CELL_RECORD, f"(E1') {ntet} != {H0075_HOLE_CELL_RECORD}: stop"

    tpl_path = OUT_DIR / H0075_HOLE_TEMPLATES[col]
    tpl = json.loads(tpl_path.read_text())
    assert tpl["Model"]["Mesh"] == HOLE_MESH_NAME, tpl["Model"]["Mesh"]
    ports = tpl["Boundaries"]["LumpedPort"]
    assert [int(p["Index"]) for p in ports] == [1, 2, 3, 4], ports
    for p in ports:
        assert p["Attributes"] == [20 + int(p["Index"])], p
        assert p["Direction"] == "+Z" and float(p["R"]) == 50.0, p
        assert int(p["Excitation"]) == int(p["Index"]), p
    assert tpl["Solver"]["Order"] == 2, tpl["Solver"]["Order"]
    samples = tpl["Solver"]["Driven"]["Samples"]
    assert len(samples) == 1 and samples[0]["Freq"] == [0.064], samples
    b = tpl["Boundaries"]
    if col == "cu":
        assert b["PEC"]["Attributes"] == [1], b["PEC"]
        assert b["Conductivity"] == [{"Attributes": [wall_attr], "Conductivity": COPPER_SIGMA,
                                      "Permeability": 1.0}], b["Conductivity"]
    else:
        assert b["PEC"]["Attributes"] == [1, wall_attr], b["PEC"]
        assert "Conductivity" not in b, b.keys()
    print(f"{tag} template {tpl_path.name}: mesh {HOLE_MESH_NAME}, ports 21-24 +Z 50 ohm, "
          f"Order 2, Freq [0.064], PEC {b['PEC']['Attributes']}, Conductivity "
          f"{json.dumps(b.get('Conductivity'))}: PASS", flush=True)

    for name, k in H0075_HOLE_SETS[col].items():
        out = json.loads(json.dumps(tpl))
        out["Model"]["Mesh"] = H0075_HOLE_MESH_NAME
        for p in out["Boundaries"]["LumpedPort"]:
            p["Excitation"] = int(p["Index"]) == k
        out["Problem"]["Output"] = f"postpro_{name}"
        targets[name].write_text(json.dumps(out, indent=2) + "\n")

        # template diff: exactly Model.Mesh, the four Excitation keys, Problem.Output
        back = json.loads(targets[name].read_text())
        assert back["Model"]["Mesh"] == H0075_HOLE_MESH_NAME
        assert back["Problem"]["Output"] == f"postpro_{name}"
        exc = [p["Excitation"] for p in back["Boundaries"]["LumpedPort"]]
        assert exc == [i == k for i in (1, 2, 3, 4)], exc
        t2 = json.loads(json.dumps(tpl))
        b2 = json.loads(json.dumps(back))
        for d in (t2, b2):
            d["Model"]["Mesh"] = None
            d["Problem"]["Output"] = None
            for p in d["Boundaries"]["LumpedPort"]:
                p["Excitation"] = None
        assert t2 == b2, f"{name}: differs from its template beyond Mesh / Excitation / Output"
        print(f"{tag} wrote {targets[name].relative_to(REPO)} from {tpl_path.name}: "
              f"excitation port {k} only, Mesh {H0075_HOLE_MESH_NAME}, Output postpro_{name}; "
              f"template diff = Mesh, Excitation, Output only: PASS", flush=True)
    print(f"{tag} done: (E1') PASS; four configs written", flush=True)


def _ours_hole_rung_s(log, header_prefix, label):
    """A copy of ``_ours_rung_s`` (not edited): the header matched by ``startswith``
    on the stripped line (exactly one), rows ``f"{label}{k}k ="``."""
    import gzip

    path = log
    if path.exists():
        lines = path.read_text().splitlines()
    else:
        with gzip.open(str(path) + ".gz", "rt") as fh:
            lines = fh.read().splitlines()
    idx = [n for n, ln in enumerate(lines) if ln.strip().startswith(header_prefix)]
    assert len(idx) == 1, f"{len(idx)} headers {header_prefix!r}"
    tok = re.compile(r"[+-]\d\.\d+e[+-]\d+[+-]\d\.\d+e[+-]\d+j")
    rows = []
    for k, ln in enumerate(lines[idx[0] + 1: idx[0] + 5], start=1):
        body = re.sub(r"^(PASSED)+", "", ln.strip()).strip()
        assert body.startswith(f"{label}{k}k ="), (k, ln)
        vals = tok.findall(body)
        assert len(vals) == 4, (k, ln)
        rows.append([complex(v) for v in vals])
    return np.array(rows), idx[0] + 1


def _hole_h0075_palace_s(col, tag):
    s = np.full((4, 4), np.nan + 0j)
    present = []
    for name, k in H0075_HOLE_SETS[col].items():
        path = OUT_DIR / f"postpro_{name}" / "port-S.csv"
        if not path.exists():
            print(f"{tag} {name}: no {path.relative_to(REPO)} (not measured)", flush=True)
            continue
        sk, header = _read_port_s(path)
        cols = [c for c in range(sk.shape[1]) if not np.all(np.isnan(sk[:, c]))]
        labels = [h for h in header if re.search(r"S\[\d+\]\[\d+\]", h)]
        print(f"{tag} {name}: S columns present {[c + 1 for c in cols]}; labels "
              f"{labels[:2]} ...", flush=True)
        if (k - 1) in cols and sk.shape[0] == 4:
            colv = sk[:, k - 1]
            print(f"{tag} {name}: kept column labelled S[i][{k}] (asserted)", flush=True)
        else:
            assert len(cols) == 1 and sk.shape[0] == 4, (name, cols, sk.shape)
            colv = sk[:, cols[0]]
            print(f"{tag} {name}: READING -- Palace labelled the single excitation "
                  f"S[i][{cols[0] + 1}], kept by position as column {k}", flush=True)
        assert not np.any(np.isnan(colv)), (name, colv)
        s[:, k - 1] = colv
        present.append(k)
    return s, present


def read_hole_h0075(col):
    """Item 454 W7 (cu) / the PEC half's read (pec): a copy of
    ``read_exact_strip_h0075_freq`` (lines 2291-2402, not edited) with the hole
    set, ours' entry-31 block, the x1 Palace output and
    ``ANS6_EXACT_STRIP_CLASS_MEANS_P2["64"][col]`` substituted.  (Tie), (P),
    (K1') / (K3') asserted (cu); (K1) / (K3) asserted (pec, with all Cu columns);
    (T5) / (C-h) / (C2-h) / (T3m-h) / spreads / price printed, never gated."""
    from tests.validation.test_ans6_exact_strip import ANS6_EXACT_STRIP_CLASS_MEANS_P2
    from tests.validation.test_port21_feed_sensitivity import KNOB_REACHES_THE_SOLVE
    from tests.validation.test_port_birdcage_four_port import _circulant_classes

    blk = H0075_HOLE_OURS[col]
    tag = f"[OPS-61 2c read {col}]"

    def means(m):
        return {c: complex(np.mean(v)) for c, v in _circulant_classes(m).items()}

    def rel(a, b):
        return {c: abs(a[c] - b[c]) / abs(b[c]) for c in b}

    def line(d):
        return "  ".join(f"{c} {v * 100:.4f}% ({v:.3e})" for c, v in d.items())

    def sci(d):
        return "  ".join(f"{c} {v:.3e}" for c, v in d.items())

    fails = []
    ours = {}
    for c in ("cu", "pec"):
        b = H0075_HOLE_OURS[c]
        s_c, hdr_line = _ours_hole_rung_s(H0075_HOLE_OURS_LOG, b["header"], b["label"])
        m_c = means(s_c)
        tie_m = max(abs(m_c[k] - b["means"][k]) for k in b["means"])
        sig_c = float(np.linalg.svd(s_c, compute_uv=False).max())
        tie_s = abs(sig_c - b["sigma_max"]) if b["sigma_max"] is not None else 0.0
        ok = tie_m <= H0075_TIE_BAND and tie_s <= H0075_TIE_BAND
        print(
            f"{tag} ours {c} parsed from {H0075_HOLE_OURS_LOG.relative_to(REPO)}:{hdr_line} "
            f"(+4 lines, label {b['label']!r}): class means max|diff| {tie_m:.3e} vs log "
            f"10-digit means; sigma_max {sig_c:.12f}"
            + (f" vs printed {b['sigma_max']} |diff| {tie_s:.3e}" if b["sigma_max"] is not None
               else " (not printed at 10 digits)")
            + f" (Tie band {H0075_TIE_BAND:.0e}): {'PASS' if ok else 'FAIL'}",
            flush=True,
        )
        for r in s_c:
            print(f"    ours {c} " + "  ".join(_fmt(z) for z in r), flush=True)
        assert ok, f"(Tie) {c}: means {tie_m:.3e} sigma {tie_s:.3e}: wrong block parsed"
        ours[c] = {"s": s_c, "means": m_c, "sigma_max": sig_c,
                   "spreads": _identities(s_c)["spreads"]}

    s, present = _hole_h0075_palace_s(col, tag)
    print(f"{tag} Palace {col} 4x4 (columns present {present}):", flush=True)
    for r in s:
        print("    Palace " + "  ".join(_fmt(z) for z in r), flush=True)
    s_ours = ours[col]["s"]
    for k in present:
        g = np.abs(s[:, k - 1] - s_ours[:, k - 1]) / np.abs(s_ours[:, k - 1])
        print(f"{tag} (T5-{col}) column {k} per entry |Palace - ours|/|ours|: "
              + "  ".join(f"{v * 100:.4f}% ({v:.3e})" for v in g), flush=True)

    if len(present) == 4:
        ident, pf = _p_check_s(f"h05x_h0075_{col}_o2 (assembled)", s, tag)
        fails += [f"(P) {x}" for x in pf]
        pal = means(s)
        t5 = rel(pal, ours[col]["means"])
        print(f"{tag} (T5-{col}) per class |Palace o2 h - ours p2 h|/|ours p2 h|: {line(t5)} "
              "(predicted <= 1e-6; printed, never gated)", flush=True)
        print(f"{tag} _class_spread per class: Palace "
              + "  ".join(f"{c} {v * 100:.6f}%" for c, v in ident["spreads"].items())
              + "; ours " + "  ".join(f"{c} {v * 100:.6f}%" for c, v in ours[col]["spreads"].items())
              + " (printed; predicted equal to the printed digit)", flush=True)
        x1_path = OUT_DIR / blk["pal_x1"] / "port-S.csv"
        pal_x1 = means(_read_port_s(x1_path)[0])
        ch = rel(pal, pal_x1)
        print(f"{tag} (C-h-{col}) Palace o2 h=0.0075 vs Palace o2 x1 ({blk['pal_x1']}) per "
              f"class: {line(ch)}; ours' h-move "
              + "  ".join(f"{c} {v * 100:.4f}%" for c, v in blk["h_move"].items())
              + " (entry-31 log :4593 / :4594; printed)", flush=True)
        ours_x1 = dict(ANS6_EXACT_STRIP_CLASS_MEANS_P2["64"][col])
        num = {c: abs(ours_x1[c] - pal[c]) for c in pal}
        den = {c: abs(ours[col]["means"][c] - pal[c]) for c in pal}
        print(
            f"{tag} (C2-h) |ours x1 p2 - Palace o2 h| / |ours h p2 - Palace o2 h| per class: "
            + "  ".join(f"{c} {num[c]:.3e}/{den[c]:.3e} = "
                        f"{(num[c] / den[c]) if den[c] > 0 else float('inf'):.4g}x" for c in pal)
            + " (predicted >= 506x vs item 229's 25x bar; printed)",
            flush=True,
        )
        if col == "cu":
            # (K1') floor asserted: Palace Cu vs ours PEC on the rung
            k1 = rel(pal, ours["pec"]["means"])
            k1_ok = min(k1.values()) > KNOB_REACHES_THE_SOLVE
            sep = {c: abs(pal[c] - ours["pec"]["means"][c]) / abs(pal[c] - ours["cu"]["means"][c])
                   if abs(pal[c] - ours["cu"]["means"][c]) > 0 else float("inf") for c in pal}
            print(
                f"{tag} (K1') |Palace Cu h - ours PEC h|/|ours PEC h| {sci(k1)}; floor > "
                f"{KNOB_REACHES_THE_SOLVE:.0e} on every class: {'PASS' if k1_ok else 'FAIL'}; "
                f"ours' Cu-PEC gap {sci(H0075_HOLE_CU_PEC_GAP_OURS)} (:4592); separation "
                f"|Palace Cu - ours PEC|/|Palace Cu - ours Cu| "
                + "  ".join(f"{c} {v:.4g}x" for c, v in sep.items())
                + " (predicted >= 279x)",
                flush=True,
            )
            if not k1_ok:
                fails.append(f"(K1') min class move {min(k1.values()):.3e}")
            # (K3') loss direction asserted
            sp = ours["pec"]["sigma_max"]
            k3_ok = ident["sigma_max"] < sp
            print(
                f"{tag} (K3') sigma_max Palace Cu h {ident['sigma_max']:.12f} < ours PEC h "
                f"{sp:.12f} (margin {sp - ident['sigma_max']:.3e}; predicted >= 3.3e-6): "
                f"{'PASS' if k3_ok else 'FAIL'}",
                flush=True,
            )
            if not k3_ok:
                fails.append(f"(K3') sigma_max cu {ident['sigma_max']!r} >= ours pec {sp!r}")
        else:
            s_cu, pres_cu = _hole_h0075_palace_s("cu", tag + " [cu]")
            if len(pres_cu) == 4:
                id_cu = _identities(s_cu)
                pal_cu = id_cu["means"]
                k1 = rel(pal_cu, pal)
                k1_ok = min(k1.values()) > KNOB_REACHES_THE_SOLVE
                print(f"{tag} (K1) Palace Cu h vs Palace PEC h class move {sci(k1)}; floor > "
                      f"{KNOB_REACHES_THE_SOLVE:.0e}: {'PASS' if k1_ok else 'FAIL'}", flush=True)
                if not k1_ok:
                    fails.append(f"(K1) min class move {min(k1.values()):.3e}")
                k3_ok = id_cu["sigma_max"] < ident["sigma_max"]
                print(f"{tag} (K3) sigma_max Palace Cu h {id_cu['sigma_max']:.12f} < Palace PEC "
                      f"h {ident['sigma_max']:.12f}: {'PASS' if k3_ok else 'FAIL'}", flush=True)
                if not k3_ok:
                    fails.append("(K3) sigma_max cu >= pec")
                d_pal = {c: pal_cu[c] - pal[c] for c in pal}
                d_ours = {c: ours["cu"]["means"][c] - ours["pec"]["means"][c] for c in pal}
                mag = {c: abs(d_pal[c] - d_ours[c]) / abs(d_ours[c]) for c in pal}
                ang = {c: float(np.degrees(np.angle(d_pal[c] / d_ours[c]))) for c in pal}
                print(f"{tag} (T3m-h) |d_pal - d_ours|/|d_ours| {sci(mag)}; arg "
                      + "  ".join(f"{c} {v:+.4f} deg" for c, v in ang.items())
                      + " (predicted <= 1e-3, |arg| < 0.1 deg; printed)", flush=True)
            else:
                print(f"{tag} (K1) / (K3) / (T3m-h): not measured (Cu columns {pres_cu})",
                      flush=True)
    else:
        print(f"{tag} (P), per-class (T5), spreads, (C-h), (C2-h), (K1'), (K3'): not measured "
              f"(columns present {present})", flush=True)

    logs = sorted((REPO / "docs" / "testing" / "logs").glob(
        f"*_OPS-61-item*-W*-palace-hole-h0075-{col}-e*.log"))
    for lg in logs:
        txt = lg.read_text()
        tot = re.findall(r"^\s*Total\s+[\d.]+\s+[\d.]+\s+([\d.]+)\s*$", txt, re.M)
        setup = re.findall(r"^\s*Setup\s+[\d.]+\s+[\d.]+\s+([\d.]+)\s*$", txt, re.M)
        its = re.findall(r"GMRES solver converged in (\d+) iteration", txt)
        peak = re.findall(r"[Pp]eak.*?([\d.]+\s*[GM]B?)", txt)
        if not tot or not setup:
            print(f"{tag} price {lg.name}: no Total/Setup line", flush=True)
            continue
        print(f"{tag} price {lg.name}: measured Total {float(tot[-1]):.1f} s, Setup "
              f"{float(setup[-1]):.1f} s, GMRES its {its}, peak {peak[-1:] if peak else 'n/a'}; "
              f"brief's projection 115-180 s", flush=True)

    assert not fails, "; ".join(fails)
    print(f"{tag} done: (Tie) PASS" + ("; (P) PASS" if len(present) == 4 else "")
          + ("; (K1') PASS; (K3') PASS" if len(present) == 4 and col == "cu" else ""),
          flush=True)


if __name__ == "__main__":
    if sys.argv[1:] == ["h0050-configs", "--f64"]:
        h0050_f64_configs()
    elif sys.argv[1:] == ["read", "--exact-strip-h0050", "--f64"]:
        read_exact_strip_h0050_f64()
    elif sys.argv[1:] == ["export", "--hole-h0075"]:
        export_hole_h0075()
    elif sys.argv[1:] == ["hole-h0075-configs", "--cu"]:
        hole_h0075_configs("cu")
    elif sys.argv[1:] == ["hole-h0075-configs", "--pec"]:
        hole_h0075_configs("pec")
    elif sys.argv[1:] == ["read", "--hole-h0075", "--cu"]:
        read_hole_h0075("cu")
    elif sys.argv[1:] == ["read", "--hole-h0075", "--pec"]:
        read_hole_h0075("pec")
    elif sys.argv[1:] == ["ans7-h001-configs", "--f128"]:
        ans7_h001_128_configs()
    elif sys.argv[1:] == ["read", "--ans7-h001", "--f128"]:
        read_ans7_h001_128()
    elif sys.argv[1:] == ["export", "--hole-h001"]:
        export_hole_h001()
    elif sys.argv[1:] == ["ans7-h001-configs"]:
        ans7_h001_configs()
    elif sys.argv[1:] == ["read", "--ans7-h001"]:
        read_ans7_h001()
    elif sys.argv[1:] == ["ans7-128-configs"]:
        ans7_128_configs()
    elif sys.argv[1:] == ["export", "--exact-strip-h0050"]:
        export_exact_strip_h0050()
    elif sys.argv[1:] == ["read", "--exact-strip-h0050"]:
        read_exact_strip_h0050()
    elif sys.argv[1:] == ["read", "--ans7-128"]:
        read_ans7_128()
    elif sys.argv[1:] == ["ans7-sweep-configs"]:
        ans7_sweep_configs()
    elif sys.argv[1:] == ["read", "--ans7-sweep"]:
        read_ans7_sweep()
    elif sys.argv[1:] == ["ans7-configs"]:
        ans7_configs()
    elif sys.argv[1:] == ["read", "--ans7-esr"]:
        read_ans7_esr()
    elif sys.argv[1:] == ["h0075-configs", "--f64"]:
        h0075_freq_configs("64")
    elif sys.argv[1:] == ["h0075-configs", "--f128"]:
        h0075_freq_configs("128")
    elif sys.argv[1:] == ["read", "--exact-strip-h0075", "--f64"]:
        read_exact_strip_h0075_freq("64")
    elif sys.argv[1:] == ["read", "--exact-strip-h0075", "--f128"]:
        read_exact_strip_h0075_freq("128")
    elif sys.argv[1:] == ["export", "--exact-strip-h0075"]:
        export_exact_strip_h0075()
    elif sys.argv[1:] == ["read", "--exact-strip-h0075"]:
        read_exact_strip_h0075()
    elif sys.argv[1:] == ["cu-configs"]:
        cu_configs()
    elif sys.argv[1:] == ["read", "--hole-cu"]:
        read_hole_cu()
    elif sys.argv[1:] == ["export", "--hole-pec"]:
        export_hole_pec()
    elif sys.argv[1:] == ["read", "--hole-pec"]:
        read_hole_pec()
    elif sys.argv[1:] == ["o2-configs"]:
        o2_configs()
    elif sys.argv[1:] == ["read", "--exact-strip-o2"]:
        read_exact_strip_o2()
    elif sys.argv[1:] == ["freq-configs"]:
        freq_configs()
    elif sys.argv[1:] == ["read", "--exact-strip-set"]:
        read_exact_strip_set()
    elif sys.argv[1:] == ["export", "--exact-strip"]:
        export(exact_strip=True)
    elif sys.argv[1:] == ["flip-config"]:
        flip_config()
    elif sys.argv[1:] == ["read", "--exact-strip"]:
        read_exact_strip()
    elif len(sys.argv) != 2 or sys.argv[1] not in ("export", "read"):
        sys.exit("usage: ops61_step1a_palace_ans4.py export [--exact-strip] | read")
    else:
        {"export": export, "read": read}[sys.argv[1]]()
