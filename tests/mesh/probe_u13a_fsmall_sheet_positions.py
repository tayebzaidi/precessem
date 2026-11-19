"""`WF-16` U13a′ W0a′ (§9 item 1000; item 993's W0a re-posed) — F-small strip positions.

Not collected (``python_files = test_*.py``).  Measurement only, asserts
nothing.  Plain ``python3`` at ``-n 1`` (a gmsh throw under ``mpiexec -n 2``
can deadlock, `GEO-34` 1d brief).

1. Builds the F-small gapped four-leg rung twice and writes each OCC model as
   STEP under ``argv[1]`` (gitignored ``/workspace/logs/u13a_prime``):
   ``_build(False, step_export_path=<dir>/generator_fsmall_whole.step)`` —
   each port box **whole** (ruling (a); the input of the fixture) — and
   ``_build(True, step_export_path=<dir>/generator_fsmall_halves.step)`` —
   item 993's halves solids, which feed the (n-split) control only.  Strip
   flag off in both (it needs ``emit_port_sheets=True``); the OCC writer drops
   the free sheets (known-issues 2026-10-07).
2. Per build prints ``port_radius``, ``port_dx, port_dy, port_dz`` and per-leg
   ``box_angle`` with ``!r``, and the generator's CAD mass of each port group
   (``100 + i``; ``200 + i`` on the halves) beside dx·dy·dz.  Pre-gate (read by
   the slot): the flag-off port_radius / dx / dz equal W0a's
   (``20261115T210354Z_WF-16-U13a-W0a-positions.log:1796``) bit for bit.
3. In a fresh gmsh session, rebuilds each leg's ``leg_sheet_strip_fraction =
   0.5`` centre rectangle with the generator's own lines and prints its placed
   corners and area — the literals ``cq_author_u13a_fsmall_four_port.py``
   restates.
"""

import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))

import gmsh  # noqa: E402

from fem_em_solver.io.mesh import _place_sheet_in_leg_frame  # noqa: E402
from tests.mesh.test_birdcage_port_sheets import _build  # noqa: E402
from tests.mesh.test_birdcage_port_tags import LEG_COUNT  # noqa: E402

F_STRIP = 0.5


def _one(label, emit_port_sheets, out):
    t0 = time.perf_counter()
    _msh, _cells, _facets, diag, _elapsed = _build(emit_port_sheets, step_export_path=out)
    t_build = time.perf_counter() - t0
    layout = diag.get("port_layout")
    port_radius = None if layout is None else float(layout["port_radius_m"])
    dx, dy, dz = (float(v) for v in diag["port_box_size_m"])
    offsets = diag.get("leg_azimuth_offsets_rad")
    theta = np.linspace(0.0, 2.0 * np.pi, LEG_COUNT, endpoint=False)
    if offsets is not None:
        theta = theta + np.asarray(offsets, dtype=float)
    print(
        f"[W0a'] {label}: emit_port_sheets={emit_port_sheets} build + STEP export "
        f"{t_build:.2f} s -> {out} bytes={os.path.getsize(out)}",
        flush=True,
    )
    print(
        f"[W0a'] {label}: port_radius={port_radius!r} port_dx={dx!r} port_dy={dy!r} "
        f"port_dz={dz!r} offsets={offsets!r} "
        f"box_angles={[float(t) for t in theta]!r}",
        flush=True,
    )
    masses = diag.get("cad_mass_by_group") or {}
    for name in sorted(masses):
        if name.startswith("port_"):
            print(
                f"[W0a'] {label}: cad_mass {name}={masses[name]!r} " f"(dx*dy*dz={dx * dy * dz!r})",
                flush=True,
            )
    return port_radius, dx, dz, theta


def main(out_dir):
    os.makedirs(out_dir, exist_ok=True)
    whole = _one("whole", False, os.path.join(out_dir, "generator_fsmall_whole.step"))
    halves = _one("halves", True, os.path.join(out_dir, "generator_fsmall_halves.step"))
    port_radius, dx, dz, theta = whole
    if port_radius is None:
        port_radius = halves[0]
        print(
            "[W0a'] flag-off diag lacks port_layout: port_radius read from the "
            "flag-on (halves) build",
            flush=True,
        )
    w0a = (0.07, 0.014, 0.008)  # W0a :1796
    print(
        f"[W0a'] pre-gate: flag-off (port_radius, dx, dz)={(port_radius, dx, dz)!r} "
        f"W0a :1796={w0a!r} bitwise_equal={(port_radius, dx, dz) == w0a}",
        flush=True,
    )
    print(f"[W0a'] strip area 0.5*dx*dz={0.5 * dx * dz!r} full dx*dz={dx * dz!r}", flush=True)

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("u13a_strip_positions")
    for idx, box_angle in enumerate(theta):
        tag = gmsh.model.occ.addRectangle(
            port_radius - 0.5 * F_STRIP * dx, -0.5 * dz, 0.0, F_STRIP * dx, dz
        )
        _place_sheet_in_leg_frame(tag, float(box_angle))
        gmsh.model.occ.synchronize()
        pts = gmsh.model.getBoundary([(2, tag)], combined=False, recursive=True)
        coords = [tuple(float(v) for v in gmsh.model.getValue(0, p, [])) for _, p in pts]
        area = float(gmsh.model.occ.getMass(2, tag))
        bb = tuple(float(v) for v in gmsh.model.occ.getBoundingBox(2, tag))
        print(
            f"[W0a'] leg {idx + 1}: box_angle={float(box_angle)!r} " f"area={area!r} bbox={bb!r}",
            flush=True,
        )
        for k, c in enumerate(coords):
            print(f"[W0a'] leg {idx + 1} corner {k}: {c!r}", flush=True)
    gmsh.finalize()


if __name__ == "__main__":
    main(sys.argv[1])
