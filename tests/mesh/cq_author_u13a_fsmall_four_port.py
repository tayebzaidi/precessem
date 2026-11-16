"""`WF-16` U13a′ W0b′ (§9 item 1000; item 993's W0b re-posed) — author the F-small four-port STEP.

Not collected (``python_files = test_*.py``).  Runs under ``cq-python``, which
clears ``PYTHONPATH`` (``docker/Dockerfile:41``), so the generator cannot be
imported and its numbers are restated as literals from W0a
(``docs/testing/logs/20261115T210354Z_WF-16-U13a-W0a-positions.log``),
re-confirmed bit for bit on the whole-box build by W0a′
(``docs/testing/logs/20261116T030547Z_WF-16-U13a-W0a-positions-whole.log``):

* W0a ``:1796`` / W0a′ ``:3274`` port_radius 0.07, port_dx 0.014, port_dz 0.008;
* W0a ``:1797`` / W0a′ ``:3275`` strip area 0.5*dx*dz = 5.6e-05;
* W0a ``:1799–1817`` / W0a′ ``:3276–3295`` the four placed
  ``leg_sheet_strip_fraction = 0.5`` strip corners.

The solids are the generator's own STEP export, re-read with
``cq.importers.importStep``: for the fixture, W0a′'s whole-box export
(``_build(False, ...)``, each port box one solid — ruling (a)); for the
(n-split) control only, its halves export (``_build(True, ...)``).  The four
strips are explicit 4-corner polygon faces (never ``makePlane``), products
``PortSheet1`` … ``PortSheet4`` of one ``cq.Assembly`` beside ``FsmallSolids``.

Usage: ``cq-python tests/mesh/cq_author_u13a_fsmall_four_port.py <out.step>
<generator.step>`` (both required; no default path).
"""

import hashlib
import os
import sys
import time

import cadquery as cq

STRIP_AREA = 5.6e-05  # W0a :1797
R_IN, R_OUT, Z_HALF = 0.0665, 0.07350000000000001, 0.004  # W0a :1799-1817
# (x, y) of the inner / outer radial corner per leg, azimuth 90 deg * (i - 1).
LEG_XY = {
    1: ((R_IN, 0.0), (R_OUT, 0.0)),
    2: ((0.0, R_IN), (0.0, R_OUT)),
    3: ((-R_IN, 0.0), (-R_OUT, 0.0)),
    4: ((0.0, -R_IN), (0.0, -R_OUT)),
}


def _strip(i):
    (x0, y0), (x1, y1) = LEG_XY[i]
    corners = [
        cq.Vector(x0, y0, -Z_HALF),
        cq.Vector(x1, y1, -Z_HALF),
        cq.Vector(x1, y1, Z_HALF),
        cq.Vector(x0, y0, Z_HALF),
    ]
    face = cq.Face.makeFromWires(cq.Wire.makePolygon(corners, close=True))
    area = face.Area()
    bb = face.BoundingBox()
    print(f"[author] PortSheet{i} area={area!r} bbox=({bb.xmin!r}, {bb.ymin!r}, "
          f"{bb.zmin!r}, {bb.xmax!r}, {bb.ymax!r}, {bb.zmax!r})", flush=True)
    assert abs(area - STRIP_AREA) <= 1e-12, (i, area)
    want = (min(x0, x1), min(y0, y1), -Z_HALF, max(x0, x1), max(y0, y1), Z_HALF)
    got = (bb.xmin, bb.ymin, bb.zmin, bb.xmax, bb.ymax, bb.zmax)
    # OCC pads bounding boxes by the shape tolerance (~1e-7): check at 1e-6.
    for g, w in zip(got, want):
        assert abs(g - w) <= 1e-6, (i, g, w)
    return face


def main(out, generator_step):
    print("cadquery", cq.__version__, flush=True)
    t0 = time.perf_counter()
    solids_wp = cq.importers.importStep(generator_step)
    solids = solids_wp.solids().vals()
    vol = sum(s.Volume() for s in solids)
    print(f"[author] importStep {generator_step}: {len(solids)} solids, "
          f"total volume={vol!r}, faces={len(solids_wp.faces().vals())}", flush=True)
    print(f"[author] imported solid count={len(solids)} (whole expected 30 by "
          f"arithmetic, halves 34; printed, not gated)", flush=True)

    asm = cq.Assembly(name="u13a_fsmall")
    asm.add(solids_wp, name="FsmallSolids")
    for i in sorted(LEG_XY):
        asm.add(cq.Workplane().add(_strip(i)), name=f"PortSheet{i}")
    if hasattr(asm, "export"):
        asm.export(out)
        route = "Assembly.export"
    else:
        asm.save(out, exportType="STEP")
        route = "Assembly.save (export absent; disclosed)"
    with open(out, "rb") as handle:
        data = handle.read()
    print(f"[author] route={route} author seconds={time.perf_counter() - t0:.2f}",
          flush=True)
    print(f"[author] wrote {out} bytes={len(data)} "
          f"sha256={hashlib.sha256(data).hexdigest()}", flush=True)
    print(f"[author] size_on_disk={os.path.getsize(out)}", flush=True)

    back = cq.importers.importStep(out)
    print(f"[author] read-back: {len(back.solids().vals())} solids, "
          f"{len(back.faces().vals())} faces", flush=True)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
