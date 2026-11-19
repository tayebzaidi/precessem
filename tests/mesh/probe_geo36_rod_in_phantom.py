"""GEO-36 step 1 (§9 item 150): rod-in-phantom meshability and price over r.

Question: does a capped rod cut as a hole from the F-small phantom mesh with
graded tips, and what does it cost, across the rod radius?

Reduced geometry (the `mesh-probe` may not touch `src/` and no rod generator
exists): the F-small phantom alone -- a cylinder R = 0.03 m, H = 0.08 m
(`tests/mesh/test_birdcage_port_tags.py:28-29`), centred on the origin, axis
z -- with a rod cut out of it as a hole.  The rod: radius r, axis z, centred
at (0.015, 0, 0), a cylindrical section of length ROD_LENGTH = 30 mm between
its two hemispherical caps (so the closed forms are pi r^2 l + 4 pi r^3 / 3
for the volume and 2 pi r l + 4 pi r^2 for the wall).

Sizing: phantom at `RESOLUTION` 0.015 (MeshSizeMax); Distance -> Threshold from
the rod wall (SizeMin r/3 -> SizeMax 0.015 over DistMin r -> DistMax 10 r);
`Ball` fields at the two cap centres (radius 2 r, VIn r/4, VOut 0.015); the
background field is the Min of the three.  The rod-wall ("cavity") surfaces
are read off the meshed volume's boundary: every boundary surface of the cut
volume whose bounding box lies inside the rod's.

Swept: r in {0.50, 0.75, ..., 2.00} mm, uniform grid, plus one phantom-only
baseline rung (no rod, no fields) so "tets above baseline" can be read.

**One fresh subprocess per rung** -- a failed gmsh build contaminates the
process (`GEO-23`), so the parent only launches children and prints what they
print.  gmsh only, `-n 1` semantics (no MPI at all).  Measurement only:
nothing asserts, nothing imports this file, nothing is re-recorded.

    python3 -u tests/mesh/probe_geo36_rod_in_phantom.py
"""

import math
import subprocess
import sys
import time

PHANTOM_RADIUS = 0.03  # m, test_birdcage_port_tags.py:28
PHANTOM_HEIGHT = 0.08  # m, test_birdcage_port_tags.py:29
RESOLUTION = 0.015  # m, test_birdcage_port_tags.py:32
ROD_CENTRE = (0.015, 0.0, 0.0)
ROD_LENGTH = 0.030  # m, cylindrical section between the cap centres
RADII_MM = [0.50, 0.75, 1.00, 1.25, 1.50, 1.75, 2.00]
BASELINE = "baseline"


def _closed_forms(r):
    big = math.pi * PHANTOM_RADIUS**2 * PHANTOM_HEIGHT
    if r is None:
        return big, None
    rod_vol = math.pi * r**2 * ROD_LENGTH + 4.0 * math.pi * r**3 / 3.0
    rod_area = 2.0 * math.pi * r * ROD_LENGTH + 4.0 * math.pi * r**2
    return big - rod_vol, rod_area


def run_rung(arg):
    """Build and mesh one rung in this (fresh) process; print one row."""
    import gmsh
    import numpy as np

    r = None if arg == BASELINE else float(arg) * 1e-3
    label = "phantom-only BASELINE (no rod)" if r is None else f"r = {float(arg):.2f} mm"
    vol_cf, area_cf = _closed_forms(r)

    gmsh.initialize(["-noenv"])
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.option.setNumber("General.NumThreads", 1)
    gmsh.logger.start()
    try:
        gmsh.model.add("geo36")
        occ = gmsh.model.occ
        phantom = occ.addCylinder(
            0.0, 0.0, -0.5 * PHANTOM_HEIGHT, 0.0, 0.0, PHANTOM_HEIGHT, PHANTOM_RADIUS
        )
        rod_bbox = None
        if r is not None:
            cx, cy, cz = ROD_CENTRE
            z0 = cz - 0.5 * ROD_LENGTH
            z1 = cz + 0.5 * ROD_LENGTH
            shaft = occ.addCylinder(cx, cy, z0, 0.0, 0.0, ROD_LENGTH, r)
            cap0 = occ.addSphere(cx, cy, z0, r)
            cap1 = occ.addSphere(cx, cy, z1, r)
            rod, _ = occ.fuse([(3, shaft)], [(3, cap0), (3, cap1)])
            out, _ = occ.cut([(3, phantom)], rod)
            occ.synchronize()
            vols = [t for d, t in out if d == 3]
            pad = 1e-6
            rod_bbox = (
                cx - r - pad,
                cy - r - pad,
                z0 - r - pad,
                cx + r + pad,
                cy + r + pad,
                z1 + r + pad,
            )
        else:
            occ.synchronize()
            vols = [phantom]
        masses = [occ.getMass(3, v) for v in vols]
        occ_mass = sum(masses)
        print(
            f"  [{label}] census: cut volumes={len(vols)} occ mass={occ_mass:.9e} "
            f"(closed form {vol_cf:.9e})",
            flush=True,
        )

        boundary = gmsh.model.getBoundary([(3, v) for v in vols], oriented=False)
        bsurfs = sorted({t for d, t in boundary if d == 2})
        rod_surfs = []
        if rod_bbox is not None:
            for s in bsurfs:
                b = gmsh.model.getBoundingBox(2, s)
                if all(b[i] >= rod_bbox[i] for i in range(3)) and all(
                    b[i + 3] <= rod_bbox[i + 3] for i in range(3)
                ):
                    rod_surfs.append(s)
        occ_area = sum(occ.getMass(2, s) for s in rod_surfs)

        gmsh.option.setNumber("Mesh.MeshSizeMax", RESOLUTION)
        if r is not None:
            gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
            gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
            gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
            f = gmsh.model.mesh.field
            dist = f.add("Distance")
            f.setNumbers(dist, "SurfacesList", rod_surfs)
            f.setNumber(dist, "Sampling", 50)
            thr = f.add("Threshold")
            f.setNumber(thr, "InField", dist)
            f.setNumber(thr, "SizeMin", r / 3.0)
            f.setNumber(thr, "SizeMax", RESOLUTION)
            f.setNumber(thr, "DistMin", r)
            f.setNumber(thr, "DistMax", 10.0 * r)
            balls = []
            for zc in (ROD_CENTRE[2] - 0.5 * ROD_LENGTH, ROD_CENTRE[2] + 0.5 * ROD_LENGTH):
                b = f.add("Ball")
                f.setNumber(b, "XCenter", ROD_CENTRE[0])
                f.setNumber(b, "YCenter", ROD_CENTRE[1])
                f.setNumber(b, "ZCenter", zc)
                f.setNumber(b, "Radius", 2.0 * r)
                f.setNumber(b, "VIn", r / 4.0)
                f.setNumber(b, "VOut", RESOLUTION)
                balls.append(b)
            mn = f.add("Min")
            f.setNumbers(mn, "FieldsList", [thr, *balls])
            f.setAsBackgroundMesh(mn)

        t0 = time.time()
        gmsh.model.mesh.generate(3)
        gsec = time.time() - t0

        tet_tags, _ = gmsh.model.mesh.getElementsByType(4)
        ntet = len(tet_tags)
        q = np.asarray(gmsh.model.mesh.getElementQualities(tet_tags, "minSICN"))
        qmin = float(q.min())
        qp1 = float(np.percentile(q, 1.0))

        mesh_area = 0.0
        cy_acc = cz_acc = 0.0
        all_tags, all_xyz, _ = gmsh.model.mesh.getNodes()
        lookup = dict(zip(np.asarray(all_tags).tolist(), range(len(all_tags))))
        all_xyz = np.asarray(all_xyz).reshape(-1, 3)
        for s in rod_surfs:
            _, ntags = gmsh.model.mesh.getElementsByType(2, s)
            if len(ntags) == 0:
                continue
            idx = np.array([lookup[int(n)] for n in ntags]).reshape(-1, 3)
            tri = all_xyz[idx]
            a = 0.5 * np.linalg.norm(np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]), axis=1)
            c = tri.mean(axis=1)
            mesh_area += a.sum()
            cy_acc += (a * c[:, 1]).sum()
            cz_acc += (a * c[:, 2]).sum()

        row = (
            f"ROW {label:<30s} OK    tets={ntet:>8d}  gmsh={gsec:7.2f} s  "
            f"SICN min={qmin:.4f} p1={qp1:.4f}  "
            f"occ mass rel err={(occ_mass - vol_cf) / vol_cf:+.3e}"
        )
        if r is not None:
            row += (
                f"  rod surfs={rod_surfs}"
                f"  wall area occ rel err={(occ_area - area_cf) / area_cf:+.3e}"
                f"  mesh rel err={(mesh_area - area_cf) / area_cf:+.3e}"
                f"  facet centroid y={cy_acc / mesh_area:+.3e} z={cz_acc / mesh_area:+.3e}"
            )
        print(row, flush=True)
    except Exception as exc:  # gmsh raises a bare Exception carrying its text
        print(f"ROW {label:<30s} FAIL  {str(exc).strip()}", flush=True)
        for line in gmsh.logger.get():
            if any(k in line for k in ("overlapping", "MeshAdapt", "equivalent", "Error", "rror")):
                print(f"  [{label}] gmsh: {line}", flush=True)
    finally:
        try:
            gmsh.logger.stop()
            gmsh.finalize()
        except Exception:
            pass


def main():
    print("GEO-36 step 1: capped rod cut from the F-small phantom, swept over r")
    print(
        f"  phantom R={PHANTOM_RADIUS} m H={PHANTOM_HEIGHT} m; RESOLUTION={RESOLUTION} m; "
        f"rod centre={ROD_CENTRE}, cylindrical length {ROD_LENGTH} m + 2 hemispherical caps"
    )
    print("  one fresh python subprocess per rung; measurement only\n", flush=True)
    for arg in [BASELINE] + [f"{r:.2f}" for r in RADII_MM]:
        t0 = time.time()
        proc = subprocess.run(
            [sys.executable, "-u", __file__, "--rung", arg],
            capture_output=True,
            text=True,
        )
        wall = time.time() - t0
        sys.stdout.write(proc.stdout)
        if proc.returncode != 0:
            print(f"  [rung {arg}] subprocess exit {proc.returncode}; stderr tail:")
            for line in proc.stderr.strip().splitlines()[-5:]:
                print(f"    {line}")
        print(f"  [rung {arg}] subprocess wall {wall:.2f} s", flush=True)
    print("\nMeasurement only -- no assertion, nothing re-recorded.")


if __name__ == "__main__":
    if len(sys.argv) == 3 and sys.argv[1] == "--rung":
        run_rung(sys.argv[2])
    else:
        main()
