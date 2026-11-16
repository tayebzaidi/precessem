"""Measurement-only probe (§9 item 182, `TH-20` step 0a re-pose, chain H13-0a):
where does the quarter-vs-full cavity-wall area gap sit, and does it follow the
facet-count (triangulation) model or a fixed face mismatch (cutter defect)?

Asserts nothing, imported by nothing.  One fresh process per rung, ``-n 1``:

    python3 -u tests/mesh/probe_th20_quarter_wall_area.py --conductor-resolution 1.6e-3

Build: ``tests/validation/test_th20_quarter_f_human_mesh._kwargs()`` plus
``conductor_resolution=h``; the full hole mesh first, then
``symmetry_quadrant=True`` (`GEO-23`).  Wall = exterior tag 401 on both.

Restriction of the full wall to a quadrant is **exact triangle clipping**
(Sutherland–Hodgman against s_x·x ≥ 0, then s_y·y ≥ 0, clipped on sign so the
leg at x ≈ 9.2e-18 m is handled by the clip, not by a tolerance), for all four
quadrants Qᵢ.  Q₁ = {x ≥ 0, y ≥ 0} is the kept quadrant.

Bands, assigned by the centroid of each clipped piece (quarter facets lie in
Q₁ already and are classified by their own centroid).  Precedence R > N > S >
F > O:
  R  neighbour cell tag in the kept ring-port set 117–120, 133–136, 217–220,
     233–236;
  N  min(x, y) ≤ 2h;
  S  the rest of the bisected sectors, θ < π/16 or θ > 7π/16;
  F  π/16 ≤ θ ≤ 7π/16, air (tag 2) neighbour;
  O  anything else (printed so nothing is silently dropped).

Per band: ΔA_b = A_q,b − A_f∩Q₁,b (m²), ΔA_b/A_f∩Q₁,b, ΔN_b/N_b with N the
facet (piece) count, and the band's facet-model share
model_b = |δ_f| · A_f∩Q₁,b · ΔN_b/N_b (m²).

Whole wall: A_f, A_q (ds-assembled), δ = A/A_CAD − 1 with the generator's CAD
cavity area; r = 4A_q/A_f − 1; model = |δ_f|·(4N_q/N_f − 1); ρ = r/model;
ΔA = A_q − A_f∩Q₁ (m²).

Negative controls (§9 rule (e)), all predicted and printed, none asserted:
  P1  Σᵢ A_f∩Qᵢ / A_f − 1 (A_f the facet-geometry sum), predicted ≤ 1e-12.
  P2  max_i |4 A_f∩Qᵢ / A_f − 1|, the full mesh's own spread (no cutter);
      magnitude unmeasured, printed.
  P3  planted missing face: A_q minus tag 117's ring-port wall
      (predicted ≈ 4.9 / 4.8 / 4.5e-5 m² at 1.6 / 2.4 / 3.2 mm); the decision
      rule below, applied to the P3 series, must read CUTTER DEFECT.

Decision rule (pre-stated, item 182; evaluated over the rungs present, each
rung's record kept as JSON under the gitignored /workspace/logs/th20_probe/):
  TRIANGULATION iff all three: ρ ∈ [0.5, 2] at every rung; ΔA changes sign or
    (max − min)/|mean| > 0.30 across the ladder; band R agrees within its own
    facet-model share, operationalised here as
    |ΔA_R − model_R| ≤ 0.5·|model_R| + 1e-8·A_R at every rung.
  CUTTER DEFECT iff either: ΔA has one sign and max |ΔA_i/mean − 1| ≤ 0.20 at
    every rung; or at every rung > 80 % of (ΔA − model·A_f/4) sits in band N
    while |ΔN_N/N_N| ≤ 1e-3 ("ΔN_N ≈ 0").
  Mixed otherwise (or both).  With fewer than three rungs the reading is
  printed as provisional.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import resource
import time

import numpy as np
from mpi4py import MPI

CAVITY_TAG = 401
AIR = 2
R_TAGS = frozenset(
    list(range(117, 121)) + list(range(133, 137)) + list(range(217, 221)) + list(range(233, 237))
)
P3_TAG = 117
QUADS = {1: (1.0, 1.0), 2: (-1.0, 1.0), 3: (-1.0, -1.0), 4: (1.0, -1.0)}
BANDS = ("R", "N", "S", "F", "O")
RECORD_DIR = "/workspace/logs/th20_probe"


def _clip(poly, k, s):
    """Sutherland–Hodgman: keep s·p[k] ≥ 0 (on sign)."""
    out = []
    n = len(poly)
    for i in range(n):
        p, q = poly[i - 1], poly[i]
        fp, fq = s * p[k], s * q[k]
        pin, qin = fp >= 0.0, fq >= 0.0
        if qin:
            if not pin:
                out.append(p + (q - p) * (fp / (fp - fq)))
            out.append(q)
        elif pin:
            out.append(p + (q - p) * (fp / (fp - fq)))
    return out


def _poly_area_centroid(poly):
    p0 = poly[0]
    a_tot = 0.0
    c = np.zeros(3)
    for i in range(1, len(poly) - 1):
        a = 0.5 * float(np.linalg.norm(np.cross(poly[i] - p0, poly[i + 1] - p0)))
        a_tot += a
        c += a * (p0 + poly[i] + poly[i + 1]) / 3.0
    return a_tot, (c / a_tot if a_tot > 0 else np.mean(poly, axis=0))


def _band(c, nbr, h):
    if int(nbr) in R_TAGS:
        return "R"
    x, y = float(c[0]), float(c[1])
    if min(x, y) <= 2.0 * h:
        return "N"
    th = math.atan2(y, x)
    if th < math.pi / 16 or th > 7 * math.pi / 16:
        return "S"
    if int(nbr) == AIR:
        return "F"
    return "O"


def _wall(msh, facets, cell_value):
    """Owned tag-401 triangles (n, 3, 3) and their neighbour cell tag."""
    import dolfinx

    tdim = msh.topology.dim
    fdim = tdim - 1
    msh.topology.create_connectivity(fdim, tdim)
    n_owned = msh.topology.index_map(fdim).size_local
    idx = np.asarray(facets.indices)[np.asarray(facets.values) == CAVITY_TAG]
    idx = idx[idx < n_owned]
    f2c = msh.topology.connectivity(fdim, tdim)
    nbr = cell_value[f2c.array[f2c.offsets[idx]]]
    g = dolfinx.cpp.mesh.entities_to_geometry(msh._cpp_object, fdim, idx.astype(np.int32), False)
    tri = msh.geometry.x[g][:, :3, :]
    return tri, nbr


def _tri_area(tri):
    return 0.5 * np.linalg.norm(np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]), axis=1)


def _quadrant_pieces(tri, nbr, q):
    """Exact clip of every triangle to quadrant q: list of (area, centroid, nbr)."""
    sx, sy = QUADS[q]
    vx, vy = sx * tri[:, :, 0], sy * tri[:, :, 1]
    inside = np.all(vx >= 0.0, axis=1) & np.all(vy >= 0.0, axis=1)
    outside = np.all(vx < 0.0, axis=1) | np.all(vy < 0.0, axis=1)
    whole = np.flatnonzero(inside)
    a = _tri_area(tri[whole])
    cen = tri[whole].mean(axis=1)
    pieces_a = list(a)
    pieces_c = list(cen)
    pieces_n = list(nbr[whole])
    n_straddle = 0
    for i in np.flatnonzero(~inside & ~outside):
        poly = [tri[i, 0].copy(), tri[i, 1].copy(), tri[i, 2].copy()]
        poly = _clip(poly, 0, sx)
        if len(poly) >= 3:
            poly = _clip(poly, 1, sy)
        if len(poly) < 3:
            continue
        ar, c = _poly_area_centroid(poly)
        if ar <= 0.0:
            continue
        n_straddle += 1
        pieces_a.append(ar)
        pieces_c.append(c)
        pieces_n.append(nbr[i])
    return np.asarray(pieces_a), np.asarray(pieces_c).reshape(-1, 3), np.asarray(pieces_n), n_straddle


def _band_sums(areas, cents, nbrs, h):
    out = {b: [0.0, 0] for b in BANDS}
    for a, c, n in zip(areas, cents, nbrs):
        b = _band(c, n, h)
        out[b][0] += float(a)
        out[b][1] += 1
    return out


def _build(kwargs, comm):
    from fem_em_solver.io.mesh import MeshGenerator

    t0 = time.perf_counter()
    msh, cells, facets, diag = MeshGenerator.birdcage_port_domain(
        comm=comm, return_diagnostics=True, **kwargs)
    t_build = time.perf_counter() - t0
    tdim = msh.topology.dim
    msh.topology.create_entities(tdim - 1)
    msh.topology.create_connectivity(tdim - 1, tdim)
    return msh, cells, facets, diag, t_build


def _readout(records, tag, key_da, key_rho, key_r_ok, key_n_ok):
    hs = sorted(records)
    da = np.array([records[h][key_da] for h in hs])
    rho = np.array([records[h][key_rho] for h in hs])
    mean = float(np.mean(da))
    sign_change = bool(np.any(da > 0) and np.any(da < 0))
    spread = float((da.max() - da.min()) / abs(mean)) if mean != 0 else float("inf")
    t1 = bool(np.all((rho >= 0.5) & (rho <= 2.0)))
    t2 = sign_change or spread > 0.30
    t3 = bool(all(records[h][key_r_ok] for h in hs))
    one_sign = bool(np.all(da > 0) or np.all(da < 0))
    dev = float(np.max(np.abs(da / mean - 1.0))) if mean != 0 else float("inf")
    d1 = one_sign and dev <= 0.20
    d2 = bool(all(records[h][key_n_ok] for h in hs))
    tri = t1 and t2 and t3
    cut = d1 or d2
    verdict = "TRIANGULATION" if (tri and not cut) else "CUTTER DEFECT" if (cut and not tri) else "Mixed"
    prov = "" if len(hs) >= 3 else f" (PROVISIONAL: {len(hs)} of 3 rungs)"
    print(f"{tag} rungs {['%.1e' % h for h in hs]} dA {['%+.4e' % v for v in da]} "
          f"rho {['%+.3f' % v for v in rho]}", flush=True)
    print(f"{tag} TRIANGULATION criteria: rho in [0.5,2] all rungs {t1}; dA sign change "
          f"{sign_change} or spread (max-min)/|mean| {spread:.3f} > 0.30 -> {t2}; band R "
          f"within facet-model share all rungs {t3}  => {tri}", flush=True)
    print(f"{tag} CUTTER DEFECT criteria: one sign {one_sign} and max|dA/mean-1| {dev:.3f} "
          f"<= 0.20 -> {d1}; >80% of (dA-model) in band N with dN_N~0 all rungs {d2}  "
          f"=> {cut}", flush=True)
    print(f"{tag} READING: {verdict}{prov}", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--conductor-resolution", type=float, required=True)
    args = ap.parse_args()
    h = args.conductor_resolution

    from tests.mesh.probe_two_torus_conductor_hole import _exterior_facet_area
    from tests.mesh.probe_wf12_cavity_wall_area import _cell_value, _class_of, _owned_facet_split
    from tests.validation.test_th20_quarter_f_human_mesh import _kwargs

    comm = MPI.COMM_WORLD
    tag = f"[TH-20 probe h={h:.2e}]"
    if comm.size != 1:
        if comm.rank == 0:
            print(f"{tag} clipping is rank-local; run at -n 1 (got -n {comm.size})", flush=True)
        return

    meshes = {}
    for name, quadrant in (("full", False), ("quarter", True)):
        kwargs = _kwargs()
        kwargs["conductor_resolution"] = h
        if quadrant:
            kwargs["symmetry_quadrant"] = True
        try:
            msh, cells, facets, diag, t_build = _build(kwargs, comm)
        except Exception as exc:  # gmsh FAIL is recorded verbatim, not hidden
            print(f"{tag} {name}: BUILD FAIL {type(exc).__name__}: {exc}", flush=True)
            return
        tdim = msh.topology.dim
        n_cells = int(msh.topology.index_map(tdim).size_global)
        rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0**2
        cv = _cell_value(msh, cells)
        tri, nbr = _wall(msh, facets, cv)
        A = _exterior_facet_area(msh, facets, CAVITY_TAG, comm)
        cad = float(diag["conductor_cavity_cad"]["area_m2"])
        split = _owned_facet_split(
            msh,
            np.asarray(facets.indices)[np.asarray(facets.values) == CAVITY_TAG],
            nbr, comm)
        meshes[name] = dict(tri=tri, nbr=nbr, A=A, cad=cad, n_cells=n_cells, N=len(tri))
        print(f"{tag} {name}: cells {n_cells}; build {t_build:.1f} s; gmsh mesh "
              f"{diag.get('mesh_wall_time_s', float('nan')):.1f} s; peak RSS {rss:.3f} GiB; "
              f"401 facets {len(tri)}; A {A:.10e} m^2; A_CAD {cad:.10e}; "
              f"delta {A / cad - 1.0:+.6e}", flush=True)
        for c in sorted(split):
            print(f"{tag}   {name} class {c}: area {split[c][0]:.10e} m^2, facets {split[c][1]}",
                  flush=True)
        del msh, cells, facets

    f, q = meshes["full"], meshes["quarter"]
    A_f, A_q = f["A"], q["A"]
    d_f, d_q = A_f / f["cad"] - 1.0, A_q / q["cad"] - 1.0
    r = 4.0 * A_q / A_f - 1.0
    model = abs(d_f) * (4.0 * q["N"] / f["N"] - 1.0)
    rho = r / model if model != 0 else float("nan")

    A_f_geo = float(_tri_area(f["tri"]).sum())
    A_q_geo = float(_tri_area(q["tri"]).sum())
    quad = {}
    for k in QUADS:
        a, c, n, ns = _quadrant_pieces(f["tri"], f["nbr"], k)
        quad[k] = dict(A=float(a.sum()), a=a, c=c, n=n, straddle=ns)
    dA = A_q - quad[1]["A"]

    print(f"{tag} A_f {A_f:.10e} (geo {A_f_geo:.10e}) A_q {A_q:.10e} (geo {A_q_geo:.10e}); "
          f"delta_f {d_f:+.6e} delta_q {d_q:+.6e}", flush=True)
    print(f"{tag} N_f {f['N']} N_q {q['N']} 4N_q/N_f-1 {4.0 * q['N'] / f['N'] - 1.0:+.6e}; "
          f"4q/f-1 {r:+.6e}; model |delta_f|(4N_q/N_f-1) {model:+.6e}; rho {rho:+.4f}", flush=True)
    for k in QUADS:
        print(f"{tag} clip Q{k}: A_f∩Q {quad[k]['A']:.10e} pieces {len(quad[k]['a'])} "
              f"(straddling triangles clipped {quad[k]['straddle']}); 4A/A_f-1 "
              f"{4.0 * quad[k]['A'] / A_f_geo - 1.0:+.6e}", flush=True)
    print(f"{tag} dA = A_q - A_f∩Q1 = {dA:+.6e} m^2; dA/(A_f/4) {dA / (A_f / 4):+.6e}; "
          f"model in m^2 {model * A_f / 4:+.6e}", flush=True)

    bq = _band_sums(_tri_area(q["tri"]), q["tri"].mean(axis=1), q["nbr"], h)
    bf = _band_sums(quad[1]["a"], quad[1]["c"], quad[1]["n"], h)
    bands = {}
    for b in BANDS:
        aq, nq = bq[b]
        af, nf = bf[b]
        da_b = aq - af
        dn = (nq - nf) / nf if nf else float("nan")
        mod_b = abs(d_f) * af * dn if nf else 0.0
        bands[b] = dict(dA=da_b, A=af, dN=dn, model=mod_b, Nq=nq, Nf=nf)
        print(f"{tag} band {b}: A_q {aq:.10e} A_f∩Q1 {af:.10e} dA {da_b:+.6e} m^2 dA/A "
              f"{(da_b / af if af else float('nan')):+.6e}; N_q {nq} N_f∩Q1 {nf} dN/N "
              f"{dn:+.6e}; facet-model share {mod_b:+.6e} m^2", flush=True)
    sum_band = sum(v["dA"] for v in bands.values())
    print(f"{tag} sum of band dA {sum_band:+.6e} vs dA {dA:+.6e} (|diff| {abs(sum_band - dA):.2e})",
          flush=True)
    other_ring = sorted({int(t) for t in quad[1]["n"]
                         if _class_of(int(t)).startswith("ring") and int(t) not in R_TAGS})
    print(f"{tag} Q1 pieces with ring-port tags outside the kept R set: {other_ring}", flush=True)

    R = bands["R"]
    r_ok = abs(R["dA"] - R["model"]) <= 0.5 * abs(R["model"]) + 1e-8 * R["A"]
    resid = dA - model * A_f / 4.0
    N = bands["N"]
    frac_N = (N["dA"] - N["model"]) / resid if resid != 0 else float("nan")
    n_ok = bool(frac_N > 0.8 and abs(N["dN"]) <= 1e-3)
    print(f"{tag} band R: |dA_R - model_R| {abs(R['dA'] - R['model']):.3e} vs "
          f"0.5|model_R| + 1e-8 A_R {0.5 * abs(R['model']) + 1e-8 * R['A']:.3e} -> agrees {r_ok}",
          flush=True)
    print(f"{tag} band N: (dA_N - model_N)/(dA - model) {frac_N:+.4f}; dN_N/N_N {N['dN']:+.3e} "
          f"-> defect-criterion (b) at this rung {n_ok}", flush=True)

    # ---- negative controls
    p1 = sum(quad[k]["A"] for k in QUADS) / A_f_geo - 1.0
    p2 = max(abs(4.0 * quad[k]["A"] / A_f_geo - 1.0) for k in QUADS)
    a117 = float(_tri_area(q["tri"][q["nbr"] == P3_TAG]).sum())
    n117 = int(np.count_nonzero(q["nbr"] == P3_TAG))
    A_q3 = A_q - a117
    r3 = 4.0 * A_q3 / A_f - 1.0
    rho3 = r3 / model if model != 0 else float("nan")
    dA3 = A_q3 - quad[1]["A"]
    print(f"{tag} (P1) sum_i A_f∩Qi/A_f - 1 = {p1:+.3e} (predicted |.| <= 1e-12)", flush=True)
    print(f"{tag} (P2) max_i |4 A_f∩Qi/A_f - 1| = {p2:.6e} (no cutter; magnitude unmeasured)",
          flush=True)
    print(f"{tag} (P3) planted missing face: tag-117 wall {a117:.6e} m^2 ({n117} facets; "
          f"predicted ~4.9/4.8/4.5e-5); dA_P3 {dA3:+.6e} m^2 ({abs(dA3) / abs(dA) if dA else float('nan'):.2f}x |dA|); "
          f"4q'/f-1 {r3:+.6e}; rho_P3 {rho3:+.4f}", flush=True)

    rec = dict(h=h, dA=dA, rho=rho, r_ok=bool(r_ok), n_ok=n_ok,
               dA_P3=dA3, rho_P3=rho3, r_ok_P3=False, n_ok_P3=False,
               n_full=f["n_cells"], n_quarter=q["n_cells"], N_f=f["N"], N_q=q["N"],
               A_f=A_f, A_q=A_q)
    os.makedirs(RECORD_DIR, exist_ok=True)
    with open(os.path.join(RECORD_DIR, f"rung_{h:.1e}.json"), "w") as fh:
        json.dump(rec, fh)
    records = {}
    for fn in sorted(os.listdir(RECORD_DIR)):
        if fn.startswith("rung_") and fn.endswith(".json"):
            with open(os.path.join(RECORD_DIR, fn)) as fh:
                d = json.load(fh)
            records[float(d["h"])] = d
    print(f"{tag} decision rule over rungs on record ({RECORD_DIR}):", flush=True)
    _readout(records, f"{tag} [measured]", "dA", "rho", "r_ok", "n_ok")
    # P3: band R is the planted face, so it cannot agree with its facet model;
    # criterion (b) is not the planted mechanism — both set False on purpose.
    _readout(records, f"{tag} [P3 control, must read CUTTER DEFECT]", "dA_P3", "rho_P3",
             "r_ok_P3", "n_ok_P3")


if __name__ == "__main__":
    main()
