"""POST-8 step 1 — `ans:2`'s pointwise sample at one width, and the per-cell
N1curl values behind it (PROJECT_PLAN §9 item 54; known-issues 2026-09-27).

Run at ``-n 2`` and ``-n 4`` on one commit (complex build,
``FEM_EM_REQUIRE_COMPLEX=1``).  Measures and fixes nothing: no band, no record,
no sampling rule (the weekly's call).

What it does
------------
* Loads ``ans:2``'s script as a module (importlib, as that script loads its own
  dependencies) and calls its pieces — ``_build_mass_averaged`` on the unset
  resolution knob (0.0025 m record rung), ``_pointwise_sar``, ``_phantom_power``
  — **never** ``main()``, so no tracked ``metrics.json`` / ``COMPARISON.md`` and
  no private comparison is written.  The probe writes nothing under
  ``examples/``.
* **Anchor (asserted, exit non-zero on failure):**
  (1) the P1 pointwise SAR row reproduces *this width's* logged row to 1e-9
      relative (``-n 2``: ``20260926T140900Z_EX-57.log:1909``; ``-n 4``:
      ``20260925T110348Z_OPS-63.log:1909`` = tracked ``metrics.json:49–94``);
  (2) whole-phantom P1 power is within 1e-9 of ``metrics.json:100`` — the
      reduced quantity that must stay still where the sample moves (the
      negative control for (1)).
* **Printed (the mechanism reading):** at each named point, every incident
  cell (rank-local ``bb_tree`` collisions on *owned* cells, global cell
  indices, gathered to rank 0) with the raw N1curl ``E`` for drives P1…P4
  evaluated in that cell (a per-cell ``eval`` by construction — the probe's
  measurement, not a sampling route), beside the CG1 value that
  ``_pointwise_sar`` samples.  For P1 also the CG1 nodal values at the
  vertices of each incident cell next to that cell's own N1curl value at the
  vertex (the interpolation-overwrite hypothesis, read directly).
* **Census (printed):** every file under ``examples/`` and ``tests/`` that
  point-samples (``evaluate_vector_field_parallel`` / ``.eval(``), with a
  *text heuristic* for whether the same file builds an N1curl field and a
  Lagrange interpolant.  Classification gated / printed / exported is by
  heuristic (``tests/`` ⇒ gated, ``examples/`` ⇒ printed unless the site is in
  an export helper) and is flagged as such for the weekly to confirm.
"""

from __future__ import annotations

import importlib.util
import re
import sys
import time
from pathlib import Path

import numpy as np
from mpi4py import MPI

from dolfinx import default_scalar_type, fem, geometry

REPO = Path(__file__).resolve().parents[2]
ANS2 = (
    REPO
    / "examples"
    / "ansys_benchmarks"
    / "ans2_birdcage_coil_driven_sar_10MHz"
    / "02_birdcage_coil_driven_sar_10MHz.py"
)

#: Logged P1 rows (W/kg), per MPI width.  Printed to 10 significant figures,
#: so the 1e-9 relative anchor carries ≤ 1.1e-10 of rounding.
LOGGED_P1_ROW = {
    2: [4.858355042e-07, 6.058341318e-09, 5.710133879e-08, 8.049438449e-09],
    4: [5.878191741e-07, 1.045326320e-08, 6.558122284e-08, 9.290633415e-09],
}
LOGGED_P1_SOURCE = {
    2: "20260926T140900Z_EX-57.log:1909",
    4: "20260925T110348Z_OPS-63.log:1909 (tracked metrics.json:49-94)",
}
#: ``metrics.json:100`` (written at -n 4).
TRACKED_P1_POWER_W = 5.503204728560125e-08
ROW_RTOL = 1e-9
POWER_RTOL = 1e-9


def _load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def _per_cell(built, points, comm, v_cg_funcs):
    """Every owned incident cell of each point, with N1curl and CG1 values."""
    msh = built["mesh"]
    tdim = msh.topology.dim
    imap = msh.topology.index_map(tdim)
    n_owned = imap.size_local
    tree = geometry.bb_tree(msh, tdim)
    cand = geometry.compute_collisions_points(tree, points)
    coll = geometry.compute_colliding_cells(msh, cand, points)
    geom_x = msh.geometry.x
    geom_dofmap = msh.geometry.dofmap
    records = []
    for i in range(points.shape[0]):
        for c in coll.links(i):
            c = int(c)
            if c >= n_owned:
                continue
            gidx = int(imap.local_to_global(np.array([c], dtype=np.int32))[0])
            cell = np.array([c], dtype=np.int32)
            x = points[i : i + 1]
            n1 = []
            cg = []
            for k in range(4):
                e = built["solves"][f"P{k + 1}"]["fields"].e_complex
                n1.append(np.asarray(e.eval(x, cell)).reshape(3))
                cg.append(np.asarray(v_cg_funcs[k].eval(x, cell)).reshape(3))
            # P1 vertex readout: CG1 nodal value vs this cell's N1curl value
            verts = []
            e1 = built["solves"]["P1"]["fields"].e_complex
            for gd in geom_dofmap[c]:
                xv = geom_x[gd : gd + 1, :3].copy()
                n1v = np.asarray(e1.eval(xv, cell)).reshape(3)
                cgv = np.asarray(v_cg_funcs[0].eval(xv, cell)).reshape(3)
                verts.append((xv.reshape(3), n1v, cgv))
            records.append(
                {
                    "point": i,
                    "cell": gidx,
                    "rank": comm.rank,
                    "n1": n1,
                    "cg": cg,
                    "verts": verts,
                }
            )
    return records


def _census():
    pat = re.compile(r"evaluate_vector_field_parallel\(|\.eval\(")
    rows = []
    for root in ("examples", "tests"):
        for f in sorted((REPO / root).rglob("*.py")):
            text = f.read_text(errors="replace")
            lines = text.splitlines()
            hits = [
                n + 1
                for n, ln in enumerate(lines)
                if pat.search(ln) and not ln.lstrip().startswith("#")
            ]
            if not hits:
                continue
            n1 = bool(re.search(r"N1curl|Nedelec|e_complex|build_.*sweep|_solve_driven", text))
            lag = bool(re.search(r"\"Lagrange\"|'Lagrange'|\"CG\"|\"DG\"", text))
            interp = ".interpolate(" in text
            candidate = n1 and lag and interp
            klass = "gated (heuristic)" if root == "tests" else "printed (heuristic)"
            rows.append((str(f.relative_to(REPO)), hits, n1, lag, interp, candidate, klass))
    return rows


def main() -> int:
    comm = MPI.COMM_WORLD
    t0 = time.perf_counter()
    if not np.issubdtype(np.dtype(default_scalar_type), np.complexfloating):
        raise RuntimeError("needs the complex DolfinX build")
    width = comm.size
    if width not in LOGGED_P1_ROW:
        raise RuntimeError(f"POST-8 step 1 runs at -n 2 or -n 4, not -n {width}")

    ans2 = _load(ANS2, "post8_ans2_module")
    if ans2.RESOLUTION_OVERRIDDEN:
        raise RuntimeError("ans:2 resolution knob must be unset (record rung)")

    built = ans2._build_mass_averaged(phantom_resolution=ans2.PHANTOM_RESOLUTION_M)
    r = ans2.CENTRE_RADIUS_M
    points = np.array(
        [[+r, 0.0, 0.0], [0.0, +r, 0.0], [-r, 0.0, 0.0], [0.0, -r, 0.0]],
        dtype=np.float64,
    )
    pointwise = ans2._pointwise_sar(built, points, comm)
    powers = ans2._phantom_power(built)

    # CG1 interpolants rebuilt exactly as _pointwise_sar builds them
    v_cg = fem.functionspace(built["mesh"], ("Lagrange", 1, (3,)))
    cg_funcs = []
    for k in range(4):
        f = fem.Function(v_cg)
        f.interpolate(built["solves"][f"P{k + 1}"]["fields"].e_complex)
        f.x.scatter_forward()
        cg_funcs.append(f)
    records = _per_cell(built, points, comm, cg_funcs)
    gathered = comm.gather(records, root=0)

    sig = float(ans2.SALINE_SIGMA)
    rho = float(ans2.PHANTOM_RHO_KG_PER_M3)

    def sar(v):
        return sig * float(np.sum(np.abs(v) ** 2)) / (2.0 * rho)

    failed = False
    if comm.rank == 0:
        print("=" * 78, flush=True)
        print(
            f"[POST-8] step 1 — ans:2 pointwise sample at mpiexec -n {width}; "
            f"cells {built['cells']}, phantom cells {built['phantom_cells']}, "
            f"phantom_resolution {ans2.PHANTOM_RESOLUTION_M} m",
            flush=True,
        )
        # --- anchor (1)
        row = pointwise[0]["sar_w_per_kg"]
        ref = LOGGED_P1_ROW[width]
        rel = [abs(a / b - 1.0) for a, b in zip(row, ref)]
        ok1 = max(rel) <= ROW_RTOL
        print(f"[POST-8] ANCHOR (1) -n {width}: P1 row vs {LOGGED_P1_SOURCE[width]}")
        for j in range(4):
            print(f"         point {j}: {row[j]:.12e} vs logged {ref[j]:.9e}  rel {rel[j]:.3e}")
        print(f"[POST-8] ANCHOR (1) max rel {max(rel):.3e} <= {ROW_RTOL:g}: {'PASS' if ok1 else 'FAIL'}")
        other = 4 if width == 2 else 2
        print(
            f"[POST-8]   (printed) same row vs the -n {other} logged row: rel "
            + "  ".join(f"{a / b - 1.0:+.4e}" for a, b in zip(row, LOGGED_P1_ROW[other]))
        )
        for k in range(4):
            print(
                f"         drive P{k + 1} row (-n {width}): "
                + "  ".join(f"{v:.9e}" for v in pointwise[k]["sar_w_per_kg"])
            )
        # --- anchor (2)
        p1 = powers[0]["dissipated_power_w"]
        rel2 = abs(p1 / TRACKED_P1_POWER_W - 1.0)
        ok2 = rel2 <= POWER_RTOL
        print(
            f"[POST-8] ANCHOR (2) -n {width}: P1 phantom power {p1:.15e} W vs "
            f"metrics.json:100 {TRACKED_P1_POWER_W:.15e}  rel {rel2:.3e} <= "
            f"{POWER_RTOL:g}: {'PASS' if ok2 else 'FAIL'}"
        )
        failed = not (ok1 and ok2)

        # --- per-cell readout
        recs = [r_ for part in gathered for r_ in part]
        recs.sort(key=lambda d: (d["point"], d["cell"]))
        print(
            f"\n[POST-8] PER-CELL READOUT (-n {width}): owned incident cells per point; "
            "SAR = sigma|E|^2/(2 rho) from the raw N1curl E evaluated in that cell "
            "vs the CG1 sample _pointwise_sar reported",
            flush=True,
        )
        for i in range(4):
            mine = [d for d in recs if d["point"] == i]
            print(f"  point {i} {points[i].tolist()}: {len(mine)} incident cell(s)")
            for k in range(4):
                cg_sample = pointwise[k]["sar_w_per_kg"][i]
                n1_sars = [sar(d["n1"][k]) for d in mine]
                lo, hi = (min(n1_sars), max(n1_sars)) if n1_sars else (np.nan, np.nan)
                inside = lo <= cg_sample <= hi if n1_sars else False
                print(
                    f"    P{k + 1}: CG1 sample {cg_sample:.9e}; N1curl per-cell "
                    f"[{lo:.9e}, {hi:.9e}]; CG1 inside per-cell range: {inside}"
                )
                for d in mine:
                    e = d["n1"][k]
                    print(
                        f"        cell {d['cell']:>7d} (rank {d['rank']}) N1curl SAR "
                        f"{sar(e):.9e}  CG1-in-cell SAR {sar(d['cg'][k]):.9e}  "
                        f"E_n1 = [{', '.join(f'{z.real:+.4e}{z.imag:+.4e}j' for z in e)}]"
                    )
            for d in mine:
                print(f"      P1 vertices of cell {d['cell']}: |E_n1(cell)|^2 vs |E_cg1 node|^2")
                for xv, n1v, cgv in d["verts"]:
                    a = float(np.sum(np.abs(n1v) ** 2))
                    b = float(np.sum(np.abs(cgv) ** 2))
                    print(
                        f"        x=({xv[0]:+.5f},{xv[1]:+.5f},{xv[2]:+.5f})  "
                        f"{a:.6e}  {b:.6e}  ratio {b / a if a else np.nan:.4f}"
                    )

        # --- census
        print("\n[POST-8] CENSUS — point-sampling sites (text heuristic; weekly to confirm)")
        print("  file | lines | N1curl-ish | Lagrange | interpolate | CANDIDATE | class")
        cands = 0
        for f, hits, n1, lag, interp, cand, klass in _census():
            cands += cand
            print(
                f"  {f} | {','.join(map(str, hits[:8]))}{'...' if len(hits) > 8 else ''} "
                f"| {n1} | {lag} | {interp} | {'CANDIDATE' if cand else '-'} | {klass}"
            )
        print(f"[POST-8] census: {cands} candidate file(s) (N1curl + Lagrange + interpolate + point sample)")
        print(f"[POST-8] elapsed {time.perf_counter() - t0:.1f} s at -n {width}", flush=True)
        print(f"[POST-8] RESULT: {'FAIL' if failed else 'PASS'}", flush=True)
    failed = comm.bcast(failed, root=0)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
