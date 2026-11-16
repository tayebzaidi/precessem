"""Example (`EX-102`): refine the mesh at an implant tip; the CAD does not move.

For the Phase 7 user who wants a finer mesh at a lead or rod tip without
touching the geometry: ``implant_rod["tip_size"]`` [m] is the cap-``Ball``
``VIn`` of the rod generator (unset = ``radius / 4``). The fixture is `GEO-36`
step 4's, **imported**: the four-port F-small birdcage with a 1 mm-radius,
32 mm tip-to-tip capped rod cut from the phantom, at the 1 g phantom rung. It
is built three times, coarse to fine (`GEO-23`): ``tip_size`` unset, ``r/(4
sqrt 2)`` and ``r/8``. No solve.

Everything gated is imported from the gate module
``tests/mesh/test_birdcage_implant_rod_tip.py`` (the `ANS-1` rule): bands,
records, rungs, fixture keywords and the summary helper ``_summarise`` (which
prints the gate's own ``[GEO-36 step 4]`` lines and raises at the runaway
guard). Two inputs are restated, each cited where used: the refusal values of
(c1) and the conductor cell tag.

Asserted per rung: (N) the unset cell count against ``ROD_1G_RUNG_CELL_RECORD``
within ``CELL_COUNT_BAND``; (A0) CAD masses, OCC wall area and wall centroid
equal the unset build's, ``tip_size_m`` equals the request; (A1) / (A3) within
``CAD_ANALYTIC_BAND``; (A4) in ``TERMINAL_AREA_BAND``; (A5) the cell census and
the tag-402 facet census; (A6) the wall centroid within ``EXACT * ROD_LENGTH``;
(D) N1curl degree 2 = 2E + 2F and degree 1 = E; (G) N_ball strictly increasing.
Control (c1): two out-of-range ``tip_size`` values are refused before gmsh
starts. (c2) the refinement is local (air count flat, phantom count grows) is
printed, never asserted.

**Not claimed:** a solve, eta_rod or tip SAR, an h-convergence statement, the
off (no-rod) build's A2, a ball-radius knob, the kept-cells build.

Run it through the example runner::

    ./scripts/run_examples.sh -e mesh:14

Outputs (full filenames, ``examples/meshing/paraview_output/``):
``meshing_14_implant_tip_refinement_unset_combined.xdmf``,
``meshing_14_implant_tip_refinement_r4rt2_combined.xdmf``,
``meshing_14_implant_tip_refinement_r8_combined.xdmf``.
"""

from __future__ import annotations

import re
import sys
import time
from pathlib import Path

import gmsh
import numpy as np
from mpi4py import MPI

import dolfinx
from dolfinx import fem

_REPO_ROOT = Path(__file__).resolve().parents[2]
# The runner puts only ``src`` on PYTHONPATH; the repo root goes on sys.path so
# the gate's constants and helpers can be imported rather than restated.
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from fem_em_solver.io.paraview_utils import (  # noqa: E402
    adopt_host_ownership,
    write_xdmf_with_tags,
)
from fem_em_solver.post.setup_figure import (  # noqa: E402
    figures_enabled,
    write_setup_figure,
)

from tests.mesh.test_birdcage_conductor_hole import (  # noqa: E402
    CAD_ANALYTIC_BAND,
    EXACT,
)
from tests.mesh.test_birdcage_implant_rod import (  # noqa: E402
    ROD,
    ROD_CENTRE,
    ROD_LENGTH,
    ROD_RADIUS,
    _build_rod,
    _phantom_minus_rod,
    _rod_area,
)
from tests.mesh.test_birdcage_implant_rod_tip import (  # noqa: E402
    AIR_CELL_TAG,
    FIXTURE_KW,
    PHANTOM_CELL_TAG,
    PLANS,
    ROD_1G_RUNG_CELL_RECORD,
    ROD_CELL_GUARD_W2_PRINTED,
    RUNAWAY_CELL_GUARD,
    T1_RADIUS,
    TIP_R4,
    _cell_geometry,
    _mean_edge,
    _summarise,
)
from tests.mesh.test_birdcage_leg_gaps import TERMINAL_AREA_BAND  # noqa: E402
from tests.mesh.test_birdcage_port_sheet_prerequisite import CELL_COUNT_BAND  # noqa: E402
from tests.mesh.test_birdcage_port_sheets import PORT_LOWER, PORT_UPPER  # noqa: E402
from tests.mesh.test_birdcage_port_tags import LEG_COUNT  # noqa: E402

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"  # committed, EX-57
BASENAME = "meshing_14_implant_tip_refinement"
#: filename-safe rung labels (``/`` and ``(`` cannot go in a filename)
FILE_LABEL = {"unset": "unset", "r/(4rt2)": "r4rt2", "r/8": "r8"}
#: the conductor cell tag, as the gate builds its expected set
#: (`test_birdcage_implant_rod_tip.py:318-320`)
CONDUCTOR_CELL_TAG = 1
#: (c1) refusal inputs, the gate's own (`test_birdcage_implant_rod_tip.py:265`):
#: coarser than the r/4 cap, and zero
REFUSED_TIP_SIZES = (ROD_RADIUS / 3.0, 0.0)


def _say(msg):
    if MPI.COMM_WORLD.rank == 0:
        print(msg, flush=True)


def _ball_mask(msh, cell_tags, ball_radius):
    """Owned-cell mask of the N_ball set (the gate's `_ball_and_t1` definition)."""
    n_owned = msh.topology.index_map(msh.topology.dim).size_local
    centroids = _cell_geometry(msh).mean(axis=1)
    tag_of = np.zeros(n_owned, dtype=np.int64)
    idx = np.asarray(cell_tags.indices)
    vals = np.asarray(cell_tags.values)
    keep = idx < n_owned
    tag_of[idx[keep]] = vals[keep]
    rod_c = np.asarray(ROD["centre"], dtype=float)
    rod_a = np.asarray(ROD["axis"], dtype=float)
    mask = np.zeros(n_owned, dtype=bool)
    for sign in (-1.0, 1.0):
        cap = rod_c + sign * 0.5 * ROD_LENGTH * rod_a
        mask |= np.linalg.norm(centroids - cap, axis=1) <= ball_radius
    return mask & (tag_of == PHANTOM_CELL_TAG)


def _write_rung_xdmf(label, out, ball_radius, n_ball_gate, comm):
    """Combined XDMF for one rung: mean_edge and in_tip_ball (DG0) + facet tags."""
    msh, cells, facets = out[0], out[1], out[2]
    tdim = msh.topology.dim
    n_owned = msh.topology.index_map(tdim).size_local
    mask = _ball_mask(msh, cells, ball_radius)
    # tie the script's own mask to the gate's reduced N_ball
    n_ball_here = comm.allreduce(int(np.count_nonzero(mask)), op=MPI.SUM)
    assert n_ball_here == n_ball_gate, (
        f"in_tip_ball field {n_ball_here} != _summarise N_ball {n_ball_gate}"
    )
    v0 = fem.functionspace(msh, ("DG", 0))
    mean_edge = fem.Function(v0, name="mean_edge")
    mean_edge.x.array[:] = 0.0
    mean_edge.x.array[:n_owned] = _mean_edge(_cell_geometry(msh))
    mean_edge.x.scatter_forward()
    in_ball = fem.Function(v0, name="in_tip_ball")
    in_ball.x.array[:] = 0.0
    in_ball.x.array[:n_owned] = mask.astype(float)
    in_ball.x.scatter_forward()
    path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_{FILE_LABEL[label]}_combined",
        msh,
        cells,
        {"mean_edge": mean_edge, "in_tip_ball": in_ball},
        comm=comm,
        facet_tags=facets,
    )
    return path


def _display_submesh_figure(out, ball_radius, comm):
    """Setup figure on a display-only phantom submesh (never enters an assertion)."""
    msh, cells = out[0], out[1]
    tdim = msh.topology.dim
    imap = msh.topology.index_map(tdim)
    n_parent = imap.size_local + imap.num_ghosts
    parent_values = np.full(n_parent, -1, dtype=np.int32)
    parent_values[cells.indices] = cells.values
    keep = np.flatnonzero(parent_values == PHANTOM_CELL_TAG).astype(np.int32)
    sub, _emap, _vmap, _gmap = dolfinx.mesh.create_submesh(msh, tdim, keep)
    n_sub_owned = sub.topology.index_map(tdim).size_local
    mids = dolfinx.mesh.compute_midpoints(sub, tdim, np.arange(n_sub_owned, dtype=np.int32))
    rod_c = np.asarray(ROD["centre"], dtype=float)
    rod_a = np.asarray(ROD["axis"], dtype=float)
    near = np.zeros(n_sub_owned, dtype=bool)
    for sign in (-1.0, 1.0):
        cap = rod_c + sign * 0.5 * ROD_LENGTH * rod_a
        near |= np.linalg.norm(mids - cap, axis=1) <= ball_radius
    display = dolfinx.mesh.meshtags(
        sub,
        tdim,
        np.arange(n_sub_owned, dtype=np.int32),
        np.where(near, 2, 1).astype(np.int32),
    )
    write_setup_figure(
        sub,
        display,
        FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names={
            1: "phantom (gel), outside the tip balls",
            2: "tip-ball cells, tip_size r/8",
        },
        translucent_tags=(1,),
        slice_normal=(0.0, 1.0, 0.0),
        slice_origin=ROD_CENTRE,
        title=("mesh:14 -- rod tip refinement: phantom cells within 2r of the caps "
               "at tip_size r/8 (GEO-36 step 4)"),
        comm=comm,
    )


def main() -> None:
    comm = MPI.COMM_WORLD
    rank = comm.rank
    started = time.perf_counter()
    OUTPUT_DIR.mkdir(exist_ok=True)

    _say("=" * 72)
    _say("EX-102 -- refine the mesh at an implant tip: the CAD does not move")
    _say("=" * 72)
    _say(f"[fixture] GEO-36 step 4: F-small birdcage + rod {ROD}, FIXTURE_KW={FIXTURE_KW}; "
         f"ranks={comm.size}")

    # ---- (c1) out-of-range tip_size refused before gmsh starts (every rank) --
    refusals = []
    for value in REFUSED_TIP_SIZES:
        t0 = time.perf_counter()
        raised = None
        try:
            _build_rod(implant_rod=dict(ROD, tip_size=value), **FIXTURE_KW)
        except ValueError as e:
            raised = str(e)
        refusals.append((value, raised, time.perf_counter() - t0))
    gmsh_untouched = not gmsh.isInitialized()
    for value, raised, dt in refusals:
        _say(f"[c1] tip_size={value!r}: refused={raised is not None} in {dt:.3f} s "
             f"({(raised or '')[:90]!r})")
    _say(f"[c1] gmsh initialised after the refusals: {not gmsh_untouched}")

    # ---- the three rungs, coarse to fine (GEO-23) ----------------------------
    rungs = PLANS["rungs"]
    summaries, files = [], []
    kept = None
    for i, (label, tip) in enumerate(rungs):
        out = _build_rod(implant_rod=ROD if tip is None else dict(ROD, tip_size=tip),
                         **FIXTURE_KW)
        s = _summarise(label, tip, out, comm)
        summaries.append(s)
        ball_radius = s["diag"]["implant_rod_cad"]["tip_ball_radius_m"]
        files.append(_write_rung_xdmf(label, out, ball_radius, s["n_ball"], comm))
        if i == len(rungs) - 1:
            kept = out
        del out

    adopt_host_ownership(OUTPUT_DIR, comm=comm)

    # ---- setup figure on the finest rung ------------------------------------
    if figures_enabled():
        _display_submesh_figure(
            kept, summaries[-1]["diag"]["implant_rod_cad"]["tip_ball_radius_m"], comm
        )
    del kept

    # ---- values ----------------------------------------------------------------
    unset = summaries[0]
    ref_cad = unset["diag"]["implant_rod_cad"]
    ref_mass = unset["diag"]["cad_mass_by_group"]
    ratio_n = unset["n_cells"] / ROD_1G_RUNG_CELL_RECORD
    target = _phantom_minus_rod()
    low, high = TERMINAL_AREA_BAND
    bound = EXACT * ROD_LENGTH
    expected = {CONDUCTOR_CELL_TAG, AIR_CELL_TAG, PHANTOM_CELL_TAG}
    expected |= {PORT_LOWER + k for k in range(1, LEG_COUNT + 1)}
    expected |= {PORT_UPPER + k for k in range(1, LEG_COUNT + 1)}
    failures = []

    def check(ok, what):
        if not ok:
            failures.append(what)

    _say(f"\n[N] unset {unset['n_cells']} cells / ROD_1G_RUNG_CELL_RECORD "
         f"{ROD_1G_RUNG_CELL_RECORD} = {ratio_n:.6f} (band {CELL_COUNT_BAND})")
    check(abs(ratio_n - 1.0) < CELL_COUNT_BAND, f"(N) ratio {ratio_n:.6f}")

    for s in summaries:
        lab = s["label"]
        cad = s["diag"]["implant_rod_cad"]
        want = TIP_R4 if s["request"] is None else s["request"]
        a0 = (s["diag"]["cad_mass_by_group"] == ref_mass
              and cad["wall_area_m2"] == ref_cad["wall_area_m2"]
              and cad["wall_centroid_m"] == ref_cad["wall_centroid_m"]
              and cad["tip_size_m"] == want)
        a1 = abs(s["diag"]["cad_mass_by_group"]["phantom"] / target - 1.0)
        a3 = abs(cad["wall_area_m2"] / _rod_area() - 1.0)
        a4 = s["area402"] / _rod_area()
        xb, yb, zb = cad["wall_centroid_m"]
        a5 = (set(s["census"]) == expected
              and sum(s["census"].values()) == s["n_cells"]
              and s["n402"] > 0 and s["n402_int"] == 0
              and s["n402_ext"] == s["n402"] and s["n402_phantom"] == s["n402"])
        a6 = (abs(xb - ROD_CENTRE[0]) <= bound and abs(yb) <= bound and abs(zb) <= bound)
        d = s["dofs"]
        dok = d[2] == 2 * d["E"] + 2 * d["F"] and d[1] == d["E"]
        _say(f"[{lab}] (A0) {a0} (tip_size_m {cad['tip_size_m']!r}, wall area "
             f"{cad['wall_area_m2']!r}); (A1) {a1:.3e} (A3) {a3:.3e} (band "
             f"{CAD_ANALYTIC_BAND:.0e}); (A4) {a4:.6f} (band {TERMINAL_AREA_BAND}); (A5) "
             f"{a5} (tag 402 {s['n402']}); (A6) |x-xc| {abs(xb - ROD_CENTRE[0]):.3e} y "
             f"{yb:+.3e} z {zb:+.3e} (bound {bound:.1e}); (D) d2 {d[2]} = 2E+2F "
             f"{2 * d['E'] + 2 * d['F']}, d1 {d[1]} = E {d['E']}: {dok}")
        check(a0, f"(A0) {lab}")
        check(a1 < CAD_ANALYTIC_BAND, f"(A1) {lab} {a1:.3e}")
        check(a3 < CAD_ANALYTIC_BAND, f"(A3) {lab} {a3:.3e}")
        check(low < a4 < high, f"(A4) {lab} {a4:.6f}")
        check(a5, f"(A5) {lab}")
        check(a6, f"(A6) {lab}")
        check(dok, f"(D) {lab}")
        check(s["n_cells"] < RUNAWAY_CELL_GUARD, f"guard {lab}")

    # (G) strictly increasing N_ball
    balls = [s["n_ball"] for s in summaries]
    _say(f"[G] N_ball {balls}: strictly increasing {balls[0] < balls[1] < balls[2]}")
    check(balls[0] < balls[1] < balls[2], f"(G) {balls}")
    # (c1)
    c1_ok = (len(refusals) == 2 and gmsh_untouched
             and all(r is not None and re.search("implant_rod", r) for _, r, _ in refusals))
    _say(f"[c1] both refused with a message naming implant_rod, gmsh untouched: {c1_ok}")
    check(c1_ok, "(c1)")

    # ---- printed, not asserted ----------------------------------------------
    _say("\n[headline] N_ball / N_ball(unset) beside the cubic model (TIP_R4 / h)^3:")
    for s in summaries[1:]:
        h = s["request"]
        _say(f"  {s['label']}: {s['n_ball'] / unset['n_ball']:.4f} vs "
             f"{(TIP_R4 / h) ** 3:.3f}; cells / unset {s['n_cells'] / unset['n_cells']:.4f}")
    _say("[cost] rung | cells | N_ball | E | F | d1 | d2 | build s | ball edge / VIn | "
         f"near-T1 ({T1_RADIUS:.1e} m) cells, mean edge")
    for s in summaries:
        vin = s["diag"]["implant_rod_cad"]["tip_size_m"]
        _say(f"  {s['label']} | {s['n_cells']} | {s['n_ball']} | {s['dofs']['E']} | "
             f"{s['dofs']['F']} | {s['dofs'][1]} | {s['dofs'][2]} | {s['t']:.2f} | "
             f"{s['ball_edge'] / vin:.4f} | {s['n_t1']}, {s['t1_edge']:.6e} m")
    _say(f"[guard] r/8 cells {summaries[-1]['n_cells']} vs ROD_CELL_GUARD_W2_PRINTED "
         f"{ROD_CELL_GUARD_W2_PRINTED}: exceeds {summaries[-1]['n_cells'] > ROD_CELL_GUARD_W2_PRINTED}")
    _say("[c2] locality (printed, not asserted): air / phantom cells per rung:")
    for s in summaries:
        _say(f"  {s['label']}: air {s['census'].get(AIR_CELL_TAG)} phantom "
             f"{s['census'].get(PHANTOM_CELL_TAG)}; census {s['census']}")
    for f in files:
        _say(f"[paraview] wrote {f}")
    _say(f"[time] total {time.perf_counter() - started:.1f} s")

    # ---- asserts, identical on every rank (reduced values only) ---------------
    assert not failures, f"EX-102 assertions missed: {failures}"
    _say("\nAll identities and controls hold.")


if __name__ == "__main__":
    main()
