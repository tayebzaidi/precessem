"""`IMP-3` step 7 part 2 (§9 item 305, chain I15): the r/(4√2) rung's rod arm at degree 2.

I14's exact-rectangle hole-rod fixture on the 1 g rung (``phantom_resolution =
PHANTOM_RESOLUTION_1G_RUNG``, ``leg_sheet_strip_fraction = GATED_WIDTH_FRACTION``,
64 MHz, κ off) built with ``implant_rod = dict(ROD, tip_size=TIP_R4S2)``.  Only
the rod arm is solved; the off arm is I14 W2's record, read back from
``logs/imp3/ours_imp3g_d2off.json``.  One env knob ``FEM_EM_IMP3_STEP7`` in
{r4s2, r8, r8build, b4, b4build}; one test dispatches on it and skips otherwise.

Step 10a (§9 item 532, §10 I30a): ``b4`` solves the rod arm at
``tip_ball_radius = TIP_BALL_4MM`` with ``tip_size`` unset (the r/4 grain) —
the `xl` entry 43 window, never in a slot — and ``b4build`` builds that mesh
only.  Both dispatch inside the same disclosed ``ROD_CELL_GUARD_XL`` patch as
r8 / r8build (entry 43's licence).  Asserted there: (Kb4) radius and grain
bitwise and cells == 351 478; (T1in) item 464's near-T1-inside-the-ball
predicate; (G1x); (N1); (G3); (Rb), (Rb8), (Rb-r8) from the stored JSONs.
``b4`` writes its JSON before the readout; ``b4build`` dry-runs the readout
with the r/8 JSON as a labelled stand-in.

Part 3a (§9 item 327): ``r8`` solves the rod arm at ``tip_size = TIP_R8`` (the
`xl` entry 33 window; never in a slot), ``r8build`` builds that mesh only.  On
those two paths only, ``test_imp3_tip_ratio_1g.ROD_CELL_GUARD_W2`` (the name
``_arm_d``'s (G1) resolves at call time) is monkeypatched to
``ROD_CELL_GUARD_XL`` — the disclosed raise `xl` entry 33's price paragraph
licenses, printed, undone at teardown.  Asserted there: (Kc8) ``tip_size_m`` ==
TIP_R8 bitwise and cells > item 305's 276 088; (G1x) the patch is in place and
cells < it; (N1) cells > the unpatched 320 000; (G3); (Rb8) η_rod^{r/(4√2)}(T1)
P1 recomputed from ``ours_imp3h_r4s2.json`` string-equal to "3.556262", sizes
equal.  (N2) every call checks the unpatched guard reads 320 000 before
dispatch, so the r4s2 path never runs inside the patch.

Asserted (bands imported by name through I14's helpers, never restated):
    (G1) rod cells < ROD_CELL_GUARD_W2 (inside ``_arm_d``, before the solve);
    (G2) N_i areas, (i), (ii), (PP), (M), (E), (V), (D) via ``_identities``;
    (G3) tag-402 facets > 0 and all exterior;
    (Kc) ``tip_size_m`` == TIP_R4S2 bitwise and rod cells > ROD_1G_RUNG_CELL_RECORD;
    (Rb) η_rod^{r/4}(T1), P1, recomputed from I14's JSONs, string-equal at
         6 d.p. to "3.482328"; the JSONs' cells / dofs equal I14's as integers.

Printed, never gated: η_rod at r/(4√2) beside r/4 at six points, the first
ladder step at T1 (no step ratio: that needs r/8), (C-h), (C2), mesh identity
against GEO-36 W2, the near-T1 and T1-cell edges, N_ball, Q per tip and drive,
the closed form, prices.  No convergence and no accuracy claim (`ANS-10`).
"""

from __future__ import annotations

import gc
import json
import os
import time
from pathlib import Path

import numpy as np
import pytest
from mpi4py import MPI

import tests.validation.test_imp3_tip_ratio_1g as _ratio
from fem_em_solver.post.evaluation import point_owner_cells
from fem_em_solver.utils.analytical import pec_prolate_spheroid_field
from tests.complex_mode import complex_only
from tests.mesh.test_birdcage_implant_rod import ROD, ROD_CENTRE, ROD_RADIUS, _build_rod
from tests.mesh.test_birdcage_implant_rod_ball import TIP_BALL_4MM, TOP_CAP, _t1_geometry
from tests.mesh.test_birdcage_implant_rod_tip import (
    FIXTURE_KW,
    ROD_1G_RUNG_CELL_RECORD,
    T1,
    T1_RADIUS,
    TIP_R4,
    TIP_R4S2,
    TIP_R8,
    _ball_and_t1,
    _cell_geometry,
    _mean_edge,
)
from tests.validation.test_imp3_palace_fields import PNAMES, PTS
from tests.validation.test_imp3_rod_birdcage import (
    DRIVES,
    PORT_IDS,
    ROD_CELL_GUARD_W2,
    TIPS,
    _outer_box_tags,
    _peak_rss_gib,
)
from tests.validation.test_imp3_tip_ratio_1g import (
    _arm_d,
    _cplx,
    _identities,
    _ser,
)

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "logs" / "imp3"
KNOB = os.environ.get("FEM_EM_IMP3_STEP7")
KNOBS = ("r4s2", "r8", "r8build", "b4", "b4build")
TAG = "[IMP-3 7]"

# I14 records (`20261006T173320Z_IMP-3-item284-W4-d1-readout.log:3894, 3905`).
I14_ETA_O2_T1_P1 = "3.482328"
I14_Q_O2_P1_TOP = 4.463326
# I14's p-move and remesh move at T1, P1 (`IMP-3.md:950`), printed for scale.
I14_PMOVE_T1 = -5.995e-02
I14_REMESH_MOVE_T1 = -3.122e-02
# I14 W2 / W3 JSON sizes (off / rod at r/4), asserted under (Rb).
I14_OFF_CELLS, I14_OFF_DOFS = 201_072, 1_284_856
I14_ROD_CELLS, I14_ROD_DOFS = 252_118, 1_620_748
# GEO-36 step 4 W2 (`20261007T170700Z_GEO-36-step4-W2-rungs.log:5612-5613, 7516`);
# printed comparands, not bands.
GEO36_R4S2_CELLS = 276_088
GEO36_R4S2_D2_DOFS = 1_775_954
GEO36_R4S2_N402 = 6_256
GEO36_R4S2_NEAR_T1 = (19, 6.422237e-04)
GEO36_UNSET_NEAR_T1_EDGE = 6.983224e-04
GEO36_R4S2_N_BALL, GEO36_R4S2_BALL_EDGE = 32_826, 2.433e-4
# GEO-36 price row (`GEO-36.md:233`).
GEO36_PRICE = "sweep 318.8-322.9 s, window 491.3-495.4 s, 57.04-57.77 GiB"
PRED_CH_FLOOR = 1.0e-3
PRED_C2_FLOOR = 3.0

# --- part 3a (§9 item 327): the r/8 rung --------------------------------------
# `xl` resource guard, NOT a physics band: the 128 G-container sizing logic of
# ROD_CELL_GUARD_W2 re-applied to the `xl` window.  Derived from item 305 W1's
# measured 53.85 GiB summed peak RSS at -n 8 for 276 088 cells
# (`20261008T033458Z_IMP-3-item305-W1-d2rod-r4s2.log`): at the guard (x1.449
# cells) ≈ 78-89.5 GiB at e = 1.0 / 1.37; r/8 (342 072, x1.239) ≈ 66.7-72.2 GiB,
# at 0.855 of the guard; r/16 (≈ 1 M cells, `GEO-36.md:236`) still trips it.
# Applied by monkeypatch on the r8 / r8build paths only (entry 33's licence).
ROD_CELL_GUARD_XL = 400_000
# GEO-36 step 4 W2, r/8 rung (`20261007T170700Z_GEO-36-step4-W2-rungs.log`):
# cells and d2 dofs :7500 (d2 also :7514); tag 402 facets, N_ball / mean ball
# edge and near-T1 (0.5 mm) count / mean edge :7501.  Printed comparands.
GEO36_R8_CELLS = 342_072
GEO36_R8_D2_DOFS = 2_200_500
GEO36_R8_N402 = 8_790
GEO36_R8_NEAR_T1 = (26, 5.081837e-04)
GEO36_R8_N_BALL, GEO36_R8_BALL_EDGE = 93_132, 1.722398e-04
# Item 305 W1 records (`20261008T033458Z_IMP-3-item305-W1-d2rod-r4s2.log`,
# `IMP-3.md:1111-1114`): (Rb8) string and sizes asserted, the rest printed.
ITEM305_ETA_T1_P1 = "3.556262"
ITEM305_ROD_CELLS, ITEM305_ROD_DOFS = 276_088, 1_775_954
ITEM305_PRICE = "sweep 267.03 s, body 383.2 s, 53.85 GiB summed peak RSS at -n 8"
# Near-T1 edge ratio h(r/(4√2)) / h(r/8) from GEO-36 W2 (6.422237e-04 / 5.081837e-04).
PRED_EDGE_RATIO_R4S2_R8 = 1.264
# `xl` entry 33's readout threshold on |Δ₂| / |Δ₁| (printed; the weekly reads it).
LADDER_RATIO_THRESHOLD = 0.7

# --- step 10a (§9 item 532, §10 I30a): the b4 knob ------------------------------
assert TIP_BALL_4MM == 4.0 * ROD_RADIUS  # import-time identity
# Item 464 W2, the same 4 mm / r/4-grain build at -n 2
# (`20261015T123700Z_IMP-3-I17-W2-n2.log`): cells :5648 (asserted; identical at
# -n 1, `…I17-W1-n1.log`, `IMP-3.md:1703`), d2 dofs / E / F :5648, :5661;
# tag-402 facets, N_ball / mean ball edge, near-T1 count / mean edge :5649;
# near-T1 vertex max and the T1 cell (original index, mean edge, max vertex
# distance) :5650; build seconds :5648.  Printed comparands except the cells.
I17_B4_CELLS = 351_478
I17_B4_D2_DOFS, I17_B4_E, I17_B4_F = 2_251_660, 418_532, 707_298
I17_B4_N402 = 5_366
I17_B4_N_BALL, I17_B4_BALL_EDGE = 111_415, 3.489008e-04
I17_B4_NEAR_T1 = (130, 3.383011e-04)
I17_B4_T1_VMAX = 3.151970e-03
I17_B4_T1_CELL = (227_010, 3.282035e-04, 2.768817e-03)
I17_B4_BUILD_S = 77.98
# The 2 mm mesh's near-T1 vertex max (`…I17-W2-n2.log:3755`), the (T1in-blind) control.
I17_2MM_T1_VMAX = 3.068100e-03
# r/8 rung η_rod(T1), P1 (`20261014T070008Z_IMP-3-step7-I15-r8.log:2643`); sizes
# asserted against GEO36_R8_CELLS / GEO36_R8_D2_DOFS (the same integers).
R8_ETA_T1_P1 = "3.530461"
ENTRY33_PRICE = (
    "build 93.47 s, sweep 278.25 s, SARs 24.57 s, body 402.2 s, 82.47 GiB "
    "summed at -n 16 (`…I15-r8.log:2663`); window 511 s, memory.peak "
    "86 486 282 240 B = 80.55 GiB (`:3773`, `:3778`)"
)
ENTRY43_PRICE_PRED = "500–700 s, 75–100 GiB"
# Entry 43's pre-registered decision interval: printed IN / OUT, never asserted,
# never a band.
ENTRY43_READOUT_INTERVAL = (3.48, 3.56)


def _say(msg):
    if MPI.COMM_WORLD.rank == 0:
        print(msg, flush=True)


def _eta(e_rod, e_off):
    """η(p, port) = ‖E_rod(p)‖ / ‖E_off(p)‖, Euclidean norm over components."""
    return np.linalg.norm(e_rod, axis=1) / np.linalg.norm(e_off, axis=1)


def _r4s2():
    comm = MPI.COMM_WORLD
    t0 = time.perf_counter()
    off = json.loads((OUT / "ours_imp3g_d2off.json").read_text())
    rod4 = json.loads((OUT / "ours_imp3g_d2rod.json").read_text())
    target = OUT / "ours_imp3h_r4s2.json"

    built = _build_rod(implant_rod=dict(ROD, tip_size=TIP_R4S2), **FIXTURE_KW)
    cad = built[3]["implant_rod_cad"]
    _say(
        f"\n{TAG} build: tip_size_m {cad['tip_size_m']!r} (request {TIP_R4S2!r}), "
        f"tip_ball_radius_m {cad['tip_ball_radius_m']!r}, build {built[4]:.2f} s"
    )
    a = _arm_d(built, "degree 2 rod r/(4√2)", True, 2)
    fails = _identities(a, True)
    msh = a["case"]["mesh"]

    # (G3) tag 402 > 0 and all exterior, as I14's `_census`.
    _, census = _outer_box_tags(msh, built[2], comm)
    _say(
        f"{TAG} (G3) 402 facets {census['wall']} (exterior {census['wall_exterior']}); "
        f"exterior {census['exterior']}"
    )
    if not (census["wall"] > 0 and census["wall_exterior"] == census["wall"]):
        fails.append(("G3", census))

    # (Kc) the knob reached gmsh.
    kc_tip = cad["tip_size_m"] == TIP_R4S2
    kc_cells = a["cells"] > ROD_1G_RUNG_CELL_RECORD
    _say(
        f"{TAG} (Kc) tip_size_m == TIP_R4S2 bitwise: {kc_tip}; cells {a['cells']} > "
        f"ROD_1G_RUNG_CELL_RECORD {ROD_1G_RUNG_CELL_RECORD}: {kc_cells} "
        f"(ratio {a['cells'] / ROD_1G_RUNG_CELL_RECORD:.4f})"
    )
    if not (kc_tip and kc_cells):
        fails.append(("Kc", cad["tip_size_m"], a["cells"]))

    # Local mesh at T1 (collective; every rank, outside any rank branch).
    n_ball, ball_sum, n_t1, t1_sum = _ball_and_t1(msh, built[1], cad["tip_ball_radius_m"], comm)
    ball_edge = ball_sum / n_ball if n_ball else float("nan")
    t1_edge = t1_sum / n_t1 if n_t1 else float("nan")
    owner_orig, owner_local = point_owner_cells(msh, T1[None, :], comm)
    local = 0.0
    if owner_local[0] >= 0:
        local = float(_mean_edge(_cell_geometry(msh)[[int(owner_local[0])]])[0])
    t1_cell_edge = comm.allreduce(local, op=MPI.MAX)
    rss = _peak_rss_gib()

    # Mesh identity against GEO-36 W2 (printed, predicted).
    _say(
        f"{TAG} mesh identity vs GEO-36 W2: cells {a['cells']} / {GEO36_R4S2_CELLS} "
        f"({a['cells'] == GEO36_R4S2_CELLS}); 402 facets {census['wall']} / {GEO36_R4S2_N402} "
        f"({census['wall'] == GEO36_R4S2_N402}); d2 dofs {a['ndofs']} / {GEO36_R4S2_D2_DOFS} "
        f"({a['ndofs'] == GEO36_R4S2_D2_DOFS}); near-T1 ({T1_RADIUS:.1e} m) {n_t1} cells, mean "
        f"edge {t1_edge:.6e} m / {GEO36_R4S2_NEAR_T1} "
        f"({n_t1 == GEO36_R4S2_NEAR_T1[0] and f'{t1_edge:.6e}' == f'{GEO36_R4S2_NEAR_T1[1]:.6e}'})"
    )
    _say(
        f"{TAG} near-T1 edge r/4 {GEO36_UNSET_NEAR_T1_EDGE:.6e} -> r/(4√2) {t1_edge:.6e} m "
        f"({t1_edge / GEO36_UNSET_NEAR_T1_EDGE - 1:+.4%}); local edge ratio h(r/4)/h(r/(4√2)) "
        f"{GEO36_UNSET_NEAR_T1_EDGE / t1_edge:.4f} beside the requested VIn ratio "
        f"{np.sqrt(2.0):.4f} (VIn, not the step at T1)"
    )
    _say(
        f"{TAG} T1 cell: original index {int(owner_orig[0])}, mean edge {t1_cell_edge:.6e} m "
        f"(first record at this rung; T1 {T1.tolist()} vs PTS[T1] "
        f"{PTS[PNAMES.index('T1')].tolist()})"
    )
    _say(
        f"{TAG} N_ball {n_ball}, mean ball edge {ball_edge:.6e} m beside GEO-36 W2 "
        f"{GEO36_R4S2_N_BALL} / {GEO36_R4S2_BALL_EDGE:.3e} m"
    )

    # (Rb) the off-arm reuse is wired.
    iT1, iF1, iF2 = PNAMES.index("T1"), PNAMES.index("F1"), PNAMES.index("F2")
    e_off = _cplx(off["E"])
    eta4 = _eta(_cplx(rod4["E"]), e_off)
    rb_str = f"{eta4[iT1, 0]:.6f}"
    rb_sizes = (off["cells"], off["ndofs"], rod4["cells"], rod4["ndofs"]) == (
        I14_OFF_CELLS,
        I14_OFF_DOFS,
        I14_ROD_CELLS,
        I14_ROD_DOFS,
    )
    rb = rb_str == I14_ETA_O2_T1_P1 and rb_sizes
    _say(
        f"{TAG} (Rb) η_rod^(r/4)(T1) P1 from I14 JSONs {rb_str} vs {I14_ETA_O2_T1_P1}: "
        f"{rb_str == I14_ETA_O2_T1_P1}; sizes off {off['cells']}/{off['ndofs']} rod "
        f"{rod4['cells']}/{rod4['ndofs']} vs {I14_OFF_CELLS}/{I14_OFF_DOFS} "
        f"{I14_ROD_CELLS}/{I14_ROD_DOFS}: {rb_sizes}"
    )
    if not rb:
        fails.append(("Rb", rb_str, off["cells"], off["ndofs"], rod4["cells"], rod4["ndofs"]))

    if rb and kc_tip and kc_cells:
        eta = _eta(a["E"], e_off)
        for j, pid in enumerate(PORT_IDS):
            _say(
                f"{TAG} η_rod {pid} r/4 / r/(4√2) / move: "
                + "  ".join(
                    f"{n} {eta4[k, j]:.6f}/{eta[k, j]:.6f}/{eta[k, j] / eta4[k, j] - 1:+.6e}"
                    for k, n in enumerate(PNAMES)
                )
            )
        move = eta[iT1, 0] / eta4[iT1, 0] - 1
        _say(
            f"{TAG} first ladder step at T1, P1: η_rod(T1) r/4 {eta4[iT1, 0]:.6f} -> "
            f"r/(4√2) {eta[iT1, 0]:.6f}, move {move:+.6e} ({move:+.4%}); beside I14 p-move "
            f"{I14_PMOVE_T1:+.3%} and remesh move {I14_REMESH_MOVE_T1:+.3%} (scale only; no "
            f"step ratio: it needs r/8)"
        )
        for n in ("F1", "F2"):
            k = PNAMES.index(n)
            _say(f"{TAG} remesh floor {n} P1 move {eta[k, 0] / eta4[k, 0] - 1:+.6e}")
        _say(
            f"{TAG} (C-h) |move| {abs(move):.6e} (predicted >= {PRED_CH_FLOOR:g}, ceiling ≈ 1): "
            f"{'IN' if abs(move) >= PRED_CH_FLOOR else 'OUT'}"
        )
        for lab, e in (("r/4", eta4), ("r/(4√2)", eta)):
            den = max(abs(e[iF1, 0] - 1), abs(e[iF2, 0] - 1))
            c2 = abs(e[iT1, 0] - 1) / den
            _say(
                f"{TAG} (C2) {lab}: |η(T1)-1| / max_F |η(F)-1| = {abs(e[iT1, 0] - 1):.6f} / "
                f"{den:.6f} = {c2:.3f} (predicted >= {PRED_C2_FLOOR:g}: "
                f"{'IN' if c2 >= PRED_C2_FLOOR else 'OUT'})"
            )
        for d in DRIVES:
            for tip in TIPS:
                s_off = off["sar"][f"{tip}|{d}"]
                q4 = rod4["sar"][f"{tip}|{d}"] / s_off
                qh = a["sar"][(tip, d)] / s_off
                _say(
                    f"{TAG} Q {d:>3s} {tip:>6s}: r/4 {q4:.6f}  r/(4√2) {qh:.6f}  move "
                    f"{qh / q4 - 1:+.6e}"
                )
        _say(
            f"{TAG} Q r/4 P1 top of record {I14_Q_O2_P1_TOP}; rod-arm 1 g masses r/4 "
            f"{rod4['mass']} r/(4√2) {a['mass']}"
        )
        cf = np.linalg.norm(
            pec_prolate_spheroid_field(
                PTS - np.array(ROD_CENTRE), semi_major=0.016, semi_minor=0.001, e0=1.0
            ),
            axis=1,
        )
        _say(
            f"{TAG} closed form |E_cf|/E0 vs η_rod(p)/η_rod(F1) (r/(4√2), P1), ratio: "
            + "  ".join(
                f"{PNAMES[k]} {cf[k]:.4f}/{eta[k, 0] / eta[iF1, 0]:.4f}="
                f"{cf[k] / (eta[k, 0] / eta[iF1, 0]):.3f}"
                for k in (PNAMES.index(n) for n in ("T1", "T2", "T3", "S1"))
            )
        )
    else:
        _say(f"{TAG} no η / Q reading printed: (Rb) {rb}, (Kc) {kc_tip and kc_cells}")

    _say(
        f"{TAG} price: cells {a['cells']}, E {a['edges']}, F {a['faces']}, dofs {a['ndofs']}, "
        f"build {built[4]:.2f} s, sweep {a['t_sweep']:.2f} s, SARs {a['t_sar']:.2f} s, arm "
        f"{a['t_arm']:.1f} s, window body {time.perf_counter() - t0:.1f} s, summed peak RSS "
        f"{rss:.2f} GiB at -n {comm.size}; beside GEO-36 price row {GEO36_PRICE}"
    )

    payload = _ser(a)
    payload.update(
        {
            "rss_gib": rss,
            "ranks": comm.size,
            "tip_size_m": cad["tip_size_m"],
            "near_t1": [n_t1, t1_edge],
            "t1_cell": [int(owner_orig[0]), t1_cell_edge],
            "n_ball": [n_ball, ball_edge],
        }
    )
    if comm.rank == 0:
        with open(target, "x") as fh:
            json.dump(payload, fh)
    comm.Barrier()
    _say(f"{TAG} anchors failed (all): {fails if fails else 'none'}")
    del a, built, msh
    gc.collect()
    assert not fails, fails


def _r8(build_only):
    """Part 3a: the r/8 rung.  ``build_only`` (r8build) never calls ``_arm_d``."""
    comm = MPI.COMM_WORLD
    t0 = time.perf_counter()
    target = OUT / "ours_imp3h_r8.json"
    if not build_only:
        # mode "x" below: fail before the build, not after the sweep.
        assert not target.exists(), f"{target} exists; the r8 window writes it with mode 'x'"
    off = json.loads((OUT / "ours_imp3g_d2off.json").read_text())
    rod4 = json.loads((OUT / "ours_imp3g_d2rod.json").read_text())
    r4s2 = json.loads((OUT / "ours_imp3h_r4s2.json").read_text())
    fails = []

    # (Rb8) the r/(4√2) reuse is wired.
    iT1 = PNAMES.index("T1")
    e_off = _cplx(off["E"])
    eta4 = _eta(_cplx(rod4["E"]), e_off)
    eta4s2 = _eta(_cplx(r4s2["E"]), e_off)
    rb_str = f"{eta4s2[iT1, 0]:.6f}"
    rb_sizes = (r4s2["cells"], r4s2["ndofs"], off["cells"], off["ndofs"]) == (
        ITEM305_ROD_CELLS,
        ITEM305_ROD_DOFS,
        I14_OFF_CELLS,
        I14_OFF_DOFS,
    )
    rb = rb_str == ITEM305_ETA_T1_P1 and rb_sizes
    _say(
        f"\n{TAG} (Rb8) η_rod^(r/(4√2))(T1) P1 from ours_imp3h_r4s2.json + off {rb_str} vs "
        f"{ITEM305_ETA_T1_P1}: {rb_str == ITEM305_ETA_T1_P1}; sizes r4s2 {r4s2['cells']}/"
        f"{r4s2['ndofs']} off {off['cells']}/{off['ndofs']} vs {ITEM305_ROD_CELLS}/"
        f"{ITEM305_ROD_DOFS} {I14_OFF_CELLS}/{I14_OFF_DOFS}: {rb_sizes}"
    )
    if not rb:
        fails.append(("Rb8", rb_str, r4s2["cells"], r4s2["ndofs"], off["cells"], off["ndofs"]))

    built = _build_rod(implant_rod=dict(ROD, tip_size=TIP_R8), **FIXTURE_KW)
    cad = built[3]["implant_rod_cad"]
    msh = built[0]
    cells = int(msh.topology.index_map(msh.topology.dim).size_global)
    _say(
        f"{TAG} build r/8: tip_size_m {cad['tip_size_m']!r} (request {TIP_R8!r}), "
        f"tip_ball_radius_m {cad['tip_ball_radius_m']!r}, {cells} cells, build {built[4]:.2f} s"
    )

    # (Kc8) the knob reached gmsh.
    kc_tip = cad["tip_size_m"] == TIP_R8
    kc_cells = cells > ITEM305_ROD_CELLS
    _say(
        f"{TAG} (Kc8) tip_size_m == TIP_R8 bitwise: {kc_tip}; cells {cells} > item 305 "
        f"{ITEM305_ROD_CELLS}: {kc_cells} (ratio {cells / ITEM305_ROD_CELLS:.4f})"
    )
    if not (kc_tip and kc_cells):
        fails.append(("Kc8", cad["tip_size_m"], cells))
    # (G1x) inside the patch; (N1) the unpatched (G1) would trip.
    patched = _ratio.ROD_CELL_GUARD_W2
    g1x = patched == ROD_CELL_GUARD_XL and cells < ROD_CELL_GUARD_XL
    n1 = cells > ROD_CELL_GUARD_W2
    _say(
        f"{TAG} DISCLOSURE: test_imp3_tip_ratio_1g.ROD_CELL_GUARD_W2 monkeypatched to "
        f"{patched} (ROD_CELL_GUARD_XL, an `xl` resource guard, not a band; original "
        f"{ROD_CELL_GUARD_W2}) on this path only; measured cells {cells}"
    )
    _say(
        f"{TAG} (G1x) patched == ROD_CELL_GUARD_XL {ROD_CELL_GUARD_XL}: "
        f"{patched == ROD_CELL_GUARD_XL}; cells {cells} < it: {cells < ROD_CELL_GUARD_XL} "
        f"({cells / ROD_CELL_GUARD_XL:.3f} of the guard)"
    )
    _say(
        f"{TAG} (N1) cells {cells} > unpatched ROD_CELL_GUARD_W2 {ROD_CELL_GUARD_W2}: {n1} "
        f"(ratio {cells / ROD_CELL_GUARD_W2:.3f}; predicted 1.069)"
    )
    if not g1x:
        fails.append(("G1x", patched, cells))
    if not n1:
        fails.append(("N1", cells, ROD_CELL_GUARD_W2))
    # A wrong mesh never reaches the `xl` solve; an (Rb8) red alone still
    # measures the identities (item 305's (Rb) rule) and prints no η / Q.
    mesh_fails = [f for f in fails if f[0] != "Rb8"]
    if mesh_fails and not build_only:
        _say(f"{TAG} mesh anchors failed before the solve: {mesh_fails}")
        assert not fails, fails

    a = None
    if not build_only:
        a = _arm_d(built, "degree 2 rod r/8", True, 2)
        fails += _identities(a, True)
        msh = a["case"]["mesh"]

    # (G3) tag 402 > 0 and all exterior.
    _, census = _outer_box_tags(msh, built[2], comm)
    _say(
        f"{TAG} (G3) 402 facets {census['wall']} (exterior {census['wall_exterior']}); "
        f"exterior {census['exterior']}"
    )
    if not (census["wall"] > 0 and census["wall_exterior"] == census["wall"]):
        fails.append(("G3", census))

    # Local mesh at T1 (collective; every rank, outside any rank branch).
    n_ball, ball_sum, n_t1, t1_sum = _ball_and_t1(msh, built[1], cad["tip_ball_radius_m"], comm)
    ball_edge = ball_sum / n_ball if n_ball else float("nan")
    t1_edge = t1_sum / n_t1 if n_t1 else float("nan")
    owner_orig, owner_local = point_owner_cells(msh, T1[None, :], comm)
    local = 0.0
    if owner_local[0] >= 0:
        local = float(_mean_edge(_cell_geometry(msh)[[int(owner_local[0])]])[0])
    t1_cell_edge = comm.allreduce(local, op=MPI.MAX)
    if build_only:
        n_e, n_f = _ratio._edge_face_counts(msh)
        ndofs = 2 * n_e + 2 * n_f
    else:
        n_e, n_f, ndofs = a["edges"], a["faces"], a["ndofs"]
    rss = _peak_rss_gib()

    _say(
        f"{TAG} mesh identity vs GEO-36 W2 r/8: cells {cells} / {GEO36_R8_CELLS} "
        f"({cells == GEO36_R8_CELLS}); 402 facets {census['wall']} / {GEO36_R8_N402} "
        f"({census['wall'] == GEO36_R8_N402}); d2 dofs {ndofs} (E {n_e}, F {n_f}"
        f"{', 2E + 2F' if build_only else ''}) / {GEO36_R8_D2_DOFS} "
        f"({ndofs == GEO36_R8_D2_DOFS}); near-T1 ({T1_RADIUS:.1e} m) {n_t1} cells, mean "
        f"edge {t1_edge:.6e} m / {GEO36_R8_NEAR_T1} "
        f"({n_t1 == GEO36_R8_NEAR_T1[0] and f'{t1_edge:.6e}' == f'{GEO36_R8_NEAR_T1[1]:.6e}'})"
    )
    h4s2 = float(r4s2["near_t1"][1])
    _say(
        f"{TAG} near-T1 edge r/(4√2) {h4s2:.6e} -> r/8 {t1_edge:.6e} m "
        f"({t1_edge / h4s2 - 1:+.4%}); local edge ratio h(r/(4√2))/h(r/8) "
        f"{h4s2 / t1_edge:.4f} beside GEO-36's {PRED_EDGE_RATIO_R4S2_R8} and the requested VIn "
        f"ratio {np.sqrt(2.0):.4f} (VIn, not the step at T1)"
    )
    _say(
        f"{TAG} T1 cell: original index {int(owner_orig[0])}, mean edge {t1_cell_edge:.6e} m "
        f"(r/(4√2) {r4s2['t1_cell']})"
    )
    _say(
        f"{TAG} N_ball {n_ball}, mean ball edge {ball_edge:.6e} m beside GEO-36 W2 "
        f"{GEO36_R8_N_BALL} / {GEO36_R8_BALL_EDGE:.6e} m"
    )

    if a is not None and rb and kc_tip and kc_cells:
        eta = _eta(a["E"], e_off)
        for j, pid in enumerate(PORT_IDS):
            _say(
                f"{TAG} η_rod {pid} r/4 / r/(4√2) / r/8: "
                + "  ".join(
                    f"{n} {eta4[k, j]:.6f}/{eta4s2[k, j]:.6f}/{eta[k, j]:.6f}"
                    for k, n in enumerate(PNAMES)
                )
            )
        d1 = eta4s2[iT1, 0] - eta4[iT1, 0]
        d2 = eta[iT1, 0] - eta4s2[iT1, 0]
        ratio = abs(d2) / abs(d1)
        _say(
            f"{TAG} ladder at T1, P1: η_rod(T1) r/4 {eta4[iT1, 0]:.6f} -> r/(4√2) "
            f"{eta4s2[iT1, 0]:.6f} -> r/8 {eta[iT1, 0]:.6f}; Δ₁ {d1:+.6f}, Δ₂ {d2:+.6f}; signs "
            f"{'same' if np.sign(d1) == np.sign(d2) else 'opposite'}; |Δ₂|/|Δ₁| {ratio:.4f} "
            f"beside entry 33's {LADDER_RATIO_THRESHOLD} (0.7|Δ₁| = "
            f"{LADDER_RATIO_THRESHOLD * abs(d1):.6f}); printed, never asserted; unequal local "
            f"steps, so no observed order"
        )
        for n in ("F1", "F2"):
            k = PNAMES.index(n)
            _say(
                f"{TAG} remesh floor {n} P1 move r/(4√2) -> r/8 "
                f"{eta[k, 0] / eta4s2[k, 0] - 1:+.6e}"
            )
        for d in DRIVES:
            for tip in TIPS:
                s_off = off["sar"][f"{tip}|{d}"]
                q4 = rod4["sar"][f"{tip}|{d}"] / s_off
                q4s2 = r4s2["sar"][f"{tip}|{d}"] / s_off
                q8 = a["sar"][(tip, d)] / s_off
                _say(
                    f"{TAG} Q {d:>3s} {tip:>6s}: r/4 {q4:.6f}  r/(4√2) {q4s2:.6f}  r/8 " f"{q8:.6f}"
                )
        _say(
            f"{TAG} rod-arm 1 g masses r/4 {rod4['mass']} r/(4√2) {r4s2['mass']} r/8 "
            f"{a['mass']}"
        )
    elif a is not None:
        _say(f"{TAG} no η / Q reading printed: (Rb8) {rb}, (Kc8) {kc_tip and kc_cells}")

    if a is not None:
        _say(
            f"{TAG} price: cells {a['cells']}, E {a['edges']}, F {a['faces']}, dofs "
            f"{a['ndofs']}, build {built[4]:.2f} s, sweep {a['t_sweep']:.2f} s, SARs "
            f"{a['t_sar']:.2f} s, arm {a['t_arm']:.1f} s, window body "
            f"{time.perf_counter() - t0:.1f} s, summed peak RSS {rss:.2f} GiB at -n "
            f"{comm.size}; beside item 305 {ITEM305_PRICE}"
        )
        payload = _ser(a)
        payload.update(
            {
                "rss_gib": rss,
                "ranks": comm.size,
                "tip_size_m": cad["tip_size_m"],
                "near_t1": [n_t1, t1_edge],
                "t1_cell": [int(owner_orig[0]), t1_cell_edge],
                "n_ball": [n_ball, ball_edge],
                "guard_patched": patched,
            }
        )
        if comm.rank == 0:
            with open(target, "x") as fh:
                json.dump(payload, fh)
        comm.Barrier()
    else:
        _say(
            f"{TAG} r8build price: build {built[4]:.2f} s, body "
            f"{time.perf_counter() - t0:.1f} s, summed peak RSS {rss:.2f} GiB at -n "
            f"{comm.size}; no solve, no JSON"
        )
    _say(f"{TAG} anchors failed (all): {fails if fails else 'none'}")
    del a, built, msh
    gc.collect()
    assert not fails, fails


def _b4_readout(eta_new, sar_new, mass_new, label, refs, h_new, prefix=""):
    """Item 532 readout (printed, never asserted).  ``sar_new`` keyed "tip|drive"."""
    off, rod4, r4s2, r8 = refs
    iT1 = PNAMES.index("T1")
    e_off = _cplx(off["E"])
    eta4 = _eta(_cplx(rod4["E"]), e_off)
    eta4s2 = _eta(_cplx(r4s2["E"]), e_off)
    eta8 = _eta(_cplx(r8["E"]), e_off)
    p = f"{TAG} {prefix}"
    for j, pid in enumerate(PORT_IDS):
        _say(
            f"{p}η_rod {pid} r/4 / r/(4√2) / r/8 / {label}: "
            + "  ".join(
                f"{n} {eta4[k, j]:.6f}/{eta4s2[k, j]:.6f}/{eta8[k, j]:.6f}/{eta_new[k, j]:.6f}"
                for k, n in enumerate(PNAMES)
            )
        )
    v = eta_new[iT1, 0]
    lo, hi = ENTRY43_READOUT_INTERVAL
    move = v / eta4[iT1, 0] - 1
    _say(
        f"{p}η_rod(T1) P1: r/4 {eta4[iT1, 0]:.6f} (record {I14_ETA_O2_T1_P1}), r/(4√2) "
        f"{eta4s2[iT1, 0]:.6f} (record {ITEM305_ETA_T1_P1}), r/8 {eta8[iT1, 0]:.6f} (record "
        f"{R8_ETA_T1_P1}), {label} {v:.6f}; move vs r/4 {move:+.6e} ({move:+.4%}); entry 43 "
        f"decision interval [{lo}, {hi}]: {'IN' if lo <= v <= hi else 'OUT'} (printed, never "
        f"asserted, not a band)"
    )
    for n in ("F1", "F2"):
        k = PNAMES.index(n)
        _say(f"{p}remesh floor {n} P1 move r/4 -> {label} {eta_new[k, 0] / eta4[k, 0] - 1:+.6e}")
    for d in DRIVES:
        for tip in TIPS:
            key = f"{tip}|{d}"
            s_off = off["sar"][key]
            _say(
                f"{p}Q {d:>3s} {tip:>6s}: r/4 {rod4['sar'][key] / s_off:.6f}  r/(4√2) "
                f"{r4s2['sar'][key] / s_off:.6f}  r/8 {r8['sar'][key] / s_off:.6f}  {label} "
                f"{sar_new[key] / s_off:.6f}"
            )
    _say(
        f"{p}rod-arm 1 g masses r/4 {rod4['mass']} r/(4√2) {r4s2['mass']} r/8 {r8['mass']} "
        f"{label} {mass_new}"
    )
    _say(
        f"{p}near-T1 edge ratio h(r/4 unset) {GEO36_UNSET_NEAR_T1_EDGE:.6e} / h({label} mesh) "
        f"{h_new:.6e} = {GEO36_UNSET_NEAR_T1_EDGE / h_new:.4f}"
    )


def _b4(build_only):
    """Step 10a: the 4 mm tip ball at the r/4 grain.  ``b4build`` never calls ``_arm_d``."""
    comm = MPI.COMM_WORLD
    t0 = time.perf_counter()
    target = OUT / "ours_imp3h_b4.json"
    if not build_only:
        # mode "x" below: fail before the build, not after the sweep.
        assert not target.exists(), f"{target} exists; the b4 window writes it with mode 'x'"
    off = json.loads((OUT / "ours_imp3g_d2off.json").read_text())
    rod4 = json.loads((OUT / "ours_imp3g_d2rod.json").read_text())
    r4s2 = json.loads((OUT / "ours_imp3h_r4s2.json").read_text())
    r8 = json.loads((OUT / "ours_imp3h_r8.json").read_text())
    refs = (off, rod4, r4s2, r8)
    fails = []

    # (Rb), (Rb8), (Rb-r8): the stored-reading paths, on both knobs.
    iT1 = PNAMES.index("T1")
    e_off = _cplx(off["E"])
    for lab, js, want, sizes in (
        ("Rb", rod4, I14_ETA_O2_T1_P1, (I14_ROD_CELLS, I14_ROD_DOFS)),
        ("Rb8", r4s2, ITEM305_ETA_T1_P1, (ITEM305_ROD_CELLS, ITEM305_ROD_DOFS)),
        ("Rb-r8", r8, R8_ETA_T1_P1, (GEO36_R8_CELLS, GEO36_R8_D2_DOFS)),
    ):
        s = f"{_eta(_cplx(js['E']), e_off)[iT1, 0]:.6f}"
        got = (js["cells"], js["ndofs"], off["cells"], off["ndofs"])
        ok_sizes = all(isinstance(x, int) for x in got) and got == (
            *sizes,
            I14_OFF_CELLS,
            I14_OFF_DOFS,
        )
        _say(
            f"{TAG} ({lab}) η_rod(T1) P1 {s} vs {want}: {s == want}; sizes "
            f"{got[0]}/{got[1]} off {got[2]}/{got[3]} vs {sizes[0]}/{sizes[1]} "
            f"{I14_OFF_CELLS}/{I14_OFF_DOFS}: {ok_sizes}"
        )
        if not (s == want and ok_sizes):
            fails.append((lab, s, got))
    rb_all = not fails

    built = _build_rod(implant_rod=dict(ROD, tip_ball_radius=TIP_BALL_4MM), **FIXTURE_KW)
    cad = built[3]["implant_rod_cad"]
    msh = built[0]
    cells = int(msh.topology.index_map(msh.topology.dim).size_global)
    r_ball = cad["tip_ball_radius_m"]
    _say(
        f"{TAG} build b4: tip_ball_radius_m {r_ball!r} (request {TIP_BALL_4MM!r}), "
        f"tip_size_m {cad['tip_size_m']!r} (TIP_R4 {TIP_R4!r}), {cells} cells, build "
        f"{built[4]:.2f} s beside item 464 {I17_B4_BUILD_S} s"
    )

    # (Kb4) the knob and the grain reached gmsh; the cell count is item 464's.
    kb_r = r_ball == TIP_BALL_4MM
    kb_g = cad["tip_size_m"] == TIP_R4
    kb_c = cells == I17_B4_CELLS
    _say(
        f"{TAG} (Kb4) tip_ball_radius_m == TIP_BALL_4MM bitwise: {kb_r}; tip_size_m == "
        f"TIP_R4 bitwise: {kb_g}; cells {cells} == I17_B4_CELLS {I17_B4_CELLS}: {kb_c} "
        f"(Δ {cells - I17_B4_CELLS:+d}; the dict(ROD) blind control builds 252 118)"
    )
    if not (kb_r and kb_g and kb_c):
        fails.append(("Kb4", r_ball, cad["tip_size_m"], cells))
    # (G1x) inside the patch; (N1) the unpatched (G1) would trip.
    patched = _ratio.ROD_CELL_GUARD_W2
    g1x = patched == ROD_CELL_GUARD_XL and cells < ROD_CELL_GUARD_XL
    n1 = cells > ROD_CELL_GUARD_W2
    _say(
        f"{TAG} DISCLOSURE: test_imp3_tip_ratio_1g.ROD_CELL_GUARD_W2 monkeypatched to "
        f"{patched} (ROD_CELL_GUARD_XL, an `xl` resource guard, not a band; original "
        f"{ROD_CELL_GUARD_W2}) on this path only; measured cells {cells}"
    )
    _say(
        f"{TAG} (G1x) patched == ROD_CELL_GUARD_XL {ROD_CELL_GUARD_XL}: "
        f"{patched == ROD_CELL_GUARD_XL}; cells {cells} < it: {cells < ROD_CELL_GUARD_XL} "
        f"({cells / ROD_CELL_GUARD_XL:.3f} of the guard; predicted 0.879)"
    )
    _say(
        f"{TAG} (N1) cells {cells} > unpatched ROD_CELL_GUARD_W2 {ROD_CELL_GUARD_W2}: {n1} "
        f"(ratio {cells / ROD_CELL_GUARD_W2:.3f}; predicted 1.098)"
    )
    if not g1x:
        fails.append(("G1x", patched, cells))
    if not n1:
        fails.append(("N1", cells, ROD_CELL_GUARD_W2))

    # (T1in) item 464's predicate (collective; every rank, outside any rank branch).
    tg = _t1_geometry(msh, comm)
    t1in = (
        tg["n_t1"] > 0
        and tg["t1_vmax"] <= r_ball
        and tg["owner"] >= 0
        and tg["owner_vmax"] <= r_ball
    )
    _say(
        f"{TAG} (T1in) near-T1 n {tg['n_t1']} > 0, vertex max {tg['t1_vmax']:.6e} <= "
        f"{r_ball:.6e}, T1 owner {tg['owner']} >= 0, owner vertex max "
        f"{tg['owner_vmax']:.6e} <= it: {t1in}; beside item 464 {I17_B4_NEAR_T1[0]} / "
        f"{I17_B4_T1_VMAX:.6e} / {I17_B4_T1_CELL[0]} / {I17_B4_T1_CELL[2]:.6e}; 2 mm control "
        f"vertex max {I17_2MM_T1_VMAX:.6e} > 2 mm (T1in-blind); TOP_CAP {TOP_CAP.tolist()}"
    )
    if not t1in:
        fails.append(("T1in", tg))

    # A wrong mesh never reaches the `xl` solve.
    mesh_fails = [f for f in fails if f[0] not in ("Rb", "Rb8", "Rb-r8")]
    if mesh_fails and not build_only:
        _say(f"{TAG} mesh anchors failed before the solve: {mesh_fails}")
        assert not fails, fails

    a = None
    if not build_only:
        a = _arm_d(built, "degree 2 rod b4", True, 2)
        fails += _identities(a, True)
        msh = a["case"]["mesh"]

    # (G3) tag 402 > 0 and all exterior.
    _, census = _outer_box_tags(msh, built[2], comm)
    _say(
        f"{TAG} (G3) 402 facets {census['wall']} (exterior {census['wall_exterior']}); "
        f"exterior {census['exterior']}"
    )
    if not (census["wall"] > 0 and census["wall_exterior"] == census["wall"]):
        fails.append(("G3", census))

    # Local mesh at T1 (collective; every rank, outside any rank branch).
    n_ball, ball_sum, n_t1, t1_sum = _ball_and_t1(msh, built[1], r_ball, comm)
    ball_edge = ball_sum / n_ball if n_ball else float("nan")
    t1_edge = t1_sum / n_t1 if n_t1 else float("nan")
    if a is None:
        tg_post = tg
        n_e, n_f = _ratio._edge_face_counts(msh)
        ndofs = 2 * n_e + 2 * n_f
    else:
        tg_post = _t1_geometry(msh, comm)
        n_e, n_f, ndofs = a["edges"], a["faces"], a["ndofs"]
    rss = _peak_rss_gib()

    _say(
        f"{TAG} mesh identity vs item 464 W2: cells {cells} / {I17_B4_CELLS} "
        f"({cells == I17_B4_CELLS}); 402 facets {census['wall']} / {I17_B4_N402} "
        f"({census['wall'] == I17_B4_N402}); d2 dofs {ndofs} (E {n_e}, F {n_f}"
        f"{', 2E + 2F' if a is None else ''}) / {I17_B4_D2_DOFS} (E {I17_B4_E}, F {I17_B4_F}) "
        f"({ndofs == I17_B4_D2_DOFS}); near-T1 ({T1_RADIUS:.1e} m) {n_t1} cells, mean edge "
        f"{t1_edge:.6e} m / {I17_B4_NEAR_T1} "
        f"({n_t1 == I17_B4_NEAR_T1[0] and f'{t1_edge:.6e}' == f'{I17_B4_NEAR_T1[1]:.6e}'})"
    )
    _say(
        f"{TAG} N_ball {n_ball}, mean ball edge {ball_edge:.6e} m / {I17_B4_N_BALL} / "
        f"{I17_B4_BALL_EDGE:.6e} ({n_ball == I17_B4_N_BALL and f'{ball_edge:.6e}' == f'{I17_B4_BALL_EDGE:.6e}'})"
    )
    t1c = (tg_post["owner"], tg_post["owner_edge"], tg_post["owner_vmax"])
    _say(
        f"{TAG} T1 cell (original index, mean edge, max vertex distance) "
        f"({t1c[0]}, {t1c[1]:.6e}, {t1c[2]:.6e}) / ({I17_B4_T1_CELL[0]}, "
        f"{I17_B4_T1_CELL[1]:.6e}, {I17_B4_T1_CELL[2]:.6e}) "
        f"({t1c[0] == I17_B4_T1_CELL[0] and f'{t1c[1]:.6e}' == f'{I17_B4_T1_CELL[1]:.6e}' and f'{t1c[2]:.6e}' == f'{I17_B4_T1_CELL[2]:.6e}'}); "
        f"near-T1 vertex max {tg_post['t1_vmax']:.6e} / {I17_B4_T1_VMAX:.6e} "
        f"({format(tg_post['t1_vmax'], '.6e') == format(I17_B4_T1_VMAX, '.6e')})"
    )

    if a is not None:
        payload = _ser(a)
        payload.update(
            {
                "rss_gib": rss,
                "ranks": comm.size,
                "tip_size_m": cad["tip_size_m"],
                "tip_ball_radius_m": r_ball,
                "near_t1": [n_t1, t1_edge],
                "t1_cell": [t1c[0], t1c[1]],
                "n_ball": [n_ball, ball_edge],
                "t1_vmax": tg_post["t1_vmax"],
                "guard_patched": patched,
            }
        )
        if comm.rank == 0:
            with open(target, "x") as fh:
                json.dump(payload, fh)
        comm.Barrier()
        if rb_all and not mesh_fails:
            _b4_readout(
                _eta(a["E"], e_off),
                {f"{t}|{d}": v for (t, d), v in a["sar"].items()},
                a["mass"],
                "b4",
                refs,
                t1_edge,
            )
        else:
            _say(f"{TAG} no η / Q reading printed: (Rb*) {rb_all}, mesh {mesh_fails}")
        _say(
            f"{TAG} price: cells {a['cells']}, E {a['edges']}, F {a['faces']}, dofs "
            f"{a['ndofs']}, build {built[4]:.2f} s, sweep {a['t_sweep']:.2f} s, SARs "
            f"{a['t_sar']:.2f} s, arm {a['t_arm']:.1f} s, window body "
            f"{time.perf_counter() - t0:.1f} s, summed peak RSS {rss:.2f} GiB at -n "
            f"{comm.size}; beside entry 33 {ENTRY33_PRICE}; entry 43 predicted "
            f"{ENTRY43_PRICE_PRED}"
        )
    else:
        if rb_all:
            _b4_readout(
                _eta(_cplx(r8["E"]), e_off),
                r8["sar"],
                r8["mass"],
                "stand-in",
                refs,
                t1_edge,
                prefix="DRY-RUN (r/8 JSON stand-in — not a reading) ",
            )
        else:
            _say(f"{TAG} no dry-run readout printed: (Rb*) red")
        _say(
            f"{TAG} b4build price: build {built[4]:.2f} s (item 464 {I17_B4_BUILD_S} s), body "
            f"{time.perf_counter() - t0:.1f} s, summed peak RSS {rss:.2f} GiB at -n "
            f"{comm.size}; no solve, no JSON"
        )
    _say(f"{TAG} anchors failed (all): {fails if fails else 'none'}")
    del a, built, msh
    gc.collect()
    assert not fails, fails


@complex_only
def test_imp3_tip_ladder_1g(monkeypatch):
    # (N2) the unpatched guard on every call, before dispatch: r4s2 never runs
    # inside the r8 patch (pytest undoes it at teardown).
    assert _ratio.ROD_CELL_GUARD_W2 == ROD_CELL_GUARD_W2 == 320_000, (
        _ratio.ROD_CELL_GUARD_W2,
        ROD_CELL_GUARD_W2,
    )
    if KNOB not in KNOBS:
        pytest.skip(f"FEM_EM_IMP3_STEP7 must be one of {KNOBS}")
    if KNOB == "r4s2":
        _r4s2()
        return
    monkeypatch.setattr(_ratio, "ROD_CELL_GUARD_W2", ROD_CELL_GUARD_XL)
    if KNOB in ("b4", "b4build"):
        _b4(build_only=KNOB == "b4build")
        return
    _r8(build_only=KNOB == "r8build")
