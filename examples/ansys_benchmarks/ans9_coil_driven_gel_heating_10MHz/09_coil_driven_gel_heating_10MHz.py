"""`ANS-9` step 1d — the runnable half of the coil-driven gel-heating SPEC (`ans:9`).

Our side of ``SPEC.md`` rows 1-7.  **Primary column (step 1d, §9 item 410):**
drive P1 of the exact-rectangle F-small birdcage at 10 MHz, κ off, through
route B (`GEO-27` step 2's builder, item 396's h0025 window):
``_build(True, phantom_resolution=0.0025, leg_sheet_strip_fraction=…)`` →
``build_four_port_sweep(built=…, build_only=True, degree=2)`` → one-factorisation
sweep with ``keep_fields=True``; N1curl degree 2 EM (201 072 cells), P1 thermal on
the 58 960-cell phantom submesh, G = 1 W / P₁, CN Δt = 31.35 s with every step
kept, the first Rannacher half-step and one BE calorimetric step.

**Degree-1 column (carried, not re-solved):** ``D1_RAGGED_COLUMN`` — N1curl
degree 1, ragged pre-10-03 sheets, `MAT-15` step 1 fixture, 2026-10-04 run
(``20261004T142653Z_ANS-9-step1a-W1.log``), adjudicated 2026-10-10.

Anchors (asserted, every band imported): (A0) SPEC == imports; (P1W), (G2), (i),
(R) on the EM half; (R1) row 1 = 1 W; (E) every run's balance residual; (R3)
row 3's calorimetric slope; (R6) the energy balance with row 6 built as an
independent per-(half-)step surface sum.  Negative controls (rule (e)): (d6)
doubled source, (s6) surface term dropped.  Printed: V_sub beside the 58 866-cell
record, the reproduction against item 396's W3, the degree-2 / degree-1 ratio.

Run through the example runner (heavy tier, ≈ 6 min)::

    ./scripts/run_examples.sh -e ans:9 -n 8 -t 590

All compute is inside ``main()``: `MAT-19` and the step 1b / 1c modules import
this file by path for ``SPEC_THERMAL``.  No AED value is read or written here
(CLAUDE.md hard rule).
"""

from __future__ import annotations

import gc
import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np
import ufl
from dolfinx import default_scalar_type, fem
from mpi4py import MPI

_REPO_ROOT = Path(__file__).resolve().parents[3]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from fem_em_solver.io.paraview_utils import adopt_host_ownership, write_xdmf_with_tags  # noqa: E402
from fem_em_solver.post.evaluation import evaluate_vector_field_parallel  # noqa: E402
from fem_em_solver.post.setup_figure import write_setup_figure  # noqa: E402

from tests.mesh.test_birdcage_port_sheets import PORT_LOWER, PORT_UPPER  # noqa: E402
from tests.mesh.test_birdcage_port_tags import LEG_COUNT, PHANTOM_HEIGHT  # noqa: E402
from tests.validation import test_mat10_sar_to_heat_sphere as m10  # noqa: E402
from tests.validation import test_mat12_transient_sar_heating_sphere as m12  # noqa: E402
from tests.validation import test_mat15_coil_driven_heating_phantom as m15  # noqa: E402
from tests.validation.test_mat15_coil_driven_heating_phantom import (  # noqa: E402
    A,
    B_BAND,
    D_CONTROL_FLOOR,
    DT_1,
    E_BAND,
    H_GEL,
    K_GEL,
    N_MAIN,
    PHANTOM_CELL_TAG,
    RHO_C,
    ROBIN_TAG,
    S2_BAND,
    T_INF,
    T_STAR,
    V_SUB_RECORD_M3,
    _be_one,
    _transient,
)

CASE_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = CASE_DIR / "paraview_output"
BASENAME = "ans9_coil_driven_gel_heating_10mhz"  # no dots (write_xdmf_with_tags suffix rule)
FIGURE_DIR = CASE_DIR / "figures"  # committed, EX-57
FIGURE_BASENAME = "ans9_coil_driven_gel_heating_10MHz"
ANS2_DIR = _REPO_ROOT / "examples" / "ansys_benchmarks" / "ans2_birdcage_coil_driven_sar_10MHz"

#: The SPEC's typed numbers (SPEC.md "Thermal problem" / "Quantities to export").
#: Used ONLY in the (A0) structure check; every computation uses the imports.
SPEC_THERMAL = {
    "k": 0.6,
    "rho_c": 4.18e6,
    "h": 10.0,
    "t_inf": 20.0,
    "t_init": 20.0,
    "duration": 1254.0,
    "dt_bound": 31.35,
    "a": 0.030,
    "half_height": 0.040,
    "times": (313.5, 627.0, 940.5, 1254.0),
    "probes": {
        "A": (0.0, 0.0, 0.0),
        "B": (0.020, 0.0, 0.0),
        "C": (0.0, 0.020, 0.0),
        "D": (-0.020, 0.0, 0.0),
        "E": (0.0, 0.0, 0.030),
    },
}
SNAP_STEPS = (10, 20, 30, 40)

#: The degree-1 column, carried as a labelled record (not re-solved): every value
#: of the tracked ``metrics.json`` at HEAD 59809ba7, written by the 2026-10-04 run
#: (docs/testing/logs/20261004T142653Z_ANS-9-step1a-W1.log).  Its (R0), (R6), (d6)
#: and (s6) were asserted then and are not re-asserted on this column now.
D1_RAGGED_LABEL = ("N1curl degree 1, ragged pre-10-03 sheets, MAT-15 step 1 fixture, "
                   "2026-10-04 run, adjudicated 2026-10-10")
D1_RAGGED_COLUMN = {
    "label": D1_RAGGED_LABEL,
    "provenance_log": "docs/testing/logs/20261004T142653Z_ANS-9-step1a-W1.log",
    "dolfinx_version": "0.11.0.post0",
    "em_element": "N1curl degree 1",
    "thermal_element": "Lagrange P1",
    "em_cells": 199920,
    "thermal_cells": 58866,
    "mpi_ranks": 2,
    "row1_power_w": 0.9999999999999971,
    "row2_mean_rise_k": {
        "313.5": 0.31809585128727264,
        "627.0": 0.614314551495653,
        "940.5": 0.892104352771447,
        "1254.0": 1.153313650151018,
    },
    "row3_slope_k_per_s": 0.0010585632847679603,
    "row3_ideal_k_per_s": 0.0010585693760630137,
    "row4_probe_rise_k": {
        "A": 0.9057980632357207,
        "B": 3.5648535011976605,
        "C": 0.6853046677515235,
        "D": 0.6974435308905691,
        "E": 0.9586388205588361,
    },
    "row5_max_rise_k": 4.457206255150254,
    "row5_max_location_m": [0.027450583355725953, 0.0015640398390111719, -0.0002930574018659884],
    "row6_surface_heat_j": 164.49781315196608,
    "row7_incident_power_w": 45428.07551777395,
    "row7_supplied_power_w": 124582.6616703806,
    "row7_accepted_power_w": 38253.74370012777,
    "first_half_step_slope_k_per_s": 0.001053210730618553,
    "anchors": {
        "R1": 2.886579864025407e-15,
        "V_sub": 3.447242491461111e-13,
        "R0": 2.524227680444612e-15,
        "R3": 5.75427098892245e-06,
        "R6": 7.09409493757295e-15,
        "d6": 1.0000000000000138,
        "s6": 0.13117847938752605,
    },
    "wall_time_s": 172.23091343790293,
    "thermal_reruns_s": 4.180006755050272,
}

#: Item 396's W3 (h0025-d2x, same fixture and route), printed reproduction
#: comparands only: docs/testing/logs/20261012T100436Z_ANS-9-step1c-W3-h0025.log
W3_G = 1.833792433645e07  # :1986
W3_SAR = {"A": 3.512880216e-01, "B": 1.749800995e01, "C": 5.395029877e-03,  # :1987-1991
          "D": 2.979742845e00, "E": 2.656157242e00}
W3_ROW2 = {10: 0.318123, 20: 0.614370, 30: 0.892184, 40: 1.153415}  # :1994
W3_ROW4 = {"A": 0.889137, "B": 3.514465, "C": 0.683153, "D": 0.756873, "E": 0.959953}  # :1997
W3_R3 = 5.6446e-06  # :1996

COLUMN_FIXTURE = "route B, exact rectangles, κ off, h_p 0.0025 m"


def _load_example(path: Path, module_name: str):
    spec = importlib.util.spec_from_file_location(module_name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


def _rel(a, b):
    return abs(a - b) / abs(b) if b != 0 else abs(a)


def _row_cells(m):
    """Row-table value strings for one metrics-shaped column (rows 1-7)."""
    t = m["row2_mean_rise_k"]
    p = m["row4_probe_rise_k"]
    xm = m["row5_max_location_m"]
    out = [("1", "power deposited, W", f"{m['row1_power_w']:.12f}")]
    out += [("2", f"<theta>({k} s), K", f"{t[k]:.6f}") for k in t]
    out.append(("3", "initial slope, K/s", f"{m['row3_slope_k_per_s']:.6e}"))
    out += [("4", f"theta probe {k} at 1254 s, K", f"{p[k]:.6f}") for k in p]
    out.append(("5", "max theta at 1254 s, K (location m)",
                f"{m['row5_max_rise_k']:.6f} ({xm[0]:+.4f}, {xm[1]:+.4f}, {xm[2]:+.4f})"))
    out.append(("6", "surface heat 0-1254 s, J", f"{m['row6_surface_heat_j']:.4f}"))
    out.append(("7", "incident power at P1 for 1 W, W", f"{m['row7_incident_power_w']:.4f}"))
    return out


_VERDICT_MARK = b"## Verdict"


def _verdict_tail(data: bytes) -> bytes:
    """Everything from the first line starting ``## Verdict``, byte for byte (or b"")."""
    if data.startswith(_VERDICT_MARK):
        return data
    idx = data.find(b"\n" + _VERDICT_MARK)
    return data[idx + 1:] if idx >= 0 else b""


def _write_comparison(case_dir: Path, metrics: dict) -> Path:
    """Re-render ``COMPARISON.md``: the generated head from ``metrics`` (degree-2
    primary column) and ``D1_RAGGED_COLUMN``, then the existing file's hand-added
    ``## Verdict`` section re-appended byte for byte (item 400's pattern)."""
    path = Path(case_dir) / "COMPARISON.md"
    tail = _verdict_tail(path.read_bytes()) if path.exists() else b""
    prim = _row_cells(metrics)
    d1 = _row_cells(D1_RAGGED_COLUMN)
    lines = [
        "# ANS-9 — comparison (our side; AED columns blank)",
        "",
        "Generated by `09_coil_driven_gel_heating_10MHz.py` (`ans:9`). Drive P1, P1 thermal, "
        "CN θ = 0.5, Δt = 31.35 s, 40 steps, Rannacher start. **Primary column (step 1d, "
        f"§9 item 410):** {metrics.get('em_element', 'N1curl degree 2')}, "
        f"{metrics.get('fixture', COLUMN_FIXTURE)}, {metrics.get('em_cells', '—')} EM / "
        f"{metrics.get('thermal_cells', '—')} thermal cells. **Degree-1 column (carried, "
        f"not re-solved):** {D1_RAGGED_LABEL}.",
        "",
        "| Row | Quantity | Ours (degree 2, exact, h_p 0.0025) | Ours (degree 1, ragged; 2026-10-04) "
        "| AED | ours vs AED |",
        "|---|---|---|---|---|---|",
        *[f"| {r} | {q} | {v} | {d1v} | — | — |" for (r, q, v), (_, _, d1v) in zip(prim, d1)],
        "",
        "## Readout (pre-registered, copied verbatim from SPEC.md)",
        "",
        "- Rows 2, 3 and 6 are integral quantities: agreement within **5 %** (the",
        "  level `ANS-2`'s phantom power integral was adjudicated at) is AGREE.",
        "- Row 4 is pointwise on a source with a hot region near the driven leg:",
        "  within **10 %** at probes A–D, with B / C / D ordered the same way in both",
        "  codes, is AGREE; probe E is a record.",
        "- Row 5 is a record (a maximum on a mapped field is mesh-dependent).",
        "- Row 7 against ours measures the feed, not the thermal solve; it is",
        "  `ANS-2`'s question and is recorded, not adjudicated here.",
        "- A disagreement on rows 2 / 3 with row 1 at 1 W on both sides is a finding",
        "  about one code's thermal solve or source mapping; a known-issues entry and",
        "  a diagnosis chunk follow.",
        "",
    ]
    head = "\n".join(lines)
    if tail:
        head += "\n"
    path.write_bytes(head.encode("utf-8") + tail)
    return path


def _draw_setup_figure(msh, cell_tags, comm):
    """`EX-57` setup figure: geometry only, opt-in via ``FEM_EM_SETUP_FIGURES=1``.

    Collective, a no-op otherwise.  Not regenerated by step 1d (the committed
    figure is the 2026-10-04 fixture's geometry; disclosed in the guide).
    """
    from tests.validation.test_port_birdcage_lumped_column import CONDUCTOR_CELL_TAG

    names = {CONDUCTOR_CELL_TAG: "conductor", 2: "air", PHANTOM_CELL_TAG: "phantom"}
    for i in range(1, LEG_COUNT + 1):
        names[PORT_LOWER + i] = f"leg L{i} port lower"
        names[PORT_UPPER + i] = f"leg L{i} port upper"
    write_setup_figure(
        msh,
        cell_tags,
        FIGURE_DIR / f"{FIGURE_BASENAME}_setup.png",
        region_names=names,
        hide_tags=(2,),
        translucent_tags=(PHANTOM_CELL_TAG,),
        slice_normal=(0.0, 0.0, 1.0),
        title="ans:9 -- coil-driven gel heating, F-small birdcage + phantom, 10 MHz",
        comm=comm,
    )


def _thermal_rows(sub, ft, q_h, v_sub, comm, say, check):  # noqa: C901
    """Rows 1-6, (R1), (E), (R3), (R6), (d6), (s6) for one DG0 source on ``sub``.

    The step 1a row code, unchanged in its schedule: the snapshot CN run, the
    first Rannacher half-step as its own run, one BE dt_1 step, the doubled
    source.  Every reduction is allreduced on every rank.
    """
    dt_main = T_STAR / N_MAIN
    ds = ufl.Measure("ds", domain=sub, subdomain_data=ft)(ROBIN_TAG)

    def allsum(form):
        return float(np.real(comm.allreduce(fem.assemble_scalar(fem.form(form)), op=MPI.SUM)))

    # ---- row 1 -------------------------------------------------------------
    int_q = allsum(q_h * ufl.dx(domain=sub))
    r1 = abs(int_q - 1.0)
    say(f"\nROW1 int_sub q_h = {int_q:.15e} W  |.-1| = {r1:.3e}  (predicted <= 2e-14)")
    check("(R1) row 1 = 1 W", r1 <= B_BAND, f"{r1:.3e} <= {B_BAND:.0e}")

    def snap_run(q):
        return _transient(sub, ft, q, dt=dt_main, n_steps=N_MAIN, theta=0.5,
                          startup_half_steps=2, snapshot_steps=range(N_MAIN + 1))

    def half_run(q):
        return _transient(sub, ft, q, dt=dt_main / 2.0, n_steps=1, theta=1.0)

    def surface_sum(res, res_half):
        """Row 6 independently: the scheme's theta-weighted conv_form per (half-)step."""
        th = {n: res.snapshots[n] - T_INF for n in range(N_MAIN + 1)}
        th_half = res_half.temperature - T_INF
        terms = [(dt_main / 2.0) * allsum(H_GEL * th_half * ds),
                 (dt_main / 2.0) * allsum(H_GEL * th[1] * ds)]
        for n in range(2, N_MAIN + 1):
            terms.append(dt_main * allsum(H_GEL * 0.5 * (th[n] + th[n - 1]) * ds))
        return float(np.sum(terms)), terms

    t1 = time.perf_counter()
    res = snap_run(q_h)
    res_half = half_run(q_h)
    res_be = _be_one(sub, ft, q_h)
    t_thermal = time.perf_counter() - t1

    e_runs = {}
    for name, rr in (("snapshot CN", res), ("half-step BE", res_half), ("BE dt_1", res_be)):
        e_runs[name] = rr.max_balance_residual
        check(f"(E) {name} max_balance", rr.max_balance_residual <= E_BAND,
              f"{rr.max_balance_residual:.3e} <= {E_BAND:.0e} (predicted <= 1e-11)")

    # ---- row 2 -------------------------------------------------------------
    mean_th = {n: allsum((res.snapshots[n] - T_INF) * ufl.dx(domain=sub)) / v_sub for n in SNAP_STEPS}
    bound = {n: n * dt_main / (RHO_C * v_sub) for n in SNAP_STEPS}
    say("\nROW2 <theta>(t) K: " + "  ".join(
        f"t={n * dt_main:.1f}s {mean_th[n]:.6f} (no-loss {bound[n]:.4f})" for n in SNAP_STEPS)
        + "  (predicted 0.318123 / 0.614370 / 0.892184 / 1.153415)")

    # ---- row 3 -------------------------------------------------------------
    thbar1 = allsum((res_be.temperature - T_INF) * ufl.dx(domain=sub)) / v_sub
    slope = thbar1 / DT_1
    ideal = 1.0 / (RHO_C * v_sub)
    r3 = abs(RHO_C * thbar1 / (DT_1 * int_q / v_sub) - 1.0)
    say(f"\nROW3 slope = {slope:.9e} K/s  vs 1 W/(rho c V_sub) = {ideal:.9e} K/s  "
        f"(predicted 1.058563e-3 / 1.058569e-3)")
    check("(R3) calorimetric slope", r3 <= S2_BAND, f"{r3:.4e} <= {S2_BAND:.0e} (predicted 5.6446e-06)")
    thbar_half = allsum((res_half.temperature - T_INF) * ufl.dx(domain=sub)) / v_sub
    slope_half = thbar_half / (dt_main / 2.0)
    deficit_half = 1.0 - RHO_C * thbar_half / ((dt_main / 2.0) * int_q / v_sub)
    say(f"PRINT main-run first half-step slope = {slope_half:.9e} K/s, deficit {deficit_half:.4e} "
        f"(not row 3)")

    # ---- row 4 -------------------------------------------------------------
    names = list(SPEC_THERMAL["probes"])
    pts = np.array([SPEC_THERMAL["probes"][k] for k in names], dtype=np.float64)
    th40_f = fem.Function(res.snapshots[N_MAIN].function_space, name="theta_40")
    th40_f.x.array[:] = res.snapshots[N_MAIN].x.array - T_INF
    th40_f.x.scatter_forward()
    # owner="lowest_original_cell" raises on a create_submesh mesh ("Mixed topology
    # unsupported.", 20261004T142026Z_ANS-9-step1a-W1.log:1977); theta is continuous
    # P1, the default collective path is used with identical points on every rank.
    vals, valid = evaluate_vector_field_parallel(th40_f, pts, comm)
    probe = {k: float(np.real(np.ravel(vals[i])[0])) for i, k in enumerate(names)}
    check("(A0) all five probes valid", bool(np.all(valid)), f"valid = {list(map(bool, valid))}")
    say("\nROW4 theta(1254 s) K: " + "  ".join(f"{k} {probe[k]:.6f}" for k in names)
        + f"  order B>C>D: {probe['B'] > probe['C'] > probe['D']}  B>D>C: "
        f"{probe['B'] > probe['D'] > probe['C']}  (predicted A 0.889137 B 3.514465 C 0.683153 "
        f"D 0.756873 E 0.959953)")

    # ---- row 5 -------------------------------------------------------------
    V = th40_f.function_space
    n_own = V.dofmap.index_map.size_local
    xyz = V.tabulate_dof_coordinates()[:n_own]
    loc_vals = np.real(th40_f.x.array[:n_own])
    if n_own:
        i = int(np.argmax(loc_vals))
        cand = (float(loc_vals[i]), tuple(float(c) for c in xyz[i]))
    else:
        cand = (-np.inf, (np.nan, np.nan, np.nan))
    allc = comm.allgather(cand)
    tmax, xmax = max(allc, key=lambda c: c[0])
    rmax = float(np.hypot(xmax[0], xmax[1]))
    azmax = float(np.degrees(np.arctan2(xmax[1], xmax[0])))
    say(f"\nROW5 max theta(1254 s) = {tmax:.6f} K at ({xmax[0]:+.5f}, {xmax[1]:+.5f}, {xmax[2]:+.5f}) m; "
        f"r = {rmax:.4f} m, azimuth {azmax:+.1f} deg (predicted 4.3-4.5 K, r in [0.020, 0.030], |az| <= 30)")

    # ---- row 6 + (R6) ------------------------------------------------------
    th40 = res.snapshots[N_MAIN] - T_INF
    row6, terms = surface_sum(res, res_half)
    stored = RHO_C * allsum(th40 * ufl.dx(domain=sub))
    deposited = T_STAR * int_q
    r6 = abs(deposited - stored - row6) / deposited
    say(f"\nROW6 surface heat 0-1254 s = {row6:.6f} J (predicted 164.40 J); stored {stored:.6f} J; "
        f"deposited {deposited:.6f} J; surface fraction {row6 / deposited:.4f} (predicted 0.1311)")
    check("(R6) energy balance", r6 <= E_BAND, f"{r6:.3e} <= {E_BAND:.0e} (predicted <= 1e-12)")
    if r6 > E_BAND:
        say("ROW6 per-step terms J: " + " ".join(f"{t:.6e}" for t in terms))

    # ---- controls -----------------------------------------------------------
    q_two = fem.Function(q_h.function_space, name="q_two")
    q_two.x.array[:] = 2.0 * q_h.x.array
    q_two.x.scatter_forward()
    res_d = snap_run(q_two)
    res_dh = half_run(q_two)
    row6_d, _ = surface_sum(res_d, res_dh)
    stored_d = RHO_C * allsum((res_d.snapshots[N_MAIN] - T_INF) * ufl.dx(domain=sub))
    d6 = abs(deposited - stored_d - row6_d) / deposited
    check("(d6) doubled source misses (R6)", d6 > D_CONTROL_FLOOR,
          f"{d6:.6f} > {D_CONTROL_FLOOR:.0e} (predicted 1.000)")
    s6 = abs(deposited - stored - 0.0) / deposited
    check("(s6) surface term dropped misses (R6)", s6 > E_BAND,
          f"{s6:.4f} > {E_BAND:.0e} (predicted 0.1311)")

    return {"res": res, "int_q": int_q, "r1": r1, "e_runs": e_runs, "mean_th": mean_th,
            "slope": slope, "ideal": ideal, "r3": r3, "slope_half": slope_half,
            "deficit_half": deficit_half, "probe": probe, "tmax": tmax, "xmax": xmax,
            "row6": row6, "r6": r6, "d6": d6, "s6": s6, "t_thermal": t_thermal}


def main() -> None:  # noqa: C901
    comm = MPI.COMM_WORLD
    rank0 = comm.rank == 0
    t_wall = time.perf_counter()

    if not np.issubdtype(np.dtype(default_scalar_type), np.complexfloating):
        raise RuntimeError(
            "this example needs the complex DolfinX build "
            "(source /usr/local/bin/dolfinx-complex-mode); the runner's `ans:` "
            "group sources it automatically"
        )

    # Route B and its bands, imported by name inside main() (the top level stays
    # import-only for MAT-19 / 1b / 1c, which load this file by path).
    from tests.mesh.test_birdcage_port_sheet_prerequisite import CELL_COUNT_BAND
    from tests.mesh.test_birdcage_port_sheets import _build
    from tests.validation.test_ans9_far_side_sar import _column_readout
    from tests.validation.test_ans9_far_side_sar_refined import _recip, _sweep
    from tests.validation.test_birdcage_sar_mass_averaged import (
        EXACT_IDENTITY_RTOL,
        _global_tag_cell_count,
    )
    from tests.validation.test_geo27_one_gram_rung_degree2 import RUNGS
    from tests.validation.test_imp3_palace_twin import _sheet_areas
    from tests.validation.test_imp3_rod_birdcage import _peak_rss_gib
    from tests.validation.test_imp3_tip_ladder_1g import I14_OFF_CELLS, I14_OFF_DOFS
    from tests.validation.test_imp3_tip_ratio_1g import NOMINAL_STRIP_AREA
    from tests.validation.test_mat10_sar_to_heat_sphere import _robin_tags, _to_sub
    from tests.validation.test_mat15_coil_driven_heating_phantom import (
        _cell_average_dg0_tag,
        _tag_submesh,
    )
    from tests.validation.test_ops61_exact_strip_ports import STRIP_BAND
    from tests.validation.test_port_birdcage_four_port import build_four_port_sweep
    from tests.validation.test_port_lumped_narrowed_sheet import GATED_WIDTH_FRACTION
    from tests.validation.test_port_lumped_sheet_sweep import RECIPROCITY_BAND

    def say(msg):
        if rank0:
            print(msg, flush=True)

    failures = []

    def check(name, ok, detail):
        say(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}")
        if not ok:
            failures.append(name)

    say("=" * 72)
    say("ANS-9 step 1d — coil-driven gel heating, 10 MHz, drive P1 (our side), degree 2")
    say("=" * 72)

    # ---- (A0) structure ----------------------------------------------------
    dt_main = T_STAR / N_MAIN
    pairs = {
        "k": (SPEC_THERMAL["k"], K_GEL),
        "rho_c": (SPEC_THERMAL["rho_c"], RHO_C),
        "h": (SPEC_THERMAL["h"], H_GEL),
        "t_inf": (SPEC_THERMAL["t_inf"], T_INF),
        "t_init": (SPEC_THERMAL["t_init"], T_INF),
        "duration": (SPEC_THERMAL["duration"], T_STAR),
        "a": (SPEC_THERMAL["a"], A),
        "half_height": (SPEC_THERMAL["half_height"], PHANTOM_HEIGHT / 2.0),
    }
    for i, frac in enumerate((0.25, 0.5, 0.75, 1.0)):
        pairs[f"time_{i}"] = (SPEC_THERMAL["times"][i], T_STAR * frac)
        pairs[f"time_{i}_step"] = (SNAP_STEPS[i] * dt_main, T_STAR * frac)
    for name, (spec_v, imp_v) in pairs.items():
        r = _rel(spec_v, imp_v)
        check(f"(A0) SPEC {name}", r <= 1e-12, f"SPEC {spec_v!r} vs import {imp_v!r} rel {r:.1e}")
    check("(A0) dt <= SPEC bound", dt_main <= SPEC_THERMAL["dt_bound"] * (1 + 1e-12),
          f"dt = {dt_main!r} vs {SPEC_THERMAL['dt_bound']}")
    structural = {
        "m15.E_BAND is m12.E_BAND": m15.E_BAND is m12.E_BAND,
        "m15.S2_BAND is m12.S2_BAND": m15.S2_BAND is m12.S2_BAND,
        "m15.B_BAND is m10.B_BAND": m15.B_BAND is m10.B_BAND,
    }
    for name, ok in structural.items():
        check(f"(A0) {name}", ok, str(ok))
    if failures:
        raise SystemExit(f"(A0) red — SPEC and MAT-15 disagree: {failures}")

    # ---- the EM half: route B, h0025, degree 2 -----------------------------
    h_p = RUNGS["h0025"]
    comm.Barrier()
    t0 = time.perf_counter()
    built = _build(True, phantom_resolution=h_p, leg_sheet_strip_fraction=GATED_WIDTH_FRACTION)
    case = build_four_port_sweep(built=built, build_only=True, degree=2)
    comm.Barrier()
    t_build = time.perf_counter() - t0
    msh, ct = case["mesh"], case["cell_tags"]
    cells = int(case["cells"])
    tag3 = _global_tag_cell_count(msh, ct, PHANTOM_CELL_TAG, comm)
    areas = _sheet_areas(msh, case["facet_tags"], comm)
    say(f"\n[EM] route B h_p = {h_p} m: {cells} cells ({tag3} tag-3; predicted 201072 / 58960), "
        f"build {t_build:.2f} s")
    t_fig = time.perf_counter()
    _draw_setup_figure(msh, ct, comm)
    t_fig = time.perf_counter() - t_fig  # excluded from the printed wall time
    check("(G2) four strip areas", len(areas) == 4, f"{len(areas)} sheets")
    for k, a in areas.items():
        rr = abs(a / NOMINAL_STRIP_AREA - 1.0)
        check(f"(G2) area {k}", rr <= STRIP_BAND, f"{a:.12e} vs {NOMINAL_STRIP_AREA:.6e} rel {rr:.3e} "
              f"<= {STRIP_BAND:.0e}")

    comm.Barrier()
    t0 = time.perf_counter()
    r2 = _sweep(case, 2)
    comm.Barrier()
    t_sweep = time.perf_counter() - t0
    s2 = np.asarray(r2.s_matrix, dtype=np.complex128)
    recip = _recip(s2)
    v2 = r2.fields["P1"].e_complex.function_space
    dofs2 = int(v2.dofmap.index_map.size_global * v2.dofmap.index_map_bs)
    say(f"[EM] degree-2 sweep {t_sweep:.2f} s, port_ids {tuple(r2.port_ids)}, d2 dofs {dofs2}")
    check("(A0) r2.port_ids[0] == 'P1'", r2.port_ids[0] == "P1", f"port_ids[0] = {r2.port_ids[0]}")
    check("(i) reciprocity", recip <= RECIPROCITY_BAND,
          f"||S-S^T||/||S|| = {recip:.3e} <= {RECIPROCITY_BAND:.0e} (predicted 1.3e-14)")
    cell_rel = abs(cells / I14_OFF_CELLS - 1.0)
    check("(R) cells vs I14_OFF_CELLS", cell_rel <= CELL_COUNT_BAND,
          f"{cells} vs {I14_OFF_CELLS} rel {cell_rel:.3e} <= {CELL_COUNT_BAND}")
    check("(R) d2 dofs == I14_OFF_DOFS", dofs2 == I14_OFF_DOFS, f"{dofs2} vs {I14_OFF_DOFS}")

    f2 = r2.fields["P1"]
    rd = _column_readout("h0025-d2x", f2.e_complex, f2.sigma_field, ct, comm)
    G = rd["G"]
    check("(P1W) all five points valid", all(rd["valid"]), f"valid = {rd['valid']}")
    check("(P1W) independent phantom integral = 1 W", rd["p1w"] <= EXACT_IDENTITY_RTOL,
          f"{rd['p1w']:.3e} <= {EXACT_IDENTITY_RTOL:.0e} (predicted <= 1e-14)")
    say(f"[repro W3] G rel {_rel(G, W3_G):.3e}; SAR " + "  ".join(
        f"{k} {_rel(rd['sar'][k], W3_SAR[k]):.3e}" for k in W3_SAR) + " (predicted <= 1e-6, printed)")

    sub, parent_ids = _tag_submesh(msh, ct, PHANTOM_CELL_TAG)
    ft = _robin_tags(sub)
    Wp, q_par = _cell_average_dg0_tag(msh, ct, PHANTOM_CELL_TAG,
                                      G * 0.5 * f2.sigma_field * ufl.inner(f2.e_complex, f2.e_complex))
    q_h = _to_sub(Wp, q_par, sub, parent_ids)

    # Row 7's current: the package's sheet_terminal_current from the same sweep.
    spec_p1 = next(s for s in case["specs"] if s.port_id == "P1")
    solved = {"driven": "P1", "source_voltage_v": complex(spec_p1.drive_voltage_v),
              "currents": {"P1": complex(r2.excitation_results["P1"].responses["P1"].current_a)}}
    del f2, r2, Wp, q_par
    gc.collect()

    one = fem.Constant(sub, default_scalar_type(1.0))
    v_sub = float(np.real(comm.allreduce(fem.assemble_scalar(fem.form(one * ufl.dx(domain=sub))),
                                         op=MPI.SUM)))
    thermal_cells = int(sub.topology.index_map(3).size_global)
    say(f"PRINT V_sub = {v_sub:.13e} m^3 vs 58 866-cell record {V_SUB_RECORD_M3:.13e} "
        f"rel {_rel(v_sub, V_SUB_RECORD_M3):.3e} (printed; thermal cells {thermal_cells})")

    # ---- the thermal half: rows 1-6 ---------------------------------------
    tr = _thermal_rows(sub, ft, q_h, v_sub, comm, say, check)
    mean_th, probe = tr["mean_th"], tr["probe"]
    say("[repro W3] row 2 " + "  ".join(f"n={n} {_rel(mean_th[n], W3_ROW2[n]):.3e}" for n in SNAP_STEPS)
        + "; row 4 " + "  ".join(f"{k} {_rel(probe[k], W3_ROW4[k]):.3e}" for k in W3_ROW4)
        + f"; (R3) {tr['r3']:.4e} vs W3 {W3_R3:.4e} (printed; W3 values carry 6-7 digits)")

    # ---- row 7 -------------------------------------------------------------
    ans2 = _load_example(ANS2_DIR / "02_birdcage_coil_driven_sar_10MHz.py",
                         "ans2_birdcage_coil_driven_sar_10mhz_for_ans9")
    exc = ans2._excitation(solved)
    row7 = G * exc["incident_power_w"]
    say(f"\nROW7 incident at P1 for 1 W = {row7:.6f} W (predicted 45844.81); supplied "
        f"{G * exc['supplied_power_re_w']:.6f} W, accepted {G * exc['accepted_power_re_w']:.6f} W")

    # ---- degree 2 / degree 1 ----------------------------------------------
    d1 = D1_RAGGED_COLUMN
    say("\nRATIO degree 2 / degree 1 (ragged, 2026-10-04): row 2 " + "  ".join(
        f"{k} {mean_th[n] / d1['row2_mean_rise_k'][k]:.6f}"
        for n, k in zip(SNAP_STEPS, d1["row2_mean_rise_k"]))
        + f"; row 3 {tr['slope'] / d1['row3_slope_k_per_s']:.6f}; row 4 "
        + "  ".join(f"{k} {probe[k] / d1['row4_probe_rise_k'][k]:.6f}" for k in probe)
        + f"; row 5 {tr['tmax'] / d1['row5_max_rise_k']:.6f}; row 6 "
        f"{tr['row6'] / d1['row6_surface_heat_j']:.6f}; row 7 {row7 / d1['row7_incident_power_w']:.6f}")

    wall = time.perf_counter() - t_wall - t_fig
    rss = _peak_rss_gib()
    say(f"PRICE build {t_build:.2f} s sweep {t_sweep:.2f} s thermal {tr['t_thermal']:.2f} s "
        f"body {wall:.2f} s peakRSS {rss:.2f} GiB ranks {comm.size}")

    # ---- artifacts ----------------------------------------------------------
    metrics = {
        "case": "ANS-9 step 1d — coil-driven gel heating, 10 MHz, drive P1 (our side)",
        "solver": "fem-em-solver (FEniCSx)",
        "dolfinx_version": __import__("dolfinx").__version__,
        "em_element": "N1curl degree 2",
        "thermal_element": "Lagrange P1",
        "fixture": COLUMN_FIXTURE,
        "em_cells": cells,
        "em_dofs": dofs2,
        "thermal_cells": thermal_cells,
        "scheme": "CN theta = 0.5, dt = 31.35 s, 40 steps, Rannacher start of 2 BE half-steps",
        "mpi_ranks": comm.size,
        "row1_power_w": tr["int_q"],
        "row2_mean_rise_k": {f"{n * dt_main:.1f}": mean_th[n] for n in SNAP_STEPS},
        "row3_slope_k_per_s": tr["slope"],
        "row3_ideal_k_per_s": tr["ideal"],
        "row4_probe_rise_k": probe,
        "row5_max_rise_k": tr["tmax"],
        "row5_max_location_m": list(tr["xmax"]),
        "row6_surface_heat_j": tr["row6"],
        "row7_incident_power_w": row7,
        "row7_supplied_power_w": G * exc["supplied_power_re_w"],
        "row7_accepted_power_w": G * exc["accepted_power_re_w"],
        "first_half_step_slope_k_per_s": tr["slope_half"],
        "point_sar_w_per_kg": rd["sar"],
        "anchors": {"P1W": rd["p1w"], "reciprocity": recip, "R1": tr["r1"],
                    "V_sub_vs_58866_record": _rel(v_sub, V_SUB_RECORD_M3), "R3": tr["r3"],
                    "R6": tr["r6"], "E": tr["e_runs"], "d6": tr["d6"], "s6": tr["s6"]},
        "wall_time_s": wall,
        "thermal_reruns_s": tr["t_thermal"],
        "columns": {"d1_ragged": D1_RAGGED_COLUMN},
    }
    if rank0:
        (CASE_DIR / "metrics.json").write_text(json.dumps(metrics, indent=2) + "\n", encoding="utf-8")
        _write_comparison(CASE_DIR, metrics)

    res = tr["res"]
    fields = {}
    for n in SNAP_STEPS:
        nm = f"theta_t{int(n * dt_main + 1e-6):04d}"  # 0313, 0627, 0940, 1254 (truncated s)
        f = fem.Function(res.snapshots[n].function_space, name=nm)
        f.x.array[:] = res.snapshots[n].x.array - T_INF
        f.x.scatter_forward()
        fields[nm] = f
    q_h.name = "q_h"
    fields["q_h"] = q_h
    OUTPUT_DIR.mkdir(exist_ok=True)
    xdmf_path, _ = write_xdmf_with_tags(OUTPUT_DIR / f"{BASENAME}_combined", sub, None, fields,
                                        comm=comm)
    adopt_host_ownership(OUTPUT_DIR, comm=comm)
    # write_xdmf_with_tags returns the path on rank 0 only (None elsewhere).
    xdmf_name = xdmf_path.name if xdmf_path is not None else ""
    say(f"\n[output] {xdmf_name}: P1 " + ", ".join(k for k in fields if k != "q_h")
        + " (K) and DG0 q_h (W/m^3); metrics.json, COMPARISON.md")
    say(f"ANS-9 wall time {wall:.1f} s (thermal re-runs {tr['t_thermal']:.1f} s)")

    if failures:
        say(f"ANS-9 anchors: FAIL {failures}")
        raise SystemExit(1)
    say("ANS-9 anchors and controls: PASS")


if __name__ == "__main__":
    main()
