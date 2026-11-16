"""`POST-8` step 2 — the point-owner rule gated on `ans:2`'s named points.

PROJECT_PLAN §9 item 177; brief ``docs/planning/chunks/POST-8.md``.  Width is
``comm.size``; the module runs at ``-n 1``, ``-n 2`` and ``-n 4`` on one commit.

Construction: the `ans:2` fixture (``build_four_port_sweep`` on the 0.0025 m
1 g rung, four single-drive ``_solve_driven`` solves P1…P4 in the order
``_build_mass_averaged`` derives them) without its 21 ball integrals.
`ans:2`'s script is loaded by importlib only for ``_pointwise_sar`` and
``_phantom_power`` — **never** ``main()``, which writes tracked files.

Anchors (asserted):

* **A0** global cells == ``ONE_GRAM_RUNG_CELL_RECORD`` exactly; P1 phantom power
  within ``POWER_RTOL`` of ``TRACKED_P1_POWER_W`` (imported from the step-1
  probe, never restated).
* **A1** each named point has exactly one owned containing cell globally.
* **A2 (the gate)** the 16 raw N1curl SARs σ|E|²/(2ρ), evaluated through
  ``evaluate_vector_field_parallel(..., owner="lowest_original_cell")``, each
  lie within 1e-9 relative of ``STEP1_N1CURL_SAR`` — the mean of step 1's
  ``-n 2`` and ``-n 4`` per-cell prints (10 significant figures; comparand
  budget ≤ 5e-10).  Because every width is held to the *same* comparand, each
  window's assertion bounds the cross-width pairwise difference at ≤ 2e-9;
  the measured pairwise maximum across the three logs is a record, not an
  assertion.
* **A3** on the owning rank the rule's value equals
  ``e_complex.eval(x, [owner_local_cell])`` bitwise; all ranks hold the same
  16 values bitwise.

Negative control, the CG1 route (``ans2._pointwise_sar``, default
``owner=None``): **C1a** (``-n 2`` / ``-n 4``) the P1 CG1 row reproduces this
width's logged row to ``ROW_RTOL``; **C1b** (``-n 2`` / ``-n 4``) the same row
against the other width's logged row has min |rel| ≥ 0.10.  Printed: C1c
(``-n 1`` CG1 row beside both logged rows, predicted ≥ 1 % from each), C1d
(CG1 / raw ratio per point and drive), C2 (owner ``local_to_global`` vs
``original_cell_index``).
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest
from mpi4py import MPI

from dolfinx import geometry

from fem_em_solver.post.evaluation import (
    evaluate_vector_field_parallel,
    point_owner_cells,
)

from tests.complex_mode import complex_only
from tests.mesh.test_birdcage_port_sheets import SHEET_IFACE
from tests.validation.test_birdcage_b1_plus_map import _solve_driven
from tests.validation.test_birdcage_b1_quadrature import QUADRATURE_STEP_DEG, _port_index
from tests.validation.test_birdcage_sar_1g_rung import (
    ONE_GRAM_RUNG_CELL_RECORD,
    PHANTOM_RESOLUTION_1G_RUNG,
)
from tests.validation.test_port_birdcage_four_port import build_four_port_sweep

REPO = Path(__file__).resolve().parents[2]
ANS2_PATH = (
    REPO
    / "examples"
    / "ansys_benchmarks"
    / "ans2_birdcage_coil_driven_sar_10MHz"
    / "02_birdcage_coil_driven_sar_10MHz.py"
)
PROBE_PATH = REPO / "scripts" / "probes" / "post8_pointwise_partition.py"

#: Step 1's per-cell raw N1curl SAR prints (W/kg), [drive P1..P4][point +x, +y,
#: -x, -y].  ``-n 2``: ``20260927T112115Z_POST-8-step1-n2.log:1918–1966``;
#: ``-n 4``: ``20260927T112501Z_POST-8-step1-n4.log:1920–1968``.  The two agree
#: in 8 of 16 entries and differ by one unit in the 10th figure in the others.
_STEP1_N2 = np.array(
    [
        [5.514140460e-07, 8.056036637e-09, 6.147781056e-08, 1.000625322e-08],
        [5.413193615e-09, 5.365427637e-07, 7.538407567e-09, 6.062393526e-08],
        [6.887651145e-08, 7.937373888e-09, 5.038331027e-07, 7.475028492e-09],
        [9.595073748e-09, 6.651644441e-08, 9.294666807e-09, 4.980486020e-07],
    ]
)
_STEP1_N4 = np.array(
    [
        [5.514140460e-07, 8.056036636e-09, 6.147781056e-08, 1.000625321e-08],
        [5.413193614e-09, 5.365427637e-07, 7.538407566e-09, 6.062393526e-08],
        [6.887651145e-08, 7.937373887e-09, 5.038331027e-07, 7.475028491e-09],
        [9.595073747e-09, 6.651644441e-08, 9.294666808e-09, 4.980486020e-07],
    ]
)
#: The comparand: the mean of the two prints (budget ≤ 5e-10 relative).
STEP1_N1CURL_SAR = 0.5 * (_STEP1_N2 + _STEP1_N4)
#: The gate, as ruled (09-30 weekly).
A2_RTOL = 1e-9
#: C1b's separation floor (measured minimum 0.1293 in step 1).
C1B_MIN_SEPARATION = 0.10

pytestmark = complex_only


def _load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def case():
    comm = MPI.COMM_WORLD
    ans2 = _load(ANS2_PATH, "post8_step2_ans2_module")
    probe = _load(PROBE_PATH, "post8_step2_probe_module")
    assert not ans2.RESOLUTION_OVERRIDDEN

    sweep = build_four_port_sweep(phantom_resolution=PHANTOM_RESOLUTION_1G_RUNG)
    azimuths = {
        f"P{s['tag'] - SHEET_IFACE}": float(s["azimuth_deg"]) for s in sweep["sheets"]
    }
    order = sorted(azimuths)
    indices = {
        pid: _port_index(azimuths[pid], azimuths["P1"], QUADRATURE_STEP_DEG)
        for pid in order
    }
    solves = {pid: _solve_driven(sweep, pid) for pid in order}
    built = {
        "mesh": sweep["mesh"],
        "cell_tags": sweep["cell_tags"],
        "cells": int(sweep["cells"]),
        "solves": solves,
    }
    r = ans2.CENTRE_RADIUS_M
    points = np.array(
        [[+r, 0.0, 0.0], [0.0, +r, 0.0], [-r, 0.0, 0.0], [0.0, -r, 0.0]],
        dtype=np.float64,
    )
    if comm.rank == 0:
        print(
            f"\n[POST-8 step 2] -n {comm.size}: cells {built['cells']}; drives {order}; "
            f"port indices {indices}; points r = {r} m",
            flush=True,
        )
    return {"ans2": ans2, "probe": probe, "built": built, "points": points, "comm": comm}


def _sar(ans2, e):
    return float(ans2.SALINE_SIGMA) * float(np.sum(np.abs(e) ** 2)) / (
        2.0 * float(ans2.PHANTOM_RHO_KG_PER_M3)
    )


def test_a0_a1_mesh_power_and_single_cell(case):
    comm, built, ans2, probe = case["comm"], case["built"], case["ans2"], case["probe"]
    points = case["points"]
    powers = ans2._phantom_power(built)
    p1 = powers[0]["dissipated_power_w"]
    rel = abs(p1 / probe.TRACKED_P1_POWER_W - 1.0)

    msh = built["mesh"]
    n_owned = msh.topology.index_map(3).size_local
    tree = geometry.bb_tree(msh, 3)
    coll = geometry.compute_colliding_cells(
        msh, geometry.compute_collisions_points(tree, points), points
    )
    local = np.array(
        [int(np.sum(np.asarray(coll.links(i)) < n_owned)) for i in range(4)],
        dtype=np.int64,
    )
    total = np.empty_like(local)
    comm.Allreduce(local, total, op=MPI.SUM)
    if comm.rank == 0:
        print(
            f"\n[POST-8 A0] -n {comm.size}: cells {built['cells']} vs record "
            f"{ONE_GRAM_RUNG_CELL_RECORD}; P1 phantom power {p1:.17e} W vs "
            f"{probe.TRACKED_P1_POWER_W:.17e}  rel {rel:.3e} (<= {probe.POWER_RTOL:g})"
        )
        print(f"[POST-8 A1] -n {comm.size}: owned containing cells per point {total.tolist()}")
    assert built["cells"] == ONE_GRAM_RUNG_CELL_RECORD
    assert rel <= probe.POWER_RTOL
    assert total.tolist() == [1, 1, 1, 1]


def test_a2_a3_rule_gate_and_cg1_negative_control(case):
    comm, built, ans2, probe = case["comm"], case["built"], case["ans2"], case["probe"]
    points = case["points"]
    width = comm.size
    msh = built["mesh"]
    imap = msh.topology.index_map(3)

    owner_original, owner_local = point_owner_cells(msh, points, comm)
    l2g_local = np.full(4, -1, dtype=np.int64)
    for i in range(4):
        if owner_local[i] >= 0:
            l2g_local[i] = int(imap.local_to_global(np.array([owner_local[i]], dtype=np.int32))[0])
    l2g = np.empty_like(l2g_local)
    comm.Allreduce(l2g_local, l2g, op=MPI.MAX)

    sar = np.zeros((4, 4))
    e_all = []
    a3_local_ok = True
    for k in range(4):
        e = built["solves"][f"P{k + 1}"]["fields"].e_complex
        vals, valid = evaluate_vector_field_parallel(
            e, points, comm, owner="lowest_original_cell"
        )
        assert bool(np.all(valid))
        vals = np.asarray(vals)
        e_all.append(vals)
        for i in range(4):
            sar[k, i] = _sar(ans2, vals[i])
            if owner_local[i] >= 0:
                direct = np.asarray(
                    e.eval(points[i : i + 1], np.array([owner_local[i]], dtype=np.int32))
                ).reshape(-1)
                if direct.tobytes() != np.ascontiguousarray(vals[i]).tobytes():
                    a3_local_ok = False
    a3_ok = comm.allreduce(int(a3_local_ok), op=MPI.MIN) == 1
    e_bytes = np.ascontiguousarray(np.stack(e_all)).tobytes()
    a3_bitwise_ranks = all(b == e_bytes for b in comm.allgather(e_bytes))

    rel = np.abs(sar / STEP1_N1CURL_SAR - 1.0)
    budget = np.abs(_STEP1_N2 - _STEP1_N4) / STEP1_N1CURL_SAR

    # Negative control: the CG1 route at default owner=None.
    pointwise = ans2._pointwise_sar(built, points, comm)
    cg = np.array([pointwise[k]["sar_w_per_kg"] for k in range(4)])
    row = cg[0]

    if comm.rank == 0:
        print(f"\n[POST-8 A2] -n {width}: raw N1curl SAR under owner='lowest_original_cell'")
        for i in range(4):
            print(
                f"  point {i} {points[i].tolist()}: owner original_cell_index "
                f"{int(owner_original[i])}, local_to_global {int(l2g[i])}"
            )
        for k in range(4):
            for i in range(4):
                e = e_all[k][i]
                print(
                    f"  P{k + 1} pt{i}: SAR {sar[k, i]:.17e}  ref {STEP1_N1CURL_SAR[k, i]:.10e}  "
                    f"rel {rel[k, i]:.3e}  E = [{', '.join(f'{z.real:.17e}{z.imag:+.17e}j' for z in e)}]"
                )
        print(
            f"[POST-8 A2] -n {width}: max rel {rel.max():.3e} (<= {A2_RTOL:g}); "
            f"comparand budget max {budget.max():.3e}"
        )
        print(
            f"[POST-8 A3] -n {width}: owner-rank eval bitwise {a3_ok}; "
            f"all ranks bitwise {a3_bitwise_ranks}"
        )
        print(f"[POST-8 C1] -n {width}: P1 CG1 row " + "  ".join(f"{v:.12e}" for v in row))
        for w in (2, 4):
            ref = np.array(probe.LOGGED_P1_ROW[w])
            print(
                f"  vs logged -n {w} ({probe.LOGGED_P1_SOURCE[w]}): rel "
                + "  ".join(f"{a / b - 1.0:+.4e}" for a, b in zip(row, ref))
            )
        print(f"[POST-8 C1d] -n {width}: CG1 / raw ratio [drive][point]")
        for k in range(4):
            print(f"  P{k + 1}: " + "  ".join(f"{c / s:.4f}" for c, s in zip(cg[k], sar[k])))
        print(
            f"[POST-8 C2] -n {width}: original {owner_original.tolist()} "
            f"local_to_global {l2g.tolist()}"
        )

    # A2 — the gate
    assert float(rel.max()) <= A2_RTOL, rel
    # A3
    assert a3_ok
    assert a3_bitwise_ranks

    # C1a / C1b at the logged widths only
    if width in (2, 4):
        own = np.array(probe.LOGGED_P1_ROW[width])
        other = np.array(probe.LOGGED_P1_ROW[4 if width == 2 else 2])
        c1a = float(np.max(np.abs(row / own - 1.0)))
        c1b = float(np.min(np.abs(row / other - 1.0)))
        if comm.rank == 0:
            print(
                f"[POST-8 C1a] -n {width}: max rel vs own logged row {c1a:.3e} "
                f"(<= {probe.ROW_RTOL:g}); [C1b] min |rel| vs other width {c1b:.4f} "
                f"(>= {C1B_MIN_SEPARATION})"
            )
        assert c1a <= probe.ROW_RTOL
        assert c1b >= C1B_MIN_SEPARATION
