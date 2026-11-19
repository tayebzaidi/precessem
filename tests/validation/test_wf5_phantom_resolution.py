"""`WF-5` step 3 (§9 item 143, chain T16): the phantom's 1/Q under ``phantom_resolution``.

Step 1's fixture (the tuned F-small birdcage on `ANS-7`'s Route A, κ held at the
tracked ``metrics.json`` value) rebuilt with the phantom refined through
``th14._sheets_build(True, as_hole=True, phantom_resolution=h_p)`` — one rung per
window, chosen by ``FEM_EM_WF5_PHANTOM_RESOLUTION`` (unset ⇒ skipped by name).

Asserted:

* (R1) the knob reached the mesh: global cells ≠ the default 80 181
  (``20261001T141511Z_WF-5-step1-W2.log:959``);
* (R2) per column |P_loss − (p_src − p_sheet(P1))|/P_loss ≤ ``DISCRETE_IDENTITY_RTOL``
  — an exact discrete identity on any mesh;
* (R3) exactly one − → + crossing of Im Y_in on 63.7 … 64.3 MHz (``_rows_on_grid``);
* control (binary, W2 only): global tag-3 cells at 0.0025 exceed those at 0.005
  (read from W1's JSON in the gitignored ``logs/wf5/``), backed by `GEO-27`'s
  monotone ladder.

Printed (the record for T13): cells, phantom cells, f₀, Q_slope, Q_field and
|Q_field − Q_slope|/Q_slope beside B_Q 3e-3 (*predicted* only — never asserted on
a new mesh), P_ph/(ω₀W) and its relative change from step 1's default-mesh value,
in-run κ beside the tracked κ.  Two rungs, degree 1, a record — no h-convergence
claim.  Tests only, complex mode.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pytest
from mpi4py import MPI

from tests.complex_mode import complex_only
from tests.validation.test_wf5_loading import STEP1_PI_PH, UNLOADED_MHZ, _rows_on_grid, _w_tot
from tests.validation.test_wf5_q_two_ways import (
    ANS7_DIR,
    ANS7_SCRIPT,
    COLUMNS,
    _budget,
    _load,
    _y_in,
)

_ENV = "FEM_EM_WF5_PHANTOM_RESOLUTION"
_REPO_ROOT = Path(__file__).resolve().parents[2]
JSON_DIR = _REPO_ROOT / "logs" / "wf5"
DEFAULT_CELLS = 80181  # 20261001T141511Z_WF-5-step1-W2.log:959
CONTROL_COARSE = 0.005
CONTROL_FINE = 0.0025
# *Predicted, printed only* (§9 item 143): step 1's B_Q is not carried to a new mesh.
PREDICTED_B_Q = 3.0e-3
PREDICTED_CELLS = {0.005: 94_000, 0.0025: 164_000}
PREDICTED_SHARE_CHANGE = 2.0e-2
SHARE_STOP = 5.0e-2

_H_RAW = os.environ.get(_ENV)

pytestmark = pytest.mark.skipif(_H_RAW is None, reason=f"{_ENV} unset (one rung per window)")


def _h():
    return float(_H_RAW)


@pytest.fixture(scope="module")
def ans7():
    return _load(ANS7_SCRIPT, "ans7_driven_sweep_wf5")


def _global_tag_cells(msh, cell_tags, tag, comm):
    n_local = msh.topology.index_map(msh.topology.dim).size_local
    owned = cell_tags.indices < n_local
    return int(comm.allreduce(int(np.count_nonzero(cell_tags.values[owned] == tag)), op=MPI.SUM))


@pytest.fixture(scope="module")
def fixture_build(ans7):
    """`test_wf5_q_two_ways.py:99–117` with ``phantom_resolution=h_p``; κ tracked."""
    from tests.validation import test_th14_birdcage_copper as th14

    comm = MPI.COMM_WORLD
    h_p = _h()
    tracked = json.loads((ANS7_DIR / "metrics.json").read_text())
    kappa = float(tracked["kappa"])
    msh, cell_tags, facet_tags, _diag, t_mesh = th14._sheets_build(
        True, as_hole=True, phantom_resolution=h_p
    )
    tdim = msh.topology.dim
    msh.topology.create_connectivity(tdim - 1, tdim)
    msh.topology.create_entity_permutations()
    outer_tags, _census = th14._outer_box_tags(msh, facet_tags, comm)
    tags_k, port_defs_k, specs_k = th14._build_ports(
        msh, cell_tags, comm, width_correction_kappa=kappa
    )
    geo_k = (msh, cell_tags, facet_tags, outer_tags, tags_k, port_defs_k)
    n_cells = int(msh.topology.index_map(tdim).size_global)
    n_ph = _global_tag_cells(msh, cell_tags, th14.PHANTOM_CELL_TAG, comm)
    # In-run κ (printed only): uncorrected specs, `_stage_a:696,724`.
    tags_u, port_defs_u, specs_unc = th14._build_ports(msh, cell_tags, comm)
    geo_u = (msh, cell_tags, facet_tags, outer_tags, tags_u, port_defs_u)
    kappa_run = ans7._kappa_in_run(th14, geo_u, specs_unc, comm)
    if comm.rank == 0:
        print(
            f"[wf5-3] phantom_resolution {h_p}: hole mesh {n_cells} cells (default "
            f"{DEFAULT_CELLS}; PREDICTED ~ {PREDICTED_CELLS.get(h_p)}), phantom (tag "
            f"{th14.PHANTOM_CELL_TAG}) cells {n_ph}, mesh {t_mesh:.1f} s",
            flush=True,
        )
        print(
            f"[wf5-3] kappa tracked {kappa:.9e} (held); kappa in-run {kappa_run:.9e} "
            f"(rel {abs(kappa_run - kappa) / abs(kappa):.4e}, printed only)",
            flush=True,
        )
    return {
        "th14": th14,
        "geo": geo_k,
        "specs": specs_k,
        "tracked": tracked,
        "comm": comm,
        "h_p": h_p,
        "n_cells": n_cells,
        "n_ph": n_ph,
        "kappa": kappa,
        "kappa_run": kappa_run,
    }


@pytest.fixture(scope="module")
def sweep(ans7, fixture_build):
    fb = fixture_build
    out = {}
    for col in COLUMNS:
        y = []
        for f_mhz in UNLOADED_MHZ:
            s11, _ = ans7._route_a_solve(fb["th14"], fb["geo"], fb["specs"], f_mhz * 1e6, col)
            y.append(_y_in(ans7, s11, f_mhz * 1e6))
        out[col] = np.array(y)
    return out


@complex_only
def test_phantom_share_record(ans7, fixture_build, sweep):
    """R1–R3 (+ the W2-only cell-count control); the T13 record printed."""
    from tests.validation.test_birdcage_power_identity import DISCRETE_IDENTITY_RTOL

    fb = fixture_build
    comm = fb["comm"]
    rank0 = comm.rank == 0
    h_p = fb["h_p"]
    fails = []
    # (R1)
    if rank0:
        print(
            f"[wf5-3] (R1) cells {fb['n_cells']} != default {DEFAULT_CELLS} (ASSERTED)", flush=True
        )
    if fb["n_cells"] == DEFAULT_CELLS:
        fails.append(("R1", fb["n_cells"]))
    record = {
        "phantom_resolution": h_p,
        "n_cells": fb["n_cells"],
        "n_phantom_cells": fb["n_ph"],
        "kappa_tracked": fb["kappa"],
        "kappa_in_run": fb["kappa_run"],
        "columns": {},
    }
    for col in COLUMNS:
        rows = _rows_on_grid(UNLOADED_MHZ, sweep[col])  # (R3) asserted inside
        f0 = rows["f0_hz"]
        b = _budget(ans7, fb, f0, col)
        omega = b["omega"]
        w_tot = _w_tot(b)
        p_loss = b["p_phantom"] + b["p_surf"]
        inv_q_field = p_loss / (omega * w_tot)
        q_field, q_slope = 1.0 / inv_q_field, 1.0 / rows["inv_q_sec"]
        miss = abs(q_field - q_slope) / q_slope
        id_rel = abs(p_loss - (b["p_src"] - b["p_sheet1"])) / p_loss
        share = b["p_phantom"] / (omega * w_tot)
        d_share = (share - STEP1_PI_PH[col]) / STEP1_PI_PH[col]
        if rank0:
            print(
                f"[wf5-3] {col}: (R3) one crossing in {rows['bracket_mhz']} MHz; f0 "
                f"{f0:.6f} Hz (cubic {rows['f0_cubic_hz']:.6f}); Q_slope {q_slope:.6f}, "
                f"Q_field {q_field:.6f}; |Q_field - Q_slope|/Q_slope = {miss:.4e} "
                f"(B_Q {PREDICTED_B_Q:.0e} PREDICTED only, not asserted on this mesh)",
                flush=True,
            )
            print(
                f"[wf5-3] {col}: secant G0 {rows['g0']:.9e} S, dB/df {rows['dbdf']:.9e} S/Hz; "
                f"cubic Q {1.0 / rows['inv_q_cubic']:.6f}, e_sec "
                f"{abs(rows['inv_q_cubic'] - rows['inv_q_sec']) / rows['inv_q_sec']:.4e}; "
                f"|Q_field - Q_cubic|/Q_cubic = "
                f"{abs(q_field * rows['inv_q_cubic'] - 1.0):.4e} (printed)",
                flush=True,
            )
            print(
                f"[wf5-3] {col}: (R2) P_phantom {b['p_phantom']:.9e} W, P_surf "
                f"{b['p_surf']:.9e} W, P_loss {p_loss:.9e} W; p_src {b['p_src']:.9e} W, "
                f"p_sheet(P1) {b['p_sheet1']:.9e} W; identity rel {id_rel:.3e} (ASSERTED <= "
                f"{DISCRETE_IDENTITY_RTOL:.0e})",
                flush=True,
            )
            print(
                f"[wf5-3] {col}: phantom share P_ph/(w0 W) {share:.6e} (default mesh "
                f"{STEP1_PI_PH[col]:.6e}); relative change {d_share:+.4e} (PREDICTED |.| <= "
                f"{PREDICTED_SHARE_CHANGE:.0e}; > {SHARE_STOP:.0e} is T13's finding, recorded)",
                flush=True,
            )
        if not id_rel <= DISCRETE_IDENTITY_RTOL:
            fails.append((col, "R2", id_rel))
        record["columns"][col] = {
            "f0_hz": f0,
            "q_slope": q_slope,
            "q_field": q_field,
            "q_miss": miss,
            "identity_rel": id_rel,
            "phantom_share": share,
            "phantom_share_change": d_share,
        }
    if rank0:
        JSON_DIR.mkdir(parents=True, exist_ok=True)
        (JSON_DIR / f"phantom_resolution_{h_p}.json").write_text(json.dumps(record, indent=2))
    comm.Barrier()
    # Control (binary, W2 only): tag-3 cells grow under refinement.
    if h_p == CONTROL_FINE:
        coarse = json.loads((JSON_DIR / f"phantom_resolution_{CONTROL_COARSE}.json").read_text())
        n_c = int(coarse["n_phantom_cells"])
        if rank0:
            print(
                f"[wf5-3] control: phantom cells {fb['n_ph']} at {CONTROL_FINE} > {n_c} at "
                f"{CONTROL_COARSE} (ASSERTED)",
                flush=True,
            )
        if not fb["n_ph"] > n_c:
            fails.append(("control", fb["n_ph"], n_c))
    assert not fails, f"WF-5 step 3 gates: {fails}"
