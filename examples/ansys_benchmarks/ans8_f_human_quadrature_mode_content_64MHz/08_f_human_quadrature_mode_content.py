"""`ANS-8` step 1 runnable half (§9 item 141, §10 chain H): the F-human
birdcage's mode content under the two-port quadrature drive, our side.

SPEC (``SPEC.md`` in this directory, v1 2026-10-01): F-human (16 legs, 32
ring gaps, the `WF-7` phantom) at 64 MHz, degree 1, σ = 800 S/m coil solved
inside, Route B — the 32×32 S of the 32 gap sheets, then circuit algebra:
C = 113.11 pF on every gap, the feeds P17 / P21 as Norton sources (V_f / 50
in parallel with 1/50 + jωC), drives (i) P17 alone, (ii) ccw ``(1, −j)``,
(iii) cw ``(1, +j)``.  Rows 1–3 come from the single shared
``fem_em_solver.ports.ring_modes.ring_mode_rows``, which the private writer
applies to AED's exported 32×32 as well.

Stages (``FEM_EM_ANS8_STAGE``), results in the gitignored
``/workspace/logs/ans8/`` (``FEM_EM_ANS8_RESULTS_DIR`` overrides):

* ``S`` (``-n 8``) — builds the 32×32 exactly as ``test_wf10_step1c_c_sweep``
  W1 does (``tests/validation/test_wf10_f_human_tuned_circuit.py:452–495``,
  imported pieces only: ``_build_f_human``, ``_ring_quadrature_slots``,
  ``_solve_one_drive``, ``terminal_form_deficit``, ``_respec``): κ from the
  uncorrected P17 solve, then the κ-corrected 32 columns.  Asserts (G1) the
  cell count within ``CELL_COUNT_BAND``, reciprocity ≤ ``RECIPROCITY_BAND``,
  σ_max ≤ ``MATRIX_PASSIVITY_CEILING``; rank 0 saves ``s32.npz``.
* ``R`` (``-n 1``) — loads that npz (missing ⇒ stop; never recomputes),
  re-checks G1, asserts (G2) the measured slot map equals the SPEC naming
  map, (G3) the matrix against item 128's held npz if present, the anchors
  A1–A6 (item 133's circuit records, `20261001T140319Z_WF-10-step1d-W1.log`)
  at ``C_MODE_F`` and the **cw** sense, and the negative control at C₂; then
  writes the SPEC rows at C = 1.1311e-10 F (``metrics.json``,
  ``COMPARISON.md``, AED column blank) and, unless
  ``FEM_EM_ANS_SKIP_PRIVATE=1``, the private table.
* SPEC v2 (§9 item 161), stage R after A1–A6: rows 2b (every gap at the
  matrix's own row-1 C_mode), 2c (the opposite pair and the four-port ccw at
  both C) and 4a (f₀ = f·√(C_mode / C)) through the shared
  ``ring_mode_rows_v2`` — anchors B0–B5, C1–C4, F1/F2 and two negative
  controls; ``metrics.json`` gains ``rows_v2`` with v1's ``rows`` unchanged.

Not claimed: any AED value or verdict, the sweep interpolation of row 4,
Route A, the loss partition, copper, the degree-2 32×32.
"""

from __future__ import annotations

import json
import os
import resource
import sys
import time
from pathlib import Path

import numpy as np

_REPO_ROOT = Path(__file__).resolve().parents[3]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

CASE_DIR = Path(__file__).resolve().parent
FIGURE_DIR = CASE_DIR / "figures"
FIGURE_BASENAME = "ans8_f_human_quadrature_mode_content_64MHz"
STAGE_ENV = "FEM_EM_ANS8_STAGE"
RESULTS_DIR_ENV = "FEM_EM_ANS8_RESULTS_DIR"
DEFAULT_RESULTS_DIR = "/workspace/logs/ans8"
NPZ_NAME = "s32.npz"
#: `ANS-4`'s shared knob: ``=1`` skips the private write even if the JSON exists.
SKIP_PRIVATE_ENV = "FEM_EM_ANS_SKIP_PRIVATE"
AED_RESULTS_PATH = CASE_DIR / "aed_results" / "ans8_aed_results.json"

#: SPEC: C on all 32 gaps (exactly 113.11 pF) and the feeds.
SPEC_C_F = 1.1311e-10
SPEC_FEED_NAMES = ("P17", "P21")
SPEC_DRIVES = {"p17_alone": (1.0, 0.0), "ccw": (1.0, -1.0j), "cw": (1.0, 1.0j)}
QUADRATURE_SLOT_OFFSET = 4

# Item 133's circuit records (`20261001T140319Z_WF-10-step1d-W1.log`), the
# anchors' literals.  C_MODE_F itself is imported (`test_wf10_…:310`).
C2_RECORD_F = 7.545440170601e-11                      # :10665
OFF_BLOCK_RECORD = 4.182887e-03                       # :10652
OFF_BLOCK_RTOL = 1.0e-6
SHARE_ATOL = 1.0e-6                                   # 6-digit printed shares
FOLD_ATOL = 1.0e-4                                    # 4-digit printed folds
A4_CW_SHARES = (0.790993, 0.952414)                   # :10666 (feed ring, other)
A4_CW_FOLDS = {"bottom": (0.1734, 0.7910, 0.0228, 0.0054, 0.0029, 0.0016, 0.0013, 0.0011, 0.0006),
               "top": (0.0004, 0.9524, 0.0299, 0.0073, 0.0039, 0.0021, 0.0017, 0.0015, 0.0008)}
A5_P17_SHARES = (0.789857, 0.951115)                  # :10660
A5_P17_FOLDS = {"bottom": (0.1735, 0.7899, 0.0227, 0.0062, 0.0029, 0.0018, 0.0013, 0.0011, 0.0005),
                "top": (0.0004, 0.9511, 0.0298, 0.0083, 0.0039, 0.0024, 0.0017, 0.0015, 0.0007)}
A6_OPPOSITE_SHARES = (0.988772, 0.987500)             # :10669
A6_FOUR_PORT_SHARES = (0.989921, 0.988704)            # :10672
NEGATIVE_CONTROL_RECORD = 0.099704                    # :10665 (parallel drive at C_2)
NEGATIVE_CONTROL_CEILING = 0.5
KAPPA_RECORD = 4.590991e-03                           # :10644, printed beside
PREDICTED_CCW_CW_GAP = 2.0e-3                         # SPEC, printed only

# SPEC v2 (§9 item 161): row 2c's symmetric-drive controls, by SPEC port name
# (mapped through ``port_map`` as SPEC_FEED_NAMES is).  The opposite pair is
# P17 / P25 alone (P21 a plain C gap), as A6's call; the four-port is **ccw**.
SPEC_CONTROLS = {"opposite": (("P17", "P25"), (1.0, -1.0)),
                 "four_ccw": (("P17", "P21", "P25", "P29"), (1.0, -1.0j, -1.0, 1.0j))}
# v2 anchors on our matrix (R = `20261001T184403Z_ANS-8-step1-W2-R.log`).
B3_ATOL = 1.0e-9                                      # own C vs C_MODE_F blocks (≈ 4e-13 expected)
B4_CCW_CW_GAP = 1.912e-3                              # R:55
B4_ATOL = 1.0e-6
B5_CCW_RANGE = (0.789075, 0.789087)                   # predicted, printed
B5_ROW2_MINUS_2B_FEED = (-2.0e-5, -1.0e-5)            # predicted, printed
B5_ROW2_MINUS_2B_OTHER = (0.0, 5.0e-6)                # predicted, printed
CONTROL_LEAK_BAR = 5.0e-5                             # item 141's bar; R:52 measured 9.46e-7
C2_PREDICTED_SHARE_ATOL = 5.0e-3                      # printed, vs A6-cw
C3_PREDICTED_SHARE_ATOL = 1.0e-4                      # printed
F1_RTOL = 1.0e-9
GAUGE_TWO_PORT_LEAK_FLOOR = 0.1                       # 0.1734 at C_mode, R:45

# SPEC v3 (§9 item 369, §10 H21): stage V3 on H9′'s tuned copper F-human.
V3_JSON_NAME = "v3-rows.json"
V3_SPEC_LABEL = "SPEC.md v3 (2026-10-10)"
#: Comparands only (logged readings, cited; never published as a row).
#: 358's mean |B₁⁺|_ccw over the masked centre sample
#: (`20261010T170356Z_WF-10-step2-W1.log:6339`), 7-digit print.
MEAN_B1P_CCW_H9P_RECORD_T = 2.884091e-07
#: 359's SAR10g at its c* (r = 45 mm, z_c, az₁₇)
#: (`20261010T171705Z_WF-10-step2s-W1.log:6316`).
SAR10G_CSTAR_H9P_RECORD_W_PER_KG = 4.852895959e-03
V3_CSTAR_LABEL_RECORD = "r=45 dz=+0 az17+0"
#: §9 item 386 amendment 2 — (Rep-P) / (Rep-Part) comparands: item 369's V3
#: window, same code / mesh / fixture (`20261011T141650Z_ANS-8-v3-V3.log:6194`,
#: P_acc,SPEC 13-digit print; `:6195`, the partition terms, 10-digit print).
P_ACC_SPEC_H21_RECORD_W = 3.057629563251e-02
PARTITION_H21_RECORD_W = {"p_rc30": 2.024950049e-02, "p_esr_feeds": 1.265338043e-03,
                          "p_phantom": 8.197883283e-03, "p_surf": 6.092761976e-04}
#: Amendment 3's pre-registered predictions (the 10:30 review's pathologist
#: arithmetic on `20261011T142050Z_ANS-8-v3-diag.log:37–39`), printed only.
V3_GAP_X_PREDICTED_W = {"P17": 1.32932e-4, "P21": 1.21365e-4}
V3_CS_OVER_X_WINDOW = (0.99, 1.01)
#: SPEC v3 probes (x, y) [m] at z = z_c: O, P₁ … P₄.
V3_PROBES_XY = (("O", 0.0, 0.0), ("P1", 0.045, 0.0), ("P2", 0.0, 0.045),
                ("P3", -0.045, 0.0), ("P4", 0.0, -0.045))
# Predictions (brief, drafter's unlogged arithmetic), printed only.
V3_PREDICTED = {
    "p_acc_spec_w": 3.03e-2, "frac_phantom": 0.270, "frac_copper": 0.020, "frac_esr": 0.710,
    "share_k1_feed": 0.9938, "leak_k0_max": 5e-3, "b1p_o_per_sqrt_w": 1.66e-6,
    "sar10g_per_w": 0.160, "purity_o": 21.5, "esr_blind": 0.71, "feed_esr_blind": 4.2e-2,
    "p361_in_pacc1": 4.2e-2,
}


def _results_dir() -> Path:
    return Path(os.environ.get(RESULTS_DIR_ENV, DEFAULT_RESULTS_DIR))


def _ring_modes():
    from fem_em_solver.ports import ring_modes
    return ring_modes


def _spec_port_names():
    return [f"P{n}" for n in range(17, 49)]


def _rss_summed_gib(comm) -> float:
    from mpi4py import MPI
    local = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss) * 1024.0 / 2**30
    return float(comm.allreduce(local, op=MPI.SUM))


def _draw_setup_figure(ctx, comm):
    """`EX-57` setup figure: geometry only (air hidden, phantom translucent).

    Opt-in via ``FEM_EM_SETUP_FIGURES=1``; collective, a no-op otherwise.  The
    32 ring-gap port boxes are the cell tags ``PORT_LOWER + i`` / ``PORT_UPPER
    + i``; the slice is the coil mid-plane (z = mesh centre), no half-space clip.
    """
    from fem_em_solver.post.setup_figure import write_setup_figure
    from tests.mesh.test_birdcage_port_sheets import PORT_LOWER, PORT_UPPER
    from tests.validation.test_port_birdcage_lumped_column import (
        CONDUCTOR_CELL_TAG,
        PHANTOM_CELL_TAG,
    )

    region_names = {2: "air", CONDUCTOR_CELL_TAG: "copper conductor",
                    PHANTOM_CELL_TAG: "phantom"}
    for spec in ctx["specs"]:
        i = int(spec.port_id[1:])
        region_names[PORT_LOWER + i] = f"port {spec.port_id}"
        region_names[PORT_UPPER + i] = f"port {spec.port_id}"
    for tag in range(PORT_LOWER + 1, PORT_LOWER + 17):  # the 16 leg-gap boxes
        region_names.setdefault(tag, "leg gap box (unported)")
    write_setup_figure(
        ctx["msh"],
        ctx["cell_tags"],
        FIGURE_DIR / f"{FIGURE_BASENAME}_setup.png",
        region_names=region_names,
        hide_tags=(2,),
        translucent_tags=(PHANTOM_CELL_TAG,),
        slice_normal=(0.0, 0.0, 1.0),
        title="ans:8 -- F-human birdcage, 32 ring-gap ports",
        window_size=(1400, 630),
        comm=comm,
    )


# -- stage S ----------------------------------------------------------------------
def _stage_s(comm, outdir):
    from tests.mesh.test_birdcage_port_sheet_prerequisite import CELL_COUNT_BAND
    from tests.validation.test_birdcage_f_human_rung import F_HUMAN_BRANCH_B_CELL_RECORD
    from tests.validation.test_port_birdcage_ring_column import (
        _reciprocity_ratio,
        _solve_one_drive,
    )
    from tests.validation.test_port_birdcage_ring_matrix import MATRIX_PASSIVITY_CEILING
    from tests.validation.test_port_drive_superposition import _ring_quadrature_slots
    from tests.validation.test_port_lumped_sheet_sweep import RECIPROCITY_BAND
    from tests.validation.test_wf10_f_human_tuned_circuit import _respec
    from tests.validation.test_wf7_f_human_b1_quadrature import (
        DEGREE,
        FREQUENCY_HZ,
        _build_f_human,
    )
    from fem_em_solver.ports.shares import terminal_form_deficit

    t_window = time.perf_counter()

    def say(msg):
        if comm.rank == 0:
            print(f"[ANS-8 S] {msg}", flush=True)

    say(f"start: -n {comm.size}, f = {FREQUENCY_HZ:.3e} Hz, degree {DEGREE}")
    built = _build_f_human(comm, emit_leg_ports=True)  # GEO-35 step 2b: unpin with the re-record
    n_cells = int(built["cells"])
    rel_cells = n_cells / F_HUMAN_BRANCH_B_CELL_RECORD - 1.0
    say(f"(G1) cells {n_cells} vs record {F_HUMAN_BRANCH_B_CELL_RECORD}: relative "
        f"{rel_cells:.3e} (band {CELL_COUNT_BAND}, ASSERTED); t={time.perf_counter() - t_window:.1f} s")
    if not abs(rel_cells) < CELL_COUNT_BAND:
        raise SystemExit(f"(G1) cells {n_cells}: relative {rel_cells:.3e} outside band — stop")
    ctx = built["ctx"]
    _draw_setup_figure(ctx, comm)
    base_specs = list(ctx["specs"])
    _az_ref, slots = _ring_quadrature_slots(built["sheets"])
    ref_pid = "P17"
    ref_k, ref_sign = slots[ref_pid]
    feed2 = next(pid for pid, (k, sg) in slots.items()
                 if sg == ref_sign and k == (ref_k + QUADRATURE_SLOT_OFFSET) % 16)
    others = sorted((pid for pid in slots if pid not in (ref_pid, feed2)),
                    key=lambda p: int(p[1:]))
    order = [ref_pid, feed2] + others
    say(f"feeds: {ref_pid} (slot {ref_k}, ring {ref_sign:+.0f}) and {feed2}")

    col = _solve_one_drive(ctx, ref_pid, reuse_factorization=False)
    fields = col.pop("fields")
    sheets = [spec.sheet(driven=(spec.port_id == ref_pid)) for spec in base_specs]
    kappa = float(terminal_form_deficit(ctx["msh"], ctx["tags_f"], sheets,
                                        fields.e_complex, comm)["pooled"])
    del fields
    say(f"kappa = {kappa:.6e} (item 133 printed {KAPPA_RECORD:.6e}); "
        f"t={time.perf_counter() - t_window:.1f} s")
    ctx32 = dict(ctx)
    ctx32["specs"] = _respec(base_specs, kappa)
    columns = {}
    for n_drive, pid in enumerate(order):
        c = _solve_one_drive(ctx32, pid, reuse_factorization=(n_drive > 0))
        c.pop("fields")
        columns[pid] = c
    kinds = [columns[p]["solve_kind"] for p in order]
    s32 = np.array([[columns[pj]["s_column"][pi] for pj in order] for pi in order],
                   dtype=np.complex128)
    ratio32 = _reciprocity_ratio(s32)
    sigma32 = float(np.linalg.norm(s32, 2))
    say(f"32x32: {kinds.count('solve')} factorising, {kinds.count('held')} held; (G1) "
        f"reciprocity {ratio32:.6e} (band {RECIPROCITY_BAND:.0e}, ASSERTED); sigma_max "
        f"{sigma32:.9f} (ceiling {MATRIX_PASSIVITY_CEILING}, ASSERTED); "
        f"t={time.perf_counter() - t_window:.1f} s")
    if not (ratio32 <= RECIPROCITY_BAND and sigma32 <= MATRIX_PASSIVITY_CEILING):
        raise SystemExit("(G1) reciprocity / passivity red — stop")
    if comm.rank == 0:
        outdir.mkdir(parents=True, exist_ok=True)
        np.savez(outdir / NPZ_NAME, s32=s32, kappa=kappa, order=np.array(order),
                 slot_k=np.array([slots[p][0] for p in order]),
                 slot_sign=np.array([slots[p][1] for p in order]),
                 cells=n_cells, freq_hz=FREQUENCY_HZ, degree=DEGREE,
                 reciprocity=ratio32, sigma_max=sigma32)
        say(f"saved {outdir / NPZ_NAME}")
    comm.Barrier()
    rss = _rss_summed_gib(comm)
    say(f"window elapsed {time.perf_counter() - t_window:.1f} s at -n {comm.size}; "
        f"summed ru_maxrss {rss:.3f} GiB")


# -- SPEC v3 (stage V3) -----------------------------------------------------------
def p_acc_line(w, cur_feeds, z_th, z_rc):
    """SPEC v3's accepted power at the 50 Ω line cut, ``½ Re Σ V_port conj(I_line)``.

    Each feed is a 50 Ω source in parallel with its own gap's ``Z_RC``, reduced
    to the Thevenin pair ``(w_i, Z_th)``: ``V_port = w − Z_th·I`` and the line
    current ``I_line = I + V_port / Z_RC`` (premise 1 of §9 item 369).
    Returns ``(P_acc, V_port, I_line)``.
    """
    w = np.asarray(w, dtype=np.complex128)
    i = np.asarray(cur_feeds, dtype=np.complex128)
    v = w - z_th * i
    i_line = i + v / z_rc
    return float(0.5 * np.real(np.sum(v * np.conj(i_line)))), v, i_line


def p_acc_wave(s2, w, v_s, z0):
    """The wave form ``(‖v‖² − ‖S v‖²) / (8 z0)`` with ``v = w·V_s``."""
    v = np.asarray(w, dtype=np.complex128) * complex(v_s)
    b = np.asarray(s2, dtype=np.complex128) @ v
    return float((np.sum(np.abs(v) ** 2) - np.sum(np.abs(b) ** 2)) / (8.0 * float(z0)))


def row5_fractions(p_acc, p_rc30, p_esr_feeds, p_phantom, p_surf, p_rest):
    """Row 5: phantom, copper and ESR fractions of ``P_acc`` (plus the rest)."""
    out = {"phantom": p_phantom / p_acc, "copper": p_surf / p_acc,
           "esr": (p_rc30 + p_esr_feeds) / p_acc, "rest": p_rest / p_acc}
    out["sum"] = out["phantom"] + out["copper"] + out["esr"] + out["rest"]
    return out


# -- §9 item 576 (§10 H30a): the v3 record route -- no builder, no field ----------
#: Rows the node record cannot carry (no field on the probes or the feed sheets).
V3_NOT_IN_RECORD = ("row3", "row6", "cs_s")


def _complex_vec(pairs):
    return np.array([complex(a, b) for a, b in pairs], dtype=np.complex128)


def record_port_state(s_node, w, z_th, v_s):
    """Feed voltages and currents under the Thevenin drive ``w`` from the node S.

    ``s_node`` is the current route's ``S_ij = 2V_i/V_s − δ_ij`` (unit drive of
    feed j through ``Z_th``), so by linearity ``V_port = ½ V_s (𝟙 + S) w`` and
    ``I = (w − V_port) / Z_th``.  Returns ``(V_port, I_feed)``.
    """
    s = np.asarray(s_node, dtype=np.complex128)
    w = np.asarray(w, dtype=np.complex128)
    v = 0.5 * complex(v_s) * ((np.eye(s.shape[0]) + s) @ w)
    i = (w - v) / complex(z_th)
    return v, i


def v3_record_rows(rec, sar10g_cstar_w_per_kg=None):
    """SPEC v3 rows 4–5, the gap and X_s / ΣX_s / P_acc from a node record (ccw).

    ``_stage_v3``'s formulas (power and row 5, the mechanism's X_s / terminal /
    Norton terms) with the record route's V / I in place of the field currents.
    Rows 3 / 6 and CS_s need a field and are named in ``not_in_record``.
    """
    from tests.validation.test_port_package_sparameters import REFERENCE_IMPEDANCE_OHM
    from tests.validation.test_wf10_f_human_copper_esr_root import (
        ESR_OHM,
        _source_voltage,
        _thevenin_feed_impedance,
        _z_rc,
    )

    assert rec["r_ohm"] == ESR_OHM, rec["r_ohm"]
    feeds = list(rec["feeds"])
    assert feeds == list(rec["order"][:2]), (feeds, rec["order"][:2])
    z0 = float(REFERENCE_IMPEDANCE_OHM)
    z_rc = _z_rc(rec["c_f"], rec["freq_hz"], rec["r_ohm"])
    z_th = _thevenin_feed_impedance(z_rc, z0)
    v_s = _source_voltage(z_rc, z0)
    s_node = np.array([[complex(a, b) for a, b in row] for row in rec["s_node"]],
                      dtype=np.complex128)
    assert s_node.shape == (2, 2), s_node.shape
    w = _complex_vec(rec["w"]["ccw"])
    v_port, i_feed = record_port_state(s_node, w, z_th, v_s)
    p_acc, v_line, _i_line = p_acc_line(w, i_feed, z_th, z_rc)
    p_wave = p_acc_wave(s_node, w, v_s, z0)
    part = rec["partitions"]["ccw"]
    by_id = part["p_sheet_by_id"]
    p_rc30 = float(sum(v for pid, v in by_id.items() if pid not in feeds))
    p_esr_feeds = float(sum(0.5 * rec["r_ohm"] * abs(v / z_rc) ** 2 for v in v_line))
    p_ph, p_surf, p_rest, p_src = (float(part[k]) for k in ("p_phantom", "p_surf", "p_rest", "p_src"))
    p_361 = p_src - sum(float(by_id[pid]) for pid in feeds)
    p_sum = p_rc30 + p_esr_feeds + p_ph + p_surf + p_rest
    r_acc1 = abs(p_acc - p_wave) / abs(p_wave)
    r_acc2 = abs(p_acc - p_sum) / abs(p_acc)
    g_th = 0.5 * float(np.real(1.0 / z_th))
    mech, x_sum = {}, 0.0
    for k, pid in enumerate(feeds):
        p_field = float(by_id[pid])
        lumped = g_th * abs(v_port[k]) ** 2
        mech[pid] = {"x_w": p_field - lumped, "lumped_norton_w": lumped, "p_sheet_field_w": p_field}
        x_sum += p_field - lumped
    terminal = float(0.5 * np.real(np.sum(v_port * np.conj(i_feed))))
    norton = terminal + sum(g_th * abs(v) ** 2 for v in v_port)
    out = {
        "degree": rec.get("degree"), "c_in_pf": rec["c_in_pf"], "node_knob": rec.get("node_knob"),
        "feeds": feeds, "v_port": v_port, "i_feed": i_feed,
        "power_w": {"p_acc_spec": p_acc, "p_acc_wave": p_wave, "p_acc_item361": p_361,
                    "p_src": p_src, "p_rc30": p_rc30, "p_esr_feeds": p_esr_feeds,
                    "p_phantom": p_ph, "p_surf": p_surf, "p_rest": p_rest},
        "p_sum": p_sum, "r_acc1": r_acc1, "r_acc2": r_acc2,
        "row5": row5_fractions(p_acc, p_rc30, p_esr_feeds, p_ph, p_surf, p_rest),
        "row5_partition": row5_fractions(p_sum, p_rc30, p_esr_feeds, p_ph, p_surf, p_rest),
        "gap": {"line_cut_minus_partition_w": p_acc - p_sum, "mechanism": mech,
                "terminal_w": terminal, "norton_lumped_w": norton},
        "sum_x_w": x_sum, "sum_x_over_p_acc": x_sum / p_acc,
        "terminal_w": terminal, "norton_lumped_w": norton,
        "negative_controls": {
            "item361_p_acc_in_p_acc_1": abs(p_361 - p_wave) / abs(p_wave),
            "esr_blind": abs(p_acc - (p_ph + p_surf + p_rest)) / abs(p_acc),
            "feed_esr_blind": abs(p_acc - (p_sum - p_esr_feeds)) / abs(p_acc)},
        "not_in_record": list(V3_NOT_IN_RECORD),
    }
    if sar10g_cstar_w_per_kg is not None:
        sar = float(sar10g_cstar_w_per_kg)
        out["row4"] = {"sar10g_w_per_kg": sar, "sar10g_per_w": sar / p_acc,
                       "sar10g_per_w_partition": sar / p_sum}
    return out


_LOG_NUM = r"([-+]?\d+\.\d+(?:e[-+]?\d+)?)"


def _read_log(path):
    """A harness log's text, from ``path`` or, once retention has gzipped it,
    ``path`` + ``.gz`` (`OPS-75`; an example cannot import `tests/log_records.py`)."""
    import gzip

    path = str(path)
    if not os.path.exists(path) and os.path.exists(f"{path}.gz"):
        path = f"{path}.gz"
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt", errors="replace") as fh:
        return fh.read()


def node_log_readings(log_path, node_basename):
    """A node window's harness log: the last ``wrote …/<basename>`` line, the last
    ccw ``(P_w)`` line before it, and every ``(Lin)`` reading (both senses)."""
    import re

    wrote_re = re.compile(r"^\[WF-10 step2\] wrote (\S*/" + re.escape(node_basename) + r")\s*$")
    pw_re = re.compile(r"^\[WF-10 step2\] \[ccw\] w = .*\(P_w\) P_src " + _LOG_NUM + r" sumP_sheet "
                       + _LOG_NUM + r" P_phantom " + _LOG_NUM + r" P_rest " + _LOG_NUM + r" P_surf "
                       + _LOG_NUM + r" W")
    lin_re = re.compile(r"^\[WF-10 step2\] \[(ccw|cw)\] \(Lin\) currents rel " + _LOG_NUM
                        + r" .*Y_in current route vs _active_y_in " + _LOG_NUM)
    lines = _read_log(log_path).splitlines()
    wrote = [(k, wrote_re.match(ln)) for k, ln in enumerate(lines)]
    wrote = [(k, m) for k, m in wrote if m]
    assert wrote, f"{log_path}: no '[WF-10 step2] wrote …/{node_basename}' line"
    k_w, m_w = wrote[-1]
    pw = [pw_re.match(ln) for ln in lines[:k_w]]
    pw = [m for m in pw if m]
    assert pw, f"{log_path}: no ccw (P_w) line before line {k_w + 1}"
    p = [float(x) for x in pw[-1].groups()]
    lin = [(m.group(1), float(m.group(2)), float(m.group(3)))
           for m in (lin_re.match(ln) for ln in lines) if m]
    return {"wrote_path": m_w.group(1), "wrote_line": k_w + 1, "n_wrote": len(wrote),
            "p_w": {"p_src": p[0], "p_sheets": p[1], "p_phantom": p[2], "p_rest": p[3],
                    "p_surf": p[4]},
            "lin": lin}


def sheet_variance_split(area, j_drive, l_tan, t_drive, y_s):
    """Amendment 3's split of a sheet's field loss above its lumped value.

    ``area`` = ∫1 dA, ``j_drive`` = ∫E·ĥ dA, ``l_tan`` = ∫|n̂×E|² dA,
    ``t_drive`` = ∫|E·ĥ|² dA.  Returns ``(cs, transverse, variance)`` with
    ``cs = ½ Re(Y_s)(L − |J|²/A)`` = ``transverse + variance``, where
    ``transverse = ½ Re(Y_s)(L − T)`` (the in-plane component across ĥ) and
    ``variance = ½ Re(Y_s)(T − |J|²/A)`` (the drive component's spread).
    """
    g = 0.5 * float(np.real(y_s))
    a = float(np.real(area))
    mean_sq = abs(complex(j_drive)) ** 2 / a
    transverse = g * (float(np.real(l_tan)) - float(np.real(t_drive)))
    variance = g * (float(np.real(t_drive)) - mean_sq)
    return g * (float(np.real(l_tan)) - mean_sq), transverse, variance


def ring_spectra(currents, port_map):
    """``{ring: k_spectrum}`` of 32 currents placed by ``port_map`` ``[(slot, ring)]``."""
    rm = _ring_modes()
    out = {}
    for ring in rm.RINGS:
        vec = np.zeros(rm.N_SLOTS, dtype=np.complex128)
        for idx, (j, r) in enumerate(port_map):
            if r == ring:
                vec[int(j)] = currents[idx]
        out[ring] = rm.k_spectrum(vec)
    return out


def _stage_v3(comm, outdir):
    """SPEC v3 rows 1–6, our side, on H9′'s tuned copper F-human (degree 1)."""
    from mpi4py import MPI

    import ufl
    from dolfinx import default_scalar_type, fem

    from fem_em_solver.ports.circuit import s_to_z
    from fem_em_solver.ports.lumped import _sheet_measure
    from fem_em_solver.post.faraday import magnetic_flux_density_from_e, project_to_cg1
    from fem_em_solver.post.sar import averaging_ball_radius, build_density_field, mass_averaged_sar
    from tests.mesh.helpers import global_cell_tag_set
    from tests.mesh.test_birdcage_port_sheet_prerequisite import CELL_COUNT_BAND
    from tests.validation.test_birdcage_b1_plus_map import (
        PHANTOM_RHO_KG_PER_M3,
        _rotate_z,
        _sample_points,
    )
    from tests.validation.test_birdcage_b1_quadrature import _mirror_xy, _read_senses
    from tests.validation.test_birdcage_power_identity import DISCRETE_IDENTITY_RTOL, _reduce
    from tests.validation.test_mass_averaged_sar_standard_masses import (
        QUADRATURE_DEGREE,
        TEN_GRAM_KG,
    )
    from tests.validation.test_port_birdcage_lumped_column import PHANTOM_CELL_TAG
    from tests.validation.test_port_birdcage_ring_matrix import REUSE_REPRODUCTION_RTOL
    from tests.validation.test_port_circuit_layer_field import _residual
    from tests.validation.test_port_package_sparameters import REFERENCE_IMPEDANCE_OHM
    from tests.validation.test_wf7_f_human_b1_quadrature import FREQUENCY_HZ
    from tests.validation.test_wf7_f_human_sar10g_hotspot import (
        CONTAINMENT_MARGIN_M,
        _phantom_geometry,
    )
    from tests.validation.test_wf10_f_human_copper_esr_root import (
        ESR_OHM,
        NPZ_PATH,
        _parallel_currents,
        _parallel_loads,
        _source_voltage,
        _thevenin_feed_impedance,
        _z_rc,
    )
    from tests.validation.test_wf10_f_human_tuned_b1_sar import (
        C_IN_F_HUMAN_PF,
        DEGREE_ENV,
        NODE_JSON,
        build_tuned_quadrature,
    )
    from tests.validation.test_wf10_f_human_tuned_circuit import _respec
    from tests.validation.test_wf10_f_human_tuned_sar10g import _containment, _lattice
    from tests.validation.test_wf12_f_human_copper_matrix import BLIND_FACTOR, HOLE_CELL_RECORD

    t_all = time.perf_counter()
    rank0 = comm.rank == 0
    fails = []

    def say(msg):
        if rank0:
            print(f"[ANS-8 v3] {msg}", flush=True)

    def check(name, ok, msg):
        say(f"{name}: {'PASS' if ok else 'FAIL'} — {msg}")
        if not ok:
            fails.append(name)

    assert os.environ.get(DEGREE_ENV, "").strip() == "", f"{DEGREE_ENV} must be unset (degree 1)"
    say(f"start: -n {comm.size}; predictions (printed only): {V3_PREDICTED}")

    # 1. The node record, read BEFORE the builder (which rewrites it).
    for path in (NODE_JSON, NPZ_PATH):
        assert os.path.exists(path), f"{path} absent: item 369 stops (BLOCKED); never re-solve"
    rec = None
    if rank0:
        with open(NODE_JSON) as fh:
            rec = json.load(fh)
    rec = comm.bcast(rec, root=0)
    assert rec["c_in_pf"] == C_IN_F_HUMAN_PF, rec["c_in_pf"]
    s_rec = np.array([[complex(a, b) for a, b in row] for row in rec["s_node"]], dtype=np.complex128)
    y_rec = [complex(a, b) for a, b in rec["y_in"]["ccw"]]

    # 2. The builder.
    out = build_tuned_quadrature(comm, with_z50=False)
    ctx, order, feeds, times = out["ctx"], out["order"], out["feeds"], out["times"]
    msh, cell_tags, sigma, omega = ctx["msh"], ctx["cell_tags"], out["sigma"], out["omega"]
    t_build = time.perf_counter() - t_all
    s_node = np.asarray(out["s_node"], dtype=np.complex128)
    z0 = float(REFERENCE_IMPEDANCE_OHM)
    c_f = C_IN_F_HUMAN_PF * 1e-12
    z_rc = _z_rc(c_f, FREQUENCY_HZ)
    z_th = _thevenin_feed_impedance(z_rc, z0)
    v_s = _source_voltage(z_rc, z0)
    w = np.array(out["w"]["ccw"], dtype=np.complex128)

    r_srep = _residual(s_node, s_rec)
    check("(S_rep)", r_srep <= REUSE_REPRODUCTION_RTOL,
          f"s_node vs {os.path.basename(NODE_JSON)} read before the builder: {r_srep:.3e} (band "
          f"{REUSE_REPRODUCTION_RTOL:g}; 359 read 6.027e-13, 361 2.705e-13)")
    anc = out["anchors"]
    cell_rel = out["cells"] / HOLE_CELL_RECORD - 1.0
    check("(A0)", abs(cell_rel) <= CELL_COUNT_BAND,
          f"{out['cells']} cells vs {HOLE_CELL_RECORD}, rel {cell_rel:.3e} (band {CELL_COUNT_BAND})")
    check("(P_RC)", anc["p_rc_node"] <= DISCRETE_IDENTITY_RTOL,
          f"{anc['p_rc_node']:.3e} (band {DISCRETE_IDENTITY_RTOL:g})")
    check("NEG node (asserted)", anc["blind_node"] >= BLIND_FACTOR * DISCRETE_IDENTITY_RTOL,
          f"P_surf-blind miss {anc['blind_node']:.4e} (bar {BLIND_FACTOR:.0f} x "
          f"{DISCRETE_IDENTITY_RTOL:g}; 359 / 361 read 7.0125e-3)")
    rel_w, blind_w = anc["p_w"]["ccw"]
    check("(P_w) ccw", rel_w <= DISCRETE_IDENTITY_RTOL, f"{rel_w:.3e} (band {DISCRETE_IDENTITY_RTOL:g})")
    check("NEG (P_w) ccw P_surf-blind (asserted)", blind_w >= BLIND_FACTOR * DISCRETE_IDENTITY_RTOL,
          f"{blind_w:.4e} (bar {BLIND_FACTOR:.0f}x; 358 read 7.0397e-3)")

    # 3. Row 1.
    y_ccw = [complex(y) for y in out["y_in"]["ccw"]]
    r_y = max(abs(a - b) for a, b in zip(y_ccw, y_rec)) / max(abs(b) for b in y_rec)
    for n, a, b in zip(feeds, y_ccw, y_rec):
        say(f"row 1: Y_in ccw {n}: {a.real:+.9e} {a.imag:+.9e}j S (JSON {b.real:+.9e} "
            f"{b.imag:+.9e}j)")
    check("(Rep-Y)", r_y <= REUSE_REPRODUCTION_RTOL,
          f"{r_y:.3e} relative to |Y_in| (band {REUSE_REPRODUCTION_RTOL:g})")

    # 4. Row 2: field route (ours) and the circuit route (printed beside).
    port_map = _ring_modes().ring_port_map_from_names(order)
    measured = [(int(out["slots"][p][0]), "top" if float(out["slots"][p][1]) > 0 else "bottom")
                for p in order]
    check("(G2-v3) names map = measured slot map", port_map == measured,
          f"ring_port_map_from_names(order) equals the builder's slots: {port_map == measured}")
    row2 = ring_spectra(np.asarray(out["cur"]["ccw"]), port_map)
    with np.load(NPZ_PATH, allow_pickle=False) as npz:
        s32_cu = np.array(npz["s32"], dtype=np.complex128)
        order_cu = [str(p) for p in npz["order"]]
    assert order_cu == order, "npz order differs"
    z32 = s_to_z(s32_cu, z0)
    i_circ = _parallel_currents(z32, _parallel_loads(32, [0, 1], z_rc, z_th),
                                np.concatenate([w, np.zeros(30, dtype=np.complex128)]))
    row2_circ = ring_spectra(i_circ, port_map)
    for ring in ("bottom", "top"):
        f, c = row2[ring], row2_circ[ring]
        say(f"row 2 ({ring}{', feed ring' if ring == 'bottom' else ''}): share_k1 "
            f"{f['share_k1']:.6f}, leak_k0 {f['leak_k0']:.6e}, |k| folds "
            + " ".join(f"{x:.4f}" for x in f["abs_k_folds"])
            + f"  || circuit-only (h12c npz; H12c identity red at 1.58; not our row): share_k1 "
            f"{c['share_k1']:.6f}, leak_k0 {c['leak_k0']:.6e}")

    # 5. Power.
    cur = np.asarray(out["cur"]["ccw"], dtype=np.complex128)
    p_acc, v_port, _i_line = p_acc_line(w, cur[:2], z_th, z_rc)
    p_wave = p_acc_wave(s_node, w, v_s, z0)
    part = out["parts"]["ccw"]
    p_rc30 = float(sum(v for pid, v in part["p_sheet_by_id"].items() if pid not in feeds))
    p_esr_feeds = float(sum(0.5 * ESR_OHM * abs(v / z_rc) ** 2 for v in v_port))
    p_ph, p_surf, p_rest, p_src = part["p_phantom"], part["p_surf"], part["p_rest"], part["p_src"]
    p_361 = p_src - sum(part["p_sheet_by_id"][pid] for pid in feeds)
    r_acc1 = abs(p_acc - p_wave) / abs(p_wave)
    p_sum = p_rc30 + p_esr_feeds + p_ph + p_surf + p_rest
    r_acc2 = abs(p_acc - p_sum) / abs(p_acc)
    fr = row5_fractions(p_acc, p_rc30, p_esr_feeds, p_ph, p_surf, p_rest)
    say(f"power: P_acc,SPEC {p_acc:.12e} W (PREDICTED ≈ {V3_PREDICTED['p_acc_spec_w']:g}); wave form "
        f"{p_wave:.12e} W; item 361's P_acc {p_361:.12e} W; P_src {p_src:.12e} W")
    say(f"partition: P_RC,30 {p_rc30:.9e}, P_ESR,feeds {p_esr_feeds:.9e}, P_phantom {p_ph:.9e}, "
        f"P_surf {p_surf:.9e}, P_rest {p_rest:.9e} W")
    check("(P_acc-1) line cut = wave form", r_acc1 <= DISCRETE_IDENTITY_RTOL,
          f"{r_acc1:.3e} (band {DISCRETE_IDENTITY_RTOL:g}; PREDICTED <= 1e-12)")
    # (P_acc-2) withdrawn (10:30 review ruling (b), §9 item 386 amendment 1):
    # its two numbers are printed as the gap, never asserted.
    gap = p_acc - p_sum
    say(f"gap (printed, (P_acc-2) withdrawn): line cut − partition = {gap:.9e} W; "
        f"{r_acc2:.4e} of the line cut, {abs(gap) / abs(p_361):.4e} of the field net "
        f"(item 361's P_acc) — feed-sheet field-vs-terminal gap, degree 1 (known-issues 2026-10-11)")
    r_rep_p = abs(p_acc / P_ACC_SPEC_H21_RECORD_W - 1.0)
    check("(Rep-P)", r_rep_p <= REUSE_REPRODUCTION_RTOL,
          f"P_acc,SPEC {p_acc:.12e} vs item 369's {P_ACC_SPEC_H21_RECORD_W:.12e} W: rel {r_rep_p:.3e} "
          f"(band {REUSE_REPRODUCTION_RTOL:g})")
    part_now = {"p_rc30": p_rc30, "p_esr_feeds": p_esr_feeds, "p_phantom": p_ph, "p_surf": p_surf}
    r_rep_part = {k: abs(part_now[k] / v - 1.0) for k, v in PARTITION_H21_RECORD_W.items()}
    check("(Rep-Part)", max(r_rep_part.values()) <= REUSE_REPRODUCTION_RTOL,
          ", ".join(f"{k} rel {v:.3e}" for k, v in r_rep_part.items())
          + f" vs item 369's terms (band {REUSE_REPRODUCTION_RTOL:g}; 10-digit print <= 1e-9)")
    nc_esr = abs(p_acc - (p_ph + p_surf + p_rest)) / abs(p_acc)
    nc_feed = abs(p_acc - (p_sum - p_esr_feeds)) / abs(p_acc)
    nc_361 = abs(p_361 - p_wave) / abs(p_wave)
    check("NC item 361's P_acc in (P_acc-1) (asserted)", nc_361 >= BLIND_FACTOR * DISCRETE_IDENTITY_RTOL,
          f"miss {nc_361:.4e} (bar {BLIND_FACTOR:.0f} x {DISCRETE_IDENTITY_RTOL:g}; item 369 read "
          f"4.9700e-02)")
    say(f"NC (printed): ESR-blind miss {nc_esr:.4e} (item 369 read 0.71196); feed-ESR-blind "
        f"{nc_feed:.4e}")
    fr_part = row5_fractions(p_sum, p_rc30, p_esr_feeds, p_ph, p_surf, p_rest)
    say(f"row 5 (of P_acc,SPEC): phantom {fr['phantom']:.6f}, copper {fr['copper']:.6f}, ESR "
        f"{fr['esr']:.6f}, rest {fr['rest']:.3e}; sum {fr['sum']:.6f}, residual {1.0 - fr['sum']:.6f} "
        f"= feed-sheet field-vs-terminal gap, degree 1 (known-issues 2026-10-11)")
    say(f"row 5 (of the field partition, printed beside): phantom {fr_part['phantom']:.6f}, copper "
        f"{fr_part['copper']:.6f}, ESR {fr_part['esr']:.6f}; sum {fr_part['sum']:.15f}")

    # Amendment 3: the gap's mechanism on the two feed sheets — printed,
    # pre-registered, never asserted.  Every rank enters every _reduce.
    msh_f, tags_f = ctx["msh"], ctx["tags_f"]
    specs_f = {sp.port_id: sp for sp in _respec(list(ctx["specs"]), None,
                                                {pid: z_th for pid in feeds})}
    e_ccw = out["e"]["ccw"]
    one = fem.Constant(msh_f, default_scalar_type(1.0))
    n_hat = ufl.FacetNormal(msh_f)("+")
    e_r = e_ccw("+")
    mech, x_sum = {}, 0.0
    t0 = time.perf_counter()
    for i, pid in enumerate(feeds):
        sh = specs_f[pid].sheet(driven=False)
        assert sh.interior, pid
        ds_s = _sheet_measure(msh_f, tags_f, sh)
        h_hat = ufl.as_vector([default_scalar_type(c) for c in sh.unit_drive()])
        e_h = ufl.dot(e_r, h_hat)
        nxe = ufl.cross(n_hat, e_r)
        a_s = _reduce(comm, one * ds_s)
        j_s = _reduce(comm, e_h * ds_s)
        l_s = _reduce(comm, ufl.inner(nxe, nxe) * ds_s)
        t_s = _reduce(comm, ufl.inner(e_h, e_h) * ds_s)
        y_s = 1.0 / complex(sh.sheet_resistivity)
        cs, trans, var = sheet_variance_split(a_s, j_s, l_s, t_s, y_s)
        p_field = float(part["p_sheet_by_id"][pid])
        lumped = 0.5 * float(np.real(1.0 / z_th)) * abs(v_port[i]) ** 2
        x_s = p_field - lumped
        x_sum += x_s
        p_l = 0.5 * float(np.real(y_s)) * float(np.real(l_s))
        mech[pid] = {"area_m2": float(np.real(a_s)), "j_drive": [float(np.real(j_s)), float(np.imag(j_s))],
                     "l_tan": float(np.real(l_s)), "t_drive": float(np.real(t_s)), "cs_w": cs,
                     "transverse_w": trans, "variance_w": var, "x_w": x_s, "cs_over_x": cs / x_s,
                     "p_sheet_field_w": p_field, "half_re_y_l_w": p_l, "lumped_norton_w": lumped}
        say(f"mechanism {pid}: A_s {np.real(a_s):.9e} m^2 (w*h {sh.sheet_width_m * sh.gap_height_m:.9e}); "
            f"J_s {complex(j_s):.9e}; L_s {np.real(l_s):.9e}; T_s {np.real(t_s):.9e}")
        say(f"mechanism {pid}: CS_s {cs:.9e} W vs X_s {x_s:.9e} W (PREDICTED "
            f"{V3_GAP_X_PREDICTED_W.get(pid, float('nan')):.5e}): CS_s/X_s {cs / x_s:.9f} (window "
            f"{V3_CS_OVER_X_WINDOW}, printed, never asserted); split transverse {trans:.9e} W + "
            f"drive-variance {var:.9e} W; ½Re(Y_s)L_s {p_l:.9e} vs P_sheet,field {p_field:.9e} W")
    terminal = float(0.5 * np.real(np.sum(v_port * np.conj(cur[:2]))))
    norton = terminal + sum(0.5 * float(np.real(1.0 / z_th)) * abs(v) ** 2 for v in v_port)
    say(f"mechanism: P_src {p_src:.12e} W vs Norton lumped (terminal {terminal:.12e} + Σ½Re(1/Z_th)|V|²) "
        f"{norton:.12e} W: rel {abs(p_src / norton - 1.0):.3e} (PREDICTED <= 1e-12, printed); "
        f"ΣX_s {x_sum:.9e} W vs field net − terminal {p_361 - terminal:.9e} W; "
        f"{time.perf_counter() - t0:.2f} s")

    # 6. Rows 3 and 6.
    comm.Barrier()
    t0 = time.perf_counter()
    cg1 = {}
    for sense in ("ccw", "cw"):
        b_dg = magnetic_flux_density_from_e(out["e"][sense], omega, name=f"B_{sense}")
        cg1[sense] = project_to_cg1(b_dg, name=f"B_{sense}_cg1")
        del b_dg
    centre, radius, half_height, _lo, _hi = _phantom_geometry(msh, cell_tags, PHANTOM_CELL_TAG, comm)
    z_c = float(centre[2])
    probes = np.array([[x, y, z_c] for _n, x, y in V3_PROBES_XY], dtype=float)
    pp, pm, pv = _read_senses(cg1["ccw"], probes)
    # (Rep-B): 358's masked centre sample, mask over the same images.
    points = _sample_points({"mesh": msh, "cell_tags": cell_tags})
    az_m, az17 = float(out["az_mirror"]), float(out["az_feeds"][0])
    ccw_x = _read_senses(cg1["ccw"], points)
    cw_x = _read_senses(cg1["cw"], points)
    cw_mx = _read_senses(cg1["cw"], _mirror_xy(points, az_m))
    cw_m17 = _read_senses(cg1["cw"], _mirror_xy(points, az17))
    rot = [_read_senses(cg1["ccw"], _rotate_z(points, np.radians(m * 22.5))) for m in range(1, 16)]
    mask = np.logical_and.reduce([ccw_x[2], cw_x[2], cw_mx[2], cw_m17[2]] + [r[2] for r in rot])
    mean_b = float(np.mean(ccw_x[0][mask]))
    comm.Barrier()
    t_b = time.perf_counter() - t0
    r_b = abs(mean_b / MEAN_B1P_CCW_H9P_RECORD_T - 1.0)
    check("(Rep-B)", r_b <= REUSE_REPRODUCTION_RTOL,
          f"mean|B1+|_ccw {mean_b:.9e} T over {int(np.count_nonzero(mask))} masked points vs 358's "
          f"{MEAN_B1P_CCW_H9P_RECORD_T:.6e} T: rel {r_b:.3e} (band {REUSE_REPRODUCTION_RTOL:g}; "
          f"7-digit rounding <= 1.8e-7)")
    check("probes valid", bool(np.all(pv)), f"{[n for n, *_x in V3_PROBES_XY]} valid {pv.tolist()}")
    row3 = {}
    for k, (n, x, y) in enumerate(V3_PROBES_XY):
        row3[n] = {"xyz_m": [x, y, z_c], "b1p_t": float(pp[k]), "b1m_t": float(pm[k]),
                   "b1p_per_sqrt_w": float(pp[k] / np.sqrt(p_acc)),
                   "b1p_per_sqrt_w_partition": float(pp[k] / np.sqrt(p_sum))}
        say(f"row 3: {n} ({x:+.3f}, {y:+.3f}, {z_c:+.3e}) m: |B1+| {pp[k]:.6e} T, |B1-| {pm[k]:.6e} "
            f"T; |B1+|/sqrt(P_acc) {row3[n]['b1p_per_sqrt_w']:.6e} T/sqrt(W) (per partition "
            f"{row3[n]['b1p_per_sqrt_w_partition']:.6e})")
    purity = float(pp[0] / pm[0])
    check("(Sense)", purity > 1.0, f"|B1+|/|B1-| at O {purity:.4f} > 1 (row 6; PREDICTED ≈ "
          f"{V3_PREDICTED['purity_o']})")

    # 7. Row 4: 359's 18-point lattice under ccw.
    a10 = averaging_ball_radius(mass_kg=TEN_GRAM_KG, rho=PHANTOM_RHO_KG_PER_M3)
    all_tags = global_cell_tag_set(msh, cell_tags)
    rho_field = build_density_field(
        msh, PHANTOM_RHO_KG_PER_M3, cell_tags=cell_tags,
        density_map={int(t): 0.0 for t in sorted(all_tags) if int(t) != PHANTOM_CELL_TAG})
    comm.Barrier()
    t0 = time.perf_counter()
    lattice = _lattice(az17, z_c)
    sar, cont, coords = {}, [], {}
    for label, _r, _dz, _d, c in lattice:
        c = tuple(float(v) for v in c)
        coords[label] = c
        cont.append(_containment(c, centre, radius, half_height, a10))
        r = mass_averaged_sar(out["e"]["ccw"], sigma=sigma, rho=rho_field, center=c, radius=a10,
                              comm=comm, quadrature_degree=QUADRATURE_DEGREE)
        sar[label] = float(r["averaged_sar_w_per_kg"])
        say(f"  ball {label:<22s} ({c[0] * 1e3:8.3f}, {c[1] * 1e3:8.3f}, {c[2] * 1e3:8.3f}) mm: "
            f"SAR10g {sar[label]:.9e} W/kg, mass {r['mass_kg'] * 1e3:.6f} g")
    comm.Barrier()
    t_sar = time.perf_counter() - t0
    key_star = max(sar, key=sar.get)
    c_star = coords[key_star]
    r_s = abs(sar[V3_CSTAR_LABEL_RECORD] / SAR10G_CSTAR_H9P_RECORD_W_PER_KG - 1.0)
    check("(Rep-S)", r_s <= REUSE_REPRODUCTION_RTOL and key_star == V3_CSTAR_LABEL_RECORD,
          f"SAR10g({V3_CSTAR_LABEL_RECORD}) {sar[V3_CSTAR_LABEL_RECORD]:.9e} vs 359's "
          f"{SAR10G_CSTAR_H9P_RECORD_W_PER_KG:.9e} W/kg: rel {r_s:.3e} (band "
          f"{REUSE_REPRODUCTION_RTOL:g}); c* = {key_star}")
    worst = max(cont)
    check("(iv) containment", worst <= 0.0,
          f"worst {worst * 1e3:.4f} mm over {len(cont)} balls (margin {CONTAINMENT_MARGIN_M * 1e3:g} "
          f"mm inside; <= 0)")
    say(f"row 4: c* (lattice argmax, sphere a = {a10 * 1e3:.4f} mm) = ({c_star[0]:.6f}, "
        f"{c_star[1]:.6f}, {c_star[2]:.6f}) m, SAR10g {sar[key_star]:.9e} W/kg; per W of P_acc "
        f"{sar[key_star] / p_acc:.6e} W/kg/W (PREDICTED ≈ {V3_PREDICTED['sar10g_per_w']}); per W of "
        f"the partition {sar[key_star] / p_sum:.6e}")
    say(f"row 6: |B1+|/|B1-| at O {purity:.4f}")

    rss = _rss_summed_gib(comm)
    say(f"price: builder {t_build:.1f} s (build {times['build']:.1f}, node pair "
        f"{times['node_pair']:.1f}, superposition {times['superposition']:.1f}); B + points {t_b:.1f} s; "
        f"18 balls {t_sar:.1f} s; total {time.perf_counter() - t_all:.1f} s at -n {comm.size}; "
        f"summed RSS {rss:.2f} GiB")

    if fails:
        say(f"RED: {fails} — stop; no v3 rows written")
        return 1
    rows = {
        "spec": V3_SPEC_LABEL, "degree": 1, "c_in_pf": C_IN_F_HUMAN_PF, "r_ohm": ESR_OHM,
        "freq_hz": float(FREQUENCY_HZ), "cells": int(out["cells"]), "drive": "ccw (1, -1j), P17 / P21",
        "feeds": list(feeds),
        "row1_y_in_ccw_s": {n: [y.real, y.imag] for n, y in zip(feeds, y_ccw)},
        "row2_field": row2, "row2_circuit_only": row2_circ,
        "row3": row3, "row4": {"c_star_m": list(c_star), "c_star_label": key_star,
                               "sar10g_w_per_kg": sar[key_star],
                               "sar10g_per_w": sar[key_star] / p_acc,
                               "averaging_volume": "sphere", "radius_m": a10,
                               "lattice_w_per_kg": sar},
        "row5": fr, "row5_partition": fr_part, "row6_purity_o": purity,
        "per_p_acc_spec": {"p_w": p_acc, "row3_b1p_per_sqrt_w": {n: row3[n]["b1p_per_sqrt_w"]
                                                                for n in row3},
                           "row4_sar10g_per_w": sar[key_star] / p_acc, "row5": fr},
        "per_partition": {"p_w": p_sum, "row3_b1p_per_sqrt_w": {
            n: row3[n]["b1p_per_sqrt_w_partition"] for n in row3},
                          "row4_sar10g_per_w": sar[key_star] / p_sum, "row5": fr_part},
        "gap": {"line_cut_minus_partition_w": gap, "rel_line_cut": r_acc2,
                "rel_field_net": abs(gap) / abs(p_361),
                "label": "feed-sheet field-vs-terminal gap, degree 1 (known-issues 2026-10-11)",
                "mechanism": mech, "terminal_w": terminal, "norton_lumped_w": norton},
        "power_w": {"p_acc_spec": p_acc, "p_acc_wave": p_wave, "p_acc_item361": p_361,
                    "p_src": p_src, "p_rc30": p_rc30, "p_esr_feeds": p_esr_feeds,
                    "p_phantom": p_ph, "p_surf": p_surf, "p_rest": p_rest},
        "anchors": {"s_rep": r_srep, "rep_y": r_y, "rep_b": r_b, "rep_s": r_s,
                    "p_acc_1": r_acc1, "rep_p": r_rep_p, "rep_part": r_rep_part,
                    "p_rc": anc["p_rc_node"],
                    "neg_node": anc["blind_node"], "p_w_ccw": rel_w, "neg_p_w_ccw": blind_w,
                    "containment_worst_m": worst},
        "negative_controls_printed": {"esr_blind": nc_esr, "feed_esr_blind": nc_feed,
                                      "item361_p_acc_in_p_acc_1": nc_361},
    }
    if rank0:
        outdir.mkdir(parents=True, exist_ok=True)
        (outdir / V3_JSON_NAME).write_text(json.dumps(rows, indent=1, sort_keys=True) + "\n")
        print(f"[ANS-8 v3] wrote {outdir / V3_JSON_NAME}", flush=True)
    comm.Barrier()
    return 0


def check_v3_power(rows, rtol):
    """Stage R's re-check of (P_acc-1) from the stored scalars; the line-cut −
    partition gap ``r2`` is printed only ((P_acc-2) withdrawn, §9 item 386)."""
    p = rows["power_w"]
    r1 = abs(p["p_acc_spec"] - p["p_acc_wave"]) / abs(p["p_acc_wave"])
    s = p["p_rc30"] + p["p_esr_feeds"] + p["p_phantom"] + p["p_surf"] + p["p_rest"]
    r2 = abs(p["p_acc_spec"] - s) / abs(p["p_acc_spec"])
    return r1, r2, r1 <= rtol


# -- stage R ----------------------------------------------------------------------
def _slot_port_map(saved):
    """``[(slot, ring)]`` from the fixture's measured slots (ring sign −1 = bottom)."""
    return [(int(k), "top" if float(sg) > 0 else "bottom")
            for k, sg in zip(saved["slot_k"], saved["slot_sign"])]


def anchor_rows(s32, port_map, freq_hz, c_mode_f):
    """The anchor computations (A1–A6, negative control), all through ``ring_mode_rows``."""
    rows = _ring_modes().ring_mode_rows
    b = lambda j: (j, "bottom")  # noqa: E731
    two = rows(s32, port_map=port_map, freq_hz=freq_hz, c_f=c_mode_f, feeds=[b(0), b(4)],
               drives={"p17_alone": (1.0, 0.0), "cw": (1.0, 1.0j), "ccw": (1.0, -1.0j)})
    opp = rows(s32, port_map=port_map, freq_hz=freq_hz, c_f=c_mode_f, feeds=[b(0), b(8)],
               drives={"opposite": (1.0, -1.0)})
    four = rows(s32, port_map=port_map, freq_hz=freq_hz, c_f=c_mode_f,
                feeds=[b(0), b(4), b(8), b(12)], drives={"four_cw": (1.0, 1.0j, -1.0, -1.0j)})
    neg = rows(s32, port_map=port_map, freq_hz=freq_hz, c_f=C2_RECORD_F, feeds=[b(0), b(4)],
               drives={"p17_alone": (1.0, 0.0)})
    return two, opp, four, neg


def _shares(block):
    return block["rings"]["bottom"]["share_k1"], block["rings"]["top"]["share_k1"]


def spec_controls(port_map, names):
    """SPEC_CONTROLS as ``{label: ([(slot, ring)], voltages)}`` through ``port_map``."""
    names = list(names)
    return {label: ([port_map[names.index(n)] for n in pnames], volts)
            for label, (pnames, volts) in SPEC_CONTROLS.items()}


def rows_v2(s32, port_map, names, freq_hz, c_spec_f=SPEC_C_F):
    """SPEC v2 rows through the shared ``ring_mode_rows_v2`` (both codes)."""
    names = list(names)
    feeds = [port_map[names.index(n)] for n in SPEC_FEED_NAMES]
    return _ring_modes().ring_mode_rows_v2(
        s32, port_map=port_map, freq_hz=freq_hz, c_spec_f=c_spec_f, feeds=feeds,
        drives=SPEC_DRIVES, controls=spec_controls(port_map, names))


def _even_max(block):
    return max(block["rings"][r]["abs_k_folds"][k] for r in ("bottom", "top") for k in (0, 2, 4, 6, 8))


def _leak_max(block):
    return max(block["rings"][r]["leak_k0"] for r in ("bottom", "top"))


def _stage_r(outdir):
    from tests.mesh.test_birdcage_port_sheet_prerequisite import CELL_COUNT_BAND
    from tests.validation.test_birdcage_f_human_rung import F_HUMAN_BRANCH_B_CELL_RECORD
    from tests.validation.test_port_birdcage_four_port import LEG_D0_REPRODUCTION_BAND
    from tests.validation.test_port_birdcage_ring_column import _reciprocity_ratio
    from tests.validation.test_port_birdcage_ring_matrix import MATRIX_PASSIVITY_CEILING
    from tests.validation.test_port_lumped_sheet_sweep import RECIPROCITY_BAND
    from tests.validation.test_wf10_f_human_tuned_circuit import C_MODE_F, STEP1C_NPZ

    t0 = time.perf_counter()
    fails = []

    def say(msg):
        print(f"[ANS-8 R] {msg}", flush=True)

    def check(label, ok, msg):
        say(f"{label} {msg} — {'green' if ok else 'RED'}")
        if not ok:
            fails.append(label)

    path = outdir / NPZ_NAME
    if not path.exists():
        say(f"{path} missing: stage R never recomputes the 32x32 — run stage S first; stop")
        return 1
    saved = np.load(path, allow_pickle=False)
    s32 = np.array(saved["s32"], dtype=np.complex128)
    order = [str(p) for p in saved["order"]]
    freq_hz = float(saved["freq_hz"])
    kappa = float(saved["kappa"])
    n_cells = int(saved["cells"])
    say(f"loaded {path}: order {order[:3]}…, f = {freq_hz:.3e} Hz, degree {int(saved['degree'])}, "
        f"kappa {kappa:.6e} beside item 133's {KAPPA_RECORD:.6e} (printed)")

    # (G1) re-checked from the npz.
    rel_cells = n_cells / F_HUMAN_BRANCH_B_CELL_RECORD - 1.0
    check("(G1) cells", abs(rel_cells) < CELL_COUNT_BAND,
          f"{n_cells} relative {rel_cells:.3e} (band {CELL_COUNT_BAND})")
    rec = _reciprocity_ratio(s32)
    sig = float(np.linalg.norm(s32, 2))
    check("(G1) reciprocity", rec <= RECIPROCITY_BAND, f"{rec:.6e} (band {RECIPROCITY_BAND:.0e})")
    check("(G1) sigma_max", sig <= MATRIX_PASSIVITY_CEILING,
          f"{sig:.9f} (ceiling {MATRIX_PASSIVITY_CEILING})")

    # (G2) measured slot map == SPEC naming map.
    rm = _ring_modes()
    measured = _slot_port_map(saved)
    named = rm.ring_port_map_from_names(order)
    mism = [(p, m, n) for p, m, n in zip(order, measured, named) if m != n]
    check("(G2) slot map", not mism, f"measured vs SPEC naming: {len(mism)} mismatches {mism[:3]}")
    if mism:
        say("(G2) red: an operator finding (port naming) — stop before rows")
        return 1
    port_map = named

    # (G3) item 128's held npz, if present.
    if os.path.exists(STEP1C_NPZ):
        held = np.load(STEP1C_NPZ, allow_pickle=False)
        held_order = [str(p) for p in held["order"]]
        same = held_order == order
        hs = np.array(held["s32"], dtype=np.complex128)
        g3 = float(np.max(np.abs(s32 - hs) / np.abs(hs))) if same else float("inf")
        check("(G3) held npz", same and g3 <= LEG_D0_REPRODUCTION_BAND,
              f"{STEP1C_NPZ}: same order {same}; entrywise {g3:.3e} (band "
              f"{LEG_D0_REPRODUCTION_BAND:.0e}; item 128 A3 3.804e-13)")
        if not (same and g3 <= LEG_D0_REPRODUCTION_BAND):
            say("(G3) red: the rebuild differs from item 128's matrix — stop before rows")
            return 1
    else:
        say(f"(G3) {STEP1C_NPZ} absent — skipped (printed)")

    # Anchors at C_MODE_F, cw.
    two, opp, four, neg = anchor_rows(s32, port_map, freq_hz, C_MODE_F)
    a1 = abs(two["c_mode_f"] / C_MODE_F - 1.0)
    check("(A1) C_mode", a1 <= LEG_D0_REPRODUCTION_BAND,
          f"{two['c_mode_f']:.12e} F vs {C_MODE_F:.12e}: relative {a1:.3e} "
          f"(band {LEG_D0_REPRODUCTION_BAND:.0e})")
    a2 = abs(two["c_k2_f"] / C2_RECORD_F - 1.0)
    check("(A2) C_2", a2 <= LEG_D0_REPRODUCTION_BAND,
          f"{two['c_k2_f']:.12e} F vs {C2_RECORD_F:.12e}: relative {a2:.3e} "
          f"(band {LEG_D0_REPRODUCTION_BAND:.0e})")
    a3 = abs(two["off_block"] / OFF_BLOCK_RECORD - 1.0)
    check("(A3) off-block", a3 <= OFF_BLOCK_RTOL,
          f"{two['off_block']:.6e} vs {OFF_BLOCK_RECORD:.6e}: relative {a3:.3e} (band {OFF_BLOCK_RTOL:.0e})")

    def shares_and_folds(label, block, rec_shares, rec_folds):
        sh = _shares(block)
        d = max(abs(sh[0] - rec_shares[0]), abs(sh[1] - rec_shares[1]))
        check(f"{label} shares", d <= SHARE_ATOL,
              f"feed ring {sh[0]:.6f} other {sh[1]:.6f} vs {rec_shares[0]:.6f} / "
              f"{rec_shares[1]:.6f}: max |d| {d:.2e} (band {SHARE_ATOL:.0e})")
        if rec_folds is not None:
            for ring in ("bottom", "top"):
                folds = block["rings"][ring]["abs_k_folds"]
                df = max(abs(a - b) for a, b in zip(folds, rec_folds[ring]))
                check(f"{label} |k|-folds {ring}", df <= FOLD_ATOL,
                      " ".join(f"{x:.4f}" for x in folds) + f": max |d| {df:.2e} (band {FOLD_ATOL:.0e})")

    shares_and_folds("(A4) cw two-port", two["drives"]["cw"], A4_CW_SHARES, A4_CW_FOLDS)
    shares_and_folds("(A5) P17 alone", two["drives"]["p17_alone"], A5_P17_SHARES, A5_P17_FOLDS)
    shares_and_folds("(A6) opposite pair", opp["drives"]["opposite"], A6_OPPOSITE_SHARES, None)
    shares_and_folds("(A6) four-port", four["drives"]["four_cw"], A6_FOUR_PORT_SHARES, None)
    for label, blk in (("opposite pair", opp["drives"]["opposite"]),
                       ("four-port", four["drives"]["four_cw"])):
        even = max(blk["rings"][r]["abs_k_folds"][k] for r in ("bottom", "top") for k in (0, 2, 4, 6, 8))
        check(f"(A6) {label} even-|k|", even <= FOLD_ATOL / 2, f"max {even:.2e} (printed 0.0000)")

    nc = neg["drives"]["p17_alone"]["rings"]["bottom"]["share_k1"]
    check("(NC) C_2 P17-alone feed-ring share", nc < NEGATIVE_CONTROL_CEILING,
          f"{nc:.6f} (< {NEGATIVE_CONTROL_CEILING}; item 133 measured {NEGATIVE_CONTROL_RECORD:.6f}; "
          f"{A5_P17_SHARES[0]:.6f} at C_mode)")
    gap = abs(_shares(two["drives"]["ccw"])[0] - _shares(two["drives"]["cw"])[0])
    say(f"(predicted, printed) |ccw - cw| feed-ring share at C_mode {gap:.3e} "
        f"(SPEC ≈ {PREDICTED_CCW_CW_GAP:.0e})")

    # ---- SPEC v2 (§9 item 161): rows 2b / 2c / 4a through ring_mode_rows_v2 ----
    v2 = rows_v2(s32, port_map, order, freq_hz)
    fixed, own = v2["fixed_c"], v2["own_c_mode"]
    c_own = own["c_f"]
    b0_rel = abs(c_own / C_MODE_F - 1.0)
    check("(B0) own C_mode", c_own == fixed["c_mode_f"] and b0_rel <= LEG_D0_REPRODUCTION_BAND,
          f"own_c_mode.c_f {c_own:.12e} bitwise == fixed_c.c_mode_f: {c_own == fixed['c_mode_f']}; "
          f"vs C_MODE_F relative {b0_rel:.3e} (band {LEG_D0_REPRODUCTION_BAND:.0e}); "
          f"vs SPEC C {SPEC_C_F:.5e}: {SPEC_C_F / c_own - 1.0:+.4e}; anti(k = 1) {v2['anti_k1']:.6f}")
    shares_and_folds("(B1) 2b P17 alone", own["drives"]["p17_alone"], A5_P17_SHARES, A5_P17_FOLDS)
    shares_and_folds("(B2) 2b cw", own["drives"]["cw"], A4_CW_SHARES, A4_CW_FOLDS)
    b3 = 0.0
    for drv in ("p17_alone", "cw"):
        for ring in ("bottom", "top"):
            a, b = own["drives"][drv]["rings"][ring], two["drives"][drv]["rings"][ring]
            b3 = max(b3, abs(a["share_k1"] - b["share_k1"]), abs(a["leak_k0"] - b["leak_k0"]),
                     max(abs(x - y) for x, y in zip(a["fractions"], b["fractions"])))
    check("(B3) 2b vs C_MODE_F blocks", b3 <= B3_ATOL,
          f"P17-alone and cw, every share / leak / fraction: max |d| {b3:.3e} (band {B3_ATOL:.0e}; "
          f"expected ≈ 4e-13)")
    sh_ccw = own["drives"]["ccw"]["rings"]["bottom"]["share_k1"]
    sh_cw = own["drives"]["cw"]["rings"]["bottom"]["share_k1"]
    b4 = abs(sh_ccw - sh_cw)
    check("(B4) 2b |ccw - cw| feed ring", abs(b4 - B4_CCW_CW_GAP) <= B4_ATOL and sh_ccw < sh_cw,
          f"{b4:.6e} vs {B4_CCW_CW_GAP:.3e} ± {B4_ATOL:.0e}; ccw {sh_ccw:.6f} < cw {sh_cw:.6f}: {sh_ccw < sh_cw}")
    in5 = B5_CCW_RANGE[0] <= sh_ccw <= B5_CCW_RANGE[1]
    say(f"(B5, predicted, printed) 2b ccw feed-ring share {sh_ccw:.6f} in {B5_CCW_RANGE}: {in5}")
    for drv in SPEC_DRIVES:
        d_feed = (fixed["drives"][drv]["rings"]["bottom"]["share_k1"]
                  - own["drives"][drv]["rings"]["bottom"]["share_k1"])
        d_other = (fixed["drives"][drv]["rings"]["top"]["share_k1"]
                   - own["drives"][drv]["rings"]["top"]["share_k1"])
        say(f"(B5, predicted, printed) row 2 - 2b {drv}: feed ring {d_feed:+.3e} in "
            f"{B5_ROW2_MINUS_2B_FEED}: {B5_ROW2_MINUS_2B_FEED[0] <= d_feed <= B5_ROW2_MINUS_2B_FEED[1]}; "
            f"other ring {d_other:+.3e} in {B5_ROW2_MINUS_2B_OTHER}: "
            f"{B5_ROW2_MINUS_2B_OTHER[0] <= d_other <= B5_ROW2_MINUS_2B_OTHER[1]}")

    ctl_own, ctl_fix = v2["controls"]["own_c_mode"], v2["controls"]["fixed_c"]
    opp_own = ctl_own["opposite"]["drives"]["opposite"]
    four_own = ctl_own["four_ccw"]["drives"]["four_ccw"]
    shares_and_folds("(C1) 2c opposite pair, own C", opp_own, A6_OPPOSITE_SHARES, None)
    c1 = max(_leak_max(opp_own), _even_max(opp_own))
    check("(C1) opposite pair leak / even-|k|, own C", c1 <= CONTROL_LEAK_BAR,
          f"max {c1:.3e} (bar {CONTROL_LEAK_BAR:.0e})")
    c2 = max(_leak_max(four_own), _even_max(four_own))
    check("(C2) four-port ccw leak / even-|k|, own C", c2 <= CONTROL_LEAK_BAR,
          f"max {c2:.3e} (bar {CONTROL_LEAK_BAR:.0e})")
    sh4 = _shares(four_own)
    d4 = max(abs(sh4[0] - A6_FOUR_PORT_SHARES[0]), abs(sh4[1] - A6_FOUR_PORT_SHARES[1]))
    wrong = max(four_own["rings"][r]["fractions"][k] for r in ("bottom", "top") for k in (1, 5, 9, 13))
    say(f"(C2, predicted, printed) four-port ccw shares {sh4[0]:.6f} / {sh4[1]:.6f}, max |d| vs A6-cw "
        f"{d4:.3e} (within {C2_PREDICTED_SHARE_ATOL:.0e}: {d4 <= C2_PREDICTED_SHARE_ATOL}); "
        f"wrong-sense bins {{1, 5, 9, 13}} max {wrong:.3e}")
    c3 = max(_leak_max(ctl_fix[lbl]["drives"][lbl]) for lbl in SPEC_CONTROLS)
    check("(C3) controls leak at SPEC C", c3 <= CONTROL_LEAK_BAR,
          f"max {c3:.3e} (bar {CONTROL_LEAK_BAR:.0e})")
    d3 = max(abs(a - b) for lbl in SPEC_CONTROLS
             for a, b in zip(_shares(ctl_fix[lbl]["drives"][lbl]), _shares(ctl_own[lbl]["drives"][lbl])))
    say(f"(C3, predicted, printed) control shares SPEC C vs own C: max |d| {d3:.3e} "
        f"(within {C3_PREDICTED_SHARE_ATOL:.0e}: {d3 <= C3_PREDICTED_SHARE_ATOL})")
    for key, ctl, rows_c in (("own C", ctl_own, own), ("SPEC C", ctl_fix, fixed)):
        ctl_min = min(min(_shares(ctl[lbl]["drives"][lbl])) for lbl in SPEC_CONTROLS)
        two_max = max(max(_shares(rows_c["drives"][d])) for d in SPEC_DRIVES)
        ratio = (_shares(ctl["opposite"]["drives"]["opposite"])[0]
                 / rows_c["drives"]["ccw"]["rings"]["bottom"]["share_k1"])
        check(f"(C4) controls > two-port, {key}", ctl_min > two_max,
              f"min control share {ctl_min:.6f} > max two-port share {two_max:.6f}; (printed) "
              f"opposite / ccw feed ring {ratio:.4f} (predicted ≈ 1.253, ceiling "
              f"{1.0 / rows_c['drives']['ccw']['rings']['bottom']['share_k1']:.4f})")

    f0_closed = freq_hz * np.sqrt(C_MODE_F / SPEC_C_F)
    f1 = abs(v2["f0_one_point_hz"] / f0_closed - 1.0)
    check("(F1) 4a one-point f0", f1 <= F1_RTOL and v2["f0_one_point_hz"] < freq_hz,
          f"{v2['f0_one_point_hz']:.2f} Hz vs f·√(C_MODE_F/SPEC_C) {f0_closed:.2f} Hz: relative "
          f"{f1:.3e} (band {F1_RTOL:.0e}); f0 - f = {v2['f0_one_point_hz'] - freq_hz:+.2f} Hz")
    delta = 1.0 - v2["f0_one_point_hz"] / freq_hz
    say(f"(F2, printed) one-point first-order error (1 - p)·δ/2 with δ = {2.0 * delta:+.3e} "
        f"(SPEC C / C_own - 1), p = d ln Im λ₁ / d ln ω unmeasured")

    # Negative controls (asserted).
    v2_c2 = rows_v2(s32, port_map, order, freq_hz, c_spec_f=C2_RECORD_F)
    nc2 = v2_c2["fixed_c"]["drives"]["p17_alone"]["rings"]["bottom"]["share_k1"]
    same_own = json.dumps(v2_c2["own_c_mode"], sort_keys=True) == json.dumps(own, sort_keys=True)
    check("(NC v2) C_2 through ring_mode_rows_v2", abs(nc2 - NEGATIVE_CONTROL_RECORD) <= SHARE_ATOL and same_own,
          f"fixed_c P17-alone feed ring {nc2:.6f} vs {NEGATIVE_CONTROL_RECORD:.6f} ± {SHARE_ATOL:.0e}; "
          f"own_c_mode bitwise equal to the SPEC-C call: {same_own}; separation "
          f"{own['drives']['p17_alone']['rings']['bottom']['share_k1'] / nc2:.2f}x (ceiling 1/{nc2:.6f} = "
          f"{1.0 / nc2:.2f}x)")
    gauge = own["drives"]["ccw"]["rings"]["bottom"]["leak_k0"]
    ctl_leak = max(c1, c2, c3)
    check("(NC v2) gauge sensitivity", gauge > GAUGE_TWO_PORT_LEAK_FLOOR,
          f"2b two-port ccw leak {gauge:.4f} > {GAUGE_TWO_PORT_LEAK_FLOOR} against controls ≤ "
          f"{ctl_leak:.3e}; (printed) separation {gauge / ctl_leak:.3e} (predicted ≈ 1.8e5)")

    if fails:
        say(f"RED: {fails} — stop; no metrics written")
        return 1

    # Published rows: the SPEC's C and drives.
    feeds = [port_map[order.index(n)] for n in SPEC_FEED_NAMES]
    rows = rm.ring_mode_rows(s32, port_map=port_map, freq_hz=freq_hz, c_f=SPEC_C_F,
                             feeds=feeds, drives=SPEC_DRIVES)
    metrics = {
        "case": "ANS-8 step 1 runnable half (ans:8)",
        "spec": "SPEC.md v2 (2026-10-01)",
        "fixture": "F-human birdcage_port_domain, 16 legs, 32 ring gaps, sigma_coil 800 S/m solved inside, degree 1",
        "cells": n_cells, "kappa": kappa, "freq_hz": freq_hz,
        "s32_reciprocity": rec, "s32_sigma_max": sig,
        "anchors_c_mode_f": C_MODE_F,
        "rows": rows,
        "rows_v2": v2,
        "loss_partition": None,
        "loss_partition_note": "not computed on our side in step 1 (scope); left blank",
        "aed": None,
    }
    # SPEC v3 (§9 item 369): stage V3's rows, if present, re-checked and added.
    v3_path = outdir / V3_JSON_NAME
    if v3_path.exists():
        from tests.validation.test_birdcage_power_identity import DISCRETE_IDENTITY_RTOL
        v3 = json.loads(v3_path.read_text())
        r1, r2, ok = check_v3_power(v3, DISCRETE_IDENTITY_RTOL)
        say(f"(v3) {v3_path.name}: (P_acc-1) {r1:.3e} from the stored scalars (band "
            f"{DISCRETE_IDENTITY_RTOL:g}): {'PASS' if ok else 'FAIL'}; gap (printed) {r2:.4e}")
        if not ok:
            say("RED: (v3) power re-check — stop; no metrics written")
            return 1
        metrics["rows_v3"] = v3
        metrics["spec"] = V3_SPEC_LABEL
    else:
        say("v3 rows absent")
    write_outputs(metrics)
    say(f"wrote {CASE_DIR / 'metrics.json'} and {CASE_DIR / 'COMPARISON.md'}; "
        f"private: {_write_private_if_available(metrics)}; elapsed {time.perf_counter() - t0:.2f} s")
    return 0


def _fmt_share_table(rows, title):
    lines = [f"### {title}", "",
             "| Drive | ring | k = ±1 share | k = 0 leak | \\|k\\| folds 0…8 |", "|---|---|---|---|---|"]
    for label, blk in rows["drives"].items():
        for ring in ("bottom", "top"):
            r = blk["rings"][ring]
            lines.append(f"| {label} | {ring} | {r['share_k1']:.6f} | {r['leak_k0']:.6f} | "
                         + " ".join(f"{x:.4f}" for x in r["abs_k_folds"]) + " |")
    return lines


def comparison_text(metrics) -> str:
    rows = metrics["rows"]
    lines = ["# ANS-8 step 1 — comparison (our side; AED column blank)", "",
             "Rows from `fem_em_solver.ports.ring_modes.ring_mode_rows` on our 32×32 "
             f"(F-human, 64 MHz, degree 1, κ = {metrics['kappa']:.6e}, {metrics['cells']} cells), "
             f"C = {rows['c_f']:.4e} F on every gap, feeds P17 / P21 (Norton, 50 Ω). "
             "AED values are never written here (PROJECT_PLAN §5.4).", "",
             "| Row | ours | AED |", "|---|---|---|",
             f"| (1) C_mode from k = 1 [pF] | {rows['c_mode_f'] * 1e12:.4f} | |",
             f"| (1) C from k = 2 [pF] | {rows['c_k2_f'] * 1e12:.4f} | |",
             f"| (1) off-block norm | {rows['off_block']:.6e} | |"]
    v2 = metrics.get("rows_v2")
    if v2 is not None:
        lines.append(f"| (1) k = 1 eigenpair ring-antisymmetric fraction (v2) | {v2['anti_k1']:.6f} | |")
        lines.append(f"| (2) rows below at C = {rows['c_f'] * 1e12:.2f} pF: detuning record (SPEC v2) | | |")
    for label in ("p17_alone", "ccw", "cw"):
        blk = rows["drives"][label]
        lines.append(f"| (2) {label} feed-ring (bottom) k = ±1 share | "
                     f"{blk['rings']['bottom']['share_k1']:.6f} | |")
        lines.append(f"| (2) {label} feed-ring (bottom) k = 0 leak | "
                     f"{blk['rings']['bottom']['leak_k0']:.6f} | |")
        lines.append(f"| (2) {label} other ring (top) k = ±1 share | "
                     f"{blk['rings']['top']['share_k1']:.6f} | |")
    for n, f in zip(SPEC_FEED_NAMES, rows["drives"]["ccw"]["feeds"]):
        lines.append(f"| (3) ccw V/I at {n} [Ω] | {f['z_port_ohm'][0]:+.6e} {f['z_port_ohm'][1]:+.6e}j | |")
    y = rows["y_feed_2x2_s"]
    lines.append("| (3) feed 2×2 Y [S], others in C | " + "; ".join(
        ", ".join(f"{e[0]:+.6e}{e[1]:+.6e}j" for e in row) for row in y) + " | |")
    if v2 is not None:
        lines += _v2_table_rows(v2)
    lines += ["", "Loss partition: not computed on our side in step 1 (left blank).", ""]
    lines += _fmt_share_table(rows, "Per-ring k-spectra (ours)")
    if v2 is not None:
        lines += [""] + _fmt_share_table(
            v2["own_c_mode"], "Per-ring k-spectra at own C_mode (2b, ours, v2)")
    lines += ["", "Not claimed: any AED value, verdict or agreement level; sweep row 4; Route A; "
              "copper. The weekly adjudicates share and leak at a stated relative level."]
    if v2 is not None:
        lines += ["", "Not claimed (v2): the row-4 sweep interpolation (4a is the one-point form, "
                  "λ₁ taken ∝ ω, its first-order error (1 − p)·δ/2 unmeasured); the degree-2 32×32; "
                  "copper; P_coil; Route A; adjudication."]
    v3 = metrics.get("rows_v3")
    if v3 is not None:
        lines += [""] + _v3_section(v3)
    d2 = metrics.get("rows_v3_d2")
    if d2 is not None and v3 is not None:  # the block follows the v3 section, never alone
        lines += [""] + _v3_d2_section(d2, v3)
    return "\n".join(lines) + "\n"


def _v3_d2_section(d2, v3=None):
    """SPEC v4's degree-2 column at the root (§10 H31 Part B, §9 item 597), from H30's
    node-record rows.  Pure; insertions only.  Every figure is in the H30 log's own
    format spec (P_acc ``.12e``, row 5 ``.6f``, X_s ``.9e``, ΣX_s / P_acc ``.6e``,
    row 4 ``.7f``); the degree-1 value from ``rows_v3`` sits beside it, labelled."""
    p = d2["power_w"]
    f = d2["row5"]
    mech = d2["gap"]["mechanism"]

    def _d1(text):
        return "" if text is None else f" (degree 1 record: {text})"

    # The `.get` guards keep a synthetic ``rows_v3`` without ``gap`` (test_ans8_v3_rows'
    # one) renderable; its degree-1 X_s cells are then left without a beside value.
    d1 = {}
    if v3 is not None:
        f1 = v3["row5"]
        p1 = v3["power_w"]["p_acc_spec"]
        d1 = {"p_acc": f"{p1:.12e}", "row4": f"{v3['row4']['sar10g_per_w']:.7f}",
              "phantom": f"{f1['phantom']:.6f}", "copper": f"{f1['copper']:.6f}",
              "esr": f"{f1['esr']:.6f}", "sum": f"{f1['sum']:.6f}"}
        m1 = v3.get("gap", {}).get("mechanism")
        if m1 is not None:
            d1["sx"] = f"{sum(m1[n]['x_w'] for n in v3['feeds']) / p1:.6e}"
            for n in v3["feeds"]:
                d1[f"x_{n}"] = f"{m1[n]['x_w']:.9e}"
    lines = ["## Step 2 at the degree-2 root (SPEC v4: rows 4–5 from the node record, degree 2; "
             "ours; AED column blank)", "",
             f"This block is the degree-2 column at the root, labelled per SPEC v4: degree "
             f"{d2['degree']}, C_in = {d2['c_in_pf']} pF on all 32 gaps, ccw (1, −j) at "
             f"{', '.join(d2['feeds'])}, 64 MHz, rows from the solved node record (knob "
             f"{d2['node_knob']}). P_acc,SPEC is the v3 section's line cut. "
             "AED values are never written here (PROJECT_PLAN §5.4).", "",
             "| Row | ours | AED |", "|---|---|---|",
             f"| P_acc,SPEC [W] | {p['p_acc_spec']:.12e}{_d1(d1.get('p_acc'))} | |",
             f"| (4) peak 10 g SAR per W of P_acc (sphere, lattice argmax) [W/kg/W] | "
             f"{d2['row4']['sar10g_per_w']:.7f}{_d1(d1.get('row4'))} | |"]
    for key, lab in (("phantom", "phantom"), ("copper", "copper (Leontovich)"),
                     ("esr", "ESR (30 gaps + 2 feed gaps)"), ("sum", "sum of the")):
        lines.append(f"| (5) {lab} fraction{'s' if key == 'sum' else ''} of P_acc | "
                     f"{f[key]:.6f}{_d1(d1.get(key))} | |")
    for n in d2["feeds"]:
        lines.append(f"| X_s at feed sheet {n} [W] | {mech[n]['x_w']:.9e}{_d1(d1.get('x_' + n))} | |")
    lines.append(f"| ΣX_s / P_acc | {d2['sum_x_over_p_acc']:.6e}{_d1(d1.get('sx'))} | |")
    row2 = d2.get("row2_circuit_route_h25")  # H35 (§9 item 623), rendered since §9 item 740
    if row2 is not None:
        log, cited = row2["source"].rsplit(":", 1)  # the page's line ranges carry en dashes
        lossless, esr = row2["lossless"], row2["esr_0p2_ohm"]
        lines.append(f"| row 2 (degree 2, **circuit route — Norton feeds on entry 36's saved "
                     f"32×32 at C_mode(d2, Cu) = {row2['c_mode_pf']} pF; in-model row not in "
                     f"record**) | ccw feed-ring share {lossless['ccw_feed_share']} / k = 0 leak "
                     f"{lossless['ccw_feed_leak_k0']} (lossless); {esr['ccw_feed_share']} / "
                     f"{esr['ccw_feed_leak_k0']} (r = 0.2 Ω, the production convention); source "
                     f"`{log}:{cited.replace('-', '–')}` | |")
    lines += ["", "Rows 1–3, row 6 and CS_s are not in the node record at degree 2 (§10 H30 as "
              "narrowed by the 2026-10-21 weekly).", "",
              "Not claimed (degree 2): any AED value, verdict or agreement level; rows 1–3 / 6 / "
              "CS_s; an order-convergence statement; f₀, Q, match or tuning accuracy; compliance."]
    if row2 is not None:
        lines += ["", f"Row 2 is now carried at degree 2 from the circuit route (§10 H35), "
                  f"labelled as that route at H25's own C_mode(d2, Cu) = {row2['c_mode_pf']} pF, "
                  f"not at the {d2['c_in_pf']} pF root; rows 1, 3, 6 and CS_s stay not in record."]
    return lines


def _v3_section(v3):
    """SPEC v3's table (ours; AED cells empty).  Insertions only: every step-1
    line above stays as it was."""
    r4 = v3["row4"]
    lines = ["## Step 2 (SPEC v3) — tuned two-port drive with ESR, degree 1 (ours; AED column blank)", "",
             "The tables above are **lossless (step 1)** (σ = 800 S/m coil, C = 113.11 pF, no ESR). "
             f"This section: H9′'s copper F-human (hole mesh, Leontovich, κ off), degree 1, the 0.015 m "
             f"phantom, C_in = {v3['c_in_pf']} pF and r = {v3['r_ohm']} Ω on all 32 gaps, ccw (1, −j) "
             f"at {', '.join(v3['feeds'])}, 64 MHz, {v3['cells']} cells. "
             "AED values are never written here (PROJECT_PLAN §5.4).", "",
             _v3_p_acc_definition(v3), "",
             "| Row | ours | AED |", "|---|---|---|"]
    for n, (re_, im_) in v3["row1_y_in_ccw_s"].items():
        lines.append(f"| (1) Y_in ccw at {n} [S] | {re_:+.6e} {im_:+.6e}j | |")
    for ring, lab in (("bottom", "feed ring (bottom)"), ("top", "other ring (top)")):
        r = v3["row2_field"][ring]
        lines.append(f"| (2) {lab} k = ±1 share (in-model sheet currents) | {r['share_k1']:.6f} | |")
        lines.append(f"| (2) {lab} k = 0 leak | {r['leak_k0']:.6e} | |")
    pp_ = v3.get("per_partition")

    def _beside(val):
        return "" if pp_ is None else f" (per W of the field partition: {val:.6e})"

    for n, r in v3["row3"].items():
        lines.append(f"| (3) \\|B₁⁺\\| at {n} per √W [T/√W] | {r['b1p_per_sqrt_w']:.6e}"
                     + (_beside(pp_["row3_b1p_per_sqrt_w"][n]) if pp_ else "") + " | |")
    c = r4["c_star_m"]
    lines.append(f"| (4) peak 10 g SAR per W (sphere, lattice argmax) [W/kg/W] | {r4['sar10g_per_w']:.6e}"
                 + (_beside(pp_["row4_sar10g_per_w"]) if pp_ else "") + " | |")
    lines.append(f"| (4) its location (x, y, z) [m] | ({c[0]:.6f}, {c[1]:.6f}, {c[2]:.6f}) | |")
    f = v3["row5"]
    fp = v3.get("row5_partition")

    def _frac(key):
        return "" if fp is None else f" (of the field partition: {fp[key]:.6f})"

    lines.append(f"| (5) phantom fraction of P_acc | {f['phantom']:.6f}{_frac('phantom')} | |")
    lines.append(f"| (5) copper (Leontovich) fraction of P_acc | {f['copper']:.6f}{_frac('copper')} | |")
    lines.append(f"| (5) ESR fraction of P_acc (30 gaps + 2 feed gaps) | {f['esr']:.6f}{_frac('esr')} | |")
    if "gap" in v3:
        lines.append(f"| (5) sum of the fractions of P_acc; residual = feed-sheet Cauchy–Schwarz excess "
                     f"(documented sheet property, weekly 2026-10-12; known-issues 2026-10-11) | "
                     f"{f['sum']:.6f}; residual "
                     f"{1.0 - f['sum']:.6f} | |")
    lines.append(f"| (6) \\|B₁⁺\\| / \\|B₁⁻\\| at O | {v3['row6_purity_o']:.4f} | |")
    lines += ["", "Row 2 is `ring_modes.k_spectrum` on the field's sheet currents; the circuit "
              "reduction of the stored 32×32 is not our row (H12c identity red at 1.58). c* is a "
              "lattice argmax (18 points), not a global maximum.", "",
              "Not claimed (v3): any AED value, verdict or agreement level; the degree-2 column; f₀, "
              "Q, match or tuning accuracy; c* as a global maximum; homogeneity; compliance."]
    return lines


def _v3_p_acc_definition(v3):
    """SPEC v3's accepted-power definition (`SPEC.md` "Accepted power, defined";
    weekly ruling 2026-10-12, §9 item 411).  Every figure is rendered from the
    stored floats; the `.get` guards keep a ``rows_v3`` without ``gap`` /
    ``p_acc_wave`` (the unit test's synthetic one) renderable."""
    p = v3["power_w"]
    txt = ("**Accepted power (SPEC v3, weekly ruling 2026-10-12).** Per-W / per-√W figures are "
           "per W of P_acc, the sum over the two feed ports of the accepted power at the 50 Ω "
           "lumped port: incident − reflected at the port plane, HFSS's default per-port "
           "accepted power. Ours is the line cut ½ Re Σ_feeds V_port·conj(I_line) = "
           f"{p['p_acc_spec']:.6e} W at unit Thevenin drive")
    wave = p.get("p_acc_wave")
    if wave is not None:
        txt += (f" (relative difference from the wave form (‖v_s‖² − ‖S v_s‖²) / (8 z0): "
                f"{abs(p['p_acc_spec'] - wave) / abs(wave):.1e})")
    txt += (f". The feed gaps' capacitor and {v3['r_ohm']} Ω ESR are load inside P_acc: their "
            "dissipation is part of row 5's ESR fraction.")
    gap = v3.get("gap")
    if gap is not None:
        mech = gap.get("mechanism", {})
        cs = ", ".join(f"{mech[n]['cs_over_x']:.9f} at {n}" for n in v3["feeds"] if n in mech)
        txt += (f" The field's feed-sheet partition exceeds the line cut by "
                f"{100.0 * gap['rel_line_cut']:.2f} %: the discrete feed sheet's Cauchy–Schwarz "
                f"excess (CS / X = {cs}), a documented sheet property (weekly ruling 2026-10-12; "
                "known-issues 2026-10-11), not a defect and not a quantity the AED column is "
                "asked for.")
    txt += " The per-partition values are printed beside each row."
    return txt


def _v2_table_rows(v2):
    own = v2["own_c_mode"]
    lines = [f"| (2b) own C_mode [pF] (from row 1, every gap incl. feeds) | {v2['c_own_f'] * 1e12:.6f} | |"]
    for label in ("p17_alone", "ccw", "cw"):
        r = own["drives"][label]["rings"]
        lines.append(f"| (2b) {label} feed-ring (bottom) k = ±1 share | {r['bottom']['share_k1']:.6f} | |")
        lines.append(f"| (2b) {label} feed-ring (bottom) k = 0 leak | {r['bottom']['leak_k0']:.6f} | |")
        lines.append(f"| (2b) {label} other ring (top) k = ±1 share | {r['top']['share_k1']:.6f} | |")
    for n, f in zip(SPEC_FEED_NAMES, own["drives"]["ccw"]["feeds"]):
        lines.append(f"| (3b) ccw V/I at {n}, own C [Ω] | {f['z_port_ohm'][0]:+.6e} {f['z_port_ohm'][1]:+.6e}j | |")
    lines.append("| (3b) feed 2×2 Y [S], others in own C | " + "; ".join(
        ", ".join(f"{e[0]:+.6e}{e[1]:+.6e}j" for e in row) for row in own["y_feed_2x2_s"]) + " | |")
    for key, cname in (("fixed_c", f"C = {v2['c_spec_f'] * 1e12:.2f} pF"), ("own_c_mode", "own C_mode")):
        for label in SPEC_CONTROLS:
            r = v2["controls"][key][label]["drives"][label]["rings"]
            lines.append(f"| (2c) {label}, {cname}: feed ring share / leak; other ring share / leak | "
                         f"{r['bottom']['share_k1']:.6f} / {r['bottom']['leak_k0']:.3e}; "
                         f"{r['top']['share_k1']:.6f} / {r['top']['leak_k0']:.3e} | |")
    lines.append(f"| (4a) f₀ = f·√(C_mode / C) [Hz] (one-point, λ₁ ∝ ω) | {v2['f0_one_point_hz']:.2f} | |")
    return lines


def write_outputs(metrics, case_dir=None):
    case_dir = CASE_DIR if case_dir is None else Path(case_dir)
    (case_dir / "metrics.json").write_text(json.dumps(metrics, indent=1, sort_keys=True) + "\n")
    (case_dir / "COMPARISON.md").write_text(comparison_text(metrics))


# -- private writer ---------------------------------------------------------------
def _parse_complex_matrix(obj):
    """S from the AED JSON: ``{"re": [[…]], "im": [[…]]}``, rows of ``{"re", "im"}``
    or ``[re, im]``, or plain numbers."""
    if isinstance(obj, dict) and "re" in obj and "im" in obj:
        return np.array(obj["re"], dtype=float) + 1j * np.array(obj["im"], dtype=float)

    def one(e):
        if isinstance(e, dict):
            return complex(e["re"], e["im"])
        if isinstance(e, (list, tuple)) and len(e) == 2:
            return complex(e[0], e[1])
        return complex(e)

    return np.array([[one(e) for e in row] for row in obj], dtype=np.complex128)


def _aed_designs(aed):
    """Every ``(name, block)`` in the AED JSON carrying an ``"S"`` 32×32."""
    out = []

    def walk(prefix, node):
        if not isinstance(node, dict):
            return
        if "S" in node:
            out.append((prefix or "aed", node))
        for k, v in node.items():
            if not str(k).startswith("_") and k != "S":
                walk(f"{prefix}/{k}" if prefix else str(k), v)

    walk("", aed)
    return out


def _load_aed_results(path=None):
    path = AED_RESULTS_PATH if path is None else Path(path)
    if not path.exists():
        return None
    return json.loads(path.read_text())


def aed_rows(block, freq_hz):
    """AED's rows through the same ``ring_mode_rows_v2`` as ours — 2b at AED's
    **own** row-1 C_mode, wired inside the shared function (never ours)."""
    rm = _ring_modes()
    s = _parse_complex_matrix(block["S"])
    names = block.get("port_names") or block.get("ports") or _spec_port_names()
    port_map = rm.ring_port_map_from_names(names)
    return rows_v2(s, port_map, names, float(block.get("freq_hz", freq_hz)))


PRIVATE_HEADER = "# ANS-8 step 1 — PRIVATE comparison (AED values; gitignored, never commit)"


def _private_keys():
    """``(label, accessor on a rows_v2 dict, absolute gap?)`` per private column."""
    keys = [("C_mode [pF]", lambda r: r["fixed_c"]["c_mode_f"] * 1e12, False),
            ("C_k2 [pF]", lambda r: r["fixed_c"]["c_k2_f"] * 1e12, False),
            ("off-block", lambda r: r["fixed_c"]["off_block"], False)]
    for blk, tag in (("fixed_c", ""), ("own_c_mode", "2b ")):
        for label in SPEC_DRIVES:
            for ring, q, qn in (("bottom", "share_k1", "share bottom"), ("bottom", "leak_k0", "leak bottom"),
                                ("top", "share_k1", "share top")):
                keys.append((f"{tag}{label} {qn}",
                             lambda r, b=blk, l=label, g=ring, qq=q: r[b]["drives"][l]["rings"][g][qq], False))
    for blk, tag in (("fixed_c", "C"), ("own_c_mode", "own C")):
        for label in SPEC_CONTROLS:
            for ring in ("bottom", "top"):
                keys.append((f"2c {label} @{tag} share {ring}",
                             lambda r, b=blk, l=label, g=ring:
                             r["controls"][b][l]["drives"][l]["rings"][g]["share_k1"], False))
                # Control leaks: absolute gaps (≈ 0 by symmetry; a relative gap divides by zero).
                keys.append((f"2c {label} @{tag} leak {ring}",
                             lambda r, b=blk, l=label, g=ring:
                             r["controls"][b][l]["drives"][l]["rings"][g]["leak_k0"], True))
    keys.append(("4a f0 [MHz]", lambda r: r["f0_one_point_hz"] / 1e6, False))
    return keys


def _private_comparison_text(metrics, aed) -> str:
    ours = metrics["rows_v2"]
    keys = _private_keys()
    lines = [PRIVATE_HEADER, "",
             "| Source | " + " | ".join(k for k, _, _ in keys) + " |",
             "|---|" + "---|" * len(keys)]
    lines.append("| ours | " + " | ".join(f"{fn(ours):.6f}" for _, fn, _ in keys) + " |")
    gaps = []
    for name, block in _aed_designs(aed):
        r = aed_rows(block, metrics["freq_hz"])
        lines.append(f"| AED {name} | " + " | ".join(f"{fn(r):.6f}" for _, fn, _ in keys) + " |")
        gaps.append(f"| {name} | " + " | ".join(
            (f"{fn(ours) - fn(r):+.3e} (abs)" if absolute else f"{(fn(ours) - fn(r)) / fn(r):+.3e}")
            for _, fn, absolute in keys) + " |")
    lines += ["", "## Gaps (ours − AED) / AED; control leaks as absolute ours − AED", "",
              "| AED | " + " | ".join(k for k, _, _ in keys) + " |", "|---|" + "---|" * len(keys)]
    lines += gaps
    return "\n".join(lines) + "\n"


def _write_private_if_available(metrics) -> str:
    """`ANS-6`'s pattern: knob first, then the JSON's presence; then the clobber
    guard — an existing private file whose first line is not this writer's
    header is kept byte-identical (only line 1 is read)."""
    if os.environ.get(SKIP_PRIVATE_ENV) == "1":
        return f"skipped ({SKIP_PRIVATE_ENV}=1) — private table not written"
    aed = _load_aed_results()
    if aed is None:
        return "no aed_results/ on this box — private table not written"
    path = CASE_DIR / "COMPARISON_private.md"
    if path.exists():
        with path.open(encoding="utf-8", errors="replace") as fh:
            first = fh.readline().rstrip("\n")
        if first != PRIVATE_HEADER:
            return f"{path.name} kept — its first line is not this writer's header (never overwritten)"
    path.write_text(_private_comparison_text(metrics, aed))
    return f"{path.name} — AED columns filled (untracked, never commit)"


def main() -> int:
    stage = os.environ.get(STAGE_ENV, "").strip().upper()
    outdir = _results_dir()
    if stage == "S":
        from mpi4py import MPI
        _stage_s(MPI.COMM_WORLD, outdir)
        return 0
    if stage == "R":
        return _stage_r(outdir)
    if stage == "V3":
        from mpi4py import MPI
        return _stage_v3(MPI.COMM_WORLD, outdir)
    print(f"[ANS-8] set {STAGE_ENV}=S (heavy, -n 8) then {STAGE_ENV}=R (-n 1); SPEC v3: "
          f"{STAGE_ENV}=V3 (heavy, -n 8) before R", flush=True)
    return 2


if __name__ == "__main__":
    sys.exit(main())
