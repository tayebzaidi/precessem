"""`WF-10` step 1 (§9 item 52, chain H8): tune F-human in the circuit layer.

From `WF-7` step 1's F-human rung (branch B, longitudinal ring sheets, 64 MHz,
degree 1) the 32×32 is recomputed in-window with κ-corrected sheet widths
(κ = ``terminal_form_deficit``'s pooled ``C/terminal − 1`` on the P17 column,
as ``test_wf7_f_human_ring_matrix.py`` computes it).  Two feed ports — P17
and the port 90° from it on the same ring — are kept at 50 Ω; the other 30
ring gaps terminate in ``C``, and ``PORT-15`` step 3's ``tuning_sweep`` /
``select_c_tuned`` find ``C_tuned``.

**Asserted (imported, never re-banded):** the cell count within
``CELL_COUNT_BAND``; ``|Im Z_in|/|Z_in| ≤ TUNING_IM_Z_RTOL`` at ``C_tuned``;
the circuit-reduced 2×2 equals one in-model solve with `PORT-14` capacitor
sheets at ``C_tuned`` on the 30 gaps within ``REDUCTION_BAND``; reciprocity
≤ ``RECIPROCITY_BAND`` on the 32×32 and both 2×2s; ``MATRIX_PASSIVITY_CEILING``
on the 32×32 and ``COLUMN_PASSIVITY_CEILING`` on the 2×2s.

**Negative control (predicted, printed — §9 rule (e)):** the 0.5·C_tuned
circuit prediction against the in-model 2×2 at ``C_tuned``, beside the
predicted ≥ 10× the band.  No prior F-human measurement of that comparison
exists, so it is not asserted here.

**Printed:** κ; every zero/pole of the sweep with the k = ±1 DFT share of the
16 slot currents per ring (circuit-side currents, feed 1 driven at 1 V through
50 Ω, feed 2 in 50 Ω); ``R_in``; feed isolation |S21| (circuit and in-model);
``cond(I − Γ S_tt)``; a capacitor-sensitivity Q proxy.

**Traps (item text):** ``_solve_one_drive`` hard-codes a 50 Ω ``z_p`` for V
and for its power residual, so on the C sheets neither is printed; the kept
(feed) ports are 50 Ω sheets, so the kept-port S is correct.

Heavy, ``-n 8``, gated by ``FEM_EM_WF10_STEP1=1``.
"""

from __future__ import annotations

import os
import resource
import time

import numpy as np
import pytest
from mpi4py import MPI

from tests.complex_mode import complex_only

STEP1_ENV = "FEM_EM_WF10_STEP1"
NEGATIVE_CONTROL_FACTOR = 0.5
NEGATIVE_CONTROL_PREDICTED_BAND_MULTIPLE = 10.0
# κ beside `WF-7` step 1's pooled range (printed only).
WF7_STEP1_POOLED_RANGE = "4.591e-3 .. 4.689e-3"
QUADRATURE_SLOT_OFFSET = 4  # 90 deg on the 22.5 deg ring grid


def _on() -> bool:
    return os.environ.get(STEP1_ENV, "").strip() == "1"


def _rss_summed_gib(comm) -> float:
    local = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss) * 1024.0 / 2**30
    return float(comm.allreduce(local, op=MPI.SUM))


def _respec(specs, kappa, impedances=None):
    """``specs`` rebuilt with ``width_correction_kappa=kappa`` (and per-port Z)."""
    from fem_em_solver.ports.lumped import LumpedSheetPortSpec

    impedances = impedances or {}
    return [
        LumpedSheetPortSpec(
            port_id=spec.port_id,
            facet_tag=spec.facet_tag,
            port_impedance_ohm=impedances.get(spec.port_id, spec.port_impedance_ohm),
            gap_height_m=spec.gap_height_m,
            sheet_width_m=spec.sheet_width_m,
            drive_direction=spec.drive_direction,
            drive_voltage_v=spec.drive_voltage_v,
            interior=spec.interior,
            width_correction_kappa=kappa,
        )
        for spec in specs
    ]


def _circuit_currents(s, frequency_hz, c_f, n_feeds, z0):
    """Port currents: feed 0 driven 1 V through z0, other feeds z0, rest in C."""
    from fem_em_solver.ports.circuit import s_to_z
    from tests.validation.test_port_circuit_layer_field import _capacitor_impedance

    z = s_to_z(np.asarray(s, dtype=np.complex128), z0)
    n = z.shape[0]
    z_load = np.full(n, _capacitor_impedance(frequency_hz, c_f), dtype=np.complex128)
    z_load[:n_feeds] = z0
    rhs = np.zeros(n, dtype=np.complex128)
    rhs[0] = 1.0
    return np.linalg.solve(z + np.diag(z_load), rhs)


def _k1_shares(currents, order, slots):
    """``{ring sign: k=±1 share}`` of the 16 slot currents per ring (Parseval)."""
    out = {}
    for sign in (1.0, -1.0):
        vec = np.zeros(16, dtype=np.complex128)
        for idx, pid in enumerate(order):
            k, sg = slots[pid]
            if sg == sign:
                vec[k] = currents[idx]
        spec = np.fft.fft(vec)
        total = float(np.sum(np.abs(spec) ** 2))
        out[sign] = (
            float((abs(spec[1]) ** 2 + abs(spec[15]) ** 2) / total) if total else float("nan")
        )
    return out


@complex_only
def test_wf10_step1_f_human_tuned_in_the_circuit_layer_matches_in_model():
    if not _on():
        pytest.skip(f"{STEP1_ENV} unset: heavy -n 8 F-human window, environment-selected")

    from fem_em_solver.ports.shares import terminal_form_deficit
    from tests.mesh.test_birdcage_port_sheet_prerequisite import CELL_COUNT_BAND
    from tests.validation.test_birdcage_f_human_rung import F_HUMAN_BRANCH_B_CELL_RECORD
    from tests.validation.test_port_birdcage_ring_column import (
        COLUMN_PASSIVITY_CEILING,
        _reciprocity_ratio,
        _solve_one_drive,
    )
    from tests.validation.test_port_birdcage_ring_matrix import MATRIX_PASSIVITY_CEILING
    from tests.validation.test_port_circuit_layer_field import (
        TUNING_IM_Z_RTOL,
        _capacitor_impedance,
        _residual,
        select_c_tuned,
        tuned_input,
        tuning_sweep,
    )
    from tests.validation.test_port_drive_superposition import _ring_quadrature_slots
    from tests.validation.test_port_lumped_rlc_termination import REDUCTION_BAND
    from tests.validation.test_port_lumped_sheet_sweep import RECIPROCITY_BAND
    from tests.validation.test_port_package_sparameters import REFERENCE_IMPEDANCE_OHM
    from tests.validation.test_wf7_f_human_b1_quadrature import (
        DEGREE,
        FREQUENCY_HZ,
        _build_f_human,
    )

    comm = MPI.COMM_WORLD
    z0 = float(REFERENCE_IMPEDANCE_OHM)
    t_window = time.perf_counter()

    def say(msg: str) -> None:
        if comm.rank == 0:
            print(f"[WF-10 step1] {msg}", flush=True)

    say(f"start: -n {comm.size}, f = {FREQUENCY_HZ:.3e} Hz, degree {DEGREE}")

    # ---- fixture -----------------------------------------------------------
    t0 = time.perf_counter()
    built = _build_f_human(comm, emit_leg_ports=True)  # GEO-35 step 2b: unpin with the re-record
    build_s = time.perf_counter() - t0
    n_cells = built["cells"]
    rel_cells = n_cells / F_HUMAN_BRANCH_B_CELL_RECORD - 1.0
    say(
        f"(i) cells {n_cells} vs record {F_HUMAN_BRANCH_B_CELL_RECORD}: relative "
        f"{rel_cells:.3e} (band {CELL_COUNT_BAND}, ASSERTED); build {build_s:.2f} s"
    )
    ctx = built["ctx"]
    base_specs = list(ctx["specs"])

    az_ref, slots = _ring_quadrature_slots(built["sheets"])
    ref_pid = "P17"
    ref_k, ref_sign = slots[ref_pid]
    feed2 = next(
        pid
        for pid, (k, sg) in slots.items()
        if sg == ref_sign and k == (ref_k + QUADRATURE_SLOT_OFFSET) % 16
    )
    others = sorted((pid for pid in slots if pid not in (ref_pid, feed2)), key=lambda p: int(p[1:]))
    order = [ref_pid, feed2] + others
    assert len(order) == 32
    say(
        f"feeds: {ref_pid} (slot {ref_k}, ring {ref_sign:+.0f}; az_ref {az_ref:.6f} deg) and "
        f"{feed2} (slot {slots[feed2][0]}); 30 ports in C"
    )

    # ---- κ from the P17 column at the uncorrected width ------------------------
    col = _solve_one_drive(ctx, ref_pid, reuse_factorization=False)
    fields = col.pop("fields")
    t_id = time.perf_counter()
    sheets = [spec.sheet(driven=(spec.port_id == ref_pid)) for spec in base_specs]
    kappa = float(
        terminal_form_deficit(ctx["msh"], ctx["tags_f"], sheets, fields.e_complex, comm)["pooled"]
    )
    del fields
    say(
        f"kappa (pooled C/terminal - 1, P17 column, uncorrected width) = {kappa:.6e} beside "
        f"WF-7 step 1's pooled {WF7_STEP1_POOLED_RANGE} (printed); identity eval "
        f"{time.perf_counter() - t_id:.2f} s; t={time.perf_counter() - t_window:.1f} s"
    )

    # ---- κ-corrected 32×32 -------------------------------------------------
    ctx32 = dict(ctx)
    ctx32["specs"] = _respec(base_specs, kappa)
    columns = {}
    for n_drive, pid in enumerate(order):
        c = _solve_one_drive(ctx32, pid, reuse_factorization=(n_drive > 0))
        c.pop("fields")
        columns[pid] = c
    kinds = [columns[p]["solve_kind"] for p in order]
    s32 = np.array(
        [[columns[pj]["s_column"][pi] for pj in order] for pi in order], dtype=np.complex128
    )
    ratio32 = _reciprocity_ratio(s32)
    sigma32 = float(np.linalg.norm(s32, 2))
    say(
        f"32x32 (kappa-corrected): {kinds.count('solve')} factorising, {kinds.count('held')} "
        f"held; reciprocity {ratio32:.6e} (band {RECIPROCITY_BAND:.0e}, ASSERTED); sigma_max "
        f"{sigma32:.9f} (ceiling {MATRIX_PASSIVITY_CEILING}, ASSERTED); S_driven(P17) "
        f"{s32[0, 0]:.6f} (printed, width-tagged); t={time.perf_counter() - t_window:.1f} s"
    )

    # ---- the circuit sweep ---------------------------------------------------
    terminated = tuple(range(2, 32))
    _grid, _vals, roots = tuning_sweep(s32, FREQUENCY_HZ, terminated=terminated)
    tuned = select_c_tuned(roots)
    say(f"tuning sweep: {len(roots)} sign change(s) of Im Z_in")
    for r in roots:
        cur = _circuit_currents(s32, FREQUENCY_HZ, r["c_f"], 2, z0)
        sh = _k1_shares(cur, order, slots)
        say(
            f"    C = {r['c_f']:.12e} F  Z_in = {r['z'].real:+.6e} {r['z'].imag:+.6e}j Ohm  "
            f"|S11| = {abs(r['s11']):.6f}  |Im Z|/|Z| = {r['im_rel']:.3e}  "
            f"{'ZERO' if r['zero'] else 'pole'}  k=+-1 share: ring+ {sh[1.0]:.6f} "
            f"ring- {sh[-1.0]:.6f}"
        )
    assert tuned is not None, "no zero of Im Z_in on the grid — the finding; never widen in-slot"
    c_t = float(tuned["c_f"])
    _, pred2 = tuned_input(s32, FREQUENCY_HZ, c_t, terminated=terminated)
    _, pred_ctl = tuned_input(
        s32, FREQUENCY_HZ, NEGATIVE_CONTROL_FACTOR * c_t, terminated=terminated
    )
    z_c = _capacitor_impedance(FREQUENCY_HZ, c_t)
    gamma = (z_c - z0) / (z_c + z0)
    s_tt = s32[2:, 2:]
    cond = float(np.linalg.cond(np.eye(30) - gamma * s_tt))
    r_in = float(tuned["z"].real)
    # Capacitor-sensitivity Q proxy: C·dX_in/dC / R (printed only, no closed form).
    dc = 1.0e-6 * c_t
    x_hi = tuned_input(s32, FREQUENCY_HZ, c_t + dc, terminated=terminated)[0].imag
    x_lo = tuned_input(s32, FREQUENCY_HZ, c_t - dc, terminated=terminated)[0].imag
    c_dxdc = c_t * (x_hi - x_lo) / (2.0 * dc)
    say(
        f"C_tuned = {c_t:.12e} F (Z_C = {z_c:.6e} Ohm): |Im Z_in|/|Z_in| = "
        f"{tuned['im_rel']:.3e} (ASSERTED <= {TUNING_IM_Z_RTOL:.0e}); R_in = {r_in:.6e} Ohm; "
        f"circuit |S21| = {abs(pred2[1, 0]):.6e}; cond(I - Gamma S_tt) = {cond:.6e}; "
        f"C dX/dC = {c_dxdc:+.6e} Ohm -> Q proxy unloaded {abs(c_dxdc) / r_in:.3f}, loaded "
        f"{abs(c_dxdc) / (r_in + z0):.3f} (printed)"
    )

    # ---- in-model 2×2 with PORT-14 C sheets on the 30 gaps -------------------
    ctx_c = dict(ctx)
    ctx_c["specs"] = _respec(base_specs, kappa, {pid: z_c for pid in others})
    meas = {}
    t_im = time.perf_counter()
    for n_drive, pid in enumerate(order[:2]):
        c = _solve_one_drive(ctx_c, pid, reuse_factorization=(n_drive > 0))
        c.pop("fields")
        meas[pid] = c
    s2 = np.array(
        [[meas[pj]["s_column"][pi] for pj in order[:2]] for pi in order[:2]], dtype=np.complex128
    )
    r_id = _residual(s2, pred2)
    r_ctl = _residual(s2, pred_ctl)
    rec_pred, rec_meas = _reciprocity_ratio(pred2), _reciprocity_ratio(s2)
    colmax_pred = float(np.max(np.sum(np.abs(pred2) ** 2, axis=0)))
    colmax_meas = float(np.max(np.sum(np.abs(s2) ** 2, axis=0)))
    say(
        f"in-model 2x2 ({[meas[p]['solve_kind'] for p in order[:2]]}, "
        f"{time.perf_counter() - t_im:.2f} s):"
    )
    for name, m in (("PRED", pred2), ("MEAS", s2)):
        say(f"    {name} " + "; ".join(", ".join(f"{complex(v):.10e}" for v in row) for row in m))
    say(
        f"(identity) residual ||MEAS - PRED||/||PRED|| = {r_id:.6e} (band {REDUCTION_BAND:.0e}, "
        f"ASSERTED); in-model |S21| = {abs(s2[1, 0]):.6e}, |S11| = {abs(s2[0, 0]):.6e}"
    )
    say(
        f"reciprocity 2x2: circuit {rec_pred:.6e}, in-model {rec_meas:.6e} (band "
        f"{RECIPROCITY_BAND:.0e}, ASSERTED); max column sum|S|^2 circuit {colmax_pred:.9f}, "
        f"in-model {colmax_meas:.9f} (ceiling {COLUMN_PASSIVITY_CEILING}, ASSERTED)"
    )
    say(
        f"negative control (PREDICTED >= {NEGATIVE_CONTROL_PREDICTED_BAND_MULTIPLE:g}x band, "
        f"printed, never asserted — rule (e)): {NEGATIVE_CONTROL_FACTOR} x C_tuned prediction "
        f"vs in-model residual {r_ctl:.6e} = {r_ctl / REDUCTION_BAND:.2f}x band "
        f"({'held' if r_ctl >= NEGATIVE_CONTROL_PREDICTED_BAND_MULTIPLE * REDUCTION_BAND else 'NOT held'})"
    )
    say(
        f"window elapsed {time.perf_counter() - t_window:.1f} s at -n {comm.size}; summed "
        f"ru_maxrss {_rss_summed_gib(comm):.3f} GiB"
    )

    # ---- asserts, in the item's order ----------------------------------------
    assert abs(rel_cells) < CELL_COUNT_BAND, f"cells {n_cells}, relative {rel_cells:.3e}"
    assert tuned["im_rel"] <= TUNING_IM_Z_RTOL
    assert (
        r_id <= REDUCTION_BAND
    ), f"identity residual {r_id:.6e} > {REDUCTION_BAND:.0e}: known-issues, never re-band"
    assert ratio32 <= RECIPROCITY_BAND
    assert rec_pred <= RECIPROCITY_BAND
    assert rec_meas <= RECIPROCITY_BAND
    assert sigma32 <= MATRIX_PASSIVITY_CEILING
    assert colmax_pred <= COLUMN_PASSIVITY_CEILING
    assert colmax_meas <= COLUMN_PASSIVITY_CEILING


# ===========================================================================
# `WF-10` step 1c (§9 item 128, chain H8a): the C-sweep C_tuned -> C_mode on
# H3's 32×32, the circuit <-> in-model identity per point.  Executes the
# "Resolves with" of known-issues 2026-09-28 (H8 0.070 % at C_tuned, item 62
# 0.349 % at C_mode — spike or trend?).  H8's configuration unchanged: P17 + P21
# in 50 Ω, 30 gaps in C.  W1 (`FEM_EM_WF10_STEP1C=1`) builds the κ-corrected
# 32×32, saves it to an npz, asserts the anchors (A1)/(A2), the negative
# control, and prints the dense circuit scan; W2 (`FEM_EM_WF10_STEP1C_W2=1`
# as well) loads the npz, asserts the guard (A3), and solves the 4 rule
# points.  The reading rule is printed only — no band moves in a slot.
# ===========================================================================

STEP1C_ENV = "FEM_EM_WF10_STEP1C"
STEP1C_W2_ENV = "FEM_EM_WF10_STEP1C_W2"
STEP1C_NPZ = "/workspace/logs/wf10/h8a-held-s32.npz"
# H8's C_tuned (`20260927T095735Z_WF-10-step1.log:10605`) and item 62's C_mode
# (`20260929T020339Z_WF-10-step1b.log:10642`), entered as literals.
C_TUNED_F = 1.973268440567e-11
C_MODE_F = 1.131080769863e-10
# (A1) H8's in-model MEAS at C_tuned, `…WF-10-step1.log:10609`; residual `:10610`.
A1_MEAS = np.array(
    [
        [-5.7176667839e-01 + 3.4507454728e-05j, -1.4806852179e-02 + 7.9149605593e-02j],
        [-1.4806852179e-02 + 7.9149605593e-02j, -5.7179555368e-01 - 3.4742542655e-05j],
    ],
    dtype=np.complex128,
)
A1_RESIDUAL_RECORD = 6.992446e-04
# (A2) item 62's in-model MEAS at C_mode, `…WF-10-step1b.log:10647`; residual `:10648`.
A2_MEAS = np.array(
    [
        [-2.0700978050e-01 + 3.9255129634e-01j, -7.2474140243e-02 - 1.2747203159e-01j],
        [-7.2474140243e-02 - 1.2747203159e-01j, -2.0774736260e-01 + 3.9176561124e-01j],
    ],
    dtype=np.complex128,
)
A2_RESIDUAL_RECORD = 3.486500e-03
# H8's cond(I - Γ S_tt) at C_tuned, `…step1.log:10606` (printed beside).
COND_AT_C_TUNED_RECORD = 4.608715
# Negative control: H8's measured 1364× (`…step1.log:10612`), asserted ≥ 10×.
STEP1C_NEGATIVE_CONTROL_MULTIPLE = 10.0
# The all-C anti-mode capacitances k = 2…8, `…step1b.log:10625–10637` (printed).
ANTI_MODE_C_PF = {2: 75.45, 3: 56.17, 4: 45.52, 5: 39.23, 6: 35.51, 7: 33.52, 8: 32.89}
SCAN_POINTS = 2001
SPIKE_RATIO = 3.0
TREND_RATIO_CEILING = 2.0


def _step1c_mode():
    if os.environ.get(STEP1C_W2_ENV, "").strip() == "1":
        return "W2"
    if os.environ.get(STEP1C_ENV, "").strip() == "1":
        return "W1"
    return None


def _entrywise_rel(measured, literal):
    """max_ij |M_ij − L_ij| / |L_ij| — "entrywise" for the anchors."""
    return float(np.max(np.abs(measured - literal) / np.abs(literal)))


def _circuit_scan(s32, frequency_hz, grid, order, slots, z0):
    """``cond(I − Γ S_tt)``, ``|S11|``, circuit k=±1 shares per ring on ``grid``."""
    from tests.validation.test_port_circuit_layer_field import (
        _capacitor_impedance,
        tuned_input,
    )

    s_tt = s32[2:, 2:]
    eye = np.eye(s_tt.shape[0])
    conds, s11s, shp, shm = [], [], [], []
    for c in grid:
        z_c = _capacitor_impedance(frequency_hz, c)
        gamma = (z_c - z0) / (z_c + z0)
        conds.append(float(np.linalg.cond(eye - gamma * s_tt)))
        _, red = tuned_input(s32, frequency_hz, c, terminated=tuple(range(2, 32)))
        s11s.append(abs(complex(red[0, 0])))
        sh = _k1_shares(_circuit_currents(s32, frequency_hz, c, 2, z0), order, slots)
        shp.append(sh[1.0])
        shm.append(sh[-1.0])
    return np.array(conds), np.array(s11s), np.array(shp), np.array(shm)


def _rule_points(grid, conds):
    """On = the two largest interior local maxima of cond (sorted by C); off =
    the cond minimum strictly between them, and strictly between the upper one
    and C_mode (the grid's last point excluded)."""
    n = len(grid)
    maxima = [i for i in range(1, n - 1) if conds[i] > conds[i - 1] and conds[i] >= conds[i + 1]]
    if len(maxima) < 2:
        return maxima, None
    top = sorted(sorted(maxima, key=lambda i: conds[i], reverse=True)[:2])
    lo, hi = top
    off1 = lo + 1 + int(np.argmin(conds[lo + 1 : hi]))
    off2 = hi + 1 + int(np.argmin(conds[hi + 1 : n - 1])) if hi + 1 < n - 1 else None
    return maxima, {"on_lo": lo, "on_hi": hi, "off1": off1, "off2": off2}


def _in_model_pair(ctx, base_specs, kappa, others, feeds, z_c):
    """In-model 2×2 with PORT-14 C sheets at ``z_c`` on the 30 gaps.  The first
    drive refactorises (each C changes the bilinear form); only the two feeds'
    ``s_column`` is read (C sheets' V / power read 50 Ω)."""
    from tests.validation.test_port_birdcage_ring_column import _solve_one_drive

    ctx_c = dict(ctx)
    ctx_c["specs"] = _respec(base_specs, kappa, {pid: z_c for pid in others})
    meas, kinds = {}, []
    for n_drive, pid in enumerate(feeds):
        c = _solve_one_drive(ctx_c, pid, reuse_factorization=(n_drive > 0))
        c.pop("fields")
        meas[pid] = c
        kinds.append(c["solve_kind"])
    s2 = np.array([[meas[pj]["s_column"][pi] for pj in feeds] for pi in feeds], dtype=np.complex128)
    return s2, kinds


@complex_only
def test_wf10_step1c_c_sweep():
    mode = _step1c_mode()
    if mode is None:
        pytest.skip(f"{STEP1C_ENV} unset: heavy -n 8 F-human C-sweep, environment-selected")

    from fem_em_solver.ports.shares import terminal_form_deficit
    from tests.mesh.test_birdcage_port_sheet_prerequisite import CELL_COUNT_BAND
    from tests.validation.test_birdcage_f_human_rung import F_HUMAN_BRANCH_B_CELL_RECORD
    from tests.validation.test_port_birdcage_four_port import LEG_D0_REPRODUCTION_BAND
    from tests.validation.test_port_birdcage_ring_column import (
        COLUMN_PASSIVITY_CEILING,
        _reciprocity_ratio,
        _solve_one_drive,
    )
    from tests.validation.test_port_birdcage_ring_matrix import MATRIX_PASSIVITY_CEILING
    from tests.validation.test_port_circuit_layer_field import (
        _capacitor_impedance,
        _residual,
        tuned_input,
    )
    from tests.validation.test_port_drive_superposition import _ring_quadrature_slots
    from tests.validation.test_port_lumped_rlc_termination import REDUCTION_BAND
    from tests.validation.test_port_lumped_sheet_sweep import RECIPROCITY_BAND
    from tests.validation.test_port_package_sparameters import REFERENCE_IMPEDANCE_OHM
    from tests.validation.test_wf7_f_human_b1_quadrature import (
        DEGREE,
        FREQUENCY_HZ,
        _build_f_human,
    )

    comm = MPI.COMM_WORLD
    z0 = float(REFERENCE_IMPEDANCE_OHM)
    terminated = tuple(range(2, 32))
    t_window = time.perf_counter()

    def say(msg: str) -> None:
        if comm.rank == 0:
            print(f"[WF-10 step1c {mode}] {msg}", flush=True)

    say(f"start: -n {comm.size}, f = {FREQUENCY_HZ:.3e} Hz, degree {DEGREE}")

    if mode == "W2":
        comm.Barrier()
        if not os.path.exists(STEP1C_NPZ):
            pytest.fail(f"{STEP1C_NPZ} missing: W2 never recomputes the 32x32 — stop")
        saved = np.load(STEP1C_NPZ, allow_pickle=False)
        say(f"loaded {STEP1C_NPZ}")

    # ---- fixture -----------------------------------------------------------
    t0 = time.perf_counter()
    built = _build_f_human(comm, emit_leg_ports=True)  # GEO-35 step 2b: unpin with the re-record
    n_cells = built["cells"]
    rel_cells = n_cells / F_HUMAN_BRANCH_B_CELL_RECORD - 1.0
    say(
        f"cells {n_cells} vs record {F_HUMAN_BRANCH_B_CELL_RECORD}: relative {rel_cells:.3e} "
        f"(band {CELL_COUNT_BAND}, ASSERTED); build {time.perf_counter() - t0:.2f} s"
    )
    assert abs(rel_cells) < CELL_COUNT_BAND
    ctx = built["ctx"]
    base_specs = list(ctx["specs"])
    az_ref, slots = _ring_quadrature_slots(built["sheets"])
    ref_pid = "P17"
    ref_k, ref_sign = slots[ref_pid]
    feed2 = next(
        pid
        for pid, (k, sg) in slots.items()
        if sg == ref_sign and k == (ref_k + QUADRATURE_SLOT_OFFSET) % 16
    )
    others = sorted((pid for pid in slots if pid not in (ref_pid, feed2)), key=lambda p: int(p[1:]))
    order = [ref_pid, feed2] + others
    feeds = order[:2]
    say(f"feeds: {ref_pid} and {feed2}; 30 ports in C")

    checks = []  # (label, pred, meas) for reciprocity / passivity

    if mode == "W1":
        # ---- κ, then the κ-corrected 32×32 (step 1's route verbatim) -------
        col = _solve_one_drive(ctx, ref_pid, reuse_factorization=False)
        fields = col.pop("fields")
        sheets = [spec.sheet(driven=(spec.port_id == ref_pid)) for spec in base_specs]
        kappa = float(
            terminal_form_deficit(ctx["msh"], ctx["tags_f"], sheets, fields.e_complex, comm)[
                "pooled"
            ]
        )
        del fields
        say(f"kappa = {kappa:.6e}; t={time.perf_counter() - t_window:.1f} s")
        ctx32 = dict(ctx)
        ctx32["specs"] = _respec(base_specs, kappa)
        columns = {}
        for n_drive, pid in enumerate(order):
            c = _solve_one_drive(ctx32, pid, reuse_factorization=(n_drive > 0))
            c.pop("fields")
            columns[pid] = c
        kinds = [columns[p]["solve_kind"] for p in order]
        s32 = np.array(
            [[columns[pj]["s_column"][pi] for pj in order] for pi in order], dtype=np.complex128
        )
        say(
            f"32x32: {kinds.count('solve')} factorising, {kinds.count('held')} held; "
            f"t={time.perf_counter() - t_window:.1f} s"
        )
    else:
        kappa = float(saved["kappa"])
        s32 = np.array(saved["s32"], dtype=np.complex128)
        saved_order = [str(p) for p in saved["order"]]
        assert saved_order == order, f"port order differs from the npz: {saved_order}"

    ratio32 = _reciprocity_ratio(s32)
    sigma32 = float(np.linalg.norm(s32, 2))
    say(
        f"32x32 reciprocity {ratio32:.6e} (band {RECIPROCITY_BAND:.0e}, ASSERTED); sigma_max "
        f"{sigma32:.9f} (ceiling {MATRIX_PASSIVITY_CEILING}, ASSERTED)"
    )
    assert ratio32 <= RECIPROCITY_BAND
    assert sigma32 <= MATRIX_PASSIVITY_CEILING

    def pair_at(c_f, label):
        z_c = _capacitor_impedance(FREQUENCY_HZ, c_f)
        t_p = time.perf_counter()
        meas, kinds = _in_model_pair(ctx, base_specs, kappa, others, feeds, z_c)
        _, pred = tuned_input(s32, FREQUENCY_HZ, c_f, terminated=terminated)
        r = _residual(meas, pred)
        gamma = (z_c - z0) / (z_c + z0)
        cond = float(np.linalg.cond(np.eye(30) - gamma * s32[2:, 2:]))
        say(
            f"{label}: C = {c_f:.12e} F, {kinds}, {time.perf_counter() - t_p:.2f} s; "
            f"t={time.perf_counter() - t_window:.1f} s"
        )
        for name, m in (("PRED", pred), ("MEAS", meas)):
            say(
                f"    {name} "
                + "; ".join(", ".join(f"{complex(v):.10e}" for v in row) for row in m)
            )
        say(f"    residual {r:.6e}  cond {cond:.6e}  residual/cond {r / cond:.6e}")
        checks.append((label, pred, meas))
        return meas, pred, r, cond

    # ---- dense circuit scan (the rule points come from it in both windows) --
    grid = np.geomspace(C_TUNED_F, C_MODE_F, SCAN_POINTS)
    conds, s11s, shp, shm = _circuit_scan(s32, FREQUENCY_HZ, grid, order, slots, z0)
    maxima, pts = _rule_points(grid, conds)
    say(
        f"scan: {SCAN_POINTS} log-spaced C in [C_tuned, C_mode]; cond at C_tuned "
        f"{conds[0]:.6e} beside H8's {COND_AT_C_TUNED_RECORD} (printed); {len(maxima)} interior "
        f"local maxima of cond"
    )
    if mode == "W1":
        say("    i     C [pF]        cond(I-GS_tt)   |S11|      share ring+  share ring-")
        extrema = set(maxima)
        extrema |= {
            i
            for i in range(1, SCAN_POINTS - 1)
            if conds[i] < conds[i - 1] and conds[i] <= conds[i + 1]
        }
        for i in range(SCAN_POINTS):
            if i % 10 == 0 or i in extrema or i == SCAN_POINTS - 1:
                tag = " MAX" if i in maxima else (" min" if i in extrema else "")
                say(
                    f"    {i:4d}  {grid[i] * 1e12:12.6f}  {conds[i]:14.6e}  {s11s[i]:.6f}  "
                    f"{shp[i]:.6f}     {shm[i]:.6f}{tag}"
                )
    say(
        "    all-C anti-mode capacitances (printed beside, …step1b.log:10625–10637): "
        + ", ".join(f"k={k} {v} pF" for k, v in ANTI_MODE_C_PF.items())
    )
    if pts is not None:
        for key in ("on_lo", "on_hi", "off1", "off2"):
            i = pts[key]
            if i is None:
                say(f"    rule point {key}: NONE")
            else:
                say(
                    f"    rule point {key}: i={i} C = {grid[i]:.12e} F ({grid[i] * 1e12:.4f} pF) "
                    f"cond {conds[i]:.6e} |S11| {s11s[i]:.6f} shares +{shp[i]:.6f} "
                    f"-{shm[i]:.6f}"
                )
    else:
        say(f"    fewer than 2 interior local maxima of cond ({maxima}): no rule points")

    if mode == "W1":
        meas_t, pred_t, r_t, cond_t = pair_at(C_TUNED_F, "(A1) C_tuned")
        a1 = _entrywise_rel(meas_t, A1_MEAS)
        say(
            f"(A1) MEAS vs …step1.log:10609 entrywise {a1:.3e} (band "
            f"{LEG_D0_REPRODUCTION_BAND:.0e}, ASSERTED); residual {r_t:.6e} beside "
            f"{A1_RESIDUAL_RECORD:.6e}"
        )
        _, pred_ctl = tuned_input(
            s32, FREQUENCY_HZ, NEGATIVE_CONTROL_FACTOR * C_TUNED_F, terminated=terminated
        )
        r_ctl = _residual(meas_t, pred_ctl)
        say(
            f"negative control: {NEGATIVE_CONTROL_FACTOR} x C_tuned prediction vs in-model at "
            f"C_tuned residual {r_ctl:.6e} = {r_ctl / REDUCTION_BAND:.2f}x band (ASSERTED >= "
            f"{STEP1C_NEGATIVE_CONTROL_MULTIPLE:g}x; H8 measured 1364x)"
        )
        meas_m, pred_m, r_m, cond_m = pair_at(C_MODE_F, "(A2) C_mode")
        a2 = _entrywise_rel(meas_m, A2_MEAS)
        say(
            f"(A2) MEAS vs …step1b.log:10647 entrywise {a2:.3e} (band "
            f"{LEG_D0_REPRODUCTION_BAND:.0e}, ASSERTED); residual {r_m:.6e} beside "
            f"{A2_RESIDUAL_RECORD:.6e}"
        )
        if comm.rank == 0:
            os.makedirs(os.path.dirname(STEP1C_NPZ), exist_ok=True)
            np.savez(
                STEP1C_NPZ,
                s32=s32,
                kappa=kappa,
                order=np.array(order),
                slot_k=np.array([slots[p][0] for p in order]),
                slot_sign=np.array([slots[p][1] for p in order]),
                cells=n_cells,
                meas_tuned=meas_t,
                r_tuned=r_t,
                cond_tuned=cond_t,
                meas_mode=meas_m,
                r_mode=r_m,
                cond_mode=cond_m,
                c_tuned=C_TUNED_F,
                c_mode=C_MODE_F,
            )
            say(f"saved {STEP1C_NPZ}")
        comm.Barrier()
    else:
        cells_saved = int(saved["cells"])
        say(f"(A3) cells {n_cells} vs npz {cells_saved} (ASSERTED equal)")
        assert n_cells == cells_saved, f"(A3) cells {n_cells} != npz {cells_saved}"
        meas_t, pred_t, r_t, cond_t = pair_at(C_TUNED_F, "(A3) guard C_tuned")
        a3 = _entrywise_rel(meas_t, np.array(saved["meas_tuned"]))
        say(
            f"(A3) re-solve vs W1 MEAS entrywise {a3:.3e} (band {LEG_D0_REPRODUCTION_BAND:.0e}, "
            f"ASSERTED)"
        )
        assert a3 <= LEG_D0_REPRODUCTION_BAND, f"(A3) {a3:.3e}: npz carry-over broken"
        r_m, cond_m = float(saved["r_mode"]), float(saved["cond_mode"])
        assert pts is not None, f"fewer than 2 interior local maxima of cond: {maxima}"
        res = {"tuned": (C_TUNED_F, r_t, cond_t), "mode": (C_MODE_F, r_m, cond_m)}
        for key in ("on_lo", "off1", "on_hi", "off2"):
            if pts[key] is None:
                continue
            _, _, r, cond = pair_at(float(grid[pts[key]]), f"rule {key}")
            res[key] = (float(grid[pts[key]]), r, cond)

        # ---- reading rule (printed only; no band moves in a slot) -----------
        say("reading table (C [pF], residual, cond, residual/cond):")
        for key in sorted(res, key=lambda k: res[k][0]):
            c, r, cd = res[key]
            say(f"    {key:6s} {c * 1e12:10.4f}  {r:.6e}  {cd:.6e}  {r / cd:.6e}")
        # Bracketing: on_lo between the C_tuned anchor and off1; on_hi between
        # off1 and off2 (C_mode if off2 is absent).
        brackets = {
            "on_lo": ("tuned", "off1"),
            "on_hi": ("off1", "off2" if "off2" in res else "mode"),
        }
        ratios = {}
        for on, (a, b) in brackets.items():
            ratios[on] = (res[on][1] / res[a][1], res[on][1] / res[b][1])
            say(f"    on/off {on}: /{a} {ratios[on][0]:.3f}, /{b} {ratios[on][1]:.3f}")
        spike = all(min(v) >= SPIKE_RATIO for v in ratios.values())
        trend_keys = sorted(
            (k for k in ("tuned", "off1", "off2", "mode") if k in res), key=lambda k: res[k][0]
        )
        trend_r = [res[k][1] for k in trend_keys]
        d = np.diff(trend_r)
        monotone = bool(np.all(d > 0) or np.all(d < 0))
        trend = monotone and all(max(v) < TREND_RATIO_CEILING for v in ratios.values())
        verdict = "spike" if spike else ("trend" if trend else "mixed")
        cs = np.array([res[k][0] for k in res])
        rs = np.array([res[k][1] for k in res])
        slope_all = float(np.polyfit(np.log(cs), np.log(rs), 1)[0])
        slope_trend = float(
            np.polyfit(np.log([res[k][0] for k in trend_keys]), np.log(trend_r), 1)[0]
        )
        say(
            f"reading (printed only): monotone over anchors+off {monotone} ({trend_keys}); "
            f"classification = {verdict.upper()}; log-log slope residual vs C: all 6 points "
            f"{slope_all:+.4f}, anchors+off-points {slope_trend:+.4f}"
        )

    for label, pred, meas in checks:
        rp, rm = _reciprocity_ratio(pred), _reciprocity_ratio(meas)
        cp = float(np.max(np.sum(np.abs(pred) ** 2, axis=0)))
        cm = float(np.max(np.sum(np.abs(meas) ** 2, axis=0)))
        say(
            f"{label}: reciprocity circuit {rp:.3e} in-model {rm:.3e} (band "
            f"{RECIPROCITY_BAND:.0e}, ASSERTED); col sum|S|^2 circuit {cp:.9f} in-model "
            f"{cm:.9f} (ceiling {COLUMN_PASSIVITY_CEILING}, ASSERTED)"
        )
    say(
        f"window elapsed {time.perf_counter() - t_window:.1f} s at -n {comm.size}; summed "
        f"ru_maxrss {_rss_summed_gib(comm):.3f} GiB"
    )

    if mode == "W1":
        assert a1 <= LEG_D0_REPRODUCTION_BAND, f"(A1) {a1:.3e}: reproducibility finding — stop"
        assert a2 <= LEG_D0_REPRODUCTION_BAND, f"(A2) {a2:.3e}: reproducibility finding — stop"
        assert r_ctl >= STEP1C_NEGATIVE_CONTROL_MULTIPLE * REDUCTION_BAND
    for label, pred, meas in checks:
        assert _reciprocity_ratio(pred) <= RECIPROCITY_BAND, label
        assert _reciprocity_ratio(meas) <= RECIPROCITY_BAND, label
        assert float(np.max(np.sum(np.abs(pred) ** 2, axis=0))) <= COLUMN_PASSIVITY_CEILING, label
        assert float(np.max(np.sum(np.abs(meas) ** 2, axis=0))) <= COLUMN_PASSIVITY_CEILING, label
