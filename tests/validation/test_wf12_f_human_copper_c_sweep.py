"""`WF-12` step 2′ (§10 chain H, H12b′, §9 item 175) — H8a's C-sweep
construction (`test_wf10_f_human_tuned_circuit.py:292–650`) on H12b's copper
F-human hole fixture (`test_wf12_f_human_copper_matrix.py`).

Selected by environment only: ``FEM_EM_WF12_STEP2P=1`` runs W1;
``FEM_EM_WF12_STEP2P_W2=1`` as well runs W2.  Brief:
docs/planning/chunks/WF-12.md "### §9 item 175".

Configuration: H8a's — P17 and P21 in 50 Ω, the other 30 ring gaps in C,
**uncorrected** sheets (exact CAD rectangles, no κ width correction — operator
ruling 2026-10-02, PROJECT_PLAN §9 "Binding ruling — sheets and κ"; the 12:00
slot's κ-corrected run went red at (M) with σ_max 1.000002129, parked on
`attempt/WF-12-step2p-20261002T172400Z`), copper (Leontovich, σ = 5.8e7) cavity wall first in every
factorising call.  Every in-model pair refactorises at P17 (never a factor held
across a C point or the κ change); only P21 is held.

Anchors (asserted, every band imported): (A0) hole cells + (iv) census;
(R) uncorrected Cu P17 reproduces H12b's P_src / P_surf; (P) the power
identity on that drive; (M) the uncorrected Cu 32×32 reciprocity / σ_max;
(C) per-C circuit and in-model 2×2 reciprocity / column passivity;
(A3, W2) guard MEAS at C_tuned equals W1's npz entrywise.
Negative control asserted: identity with P_surf dropped misses ≥ 10× band.
The per-C identity residual is a record, printed only, with the
pre-registered reading rule (printed only).

κ = terminal_form_deficit(...)["pooled"] on the uncorrected P17 drive is
computed and printed as a diagnostic only, never applied (ruling 2026-10-02).
The σ = 800 comparand residuals (`…step1c-W2.log`) were κ-corrected, so
ρ = r_Cu / r_800 compares copper-uncorrected against σ = 800 κ-corrected; the
caveat is printed and the reading rule prints as pre-registered.

Scope: degree 1, 64 MHz, H8a's two-feed configuration, this one fixture.
"""

from __future__ import annotations

import os
import time

import numpy as np
import pytest
from mpi4py import MPI

from tests.complex_mode import complex_only

STEP2P_ENV = "FEM_EM_WF12_STEP2P"
STEP2P_W2_ENV = "FEM_EM_WF12_STEP2P_W2"
STEP2P_NPZ = "/workspace/logs/wf12/h12b-prime-held-s32.npz"
H8A_NPZ = "/workspace/logs/wf10/h8a-held-s32.npz"
DRIVEN = "P17"

# (A0) record and (R) records: H12b W2 (`20261002T094751Z_WF-12-step2-W2.log:6352,6387`)
# are imported from H12b's module (HOLE_CELL_RECORD, P17_P_SRC_RECORD_W,
# P17_P_SURF_RECORD_W); S_driven(P17) printed beside `:6411`.
S_DRIVEN_H12B = "0.383271+0.391623j"

# σ = 800 literals, `20261001T110957Z_WF-10-step1c-W2.log` (C at :10625–10637,
# residual and cond at :10642–10647).  Keys: label -> (C [F] or None for the
# imported constant, residual, cond).
SIGMA800 = {
    "tuned": (None, 6.992446e-04, 4.608715),
    "on_lo": (7.569773860279e-11, 2.016987e-03, 23.12163),
    "off1": (8.834726608692e-11, 3.576879e-03, 12.51408),
    "on_hi": (1.050184374462e-10, 6.106169e-03, 21.35226),
    "off2": (1.130093729141e-10, 3.513256e-03, 18.61198),
    "mode": (None, 3.486500e-03, 18.56031),
}
H8B_RESIDUAL = 2.399380e-03  # `20261001T140319Z_WF-10-step1d-W1.log:10658` — NOT the comparand
SIGMA800_SLOPE = 1.0332  # `…step1c-W2.log:10650`
SIGMA800_KAPPA = 4.590991e-03  # `…step1c-W1.log:10611`
SIGMA800_C_TUNED_PF = 19.7327
SIGMA800_COND_MAX_PF = (75.70, 105.02)
SIGMA800_NEG_MULTIPLE = 1364.0  # `…step1c-W1.log:10843`
SIGMA800_ROUTE_PRECEDENT = 9.1e-15  # `…step1d-W1.log:10654`
H8B_TRANSFER = 0.688

# Predictions written first (§9 item 175), printed only.
PRED_KAPPA = (2.3e-3, 9.2e-3)
PRED_RECIP = 1e-14
PRED_SMAX_MARGIN = (1e-7, 1e-4)
PRED_C_TUNED_PF = (17.0, 23.0)
PRED_COND_MAX_RTOL = 0.05
PRED_NEG = (300.0, 3000.0)
PRED_RHO = (0.5, 2.0)
PRED_R_MODE = (1.7e-3, 7e-3)
PRED_SLOPE = (1.03, 0.3)
PRED_PAIR_S = (35.0, 55.0)
PRED_RSS_GIB = 13.0
# Reading rule thresholds (pre-registered, printed only).
RULE_CONDUCTOR_R_MODE = 1e-3
RULE_CONDUCTOR_RHO = 0.5
RULE_CONDUCTOR_C_MIN_F = 75e-12
RULE_TRACKS_RHO = (0.5, 2.0)
RULE_TRACKS_SLOPE = (1.03, 0.3)
ON_POINT_C_RTOL = 0.005
W2_PRICE_CEILING_S = 540.0
CURRENT_ROUTE_CONVERSION_BAND = 1e-12  # printed beside, not asserted (brief)
# Operator ruling 2026-10-02 (PROJECT_PLAN §9 binding rulings): port / lumped
# sheets are exact CAD rectangles and the κ width correction is NOT applied to
# them.  Every spec built here (32×32, in-model pairs) uses this value; κ is a
# printed diagnostic only.
KAPPA_APPLIED = None
H12B_UNCORRECTED_SIGMA_MAX = 0.999993339  # H12b uncorrected Cu 32×32 (ring order), printed beside


def _mode():
    if os.environ.get(STEP2P_W2_ENV, "").strip() == "1":
        return "W2"
    if os.environ.get(STEP2P_ENV, "").strip() == "1":
        return "W1"
    return None


def _drive(ctx, built, specs, driven_id, *, held, wall_term, z_s, partition, keep_e):
    """Module-level local copy of H12b's closure ``_drive``
    (`test_wf12_f_human_copper_matrix.py:170–217`), taking ``specs`` as an
    argument; ``partition=False`` skips the volume / surface / sheet power
    forms; ``keep_e=True`` returns ``e`` (for κ) before ``del fields``.  The
    wall term goes first in the factorising call's bilinear list."""
    import ufl
    from dolfinx import fem

    from fem_em_solver.ports.lumped import (
        lumped_port_bilinear_term,
        lumped_port_linear_term,
        sheet_terminal_current,
    )
    from tests.validation.test_birdcage_power_identity import (
        _sheet_field_dissipation_w,
        _source_power_w,
    )
    from tests.validation.test_port_birdcage_four_port import TERMINATED_PORT_IMPEDANCE_OHM
    from tests.validation.test_port_birdcage_lumped_column import PHANTOM_CELL_TAG
    from tests.validation.test_th14_birdcage_copper import _surface_loss_w

    comm = MPI.COMM_WORLD
    msh, tags_f, omega, cell_tags = ctx["msh"], ctx["tags_f"], ctx["omega"], ctx["cell_tags"]
    solver = ctx["solver"]
    z_p = complex(TERMINATED_PORT_IMPEDANCE_OHM)
    sheets = [sp.sheet(driven=(sp.port_id == driven_id)) for sp in specs]
    dsheet = next(s for s in sheets if s.port_id == driven_id)
    linear = [
        lambda test, _s=dsheet: lumped_port_linear_term(
            msh, tags_f, _s, test, omega_rad_per_s=omega
        )
    ]
    comm.Barrier()
    t0 = time.perf_counter()
    if held:
        fields = solver.solve_with_held_factorization(extra_linear_terms=linear)
    else:
        bil = ([wall_term] if wall_term is not None else []) + [
            (
                lambda trial, test, _s=s: lumped_port_bilinear_term(
                    msh, tags_f, _s, trial, test, omega_rad_per_s=omega
                )
            )
            for s in sheets
        ]
        fields = solver.solve(
            current_density=None,
            project_source=False,
            extra_bilinear_terms=bil,
            extra_linear_terms=linear,
        )
    comm.Barrier()
    t_solve = time.perf_counter() - t0
    t1 = time.perf_counter()
    e = fields.e_complex
    currents = {s.port_id: sheet_terminal_current(msh, tags_f, s, e, comm) for s in sheets}
    v_src = complex(dsheet.source_voltage_v)
    s_col = {}
    for s in sheets:
        v = complex(s.source_voltage_v) - currents[s.port_id] * z_p
        s_col[s.port_id] = (v - z_p * currents[s.port_id]) / v_src
    out = {"s": s_col, "kind": "held" if held else "solve", "t_solve": t_solve, "sheets": sheets}
    if partition:
        sigma_f = fields.sigma_field
        dx = ufl.Measure("dx", domain=msh, subdomain_data=cell_tags)
        acct = {"mesh": msh, "facet_tags": tags_f}

        def _vol(measure):
            form = fem.form(0.5 * sigma_f * ufl.inner(e, e) * measure)
            return float(np.real(comm.allreduce(fem.assemble_scalar(form), op=MPI.SUM)))

        out["p_sheets"] = float(sum(_sheet_field_dissipation_w(acct, s, e, omega) for s in sheets))
        out["p_src"] = float(_source_power_w(acct, dsheet, e, omega))
        p_total = _vol(dx)
        out["p_phantom"] = _vol(dx(PHANTOM_CELL_TAG))
        out["p_rest"] = p_total - out["p_phantom"]
        out["p_surf"] = _surface_loss_w(msh, built["facets"], z_s, e, comm)
    if keep_e:
        out["e"] = e
    del fields
    out["t_post"] = time.perf_counter() - t1
    return out


@complex_only
def test_wf12_step2p_copper_c_sweep():
    mode = _mode()
    if mode is None:
        pytest.skip(f"{STEP2P_ENV} unset: heavy -n 8 copper C-sweep, environment-selected")

    from fem_em_solver.core.cavity import surface_resistance_ohm
    from fem_em_solver.ports.shares import terminal_form_deficit
    from tests.mesh.test_birdcage_port_sheet_prerequisite import CELL_COUNT_BAND
    from tests.validation.test_birdcage_f_human_rung import _ring_ports
    from tests.validation.test_birdcage_power_identity import DISCRETE_IDENTITY_RTOL
    from tests.validation.test_port_birdcage_four_port import LEG_D0_REPRODUCTION_BAND
    from tests.validation.test_port_birdcage_ring_column import (
        COLUMN_PASSIVITY_CEILING,
        _reciprocity_ratio,
    )
    from tests.validation.test_port_birdcage_ring_matrix import (
        MATRIX_PASSIVITY_CEILING,
        REUSE_REPRODUCTION_RTOL,
    )
    from tests.validation.test_port_circuit_layer_field import (
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
    from tests.validation.test_th14_birdcage_copper import _leontovich_term
    from tests.validation.test_wf7_f_human_b1_quadrature import (
        DEGREE,
        FREQUENCY_HZ,
        _build_f_human,
    )
    from tests.validation.test_wf10_f_human_tuned_circuit import (
        ANTI_MODE_C_PF,
        C_MODE_F,
        C_TUNED_F,
        NEGATIVE_CONTROL_FACTOR,
        QUADRATURE_SLOT_OFFSET,
        SCAN_POINTS,
        _circuit_currents,
        _circuit_scan,
        _entrywise_rel,
        _k1_shares,
        _respec,
        _rule_points,
    )
    from tests.validation.test_wf12_f_human_copper_matrix import (
        BLIND_FACTOR,
        COPPER_SIGMA,
        HOLE_CELL_RECORD,
        P17_P_SRC_RECORD_W,
        P17_P_SURF_RECORD_W,
    )
    from tests.validation.test_wf12_f_human_hole import (
        _census_checks,
        _owned_census,
        _rss_gib,
    )

    comm = MPI.COMM_WORLD
    z0 = float(REFERENCE_IMPEDANCE_OHM)
    terminated = tuple(range(2, 32))
    t_window = time.perf_counter()
    sig800 = dict(SIGMA800)
    sig800["tuned"] = (C_TUNED_F,) + SIGMA800["tuned"][1:]
    sig800["mode"] = (C_MODE_F,) + SIGMA800["mode"][1:]

    def say(msg):
        if comm.rank == 0:
            print(f"[WF-12 step2p {mode}] {msg}", flush=True)

    say(f"start: -n {comm.size}, f = {FREQUENCY_HZ:.3e} Hz, degree {DEGREE}")
    saved = None
    if mode == "W2":
        comm.Barrier()
        if not os.path.exists(STEP2P_NPZ):
            pytest.fail(f"{STEP2P_NPZ} missing: W2 never recomputes it — stop")
        saved = np.load(STEP2P_NPZ, allow_pickle=False)
        say(f"loaded {STEP2P_NPZ}")

    # ---- fixture, (A0), (iv) -------------------------------------------------
    t0 = time.perf_counter()
    built = _build_f_human(comm, as_hole=True)
    t_build = time.perf_counter() - t0
    ctx = built["ctx"]
    msh, cell_tags, omega = ctx["msh"], ctx["cell_tags"], ctx["omega"]
    census = _owned_census(msh, cell_tags, comm)
    census_fails = _census_checks(census, built["cells"], _ring_ports())
    oc = built["outer_census"]
    n_cells = built["cells"]
    cell_rel = n_cells / HOLE_CELL_RECORD - 1.0
    say(
        f"(A0) hole cells {n_cells} vs record {HOLE_CELL_RECORD}, rel {cell_rel:.3e} "
        f"(band {CELL_COUNT_BAND}, ASSERTED; solid sigma=800 route 507266 — binary control); "
        f"build {t_build:.1f} s; (iv) failures {census_fails or 'none'}; outer census {oc}"
    )
    assert abs(cell_rel) <= CELL_COUNT_BAND, cell_rel
    assert not census_fails, census_fails
    assert oc["tag_collision"] == 0 and oc["cavity"] > 0, oc
    assert (
        oc["cavity_exterior"] == oc["cavity"] and oc["outer"] + oc["cavity"] == oc["exterior"]
    ), oc

    z_s = (1.0 + 1.0j) * surface_resistance_ohm(omega, COPPER_SIGMA)
    wall = _leontovich_term(msh, built["facets"], z_s, omega)
    base_specs = list(ctx["specs"])
    az_ref, slots = _ring_quadrature_slots(built["sheets"])
    ref_k, ref_sign = slots[DRIVEN]
    feed2 = next(
        pid
        for pid, (k, sg) in slots.items()
        if sg == ref_sign and k == (ref_k + QUADRATURE_SLOT_OFFSET) % 16
    )
    others = sorted((pid for pid in slots if pid not in (DRIVEN, feed2)), key=lambda p: int(p[1:]))
    order = [DRIVEN, feed2] + others
    feeds = order[:2]
    assert len(order) == 32
    say(f"feeds: {DRIVEN} and {feed2} (H8a's order); 30 ports in C")

    checks = []
    pair_times = []

    def pair_at(c_f, label, kappa, s32):
        """P17 factorises (wall + 30 C sheets + 2 feeds at 50 Ω, uncorrected
        widths — κ not applied, ruling 2026-10-02), P21 held; only the feeds'
        S column is read (partition=False).  ``kappa`` is the applied value,
        always KAPPA_APPLIED (None)."""
        z_c = _capacitor_impedance(FREQUENCY_HZ, c_f)
        specs_c = _respec(base_specs, kappa, {pid: z_c for pid in others})
        t_p = time.perf_counter()
        cols, kinds = {}, []
        for n_drive, pid in enumerate(feeds):
            d = _drive(
                ctx,
                built,
                specs_c,
                pid,
                held=(n_drive > 0),
                wall_term=wall,
                z_s=z_s,
                partition=False,
                keep_e=False,
            )
            cols[pid] = d["s"]
            kinds.append(d["kind"])
        meas = np.array([[cols[pj][pi] for pj in feeds] for pi in feeds], dtype=np.complex128)
        dt = time.perf_counter() - t_p
        pair_times.append(dt)
        _, pred = tuned_input(s32, FREQUENCY_HZ, c_f, terminated=terminated)
        r = _residual(meas, pred)
        gamma = (z_c - z0) / (z_c + z0)
        cond = float(np.linalg.cond(np.eye(30) - gamma * s32[2:, 2:]))
        cur = _circuit_currents(s32, FREQUENCY_HZ, c_f, 2, z0)
        sh = _k1_shares(cur, order, slots)
        route = np.array([(1.0 if i == 0 else 0.0) - 2.0 * z0 * cur[i] for i in range(2)])
        route_diff = float(np.max(np.abs(pred[:, 0] - route)))
        say(
            f"{label}: C = {c_f:.12e} F ({c_f * 1e12:.4f} pF), {kinds}, pair {dt:.2f} s "
            f"(PREDICTED {PRED_PAIR_S} s); t={time.perf_counter() - t_window:.1f} s"
        )
        for name, m in (("PRED", pred), ("MEAS", meas)):
            say(
                f"    {name} "
                + "; ".join(", ".join(f"{complex(v):.10e}" for v in row) for row in m)
            )
        say(
            f"    residual {r:.6e}  cond {cond:.6e}  residual/cond {r / cond:.6e}  "
            f"k=+-1 shares ring+ {sh[1.0]:.6f} ring- {sh[-1.0]:.6f}"
        )
        say(
            f"    circuit routes (printed): max|tuned_input col0 - (delta - 2 z0 I)| "
            f"{route_diff:.3e} (beside {CURRENT_ROUTE_CONVERSION_BAND:g}; sigma=800 "
            f"{SIGMA800_ROUTE_PRECEDENT:g})"
        )
        checks.append((label, pred, meas))
        return meas, pred, r, cond, sh

    if mode == "W1":
        # ---- uncorrected Cu P17: (R), (P), NEG, κ ----------------------------
        d = _drive(
            ctx,
            built,
            base_specs,
            DRIVEN,
            held=False,
            wall_term=wall,
            z_s=z_s,
            partition=True,
            keep_e=True,
        )
        p_vol = d["p_phantom"] + d["p_rest"]
        rel_p = abs(d["p_src"] - d["p_sheets"] - p_vol - d["p_surf"]) / abs(d["p_src"])
        blind = abs(d["p_src"] - d["p_sheets"] - p_vol) / abs(d["p_src"])
        r_src = abs(d["p_src"] / P17_P_SRC_RECORD_W - 1.0)
        r_surf = abs(d["p_surf"] / P17_P_SURF_RECORD_W - 1.0)
        t_k = time.perf_counter()
        kappa = float(
            terminal_form_deficit(msh, ctx["tags_f"], d["sheets"], d["e"], comm)["pooled"]
        )
        del d["e"]
        say(
            f"Cu P17 uncorrected: {d['kind']} {d['t_solve']:.2f} s + post {d['t_post']:.2f} s; "
            f"S_driven {d['s'][DRIVEN]:.6f} beside H12b {S_DRIVEN_H12B} (printed)"
        )
        say(
            f"(R) P_src {d['p_src']:.9e} vs {P17_P_SRC_RECORD_W:.9e} rel {r_src:.3e}; P_surf "
            f"{d['p_surf']:.9e} vs {P17_P_SURF_RECORD_W:.9e} rel {r_surf:.3e} (band "
            f"{REUSE_REPRODUCTION_RTOL:g}, ASSERTED)"
        )
        say(
            f"(P) rel {rel_p:.3e} (band {DISCRETE_IDENTITY_RTOL:g}, ASSERTED); NEG blind "
            f"{blind:.4e} = {blind / DISCRETE_IDENTITY_RTOL:.1f}x band (bar {BLIND_FACTOR:.0f}x, "
            f"ASSERTED; ceiling ≈124x)"
        )
        say(
            f"kappa_Cu = {kappa:.6e} beside sigma=800 {SIGMA800_KAPPA:.6e} (PREDICTED "
            f"{PRED_KAPPA}); eval {time.perf_counter() - t_k:.2f} s — DIAGNOSTIC ONLY, not "
            f"applied (operator ruling 2026-10-02: exact-rectangle sheets, kappa off)"
        )
        assert r_src <= REUSE_REPRODUCTION_RTOL, r_src
        assert r_surf <= REUSE_REPRODUCTION_RTOL, r_surf
        assert rel_p <= DISCRETE_IDENTITY_RTOL, rel_p
        assert blind >= BLIND_FACTOR * DISCRETE_IDENTITY_RTOL, blind
        del d

        # ---- uncorrected Cu 32×32 in H8a's order (κ not applied, 2026-10-02) --
        specs32 = _respec(base_specs, KAPPA_APPLIED)
        cols, kinds, cad = {}, [], []
        for n_drive, pid in enumerate(order):
            d = _drive(
                ctx,
                built,
                specs32,
                pid,
                held=(n_drive > 0),
                wall_term=wall,
                z_s=z_s,
                partition=False,
                keep_e=False,
            )
            cols[pid] = d["s"]
            kinds.append(d["kind"])
            cad.append(d["t_solve"] + d["t_post"])
        s32 = np.array([[cols[pj][pi] for pj in order] for pi in order], dtype=np.complex128)
        say(
            f"32x32 (uncorrected widths, kappa not applied): {kinds.count('solve')} factorising, {kinds.count('held')} held; first "
            f"{cad[0]:.2f} s, held max {max(cad[1:]):.2f} s; t={time.perf_counter() - t_window:.1f} s"
        )
        assert kinds.count("solve") == 1 and kinds.count("held") == 31, kinds
    else:
        kappa = float(saved["kappa"])
        s32 = np.array(saved["s32"], dtype=np.complex128)
        saved_order = [str(p) for p in saved["order"]]
        cells_saved = int(saved["cells"])
        say(
            f"npz: cells {cells_saved}, kappa {kappa:.6e} (diagnostic only, not applied; "
            f"ASSERTED cells and order equal)"
        )
        assert n_cells == cells_saved, (n_cells, cells_saved)
        assert saved_order == order, saved_order

    # ---- (M) ---------------------------------------------------------------
    ratio32 = _reciprocity_ratio(s32)
    sigma32 = float(np.linalg.norm(s32, 2))
    eye32 = np.eye(32)
    say(
        f"(M) 32x32 reciprocity {ratio32:.6e} (band {RECIPROCITY_BAND:.0e}, ASSERTED; PREDICTED "
        f"≈{PRED_RECIP:.0e}); sigma_max {sigma32:.12f} (ceiling {MATRIX_PASSIVITY_CEILING}, "
        f"ASSERTED), 1-sigma_max {1.0 - sigma32:.3e} (PREDICTED {PRED_SMAX_MARGIN}); "
        f"cond(I-S32) {np.linalg.cond(eye32 - s32):.6e}, cond(I+S32) "
        f"{np.linalg.cond(eye32 + s32):.6e} (printed); S_driven(P17) {s32[0, 0]:.6f}; "
        f"beside H12b uncorrected sigma_max {H12B_UNCORRECTED_SIGMA_MAX} (ring order, printed)"
    )
    if os.path.exists(H8A_NPZ):
        h8a = np.load(H8A_NPZ, allow_pickle=False)
        if [str(p) for p in h8a["order"]] == order:
            s800 = np.array(h8a["s32"], dtype=np.complex128)
            say(
                f"||S_Cu - S_800||_F/||S_800||_F = "
                f"{np.linalg.norm(s32 - s800) / np.linalg.norm(s800):.6e} (printed)"
            )
        else:
            say("H8a npz order differs; S_Cu vs S_800 not printed")
    else:
        say(f"{H8A_NPZ} absent (not a failure)")
    assert ratio32 <= RECIPROCITY_BAND, ratio32
    assert sigma32 <= MATRIX_PASSIVITY_CEILING, sigma32

    # ---- circuit sweep and dense scan ---------------------------------------
    _g, _v, roots = tuning_sweep(s32, FREQUENCY_HZ, terminated=terminated)
    tuned = select_c_tuned(roots)
    say(
        f"tuning sweep: {len(roots)} sign change(s) of Im Z_in; copper select_c_tuned "
        + (f"{tuned['c_f'] * 1e12:.4f} pF" if tuned is not None else "NONE")
        + f" beside sigma=800 {SIGMA800_C_TUNED_PF} pF (PREDICTED {PRED_C_TUNED_PF} pF, printed)"
    )
    grid = np.geomspace(C_TUNED_F, C_MODE_F, SCAN_POINTS)
    conds, s11s, shp, shm = _circuit_scan(s32, FREQUENCY_HZ, grid, order, slots, z0)
    maxima, pts = _rule_points(grid, conds)
    say(
        f"scan: {SCAN_POINTS} log-spaced C in [C_tuned, C_mode]; {len(maxima)} interior local "
        f"maxima of cond; cond at C_tuned {conds[0]:.6e}"
    )
    if mode == "W1":
        extrema = set(maxima) | {
            i
            for i in range(1, SCAN_POINTS - 1)
            if conds[i] < conds[i - 1] and conds[i] <= conds[i + 1]
        }
        say("    i     C [pF]        cond(I-GS_tt)   |S11|      share ring+  share ring-")
        for i in range(SCAN_POINTS):
            if i % 20 == 0 or i in extrema or i == SCAN_POINTS - 1:
                tag = " MAX" if i in maxima else (" min" if i in extrema else "")
                say(
                    f"    {i:4d}  {grid[i] * 1e12:12.6f}  {conds[i]:14.6e}  {s11s[i]:.6f}  "
                    f"{shp[i]:.6f}     {shm[i]:.6f}{tag}"
                )
    say(
        "    all-C anti-mode capacitances (sigma=800, printed): "
        + ", ".join(f"k={k} {v} pF" for k, v in ANTI_MODE_C_PF.items())
    )
    cu_pts = {}
    if pts is not None:
        for key in ("on_lo", "on_hi", "off1", "off2"):
            i = pts[key]
            if i is None:
                say(f"    copper rule point {key}: NONE")
                continue
            cu_pts[key] = float(grid[i])
            say(
                f"    copper rule point {key}: C = {grid[i]:.12e} F ({grid[i] * 1e12:.4f} pF) "
                f"beside sigma=800 {sig800[key][0] * 1e12:.4f} pF; cond {conds[i]:.6e}"
            )
        for key, ref in zip(("on_lo", "on_hi"), SIGMA800_COND_MAX_PF):
            if key in cu_pts:
                say(
                    f"    cond max {key}: {cu_pts[key] * 1e12:.4f} pF vs {ref} pF rel "
                    f"{cu_pts[key] * 1e12 / ref - 1.0:+.4f} (PREDICTED ±{PRED_COND_MAX_RTOL})"
                )
    else:
        say(f"    fewer than 2 interior local maxima of cond ({maxima}): no copper rule points")

    if mode == "W1":
        meas_t, pred_t, r_t, cond_t, _ = pair_at(C_TUNED_F, "C_tuned", KAPPA_APPLIED, s32)
        _, pred_ctl = tuned_input(
            s32, FREQUENCY_HZ, NEGATIVE_CONTROL_FACTOR * C_TUNED_F, terminated=terminated
        )
        r_ctl = _residual(meas_t, pred_ctl)
        mult = r_ctl / REDUCTION_BAND
        say(
            f"negative control (PREDICTED {PRED_NEG}x, printed, not asserted): "
            f"{NEGATIVE_CONTROL_FACTOR} x C_tuned prediction vs MEAS at C_tuned residual "
            f"{r_ctl:.6e} = {mult:.1f}x REDUCTION_BAND (>= 10x: "
            f"{'held' if mult >= 10.0 else 'NOT held'}; sigma=800 {SIGMA800_NEG_MULTIPLE:.0f}x)"
        )
        meas_m, pred_m, r_m, cond_m, _ = pair_at(C_MODE_F, "C_mode", KAPPA_APPLIED, s32)
        say(
            f"W1 residuals (printed): C_tuned r_Cu {r_t:.6e} vs r_800 {sig800['tuned'][1]:.6e}; "
            f"C_mode r_Cu {r_m:.6e} vs r_800 {sig800['mode'][1]:.6e} (rho "
            f"{r_m / sig800['mode'][1]:.4f}; PREDICTED r_Cu(C_mode) in {PRED_R_MODE}); CAVEAT: "
            f"r_Cu at uncorrected widths, r_800 kappa-corrected"
        )
        if comm.rank == 0:
            os.makedirs(os.path.dirname(STEP2P_NPZ), exist_ok=True)
            np.savez(
                STEP2P_NPZ,
                s32=s32,
                kappa=kappa,
                order=np.array(order),
                cells=n_cells,
                meas_tuned=meas_t,
                r_tuned=r_t,
                cond_tuned=cond_t,
                meas_mode=meas_m,
                r_mode=r_m,
                cond_mode=cond_m,
            )
            say(f"saved {STEP2P_NPZ}")
        comm.Barrier()
    else:
        meas_t, pred_t, r_t, cond_t, _ = pair_at(
            C_TUNED_F, "(A3) guard C_tuned", KAPPA_APPLIED, s32
        )
        a3 = _entrywise_rel(meas_t, np.array(saved["meas_tuned"]))
        say(
            f"(A3) guard MEAS vs W1 npz entrywise {a3:.3e} (band {LEG_D0_REPRODUCTION_BAND:.0e}, "
            f"ASSERTED)"
        )
        assert a3 <= LEG_D0_REPRODUCTION_BAND, a3
        res = {
            "tuned": (C_TUNED_F, r_t, cond_t),
            "mode": (C_MODE_F, float(saved["r_mode"]), float(saved["cond_mode"])),
        }
        for key in ("on_lo", "off1", "on_hi", "off2"):
            _, _, r, cond, _ = pair_at(sig800[key][0], f"sigma=800 rule {key}", KAPPA_APPLIED, s32)
            res[key] = (sig800[key][0], r, cond)
        guard_s = pair_times[0]
        extra = [
            k
            for k in ("on_lo", "on_hi")
            if k in cu_pts and abs(cu_pts[k] / sig800[k][0] - 1.0) > ON_POINT_C_RTOL
        ]
        projection = t_build + (1 + 4 + len(extra)) * guard_s + 50.0
        cu_own = {}
        if extra and projection <= W2_PRICE_CEILING_S:
            for k in extra:
                _, _, r, cond, _ = pair_at(cu_pts[k], f"copper own {k}", KAPPA_APPLIED, s32)
                cu_own[k] = (cu_pts[k], r, cond)
        elif extra:
            say(f"copper own {extra}: deferred, projection {projection:.0f} s")
        else:
            say("copper on_lo / on_hi within 0.5 % of sigma=800's (or absent): no extra point")

        # ---- reading table and rule (printed only) --------------------------
        keys = ["tuned", "on_lo", "off1", "on_hi", "off2", "mode"]
        say(
            "CAVEAT: r_Cu is copper at UNCORRECTED widths (kappa not applied, ruling "
            "2026-10-02); r_800 / cond_800 are sigma=800 KAPPA-CORRECTED (step1c-W2), so rho "
            "compares across that difference; rule printed as pre-registered, not refit"
        )
        say(
            "reading table: key  C[pF]  r_Cu  r_800  rho=r_Cu/r_800  cond_Cu  r_Cu/cond_Cu  cond_800"
        )
        rho = {}
        for k in keys:
            c, r, cd = res[k]
            rho[k] = r / sig800[k][1]
            say(
                f"    {k:6s} {c * 1e12:10.4f}  {r:.6e}  {sig800[k][1]:.6e}  {rho[k]:.4f}  "
                f"{cd:.6e}  {r / cd:.6e}  {sig800[k][2]:.6e}"
            )
        for k, (c, r, cd) in cu_own.items():
            say(f"    copper own {k}: {c * 1e12:.4f} pF  r_Cu {r:.6e}  cond {cd:.6e}")
        cs = np.array([res[k][0] for k in keys])
        rs = np.array([res[k][1] for k in keys])
        slope_cu = float(np.polyfit(np.log(cs), np.log(rs), 1)[0])
        r_mode = res["mode"][1]
        big = [k for k in keys if res[k][0] >= RULE_CONDUCTOR_C_MIN_F]
        conductor = r_mode <= RULE_CONDUCTOR_R_MODE and all(
            rho[k] <= RULE_CONDUCTOR_RHO for k in big
        )
        tracks = (
            all(RULE_TRACKS_RHO[0] <= rho[k] <= RULE_TRACKS_RHO[1] for k in keys)
            and abs(slope_cu - RULE_TRACKS_SLOPE[0]) <= RULE_TRACKS_SLOPE[1]
        )
        verdict = "conductor" if conductor else ("tracks" if tracks else "mixed")
        say(
            f"log-log slope r_Cu vs C (6 points) {slope_cu:+.4f} beside sigma=800 "
            f"{SIGMA800_SLOPE:+.4f} (PREDICTED {PRED_SLOPE[0]} ± {PRED_SLOPE[1]})"
        )
        say(
            f"rho range {min(rho.values()):.4f} .. {max(rho.values()):.4f} (PREDICTED "
            f"{PRED_RHO}); r_Cu(C_mode) {r_mode:.6e} (PREDICTED {PRED_R_MODE}); C >= 75 pF "
            f"points {big}"
        )
        say(f"reading rule (pre-registered, printed only): classification = {verdict.upper()}")
        say(
            f"H8b-configuration estimate {H8B_TRANSFER} x r_Cu(C_mode) = "
            f"{H8B_TRANSFER * r_mode:.6e} — assumed transfer, not measured (H8b sigma=800 "
            f"{H8B_RESIDUAL:.6e}, not the comparand)"
        )

    for label, pred, meas in checks:
        rp, rm = _reciprocity_ratio(pred), _reciprocity_ratio(meas)
        cp = float(np.max(np.sum(np.abs(pred) ** 2, axis=0)))
        cm = float(np.max(np.sum(np.abs(meas) ** 2, axis=0)))
        say(
            f"(C) {label}: reciprocity circuit {rp:.3e} in-model {rm:.3e} (band "
            f"{RECIPROCITY_BAND:.0e}, ASSERTED); col sum|S|^2 circuit {cp:.9f} in-model "
            f"{cm:.9f} (ceiling {COLUMN_PASSIVITY_CEILING}, ASSERTED)"
        )
    _m, rss_per, rss_sum = _rss_gib(comm)
    say(
        f"price: pairs {', '.join(f'{t:.1f}' for t in pair_times)} s; window "
        f"{time.perf_counter() - t_window:.1f} s at -n {comm.size}; RSS per rank "
        f"{', '.join(f'{r:.2f}' for r in (rss_per or []))} GiB, summed {rss_sum:.2f} GiB "
        f"(PREDICTED <= {PRED_RSS_GIB} GiB)"
    )
    for label, pred, meas in checks:
        assert _reciprocity_ratio(pred) <= RECIPROCITY_BAND, label
        assert _reciprocity_ratio(meas) <= RECIPROCITY_BAND, label
        assert float(np.max(np.sum(np.abs(pred) ** 2, axis=0))) <= COLUMN_PASSIVITY_CEILING, label
        assert float(np.max(np.sum(np.abs(meas) ** 2, axis=0))) <= COLUMN_PASSIVITY_CEILING, label
