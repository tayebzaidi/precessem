"""`WF-12` step 2 (§10 chain H12b, §9 item 169) — the copper 32×32 on the
F-human hole mesh, with one PEC column.

Selected by environment only (``FEM_EM_WF12_STEP2=1``), never by ``-k``.
Knobs: ``FEM_EM_WF12_STEP2_DRIVES=N`` drives P17 onward for N drives (unset =
32; below 32 the matrix gates print SKIPPED); ``FEM_EM_WF12_STEP2_PEC=0`` skips
the PEC column.

Fixture: ``_build_f_human(comm, as_hole=True)`` unchanged (H12a).  Every Cu
drive terminates all 32 ring sheets at 50 Ω with `TH-14`'s Leontovich copper
wall (σ = 5.8e7) on the cavity, first in the factorising call's bilinear list
(H12a's order); later drives back-substitute through the held factor.

Anchors (asserted, every band imported):
  (A0)  hole cells = 369 387 (`20261002T020713Z_WF-12-step1-W2.log:6352`)
        within ``CELL_COUNT_BAND``; (iv) under reduction as in step 1;
  (R)   Cu P17 reproduces P_src 6.916352890e-03 W and P_surf 8.551752204e-07 W
        (`…step1-W2.log:6354`) within ``REUSE_REPRODUCTION_RTOL``;
  (P)   every Cu drive: P_src = ΣP_sheet + P_phantom + P_(Ω∖phantom) + P_surf
        within ``DISCRETE_IDENTITY_RTOL``; exactly one factorising Cu solve;
  (ii)  ‖S − Sᵀ‖_F/‖S‖_F ≤ ``RECIPROCITY_BAND``;
  (iii) σ_max and every column Σ|S|² ≤ ``MATRIX_PASSIVITY_CEILING``;
  (iv′) ``N_SYMMETRY_CLASSES`` classes, azimuth residual < ``AZIMUTH_MATCH_DEG``,
        worst class spread ≤ ``OPPOSITE_SPREAD_BAND``;
  (PEC) the PEC P17 identity closes without P_surf within the identity band,
        and P_surf evaluated on the PEC field ≤ 1e-12 × the Cu value.
Negative control (asserted, rule (e)): on every Cu drive the identity with
P_surf dropped misses by ≥ 10× ``DISCRETE_IDENTITY_RTOL`` (step 1 measured
124× on P17, `…step1-W2.log:6355`).
Predicted, printed only: Cu ≠ PEC on P17; the 1.01 column-scale lift; the
coil-loss share P_surf/(P_src − ΣP_sheet) per drive; spreads; S_driven(P17).

Scope: degree 1, 64 MHz, 50 Ω, this one fixture.  Not a PEC 32×32, not
|B1+| / hotspot (H9), not C_mode (H12c), not loss against AED, not degree 2.
"""

from __future__ import annotations

import os
import time

import numpy as np
import pytest
from mpi4py import MPI

from tests.complex_mode import complex_only

STEP2_ENV = "FEM_EM_WF12_STEP2"
DRIVES_ENV = "FEM_EM_WF12_STEP2_DRIVES"
PEC_ENV = "FEM_EM_WF12_STEP2_PEC"
DRIVEN = "P17"
COPPER_SIGMA = 5.8e7

# (A0) and (R) records, H12a W2 (`20261002T020713Z_WF-12-step1-W2.log:6352,6354`).
HOLE_CELL_RECORD = 369_387
P17_P_SRC_RECORD_W = 6.916352890e-03
P17_P_SURF_RECORD_W = 8.551752204e-07
# (PEC): n×E ≡ 0 on pinned facets, so P_surf on the PEC field is round-off.
PEC_SURF_RATIO_CEILING = 1.0e-12
# Negative control bar (rule (e), asserted): ceiling is P_surf/P_src ≈ 124×.
BLIND_FACTOR = 10.0

# Predictions written first (§9 item 169), printed only.
PREDICTED_COIL_SHARE = 0.370
SIGMA800_COIL_SHARE = 0.9945
PREDICTED_RECIPROCITY = 1e-14
PREDICTED_SPREAD = (0.025, 0.045)
PREDICTED_COLSUM = (0.985, 0.9991)
PREDICTED_SMAX_MARGIN = (1e-8, 1e-4)
PREDICTED_DSUM = (1.5e-4, 7e-4)  # central 3.4e-4
PREDICTED_MAXDS = (1e-5, 1e-2)
PREDICTED_LIFT = 1e6
S_DRIVEN_STEP0_RECORD = "0.407423+0.344417j"
W2_PRICE_CEILING_S = 540.0


def _knobs():
    raw = os.environ.get(DRIVES_ENV, "").strip()
    n = int(raw) if raw else 32
    if not 1 <= n <= 32:
        raise ValueError(f"{DRIVES_ENV} must be in 1..32, got {raw!r}")
    pec = os.environ.get(PEC_ENV, "").strip() != "0"
    return n, pec


def _run_step2(n_drives=32, run_pec=True, *, field_hook=None):
    """Body of the step-2 gate, lifted to module level (rule (a), `EX-80`, §9
    item 179) so `examples/mri/06_*.py` can call it.  ``field_hook`` (default
    None: every print and assert is unchanged) is called inside ``_drive`` as
    ``field_hook(kind, driven_id, e, sigma_f, z_s)`` before the field is freed."""
    import ufl
    from dolfinx import fem

    from fem_em_solver.core import TimeHarmonicSolver
    from fem_em_solver.core.cavity import surface_resistance_ohm
    from fem_em_solver.ports.lumped import (
        lumped_port_bilinear_term,
        lumped_port_linear_term,
        sheet_terminal_current,
    )
    from tests.mesh.test_birdcage_port_sheet_prerequisite import CELL_COUNT_BAND
    from tests.validation.test_birdcage_f_human_rung import _ring_ports
    from tests.validation.test_birdcage_power_identity import (
        DISCRETE_IDENTITY_RTOL,
        _sheet_field_dissipation_w,
        _source_power_w,
    )
    from tests.validation.test_port_birdcage_four_port import TERMINATED_PORT_IMPEDANCE_OHM
    from tests.validation.test_port_birdcage_lumped_column import PHANTOM_CELL_TAG
    from tests.validation.test_port_birdcage_ring_column import (
        AZIMUTH_MATCH_DEG,
        CONTROL_COLUMN_SCALE,
        OPPOSITE_SPREAD_BAND,
        _reciprocity_ratio,
    )
    from tests.validation.test_port_birdcage_ring_matrix import (
        MATRIX_PASSIVITY_CEILING,
        N_AZIMUTH_CLASSES,
        N_SYMMETRY_CLASSES,
        REUSE_REPRODUCTION_RTOL,
        _azimuth_class,
    )
    from tests.validation.test_port_lumped_sheet_sweep import RECIPROCITY_BAND
    from tests.validation.test_th14_birdcage_copper import (
        _leontovich_term,
        _problem,
        _surface_loss_w,
    )
    from tests.validation.test_wf7_f_human_b1_quadrature import DEGREE, _build_f_human
    from tests.validation.test_wf12_f_human_hole import (
        _census_checks,
        _owned_census,
        _rss_gib,
    )

    comm = MPI.COMM_WORLD
    t_window = time.perf_counter()

    def say(msg):
        if comm.rank == 0:
            print(f"[WF-12 step2] {msg}", flush=True)

    say(f"start: -n {comm.size}, degree {DEGREE}, drives {n_drives} from {DRIVEN}, "
        f"PEC column {'on' if run_pec else 'off'}")

    # ---- fixture, (A0), (iv) ---------------------------------------------
    t0 = time.perf_counter()
    built = _build_f_human(comm, as_hole=True)
    t_build = time.perf_counter() - t0
    ctx = built["ctx"]
    msh, tags_f, omega, cell_tags = ctx["msh"], ctx["tags_f"], ctx["omega"], ctx["cell_tags"]
    census = _owned_census(msh, cell_tags, comm)
    census_fails = _census_checks(census, built["cells"], _ring_ports())
    oc = built["outer_census"]
    cell_rel = built["cells"] / HOLE_CELL_RECORD - 1.0
    say(f"(A0) hole cells {built['cells']} vs record {HOLE_CELL_RECORD}, rel {cell_rel:.3e} "
        f"(band {CELL_COUNT_BAND}, ASSERTED); build {t_build:.1f} s; (iv) failures "
        f"{census_fails or 'none'}; outer census {oc}")
    assert abs(cell_rel) <= CELL_COUNT_BAND, cell_rel
    assert not census_fails, census_fails
    assert oc["tag_collision"] == 0 and oc["cavity"] > 0, oc
    assert oc["cavity_exterior"] == oc["cavity"] and oc["outer"] + oc["cavity"] == oc["exterior"], oc

    z_s = (1.0 + 1.0j) * surface_resistance_ohm(omega, COPPER_SIGMA)
    wall = _leontovich_term(msh, built["facets"], z_s, omega)
    z_p = complex(TERMINATED_PORT_IMPEDANCE_OHM)
    dx = ufl.Measure("dx", domain=msh, subdomain_data=cell_tags)
    acct = {"mesh": msh, "facet_tags": tags_f}

    def _drive(solver, driven_id, *, held, wall_term):
        """Local copy of `_solve_one_drive` (`test_port_birdcage_ring_column.py:
        223–293`) with the wall term first in the factorising call, plus H12a's
        (P) block (`test_wf12_f_human_hole.py:464–481`)."""
        sheets = [sp.sheet(driven=(sp.port_id == driven_id)) for sp in ctx["specs"]]
        dsheet = next(s for s in sheets if s.port_id == driven_id)
        linear = [lambda test, _s=dsheet: lumped_port_linear_term(
            msh, tags_f, _s, test, omega_rad_per_s=omega)]
        comm.Barrier()
        t0 = time.perf_counter()
        if held:
            fields = solver.solve_with_held_factorization(extra_linear_terms=linear)
        else:
            bil = ([wall_term] if wall_term is not None else []) + [
                (lambda trial, test, _s=s: lumped_port_bilinear_term(
                    msh, tags_f, _s, trial, test, omega_rad_per_s=omega))
                for s in sheets
            ]
            fields = solver.solve(current_density=None, project_source=False,
                                  extra_bilinear_terms=bil, extra_linear_terms=linear)
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
        sigma_f = fields.sigma_field

        def _vol(measure):
            form = fem.form(0.5 * sigma_f * ufl.inner(e, e) * measure)
            return float(np.real(comm.allreduce(fem.assemble_scalar(form), op=MPI.SUM)))

        p_sheets = float(sum(_sheet_field_dissipation_w(acct, s, e, omega) for s in sheets))
        p_src = float(_source_power_w(acct, dsheet, e, omega))
        p_total = _vol(dx)
        p_phantom = _vol(dx(PHANTOM_CELL_TAG))
        p_surf = _surface_loss_w(msh, built["facets"], z_s, e, comm)
        if field_hook is not None:
            field_hook("cu" if wall_term is not None else "pec", driven_id, e, sigma_f,
                       z_s if wall_term is not None else None)
        del fields
        t_post = time.perf_counter() - t1
        return {
            "s": s_col, "kind": "held" if held else "solve", "t_solve": t_solve,
            "t_post": t_post, "p_src": p_src, "p_sheets": p_sheets,
            "p_phantom": p_phantom, "p_rest": p_total - p_phantom, "p_surf": p_surf,
        }

    # ---- Cu drives -------------------------------------------------------
    ring_ports = _ring_ports()
    all_ids = [f"P{i}" for i in ring_ports]
    assert all_ids.index(DRIVEN) == 0, all_ids[:3]
    driven_ids = all_ids[:n_drives]
    cols = {}
    for k, pid in enumerate(driven_ids):
        d = _drive(ctx["solver"], pid, held=(k > 0), wall_term=wall)
        p_vol = d["p_phantom"] + d["p_rest"]
        d["rel"] = abs(d["p_src"] - d["p_sheets"] - p_vol - d["p_surf"]) / abs(d["p_src"])
        d["blind"] = abs(d["p_src"] - d["p_sheets"] - p_vol) / abs(d["p_src"])
        d["share"] = d["p_surf"] / (d["p_src"] - d["p_sheets"])
        cols[pid] = d
        say(f"Cu drive {k + 1:2d}/{n_drives} {pid}: {d['kind']} {d['t_solve']:.2f} s + post "
            f"{d['t_post']:.2f} s; S_driven {d['s'][pid]:.6f}; sum|S|^2 "
            f"{sum(abs(v) ** 2 for v in d['s'].values()):.6f} | P_src {d['p_src']:.9e} W "
            f"sumP_sheet {d['p_sheets']:.9e} P_phantom {d['p_phantom']:.9e} P_rest "
            f"{d['p_rest']:.9e} P_surf {d['p_surf']:.9e} W; (P) rel {d['rel']:.3e}; NEG blind "
            f"{d['blind']:.4e} = {d['blind'] / DISCRETE_IDENTITY_RTOL:.1f}x band; coil share "
            f"P_surf/(P_src-sumP_sheet) {d['share']:.4f} (PREDICTED ≈{PREDICTED_COIL_SHARE}; "
            f"sigma=800 {SIGMA800_COIL_SHARE}) t={time.perf_counter() - t_window:.1f} s")
    kinds = [cols[p]["kind"] for p in driven_ids]
    rss_mine, rss_per, rss_sum = _rss_gib(comm)
    held_cad = [cols[p]["t_solve"] + cols[p]["t_post"] for p in driven_ids[1:]]
    first = cols[DRIVEN]
    t_first = first["t_solve"] + first["t_post"]
    max_cad = max(held_cad) if held_cad else float("nan")
    say(f"solve kinds: {kinds.count('solve')} factorising, {kinds.count('held')} held; "
        f"held cadence (solve + post) max {max_cad:.2f} s, mean "
        f"{(np.mean(held_cad) if held_cad else float('nan')):.2f} s (proxy 3.75 s at sigma=800); "
        f"RSS per rank {', '.join(f'{r:.2f}' for r in (rss_per or []))} GiB, summed {rss_sum:.2f} GiB")
    if held_cad:
        projection = t_build + t_first + 31 * max_cad + t_first + 50.0
        projection_nopec = t_build + t_first + 31 * max_cad + 50.0
        say(f"PROJECTED W2: build {t_build:.1f} + Cu factor {t_first:.1f} + 31 x {max_cad:.2f} "
            f"+ PEC {t_first:.1f} + 50 = {projection:.0f} s (Cu-only {projection_nopec:.0f} s; "
            f"ceiling {W2_PRICE_CEILING_S:.0f} s; predicted W2 320-420 s)")

    # (R), (P), NEG
    r_src = abs(first["p_src"] / P17_P_SRC_RECORD_W - 1.0)
    r_surf = abs(first["p_surf"] / P17_P_SURF_RECORD_W - 1.0)
    say(f"(R) P17 P_src {first['p_src']:.9e} vs {P17_P_SRC_RECORD_W:.9e} rel {r_src:.3e}; "
        f"P_surf {first['p_surf']:.9e} vs {P17_P_SURF_RECORD_W:.9e} rel {r_surf:.3e} "
        f"(band {REUSE_REPRODUCTION_RTOL:g}, ASSERTED)")
    worst_rel = max(cols[p]["rel"] for p in driven_ids)
    min_blind = min(cols[p]["blind"] for p in driven_ids)
    shares = [cols[p]["share"] for p in driven_ids]
    say(f"(P) worst rel over {n_drives} Cu drives {worst_rel:.3e} (band "
        f"{DISCRETE_IDENTITY_RTOL:g}, ASSERTED); NEG min blind miss {min_blind:.4e} = "
        f"{min_blind / DISCRETE_IDENTITY_RTOL:.1f}x the band (bar {BLIND_FACTOR:.0f}x, ASSERTED; "
        f"ceiling ≈124x); coil share min {min(shares):.4f} max {max(shares):.4f}")
    assert kinds.count("solve") == 1 and kinds.count("held") == n_drives - 1, kinds
    assert r_src <= REUSE_REPRODUCTION_RTOL, r_src
    assert r_surf <= REUSE_REPRODUCTION_RTOL, r_surf
    assert worst_rel <= DISCRETE_IDENTITY_RTOL, worst_rel
    assert min_blind >= BLIND_FACTOR * DISCRETE_IDENTITY_RTOL, min_blind

    # ---- matrix gates ----------------------------------------------------
    matrix = None
    pec = None
    if n_drives == 32:
        s = np.array([[cols[pj]["s"][pi] for pj in driven_ids] for pi in driven_ids],
                     dtype=complex)
        ratio = _reciprocity_ratio(s)
        control = s.copy()
        control[:, 0] *= CONTROL_COLUMN_SCALE
        control_ratio = _reciprocity_ratio(control)
        sigma_max = float(np.linalg.norm(s, 2))
        col_norms = np.sum(np.abs(s) ** 2, axis=0)
        # Class table lifted from `test_wf7_f_human_ring_matrix.py:412–461`,
        # on `built["sheets"]`.
        by_ord = {sh["ordinal"]: sh for sh in built["sheets"]}
        az = np.array([by_ord[i]["azimuth_deg"] for i in ring_ports])
        zc = np.array([by_ord[i]["z"] for i in ring_ports])
        n = len(ring_ports)
        classes = {}
        worst_residual_deg = 0.0
        for i in range(n):
            for j in range(n):
                steps, residual_deg = _azimuth_class(az[i], az[j])
                worst_residual_deg = max(worst_residual_deg, residual_deg)
                classes.setdefault((bool(zc[i] * zc[j] > 0.0), steps), []).append(
                    float(abs(s[i, j])))
        spreads = {}
        for key, mags in classes.items():
            m = np.array(mags)
            spreads[key] = (float((m.max() - m.min()) / m.mean()), len(m), float(m.mean()))
        worst = max(spreads, key=lambda key: spreads[key][0])
        say(f"(ii) reciprocity {ratio:.6e} (band {RECIPROCITY_BAND:.0e}, ASSERTED; PREDICTED "
            f"≈{PREDICTED_RECIPROCITY:.0e}, sigma=800 1.29e-14)")
        say(f"(iii) sigma_max {sigma_max:.12f} (ceiling {MATRIX_PASSIVITY_CEILING}, ASSERTED), "
            f"1-sigma_max {MATRIX_PASSIVITY_CEILING - sigma_max:+.3e} (PREDICTED "
            f"{PREDICTED_SMAX_MARGIN}); column sum|S|^2 min {col_norms.min():.6f} max "
            f"{col_norms.max():.6f} (ASSERTED <= 1; PREDICTED {PREDICTED_COLSUM})")
        say(f"(iv') classes {len(classes)} (expected {N_SYMMETRY_CLASSES}, ASSERTED); worst "
            f"azimuth residual {worst_residual_deg:.3e} deg (< {AZIMUTH_MATCH_DEG:.0e}, ASSERTED)")
        for same_ring in (True, False):
            for steps in range(N_AZIMUTH_CLASSES):
                if (same_ring, steps) not in spreads:
                    say(f"     {'same ' if same_ring else 'other'} ring, {steps} steps: EMPTY")
                    continue
                sp, cnt, mean = spreads[(same_ring, steps)]
                say(f"     {'same ' if same_ring else 'other'} ring, {steps} steps  n = {cnt:4d}"
                    f"  mean |S| = {mean:.9e}  spread {sp * 100:8.4f}%"
                    + ("   <-- worst" if (same_ring, steps) == worst else ""))
        say(f"     worst class spread {spreads[worst][0]:.4e} (band {OPPOSITE_SPREAD_BAND}, "
            f"ASSERTED; PREDICTED {PREDICTED_SPREAD}, sigma=800 3.70%)")
        say(f"S_driven(P17) {cols[DRIVEN]['s'][DRIVEN]:.6f} beside sigma=800 record "
            f"{S_DRIVEN_STEP0_RECORD} (printed)")
        say(f"control (predicted, printed): column {driven_ids[0]} x {CONTROL_COLUMN_SCALE}: "
            f"ratio {control_ratio:.6e} = {control_ratio / max(ratio, 1e-300):.3e}x the measured "
            f"(PREDICTED > {PREDICTED_LIFT:.0e}x; sigma=800 1.6e11x)")
        assert ratio <= RECIPROCITY_BAND, ratio
        assert sigma_max <= MATRIX_PASSIVITY_CEILING, sigma_max
        assert float(col_norms.max()) <= MATRIX_PASSIVITY_CEILING, col_norms.max()
        assert len(classes) == N_SYMMETRY_CLASSES, len(classes)
        assert worst_residual_deg < AZIMUTH_MATCH_DEG, worst_residual_deg
        assert spreads[worst][0] <= OPPOSITE_SPREAD_BAND, spreads[worst]
        matrix = {"s": s, "ratio": ratio, "control_ratio": control_ratio,
                  "sigma_max": sigma_max, "col_norms": col_norms, "n_classes": len(classes),
                  "worst_residual_deg": worst_residual_deg, "spreads": spreads, "worst": worst}
    else:
        say(f"(ii), (iii), (iv') and the column-scale control: SKIPPED "
            f"({DRIVES_ENV}={n_drives} < 32)")

    # ---- PEC column ------------------------------------------------------
    if run_pec:
        del ctx["solver"]
        t0 = time.perf_counter()
        pec_solver = TimeHarmonicSolver(
            _problem(msh, cell_tags, 64e6, built["facets"], None), degree=DEGREE)
        d = _drive(pec_solver, DRIVEN, held=False, wall_term=None)
        t_pec = time.perf_counter() - t0
        p_vol = d["p_phantom"] + d["p_rest"]
        rel_pec = abs(d["p_src"] - d["p_sheets"] - p_vol) / abs(d["p_src"])
        surf_ratio = d["p_surf"] / first["p_surf"]
        sum_cu = sum(abs(v) ** 2 for v in first["s"].values())
        sum_pec = sum(abs(v) ** 2 for v in d["s"].values())
        max_ds = max(abs(d["s"][p] - first["s"][p]) for p in first["s"])
        say(f"(PEC) P17 fresh solver {t_pec:.1f} s (solve {d['t_solve']:.1f} s): P_src "
            f"{d['p_src']:.9e} sumP_sheet {d['p_sheets']:.9e} P_phantom {d['p_phantom']:.9e} "
            f"P_rest {d['p_rest']:.9e} W; identity w/o P_surf rel {rel_pec:.3e} (band "
            f"{DISCRETE_IDENTITY_RTOL:g}, ASSERTED); P_surf(PEC field) {d['p_surf']:.3e} W = "
            f"{surf_ratio:.3e} x Cu (ceiling {PEC_SURF_RATIO_CEILING:g}, ASSERTED)")
        say(f"Cu != PEC (predicted, printed): sum|S_i,17|^2 PEC {sum_pec:.9f} - Cu {sum_cu:.9f} "
            f"= {sum_pec - sum_cu:.4e} (PREDICTED ≈3.4e-4 in {PREDICTED_DSUM}); max_i|dS| "
            f"{max_ds:.4e} (PREDICTED in {PREDICTED_MAXDS}); S_driven PEC {d['s'][DRIVEN]:.6f}")
        assert rel_pec <= DISCRETE_IDENTITY_RTOL, rel_pec
        assert surf_ratio <= PEC_SURF_RATIO_CEILING, surf_ratio
        pec = {"d": d, "rel_pec": rel_pec, "surf_ratio": surf_ratio,
               "dsum": sum_pec - sum_cu, "max_ds": max_ds}
    else:
        say(f"(PEC) SKIPPED ({PEC_ENV}=0)")
    _m, rss_per, rss_sum = _rss_gib(comm)
    say(f"window {time.perf_counter() - t_window:.1f} s; RSS summed {rss_sum:.2f} GiB")
    return {
        "cells": built["cells"], "t_build": t_build, "cols": cols, "kinds": kinds,
        "r_src": r_src, "r_surf": r_surf, "worst_rel": worst_rel, "min_blind": min_blind,
        "shares": shares, "matrix": matrix, "pec": pec, "msh": msh, "cell_tags": cell_tags,
        "facets": built["facets"], "specs": ctx["specs"],
    }


@complex_only
def test_wf12_step2_copper_matrix():
    if os.environ.get(STEP2_ENV, "").strip() != "1":
        pytest.skip(f"{STEP2_ENV} unset: WF-12 step 2 is selected by environment only")
    _run_step2(*_knobs())
