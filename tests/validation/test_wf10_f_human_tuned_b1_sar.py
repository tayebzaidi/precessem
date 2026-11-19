"""`WF-10` step 2 part 1 (§10 chain H9′, §9 item 358) — the tuned two-port
co-rotating ``|B₁⁺|`` on the copper F-human hole mesh at C_in,F-human with
r = 0.2 Ω on all 32 gaps, gated by the feed-bisector mirror identity and the
power identities.

Brief: ``docs/planning/chunks/WF-10.md`` "§9 item 358" and the 2026-10-10 10:30
review rulings beneath item 359 (which outrank the brief): the drive sense comes
from the convention code (``quadrature_phase_weights``, ``ccw`` = (1, −1j) at
slots 0 / 4); the node stays C_in,F-human = 118.2437 pF; the asserted set is the
feed-bisector mirror at ``C4_COVARIANCE_BAND`` plus (P_RC), (P_w), (Lin), (Σ_P);
C16 is printed, never asserted.

Two tests:

* ``test_wf10_step2_sense_coding`` (S0) — numpy only, always runs.
* ``test_wf10_step2_tuned_b1_quadrature`` — selected by ``FEM_EM_WF10_STEP2=1``
  (never ``-k``); heavy, ``-n 8``.

``build_tuned_quadrature`` is module-level so item 359 imports it; its returned
keys are fixed here.

Step 4a (§9 item 368, §10 H20a): ``FEM_EM_WF10_DEGREE`` — a degree knob on
``build_tuned_quadrature``.  Unset or empty ⇒ the degree-1 path, untouched;
``"2"`` ⇒ the degree-1 solver is deleted and a ``TimeHarmonicSolver(problem,
degree=2)`` built on the same problem (the `WF-12` step 4 precedent); anything
else ⇒ ``ValueError``.  (Deg) is asserted inside ``_apply_degree``.  At degree
2 the builder writes ``_node_json_path(2)`` (``…-d2.json``) and the step 2 test
runs ``with_z50=False`` ((R) / (Z50) are degree-1 comparands, printed only).
Two more tests: ``test_wf10_step4a_degree_knob_coding`` (K0, numpy only) and
``test_wf10_step4a_build_only`` (K1, ``FEM_EM_WF10_STEP4A=1``, build-only).

Scope: degree 1, 64 MHz, κ off, one copper hole fixture, one phantom; not f₀ /
Q / match / S₁₁ / tuning accuracy, not homogeneity, not absolute ``|B₁⁺|`` or
SAR, not C16 invariance, not compliance, not degree 2, not AED, not the 10 g
hotspot (item 359).
"""

from __future__ import annotations

import dataclasses
import json
import os
import time

import numpy as np
import pytest
from mpi4py import MPI

from tests.complex_mode import complex_only
from tests.validation.test_birdcage_b1_quadrature import quadrature_phase_weights
from tests.validation.test_wf10_f_human_copper_esr_root import (
    ESR_OHM,
    JSON_DIR,
    NPZ_PATH,
    SENSES,
    Y0_CODING_BAND,
    _active_y_in,
    _current_route_s,
    _source_voltage,
    _thevenin_feed_impedance,
    _z_rc,
)

STEP2_ENV = "FEM_EM_WF10_STEP2"
DRIVEN = "P17"
#: C_in,F-human (item 354, ruling (2); 10:30 ruling (1)) — the one constant
#: item 359 and `EX-75` import.  Never retyped.
C_IN_F_HUMAN_PF = 118.2437
#: The full-precision secant root (`20261010T094215Z_WF-10-step1e2-W2.log:76`),
#: printed beside the constant only.
C_IN_SECANT_PF = 118.243680
#: The co-rotating (1, −1j) P17 root (`…W2.log:72`), printed only.
C_IN_CCW_ROOT_PF = 118.241921
#: Step 4c (§9 item 511, §10 H28): the degree-2 P17 (1, +1j) root to 4
#: decimals (`20261017T070008Z_WF-10-step1e-H26.log:9000`), the same precedent
#: as ``C_IN_F_HUMAN_PF``.  Selected only by ``FEM_EM_WF10_NODE=d2root`` with
#: ``FEM_EM_WF10_DEGREE=2``.
C_IN_F_HUMAN_D2_PF = 115.5209
#: The full-precision degree-2 secant root, printed only.
C_IN_SECANT_D2_PF = 115.520946
NODE_ENV = "FEM_EM_WF10_NODE"
#: The build-only selector of step 4c's (K1n); never set in `xl` entry 41.
STEP4C_ENV = "FEM_EM_WF10_STEP4C_BUILD"
# Predictions written first (brief), printed only.
PREDICTED_P_W = 1.0e-10
PREDICTED_LIN = 1.0e-13
PREDICTED_SIGMA_P = 1.0e-12
PREDICTED_SENSE_RATIO = 10.0
PREDICTED_MIRROR = (0.2e-2, 1.5e-2)
PREDICTED_MISPAIRED = 0.90
PREDICTED_BLIND_PW = 7.0e-3
PREDICTED_DZ = 0.233
PREDICTED_C16 = 3.0e-2
PREDICTED_LAG_SHARE = 0.9
H17_K1_SHARE = 0.9940
PREDICTED_Y_IN = {
    ("P17", "(1,+1j)"): "Im ≈ +2.7e-6 S (|·| <= 1e-5)",
    ("P17", "(1,-1j)"): "Im ≈ +9e-6 S",
    ("P21", "(1,+1j)"): "Im ≈ -1.9e-4 S",
    ("P21", "(1,-1j)"): "Im ≈ -2.2e-4 S",
}
NODE_JSON = os.path.join(JSON_DIR, f"h9p-node-{C_IN_F_HUMAN_PF:.4f}.json")

# ---- step 4a (§9 item 368, §10 H20a): the degree knob ----------------------------
DEGREE_ENV = "FEM_EM_WF10_DEGREE"
STEP4A_ENV = "FEM_EM_WF10_STEP4A"
#: Hole-mesh topology records (`20261008T093243Z_WF-12-step4a-W1.log:6221–6222`,
#: same mesh): E, F and the degree-2 count 2E + 2F; printed beside, never a band.
E_RECORD = 500_568
F_RECORD = 784_444
N_DOFS_D2_RECORD = 2_570_024
#: `WF-12` 4a's summed RSS at build (`…093243Z…:6223`), printed beside.
RSS_D2_BUILD_RECORD_GIB = 1.13
PREDICTED_RSS_D2_GIB = 2.0
#: 358's W1 readings (`20261010T170356Z_WF-10-step2-W1.log:6296–6345`), printed
#: beside the step 2 re-gate (W2).
M_B_RECORD = 1.9647e-03
MEAN_B1P_RECORD_T = 2.884091e-07
#: At degree 2: which of step 2's gates are degree-1 record comparands (printed,
#: not run) and which stay asserted (the review's ruling 2 set).
DEGREE2_RECORD_GATES = {
    "printed": ("(R)", "(Z50)", "(S_rep)"),
    "asserted": (
        "(A0)",
        "(iv)",
        "(Deg)",
        "(P_RC)",
        "NEG node",
        "(P_w)",
        "(Lin)",
        "(Sigma_P)",
        "(Sense)",
        "(M_B)",
    ),
}


def _degree_knob():
    """``FEM_EM_WF10_DEGREE`` read at call time: unset / empty ⇒ ``None`` (the
    degree-1 path), ``"2"`` ⇒ 2, anything else (``"1"`` included) ⇒ ValueError."""
    raw = os.environ.get(DEGREE_ENV, "").strip()
    if raw == "":
        return None
    if raw == "2":
        return 2
    raise ValueError(f"{DEGREE_ENV}={raw!r}: only unset/empty (degree 1) or '2' are accepted")


def _node_pf():
    """Step 4c: ``FEM_EM_WF10_NODE`` read at call time.  Unset / empty ⇒
    ``C_IN_F_HUMAN_PF`` (the same object); ``"d2root"`` ⇒ ``C_IN_F_HUMAN_D2_PF``,
    only with ``FEM_EM_WF10_DEGREE=2`` (the degree-1 path has one node); anything
    else ⇒ ValueError.  Never parses a float from the environment."""
    raw = os.environ.get(NODE_ENV, "").strip()
    if raw == "":
        return C_IN_F_HUMAN_PF
    if raw == "d2root":
        if _degree_knob() != 2:
            raise ValueError(
                f"{NODE_ENV}=d2root needs {DEGREE_ENV}=2: the degree-1 path has "
                f"one node ({C_IN_F_HUMAN_PF} pF)"
            )
        return C_IN_F_HUMAN_D2_PF
    raise ValueError(f"{NODE_ENV}={raw!r}: only unset/empty or 'd2root' are accepted")


def _node_json_path(degree, c_pf=None):
    """``NODE_JSON`` for the degree-1 path; a ``-d2`` sibling at degree 2 so a
    degree-2 run never overwrites the degree-1 record items 359 / 361 read.

    Step 4c: ``c_pf`` (default ``None``) names the node; ``None`` or
    ``C_IN_F_HUMAN_PF`` give the paths above character for character, another
    node at degree 2 gives ``h9p-node-<c:.4f>-d2.json`` under ``JSON_DIR`` (read
    at call time); another node at degree 1 ⇒ ValueError."""
    if c_pf is not None and c_pf != C_IN_F_HUMAN_PF:
        if degree == 2:
            return os.path.join(JSON_DIR, f"h9p-node-{c_pf:.4f}-d2.json")
        raise ValueError((degree, c_pf))
    if degree is None:
        return NODE_JSON
    if degree == 2:
        return os.path.join(JSON_DIR, f"h9p-node-{C_IN_F_HUMAN_PF:.4f}-d2.json")
    raise ValueError(degree)


def _apply_degree(ctx, comm):
    """Apply the knob to ``ctx`` and assert (Deg).  Returns ``(solver, degree,
    counts)``; on the unset path ``ctx`` is not touched."""
    knob = _degree_knob()
    if knob is None:
        solver, degree = ctx["solver"], 1
    else:
        from fem_em_solver.core import TimeHarmonicSolver

        problem = ctx["solver"].problem
        del ctx["solver"]
        ctx["solver"] = TimeHarmonicSolver(problem, degree=2)
        solver, degree = ctx["solver"], 2
    msh = solver.mesh
    v_space = solver.function_space()
    imap = v_space.dofmap.index_map
    n_dofs = int(imap.size_global * v_space.dofmap.index_map_bs)
    for dim in (1, 2):
        if msh.topology.index_map(dim) is None:
            msh.topology.create_entities(dim)
    n_edges = int(msh.topology.index_map(1).size_global)
    n_faces = int(msh.topology.index_map(2).size_global)
    expected = n_edges if degree == 1 else 2 * n_edges + 2 * n_faces
    if comm.rank == 0:
        print(
            f"[WF-10 step4a] (Deg) knob {knob}; solver.degree {solver.degree}; global unknowns "
            f"{n_dofs} vs {'E' if degree == 1 else '2E+2F'} = {expected} (E {n_edges}, F "
            f"{n_faces}; records E {E_RECORD}, F {F_RECORD}, 2E+2F {N_DOFS_D2_RECORD}; "
            f"ASSERTED exact)",
            flush=True,
        )
    assert solver.degree == degree, (solver.degree, degree)
    assert n_dofs == expected, (n_dofs, expected, n_edges, n_faces)
    return solver, degree, {"n_dofs": n_dofs, "n_edges": n_edges, "n_faces": n_faces}


def _signed_k_shares(currents, order, slots):
    """``{ring sign: {"lead": |X₁|²/Σ, "lag": |X₁₅|²/Σ}}`` of the 16 slot currents.

    `_k1_shares`' Parseval (`tests/validation/test_wf10_f_human_tuned_circuit.py:99`)
    split by sign: a lagging ``e^{−j2πk/16}`` pattern lands in bin 15.
    """
    out = {}
    for sign in (1.0, -1.0):
        vec = np.zeros(16, dtype=np.complex128)
        for idx, pid in enumerate(order):
            k, sg = slots[pid]
            if sg == sign:
                vec[k] = currents[idx]
        spec = np.fft.fft(vec)
        total = float(np.sum(np.abs(spec) ** 2))
        if total:
            out[sign] = {
                "lead": float(abs(spec[1]) ** 2 / total),
                "lag": float(abs(spec[15]) ** 2 / total),
            }
        else:
            out[sign] = {"lead": float("nan"), "lag": float("nan")}
    return out


def _sense_weights():
    """``{sense: (w₁₇, w₂₁)}`` from the convention code at slots 0 / 4 (k = 0, 1),
    normalised to the first entry."""
    out = {}
    for sense in ("ccw", "cw"):
        w = quadrature_phase_weights([0, 1], sense)
        out[sense] = tuple(complex(x) for x in (w / w[0]))
    return out


def test_wf10_step2_sense_coding():
    """(S0): convention code, FFT bins, and the signed share helper."""
    w = _sense_weights()
    d_ccw = max(abs(np.asarray(w["ccw"]) - np.array([1.0, -1.0j])))
    d_ref = max(abs(np.asarray(w["ccw"]) - np.asarray(SENSES["(1,-1j)"])))
    print(
        f"[WF-10 step2 S0] w_ccw {w['ccw']}, w_cw {w['cw']}; |w_ccw - (1,-1j)| {d_ccw:.3e}, "
        f"vs SENSES['(1,-1j)'] {d_ref:.3e} (band 1e-15, ASSERTED)"
    )
    assert d_ccw <= 1e-15 and d_ref <= 1e-15
    bins = {}
    for sense in ("ccw", "cw"):
        x = np.zeros(16, dtype=np.complex128)
        x[0], x[4] = 1.0, w[sense][1]
        spec = np.fft.fft(x)
        bins[sense] = (abs(spec[1]), abs(spec[15]))
        print(f"[WF-10 step2 S0] {sense}: |X1| {bins[sense][0]:.3e}, |X15| {bins[sense][1]:.3e}")
    assert bins["ccw"][0] <= 1e-12 and abs(bins["ccw"][1] - 2.0) <= 1e-12
    assert bins["cw"][1] <= 1e-12 and abs(bins["cw"][0] - 2.0) <= 1e-12
    # Planted patterns (asserted negative/positive controls).
    order = [f"P{i}" for i in range(32)]
    slots = {f"P{i}": (i % 16, 1.0 if i < 16 else -1.0) for i in range(32)}
    k = np.array([slots[p][0] for p in order], dtype=float)
    lag = _signed_k_shares(np.exp(-2j * np.pi * k / 16.0), order, slots)
    lead = _signed_k_shares(np.exp(+2j * np.pi * k / 16.0), order, slots)
    print(f"[WF-10 step2 S0] planted lag pattern -> {lag}; planted lead pattern -> {lead}")
    for ring in (1.0, -1.0):
        assert abs(lag[ring]["lag"] - 1.0) <= 1e-12 and lag[ring]["lead"] <= 1e-12
        assert abs(lead[ring]["lead"] - 1.0) <= 1e-12 and lead[ring]["lag"] <= 1e-12


def build_tuned_quadrature(comm, *, with_z50=True, build_only=False):
    """Fixture, (A0)/(iv), optional Z50 anchors, the node pair at C_in, and the
    two sense superpositions with their partitions.  Never a ``test_*``.

    Returned keys (fixed for item 359): ``built``, ``ctx``, ``order``,
    ``slots``, ``feeds``, ``az_feeds``, ``az_mirror``, ``w``, ``e``, ``e17``,
    ``e21``, ``sigma``, ``s_node``, ``cur``, ``y_in``, ``parts``,
    ``p_rc_node``, ``blind_node``, ``times``, ``cells`` — plus ``anchors``
    (the asserted scalars), ``say`` and ``degree``.

    ``build_only=True`` (step 4a) returns ``built``, ``ctx``, ``cells``,
    ``degree``, ``n_dofs``, ``n_edges``, ``n_faces``, ``times`` and
    ``solver_built`` (the solver ``_build_f_human`` created; ``None`` at degree
    2) right after (Deg): no form compiled, nothing assembled or factorised.
    """
    import ufl
    from dolfinx import fem

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
    from tests.validation.test_port_birdcage_ring_column import AZIMUTH_MATCH_DEG
    from tests.validation.test_port_birdcage_ring_matrix import REUSE_REPRODUCTION_RTOL
    from tests.validation.test_port_circuit_layer_field import _residual
    from tests.validation.test_port_drive_superposition import _ring_quadrature_slots
    from tests.validation.test_port_package_sparameters import REFERENCE_IMPEDANCE_OHM
    from tests.validation.test_th14_birdcage_copper import _leontovich_term, _surface_loss_w
    from tests.validation.test_wf7_f_human_b1_quadrature import (
        DEGREE,
        FREQUENCY_HZ,
        _build_f_human,
    )
    from tests.validation.test_wf10_f_human_tuned_circuit import (
        QUADRATURE_SLOT_OFFSET,
        _respec,
    )
    from tests.validation.test_wf12_f_human_copper_matrix import (
        BLIND_FACTOR,
        COPPER_SIGMA,
        HOLE_CELL_RECORD,
        P17_P_SRC_RECORD_W,
        P17_P_SURF_RECORD_W,
    )
    from tests.validation.test_wf12_f_human_hole import _census_checks, _owned_census

    z0 = float(REFERENCE_IMPEDANCE_OHM)
    t_window = time.perf_counter()
    times = {}

    def say(msg):
        if comm.rank == 0:
            print(f"[WF-10 step2] {msg}", flush=True)

    # ---- the stored comparand (ruling (5): absent ⇒ stop, never re-solve) -----
    assert os.path.exists(
        NPZ_PATH
    ), f"{NPZ_PATH} absent: item 358 stops (BLOCKED, rule (d)); never re-solve"
    with np.load(NPZ_PATH, allow_pickle=False) as npz:
        s32_npz = np.array(npz["s32"], dtype=np.complex128)
        order_npz = [str(p) for p in npz["order"]]
        cells_npz = int(npz["cells"])
    c_pf = _node_pf()
    node_knob = os.environ.get(NODE_ENV, "").strip()
    c_f = c_pf * 1e-12
    knob = _degree_knob()
    eff_degree = DEGREE if knob is None else knob
    node_json = _node_json_path(knob, c_pf=c_pf)
    say(
        f"start: -n {comm.size}, f = {FREQUENCY_HZ:.3e} Hz, degree {eff_degree} ({DEGREE_ENV} "
        f"{knob if knob is not None else 'unset'}; fixture DEGREE {DEGREE}), kappa off, r = "
        f"{ESR_OHM} Ohm on all 32 gaps; node C_in,F-human {C_IN_F_HUMAN_PF} pF (secant "
        f"{C_IN_SECANT_PF:.6f} pF, (1,-1j) root {C_IN_CCW_ROOT_PF:.6f} pF; printed)"
    )
    say(
        f"step 4c node: {NODE_ENV} {node_knob or 'unset'} -> effective node {c_pf} pF "
        f"(degree-2 root {C_IN_F_HUMAN_D2_PF} pF, secant {C_IN_SECANT_D2_PF:.6f} pF); node "
        f"record {os.path.basename(node_json)}"
    )

    # ---- fixture, (A0), (iv) ----------------------------------------------
    t0 = time.perf_counter()
    built = _build_f_human(comm, as_hole=True)
    times["build"] = time.perf_counter() - t0
    ctx = built["ctx"]
    msh, tags_f, omega, cell_tags = ctx["msh"], ctx["tags_f"], ctx["omega"], ctx["cell_tags"]
    assert abs(omega / (2.0 * np.pi * float(FREQUENCY_HZ)) - 1.0) < 1e-12, omega
    census = _owned_census(msh, cell_tags, comm)
    census_fails = _census_checks(census, built["cells"], _ring_ports())
    cell_rel = built["cells"] / HOLE_CELL_RECORD - 1.0
    say(
        f"(A0) hole cells {built['cells']} vs record {HOLE_CELL_RECORD}, rel {cell_rel:.3e} "
        f"(band {CELL_COUNT_BAND}, ASSERTED); npz cells {cells_npz}; build {times['build']:.1f} s; "
        f"(iv) failures {census_fails or 'none'}"
    )
    assert abs(cell_rel) <= CELL_COUNT_BAND, cell_rel
    assert built["cells"] == cells_npz, (built["cells"], cells_npz)
    assert not census_fails, census_fails

    # ---- step 4a: the degree knob, (Deg) asserted inside ----------------------
    solver_built = ctx["solver"]
    _solver, degree, counts = _apply_degree(ctx, comm)
    del _solver
    if degree != 1:
        solver_built = None  # never hold two solvers (`WF-12` step 4's trap)
    if build_only:
        say(
            f"build_only: return after (A0)/(iv)/(Deg) at degree {degree}; no form compiled, "
            f"nothing assembled; build {times['build']:.1f} s"
        )
        return {
            "built": built,
            "ctx": ctx,
            "cells": built["cells"],
            "degree": degree,
            "n_dofs": counts["n_dofs"],
            "n_edges": counts["n_edges"],
            "n_faces": counts["n_faces"],
            "times": times,
            "solver_built": solver_built,
            "c_in_pf": c_pf,
            "node_json": node_json,
        }
    del solver_built

    z_s = (1.0 + 1.0j) * surface_resistance_ohm(omega, COPPER_SIGMA)
    wall = _leontovich_term(msh, built["facets"], z_s, omega)
    z_p = complex(TERMINATED_PORT_IMPEDANCE_OHM)
    dx = ufl.Measure("dx", domain=msh, subdomain_data=cell_tags)
    acct = {"mesh": msh, "facet_tags": tags_f}
    solver = ctx["solver"]
    base_specs = list(ctx["specs"])
    keep = {}

    def _vol(sigma_f, e, measure):
        form = fem.form(0.5 * sigma_f * ufl.inner(e, e) * measure)
        return float(np.real(comm.allreduce(fem.assemble_scalar(form), op=MPI.SUM)))

    # Lifted from `main` at this item's base (524bd617),
    # tests/validation/test_wf10_f_human_copper_esr_root.py lines 339–386 (the
    # `_drive` closure), verbatim except: the solved field is copied into a fresh
    # Function and kept when ``keep_as`` is given (the held solve may reuse the
    # solver's Function), and the σ field is kept once.
    def _drive(specs, driven_id, *, held, partition, keep_as=None):
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
            bil = [wall] + [
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
        out = {
            "s": s_col,
            "currents": currents,
            "kind": "held" if held else "solve",
            "t_solve": t_solve,
        }
        if "sigma" not in keep:
            keep["sigma"] = fields.sigma_field
        if partition:
            sigma_f = fields.sigma_field
            by_id = {
                s.port_id: float(_sheet_field_dissipation_w(acct, s, e, omega)) for s in sheets
            }
            out["p_sheet_by_id"] = by_id
            out["p_sheets"] = float(sum(by_id.values()))
            out["p_src"] = float(_source_power_w(acct, dsheet, e, omega))
            p_total = _vol(sigma_f, e, dx)
            out["p_phantom"] = _vol(sigma_f, e, dx(PHANTOM_CELL_TAG))
            out["p_rest"] = p_total - out["p_phantom"]
            out["p_surf"] = _surface_loss_w(msh, built["facets"], z_s, e, comm)
        if keep_as is not None:
            cp = fem.Function(e.function_space, name=keep_as)
            cp.x.array[:] = e.x.array
            cp.x.scatter_forward()
            keep[keep_as] = cp
        del fields, e
        out["t_post"] = time.perf_counter() - t1
        return out

    # ---- order (H17's: [P17, feed2] + rest by number) -------------------------
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
    feeds_idx = [0, 1]
    assert len(order) == 32 and sorted(order) == sorted(sp.port_id for sp in base_specs)
    assert order_npz == order, "npz order differs from this run's order"
    az_by = {f"P{s['ordinal']}": float(s["azimuth_deg"]) for s in built["sheets"]}
    az_feeds = (az_by[feeds[0]], az_by[feeds[1]])
    d_az = (az_feeds[1] - az_feeds[0]) % 360.0
    az_mirror = az_feeds[0] + 45.0
    say(
        f"feeds: {feeds[0]} (slot {ref_k}, ring {ref_sign:+.0f}, az {az_feeds[0]:.9f} deg) and "
        f"{feeds[1]} (slot {slots[feeds[1]][0]}, az {az_feeds[1]:.9f} deg); az21 - az17 = "
        f"{d_az:.9f} deg (90 within {AZIMUTH_MATCH_DEG:g}, ASSERTED); mirror plane az "
        f"{az_mirror:.9f} deg; az_ref {az_ref:.9f} deg"
    )
    assert abs(d_az - 90.0) < AZIMUTH_MATCH_DEG, d_az
    assert slots[feeds[0]][0] == 0 and slots[feeds[1]][0] == 4, (slots[feeds[0]], slots[feeds[1]])

    # Lifted `_pair` (H17 module lines 406–418), verbatim except that the field
    # of each drive is kept when ``keep_fields`` is set.
    def _pair(imp_x, z_feed, v_feed, keep_fields=False):
        specs_x = _respec(base_specs, None, imp_x)
        assert all(sp.width_correction_kappa is None for sp in specs_x)
        cur, knd, tt, parts = [], [], [], []
        for n_drive, pid in enumerate(feeds):
            d = _drive(
                specs_x,
                pid,
                held=(n_drive > 0),
                partition=(n_drive == 0),
                keep_as=(f"E_{pid}" if keep_fields else None),
            )
            knd.append(d["kind"])
            tt.append(d["t_solve"] + d["t_post"])
            cur.append(np.array([d["currents"][p] for p in order], dtype=np.complex128))
            if n_drive == 0:
                parts.append({k: v for k, v in d.items() if k.startswith("p_")})
            del d
        return _current_route_s(cur, feeds_idx, z_feed, v_feed), cur, knd, tt, parts[0], specs_x

    def _partition(p):
        p_vol = p["p_phantom"] + p["p_rest"]
        rel = abs(p["p_src"] - p["p_sheets"] - p_vol - p["p_surf"]) / abs(p["p_src"])
        blind = abs(p["p_src"] - p["p_sheets"] - p_vol) / abs(p["p_src"])
        return rel, blind

    anchors = {}
    # ---- the Z50 pair: (R), (P50), NEG, (Z50) — H17's, unchanged --------------
    if with_z50:
        s_50, _c50, kinds_50, t_50, p50, _sp = _pair({pid: z0 for pid in order}, z0, 1.0)
        del _c50, _sp
        rel_p50, blind50 = _partition(p50)
        r_src = abs(p50["p_src"] / P17_P_SRC_RECORD_W - 1.0)
        r_surf = abs(p50["p_surf"] / P17_P_SURF_RECORD_W - 1.0)
        r_z50 = _residual(s_50, s32_npz[:2, :2])
        times["z50_pair"] = float(sum(t_50))
        say(f"Z50 pair {kinds_50}, {', '.join(f'{t:.2f}' for t in t_50)} s")
        say(
            f"(R) P_src rel {r_src:.3e}; P_surf rel {r_surf:.3e} (band {REUSE_REPRODUCTION_RTOL:g}, "
            f"ASSERTED; measured 2.109e-11 / 1.509e-11)"
        )
        say(
            f"(P50) rel {rel_p50:.3e} (band {DISCRETE_IDENTITY_RTOL:g}, ASSERTED; measured "
            f"2.417e-13)"
        )
        say(
            f"NEG P_surf-blind miss {blind50:.4e} = {blind50 / DISCRETE_IDENTITY_RTOL:.1f}x the band "
            f"(bar {BLIND_FACTOR:.0f}x, ASSERTED; measured 123.6x)"
        )
        say(
            f"(Z50) residual {r_z50:.6e} (band {REUSE_REPRODUCTION_RTOL:g}, ASSERTED; measured "
            f"2.981e-14); t={time.perf_counter() - t_window:.1f} s"
        )
        assert kinds_50 == ["solve", "held"], kinds_50
        assert r_src <= REUSE_REPRODUCTION_RTOL, r_src
        assert r_surf <= REUSE_REPRODUCTION_RTOL, r_surf
        assert rel_p50 <= DISCRETE_IDENTITY_RTOL, rel_p50
        assert blind50 >= BLIND_FACTOR * DISCRETE_IDENTITY_RTOL, blind50
        assert r_z50 <= REUSE_REPRODUCTION_RTOL, f"(Z50) closure {r_z50:.6e}"
        anchors.update(r_src=r_src, r_surf=r_surf, rel_p50=rel_p50, blind50=blind50, r_z50=r_z50)

    # ---- the node pair at C_in -------------------------------------------------
    z_rc = _z_rc(c_f, FREQUENCY_HZ)
    z_th = _thevenin_feed_impedance(z_rc, z0)
    v_s = _source_voltage(z_rc, z0)
    say(f"node C = {c_f:.12e} F: Z_RC {z_rc:.9e} Ohm, Z_th {z_th:.9e} Ohm, V_s {v_s:.9e} V")
    imp = {pid: z_rc for pid in others}
    imp.update({pid: z_th for pid in feeds})
    s_node, cur_node, kinds_n, t_n, pn, specs_n = _pair(imp, z_th, v_s, keep_fields=True)
    times["node_pair"] = float(sum(t_n))
    assert kinds_n == ["solve", "held"], kinds_n
    p_rc_node, blind_node = _partition(pn)
    p_rc_sheets = sum(v for pid, v in pn["p_sheet_by_id"].items() if pid not in feeds)
    p_vol_n = pn["p_phantom"] + pn["p_rest"]
    dz = abs(pn["p_src"] - (pn["p_sheets"] - p_rc_sheets) - p_vol_n - pn["p_surf"]) / abs(
        pn["p_src"]
    )
    say(f"node pair {kinds_n}, {', '.join(f'{t:.2f}' for t in t_n)} s")
    say(
        f"(P_RC) P17 drive: P_src {pn['p_src']:.9e} sumP_sheet {pn['p_sheets']:.9e} (30 RC "
        f"{p_rc_sheets:.9e}) P_phantom {pn['p_phantom']:.9e} P_rest {pn['p_rest']:.9e} P_surf "
        f"{pn['p_surf']:.9e} W; rel {p_rc_node:.3e} (band {DISCRETE_IDENTITY_RTOL:g}, ASSERTED; "
        f"flanking 7.878e-13 / 7.726e-13)"
    )
    say(
        f"NEG at the node: P_surf-blind miss {blind_node:.4e} = "
        f"{blind_node / DISCRETE_IDENTITY_RTOL:.1f}x the band (bar {BLIND_FACTOR:.0f}x, ASSERTED)"
    )
    say(
        f"(Dz) partition without the 30 RC sheets misses {dz:.4e} (PREDICTED ≈ {PREDICTED_DZ}, "
        f"printed)"
    )
    assert p_rc_node <= DISCRETE_IDENTITY_RTOL, f"(P_RC) {p_rc_node:.3e} — never re-band"
    assert blind_node >= BLIND_FACTOR * DISCRETE_IDENTITY_RTOL, blind_node
    anchors.update(p_rc_node=p_rc_node, blind_node=blind_node, dz=dz)

    e17, e21 = keep.pop(f"E_{feeds[0]}"), keep.pop(f"E_{feeds[1]}")
    sigma = keep["sigma"]
    # Shape, mesh and element checks — never object identity (POST-6 step 1's trap).
    assert e17.function_space.mesh is e21.function_space.mesh is msh
    assert e17.x.array.shape == e21.x.array.shape
    assert str(e17.function_space.ufl_element()) == str(e21.function_space.ufl_element())
    assert (
        e17.function_space.dofmap.index_map.size_local
        == e21.function_space.dofmap.index_map.size_local
    )

    # ---- the two senses: superposition, partition, (Lin), Y_in ---------------
    w = _sense_weights()
    specs_by = {sp.port_id: sp for sp in specs_n}
    e_s, cur_s, y_in, parts, lin = {}, {}, {}, {}, {}
    t0 = time.perf_counter()
    for sense in ("ccw", "cw"):
        w17, w21 = w[sense]
        ef = fem.Function(e17.function_space, name=f"E_{sense}")
        ef.x.array[:] = w17 * e17.x.array + w21 * e21.x.array
        ef.x.scatter_forward()
        e_s[sense] = ef
        wf = {feeds[0]: w17, feeds[1]: w21}
        cur_s[sense] = w17 * cur_node[0] + w21 * cur_node[1]
        sheets = {
            pid: (
                dataclasses.replace(
                    specs_by[pid].sheet(driven=True), source_voltage_v=complex(wf[pid])
                )
                if pid in feeds
                else specs_by[pid].sheet(driven=False)
            )
            for pid in order
        }
        cur_read = np.array(
            [sheet_terminal_current(msh, tags_f, sheets[pid], ef, comm) for pid in order],
            dtype=np.complex128,
        )
        r_lin = float(np.max(np.abs(cur_read - cur_s[sense])) / np.max(np.abs(cur_s[sense])))
        y_cur = [
            complex(1.0 / z_rc + cur_read[i] / (wf[pid] - z_th * cur_read[i]))
            for i, pid in enumerate(feeds)
        ]
        y_ref = _active_y_in(s_node, (w17, w21), z0)
        r_y = float(max(abs(a - b) for a, b in zip(y_cur, y_ref)) / max(abs(b) for b in y_ref))
        y_in[sense] = y_ref
        p = {
            "p_src": float(sum(_source_power_w(acct, sheets[pid], ef, omega) for pid in feeds)),
            "p_sheet_by_id": {
                pid: float(_sheet_field_dissipation_w(acct, sheets[pid], ef, omega))
                for pid in order
            },
        }
        p["p_sheets"] = float(sum(p["p_sheet_by_id"].values()))
        p_total = _vol(sigma, ef, dx)
        p["p_phantom"] = _vol(sigma, ef, dx(PHANTOM_CELL_TAG))
        p["p_rest"] = p_total - p["p_phantom"]
        p["p_surf"] = _surface_loss_w(msh, built["facets"], z_s, ef, comm)
        rel_w, blind_w = _partition(p)
        p["rel"], p["blind"] = rel_w, blind_w
        parts[sense] = p
        lin[sense] = (r_lin, r_y)
        say(
            f"[{sense}] w = ({w17:.3f}, {w21:.3f}); (P_w) P_src {p['p_src']:.9e} sumP_sheet "
            f"{p['p_sheets']:.9e} P_phantom {p['p_phantom']:.9e} P_rest {p['p_rest']:.9e} P_surf "
            f"{p['p_surf']:.9e} W; rel {rel_w:.3e} (band {DISCRETE_IDENTITY_RTOL:g}, ASSERTED; "
            f"PREDICTED <= {PREDICTED_P_W:g}); P_surf-blind miss {blind_w:.4e} (PREDICTED ≈ "
            f"{PREDICTED_BLIND_PW:g}, printed)"
        )
        say(
            f"[{sense}] (Lin) currents rel {r_lin:.3e} (band 1e-12 UNIT_WEIGHT_CONTROL_RTOL, "
            f"ASSERTED; PREDICTED <= {PREDICTED_LIN:g}); Y_in current route vs _active_y_in "
            f"{r_y:.3e} (band {Y0_CODING_BAND:g}, ASSERTED)"
        )
    times["superposition"] = time.perf_counter() - t0
    anchors["lin"] = lin
    anchors["p_w"] = {s: (parts[s]["rel"], parts[s]["blind"]) for s in parts}

    rec = {
        "c_f": c_f,
        "c_in_pf": c_pf,
        "r_ohm": ESR_OHM,
        "freq_hz": float(FREQUENCY_HZ),
        "order": order,
        "feeds": feeds,
        "az_feeds": list(az_feeds),
        "az_mirror": az_mirror,
        "s_node": [[[complex(v).real, complex(v).imag] for v in row] for row in s_node],
        "w": {s: [[complex(x).real, complex(x).imag] for x in w[s]] for s in w},
        "y_in": {s: [[y.real, y.imag] for y in ys] for s, ys in y_in.items()},
        "partitions": {s: parts[s] for s in parts},
        "partition_p17_node": dict(pn),
        "p_rc_node": p_rc_node,
        "blind_node": blind_node,
        "times_s": dict(times),
        "degree": degree,
        "node_knob": node_knob,
    }
    if comm.rank == 0:
        os.makedirs(JSON_DIR, exist_ok=True)
        with open(node_json, "w") as fh:
            json.dump(rec, fh, indent=1)
        print(f"[WF-10 step2] wrote {node_json}", flush=True)
    comm.Barrier()

    return {
        "built": built,
        "ctx": ctx,
        "order": order,
        "slots": slots,
        "feeds": feeds,
        "az_feeds": az_feeds,
        "az_mirror": az_mirror,
        "w": w,
        "e": e_s,
        "e17": e17,
        "e21": e21,
        "sigma": sigma,
        "s_node": s_node,
        "cur": cur_s,
        "y_in": y_in,
        "parts": parts,
        "p_rc_node": p_rc_node,
        "blind_node": blind_node,
        "times": times,
        "cells": built["cells"],
        "anchors": anchors,
        "say": say,
        "omega": omega,
        "t_window": t_window,
        "degree": degree,
        "c_in_pf": c_pf,
        "node_json": node_json,
    }


@complex_only
def test_wf10_step2_tuned_b1_quadrature():
    if os.environ.get(STEP2_ENV, "").strip() != "1":
        pytest.skip(f"{STEP2_ENV} unset: heavy -n 8 F-human window, environment-selected")

    from fem_em_solver.post.faraday import magnetic_flux_density_from_e, project_to_cg1
    from fem_em_solver.post.sar import mean_sar
    from tests.validation.test_birdcage_b1_plus_map import (
        C4_COVARIANCE_BAND,
        PHANTOM_RHO_KG_PER_M3,
        _relative_l2,
        _rotate_z,
        _sample_points,
    )
    from tests.validation.test_birdcage_b1_quadrature import _mirror_xy, _read_senses
    from tests.validation.test_birdcage_power_identity import DISCRETE_IDENTITY_RTOL
    from tests.validation.test_birdcage_sar_mass_averaged import EXACT_IDENTITY_RTOL
    from tests.validation.test_port_birdcage_lumped_column import PHANTOM_CELL_TAG
    from tests.validation.test_port_birdcage_ring_matrix import REUSE_REPRODUCTION_RTOL
    from tests.validation.test_port_circuit_layer_field import _residual
    from tests.validation.test_port_drive_superposition import UNIT_WEIGHT_CONTROL_RTOL
    from tests.validation.test_wf10_f_human_tuned_circuit import _k1_shares
    from tests.validation.test_wf12_f_human_hole import _rss_gib

    comm = MPI.COMM_WORLD
    knob = _degree_knob()
    # ---- step 4a (Bit): the degree-1 node record, read BEFORE the builder ------
    # (which rewrites it), exactly as item 359 does.  Degree-1 comparand only.
    s_rec = None
    if knob is None:
        assert os.path.exists(NODE_JSON), f"{NODE_JSON} absent: the (S_rep) re-gate cannot run"
        rec = None
        if comm.rank == 0:
            with open(NODE_JSON) as fh:
                rec = json.load(fh)
        rec = comm.bcast(rec, root=0)
        assert abs(rec["c_in_pf"] - C_IN_F_HUMAN_PF) == 0.0, rec["c_in_pf"]
        s_rec = np.array(
            [[complex(a, b) for a, b in row] for row in rec["s_node"]], dtype=np.complex128
        )
    out = build_tuned_quadrature(comm, with_z50=(knob is None))
    say, omega, ctx = out["say"], out["omega"], out["ctx"]
    order, slots, feeds = out["order"], out["slots"], out["feeds"]
    times = out["times"]
    if knob is None:
        r_srep = _residual(np.asarray(out["s_node"], dtype=np.complex128), s_rec)
        say(
            f"(S_rep) s_node vs {os.path.basename(NODE_JSON)} read before the builder: residual "
            f"{r_srep:.3e} (band {REUSE_REPRODUCTION_RTOL:g}, ASSERTED)"
        )
        assert r_srep <= REUSE_REPRODUCTION_RTOL, f"(S_rep) {r_srep:.3e}"
    else:
        say(
            f"(R) / (Z50): degree-1 records, not run at degree 2 (effective node "
            f"{out['c_in_pf']} pF, record {os.path.basename(out['node_json'])}); (S_rep) likewise "
            f"(printed: {DEGREE2_RECORD_GATES['printed']}; asserted: "
            f"{DEGREE2_RECORD_GATES['asserted']})"
        )

    # ---- Y_in, four readings ---------------------------------------------------
    for sense, label in (("cw", "(1,+1j)"), ("ccw", "(1,-1j)")):
        for i, pid in enumerate(feeds):
            y = out["y_in"][sense][i]
            say(
                f"Y_in a = {label} [{sense}] {pid}: {y.real:+.9e} {y.imag:+.9e}j S (PREDICTED "
                f"{PREDICTED_Y_IN[(pid, label)]}; reading)"
            )
    yc, yw = out["y_in"]["ccw"], out["y_in"]["cw"]
    say(
        f"M-paired Y_in spreads: |Y17(ccw) - Y21(cw)|/|Y17(ccw)| {abs(yc[0] - yw[1]) / abs(yc[0]):.4e}; "
        f"|Y21(ccw) - Y17(cw)|/|Y21(ccw)| {abs(yc[1] - yw[0]) / abs(yc[1]):.4e} (PREDICTED "
        f"2.1-2.3 %, reading)"
    )

    # ---- signed k-spectrum -----------------------------------------------------
    for sense in ("ccw", "cw"):
        sh = _signed_k_shares(out["cur"][sense], order, slots)
        un = _k1_shares(out["cur"][sense], order, slots)
        say(
            f"[{sense}] signed k shares: ring+ lead {sh[1.0]['lead']:.6f} lag {sh[1.0]['lag']:.6f}; "
            f"ring- lead {sh[-1.0]['lead']:.6f} lag {sh[-1.0]['lag']:.6f}; unsigned +-1 ring+ "
            f"{un[1.0]:.6f} ring- {un[-1.0]:.6f} (H17 ≈ {H17_K1_SHARE}; PREDICTED lag share under "
            f"ccw >= {PREDICTED_LAG_SHARE}; reading)"
        )

    # ---- B, CG1 projection, point reads ------------------------------------------
    comm.Barrier()
    t0 = time.perf_counter()
    cg1 = {}
    for sense in ("ccw", "cw"):
        b_dg = magnetic_flux_density_from_e(out["e"][sense], omega, name=f"B_{sense}")
        cg1[sense] = project_to_cg1(b_dg, name=f"B_{sense}_cg1")
        del b_dg
    comm.Barrier()
    t_proj = time.perf_counter() - t0
    sweep = {"mesh": ctx["msh"], "cell_tags": ctx["cell_tags"]}
    points = _sample_points(sweep)
    az_m = out["az_mirror"]
    az17 = out["az_feeds"][0]
    comm.Barrier()
    t0 = time.perf_counter()
    ccw_x = _read_senses(cg1["ccw"], points)
    cw_x = _read_senses(cg1["cw"], points)
    cw_mx = _read_senses(cg1["cw"], _mirror_xy(points, az_m))
    cw_m17 = _read_senses(cg1["cw"], _mirror_xy(points, az17))
    rot = {
        m: _read_senses(cg1["ccw"], _rotate_z(points, np.radians(m * 22.5))) for m in range(1, 16)
    }
    comm.Barrier()
    t_eval = time.perf_counter() - t0
    mask = np.logical_and.reduce(
        [ccw_x[2], cw_x[2], cw_mx[2], cw_m17[2]] + [r[2] for r in rot.values()]
    )
    n_pts, n_mask = int(points.shape[0]), int(np.count_nonzero(mask))
    plus_x, minus_x = ccw_x[0], ccw_x[1]
    mirror = _relative_l2(cw_mx[1], plus_x, mask)
    mispaired = _relative_l2(cw_mx[0], plus_x, mask)
    mispaired_pred = _relative_l2(minus_x, plus_x, mask)
    wrong_plane = _relative_l2(cw_m17[1], plus_x, mask)
    c16 = {m: _relative_l2(rot[m][0], plus_x, mask) for m in rot}
    worst_m = max(c16, key=c16.get)
    sense_ratio = float(np.mean(plus_x[mask]) / np.mean(minus_x[mask]))
    cw_ratio = float(np.mean(cw_x[1][mask]) / np.mean(cw_x[0][mask]))
    r_pts = np.hypot(points[:, 0], points[:, 1]) + np.abs(points[:, 2])
    r_pts = np.where(mask, r_pts, np.inf)
    ic = int(np.argmin(r_pts))
    mean_b = float(np.mean(plus_x[mask]))
    cv_b = float(np.std(plus_x[mask]) / mean_b)
    p_src_ccw = out["parts"]["ccw"]["p_src"]
    say(
        f"B: two CG1 projections {t_proj:.2f} s; {n_pts} sample points, {n_mask} valid on every "
        f"image; evaluations {t_eval:.2f} s"
    )
    say(
        f"(Sense) mean|B1+|_ccw / mean|B1-|_ccw = {sense_ratio:.4f} (> 1 ASSERTED; PREDICTED >= "
        f"{PREDICTED_SENSE_RATIO:g}); cw: mean|B1-|/mean|B1+| = {cw_ratio:.4f} (printed)"
    )
    say(
        f"centre purity at ({points[ic, 0]:.4e}, {points[ic, 1]:.4e}, {points[ic, 2]:.4e}) m: ccw "
        f"|B1+|/|B1-| {plus_x[ic] / minus_x[ic]:.4f}; cw |B1-|/|B1+| {cw_x[1][ic] / cw_x[0][ic]:.4f} "
        f"(reading)"
    )
    say(
        f"(M_B) feed-bisector mirror (az {az_m:.6f} deg): |B1-|_cw(Mx) vs |B1+|_ccw(x) rel L2 "
        f"{mirror:.4e} (band {C4_COVARIANCE_BAND:g}, ASSERTED; PREDICTED {PREDICTED_MIRROR[0]:g}-"
        f"{PREDICTED_MIRROR[1]:g})"
    )
    say(
        f"mis-paired control |B1+|_cw(Mx) vs |B1+|_ccw(x) {mispaired:.4e} beside its prediction "
        f"|B1-|_ccw(x) vs |B1+|_ccw(x) {mispaired_pred:.4e} (PREDICTED >= {PREDICTED_MISPAIRED:g}, "
        f"printed)"
    )
    say(f"wrong-plane mirror (az17 {az17:.6f} deg): {wrong_plane:.4e} (printed, reading)")
    say(
        "C16 readings |B1+|_ccw(R_m x) vs |B1+|_ccw(x): "
        + ", ".join(f"m={m}: {v:.4e}" for m, v in c16.items())
        + f"; worst m={worst_m} {c16[worst_m]:.4e} (PREDICTED <= {PREDICTED_C16:g}; printed, "
        f"never asserted — the drive is not C16-symmetric)"
    )
    say(
        f"step 4a records (printed beside): (M_B) {mirror:.4e} vs 358's {M_B_RECORD:.4e}; mean"
        f"|B1+|_ccw {mean_b:.6e} T vs 358's {MEAN_B1P_RECORD_T:.6e} T (degree {out['degree']})"
    )
    say(
        f"mean|B1+|_ccw {mean_b:.6e} T, CV {cv_b:.4e} (reading, no homogeneity claim); "
        f"mean|B1+|_ccw/sqrt(P_src) {mean_b / np.sqrt(p_src_ccw):.6e} T/sqrt(W)"
    )

    # ---- SAR, (Σ_P) ----------------------------------------------------------
    sar = mean_sar(
        out["e"]["ccw"],
        sigma=out["sigma"],
        rho=PHANTOM_RHO_KG_PER_M3,
        cell_tags=ctx["cell_tags"],
        subdomain_ids=PHANTOM_CELL_TAG,
        comm=comm,
    )
    p_ph = out["parts"]["ccw"]["p_phantom"]
    r_sigma = abs(sar["dissipated_power_w"] - p_ph) / abs(p_ph)
    say(
        f"(Sigma_P) P_phantom {p_ph:.12e} W vs mean_sar dissipated {sar['dissipated_power_w']:.12e} "
        f"W, rel {r_sigma:.3e} (band {EXACT_IDENTITY_RTOL:g}, ASSERTED; PREDICTED <= "
        f"{PREDICTED_SIGMA_P:g})"
    )
    say(
        f"SAR (ccw): whole-phantom mean {sar['mean_sar_w_per_kg']:.6e} W/kg, per W of P_src "
        f"{sar['mean_sar_w_per_kg'] / p_src_ccw:.6e} W/kg/W; P_phantom/P_src {p_ph / p_src_ccw:.6e}"
    )

    _m, rss_per, rss_sum = _rss_gib(comm)
    say(
        f"price: build {times['build']:.1f} s, Z50 pair {times.get('z50_pair', 0.0):.1f} s, node "
        f"pair {times['node_pair']:.1f} s, superposition + partitions {times['superposition']:.1f} s, "
        f"B projections {t_proj:.1f} s + evaluations {t_eval:.1f} s; RSS summed {rss_sum:.2f} GiB; "
        f"window elapsed {time.perf_counter() - out['t_window']:.1f} s at -n {comm.size}"
    )

    # ---- asserts ---------------------------------------------------------------
    for sense in ("ccw", "cw"):
        r_lin, r_y = out["anchors"]["lin"][sense]
        rel_w, _b = out["anchors"]["p_w"][sense]
        assert rel_w <= DISCRETE_IDENTITY_RTOL, f"(P_w) {sense} {rel_w:.3e} — never re-band"
        assert r_lin <= UNIT_WEIGHT_CONTROL_RTOL, f"(Lin) currents {sense} {r_lin:.3e}"
        assert r_y <= Y0_CODING_BAND, f"(Lin) Y_in {sense} {r_y:.3e}"
    assert r_sigma <= EXACT_IDENTITY_RTOL, f"(Sigma_P) {r_sigma:.3e}"
    assert n_mask > 0, "no sample point valid on every image"
    assert sense_ratio > 1.0, f"(Sense) red: B1- dominant under ccw ({sense_ratio:.4f})"
    assert mirror <= C4_COVARIANCE_BAND, f"(M_B) mirror {mirror:.4e} — never re-band"


def test_wf10_step4a_degree_knob_coding(monkeypatch):
    """(K0) step 4a: the knob parser, the node-JSON path, the degree-2 gate table."""
    monkeypatch.delenv(DEGREE_ENV, raising=False)
    assert _degree_knob() is None
    for raw, want in (("", None), ("  ", None), ("2", 2), (" 2 ", 2)):
        monkeypatch.setenv(DEGREE_ENV, raw)
        assert _degree_knob() == want, (raw, want)
    for raw in ("1", "3", "x", "2.0"):
        monkeypatch.setenv(DEGREE_ENV, raw)
        with pytest.raises(ValueError):
            _degree_knob()
    monkeypatch.delenv(DEGREE_ENV, raising=False)
    p1, p2 = _node_json_path(None), _node_json_path(2)
    print(f"[WF-10 step4a K0] node JSON: unset -> {p1}; degree 2 -> {p2}")
    assert p1 == NODE_JSON
    assert p2.endswith("-d2.json") and os.path.dirname(p2) == JSON_DIR and p2 != NODE_JSON
    with pytest.raises(ValueError):
        _node_json_path(1)
    printed, asserted = DEGREE2_RECORD_GATES["printed"], DEGREE2_RECORD_GATES["asserted"]
    print(f"[WF-10 step4a K0] degree-2 gates: printed {printed}; asserted {asserted}")
    assert "(R)" in printed and "(Z50)" in printed
    assert not set(printed) & set(asserted)
    for gate in (
        "(A0)",
        "(iv)",
        "(Deg)",
        "(P_RC)",
        "(P_w)",
        "(Lin)",
        "(Sigma_P)",
        "(Sense)",
        "(M_B)",
    ):
        assert gate in asserted, gate
    # The records' own arithmetic (exact integers, `…093243Z…:6221–6222`).
    assert N_DOFS_D2_RECORD == 2 * E_RECORD + 2 * F_RECORD


@complex_only
def test_wf10_step4a_build_only(monkeypatch):
    """(K1) step 4a: one build, both paths' unknown counts; no form compiled."""
    if os.environ.get(STEP4A_ENV, "").strip() != "1":
        pytest.skip(f"{STEP4A_ENV} unset: build-only -n 2 F-human window, environment-selected")
    from tests.validation.test_wf12_f_human_hole import _rss_gib

    comm = MPI.COMM_WORLD

    def say(msg):
        if comm.rank == 0:
            print(f"[WF-10 step4a K1] {msg}", flush=True)

    monkeypatch.delenv(DEGREE_ENV, raising=False)
    monkeypatch.delenv(NODE_ENV, raising=False)  # step 4c: the knob needs degree 2
    t0 = time.perf_counter()
    out = build_tuned_quadrature(comm, build_only=True)
    t_build = time.perf_counter() - t0
    ctx = out["ctx"]
    n1, e, f = out["n_dofs"], out["n_edges"], out["n_faces"]
    same = ctx["solver"] is out["solver_built"]
    say(
        f"(Deg-1) unset: degree {out['degree']}, solver is _build_f_human's: {same}; unknowns "
        f"{n1} vs E {e} (record E {E_RECORD}; ASSERTED exact); F {f} (record {F_RECORD}); "
        f"build {t_build:.1f} s"
    )
    assert out["degree"] == 1 and ctx["solver"].degree == 1
    assert same, "(Deg-1) the unset path rebuilt the solver: the knob leaks"
    assert n1 == e, (n1, e)
    out.pop("solver_built")

    monkeypatch.setenv(DEGREE_ENV, "2")
    t0 = time.perf_counter()
    solver2, deg2, counts2 = _apply_degree(ctx, comm)
    t_d2 = time.perf_counter() - t0
    n2 = counts2["n_dofs"]
    assert counts2["n_edges"] == e and counts2["n_faces"] == f
    _m, _per, rss_sum = _rss_gib(comm)
    say(
        f"(Deg-2) knob 2: degree {deg2} (solver.degree {solver2.degree}); unknowns {n2} vs "
        f"2E+2F = {2 * e + 2 * f} (record {N_DOFS_D2_RECORD}; ASSERTED exact); degree-2 space "
        f"{t_d2:.1f} s"
    )
    say(
        f"NEG (asserted, rule (e)): n(2) {n2} != n(1) {n1}; n(2)/n(1) = {n2 / n1:.6f} vs "
        f"(2E+2F)/E = {(2 * e + 2 * f) / e:.6f} (exact integer identity ASSERTED)"
    )
    say(
        f"RSS after the degree-2 space, summed {rss_sum:.2f} GiB (PREDICTED <= "
        f"{PREDICTED_RSS_D2_GIB:g}; WF-12 4a at build {RSS_D2_BUILD_RECORD_GIB} GiB; printed)"
    )
    assert deg2 == 2 and solver2.degree == 2 and ctx["solver"] is solver2
    assert n2 == 2 * e + 2 * f, (n2, e, f)
    assert n2 != n1
    assert n2 * e == n1 * (2 * e + 2 * f), (n1, n2, e, f)


def test_wf10_step4c_node_knob_coding(monkeypatch):
    """(N0) step 4c (§9 item 511, H28): the node knob, the node-record path and
    the series-RC closed form (N2) at both nodes.  numpy only."""
    from tests.validation.test_wf7_f_human_b1_quadrature import FREQUENCY_HZ

    monkeypatch.delenv(NODE_ENV, raising=False)
    monkeypatch.delenv(DEGREE_ENV, raising=False)
    assert _node_pf() is C_IN_F_HUMAN_PF
    monkeypatch.setenv(NODE_ENV, "")
    assert _node_pf() is C_IN_F_HUMAN_PF
    monkeypatch.setenv(NODE_ENV, "d2root")
    with pytest.raises(ValueError):
        _node_pf()  # degree knob unset: the degree-1 path has one node
    monkeypatch.setenv(DEGREE_ENV, "2")
    assert _node_pf() == 115.5209 and _node_pf() is C_IN_F_HUMAN_D2_PF
    for raw in ("115.5209", "x", "d2"):
        monkeypatch.setenv(NODE_ENV, raw)
        with pytest.raises(ValueError):
            _node_pf()
    monkeypatch.delenv(NODE_ENV, raising=False)
    monkeypatch.delenv(DEGREE_ENV, raising=False)

    p2 = _node_json_path(2)
    p2n = _node_json_path(2, c_pf=C_IN_F_HUMAN_D2_PF)
    assert _node_json_path(2, c_pf=C_IN_F_HUMAN_PF) == p2
    assert _node_json_path(None, c_pf=C_IN_F_HUMAN_PF) == _node_json_path(None) == NODE_JSON
    assert p2n.endswith("h9p-node-115.5209-d2.json") and os.path.dirname(p2n) == JSON_DIR
    assert p2n != p2 and p2n != NODE_JSON
    with pytest.raises(ValueError):
        _node_json_path(None, c_pf=C_IN_F_HUMAN_D2_PF)
    print(
        f"[WF-10 step4c N0] knob unset -> {C_IN_F_HUMAN_PF} pF; d2root + degree 2 -> "
        f"{C_IN_F_HUMAN_D2_PF} pF (secant {C_IN_SECANT_D2_PF:.6f}); paths: {p2} / {p2n} "
        f"(ASSERTED)"
    )

    # (N2) the closed form Z_RC = r - j/(2 pi f C), and the asserted separation.
    im = {}
    for c_pf in (C_IN_F_HUMAN_PF, C_IN_F_HUMAN_D2_PF):
        c_f = c_pf * 1e-12
        z = complex(_z_rc(c_f, FREQUENCY_HZ))
        ref = ESR_OHM - 1j / (2.0 * np.pi * float(FREQUENCY_HZ) * c_f)
        rel = abs(z - ref) / abs(ref)
        im[c_pf] = z.imag
        print(
            f"[WF-10 step4c N2] C {c_pf} pF: Z_RC {z:.12e} vs closed form {ref:.12e}, rel "
            f"{rel:.3e} (band 1e-12, ASSERTED)"
        )
        assert rel <= 1e-12, (c_pf, rel)
    ratio = im[C_IN_F_HUMAN_D2_PF] / im[C_IN_F_HUMAN_PF]
    want = C_IN_F_HUMAN_PF / C_IN_F_HUMAN_D2_PF
    print(
        f"[WF-10 step4c N2] NEG separation Im Z_RC(d2) / Im Z_RC(d1) = {ratio:.6f} vs "
        f"118.2437/115.5209 = {want:.6f} (ASSERTED to 1e-12 rel; a knob-blind builder gives 1)"
    )
    assert abs(ratio / want - 1.0) <= 1e-12, (ratio, want)
    assert abs(ratio - 1.0) > 1e-2, ratio


@complex_only
def test_wf10_step4c_build_only():
    """(K1n) step 4c: the builder sees ``FEM_EM_WF10_NODE=d2root`` through
    ``mpiexec`` at degree 2; build-only, no node record written.  Reads the
    window's command-line environment (no monkeypatch)."""
    if os.environ.get(STEP4C_ENV, "").strip() != "1":
        pytest.skip(f"{STEP4C_ENV} unset: build-only -n 2 F-human window, environment-selected")
    from tests.validation.test_wf12_f_human_hole import _rss_gib

    comm = MPI.COMM_WORLD

    def say(msg):
        if comm.rank == 0:
            print(f"[WF-10 step4c K1n] {msg}", flush=True)

    target = os.path.join(JSON_DIR, f"h9p-node-{C_IN_F_HUMAN_D2_PF:.4f}-d2.json")
    existed_before = os.path.exists(target)
    t0 = time.perf_counter()
    out = build_tuned_quadrature(comm, with_z50=False, build_only=True)
    t_build = time.perf_counter() - t0
    n2, e, f = out["n_dofs"], out["n_edges"], out["n_faces"]
    pfs = comm.allgather(_node_pf())
    exists_after = os.path.exists(out["node_json"])
    _m, _per, rss_sum = _rss_gib(comm)
    say(
        f"(Deg-2) degree {out['degree']}; unknowns {n2} vs 2E+2F = {2 * e + 2 * f} (record "
        f"{N_DOFS_D2_RECORD}; ASSERTED exact); cells {out['cells']}"
    )
    say(
        f"node: c_in_pf {out['c_in_pf']} (want {C_IN_F_HUMAN_D2_PF}, ASSERTED exact; a knob-blind "
        f"builder returns {C_IN_F_HUMAN_PF}); node_json {os.path.basename(out['node_json'])} "
        f"(ASSERTED); per-rank _node_pf {pfs}; file existed before {existed_before}, after "
        f"{exists_after} (absent ASSERTED)"
    )
    say(
        f"build {t_build:.1f} s; RSS summed {rss_sum:.2f} GiB (PREDICTED <= "
        f"{PREDICTED_RSS_D2_GIB:g}; H26a W1 110.0 s / {RSS_D2_BUILD_RECORD_GIB} GiB; printed)"
    )
    assert out["degree"] == 2
    assert n2 == 2 * e + 2 * f == N_DOFS_D2_RECORD, (n2, e, f)
    assert out["c_in_pf"] == C_IN_F_HUMAN_D2_PF, out["c_in_pf"]
    assert out["c_in_pf"] != C_IN_F_HUMAN_PF
    assert os.path.basename(out["node_json"]) == "h9p-node-115.5209-d2.json", out["node_json"]
    assert out["node_json"] != _node_json_path(2)
    assert len(set(pfs)) == 1 and pfs[0] == C_IN_F_HUMAN_D2_PF, pfs
    assert not exists_after, f"{out['node_json']} exists after a build-only call"
