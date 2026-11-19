"""`WF-10` step 1e part 1 (§10 chain H17 part 1, §9 item 334) — a C-node module on
the copper F-human hole mesh with series-RC sheets (r = 0.2 Ω) on all 32 gaps
and the two-port 90° feed, reading the sign of Im ``Y_in`` per node.

Brief: ``docs/planning/chunks/WF-10.md`` "§9 item 334" and its 2026-10-08 18:00
review rulings (which outrank the brief): drive sense never asserted; the
stored ``h12c-cu-s32-a2.npz`` is the (Z50) comparand (absent ⇒ stop, never
re-solve); the smoke node is copper's lossless ``C_mode`` 118.1186 pF; the
"dropping r moves Y_in" control is **predicted** (printed, never asserted);
the RC-node power partition is the anchor, asserted at
``DISCRETE_IDENTITY_RTOL``.

Two tests:

* ``test_wf10_step1e_y_in_coding`` (Y0) — numpy only, always runs: the S-route
  ``Y_in`` (feeds' 2×2 S at 50 Ω with the gap RC inside, then Γ_act) equals the
  current-route ``1/Z_RC + I_f/V_f`` within 1e-12 on a seeded random passive Z.
* ``test_wf10_step1e_copper_esr_nodes`` — selected by ``FEM_EM_WF10_STEP1E=1``
  (never ``-k``); heavy, ``-n 8``.  Node list from
  ``FEM_EM_WF10_STEP1E_NODES_PF`` (comma-separated pF; unset ⇒ the npz
  ``c_mode``).

``Y_in`` (defined here; §10 does not): the active input admittance at P17
under incident waves a = (1, +1j) at 50 Ω: Γ_act = S₁₁ + S₁₂·(+1j),
``Y_in = (1 − Γ_act)/(50(1 + Γ_act))`` — the source's view of
(network ∥ the feed gap's r + C), i.e. ``ring_mode_rows``' ``y_in_with_c_s``
(`src/fem_em_solver/ports/ring_modes.py:163`) with 1/Z_RC for jωC.  The sign of
Im ``Y_in`` is a reading, never asserted.

Scope: degree 1, 64 MHz, κ off, one copper fixture, one C node per pair.  Not a
root, not C_in, not tuning / f₀ / Q / match, not the circuit ↔ field identity,
not |B₁⁺| / hotspot, not degree 2, not AED.
"""

from __future__ import annotations

import json
import os
import time

import numpy as np
import pytest
from mpi4py import MPI

from tests.complex_mode import complex_only

STEP1E_ENV = "FEM_EM_WF10_STEP1E"
NODES_ENV = "FEM_EM_WF10_STEP1E_NODES_PF"
DRIVEN = "P17"
NPZ_PATH = "/workspace/logs/wf12/h12c-cu-s32-a2.npz"
JSON_DIR = "/workspace/logs/wf10"
#: Capacitor ESR (§9 binding ruling "capacitor ESR"): r = 0.2 Ω on every gap.
ESR_OHM = 0.2
#: (Y0) coding band (brief: 1e-12).
Y0_CODING_BAND = 1.0e-12
#: Ruling (4): predicted, printed only — dropping r moves Y_in by ≥ 1e-6 rel.
PREDICTED_DROP_R_MOVE_FLOOR = 1.0e-6
#: The smoke node (ruling (3)): copper's lossless C_mode, the npz value.
C_MODE_CU_PF = 118.1186
# Predictions written first (§9 item 334), printed only.
PREDICTED_Z50 = 1.0e-10
PREDICTED_P_RC = 1.0e-10
PREDICTED_DZ_ORDER = 1.0e-1
Q_C_PREDICTED = 105.0
LOSSLESS_K1_SHARE = 0.99886

# ---- §9 item 459 (§10 H26a): the degree knob on step 1e's root module ------------
# The knob is H20a's (`test_wf10_f_human_tuned_b1_sar.DEGREE_ENV` /
# `_degree_knob` / `_apply_degree`), imported inside function bodies only:
# that module imports this one at module level, so a module-level import here
# would be circular.
#: Build-only selector (W1 of item 459 only; never set in `xl` entry 39).
BUILD_ENV = "FEM_EM_WF10_STEP1E_BUILD"
#: At degree 2: the degree-1 record comparands (printed only) and the gates that
#: stay asserted.
STEP1E_DEGREE2_GATES = {
    "printed": ("(R)", "(Z50)"),
    "asserted": ("(A0)", "(iv)", "(Deg)", "kinds", "(P50)", "NEG", "(P_RC)"),
}
#: Entry 32's degree-2 same-comparison readings (P17 driven, all 32 sheets at
#: 50 Ω, copper hole; `20261013T070008Z_WF-12-step4-H16.log:6711–6712`), printed
#: beside (P50) / NEG at degree 2, never a band.
ENTRY32_P_D2 = 1.468e-13
ENTRY32_NEG_D2 = 119.0


def _h17_node_json_path(c_pf, degree):
    """The step 1e node JSON: the degree-1 name unchanged (``degree is None``),
    a ``-d2`` sibling at degree 2 so a degree-2 run never overwrites a degree-1
    record; anything else raises."""
    if degree is None:
        return os.path.join(JSON_DIR, f"h17-node-{c_pf:.4f}.json")
    if degree == 2:
        return os.path.join(JSON_DIR, f"h17-node-{c_pf:.4f}-d2.json")
    raise ValueError(degree)


def _build_root_fixture(comm):
    """Step 1e's fixture: build, (A0) / (iv), then H20a's knob (``_apply_degree``,
    (Deg) asserted inside).  Compiles no form and assembles nothing.  At degree 2
    ``solver_built`` is ``None`` so two solvers are never held."""
    from tests.mesh.test_birdcage_port_sheet_prerequisite import CELL_COUNT_BAND
    from tests.validation.test_birdcage_f_human_rung import _ring_ports
    from tests.validation.test_wf7_f_human_b1_quadrature import FREQUENCY_HZ, _build_f_human
    from tests.validation.test_wf10_f_human_tuned_b1_sar import _apply_degree
    from tests.validation.test_wf12_f_human_copper_matrix import HOLE_CELL_RECORD
    from tests.validation.test_wf12_f_human_hole import _census_checks, _owned_census

    def say(msg):
        if comm.rank == 0:
            print(f"[WF-10 step1e] {msg}", flush=True)

    t0 = time.perf_counter()
    built = _build_f_human(comm, as_hole=True)
    t_build = time.perf_counter() - t0
    ctx = built["ctx"]
    msh, omega, cell_tags = ctx["msh"], ctx["omega"], ctx["cell_tags"]
    assert abs(omega / (2.0 * np.pi * float(FREQUENCY_HZ)) - 1.0) < 1e-12, omega
    census = _owned_census(msh, cell_tags, comm)
    census_fails = _census_checks(census, built["cells"], _ring_ports())
    cell_rel = built["cells"] / HOLE_CELL_RECORD - 1.0
    say(
        f"(A0) hole cells {built['cells']} vs record {HOLE_CELL_RECORD}, rel {cell_rel:.3e} "
        f"(band {CELL_COUNT_BAND}, ASSERTED); build {t_build:.1f} s; "
        f"(iv) failures {census_fails or 'none'}"
    )
    assert abs(cell_rel) <= CELL_COUNT_BAND, cell_rel
    assert not census_fails, census_fails
    solver_built = ctx["solver"]
    _solver, degree, counts = _apply_degree(ctx, comm)
    del _solver
    if degree != 1:
        solver_built = None  # never hold two solvers (`WF-12` step 4's trap)
    return {
        "built": built,
        "ctx": ctx,
        "cells": built["cells"],
        "t_build": t_build,
        "degree": degree,
        "counts": counts,
        "solver_built": solver_built,
    }


# ============================================================================
# Lifted verbatim from `attempt/WF-12-step3-a2-20261006T1702Z` (tip b3cafbfd),
# tests/validation/test_wf12_f_human_copper_c_mode.py lines 130–166
# (`_thevenin_feed_impedance` … `_parallel_loads`).  Copied, never imported
# from the branch.  ``z_c`` may be any complex gap impedance (here Z_RC).
# ============================================================================


def _thevenin_feed_impedance(z_c, z0):
    """Z_th of a z0 source in parallel with the gap's Z_C."""
    return z0 * z_c / (z0 + z_c)


def _source_voltage(z_c, z0):
    """The z0 source's open-circuit V_s for V_th = 1 V across Z_C."""
    return (z0 + z_c) / z_c


def _current_route_s(currents_by_drive, feeds_idx, z_th, v_s):
    """``S`` (len(feeds)²) from port currents: column j from the drive of feed j.

    ``currents_by_drive[j]`` is the full current vector with feed j driven at
    V_th = 1 V through Z_th; V_i = δ_ij − Z_th I_i; S_ij = 2V_i/V_s − δ_ij.
    """
    n = len(feeds_idx)
    s = np.zeros((n, n), dtype=np.complex128)
    for j in range(n):
        cur = currents_by_drive[j]
        for i, idx in enumerate(feeds_idx):
            d = 1.0 if i == j else 0.0
            v = d - z_th * cur[idx]
            s[i, j] = 2.0 * v / v_s - d
    return s


def _parallel_currents(z32, z_load, drive):
    """``(Z + diag Z_load)⁻¹ v`` — circuit-side port currents for Thevenin ``v``."""
    return np.linalg.solve(
        np.asarray(z32) + np.diag(z_load), np.asarray(drive, dtype=np.complex128)
    )


def _parallel_loads(n, feeds_idx, z_c, z_feed):
    z_load = np.full(n, z_c, dtype=np.complex128)
    for i in feeds_idx:
        z_load[i] = z_feed
    return z_load


# ============================================================================
# End of the lifted helper block.
# ============================================================================


def _z_rc(c_f, frequency_hz=64.0e6, r_ohm=ESR_OHM):
    from fem_em_solver.ports.lumped import series_rlc_impedance

    return series_rlc_impedance(float(frequency_hz), r_ohm=float(r_ohm), l_h=0.0, c_f=float(c_f))


def _y_from_gamma(gamma, z0):
    return (1.0 - gamma) / (z0 * (1.0 + gamma))


def _active_y_in(s2, a, z0):
    """``[Y_in at feed 0, Y_in at feed 1]`` under incident waves ``a`` (2-vector)."""
    a = np.asarray(a, dtype=np.complex128)
    b = np.asarray(s2, dtype=np.complex128) @ a
    return [complex(_y_from_gamma(b[i] / a[i], z0)) for i in range(2)]


def _s_route_2x2(s_full, feeds_idx, z_load_others, z_shunt_feed, z0):
    """Feeds' 2×2 S at z0 with every other port terminated in ``z_load_others``
    (S-domain, `reduce_terminated_ports`) and the feed gap's ``z_shunt_feed``
    across each feed (Y-domain shunt)."""
    from fem_em_solver.ports.circuit import reduce_terminated_ports, s_to_z, z_to_s

    n = s_full.shape[0]
    assert list(feeds_idx) == [0, 1]
    s_red = reduce_terminated_ports(s_full, z0, {k: z_load_others for k in range(2, n)})
    y2 = np.linalg.inv(s_to_z(s_red, z0)) + np.eye(2) / z_shunt_feed
    return z_to_s(np.linalg.inv(y2), z0)


def _current_route_y_in(z_full, feeds_idx, z_load_others, z_shunt_feed, a, z0):
    """``1/Z_shunt + I_f/V_f`` at each feed, the drive's Thevenin vector ∝ a."""
    n = z_full.shape[0]
    z_th = _thevenin_feed_impedance(z_shunt_feed, z0)
    z_load = _parallel_loads(n, feeds_idx, z_load_others, z_th)
    v_th = np.zeros(n, dtype=np.complex128)
    for i, w in zip(feeds_idx, a):
        v_th[i] = w
    cur = _parallel_currents(z_full, z_load, v_th)
    out = []
    for i in feeds_idx:
        v_port = v_th[i] - z_th * cur[i]
        out.append(complex(1.0 / z_shunt_feed + cur[i] / v_port))
    return out


def _circuit_current_route_s(z_full, feeds_idx, z_load_others, z_shunt_feed, z0):
    """The field route's exact circuit mirror: unit Thevenin drives through Z_th,
    then `_current_route_s` with V_s."""
    n = z_full.shape[0]
    z_th = _thevenin_feed_impedance(z_shunt_feed, z0)
    v_s = _source_voltage(z_shunt_feed, z0)
    z_load = _parallel_loads(n, feeds_idx, z_load_others, z_th)
    curs = []
    for j in feeds_idx:
        e = np.zeros(n, dtype=np.complex128)
        e[j] = 1.0
        curs.append(_parallel_currents(z_full, z_load, e))
    return _current_route_s(curs, feeds_idx, z_th, v_s)


def _rel(a, b):
    a, b = np.asarray(a, dtype=np.complex128), np.asarray(b, dtype=np.complex128)
    return float(np.max(np.abs(a - b)) / np.max(np.abs(b)))


SENSES = {"(1,+1j)": (1.0, 1.0j), "(1,-1j)": (1.0, -1.0j)}


def test_wf10_step1e_y_in_coding():
    """(Y0): S-route Y_in = current-route Y_RC + I_f/V_f within 1e-12."""
    from fem_em_solver.ports.circuit import z_to_s

    z0 = 50.0
    f = 64.0e6
    rng = np.random.default_rng(20261008)
    n = 32
    b = rng.standard_normal((n, n))
    r_part = 5.0 * (b @ b.T) / n + 0.5 * np.eye(n)  # Re Z symmetric positive definite
    x = rng.standard_normal((n, n))
    x_part = 40.0 * (x + x.T)
    z = r_part + 1j * x_part
    assert np.min(np.linalg.eigvalsh(0.5 * (z + z.conj().T))) > 0.0  # passive
    s = z_to_s(z, z0)
    c_f = C_MODE_CU_PF * 1e-12
    z_rc = _z_rc(c_f, f)
    z_c = _z_rc(c_f, f, r_ohm=0.0)
    feeds_idx = [0, 1]

    s_rt = _s_route_2x2(s, feeds_idx, z_rc, z_rc, z0)
    s_cm = _circuit_current_route_s(z, feeds_idx, z_rc, z_rc, z0)
    r_s = _rel(s_cm, s_rt)
    print(
        f"[WF-10 step1e Y0] Z_RC = {z_rc:.9e} Ohm (r {ESR_OHM}, C {C_MODE_CU_PF} pF); "
        f"S-route vs current-route-S {r_s:.3e} (band {Y0_CODING_BAND:g}, ASSERTED)"
    )
    worst = 0.0
    for label, a in SENSES.items():
        y_s = _active_y_in(s_rt, a, z0)
        y_c = _current_route_y_in(z, feeds_idx, z_rc, z_rc, a, z0)
        y_m = _active_y_in(s_cm, a, z0)
        r = max(_rel(y_s, y_c), _rel(y_m, y_c))
        worst = max(worst, r)
        print(
            f"[WF-10 step1e Y0] a = {label}: Y_in S-route {y_s[0]:.9e} / {y_s[1]:.9e} S; "
            f"current route {y_c[0]:.9e} / {y_c[1]:.9e} S; rel {r:.3e}"
        )
    # Ruling (4): predicted, printed only — dropping r on the 30 loads only.
    a = SENSES["(1,+1j)"]
    y_rc = _current_route_y_in(z, feeds_idx, z_rc, z_rc, a, z0)
    y_c0 = _current_route_y_in(z, feeds_idx, z_c, z_rc, a, z0)
    move = abs(y_c0[0] - y_rc[0]) / abs(y_rc[0])
    per_load = abs(1.0 - 1.0 / (1.0 + 1j * 2 * np.pi * f * c_f * ESR_OHM))
    print(
        f"[WF-10 step1e Y0] drop-r control (PREDICTED >= {PREDICTED_DROP_R_MOVE_FLOOR:g}, "
        f"printed, never asserted): measured rel move {move:.3e}; closed-form per-load "
        f"|1 - 1/(1 + jwCr)| = {per_load:.3e}"
    )
    assert r_s <= Y0_CODING_BAND, f"(Y0) S routes differ {r_s:.3e}"
    assert worst <= Y0_CODING_BAND, f"(Y0) Y_in routes differ {worst:.3e}"


@complex_only
def test_wf10_step1e_copper_esr_nodes():
    if os.environ.get(STEP1E_ENV, "").strip() != "1":
        pytest.skip(f"{STEP1E_ENV} unset: heavy -n 8 F-human window, environment-selected")

    import ufl
    from dolfinx import fem

    from fem_em_solver.core.cavity import surface_resistance_ohm
    from fem_em_solver.ports.circuit import s_to_z
    from fem_em_solver.ports.lumped import (
        lumped_port_bilinear_term,
        lumped_port_linear_term,
        sheet_terminal_current,
    )
    from tests.mesh.test_birdcage_port_sheet_prerequisite import CELL_COUNT_BAND
    from tests.validation.test_birdcage_power_identity import (
        DISCRETE_IDENTITY_RTOL,
        _sheet_field_dissipation_w,
        _source_power_w,
    )
    from tests.validation.test_port_birdcage_four_port import TERMINATED_PORT_IMPEDANCE_OHM
    from tests.validation.test_port_birdcage_lumped_column import PHANTOM_CELL_TAG
    from tests.validation.test_port_birdcage_ring_matrix import REUSE_REPRODUCTION_RTOL
    from tests.validation.test_port_circuit_layer_field import _residual
    from tests.validation.test_port_drive_superposition import _ring_quadrature_slots
    from tests.validation.test_port_package_sparameters import REFERENCE_IMPEDANCE_OHM
    from tests.validation.test_th14_birdcage_copper import _leontovich_term, _surface_loss_w
    from tests.validation.test_wf7_f_human_b1_quadrature import DEGREE, FREQUENCY_HZ
    from tests.validation.test_wf10_f_human_tuned_circuit import (
        QUADRATURE_SLOT_OFFSET,
        _k1_shares,
        _respec,
    )
    from tests.validation.test_wf12_f_human_copper_matrix import (
        BLIND_FACTOR,
        COPPER_SIGMA,
        HOLE_CELL_RECORD,
        P17_P_SRC_RECORD_W,
        P17_P_SURF_RECORD_W,
    )
    from tests.validation.test_wf12_f_human_hole import _rss_gib

    comm = MPI.COMM_WORLD
    z0 = float(REFERENCE_IMPEDANCE_OHM)
    t_window = time.perf_counter()

    def say(msg):
        if comm.rank == 0:
            print(f"[WF-10 step1e] {msg}", flush=True)

    # ---- the stored comparand (ruling (2): absent ⇒ stop, never re-solve) ----
    assert os.path.exists(
        NPZ_PATH
    ), f"{NPZ_PATH} absent: item 334 stops (BLOCKED, rule (d)); never re-solve"
    with np.load(NPZ_PATH, allow_pickle=False) as npz:
        s32_npz = np.array(npz["s32"], dtype=np.complex128)
        order_npz = [str(p) for p in npz["order"]]
        cells_npz = int(npz["cells"])
        c_mode_npz = float(npz["c_mode"])
        meas50_npz = np.array(npz["meas50"], dtype=np.complex128)
    nodes_env = os.environ.get(NODES_ENV, "").strip()
    nodes_f = (
        [float(x) * 1e-12 for x in nodes_env.split(",") if x.strip()] if nodes_env else [c_mode_npz]
    )
    from tests.validation.test_wf10_f_human_tuned_b1_sar import DEGREE_ENV, _degree_knob

    knob = _degree_knob()
    d2 = knob is not None
    say(
        f"start: -n {comm.size}, f = {FREQUENCY_HZ:.3e} Hz, degree "
        f"{DEGREE if knob is None else knob} ({DEGREE_ENV} {knob if d2 else 'unset'}; fixture "
        f"DEGREE {DEGREE}), kappa off, "
        f"r = {ESR_OHM} Ohm on all 32 gaps; npz c_mode {c_mode_npz:.12e} F; nodes "
        f"{', '.join(f'{c * 1e12:.4f}' for c in nodes_f)} pF"
    )

    # ---- fixture, (A0), (iv), (Deg) (item 459: `_build_root_fixture`) ---------
    fixture = _build_root_fixture(comm)
    built, t_build = fixture["built"], fixture["t_build"]
    assert fixture["degree"] == (knob or 1), (fixture["degree"], knob)
    del fixture
    ctx = built["ctx"]
    msh, tags_f, omega, cell_tags = ctx["msh"], ctx["tags_f"], ctx["omega"], ctx["cell_tags"]
    say(
        f"(A0) npz cells {cells_npz} (ASSERTED equal); census {CELL_COUNT_BAND}-banded above; "
        f"hole record {HOLE_CELL_RECORD}"
    )
    assert built["cells"] == cells_npz, (built["cells"], cells_npz)

    z_s = (1.0 + 1.0j) * surface_resistance_ohm(omega, COPPER_SIGMA)
    wall = _leontovich_term(msh, built["facets"], z_s, omega)
    z_p = complex(TERMINATED_PORT_IMPEDANCE_OHM)
    dx = ufl.Measure("dx", domain=msh, subdomain_data=cell_tags)
    acct = {"mesh": msh, "facet_tags": tags_f}
    solver = ctx["solver"]
    base_specs = list(ctx["specs"])

    # Lifted from `attempt/WF-12-step3-a2-20261006T1702Z` (tip b3cafbfd),
    # test_wf12_f_human_copper_c_mode.py lines 314–365 (the `_drive` closure),
    # verbatim except: the κ branch (`keep_e`) dropped (κ is never computed
    # here), and the partition also keeps the per-sheet dissipation by id (for
    # the printed Dz control).  ``s`` uses a 50 Ω ``z_p``: valid only where
    # every sheet is 50 Ω — at a node only ``currents`` is read.
    def _drive(specs, driven_id, *, held, partition):
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
        if partition:
            sigma_f = fields.sigma_field

            # F821: pyflakes sees the `del` at the end of this function; the closure runs before it.
            def _vol(measure):
                form = fem.form(0.5 * sigma_f * ufl.inner(e, e) * measure)  # noqa: F821
                return float(np.real(comm.allreduce(fem.assemble_scalar(form), op=MPI.SUM)))

            by_id = {
                s.port_id: float(_sheet_field_dissipation_w(acct, s, e, omega)) for s in sheets
            }
            out["p_sheet_by_id"] = by_id
            out["p_sheets"] = float(sum(by_id.values()))
            out["p_src"] = float(_source_power_w(acct, dsheet, e, omega))
            p_total = _vol(dx)
            out["p_phantom"] = _vol(dx(PHANTOM_CELL_TAG))
            out["p_rest"] = p_total - out["p_phantom"]
            out["p_surf"] = _surface_loss_w(msh, built["facets"], z_s, e, comm)
        del fields, e
        out["t_post"] = time.perf_counter() - t1
        return out

    # ---- order (item 285's: [P17, feed2] + rest by number) --------------------
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
    say(
        f"feeds: {DRIVEN} (slot {ref_k}, ring {ref_sign:+.0f}) and {feed2} (slot "
        f"{slots[feed2][0]}); npz order[:2] {order_npz[:2]}"
    )
    assert order_npz[:2] == [DRIVEN, feed2], order_npz[:2]
    assert order_npz == order, "npz order differs from this run's order"

    # Lifted `_pair` (branch lines 471–482), verbatim except the P17 drive's
    # partition switch (True here: (R)/(P50) on the Z50 pair, (P_RC) at a node).
    def _pair(imp_x, z_feed, v_feed):
        specs_x = _respec(base_specs, None, imp_x)
        assert all(sp.width_correction_kappa is None for sp in specs_x)
        cur, knd, tt, parts = [], [], [], []
        for n_drive, pid in enumerate(feeds):
            d = _drive(specs_x, pid, held=(n_drive > 0), partition=(n_drive == 0))
            knd.append(d["kind"])
            tt.append(d["t_solve"] + d["t_post"])
            cur.append(np.array([d["currents"][p] for p in order], dtype=np.complex128))
            if n_drive == 0:
                parts.append({k: v for k, v in d.items() if k.startswith("p_")})
            del d
        return _current_route_s(cur, feeds_idx, z_feed, v_feed), cur, knd, tt, parts[0]

    def _show(name, m):
        say(f"    {name} " + "; ".join(", ".join(f"{complex(v):.10e}" for v in row) for row in m))

    def _partition(p):
        p_vol = p["p_phantom"] + p["p_rest"]
        rel = abs(p["p_src"] - p["p_sheets"] - p_vol - p["p_surf"]) / abs(p["p_src"])
        blind = abs(p["p_src"] - p["p_sheets"] - p_vol) / abs(p["p_src"])
        return rel, blind

    # ---- the Z50 pair: (R), (P50), NEG, (Z50) ---------------------------------
    s_50, _c50, kinds_50, t_50, p50 = _pair({pid: z0 for pid in order}, z0, 1.0)
    del _c50
    rel_p50, blind50 = _partition(p50)
    r_src = abs(p50["p_src"] / P17_P_SRC_RECORD_W - 1.0)
    r_surf = abs(p50["p_surf"] / P17_P_SURF_RECORD_W - 1.0)
    say(f"Z50 pair {kinds_50}, {', '.join(f'{t:.2f}' for t in t_50)} s")
    d1_only = "degree-1 record comparand, d1 -> d2 move, printed only (H26a)" if d2 else "ASSERTED"
    say(
        f"(R) P_src {p50['p_src']:.9e} vs {P17_P_SRC_RECORD_W:.9e} rel {r_src:.3e}; P_surf "
        f"{p50['p_surf']:.9e} vs {P17_P_SURF_RECORD_W:.9e} rel {r_surf:.3e} (band "
        f"{REUSE_REPRODUCTION_RTOL:g}, {d1_only}; measured 2.109e-11 / 1.506e-11)"
    )
    say(
        f"(P50) P_src {p50['p_src']:.9e} sumP_sheet {p50['p_sheets']:.9e} P_phantom "
        f"{p50['p_phantom']:.9e} P_rest {p50['p_rest']:.9e} P_surf {p50['p_surf']:.9e} W; rel "
        f"{rel_p50:.3e} (band {DISCRETE_IDENTITY_RTOL:g}, ASSERTED; measured 2.420e-13"
        + (f"; entry 32 degree 2 {ENTRY32_P_D2:.3e}, beside)" if d2 else ")")
    )
    say(
        f"NEG P_surf-blind miss {blind50:.4e} = {blind50 / DISCRETE_IDENTITY_RTOL:.1f}x the band "
        f"(bar {BLIND_FACTOR:.0f}x, ASSERTED; measured 123.6x / 123.3x"
        + (f"; entry 32 degree 2 {ENTRY32_NEG_D2:.1f}x, beside)" if d2 else ")")
    )
    s32_ff = s32_npz[:2, :2]
    r_z50 = _residual(s_50, s32_ff)
    r_z50_meas = _residual(s_50, meas50_npz)
    _show("S32npz[:2,:2]", s32_ff)
    _show("MEAS50       ", s_50)
    say(
        f"(Z50) ||MEAS50 - S32npz[:2,:2]||/||S32npz[:2,:2]|| = {r_z50:.6e} (band "
        f"{REUSE_REPRODUCTION_RTOL:g}, {d1_only}; PREDICTED <= {PREDICTED_Z50:g}); vs npz meas50 "
        f"{r_z50_meas:.6e} (printed); t={time.perf_counter() - t_window:.1f} s"
    )

    # ---- per node --------------------------------------------------------------
    z32_npz = s_to_z(s32_npz, z0)
    node_results = []
    node_pair_s = []
    for c_f in nodes_f:
        z_rc = _z_rc(c_f, FREQUENCY_HZ)
        z_c = _z_rc(c_f, FREQUENCY_HZ, r_ohm=0.0)
        z_th = _thevenin_feed_impedance(z_rc, z0)
        v_s = _source_voltage(z_rc, z0)
        say(
            f"node C = {c_f:.12e} F ({c_f * 1e12:.4f} pF): Z_RC {z_rc:.9e} Ohm, Z_th "
            f"{z_th:.9e} Ohm, V_s {v_s:.9e} V"
        )
        imp = {pid: z_rc for pid in others}
        imp.update({pid: z_th for pid in feeds})
        s_node, cur_node, kinds_n, t_n, pn = _pair(imp, z_th, v_s)
        _show("S_node", s_node)
        rel_prc, blind_n = _partition(pn)
        p_rc = sum(v for pid, v in pn["p_sheet_by_id"].items() if pid not in feeds)
        p_feeds = sum(v for pid, v in pn["p_sheet_by_id"].items() if pid in feeds)
        p_vol = pn["p_phantom"] + pn["p_rest"]
        dz = abs(pn["p_src"] - (pn["p_sheets"] - p_rc) - p_vol - pn["p_surf"]) / abs(pn["p_src"])
        say(
            f"(P_RC) P_src {pn['p_src']:.9e} sumP_sheet {pn['p_sheets']:.9e} (30 RC "
            f"{p_rc:.9e}, 2 feeds {p_feeds:.9e}) P_phantom {pn['p_phantom']:.9e} P_rest "
            f"{pn['p_rest']:.9e} P_surf {pn['p_surf']:.9e} W; rel {rel_prc:.3e} (band "
            f"{DISCRETE_IDENTITY_RTOL:g}, ASSERTED; PREDICTED <= {PREDICTED_P_RC:g})"
        )
        say(
            f"(Dz) partition without the 30 RC sheets misses {dz:.4e} (PREDICTED order "
            f"{PREDICTED_DZ_ORDER:g}, printed); P_surf-blind miss at the node {blind_n:.4e} "
            f"(printed)"
        )
        q_frac = p_rc / (p_rc + pn["p_surf"] + pn["p_phantom"])
        say(
            f"P_RC/(P_RC + P_surf + P_phantom) = {q_frac:.6e} (beside Q_C = 1/(wCr) = "
            f"{1.0 / (omega * c_f * ESR_OHM):.2f}, predicted ≈ {Q_C_PREDICTED:g}; printed)"
        )
        y_in = {}
        for label, a in SENSES.items():
            y_in[label] = _active_y_in(s_node, a, z0)
            y17, y21 = y_in[label]
            spread = abs(y17 - y21) / abs(y17)
            say(
                f"Y_in a = {label}: {DRIVEN} {y17.real:+.9e} {y17.imag:+.9e}j S (sign Im "
                f"{'+' if y17.imag > 0 else '-' if y17.imag < 0 else '0'}); {feed2} "
                f"{y21.real:+.9e} {y21.imag:+.9e}j S (sign Im "
                f"{'+' if y21.imag > 0 else '-' if y21.imag < 0 else '0'}); C4/mirror spread "
                f"{spread:.4e} (READING, never asserted)"
            )
        for lbl, z_load in (("Z_RC (r=0.2)", z_rc), ("Z_C (r=0)", z_c)):
            yc = _current_route_y_in(z32_npz, feeds_idx, z_load, z_load, SENSES["(1,+1j)"], z0)
            say(
                f"circuit-only ({'degree-1 npz' if d2 else 'npz'} s32, local loads {lbl}; not a prediction of the in-model "
                f"sign — H12c measured the circuit <-> field identity at 1.58 here): Y_in "
                f"(1,+1j) {DRIVEN} {yc[0]:.9e} S, {feed2} {yc[1]:.9e} S"
            )
        sh = _k1_shares(cur_node[0], order, slots)
        say(
            f"k=+-1 share, P17 column at the node: feed ring {sh[ref_sign]:.6f} (ring+ "
            f"{sh[1.0]:.6f} ring- {sh[-1.0]:.6f}) beside lossless {LOSSLESS_K1_SHARE}"
        )
        say(
            f"node pair {kinds_n}, {', '.join(f'{t:.2f}' for t in t_n)} s; "
            f"t={time.perf_counter() - t_window:.1f} s"
        )
        rec = {
            "c_f": c_f,
            "r_ohm": ESR_OHM,
            "freq_hz": float(FREQUENCY_HZ),
            "order": order,
            "degree": knob or 1,
            "s_node": [[[complex(v).real, complex(v).imag] for v in row] for row in s_node],
            "y_in": {lbl: [[y.real, y.imag] for y in ys] for lbl, ys in y_in.items()},
            "partition": {k: v for k, v in pn.items()},
            "rel_p_rc": rel_prc,
            "dz": dz,
            "k1_share_feed_ring": sh[ref_sign],
            "times_s": {"pair": t_n, "z50_pair": t_50, "build": t_build},
        }
        if comm.rank == 0:
            os.makedirs(JSON_DIR, exist_ok=True)
            path = _h17_node_json_path(c_f * 1e12, knob)
            with open(path, "w") as fh:
                json.dump(rec, fh, indent=1)
            print(f"[WF-10 step1e] wrote {path}", flush=True)
        comm.Barrier()
        node_results.append((c_f, kinds_n, rel_prc))
        node_pair_s.append(sum(t_n))

    _m, rss_per, rss_sum = _rss_gib(comm)
    say(
        f"price: build {t_build:.1f} s, Z50 pair {sum(t_50):.1f} s, node pairs "
        f"{', '.join(f'{t:.1f}' for t in node_pair_s)} s ({len(node_results)} node(s)); RSS per rank "
        f"{', '.join(f'{r:.2f}' for r in (rss_per or []))} GiB, summed {rss_sum:.2f} GiB; "
        f"window elapsed {time.perf_counter() - t_window:.1f} s at -n {comm.size}"
    )

    # ---- asserts: anchors first -----------------------------------------------
    assert kinds_50 == ["solve", "held"], kinds_50
    if not d2:  # (R): degree-1 records, printed only at degree 2 (H26a)
        assert r_src <= REUSE_REPRODUCTION_RTOL, r_src
        assert r_surf <= REUSE_REPRODUCTION_RTOL, r_surf
    assert rel_p50 <= DISCRETE_IDENTITY_RTOL, rel_p50
    assert blind50 >= BLIND_FACTOR * DISCRETE_IDENTITY_RTOL, blind50
    if not d2:  # (Z50): the degree-1 npz, printed only at degree 2 (H26a)
        assert r_z50 <= REUSE_REPRODUCTION_RTOL, f"(Z50) closure {r_z50:.6e}"
    for c_f, kinds_n, rel_prc in node_results:
        assert kinds_n == ["solve", "held"], kinds_n
        assert (
            rel_prc <= DISCRETE_IDENTITY_RTOL
        ), f"(P_RC) at {c_f * 1e12:.4f} pF: {rel_prc:.3e} — never re-band"


# ============================================================================
# `WF-10` step 1e part 2 (§10 H17 part 2, §9 item 354): the root reader.
# Brief: docs/planning/chunks/WF-10.md "§9 item 354" and its 2026-10-10 03:00
# review rulings: C_in,F-human is the P17 (1, +1j) reading, the four-reading
# spread printed beside it; the root is a record, never asserted.
# ============================================================================

ROOT_ENV = "FEM_EM_WF10_STEP1E_ROOT"
#: Coding band for the secant on a linear B(C) (brief: 1e-12).
SECANT_CODING_BAND = 1.0e-12
#: Prediction written first (brief): slope ≈ 16ω S/F (≈ 6.43e-3 S/pF).
PREDICTED_SLOPE_FACTOR = 16.0
#: Part 1's node-118.1186 readings (`20261009T003325Z_WF-10-step1e-W1.log:6322–6323`),
#: printed beside W1a's re-measurement (predicted ≤ 1e-9 rel), never asserted.
PART1_Y_IN = {
    "(1,+1j)": (1.215946766e-02 - 8.803805802e-04j, 1.054684737e-02 - 8.466531661e-04j),
    "(1,-1j)": (1.042881556e-02 - 6.442348932e-04j, 1.228846537e-02 - 1.108258906e-03j),
}
PREDICTED_PART1_REPRO = 1.0e-9


def _sign_change_roots(c_nodes, b_vals):
    """Every sign change of ``b`` over the ascending grid ``c``: a list of
    ``(c_a, c_b, secant_root)``; an exact zero at a node is its own root."""
    c = np.asarray(c_nodes, dtype=float)
    b = np.asarray(b_vals, dtype=float)
    assert np.all(np.diff(c) > 0.0), c
    out = []
    for i in range(len(c) - 1):
        ba, bb = b[i], b[i + 1]
        if ba == 0.0:
            out.append((c[i], c[i], c[i]))
        elif ba * bb < 0.0:
            out.append((c[i], c[i + 1], c[i] - ba * (c[i + 1] - c[i]) / (bb - ba)))
    if len(b) and b[-1] == 0.0:
        out.append((c[-1], c[-1], c[-1]))
    return out


def test_wf10_step1e_bracket_coding():
    """Ungated, synthetic: the secant recovers a planted root of a linear B(C)
    ≤ 1e-12 rel; a no-sign-change grid returns none; two changes return both."""
    grid = np.array([118.1186, 118.20, 118.30, 118.60, 121.00]) * 1e-12
    slope = PREDICTED_SLOPE_FACTOR * 2.0 * np.pi * 64.0e6
    worst = 0.0
    for c_root in (118.255e-12, 118.15e-12, 120.5e-12):
        b = slope * (grid - c_root)
        roots = _sign_change_roots(grid, b)
        assert len(roots) == 1, roots
        c_a, c_b, r = roots[0]
        assert c_a < c_root < c_b, (c_a, c_root, c_b)
        worst = max(worst, abs(r / c_root - 1.0))
    print(
        f"[WF-10 step1e2 coding] planted-root secant worst rel {worst:.3e} "
        f"(band {SECANT_CODING_BAND:g}, ASSERTED)"
    )
    assert worst <= SECANT_CODING_BAND, worst
    # Negative control (asserted): no sign change ⇒ no root.
    none = _sign_change_roots(grid, slope * (grid - 130.0e-12))
    print(f"[WF-10 step1e2 coding] no-sign-change grid returns {none} (ASSERTED empty)")
    assert none == []
    # Two changes: B = (C − a)(C − b), both bracketed and returned.
    a_r, b_r = 118.25e-12, 119.5e-12
    two = _sign_change_roots(grid, (grid - a_r) * (grid - b_r) * 1e24)
    print(
        f"[WF-10 step1e2 coding] two-change grid returns "
        f"{[(x * 1e12, y * 1e12) for x, y, _ in two]} pF brackets (ASSERTED 2)"
    )
    assert len(two) == 2
    assert two[0][0] < a_r < two[0][1] and two[1][0] < b_r < two[1][1], two


def test_wf10_step1e_root_bracket():
    """Numpy only, ``-n 1``: read the five node JSONs and bracket Im Y_in = 0."""
    if os.environ.get(ROOT_ENV, "").strip() != "1":
        pytest.skip(f"{ROOT_ENV} unset: the root reader runs after W1a/W1b")
    nodes_env = os.environ.get(NODES_ENV, "").strip()
    assert nodes_env, f"{NODES_ENV} must list the five nodes"
    nodes_pf = sorted(float(x) for x in nodes_env.split(",") if x.strip())
    from tests.validation.test_wf10_f_human_tuned_b1_sar import _degree_knob

    knob = _degree_knob()
    recs = []
    for c_pf in nodes_pf:
        path = _h17_node_json_path(c_pf, knob)
        assert os.path.exists(path), f"{path} absent"
        with open(path) as fh:
            rec = json.load(fh)
        assert abs(rec["c_f"] * 1e12 / c_pf - 1.0) < 1e-12, (rec["c_f"], c_pf)
        assert rec["r_ohm"] == ESR_OHM, rec["r_ohm"]
        assert rec.get("degree", 1) == (knob or 1), (rec.get("degree"), knob)
        recs.append(rec)
    print(
        f"[WF-10 step1e2 root] degree {knob or 1} (FEM_EM_WF10_DEGREE "
        f"{knob if knob is not None else 'unset'}); node JSONs "
        f"{'h17-node-<pF>-d2.json' if knob else 'h17-node-<pF>.json'}"
    )
    c = np.array([r["c_f"] for r in recs])
    omega = 2.0 * np.pi * float(recs[0]["freq_hz"])
    order = recs[0]["order"]
    feed_names = order[:2]
    print(
        f"[WF-10 step1e2 root] nodes {', '.join(f'{x * 1e12:.4f}' for x in c)} pF; feeds "
        f"{feed_names}; r = {ESR_OHM} Ohm; predicted slope 16w = "
        f"{PREDICTED_SLOPE_FACTOR * omega * 1e-12:.4e} S/pF; feed gap's own w = "
        f"{omega * 1e-12:.4e} S/pF"
    )
    # Part 1 reproduction at 118.1186 (printed, never asserted; degree 1 only —
    # at degree 2 it would compare a degree-2 node to a degree-1 log).
    for rec in (recs if knob is None else []):
        if abs(rec["c_f"] * 1e12 - C_MODE_CU_PF) < 1e-9:
            for lbl, (p17, p21) in PART1_Y_IN.items():
                y = [complex(*v) for v in rec["y_in"][lbl]]
                rr = max(abs(y[0] - p17) / abs(p17), abs(y[1] - p21) / abs(p21))
                print(
                    f"[WF-10 step1e2 root] node 118.1186 a = {lbl}: re-measured {y[0]:.9e} / "
                    f"{y[1]:.9e} S vs part 1 {p17:.9e} / {p21:.9e} S; rel {rr:.3e} "
                    f"(PREDICTED <= {PREDICTED_PART1_REPRO:g}, printed)"
                )
    roots = {}
    for lbl in SENSES:
        for k, name in enumerate(feed_names):
            y = np.array([complex(*r["y_in"][lbl][k]) for r in recs])
            b = y.imag
            print(
                f"[WF-10 step1e2 root] {name} a = {lbl}: Im Y_in "
                + ", ".join(f"{x * 1e12:.4f}:{v:+.6e}" for x, v in zip(c, b))
                + " S; Re Y_in "
                + ", ".join(f"{v:.6e}" for v in y.real)
                + " S"
            )
            slopes = np.diff(b) / np.diff(c)
            print(
                f"[WF-10 step1e2 root] {name} a = {lbl}: interval slopes "
                + ", ".join(f"{s * 1e-12:+.4e}" for s in slopes)
                + " S/pF (= "
                + ", ".join(f"{s / (PREDICTED_SLOPE_FACTOR * omega):.3f}" for s in slopes)
                + " x 16w)"
            )
            found = _sign_change_roots(c, b)
            if not found:
                i_near = int(np.argmin(np.abs(b)))
                print(
                    f"[WF-10 step1e2 root] {name} a = {lbl}: NO sign change in the grid; "
                    f"nearest node {c[i_near] * 1e12:.4f} pF (Im {b[i_near]:+.6e} S) — a record"
                )
            for c_a, c_b, r in found:
                print(
                    f"[WF-10 step1e2 root] {name} a = {lbl}: sign change bracket "
                    f"[{c_a * 1e12:.4f}, {c_b * 1e12:.4f}] pF, secant root {r * 1e12:.6f} pF, "
                    f"resolution {(c_b - c_a) * 1e12:.4f} pF ({(c_b - c_a) / r:.3e} rel)"
                )
            roots[(name, lbl)] = found
    key = (feed_names[0], "(1,+1j)")
    all_first = [v[0][2] for v in roots.values() if v]
    if roots[key]:
        c_in = roots[key][0][2]
        spread = (max(all_first) - min(all_first)) if all_first else float("nan")
        print(
            f"[WF-10 step1e2 root] C_in,F-human ({key[0]} (1,+1j), ruling (2)) = "
            f"{c_in * 1e12:.6f} pF, bracket [{roots[key][0][0] * 1e12:.4f}, "
            f"{roots[key][0][1] * 1e12:.4f}] pF; four-reading roots "
            + ", ".join(
                f"{n} {l}: " + (f"{v[0][2] * 1e12:.6f}" if v else "none")
                for (n, l), v in roots.items()
            )
            + f" pF; spread {spread * 1e12:.6f} pF ({spread / c_in:.3e} rel) — RECORD, "
            f"never asserted"
        )
    else:
        print(f"[WF-10 step1e2 root] C_in,F-human: no sign change on {key} in the grid — a record")


# ============================================================================
# §9 item 459 (§10 H26a): the degree knob's tests (H20a's K0 / K1, adapted).
# ============================================================================


def test_wf10_step1e_degree_knob_coding():
    """(K0) H26a: the node-JSON path coding and the degree-2 gate table (numpy only;
    the knob's parser is H20a's and is tested there)."""
    for c_pf in (118.2, 118.1186, 121.0):
        old = os.path.join(JSON_DIR, f"h17-node-{c_pf:.4f}.json")
        assert _h17_node_json_path(c_pf, None) == old, (c_pf, old)
    p1, p2 = _h17_node_json_path(118.2, None), _h17_node_json_path(118.2, 2)
    print(f"[WF-10 H26a K0] node JSON: unset -> {p1}; degree 2 -> {p2}")
    assert p2.endswith("-d2.json") and os.path.dirname(p2) == JSON_DIR and p2 != p1
    with pytest.raises(ValueError):
        _h17_node_json_path(118.2, 1)
    printed, asserted = STEP1E_DEGREE2_GATES["printed"], STEP1E_DEGREE2_GATES["asserted"]
    print(f"[WF-10 H26a K0] degree-2 gates: printed {printed}; asserted {asserted}")
    assert "(R)" in printed and "(Z50)" in printed
    assert not set(printed) & set(asserted)
    for gate in ("(A0)", "(Deg)", "(P50)", "NEG", "(P_RC)"):
        assert gate in asserted, gate


@complex_only
def test_wf10_step1e_build_only(monkeypatch):
    """(K1) H26a: one build of step 1e's fixture, both paths' unknown counts; no
    form compiled."""
    if os.environ.get(BUILD_ENV, "").strip() != "1":
        pytest.skip(f"{BUILD_ENV} unset: build-only -n 2 F-human window, environment-selected")
    from tests.validation.test_wf10_f_human_tuned_b1_sar import (
        DEGREE_ENV,
        E_RECORD,
        F_RECORD,
        N_DOFS_D2_RECORD,
        PREDICTED_RSS_D2_GIB,
        RSS_D2_BUILD_RECORD_GIB,
        _apply_degree,
    )
    from tests.validation.test_wf12_f_human_hole import _rss_gib

    comm = MPI.COMM_WORLD

    def say(msg):
        if comm.rank == 0:
            print(f"[WF-10 H26a K1] {msg}", flush=True)

    monkeypatch.delenv(DEGREE_ENV, raising=False)
    t0 = time.perf_counter()
    out = _build_root_fixture(comm)
    t_build = time.perf_counter() - t0
    ctx = out["ctx"]
    n1, e, f = (out["counts"][k] for k in ("n_dofs", "n_edges", "n_faces"))
    same = ctx["solver"] is out["solver_built"]
    say(
        f"(Deg-1) unset: degree {out['degree']}, solver is _build_f_human's: {same}; unknowns "
        f"{n1} vs E {e} (record E {E_RECORD}; ASSERTED exact); F {f} (record {F_RECORD}); "
        f"cells {out['cells']}; build {t_build:.1f} s"
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
        f"{PREDICTED_RSS_D2_GIB:g}; H20a 1.10 GiB; WF-12 4a at build {RSS_D2_BUILD_RECORD_GIB} "
        f"GiB; printed)"
    )
    assert deg2 == 2 and solver2.degree == 2 and ctx["solver"] is solver2
    assert n2 == 2 * e + 2 * f, (n2, e, f)
    assert n2 != n1
    assert n2 * e == n1 * (2 * e + 2 * f), (n1, n2, e, f)
