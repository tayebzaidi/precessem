"""Example (`EX-128`): a sigma = 800 S/m rod in the tuned copper F-human's gel phantom.

`IMP-3` step 8 (`tests/validation/test_imp3_f_human_tuned_rod.py`, item 381, landed
`0e7d3abd`) puts one capped rod (r = 1 mm, l = 30 mm, axis z, centre (0.015, 0, 0)) into
I21a's kept-rod copper (Leontovich) F-human hole mesh and switches it by *material only*:
gel (the *without* arm) or sigma = 800 S/m (the *with* arm).  Both arms run on ONE mesh, the
two-port 90 degree drive at the rod-free root C_in = 118.2437 pF.  This example repeats the
gate's body verbatim (so every gate assertion runs, in the gate's print format), then adds
what the gate returns nothing of: the cellwise local-SAR ratio map, the rod's own loss map,
the near-tip / shaft / far ratio readings, and the ParaView outputs.

Provenance of the repeated body: `tests/validation/test_imp3_f_human_tuned_rod.py` at
`0e7d3abd`, `_run_rod` lines 340-847 (imports 341-408 left inside the function), with the
`_arm` closure (548-689) whole.  Three changes only, each marked: the function is named
`_main` and does not return after the asserts; `_say` / `_in` and every helper are imported
from the gate module; `mass_averaged_sar` is wrapped (results untouched) so the tip-ball
results can be read back for the CSV and the record checks.

Not claimed: an accurate tip ratio, a PEC / insulated / real lead, ISO 10974 or ASTM F2182,
heating, absolute / tuned-coil / pointwise SAR, compliance, degree 2 or h-convergence, f0 / Q
/ match / a re-tune with the rod in, AED or Palace.  Degree 1, one fixture, one sigma rod.
"""

from __future__ import annotations

import dataclasses  # noqa: F401  (used by the repeated body)
import sys
import time
from pathlib import Path

import numpy as np
import ufl
from dolfinx import default_scalar_type, fem, la
from mpi4py import MPI

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

import tests.validation.test_imp3_f_human_tuned_rod as g8  # noqa: E402
from tests.validation.test_imp3_f_human_tuned_rod import (  # noqa: E402,F401
    C_STAR_LABEL,
    PRED_COVERAGE,
    PRED_D_SAR10G,
    PRED_ETA_Z_CCW,
    PRED_MASS_1G_G,
    PRED_NEG,
    PRED_P_W,
    PRED_Q_TOP_CCW,
    PRED_C1,
    PRED_C2,
    PRED_REMESH,
    PRED_ROD_S_MOVE,
    PRED_SIGMA_P,
    PRED_T3,
    SAR10G_C_STAR_359,
    _ball_mean_closed_form_tag3,
    _build_f_human_rod,
    _c_star,
    _in,
    _load,
    _read_record,
    _say,
)

from fem_em_solver.io.paraview_utils import adopt_host_ownership, write_xdmf_with_tags  # noqa: E402
from fem_em_solver.post import sar as _sar_mod  # noqa: E402
from fem_em_solver.post.setup_figure import write_setup_figure  # noqa: E402

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"  # committed, EX-57
BASENAME = "implants_10_sigma_rod_in_tuned_f_human_tip_sar"

# ---- W1b records (20261011T200932Z_IMP-3-step8-W1b-rod.log), tolerance 10^-(d-2) ------
# (value, printed significant digits d); complex records are compared as |a-b|/|b|.
W1B = "20261011T200932Z_IMP-3-step8-W1b-rod.log"
REC = {
    # :6350
    "cells": 436541, "n403": 15025, "dofs": 579178,
    # :6352
    "n_lead": 15025, "n_ball": 7974, "v_ball": (1.080564e-07, 7),
    # :6354 (node pair, without) / :6357 (with)
    "node_o": {"p_src": (4.336252680e-02, 10), "p_sheets": (3.896210440e-02, 10),
               "p_phantom": (4.096215247e-03, 10), "p_surf": (3.042036904e-04, 10)},
    "node_w": {"p_src": (4.336218997e-02, 10), "p_sheets": (3.896151148e-02, 10),
               "p_phantom": (4.096331988e-03, 10), "p_surf": (3.041994709e-04, 10)},
    "neg_o": (7.0154e-03, 5), "neg_w": (7.0153e-03, 5),
    # :6355 / :6356 / :6359 / :6360 : P_src, sumP_sheet, P_phantom, P_surf, blind
    "pw": {("o", "ccw"): (8.664505689e-02, 7.784207601e-02, 8.193837110e-03, 6.090942816e-04, 7.0298e-03),
           ("o", "cw"): (8.664505689e-02, 7.783192783e-02, 8.204497083e-03, 6.085900205e-04, 7.0239e-03),
           ("w", "ccw"): (8.664028376e-02, 7.783302386e-02, 8.197568987e-03, 6.090099488e-04, 7.0292e-03),
           ("w", "cw"): (8.664028376e-02, 7.782340442e-02, 8.207742087e-03, 6.085282764e-04, 7.0236e-03)},
    # :6361-6363
    "t3": {"P17": dict(qo=1.6524391923e-01 + 2.2646778905e-01j, qw=3.9574634613e-01 + 4.8522179881e-01j,
                       dq=2.3050242690e-01 + 2.5875400976e-01j, ratio=2.233475, c1=1.9933, c2=0.9306),
           "P21": dict(qo=5.1893178432e-01 + 1.1480394572e+00j, qw=1.3221407802e+00 + 2.5900322174e+00j,
                       dq=8.0320899585e-01 + 1.4419927603e+00j, ratio=2.308146, c1=1.9932, c2=0.9305)},
    "t3_ccw_dq": 1.6724951872e+00 - 5.4445498609e-01j, "t3_ccw_ratio": 2.304690,
    # :6364-6365
    "r_cov_mm": (108.143, 6), "c3": (8.3070e-05, 5),
    # :6366-6371 : sar_o, sar_w, Q, mass_g, eta, M_cf, est, Q/est
    "tips": {("top", "ccw"): (5.774083444e-04, 2.196670608e-03, 3.804362, 0.980607, 0.918468, 2.905124, 2.749796, 1.383507),
             ("top", "P17"): (8.345568409e-05, 1.482919877e-04, 1.776895, 0.980607, 0.459361, 2.905124, 1.875140, 0.947607),
             ("top", "P21"): (4.619880613e-04, 1.887187454e-03, 4.084927, 0.980607, 0.954367, 2.905124, 2.818188, 1.449487),
             ("bottom", "ccw"): (6.111085336e-04, 2.225663836e-03, 3.642011, 0.982115, 0.918468, 2.903149, 2.747982, 1.325340),
             ("bottom", "P17"): (9.617889423e-05, 1.584905128e-04, 1.647872, 0.982115, 0.459361, 2.903149, 1.874232, 0.879225),
             ("bottom", "P21"): (4.716739884e-04, 1.893382869e-03, 4.014177, 0.982115, 0.954367, 2.903149, 2.816303, 1.425336)},
    # :6372-6375
    "a1g_mm": (6.2035, 5), "s10_o": (4.877452318e-03, 10), "s10_w": (4.796853846e-03, 10),
    "dsar": (-1.6525e-02, 5), "remesh": (1.2408e-02, 5), "rod_move": (3.5593e-04, 5),
    "p_rod": (6.809682718e-07, 10), "p_rod_over_ph": (8.3070e-05, 5),
}

_CHK = {"fail": [], "groups": {}}


def _chk(group, name, got, rec, d=None):
    """Exact for int / str / list; otherwise relative at 10^-(d-2)."""
    if d is None:
        ok = got == rec
        rel = 0.0 if ok else float("inf")
    else:
        rel = float(abs(got - rec) / abs(rec))
        ok = rel <= 10.0 ** (-(d - 2))
    g = _CHK["groups"].setdefault(group, [0, 0.0])
    g[0] += 1
    g[1] = max(g[1], rel)
    if not ok:
        _CHK["fail"].append((group, name, str(got)[:40], str(rec)[:40], rel))


# ---- result recorder: wraps mass_averaged_sar, results untouched (change 3) -------------
CALLS = []
_orig_mass_averaged_sar = _sar_mod.mass_averaged_sar


def _recording_mass_averaged_sar(*args, **kwargs):
    r = _orig_mass_averaged_sar(*args, **kwargs)
    CALLS.append((kwargs.get("center"), kwargs.get("radius"), r))
    return r


def _dg0_vec(form):
    """Assemble a DG0-test-function linear form to a per-cell real vector (`mri:17`'s)."""
    b = fem.assemble_vector(fem.form(form))
    b.scatter_reverse(la.InsertMode.add)
    b.scatter_forward()
    return np.real(b.array)


def _main(comm):  # change 1: was `_run_rod`; execution continues after the asserts
    rank0 = comm.rank == 0

    def P(msg):
        if rank0:
            print(msg, flush=True)

    # ---- structural start-up asserts ------------------------------------------------
    assert comm.size == 8, comm.size
    if not np.issubdtype(np.dtype(default_scalar_type), np.complexfloating):
        raise RuntimeError("complex build required (source dolfinx-complex-mode)")
    _sar_mod.mass_averaged_sar = _recording_mass_averaged_sar
    # ---- BEGIN verbatim body: tests/validation/test_imp3_f_human_tuned_rod.py @ 0e7d3abd,
    # ---- lines 341-847 of `_run_rod` (imports 341-408 kept inside the function) ----
    import ufl
    from dolfinx import fem

    from fem_em_solver.core import HomogeneousMaterial, TimeHarmonicSolver
    from fem_em_solver.core.cavity import surface_resistance_ohm
    from fem_em_solver.core.time_harmonic import build_material_fields
    from fem_em_solver.io.mesh import BIRDCAGE_IMPLANT_CELL_TAG
    from fem_em_solver.ports.lumped import (
        lumped_port_bilinear_term,
        lumped_port_linear_term,
        sheet_terminal_current,
    )
    from fem_em_solver.post.evaluation import evaluate_vector_field_parallel
    from fem_em_solver.post.sar import (
        averaging_ball_radius,
        build_density_field,
        mass_averaged_sar,
        mean_sar,
    )

    from tests.mesh.helpers import global_cell_tag_set
    from tests.mesh.test_birdcage_implant_rod import ROD_CENTRE, _cell_values_with_ghosts
    from tests.mesh.test_birdcage_implant_rod_f_human import (
        F_HUMAN_KEPT_ROD_CELL_GUARD,
        F_HUMAN_KEPT_ROD_CELL_RECORD,
        ROD_FH,
    )
    from tests.mesh.test_birdcage_port_sheet_prerequisite import CELL_COUNT_BAND
    from tests.validation.test_birdcage_b1_plus_map import PHANTOM_RHO_KG_PER_M3
    from tests.validation.test_birdcage_f_human_rung import F_HUMAN_RING_RADIUS, _ring_ports
    from tests.validation.test_birdcage_power_identity import (
        DISCRETE_IDENTITY_RTOL,
        _sheet_field_dissipation_w,
        _source_power_w,
    )
    from tests.validation.test_birdcage_sar_mass_averaged import EXACT_IDENTITY_RTOL
    from tests.validation.test_imp3_rod_birdcage import AXIS_SAMPLES, PRED_Q_OVER_ESTIMATE, TIPS
    from tests.validation.test_mass_averaged_sar_standard_masses import (
        ONE_GRAM_KG,
        QUADRATURE_DEGREE,
        TEN_GRAM_KG,
    )
    from tests.validation.test_port_birdcage_four_port import TERMINATED_PORT_IMPEDANCE_OHM
    from tests.validation.test_port_birdcage_lumped_column import PHANTOM_CELL_TAG
    from tests.validation.test_port_birdcage_ring_column import AZIMUTH_MATCH_DEG
    from tests.validation.test_port_circuit_layer_field import _residual
    from tests.validation.test_port_drive_superposition import _ring_quadrature_slots
    from tests.validation.test_port_gap_voltage_impedance import SIGMA_WIRE_S_PER_M
    from tests.validation.test_port_package_sparameters import REFERENCE_IMPEDANCE_OHM
    from tests.validation.test_th14_birdcage_copper import _leontovich_term, _surface_loss_w
    from tests.validation.test_wf7_f_human_b1_quadrature import DEGREE, FREQUENCY_HZ
    from tests.validation.test_wf7_f_human_sar10g_hotspot import (
        CONDUCTOR_CLEARANCE_M,
        COVERAGE_BALL_PAD_M,
        _phantom_geometry,
    )
    from tests.validation.test_wf9_transfer_function_identity import BALL_RADIUS, _free
    from tests.validation.test_wf10_f_human_copper_esr_root import (
        ESR_OHM,
        _current_route_s,
        _source_voltage,
        _thevenin_feed_impedance,
        _z_rc,
    )
    from tests.validation.test_wf10_f_human_tuned_b1_sar import C_IN_F_HUMAN_PF, _sense_weights
    from tests.validation.test_wf10_f_human_tuned_circuit import QUADRATURE_SLOT_OFFSET, _respec
    from tests.validation.test_wf12_f_human_copper_matrix import BLIND_FACTOR, COPPER_SIGMA
    from tests.validation.test_wf12_f_human_hole import _census_checks, _owned_census, _rss_gib

    t_w = time.perf_counter()
    times = {}
    s_rec = _read_record(comm)
    load = _load()

    # ---- build, (A0r) ------------------------------------------------------------
    t0 = time.perf_counter()
    built = _build_f_human_rod(comm, ROD_FH)
    times["build"] = time.perf_counter() - t0
    ctx = built["ctx"]
    msh, tags_f, omega, cell_tags = ctx["msh"], ctx["tags_f"], ctx["omega"], ctx["cell_tags"]
    n_cells = int(built["cells"])
    assert n_cells < F_HUMAN_KEPT_ROD_CELL_GUARD, n_cells
    assert abs(omega / (2.0 * np.pi * float(FREQUENCY_HZ)) - 1.0) < 1e-12
    census = _owned_census(msh, cell_tags, comm)
    census_fails = _census_checks(census, n_cells, _ring_ports())
    n403 = int(census.get(BIRDCAGE_IMPLANT_CELL_TAG, 0))
    cell_rel = n_cells / F_HUMAN_KEPT_ROD_CELL_RECORD - 1.0
    msh.topology.create_entities(1)
    n_edges = int(msh.topology.index_map(1).size_global)
    _say(f"[rod] (A0r) cells {n_cells} vs F_HUMAN_KEPT_ROD_CELL_RECORD "
         f"{F_HUMAN_KEPT_ROD_CELL_RECORD} rel {cell_rel:.3e} (band {CELL_COUNT_BAND:g}, ASSERTED; "
         f"guard {F_HUMAN_KEPT_ROD_CELL_GUARD}); census failures {census_fails or 'none'}; tag-403 "
         f"cells {n403} (PREDICTED ≈ 15 000); E = degree-1 dofs {n_edges}; build "
         f"{times['build']:.1f} s")
    assert abs(cell_rel) <= CELL_COUNT_BAND, cell_rel
    assert not census_fails, census_fails
    assert n403 > 0

    # ---- problems, (M0) ---------------------------------------------------------
    p_without = built["problem"]
    mmap = dict(p_without.material_map)
    gel = mmap[PHANTOM_CELL_TAG]
    p_with = dataclasses.replace(p_without, material_map={
        **mmap, BIRDCAGE_IMPLANT_CELL_TAG: HomogeneousMaterial(
            sigma=SIGMA_WIRE_S_PER_M, epsilon_r=gel.epsilon_r, mu_r=gel.mu_r)})
    delta_sigma = float(SIGMA_WIRE_S_PER_M - gel.sigma)

    tdim = msh.topology.dim
    imap = msh.topology.index_map(tdim)
    n_owned = imap.size_local
    cells = np.arange(n_owned + imap.num_ghosts, dtype=np.int32)
    mid = msh.geometry.x[msh.geometry.dofmap[cells]].mean(axis=1)
    tag = _cell_values_with_ghosts(msh, cell_tags)
    lead = tag == BIRDCAGE_IMPLANT_CELL_TAG
    gelc = tag == PHANTOM_CELL_TAG
    ballc = gelc & (np.linalg.norm(mid - np.asarray(TIPS["top"]), axis=1) < BALL_RADIUS)

    def _minmax(problem):
        sigma_f, eps_f = build_material_fields(
            msh, problem.material, cell_tags=problem.cell_tags,
            material_map=problem.material_map, phantom_material=problem.phantom_material,
            phantom_tag=problem.phantom_tag)
        owned_lead = np.flatnonzero(lead[:n_owned])
        res = {}
        for name, f in (("sigma", sigma_f), ("eps", eps_f)):
            dofs = f.function_space.dofmap.list[:, 0]
            vals = np.real(f.x.array[dofs[owned_lead]])
            lo = float(vals.min()) if vals.size else np.inf
            hi = float(vals.max()) if vals.size else -np.inf
            res[name] = (comm.allreduce(lo, op=MPI.MIN), comm.allreduce(hi, op=MPI.MAX))
        return res

    m_with, m_without = _minmax(p_with), _minmax(p_without)
    _say(f"[rod] (M0) with: {m_with}; without: {m_without}; gel sigma {gel.sigma} eps_r "
         f"{gel.epsilon_r}; delta_sigma {delta_sigma:.9f} (ASSERTED min = max)")
    assert m_with["sigma"] == (SIGMA_WIRE_S_PER_M, SIGMA_WIRE_S_PER_M), m_with
    assert m_with["eps"] == (gel.epsilon_r, gel.epsilon_r), m_with
    assert m_without["sigma"] == (gel.sigma, gel.sigma), m_without
    assert m_without["eps"] == (gel.epsilon_r, gel.epsilon_r), m_without

    # ---- masks ------------------------------------------------------------------
    q0 = fem.functionspace(msh, ("DG", 0))
    dof_of_cell = q0.dofmap.list[:, 0]

    def mask(flags, value=1.0):
        f = fem.Function(q0)
        f.x.array[:] = 0.0
        f.x.array[dof_of_cell[cells]] = value * flags.astype(float)
        f.x.scatter_forward()
        return f

    def owned_count(flags):
        return comm.allreduce(int(np.count_nonzero(flags[:n_owned])), op=MPI.SUM)

    chi_lead, chi_b, sigma3 = mask(lead), mask(ballc), mask(gelc, float(gel.sigma))
    n_lead_c, n_ball_c = owned_count(lead), owned_count(ballc)
    disjoint = comm.allreduce(bool(not np.any(ballc & lead)), op=MPI.LAND)
    dxx = ufl.dx(domain=msh)

    def _scalar(expr):
        return comm.allreduce(fem.assemble_scalar(fem.form(expr)), op=MPI.SUM)

    v_ball = float(np.real(_scalar(chi_b * dxx)))
    _say(f"[rod] masks: lead {n_lead_c}, ball {n_ball_c} cells (V_b {v_ball:.6e} m^3), disjoint "
         f"{disjoint} (ASSERTED)")
    assert n_lead_c > 0 and n_ball_c > 0 and disjoint

    # ---- shared coil terms -------------------------------------------------------
    z0 = float(REFERENCE_IMPEDANCE_OHM)
    z_s = (1.0 + 1.0j) * surface_resistance_ohm(omega, COPPER_SIGMA)
    wall = _leontovich_term(msh, built["facets"], z_s, omega)
    z_p = complex(TERMINATED_PORT_IMPEDANCE_OHM)
    dx = ufl.Measure("dx", domain=msh, subdomain_data=cell_tags)
    acct = {"mesh": msh, "facet_tags": tags_f}
    base_specs = list(ctx["specs"])
    c_f = C_IN_F_HUMAN_PF * 1e-12

    def _vol(sigma_f, e, measure):
        form = fem.form(0.5 * sigma_f * ufl.inner(e, e) * measure)
        return float(np.real(comm.allreduce(fem.assemble_scalar(form), op=MPI.SUM)))

    # Order block lifted from `main` at 232a6fe2,
    # tests/validation/test_wf10_f_human_tuned_b1_sar.py lines 416–436, verbatim.
    az_ref, slots = _ring_quadrature_slots(built["sheets"])
    ref_k, ref_sign = slots["P17"]
    feed2 = next(pid for pid, (k, sg) in slots.items()
                 if sg == ref_sign and k == (ref_k + QUADRATURE_SLOT_OFFSET) % 16)
    others = sorted((pid for pid in slots if pid not in ("P17", feed2)),
                    key=lambda p: int(p[1:]))
    order = ["P17", feed2] + others
    feeds = order[:2]
    feeds_idx = [0, 1]
    assert len(order) == 32 and sorted(order) == sorted(sp.port_id for sp in base_specs)
    az_by = {f"P{s['ordinal']}": float(s["azimuth_deg"]) for s in built["sheets"]}
    az_feeds = (az_by[feeds[0]], az_by[feeds[1]])
    d_az = (az_feeds[1] - az_feeds[0]) % 360.0
    assert abs(d_az - 90.0) < AZIMUTH_MATCH_DEG, d_az
    assert slots[feeds[0]][0] == 0 and slots[feeds[1]][0] == 4
    _say(f"[rod] feeds {feeds}, az {az_feeds}, d_az {d_az:.9f} deg")

    z_rc = _z_rc(c_f, FREQUENCY_HZ)
    z_th = _thevenin_feed_impedance(z_rc, z0)
    v_s = _source_voltage(z_rc, z0)
    imp = {pid: z_rc for pid in others}
    imp.update({pid: z_th for pid in feeds})
    w_sense = _sense_weights()

    def _arm(problem, label, solver, *, ell=None):
        keep = {}

        # Lifted (L1) from `main` at 232a6fe2,
        # tests/validation/test_wf10_f_human_tuned_b1_sar.py lines 364–413 (`_drive`)
        # verbatim except that ``solver`` is the arm's solver, not ``ctx["solver"]``.
        def _drive(specs, driven_id, *, held, partition, keep_as=None):
            sheets = [sp.sheet(driven=(sp.port_id == driven_id)) for sp in specs]
            dsheet = next(s for s in sheets if s.port_id == driven_id)
            linear = [lambda test, _s=dsheet: lumped_port_linear_term(
                msh, tags_f, _s, test, omega_rad_per_s=omega)]
            comm.Barrier()
            t0 = time.perf_counter()
            if held:
                fields = solver.solve_with_held_factorization(extra_linear_terms=linear)
            else:
                bil = [wall] + [
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
            currents = {s.port_id: sheet_terminal_current(msh, tags_f, s, e, comm)
                        for s in sheets}
            v_src = complex(dsheet.source_voltage_v)
            s_col = {}
            for s in sheets:
                v = complex(s.source_voltage_v) - currents[s.port_id] * z_p
                s_col[s.port_id] = (v - z_p * currents[s.port_id]) / v_src
            out = {"s": s_col, "currents": currents, "kind": "held" if held else "solve",
                   "t_solve": t_solve}
            if "sigma" not in keep:
                keep["sigma"] = fields.sigma_field
            if partition:
                sigma_f = fields.sigma_field
                by_id = {s.port_id: float(_sheet_field_dissipation_w(acct, s, e, omega))
                         for s in sheets}
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

        # Lifted `_pair` (same file, lines 440–453), verbatim.
        def _pair(imp_x, z_feed, v_feed, keep_fields=False):
            specs_x = _respec(base_specs, None, imp_x)
            assert all(sp.width_correction_kappa is None for sp in specs_x)
            cur, knd, tt, parts = [], [], [], []
            for n_drive, pid in enumerate(feeds):
                d = _drive(specs_x, pid, held=(n_drive > 0), partition=(n_drive == 0),
                           keep_as=(f"E_{pid}" if keep_fields else None))
                knd.append(d["kind"])
                tt.append(d["t_solve"] + d["t_post"])
                cur.append(np.array([d["currents"][p] for p in order], dtype=np.complex128))
                if n_drive == 0:
                    parts.append({k: v for k, v in d.items() if k.startswith("p_")})
                del d
            return (_current_route_s(cur, feeds_idx, z_feed, v_feed), cur, knd, tt, parts[0],
                    specs_x)

        # Lifted `_partition` (same file, lines 455–459), verbatim.
        def _partition(p):
            p_vol = p["p_phantom"] + p["p_rest"]
            rel = abs(p["p_src"] - p["p_sheets"] - p_vol - p["p_surf"]) / abs(p["p_src"])
            blind = abs(p["p_src"] - p["p_sheets"] - p_vol) / abs(p["p_src"])
            return rel, blind

        t0 = time.perf_counter()
        s_node, cur_node, kinds, t_n, pn, specs_n = _pair(imp, z_th, v_s, keep_fields=True)
        assert kinds == ["solve", "held"], kinds
        p_rc, blind = _partition(pn)
        _say(f"[rod/{label}] node pair {kinds}, {', '.join(f'{t:.2f}' for t in t_n)} s; (P_RC) "
             f"P_src {pn['p_src']:.9e} sumP_sheet {pn['p_sheets']:.9e} P_phantom "
             f"{pn['p_phantom']:.9e} P_rest {pn['p_rest']:.9e} P_surf {pn['p_surf']:.9e} W; rel "
             f"{p_rc:.3e} (band {DISCRETE_IDENTITY_RTOL:g}, ASSERTED; PREDICTED <= {PRED_P_W:g}); "
             f"node NEG {blind:.4e} = {blind / DISCRETE_IDENTITY_RTOL:.1f}x the band (bar "
             f"{BLIND_FACTOR:.0f}x; PRINTED per ruling (3); PREDICTED ≈ {PRED_NEG:g})")
        w_field = None
        if ell is not None:
            t1 = time.perf_counter()
            f_ell = solver.solve_with_held_factorization(extra_linear_terms=[ell])
            w_field = fem.Function(f_ell.e_complex.function_space, name="w")
            w_field.x.array[:] = f_ell.e_complex.x.array
            w_field.x.scatter_forward()
            del f_ell
            _say(f"[rod/{label}] ell held on the P17 node factorisation: "
                 f"{time.perf_counter() - t1:.2f} s")
        e17, e21 = keep.pop(f"E_{feeds[0]}"), keep.pop(f"E_{feeds[1]}")
        sigma_solver = keep["sigma"]
        specs_by = {sp.port_id: sp for sp in specs_n}
        e_s, parts, lin = {}, {}, {}
        for sense in ("ccw", "cw"):
            w17, w21 = w_sense[sense]
            ef = fem.Function(e17.function_space, name=f"E_{sense}")
            ef.x.array[:] = w17 * e17.x.array + w21 * e21.x.array
            ef.x.scatter_forward()
            e_s[sense] = ef
            wf = {feeds[0]: w17, feeds[1]: w21}
            cur_lin = w17 * cur_node[0] + w21 * cur_node[1]
            sheets = {pid: (dataclasses.replace(specs_by[pid].sheet(driven=True),
                                                source_voltage_v=complex(wf[pid]))
                            if pid in feeds else specs_by[pid].sheet(driven=False))
                      for pid in order}
            cur_read = np.array([sheet_terminal_current(msh, tags_f, sheets[pid], ef, comm)
                                 for pid in order], dtype=np.complex128)
            r_lin = float(np.max(np.abs(cur_read - cur_lin)) / np.max(np.abs(cur_lin)))
            p = {"p_src": float(sum(_source_power_w(acct, sheets[pid], ef, omega)
                                    for pid in feeds))}
            p["p_sheets"] = float(sum(_sheet_field_dissipation_w(acct, sheets[pid], ef, omega)
                                      for pid in order))
            p_total = _vol(sigma_solver, ef, dx)
            p["p_phantom"] = _vol(sigma_solver, ef, dx(PHANTOM_CELL_TAG))
            p["p_rest"] = p_total - p["p_phantom"]
            p["p_surf"] = _surface_loss_w(msh, built["facets"], z_s, ef, comm)
            p["p_rod"] = _vol(sigma_solver, ef, dx(BIRDCAGE_IMPLANT_CELL_TAG))
            p["rel"], p["blind"] = _partition(p)
            parts[sense] = p
            lin[sense] = r_lin
            _say(f"[rod/{label}] [{sense}] (P_w) P_src {p['p_src']:.9e} sumP_sheet "
                 f"{p['p_sheets']:.9e} P_phantom {p['p_phantom']:.9e} P_rest {p['p_rest']:.9e} "
                 f"(P_rod {p['p_rod']:.9e}) P_surf {p['p_surf']:.9e} W; rel {p['rel']:.3e} (band "
                 f"{DISCRETE_IDENTITY_RTOL:g}, ASSERTED); blind {p['blind']:.4e}; (Lin) "
                 f"{r_lin:.3e} (printed)")
        times[f"{label}_arm"] = time.perf_counter() - t0
        return {"s_node": np.asarray(s_node, dtype=np.complex128), "e": e_s, "e17": e17,
                "e21": e21, "sigma": sigma_solver, "pn": pn, "p_rc": p_rc, "blind": blind,
                "parts": parts, "lin": lin, "w": w_field}

    # ---- the without arm ------------------------------------------------------------
    solver = ctx.pop("solver")
    arm_o = _arm(p_without, "without", solver)
    _free(solver)
    del solver
    # ---- the with arm ---------------------------------------------------------------
    zhat = ufl.as_vector((0.0, 0.0, 1.0))
    ell = lambda v: ufl.inner((chi_b / v_ball) * zhat, v) * ufl.dx  # noqa: E731
    solver = TimeHarmonicSolver(p_with, degree=DEGREE)
    arm_w = _arm(p_with, "with", solver, ell=ell)
    _free(solver)
    del solver
    w = arm_w["w"]

    # ---- (T3) + controls -----------------------------------------------------------
    def q(e):
        return complex(_scalar(chi_b * e[2] * dxx)) / v_ball

    t3 = {}
    for d, key in (("P17", "e17"), ("P21", "e21")):
        eo, ew = arm_o[key], arm_w[key]
        qo, qw = q(eo), q(ew)
        dq = qw - qo
        fold = load * delta_sigma * complex(_scalar(chi_lead * ufl.dot(w, eo) * dxx))
        c1 = load * delta_sigma * complex(_scalar(chi_lead * ufl.dot(ufl.conj(w), eo) * dxx))
        c2 = load * delta_sigma * complex(_scalar(chi_lead * ufl.dot(w, ew) * dxx))
        rel = abs(fold - dq) / abs(dq)
        f1, f2 = abs(c1 - dq) / abs(dq), abs(c2 - dq) / abs(dq)
        t3[d] = dict(qo=qo, qw=qw, dq=dq, fold=fold, rel=rel, c1=f1, c2=f2)
        _say(f"[rod] (T3) {d}: q_without {qo:.10e} q_with {qw:.10e} dq {dq:.10e} fold "
             f"{fold:.10e}; rel {rel:.3e} (band {DISCRETE_IDENTITY_RTOL:g}, ASSERTED; PREDICTED <= "
             f"{PRED_T3:g}); |q_with/q_without| {abs(qw / qo):.6f}; C1 (conj w) miss {f1:.4f} "
             f"(PREDICTED {PRED_C1} {_in(f1, PRED_C1)}; wired ASSERTED); C2 (E_with) miss "
             f"{f2:.4f} (PREDICTED {PRED_C2} {_in(f2, PRED_C2)}; wired ASSERTED)")
    wc = w_sense["ccw"]
    dq_ccw = wc[0] * t3["P17"]["dq"] + wc[1] * t3["P21"]["dq"]
    fold_ccw = wc[0] * t3["P17"]["fold"] + wc[1] * t3["P21"]["fold"]
    q_ccw_o = wc[0] * t3["P17"]["qo"] + wc[1] * t3["P21"]["qo"]
    q_ccw_w = wc[0] * t3["P17"]["qw"] + wc[1] * t3["P21"]["qw"]
    _say(f"[rod] (T3) ccw by linearity: dq {dq_ccw:.10e} fold {fold_ccw:.10e} rel "
         f"{abs(fold_ccw - dq_ccw) / abs(dq_ccw):.3e}; |q_with/q_without| "
         f"{abs(q_ccw_w / q_ccw_o):.6f} (printed)")

    # ---- (Σ_P), coverage --------------------------------------------------------------
    centre, radius, half_height, _lo, _hi = _phantom_geometry(msh, cell_tags, PHANTOM_CELL_TAG,
                                                              comm)
    all_tags = global_cell_tag_set(msh, cell_tags)
    rho = build_density_field(
        msh, PHANTOM_RHO_KG_PER_M3, cell_tags=cell_tags,
        density_map={int(t): 0.0 for t in sorted(all_tags) if int(t) != PHANTOM_CELL_TAG})
    coverage_radius = float(np.hypot(radius, half_height) + COVERAGE_BALL_PAD_M)
    conductor_limit = F_HUMAN_RING_RADIUS - CONDUCTOR_CLEARANCE_M
    assert coverage_radius < conductor_limit, (coverage_radius, conductor_limit)
    cen = tuple(float(v) for v in centre)
    sp_cov = {}
    for label, arm in (("without", arm_o), ("with", arm_w)):
        e = arm["e"]["ccw"]
        ms = mean_sar(e, sigma=arm["sigma"], rho=PHANTOM_RHO_KG_PER_M3, cell_tags=cell_tags,
                      subdomain_ids=PHANTOM_CELL_TAG, comm=comm)
        r_sp = abs(ms["dissipated_power_w"] - arm["parts"]["ccw"]["p_phantom"]) / abs(
            arm["parts"]["ccw"]["p_phantom"])
        whole = mass_averaged_sar(e, sigma=sigma3, rho=rho, center=cen, radius=coverage_radius,
                                  comm=comm, quadrature_degree=QUADRATURE_DEGREE)
        tagged = mean_sar(e, sigma=sigma3, rho=PHANTOM_RHO_KG_PER_M3, cell_tags=cell_tags,
                          subdomain_ids=PHANTOM_CELL_TAG, comm=comm)
        pmiss = abs(whole["dissipated_power_w"] / tagged["dissipated_power_w"] - 1.0)
        mmiss = abs(whole["mass_kg"] / (PHANTOM_RHO_KG_PER_M3 * tagged["volume_m3"]) - 1.0)
        c3 = None
        if label == "with":
            whole_s = mass_averaged_sar(e, sigma=arm["sigma"], rho=rho, center=cen,
                                        radius=coverage_radius, comm=comm,
                                        quadrature_degree=QUADRATURE_DEGREE)
            c3 = abs(whole_s["dissipated_power_w"] / tagged["dissipated_power_w"] - 1.0)
        sp_cov[label] = (r_sp, pmiss, mmiss)
        _say(f"[rod/{label}] (Sigma_P) mean_sar(solver sigma, tag 3) vs partition P_phantom: rel "
             f"{r_sp:.3e} (band {EXACT_IDENTITY_RTOL:g}, ASSERTED; PREDICTED <= {PRED_SIGMA_P:g}); "
             f"coverage (sigma3, R_cov {coverage_radius * 1e3:.3f} mm < {conductor_limit * 1e3:.3f}"
             f" mm): power miss {pmiss:.3e}, mass miss {mmiss:.3e} (band {EXACT_IDENTITY_RTOL:g}, "
             f"ASSERTED; PREDICTED <= {PRED_COVERAGE:g})"
             + (f"; C3 solver-sigma whole ball miss {c3:.4e} vs P_rod/P_phantom "
                f"{arm['parts']['ccw']['p_rod'] / arm['parts']['ccw']['p_phantom']:.4e} "
                f"(predicted, printed)" if c3 is not None else ""))

    # ---- readings: tip 1 g ratio, η_z, M_cf, estimate ---------------------------------
    a1 = averaging_ball_radius(mass_kg=ONE_GRAM_KG, rho=PHANTOM_RHO_KG_PER_M3)
    a10 = averaging_ball_radius(mass_kg=TEN_GRAM_KG, rho=PHANTOM_RHO_KG_PER_M3)
    zs = np.linspace(-0.016, 0.016, AXIS_SAMPLES)
    pts = np.column_stack([np.full_like(zs, ROD_CENTRE[0]), np.zeros_like(zs), zs])
    fields_o = {"ccw": arm_o["e"]["ccw"], "P17": arm_o["e17"], "P21": arm_o["e21"]}
    fields_w = {"ccw": arm_w["e"]["ccw"], "P17": arm_w["e17"], "P21": arm_w["e21"]}
    eta = {}
    for d, f in fields_o.items():
        vals, ok = evaluate_vector_field_parallel(f, pts, owner="lowest_original_cell")
        vals = comm.bcast(vals if comm.rank == 0 else None, root=0)
        ok = comm.bcast(ok if comm.rank == 0 else None, root=0)
        assert bool(np.all(ok)), "axis sample off-mesh"
        ev = np.asarray(vals).reshape(len(pts), 3)
        eta[d] = float(np.sum(np.abs(ev[:, 2]) ** 2) / np.sum(np.abs(ev) ** 2))
    m_cf = {tip: _ball_mean_closed_form_tag3(msh, cell_tags, TIPS[tip], a1) for tip in TIPS}
    t0 = time.perf_counter()
    for tip in ("top", "bottom"):
        c = tuple(float(v) for v in TIPS[tip])
        for d in ("ccw", "P17", "P21"):
            ro = mass_averaged_sar(fields_o[d], sigma=sigma3, rho=rho, center=c, radius=a1,
                                   comm=comm, quadrature_degree=QUADRATURE_DEGREE)
            rw = mass_averaged_sar(fields_w[d], sigma=sigma3, rho=rho, center=c, radius=a1,
                                   comm=comm, quadrature_degree=QUADRATURE_DEGREE)
            qq = rw["averaged_sar_w_per_kg"] / ro["averaged_sar_w_per_kg"]
            est = eta[d] * m_cf[tip][0] + (1.0 - eta[d])
            extra = (f" (PREDICTED {PRED_Q_TOP_CCW} {_in(qq, PRED_Q_TOP_CCW)})"
                     if (tip, d) == ("top", "ccw") else "")
            _say(f"[rod] tip {tip:<6s} {d:<4s}: SAR1g without {ro['averaged_sar_w_per_kg']:.9e} "
                 f"with {rw['averaged_sar_w_per_kg']:.9e} W/kg; Q {qq:.6f}{extra}; masses "
                 f"{ro['mass_kg'] * 1e3:.6f} / {rw['mass_kg'] * 1e3:.6f} g (PREDICTED ≈ "
                 f"{PRED_MASS_1G_G}); eta_z {eta[d]:.6f} M_cf {m_cf[tip][0]:.6f} est {est:.6f} "
                 f"Q/est {qq / est:.6f} (PREDICTED {PRED_Q_OVER_ESTIMATE} "
                 f"{_in(qq / est, PRED_Q_OVER_ESTIMATE)}); reading")
    times["tip_balls"] = time.perf_counter() - t0
    _say(f"[rod] eta_z ccw {eta['ccw']:.6f} (PREDICTED {PRED_ETA_Z_CCW} "
         f"{_in(eta['ccw'], PRED_ETA_Z_CCW)}); a1g {a1 * 1e3:.4f} mm")

    # ---- readings: SAR10g(c*), residuals, P_rod -----------------------------------------
    c_st = _c_star(az_feeds[0], float(centre[2]))
    s10 = {}
    for label, f in (("without", fields_o["ccw"]), ("with", fields_w["ccw"])):
        s10[label] = float(mass_averaged_sar(
            f, sigma=sigma3, rho=rho, center=c_st, radius=a10, comm=comm,
            quadrature_degree=QUADRATURE_DEGREE)["averaged_sar_w_per_kg"])
    dsar = s10["with"] / s10["without"] - 1.0
    _say(f"[rod] SAR10g(c* {C_STAR_LABEL}) ccw sigma3: without {s10['without']:.9e} (359 "
         f"{SAR10G_C_STAR_359:.9e}, remesh rel {s10['without'] / SAR10G_C_STAR_359 - 1.0:.3e}) with "
         f"{s10['with']:.9e} W/kg; with/without - 1 = {dsar:.4e} (PREDICTED |.| in {PRED_D_SAR10G} "
         f"{_in(abs(dsar), PRED_D_SAR10G)}; reading)")
    r_remesh = _residual(arm_o["s_node"], s_rec)
    r_rod = _residual(arm_w["s_node"], arm_o["s_node"])
    _say(f"[rod] remesh move _residual(s_without, s_rec) {r_remesh:.4e} (PREDICTED {PRED_REMESH} "
         f"{_in(r_remesh, PRED_REMESH)}); rod move _residual(s_with, s_without) {r_rod:.4e} "
         f"(PREDICTED <= {PRED_ROD_S_MOVE:g}); readings")
    pw = arm_w["parts"]["ccw"]
    _say(f"[rod] P_rod (with, ccw) {pw['p_rod']:.9e} W; P_rod/P_phantom "
         f"{pw['p_rod'] / pw['p_phantom']:.4e}; ESR {ESR_OHM} Ohm, C_in {C_IN_F_HUMAN_PF} pF")
    _m, _p, rss = _rss_gib(comm)
    _say(f"[rod] times {dict((k, round(v, 1)) for k, v in times.items())}; window "
         f"{time.perf_counter() - t_w:.1f} s; summed RSS {rss:.2f} GiB; -n {comm.size}")

    # ---- asserts (after every print) ----------------------------------------------------
    for arm in (arm_o, arm_w):
        assert arm["p_rc"] <= DISCRETE_IDENTITY_RTOL, arm["p_rc"]
        for s in ("ccw", "cw"):
            assert arm["parts"][s]["rel"] <= DISCRETE_IDENTITY_RTOL, (s, arm["parts"][s]["rel"])
    for d in ("P17", "P21"):
        assert t3[d]["rel"] <= DISCRETE_IDENTITY_RTOL, (d, t3[d]["rel"])
        assert t3[d]["c1"] > DISCRETE_IDENTITY_RTOL and t3[d]["c2"] > DISCRETE_IDENTITY_RTOL
    for label, (r_sp, pmiss, mmiss) in sp_cov.items():
        assert r_sp <= EXACT_IDENTITY_RTOL, (label, r_sp)
        assert pmiss <= EXACT_IDENTITY_RTOL, (label, pmiss)
        assert mmiss <= EXACT_IDENTITY_RTOL, (label, mmiss)
    # ---- END verbatim body (line 847) ----

    # =================================================================================
    # ---- the example's own section (everything below is new; the gate's asserts held)
    # =================================================================================
    _sar_mod.mass_averaged_sar = _orig_mass_averaged_sar
    from tests.mesh.test_birdcage_implant_rod import ROD
    from tests.mesh.test_birdcage_port_sheets import PORT_LOWER, PORT_UPPER
    from tests.validation.test_port_birdcage_lumped_column import CONDUCTOR_CELL_TAG

    assert ROD_FH == dict(ROD, keep_cells=True)
    top_expected = tuple(np.asarray(ROD_CENTRE, dtype=float)
                         + np.array([0.0, 0.0, ROD["length"] / 2.0 + ROD["radius"]]))
    assert top_expected == tuple(float(x) for x in TIPS["top"]), (top_expected, TIPS["top"])

    # ---- record checks (the W1b log, tolerance 10^-(d-2)) -------------------------------
    _chk("A0r", "cells", n_cells, REC["cells"])
    _chk("A0r", "n403", n403, REC["n403"])
    _chk("A0r", "dofs", n_edges, REC["dofs"])
    _chk("masks", "lead", n_lead_c, REC["n_lead"])
    _chk("masks", "ball", n_ball_c, REC["n_ball"])
    _chk("masks", "v_ball", v_ball, *REC["v_ball"])
    _chk("feeds", "feeds", feeds, ["P17", "P21"])
    for key, arm in (("node_o", arm_o), ("node_w", arm_w)):
        for k, (val, d) in REC[key].items():
            _chk(key, k, arm["pn"][k], val, d)
    _chk("node_neg", "o", arm_o["blind"], *REC["neg_o"])
    _chk("node_neg", "w", arm_w["blind"], *REC["neg_w"])
    for (ak, sense), vals in REC["pw"].items():
        pp = (arm_o if ak == "o" else arm_w)["parts"][sense]
        for k, val, d in zip(("p_src", "p_sheets", "p_phantom", "p_surf", "blind"), vals,
                             (10, 10, 10, 10, 5)):
            _chk("P_w", f"{ak}/{sense}/{k}", pp[k], val, d)
    for dd, rr in REC["t3"].items():
        for k in ("qo", "qw", "dq"):
            _chk("T3", f"{dd}/{k}", t3[dd][k], rr[k], 11)
        _chk("T3", f"{dd}/fold", t3[dd]["fold"], rr["dq"], 11)
        _chk("T3", f"{dd}/ratio", abs(t3[dd]["qw"] / t3[dd]["qo"]), rr["ratio"], 7)
        _chk("T3", f"{dd}/c1", t3[dd]["c1"], rr["c1"], 5)
        _chk("T3", f"{dd}/c2", t3[dd]["c2"], rr["c2"], 5)
    _chk("T3", "ccw/dq", dq_ccw, REC["t3_ccw_dq"], 11)
    _chk("T3", "ccw/fold", fold_ccw, REC["t3_ccw_dq"], 11)
    _chk("T3", "ccw/ratio", abs(q_ccw_w / q_ccw_o), REC["t3_ccw_ratio"], 7)
    _chk("coverage", "R_cov_mm", coverage_radius * 1e3, *REC["r_cov_mm"])
    _chk("coverage", "C3", c3, *REC["c3"])
    tip_calls = [c for c in CALLS if c[1] is not None and c[1] == a1]
    assert len(tip_calls) == 12, len(tip_calls)
    tip_rows = []
    n_call = 0
    for tip in ("top", "bottom"):
        for dd in ("ccw", "P17", "P21"):
            ro, rw = tip_calls[n_call][2], tip_calls[n_call + 1][2]
            n_call += 2
            so, sw = float(ro["averaged_sar_w_per_kg"]), float(rw["averaged_sar_w_per_kg"])
            qq = sw / so
            est = eta[dd] * m_cf[tip][0] + (1.0 - eta[dd])
            rec = REC["tips"][(tip, dd)]
            _chk("tips", f"{tip}/{dd}/sar_o", so, rec[0], 10)
            _chk("tips", f"{tip}/{dd}/sar_w", sw, rec[1], 10)
            _chk("tips", f"{tip}/{dd}/Q", qq, rec[2], 7)
            _chk("tips", f"{tip}/{dd}/mass", float(rw["mass_kg"]) * 1e3, rec[3], 6)
            _chk("tips", f"{tip}/{dd}/eta", eta[dd], rec[4], 6)
            _chk("tips", f"{tip}/{dd}/Mcf", m_cf[tip][0], rec[5], 7)
            _chk("tips", f"{tip}/{dd}/est", est, rec[6], 7)
            _chk("tips", f"{tip}/{dd}/Qest", qq / est, rec[7], 7)
            tip_rows.append((tip, dd, so, sw, qq, float(rw["mass_kg"]) * 1e3, eta[dd],
                             m_cf[tip][0], est, qq / est))
    _chk("readings", "a1g_mm", a1 * 1e3, *REC["a1g_mm"])
    _chk("readings", "s10_o", s10["without"], *REC["s10_o"])
    _chk("readings", "s10_w", s10["with"], *REC["s10_w"])
    _chk("readings", "dsar", dsar, *REC["dsar"])
    _chk("readings", "remesh", r_remesh, *REC["remesh"])
    _chk("readings", "rod_move", r_rod, *REC["rod_move"])
    _chk("readings", "p_rod", pw["p_rod"], *REC["p_rod"])
    _chk("readings", "p_rod/p_ph", pw["p_rod"] / pw["p_phantom"], *REC["p_rod_over_ph"])
    for grp, (cnt, mx) in _CHK["groups"].items():
        P(f"[EX-128] REC {grp}: {cnt} records vs {W1B}, max relative {mx:.3g}")
    for f in _CHK["fail"]:
        P(f"[EX-128] REC FAIL {f}")
    assert not _CHK["fail"], _CHK["fail"][:5]

    # ---- newly asserted negative controls (W1b is the prior, ruling (4)) ----------------
    neg_bar = BLIND_FACTOR * DISCRETE_IDENTITY_RTOL
    for lab, arm in (("without", arm_o), ("with", arm_w)):
        P(f"[EX-128] node NEG [{lab}] {arm['blind']:.5g} >= {neg_bar:.3g} "
          f"({arm['blind'] / neg_bar:.4g}x the bar), ASSERTED")
        assert arm["blind"] >= neg_bar, (lab, arm["blind"])
    P(f"[EX-128] C3 solver-sigma whole ball miss {c3:.5g} > {EXACT_IDENTITY_RTOL:g}, ASSERTED; "
      f"P_rod/P_phantom {pw['p_rod'] / pw['p_phantom']:.5g} (printed)")
    assert c3 > EXACT_IDENTITY_RTOL, c3
    q_all = [r[4] for r in tip_rows]
    P("[EX-128] (Cq) Q on the six tip rows " + ", ".join(f"{x:.5g}" for x in q_all)
      + f"; min {min(q_all):.5g} > 1, ASSERTED")
    assert all(x > 1.0 for x in q_all), q_all

    # ---- DG0 maps ------------------------------------------------------------------------
    own = np.arange(len(cells)) < n_owned
    d_of_c = dof_of_cell[cells]
    v_t = ufl.TestFunction(q0)
    vol = _dg0_vec(ufl.conj(v_t) * ufl.dx(domain=msh))

    def _loss(arm, measure):
        e = arm["e"]["ccw"]
        return _dg0_vec(0.5 * arm["sigma"] * ufl.inner(e, e) * ufl.conj(v_t) * measure)

    loss_o = _loss(arm_o, dx(PHANTOM_CELL_TAG))
    loss_w = _loss(arm_w, dx(PHANTOM_CELL_TAG))
    loss_rod = _loss(arm_w, dx(BIRDCAGE_IMPLANT_CELL_TAG))

    def _owned_sum(vec):
        return float(comm.allreduce(float(np.sum(vec[d_of_c[own]])), op=MPI.SUM))

    for lab, vec, arm in (("(L)_0", loss_o, arm_o), ("(L)_w", loss_w, arm_w)):
        s_l, p_ph = _owned_sum(vec), arm["parts"]["ccw"]["p_phantom"]
        r_l = abs(s_l / p_ph - 1.0)
        P(f"[EX-128] {lab} DG0 phantom-loss map sums to the partition P_phantom: rel {r_l:.3g} "
          f"(band {EXACT_IDENTITY_RTOL:g}, ASSERTED)")
        assert r_l <= EXACT_IDENTITY_RTOL, (lab, r_l)
    s_r = _owned_sum(loss_rod)
    r_r = abs(s_r / pw["p_rod"] - 1.0)
    P(f"[EX-128] (L-rod) DG0 rod-loss map sums to P_rod: rel {r_r:.3g} (band "
      f"{EXACT_IDENTITY_RTOL:g}, ASSERTED)")
    assert r_r <= EXACT_IDENTITY_RTOL, r_r

    p_src_o = arm_o["parts"]["ccw"]["p_src"]
    p_src_w = arm_w["parts"]["ccw"]["p_src"]
    sar_o = np.zeros_like(loss_o)
    sar_w = np.zeros_like(loss_w)
    np.divide(loss_o, PHANTOM_RHO_KG_PER_M3 * vol, out=sar_o,
              where=(vol > 0.0) & (loss_o > 0.0))
    np.divide(loss_w, PHANTOM_RHO_KG_PER_M3 * vol, out=sar_w,
              where=(vol > 0.0) & (loss_w > 0.0))
    sar_o /= p_src_o
    sar_w /= p_src_w
    ratio = np.zeros_like(sar_o)
    np.divide(sar_w, sar_o, out=ratio, where=sar_o > 0.0)
    rod_dens = np.zeros_like(loss_rod)
    np.divide(loss_rod, vol, out=rod_dens, where=(vol > 0.0) & (loss_rod > 0.0))
    rod_dens /= p_src_w

    # ---- (B) masks and the three printed readings (cell-ordered numpy) --------------------
    ratio_c, vol_c, so_c = ratio[d_of_c], vol[d_of_c], sar_o[d_of_c]
    mid_c = mid[cells]
    gel_ok = gelc[cells] & own & (so_c > 0.0)
    rod_c = np.hypot(mid_c[:, 0] - ROD_CENTRE[0], mid_c[:, 1] - ROD_CENTRE[1])
    balls = {}
    for tip in ("top", "bottom"):
        dist = np.linalg.norm(mid_c - np.asarray(TIPS[tip]), axis=1)
        b1 = gelc[cells] & (dist < a1)
        n_b = comm.allreduce(int(np.count_nonzero(b1 & own)), op=MPI.SUM)
        n_b403 = comm.allreduce(int(np.count_nonzero(b1 & own & lead[cells])), op=MPI.SUM)
        P(f"[EX-128] (B) 1 g ball {tip}: {n_b} gel cells, {n_b403} tag-403 cells (ASSERTED > 0, 0)")
        assert n_b > 0 and n_b403 == 0, (tip, n_b, n_b403)
        balls[tip] = b1

    def _reading(label, sel, predicted, holds):
        # item 512: the IN / OUT word is built here, after the reductions, from the predicate
        # `holds(mean, max)`; an f-string argument would be evaluated at the call, before the
        # reading exists (UnboundLocalError on the first call, a stale verdict after).
        sel = sel & gel_ok
        cnt = comm.allreduce(int(np.count_nonzero(sel)), op=MPI.SUM)
        wv = comm.allreduce(float(np.sum(ratio_c[sel] * vol_c[sel])), op=MPI.SUM)
        vv = comm.allreduce(float(np.sum(vol_c[sel])), op=MPI.SUM)
        mx = comm.allreduce(float(np.max(ratio_c[sel])) if np.any(sel) else -np.inf, op=MPI.MAX)
        mean = wv / vv if vv > 0.0 else float("nan")
        verdict = "IN" if holds(mean, mx) else "OUT"
        P(f"[EX-128] reading {label}: {cnt} cells, volume-weighted mean ratio {mean:.5g}, max "
          f"{mx:.5g}; {predicted}: {verdict} (printed, no band)")
        return cnt, mean, mx

    q_top_ccw = tip_rows[0][4]
    for tip in ("top", "bottom"):
        near = np.linalg.norm(mid_c - np.asarray(TIPS[tip]), axis=1) < 1.0e-3
        _reading(f"tip-near {tip} (<1 mm)", near,
                 f"PREDICTED max > Q_top,ccw {q_top_ccw:.5g}",
                 lambda m, x: x > q_top_ccw)
    shaft = (rod_c > 1.0e-3) & (rod_c < 2.0e-3) & (np.abs(mid_c[:, 2]) < 5.0e-3)
    _reading("shaft shadow (1-2 mm from axis, |z| < 5 mm)", shaft,
             "PREDICTED mean in [0, 1)", lambda m, x: 0.0 <= m < 1.0)
    far = rod_c > 30.0e-3
    _reading("far (> 30 mm from axis)", far,
             "PREDICTED mean in [0.97, 1.03]", lambda m, x: 0.97 <= m <= 1.03)

    # ---- outputs (after every assertion) ----------------------------------------------------
    names = {2: "air", PHANTOM_CELL_TAG: "phantom", CONDUCTOR_CELL_TAG: "conductor (hole route: void)",
             BIRDCAGE_IMPLANT_CELL_TAG: "rod sigma = 800 S/m (tag 403)"}
    for sp in base_specs:
        n = int(sp.port_id[1:])
        if sp.port_id in feeds:
            label = f"feed {sp.port_id} ({int(round(90 * feeds.index(sp.port_id)))} deg)"
        else:
            label = "C_in + 0.2 Ohm gap"
        names[PORT_LOWER + n] = label
        names[PORT_UPPER + n] = label
    write_setup_figure(
        msh,
        cell_tags,
        FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names=names,
        hide_tags=(2,),
        translucent_tags=(PHANTOM_CELL_TAG,),
        slice_normal=(0.0, 1.0, 0.0),
        slice_origin=(0.0, 0.0, 0.0),
        title="imp:10 -- sigma = 800 S/m rod in the tuned copper F-human gel phantom, two-port "
              "90 deg drive (IMP-3 step 8)",
        window_size=(1400, 630),
        comm=comm,
    )

    def named(name, cell_values):
        f = fem.Function(q0, name=name)
        f.x.array[:] = 0.0
        f.x.array[d_of_c] = np.real(cell_values)
        f.x.scatter_forward()
        return f

    def named_dof(name, dof_values):
        f = fem.Function(q0, name=name)
        f.x.array[:] = dof_values
        f.x.scatter_forward()
        return f

    export = {
        "SAR_local_without_ccw_per_W": named_dof("SAR_local_without_ccw_per_W", sar_o),
        "SAR_local_with_ccw_per_W": named_dof("SAR_local_with_ccw_per_W", sar_w),
        "SAR_local_with_over_without_ccw": named_dof("SAR_local_with_over_without_ccw", ratio),
        "rod_loss_density_with_ccw_per_W": named_dof("rod_loss_density_with_ccw_per_W", rod_dens),
        "one_gram_ball_top": named("one_gram_ball_top", balls["top"].astype(float)),
        "one_gram_ball_bottom": named("one_gram_ball_bottom", balls["bottom"].astype(float)),
        "tf_probe_ball": named("tf_probe_ball", ballc[cells].astype(float)),
    }
    OUTPUT_DIR.mkdir(exist_ok=True)
    write_xdmf_with_tags(OUTPUT_DIR / f"{BASENAME}_combined", msh, cell_tags, export, comm=comm)
    adopt_host_ownership(OUTPUT_DIR, comm=comm)
    xdmf_name = f"{BASENAME}_combined.xdmf"
    csv_path = OUTPUT_DIR / f"{BASENAME}_tip_rows.csv"
    if rank0:
        with open(csv_path, "w") as fh:
            fh.write("tip,drive,SAR1g_without_W_per_kg,SAR1g_with_W_per_kg,Q,mass_g,eta_z,M_cf,"
                     "estimate,Q_over_estimate\n")
            for r in tip_rows:
                fh.write(f"{r[0]},{r[1]}," + ",".join(f"{x:.6g}" for x in r[2:]) + "\n")
    comm.Barrier()
    P(f"[EX-128] wrote {xdmf_name}, {csv_path.name}, {BASENAME}_setup.png (flagged)")
    P("[EX-128] degree 1, one sigma rod, tuned copper F-human at the rod-free root: SAR values "
      "are not an absolute, a tip-ratio accuracy or a compliance figure")
    P(f"[EX-128] example window {time.perf_counter() - t_w:.1f} s; summed RSS {rss:.2f} GiB")
    P("[EX-128] All assertions hold")


def main() -> None:
    _main(MPI.COMM_WORLD)


if __name__ == "__main__":
    main()
