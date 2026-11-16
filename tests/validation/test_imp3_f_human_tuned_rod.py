"""`IMP-3` step 8 (§9 item 381, §10 I21b): the σ = 800 S/m rod under the tuned
copper F-human's two-port drive.

Brief: ``docs/planning/chunks/IMP-3.md`` "§9 item 381" and the 2026-10-11 03:00
review ruling beneath it (which outranks the brief):

* (1) the reproduction is posed on item 358's own mesh (``repro``: a separate
  ``build_tuned_quadrature(comm, with_z50=False)`` build); the rod-mesh
  rod-as-gel residual against the record is printed as the remesh move;
* (2) ``_build_f_human`` / ``_f_human_kwargs`` are **not edited**: the
  rod-as-gel mapping lives here (``_rod_kwargs`` / ``_build_f_human_rod``).
  The wrapper swaps ``_f_human_kwargs`` in the `WF-7` module namespace for the
  duration of one ``_build_f_human(comm, as_hole=True)`` call (disclosed; the
  original is restored in ``finally``), then replaces the problem's material
  map with tag 403 → the gel and rebuilds the solver on it;
* (3) W1b's node NEG has no same-fixture prior: printed, not asserted;
* (4) SAR balls and coverage use the gel σ masked to tag 3 (σ₃).

Knob ``FEM_EM_IMP3_STEP8`` ∈ {``repro``, ``rod``}; unset ⇒ the gated test
skips; anything else ⇒ ``ValueError``.

Asserted (bands imported by name): repro — (S_rep), (A0), (P_RC), node NEG,
(P_w) ccw / cw; rod — (A0r), (M0), mask disjointness, (P_RC) per arm, (P_w)
per arm, (T3) on P17 and P21, (Σ_P) per arm, coverage power and mass per arm,
controls C1 / C2 wired.  Everything else (tip 1 g Q, η_z, M_cf, the `IMP-2`
estimate, ΔSAR10g(c*), remesh move, P_rod, masses) is printed.

Not claimed: an accurate tip ratio, a PEC / insulated / real lead, ISO 10974,
heating, absolute SAR or compliance, degree 2 or h-convergence, f₀ / Q / match
with the rod in, C16, AED or Palace.
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

STEP8_ENV = "FEM_EM_IMP3_STEP8"
STEP8_VALUES = ("repro", "rod")

#: 359's SAR10g(c*) under ccw (`20261010T171705Z_WF-10-step2s-W1.log:6372`).
SAR10G_C_STAR_359 = 4.852895959e-03
C_STAR_LABEL = "r=45 dz=+0 az17+0"

# ---- predictions, written first (brief), printed only -----------------------
PRED_S_REP = 1e-12
PRED_P_W = 1e-12
PRED_T3 = 5e-12
PRED_SIGMA_P = 1e-14
PRED_COVERAGE = 1e-13
PRED_Q_TOP_CCW = (1.5, 8.0)
PRED_ETA_Z_CCW = (0.2, 1.0)
PRED_D_SAR10G = (1e-4, 3e-2)
PRED_REMESH = (1e-5, 1e-2)
PRED_ROD_S_MOVE = 1e-3
PRED_C1 = (0.10, 2.00)
PRED_C2 = (0.80, 1.05)
PRED_MASS_1G_G = 0.98
PRED_NEG = 7.0e-3


def _say(msg):
    if MPI.COMM_WORLD.rank == 0:
        print(f"[IMP-3 step8] {msg}", flush=True)


def _in(x, band):
    return "IN" if band[0] <= x <= band[1] else "OUT"


def _knob():
    raw = os.environ.get(STEP8_ENV, "").strip()
    if raw == "":
        return None
    if raw in STEP8_VALUES:
        return raw
    raise ValueError(f"{STEP8_ENV}={raw!r}: only unset, 'repro' or 'rod'")


def _rod_kwargs(as_hole=True, implant_rod=None):
    """The F-human generator kwargs with an optional kept rod (ruling (2): a
    wrapper, `WF-7`'s ``_f_human_kwargs`` untouched).  ``None`` ⇒ exactly
    ``_f_human_kwargs(as_hole=as_hole)``; a rod requires ``as_hole=True``."""
    from tests.validation.test_wf7_f_human_b1_quadrature import _f_human_kwargs

    if implant_rod is not None and not as_hole:
        raise ValueError("IMP-3 step 8: implant_rod requires as_hole=True (solid route not gated)")
    kwargs = _f_human_kwargs(as_hole=as_hole)
    if implant_rod is not None:
        kwargs["implant_rod"] = dict(implant_rod)
    return kwargs


def _load():
    from fem_em_solver.utils.constants import MU_0

    from tests.validation.test_wf7_f_human_b1_quadrature import FREQUENCY_HZ

    return -1j * (2.0 * np.pi * float(FREQUENCY_HZ)) * MU_0


def _build_f_human_rod(comm, implant_rod):
    """``_build_f_human(comm, as_hole=True)`` on the kept-rod mesh, with tag 403
    mapped to the gel (the *without* arm) and ``problem`` added to the dict."""
    from fem_em_solver.core import TimeHarmonicSolver
    from fem_em_solver.io.mesh import BIRDCAGE_IMPLANT_CELL_TAG

    from tests.validation import test_wf7_f_human_b1_quadrature as wf7
    from tests.validation.test_port_birdcage_lumped_column import PHANTOM_CELL_TAG

    if implant_rod is None:
        raise ValueError("_build_f_human_rod: a rod is required (use _build_f_human otherwise)")
    orig = wf7._f_human_kwargs

    def _patched(phantom_resolution=None, as_hole=False, emit_leg_ports=None):
        k = orig(phantom_resolution=phantom_resolution, as_hole=as_hole,
                 emit_leg_ports=emit_leg_ports)
        if not as_hole:
            raise ValueError("implant_rod requires as_hole=True")
        k["implant_rod"] = dict(implant_rod)
        return k

    wf7._f_human_kwargs = _patched
    try:
        built = wf7._build_f_human(comm, as_hole=True)
    finally:
        wf7._f_human_kwargs = orig
    ctx = built["ctx"]
    old = ctx.pop("solver")
    problem = old.problem
    del old
    mmap = dict(problem.material_map)
    gel = mmap[PHANTOM_CELL_TAG]
    p_without = dataclasses.replace(
        problem, material_map={**mmap, BIRDCAGE_IMPLANT_CELL_TAG: gel})
    ctx["solver"] = TimeHarmonicSolver(p_without, degree=wf7.DEGREE)
    built["problem"] = p_without
    return built


def test_imp3_s8_coding():
    """Numpy only: knob, kwargs identity, ROD_FH, TIPS, LOAD."""
    from tests.mesh.test_birdcage_implant_rod import ROD, ROD_CENTRE
    from tests.mesh.test_birdcage_implant_rod_f_human import ROD_FH
    from tests.validation.test_imp3_rod_birdcage import TIPS
    from tests.validation.test_wf7_f_human_b1_quadrature import _f_human_kwargs
    from tests.validation.test_wf9_lead_in_coil_phantom import LOAD as WF9_LOAD

    saved = os.environ.get(STEP8_ENV)
    try:
        os.environ.pop(STEP8_ENV, None)
        assert _knob() is None
        for v in STEP8_VALUES:
            os.environ[STEP8_ENV] = v
            assert _knob() == v
        for bad in ("1", "Rod", "both"):
            os.environ[STEP8_ENV] = bad
            with pytest.raises(ValueError):
                _knob()
    finally:
        if saved is None:
            os.environ.pop(STEP8_ENV, None)
        else:
            os.environ[STEP8_ENV] = saved
    for hole in (True, False):
        a = _f_human_kwargs(as_hole=hole)
        b = _rod_kwargs(as_hole=hole, implant_rod=None)
        assert a == b and "implant_rod" not in a and "implant_rod" not in b
    k = _rod_kwargs(as_hole=True, implant_rod=ROD_FH)
    assert k["implant_rod"] == ROD_FH and {kk: v for kk, v in k.items() if kk != "implant_rod"} \
        == _f_human_kwargs(as_hole=True)
    with pytest.raises(ValueError):
        _rod_kwargs(as_hole=False, implant_rod=ROD_FH)
    assert ROD_FH == dict(ROD, keep_cells=True)
    top = tuple(np.asarray(ROD_CENTRE, dtype=float)
                + np.array([0.0, 0.0, ROD["length"] / 2.0 + ROD["radius"]]))
    assert top == tuple(float(x) for x in TIPS["top"]), (top, TIPS["top"])
    load = _load()
    rel = abs(load / WF9_LOAD - 1.0)
    print(f"[IMP-3 step8 coding] LOAD {load:.15e} vs WF-9 {WF9_LOAD:.15e}: rel {rel:.3e} "
          f"(band 1e-15, ASSERTED); TIPS top {TIPS['top']} == ROD_CENTRE + (0,0,l/2+r) {top}")
    assert rel <= 1e-15


# --------------------------------------------------------------------------------
# Lifted (L2): `_ball_mean_closed_form` from `main` at 232a6fe2,
# tests/validation/test_imp3_rod_birdcage.py lines 229–272, with one added
# filter — only owned cells whose tag is PHANTOM_CELL_TAG contribute (the
# kept rod's 403 cells exist here and are excluded the way F-small's hole did).
def _ball_mean_closed_form_tag3(msh, cell_tags, centre, radius):
    import basix

    from fem_em_solver.utils.analytical import pec_prolate_spheroid_field

    from tests.mesh.test_birdcage_implant_rod import ROD_CENTRE, _cell_values_with_ghosts
    from tests.validation.test_imp3_rod_birdcage import (
        CF_CHUNK,
        QDEG_CF,
        SPHEROID_SEMI_MAJOR,
        SPHEROID_SEMI_MINOR,
    )
    from tests.validation.test_port_birdcage_lumped_column import PHANTOM_CELL_TAG

    comm = msh.comm
    tdim = msh.topology.dim
    n_owned = msh.topology.index_map(tdim).size_local
    dmap = msh.geometry.dofmap
    if dmap.shape[1] != 4:
        raise RuntimeError(f"expected affine tets, geometry dofmap {dmap.shape}")
    xg = msh.geometry.x
    c = np.asarray(centre, dtype=float)
    rc = np.asarray(ROD_CENTRE, dtype=float)
    pts_ref, wts = basix.make_quadrature(basix.CellType.tetrahedron, QDEG_CF)
    tag = _cell_values_with_ghosts(msh, cell_tags)
    cells = np.arange(n_owned, dtype=np.int32)
    cells = cells[tag[:n_owned] == PHANTOM_CELL_TAG]
    vmin = np.min(np.linalg.norm(xg[dmap[cells]] - c, axis=2), axis=1)
    cells = cells[vmin < radius + 0.01]
    num = 0.0
    vol = 0.0
    for start in range(0, cells.size, CF_CHUNK):
        chunk = cells[start:start + CF_CHUNK]
        v = xg[dmap[chunk]]
        jac = np.stack([v[:, 1] - v[:, 0], v[:, 2] - v[:, 0], v[:, 3] - v[:, 0]], axis=2)
        det = np.abs(np.linalg.det(jac))
        phys = v[:, 0][:, None, :] + np.einsum("cij,qj->cqi", jac, pts_ref)
        chi = (np.sum((phys - c) ** 2, axis=2) < radius**2).astype(float)
        if not np.any(chi):
            continue
        e = pec_prolate_spheroid_field(
            phys.reshape(-1, 3) - rc,
            semi_major=SPHEROID_SEMI_MAJOR,
            semi_minor=SPHEROID_SEMI_MINOR,
            e0=1.0,
            interior="continue",
        ).reshape(len(chunk), len(wts), 3)
        e2 = np.sum(np.abs(e) ** 2, axis=2)
        num += float(np.sum(det[:, None] * wts[None, :] * chi * e2))
        vol += float(np.sum(det[:, None] * wts[None, :] * chi))
    num = comm.allreduce(num, op=MPI.SUM)
    vol = comm.allreduce(vol, op=MPI.SUM)
    return num / vol, vol


def _read_record(comm):
    from tests.validation.test_wf10_f_human_tuned_b1_sar import C_IN_F_HUMAN_PF, NODE_JSON

    assert os.path.exists(NODE_JSON), f"{NODE_JSON} absent: item 381 BLOCKED (never re-run 358)"
    rec = None
    if comm.rank == 0:
        with open(NODE_JSON) as fh:
            rec = json.load(fh)
    rec = comm.bcast(rec, root=0)
    assert abs(rec["c_in_pf"] - C_IN_F_HUMAN_PF) == 0.0, rec["c_in_pf"]
    return np.array([[complex(a, b) for a, b in row] for row in rec["s_node"]],
                    dtype=np.complex128)


def _c_star(az17, z_c):
    from tests.validation.test_wf10_f_human_tuned_sar10g import _lattice

    for label, *_x, c in _lattice(az17, z_c):
        if label == C_STAR_LABEL:
            return tuple(float(v) for v in c)
    raise KeyError(C_STAR_LABEL)


def _run_repro(comm):
    from fem_em_solver.post.sar import averaging_ball_radius, build_density_field, mass_averaged_sar

    from tests.mesh.helpers import global_cell_tag_set
    from tests.mesh.test_birdcage_port_sheet_prerequisite import CELL_COUNT_BAND
    from tests.validation.test_birdcage_b1_plus_map import PHANTOM_RHO_KG_PER_M3
    from tests.validation.test_birdcage_power_identity import DISCRETE_IDENTITY_RTOL
    from tests.validation.test_mass_averaged_sar_standard_masses import (
        QUADRATURE_DEGREE,
        TEN_GRAM_KG,
    )
    from tests.validation.test_port_birdcage_lumped_column import PHANTOM_CELL_TAG
    from tests.validation.test_port_birdcage_ring_matrix import REUSE_REPRODUCTION_RTOL
    from tests.validation.test_port_circuit_layer_field import _residual
    from tests.validation.test_wf7_f_human_sar10g_hotspot import _phantom_geometry
    from tests.validation.test_wf10_f_human_tuned_b1_sar import build_tuned_quadrature
    from tests.validation.test_wf12_f_human_copper_matrix import BLIND_FACTOR, HOLE_CELL_RECORD
    from tests.validation.test_wf12_f_human_hole import _rss_gib

    t_w = time.perf_counter()
    s_rec = _read_record(comm)
    out = build_tuned_quadrature(comm, with_z50=False)
    s_node = np.asarray(out["s_node"], dtype=np.complex128)
    r_srep = _residual(s_node, s_rec)
    a = out["anchors"]
    cells = out["cells"]
    cell_rel = cells / HOLE_CELL_RECORD - 1.0
    _say(f"[repro] (S_rep) s_node vs NODE_JSON read before the builder: residual {r_srep:.6e} "
         f"(band {REUSE_REPRODUCTION_RTOL:g}, ASSERTED; PREDICTED <= {PRED_S_REP:g})")
    _say(f"[repro] (A0) cells {cells} vs {HOLE_CELL_RECORD} rel {cell_rel:.3e} (band "
         f"{CELL_COUNT_BAND:g}, ASSERTED; census asserted clean by the builder)")
    _say(f"[repro] (P_RC) {a['p_rc_node']:.3e} (band {DISCRETE_IDENTITY_RTOL:g}, ASSERTED); node "
         f"NEG {a['blind_node']:.4e} = {a['blind_node'] / DISCRETE_IDENTITY_RTOL:.1f}x (bar "
         f"{BLIND_FACTOR:.0f}x, ASSERTED; 358 read 7.0125e-3)")
    for s, (rel, blind) in a["p_w"].items():
        _say(f"[repro] (P_w) [{s}] rel {rel:.3e} (band {DISCRETE_IDENTITY_RTOL:g}, ASSERTED); "
             f"blind {blind:.4e}; (Lin) {a['lin'][s][0]:.3e} (printed)")
    # Reading: SAR10g(c*) under ccw on this mesh with the builder's σ (as 359).
    ctx = out["ctx"]
    msh, cell_tags = ctx["msh"], ctx["cell_tags"]
    centre, _r, _h, _lo, _hi = _phantom_geometry(msh, cell_tags, PHANTOM_CELL_TAG, comm)
    c = _c_star(float(out["az_feeds"][0]), float(centre[2]))
    all_tags = global_cell_tag_set(msh, cell_tags)
    rho = build_density_field(
        msh, PHANTOM_RHO_KG_PER_M3, cell_tags=cell_tags,
        density_map={int(t): 0.0 for t in sorted(all_tags) if int(t) != PHANTOM_CELL_TAG})
    a10 = averaging_ball_radius(mass_kg=TEN_GRAM_KG, rho=PHANTOM_RHO_KG_PER_M3)
    r = mass_averaged_sar(out["e"]["ccw"], sigma=out["sigma"], rho=rho, center=c, radius=a10,
                          comm=comm, quadrature_degree=QUADRATURE_DEGREE)
    sar = float(r["averaged_sar_w_per_kg"])
    _say(f"[repro] SAR10g(c* {C_STAR_LABEL} at {tuple(round(x * 1e3, 3) for x in c)} mm) ccw "
         f"{sar:.9e} W/kg vs 359 {SAR10G_C_STAR_359:.9e}: rel {sar / SAR10G_C_STAR_359 - 1.0:.3e} "
         f"(PREDICTED <= 1e-9, reading)")
    _m, _p, rss = _rss_gib(comm)
    _say(f"[repro] times {out['times']}; window {time.perf_counter() - t_w:.1f} s; summed RSS "
         f"{rss:.2f} GiB")
    assert r_srep <= REUSE_REPRODUCTION_RTOL, f"(S_rep) {r_srep:.3e}"
    assert abs(cell_rel) <= CELL_COUNT_BAND, cell_rel
    assert a["p_rc_node"] <= DISCRETE_IDENTITY_RTOL
    assert a["blind_node"] >= BLIND_FACTOR * DISCRETE_IDENTITY_RTOL
    for s, (rel, _b) in a["p_w"].items():
        assert rel <= DISCRETE_IDENTITY_RTOL, (s, rel)


def _run_rod(comm):
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


@complex_only
def test_imp3_s8_tuned_rod():
    knob = _knob()
    if knob is None:
        pytest.skip(f"{STEP8_ENV} unset: heavy -n 8 F-human windows, environment-selected")
    comm = MPI.COMM_WORLD
    if knob == "repro":
        _run_repro(comm)
    else:
        _run_rod(comm)
