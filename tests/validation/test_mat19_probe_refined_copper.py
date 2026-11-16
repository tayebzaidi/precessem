"""`MAT-19` step 3 re-pose / P10′ -- the probe on a phantom-refined tuned copper fixture (§9 item 314).

Brief: ``docs/planning/chunks/MAT-19.md`` (§9 item 314, chain P10′), with the
2026-10-08 03:00 review rulings.  Attempt 2 of 4 on `MAT-19` step 3 (attempt 1,
item 265, is ``test_mat19_probe_tuned_copper.py``, unedited; it is the template).

Fixture: `ANS-7`'s exact-rectangle hole mesh (κ off) with the phantom at
``PHANTOM_RES`` = 0.005 m, degree 2, Leontovich copper.  The phantom knob reaches
the mesh through a LOCAL copy of ``_exact_geo`` (``_geo_p``) and a LOCAL copy of
``_build_ladder``'s Cu 4×4 sweep (``_z4_cu``) -- ruling (1): neither
``test_th14_birdcage_copper.py`` nor ``test_ans7_exact_strip.py`` is edited.  The
local 4×4 runs the Cu column only (no PEC column, no field solve: §9 item 310's
flags on the h rung, disclosed; the PEC column feeds no anchor here).

Stages (opt-in ``FEM_EM_MAT19_STEP3P``; unset ⇒ skip, anything else raises):
  ``build``  (W0, -n 8): default-resolution local build cells == ``_build_ladder``'s
             (ruling 1, build only); refined cells / tag-3 cells printed;
             tag-3 > 3 000 asserted; cells > 130 000 is a RESOURCE STOP (ruling 2).
  ``cin``    (W1, -n 8): Cu 4×4 at 64 MHz, C_mode (TH-17 construction), C_in by
             ``_c_in`` (secant, |g| <= 1e-15 S); writes ``CIN_JSON`` (mode "x").
  ``series`` (W2 -n 8 / W3 -n 4): reads ``CIN_JSON``; at -n 8 three Route A drives
             at C_in (63.9 / 64.0 / 64.1 MHz), in-model zero printed; at the 64.0
             node (every width) (AB′) and (P), G = 1 W / P_acc, P1 thermal 40
             steps, the probe series (item 258 / 265 construction); -n 8 writes
             ``NPZ``, other widths check (W) against it.

Asserted: (A0) tie on the default build, tag-3 > 3 000, snapshot 0 = 0.0, (M);
(Φ); (Src) (AB′) <= ``REDUCTION_BAND`` and (P) <= ``DISCRETE_IDENTITY_RTOL``;
(A) <= ``B_BAND``; (U4) bitwise; (W) <= ``WIDTH_BAND``; (Σ) <=
``C4_COVARIANCE_BAND``.  Negative controls (rule (e), predicted factors, miss
only): (U3c) iff G_w > 0, (C1″), (f″).
"""

from __future__ import annotations

import dataclasses
import hashlib
import json
import os
import resource
import time
from pathlib import Path

import numpy as np
import pytest
import ufl
from dolfinx import default_scalar_type, fem
from mpi4py import MPI

from fem_em_solver.core.cavity import surface_resistance_ohm
from fem_em_solver.ports.sparameters import run_n_port_sparameter_sweep
from fem_em_solver.thermal import probe_series
from fem_em_solver.thermal import probe as _probe_mod
from fem_em_solver.thermal.probe import _ball_average_on_cells, _lattice

from tests.complex_mode import complex_only
from tests.mesh.test_birdcage_leg_offset import _sheet_azimuth_deg
from tests.mesh.test_birdcage_port_sheets import SHEET_IFACE
from tests.mesh.test_birdcage_port_tags import PHANTOM_HEIGHT
from tests.validation import test_th14_birdcage_copper as th14
from tests.validation.test_ans7_exact_strip import (
    GATED_WIDTH_FRACTION,
    REDUCTION_BAND,
    _assert_tie,
    _c_in,
    _c_mode,
    _fixture_tie,
    _route_a_solve,
    _s11_from_z,
    parallel_rows,
    route_b_z,
    z_from_s,
)
from tests.validation.test_ans7_sparse_route_b import _ANS7
from tests.validation.test_birdcage_power_identity import (
    DISCRETE_IDENTITY_RTOL,
    _sheet_field_dissipation_w,
    _source_power_w,
)
from tests.validation.test_mat10_sar_to_heat_sphere import _robin_tags, _to_sub
from tests.validation.test_mat15_coil_driven_heating_phantom import (
    A,
    B_BAND,
    C4_COVARIANCE_BAND,
    DT_1,
    F_CONTROL_FLOOR,
    N_MAIN,
    RHO_C,
    T_INF,
    T_STAR,
    _allsum,
    _be_one,
    _cell_average_dg0_tag,
    _say,
    _tag_submesh,
    _transient,
)
from tests.validation.test_mat15_tuned_copper_heating import PHI_BAND
from tests.validation.test_mat19_probe_coil_phantom import (
    A_PROBE,
    _lattice_q_mean,
    _probe_set,
)
from tests.validation.test_mat19_probe_tuned_copper import (
    MEAN_THETA_236_K,
    STEP2_N40_PER_DEP_W,
    _point_values,
)
from tests.validation.test_mat19_virtual_temperature_probe import WIDTH_BAND
from tests.validation.test_port_lumped_sheet_sweep import RECIPROCITY_BAND

KNOB = "FEM_EM_MAT19_STEP3P"
STAGES = ("build", "cin", "series")
PHANTOM_RES = 0.005
DRIVE_HZ = 64.0e6
TUNE_MHZ = (63.9, 64.0, 64.1)
# Resource stops (ruling 2) -- never physics pass/fail.
CELL_CAP = 130_000
DRIVE_CLOCK_GATE_S = 240.0
TAG3_FLOOR = 3_000  # (A0) precondition: the phantom *is* refined (item 265 had 537)
C_IN_X1_PF = 63.590091  # the x1 C_in (T21′), printed only
CIN_JSON = "/workspace/logs/mat19/step3p-cin.json"
NPZ = "/workspace/logs/mat19/step3p-n8.npz"
NPZ_WIDTH = 8
DT = T_STAR / N_MAIN


def _stage():
    raw = os.environ.get(KNOB, "").strip()
    if not raw:
        return None
    if raw not in STAGES:
        raise ValueError(f"{KNOB} must be unset or one of {STAGES}, got {raw!r}")
    return raw


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _geo_p(comm, phantom_res):
    """LOCAL copy of ``_exact_geo``'s body (``test_ans7_exact_strip.py:409–425``)
    with ``phantom_resolution`` forwarded to ``_sheets_build`` (``None`` forwards
    nothing, so the default build is ``_exact_geo``'s)."""
    extra = {"leg_sheet_strip_fraction": GATED_WIDTH_FRACTION}
    if phantom_res is not None:
        extra["phantom_resolution"] = float(phantom_res)
    msh, cell_tags, facet_tags, _d, _t = th14._sheets_build(True, as_hole=True, **extra)
    tdim = msh.topology.dim
    msh.topology.create_connectivity(tdim - 1, tdim)
    msh.topology.create_entity_permutations()
    cells = int(msh.topology.index_map(tdim).size_global)
    outer_tags, _census = th14._outer_box_tags(msh, facet_tags, comm)
    tags_f, port_defs, specs = th14._build_ports(msh, cell_tags, comm,
                                                 width_correction_kappa=None)
    n_own = msh.topology.index_map(tdim).size_local
    own = cell_tags.indices < n_own
    tag3 = int(comm.allreduce(int(np.count_nonzero(cell_tags.values[own] == th14.PHANTOM_CELL_TAG)),
                              op=MPI.SUM))
    return (msh, cell_tags, facet_tags, outer_tags, tags_f, port_defs), specs, cells, tag3


def _z4_cu(geo, specs, f_hz, degree):
    """LOCAL copy of ``_build_ladder``'s Cu sigma sweep (``test_th14_birdcage_copper.py:345–356``),
    one frequency, pec_column=False, field_solve=False."""
    msh, cell_tags, facet_tags, outer_tags, tags_f, port_defs = geo
    omega = 2.0 * np.pi * f_hz
    imp_problem = th14._problem(msh, cell_tags, f_hz, outer_tags, (th14.OUTER_BOX_TAG,))
    z_s = (1.0 + 1.0j) * surface_resistance_ohm(omega, th14.COPPER_SIGMA)
    term = th14._leontovich_term(msh, facet_tags, z_s, omega)
    return th14._sweep_record(run_n_port_sparameter_sweep(
        imp_problem, port_defs, lumped_sheet_ports=specs, lumped_sheet_facet_tags=tags_f,
        extra_bilinear_terms=[term], degree=degree))


def _build_stage(comm):
    t0 = time.perf_counter()
    lad = th14._build_ladder(2, sigmas=(th14.COPPER_SIGMA,), frequencies_hz=[],
                             width_correction_kappa=None, measure_kappa=False,
                             pec_column=False, field_solve=False,
                             leg_sheet_strip_fraction=GATED_WIDTH_FRACTION)
    t_lad = time.perf_counter() - t0
    t1 = time.perf_counter()
    _g0, _s0, cells0, tag3_0 = _geo_p(comm, None)
    del _g0, _s0
    t_def = time.perf_counter() - t1
    t2 = time.perf_counter()
    _g, _s, cells, tag3 = _geo_p(comm, PHANTOM_RES)
    t_ref = time.perf_counter() - t2
    _say(f"\nMAT19-3p (W0) _build_ladder cells {lad['cells']} ({t_lad:.1f} s); local default build "
         f"cells {cells0} tag-3 {tag3_0} ({t_def:.1f} s); refined (phantom {PHANTOM_RES} m) cells "
         f"{cells} tag-3 {tag3} ({t_ref:.1f} s) (PREDICTED tag-3 6 000-10 000, cells 90 000-115 000; "
         f"resource stop > {CELL_CAP}); -n {comm.size}")
    return dict(stage="build", lad_cells=int(lad["cells"]), cells0=cells0, tag3_0=tag3_0,
                cells=cells, tag3=tag3, tie=_fixture_tie())


def _cin_stage(comm):
    t0 = time.perf_counter()
    geo, specs, cells, tag3 = _geo_p(comm, PHANTOM_RES)
    t_mesh = time.perf_counter() - t0
    t1 = time.perf_counter()
    sw = _z4_cu(geo, specs, DRIVE_HZ, 2)
    comm.Barrier()
    t_sweep = time.perf_counter() - t1
    z64 = z_from_s(np.asarray(sw["s"], dtype=np.complex128))
    omega = 2.0 * np.pi * DRIVE_HZ
    cm = _c_mode(z64, omega)
    c_in, g_in, it = _c_in(z64, DRIVE_HZ, cm["c_mode"])
    sig_max = float(np.max(sw["sigma"]))
    _say(f"\nMAT19-3p (W1) cells {cells} tag-3 {tag3}; mesh {t_mesh:.1f} s, Cu 4x4 p2 sweep "
         f"{t_sweep:.1f} s; recip {sw['reciprocity']:.3e} (band {RECIPROCITY_BAND:.0e}); sigma_max "
         f"{sig_max:.12f}; off-diag {cm['off']:.3e}; z1 {cm['z1']:.12e}; |z1-z3|/|z1| {cm['split']:.3e}")
    _say(f"MAT19-3p (W1) C_mode,Cu = {cm['c_mode']!r} F; C_in = {c_in!r} F = {c_in * 1e12:.6f} pF "
         f"(|g| {abs(g_in):.3e} S after {it} it); C_in/C_in(x1 {C_IN_X1_PF} pF) - 1 = "
         f"{c_in * 1e12 / C_IN_X1_PF - 1.0:+.4e} (PREDICTED within ~1e-3, printed)")
    if comm.rank == 0:
        os.makedirs(os.path.dirname(CIN_JSON), exist_ok=True)
        with open(CIN_JSON, "x") as fh:
            json.dump({"c_in": c_in, "c_mode": cm["c_mode"], "g_in": g_in, "iter": it,
                       "z64_re": np.real(z64).tolist(), "z64_im": np.imag(z64).tolist(),
                       "cells": cells, "tag3": tag3, "degree": 2, "phantom_res": PHANTOM_RES,
                       "f_hz": DRIVE_HZ, "ranks": comm.size, "reciprocity": sw["reciprocity"],
                       "t_sweep_s": t_sweep, "t_mesh_s": t_mesh}, fh)
        print(f"MAT19-3p (W1) wrote {CIN_JSON}", flush=True)
    comm.Barrier()
    return dict(stage="cin", cells=cells, tag3=tag3, g_in=g_in, recip=sw["reciprocity"],
                sig_max=sig_max)


def _series_stage(comm):  # noqa: C901
    if not os.path.exists(CIN_JSON):
        pytest.fail(f"{CIN_JSON} absent: run W1 (stage cin) first")
    cj = json.loads(Path(CIN_JSON).read_text())
    c_in = float(cj["c_in"])
    z64 = np.array(cj["z64_re"]) + 1j * np.array(cj["z64_im"])
    clock = {}
    t0 = time.perf_counter()
    geo, specs, cells, tag3 = _geo_p(comm, PHANTOM_RES)
    msh, cell_tags, facet_tags, _outer, tags_f, _pdefs = geo
    if cells != cj["cells"] or tag3 != cj["tag3"]:
        pytest.fail(f"mesh differs from W1's: cells {cells}/{cj['cells']} tag-3 {tag3}/{cj['tag3']}")
    phi1 = _sheet_azimuth_deg(msh, tags_f, SHEET_IFACE + 1, comm)
    sin_phi = abs(np.sin(np.radians(phi1)))
    clock["mesh"] = time.perf_counter() - t0
    _say(f"\nMAT19-3p C_in = {c_in!r} F (from {CIN_JSON}); cells {cells} tag-3 {tag3}; (Phi) phi_1 = "
         f"{phi1:.12f} deg |sin| = {sin_phi:.3e} (bound {PHI_BAND:.0e}); -n {comm.size}")

    drives = {}
    fields = None
    mhz_list = TUNE_MHZ if comm.size == NPZ_WIDTH else (64.0,)
    for mhz in mhz_list:
        f = mhz * 1e6
        t1 = time.perf_counter()
        s_a, fl = _route_a_solve(th14, geo, specs, f, "cu", c_f=c_in,
                                 return_fields=(abs(mhz - 64.0) < 1e-9), degree=2)
        comm.Barrier()
        dt_d = time.perf_counter() - t1
        z_in = _ANS7.add_c1(_ANS7._z_from_s11(s_a), 2.0 * np.pi * f, c_in)
        drives[mhz] = dict(s_a=s_a, z_in=z_in, im_y=float((1.0 / complex(z_in)).imag), t=dt_d)
        if fl is not None:
            fields = fl
        _say(f"MAT19-3p drive {mhz:.1f} MHz at C_in: S11_A {complex(s_a):.12e}; Im Y_in,A "
             f"{drives[mhz]['im_y']:+.9e} S; drive {dt_d:.1f} s")
    clock["EM drive 64.0"] = drives[64.0]["t"]
    s_b64 = _s11_from_z(route_b_z(z64, DRIVE_HZ, c_in)[1])
    ab = float(abs(drives[64.0]["s_a"] - s_b64) / abs(s_b64))
    _say(f"MAT19-3p (Src/AB') circuit (solved Z64, C_in) vs field at 64.0 MHz: S11_B "
         f"{complex(s_b64):.12e}; |S_A - S_B|/|S_B| = {ab:.3e} (band {REDUCTION_BAND:.0e})")
    f0_rel = None
    if len(drives) == 3:
        rows = parallel_rows(np.array(TUNE_MHZ) * 1e6, np.array([drives[m]["z_in"] for m in TUNE_MHZ]))
        f0 = rows.get("f0_hz")
        f0_rel = (f0 / DRIVE_HZ - 1.0) if f0 else float("nan")
        _say(f"MAT19-3p in-model zero: up-crossings {rows['crossings_up_hz']}; f0,A/64 MHz - 1 = "
             f"{f0_rel:+.6e}; grid resolution (the 64.0 node's AB' residual) {ab:.3e} -> "
             f"|f0 rel|/res = {abs(f0_rel) / ab:.3f} (printed, never asserted)")

    t2 = time.perf_counter()
    omega = 2.0 * np.pi * DRIVE_HZ
    z_c = _ANS7.z_cap(omega, c_in)
    sp_list = [dataclasses.replace(s, port_impedance_ohm=z_c)
               if i in _ANS7.TUNING_TERMINATED_INDICES else s for i, s in enumerate(specs)]
    e = fields.e_complex
    sigma_f = fields.sigma_field
    dx = ufl.Measure("dx", domain=msh, subdomain_data=cell_tags)

    def _vol(measure):
        form = fem.form(0.5 * sigma_f * ufl.inner(e, e) * measure)
        return float(np.real(comm.allreduce(fem.assemble_scalar(form), op=MPI.SUM)))

    acct = {"mesh": msh, "facet_tags": tags_f}
    sheet_objs = {sp.port_id: sp.sheet(driven=(sp.port_id == th14.DRIVEN)) for sp in sp_list}
    p_sheets = sum(float(_sheet_field_dissipation_w(acct, sh, e, omega))
                   for sh in sheet_objs.values())
    p_src = float(_source_power_w(acct, sheet_objs[th14.DRIVEN], e, omega))
    p_total = _vol(dx)
    p_phantom = _vol(dx(th14.PHANTOM_CELL_TAG))
    z_s = (1.0 + 1.0j) * surface_resistance_ohm(omega, th14.COPPER_SIGMA)
    p_surf = th14._surface_loss_w(msh, facet_tags, z_s, e, comm)
    p_non = p_total - p_phantom
    p_rel = abs(p_src - p_sheets - p_phantom - p_non - p_surf) / abs(p_src)
    p_acc = p_src - p_sheets
    G = 1.0 / p_acc
    eta = p_phantom / p_acc
    _say(f"MAT19-3p (Src/P) identity rel {p_rel:.3e} (bound {DISCRETE_IDENTITY_RTOL:.0e}); P_acc = "
         f"{p_acc:.15e} W  G = {G:.12e} 1/W  eta = {eta:.9f} (item 236's 0.5777 on x1, reference)")

    sub, parent_ids = _tag_submesh(msh, cell_tags, th14.PHANTOM_CELL_TAG)
    ft = _robin_tags(sub)
    n_sub = sub.topology.index_map(3).size_global
    one = fem.Constant(sub, default_scalar_type(1.0))
    v_sub = float(np.real(_allsum(comm, one * ufl.dx(domain=sub))))
    Wp, q_par = _cell_average_dg0_tag(msh, cell_tags, th14.PHANTOM_CELL_TAG,
                                      G * 0.5 * sigma_f * ufl.inner(e, e))
    q_h = _to_sub(Wp, q_par, sub, parent_ids)
    clock["power+source transfer"] = time.perf_counter() - t2

    t3 = time.perf_counter()
    res = _transient(sub, ft, q_h, dt=DT, n_steps=N_MAIN, theta=0.5, startup_half_steps=2,
                     snapshot_steps=range(N_MAIN + 1))
    res_be = _be_one(sub, ft, q_h)
    clock["transient + BE"] = time.perf_counter() - t3
    V = res.snapshots[N_MAIN].function_space
    th = {}
    for n in range(N_MAIN + 1):
        f_ = fem.Function(V)
        f_.x.array[:] = res.snapshots[n].x.array - T_INF
        f_.x.scatter_forward()
        th[n] = f_
    th_be = fem.Function(V)
    th_be.x.array[:] = res_be.temperature.x.array - T_INF
    th_be.x.scatter_forward()

    names, order, probes, centres = _probe_set()
    t4 = time.perf_counter()
    lam = _lattice(3)
    s = probe_series(th, DT, probes)
    s_be = probe_series({1: th_be}, DT_1, probes)
    qbar = {k: _lattice_q_mean(V, q_h, p, lam)[0] for k, p in zip(order, probes)}
    point = _point_values(th[N_MAIN], centres, names)
    clock["probe series"] = time.perf_counter() - t4
    means = np.array([float(np.real(_allsum(comm, th[n] * ufl.dx(domain=sub)))) / v_sub
                      for n in range(N_MAIN + 1)])
    mean_be = float(np.real(_allsum(comm, th_be * ufl.dx(domain=sub)))) / v_sub

    imap = sub.topology.index_map(3)
    part = comm.allgather((int(imap.size_local), int(imap.num_ghosts)))
    o_w = sum(1 for so, _g in part if so > 0)
    g_w = sum(g for _so, g in part)
    _say(f"MAT19-3p (part) -n {comm.size}: sub owned {[p[0] for p in part]} ghosts "
         f"{[p[1] for p in part]} O_w {o_w} G_w {g_w}")
    ghost = None
    if comm.size > 1:
        allc = np.arange(imap.size_local + imap.num_ghosts, dtype=np.int32)
        gval, gvol = _ball_average_on_cells(th[N_MAIN], probes[order.index("Omega")], allc)
        ghost = (abs(gvol / v_sub - 1.0), abs(gval / means[N_MAIN] - 1.0))
    allrows = comm.allgather(s.values[N_MAIN].copy())
    u4 = all(np.array_equal(allrows[0], r) for r in allrows)
    rss = comm.allreduce(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, op=MPI.SUM)
    hashes = (_sha(Path(__file__).resolve()), _sha(_probe_mod.__file__))
    clock["window"] = time.perf_counter() - t0
    return dict(stage="series", cells=cells, tag3=tag3, sin_phi=sin_phi, ab=ab, p_rel=p_rel,
                f0_rel=f0_rel, eta=eta, G=G, n_sub=n_sub, v_sub=v_sub, s=s, s_be=s_be,
                qbar=qbar, point=point, means=means, mean_be=mean_be, ghost=ghost, u4=u4,
                names=names, order=order, centres=centres, hashes=hashes, clock=clock,
                o_w=o_w, g_w=g_w, rss_kib=rss, drive_t=drives[64.0]["t"])


@pytest.fixture(scope="module")
def run():
    st = _stage()
    if st is None:
        pytest.skip(f"{KNOB} unset (heavy windows)")
    comm = MPI.COMM_WORLD
    out = {"build": _build_stage, "cin": _cin_stage, "series": _series_stage}[st](comm)
    out["comm"] = comm
    return out


@complex_only
def test_mat19_step3p(run):  # noqa: C901
    comm = run["comm"]
    if run["stage"] == "build":
        assert run["cells0"] == run["lad_cells"], (
            f"(W0) local default build {run['cells0']} != _build_ladder {run['lad_cells']}")
        _assert_tie(run["cells0"], run["tie"])
        assert run["tag3"] > TAG3_FLOOR, f"(A0) tag-3 {run['tag3']} <= {TAG3_FLOOR}"
        if run["cells"] > CELL_CAP:
            pytest.fail(f"RESOURCE STOP (ruling 2): refined cells {run['cells']} > {CELL_CAP}; "
                        f"park before W1")
        return
    if run["stage"] == "cin":
        assert abs(run["g_in"]) <= 1e-15, f"(C_in) |g| {abs(run['g_in']):.3e}"
        assert run["recip"] <= RECIPROCITY_BAND, run["recip"]
        assert run["tag3"] > TAG3_FLOOR
        return

    s, s_be = run["s"], run["s_be"]
    order, names, centres = run["order"], run["names"], run["centres"]
    v = s.values
    j = {k: order.index(k) for k in order}
    om = "Omega"
    n_sub = run["n_sub"]
    _say(f"\nMAT19-3p width -n {comm.size}; dt = {DT:.4f} s; t* = {T_STAR:.3f} s; sub cells {n_sub} "
         f"(tag-3 {run['tag3']})")
    _say("MAT19-3p ball V_h/(4 pi a^3/3) - 1: " + "  ".join(
        f"{k} {s.volumes[j[k]] / s.nominal_volumes[j[k]] - 1.0:+.3e}" for k in order if k != om))
    m_ok = all(np.hypot(centres[k][0], centres[k][1]) + A_PROBE <= A
               and abs(centres[k][2]) + A_PROBE <= PHANTOM_HEIGHT / 2.0
               for k in order if k != om)
    PO = v[:, j[om]]
    relA = np.abs(PO[1:] - run["means"][1:]) / np.abs(run["means"][1:])
    relA_be = abs(s_be.values[0, j[om]] - run["mean_be"]) / abs(run["mean_be"])
    volA = abs(s.volumes[j[om]] / run["v_sub"] - 1.0)
    _say(f"MAT19-3p (A) P_Omega vs int theta/V_sub: max n=1..40 {np.max(relA):.3e}, BE "
         f"{relA_be:.3e}; V_Omega/V_sub - 1 = {volA:.3e} (bound {B_BAND:.0e}); (M) {m_ok}")
    _say(f"MAT19-3p <theta>(t*) per accepted watt = {run['means'][N_MAIN]:.6f} K (item 236's "
         f"{MEAN_THETA_236_K} on x1, reference only); eta = {run['eta']:.9f}")
    _say("MAT19-3p series theta_bar [K per accepted W] (n: " + " ".join(order) + "):")
    for n in (1, 10, 20, 30, 40):
        _say(f"MAT19-3p   n={n:2d} t={n * DT:7.2f}  " + " ".join(f"{v[n, j[k]]:.6e}" for k in order))
    _say("MAT19-3p (P4') ball(t*)/point(t*) - 1: " + "  ".join(
        f"{k} {v[N_MAIN, j[k]] / run['point'][k][0] - 1.0:+.4e}" for k in names)
        + " (predicted |.| <= 2 %, printed)")
    keys_r = names + ["C'", om]
    _say("MAT19-3p ratio (i) [MIXED: 64 MHz per ACCEPTED W / 10 MHz per DEPOSITED W, step 2 "
         "20261006T005727Z_MAT-19-step2-a2-W1.log:1970]: " + "  ".join(
             f"{k} {v[N_MAIN, j[k]] / STEP2_N40_PER_DEP_W[k]:.6f}" for k in keys_r))
    _say("MAT19-3p ratio (ii) [per DEPOSITED W at both; 64 MHz rise / eta]: " + "  ".join(
        f"{k} {v[N_MAIN, j[k]] / run['eta'] / STEP2_N40_PER_DEP_W[k]:.6f}" for k in keys_r))
    c1 = abs(v[N_MAIN, j["B+"]] - v[N_MAIN, j["B"]]) / v[N_MAIN, j["B"]]
    _say(f"MAT19-3p (C1'') |B+ - B|/B = {c1:.4e} (predicted [1e-3, 1e-2]; asserted > {WIDTH_BAND:.0e})")
    sig = abs(v[N_MAIN, j["C"]] - v[N_MAIN, j["C'"]]) / v[N_MAIN, j["C"]]
    cp = "C'"
    qc, qcp = run["qbar"]["C"], run["qbar"][cp]
    _say(f"MAT19-3p (Sigma) |C - C'|/C at t* = {sig:.4e} (bound {C4_COVARIANCE_BAND:.0e}; predicted "
         f"[2e-3, 3e-2]); q_h ball pair C {qc:.6e} C' {qcp:.6e} rel {abs(qc - qcp) / abs(qc):.4e} "
         f"(predicted <= 5 %); ball volumes C {s.volumes[j['C']]:.6e} C' {s.volumes[j[cp]]:.6e}")
    ff = abs(v[N_MAIN, j["C"]] - v[N_MAIN, j["B"]]) / v[N_MAIN, j["B"]]
    _say(f"MAT19-3p (f'') |C - B|/B = {ff:.4f} (predicted [0.3, 0.8]; asserted > {F_CONTROL_FLOOR:.0e})")
    _say("MAT19-3p (S) rho c theta_bar(theta^1)/(dt_1 q_bar) - 1: " + "  ".join(
        f"{k} {RHO_C * s_be.values[0, j[k]] / (DT_1 * run['qbar'][k]) - 1.0:+.4e}" for k in order))
    if run["ghost"] is not None:
        ratio = run["ghost"][0] / (run["g_w"] / n_sub) if run["g_w"] else float("nan")
        _say(f"MAT19-3p (U3c) ghost-included Omega: volume miss {run['ghost'][0]:.4e}, value miss "
             f"{run['ghost'][1]:.4e}; G_w/N_sub {run['g_w'] / n_sub:.4e}; factor {ratio:.3f} "
             f"(predicted x[0.5, 2]; asserted miss iff G_w > 0)")
    _say(f"MAT19-3p (U4) bitwise rank agreement at n=40: {run['u4']}; O_w {run['o_w']}")

    w = None
    if comm.size == NPZ_WIDTH:
        if comm.rank == 0:
            np.savez(NPZ, series=v, volumes=s.volumes, cells=np.array([n_sub]),
                     sha_module=run["hashes"][0], sha_probe=run["hashes"][1])
            print(f"MAT19-3p (W) wrote {NPZ} (sha module {run['hashes'][0][:12]})", flush=True)
        comm.Barrier()
    else:
        if not os.path.exists(NPZ):
            pytest.fail(f"(W) {NPZ} missing: run W2 (-n 8) first")
        z = np.load(NPZ)
        if str(z["sha_module"]) != run["hashes"][0] or str(z["sha_probe"]) != run["hashes"][1]:
            pytest.fail("(W) npz hash differs from this module / probe.py")
        cells_eq = int(z["cells"][0]) == n_sub
        rel = np.abs(v[1:] / z["series"][1:] - 1.0)
        w = float(np.max(rel))
        _say(f"MAT19-3p (W) -n {comm.size} vs -n {NPZ_WIDTH}: cells equal {cells_eq}; series max rel "
             f"{w:.3e} (band {WIDTH_BAND:.0e})")
        gate = run["drive_t"] > DRIVE_CLOCK_GATE_S
        _say(f"MAT19-3p time gate: drive clock {run['drive_t']:.1f} s vs {DRIVE_CLOCK_GATE_S:.0f} s -> "
             f"{'-n 2 NOT run (resource stop)' if gate else '-n 2 (W4) may run'}")
    _say("MAT19-3p phase clocks: " + "  ".join(f"{k} {t:.1f} s" for k, t in run["clock"].items())
         + f"; summed ru_maxrss {run['rss_kib'] / 2**20:.2f} GiB")

    assert run["tag3"] > TAG3_FLOOR
    assert run["sin_phi"] <= PHI_BAND, f"(Phi) {run['sin_phi']:.3e}"
    assert m_ok, "(M) a probe ball leaves the phantom"
    assert np.all(v[0] == 0.0) and run["means"][0] == 0.0
    assert run["ab"] <= REDUCTION_BAND, f"(Src/AB') {run['ab']:.3e}"
    assert run["p_rel"] <= DISCRETE_IDENTITY_RTOL, f"(Src/P) {run['p_rel']:.3e}"
    assert np.max(relA) <= B_BAND, f"(A) {np.max(relA):.3e}"
    assert relA_be <= B_BAND and volA <= B_BAND, (relA_be, volA)
    if comm.size > 1:
        assert run["u4"], "(U4) ranks disagree bitwise"
    if run["ghost"] is not None and run["g_w"] > 0:
        assert run["ghost"][0] > B_BAND, "(U3c) ghost-included volume did not miss"
    assert c1 > WIDTH_BAND, f"(C1'') {c1:.3e}"
    assert ff > F_CONTROL_FLOOR, f"(f'') {ff:.4f}"
    if w is not None:
        assert cells_eq, "(W) cells differ"
        assert w <= WIDTH_BAND, f"(W) {w:.3e} > {WIDTH_BAND:.0e}"
    assert sig <= C4_COVARIANCE_BAND, f"(Sigma) {sig:.4e}"
