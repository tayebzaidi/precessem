"""`IMP-3` step 9 (§9 item 428, §10 I25): a PEC oblate disc in the coil-driven phantom.

Step 1's fixture (item 174) with `GEO-36` step 7's disc in place of the rod:
the loaded F-small 4-leg birdcage at 64 MHz (untuned sigma = 800 S/m solid
coil, 50 Ohm lumped sheets, degree 1), swept with and without an exactly
revolved oblate spheroid (a = 6 mm, c = 1.2 mm, aspect 5, symmetry axis x-hat,
centre (0.015, 0, 0)) cut out of the phantom as a hole whose wall (facet tag
402) is pinned by the default exterior PEC.

Asserted (every band imported, never restated):

    (A0) census within CELL_COUNT_BAND (W1: 116 085 / SPHEROID_HOLE_CELL_RECORD;
         W2: 199 920, disc arm guarded at ROD_CELL_GUARD_W2 before its solve);
    (R0) W1 without-disc Z reproduces the 64 MHz record within
         REUSE_REPRODUCTION_RTOL, entrywise;
    (i)  ||S - S^T|| / ||S|| <= RECIPROCITY_BAND, both arms;
    (ii) sigma_max(S) <= 1 + PASSIVITY_SIGMA_TOLERANCE and every column
         power <= 1, both arms;
    (P)  P_src,exact = P_vol + P_sheet,exact within DISCRETE_IDENTITY_RTOL,
         every drive of both arms;
    (M)  the five sigma_y (P2 <-> P4) Z pairs within ADJACENT_SPREAD_BAND,
         disc arm;
    (E)  n x E = 0 on the tag-402 dofs (count > 0, max |E| <= PEC_WALL_MAX_E);
    (S)  pooled off-diagonal spread >= POOLED_SEPARATION_FLOOR x the worst
         class spread, without-disc arm;
    containment of the 1 g ball;
    (C1, W1) the disc wall left natural: (P) still closes, and the knob is
         wired (max |E| on the 402 dofs > PEC_WALL_MAX_E).

The estimate (printed, never asserted).  Near a rim point the background
field E_bg is taken uniform over the 1 g ball.  In the local frame
(local x = world z toward the rim, local y = world y, local z = world x, the
disc normal) write E_bg = sum_i b_i e_i, with e_i the quasi-static PEC field
for a unit drive along local axis i (linear in the drive).  Then
<|E|^2> = sum_ij Re(b_i conj(b_j)) <e_i . e_j> (the e_i are real).  The ball
centred on a rim point (local (a, 0, 0)) is invariant under y_l -> -y_l and
under z_l -> -z_l (the disc's mid-plane), and so are the body and its gel
mask.  Under the mirror M_y, the field of a drive along x or z is M_y-even
(e(M_y p) = M_y e(p)) while the field of the y drive is M_y-odd
(e(M_y p) = -M_y e(p)); M_y is orthogonal, so e_x . e_y and e_z . e_y are odd
functions of y_l and average to zero over the ball.  The same argument with
M_z kills e_x . e_z.  Every cross term vanishes, leaving the diagonal.
Hence Q_cf = <|E|^2> / |E_bg|^2 = eta_par M_par + eta_t M_t + eta_n M_n, with
eta_i = |b_i|^2 / |b|^2 (sampled on the disc's in-plane diameter in the
without-disc arm) and M_i = <|e_i|^2> over the ball's gel points.  M_par is
`pec_oblate_spheroid_field(drive="inplane")` at p_l; M_t is the same closed
form at the 90 deg-rotated point (y_l, -x_l, z_l) (|E| is rotation
invariant); M_n is `drive="axial"` at p_l.  Step 1's eta_z M_cf + (1 - eta_z)
is the rod's two-way version of the same identity.

Printed beside predictions written before the first window: the rim 1 g Q,
the three eta's, M_par / M_t / M_n, Q / estimate, the ball masses, rim
fidelity, C2 and C3.  The ratio's accuracy is `ANS-10`'s / (7c)'s.

Scope: degree 1, 64 MHz, the untuned solid-coil 4x4, one aspect-5 disc with
normal x-hat at (0.015, 0, 0).  No absolute or pointwise SAR, heating,
compliance, convergence, tuned or copper coil, F-human, other placements or a
square-edged plate.
"""

from __future__ import annotations

import gc
import math
import time
from dataclasses import replace

import basix
import numpy as np
from mpi4py import MPI

from fem_em_solver.ports.sparameters import run_n_port_sparameter_sweep
from fem_em_solver.ports.superposition import superpose_drives
from fem_em_solver.post.evaluation import evaluate_vector_field_parallel
from fem_em_solver.post.sar import (
    averaging_ball_radius,
    build_density_field,
    mass_averaged_sar,
)
from fem_em_solver.utils.analytical import (
    oblate_spheroid_depolarisation,
    pec_oblate_spheroid_field,
)
from fem_em_solver.utils.constants import EPSILON_0, MU_0

from tests.complex_mode import complex_only
from tests.mesh.helpers import global_cell_tag_set
from tests.mesh.test_birdcage_conductor_hole import SOLID_CELL_RECORD
from tests.mesh.test_birdcage_implant_spheroid import (
    RIM_POINTS,
    SPHEROID_CENTRE,
    SPHEROID_HOLE_CELL_RECORD,
    SPHEROID_SEMI_MAJOR,
    SPHEROID_SEMI_MINOR,
    _build_spheroid,
    _inside_true_spheroid,
)
from tests.mesh.test_birdcage_port_sheet_prerequisite import CELL_COUNT_BAND
from tests.mesh.test_birdcage_port_sheets import SHEET_IFACE, _build
from tests.mesh.test_birdcage_port_tags import PHANTOM_HEIGHT, PHANTOM_RADIUS
from tests.validation.test_birdcage_b1_plus_map import C4_COVARIANCE_BAND, _solve_driven
from tests.validation.test_birdcage_b1_quadrature import (
    QUADRATURE_STEP_DEG,
    _port_index,
    quadrature_phase_weights,
)
from tests.validation.test_birdcage_power_identity import (
    DISCRETE_IDENTITY_RTOL,
    STEP2D_FREQUENCY_HZ,
    _exact_shares,
)
from tests.validation.test_birdcage_sar_1g_rung import (
    ONE_GRAM_RUNG_CELL_RECORD,
    PHANTOM_RESOLUTION_1G_RUNG,
)
from tests.validation.test_birdcage_sar_mass_averaged import PHANTOM_RHO_KG_PER_M3
from tests.validation.test_imp3_rod_birdcage import (
    AXIS_SAMPLES,
    CF_CHUNK,
    DRIVES,
    IMP3_FLAG_OFF_Z_64MHZ,
    MIRROR_PAIRS,
    PEC_WALL_MAX_E,
    PORT_IDS,
    QDEG_CF,
    ROD_CELL_GUARD_W2,
    SIGMA_X_PAIRS,
    _in,
    _matrix_quantities,
    _outer_box_tags,
    _pair_spread,
    _peak_rss_gib,
    _say,
    _wall_dofs,
    _wall_max_e,
)
from tests.validation.test_lossy_sphere_fullwave import SALINE_EPSILON_R, SALINE_SIGMA
from tests.validation.test_mass_averaged_sar_standard_masses import (
    ONE_GRAM_KG,
    QUADRATURE_DEGREE,
)
from tests.validation.test_port_birdcage_four_port import (
    PASSIVITY_SIGMA_TOLERANCE,
    POOLED_SEPARATION_FLOOR,
    build_four_port_sweep,
)
from tests.validation.test_port_birdcage_lumped_column import (
    ADJACENT_SPREAD_BAND,
    PHANTOM_CELL_TAG,
)
from tests.validation.test_port_birdcage_ring_matrix import REUSE_REPRODUCTION_RTOL
from tests.validation.test_port_lumped_sheet_sweep import RECIPROCITY_BAND
from tests.validation.test_th14_birdcage_copper import OUTER_BOX_TAG

assert len(SIGMA_X_PAIRS) == 5  # printed through _matrix_quantities' "sigma_x"

# --- predictions, written before the first window (printed, never asserted) --
PRED_DISC_CELLS_W2 = (235_000, 300_000)
PRED_MIRROR_SPREAD_DISC = (5.0e-5, 1.5e-3)
PRED_Q_P1 = (1.3, 2.6)
PRED_Q_P24 = (0.9, 2.0)
PRED_Q_OVER_ESTIMATE = (0.7, 1.6)
PRED_ETA_PAR_P1 = 0.99
PRED_ETA_PAR_P24 = (0.3, 0.6)
PRED_EST_P1 = 1.64
PRED_EST_P24 = (1.2, 1.3)
PRED_M = {"par": 1.650, "t": 0.880, "n": 1.068}
PRED_GEL_FRACTION = 0.923
PRED_DP_PHANTOM = (0.005, 0.10)
ROD_D_COLPOWER = -3.55e-5
ROD_MAX_DS = 1.32e-4
PRED_C3_PAIR_MAX = 0.01
PRED_C3_P1P3 = (0.50, 0.95)
PRED_C1_NATURAL_OVER_OFF = (0.7, 1.3)
PRED_C1_PEC_OVER_NATURAL = (1.1, 3.0)
L_X_RECORD = 0.124758


def _rim_fidelity(msh, cell_tags, comm):
    """Owned gel cells within 2 rho_r of each rim point, and how many midpoints
    lie inside the true spheroid (reduced)."""
    import dolfinx

    tdim = msh.topology.dim
    n_owned = msh.topology.index_map(tdim).size_local
    idx = np.asarray(cell_tags.indices)
    val = np.asarray(cell_tags.values)
    gel = idx[(idx < n_owned) & (val == PHANTOM_CELL_TAG)].astype(np.int32)
    mid = dolfinx.mesh.compute_midpoints(msh, tdim, gel)
    radius = 2.0 * SPHEROID_SEMI_MINOR**2 / SPHEROID_SEMI_MAJOR
    out = {}
    for name, p in RIM_POINTS.items():
        near = np.linalg.norm(mid - np.asarray(p), axis=1) < radius
        s = _inside_true_spheroid(mid[near]) if np.any(near) else np.zeros(0)
        out[name] = (
            comm.allreduce(int(near.sum()), op=MPI.SUM),
            comm.allreduce(int(np.count_nonzero(s < 1.0)), op=MPI.SUM),
        )
    return out


def _disc_ball_means(msh, cell_tags, centre, radius):
    """<|E_cf|^2>/E0^2 over the ball at ``centre`` on owned gel cells (numpy, deg 8).

    Three closed forms in the disc's local frame d = p - SPHEROID_CENTRE,
    p_l = (d_z, d_y, d_x): M_par (in-plane drive along local x, toward the
    rim), M_t (in-plane drive along local y, evaluated at the rotated point
    (y_l, -x_l, z_l)), M_n (axial drive, along the normal).  Returns
    ({"par", "t", "n"}, ball volume, quadrature points inside the true
    spheroid, non-finite closed-form points dropped).
    """
    comm = msh.comm
    tdim = msh.topology.dim
    n_owned = msh.topology.index_map(tdim).size_local
    dmap = msh.geometry.dofmap
    if dmap.shape[1] != 4:
        raise RuntimeError(f"expected affine tets, geometry dofmap {dmap.shape}")
    xg = msh.geometry.x
    c = np.asarray(centre, dtype=float)
    sc = np.asarray(SPHEROID_CENTRE, dtype=float)
    pts_ref, wts = basix.make_quadrature(basix.CellType.tetrahedron, QDEG_CF)
    idx = np.asarray(cell_tags.indices)
    val = np.asarray(cell_tags.values)
    cells = np.sort(idx[(idx < n_owned) & (val == PHANTOM_CELL_TAG)].astype(np.int32))
    if cells.size:
        vmin = np.min(np.linalg.norm(xg[dmap[cells]] - c, axis=2), axis=1)
        cells = cells[vmin < radius + 0.01]
    num = {"par": 0.0, "t": 0.0, "n": 0.0}
    vol = 0.0
    n_inside = 0
    n_bad = 0
    kw = dict(semi_major=SPHEROID_SEMI_MAJOR, semi_minor=SPHEROID_SEMI_MINOR, e0=1.0,
              interior="continue")
    for start in range(0, cells.size, CF_CHUNK):
        chunk = cells[start:start + CF_CHUNK]
        v = xg[dmap[chunk]]
        jac = np.stack([v[:, 1] - v[:, 0], v[:, 2] - v[:, 0], v[:, 3] - v[:, 0]], axis=2)
        det = np.abs(np.linalg.det(jac))
        phys = v[:, 0][:, None, :] + np.einsum("cij,qj->cqi", jac, pts_ref)
        chi = (np.sum((phys - c) ** 2, axis=2) < radius**2)
        if not np.any(chi):
            continue
        flat = phys.reshape(-1, 3)
        d = flat - sc
        p_l = np.column_stack([d[:, 2], d[:, 1], d[:, 0]])
        p_rot = np.column_stack([p_l[:, 1], -p_l[:, 0], p_l[:, 2]])
        chi_f = chi.reshape(-1)
        n_inside += int(np.count_nonzero(chi_f & (_inside_true_spheroid(flat) < 1.0)))
        with np.errstate(all="ignore"):
            fields = {
                "par": pec_oblate_spheroid_field(p_l, drive="inplane", **kw),
                "t": pec_oblate_spheroid_field(p_rot, drive="inplane", **kw),
                "n": pec_oblate_spheroid_field(p_l, drive="axial", **kw),
            }
        e2 = {k: np.sum(np.abs(f) ** 2, axis=1) for k, f in fields.items()}
        good = np.ones(flat.shape[0], dtype=bool)
        for k in e2:
            good &= np.isfinite(e2[k])
        n_bad += int(np.count_nonzero(chi_f & ~good))
        w = (det[:, None] * wts[None, :]).reshape(-1) * (chi_f & good)
        for k in e2:
            num[k] += float(np.sum(w * np.where(good, e2[k], 0.0)))
        vol += float(np.sum(w))
    vol = comm.allreduce(vol, op=MPI.SUM)
    means = {k: comm.allreduce(v, op=MPI.SUM) / vol for k, v in num.items()}
    return (means, vol, comm.allreduce(n_inside, op=MPI.SUM),
            comm.allreduce(n_bad, op=MPI.SUM))


def _arm_disc(built, label, disc):
    """Sweep one arm; return scalars and, for the disc arm, the mesh handles."""
    comm = MPI.COMM_WORLD
    t_arm = time.perf_counter()
    case = build_four_port_sweep(
        frequency_hz=STEP2D_FREQUENCY_HZ, built=built, build_only=True
    )
    msh = case["mesh"]
    ncells = int(case["cells"])
    _say(
        f"\n[IMP-3 {label}] {ncells} cells, mesh wall {built[3]['mesh_wall_time_s']:.2f} s, "
        f"build {built[4]:.2f} s; sheets "
        + ", ".join(
            f"P{s['tag'] - SHEET_IFACE} az {s['azimuth_deg']:.3f} deg {s['facets']} facets"
            for s in case["sheets"]
        )
    )
    comm.Barrier()
    t0 = time.perf_counter()
    result = run_n_port_sparameter_sweep(
        case["problem"],
        case["port_defs"],
        lumped_sheet_ports=case["specs"],
        lumped_sheet_facet_tags=case["facet_tags"],
        keep_fields=True,
        degree=case["degree"],
    )
    comm.Barrier()
    t_sweep = time.perf_counter() - t0
    assert tuple(result.port_ids) == PORT_IDS, result.port_ids
    mq = _matrix_quantities(result)

    omega = 2.0 * np.pi * float(STEP2D_FREQUENCY_HZ)
    shares = {
        pid: _exact_shares(case, {"driven": pid, "omega": omega, "fields": result.fields[pid]})
        for pid in PORT_IDS
    }

    azimuths = {f"P{s['tag'] - SHEET_IFACE}": float(s["azimuth_deg"]) for s in case["sheets"]}
    ks = np.array(
        [_port_index(azimuths[p], azimuths["P1"], QUADRATURE_STEP_DEG) for p in PORT_IDS],
        dtype=float,
    )
    weights = {sense: quadrature_phase_weights(ks, sense) for sense in ("ccw", "cw")}
    superposed = {sense: superpose_drives(result, weights[sense]) for sense in weights}
    e_of = {pid: result.fields[pid].e_complex for pid in PORT_IDS}
    e_of.update({sense: superposed[sense].e_complex for sense in superposed})
    sigma_field = result.fields["P1"].sigma_field

    cell_tags = case["cell_tags"]
    all_tags = global_cell_tag_set(msh, cell_tags)
    rho_field = build_density_field(
        msh,
        PHANTOM_RHO_KG_PER_M3,
        cell_tags=cell_tags,
        density_map={int(t): 0.0 for t in sorted(all_tags) if int(t) != PHANTOM_CELL_TAG},
    )
    a_1g = averaging_ball_radius(mass_kg=ONE_GRAM_KG, rho=PHANTOM_RHO_KG_PER_M3)
    t1 = time.perf_counter()
    sar = {}
    mass = {}
    for rim, centre in RIM_POINTS.items():
        for d in DRIVES:
            out = mass_averaged_sar(
                e_of[d], sigma=sigma_field, rho=rho_field, center=centre,
                radius=a_1g, quadrature_degree=QUADRATURE_DEGREE,
            )
            sar[(rim, d)] = float(out["averaged_sar_w_per_kg"])
            mass[rim] = float(out["mass_kg"])
    t_sar = time.perf_counter() - t1

    out = {
        "label": label,
        "cells": ncells,
        "mq": mq,
        "shares": shares,
        "sar": sar,
        "mass": mass,
        "a_1g": a_1g,
        "t_sweep": t_sweep,
        "t_sar": t_sar,
        "rim_fidelity": _rim_fidelity(msh, built[1], comm),
    }

    if not disc:
        # Background on the disc's in-plane diameter, without-disc arm only.
        zs = np.linspace(-SPHEROID_SEMI_MAJOR, SPHEROID_SEMI_MAJOR, AXIS_SAMPLES)
        pts = np.column_stack(
            [np.full_like(zs, SPHEROID_CENTRE[0]), np.full_like(zs, SPHEROID_CENTRE[1]), zs]
        )
        diam = {}
        for pid in PORT_IDS:
            f = result.fields[pid]
            re, ok_r = evaluate_vector_field_parallel(f.e_real, pts)
            im, ok_i = evaluate_vector_field_parallel(f.e_imag, pts)
            assert bool(np.all(ok_r)) and bool(np.all(ok_i)), "diameter sample off-mesh"
            diam[pid] = np.real(re) + 1j * np.real(im)
        for sense in ("ccw", "cw"):
            diam[sense] = sum(w * diam[p] for w, p in zip(weights[sense], PORT_IDS))
        out["eta"] = {}
        for d in DRIVES:
            tot = float(np.sum(np.abs(diam[d]) ** 2))
            out["eta"][d] = {
                "par": float(np.sum(np.abs(diam[d][:, 2]) ** 2)) / tot,
                "t": float(np.sum(np.abs(diam[d][:, 1]) ** 2)) / tot,
                "n": float(np.sum(np.abs(diam[d][:, 0]) ** 2)) / tot,
            }
    else:
        v_space = e_of["P1"].function_space
        owned, n_wall = _wall_dofs(msh, built[2], v_space)
        out["wall_dofs"] = n_wall
        out["wall_max"] = {pid: _wall_max_e(e_of[pid], owned, comm) for pid in PORT_IDS}
        out["case"] = case
        out["built"] = built
        out["owned_wall"] = owned
        out["rho_field"] = rho_field
    out["t_arm"] = time.perf_counter() - t_arm
    del superposed, e_of, result
    gc.collect()
    return out


def _print_arm(a):
    mq = a["mq"]
    _say(f"[IMP-3 {a['label']}] Z (Ohm), column k = port k driven:")
    for i in range(4):
        _say("    Z_%dk = " % (i + 1) + "  ".join(f"{v:+.9e}" for v in mq["z"][i]))
    _say(f"[IMP-3 {a['label']}] S (z0 = 50 Ohm):")
    for i in range(4):
        _say("    S_%dk = " % (i + 1) + "  ".join(f"{v:+.9e}" for v in mq["s"][i]))
    _say(
        f"    (i) ||S-S^T||/||S|| = {mq['recip']:.3e} (band {RECIPROCITY_BAND:g});  "
        f"(ii) sigma(S) = {', '.join(f'{x:.9f}' for x in mq['sigma'])} "
        f"(<= 1 + {PASSIVITY_SIGMA_TOLERANCE:g}); column powers "
        f"{', '.join(f'{x:.9f}' for x in mq['colpow'])}"
    )
    _say(
        "    class spreads: "
        + "  ".join(f"{k} {v:.4%}" for k, v in mq["spreads"].items())
        + f"; pooled {mq['pooled']:.4%}; separation {mq['separation']:.1f}x "
        f"(floor {POOLED_SEPARATION_FLOOR:g})"
    )
    _say(
        "    (M) sigma_y pairs (Z22,Z44) (Z12,Z14) (Z21,Z41) (Z32,Z34) (Z23,Z43): "
        + " / ".join(f"{v:.4%}" for v in mq["mirror"])
        + f" (band {ADJACENT_SPREAD_BAND:.2%})"
    )
    _say(
        "    (C2) sigma_x pairs (Z11,Z33) (Z21,Z23) (Z12,Z32) (Z41,Z43) (Z14,Z34): "
        + " / ".join(f"{v:.4%}" for v in mq["sigma_x"])
    )
    for pid in PORT_IDS:
        ex = a["shares"][pid]
        resid = abs(ex["p_src"] - ex["p_vol"] - ex["sheet_field_total"]) / abs(ex["p_src"])
        _say(
            f"    (P) {pid}: P_src {ex['p_src']:.9e} W  P_vol {ex['p_vol']:.9e} W  "
            f"P_sheet {ex['sheet_field_total']:.9e} W  rel {resid:.3e} "
            f"(band {DISCRETE_IDENTITY_RTOL:g}); phantom {ex['phantom']:.9e} W "
            f"({ex['phantom'] / ex['p_src']:.3e} of P_src)"
        )
    _say(
        f"    1 g ball a = {a['a_1g'] * 1e3:.4f} mm; meshed masses "
        + ", ".join(f"{t} {m * 1e3:.6f} g" for t, m in a["mass"].items())
    )
    for rim in RIM_POINTS:
        _say(
            f"    SAR_1g at {rim} rim point (W/kg): "
            + "  ".join(f"{d} {a['sar'][(rim, d)]:.6e}" for d in DRIVES)
        )
    _say(
        "    rim fidelity (gel cells within 2 rho_r = 0.48 mm, of which midpoint inside "
        "the true spheroid): "
        + ", ".join(f"{r} {n} / {i}" for r, (n, i) in a["rim_fidelity"].items())
    )
    _say(
        f"    times: arm {a['t_arm']:.1f} s, sweep {a['t_sweep']:.2f} s, "
        f"1 g SARs {a['t_sar']:.2f} s"
    )


def _c3(a):
    s = a["sar"]
    worst = 0.0
    _say(
        f"    (C3) field-level mirror, {a['label']} (C4_COVARIANCE_BAND "
        f"{C4_COVARIANCE_BAND:.0%} shown as a reference only):"
    )
    for rim in RIM_POINTS:
        p24 = _pair_spread(s[(rim, "P2")], s[(rim, "P4")])
        ccw_cw = _pair_spread(s[(rim, "ccw")], s[(rim, "cw")])
        x, y = s[(rim, "P1")], s[(rim, "P3")]
        p13 = abs(x - y) / max(abs(x), abs(y))
        worst = max(worst, p24, ccw_cw)
        _say(
            f"        {rim}: P2/P4 {p24:.4%}, ccw/cw {ccw_cw:.4%} (each predicted <= "
            f"{PRED_C3_PAIR_MAX:.0%}); P1/P3 |P1-P3|/max {p13:.4%} (predicted "
            f"[{PRED_C3_P1P3[0]:.0%}, {PRED_C3_P1P3[1]:.0%}], {_in(p13, PRED_C3_P1P3)})"
        )
    for d in DRIVES:
        tb = _pair_spread(s[("top", d)], s[("bottom", d)])
        worst = max(worst, tb)
        _say(f"        top/bottom {d}: {tb:.4%}")
    _say(
        f"        worst kept-mirror pair {worst:.4%} -> "
        f"{'<=' if worst <= PRED_C3_PAIR_MAX else '>'} prediction"
    )


def _case_disc(phantom_resolution):
    """Without-disc arm first, freed, then the disc arm (GEO-23 order)."""
    comm = MPI.COMM_WORLD
    t_case = time.perf_counter()
    window = "W1" if phantom_resolution is None else "W2"
    kw = {} if phantom_resolution is None else {"phantom_resolution": phantom_resolution}

    built_off = _build(True, **kw)
    off = _arm_disc(built_off, f"{window} without disc", disc=False)
    del built_off
    gc.collect()

    built_disc = _build_spheroid(**kw)
    disc_cells = int(built_disc[0].topology.index_map(3).size_global)
    _say(f"\n[IMP-3 {window}] disc-arm build: {disc_cells} cells (guard {ROD_CELL_GUARD_W2})")
    assert disc_cells < ROD_CELL_GUARD_W2, (disc_cells, ROD_CELL_GUARD_W2)
    disc = _arm_disc(built_disc, f"{window} with disc", disc=True)

    t0 = time.perf_counter()
    m = {}
    for rim, centre in RIM_POINTS.items():
        m[rim] = _disc_ball_means(disc["case"]["mesh"], built_disc[1], centre, off["a_1g"])
    t_cf = time.perf_counter() - t0
    e = math.sqrt(1.0 - (SPHEROID_SEMI_MINOR / SPHEROID_SEMI_MAJOR) ** 2)
    lx, lz = oblate_spheroid_depolarisation(e)
    case = {"window": window, "off": off, "disc": disc, "disc_cells": disc_cells,
            "m": m, "lx": lx, "lz": lz, "t_cf": t_cf}
    case["t_case"] = time.perf_counter() - t_case
    return case


def _report(case):
    off, disc = case["off"], case["disc"]
    w = case["window"]
    _print_arm(off)
    _print_arm(disc)
    _say(
        f"    (E) tag-402 dofs (owned, reduced) {disc['wall_dofs']}; max |E| per drive "
        + ", ".join(f"{p} {v:.3e}" for p, v in disc["wall_max"].items())
        + f" (bound {PEC_WALL_MAX_E:g})"
    )
    mx = max(disc["mq"]["mirror"])
    _say(
        f"    (S) disc-arm separation {disc['mq']['separation']:.1f}x (printed); disc-arm "
        f"mirror pairs max {mx:.4%} (predicted [{PRED_MIRROR_SPREAD_DISC[0]:.3%}, "
        f"{PRED_MIRROR_SPREAD_DISC[1]:.2%}], {_in(mx, PRED_MIRROR_SPREAD_DISC)})"
    )
    ds = disc["mq"]["s"] - off["mq"]["s"]
    dcol = disc["mq"]["colpow"][0] - off["mq"]["colpow"][0]
    _say(f"[IMP-3 {w}] (C2) what the 4x4 sees of the disc (predicted, printed):")
    for pid in PORT_IDS:
        p0 = off["shares"][pid]["phantom"]
        p1 = disc["shares"][pid]["phantom"]
        tag = (f"predicted [{PRED_DP_PHANTOM[0]:.1%}, {PRED_DP_PHANTOM[1]:.0%}], "
               f"{_in((p1 - p0) / p0, PRED_DP_PHANTOM)}" if pid == "P1" else "printed")
        _say(
            f"        {pid}: P_phantom {p0:.9e} -> {p1:.9e} W, dP/P {(p1 - p0) / p0:+.4%} "
            f"({tag})"
        )
    _say(
        f"        d sum_i |S_i1|^2 = {dcol:+.3e} (rod {ROD_D_COLPOWER:+.2e}); max |dS| = "
        f"{float(np.max(np.abs(ds))):.3e} (rod {ROD_MAX_DS:.2e})"
    )
    _c3(off)
    _c3(disc)
    f_hz = float(STEP2D_FREQUENCY_HZ)
    k_mag = 2.0 * np.pi * f_hz * math.sqrt(
        MU_0 * EPSILON_0 * abs(SALINE_EPSILON_R - 1j * SALINE_SIGMA
                               / (2.0 * np.pi * f_hz * EPSILON_0))
    )
    _say(
        f"[IMP-3 {w}] rim 1 g ratio Q = SAR_1g,rim(disc)/SAR_1g,rim(without); "
        f"L_x = {case['lx']:.6f} (record {L_X_RECORD}), L_z = {case['lz']:.6f}, "
        f"(1/L_x)^2 = {1.0 / case['lx'] ** 2:.2f}; |k| = {k_mag:.3f} 1/m, "
        f"(|k| a)^2 = {(k_mag * SPHEROID_SEMI_MAJOR) ** 2:.4%}; closed-form means "
        f"{case['t_cf']:.1f} s"
    )
    _say(
        f"    ball masses: without {off['mass']['top'] * 1e3:.6f} / "
        f"{off['mass']['bottom'] * 1e3:.6f} g, with disc {disc['mass']['top'] * 1e3:.6f} / "
        f"{disc['mass']['bottom'] * 1e3:.6f} g (predicted 1.000 / {PRED_GEL_FRACTION} g)"
    )
    for rim in RIM_POINTS:
        means, vol, n_in, n_bad = case["m"][rim]
        vball = 4.0 / 3.0 * math.pi * off["a_1g"] ** 3
        _say(
            f"    {rim}: M_par {means['par']:.6f} (pred {PRED_M['par']}), M_t "
            f"{means['t']:.6f} (pred {PRED_M['t']}), M_n {means['n']:.6f} (pred "
            f"{PRED_M['n']}); gel ball volume / 4/3 pi a^3 = {vol / vball:.6f} (pred "
            f"{PRED_GEL_FRACTION}); closed-form quadrature points inside the body {n_in}, "
            f"non-finite dropped {n_bad}"
        )
    for d in DRIVES:
        eta = off["eta"][d]
        if d == "P1":
            eta_tag = f"pred ~{PRED_ETA_PAR_P1}"
            q_band = PRED_Q_P1
        elif d in ("P2", "P4"):
            eta_tag = f"pred [{PRED_ETA_PAR_P24[0]}, {PRED_ETA_PAR_P24[1]}], " \
                      f"{_in(eta['par'], PRED_ETA_PAR_P24)}"
            q_band = PRED_Q_P24
        else:
            eta_tag = "no prediction"
            q_band = None
        for rim in RIM_POINTS:
            means = case["m"][rim][0]
            q = disc["sar"][(rim, d)] / off["sar"][(rim, d)]
            est = eta["par"] * means["par"] + eta["t"] * means["t"] + eta["n"] * means["n"]
            qtag = (f"pred [{q_band[0]}, {q_band[1]}], {_in(q, q_band)}" if q_band
                    else "no prediction")
            _say(
                f"    {d:>3s} {rim:>6s}: Q = {q:.6f} ({qtag}); eta_par {eta['par']:.6f} "
                f"({eta_tag}), eta_t {eta['t']:.6f}, eta_n {eta['n']:.6f}; estimate "
                f"{est:.6f}; Q/estimate {q / est:.6f} ({_in(q / est, PRED_Q_OVER_ESTIMATE)} "
                f"[{PRED_Q_OVER_ESTIMATE[0]}, {PRED_Q_OVER_ESTIMATE[1]}])"
            )
    _say(
        f"    estimate predictions: P1 ~{PRED_EST_P1}, P2/P4 [{PRED_EST_P24[0]}, "
        f"{PRED_EST_P24[1]}]"
    )
    _say(
        f"[IMP-3 {w}] cells without {off['cells']} / with {disc['cells']}; case "
        f"{case['t_case']:.1f} s; summed peak RSS {_peak_rss_gib():.2f} GiB at -n "
        f"{MPI.COMM_WORLD.size}"
    )


def _assert_identities_disc(case):
    off, disc = case["off"], case["disc"]
    for arm in (off, disc):
        mq = arm["mq"]
        assert mq["recip"] <= RECIPROCITY_BAND, (arm["label"], mq["recip"])
        assert float(mq["sigma"][0]) <= 1.0 + PASSIVITY_SIGMA_TOLERANCE, (arm["label"], mq["sigma"])
        assert np.all(mq["colpow"] <= 1.0), (arm["label"], mq["colpow"])
        for pid in PORT_IDS:
            ex = arm["shares"][pid]
            resid = abs(ex["p_src"] - ex["p_vol"] - ex["sheet_field_total"]) / abs(ex["p_src"])
            assert resid <= DISCRETE_IDENTITY_RTOL, (arm["label"], pid, resid)
    for v in disc["mq"]["mirror"]:
        assert v <= ADJACENT_SPREAD_BAND, disc["mq"]["mirror"]
    assert disc["wall_dofs"] > 0
    for pid, v in disc["wall_max"].items():
        assert v <= PEC_WALL_MAX_E, (pid, v)
    assert off["mq"]["separation"] >= POOLED_SEPARATION_FLOOR, off["mq"]["separation"]
    a = off["a_1g"]
    assert SPHEROID_CENTRE[0] + a < PHANTOM_RADIUS, (SPHEROID_CENTRE[0] + a, PHANTOM_RADIUS)
    assert SPHEROID_SEMI_MAJOR + a < 0.5 * PHANTOM_HEIGHT, (SPHEROID_SEMI_MAJOR + a,
                                                            0.5 * PHANTOM_HEIGHT)


def _control_c1_disc(case):
    """(C1) the disc wall left natural, P1 only, on a fresh solver."""
    comm = MPI.COMM_WORLD
    disc = case["disc"]
    sweep = disc["case"]
    msh = sweep["mesh"]
    outer, census = _outer_box_tags(msh, disc["built"][2], comm)
    _say(
        f"\n[IMP-3 W1] (C1) census: exterior {census['exterior']} = outer {census['outer']} "
        f"+ 402 {census['wall_exterior']} (402 total {census['wall']}); tag 499 "
        f"pre-existing on {census['tag_collision']} facets"
    )
    assert census["tag_collision"] == 0, census
    assert census["exterior"] == census["outer"] + census["wall_exterior"], census
    assert census["wall_exterior"] == census["wall"] > 0, census
    t0 = time.perf_counter()
    problem = replace(sweep["problem"], facet_tags=outer, pec_facet_tags=(OUTER_BOX_TAG,))
    nat_case = dict(sweep, problem=problem)
    solved = _solve_driven(nat_case, "P1")
    ex = _exact_shares(nat_case, solved)
    resid = abs(ex["p_src"] - ex["p_vol"] - ex["sheet_field_total"]) / abs(ex["p_src"])
    wall_max = _wall_max_e(solved["fields"].e_complex, disc["owned_wall"], comm)
    a_1g = case["off"]["a_1g"]
    sar_nat = {}
    for rim, centre in RIM_POINTS.items():
        out = mass_averaged_sar(
            solved["fields"].e_complex, sigma=solved["fields"].sigma_field,
            rho=disc["rho_field"], center=centre, radius=a_1g,
            quadrature_degree=QUADRATURE_DEGREE,
        )
        sar_nat[rim] = float(out["averaged_sar_w_per_kg"])
    _say(
        f"    (C1) natural wall, P1: (P) rel {resid:.3e} (band {DISCRETE_IDENTITY_RTOL:g}); "
        f"max |E| on 402 dofs {wall_max:.3e} (> {PEC_WALL_MAX_E:g} = knob wired); "
        f"phantom {ex['phantom']:.9e} W; {time.perf_counter() - t0:.1f} s"
    )
    for rim in RIM_POINTS:
        r_off = sar_nat[rim] / case["off"]["sar"][(rim, "P1")]
        r_pec = disc["sar"][(rim, "P1")] / sar_nat[rim]
        _say(
            f"        {rim}: SAR_1g natural {sar_nat[rim]:.6e} W/kg; natural/without "
            f"{r_off:.6f} (predicted [{PRED_C1_NATURAL_OVER_OFF[0]:g}, "
            f"{PRED_C1_NATURAL_OVER_OFF[1]:g}], {_in(r_off, PRED_C1_NATURAL_OVER_OFF)}); "
            f"PEC/natural {r_pec:.6f} (predicted [{PRED_C1_PEC_OVER_NATURAL[0]:g}, "
            f"{PRED_C1_PEC_OVER_NATURAL[1]:g}], {_in(r_pec, PRED_C1_PEC_OVER_NATURAL)})"
        )
    return {"resid": resid, "wall_max": wall_max}


@complex_only
def test_imp3_disc_default_phantom():
    """W1: the record meshes, (A0), (R0), identities, C1."""
    case = _case_disc(None)
    off, disc = case["off"], case["disc"]
    rel = np.abs(off["mq"]["z"] - IMP3_FLAG_OFF_Z_64MHZ) / np.abs(IMP3_FLAG_OFF_Z_64MHZ)
    _report(case)
    _say(
        f"[IMP-3 W1] (A0) without {off['cells']} vs {SOLID_CELL_RECORD}, with "
        f"{disc['cells']} vs {SPHEROID_HOLE_CELL_RECORD} (band {CELL_COUNT_BAND:.0%}); "
        f"(R0) max entrywise |Z - Z_rec|/|Z_rec| = {float(rel.max()):.3e} "
        f"(band {REUSE_REPRODUCTION_RTOL:g})"
    )
    c1 = _control_c1_disc(case)
    assert abs(off["cells"] / SOLID_CELL_RECORD - 1.0) <= CELL_COUNT_BAND
    assert abs(disc["cells"] / SPHEROID_HOLE_CELL_RECORD - 1.0) <= CELL_COUNT_BAND
    assert float(rel.max()) <= REUSE_REPRODUCTION_RTOL, rel
    _assert_identities_disc(case)
    assert c1["resid"] <= DISCRETE_IDENTITY_RTOL, c1
    assert c1["wall_max"] > PEC_WALL_MAX_E, c1


@complex_only
def test_imp3_disc_one_gram_rung():
    """W2: the 1 g rung, identities, the rim ratio of record."""
    case = _case_disc(PHANTOM_RESOLUTION_1G_RUNG)
    off, disc = case["off"], case["disc"]
    _report(case)
    _say(
        f"[IMP-3 W2] (A0) without {off['cells']} vs {ONE_GRAM_RUNG_CELL_RECORD} (band "
        f"{CELL_COUNT_BAND:.0%}); with {disc['cells']} — first record, predicted "
        f"[{PRED_DISC_CELLS_W2[0]}, {PRED_DISC_CELLS_W2[1]}] "
        f"({_in(disc['cells'], PRED_DISC_CELLS_W2)}), guard {ROD_CELL_GUARD_W2}"
    )
    assert abs(off["cells"] / ONE_GRAM_RUNG_CELL_RECORD - 1.0) <= CELL_COUNT_BAND
    assert disc["cells"] < ROD_CELL_GUARD_W2
    _assert_identities_disc(case)
