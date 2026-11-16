"""`IMP-3` step 1 (§9 item 174, chain I6): a capped PEC rod in the coil-driven phantom.

The loaded F-small 4-leg birdcage at 64 MHz (PORT-11 / PORT-16's untuned
σ = 800 S/m solid coil, 50 Ohm lumped sheets), swept with and without
`GEO-36`'s rod (r = 1 mm, 30 mm between cap centres, axis z, centre
(0.015, 0, 0) on leg 1's azimuth) cut out of the phantom as a hole whose wall
(facet tag 402) is pinned by the default exterior PEC.

Asserted (every band imported, never restated):

    (A0) census within CELL_COUNT_BAND (W1: 116 085 / 167 371; W2: 199 920,
         rod arm guarded at 320 000 before its solve);
    (R0) W1 without-rod Z reproduces the 64 MHz record within
         REUSE_REPRODUCTION_RTOL, entrywise;
    (i)  ||S - S^T|| / ||S|| <= RECIPROCITY_BAND, both arms;
    (ii) sigma_max(S) <= 1 + PASSIVITY_SIGMA_TOLERANCE and every column
         power <= 1, both arms;
    (P)  P_src,exact = P_vol + P_sheet,exact within DISCRETE_IDENTITY_RTOL,
         every drive of both arms;
    (M)  the five sigma_y (P2 <-> P4) Z pairs within ADJACENT_SPREAD_BAND,
         rod arm;
    (E)  n x E = 0 on the tag-402 dofs (count > 0, max |E| <= 1e-12), rod arm;
    (S)  pooled off-diagonal spread >= POOLED_SEPARATION_FLOOR x the worst
         class spread, without-rod arm;
    (C1, W1) the rod wall left natural: (P) still closes, and the knob is
         wired (max |E| on the 402 dofs > 1e-12).

Printed beside predictions written before the first window: the tip 1 g SAR
ratio Q, the ball masses, the closed-form M_cf, L_z, eta_z, the estimate,
the quasi-static scale, C2 (what the 4x4 sees) and C3 (field-level mirror).
None of those is asserted; the ratio's accuracy is `ANS-10`'s.

Scope: degree 1, 64 MHz, the untuned solid-coil 4x4, one axial rod on leg 1's
azimuth.  No absolute SAR, heating, compliance, convergence, tuned or copper
coil, F-human, other placements, ISO 10974 or ASTM F2182.
"""

from __future__ import annotations

import gc
import math
import resource
import time
from dataclasses import replace

import basix
import dolfinx
import numpy as np
import pytest
from dolfinx import fem
from mpi4py import MPI

from fem_em_solver.io.mesh import BIRDCAGE_IMPLANT_SURFACE_TAG
from fem_em_solver.ports.sparameters import run_n_port_sparameter_sweep
from fem_em_solver.ports.superposition import superpose_drives
from fem_em_solver.post.evaluation import evaluate_vector_field_parallel
from fem_em_solver.post.sar import (
    averaging_ball_radius,
    build_density_field,
    mass_averaged_sar,
)
from fem_em_solver.utils.analytical import (
    pec_prolate_spheroid_field,
    prolate_spheroid_depolarisation,
)
from fem_em_solver.utils.constants import EPSILON_0, MU_0

from tests.complex_mode import complex_only
from tests.mesh.helpers import global_cell_tag_set
from tests.mesh.test_birdcage_conductor_hole import SOLID_CELL_RECORD
from tests.mesh.test_birdcage_implant_rod import ROD_CENTRE, _build_rod
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
from tests.validation.test_lossy_sphere_fullwave import SALINE_EPSILON_R, SALINE_SIGMA
from tests.validation.test_mass_averaged_sar_standard_masses import (
    ONE_GRAM_KG,
    QUADRATURE_DEGREE,
)
from tests.validation.test_port_birdcage_four_port import (
    PASSIVITY_SIGMA_TOLERANCE,
    POOLED_SEPARATION_FLOOR,
    _circulant_classes,
    _class_spread,
    build_four_port_sweep,
)
from tests.validation.test_port_birdcage_lumped_column import (
    ADJACENT_SPREAD_BAND,
    PHANTOM_CELL_TAG,
)
from tests.validation.test_port_birdcage_ring_matrix import REUSE_REPRODUCTION_RTOL
from tests.validation.test_port_lumped_sheet_sweep import RECIPROCITY_BAND
from tests.validation.test_th14_birdcage_copper import OUTER_BOX_TAG

# --- records (version-tagged; the (1*) licence) -----------------------------
# W1 rod-arm census: `20261002T110559Z_GEO-36-step2-W2.log:3612`.
ROD_CELL_RECORD_W1 = 167_371
# W2 rod-arm guard and first-record prediction (§9 item 174 (A0)).
ROD_CELL_GUARD_W2 = 320_000
ROD_CELL_PREDICTED_W2 = (225_000, 275_000)

# (R0) the 64 MHz record, `20260912T123236Z_PORT-14-step2d-64mhz.log.gz:1866-1869`
# (digit-identical at `20261001T171613Z_PORT-21-step5-W2.log:3670-3673`).
IMP3_FLAG_OFF_Z_64MHZ = np.array(
    [
        [+2.647082952e01 + 4.646185233e01j, +1.877248178e01 + 6.880903941e-01j,
         +1.429887332e01 - 4.749644673e00j, +1.877380478e01 + 6.929877997e-01j],
        [+1.877079735e01 + 6.864775531e-01j, +2.647408805e01 + 4.647247309e01j,
         +1.877301423e01 + 6.893105431e-01j, +1.429704807e01 - 4.748654779e00j],
        [+1.429428638e01 - 4.749063864e00j, +1.876908728e01 + 6.896008822e-01j,
         +2.649109802e01 + 4.648038833e01j, +1.876465515e01 + 6.882647066e-01j],
        [+1.877383419e01 + 6.947656906e-01j, +1.429928713e01 - 4.746587865e00j,
         +1.877029497e01 + 6.913645278e-01j, +2.647483350e01 + 4.645438954e01j],
    ],
    dtype=np.complex128,
)

# (E) `IMP-2` (iii)'s bound (`test_imp2_spheroid_fem.py:446`).
PEC_WALL_MAX_E = 1.0e-12

# Equivalent spheroid of the capped rod: tip to tip 32 mm, r = 1 mm (aspect 16).
SPHEROID_SEMI_MAJOR = 0.016
SPHEROID_SEMI_MINOR = 0.001
TIPS = {"top": (0.015, 0.0, +0.016), "bottom": (0.015, 0.0, -0.016)}
QDEG_CF = 8
CF_CHUNK = 5000
AXIS_SAMPLES = 9

PORT_IDS = ("P1", "P2", "P3", "P4")
DRIVES = PORT_IDS + ("ccw", "cw")
# sigma_y swaps P2 <-> P4 (0-based 1 <-> 3).
MIRROR_PAIRS = (((1, 1), (3, 3)), ((0, 1), (0, 3)), ((1, 0), (3, 0)),
                ((2, 1), (2, 3)), ((1, 2), (3, 2)))
# sigma_x (broken) swaps P1 <-> P3 (0-based 0 <-> 2) — C2, printed.
SIGMA_X_PAIRS = (((0, 0), (2, 2)), ((1, 0), (1, 2)), ((0, 1), (2, 1)),
                 ((3, 0), (3, 2)), ((0, 3), (2, 3)))

# --- predictions, written before the first window (printed, never asserted) --
PRED_Q = (2.0, 10.0)
PRED_MASS_ROD_G = 0.9817
PRED_M_CF = (2.5, 6.0)
PRED_ETA_Z = (0.5, 1.0)
PRED_Q_OVER_ESTIMATE = (0.6, 2.5)
PRED_MIRROR_SPREAD_ROD = (1.0e-4, 1.5e-3)
PRED_SEPARATION_ROD = (300.0, 1000.0)
PRED_DP_PHANTOM = (0.01, 0.15)
PRED_D_COLPOWER = (5.0e-6, 1.0e-4)
PRED_MAX_DS = (1.0e-6, 1.0e-3)
PRED_C3_PAIR_MAX = 0.01
PRED_C3_P1P3 = (0.50, 0.95)
PRED_C1_NATURAL_OVER_OFF = (0.7, 1.3)
PRED_C1_PEC_OVER_NATURAL = (2.0, 10.0)


def _rank0():
    return MPI.COMM_WORLD.rank == 0


def _say(msg):
    if _rank0():
        print(msg, flush=True)


def _in(x, band):
    return "IN" if band[0] <= x <= band[1] else "OUT"


def _owned_count(a, n_owned, comm):
    return comm.allreduce(int(np.count_nonzero(np.asarray(a) < n_owned)), op=MPI.SUM)


def _wall_dofs(msh, build_facet_tags, v_space):
    """Local copy of `IMP-1`'s ``_cavity_dofs`` on tag 402 (owned count reduced)."""
    fdim = msh.topology.dim - 1
    msh.topology.create_connectivity(fdim, msh.topology.dim)
    facets = np.asarray(build_facet_tags.find(BIRDCAGE_IMPLANT_SURFACE_TAG), dtype=np.int32)
    dofs = fem.locate_dofs_topological(v_space, fdim, facets)
    n_owned = v_space.dofmap.index_map.size_local * v_space.dofmap.index_map_bs
    owned = dofs[dofs < n_owned]
    return owned, msh.comm.allreduce(int(owned.size), op=MPI.SUM)


def _wall_max_e(e_complex, owned_dofs, comm):
    vals = np.abs(np.asarray(e_complex.x.array)[owned_dofs])
    local = float(vals.max()) if vals.size else 0.0
    return comm.allreduce(local, op=MPI.MAX)


def _outer_box_tags(msh, build_facet_tags, comm):
    """Local copy of `TH-14`'s ``_outer_box_tags``: exterior minus 402, tag 499."""
    fdim = msh.topology.dim - 1
    msh.topology.create_connectivity(fdim, msh.topology.dim)
    n_owned = msh.topology.index_map(fdim).size_local
    exterior = np.asarray(dolfinx.mesh.exterior_facet_indices(msh.topology), dtype=np.int32)
    wall = np.asarray(build_facet_tags.find(BIRDCAGE_IMPLANT_SURFACE_TAG), dtype=np.int32)
    outer = np.setdiff1d(exterior, wall).astype(np.int32)
    tags = dolfinx.mesh.meshtags(
        msh, fdim, outer, np.full(outer.size, OUTER_BOX_TAG, dtype=np.int32)
    )
    census = {
        "exterior": _owned_count(exterior, n_owned, comm),
        "wall": _owned_count(wall, n_owned, comm),
        "wall_exterior": _owned_count(np.intersect1d(wall, exterior), n_owned, comm),
        "outer": _owned_count(outer, n_owned, comm),
        "tag_collision": comm.allreduce(
            int(np.count_nonzero(np.asarray(build_facet_tags.values) == OUTER_BOX_TAG)),
            op=MPI.SUM,
        ),
    }
    return tags, census


def _ball_mean_closed_form(msh, centre, radius):
    """<|E_cf|^2>/E0^2 over the ball at ``centre`` on owned FEM cells (numpy, deg 8).

    Modelled on `IMP-2`'s ``_numpy_quadrature_mean`` with a ball indicator at
    the quadrature points; only cells that exist in this mesh contribute, so
    the rod's hole is excluded the way the FEM ball excludes it.
    """
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
    cells = np.arange(n_owned, dtype=np.int32)
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


def _matrix_quantities(result):
    z = np.asarray(result.z_matrix, dtype=np.complex128)
    s = np.asarray(result.s_matrix, dtype=np.complex128)
    classes = _circulant_classes(z)
    spreads = {k: _class_spread(v) for k, v in classes.items()}
    pooled = _class_spread(np.concatenate([classes["adjacent"], classes["opposite"]]))
    return {
        "z": z,
        "s": s,
        "sigma": np.linalg.svd(s, compute_uv=False),
        "colpow": np.sum(np.abs(s) ** 2, axis=0),
        "recip": float(np.linalg.norm(s - s.T) / np.linalg.norm(s)),
        "spreads": spreads,
        "pooled": pooled,
        "separation": pooled / max(spreads.values()),
        "mirror": [_class_spread([z[a], z[b]]) for a, b in MIRROR_PAIRS],
        "sigma_x": [_class_spread([z[a], z[b]]) for a, b in SIGMA_X_PAIRS],
    }


def _arm(built, label, rod, *, field_hook=None):
    """Sweep one arm (steps 1-5, 9); return scalars and, for the rod arm, the mesh."""
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
    for tip, centre in TIPS.items():
        for d in DRIVES:
            out = mass_averaged_sar(
                e_of[d], sigma=sigma_field, rho=rho_field, center=centre,
                radius=a_1g, quadrature_degree=QUADRATURE_DEGREE,
            )
            sar[(tip, d)] = float(out["averaged_sar_w_per_kg"])
            mass[tip] = float(out["mass_kg"])
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
    }

    if not rod:
        # Step 6: background field on the rod axis, without-rod arm only.
        zs = np.linspace(-0.016, 0.016, AXIS_SAMPLES)
        pts = np.column_stack([np.full_like(zs, ROD_CENTRE[0]), np.zeros_like(zs), zs])
        axis_e = {}
        for pid in PORT_IDS:
            f = result.fields[pid]
            re, ok_r = evaluate_vector_field_parallel(f.e_real, pts)
            im, ok_i = evaluate_vector_field_parallel(f.e_imag, pts)
            assert bool(np.all(ok_r)) and bool(np.all(ok_i)), "axis sample off-mesh"
            axis_e[pid] = np.real(re) + 1j * np.real(im)
        for sense in ("ccw", "cw"):
            axis_e[sense] = sum(w * axis_e[p] for w, p in zip(weights[sense], PORT_IDS))
        out["eta_z"] = {
            d: float(np.sum(np.abs(axis_e[d][:, 2]) ** 2) / np.sum(np.abs(axis_e[d]) ** 2))
            for d in DRIVES
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
        out["sigma_field"] = sigma_field
    if field_hook is not None:
        field_hook(label, rod, msh, cell_tags, e_of, sigma_field, rho_field)
    out["t_arm"] = time.perf_counter() - t_arm
    # Step 7: drop the solved fields and the solver.
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
    for tip in TIPS:
        _say(
            f"    SAR_1g at {tip} tip (W/kg): "
            + "  ".join(f"{d} {a['sar'][(tip, d)]:.6e}" for d in DRIVES)
        )
    _say(
        f"    times: arm {a['t_arm']:.1f} s, sweep {a['t_sweep']:.2f} s, "
        f"1 g SARs {a['t_sar']:.2f} s"
    )


def _pair_spread(x, y):
    return abs(x - y) / (0.5 * (abs(x) + abs(y)))


def _c3(a):
    s = a["sar"]
    lines = []
    worst = 0.0
    for tip in TIPS:
        p24 = _pair_spread(s[(tip, "P2")], s[(tip, "P4")])
        ccw_cw = _pair_spread(s[(tip, "ccw")], s[(tip, "cw")])
        p13 = _pair_spread(s[(tip, "P1")], s[(tip, "P3")])
        worst = max(worst, p24, ccw_cw)
        lines.append(
            f"        {tip}: P2/P4 {p24:.4%}, ccw/cw {ccw_cw:.4%} (each predicted <= "
            f"{PRED_C3_PAIR_MAX:.0%}); P1/P3 {p13:.4%} (predicted "
            f"[{PRED_C3_P1P3[0]:.0%}, {PRED_C3_P1P3[1]:.0%}])"
        )
    for d in DRIVES:
        tb = _pair_spread(s[("top", d)], s[("bottom", d)])
        worst = max(worst, tb)
        lines.append(f"        top/bottom {d}: {tb:.4%}")
    _say(
        f"    (C3) field-level mirror, {a['label']} (C4_COVARIANCE_BAND "
        f"{C4_COVARIANCE_BAND:.0%} shown as a reference only):"
    )
    for line in lines:
        _say(line)
    _say(f"        worst kept-mirror pair {worst:.4%} -> {'<=' if worst <= PRED_C3_PAIR_MAX else '>'} prediction")


def _peak_rss_gib():
    kib = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return MPI.COMM_WORLD.allreduce(float(kib), op=MPI.SUM) / 1024.0**2


def _case(phantom_resolution, *, field_hook=None):
    """Without-rod arm first, freed, then the rod arm (§9 item 174 steps 1-11)."""
    comm = MPI.COMM_WORLD
    t_case = time.perf_counter()
    window = "W1" if phantom_resolution is None else "W2"
    kw = {} if phantom_resolution is None else {"phantom_resolution": phantom_resolution}

    built_off = _build(True, **kw)
    off = _arm(built_off, f"{window} without rod", rod=False, field_hook=field_hook)
    del built_off
    gc.collect()

    built_rod = _build_rod(**kw)
    rod_cells = int(built_rod[0].topology.index_map(3).size_global)
    _say(f"\n[IMP-3 {window}] rod-arm build: {rod_cells} cells (guard {ROD_CELL_GUARD_W2})")
    if phantom_resolution is not None:
        assert rod_cells < ROD_CELL_GUARD_W2, (rod_cells, ROD_CELL_GUARD_W2)
    rod = _arm(built_rod, f"{window} with rod", rod=True, field_hook=field_hook)

    # Closed form on the rod-arm mesh.
    lz = prolate_spheroid_depolarisation(
        math.sqrt(1.0 - (SPHEROID_SEMI_MINOR / SPHEROID_SEMI_MAJOR) ** 2)
    )
    m_cf = {}
    for tip, centre in TIPS.items():
        m_cf[tip], _vol = _ball_mean_closed_form(rod["case"]["mesh"], centre, off["a_1g"])

    case = {
        "window": window,
        "off": off,
        "rod": rod,
        "rod_cells": rod_cells,
        "lz": lz,
        "m_cf": m_cf,
    }
    case["t_case"] = time.perf_counter() - t_case
    return case


def _report(case):
    off, rod = case["off"], case["rod"]
    w = case["window"]
    _print_arm(off)
    _print_arm(rod)
    _say(
        f"    (E) tag-402 dofs (owned, reduced) {rod['wall_dofs']}; max |E| per drive "
        + ", ".join(f"{p} {v:.3e}" for p, v in rod["wall_max"].items())
        + f" (bound {PEC_WALL_MAX_E:g})"
    )
    _say(
        f"    (S) rod-arm separation {rod['mq']['separation']:.1f}x (predicted "
        f"[{PRED_SEPARATION_ROD[0]:g}, {PRED_SEPARATION_ROD[1]:g}]x, "
        f"{_in(rod['mq']['separation'], PRED_SEPARATION_ROD)}); rod-arm mirror pairs "
        f"max {max(rod['mq']['mirror']):.4%} (predicted [{PRED_MIRROR_SPREAD_ROD[0]:.2%}, "
        f"{PRED_MIRROR_SPREAD_ROD[1]:.2%}])"
    )
    # C2
    ds = rod["mq"]["s"] - off["mq"]["s"]
    dcol = rod["mq"]["colpow"][0] - off["mq"]["colpow"][0]
    _say(f"[IMP-3 {w}] (C2) what the 4x4 sees of the rod (predicted, printed):")
    for pid in PORT_IDS:
        p0 = off["shares"][pid]["phantom"]
        p1 = rod["shares"][pid]["phantom"]
        _say(
            f"        {pid}: P_phantom {p0:.9e} -> {p1:.9e} W, dP/P {(p1 - p0) / p0:+.4%} "
            f"(predicted |dP/P| in [{PRED_DP_PHANTOM[0]:.0%}, {PRED_DP_PHANTOM[1]:.0%}], "
            f"{_in(abs((p1 - p0) / p0), PRED_DP_PHANTOM)})"
        )
    _say(
        f"        d sum_i |S_i1|^2 = {dcol:+.3e} (predicted |.| in [{PRED_D_COLPOWER[0]:g}, "
        f"{PRED_D_COLPOWER[1]:g}], {_in(abs(dcol), PRED_D_COLPOWER)}); max |dS| = "
        f"{np.max(np.abs(ds)):.3e} (predicted [{PRED_MAX_DS[0]:g}, {PRED_MAX_DS[1]:g}], "
        f"{_in(float(np.max(np.abs(ds))), PRED_MAX_DS)})"
    )
    _c3(off)
    _c3(rod)
    # Ratio of record and the closed form.
    k_mag = 2.0 * np.pi * float(STEP2D_FREQUENCY_HZ) * math.sqrt(
        MU_0 * EPSILON_0 * abs(SALINE_EPSILON_R - 1j * SALINE_SIGMA
                               / (2.0 * np.pi * float(STEP2D_FREQUENCY_HZ) * EPSILON_0))
    )
    c_eq = math.sqrt(SPHEROID_SEMI_MAJOR**2 - SPHEROID_SEMI_MINOR**2)
    _say(
        f"[IMP-3 {w}] ratio Q = SAR_1g,tip(rod)/SAR_1g,tip(without) (predicted "
        f"[{PRED_Q[0]:g}, {PRED_Q[1]:g}]); L_z(16) = {case['lz']:.6f}; "
        f"|k| = {k_mag:.3f} 1/m, (|k| c_eq)^2 = {(k_mag * c_eq) ** 2:.4%}"
    )
    _say(
        f"    ball masses: without {off['mass']['top'] * 1e3:.6f} / "
        f"{off['mass']['bottom'] * 1e3:.6f} g, with rod {rod['mass']['top'] * 1e3:.6f} / "
        f"{rod['mass']['bottom'] * 1e3:.6f} g (predicted 1.000 / {PRED_MASS_ROD_G} g)"
    )
    for tip in TIPS:
        _say(
            f"    M_cf({tip}) = {case['m_cf'][tip]:.6f} (predicted [{PRED_M_CF[0]:g}, "
            f"{PRED_M_CF[1]:g}], {_in(case['m_cf'][tip], PRED_M_CF)})"
        )
    for d in DRIVES:
        eta = off["eta_z"][d]
        for tip in TIPS:
            q = rod["sar"][(tip, d)] / off["sar"][(tip, d)]
            est = eta * case["m_cf"][tip] + (1.0 - eta)
            _say(
                f"    {d:>3s} {tip:>6s}: Q = {q:.6f} ({_in(q, PRED_Q)}); eta_z = {eta:.6f} "
                f"({_in(eta, PRED_ETA_Z)}); estimate {est:.6f}; Q/estimate "
                f"{q / est:.6f} ({_in(q / est, PRED_Q_OVER_ESTIMATE)})"
            )
    _say(
        f"[IMP-3 {w}] cells without {off['cells']} / with {rod['cells']}; case "
        f"{case['t_case']:.1f} s; summed peak RSS {_peak_rss_gib():.2f} GiB at -n "
        f"{MPI.COMM_WORLD.size}"
    )


def _assert_identities(case):
    off, rod = case["off"], case["rod"]
    for arm in (off, rod):
        mq = arm["mq"]
        assert mq["recip"] <= RECIPROCITY_BAND, (arm["label"], mq["recip"])
        assert float(mq["sigma"][0]) <= 1.0 + PASSIVITY_SIGMA_TOLERANCE, (arm["label"], mq["sigma"])
        assert np.all(mq["colpow"] <= 1.0), (arm["label"], mq["colpow"])
        for pid in PORT_IDS:
            ex = arm["shares"][pid]
            resid = abs(ex["p_src"] - ex["p_vol"] - ex["sheet_field_total"]) / abs(ex["p_src"])
            assert resid <= DISCRETE_IDENTITY_RTOL, (arm["label"], pid, resid)
    for v in rod["mq"]["mirror"]:
        assert v <= ADJACENT_SPREAD_BAND, rod["mq"]["mirror"]
    assert rod["wall_dofs"] > 0
    for pid, v in rod["wall_max"].items():
        assert v <= PEC_WALL_MAX_E, (pid, v)
    assert off["mq"]["separation"] >= POOLED_SEPARATION_FLOOR, off["mq"]["separation"]
    a = off["a_1g"]
    assert ROD_CENTRE[0] + a < PHANTOM_RADIUS, (ROD_CENTRE[0] + a, PHANTOM_RADIUS)
    assert 0.016 + a < 0.5 * PHANTOM_HEIGHT, (0.016 + a, 0.5 * PHANTOM_HEIGHT)


def _control_c1(case):
    """(C1) the rod wall left natural, P1 only, on a fresh solver."""
    comm = MPI.COMM_WORLD
    rod = case["rod"]
    sweep = rod["case"]
    msh = sweep["mesh"]
    outer, census = _outer_box_tags(msh, rod["built"][2], comm)
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
    wall_max = _wall_max_e(solved["fields"].e_complex, rod["owned_wall"], comm)
    a_1g = case["off"]["a_1g"]
    sar_nat = {}
    for tip, centre in TIPS.items():
        out = mass_averaged_sar(
            solved["fields"].e_complex, sigma=solved["fields"].sigma_field,
            rho=rod["rho_field"], center=centre, radius=a_1g,
            quadrature_degree=QUADRATURE_DEGREE,
        )
        sar_nat[tip] = float(out["averaged_sar_w_per_kg"])
    _say(
        f"    (C1) natural wall, P1: (P) rel {resid:.3e} (band {DISCRETE_IDENTITY_RTOL:g}); "
        f"max |E| on 402 dofs {wall_max:.3e} (> {PEC_WALL_MAX_E:g} = knob wired); "
        f"phantom {ex['phantom']:.9e} W; {time.perf_counter() - t0:.1f} s"
    )
    for tip in TIPS:
        r_off = sar_nat[tip] / case["off"]["sar"][(tip, "P1")]
        r_pec = rod["sar"][(tip, "P1")] / sar_nat[tip]
        _say(
            f"        {tip}: SAR_1g natural {sar_nat[tip]:.6e} W/kg; natural/without "
            f"{r_off:.6f} (predicted [{PRED_C1_NATURAL_OVER_OFF[0]:g}, "
            f"{PRED_C1_NATURAL_OVER_OFF[1]:g}], {_in(r_off, PRED_C1_NATURAL_OVER_OFF)}); "
            f"PEC/natural {r_pec:.6f} (predicted [{PRED_C1_PEC_OVER_NATURAL[0]:g}, "
            f"{PRED_C1_PEC_OVER_NATURAL[1]:g}], {_in(r_pec, PRED_C1_PEC_OVER_NATURAL)})"
        )
    return {"resid": resid, "wall_max": wall_max}


@complex_only
def test_imp3_default_phantom():
    """W1 (the 🟡 window): the record meshes, (R0), identities, C1."""
    case = _case(None)
    off, rod = case["off"], case["rod"]
    rel = np.abs(off["mq"]["z"] - IMP3_FLAG_OFF_Z_64MHZ) / np.abs(IMP3_FLAG_OFF_Z_64MHZ)
    _report(case)
    _say(
        f"[IMP-3 W1] (A0) without {off['cells']} vs {SOLID_CELL_RECORD}, with "
        f"{rod['cells']} vs {ROD_CELL_RECORD_W1} (band {CELL_COUNT_BAND:.0%}); "
        f"(R0) max entrywise |Z - Z_rec|/|Z_rec| = {float(rel.max()):.3e} "
        f"(band {REUSE_REPRODUCTION_RTOL:g})"
    )
    c1 = _control_c1(case)
    assert abs(off["cells"] / SOLID_CELL_RECORD - 1.0) <= CELL_COUNT_BAND
    assert abs(rod["cells"] / ROD_CELL_RECORD_W1 - 1.0) <= CELL_COUNT_BAND
    assert float(rel.max()) <= REUSE_REPRODUCTION_RTOL, rel
    _assert_identities(case)
    assert c1["resid"] <= DISCRETE_IDENTITY_RTOL, c1
    assert c1["wall_max"] > PEC_WALL_MAX_E, c1


def _run_one_gram_rung(*, field_hook=None):
    """W2 body (the 1 g rung, identities, the ratio of record); returns ``case``.

    ``field_hook`` (default None, prints unchanged) is called per arm just before
    the solved fields are freed; EX-83 uses it to keep DG0 maps (§9 item 191).
    """
    case = _case(PHANTOM_RESOLUTION_1G_RUNG, field_hook=field_hook)
    off, rod = case["off"], case["rod"]
    _report(case)
    _say(
        f"[IMP-3 W2] (A0) without {off['cells']} vs {ONE_GRAM_RUNG_CELL_RECORD} (band "
        f"{CELL_COUNT_BAND:.0%}); with {rod['cells']} — first record, predicted "
        f"[{ROD_CELL_PREDICTED_W2[0]}, {ROD_CELL_PREDICTED_W2[1]}] "
        f"({_in(rod['cells'], ROD_CELL_PREDICTED_W2)}), guard {ROD_CELL_GUARD_W2}"
    )
    assert abs(off["cells"] / ONE_GRAM_RUNG_CELL_RECORD - 1.0) <= CELL_COUNT_BAND
    assert rod["cells"] < ROD_CELL_GUARD_W2
    _assert_identities(case)
    return case


@complex_only
def test_imp3_one_gram_rung():
    """W2 (the ✅ window): the 1 g rung, identities, the ratio of record."""
    _run_one_gram_rung()
