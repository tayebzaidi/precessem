"""Example (`EX-122`): a PEC oblate disc in the coil-driven phantom -- rim 1 g SAR, PEC vs natural wall.

`imp:8` puts this disc (a = 6 mm, c = 1.2 mm, aspect 5, normal x-hat) in an
imposed quasi-static field with no coil; `imp:3` puts a rod in the coil. This
example puts the disc in the coil: the loaded F-small 4-leg birdcage at 64 MHz
(untuned sigma = 800 S/m solid coil, 50 Ohm lumped sheets, degree 1), through
`IMP-3` step 9's own builders on the W1 record meshes.

**Import, never restate.** ``_case_disc(None)`` (both arms), ``_report`` (every
gate line, in W1's format) and ``_control_c1_disc`` are called unchanged from
``tests/validation/test_imp3_disc_birdcage.py``. Every band and record is
taken through that module's namespace (``g9``). No ``test_*`` name is used.

**Anchors** (asserted at the imported bands): (A0) cell counts, (R0) the
without-disc Z against the 64 MHz record, (i) reciprocity, (ii) passivity,
(P) the power identity on every drive of both arms, (M) the disc-arm mirror
pairs, (E) n x E = 0 on the tag-402 dofs, (S) the without-disc separation,
containment of the 1 g ball. Example-side: (F-pec) a fresh P1 PEC solve closes
(P) and (E); (L) the DG0 phantom-loss map sums to that solve's phantom power,
for the PEC and the natural wall; (N) the natural-wall problem closes (P) and
its census reads exterior = outer + 402.

**Negative controls** (asserted; each backed by a W1 measurement of the same
comparison): (C1) the gate's, called unchanged; (C1-r) PEC / natural 1 g SAR
> 1 at both rim points; (Cq) Q = SAR_1g(disc) / SAR_1g(without) > 1 on P1, P3,
ccw, cw at both rim points. P2 / P4 are printed only (their rim SARs are
mesh-noisy on the W1 mesh).

**Printed only.** The rim-near PEC / natural local-SAR ratio, a rim |E|
profile beside the shape of the closed form, and times.

**Scope.** One aspect-5 PEC disc, normal x-hat, centre (0.015, 0, 0), in the
untuned solid-coil F-small 4x4 with 50 Ohm sheets, 64 MHz, degree 1, W1's
default phantom, ``-n 2``. W1 is the coarse record mesh; W2 (the 1 g rung) is
the ratio of record and is cited, not run. Not claimed: the accuracy of the rim
ratio (`ANS-10`), absolute or pointwise SAR, convergence, tissue, heating,
ISO 10974 / F2182, a tuned or copper coil, a square-edged plate, an external
check.

Run it through the example runner (the ``imp:`` group sources the complex
build)::

    ./scripts/run_examples.sh -e imp:9 -n 2 -t 590

Output lands in ``examples/implants/paraview_output/``:
``implants_09_pec_disc_in_birdcage_rim_sar_combined.xdmf`` and
``implants_09_pec_disc_in_birdcage_rim_sar_rim_profile.csv``.
"""

from __future__ import annotations

import sys
import time
from dataclasses import replace
from pathlib import Path

import numpy as np
import ufl
from mpi4py import MPI

import dolfinx
from dolfinx import default_scalar_type, fem, la

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from fem_em_solver.io.paraview_utils import (  # noqa: E402
    adopt_host_ownership,
    write_xdmf_with_tags,
)
from fem_em_solver.post.evaluation import evaluate_vector_field_parallel  # noqa: E402
from fem_em_solver.post.sar import mass_averaged_sar  # noqa: E402
from fem_em_solver.post.setup_figure import write_setup_figure  # noqa: E402
from fem_em_solver.utils.analytical import pec_oblate_spheroid_field  # noqa: E402

import tests.validation.test_imp3_disc_birdcage as g9  # noqa: E402

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"
BASENAME = "implants_09_pec_disc_in_birdcage_rim_sar"
W1_LOG = "20261013T213344Z_IMP-3-step9-W1.log"

# Records, printed beside the readings and never asserted by equality.
REC = {
    "cells": ("116 085 / 150 140", f"{W1_LOG}:3764"),
    "R0": (2.983e-10, f"{W1_LOG}:3765"),
    "P_disc": (2.988e-15, f"{W1_LOG}:3710"),
    "P_phantom_disc": (1.562944066e-06, f"{W1_LOG}:3710"),
    "E": ("6 669 dofs / 0.000e+00", f"{W1_LOG}:3719"),
    "M": (0.0425, f"{W1_LOG}:3720"),
    "S": (671.1, f"{W1_LOG}:3684"),
    "Q_P1": ((1.637623, 1.646293), f"{W1_LOG}:3751-3752"),
    "Q_min": (1.531498, f"{W1_LOG}:3762"),
    "C1_P": (3.154e-15, f"{W1_LOG}:3768"),
    "C1_wall": (2.627e-04, f"{W1_LOG}:3768"),
    "C1_pec_over_nat": ((1.536348, 1.551662), f"{W1_LOG}:3769-3770"),
    "W2_Q_P1": ((1.592094, 1.593359), "20261013T213816Z_IMP-3-step9-W2.log:3862-3863"),
    "W2_Q_over_est": ((0.961546, 0.964433), "20261013T213816Z_IMP-3-step9-W2.log:3862-3863"),
}

# Predictions, written before window 3 (printed, never asserted).
PRED_FRESH_REPRO = 1.0e-8
PRED_MAP_RATIO_MAX_MIN = 1.54  # > the 1 g ratio
PRED_PROFILE_FIRST_PEC_ABOVE_NAT = True
PRED_PROFILE_VS_CF = (0.5, 2.0)

NEAR_RADIUS_M = 1.0e-3
N_PROFILE = 12
PROFILE_Z0_OFF_M = 0.25e-3
PROFILE_Z1_OFF_M = 6.0e-3
CQ_DRIVES = ("P1", "P3", "ccw", "cw")

PHANTOM_CELL_TAG = g9.PHANTOM_CELL_TAG
RIM_POINTS = g9.RIM_POINTS
DRIVES = g9.DRIVES


def _cellwise(msh, form_expr):
    """Per-cell integral of ``form_expr`` against a DG0 test function (owned+ghost)."""
    dg0 = fem.functionspace(msh, ("DG", 0))
    v = ufl.TestFunction(dg0)
    b = fem.assemble_vector(fem.form(form_expr(v)))
    b.scatter_reverse(la.InsertMode.add)
    b.scatter_forward()
    return dg0, np.real(np.array(b.array))


def _ball_mask(msh, dg0, centre, radius):
    """DG0 0/1 field: cell midpoint inside the ball (numpy, no UFL comparison)."""
    gdm = msh.geometry.dofmap
    xg = msh.geometry.x
    mid = xg[gdm].mean(axis=1)
    inside = np.sum((mid - np.asarray(centre)) ** 2, axis=1) < radius**2
    f = fem.Function(dg0)
    f.x.array[:] = 0.0
    f.x.array[dg0.dofmap.list[:, 0]] = inside.astype(float)
    f.x.scatter_forward()
    return f


def _maps(msh, cell_tags, solved):
    """DG0 maps of one solved P1 field: loss vector, |E|, local SAR, owned loss sum."""
    comm = msh.comm
    e = solved["fields"].e_complex
    sig = ufl.real(solved["fields"].sigma_field)
    dx = ufl.Measure("dx", domain=msh, subdomain_data=cell_tags)
    dxa = ufl.Measure("dx", domain=msh)
    n_loc = msh.topology.index_map(msh.topology.dim).size_local
    dg0, vol = _cellwise(msh, lambda v: ufl.conj(v) * dxa)
    _, loss = _cellwise(
        msh, lambda v: 0.5 * sig * ufl.inner(e, e) * ufl.conj(v) * dx(PHANTOM_CELL_TAG)
    )
    _, e2 = _cellwise(msh, lambda v: ufl.inner(e, e) * ufl.conj(v) * dxa)
    loss_sum = float(comm.allreduce(float(np.sum(loss[:n_loc])), op=MPI.SUM))
    return {
        "dg0": dg0,
        "vol": vol,
        "loss": loss,
        "emag": np.sqrt(np.maximum(e2 / vol, 0.0)),
        "lsar": loss / (g9.PHANTOM_RHO_KG_PER_M3 * vol),
        "loss_sum": loss_sum,
    }


def _profile_points():
    a = g9.SPHEROID_SEMI_MAJOR
    c = g9.SPHEROID_CENTRE
    z = np.linspace(c[2] + a + PROFILE_Z0_OFF_M, c[2] + a + PROFILE_Z1_OFF_M, N_PROFILE)
    return np.column_stack([np.full(N_PROFILE, c[0]), np.full(N_PROFILE, c[1]), z])


def _sar_1g(solved, rho_field, a_1g):
    out = {}
    for rim, centre in RIM_POINTS.items():
        r = mass_averaged_sar(
            solved["fields"].e_complex, sigma=solved["fields"].sigma_field,
            rho=rho_field, center=centre, radius=a_1g,
            quadrature_degree=g9.QUADRATURE_DEGREE,
        )
        out[rim] = float(r["averaged_sar_w_per_kg"])
    return out


def main() -> None:
    comm = MPI.COMM_WORLD
    rank0 = comm.rank == 0
    started = time.perf_counter()
    say = g9._say
    timers = {}

    if not np.issubdtype(np.dtype(default_scalar_type), np.complexfloating):
        raise RuntimeError(
            "this example needs the complex DolfinX build "
            "(source /usr/local/bin/dolfinx-complex-mode); the runner's `imp:` "
            "group sources it automatically"
        )
    assert comm.size == 2, comm.size
    assert len(RIM_POINTS) == 2, RIM_POINTS
    assert abs(g9.SPHEROID_SEMI_MAJOR / g9.SPHEROID_SEMI_MINOR - 5.0) <= 1e-12

    say("=" * 72)
    say("EX-122 -- a PEC oblate disc in the coil-driven phantom: rim 1 g SAR, PEC vs natural")
    say("=" * 72)
    say(
        "\n[scope] W1 record mesh (the coarse mesh; W2 is the ratio of record, cited not run), "
        "degree 1, 64 MHz, untuned solid coil, 50 Ohm sheets, one aspect-5 disc. No accuracy of "
        "the rim ratio (ANS-10), no absolute SAR, no heating, no compliance."
    )

    # ---- the gate's own body and report, unchanged ------------------------------
    t0 = time.perf_counter()
    case = g9._case_disc(None)
    timers["gate case"] = time.perf_counter() - t0
    off, disc = case["off"], case["disc"]
    rel_z = np.abs(off["mq"]["z"] - g9.IMP3_FLAG_OFF_Z_64MHZ) / np.abs(g9.IMP3_FLAG_OFF_Z_64MHZ)
    g9._report(case)
    say(
        f"[EX-122] (A0) without {off['cells']} vs {g9.SOLID_CELL_RECORD}, with "
        f"{disc['cells']} vs {g9.SPHEROID_HOLE_CELL_RECORD} (band {g9.CELL_COUNT_BAND:.0%}); "
        f"record {REC['cells'][0]} ({REC['cells'][1]}); (R0) max entrywise "
        f"{float(rel_z.max()):.3e} (band {g9.REUSE_REPRODUCTION_RTOL:g}; record "
        f"{REC['R0'][0]:.3e}, {REC['R0'][1]})"
    )

    # ---- (C1), the gate's, unchanged --------------------------------------------
    t0 = time.perf_counter()
    c1 = g9._control_c1_disc(case)
    timers["C1"] = time.perf_counter() - t0

    # ---- anchors, asserted at the imported bands --------------------------------
    say("\n[anchors] asserted at the imported bands (digits beside "
        f"{W1_LOG}):")
    a_cells_off = abs(off["cells"] / g9.SOLID_CELL_RECORD - 1.0)
    a_cells_disc = abs(disc["cells"] / g9.SPHEROID_HOLE_CELL_RECORD - 1.0)
    say(f"  (A0) |cells/record - 1|: without {a_cells_off:.3e}, with {a_cells_disc:.3e} "
        f"(band {g9.CELL_COUNT_BAND:.0%})")
    assert a_cells_off <= g9.CELL_COUNT_BAND, a_cells_off
    assert a_cells_disc <= g9.CELL_COUNT_BAND, a_cells_disc
    say(f"  (R0) {float(rel_z.max()):.3e} (band {g9.REUSE_REPRODUCTION_RTOL:g})")
    assert float(rel_z.max()) <= g9.REUSE_REPRODUCTION_RTOL, rel_z
    for arm in (off, disc):
        mq = arm["mq"]
        say(f"  (i)  {arm['label']}: {mq['recip']:.3e} (band {g9.RECIPROCITY_BAND:g}); "
            f"(ii) sigma_max {float(mq['sigma'][0]):.9f} (<= 1 + {g9.PASSIVITY_SIGMA_TOLERANCE:g}), "
            f"max column power {float(np.max(mq['colpow'])):.9f} (<= 1)")
        assert mq["recip"] <= g9.RECIPROCITY_BAND, mq["recip"]
        assert float(mq["sigma"][0]) <= 1.0 + g9.PASSIVITY_SIGMA_TOLERANCE, mq["sigma"]
        assert np.all(mq["colpow"] <= 1.0), mq["colpow"]
        worst = 0.0
        for pid in g9.PORT_IDS:
            ex = arm["shares"][pid]
            resid = abs(ex["p_src"] - ex["p_vol"] - ex["sheet_field_total"]) / abs(ex["p_src"])
            worst = max(worst, resid)
            assert resid <= g9.DISCRETE_IDENTITY_RTOL, (arm["label"], pid, resid)
        say(f"  (P)  {arm['label']}: worst over drives {worst:.3e} (band "
            f"{g9.DISCRETE_IDENTITY_RTOL:g}; W1 <= 7.0e-15; disc P1 record "
            f"{REC['P_disc'][0]:.3e}, {REC['P_disc'][1]})")
    mx = max(disc["mq"]["mirror"])
    say(f"  (M)  disc-arm sigma_y pairs max {mx:.4%} (band {g9.ADJACENT_SPREAD_BAND:.2%}; "
        f"record {REC['M'][0]:.4f} %, {REC['M'][1]})")
    for v in disc["mq"]["mirror"]:
        assert v <= g9.ADJACENT_SPREAD_BAND, disc["mq"]["mirror"]
    say(f"  (E)  tag-402 dofs {disc['wall_dofs']}; max |E| "
        + ", ".join(f"{p} {v:.3e}" for p, v in disc["wall_max"].items())
        + f" (bound {g9.PEC_WALL_MAX_E:g}; record {REC['E'][0]}, {REC['E'][1]})")
    assert disc["wall_dofs"] > 0
    for pid, v in disc["wall_max"].items():
        assert v <= g9.PEC_WALL_MAX_E, (pid, v)
    say(f"  (S)  without-disc separation {off['mq']['separation']:.1f}x (floor "
        f"{g9.POOLED_SEPARATION_FLOOR:g}; record {REC['S'][0]}x, {REC['S'][1]})")
    assert off["mq"]["separation"] >= g9.POOLED_SEPARATION_FLOOR, off["mq"]["separation"]
    a_1g = off["a_1g"]
    say(f"  containment: x0 + a_1g = {g9.SPHEROID_CENTRE[0] + a_1g:.6f} < R = "
        f"{g9.PHANTOM_RADIUS:.6f}; a + a_1g = {g9.SPHEROID_SEMI_MAJOR + a_1g:.6f} < H/2 = "
        f"{0.5 * g9.PHANTOM_HEIGHT:.6f}")
    assert g9.SPHEROID_CENTRE[0] + a_1g < g9.PHANTOM_RADIUS
    assert g9.SPHEROID_SEMI_MAJOR + a_1g < 0.5 * g9.PHANTOM_HEIGHT

    # ---- (F-pec): a fresh PEC solve ---------------------------------------------
    sweep = disc["case"]
    msh = sweep["mesh"]
    cell_tags = sweep["cell_tags"]
    rho_field = disc["rho_field"]
    t0 = time.perf_counter()
    pec = g9._solve_driven(sweep, "P1")
    ex_pec = g9._exact_shares(sweep, pec)
    resid_pec = abs(ex_pec["p_src"] - ex_pec["p_vol"] - ex_pec["sheet_field_total"]) / abs(ex_pec["p_src"])
    wall_pec = g9._wall_max_e(pec["fields"].e_complex, disc["owned_wall"], comm)
    sar_pec = _sar_1g(pec, rho_field, a_1g)
    timers["fresh PEC solve + shares + SAR"] = time.perf_counter() - t0
    say("\n[example-side, asserted] (F-pec) fresh P1 PEC solve on a fresh solver:")
    say(f"  (P) rel {resid_pec:.3e} (band {g9.DISCRETE_IDENTITY_RTOL:g}); (E) max |E| on 402 dofs "
        f"{wall_pec:.3e} (bound {g9.PEC_WALL_MAX_E:g})")
    assert resid_pec <= g9.DISCRETE_IDENTITY_RTOL, resid_pec
    assert wall_pec <= g9.PEC_WALL_MAX_E, wall_pec
    p_ref = disc["shares"]["P1"]["phantom"]
    say(f"  printed only: |P_phantom,fresh / sweep - 1| = {abs(ex_pec['phantom'] / p_ref - 1.0):.3e} "
        f"(predicted <= {PRED_FRESH_REPRO:g}); fresh {ex_pec['phantom']:.9e} W vs sweep {p_ref:.9e} W "
        f"(record {REC['P_phantom_disc'][0]:.9e} W, {REC['P_phantom_disc'][1]})")
    for rim in RIM_POINTS:
        rr = abs(sar_pec[rim] / disc["sar"][(rim, "P1")] - 1.0)
        say(f"  printed only: {rim} |SAR_1g,fresh / sweep - 1| = {rr:.3e} (predicted <= "
            f"{PRED_FRESH_REPRO:g}); fresh {sar_pec[rim]:.6e} W/kg")

    # ---- (N): the natural-wall solve --------------------------------------------
    t0 = time.perf_counter()
    outer, census = g9._outer_box_tags(msh, disc["built"][2], comm)
    problem = replace(sweep["problem"], facet_tags=outer, pec_facet_tags=(g9.OUTER_BOX_TAG,))
    nat_case = dict(sweep, problem=problem)
    nat = g9._solve_driven(nat_case, "P1")
    ex_nat = g9._exact_shares(nat_case, nat)
    resid_nat = abs(ex_nat["p_src"] - ex_nat["p_vol"] - ex_nat["sheet_field_total"]) / abs(ex_nat["p_src"])
    sar_nat = _sar_1g(nat, rho_field, a_1g)
    timers["natural solve + shares + SAR"] = time.perf_counter() - t0
    say("\n[example-side, asserted] (N) the natural-wall solve:")
    say(f"  census: exterior {census['exterior']} = outer {census['outer']} + 402 "
        f"{census['wall_exterior']} (402 total {census['wall']}); tag 499 pre-existing on "
        f"{census['tag_collision']} facets")
    assert census["tag_collision"] == 0, census
    assert census["exterior"] == census["outer"] + census["wall_exterior"], census
    assert census["wall_exterior"] == census["wall"] > 0, census
    say(f"  (P) rel {resid_nat:.3e} (band {g9.DISCRETE_IDENTITY_RTOL:g}); C1 resid {c1['resid']:.3e}")
    assert resid_nat <= g9.DISCRETE_IDENTITY_RTOL, resid_nat

    # ---- (L): loss-map sums ------------------------------------------------------
    t0 = time.perf_counter()
    m_pec = _maps(msh, cell_tags, pec)
    m_nat = _maps(msh, cell_tags, nat)
    timers["DG0 maps"] = time.perf_counter() - t0
    say("\n[example-side, asserted] (L) owned DG0 phantom-loss sums vs the solve's phantom power:")
    for name, mp, ex in (("PEC", m_pec, ex_pec), ("natural", m_nat, ex_nat)):
        rel = abs(mp["loss_sum"] - ex["phantom"]) / abs(ex["phantom"])
        say(f"  {name}: map {mp['loss_sum']:.9e} W vs shares {ex['phantom']:.9e} W, rel {rel:.3e} "
            f"(band {g9.DISCRETE_IDENTITY_RTOL:g})")
        assert rel <= g9.DISCRETE_IDENTITY_RTOL, (name, rel)

    # ---- negative controls --------------------------------------------------------
    say("\n[control C1] the gate's, called unchanged: "
        f"resid {c1['resid']:.3e} <= {g9.DISCRETE_IDENTITY_RTOL:g}; wall max |E| "
        f"{c1['wall_max']:.3e} > {g9.PEC_WALL_MAX_E:g} (record {REC['C1_wall'][0]:.3e}, "
        f"{REC['C1_wall'][1]}) -- ASSERTED")
    assert c1["resid"] <= g9.DISCRETE_IDENTITY_RTOL, c1
    assert c1["wall_max"] > g9.PEC_WALL_MAX_E, c1

    say("\n[control C1-r] SAR_1g(PEC, P1) / SAR_1g(natural, P1) > 1 at both rim points -- ASSERTED:")
    for k, rim in enumerate(RIM_POINTS):
        ratio = sar_pec[rim] / sar_nat[rim]
        rec = REC["C1_pec_over_nat"][0][k]
        say(f"  {rim}: natural {sar_nat[rim]:.6e} W/kg, PEC {sar_pec[rim]:.6e} W/kg, ratio "
            f"{ratio:.6f} (gate prediction [{g9.PRED_C1_PEC_OVER_NATURAL[0]:g}, "
            f"{g9.PRED_C1_PEC_OVER_NATURAL[1]:g}], "
            f"{g9._in(ratio, g9.PRED_C1_PEC_OVER_NATURAL)}; record {rec:.6f}, "
            f"{REC['C1_pec_over_nat'][1]}; fresh-vs-record {ratio / rec - 1.0:+.2e})")
        assert ratio > 1.0, (rim, ratio)

    say("\n[control Cq] Q = SAR_1g,rim(disc) / SAR_1g,rim(without); P1, P3, ccw, cw > 1 ASSERTED, "
        "P2 / P4 printed only:")
    qmin = None
    for d in DRIVES:
        row = []
        for rim in RIM_POINTS:
            q = disc["sar"][(rim, d)] / off["sar"][(rim, d)]
            row.append(f"{rim} {q:.6f}")
            if d in CQ_DRIVES:
                assert q > 1.0, (d, rim, q)
                if qmin is None or q < qmin[0]:
                    qmin = (q, d, rim)
        say(f"  {d:>3s}: " + "  ".join(row) + ("" if d in CQ_DRIVES else "  (printed only)"))
    say(f"  min Q over the 8 asserted pairs {qmin[0]:.6f} ({qmin[1]} {qmin[2]}); record minimum "
        f"{REC['Q_min'][0]:.6f} ({REC['Q_min'][1]}); P1 record {REC['Q_P1'][0]} "
        f"({REC['Q_P1'][1]})")
    say(f"  ratio of record (W2, cited not run): P1 Q {REC['W2_Q_P1'][0]}, Q/estimate "
        f"{REC['W2_Q_over_est'][0]} ({REC['W2_Q_P1'][1]})")

    # ---- printed: rim-near map ratio ---------------------------------------------
    say("\n[printed, no band] rim-near map: gel cells with midpoint within "
        f"{NEAR_RADIUS_M * 1e3:.1f} mm of a rim point, local SAR PEC / natural:")
    tdim = msh.topology.dim
    n_owned = msh.topology.index_map(tdim).size_local
    idx = np.asarray(cell_tags.indices)
    val = np.asarray(cell_tags.values)
    gel = np.sort(idx[(idx < n_owned) & (val == PHANTOM_CELL_TAG)].astype(np.int32))
    mid = dolfinx.mesh.compute_midpoints(msh, tdim, gel)
    dofs = m_pec["dg0"].dofmap.list[gel, 0]
    nat_loss = m_nat["loss"][dofs]
    ok_nat = nat_loss > 0.0
    ratio_gel = np.zeros(gel.size)
    ratio_gel[ok_nat] = m_pec["lsar"][dofs][ok_nat] / m_nat["lsar"][dofs][ok_nat]
    map_max = {}
    for rim, centre in RIM_POINTS.items():
        near = (np.linalg.norm(mid - np.asarray(centre), axis=1) < NEAR_RADIUS_M) & ok_nat
        n = comm.allreduce(int(near.sum()), op=MPI.SUM)
        vloc = m_pec["vol"][dofs][near]
        wsum = comm.allreduce(float(np.sum(ratio_gel[near] * vloc)), op=MPI.SUM)
        vsum = comm.allreduce(float(np.sum(vloc)), op=MPI.SUM)
        mloc = float(ratio_gel[near].max()) if near.any() else 0.0
        mmax = comm.allreduce(mloc, op=MPI.MAX)
        map_max[rim] = mmax
        say(f"  {rim}: {n} cells, max ratio {mmax:.4f}, volume-weighted mean {wsum / vsum:.4f} "
            f"(predicted max > {PRED_MAP_RATIO_MAX_MIN}, {'IN' if mmax > PRED_MAP_RATIO_MAX_MIN else 'OUT'}; "
            "no ceiling claimed)")

    # ---- printed: rim profile -----------------------------------------------------
    pts = _profile_points()
    prof = {}
    for name, sol in (("pec", pec), ("natural", nat)):
        vals, ok = evaluate_vector_field_parallel(sol["fields"].e_complex, pts)
        assert bool(np.all(ok)), (name, ok)
        prof[name] = np.sqrt(np.sum(np.abs(np.asarray(vals)) ** 2, axis=1))
    d = pts - np.asarray(g9.SPHEROID_CENTRE)
    p_l = np.column_stack([d[:, 2], d[:, 1], d[:, 0]])
    cf = pec_oblate_spheroid_field(
        p_l, drive="inplane", semi_major=g9.SPHEROID_SEMI_MAJOR,
        semi_minor=g9.SPHEROID_SEMI_MINOR, e0=1.0, interior="continue",
    )
    cf = np.sqrt(np.sum(np.abs(np.asarray(cf)) ** 2, axis=1))
    pec_n = prof["pec"] / prof["pec"][-1]
    nat_n = prof["natural"] / prof["natural"][-1]
    cf_n = cf / cf[-1]
    say("\n[printed, no band] |E| along (0.015, 0, z) above the top rim, normalised at the far point "
        "(shape comparator only):")
    say("   z (mm)   PEC/far   natural/far   closed-form/far   PEC/cf")
    for i in range(N_PROFILE):
        say(f"  {1e3 * pts[i, 2]:7.3f}  {pec_n[i]:8.4f}  {nat_n[i]:11.4f}  {cf_n[i]:15.4f}  "
            f"{pec_n[i] / cf_n[i]:7.4f}")
    r0 = float(pec_n[0] / cf_n[0])
    say(f"  first point: PEC/far {pec_n[0]:.4f} vs natural/far {nat_n[0]:.4f} (predicted PEC above "
        f"natural: {'IN' if pec_n[0] > nat_n[0] else 'OUT'}); PEC/closed-form {r0:.4f} (predicted "
        f"[{PRED_PROFILE_VS_CF[0]:g}, {PRED_PROFILE_VS_CF[1]:g}], {g9._in(r0, PRED_PROFILE_VS_CF)})")
    OUTPUT_DIR.mkdir(exist_ok=True)
    csv_path = OUTPUT_DIR / f"{BASENAME}_rim_profile.csv"
    if rank0:
        np.savetxt(
            csv_path,
            np.column_stack([pts[:, 2], prof["pec"], prof["natural"], cf_n]),
            delimiter=",",
            header="z_m,E_pec_P1,E_natural_P1,closed_form_shape_normalised",
            comments="",
        )
    comm.Barrier()

    # ---- ParaView -----------------------------------------------------------------
    t0 = time.perf_counter()
    dg0 = m_pec["dg0"]

    def _fn(name, arr):
        f = fem.Function(dg0, name=name)
        f.x.array[:] = arr
        f.x.scatter_forward()
        return f

    ratio_full = np.zeros_like(m_pec["lsar"])
    pos = m_nat["loss"] > 0.0
    ratio_full[pos] = m_pec["lsar"][pos] / m_nat["lsar"][pos]
    fields = {
        "E_magnitude_pec_P1": _fn("E_magnitude_pec_P1", m_pec["emag"]),
        "local_sar_pec_P1": _fn("local_sar_pec_P1", m_pec["lsar"]),
        "E_magnitude_natural_P1": _fn("E_magnitude_natural_P1", m_nat["emag"]),
        "local_sar_natural_P1": _fn("local_sar_natural_P1", m_nat["lsar"]),
        "local_sar_pec_over_natural_P1": _fn("local_sar_pec_over_natural_P1", ratio_full),
    }
    for name, rim in (("one_gram_ball_top", "top"), ("one_gram_ball_bottom", "bottom")):
        b = _ball_mask(msh, dg0, RIM_POINTS[rim], a_1g)
        b.name = name
        fields[name] = b
    xdmf_path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_combined", msh, cell_tags, fields,
        comm=comm, facet_tags=disc["built"][2],
    )
    adopt_host_ownership(OUTPUT_DIR, comm=comm)
    timers["XDMF"] = time.perf_counter() - t0

    # ---- setup figure (no-op unless FEM_EM_SETUP_FIGURES=1), after the timers ------
    t0 = time.perf_counter()
    port_names = {100 + i: f"port P{i} box (lower)" for i in range(1, 5)}
    port_names.update({200 + i: f"port P{i} box (upper)" for i in range(1, 5)})
    write_setup_figure(
        msh,
        cell_tags,
        FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names={1: "conductor", 2: "air", 3: "phantom (saline)", **port_names},
        hide_tags=(2,),
        translucent_tags=(3,),
        slice_normal=(0.0, 1.0, 0.0),
        slice_origin=(0.0, 0.0, 0.0),
        title="imp:9 -- PEC oblate disc (a = 6 mm, aspect 5) in the driven birdcage phantom "
              "(IMP-3 step 9, W1 mesh)",
        window_size=(1400, 630),
        comm=comm,
    )
    timers["setup figure"] = time.perf_counter() - t0

    say("\n[times] " + ", ".join(f"{k} {v:.1f} s" for k, v in timers.items())
        + f"; summed peak RSS {g9._peak_rss_gib():.2f} GiB at -n {comm.size}")
    say(f"[paraview] wrote {xdmf_path} and {csv_path}")
    say(f"\nAll assertions hold. Total elapsed {time.perf_counter() - started:.1f} s.")


if __name__ == "__main__":
    main()
