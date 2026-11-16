"""`IMP-6` step 2e (§9 item 466, §10 chain I29): the body grain — the last
step-2 sub-step.  No solve.

For each of I28's three rim rungs (rho_r/3, /4, /5 on ``oblate_rim="edge"``)
the meshed body (a hole) is measured by its **second moments**: the spheroid
with the same per-volume moments has a_I = sqrt(5 (I_xx + I_yy) / (2 V_h)) and
c_I = sqrt(5 I_zz / V_h) (exact on a true spheroid, since ∫x² dV = V a²/5).
The moments are computed twice — box closed form minus the exact P1-tet
moments of the air cells (volume path), and as divergence-theorem surface
integrals over the cavity facets (surface path) — and the two must agree
(Mv).  The node extents a_v, c_v are printed as the control (every wall node
lies on the CAD surface on ``edge``, so they read a and c by construction).
The rim-ball comparand is re-evaluated on the nominal mask cells at the
measured axes; I28's solved Φ_h (s = 0.1) gives the body-corrected miss r^I.

Anchors: (Ce) cells, (V0) comparand reproduction, (Vr) the reading path,
(V1) the deficit, (Mv) two-path divergence identity, (B) |r^I| ≤ MAT4_BAND,
non-finite count 0.  Everything else (the E3 reading rule, the predicted
ranges) is printed, never asserted.  The step-2 family freezes after this
item whatever it reads (10-14 weekly).
"""

from __future__ import annotations

import time

import basix
import numpy as np
from mpi4py import MPI

from fem_em_solver.utils.analytical import (
    oblate_spheroid_depolarisation,
    pec_oblate_spheroid_field,
)

from tests.complex_mode import complex_only
from tests.validation.test_imp1_pec_sphere_lossy_host import CAVITY_TAG, E0, MAT4_BAND
from tests.validation.test_imp6_oblate_fem import (  # private names only (never test_*)
    BOX_HALF_EXTENTS,
    SEMI_MAJOR,
    _body_volume_deficit,
    _ncells,
    _numpy_quadrature_mean,
    _owned_masked_cells,
    _rho_r,
    _rim_mask,
    _semi_minor,
)
from tests.validation.test_imp6_oblate_rim_ladder import (  # private names only
    ASPECT,
    DIVISORS,
    ITEM426_RUNGS,
    ITEM440_R1,
    ITEM450_REPRO_BAR,
    _mesh_rung,
    _mesh_seconds,
)

# I28's solved Phi_h at s = 0.1 and its printed miss on these exact meshes:
# rho_r/3: 20261014T183329Z_IMP-6-step2d-W1a.log:173;
# rho_r/4: 20261014T183202Z_IMP-6-step2d-W1b.log:181;
# rho_r/5: 20261014T183415Z_IMP-6-step2d-W1c.log:194.
I28_PHI_H = {3.0: 22.781980, 4.0: 22.852453, 5.0: 22.888993}
I28_MISS_STR = {3.0: "+0.3900%", 4.0: "+0.6506%", 5.0: "+0.7762%"}
# I28's census (20261014T183126Z_IMP-6-step2d-W0r.log:173, :292, :424).
I28_DEFICIT_STR = {3.0: "-5.147e-03", 4.0: "-5.148e-03", 5.0: "-5.088e-03"}
I28_CELLS = {3.0: 33_334, 4.0: 53_959, 5.0: 89_530}
# I28's (Wk') / (Wk) on these meshes (…W1a.log:179, …W1b.log:188, …W1c.log:199), printed beside.
I28_WK = {3.0: "1.325e-8", 4.0: "1.199e-14", 5.0: "6.947e-10"}
I28_R_INF = "+0.893 %"
I28_DELTA = {"34": 0.00261, "45": 0.00126}
# Instrument tolerances (Mv), set >> round-off; box cancellation ~1 200x on I_zz.
MV_VOL_TOL = 1e-9
MV_MOM_TOL = 1e-8
V1_TOL = 1e-12
V0_R4_TOL = 1e-12
E3_QUARTER = 0.25  # item 450's (Fix) quarter; predicted threshold, printed only


def _volume_moments(msh):
    """(V, first moments (3,), second moments diag (3,)) of the hole: box minus air tets; reduced."""
    n_owned = msh.topology.index_map(msh.topology.dim).size_local
    v = msh.geometry.x[msh.geometry.dofmap[np.arange(n_owned)]]  # (n, 4, 3)
    jac = np.stack([v[:, 1] - v[:, 0], v[:, 2] - v[:, 0], v[:, 3] - v[:, 0]], axis=2)
    vol_t = np.abs(np.linalg.det(jac)) / 6.0
    s = np.sum(v, axis=1)  # (n, 3)
    first_t = vol_t[:, None] * s / 4.0
    second_t = (vol_t[:, None] / 20.0) * (np.sum(v * v, axis=1) + s * s)
    loc = np.concatenate([[np.sum(vol_t)], np.sum(first_t, axis=0), np.sum(second_t, axis=0)])
    air = np.zeros(7)
    msh.comm.Allreduce(loc, air, op=MPI.SUM)
    bx, by, bz = BOX_HALF_EXTENTS
    box_v = 8.0 * bx * by * bz
    box_second = box_v * np.array([bx * bx, by * by, bz * bz]) / 3.0
    return box_v - air[0], -air[1:4], box_second - air[4:7]


def _surface_moments(msh, facet_tags):
    """(V_s, I_ii,s (3,)) by the divergence theorem over owned cavity facets, normals
    oriented into the single adjacent air cell (outward from the hole); reduced."""
    tdim = msh.topology.dim
    msh.topology.create_connectivity(tdim - 1, tdim)
    msh.topology.create_connectivity(tdim - 1, 0)
    f2c = msh.topology.connectivity(tdim - 1, tdim)
    n_owned_f = msh.topology.index_map(tdim - 1).size_local
    facets = facet_tags.find(CAVITY_TAG)
    facets = facets[facets < n_owned_f]
    from dolfinx.mesh import entities_to_geometry

    gx = msh.geometry.x
    if len(facets):
        fgeo = entities_to_geometry(msh, tdim - 1, facets, False)
        fv = gx[fgeo]  # (nf, 3, 3)
        adj = np.array([f2c.links(f)[0] for f in facets], dtype=np.int32)
        nadj = np.array([len(f2c.links(f)) for f in facets])
        cv = gx[msh.geometry.dofmap[adj]]  # (nf, 4, 3)
        ccent = np.mean(cv, axis=1)
        e1 = fv[:, 1] - fv[:, 0]
        e2 = fv[:, 2] - fv[:, 0]
        cr = np.cross(e1, e2)  # |cr| = 2 area
        fcent = np.mean(fv, axis=1)
        sign = np.sign(np.sum(cr * (ccent - fcent), axis=1))
        cr = cr * sign[:, None]  # area-weighted normal (x2), into the air cell
        pts, wts = basix.make_quadrature(basix.CellType.triangle, 3)
        phys = fv[:, 0][:, None, :] + pts[None, :, 0:1] * e1[:, None, :] + pts[None, :, 1:2] * e2[:, None, :]
        # ∫ f n_i dS = Σ_q w_q f(x_q) cr_i  (reference area 1/2 ↔ |cr|/2)
        vs = float(np.sum(np.sum(wts[None, :] * phys[:, :, 0], axis=1) * cr[:, 0]))
        ii = np.array([
            float(np.sum(np.sum(wts[None, :] * phys[:, :, i] ** 3 / 3.0, axis=1) * cr[:, i]))
            for i in range(3)
        ])
        bad_adj = int(np.sum(nadj != 1))
        nzero = int(np.sum(sign == 0))
        # wall nodes: extents
        rho = np.hypot(fv[..., 0], fv[..., 1])
        a_v_loc, c_v_loc = float(np.max(rho)), float(np.max(np.abs(fv[..., 2])))
    else:
        vs, ii, bad_adj, nzero, a_v_loc, c_v_loc = 0.0, np.zeros(3), 0, 0, 0.0, 0.0
    loc = np.concatenate([[vs], ii, [float(len(facets)), float(bad_adj), float(nzero)]])
    tot = np.zeros(7)
    msh.comm.Allreduce(loc, tot, op=MPI.SUM)
    a_v = msh.comm.allreduce(a_v_loc, op=MPI.MAX)
    c_v = msh.comm.allreduce(c_v_loc, op=MPI.MAX)
    return tot[0], tot[1:4], int(tot[4]), int(tot[5]), int(tot[6]), a_v, c_v


def _phi_cf(msh, cells, a_, c_):
    """Rim-ball comparand at axes (a_, c_) on ``cells``, plus the reduced non-finite |E|² count."""
    bad = [0]

    def fn(p):
        e = pec_oblate_spheroid_field(p, semi_major=a_, semi_minor=c_, drive="inplane",
                                      e0=E0, interior="continue")
        bad[0] += int(np.sum(~np.isfinite(np.sum(np.abs(e) ** 2, axis=1))))
        return e

    phi = _numpy_quadrature_mean(msh, cells, fn)
    nbad = msh.comm.allreduce(bad[0], op=MPI.SUM)
    return phi, nbad


def _lx(a_, c_):
    e = np.sqrt(1.0 - (c_ / a_) ** 2)
    return oblate_spheroid_depolarisation(e)


@complex_only
def test_body_grain():
    comm = MPI.COMM_WORLD
    a = SEMI_MAJOR
    c = _semi_minor(ASPECT)
    v_nom = 4.0 * np.pi * a * a * c / 3.0
    lx_nom, lz_nom = _lx(a, c)
    out = {}
    fails = []
    for d in DIVISORS:
        t0 = time.perf_counter()
        msh, cell_tags, facet_tags = _mesh_rung(d)
        cells_n = _ncells(msh)
        deficit = _body_volume_deficit(msh, ASPECT)
        vh, m1, ih = _volume_moments(msh)
        vs, is_, nfac, bad_adj, nzero, a_v, c_v = _surface_moments(msh, facet_tags)
        a_i = float(np.sqrt(5.0 * (ih[0] + ih[1]) / (2.0 * vh)))
        c_i = float(np.sqrt(5.0 * ih[2] / vh))
        cent = m1 / vh
        lx_i, lz_i = _lx(a_i, c_i)
        chi = _rim_mask(msh, ASPECT)
        cells = _owned_masked_cells(msh, chi)
        n_mask = msh.comm.allreduce(len(cells), op=MPI.SUM)
        phi_nom, bad_nom = _phi_cf(msh, cells, a, c)
        phi_i, bad_i = _phi_cf(msh, cells, a_i, c_i)
        r = I28_PHI_H[d] / phi_nom - 1.0
        r_i = I28_PHI_H[d] / phi_i - 1.0
        dt = time.perf_counter() - t0
        v0 = abs(phi_nom / ITEM426_RUNGS[d]["phi_cf_h"] - 1.0)
        v0_r4 = abs(phi_nom / ITEM440_R1["phi_cf_h8"] - 1.0) if d == 4.0 else None
        v1 = abs(vh / v_nom - 1.0 - deficit)
        mv_v = abs(vs / vh - 1.0)
        mv_i = [abs(is_[i] / ih[i] - 1.0) for i in range(3)]
        out[d] = dict(r=r, r_i=r_i, a_i=a_i, c_i=c_i)
        if comm.rank == 0:
            print(f"\n[item 466] rung rho_r/{d:g} (rho_r = {_rho_r(ASPECT) * 1e3:.4f} mm): "
                  f"cells {cells_n} (record {I28_CELLS[d]}), mesh {_mesh_seconds(d):.1f} s, "
                  f"rung {dt:.1f} s, mask cells {n_mask}, cavity facets {nfac} "
                  f"(adjacency != 1: {bad_adj}, zero-sign: {nzero})")
            print(f"  (V1) deficit {deficit:.3e} (record {I28_DEFICIT_STR[d]}); "
                  f"|V_h/V - 1 - deficit| = {v1:.3e}")
            print(f"  (Mv) V_h = {vh:.12e}, V_s = {vs:.12e}, |V_s/V_h - 1| = {mv_v:.3e}; "
                  f"I_ii,h = {ih[0]:.9e} / {ih[1]:.9e} / {ih[2]:.9e}; "
                  f"|I_s/I_h - 1| = {mv_i[0]:.3e} / {mv_i[1]:.3e} / {mv_i[2]:.3e}")
            print(f"  centroid (m) = ({cent[0]:.3e}, {cent[1]:.3e}, {cent[2]:.3e})")
            print(f"  (N-ext, control) a_v/a - 1 = {a_v / a - 1:+.3e}, c_v/c - 1 = {c_v / c - 1:+.3e} "
                  f"(predicted <= 1e-9: extents carry nothing)")
            print(f"  moment body: a_I/a - 1 = {a_i / a - 1:+.4e}, c_I/c - 1 = {c_i / c - 1:+.4e} "
                  f"(predicted |.| in [1e-4, 5e-3], c more negative); c_I/a_I = {c_i / a_i:.6f} "
                  f"(nominal {c / a:.6f})")
            print(f"  L_x: {lx_nom:.8f} -> {lx_i:.8f} (1/L_x {1 / lx_nom:.6f} -> {1 / lx_i:.6f}); "
                  f"L_z {lz_nom:.8f} -> {lz_i:.8f}; (L_x/L_x^I)^2 - 1 = {(lx_nom / lx_i) ** 2 - 1:+.4%}")
            print(f"  (V0) Phi_cf(a, c) = {phi_nom:.12f} vs item 426 {ITEM426_RUNGS[d]['phi_cf_h']:.6f}: "
                  f"|ratio - 1| = {v0:.3e} (bar {ITEM450_REPRO_BAR:g})"
                  + (f"; vs item 440 {ITEM440_R1['phi_cf_h8']:.12f}: {v0_r4:.3e} (bar 1e-12)"
                     if v0_r4 is not None else "")
                  + f"; I28 (Wk')/(Wk) {I28_WK[d]}")
            print(f"  Phi_cf(a_I, c_I) = {phi_i:.12f}; Phi_cf(a_I,c_I)/Phi_cf(a,c) - 1 = "
                  f"{phi_i / phi_nom - 1:+.4%} (predicted |.| in [0.05 %, 3 %], sign not predicted)")
            print(f"  non-finite |E|^2: nominal {bad_nom}, measured {bad_i}")
            print(f"  (Vr) r_d = {r:+.4%} (I28 {I28_MISS_STR[d]}); body-corrected r_d^I = {r_i:+.4%}"
                  + ("  [no correction applied]" if r_i == r else ""))
        # asserts outside the rank-0 branch (all values reduced / rank-identical)
        checks = [
            ("Ce", cells_n == I28_CELLS[d]),
            ("V0", v0 <= ITEM450_REPRO_BAR),
            ("V0-r4", v0_r4 is None or v0_r4 <= V0_R4_TOL),
            ("Vr", f"{r:+.4%}" == I28_MISS_STR[d]),
            ("V1-str", f"{deficit:.3e}" == I28_DEFICIT_STR[d]),
            ("V1", v1 <= V1_TOL),
            ("Mv-V", mv_v <= MV_VOL_TOL),
            ("Mv-I", max(mv_i) <= MV_MOM_TOL),
            ("B", abs(r_i) <= MAT4_BAND),
            ("finite", bad_nom == 0 and bad_i == 0),
        ]
        fails += [f"rho_r/{d:g}:{name}" for name, ok in checks if not ok]

    d3, d4, d5 = DIVISORS
    del34 = out[d4]["r_i"] - out[d3]["r_i"]
    del45 = out[d5]["r_i"] - out[d4]["r_i"]
    q = del45 / del34 if del34 != 0 else float("nan")
    r_inf = out[d5]["r_i"] + del45 * q / (1 - q) if 0 < q < 1 else float("nan")
    level = (all(abs(out[d]["r_i"]) <= abs(out[d]["r"]) for d in DIVISORS)
             and abs(out[d5]["r_i"]) <= E3_QUARTER * abs(out[d5]["r"]))
    partial = all(abs(out[d]["r_i"]) < abs(out[d]["r"]) for d in DIVISORS)
    verdict = "E3-level" if level else ("E3-partial" if partial else "E3-not")
    if comm.rank == 0:
        print("\n[item 466] summary (printed, never asserted):")
        print("  r_d (I28)  : " + " / ".join(f"{out[d]['r']:+.4%}" for d in DIVISORS))
        print("  r_d^I      : " + " / ".join(f"{out[d]['r_i']:+.4%}" for d in DIVISORS))
        print("  a_I/a - 1  : " + " / ".join(f"{out[d]['a_i'] / a - 1:+.4e}" for d in DIVISORS))
        print("  c_I/c - 1  : " + " / ".join(f"{out[d]['c_i'] / c - 1:+.4e}" for d in DIVISORS))
        print(f"  Delta^I_34 = {del34 * 100:+.4f} pp (predicted within +-0.05 pp of "
              f"{I28_DELTA['34'] * 100:+.3f}); Delta^I_45 = {del45 * 100:+.4f} pp (predicted within "
              f"+-0.05 pp of {I28_DELTA['45'] * 100:+.3f}); q^I = {q:.4f}; r_inf^I = "
              f"{r_inf:+.4%} (I28 {I28_R_INF})")
        print(f"  reading rule: {verdict} (predicted: not E3-level). "
              f"The step-2 family is frozen (10-14 weekly).")
        print(f"  failures: {fails}")
    assert not fails, f"item 466 anchors red: {fails}"
