"""Example (`EX-112`, `mat:19`): the F2182-sized gel box -- a 90 mm slab to ≈ 3 % at its
centre, cooled ≈ 80 % through its two large faces, and hours from steady.

`MAT-21` step 1 gates a P1 FEM of the 650 x 420 x 90 mm gel box (uniform q, Robin
h = 10 on every face) against step 0's triple Robin eigen-series.  Here the gate's own
m = 6 rung is re-run through its module-level helpers (unchanged; its pytest fixtures
cannot be called) and three readings no example has are printed beside it: the 1-D slab
estimate against the box along the three centre lines, the heat leaving through each
face pair against the series' face integrals, and the centre history against q t / rho c.

Everything -- every band, every helper -- is imported from
``tests/validation/test_mat21_box_fem.py`` (``m1``) and
``tests/validation/test_mat21_box_closed_form.py`` (``m0``).  Records live only in the
example-local ``REC`` dict, each entry commented with its ``log:line``; they are printed
beside the readings and never asserted by equality.

Asserted (imported bands only): (S1) interior, (S2) Gauss, (T) interior, (E-cn) balance;
(F0) six face tags partition the Robin facets (integer counts); (F1) sum of the six face
powers = q V_h within ``S2_GAUSS_BAND``.  Controls (asserted as misses): (c1) comparand
at 2h, (c2) Dirichlet z-roots at t*, (F1-c) the partition with the -z face dropped.
Not claimed: F2182 compliance, a coil or SAR source, the duty cycle, the LSQ rate,
degree 2, an external check.

Run through the example runner (the ``mat:`` group sources the complex build)::

    ./scripts/run_examples.sh -e mat:19 -n 2 -t 420

Output: ``examples/materials/paraview_output/materials_19_f2182_gel_box_slab_and_faces_combined.xdmf``
and ``..._centre_lines.csv``.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import ufl
from dolfinx import default_scalar_type, fem
from dolfinx import mesh as dmesh
from mpi4py import MPI

from fem_em_solver.io.paraview_utils import adopt_host_ownership, write_xdmf_with_tags
from fem_em_solver.post.evaluation import evaluate_vector_field_parallel
from fem_em_solver.post.setup_figure import write_setup_figure
from fem_em_solver.thermal import solve_transient_conduction

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

import tests.validation.test_mat12_transient_sar_heating_sphere as m12  # noqa: E402
import tests.validation.test_mat13_transient_perfusion_sphere as m13  # noqa: E402
import tests.validation.test_mat15_finite_cylinder_fem as m15b  # noqa: E402
import tests.validation.test_mat21_box_closed_form as m0  # noqa: E402
import tests.validation.test_mat21_box_fem as m1  # noqa: E402

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"  # committed, EX-57
BASENAME = "materials_19_f2182_gel_box_slab_and_faces"

W1 = "20261013T142646Z_MAT-21-step1r-W1-steady.log"
W2 = "20261013T142900Z_MAT-21-step1r-W2-transient.log"
W0 = "20261013T004329Z_MAT-21-step0-W1.log"

# Records, printed beside readings only (EX-89 item 332 (Q1)); never asserted by equality.
REC = {
    "cells": 269568,          # W1 :106
    "whole": 1.547121e-03,    # W1 :106
    "gauss": 1.682e-12,       # W1 :106 (round-off floor; band assert only)
    "s1": 1.417055e-03,       # W1 :108
    "th0": 12.005487,         # W1 :108
    "th0_series": 12.003947,  # W1 :108
    "c1": 5.9734e-01,         # W1 :110
    "t_int": 9.071528e-04,    # W2 :92
    "t_whole": 1.417652e-03,  # W2 :92
    "balance": 2.809e-12,     # W2 :92 (round-off floor; band assert only)
    "c2": 3.6135e-01,         # W2 :93
    "gam1": (1.330018, 1.232272, 0.771359),  # W0 :91
    "tau1": 20193.0,          # W0 :94
    "t_C": 14107.5,           # W0 :94
    "frac_tstar": 0.1113,     # W0 :80 (1.335697 / 12.003947)
    "frac_tC": 0.4592,        # W0 :80 (5.511710 / 12.003947)
}
# Predictions from stdlib arithmetic while drafting (brief item 438), printed only.
PRED_FACE = {"x": 0.078249, "y": 0.126245, "z": 0.795501}
PRED_MEAN = {"x": 5.0862, "y": 5.3023, "z": 7.1595}
PRED_F1C = 0.39775

A, B, C = m1.A_BOX, m1.B_BOX, m1.C_BOX
LHALF = (A, B, C)
Q, K, H, T_INF = m1.Q, m1.K_GEL, m1.H_GEL, m1.T_INF
T_STAR, N_MAIN = m1.T_STAR, m1.N_MAIN
# (axis, sign) of the six example-local face tags: 1/2 = ±x, 3/4 = ±y, 5/6 = ±z.
FACES = {1: (0, +1), 2: (0, -1), 3: (1, +1), 4: (1, -1), 5: (2, +1), 6: (2, -1)}
AXNAME = "xyz"


def _slab(z):
    """1-D Robin slab limit (step 0's (R5) formula): q (C^2 - z^2) / (2 k) + q C / h."""
    return Q * (C**2 - z**2) / (2.0 * K) + Q * C / H


def _say(msg):
    m1._say(msg)


def _at(T, pts, offset):
    """T - offset at points, real; raises if any point is not found on any rank."""
    vals, ok = evaluate_vector_field_parallel(T, pts)
    out = np.empty(pts.shape[0])
    for i in range(pts.shape[0]):
        if not bool(ok[i]):
            raise RuntimeError(f"point {pts[i]} not found on any rank")
        out[i] = float(np.real(np.ravel(vals[i])[0])) - offset
    return out


def _series_at_points(fn, x, y, z, chunk=20000):
    out = np.empty(x.shape[0])
    for s in range(0, x.shape[0], chunk):
        sl = slice(s, s + chunk)
        out[sl] = np.real(fn(x[sl], y[sl], z[sl]))
    return out


def _face_quadrature():
    """Series face integrals h * int theta dA per face (128-pt Gauss-Legendre, (R3) locally)."""
    xg, wg = np.polynomial.legendre.leggauss(128)
    fn = m1._steady_ex()
    out = {}
    for d in range(3):
        o = [i for i in range(3) if i != d]
        u, v = np.meshgrid(xg, xg, indexing="ij")
        wuv = np.outer(wg, wg).ravel() * LHALF[o[0]] * LHALF[o[1]]
        for sgn in (1.0, -1.0):
            pts = [None, None, None]
            pts[d] = np.full(u.size, sgn * LHALF[d])
            pts[o[0]] = u.ravel() * LHALF[o[0]]
            pts[o[1]] = v.ravel() * LHALF[o[1]]
            th = np.real(fn(*pts))
            out[(d, sgn)] = (H * float(np.sum(wuv * th)), float(np.sum(wuv)))
    return out


def main() -> None:
    comm = MPI.COMM_WORLD
    rank0 = comm.rank == 0

    def say(msg=""):
        if rank0:
            print(msg, flush=True)

    if not np.issubdtype(np.dtype(default_scalar_type), np.complexfloating):
        raise RuntimeError(
            "this example needs the complex DolfinX build "
            "(source /usr/local/bin/dolfinx-complex-mode); the runner's `mat:` "
            "group sources it automatically"
        )

    structural = {
        "m1.INTERIOR_BAND is m13.INTERIOR_BAND": m1.INTERIOR_BAND is m13.INTERIOR_BAND,
        "m1.E_BAND is m12.E_BAND": m1.E_BAND is m12.E_BAND,
        "m1.S2_GAUSS_BAND is m15b.S2_GAUSS_BAND": m1.S2_GAUSS_BAND is m15b.S2_GAUSS_BAND,
        "m1.HALF is m0.HALF": m1.HALF is m0.HALF,
        "m1.RUNGS[-1] == 6": m1.RUNGS[-1] == 6,
        "m1.T_STAR == 0.2 * m1.T_C": m1.T_STAR == 0.2 * m1.T_C,
        "m1.K_GEL == m1.K and m1.H_GEL == m1.H": m1.K_GEL == m1.K and m1.H_GEL == m1.H,
    }
    for name, ok in structural.items():
        assert ok, f"structural import check failed: {name}"
    say("=" * 72)
    say("EX-112 -- the F2182-sized gel box: slab at the centre, cooled through its large faces")
    say("=" * 72)
    say("[structural] " + "; ".join(f"{k}: {v}" for k, v in structural.items()))
    say(f"[records] steady beside {W1}:106-110; transient beside {W2}:92-93; closed form beside {W0}")
    bands = (f"INTERIOR_BAND {m1.INTERIOR_BAND:.0e}, S2_GAUSS_BAND {m1.S2_GAUSS_BAND:.0e}, "
             f"E_BAND {m1.E_BAND:.0e}")
    say(f"[bands] imported: {bands}")

    clk = m1._Clock("EX-112")

    # ---- (1) anchor: steady, the gate's own helpers --------------------------------
    msh, ft, ns, n, v_h, st = m1._rung(6, clk)
    tdim = msh.topology.dim
    fdim = tdim - 1
    imap_c = msh.topology.index_map(tdim)
    n_cells_all = imap_c.size_local + imap_c.num_ghosts

    cells = np.arange(n_cells_all, dtype=np.int32)
    cell_tags = dmesh.meshtags(msh, tdim, cells, np.full(n_cells_all, 1, dtype=np.int32))
    write_setup_figure(
        msh,
        cell_tags,
        FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names={1: "gel (650 x 420 x 90 mm)"},
        title="mat:19 -- F2182-sized gel box (650 x 420 x 90 mm), Robin on all six faces, m = 6 rung: x-z slice",
        slice_normal=(0.0, 1.0, 0.0),
        slice_origin=(0.0, 0.0, 0.0),
        clip_normal=(0.0, -1.0, 0.0),
        comm=comm,
    )
    clk.lap("setup figure (no-op unless FEM_EM_SETUP_FIGURES=1)")

    th = st.temperature - T_INF
    ex, npts = m1._exact(msh, m1._steady_ex())
    clk.lap(f"steady series interpolation ({npts} P2 points)")
    chi = m1._interior_chi(msh)
    e_whole = m1._rel(msh, th, ex)
    s1 = m1._rel(msh, th, ex, chi)
    ds_r = ufl.Measure("ds", domain=msh, subdomain_data=ft)(m1.ROBIN_TAG)
    conv = float(np.real(m1._allsum(comm, H * th * ds_r)))
    gauss = abs(conv - Q * v_h) / (Q * v_h)
    th_c = float(_at(st.temperature, np.zeros((1, 3)), T_INF)[0])
    th_c_ex = float(np.real(m1._steady_ex()(np.zeros(1), np.zeros(1), np.zeros(1))[0]))
    clk.lap("steady anchors")

    ex_c1, _ = m1._exact(msh, m1._steady_ex(h=2.0 * m1.H))
    c1 = m1._rel(msh, th, ex_c1, chi)
    del ex_c1
    clk.lap("control (c1)")

    say()
    say(f"[anchor] cells {n} (record {REC['cells']}, W1:106); whole L2 {e_whole:.6e} "
        f"(record {REC['whole']:.6e}, W1:106)")
    say(f"[anchor] (S1) interior {s1:.6e} <= {m1.INTERIOR_BAND:.0e} (record {REC['s1']:.6e}, W1:108)")
    say(f"[anchor] (S2) Gauss {gauss:.3e} <= {m1.S2_GAUSS_BAND:.0e} (record {REC['gauss']:.3e}, W1:106; "
        "round-off floor, band only)")
    say(f"[anchor] theta_h(0) {th_c:.6f} K beside series {th_c_ex:.6f} K (records {REC['th0']:.6f} / "
        f"{REC['th0_series']:.6f}, W1:108)")
    say(f"[control] (c1) comparand at 2h: interior miss {c1:.4e} > {m1.INTERIOR_BAND:.0e} "
        f"(record {REC['c1']:.4e}, W1:110)")
    assert s1 <= m1.INTERIOR_BAND, f"(S1) interior {s1:.4e}"
    assert gauss <= m1.S2_GAUSS_BAND, f"(S2) Gauss {gauss:.3e}"
    assert c1 > m1.INTERIOR_BAND, f"(c1) did not miss: {c1:.4e}"

    # ---- (F0) face tags partition the Robin boundary --------------------------------
    msh.topology.create_connectivity(fdim, tdim)
    size_f = msh.topology.index_map(fdim).size_local
    idx_all, tag_all = [], []
    for tag, (d, sgn) in FACES.items():
        f = dmesh.locate_entities_boundary(
            msh, fdim,
            lambda x, d=d, sgn=sgn: np.isclose(x[d], sgn * LHALF[d], rtol=0, atol=1e-9))
        idx_all.append(np.asarray(f, dtype=np.int32))
        tag_all.append(np.full(f.size, tag, dtype=np.int32))
    idx = np.concatenate(idx_all) if idx_all else np.zeros(0, dtype=np.int32)
    tg = np.concatenate(tag_all) if tag_all else np.zeros(0, dtype=np.int32)
    order = np.argsort(idx, kind="stable")
    face_tags = dmesh.meshtags(msh, fdim, idx[order], tg[order])
    face_tags.name = "face_tags"

    ndiv = tuple(ns)  # (78, 48, 12); 2 * N_a * N_b triangles per face
    expected = {1: ndiv[1] * ndiv[2] * 2, 3: ndiv[0] * ndiv[2] * 2, 5: ndiv[0] * ndiv[1] * 2}
    expected_counts = {1: expected[1], 2: expected[1], 3: expected[3], 4: expected[3],
                       5: expected[5], 6: expected[5]}
    counts = {}
    for tag in FACES:
        loc = face_tags.find(tag)
        counts[tag] = comm.allreduce(int(np.count_nonzero(loc < size_f)), op=MPI.SUM)
    robin_loc = ft.find(m1.ROBIN_TAG)
    n_robin = comm.allreduce(int(np.count_nonzero(robin_loc < size_f)), op=MPI.SUM)
    one = fem.Constant(msh, default_scalar_type(1.0))
    ds = ufl.Measure("ds", domain=msh, subdomain_data=face_tags)
    areas = {tag: float(np.real(m1._allsum(comm, one * ds(tag)))) for tag in FACES}
    exact_area = {1: 4 * B * C, 2: 4 * B * C, 3: 4 * A * C, 4: 4 * A * C, 5: 4 * A * B, 6: 4 * A * B}
    say()
    for tag, (d, sgn) in FACES.items():
        say(f"[F0] tag {tag} ({'+' if sgn > 0 else '-'}{AXNAME[d]}): owned facets {counts[tag]} "
            f"(expected {expected_counts[tag]}); area {areas[tag]:.9e} m^2 (exact {exact_area[tag]:.9e}, "
            f"rel {areas[tag] / exact_area[tag] - 1.0:+.2e})")
    say(f"[F0] sum of tagged facets {sum(counts.values())} vs Robin facets {n_robin} (expected 21024)")
    assert counts == expected_counts, f"(F0) per-tag counts {counts} vs {expected_counts}"
    assert sum(counts.values()) == 21024 == n_robin, (counts, n_robin)

    # ---- (F1) the six face powers sum to the gate's Gauss ---------------------------
    powers = {tag: float(np.real(m1._allsum(comm, H * th * ds(tag)))) for tag in FACES}
    qv = Q * v_h
    sum_p = sum(powers.values())
    f1 = abs(sum_p / qv - 1.0)
    sum_drop = sum(p for tag, p in powers.items() if tag != 6)
    f1c = abs(sum_drop / qv - 1.0)
    say(f"[F1] |sum P_f / (q V_h) - 1| = {f1:.3e} <= {m1.S2_GAUSS_BAND:.0e} (equal to the (S2) reading "
        f"{gauss:.3e} to ~1e-13 predicted)")
    say(f"[control] (F1-c) -z face dropped: |sum_(f!=6) P_f / (q V_h) - 1| = {f1c:.5f} > "
        f"{m1.S2_GAUSS_BAND:.0e} (predicted {PRED_F1C}, range [0.39, 0.41])")
    assert f1 <= m1.S2_GAUSS_BAND, f"(F1) {f1:.3e}"
    assert f1c > m1.S2_GAUSS_BAND, f"(F1-c) did not miss: {f1c:.3e}"

    # series face integrals and the pair table (printed)
    quad = _face_quadrature()
    clk.lap("series face quadrature")
    v_box = m1.V_BOX
    say()
    say("[faces] pair | FEM P/(qV_h) | series P/(qV) | rel | prediction | area share | FEM mean theta | "
        "series mean theta | prediction (K)")
    area_tot = sum(areas.values())
    for d, (t_pos, t_neg) in enumerate(((1, 2), (3, 4), (5, 6))):
        fem_p = powers[t_pos] + powers[t_neg]
        fem_a = areas[t_pos] + areas[t_neg]
        s_p = quad[(d, 1.0)][0] + quad[(d, -1.0)][0]
        s_a = quad[(d, 1.0)][1] + quad[(d, -1.0)][1]
        ax = AXNAME[d]
        say(f"  ±{ax}  {fem_p / qv:.6f} | {s_p / (Q * v_box):.6f} | {fem_p / qv / (s_p / (Q * v_box)) - 1.0:+.3e} "
            f"| {PRED_FACE[ax]:.6f} | {fem_a / area_tot:.4f} | {fem_p / (H * fem_a):.4f} | "
            f"{s_p / (H * s_a):.4f} | {PRED_MEAN[ax]:.4f}")
    s_all = sum(v[0] for v in quad.values())
    a_all = sum(v[1] for v in quad.values())
    say(f"[faces] series sum {s_all / (Q * v_box):.7f} (= prod S_d(40), step 0 deficit 4.6708e-06, {W0}:67); "
        f"all-face mean theta {s_all / (H * a_all):.4f} K (record 6.6531, {W0}:68)")

    # ---- slab against box: centre lines (printed) -----------------------------------
    fn_s = m1._steady_ex()
    say()
    say("[slab] centre lines: s/L | FEM theta | series theta | theta_slab | slab/series - 1")
    fracs = (0.0, 0.25, 0.5, 0.7, 0.9, 0.995)
    pred_slab = {"z": {0.0: 3.09, 0.5: 3.07, 0.9: 3.05, 0.995: 3.05},
                 "x": {0.25: 3.5, 0.5: 5.7, 0.7: 12.7, 0.9: 42.0, 0.995: 101.0},
                 "y": {0.5: 9.5, 0.7: 20.7, 0.9: 53.0, 0.995: 98.0}}
    for d in range(3):
        ax = AXNAME[d]
        pts = np.zeros((len(fracs), 3))
        pts[:, d] = np.array(fracs) * LHALF[d]
        t_f = _at(st.temperature, pts, T_INF)
        t_s = np.real(fn_s(pts[:, 0], pts[:, 1], pts[:, 2]))
        t_l = _slab(pts[:, 2])
        for i, fr in enumerate(fracs):
            p = pred_slab[ax].get(fr)
            pr = f"  (predicted ≈ {p} %)" if p is not None else ""
            say(f"  {ax} {fr:5.3f}  {t_f[i]:.6f} | {t_s[i]:.6f} | {t_l[i]:.6f} | "
                f"{100.0 * (t_l[i] / t_s[i] - 1.0):+.2f} %{pr}")
    # where the slab estimate leaves 5 % and 10 % along x and y (series only)
    for d, pr in ((0, 0.64), (1, 0.51)):
        s = np.linspace(0.0, 0.995, 1991) * LHALF[d]
        pts = np.zeros((s.size, 3))
        pts[:, d] = s
        r = _slab(pts[:, 2]) / np.real(fn_s(pts[:, 0], pts[:, 1], pts[:, 2])) - 1.0
        cross = {}
        for thr in (0.05, 0.10):
            j = int(np.argmax(r > thr)) if np.any(r > thr) else -1
            if j > 0:
                w = (thr - r[j - 1]) / (r[j] - r[j - 1])
                cross[thr] = (s[j - 1] + w * (s[j] - s[j - 1])) / LHALF[d]
            else:
                cross[thr] = float("nan")
        say(f"[slab] slab/series - 1 crosses 5 % at |{AXNAME[d]}|/{'AB'[d]} = {cross[0.05]:.3f}, 10 % at "
            f"{cross[0.10]:.3f} (10 % predicted ≈ {pr})")

    # the CSV lines: 101 points each
    sline = np.linspace(-0.995, 0.995, 101)
    csv_rows = []
    for d in range(3):
        pts = np.zeros((sline.size, 3))
        pts[:, d] = sline * LHALF[d]
        t_f = _at(st.temperature, pts, T_INF)
        t_s = np.real(fn_s(pts[:, 0], pts[:, 1], pts[:, 2]))
        t_l = _slab(pts[:, 2])
        for i in range(sline.size):
            csv_rows.append((AXNAME[d], sline[i], pts[i, 0], pts[i, 1], pts[i, 2], t_f[i], t_s[i],
                             t_l[i], t_l[i] / t_s[i] - 1.0))
    clk.lap("centre lines")

    # volume shares from the series at cell midpoints (numpy masks)
    mids = dmesh.compute_midpoints(msh, tdim, np.arange(n_cells_all, dtype=np.int32))
    ser_mid = _series_at_points(fn_s, mids[:, 0], mids[:, 1], mids[:, 2])
    ratio_mid = _slab(mids[:, 2]) / ser_mid - 1.0
    V0 = fem.functionspace(msh, ("DG", 0))
    shares = {}
    mask10 = None
    for thr, key in ((0.05, "5"), (0.10, "10")):
        mf = fem.Function(V0, name=f"slab_within_{key}pct")
        mf.x.array[:] = (np.abs(ratio_mid) <= thr).astype(float)
        mf.x.scatter_forward()
        shares[key] = float(np.real(m1._allsum(comm, mf * ufl.dx(domain=msh)))) / v_h
        if key == "10":
            mask10 = mf
    del ser_mid
    clk.lap("series at cell midpoints, volume shares")
    say(f"[slab] volume share with |slab/series - 1| <= 5 %: {shares['5']:.4f} (predicted [0.06, 0.16]); "
        f"<= 10 %: {shares['10']:.4f} (predicted [0.20, 0.35])")

    # ---- transient: the gate's (T) call, copied verbatim -----------------------------
    res = solve_transient_conduction(
        msh, ft, k=K, rho_c=m1.RHO_C, q=Q, h=H, t_inf=T_INF, t_init=T_INF,
        robin_tag=m1.ROBIN_TAG, dt=T_STAR / N_MAIN, n_steps=N_MAIN, theta=0.5,
        startup_half_steps=2, degree=1, snapshot_steps=range(0, N_MAIN + 1, 5))
    clk.lap(f"CN N={N_MAIN}")
    th_t = res.temperature - T_INF
    ex_t, npts_t = m1._exact(msh, m1._trans_ex(T_STAR))
    clk.lap(f"transient series interpolation at t* ({npts_t} P2 points)")
    t_int = m1._rel(msh, th_t, ex_t, chi)
    t_whole = m1._rel(msh, th_t, ex_t)
    bal = res.max_balance_residual
    ex_c2, _ = m1._exact(
        msh, m1._trans_ex(T_STAR, roots=(None, None, (np.arange(1, 41) - 0.5) * np.pi)))
    c2 = m1._rel(msh, th_t, ex_c2, chi)
    del ex_c2
    clk.lap("control (c2)")
    say()
    say(f"[anchor] (T) CN N={N_MAIN} at t* = {T_STAR:.1f} s: interior {t_int:.6e} <= {m1.INTERIOR_BAND:.0e} "
        f"(record {REC['t_int']:.6e}, W2:92); whole {t_whole:.6e} (record {REC['t_whole']:.6e}, W2:92)")
    say(f"[anchor] (E-cn) balance {bal:.3e} <= {m1.E_BAND:.0e} (record {REC['balance']:.3e}, W2:92; "
        "round-off floor, band only)")
    say(f"[control] (c2) Dirichlet z-roots at t*: interior miss {c2:.4e} > {m1.INTERIOR_BAND:.0e} "
        f"(record {REC['c2']:.4e}, W2:93)")
    assert t_int <= m1.INTERIOR_BAND, f"(T) interior {t_int:.4e}"
    assert bal <= m1.E_BAND, f"(E-cn) balance {bal:.3e}"
    assert c2 > m1.INTERIOR_BAND, f"(c2) did not miss: {c2:.4e}"

    # covered-record agreement flags (printed, not asserted)
    covered = [
        ("cells", f"{n}", f"{REC['cells']}"),
        ("whole L2", f"{e_whole:.6e}", f"{REC['whole']:.6e}"),
        ("(S1) interior", f"{s1:.6e}", f"{REC['s1']:.6e}"),
        ("theta_h(0)", f"{th_c:.6f}", f"{REC['th0']:.6f}"),
        ("(c1)", f"{c1:.4e}", f"{REC['c1']:.4e}"),
        ("(T) interior", f"{t_int:.6e}", f"{REC['t_int']:.6e}"),
        ("(T) whole", f"{t_whole:.6e}", f"{REC['t_whole']:.6e}"),
        ("(c2)", f"{c2:.4e}", f"{REC['c2']:.4e}"),
    ]
    say("[records] " + "; ".join(f"{nm}: {got} vs {rec} {'same' if got == rec else 'DIFFERS'}"
                                 for nm, got, rec in covered))

    # ---- hours from steady: centre history (printed) ---------------------------------
    say()
    say("[history] step | t (s) | FEM theta(0,t) | series theta(0,t) | q t / rho c | FEM/(q t/rho c) | FEM/series - 1")
    snap_steps = sorted(res.snapshots)
    origin = np.zeros((1, 3))
    hist = []
    for stp in snap_steps:
        t = stp * T_STAR / N_MAIN
        f_c = float(_at(res.snapshots[stp], origin, T_INF)[0])
        s_c = 0.0 if stp == 0 else float(np.real(m1._trans_ex(t)(np.zeros(1), np.zeros(1), np.zeros(1))[0]))
        ad = Q * t / m1.RHO_C
        hist.append((stp, t, f_c, s_c, ad))
        r1 = f"{f_c / ad:.6f}" if ad > 0 else "n/a"
        r2 = f"{f_c / s_c - 1.0:+.3e}" if s_c != 0 else "n/a"
        say(f"  {stp:3d} | {t:9.1f} | {f_c:.6f} | {s_c:.6f} | {ad:.6f} | {r1} | {r2}")
    say(f"[history] predicted: FEM ≈ series to ≲ 0.1 %; centre/(q t/rho c) falls from 1 to 0.989 at t*")

    gam1 = [float(m0._axis_data(L, 1)[0][0]) for L in LHALF]
    tau1 = 1.0 / (m0.ALPHA * sum((g / L) ** 2 for g, L in zip(gam1, LHALF)))
    say(f"[tau] gamma_1 per axis {gam1[0]:.6f} / {gam1[1]:.6f} / {gam1[2]:.6f} "
        f"(records {REC['gam1'][0]} / {REC['gam1'][1]} / {REC['gam1'][2]}, {W0}:91); "
        f"tau_1 = {tau1:.1f} s = {tau1 / 3600.0:.2f} h (record {REC['tau1']:.0f} s, {W0}:94)")
    th_ss0 = th_c_ex
    for lab, t in (("t*", T_STAR), ("t_C", m1.T_C), ("tau_1", tau1), ("3 tau_1", 3.0 * tau1)):
        v = float(np.real(m1._trans_ex(t)(np.zeros(1), np.zeros(1), np.zeros(1))[0]))
        extra = ""
        if lab == "t*":
            extra = f" (record {REC['frac_tstar']})"
        elif lab == "t_C":
            extra = f" (record {REC['frac_tC']})"
        say(f"[tau] series theta(0, {lab} = {t:.1f} s) = {v:.6f} K = {v / th_ss0:.4f} of theta_ss(0){extra}")
    clk.lap("history")

    # ---- outputs: CSV, then XDMF, then adopt_host_ownership ---------------------------
    OUTPUT_DIR.mkdir(exist_ok=True)
    if rank0:
        csv_path = OUTPUT_DIR / f"{BASENAME}_centre_lines.csv"
        with open(csv_path, "w") as fh:
            fh.write("line,s_over_L,x_m,y_m,z_m,theta_fem_K,theta_series_K,theta_slab_K,"
                     "slab_over_series_minus_1\n")
            for row in csv_rows:
                fh.write(f"{row[0]},{row[1]:.6f},{row[2]:.9e},{row[3]:.9e},{row[4]:.9e},"
                         f"{row[5]:.9e},{row[6]:.9e},{row[7]:.9e},{row[8]:.9e}\n")
    comm.Barrier()

    V1 = st.temperature.function_space
    zc = V1.tabulate_dof_coordinates()[:, 2]
    f_steady = fem.Function(V1, name="theta_steady_fem")
    f_steady.x.array[:] = np.real(st.temperature.x.array) - T_INF
    f_series = fem.Function(V1, name="theta_steady_series")
    f_series.interpolate(ex)
    f_series.x.array[:] = np.real(f_series.x.array)
    f_diff = fem.Function(V1, name="theta_steady_fem_minus_series")
    f_diff.x.array[:] = np.real(f_steady.x.array) - np.real(f_series.x.array)
    f_slab = fem.Function(V1, name="slab_over_series_minus_1")
    f_slab.x.array[:] = _slab(zc) / np.real(f_series.x.array) - 1.0
    f_tf = fem.Function(V1, name="theta_tstar_fem")
    f_tf.x.array[:] = np.real(res.temperature.x.array) - T_INF
    f_ts = fem.Function(V1, name="theta_tstar_series")
    f_ts.interpolate(ex_t)
    f_ts.x.array[:] = np.real(f_ts.x.array)
    fields = {}
    for f in (f_steady, f_series, f_diff, f_slab, f_tf, f_ts, mask10):
        f.x.scatter_forward()
        fields[f.name] = f
    xdmf_path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_combined", msh, cell_tags, fields, comm=comm, facet_tags=face_tags)
    adopt_host_ownership(OUTPUT_DIR, comm=comm)
    clk.lap("outputs")
    if rank0:
        print(f"\n[output] {xdmf_path.name}: CG1 theta_steady_fem, theta_steady_series, "
              "theta_steady_fem_minus_series, slab_over_series_minus_1, theta_tstar_fem, "
              f"theta_tstar_series (K, t* = {T_STAR:.1f} s), DG0 slab_within_10pct, six face tags")
        print(f"[output] {BASENAME}_centre_lines.csv ({len(csv_rows)} rows)")
        print("[scope] one rung (m = 6), P1, uniform q = 2000 W/m^3, Robin h = 10 on all six faces; "
              "no F2182 compliance, coil, SAR source or duty cycle claimed")
        print("All assertions hold")


if __name__ == "__main__":
    main()
