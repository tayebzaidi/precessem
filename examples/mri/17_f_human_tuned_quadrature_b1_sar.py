"""Example (`EX-109`): the tuned copper F-human birdcage under the clinical
two-port 90 degree drive -- ``|B1+|`` and the 10 g SAR hotspot.

`WF-10` step 2 (`tests/validation/test_wf10_f_human_tuned_b1_sar.py`, item 358,
and `tests/validation/test_wf10_f_human_tuned_sar10g.py`, item 359) drives the
copper hole mesh (369 387 cells, Leontovich wall, kappa off, degree 1, 64 MHz,
r = 0.2 Ohm on all 32 ring gaps) at C_in,F-human = 118.2437 pF with feeds P17
and P21 (slots 0 / 4), then reads ``|B1+|`` and the 10 g lattice.

``build_tuned_quadrature`` is imported (one build serves both readings; no Z50
pair).  The two ``test_*`` bodies are never imported: this script repeats their
*sequence* (provenance: 358 body lines 540-613 and 359 body lines 173-289 at
commit ``7ef44c35``) through the same imported helpers and bands; the printed
rows are compared at the printed digits with the two W1 logs
(``20261010T170356Z_WF-10-step2-W1.log:6308-6341`` and
``20261010T171705Z_WF-10-step2s-W1.log:6312-6383``) in the example-local
``REC`` tables below (EX-89 item 332 (Q1)).  A digit mismatch is reported, never
a stop; every asserted band is an imported one.

Negative controls (asserted): (N1) the mis-paired sense misses the mirror by
more than ``C4_COVARIANCE_BAND``; (N2) the single feed P17 without the partner
swap misses the mirror by more than the band; the node ``P_surf``-blind miss is
asserted inside the builder.

Scope: degree 1, 64 MHz, kappa off, one copper hole fixture, one phantom.  Not
the four-port drive (`EX-75`), not f0 / Q / a match / S11 / Z_in, not
homogeneity, not an absolute |B1+| or SAR, not C16 invariance, not a global
SAR maximum, not compliance, not degree 2, not AED.

Run through the example runner (the ``mri:`` group sources the complex build)::

    ./scripts/run_examples.sh -e mri:17 -n 8 -t 590

Output: ``examples/mri/paraview_output/mri_17_f_human_tuned_quadrature_b1_sar_combined.xdmf``
and ``..._lattice.csv``.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import ufl
from dolfinx import default_scalar_type, fem, la
from dolfinx.mesh import compute_midpoints
from mpi4py import MPI

from fem_em_solver.io.paraview_utils import adopt_host_ownership, write_xdmf_with_tags
from fem_em_solver.post.faraday import (
    b1_minus,
    b1_plus,
    magnetic_flux_density_from_e,
    project_to_cg1,
)
from fem_em_solver.post.sar import (
    averaging_ball_radius,
    build_density_field,
    mass_averaged_sar,
    mean_sar,
)
from fem_em_solver.post.setup_figure import write_setup_figure

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

import os  # noqa: E402

from tests.mesh.helpers import global_cell_tag_set  # noqa: E402
from tests.mesh.test_birdcage_port_sheet_prerequisite import CELL_COUNT_BAND  # noqa: E402
from tests.mesh.test_birdcage_port_sheets import PORT_LOWER, PORT_UPPER  # noqa: E402
from tests.validation.test_birdcage_b1_plus_map import (  # noqa: E402
    C4_COVARIANCE_BAND,
    PHANTOM_RHO_KG_PER_M3,
    _relative_l2,
    _rotate_z,
    _sample_points,
)
from tests.validation.test_birdcage_b1_quadrature import _mirror_xy, _read_senses  # noqa: E402
from tests.validation.test_birdcage_f_human_rung import F_HUMAN_RING_RADIUS  # noqa: E402
from tests.validation.test_birdcage_power_identity import DISCRETE_IDENTITY_RTOL  # noqa: E402
from tests.validation.test_birdcage_sar_mass_averaged import EXACT_IDENTITY_RTOL  # noqa: E402
from tests.validation.test_mass_averaged_sar_standard_masses import (  # noqa: E402
    QUADRATURE_DEGREE,
    TEN_GRAM_KG,
)
from tests.validation.test_port_birdcage_lumped_column import (  # noqa: E402
    CONDUCTOR_CELL_TAG,
    PHANTOM_CELL_TAG,
)
from tests.validation.test_port_birdcage_ring_matrix import REUSE_REPRODUCTION_RTOL  # noqa: E402
from tests.validation.test_port_circuit_layer_field import _residual  # noqa: E402
from tests.validation.test_port_drive_superposition import UNIT_WEIGHT_CONTROL_RTOL  # noqa: E402
from tests.validation.test_wf7_f_human_sar10g_hotspot import (  # noqa: E402
    CONDUCTOR_CLEARANCE_M,
    COVERAGE_BALL_PAD_M,
    FIRST_BALL_STOP_S,
    LATTICE_R_M,
    _phantom_geometry,
)
from tests.validation.test_wf10_f_human_copper_esr_root import Y0_CODING_BAND  # noqa: E402
from tests.validation.test_wf10_f_human_tuned_b1_sar import (  # noqa: E402
    C_IN_F_HUMAN_PF,
    C_IN_SECANT_PF,
    NODE_JSON,
    PREDICTED_C16,
    PREDICTED_MIRROR as B1_PREDICTED_MIRROR,
    PREDICTED_MISPAIRED,
    PREDICTED_SENSE_RATIO,
    _signed_k_shares,
    build_tuned_quadrature,
)
from tests.validation.test_wf10_f_human_tuned_circuit import _k1_shares  # noqa: E402
from tests.validation.test_wf10_f_human_tuned_sar10g import (  # noqa: E402
    PREDICTED_AZ_SPREAD,
    PREDICTED_CONTROL_MISS,
    PREDICTED_MIRROR as SAR_PREDICTED_MIRROR,
    _containment,
    _lattice,
    _plane_distance,
)
from tests.validation.test_wf12_f_human_copper_matrix import HOLE_CELL_RECORD  # noqa: E402
from tests.validation.test_wf12_f_human_hole import _rss_gib  # noqa: E402

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"  # committed, EX-57
BASENAME = "mri_17_f_human_tuned_quadrature_b1_sar"
W1_B1 = "20261010T170356Z_WF-10-step2-W1.log"
W1_SAR = "20261010T171705Z_WF-10-step2s-W1.log"

# ---- the W1 records, restated as EX-89 item 238's REC was: each commented with
# its log:line, compared at the printed precision (a mismatch is reported). ----
REC = {
    # W1_B1 :6325-6328 (Y_in), formatted "%+.9e" of re, im
    "Y17_cw": ("+1.172678791e-02", "+2.783852939e-06"),
    "Y21_cw": ("+1.021104330e-02", "-1.869819907e-04"),
    "Y17_ccw": ("+1.010006039e-02", "+1.017991223e-05"),
    "Y21_ccw": ("+1.184992763e-02", "-2.164134517e-04"),
    # W1_B1 :6329
    "Ypair_a": "2.2401e-02", "Ypair_b": "2.1213e-02",
    # W1_B1 :6330-6331 (signed k shares "%.6f": ring+ lead, lag, ring- lead, lag, unsigned ring+, ring-)
    "k_ccw": ("0.002182", "0.996265", "0.002173", "0.991672", "0.998448", "0.993845"),
    "k_cw": ("0.996794", "0.002035", "0.992216", "0.002026", "0.998829", "0.994241"),
    # W1_B1 :6332-6339
    "n_pts": (41, 41),
    "sense": "21.5155", "sense_cw": "22.2551",
    "purity_ccw": "21.5146", "purity_cw": "20.1394",
    "mirror_b": "1.9647e-03", "mispaired": "9.5476e-01", "mispaired_pred": "9.5353e-01",
    "wrong_plane": "5.2887e-03", "c16_worst": ("8", "8.0260e-03"),
    "mean_b": "2.884091e-07", "cv_b": "7.3397e-03", "b_per_sqrt_w": "9.803471e-07",
    # W1_B1 :6341
    "sar_mean": "3.708285e-03", "sar_mean_per_w": "4.284647e-02", "pph_over_psrc": "9.472043e-02",
    # W1_SAR :6354-6371 (ccw, cw(M) "%.9e", rel "%.4f" %)
    "lattice": [
        ("r=15 dz=+0 az17+0", "6.830007305e-04", "6.606850198e-04", "3.2673"),
        ("r=30 dz=+0 az17+0", "2.230610345e-03", "2.138875190e-03", "4.1126"),
        ("r=45 dz=+0 az17+0", "4.852895959e-03", "4.821772251e-03", "0.6413"),
        ("r=15 dz=+0 az17+90", "5.781188640e-04", "5.968531327e-04", "3.2406"),
        ("r=30 dz=+0 az17+90", "1.803115212e-03", "1.887962586e-03", "4.7056"),
        ("r=45 dz=+0 az17+90", "4.026865820e-03", "4.073878924e-03", "1.1675"),
        ("r=15 dz=+0 az17+180", "6.587520880e-04", "6.671098227e-04", "1.2687"),
        ("r=30 dz=+0 az17+180", "2.230637254e-03", "2.224615414e-03", "0.2700"),
        ("r=45 dz=+0 az17+180", "4.816772104e-03", "4.848794891e-03", "0.6648"),
        ("r=15 dz=+0 az17+270", "5.882637672e-04", "5.803609204e-04", "1.3434"),
        ("r=30 dz=+0 az17+270", "1.879465733e-03", "1.896434461e-03", "0.9028"),
        ("r=45 dz=+0 az17+270", "4.063626880e-03", "4.053371643e-03", "0.2524"),
        ("r=15 dz=-50 az17+0", "1.635320677e-03", "1.636807026e-03", "0.0909"),
        ("r=30 dz=-50 az17+0", "2.243679310e-03", "2.244344497e-03", "0.0296"),
        ("r=45 dz=-50 az17+0", "3.414395055e-03", "3.391893557e-03", "0.6590"),
        ("r=15 dz=+50 az17+0", "1.578408674e-03", "1.587542278e-03", "0.5787"),
        ("r=30 dz=+50 az17+0", "2.091815655e-03", "2.096323755e-03", "0.2155"),
        ("r=45 dz=+50 az17+0", "3.209282450e-03", "3.198000962e-03", "0.3515"),
    ],
    # W1_SAR :6372-6376, :6383
    "c_star": "r=45 dz=+0 az17+0",
    "sar_star": "4.852895959e-03",
    "mirror_star": "0.6413", "mirror_l2": "1.2789",
    "spreads": ("16.7267", "20.9986", "18.6041"),
    "ctl_miss": "3674.7755", "ctl_partner": "7.3845",
    "ctl17_c": "1.144006e-04", "ctl17_m": "4.318366e-03", "ctl21_m": "1.059527e-04",
    "sar_per_w": ("5.607159e-02", "5.919694e-01"),
    "worst_containment_mm": "-4.8813", "n_balls": 39,
    "coverage_mm": "108.1429",
}


def _dg0_vec(form):
    """Assemble a DG0-test-function linear form to a per-cell real vector
    (``mri:6``'s ``_dg0_vec``, re-written here)."""
    b = fem.assemble_vector(fem.form(form))
    b.scatter_reverse(la.InsertMode.add)
    b.scatter_forward()
    return np.real(b.array)


def main() -> None:
    comm = MPI.COMM_WORLD
    rank0 = comm.rank == 0
    started = time.perf_counter()

    if not np.issubdtype(np.dtype(default_scalar_type), np.complexfloating):
        raise RuntimeError(
            "this example needs the complex DolfinX build "
            "(source /usr/local/bin/dolfinx-complex-mode); the runner's `mri:` "
            "group sources it automatically"
        )

    if rank0:
        print("=" * 78)
        print("EX-109 -- the tuned copper F-human birdcage, two-port 90 deg drive: "
              "|B1+| and the 10 g hotspot (WF-10 step 2)")
        print("=" * 78)
        print(f"[records] every printed row is beside {W1_B1}:6308-6341 and "
              f"{W1_SAR}:6312-6383", flush=True)

    mism = []  # (label, got, want)

    def chk(label, got, want):
        if got != want:
            mism.append((label, got, want))
        return got == want

    # ---- the stored node record, read BEFORE the builder (it rewrites the file) --
    assert os.path.exists(NODE_JSON), f"{NODE_JSON} absent: stop (never re-solve)"
    node_rec = None
    if rank0:
        with open(NODE_JSON) as fh:
            node_rec = json.load(fh)
    node_rec = comm.bcast(node_rec, root=0)
    assert abs(node_rec["c_in_pf"] - C_IN_F_HUMAN_PF) == 0.0, node_rec["c_in_pf"]
    s_rec = np.array([[complex(a, b) for a, b in row] for row in node_rec["s_node"]],
                     dtype=np.complex128)

    # ---- one build serves both readings ------------------------------------------
    out = build_tuned_quadrature(comm, with_z50=False)
    say, ctx, times, omega = out["say"], out["ctx"], out["times"], out["omega"]
    msh, cell_tags, sigma = ctx["msh"], ctx["cell_tags"], out["sigma"]
    order, slots, feeds = out["order"], out["slots"], out["feeds"]
    parts = out["parts"]
    centre, radius, half_height, lo, hi = _phantom_geometry(msh, cell_tags, PHANTOM_CELL_TAG, comm)
    z_c = float(centre[2])

    # EX-57 setup figure -- right after the fixture returns the built mesh.
    names = {2: "air", PHANTOM_CELL_TAG: "phantom", CONDUCTOR_CELL_TAG: "conductor (hole route: void)"}
    for sp in ctx["specs"]:
        n = int(sp.port_id[1:])
        label = {"P17": "feed P17 (0 deg)", "P21": "feed P21 (90 deg)"}.get(
            sp.port_id, "C_in + 0.2 Ohm gap")
        names[PORT_LOWER + n] = label
        names[PORT_UPPER + n] = label
    write_setup_figure(
        msh,
        cell_tags,
        FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names=names,
        hide_tags=(2,),
        translucent_tags=(PHANTOM_CELL_TAG,),
        slice_normal=(0.0, 0.0, 1.0),
        slice_origin=(0.0, 0.0, z_c),
        title="mri:17 -- tuned copper F-human birdcage, two-port 90 deg feed at C_in 118.2437 pF, "
              "r = 0.2 Ohm (WF-10 H9')",
        window_size=(1400, 630),
        comm=comm,
    )

    # ---- (A0), (S_rep) -----------------------------------------------------------
    cell_rel = out["cells"] / HOLE_CELL_RECORD - 1.0
    s_node = np.asarray(out["s_node"], dtype=np.complex128)
    r_srep = _residual(s_node, s_rec)
    if rank0:
        print(f"\n[anchor] (A0) cells {out['cells']} vs record {HOLE_CELL_RECORD}, rel "
              f"{cell_rel:.3e} (band {CELL_COUNT_BAND}); node C_in {C_IN_F_HUMAN_PF} pF "
              f"(secant {C_IN_SECANT_PF:.6f}); feeds {feeds}, az "
              f"{out['az_feeds'][0]:.6f} / {out['az_feeds'][1]:.6f} deg, mirror plane "
              f"{out['az_mirror']:.6f} deg  (W1: 369387, 0; 11.25 / 101.25; 56.25)")
        print(f"[anchor] (S_rep) s_node vs the stored node JSON: {r_srep:.6e} <= "
              f"{REUSE_REPRODUCTION_RTOL:g}  (W1_SAR :6312: 6.027e-13)", flush=True)

    # ---- Y_in and the signed k spectrum ------------------------------------------
    y_rows = {"cw": ("Y17_cw", "Y21_cw"), "ccw": ("Y17_ccw", "Y21_ccw")}
    if rank0:
        print("\n[Y_in] (readings)")
    for sense in ("cw", "ccw"):
        for i, pid in enumerate(feeds):
            y = out["y_in"][sense][i]
            key = y_rows[sense][i]
            ok = chk(key, (f"{y.real:+.9e}", f"{y.imag:+.9e}"), REC[key])
            if rank0:
                print(f"   [{sense}] {pid}: {y.real:+.9e} {y.imag:+.9e}j S   {'= W1' if ok else 'DIFFERS'}")
    yc, yw = out["y_in"]["ccw"], out["y_in"]["cw"]
    sp_a = abs(yc[0] - yw[1]) / abs(yc[0])
    sp_b = abs(yc[1] - yw[0]) / abs(yc[1])
    chk("Ypair_a", f"{sp_a:.4e}", REC["Ypair_a"])
    chk("Ypair_b", f"{sp_b:.4e}", REC["Ypair_b"])
    if rank0:
        print(f"[Y_in] M-paired spreads {sp_a:.4e}, {sp_b:.4e}  (W1: {REC['Ypair_a']}, {REC['Ypair_b']})")
    for sense in ("ccw", "cw"):
        sh = _signed_k_shares(out["cur"][sense], order, slots)
        un = _k1_shares(out["cur"][sense], order, slots)
        got = (f"{sh[1.0]['lead']:.6f}", f"{sh[1.0]['lag']:.6f}", f"{sh[-1.0]['lead']:.6f}",
               f"{sh[-1.0]['lag']:.6f}", f"{un[1.0]:.6f}", f"{un[-1.0]:.6f}")
        ok = chk(f"k_{sense}", got, REC[f"k_{sense}"])
        if rank0:
            print(f"[k] [{sense}] ring+ lead/lag {got[0]}/{got[1]}, ring- lead/lag {got[2]}/{got[3]}, "
                  f"unsigned ring+ {got[4]} ring- {got[5]}   {'= W1' if ok else 'DIFFERS'}")

    # ---- B, CG1 projection, point reads -------------------------------------------
    comm.Barrier()
    b_dg = {}
    cg1 = {}
    for sense in ("ccw", "cw"):
        b_dg[sense] = magnetic_flux_density_from_e(out["e"][sense], omega, name=f"B_{sense}")
        cg1[sense] = project_to_cg1(b_dg[sense], name=f"B_{sense}_cg1")
    comm.Barrier()
    sweep = {"mesh": msh, "cell_tags": cell_tags}
    points = _sample_points(sweep)
    az_m = out["az_mirror"]
    az17 = out["az_feeds"][0]
    ccw_x = _read_senses(cg1["ccw"], points)
    cw_x = _read_senses(cg1["cw"], points)
    cw_mx = _read_senses(cg1["cw"], _mirror_xy(points, az_m))
    cw_m17 = _read_senses(cg1["cw"], _mirror_xy(points, az17))
    rot = {m: _read_senses(cg1["ccw"], _rotate_z(points, np.radians(m * 22.5)))
           for m in range(1, 16)}
    mask = np.logical_and.reduce([ccw_x[2], cw_x[2], cw_mx[2], cw_m17[2]]
                                 + [r[2] for r in rot.values()])
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
    r_pts = np.where(mask, np.hypot(points[:, 0], points[:, 1]) + np.abs(points[:, 2]), np.inf)
    ic = int(np.argmin(r_pts))
    mean_b = float(np.mean(plus_x[mask]))
    cv_b = float(np.std(plus_x[mask]) / mean_b)
    p_src_ccw = parts["ccw"]["p_src"]
    pur_ccw = f"{plus_x[ic] / minus_x[ic]:.4f}"
    pur_cw = f"{cw_x[1][ic] / cw_x[0][ic]:.4f}"
    chk("n_pts", (n_pts, n_mask), REC["n_pts"])
    chk("sense", f"{sense_ratio:.4f}", REC["sense"])
    chk("sense_cw", f"{cw_ratio:.4f}", REC["sense_cw"])
    chk("purity_ccw", pur_ccw, REC["purity_ccw"])
    chk("purity_cw", pur_cw, REC["purity_cw"])
    chk("mirror_b", f"{mirror:.4e}", REC["mirror_b"])
    chk("mispaired", f"{mispaired:.4e}", REC["mispaired"])
    chk("mispaired_pred", f"{mispaired_pred:.4e}", REC["mispaired_pred"])
    chk("wrong_plane", f"{wrong_plane:.4e}", REC["wrong_plane"])
    chk("c16_worst", (str(worst_m), f"{c16[worst_m]:.4e}"), REC["c16_worst"])
    chk("mean_b", f"{mean_b:.6e}", REC["mean_b"])
    chk("cv_b", f"{cv_b:.4e}", REC["cv_b"])
    chk("b_per_sqrt_w", f"{mean_b / np.sqrt(p_src_ccw):.6e}", REC["b_per_sqrt_w"])
    if rank0:
        print(f"\n[B1+] {n_pts} sample points, {n_mask} valid on every image (W1: 41, 41)")
        print(f"[B1+] (Sense) mean|B1+|_ccw / mean|B1-|_ccw = {sense_ratio:.4f} (> 1 ASSERTED, "
              f"358's bar; PREDICTED >= {PREDICTED_SENSE_RATIO:g}; W1 {REC['sense']}); cw: "
              f"mean|B1-|/mean|B1+| = {cw_ratio:.4f} (W1 {REC['sense_cw']})")
        print(f"[B1+] centre purity at ({points[ic, 0]:.4e}, {points[ic, 1]:.4e}, "
              f"{points[ic, 2]:.4e}) m: ccw |B1+|/|B1-| {pur_ccw}; cw |B1-|/|B1+| {pur_cw}  "
              f"(W1 {REC['purity_ccw']}, {REC['purity_cw']})")
        print(f"[B1+] (M_B) feed-bisector mirror (az {az_m:.6f} deg): |B1-|_cw(Mx) vs "
              f"|B1+|_ccw(x) rel L2 {mirror:.4e} (<= {C4_COVARIANCE_BAND:g} ASSERTED; PREDICTED "
              f"{B1_PREDICTED_MIRROR[0]:g}-{B1_PREDICTED_MIRROR[1]:g}; W1 {REC['mirror_b']})")
        print(f"[B1+] (N1, ASSERTED > band) mis-paired |B1+|_cw(Mx) vs |B1+|_ccw(x) {mispaired:.4e} "
              f"= {mispaired / C4_COVARIANCE_BAND:.1f}x the band (W1 {REC['mispaired']}); "
              f"|B1-|_ccw(x) vs |B1+|_ccw(x) {mispaired_pred:.4e} (PREDICTED >= "
              f"{PREDICTED_MISPAIRED:g}, printed)")
        print(f"[B1+] wrong-plane mirror (az17 {az17:.6f} deg): {wrong_plane:.4e} (printed; W1 "
              f"{REC['wrong_plane']})")
        print("[B1+] C16 readings (printed, never asserted -- the drive is not C16-symmetric): "
              + ", ".join(f"m={m}: {v:.4e}" for m, v in c16.items())
              + f"; worst m={worst_m} {c16[worst_m]:.4e} (PREDICTED <= {PREDICTED_C16:g}; W1 "
              f"m=8 {REC['c16_worst'][1]})")
        print(f"[B1+] mean|B1+|_ccw {mean_b:.6e} T, CV {cv_b:.4e} (reading, no homogeneity claim); "
              f"/sqrt(P_src) {mean_b / np.sqrt(p_src_ccw):.6e} T/sqrt(W)  (W1 {REC['mean_b']}, "
              f"{REC['cv_b']}, {REC['b_per_sqrt_w']})", flush=True)

    # ---- SAR on the tag, (Sigma_P), (X1) ------------------------------------------
    sar = mean_sar(out["e"]["ccw"], sigma=sigma, rho=PHANTOM_RHO_KG_PER_M3,
                   cell_tags=cell_tags, subdomain_ids=PHANTOM_CELL_TAG, comm=comm)
    p_ph = parts["ccw"]["p_phantom"]
    r_sigma = abs(sar["dissipated_power_w"] - p_ph) / abs(p_ph)

    dg0 = fem.functionspace(msh, ("DG", 0))
    v = ufl.TestFunction(dg0)
    dx = ufl.Measure("dx", domain=msh, subdomain_data=cell_tags)
    n_loc = dg0.dofmap.index_map.size_local
    e_ccw = out["e"]["ccw"]
    vol = _dg0_vec(ufl.conj(v) * ufl.dx(domain=msh))
    loss = _dg0_vec(0.5 * sigma * ufl.inner(e_ccw, e_ccw) * ufl.conj(v) * dx(PHANTOM_CELL_TAG))
    x1_sum = float(comm.allreduce(np.sum(loss[:n_loc]), op=MPI.SUM))
    r_x1 = abs(x1_sum / p_ph - 1.0)
    chk("sar_mean", f"{sar['mean_sar_w_per_kg']:.6e}", REC["sar_mean"])
    chk("sar_mean_per_w", f"{sar['mean_sar_w_per_kg'] / p_src_ccw:.6e}", REC["sar_mean_per_w"])
    chk("pph_over_psrc", f"{p_ph / p_src_ccw:.6e}", REC["pph_over_psrc"])
    if rank0:
        print(f"\n[SAR] (Sigma_P) P_phantom {p_ph:.12e} W vs mean_sar dissipated "
              f"{sar['dissipated_power_w']:.12e} W, rel {r_sigma:.3e} (<= {EXACT_IDENTITY_RTOL:g} "
              f"ASSERTED; W1 0.000e+00)")
        print(f"[SAR] (X1) per-cell phantom loss vector summed over owned cells {x1_sum:.12e} W vs "
              f"P_phantom, rel {r_x1:.3e} (<= {EXACT_IDENTITY_RTOL:g} ASSERTED, example-side)")
        print(f"[SAR] whole-phantom mean {sar['mean_sar_w_per_kg']:.6e} W/kg, per W of P_src "
              f"{sar['mean_sar_w_per_kg'] / p_src_ccw:.6e}; P_phantom/P_src {p_ph / p_src_ccw:.6e}  "
              f"(W1 {REC['sar_mean']}, {REC['sar_mean_per_w']}, {REC['pph_over_psrc']})", flush=True)

    # ---- the 10 g lattice, repeating 359's sequence ---------------------------------
    a10 = averaging_ball_radius(mass_kg=TEN_GRAM_KG, rho=PHANTOM_RHO_KG_PER_M3)
    all_tags = global_cell_tag_set(msh, cell_tags)
    rho_field = build_density_field(
        msh, PHANTOM_RHO_KG_PER_M3, cell_tags=cell_tags,
        density_map={int(t): 0.0 for t in sorted(all_tags) if int(t) != PHANTOM_CELL_TAG})
    fields = {"ccw": out["e"]["ccw"], "cw": out["e"]["cw"], "P17": out["e17"], "P21": out["e21"]}
    containment, walls = [], []

    def ball(key, c, label):
        c = tuple(float(x) for x in c)
        containment.append((label, c, _containment(c, centre, radius, half_height, a10)))
        t = time.perf_counter()
        r = mass_averaged_sar(fields[key], sigma=sigma, rho=rho_field, center=c, radius=a10,
                              comm=comm, quadrature_degree=QUADRATURE_DEGREE)
        r["wall_s"] = comm.allreduce(time.perf_counter() - t, op=MPI.MAX)
        walls.append(r["wall_s"])
        return r

    comm.Barrier()
    t_sar = time.perf_counter()
    lattice = _lattice(az17, z_c)
    ccw, reduced = {}, False
    for idx, (label, _r, dz, _d, c) in enumerate(lattice):
        if reduced and dz != 0.0:
            continue
        ccw[label] = ball("ccw", c, label)
        if idx == 0 and ccw[label]["wall_s"] > FIRST_BALL_STOP_S:
            reduced = True
            if rank0:
                print(f"[SAR] REDUCED LATTICE: first ball {ccw[label]['wall_s']:.2f} s > "
                      f"{FIRST_BALL_STOP_S:.0f} s -- z = z_c rows only (pre-registered, not a red)")
    coords = {label: np.array(c) for label, *_x, c in lattice}
    used = [label for label, *_x in lattice if label in ccw]
    key_star = max(used, key=lambda k: ccw[k]["averaged_sar_w_per_kg"])
    sar_star = ccw[key_star]["averaged_sar_w_per_kg"]
    c_star = coords[key_star]
    cw = {}
    for label in used:
        mc = _mirror_xy(coords[label][None, :], az_m)[0]
        cw[label] = ball("cw", mc, f"M({label})")
    mc_star = _mirror_xy(c_star[None, :], az_m)[0]
    ctl17_c = ball("P17", c_star, "control P17 c*")["averaged_sar_w_per_kg"]
    ctl17_m = ball("P17", mc_star, "control P17 Mc*")["averaged_sar_w_per_kg"]
    ctl21_m = ball("P21", mc_star, "control P21 Mc*")["averaged_sar_w_per_kg"]
    coverage_radius = float(np.hypot(radius, half_height) + COVERAGE_BALL_PAD_M)
    conductor_limit = F_HUMAN_RING_RADIUS - CONDUCTOR_CLEARANCE_M
    whole = mass_averaged_sar(fields["ccw"], sigma=sigma, rho=rho_field,
                              center=tuple(float(x) for x in centre), radius=coverage_radius,
                              comm=comm, quadrature_degree=QUADRATURE_DEGREE)
    t_sar = comm.allreduce(time.perf_counter() - t_sar, op=MPI.MAX)
    power_miss = abs(whole["dissipated_power_w"] / sar["dissipated_power_w"] - 1.0)
    mass_miss = abs(whole["mass_kg"] / (PHANTOM_RHO_KG_PER_M3 * sar["volume_m3"]) - 1.0)

    sar_c = np.array([ccw[k]["averaged_sar_w_per_kg"] for k in used])
    sar_mm = np.array([cw[k]["averaged_sar_w_per_kg"] for k in used])
    mask18 = np.ones(len(used), dtype=bool)
    mirror_star = abs(cw[key_star]["averaged_sar_w_per_kg"] - sar_star) / sar_star
    mirror_l2 = _relative_l2(sar_mm, sar_c, mask18)
    ctl_miss = abs(ctl17_m - ctl17_c) / ctl17_c
    ctl_partner = abs(ctl21_m - ctl17_c) / ctl17_c
    spreads = {}
    for r in LATTICE_R_M:
        vals = [ccw[lab]["averaged_sar_w_per_kg"] for lab, rr, dz, _d, _c in lattice
                if rr == r and dz == 0.0]
        spreads[r] = (max(vals) - min(vals)) / np.mean(vals)
    worst_cont = max(x[2] for x in containment)

    rec_lat = {row[0]: row[1:] for row in REC["lattice"]}
    n_lat_equal = 0
    for k in used:
        d = abs(cw[k]["averaged_sar_w_per_kg"] / ccw[k]["averaged_sar_w_per_kg"] - 1.0)
        got = (f"{ccw[k]['averaged_sar_w_per_kg']:.9e}", f"{cw[k]['averaged_sar_w_per_kg']:.9e}",
               f"{d * 100:.4f}")
        n_lat_equal += int(chk(f"lattice {k}", got, rec_lat[k]))
    chk("c_star", key_star, REC["c_star"])
    chk("sar_star", f"{sar_star:.9e}", REC["sar_star"])
    chk("mirror_star", f"{mirror_star * 100:.4f}", REC["mirror_star"])
    chk("mirror_l2", f"{mirror_l2 * 100:.4f}", REC["mirror_l2"])
    chk("spreads", tuple(f"{v * 100:.4f}" for v in spreads.values()), REC["spreads"])
    chk("ctl_miss", f"{ctl_miss * 100:.4f}", REC["ctl_miss"])
    chk("ctl_partner", f"{ctl_partner * 100:.4f}", REC["ctl_partner"])
    chk("ctl17_c", f"{ctl17_c:.6e}", REC["ctl17_c"])
    chk("ctl17_m", f"{ctl17_m:.6e}", REC["ctl17_m"])
    chk("ctl21_m", f"{ctl21_m:.6e}", REC["ctl21_m"])
    sar_pw = (f"{sar_star / p_src_ccw:.6e}", f"{sar_star / p_ph:.6e}")
    chk("sar_per_w", sar_pw, REC["sar_per_w"])
    chk("worst_containment", f"{worst_cont * 1e3:.4f}", REC["worst_containment_mm"])
    chk("n_balls", len(containment), REC["n_balls"])
    chk("coverage", f"{coverage_radius * 1e3:.4f}", REC["coverage_mm"])
    _m, _rss_per, rss_sum = _rss_gib(comm)

    if rank0:
        print(f"\n[10 g] ccw table ({len(used)} points{', REDUCED' if reduced else ''}); "
              f"cw at the feed-bisector image; W1 equal on {n_lat_equal} / {len(used)} rows:")
        for k in used:
            d = abs(cw[k]["averaged_sar_w_per_kg"] / ccw[k]["averaged_sar_w_per_kg"] - 1.0)
            print(f"   {k:<26s} ccw {ccw[k]['averaged_sar_w_per_kg']:.9e} cw(M) "
                  f"{cw[k]['averaged_sar_w_per_kg']:.9e} W/kg, rel {d * 100:.4f}%"
                  f"{'   <- c*' if k == key_star else ''}")
        print(f"[10 g] c* = {key_star} ({c_star[0] * 1e3:.3f}, {c_star[1] * 1e3:.3f}, "
              f"{c_star[2] * 1e3:.3f}) mm, SAR10g(c*) {sar_star:.9e} W/kg (W1 {REC['c_star']}, "
              f"{REC['sar_star']}); a lattice argmax, not a global maximum")
        print(f"[10 g] (iii-a) mirror at c*: {mirror_star * 100:.4f}% (<= "
              f"{C4_COVARIANCE_BAND * 100:.1f}% ASSERTED; PREDICTED {SAR_PREDICTED_MIRROR[0] * 100:g}-"
              f"{SAR_PREDICTED_MIRROR[1] * 100:g}%; W1 {REC['mirror_star']}%)")
        print(f"[10 g] (iii-b) lattice rel L2 over {len(used)} points: {mirror_l2 * 100:.4f}% (<= "
              f"{C4_COVARIANCE_BAND * 100:.1f}% ASSERTED; W1 {REC['mirror_l2']}%)")
        print(f"[10 g] (N2, ASSERTED > band) single feed P17 without the partner swap: "
              f"SAR10g(c*) {ctl17_c:.6e}, SAR10g(Mc*) {ctl17_m:.6e} W/kg, miss {ctl_miss * 100:.4f}% "
              f"= {ctl_miss / C4_COVARIANCE_BAND:.0f}x the band (W1 {REC['ctl_miss']}%; "
              f"PREDICTED >= {PREDICTED_CONTROL_MISS * 100:g}%)")
        print(f"[10 g] M-partner P21(Mc*) {ctl21_m:.6e} vs P17(c*): {ctl_partner * 100:.4f}% "
              f"(printed; W1 {REC['ctl_partner']}%, OUT of 359's predicted 5%, a record)")
        print("[10 g] azimuthal spread of the four z_c rows (printed, a record; PREDICTED <= "
              f"{PREDICTED_AZ_SPREAD * 100:g}%): "
              + ", ".join(f"r={r * 1e3:.0f} mm {v * 100:.4f}%" for r, v in spreads.items())
              + f"  (W1 {', '.join(REC['spreads'])})")
        print(f"[10 g] (i) coverage ball {coverage_radius * 1e3:.4f} mm (< conductor limit "
              f"{conductor_limit * 1e3:.1f} mm ASSERTED): power rel {power_miss:.3e}; mass rel "
              f"{mass_miss:.3e} (<= {EXACT_IDENTITY_RTOL:g} ASSERTED; W1 2.442e-15 / 6.994e-15)")
        print(f"[10 g] (iv) worst containment {worst_cont * 1e3:.4f} mm over {len(containment)} balls "
              f"(<= 0 ASSERTED; W1 {REC['worst_containment_mm']} mm over {REC['n_balls']})")
        print(f"[10 g] SAR per watt (ccw): SAR10g(c*)/P_src {sar_pw[0]}, /P_phantom {sar_pw[1]} W/kg/W "
              f"(W1 {REC['sar_per_w'][0]}, {REC['sar_per_w'][1]}); peak-to-mean "
              f"{sar_star / sar['mean_sar_w_per_kg']:.4f}", flush=True)

    # ---- asserts on the imported bands (before the artifacts; rank-consistent values) --
    for sense in ("ccw", "cw"):
        r_lin, r_y = out["anchors"]["lin"][sense]
        rel_w, blind_w = out["anchors"]["p_w"][sense]
        if rank0:
            print(f"[anchor] [{sense}] (P_w) {rel_w:.3e} <= {DISCRETE_IDENTITY_RTOL:g}, blind "
                  f"{blind_w:.4e}; (Lin) currents {r_lin:.3e} <= {UNIT_WEIGHT_CONTROL_RTOL:g}, Y_in "
                  f"{r_y:.3e} <= {Y0_CODING_BAND:g}  (W1: 5.065e-13 / 5.893e-13)")
        assert rel_w <= DISCRETE_IDENTITY_RTOL, f"(P_w) {sense} {rel_w:.3e}"
        assert r_lin <= UNIT_WEIGHT_CONTROL_RTOL, f"(Lin) {sense} {r_lin:.3e}"
        assert r_y <= Y0_CODING_BAND, f"(Lin) Y_in {sense} {r_y:.3e}"
    if rank0:
        print(f"[anchor] (P_RC) {out['p_rc_node']:.3e} <= {DISCRETE_IDENTITY_RTOL:g}; node "
              f"P_surf-blind miss {out['blind_node']:.4e} = "
              f"{out['blind_node'] / DISCRETE_IDENTITY_RTOL:.1f}x the band (asserted in the "
              f"builder; W1 7.830e-13, 7012.5x)")
    assert abs(cell_rel) <= CELL_COUNT_BAND and out["cells"] == HOLE_CELL_RECORD, out["cells"]
    assert r_srep <= REUSE_REPRODUCTION_RTOL, f"(S_rep) {r_srep:.6e}"
    assert r_sigma <= EXACT_IDENTITY_RTOL, f"(Sigma_P) {r_sigma:.3e}"
    assert r_x1 <= EXACT_IDENTITY_RTOL, f"(X1) {r_x1:.3e}"
    assert n_mask > 0, "no sample point valid on every image"
    assert sense_ratio > 1.0, f"(Sense) red: B1- dominant under ccw ({sense_ratio:.4f})"
    assert mirror <= C4_COVARIANCE_BAND, f"(M_B) mirror {mirror:.4e}"
    assert mispaired > C4_COVARIANCE_BAND, f"(N1) mis-paired control inside the band {mispaired:.4e}"
    assert coverage_radius < conductor_limit, (coverage_radius, conductor_limit)
    assert power_miss <= EXACT_IDENTITY_RTOL, f"(i) power {power_miss:.3e}"
    assert mass_miss <= EXACT_IDENTITY_RTOL, f"(i) mass {mass_miss:.3e}"
    assert worst_cont <= 0.0, f"(iv) containment {worst_cont:.4e} m"
    assert mirror_star <= C4_COVARIANCE_BAND, f"(iii-a) {mirror_star:.4e}"
    assert mirror_l2 <= C4_COVARIANCE_BAND, f"(iii-b) {mirror_l2:.4e}"
    assert ctl_miss > C4_COVARIANCE_BAND, f"(N2) single-feed control inside the band {ctl_miss:.4e}"

    # ---- ParaView maps: DG0, real, named once ------------------------------------
    def named(name, values):
        f = fem.Function(dg0, name=name)
        f.x.array[:] = values
        f.x.scatter_forward()
        return f

    def real_dg0(fn, name):
        f = fem.Function(fn.function_space, name=name)
        f.x.array[:] = np.real(fn.x.array)
        f.x.scatter_forward()
        return f

    bp_ccw = real_dg0(b1_plus(b_dg["ccw"]), "B1_plus_ccw")
    bm_ccw = real_dg0(b1_minus(b_dg["ccw"]), "B1_minus_ccw")
    bm_cw = real_dg0(b1_minus(b_dg["cw"]), "B1_minus_cw")
    purity = np.zeros_like(bp_ccw.x.array)
    np.divide(bp_ccw.x.array, bm_ccw.x.array, out=purity, where=bm_ccw.x.array > 0.0)
    e_sq = fem.Function(dg0, name="e_sq_ccw")
    e_sq.interpolate(fem.Expression(ufl.inner(e_ccw, e_ccw), dg0.element.interpolation_points))
    e_mag = np.sqrt(np.maximum(np.real(e_sq.x.array), 0.0))
    sar_local = np.zeros_like(loss)
    np.divide(loss, PHANTOM_RHO_KG_PER_M3 * vol, out=sar_local, where=(vol > 0.0) & (loss > 0.0))

    tdim = msh.topology.dim
    imap = msh.topology.index_map(tdim)
    n_cells = imap.size_local + imap.num_ghosts
    mids = compute_midpoints(msh, tdim, np.arange(n_cells, dtype=np.int32))
    cell_dof = np.asarray(dg0.dofmap.list).reshape(-1)[:n_cells]
    lat_marker = np.zeros_like(loss)
    c_marker = np.zeros_like(loss)
    side = np.ones_like(loss)
    nrm = np.array([-np.sin(np.radians(az_m)), np.cos(np.radians(az_m)), 0.0])
    side[cell_dof] = np.where(mids @ nrm < 0.0, -1.0, 1.0)
    for k, label in enumerate(label for label, *_x in lattice):
        if label not in ccw:
            continue
        inside = np.linalg.norm(mids - coords[label][None, :], axis=1) < a10
        lat_marker[cell_dof[inside]] = float(k + 1)
    c_marker[cell_dof[np.linalg.norm(mids - c_star[None, :], axis=1) < a10]] = 1.0

    export = {
        "B1_plus_ccw": bp_ccw,
        "B1_minus_ccw": bm_ccw,
        "B1_minus_cw": bm_cw,
        "B1_purity_ccw": named("B1_purity_ccw", purity),
        "E_magnitude_ccw": named("E_magnitude_ccw", e_mag),
        "phantom_loss_per_cell_W_ccw": named("phantom_loss_per_cell_W_ccw", loss),
        "SAR_local_ccw": named("SAR_local_ccw", sar_local),
        "LatticeBall": named("LatticeBall", lat_marker),
        "CStar": named("CStar", c_marker),
        "MirrorSide": named("MirrorSide", side),
    }
    OUTPUT_DIR.mkdir(exist_ok=True)
    xdmf_path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_combined", msh, cell_tags, export, comm=comm
    )
    adopt_host_ownership(OUTPUT_DIR, comm=comm)

    csv_path = OUTPUT_DIR / f"{BASENAME}_lattice.csv"
    if rank0:
        with open(csv_path, "w") as fh:
            fh.write("label,r_m,dz_m,d_az_deg,x_m,y_m,z_m,SAR10g_ccw,SAR10g_cw_M,rel,mass_kg,"
                     "wall_s,is_c_star\n")
            for label, r_, dz_, daz_, c_ in lattice:
                if label not in ccw:
                    continue
                rel = abs(cw[label]["averaged_sar_w_per_kg"] / ccw[label]["averaged_sar_w_per_kg"] - 1.0)
                fh.write(f"\"{label}\",{r_:.6e},{dz_:.6e},{daz_:g},{c_[0]:.9e},{c_[1]:.9e},"
                         f"{c_[2]:.9e},{ccw[label]['averaged_sar_w_per_kg']:.9e},"
                         f"{cw[label]['averaged_sar_w_per_kg']:.9e},{rel:.9e},"
                         f"{ccw[label]['mass_kg']:.9e},{ccw[label]['wall_s']:.3f},"
                         f"{int(label == key_star)}\n")

    if rank0:
        print(f"\n[records] {len(mism)} printed record row(s) differ from the W1 digits"
              + ("" if not mism else ":"))
        for label, got, want in mism:
            print(f"   DIFFERS {label}: got {got}, W1 {want}")
        print(f"[price] {len(walls)} 10 g balls + coverage, per ball mean {np.mean(walls):.2f} s "
              f"max {np.max(walls):.2f} s; SAR phase {t_sar:.1f} s; build {times['build']:.1f} s, "
              f"node pair {times['node_pair']:.1f} s, superposition {times['superposition']:.1f} s; "
              f"RSS summed {rss_sum:.2f} GiB; -n {comm.size}")
        print(f"\n[output] {xdmf_path.name}: DG0 " + ", ".join(export) + ", CellTags")
        print(f"[output] {csv_path.name}")
        print("[paraview] slice z = z_c on `B1_plus_ccw`; threshold `LatticeBall` >= 1 coloured by "
              "`SAR_local_ccw`; compare `B1_plus_ccw` with `B1_minus_cw` across `MirrorSide`.")
        print(f"\nEX-109 anchors, identities and both negative controls: PASS. Total elapsed "
              f"{time.perf_counter() - started:.1f} s.", flush=True)


if __name__ == "__main__":
    main()
