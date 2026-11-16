"""Example (`EX-75`): the tuned copper F-human birdcage under the four-port
0 / 90 / 180 / 270 degree quadrature drive, beside the clinical two-port drive.

`WF-10` step 2 (`tests/validation/test_wf10_f_human_tuned_b1_sar.py`, item 358)
drives the copper hole mesh (369 387 cells, Leontovich wall, kappa off,
degree 1, 64 MHz, r = 0.2 Ohm on all 32 ring gaps) at C_in,F-human with the
two feeds P17 / P21 (slots 0 / 4).  This example builds that fixture once
(``build_tuned_quadrature``, imported, ``with_z50=False``), re-asserts the
two-port gates, then loads the same C with four Thevenin feeds at slots
0 / 4 / 8 / 12 on the feed ring (the other 28 gaps at Z_RC): one factorising
solve plus three held solves, superposed with
``quadrature_phase_weights([0, 1, 2, 3], "ccw")``.

Provenance: the four-port solve repeats the sequence of 358's ``_drive``
closure (lines 364-413 of the gate module) and the sense superposition
(lines 524-571) through imported helpers; the two ``test_*`` bodies are never
imported.  The four-port drive runs UNTUNED at the two-port root C: its figures
are printed, not gated.  Only the construction identities (X4-P) / (X4-L) /
(X4-K) of the script's own solve are asserted (review ruling 2026-10-11).
The two-port rows are printed beside ``20261010T170356Z_WF-10-step2-W1.log``
:6308-6341 in the example-local ``REC`` table (EX-89 item 332 (Q1)); a digit
mismatch is reported, never a stop.

Negative controls: (N1, asserted) the mis-paired sense misses the mirror by
more than ``C4_COVARIANCE_BAND``; (W0, asserted) the two-port weights feed
k = 0 (FFT bin 0 = sqrt 2) and the four-port weights do not.

Scope: degree 1, 64 MHz, kappa off, one copper hole fixture, one phantom,
C fixed at the two-port root.  Not a four-port tuning, not f0 / Q / match /
S11 / Z_in, not homogeneity, not an absolute |B1+|, not C16 or C4 invariance
as a gate, not 10 g SAR, not degree 2, not AED, not compliance.

Run through the example runner (the ``mri:`` group sources the complex build)::

    ./scripts/run_examples.sh -e mri:18 -n 8 -t 590

Output: ``examples/mri/paraview_output/mri_18_f_human_four_port_quadrature_combined.xdmf``
and ``mri_18_f_human_four_port_quadrature_kspectrum.csv``.
"""

from __future__ import annotations

import dataclasses
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import ufl
from dolfinx import default_scalar_type, fem
from mpi4py import MPI

from fem_em_solver.core.cavity import surface_resistance_ohm
from fem_em_solver.io.paraview_utils import adopt_host_ownership, write_xdmf_with_tags
from fem_em_solver.ports.lumped import (
    lumped_port_bilinear_term,
    lumped_port_linear_term,
    sheet_terminal_current,
)
from fem_em_solver.post.faraday import (
    b1_minus,
    b1_plus,
    magnetic_flux_density_from_e,
    project_to_cg1,
)
from fem_em_solver.post.sar import mean_sar
from fem_em_solver.post.setup_figure import write_setup_figure

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from tests.mesh.test_birdcage_port_sheet_prerequisite import CELL_COUNT_BAND  # noqa: E402
from tests.mesh.test_birdcage_port_sheets import PORT_LOWER, PORT_UPPER  # noqa: E402
from tests.validation.test_birdcage_b1_plus_map import (  # noqa: E402
    C4_COVARIANCE_BAND,
    PHANTOM_RHO_KG_PER_M3,
    _relative_l2,
    _rotate_z,
    _sample_points,
)
from tests.validation.test_birdcage_b1_quadrature import (  # noqa: E402
    _mirror_xy,
    _read_senses,
    quadrature_phase_weights,
)
from tests.validation.test_birdcage_power_identity import (  # noqa: E402
    DISCRETE_IDENTITY_RTOL,
    _sheet_field_dissipation_w,
    _source_power_w,
)
from tests.validation.test_birdcage_sar_mass_averaged import EXACT_IDENTITY_RTOL  # noqa: E402
from tests.validation.test_port_birdcage_lumped_column import (  # noqa: E402
    CONDUCTOR_CELL_TAG,
    PHANTOM_CELL_TAG,
)
from tests.validation.test_port_birdcage_ring_column import AZIMUTH_MATCH_DEG  # noqa: E402
from tests.validation.test_port_birdcage_ring_matrix import REUSE_REPRODUCTION_RTOL  # noqa: E402
from tests.validation.test_port_circuit_layer_field import _residual  # noqa: E402
from tests.validation.test_port_drive_superposition import (  # noqa: E402
    UNIT_WEIGHT_CONTROL_RTOL,
    _ring_quadrature_slots,
)
from tests.validation.test_port_package_sparameters import REFERENCE_IMPEDANCE_OHM  # noqa: E402
from tests.validation.test_th14_birdcage_copper import _leontovich_term, _surface_loss_w  # noqa: E402
from tests.validation.test_wf7_f_human_b1_quadrature import FREQUENCY_HZ  # noqa: E402
from tests.validation.test_wf7_f_human_sar10g_hotspot import _phantom_geometry  # noqa: E402
from tests.validation.test_wf10_f_human_copper_esr_root import (  # noqa: E402
    Y0_CODING_BAND,
    _source_voltage,
    _thevenin_feed_impedance,
    _z_rc,
)
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
from tests.validation.test_wf10_f_human_tuned_circuit import (  # noqa: E402
    QUADRATURE_SLOT_OFFSET,
    _respec,
)
from tests.validation.test_wf12_f_human_copper_matrix import (  # noqa: E402
    COPPER_SIGMA,
    HOLE_CELL_RECORD,
)
from tests.validation.test_wf12_f_human_hole import _rss_gib  # noqa: E402

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"  # committed, EX-57
BASENAME = "mri_18_f_human_four_port_quadrature"
W1_B1 = "20261010T170356Z_WF-10-step2-W1.log"

# ---- the W1 records for the two-port drive, restated as EX-89 item 238's REC was:
# each commented with its log:line, compared at the printed precision. ----
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
}


def _spectrum_share(vec16):
    """Power share |X_b|^2 / sum over the 16 FFT bins of a 16-slot vector."""
    spec = np.fft.fft(vec16)
    p = np.abs(spec) ** 2
    tot = float(np.sum(p))
    return p / tot if tot else np.full(16, np.nan)


def _ring_vector(currents, order, slots, sign):
    vec = np.zeros(16, dtype=np.complex128)
    for idx, pid in enumerate(order):
        k, sg = slots[pid]
        if sg == sign:
            vec[k] = currents[idx]
    return vec


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
        print("EX-75 -- the tuned copper F-human birdcage: four-port 0/90/180/270 deg "
              "quadrature drive beside the two-port drive (WF-10 H9')")
        print("=" * 78)
        print(f"[records] two-port rows are beside {W1_B1}:6308-6341; the four-port "
              "figures are printed, not gated", flush=True)

    mism = []  # (label, got, want)

    def chk(label, got, want):
        if got != want:
            mism.append((label, got, want))
        return got == want

    # ---- (W0) drive-pattern coding: ungated numpy, asserted -------------------------
    w4 = quadrature_phase_weights([0, 1, 2, 3], "ccw")
    w4 = w4 / w4[0]
    w2 = quadrature_phase_weights([0, 1], "ccw")
    w2 = w2 / w2[0]
    x4 = np.zeros(16, dtype=np.complex128)
    for n, k in enumerate((0, 4, 8, 12)):
        x4[k] = w4[n]
    x2 = np.zeros(16, dtype=np.complex128)
    x2[0], x2[QUADRATURE_SLOT_OFFSET] = w2[0], w2[1]
    f4, f2 = np.fft.fft(x4), np.fft.fft(x2)
    if rank0:
        print(f"\n[W0] four-port ccw weights {np.round(w4, 12)}: FFT bin 0 {abs(f4[0]):.3e}, bin 1 "
              f"{abs(f4[1]):.3e}, bin 15 {abs(f4[15]):.6f} (ASSERTED 0, 0, 4; the k = -1 alias "
              f"family is bins 3 / 7 / 11 / 15)")
        print(f"[W0] two-port ccw weights {np.round(w2, 12)}: FFT bin 0 {abs(f2[0]):.6f} "
              f"(ASSERTED sqrt(2) = {np.sqrt(2.0):.6f}; the k = 0 term the two-port drive feeds "
              f"and the four-port drive does not)", flush=True)
    assert abs(f4[0]) <= 1e-12 and abs(f4[1]) <= 1e-12 and abs(abs(f4[15]) - 4.0) <= 1e-12
    assert abs(abs(f2[0]) - np.sqrt(2.0)) <= 1e-12 and abs(f2[0]) > 1e-12  # negative control

    # ---- the stored node record, read BEFORE the builder (it rewrites the file) -----
    assert os.path.exists(NODE_JSON), f"{NODE_JSON} absent: stop (never re-solve)"
    node_rec = None
    if rank0:
        with open(NODE_JSON) as fh:
            node_rec = json.load(fh)
    node_rec = comm.bcast(node_rec, root=0)
    assert abs(node_rec["c_in_pf"] - C_IN_F_HUMAN_PF) == 0.0, node_rec["c_in_pf"]
    s_rec = np.array([[complex(a, b) for a, b in row] for row in node_rec["s_node"]],
                     dtype=np.complex128)

    # ---- one build serves both drives ---------------------------------------------
    out = build_tuned_quadrature(comm, with_z50=False)
    omega, ctx, times = out["omega"], out["ctx"], out["times"]
    built = out["built"]
    msh, cell_tags, sigma = ctx["msh"], ctx["cell_tags"], out["sigma"]
    tags_f = ctx["tags_f"]
    order, slots, feeds2 = out["order"], out["slots"], out["feeds"]
    parts = out["parts"]
    centre, radius, half_height, _lo, _hi = _phantom_geometry(msh, cell_tags, PHANTOM_CELL_TAG, comm)
    z_c = float(centre[2])

    # ---- the four feeds: slots 0 / 4 / 8 / 12 on the feed ring ---------------------
    ref_k, ref_sign = slots[order[0]]
    feeds4 = []
    for k in (0, 4, 8, 12):
        match = [pid for pid, (kk, sg) in slots.items() if sg == ref_sign and kk == (ref_k + k) % 16]
        assert len(match) == 1, (k, match)
        feeds4.append(match[0])
    az_by = {f"P{s['ordinal']}": float(s["azimuth_deg"]) for s in built["sheets"]}
    az4 = [az_by[p] for p in feeds4]
    d_az4 = [((az4[(i + 1) % 4] - az4[i]) % 360.0) for i in range(4)]
    assert all(abs(d - 90.0) < AZIMUTH_MATCH_DEG for d in d_az4), d_az4
    assert feeds4[0] == feeds2[0] and feeds4[1] == feeds2[1], (feeds4, feeds2)

    # EX-57 setup figure -- right after the fixture returns the built mesh.
    names = {2: "air", PHANTOM_CELL_TAG: "phantom", CONDUCTOR_CELL_TAG: "conductor (hole route: void)"}
    for sp in ctx["specs"]:
        n = int(sp.port_id[1:])
        if sp.port_id in feeds4:
            label = f"feed {sp.port_id} ({int(round(90 * feeds4.index(sp.port_id)))} deg)"
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
        slice_normal=(0.0, 0.0, 1.0),
        slice_origin=(0.0, 0.0, z_c),
        title="mri:18 -- tuned copper F-human birdcage, four feeds 0 / 90 / 180 / 270 deg at "
              "C_in 118.2437 pF, r = 0.2 Ohm (WF-10 H9')",
        window_size=(1400, 630),
        comm=comm,
    )

    # ---- (A0), (S_rep) -------------------------------------------------------------
    cell_rel = out["cells"] / HOLE_CELL_RECORD - 1.0
    s_node = np.asarray(out["s_node"], dtype=np.complex128)
    r_srep = _residual(s_node, s_rec)
    if rank0:
        print(f"\n[anchor] (A0) cells {out['cells']} vs record {HOLE_CELL_RECORD}, rel "
              f"{cell_rel:.3e} (band {CELL_COUNT_BAND}); node C_in {C_IN_F_HUMAN_PF} pF "
              f"(secant {C_IN_SECANT_PF:.6f}); two-port feeds {feeds2}, az "
              f"{out['az_feeds'][0]:.6f} / {out['az_feeds'][1]:.6f} deg, mirror plane "
              f"{out['az_mirror']:.6f} deg  (W1: 369387, 0; 11.25 / 101.25; 56.25)")
        print(f"[anchor] four-port feeds {feeds4}, az "
              + ", ".join(f"{a:.6f}" for a in az4) + f" deg, successive spacings "
              + ", ".join(f"{d:.6f}" for d in d_az4) + f" (90 within {AZIMUTH_MATCH_DEG:g}, ASSERTED)")
        print(f"[anchor] (S_rep) s_node vs the stored node JSON: {r_srep:.6e} <= "
              f"{REUSE_REPRODUCTION_RTOL:g}  (W1_SAR :6312: 6.027e-13)", flush=True)

    # ---- two-port: Y_in and the signed k spectrum ------------------------------------
    y_rows = {"cw": ("Y17_cw", "Y21_cw"), "ccw": ("Y17_ccw", "Y21_ccw")}
    if rank0:
        print("\n[Y_in] (two-port readings)")
    for sense in ("cw", "ccw"):
        for i, pid in enumerate(feeds2):
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
    k2 = {}
    for sense in ("ccw", "cw"):
        sh = _signed_k_shares(out["cur"][sense], order, slots)
        k2[sense] = sh
        spec_f = _spectrum_share(_ring_vector(out["cur"][sense], order, slots, ref_sign))
        un = {s: sh[s]["lead"] + sh[s]["lag"] for s in (1.0, -1.0)}
        got = (f"{sh[1.0]['lead']:.6f}", f"{sh[1.0]['lag']:.6f}", f"{sh[-1.0]['lead']:.6f}",
               f"{sh[-1.0]['lag']:.6f}", f"{un[1.0]:.6f}", f"{un[-1.0]:.6f}")
        ok = chk(f"k_{sense}", got, REC[f"k_{sense}"])
        if rank0:
            print(f"[k] [{sense}] ring+ lead/lag {got[0]}/{got[1]}, ring- lead/lag {got[2]}/{got[3]}, "
                  f"unsigned ring+ {got[4]} ring- {got[5]}   {'= W1' if ok else 'DIFFERS'}; "
                  f"feed-ring k = 0 share {spec_f[0]:.3e}")

    # ---- two-port B, CG1 projection, point reads (358's sequence) ----------------------
    comm.Barrier()
    b_dg, cg1 = {}, {}
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
        print(f"[B1+] (Sense) mean|B1+|_ccw / mean|B1-|_ccw = {sense_ratio:.4f} (> 1 ASSERTED; "
              f"PREDICTED >= {PREDICTED_SENSE_RATIO:g}; W1 {REC['sense']}); cw: "
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
        print("[B1+] two-port C16 readings (printed, never asserted): "
              + ", ".join(f"m={m}: {v:.4e}" for m, v in c16.items())
              + f"; worst m={worst_m} {c16[worst_m]:.4e} (PREDICTED <= {PREDICTED_C16:g}; W1 "
              f"m=8 {REC['c16_worst'][1]})")
        print(f"[B1+] two-port mean|B1+|_ccw {mean_b:.6e} T, CV {cv_b:.4e} (reading, no homogeneity "
              f"claim); /sqrt(P_src) {mean_b / np.sqrt(p_src_ccw):.6e} T/sqrt(W)  (W1 "
              f"{REC['mean_b']}, {REC['cv_b']}, {REC['b_per_sqrt_w']})", flush=True)

    # ---- (Sigma_P) on the two-port field --------------------------------------------
    sar = mean_sar(out["e"]["ccw"], sigma=sigma, rho=PHANTOM_RHO_KG_PER_M3,
                   cell_tags=cell_tags, subdomain_ids=PHANTOM_CELL_TAG, comm=comm)
    p_ph = parts["ccw"]["p_phantom"]
    r_sigma = abs(sar["dissipated_power_w"] - p_ph) / abs(p_ph)
    chk("sar_mean", f"{sar['mean_sar_w_per_kg']:.6e}", REC["sar_mean"])
    chk("sar_mean_per_w", f"{sar['mean_sar_w_per_kg'] / p_src_ccw:.6e}", REC["sar_mean_per_w"])
    chk("pph_over_psrc", f"{p_ph / p_src_ccw:.6e}", REC["pph_over_psrc"])
    if rank0:
        print(f"\n[SAR] (Sigma_P) P_phantom {p_ph:.12e} W vs mean_sar dissipated "
              f"{sar['dissipated_power_w']:.12e} W, rel {r_sigma:.3e} (<= {EXACT_IDENTITY_RTOL:g} "
              f"ASSERTED; W1 0.000e+00); whole-phantom mean {sar['mean_sar_w_per_kg']:.6e} W/kg, "
              f"per W of P_src {sar['mean_sar_w_per_kg'] / p_src_ccw:.6e}; P_phantom/P_src "
              f"{p_ph / p_src_ccw:.6e}  (W1 {REC['sar_mean']}, {REC['sar_mean_per_w']}, "
              f"{REC['pph_over_psrc']})", flush=True)

    # ---- the two-port anchors, asserted at their imported bands -----------------------
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
              f"builder; W1 7.830e-13, 7012.5x)", flush=True)
    assert abs(cell_rel) <= CELL_COUNT_BAND and out["cells"] == HOLE_CELL_RECORD, out["cells"]
    assert r_srep <= REUSE_REPRODUCTION_RTOL, f"(S_rep) {r_srep:.6e}"
    assert r_sigma <= EXACT_IDENTITY_RTOL, f"(Sigma_P) {r_sigma:.3e}"
    assert n_mask > 0, "no sample point valid on every image"
    assert sense_ratio > 1.0, f"(Sense) red: B1- dominant under ccw ({sense_ratio:.4f})"
    assert mirror <= C4_COVARIANCE_BAND, f"(M_B) mirror {mirror:.4e}"
    assert mispaired > C4_COVARIANCE_BAND, f"(N1) mis-paired control inside the band {mispaired:.4e}"

    # ---- keep only what the maps need; free the rest before the four-port solves -------
    bp2 = b1_plus(b_dg["ccw"])
    bm2 = b1_minus(b_dg["ccw"])
    dg0 = fem.functionspace(msh, ("DG", 0))

    def real_dg0(fn, name):
        f = fem.Function(fn.function_space, name=name)
        f.x.array[:] = np.real(fn.x.array)
        f.x.scatter_forward()
        return f

    bplus_2 = real_dg0(bp2, "B1_plus_2port")
    bminus_2 = real_dg0(bm2, "B1_minus_2port")
    e_ccw2 = out["e"]["ccw"]
    e_sq2 = fem.Function(dg0, name="e_sq_2port")
    e_sq2.interpolate(fem.Expression(ufl.inner(e_ccw2, e_ccw2), dg0.element.interpolation_points))
    emag_2 = np.sqrt(np.maximum(np.real(e_sq2.x.array), 0.0))
    del bp2, bm2, e_sq2, cg1, b_dg
    del out["e"]["cw"]
    del out["e17"], out["e21"]
    cur2_ccw = out["cur"]["ccw"]
    p_src2 = float(p_src_ccw)

    # ---- the four-port drive: 358's _drive sequence at the base commit ------------------
    # (gate module lines 346-413): wall first in the factorising call; the first solve
    # factorises (never hold a factor across an impedance change, PORT-19), the next
    # three are held.
    z0 = float(REFERENCE_IMPEDANCE_OHM)
    c_f = C_IN_F_HUMAN_PF * 1e-12
    z_rc = _z_rc(c_f, FREQUENCY_HZ)
    z_th = _thevenin_feed_impedance(z_rc, z0)
    v_s = _source_voltage(z_rc, z0)
    z_s = (1.0 + 1.0j) * surface_resistance_ohm(omega, COPPER_SIGMA)
    wall = _leontovich_term(msh, built["facets"], z_s, omega)
    solver = ctx["solver"]
    base_specs = list(ctx["specs"])
    imp4 = {pid: z_rc for pid in order}
    imp4.update({pid: z_th for pid in feeds4})
    specs4 = _respec(base_specs, None, imp4)
    assert all(sp.width_correction_kappa is None for sp in specs4)
    specs4_by = {sp.port_id: sp for sp in specs4}
    acct = {"mesh": msh, "facet_tags": tags_f}
    if rank0:
        print(f"\n[four-port] node C = {c_f:.12e} F: Z_RC {z_rc:.9e} Ohm, Z_th {z_th:.9e} Ohm, "
              f"V_s {v_s:.9e} V; feeds {feeds4} at Z_th, the other 28 gaps at Z_RC", flush=True)

    e4, cur4, kinds4, t4 = [], [], [], []
    for j, pid_j in enumerate(feeds4):
        held = j > 0
        sheets = [sp.sheet(driven=(sp.port_id == pid_j)) for sp in specs4]
        dsheet = next(s for s in sheets if s.port_id == pid_j)
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
        t4.append(time.perf_counter() - t0)
        kinds4.append("held" if held else "solve")
        e = fields.e_complex
        currents = {s.port_id: sheet_terminal_current(msh, tags_f, s, e, comm) for s in sheets}
        cur4.append(np.array([currents[p] for p in order], dtype=np.complex128))
        cp = fem.Function(e.function_space, name=f"E_{pid_j}_4port")
        cp.x.array[:] = e.x.array
        cp.x.scatter_forward()
        e4.append(cp)
        del fields, e, sheets, dsheet, linear
    times["four_port"] = float(sum(t4))
    if rank0:
        print(f"[four-port] solves {kinds4}, " + ", ".join(f"{t:.2f}" for t in t4) + " s", flush=True)

    # ---- the superposition, (X4-P), (X4-L) --------------------------------------------
    ef4 = fem.Function(e4[0].function_space, name="E_4port")
    ef4.x.array[:] = sum(complex(w4[j]) * e4[j].x.array for j in range(4))
    ef4.x.scatter_forward()
    del e4
    cur4_sup = sum(complex(w4[j]) * cur4[j] for j in range(4))
    wf = {pid: complex(w4[j]) for j, pid in enumerate(feeds4)}
    sheets_w = {pid: (dataclasses.replace(specs4_by[pid].sheet(driven=True),
                                          source_voltage_v=wf[pid])
                      if pid in wf else specs4_by[pid].sheet(driven=False))
                for pid in order}
    cur_read = np.array([sheet_terminal_current(msh, tags_f, sheets_w[pid], ef4, comm)
                         for pid in order], dtype=np.complex128)
    r_x4l = float(np.max(np.abs(cur_read - cur4_sup)) / np.max(np.abs(cur4_sup)))
    dx = ufl.Measure("dx", domain=msh, subdomain_data=cell_tags)

    def _vol(e, measure):
        form = fem.form(0.5 * sigma * ufl.inner(e, e) * measure)
        return float(np.real(comm.allreduce(fem.assemble_scalar(form), op=MPI.SUM)))

    p4_src = float(sum(_source_power_w(acct, sheets_w[pid], ef4, omega) for pid in feeds4))
    p4_sheet_by = {pid: float(_sheet_field_dissipation_w(acct, sheets_w[pid], ef4, omega))
                   for pid in order}
    p4_sheets = float(sum(p4_sheet_by.values()))
    p4_total = _vol(ef4, dx)
    p4_phantom = _vol(ef4, dx(PHANTOM_CELL_TAG))
    p4_rest = p4_total - p4_phantom
    p4_surf = _surface_loss_w(msh, built["facets"], z_s, ef4, comm)
    r_x4p = abs(p4_src - p4_sheets - p4_total - p4_surf) / abs(p4_src)
    y4 = [complex(1.0 / z_rc + cur_read[order.index(pid)]
                  / (wf[pid] - z_th * cur_read[order.index(pid)])) for pid in feeds4]
    if rank0:
        print(f"[four-port] (X4-P) P_src {p4_src:.9e} sumP_sheet {p4_sheets:.9e} P_phantom "
              f"{p4_phantom:.9e} P_rest {p4_rest:.9e} P_surf {p4_surf:.9e} W; rel {r_x4p:.3e} "
              f"(<= {DISCRETE_IDENTITY_RTOL:g} ASSERTED, example-side; PREDICTED <= 1e-12)")
        print(f"[four-port] (X4-L) 32 terminal currents off the superposed field vs sum w_i I_i: "
              f"{r_x4l:.3e} (<= {UNIT_WEIGHT_CONTROL_RTOL:g} ASSERTED, example-side; PREDICTED "
              f"<= 1e-13)")
        print(f"[four-port] (X4-K) kinds {kinds4} (ASSERTED ['solve', 'held', 'held', 'held'])")
        print("[four-port] Y_in at the four feeds, untuned at the two-port root C (printed): "
              + "; ".join(f"{pid} {y.real:+.6e} {y.imag:+.6e}j S" for pid, y in zip(feeds4, y4)),
              flush=True)
    assert kinds4 == ["solve", "held", "held", "held"], kinds4
    assert r_x4p <= DISCRETE_IDENTITY_RTOL, f"(X4-P) {r_x4p:.3e}"
    assert r_x4l <= UNIT_WEIGHT_CONTROL_RTOL, f"(X4-L) {r_x4l:.3e}"

    # ---- k spectra, both drives (printed) ------------------------------------------------
    k4 = _signed_k_shares(cur4_sup, order, slots)
    spec2 = {s: _spectrum_share(_ring_vector(cur2_ccw, order, slots, s)) for s in (1.0, -1.0)}
    spec4 = {s: _spectrum_share(_ring_vector(cur4_sup, order, slots, s)) for s in (1.0, -1.0)}
    csv_k = OUTPUT_DIR / f"{BASENAME}_kspectrum.csv"
    if rank0:
        OUTPUT_DIR.mkdir(exist_ok=True)
        with open(csv_k, "w") as fh:
            fh.write("drive,ring," + ",".join(f"bin{b}" for b in range(16))
                     + ",k0_share,lead_k1_share,lag_k15_share\n")
            for name, spec in (("2port_ccw", spec2), ("4port_ccw", spec4)):
                for s in (1.0, -1.0):
                    row = spec[s]
                    fh.write(f"{name},{s:+.0f}," + ",".join(f"{v:.9e}" for v in row)
                             + f",{row[0]:.9e},{row[1]:.9e},{row[15]:.9e}\n")
        k2c = k2["ccw"]
        print(f"\n[k] feed ring (sign {ref_sign:+.0f}) k-spectrum, power shares (printed):")
        print(f"   two-port : k = 0 {spec2[ref_sign][0]:.4e}, k = +1 (lead) {spec2[ref_sign][1]:.6f}, "
              f"k = -1 (lag, bin 15) {spec2[ref_sign][15]:.6f}  (W1 lag ring+/ring- "
              f"{REC['k_ccw'][1]} / {REC['k_ccw'][3]}; {k2c[ref_sign]['lag']:.6f})")
        print(f"   four-port: k = 0 {spec4[ref_sign][0]:.4e} (PREDICTED <= 1e-4), k = +1 (lead) "
              f"{spec4[ref_sign][1]:.6f}, k = -1 (lag, bin 15) {spec4[ref_sign][15]:.6f} "
              f"(PREDICTED >= two-port's {REC['k_ccw'][1]} / {REC['k_ccw'][3]}); other-ring "
              f"k = 0 {spec4[-ref_sign][0]:.4e}, lag {spec4[-ref_sign][15]:.6f}; unsigned shares "
              f"ring+ {k4[1.0]['lead'] + k4[1.0]['lag']:.6f}, ring- "
              f"{k4[-1.0]['lead'] + k4[-1.0]['lag']:.6f}")
        print("   why no k = 0: the four weights sum to zero (W0 bin 0) and the loading is "
              "C4-symmetric, so the 16-slot drive has only bins 3 / 7 / 11 / 15 (k = -1 family) "
              "plus mesh-asymmetry leakage.", flush=True)

    # ---- four-port B, CG1, point reads ----------------------------------------------------
    comm.Barrier()
    b4_dg = magnetic_flux_density_from_e(ef4, omega, name="B_4port")
    b4_cg1 = project_to_cg1(b4_dg, name="B_4port_cg1")
    comm.Barrier()
    x4r = _read_senses(b4_cg1, points)
    rot4 = {m: _read_senses(b4_cg1, _rotate_z(points, np.radians(m * 22.5))) for m in range(1, 16)}
    mask4 = np.logical_and.reduce([x4r[2]] + [r[2] for r in rot4.values()])
    plus4 = x4r[0]
    c16_4 = {m: _relative_l2(rot4[m][0], plus4, mask4) for m in rot4}
    worst4 = max(c16_4, key=c16_4.get)
    c4_spread = max(c16_4[m] for m in (4, 8, 12))
    mask_c = np.logical_and(mask, mask4)
    mean4 = float(np.mean(plus4[mask_c]))
    mean2c = float(np.mean(plus_x[mask_c]))
    bps2 = mean2c / np.sqrt(p_src2)
    bps4 = mean4 / np.sqrt(p4_src)
    sense4 = float(np.mean(x4r[0][mask4]) / np.mean(x4r[1][mask4]))
    if rank0:
        print(f"\n[B1+] four-port: {int(np.count_nonzero(mask4))} valid points; (Sense) mean|B1+|/"
              f"mean|B1-| {sense4:.4f} (printed); mean|B1+|_ccw {mean4:.6e} T (common mask "
              f"{int(np.count_nonzero(mask_c))} points)")
        print(f"[B1+] four-port C4 spread of |B1+|_ccw at 90 deg rotations (m = 4 / 8 / 12): "
              + " / ".join(f"{c16_4[m]:.4e}" for m in (4, 8, 12))
              + f"; max {c4_spread:.4e} (PREDICTED <= 1e-2, printed)")
        print("[B1+] four-port C16 readings (printed): "
              + ", ".join(f"m={m}: {v:.4e}" for m, v in c16_4.items())
              + f"; worst m={worst4} {c16_4[worst4]:.4e}; two-port worst m={worst_m} "
              f"{c16[worst_m]:.4e}  (both PREDICTED worst at a non-multiple of the symmetry)")
        print(f"[B1+] |B1+|_ccw/sqrt(P_src) on the common mask: two-port {bps2:.6e}, four-port "
              f"{bps4:.6e} T/sqrt(W); four-port / two-port {bps4 / bps2:.4f} (PREDICTED 1 +- 0.1, "
              f"printed); P_src two-port {p_src2:.9e} W, four-port {p4_src:.9e} W", flush=True)

    # ---- maps on one colour scale ---------------------------------------------------------
    def named(name, values):
        f = fem.Function(dg0, name=name)
        f.x.array[:] = values
        f.x.scatter_forward()
        return f

    bplus_4 = real_dg0(b1_plus(b4_dg), "B1_plus_4port")
    bminus_4 = real_dg0(b1_minus(b4_dg), "B1_minus_4port")
    e_sq4 = fem.Function(dg0, name="e_sq_4port")
    e_sq4.interpolate(fem.Expression(ufl.inner(ef4, ef4), dg0.element.interpolation_points))
    emag_4 = np.sqrt(np.maximum(np.real(e_sq4.x.array), 0.0))
    n_loc = dg0.dofmap.index_map.size_local
    lo_l = float(min(np.min(bplus_2.x.array[:n_loc]), np.min(bplus_4.x.array[:n_loc])))
    hi_l = float(max(np.max(bplus_2.x.array[:n_loc]), np.max(bplus_4.x.array[:n_loc])))
    lo = float(comm.allreduce(lo_l, op=MPI.MIN))
    hi = float(comm.allreduce(hi_l, op=MPI.MAX))

    tdim = msh.topology.dim
    imap = msh.topology.index_map(tdim)
    n_cells = imap.size_local + imap.num_ghosts
    cell_dof = np.asarray(dg0.dofmap.list).reshape(-1)[:n_cells]
    marker4 = np.zeros_like(bplus_2.x.array)
    marker2 = np.zeros_like(bplus_2.x.array)
    ct_idx = np.asarray(cell_tags.indices)
    ct_val = np.asarray(cell_tags.values)
    for j, pid in enumerate(feeds4):
        n = int(pid[1:])
        sel = ct_idx[(ct_val == PORT_LOWER + n) | (ct_val == PORT_UPPER + n)]
        sel = sel[sel < n_cells]
        marker4[cell_dof[sel]] = float(j + 1)
        if j < 2:
            marker2[cell_dof[sel]] = float(j + 1)

    export = {
        "B1_plus_2port": bplus_2,
        "B1_minus_2port": bminus_2,
        "B1_plus_4port": bplus_4,
        "B1_minus_4port": bminus_4,
        "E_magnitude_2port": named("E_magnitude_2port", emag_2),
        "E_magnitude_4port": named("E_magnitude_4port", emag_4),
        "FeedMarker": named("FeedMarker", marker4),
        "FeedMarker2": named("FeedMarker2", marker2),
    }
    OUTPUT_DIR.mkdir(exist_ok=True)
    xdmf_path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_combined", msh, cell_tags, export, comm=comm
    )
    adopt_host_ownership(OUTPUT_DIR, comm=comm)
    _m, _rss_per, rss_sum = _rss_gib(comm)

    if rank0:
        print(f"\n[records] {len(mism)} printed two-port record row(s) differ from the W1 digits"
              + ("" if not mism else ":"))
        for label, got, want in mism:
            print(f"   DIFFERS {label}: got {got}, W1 {want}")
        print(f"[price] build {times['build']:.1f} s, node pair {times['node_pair']:.1f} s, "
              f"superposition {times['superposition']:.1f} s, four-port solves "
              f"{times['four_port']:.1f} s; RSS summed {rss_sum:.2f} GiB; -n {comm.size}")
        print(f"\n[output] {xdmf_path.name}: DG0 " + ", ".join(export) + ", CellTags")
        print(f"[output] {csv_k.name}")
        print(f"[paraview] ONE COLOUR SCALE for B1_plus_2port and B1_plus_4port: "
              f"[{lo:.6e}, {hi:.6e}] T (allreduced min / max over owned cells). Slice z = z_c; "
              "compare B1_plus_2port with B1_plus_4port; threshold FeedMarker >= 1 to see the feeds.")
        print(f"\nEX-75 two-port anchors re-asserted, four-port construction identities and both "
              f"negative controls: PASS. Total elapsed {time.perf_counter() - started:.1f} s.",
              flush=True)


if __name__ == "__main__":
    main()
