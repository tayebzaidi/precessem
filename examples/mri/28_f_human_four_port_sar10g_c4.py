"""Example (`EX-127`): the 10 g lattice under the four-port drive on the tuned
copper F-human birdcage -- C4 holds at 10 g, the clinical two-port misses it,
and the spread per radius reads drive anisotropy (`WF-10` step 5).

`WF-10` step 5 (`tests/validation/test_wf10_f_human_four_port_sar10g.py`, item
422, commit ``a11ef456``) drives the copper hole mesh (369 387 cells, Leontovich
wall, kappa off, degree 1, 64 MHz, r = 0.2 Ohm on all 32 ring gaps) at the
two-port root C_in = 118.2437 pF with `EX-75`'s four feeds (P17 / P21 / P25 /
P29, slots 0 / 4 / 8 / 12, ccw quadrature phases), then reads the 18-point 10 g
lattice under that drive and the 12 z_c points under the two-port ccw drive.

Every record, band and helper is imported from the gate module (the ANS-1
rule).  The gate's ``test_*`` bodies are never imported: this script repeats
the body's sequence (provenance: ``test_wf10_step5_four_port_sar10g``, lines
342-565 at ``a11ef456``) through ``_four_port_drive``, ``_lattice``,
``_zc_rows``, ``_c4_misses``, ``_spreads`` and ``_containment``.  The printed
records are compared with the gate run of record
(``20261013T141801Z_WF-10-step5-W1.log``) in the example-local ``REC`` table at
a relative tolerance of 10**-(d-2), d the printed significant digits; a record
off its tolerance stops the run.  Every asserted band is an imported one.

Example-side identities (asserted, exact): (O) the 12 z_c balls are one C4
orbit per radius; (X1)_4 / (X1)_2 the per-cell DG0 phantom loss maps sum to the
phantom power.  Negative control (asserted): the two-port ccw table misses C4
by more than ``C4_COVARIANCE_BAND``.  Printed only: the 180 degree blindness
(m = 2 alone passes both drives).

Scope: degree 1, 64 MHz, kappa off, one copper hole fixture, the 0.015 m
phantom; the four-port drive runs UNTUNED at the two-port root (as in
`EX-75`): per-W values compare shapes, not tuned-coil SARs.  No tuning, f0, Q,
match, S11 / Z_in, absolute SAR, compliance, h-statement, degree 2, AED or
external claim.

Run through the example runner (the ``mri:`` group sources the complex build)::

    ./scripts/run_examples.sh -e mri:28 -n 8 -t 590

Output (gitignored): ``examples/mri/paraview_output/mri_28_f_human_four_port_sar10g_c4_combined.xdmf``
and ``mri_28_f_human_four_port_sar10g_c4_zc_rows.csv``.
"""

from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import ufl
from dolfinx import default_scalar_type, fem, la
from dolfinx.mesh import compute_midpoints
from mpi4py import MPI

from fem_em_solver.io.paraview_utils import adopt_host_ownership, write_xdmf_with_tags
from fem_em_solver.post.setup_figure import write_setup_figure

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

import tests.validation.test_birdcage_b1_plus_map as _b1map  # noqa: E402
import tests.validation.test_port_birdcage_ring_matrix as _ringmatrix  # noqa: E402
import tests.validation.test_wf10_f_human_four_port_sar10g as g5  # noqa: E402

# ---- the gate body's own import lines (``:308-337``), copied verbatim --------------
from fem_em_solver.post.sar import (  # noqa: E402
    averaging_ball_radius,
    build_density_field,
    mass_averaged_sar,
    mean_sar,
)
from tests.mesh.helpers import global_cell_tag_set  # noqa: E402
from tests.validation.test_birdcage_b1_plus_map import PHANTOM_RHO_KG_PER_M3  # noqa: E402
from tests.validation.test_birdcage_f_human_rung import F_HUMAN_RING_RADIUS  # noqa: E402
from tests.validation.test_birdcage_power_identity import DISCRETE_IDENTITY_RTOL  # noqa: E402
from tests.validation.test_birdcage_sar_mass_averaged import EXACT_IDENTITY_RTOL  # noqa: E402
from tests.validation.test_mass_averaged_sar_standard_masses import (  # noqa: E402
    QUADRATURE_DEGREE,
    TEN_GRAM_KG,
)
from tests.validation.test_port_birdcage_lumped_column import PHANTOM_CELL_TAG  # noqa: E402
from tests.validation.test_port_birdcage_ring_column import AZIMUTH_MATCH_DEG  # noqa: E402
from tests.validation.test_port_birdcage_ring_matrix import REUSE_REPRODUCTION_RTOL  # noqa: E402
from tests.validation.test_port_circuit_layer_field import _residual  # noqa: E402
from tests.validation.test_port_drive_superposition import UNIT_WEIGHT_CONTROL_RTOL  # noqa: E402
from tests.validation.test_wf7_f_human_sar10g_hotspot import (  # noqa: E402
    CONDUCTOR_CLEARANCE_M,
    COVERAGE_BALL_PAD_M,
    FIRST_BALL_STOP_S,
    _phantom_geometry,
)
from tests.validation.test_wf10_f_human_copper_esr_root import Y0_CODING_BAND  # noqa: E402
from tests.validation.test_wf10_f_human_tuned_b1_sar import build_tuned_quadrature  # noqa: E402
from tests.validation.test_wf12_f_human_hole import _rss_gib  # noqa: E402

# ---- as ``mri:18`` imports them --------------------------------------------------------
from tests.mesh.test_birdcage_port_sheets import PORT_LOWER, PORT_UPPER  # noqa: E402
from tests.validation.test_port_birdcage_lumped_column import CONDUCTOR_CELL_TAG  # noqa: E402

# ---- names taken from the gate module's namespace (never restated) --------------------
_four_port_drive = g5._four_port_drive
_c4_misses = g5._c4_misses
_spreads = g5._spreads
_zc_rows = g5._zc_rows
_spectrum_share = g5._spectrum_share
_ring_vector = g5._ring_vector
_rotate_z = g5._rotate_z
_lattice = g5._lattice
_containment = g5._containment
_degree_knob = g5._degree_knob
_signed_k_shares = g5._signed_k_shares
TWO_PORT_ZC_RECORD = g5.TWO_PORT_ZC_RECORD
TWO_PORT_SPREADS = g5.TWO_PORT_SPREADS
TWO_PORT_C4_MISSES = g5.TWO_PORT_C4_MISSES
VERDICT_FRACTION = g5.VERDICT_FRACTION
PREDICTED_C4 = g5.PREDICTED_C4
PREDICTED_SPREAD4 = g5.PREDICTED_SPREAD4
PREDICTED_TWO_PORT_REPRO = g5.PREDICTED_TWO_PORT_REPRO
C4_COVARIANCE_BAND = g5.C4_COVARIANCE_BAND
LATTICE_R_M = g5.LATTICE_R_M
C_IN_F_HUMAN_PF = g5.C_IN_F_HUMAN_PF
NODE_JSON = g5.NODE_JSON

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"  # committed, EX-57
BASENAME = "mri_28_f_human_four_port_sar10g_c4"
W1 = "20261013T141801Z_WF-10-step5-W1.log"

# ---- the gate run of record, restated as EX-89 item 238's REC was: each entry
# commented with its W1 line; compared at 10**-(d-2), d the printed digits. -----------
REC = {
    "cells": 369387,                                  # :6324 (exact)
    # :6376-6393, "%.9e" W/kg, d = 10 -> 1e-8
    "sar4": {
        "r=15 dz=+0 az17+0": "9.391419740e-04",       # :6376
        "r=30 dz=+0 az17+0": "3.020877141e-03",
        "r=45 dz=+0 az17+0": "6.546455227e-03",
        "r=15 dz=+0 az17+90": "9.112890186e-04",
        "r=30 dz=+0 az17+90": "2.898964058e-03",
        "r=45 dz=+0 az17+90": "6.504834035e-03",
        "r=15 dz=+0 az17+180": "9.095721163e-04",
        "r=30 dz=+0 az17+180": "3.026573257e-03",
        "r=45 dz=+0 az17+180": "6.503503968e-03",
        "r=15 dz=+0 az17+270": "9.221820473e-04",
        "r=30 dz=+0 az17+270": "3.016641875e-03",
        "r=45 dz=+0 az17+270": "6.554475109e-03",     # :6387
        "r=15 dz=-50 az17+0": "2.333015501e-03",      # :6388
        "r=30 dz=-50 az17+0": "3.099452879e-03",
        "r=45 dz=-50 az17+0": "4.626128485e-03",
        "r=15 dz=+50 az17+0": "2.315224069e-03",
        "r=30 dz=+50 az17+0": "3.017352834e-03",
        "r=45 dz=+50 az17+0": "4.530598061e-03",      # :6393
    },
    # ball masses, g, "%.6f" (7 significant digits) -> 1e-5; :6339-6368 (the 4p rows)
    "mass_g": {
        "r=15 dz=+0 az17+0": "9.994606", "r=30 dz=+0 az17+0": "10.028757",
        "r=45 dz=+0 az17+0": "10.010236", "r=15 dz=+0 az17+90": "10.050608",
        "r=30 dz=+0 az17+90": "10.098068", "r=45 dz=+0 az17+90": "10.009283",
        "r=15 dz=+0 az17+180": "10.024962", "r=30 dz=+0 az17+180": "9.945828",
        "r=45 dz=+0 az17+180": "10.032923", "r=15 dz=+0 az17+270": "9.976829",
        "r=30 dz=+0 az17+270": "9.926530", "r=45 dz=+0 az17+270": "10.032609",
        "r=15 dz=-50 az17+0": "10.004065", "r=30 dz=-50 az17+0": "9.991901",
        "r=45 dz=-50 az17+0": "9.976512", "r=15 dz=+50 az17+0": "9.986379",
        "r=30 dz=+50 az17+0": "10.064920", "r=45 dz=+50 az17+0": "9.994204",
    },
    # :6395 (C4_4) per m, percent: (string, d)
    "c4_4": {1: ("1.8127", 5), 2: ("0.7215", 4), 3: ("0.2645", 4)},
    # :6396 control per m, percent
    "c4_2": {1: ("17.3831", 6), 2: ("0.8080", 4), 3: ("16.1396", 6)},
    # :6398 spread per radius, percent: four-port (5 / 5 / 4), two-port in-window (6)
    "sp4": (("3.2122", 5), ("4.2668", 5), ("0.7809", 4)),
    "sp2": (("16.7267", 6), ("20.9986", 6), ("18.6041", 6)),
    "verdict": "drive anisotropy",                    # :6399
    "below": [True, True, True],                      # :6399
    "c_star4": "r=45 dz=+0 az17+270",                 # :6394
    "sar_star4": "6.554475109e-03",                   # :6394, d = 10
    "p_src4": "2.102028870e-01",                      # :6371, d = 10
    "p_phantom4": "1.204416570e-02",                  # :6371, d = 10
    "k0": ("4.5411e-07", 5),                          # :6374
    "k_bin15": ("0.999450", 6),                       # :6374
    "per_watt": ("3.118166e-02", 7),                  # :6404
    "peak_to_mean": ("1.2031", 5),                    # :6404
    "pph_over_psrc": ("5.729781e-02", 7),             # :6404
}
# Excluded from record equality (round-off readings, printed at their bands only):
# (S_rep) :6337, (X4-P) :6371, (X4-L) :6372, (Sigma_P)4 :6401, (i)4 :6400, the
# two-port (P_w) / (Lin) :6402, the in-window-vs-359 reproduction :6397.


def _dg0_vec(form):
    """Assemble a DG0-test-function linear form to a per-cell real vector
    (copied from ``mri:17``, ``examples/mri/17_f_human_tuned_quadrature_b1_sar.py``)."""
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

    def P(msg: str = "") -> None:
        if rank0:
            print(msg, flush=True)

    P("=" * 78)
    P("EX-127 -- the 10 g lattice under the four-port drive, tuned copper F-human "
      "(WF-10 step 5)")
    P("=" * 78)
    P("degree 1, untuned four-port at the two-port root: per-W SAR values are not an "
      "absolute or tuned-coil SAR")
    P(f"[records] every record row is compared with {W1} (REC table) at 10**-(d-2)")

    # ---- structural asserts ---------------------------------------------------------
    assert comm.size == 8, f"this window is -n 8 (got {comm.size})"
    assert _degree_knob() is None, "FEM_EM_WF10_DEGREE must be unset (degree 1 only)"
    assert os.path.exists(NODE_JSON), "NODE_JSON absent: stop (BLOCKED, never a re-solve)"
    assert C4_COVARIANCE_BAND is _b1map.C4_COVARIANCE_BAND
    assert REUSE_REPRODUCTION_RTOL is _ringmatrix.REUSE_REPRODUCTION_RTOL

    # ---- NODE_JSON, read BEFORE the builder (which rewrites it) ---------------------
    rec = None
    if rank0:
        with open(NODE_JSON) as fh:
            rec = json.load(fh)
    rec = comm.bcast(rec, root=0)
    assert abs(rec["c_in_pf"] - C_IN_F_HUMAN_PF) == 0.0, rec["c_in_pf"]
    s_rec = np.array([[complex(a, b) for a, b in row] for row in rec["s_node"]],
                     dtype=np.complex128)

    # ---- build + (S_rep) -------------------------------------------------------------
    out = build_tuned_quadrature(comm, with_z50=False)
    ctx, times, omega = out["ctx"], out["times"], out["omega"]
    msh, cell_tags, sigma = ctx["msh"], ctx["cell_tags"], out["sigma"]
    s_node = np.asarray(out["s_node"], dtype=np.complex128)
    r_srep = _residual(s_node, s_rec)

    # ---- geometry ---------------------------------------------------------------------
    centre, radius, half_height, lo, hi = _phantom_geometry(msh, cell_tags, PHANTOM_CELL_TAG, comm)
    z_c = float(centre[2])
    a10 = averaging_ball_radius(mass_kg=TEN_GRAM_KG, rho=PHANTOM_RHO_KG_PER_M3)
    az17 = float(out["az_feeds"][0])
    P(f"[geometry] z_c {z_c:.4e} m; R {radius * 1e3:.4f} mm, H/2 {half_height * 1e3:.4f} mm; "
      f"10 g a = {a10 * 1e3:.4f} mm; az17 {az17:.4f} deg; cells {out['cells']}")
    all_tags = global_cell_tag_set(msh, cell_tags)
    rho_field = build_density_field(
        msh, PHANTOM_RHO_KG_PER_M3, cell_tags=cell_tags,
        density_map={int(t): 0.0 for t in sorted(all_tags) if int(t) != PHANTOM_CELL_TAG})

    # ---- memory: keep only ccw (gate :381-382) ------------------------------------------
    e_ccw = out["e"]["ccw"]
    del out["e"]["cw"], out["e17"], out["e21"]

    # ---- four-port drive + (Sigma_P)4 ----------------------------------------------------
    d4 = _four_port_drive(out, comm)
    ef4 = d4["ef4"]
    tagged4 = mean_sar(ef4, sigma=sigma, rho=PHANTOM_RHO_KG_PER_M3,
                       cell_tags=cell_tags, subdomain_ids=PHANTOM_CELL_TAG, comm=comm)
    p4_ph = d4["parts"]["p_phantom"]
    p4_src = d4["parts"]["p_src"]
    p2_src = out["parts"]["ccw"]["p_src"]
    p2_ph = out["parts"]["ccw"]["p_phantom"]
    r_sigma4 = abs(tagged4["dissipated_power_w"] - p4_ph) / abs(p4_ph)
    az4 = d4["az4"]
    d_az4 = [((az4[(i + 1) % 4] - az4[i]) % 360.0) for i in range(4)]
    az_ok = all(abs(d - 90.0) < AZIMUTH_MATCH_DEG for d in d_az4)
    feeds_ok = d4["feeds4"][0] == d4["feeds2"][0] and d4["feeds4"][1] == d4["feeds2"][1]
    ref_sign = d4["ref_sign"]
    spec4 = {s: _spectrum_share(_ring_vector(d4["cur4"], out["order"], out["slots"], s))
             for s in (1.0, -1.0)}
    k4 = _signed_k_shares(d4["cur4"], out["order"], out["slots"])

    # ---- EX-57 setup figure: after the four feeds are known --------------------------------
    t_fig = time.perf_counter()
    names = {2: "air", PHANTOM_CELL_TAG: "phantom", CONDUCTOR_CELL_TAG: "conductor (hole route: void)"}
    for sp in ctx["specs"]:
        n = int(sp.port_id[1:])
        if sp.port_id in d4["feeds4"]:
            deg = 90 * d4["feeds4"].index(sp.port_id)
            label = f"feed {sp.port_id} ({deg} deg)"
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
        title="mri:28 -- tuned copper F-human, four feeds 0 / 90 / 180 / 270 deg at the "
              "two-port root: 10 g C4 orbit (WF-10 step 5)",
        window_size=(1400, 630),
        comm=comm,
    )
    t_fig = comm.allreduce(time.perf_counter() - t_fig, op=MPI.MAX)

    # ---- the 10 g lattice (gate :401-438) ---------------------------------------------------
    fields = {"4p": ef4, "ccw": e_ccw}
    containment, walls, masses = [], [], {}

    def ball(key, c, label):
        c = tuple(float(v) for v in c)
        containment.append((label, c, _containment(c, centre, radius, half_height, a10)))
        t = time.perf_counter()
        r = mass_averaged_sar(fields[key], sigma=sigma, rho=rho_field, center=c, radius=a10,
                              comm=comm, quadrature_degree=QUADRATURE_DEGREE)
        r["wall_s"] = comm.allreduce(time.perf_counter() - t, op=MPI.MAX)
        walls.append(r["wall_s"])
        if key == "4p":
            masses[label] = r["mass_kg"] * 1e3
        return r

    comm.Barrier()
    t_sar = time.perf_counter()
    lattice = _lattice(az17, z_c)
    sar4 = {}
    reduced = False
    for idx, (label, _r, dz, _d, c) in enumerate(lattice):
        if reduced and dz != 0.0:
            continue
        sar4[label] = ball("4p", c, label)["averaged_sar_w_per_kg"]
        if idx == 0 and walls[-1] > FIRST_BALL_STOP_S:
            reduced = True
            P(f"REDUCED LATTICE: first ball {walls[-1]:.2f} s > {FIRST_BALL_STOP_S:.0f} s -- "
              f"z = z_c rows only (pre-registered, not a red)")
    used4 = [label for label, *_x in lattice if label in sar4]
    key4 = max(used4, key=lambda k: sar4[k])
    coords = {label: np.array(c) for label, *_x, c in lattice}
    sar2 = {}
    for label, _r, dz, _d, c in lattice:
        if dz == 0.0:
            sar2[label] = ball("ccw", c, label)["averaged_sar_w_per_kg"]

    # ---- coverage ------------------------------------------------------------------------------
    coverage_radius = float(np.hypot(radius, half_height) + COVERAGE_BALL_PAD_M)
    conductor_limit = F_HUMAN_RING_RADIUS - CONDUCTOR_CLEARANCE_M
    t = time.perf_counter()
    whole4 = mass_averaged_sar(ef4, sigma=sigma, rho=rho_field,
                               center=tuple(float(v) for v in centre), radius=coverage_radius,
                               comm=comm, quadrature_degree=QUADRATURE_DEGREE)
    t_cov = comm.allreduce(time.perf_counter() - t, op=MPI.MAX)
    t_sar = comm.allreduce(time.perf_counter() - t_sar, op=MPI.MAX)
    power_miss = abs(whole4["dissipated_power_w"] / tagged4["dissipated_power_w"] - 1.0)
    mass_miss = abs(whole4["mass_kg"] / (PHANTOM_RHO_KG_PER_M3 * tagged4["volume_m3"]) - 1.0)

    # ---- readings --------------------------------------------------------------------------------
    rows4 = _zc_rows(sar4, lattice)
    rows2 = _zc_rows(sar2, lattice)
    c4_4 = _c4_misses(rows4)
    c4_2 = _c4_misses(rows2)
    c4_max4 = max(c4_4.values())
    c4_max2 = max(c4_2.values())
    sp4 = _spreads(rows4)
    sp2 = _spreads(rows2)
    sp359 = np.array(TWO_PORT_SPREADS)
    repro2 = float(np.max(np.abs(rows2 / TWO_PORT_ZC_RECORD - 1.0)))
    below = [bool(sp4[i] <= VERDICT_FRACTION * _spreads(TWO_PORT_ZC_RECORD)[i]) for i in range(3)]
    verdict = ("drive anisotropy" if all(below)
               else "lattice resolution next suspect" if not any(below) else "mixed")
    worst_cont = max(x[2] for x in containment)
    mean4 = tagged4["mean_sar_w_per_kg"]
    sar_star4 = sar4[key4]
    c_star4 = coords[key4]
    ratio_sp = [float(sp4[i] / sp2[i]) for i in range(3)]

    # ---- (O): the 12 z_c balls are one C4 orbit per radius ---------------------------------------
    zc_pts = np.array([c for _l, _r, dz, _d, c in lattice if dz == 0.0], dtype=float).reshape(4, 3, 3)
    r_orbit = max(float(np.max(np.abs(_rotate_z(zc_pts[0], m * np.pi / 2.0) - zc_pts[m])))
                  for m in (1, 2, 3)) / float(max(LATTICE_R_M))

    # ---- (X1)_4 / (X1)_2: DG0 loss maps sum to the phantom power ---------------------------------
    dg0 = fem.functionspace(msh, ("DG", 0))
    v = ufl.TestFunction(dg0)
    dx = ufl.Measure("dx", domain=msh, subdomain_data=cell_tags)
    n_loc = dg0.dofmap.index_map.size_local
    t_dg = time.perf_counter()
    vol = _dg0_vec(ufl.conj(v) * ufl.dx(domain=msh))
    loss4 = _dg0_vec(0.5 * sigma * ufl.inner(ef4, ef4) * ufl.conj(v) * dx(PHANTOM_CELL_TAG))
    loss2 = _dg0_vec(0.5 * sigma * ufl.inner(e_ccw, e_ccw) * ufl.conj(v) * dx(PHANTOM_CELL_TAG))
    t_dg = comm.allreduce(time.perf_counter() - t_dg, op=MPI.MAX)
    x1_4 = float(comm.allreduce(np.sum(loss4[:n_loc]), op=MPI.SUM))
    x1_2 = float(comm.allreduce(np.sum(loss2[:n_loc]), op=MPI.SUM))
    r_x1_4 = abs(x1_4 / p4_ph - 1.0)
    r_x1_2 = abs(x1_2 / p2_ph - 1.0)

    _m, _rss_per, rss_sum = _rss_gib(comm)

    # ---- record comparison (relative, 10**-(d-2)); a miss stops the run -----------------------
    recs = []  # (label, got, want_str, rel, tol)

    def rec_check(label, got, want, d):
        w = float(want)
        tol = 10.0 ** (-(d - 2))
        rel = abs(float(got) / w - 1.0)
        recs.append((label, float(got), want, rel, tol))

    rec_exact_bad = []
    if out["cells"] != REC["cells"]:
        rec_exact_bad.append(("cells", out["cells"], REC["cells"]))
    for lab in used4:
        rec_check(f"sar4 {lab}", sar4[lab], REC["sar4"][lab], 10)
        rec_check(f"mass {lab}", masses[lab], REC["mass_g"][lab], 7)
    for i, lab in enumerate([x[0] for x in lattice if x[2] == 0.0]):
        rec_check(f"sar2 {lab}", rows2.reshape(-1)[i], f"{TWO_PORT_ZC_RECORD.reshape(-1)[i]:.9e}", 10)
    for m in (1, 2, 3):
        rec_check(f"c4_4 m={m} %", c4_4[m] * 100, *REC["c4_4"][m])
        rec_check(f"c4_2 m={m} %", c4_2[m] * 100, *REC["c4_2"][m])
    for i in range(3):
        rec_check(f"spread4 r={LATTICE_R_M[i] * 1e3:.0f} %", sp4[i] * 100, *REC["sp4"][i])
        rec_check(f"spread2 r={LATTICE_R_M[i] * 1e3:.0f} %", sp2[i] * 100, *REC["sp2"][i])
    rec_check("sar_star4", sar_star4, REC["sar_star4"], 10)
    rec_check("p_src4", p4_src, REC["p_src4"], 10)
    rec_check("p_phantom4", p4_ph, REC["p_phantom4"], 10)
    rec_check("k=0 share", spec4[ref_sign][0], *REC["k0"])
    rec_check("k=-1 (bin 15) share", spec4[ref_sign][15], *REC["k_bin15"])
    rec_check("SAR10g(c*4)/P_src,4", sar_star4 / p4_src, *REC["per_watt"])
    rec_check("peak-to-mean", sar_star4 / mean4, *REC["peak_to_mean"])
    rec_check("P_phantom/P_src,4", p4_ph / p4_src, *REC["pph_over_psrc"])
    rec_bad = [r for r in recs if r[3] > r[4]]
    if key4 != REC["c_star4"]:
        rec_exact_bad.append(("c_star4", key4, REC["c_star4"]))
    if verdict != REC["verdict"] or below != REC["below"]:
        rec_exact_bad.append(("verdict", (verdict, below), (REC["verdict"], REC["below"])))

    # ---- printed output -----------------------------------------------------------------------------
    P(f"\n[anchor] (S_rep) s_node vs NODE_JSON at {C_IN_F_HUMAN_PF} pF: {r_srep:.4e} <= "
      f"{REUSE_REPRODUCTION_RTOL:g} (round-off reading, W1 :6337 2.578e-13)")
    P(f"[four-port] feeds {d4['feeds4']} (two-port {d4['feeds2']}), azimuth steps "
      + " / ".join(f"{d:.6f}" for d in d_az4) + f" deg (90 within {AZIMUTH_MATCH_DEG:g}); "
      f"kappa off {d4['kappa_off']}; kinds {d4['kinds4']}")
    P(f"[four-port] (X4-P) rel {d4['r_x4p']:.4e} (band {DISCRETE_IDENTITY_RTOL:g}; W1 :6371 3.286e-13); "
      f"(X4-L) {d4['r_x4l']:.4e} (band {UNIT_WEIGHT_CONTROL_RTOL:g}; W1 :6372 6.920e-16)")
    P(f"[four-port] P_src,4 {p4_src:.9e} W, P_phantom,4 {p4_ph:.9e} W (W1 :6371 {REC['p_src4']}, "
      f"{REC['p_phantom4']})")
    P("[four-port] Y_in four feeds (printed): "
      + "; ".join(f"{pid} {y.real:+.5e} {y.imag:+.5e}j S" for pid, y in zip(d4["feeds4"], d4["y4"])))
    P(f"[four-port] feed-ring k shares: k = 0 {spec4[ref_sign][0]:.4e}, k = -1 (bin 15) "
      f"{spec4[ref_sign][15]:.6f} (W1 :6374 {REC['k0'][0]}, {REC['k_bin15'][0]}); signed lag "
      f"ring+ {k4[1.0]['lag']:.6f}, ring- {k4[-1.0]['lag']:.6f}")
    P(f"\n[10 g] four-port table ({len(used4)} points{', REDUCED' if reduced else ''}); "
      f"two-port ccw beside it on the z_c rows (record reprints, W/kg):")
    for k in used4:
        extra = f"   two-port ccw {sar2[k]:.9e}" if k in sar2 else ""
        P(f"   {k:<24s} 4p {sar4[k]:.9e}{extra}{'   <- c*4' if k == key4 else ''}")
    if reduced:
        P("[10 g] the six off-plane rows were NOT read (reduced lattice)")
    P(f"[10 g] c*4 = {key4} ({c_star4[0] * 1e3:.3f}, {c_star4[1] * 1e3:.3f}, {c_star4[2] * 1e3:.3f}) mm, "
      f"SAR10g(c*4) {sar_star4:.9e} W/kg (PREDICTED r = 45 mm, z_c; W1 :6394 {REC['c_star4']}); "
      f"a lattice argmax, not a global maximum")
    P("[10 g] (C4_4) four-port, row m vs row 0, rel L2 over 3 radii: "
      + " / ".join(f"m={m} {v * 100:.4f}%" for m, v in c4_4.items())
      + f"; max {c4_max4 * 100:.4f}% (band {C4_COVARIANCE_BAND * 100:.1f}%, ASSERTED; "
      f"PREDICTED {PREDICTED_C4[0] * 100:g}-{PREDICTED_C4[1] * 100:g}% per m; "
      f"W1 :6395 1.8127 / 0.7215 / 0.2645)")
    P("[10 g] control (ASSERTED > band): two-port ccw, same 12 balls: "
      + " / ".join(f"m={m} {v * 100:.4f}%" for m, v in c4_2.items())
      + f"; max {c4_max2 * 100:.4f}% = {c4_max2 / C4_COVARIANCE_BAND:.1f}x the band "
      f"(W1 :6396 17.3831 / 0.8080 / 16.1396; 359 {TWO_PORT_C4_MISSES})")
    blind4 = c4_4[2] <= C4_COVARIANCE_BAND
    blind2 = c4_2[2] <= C4_COVARIANCE_BAND
    P(f"[10 g] 180 degree blindness, m = 2 alone: four-port {c4_4[2] * 100:.4f}%, two-port "
      f"{c4_2[2] * 100:.4f}% -- both inside the band: {blind4 and blind2} "
      f"(PREDICTED inside band (printed only); PREDICTED 0.72 / 0.81%; the gate takes the max over m = 1, 2, 3)")
    P(f"[10 g] two-port in-window vs 359's logged table: max rel {repro2:.3e} "
      f"(printed; PREDICTED <= {PREDICTED_TWO_PORT_REPRO:g}; W1 :6397 1.912e-10)")
    P("[10 g] spread per radius (the record): "
      + "; ".join(f"r={r * 1e3:.0f} mm four-port {sp4[i] * 100:.4f}%, two-port in-window "
                  f"{sp2[i] * 100:.4f}%, 359 {sp359[i] * 100:.1f}% (1/3: {sp359[i] * 100 / 3:.2f}%)"
                  for i, r in enumerate(LATTICE_R_M))
      + f" (PREDICTED four-port {PREDICTED_SPREAD4[0] * 100:g}-{PREDICTED_SPREAD4[1] * 100:g}%)")
    P("[10 g] spread four-port / two-port per radius: " + " / ".join(f"{x:.3f}" for x in ratio_sp)
      + " (PREDICTED 0.19 / 0.20 / 0.04, printed)")
    P(f"[10 g] VERDICT (printed; PREDICTED drive anisotropy): {verdict} "
      f"(per radius <= 1/3 of 359: {below})")
    P(f"[10 g] (i)4 coverage ball {coverage_radius * 1e3:.4f} mm (< conductor limit "
      f"{conductor_limit * 1e3:.1f} mm); power rel {power_miss:.3e}, mass rel {mass_miss:.3e} "
      f"(band {EXACT_IDENTITY_RTOL:g}); (Sigma_P)4 rel {r_sigma4:.3e}; (iv) worst containment "
      f"{worst_cont * 1e3:.4f} mm over {len(containment)} balls (<= 0)")
    P(f"[10 g] per watt (untuned four-port at the two-port root, degree 1): SAR10g(c*4)/P_src,4 "
      f"{sar_star4 / p4_src:.6e} W/kg/W (two-port 5.607e-02); peak-to-mean {sar_star4 / mean4:.4f} "
      f"(two-port 1.309); P_phantom/P_src,4 {p4_ph / p4_src:.6e} (two-port 0.0947)")
    two = out["anchors"]
    for sense in ("ccw", "cw"):
        rel_w, _b = two["p_w"][sense]
        r_lin, r_y = two["lin"][sense]
        P(f"[anchor] two-port [{sense}] (P_w) {rel_w:.3e} <= {DISCRETE_IDENTITY_RTOL:g}; (Lin) "
          f"{r_lin:.3e} <= {UNIT_WEIGHT_CONTROL_RTOL:g}, Y_in {r_y:.3e} <= {Y0_CODING_BAND:g}")
    P(f"[anchor] (P_RC) {out['p_rc_node']:.3e} (asserted inside the builder)")
    P(f"[example] (O) z_c orbit: max |R(m*90) row0 - row m| / max(LATTICE_R_M) {r_orbit:.3e} "
      f"<= {EXACT_IDENTITY_RTOL:g} (PREDICTED <= 1e-15)")
    P(f"[example] (X1)_4 DG0 loss map sum rel {r_x1_4:.3e}; (X1)_2 {r_x1_2:.3e} "
      f"(<= {EXACT_IDENTITY_RTOL:g} ASSERTED)")
    n_ok = len(recs) - len(rec_bad)
    P(f"\n[records] {n_ok} / {len(recs)} relative records inside 10**-(d-2) against {W1}; "
      f"max rel {max(r[3] for r in recs):.3e}; exact: cells {out['cells']}, c*4 {key4!r}, verdict {verdict!r}")
    for lab, got, want, rel, tol in rec_bad:
        P(f"   OFF TOLERANCE {lab}: got {got:.9e}, W1 {want}, rel {rel:.3e} > {tol:.0e}")
    for lab, got, want in rec_exact_bad:
        P(f"   EXACT MISMATCH {lab}: got {got!r}, W1 {want!r}")
    P(f"[price] build {times['build']:.1f} s, node pair {times['node_pair']:.1f} s, superposition "
      f"{times['superposition']:.1f} s; four-port solves {sum(d4['t4']):.1f} s, partitions "
      f"{d4['t_sup']:.1f} s; {len(walls)} 10 g balls (per ball mean {np.mean(walls):.2f} s max "
      f"{np.max(walls):.2f} s) + coverage {t_cov:.2f} s, SAR phase {t_sar:.1f} s; DG0 maps "
      f"{t_dg:.1f} s; setup figure {t_fig:.1f} s; RSS summed {rss_sum:.2f} GiB; -n {comm.size}")

    # ---- asserts (imported bands only), before the artifacts -----------------------------------
    assert not rec_exact_bad, rec_exact_bad
    assert not rec_bad, [(r[0], r[3], r[4]) for r in rec_bad]
    assert r_srep <= REUSE_REPRODUCTION_RTOL, f"(S_rep) {r_srep:.6e}"
    for sense in ("ccw", "cw"):
        rel_w, _b = out["anchors"]["p_w"][sense]
        r_lin, r_y = out["anchors"]["lin"][sense]
        assert rel_w <= DISCRETE_IDENTITY_RTOL, f"(P_w) {sense} {rel_w:.3e}"
        assert r_lin <= UNIT_WEIGHT_CONTROL_RTOL, f"(Lin) {sense} {r_lin:.3e}"
        assert r_y <= Y0_CODING_BAND, f"(Lin) Y_in {sense} {r_y:.3e}"
    assert d4["kappa_off"], "kappa must be off on every sheet"
    assert az_ok, d_az4
    assert feeds_ok, (d4["feeds4"], d4["feeds2"])
    assert d4["kinds4"] == ["solve", "held", "held", "held"], d4["kinds4"]
    assert d4["r_x4p"] <= DISCRETE_IDENTITY_RTOL, f"(X4-P) {d4['r_x4p']:.3e}"
    assert d4["r_x4l"] <= UNIT_WEIGHT_CONTROL_RTOL, f"(X4-L) {d4['r_x4l']:.3e}"
    assert r_sigma4 <= EXACT_IDENTITY_RTOL, f"(Sigma_P)4 {r_sigma4:.3e}"
    assert coverage_radius < conductor_limit, (coverage_radius, conductor_limit)
    assert power_miss <= EXACT_IDENTITY_RTOL, f"(i)4 power {power_miss:.3e}"
    assert mass_miss <= EXACT_IDENTITY_RTOL, f"(i)4 mass {mass_miss:.3e}"
    assert worst_cont <= 0.0, f"(iv) containment {worst_cont:.4e} m"
    assert r_orbit <= EXACT_IDENTITY_RTOL, f"(O) {r_orbit:.3e}"
    assert r_x1_4 <= EXACT_IDENTITY_RTOL, f"(X1)_4 {r_x1_4:.3e}"
    assert r_x1_2 <= EXACT_IDENTITY_RTOL, f"(X1)_2 {r_x1_2:.3e}"
    assert c4_max2 > C4_COVARIANCE_BAND, f"control: two-port C4 max {c4_max2:.4e} inside band"
    assert c4_max4 <= C4_COVARIANCE_BAND, f"(C4_4) {c4_4} -- never re-band"

    # ---- ParaView maps: DG0, real, named once -----------------------------------------------------
    def named(name, values):
        f = fem.Function(dg0, name=name)
        f.x.array[:] = values
        f.x.scatter_forward()
        return f

    sar_loc4 = np.zeros_like(loss4)
    sar_loc2 = np.zeros_like(loss2)
    np.divide(loss4, PHANTOM_RHO_KG_PER_M3 * vol, out=sar_loc4, where=(vol > 0.0) & (loss4 > 0.0))
    np.divide(loss2, PHANTOM_RHO_KG_PER_M3 * vol, out=sar_loc2, where=(vol > 0.0) & (loss2 > 0.0))
    sar_loc4 = sar_loc4 / p4_src
    sar_loc2 = sar_loc2 / p2_src
    ratio = np.zeros_like(sar_loc4)
    np.divide(sar_loc4, sar_loc2, out=ratio, where=(sar_loc4 > 0.0) & (sar_loc2 > 0.0))

    tdim = msh.topology.dim
    imap = msh.topology.index_map(tdim)
    n_cells = imap.size_local + imap.num_ghosts
    mids = compute_midpoints(msh, tdim, np.arange(n_cells, dtype=np.int32))
    cell_dof = np.asarray(dg0.dofmap.list).reshape(-1)[:n_cells]
    zc_labels = [x[0] for x in lattice if x[2] == 0.0]
    zc_ball = np.zeros_like(loss4)
    c4_orbit = np.zeros_like(loss4)
    for k, label in enumerate(zc_labels):
        m_idx, r_idx = divmod(k, 3)
        inside = np.linalg.norm(mids - coords[label][None, :], axis=1) < a10
        zc_ball[cell_dof[inside]] = float(1 + 3 * m_idx + r_idx)
        c4_orbit[cell_dof[inside]] = float(r_idx + 1)

    export = {
        "SAR_local_4port_per_W": named("SAR_local_4port_per_W", sar_loc4),
        "SAR_local_2port_ccw_per_W": named("SAR_local_2port_ccw_per_W", sar_loc2),
        "SAR_local_4p_over_2p_per_W": named("SAR_local_4p_over_2p_per_W", ratio),
        "ZcBall": named("ZcBall", zc_ball),
        "C4Orbit": named("C4Orbit", c4_orbit),
    }
    OUTPUT_DIR.mkdir(exist_ok=True)
    xdmf_path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_combined", msh, cell_tags, export, comm=comm
    )
    adopt_host_ownership(OUTPUT_DIR, comm=comm)

    csv_path = OUTPUT_DIR / f"{BASENAME}_zc_rows.csv"
    if rank0:
        with open(csv_path, "w") as fh:
            fh.write("label,r_m,az_offset_deg,x_m,y_m,z_m,SAR10g_4port_per_W,SAR10g_2port_ccw_per_W,"
                     "mass_kg\n")
            for label, r_, dz_, daz_, c_ in lattice:
                if dz_ != 0.0:
                    continue
                fh.write(f"\"{label}\",{r_:.6e},{daz_:g},{c_[0]:.6e},{c_[1]:.6e},{c_[2]:.6e},"
                         f"{sar4[label] / p4_src:.6e},{sar2[label] / p2_src:.6e},"
                         f"{masses[label] * 1e-3:.6e}\n")

    if rank0:  # xdmf_path is None on ranks != 0 (item 506)
        P(f"\n[output] {xdmf_path.name}: DG0 " + ", ".join(export) + ", CellTags")
    P(f"[output] {csv_path.name}")
    P("[paraview] threshold `ZcBall` >= 1 coloured by `SAR_local_4port_per_W`, then by "
      "`SAR_local_2port_ccw_per_W` on one colour scale; slice z = z_c on "
      "`SAR_local_4p_over_2p_per_W`; group by `C4Orbit`.")
    P(f"\nEX-127 all assertions hold: gate anchors, (O), (X1)_4 / (X1)_2, the asserted control, "
      f"{len(recs)} records. Total elapsed {time.perf_counter() - started:.1f} s.")


if __name__ == "__main__":
    main()
