"""`WF-16` step 9 part B (§9 item 382b, chain U36): ``Study.b1_plus_at`` and
``Study.sar(mass_g=10)`` on the tuned copper F-human Study.

The Study of part A (``FHumanBirdcageGeometry(as_hole=True)``,
``TwoPortQuadratureDrive(C_in, esr_ohm=0.2)``, ``LeontovichWall()``) is read
through the public B / SAR API and reproduces H9′'s records — item 358's mean
``|B₁⁺|`` (``20261010T170356Z_WF-10-step2-W1.log``) and item 359's 18-point
10 g lattice (``20261010T171705Z_WF-10-step2s-W1.log``) — at
``REUSE_REPRODUCTION_RTOL``; the feed-bisector mirror gates at
``C4_COVARIANCE_BAND``.  Every band is imported; none is set here.  The Study
never calls ``build_tuned_quadrature``.

Anchors (asserted):
  (B)     mean over the 41 sample points of ``b1_plus_at(·, drives=("ccw",))``
          == ``H9P_MEAN_B1_PLUS_CCW_RECORD_T`` within ``REUSE_REPRODUCTION_RTOL``.
  (Sense) mean|B₁⁺|_ccw / mean|B₁⁻|_ccw > 1 (record 21.5155).
  (M_B)   |B₁⁻|_cw(Mx) vs |B₁⁺|_ccw(x) rel L2 <= ``C4_COVARIANCE_BAND``.
  (S10)   the 18 ccw and 18 cw(M) 10 g values == the records within
          ``REUSE_REPRODUCTION_RTOL``; c* == 359's, SAR10g(c*) == the record.
  (iii-a) mirror at c*; (iii-b) lattice rel L2 — both <= ``C4_COVARIANCE_BAND``.

Negative control (rule (e), asserted, backed): (N1) the mis-paired sense
|B₁⁺|_cw(Mx) vs |B₁⁺|_ccw(x) > ``C4_COVARIANCE_BAND`` (record 9.5476e-01,
19.1x); the factor is printed, not asserted.

Opt-in: ``FEM_EM_WF16_S9B=1`` (unset ⇒ the solve node skips; any other value
raises); width 8 asserted when set.  Scope: degree 1, 64 MHz, κ off, one copper
hole fixture, one phantom.  Not claimed: homogeneity, an absolute |B₁⁺| or SAR,
a global SAR maximum (c* is a lattice argmax), compliance, degree 2, C16.
"""

from __future__ import annotations

import os
import time

import numpy as np
import pytest
from mpi4py import MPI

from fem_em_solver import FHumanBirdcageGeometry, TwoPortQuadratureDrive
from fem_em_solver.study import _TUNED_SENSES, SAR_QUADRATURE_DEGREE
from tests.complex_mode import complex_only
from tests.package.test_wf16_study_tuned_f_human import _c_f, _study
from tests.validation.test_birdcage_b1_plus_map import (
    C4_COVARIANCE_BAND,
    PHANTOM_RHO_KG_PER_M3,
    _relative_l2,
    _sample_points,
)
from tests.validation.test_birdcage_b1_quadrature import _mirror_xy, _read_senses
from tests.validation.test_mass_averaged_sar_standard_masses import QUADRATURE_DEGREE
from tests.validation.test_port_birdcage_lumped_column import PHANTOM_CELL_TAG
from tests.validation.test_port_birdcage_ring_matrix import REUSE_REPRODUCTION_RTOL
from tests.validation.test_wf7_f_human_sar10g_hotspot import _phantom_geometry
from tests.validation.test_wf10_f_human_copper_esr_root import ESR_OHM
from tests.validation.test_wf10_f_human_tuned_b1_sar import NODE_JSON
from tests.validation.test_wf10_f_human_tuned_sar10g import _lattice

S9B_ENV = "FEM_EM_WF16_S9B"
RECORD_RANK_WIDTH = 8  # the width H9′'s records were taken at

# ---- record constants (printed comparands, transcribed; never a band) --------
#: Item 358's mean |B₁⁺|_ccw over the 41 sample points
#: (`20261010T170356Z_WF-10-step2-W1.log:6339`, printed at .6e).
H9P_MEAN_B1_PLUS_CCW_RECORD_T = 2.884091e-07
#: Item 358's (Sense) ratio (`…step2-W1.log:6333`), printed beside.
H9P_SENSE_RATIO_RECORD = 21.5155
#: Item 358's (M_B) mirror and mis-paired control (`…step2-W1.log:6335–6336`).
H9P_M_B_RECORD = 1.9647e-03
H9P_N1_RECORD = 9.5476e-01
#: Item 359's SAR10g(c*) (`20261010T171705Z_WF-10-step2s-W1.log:6372`).
H9P_SAR10G_CSTAR_RECORD_W_PER_KG = 4.852895959e-03
#: Item 359's c* label (`…step2s-W1.log:6372`).
H9P_CSTAR_LABEL = "r=45 dz=+0 az17+0"
#: Item 359's (iii-a) / (iii-b) readings (`…step2s-W1.log:6373–6374`), printed beside.
H9P_MIRROR_CSTAR_RECORD = 0.6413e-2
H9P_MIRROR_L2_RECORD = 1.2789e-2
#: Item 359's 18-row table (`…step2s-W1.log:6354–6371`): (label, ccw, cw(M)) W/kg.
H9P_SAR10G_TABLE_RECORD = (
    ("r=15 dz=+0 az17+0", 6.830007305e-04, 6.606850198e-04),
    ("r=30 dz=+0 az17+0", 2.230610345e-03, 2.138875190e-03),
    ("r=45 dz=+0 az17+0", 4.852895959e-03, 4.821772251e-03),
    ("r=15 dz=+0 az17+90", 5.781188640e-04, 5.968531327e-04),
    ("r=30 dz=+0 az17+90", 1.803115212e-03, 1.887962586e-03),
    ("r=45 dz=+0 az17+90", 4.026865820e-03, 4.073878924e-03),
    ("r=15 dz=+0 az17+180", 6.587520880e-04, 6.671098227e-04),
    ("r=30 dz=+0 az17+180", 2.230637254e-03, 2.224615414e-03),
    ("r=45 dz=+0 az17+180", 4.816772104e-03, 4.848794891e-03),
    ("r=15 dz=+0 az17+270", 5.882637672e-04, 5.803609204e-04),
    ("r=30 dz=+0 az17+270", 1.879465733e-03, 1.896434461e-03),
    ("r=45 dz=+0 az17+270", 4.063626880e-03, 4.053371643e-03),
    ("r=15 dz=-50 az17+0", 1.635320677e-03, 1.636807026e-03),
    ("r=30 dz=-50 az17+0", 2.243679310e-03, 2.244344497e-03),
    ("r=45 dz=-50 az17+0", 3.414395055e-03, 3.391893557e-03),
    ("r=15 dz=+50 az17+0", 1.578408674e-03, 1.587542278e-03),
    ("r=30 dz=+50 az17+0", 2.091815655e-03, 2.096323755e-03),
    ("r=45 dz=+50 az17+0", 3.209282450e-03, 3.198000962e-03),
)
#: 358's sample-point count (`…step2-W1.log:6332`: "41 sample points, 41 valid").
H9P_N_SAMPLE_POINTS = 41


def _knob() -> bool:
    raw = os.environ.get(S9B_ENV, "").strip()
    if raw == "":
        return False
    if raw == "1":
        return True
    raise ValueError(f"{S9B_ENV}={raw!r}: only unset/empty or '1' are accepted")


def _p0(msg):
    if MPI.COMM_WORLD.rank == 0:
        print(f"[WF-16 s9b] {msg}", flush=True)


def test_wf16_s9b_hooks_and_records_no_mesh():
    """W0: the hooks' refusals before any mesh; the record constants' coding."""
    assert SAR_QUADRATURE_DEGREE == QUADRATURE_DEGREE
    assert _TUNED_SENSES == ("ccw", "cw")
    study = _study()
    assert study.density_kg_m3 == PHANTOM_RHO_KG_PER_M3
    assert int(study.geometry.phantom_tag) == PHANTOM_CELL_TAG
    for bad in (
        {"mass_g": 0},
        {"mass_g": 10.0, "radius_m": 0.01},
        {"mass_g": None},
        {"mass_g": 10.0, "drives": ("P17",)},
    ):
        with pytest.raises(ValueError):
            study.sar(centres=[(0.0, 0.0, 0.0)], **bad)
    with pytest.raises(ValueError):
        study.b1_plus(drives=("P1",))
    with pytest.raises(ValueError):
        study.b1_plus_at(np.zeros((1, 3)), drives=("x",))
    # The solve's own refusals, reached through the hooks, before any mesh.
    solid = _study(geom=FHumanBirdcageGeometry())
    with pytest.raises((TypeError, ValueError)):
        solid.b1_plus()
    with pytest.raises((TypeError, ValueError)):
        solid.sar(centres=[(0.0, 0.0, 0.0)])
    bare = _study(conductor=None)
    with pytest.raises((TypeError, ValueError)):
        bare.b1_plus()
    for s in (study, solid, bare):
        assert getattr(s, "result", None) is None, "a refusal built a mesh"
        assert getattr(s, "_b1_cg1", None) is None
    # Record coding.
    labels = [row[0] for row in H9P_SAR10G_TABLE_RECORD]
    assert labels == [lab for lab, *_x in _lattice(11.25, 0.0)]
    best = max(H9P_SAR10G_TABLE_RECORD, key=lambda row: row[1])
    assert best[0] == H9P_CSTAR_LABEL and best[1] == H9P_SAR10G_CSTAR_RECORD_W_PER_KG
    sar_c = np.array([row[1] for row in H9P_SAR10G_TABLE_RECORD])
    sar_m = np.array([row[2] for row in H9P_SAR10G_TABLE_RECORD])
    l2 = _relative_l2(sar_m, sar_c, np.ones(sar_c.size, dtype=bool))
    star = abs(best[2] - best[1]) / best[1]
    _p0(
        f"(R0) table re-reduced: (iii-a) {star:.4e} vs record {H9P_MIRROR_CSTAR_RECORD:.4e}; "
        f"(iii-b) {l2:.4e} vs record {H9P_MIRROR_L2_RECORD:.4e}"
    )
    assert abs(star - H9P_MIRROR_CSTAR_RECORD) <= 5e-7
    assert abs(l2 - H9P_MIRROR_L2_RECORD) <= 5e-7
    assert isinstance(study.excitation, TwoPortQuadratureDrive)
    assert study.excitation.c_in_f == _c_f() and study.excitation.esr_ohm == ESR_OHM


@complex_only
def test_wf16_s9b_tuned_b1_sar():
    if not _knob():
        pytest.skip(f"{S9B_ENV} unset: heavy -n 8 F-human window, environment-selected")
    comm = MPI.COMM_WORLD
    assert comm.size == RECORD_RANK_WIDTH, f"records are -n {RECORD_RANK_WIDTH}; got {comm.size}"
    assert os.path.exists(NODE_JSON), f"{NODE_JSON} absent: item 382b stops (BLOCKED)"
    t_window = time.perf_counter()
    study = _study()

    # ---- B -------------------------------------------------------------------------
    comm.Barrier()
    t0 = time.perf_counter()
    b1 = study.b1_plus(comm=comm)  # solves first (study.result is None)
    comm.Barrier()
    t_b = time.perf_counter() - t0
    res = study.result
    assert set(b1) == {"ccw", "cw"}
    points = _sample_points({"mesh": res.mesh, "cell_tags": res.cell_tags})
    az_m = float(res.az_mirror)
    vals, valid = study.b1_plus_at(points, drives=("ccw",), comm=comm)["ccw"]
    ccw_x = _read_senses(b1["ccw"], points)
    cw_mx = _read_senses(b1["cw"], _mirror_xy(points, az_m))
    mask = np.logical_and.reduce([valid, ccw_x[2], cw_mx[2]])
    n_pts, n_mask = int(points.shape[0]), int(np.count_nonzero(mask))
    plus_x, minus_x = ccw_x[0], ccw_x[1]
    api_vs_read = float(np.max(np.abs(vals[mask] - plus_x[mask])))
    mean_b = float(np.mean(vals[mask]))
    r_b = abs(mean_b - H9P_MEAN_B1_PLUS_CCW_RECORD_T) / H9P_MEAN_B1_PLUS_CCW_RECORD_T
    sense_ratio = float(np.mean(plus_x[mask]) / np.mean(minus_x[mask]))
    mirror = _relative_l2(cw_mx[1], plus_x, mask)
    mispaired = _relative_l2(cw_mx[0], plus_x, mask)
    _p0(
        f"B: solve + two CG1 projections {t_b:.1f} s (build {res.times.get('build', 0.0):.1f} s, "
        f"node pair {res.times.get('node_pair', 0.0):.1f} s); {res.cells} cells; {n_pts} sample "
        f"points, {n_mask} valid on every image; b1_plus_at vs _read_senses max |d| "
        f"{api_vs_read:.3e} T"
    )
    _p0(
        f"(B) mean|B1+|_ccw {mean_b:.9e} T vs record {H9P_MEAN_B1_PLUS_CCW_RECORD_T:.6e} T: rel "
        f"{r_b:.3e} (band {REUSE_REPRODUCTION_RTOL:g}, ASSERTED)"
    )
    _p0(
        f"(Sense) mean|B1+|_ccw / mean|B1-|_ccw {sense_ratio:.4f} (> 1 ASSERTED; record "
        f"{H9P_SENSE_RATIO_RECORD})"
    )
    _p0(
        f"(M_B) feed-bisector mirror (az {az_m:.6f} deg) rel L2 {mirror:.4e} (band "
        f"{C4_COVARIANCE_BAND:g}, ASSERTED; record {H9P_M_B_RECORD:.4e})"
    )
    _p0(
        f"(N1) mis-paired |B1+|_cw(Mx) vs |B1+|_ccw(x) {mispaired:.4e} = "
        f"{mispaired / C4_COVARIANCE_BAND:.1f}x band (> band ASSERTED; factor printed; record "
        f"{H9P_N1_RECORD:.4e})"
    )

    # ---- SAR, 10 g lattice ------------------------------------------------------------
    centre, _radius, _hh, _lo, _hi = _phantom_geometry(
        res.mesh, res.cell_tags, PHANTOM_CELL_TAG, comm
    )
    z_c = float(centre[2])
    az17 = float(res.az_feeds[0])
    lattice = _lattice(az17, z_c)
    labels = [lab for lab, *_x in lattice]
    centres = [tuple(float(v) for v in c) for *_x, c in lattice]
    mirrored = [tuple(float(v) for v in _mirror_xy(np.array([c]), az_m)[0]) for c in centres]
    comm.Barrier()
    t0 = time.perf_counter()
    ccw = study.sar(mass_g=10.0, centres=centres, drives=("ccw",), comm=comm)["ccw"]
    cw = study.sar(mass_g=10.0, centres=mirrored, drives=("cw",), comm=comm)["cw"]
    comm.Barrier()
    t_sar = time.perf_counter() - t0
    sar_c = np.array([d["averaged_sar_w_per_kg"] for d in ccw])
    sar_m = np.array([d["averaged_sar_w_per_kg"] for d in cw])
    rec = {row[0]: row[1:] for row in H9P_SAR10G_TABLE_RECORD}
    r_ccw = np.array([abs(sar_c[i] / rec[lab][0] - 1.0) for i, lab in enumerate(labels)])
    r_cw = np.array([abs(sar_m[i] / rec[lab][1] - 1.0) for i, lab in enumerate(labels)])
    i_star = int(np.argmax(sar_c))
    key_star = labels[i_star]
    sar_star = float(sar_c[i_star])
    r_star = abs(sar_star - H9P_SAR10G_CSTAR_RECORD_W_PER_KG) / H9P_SAR10G_CSTAR_RECORD_W_PER_KG
    mirror_star = abs(sar_m[i_star] - sar_star) / sar_star
    mirror_l2 = _relative_l2(sar_m, sar_c, np.ones(sar_c.size, dtype=bool))
    for i, lab in enumerate(labels):
        _p0(
            f"  {lab:<22s} ccw {sar_c[i]:.9e} (rec {rec[lab][0]:.9e}, rel {r_ccw[i]:.2e}) cw(M) "
            f"{sar_m[i]:.9e} (rec {rec[lab][1]:.9e}, rel {r_cw[i]:.2e}) mass "
            f"{ccw[i]['mass_kg'] * 1e3:.6f} g{'   <- c*' if i == i_star else ''}"
        )
    _p0(
        f"(S10) max rel vs records: ccw {r_ccw.max():.3e}, cw(M) {r_cw.max():.3e} (band "
        f"{REUSE_REPRODUCTION_RTOL:g}, ASSERTED); c* = {key_star} (record {H9P_CSTAR_LABEL}), "
        f"SAR10g(c*) {sar_star:.9e} W/kg, rel {r_star:.3e}; z_c {z_c:.6e} m, az17 {az17:.9f} deg"
    )
    _p0(
        f"(iii-a) mirror at c* {mirror_star * 100:.4f}% (record {H9P_MIRROR_CSTAR_RECORD * 100:.4f}%); "
        f"(iii-b) lattice rel L2 {mirror_l2 * 100:.4f}% (record {H9P_MIRROR_L2_RECORD * 100:.4f}%) "
        f"(band {C4_COVARIANCE_BAND * 100:.1f}%, ASSERTED)"
    )
    _p0(
        f"price: B {t_b:.1f} s, 36 balls {t_sar:.1f} s; window "
        f"{time.perf_counter() - t_window:.1f} s at -n {comm.size}"
    )

    assert n_pts == H9P_N_SAMPLE_POINTS and n_mask == n_pts, (n_pts, n_mask)
    assert r_b <= REUSE_REPRODUCTION_RTOL, f"(B) {r_b:.3e}"
    assert sense_ratio > 1.0, f"(Sense) {sense_ratio:.4f}"
    assert mirror <= C4_COVARIANCE_BAND, f"(M_B) {mirror:.4e} — never re-band"
    assert mispaired > C4_COVARIANCE_BAND, f"(N1) blind comparator {mispaired:.4e}"
    assert r_ccw.max() <= REUSE_REPRODUCTION_RTOL, f"(S10) ccw {r_ccw.max():.3e}"
    assert r_cw.max() <= REUSE_REPRODUCTION_RTOL, f"(S10) cw {r_cw.max():.3e}"
    assert key_star == H9P_CSTAR_LABEL, key_star
    assert r_star <= REUSE_REPRODUCTION_RTOL, f"(S10) c* {r_star:.3e}"
    assert mirror_star <= C4_COVARIANCE_BAND, f"(iii-a) {mirror_star:.4e} — never re-band"
    assert mirror_l2 <= C4_COVARIANCE_BAND, f"(iii-b) {mirror_l2:.4e} — never re-band"
