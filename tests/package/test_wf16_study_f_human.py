"""`WF-16` step 8 (§9 item 355, chain U27): the port ``Study`` on F-human.

``Study(FHumanBirdcageGeometry(), gel, LumpedPortDrive(), frequency_hz=64e6,
degree=1)`` reproduces `WF-7` step 1's degree-1 F-human 32x32 (`GEO-35`
step 2a's leg-box-free fixture, re-recorded by §9 item 351) through the
package's matched-drive S route.

Anchors (asserted, every band imported, none set here):
  (W0) no mesh: the geometry's generator kwargs equal `WF-7`'s
       ``_f_human_kwargs()`` (and ``emit_leg_ports=True`` equals the boxed,
       no-key dict); every ring port's drive direction equals
       ``_ring_gap_frame(i, 16)[0]``; tags, conductor sigma and port impedance
       equal `WF-7`'s imports; the other Study methods and a
       ``LeontovichWall`` refuse with ``TypeError`` before any mesh.
  (F)  cells == ``WF7_CELL_RECORD`` and within ``CELL_COUNT_BAND`` of
       ``F_HUMAN_BRANCH_B_CELL_RECORD``; (N') leg census 0; 32 non-empty sheets.
  (R)  sigma_max and the 18 C16 x mirror class means against
       ``SIGMA_MAX_D1_RECORD`` / ``CLASS_MEAN_D1_RECORD`` at
       ``LEG_D0_REPRODUCTION_BAND`` -- at ``-n 8`` only (the record width).
  (G)  `WF-7`'s gates at their bands: reciprocity, passivity (sigma_max and
       column norms), C16 class count / azimuth residual / worst spread, the
       column-scale control lift.

Negative controls (rule (e)):
  (c1) asserted: the matrix misses the ``*_BOXED`` records by > the band
       (sigma_max and worst class mean each).
  (c2) ``FEM_EM_WF7_LEG_PORTS=boxed`` -> ``emit_leg_ports=True``: cells ==
       ``WF7_BOXED_CELL_RECORD``, the ``*_BOXED`` records reproduced at the
       band, the leg-box-free records missed by > the band (asserted).

Printed only: ``S_driven(P17)`` beside the step-0 record, the chord, the
entrywise max relative difference against item 351's saved
``wf7-s32-degree1-n8-noleg.npz`` if present, seconds, RSS.

Scope: degree 1, 64 MHz, sigma = 800 solid, one fixture, S only -- no fields,
so `WF-7`'s exact identity (v) and P17 ``P_src`` are out of scope.
"""

from __future__ import annotations

import os
import time
import unittest.mock

import numpy as np
import pytest
from mpi4py import MPI

from fem_em_solver import (
    FHumanBirdcageGeometry,
    LeontovichWall,
    LumpedPortDrive,
    QuasiStaticSphereField,
    RobinCooling,
    Study,
)
from fem_em_solver.materials import GelledSalinePhantomMaterial

from tests.complex_mode import complex_only
from tests.mesh.test_birdcage_leg_offset import _sheet_azimuth_deg
from tests.mesh.test_birdcage_port_sheet_prerequisite import CELL_COUNT_BAND
from tests.mesh.test_birdcage_port_sheets import PORT_LOWER, PORT_UPPER, SHEET_IFACE
from tests.mesh.test_birdcage_ring_gaps_scaleup import _ring_gap_frame
from tests.mesh.test_geo35_leg_port_boxes import _owned_census
from tests.validation.test_birdcage_f_human_rung import (
    F_HUMAN_BRANCH_B_CELL_RECORD,
    LEG_COUNT,
)
from tests.validation.test_lossy_sphere_fullwave import SALINE_EPSILON_R, SALINE_SIGMA
from tests.validation.test_port_birdcage_four_port import (
    LEG_D0_REPRODUCTION_BAND,
    TERMINATED_PORT_IMPEDANCE_OHM,
)
from tests.validation.test_port_birdcage_lumped_column import (
    CONDUCTOR_CELL_TAG,
    PHANTOM_CELL_TAG,
)
from tests.validation.test_port_birdcage_ring_column import (
    AZIMUTH_MATCH_DEG,
    CONTROL_COLUMN_SCALE,
    OPPOSITE_SPREAD_BAND,
    _reciprocity_ratio,
)
from tests.validation.test_port_birdcage_ring_matrix import (
    MATRIX_PASSIVITY_CEILING,
    N_SYMMETRY_CLASSES,
)
from tests.validation.test_port_gap_voltage_impedance import SIGMA_WIRE_S_PER_M
from tests.validation.test_port_lumped_sheet_sweep import RECIPROCITY_BAND
from tests.validation.test_wf7_f_human_b1_quadrature import _f_human_kwargs
from tests.validation.test_wf7_f_human_ring_matrix import (
    CLASS_MEAN_D1_RECORD,
    CLASS_MEAN_D1_RECORD_BOXED,
    CONTROL_RATIO_LIFT,
    DEGREE,
    FREQUENCY_HZ,
    S_DRIVEN_STEP0_RECORD,
    S_DRIVEN_STEP0_RECORD_BOXED,
    SIGMA_MAX_D1_RECORD,
    SIGMA_MAX_D1_RECORD_BOXED,
    WF7_BOXED_CELL_RECORD,
    WF7_CELL_RECORD,
    _boxed_knob,
    _class_means,
    _rss_summed_gib,
    _sheet_centre_z,
)

RECORD_RANK_WIDTH = 8  # the width item 351 recorded the 32x32 at
CHORD_PRINT_REFERENCE_M = 1.716154e-02  # WF-7 step 1's printed chord, printed beside
NPZ_NOLEG = "/workspace/logs/wf7/wf7-s32-degree1-n8-noleg.npz"
NPZ_BOXED = "/workspace/logs/wf7/wf7-s32-degree1-n8.npz"


def _phantom():
    return GelledSalinePhantomMaterial(SALINE_SIGMA, SALINE_EPSILON_R, FREQUENCY_HZ)


def _record_rel(s, az, zc, sigma_rec, class_rec):
    """``(sigma_max rel, worst class-mean rel, per-class rels, spreads, residual)``."""
    _classes, spreads, residual = _class_means(s, az, zc)
    sigma = float(np.linalg.norm(s, 2))
    rels = {k: abs(spreads[k][2] / class_rec[k] - 1.0) for k in class_rec}
    return abs(sigma / sigma_rec - 1.0), max(rels.values()), rels, spreads, residual


def test_wf16_f_human_w0_kwargs_tags_refusals():
    """(W0): no mesh anywhere in this node."""
    geom = FHumanBirdcageGeometry()
    assert geom.generator_kwargs() == _f_human_kwargs()
    assert FHumanBirdcageGeometry(emit_leg_ports=True).generator_kwargs() == \
        _f_human_kwargs(emit_leg_ports=True)
    assert "emit_leg_ports" not in FHumanBirdcageGeometry(
        emit_leg_ports=True).generator_kwargs()
    assert geom.leg_count == LEG_COUNT
    assert geom.ring_ordinals() == list(range(LEG_COUNT + 1, 3 * LEG_COUNT + 1))
    for i in geom.ring_ordinals():
        ref = np.asarray(_ring_gap_frame(i, LEG_COUNT)[0], dtype=float)
        assert np.array_equal(np.asarray(geom.drive_direction(i)), ref), i
    assert (geom.sheet_iface, geom.port_lower, geom.port_upper) == (
        SHEET_IFACE, PORT_LOWER, PORT_UPPER)
    assert geom.conductor_tag == CONDUCTOR_CELL_TAG
    assert geom.phantom_tag == PHANTOM_CELL_TAG
    assert geom.conductor_sigma_s_per_m == SIGMA_WIRE_S_PER_M
    assert LumpedPortDrive().port_impedance_ohm == TERMINATED_PORT_IMPEDANCE_OHM
    assert LumpedPortDrive().z0_ohm == TERMINATED_PORT_IMPEDANCE_OHM
    with pytest.raises(dataclasses_frozen_error()):
        geom.ring_radius = 0.2

    st = Study(geom, _phantom(), LumpedPortDrive(), frequency_hz=FREQUENCY_HZ,
               degree=DEGREE)
    with pytest.raises(TypeError):
        st.mean_sar()
    with pytest.raises(TypeError):
        st.temperature(RobinCooling(0.5, 10.0, 20.0))
    # estimate is priced since WF-13 step 10 (§9 item 562): the refusal is lifted; the branch reaches the build, nothing meshes
    with unittest.mock.patch.object(FHumanBirdcageGeometry, "build", side_effect=_BuildReached):
        with pytest.raises(_BuildReached):
            st.estimate()
    with pytest.raises(TypeError):
        st.sar(centres=[(0.0, 0.0, 0.0)])
    with pytest.raises(TypeError):
        st.b1_plus()
    with pytest.raises(TypeError):
        st.solve_sweep([FREQUENCY_HZ])
    with pytest.raises(TypeError):
        Study(geom, _phantom(), LumpedPortDrive(), frequency_hz=FREQUENCY_HZ,
              degree=DEGREE, conductor=LeontovichWall()).solve()
    with pytest.raises(TypeError):
        Study(geom, _phantom(), QuasiStaticSphereField(), frequency_hz=FREQUENCY_HZ,
              degree=DEGREE).solve()


class _BuildReached(Exception):
    """FH1: §9 item 562 — raised by the patched F-human build (nothing meshes)."""


def dataclasses_frozen_error():
    import dataclasses

    return dataclasses.FrozenInstanceError


@complex_only
def test_wf16_f_human_s_reproduces_wf7():
    """(F), (R), (G), (c1); (c2) under ``FEM_EM_WF7_LEG_PORTS=boxed``."""
    comm = MPI.COMM_WORLD
    if comm.size != RECORD_RANK_WIDTH:
        pytest.skip(f"(R) is asserted at -n {RECORD_RANK_WIDTH} only (the record "
                    f"width); this run is -n {comm.size}")
    boxed = _boxed_knob()

    def say(msg):
        if comm.rank == 0:
            print(f"[WF-16 s8] {msg}", flush=True)

    t_all = time.perf_counter()
    geom = FHumanBirdcageGeometry(emit_leg_ports=boxed)
    say(f"start: -n {comm.size}, f = {FREQUENCY_HZ:.3e} Hz, degree {DEGREE}, "
        f"route {'BOXED (c2)' if boxed else 'leg-box-free default'}")
    st = Study(geom, _phantom(), LumpedPortDrive(), frequency_hz=FREQUENCY_HZ,
               degree=DEGREE)
    res = st.solve(comm)
    say(f"solve: mesh {res.mesh_time:.2f} s, sweep {res.sweep_time:.2f} s, total "
        f"{time.perf_counter() - t_all:.2f} s; summed ru_maxrss "
        f"{_rss_summed_gib(comm):.3f} GiB")

    # ---- (F) ---------------------------------------------------------------
    msh, cell_tags, tags_f = res.mesh, res.cell_tags, res.facet_tags
    cells = int(res.cells)
    rel_b = cells / F_HUMAN_BRANCH_B_CELL_RECORD - 1.0
    leg_tags = ([PORT_LOWER + i for i in range(1, LEG_COUNT + 1)]
                + [PORT_UPPER + i for i in range(1, LEG_COUNT + 1)])
    census = sum(_owned_census(msh, cell_tags, leg_tags, comm).values())
    facets = [s["facets"] for s in res.sheets]
    cell_rec = WF7_BOXED_CELL_RECORD if boxed else WF7_CELL_RECORD
    say(f"(F) cells {cells} vs {'WF7_BOXED_CELL_RECORD' if boxed else 'WF7_CELL_RECORD'}"
        f" {cell_rec} (== ASSERTED); vs F_HUMAN_BRANCH_B {F_HUMAN_BRANCH_B_CELL_RECORD} "
        f"rel {rel_b:+.4e} (band {CELL_COUNT_BAND}, ASSERTED)")
    say(f"(N') leg census {census} ({'printed: boxed' if boxed else '== 0 ASSERTED'}); "
        f"sheets {len(facets)}, facets min {min(facets)} max {max(facets)} (> 0 ASSERTED)")
    chord = res.sheets[0]["h"]
    say(f"chord {chord:.6e} m beside {CHORD_PRINT_REFERENCE_M:.6e} m (printed); sheet "
        f"areas {min(s['area'] for s in res.sheets):.6e}..{max(s['area'] for s in res.sheets):.6e}")
    assert cells == cell_rec
    assert abs(rel_b) <= CELL_COUNT_BAND
    if not boxed:
        assert census == 0
    assert len(facets) == 2 * LEG_COUNT and min(facets) > 0

    # ---- the matrix, in WF-7's driven order (P17..P48) ----------------------
    s = np.asarray(res.s_matrix, dtype=complex)
    ids = [d.port_id for d in res.port_definitions]
    assert ids == [f"P{i}" for i in geom.ring_ordinals()]
    az = np.array([_sheet_azimuth_deg(msh, tags_f, SHEET_IFACE + i, comm)
                   for i in geom.ring_ordinals()])
    zc = np.array([_sheet_centre_z(msh, tags_f, SHEET_IFACE + i, comm)
                   for i in geom.ring_ordinals()])

    # ---- (G) ---------------------------------------------------------------
    ratio = _reciprocity_ratio(s)
    control = s.copy()
    control[:, 0] *= CONTROL_COLUMN_SCALE
    lift = _reciprocity_ratio(control) / max(ratio, 1e-300)
    sigma_max = float(np.linalg.norm(s, 2))
    col_norms = np.sum(np.abs(s) ** 2, axis=0)
    classes, spreads, residual = _class_means(s, az, zc)
    worst = max(spreads, key=lambda k: spreads[k][0])
    say(f"(G) reciprocity {ratio:.6e} (band {RECIPROCITY_BAND:.0e}); sigma_max "
        f"{sigma_max:.9f} margin {MATRIX_PASSIVITY_CEILING - sigma_max:+.3e}; column "
        f"norms {col_norms.min():.6f}..{col_norms.max():.6f} (ceiling "
        f"{MATRIX_PASSIVITY_CEILING}); classes {len(classes)} (== {N_SYMMETRY_CLASSES}); "
        f"azimuth residual {residual:.3e} deg (< {AZIMUTH_MATCH_DEG:.0e}); worst spread "
        f"{spreads[worst][0]:.4e} {worst} (band {OPPOSITE_SPREAD_BAND}); control lift "
        f"{lift:.3e}x (bar {CONTROL_RATIO_LIFT:.0e}x) -- all ASSERTED")

    # ---- (R) and the controls ------------------------------------------------
    own_sig, own_cls = ((SIGMA_MAX_D1_RECORD_BOXED, CLASS_MEAN_D1_RECORD_BOXED) if boxed
                        else (SIGMA_MAX_D1_RECORD, CLASS_MEAN_D1_RECORD))
    oth_sig, oth_cls = ((SIGMA_MAX_D1_RECORD, CLASS_MEAN_D1_RECORD) if boxed
                        else (SIGMA_MAX_D1_RECORD_BOXED, CLASS_MEAN_D1_RECORD_BOXED))
    r_sig, r_cls, r_rels, _sp, _r = _record_rel(s, az, zc, own_sig, own_cls)
    c_sig, c_cls, c_rels, _sp, _r = _record_rel(s, az, zc, oth_sig, oth_cls)
    lab = "BOXED records" if boxed else "leg-box-free records"
    say(f"(R) vs {lab}: sigma_max {sigma_max:.9f} vs {own_sig:.9f} rel {r_sig:.3e}; "
        f"worst class-mean rel {r_cls:.3e} (band {LEG_D0_REPRODUCTION_BAND:g}, ASSERTED)")
    for k in own_cls:
        say(f"     {'same ' if k[0] else 'other'} ring, {k[1]} steps  mean |S| "
            f"{spreads[k][2]:.12e} vs {own_cls[k]:.9e} rel {r_rels[k]:.3e}  | other "
            f"record rel {c_rels[k]:.3e}")
    say(f"({'c2' if boxed else 'c1'}) vs the other fixture's records: sigma_max rel "
        f"{c_sig:.3e} = {c_sig / LEG_D0_REPRODUCTION_BAND:.1f}x the band; worst class "
        f"rel {c_cls:.3e} = {c_cls / LEG_D0_REPRODUCTION_BAND:.3e}x (> 1x each ASSERTED)")
    s17 = complex(s[ids.index("P17"), ids.index("P17")])
    s_rec = S_DRIVEN_STEP0_RECORD_BOXED if boxed else S_DRIVEN_STEP0_RECORD
    say(f"S_driven(P17) {s17:.6f} beside {s_rec} (printed)")
    npz = NPZ_BOXED if boxed else NPZ_NOLEG
    if os.path.exists(npz):
        with np.load(npz, allow_pickle=False) as z:
            ref = z["s32"]
            order = [str(x) for x in z["order"]]
        same = order == ids
        if same:
            ent = float(np.max(np.abs(s - ref)) / np.max(np.abs(ref)))
            ent_rel = float(np.max(np.abs(s - ref) / np.abs(ref)))
            say(f"npz {npz}: order equal; max|dS|/max|S| {ent:.3e}, entrywise max rel "
                f"{ent_rel:.3e} (printed only)")
        else:
            say(f"npz {npz}: order differs (printed only)")
    else:
        say(f"npz {npz}: absent (printed)")
    say(f"window {time.perf_counter() - t_all:.2f} s")

    assert ratio <= RECIPROCITY_BAND
    assert sigma_max <= MATRIX_PASSIVITY_CEILING
    assert float(col_norms.max()) <= MATRIX_PASSIVITY_CEILING
    assert len(classes) == N_SYMMETRY_CLASSES
    assert residual < AZIMUTH_MATCH_DEG
    assert spreads[worst][0] <= OPPOSITE_SPREAD_BAND
    assert lift > CONTROL_RATIO_LIFT
    assert r_sig <= LEG_D0_REPRODUCTION_BAND
    assert r_cls <= LEG_D0_REPRODUCTION_BAND
    assert c_sig > LEG_D0_REPRODUCTION_BAND
    assert c_cls > LEG_D0_REPRODUCTION_BAND
