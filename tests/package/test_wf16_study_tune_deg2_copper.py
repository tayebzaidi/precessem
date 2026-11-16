"""`WF-16` step 6(d): ``Study.tune(k=1)`` at degree 2 on copper reproduces `ANS-7`'s copper mode capacitance.

``Study.tune(k=1)`` on the exact-rectangle copper (Leontovich) hole F-small Study
(``FourPortBirdcageGeometry(as_hole=True, strip_fraction=GATED_WIDTH_FRACTION)``,
``conductor=LeontovichWall(COPPER_SIGMA)``, the Study's default degree 2, 64 MHz,
``-n 8``, κ off) returns `TH-17`'s ``C_mode = 1/(omega Im z_1)``; it is compared
with `ANS-7`'s same-fixture degree-2 copper reading within
``LEG_D0_REPRODUCTION_BAND``, and the class means are tied to `ANS-6`'s degree-2
``"cu"`` record.

Anchors: (R0'') the Study's class means equal the registered degree-2 ``"cu"``
column at ``LEG_D0_REPRODUCTION_BAND`` and its cell count equals the hole
record; (I) the API equals `TH-17`'s own modal route on the same S at
``EXACT_IDENTITY_RTOL``; (Im) ``Im z_1 > 0``; (a3) the comparand.  Negative
controls (rule (e); each the same comparison already measured on the same
fixture, only the side of the band asserted): (c1) ``k = 2`` misses the copper
comparand by more than ``REDUCTION_BAND``; (c-pec) the copper Study misses the
degree-2 PEC comparands (``C_mode`` and class means) by more than
``NEGATIVE_CONTROL_FLOOR``.  Lossless capacitor, one mesh, no h-statement.
Every band is imported.
"""

from __future__ import annotations

import time

import numpy as np
import pytest
from mpi4py import MPI

import fem_em_solver.study as study_mod
from fem_em_solver import LeontovichWall
from fem_em_solver.materials import GelledSalinePhantomMaterial
from fem_em_solver.ports.circuit import s_to_z
from fem_em_solver.study import LumpedPortDrive, Study
from fem_em_solver.utils.instrumentation import report_peak_rss

from tests.complex_mode import complex_only
from tests.package.test_wf16_study_copper import NEGATIVE_CONTROL_FLOOR, _geometry, _rel_means
from tests.package.test_wf16_study_tune import TH17_STEP3_C_MODE_RECORD_F
from tests.package.test_wf16_study_tune_deg2 import ANS7_T21_C_MODE_PEC_P2_F
from tests.validation import test_th14_birdcage_copper as th14
from tests.validation.test_ans6_exact_strip import (
    ANS6_EXACT_STRIP_CLASS_MEANS_P2,
    ANS6_EXACT_STRIP_HOLE_CELL_RECORD,
    ANS6_EXACT_STRIP_P2_RECORD_RANK_WIDTH,
    CLASSES,
    _class_means,
)
from tests.validation.test_birdcage_sar_mass_averaged import EXACT_IDENTITY_RTOL
from tests.validation.test_lossy_sphere_fullwave import (
    FREQUENCY_64_HZ,
    SALINE_EPSILON_R,
    SALINE_SIGMA,
)
from tests.validation.test_port_birdcage_four_port import LEG_D0_REPRODUCTION_BAND
from tests.validation.test_port_lumped_rlc_termination import REDUCTION_BAND
from tests.validation.test_port_package_sparameters import REFERENCE_IMPEDANCE_OHM
from tests.validation.test_th17_birdcage_eigenmodes import _modal_z

KEY = "64"
# ANS-7 step 7 (T21) k = 1 copper mode capacitance: degree 2, -n 8, 80 798
# cells, exact-rectangle hole, kappa off, Cu sigma = 5.8e7 S/m solved directly
# (20261004T140347Z_ANS-7-item217-W1-cmode.log:2084, :2088; the json MAT-15
# step 3 reads carries the same repr, 20261005T040135Z_MAT-15-step3-W1.log:215).
# A comparand, never a band; its repr is exact, so there is no print floor.
ANS7_T21_C_MODE_CU_P2_F = 6.358914053477356e-11
# Same log line (:2084), same fixture: Im z_2 of the Cu column (printed .9e,
# floor 2e-10).  Used only for the (c1) predicted factor, never as a comparand.
ANS7_T21_Z2_CU_P2_IMAG_OHM = 2.628994510e+01
# WF-16 step 6(c)'s API PEC reading, same fixture, degree 2, -n 8
# (20261014T111026Z_WF-16-s6c-W1.log:1013, :1036).  Printed only.
WF16_S6C_API_C_MODE_PEC_F = 6.361715963781403e-11


def _p0(*args):
    if MPI.COMM_WORLD.rank == 0:
        print(*args, flush=True)


@pytest.fixture(scope="module")
def tuned():
    comm = MPI.COMM_WORLD
    assert comm.size == ANS6_EXACT_STRIP_P2_RECORD_RANK_WIDTH, (
        f"asserted at -n {ANS6_EXACT_STRIP_P2_RECORD_RANK_WIDTH}; got -n {comm.size}")
    phantom = GelledSalinePhantomMaterial(SALINE_SIGMA, SALINE_EPSILON_R, FREQUENCY_64_HZ)
    study = Study(_geometry(), phantom, LumpedPortDrive(), frequency_hz=FREQUENCY_64_HZ,
                  conductor=LeontovichWall(study_mod.COPPER_SIGMA))
    assert study.degree == 2
    assert isinstance(study.conductor, LeontovichWall)
    assert study.conductor.sigma_s_per_m == th14.COPPER_SIGMA
    t0 = time.perf_counter()
    res = study.solve(comm)
    t_solve = time.perf_counter() - t0
    t1 = study.tune(k=1, comm=comm)
    t2 = study.tune(k=2, comm=comm)
    assert study.result is res
    s = np.asarray(res.s_matrix, dtype=np.complex128)
    _p0(f"\n[WF-16 s6d] cu deg 2: cells {res.cells}, mesh {res.mesh_time:.1f} s, sweep "
        f"{res.sweep_time:.1f} s, total {t_solve:.1f} s, -n {comm.size}")
    _p0("    z_k = " + ", ".join(f"{z.real:+.9e}{z.imag:+.9e}j" for z in t1.z_k))
    for k, t in ((1, t1), (2, t2)):
        c = t.c_mode_farad
        _p0(f"    C_{k} = {c:.15e} F; vs ANS-7 T21 Cu {c / ANS7_T21_C_MODE_CU_P2_F - 1:+.6e}; "
            f"vs ANS-7 T21 PEC {c / ANS7_T21_C_MODE_PEC_P2_F - 1:+.6e}; "
            f"vs WF-16 s6c API PEC {c / WF16_S6C_API_C_MODE_PEC_F - 1:+.6e}; "
            f"vs TH-17 deg-1 record {c / TH17_STEP3_C_MODE_RECORD_F - 1:+.6e}")
    report_peak_rss(comm, label="[WF-16 s6d] summed peak ru_maxrss")
    return {"study": study, "res": res, "s": s, "means": _class_means(s), "t1": t1, "t2": t2}


@complex_only
def test_wf16_tune_p2cu_r0_fixture_tie(tuned):
    """(R0''): the degree-2 copper Study is the registered exact-rectangle hole fixture."""
    assert tuned["res"].cells == ANS6_EXACT_STRIP_HOLE_CELL_RECORD
    rel = _rel_means(tuned["means"], ANS6_EXACT_STRIP_CLASS_MEANS_P2[KEY]["cu"])
    _p0(f"\n[WF-16 s6d R0''] class-mean rel vs P2 'cu' record: "
        + ", ".join(f"{c} {rel[c]:.3e}" for c in CLASSES)
        + f" (band {LEG_D0_REPRODUCTION_BAND:.1e}; predicted <= 3e-10)")
    for c in CLASSES:
        assert rel[c] <= LEG_D0_REPRODUCTION_BAND, (c, rel[c])


@complex_only
def test_wf16_tune_p2cu_i_api_equals_modal_route(tuned):
    """(I) and (Im): the API equals `TH-17`'s modal route on the same S."""
    t1 = tuned["t1"]
    zk, resid = _modal_z(s_to_z(tuned["s"], REFERENCE_IMPEDANCE_OHM))
    dz = float(np.max(np.abs(np.asarray(t1.z_k) - zk) / np.abs(zk)))
    dr = abs(t1.off_diagonal_residual - resid) / max(abs(resid), np.finfo(float).tiny)
    z1, z3 = t1.z_k[1], t1.z_k[3]
    _p0(f"\n[WF-16 s6d I] max rel |z_k(api) - z_k(TH-17)| = {dz:.3e}; residual rel "
        f"{dr:.3e} (api {t1.off_diagonal_residual:.3e}, T21 read 4.561e-06); "
        f"|z1 - z3|/|z1| = {abs(z1 - z3) / abs(z1):.3e} (reciprocity, printed; "
        f"T21 read 8.369e-16)")
    assert dz <= EXACT_IDENTITY_RTOL
    assert dr <= EXACT_IDENTITY_RTOL
    assert t1.k == 1 and t1.frequency_hz == FREQUENCY_64_HZ
    assert z1.imag > 0.0
    omega = 2.0 * np.pi * FREQUENCY_64_HZ
    assert abs(t1.c_mode_farad * omega * zk[1].imag - 1.0) <= EXACT_IDENTITY_RTOL


@complex_only
def test_wf16_tune_p2cu_a_c_mode_equals_ans7_t21_cu(tuned):
    """(a3): C_mode within LEG_D0_REPRODUCTION_BAND of `ANS-7`'s degree-2 copper reading."""
    c = tuned["t1"].c_mode_farad
    gap = c / ANS7_T21_C_MODE_CU_P2_F - 1.0
    _p0(f"\n[WF-16 s6d a3] C_mode = {c!r} F; comparand {ANS7_T21_C_MODE_CU_P2_F!r} F; "
        f"gap {gap:+.3e} (predicted <= 3e-10, printed); band {LEG_D0_REPRODUCTION_BAND:.1e}")
    assert abs(gap) <= LEG_D0_REPRODUCTION_BAND


@complex_only
def test_wf16_tune_p2cu_c1_mode_index_control(tuned):
    """(c1): k = 2 on the same Study misses the copper comparand (predicted +48.7536 %)."""
    c2 = tuned["t2"].c_mode_farad
    gap = c2 / ANS7_T21_C_MODE_CU_P2_F - 1.0
    pred = 1.0 / (2.0 * np.pi * FREQUENCY_64_HZ * ANS7_T21_Z2_CU_P2_IMAG_OHM)
    _p0(f"\n[WF-16 s6d c1] k = 2 C = {c2:.9e} F, gap {gap:+.6e} "
        f"({abs(gap) / LEG_D0_REPRODUCTION_BAND:.3e}x LEG_D0, "
        f"{abs(gap) / REDUCTION_BAND:.1f}x REDUCTION_BAND); predicted from T21 Cu z_2 "
        f"{pred:.6e} F ({pred / ANS7_T21_C_MODE_CU_P2_F - 1:+.6e})")
    assert abs(gap) > LEG_D0_REPRODUCTION_BAND
    assert abs(gap) > REDUCTION_BAND


@complex_only
def test_wf16_tune_p2cu_cpec_conductor_control(tuned):
    """(c-pec): the copper Study misses the degree-2 PEC comparands (predicted -4.4043e-4)."""
    c = tuned["t1"].c_mode_farad
    gap = c / ANS7_T21_C_MODE_PEC_P2_F - 1.0
    gap_api = c / WF16_S6C_API_C_MODE_PEC_F - 1.0
    rel = _rel_means(tuned["means"], ANS6_EXACT_STRIP_CLASS_MEANS_P2[KEY]["pec"])
    _p0(f"\n[WF-16 s6d c-pec] C_1 vs ANS-7 T21 PEC gap {gap:+.6e} "
        f"({abs(gap) / LEG_D0_REPRODUCTION_BAND:.3e}x LEG_D0, "
        f"{abs(gap) / NEGATIVE_CONTROL_FLOOR:.1f}x floor; predicted -4.4043e-4); "
        f"vs WF-16 s6c API PEC {gap_api:+.6e} (printed); class-mean rel vs P2 'pec' record: "
        + ", ".join(f"{k} {rel[k]:.3e}" for k in CLASSES)
        + f" (floor {NEGATIVE_CONTROL_FLOOR:.1e})")
    assert abs(gap) > LEG_D0_REPRODUCTION_BAND
    assert abs(gap) > NEGATIVE_CONTROL_FLOOR
    for k in CLASSES:
        assert rel[k] > NEGATIVE_CONTROL_FLOOR, (k, rel[k])
