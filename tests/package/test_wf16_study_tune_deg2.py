"""`WF-16` step 6(c): ``Study.tune(k=1)`` at degree 2 reproduces `ANS-7`'s mode capacitance.

``Study.tune(k=1)`` on the exact-rectangle PEC hole F-small Study
(``FourPortBirdcageGeometry(as_hole=True, strip_fraction=GATED_WIDTH_FRACTION)``,
``conductor=None``, the Study's default degree 2, 64 MHz, ``-n 8``, κ off)
returns ``C_mode = 1/(omega Im z_1)``; it is compared with `ANS-7`'s
same-fixture degree-2 PEC reading within ``LEG_D0_REPRODUCTION_BAND``, and the
degree-1 `TH-17` record is the order control.

Anchors: (R0') the Study's class means equal the registered degree-2 ``"pec"``
column at ``LEG_D0_REPRODUCTION_BAND`` and its cell count equals the hole
record; (I) the API equals `TH-17`'s own modal route on the same S at
``EXACT_IDENTITY_RTOL``; (Im) ``Im z_1 > 0``; (a2) the comparand.  Negative
controls (rule (e); factors predicted from measured readings, only the side of
the band asserted): (c1) ``k = 2`` misses the comparand by more than
``REDUCTION_BAND``; (c-deg) the degree-1 `TH-17` record is missed by more than
``REDUCTION_BAND``.  Lossless, one mesh, no h-statement.  Every band is imported.
"""

from __future__ import annotations

import time

import numpy as np
import pytest
from mpi4py import MPI

from fem_em_solver.materials import GelledSalinePhantomMaterial
from fem_em_solver.ports.circuit import s_to_z
from fem_em_solver.study import LumpedPortDrive, Study
from fem_em_solver.utils.instrumentation import report_peak_rss
from tests.complex_mode import complex_only
from tests.package.test_wf16_study_copper import _geometry, _rel_means
from tests.package.test_wf16_study_tune import TH17_STEP3_C_MODE_RECORD_F
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
# ANS-7 step 7 (T21) k = 1 PEC mode capacitance: degree 2, -n 8, 80 798 cells,
# exact-rectangle PEC hole, kappa off, PEC column solved directly
# (20261004T140347Z_ANS-7-item217-W1-cmode.log:2086).  A comparand, never a
# band; its repr is exact, so there is no print floor.
ANS7_T21_C_MODE_PEC_P2_F = 6.361715963780996e-11
# Same log line: Im z_2 of the PEC column (printed .9e, floor 2e-10).  Used only
# for the (c1) predicted factor, never as a comparand.
ANS7_T21_Z2_PEC_P2_IMAG_OHM = 2.627824708e01
# Item 429's degree-1 API reading (20261014T004318Z_WF-16-s6-W1.log); printed only.
WF16_S6_DEG1_API_C_MODE_F = 6.626796421e-11


def _p0(*args):
    if MPI.COMM_WORLD.rank == 0:
        print(*args, flush=True)


@pytest.fixture(scope="module")
def tuned():
    comm = MPI.COMM_WORLD
    assert (
        comm.size == ANS6_EXACT_STRIP_P2_RECORD_RANK_WIDTH
    ), f"asserted at -n {ANS6_EXACT_STRIP_P2_RECORD_RANK_WIDTH}; got -n {comm.size}"
    phantom = GelledSalinePhantomMaterial(SALINE_SIGMA, SALINE_EPSILON_R, FREQUENCY_64_HZ)
    study = Study(_geometry(), phantom, LumpedPortDrive(), frequency_hz=FREQUENCY_64_HZ)
    assert study.degree == 2
    assert study.conductor is None
    t0 = time.perf_counter()
    res = study.solve(comm)
    t_solve = time.perf_counter() - t0
    t1 = study.tune(k=1, comm=comm)
    t2 = study.tune(k=2, comm=comm)
    assert study.result is res
    s = np.asarray(res.s_matrix, dtype=np.complex128)
    _p0(
        f"\n[WF-16 s6c] pec deg 2: cells {res.cells}, mesh {res.mesh_time:.1f} s, sweep "
        f"{res.sweep_time:.1f} s, total {t_solve:.1f} s, -n {comm.size}"
    )
    _p0("    z_k = " + ", ".join(f"{z.real:+.9e}{z.imag:+.9e}j" for z in t1.z_k))
    for k, t in ((1, t1), (2, t2)):
        c = t.c_mode_farad
        _p0(
            f"    C_{k} = {c:.15e} F; vs ANS-7 T21 {c / ANS7_T21_C_MODE_PEC_P2_F - 1:+.6e}; "
            f"vs TH-17 deg-1 record {c / TH17_STEP3_C_MODE_RECORD_F - 1:+.6e}"
        )
    _p0(
        f"    deg-2 / deg-1 API C_1 ratio {t1.c_mode_farad / WF16_S6_DEG1_API_C_MODE_F:.9f} "
        f"(move {t1.c_mode_farad / WF16_S6_DEG1_API_C_MODE_F - 1:+.6e}, predicted -4.0001 %, "
        f"printed)"
    )
    report_peak_rss(comm, label="[WF-16 s6c] summed peak ru_maxrss")
    return {"study": study, "res": res, "s": s, "means": _class_means(s), "t1": t1, "t2": t2}


@complex_only
def test_wf16_tune_p2_r0_fixture_tie(tuned):
    """(R0'): the degree-2 PEC Study is the registered exact-rectangle hole fixture."""
    assert tuned["res"].cells == ANS6_EXACT_STRIP_HOLE_CELL_RECORD
    rel = _rel_means(tuned["means"], ANS6_EXACT_STRIP_CLASS_MEANS_P2[KEY]["pec"])
    _p0(
        "\n[WF-16 s6c R0'] class-mean rel vs P2 'pec' record: "
        + ", ".join(f"{c} {rel[c]:.3e}" for c in CLASSES)
        + f" (band {LEG_D0_REPRODUCTION_BAND:.1e}; predicted <= 3e-10)"
    )
    for c in CLASSES:
        assert rel[c] <= LEG_D0_REPRODUCTION_BAND, (c, rel[c])


@complex_only
def test_wf16_tune_p2_i_api_equals_modal_route(tuned):
    """(I) and (Im): the API equals `TH-17`'s modal route on the same S."""
    t1 = tuned["t1"]
    zk, resid = _modal_z(s_to_z(tuned["s"], REFERENCE_IMPEDANCE_OHM))
    dz = float(np.max(np.abs(np.asarray(t1.z_k) - zk) / np.abs(zk)))
    dr = abs(t1.off_diagonal_residual - resid) / max(abs(resid), np.finfo(float).tiny)
    z1, z3 = t1.z_k[1], t1.z_k[3]
    _p0(
        f"\n[WF-16 s6c I] max rel |z_k(api) - z_k(TH-17)| = {dz:.3e}; residual rel "
        f"{dr:.3e} (api {t1.off_diagonal_residual:.3e}, item 217 read 4.550e-06); "
        f"|z1 - z3|/|z1| = {abs(z1 - z3) / abs(z1):.3e} (reciprocity, printed)"
    )
    assert dz <= EXACT_IDENTITY_RTOL
    assert dr <= EXACT_IDENTITY_RTOL
    assert t1.k == 1 and t1.frequency_hz == FREQUENCY_64_HZ
    assert z1.imag > 0.0
    omega = 2.0 * np.pi * FREQUENCY_64_HZ
    assert abs(t1.c_mode_farad * omega * zk[1].imag - 1.0) <= EXACT_IDENTITY_RTOL


@complex_only
def test_wf16_tune_p2_a_c_mode_equals_ans7_t21(tuned):
    """(a2): C_mode within LEG_D0_REPRODUCTION_BAND of `ANS-7`'s degree-2 reading."""
    c = tuned["t1"].c_mode_farad
    gap = c / ANS7_T21_C_MODE_PEC_P2_F - 1.0
    _p0(
        f"\n[WF-16 s6c a2] C_mode = {c:.15e} F; comparand {ANS7_T21_C_MODE_PEC_P2_F!r} F; "
        f"gap {gap:+.3e} (predicted <= 3e-10, printed); band {LEG_D0_REPRODUCTION_BAND:.1e}"
    )
    assert abs(gap) <= LEG_D0_REPRODUCTION_BAND


@complex_only
def test_wf16_tune_p2_c1_mode_index_control(tuned):
    """(c1): k = 2 on the same Study misses the comparand (predicted +48.754 %)."""
    c2 = tuned["t2"].c_mode_farad
    gap = c2 / ANS7_T21_C_MODE_PEC_P2_F - 1.0
    pred = 1.0 / (2.0 * np.pi * FREQUENCY_64_HZ * ANS7_T21_Z2_PEC_P2_IMAG_OHM)
    _p0(
        f"\n[WF-16 s6c c1] k = 2 C = {c2:.9e} F, gap {gap:+.6e} "
        f"({abs(gap) / LEG_D0_REPRODUCTION_BAND:.3e}x LEG_D0, "
        f"{abs(gap) / REDUCTION_BAND:.1f}x REDUCTION_BAND); predicted from T21 z_2 "
        f"{pred:.6e} F ({pred / ANS7_T21_C_MODE_PEC_P2_F - 1:+.6e})"
    )
    assert abs(gap) > LEG_D0_REPRODUCTION_BAND
    assert abs(gap) > REDUCTION_BAND


@complex_only
def test_wf16_tune_p2_cdeg_order_control(tuned):
    """(c-deg): the degree-2 C_mode misses `TH-17`'s degree-1 record (predicted -4.0287 %)."""
    c = tuned["t1"].c_mode_farad
    gap = c / TH17_STEP3_C_MODE_RECORD_F - 1.0
    _p0(
        f"\n[WF-16 s6c c-deg] C_1 vs TH-17 deg-1 record gap {gap:+.6e} "
        f"({abs(gap) / REDUCTION_BAND:.1f}x band; predicted -4.0287 %)"
    )
    assert abs(gap) > REDUCTION_BAND
