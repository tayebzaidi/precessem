"""`WF-16` step 6(a): ``Study.tune(k=1)`` reproduces `TH-17`'s mode capacitance.

``Study.tune(k=1)`` on the exact-rectangle PEC hole F-small Study
(``FourPortBirdcageGeometry(as_hole=True, strip_fraction=GATED_WIDTH_FRACTION)``,
``conductor=None``, degree 1, 64 MHz, ``-n 2``, κ off) returns
``C_mode = 1/(omega Im z_1)`` with ``z_k = diag(F Z F^H)``; it is compared with
`TH-17`'s step-3 record (degree 1, ragged-sheet PEC hole) within
``REDUCTION_BAND``.

Anchors: (R0) the Study's class means equal the registered ``"pec"`` column at
``LEG_D0_REPRODUCTION_BAND`` and its cell count equals the hole record; (I) the
API equals `TH-17`'s own modal route on the same S at ``EXACT_IDENTITY_RTOL``;
(Im) ``Im z_1 > 0``; (a) the record comparison.  Negative controls (rule (e);
factors predicted from registered records, only the sign asserted): (c1)
``k = 2`` and (c2) the solid sigma = 800 F-small Study, each missing the record
by more than ``REDUCTION_BAND``; (c3) refusals with no mesh.  Lossless: a
uniform capacitor ESR shifts ``Re z_k`` only.  Every band is imported.
"""

from __future__ import annotations

import time
from types import SimpleNamespace

import numpy as np
import pytest
from mpi4py import MPI

import fem_em_solver.study as study_mod
from fem_em_solver.materials import GelledSalinePhantomMaterial
from fem_em_solver.ports.circuit import s_to_z
from fem_em_solver.study import (
    CapacitorSheets,
    FHumanBirdcageGeometry,
    FourPortBirdcageGeometry,
    LeontovichWall,
    LumpedPortDrive,
    QuasiStaticSphereField,
    SphereInBoxGeometry,
    Study,
    TwoPortQuadratureDrive,
)

from tests.complex_mode import complex_only
from tests.package.test_wf16_study_copper import _rel_means
from tests.package.test_wf16_study_four_port import _record, _rel
from tests.package.test_wf1_yaml_copper import _py_copper
from tests.validation.test_ans6_exact_strip import (
    ANS6_EXACT_STRIP_CLASS_MEANS_P1,
    ANS6_EXACT_STRIP_HOLE_CELL_RECORD,
    ANS6_EXACT_STRIP_RECORD_RANK_WIDTH,
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
# TH-17 step 3's k = 1 mode capacitance, degree 1, item 51's ragged-sheet PEC
# hole mesh (TH-17.md step 3 record; 20260928T095346Z_TH-17-step3.log.gz,
# -n 4).  A comparand, never a band; its .9e print floor (<= 5e-10 relative)
# sits far under REDUCTION_BAND.
TH17_STEP3_C_MODE_RECORD_F = 6.628769508e-11


def _p0(*args):
    if MPI.COMM_WORLD.rank == 0:
        print(*args, flush=True)


def _pec_circulant_s():
    m = ANS6_EXACT_STRIP_CLASS_MEANS_P1[KEY]["pec"]
    cls = (m["self"], m["adjacent"], m["opposite"], m["adjacent"])
    return np.array([[cls[(j - i) % 4] for j in range(4)] for i in range(4)],
                    dtype=np.complex128)


def _solid_study():
    phantom = GelledSalinePhantomMaterial(SALINE_SIGMA, SALINE_EPSILON_R, FREQUENCY_64_HZ)
    return Study(FourPortBirdcageGeometry.f_small(), phantom, LumpedPortDrive(),
                 frequency_hz=FREQUENCY_64_HZ, degree=1)


def test_wf16_tune_refusals(monkeypatch):
    """(c3): refusals by construction, no mesh is built."""
    calls = {"n": 0}

    def stub(self, comm=MPI.COMM_WORLD):
        calls["n"] += 1
        raise AssertionError("build must not be called by a refused tune")

    for cls in (SphereInBoxGeometry, FHumanBirdcageGeometry, FourPortBirdcageGeometry):
        monkeypatch.setattr(cls, "build", stub)
    phantom = GelledSalinePhantomMaterial(SALINE_SIGMA, SALINE_EPSILON_R, FREQUENCY_64_HZ)
    sphere = Study(SphereInBoxGeometry(0.05, 0.15, 0.01, 0.03), phantom,
                   QuasiStaticSphereField(), frequency_hz=FREQUENCY_64_HZ, degree=1)
    with pytest.raises(TypeError):
        sphere.tune()
    human = Study(FHumanBirdcageGeometry(as_hole=True), phantom,
                  TwoPortQuadratureDrive(118.0e-12, esr_ohm=0.2),
                  frequency_hz=FREQUENCY_64_HZ, degree=1,
                  conductor=LeontovichWall(study_mod.COPPER_SIGMA))
    with pytest.raises(TypeError):
        human.tune()
    caps = Study(FourPortBirdcageGeometry.f_small(), phantom,
                 CapacitorSheets(15.57e-12, on=("P2",)), frequency_hz=FREQUENCY_64_HZ,
                 degree=1)
    with pytest.raises(TypeError):
        caps.tune()
    pec = _py_copper("pec")
    assert pec.degree == 1 and pec.conductor is None
    with pytest.raises(ValueError):
        pec.tune(k=0)  # no result yet: refused before any solve
    for bad in (3, -1, 4):
        with pytest.raises(ValueError):
            pec.tune(k=bad)
    with pytest.raises(TypeError):
        pec.tune(k=1.0)
    assert calls["n"] == 0 and pec.result is None
    pec.result = SimpleNamespace(s_matrix=_pec_circulant_s())
    with pytest.raises(ValueError):
        pec.tune(k=0)  # Im z_0 < 0 on the PEC circulant
    t1 = pec.tune(k=1)
    zk, _ = _modal_z(s_to_z(_pec_circulant_s(), REFERENCE_IMPEDANCE_OHM))
    _p0(f"\n[WF-16 s6 W0] circulant z_k = {zk}; circulant C_1 = {t1.c_mode_farad:.9e} F, "
        f"vs record {t1.c_mode_farad / TH17_STEP3_C_MODE_RECORD_F - 1:+.3e} (printed)")
    assert zk[0].imag < 0.0
    assert calls["n"] == 0


@pytest.fixture(scope="module")
def tuned():
    comm = MPI.COMM_WORLD
    assert comm.size == ANS6_EXACT_STRIP_RECORD_RANK_WIDTH, (
        f"asserted at -n {ANS6_EXACT_STRIP_RECORD_RANK_WIDTH}; got -n {comm.size}")
    out = {}
    for name, make in (("pec", lambda: _py_copper("pec")), ("solid", _solid_study)):
        study = make()
        assert study.degree == 1
        t0 = time.perf_counter()
        res = study.solve(comm)
        t_solve = time.perf_counter() - t0
        s = np.asarray(res.s_matrix, dtype=np.complex128)
        out[name] = {"study": study, "res": res, "s": s, "means": _class_means(s),
                     "t1": study.tune(k=1, comm=comm), "t2": study.tune(k=2, comm=comm),
                     "t": t_solve}
        _p0(f"\n[WF-16 s6] {name}: cells {res.cells}, mesh {res.mesh_time:.1f} s, sweep "
            f"{res.sweep_time:.1f} s, total {t_solve:.1f} s, -n {comm.size}")
        _p0(f"    {name} z_k = " + ", ".join(
            f"{z.real:+.6e}{z.imag:+.6e}j" for z in out[name]["t1"].z_k))
        for k in (1, 2):
            c = out[name][f"t{k}"].c_mode_farad
            _p0(f"    {name} C_{k} = {c:.9e} F, vs TH-17 record "
                f"{c / TH17_STEP3_C_MODE_RECORD_F - 1:+.6e}")
    return out


@complex_only
def test_wf16_tune_r0_fixture_tie(tuned):
    """(R0): the PEC Study is the registered exact-rectangle hole fixture."""
    pec = tuned["pec"]
    assert pec["res"].cells == ANS6_EXACT_STRIP_HOLE_CELL_RECORD
    rel = _rel_means(pec["means"], ANS6_EXACT_STRIP_CLASS_MEANS_P1[KEY]["pec"])
    _p0(f"\n[WF-16 s6 R0] class-mean rel vs 'pec' record: "
        + ", ".join(f"{c} {rel[c]:.3e}" for c in CLASSES))
    for c in CLASSES:
        assert rel[c] <= LEG_D0_REPRODUCTION_BAND, (c, rel[c])


@complex_only
def test_wf16_tune_i_api_equals_modal_route(tuned):
    """(I) and (Im): the API equals `TH-17`'s modal route on the same S."""
    pec = tuned["pec"]
    t1 = pec["t1"]
    zk, resid = _modal_z(s_to_z(pec["s"], REFERENCE_IMPEDANCE_OHM))
    dz = float(np.max(np.abs(np.asarray(t1.z_k) - zk) / np.abs(zk)))
    dr = abs(t1.off_diagonal_residual - resid) / max(abs(resid), np.finfo(float).tiny)
    z1, z3 = t1.z_k[1], t1.z_k[3]
    _p0(f"\n[WF-16 s6 I] max rel |z_k(api) - z_k(TH-17)| = {dz:.3e}; residual rel "
        f"{dr:.3e} (api {t1.off_diagonal_residual:.3e}); |z1 - z3|/|z1| = "
        f"{abs(z1 - z3) / abs(z1):.3e} (reciprocity, printed)")
    assert dz <= EXACT_IDENTITY_RTOL
    assert dr <= EXACT_IDENTITY_RTOL
    assert t1.k == 1 and t1.frequency_hz == FREQUENCY_64_HZ
    assert z1.imag > 0.0
    omega = 2.0 * np.pi * FREQUENCY_64_HZ
    assert abs(t1.c_mode_farad * omega * zk[1].imag - 1.0) <= EXACT_IDENTITY_RTOL


@complex_only
def test_wf16_tune_a_c_mode_equals_th17_record(tuned):
    """(a): C_mode on the exact PEC hole within REDUCTION_BAND of `TH-17`'s."""
    c = tuned["pec"]["t1"].c_mode_farad
    gap = c / TH17_STEP3_C_MODE_RECORD_F - 1.0
    _p0(f"\n[WF-16 s6 a] C_mode = {c:.9e} F; record {TH17_STEP3_C_MODE_RECORD_F:.9e} F; "
        f"gap {gap:+.6e} (predicted -2.98e-4 in [-4.5e-4, -1.5e-4], printed); "
        f"band {REDUCTION_BAND:.1e}")
    assert abs(gap) <= REDUCTION_BAND


@complex_only
def test_wf16_tune_c1_mode_index_control(tuned):
    """(c1): k = 2 on the same Study misses the record (predicted +48.9 %)."""
    c2 = tuned["pec"]["t2"].c_mode_farad
    gap = c2 / TH17_STEP3_C_MODE_RECORD_F - 1.0
    _p0(f"\n[WF-16 s6 c1] k = 2 C = {c2:.9e} F, gap {gap:+.6e} "
        f"({abs(gap) / REDUCTION_BAND:.1f}x band; predicted +48.9 %)")
    assert abs(gap) > REDUCTION_BAND


@complex_only
def test_wf16_tune_c2_wrong_fixture_control(tuned):
    """(c2): the solid sigma = 800 Study ties to its record, then misses `TH-17`'s."""
    solid = tuned["solid"]
    tie = _rel(solid["s"], _record(1))
    c = solid["t1"].c_mode_farad
    gap = c / TH17_STEP3_C_MODE_RECORD_F - 1.0
    _p0(f"\n[WF-16 s6 c2] solid S vs record rel {tie:.3e}; C_1 = {c:.9e} F, gap "
        f"{gap:+.6e} ({abs(gap) / REDUCTION_BAND:.1f}x band; predicted -6.98 %)")
    assert tie <= LEG_D0_REPRODUCTION_BAND
    assert abs(gap) > REDUCTION_BAND
