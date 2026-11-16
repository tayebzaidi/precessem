"""`WF-16` step 6(b): ``Study.tune(k=1)``'s ``C_mode`` placed on the gap sheets.

On the degree-1 exact-rectangle PEC hole F-small Study
(``FourPortBirdcageGeometry(as_hole=True, strip_fraction=GATED_WIDTH_FRACTION)``,
``conductor=None``, 64 MHz, ``-n 2``, κ off), ``CapacitorSheets(C_mode, esr_ohm=0)``
with ``C_mode`` from ``Study.tune(k=1)`` reproduces the circuit reduction of the
same Study's untuned 4x4 at ``REDUCTION_BAND``: Route A (`PORT-22`'s
construction, P1 kept, P2-P4 at ``C_mode``) and per sheet (`PORT-24`'s
construction).  The tuned modal ``Im z_1`` of the circuit-side ``Z + z_cap I``
is zero **by construction** (``z_1 + 1/(j omega C_mode)`` is `TH-17`'s loop
impedance at the C that ``tune`` defined from that z_1): it is printed as a
record, never a gate (2026-10-15 03:00 review ruling).  The identity holds at
any C, so it is not a tuning statement; the k = 2 controls show its comparator
reads C.

Anchors: (R0) the Study's cells and class means equal the ``"pec"`` hole record
at ``LEG_D0_REPRODUCTION_BAND``; (a') ``C_mode`` within ``REDUCTION_BAND`` of
`TH-17`'s record.  Negative controls (rule (e); factors predicted from item
429's printed z_k, only the band exceedance asserted): (c1) the k = 2
capacitance leaves ``Im z_1`` of order 0.33 |z_1|; (c2) Route A's S11 against
the reduction at the k = 2 capacitance.  Lossless only.  Every band is
imported; none is set here.  Mesh economy (copied from item 373's fixture,
``tests/package/test_wf16_study_capacitor_sheets.py``, accepted by the
2026-10-11 review): ``FourPortBirdcageGeometry.build`` is wrapped so the first
call builds the mesh and later equal-geometry calls return it; counts asserted.
"""

from __future__ import annotations

import time

import numpy as np
import pytest
from mpi4py import MPI

from fem_em_solver import (
    CapacitorSheets,
    FourPortBirdcageGeometry,
    LumpedPortDrive,
    Study,
)
from fem_em_solver.materials import GelledSalinePhantomMaterial
from fem_em_solver.ports.circuit import reduce_terminated_ports, s_to_z

from tests.complex_mode import complex_only
from tests.package.test_wf16_study_copper import _geometry
from tests.package.test_wf16_study_tune import TH17_STEP3_C_MODE_RECORD_F
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
from tests.validation.test_th17_birdcage_eigenmodes import _loop_z, _modal_z

PORTS = ("P1", "P2", "P3", "P4")
# Printed beside, never asserted.
ITEM429_C_MODE_F = 6.626796421e-11  # 20261014T004318Z_WF-16-s6-W1.log
ITEM429_RE_Z1_OHM = 0.02500348  # same log, z_k print
PRIOR_ROUTE_A_COPPER = 1.406e-05  # item 373, copper, C_F, r = 0.2
PRIOR_ROUTE_A_S11_COPPER = complex(-0.2760, 0.9536)
PRIOR_I4_COPPER_R0 = (3.096783e-05, 3.083764e-05, 3.088213e-05, 3.050839e-05)  # PORT-24 r = 0


def _p0(*args):
    if MPI.COMM_WORLD.rank == 0:
        print(*args, flush=True)


def _phantom():
    return GelledSalinePhantomMaterial(SALINE_SIGMA, SALINE_EPSILON_R, FREQUENCY_64_HZ)


def _study(excitation):
    return Study(_geometry(), _phantom(), excitation, frequency_hz=FREQUENCY_64_HZ, degree=1)


def _relf(a, b):
    a = np.asarray(a, dtype=np.complex128)
    b = np.asarray(b, dtype=np.complex128)
    return float(np.linalg.norm(a - b) / np.linalg.norm(b))


@pytest.fixture(scope="module")
def s6b():
    comm = MPI.COMM_WORLD
    assert comm.size == ANS6_EXACT_STRIP_RECORD_RANK_WIDTH, (
        f"asserted at -n {ANS6_EXACT_STRIP_RECORD_RANK_WIDTH}; got -n {comm.size}")
    real_build = FourPortBirdcageGeometry.build
    cache = {"key": None, "built": None, "hits": 0, "real": 0}

    def _build(self, comm_):
        if cache["built"] is not None and cache["key"] == self and comm_ is comm:
            cache["hits"] += 1
            return cache["built"]
        cache["real"] += 1
        out = real_build(self, comm_)
        cache["key"], cache["built"] = self, out
        return out

    out = {}
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(FourPortBirdcageGeometry, "build", _build)
        s0 = _study(LumpedPortDrive())
        assert s0.degree == 1 and s0.conductor is None
        t0 = time.perf_counter()
        res0 = s0.solve(comm)
        out["S0"] = {"study": s0, "res": res0, "t": time.perf_counter() - t0}
        t1 = s0.tune(k=1, comm=comm)
        t2 = s0.tune(k=2, comm=comm)
        assert s0.result is res0
        out["t1"], out["t2"] = t1, t2
        _p0(f"\n[WF-16 s6b] S0: cells {res0.cells}, {out['S0']['t']:.1f} s; z_k = " + ", ".join(
            f"{z.real:+.9e}{z.imag:+.9e}j" for z in t1.z_k))
        _p0(f"[WF-16 s6b] C_1 = {t1.c_mode_farad:.9e} F, C_2 = {t2.c_mode_farad:.9e} F")
        plan = [("S_A", CapacitorSheets(t1.c_mode_farad, 0.0, on=("P2", "P3", "P4")))]
        plan += [(f"S{i + 1}", CapacitorSheets(t1.c_mode_farad, 0.0, on=(PORTS[i],)))
                 for i in range(4)]
        plan += [("S_A2", CapacitorSheets(t2.c_mode_farad, 0.0, on=("P2", "P3", "P4")))]
        for name, exc in plan:
            st = _study(exc)
            assert st.degree == 1 and st.conductor is None
            t0 = time.perf_counter()
            res = st.solve(comm)
            out[name] = {"study": st, "res": res, "exc": exc, "t": time.perf_counter() - t0}
            _p0(f"[WF-16 s6b] {name}: cells {res.cells}, S {res.s_matrix.shape}, "
                f"{out[name]['t']:.1f} s")
    n_studies = 1 + len(plan)
    _p0(f"[WF-16 s6b] build cache: {cache['real']} real build, {cache['hits']} hits "
        f"({n_studies} Studies)")
    assert cache["real"] == 1 and cache["hits"] == n_studies - 1 == 6
    out["omega"] = s0.omega
    return out


@complex_only
def test_wf16_s6b_r0_fixture_tie(s6b):
    """(R0) and (a'): the PEC hole fixture and the in-run ``C_mode``."""
    res0 = s6b["S0"]["res"]
    assert res0.cells == ANS6_EXACT_STRIP_HOLE_CELL_RECORD
    means = _class_means(np.asarray(res0.s_matrix, dtype=np.complex128))
    want = ANS6_EXACT_STRIP_CLASS_MEANS_P1["64"]["pec"]
    rel = {c: abs(means[c] - want[c]) / abs(want[c]) for c in CLASSES}
    c1 = s6b["t1"].c_mode_farad
    a = c1 / TH17_STEP3_C_MODE_RECORD_F - 1.0
    _p0("\n[WF-16 s6b R0] class means vs 'pec' record: " + ", ".join(
        f"{c} {rel[c]:.3e}" for c in CLASSES)
        + f" (band {LEG_D0_REPRODUCTION_BAND:.0e}; item 429 5.052e-11 / 3.878e-11 / 5.875e-11)")
    _p0(f"[WF-16 s6b a'] C_mode {c1:.9e} F vs TH-17 record {a:+.6e} (band {REDUCTION_BAND:.0e}; "
        f"item 429 -2.976552e-4); vs item 429's C_mode {ITEM429_C_MODE_F:.9e} F: "
        f"{c1 / ITEM429_C_MODE_F - 1:+.3e} (printed)")
    assert max(rel.values()) <= LEG_D0_REPRODUCTION_BAND
    assert abs(a) <= REDUCTION_BAND


@complex_only
def test_wf16_s6b_route_a_identity(s6b):
    """(I-A): Route A at ``C_mode`` equals the three-sheet circuit reduction of S0."""
    s0 = s6b["S0"]["res"].s_matrix
    z0 = float(REFERENCE_IMPEDANCE_OHM)
    sa = s6b["S_A"]
    zc = sa["exc"].impedance_ohm(FREQUENCY_64_HZ)
    assert [p.port_id for p in sa["res"].port_definitions] == ["P1"]
    assert sa["res"].capacitor_impedance_ohm == zc
    assert sa["res"].s_matrix.shape == (1, 1)
    pa = reduce_terminated_ports(s0, z0, {1: zc, 2: zc, 3: zc})
    r = abs(sa["res"].s_matrix[0, 0] - pa[0, 0]) / abs(pa[0, 0])
    _p0(f"\n[WF-16 s6b I-A] S_A S11 {complex(sa['res'].s_matrix[0, 0]):.9e} vs reduction "
        f"{complex(pa[0, 0]):.9e}: rel {r:.6e} ({r / REDUCTION_BAND:.4f} x REDUCTION_BAND "
        f"{REDUCTION_BAND:.0e}); copper prior {PRIOR_ROUTE_A_COPPER:.3e}, "
        f"S11 {PRIOR_ROUTE_A_S11_COPPER} (printed)")
    assert r <= REDUCTION_BAND


@complex_only
def test_wf16_s6b_per_sheet_identity(s6b):
    """(I-4): per sheet at ``C_mode``, r = 0, the Study 3x3 equals the reduction of S0."""
    s0 = s6b["S0"]["res"].s_matrix
    z0 = float(REFERENCE_IMPEDANCE_OHM)
    worst = 0.0
    for i in range(4):
        d = s6b[f"S{i + 1}"]
        zc = d["exc"].impedance_ohm(FREQUENCY_64_HZ)
        assert [p.port_id for p in d["res"].port_definitions] == [
            p for j, p in enumerate(PORTS) if j != i]
        assert d["res"].capacitor_impedance_ohm == zc
        pred = reduce_terminated_ports(s0, z0, {i: zc})
        r = _relf(d["res"].s_matrix, pred)
        worst = max(worst, r)
        _p0(f"[WF-16 s6b I-4] sheet P{i + 1}: {r:.6e} ({r / REDUCTION_BAND:.4f} x "
            f"REDUCTION_BAND {REDUCTION_BAND:.0e}); copper r=0 prior {PRIOR_I4_COPPER_R0[i]:.6e}")
    assert worst <= REDUCTION_BAND


@complex_only
def test_wf16_s6b_tuned_modal_by_construction(s6b):
    """(T): the tuned modal z_1 -- a by-construction identity, printed as a record.

    Asserted: the sheet excitation's own impedance at ``C_mode`` turns the API's
    z_1 into `TH-17`'s ``_loop_z`` (the excitation's ``impedance_ohm`` is the
    loop's ``1/(j omega C)``).  ``Im`` of the tuned z_1 is printed, never gated
    (rule (h); 2026-10-15 03:00 ruling).
    """
    s0 = np.asarray(s6b["S0"]["res"].s_matrix, dtype=np.complex128)
    z0 = float(REFERENCE_IMPEDANCE_OHM)
    t1 = s6b["t1"]
    zc = s6b["S_A"]["exc"].impedance_ohm(FREQUENCY_64_HZ)
    zt = _modal_z(s_to_z(s0, z0) + zc * np.eye(4))[0]
    z1 = complex(t1.z_k[1])
    loop = _loop_z(z1, FREQUENCY_64_HZ, t1.c_mode_farad)
    d_loop = abs(zt[1] - loop) / abs(z1)
    im_rel = abs(zt[1].imag) / abs(z1)
    _p0(f"\n[WF-16 s6b T] tuned z_1 = {zt[1].real:+.9e}{zt[1].imag:+.9e}j; |Im|/|z_1| = "
        f"{im_rel:.3e} (by construction, record; {im_rel / REDUCTION_BAND:.2e} x REDUCTION_BAND); "
        f"Re z_1 {z1.real:.8e} (item 429 {ITEM429_RE_Z1_OHM}); |z_1 - z_3|/|z_1| = "
        f"{abs(t1.z_k[1] - t1.z_k[3]) / abs(z1):.3e} (reciprocity, printed)")
    _p0(f"[WF-16 s6b T] |zt_1 - _loop_z|/|z_1| = {d_loop:.3e} (band EXACT_IDENTITY_RTOL "
        f"{EXACT_IDENTITY_RTOL:.0e})")
    assert d_loop <= EXACT_IDENTITY_RTOL


@complex_only
def test_wf16_s6b_k2_controls(s6b):
    """(c1), (c2): the k = 2 capacitance misses; factors printed, exceedance asserted."""
    s0 = np.asarray(s6b["S0"]["res"].s_matrix, dtype=np.complex128)
    z0 = float(REFERENCE_IMPEDANCE_OHM)
    t1, t2 = s6b["t1"], s6b["t2"]
    zc = s6b["S_A"]["exc"].impedance_ohm(FREQUENCY_64_HZ)
    zc2 = s6b["S_A2"]["exc"].impedance_ohm(FREQUENCY_64_HZ)
    omega = s6b["omega"]
    _p0(f"\n[WF-16 s6b c] zc2 {complex(zc2):.9e} vs 1/(j w C_2) "
        f"{1.0 / (1j * omega * t2.c_mode_farad):.9e} (printed)")
    z1 = complex(t1.z_k[1])
    zt2 = _modal_z(s_to_z(s0, z0) + zc2 * np.eye(4))[0]
    c1 = abs(zt2[1].imag) / abs(z1)
    sa = s6b["S_A"]["res"].s_matrix[0, 0]
    pa2 = reduce_terminated_ports(s0, z0, {1: zc2, 2: zc2, 3: zc2})
    c2 = abs(sa - pa2[0, 0]) / abs(pa2[0, 0])
    sa2 = s6b["S_A2"]["res"].s_matrix[0, 0]
    own2 = abs(sa2 - pa2[0, 0]) / abs(pa2[0, 0])
    _p0(f"[WF-16 s6b c1] |Im z_1(C_2)|/|z_1| = {c1:.4e} ({c1 / REDUCTION_BAND:.1f} x band; "
        f"predicted 0.3288)")
    _p0(f"[WF-16 s6b c2] S_A(C_mode) vs reduction at C_2: {c2:.4e} "
        f"({c2 / REDUCTION_BAND:.1f} x band; predicted 0.352); reduction at C_2 S11 "
        f"{complex(pa2[0, 0]):.6e} (predicted -0.5950+0.8033j)")
    _p0(f"[WF-16 s6b c] S_A2 vs its own reduction at C_2: {own2:.3e} "
        f"(the identity at the wrong C, printed; predicted approx 1e-5)")
    zin = {}
    for name, z in (("S_A", zc), ("S_A2", zc2)):
        s11 = complex(s6b[name]["res"].s_matrix[0, 0])
        zin[name] = z0 * (1 + s11) / (1 - s11) + z
    _p0(f"[WF-16 s6b c] field-side Z_in(P1)+zc at C_mode {zin['S_A']:.6e} ohm, "
        f"at C_2 {zin['S_A2']:.6e} ohm (printed; neither separates the modes)")
    assert c1 > REDUCTION_BAND
    assert c2 > REDUCTION_BAND
