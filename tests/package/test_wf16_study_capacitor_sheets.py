"""`WF-16` step 10 (§9 item 373, chain U37): capacitor sheets with ESR in the ``Study``.

``Study(FourPortBirdcageGeometry(as_hole=True, ...), gel, CapacitorSheets(C, esr_ohm,
on=...), conductor=LeontovichWall(COPPER_SIGMA), degree=1)`` on the F-small exact
four-port at 64 MHz, ``-n 2``, reproduces `PORT-24` step 2's identities:

- (R0) the 50 ohm Study (``LumpedPortDrive``) reproduces the ``ANS-6`` Cu class
  means at ``LEG_D0_REPRODUCTION_BAND`` through the edited function;
- (I) per sheet, the Study's reduced 3x3 at r = 0.2 ohm equals
  ``reduce_terminated_ports(S0, z0, {i: Z_RC})`` at ``REDUCTION_BAND``;
- (Z) at r = 0 the Study's 3x3 equals the direct capacitor-only 3x3
  (``_terminated_three_port`` with ``_ANS7.z_cap``) at ``LEG_D0_REPRODUCTION_BAND``;
- (P) the exact power identity with the sheets' dissipation counted at
  ``DISCRETE_IDENTITY_RTOL``.

Negative control (rule (e)): (N1) planted ESR 0.4 ohm on P2 misses the 0.2 ohm
reduction by more than ``REDUCTION_BAND`` (asserted, backed by `PORT-24`'s plant
arm, 10.95x); (N2) (Z)'s comparator against the r = 0.2 Study is printed only.
Every band is imported; none is set here.  Mesh economy (disclosed): inside the
fixture, ``FourPortBirdcageGeometry.build`` is wrapped so the first call builds
the real mesh and later equal-geometry calls return it; the hit count is asserted.
"""

from __future__ import annotations

import time

import numpy as np
import pytest
import ufl
from dolfinx import fem
from mpi4py import MPI

import fem_em_solver.study as study_mod
from fem_em_solver import (
    CapacitorSheets,
    FHumanBirdcageGeometry,
    FourPortBirdcageGeometry,
    LeontovichWall,
    LumpedPortDrive,
    Study,
)
from fem_em_solver.materials import GelledSalinePhantomMaterial
from fem_em_solver.ports.circuit import reduce_terminated_ports
from fem_em_solver.ports.lumped import run_lumped_sheet_port_case

from tests.complex_mode import complex_only
from tests.package.test_wf16_study_copper import _geometry
from tests.validation import test_th14_birdcage_copper as th14
from tests.validation.test_ans6_exact_strip import (
    ANS6_EXACT_STRIP_CLASS_MEANS_P1,
    ANS6_EXACT_STRIP_HOLE_CELL_RECORD,
    ANS6_EXACT_STRIP_RECORD_RANK_WIDTH,
    CLASSES,
    _class_means,
)
from tests.validation.test_ans7_exact_strip import C_F, _ANS7
from tests.validation.test_birdcage_power_identity import (
    DISCRETE_IDENTITY_RTOL,
    _sheet_field_dissipation_w,
    _source_power_w,
)
from tests.validation.test_lossy_sphere_fullwave import (
    FREQUENCY_64_HZ,
    SALINE_EPSILON_R,
    SALINE_SIGMA,
)
from tests.validation.test_port24_series_rc_sheet import (
    DEGREE,
    ESR_OHM,
    FREQUENCY_HZ,
    PLANT_INDEX,
    PLANT_OHM,
    POWER_DRIVEN,
    _rc,
    _rel,
    _terminated_three_port,
)
from tests.validation.test_port_birdcage_four_port import LEG_D0_REPRODUCTION_BAND
from tests.validation.test_port_lumped_rlc_termination import REDUCTION_BAND
from tests.validation.test_port_package_sparameters import REFERENCE_IMPEDANCE_OHM

PORTS = ("P1", "P2", "P3", "P4")
# Priors (PORT-24 step 2, same comparison, direct route), printed beside.
PRIOR_I_R02 = (3.049022e-05, 3.036163e-05, 3.040555e-05, 3.003774e-05)
PRIOR_P_R02 = 1.011e-12
PRIOR_PLANT = 1.095344e-02


def _p0(*args):
    if MPI.COMM_WORLD.rank == 0:
        print(*args, flush=True)


def _phantom():
    return GelledSalinePhantomMaterial(SALINE_SIGMA, SALINE_EPSILON_R, FREQUENCY_64_HZ)


def _study(excitation):
    return Study(_geometry(), _phantom(), excitation, frequency_hz=FREQUENCY_64_HZ,
                 degree=1, conductor=LeontovichWall(study_mod.COPPER_SIGMA))


def test_wf16_s10_api_and_refusals():
    """(K): no solve."""
    assert FREQUENCY_HZ == FREQUENCY_64_HZ and DEGREE == 1
    cs = CapacitorSheets(C_F, ESR_OHM, on=("P2",))
    assert cs.impedance_ohm(FREQUENCY_HZ) == _rc(ESR_OHM)
    assert CapacitorSheets(C_F, on=("P1",)).drive == LumpedPortDrive()
    assert CapacitorSheets(C_F, on=["P1", "P3"]).on == ("P1", "P3")
    for kwargs in (dict(c_f=0.0, on=("P1",)), dict(c_f=-1e-12, on=("P1",)),
                   dict(c_f=float("nan"), on=("P1",)), dict(c_f=float("inf"), on=("P1",)),
                   dict(c_f=C_F, esr_ohm=-0.1, on=("P1",)),
                   dict(c_f=C_F, esr_ohm=float("nan"), on=("P1",)),
                   dict(c_f=C_F, on=()), dict(c_f=C_F, on=("P1", "P1")),
                   dict(c_f=C_F, on=("P1",), drive="LumpedPortDrive")):
        with pytest.raises(ValueError):
            CapacitorSheets(**kwargs)
    # Refused before any mesh.
    with pytest.raises(TypeError):
        Study(FHumanBirdcageGeometry(), _phantom(), CapacitorSheets(C_F, on=("P1",)),
              frequency_hz=FREQUENCY_64_HZ, degree=1).solve()
    for on in (("P5",), PORTS):
        with pytest.raises(ValueError):
            _study(CapacitorSheets(C_F, on=on)).solve()
    st = _study(CapacitorSheets(C_F, ESR_OHM, on=("P2",)))
    with pytest.raises(TypeError):
        st.solve_sweep((FREQUENCY_64_HZ,))
    with pytest.raises(TypeError):
        st.sar(centres=[(0.0, 0.0, 0.0)])
    with pytest.raises(TypeError):
        st.b1_plus()
    with pytest.raises(TypeError):
        st.mean_sar()
    assert "CapacitorSheets" in study_mod.__all__


@pytest.fixture(scope="module")
def s10():
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

    plan = [("S0", LumpedPortDrive())]
    plan += [(f"S{i + 1}(0.2)", CapacitorSheets(C_F, ESR_OHM, on=(PORTS[i],))) for i in range(4)]
    plan += [("S_plant", CapacitorSheets(C_F, PLANT_OHM, on=(PORTS[PLANT_INDEX],)))]
    plan += [(f"S{i + 1}(0)", CapacitorSheets(C_F, 0.0, on=(PORTS[i],))) for i in range(4)]
    plan += [("S_A", CapacitorSheets(C_F, ESR_OHM, on=("P2", "P3", "P4")))]
    out = {}
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(FourPortBirdcageGeometry, "build", _build)
        for name, exc in plan:
            st = _study(exc)
            t0 = time.perf_counter()
            res = st.solve(comm)
            out[name] = {"res": res, "study": st, "t": time.perf_counter() - t0}
            _p0(f"[WF-16 s10] {name}: {type(exc).__name__}, cells {res.cells}, "
                f"S {res.s_matrix.shape}, {out[name]['t']:.1f} s")
    _p0(f"[WF-16 s10] build cache: {cache['real']} real build, {cache['hits']} hits "
        f"({len(plan)} Studies)")
    assert cache["real"] == 1 and cache["hits"] == len(plan) - 1
    out["omega"] = out["S0"]["study"].omega
    return out


@complex_only
def test_wf16_s10_cu_record(s10):
    """(R0): the 50 ohm Study reproduces the Cu record through the edited function."""
    res = s10["S0"]["res"]
    assert res.cells == ANS6_EXACT_STRIP_HOLE_CELL_RECORD
    means = _class_means(res.s_matrix)
    want = ANS6_EXACT_STRIP_CLASS_MEANS_P1["64"]["cu"]
    rel = {c: abs(means[c] - want[c]) / abs(want[c]) for c in CLASSES}
    _p0("[WF-16 s10] (R0) Cu class means vs record: " + ", ".join(
        f"{c} {rel[c]:.3e}" for c in CLASSES)
        + f" (band {LEG_D0_REPRODUCTION_BAND:.0e}; prior 4.093e-11 / 8.673e-11 / 1.376e-10)")
    assert max(rel.values()) <= LEG_D0_REPRODUCTION_BAND


@complex_only
def test_wf16_s10_per_sheet_identity(s10):
    """(I): Study 3x3 at r = 0.2 vs the circuit reduction of S0, per sheet."""
    s0 = s10["S0"]["res"].s_matrix
    z0 = float(REFERENCE_IMPEDANCE_OHM)
    worst = 0.0
    for i in range(4):
        res = s10[f"S{i + 1}(0.2)"]["res"]
        assert [p.port_id for p in res.port_definitions] == [p for j, p in enumerate(PORTS) if j != i]
        assert res.capacitor_impedance_ohm == _rc(ESR_OHM)
        pred = reduce_terminated_ports(s0, z0, {i: _rc(ESR_OHM)})
        r = _rel(res.s_matrix, pred)
        worst = max(worst, r)
        _p0(f"[WF-16 s10] (I) sheet P{i + 1} r = {ESR_OHM:g}: {r:.6e} "
            f"({r / REDUCTION_BAND:.4f} x REDUCTION_BAND {REDUCTION_BAND:.0e}); "
            f"prior {PRIOR_I_R02[i]:.6e}")
    # Printed, never asserted: Route A's S11 against the three-sheet reduction.
    sa = s10["S_A"]["res"].s_matrix
    pa = reduce_terminated_ports(s0, z0, {1: _rc(ESR_OHM), 2: _rc(ESR_OHM), 3: _rc(ESR_OHM)})
    _p0(f"[WF-16 s10] S_A S11 {complex(sa[0, 0]):.9e} vs reduction {complex(pa[0, 0]):.9e}: "
        f"rel {abs(sa[0, 0] - pa[0, 0]) / abs(pa[0, 0]):.3e} (printed only)")
    assert worst <= REDUCTION_BAND


@complex_only
def test_wf16_s10_r0_equals_capacitor_only(s10):
    """(Z): the Study at r = 0 equals the direct capacitor-only 3x3 (fresh factor)."""
    omega = s10["omega"]
    zc = _ANS7.z_cap(omega, C_F)
    worst = 0.0
    for i in range(4):
        res = s10[f"S{i + 1}(0)"]["res"]
        t0 = time.perf_counter()
        cap = _terminated_three_port(res.problem, res.all_port_definitions, res.specs,
                                     res.facet_tags, res.conductor_terms, i, zc)
        r = _rel(res.s_matrix, cap)
        worst = max(worst, r)
        r02 = s10[f"S{i + 1}(0.2)"]["res"].s_matrix
        n2 = _rel(r02, cap)
        _p0(f"[WF-16 s10] (Z) sheet P{i + 1}: Study r=0 vs direct z_cap {r:.6e} "
            f"(band {LEG_D0_REPRODUCTION_BAND:.0e}); (N2) Study r=0.2 vs it {n2:.6e} "
            f"(predicted ~1e-2, printed only); {time.perf_counter() - t0:.1f} s")
    assert worst <= LEG_D0_REPRODUCTION_BAND


@complex_only
def test_wf16_s10_power_identity(s10):
    """(P): P1 driven, P2-P4 series-RC; sheets' dissipation counted.

    Accounted as `PORT-24` step 2's (P) (tests/validation/test_port24_series_rc_sheet.py,
    the power-identity block), with ``P_surf`` on the Study's conductor facets as
    step 7's ``_identity`` does.
    """
    comm = MPI.COMM_WORLD
    res = s10["S_A"]["res"]
    omega = s10["omega"]
    msh = res.mesh
    _r, fields = run_lumped_sheet_port_case(
        res.problem, res.all_port_definitions, res.specs, facet_tags=res.facet_tags,
        driven_port_id=POWER_DRIVEN, verbose=False, return_fields=True,
        extra_bilinear_terms=res.conductor_terms, degree=res.degree)
    e = fields.e_complex
    dx = ufl.Measure("dx", domain=msh, subdomain_data=res.cell_tags)
    p_vol = float(np.real(comm.allreduce(fem.assemble_scalar(fem.form(
        0.5 * fields.sigma_field * ufl.inner(e, e) * dx)), op=MPI.SUM)))
    acct = {"mesh": msh, "facet_tags": res.facet_tags}
    sheets = {s.port_id: s.sheet(driven=(s.port_id == POWER_DRIVEN)) for s in res.specs}
    p_sheet = {pid: _sheet_field_dissipation_w(acct, sh, e, omega) for pid, sh in sheets.items()}
    p_src = float(_source_power_w(acct, sheets[POWER_DRIVEN], e, omega))
    z_s = LeontovichWall(study_mod.COPPER_SIGMA).surface_impedance_ohm(omega)
    p_surf = th14._surface_loss_w(msh, res.conductor_facet_tags, z_s, e, comm)
    p_all = sum(p_sheet.values())
    p_rc = sum(v for pid, v in p_sheet.items() if pid != POWER_DRIVEN)
    resid = abs(p_src - p_all - p_vol - p_surf) / abs(p_src)
    resid_dz = abs(p_src - (p_all - p_rc) - p_vol - p_surf) / abs(p_src)
    p_rest = p_vol + p_surf
    q_c = 1.0 / (omega * C_F * ESR_OHM)
    _p0(f"[WF-16 s10] (P) P_src {p_src:.9e} W, P_sheets {p_all:.9e} W (RC "
        + ", ".join(f"{k} {v:.9e}" for k, v in p_sheet.items() if k != POWER_DRIVEN)
        + f"), P_vol {p_vol:.9e} W, P_surf {p_surf:.9e} W; residual {resid:.3e} "
        f"(<= {DISCRETE_IDENTITY_RTOL:g}; prior {PRIOR_P_R02:.3e})")
    _p0(f"[WF-16 s10] (Dz) {resid_dz:.6e} (prior 7.031342e-03); Q(r)/Q(0) "
        f"{p_rest / (p_rest + p_rc):.6e} (prior 2.922955e-01); Q_rest "
        f"{p_rc / p_rest * q_c:.6e} (prior 4.541422e+02); printed only")
    assert p_surf > 0.0 and p_rc > 0.0
    assert resid <= DISCRETE_IDENTITY_RTOL


@complex_only
def test_wf16_s10_planted_esr_is_seen(s10):
    """(N1): a planted 0.4 ohm on P2 misses the 0.2 ohm reduction (asserted)."""
    s0 = s10["S0"]["res"].s_matrix
    res = s10["S_plant"]["res"]
    pred = reduce_terminated_ports(s0, float(REFERENCE_IMPEDANCE_OHM),
                                   {PLANT_INDEX: _rc(ESR_OHM)})
    r = _rel(res.s_matrix, pred)
    _p0(f"[WF-16 s10] (N1) planted r = {PLANT_OHM:g} on P{PLANT_INDEX + 1}: {r:.6e} = "
        f"{r / REDUCTION_BAND:.3f} x REDUCTION_BAND (printed; prior {PRIOR_PLANT:.6e} = 10.953x)")
    assert r > REDUCTION_BAND
