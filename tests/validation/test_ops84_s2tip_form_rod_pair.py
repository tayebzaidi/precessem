"""OPS-84 (§9 item 860): the (S2-tip) comparand's form on a field-derived source.

I41's matched rod pair (`test_imp3_rod_heating_matched`): arm M (the kept-cells
build, tag 403 mapped to gel and merged into the phantom) and arm R (the PEC rod
hole), run in I41's order without its control arm.  Inside each arm's field hook,
for s ∈ {ccw, P1}: q_s = G·(DG0 cell average of ½σ|E_s|² on the phantom) on the
phantom submesh, one backward-Euler step from T∞ at Δt₁ = ``DT_1`` and Δt₁/10,
and the 1 g ball means at both ``TIPS`` of θ¹ − T∞ against ⟨q⟩ (the gate as
written), ⟨Πq⟩ (consistent-mass L2 projection onto the P1 temperature space)
and ⟨Π_L q⟩ (lumped), via `OPS-83`'s ``_projections``.

Asserted (bands imported, none set): (A0m) as I41; (B) ≤ ``B_BAND`` on the four
(arm, source) pairs; (E) ≤ ``E_BAND`` and the global (S2) ≤ ``S2_BAND`` on all
eight solves; the reproduction control — the ccw Δt₁ global ratio equals
`OPS-81` W2's printed 0.999994872614 (arm M) / 0.999995000025 (arm R) within
1e-11 (`20261105T183330Z_OPS-81-W2-i41.log:3683`, `:3728`); the projection
identity |⟨Πθ⟩_B/⟨θ⟩_B − 1| ≤ ``PROJECTION_ID_TOL`` at both tips on both arms.

Printed only (rule (e)): R_DG0, R_Π, R_L, their Δt-ratios, ⟨Πq⟩/⟨q⟩ − 1 and the
pre-stated reading per arm × source × tip.  Not claimed: a re-pose of the tip
gate, a band, a reading on U127's mesh or on F-clinical.
"""

from __future__ import annotations

import gc
import os
import time

import numpy as np
import pytest
import ufl
from dolfinx import default_scalar_type, fem
from mpi4py import MPI

import tests.validation.test_imp3_rod_birdcage as _rb
from fem_em_solver.post.sar import averaging_ball_radius, mean_sar

from tests.complex_mode import complex_only
from tests.mesh.test_birdcage_implant_rod import (
    KEPT_ROD_CELL_RECORD,
    _build_kept_rod,
    _build_rod,
)
from tests.validation.test_birdcage_sar_mass_averaged import PHANTOM_RHO_KG_PER_M3
from tests.validation.test_imp3_rod_birdcage import TIPS
from tests.validation.test_imp3_rod_heating_matched import (
    ROD_CELLS,
    _arm_with_403,
    _census,
    _merge_rod_into_phantom,
)
from tests.validation.test_imp3_rod_heating_phantom import (
    _outer_robin_tags,
    _parent_wall_count,
)
from tests.validation.test_mass_averaged_sar_standard_masses import (
    ONE_GRAM_KG,
    QUADRATURE_DEGREE,
)
from tests.validation.test_mat10_sar_to_heat_sphere import B_BAND, ROBIN_TAG, _allsum, _to_sub
from tests.validation.test_mat15_coil_driven_heating_phantom import (
    DT_1,
    E_BAND,
    RHO_C,
    S2_BAND,
    T_INF,
    _cell_average_dg0_tag,
    _max_imag,
    _tag_submesh,
    _transient,
)
from tests.validation.test_ops83_s2tip_form_on_stub import (
    CONST_RATIO,
    PROJECTION_ID_TOL,
    SCALING_RATIO,
    _projections,
)
from tests.validation.test_port_birdcage_lumped_column import PHANTOM_CELL_TAG

KNOB = "FEM_EM_OPS84"
SOURCES = ("ccw", "P1")
DTS = {"dt1": DT_1, "dt1/10": DT_1 / 10.0}
# `OPS-81` W2 printed ccw Δt₁ global ratios (`20261105T183330Z_OPS-81-W2-i41.log:3683`, `:3728`).
I41_S2_RATIO = {"M": 0.999994872614, "R": 0.999995000025}
REPRO_TOL = 1.0e-11  # ten units of the `.12f` print's last digit (item 860 brief)


def _say(msg):
    if MPI.COMM_WORLD.rank == 0:
        print(msg, flush=True)


def _ball_reader(sub, funcs):
    """1 g indicator-ball means of each named integrand, I39's form (`:202–216`)."""
    comm = sub.comm
    centre_c = fem.Constant(sub, np.zeros(3, dtype=default_scalar_type))
    rad_c = fem.Constant(sub, default_scalar_type(0.0))
    xs = ufl.SpatialCoordinate(sub)
    off = xs - centre_c
    inside = ufl.conditional(ufl.lt(ufl.real(ufl.dot(off, off)), ufl.real(rad_c) ** 2), 1.0, 0.0)
    dxq = ufl.dx(domain=sub, metadata={"quadrature_degree": int(QUADRATURE_DEGREE)})
    f_vol = fem.form(inside * dxq)
    forms = {k: fem.form(f * inside * dxq) for k, f in funcs.items()}

    def _read(centre, radius):
        centre_c.value[:] = np.asarray(centre, dtype=float)
        rad_c.value = radius
        vol = float(np.real(comm.allreduce(fem.assemble_scalar(f_vol), op=MPI.SUM)))
        out = {k: float(np.real(comm.allreduce(fem.assemble_scalar(f), op=MPI.SUM))) / vol
               for k, f in forms.items()}
        return vol, out

    return _read


def _reading(r):
    """The pre-stated reading for one arm × source × tip (printed, never asserted)."""
    dg0_1, pi_1 = r["DG0"]["dt1"], r["Pi"]["dt1"]
    const = CONST_RATIO[0] <= r["ratio"]["DG0"] <= CONST_RATIO[1]
    if dg0_1 <= S2_BAND:
        return "INCONCLUSIVE"
    if const and pi_1 <= S2_BAND:
        return "form ACCOUNTS"
    if const and pi_1 > S2_BAND and r["ratio"]["Pi"] >= SCALING_RATIO:
        return "form + O(dt)"
    if pi_1 > S2_BAND and r["ratio"]["Pi"] < SCALING_RATIO:
        return "RESIDUAL"
    return "INCONCLUSIVE"


def _slope_arm(arm, label, msh, cell_tags, e_of, sigma_field, G_ccw):
    comm = MPI.COMM_WORLD
    t0 = time.perf_counter()
    sub, parent_ids = _tag_submesh(msh, cell_tags, PHANTOM_CELL_TAG)
    ft, counts = _outer_robin_tags(sub)
    one = fem.Constant(sub, default_scalar_type(1.0))
    v_sub = float(np.real(_allsum(comm, one * ufl.dx(domain=sub))))
    a1 = averaging_ball_radius(mass_kg=ONE_GRAM_KG, rho=PHANTOM_RHO_KG_PER_M3)
    v_ball = 4.0 / 3.0 * np.pi * a1 ** 3
    _say(f"\n[OPS-84 {label}] sub cells {sub.topology.index_map(3).size_global}; V_sub "
         f"{v_sub:.12e} m^3; Robin facets {counts} (tag {ROBIN_TAG}); a_1g {a1 * 1e3:.4f} mm, "
         f"4/3 pi a^3 {v_ball:.6e} m^3; quadrature degree {QUADRATURE_DEGREE}")
    out = {"B": {}, "E": {}, "S2": {}, "ratio_glob": {}, "R": {}, "imag": 0.0, "pid": {},
           "G": {}, "P": {}}
    for s in SOURCES:
        P = float(np.real(mean_sar(e_of[s], sigma=sigma_field, rho=PHANTOM_RHO_KG_PER_M3,
                                   cell_tags=cell_tags, comm=comm,
                                   subdomain_ids=PHANTOM_CELL_TAG)["dissipated_power_w"]))
        if s == "ccw":
            G = (1.0 / P) if G_ccw is None else G_ccw
        else:
            G = 1.0 / P
        dep = G * P
        out["G"][s], out["P"][s] = G, P
        Wp, q_par = _cell_average_dg0_tag(msh, cell_tags, PHANTOM_CELL_TAG,
                                          G * 0.5 * sigma_field * ufl.inner(e_of[s], e_of[s]))
        q = _to_sub(Wp, q_par, sub, parent_ids)
        int_q = float(np.real(_allsum(comm, q * ufl.dx(domain=sub))))
        out["B"][s] = abs(int_q - dep) / abs(dep)
        _say(f"[OPS-84 {label}] source {s}: P_phantom {P:.12e} W, G {G:.12e} 1/W, G*P "
             f"{dep:.12e} W; (B) int_sub q {int_q:.15e}: rel {out['B'][s]:.3e} (band "
             f"{B_BAND:.0e}, ASSERTED)")
        th = {}
        V1 = None
        for dname, dt in DTS.items():
            res = _transient(sub, ft, q, dt=dt, n_steps=1, theta=1.0)
            thbar = float(np.real(_allsum(comm, (res.temperature - T_INF)
                                          * ufl.dx(domain=sub)))) / v_sub
            ratio = RHO_C * thbar / (dt * dep / v_sub)
            out["ratio_glob"][(s, dname)] = ratio
            out["S2"][(s, dname)] = abs(ratio - 1.0)
            out["E"][(s, dname)] = float(res.max_balance_residual)
            out["imag"] = max(out["imag"], _max_imag(res.temperature))
            V1 = res.temperature.function_space
            f = fem.Function(V1)
            f.x.array[:] = res.temperature.x.array - T_INF
            f.x.scatter_forward()
            th[dname] = f
            _say(f"[OPS-84 {label}] source {s} {dname} ({dt:.6e} s): global rho c thbar/(dt G P/V) "
                 f"= {ratio:.12f}, (S2) {out['S2'][(s, dname)]:.4e} (band {S2_BAND:.0e}, "
                 f"ASSERTED); (E) {out['E'][(s, dname)]:.3e} (band {E_BAND:.0e}, ASSERTED)")
        pq, plq = _projections(q, V1, f"ops84_{arm}_{s}")
        read = _ball_reader(sub, {"th_dt1": th["dt1"], "th_dt1/10": th["dt1/10"], "q": q,
                                  "pq": pq, "plq": plq})
        if s == "ccw":
            pth, _ = _projections(th["dt1"], V1, f"ops84_{arm}_id")
            read_id = _ball_reader(sub, {"th": th["dt1"], "pth": pth})
        for tip, c in TIPS.items():
            vol, mb = read(c, a1)
            r = {"DG0": {}, "Pi": {}, "L": {}, "ratio": {}, "vol": vol,
                 "pq/q": mb["pq"] / mb["q"] - 1.0}
            for dname, dt in DTS.items():
                lhs = RHO_C * mb[f"th_{dname}"] / dt
                r["DG0"][dname] = abs(lhs / mb["q"] - 1.0)
                r["Pi"][dname] = abs(lhs / mb["pq"] - 1.0)
                r["L"][dname] = abs(lhs / mb["plq"] - 1.0)
            for k in ("DG0", "Pi", "L"):
                b = r[k]["dt1/10"]
                r["ratio"][k] = r[k]["dt1"] / b if b > 0 else np.inf
            r["reading"] = _reading(r)
            out["R"][(s, tip)] = r
            _say(f"[OPS-84 {label}] {s} {tip} tip {tuple(c)}: ball vol {vol:.6e} m^3 "
                 f"({vol / v_ball:.6f} x 4/3 pi a^3); <q> {mb['q']:.12e} <Pq> {mb['pq']:.12e} "
                 f"<PLq> {mb['plq']:.12e} W/m^3; <Pq>/<q> - 1 = {r['pq/q']:+.4e}")
            if s == "ccw":
                _, mi = read_id(c, a1)
                out["pid"][tip] = abs(mi["pth"] / mi["th"] - 1.0)
                _say(f"[OPS-84 {label}] projection identity {tip}: |<P theta>/<theta> - 1| = "
                     f"{out['pid'][tip]:.3e} (tol {PROJECTION_ID_TOL:.0e}, ASSERTED)")
    _say(f"[OPS-84 {label}] max|Im theta| {out['imag']:.3e} (printed); thermal elapsed "
         f"{time.perf_counter() - t0:.1f} s")
    return out


def _run():
    comm = MPI.COMM_WORLD
    t0 = time.perf_counter()
    res = {}
    built_rod = _build_rod()
    c_rod = _census(built_rod)
    rod_wall = _parent_wall_count(built_rod[0], built_rod[2])
    built_m = _build_kept_rod()
    c_m = _census(built_m)
    _say(f"\nOPS-84 (A0m) census: arm M {c_m['cells']} cells (record {KEPT_ROD_CELL_RECORD}), "
         f"tag-3 {c_m['t3']}, tag-403 {c_m['t403']} (expected {ROD_CELLS}); arm R "
         f"{c_rod['cells']} cells, tag-3 {c_rod['t3']}; R parent tag-402 owned {rod_wall}")
    assert c_m["cells"] == KEPT_ROD_CELL_RECORD, c_m
    assert c_m["t403"] == ROD_CELLS, c_m
    assert c_m["t3"] == c_rod["t3"], (c_m, c_rod)

    def hook_m(label, rod, msh, cell_tags, e_of, sigma_field, rho_field):
        merged = _merge_rod_into_phantom(msh, cell_tags)
        res["M"] = _slope_arm("M", label, msh, merged, e_of, sigma_field, None)

    def hook_r(label, rod, msh, cell_tags, e_of, sigma_field, rho_field):
        res["R"] = _slope_arm("R", label, msh, cell_tags, e_of, sigma_field, res["M"]["G"]["ccw"])

    m = _arm_with_403(built_m, "OPS-84 arm M (rod as gel)", hook_m, map_gel=True)
    del built_m
    gc.collect()
    r = _rb._arm(built_rod, "OPS-84 arm R (with rod)", rod=True, field_hook=hook_r)
    del built_rod
    gc.collect()

    for arm in ("M", "R"):
        rr = I41_S2_RATIO[arm]
        got = res[arm]["ratio_glob"][("ccw", "dt1")]
        _say(f"OPS-84 reproduction control arm {arm}: ccw dt1 global ratio {got:.15f} vs OPS-81 W2 "
             f"{rr:.12f}: |diff| {abs(got - rr):.3e} (tol {REPRO_TOL:.0e}, ASSERTED)")
    _say("OPS-84 table: arm source tip | R_DG0(dt1) R_DG0(dt1/10) ratio | R_Pi(dt1) "
         "R_Pi(dt1/10) ratio | R_L(dt1) R_L(dt1/10) ratio | <Pq>/<q>-1 | reading")
    for arm in ("M", "R"):
        for s in SOURCES:
            for tip in TIPS:
                x = res[arm]["R"][(s, tip)]
                row = " | ".join(f"{x[k]['dt1']:.4e} {x[k]['dt1/10']:.4e} {x['ratio'][k]:.3f}"
                                 for k in ("DG0", "Pi", "L"))
                _say(f"OPS-84 table: {arm} {s} {tip} | {row} | {x['pq/q']:+.4e} | {x['reading']}")
    for arm in ("M", "R"):
        for s in SOURCES:
            for tip in TIPS:
                _say(f"OPS-84 READING (printed, never asserted): arm {arm} source {s} tip {tip}: "
                     f"{res[arm]['R'][(s, tip)]['reading']}")
    _say(f"OPS-84 timing: arms M {m['t_arm']:.1f} / R {r['t_arm']:.1f} s; total "
         f"{time.perf_counter() - t0:.1f} s at -n {comm.size}")
    return {"res": res, "c_m": c_m, "c_rod": c_rod, "M": m}


@pytest.fixture(scope="module")
def ops84():
    if os.environ.get(KNOB) != "1":
        pytest.skip(f"{KNOB}=1 not set (heavy window, ≈ 4-6 min at -n 2)")
    return _run()


# --- coding tests (W0) -------------------------------------------------------

def test_knob_refuses_without_value():
    assert os.environ.get(KNOB) in (None, "1"), f"{KNOB} must be unset or 1"


def test_projection_tags_unique():
    tags = [f"ops84_{a}_{s}" for a in ("M", "R") for s in SOURCES]
    tags += [f"ops84_{a}_id" for a in ("M", "R")]
    assert len(set(tags)) == len(tags)


def test_reading_classes():
    def mk(dg0, pi, rdg0, rpi):
        return {"DG0": {"dt1": dg0}, "Pi": {"dt1": pi}, "ratio": {"DG0": rdg0, "Pi": rpi}}
    big, small = 10 * S2_BAND, 0.1 * S2_BAND
    assert _reading(mk(small, small, 1.0, 1.0)) == "INCONCLUSIVE"
    assert _reading(mk(big, small, 1.0, 10.0)) == "form ACCOUNTS"
    assert _reading(mk(big, big, 1.0, 10.0)) == "form + O(dt)"
    assert _reading(mk(big, big, 1.0, 1.0)) == "RESIDUAL"


# --- W1 ------------------------------------------------------------------------

@complex_only
def test_a0m_census(ops84):
    assert ops84["c_m"]["cells"] == KEPT_ROD_CELL_RECORD
    assert ops84["c_m"]["t403"] == ROD_CELLS
    assert ops84["c_m"]["t3"] == ops84["c_rod"]["t3"]
    assert ops84["M"]["material_403"]["mapped"]


@complex_only
@pytest.mark.parametrize("arm", ["M", "R"])
def test_identities(ops84, arm):
    h = ops84["res"][arm]
    for s in SOURCES:
        assert h["B"][s] <= B_BAND, f"(B) {arm} {s}: {h['B'][s]:.3e}"
        for d in DTS:
            assert h["E"][(s, d)] <= E_BAND, f"(E) {arm} {s} {d}: {h['E'][(s, d)]:.3e}"
            assert h["S2"][(s, d)] <= S2_BAND, f"(S2) {arm} {s} {d}: {h['S2'][(s, d)]:.4e}"


@complex_only
@pytest.mark.parametrize("arm", ["M", "R"])
def test_reproduction_control(ops84, arm):
    got = ops84["res"][arm]["ratio_glob"][("ccw", "dt1")]
    assert abs(got - I41_S2_RATIO[arm]) <= REPRO_TOL, (arm, got)


@complex_only
@pytest.mark.parametrize("arm", ["M", "R"])
def test_projection_identity(ops84, arm):
    for tip, v in ops84["res"][arm]["pid"].items():
        assert v <= PROJECTION_ID_TOL, (arm, tip, v)
