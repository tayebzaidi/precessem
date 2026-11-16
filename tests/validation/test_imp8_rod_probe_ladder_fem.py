"""`IMP-8` step 19 (§9 item 895, §10 I81): the gel / Ti / hole bracket at the
clinical probe placement and on a 1 / 2 / 3 / 5 mm axial ladder.

Brief: docs/planning/chunks/IMP-8.md, "§9 item 895".  I75's fixture body
(``test_imp8_rod_probe_f2182_fem.probe_run``, I55's graded gel box, real build,
degree 2, CN to t* = 900 s, I55's uniform rod source, I73's body cells, I66's
``_coefficients``) with the solve call byte-identical, and three read batches,
one per call: I73's ``READ_PTS``; I75's ``PROBE_PTS`` (final field + each
snapshot), unchanged so (R75) compares call for call; and the new
``LADDER_PTS`` (z = ±17, ±18, ±20 mm = L/2 + 2, 3, 5 mm) on the final field and
each snapshot.  I75's ``probe_run`` fixture is not imported (it would make a
second solve available to collection).

Arm knob: I75's ``FEM_EM_IMP8_I75_ARM`` = ``ti`` (default) / ``gel`` / ``hole``.

Asserted (bands imported): I75's gates on every arm — (E) ≤ ``E_BAND``, (S0) ≤
``B_BAND``, cells = ``W0C_CELLS``, body cells = ``BODY_CELLS``, no non-finite
coefficient dofs, snapshots at 225 / 450 / 675 s, θ_probe(900) > θ_probe(675) at
both ±16 mm probes; I75's arm gates unchanged (gel: I73's six θ and I55's (I)
at ``INTERIOR_BAND``; ti: ``TI_THETA``; hole: θ(2, 0) moves off I55's by more
than ``CONTROL_FLOOR``); (R75) I75's printed probe reads at ±16 mm, θ(900 s) and
the D = 0.5 superposed θ(900 s), within one unit of the sixth significant digit;
the control, at each sign: θ(16) > θ(17) > θ(18) > θ(20) strictly and
1 − θ(17)/θ(16) > ``INTERIOR_BAND`` (factor 1; ceiling 1).

Printed only (🟡 record): per ladder point θ(900 s), θ at 225 / 450 / 675 s,
θ(900)/θ(675), the D = 0.5 superposed θ(900 s) (a superposition, not a switched
solve), θ(z)/θ(16), the ± mean and |difference|, and beside each the
infinite-medium line-source closed forms QP·Φ and ``_duhamel`` with FEM /
closed form.
"""

from __future__ import annotations

import time

import numpy as np
import pytest

from fem_em_solver.thermal import solve_transient_conduction

from tests.mesh.test_birdcage_implant_rod import ROD_RADIUS
from tests.validation.test_imp8_composite_cylinder_closed_form import (
    ROD_K,  # nominal
    ROD_RHO_C,  # nominal
)
from tests.validation.test_imp8_rod_conducting_body_fem import EPS_HOLE, _coefficients
from tests.validation.test_imp8_rod_line_source_closed_form import L, MM, T_PRINT, _phi
from tests.validation.test_imp8_rod_line_source_duty_closed_form import PERIOD, _duhamel
from tests.validation.test_imp8_rod_line_source_duty_fem import STEP_T2, T2
from tests.validation.test_imp8_rod_line_source_fem import (  # noqa: F401  (box: fixture)
    DEGREE,
    GATE_PTS,
    QP,
    T_INF,
    W0C_CELLS,
    _need,
    _source_and_s0,
    _theta_pts,
    box,
)
from tests.validation.test_imp8_rod_probe_f2182_fem import (
    ARM_ENV,
    ARMS,
    DUTY,
    GEL_NEW,
    PROBE_PTS,
    SNAP_STEPS,
    TI_THETA,
    _arm,
    _body_cells,
    _sixth,
    _superpose,
)
from tests.validation.test_imp8_rod_self_heated_fem import (
    BODY_CELLS,
    CONTROL_FLOOR,
    I55_THETA,
    READ_PTS,
)
from tests.validation.test_mat10_sar_to_heat_sphere import B_BAND, ROBIN_TAG, _say
from tests.validation.test_mat12_transient_sar_heating_sphere import E_BAND, RHO_C
from tests.validation.test_mat13_transient_perfusion_sphere import INTERIOR_BAND
from tests.validation.test_mat15_coil_driven_heating_phantom import H_GEL, K_GEL, N_MAIN

#: mm beyond L/2 for the new ladder rungs (rung 1, 1 mm = 16 mm, is I75's probe batch).
LADDER_OFF = (2.0, 3.0, 5.0)
LADDER_PTS = ((0.0, 17.0), (0.0, -17.0), (0.0, 18.0), (0.0, -18.0), (0.0, 20.0), (0.0, -20.0))
Z_ALL = (16.0, 17.0, 18.0, 20.0)  # mm, the four rungs (16 = I75's probe; 17 = the clinical probe)
#: I75's printed probe reads, K/W, (θ(900 s) at +16, −16; D = 0.5 at +16, −16):
#: ti `20261108T194323Z_IMP-8-item881-W1.log:120–121`, gel
#: `20261108T195853Z_IMP-8-item881-W2b.log:120–121`, hole
#: `20261108T200636Z_IMP-8-item881-W3.log:120–121`.
I75_PROBE = {
    "ti": ((1.064117e01, 1.049960e01), (2.109553e00, 2.099602e00)),
    "gel": ((9.220827e00, 9.102961e00), (1.979610e00, 1.970586e00)),
    "hole": ((8.815610e00, 8.750725e00), (1.895292e00, 1.888582e00)),
}
#: I80's (item 894) line record, continuous θ(900 s), K/W, printed beside only
#: (`20261109T043829Z_IMP-8-item894-W1.log:98–107`); the inline QP·Φ is the comparand.
I80_LINE = {16.0: 9.648662, 17.0: 6.810428, 18.0: 5.242317, 20.0: 3.427288}


def test_k0_ladder_points():
    half_mm = L / MM / 2.0
    assert abs(half_mm - 15.0) < 1e-9, half_mm
    want = []
    for off in LADDER_OFF:
        z = half_mm + off
        want += [(0.0, z), (0.0, -z)]
    assert len(LADDER_PTS) == len(want)
    for (r, z), (rw, zw) in zip(LADDER_PTS, want):
        assert r == rw == 0.0
        assert abs(z - zw) < 1e-9, (z, zw)
    assert PROBE_PTS == ((0.0, half_mm + 1.0), (0.0, -(half_mm + 1.0)))
    assert set(I75_PROBE) == set(ARMS)


def _by_z(probe, ladder, j):
    """{signed z (mm): θ at step j} from the probe and ladder batches."""
    d = {z: probe[j][i] for i, (_, z) in enumerate(PROBE_PTS)}
    d.update({z: ladder[j][i] for i, (_, z) in enumerate(LADDER_PTS)})
    return d


@pytest.fixture(scope="module")
def ladder_run(box):
    _need("transient")
    arm = _arm()
    msh, ct, ft, n = box
    q, out = _source_and_s0(msh)
    body_idx, n_body = _body_cells(msh, q)
    k, rho_c, n_bad = _coefficients(msh, body_idx, arm)
    out.update({"arm": arm, "n": n, "n_body": n_body, "n_bad": n_bad})
    _say(f"[imp8-i81] arm {arm}; ROD_K {ROD_K:g}, ROD_RHO_C {ROD_RHO_C:.6g} (nominal); K_GEL {K_GEL:g}, "
         f"RHO_C {RHO_C:.6g}; EPS_HOLE {EPS_HOLE:.0e}; ROD_RADIUS {ROD_RADIUS:g} m, L {L:g} m")
    _say(f"[imp8-i81] body cells (positive q dofs -> cells, owned, allreduced) {n_body} (BODY_CELLS "
         f"{BODY_CELLS}); source-selected cells {out['n_sel']}")
    dt = T_PRINT / N_MAIN
    out["snap_t"] = [j * dt for j in SNAP_STEPS]
    t1 = time.perf_counter()
    res = solve_transient_conduction(
        msh, ft, k=k, rho_c=rho_c, q=q, h=H_GEL, t_inf=T_INF, t_init=T_INF,
        robin_tag=ROBIN_TAG, dt=dt, n_steps=N_MAIN, theta=0.5, startup_half_steps=2,
        degree=DEGREE, snapshot_steps=SNAP_STEPS)
    out["t_call"] = time.perf_counter() - t1
    out["E"] = float(res.max_balance_residual)
    out["snap_keys"] = sorted(int(s) for s in res.snapshots)
    _say(f"[imp8-i81] CN N={N_MAIN} dt={dt:g} s: t_call {out['t_call']:.1f} s; n_factorizations "
         f"{res.n_factorizations}; ranks {msh.comm.size}; snapshots {out['snap_keys']} at t "
         f"{out['snap_t']} s (I60 T2 {T2:g} s = step {STEP_T2})")

    # 1. I73's batch, unchanged.
    th = _theta_pts(res.temperature, READ_PTS)
    out["theta"] = dict(zip(READ_PTS, th))
    out["I"] = []
    ref_tab = ({**I55_THETA, **GEL_NEW} if arm == "gel" else TI_THETA if arm == "ti" else {})
    for (r, z), v in zip(READ_PTS, th):
        cf = QP * _phi(r * MM, z * MM, T_PRINT)
        rel = abs(v / cf - 1.0)
        if (r, z) in GATE_PTS:
            out["I"].append(rel)
        ref = ref_tab.get((r, z))
        cmp = (f"ref {ref:.6e}, |diff| {abs(v - ref):.3e} (6th-digit unit {_sixth(ref):.0e})"
               if ref is not None else "no reproduction on this arm")
        _say(f"[imp8-i81] ({r:g}, {z:g}) mm t*={T_PRINT:g} s arm {arm}: theta {v:.6e} K/W; {cmp}; vs "
             f"q'/(4 pi k) Phi {cf:.6e}: rel {rel:.3e}")
    out["move"] = out["theta"][(2.0, 0.0)] / I55_THETA[(2.0, 0.0)] - 1.0
    _say(f"[imp8-i81] (2, 0) mm theta/theta_55 - 1 = {out['move']:.4e} (control floor "
         f"{CONTROL_FLOOR:.0e}, asserted on arm hole)")

    # 2. I75's probe batch on the final field and on each snapshot, unchanged.
    fields = {N_MAIN: res.temperature, **{j: res.snapshots[j] for j in SNAP_STEPS}}
    probe = {j: _theta_pts(f, PROBE_PTS) for j, f in sorted(fields.items())}
    out["probe"] = probe
    duty = [_superpose({j: probe[j][i] for j in probe}) for i in range(len(PROBE_PTS))]
    out["probe_duty"] = duty
    ref900, refD = I75_PROBE[arm]
    for i, (r, z) in enumerate(PROBE_PTS):
        _say(f"[imp8-i81] (R75) probe ({r:g}, {z:+g}) mm arm {arm}: theta(900 s) {probe[N_MAIN][i]:.6e} "
             f"vs I75 {ref900[i]:.6e} |diff| {abs(probe[N_MAIN][i] - ref900[i]):.3e} (unit "
             f"{_sixth(ref900[i]):.0e}); D=0.5 {duty[i]:.6e} vs I75 {refD[i]:.6e} |diff| "
             f"{abs(duty[i] - refD[i]):.3e} (unit {_sixth(refD[i]):.0e})")

    # 3. New: the ladder batch on the final field and on each snapshot.
    ladder = {j: _theta_pts(f, LADDER_PTS) for j, f in sorted(fields.items())}
    out["ladder"] = ladder
    by = {j: _by_z(probe, ladder, j) for j in fields}
    out["by"] = by
    dutyz = {z: _superpose({j: by[j][z] for j in by}) for z in by[N_MAIN]}
    out["dutyz"] = dutyz
    for za in Z_ALL:
        tag = (" (I75's probe)" if za == 16.0 else " (the clinical probe)" if za == 17.0 else "")
        cf = QP * _phi(0.0, za * MM, T_PRINT)
        cfd = _duhamel(0.0, za * MM, T_PRINT, DUTY)
        for z in (za, -za):
            v = by[N_MAIN][z]
            series = ", ".join(f"{j * dt:g} s {by[j][z]:.6e}" for j in SNAP_STEPS)
            ref16 = by[N_MAIN][16.0 if z > 0 else -16.0]
            _say(f"[imp8-i81] ladder (0, {z:+g}) mm{tag} arm {arm}: theta(900 s) {v:.6e} K/W; {series}; "
                 f"theta(900)/theta(675) {v / by[30][z]:.6f}; D=0.5 superposed theta(900 s) "
                 f"{dutyz[z]:.6e} K/W; theta(z)/theta(16) {v / ref16:.6f}; infinite medium, line "
                 f"source, no rod: QP*Phi {cf:.6e} (FEM/cf {v / cf:.6f}), _duhamel(D=0.5, PERIOD "
                 f"{PERIOD:g} s) {cfd:.6e} (FEM/cf {dutyz[z] / cfd:.6f})")
        vp, vm = by[N_MAIN][za], by[N_MAIN][-za]
        dp, dm = dutyz[za], dutyz[-za]
        mean = 0.5 * (vp + vm)
        dmean = 0.5 * (dp + dm)
        _say(f"[imp8-i81] ladder |z|={za:g} mm{tag} arm {arm}: +/- mean theta(900 s) {mean:.6e} K/W, "
             f"|+ - -| {abs(vp - vm):.3e} (mesh asymmetry); D=0.5 mean {dmean:.6e}, |+ - -| "
             f"{abs(dp - dm):.3e}; mean FEM/QP*Phi {mean / cf:.6f}; D mean FEM/_duhamel "
             f"{dmean / cfd:.6f}; I80 line record {I80_LINE[za]:.7g} (printed beside only)")
    for s in (1.0, -1.0):
        d = by[N_MAIN]
        _say(f"[imp8-i81] control sign {s:+g} arm {arm}: theta 16/17/18/20 {d[16 * s]:.6e} / "
             f"{d[17 * s]:.6e} / {d[18 * s]:.6e} / {d[20 * s]:.6e}; first drop 1 - theta(17)/theta(16) "
             f"{1.0 - d[17 * s] / d[16 * s]:.6f} (band > {INTERIOR_BAND:g}, < 1; infinite line 0.2941583)")
    _say(f"[imp8-i81] (E) {out['E']:.3e} (band {E_BAND:.0e}); (S0) {out['S0']:.3e} (band {B_BAND:.0e}); "
         f"cells {n} (W0C {W0C_CELLS}); non-finite coefficient dofs {n_bad}")
    return out


def test_w_preconditions(ladder_run):
    assert ladder_run["n"] == W0C_CELLS, ladder_run["n"]
    assert ladder_run["n_body"] == BODY_CELLS, ladder_run["n_body"]
    assert ladder_run["n_bad"] == 0
    assert ladder_run["snap_keys"] == list(SNAP_STEPS), ladder_run["snap_keys"]
    assert ladder_run["snap_t"] == [225.0, 450.0, 675.0], ladder_run["snap_t"]
    assert ladder_run["snap_t"][-1] == T2 and SNAP_STEPS[-1] == STEP_T2


def test_w_e_balance(ladder_run):
    assert ladder_run["E"] <= E_BAND, ladder_run["E"]


def test_w_s0_source_power(ladder_run):
    assert ladder_run["S0"] <= B_BAND, ladder_run["S0"]


def test_w_probe_rises_675_to_900(ladder_run):
    p = ladder_run["probe"]
    for i in range(len(PROBE_PTS)):
        assert p[N_MAIN][i] > p[30][i], (PROBE_PTS[i], p[N_MAIN][i], p[30][i])


def test_w_r75_reproduces_i75_probe(ladder_run):
    ref900, refD = I75_PROBE[ladder_run["arm"]]
    p900 = ladder_run["probe"][N_MAIN]
    duty = ladder_run["probe_duty"]
    for i in range(len(PROBE_PTS)):
        assert abs(p900[i] - ref900[i]) <= _sixth(ref900[i]), (PROBE_PTS[i], p900[i], ref900[i])
        assert abs(duty[i] - refD[i]) <= _sixth(refD[i]), (PROBE_PTS[i], duty[i], refD[i])


def test_w_control_ladder_falls(ladder_run):
    d = ladder_run["by"][N_MAIN]
    for s in (1.0, -1.0):
        vals = [d[z * s] for z in Z_ALL]
        assert all(np.isfinite(vals)), vals
        for a, b in zip(vals, vals[1:]):
            assert a > b, (s, vals)
        drop = 1.0 - d[17.0 * s] / d[16.0 * s]
        assert INTERIOR_BAND < drop < 1.0, (s, drop)


def test_w2_gel_reproduces_i73(ladder_run):
    if ladder_run["arm"] != "gel":
        pytest.skip("arm gel only: the anchor")
    for pt, ref in {**I55_THETA, **GEL_NEW}.items():
        v = ladder_run["theta"][pt]
        assert abs(v - ref) <= _sixth(ref), (pt, v, ref)


def test_w2_gel_interior_against_line_source(ladder_run):
    if ladder_run["arm"] != "gel":
        pytest.skip("arm gel only: I55's (I)")
    assert max(ladder_run["I"]) <= INTERIOR_BAND, ladder_run["I"]


def test_w1_ti_reproduces_i73(ladder_run):
    if ladder_run["arm"] != "ti":
        pytest.skip("arm ti only: the route control")
    for pt, ref in TI_THETA.items():
        v = ladder_run["theta"][pt]
        assert abs(v - ref) <= _sixth(ref), (pt, v, ref)


def test_w3_hole_control_moves(ladder_run):
    if ladder_run["arm"] != "hole":
        pytest.skip("arm hole only: the coefficient control")
    assert abs(ladder_run["move"]) > CONTROL_FLOOR, ladder_run["move"]
