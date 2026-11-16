"""`IMP-8` step 17 (§9 item 881, §10 I75): the F2182-shaped probe readout on the
gel / Ti / hole bracket.

Brief: docs/planning/chunks/IMP-8.md, "§9 item 881".  I73's fixture unchanged
(I55's graded gel box, real build, degree 2, CN to t* = 900 s, I55's uniform rod
source); the body cells are I73's (the positive DG0 dofs of q, as cell indices
through the DG0 dofmap, owned + ghost) and get their coefficients through I66's
``_coefficients(msh, body_idx, arm)`` (``ti`` / ``gel`` / ``hole``, ε =
``EPS_HOLE``).  Snapshots at steps 10 / 20 / 30 (225 / 450 / 675 s).

Arm knob ``FEM_EM_IMP8_I75_ARM`` = ``ti`` (default) / ``gel`` / ``hole``.

Asserted (bands imported): (E) ≤ ``E_BAND``; (S0) ≤ ``B_BAND``; cells =
``W0C_CELLS``, body cells = ``BODY_CELLS``, no non-finite coefficient dofs; arm
``gel`` reproduces I73's gel θ (I55's five plus (2, 15) mm) to one unit of the
sixth significant digit and I55's (I) at ``INTERIOR_BAND``; arm ``ti`` reproduces
I73's six Ti θ the same way; arm ``hole`` moves θ(2, 0) off I55's by more than
``CONTROL_FLOOR``; every arm, the re-posed weekly control θ_probe(900 s) >
θ_probe(675 s) at both probes (factor 1, strict).

Printed only (🟡 record, rule (e)): the probe θ at (0, ±16) mm, the D = 0.5 duty
read by superposition θ(900) − θ(675) + θ(450) − θ(225) (not a switched solve),
the infinite-medium line-source closed forms, and on arm ``gel`` the superposed
duty at ``GATE_PTS`` beside I59's ``_duhamel`` and I60's direct FEM reads.
"""

from __future__ import annotations

import math
import os
import time

import numpy as np
from dolfinx import fem
from mpi4py import MPI

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

ARM_ENV = "FEM_EM_IMP8_I75_ARM"
ARMS = ("ti", "gel", "hole")
Z_P = 16.0  # mm, the probe: on the axis, 1 mm beyond each rod end (L/2 = 15 mm)
PROBE_PTS = ((0.0, Z_P), (0.0, -Z_P))  # (r, z) mm
SNAP_STEPS = (10, 20, 30)
DUTY = 0.5
#: I73's gel point, K/W (`20261106T103928Z_IMP-8-item856-W2.log:118`).
GEL_NEW = {(2.0, 15.0): 9.482226e00}
#: I73's Ti θ at ``READ_PTS``, K/W (`20261106T103149Z_IMP-8-item856-W1.log:113–118`).
TI_THETA = dict(zip(READ_PTS, (4.583697e00, 1.098572e00, 2.774057e00, 7.112581e-01,
                               1.735615e01, 1.041661e01)))
#: I60's direct D = 0.5 FEM reads at ``GATE_PTS``, K/W
#: (`20261104T044429Z_IMP-8-I60-W1-item799.log:112–115`), printed only.
I60_DIRECT = dict(zip(GATE_PTS, (2.240893e00, 6.683594e-01, 1.366131e00, 4.289727e-01)))


def _arm():
    raw = os.environ.get(ARM_ENV, "").strip()
    if raw == "":
        return "ti"
    if raw in ARMS:
        return raw
    raise ValueError(f"{ARM_ENV} must be unset or one of {ARMS}, got {raw!r}")


def test_k0_arm_knob(monkeypatch):
    monkeypatch.delenv(ARM_ENV, raising=False)
    assert _arm() == "ti"
    for a in ARMS:
        monkeypatch.setenv(ARM_ENV, a)
        assert _arm() == a
    for bad in ("Ti", "pec", "1", "HOLE"):
        monkeypatch.setenv(ARM_ENV, bad)
        with pytest.raises(ValueError):
            _arm()


def _sixth(ref):
    return 10.0 ** (math.floor(math.log10(abs(ref))) - 5)


def _body_cells(msh, q):
    """I73's body (q's positive DG0 dofs) as cell indices, owned + ghost; the
    allreduced owned count."""
    W = q.function_space
    tdim = msh.topology.dim
    n_loc = msh.topology.index_map(tdim).size_local
    n_all = n_loc + msh.topology.index_map(tdim).num_ghosts
    dofs = np.asarray(W.dofmap.list).reshape(-1)[:n_all]
    body = np.flatnonzero(np.real(q.x.array[dofs]) > 0.0).astype(np.int32)
    n_body = msh.comm.allreduce(int(np.count_nonzero(body < n_loc)), op=MPI.SUM)
    return body, n_body


def _superpose(th):
    """D = 0.5 at PERIOD 450 s from the continuous step response, from rest:
    θ_D(900) = θ(900) − θ(675) + θ(450) − θ(225)."""
    return th[40] - th[30] + th[20] - th[10]


@pytest.fixture(scope="module")
def probe_run(box):
    _need("transient")
    arm = _arm()
    msh, ct, ft, n = box
    q, out = _source_and_s0(msh)
    body_idx, n_body = _body_cells(msh, q)
    k, rho_c, n_bad = _coefficients(msh, body_idx, arm)
    out.update({"arm": arm, "n": n, "n_body": n_body, "n_bad": n_bad})
    _say(f"[imp8-i75] arm {arm}; ROD_K {ROD_K:g}, ROD_RHO_C {ROD_RHO_C:.6g} (nominal); K_GEL {K_GEL:g}, "
         f"RHO_C {RHO_C:.6g}; EPS_HOLE {EPS_HOLE:.0e}; ROD_RADIUS {ROD_RADIUS:g} m, L {L:g} m")
    _say(f"[imp8-i75] body cells (positive q dofs -> cells, owned, allreduced) {n_body} (BODY_CELLS "
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
    _say(f"[imp8-i75] CN N={N_MAIN} dt={dt:g} s: t_call {out['t_call']:.1f} s; n_factorizations "
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
        _say(f"[imp8-i75] ({r:g}, {z:g}) mm t*={T_PRINT:g} s arm {arm}: theta {v:.6e} K/W; {cmp}; vs "
             f"q'/(4 pi k) Phi {cf:.6e}: rel {rel:.3e}")
    out["move"] = out["theta"][(2.0, 0.0)] / I55_THETA[(2.0, 0.0)] - 1.0
    _say(f"[imp8-i75] (2, 0) mm theta/theta_55 - 1 = {out['move']:.4e} (control floor "
         f"{CONTROL_FLOOR:.0e}, asserted on arm hole)")

    # 2. The probe batch on the final field and on each snapshot.
    fields = {N_MAIN: res.temperature, **{j: res.snapshots[j] for j in SNAP_STEPS}}
    probe = {j: _theta_pts(f, PROBE_PTS) for j, f in sorted(fields.items())}
    out["probe"] = probe
    duty = [_superpose({j: probe[j][i] for j in probe}) for i in range(len(PROBE_PTS))]
    for i, (r, z) in enumerate(PROBE_PTS):
        series = ", ".join(f"{j * dt:g} s {probe[j][i]:.6e}" for j in SNAP_STEPS)
        _say(f"[imp8-i75] probe ({r:g}, {z:+g}) mm arm {arm}: theta(900 s) {probe[N_MAIN][i]:.6e} K/W; "
             f"{series}; theta(900)/theta(675) {probe[N_MAIN][i] / probe[30][i]:.6f}; D=0.5 superposed "
             f"theta(900 s) {duty[i]:.6e} K/W")
    p900 = probe[N_MAIN]
    _say(f"[imp8-i75] probe arm {arm}: mean {0.5 * (p900[0] + p900[1]):.6e} K/W, |+16 - -16| "
         f"{abs(p900[0] - p900[1]):.3e} (mesh asymmetry); D=0.5 mean {0.5 * (duty[0] + duty[1]):.6e}")
    for z in (Z_P, -Z_P):
        _say(f"[imp8-i75] infinite medium, line source, no rod: (0, {z:+g}) mm QP*Phi(900 s) "
             f"{QP * _phi(0.0, z * MM, T_PRINT):.6e} K/W; _duhamel(D=0.5, PERIOD {PERIOD:g} s) "
             f"{_duhamel(0.0, z * MM, T_PRINT, DUTY):.6e} K/W")

    # 3. Arm gel: the superposed duty at GATE_PTS beside I59 and I60 (printed only).
    if arm == "gel":
        gate = {j: _theta_pts(f, GATE_PTS) for j, f in sorted(fields.items())}
        for i, (r, z) in enumerate(GATE_PTS):
            sup = _superpose({j: gate[j][i] for j in gate})
            duh = _duhamel(r * MM, z * MM, T_PRINT, DUTY)
            ref = I60_DIRECT[(r, z)]
            _say(f"[imp8-i75] gel duty ({r:g}, {z:g}) mm: superposed {sup:.6e} K/W; I59 _duhamel "
                 f"{duh:.6e} (rel {sup / duh - 1.0:+.3e}); I60 direct FEM {ref:.6e} (restart vs no "
                 f"restart: diff {sup - ref:+.3e}, rel {sup / ref - 1.0:+.3e}; printed, rule (e))")
    _say(f"[imp8-i75] (E) {out['E']:.3e} (band {E_BAND:.0e}); (S0) {out['S0']:.3e} (band {B_BAND:.0e}); "
         f"cells {n} (W0C {W0C_CELLS}); non-finite coefficient dofs {n_bad}")
    return out


def test_w_preconditions(probe_run):
    assert probe_run["n"] == W0C_CELLS, probe_run["n"]
    assert probe_run["n_body"] == BODY_CELLS, probe_run["n_body"]
    assert probe_run["n_bad"] == 0
    assert probe_run["snap_keys"] == list(SNAP_STEPS), probe_run["snap_keys"]
    assert probe_run["snap_t"] == [225.0, 450.0, 675.0], probe_run["snap_t"]
    assert probe_run["snap_t"][-1] == T2 and SNAP_STEPS[-1] == STEP_T2


def test_w_e_balance(probe_run):
    assert probe_run["E"] <= E_BAND, probe_run["E"]


def test_w_s0_source_power(probe_run):
    assert probe_run["S0"] <= B_BAND, probe_run["S0"]


def test_w_probe_rises_675_to_900(probe_run):
    """The re-posed weekly control: monotone heating under a continuous source."""
    p = probe_run["probe"]
    for i in range(len(PROBE_PTS)):
        assert p[N_MAIN][i] > p[30][i], (PROBE_PTS[i], p[N_MAIN][i], p[30][i])


def test_w2_gel_reproduces_i73(probe_run):
    if probe_run["arm"] != "gel":
        pytest.skip("arm gel only: the anchor")
    for pt, ref in {**I55_THETA, **GEL_NEW}.items():
        v = probe_run["theta"][pt]
        assert abs(v - ref) <= _sixth(ref), (pt, v, ref)


def test_w2_gel_interior_against_line_source(probe_run):
    if probe_run["arm"] != "gel":
        pytest.skip("arm gel only: I55's (I)")
    assert max(probe_run["I"]) <= INTERIOR_BAND, probe_run["I"]


def test_w1_ti_reproduces_i73(probe_run):
    if probe_run["arm"] != "ti":
        pytest.skip("arm ti only: the route control")
    for pt, ref in TI_THETA.items():
        v = probe_run["theta"][pt]
        assert abs(v - ref) <= _sixth(ref), (pt, v, ref)


def test_w3_hole_control_moves(probe_run):
    if probe_run["arm"] != "hole":
        pytest.skip("arm hole only: the coefficient control")
    assert abs(probe_run["move"]) > CONTROL_FLOOR, probe_run["move"]
