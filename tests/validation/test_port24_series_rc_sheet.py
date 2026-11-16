"""`PORT-24` step 2 (§9 item 309, chain T26 prerequisite) -- the in-model
series-RC sheet: a 0.2 Ω ESR on a lumped capacitor sheet, on the F-small exact
four-port at 64 MHz, degree 1, κ off, the copper (Leontovich) column.

The sheet is ``LumpedSheetPortSpec(port_impedance_ohm=series_rlc_impedance(f,
r_ohm=r, l_h=0, c_f=C_F))`` on the exact CAD rectangles of
``test_ans7_exact_strip._exact_geo`` (``width_correction_kappa=None``).  Nothing
here is new physics: the sheet law already takes a complex ``Z_p``; this module
asserts that an ESR-carrying capacitor sheet behaves as its circuit element.

Arms (``FEM_EM_PORT24_ARM``; unset ⇒ skip; anything else raises), one harness
window each:

* ``r02`` -- every sheet carries r = 0.2 Ω.  (I) per sheet i, the field-solved
  3×3 with sheet i terminated in its series-RC equals
  ``reduce_terminated_ports(S₄, 50, {i: series_rlc_impedance(f, r_ohm=0.2,
  l_h=0, c_f=C)})`` at ``REDUCTION_BAND`` (`PORT-14`'s band, imported).
* ``r0`` -- r = 0: (I) as above at r = 0, and (Z) each sheet's terminated 3×3
  reproduces the capacitor-only sheet's (``z_cap`` of the ``ans:7`` circuit
  layer, the sheet Route A places today) at ``LEG_D0_REPRODUCTION_BAND``.
* ``plant`` -- (Pl), predicted, asserted: sheet P2 at 0.4 Ω, the others 0.2;
  P2's identity held against the 0.2 Ω reduction must exceed
  ``REDUCTION_BAND``; the three unplanted must pass.  A planted residual inside
  the band fails here -- the identity cannot see r at this C (never a pass).

Every arm: (P) the exact discrete power identity on one P1 drive with P2–P4
terminated in the arm's series-RC, the sheets' own dissipation counted, at
``DISCRETE_IDENTITY_RTOL``.  (Dz) the same residual with the RC sheets'
dissipation omitted -- printed, never asserted.  Q(r)/Q(0) is printed against
1/Q = 1/Q_rest + 1/Q_C, Q_C = 1/(ωCr), never asserted.

Not claimed: a measured ESR, frequency dependence, any AED re-run.
"""

from __future__ import annotations

import dataclasses
import os
import time

import numpy as np
import pytest
import ufl
from mpi4py import MPI

from dolfinx import fem

from fem_em_solver.core.cavity import surface_resistance_ohm
from fem_em_solver.ports.circuit import reduce_terminated_ports
from fem_em_solver.ports.lumped import run_lumped_sheet_port_case, series_rlc_impedance
from fem_em_solver.ports.sparameters import _power_waves, run_n_port_sparameter_sweep

from tests.complex_mode import complex_only
from tests.validation.test_ans7_exact_strip import C_F, _ANS7, _exact_geo
from tests.validation.test_birdcage_power_identity import (
    DISCRETE_IDENTITY_RTOL,
    _sheet_field_dissipation_w,
    _source_power_w,
)
from tests.validation.test_port_birdcage_four_port import LEG_D0_REPRODUCTION_BAND
from tests.validation.test_port_lumped_narrowed_sheet import GATED_WIDTH_FRACTION
from tests.validation.test_port_lumped_rlc_termination import REDUCTION_BAND
from tests.validation.test_port_package_sparameters import REFERENCE_IMPEDANCE_OHM

ARM_ENV = "FEM_EM_PORT24_ARM"
ARMS = ("r02", "r0", "plant")
FREQUENCY_HZ = 64.0e6
DEGREE = 1
ESR_OHM = 0.2          # operator ESR ruling 2026-10-06 (§9)
PLANT_OHM = 0.4
PLANT_INDEX = 1        # P2: a capacitor sheet in the Route A configuration
POWER_DRIVEN = "P1"


def _arm():
    raw = os.environ.get(ARM_ENV, "").strip()
    if not raw:
        return None
    if raw not in ARMS:
        raise ValueError(f"{ARM_ENV}={raw!r}: one of {ARMS}")
    return raw


def _say(msg):
    if MPI.COMM_WORLD.rank == 0:
        print(msg, flush=True)


def _arm_r(arm):
    """Per-sheet r (index 0..3) for the arm, and the reference r of (I)."""
    if arm == "r0":
        return [0.0] * 4, 0.0
    r = [ESR_OHM] * 4
    if arm == "plant":
        r[PLANT_INDEX] = PLANT_OHM
    return r, ESR_OHM


def _rc(r):
    return series_rlc_impedance(FREQUENCY_HZ, r_ohm=r, l_h=0.0, c_f=C_F)


def _terminated_three_port(problem, port_defs, specs, tags_f, extra, index, z_term):
    """Sheet ``index`` carries ``z_term``; the 3×3 over the kept ports (ascending)."""
    sp = [dataclasses.replace(s, port_impedance_ohm=z_term) if i == index else s
          for i, s in enumerate(specs)]
    kept = [p for i, p in enumerate(port_defs) if i != index]
    z0 = float(REFERENCE_IMPEDANCE_OHM)
    s3 = np.zeros((3, 3), dtype=np.complex128)
    for col, drv in enumerate(kept):
        res = run_lumped_sheet_port_case(problem, port_defs, sp, facet_tags=tags_f,
                                         driven_port_id=drv.port_id, verbose=False,
                                         extra_bilinear_terms=extra, degree=DEGREE)
        d = res.responses[drv.port_id]
        a_d, _ = _power_waves(d.voltage_v, d.current_a, z0)
        assert abs(a_d) > 0.0
        for row, rcv in enumerate(kept):
            rr = res.responses[rcv.port_id]
            _, b = _power_waves(rr.voltage_v, rr.current_a, z0)
            s3[row, col] = b / a_d
    return s3


def _rel(a, b):
    return float(np.linalg.norm(a - b) / np.linalg.norm(b))


@pytest.fixture(scope="module")
def port24():
    arm = _arm()
    if arm is None:
        pytest.skip(f"{ARM_ENV} unset")
    return _run_port24_arm(arm)


def _run_port24_arm(arm, return_fields=False):
    """The fixture body, lifted (EX-105): one arm, solved and accounted.

    ``return_fields=True`` (default off) additionally returns the mesh, tags
    and field for the ``ports:21`` example's figure and XDMF.
    """
    from tests.validation import test_th14_birdcage_copper as th14

    comm =MPI.COMM_WORLD
    t0 = time.perf_counter()
    geo, specs, cells = _exact_geo(th14, comm, GATED_WIDTH_FRACTION)
    msh, cell_tags, facet_tags, outer_tags, tags_f, port_defs = geo
    assert all(s.width_correction_kappa is None for s in specs)
    omega = 2.0 * np.pi * FREQUENCY_HZ
    problem = th14._problem(msh, cell_tags, FREQUENCY_HZ, outer_tags, (th14.OUTER_BOX_TAG,))
    z_s = (1.0 + 1.0j) * surface_resistance_ohm(omega, th14.COPPER_SIGMA)
    extra = [th14._leontovich_term(msh, facet_tags, z_s, omega)]
    t_mesh = time.perf_counter() - t0
    _say(f"\n[PORT-24 s2] arm {arm}: exact F-small hole mesh {cells} cells "
         f"({t_mesh:.1f} s), f = {FREQUENCY_HZ:.6e} Hz, degree {DEGREE}, kappa off, Cu "
         f"Leontovich, C = {C_F!r} F, -n {comm.size}")

    t1 = time.perf_counter()
    sw = run_n_port_sparameter_sweep(problem, port_defs, lumped_sheet_ports=specs,
                                     lumped_sheet_facet_tags=tags_f,
                                     extra_bilinear_terms=extra, degree=DEGREE)
    s4 = np.asarray(sw.s_matrix, dtype=np.complex128)
    _say(f"[PORT-24 s2] 50 ohm 4x4 in {time.perf_counter() - t1:.1f} s; "
         f"||S-S^T||/||S|| = {np.linalg.norm(s4 - s4.T) / np.linalg.norm(s4):.3e}")

    r_list, r_ref = _arm_r(arm)
    z_ref = _rc(r_ref)
    z0 = float(REFERENCE_IMPEDANCE_OHM)
    ident = []
    for i in range(4):
        t2 = time.perf_counter()
        z_i = _rc(r_list[i])
        meas = _terminated_three_port(problem, port_defs, specs, tags_f, extra, i, z_i)
        pred = reduce_terminated_ports(s4, z0, {i: z_ref})
        row = {"i": i, "r": r_list[i], "z": z_i, "meas": meas, "pred": pred,
               "res": _rel(meas, pred)}
        if arm == "r0":
            zc = _ANS7.z_cap(omega, C_F)
            cap = _terminated_three_port(problem, port_defs, specs, tags_f, extra, i, zc)
            row["zc"] = zc
            row["cap_res"] = _rel(meas, cap)
        ident.append(row)
        _say(f"[PORT-24 s2] (I) sheet P{i + 1}: r = {r_list[i]:g} ohm, Z = "
             f"{z_i.real:+.9e}{z_i.imag:+.9e}j ohm vs reduction at "
             f"series_rlc_impedance(r={r_ref:g}, L=0, C) = {z_ref.real:+.9e}{z_ref.imag:+.9e}j; "
             f"||S'_meas - S'_pred||_F/||S'_pred||_F = {row['res']:.6e} "
             f"({row['res'] / REDUCTION_BAND:.4f} x REDUCTION_BAND {REDUCTION_BAND:.0e})"
             + (f"; (Z) vs capacitor-only z_cap {row['zc'].real:+.3e}{row['zc'].imag:+.9e}j: "
                f"{row['cap_res']:.6e} (band {LEG_D0_REPRODUCTION_BAND:.0e})" if arm == "r0" else "")
             + f"; {time.perf_counter() - t2:.1f} s")

    # (P) power identity: P1 driven at 50 ohm, P2-P4 in the arm's series-RC.
    t3 = time.perf_counter()
    sp = [s if i == 0 else dataclasses.replace(s, port_impedance_ohm=_rc(r_list[i]))
          for i, s in enumerate(specs)]
    _res, fields = run_lumped_sheet_port_case(
        problem, port_defs, sp, facet_tags=tags_f, driven_port_id=POWER_DRIVEN,
        verbose=False, return_fields=True, extra_bilinear_terms=extra, degree=DEGREE)
    e = fields.e_complex
    dx = ufl.Measure("dx", domain=msh, subdomain_data=cell_tags)
    p_vol = float(np.real(comm.allreduce(fem.assemble_scalar(fem.form(
        0.5 * fields.sigma_field * ufl.inner(e, e) * dx)), op=MPI.SUM)))
    acct = {"mesh": msh, "facet_tags": tags_f}
    sheets = {s.port_id: s.sheet(driven=(s.port_id == POWER_DRIVEN)) for s in sp}
    p_sheet = {pid: _sheet_field_dissipation_w(acct, sh, e, omega) for pid, sh in sheets.items()}
    p_src = float(_source_power_w(acct, sheets[POWER_DRIVEN], e, omega))
    p_surf = th14._surface_loss_w(msh, facet_tags, z_s, e, comm)
    p_rc = sum(v for pid, v in p_sheet.items() if pid != POWER_DRIVEN)
    p_all = sum(p_sheet.values())
    resid = abs(p_src - p_all - p_vol - p_surf) / abs(p_src)
    resid_dz = abs(p_src - (p_all - p_rc) - p_vol - p_surf) / abs(p_src)
    p_rest = p_vol + p_surf          # dissipation other than the RC sheets
    share = p_rc / (p_rc + p_rest)
    q_c = 1.0 / (omega * C_F * ESR_OHM)
    _say(f"[PORT-24 s2] (P) P1 drive, P2-P4 series-RC r = {r_list[1:]}: P_src {p_src:.9e} W, "
         f"P_sheets {p_all:.9e} W (P1 {p_sheet['P1']:.9e}, RC "
         + ", ".join(f"{k} {v:.9e}" for k, v in p_sheet.items() if k != 'P1')
         + f"), P_vol {p_vol:.9e} W, P_surf {p_surf:.9e} W; residual/P_src {resid:.3e} "
         f"(ASSERTED <= {DISCRETE_IDENTITY_RTOL:g}); {time.perf_counter() - t3:.1f} s")
    _say(f"[PORT-24 s2] (Dz) RC-sheet dissipation omitted: residual/P_src {resid_dz:.6e} "
         f"(printed, never asserted; PREDICTED ~ ESR loss share Q_rest/(Q_rest + Q_C)); "
         f"measured ESR share P_RC/(P_RC + P_vol + P_surf) = {share:.6e}")
    q_ratio = None
    q_rest = None
    if p_rc > 0.0:
        q_ratio = p_rest / (p_rest + p_rc)
        q_rest_over_qc = p_rc / p_rest
        q_rest = q_rest_over_qc * q_c
        _say(f"[PORT-24 s2] Q(r)/Q(0) at fixed stored energy = P_rest/(P_rest + P_RC) = "
             f"{p_rest / (p_rest + p_rc):.6e}; 1/Q = 1/Q_rest + 1/Q_C with Q_C = 1/(omega C r) "
             f"= {q_c:.6e} => Q_rest = {q_rest_over_qc * q_c:.6e} (printed only)")
    else:
        _say("[PORT-24 s2] Q(r)/Q(0) = 1 (r = 0: the RC sheets dissipate "
             f"{p_rc:.3e} W); Q_C at 0.2 ohm would be {q_c:.6e} (printed only)")
    _say(f"[PORT-24 s2] arm {arm} total {time.perf_counter() - t0:.1f} s, -n {comm.size}")
    out = {"arm": arm, "ident": ident, "resid": resid, "resid_dz": resid_dz,
           "p_rc": p_rc, "s4": s4}
    if return_fields:
        out.update({
            "cells": cells, "t_mesh": t_mesh, "omega": omega, "p_src": p_src,
            "p_sheet": p_sheet, "p_vol": p_vol, "p_surf": p_surf, "p_rest": p_rest,
            "share": share, "q_c": q_c, "q_ratio": q_ratio, "q_rest": q_rest,
            "msh": msh, "cell_tags": cell_tags, "facet_tags": facet_tags,
            "specs": sp, "e": e,
        })
    return out


@complex_only
def test_port24_series_rc_sheet(port24):
    arm = port24["arm"]
    rows = port24["ident"]
    assert len(rows) == 4
    if arm == "plant":
        planted = rows[PLANT_INDEX]
        others = [r for r in rows if r["i"] != PLANT_INDEX]
        _say(f"[PORT-24 s2] (Pl) planted P{PLANT_INDEX + 1} residual {planted['res']:.6e} = "
             f"{planted['res'] / REDUCTION_BAND:.3f} x REDUCTION_BAND (ASSERTED > band; "
             f"predicted ceiling ~5.3e-3); unplanted max "
             f"{max(r['res'] for r in others):.6e} (ASSERTED <= band)")
        assert planted["res"] > REDUCTION_BAND, (
            "(Pl): the planted 0.4 ohm sheet passes the 0.2 ohm identity -- the identity "
            f"cannot see r at this C (residual {planted['res']:.6e}); report and stop")
        for r in others:
            assert r["res"] <= REDUCTION_BAND, (r["i"], r["res"])
    else:
        for r in rows:
            assert r["res"] <= REDUCTION_BAND, (r["i"], r["res"])
        if arm == "r0":
            for r in rows:
                assert r["cap_res"] <= LEG_D0_REPRODUCTION_BAND, (r["i"], r["cap_res"])
    assert port24["resid"] <= DISCRETE_IDENTITY_RTOL, port24["resid"]
