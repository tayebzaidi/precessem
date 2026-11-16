"""Example (`EX-120`): tune, then load -- ``Study.tune(k=1)``'s ``C_mode`` placed on three
gap sheets through ``CapacitorSheets``, the reduced network equals the circuit reduction,
and the field of the loaded coil.

`WF-16` step 6(b) (item 463) gated ``CapacitorSheets`` at tune's ``C_mode`` on the
exact-rectangle PEC hole F-small birdcage (degree 1, 64 MHz, ``-n 2``, kappa off).  No
example placed capacitors in the solved model through the public ``Study`` API.  This one:

* solves the 50 ohm four-port Study, runs ``tune(k=1)`` and ``tune(k=2)``, reproduces the
  gate's anchors (R0), (a');
* feeds the in-run ``C_mode`` into ``CapacitorSheets(C, 0.0, on=("P2", "P3", "P4"))`` and
  solves the loaded (1 x 1) Study; (I-A) the reduced S11 equals the circuit reduction of
  the untuned 4 x 4 (``reduce_terminated_ports``);
* because a loaded Study keeps no fields, solves once more on the result's kept problem
  to get the loaded field: (S=) it reproduces the S the identity read; (P) the exact power
  identity with the sheets' dissipation counted;
* writes the loaded field and an open-gap control to ParaView.

Lossless capacitors only (item 463's record; the capacitor-ESR ruling of 2026-10-06 scopes
r to new fixtures, and lossless records stay valid, labelled lossless).

Every record, band and helper is imported from the gate modules (the ANS-1 rule); record
digits sit only in the example-local ``REC`` dict, each with its ``log:line``.  New values
print at <= 6 significant digits.

Negative controls (asserted): (c1) the k = 2 capacitance misses; (c2) (I-A)'s comparator
reads C; (cO) the capacitors dropped from the field solve (open gaps) move S11;
(K) the refusals.

Run it through the example runner (the ``mri:`` group sources the complex build)::

    ./scripts/run_examples.sh -e mri:24 -n 2 -t 420

Not claimed: an in-model 4 x 4 at C_mode, f0, Q, R_in, a match, absolute S11 / Z_in, ESR,
copper, degree 2, other frequencies, F-human, an external check.
"""

from __future__ import annotations

import gc
import sys
import time
from pathlib import Path

import numpy as np
import ufl
from mpi4py import MPI

from dolfinx import default_scalar_type, fem

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from fem_em_solver.io.paraview_utils import (  # noqa: E402
    adopt_host_ownership,
    write_xdmf_with_tags,
)
from fem_em_solver.ports.lumped import run_lumped_sheet_port_case  # noqa: E402
from fem_em_solver.ports.sparameters import _power_waves  # noqa: E402
from fem_em_solver.post.evaluation import evaluate_vector_field_parallel  # noqa: E402
from fem_em_solver.post.faraday import (  # noqa: E402
    magnetic_flux_density_from_e,
    project_to_cg1,
)
from fem_em_solver.post.setup_figure import write_setup_figure  # noqa: E402
from fem_em_solver.utils.instrumentation import report_peak_rss  # noqa: E402

import tests.package.test_wf16_study_tune_capacitors as g  # noqa: E402
import tests.validation.test_port_birdcage_four_port as bc4  # noqa: E402
import tests.validation.test_port_lumped_rlc_termination as rlc  # noqa: E402
from tests.validation.test_birdcage_power_identity import (  # noqa: E402
    DISCRETE_IDENTITY_RTOL,
    _sheet_field_dissipation_w,
    _source_power_w,
)
from tests.validation.test_port_birdcage_lumped_column import PHANTOM_CELL_TAG  # noqa: E402

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"
BASENAME = "mri_24_study_tune_then_capacitor_sheets"

LOADED = ("P2", "P3", "P4")

# Gate-log records (docs/testing/logs/20261015T201110Z_WF-16-s6b-W1.log).  Tolerance of
# each is 10^-(d-2) with d its printed significant digits (2026-10-14 10:30 ruling's form).
REC = {
    "cells": 80798,                                      # :939
    "c1": 6.626796421e-11, "c2": 9.873468900e-11,        # :940   (1e-8)
    "z1": complex(2.500348022e-02, 3.752636761e+01),     # :939   (1e-8, per component)
    "a_prime": -2.976552e-04,                            # :950   (1e-5)
    "i_a": 1.391811e-05,                                 # :954   (1e-5)
    "s_a": complex(-2.787323415e-01, 9.590409206e-01),   # :954   (1e-8, |dS|/|S|)
    "c1_factor": 3.2883e-01,                             # :970   (1e-3)
    "c2_factor": 3.5268e-01,                             # :971   (1e-3)
}
R0_PRIOR = (5.052e-11, 3.878e-11, 5.875e-11)             # :949, printed only
ZIN_PRIOR = complex(4.988020e-02, 1.071628e-02)          # :973, printed only


def _p0(comm, text):
    if comm.rank == 0:
        print(text, flush=True)


def _raises(exc_type, fn):
    try:
        fn()
    except exc_type:
        return True
    except Exception as err:  # wrong type is a failure
        raise AssertionError(f"expected {exc_type.__name__}, got {type(err).__name__}: {err}")
    return False


def _field_solve(res, ports, z0):
    """Item 373's (P) call (tests/package/test_wf16_study_capacitor_sheets.py:248-251)."""
    r, f = run_lumped_sheet_port_case(
        res.problem, ports, res.specs, facet_tags=res.facet_tags,
        driven_port_id="P1", verbose=False, return_fields=True,
        extra_bilinear_terms=res.conductor_terms, degree=res.degree)
    resp = r.responses["P1"]
    a, b = _power_waves(resp.voltage_v, resp.current_a, z0)
    return r, f, complex(a), complex(b)


def main() -> None:
    comm = MPI.COMM_WORLD
    started = time.perf_counter()
    wall = {}

    # ---- structural asserts -------------------------------------------------
    assert comm.size == g.ANS6_EXACT_STRIP_RECORD_RANK_WIDTH, (
        f"(W) the degree-1 record is asserted at -n {g.ANS6_EXACT_STRIP_RECORD_RANK_WIDTH}; "
        f"got -n {comm.size}")
    if not np.issubdtype(default_scalar_type, np.complexfloating):
        raise RuntimeError("complex build required (source dolfinx-complex-mode)")
    assert g.REDUCTION_BAND is rlc.REDUCTION_BAND
    assert g.LEG_D0_REPRODUCTION_BAND is bc4.LEG_D0_REPRODUCTION_BAND
    s0 = g._study(g.LumpedPortDrive())
    assert s0.degree == 1 and s0.conductor is None
    assert s0.geometry.as_hole is True and s0.result is None
    z0 = float(g.REFERENCE_IMPEDANCE_OHM)
    assert float(s0.excitation.z0_ohm) == z0
    omega = s0.omega
    band = g.REDUCTION_BAND

    # ---- (K) refusals, before any solve, no mesh ----------------------------
    assert _raises(ValueError, lambda: s0.tune(k=0))
    assert _raises(TypeError, lambda: s0.tune(k=True))
    c_const = g.TH17_STEP3_C_MODE_RECORD_F
    sk = g._study(g.CapacitorSheets(c_const, 0.0, on=LOADED))
    assert _raises(TypeError, lambda: sk.tune())
    assert _raises(TypeError, lambda: sk.b1_plus(comm=comm))
    assert _raises(ValueError, lambda: g._study(
        g.CapacitorSheets(c_const, 0.0, on=g.PORTS)).solve(comm))
    assert s0.result is None and sk.result is None

    # ---- (1) anchor: S0 and the two tune reads ------------------------------
    t0 = time.perf_counter()
    res0 = s0.solve(comm)
    wall["S0"] = time.perf_counter() - t0
    t1 = s0.tune(k=1, comm=comm)
    t2 = s0.tune(k=2, comm=comm)
    assert s0.result is res0
    s0m = np.asarray(res0.s_matrix, dtype=np.complex128)
    means = g._class_means(s0m)
    want = g.ANS6_EXACT_STRIP_CLASS_MEANS_P1["64"]["pec"]
    r0 = {c: abs(means[c] - want[c]) / abs(want[c]) for c in g.CLASSES}
    a_prime = t1.c_mode_farad / g.TH17_STEP3_C_MODE_RECORD_F - 1.0
    z1 = complex(t1.z_k[1])
    rec_c1 = abs(t1.c_mode_farad / REC["c1"] - 1.0)
    rec_c2 = abs(t2.c_mode_farad / REC["c2"] - 1.0)
    rec_z1 = (abs(z1.real / REC["z1"].real - 1.0), abs(z1.imag / REC["z1"].imag - 1.0))
    rec_ap = abs(a_prime / REC["a_prime"] - 1.0)

    # ---- load: S_A through CapacitorSheets at the in-run C_mode -------------
    exc = g.CapacitorSheets(t1.c_mode_farad, 0.0, on=LOADED)
    sA = g._study(exc)
    t0 = time.perf_counter()
    resA = sA.solve(comm)
    wall["SA"] = time.perf_counter() - t0
    zc = exc.impedance_ohm(g.FREQUENCY_64_HZ)
    assert [p.port_id for p in resA.port_definitions] == ["P1"]
    assert resA.capacitor_port_ids == LOADED
    assert resA.capacitor_impedance_ohm == zc
    assert resA.s_matrix.shape == (1, 1)
    assert resA.cells == res0.cells
    assert sA.degree == 1 and sA.conductor is None
    sa11 = complex(resA.s_matrix[0, 0])

    pA = g.reduce_terminated_ports(s0m, z0, {1: zc, 2: zc, 3: zc})
    i_a = abs(sa11 - pA[0, 0]) / abs(pA[0, 0])
    rec_ia = abs(i_a / REC["i_a"] - 1.0)
    rec_sa = abs(sa11 - REC["s_a"]) / abs(REC["s_a"])

    # setup figure (no-op unless FEM_EM_SETUP_FIGURES=1)
    local_tags = np.unique(resA.cell_tags.values)
    all_tags = sorted({int(t) for part in comm.allgather(local_tags.tolist()) for t in part})
    air_tags = tuple(t for t in all_tags if t != PHANTOM_CELL_TAG)
    names = {PHANTOM_CELL_TAG: "gel phantom"}
    for t in air_tags:
        names[t] = "air (outer box PEC)"
    write_setup_figure(
        resA.mesh, resA.cell_tags, FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names=names, hide_tags=air_tags, translucent_tags=(PHANTOM_CELL_TAG,),
        slice_normal=(0.0, 0.0, 1.0), slice_origin=(0.0, 0.0, 0.0),
        title=("mri:24 -- tune then load: CapacitorSheets(C_mode) on P2-P4 of the PEC hole "
               "F-small birdcage, P1 fed (WF-16 step 6b)"),
        comm=comm)

    # ---- negative controls (circuit side) -----------------------------------
    zc2 = g.CapacitorSheets(t2.c_mode_farad, 0.0, on=("P2",)).impedance_ohm(
        g.FREQUENCY_64_HZ)
    zt2 = g._modal_z(g.s_to_z(s0m, z0) + zc2 * np.eye(4))[0]
    c1_val = abs(zt2[1].imag) / abs(t1.z_k[1])
    pA2 = g.reduce_terminated_ports(s0m, z0, {1: zc2, 2: zc2, 3: zc2})
    c2_val = abs(sa11 - pA2[0, 0]) / abs(pA2[0, 0])
    rec_c1f = abs(c1_val / REC["c1_factor"] - 1.0)
    rec_c2f = abs(c2_val / REC["c2_factor"] - 1.0)

    # ---- the field route: one extra solve on the loaded Study's kept problem -
    assert not resA.conductor_terms
    gc.collect()
    comm.Barrier()
    t0 = time.perf_counter()
    rT, fT, aT, bT = _field_solve(resA, resA.all_port_definitions, z0)
    wall["field_tuned"] = time.perf_counter() - t0
    s_eq = abs(bT / aT - sa11) / abs(sa11)

    # (P): copy of item 373's block (test_wf16_study_capacitor_sheets.py:252-264), no P_surf
    msh = resA.mesh
    e = fT.e_complex
    dx = ufl.Measure("dx", domain=msh, subdomain_data=resA.cell_tags)
    p_vol = float(np.real(comm.allreduce(fem.assemble_scalar(fem.form(
        0.5 * fT.sigma_field * ufl.inner(e, e) * dx)), op=MPI.SUM)))
    acct = {"mesh": msh, "facet_tags": resA.facet_tags}
    sheets = {s.port_id: s.sheet(driven=(s.port_id == "P1")) for s in resA.specs}
    p_sheet = {pid: _sheet_field_dissipation_w(acct, sh, e, omega)
               for pid, sh in sheets.items()}
    p_src = float(_source_power_w(acct, sheets["P1"], e, omega))
    p_all = sum(p_sheet.values())
    p_rc = sum(v for pid, v in p_sheet.items() if pid != "P1")
    p_resid = abs(p_src - p_all - p_vol) / abs(p_src)

    # (cO): capacitors dropped -- kept ports only, so P2-P4 are open gaps
    t0 = time.perf_counter()
    rO, fO, aO, bO = _field_solve(resA, resA.port_definitions, z0)
    wall["field_open"] = time.perf_counter() - t0
    s11_open = bO / aO
    co_val = abs(s11_open - sa11) / abs(sa11)

    # ---- fields per unit incident wave --------------------------------------
    t0 = time.perf_counter()
    bt = project_to_cg1(magnetic_flux_density_from_e(fT.e_complex, omega, name="B_tuned_dg0"),
                        name="B_tuned")
    bt.x.array[:] = bt.x.array / aT
    bt.x.scatter_forward()
    bo = project_to_cg1(magnetic_flux_density_from_e(fO.e_complex, omega, name="B_open_dg0"),
                        name="B_open")
    bo.x.array[:] = bo.x.array / aO
    bo.x.scatter_forward()
    wall["projections"] = time.perf_counter() - t0

    origin = np.zeros((1, 3))
    centre = {}
    for key, fn in (("tuned", bt), ("open", bo)):
        vals, valid = evaluate_vector_field_parallel(fn, origin)
        vals = np.asarray(vals).reshape(-1, 3)
        assert bool(np.all(np.asarray(valid, dtype=bool))), f"{key}: origin not found"
        v = vals[0]
        centre[key] = {"b1p": abs(v[0] + 1j * v[1]) / 2.0, "b1m": abs(v[0] - 1j * v[1]) / 2.0,
                       "norm": float(np.linalg.norm(v))}

    space = bt.function_space
    space_s = fem.functionspace(space.mesh, ("Lagrange", 1))

    def _scalar(name, values):
        f = fem.Function(space_s, name=name)
        assert values.shape[0] == f.x.array.shape[0], (values.shape, f.x.array.shape)
        f.x.array[:] = values
        f.x.scatter_forward()
        return f

    nt = bt.x.array.reshape(-1, 3)
    no = bo.x.array.reshape(-1, 3)
    b_real = fem.Function(space, name="B_tuned_real")
    b_real.x.array[:] = np.real(bt.x.array)
    b_real.x.scatter_forward()
    fields = {
        "B_tuned_real": b_real,
        "B_tuned_abs": _scalar("B_tuned_abs", np.sqrt(np.sum(np.abs(nt) ** 2, axis=1))),
        "B1_plus_tuned": _scalar("B1_plus_tuned", np.abs(nt[:, 0] + 1j * nt[:, 1]) / 2.0),
        "B1_minus_tuned": _scalar("B1_minus_tuned", np.abs(nt[:, 0] - 1j * nt[:, 1]) / 2.0),
        "B_open_abs": _scalar("B_open_abs", np.sqrt(np.sum(np.abs(no) ** 2, axis=1))),
    }
    t0 = time.perf_counter()
    OUTPUT_DIR.mkdir(exist_ok=True)
    xdmf_path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_combined", msh, resA.cell_tags, fields, comm=comm)
    adopt_host_ownership(OUTPUT_DIR, comm=comm)
    wall["xdmf"] = time.perf_counter() - t0
    peak = report_peak_rss(comm, label="EX-120 after outputs")

    # ---- prints (no collective inside) --------------------------------------
    zin = z0 * (1.0 + sa11) / (1.0 - sa11) + zc
    zk = np.asarray(t1.z_k)
    if comm.rank == 0:
        print("=" * 72)
        print("EX-120 -- tune then load: CapacitorSheets(C_mode) on P2-P4, lossless")
        print("=" * 72)
        print(f"[study] degree 1, -n {comm.size}, f = {s0.frequency_hz:.6e} Hz, "
              f"z0 = {z0:g} ohm, omega = {omega:.6e} rad/s; lossless capacitors (r = 0)")
        print(f"[(K)] tune(k=0) ValueError, tune(k=True) TypeError, loaded tune() TypeError, "
              f"loaded b1_plus() TypeError, no-port-left solve ValueError; "
              f"results still None")
        print(f"[S0] cells {res0.cells} (record {REC['cells']}), {wall['S0']:.1f} s")
        print("[(R0)] class-mean rel vs 'pec' record: "
              + "  ".join(f"{c} {r0[c]:.3e}" for c in g.CLASSES)
              + f" (band {g.LEG_D0_REPRODUCTION_BAND:.0e}; gate prior "
              + " / ".join(f"{x:.3e}" for x in R0_PRIOR) + ", :949)")
        print("[z_k] Re, Im (ohm): " + ", ".join(f"{z.real:+.5e}{z.imag:+.5e}j" for z in zk))
        print(f"[z_k] off-diagonal residual {t1.off_diagonal_residual:.5e}")
        print(f"[C] C_1 = {t1.c_mode_farad:.6e} F (|rel to :940| {rec_c1:.2e}), "
              f"C_2 = {t2.c_mode_farad:.6e} F (|rel to :940| {rec_c2:.2e}); z_1 rel to :939 "
              f"Re {rec_z1[0]:.2e}, Im {rec_z1[1]:.2e}")
        print(f"[(a')] C_1 vs TH17_STEP3_C_MODE_RECORD_F gap {a_prime:+.6e} (gate :950 "
              f"{REC['a_prime']:+.6e}; rel {rec_ap:.2e}); band {band:.1e}")
        print(f"[S_A] cells {resA.cells}, {wall['SA']:.1f} s; capacitor ports "
              f"{resA.capacitor_port_ids}; zc = {zc.real:.5e}{zc.imag:+.5e}j ohm")
        print(f"[S_A] S11 = {sa11.real:+.6g}{sa11.imag:+.6g}j (rel to :954 {rec_sa:.2e}), "
              f"|S11| = {abs(sa11):.6g}")
        print(f"[(I-A)] |S_A - reduction|/|reduction| {i_a:.6e} ({i_a / band:.4f} x band "
              f"{band:.0e}); gate :954 {REC['i_a']:.6e} (rel {rec_ia:.2e})")
        print(f"[Zin] Z_in(P1)+zc = {zin.real:.5e}{zin.imag:+.5e}j ohm (gate :973 "
              f"{ZIN_PRIOR.real:.5e}{ZIN_PRIOR.imag:+.5e}j; a cancellation of two ~37.5 ohm "
              f"reactances, printed only; near-zero reactance seen from the feed, not a match)")
        print(f"[(c1)] |Im zt_1(C_2)|/|z_1| {c1_val:.5e} ({c1_val / band:.1f} x band; "
              f"gate :970 {REC['c1_factor']:.4e}, rel {rec_c1f:.2e})")
        print(f"[(c2)] S_A vs reduction at C_2 {c2_val:.5e} ({c2_val / band:.1f} x band; "
              f"gate :971 {REC['c2_factor']:.4e}, rel {rec_c2f:.2e})")
        print(f"[(S=)] |bT/aT - S_A|/|S_A| {s_eq:.3e} (band "
              f"{g.LEG_D0_REPRODUCTION_BAND:.0e}; predicted <= 1e-14), "
              f"field solve {wall['field_tuned']:.1f} s")
        print(f"[(P)] P_src {p_src:.6e} W, P_sheet(P1) {p_sheet['P1']:.6e} W, "
              f"sum P_sheet(P2-P4)/P_src {p_rc / p_src:.3e} (predicted <= 1e-9), "
              f"P_vol/P_src {p_vol / p_src:.3e} (predicted order 1e-3); residual "
              f"{p_resid:.3e} (band {DISCRETE_IDENTITY_RTOL:.0e}; predicted ~1e-12)")
        print(f"[(cO)] open-gap S11 = {s11_open.real:+.6g}{s11_open.imag:+.6g}j; "
              f"|S11_open - S_A|/|S_A| {co_val:.5g} ({co_val / band:.3g} x band; predicted "
              f"range [0.3, 2.0], not asserted), solve {wall['field_open']:.1f} s")
        for key in ("tuned", "open"):
            c = centre[key]
            print(f"[B {key}] centre |B| {c['norm']:.5e}, |B1+| {c['b1p']:.5e}, "
                  f"|B1-| {c['b1m']:.5e} T per unit incident wave")
        print(f"[B] |B1+|/|B1-| tuned {centre['tuned']['b1p'] / centre['tuned']['b1m']:.4g} "
              f"(predicted ~1 within a few %); |B_T(0)|/|B_O(0)| "
              f"{centre['tuned']['norm'] / centre['open']['norm']:.4g} (predicted >> 1)")
        print(f"[wall] S0 {wall['S0']:.1f}, S_A {wall['SA']:.1f}, field tuned "
              f"{wall['field_tuned']:.1f}, field open {wall['field_open']:.1f}, "
              f"projections {wall['projections']:.1f}, xdmf {wall['xdmf']:.1f} s; "
              f"summed peak RSS {peak / 2**30:.2f} GiB")
        print(f"[paraview] wrote {xdmf_path}")

    # ---- asserts -------------------------------------------------------------
    assert res0.cells == REC["cells"]
    for c in g.CLASSES:
        assert r0[c] <= g.LEG_D0_REPRODUCTION_BAND, f"(R0) {c} {r0[c]:.3e}"
    assert abs(a_prime) <= band, f"(a') {a_prime:+.6e}"
    assert rec_c1 <= 1e-8 and rec_c2 <= 1e-8, f"C records {rec_c1:.2e} {rec_c2:.2e}"
    assert rec_z1[0] <= 1e-8 and rec_z1[1] <= 1e-8, f"z_1 records {rec_z1}"
    assert rec_ap <= 1e-5, f"(a') record {rec_ap:.2e}"
    assert i_a <= band, f"(I-A) {i_a:.6e}"
    assert rec_ia <= 1e-5, f"(I-A) record {rec_ia:.2e}"
    assert rec_sa <= 1e-8, f"S_A record {rec_sa:.2e}"
    assert rec_c1f <= 1e-3 and rec_c2f <= 1e-3, f"control records {rec_c1f:.2e} {rec_c2f:.2e}"
    assert c1_val > band, f"(c1) blind {c1_val:.4e}"
    assert c2_val > band, f"(c2) blind {c2_val:.4e}"
    assert s_eq <= g.LEG_D0_REPRODUCTION_BAND, f"(S=) {s_eq:.3e}"
    assert p_vol > 0.0
    assert p_resid <= DISCRETE_IDENTITY_RTOL, f"(P) {p_resid:.3e}"
    assert co_val > band, f"(cO) blind {co_val:.4e}"

    _p0(comm, f"\nAll assertions hold. Total elapsed {time.perf_counter() - started:.1f} s.")


if __name__ == "__main__":
    main()
