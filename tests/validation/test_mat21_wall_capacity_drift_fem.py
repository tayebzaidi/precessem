"""`MAT-21` step 40 — P63's bias table re-read against the capacitive drift (§9 item 836,
§10 chain P, P69; brief `docs/planning/chunks/MAT-21.md` "§9 item 836", under the review
ruling 2026-11-05 03:00).

Sets no band and moves none; every band is imported.  Item 682's 1D-reduced walled box
(m = 4, 119 808 cells, Δt = ``T_STAR`` / ``N_MAIN``), CN with two BE start-up half-steps,
``q=Q``, Robin h = ``H`` on z = ±(C + d).  Three runs, P63's fixture body without R3: R0 from
θ = 0, source on from 0, to t\\*; R1 from θ₀ = g·z on every cell (gel and wall), to t\\*; R2 from
θ₀ = g·z with the source off on [0, T_C) and on after, to T_C + t\\*.

Comparand (P67, item 823): item 682's two-layer sourced series ``s_w`` plus the odd
composite-mode drift with the wall's heat capacity, ``_drift_cap_z`` at t\\* (R1) and at
T_C + t\\* (R2), on the unique P2 dof z (P63's ``_ZField``).  The h_eff drift ``_drift_z`` is
kept only as the reproduction of P63's (I).

Asserted (collected, item 668): (I) interior rel L2 ≤ ``INTERIOR_BAND`` on R1 and R2 against
the composite; (E) ≤ ``E_BAND`` on all three, ``source_power`` identical; (S0) on R2; (Sym′) on
R1 − R0 at ``SYM_BAR`` g C; the unique-z precondition ≤ ``PRE_BAR``.  Controls (asserted,
reproduction at the print's last digit): R0 = item 682's 2.244928e-3, and the h_eff
comparand reproduces P63's (I) R1 8.105518e-03 / R2 2.415575e-03, each within 5e-10.
Printed (records, no band): the FEM floor ceiling, and the gap table.

``K_W`` / ``RHO_C_W`` nominal PMMA, citation owed (W3).  Not claimed: the six-wall box, a rate,
a recommended wait.
"""

from __future__ import annotations

import time

import numpy as np
import pytest
from dolfinx import fem
from mpi4py import MPI

from fem_em_solver.thermal import solve_transient_conduction

from tests.validation.test_imp7_duty_cycle_fem import S0_BAND
from tests.validation.test_mat10_sar_to_heat_sphere import _say
from tests.validation.test_mat12_transient_sar_heating_sphere import E_BAND
from tests.validation.test_mat13_transient_perfusion_sphere import INTERIOR_BAND
from tests.validation.test_mat15_coil_driven_heating_phantom import H_GEL, K_GEL, N_MAIN
from tests.validation.test_mat15_finite_cylinder_closed_form import H, K, Q, T_INF
from tests.validation.test_mat15_finite_cylinder_fem import _rel
from tests.validation.test_mat21_box_closed_form import A_BOX, B_BOX, C_BOX, T_C
from tests.validation.test_mat21_initial_gradient_closed_form import G_GRAD, GC
from tests.validation.test_mat21_initial_gradient_fem import SYM_BAR
from tests.validation.test_mat21_wall_capacity_drift_closed_form import _drift_cap, _drift_cap_z
from tests.validation.test_mat21_wall_closed_form import D_WALL, RHO_C_W
from tests.validation.test_mat21_wall_fem import ROBIN, _materials, _robin_facets, _walled_box
from tests.validation.test_mat21_wall_pre_equilibration_fem import (
    DT,
    I_ADJ,
    I_P0,
    I_P5,
    ITEM_682_TSTAR,
    N_OFF,
    N_R2,
    PRE_BAR,
    READ_PTS,
    REPRO_BOUND,
    _census,
    _chi_off_on,
    _drift_z,
    _read,
    _series_z,
    _ZField,
)
from tests.validation.test_mat21_wall_transient_closed_form import T_STAR
from tests.validation.test_mat21_wall_transient_fem import RUNG, _rho_c_field

#: P63's (I) against the h_eff comparand, `20261104T060527Z_MAT-21-P63-W1b.log:66`.
P63_I1, P63_I2 = 8.105518e-03, 2.415575e-03
#: P63's I-mw (R3 vs its exact comparand), same log, printed as the FEM floor only.
P63_IMW = 2.284679e-03
#: P62's FEM − _drift_1d (the wall-capacity term), `…P63-W1b.log:76–81`, printed beside.
P62_TERM = {("R1", "P5"): 2.173640e-02, ("R1", "ADJ_PT"): 1.544980e-02,
            ("R2", "P5"): 1.803301e-03, ("R2", "ADJ_PT"): 1.617858e-03}


def _cap_z(z, t):
    return np.asarray(_drift_cap_z(np.asarray(z, dtype=float), t), dtype=float)


@pytest.fixture(scope="module")
def runs():
    comm = MPI.COMM_WORLD
    t0 = time.perf_counter()
    assert K == K_GEL and H == H_GEL, (K, K_GEL, H, H_GEL)
    msh, ns, nw = _walled_box(RUNG, (False, False, True))
    chi, kf, inter, _ct = _materials(msh)
    rcf = _rho_c_field(msh)
    ft = _robin_facets(msh, (A_BOX, B_BOX, C_BOX + D_WALL), (2,))
    _say(f"\nMAT21-40 walled box m = {RUNG} cells {tuple(ns)} N_w = {nw}: n = "
         f"{msh.topology.index_map(3).size_global}; g {G_GRAD} K/m, g C {GC} K; T_C {T_C:.3f} s, "
         f"t* {T_STAR:.3f} s, dt {DT:.4f} s; R2 {N_R2} steps (off {N_OFF}); ranks {comm.size}")
    grad = lambda x: T_INF + G_GRAD * x[2]  # noqa: E731  (every cell, gel and wall)
    kw = dict(k=kf, q=Q * chi, h=H_GEL, t_inf=T_INF, t_init=T_INF, robin_tag=ROBIN, theta=0.5,
              startup_half_steps=2, dt=DT)
    specs = {"R0": dict(rho_c=rcf, n_steps=N_MAIN),
             "R1": dict(rho_c=rcf, n_steps=N_MAIN, t_init_expr=grad),
             "R2": dict(rho_c=rcf, n_steps=N_R2, t_init_expr=grad, source_schedule=_chi_off_on)}
    out = {"E": {}, "src": {}, "th": {}, "ok": {}, "wall": {}}
    res = {}
    for key, sp in specs.items():
        tq = time.perf_counter()
        res[key] = solve_transient_conduction(msh, ft, **kw, **sp)
        out["wall"][key] = time.perf_counter() - tq
        out["E"][key] = float(res[key].max_balance_residual)
        out["src"][key] = float(res[key].source_power)
        out["th"][key], out["ok"][key] = _read(res[key].temperature, READ_PTS)
        _say(f"MAT21-40 run {key}: {out['wall'][key]:.1f} s; (E) {out['E'][key]:.3e} (E_BAND "
             f"{E_BAND:.0e}); source_power {out['src'][key]:.9e} W")
    out["s0"], tot = _census(res["R2"].source_multipliers)
    _say(f"MAT21-40 (S0) R2 sum chi dt {tot:.9f} s vs t* {T_STAR:.9f} s rel {out['s0']:.3e} "
         f"(band {S0_BAND:g})")

    tq = time.perf_counter()
    zf = _ZField()
    s_w = zf(msh, _series_z(RHO_C_W), T_STAR)
    c_0 = zf(msh, _cap_z, T_STAR)
    c_c = zf(msh, _cap_z, T_C + T_STAR)
    d_0 = zf(msh, _drift_z, T_STAR)
    d_c = zf(msh, _drift_z, T_C + T_STAR)
    out["pre"] = comm.allreduce(zf.worst, op=MPI.MAX)
    _say(f"MAT21-40 five z-only comparands on the unique dof z ({zf.n_unique} values per rank "
         f"max): {time.perf_counter() - tq:.1f} s; precondition worst rel {out['pre']:.3e} (bar "
         f"{PRE_BAR:.0e})")

    th = {key: r.temperature - T_INF for key, r in res.items()}
    out["ref"] = _rel(msh, th["R0"], s_w, inter)
    out["I1"] = _rel(msh, th["R1"], s_w + c_0, inter)
    out["I2"] = _rel(msh, th["R2"], s_w + c_c, inter)
    out["h1"] = _rel(msh, th["R1"], s_w + d_0, inter)
    out["h2"] = _rel(msh, th["R2"], s_w + d_c, inter)
    _say(f"MAT21-40 (I) vs the composite: R1 interior rel L2 {out['I1']:.6e}, R2 {out['I2']:.6e} "
         f"(bound {INTERIOR_BAND:.0e}); whole {_rel(msh, th['R1'], s_w + c_0):.6e} / "
         f"{_rel(msh, th['R2'], s_w + c_c):.6e} (printed)")
    _say(f"MAT21-40 control h_eff reproduction of P63's (I): R1 {out['h1']:.6e} vs {P63_I1:.6e} "
         f"|delta| {abs(out['h1'] - P63_I1):.3e}; R2 {out['h2']:.6e} vs {P63_I2:.6e} |delta| "
         f"{abs(out['h2'] - P63_I2):.3e} (bound {REPRO_BOUND:g})")
    _say(f"MAT21-40 control R0 interior {out['ref']:.6e} vs item 682 {ITEM_682_TSTAR:.6e}, |delta| "
         f"{abs(out['ref'] - ITEM_682_TSTAR):.3e} (bound {REPRO_BOUND:g})")
    _say(f"MAT21-40 ceiling (printed only): FEM floor R0 {out['ref']:.6e} (P63 I-mw {P63_IMW:.6e});"
         f" h_eff/floor R1 {out['h1'] / out['ref']:.3f}x, R2 {out['h2'] / out['ref']:.3f}x; "
         f"composite/floor R1 {out['I1'] / out['ref']:.3f}x, R2 {out['I2'] / out['ref']:.3f}x")
    _say(f"MAT21-40 elapsed {time.perf_counter() - t0:.1f} s")
    return out


def test_all_gates(runs):
    """Collect every gate before asserting (item 668), so one red hides no other."""
    r = runs
    src = r["src"]
    gates = {"valid_points": all(all(o) for o in r["ok"].values()),
             "precondition_unique_z": r["pre"] <= PRE_BAR,
             "I_R1_composite": r["I1"] <= INTERIOR_BAND,
             "I_R2_composite": r["I2"] <= INTERIOR_BAND,
             "E_source_power_identical": src["R0"] > 0.0 and len(set(src.values())) == 1,
             "S0_R2": r["s0"] <= S0_BAND,
             "control_repro_682": abs(r["ref"] - ITEM_682_TSTAR) <= REPRO_BOUND,
             "control_repro_P63_I1": abs(r["h1"] - P63_I1) <= REPRO_BOUND,
             "control_repro_P63_I2": abs(r["h2"] - P63_I2) <= REPRO_BOUND}
    for key, e in r["E"].items():
        gates[f"E_{key}"] = e <= E_BAND
    d = r["th"]["R1"] - r["th"]["R0"]
    s12, s34, s0 = abs(d[1] + d[2]), abs(d[3] + d[4]), abs(d[0])
    _say(f"MAT21-40 (Sym') R1 - R0: |P1 + P2|/(g C) {s12 / GC:.3e}, |P3 + P4|/(g C) "
         f"{s34 / GC:.3e}, |P0|/(g C) {s0 / GC:.3e} (bar {SYM_BAR:.0e}, ASSERTED)")
    gates["Sym_R1"] = max(s12, s34, s0) <= SYM_BAR * GC
    _say(f"MAT21-40 gates: {gates}")
    assert all(gates.values()), (gates, {k: v for k, v in r.items() if k != "th"})


def test_printed_records(runs):
    """The gap table (no band); asserts only that it was formed."""
    th = runs["th"]
    ser = _series_z(RHO_C_W)
    _say("MAT21-40 TABLE run | point | FEM drift (Rk - R0) K | composite K | FEM - composite K | "
         "P62 FEM - _drift_1d K | FEM bias (Rk - R0)/R0 | composite / _series")
    for key, t_w, tag in (("R1", 0.0, "t_w = 0"), ("R2", T_C, "t_w = T_C")):
        for i, lab in ((I_P5, "P5"), (I_ADJ, "ADJ_PT"), (I_P0, "P0")):
            p = tuple(READ_PTS[i])
            dh = th[key][i] - th["R0"][i]
            cap = _drift_cap(p, T_STAR + t_w)
            sv = float(np.ravel(ser(np.array([p[2]]), T_STAR))[0])
            term = P62_TERM.get((key, lab), float("nan"))
            cb = cap / sv if sv != 0.0 else float("nan")
            _say(f"MAT21-40 TABLE {key} {tag} | {lab} | {dh:+.6e} | {cap:+.6e} | {dh - cap:+.3e} | "
                 f"{term:+.6e} | {dh / th['R0'][i]:+.6e} | {cb:+.6e}")
            assert np.isfinite(dh / th["R0"][i])
    _say("MAT21-40 wall time per run: " + ", ".join(f"{k} {v:.1f} s"
                                                     for k, v in runs["wall"].items()))
