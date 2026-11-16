"""`MAT-21` step 27 — the walled box's Crank–Nicolson under D = 0.25 against the two-layer
slab series (§9 item 721, §10 chain P, P54; brief `docs/planning/chunks/MAT-21.md` "§9 item
721").

Sets no band and moves none; every band is imported.  Item 688's fixture by import, on one
mesh, two runs: ``_chi(PERIOD, 0.25)`` (an on-phase of 5 steps) and the D = 0.5 control.
The comparand is `two_layer_slab_duty_cycle_temperature` at each run's duty, with item 688's
side convention, in chunks into P2.  Read in one fixed batch: P37's six points plus item
713's third height, (0, 0, 0.9 C).

Asserted (collected, item 668): interior rel L2 ≤ ``INTERIOR_BAND`` at t\\* and t₂; (E) ≤
``E_BAND``; (S0) ≤ ``S0_BAND``, each run against its own D·t₂.  Controls: the D = 0.5 run
reproduces item 688's record to the print's last digit, 5e-10 absolute; the continuous and
the D = 0.5 series each miss the D = 0.25 run by > ``CTL_FACTOR`` × its interior error.
Printed (records, **no band**): the FEM θ beside item 713's series; the z = 0 five-point
spread, the lateral discretisation floor of the 1D-reduced geometry.

Not claimed: the six-wall box, a dose, lateral structure.
"""

from __future__ import annotations

import time

import numpy as np
import pytest
from dolfinx import fem

from fem_em_solver.post.evaluation import evaluate_vector_field_parallel
from fem_em_solver.thermal import solve_transient_conduction
from fem_em_solver.utils.analytical import two_layer_slab_duty_cycle_temperature

from tests.validation.test_imp7_duty_cycle_fem import S0_BAND, _chi
from tests.validation.test_mat10_sar_to_heat_sphere import _say
from tests.validation.test_mat12_transient_sar_heating_sphere import E_BAND
from tests.validation.test_mat13_transient_perfusion_sphere import INTERIOR_BAND
from tests.validation.test_mat15_coil_driven_heating_phantom import H_GEL, K_GEL, N_MAIN
from tests.validation.test_mat15_finite_cylinder_closed_form import H, K, Q, T_INF
from tests.validation.test_mat15_finite_cylinder_fem import _rel
from tests.validation.test_mat21_box_closed_form import A_BOX, B_BOX, C_BOX
from tests.validation.test_mat21_box_duty25_dose_closed_form import DUTY25
from tests.validation.test_mat21_box_probe_points_fem import PROBE_LABELS, PROBE_PTS
from tests.validation.test_mat21_wall_closed_form import D_WALL, RHO_C_W
from tests.validation.test_mat21_wall_duty_closed_form import PERIOD, _duty
from tests.validation.test_mat21_wall_duty_fem import DT
from tests.validation.test_mat21_wall_fem import ROBIN, _materials, _robin_facets, _walled_box
from tests.validation.test_mat21_wall_transient_closed_form import N_SERIES, T_STAR
from tests.validation.test_mat21_wall_transient_fem import (
    CTL_FACTOR,
    INTERP_CHUNK,
    N2,
    RUNG,
    SLAB,
    T2,
    _exact_t,
    _lam,
    _rho_c_field,
)

#: Item 688's (t*, t₂), `20261028T172612Z_MAT-21-step21b-W1.log:55`, `:57`.
ITEM_688 = {N_MAIN: 1.697155e-3, N2: 1.447368e-3}
REPRO_BOUND = 5e-10  # the print's last digit, the 18:00 ruling on I41
READ_LABELS = tuple(PROBE_LABELS) + ("(0, 0, 0.9 C)",)
READ_PTS = np.vstack([np.array(PROBE_PTS, dtype=float), [0.0, 0.0, 0.9 * C_BOX]])


def _exact_duty_at(msh, t, duty, lam):
    """Item 688's `_exact_duty_t`, whose series is bound to ``DUTY``, at any duty."""
    g = fem.Function(fem.functionspace(msh, ("Lagrange", 2)))

    def ev(x):
        n = x.shape[1]
        out = np.empty(n)
        for s in range(0, n, INTERP_CHUNK):
            out[s:s + INTERP_CHUNK] = two_layer_slab_duty_cycle_temperature(
                x[2, s:s + INTERP_CHUNK], t, q=Q, h=H, t_inf=0.0, rho_c_w=RHO_C_W, period=PERIOD,
                duty=duty, n_terms=N_SERIES, eigenvalues=lam, **SLAB)
        return out

    g.interpolate(ev)
    return g


def _census_at(mult, duty, startup=2):
    m = np.asarray(mult, dtype=float)
    dts = np.array([0.5 * DT] * startup + [DT] * (len(m) - startup))
    tot, target = float(np.dot(m, dts)), duty * T2
    return abs(tot - target) / target, tot, target


def _mids():
    """The schedule's evaluation times (`transient.py:352–357`)."""
    return [(i + 0.5) * (0.5 * DT) for i in range(2)] + [(1 + i + 0.5) * DT for i in range(N2 - 1)]


@pytest.fixture(scope="module")
def walled():
    t0 = time.perf_counter()
    assert K == K_GEL and H == H_GEL, (K, K_GEL, H, H_GEL)
    msh, ns, nw = _walled_box(RUNG, (False, False, True))
    chi, kf, inter, _ct = _materials(msh)
    rcf = _rho_c_field(msh)
    ft = _robin_facets(msh, (A_BOX, B_BOX, C_BOX + D_WALL), (2,))
    lam_w = _lam(RHO_C_W)
    _say(f"\nMAT21-27 1D-reduced walled box m = {RUNG} cells {tuple(ns)} N_w = {nw}: n = "
         f"{msh.topology.index_map(3).size_global}; t* {T_STAR:.3f} s, t2 {T2:.3f} s; P "
         f"{PERIOD:.4f} s = {PERIOD / DT:.1f} steps; dt {DT:.4f} s, {N2} CN steps")
    out = {"E": {}, "s0": {}, "e": {}, "read": {}, "ok": {}}
    res = {}
    for d in (DUTY25, 0.5):
        tq = time.perf_counter()
        r = solve_transient_conduction(
            msh, ft, k=kf, rho_c=rcf, q=Q * chi, h=H_GEL, t_inf=T_INF, t_init=T_INF,
            robin_tag=ROBIN, theta=0.5, startup_half_steps=2, dt=DT, n_steps=N2,
            snapshot_steps=(N_MAIN, N2), source_schedule=_chi(PERIOD, d))
        res[d] = r
        out["E"][d] = float(r.max_balance_residual)
        s0, tot, target = _census_at(r.source_multipliers, d)
        out["s0"][d] = s0
        _say(f"MAT21-27 D={d} CN solve {time.perf_counter() - tq:.1f} s; (E) {out['E'][d]:.3e} "
             f"(E_BAND {E_BAND:.0e}); (S0) sum chi dt {tot:.6f} s vs D t2 {target:.6f} s rel "
             f"{s0:.3e} (band {S0_BAND:g}); on {int(np.sum(np.asarray(r.source_multipliers)))}")

    for step, t in ((N_MAIN, T_STAR), (N2, T2)):
        tq = time.perf_counter()
        ex25 = _exact_duty_at(msh, t, DUTY25, lam_w)
        ex50 = _exact_duty_at(msh, t, 0.5, lam_w)
        exc = _exact_t(msh, t, RHO_C_W, lam_w)
        th25 = res[DUTY25].snapshots[step] - T_INF
        th50 = res[0.5].snapshots[step] - T_INF
        rec = {"e25": _rel(msh, th25, ex25, inter), "e50": _rel(msh, th50, ex50, inter),
               "miss_c": _rel(msh, th25, exc, inter), "miss_50": _rel(msh, th25, ex50, inter)}
        out["e"][step] = rec
        _say(f"MAT21-27 step {step} t = {t:.3f} s: D=0.25 interior rel L2 {rec['e25']:.6e} (bound "
             f"{INTERIOR_BAND:.0e}), whole {_rel(msh, th25, ex25):.6e} (printed); D=0.5 interior "
             f"{rec['e50']:.6e} vs item 688 {ITEM_688[step]:.6e}, |delta| "
             f"{abs(rec['e50'] - ITEM_688[step]):.3e} (bound {REPRO_BOUND:g}); misses of the "
             f"D=0.25 run: continuous {rec['miss_c']:.6e}, D=0.5 series {rec['miss_50']:.6e}; "
             f"{time.perf_counter() - tq:.1f} s")
        vals, ok = evaluate_vector_field_parallel(res[DUTY25].snapshots[step], READ_PTS)
        out["ok"][step] = [bool(o) for o in ok]
        out["read"][step] = np.array(
            [float(np.real(np.ravel(vals[i])[0])) - T_INF if bool(ok[i]) else np.nan
             for i in range(READ_PTS.shape[0])])
    _say(f"MAT21-27 elapsed {time.perf_counter() - t0:.1f} s")
    return out


def test_w0_schedule_arithmetic():
    """W0 (numpy): the grid and each duty's census on its own closure, with no solve."""
    assert abs(PERIOD / DT - 20.0) < 1e-12 and abs(DUTY25 * PERIOD / DT - 5.0) < 1e-12
    assert len(READ_LABELS) == READ_PTS.shape[0] == 7
    for d in (DUTY25, 0.5):
        chi = _chi(PERIOD, d)
        mult = [chi(t) for t in _mids()]
        rel, tot, target = _census_at(mult, d)
        _say(f"MAT21-27 W0 D={d}: {len(mult)} (half-)steps; sum chi dt {tot:.6f} s vs "
             f"{target:.6f} s rel {rel:.3e} (band {S0_BAND:g})")
        assert rel <= S0_BAND, (d, rel)


def test_all_gates(walled):
    """Collect every gate before asserting (item 668), so one red hides no other."""
    w = walled
    gates = {"valid_points": all(all(o) for o in w["ok"].values())}
    for d in (DUTY25, 0.5):
        gates[f"balance_D{d}"] = w["E"][d] <= E_BAND
        gates[f"S0_D{d}"] = w["s0"][d] <= S0_BAND
    for step, name in ((N_MAIN, "tstar"), (N2, "t2")):
        r = w["e"][step]
        gates[f"interior_{name}"] = r["e25"] <= INTERIOR_BAND
        gates[f"repro_688_{name}"] = abs(r["e50"] - ITEM_688[step]) <= REPRO_BOUND
        gates[f"control_continuous_{name}"] = r["miss_c"] > CTL_FACTOR * r["e25"]
        gates[f"control_D50_{name}"] = r["miss_50"] > CTL_FACTOR * r["e25"]
        _say(f"MAT21-27 (control) {name}: continuous {r['miss_c'] / r['e25']:.2f}x, D=0.5 series "
             f"{r['miss_50'] / r['e25']:.2f}x the FEM error (asserted > {CTL_FACTOR:.0f}x)")
    _say(f"MAT21-27 gates: {gates}")
    assert all(gates.values()), (gates, w)


def test_printed_six_points(walled):
    """The records (no band); asserts only that they were formed."""
    _say("MAT21-27 TABLE t | point | FEM theta K | item 713 series K | FEM - series K")
    for step, name, t in ((N_MAIN, "t*", T_STAR), (N2, "t2", T2)):
        th = walled["read"][step]
        for i, lab in enumerate(READ_LABELS):
            s = float(_duty(READ_PTS[i, 2], t, duty=DUTY25, side="left")[0])
            _say(f"MAT21-27 TABLE {name} | {lab} | {th[i]:.6f} | {s:.6f} | {th[i] - s:+.3e}")
        spread = float(np.max(th[:5]) - np.min(th[:5]))
        _say(f"MAT21-27 z = 0 five-point spread at {name}: {spread:.3e} K (record, the lateral "
             f"discretisation floor)")
        assert np.all(np.isfinite(th))
