"""`MAT-21` step 36 — the pre-equilibration FEM twin inside the walled box (§9 item 805,
§10 chain P, P63; brief `docs/planning/chunks/MAT-21.md` "§9 item 805").

Sets no band and moves none; every band is imported.  Item 682's 1D-reduced walled box
(`test_mat21_wall_transient_fem`: ``_walled_box`` / ``_materials`` / ``_rho_c_field``, m = 4,
119 808 cells, Δt = ``T_STAR`` / ``N_MAIN``), CN with two BE start-up half-steps, ``q=Q``,
Robin h = ``H`` on z = ±(C + d).  Four runs: R0 from θ = 0, source on from 0, to t\\*; R1 from
θ₀ = g·z on every cell (gel and wall), to t\\*; R2 from θ₀ = g·z with the source off on
[0, T_C) and on after, to T_C + t\\*; R3 = R1 with ρ_wc_w = ``MASSLESS`` · ``RHO_C_W``.

Comparand (P62, item 788): the z-only h_eff drift ``_drift_1d`` plus item 682's two-layer
series ``_series``; R3's series at the massless wall, where the meshed wall *is* the h_eff
resistance, so R3's comparand is exact up to discretisation.  Every comparand depends on z
only, so it is evaluated on the unique P2 dof z (rounded to 1e-12 m for grouping; each group
is evaluated at its first member's own z) and mapped back; precondition: on the first
``INTERP_CHUNK`` dofs of every rank the mapped values equal the direct evaluation to 1e-13
relative (max |diff| / max |direct|).

Asserted (collected, item 668): (I) interior rel L2 ≤ ``INTERIOR_BAND`` on R1 and R2; (I-mw)
the same on R3; (E) ≤ ``E_BAND`` on all four and ``source_power`` identical on all four; (S0)
R2's Σχ·Δt = t\\* to ``S0_BAND``; (Sym′) the point-reflection pair identity on R1 − R0 at
``SYM_BAR`` g C (W0 first checks the P2 dof set is closed under p → −p to 1e-12 m); R2 − R0's
pair sums are printed only (R2's CN switch-on vs R0's BE start-up leaves an even residue).
Controls (asserted): R0 reproduces item 682's t\\* interior record 2.244928e-3 within 5e-10;
the wall-free series (P51's comparand at H: ``_grad_ex`` + ``_trans_ex``) misses R1 at t\\* by
> ``INTERIOR_BAND``.  Printed (records, no band): the FEM drift part beside ``_drift_1d``
(the wall-capacity term), R1 − R3, the bias table beside item 788's six-face form.

``K_W`` / ``RHO_C_W`` are nominal PMMA values, citation owed (W3).  Not claimed: the six-wall
box, a rate, a recommended wait.
"""

from __future__ import annotations

import time

import numpy as np
import pytest
from dolfinx import fem
from dolfinx import mesh as dmesh
from mpi4py import MPI

from fem_em_solver.post.evaluation import evaluate_vector_field_parallel
from fem_em_solver.thermal import solve_transient_conduction

from tests.validation.test_imp7_duty_cycle_fem import S0_BAND
from tests.validation.test_mat10_sar_to_heat_sphere import _say
from tests.validation.test_mat12_transient_sar_heating_sphere import E_BAND
from tests.validation.test_mat13_transient_perfusion_sphere import INTERIOR_BAND
from tests.validation.test_mat15_coil_driven_heating_phantom import H_GEL, K_GEL, N_MAIN
from tests.validation.test_mat15_finite_cylinder_closed_form import H, K, Q, T_INF
from tests.validation.test_mat15_finite_cylinder_fem import _rel
from tests.validation.test_mat21_box_closed_form import A_BOX, B_BOX, C_BOX, T_C
from tests.validation.test_mat21_box_fem import _exact, _trans_ex
from tests.validation.test_mat21_box_probe_points_fem import PROBE_LABELS, PROBE_PTS
from tests.validation.test_mat21_initial_gradient_closed_form import G_GRAD, GC
from tests.validation.test_mat21_initial_gradient_fem import SYM_BAR, _grad_ex
from tests.validation.test_mat21_pre_equilibration_closed_form import ADJ_PT, P5
from tests.validation.test_mat21_wall_closed_form import D_WALL, RHO_C_W
from tests.validation.test_mat21_wall_fem import ROBIN, _materials, _robin_facets, _walled_box
from tests.validation.test_mat21_wall_pre_equilibration_closed_form import (
    _drift,
    _drift_1d,
    _sourced,
)
from tests.validation.test_mat21_wall_transient_closed_form import T_STAR
from tests.validation.test_mat21_wall_transient_fem import (
    INTERP_CHUNK,
    MASSLESS,
    RUNG,
    _lam,
    _rho_c_field,
    _series,
)

DT = T_STAR / N_MAIN
N_OFF = int(round(T_C / DT))
N_R2 = N_OFF + N_MAIN
#: Item 682's t* interior record, `20261028T021107Z_MAT-21-step20c-W1-fem.log:54`.
ITEM_682_TSTAR = 2.244928e-3
REPRO_BOUND = 5e-10  # the print's last digit, the 18:00 ruling on I41
PRE_BAR = 1e-13  # the brief's unique-z precondition, relative
SYM_DOF_BAR = 1e-12  # m, the brief's dof-set closure check
Z_ROUND = 12
READ_LABELS = tuple(PROBE_LABELS) + ("P5 (gate)", "ADJ_PT")
READ_PTS = np.vstack([np.array(PROBE_PTS, dtype=float), np.array(P5, dtype=float),
                      np.array(ADJ_PT, dtype=float)])
I_P5, I_ADJ, I_P0 = 6, 7, 0


def _chi_off_on(t):
    """0 on [0, T_C), 1 after (one switch, on the grid); t is always a (half-)step midpoint."""
    return 1.0 if int(np.floor(t / DT)) >= N_OFF else 0.0


def _mids(n_steps):
    """The schedule's evaluation times (`transient.py:352–357`)."""
    return [(i + 0.5) * (0.5 * DT) for i in range(2)] + [(1 + i + 0.5) * DT
                                                          for i in range(n_steps - 1)]


def _census(mult, startup=2):
    m = np.asarray(mult, dtype=float)
    dts = np.array([0.5 * DT] * startup + [DT] * (len(m) - startup))
    tot = float(np.dot(m, dts))
    return abs(tot - T_STAR) / T_STAR, tot


def _drift_z(z, t):
    return np.array([_drift_1d((0.0, 0.0, float(zi)), t) for zi in np.atleast_1d(z)])


def _series_z(rho_c_w):
    lam = _lam(rho_c_w)
    return lambda z, t: np.asarray(_series(np.asarray(z, dtype=float), t, rho_c_w, lam))


class _ZField:
    """A z-only comparand on the unique dof z, mapped back; records the precondition."""

    def __init__(self):
        self.worst = 0.0
        self.n_unique = 0

    def __call__(self, msh, fz, t):
        g = fem.Function(fem.functionspace(msh, ("Lagrange", 2)))

        def ev(x):
            z = np.ascontiguousarray(x[2])
            _u, idx, inv = np.unique(np.round(z, Z_ROUND), return_index=True,
                                     return_inverse=True)
            out = np.asarray(fz(z[idx], t), dtype=float)[np.ravel(inv)]
            k = min(INTERP_CHUNK, z.size)
            if k:
                direct = np.asarray(fz(z[:k], t), dtype=float)
                scale = max(float(np.max(np.abs(direct))), 1e-300)
                self.worst = max(self.worst, float(np.max(np.abs(out[:k] - direct))) / scale)
            self.n_unique = max(self.n_unique, int(idx.size))
            return out

        g.interpolate(ev)
        return g


def _dof_closure(msh):
    """max over P2 dofs of the distance from −p to the nearest dof (gathered, -n 1 or more)."""
    from scipy.spatial import cKDTree

    V2 = fem.functionspace(msh, ("Lagrange", 2))
    nloc = V2.dofmap.index_map.size_local
    xs = V2.tabulate_dof_coordinates()[:nloc]
    allx = msh.comm.gather(xs, root=0)
    d = None
    if msh.comm.rank == 0:
        X = np.vstack(allx)
        d = float(np.max(cKDTree(X).query(-X)[0]))
    return msh.comm.bcast(d, root=0), msh.comm.allreduce(nloc, op=MPI.SUM)


def _read(T, pts):
    vals, ok = evaluate_vector_field_parallel(T, pts)
    th = np.array([float(np.real(np.ravel(vals[i])[0])) - T_INF if bool(ok[i]) else np.nan
                   for i in range(pts.shape[0])])
    return th, [bool(o) for o in ok]


def test_w0_premises():
    """W0 (no solve): the schedule arithmetic, the dof-set closure under p → −p on item 682's
    mesh, and the unique-z precondition on a small create_box."""
    t0 = time.perf_counter()
    assert K == K_GEL and H == H_GEL, (K, K_GEL, H, H_GEL)
    assert abs(T_C / DT - N_OFF) <= 1e-12 * N_OFF, (T_C / DT, N_OFF)
    mult = [_chi_off_on(t) for t in _mids(N_R2)]
    rel, tot = _census(mult)
    _say(f"\nMAT21-36 W0 T_C/dt {T_C / DT:.15f} (N_OFF {N_OFF}); R2 {N_R2} CN steps; sum chi dt "
         f"{tot:.6f} s vs t* {T_STAR:.6f} s rel {rel:.3e} (band {S0_BAND:g}); on "
         f"{int(sum(mult))} of {len(mult)}")
    assert rel <= S0_BAND, rel
    msh, ns, nw = _walled_box(RUNG, (False, False, True))
    dmax, ndof = _dof_closure(msh)
    _say(f"MAT21-36 W0 walled box m = {RUNG} cells {tuple(ns)} N_w = {nw}: P2 dofs {ndof}; dof "
         f"set closure under p -> -p: max dist {dmax:.3e} m (bar {SYM_DOF_BAR:.0e})")
    small = dmesh.create_box(MPI.COMM_WORLD, [np.array([-A_BOX, -B_BOX, -C_BOX - D_WALL]),
                                              np.array([A_BOX, B_BOX, C_BOX + D_WALL])],
                             [4, 3, 8], cell_type=dmesh.CellType.tetrahedron)
    zf = _ZField()
    zf(small, _series_z(RHO_C_W), T_STAR)
    zf(small, _drift_z, T_STAR)
    zf(small, _drift_z, T_C + T_STAR)
    zf(small, _series_z(MASSLESS * RHO_C_W), T_STAR)
    _say(f"MAT21-36 W0 unique-z precondition on a small create_box: worst rel {zf.worst:.3e} "
         f"(bar {PRE_BAR:.0e}); unique z {zf.n_unique}; elapsed {time.perf_counter() - t0:.1f} s")
    assert dmax <= SYM_DOF_BAR, dmax
    assert zf.worst <= PRE_BAR, zf.worst


@pytest.fixture(scope="module")
def runs():
    comm = MPI.COMM_WORLD
    t0 = time.perf_counter()
    assert K == K_GEL and H == H_GEL, (K, K_GEL, H, H_GEL)
    msh, ns, nw = _walled_box(RUNG, (False, False, True))
    chi, kf, inter, _ct = _materials(msh)
    rcf = _rho_c_field(msh)
    rc0 = fem.Function(rcf.function_space)
    rc0.x.array[:] = np.where(rcf.x.array == RHO_C_W, MASSLESS * RHO_C_W, rcf.x.array)
    ft = _robin_facets(msh, (A_BOX, B_BOX, C_BOX + D_WALL), (2,))
    _say(f"\nMAT21-36 walled box m = {RUNG} cells {tuple(ns)} N_w = {nw}: n = "
         f"{msh.topology.index_map(3).size_global}; g {G_GRAD} K/m, g C {GC} K; T_C {T_C:.3f} s, "
         f"t* {T_STAR:.3f} s, dt {DT:.4f} s; R2 {N_R2} steps (off {N_OFF}); ranks {comm.size}")
    grad = lambda x: T_INF + G_GRAD * x[2]  # noqa: E731  (every cell, gel and wall)
    kw = dict(k=kf, q=Q * chi, h=H_GEL, t_inf=T_INF, t_init=T_INF, robin_tag=ROBIN, theta=0.5,
              startup_half_steps=2, dt=DT)
    specs = {"R0": dict(rho_c=rcf, n_steps=N_MAIN),
             "R1": dict(rho_c=rcf, n_steps=N_MAIN, t_init_expr=grad),
             "R2": dict(rho_c=rcf, n_steps=N_R2, t_init_expr=grad, source_schedule=_chi_off_on),
             "R3": dict(rho_c=rc0, n_steps=N_MAIN, t_init_expr=grad)}
    out = {"E": {}, "src": {}, "th": {}, "ok": {}, "wall": {}}
    res = {}
    for key, sp in specs.items():
        tq = time.perf_counter()
        res[key] = solve_transient_conduction(msh, ft, **kw, **sp)
        out["wall"][key] = time.perf_counter() - tq
        out["E"][key] = float(res[key].max_balance_residual)
        out["src"][key] = float(res[key].source_power)
        out["th"][key], out["ok"][key] = _read(res[key].temperature, READ_PTS)
        _say(f"MAT21-36 run {key}: {out['wall'][key]:.1f} s; (E) {out['E'][key]:.3e} (E_BAND "
             f"{E_BAND:.0e}); source_power {out['src'][key]:.9e} W")
    out["s0"], tot = _census(res["R2"].source_multipliers)
    _say(f"MAT21-36 (S0) R2 sum chi dt {tot:.9f} s vs t* {T_STAR:.9f} s rel {out['s0']:.3e} "
         f"(band {S0_BAND:g}); on {int(np.sum(np.asarray(res['R2'].source_multipliers)))} of "
         f"{len(res['R2'].source_multipliers)}")

    tq = time.perf_counter()
    zf = _ZField()
    s_w = zf(msh, _series_z(RHO_C_W), T_STAR)
    s_0 = zf(msh, _series_z(MASSLESS * RHO_C_W), T_STAR)
    d_0 = zf(msh, _drift_z, T_STAR)
    d_c = zf(msh, _drift_z, T_C + T_STAR)
    out["pre"] = comm.allreduce(zf.worst, op=MPI.MAX)
    _say(f"MAT21-36 four z-only comparands on the unique dof z ({zf.n_unique} values per rank "
         f"max): {time.perf_counter() - tq:.1f} s; precondition worst rel {out['pre']:.3e} (bar "
         f"{PRE_BAR:.0e})")
    tq = time.perf_counter()
    wf_d, npts = _exact(msh, _grad_ex(T_STAR))
    wf_s, _ = _exact(msh, _trans_ex(T_STAR))
    _say(f"MAT21-36 wall-free control, two 3D series interpolations ({npts} P2 points): "
         f"{time.perf_counter() - tq:.1f} s")

    th = {key: r.temperature - T_INF for key, r in res.items()}
    out["ref"] = _rel(msh, th["R0"], s_w, inter)
    out["I1"] = _rel(msh, th["R1"], s_w + d_0, inter)
    out["I2"] = _rel(msh, th["R2"], s_w + d_c, inter)
    out["Imw"] = _rel(msh, th["R3"], s_0 + d_0, inter)
    out["wf"] = _rel(msh, th["R1"], wf_d + wf_s, inter)
    _say(f"MAT21-36 (I) R1 interior rel L2 {out['I1']:.6e}, R2 {out['I2']:.6e} (bound "
         f"{INTERIOR_BAND:.0e}); whole {_rel(msh, th['R1'], s_w + d_0):.6e} / "
         f"{_rel(msh, th['R2'], s_w + d_c):.6e} (printed)")
    _say(f"MAT21-36 (I-mw) R3 interior rel L2 vs its exact comparand {out['Imw']:.6e} (bound "
         f"{INTERIOR_BAND:.0e}); whole {_rel(msh, th['R3'], s_0 + d_0):.6e} (printed)")
    _say(f"MAT21-36 control R0 interior {out['ref']:.6e} vs item 682 {ITEM_682_TSTAR:.6e}, |delta| "
         f"{abs(out['ref'] - ITEM_682_TSTAR):.3e} (bound {REPRO_BOUND:g}); wall-free series miss "
         f"of R1 {out['wf']:.4e} (asserted > {INTERIOR_BAND:.0e})")
    _say(f"MAT21-36 elapsed {time.perf_counter() - t0:.1f} s")
    return out


def test_all_gates(runs):
    """Collect every gate before asserting (item 668), so one red hides no other."""
    r = runs
    src = r["src"]
    gates = {"valid_points": all(all(o) for o in r["ok"].values()),
             "precondition_unique_z": r["pre"] <= PRE_BAR,
             "I_R1": r["I1"] <= INTERIOR_BAND,
             "I_R2": r["I2"] <= INTERIOR_BAND,
             "I_mw_R3": r["Imw"] <= INTERIOR_BAND,
             "E_source_power_identical": src["R0"] > 0.0 and len(set(src.values())) == 1,
             "S0_R2": r["s0"] <= S0_BAND,
             "control_repro_682": abs(r["ref"] - ITEM_682_TSTAR) <= REPRO_BOUND,
             "control_wall_free_misses": r["wf"] > INTERIOR_BAND}
    for key, e in r["E"].items():
        gates[f"E_{key}"] = e <= E_BAND
    ref = r["th"]["R0"]
    # (Sym') is the brief's on R1 - R0 only.  R2 - R0 is printed, not gated: R2's source
    # switches on at step N_OFF under CN while R0's starts under the two BE half-steps, so
    # R2 - R0 carries an even (z-symmetric) discretisation residue besides the odd drift.
    # Measured in the first W1 (`20261104T060402Z_MAT-21-P63-W1.log:71`, `:81`): |P0|/(g C)
    # 3.951e-05 = the P0 reading 1.778e-05 K, |P1 + P2| = 2 |P0|; a record, not a symmetry break.
    for key in ("R1", "R2"):
        d = r["th"][key] - ref
        s12, s34, s0 = abs(d[1] + d[2]), abs(d[3] + d[4]), abs(d[0])
        tag = f"bar {SYM_BAR:.0e}, ASSERTED" if key == "R1" else "printed, no band"
        _say(f"MAT21-36 (Sym') {key} - R0: |P1 + P2|/(g C) {s12 / GC:.3e}, |P3 + P4|/(g C) "
             f"{s34 / GC:.3e}, |P0|/(g C) {s0 / GC:.3e} ({tag})")
        if key == "R1":
            gates["Sym_R1"] = max(s12, s34, s0) <= SYM_BAR * GC
    _say(f"MAT21-36 gates: {gates}")
    assert all(gates.values()), (gates, {k: v for k, v in r.items() if k != "th"})


def test_printed_records(runs):
    """The records (no band); asserts only that they were formed."""
    r = runs
    th = r["th"]
    _say("MAT21-36 TABLE run | point | FEM drift (Rk - R0) K | _drift_1d K | FEM - _drift_1d K "
         "(wall-capacity term) | R1 - R3 K | FEM bias | item 788 bias (six-face) | 1D bias")
    for key, t_w, tag in (("R1", 0.0, "t_w = 0"), ("R2", T_C, "t_w = T_C")):
        for i, lab in ((I_P5, "P5"), (I_ADJ, "ADJ_PT"), (I_P0, "P0")):
            p = tuple(READ_PTS[i])
            dh = th[key][i] - th["R0"][i]
            d1 = _drift_1d(p, T_STAR + t_w)
            r13 = th["R1"][i] - th["R3"][i] if key == "R1" else float("nan")
            s788 = _sourced(p, T_STAR)
            b788 = _drift(p, T_STAR + t_w) / s788
            _say(f"MAT21-36 TABLE {key} {tag} | {lab} | {dh:+.6e} | {d1:+.6e} | {dh - d1:+.3e} | "
                 f"{r13:+.3e} | {dh / th['R0'][i]:+.6e} | {b788:+.6e} | {d1 / s788:+.6e}")
            assert np.isfinite(dh / th["R0"][i])
    for i, lab in ((I_P5, "P5"), (I_ADJ, "ADJ_PT")):
        p = tuple(READ_PTS[i])
        worst = max(abs(_drift(p, T_STAR + w) - _drift_1d(p, T_STAR + w)) / GC for w in (0.0, T_C))
        _say(f"MAT21-36 SIDE {lab}: |six-face _drift - _drift_1d|/(g C) max over t_w in "
             f"{{0, T_C}} {worst:.3e} (printed, item 788 measured <= 2.295e-06)")
    _say("MAT21-36 wall time per run: " + ", ".join(f"{k} {v:.1f} s" for k, v in r["wall"].items()))
