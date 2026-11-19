"""`MAT-21` step 26 — the pre-equilibration FEM twin of P51 on step 1's box (§9 item 720,
§10 chain P, P52; brief `docs/planning/chunks/MAT-21.md` "§9 item 720").

Sets no band and moves none; every band is imported.  P47b′'s build (`gradient_run`,
``RUNGS[-1]``, ``create_box``) with the source on, P1 CN to t\\*, ``N_MAIN`` steps, two BE
half-steps, ``q=Q`` from t = 0, three runs: the reference from uniform T∞; t_w = 0 from
T∞ + g·z; t_w = ``T_C`` from T∞ + the drift series at t = ``T_C`` on the mesh.  Real build.

Asserted: (I) interior rel L2 ≤ ``INTERIOR_BAND`` against P51's superposition
θ_drift(t\\* + t_w) + θ_sourced(t\\*), on both gradient runs; (E) ≤ ``E_BAND`` on all three,
the source branch of the balance (``source_power`` > 0, equal on all three); (Sym′) the
point-reflection pair identity on the drift part, run minus reference, because the sourced
part is even.  Controls: the reference's interior error reproduces step 1's record, to the
print's last digit; the Dirichlet-z series misses at t\\* by > ``INTERIOR_BAND``.  Printed
(records, **no band**): the bias at P37's six points beside P51's; the five z = 0 rows are
the FEM's symmetry floor.

Not claimed: a recommended wait, a rate, the walled box.
"""

from __future__ import annotations

import time

import numpy as np
import pytest
from dolfinx import mesh as dmesh
from mpi4py import MPI

from fem_em_solver.post.evaluation import evaluate_vector_field_parallel
from fem_em_solver.thermal import solve_transient_conduction
from tests.validation.test_mat10_sar_to_heat_sphere import ROBIN_TAG, _robin_tags, _say
from tests.validation.test_mat12_transient_sar_heating_sphere import E_BAND, RHO_C
from tests.validation.test_mat13_transient_perfusion_sphere import INTERIOR_BAND
from tests.validation.test_mat15_coil_driven_heating_phantom import H_GEL, K_GEL, N_MAIN, T_INF
from tests.validation.test_mat15_finite_cylinder_fem import _rel
from tests.validation.test_mat21_box_closed_form import A_BOX, B_BOX, C_BOX, T_C, Q
from tests.validation.test_mat21_box_fem import (
    BASE,
    INTERP_CHUNK,
    RUNGS,
    T_STAR,
    _exact,
    _interior_chi,
    _trans_ex,
)
from tests.validation.test_mat21_box_probe_points_fem import PROBE_LABELS, PROBE_PTS
from tests.validation.test_mat21_initial_gradient_closed_form import G_GRAD, GC
from tests.validation.test_mat21_initial_gradient_fem import DZ_ROOTS, SYM_BAR, _grad_ex
from tests.validation.test_mat21_pre_equilibration_closed_form import _drift, _sourced

#: Step 1's (T), `20261028T003134Z_MAT-21-step19-Wreg.log:92`.
STEP1_T = 9.071528e-04
STEP1_T_BOUND = 1e-10  # the print's last digit (`.6e`), the 18:00 ruling on I41


def _wait_init(t_w):
    """T∞ + the drift series at t_w, on any point set, chunked as `_exact`."""
    fn = _grad_ex(t_w)

    def ev(x):
        n = x.shape[1]
        out = np.empty(n)
        for s in range(0, n, INTERP_CHUNK):
            out[s : s + INTERP_CHUNK] = np.real(
                fn(
                    x[0, s : s + INTERP_CHUNK],
                    x[1, s : s + INTERP_CHUNK],
                    x[2, s : s + INTERP_CHUNK],
                )
            )
        return T_INF + out

    return ev


def _read(T, pts):
    vals, ok = evaluate_vector_field_parallel(T, pts)
    th = np.array(
        [
            float(np.real(np.ravel(vals[i])[0])) - T_INF if bool(ok[i]) else np.nan
            for i in range(pts.shape[0])
        ]
    )
    return th, [bool(o) for o in ok]


@pytest.fixture(scope="module")
def runs():
    comm = MPI.COMM_WORLD
    t0 = time.perf_counter()
    m = RUNGS[-1]
    ns = [b * m for b in BASE]
    msh = dmesh.create_box(
        comm,
        [np.array([-A_BOX, -B_BOX, -C_BOX]), np.array([A_BOX, B_BOX, C_BOX])],
        ns,
        cell_type=dmesh.CellType.tetrahedron,
    )
    ft = _robin_tags(msh)
    chi = _interior_chi(msh)
    pts = np.array(PROBE_PTS, dtype=float)
    _say(
        f"\nMAT21-26 rung m = {m} {tuple(ns)}: cells {msh.topology.index_map(3).size_global}; "
        f"g = {G_GRAD} K/m, g C = {GC} K; t_C = {T_C:.3f} s, t* = {T_STAR:.3f} s, q = {Q} W/m^3, "
        f"CN N = {N_MAIN} + 2 BE half-steps; ranks {comm.size}"
    )

    kw = dict(
        k=K_GEL,
        rho_c=RHO_C,
        q=Q,
        h=H_GEL,
        t_inf=T_INF,
        t_init=T_INF,
        robin_tag=ROBIN_TAG,
        dt=T_STAR / N_MAIN,
        n_steps=N_MAIN,
        theta=0.5,
        startup_half_steps=2,
        degree=1,
    )
    starts = {"ref": None, "tw0": lambda x: T_INF + G_GRAD * x[2], "twC": _wait_init(T_C)}
    out = {"E": {}, "src": {}, "th": {}, "ok": {}}
    res = {}
    for key, expr in starts.items():
        tq = time.perf_counter()
        res[key] = solve_transient_conduction(msh, ft, t_init_expr=expr, **kw)
        out["E"][key] = float(res[key].max_balance_residual)
        out["src"][key] = float(res[key].source_power)
        out["th"][key], out["ok"][key] = _read(res[key].temperature, pts)
        _say(
            f"MAT21-26 run {key}: {time.perf_counter() - tq:.1f} s; (E) {out['E'][key]:.3e} "
            f"(E_BAND {E_BAND:.0e}); source_power {out['src'][key]:.9e} W, / q V_box "
            f"{out['src'][key] / (Q * 8.0 * A_BOX * B_BOX * C_BOX) - 1.0:+.3e}"
        )

    tq = time.perf_counter()
    ex_s, npts = _exact(msh, _trans_ex(T_STAR))
    ex_0, _ = _exact(msh, _grad_ex(T_STAR))
    ex_c, _ = _exact(msh, _grad_ex(T_STAR + T_C))
    ex_dz, _ = _exact(msh, _grad_ex(T_STAR, roots_z=DZ_ROOTS))
    _say(
        f"MAT21-26 four series interpolations ({npts} P2 points): {time.perf_counter() - tq:.1f} s"
    )

    th = {key: r.temperature - T_INF for key, r in res.items()}
    out["ref_T"] = _rel(msh, th["ref"], ex_s, chi)
    out["I0"] = _rel(msh, th["tw0"], ex_s + ex_0, chi)
    out["IC"] = _rel(msh, th["twC"], ex_s + ex_c, chi)
    out["Dz"] = _rel(msh, th["tw0"], ex_s + ex_dz, chi)
    _say(
        f"MAT21-26 (I) t_w = 0: interior rel L2 {out['I0']:.6e}; t_w = T_C: {out['IC']:.6e} "
        f"(bound {INTERIOR_BAND:.0e}); whole {_rel(msh, th['tw0'], ex_s + ex_0):.6e} / "
        f"{_rel(msh, th['twC'], ex_s + ex_c):.6e} (printed)"
    )
    _say(
        f"MAT21-26 control (ref) reference interior {out['ref_T']:.6e} vs step 1 {STEP1_T:.6e}, "
        f"|delta| {abs(out['ref_T'] - STEP1_T):.3e} (bound {STEP1_T_BOUND:.0e}); (Dz) Dirichlet-z "
        f"interior miss {out['Dz']:.4e} (asserted > {INTERIOR_BAND:.0e})"
    )
    _say(f"MAT21-26 elapsed {time.perf_counter() - t0:.1f} s")
    return out


def test_w0_comparand_at_p5():
    """W0 (numpy): the imports and P51's comparand at P5, with no solve."""
    p = PROBE_PTS[5]
    d0, dc, s = _drift(p, T_STAR), _drift(p, T_STAR + T_C), _sourced(p, T_STAR)
    _say(
        f"MAT21-26 W0 P5: drift(t*) {d0:+.6e} K, drift(t* + T_C) {dc:+.6e} K, sourced "
        f"{s:.6f} K; T_C {T_C:.3f} s; {N_MAIN} steps; {len(DZ_ROOTS)} Dz roots"
    )
    assert N_MAIN == 40 and RUNGS[-1] == 6
    assert 0.0 < dc < d0 and s > 0.0


def test_i_both_gradient_runs(runs):
    assert runs["I0"] <= INTERIOR_BAND, runs["I0"]
    assert runs["IC"] <= INTERIOR_BAND, runs["IC"]


def test_e_all_three_on_the_source_branch(runs):
    src = runs["src"]
    assert src["ref"] > 0.0 and src["ref"] == src["tw0"] == src["twC"], src
    for key, e in runs["E"].items():
        assert e <= E_BAND, (key, e)


def test_sym_prime_on_the_drift_part(runs):
    assert all(all(o) for o in runs["ok"].values()), runs["ok"]
    worst = {}
    for key in ("tw0", "twC"):
        d = runs["th"][key] - runs["th"]["ref"]
        s12, s34, s0 = abs(d[1] + d[2]), abs(d[3] + d[4]), abs(d[0])
        worst[key] = max(s12, s34, s0)
        _say(
            f"MAT21-26 (Sym') {key}: |P1 + P2|/(g C) {s12 / GC:.3e}, |P3 + P4|/(g C) "
            f"{s34 / GC:.3e}, |P0|/(g C) {s0 / GC:.3e} (bar {SYM_BAR:.0e})"
        )
    for key, w in worst.items():
        assert w <= SYM_BAR * GC, (key, w)


def test_control_reference_reproduces_step_1(runs):
    assert abs(runs["ref_T"] - STEP1_T) <= STEP1_T_BOUND, runs["ref_T"]


def test_control_dirichlet_z_misses(runs):
    assert runs["Dz"] > INTERIOR_BAND, runs["Dz"]


def test_printed_bias_table(runs):
    """The record (no band); asserts only that it was formed."""
    _say(
        "MAT21-26 TABLE point | t_w | FEM drift K | P51 drift K | FEM sourced K | P51 sourced K | "
        "FEM bias | P51 bias | |FEM drift|/(g C)"
    )
    ref = runs["th"]["ref"]
    for key, t_w, tag in (("tw0", 0.0, "0"), ("twC", T_C, "T_C")):
        for i, (lab, p) in enumerate(zip(PROBE_LABELS, PROBE_PTS)):
            d_h = runs["th"][key][i] - ref[i]
            d_s, s_s = _drift(p, T_STAR + t_w), _sourced(p, T_STAR)
            _say(
                f"MAT21-26 TABLE {lab} | {tag} | {d_h:+.6e} | {d_s:+.6e} | {ref[i]:.6f} | "
                f"{s_s:.6f} | {d_h / ref[i]:+.6e} | {d_s / s_s:+.6e} | {abs(d_h) / GC:.3e}"
            )
            assert np.isfinite(d_h / ref[i])
